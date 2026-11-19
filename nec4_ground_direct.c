#include "nec4_ground_direct.h"

#include "nec2c.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#define NEC4_MU0 1.2566370621238373306e-6
#define NEC4_PI 3.14159265358979323846264338327950288

typedef struct {
  double frequency_hz;
  double complex lower_k;
  double complex upper_k;
  double rho;
  double observer_z;
  double source_z;
} direct_context;

typedef int (*density_function)(const direct_context *,double complex,
    double complex [NEC4_GROUND_COMPONENT_COUNT]);

typedef struct {
  direct_context direct;
  double complex previous_gamma[2];
  int have_previous[2];
} hankel_context;

static double complex complex_expm1(double complex value);

static double complex outgoing_gamma(double complex lambda,
    double complex wave_number)
{
  double complex value=csqrt(lambda*lambda-wave_number*wave_number);
  double tolerance=32.*DBL_EPSILON*fmax(1.,cabs(value));
  if(creal(value)<-tolerance||
     (fabs(creal(value))<=tolerance&&cimag(value)<0.)) value=-value;
  return value;
}

static int gauss_legendre(unsigned int count,double *nodes,double *weights)
{
  unsigned int index,iteration,mirror;
  if(count<12||!nodes||!weights) return 0;
  for(index=0;index<(count+1)/2;++index)
  {
    double root=cos(NEC4_PI*((double)index+.75)/((double)count+.5));
    double derivative=0.;
    for(iteration=0;iteration<32;++iteration)
    {
      unsigned int degree;
      double p0=1.,p1=root,pn;
      if(count==1) pn=p1;
      else
      {
        pn=p1;
        for(degree=2;degree<=count;++degree)
        {
          pn=((2.*degree-1.)*root*p1-(degree-1.)*p0)/degree;
          p0=p1;p1=pn;
        }
      }
      derivative=count*(root*pn-p0)/(root*root-1.);
      {
        double update=pn/derivative;
        root-=update;
        if(fabs(update)<=4.*DBL_EPSILON*fmax(1.,fabs(root))) break;
      }
    }
    if(iteration==32||!isfinite(root)||!isfinite(derivative)) return 0;
    {
      unsigned int degree;
      double p0=1.,p1=root,pn=p1;
      for(degree=2;degree<=count;++degree)
      {
        pn=((2.*degree-1.)*root*p1-(degree-1.)*p0)/degree;
        p0=p1;p1=pn;
      }
      derivative=count*(root*pn-p0)/(root*root-1.);
    }
    mirror=count-1-index;
    nodes[index]=-root;nodes[mirror]=root;
    weights[index]=2./((1.-root*root)*derivative*derivative);
    weights[mirror]=weights[index];
  }
  return 1;
}

static int density(const direct_context *context,double complex lambda,
    double complex value[NEC4_GROUND_COMPONENT_COUNT])
{
  double complex gamma1=outgoing_gamma(lambda,context->lower_k);
  double complex gamma2=outgoing_gamma(lambda,context->upper_k);
  double complex denominator=(context->lower_k*context->lower_k*gamma2+
      context->upper_k*context->upper_k*gamma1);
  double complex exponential=cexp(gamma1*context->source_z-
      gamma2*context->observer_z);
  double complex j0,j0_prime,j0_second,u,common;
  double complex argument=lambda*context->rho;
  if(cabs(denominator)==0.||cabs(gamma1+gamma2)==0.)
    return 0;
  if(context->rho==0.)
  {
    /* J0(x)=1-x^2/4+O(x^4): J0'(0)=0, J0''(0)=-1/2. */
    j0=1.;j0_prime=0.;j0_second=-.5;
  }
  else
  {
    bessel(argument,&j0,&j0_prime);
    j0_second=-j0-j0_prime/argument;
  }
  u=2.*exponential*j0*lambda/(gamma1+gamma2);
  common=I*context->frequency_hz*NEC4_MU0/2.;
  value[NEC4_GROUND_RHO_VERTICAL]=common*2.*gamma2*exponential*j0_prime*
    lambda*lambda/denominator;
  value[NEC4_GROUND_Z_VERTICAL]=-common*2.*exponential*j0*
    lambda*lambda*lambda/denominator;
  value[NEC4_GROUND_RHO_HORIZONTAL]=-common*(2.*exponential*j0_second*
    lambda*lambda*lambda/denominator+u);
  value[NEC4_GROUND_PHI_HORIZONTAL]=common*(context->rho==0.?
    -exponential*lambda*lambda*lambda/denominator+u:
    2.*exponential*j0_prime*lambda*lambda/
      (context->rho*denominator)+u);
  value[NEC4_GROUND_Z_HORIZONTAL]=common*2.*gamma1*exponential*j0_prime*
    lambda*lambda/denominator;
  return 1;
}

static double complex continuous_gamma(hankel_context *context,
    double complex lambda,double complex wave_number,unsigned int medium)
{
  double complex principal=csqrt(lambda*lambda-wave_number*wave_number);
  double complex selected;
  if(!context->have_previous[medium])
  {
    /* Appendix H's large-|lambda| continuation has
       gamma~S*lambda, S=sign(Re(lambda-k)).  Initializing against that
       asymptote avoids a sheet flip when a nearly vertical left ray moves
       by roundoff as the interface is approached. */
    double complex target=creal(lambda)<creal(wave_number)?-lambda:lambda;
    selected=cabs(-principal-target)<cabs(principal-target)?-principal:
      principal;
    context->have_previous[medium]=1;
  }
  else if(cabs(-principal-context->previous_gamma[medium])<
          cabs(principal-context->previous_gamma[medium])) selected=-principal;
  else selected=principal;
  context->previous_gamma[medium]=selected;
  return selected;
}

static int hankel_density(hankel_context *context,double complex lambda,
    double complex value[NEC4_GROUND_COMPONENT_COUNT])
{
  const direct_context *direct=&context->direct;
  double complex gamma1=continuous_gamma(context,lambda,direct->lower_k,0);
  double complex gamma2=continuous_gamma(context,lambda,direct->upper_k,1);
  double complex denominator=direct->lower_k*direct->lower_k*gamma2+
    direct->upper_k*direct->upper_k*gamma1;
  double complex gamma_sum=gamma1+gamma2;
  double complex exponential=cexp(gamma1*direct->source_z-
    gamma2*direct->observer_z);
  double complex h0,h0_prime,h0_second,u,common;
  double complex argument=lambda*direct->rho;
  if(cabs(denominator)==0.||cabs(gamma_sum)==0.||cabs(argument)==0.)
    return 0;
  hankel2(argument,&h0,&h0_prime);
  h0_second=-h0-h0_prime/argument;
  u=exponential*h0*lambda/gamma_sum;
  common=I*direct->frequency_hz*NEC4_MU0/2.;
  value[NEC4_GROUND_RHO_VERTICAL]=common*gamma2*exponential*h0_prime*
    lambda*lambda/denominator;
  value[NEC4_GROUND_Z_VERTICAL]=-common*exponential*h0*
    lambda*lambda*lambda/denominator;
  value[NEC4_GROUND_RHO_HORIZONTAL]=-common*(exponential*h0_second*
    lambda*lambda*lambda/denominator+u);
  value[NEC4_GROUND_PHI_HORIZONTAL]=common*(exponential*h0_prime*
    lambda*lambda/(direct->rho*denominator)+u);
  value[NEC4_GROUND_Z_HORIZONTAL]=common*gamma1*exponential*h0_prime*
    lambda*lambda/denominator;
  return 1;
}

static int appendix_h_hankel_density(hankel_context *context,
    double complex lambda,
    double complex value[NEC4_GROUND_COMPONENT_COUNT])
{
  const direct_context *direct=&context->direct;
  double complex k1_squared=direct->lower_k*direct->lower_k;
  double complex k2_squared=direct->upper_k*direct->upper_k;
  double complex sum_squared=k1_squared+k2_squared;
  double complex gamma1=continuous_gamma(context,lambda,direct->lower_k,0);
  double complex gamma2=continuous_gamma(context,lambda,direct->upper_k,1);
  double complex gamma_sum=gamma1+gamma2;
  double complex denominator=k1_squared*gamma2+k2_squared*gamma1;
  double complex delta,t1,t2,t3,propagation,v_coefficient,source_coefficient;
  double complex h0,h0_prime,h0_second,u,common,argument=lambda*direct->rho;
  if(cabs(argument)==0.||cabs(gamma_sum)==0.||cabs(gamma2)==0.||
     cabs(denominator)==0.||cabs(sum_squared)==0.) return 0;
  delta=(k2_squared-k1_squared)/gamma_sum;
  t1=complex_expm1(direct->source_z*delta);
  t2=gamma2*t1-delta;
  t3=gamma1*t1+delta;
  propagation=cexp(-gamma2*(direct->observer_z-direct->source_z));
  /* The full-line H2 representation contributes half of the positive-axis
     J0 density. */
  v_coefficient=(gamma2*k1_squared*t1+k2_squared*t2)/
    (gamma2*sum_squared*denominator)*propagation;
  source_coefficient=(gamma1*k2_squared*t1+k1_squared*t3)/
    (sum_squared*denominator)*propagation;
  hankel2(argument,&h0,&h0_prime);
  h0_second=-h0-h0_prime/argument;
  u=cexp(gamma1*direct->source_z-gamma2*direct->observer_z)*
    h0*lambda/gamma_sum;
  common=I*direct->frequency_hz*NEC4_MU0/2.;
  value[NEC4_GROUND_RHO_VERTICAL]=common*gamma2*v_coefficient*h0_prime*
    lambda*lambda;
  value[NEC4_GROUND_Z_VERTICAL]=-common*v_coefficient*h0*
    lambda*lambda*lambda;
  value[NEC4_GROUND_RHO_HORIZONTAL]=-common*(v_coefficient*h0_second*
    lambda*lambda*lambda+u);
  value[NEC4_GROUND_PHI_HORIZONTAL]=common*(v_coefficient*h0_prime*
    lambda*lambda/direct->rho+u);
  value[NEC4_GROUND_Z_HORIZONTAL]=common*source_coefficient*h0_prime*
    lambda*lambda;
  return 1;
}

static double complex complex_expm1(double complex value)
{
  if(cabs(value)<1.e-5)
  {
    double complex term=value,sum=value;
    unsigned int order;
    for(order=2;order<=12;++order)
    {
      term*=value/order;
      sum+=term;
    }
    return sum;
  }
  return cexp(value)-1.;
}

static int appendix_h_density(const direct_context *context,
    double complex lambda,
    double complex value[NEC4_GROUND_COMPONENT_COUNT])
{
  double complex k1_squared=context->lower_k*context->lower_k;
  double complex k2_squared=context->upper_k*context->upper_k;
  double complex sum_squared=k1_squared+k2_squared;
  double complex gamma1=outgoing_gamma(lambda,context->lower_k);
  double complex gamma2=outgoing_gamma(lambda,context->upper_k);
  double complex gamma_sum=gamma1+gamma2;
  double complex denominator=k1_squared*gamma2+k2_squared*gamma1;
  double complex delta,t1,t2,t3,propagation,v_coefficient,source_coefficient;
  double complex j0,j0_prime,j0_second,u,common,argument=lambda*context->rho;
  if(cabs(gamma_sum)==0.||cabs(gamma2)==0.||cabs(denominator)==0.||
     cabs(sum_squared)==0.) return 0;
  delta=(k2_squared-k1_squared)/gamma_sum;
  t1=complex_expm1(context->source_z*delta);
  t2=gamma2*t1-delta;
  t3=gamma1*t1+delta;
  propagation=cexp(-gamma2*(context->observer_z-context->source_z));
  v_coefficient=2.*(gamma2*k1_squared*t1+k2_squared*t2)/
    (gamma2*sum_squared*denominator)*propagation;
  source_coefficient=2.*(gamma1*k2_squared*t1+k1_squared*t3)/
    (sum_squared*denominator)*propagation;
  if(context->rho==0.)
  {
    /* The same cylindrical limits apply after Appendix-H subtraction. */
    j0=1.;j0_prime=0.;j0_second=-.5;
  }
  else
  {
    bessel(argument,&j0,&j0_prime);
    j0_second=-j0-j0_prime/argument;
  }
  u=2.*cexp(gamma1*context->source_z-gamma2*context->observer_z)*
    j0*lambda/gamma_sum;
  common=I*context->frequency_hz*NEC4_MU0/2.;
  value[NEC4_GROUND_RHO_VERTICAL]=common*gamma2*v_coefficient*j0_prime*
    lambda*lambda;
  value[NEC4_GROUND_Z_VERTICAL]=-common*v_coefficient*j0*
    lambda*lambda*lambda;
  value[NEC4_GROUND_RHO_HORIZONTAL]=-common*(v_coefficient*j0_second*
    lambda*lambda*lambda+u);
  value[NEC4_GROUND_PHI_HORIZONTAL]=common*(context->rho==0.?
    -.5*v_coefficient*lambda*lambda*lambda+u:
    v_coefficient*j0_prime*lambda*lambda/context->rho+u);
  value[NEC4_GROUND_Z_HORIZONTAL]=common*source_coefficient*j0_prime*
    lambda*lambda;
  return 1;
}

static int integrate_segment(const direct_context *context,
    double complex start,double complex end,unsigned int order,
    const double *nodes,const double *weights,
    density_function evaluate_density,
    double complex sum[NEC4_GROUND_COMPONENT_COUNT])
{
  unsigned int node,component;
  double complex midpoint=(start+end)/2.;
  double complex half_width=(end-start)/2.;
  for(node=0;node<order;++node)
  {
    double complex values[NEC4_GROUND_COMPONENT_COUNT];
    double complex lambda=midpoint+half_width*nodes[node];
    if(!evaluate_density(context,lambda,values)) return 0;
    for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
      sum[component]+=half_width*weights[node]*values[component];
  }
  return 1;
}

int nec4_ground_direct_cross_bessel(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int quadrature_order,
    double cutoff_decay,nec4_ground_direct_result *result)
{
  direct_context context;
  double height_sum,scale,panel_width,cutoff,left,right;
  double complex p0;
  double *nodes=NULL,*weights=NULL;
  unsigned int panel,panels,component;
  if(!result||!(frequency_hz>0.)||rho<0.||observer_z<0.||source_depth<0.||
     observer_z+source_depth<=0.||quadrature_order<12||cutoff_decay<30.||
     !isfinite(frequency_hz)||!isfinite(rho)||!isfinite(observer_z)||
     !isfinite(source_depth)||!(cabs(lower_k)>0.)||!(cabs(upper_k)>0.))
    return 0;
  height_sum=observer_z+source_depth;
  scale=rho==0.?1./height_sum:fmin(1./rho,1./height_sum);
  p0=scale+I*scale;
  panel_width=NEC4_PI*scale;
  cutoff=scale+cutoff_decay/height_sum;
  panels=(unsigned int)ceil((cutoff-scale)/panel_width);
  nodes=malloc(quadrature_order*sizeof(*nodes));
  weights=malloc(quadrature_order*sizeof(*weights));
  if(!nodes||!weights||!gauss_legendre(quadrature_order,nodes,weights))
  { free(nodes);free(weights);return 0; }
  context.frequency_hz=frequency_hz;
  context.lower_k=lower_k;context.upper_k=upper_k;
  context.rho=rho;context.observer_z=observer_z;context.source_z=-source_depth;
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    result->value[component]=0.;
  if(!integrate_segment(&context,0.,p0,quadrature_order,nodes,weights,density,
      result->value)) { free(nodes);free(weights);return 0; }
  left=scale;
  for(panel=0;panel<panels;++panel)
  {
    right=fmin(left+panel_width,cutoff);
    if(!integrate_segment(&context,left+I*scale,right+I*scale,
        quadrature_order,nodes,weights,density,result->value))
    { free(nodes);free(weights);return 0; }
    left=right;
  }
  free(nodes);free(weights);
  result->panels=panels+1;
  result->quadrature_order=quadrature_order;
  return 1;
}

typedef int (*hankel_density_function)(hankel_context *,double complex,
    double complex [NEC4_GROUND_COMPONENT_COUNT]);

static int integrate_hankel_segment(hankel_context *context,
    double complex start,double complex end,unsigned int order,
    const double *nodes,const double *weights,
    hankel_density_function evaluate_density,
    double complex sum[NEC4_GROUND_COMPONENT_COUNT])
{
  unsigned int node,component;
  double complex midpoint=(start+end)/2.,half_width=(end-start)/2.;
  for(node=0;node<order;++node)
  {
    double complex values[NEC4_GROUND_COMPONENT_COUNT];
    double complex lambda=midpoint+half_width*nodes[node];
    if(!evaluate_density(context,lambda,values)) return 0;
    for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
      sum[component]+=half_width*weights[node]*values[component];
  }
  return 1;
}

static int integrate_hankel_line(hankel_context *context,
    double complex start,double complex end,unsigned int panels,
    unsigned int order,const double *nodes,const double *weights,
    hankel_density_function evaluate_density,
    double complex sum[NEC4_GROUND_COMPONENT_COUNT],unsigned int *total)
{
  unsigned int panel;
  if(panels<1) panels=1;
  for(panel=0;panel<panels;++panel)
    if(!integrate_hankel_segment(context,
        start+(end-start)*(double)panel/panels,
        start+(end-start)*(double)(panel+1)/panels,
        order,nodes,weights,evaluate_density,sum)) return 0;
  *total+=panels;
  return 1;
}

static unsigned int line_panel_count(double complex start,double complex end,
    double panel_length,unsigned int minimum)
{
  double requested=ceil(cabs(end-start)/panel_length);
  unsigned int count=requested>(double)UINT_MAX?UINT_MAX:(unsigned int)requested;
  return count<minimum?minimum:count;
}

static int cross_hankel_integrate(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int panel_subdivisions,
    unsigned int quadrature_order,double cutoff_decay,
    hankel_density_function evaluate_density,
    double complex result[NEC4_GROUND_COMPONENT_COUNT],
    unsigned int *used_panels)
{
  hankel_context context;
  double height_sum,radial_distance,upper_real,d_min,d_cap;
  double panel_length,ray_parameter;
  double complex p1,p2,p3,left_end,right_end,delta;
  double fractions[4]={0.,1.};
  double *nodes=NULL,*weights=NULL;
  unsigned int fraction_count=2,index,j,total_panels=0,ray_panels;
  if(!result||!used_panels||!evaluate_density||!(frequency_hz>0.)||
     !(rho>0.)||observer_z<0.||source_depth<0.||
     !(observer_z+source_depth>0.)||
     panel_subdivisions<1||quadrature_order<12||cutoff_decay<30.||
     !isfinite(frequency_hz)||!isfinite(rho)||!isfinite(observer_z)||
     !isfinite(source_depth)||!(cabs(lower_k)>0.)||!(cabs(upper_k)>0.))
    return 0;
  if(fabs(cimag(upper_k))>128.*DBL_EPSILON*fmax(1.,cabs(upper_k))) return 0;
  height_sum=observer_z+source_depth;
  if(height_sum>.6*rho) return 0;
  radial_distance=hypot(rho,height_sum);upper_real=creal(upper_k);
  if(!(upper_real>0.)) return 0;
  d_min=upper_real*fmin(.06/(upper_real*radial_distance)+.04,
      .8/(1.+height_sum/rho));
  d_cap=creal(lower_k+upper_k)/2.-upper_real;
  d_min=fmin(d_min,d_cap);
  if(!(d_min>0.)) return 0;
  p1=upper_k+(-1.+I)*d_min;p2=upper_k+(1.+I)*d_min;
  p3=fmax(1.01*creal(lower_k),1.2*upper_real)+I*.95*cimag(lower_k);
  {
    double slope=rho/height_sum;
    double asymptotic_decay=height_sum+rho*rho/height_sum;
    ray_parameter=cutoff_decay/asymptotic_decay;
    left_end=p1-ray_parameter*(1.+I*slope);
    right_end=p3+ray_parameter*(1.-I*slope);
  }
  panel_length=NEC4_PI/(fmax(rho,height_sum)*panel_subdivisions);
  ray_panels=line_panel_count(p3,right_end,panel_length,
                              4*panel_subdivisions);
  nodes=malloc(quadrature_order*sizeof(*nodes));
  weights=malloc(quadrature_order*sizeof(*weights));
  if(!nodes||!weights||!gauss_legendre(quadrature_order,nodes,weights))
  { free(nodes);free(weights);return 0; }
  memset(&context,0,sizeof(context));
  context.direct.frequency_hz=frequency_hz;
  context.direct.lower_k=lower_k;context.direct.upper_k=upper_k;
  context.direct.rho=rho;context.direct.observer_z=observer_z;
  context.direct.source_z=-source_depth;
  for(index=0;index<NEC4_GROUND_COMPONENT_COUNT;++index)
    result[index]=0.;

  delta=p1-left_end;
  if(creal(delta)!=0.)
  {
    double fraction=-creal(left_end)/creal(delta);
    if(fraction>0.&&fraction<1.) fractions[fraction_count++]=fraction;
  }
  if(cimag(delta)!=0.)
  {
    double fraction=-cimag(left_end)/cimag(delta);
    if(fraction>0.&&fraction<1.) fractions[fraction_count++]=fraction;
  }
  for(index=0;index<fraction_count;++index)
    for(j=index+1;j<fraction_count;++j)
      if(fractions[j]<fractions[index])
      { double swap=fractions[index];fractions[index]=fractions[j];
        fractions[j]=swap; }
  {
    unsigned int unique=1;
    for(index=1;index<fraction_count;++index)
      if(fractions[index]!=fractions[unique-1])
        fractions[unique++]=fractions[index];
    fraction_count=unique;
  }
  for(index=0;index+1<fraction_count;++index)
  {
    double span=fractions[index+1]-fractions[index];
    unsigned int panels=(unsigned int)ceil(ray_panels*span);
    if(panels<2*panel_subdivisions) panels=2*panel_subdivisions;
    if(!integrate_hankel_line(&context,left_end+delta*fractions[index],
        left_end+delta*fractions[index+1],panels,quadrature_order,
        nodes,weights,evaluate_density,result,&total_panels)) goto failure;
  }
  if(!integrate_hankel_line(&context,p1,p2,
      line_panel_count(p1,p2,panel_length,2*panel_subdivisions),
      quadrature_order,nodes,weights,evaluate_density,result,
      &total_panels)) goto failure;
  if(!integrate_hankel_line(&context,p2,p3,
      line_panel_count(p2,p3,panel_length,2*panel_subdivisions),
      quadrature_order,nodes,weights,evaluate_density,result,
      &total_panels)) goto failure;
  if(!integrate_hankel_line(&context,p3,right_end,ray_panels,
      quadrature_order,nodes,weights,evaluate_density,result,
      &total_panels)) goto failure;
  free(nodes);free(weights);
  *used_panels=total_panels;
  return 1;
failure:
  free(nodes);free(weights);return 0;
}

int nec4_ground_direct_cross_hankel(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int panel_subdivisions,
    unsigned int quadrature_order,double cutoff_decay,
    nec4_ground_direct_result *result)
{
  if(!result||!cross_hankel_integrate(frequency_hz,lower_k,upper_k,rho,
      observer_z,source_depth,panel_subdivisions,quadrature_order,cutoff_decay,
      hankel_density,result->value,&result->panels)) return 0;
  result->quadrature_order=quadrature_order;
  return 1;
}

static void appendix_h_addbacks(const direct_context *context,
    const double complex first[NEC4_GROUND_COMPONENT_COUNT],
    nec4_ground_table1_result *result)
{
  double complex k1_squared=context->lower_k*context->lower_k;
  double complex k2_squared=context->upper_k*context->upper_k;
  double complex sum_squared=k1_squared+k2_squared;
  double dz=context->observer_z-context->source_z;
  double distance=hypot(context->rho,dz);
  double sine=dz/distance,cosine=context->rho/distance;
  double complex wave=cexp(-I*context->upper_k*distance);
  double complex radial=3./(distance*distance)+3.*I*context->upper_k/distance-
    k2_squared;
  double complex v0_factor=(-I*context->frequency_hz*NEC4_MU0/2.)*
    2./sum_squared*wave/(distance*distance*distance);
  double complex c1=(k1_squared-k2_squared)/sum_squared;
  double complex c2=k2_squared*(k1_squared-k2_squared)/
    (sum_squared*sum_squared);
  double complex negative=(-I*context->frequency_hz*NEC4_MU0/2.)*
    wave/distance;
  double complex positive=-negative;
  double one_minus_sine=1.-sine;
  double complex v0[NEC4_GROUND_COMPONENT_COUNT];
  double complex residual[NEC4_GROUND_COMPONENT_COUNT];
  unsigned int component;
  v0[NEC4_GROUND_RHO_VERTICAL]=v0_factor*context->rho*dz*radial;
  v0[NEC4_GROUND_Z_VERTICAL]=v0_factor*(dz*dz*radial-1.-
    I*context->upper_k*distance+k2_squared*distance*distance);
  v0[NEC4_GROUND_RHO_HORIZONTAL]=v0_factor*(context->rho*context->rho*
    radial-1.-I*context->upper_k*distance);
  v0[NEC4_GROUND_PHI_HORIZONTAL]=v0_factor*(1.+
    I*context->upper_k*distance);
  v0[NEC4_GROUND_Z_HORIZONTAL]=v0_factor*context->rho*dz*radial;
  if(context->rho==0.)
  {
    /* With c=rho/R and s=sqrt(1-c^2), (1-s)/c^2 -> 1/2 and
       (s*(1+c^2)-1)/c^2 -> 1/2.  Evaluating these limits explicitly
       avoids an arbitrary off-axis displacement in the table generator. */
    double complex transverse=-.5*c2+
      .5*c1*context->source_z/distance+1.;
    residual[NEC4_GROUND_RHO_VERTICAL]=0.;
    residual[NEC4_GROUND_Z_VERTICAL]=negative*(c2-
      c1*context->source_z/distance);
    residual[NEC4_GROUND_RHO_HORIZONTAL]=negative*transverse;
    residual[NEC4_GROUND_PHI_HORIZONTAL]=positive*transverse;
    residual[NEC4_GROUND_Z_HORIZONTAL]=0.;
  }
  else
  {
    residual[NEC4_GROUND_RHO_VERTICAL]=negative*(c2*one_minus_sine/cosine-
      c1*context->source_z*cosine/distance);
    residual[NEC4_GROUND_Z_VERTICAL]=negative*(c2-
      c1*context->source_z*sine/distance);
    residual[NEC4_GROUND_RHO_HORIZONTAL]=negative*(
      c2*(one_minus_sine/(cosine*cosine)-1.)+
      c1*context->source_z/distance*
      ((sine*(1.+cosine*cosine)-1.)/(cosine*cosine))+1.);
    residual[NEC4_GROUND_PHI_HORIZONTAL]=positive*(
      -c2*one_minus_sine/(cosine*cosine)+
      c1*context->source_z*one_minus_sine/(distance*cosine*cosine)+1.);
    residual[NEC4_GROUND_Z_HORIZONTAL]=positive*(
      k1_squared*c2/k2_squared*one_minus_sine/cosine+
      c1*context->source_z*cosine/distance);
  }
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
  {
    result->total[component]=first[component]+v0[component];
    result->finite_remainder[component]=first[component]-residual[component];
  }
}

int nec4_ground_direct_cross_table1_hankel(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int panel_subdivisions,
    unsigned int quadrature_order,double cutoff_decay,
    nec4_ground_table1_result *result)
{
  direct_context context;
  double complex first[NEC4_GROUND_COMPONENT_COUNT];
  if(!result||!cross_hankel_integrate(frequency_hz,lower_k,upper_k,rho,
      observer_z,source_depth,panel_subdivisions,quadrature_order,cutoff_decay,
      appendix_h_hankel_density,first,&result->panels)) return 0;
  context.frequency_hz=frequency_hz;
  context.lower_k=lower_k;context.upper_k=upper_k;
  context.rho=rho;context.observer_z=observer_z;
  context.source_z=-source_depth;
  appendix_h_addbacks(&context,first,result);
  result->quadrature_order=quadrature_order;
  return 1;
}

static double interface_error(
    const double complex actual[NEC4_GROUND_COMPONENT_COUNT],
    const double complex expected[NEC4_GROUND_COMPONENT_COUNT])
{
  double scale=0.,error=0.;
  unsigned int component;
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    scale=fmax(scale,cabs(expected[component]));
  if(!(scale>0.)) scale=1.;
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    error=fmax(error,cabs(actual[component]-expected[component])/scale);
  return error;
}

int nec4_ground_direct_cross_interface(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    unsigned int panel_subdivisions,unsigned int quadrature_order,
    double cutoff_decay,nec4_ground_interface_result *result)
{
  const double coarse_ratio=.002;
  nec4_ground_table1_result symmetric[3],source_only[2];
  double complex coarse[NEC4_GROUND_COMPONENT_COUNT];
  double complex path[NEC4_GROUND_COMPONENT_COUNT];
  unsigned int component,index,total_panels=0;
  if(!result||!(rho>0.)) return 0;
  for(index=0;index<3;++index)
  {
    double height=coarse_ratio*rho/(1u<<index);
    if(!nec4_ground_direct_cross_table1_hankel(frequency_hz,lower_k,upper_k,
        rho,height/2.,height/2.,panel_subdivisions,quadrature_order,
        cutoff_decay,&symmetric[index])) return 0;
    total_panels+=symmetric[index].panels;
  }
  for(index=0;index<2;++index)
  {
    double height=coarse_ratio*rho/(2u<<index);
    if(!nec4_ground_direct_cross_table1_hankel(frequency_hz,lower_k,upper_k,
        rho,0.,height,panel_subdivisions,quadrature_order,cutoff_decay,
        &source_only[index])) return 0;
    total_panels+=source_only[index].panels;
  }
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
  {
    coarse[component]=2.*symmetric[1].finite_remainder[component]-
      symmetric[0].finite_remainder[component];
    result->field.finite_remainder[component]=
      2.*symmetric[2].finite_remainder[component]-
      symmetric[1].finite_remainder[component];
    result->field.total[component]=2.*symmetric[2].total[component]-
      symmetric[1].total[component];
    path[component]=2.*source_only[1].finite_remainder[component]-
      source_only[0].finite_remainder[component];
  }
  result->refinement_error=interface_error(
      result->field.finite_remainder,coarse);
  result->approach_path_error=interface_error(
      result->field.finite_remainder,path);
  result->field.panels=total_panels;
  result->field.quadrature_order=quadrature_order;
  return 1;
}

int nec4_ground_direct_table1_restore(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,
    const double complex finite_remainder[NEC4_GROUND_COMPONENT_COUNT],
    double complex total[NEC4_GROUND_COMPONENT_COUNT])
{
  direct_context context;
  nec4_ground_table1_result addbacks;
  double complex zero[NEC4_GROUND_COMPONENT_COUNT]={0};
  unsigned int component;
  if(!finite_remainder||!total||!(frequency_hz>0.)||rho<0.||
     observer_z<0.||source_depth<0.||
     !(rho>0.||observer_z+source_depth>0.)||
     !isfinite(frequency_hz)||!isfinite(rho)||!isfinite(observer_z)||
     !isfinite(source_depth)||!(cabs(lower_k)>0.)||!(cabs(upper_k)>0.))
    return 0;
  context.frequency_hz=frequency_hz;
  context.lower_k=lower_k;context.upper_k=upper_k;
  context.rho=rho;context.observer_z=observer_z;
  context.source_z=-source_depth;
  appendix_h_addbacks(&context,zero,&addbacks);
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
  {
    if(!isfinite(creal(finite_remainder[component]))||
       !isfinite(cimag(finite_remainder[component]))) return 0;
    /* addbacks.total is V0 and addbacks.finite_remainder is -S'. */
    total[component]=finite_remainder[component]-
      addbacks.finite_remainder[component]+addbacks.total[component];
  }
  return 1;
}

int nec4_ground_direct_cross_table1(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int quadrature_order,
    double cutoff_decay,nec4_ground_table1_result *result)
{
  direct_context context;
  double height_sum,scale,panel_width,cutoff,left,right;
  double complex p0,first[NEC4_GROUND_COMPONENT_COUNT]={0};
  double *nodes=NULL,*weights=NULL;
  unsigned int panel,panels;
  if(!result||!(frequency_hz>0.)||rho<0.||observer_z<0.||source_depth<0.||
     observer_z+source_depth<=0.||quadrature_order<12||cutoff_decay<30.||
     !isfinite(frequency_hz)||!isfinite(rho)||!isfinite(observer_z)||
     !isfinite(source_depth)||!(cabs(lower_k)>0.)||!(cabs(upper_k)>0.))
    return 0;
  context.frequency_hz=frequency_hz;
  context.lower_k=lower_k;context.upper_k=upper_k;
  context.rho=rho;context.observer_z=observer_z;context.source_z=-source_depth;
  height_sum=observer_z+source_depth;
  scale=rho==0.?1./height_sum:fmin(1./rho,1./height_sum);
  p0=scale+I*scale;
  panel_width=NEC4_PI*scale;
  cutoff=scale+cutoff_decay/height_sum;
  panels=(unsigned int)ceil((cutoff-scale)/panel_width);
  nodes=malloc(quadrature_order*sizeof(*nodes));
  weights=malloc(quadrature_order*sizeof(*weights));
  if(!nodes||!weights||!gauss_legendre(quadrature_order,nodes,weights))
  { free(nodes);free(weights);return 0; }
  if(!integrate_segment(&context,0.,p0,quadrature_order,nodes,weights,
      appendix_h_density,first)) { free(nodes);free(weights);return 0; }
  left=scale;
  for(panel=0;panel<panels;++panel)
  {
    right=fmin(left+panel_width,cutoff);
    if(!integrate_segment(&context,left+I*scale,right+I*scale,
        quadrature_order,nodes,weights,appendix_h_density,first))
    { free(nodes);free(weights);return 0; }
    left=right;
  }
  free(nodes);free(weights);
  appendix_h_addbacks(&context,first,result);
  result->panels=panels+1;
  result->quadrature_order=quadrature_order;
  return 1;
}
