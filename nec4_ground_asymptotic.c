#include "nec4_ground_asymptotic.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define JET_ORDER 4
#define MU0 (4.e-7*M_PI)

typedef struct { double complex c[JET_ORDER+1]; } jet;

static double complex outgoing_beta(double complex k,double complex lambda)
{
  double complex value=csqrt(k*k-lambda*lambda);
  if(creal(value)<0. ||
     (fabs(creal(value))<1.e-14*fmax(1.,cabs(value)) && cimag(value)>0.))
    value=-value;
  return value;
}

static int saddle_point(double complex k1,double complex k2,double rho,
    double z,double depth,double complex *answer)
{
  unsigned int iteration;
  double scale;
  double complex value;
  if(!answer || rho<=0. || z<0. || depth<0.) return 0;
  if(depth==0.) {
    *answer=k2*rho/hypot(rho,z);
    return cabs(*answer)>0.;
  }
  if(z==0.) { *answer=k2; return cabs(*answer)>0.; }
  if(fabs(cimag(k1))<=1.e-14*fmax(1.,cabs(k1)) &&
     fabs(cimag(k2))<=1.e-14*fmax(1.,cabs(k2)) && creal(k2)>0.) {
    double low=0.,high=creal(k2)*(1.-2.e-14);
    for(iteration=0;iteration<80;iteration++) {
      double middle=.5*(low+high);
      double beta1=sqrt(fmax(0.,creal(k1)*creal(k1)-middle*middle));
      double beta2=sqrt(fmax(0.,creal(k2)*creal(k2)-middle*middle));
      double residual=rho-middle*(z/beta2+depth/beta1);
      if(residual>0.) low=middle; else high=middle;
    }
    *answer=.5*(low+high);
    return 1;
  }
  scale=cabs(k2/k1);
  value=k2*rho/hypot(rho,z+scale*depth);
  for(iteration=0;iteration<40;iteration++) {
    double complex beta1=outgoing_beta(k1,value);
    double complex beta2=outgoing_beta(k2,value);
    double complex residual,derivative,update;
    if(cabs(beta1)==0. || cabs(beta2)==0.) { value*=1.-1.e-8; continue; }
    residual=rho-value*(z/beta2+depth/beta1);
    derivative=-(z*k2*k2/(beta2*beta2*beta2)+
                 depth*k1*k1/(beta1*beta1*beta1));
    if(cabs(derivative)==0.) break;
    update=residual/derivative;
    value-=update;
    if(cabs(update)<=2.e-13*fmax(1.,cabs(value))) break;
  }
  {
    double complex b1=outgoing_beta(k1,value),b2=outgoing_beta(k2,value);
    double complex residual=rho-value*(z/b2+depth/b1);
    if(cabs(residual)>2.e-10*fmax(1.,rho)) return 0;
  }
  if(creal(value)>creal(k2)*(1.+1.e-10)) value=k2;
  *answer=value;
  return 1;
}

static jet jet_constant(double complex value)
{
  jet result={{0.}}; result.c[0]=value; return result;
}

static jet jet_variable(double complex value)
{
  jet result=jet_constant(value); result.c[1]=1.; return result;
}

static jet jet_add(jet a,jet b)
{
  jet result; unsigned int n;
  for(n=0;n<=JET_ORDER;n++) result.c[n]=a.c[n]+b.c[n];
  return result;
}

static jet jet_scale(jet a,double complex value)
{
  jet result; unsigned int n;
  for(n=0;n<=JET_ORDER;n++) result.c[n]=value*a.c[n];
  return result;
}

static jet jet_mul(jet a,jet b)
{
  jet result={{0.}}; unsigned int n,k;
  for(n=0;n<=JET_ORDER;n++)
    for(k=0;k<=n;k++) result.c[n]+=a.c[k]*b.c[n-k];
  return result;
}

static jet jet_inv(jet a)
{
  jet result={{0.}}; unsigned int n,k;
  result.c[0]=1./a.c[0];
  for(n=1;n<=JET_ORDER;n++) {
    for(k=1;k<=n;k++) result.c[n]+=a.c[k]*result.c[n-k];
    result.c[n]/=-a.c[0];
  }
  return result;
}

static jet jet_div(jet a,jet b) { return jet_mul(a,jet_inv(b)); }

static jet jet_sqrt_branch(jet a,double complex root)
{
  jet result={{0.}}; unsigned int n,k;
  result.c[0]=root;
  for(n=1;n<=JET_ORDER;n++) {
    double complex sum=0.;
    for(k=1;k<n;k++) sum+=result.c[k]*result.c[n-k];
    result.c[n]=(a.c[n]-sum)/(2.*root);
  }
  return result;
}

static jet beta_jet(double complex k,jet lambda)
{
  jet radicand=jet_add(jet_constant(k*k),jet_scale(jet_mul(lambda,lambda),-1.));
  return jet_sqrt_branch(radicand,outgoing_beta(k,lambda.c[0]));
}

static int build_jets(double frequency_hz,double complex k1,double complex k2,
    double rho,double z,double depth,double complex saddle,jet *phase,jet field[5])
{
  jet lambda=jet_variable(saddle);
  jet beta1=beta_jet(k1,lambda),beta2=beta_jet(k2,lambda);
  jet gamma1=jet_scale(beta1,I),gamma2=jet_scale(beta2,I);
  jet denominator=jet_add(jet_scale(gamma2,k1*k1),
                          jet_scale(gamma1,k2*k2));
  jet root_arg=jet_div(jet_constant(2./(M_PI*rho)),lambda);
  jet root=jet_sqrt_branch(root_arg,csqrt(root_arg.c[0]));
  jet correction=jet_scale(jet_inv(lambda),1./(8.*rho));
  jet h0=jet_scale(jet_mul(root,jet_add(jet_constant(1.),
                                        jet_scale(correction,I))),
                   cexp(I*M_PI/4.));
  jet h0p=jet_scale(jet_mul(root,jet_add(jet_constant(1.),
                                         jet_scale(correction,-3.*I))),
                    cexp(-I*M_PI/4.));
  jet h0pp=jet_scale(jet_mul(root,jet_add(jet_constant(1.),
                                          jet_scale(correction,-7.*I))),
                     -cexp(I*M_PI/4.));
  jet lambda2=jet_mul(lambda,lambda),lambda3=jet_mul(lambda2,lambda);
  jet q=jet_constant(I*(2.*M_PI*frequency_hz)*MU0/(4.*M_PI));
  jet u=jet_div(jet_mul(h0,lambda),jet_add(gamma1,gamma2));
  *phase=jet_add(jet_add(jet_scale(gamma2,z),jet_scale(gamma1,depth)),
                 jet_scale(lambda,I*rho));
  field[NEC4_GROUND_RHO_VERTICAL]=jet_mul(q,jet_div(
      jet_mul(jet_mul(gamma2,h0p),lambda2),denominator));
  field[NEC4_GROUND_Z_VERTICAL]=jet_scale(jet_mul(q,
      jet_div(jet_mul(h0,lambda3),denominator)),-1.);
  field[NEC4_GROUND_RHO_HORIZONTAL]=jet_scale(jet_mul(q,jet_add(
      jet_div(jet_mul(h0pp,lambda3),denominator),u)),-1.);
  field[NEC4_GROUND_PHI_HORIZONTAL]=jet_mul(q,jet_add(
      jet_scale(jet_div(jet_mul(h0p,lambda2),denominator),1./rho),u));
  field[NEC4_GROUND_Z_HORIZONTAL]=jet_mul(q,jet_div(
      jet_mul(jet_mul(gamma1,h0p),lambda2),denominator));
  return 1;
}

int nec4_ground_asymptotic_first(double frequency_hz,
    double complex k1,double complex k2,double rho,double z,double depth,
    nec4_ground_asymptotic_result *result)
{
  double complex saddle,gamma1,gamma2,denominator,curvature,spreading,psi;
  double scale;
  if(!result || frequency_hz<=0. || rho<=0. || z<0. || depth<0. ||
     !saddle_point(k1,k2,rho,z,depth,&saddle)) return 0;
  gamma1=I*outgoing_beta(k1,saddle); gamma2=I*outgoing_beta(k2,saddle);
  if(cabs(saddle)==0. || cabs(gamma1)==0. || cabs(gamma2)==0.) return 0;
  denominator=k1*k1*gamma2+k2*k2*gamma1;
  curvature=-k1*k1*gamma2*gamma2*depth/gamma1-
             k2*k2*gamma1*gamma1*z/gamma2;
  spreading=-csqrt(rho/saddle)*csqrt(curvature);
  if(cabs(denominator)==0. || cabs(spreading)==0.) return 0;
  psi=gamma1*gamma2*cexp(-gamma1*depth-gamma2*z-I*saddle*rho)/spreading;
  scale=(2.*M_PI*frequency_hz)*MU0/(2.*M_PI);
  result->value[0]=cexp(I*M_PI/4.)*scale*saddle*gamma2/denominator*psi;
  result->value[1]=cexp(-I*M_PI/4.)*scale*saddle*saddle*saddle/denominator*psi;
  result->value[2]=cexp(-I*M_PI/4.)*scale*(-saddle*saddle/denominator+
      1./(gamma1+gamma2))*psi;
  result->value[3]=cexp(-I*M_PI/4.)*scale*(I*saddle/(rho*denominator)-
      1./(gamma1+gamma2))*psi;
  result->value[4]=cexp(I*M_PI/4.)*scale*saddle*gamma1/denominator*psi;
  result->saddle=saddle;
  return 1;
}

static int second_saddle_raw(double frequency_hz,
    double complex k1,double complex k2,double rho,double z,double depth,
    nec4_ground_asymptotic_result *result)
{
  nec4_ground_asymptotic_result first;
  double complex saddle,prefactor,f2,f3,f4,probe;
  jet phase,field[5];
  unsigned int component;
  if(!result || !nec4_ground_asymptotic_first(frequency_hz,k1,k2,rho,z,depth,
                                               &first)) return 0;
  saddle=first.saddle;
  if(!build_jets(frequency_hz,k1,k2,rho,z,depth,saddle,&phase,field)) return 0;
  f2=2.*phase.c[2]; f3=6.*phase.c[3]; f4=24.*phase.c[4];
  if(cabs(f2)==0.) return 0;
  prefactor=csqrt(2.*M_PI/f2)*cexp(-phase.c[0]);
  probe=prefactor*field[NEC4_GROUND_RHO_VERTICAL].c[0];
  if(cabs(probe-first.value[0])>cabs(-probe-first.value[0])) prefactor=-prefactor;
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;component++) {
    double complex g0=field[component].c[0];
    double complex g1=field[component].c[1];
    double complex g2=2.*field[component].c[2];
    double complex q2=(g0*(5.*f3*f3-3.*f2*f4)-12.*g1*f2*f3+
                       12.*g2*f2*f2)/(24.*f2*f2*f2);
    result->value[component]=prefactor*(g0+q2);
  }
  result->saddle=saddle;
  return 1;
}

/* Appendix J equations (J-4)--(J-13) after taking z -> 0 for a source
   at or below the interface.  The ordinary saddle reaches the k2 branch
   point in this limit, so evaluating the generic derivatives loses the
   finite products of their singular factors.  Burke's Lambda2=Gamma2/z
   form preserves those products without displacing the observation point. */
static int second_interface_raw(double frequency_hz,
    double complex k1,double complex k2,double rho,double depth,
    nec4_ground_asymptotic_result *result)
{
  double complex lambda=k2,gamma1=I*outgoing_beta(k1,k2),lambda2,base;
  double complex f1,f2,denominator,dv,dv1,dv2,du1,du2;
  double complex root,l32,l52,rv,zv,rh,ph,zh,u,phase,scale;
  double complex plus=cexp(I*M_PI/4.),minus=cexp(-I*M_PI/4.);
  double common;
  if(!result || frequency_hz<=0. || rho<=0. || depth<0. ||
     cabs(gamma1)==0.) return 0;
  lambda2=lambda/(lambda*(-depth)/gamma1-I*rho);
  if(cabs(lambda2)==0.) return 0;
  base=-k2*k2*gamma1*gamma1*gamma1;
  if(cabs(base)==0.) return 0;

  f1=-3.*lambda*cpow(gamma1,2.5)*cpow(lambda2,2.5)*
      k2*k2*cpow(gamma1,5.)/(2.*cpow(base,2.5));
  f2=cpow(gamma1,4.5)*cpow(lambda2,4.5)/(2.*cpow(base,1.5));
  denominator=k2*k2*gamma1;
  dv=1./denominator;
  dv1=-lambda*k1*k1/(lambda2*denominator*denominator);
  dv2=k1*k1*k2*k2/(lambda2*lambda2*lambda2*
                    denominator*denominator);
  du1=-lambda/(gamma1*gamma1*lambda2);
  du2=k2*k2/(gamma1*gamma1*lambda2*lambda2*lambda2)+
      2.*lambda*lambda/(gamma1*gamma1*gamma1*lambda2*lambda2);

  root=csqrt(lambda); l32=lambda*root; l52=lambda*lambda*root;
  rv=f1/lambda2*(l52-3.*I*l32/(8.*rho))*dv-
     f2/(lambda2*lambda2*lambda2)*
       (lambda*l52-3.*I*l52/(8.*rho))*dv;
  zv=(l52+I*l32/(8.*rho))*(f1*dv1+f2*dv2);
  rh=(l52-7.*I*l32/(8.*rho))*(f1*dv1+f2*dv2);
  ph=(l32-3.*I*root/(8.*rho))*(f1*dv1+f2*dv2);
  zh=gamma1*ph;
  u=(root+I/(8.*rho*root))*(f1*du1+f2*du2);

  common=(2.*M_PI*frequency_hz)*MU0/(4.*M_PI);
  phase=cexp(-gamma1*depth-I*lambda*rho);
  scale=2./sqrt(rho)*phase;
  /* Principal fractional powers are opposite to the continuation along the
     outgoing steepest-descent path at this branch point. */
  result->value[NEC4_GROUND_RHO_VERTICAL]=
    -common*plus*scale*rv;
  result->value[NEC4_GROUND_Z_VERTICAL]=
    -common*minus*scale*zv;
  result->value[NEC4_GROUND_RHO_HORIZONTAL]=
    -I*common*plus*scale*(rh-u);
  result->value[NEC4_GROUND_PHI_HORIZONTAL]=
    -I*common*scale*(minus*ph/rho+plus*u);
  result->value[NEC4_GROUND_Z_HORIZONTAL]=
    -common*plus*scale*zh;
  result->saddle=k2;
  return 1;
}

static int second_axis(double frequency,double complex k1,double complex k2,
    double z,double depth,nec4_ground_asymptotic_result *result)
{
  double complex f0,f2,f4,gv0,gv2,gu0,gu2,iv,u;
  double complex q;
  if(!result || z<0. || depth<0. || z+depth<=0. || cabs(k1)==0. ||
     cabs(k2)==0. || cabs(k1+k2)==0.) return 0;
  f0=I*(k2*z+k1*depth);
  f2=-I*(z/k2+depth/k1);
  f4=-3.*I*(z/(k2*k2*k2)+depth/(k1*k1*k1));
  if(cabs(f2)==0.) return 0;
  gv0=-I/(k1*k2*(k1+k2));
  gv2=-I*(1./(k1*k1*k1)+1./(k2*k2*k2))/((k1+k2)*(k1+k2));
  gu0=-2.*I/(k1+k2);
  gu2=-2.*I*(1./k1+1./k2)/((k1+k2)*(k1+k2));
  iv=2./(f2*f2)*cexp(-f0)*(gv0+(2.*gv2*f2-gv0*f4)/(f2*f2));
  u=1./f2*cexp(-f0)*(gu0+(gu2*f2-gu0*f4/3.)/(f2*f2));
  q=I*(2.*M_PI*frequency)*MU0/(4.*M_PI);
  result->value[0]=0.;
  result->value[1]=-q*2.*iv;
  result->value[2]=-q*(-iv+u);
  result->value[3]=q*(-iv+u);
  result->value[4]=0.;
  result->saddle=0.;
  return 1;
}

static double complex saddle_exponential(double complex k1,double complex k2,
    double rho,double z,double depth,double complex saddle)
{
  return cexp(-I*outgoing_beta(k2,saddle)*z-
              I*outgoing_beta(k1,saddle)*depth-I*saddle*rho);
}

int nec4_ground_asymptotic_second_saddle(double frequency_hz,
    double complex k1,double complex k2,double rho,double z,double depth,
    nec4_ground_asymptotic_result *result)
{
  double rho_a;
  if(!result || frequency_hz<=0. || rho<0. || z<0. || depth<0.) return 0;
  /* Very near a source on the interface the generic double-precision
     derivatives alternate across the k2 branch point.  The exact z=0
     Appendix-J limit differs by less than 0.1 percent at this switch in the
     qualified high-loss case and remains continuous with the outer form. */
  if(rho>0. && ((z==0.) ||
      (depth==0. && z<=3.e-4*(2.*M_PI/cabs(k2)))))
    return second_interface_raw(frequency_hz,k1,k2,rho,depth,result);
  rho_a=.02*(z+depth);
  if(rho_a<=0.) return 0;
  if(rho==0.) return second_axis(frequency_hz,k1,k2,z,depth,result);
  if(rho>=rho_a) return second_saddle_raw(frequency_hz,k1,k2,rho,z,depth,result);
  {
    nec4_ground_asymptotic_result first,first_a,second_a,axis;
    double complex phase,phase_a,phase_0;
    unsigned int component;
    if(!nec4_ground_asymptotic_first(frequency_hz,k1,k2,rho,z,depth,&first) ||
       !nec4_ground_asymptotic_first(frequency_hz,k1,k2,rho_a,z,depth,&first_a) ||
       !second_saddle_raw(frequency_hz,k1,k2,rho_a,z,depth,&second_a) ||
       !second_axis(frequency_hz,k1,k2,z,depth,&axis)) return 0;
    phase=saddle_exponential(k1,k2,rho,z,depth,first.saddle);
    phase_a=saddle_exponential(k1,k2,rho_a,z,depth,first_a.saddle);
    phase_0=cexp(-I*(k2*z+k1*depth));
    if(cabs(phase_a)==0. || cabs(phase_0)==0.) return 0;
    for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;component++) {
      double complex higher_a=(second_a.value[component]-first_a.value[component])/
                              phase_a;
      double complex higher_0;
      double t=rho/rho_a;
      if(component==NEC4_GROUND_RHO_VERTICAL ||
         component==NEC4_GROUND_Z_HORIZONTAL)
        higher_0=0.;
      else {
        /* J-19 is the higher-order axis correction, not the complete axis
           field.  Only E_phi^H has a nonzero first-order axis limit. */
        higher_0=axis.value[component]/phase_0;
        if(component==NEC4_GROUND_PHI_HORIZONTAL) {
          double complex d=z/k2+depth/k1;
          double complex first_phi=I*(2.*M_PI*frequency_hz)*MU0/(4.*M_PI)*
            (2.*I/(k1*k2*(k1+k2)*d*d))*phase_0;
          higher_0=(axis.value[component]-first_phi)/phase_0;
        }
      }
      if(component==NEC4_GROUND_RHO_VERTICAL ||
         component==NEC4_GROUND_Z_HORIZONTAL)
        result->value[component]=first.value[component]+phase*t*higher_a;
      else
        result->value[component]=first.value[component]+phase*
          (higher_0+t*t*(higher_a-higher_0));
    }
    result->saddle=first.saddle;
    return 1;
  }
}

/* Q(z)=-sqrt(pi) z exp(z^2) erfc(z).  The pole mapping used here keeps z in
   the right half-plane.  A convergent erf series covers the coalescence
   region; farther out the scaled erfc asymptotic avoids overflow and loss of
   significance.  Long-double intermediates leave double-roundoff headroom. */
static double complex modified_q(double complex input)
{
  long double complex z=(long double complex)input;
  long double magnitude=cabsl(z);
  if(magnitude<4.L) {
    long double complex term=z,sum=z;
    unsigned int n;
    for(n=0;n<400;n++) {
      long double complex next=term*(-z*z)*(2.L*n+1.L)/
                               ((n+1.L)*(2.L*n+3.L));
      sum+=next; term=next;
      if(cabsl(next)<=8.L*LDBL_EPSILON*fmaxl(1.L,cabsl(sum))) break;
    }
    {
      long double complex erfc_value=1.L-2.L/sqrtl(M_PI)*sum;
      return (double complex)(-sqrtl(M_PI)*z*cexpl(z*z)*erfc_value);
    }
  } else {
    long double complex term=1.L,sum=1.L;
    long double previous=cabsl(term);
    unsigned int n;
    for(n=1;n<200;n++) {
      long double complex next=term*-(2.L*n-1.L)/(2.L*z*z);
      long double current=cabsl(next);
      if(current>previous) break;
      sum+=next; term=next; previous=current;
      if(current<=8.L*LDBL_EPSILON*fmaxl(1.L,cabsl(sum))) break;
    }
    return (double complex)(-sum);
  }
}

typedef struct {
  double complex pole,gamma1,gamma2,a,b,saddle_phase,saddle;
} pole_geometry;

static int build_pole_geometry(double complex k1,double complex k2,double rho,
    double z,double depth,pole_geometry *geometry)
{
  double complex sum_squared,denominator_derivative,difference;
  if(!geometry || rho<=0. ||
     !saddle_point(k1,k2,rho,z,depth,&geometry->saddle)) return 0;
  sum_squared=k1*k1+k2*k2;
  geometry->pole=k1*k2/csqrt(sum_squared);
  geometry->gamma1=I*outgoing_beta(k1,geometry->pole);
  geometry->gamma2=I*outgoing_beta(k2,geometry->pole);
  geometry->saddle_phase=I*outgoing_beta(k2,geometry->saddle)*z+
                         I*outgoing_beta(k1,geometry->saddle)*depth+
                         I*geometry->saddle*rho;
  difference=(geometry->gamma2*z+geometry->gamma1*depth+
              I*geometry->pole*rho)-geometry->saddle_phase;
  geometry->b=csqrt(difference);
  if(cimag(geometry->b)>0. ||
     (fabs(cimag(geometry->b))<1.e-14*fmax(1.,cabs(geometry->b)) &&
      creal(geometry->b)<0.)) geometry->b=-geometry->b;
  if(cabs(geometry->b)==0.) return 0;
  denominator_derivative=k1*k1*geometry->pole/geometry->gamma2+
                         k2*k2*geometry->pole/geometry->gamma1;
  if(cabs(denominator_derivative)==0.) return 0;
  geometry->a=csqrt(2.*I/(M_PI*rho))*csqrt(geometry->pole)/
              denominator_derivative;
  return 1;
}

static int pole_potential(double complex k1,double complex k2,double rho,
    double z,double depth,int second,pole_geometry *geometry,
    double complex *potential)
{
  double complex bracket;
  if(!potential || !build_pole_geometry(k1,k2,rho,z,depth,geometry)) return 0;
  bracket=modified_q(I*geometry->b)+1.;
  if(second) bracket+=1./(2.*geometry->b*geometry->b);
  *potential=sqrt(M_PI)*geometry->a/geometry->b*
             cexp(-geometry->saddle_phase)*bracket;
  return isfinite(creal(*potential)) && isfinite(cimag(*potential));
}

int nec4_ground_asymptotic_first_pole(double frequency_hz,
    double complex k1,double complex k2,double rho,double z,double depth,
    nec4_ground_asymptotic_result *result)
{
  pole_geometry geometry;
  double complex potential,sum_squared;
  double common;
  if(!result || frequency_hz<=0. || z<0. || depth<0. ||
     !pole_potential(k1,k2,rho,z,depth,0,&geometry,&potential)) return 0;
  sum_squared=k1*k1+k2*k2;
  common=(2.*M_PI*frequency_hz)*MU0/(4.*M_PI);
  result->value[0]=-I*common*k1*k2*k2*k2/sum_squared*potential;
  result->value[1]=I*common*k1*k1*k2*k2/sum_squared*potential;
  result->value[2]=result->value[1];
  result->value[3]=common/rho*k1*k2/csqrt(sum_squared)*potential;
  result->value[4]=-I*common*k1*k1*k1*k2/sum_squared*potential;
  result->saddle=geometry.saddle;
  return 1;
}

int nec4_ground_asymptotic_second_pole(double frequency_hz,
    double complex k1,double complex k2,double rho,double z,double depth,
    nec4_ground_asymptotic_result *result)
{
  pole_geometry geometry;
  double complex potential,radial_first;
  double common;
  if(!result || frequency_hz<=0. || z<0. || depth<0. ||
     !pole_potential(k1,k2,rho,z,depth,1,&geometry,&potential)) return 0;
  common=(2.*M_PI*frequency_hz)*MU0/(4.*M_PI);
  radial_first=geometry.pole-3.*I/(8.*rho);
  result->value[0]=common*radial_first*geometry.gamma2*potential;
  result->value[1]=-I*common*(geometry.pole*geometry.pole+
                             I*geometry.pole/(8.*rho))*potential;
  result->value[2]=I*common*(geometry.pole*geometry.pole-
                            7.*I*geometry.pole/(8.*rho))*potential;
  result->value[3]=common/rho*radial_first*potential;
  result->value[4]=common*radial_first*geometry.gamma1*potential;
  result->saddle=geometry.saddle;
  return 1;
}

int nec4_ground_asymptotic_first_uniform(double frequency_hz,
    double complex k1,double complex k2,double rho,double z,double depth,
    nec4_ground_asymptotic_result *result)
{
  nec4_ground_asymptotic_result saddle,pole;
  unsigned int component;
  if(!result || !nec4_ground_asymptotic_first(frequency_hz,k1,k2,rho,z,depth,
                                               &saddle) ||
     !nec4_ground_asymptotic_first_pole(frequency_hz,k1,k2,rho,z,depth,&pole))
    return 0;
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;component++)
    result->value[component]=saddle.value[component]+pole.value[component];
  result->saddle=saddle.saddle;
  return 1;
}

int nec4_ground_asymptotic_second_uniform(double frequency_hz,
    double complex k1,double complex k2,double rho,double z,double depth,
    nec4_ground_asymptotic_result *result)
{
  nec4_ground_asymptotic_result saddle,pole;
  unsigned int component;
  if(!result || !nec4_ground_asymptotic_second_saddle(frequency_hz,k1,k2,rho,
                                                       z,depth,&saddle)) return 0;
  if(rho==0.) { *result=saddle; return 1; }
  if(!nec4_ground_asymptotic_second_pole(frequency_hz,k1,k2,rho,z,depth,&pole))
    return 0;
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;component++)
    result->value[component]=saddle.value[component]+pole.value[component];
  result->saddle=saddle.saddle;
  return 1;
}
