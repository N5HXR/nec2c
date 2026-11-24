#include "junction_charge.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef JUNCTION_CHARGE_PURE
#include "nec2c.h"
#include "shared.h"
#endif

/* Positive abscissas and weights for 32-point Gauss-Legendre quadrature. */
static const double gauss_x[16]={
  .0483076656877383162,.144471961582796493,.239287362252137075,
  .331868602282127650,.421351276130635345,.506899908932229390,
  .587715757240762329,.663044266930215201,.732182118740289680,
  .794483795967942407,.849367613732569970,.896321155766052124,
  .934906075937739689,.964762255587506431,.985611511545268335,
  .997263861849481564};
static const double gauss_w[16]={
  .0965400885147278006,.0956387200792748594,.0938443990808045656,
  .0911738786957638847,.0876520930044038111,.0833119242269467552,
  .0781938957870703065,.0723457941088485062,.0658222227763618468,
  .0586840934785355471,.0509980592623761762,.0428358980222266807,
  .0342738629130214331,.0253920653092620595,.0162743947309056706,
  .00701861000947009660};

#ifndef JUNCTION_CHARGE_PURE
static double complex *end_factor;
static unsigned char *end_assigned;
static size_t endpoint_capacity;
static junction_charge_diagnostics last_diagnostics;
#endif

static double dot3(const double a[3],const double b[3])
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

static double distance_to_source(const double observer[3],
                                 const double source_direction[3],
                                 double source_s,double source_radius)
{
  double dx=observer[0]-source_s*source_direction[0];
  double dy=observer[1]-source_s*source_direction[1];
  double dz=observer[2]-source_s*source_direction[2];
  return sqrt(dx*dx+dy*dy+dz*dz+source_radius*source_radius);
}

static double complex green(double distance)
{
  double wave_number=2.*acos(-1.);
  return cexp(-I*wave_number*distance)/distance;
}

static double complex triangle_panel(const double observer[3],double length,
                                     double radius,const double direction[3],
                                     double lower,double upper)
{
  double center=.5*(lower+upper),half=.5*(upper-lower);
  double complex sum=0.;
  int q;
  for(q=0;q<16;++q)
  {
    double offset=half*gauss_x[q];
    double s1=center-offset,s2=center+offset;
    sum += gauss_w[q]*((length-s1)*green(distance_to_source(observer,direction,s1,radius))+
                       (length-s2)*green(distance_to_source(observer,direction,s2,radius)));
  }
  return half*sum;
}

static double complex triangle_adaptive(const double observer[3],double length,
                                        double radius,const double direction[3],
                                        double lower,double upper,
                                        double complex coarse,int depth)
{
  double middle=.5*(lower+upper);
  double complex left=triangle_panel(observer,length,radius,direction,lower,middle);
  double complex right=triangle_panel(observer,length,radius,direction,middle,upper);
  double complex fine=left+right;
  if(depth>=12 || cabs(fine-coarse)<=1.e-13*fmax(length,cabs(fine))) return fine;
  return triangle_adaptive(observer,length,radius,direction,lower,middle,left,depth+1)+
    triangle_adaptive(observer,length,radius,direction,middle,upper,right,depth+1);
}

static double complex triangle_potential(const double observer[3],double length,
                                         double radius,const double direction[3])
{
  double complex coarse=triangle_panel(observer,length,radius,direction,0.,length);
  return triangle_adaptive(observer,length,radius,direction,0.,length,coarse,0);
}

static double complex line_panel(const double observer[3],double radius,
                                 const double direction[3],double lower,
                                 double upper)
{
  double center=.5*(lower+upper),half=.5*(upper-lower);
  double complex sum=0.;
  int q;
  for(q=0;q<16;++q)
  {
    double offset=half*gauss_x[q];
    sum += gauss_w[q]*(green(distance_to_source(observer,direction,center-offset,radius))+
                       green(distance_to_source(observer,direction,center+offset,radius)));
  }
  return half*sum;
}

static double complex line_adaptive(const double observer[3],double radius,
                                    const double direction[3],double lower,
                                    double upper,double complex coarse,int depth)
{
  double middle=.5*(lower+upper);
  double complex left=line_panel(observer,radius,direction,lower,middle);
  double complex right=line_panel(observer,radius,direction,middle,upper);
  double complex fine=left+right;
  if(depth>=12 || cabs(fine-coarse)<=1.e-13*fmax(fabs(upper-lower),cabs(fine))) return fine;
  return line_adaptive(observer,radius,direction,lower,middle,left,depth+1)+
    line_adaptive(observer,radius,direction,middle,upper,right,depth+1);
}

/* Appendix A equation A-5: finite integral plus a half infinite line. */
static double complex tail_potential(const double observer[3],double radius,
                                     const double direction[3])
{
  double wave_number=2.*acos(-1.);
  double axial=dot3(observer,direction);
  double rho=sqrt(fmax(dot3(observer,observer)+radius*radius-axial*axial,DBL_MIN));
  double complex finite=0.,hankel=j0(wave_number*rho)-I*y0(wave_number*rho);
  if(axial!=0.)
  {
    double complex coarse=line_panel(observer,radius,direction,0.,axial);
    finite=line_adaptive(observer,radius,direction,0.,axial,coarse,0);
  }
  return finite-I*.5*acos(-1.)*hankel;
}

static int solve_dense(size_t order,double complex *a,double complex *b,
                       double *pivot_ratio)
{
  size_t i,j,k,pivot;
  double largest=0.,smallest=DBL_MAX;
  for(k=0;k<order;++k)
  {
    double magnitude;
    pivot=k;magnitude=cabs(a[k*order+k]);
    for(i=k+1;i<order;++i)
      if(cabs(a[i*order+k])>magnitude)
      { magnitude=cabs(a[i*order+k]);pivot=i; }
    if(!isfinite(magnitude) || magnitude<=DBL_EPSILON) return 0;
    if(pivot!=k)
    {
      for(j=k;j<order;++j)
      { double complex temporary=a[k*order+j];a[k*order+j]=a[pivot*order+j];a[pivot*order+j]=temporary; }
      { double complex temporary=b[k];b[k]=b[pivot];b[pivot]=temporary; }
    }
    largest=fmax(largest,magnitude);smallest=fmin(smallest,magnitude);
    for(i=k+1;i<order;++i)
    {
      double complex multiplier=a[i*order+k]/a[k*order+k];
      a[i*order+k]=0.;
      for(j=k+1;j<order;++j) a[i*order+j]-=multiplier*a[k*order+j];
      b[i]-=multiplier*b[k];
    }
  }
  for(i=order;i-- >0;)
  {
    double complex value=b[i];
    for(j=i+1;j<order;++j) value-=a[i*order+j]*b[j];
    b[i]=value/a[i*order+i];
    if(!isfinite(creal(b[i])) || !isfinite(cimag(b[i]))) return 0;
  }
  *pivot_ratio=smallest>0.?largest/smallest:INFINITY;
  return 1;
}

int junction_charge_solve(size_t count,const double *length,
                          const double *radius,const double (*direction)[3],
                          double complex *factor,double *scaled_residual,
                          double *pivot_ratio)
{
  size_t order=2*count,row,column,i,k,l;
  double complex *matrix,*original,*rhs,*original_rhs,mean=0.;
  double residual=0.;
  if(count<2 || order>SIZE_MAX/order/sizeof(double complex)) return 0;
  matrix=calloc(order*order,sizeof(*matrix));
  original=calloc(order*order,sizeof(*original));
  rhs=calloc(order,sizeof(*rhs));
  original_rhs=calloc(order,sizeof(*original_rhs));
  if(!matrix || !original || !rhs || !original_rhs)
  { free(matrix);free(original);free(rhs);free(original_rhs);return 0; }

  row=0;
  for(i=0;i<count;++i)
    for(k=0;k<2;++k)
    {
      double observer[3],s=(k==0?.5:1.)*length[i];
      observer[0]=s*direction[i][0];observer[1]=s*direction[i][1];observer[2]=s*direction[i][2];
      for(l=0;l<count;++l)
      {
        column=2*l;
        matrix[row*order+column]=triangle_potential(observer,length[l],radius[l],direction[l]);
        matrix[row*order+column+1]=length[l]*
          tail_potential(observer,radius[l],direction[l]);
      }
      rhs[row]=1.;
      ++row;
    }

  memcpy(original,matrix,order*order*sizeof(*matrix));
  memcpy(original_rhs,rhs,order*sizeof(*rhs));
  if(!solve_dense(order,matrix,rhs,pivot_ratio))
  { free(matrix);free(original);free(rhs);free(original_rhs);return 0; }
  for(i=0;i<count;++i)
  { factor[i]=(rhs[2*i]+rhs[2*i+1])*length[i];mean+=factor[i]; }
  mean/=count;
  if(!isfinite(creal(mean)) || !isfinite(cimag(mean)) || cabs(mean)<=DBL_MIN)
  { free(matrix);free(original);free(rhs);free(original_rhs);return 0; }
  for(i=0;i<count;++i)
  {
    factor[i]/=mean;
    if(!isfinite(creal(factor[i])) || !isfinite(cimag(factor[i])) || creal(factor[i])<=0.)
    { free(matrix);free(original);free(rhs);free(original_rhs);return 0; }
  }
  for(i=0;i<order;++i)
  {
    double complex value=-original_rhs[i];
    double scale=cabs(original_rhs[i]);
    for(k=0;k<order;++k)
    { value+=original[i*order+k]*rhs[k];scale+=cabs(original[i*order+k]*rhs[k]); }
    residual=fmax(residual,cabs(value)/fmax(1.,scale));
  }
  *scaled_residual=residual;
  free(matrix);free(original);free(rhs);free(original_rhs);
  return 1;
}

#ifndef JUNCTION_CHARGE_PURE
static size_t find_root(size_t *parent,size_t item)
{
  size_t root=item;
  while(parent[root]!=root) root=parent[root];
  while(parent[item]!=item)
  { size_t next=parent[item];parent[item]=root;item=next; }
  return root;
}

static void join(size_t *parent,size_t a,size_t b)
{
  a=find_root(parent,a);b=find_root(parent,b);
  if(a!=b) parent[b]=a;
}

int junction_charge_prepare(FILE *report,int strict,
                            const unsigned char *prescribed_endpoints)
{
  size_t endpoints=2*(size_t)data.n,*parent,*counts,*offsets,*cursor,*members;
  size_t i,end,root,group;
  memset(&last_diagnostics,0,sizeof(last_diagnostics));
  if(endpoints>endpoint_capacity)
  {
    double complex *new_factor=realloc(end_factor,endpoints*sizeof(*end_factor));
    unsigned char *new_assigned=realloc(end_assigned,endpoints*sizeof(*end_assigned));
    if(!new_factor || !new_assigned) abort_on_error(-1);
    end_factor=new_factor;end_assigned=new_assigned;endpoint_capacity=endpoints;
  }
  memset(end_assigned,0,endpoints*sizeof(*end_assigned));
  if(!data.n) return 1;
  parent=malloc(endpoints*sizeof(*parent));counts=calloc(endpoints,sizeof(*counts));
  offsets=calloc(endpoints+1,sizeof(*offsets));cursor=calloc(endpoints,sizeof(*cursor));
  members=malloc(endpoints*sizeof(*members));
  if(!parent||!counts||!offsets||!cursor||!members) abort_on_error(-1);
  for(i=0;i<endpoints;++i) parent[i]=i;
  for(i=0;i<(size_t)data.n;++i)
    for(end=0;end<2;++end)
    {
      int icon=end?data.icon2[i]:data.icon1[i];
      size_t other,other_end;
      if(!icon || icon>PCHCON || icon==-(int)(i+1) || icon==(int)(i+1)) continue;
      other=(size_t)(abs(icon)-1);
      other_end=(size_t)(((icon>0)?1:0) ^ (int)end);
      if(other<(size_t)data.n) join(parent,2*i+end,2*other+other_end);
    }
  for(i=0;i<endpoints;++i) counts[find_root(parent,i)]++;
  for(i=0;i<endpoints;++i) offsets[i+1]=offsets[i]+counts[i];
  memcpy(cursor,offsets,endpoints*sizeof(*cursor));
  for(i=0;i<endpoints;++i) { root=find_root(parent,i);members[cursor[root]++]=i; }

  for(root=0;root<endpoints;++root)
  {
    size_t count=counts[root];
    double *length,*radius,(*direction)[3],residual=0.,pivot_ratio=1.;
    double complex *factor;
    int equal_radii=1,all_prescribed=prescribed_endpoints!=NULL;
    int success;
    if(count<2) continue;
    ++last_diagnostics.junctions;
    if(count>last_diagnostics.maximum_order) last_diagnostics.maximum_order=count;
    length=malloc(count*sizeof(*length));radius=malloc(count*sizeof(*radius));
    factor=malloc(count*sizeof(*factor));direction=malloc(count*sizeof(*direction));
    if(!length||!radius||!factor||!direction)
      abort_on_error(-1);
    for(group=0;group<count;++group)
    {
      size_t endpoint_index=members[offsets[root]+group];
      size_t segment=endpoint_index/2;
      int which=(int)(endpoint_index&1);
      length[group]=data.si[segment];radius[group]=data.bi[segment];
      direction[group][0]=(which?-1.:1.)*data.cab[segment];
      direction[group][1]=(which?-1.:1.)*data.sab[segment];
      direction[group][2]=(which?-1.:1.)*data.salp[segment];
      if(group && fabs(radius[group]/radius[0]-1.)>1.e-12) equal_radii=0;
      if(!prescribed_endpoints || !prescribed_endpoints[endpoint_index])
        all_prescribed=0;
    }
    if(all_prescribed)
    {
      /* The NEC4 interface condition owns these endpoint ratios.  Solving a
         homogeneous local junction problem here would both be redundant and
         ill-conditioned for a dense regular radial screen. */
      success=1;
      ++last_diagnostics.prescribed;
    }
    else if(count==2 && equal_radii)
    { factor[0]=factor[1]=1.;success=1;++last_diagnostics.symmetric; }
    else
    {
      success=junction_charge_solve(count,length,radius,direction,factor,&residual,&pivot_ratio);
      if(success) ++last_diagnostics.solved;
    }
    if(all_prescribed)
    {
      /* Ratios remain owned by nec4_medium_state. */
    }
    else if(success)
    {
      for(group=0;group<count;++group)
      { size_t endpoint=members[offsets[root]+group];end_factor[endpoint]=factor[group];end_assigned[endpoint]=1; }
      last_diagnostics.maximum_scaled_residual=fmax(last_diagnostics.maximum_scaled_residual,residual);
      last_diagnostics.maximum_pivot_ratio=fmax(last_diagnostics.maximum_pivot_ratio,pivot_ratio);
    }
    else ++last_diagnostics.fallback;
    free(length);free(radius);free(factor);free(direction);
  }
  if(report && last_diagnostics.junctions)
    fprintf(report,"\n NEC4ISH JUNCTION CHARGE: %zu JUNCTIONS, %zu LOCAL SOLVES, %zu SYMMETRIC, %zu INTERFACE-PRESCRIBED, %zu FALLBACKS"
      "\n                       MAX ORDER %zu, RESIDUAL %.3E, PIVOT RATIO %.3E\n",
      last_diagnostics.junctions,last_diagnostics.solved,last_diagnostics.symmetric,
      last_diagnostics.prescribed,
      last_diagnostics.fallback,last_diagnostics.maximum_order,
      last_diagnostics.maximum_scaled_residual,last_diagnostics.maximum_pivot_ratio);
  free(parent);free(counts);free(offsets);free(cursor);free(members);
  return !(strict && last_diagnostics.fallback);
}

double complex junction_charge_factor(size_t segment,int end,double complex fallback)
{
  size_t endpoint=2*segment+(end?1:0);
  return endpoint<endpoint_capacity && end_assigned[endpoint]?end_factor[endpoint]:fallback;
}

const junction_charge_diagnostics *junction_charge_last_diagnostics(void)
{ return &last_diagnostics; }
#endif
