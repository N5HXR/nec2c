#include "vie_halfspace.h"
#include <float.h>
#include <math.h>
#include <string.h>

#define PI 3.1415926535897932384626433832795
#define MAX_ORDER 256

static int material_valid(complex double e)
{ return isfinite(creal(e)) && isfinite(cimag(e)) && creal(e)>=1 && cimag(e)<=0; }

void vie_halfspace_defaults(vie_halfspace_options *o)
{ o->relative_tolerance=2e-7; o->absolute_tolerance=1e-10; o->maximum_order=128; }

int vie_halfspace_fresnel(complex double e, double u,
                          complex double *rs, complex double *rp)
{
  complex double w,v,a,b;
  if(!material_valid(e) || !isfinite(u) || u<0 || !rs || !rp) return 0;
  if(e==1.0) { *rs=0; *rp=0; return 1; }
  w=u<=1 ? sqrt((1-u)*(1+u)):-I*sqrt((u-1)*(u+1));
  v=csqrt(e-u*u);
  if(cimag(v)>0) v=-v; /* outgoing/decaying branch, including lossless evanescence */
  a=w+v; b=e*w+v;
  if(cabs(a)==0 || cabs(b)==0) return 0;
  *rs=(w-v)/a; *rp=(e*w-v)/b;
  return isfinite(creal(*rs)) && isfinite(cimag(*rs)) &&
         isfinite(creal(*rp)) && isfinite(cimag(*rp));
}

static int gauss(int count, double *x, double *weight)
{
  int j;
  for(j=0;j<(count+1)/2;++j)
  {
    double z=cos(PI*(j+.75)/(count+.5)),derivative=0;
    int it;
    for(it=0;it<40;++it)
    {
      double p=1,q=z,previous=z;
      int m;
      for(m=2;m<=count;++m) { double r=((2*m-1)*z*q-(m-1)*p)/m; p=q; q=r; }
      derivative=count*(z*q-p)/(z*z-1);
      z-=q/derivative;
      if(fabs(z-previous)<4*DBL_EPSILON) break;
    }
    if(it==40) return 0;
    x[j]=-z; x[count-1-j]=z;
    weight[j]=weight[count-1-j]=2/((1-z*z)*derivative*derivative);
  }
  return 1;
}

static double sinc_real(double v)
{ return fabs(v)<1e-5 ? 1-v*v/6+v*v*v*v/120:sin(v)/v; }

/* Stable even for very large evanescent kz*cell height: no sinh overflow. */
static complex double z_average(double z, double size, complex double kz)
{
  complex double a=kz*size;
  if(cabs(a)<1e-4)
    return cexp(-I*kz*z)*(1-a*a/24+a*a*a*a/1920);
  return (cexp(-I*kz*(z-.5*size))-cexp(-I*kz*(z+.5*size)))/(I*a);
}

static int quadrature(const vie_cuboid *a, const vie_cuboid *b, double volume,
                       double k, complex double e, int order,
                       complex double n[3][3], complex double h[3][3])
{
  double nodes[MAX_ORDER],weights[MAX_ORDER];
  double gap=a->center[2]+b->center[2]-.5*(a->size[2]+b->size[2]);
  double scale=k*gap;
  /* Integrate v=k*gap*sqrt(u^2-1); exp(-v) controls the evanescent tail.
     Split at the lossless substrate branch and fixed decay scales. */
  double cut[7]={0,1,4,16,64,0,0};
  int cuts=5,segment,j,az,p,q;
  double branch=scale*sqrt(creal(e)-1);
  if(branch>0 && branch<64) cut[cuts++]=branch;
  for(j=1;j<cuts;++j)
  { double v=cut[j]; int m=j; while(m>0 && cut[m-1]>v) { cut[m]=cut[m-1]; --m; } cut[m]=v; }
  memset(n,0,9*sizeof(complex double)); memset(h,0,9*sizeof(complex double));
  if(!isfinite(scale) || !(scale>0) || !gauss(order,nodes,weights)) return 0;
  for(segment=-1;segment<cuts-1;++segment)
  {
    double lo=segment<0 ? 0:cut[segment];
    double hi=segment<0 ? PI/2:cut[segment+1];
    if(hi==lo) continue;
    for(j=0;j<order;++j)
    {
      double t=.5*((hi-lo)*nodes[j]+hi+lo),u;
      complex double w,radial,rs,rp,zfactor;
      double rw=.5*(hi-lo)*weights[j];
      if(segment<0) { u=sin(t); w=cos(t); radial=sin(t)*rw; }
      else { double ev=t/scale; u=hypot(1,ev); w=-I*ev; radial=I*rw/scale; }
      if(!vie_halfspace_fresnel(e,u,&rs,&rp)) return 0;
      zfactor=z_average(a->center[2],a->size[2],k*w)*
              z_average(b->center[2],b->size[2],k*w);
      for(az=0;az<2*order;++az)
      {
        double phi=PI*az/order,c=cos(phi),s=sin(phi);
        double lateral=k*u*(c*(a->center[0]-b->center[0])+s*(a->center[1]-b->center[1]));
        double shape=sinc_real(.5*k*u*c*a->size[0])*sinc_real(.5*k*u*s*a->size[1])*
                     sinc_real(.5*k*u*c*b->size[0])*sinc_real(.5*k*u*s*b->size[1]);
        complex double te[3]={-s,c,0},up[3]={w*c,w*s,-u},down[3]={-w*c,-w*s,-u};
        complex double kv[3]={u*c,u*s,w};
        complex double weight=-I*k*k*k/(8*PI*PI)*volume*(PI/order)*radial*
                               zfactor*shape*cexp(-I*lateral);
        complex double block[3][3];
        for(p=0;p<3;++p) for(q=0;q<3;++q)
        { block[p][q]=rs*te[p]*te[q]+rp*up[p]*down[q]; n[p][q]+=weight*block[p][q]; }
        for(p=0;p<3;++p) for(q=0;q<3;++q)
          h[p][q]+=weight*(kv[(p+1)%3]*block[(p+2)%3][q]-
                          kv[(p+2)%3]*block[(p+1)%3][q])/(I*k);
      }
    }
  }
  for(p=0;p<3;++p) for(q=0;q<3;++q)
    if(!isfinite(creal(n[p][q])) || !isfinite(cimag(n[p][q])) ||
       !isfinite(creal(h[p][q])) || !isfinite(cimag(h[p][q]))) return 0;
  return 1;
}

static int evaluate(const vie_cuboid *a,const vie_cuboid *b,double volume,
                     double k,complex double e,const vie_halfspace_options *requested,
                     complex double n[3][3],complex double h[3][3],vie_halfspace_report *report)
{
  vie_halfspace_options defaults;
  const vie_halfspace_options *o=requested;
  complex double oldn[3][3],oldh[3][3],newn[3][3],newh[3][3];
  int order,p,q;
  if(!isfinite(k) || !(k>0) || !material_valid(e) || !n || !h) return 0;
  vie_halfspace_defaults(&defaults); if(!o) o=&defaults;
  if(!isfinite(o->relative_tolerance) || !(o->relative_tolerance>0) ||
     !isfinite(o->absolute_tolerance) || !(o->absolute_tolerance>0) ||
     o->maximum_order<16 || o->maximum_order>MAX_ORDER) return 0;
  if(report) { report->order=0; report->electric_delta=INFINITY; report->magnetic_delta=INFINITY; }
  if(e==1.0)
  { memset(n,0,sizeof(newn)); memset(h,0,sizeof(newh));
    if(report) { report->electric_delta=0; report->magnetic_delta=0; } return 1; }
  for(order=8;order<=o->maximum_order;order*=2)
  {
    double dn=0,dh=0,sn=0,sh=0;
    if(!quadrature(a,b,volume,k,e,order,newn,newh)) return 0;
    if(order>8)
    {
      for(p=0;p<3;++p) for(q=0;q<3;++q)
      { dn=fmax(dn,cabs(newn[p][q]-oldn[p][q])); dh=fmax(dh,cabs(newh[p][q]-oldh[p][q]));
        sn=fmax(sn,cabs(newn[p][q])); sh=fmax(sh,cabs(newh[p][q])); }
      if(report) { report->order=order; report->electric_delta=dn; report->magnetic_delta=dh; }
      if(dn<=o->absolute_tolerance+o->relative_tolerance*sn &&
         dh<=o->absolute_tolerance+o->relative_tolerance*sh)
      { memcpy(n,newn,sizeof(newn)); memcpy(h,newh,sizeof(newh)); return 1; }
    }
    memcpy(oldn,newn,sizeof(newn)); memcpy(oldh,newh,sizeof(newh));
  }
  return 0;
}

static int point_valid(const double p[3])
{ return p && isfinite(p[0]) && isfinite(p[1]) && isfinite(p[2]) && p[2]>0; }

static int cuboid_valid(const vie_cuboid *a)
{
  int j;
  if(!a) return 0;
  for(j=0;j<3;++j)
    if(!isfinite(a->center[j]) || !isfinite(a->size[j]) || !(a->size[j]>0)) return 0;
  return a->center[2]-.5*a->size[2]>0;
}

int vie_halfspace_point(const double t[3],const double s[3],double k,complex double e,
                        const vie_halfspace_options *o,complex double n[3][3],
                        complex double h[3][3],vie_halfspace_report *r)
{
  vie_cuboid a={{0,0,0},{0,0,0}},b=a;
  int j;
  if(!point_valid(t) || !point_valid(s)) return 0;
  for(j=0;j<3;++j) { a.center[j]=t[j]; b.center[j]=s[j]; }
  return evaluate(&a,&b,1,k,e,o,n,h,r);
}

int vie_halfspace_cuboids(const vie_cuboid *a,const vie_cuboid *b,double k,complex double e,
                          const vie_halfspace_options *o,complex double n[3][3],
                          complex double h[3][3],vie_halfspace_report *r)
{
  double volume;
  if(!cuboid_valid(a) || !cuboid_valid(b)) return 0;
  volume=b->size[0]*b->size[1]*b->size[2];
  if(!isfinite(volume) || !(volume>0)) return 0;
  return evaluate(a,b,volume,k,e,o,n,h,r);
}

int vie_halfspace_point_cuboid(const double t[3],const vie_cuboid *b,double k,complex double e,
                              const vie_halfspace_options *o,complex double n[3][3],
                              complex double h[3][3],vie_halfspace_report *r)
{
  vie_cuboid a={{0,0,0},{0,0,0}};
  double volume;
  if(!point_valid(t) || !cuboid_valid(b)) return 0;
  memcpy(a.center,t,sizeof(a.center));
  volume=b->size[0]*b->size[1]*b->size[2];
  if(!isfinite(volume) || !(volume>0)) return 0;
  return evaluate(&a,b,volume,k,e,o,n,h,r);
}

int vie_halfspace_cuboid_point(const vie_cuboid *a,const double s[3],double k,complex double e,
                              const vie_halfspace_options *o,complex double n[3][3],
                              complex double h[3][3],vie_halfspace_report *r)
{
  vie_cuboid b={{0,0,0},{0,0,0}};
  if(!cuboid_valid(a) || !point_valid(s)) return 0;
  memcpy(b.center,s,sizeof(b.center));
  return evaluate(a,&b,1,k,e,o,n,h,r);
}

/* Parseval flux through z=0 from the downward spectrum of the complete
   source. Include wire/volume interference BEFORE taking squared magnitudes. */
static int power_quadrature(double k,complex double e,double clearance,
                            vie_halfspace_spectrum spectrum,void *context,
                            int order,double *power)
{
  double nodes[MAX_ORDER],weights[MAX_ORDER],sum=0,scale=k*clearance;
  int segment,j,az,p,count=5;
  double cuts[6]={0,1,4,16,64,0},branch=scale*sqrt(creal(e)-1);
  if(branch>0 && branch<64) cuts[count++]=branch;
  for(j=1;j<count;++j)
  { double value=cuts[j]; int m=j; while(m>0 && cuts[m-1]>value) { cuts[m]=cuts[m-1]; --m; } cuts[m]=value; }
  if(!gauss(order,nodes,weights)) return 0;
  for(segment=-1;segment<count-1;++segment)
  {
    double lo=segment<0 ? 0:cuts[segment],hi=segment<0 ? PI/2:cuts[segment+1];
    for(j=0;j<order;++j)
    {
      double t=.5*((hi-lo)*nodes[j]+hi+lo),u,jac;
      complex double w,v,den;
      double weight=.5*(hi-lo)*weights[j],fs,fp;
      if(segment<0) { u=sin(t); w=cos(t); jac=u*creal(w); }
      else { double ev=t/scale; u=hypot(1,ev); w=-I*ev; jac=ev/scale; }
      v=csqrt(e-u*u); if(cimag(v)>0) v=-v;
      /* Algebraically cancel |w|^2 before evaluating the interface flux.
         Avoid subtracting nearly equal Fresnel coefficients near PEC/grazing. */
      den=w+v; fs=4*creal(v)/creal(den*conj(den));
      den=e*w+v; fp=4*creal(v*conj(e))/creal(den*conj(den));
      for(az=0;az<2*order;++az)
      {
        double phi=PI*az/order,c=cos(phi),s=sin(phi);
        complex double source[3],es=0,ep=0;
        complex double te[3]={-s,c,0},pd[3]={-w*c,-w*s,-u};
        if(!spectrum(context,k,u,w,phi,source)) return 0;
        for(p=0;p<3;++p) { es+=te[p]*source[p]; ep+=pd[p]*source[p]; }
        sum+=weight*(PI/order)*jac*
             (fs*creal(es*conj(es))+fp*creal(ep*conj(ep)));
      }
    }
  }
  *power=sum*376.730313668*k*k/(32*PI*PI);
  return isfinite(*power) && *power>=0;
}

int vie_halfspace_power(double k,complex double e,double clearance,
                        vie_halfspace_spectrum spectrum,void *context,
                        const vie_halfspace_options *requested,double *power,
                        double *delta)
{
  vie_halfspace_options defaults;
  const vie_halfspace_options *o=requested;
  double old=0,value;
  int order;
  vie_halfspace_defaults(&defaults); if(!o) o=&defaults;
  if(!power || !spectrum || !isfinite(k) || k<=0 || !material_valid(e) ||
     !isfinite(clearance) || clearance<=0 || !isfinite(k*clearance) || k*clearance<=0 ||
     !isfinite(o->relative_tolerance) || o->relative_tolerance<=0 ||
     !isfinite(o->absolute_tolerance) || o->absolute_tolerance<=0 ||
     o->maximum_order<16 || o->maximum_order>MAX_ORDER) return 0;
  for(order=8;order<=o->maximum_order;order*=2)
  {
    if(!power_quadrature(k,e,clearance,spectrum,context,order,&value)) return 0;
    if(order>8)
    {
      if(delta) *delta=fabs(value-old);
      if(fabs(value-old)<=o->absolute_tolerance+o->relative_tolerance*fabs(value))
      { *power=value; return 1; }
    }
    old=value;
  }
  return 0;
}
