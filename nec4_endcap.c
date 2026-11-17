#include "nec4_endcap.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

static void gauss_legendre(int n,double *x,double *w)
{
  int i,j,m=(n+1)/2;
  for(i=0;i<m;++i)
  {
    double z=cos(M_PI*((double)i+.75)/((double)n+.5)),last;
    double p1=0.,p2=0.,pp=0.;
    do
    {
      p1=1.;p2=0.;
      for(j=1;j<=n;++j)
      {
        double p3=p2;
        p2=p1;
        p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.);
      last=z;
      z=last-p1/pp;
    } while(fabs(z-last)>4.*DBL_EPSILON);
    x[i]=-z;x[n-1-i]=z;
    w[i]=2./((1.-z*z)*pp*pp);w[n-1-i]=w[i];
  }
}

static void disk_quadrature(double a,double rho,double z,double k,int nr,
                            double complex *radial,double complex *axial)
{
  double xr[64],wr[64],xp[128],wp[128];
  double complex er=0.,ez=0.;
  int i,j,np=2*nr;
  gauss_legendre(nr,xr,wr);
  gauss_legendre(np,xp,wp);
  for(i=0;i<nr;++i)
  {
    double rp=.5*a*(xr[i]+1.),rw=.5*a*wr[i]*rp;
    for(j=0;j<np;++j)
    {
      double phi=M_PI*(xp[j]+1.),c=cos(phi);
      double dx=rho-rp*c,r2=dx*dx+rp*rp*(1.-c*c)+z*z;
      double r=sqrt(r2),weight=rw*M_PI*wp[j];
      double complex q=cexp(-I*k*r)*(1.+I*k*r)/(r2*r)*weight;
      er+=q*dx;
      ez+=q*z;
    }
  }
  *radial=er;*axial=ez;
}

static double scaled_delta(double complex a,double complex b)
{
  return cabs(a-b)/(1.+fmax(cabs(a),cabs(b)));
}

int nec4_endcap_disk_field(double a,double rho,double z,double k,
                           nec4_endcap_field *f)
{
  double r0;
  if(!f || !(a>0.) || rho<0. || !(k>0.) ||
     !isfinite(a) || !isfinite(rho) || !isfinite(z)) return 0;
  r0=hypot(rho,z);
  if(rho<=64.*DBL_EPSILON*fmax(a,fabs(z)))
  {
    double az=fabs(z),ra=hypot(a,z),sgn=z<0.?-1.:1.;
    if(az<=DBL_MIN) return 0;
    f->radial=0.;
    f->axial=2.*M_PI*sgn*cexp(-I*k*az)*
      (1.-az/ra*cexp(-I*k*(ra-az)));
    f->branch=NEC4_ENDCAP_AXIS;
  }
  else if(r0>1000.*a)
  {
    double r2=r0*r0,a2=a*a,r4=r2*r2,r3=r2*r0;
    double complex e=cexp(-I*k*r0),base=1.+I*k*r0;
    double complex q=I*k*k*k*r3+6.*k*k*r2-15.*base;
    f->axial=M_PI*z*a2*e/r3*(base-a2/(4.*r2)*
      (3.*base-k*k*r2)-rho*rho*a2/(8.*r4)*q);
    f->radial=M_PI*rho*a2*e/r3*(base-a2/(2.*r2)*
      (3.*base-k*k*r2)-rho*rho*a2/(8.*r4)*q);
    f->branch=NEC4_ENDCAP_FAR;
  }
  else
  {
    double complex er0,ez0,er1,ez1;
    int n=12;
    if(fabs(z)<16.*a && rho<2.*a) n=24;
    if(n>32) n=32;
    disk_quadrature(a,rho,z,k,n,&er0,&ez0);
    for(;;)
    {
      int next=2*n;
      if(next>64) next=64;
      disk_quadrature(a,rho,z,k,next,&er1,&ez1);
      if((scaled_delta(er0,er1)<=2.e-13 &&
          scaled_delta(ez0,ez1)<=2.e-13) || next==64) break;
      er0=er1;ez0=ez1;n=next;
    }
    f->radial=er1;f->axial=ez1;
    f->branch=NEC4_ENDCAP_QUADRATURE;
  }
  return isfinite(creal(f->radial)) && isfinite(cimag(f->radial)) &&
    isfinite(creal(f->axial)) && isfinite(cimag(f->axial));
}

const char *nec4_endcap_branch_name(nec4_endcap_branch branch)
{
  static const char *const names[]={"axis","quadrature","far"};
  return branch<=NEC4_ENDCAP_FAR?names[branch]:"unknown";
}
