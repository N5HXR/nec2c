/* Complex cylindrical functions shared by the legacy and NEC4 ground paths. */

#include "nec2c.h"

/* Evaluate J0(z) and its derivative. This is the original SOMNEC algorithm,
   moved unchanged apart from making its one-time coefficient initialization
   safe before parallel NEC4 table generation. */
void bessel(complex double z,complex double *j0,complex double *j0p)
{
  int k,ib;
  static _Thread_local int initialized,m[101];
  static _Thread_local double a1[25],a2[25];
  double zms;
  complex double p0z,p1z,q0z,q1z,zi,zi2,zk,cz,sz,j0x=CPLX_00,j0px=CPLX_00;
  if(!initialized)
  {
    int i;
    for(k=1;k<=25;k++)
    {
      i=k-1;
      a1[i]=-.25/(k*k);
      a2[i]=1.0/(k+1.0);
    }
    for(i=1;i<=101;i++)
    {
      double tst=1.0;
      int terms=1;
      for(k=0;k<24;k++)
      {
        terms=k+1;
        tst *= -i*a1[k];
        if(tst<1.0e-6) break;
      }
      m[i-1]=terms;
    }
    initialized=1;
  }

  zms=creal(z*conj(z));
  if(zms<=1.e-12)
  {
    *j0=CPLX_10;
    *j0p=-.5*z;
    return;
  }

  ib=0;
  if(zms<=37.21)
  {
    int iz,miz;
    if(zms>36.) ib=1;
    iz=(int)zms;
    miz=m[iz];
    *j0=CPLX_10;
    *j0p=*j0;
    zk=*j0;
    zi=z*z;
    for(k=0;k<miz;k++)
    {
      zk *= a1[k]*zi;
      *j0 += zk;
      *j0p += a2[k]*zk;
    }
    *j0p *= -.5*z;
    if(ib==0) return;
    j0x=*j0;
    j0px=*j0p;
  }

  zi=1./z;
  zi2=zi*zi;
  p0z=1.+(P20*zi2-P10)*zi2;
  p1z=1.+(P11-P21*zi2)*zi2;
  q0z=(Q20*zi2-Q10)*zi;
  q1z=(Q11-Q21*zi2)*zi;
  zk=cexp(CPLX_01*(z-POF));
  zi2=1./zk;
  cz=.5*(zk+zi2);
  sz=CPLX_01*.5*(zi2-zk);
  zk=C3*csqrt(zi);
  *j0=zk*(p0z*cz-q0z*sz);
  *j0p=-zk*(p1z*sz+q1z*cz);
  if(ib==0) return;
  zms=cos((sqrt(zms)-6.)*PI10);
  *j0=.5*(j0x*(1.+zms)+*j0*(1.-zms));
  *j0p=.5*(j0px*(1.+zms)+*j0p*(1.-zms));
}

static void hankel_kind(complex double z,double sign,complex double *h0,
    complex double *h0p)
{
  int k,ib;
  static _Thread_local int initialized,m[101];
  static _Thread_local double a1[25],a2[25],a3[25],a4[25];
  double zms;
  complex double clogz,j0,j0p,p0z,p1z,q0z,q1z,y0=CPLX_00,y0p=CPLX_00;
  complex double zi,zi2,zk;
  if(!initialized)
  {
    int i;
    double psi=-GAMMA;
    for(k=1;k<=25;k++)
    {
      i=k-1;
      a1[i]=-.25/(k*k);a2[i]=1./(k+1.);
      psi+=1./k;a3[i]=2.*psi;
      a4[i]=(2.*psi+1./(k+1.))/(k+1.);
    }
    for(i=1;i<=101;i++)
    {
      double test=1.;
      int terms=1;
      for(k=0;k<24;k++)
      {
        terms=k+1;test*=-i*a1[k];
        if(test*a3[k]<1.e-6) break;
      }
      m[i-1]=terms;
    }
    initialized=1;
  }
  zms=creal(z*conj(z));
  if(zms==0.) abort();
  ib=0;
  if(zms<=16.81)
  {
    int iz=(int)zms,miz=m[iz];
    if(zms>16.) ib=1;
    j0=CPLX_10;j0p=j0;zk=j0;
    zi=z*z;
    for(k=0;k<miz;k++)
    {
      zk*=a1[k]*zi;j0+=zk;j0p+=a2[k]*zk;
      y0+=a3[k]*zk;y0p+=a4[k]*zk;
    }
    j0p*=-.5*z;clogz=clog(.5*z);
    y0=(2.*j0*clogz-y0)/PI+C2;
    y0p=(2./z+2.*j0p*clogz+.5*y0p*z)/PI+C1*z;
    *h0=j0+sign*CPLX_01*y0;
    *h0p=j0p+sign*CPLX_01*y0p;
    if(ib==0) return;
    y0=*h0;y0p=*h0p;
  }
  zi=1./z;zi2=zi*zi;
  p0z=1.+(P20*zi2-P10)*zi2;
  p1z=1.+(P11-P21*zi2)*zi2;
  q0z=(Q20*zi2-Q10)*zi;
  q1z=(Q11-Q21*zi2)*zi;
  zk=cexp(sign*CPLX_01*(z-POF))*csqrt(zi)*C3;
  *h0=zk*(p0z+sign*CPLX_01*q0z);
  *h0p=sign*CPLX_01*zk*(p1z+sign*CPLX_01*q1z);
  if(ib==0) return;
  zms=cos((sqrt(zms)-4.)*PI10);
  *h0=.5*(y0*(1.+zms)+*h0*(1.-zms));
  *h0p=.5*(y0p*(1.+zms)+*h0p*(1.-zms));
}

/* Legacy SOMNEC convention: first-kind Hankel function and derivative. */
void hankel(complex double z,complex double *h0,complex double *h0p)
{ hankel_kind(z,1.,h0,h0p); }

/* NEC4 Figure 5-3 full-line representation uses the second kind. */
void hankel2(complex double z,complex double *h0,complex double *h0p)
{ hankel_kind(z,-1.,h0,h0p); }
