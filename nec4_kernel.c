#include "nec4_kernel.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

#define NEC4_ETA_OVER_4PI 29.97922085

/* Appendix B suggests moving the single-precision switches two decimal
   orders for double precision.  The M4 oracle freezes/refines these values. */
#define NEC4_FAR_RATIO 1000000.0
#define NEC4_FAR_KR (0.0001*6.2831853071795864769)
#define NEC4_EXACT_KR (0.0001*6.2831853071795864769)
#define NEC4_AXIS_RATIO 0.001
#define NEC4_SMALL_AXIS_RATIO 0.005

static double complex outgoing(double complex wave_number,double distance)
{ return cexp(-I*wave_number*distance); }

static void static_integrals(double delta,double z,double rho,
                             double *i1,double *i2,double *i3,
                             double *i4,double *i5)
{
  double s1=-delta-z,s2=delta-z;
  double r1=hypot(rho,s1),r2=hypot(rho,s2);
  *i1=asinh(s2/rho)-asinh(s1/rho);
  *i2=2.*delta;
  *i3=.5*rho*rho*(*i1)+.5*(s2*r2-s1*r1);
  *i4=rho*rho*(*i2)+(s2*s2*s2-s1*s1*s1)/3.;
  *i5=.75*rho*rho*(*i3)+.25*(s2*r2*r2*r2-s1*r1*r1*r1);
}

static double complex integral_series(double delta,double z,double rho,
                                      double complex wave_number)
{
  long double s1=-delta-z,s2=delta-z,rho2=(long double)rho*rho;
  long double r0=hypotl(rho,z),r1=hypotl(rho,s1),r2=hypotl(rho,s2);
  long double moment[42],rp1=1.,rp2=1.,r0power[42],coefficient=1.;
  long double complex lk=wave_number;
  long double complex alpha=-I*lk,power=1.,sum,term,
    phase=cexpl(-I*lk*r0);
  int m,j;
  moment[0]=asinhl(s2/rho)-asinhl(s1/rho); /* integral R^-1 */
  moment[1]=s2-s1;                         /* integral R^0  */
  r0power[0]=1.;
  for(m=1;m<42;++m) r0power[m]=r0power[m-1]*(-r0);
  sum=moment[0]*phase;
  for(m=1;m<40;++m)
  {
    int n=m-1;
    long double shifted=0.,binomial=1.;
    if(n>=1)
    {
      rp1*=r1;rp2*=r2;
      moment[m]=(s2*rp2-s1*rp1)/(n+1.)+
        n*rho2*moment[m-2]/(n+1.);
    }
    for(j=0;j<=m;++j)
    {
      shifted+=binomial*r0power[m-j]*moment[j];
      binomial*=((long double)(m-j))/(j+1.);
    }
    power*=alpha;
    coefficient/=m;
    term=phase*power*coefficient*shifted;
    sum+=term;
    if(m>6 && cabsl(term)<=2.*LDBL_EPSILON*(1.+cabsl(sum))) break;
  }
  return (double complex)sum;
}

/* One 16-point Gauss-Legendre panel on [a,b]. */
static double complex integral_gauss_panel(double a,double b,double z,
                                           double rho,double complex wave_number)
{
  static const double x[8]={
    .0950125098376374402,.281603550779258913,.458016777657227386,
    .617876244402643748,.755404408355003034,.865631202387831744,
    .944575023073232576,.989400934991649933};
  static const double w[8]={
    .189450610455068496,.182603415044923589,.169156519395002538,
    .149595988816576733,.124628971255533872,.095158511682492785,
    .062253523938647893,.027152459411754095};
  double middle=.5*(a+b),half=.5*(b-a);
  double complex sum=0.;
  size_t i;
  for(i=0;i<8;++i)
  {
    double q=half*x[i];
    double r1=hypot(rho,z-(middle-q)),r2=hypot(rho,z-(middle+q));
    sum+=w[i]*(outgoing(wave_number,r1)/r1+
               outgoing(wave_number,r2)/r2);
  }
  return half*sum;
}

static double complex integral_gauss6(double a,double b,double z,double rho,
                                      double complex wave_number)
{
  static const double x[3]={.238619186083196909,.661209386466264514,
    .932469514203152028};
  static const double w[3]={.467913934572691047,.360761573048138608,
    .171324492379170345};
  double middle=.5*(a+b),half=.5*(b-a);
  double complex sum=0.;
  size_t i;
  for(i=0;i<3;++i)
  {
    double q=half*x[i];
    double r1=hypot(rho,z-(middle-q)),r2=hypot(rho,z-(middle+q));
    sum+=w[i]*(outgoing(wave_number,r1)/r1+
               outgoing(wave_number,r2)/r2);
  }
  return half*sum;
}

static double complex integral_adaptive(double a,double b,double z,double rho,
                                        double complex wave_number,
                                        double complex whole,int depth)
{
  double middle=.5*(a+b);
  double complex left=integral_gauss_panel(a,middle,z,rho,wave_number);
  double complex right=integral_gauss_panel(middle,b,z,rho,wave_number);
  double complex split=left+right;
  if(depth<=0 || cabs(split-whole)<=8.*DBL_EPSILON*(1.+cabs(split)))
    return split;
  return integral_adaptive(a,middle,z,rho,wave_number,left,depth-1)+
    integral_adaptive(middle,b,z,rho,wave_number,right,depth-1);
}

static double complex scalar_integral(double delta,double z,double rho,
                                      double complex wave_number)
{
  double r0=hypot(rho,z);
  /* Appendix C's moment expansion, continued until the double-precision
     residual is negligible, is both faster and more accurate than its
     historical five-term truncation throughout this range. */
  if(cabs(wave_number)*delta<=1. && r0<=10.*delta)
    return integral_series(delta,z,rho,wave_number);
  if(cabs(wave_number)*delta<=1. && r0>10.*delta)
    return integral_gauss6(-delta,delta,z,rho,wave_number);
  {
    double complex whole=integral_gauss_panel(-delta,delta,z,rho,wave_number);
    return integral_adaptive(-delta,delta,z,rho,wave_number,whole,14);
  }
}

static void exact_axial(double delta,double z,double rho,double complex k,
                        double complex eta4pi,
                        double complex *ezs,double complex *ezc,
                        double complex *ezk)
{
  double r1=hypot(rho,z+delta),r2=hypot(rho,z-delta);
  double complex g1=outgoing(k,r1)/r1,g2=outgoing(k,r2)/r2;
  double complex integral=scalar_integral(delta,z,rho,k);
  double complex sd=csin(k*delta),cd=ccos(k*delta);
  double complex p=-I*eta4pi,q=-p;
  *ezs=q*(g2-g1)*cd;
  *ezk=p*k*integral;
  *ezc=*ezk+p*((g2+g1)*sd-k*integral);
}

static void exact_radial(double delta,double z,double rho,double complex k,
                         double complex eta4pi,
                         double complex *ers,double complex *erc)
{
  long double ld=delta,lz=z,lrho=rho;
  long double r1=hypotl(lrho,lz+ld),r2=hypotl(lrho,lz-ld);
  long double complex lk=k;
  long double complex e1=cexpl(-I*lk*r1),e2=cexpl(-I*lk*r2);
  long double complex g1=e1/r1,g2=e2/r2;
  long double complex sd=csinl(lk*ld),cd=ccosl(lk*ld);
  long double complex p=-I*(long double complex)eta4pi,q=-p;
  *ers=(double complex)(p*(((lz-ld)*g2-(lz+ld)*g1)*cd-
                           I*(e2+e1)*sd)/lrho);
  *erc=(double complex)(q*(((lz-ld)*g2+(lz+ld)*g1)*sd+
                           I*(e2-e1)*cd)/lrho);
}

static void near_axis_radial(double delta,double z,double rho,double complex k,
                             double complex eta4pi,
                             double complex *ers,double complex *erc)
{
  double az=fabs(z),sz=z<0.?-1.:1.,den=z*z-delta*delta;
  double complex em=outgoing(k,az),edm=outgoing(k,delta),
    edp=cexp(I*k*delta);
  double complex sd=csin(k*delta),cd=ccos(k*delta);
  double complex p=-I*eta4pi;
  *ers=p*rho*em*(.5*(edm/((az+delta)*(az+delta))-
    edp/((az-delta)*(az-delta)))*cd-I*k*delta/den);
  *erc=p*sz*rho*em*(.5*(edm/((az+delta)*(az+delta))+
    edp/((az-delta)*(az-delta)))*sd-k*delta/den);
}

static void small_axial(double delta,double z,double rho,double complex k,
                        double complex eta4pi,
                        double complex *ezs,double complex *ezc,
                        double complex *ezk)
{
  double r1=hypot(rho,z+delta),r2=hypot(rho,z-delta);
  double i1,i2,i3,i4,i5;
  double complex sd=csin(k*delta),cd=ccos(k*delta);
  double r1_2=r1*r1,r2_2=r2*r2,r1_3=r1_2*r1,r2_3=r2_2*r2;
  double complex k2=k*k,k3=k2*k,k4=k2*k2;
  double complex p=-I*eta4pi,q=-p,difference;
  static_integrals(delta,z,rho,&i1,&i2,&i3,&i4,&i5);
  *ezs=q*(1./r2-1./r1-.5*k2*(r2-r1)-
          2.*I*k3*z*delta/3.+k4*(r2_3-r1_3)/24.)*cd;
  *ezk=p*k*(i1-2.*I*k*delta-.5*k2*i3+I*k3*i4/6.+k4*i5/24.);
  difference=p*((1./r2+1./r1-.5*k2*(r2+r1)+
    I*k3*(r2_2+r1_2)/6.+k4*(r2_3+r1_3)/24.)*sd+
    I*k4*delta*delta*delta/3.-
    k*(i1-.5*k2*i3+I*k3*i4/6.+k4*i5/24.));
  *ezc=*ezk+difference;
}

static void small_radial(double delta,double z,double rho,double complex k,
                         double complex eta4pi,
                         double complex *ers,double complex *erc)
{
  long double ld=delta,lz=z,lrho=rho;
  long double complex lk=k;
  long double r1=hypotl(lrho,lz+ld),r2=hypotl(lrho,lz-ld);
  long double r1_3=r1*r1*r1,r2_3=r2*r2*r2;
  long double complex k2=lk*lk,k3=k2*lk,k4=k2*k2,k6=k3*k3;
  long double complex sd=csinl(lk*ld),cd=ccosl(lk*ld);
  long double complex p=-I*(long double complex)eta4pi,q=-p;
  *ers=(double complex)(p*(lz*(1./r2-1./r1-.5*k2*(r2-r1)+
      k4*(r2_3-r1_3)/24.)*cd-
    ld*(1./r2+1./r1-.5*k2*(r2+r1)+
      k4*(r2_3+r1_3)/24.)*cd-
    (lk*(r2+r1)-k3*(r2_3+r1_3)/6.)*sd+
    2.*I*k3*lrho*lrho*ld/3.)/lrho);
  *erc=(double complex)(q*(lz*(1./r2+1./r1-.5*k2*(r2+r1)+
      k4*(r2_3+r1_3)/24.)*sd-
    ld*(1./r2-1./r1-.5*k2*(r2-r1)+
      k4*(r2_3-r1_3)/24.)*sd+
    (lk*(r2-r1)-k3*(r2_3-r1_3)/6.)*cd+
    2.*I*k6*ld*ld*ld*lrho*lrho*lz/45.)/lrho);
}

static void small_axis_radial(double delta,double z,double rho,double complex k,
                              double complex eta4pi,
                              double complex *ers,double complex *erc)
{
  double az=fabs(z),sz=z<0.?-1.:1.,d2=delta*delta;
  double den=z*z-d2;
  double complex k2=k*k,k3=k2*k,k5=k3*k2;
  double complex p=-I*eta4pi,q=-p;
  *ers=q*rho*delta*(az*(2.+k2*(z*z-2.*d2))/(den*den)-
                    2.*I*k3/3.);
  *erc=p*sz*k*rho*delta*d2*((6.+k2*(z*z-2.*d2))/(3.*den*den)-
                            2.*I*k5*az/45.);
}

/* Continue the Appendix B on-axis limit through rho^5.  Evaluating the
   direct expression at two safe auxiliary radii supplies the even rho^2 and
   rho^4 coefficients without exposing the requested point to cancellation. */
static void axis_corrected_radial(double delta,double z,double rho,
                                  double complex k,double complex eta4pi,
                                  int small,double complex *ers,
                                  double complex *erc)
{
  double minimum=fmin(hypot(rho,z+delta),hypot(rho,z-delta));
  double h=(small?NEC4_SMALL_AXIS_RATIO:NEC4_AXIS_RATIO)*minimum;
  double h2=h*h,h4=h2*h2,x=rho*rho;
  double complex a_s,a_c,e1s,e1c,e2s,e2c,d1,d2,b,c;
  if(small) small_axis_radial(delta,z,1.,k,eta4pi,&a_s,&a_c);
  else near_axis_radial(delta,z,1.,k,eta4pi,&a_s,&a_c);
  exact_radial(delta,z,h,k,eta4pi,&e1s,&e1c);
  exact_radial(delta,z,2.*h,k,eta4pi,&e2s,&e2c);
  d1=e1s/h-a_s;d2=e2s/(2.*h)-a_s;
  c=(d2-4.*d1)/(12.*h4);b=(d1-c*h4)/h2;
  *ers=rho*(a_s+b*x+c*x*x);
  d1=e1c/h-a_c;d2=e2c/(2.*h)-a_c;
  c=(d2-4.*d1)/(12.*h4);b=(d1-c*h4)/h2;
  *erc=rho*(a_c+b*x+c*x*x);
}

static void far_field(double delta,double z,double rho,double complex k,
                      double complex eta4pi,int small,nec4_segment_fields *f)
{
  double r=hypot(rho,z),r2=r*r,r3=r2*r,r4=r2*r2,r5=r4*r;
  double d2=delta*delta,d3=d2*delta;
  double complex k2=k*k,k3=k2*k,k4=k2*k2,k5=k4*k;
  double complex p=-I*eta4pi,q=-p,difference;
  if(!small)
  {
    double complex e=outgoing(k,r),common=(1.+I*k*r)*e;
    f->ers=q*2.*delta*rho/r3*common;
    f->erc=p*2.*k*d3*rho*z/(3.*r5)*
      (1.+I*k*r-k2*r2/3.)*e;
    f->ezs=q*2.*delta*z/r3*common;
    f->ezk=p*k*delta/r*(2.+d2/(3.*r2)*
      ((2.*z*z-rho*rho)/r2*(1.+I*k*r)-k2*z*z))*e;
    difference=p*k*d3/(3.*r3)*
      (2.*(2.*z*z-rho*rho)/r2*(1.+I*k*r)-k2*(2.*z*z+r2))*e;
  }
  else
  {
    double complex poly=1.+k2*r2/2.-I*k3*r3/3.-k4*r4/8.+I*k5*r5/30.;
    f->ers=q*2.*rho*delta/r3*poly;
    f->erc=p*2.*k*rho*z*d3/(3.*r5)*
      (3.+k2*r2/2.+k4*r4/8.-I*k5*r5/15.);
    f->ezs=q*2.*z*delta/r3*poly;
    f->ezk=p*k*delta/(3.*r5)*(6.*r4+(2.*z*z-rho*rho)*d2-
      6.*I*k*r5-k2*r2*(3.*r4+rho*rho*d2/2.)+
      I*k3*r5*(r2+d2));
    difference=p*k*d3/(3.*r5)*(2.*(2.*z*z-rho*rho)-
      k2*r2*(z*z+2.*rho*rho)+5.*I*k3*r5/3.+
      k4*r4*(4.*z*z+3.*rho*rho)/4.+
      I*k5*r5*(13.*z*z-9.*rho*rho)/120.);
  }
  f->ezc=f->ezk+difference;
}

static int reduced_segment_field(double s,double z,double rho,
    double complex k,double complex eta4pi,nec4_segment_fields *f)
{
  double delta=.5*s,r0=hypot(rho,z),r1=hypot(rho,z+delta),
    r2=hypot(rho,z-delta),minimum=fmin(r1,r2),kr=cabs(k)*r0;
  if(!f || !(s>0.) || !(rho>0.) || !(cabs(k)>0.) || !(cabs(eta4pi)>0.) ||
     !isfinite(r0)||!isfinite(creal(k))||!isfinite(cimag(k))||
     !isfinite(creal(eta4pi))||!isfinite(cimag(eta4pi))) return 0;
  f->erk=0.;
  if(r0>NEC4_FAR_RATIO*delta)
  {
    int small=kr<=NEC4_FAR_KR;
    far_field(delta,z,rho,k,eta4pi,small,f);
    f->branch=small?NEC4_FIELD_FAR_SMALL_KR:NEC4_FIELD_FAR;
  }
  else if(kr>NEC4_EXACT_KR)
  {
    exact_axial(delta,z,rho,k,eta4pi,&f->ezs,&f->ezc,&f->ezk);
    if(rho>NEC4_AXIS_RATIO*minimum)
    { exact_radial(delta,z,rho,k,eta4pi,&f->ers,&f->erc);f->branch=NEC4_FIELD_EXACT; }
    else
    { axis_corrected_radial(delta,z,rho,k,eta4pi,0,&f->ers,&f->erc);
      f->branch=NEC4_FIELD_NEAR_AXIS; }
  }
  else
  {
    small_axial(delta,z,rho,k,eta4pi,&f->ezs,&f->ezc,&f->ezk);
    if(rho>NEC4_SMALL_AXIS_RATIO*minimum)
    { small_radial(delta,z,rho,k,eta4pi,&f->ers,&f->erc);f->branch=NEC4_FIELD_SMALL_KR; }
    else
    { axis_corrected_radial(delta,z,rho,k,eta4pi,1,&f->ers,&f->erc);
      f->branch=NEC4_FIELD_SMALL_KR_AXIS; }
  }
  return isfinite(creal(f->ezs)) && isfinite(cimag(f->ezs)) &&
    isfinite(creal(f->ers)) && isfinite(cimag(f->ers)) &&
    isfinite(creal(f->ezc)) && isfinite(cimag(f->ezc)) &&
    isfinite(creal(f->erc)) && isfinite(cimag(f->erc)) &&
    isfinite(creal(f->ezk)) && isfinite(cimag(f->ezk));
}

int nec4_reduced_segment_field_medium(double s,double z,double rho,
    double complex k,double complex eta,nec4_segment_fields *f)
{
  return reduced_segment_field(s,z,rho,k,eta/(4.*acos(-1.)),f);
}

int nec4_reduced_segment_field(double s,double z,double rho,double k,
                               nec4_segment_fields *f)
{
  /* Preserve the historical NEC4 free-space constant bit-for-bit. */
  return reduced_segment_field(s,z,rho,k,NEC4_ETA_OVER_4PI,f);
}

const char *nec4_field_branch_name(nec4_field_branch branch)
{
  static const char *const name[]={"exact","near-axis","small-kR",
    "small-kR-axis","far","far-small-kR"};
  return branch<=NEC4_FIELD_FAR_SMALL_KR?name[branch]:"unknown";
}
