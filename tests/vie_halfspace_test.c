#include "vie_halfspace.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#define PI 3.14159265358979323846
static int failures;
#define CHECK(x,m) do { if(!(x)) { fprintf(stderr,"FAIL: %s\n",m); ++failures; } } while(0)

static void test_fresnel(void)
{
  complex double rs,rp;
  CHECK(vie_halfspace_fresnel(4,0,&rs,&rp),"normal incidence");
  CHECK(cabs(rs+1.0/3)<1e-14 && cabs(rp-1.0/3)<1e-14,"normal Fresnel amplitudes");
  CHECK(vie_halfspace_fresnel(4,sqrt(.8),&rs,&rp),"Brewster incidence");
  CHECK(cabs(rp)<1e-14,"Brewster zero");
  CHECK(vie_halfspace_fresnel(4,1,&rs,&rp),"grazing incidence");
  CHECK(cabs(rs+1)<1e-14 && cabs(rp+1)<1e-14,"grazing limits");
  CHECK(vie_halfspace_fresnel(4-.5*I,1e5,&rs,&rp),"large evanescent wave number");
  CHECK(cabs(rp-(3-.5*I)/(5-.5*I))<1e-9 && cabs(rs)<1e-9,"electrostatic reflection limit");
  CHECK(!vie_halfspace_fresnel(4+.1*I,0,&rs,&rp),"reject active medium");
}

static void test_pec_point(void)
{
  double a[3]={.12,-.1,.4},b[3]={-.08,.05,.3},d[3];
  complex double n[3][3],h[3][3],g,diag,radial,gradient[3];
  double r=0,k=2.1,errn=0,errh=0;
  vie_halfspace_report report;
  int p,q;
  CHECK(vie_halfspace_point(a,b,k,1e24,NULL,n,h,&report),"PEC point limit convergence");
  for(p=0;p<3;++p) { d[p]=a[p]-(p==2 ? -b[p]:b[p]); r+=d[p]*d[p]; }
  r=sqrt(r); g=cexp(-I*k*r)/(4*PI*r);
  diag=g*(k*k-I*k/r-1/(r*r)); radial=g*(-k*k+3*I*k/r+3/(r*r));
  for(p=0;p<3;++p) gradient[p]=g*(-I*k-1/r)*d[p]/r;
  for(p=0;p<3;++p) for(q=0;q<3;++q)
  {
    complex double expected=(q==2 ? 1:-1)*((p==q ? diag:0)+radial*d[p]*d[q]/(r*r));
    complex double hm=0;
    if(p==(q+1)%3) hm=gradient[(q+2)%3];
    if(p==(q+2)%3) hm=-gradient[(q+1)%3];
    hm*=q==2 ? 1:-1;
    errn=fmax(errn,cabs(n[p][q]-expected)); errh=fmax(errh,cabs(h[p][q]-hm));
  }
  CHECK(errn<3e-9 && errh<3e-9,"PEC electric and magnetic closed-form image dyads");
  printf("PEC point maximum electric/magnetic error %.3g / %.3g, order %d\n",errn,errh,report.order);
}

static void test_cuboids(void)
{
  vie_cuboid a={{.1,-.1,.4},{.06,.04,.05}},b={{-.1,.05,.3},{.04,.03,.06}},image=b;
  complex double n[3][3],h[3][3],reverse[3][3],hr[3][3],exact[3][3];
  double va=a.size[0]*a.size[1]*a.size[2],vb=b.size[0]*b.size[1]*b.size[2],err=0;
  int p,q;
  CHECK(vie_halfspace_cuboids(&a,&b,2,4-.4*I,NULL,n,h,NULL),"finite lossy cuboids");
  CHECK(vie_halfspace_cuboids(&b,&a,2,4-.4*I,NULL,reverse,hr,NULL),"reverse lossy cuboids");
  for(p=0;p<3;++p) for(q=0;q<3;++q)
    CHECK(cabs(va*n[p][q]-vb*reverse[q][p])<1e-14,"finite-volume reciprocity");
  CHECK(vie_halfspace_cuboids(&a,&b,2,1e24,NULL,n,h,NULL),"PEC cuboid limit");
  image.center[2]=-image.center[2];
  CHECK(vie_cuboid_n(&a,&image,2,NULL,exact,NULL),"independent real-space cuboid oracle");
  for(p=0;p<3;++p) for(q=0;q<3;++q)
    err=fmax(err,cabs(n[p][q]-(q==2 ? 1:-1)*exact[p][q]));
  CHECK(err<1e-11,"spectral cuboid averages match real-space PEC image kernel");
  printf("PEC finite-cuboid maximum error %.3g\n",err);
}

static void test_mixed(void)
{
  double point[3]={.12,-.1,.4};
  vie_cuboid cell={{-.08,.05,.3},{.04,.03,.06}};
  const double nodes[3]={-.7745966692414834,0,.7745966692414834};
  const double weights[3]={5.0/9,8.0/9,5.0/9};
  complex double n[3][3],h[3][3],reverse[3][3],hr[3][3];
  complex double oracle_n[3][3]={{0}},oracle_h[3][3]={{0}},saved[3][3];
  double volume=cell.size[0]*cell.size[1]*cell.size[2],k=2.1,errn=0,errh=0;
  int p,q,ix,iy,iz;
  CHECK(vie_halfspace_point_cuboid(point,&cell,k,4-.4*I,NULL,n,h,NULL),"point testing finite lossy source");
  CHECK(vie_halfspace_cuboid_point(&cell,point,k,4-.4*I,NULL,reverse,hr,NULL),"cell testing point lossy source");
  for(p=0;p<3;++p) for(q=0;q<3;++q)
    CHECK(cabs(n[p][q]-volume*reverse[q][p])<1e-14,"mixed electric reciprocity and source normalization");
  CHECK(vie_halfspace_point_cuboid(point,&cell,k,1e24,NULL,n,h,NULL),"mixed PEC convergence");
  /* Independent real-space Gauss integration of the closed-form image dipole.
     Tests magnetic curl signs as well as the electric source-volume factor. */
  for(ix=0;ix<3;++ix) for(iy=0;iy<3;++iy) for(iz=0;iz<3;++iz)
  {
    int index[3]={ix,iy,iz};
    double d[3],radius=0,weight=volume*weights[ix]*weights[iy]*weights[iz]/8;
    complex double g,diag,radial,gradient[3];
    for(p=0;p<3;++p)
    {
      double source=cell.center[p]+.5*cell.size[p]*nodes[index[p]];
      d[p]=point[p]-(p==2 ? -source:source); radius+=d[p]*d[p];
    }
    radius=sqrt(radius); g=cexp(-I*k*radius)/(4*PI*radius);
    diag=g*(k*k-I*k/radius-1/(radius*radius));
    radial=g*(-k*k+3*I*k/radius+3/(radius*radius));
    for(p=0;p<3;++p) gradient[p]=g*(-I*k-1/radius)*d[p]/radius;
    for(p=0;p<3;++p) for(q=0;q<3;++q)
    {
      double image_sign=q==2 ? 1:-1;
      complex double hm=0;
      if(p==(q+1)%3) hm=gradient[(q+2)%3];
      if(p==(q+2)%3) hm=-gradient[(q+1)%3];
      oracle_n[p][q]+=weight*image_sign*((p==q ? diag:0)+radial*d[p]*d[q]/(radius*radius));
      oracle_h[p][q]+=weight*image_sign*hm;
    }
  }
  for(p=0;p<3;++p) for(q=0;q<3;++q)
  { errn=fmax(errn,cabs(n[p][q]-oracle_n[p][q])); errh=fmax(errh,cabs(h[p][q]-oracle_h[p][q])); }
  CHECK(errn<1e-11 && errh<1e-11,"mixed finite-source E/H match independent PEC image integration");
  memcpy(saved,n,sizeof(n)); cell.size[0]=0;
  CHECK(!vie_halfspace_point_cuboid(point,&cell,k,4,NULL,n,h,NULL),"mixed rejects degenerate source cell");
  CHECK(!vie_halfspace_cuboid_point(&cell,point,k,4,NULL,n,h,NULL),"mixed rejects degenerate testing cell");
  CHECK(memcmp(saved,n,sizeof(n))==0,"mixed invalid geometry preserves outputs");
  printf("Mixed PEC electric/magnetic integration error %.3g / %.3g\n",errn,errh);
}

static void test_independent_lossy_reference(void)
{
  /* tools/m10_halfspace_reference.py analytically integrates azimuth and
     independently integrates the resulting scalar integrals by adaptive Simpson. */
  double a[3]={0,0,.4},b[3]={0,0,.3};
  complex double n[3][3],h[3][3];
  vie_halfspace_options options;
  vie_halfspace_report report;
  CHECK(vie_halfspace_point(a,b,2.1,4-.4*I,NULL,n,h,&report),"independent lossy axial reference convergence");
  CHECK(cabs(n[0][0]-(.15390406579571347-.0049065908033000816*I))<2e-8,"independent finite-frequency Nxx");
  CHECK(cabs(n[2][2]-(.33189432434537552-.34517117974856282*I))<2e-8,"independent finite-frequency Nzz");
  CHECK(cabs(h[0][1]-(-.085453897701260606+.047022971888858248*I))<2e-8,"independent finite-frequency Hxy");
  vie_halfspace_defaults(&options); options.maximum_order=256;
  CHECK(vie_halfspace_point(a,b,2.1,4,&options,n,h,NULL),"lossless substrate branch handling");
  a[0]=3; options.maximum_order=16; options.relative_tolerance=1e-12; options.absolute_tolerance=1e-14;
  CHECK(!vie_halfspace_point(a,b,2.1,4-.4*I,&options,n,h,&report),"unresolved integral fails at maximum order");
  printf("Lossy axial reference and lossless branch checks passed; unresolved order %d rejected\n",report.order);
}

static void test_quasistatic_and_rejection(void)
{
  double a[3]={0,0,.2},b[3]={0,0,.3};
  complex double n[3][3],h[3][3],saved[3][3],r=(3-.5*I)/(5-.5*I);
  double height=.5,base=1/(4*PI*height*height*height);
  int p,q;
  CHECK(vie_halfspace_point(a,b,.001,4-.5*I,NULL,n,h,NULL),"quasistatic lossy limit");
  CHECK(cabs(n[0][0]-r*base)<2e-6 && cabs(n[2][2]-2*r*base)<4e-6,"quasistatic dielectric image");
  CHECK(cimag(n[2][2])<0,"passive quasistatic substrate loss sign");
  CHECK(vie_halfspace_point(a,b,1,1,NULL,n,h,NULL),"matched background");
  for(p=0;p<3;++p) for(q=0;q<3;++q) CHECK(n[p][q]==0 && h[p][q]==0,"matched background exactly zero");
  for(p=0;p<3;++p) for(q=0;q<3;++q) n[p][q]=7+I;
  memcpy(saved,n,sizeof(n));
  CHECK(!vie_halfspace_point(a,b,0,4,NULL,n,h,NULL),"reject static API use");
  CHECK(!vie_halfspace_point(a,b,1,4+I,NULL,n,h,NULL),"reject active substrate");
  b[2]=0;
  CHECK(!vie_halfspace_point(a,b,1,4,NULL,n,h,NULL),"reject interface source");
  CHECK(memcmp(n,saved,sizeof(n))==0,"failure preserves output");
}
static int vertical_dipole_spectrum(void *context,double k,double u,
                                    complex double w,double phi,complex double source[3])
{
  double height=*(double *)context;
  (void)u; (void)phi;
  source[0]=source[1]=0; source[2]=cexp(-I*k*w*height);
  return 1;
}

static void test_ground_power(void)
{
  double height=.35,k=2.1,power=-1,delta,upper=0;
  double point[3]={0,0,.35};
  complex double n[3][3],h[3][3],epsilon=4-.4*I;
  int j;
  CHECK(vie_halfspace_power(k,epsilon,height,vertical_dipole_spectrum,&height,NULL,&power,&delta),
        "lossy lower-half-space flux convergence");
  CHECK(vie_halfspace_point(point,point,k,epsilon,NULL,n,h,NULL),"dipole reflected self field");
  /* Independent upper-hemisphere angular integration (Simpson in cos theta).
     Total emitted power comes independently from the dipole's local E field. */
  for(j=0;j<=2048;++j)
  {
    double w=(double)j/2048,u=sqrt(1-w*w);
    complex double rs,rp,a;
    CHECK(vie_halfspace_fresnel(epsilon,u,&rs,&rp),"power angular Fresnel evaluation");
    a=cexp(I*k*w*height)+rp*cexp(-I*k*w*height);
    upper+=(j==0 || j==2048 ? 1:j%2 ? 4:2)*(1-w*w)*creal(a*conj(a));
  }
  upper*=376.730313668*k*k/(16*PI*3*2048);
  {
    double emitted=376.730313668*k*k/(12*PI)-.5*376.730313668/k*cimag(n[2][2]);
    CHECK(power>0 && fabs(emitted-upper-power)<2e-6,"independent dipole work equals air radiation plus ground flux");
    printf("Lossy dipole power: work %.12g, air %.12g, ground %.12g, defect %.3g\n",
           emitted,upper,power,fabs(emitted-upper-power));
  }
  CHECK(vie_halfspace_power(k,1,height,vertical_dipole_spectrum,&height,NULL,&power,&delta),
        "matched lower half-space transmission");
  CHECK(fabs(power-376.730313668*k*k/(24*PI))<1e-7,"matched background transmits half free-space dipole power");
  /* A huge finite dielectric develops a narrow grazing transmission layer.
     Unlike the E/H image limit, this flux integral is unresolved at max order;
     it must fail explicitly rather than certify an apparent PEC zero. */
  power=7;
  CHECK(!vie_halfspace_power(k,1e24,height,vertical_dipole_spectrum,&height,NULL,&power,&delta) && power==7,
        "unresolved grazing transmission limit fails and preserves power");
  power=7;
  CHECK(!vie_halfspace_power(k,epsilon,0,vertical_dipole_spectrum,&height,NULL,&power,&delta) && power==7,
        "invalid flux input preserves power");
}

static void test_remote_ground(void)
{
  /* Translate both localized sources away from the interface while preserving
     their direct separation. Reflected E/H, rather than plane-wave illumination
     or angular interference patterns, must tend to zero. */
  double a[3]={0,0,.5},b[3]={.1,0,.5},near_e=0,near_h=0,far_e=0,far_h=0;
  complex double n[3][3],h[3][3];
  int p,q;
  CHECK(vie_halfspace_point(a,b,.01,4-.4*I,NULL,n,h,NULL),"near remote-ground reference");
  for(p=0;p<3;++p) for(q=0;q<3;++q) {
    near_e=fmax(near_e,cabs(n[p][q]));near_h=fmax(near_h,cabs(h[p][q]));
  }
  a[2]=b[2]=4;
  CHECK(vie_halfspace_point(a,b,.01,4-.4*I,NULL,n,h,NULL),"remote-ground convergence");
  for(p=0;p<3;++p) for(q=0;q<3;++q) {
    far_e=fmax(far_e,cabs(n[p][q]));far_h=fmax(far_h,cabs(h[p][q]));
  }
  CHECK(near_e>0 && near_h>0 && far_e<near_e/100 && far_h<near_h/40,
        "localized reflected E/H vanishes as ground recedes");
  printf("Remote-ground electric/magnetic ratios %.6g / %.6g\n",far_e/near_e,far_h/near_h);
}

int main(void)
{
  test_ground_power(); test_fresnel(); test_pec_point(); test_cuboids(); test_mixed(); test_independent_lossy_reference(); test_quasistatic_and_rejection(); test_remote_ground();
  if(failures) return 1;
  puts("Finite-ground reflected E/H operator checks passed"); return 0;
}
