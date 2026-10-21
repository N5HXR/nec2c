#include "dielectric.h"
#include "shared.h"
#include "vie_ground.h"
#include "vie_halfspace.h"

#include <math.h>
#include <string.h>

/* Called before the coupled solver normalizes geometry to wavelengths. */
int dielectric_ground_validate(const dielectric_model *model,
                               const dielectric_mesh_result *mesh, int iterative,
                               const dielectric_plane_wave *wave,
                               const dielectric_pattern_request *pattern,
                               dielectric_diagnostic *diagnostic)
{
  const char *reason = NULL;
  size_t i;
  if(gnd.ksymp == 1)
    return 1;
  if(gnd.iperf != 1 && gnd.iperf != 2)
    reason = "GN 0 dielectric coupling requires a separate Fresnel approximation implementation; use GN 1 or GN 2";
  else if(gnd.iperf == 2 && data.m != 0)
    reason = "GN 2 dielectric coupling currently supports wires only; NEC patch ground uses a different approximation";
  else if(gnd.iperf == 2 && (!isfinite(save.epsr) || save.epsr<1 || !isfinite(save.sig)))
    reason = "GN 2 requires passive finite ground with relative permittivity >= 1";
  else if(gnd.nradl != 0 || gnd.ifar > 0)
    reason = "dielectric ground supports only infinite planar ground without screens or cliff ground";
  else if(iterative == 2 && gnd.iperf != 1)
    reason = "GN 2 planar ground requires -s dense or -s gmres; real-ground FFT is not supported";
  else if(wave != NULL && (!isfinite(wave->theta_degrees) ||
          wave->theta_degrees < 0.0 || wave->theta_degrees > 90.0))
    reason = "ground incident plane wave requires theta between 0 and 90 degrees";
  if(reason == NULL && wave != NULL && (wave->type==2 || wave->type==3) &&
     (!isfinite(wave->axial_ratio) || !(wave->axial_ratio>0.0)))
    reason = "ground elliptical EX requires an explicit positive axial ratio in F6";
  for(i=0; reason == NULL && i<model->box_count; ++i)
    if(!isfinite(model->boxes[i].minimum[2]) ||
       (gnd.iperf==1 ? model->boxes[i].minimum[2]<0.0 : !(model->boxes[i].minimum[2]>0.0)))
      reason = gnd.iperf==1 ? "PEC dielectric boxes require z >= 0, including vacuum boxes" :
        "real ground requires strictly positive clearance for every dielectric box, including vacuum boxes";
  for(i=0; reason == NULL && i<mesh->cell_count; ++i)
  {
    double bottom=mesh->cells[i].center[2]-0.5*mesh->cells[i].size[2];
    if(!isfinite(bottom) || (gnd.iperf==1 ? bottom<0.0 : !(bottom>0.0)))
      reason = gnd.iperf==1 ? "PEC complete volume cells require z >= 0" :
        "real ground requires strictly positive clearance for every complete volume cell";
  }
  for(i=0; reason == NULL && i<(size_t)data.n; ++i)
    if(!(data.z[i]-0.5*data.si[i]*fabs(data.salp[i])-data.bi[i] > 0.0))
      reason = "planar ground requires positive wire endpoint/radius clearance; wire/ground contact is unsupported";
  for(i=0; reason == NULL && i<(size_t)data.m; ++i)
    /* Same effective patch-size scale used by the conductor/cell gate.
       A full size clearance is conservative for the NEC patch abstraction. */
    if(!(data.pz[i]-sqrt(data.pbi[i]) > 0.0))
      reason = "planar ground requires patch centroid height greater than sqrt(area)";
  if(reason == NULL && pattern != NULL && pattern->enabled)
  {
    double first=pattern->theta_start_degrees;
    double last=first+(pattern->theta_count ? pattern->theta_count-1:0)*
                              pattern->theta_step_degrees;
    if(!isfinite(first) || !isfinite(last) || fmin(first,last)<0 || fmax(first,last)>90)
      reason = "ground radiation pattern requires upper-hemisphere theta between 0 and 90 degrees";
  }
  if(reason == NULL && pattern != NULL && pattern->near_enabled)
  {
    int iz;
    double rfirst=fpat.xnr, rlast=rfirst+(fpat.nrx-1)*fpat.dxnr;
    if(fpat.nrx<=0 || fpat.nry<=0 || fpat.nrz<=0)
      reason="ground near-field counts must be positive";
    for(iz=0; reason == NULL && iz<fpat.nrz; ++iz)
    {
      double z=fpat.znr+iz*fpat.dznr;
      if(fpat.near != 0)
      {
        double cosine=cos(z*3.14159265358979323846/180.0);
        z=fmin(rfirst*cosine,rlast*cosine);
      }
      if(!isfinite(z) || z<0.0 || (gnd.iperf==2 && z==0.0))
        reason=gnd.iperf==2 ? "GN 2 near-field observations require z > 0":"ground near-field observations require z >= 0";
    }
  }
  if(reason == NULL)
    return 1;
  memset(diagnostic,0,sizeof(*diagnostic));
  snprintf(diagnostic->resource,sizeof(diagnostic->resource),"%s",
           model->resource_path ? model->resource_path:"<dielectric>");
  snprintf(diagnostic->record,sizeof(diagnostic->record),"GN");
  snprintf(diagnostic->message,sizeof(diagnostic->message),"%s",reason);
  return 0;
}

/* Preserve GN's negative-conductivity convention (specified imaginary epsilon)
   without mutating save.sig across frequency sweeps. Use NEC's constants so
   native Sommerfeld tables and the volume spectrum see identical material. */
complex double dielectric_ground_epsilon(double frequency_mhz)
{
  return save.epsr + I*(save.sig<0 ? save.sig : -save.sig*(CVEL/frequency_mhz)*59.96);
}

void dielectric_ground_initialize(double frequency_mhz)
{
  complex double epsilon;
  if(gnd.ksymp==1) return;
  gnd.frati=1;
  if(gnd.iperf==1) return;
  epsilon=dielectric_ground_epsilon(frequency_mhz);
  gnd.zrati=1/csqrt(epsilon);
  gwav.u=gnd.zrati; gwav.u2=gwav.u*gwav.u;
  if(gnd.iperf==2)
  {
    if(data.n!=0 && epsilon!=1.0) somnec(save.epsr,save.sig,frequency_mhz);
    gnd.frati=(epsilon-1)/(epsilon+1);
  }
}

int dielectric_background_n(const vie_cuboid *a,const vie_cuboid *b,double k,
                           const vie_kernel_options *options,complex double n[3][3])
{
  complex double reflected[3][3],h[3][3];
  int p,q;
  if(gnd.ksymp!=1 && gnd.iperf==1)
    return vie_cuboid_pec_n(a,b,k,0,options,n,NULL);
  if(!vie_cuboid_n(a,b,k,options,n,NULL)) return 0;
  if(gnd.ksymp==1) return 1;
  if(!vie_halfspace_cuboids(a,b,k,dielectric_ground_epsilon(k*299792458.0/(2*3.14159265358979323846*1e6)),
                            NULL,reflected,h,NULL)) return 0;
  for(p=0;p<3;++p) for(q=0;q<3;++q) n[p][q]+=reflected[p][q];
  return 1;
}

/* Reflection of a transverse plane wave; downward incidence or the downward
   source direction associated with an upward far-field ray. */
void dielectric_ground_reflect(const double down[3],const complex double field[3],
                              double k,complex double reflected[3])
{
  double u=hypot(down[0],down[1]),c=u>1e-14 ? down[0]/u:1,s=u>1e-14 ? down[1]/u:0;
  double w=fabs(down[2]),te[3]={-s,c,0},pd[3]={-w*c,-w*s,-u},pu[3]={w*c,w*s,-u};
  complex double rs=-1,rp=1,es=0,ep=0;
  int p;
  if(gnd.iperf!=1)
    vie_halfspace_fresnel(dielectric_ground_epsilon(k*299792458.0/(2*3.14159265358979323846*1e6)),fmin(u,1),&rs,&rp);
  for(p=0;p<3;++p) { es+=te[p]*field[p]; ep+=pd[p]*field[p]; }
  for(p=0;p<3;++p) reflected[p]=rs*te[p]*es+rp*pu[p]*ep;
}

typedef struct {
  const dielectric_mesh_result *mesh;
  const complex double *current;
  int conductor;
} ground_sources;

static double ground_sinc(double x)
{ return fabs(x)<1e-5 ? 1-x*x/6+x*x*x*x/120:sin(x)/x; }

static complex double ground_z_average(double z,double size,complex double kz)
{
  complex double a=kz*size;
  if(cabs(a)<1e-4) return cexp(-I*kz*z)*(1-a*a/24+a*a*a*a/1920);
  return (cexp(-I*kz*(z-.5*size))-cexp(-I*kz*(z+.5*size)))/(I*a);
}

/* Integral of exp(center_phase+i q ell) over a segment. Combine endpoint
   exponents before evaluating them so evanescent factors never overflow. */
static complex double ground_segment_integral(complex double center_phase,complex double q,double length)
{
  complex double a=q*length;
  if(cabs(a)<1e-4) return cexp(center_phase)*length*(1-a*a/24+a*a*a*a/1920);
  return (cexp(center_phase+I*a/2)-cexp(center_phase-I*a/2))/(I*q);
}

static int ground_spectrum(void *context,double k,double u,complex double w,
                           double phi,complex double source[3])
{
  ground_sources *sources=context;
  double c=cos(phi),s=sin(phi);
  size_t cell;
  int p,segment;
  source[0]=source[1]=source[2]=0;
  for(cell=0;cell<sources->mesh->cell_count;++cell)
  {
    const dielectric_cell *a=sources->mesh->cells+cell;
    double volume=a->size[0]*a->size[1]*a->size[2];
    complex double phase=volume*ground_sinc(.5*k*u*c*a->size[0])*ground_sinc(.5*k*u*s*a->size[1])*
      ground_z_average(a->center[2],a->size[2],k*w)*cexp(I*k*u*(c*a->center[0]+s*a->center[1]));
    for(p=0;p<3;++p) source[p]+=phase*sources->current[3*cell+p];
  }
  if(sources->conductor)
    for(segment=0;segment<data.n;++segment)
    {
      double lambda=data.wlam,length=data.si[segment];
      double direction[3]={data.cab[segment],data.sab[segment],data.salp[segment]};
      complex double a=crnt.air[segment]+I*crnt.aii[segment];
      complex double b=crnt.bir[segment]+I*crnt.bii[segment];
      complex double cc=crnt.cir[segment]+I*crnt.cii[segment];
      complex double phase=I*k*lambda*(u*(c*data.x[segment]+s*data.y[segment])-w*data.z[segment]);
      complex double q=k*lambda*(u*(c*direction[0]+s*direction[1])-w*direction[2]);
      complex double base=ground_segment_integral(phase,q,length);
      complex double plus=ground_segment_integral(phase,q+TP,length);
      complex double minus=ground_segment_integral(phase,q-TP,length);
      complex double moment=lambda*lambda*(a*base+b*(plus-minus)/(2*I)+cc*(plus+minus)/2);
      for(p=0;p<3;++p) source[p]+=moment*direction[p];
    }
  if(sources->conductor)
    for(segment=0;segment<data.m;++segment)
    {
      /* Native patches carry Cartesian sheet current (A/m), with the same
         centroid/area approximation used by NEC's far-field transform. */
      double lambda=data.wlam;
      complex double moment=data.pbi[segment]*lambda*lambda*
        cexp(I*k*lambda*(u*(c*data.px[segment]+s*data.py[segment])-w*data.pz[segment]));
      for(p=0;p<3;++p) source[p]+=moment*crnt.cur[data.n+3*segment+p];
    }
  return 1;
}

int dielectric_ground_power(const dielectric_mesh_result *mesh,const complex double *current,
                             double k,int conductor,double *power,double *delta)
{
  ground_sources sources={mesh,current,conductor};
  double clearance=INFINITY;
  size_t cell;
  int segment;
  vie_halfspace_options options;
  for(cell=0;cell<mesh->cell_count;++cell)
    clearance=fmin(clearance,mesh->cells[cell].center[2]-.5*mesh->cells[cell].size[2]);
  if(conductor)
    for(segment=0;segment<data.n;++segment)
      clearance=fmin(clearance,data.wlam*(data.z[segment]-.5*data.si[segment]*fabs(data.salp[segment])));
  if(!isfinite(clearance)) { *power=0; if(delta) *delta=0; return 1; }
  vie_halfspace_defaults(&options); options.relative_tolerance=2e-5; options.absolute_tolerance=1e-16;
  return vie_halfspace_power(k,dielectric_ground_epsilon(k*299792458.0/(TP*1e6)),
                             clearance,ground_spectrum,&sources,&options,power,delta);
}

double dielectric_ground_extinction(const dielectric_mesh_result *mesh,const complex double *current,
                                    double k,const dielectric_plane_wave *wave,int conductor)
{
  ground_sources sources={mesh,current,conductor};
  double direction[3],u,phi,power=0;
  complex double field[3],reflected[3],moment[3];
  int p;
  dielectric_plane_wave_vectors(wave,direction,field);
  memset(reflected,0,sizeof(reflected));
  if(gnd.ksymp!=1) dielectric_ground_reflect(direction,field,k,reflected);
  u=hypot(direction[0],direction[1]); phi=atan2(direction[1],direction[0]);
  ground_spectrum(&sources,k,u,-direction[2],phi,moment);
  for(p=0;p<3;++p) power+=.5*creal(field[p]*conj(moment[p]));
  ground_spectrum(&sources,k,u,direction[2],phi,moment);
  for(p=0;p<3;++p) power+=.5*creal(reflected[p]*conj(moment[p]));
  return power;
}
