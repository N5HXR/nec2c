#include "nec4_ground_matrix.h"
#include "nec4_ground_direct.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define NEC4_C0 299792458.0
#define NEC4_MU0 1.25663706212e-6

static int finite_complex(double complex value)
{ return isfinite(creal(value))&&isfinite(cimag(value)); }

int nec4_ground_matrix_table_cross(void *opaque,double rho,double observer_z,
    double source_depth,double complex field[5])
{
  nec4_ground_matrix_tables *tables=opaque;
  if(!tables||!field) return 0;
  for(int kind=0;kind<3;++kind)
    if(tables->table[kind]&&
       nec4_ground_table_evaluate(tables->table[kind],
         tables->lower_wave_number_per_m,tables->upper_wave_number_per_m,
         rho,observer_z,source_depth,field)) return 1;
  return 0;
}

static int homogeneous(double frequency_hz,double complex wave_number,
    double rho,double source_z,double observer_z,double complex value[5])
{
  double dz=observer_z-source_z,distance=hypot(rho,dz),omega;
  double complex common,radial,transverse,mixed;
  if(!(frequency_hz>0.) || !(distance>0.) || cabs(wave_number)==0.) return 0;
  omega=2.*acos(-1.)*frequency_hz;
  common=-I*omega*NEC4_MU0/(4.*acos(-1.)*wave_number*wave_number)*
    cexp(-I*wave_number*distance)/(distance*distance*distance);
  radial=3.+3.*I*wave_number*distance-wave_number*wave_number*distance*distance;
  transverse=1.+I*wave_number*distance-wave_number*wave_number*distance*distance;
  mixed=rho*dz/(distance*distance)*radial;
  value[NEC4_GROUND_RHO_VERTICAL]=common*mixed;
  value[NEC4_GROUND_Z_VERTICAL]=common*
    (dz*dz/(distance*distance)*radial-transverse);
  value[NEC4_GROUND_RHO_HORIZONTAL]=common*
    (rho*rho/(distance*distance)*radial-transverse);
  value[NEC4_GROUND_PHI_HORIZONTAL]=common*transverse;
  value[NEC4_GROUND_Z_HORIZONTAL]=common*mixed;
  for(int i=0;i<5;++i) if(!finite_complex(value[i])) return 0;
  return 1;
}

int nec4_ground_matrix_homogeneous_components(double frequency,
    double complex wave_number,double rho,double source_z,double observer_z,
    double complex value[5])
{ return homogeneous(frequency,wave_number,rho,source_z,observer_z,value); }

static int canonical(const nec4_ground_matrix_context *context,double rho,
    double observer_z,double source_depth,double complex value[5])
{
  return context&&context->cross&&rho>=0.&&observer_z>=0.&&source_depth>=0.&&
    context->cross(context->cross_context,rho,observer_z,source_depth,value);
}

static int pairing(const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    double rho,double source_z,double observer_z,double complex value[5])
{
  double complex base[5],direct[5],image[5];
  double complex ratio=context->lower_epsilon_relative;
  if(source_medium==NEC4_MEDIUM_LOWER&&observer_medium==NEC4_MEDIUM_UPPER)
    return canonical(context,rho,fabs(observer_z),fabs(source_z),value);
  if(source_medium==NEC4_MEDIUM_UPPER&&observer_medium==NEC4_MEDIUM_LOWER)
  {
    if(!canonical(context,rho,fabs(source_z),fabs(observer_z),base)) return 0;
    value[0]=-base[4];value[1]=base[1];value[2]=base[2];
    value[3]=base[3];value[4]=-base[0];
    return 1;
  }
  if(source_medium==NEC4_MEDIUM_LOWER&&observer_medium==NEC4_MEDIUM_LOWER)
  {
    double a=fabs(source_z),b=fabs(observer_z);
    if(!canonical(context,rho,0.,a+b,base) ||
       !homogeneous(context->frequency_hz,context->lower_wave_number_per_m,
         rho,-a,-b,direct) ||
       !homogeneous(context->frequency_hz,context->lower_wave_number_per_m,
         rho,-a,b,image)) return 0;
    value[0]=base[0]+direct[0]-image[0];
    value[1]=base[1]/ratio+direct[1]-image[1];
    value[2]=base[2]+direct[2]-image[2];
    value[3]=base[3]+direct[3]-image[3];
    value[4]=base[4]/ratio+direct[4]-image[4];
    return 1;
  }
  if(source_medium==NEC4_MEDIUM_UPPER&&observer_medium==NEC4_MEDIUM_UPPER)
  {
    double a=fabs(source_z),b=fabs(observer_z);
    if(!canonical(context,rho,a+b,0.,base)) {
      fprintf(stderr,"nec2c: NEC4 upper pairing canonical term failed\n");
      return 0;
    }
    if(!homogeneous(context->frequency_hz,context->upper_wave_number_per_m,
         rho,a,b,direct)) {
      fprintf(stderr,"nec2c: NEC4 upper pairing direct term failed\n");
      return 0;
    }
    if(!homogeneous(context->frequency_hz,context->upper_wave_number_per_m,
         rho,a,-b,image)) {
      fprintf(stderr,"nec2c: NEC4 upper pairing image term failed\n");
      return 0;
    }
    value[0]=-base[4]+direct[0]-image[0];
    value[1]=ratio*base[1]+direct[1]-image[1];
    value[2]=base[2]+direct[2]-image[2];
    value[3]=base[3]+direct[3]-image[3];
    value[4]=-ratio*base[0]+direct[4]-image[4];
    return 1;
  }
  return 0;
}

static void assemble(const double complex c[5],double cosine,double sine,
    const double direction[3],double complex field[3])
{
  double complex xx=cosine*cosine*c[2]-sine*sine*c[3];
  double complex xy=sine*cosine*(c[2]+c[3]);
  double complex yy=sine*sine*c[2]-cosine*cosine*c[3];
  field[0]=xx*direction[0]+xy*direction[1]+cosine*c[0]*direction[2];
  field[1]=xy*direction[0]+yy*direction[1]+sine*c[0]*direction[2];
  field[2]=cosine*c[4]*direction[0]+sine*c[4]*direction[1]+c[1]*direction[2];
}

int nec4_ground_matrix_point(const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    const double source[3],const double observer[3],
    const double direction[3],double complex field[3])
{
  double dx=observer[0]-source[0],dy=observer[1]-source[1];
  double rho=hypot(dx,dy),cosine=1.,sine=0.;
  double complex components[5];
  if(!context||!source||!observer||!direction||!field) return 0;
  if(rho>0.) { cosine=dx/rho;sine=dy/rho; }
  if(!pairing(context,source_medium,observer_medium,rho,source[2],observer[2],
      components)) return 0;
  assemble(components,cosine,sine,direction,field);
  return finite_complex(field[0])&&finite_complex(field[1])&&finite_complex(field[2]);
}

static int correction_point(const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    const double source[3],const double observer[3],const double direction[3],
    double complex field[3])
{
  double dx=observer[0]-source[0],dy=observer[1]-source[1],rho=hypot(dx,dy);
  double cosine=rho>0.?dx/rho:1.,sine=rho>0.?dy/rho:0.;
  double complex free_components[5],free_field[3];
  if(!nec4_ground_matrix_point(context,source_medium,observer_medium,source,
      observer,direction,field)) {
    fprintf(stderr,"nec2c: NEC4 complete matrix point failed at "
        "rho=%.17g source_z=%.17g observer_z=%.17g\n",
        rho,source[2],observer[2]);
    return 0;
  }
  if(!homogeneous(context->frequency_hz,
       source_medium==NEC4_MEDIUM_LOWER?context->lower_wave_number_per_m:
       context->upper_wave_number_per_m,rho,
       source[2],observer[2],free_components)) {
    fprintf(stderr,"nec2c: NEC4 free-space subtraction failed at "
        "rho=%.17g source_z=%.17g observer_z=%.17g\n",
        rho,source[2],observer[2]);
    return 0;
  }
  assemble(free_components,cosine,sine,direction,free_field);
  for(int i=0;i<3;++i) field[i]-=free_field[i];
  return 1;
}

static int segment_integral(
    const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    const double center[3],const double direction[3],double length,
    const double observer[3],double radius,int subtract_free,
    nec4_ground_segment_field *result)
{
  /* 16-point Gauss-Legendre is deterministic and resolves the smooth ground
     correction.  The singular direct wire term is deliberately not here. */
  static const double x[8]={
    .0950125098376374402,.281603550779258913,.458016777657227386,
    .617876244402643748,.755404408355003034,.865631202387831744,
    .944575023073232576,.989400934991649933};
  static const double w[8]={
    .189450610455068496,.182603415044923589,.169156519395002538,
    .149595988816576733,.124628971255533872,.095158511682492785,
    .062253523938647893,.027152459411754095};
  double shifted[3]={observer[0],observer[1],observer[2]};
  int panels=1;
  if(!context||!center||!direction||!observer||!result||!(length>0.)||radius<0.)
    return 0;
  memset(result,0,sizeof(*result));
  /* Match on the wire surface when the observer lies on the source axis. */
  {
    double vx=observer[0]-center[0],vy=observer[1]-center[1],
      vz=observer[2]-center[2];
    double projection=vx*direction[0]+vy*direction[1]+vz*direction[2];
    double rx=vx-projection*direction[0],ry=vy-projection*direction[1],
      rz=vz-projection*direction[2],r=hypot(hypot(rx,ry),rz);
    if(radius>0.&&r<.5*radius)
    {
      double nx=-direction[1],ny=direction[0],nz=0.,n=hypot(nx,ny);
      if(n<1.e-12) { nx=1.;ny=0.;nz=0.;n=1.; }
      shifted[0]+=radius*nx/n;shifted[1]+=radius*ny/n;shifted[2]+=radius*nz/n;
    }
  }
  /* A point-current representation at the wire radius has a narrow but
     finite near-self peak.  One panel is sufficient for smooth ground
     corrections and separated segments; subdivide only when the observer's
     projection lies on the source and its transverse distance is small. */
  {
    double vx=shifted[0]-center[0],vy=shifted[1]-center[1],
      vz=shifted[2]-center[2];
    double axial=vx*direction[0]+vy*direction[1]+vz*direction[2];
    double tx=vx-axial*direction[0],ty=vy-axial*direction[1],
      tz=vz-axial*direction[2],transverse=hypot(hypot(tx,ty),tz);
    if(fabs(axial)<=.5*length&&transverse<length)
    {
      double scale=fmax(radius,transverse);
      if(scale>0.) panels=(int)ceil(length/(4.*scale));
      if(panels<1) panels=1;
      if(panels>128) panels=128;
    }
  }
  for(int panel=0;panel<panels;++panel)
    for(int q=0;q<8;++q) for(int sign=-1;sign<=1;sign+=2)
  {
    double panel_length=length/panels;
    double middle=-.5*length+(panel+.5)*panel_length;
    double half=.5*panel_length;
    double t=middle+half*sign*x[q],source[3];
    double complex source_wave_number=
      source_medium==NEC4_MEDIUM_LOWER?
      context->lower_wave_number_per_m:context->upper_wave_number_per_m;
    double complex basis_sine=csin(source_wave_number*t);
    double complex basis_cosine=ccos(source_wave_number*t);
    double complex value[3];
    for(int k=0;k<3;++k) source[k]=center[k]+t*direction[k];
    if(subtract_free)
    {
      if(!correction_point(context,source_medium,observer_medium,source,
          shifted,direction,value)) return 0;
    }
    else if(!nec4_ground_matrix_point(context,source_medium,observer_medium,
        source,shifted,direction,value)) return 0;
    for(int k=0;k<3;++k)
    {
      double complex weighted=half*w[q]*value[k];
      result->constant[k]+=weighted;
      result->sine[k]+=weighted*basis_sine;
      result->cosine[k]+=weighted*basis_cosine;
    }
  }
  return 1;
}

int nec4_ground_matrix_segment_correction(
    const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    const double center[3],const double direction[3],double length,
    const double observer[3],double radius,nec4_ground_segment_field *result)
{
  return segment_integral(context,source_medium,observer_medium,center,
    direction,length,observer,radius,1,result);
}

int nec4_ground_matrix_segment_nec_correction(
    const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    const double center[3],const double direction[3],double length,
    const double observer[3],double radius,double normalization_length,
    nec4_ground_segment_field *result)
{
  if(!(normalization_length>0.)||!isfinite(normalization_length)||
     !nec4_ground_matrix_segment_correction(context,source_medium,
       observer_medium,center,direction,length,observer,radius,result)) return 0;
  for(int component=0;component<3;++component)
  {
    result->constant[component]*=normalization_length;
    result->sine[component]*=normalization_length;
    result->cosine[component]*=normalization_length;
  }
  return 1;
}

int nec4_ground_matrix_segment_nec_complete(
    const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    const double center[3],const double direction[3],double length,
    const double observer[3],double radius,double normalization_length,
    nec4_ground_segment_field *result)
{
  if(!(normalization_length>0.)||!isfinite(normalization_length)||
     !segment_integral(context,source_medium,observer_medium,center,direction,
       length,observer,radius,0,result)) return 0;
  for(int component=0;component<3;++component)
  {
    result->constant[component]*=normalization_length;
    result->sine[component]*=normalization_length;
    result->cosine[component]*=normalization_length;
  }
  return 1;
}

int nec4_ground_matrix_endpoint_charge_nec(
    const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,const double endpoint[3],
    const double observer[3],double normalization_length,int end_sign,
    double complex field[3])
{
  const double eps0=8.8541878128e-12;
  double dx,dy,dz,r,omega;
  double complex epsilon,k,common;
  if(!context||!endpoint||!observer||!field||!(normalization_length>0.)||
     (end_sign!=1&&end_sign!=-1)||!(context->frequency_hz>0.)) return 0;
  dx=observer[0]-endpoint[0];dy=observer[1]-endpoint[1];
  dz=observer[2]-endpoint[2];r=hypot(hypot(dx,dy),dz);
  if(!(r>0.)) return 0;
  omega=2.*acos(-1.)*context->frequency_hz;
  if(source_medium==NEC4_MEDIUM_LOWER)
  { epsilon=context->lower_epsilon_relative;k=context->lower_wave_number_per_m; }
  else { epsilon=1.;k=context->upper_wave_number_per_m; }
  common=end_sign*I*normalization_length/omega/(4.*acos(-1.)*eps0*epsilon)*
    cexp(-I*k*r)*(1.+I*k*r)/(r*r*r);
  field[0]=common*dx;field[1]=common*dy;field[2]=common*dz;
  return finite_complex(field[0])&&finite_complex(field[1])&&
    finite_complex(field[2]);
}
