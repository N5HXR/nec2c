#include "nec2c.h"
#include "shared.h"
#include "wire_operator.h"
#include "junction_charge.h"
#include "nec4_endcap.h"

#include <float.h>

typedef struct
{
  const char *name;
  void (*sbf)(int, int, double complex *, double complex *, double complex *);
  void (*tbf)(int, int);
  void (*trio)(int);
  void (*efld)(double, double, double, double, int);
  void (*hsfld)(double, double, double, double);
} wire_operator_vtable;

static const wire_operator_vtable nec2_operator = {
  "nec2", sbf, tbf, trio, efld, hsfld
};

static const wire_operator_vtable nec4_operator = {
  "nec4ish", sbf, tbf, trio, efld, hsfld
};

static const wire_operator_vtable *active_operator = &nec2_operator;
static const char *requested_engine = "nec2";
static const char *requested_ground_engine = "nec2";
static const char *resolved_ground_engine = "nec2";
static const char *selection_reason = "not-resolved";
static const char *model_check_policy = "warn";
static wire_operator_validity last_validity;
static size_t model_patches, model_coated_domains;
static int model_dielectric_loaded;
static const nec4_medium_state *medium_state;
static _Thread_local int ground_direct_only;
static const nec4_ground_matrix_context *ground_matrix_context;
static _Thread_local size_t matrix_source_segment=(size_t)-1;

static double dot3(double ax,double ay,double az,
                   double bx,double by,double bz)
{ return ax*bx+ay*by+az*bz; }

static int endcap_disk_global_image(size_t segment,int end,int image,
    double xi,double yi,double zi,double complex *ex,
    double complex *ey,double complex *ez)
{
  double cx=end?data.x2[segment]:data.x1[segment];
  double cy=end?data.y2[segment]:data.y1[segment];
  double cz=end?data.z2[segment]:data.z1[segment];
  double cab=data.cab[segment],sab=data.sab[segment],salp=data.salp[segment];
  if(image) { cz=-cz;salp=-salp; }
  double dx=xi-cx,dy=yi-cy,dz=zi-cz;
  double axial=dot3(dx,dy,dz,cab,sab,salp);
  double rx=dx-axial*cab;
  double ry=dy-axial*sab;
  double rz=dz-axial*salp;
  double radial=hypot(hypot(rx,ry),rz);
  nec4_endcap_field field;
  if(!nec4_endcap_disk_field(data.bi[segment],radial,axial,
      wire_operator_segment_wave_number(segment),&field))
    return 0;
  if(radial>0.)
  {
    rx/=radial;ry/=radial;rz/=radial;
  }
  else rx=ry=rz=0.;
  *ex=field.axial*cab+field.radial*rx;
  *ey=field.axial*sab+field.radial*ry;
  *ez=field.axial*salp+field.radial*rz;
  return 1;
}

static void reflect_cap_field(size_t segment,int end,double xi,double yi,
    double zi,double complex *ex,double complex *ey,double complex *ez)
{
  double cx=end?data.x2[segment]:data.x1[segment];
  double cy=end?data.y2[segment]:data.y1[segment];
  double cz=end?data.z2[segment]:data.z1[segment];
  double xij=xi-cx,yij=yi-cy,zij=zi+cz;
  double r=hypot(hypot(xij,yij),zij),xymag=hypot(xij,yij);
  double px,py,cth;
  complex double refs,refps,zrsin,zratx=gnd.zrati,epy,epx;
  if(gnd.iperf>0) return;
  /* This is the same Fresnel dyadic used by efld() for GN 0. */
  if(xymag<=1.e-6) { px=py=0.;cth=1.;zrsin=CPLX_10; }
  else
  {
    px=-yij/xymag;py=xij/xymag;cth=zij/r;
    zrsin=csqrt(1.-zratx*zratx*(1.-cth*cth));
  }
  refs=(cth-zratx*zrsin)/(cth+zratx*zrsin);
  refps=-(zratx*cth-zrsin)/(zratx*cth+zrsin)-refs;
  epy=px*(*ex)+py*(*ey);epx=px*epy;epy=py*epy;
  *ex=refs*(*ex)+refps*epx;
  *ey=refs*(*ey)+refps*epy;
  *ez=refs*(*ez);
}

static int endcap_disk_global(size_t segment,int end,double xi,double yi,
                              double zi,double complex *ex,
                              double complex *ey,double complex *ez)
{
  complex double ix,iy,iz;
  if(!endcap_disk_global_image(segment,end,0,xi,yi,zi,ex,ey,ez)) return 0;
  if(gnd.ksymp<=1) return 1;
  if(!endcap_disk_global_image(segment,end,1,xi,yi,zi,&ix,&iy,&iz)) return 0;
  reflect_cap_field(segment,end,xi,yi,zi,&ix,&iy,&iz);
  *ex-=gnd.frati*ix;*ey-=gnd.frati*iy;*ez-=gnd.frati*iz;
  return 1;
}

int wire_operator_endcap_current_field(size_t segment,int end,
    double xi,double yi,double zi,double complex *ex,
    double complex *ey,double complex *ez)
{
  double complex scale,k,eta=ETA;
  int lower=0;
  if(!wire_operator_nec4_basis() || segment>=(size_t)data.n ||
     (end!=0 && end!=1)) return 0;
  if(!endcap_disk_global(segment,end,xi,yi,zi,ex,ey,ez)) return 0;
  /* Equations (4-35), (4-37), generalized to the segment medium.  With
     lengths normalized by lambda0 the free-space expression reduces to
     I*eta0/(8*pi^3*a^2). */
  k=wire_operator_segment_wave_number(segment);
  if(wire_operator_nec4_ground()&&medium_state&&
     segment<medium_state->segment_count&&
     medium_state->segment_medium[segment]==NEC4_MEDIUM_LOWER)
  { eta=medium_state->intrinsic_impedance_ohm;lower=1; }
  if(lower)
    scale=I*(end?-1.:1.)*eta/
      (4.*PI*PI*k*data.bi[segment]*data.bi[segment]);
  else
    /* Keep the accepted free-space M5 arithmetic bit-for-bit. */
    scale=I*(end?-1.:1.)*ETA/
      (8.*PI*PI*PI*data.bi[segment]*data.bi[segment]);
  *ex*=scale;*ey*=scale;*ez*=scale;
  return 1;
}

int wire_operator_gap_voltage_field(size_t segment,double xi,double yi,
    double zi,double complex *ex,double complex *ey,double complex *ez)
{
  double complex nx,ny,nz,px,py,pz;
  double scale;
  if(!wire_operator_nec4_basis() || segment>=(size_t)data.n) return 0;
  if(!endcap_disk_global(segment,0,xi,yi,zi,&nx,&ny,&nz) ||
     !endcap_disk_global(segment,1,xi,yi,zi,&px,&py,&pz)) return 0;
  /* Equations (4-39), (4-40): positive cap at end 2, negative at end 1. */
  scale=1./(4.*PI*data.si[segment]);
  *ex=(px-nx)*scale;*ey=(py-ny)*scale;*ez=(pz-nz)*scale;
  return 1;
}

static int same_point(double ax,double ay,double az,
                      double bx,double by,double bz)
{
  double scale=1.+fmax(fmax(fabs(ax),fabs(ay)),fabs(az));
  return hypot(hypot(ax-bx,ay-by),az-bz)<=64.*DBL_EPSILON*scale;
}

static int connected_pair(size_t i,size_t j)
{
  return same_point(data.x1[i],data.y1[i],data.z1[i],
                    data.x1[j],data.y1[j],data.z1[j]) ||
    same_point(data.x1[i],data.y1[i],data.z1[i],
               data.x2[j],data.y2[j],data.z2[j]) ||
    same_point(data.x2[i],data.y2[i],data.z2[i],
               data.x1[j],data.y1[j],data.z1[j]) ||
    same_point(data.x2[i],data.y2[i],data.z2[i],
               data.x2[j],data.y2[j],data.z2[j]);
}

static double point_segment_distance(size_t observer,size_t source)
{
  double dx=data.x[observer]-data.x[source];
  double dy=data.y[observer]-data.y[source];
  double dz=data.z[observer]-data.z[source];
  double axial=dot3(dx,dy,dz,data.cab[source],data.sab[source],
                    data.salp[source]);
  double half=.5*data.si[source];
  if(axial>half) axial=half;
  else if(axial<-half) axial=-half;
  dx-=axial*data.cab[source];
  dy-=axial*data.sab[source];
  dz-=axial*data.salp[source];
  return hypot(hypot(dx,dy),dz);
}

static int prepare_validity(FILE *report)
{
  size_t i,j;
  memset(&last_validity,0,sizeof(last_validity));
  last_validity.minimum_length_radius_ratio=DBL_MAX;
  for(i=0;i<(size_t)data.n;++i)
  {
    double ratio=data.si[i]/data.bi[i];
    if(ratio<last_validity.minimum_length_radius_ratio)
      last_validity.minimum_length_radius_ratio=ratio;
    if(ratio<2.) ++last_validity.short_thick_segments;
    for(j=i+1;j<(size_t)data.n;++j)
    {
      double alignment=fabs(dot3(data.cab[i],data.sab[i],data.salp[i],
                                 data.cab[j],data.sab[j],data.salp[j]));
      double dij=point_segment_distance(i,j),dji=point_segment_distance(j,i);
      int connected=connected_pair(i,j);
      if(connected && fabs(data.bi[i]/data.bi[j]-1.)>1.e-6)
        ++last_validity.radius_step_junctions;
      if(!connected && alignment>.95)
      {
        double dx=data.x[j]-data.x[i],dy=data.y[j]-data.y[i],
          dz=data.z[j]-data.z[i];
        double axial=fabs(dot3(dx,dy,dz,data.cab[i],data.sab[i],data.salp[i]));
        double transverse=sqrt(fmax(0.,dx*dx+dy*dy+dz*dz-axial*axial));
        if(axial<.5*(data.si[i]+data.si[j]) &&
           transverse<4.*fmax(data.bi[i],data.bi[j]))
          ++last_validity.close_parallel_pairs;
      }
      if(!connected && (dij<data.bi[j] || dji<data.bi[i]))
        ++last_validity.buried_match_points;
      if(connected && alignment<.999999 &&
         (dij<data.bi[j] || dji<data.bi[i]))
        ++last_validity.overlapping_bends;
    }
  }
  for(i=0;i<(size_t)vsorc.nsant;++i)
  {
    size_t source=(size_t)(vsorc.isant[i]-1),j;
    if(source>=(size_t)data.n) continue;
    for(j=0;j<(size_t)data.n;++j)
      if(j!=source && connected_pair(source,j))
      {
        double ratio=data.si[source]/data.si[j];
        if(ratio<1.) ratio=1./ratio;
        if(ratio>2.) {++last_validity.source_length_mismatches;break;}
      }
  }
  last_validity.legacy_discontinuity_sources=(size_t)vsorc.nvqd;
  if(last_validity.minimum_length_radius_ratio==DBL_MAX)
    last_validity.minimum_length_radius_ratio=0.;
  if(report)
    fprintf(report,"\n  NEC4ISH WIRE VALIDITY: %lu SEGMENTS BELOW L/A=2, "
            "%lu CLOSE PARALLEL PAIRS, %lu BURIED MATCH POINTS, "
            "%lu OVERLAPPING BENDS, %lu SOURCE LENGTH MISMATCHES, "
            "%lu RADIUS-STEP JUNCTIONS, %lu LEGACY EX5 SOURCES; "
            "MINIMUM L/A %.6G\n",
            (unsigned long)last_validity.short_thick_segments,
            (unsigned long)last_validity.close_parallel_pairs,
            (unsigned long)last_validity.buried_match_points,
            (unsigned long)last_validity.overlapping_bends,
            (unsigned long)last_validity.source_length_mismatches,
            (unsigned long)last_validity.radius_step_junctions,
            (unsigned long)last_validity.legacy_discontinuity_sources,
            last_validity.minimum_length_radius_ratio);
  return !(last_validity.short_thick_segments ||
           last_validity.close_parallel_pairs ||
           last_validity.buried_match_points ||
           last_validity.overlapping_bends ||
           last_validity.source_length_mismatches);
}

int wire_operator_set_engine(const char *name)
{
  if(!name) return 0;
  if(!strcmp(name, "auto"))
  {
    requested_engine = "auto";
    active_operator = &nec2_operator;
    return 1;
  }
  if(!strcmp(name, "nec2"))
  {
    requested_engine = "nec2";
    active_operator = &nec2_operator;
    return 1;
  }
  if(!strcmp(name, "nec4ish")||!strcmp(name, "nec4"))
  {
    requested_engine = "nec4ish";
    active_operator = &nec4_operator;
    return 1;
  }
  return 0;
}

int wire_operator_set_ground_engine(const char *name)
{
  if(!name) return 0;
  if(!strcmp(name, "auto"))
  {
    requested_ground_engine = "auto";
    resolved_ground_engine = "nec2";
    return 1;
  }
  if(!strcmp(name, "nec2"))
  {
    requested_ground_engine = "nec2";
    resolved_ground_engine = "nec2";
    return 1;
  }
  if(!strcmp(name, "nec4ish")||!strcmp(name, "nec4"))
  {
    requested_ground_engine = "nec4ish";
    resolved_ground_engine = "nec4ish";
    return 1;
  }
  return 0;
}

int wire_operator_resolve_auto(void)
{
  int automatic=!strcmp(requested_engine,"auto")||
    !strcmp(requested_ground_engine,"auto");
  int wire_unsupported=model_coated_domains||model_dielectric_loaded;
  int finite_sommerfeld=gnd.ksymp==2&&gnd.iperf==2;
  int interface_wire=0;
  for(int segment=0;segment<data.n&&!interface_wire;++segment)
    interface_wire=data.z1[segment]<0.||data.z2[segment]<0.;
  int ground_unsupported=wire_unsupported||model_patches||!finite_sommerfeld||
    !interface_wire;

  if(!strcmp(requested_engine,"auto"))
  {
    /* An explicit NEC4ish ground request also requests its required wire
       basis.  Unsupported combinations then fail with the ordinary,
       actionable compatibility diagnostic rather than silently changing an
       explicit choice. */
    if(!strcmp(requested_ground_engine,"nec4ish")||!wire_unsupported)
      active_operator=&nec4_operator;
    else
      active_operator=&nec2_operator;
  }

  if(!strcmp(requested_ground_engine,"auto"))
    resolved_ground_engine=active_operator==&nec4_operator&&
      !ground_unsupported?"nec4ish":"nec2";

  if(!automatic)
    selection_reason="explicit-selection";
  else if(wire_unsupported&&active_operator==&nec2_operator)
    selection_reason=model_coated_domains&&model_dielectric_loaded?
      "coated-wire-and-dielectric-vie":
      (model_coated_domains?"coated-wire":"dielectric-vie");
  else if(model_patches&&!strcmp(resolved_ground_engine,"nec2")&&
          !strcmp(requested_ground_engine,"auto"))
    selection_reason="surface-patch-ground-fallback";
  else if(finite_sommerfeld&&!interface_wire&&
          !strcmp(resolved_ground_engine,"nec2")&&
          !strcmp(requested_ground_engine,"auto"))
    selection_reason="above-ground-optimized-ground";
  else if(!finite_sommerfeld&&!strcmp(resolved_ground_engine,"nec2")&&
          !strcmp(requested_ground_engine,"auto"))
    selection_reason="classic-ground-mode";
  else
    selection_reason="supported-nec4ish-model";
  return 1;
}

int wire_operator_set_model_check(const char *name)
{
  if(!name) return 0;
  if(strcmp(name, "warn") && strcmp(name, "strict") && strcmp(name, "off"))
    return 0;
  model_check_policy = !strcmp(name, "strict") ? "strict" :
    (!strcmp(name, "off") ? "off" : "warn");
  return 1;
}

void wire_operator_set_model_features(size_t patches,size_t coated_domains,
                                      int dielectric_loaded)
{
  model_patches=patches;
  model_coated_domains=coated_domains;
  model_dielectric_loaded=dielectric_loaded!=0;
}

void wire_operator_set_medium_state(const nec4_medium_state *state)
{ medium_state=state; }

void wire_operator_set_ground_matrix_context(
    const nec4_ground_matrix_context *context)
{ ground_matrix_context=context; }

void wire_operator_set_source_segment(size_t segment)
{ matrix_source_segment=segment; }

int wire_operator_ground_direct_only(void)
{ return ground_direct_only; }

const char *wire_operator_requested_engine(void) { return requested_engine; }
const char *wire_operator_resolved_engine(void) { return active_operator->name; }
const char *wire_operator_requested_ground_engine(void)
{
  return requested_ground_engine;
}
const char *wire_operator_resolved_ground_engine(void)
{
  return resolved_ground_engine;
}
const char *wire_operator_selection_reason(void) { return selection_reason; }
const char *wire_operator_model_check(void) { return model_check_policy; }
int wire_operator_nec4_basis(void) { return active_operator == &nec4_operator; }
int wire_operator_nec4_ground(void)
{ return !strcmp(resolved_ground_engine,"nec4ish"); }

int wire_operator_prepare_frequency(FILE *report)
{
  if(!wire_operator_nec4_basis()) return 1;
  {
    int valid=prepare_validity(!strcmp(model_check_policy,"off")?NULL:report);
    if(model_coated_domains || model_dielectric_loaded)
    {
      if(report) fprintf(report,"\n  NEC4ISH WIRE COMBINATION ERROR: %s%s%s has not "
          "been qualified with the NEC4ish wire operator; select --wire-engine nec2.\n",
          model_coated_domains?"LD 7 coated wire":"",
          model_coated_domains&&model_dielectric_loaded?" and ":"",
          model_dielectric_loaded?"dielectric VIE":"");
      return 0;
    }
    if(report && model_patches)
      fprintf(report,"  NEC4ISH WIRE COMBINATION: %lu surface patches enabled.\n",
              (unsigned long)model_patches);
    if(report && gnd.ksymp>1 && gnd.iperf==2)
      fprintf(report,"  NEC4ISH GROUND NOTICE: flat-cap direct and image fields are "
          "included; the cap-specific Sommerfeld correction remains provisional.\n");
    if(report && gnd.nradl)
      fprintf(report,"  NEC4ISH GROUND NOTICE: radial-screen reflection of flat-cap "
          "fields uses the homogeneous-ground coefficient.\n");
    if(!junction_charge_prepare(report,!strcmp(model_check_policy,"strict"),
        medium_state?medium_state->endpoint_factor_assigned:NULL))
      return 0;
    return strcmp(model_check_policy,"strict") || valid;
  }
}

const wire_operator_validity *wire_operator_last_validity(void)
{ return &last_validity; }

double complex wire_operator_charge_factor(size_t segment,int end,
                                           double complex fallback)
{
  double complex interface_factor;
  if(!wire_operator_nec4_basis()) return fallback;
  if(medium_state&&nec4_medium_endpoint_factor(medium_state,segment,end,
      &interface_factor)) return interface_factor;
  return junction_charge_factor(segment,end,fallback);
}

double complex wire_operator_segment_wave_number(size_t segment)
{
  if(wire_operator_nec4_ground()&&medium_state&&
     segment<medium_state->segment_count&&
     medium_state->segment_medium[segment]==NEC4_MEDIUM_LOWER)
    /* Geometry is normalized with NEC's historical free-space wavelength.
       Use the exact medium/free-space ratio here, rather than multiplying a
       modern-SI k by that historical wavelength.  Besides avoiding a small
       mixed-constant scale error, this makes an er=1 invisible interface an
       exact algebraic null. */
    return TP*medium_state->refractive_index;
  return TP;
}

double complex wire_operator_source_wave_number(void)
{ return wire_operator_segment_wave_number(matrix_source_segment); }

double complex wire_operator_source_intrinsic_impedance(void)
{
  if(wire_operator_nec4_ground()&&medium_state&&
     matrix_source_segment<medium_state->segment_count&&
     medium_state->segment_medium[matrix_source_segment]==NEC4_MEDIUM_LOWER)
    return ETA/medium_state->refractive_index;
  return ETA;
}

double complex wire_operator_segment_phase(size_t segment,
    double offset_wavelengths)
{ return wire_operator_segment_wave_number(segment)*offset_wavelengths; }

double complex wire_operator_end_parameter(size_t segment,double radius)
{
  double complex ka=.5*wire_operator_segment_wave_number(segment)*radius;
  if(wire_operator_nec4_basis()) return ka;
  return ka*(1.-.5*ka*ka)/(1.-ka*ka);
}

double complex wire_operator_cosine_shape(double complex angle)
{
  if(!wire_operator_nec4_basis()) return cos(creal(angle));
  if(cimag(angle)==0.)
  {
    double sine=sin(.5*creal(angle));
    return -2.*sine*sine;
  }
  return -2.0*csin(0.5*angle)*csin(0.5*angle);
}

double complex wire_operator_center_basis(double complex constant,
                                          double complex cosine)
{
  return wire_operator_nec4_basis() ? constant : constant+cosine;
}

double wire_operator_cosine_integral(double constant_integral,
                                     double cosine_integral)
{
  return wire_operator_nec4_basis() ? cosine_integral-constant_integral :
    cosine_integral;
}

double complex wire_operator_current(double complex constant,
                                     double complex sine,
                                     double complex cosine,
                                     size_t segment,double offset_wavelengths)
{
  double complex angle=wire_operator_segment_phase(segment,offset_wavelengths);
  double complex sine_shape=cimag(angle)==0.?sin(creal(angle)):csin(angle);
  return constant+sine*sine_shape+cosine*wire_operator_cosine_shape(angle);
}

void wire_operator_sbf(int i, int is, double complex *aa,
                       double complex *bb, double complex *cc)
{
  active_operator->sbf(i, is, aa, bb, cc);
}

void wire_operator_tbf(int i, int icap)
{
  active_operator->tbf(i, icap);
}

void wire_operator_trio(int j)
{
  active_operator->trio(j);
}

static void wire_operator_efld_impl(double xi,double yi,double zi,double ai,
                                   int ij,int point_observer)
{
  if(wire_operator_nec4_ground()&&ground_matrix_context&&medium_state&&
     matrix_source_segment<medium_state->segment_count)
  {
    size_t observer_segment=(size_t)((long)matrix_source_segment+(long)ij);
    nec4_segment_medium source_medium=
      medium_state->segment_medium[matrix_source_segment];
    nec4_segment_medium observer_medium=point_observer?
      (zi*data.wlam>=-medium_state->interface_tolerance_m?
        NEC4_MEDIUM_UPPER:NEC4_MEDIUM_LOWER):
      (observer_segment<medium_state->segment_count?
        medium_state->segment_medium[observer_segment]:NEC4_MEDIUM_UPPER);
    nec4_ground_segment_field correction;
    double center[3]={dataj.xj*data.wlam,dataj.yj*data.wlam,
      dataj.zj*data.wlam};
    double direction[3]={dataj.cabj,dataj.sabj,dataj.salpj};
    double observer[3]={xi*data.wlam,yi*data.wlam,zi*data.wlam};
    /* Matrix matching retains the reduced finite-radius direct kernel and
       adds complete-interface minus homogeneous-medium field.  For an
       arbitrary field point, use the complete interface field directly:
       subtracting one far-field implementation and adding another loses the
       grazing direct/image cancellation after many wavelengths. */
    if(point_observer)
    {
      if(!nec4_ground_matrix_segment_nec_complete(ground_matrix_context,
          source_medium,observer_medium,center,direction,dataj.s*data.wlam,
          observer,ai*data.wlam,data.wlam,&correction))
      {
        fprintf(stderr,"nec2c: NEC4ish complete point-field evaluation failed "
            "(source=%lu source_xyz=%.17g,%.17g,%.17g "
            "observer_xyz=%.17g,%.17g,%.17g)\n",
            (unsigned long)matrix_source_segment,
            center[0],center[1],center[2],observer[0],observer[1],observer[2]);
        stop(-1);
      }
      dataj.exk=correction.constant[0];
      dataj.eyk=correction.constant[1];
      dataj.ezk=correction.constant[2];
      dataj.exs=correction.sine[0];
      dataj.eys=correction.sine[1];
      dataj.ezs=correction.sine[2];
      dataj.exc=correction.cosine[0];
      dataj.eyc=correction.cosine[1];
      dataj.ezc=correction.cosine[2];
      if(wire_operator_nec4_basis())
      {
        dataj.exc-=dataj.exk;
        dataj.eyc-=dataj.eyk;
        dataj.ezc-=dataj.ezk;
      }
      return;
    }
    ground_direct_only=1;
    active_operator->efld(xi,yi,zi,ai,ij);
    ground_direct_only=0;
    if(observer_segment>=medium_state->segment_count||
       !nec4_ground_matrix_segment_nec_correction(ground_matrix_context,
         source_medium,observer_medium,center,direction,
         dataj.s*data.wlam,observer,ai*data.wlam,data.wlam,&correction))
    {
      fprintf(stderr,"nec2c: NEC4ish ground matrix evaluation failed "
          "(source=%lu observer=%lu source_xyz=%.17g,%.17g,%.17g "
          "observer_xyz=%.17g,%.17g,%.17g)\n",
          (unsigned long)matrix_source_segment,
          point_observer?(unsigned long)medium_state->segment_count:
            (unsigned long)observer_segment,
          center[0],center[1],center[2],
          observer[0],observer[1],observer[2]);
      stop(-1);
    }
    /* The matrix adapter has converted physical E/I in ohm/m to NEC's ohm
       normalization using data.wlam. */
    dataj.exk+=correction.constant[0];
    dataj.eyk+=correction.constant[1];
    dataj.ezk+=correction.constant[2];
    dataj.exs+=correction.sine[0];
    dataj.eys+=correction.sine[1];
    dataj.ezs+=correction.sine[2];
    dataj.exc+=correction.cosine[0];
    dataj.eyc+=correction.cosine[1];
    dataj.ezc+=correction.cosine[2];
    /* Reduced segment fields omit artificial endpoint charges because they
       cancel at ordinary same-medium junctions.  At an interface crossing
       the two charges occupy different media and section 4.1 requires both
       endpoint fields to be restored.  The interface correction already
       contains complete-minus-homogeneous endpoint fields; add the omitted
       homogeneous part here. */
    for(int end=0;end<2;++end)
    {
      size_t endpoint_index=2*matrix_source_segment+(size_t)end;
      if(medium_state->endpoint_factor_assigned[endpoint_index])
      {
        double endpoint[3]={
          (end?data.x2[matrix_source_segment]:data.x1[matrix_source_segment])*
            data.wlam,
          (end?data.y2[matrix_source_segment]:data.y1[matrix_source_segment])*
            data.wlam,
          (end?data.z2[matrix_source_segment]:data.z1[matrix_source_segment])*
            data.wlam};
        double complex charge[3];
        double offset=(end?1.:-1.)*.5*dataj.s;
        double complex phase=wire_operator_segment_phase(
          matrix_source_segment,offset);
        if(!nec4_ground_matrix_endpoint_charge_nec(ground_matrix_context,
            source_medium,endpoint,observer,data.wlam,end? -1:1,charge))
          stop(-1);
        dataj.exk+=charge[0];dataj.eyk+=charge[1];dataj.ezk+=charge[2];
        dataj.exs+=charge[0]*csin(phase);
        dataj.eys+=charge[1]*csin(phase);
        dataj.ezs+=charge[2]*csin(phase);
        dataj.exc+=charge[0]*ccos(phase);
        dataj.eyc+=charge[1]*ccos(phase);
        dataj.ezc+=charge[2]*ccos(phase);
      }
    }
  }
  else active_operator->efld(xi, yi, zi, ai, ij);
  if(wire_operator_nec4_basis())
  {
    dataj.exc -= dataj.exk;
    dataj.eyc -= dataj.eyk;
    dataj.ezc -= dataj.ezk;
  }
}

void wire_operator_efld(double xi,double yi,double zi,double ai,int ij)
{ wire_operator_efld_impl(xi,yi,zi,ai,ij,0); }

void wire_operator_efld_point(double xi,double yi,double zi,double ai)
{
  /* nefld historically evaluates an arbitrary observation point with ij=1;
     retain that basis/end-condition selector for every wire engine. */
  wire_operator_efld_impl(xi,yi,zi,ai,1,1);
}

void wire_operator_hsfld(double xi, double yi, double zi, double ai)
{
  if(wire_operator_nec4_ground()&&ground_matrix_context&&
     ground_matrix_context->cross_h&&medium_state&&
     matrix_source_segment<medium_state->segment_count&&
     medium_state->segment_medium[matrix_source_segment]==NEC4_MEDIUM_LOWER&&
     zi*data.wlam>=-medium_state->interface_tolerance_m)
  {
    nec4_ground_segment_field complete;
    double center[3]={dataj.xj*data.wlam,dataj.yj*data.wlam,
      dataj.zj*data.wlam};
    double direction[3]={dataj.cabj,dataj.sabj,dataj.salpj};
    double observer[3]={xi*data.wlam,yi*data.wlam,zi*data.wlam};
    (void)ai;
    if(!nec4_ground_matrix_h_segment_nec_cross(ground_matrix_context,center,
        direction,dataj.s*data.wlam,observer,data.wlam,&complete))
    {
      fprintf(stderr,"nec2c: NEC4ish direct H evaluation failed "
          "(source=%lu observer_xyz=%.17g,%.17g,%.17g)\n",
          (unsigned long)matrix_source_segment,observer[0],observer[1],
          observer[2]);
      stop(-1);
    }
    dataj.exk=complete.constant[0];dataj.eyk=complete.constant[1];
    dataj.ezk=complete.constant[2];dataj.exs=complete.sine[0];
    dataj.eys=complete.sine[1];dataj.ezs=complete.sine[2];
    dataj.exc=complete.cosine[0];dataj.eyc=complete.cosine[1];
    dataj.ezc=complete.cosine[2];
  }
  else active_operator->hsfld(xi, yi, zi, ai);
  if(wire_operator_nec4_basis())
  {
    dataj.exc -= dataj.exk;
    dataj.eyc -= dataj.eyk;
    dataj.ezc -= dataj.ezk;
  }
}
