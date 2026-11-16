#include "nec2c.h"
#include "shared.h"
#include "wire_operator.h"
#include "junction_charge.h"

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
  "nec4", sbf, tbf, trio, efld, hsfld
};

static const wire_operator_vtable *active_operator = &nec2_operator;
static const char *requested_engine = "auto";
static const char *requested_ground_engine = "auto";
static const char *resolved_ground_engine = "nec2";
static const char *model_check_policy = "warn";
static wire_operator_validity last_validity;

static double dot3(double ax,double ay,double az,
                   double bx,double by,double bz)
{ return ax*bx+ay*by+az*bz; }

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
  if(last_validity.minimum_length_radius_ratio==DBL_MAX)
    last_validity.minimum_length_radius_ratio=0.;
  if(report)
    fprintf(report,"\n  NEC4 WIRE VALIDITY: %lu SEGMENTS BELOW L/A=2, "
            "%lu CLOSE PARALLEL PAIRS, %lu BURIED MATCH POINTS, "
            "%lu OVERLAPPING BENDS; MINIMUM L/A %.6G\n",
            (unsigned long)last_validity.short_thick_segments,
            (unsigned long)last_validity.close_parallel_pairs,
            (unsigned long)last_validity.buried_match_points,
            (unsigned long)last_validity.overlapping_bends,
            last_validity.minimum_length_radius_ratio);
  return !(last_validity.short_thick_segments ||
           last_validity.close_parallel_pairs ||
           last_validity.buried_match_points ||
           last_validity.overlapping_bends);
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
  if(!strcmp(name, "nec4"))
  {
    requested_engine = "nec4";
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
  return 0;
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
const char *wire_operator_model_check(void) { return model_check_policy; }
int wire_operator_nec4_basis(void) { return active_operator == &nec4_operator; }

int wire_operator_prepare_frequency(FILE *report)
{
  if(!wire_operator_nec4_basis()) return 1;
  {
    int valid=prepare_validity(!strcmp(model_check_policy,"off")?NULL:report);
    if(!junction_charge_prepare(report,!strcmp(model_check_policy,"strict")))
      return 0;
    return strcmp(model_check_policy,"strict") || valid;
  }
}

const wire_operator_validity *wire_operator_last_validity(void)
{ return &last_validity; }

double complex wire_operator_charge_factor(size_t segment,int end,
                                           double complex fallback)
{
  if(!wire_operator_nec4_basis()) return fallback;
  return junction_charge_factor(segment,end,fallback);
}

double wire_operator_cosine_shape(double angle)
{
  if(!wire_operator_nec4_basis()) return cos(angle);
  return -2.0*sin(0.5*angle)*sin(0.5*angle);
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
                                     double angle)
{
  return constant+sine*sin(angle)+cosine*wire_operator_cosine_shape(angle);
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

void wire_operator_efld(double xi, double yi, double zi, double ai, int ij)
{
  active_operator->efld(xi, yi, zi, ai, ij);
  if(wire_operator_nec4_basis())
  {
    dataj.exc -= dataj.exk;
    dataj.eyc -= dataj.eyk;
    dataj.ezc -= dataj.ezk;
  }
}

void wire_operator_hsfld(double xi, double yi, double zi, double ai)
{
  active_operator->hsfld(xi, yi, zi, ai);
  if(wire_operator_nec4_basis())
  {
    dataj.exc -= dataj.exk;
    dataj.eyc -= dataj.eyk;
    dataj.ezc -= dataj.ezk;
  }
}
