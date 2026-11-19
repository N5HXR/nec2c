#include "nec4_medium.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NEC4_C0 299792458.0
#define NEC4_ETA0 376.730313668
#define NEC4_EPS0 8.8541878128e-12

static int side(double z,double tolerance)
{ return z>tolerance?1:(z< -tolerance?-1:0); }

static double endpoint(const nec4_medium_geometry *g,size_t segment,int end,int axis)
{
  const double *v;
  if(axis==0) v=end?g->x2:g->x1;
  else if(axis==1) v=end?g->y2:g->y1;
  else v=end?g->z2:g->z1;
  return v[segment];
}

static int same_endpoint(const nec4_medium_geometry *g,size_t a,int ae,
                         size_t b,int be,double tolerance)
{
  double dx=endpoint(g,a,ae,0)-endpoint(g,b,be,0);
  double dy=endpoint(g,a,ae,1)-endpoint(g,b,be,1);
  double dz=endpoint(g,a,ae,2)-endpoint(g,b,be,2);
  return hypot(hypot(dx,dy),dz)<=tolerance;
}

void nec4_medium_state_init(nec4_medium_state *state)
{ memset(state,0,sizeof(*state)); }

void nec4_medium_state_free(nec4_medium_state *state)
{
  free(state->segment_medium);
  free(state->endpoint_factor);
  free(state->endpoint_factor_assigned);
  memset(state,0,sizeof(*state));
}

double complex nec4_medium_epsilon(double er,double loss,double frequency_mhz,
                                   int loss_is_epsilon)
{
  const double pi=acos(-1.);
  if(loss_is_epsilon) return er-I*fabs(loss);
  return er-I*loss/(2.*pi*frequency_mhz*1.e6*NEC4_EPS0);
}

int nec4_medium_classify(const nec4_medium_geometry *g,double er,double loss,
    double frequency_mhz,int loss_is_epsilon,nec4_medium_state *state)
{
  size_t i,j;double scale=1.,tolerance;
  const double pi=acos(-1.);
  nec4_medium_state_free(state);
  if(!g || !g->count || !g->x1 || !g->y1 || !g->z1 || !g->x2 || !g->y2 ||
     !g->z2 || !(frequency_mhz>0.) || !(er>0.) || (!loss_is_epsilon&&loss<0.))
    return 0;
  state->segment_medium=calloc(g->count,sizeof(*state->segment_medium));
  state->endpoint_factor=calloc(2*g->count,sizeof(*state->endpoint_factor));
  state->endpoint_factor_assigned=calloc(2*g->count,
      sizeof(*state->endpoint_factor_assigned));
  if(!state->segment_medium||!state->endpoint_factor||
     !state->endpoint_factor_assigned) { nec4_medium_state_free(state);return 0; }
  for(i=0;i<g->count;++i)
  {
    scale=fmax(scale,fabs(g->x1[i]));scale=fmax(scale,fabs(g->y1[i]));
    scale=fmax(scale,fabs(g->z1[i]));scale=fmax(scale,fabs(g->x2[i]));
    scale=fmax(scale,fabs(g->y2[i]));scale=fmax(scale,fabs(g->z2[i]));
  }
  tolerance=256.*DBL_EPSILON*scale;
  state->coordinate_scale_m=scale;state->interface_tolerance_m=tolerance;
  state->segment_count=g->count;
  state->frequency_mhz=frequency_mhz;state->relative_permittivity=er;
  state->conductivity_s_per_m=loss_is_epsilon?0.:loss;
  state->loss_is_epsilon=loss_is_epsilon;
  state->epsilon_relative=nec4_medium_epsilon(er,loss,frequency_mhz,loss_is_epsilon);
  state->refractive_index=csqrt(state->epsilon_relative);
  state->wave_number_per_m=2.*pi*frequency_mhz*1.e6/NEC4_C0*
    state->refractive_index;
  state->intrinsic_impedance_ohm=NEC4_ETA0/state->refractive_index;
  for(i=0;i<g->count;++i)
  {
    int a=side(g->z1[i],tolerance),b=side(g->z2[i],tolerance);
    if(a*b<0) { ++state->unsplit_crossings;continue; }
    if(!a&&!b) { ++state->in_plane_segments;continue; }
    state->segment_medium[i]=(a<0||b<0)?NEC4_MEDIUM_LOWER:NEC4_MEDIUM_UPPER;
    if(state->segment_medium[i]==NEC4_MEDIUM_LOWER) ++state->lower_segments;
    else ++state->upper_segments;
    if(!a) ++state->interface_endpoints;
    if(!b) ++state->interface_endpoints;
  }
  /* Count unique interface points joining at least one segment in each medium. */
  for(i=0;i<g->count;++i) for(int ie=0;ie<2;++ie)
  {
    int seen=0,upper_count=0,lower_count=0;
    if(side(endpoint(g,i,ie,2),tolerance)) continue;
    for(j=0;j<i&&!seen;++j) for(int je=0;je<2;++je)
      if(!side(endpoint(g,j,je,2),tolerance)&&
         same_endpoint(g,i,ie,j,je,tolerance)) { seen=1;break; }
    if(seen) continue;
    for(j=0;j<g->count;++j) for(int je=0;je<2;++je)
      if(!side(endpoint(g,j,je,2),tolerance)&&
         same_endpoint(g,i,ie,j,je,tolerance))
      {
        upper_count+=state->segment_medium[j]==NEC4_MEDIUM_UPPER;
        lower_count+=state->segment_medium[j]==NEC4_MEDIUM_LOWER;
      }
    if(upper_count==1&&lower_count==1)
    {
      double worst=1.;
      ++state->crossing_junctions;
      for(j=0;j<g->count;++j) for(int je=0;je<2;++je)
        if(!side(endpoint(g,j,je,2),tolerance)&&
           same_endpoint(g,i,ie,j,je,tolerance))
        {
          double dx=g->x2[j]-g->x1[j],dy=g->y2[j]-g->y1[j],
            dz=g->z2[j]-g->z1[j],length=hypot(hypot(dx,dy),dz);
          if(length>0.) worst=fmin(worst,fabs(dz)/length);
        }
      /* Audit the highly oblique region explicitly mentioned in the public
         guide.  This is not a rejection criterion: the guide reports that
         the normal-crossing condition was used beyond 60 degrees from the
         normal with apparently good results. */
      if(worst<.5) ++state->oblique_crossings;
      for(j=0;j<g->count;++j) for(int je=0;je<2;++je)
        if(!side(endpoint(g,j,je,2),tolerance)&&
           same_endpoint(g,i,ie,j,je,tolerance))
        {
          size_t ep=2*j+(size_t)je;
          state->endpoint_factor[ep]=state->segment_medium[j]==NEC4_MEDIUM_LOWER?
            state->epsilon_relative:1.;
          state->endpoint_factor_assigned[ep]=1;
        }
    }
    else ++state->unsupported_interface_junctions;
  }
  state->valid=state->unsplit_crossings==0&&state->in_plane_segments==0&&
    state->unsupported_interface_junctions==0;
  return state->valid;
}

int nec4_medium_endpoint_factor(const nec4_medium_state *state,size_t segment,
    int end,double complex *factor)
{
  size_t endpoint_index=2*segment+(end?1:0);
  if(!state||!factor||segment>=state->segment_count||
     end<0||end>1||!state->endpoint_factor_assigned[endpoint_index]) return 0;
  *factor=state->endpoint_factor[endpoint_index];
  return 1;
}

int nec4_medium_update_frequency(nec4_medium_state *state,double frequency_mhz,
    double complex epsilon_relative)
{
  const double pi=acos(-1.);
  if(!state||!state->valid||!(frequency_mhz>0.)||
     !(creal(epsilon_relative)>0.)||cimag(epsilon_relative)>0.||
     !isfinite(frequency_mhz)||!isfinite(creal(epsilon_relative))||
     !isfinite(cimag(epsilon_relative))) return 0;
  state->frequency_mhz=frequency_mhz;
  state->relative_permittivity=creal(epsilon_relative);
  state->epsilon_relative=epsilon_relative;
  state->conductivity_s_per_m=-cimag(epsilon_relative)*2.*pi*
    frequency_mhz*1.e6*NEC4_EPS0;
  state->refractive_index=csqrt(epsilon_relative);
  if(creal(state->refractive_index)<0.||
     (creal(state->refractive_index)==0.&&cimag(state->refractive_index)>0.))
    state->refractive_index=-state->refractive_index;
  state->wave_number_per_m=2.*pi*frequency_mhz*1.e6/NEC4_C0*
    state->refractive_index;
  state->intrinsic_impedance_ohm=NEC4_ETA0/state->refractive_index;
  for(size_t segment=0;segment<state->segment_count;++segment)
    for(int end=0;end<2;++end)
    {
      size_t index=2*segment+(size_t)end;
      if(state->endpoint_factor_assigned[index])
        state->endpoint_factor[index]=
          state->segment_medium[segment]==NEC4_MEDIUM_LOWER?
          epsilon_relative:1.;
    }
  return isfinite(creal(state->wave_number_per_m))&&
    isfinite(cimag(state->wave_number_per_m))&&
    isfinite(creal(state->intrinsic_impedance_ohm))&&
    isfinite(cimag(state->intrinsic_impedance_ohm));
}

void nec4_medium_print(FILE *stream,const nec4_medium_state *s)
{
  if(!stream||!s) return;
  fprintf(stream,"\n  NEC4 INTERFACE PREFLIGHT: %lu UPPER, %lu LOWER, "
      "%lu INTERFACE ENDPOINTS, %lu CROSSING JUNCTIONS"
      "\n                            %lu UNSUPPORTED INTERFACE JUNCTIONS, "
      "%lu UNSPLIT CROSSINGS, %lu IN-PLANE "
      "SEGMENTS, %lu OBLIQUE CROSSINGS; TOLERANCE %.6E M"
      "\n  LOWER MEDIUM: EPSILON_R %.9G%+.9GJ, N %.9G%+.9GJ, "
      "K %.9G%+.9GJ 1/M, ETA %.9G%+.9GJ OHM\n",
      (unsigned long)s->upper_segments,(unsigned long)s->lower_segments,
      (unsigned long)s->interface_endpoints,(unsigned long)s->crossing_junctions,
      (unsigned long)s->unsupported_interface_junctions,
      (unsigned long)s->unsplit_crossings,(unsigned long)s->in_plane_segments,
      (unsigned long)s->oblique_crossings,s->interface_tolerance_m,
      creal(s->epsilon_relative),cimag(s->epsilon_relative),
      creal(s->refractive_index),cimag(s->refractive_index),
      creal(s->wave_number_per_m),cimag(s->wave_number_per_m),
      creal(s->intrinsic_impedance_ohm),cimag(s->intrinsic_impedance_ohm));
}
