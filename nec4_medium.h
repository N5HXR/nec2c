#ifndef NEC2C_NEC4_MEDIUM_H
#define NEC2C_NEC4_MEDIUM_H 1

#include <complex.h>
#include <stddef.h>
#include <stdio.h>

typedef enum {
  NEC4_MEDIUM_INVALID=0,
  NEC4_MEDIUM_UPPER=1,
  NEC4_MEDIUM_LOWER=2
} nec4_segment_medium;

typedef struct {
  size_t count;
  const double *x1,*y1,*z1,*x2,*y2,*z2;
} nec4_medium_geometry;

typedef struct {
  int valid;
  size_t segment_count;
  size_t upper_segments;
  size_t lower_segments;
  size_t interface_endpoints;
  size_t crossing_junctions;
  size_t unsupported_interface_junctions;
  size_t unsplit_crossings;
  size_t in_plane_segments;
  size_t oblique_crossings;
  double coordinate_scale_m;
  double interface_tolerance_m;
  double frequency_mhz;
  double relative_permittivity;
  double conductivity_s_per_m;
  int loss_is_epsilon;
  double complex epsilon_relative;
  double complex refractive_index;
  double complex wave_number_per_m;
  double complex intrinsic_impedance_ohm;
  nec4_segment_medium *segment_medium;
  double complex *endpoint_factor;
  unsigned char *endpoint_factor_assigned;
} nec4_medium_state;

void nec4_medium_state_init(nec4_medium_state *state);
void nec4_medium_state_free(nec4_medium_state *state);
double complex nec4_medium_epsilon(double relative_permittivity,
    double conductivity_or_epsilon_loss,double frequency_mhz,int loss_is_epsilon);
int nec4_medium_classify(const nec4_medium_geometry *geometry,
    double relative_permittivity,double conductivity_or_epsilon_loss,
    double frequency_mhz,int loss_is_epsilon,nec4_medium_state *state);
int nec4_medium_endpoint_factor(const nec4_medium_state *state,size_t segment,int end,
    double complex *factor);
void nec4_medium_print(FILE *stream,const nec4_medium_state *state);

#endif
