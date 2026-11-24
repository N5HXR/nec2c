#ifndef NEC2C_NEC4_GROUND_MATRIX_H
#define NEC2C_NEC4_GROUND_MATRIX_H 1

#include "nec4_medium.h"
#include "nec4_ground_table.h"

#include <complex.h>

typedef int (*nec4_ground_cross_evaluator)(void *context,double rho_m,
    double observer_z_m,double source_depth_m,
    double complex field[5]);

typedef int (*nec4_ground_cross_h_evaluator)(void *context,double rho_m,
    double observer_z_m,double source_depth_m,
    double complex field[NEC4_GROUND_H_COMPONENT_COUNT]);

typedef struct {
  double frequency_hz;
  double complex lower_epsilon_relative;
  double complex lower_wave_number_per_m;
  double complex upper_wave_number_per_m;
  nec4_ground_cross_evaluator cross;
  nec4_ground_cross_h_evaluator cross_h;
  void *cross_context;
} nec4_ground_matrix_context;

typedef struct {
  double complex constant[3];
  double complex sine[3];
  double complex cosine[3];
} nec4_ground_segment_field;

/* Complete lower-medium source to upper-medium observer H field.  M11D
   independently qualified this canonical transmission pairing; same-medium
   H stays on the legacy path until its transformations are derived. */
int nec4_ground_matrix_h_point_cross(
    const nec4_ground_matrix_context *context,const double source_m[3],
    const double observer_m[3],const double source_direction[3],
    double complex field[3]);
int nec4_ground_matrix_h_segment_nec_cross(
    const nec4_ground_matrix_context *context,const double source_center_m[3],
    const double source_direction[3],double source_length_m,
    const double observer_m[3],double normalization_length_m,
    nec4_ground_segment_field *field);

typedef struct {
  const nec4_ground_table *table[3];
  double complex lower_wave_number_per_m;
  double complex upper_wave_number_per_m;
} nec4_ground_matrix_tables;

int nec4_ground_matrix_homogeneous_components(double frequency_hz,
    double complex wave_number_per_m,double rho_m,double source_z_m,
    double observer_z_m,double complex field[5]);

/* Canonical callback for a complete three-table M9 bank.  Shared depth
   boundaries belong to the lower-numbered table. */
int nec4_ground_matrix_table_cross(void *context,double rho_m,
    double observer_z_m,double source_depth_m,double complex field[5]);

/* Appendix-G completion of the canonical lower-source/upper-observer field.
   source_direction is a unit current direction.  The returned field is due
   to a unit electric-current element (A m). */
int nec4_ground_matrix_point(const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    const double source_m[3],const double observer_m[3],
    const double source_direction[3],double complex field[3]);

/* Integrate the difference between the complete interface dyadic and the
   source-medium homogeneous dyadic over one straight source segment.  This is
   the correction added to NEC's finite-radius direct wire field. */
int nec4_ground_matrix_segment_correction(
    const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    const double source_center_m[3],const double source_direction[3],
    double source_length_m,const double observer_m[3],
    double regularization_radius_m,nec4_ground_segment_field *field);

/* Convert the physical E/I segment correction (ohm/m) to NEC's matrix
   normalization (ohm).  NEC solves for current per metre and multiplies the
   solved coefficient by wavelength before reporting amperes. */
int nec4_ground_matrix_segment_nec_correction(
    const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    const double source_center_m[3],const double source_direction[3],
    double source_length_m,const double observer_m[3],
    double regularization_radius_m,double normalization_length_m,
    nec4_ground_segment_field *field);

/* Complete interface field for a source segment, in NEC matrix units.  This
   path is required for lower-medium sources whose complex current basis and
   homogeneous direct kernel are not the free-space NEC kernel. */
int nec4_ground_matrix_segment_nec_complete(
    const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    const double source_center_m[3],const double source_direction[3],
    double source_length_m,const double observer_m[3],
    double regularization_radius_m,double normalization_length_m,
    nec4_ground_segment_field *field);

/* Homogeneous-medium point-charge field omitted by the reduced segment
   kernel.  The sign is +1 at source end 1 and -1 at source end 2 for the
   exp(+j omega t) convention.  Returned units are NEC matrix ohms for unit
   endpoint current shape and current-per-metre matrix unknown. */
int nec4_ground_matrix_endpoint_charge_nec(
    const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,const double endpoint_m[3],
    const double observer_m[3],double normalization_length_m,int end_sign,
    double complex field[3]);

#endif
