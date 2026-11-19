#ifndef NEC2C_NEC4_GROUND_MATRIX_H
#define NEC2C_NEC4_GROUND_MATRIX_H 1

#include "nec4_medium.h"
#include "nec4_ground_table.h"

#include <complex.h>

typedef int (*nec4_ground_cross_evaluator)(void *context,double rho_m,
    double observer_z_m,double source_depth_m,
    double complex field[5]);

typedef struct {
  double frequency_hz;
  double complex lower_epsilon_relative;
  double complex lower_wave_number_per_m;
  double complex upper_wave_number_per_m;
  nec4_ground_cross_evaluator cross;
  void *cross_context;
} nec4_ground_matrix_context;

typedef struct {
  double complex constant[3];
  double complex sine[3];
  double complex cosine[3];
} nec4_ground_segment_field;

typedef struct {
  const nec4_ground_table *table[3];
  double complex lower_wave_number_per_m;
  double complex upper_wave_number_per_m;
} nec4_ground_matrix_tables;

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
   upper-medium homogeneous dyadic over one straight source segment.  This is
   the correction added to NEC's finite-radius direct wire field. */
int nec4_ground_matrix_segment_correction(
    const nec4_ground_matrix_context *context,
    nec4_segment_medium source_medium,nec4_segment_medium observer_medium,
    const double source_center_m[3],const double source_direction[3],
    double source_length_m,const double observer_m[3],
    double regularization_radius_m,nec4_ground_segment_field *field);

#endif
