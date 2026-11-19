#ifndef NEC2C_WIRE_OPERATOR_H
#define NEC2C_WIRE_OPERATOR_H 1

#include <complex.h>
#include <stdio.h>
#include "nec4_medium.h"
#include "nec4_ground_matrix.h"

/*
 * Public engine-selection boundary for wire and ground physics.
 *
 * M1 deliberately exposes only the unchanged NEC2 implementation.  Later
 * milestones add implementations behind this interface without changing the
 * classic entry points or their data representation.
 */

int wire_operator_set_engine(const char *name);
int wire_operator_set_ground_engine(const char *name);
int wire_operator_set_model_check(const char *name);
void wire_operator_set_model_features(size_t patches, size_t coated_domains,
                                      int dielectric_loaded);
void wire_operator_set_medium_state(const nec4_medium_state *state);
void wire_operator_set_ground_matrix_context(
    const nec4_ground_matrix_context *context);
void wire_operator_set_source_segment(size_t segment);
int wire_operator_ground_direct_only(void);

const char *wire_operator_requested_engine(void);
const char *wire_operator_resolved_engine(void);
const char *wire_operator_requested_ground_engine(void);
const char *wire_operator_resolved_ground_engine(void);
const char *wire_operator_model_check(void);
int wire_operator_nec4_basis(void);
int wire_operator_nec4_ground(void);

typedef struct
{
  size_t short_thick_segments;
  size_t close_parallel_pairs;
  size_t buried_match_points;
  size_t overlapping_bends;
  size_t source_length_mismatches;
  size_t radius_step_junctions;
  size_t legacy_discontinuity_sources;
  double minimum_length_radius_ratio;
} wire_operator_validity;

int wire_operator_prepare_frequency(FILE *report);
const wire_operator_validity *wire_operator_last_validity(void);
double complex wire_operator_charge_factor(size_t segment, int end,
                                           double complex fallback);
double complex wire_operator_segment_wave_number(size_t segment);
double complex wire_operator_source_wave_number(void);
double complex wire_operator_source_intrinsic_impedance(void);
double complex wire_operator_segment_phase(size_t segment,
                                            double offset_wavelengths);
double complex wire_operator_end_parameter(size_t segment,double radius);

/* Representation-aware operations for A + B sin(x) + C f_cos(x). */
double complex wire_operator_cosine_shape(double complex angle);
double complex wire_operator_center_basis(double complex constant,
                                          double complex cosine);
double wire_operator_cosine_integral(double constant_integral,
                                     double cosine_integral);
double complex wire_operator_current(double complex constant,
                                     double complex sine,
                                     double complex cosine,
                                     size_t segment,double offset_wavelengths);

void wire_operator_sbf(int i, int is, double complex *aa,
                       double complex *bb, double complex *cc);
void wire_operator_tbf(int i, int icap);
void wire_operator_trio(int j);
void wire_operator_efld(double xi, double yi, double zi, double ai, int ij);
void wire_operator_efld_point(double xi,double yi,double zi,double ai);
void wire_operator_hsfld(double xi, double yi, double zi, double ai);

/* NEC4 flat-cap fields in wavelength-normalized coordinates. */
int wire_operator_endcap_current_field(size_t segment, int end,
    double xi, double yi, double zi, double complex *ex,
    double complex *ey, double complex *ez);
int wire_operator_gap_voltage_field(size_t segment,
    double xi, double yi, double zi, double complex *ex,
    double complex *ey, double complex *ez);

#endif
