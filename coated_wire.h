#ifndef COATED_WIRE_H
#define COATED_WIRE_H 1

#include "nec2c.h"

#include <stddef.h>
#include <stdio.h>

typedef struct
{
  int tag;
  int first;
  int last;
  double epsilon_r;
  double outer_radius_m;
  double loss_tangent;
  double conductivity_s_per_m;
  double conductor_radius_m;
  size_t segment_count;
} coated_wire_domain;

typedef struct
{
  coated_wire_domain *domains;
  size_t domain_count;
  int *segment_domain;
  size_t segment_count;
  int finalized;
  int operator_enabled;
  int experimental_junctions_enabled;
  int experimental_junction_count;
  int experimental_geometry;
  int experimental_geometry_warning_emitted;
  int qualified_bend_count;
  double qualified_bend_maximum_turn_deg;
  double qualified_bend_minimum_chord_ratio;
} coated_wire_model;

typedef struct
{
  const char *code;
  char message[192];
  int tag;
  int segment;
  double value;
  int value_defined;
} coated_wire_diagnostic;

void coated_wire_model_init(coated_wire_model *model);
void coated_wire_model_clear(coated_wire_model *model);
void coated_wire_allow_experimental_junctions(coated_wire_model *model,
    int enabled);
int coated_wire_add(coated_wire_model *model, int tag, int first, int last,
    double epsilon_r, double outer_radius_m, double loss_tangent,
    double conductivity_s_per_m, coated_wire_diagnostic *diagnostic);
int coated_wire_finalize(coated_wire_model *model, const data_t *geometry,
    coated_wire_diagnostic *diagnostic);
int coated_wire_validate_geometry(coated_wire_model *model,
    const data_t *geometry, coated_wire_diagnostic *diagnostic);
int coated_wire_domain_equal(const coated_wire_domain *left,
    const coated_wire_domain *right);
void coated_wire_activate(const coated_wire_model *model, double wavelength_m);
complex double coated_wire_relative_permittivity(
    const coated_wire_domain *domain, double frequency_mhz);
int coated_wire_has_loss(void);
double coated_wire_absorbed_power(void);
double coated_wire_radiated_power(size_t sample_count);
int coated_wire_segment_active(int segment);
void coated_wire_field_correction(int source_segment, int observation_segment,
    complex double *sine_field, complex double *cosine_field);
void coated_wire_print_summary(FILE *stream, const coated_wire_model *model);
void coated_wire_print_diagnostic(FILE *stream,
    const coated_wire_diagnostic *diagnostic);

#endif
