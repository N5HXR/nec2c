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
int coated_wire_add(coated_wire_model *model, int tag, int first, int last,
    double epsilon_r, double outer_radius_m, double loss_tangent,
    double conductivity_s_per_m, coated_wire_diagnostic *diagnostic);
int coated_wire_finalize(coated_wire_model *model, const data_t *geometry,
    coated_wire_diagnostic *diagnostic);
int coated_wire_domain_equal(const coated_wire_domain *left,
    const coated_wire_domain *right);
void coated_wire_print_summary(FILE *stream, const coated_wire_model *model);
void coated_wire_print_diagnostic(FILE *stream,
    const coated_wire_diagnostic *diagnostic);

#endif
