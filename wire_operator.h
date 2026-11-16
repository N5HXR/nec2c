#ifndef NEC2C_WIRE_OPERATOR_H
#define NEC2C_WIRE_OPERATOR_H 1

#include <complex.h>

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

const char *wire_operator_requested_engine(void);
const char *wire_operator_resolved_engine(void);
const char *wire_operator_requested_ground_engine(void);
const char *wire_operator_resolved_ground_engine(void);
const char *wire_operator_model_check(void);
int wire_operator_nec4_basis(void);

/* Representation-aware operations for A + B sin(x) + C f_cos(x). */
double wire_operator_cosine_shape(double angle);
double wire_operator_center_basis(double constant, double cosine);
double wire_operator_cosine_integral(double constant_integral,
                                     double cosine_integral);
double complex wire_operator_current(double complex constant,
                                     double complex sine,
                                     double complex cosine,
                                     double angle);

void wire_operator_sbf(int i, int is, double *aa, double *bb, double *cc);
void wire_operator_tbf(int i, int icap);
void wire_operator_trio(int j);
void wire_operator_efld(double xi, double yi, double zi, double ai, int ij);
void wire_operator_hsfld(double xi, double yi, double zi, double ai);

#endif
