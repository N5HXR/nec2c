#ifndef NEC2C_JUNCTION_CHARGE_H
#define NEC2C_JUNCTION_CHARGE_H 1

#include <stddef.h>
#include <stdio.h>
#include <complex.h>

typedef struct
{
  size_t junctions;
  size_t solved;
  size_t symmetric;
  size_t prescribed;
  size_t fallback;
  size_t maximum_order;
  double maximum_scaled_residual;
  double maximum_pivot_ratio;
} junction_charge_diagnostics;

/*
 * Section 3.4 / Appendix A local frequency-domain junction solution. Directions
 * point away from the junction.  Factors are normalized to unit arithmetic
 * mean because only their ratios enter the wire basis.
 */
int junction_charge_solve(size_t count, const double *length,
                          const double *radius, const double (*direction)[3],
                          double complex *factor, double *scaled_residual,
                          double *pivot_ratio);

/* Build the endpoint graph and prepare factors for the current frequency. */
int junction_charge_prepare(FILE *report, int strict,
                            const unsigned char *prescribed_endpoints);
double complex junction_charge_factor(size_t segment, int end,
                                      double complex fallback);
const junction_charge_diagnostics *junction_charge_last_diagnostics(void);

#endif
