#ifndef NEC2C_NEC4_GROUND_ASYMPTOTIC_H
#define NEC2C_NEC4_GROUND_ASYMPTOTIC_H 1

#include "nec4_ground_direct.h"

typedef struct {
  double complex value[NEC4_GROUND_COMPONENT_COUNT];
  double complex saddle;
} nec4_ground_asymptotic_result;

/* Section 5.1.3 equations (5-21)--(5-23), excluding the pole term. */
int nec4_ground_asymptotic_first(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,
    nec4_ground_asymptotic_result *result);

/* Appendix-J second-order saddle expansion with the printed two-term Hankel
   amplitudes.  This deliberately excludes equations (5-26)/(5-28), allowing
   the surface-wave correction to be qualified independently. */
int nec4_ground_asymptotic_second_saddle(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,
    nec4_ground_asymptotic_result *result);

/* Equations (5-25)--(5-28), kept separate from the saddle contribution so
   callers and tests can observe the surface-wave correction explicitly. */
int nec4_ground_asymptotic_first_pole(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,
    nec4_ground_asymptotic_result *result);
int nec4_ground_asymptotic_second_pole(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,
    nec4_ground_asymptotic_result *result);

/* Complete first- and second-order uniform outer fields.  At rho=0 the pole
   representation is inapplicable and the analytic Appendix-J axis value is
   returned by the second-order entry point. */
int nec4_ground_asymptotic_first_uniform(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,
    nec4_ground_asymptotic_result *result);
int nec4_ground_asymptotic_second_uniform(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,
    nec4_ground_asymptotic_result *result);

#endif
