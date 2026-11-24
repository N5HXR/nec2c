#ifndef NEC2C_NEC4_GROUND_DISPATCH_H
#define NEC2C_NEC4_GROUND_DISPATCH_H 1

#include <complex.h>

typedef enum {
  NEC4_GROUND_LINEAR_TABLE=1,
  NEC4_GROUND_MODEL_FIT=2,
  NEC4_GROUND_ASYMPTOTIC=3
} nec4_ground_method;

typedef struct {
  nec4_ground_method method;
  int model_region;
  int used_direct_fallback;
  int grazing_risk;
  double grazing_angle_degrees;
  double separation_wavelengths;
} nec4_ground_dispatch;

int nec4_ground_route(double rho,double observer_z,double source_z,
    double complex k_lower,double complex k_upper,double wavelength0,
    nec4_ground_dispatch *route);
const char *nec4_ground_method_name(nec4_ground_method method);

#endif
