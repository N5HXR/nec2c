#include "nec4_ground_dispatch.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static nec4_ground_dispatch route(double rho,double z,double depth)
{
  const double wavelength=1.;
  /* |k_upper/k_lower|=0.25 gives Lg=0.25 and Lmax=1 wavelength. */
  double complex upper=2.*M_PI;
  double complex lower=4.*upper;
  nec4_ground_dispatch result;
  assert(nec4_ground_route(rho,z,-depth,lower,upper,wavelength,&result));
  return result;
}

int main(void)
{
  nec4_ground_dispatch r;
  r=route(.25,.25,.1);assert(r.method==NEC4_GROUND_LINEAR_TABLE);
  r=route(.2500001,.25,.1);assert(r.method==NEC4_GROUND_MODEL_FIT);
  r=route(.1,.099,.01);assert(r.method==NEC4_GROUND_LINEAR_TABLE);
  r=route(.5,.25-.000001,.01);assert(r.model_region==2);
  r=route(.5,.25,.01);assert(r.model_region==1);
  r=route(.600001,.6,.01);assert(r.model_region==3);
  r=route(.600001,.600001,.01);assert(r.model_region==1);
  r=route(3.,2.,.2);assert(r.method==NEC4_GROUND_MODEL_FIT);
  r=route(3.000001,.2,.01);assert(r.method==NEC4_GROUND_ASYMPTOTIC);
  r=route(4.,.02,.02);assert(r.grazing_risk);
  r=route(4.,.2,.02);assert(r.grazing_risk);
  r=route(4.,.8,.8);assert(!r.grazing_risk);
  assert(!strcmp(nec4_ground_method_name(NEC4_GROUND_ASYMPTOTIC),
                 "asymptotic"));
  assert(!nec4_ground_route(-1.,0.,0.,1.,1.,1.,&r));
  return 0;
}
