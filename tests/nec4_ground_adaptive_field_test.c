#define _POSIX_C_SOURCE 200809L

#include "nec4_ground_adaptive_e.h"
#include "nec4_ground_adaptive_h.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main(void)
{
  const double c0=299792458.,frequency=c0,wavelength=c0/frequency;
  const double complex upper_k=2.*acos(-1.)/wavelength;
  const double complex lower_k=2.*upper_k;
  const nec4_ground_table_direct_options options={3,48,80.};
  const nec4_ground_adaptive_bounds bounds={
    {2.9*wavelength,.02*wavelength,.01*wavelength},
    {3.1*wavelength,.03*wavelength,.015*wavelength}};
  char epath[]="/tmp/nec2c-adaptive-e-XXXXXX";
  char hpath[]="/tmp/nec2c-adaptive-h-XXXXXX";
  int efd=mkstemp(epath),hfd=mkstemp(hpath);
  nec4_ground_adaptive_e e,eloaded;
  nec4_ground_adaptive_h h,hloaded;
  double complex ev[5],ev2[5],hv[4],hv2[4];
  assert(efd>=0&&hfd>=0);close(efd);close(hfd);unlink(epath);unlink(hpath);
  nec4_ground_adaptive_e_init(&e);nec4_ground_adaptive_e_init(&eloaded);
  nec4_ground_adaptive_h_init(&h);nec4_ground_adaptive_h_init(&hloaded);
  assert(nec4_ground_adaptive_e_build_anisotropic(&e,frequency,lower_k,
    upper_k,&bounds,1.e-3,8,128,4,&options));
  assert(nec4_ground_adaptive_h_build_anisotropic(&h,frequency,lower_k,
    upper_k,&bounds,2.e-3,8,128,4,&options));
  assert(nec4_ground_adaptive_e_store(&e,epath));
  assert(nec4_ground_adaptive_h_store(&h,hpath));
  assert(nec4_ground_adaptive_e_cache_key(frequency,lower_k,upper_k,&bounds,
    1.e-3,8,128,1,&options)!=
    nec4_ground_adaptive_h_cache_key(frequency,lower_k,upper_k,&bounds,
    2.e-3,8,128,&options));
  assert(nec4_ground_adaptive_e_load_anisotropic(&eloaded,epath,frequency,
    lower_k,upper_k,&bounds,1.e-3,8,128,&options));
  assert(nec4_ground_adaptive_h_load(&hloaded,hpath,frequency,lower_k,upper_k,
    &bounds,2.e-3,8,128,&options));
  assert(nec4_ground_adaptive_e_evaluate(&e,3.*wavelength,.025*wavelength,
    .0125*wavelength,ev));
  assert(nec4_ground_adaptive_e_evaluate(&eloaded,3.*wavelength,
    .025*wavelength,.0125*wavelength,ev2));
  assert(nec4_ground_adaptive_h_evaluate(&h,3.*wavelength,.025*wavelength,
    .0125*wavelength,hv));
  assert(nec4_ground_adaptive_h_evaluate(&hloaded,3.*wavelength,
    .025*wavelength,.0125*wavelength,hv2));
  assert(!memcmp(ev,ev2,sizeof(ev))&&!memcmp(hv,hv2,sizeof(hv)));
  /* E/H identity and tensor/directional identity are not interchangeable. */
  assert(!nec4_ground_adaptive_e_load(&eloaded,epath,frequency,lower_k,
    upper_k,&bounds,1.e-3,8,128,&options));
  assert(!nec4_ground_adaptive_h_load(&hloaded,epath,frequency,lower_k,
    upper_k,&bounds,2.e-3,8,128,&options));
  nec4_ground_adaptive_e_destroy(&e);nec4_ground_adaptive_e_destroy(&eloaded);
  nec4_ground_adaptive_h_destroy(&h);nec4_ground_adaptive_h_destroy(&hloaded);
  unlink(epath);unlink(hpath);
  puts("NEC4ish adaptive E/H cache identity and directional topology passed");
  return 0;
}
