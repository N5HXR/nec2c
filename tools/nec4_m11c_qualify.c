#define _POSIX_C_SOURCE 200809L

#include "nec4_ground_adaptive_e.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define C0 299792458.0

static double seconds(void)
{
  struct timespec value;clock_gettime(CLOCK_MONOTONIC,&value);
  return value.tv_sec+1.e-9*value.tv_nsec;
}

static double halton(unsigned int index,unsigned int base)
{
  double fraction=1.,value=0.;
  while(index) { fraction/=base;value+=fraction*(index%base);index/=base; }
  return value;
}

static double scaled_error(const double complex actual[5],
    const double complex expected[5])
{
  double difference=0.,scale=0.;
  for(unsigned int component=0;component<5;++component) {
    difference=fmax(difference,cabs(actual[component]-expected[component]));
    scale=fmax(scale,cabs(expected[component]));
  }
  return difference/fmax(scale,1.e-300);
}

static nec4_ground_adaptive_bounds bounds_for(const char *region,double wl,
    double lg)
{
  nec4_ground_adaptive_bounds bounds={{0},{0}};
  if(!strcmp(region,"switch-volume")) {
    bounds.lower[0]=2.8*wl;bounds.upper[0]=3.2*wl;
    bounds.lower[1]=.015*wl;bounds.upper[1]=.035*wl;
    bounds.lower[2]=.002*lg;bounds.upper[2]=.05*lg;
  } else {
    bounds.lower[0]=3.*wl;bounds.upper[0]=5.*wl;
    bounds.lower[1]=.01*wl;bounds.upper[1]=.4*wl;
    bounds.lower[2]=bounds.upper[2]=.02*lg;
  }
  return bounds;
}

int main(int argc,char **argv)
{
  const char *region=argc>1?argv[1]:"grazing-sheet";
  int workers=argc>2?atoi(argv[2]):16;
  size_t directional_budget=argc>3?(size_t)strtoull(argv[3],NULL,10):4096;
  const double frequency=C0,epsilon_r=81.,wl=C0/frequency;
  const double complex upper_k=2.*M_PI/wl,lower_k=upper_k*sqrt(epsilon_r);
  const double lg=fmin(cabs(upper_k/lower_k),.3)*wl;
  const nec4_ground_adaptive_bounds bounds=bounds_for(region,wl,lg);
  const nec4_ground_table_direct_options options={5,80,120.};
  nec4_ground_adaptive_e baseline,directional;
  double start,baseline_seconds,directional_seconds;
  double baseline_worst=0.,directional_worst=0.,maximum_cancellation=0.;
  unsigned int baseline_fallback=0,directional_fallback=0,both_fallback=0;
  if((strcmp(region,"switch-volume")&&strcmp(region,"grazing-sheet"))||
     workers<1||directional_budget<1) return 2;
  nec4_ground_adaptive_e_init(&baseline);
  nec4_ground_adaptive_e_init(&directional);
  start=seconds();
  if(!nec4_ground_adaptive_e_build(&baseline,frequency,lower_k,upper_k,&bounds,
      5.e-5,!strcmp(region,"switch-volume")?4:6,512,workers,&options))return 1;
  baseline_seconds=seconds()-start;start=seconds();
  if(!nec4_ground_adaptive_e_build_anisotropic(&directional,frequency,
      lower_k,upper_k,&bounds,5.e-5,30,directional_budget,workers,&options))
    return 1;
  directional_seconds=seconds()-start;
  for(unsigned int index=1;index<=256;++index) {
    double point[3];double complex expected[5],actual[5];
    nec4_ground_direct_split_result direct;
    int baseline_ok,directional_ok;
    point[0]=bounds.lower[0]+halton(index,2)*(bounds.upper[0]-bounds.lower[0]);
    point[1]=bounds.lower[1]+halton(index,3)*(bounds.upper[1]-bounds.lower[1]);
    point[2]=bounds.lower[2]+halton(index,5)*(bounds.upper[2]-bounds.lower[2]);
    if(!nec4_ground_direct_cross_hankel_split(frequency,lower_k,upper_k,
        point[0],point[1],point[2],5,80,120.,&direct)) return 1;
    memcpy(expected,direct.value,sizeof(expected));
    for(unsigned int component=0;component<5;++component)
      maximum_cancellation=fmax(maximum_cancellation,
        fmax(cabs(direct.lower_spectral[component]),
             cabs(direct.upper_spectral[component]))/
        fmax(cabs(direct.value[component]),1.e-300));
    baseline_ok=nec4_ground_adaptive_e_evaluate(&baseline,point[0],point[1],
      point[2],actual);
    if(baseline_ok) baseline_worst=fmax(baseline_worst,
      scaled_error(actual,expected));
    else baseline_fallback++;
    directional_ok=nec4_ground_adaptive_e_evaluate(&directional,point[0],
      point[1],point[2],actual);
    if(directional_ok) directional_worst=fmax(directional_worst,
      scaled_error(actual,expected));
    else directional_fallback++;
    if(!baseline_ok&&!directional_ok)both_fallback++;
  }
  printf("region,baseline_cells,baseline_unresolved,baseline_fallback,"
    "baseline_worst,directional_budget,directional_cells,"
    "directional_unresolved,directional_fallback,directional_worst,"
    "both_fallback,maximum_cancellation,baseline_seconds,"
    "directional_seconds\n");
  printf("%s,%zu,%zu,%u,%.17g,%zu,%zu,%zu,%u,%.17g,%u,%.17g,%.9g,%.9g\n",
    region,baseline.table.stats.cells,baseline.table.stats.unresolved_leaves,
    baseline_fallback,baseline_worst,directional_budget,
    directional.table.stats.cells,directional.table.stats.unresolved_leaves,
    directional_fallback,directional_worst,both_fallback,
    maximum_cancellation,baseline_seconds,directional_seconds);
  nec4_ground_adaptive_e_destroy(&directional);
  nec4_ground_adaptive_e_destroy(&baseline);return 0;
}
