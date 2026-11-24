#define _POSIX_C_SOURCE 200809L

#include "nec4_ground_adaptive_h.h"

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

static double scaled_error(
    const double complex actual[NEC4_GROUND_H_COMPONENT_COUNT],
    const double complex expected[NEC4_GROUND_H_COMPONENT_COUNT])
{
  double difference=0.,scale=0.;
  for(unsigned int component=0;
      component<NEC4_GROUND_H_COMPONENT_COUNT;++component) {
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
  size_t budget=argc>3?(size_t)strtoull(argv[3],NULL,10):4096;
  const double frequency=C0,epsilon_r=81.,wl=C0/frequency;
  const double complex upper_k=2.*M_PI/wl,lower_k=upper_k*sqrt(epsilon_r);
  const double lg=fmin(cabs(upper_k/lower_k),.3)*wl;
  const nec4_ground_adaptive_bounds bounds=bounds_for(region,wl,lg);
  const nec4_ground_table_direct_options options={5,80,120.};
  nec4_ground_adaptive_h table;
  double start,build_seconds,worst=0.;unsigned int fallback=0;
  if((strcmp(region,"switch-volume")&&strcmp(region,"grazing-sheet"))||
     workers<1||budget<1) return 2;
  nec4_ground_adaptive_h_init(&table);start=seconds();
  if(!nec4_ground_adaptive_h_build_anisotropic(&table,frequency,lower_k,
      upper_k,&bounds,3.e-4,30,budget,workers,&options)) return 1;
  build_seconds=seconds()-start;
  for(unsigned int index=1;index<=128;++index) {
    double point[3];double complex expected[NEC4_GROUND_H_COMPONENT_COUNT];
    double complex actual[NEC4_GROUND_H_COMPONENT_COUNT];
    nec4_ground_direct_h_result direct;
    point[0]=bounds.lower[0]+halton(index,2)*(bounds.upper[0]-bounds.lower[0]);
    point[1]=bounds.lower[1]+halton(index,3)*(bounds.upper[1]-bounds.lower[1]);
    point[2]=bounds.lower[2]+halton(index,5)*(bounds.upper[2]-bounds.lower[2]);
    if(!nec4_ground_direct_cross_h_hankel(frequency,lower_k,upper_k,
        point[0],point[1],point[2],5,80,120.,&direct)) return 1;
    memcpy(expected,direct.value,sizeof(expected));
    if(nec4_ground_adaptive_h_evaluate(&table,point[0],point[1],point[2],
        actual)) worst=fmax(worst,scaled_error(actual,expected));
    else fallback++;
  }
  printf("region,cells,unresolved,fallback,worst,build_seconds\n");
  printf("%s,%zu,%zu,%u,%.17g,%.9g\n",region,table.table.stats.cells,
    table.table.stats.unresolved_leaves,fallback,worst,build_seconds);
  nec4_ground_adaptive_h_destroy(&table);return 0;
}
