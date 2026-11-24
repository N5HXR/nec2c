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
#define EPS0 8.8541878128e-12

typedef struct { const char *name;double frequency,epsilon_r,conductivity; } material;
typedef struct { const char *name;unsigned int maximum_depth; } region;

static const material materials[]={
  {"average-soil",7.1e6,13.,.005},
  {"low-loss-water",5.e6,81.,.01},
  {"seawater",5.e6,81.,4.},
  {"lossless-water",299792458.,81.,0.}
};
static const region regions[]={
  {"switch-volume",4},
  {"interface-line",10},
  {"grazing-sheet",6},
  {"axis-volume",5}
};

static double seconds(void)
{
  struct timespec value;
  clock_gettime(CLOCK_MONOTONIC,&value);
  return value.tv_sec+1.e-9*value.tv_nsec;
}

static double halton(unsigned int index,unsigned int base)
{
  double fraction=1.,value=0.;
  while(index) { fraction/=base;value+=fraction*(index%base);index/=base; }
  return value;
}

static int direct_fine(double frequency,double complex lower_k,
    double complex upper_k,double rho,double z,double depth,double complex value[5])
{
  nec4_ground_direct_result direct;nec4_ground_interface_result interface;
  if(z+depth==0.) {
    if(!nec4_ground_direct_cross_interface(frequency,lower_k,upper_k,rho,
                                            5,80,120.,&interface)) return 0;
    memcpy(value,interface.field.total,sizeof(interface.field.total));return 1;
  }
  if(rho>0.&&z+depth<=.5*rho) {
    if(!nec4_ground_direct_cross_hankel(frequency,lower_k,upper_k,rho,z,depth,
                                        5,80,120.,&direct)) return 0;
  } else if(!nec4_ground_direct_cross_bessel(frequency,lower_k,upper_k,rho,z,
                                              depth,80,120.,&direct)) return 0;
  memcpy(value,direct.value,sizeof(direct.value));return 1;
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

static nec4_ground_adaptive_bounds bounds_for(const char *name,double wavelength,
    double lg)
{
  nec4_ground_adaptive_bounds bounds={{0}};
  if(!strcmp(name,"switch-volume")) {
    bounds.lower[0]=2.8*wavelength;bounds.upper[0]=3.2*wavelength;
    bounds.lower[1]=.015*wavelength;bounds.upper[1]=.035*wavelength;
    bounds.lower[2]=0.;bounds.upper[2]=.05*lg;
  } else if(!strcmp(name,"interface-line")) {
    bounds.lower[0]=.5*wavelength;bounds.upper[0]=4.*wavelength;
  } else if(!strcmp(name,"grazing-sheet")) {
    bounds.lower[0]=3.*wavelength;bounds.upper[0]=5.*wavelength;
    bounds.lower[1]=.01*wavelength;bounds.upper[1]=.4*wavelength;
    bounds.lower[2]=bounds.upper[2]=.02*lg;
  } else {
    bounds.lower[0]=0.;bounds.upper[0]=.5*wavelength;
    bounds.lower[1]=2.8*wavelength;bounds.upper[1]=3.2*wavelength;
    bounds.lower[2]=0.;bounds.upper[2]=.2*lg;
  }
  return bounds;
}

static void usage(const char *program)
{
  fprintf(stderr,"usage: %s [--material NAME] [--region NAME] [--workers N] "
                 "[--output PATH]\n",program);
}

int main(int argc,char **argv)
{
  const char *material_filter=NULL,*region_filter=NULL,*output_path=NULL;
  int workers=16;FILE *output=stdout;
  for(int argument=1;argument<argc;++argument) {
    if(!strcmp(argv[argument],"--material")&&argument+1<argc)
      material_filter=argv[++argument];
    else if(!strcmp(argv[argument],"--region")&&argument+1<argc)
      region_filter=argv[++argument];
    else if(!strcmp(argv[argument],"--workers")&&argument+1<argc)
      workers=(int)strtol(argv[++argument],NULL,10);
    else if(!strcmp(argv[argument],"--output")&&argument+1<argc)
      output_path=argv[++argument];
    else { usage(argv[0]);return 2; }
  }
  if(workers<1) { usage(argv[0]);return 2; }
  if(output_path) { output=fopen(output_path,"w");if(!output){perror(output_path);return 1;} }
  fprintf(output,"material,region,cells,leaves,accepted_leaves,unresolved_leaves,"
    "deepest_level,fitted_samples,withheld_samples,worst_leaf_error,"
    "worst_encountered_error,independent_worst_error,build_seconds,"
    "evaluation_seconds,direct_failures,adaptive_fallbacks\n");
  for(size_t m=0;m<sizeof(materials)/sizeof(materials[0]);++m) {
    const material *medium=&materials[m];
    double omega,wavelength,lg;double complex epsilon,upper_k,lower_k;
    if(material_filter&&strcmp(material_filter,medium->name)) continue;
    omega=2.*M_PI*medium->frequency;wavelength=C0/medium->frequency;
    epsilon=medium->epsilon_r-I*medium->conductivity/(omega*EPS0);
    upper_k=2.*M_PI/wavelength;lower_k=upper_k*csqrt(epsilon);
    lg=fmin(cabs(upper_k/lower_k),.3)*wavelength;
    for(size_t r=0;r<sizeof(regions)/sizeof(regions[0]);++r) {
      const region *area=&regions[r];
      nec4_ground_adaptive_bounds bounds;
      nec4_ground_adaptive_e table;
      const nec4_ground_table_direct_options fine={5,80,120.};
      double before,build_time,evaluation_time,worst=0.;
      int direct_failures=0,adaptive_fallbacks=0;
      if(region_filter&&strcmp(region_filter,area->name)) continue;
      bounds=bounds_for(area->name,wavelength,lg);
      nec4_ground_adaptive_e_init(&table);before=seconds();
      if(!nec4_ground_adaptive_e_build(&table,medium->frequency,lower_k,upper_k,
          &bounds,5.e-5,area->maximum_depth,
          !strcmp(area->name,"interface-line")?2048:512,workers,&fine)) {
        fprintf(stderr,"adaptive build failed: %s/%s\n",medium->name,area->name);
        return 1;
      }
      build_time=seconds()-before;before=seconds();
      for(unsigned int index=1;index<=128;++index) {
        const unsigned int bases[3]={2,3,5};double point[3];
        double complex expected[5],actual[5];
        for(unsigned int axis=0;axis<3;++axis)
          point[axis]=bounds.lower[axis]+(bounds.upper[axis]-bounds.lower[axis])*
                      halton(index+137,bases[axis]);
        if(!direct_fine(medium->frequency,lower_k,upper_k,point[0],point[1],
                        point[2],expected)) direct_failures++;
        else if(!nec4_ground_adaptive_e_evaluate(&table,point[0],point[1],
                                                  point[2],actual))
          adaptive_fallbacks++;
        else worst=fmax(worst,scaled_error(actual,expected));
      }
      evaluation_time=seconds()-before;
      fprintf(output,"%s,%s,%zu,%zu,%zu,%zu,%u,%zu,%zu,%.17g,%.17g,"
        "%.17g,%.9g,%.9g,%d,%d\n",medium->name,area->name,
        table.table.stats.cells,table.table.stats.leaves,
        table.table.stats.accepted_leaves,table.table.stats.unresolved_leaves,
        table.table.stats.deepest_level,table.table.stats.fitted_samples,
        table.table.stats.withheld_samples,
        table.table.stats.worst_withheld_scaled_error,
        table.table.stats.worst_encountered_scaled_error,worst,build_time,
        evaluation_time,direct_failures,adaptive_fallbacks);
      fflush(output);nec4_ground_adaptive_e_destroy(&table);
    }
  }
  if(output_path&&fclose(output)) { perror(output_path);return 1; }
  return 0;
}
