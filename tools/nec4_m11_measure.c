#define _POSIX_C_SOURCE 200809L

#include "nec4_ground_asymptotic.h"
#include "nec4_ground_direct.h"
#include "nec4_ground_dispatch.h"

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

typedef struct {
  const char *name;
  double frequency_hz,epsilon_r,conductivity;
} material;

typedef struct {
  const char *family;
  double rho_lambda,z_lambda,depth_lg;
} witness;

static const material materials[] = {
  {"average-soil",7.1e6,13.,.005},
  {"low-loss-water",5.e6,81.,.01},
  {"seawater",5.e6,81.,4.},
  {"lossless-water",299792458.,81.,0.}
};

/* These points intentionally concentrate on the original-NEC4 weak region.
   Observer height is in free-space wavelengths; source depth is in Lg. */
static const witness witnesses[] = {
  {"switch",2.90,.025,.025},
  {"switch",3.00,.025,.025},
  {"switch",3.01,.025,.025},
  {"switch",3.10,.025,.025},
  {"grazing-angle",4.00,.069816,.02}, /* approximately 1 degree */
  {"grazing-angle",4.00,.349955,.02}, /* approximately 5 degrees */
  {"grazing-angle",4.00,.705308,.02}, /* approximately 10 degrees */
  {"grazing-angle",4.00,1.071797,.02},/* approximately 15 degrees */
  {"grazing-angle",4.00,1.455881,.02},/* approximately 20 degrees */
  {"grazing-range",4.00,.02,.02},
  {"grazing-range",6.00,.02,.02},
  {"grazing-range",10.00,.02,.02},
  {"interface",.50,0.,0.},
  {"interface",1.00,0.,0.},
  {"interface",2.00,0.,0.},
  {"interface",3.00,0.,0.},
  {"interface",3.10,0.,0.},
  {"axis",.01,3.00,.10},
  {"axis",.10,3.00,.10},
  {"axis",.50,3.00,.10}
};

static double seconds(void)
{
  struct timespec now;
  if(clock_gettime(CLOCK_MONOTONIC,&now)) return 0.;
  return now.tv_sec+1.e-9*now.tv_nsec;
}

static double scaled_error(const double complex a[5],
    const double complex b[5])
{
  double numerator=0.,scale=0.;
  for(int component=0;component<5;++component) {
    numerator=fmax(numerator,cabs(a[component]-b[component]));
    scale=fmax(scale,cabs(b[component]));
  }
  return numerator/fmax(scale,1.e-300);
}

static int direct_selected(double frequency,double complex lower_k,
    double complex upper_k,double rho,double z,double depth,
    unsigned int panels,unsigned int order,double cutoff,
    double complex value[5])
{
  nec4_ground_direct_result direct;
  nec4_ground_interface_result interface;
  if(z+depth==0.) {
    if(!nec4_ground_direct_cross_interface(frequency,lower_k,upper_k,rho,
        panels,order,cutoff,&interface)) return 0;
    memcpy(value,interface.field.total,sizeof(interface.field.total));
    return 1;
  }
  if(rho>0.&&z+depth<=.5*rho) {
    if(!nec4_ground_direct_cross_hankel(frequency,lower_k,upper_k,rho,z,depth,
        panels,order,cutoff,&direct)) return 0;
  } else if(!nec4_ground_direct_cross_bessel(frequency,lower_k,upper_k,rho,z,
      depth,order,cutoff,&direct)) return 0;
  memcpy(value,direct.value,sizeof(direct.value));
  return 1;
}

static int direct_alternate(double frequency,double complex lower_k,
    double complex upper_k,double rho,double z,double depth,
    unsigned int panels,unsigned int order,double cutoff,
    double complex value[5])
{
  nec4_ground_direct_result direct;
  if(z+depth==0.) return 0;
  if(rho>0.&&z+depth<=.5*rho) {
    if(!nec4_ground_direct_cross_bessel(frequency,lower_k,upper_k,rho,z,depth,
        order,cutoff,&direct)) return 0;
  } else if(!nec4_ground_direct_cross_hankel(frequency,lower_k,upper_k,rho,z,
      depth,panels,order,cutoff,&direct)) return 0;
  memcpy(value,direct.value,sizeof(direct.value));
  return 1;
}

static void usage(const char *program)
{
  fprintf(stderr,"usage: %s [--material NAME] [--output PATH]\n",program);
}

int main(int argc,char **argv)
{
  const char *material_filter=NULL,*output_path=NULL;
  FILE *output=stdout;
  int argument;
  for(argument=1;argument<argc;++argument) {
    if(!strcmp(argv[argument],"--material")&&argument+1<argc)
      material_filter=argv[++argument];
    else if(!strcmp(argv[argument],"--output")&&argument+1<argc)
      output_path=argv[++argument];
    else { usage(argv[0]);return 2; }
  }
  if(output_path) {
    output=fopen(output_path,"w");
    if(!output) { perror(output_path);return 1; }
  }
  fprintf(output,"material,family,rho_lambda,z_lambda,depth_lg,"
      "method,grazing_angle_deg,grazing_risk,medium_ok,fine_ok,alternate_ok,"
      "medium_fine_error,alternate_fine_error,production_fine_error,"
      "medium_seconds,fine_seconds,alternate_seconds");
  for(int component=0;component<5;++component)
    fprintf(output,",fine_c%d_real,fine_c%d_imag",component,component);
  fputc('\n',output);
  for(size_t m=0;m<sizeof(materials)/sizeof(materials[0]);++m) {
    const material *medium=&materials[m];
    double omega,wavelength,lg;
    double complex epsilon,upper_k,lower_k;
    if(material_filter&&strcmp(material_filter,medium->name)) continue;
    omega=2.*M_PI*medium->frequency_hz;wavelength=C0/medium->frequency_hz;
    epsilon=medium->epsilon_r-I*medium->conductivity/(omega*EPS0);
    upper_k=2.*M_PI/wavelength;lower_k=upper_k*csqrt(epsilon);
    lg=fmin(cabs(upper_k/lower_k),.3)*wavelength;
    for(size_t p=0;p<sizeof(witnesses)/sizeof(witnesses[0]);++p) {
      const witness *point=&witnesses[p];
      double rho=point->rho_lambda*wavelength,z=point->z_lambda*wavelength;
      double depth=point->depth_lg*lg;
      double complex coarse[5]={0},fine[5]={0},alternate[5]={0},production[5]={0};
      nec4_ground_dispatch route;
      nec4_ground_asymptotic_result asymptotic;
      double before,coarse_time,fine_time,alternate_time;
      int coarse_ok,fine_ok,alternate_ok,production_ok=0;
      before=seconds();
      coarse_ok=direct_selected(medium->frequency_hz,lower_k,upper_k,rho,z,
          depth,3,48,80.,coarse);
      coarse_time=seconds()-before;
      before=seconds();
      fine_ok=direct_selected(medium->frequency_hz,lower_k,upper_k,rho,z,
          depth,5,80,120.,fine);
      fine_time=seconds()-before;
      before=seconds();
      alternate_ok=direct_alternate(medium->frequency_hz,lower_k,upper_k,rho,z,
          depth,5,80,120.,alternate);
      alternate_time=seconds()-before;
      if(!nec4_ground_route(rho,z,-depth,lower_k,upper_k,wavelength,&route))
        return 1;
      if(route.method==NEC4_GROUND_ASYMPTOTIC) {
        production_ok=nec4_ground_asymptotic_second_uniform(
          medium->frequency_hz,lower_k,upper_k,rho,z,depth,&asymptotic);
        if(production_ok) memcpy(production,asymptotic.value,sizeof(production));
      } else if(coarse_ok) {
        /* All non-asymptotic witnesses in this harness lie in the accepted
           low-height direct strip.  This represents current production. */
        memcpy(production,coarse,sizeof(production));production_ok=1;
      }
      fprintf(output,"%s,%s,%.17g,%.17g,%.17g,%s,%.17g,%d,%d,%d,%d,",
          medium->name,point->family,point->rho_lambda,point->z_lambda,
          point->depth_lg,nec4_ground_method_name(route.method),
          route.grazing_angle_degrees,route.grazing_risk,coarse_ok,fine_ok,
          alternate_ok);
      fprintf(output,"%.17g,%.17g,%.17g,%.9g,%.9g,%.9g",
          coarse_ok&&fine_ok?scaled_error(coarse,fine):NAN,
          alternate_ok&&fine_ok?scaled_error(alternate,fine):NAN,
          production_ok&&fine_ok?scaled_error(production,fine):NAN,
          coarse_time,fine_time,alternate_time);
      for(int component=0;component<5;++component)
        fprintf(output,",%.17g,%.17g",creal(fine[component]),cimag(fine[component]));
      fputc('\n',output);fflush(output);
    }
  }
  if(material_filter) {
    int found=0;
    for(size_t m=0;m<sizeof(materials)/sizeof(materials[0]);++m)
      if(!strcmp(material_filter,materials[m].name)) found=1;
    if(!found) { fprintf(stderr,"unknown material: %s\n",material_filter);return 2; }
  }
  if(output_path&&fclose(output)) { perror(output_path);return 1; }
  return 0;
}
