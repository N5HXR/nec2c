#define _POSIX_C_SOURCE 200809L

#include "nec4_ground_direct.h"

#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define C0 299792458.0
#define EPS0 8.8541878128e-12

typedef struct {
  double rho,z,depth;
  double complex value[NEC4_GROUND_COMPONENT_COUNT];
  int valid,heldout;
} sample;

static double halton(unsigned int index,unsigned int base)
{
  double fraction=1.,value=0.;
  while(index) { fraction/=base;value+=fraction*(index%base);index/=base; }
  return value;
}

static void usage(const char *program)
{
  fprintf(stderr,"usage: %s --epsilon ER --conductivity S/M [--frequency HZ] "
      "[--workers N] [--level N] [--rho-level N] [--z-level N] "
      "[--depth-level N] [--output PATH]\n",program);
}

static int evaluate(double frequency,double complex k1,double complex k2,
    sample *point)
{
  nec4_ground_direct_result direct;
  nec4_ground_interface_result interface;
  if(point->rho==0.&&point->z==0.&&point->depth==0.) {
    memset(point->value,0,sizeof(point->value));return 1;
  }
  if(point->z+point->depth==0.) {
    if(!nec4_ground_direct_cross_interface(frequency,k1,k2,point->rho,
                                           3,48,80.,&interface)) return 0;
    memcpy(point->value,interface.field.total,sizeof(point->value));return 1;
  }
  if(point->rho>0.&&point->z+point->depth<=.5*point->rho) {
    if(!nec4_ground_direct_cross_hankel(frequency,k1,k2,point->rho,point->z,
                                        point->depth,3,48,80.,&direct)) return 0;
  } else if(!nec4_ground_direct_cross_bessel(frequency,k1,k2,point->rho,
                                              point->z,point->depth,48,80.,
                                              &direct)) return 0;
  memcpy(point->value,direct.value,sizeof(point->value));return 1;
}

int main(int argc,char **argv)
{
  double frequency=7.1e6,epsilon_r=NAN,conductivity=NAN;
  int workers=0;
  int level=0;
  int rho_level=-1,z_level=-1,depth_level=-1;
  const char *output_path=NULL;
  FILE *output=stdout;
  int nr,nz,nd,multiplier;
  size_t train_count,heldout_count,count;
  sample *points;
  double omega,wavelength,lg;
  double complex epsilon,k1,k2;
  size_t cursor=0,index;
  int argument;
  for(argument=1;argument<argc;argument++) {
    if(!strcmp(argv[argument],"--epsilon")&&argument+1<argc)
      epsilon_r=strtod(argv[++argument],NULL);
    else if(!strcmp(argv[argument],"--conductivity")&&argument+1<argc)
      conductivity=strtod(argv[++argument],NULL);
    else if(!strcmp(argv[argument],"--frequency")&&argument+1<argc)
      frequency=strtod(argv[++argument],NULL);
    else if(!strcmp(argv[argument],"--workers")&&argument+1<argc)
      workers=(int)strtol(argv[++argument],NULL,10);
    else if(!strcmp(argv[argument],"--level")&&argument+1<argc)
      level=(int)strtol(argv[++argument],NULL,10);
    else if(!strcmp(argv[argument],"--rho-level")&&argument+1<argc)
      rho_level=(int)strtol(argv[++argument],NULL,10);
    else if(!strcmp(argv[argument],"--z-level")&&argument+1<argc)
      z_level=(int)strtol(argv[++argument],NULL,10);
    else if(!strcmp(argv[argument],"--depth-level")&&argument+1<argc)
      depth_level=(int)strtol(argv[++argument],NULL,10);
    else if(!strcmp(argv[argument],"--output")&&argument+1<argc)
      output_path=argv[++argument];
    else { usage(argv[0]);return 2; }
  }
  if(rho_level<0) rho_level=level;
  if(z_level<0) z_level=level;
  if(depth_level<0) depth_level=level;
  if(!(frequency>0.)||!(epsilon_r>0.)||conductivity<0.||level<0||level>3||
     rho_level<0||rho_level>3||z_level<0||z_level>4||
     depth_level<0||depth_level>3||
     !isfinite(frequency)||!isfinite(epsilon_r)||!isfinite(conductivity)) {
    usage(argv[0]);return 2;
  }
  omega=2.*M_PI*frequency;wavelength=C0/frequency;k2=2.*M_PI/wavelength;
  epsilon=epsilon_r-I*conductivity/(omega*EPS0);k1=k2*csqrt(epsilon);
  lg=fmin(cabs(k2/k1),.3)*wavelength;
  multiplier=1<<rho_level;nr=12*multiplier+1;
  multiplier=1<<z_level;nz=8*multiplier+1;
  multiplier=1<<depth_level;nd=4*multiplier+1;
  train_count=(size_t)nr*nz*nd;
  heldout_count=4096;
  count=train_count+heldout_count;
  points=calloc(count,sizeof(*points));
  if(!points) { perror("calloc");return 1; }
  for(int ir=0;ir<nr;ir++) for(int iz=0;iz<nz;iz++) for(int id=0;id<nd;id++) {
    points[cursor].rho=3.*wavelength*ir/(nr-1);
    points[cursor].z=2.*wavelength*iz/(nz-1);
    points[cursor].depth=.8*lg*id/(nd-1);
    cursor++;
  }
  for(unsigned int heldout=1;heldout<=heldout_count;heldout++) {
    points[cursor].rho=3.*wavelength*halton(heldout,2);
    points[cursor].z=2.*wavelength*halton(heldout,3);
    points[cursor].depth=.8*lg*halton(heldout,5);
    points[cursor].heldout=1;cursor++;
  }
  if(cursor!=count) { free(points);return 1; }
#ifdef _OPENMP
  if(workers>0) omp_set_num_threads(workers);
#pragma omp parallel for schedule(dynamic,1)
#endif
  for(index=0;index<count;index++)
    points[index].valid=evaluate(frequency,k1,k2,&points[index]);
  for(index=0;index<count;index++) if(!points[index].valid) {
    fprintf(stderr,"oracle failed at %.17g %.17g %.17g\n",
            points[index].rho,points[index].z,points[index].depth);
    free(points);return 1;
  }
  if(output_path) {
    output=fopen(output_path,"w");
    if(!output) { perror(output_path);free(points);return 1; }
  }
  fprintf(output,"# frequency=%.17g epsilon=%.17g conductivity=%.17g "
      "wavelength=%.17g lg=%.17g nr=%d nz=%d nd=%d\n",
      frequency,epsilon_r,conductivity,wavelength,lg,nr,nz,nd);
  fprintf(output,"kind,rho,z,depth");
  for(int component=0;component<NEC4_GROUND_COMPONENT_COUNT;component++)
    fprintf(output,",c%d_real,c%d_imag",component,component);
  fputc('\n',output);
  for(index=0;index<count;index++) {
    fprintf(output,"%s,%.17g,%.17g,%.17g",
            points[index].heldout?"heldout":"train",points[index].rho,
            points[index].z,points[index].depth);
    for(int component=0;component<NEC4_GROUND_COMPONENT_COUNT;component++)
      fprintf(output,",%.17g,%.17g",creal(points[index].value[component]),
              cimag(points[index].value[component]));
    fputc('\n',output);
  }
  if(output_path&&fclose(output)) { perror(output_path);free(points);return 1; }
  free(points);return 0;
}
