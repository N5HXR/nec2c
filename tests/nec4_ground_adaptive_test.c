#define _POSIX_C_SOURCE 200809L

#include "nec4_ground_adaptive.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct { int oscillatory; } fixture;

static int sample(void *opaque,double x,double y,double z,double complex value[5])
{
  fixture *context=opaque;
  for(unsigned int component=0;component<5;++component) {
    double scale=component+1.;
    double real=1.+scale*x+(.5+scale)*y+(.25+scale)*z+
      .2*x*y-.1*y*z+.05*x*z+.03*x*y*z+
      .02*x*x*x-.015*y*y*y+.01*z*z*z;
    double imaginary=-.4+.1*scale*x-.2*y+.3*z+.04*x*x*y;
    if(context->oscillatory) {
      real+=.05*sin((3.+component)*x*y+.7*z);
      imaginary+=.03*cos(2.3*x+(1.+component)*y*z);
    }
    value[component]=real+I*imaginary;
  }
  return 1;
}

static double error(const double complex actual[5],const double complex expected[5])
{
  double difference=0.,scale=0.;
  for(unsigned int component=0;component<5;++component) {
    difference=fmax(difference,cabs(actual[component]-expected[component]));
    scale=fmax(scale,cabs(expected[component]));
  }
  return difference/scale;
}

int main(void)
{
  nec4_ground_adaptive serial,parallel,curve,limited,loaded;
  nec4_ground_adaptive_bounds volume={{0.,0.,0.},{1.,1.,1.}};
  nec4_ground_adaptive_bounds line={{0.,.25,.1},{1.,.25,.1}};
  fixture cubic={0},oscillatory={1};
  double worst=0.;
  char cache[]="/tmp/nec4-adaptive-cache-XXXXXX";
  const char identity[]="synthetic-adaptive-v1";
  nec4_ground_adaptive_init(&serial);nec4_ground_adaptive_init(&parallel);
  nec4_ground_adaptive_init(&curve);
  assert(nec4_ground_adaptive_build(&serial,&volume,1.e-12,2,0,1,sample,
                                     &cubic));
  assert(serial.stats.cells==1&&serial.stats.leaves==1&&
         serial.stats.accepted_leaves==1&&serial.stats.unresolved_leaves==0);
  assert(serial.stats.fitted_samples==64&&serial.stats.withheld_samples==27);
  assert(serial.stats.worst_withheld_scaled_error<1.e-12);
  assert(nec4_ground_adaptive_build(&parallel,&volume,1.e-4,5,0,4,sample,
                                     &oscillatory));
  assert(nec4_ground_adaptive_build(&serial,&volume,1.e-4,5,0,1,sample,
                                     &oscillatory));
  assert(!memcmp(&serial.stats,&parallel.stats,sizeof(serial.stats)));
  assert(serial.stats.cells>1&&serial.stats.unresolved_leaves==0);
  for(unsigned int index=1;index<=343;++index) {
    double x=fmod(index*.6180339887498948,1.);
    double y=fmod(index*.4142135623730950,1.);
    double z=fmod(index*.7320508075688772,1.);
    double complex expected[5],one[5],many[5];
    assert(sample(&oscillatory,x,y,z,expected));
    assert(nec4_ground_adaptive_evaluate(&serial,x,y,z,one));
    assert(nec4_ground_adaptive_evaluate(&parallel,x,y,z,many));
    assert(!memcmp(one,many,sizeof(one)));
    worst=fmax(worst,error(one,expected));
  }
  assert(worst<1.e-4);
  assert(!nec4_ground_adaptive_evaluate(&serial,1.01,.5,.5,
                                         (double complex[5]){0}));
  nec4_ground_adaptive_init(&limited);
  assert(nec4_ground_adaptive_build(&limited,&volume,1.e-8,5,2,4,sample,
                                     &oscillatory));
  assert(limited.stats.cells==2&&limited.stats.unresolved_leaves>0);
  assert(!nec4_ground_adaptive_evaluate(&limited,.9,.9,.9,
                                         (double complex[5]){0}));
  assert(nec4_ground_adaptive_build(&curve,&line,1.e-4,5,0,4,sample,
                                     &oscillatory));
  assert(curve.stats.fitted_samples>=4&&curve.stats.withheld_samples>=3&&
         curve.stats.unresolved_leaves==0);
  for(unsigned int index=0;index<=100;++index) {
    double x=index/100.;double complex expected[5],actual[5];
    assert(sample(&oscillatory,x,.25,.1,expected));
    assert(nec4_ground_adaptive_evaluate(&curve,x,.25,.1,actual));
    assert(error(actual,expected)<1.e-4);
  }
  {
    int descriptor=mkstemp(cache);double complex original[5],restored[5];
    FILE *stream;
    assert(descriptor>=0);close(descriptor);assert(!unlink(cache));
    assert(nec4_ground_adaptive_store(&serial,cache,identity,sizeof(identity)));
    nec4_ground_adaptive_init(&loaded);
    assert(!nec4_ground_adaptive_load(&loaded,cache,"wrong",6));
    assert(nec4_ground_adaptive_load(&loaded,cache,identity,sizeof(identity)));
    assert(!memcmp(&serial.stats,&loaded.stats,sizeof(serial.stats)));
    assert(nec4_ground_adaptive_evaluate(&serial,.37,.41,.23,original));
    assert(nec4_ground_adaptive_evaluate(&loaded,.37,.41,.23,restored));
    assert(!memcmp(original,restored,sizeof(original)));
    assert(!nec4_ground_adaptive_load(&loaded,cache,"wrong",6));
    stream=fopen(cache,"ab");assert(stream);assert(fputc(0,stream)!=EOF);
    assert(!fclose(stream));
    assert(!nec4_ground_adaptive_load(&loaded,cache,identity,sizeof(identity)));
    assert(!unlink(cache));nec4_ground_adaptive_destroy(&loaded);
  }
  printf("NEC4 adaptive cubic cells passed: %zu cells, depth %u, "
         "worst heldout %.3e, independent %.3e\n",serial.stats.cells,
         serial.stats.deepest_level,serial.stats.worst_withheld_scaled_error,
         worst);
  nec4_ground_adaptive_destroy(&curve);nec4_ground_adaptive_destroy(&parallel);
  nec4_ground_adaptive_destroy(&limited);
  nec4_ground_adaptive_destroy(&serial);
  return 0;
}
