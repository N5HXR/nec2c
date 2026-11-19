#define _POSIX_C_SOURCE 200809L

#include "nec4_ground_runtime.h"

#include <assert.h>
#include <complex.h>
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int sample(void *opaque,double rho,double observer,double depth,
    double complex value[5])
{
  double bias=*(double *)opaque;
  for(int component=0;component<5;++component)
    value[component]=bias+component+.2*rho-.3*observer+.4*depth+
      I*(-bias+.1*component-.5*rho+.6*observer-.7*depth);
  return 1;
}

static void remove_directory(const char *path)
{
  DIR *directory=opendir(path);struct dirent *entry;char file[4096];
  assert(directory);
  while((entry=readdir(directory)))
  {
    if(!strcmp(entry->d_name,".")||!strcmp(entry->d_name,"..")) continue;
    assert(snprintf(file,sizeof(file),"%s/%s",path,entry->d_name)>0);
    assert(!unlink(file));
  }
  closedir(directory);assert(!rmdir(path));
}

int main(void)
{
  char directory[]="/tmp/nec2c-nec4-runtime-XXXXXX";
  const double frequency=7.1e6,wavelength=299792458./frequency;
  const double complex epsilon=13.-.1*I;
  const double complex upper_k=2.*acos(-1.)/wavelength;
  const double complex lower_k=upper_k*csqrt(epsilon);
  double bias=1.25;
  nec4_ground_runtime first,second,changed;
  nec4_ground_runtime_options options={
    .cache_directory=directory,.workers=4,.sampler=sample,
    .sampler_context=&bias,.direct_options=NULL};
  assert(mkdtemp(directory));
  nec4_ground_runtime_init(&first);nec4_ground_runtime_init(&second);
  nec4_ground_runtime_init(&changed);
  assert(nec4_ground_runtime_prepare(&first,frequency,epsilon,lower_k,upper_k,
      &options));
  assert(first.ready&&first.cache_hits==0&&first.cache_misses==4&&
      first.generated_tables==4);
  assert(fabs(first.wavelength0_m-wavelength)<1.e-14*wavelength);
  assert(fabs(first.lg_m-fmin(cabs(upper_k/lower_k),.3)*wavelength)<
      1.e-14*wavelength);
  assert(fabs(first.lmax_m-fmin(15.*cabs(upper_k/lower_k),1.)*wavelength)<
      1.e-14*wavelength);
  assert(nec4_ground_runtime_matrix_context(&first));
  {
    double complex field[5];nec4_ground_dispatch route;
    nec4_ground_direct_result direct;
    nec4_ground_interface_result interface;
    assert(nec4_ground_runtime_evaluate(&first,.1,.1,.1,field,&route));
    assert(route.method==NEC4_GROUND_LINEAR_TABLE);
    assert(nec4_ground_runtime_evaluate(&first,first.wavelength0_m,
      .5*first.wavelength0_m,.5*first.lg_m,field,&route));
    assert(route.method==NEC4_GROUND_MODEL_FIT);
    assert(!route.used_direct_fallback);
    {
      double rho=first.wavelength0_m,z=.5*first.wavelength0_m,
        depth=.5*first.lg_m,radius=hypot(rho,z+depth);
      double complex factor=radius*radius*cexp(I*(upper_k*hypot(rho,z)+
                                                  lower_k*depth));
      double complex expected[5];
      assert(sample(&bias,rho,z,depth,expected));
      for(int component=0;component<5;++component)
        assert(cabs(field[component]*factor-expected[component])<2.e-12);
    }
    assert(nec4_ground_runtime_evaluate(&first,first.wavelength0_m,
      .05*first.lg_m,.05*first.lg_m,field,&route));
    assert(route.method==NEC4_GROUND_MODEL_FIT&&route.used_direct_fallback);
    assert(nec4_ground_direct_cross_hankel(frequency,lower_k,upper_k,
      first.wavelength0_m,.05*first.lg_m,.05*first.lg_m,3,48,80.,&direct));
    assert(!memcmp(field,direct.value,sizeof(field)));
    assert(nec4_ground_runtime_evaluate(&first,.5*first.wavelength0_m,
      0.,0.,field,&route));
    assert(route.method==NEC4_GROUND_MODEL_FIT&&route.used_direct_fallback);
    assert(nec4_ground_direct_cross_interface(frequency,lower_k,upper_k,
      .5*first.wavelength0_m,3,48,80.,&interface));
    assert(!memcmp(field,interface.field.total,sizeof(field)));
    assert(nec4_ground_runtime_evaluate(&first,4.*first.wavelength0_m,
      .8*first.wavelength0_m,.8*first.wavelength0_m,field,&route));
    assert(route.method==NEC4_GROUND_ASYMPTOTIC);
    for(int component=0;component<5;++component)
      assert(isfinite(creal(field[component]))&&isfinite(cimag(field[component])));
  }
  assert(nec4_ground_runtime_prepare(&second,frequency,epsilon,lower_k,upper_k,
      &options));
  assert(second.ready&&second.cache_hits==4&&second.cache_misses==0&&
      second.generated_tables==0&&second.payload_bytes==first.payload_bytes);
  for(int table=0;table<3;++table)
    assert(!memcmp(first.table[table].values,second.table[table].values,
      first.table[table].value_count*sizeof(*first.table[table].values)));
  {
    const nec4_ground_matrix_context *a=nec4_ground_runtime_matrix_context(&first);
    const nec4_ground_matrix_context *b=nec4_ground_runtime_matrix_context(&second);
    double complex left[5],right[5];
    assert(a&&b&&a->cross(a->cross_context,.2*first.lg_m,.3*first.lg_m,
      .7*first.lg_m,left));
    assert(b->cross(b->cross_context,.2*second.lg_m,.3*second.lg_m,
      .7*second.lg_m,right));
    assert(!memcmp(left,right,sizeof(left)));
  }
  assert(nec4_ground_runtime_prepare(&changed,frequency*1.01,epsilon,
      lower_k*1.01,upper_k*1.01,&options));
  assert(changed.cache_hits==0&&changed.generated_tables==4);
  printf("NEC4 frequency table runtime passed: %zu bytes, first 0/4 hits, "
      "repeat 4/4 hits, changed frequency 0/4 hits\n",first.payload_bytes);
  nec4_ground_runtime_destroy(&changed);nec4_ground_runtime_destroy(&second);
  nec4_ground_runtime_destroy(&first);remove_directory(directory);
  return 0;
}
