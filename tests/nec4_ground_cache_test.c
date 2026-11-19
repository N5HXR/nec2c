#define _POSIX_C_SOURCE 200809L

#include "nec4_ground_cache.h"

#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

int main(void)
{
  char path[]="/tmp/nec2c-nec4-cache-XXXXXX";
  int descriptor=mkstemp(path);
  FILE *stream;
  long offset;
  unsigned char byte;
  nec4_ground_cache_key key={
    .equation_version=NEC4_GROUND_EQUATION_VERSION,
    .precision_bits=53,.table_kind=1,.component_count=1,
    .frequency_hz=7.1e6,.epsilon_real=13.,.epsilon_imag=-.0123,
    .wavelength0_m=299792458./7.1e6,.lg_m=2.5,.lmax_m=12.5,
    .rho_start_m=0.,.rho_step_m=.5,
    .observer_start_m=0.,.observer_step_m=.5,
    .source_depth_start_m=0.,.source_depth_step_m=.5,
    .rho_count=2,.observer_count=2,.source_depth_count=2,
    .value_count=8};
  nec4_ground_cache_key wrong;
  double complex values[8]={1.+2.*I,-3.+4.*I,5.5-.25*I,-0.0+0.0*I,
    8.-9.*I,10.+11.*I,-12.-13.*I,14.+15.*I};
  double complex loaded[8]={0};
  assert(descriptor>=0);close(descriptor);unlink(path);
  assert(nec4_ground_cache_load(path,&key,loaded)==NEC4_CACHE_NOT_FOUND);
  assert(nec4_ground_cache_store(path,&key,values)==NEC4_CACHE_OK);
  assert(nec4_ground_cache_load(path,&key,loaded)==NEC4_CACHE_OK);
  assert(!memcmp(values,loaded,sizeof(values)));
  wrong=key;wrong.frequency_hz=7.2e6;
  assert(nec4_ground_cache_load(path,&wrong,loaded)==NEC4_CACHE_KEY_MISMATCH);
  wrong=key;wrong.rho_step_m=.25;
  assert(nec4_ground_cache_load(path,&wrong,loaded)==NEC4_CACHE_KEY_MISMATCH);
  wrong=key;wrong.rho_count=3;
  assert(nec4_ground_cache_load(path,&wrong,loaded)==NEC4_CACHE_FORMAT_ERROR);
  stream=fopen(path,"r+b");assert(stream);
  assert(!fseek(stream,-(long)sizeof(double),SEEK_END));
  offset=ftell(stream);assert(offset>=0);
  assert(fread(&byte,1,1,stream)==1);byte^=0x80;
  assert(!fseek(stream,offset,SEEK_SET));assert(fwrite(&byte,1,1,stream)==1);
  assert(!fclose(stream));
  assert(nec4_ground_cache_load(path,&key,loaded)==NEC4_CACHE_CHECKSUM_ERROR);

  assert(nec4_ground_cache_store(path,&key,values)==NEC4_CACHE_OK);
  stream=fopen(path,"ab");assert(stream);assert(fputc(0x5a,stream)!=EOF);
  assert(!fclose(stream));
  assert(nec4_ground_cache_load(path,&key,loaded)==NEC4_CACHE_FORMAT_ERROR);

  assert(nec4_ground_cache_store(path,&key,values)==NEC4_CACHE_OK);
  {
    struct stat status;
    assert(!stat(path,&status));assert(status.st_size>1);
    assert(!truncate(path,status.st_size-1));
  }
  assert(nec4_ground_cache_load(path,&key,loaded)==NEC4_CACHE_FORMAT_ERROR);

  wrong=key;wrong.value_count--;
  assert(nec4_ground_cache_store(path,&wrong,values)==NEC4_CACHE_FORMAT_ERROR);
  unlink(path);
  return 0;
}
