#define _POSIX_C_SOURCE 200809L

#include "nec4_ground_cache.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define CACHE_ENDIAN_MARKER 0x01020304u
#define CACHE_MAGIC "N4GTBL1\0"
#define CACHE_MAGIC_SIZE 8u
#define FNV_OFFSET UINT64_C(14695981039346656037)
#define FNV_PRIME UINT64_C(1099511628211)

static uint64_t checksum_bytes(uint64_t checksum,const void *data,size_t size)
{
  const unsigned char *bytes=data;
  size_t index;
  for(index=0;index<size;++index)
  {
    checksum^=bytes[index];
    checksum*=FNV_PRIME;
  }
  return checksum;
}

static uint64_t checksum_key(uint64_t checksum,const nec4_ground_cache_key *key)
{
  checksum=checksum_bytes(checksum,&key->equation_version,
                          sizeof(key->equation_version));
  checksum=checksum_bytes(checksum,&key->precision_bits,
                          sizeof(key->precision_bits));
  checksum=checksum_bytes(checksum,&key->table_kind,sizeof(key->table_kind));
  checksum=checksum_bytes(checksum,&key->component_count,
                          sizeof(key->component_count));
  checksum=checksum_bytes(checksum,&key->frequency_hz,sizeof(key->frequency_hz));
  checksum=checksum_bytes(checksum,&key->epsilon_real,sizeof(key->epsilon_real));
  checksum=checksum_bytes(checksum,&key->epsilon_imag,sizeof(key->epsilon_imag));
  checksum=checksum_bytes(checksum,&key->wavelength0_m,
                          sizeof(key->wavelength0_m));
  checksum=checksum_bytes(checksum,&key->lg_m,sizeof(key->lg_m));
  checksum=checksum_bytes(checksum,&key->lmax_m,sizeof(key->lmax_m));
  checksum=checksum_bytes(checksum,&key->rho_start_m,sizeof(key->rho_start_m));
  checksum=checksum_bytes(checksum,&key->rho_step_m,sizeof(key->rho_step_m));
  checksum=checksum_bytes(checksum,&key->observer_start_m,
                          sizeof(key->observer_start_m));
  checksum=checksum_bytes(checksum,&key->observer_step_m,
                          sizeof(key->observer_step_m));
  checksum=checksum_bytes(checksum,&key->source_depth_start_m,
                          sizeof(key->source_depth_start_m));
  checksum=checksum_bytes(checksum,&key->source_depth_step_m,
                          sizeof(key->source_depth_step_m));
  checksum=checksum_bytes(checksum,&key->rho_count,sizeof(key->rho_count));
  checksum=checksum_bytes(checksum,&key->observer_count,
                          sizeof(key->observer_count));
  checksum=checksum_bytes(checksum,&key->source_depth_count,
                          sizeof(key->source_depth_count));
  return checksum_bytes(checksum,&key->value_count,sizeof(key->value_count));
}

static int multiply_u64(uint64_t a,uint64_t b,uint64_t *result)
{
  if(a&&b>UINT64_MAX/a) return 0;
  *result=a*b;
  return 1;
}

static int valid_key(const nec4_ground_cache_key *key)
{
  uint64_t points,values;
  if(!key||key->equation_version==0||key->precision_bits==0||
     key->table_kind==0||key->component_count==0||
     !(key->frequency_hz>0.)||!(key->wavelength0_m>0.)||
     !(key->lg_m>0.)||!(key->lmax_m>=key->lg_m)||
     !isfinite(key->frequency_hz)||!isfinite(key->epsilon_real)||
     !isfinite(key->epsilon_imag)||!isfinite(key->wavelength0_m)||
     !isfinite(key->lg_m)||!isfinite(key->lmax_m)||
     !isfinite(key->rho_start_m)||!isfinite(key->rho_step_m)||
     !isfinite(key->observer_start_m)||!isfinite(key->observer_step_m)||
     !isfinite(key->source_depth_start_m)||
     !isfinite(key->source_depth_step_m)||
     key->rho_count<2||key->observer_count<2||key->source_depth_count<2||
     !(key->rho_step_m>0.)||!(key->observer_step_m>0.)||
     !(key->source_depth_step_m>0.)) return 0;
  return multiply_u64(key->rho_count,key->observer_count,&points)&&
    multiply_u64(points,key->source_depth_count,&points)&&
    multiply_u64(points,key->component_count,&values)&&
    values==key->value_count&&values<=SIZE_MAX/(2u*sizeof(double));
}

uint64_t nec4_ground_cache_key_hash(const nec4_ground_cache_key *key)
{
  if(!valid_key(key)) return 0;
  return checksum_key(FNV_OFFSET,key);
}

static int write_exact(FILE *stream,const void *value,size_t size)
{ return fwrite(value,1,size,stream)==size; }

static int read_exact(FILE *stream,void *value,size_t size)
{ return fread(value,1,size,stream)==size; }

static int write_key(FILE *stream,const nec4_ground_cache_key *key)
{
  return write_exact(stream,&key->equation_version,sizeof(key->equation_version))&&
    write_exact(stream,&key->precision_bits,sizeof(key->precision_bits))&&
    write_exact(stream,&key->table_kind,sizeof(key->table_kind))&&
    write_exact(stream,&key->component_count,sizeof(key->component_count))&&
    write_exact(stream,&key->frequency_hz,sizeof(key->frequency_hz))&&
    write_exact(stream,&key->epsilon_real,sizeof(key->epsilon_real))&&
    write_exact(stream,&key->epsilon_imag,sizeof(key->epsilon_imag))&&
    write_exact(stream,&key->wavelength0_m,sizeof(key->wavelength0_m))&&
    write_exact(stream,&key->lg_m,sizeof(key->lg_m))&&
    write_exact(stream,&key->lmax_m,sizeof(key->lmax_m))&&
    write_exact(stream,&key->rho_start_m,sizeof(key->rho_start_m))&&
    write_exact(stream,&key->rho_step_m,sizeof(key->rho_step_m))&&
    write_exact(stream,&key->observer_start_m,sizeof(key->observer_start_m))&&
    write_exact(stream,&key->observer_step_m,sizeof(key->observer_step_m))&&
    write_exact(stream,&key->source_depth_start_m,
                sizeof(key->source_depth_start_m))&&
    write_exact(stream,&key->source_depth_step_m,
                sizeof(key->source_depth_step_m))&&
    write_exact(stream,&key->rho_count,sizeof(key->rho_count))&&
    write_exact(stream,&key->observer_count,sizeof(key->observer_count))&&
    write_exact(stream,&key->source_depth_count,
                sizeof(key->source_depth_count))&&
    write_exact(stream,&key->value_count,sizeof(key->value_count));
}

static int read_key(FILE *stream,nec4_ground_cache_key *key)
{
  return read_exact(stream,&key->equation_version,sizeof(key->equation_version))&&
    read_exact(stream,&key->precision_bits,sizeof(key->precision_bits))&&
    read_exact(stream,&key->table_kind,sizeof(key->table_kind))&&
    read_exact(stream,&key->component_count,sizeof(key->component_count))&&
    read_exact(stream,&key->frequency_hz,sizeof(key->frequency_hz))&&
    read_exact(stream,&key->epsilon_real,sizeof(key->epsilon_real))&&
    read_exact(stream,&key->epsilon_imag,sizeof(key->epsilon_imag))&&
    read_exact(stream,&key->wavelength0_m,sizeof(key->wavelength0_m))&&
    read_exact(stream,&key->lg_m,sizeof(key->lg_m))&&
    read_exact(stream,&key->lmax_m,sizeof(key->lmax_m))&&
    read_exact(stream,&key->rho_start_m,sizeof(key->rho_start_m))&&
    read_exact(stream,&key->rho_step_m,sizeof(key->rho_step_m))&&
    read_exact(stream,&key->observer_start_m,sizeof(key->observer_start_m))&&
    read_exact(stream,&key->observer_step_m,sizeof(key->observer_step_m))&&
    read_exact(stream,&key->source_depth_start_m,
               sizeof(key->source_depth_start_m))&&
    read_exact(stream,&key->source_depth_step_m,
               sizeof(key->source_depth_step_m))&&
    read_exact(stream,&key->rho_count,sizeof(key->rho_count))&&
    read_exact(stream,&key->observer_count,sizeof(key->observer_count))&&
    read_exact(stream,&key->source_depth_count,
               sizeof(key->source_depth_count))&&
    read_exact(stream,&key->value_count,sizeof(key->value_count));
}

static int keys_equal(const nec4_ground_cache_key *a,
                      const nec4_ground_cache_key *b)
{
  return a->equation_version==b->equation_version&&
    a->precision_bits==b->precision_bits&&
    a->table_kind==b->table_kind&&
    a->component_count==b->component_count&&
    a->rho_count==b->rho_count&&
    a->observer_count==b->observer_count&&
    a->source_depth_count==b->source_depth_count&&
    a->value_count==b->value_count&&
    !memcmp(&a->frequency_hz,&b->frequency_hz,sizeof(double))&&
    !memcmp(&a->epsilon_real,&b->epsilon_real,sizeof(double))&&
    !memcmp(&a->epsilon_imag,&b->epsilon_imag,sizeof(double))&&
    !memcmp(&a->wavelength0_m,&b->wavelength0_m,sizeof(double))&&
    !memcmp(&a->lg_m,&b->lg_m,sizeof(double))&&
    !memcmp(&a->lmax_m,&b->lmax_m,sizeof(double))&&
    !memcmp(&a->rho_start_m,&b->rho_start_m,sizeof(double))&&
    !memcmp(&a->rho_step_m,&b->rho_step_m,sizeof(double))&&
    !memcmp(&a->observer_start_m,&b->observer_start_m,sizeof(double))&&
    !memcmp(&a->observer_step_m,&b->observer_step_m,sizeof(double))&&
    !memcmp(&a->source_depth_start_m,&b->source_depth_start_m,sizeof(double))&&
    !memcmp(&a->source_depth_step_m,&b->source_depth_step_m,sizeof(double));
}

nec4_ground_cache_status nec4_ground_cache_store(const char *path,
    const nec4_ground_cache_key *key,const double complex *values)
{
  FILE *stream=NULL;
  char *temporary=NULL;
  size_t length,index;
  int descriptor=-1,failed=0;
  uint32_t format=NEC4_GROUND_CACHE_FORMAT_VERSION;
  uint32_t endian=CACHE_ENDIAN_MARKER;
  uint32_t double_size=(uint32_t)sizeof(double);
  uint64_t checksum=FNV_OFFSET;
  if(!path||!*path||!key||(!values&&key->value_count)) return NEC4_CACHE_IO_ERROR;
  if(!valid_key(key)) return NEC4_CACHE_FORMAT_ERROR;
  checksum=checksum_key(checksum,key);
  for(index=0;index<(size_t)key->value_count;++index)
  {
    double pair[2]={creal(values[index]),cimag(values[index])};
    checksum=checksum_bytes(checksum,pair,sizeof(pair));
  }
  length=strlen(path);
  temporary=malloc(length+12);
  if(!temporary) return NEC4_CACHE_IO_ERROR;
  memcpy(temporary,path,length);
  memcpy(temporary+length,".tmp.XXXXXX",12);
  descriptor=mkstemp(temporary);
  if(descriptor<0) { free(temporary);return NEC4_CACHE_IO_ERROR; }
  stream=fdopen(descriptor,"wb");
  if(!stream) { close(descriptor);unlink(temporary);free(temporary);
    return NEC4_CACHE_IO_ERROR; }
  failed=!write_exact(stream,CACHE_MAGIC,CACHE_MAGIC_SIZE)||
    !write_exact(stream,&format,sizeof(format))||
    !write_exact(stream,&endian,sizeof(endian))||
    !write_exact(stream,&double_size,sizeof(double_size))||
    !write_key(stream,key)||!write_exact(stream,&checksum,sizeof(checksum));
  for(index=0;!failed&&index<(size_t)key->value_count;++index)
  {
    double pair[2]={creal(values[index]),cimag(values[index])};
    failed=!write_exact(stream,pair,sizeof(pair));
  }
  if(!failed&&(fflush(stream)||fsync(fileno(stream)))) failed=1;
  if(fclose(stream)) failed=1;
  if(!failed&&rename(temporary,path)) failed=1;
  if(failed) unlink(temporary);
  free(temporary);
  return failed?NEC4_CACHE_IO_ERROR:NEC4_CACHE_OK;
}

nec4_ground_cache_status nec4_ground_cache_load(const char *path,
    const nec4_ground_cache_key *expected,double complex *values)
{
  FILE *stream;
  char magic[CACHE_MAGIC_SIZE];
  uint32_t format,endian,double_size;
  uint64_t stored_checksum,checksum;
  nec4_ground_cache_key actual;
  size_t index;
  int trailing;
  if(!path||!expected||(!values&&expected->value_count)) return NEC4_CACHE_IO_ERROR;
  if(!valid_key(expected)) return NEC4_CACHE_FORMAT_ERROR;
  stream=fopen(path,"rb");
  if(!stream) return errno==ENOENT?NEC4_CACHE_NOT_FOUND:NEC4_CACHE_IO_ERROR;
  if(!read_exact(stream,magic,sizeof(magic))||
     !read_exact(stream,&format,sizeof(format))||
     !read_exact(stream,&endian,sizeof(endian))||
     !read_exact(stream,&double_size,sizeof(double_size))||
     !read_key(stream,&actual)||
     !read_exact(stream,&stored_checksum,sizeof(stored_checksum)))
  { fclose(stream);return NEC4_CACHE_FORMAT_ERROR; }
  if(memcmp(magic,CACHE_MAGIC,sizeof(magic))||
     format!=NEC4_GROUND_CACHE_FORMAT_VERSION||
     endian!=CACHE_ENDIAN_MARKER||double_size!=sizeof(double))
  { fclose(stream);return NEC4_CACHE_FORMAT_ERROR; }
  if(!keys_equal(&actual,expected))
  { fclose(stream);return NEC4_CACHE_KEY_MISMATCH; }
  if(!valid_key(&actual))
  { fclose(stream);return NEC4_CACHE_FORMAT_ERROR; }
  checksum=checksum_key(FNV_OFFSET,&actual);
  for(index=0;index<(size_t)actual.value_count;++index)
  {
    double pair[2];
    if(!read_exact(stream,pair,sizeof(pair)))
    { fclose(stream);return NEC4_CACHE_FORMAT_ERROR; }
    values[index]=pair[0]+I*pair[1];
    checksum=checksum_bytes(checksum,pair,sizeof(pair));
  }
  trailing=fgetc(stream);
  fclose(stream);
  if(trailing!=EOF) return NEC4_CACHE_FORMAT_ERROR;
  return checksum==stored_checksum?NEC4_CACHE_OK:NEC4_CACHE_CHECKSUM_ERROR;
}

const char *nec4_ground_cache_status_name(nec4_ground_cache_status status)
{
  switch(status)
  {
    case NEC4_CACHE_OK:return "ok";
    case NEC4_CACHE_NOT_FOUND:return "not-found";
    case NEC4_CACHE_KEY_MISMATCH:return "key-mismatch";
    case NEC4_CACHE_FORMAT_ERROR:return "format-error";
    case NEC4_CACHE_CHECKSUM_ERROR:return "checksum-error";
    case NEC4_CACHE_IO_ERROR:return "io-error";
    default:return "invalid";
  }
}
