#ifndef NEC2C_NEC4_GROUND_CACHE_H
#define NEC2C_NEC4_GROUND_CACHE_H 1

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define NEC4_GROUND_CACHE_FORMAT_VERSION 2u
#define NEC4_GROUND_EQUATION_VERSION 0x00040009u

typedef struct {
  uint32_t equation_version;
  uint32_t precision_bits;
  uint32_t table_kind;
  uint32_t component_count;
  double frequency_hz;
  double epsilon_real;
  double epsilon_imag;
  double wavelength0_m;
  double lg_m;
  double lmax_m;
  double rho_start_m;
  double rho_step_m;
  double observer_start_m;
  double observer_step_m;
  double source_depth_start_m;
  double source_depth_step_m;
  uint64_t rho_count;
  uint64_t observer_count;
  uint64_t source_depth_count;
  uint64_t value_count;
} nec4_ground_cache_key;

typedef enum {
  NEC4_CACHE_OK=0,
  NEC4_CACHE_NOT_FOUND=1,
  NEC4_CACHE_KEY_MISMATCH=2,
  NEC4_CACHE_FORMAT_ERROR=3,
  NEC4_CACHE_CHECKSUM_ERROR=4,
  NEC4_CACHE_IO_ERROR=5
} nec4_ground_cache_status;

nec4_ground_cache_status nec4_ground_cache_store(const char *path,
    const nec4_ground_cache_key *key,const double complex *values);
nec4_ground_cache_status nec4_ground_cache_load(const char *path,
    const nec4_ground_cache_key *expected,double complex *values);
const char *nec4_ground_cache_status_name(nec4_ground_cache_status status);
uint64_t nec4_ground_cache_key_hash(const nec4_ground_cache_key *key);

#endif
