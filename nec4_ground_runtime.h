#ifndef NEC2C_NEC4_GROUND_RUNTIME_H
#define NEC2C_NEC4_GROUND_RUNTIME_H 1

#include "nec4_ground_matrix.h"
#include "nec4_ground_dispatch.h"
#include "nec4_ground_model_grid.h"
#include "nec4_ground_adaptive_e.h"
#include "nec4_ground_adaptive_h.h"

#include <complex.h>
#include <stddef.h>

typedef enum {
  NEC4_GROUND_MODEL_GRID=0,
  NEC4_GROUND_MODEL_DIRECT=1,
  NEC4_GROUND_MODEL_ADAPTIVE=2
} nec4_ground_model;

enum { NEC4_GROUND_ADAPTIVE_REGION_COUNT=3 };

typedef struct {
  const char *cache_directory;
  int workers;
  nec4_ground_model model;
  nec4_ground_table_sampler sampler;
  void *sampler_context;
  const nec4_ground_table_direct_options *direct_options;
} nec4_ground_runtime_options;

typedef struct {
  nec4_ground_table table[3];
  nec4_ground_table model_grid;
  nec4_ground_adaptive_e adaptive_e[NEC4_GROUND_ADAPTIVE_REGION_COUNT];
  unsigned char adaptive_e_ready[NEC4_GROUND_ADAPTIVE_REGION_COUNT];
  nec4_ground_adaptive_h adaptive_h[NEC4_GROUND_ADAPTIVE_REGION_COUNT];
  unsigned char adaptive_h_ready[NEC4_GROUND_ADAPTIVE_REGION_COUNT];
  nec4_ground_matrix_tables matrix_tables;
  nec4_ground_matrix_context matrix_context;
  double wavelength0_m,lg_m,lmax_m;
  nec4_ground_table_direct_options direct_options;
  size_t cache_hits,cache_misses,generated_tables,payload_bytes;
  nec4_ground_model model;
  int workers;
  size_t adaptive_hits,direct_fallbacks,asymptotic_hits,linear_hits;
  size_t h_adaptive_hits,h_direct_fallbacks;
  int grazing_seen;
  int grazing_reported;
  int direct_strip_seen;
  int identical_media;
  int ready;
} nec4_ground_runtime;

void nec4_ground_runtime_init(nec4_ground_runtime *runtime);
void nec4_ground_runtime_destroy(nec4_ground_runtime *runtime);

/* Prepare the selected linear/grid/adaptive table bank for one
   frequency/material pair.  A sampler in options is a test/qualification
   hook; production callers leave it NULL to use the accepted direct
   evaluator. */
int nec4_ground_runtime_prepare(nec4_ground_runtime *runtime,
    double frequency_hz,double complex lower_epsilon_relative,
    double complex lower_wave_number_per_m,
    double complex upper_wave_number_per_m,
    const nec4_ground_runtime_options *options);

const nec4_ground_matrix_context *nec4_ground_runtime_matrix_context(
    const nec4_ground_runtime *runtime);
int nec4_ground_runtime_grazing_seen(const nec4_ground_runtime *runtime);
int nec4_ground_runtime_claim_grazing_warning(nec4_ground_runtime *runtime);
int nec4_ground_runtime_direct_strip_seen(const nec4_ground_runtime *runtime);
const char *nec4_ground_model_name(nec4_ground_model model);

/* Evaluate the accepted production regions.  Adaptive mode uses bounded
   directional tables only in qualified regions and retains the direct contour
   as its correctness-first fallback. */
int nec4_ground_runtime_evaluate(const nec4_ground_runtime *runtime,
    double rho_m,double observer_z_m,double source_depth_m,
    double complex field[NEC4_GROUND_COMPONENT_COUNT],
    nec4_ground_dispatch *route);
int nec4_ground_runtime_evaluate_h(const nec4_ground_runtime *runtime,
    double rho_m,double observer_z_m,double source_depth_m,
    double complex field[NEC4_GROUND_H_COMPONENT_COUNT],
    nec4_ground_dispatch *route);

#endif
