#ifndef NEC2C_NEC4_GROUND_TABLE_H
#define NEC2C_NEC4_GROUND_TABLE_H 1

#include "nec4_ground_cache.h"
#include "nec4_ground_direct.h"

#include <complex.h>
#include <stddef.h>

typedef struct {
  unsigned int kind;
  double frequency_hz;
  double complex epsilon_relative;
  double wavelength0_m;
  double lg_m;
  double lmax_m;
  double rho_start_m,rho_step_m;
  double observer_start_m,observer_step_m;
  double source_depth_start_m,source_depth_step_m;
  size_t rho_count,observer_count,source_depth_count,value_count;
  double complex *values;
} nec4_ground_table;

typedef int (*nec4_ground_table_sampler)(void *context,double rho_m,
    double observer_z_m,double source_depth_m,
    double complex values[NEC4_GROUND_COMPONENT_COUNT]);

typedef struct {
  unsigned int panel_subdivisions;
  unsigned int quadrature_order;
  double cutoff_decay;
} nec4_ground_table_direct_options;

/* Allocate one of the three Section 5.1.3 direct tables using the qualified
   M9 grid. Values are transformed/stored fields supplied by the callback. */
int nec4_ground_table_init(nec4_ground_table *table,unsigned int kind,
    double frequency_hz,double complex epsilon_relative,double wavelength0_m,
    double lg_m,double lmax_m);
void nec4_ground_table_destroy(nec4_ground_table *table);

/* Sampling is deterministic: each point owns one fixed contiguous block and
   no reduction depends on worker count. The sampler and context must be
   reentrant. The singular table-1 origin is stored as its analytic transformed
   zero without invoking the sampler. workers<=0 uses the OpenMP default, or
   one worker without it. */
int nec4_ground_table_fill(nec4_ground_table *table,
    nec4_ground_table_sampler sampler,void *context,int workers);

/* Generate transformed production values from the accepted Section 5.1.3 /
   Appendix-H direct evaluator.  A NULL options pointer selects the qualified
   production defaults. */
int nec4_ground_table_fill_direct(nec4_ground_table *table,
    double complex lower_k,double complex upper_k,
    const nec4_ground_table_direct_options *options,int workers);

/* Interpolate stored fields. Analytic table-1 add-backs and F1/F2/F3 restore
   factors belong to the physical evaluator above this representation layer. */
int nec4_ground_table_interpolate(const nec4_ground_table *table,double rho_m,
    double observer_z_m,double source_depth_m,
    double complex values[NEC4_GROUND_COMPONENT_COUNT]);

/* Interpolate, undo F1/F2/F3, and for table 1 reconstruct the physical field
   by applying the Appendix-H analytic add-backs. */
int nec4_ground_table_evaluate(const nec4_ground_table *table,
    double complex lower_k,double complex upper_k,double rho_m,
    double observer_z_m,double source_depth_m,
    double complex values[NEC4_GROUND_COMPONENT_COUNT]);

void nec4_ground_table_cache_key(const nec4_ground_table *table,
    nec4_ground_cache_key *key);
nec4_ground_cache_status nec4_ground_table_store(
    const nec4_ground_table *table,const char *path);
nec4_ground_cache_status nec4_ground_table_load(
    nec4_ground_table *table,const char *path);

#endif
