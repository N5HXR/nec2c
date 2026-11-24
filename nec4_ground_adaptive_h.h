#ifndef NEC2C_NEC4_GROUND_ADAPTIVE_H_FIELD_H
#define NEC2C_NEC4_GROUND_ADAPTIVE_H_FIELD_H 1

#include "nec4_ground_adaptive.h"
#include "nec4_ground_table.h"

#include <stdint.h>

typedef struct {
  nec4_ground_adaptive table;
  double frequency_hz;
  double complex lower_wave_number_per_m;
  double complex upper_wave_number_per_m;
  nec4_ground_table_direct_options direct_options;
} nec4_ground_adaptive_h;

void nec4_ground_adaptive_h_init(nec4_ground_adaptive_h *table);
void nec4_ground_adaptive_h_destroy(nec4_ground_adaptive_h *table);

/* Build an independently sampled adaptive H representation.  The same F2
   phase and inverse-square storage normalization used by adaptive E is used
   here, but every sample comes from the direct magnetic spectral kernels. */
int nec4_ground_adaptive_h_build_anisotropic(nec4_ground_adaptive_h *table,
    double frequency_hz,double complex lower_wave_number_per_m,
    double complex upper_wave_number_per_m,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,int workers,
    const nec4_ground_table_direct_options *direct_options);

int nec4_ground_adaptive_h_evaluate(const nec4_ground_adaptive_h *table,
    double rho_m,double observer_z_m,double source_depth_m,
    double complex values[NEC4_GROUND_H_COMPONENT_COUNT]);

int nec4_ground_adaptive_h_store(const nec4_ground_adaptive_h *table,
    const char *path);
int nec4_ground_adaptive_h_load(nec4_ground_adaptive_h *table,const char *path,
    double frequency_hz,double complex lower_wave_number_per_m,
    double complex upper_wave_number_per_m,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,
    const nec4_ground_table_direct_options *direct_options);

uint64_t nec4_ground_adaptive_h_cache_key(double frequency_hz,
    double complex lower_wave_number_per_m,
    double complex upper_wave_number_per_m,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,
    const nec4_ground_table_direct_options *direct_options);

#endif
