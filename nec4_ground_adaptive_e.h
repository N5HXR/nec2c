#ifndef NEC2C_NEC4_GROUND_ADAPTIVE_E_H
#define NEC2C_NEC4_GROUND_ADAPTIVE_E_H 1

#include "nec4_ground_adaptive.h"
#include "nec4_ground_table.h"

#include <stdint.h>

typedef struct {
  nec4_ground_adaptive table;
  double frequency_hz;
  double complex lower_wave_number_per_m;
  double complex upper_wave_number_per_m;
  nec4_ground_table_direct_options direct_options;
} nec4_ground_adaptive_e;

void nec4_ground_adaptive_e_init(nec4_ground_adaptive_e *table);
void nec4_ground_adaptive_e_destroy(nec4_ground_adaptive_e *table);

/* Build an adaptive total-E representation after removing the accepted F2
   ray phase and inverse-square spreading factor.  Bounds use metres. */
int nec4_ground_adaptive_e_build(nec4_ground_adaptive_e *table,
    double frequency_hz,double complex lower_wave_number_per_m,
    double complex upper_wave_number_per_m,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,int workers,
    const nec4_ground_table_direct_options *direct_options);

int nec4_ground_adaptive_e_build_anisotropic(nec4_ground_adaptive_e *table,
    double frequency_hz,double complex lower_wave_number_per_m,
    double complex upper_wave_number_per_m,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,int workers,
    const nec4_ground_table_direct_options *direct_options);

int nec4_ground_adaptive_e_evaluate(const nec4_ground_adaptive_e *table,
    double rho_m,double observer_z_m,double source_depth_m,
    double complex values[NEC4_GROUND_COMPONENT_COUNT]);

int nec4_ground_adaptive_e_store(const nec4_ground_adaptive_e *table,
    const char *path);
int nec4_ground_adaptive_e_load(nec4_ground_adaptive_e *table,const char *path,
    double frequency_hz,double complex lower_wave_number_per_m,
    double complex upper_wave_number_per_m,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,
    const nec4_ground_table_direct_options *direct_options);

int nec4_ground_adaptive_e_load_anisotropic(nec4_ground_adaptive_e *table,
    const char *path,double frequency_hz,
    double complex lower_wave_number_per_m,
    double complex upper_wave_number_per_m,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,
    const nec4_ground_table_direct_options *direct_options);

/* Stable cache filename key.  Store/load still verify the complete identity
   and payload checksum; this hash is only used to avoid filename collisions. */
uint64_t nec4_ground_adaptive_e_cache_key(double frequency_hz,
    double complex lower_wave_number_per_m,
    double complex upper_wave_number_per_m,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,int anisotropic,
    const nec4_ground_table_direct_options *direct_options);

#endif
