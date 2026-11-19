#ifndef NEC2C_NEC4_GROUND_MODEL_GRID_H
#define NEC2C_NEC4_GROUND_MODEL_GRID_H 1

#include "nec4_ground_table.h"

/* Independently developed, accuracy-qualified replacement for the missing
   tuned equation-(5-10) coefficients.  Values use the accepted F2 transform
   on a level-3 tensor grid; z below 0.16 lambda remains direct-contour. */
int nec4_ground_model_grid_init(nec4_ground_table *grid,double frequency_hz,
    double complex epsilon_relative,double wavelength0_m,double lg_m,
    double lmax_m);
int nec4_ground_model_grid_fill_direct(nec4_ground_table *grid,
    double complex lower_k,double complex upper_k,
    const nec4_ground_table_direct_options *options,int workers);
int nec4_ground_model_grid_evaluate(const nec4_ground_table *grid,
    double complex lower_k,double complex upper_k,double rho_m,
    double observer_z_m,double source_depth_m,
    double complex values[NEC4_GROUND_COMPONENT_COUNT]);

#endif
