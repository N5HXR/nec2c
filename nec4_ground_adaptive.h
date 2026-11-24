#ifndef NEC2C_NEC4_GROUND_ADAPTIVE_H
#define NEC2C_NEC4_GROUND_ADAPTIVE_H 1

#include "nec4_ground_direct.h"

#include <complex.h>
#include <stddef.h>

typedef int (*nec4_ground_adaptive_sampler)(void *context,double rho_m,
    double observer_z_m,double source_depth_m,
    double complex values[NEC4_GROUND_COMPONENT_COUNT]);

typedef struct {
  double lower[3];
  double upper[3];
} nec4_ground_adaptive_bounds;

typedef struct {
  size_t cells,leaves,accepted_leaves,unresolved_leaves;
  size_t fitted_samples,withheld_samples;
  unsigned int deepest_level;
  double worst_withheld_scaled_error;
  double worst_encountered_scaled_error;
} nec4_ground_adaptive_stats;

typedef struct nec4_ground_adaptive_node nec4_ground_adaptive_node;

typedef struct {
  nec4_ground_adaptive_node *root;
  nec4_ground_adaptive_bounds bounds;
  nec4_ground_adaptive_stats stats;
  double tolerance;
  unsigned int maximum_depth;
  size_t maximum_cells;
  int workers;
  int anisotropic;
} nec4_ground_adaptive;

void nec4_ground_adaptive_init(nec4_ground_adaptive *table);
void nec4_ground_adaptive_destroy(nec4_ground_adaptive *table);

/* Build deterministic tensor-cubic cells.  Fit nodes lie at 0, 1/3, 2/3,
   and 1 on each active axis.  Error is measured only at withheld 1/6, 1/2,
   and 5/6 tensor witnesses.  A zero-width axis is retained exactly and does
   not create duplicate samples or children. */
int nec4_ground_adaptive_build(nec4_ground_adaptive *table,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,int workers,
    nec4_ground_adaptive_sampler sampler,void *sampler_context);

/* Accuracy-first k-d refinement.  Each rejected cell is bisected only along
   the active axis with the largest center-line withheld error.  This avoids
   spending 4x/8x cells when phase variation is strongly directional.  The
   versioned cache records every split axis, so reconstructed trees retain
   the exact deterministic topology. */
int nec4_ground_adaptive_build_anisotropic(nec4_ground_adaptive *table,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,int workers,
    nec4_ground_adaptive_sampler sampler,void *sampler_context);

/* Returns zero when the coordinate is outside the table or its leaf did not
   meet the requested tolerance at maximum_depth.  The caller can then use the
   direct-contour fallback without silently accepting an unresolved cell. */
int nec4_ground_adaptive_evaluate(const nec4_ground_adaptive *table,
    double rho_m,double observer_z_m,double source_depth_m,
    double complex values[NEC4_GROUND_COMPONENT_COUNT]);

/* Versioned, checksummed, atomic persistent cache.  identity is an exact
   caller-defined key containing material, frequency, representation, and
   direct-sampler settings. */
int nec4_ground_adaptive_store(const nec4_ground_adaptive *table,
    const char *path,const void *identity,size_t identity_size);
int nec4_ground_adaptive_load(nec4_ground_adaptive *table,const char *path,
    const void *expected_identity,size_t identity_size);

#endif
