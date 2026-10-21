#ifndef VIE_FFT_H
#define VIE_FFT_H
#include "vie_iterative.h"
/* Experimental volume-only backend: equal cubes on a common integer lattice.
   Coordinates are nonnegative grid indices, cell-major xyz. No mesh inference.
   Plan owns mutable FFT workspace: apply is not reentrant. */
typedef struct vie_fft_plan vie_fft_plan;
/* Budget covers all plan-owned allocations, not caller/GMRES storage.
   On failure *out is NULL. k0*pitch <= 1 is the certified kernel domain. */
int vie_fft_create(vie_fft_plan **out, const size_t shape[3], size_t cells,
                   const size_t *xyz, double pitch, double k0, complex double epsilon_r,
                   size_t byte_limit);
/* Assign passive scalar permittivity per cell after creation. Geometry/kernel
   and the plan memory budget remain unchanged; updates the local diagonal. */
int vie_fft_set_permittivity(vie_fft_plan *plan, size_t cell, complex double epsilon_r);
void vie_fft_destroy(vie_fft_plan *plan);
size_t vie_fft_bytes(const vie_fft_plan *plan);
void vie_fft_operator(vie_fft_plan *plan, vie_linear_operator *op);
#endif
