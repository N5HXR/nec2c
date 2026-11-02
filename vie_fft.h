#ifndef VIE_FFT_H
#define VIE_FFT_H
#include "vie_iterative.h"
/* Experimental volume-only backend: equal cubes on a common integer lattice.
   Coordinates are nonnegative grid indices, cell-major xyz. No mesh inference.
   Plan owns mutable FFT workspace: apply is not reentrant. */
typedef struct vie_fft_plan vie_fft_plan;
int vie_fft_backend_set(const char *name);
int vie_fft_fftw_available(void);
const char *vie_fft_backend_requested_name(void);
const char *vie_fft_backend_name(const vie_fft_plan *plan);
/* Backend used by the most recently created plan in this process, or NULL
   before any plan has been created.  This lets forked sweep workers report
   their actual auto-selection to the parent timing document. */
const char *vie_fft_backend_last_used_name(void);
void vie_fft_backend_clear_last_used(void);
/* Budget covers application-owned plan allocations, not caller/GMRES storage.
   Optional external libraries may retain small opaque plan metadata; peak RSS
   is measured separately. On failure *out is NULL. k0*pitch <= 1 is the
   certified kernel domain. */
int vie_fft_create(vie_fft_plan **out, const size_t shape[3], size_t cells,
                   const size_t *xyz, double pitch, double k0, complex double epsilon_r,
                   size_t byte_limit);
/* Assign passive scalar permittivity per cell after creation. Geometry/kernel
   and the plan memory budget remain unchanged; updates the local diagonal. */
int vie_fft_set_permittivity(vie_fft_plan *plan, size_t cell, complex double epsilon_r);
void vie_fft_destroy(vie_fft_plan *plan);
size_t vie_fft_bytes(const vie_fft_plan *plan);
/* Immutable thread budget captured when the plan is created. */
size_t vie_fft_threads(const vie_fft_plan *plan);
int vie_fft_parallel_available(void);
void vie_fft_operator(vie_fft_plan *plan, vie_linear_operator *op);
#endif
