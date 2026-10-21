#ifndef DIELECTRIC_FFT_H
#define DIELECTRIC_FFT_H
#include "dielectric.h"
#include "vie_fft.h"
#define DIELECTRIC_FFT_BUDGET ((size_t)512*1024*1024)
typedef struct dielectric_fft_system dielectric_fft_system;
extern vie_gmres_options dielectric_gmres_settings;
/* Equal cubes on a common lattice; cell-local passive scalar materials. extra_bytes reserves caller-owned solver storage.
   Includes GMRES, temporary lattice map and geometry in the 512 MiB cap. */
int dielectric_fft_prepare(dielectric_fft_system **out,
    const dielectric_model *model,const dielectric_mesh_result *mesh,
    size_t conductors,double k0,size_t extra_bytes,dielectric_diagnostic *diagnostic);
/* PEC image embedding; bottom gap must be an integer multiple of half pitch.
   Only physical cells are unknowns. Image sources carry diag(-1,-1,+1) J.
   Uses the free-space self-block preconditioner as an approximation. */
int dielectric_fft_prepare_pec(dielectric_fft_system **out,
    const dielectric_model *model,const dielectric_mesh_result *mesh,
    size_t conductors,double k0,size_t extra_bytes,dielectric_diagnostic *diagnostic);
void dielectric_fft_free(dielectric_fft_system *system);
/* Only conductor/conductor and cross blocks may be assembled through this API. */
complex double *dielectric_fft_entry(dielectric_fft_system *system,size_t row,size_t column);
void dielectric_fft_get_operator(dielectric_fft_system *system,vie_linear_operator *op);
size_t dielectric_fft_storage(const dielectric_fft_system *system);
#endif
