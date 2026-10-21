#ifndef VIE_HALFSPACE_H
#define VIE_HALFSPACE_H 1
#include "vie_kernel.h"

typedef struct {
  double relative_tolerance, absolute_tolerance;
  int maximum_order;
} vie_halfspace_options;
typedef struct {
  int order;
  double electric_delta, magnetic_delta;
} vie_halfspace_report;

/* exp(+i omega t), passive nonmagnetic lower half-space, epsilon = er-i loss.
 * These functions return REFLECTED terms only. N/(i omega epsilon0) maps an
 * electric current moment to E; H maps the moment directly to magnetic field.
 * Point outputs have units m^-3 and m^-2. Cuboid outputs integrate source
 * volume and average testing volume (dimensionless N, H in metres).
 * Positive source/test clearance, k0>0, er>=1 and loss>=0 are required.
 * Maximum-order failure is explicit and leaves output matrices unchanged.
 * Success means successive-quadrature agreement, not a rigorous error bound.
 */
void vie_halfspace_defaults(vie_halfspace_options *options);
int vie_halfspace_fresnel(complex double epsilon, double transverse_k_over_k0,
                          complex double *rs, complex double *rp);
int vie_halfspace_point(const double testing[3], const double source[3],
                        double k0, complex double epsilon,
                        const vie_halfspace_options *options,
                        complex double n[3][3], complex double h[3][3],
                        vie_halfspace_report *report);
int vie_halfspace_cuboids(const vie_cuboid *testing, const vie_cuboid *source,
                          double k0, complex double epsilon,
                          const vie_halfspace_options *options,
                          complex double n[3][3], complex double h[3][3],
                          vie_halfspace_report *report);
/* Mixed interactions: point testing integrates the finite source cell;
 * cuboid testing averages the response to a unit point current moment.
 * Their N matrices obey N_point_cell = source_volume * transpose(N_cell_point)
 * with source and testing positions exchanged. H is not transposed this way.
 */
int vie_halfspace_point_cuboid(const double testing[3], const vie_cuboid *source,
                              double k0, complex double epsilon,
                              const vie_halfspace_options *options,
                              complex double n[3][3], complex double h[3][3],
                              vie_halfspace_report *report);
int vie_halfspace_cuboid_point(const vie_cuboid *testing, const double source[3],
                              double k0, complex double epsilon,
                              const vie_halfspace_options *options,
                              complex double n[3][3], complex double h[3][3],
                              vie_halfspace_report *report);
/* Downward source spectrum: integral J(r) exp(+i k u rho.r_xy-i k w z) dv.
 * The callback includes every source so flux retains all interference terms.
 * clearance is a positive LOWER bound on source height. Flux is delivered to
 * the lower half-space (absorption for lossy ground, transmission if lossless).
 * k and source coordinates are SI; outputs are watts. Same convergence caveat
 * as the field integral applies. Failure leaves *power unchanged.
 */
typedef int (*vie_halfspace_spectrum)(void *context,double k,double u,
                                    complex double w,double phi,complex double source[3]);
int vie_halfspace_power(double k,complex double epsilon,double clearance,
                        vie_halfspace_spectrum spectrum,void *context,
                        const vie_halfspace_options *options,double *power,double *delta);
#endif
