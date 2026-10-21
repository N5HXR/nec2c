#ifndef VIE_GROUND_H
#define VIE_GROUND_H 1

#include "vie_kernel.h"

typedef struct
{
  vie_kernel_report direct;
  vie_kernel_report reflected;
} vie_pec_report;

/* Volume-averaged electric-current N operator above the PEC plane z=0.
 * Both complete cuboids must lie in z>=0 with clearance at least minimum_gap
 * (metres). minimum_gap=0 permits face contact with PEC: the reflected
 * interaction then uses the existing touching-cuboid kernel. Buried volume
 * is rejected. The image current is diag(-1,-1,+1) J. The deck-level contact
 * policy is validated separately by the dielectric integration layer.
 * No real-ground approximation is made by this interface. On failure the
 * output matrix is unchanged. Reports retain the underlying kernel meaning;
 * in particular far-order evaluation is not an adaptive error certificate.
 */
int vie_cuboid_pec_n(const vie_cuboid *testing, const vie_cuboid *source,
                     double k0, double minimum_gap,
                     const vie_kernel_options *options,
                     complex double result[3][3], vie_pec_report *report);

#endif
