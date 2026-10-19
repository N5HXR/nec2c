#ifndef VIE_KERNEL_H
#define VIE_KERNEL_H 1

#include <complex.h>

typedef struct
{
  double center[3];
  double size[3];
} vie_cuboid;

typedef enum
{
  VIE_RELATION_SELF,
  VIE_RELATION_FACE,
  VIE_RELATION_EDGE,
  VIE_RELATION_CORNER,
  VIE_RELATION_NEAR,
  VIE_RELATION_FAR
} vie_relation;

typedef struct
{
  double relative_tolerance;
  double absolute_tolerance;
  int maximum_order;
  int far_order;
} vie_kernel_options;

typedef struct
{
  vie_relation relation;
  int quadrature_order;
  double estimated_error;
  double estimated_absolute_error;
  int worst_row;
  int worst_column;
} vie_kernel_report;

/*
 * Return the volume-averaged Galerkin matrix of
 *
 *   N u = curl curl integral_V g(r-r') u(r') dV',
 *   g(R) = exp(-i k0 R) / (4 pi R),
 *
 * for constant Cartesian source and test functions on axis-aligned cuboids.
 * The first cuboid is the testing cell and sets the averaging volume.  The
 * result is dimensionless.  k0 is in rad/m and must be nonnegative.  M4's
 * certified domain requires k0 times either cuboid's longest edge <= 1.
 */
int vie_cuboid_n(const vie_cuboid *testing, const vie_cuboid *source,
                 double k0, const vie_kernel_options *options,
                 complex double result[3][3], vie_kernel_report *report);

void vie_kernel_default_options(vie_kernel_options *options);
const char *vie_relation_name(vie_relation relation);

#endif
