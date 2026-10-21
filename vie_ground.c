#include "vie_ground.h"

#include <math.h>
#include <string.h>

static int above_ground(const vie_cuboid *cell, double minimum_gap)
{
  double gap;
  int axis;
  if(cell == NULL)
    return 0;
  for(axis = 0; axis < 3; ++axis)
    if(!isfinite(cell->center[axis]) || !isfinite(cell->size[axis]) ||
       !(cell->size[axis] > 0.0))
      return 0;
  gap = cell->center[2] - 0.5*cell->size[2];
  return isfinite(gap) && gap >= 0.0 && gap >= minimum_gap;
}

int vie_cuboid_pec_n(const vie_cuboid *testing, const vie_cuboid *source,
                     double k0, double minimum_gap,
                     const vie_kernel_options *options,
                     complex double result[3][3], vie_pec_report *report)
{
  vie_cuboid image;
  complex double direct[3][3], reflected[3][3];
  vie_pec_report local_report;
  int row, column;
  if(result == NULL || !isfinite(minimum_gap) || minimum_gap < 0.0 ||
     !above_ground(testing, minimum_gap) || !above_ground(source, minimum_gap))
    return 0;
  image = *source;
  image.center[2] = -image.center[2];
  memset(&local_report, 0, sizeof(local_report));
  if(!vie_cuboid_n(testing, source, k0, options, direct, &local_report.direct) ||
     !vie_cuboid_n(testing, &image, k0, options, reflected,
                   &local_report.reflected))
    return 0;
  /* Reflection acts on the SOURCE component (columns), not the test rows. */
  for(row = 0; row < 3; ++row)
    for(column = 0; column < 3; ++column)
    {
      direct[row][column] += (column == 2 ? 1.0 : -1.0)*reflected[row][column];
      if(!isfinite(creal(direct[row][column])) ||
         !isfinite(cimag(direct[row][column])))
        return 0;
    }
  memcpy(result, direct, sizeof(direct));
  if(report != NULL)
    *report = local_report;
  return 1;
}
