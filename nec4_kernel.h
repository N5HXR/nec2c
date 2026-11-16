#ifndef NEC2C_NEC4_KERNEL_H
#define NEC2C_NEC4_KERNEL_H 1

#include <complex.h>

typedef enum
{
  NEC4_FIELD_EXACT,
  NEC4_FIELD_NEAR_AXIS,
  NEC4_FIELD_SMALL_KR,
  NEC4_FIELD_SMALL_KR_AXIS,
  NEC4_FIELD_FAR,
  NEC4_FIELD_FAR_SMALL_KR
} nec4_field_branch;

typedef struct
{
  double complex ezs, ers;
  double complex ezc, erc;
  double complex ezk, erk;
  nec4_field_branch branch;
} nec4_segment_fields;

/*
 * Documented NEC4 reduced thin-wire field, theory section 4.1 and Appendix B.
 * The source segment is centered at zero, has full length s, and rho is the
 * reduced-kernel radial distance from its surface filament to the match axis.
 */
int nec4_reduced_segment_field(double s, double z, double rho,
                               double wave_number,
                               nec4_segment_fields *field);

const char *nec4_field_branch_name(nec4_field_branch branch);

#endif
