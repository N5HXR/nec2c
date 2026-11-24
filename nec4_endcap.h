#ifndef NEC2C_NEC4_ENDCAP_H
#define NEC2C_NEC4_ENDCAP_H 1

#include <complex.h>

typedef enum
{
  NEC4_ENDCAP_AXIS,
  NEC4_ENDCAP_QUADRATURE,
  NEC4_ENDCAP_FAR
} nec4_endcap_branch;

typedef struct
{
  double complex radial;
  double complex axial;
  nec4_endcap_branch branch;
} nec4_endcap_field;

/*
 * Return D = -grad integral_A exp(-j k R)/R dA for a flat disk of radius a.
 * The disk is centered at the origin in the z=0 plane.  The observation
 * point is (rho,0,z).  Multiplication by charge density/(4*pi*epsilon)
 * gives the electric field of the disk (NEC4 theory (4-37), Appendix E).
 */
int nec4_endcap_disk_field(double radius, double rho, double z,
                           double complex wave_number,
                           nec4_endcap_field *field);

const char *nec4_endcap_branch_name(nec4_endcap_branch branch);

#endif
