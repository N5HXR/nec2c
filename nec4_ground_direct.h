#ifndef NEC2C_NEC4_GROUND_DIRECT_H
#define NEC2C_NEC4_GROUND_DIRECT_H 1

#include <complex.h>

enum { NEC4_GROUND_COMPONENT_COUNT=5 };

typedef enum {
  NEC4_GROUND_RHO_VERTICAL=0,
  NEC4_GROUND_Z_VERTICAL=1,
  NEC4_GROUND_RHO_HORIZONTAL=2,
  NEC4_GROUND_PHI_HORIZONTAL=3,
  NEC4_GROUND_Z_HORIZONTAL=4
} nec4_ground_component;

typedef struct {
  double complex value[NEC4_GROUND_COMPONENT_COUNT];
  unsigned int panels;
  unsigned int quadrature_order;
} nec4_ground_direct_result;

enum { NEC4_GROUND_H_COMPONENT_COUNT=4 };

typedef enum {
  NEC4_GROUND_H_PHI_VERTICAL=0,
  NEC4_GROUND_H_RHO_HORIZONTAL=1,
  NEC4_GROUND_H_PHI_HORIZONTAL=2,
  NEC4_GROUND_H_Z_HORIZONTAL=3
} nec4_ground_h_component;

/* Magnetic-field kernels use the same cylindrical dyadic convention as the
   electric kernels.  For a horizontal x-directed source,

     H_rho = H_RHO_HORIZONTAL sin(phi)
     H_phi = H_PHI_HORIZONTAL cos(phi)
     H_z   = H_Z_HORIZONTAL sin(phi).

   The values are derived at the spectral-integrand level from
   curl(E)=-j omega mu H; they are not numerical derivatives of E. */
typedef struct {
  double complex value[NEC4_GROUND_H_COMPONENT_COUNT];
  unsigned int panels;
  unsigned int quadrature_order;
} nec4_ground_direct_h_result;

typedef struct {
  nec4_ground_direct_h_result field;
  double refinement_error;
  double approach_path_error;
} nec4_ground_h_interface_result;

/* Exact partition of the documented Figure 5-3a contour at P2, on the
   upper-medium branch transition.  The two arrays are numerical contour
   contributions, not a claim to reproduce NEC-4.2 GN3's unpublished
   branch-cut representation.  Their ordered sum is value[]. */
typedef struct {
  double complex value[NEC4_GROUND_COMPONENT_COUNT];
  double complex lower_spectral[NEC4_GROUND_COMPONENT_COUNT];
  double complex upper_spectral[NEC4_GROUND_COMPONENT_COUNT];
  unsigned int panels;
  unsigned int quadrature_order;
} nec4_ground_direct_split_result;

typedef struct {
  double complex total[NEC4_GROUND_COMPONENT_COUNT];
  double complex finite_remainder[NEC4_GROUND_COMPONENT_COUNT];
  unsigned int panels;
  unsigned int quadrature_order;
} nec4_ground_table1_result;

typedef struct {
  nec4_ground_table1_result field;
  double refinement_error;
  double approach_path_error;
} nec4_ground_interface_result;

/* Direct Section 5.1 cross-interface Bessel-contour evaluator.
   lower_k and upper_k use the exp(+j omega t) passive-medium convention.
   source_depth and observer_z are nonnegative distances from z=0. */
int nec4_ground_direct_cross_bessel(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int quadrature_order,
    double cutoff_decay,nec4_ground_direct_result *result);

int nec4_ground_direct_cross_h_bessel(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int quadrature_order,
    double cutoff_decay,nec4_ground_direct_h_result *result);

/* Figure 5-3a full-line Hankel-contour evaluator.  Both endpoints cannot be
   exactly on the interface; use the physical-limit entry point below. */
int nec4_ground_direct_cross_hankel(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int panel_subdivisions,
    unsigned int quadrature_order,double cutoff_decay,
    nec4_ground_direct_result *result);

int nec4_ground_direct_cross_h_hankel(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int panel_subdivisions,
    unsigned int quadrature_order,double cutoff_decay,
    nec4_ground_direct_h_result *result);

/* Physical z=z'=0 magnetic-field limit for rho>0, independently checked
   along symmetric and source-only approach paths. */
int nec4_ground_direct_cross_h_interface(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    unsigned int panel_subdivisions,unsigned int quadrature_order,
    double cutoff_decay,nec4_ground_h_interface_result *result);

int nec4_ground_direct_cross_hankel_split(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int panel_subdivisions,
    unsigned int quadrature_order,double cutoff_decay,
    nec4_ground_direct_split_result *result);

/* Appendix-H two-stage subtraction used to generate linear table 1. */
int nec4_ground_direct_cross_table1(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int quadrature_order,
    double cutoff_decay,nec4_ground_table1_result *result);

/* Figure 5-3a evaluation of the Appendix-H finite table-1 quantity at a
   positive total endpoint height. */
int nec4_ground_direct_cross_table1_hankel(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,unsigned int panel_subdivisions,
    unsigned int quadrature_order,double cutoff_decay,
    nec4_ground_table1_result *result);

/* Physical z=z'=0 table-1 limit for rho>0.  Richardson extrapolation is
   independently checked along symmetric and source-only approaches. */
int nec4_ground_direct_cross_interface(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    unsigned int panel_subdivisions,unsigned int quadrature_order,
    double cutoff_decay,nec4_ground_interface_result *result);

/* Reconstruct the physical table-1 field from Appendix H's finite E'' value
   by adding the residual 1/R field and V0 field in that order. */
int nec4_ground_direct_table1_restore(double frequency_hz,
    double complex lower_k,double complex upper_k,double rho,
    double observer_z,double source_depth,
    const double complex finite_remainder[NEC4_GROUND_COMPONENT_COUNT],
    double complex total[NEC4_GROUND_COMPONENT_COUNT]);

#endif
