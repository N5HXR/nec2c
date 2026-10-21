#ifndef VIE_ITERATIVE_H
#define VIE_ITERATIVE_H

#include <complex.h>
#include <stddef.h>

/* Callbacks write all n output entries, accept distinct input/output arrays,
   and return zero on failure. A NULL preconditioner means identity. */
typedef int (*vie_vector_map)(void *context, const complex double *input,
                              complex double *output);
typedef struct
{
  size_t order;
  vie_vector_map apply;
  void *context;
  vie_vector_map precondition;
  void *precondition_context;
} vie_linear_operator;

typedef struct
{
  size_t restart;
  size_t maximum_iterations;
  double relative_tolerance;
  double absolute_tolerance;
} vie_gmres_options;

typedef enum
{
  VIE_ITERATIVE_CONVERGED,
  VIE_ITERATIVE_ITERATION_LIMIT,
  VIE_ITERATIVE_BREAKDOWN,
  VIE_ITERATIVE_CALLBACK_FAILED,
  VIE_ITERATIVE_INVALID_INPUT,
  VIE_ITERATIVE_ALLOCATION_FAILED,
  VIE_ITERATIVE_NONFINITE
} vie_iterative_status;

typedef struct
{
  vie_iterative_status status;
  size_t iterations;
  size_t operator_applications;
  size_t restarts;
  double residual_norm;
  double relative_residual;
} vie_gmres_report;

void vie_gmres_default_options(vie_gmres_options *options);
const char *vie_iterative_status_name(vie_iterative_status status);
/* Right-preconditioned restarted complex GMRES. x is an initial guess and is
   updated even on iteration-limit return. Success requires a recomputed
   unpreconditioned ||b-Ax|| <= max(atol, rtol ||b||), never just its estimate.
   RHS and x must not alias. Storage is O(n*restart + restart^2). */
int vie_gmres(const vie_linear_operator *operator, const complex double *rhs,
               complex double *x, const vie_gmres_options *options,
               vie_gmres_report *report);
/* Reference adapter only: retains the dense matrix and uses diagonal Jacobi.
   This is for validation of iteration, not an accelerated production backend. */
int vie_dense_gmres(const complex double *matrix, size_t order,
                     const complex double *rhs, complex double *x,
                     const vie_gmres_options *options, vie_gmres_report *report);

#endif
