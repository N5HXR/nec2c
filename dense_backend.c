#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dense_backend.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define LAPACK_AUTO_MIN_ORDER 128U

static dense_backend_kind requested_backend = DENSE_BACKEND_AUTO;

#ifdef HAVE_LAPACK
#include <lapacke.h>
#endif

int dense_backend_set(const char *name)
{
  if(name == NULL || strcmp(name, "auto") == 0)
    requested_backend = DENSE_BACKEND_AUTO;
  else if(strcmp(name, "builtin") == 0)
    requested_backend = DENSE_BACKEND_BUILTIN;
  else if(strcmp(name, "lapack") == 0)
  {
    if(!dense_backend_lapack_available()) return 0;
    requested_backend = DENSE_BACKEND_LAPACK;
  }
  else return 0;
  return 1;
}

int dense_backend_lapack_available(void)
{
#ifdef HAVE_LAPACK
  return 1;
#else
  return 0;
#endif
}

const char *dense_backend_requested_name(void)
{
  switch(requested_backend)
  {
    case DENSE_BACKEND_BUILTIN: return "builtin";
    case DENSE_BACKEND_LAPACK: return "lapack";
    default: return "auto";
  }
}

static dense_backend_kind backend_for_order(size_t order)
{
  if(requested_backend == DENSE_BACKEND_BUILTIN) return DENSE_BACKEND_BUILTIN;
  if(requested_backend == DENSE_BACKEND_LAPACK) return DENSE_BACKEND_LAPACK;
  if(dense_backend_lapack_available() && order >= LAPACK_AUTO_MIN_ORDER)
    return DENSE_BACKEND_LAPACK;
  return DENSE_BACKEND_BUILTIN;
}

const char *dense_backend_name_for_order(size_t order)
{
  return backend_for_order(order) == DENSE_BACKEND_LAPACK ? "lapack" : "builtin";
}

const char *dense_backend_solver_name_for_order(size_t order)
{
  return backend_for_order(order) == DENSE_BACKEND_LAPACK
      ? "dense-lapack" : "dense-builtin";
}

const char *dense_backend_classic_solver_name_for_order(size_t order)
{
  return backend_for_order(order) == DENSE_BACKEND_LAPACK
      ? "legacy-nec-direct-lapack" : "legacy-nec-direct-builtin";
}

size_t dense_backend_lapack_min_order(void) { return LAPACK_AUTO_MIN_ORDER; }

#ifdef HAVE_LAPACK
static int lapack_dimensions(size_t order, size_t leading_dimension,
                             int *n, int *lda)
{
  if(order == 0 || order > (size_t)INT_MAX ||
     leading_dimension < order || leading_dimension > (size_t)INT_MAX)
    return 0;
  *n = (int)order;
  *lda = (int)leading_dimension;
  return 1;
}
#endif

int dense_backend_factor_classic(size_t order, complex double *matrix,
                                 size_t leading_dimension, int *pivots)
{
#ifdef HAVE_LAPACK
  int n, lda;
  size_t index;
  lapack_int info;
  lapack_int *lapack_pivots;
  if(backend_for_order(order) != DENSE_BACKEND_LAPACK) return -1;
  if(!matrix || !pivots || !lapack_dimensions(order, leading_dimension, &n, &lda))
    return 0;
  lapack_pivots = malloc(order * sizeof(*lapack_pivots));
  if(!lapack_pivots) return 0;
  info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n,
      (lapack_complex_double *)matrix, lda, lapack_pivots);
  if(info == 0)
    for(index = 0; index < order; ++index) pivots[index] = (int)lapack_pivots[index];
  free(lapack_pivots);
  return info == 0;
#else
  (void)order; (void)matrix; (void)leading_dimension; (void)pivots;
  return -1;
#endif
}

int dense_backend_solve_classic(size_t order, const complex double *matrix,
                                size_t leading_dimension, const int *pivots,
                                complex double *rhs)
{
  return dense_backend_solve_classic_many(order, matrix, leading_dimension,
                                           pivots, rhs, order, 1);
}

int dense_backend_solve_classic_many(size_t order,
                                     const complex double *matrix,
                                     size_t leading_dimension,
                                     const int *pivots,
                                     complex double *rhs,
                                     size_t rhs_stride, size_t rhs_count)
{
#ifdef HAVE_LAPACK
  int n, lda, ldb, nrhs;
  size_t index;
  lapack_int info;
  lapack_int *lapack_pivots;
  if(backend_for_order(order) != DENSE_BACKEND_LAPACK) return -1;
  if(rhs_count == 0) return 1;
  if(!matrix || !pivots || !rhs ||
     !lapack_dimensions(order, leading_dimension, &n, &lda) ||
     rhs_count > (size_t)INT_MAX ||
     rhs_stride < order || rhs_stride > (size_t)INT_MAX) return 0;
  ldb = (int)rhs_stride;
  nrhs = (int)rhs_count;
  lapack_pivots = malloc(order * sizeof(*lapack_pivots));
  if(!lapack_pivots) return 0;
  for(index = 0; index < order; ++index) lapack_pivots[index] = (lapack_int)pivots[index];
  info = LAPACKE_zgetrs(LAPACK_COL_MAJOR, 'T', n, nrhs,
      (const lapack_complex_double *)matrix, lda, lapack_pivots,
      (lapack_complex_double *)rhs, ldb);
  free(lapack_pivots);
  return info == 0;
#else
  (void)order; (void)matrix; (void)leading_dimension; (void)pivots;
  (void)rhs; (void)rhs_stride; (void)rhs_count;
  return -1;
#endif
}

static int builtin_row_major_solve(complex double *a, complex double *b, size_t n)
{
  size_t k, i, j;
  for(k = 0; k < n; ++k)
  {
    size_t pivot = k;
    double largest = cabs(a[k*n+k]);
    for(i = k + 1; i < n; ++i)
    {
      double candidate = cabs(a[i*n+k]);
      if(candidate > largest) { largest = candidate; pivot = i; }
    }
    if(!(largest > 64.0 * DBL_EPSILON)) return 0;
    if(pivot != k)
    {
      for(j = k; j < n; ++j)
      {
        complex double temporary = a[k*n+j];
        a[k*n+j] = a[pivot*n+j];
        a[pivot*n+j] = temporary;
      }
      { complex double temporary = b[k]; b[k] = b[pivot]; b[pivot] = temporary; }
    }
    for(i = k + 1; i < n; ++i)
    {
      complex double multiplier = a[i*n+k] / a[k*n+k];
      a[i*n+k] = multiplier;
      for(j = k + 1; j < n; ++j) a[i*n+j] -= multiplier * a[k*n+j];
      b[i] -= multiplier * b[k];
    }
  }
  for(i = n; i-- > 0; )
  {
    complex double sum = b[i];
    for(j = i + 1; j < n; ++j) sum -= a[i*n+j] * b[j];
    b[i] = sum / a[i*n+i];
  }
  return 1;
}

int dense_backend_solve_row_major(complex double *matrix,
                                  complex double *rhs, size_t order)
{
  if(backend_for_order(order) == DENSE_BACKEND_BUILTIN)
    return builtin_row_major_solve(matrix, rhs, order);
#ifdef HAVE_LAPACK
  {
    int n, lda, nrhs = 1;
    lapack_int info, *pivots;
    if(!lapack_dimensions(order, order, &n, &lda)) return 0;
    pivots = malloc(order * sizeof(*pivots));
    if(!pivots) return 0;
    info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n,
        (lapack_complex_double *)matrix, lda, pivots);
    if(info == 0) info = LAPACKE_zgetrs(LAPACK_COL_MAJOR, 'T', n, nrhs,
        (lapack_complex_double *)matrix, lda, pivots,
        (lapack_complex_double *)rhs, n);
    free(pivots);
    return info == 0;
  }
#else
  return 0;
#endif
}
