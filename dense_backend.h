#ifndef DENSE_BACKEND_H
#define DENSE_BACKEND_H 1

#include <complex.h>
#include <stddef.h>

typedef enum
{
  DENSE_BACKEND_AUTO = 0,
  DENSE_BACKEND_BUILTIN = 1,
  DENSE_BACKEND_LAPACK = 2
} dense_backend_kind;

int dense_backend_set(const char *name);
int dense_backend_lapack_available(void);
const char *dense_backend_requested_name(void);
const char *dense_backend_name_for_order(size_t order);
size_t dense_backend_lapack_min_order(void);

/* Classic NEC stores the algebraic matrix transposed. LAPACK factors that
   storage directly and solves with TRANS='T', avoiding a physical transpose. */
int dense_backend_factor_classic(size_t order, complex double *matrix,
                                 size_t leading_dimension, int *pivots);
int dense_backend_solve_classic(size_t order, const complex double *matrix,
                                size_t leading_dimension, const int *pivots,
                                complex double *rhs);
int dense_backend_solve_classic_many(size_t order,
                                     const complex double *matrix,
                                     size_t leading_dimension,
                                     const int *pivots,
                                     complex double *rhs,
                                     size_t rhs_stride, size_t rhs_count);

/* Dense VIE matrices are row-major, which is also column-major storage of the
   transpose. The LAPACK path uses the same transpose-solve convention. */
int dense_backend_solve_row_major(complex double *matrix,
                                  complex double *rhs, size_t order);

#endif
