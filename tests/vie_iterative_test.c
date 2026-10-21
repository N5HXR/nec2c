#include "vie_iterative.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static int failures;
#define CHECK(condition, message) do { if(!(condition)) { \
  fprintf(stderr, "FAIL: %s (line %d)\n", message, __LINE__); ++failures; } } while(0)

static int identity(void *context, const complex double *x, complex double *y)
{
  size_t i, n = *(size_t *)context;
  for(i = 0; i < n; ++i) y[i] = x[i];
  return 1;
}

static int zero(void *context, const complex double *x, complex double *y)
{
  size_t i, n = *(size_t *)context;
  (void)x;
  for(i = 0; i < n; ++i) y[i] = 0.0;
  return 1;
}

static int fail_callback(void *context, const complex double *x, complex double *y)
{
  (void)context; (void)x; (void)y;
  return 0;
}

static int nan_callback(void *context, const complex double *x, complex double *y)
{
  identity(context, x, y); y[0] = NAN;
  return 1;
}

/* A callback without a stored matrix: asymmetric complex tridiagonal stencil. */
static int stencil(void *context, const complex double *x, complex double *y)
{
  size_t i, n = *(size_t *)context;
  for(i = 0; i < n; ++i)
  {
    y[i] = (3.0+0.2*I)*x[i];
    if(i > 0) y[i] += (0.3-0.1*I)*x[i-1];
    if(i+1 < n) y[i] += (-0.5+0.3*I)*x[i+1];
  }
  return 1;
}

static void test_matrix_free(void)
{
  size_t n = 40, i;
  complex double expected[40], b[40], x[40] = {0};
  vie_linear_operator op = {40, stencil, &n, NULL, NULL};
  vie_gmres_options options;
  vie_gmres_report report;
  for(i = 0; i < n; ++i) expected[i] = sin((double)i)+I*cos(0.4*(double)i);
  stencil(&n, expected, b);
  vie_gmres_default_options(&options); options.restart = 3;
  CHECK(vie_gmres(&op, b, x, &options, &report), "matrix-free stencil converges");
  CHECK(report.restarts > 0, "restart path exercised");
  CHECK(report.relative_residual <= options.relative_tolerance, "true residual meets target");
  for(i = 0; i < n; ++i) CHECK(cabs(x[i]-expected[i]) < 1e-10, "matrix-free solution accuracy");
  CHECK(vie_gmres(&op, b, x, &options, &report) && report.iterations == 0,
        "already converged initial guess");
  memset(x, 0, sizeof(x)); options.maximum_iterations = 1;
  CHECK(!vie_gmres(&op, b, x, &options, &report) &&
        report.status == VIE_ITERATIVE_ITERATION_LIMIT && report.iterations == 1 &&
        report.relative_residual > options.relative_tolerance,
        "iteration limit is explicit and retains measured residual");
}

static void test_dense_complex(void)
{
  enum { N = 9 };
  complex double a[N*N] = {0}, expected[N], b[N] = {0}, x[N] = {0};
  vie_gmres_report report;
  size_t i, j;
  for(i = 0; i < N; ++i)
  {
    a[i*N+i] = (double)(i+1)*(2.0+I);
    if(i+1 < N) a[i*N+i+1] = 0.6-0.8*I;
    if(i > 0) a[i*N+i-1] = -0.4+0.1*I;
    expected[i] = (double)(i+1)-I*(double)(i%3);
  }
  for(i = 0; i < N; ++i)
    for(j = 0; j < N; ++j) b[i] += a[i*N+j]*expected[j];
  CHECK(vie_dense_gmres(a, N, b, x, NULL, &report), "complex Jacobi-preconditioned solve");
  for(i = 0; i < N; ++i) CHECK(cabs(x[i]-expected[i]) < 1e-10, "dense complex solution accuracy");
  CHECK(report.relative_residual <= 1e-12, "dense true residual");
}

static void test_edge_cases(void)
{
  size_t n = 3;
  complex double b[3] = {1+I, 2-I, -3+2*I}, x[3] = {0};
  vie_linear_operator op = {3, identity, &n, NULL, NULL};
  vie_gmres_options options;
  vie_gmres_report report;
  vie_gmres_default_options(&options);
  CHECK(vie_gmres(&op, b, x, NULL, &report) && report.iterations == 1,
        "identity happy breakdown succeeds");
  memset(b, 0, sizeof(b)); memset(x, 0, sizeof(x));
  CHECK(vie_gmres(&op, b, x, NULL, &report) && report.iterations == 0 &&
        report.relative_residual == 0.0, "zero RHS exact solution");
  x[0] = 1e-18; options.absolute_tolerance = 1e-15;
  CHECK(vie_gmres(&op, b, x, &options, &report) && report.iterations == 0,
        "absolute tolerance on zero RHS");
  b[0] = 1.0; memset(x, 0, sizeof(x)); op.apply = zero;
  CHECK(!vie_gmres(&op, b, x, NULL, &report) && report.status == VIE_ITERATIVE_BREAKDOWN,
        "inconsistent singular operator cannot claim convergence");
  op.apply = fail_callback;
  CHECK(!vie_gmres(&op, b, x, NULL, &report) && report.status == VIE_ITERATIVE_CALLBACK_FAILED,
        "operator failure propagates");
  op.apply = identity; op.precondition = fail_callback;
  CHECK(!vie_gmres(&op, b, x, NULL, &report) && report.status == VIE_ITERATIVE_CALLBACK_FAILED,
        "preconditioner failure propagates");
  op.precondition = nan_callback; op.precondition_context = &n;
  CHECK(!vie_gmres(&op, b, x, NULL, &report) && report.status == VIE_ITERATIVE_NONFINITE,
        "nonfinite preconditioner rejected");
  op.precondition = NULL; op.apply = nan_callback;
  CHECK(!vie_gmres(&op, b, x, NULL, &report) && report.status == VIE_ITERATIVE_NONFINITE,
        "nonfinite operator rejected");
  op.apply = identity; options.restart = 0;
  CHECK(!vie_gmres(&op, b, x, &options, &report) && report.status == VIE_ITERATIVE_INVALID_INPUT,
        "invalid restart rejected");
  CHECK(!vie_gmres(&op, b, b, NULL, &report) && report.status == VIE_ITERATIVE_INVALID_INPUT,
        "RHS alias rejected");
  op.order = SIZE_MAX;
  CHECK(!vie_gmres(&op, b, x, NULL, &report) && report.status == VIE_ITERATIVE_ALLOCATION_FAILED,
        "workspace size overflow rejected before reading vectors");
  CHECK(!vie_dense_gmres(NULL, 0, b, x, NULL, &report), "invalid dense operator rejected");
}

int main(void)
{
  test_matrix_free(); test_dense_complex(); test_edge_cases();
  if(failures) return EXIT_FAILURE;
  puts("M9 iterative operator, complex GMRES, and failure-path validation passed");
  return 0;
}
