#include "vie_iterative.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static double norm(size_t n, const complex double *x)
{
  size_t i;
  double result = 0.0;
  for(i = 0; i < n; ++i)
  {
    if(!isfinite(creal(x[i])) || !isfinite(cimag(x[i]))) return NAN;
    result = hypot(result, cabs(x[i]));
  }
  return result;
}

static void *allocate(size_t n, size_t m, size_t bytes)
{
  if(m == 0 || n > SIZE_MAX/m || n*m > SIZE_MAX/bytes) return NULL;
  return calloc(n*m, bytes);
}

void vie_gmres_default_options(vie_gmres_options *o)
{
  o->restart = 30;
  o->maximum_iterations = 300;
  o->relative_tolerance = 1e-12;
  o->absolute_tolerance = 0.0;
}

const char *vie_iterative_status_name(vie_iterative_status s)
{
  switch(s)
  {
    case VIE_ITERATIVE_CONVERGED: return "converged";
    case VIE_ITERATIVE_ITERATION_LIMIT: return "iteration limit";
    case VIE_ITERATIVE_BREAKDOWN: return "Krylov breakdown";
    case VIE_ITERATIVE_CALLBACK_FAILED: return "operator/preconditioner failure";
    case VIE_ITERATIVE_INVALID_INPUT: return "invalid input";
    case VIE_ITERATIVE_ALLOCATION_FAILED: return "workspace allocation/size failure";
    case VIE_ITERATIVE_NONFINITE: return "nonfinite arithmetic";
  }
  return "unknown status";
}

static int residual(const vie_linear_operator *op, const complex double *b,
                    const complex double *x, complex double *r,
                    vie_gmres_report *report)
{
  size_t i;
  ++report->operator_applications;
  if(!op->apply(op->context, x, r))
  { report->status = VIE_ITERATIVE_CALLBACK_FAILED; return 0; }
  for(i = 0; i < op->order; ++i) r[i] = b[i]-r[i];
  report->residual_norm = norm(op->order, r);
  if(!isfinite(report->residual_norm))
  { report->status = VIE_ITERATIVE_NONFINITE; return 0; }
  return 1;
}

/* Complex unitary Givens rotation [c s; -conj(s) c]. */
static void rotation(complex double a, complex double b,
                      double *c, complex double *s)
{
  double aa = cabs(a), bb = cabs(b);
  if(bb == 0.0) { *c = 1.0; *s = 0.0; }
  else if(aa == 0.0) { *c = 0.0; *s = conj(b)/bb; }
  else
  {
    double length = hypot(aa, bb);
    *c = aa/length;
    *s = (a/aa)*(conj(b)/length);
  }
}

int vie_gmres(const vie_linear_operator *op, const complex double *b,
               complex double *x, const vie_gmres_options *options,
               vie_gmres_report *report)
{
  vie_gmres_report local;
  vie_gmres_options defaults;
  complex double *v = NULL, *z = NULL, *h = NULL, *r = NULL, *w = NULL;
  complex double *g = NULL, *y = NULL, *s = NULL;
  double *c = NULL, bnorm, target;
  size_t n, m, i, j, k;
  int success = 0;
  if(report == NULL) report = &local;
  memset(report, 0, sizeof(*report));
  report->status = VIE_ITERATIVE_INVALID_INPUT;
  report->residual_norm = report->relative_residual = INFINITY;
  vie_gmres_default_options(&defaults);
  if(options == NULL) options = &defaults;
  if(op == NULL || op->apply == NULL || op->order == 0 || b == NULL ||
     x == NULL || b == x || options->restart == 0 ||
     options->maximum_iterations == 0 ||
     !isfinite(options->relative_tolerance) || options->relative_tolerance < 0.0 ||
     !isfinite(options->absolute_tolerance) || options->absolute_tolerance < 0.0 ||
     (options->relative_tolerance == 0.0 && options->absolute_tolerance == 0.0))
    return 0;
  n = op->order;
  m = options->restart;
  if(m > n) m = n;
  if(m > options->maximum_iterations) m = options->maximum_iterations;
  if(m == SIZE_MAX) { report->status = VIE_ITERATIVE_ALLOCATION_FAILED; return 0; }
  v = allocate(n, m+1, sizeof(*v)); z = allocate(n, m, sizeof(*z));
  h = allocate(m+1, m, sizeof(*h)); r = allocate(n, 1, sizeof(*r));
  w = allocate(n, 1, sizeof(*w)); g = allocate(m+1, 1, sizeof(*g));
  y = allocate(m, 1, sizeof(*y)); s = allocate(m, 1, sizeof(*s));
  c = allocate(m, 1, sizeof(*c));
  if(!v || !z || !h || !r || !w || !g || !y || !s || !c)
  { report->status = VIE_ITERATIVE_ALLOCATION_FAILED; goto done; }
  bnorm = norm(n, b);
  target = fmax(options->absolute_tolerance, options->relative_tolerance*bnorm);
  if(!isfinite(bnorm) || !isfinite(norm(n, x)) || !isfinite(target))
  { report->status = VIE_ITERATIVE_NONFINITE; goto done; }
  if(!residual(op, b, x, r, report)) goto done;
  for(;;)
  {
    size_t used = 0;
    int breakdown = 0;
    double beta = report->residual_norm;
    report->relative_residual = bnorm > 0.0 ? beta/bnorm : beta;
    if(beta <= target)
    { report->status = VIE_ITERATIVE_CONVERGED; success = 1; goto done; }
    if(report->iterations >= options->maximum_iterations)
    { report->status = VIE_ITERATIVE_ITERATION_LIMIT; goto done; }
    memset(h, 0, (m+1)*m*sizeof(*h));
    memset(g, 0, (m+1)*sizeof(*g)); g[0] = beta;
    for(i = 0; i < n; ++i) v[i] = r[i]/beta;
    for(j = 0; j < m && report->iterations < options->maximum_iterations; ++j)
    {
      double next, before;
      int pass;
      if(op->precondition)
      {
        if(!op->precondition(op->precondition_context, v+j*n, z+j*n))
        { report->status = VIE_ITERATIVE_CALLBACK_FAILED; goto done; }
      }
      else memcpy(z+j*n, v+j*n, n*sizeof(*z));
      if(!isfinite(norm(n, z+j*n)))
      { report->status = VIE_ITERATIVE_NONFINITE; goto done; }
      ++report->operator_applications;
      if(!op->apply(op->context, z+j*n, w))
      { report->status = VIE_ITERATIVE_CALLBACK_FAILED; goto done; }
      before = norm(n, w);
      if(!isfinite(before)) { report->status = VIE_ITERATIVE_NONFINITE; goto done; }
      /* Two modified Gram-Schmidt passes reduce loss of orthogonality. */
      for(pass = 0; pass < 2; ++pass)
        for(k = 0; k <= j; ++k)
        {
          complex double projection = 0.0;
          for(i = 0; i < n; ++i) projection += conj(v[k*n+i])*w[i];
          h[k*m+j] += projection;
          for(i = 0; i < n; ++i) w[i] -= projection*v[k*n+i];
        }
      next = norm(n, w);
      if(!isfinite(next)) { report->status = VIE_ITERATIVE_NONFINITE; goto done; }
      h[(j+1)*m+j] = next;
      breakdown = next <= 32*DBL_EPSILON*before;
      if(!breakdown)
        for(i = 0; i < n; ++i) v[(j+1)*n+i] = w[i]/next;
      for(k = 0; k < j; ++k)
      {
        complex double first = h[k*m+j], second = h[(k+1)*m+j];
        h[k*m+j] = c[k]*first+s[k]*second;
        h[(k+1)*m+j] = -conj(s[k])*first+c[k]*second;
      }
      rotation(h[j*m+j], h[(j+1)*m+j], &c[j], &s[j]);
      h[j*m+j] = c[j]*h[j*m+j]+s[j]*h[(j+1)*m+j];
      h[(j+1)*m+j] = 0.0;
      g[j+1] = -conj(s[j])*g[j]; g[j] *= c[j];
      ++report->iterations; used = j+1;
      if(breakdown || cabs(g[j+1]) <= target) break;
    }
    for(i = used; i-- > 0; )
    {
      complex double value = g[i];
      for(k = i+1; k < used; ++k) value -= h[i*m+k]*y[k];
      if(cabs(h[i*m+i]) == 0.0)
      { report->status = VIE_ITERATIVE_BREAKDOWN; goto done; }
      y[i] = value/h[i*m+i];
    }
    for(j = 0; j < used; ++j)
      for(i = 0; i < n; ++i) x[i] += z[j*n+i]*y[j];
    if(!residual(op, b, x, r, report)) goto done;
    report->relative_residual = bnorm > 0.0 ? report->residual_norm/bnorm : report->residual_norm;
    if(breakdown && report->residual_norm > target)
    { report->status = VIE_ITERATIVE_BREAKDOWN; goto done; }
    if(report->residual_norm > target && report->iterations < options->maximum_iterations)
      ++report->restarts;
  }
done:
  free(c); free(s); free(y); free(g); free(w); free(r); free(h); free(z); free(v);
  return success;
}

typedef struct { const complex double *matrix; size_t n; } dense_context;

static int dense_apply(void *context, const complex double *x, complex double *y)
{
  const dense_context *a = context;
  size_t i, j;
  for(i = 0; i < a->n; ++i)
  {
    y[i] = 0.0;
    for(j = 0; j < a->n; ++j) y[i] += a->matrix[i*a->n+j]*x[j];
  }
  return 1;
}

static int diagonal_apply(void *context, const complex double *x, complex double *y)
{
  const dense_context *a = context;
  size_t i;
  for(i = 0; i < a->n; ++i)
  {
    complex double diagonal = a->matrix[i*a->n+i];
    if(cabs(diagonal) == 0.0) return 0;
    y[i] = x[i]/diagonal;
  }
  return 1;
}

int vie_dense_gmres(const complex double *matrix, size_t n,
                     const complex double *b, complex double *x,
                     const vie_gmres_options *options, vie_gmres_report *report)
{
  dense_context context = {matrix, n};
  vie_linear_operator op = {n, dense_apply, &context, diagonal_apply, &context};
  if(matrix == NULL || n == 0 || n > SIZE_MAX/n || n*n > SIZE_MAX/sizeof(*matrix))
  {
    if(report)
    {
      memset(report, 0, sizeof(*report)); report->status = VIE_ITERATIVE_INVALID_INPUT;
      report->residual_norm = report->relative_residual = INFINITY;
    }
    return 0;
  }
  return vie_gmres(&op, b, x, options, report);
}
