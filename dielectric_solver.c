#include "dielectric.h"
#include "vie_kernel.h"
#include "vie_ground.h"
#include "vie_halfspace.h"
#include "shared.h"
#include "vie_iterative.h"
#include "dielectric_fft.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define VIE_EPSILON_0 8.8541878128e-12
#define VIE_MU_0 1.25663706212e-6
#define VIE_C0 299792458.0
#define VIE_PI 3.141592653589793238462643383279502884

static const dielectric_material *material_by_id(const dielectric_model *model,
                                                  int id)
{
  size_t i;
  for(i = 0; i < model->material_count; ++i)
    if(model->materials[i].id == id)
      return &model->materials[i];
  return NULL;
}

static void solver_error(dielectric_diagnostic *d, const dielectric_model *model,
                         const char *record, const char *value,
                         const char *message)
{
  memset(d, 0, sizeof(*d));
  snprintf(d->resource, sizeof(d->resource), "%s",
           model->resource_path != NULL ? model->resource_path : "<dielectric>");
  snprintf(d->record, sizeof(d->record), "%s", record);
  snprintf(d->value, sizeof(d->value), "%s", value);
  snprintf(d->message, sizeof(d->message), "%s", message);
}

static double sinc_real(double x)
{
  if(fabs(x) < 1.0e-5)
  {
    double x2 = x * x;
    return 1.0 - x2 / 6.0 + x2 * x2 / 120.0;
  }
  return sin(x) / x;
}

static void cross3(const complex double a[3], const double b[3],
                   complex double out[3])
{
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

void dielectric_plane_wave_vectors(const dielectric_plane_wave *wave,
                               double direction[3], complex double field[3])
{
  double theta = wave->theta_degrees * VIE_PI / 180.0;
  double phi = wave->phi_degrees * VIE_PI / 180.0;
  double eta = wave->eta_degrees * VIE_PI / 180.0;
  double ct = cos(theta), st = sin(theta);
  double cp = cos(phi), sp = sin(phi);
  double ce = cos(eta), se = sin(eta);
  double p[3], q[3];
  p[0] = ct * cp * ce - sp * se;
  p[1] = ct * sp * ce + cp * se;
  p[2] = -st * ce;
  direction[0] = -st * cp;
  direction[1] = -st * sp;
  direction[2] = -ct;
  q[0] = direction[1] * p[2] - direction[2] * p[1];
  q[1] = direction[2] * p[0] - direction[0] * p[2];
  q[2] = direction[0] * p[1] - direction[1] * p[0];
  field[0] = p[0]; field[1] = p[1]; field[2] = p[2];
  if(wave->type == 2 || wave->type == 3)
  {
    double sign = wave->type == 2 ? -1.0 : 1.0;
    double ratio = wave->axial_ratio == 0.0 ? 1.0 : wave->axial_ratio;
    field[0] += I * sign * ratio * q[0];
    field[1] += I * sign * ratio * q[1];
    field[2] += I * sign * ratio * q[2];
  }
}

static complex double cell_phase_average(const dielectric_cell *cell,
                                         const double direction[3], double k0,
                                         double phase_sign)
{
  double dot = direction[0] * cell->center[0] +
               direction[1] * cell->center[1] +
               direction[2] * cell->center[2];
  double average = 1.0;
  int axis;
  for(axis = 0; axis < 3; ++axis)
    average *= sinc_real(0.5 * k0 * direction[axis] * cell->size[axis]);
  return average * cexp(I * phase_sign * k0 * dot);
}

static complex double background_component(const dielectric_cell *cell,
                       const double direction[3], const complex double field[3],
                       double k0, int component)
{
  complex double value=field[component]*cell_phase_average(cell,direction,k0,-1.0);
  if(gnd.ksymp != 1)
  {
    double reflected[3]={direction[0],direction[1],-direction[2]};
    complex double reflected_field[3];
    dielectric_ground_reflect(direction,field,k0,reflected_field);
    value+=reflected_field[component]*cell_phase_average(cell,reflected,k0,-1.0);
  }
  return value;
}

static int lu_solve(complex double *a, complex double *b, size_t n)
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
    if(!(largest > 64.0 * DBL_EPSILON))
      return 0;
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
      for(j = k + 1; j < n; ++j)
        a[i*n+j] -= multiplier * a[k*n+j];
      b[i] -= multiplier * b[k];
    }
  }
  for(i = n; i-- > 0; )
  {
    complex double sum = b[i];
    for(j = i + 1; j < n; ++j)
      sum -= a[i*n+j] * b[j];
    b[i] = sum / a[i*n+i];
  }
  return 1;
}

void dielectric_radiation_amplitude(const dielectric_mesh_result *mesh,
                                const complex double *current, double k0,
                                const double direction[3],
                                complex double amplitude[3])
{
  complex double integral[3] = {0.0, 0.0, 0.0};
  complex double crossed[3], transverse[3];
  size_t cell;
  int component;
  for(cell = 0; cell < mesh->cell_count; ++cell)
  {
    double volume = mesh->cells[cell].size[0] * mesh->cells[cell].size[1] *
                    mesh->cells[cell].size[2];
    complex double weight = volume * cell_phase_average(&mesh->cells[cell],
                                                         direction, k0, 1.0);
    for(component = 0; component < 3; ++component)
      integral[component] += current[3*cell+component] * weight;
    if(gnd.ksymp != 1)
    {
      dielectric_cell image=mesh->cells[cell];
      image.center[2]=-image.center[2];
      weight=volume*cell_phase_average(&image,direction,k0,1.0);
      double down[3]={direction[0],direction[1],-direction[2]};
      complex double reflected[3];
      /* Only the transverse source component radiates. */
      dielectric_ground_reflect(down,current+3*cell,k0,reflected);
      for(component=0;component<3;++component) integral[component]+=reflected[component]*weight;
    }
  }
  cross3(integral, direction, crossed);
  cross3(crossed, direction, transverse);
  for(component = 0; component < 3; ++component)
    amplitude[component] = I * k0 * VIE_C0 * VIE_MU_0 *
                           transverse[component] / (4.0 * VIE_PI);
}

static double scattered_power(const dielectric_mesh_result *mesh,
                              const complex double *current, double k0)
{
  const size_t count = 974;
  const double golden = VIE_PI * (3.0 - sqrt(5.0));
  const double eta0 = sqrt(VIE_MU_0 / VIE_EPSILON_0);
  double sum = 0.0;
  size_t i;
  for(i = 0; i < count; ++i)
  {
    double z = 1.0 - (gnd.ksymp == 1 ? 2.0:1.0) * ((double)i + 0.5) / (double)count;
    double radial = sqrt(fmax(0.0, 1.0 - z*z));
    double azimuth = golden * (double)i;
    double direction[3] = { radial*cos(azimuth), radial*sin(azimuth), z };
    complex double a[3];
    dielectric_radiation_amplitude(mesh, current, k0, direction, a);
    sum += creal(a[0]*conj(a[0]) + a[1]*conj(a[1]) + a[2]*conj(a[2]));
  }
  return ((gnd.ksymp == 1 ? 4.0:2.0) * VIE_PI / (double)count) * sum / (2.0 * eta0);
}

static void print_pattern(const dielectric_mesh_result *mesh,
                          const complex double *current, double k0,
                          const dielectric_pattern_request *request,
                          double incident_field_squared, FILE *stream)
{
  size_t ip, it;
  if(request == NULL || !request->enabled)
    return;
  fprintf(stream, "\n\n  DIELECTRIC SCATTERED FAR FIELD (R TIMES E, V)"
                  "\n  THETA      PHI          ETHETA(RE,IM)              EPHI(RE,IM)"
                  "        BISTATIC RCS (M^2)");
  for(ip = 0; ip < request->phi_count; ++ip)
    for(it = 0; it < request->theta_count; ++it)
    {
      double td = request->theta_start_degrees + request->theta_step_degrees*(double)it;
      double pd = request->phi_start_degrees + request->phi_step_degrees*(double)ip;
      double t = td * VIE_PI / 180.0, p = pd * VIE_PI / 180.0;
      double direction[3] = {sin(t)*cos(p), sin(t)*sin(p), cos(t)};
      double etheta[3] = {cos(t)*cos(p), cos(t)*sin(p), -sin(t)};
      double ephi[3] = {-sin(p), cos(p), 0.0};
      complex double a[3], at, ap;
      double rcs;
      dielectric_radiation_amplitude(mesh, current, k0, direction, a);
      at = a[0]*etheta[0] + a[1]*etheta[1] + a[2]*etheta[2];
      ap = a[0]*ephi[0] + a[1]*ephi[1] + a[2]*ephi[2];
      rcs = 4.0 * VIE_PI * (creal(at*conj(at)) + creal(ap*conj(ap))) /
            incident_field_squared;
      fprintf(stream, "\n %8.3f %8.3f  %19.12E %19.12E  %19.12E %19.12E  %19.12E",
              td, pd, creal(at), cimag(at), creal(ap), cimag(ap), rcs);
    }
}

int dielectric_solve_plane_wave(const dielectric_model *model,
                                const dielectric_mesh_result *mesh,
                                double frequency_mhz, int iterative,
                                const dielectric_plane_wave *wave,
                                const dielectric_pattern_request *pattern,
                                FILE *stream,
                                dielectric_diagnostic *diagnostic)
{
  size_t cells = mesh->cell_count, unknowns = 3*cells, row, column, cell;
  complex double *matrix = NULL, *original = NULL, *rhs = NULL, *rhs_original = NULL;
  complex double *chi_values = NULL;
  complex double polarization[3] = {0.0, 0.0, 0.0};
  complex double incident_vector[3];
  double propagation[3], frequency_hz, omega, k0, rhs_norm = 0.0, residual = 0.0;
  double pext = 0.0, psca, absorbed;
  vie_kernel_options kernel_options;
  dielectric_fft_system *accelerated=NULL;
  if(wave == NULL || stream == NULL)
  {
    solver_error(diagnostic, model, "SOLVER", "input", "missing plane wave or output stream");
    return 0;
  }
  if(wave->type < 1 || wave->type > 3)
  {
    solver_error(diagnostic, model, "EX", "type", "DF-only solving requires plane-wave EX type 1, 2, or 3");
    return 0;
  }
  if(!dielectric_ground_validate(model,mesh,iterative,wave,pattern,diagnostic))
    return 0;
  dielectric_ground_initialize(frequency_mhz);
  if(gnd.ksymp != 1)
    fprintf(stream,"\n\n  DIELECTRIC BACKGROUND: %s",gnd.iperf==1 ?
      "PERFECT GROUND (GN 1; DIELECTRIC Z >= 0)":"SOMMERFELD GROUND (GN 2; POSITIVE GAP)");
  if(cells == 0)
  {
    fprintf(stream, "\n\n  ---------- DIELECTRIC VOLUME INTEGRAL SOLUTION ----------"
                    "\n  FREQUENCY (MHZ): %12.6f"
                    "\n  VOLUME CELLS:               0"
                    "\n  COMPLEX UNKNOWNS:           0"
                    "\n  RELATIVE RESIDUAL: 0.0000E+00"
                    "\n  VACUUM IDENTITY: NO POLARIZATION CURRENT OR SCATTERED FIELD"
                    "\n\n  TOTAL ELECTRIC DIPOLE MOMENT (C M)"
                    "\n       PX(RE,IM)                 PY(RE,IM)                 PZ(RE,IM)"
                    "\n   0.0000E+00  0.0000E+00   0.0000E+00  0.0000E+00   0.0000E+00  0.0000E+00"
                    "\n\n  DIELECTRIC POWER CHECK"
                    "\n  EXTINCTION POWER (W):  0.00000E+00"
                    "\n  SCATTERED POWER  (W):  0.00000E+00"
                    "\n  RELATIVE IMBALANCE:    0.00000E+00", frequency_mhz);
    double ratio = wave->axial_ratio == 0.0 ? 1.0 : wave->axial_ratio;
    double field_squared = wave->type == 1 ? 1.0 : 1.0 + ratio*ratio;
    print_pattern(mesh, NULL, 2.0*VIE_PI*frequency_mhz*1.0e6/VIE_C0,
                  pattern, field_squared, stream);
    if(pattern != NULL && pattern->near_enabled &&
       !dielectric_print_near_fields(mesh, NULL,
           2.0*VIE_PI*frequency_mhz*1.0e6/VIE_C0, 0, stream))
    {
      solver_error(diagnostic, model, "OUTPUT", "NE/NH", "near-field counts must be positive");
      return 0;
    }
    return 1;
  }
  if(unknowns > SIZE_MAX/unknowns || unknowns*unknowns > SIZE_MAX/sizeof(*matrix))
  {
    solver_error(diagnostic, model, "SOLVER", "matrix", "dense matrix size overflow");
    return 0;
  }
  if(iterative==2) {
    if(!(gnd.ksymp!=1 && gnd.iperf==1 ? dielectric_fft_prepare_pec : dielectric_fft_prepare)(&accelerated,model,mesh,0,
        2.0*VIE_PI*frequency_mhz*1e6/VIE_C0,
        3*unknowns*sizeof(complex double)+cells*sizeof(double),diagnostic)) goto failure;
    fprintf(stream,"\n\n  FFT ACCOUNTED SOLVER BYTES: %zu",dielectric_fft_storage(accelerated));
    fflush(stream);
  } else {
    matrix = calloc(unknowns*unknowns, sizeof(*matrix));
    original = calloc(unknowns*unknowns, sizeof(*original));
  }
  rhs = calloc(unknowns, sizeof(*rhs));
  rhs_original = calloc(unknowns, sizeof(*rhs_original));
  chi_values = calloc(cells, sizeof(*chi_values));
  if((!accelerated && (matrix == NULL || original == NULL)) || rhs == NULL || rhs_original == NULL ||
     chi_values == NULL)
  {
    solver_error(diagnostic, model, "SOLVER", "memory", "unable to allocate dense VIE system");
    goto failure;
  }
  frequency_hz = frequency_mhz * 1.0e6;
  omega = 2.0 * VIE_PI * frequency_hz;
  k0 = omega / VIE_C0;
  dielectric_plane_wave_vectors(wave, propagation, incident_vector);
  vie_kernel_default_options(&kernel_options);
  /* Matrix entries only need to be well below the volume-discretization
     error.  The stricter M4 defaults are retained for kernel certification;
     this tolerance also avoids rejecting symmetry-limited cancellation in
     electrically larger face-adjacent pairs. */
  kernel_options.relative_tolerance = 5.0e-7;
  kernel_options.absolute_tolerance = 5.0e-10;

  for(cell = 0; cell < cells; ++cell)
  {
    const dielectric_material *material = material_by_id(model, mesh->cells[cell].material_id);
    complex double chi, epsilon_r;
    int component;
    if(material == NULL) { solver_error(diagnostic, model, "SOLVER", "material", "mesh references an unknown material"); goto failure; }
    epsilon_r = dielectric_material_epsilon(material,k0);
    chi = epsilon_r - 1.0;
    chi_values[cell] = chi;
    for(component = 0; component < 3; ++component)
      rhs[3*cell+component] = I*omega*VIE_EPSILON_0*chi*background_component(&mesh->cells[cell],propagation,incident_vector,k0,component);
    for(component = 0; component < 3; ++component)
      if(!accelerated) matrix[(3*cell+component)*unknowns + 3*cell+component] = epsilon_r;
  }
  for(cell = 0; !accelerated && cell < cells; ++cell)
  {
    for(column = cell; column < cells; ++column)
    {
      vie_cuboid testing, source;
      complex double nblock[3][3];
      size_t p, q;
      memcpy(testing.center, mesh->cells[cell].center, sizeof(testing.center));
      memcpy(testing.size, mesh->cells[cell].size, sizeof(testing.size));
      memcpy(source.center, mesh->cells[column].center, sizeof(source.center));
      memcpy(source.size, mesh->cells[column].size, sizeof(source.size));
      if(!dielectric_background_n(&testing,&source,k0,&kernel_options,nblock))
      {
        char pair[192];
        snprintf(pair, sizeof(pair),"cells=%zu,%zu k0=%.17g",cell+1,column+1,k0);
        solver_error(diagnostic, model, "SOLVER", pair,
                     "direct or reflected cuboid interaction quadrature failed");
        goto failure;
      }
      for(p = 0; p < 3; ++p)
        for(q = 0; q < 3; ++q)
        {
          matrix[(3*cell+p)*unknowns + 3*column+q] -= chi_values[cell]*nblock[p][q];
          if(column != cell)
          {
            double testing_volume = mesh->cells[cell].size[0]*mesh->cells[cell].size[1]*mesh->cells[cell].size[2];
            double source_volume = mesh->cells[column].size[0]*mesh->cells[column].size[1]*mesh->cells[column].size[2];
            matrix[(3*column+q)*unknowns + 3*cell+p] -=
                chi_values[column]*(testing_volume/source_volume)*nblock[p][q];
          }
        }
    }
  }
  if(!accelerated) memcpy(original, matrix, unknowns*unknowns*sizeof(*matrix));
  memcpy(rhs_original, rhs, unknowns*sizeof(*rhs));
  if(iterative)
  {
    vie_gmres_report report;
    memset(rhs, 0, unknowns*sizeof(*rhs));
    int converged;
    if(accelerated) {
      vie_linear_operator op;dielectric_fft_get_operator(accelerated,&op);
      converged=vie_gmres(&op,rhs_original,rhs,&dielectric_gmres_settings,&report);
    } else converged=vie_dense_gmres(original,unknowns,rhs_original,rhs,&dielectric_gmres_settings,&report);
    residual=report.relative_residual;
    fprintf(stream, "\n\n  SOLVER: %s"
                    "\n  GMRES STATUS: %s"
                    "\n  GMRES ITERATIONS: %zu"
                    "\n  GMRES RESTARTS: %zu"
                    "\n  GMRES OPERATOR APPLICATIONS: %zu"
                    "\n  GMRES TRUE RELATIVE RESIDUAL: %.12E",
            accelerated ? "FFT-OPERATOR GMRES (UNIFORM CUBES; 512 MIB CAP)":"DENSE-OPERATOR GMRES (EXPERIMENTAL; NO ACCELERATION)",
            vie_iterative_status_name(report.status),
            report.iterations, report.restarts, report.operator_applications,
            report.relative_residual);
    if(!converged)
    {
      solver_error(diagnostic, model, "SOLVER", "GMRES", vie_iterative_status_name(report.status));
      goto failure;
    }
  }
  else
  if(!lu_solve(matrix, rhs, unknowns))
  {
    solver_error(diagnostic, model, "SOLVER", "matrix", "singular or numerically rank-deficient VIE system");
    goto failure;
  }
  if(!accelerated) {
  for(row = 0; row < unknowns; ++row)
  {
    complex double value = -rhs_original[row];
    rhs_norm += creal(rhs_original[row]*conj(rhs_original[row]));
    for(column = 0; column < unknowns; ++column)
      value += original[row*unknowns+column]*rhs[column];
    residual += creal(value*conj(value));
  }
  residual = sqrt(residual / fmax(rhs_norm, DBL_MIN));
  }

  fprintf(stream, "\n\n  ---------- DIELECTRIC VOLUME INTEGRAL SOLUTION ----------"
                  "\n  FREQUENCY (MHZ): %12.6f"
                  "\n  VOLUME CELLS:    %12zu"
                  "\n  COMPLEX UNKNOWNS:%12zu"
                  "\n  RELATIVE RESIDUAL:%11.4E",
          frequency_mhz, cells, unknowns, residual);
  fprintf(stream, "\n\n  CELL  REGION MATERIAL       X          Y          Z"
                  "             JX(RE,IM)                 JY(RE,IM)                 JZ(RE,IM)  A/M^2");
  for(cell = 0; cell < cells; ++cell)
  {
    double volume = mesh->cells[cell].size[0]*mesh->cells[cell].size[1]*mesh->cells[cell].size[2];
    int component;
    for(component = 0; component < 3; ++component)
    {
      polarization[component] += rhs[3*cell+component]*volume/(I*omega);
      pext += 0.5*volume*creal(background_component(&mesh->cells[cell],propagation,incident_vector,k0,component)*conj(rhs[3*cell+component]));
    }
    fprintf(stream, "\n %5zu %7d %8d %10.4E %10.4E %10.4E"
                    "  %19.12E %19.12E  %19.12E %19.12E  %19.12E %19.12E",
            cell+1, mesh->cells[cell].region_id, mesh->cells[cell].material_id,
            mesh->cells[cell].center[0], mesh->cells[cell].center[1], mesh->cells[cell].center[2],
            creal(rhs[3*cell]), cimag(rhs[3*cell]), creal(rhs[3*cell+1]), cimag(rhs[3*cell+1]),
            creal(rhs[3*cell+2]), cimag(rhs[3*cell+2]));
  }
  fprintf(stream, "\n\n  TEST DIAGNOSTIC: CELL-AVERAGED INTERNAL ELECTRIC FIELD (V/M)"
                  "\n  CELL             EX(RE,IM)                 EY(RE,IM)                 EZ(RE,IM)");
  for(cell = 0; cell < cells; ++cell)
  {
    const dielectric_material *material = material_by_id(model, mesh->cells[cell].material_id);
    complex double chi = dielectric_material_epsilon(material,k0) - 1.0;
    complex double divisor = I*omega*VIE_EPSILON_0*chi;
    fprintf(stream, "\n %5zu  %11.4E %11.4E  %11.4E %11.4E  %11.4E %11.4E",
            cell+1, creal(rhs[3*cell]/divisor), cimag(rhs[3*cell]/divisor),
            creal(rhs[3*cell+1]/divisor), cimag(rhs[3*cell+1]/divisor),
            creal(rhs[3*cell+2]/divisor), cimag(rhs[3*cell+2]/divisor));
  }
  psca = scattered_power(mesh, rhs, k0);
  absorbed = dielectric_absorbed_power(model,mesh,rhs,k0);
  fprintf(stream,"\n\n  DIELECTRIC ABSORBED POWER (W): %19.12E",absorbed);
  if(gnd.ksymp!=1 && gnd.iperf==2)
  {
    double ground_power,ground_delta;
    if(!dielectric_ground_power(mesh,rhs,k0,0,&ground_power,&ground_delta))
    { solver_error(diagnostic,model,"SOLVER","ground power","lower-half-space flux quadrature did not converge"); goto failure; }
    fprintf(stream,"\n\n  REAL-GROUND SCATTERED POWER (W)"
                   "\n  EXTINCTION POWER: %19.12E"
                   "\n  UPPER-HEMISPHERE SCATTERED POWER: %19.12E"
                   "\n  LOWER-HALF-SPACE SCATTERED POWER: %19.12E"
                   "\n  RELATIVE IMBALANCE: %19.12E"
                   "\n  GROUND FLUX QUADRATURE CHANGE (W): %19.12E",
            pext,psca,ground_power,fabs(pext-psca-ground_power-absorbed)/fmax(fabs(pext),DBL_MIN),ground_delta);
  }
  else
  fprintf(stream, "\n\n  TOTAL ELECTRIC DIPOLE MOMENT (C M)"
                  "\n       PX(RE,IM)                 PY(RE,IM)                 PZ(RE,IM)"
                  "\n  %11.4E %11.4E  %11.4E %11.4E  %11.4E %11.4E"
                  "\n\n  DIELECTRIC POWER CHECK"
                  "\n  EXTINCTION POWER (W): %12.5E"
                  "\n  SCATTERED POWER  (W): %12.5E"
                  "\n  RELATIVE IMBALANCE:   %12.5E",
          creal(polarization[0]), cimag(polarization[0]),
          creal(polarization[1]), cimag(polarization[1]),
          creal(polarization[2]), cimag(polarization[2]), pext, psca,
          fabs(pext-psca-absorbed)/fmax(fmax(fabs(pext), fabs(psca)), DBL_MIN));
  {
    double incident_squared = creal(incident_vector[0]*conj(incident_vector[0]) +
                                     incident_vector[1]*conj(incident_vector[1]) +
                                     incident_vector[2]*conj(incident_vector[2]));
    print_pattern(mesh, rhs, k0, pattern, incident_squared, stream);
  }
  if(pattern != NULL && pattern->near_enabled &&
     !dielectric_print_near_fields(mesh, rhs, k0, 0, stream))
  {
    solver_error(diagnostic, model, "OUTPUT", "NE/NH",
        "near-field evaluation requires positive counts, at least one cell longest-edge clearance from every dielectric cell, and converged background quadrature");
    goto failure;
  }
  dielectric_fft_free(accelerated);
  free(chi_values);
  free(matrix); free(original); free(rhs_original); free(rhs);
  return 1;

failure:
  dielectric_fft_free(accelerated);
  free(chi_values);
  free(matrix); free(original); free(rhs_original); free(rhs);
  return 0;
}
