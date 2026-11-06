#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "coated_wire.h"
#include "shared.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static const coated_wire_model *active_model;
static double active_wavelength_m;

typedef struct
{
  complex double cosine_charge;
  complex double sine_charge;
} coating_integral;

static const double gauss8_x[4] = {
  0.183434642495649805, 0.525532409916328991,
  0.796666477413626740, 0.960289856497536232
};
static const double gauss8_w[4] = {
  0.362683783378361983, 0.313706645877887287,
  0.222381034453374471, 0.101228536290376259
};
static const double gauss16_x[8] = {
  0.0950125098376374402, 0.281603550779258913,
  0.458016777657227386, 0.617876244402643748,
  0.755404408355003034, 0.865631202387831744,
  0.944575023073232576, 0.989400934991649933
};
static const double gauss16_w[8] = {
  0.189450610455068496, 0.182603415044923589,
  0.169156519395002538, 0.149595988816576732,
  0.124628971255533872, 0.0951585116824927848,
  0.0622535239386478929, 0.0271524594117540949
};

static int fail(coated_wire_diagnostic *diagnostic, const char *code,
    int tag, int segment, const char *message)
{
  if(diagnostic)
  {
    memset(diagnostic, 0, sizeof(*diagnostic));
    diagnostic->code = code;
    diagnostic->tag = tag;
    diagnostic->segment = segment;
    snprintf(diagnostic->message, sizeof(diagnostic->message), "%s", message);
  }
  return 0;
}

static int fail_value(coated_wire_diagnostic *diagnostic, const char *code,
    int tag, const char *message, double value)
{
  fail(diagnostic, code, tag, 0, message);
  if(diagnostic)
  {
    diagnostic->value = value;
    diagnostic->value_defined = 1;
  }
  return 0;
}

void coated_wire_model_init(coated_wire_model *model)
{
  if(model) memset(model, 0, sizeof(*model));
}

void coated_wire_model_clear(coated_wire_model *model)
{
  int experimental_junctions_enabled;
  if(!model) return;
  experimental_junctions_enabled = model->experimental_junctions_enabled;
  free(model->domains);
  free(model->segment_domain);
  memset(model, 0, sizeof(*model));
  model->experimental_junctions_enabled = experimental_junctions_enabled;
}

void coated_wire_allow_experimental_junctions(coated_wire_model *model,
    int enabled)
{
  if(model) model->experimental_junctions_enabled = enabled != 0;
}

int coated_wire_add(coated_wire_model *model, int tag, int first, int last,
    double epsilon_r, double outer_radius_m, double loss_tangent,
    double conductivity_s_per_m, coated_wire_diagnostic *diagnostic)
{
  coated_wire_domain *grown;
  size_t i;
  if(!model)
    return fail(diagnostic, "COATED_WIRE_INTERNAL_ERROR", tag, 0,
        "coated-wire model is unavailable");
  if(tag <= 0)
    return fail(diagnostic, "COATED_WIRE_INVALID_TAG", tag, 0,
        "LD 7 requires a positive, existing wire tag");
  if(first < 0 || last < 0 || ((first == 0) != (last == 0)) ||
      (first != 0 && last < first))
    return fail(diagnostic, "COATED_WIRE_INVALID_RANGE", tag, 0,
        "LD 7 requires 0 0 or an ordered positive tag-relative range");
  if(!isfinite(epsilon_r) || epsilon_r < 1.0)
    return fail_value(diagnostic, "COATED_WIRE_INVALID_ER", tag,
        "LD 7 ER must be finite and at least 1", epsilon_r);
  if(!isfinite(outer_radius_m) || outer_radius_m <= 0.0)
    return fail_value(diagnostic, "COATED_WIRE_INVALID_OUTER_RADIUS", tag,
        "LD 7 OUTER_RADIUS must be finite and positive", outer_radius_m);
  if(!isfinite(loss_tangent) || loss_tangent < 0.0)
    return fail_value(diagnostic, "COATED_WIRE_INVALID_LOSS_TANGENT", tag,
        "LD 7 TAN_DELTA must be finite and nonnegative", loss_tangent);
  if(!isfinite(conductivity_s_per_m) || conductivity_s_per_m < 0.0)
    return fail_value(diagnostic, "COATED_WIRE_INVALID_CONDUCTIVITY", tag,
        "LD 7 SIGMA_DC must be finite and nonnegative", conductivity_s_per_m);
  for(i = 0; i < model->domain_count; ++i)
    if(model->domains[i].tag == tag)
      return fail(diagnostic, "COATED_WIRE_DUPLICATE_TAG", tag, 0,
          "a wire tag may have only one LD 7 coating designation");
  grown = realloc(model->domains,
      (model->domain_count + 1) * sizeof(*model->domains));
  if(!grown)
    return fail(diagnostic, "COATED_WIRE_ALLOCATION_FAILED", tag, 0,
        "unable to allocate coated-wire metadata");
  model->domains = grown;
  memset(&model->domains[model->domain_count], 0,
      sizeof(model->domains[model->domain_count]));
  model->domains[model->domain_count].tag = tag;
  model->domains[model->domain_count].first = first;
  model->domains[model->domain_count].last = last;
  model->domains[model->domain_count].epsilon_r = epsilon_r;
  model->domains[model->domain_count].outer_radius_m = outer_radius_m;
  model->domains[model->domain_count].loss_tangent = loss_tangent;
  model->domains[model->domain_count].conductivity_s_per_m = conductivity_s_per_m;
  ++model->domain_count;
  model->finalized = 0;
  return 1;
}

int coated_wire_domain_equal(const coated_wire_domain *left,
    const coated_wire_domain *right)
{
  return left && right &&
      left->epsilon_r == right->epsilon_r &&
      left->outer_radius_m == right->outer_radius_m &&
      left->loss_tangent == right->loss_tangent &&
      left->conductivity_s_per_m == right->conductivity_s_per_m &&
      left->conductor_radius_m == right->conductor_radius_m;
}

static int connected_segment(int connection, int self, int segment_count)
{
  int absolute;
  if(connection == 0 || connection > PCHCON) return -1;
  absolute = connection < 0 ? -connection : connection;
  if(absolute < 1 || absolute > segment_count) return -1;
  if(absolute - 1 == self) return -2; /* NEC ground-plane self connection. */
  return absolute - 1;
}

int coated_wire_finalize(coated_wire_model *model, const data_t *geometry,
    coated_wire_diagnostic *diagnostic)
{
  size_t d;
  int i;
  if(!model || !geometry)
    return fail(diagnostic, "COATED_WIRE_INTERNAL_ERROR", 0, 0,
        "coated-wire geometry is unavailable");
  if(model->finalized) return 1;
  model->experimental_junction_count = 0;
  free(model->segment_domain);
  model->segment_domain = NULL;
  model->segment_count = (size_t)geometry->n;
  if(geometry->n)
  {
    model->segment_domain = malloc((size_t)geometry->n * sizeof(int));
    if(!model->segment_domain)
      return fail(diagnostic, "COATED_WIRE_ALLOCATION_FAILED", 0, 0,
          "unable to allocate coated-wire segment lookup");
    for(i = 0; i < geometry->n; ++i) model->segment_domain[i] = -1;
  }

  for(d = 0; d < model->domain_count; ++d)
  {
    coated_wire_domain *domain = &model->domains[d];
    int count = 0;
    double radius = 0.0;
    for(i = 0; i < geometry->n; ++i)
      if(geometry->itag[i] == domain->tag)
      {
        ++count;
        if(radius == 0.0) radius = geometry->bi[i];
        else if(geometry->bi[i] != radius)
          return fail(diagnostic, "COATED_WIRE_MIXED_CONDUCTOR_RADIUS",
              domain->tag, i + 1,
              "all segments in an LD 7 tag must have one conductor radius");
      }
    if(count == 0)
      return fail(diagnostic, "COATED_WIRE_TAG_NOT_FOUND", domain->tag, 0,
          "LD 7 references a wire tag that does not exist");
    if(!((domain->first == 0 && domain->last == 0) ||
         (domain->first == 1 && domain->last == count)))
      return fail(diagnostic, "COATED_WIRE_PARTIAL_COVERAGE_UNSUPPORTED",
          domain->tag, 0,
          "LD 7 must cover the complete tag; partial coating is not supported");
    if(domain->outer_radius_m < radius)
      return fail_value(diagnostic, "COATED_WIRE_OUTER_RADIUS_TOO_SMALL",
          domain->tag,
          "LD 7 OUTER_RADIUS must be at least the selected GW conductor radius",
          domain->outer_radius_m);
    domain->first = 1;
    domain->last = count;
    domain->conductor_radius_m = radius;
    domain->segment_count = (size_t)count;
    for(i = 0; i < geometry->n; ++i)
      if(geometry->itag[i] == domain->tag)
      {
        if(model->segment_domain[i] != -1)
          return fail(diagnostic, "COATED_WIRE_SEGMENT_OVERLAP", domain->tag,
              i + 1, "a wire segment is selected by more than one coating");
        model->segment_domain[i] = (int)d;
      }
  }

  for(i = 0; i < geometry->n; ++i)
  {
    int endpoint;
    int own = model->segment_domain ? model->segment_domain[i] : -1;
    if(own < 0) continue;
    for(endpoint = 0; endpoint < 2; ++endpoint)
    {
      int connection = endpoint ? geometry->icon2[i] : geometry->icon1[i];
      int other;
      if(connection > PCHCON)
        return fail(diagnostic, "COATED_WIRE_PATCH_JUNCTION_UNSUPPORTED",
            geometry->itag[i], i + 1,
            "a coated wire may not terminate on a surface patch in this release");
      other = connected_segment(connection, i, geometry->n);
      if(other == -2)
        return fail(diagnostic, "COATED_WIRE_GROUND_JUNCTION_UNSUPPORTED",
            geometry->itag[i], i + 1,
            "a coated wire may not terminate directly on a ground plane in this release");
      if(other >= 0)
      {
        int adjacent = model->segment_domain[other];
        if(adjacent < 0)
        {
          if(!model->experimental_junctions_enabled)
            return fail(diagnostic, "COATED_WIRE_BARE_JUNCTION_UNSUPPORTED",
                geometry->itag[i], i + 1,
                "a coated wire may not contact an uncoated wire in this release");
          ++model->experimental_junction_count;
        }
        if(!coated_wire_domain_equal(&model->domains[own],
                                     adjacent < 0 ? NULL :
                                     &model->domains[adjacent]))
        {
          if(adjacent >= 0 && !model->experimental_junctions_enabled)
            return fail(diagnostic, "COATED_WIRE_JUNCTION_MISMATCH",
                geometry->itag[i], i + 1,
                "contacting coated wires must have identical coating parameters");
          if(adjacent >= 0 && i < other) ++model->experimental_junction_count;
        }
      }
    }
  }
  model->finalized = 1;
  return 1;
}

int coated_wire_validate_isolated(coated_wire_model *model,
    const data_t *geometry, coated_wire_diagnostic *diagnostic)
{
  const coated_wire_domain *domain;
  double ux, uy, uz, x0, y0, z0, scale = 1.0;
  int i, open_ends = 0;

  if(!model || !geometry || !model->finalized)
    return fail(diagnostic, "COATED_WIRE_INTERNAL_ERROR", 0, 0,
        "coated-wire geometry has not been finalized");
  if(model->domain_count == 0 ||
      model->segment_count != (size_t)geometry->n ||
      !model->segment_domain || geometry->m != 0 || geometry->ipsym != 0)
    return fail(diagnostic, "COATED_WIRE_M4_GEOMETRY_UNSUPPORTED", 0, 0,
        "the coated-wire operator requires an unsymmetrized isolated coating covering every wire segment and no surface patches");
  domain = &model->domains[0];
  if(geometry->n < 1)
    return fail(diagnostic, "COATED_WIRE_M4_GEOMETRY_UNSUPPORTED", domain->tag, 0,
        "M4 requires at least one coated wire segment");

  ux = geometry->cab[0];
  uy = geometry->sab[0];
  uz = geometry->salp[0];
  x0 = geometry->x[0];
  y0 = geometry->y[0];
  z0 = geometry->z[0];
  for(i = 0; i < geometry->n; ++i)
  {
    double dot = geometry->cab[i]*ux + geometry->sab[i]*uy +
        geometry->salp[i]*uz;
    double dx = geometry->x[i]-x0, dy = geometry->y[i]-y0,
        dz = geometry->z[i]-z0;
    double axial = dx*ux + dy*uy + dz*uz;
    double radial = sqrt(fmax(0.0, dx*dx+dy*dy+dz*dz-axial*axial));
    int endpoint;
    if(model->segment_domain[i] < 0 &&
       !model->experimental_junctions_enabled)
      return fail(diagnostic, "COATED_WIRE_M4_GEOMETRY_UNSUPPORTED",
          geometry->itag[i], i+1,
          "the coated-wire operator requires coating on every wire segment");
    if(model->segment_domain[i] >= 0 &&
       !coated_wire_domain_equal(domain,
          &model->domains[model->segment_domain[i]]) &&
       !model->experimental_junctions_enabled)
      return fail(diagnostic, "COATED_WIRE_JUNCTION_MISMATCH",
          geometry->itag[i], i+1,
          "one continuous coated wire must use identical coating parameters");
    if(geometry->si[i] > scale) scale = geometry->si[i];
    if(fabs(fabs(dot)-1.0) > 1.0e-10 || radial > scale*1.0e-9)
      return fail(diagnostic, "COATED_WIRE_M4_GEOMETRY_UNSUPPORTED",
          domain->tag, i+1,
          "the coated-wire operator supports one straight collinear coated wire; bends are deferred");
    for(endpoint = 0; endpoint < 2; ++endpoint)
    {
      int connection = endpoint ? geometry->icon2[i] : geometry->icon1[i];
      if(connection == 0) ++open_ends;
    }
  }
  if(open_ends != 2)
    return fail(diagnostic, "COATED_WIRE_M4_GEOMETRY_UNSUPPORTED",
        domain->tag, 0,
        "the coated-wire operator requires one connected coated wire with exactly two free ends");
  model->operator_enabled = 1;
  return 1;
}

void coated_wire_activate(const coated_wire_model *model, double wavelength_m)
{
  active_model = model && model->finalized && model->domain_count &&
      model->operator_enabled ? model : NULL;
  active_wavelength_m = wavelength_m;
}

complex double coated_wire_relative_permittivity(
    const coated_wire_domain *domain, double frequency_mhz)
{
  const double epsilon0 = 8.8541878128e-12;
  double loss;
  if(!domain || frequency_mhz <= 0.0) return NAN + I*NAN;
  loss = domain->epsilon_r*domain->loss_tangent +
      domain->conductivity_s_per_m/
      (TP*frequency_mhz*1.0e6*epsilon0);
  return cmplx(domain->epsilon_r, -loss);
}

int coated_wire_has_loss(void)
{
  size_t domain;
  if(!active_model) return 0;
  for(domain = 0; domain < active_model->domain_count; ++domain)
    if(active_model->domains[domain].loss_tangent != 0.0 ||
        active_model->domains[domain].conductivity_s_per_m != 0.0)
      return 1;
  return 0;
}

int coated_wire_segment_active(int segment)
{
  return active_model && segment >= 1 &&
      (size_t)segment <= active_model->segment_count &&
      active_model->segment_domain[segment-1] >= 0;
}

static coating_integral coating_integrand(int source, int observation,
    double position, double inner_radius, double outer_radius)
{
  coating_integral value;
  double sx = data.x[source] + data.cab[source]*position;
  double sy = data.y[source] + data.sab[source]*position;
  double sz = data.z[source] + data.salp[source]*position;
  double dx = data.x[observation]-sx;
  double dy = data.y[observation]-sy;
  double dz = data.z[observation]-sz;
  double projection = dx*data.cab[observation] +
      dy*data.sab[observation] + dz*data.salp[observation];
  double distance2 = dx*dx + dy*dy + dz*dz;
  double ra2 = distance2 + inner_radius*inner_radius;
  double rb2 = distance2 + outer_radius*outer_radius;
  double ra = sqrt(ra2), rb = sqrt(rb2);
  complex double ga = -cmplx(1.0, TP*ra)*
      cmplx(cos(TP*ra), -sin(TP*ra))/(ra2*ra);
  complex double gb = -cmplx(1.0, TP*rb)*
      cmplx(cos(TP*rb), -sin(TP*rb))/(rb2*rb);
  complex double projected_gradient = projection*(gb-ga);

  value.cosine_charge = cos(TP*position)*projected_gradient;
  value.sine_charge = sin(TP*position)*projected_gradient;
  return value;
}

static coating_integral coating_gauss(int source, int observation,
    double inner_radius, double outer_radius, double lower, double upper,
    const double *nodes, const double *weights, size_t count)
{
  coating_integral sum = { CPLX_00, CPLX_00 };
  double midpoint = 0.5*(lower+upper), half = 0.5*(upper-lower);
  size_t i;
  for(i = 0; i < count; ++i)
  {
    coating_integral left = coating_integrand(source, observation,
        midpoint-half*nodes[i], inner_radius, outer_radius);
    coating_integral right = coating_integrand(source, observation,
        midpoint+half*nodes[i], inner_radius, outer_radius);
    sum.cosine_charge += weights[i]*(left.cosine_charge+right.cosine_charge);
    sum.sine_charge += weights[i]*(left.sine_charge+right.sine_charge);
  }
  sum.cosine_charge *= half;
  sum.sine_charge *= half;
  return sum;
}

static coating_integral coating_adaptive(int source, int observation,
    double inner_radius, double outer_radius, double lower, double upper,
    int depth)
{
  coating_integral q8 = coating_gauss(source, observation, inner_radius,
      outer_radius, lower, upper, gauss8_x, gauss8_w, 4);
  coating_integral q16 = coating_gauss(source, observation, inner_radius,
      outer_radius, lower, upper, gauss16_x, gauss16_w, 8);
  double error = fmax(cabs(q16.cosine_charge-q8.cosine_charge),
      cabs(q16.sine_charge-q8.sine_charge));
  double scale = fmax(cabs(q16.cosine_charge), cabs(q16.sine_charge));
  if(depth >= 14 || error <= 2.0e-12*(1.0+scale)) return q16;
  {
    double midpoint = 0.5*(lower+upper);
    coating_integral left = coating_adaptive(source, observation,
        inner_radius, outer_radius, lower, midpoint, depth+1);
    coating_integral right = coating_adaptive(source, observation,
        inner_radius, outer_radius, midpoint, upper, depth+1);
    left.cosine_charge += right.cosine_charge;
    left.sine_charge += right.sine_charge;
    return left;
  }
}

void coated_wire_field_correction(int source_segment, int observation_segment,
    complex double *sine_field, complex double *cosine_field)
{
  const coated_wire_domain *domain;
  coating_integral integral;
  complex double epsilon, alpha;
  double inner_radius, outer_radius, half_length, frequency_mhz;
  int source = source_segment-1, observation = observation_segment-1;

  if(!active_model || !sine_field || !cosine_field || source < 0 ||
      observation < 0 || (size_t)source >= active_model->segment_count ||
      (size_t)observation >= active_model->segment_count ||
      active_model->segment_domain[source] < 0 || active_wavelength_m <= 0.0)
    return;
  domain = &active_model->domains[active_model->segment_domain[source]];
  frequency_mhz = CVEL/active_wavelength_m;
  epsilon = coated_wire_relative_permittivity(domain, frequency_mhz);
  alpha = 1.0-1.0/epsilon;
  inner_radius = data.bi[source];
  outer_radius = domain->outer_radius_m/active_wavelength_m;
  if(alpha == CPLX_00 ||
      domain->outer_radius_m == domain->conductor_radius_m ||
      outer_radius == inner_radius) return;
  half_length = 0.5*data.si[source];
  integral = coating_adaptive(source, observation, inner_radius, outer_radius,
      -half_length, half_length, 0);
  *sine_field += alpha*CONST3*integral.cosine_charge;
  *cosine_field -= alpha*CONST3*integral.sine_charge;
}

double coated_wire_absorbed_power(void)
{
  const double epsilon0 = 8.8541878128e-12;
  double frequency_mhz, frequency_hz, omega, total = 0.0;
  size_t segment;
  if(!coated_wire_has_loss() || active_wavelength_m <= 0.0) return 0.0;
  frequency_mhz = CVEL/active_wavelength_m;
  frequency_hz = frequency_mhz*1.0e6;
  omega = TP*frequency_hz;
  for(segment = 0; segment < active_model->segment_count; ++segment)
  {
    int domain_index = active_model->segment_domain[segment];
    const coated_wire_domain *domain;
    complex double epsilon, b, c;
    double half, oscillatory, cosine_norm, sine_norm, charge_integral;
    if(domain_index < 0) continue;
    domain = &active_model->domains[domain_index];
    epsilon = coated_wire_relative_permittivity(domain, frequency_mhz);
    if(cimag(epsilon) == 0.0 ||
        domain->outer_radius_m == domain->conductor_radius_m) continue;
    b = cmplx(crnt.bir[segment], crnt.bii[segment]);
    c = cmplx(crnt.cir[segment], crnt.cii[segment]);
    half = 0.5*data.si[segment];
    oscillatory = sin(2.0*TP*half)/(2.0*TP);
    cosine_norm = half+oscillatory;
    sine_norm = half-oscillatory;
    charge_integral = active_wavelength_m/(frequency_hz*frequency_hz)*
        (creal(b*conj(b))*cosine_norm + creal(c*conj(c))*sine_norm);
    total += omega*(-cimag(epsilon))*charge_integral*
        log(domain->outer_radius_m/domain->conductor_radius_m)/
        (4.0*PI*epsilon0*creal(epsilon*conj(epsilon)));
  }
  return total;
}

double coated_wire_radiated_power(size_t sample_count)
{
  double sum = 0.0;
  size_t sample;
  if(!active_model || sample_count == 0) return NAN;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum) schedule(static) if(sample_count >= 1024)
#endif
  for(sample = 0; sample < sample_count; ++sample)
  {
    double z = 1.0-2.0*((double)sample+0.5)/(double)sample_count;
    double phi = PI*(3.0-sqrt(5.0))*(double)sample;
    complex double et, ep;
    ffld(acos(z), phi, &et, &ep);
    et *= active_wavelength_m;
    ep *= active_wavelength_m;
    sum += creal(et*conj(et)+ep*conj(ep));
  }
  return sum*(4.0*PI/(double)sample_count)/(2.0*ETA);
}

void coated_wire_print_summary(FILE *stream, const coated_wire_model *model)
{
  size_t i;
  if(!stream || !model || !model->domain_count) return;
  fprintf(stream, "\n\n\n"
      "                         -------- COATED WIRES --------\n"
      "   TAG   SEGMENTS   WIRE RADIUS   OUTER RADIUS       ER"
      "    TAN DELTA    SIGMA DC\n");
  for(i = 0; i < model->domain_count; ++i)
  {
    const coated_wire_domain *d = &model->domains[i];
    fprintf(stream, " %5d %4d THROUGH %-4d %12.5E %12.5E %9.5g"
        " %12.5E %12.5E\n", d->tag, d->first, d->last,
        d->conductor_radius_m, d->outer_radius_m, d->epsilon_r,
        d->loss_tangent, d->conductivity_s_per_m);
  }
  if(model->experimental_junction_count)
    fprintf(stream, "  STATUS: EXPERIMENTAL COATING JUNCTION OPERATOR REQUESTED\n");
  else
    fprintf(stream, "  STATUS: VALIDATED; STRAIGHT UNIFORM-COATING OPERATOR ENABLED\n");
}

void coated_wire_print_diagnostic(FILE *stream,
    const coated_wire_diagnostic *diagnostic)
{
  if(!stream || !diagnostic) return;
  fprintf(stream, "\n\n  COATED WIRE ERROR [%s]: %s",
      diagnostic->code ? diagnostic->code : "COATED_WIRE_ERROR",
      diagnostic->message);
  if(diagnostic->tag) fprintf(stream, " [tag: %d]", diagnostic->tag);
  if(diagnostic->segment) fprintf(stream, " [segment: %d]", diagnostic->segment);
  if(diagnostic->value_defined) fprintf(stream, " [value: %.17g]", diagnostic->value);
}
