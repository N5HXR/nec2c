#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "coated_wire.h"
#include "shared.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

static const coated_wire_model *active_model;
static double active_wavelength_m;
static const double coated_wire_m8_warning_threshold = 0.12;

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

static int connected_endpoint(const data_t *geometry, int segment, int other)
{
  if(connected_segment(geometry->icon1[segment], segment, geometry->n) == other)
    return 0;
  if(connected_segment(geometry->icon2[segment], segment, geometry->n) == other)
    return 1;
  return -1;
}

static double segment_alignment(const data_t *geometry, int left, int right)
{
  return geometry->cab[left]*geometry->cab[right] +
      geometry->sab[left]*geometry->sab[right] +
      geometry->salp[left]*geometry->salp[right];
}

/* Length of the straight geometric chord from a bend to the next bend/end.
 * NEC subdivides each GW card into MoM segments, so the validated Y=4 rule
 * must apply to the whole collinear run rather than an individual unknown. */
static double straight_chord_length(const data_t *geometry, int segment,
    int bend_endpoint)
{
  double length = 0.0;
  int current = segment, endpoint = bend_endpoint, visited;
  for(visited = 0; visited < geometry->n; ++visited)
  {
    int next, next_endpoint;
    /* Segment centers/lengths may already be frequency-normalized when XQ is
     * dispatched.  Endpoint arrays retain the input geometry in meters. */
    {
      double dx = geometry->x2[current]-geometry->x1[current];
      double dy = geometry->y2[current]-geometry->y1[current];
      double dz = geometry->z2[current]-geometry->z1[current];
      length += sqrt(dx*dx+dy*dy+dz*dz);
    }
    endpoint = 1-endpoint;
    next = connected_segment(endpoint ? geometry->icon2[current] :
        geometry->icon1[current], current, geometry->n);
    if(next < 0 || fabs(segment_alignment(geometry, current, next)) < 0.999999)
      break;
    next_endpoint = connected_endpoint(geometry, next, current);
    if(next_endpoint < 0) break;
    current = next;
    endpoint = next_endpoint;
  }
  return length;
}

static int wire_component_count(const data_t *geometry)
{
  unsigned char *visited;
  int *stack, components = 0, start;
  if(geometry->n == 0) return 0;
  visited = calloc((size_t)geometry->n, sizeof(*visited));
  stack = malloc((size_t)geometry->n*sizeof(*stack));
  if(!visited || !stack)
  {
    free(visited);
    free(stack);
    return -1;
  }
  for(start = 0; start < geometry->n; ++start)
    if(!visited[start])
    {
      int top = 0;
      ++components;
      visited[start] = 1;
      stack[top++] = start;
      while(top)
      {
        int segment = stack[--top], endpoint;
        for(endpoint = 0; endpoint < 2; ++endpoint)
        {
          int next = connected_segment(endpoint ? geometry->icon2[segment] :
              geometry->icon1[segment], segment, geometry->n);
          if(next >= 0 && !visited[next])
          {
            visited[next] = 1;
            stack[top++] = next;
          }
        }
      }
    }
  free(visited);
  free(stack);
  return components;
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
  model->experimental_geometry = 0;
  model->experimental_geometry_warning_emitted = 0;
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

int coated_wire_validate_geometry(coated_wire_model *model,
    const data_t *geometry, coated_wire_diagnostic *diagnostic)
{
  const coated_wire_domain *domain;
  int i, open_ends = 0, components, unqualified_geometry = 0;

  if(!model || !geometry || !model->finalized)
    return fail(diagnostic, "COATED_WIRE_INTERNAL_ERROR", 0, 0,
        "coated-wire geometry has not been finalized");
  if(model->domain_count == 0 ||
      model->segment_count != (size_t)geometry->n ||
      !model->segment_domain || geometry->m != 0 || geometry->ipsym != 0)
    return fail(diagnostic, "COATED_WIRE_M4_GEOMETRY_UNSUPPORTED", 0, 0,
        "the coated-wire operator requires unsymmetrized, uniformly coated wire geometry and no surface patches");
  domain = &model->domains[0];
  if(geometry->n < 1)
    return fail(diagnostic, "COATED_WIRE_M4_GEOMETRY_UNSUPPORTED", domain->tag, 0,
        "the coated-wire operator requires at least one coated wire segment");
  model->operator_enabled = 0;
  model->experimental_geometry = 0;
  model->qualified_bend_count = 0;
  model->qualified_bend_maximum_turn_deg = 0.0;
  model->qualified_bend_minimum_chord_ratio = DBL_MAX;
  for(i = 0; i < geometry->n; ++i)
  {
    int connections[2] = { geometry->icon1[i], geometry->icon2[i] };
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
    for(endpoint = 0; endpoint < 2; ++endpoint)
    {
      int other = connected_segment(connections[endpoint], i, geometry->n);
      if(connections[endpoint] == 0) ++open_ends;
      if(other >= 0)
      {
        double alignment = segment_alignment(geometry, i, other);
        if(fabs(alignment) < 0.999999 && i < other)
        {
          int other_endpoint = connected_endpoint(geometry, other, i);
          double away_dot, turn_cos, turn_deg, left_length, right_length;
          double outer_radius, minimum_ratio;
          if(model->segment_domain[i] < 0 ||
             model->segment_domain[other] < 0 ||
             !coated_wire_domain_equal(
                 &model->domains[model->segment_domain[i]],
                 &model->domains[model->segment_domain[other]]))
          {
            unqualified_geometry = 1;
            continue;
          }
          if(other_endpoint < 0)
          {
            unqualified_geometry = 1;
            continue;
          }
          /* Direction vectors point from endpoint 1 to endpoint 2.  Reverse
           * the vector when the shared junction is endpoint 2 so both vectors
           * point away from the junction. */
          away_dot = alignment * (endpoint == other_endpoint ? 1.0 : -1.0);
          turn_cos = fmax(-1.0, fmin(1.0, -away_dot));
          turn_deg = acos(turn_cos)*180.0/PI;
          left_length = straight_chord_length(geometry, i, endpoint);
          right_length = straight_chord_length(geometry, other,
              other_endpoint);
          outer_radius = model->domains[model->segment_domain[i]].outer_radius_m;
          minimum_ratio = fmin(left_length, right_length)/outer_radius;
          if(turn_deg > model->qualified_bend_maximum_turn_deg)
            model->qualified_bend_maximum_turn_deg = turn_deg;
          if(minimum_ratio < model->qualified_bend_minimum_chord_ratio)
            model->qualified_bend_minimum_chord_ratio = minimum_ratio;
          if(turn_cos + 1.0e-12 < sqrt(0.5) ||
              minimum_ratio + 1.0e-12 < 4.0)
          {
            unqualified_geometry = 1;
            if(!model->experimental_junctions_enabled)
              return fail(diagnostic, "COATED_WIRE_BEND_UNVALIDATED",
                  geometry->itag[i], i+1,
                  "coated-wire bend exceeds the validated 45-degree maximum turn or four-outer-radius minimum chord; use --coated-wire-junctions experimental to override");
          }
          else
            ++model->qualified_bend_count;
        }
      }
    }
  }
  if(open_ends != 2)
    unqualified_geometry = 1;
  components = wire_component_count(geometry);
  if(components < 0)
    return fail(diagnostic, "COATED_WIRE_ALLOCATION_FAILED", domain->tag, 0,
        "unable to allocate coated-wire geometry traversal");
  if(components != 1) unqualified_geometry = 1;
  if(unqualified_geometry && !model->experimental_junctions_enabled)
    return fail(diagnostic, "COATED_WIRE_GENERAL_GEOMETRY_UNVALIDATED",
        domain->tag, 0,
        "disconnected or non-collinear coated-wire geometry requires --coated-wire-junctions experimental because M7 lacks a converged full-wave reference");
  if(unqualified_geometry) model->experimental_geometry = 1;
  if(model->qualified_bend_count == 0)
    model->qualified_bend_minimum_chord_ratio = 0.0;
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

double coated_wire_validity_metric(const coated_wire_domain *domain,
    double frequency_mhz)
{
  if(!domain || frequency_mhz <= 0.0) return NAN;
  /* NEC uses MHz and CVEL in metres per microsecond, so f/CVEL is 1/metre. */
  return sqrt(domain->epsilon_r)*TP*frequency_mhz/CVEL*
      domain->outer_radius_m;
}

double coated_wire_validity_warning_threshold(void)
{
  return coated_wire_m8_warning_threshold;
}

int coated_wire_check_validity(coated_wire_model *model, double frequency_mhz,
    coated_wire_diagnostic *diagnostic)
{
  double maximum = 0.0;
  size_t domain;
  if(!model || !model->operator_enabled || model->validity_warning_emitted)
    return 0;
  for(domain = 0; domain < model->domain_count; ++domain)
  {
    if(model->domains[domain].outer_radius_m <=
       model->domains[domain].conductor_radius_m)
      continue;
    double metric = coated_wire_validity_metric(&model->domains[domain],
        frequency_mhz);
    if(metric > maximum) maximum = metric;
  }
  if(maximum <= coated_wire_m8_warning_threshold) return 0;
  model->validity_warning_emitted = 1;
  fail_value(diagnostic, "COATED_WIRE_VALIDITY_ENVELOPE_EXCEEDED", 0,
      "sqrt(ER)*k0*outer_radius exceeds the M8 0.12 reference envelope; continuing with reduced-operator results",
      maximum);
  return 1;
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
  if(model->experimental_junctions_enabled)
    fprintf(stream, "  STATUS: EXPERIMENTAL COATING GEOMETRY/JUNCTION MODE REQUESTED\n");
  else
    fprintf(stream, "  STATUS: VALIDATED; STRAIGHT OR QUALIFIED-BEND UNIFORM-COATING OPERATOR ENABLED\n");
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
