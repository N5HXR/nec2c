#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "coated_wire.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

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
  if(!model) return;
  free(model->domains);
  free(model->segment_domain);
  memset(model, 0, sizeof(*model));
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
          return fail(diagnostic, "COATED_WIRE_BARE_JUNCTION_UNSUPPORTED",
              geometry->itag[i], i + 1,
              "a coated wire may not contact an uncoated wire in this release");
        if(!coated_wire_domain_equal(&model->domains[own],
                                     &model->domains[adjacent]))
          return fail(diagnostic, "COATED_WIRE_JUNCTION_MISMATCH",
              geometry->itag[i], i + 1,
              "contacting coated wires must have identical coating parameters");
      }
    }
  }
  model->finalized = 1;
  return 1;
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
  fprintf(stream, "  STATUS: VALIDATED; NUMERICAL OPERATOR DEFERRED TO M3\n");
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
