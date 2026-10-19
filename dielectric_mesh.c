#include "nec2c.h"
#include "dielectric.h"

#include <float.h>
#include <stdint.h>

#define MESH_TOLERANCE 1.0e-12

static void mesh_diagnostic(dielectric_diagnostic *diagnostic,
                            const dielectric_model *model, size_t line,
                            const char *record, const char *value,
                            const char *message)
{
  snprintf(diagnostic->resource, sizeof(diagnostic->resource), "%s",
           model->resource_path != NULL ? model->resource_path : "<dielectric>");
  snprintf(diagnostic->record, sizeof(diagnostic->record), "%s", record);
  snprintf(diagnostic->value, sizeof(diagnostic->value), "%s", value);
  snprintf(diagnostic->message, sizeof(diagnostic->message), "%s", message);
  diagnostic->line = line;
}

static const dielectric_region *find_region(const dielectric_model *model, int id)
{
  size_t index;
  for(index = 0; index < model->region_count; index++)
    if(model->regions[index].id == id)
      return &model->regions[index];
  return NULL;
}

static const dielectric_material *find_material(const dielectric_model *model, int id)
{
  size_t index;
  for(index = 0; index < model->material_count; index++)
    if(model->materials[index].id == id)
      return &model->materials[index];
  return NULL;
}

static int checked_multiply(size_t left, size_t right, size_t *result)
{
  if((left != 0) && (right > SIZE_MAX / left))
    return 0;
  *result = left * right;
  return 1;
}

static int checked_add(size_t left, size_t right, size_t *result)
{
  if(right > SIZE_MAX - left)
    return 0;
  *result = left + right;
  return 1;
}

static int division_count(double length, double maximum_edge, size_t *count)
{
  double ratio = length / maximum_edge;
  if(!isfinite(ratio) || (ratio > (double)SIZE_MAX))
    return 0;
  *count = (size_t)ceil(ratio * (1.0 - 4.0 * DBL_EPSILON));
  if(*count == 0)
    *count = 1;
  return 1;
}

static int box_cell_count(const size_t divisions[3], size_t *count)
{
  size_t xy;
  return checked_multiply(divisions[0], divisions[1], &xy) &&
         checked_multiply(xy, divisions[2], count);
}

static int intervals_overlap_interior(double amin, double amax,
                                      double bmin, double bmax)
{
  return fmax(amin, bmin) < fmin(amax, bmax);
}

static int intervals_touch_or_overlap(double amin, double amax,
                                      double bmin, double bmax)
{
  return fmax(amin, bmin) <= fmin(amax, bmax);
}

static int validate_box_relationships(const dielectric_model *model,
                                      dielectric_diagnostic *diagnostic)
{
  size_t left, right;
  for(left = 0; left < model->box_count; left++)
  {
    const dielectric_region *left_region = find_region(model, model->boxes[left].region_id);
    if(left_region == NULL)
    {
      mesh_diagnostic(diagnostic, model, model->boxes[left].source_line,
                      "BOX", "region", "internal unresolved region reference");
      return 0;
    }
    for(right = left + 1; right < model->box_count; right++)
    {
      const dielectric_region *right_region = find_region(model, model->boxes[right].region_id);
      int interior = 1;
      int closed = 1;
      int axis;
      if(right_region == NULL)
      {
        mesh_diagnostic(diagnostic, model, model->boxes[right].source_line,
                        "BOX", "region", "internal unresolved region reference");
        return 0;
      }
      for(axis = 0; axis < 3; axis++)
      {
        interior &= intervals_overlap_interior(
            model->boxes[left].minimum[axis], model->boxes[left].maximum[axis],
            model->boxes[right].minimum[axis], model->boxes[right].maximum[axis]);
        closed &= intervals_touch_or_overlap(
            model->boxes[left].minimum[axis], model->boxes[left].maximum[axis],
            model->boxes[right].minimum[axis], model->boxes[right].maximum[axis]);
      }
      if(interior)
      {
        char value[64];
        snprintf(value, sizeof(value), "%d,%d", model->boxes[left].region_id,
                 model->boxes[right].region_id);
        mesh_diagnostic(diagnostic, model, model->boxes[right].source_line,
                        "BOX", value, "physical box interiors overlap");
        return 0;
      }
      if(closed && (left_region->material_id != right_region->material_id))
      {
        char value[64];
        snprintf(value, sizeof(value), "%d,%d", left_region->material_id,
                 right_region->material_id);
        mesh_diagnostic(diagnostic, model, model->boxes[right].source_line,
                        "BOX", value,
                        "touching boxes with different material IDs are not supported");
        return 0;
      }
    }
  }
  return 1;
}

static double point_box_distance_squared(const double point[3],
                                         const dielectric_cell *cell)
{
  double distance_squared = 0.0;
  int axis;
  for(axis = 0; axis < 3; axis++)
  {
    double minimum = cell->center[axis] - 0.5 * cell->size[axis];
    double maximum = cell->center[axis] + 0.5 * cell->size[axis];
    double delta = point[axis] < minimum ? minimum - point[axis] :
                   point[axis] > maximum ? point[axis] - maximum : 0.0;
    distance_squared += delta * delta;
  }
  return distance_squared;
}

static double segment_box_distance_squared(const double start[3],
                                           const double end[3],
                                           const dielectric_cell *cell)
{
  double breaks[8];
  size_t break_count = 2;
  size_t index, interval;
  double best = DBL_MAX;
  int axis;
  breaks[0] = 0.0;
  breaks[1] = 1.0;
  for(axis = 0; axis < 3; axis++)
  {
    double direction = end[axis] - start[axis];
    double minimum = cell->center[axis] - 0.5 * cell->size[axis];
    double maximum = cell->center[axis] + 0.5 * cell->size[axis];
    if(direction != 0.0)
    {
      double first = (minimum - start[axis]) / direction;
      double second = (maximum - start[axis]) / direction;
      if((first > 0.0) && (first < 1.0))
        breaks[break_count++] = first;
      if((second > 0.0) && (second < 1.0))
        breaks[break_count++] = second;
    }
  }
  for(index = 1; index < break_count; index++)
  {
    double value = breaks[index];
    size_t cursor = index;
    while((cursor > 0) && (breaks[cursor - 1] > value))
    {
      breaks[cursor] = breaks[cursor - 1];
      cursor--;
    }
    breaks[cursor] = value;
  }
  for(interval = 0; interval + 1 < break_count; interval++)
  {
    double lower = breaks[interval];
    double upper = breaks[interval + 1];
    double middle = 0.5 * (lower + upper);
    double quadratic = 0.0;
    double linear = 0.0;
    double constant = 0.0;
    double candidate;
    for(axis = 0; axis < 3; axis++)
    {
      double direction = end[axis] - start[axis];
      double at_middle = start[axis] + direction * middle;
      double minimum = cell->center[axis] - 0.5 * cell->size[axis];
      double maximum = cell->center[axis] + 0.5 * cell->size[axis];
      double bound;
      if((at_middle >= minimum) && (at_middle <= maximum))
        continue;
      bound = at_middle < minimum ? minimum : maximum;
      quadratic += direction * direction;
      linear += 2.0 * direction * (start[axis] - bound);
      constant += (start[axis] - bound) * (start[axis] - bound);
    }
    candidate = quadratic > 0.0 ? -linear / (2.0 * quadratic) : lower;
    if(candidate < lower)
      candidate = lower;
    if(candidate > upper)
      candidate = upper;
    {
      double value = quadratic * candidate * candidate + linear * candidate + constant;
      if(value < best)
        best = value;
    }
  }
  return best < 0.0 && best > -MESH_TOLERANCE ? 0.0 : best;
}

static int evaluate_clearance(const data_t *geometry,
                              dielectric_mesh_result *mesh,
                              dielectric_diagnostic *diagnostic,
                              const dielectric_model *model)
{
  size_t cell_index;
  double worst_margin = DBL_MAX;
  double worst_gap = 0.0;
  double worst_required = 0.0;
  mesh->has_conductors = (mesh->cell_count > 0) &&
                         ((geometry->n > 0) || (geometry->m > 0));
  mesh->minimum_conductor_gap = DBL_MAX;
  mesh->required_conductor_gap = 0.0;
  if(!mesh->has_conductors)
    return 1;
  for(cell_index = 0; cell_index < mesh->cell_count; cell_index++)
  {
    const dielectric_cell *cell = &mesh->cells[cell_index];
    double cell_longest = fmax(cell->size[0], fmax(cell->size[1], cell->size[2]));
    int conductor;
    for(conductor = 0; conductor < geometry->n; conductor++)
    {
      double start[3], end[3];
      double gap, required;
      start[0] = geometry->x[conductor] - 0.5 * geometry->si[conductor] * geometry->cab[conductor];
      start[1] = geometry->y[conductor] - 0.5 * geometry->si[conductor] * geometry->sab[conductor];
      start[2] = geometry->z[conductor] - 0.5 * geometry->si[conductor] * geometry->salp[conductor];
      end[0] = geometry->x[conductor] + 0.5 * geometry->si[conductor] * geometry->cab[conductor];
      end[1] = geometry->y[conductor] + 0.5 * geometry->si[conductor] * geometry->sab[conductor];
      end[2] = geometry->z[conductor] + 0.5 * geometry->si[conductor] * geometry->salp[conductor];
      gap = sqrt(segment_box_distance_squared(start, end, cell)) - geometry->bi[conductor];
      if(gap < 0.0)
        gap = 0.0;
      required = fmax(cell_longest, geometry->si[conductor]);
      if(gap - required < worst_margin)
      {
        worst_margin = gap - required;
        worst_gap = gap;
        worst_required = required;
      }
      if(gap < mesh->minimum_conductor_gap)
      {
        mesh->minimum_conductor_gap = gap;
        mesh->required_conductor_gap = required;
      }
    }
    for(conductor = 0; conductor < geometry->m; conductor++)
    {
      double point[3] = {geometry->px[conductor], geometry->py[conductor],
                         geometry->pz[conductor]};
      double patch_radius = sqrt(geometry->pbi[conductor] / PI);
      double gap = sqrt(point_box_distance_squared(point, cell)) - patch_radius;
      double required = fmax(cell_longest, sqrt(geometry->pbi[conductor]));
      if(gap < 0.0)
        gap = 0.0;
      if(gap - required < worst_margin)
      {
        worst_margin = gap - required;
        worst_gap = gap;
        worst_required = required;
      }
      if(gap < mesh->minimum_conductor_gap)
      {
        mesh->minimum_conductor_gap = gap;
        mesh->required_conductor_gap = required;
      }
    }
  }
  if(worst_margin < -MESH_TOLERANCE)
  {
    char value[128];
    snprintf(value, sizeof(value), "measured=%.9g required=%.9g",
             worst_gap, worst_required);
    mesh_diagnostic(diagnostic, model, 0, "CLEARANCE", value,
                    "dielectric-to-conductor clearance is insufficient");
    return 0;
  }
  return 1;
}

static int compare_doubles(const void *left, const void *right)
{
  double a = *(const double *)left;
  double b = *(const double *)right;
  return (a > b) - (a < b);
}

void dielectric_mesh_init(dielectric_mesh_result *mesh)
{
  memset(mesh, 0, sizeof(*mesh));
}

void dielectric_mesh_free(dielectric_mesh_result *mesh)
{
  free(mesh->regions);
  free(mesh->boxes);
  free(mesh->cells);
  dielectric_mesh_init(mesh);
}

int dielectric_mesh_build(const dielectric_model *model, double frequency_mhz,
                          const data_t *geometry,
                          dielectric_mesh_result *mesh,
                          dielectric_diagnostic *diagnostic)
{
  dielectric_mesh_result candidate;
  double k0;
  size_t box_index;
  size_t total_cells = 0;
  dielectric_mesh_init(&candidate);
  candidate.frequency_mhz = frequency_mhz;
  if(!isfinite(frequency_mhz) || (frequency_mhz <= 0.0))
  {
    mesh_diagnostic(diagnostic, model, 0, "FR", "frequency",
                    "highest requested frequency must be finite and positive");
    return 0;
  }
  if(model->max_kh > 1.0 + MESH_TOLERANCE)
  {
    char value[64];
    snprintf(value, sizeof(value), "%.17g", model->max_kh);
    mesh_diagnostic(diagnostic, model, 0, "MESH", value,
                    "MAX_KH greater than 1 is not permitted");
    return 0;
  }
  candidate.accuracy_warning = model->max_kh > 0.5 + MESH_TOLERANCE;
  if(!validate_box_relationships(model, diagnostic))
    return 0;
  candidate.regions = calloc(model->region_count, sizeof(*candidate.regions));
  if(candidate.regions == NULL)
  {
    mesh_diagnostic(diagnostic, model, 0, "MESH", "regions",
                    "memory allocation failed");
    return 0;
  }
  candidate.box_count = model->box_count;
  candidate.boxes = calloc(candidate.box_count, sizeof(*candidate.boxes));
  if(candidate.boxes == NULL)
  {
    mesh_diagnostic(diagnostic, model, 0, "MESH", "", "memory allocation failed");
    return 0;
  }
  k0 = TP * frequency_mhz / CVEL;
  for(box_index = 0; box_index < model->region_count; box_index++)
  {
    const dielectric_region *region = &model->regions[box_index];
    const dielectric_material *material = find_material(model, region->material_id);
    double minimum[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double maximum[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    size_t physical_box;
    int has_box = 0;
    int axis;
    if(material == NULL)
    {
      mesh_diagnostic(diagnostic, model, region->source_line, "REGION", "material",
                      "internal unresolved material reference");
      goto failure;
    }
    for(physical_box = 0; physical_box < model->box_count; physical_box++)
    {
      const dielectric_box *box = &model->boxes[physical_box];
      if(box->region_id != region->id)
        continue;
      has_box = 1;
      for(axis = 0; axis < 3; axis++)
      {
        minimum[axis] = fmin(minimum[axis], box->minimum[axis]);
        maximum[axis] = fmax(maximum[axis], box->maximum[axis]);
      }
    }
    if(has_box)
    {
      dielectric_region_mesh *region_mesh =
          &candidate.regions[candidate.region_count++];
      double dx = maximum[0] - minimum[0];
      double dy = maximum[1] - minimum[1];
      double dz = maximum[2] - minimum[2];
      region_mesh->region_id = region->id;
      region_mesh->material_id = material->id;
      region_mesh->diameter = sqrt(dx * dx + dy * dy + dz * dz);
      region_mesh->electrical_diameter =
          k0 * sqrt(material->epsilon_r) * region_mesh->diameter;
      if(!isfinite(region_mesh->electrical_diameter))
      {
        mesh_diagnostic(diagnostic, model, region->source_line, "REGION", "diameter",
                        "region electrical diameter is not representable");
        goto failure;
      }
    }
  }
  for(box_index = 0; box_index < model->box_count; box_index++)
  {
    const dielectric_box *box = &model->boxes[box_index];
    const dielectric_region *region = find_region(model, box->region_id);
    const dielectric_material *material;
    dielectric_box_mesh *box_mesh = &candidate.boxes[box_index];
    double lengths[3] = {0.0, 0.0, 0.0};
    double maximum_edge;
    size_t box_cells;
    int axis;
    if(region == NULL)
    {
      mesh_diagnostic(diagnostic, model, box->source_line, "BOX", "region",
                      "internal unresolved region reference");
      goto failure;
    }
    material = find_material(model, region->material_id);
    if(material == NULL)
    {
      mesh_diagnostic(diagnostic, model, region->source_line, "REGION", "material",
                      "internal unresolved material reference");
      goto failure;
    }
    box_mesh->region_id = region->id;
    box_mesh->material_id = material->id;
    for(axis = 0; axis < 3; axis++)
    {
      lengths[axis] = box->maximum[axis] - box->minimum[axis];
      if(!isfinite(lengths[axis]) || (lengths[axis] <= 0.0))
      {
        mesh_diagnostic(diagnostic, model, box->source_line, "BOX", "bounds",
                        "box edge length is not finite and positive");
        goto failure;
      }
    }
    box_mesh->physical_volume = lengths[0] * lengths[1] * lengths[2];
    if(!isfinite(box_mesh->physical_volume) || (box_mesh->physical_volume <= 0.0) ||
       !isfinite(candidate.declared_volume + box_mesh->physical_volume))
    {
      mesh_diagnostic(diagnostic, model, box->source_line, "BOX", "volume",
                      "box or accumulated physical volume is not representable");
      goto failure;
    }
    candidate.declared_volume += box_mesh->physical_volume;
    if(material->epsilon_r == 1.0)
    {
      candidate.vacuum_box_count++;
      candidate.vacuum_volume += box_mesh->physical_volume;
      continue;
    }
    candidate.active_volume += box_mesh->physical_volume;
    maximum_edge = model->max_kh / (k0 * sqrt(material->epsilon_r));
    if(!isfinite(maximum_edge) || (maximum_edge <= 0.0))
    {
      mesh_diagnostic(diagnostic, model, box->source_line, "MESH", "cell edge",
                      "electrical cell-size limit is not representable");
      goto failure;
    }
    for(axis = 0; axis < 3; axis++)
      if(!division_count(lengths[axis], maximum_edge, &box_mesh->divisions[axis]))
      {
        mesh_diagnostic(diagnostic, model, box->source_line, "BOX", "division count",
                        "required subdivision count is not representable");
        goto failure;
      }
    while(1)
    {
      double edge[3] = {0.0, 0.0, 0.0};
      double smallest, largest;
      int largest_axis = 0;
      for(axis = 0; axis < 3; axis++)
        edge[axis] = lengths[axis] / (double)box_mesh->divisions[axis];
      smallest = fmin(edge[0], fmin(edge[1], edge[2]));
      largest = fmax(edge[0], fmax(edge[1], edge[2]));
      if(largest <= model->max_aspect * smallest * (1.0 + MESH_TOLERANCE))
        break;
      for(axis = 1; axis < 3; axis++)
        if(edge[axis] > edge[largest_axis])
          largest_axis = axis;
      if(box_mesh->divisions[largest_axis] == SIZE_MAX)
      {
        mesh_diagnostic(diagnostic, model, box->source_line, "BOX", "division count",
                        "aspect-ratio refinement overflowed");
        goto failure;
      }
      box_mesh->divisions[largest_axis]++;
      if(!box_cell_count(box_mesh->divisions, &box_cells) ||
         (box_cells > model->max_cells))
      {
        mesh_diagnostic(diagnostic, model, box->source_line, "MESH", "MAX_CELLS",
                        "aspect-ratio refinement exceeds the model cell limit");
        goto failure;
      }
    }
    if(!box_cell_count(box_mesh->divisions, &box_cells) ||
       !checked_add(total_cells, box_cells, &total_cells) ||
       (total_cells > model->max_cells))
    {
      char value[64];
      snprintf(value, sizeof(value), "%lu", (unsigned long)model->max_cells);
      mesh_diagnostic(diagnostic, model, box->source_line, "MESH", value,
                      "required cells exceed MAX_CELLS before cell allocation");
      goto failure;
    }
    box_mesh->cell_offset = total_cells - box_cells;
    box_mesh->cell_count = box_cells;
    if((box_mesh->divisions[0] == 1) || (box_mesh->divisions[1] == 1) ||
       (box_mesh->divisions[2] == 1))
      candidate.feature_warning_count++;
  }
  candidate.cell_count = total_cells;
  if(total_cells > 0)
  {
    candidate.cells = calloc(total_cells, sizeof(*candidate.cells));
    if(candidate.cells == NULL)
    {
      mesh_diagnostic(diagnostic, model, 0, "MESH", "cells", "memory allocation failed");
      goto failure;
    }
  }
  for(box_index = 0; box_index < model->box_count; box_index++)
  {
    const dielectric_box *box = &model->boxes[box_index];
    dielectric_box_mesh *box_mesh = &candidate.boxes[box_index];
    size_t ix, iy, iz, cell_index = box_mesh->cell_offset;
    double edge[3];
    if(box_mesh->cell_count == 0)
      continue;
    for(ix = 0; ix < 3; ix++)
      edge[ix] = (box->maximum[ix] - box->minimum[ix]) /
                 (double)box_mesh->divisions[ix];
    for(ix = 0; ix < box_mesh->divisions[0]; ix++)
      for(iy = 0; iy < box_mesh->divisions[1]; iy++)
        for(iz = 0; iz < box_mesh->divisions[2]; iz++)
        {
          dielectric_cell *cell = &candidate.cells[cell_index++];
          double longest, shortest;
          cell->region_id = box_mesh->region_id;
          cell->material_id = box_mesh->material_id;
          cell->size[0] = edge[0];
          cell->size[1] = edge[1];
          cell->size[2] = edge[2];
          cell->center[0] = box->minimum[0] + ((double)ix + 0.5) * edge[0];
          cell->center[1] = box->minimum[1] + ((double)iy + 0.5) * edge[1];
          cell->center[2] = box->minimum[2] + ((double)iz + 0.5) * edge[2];
          longest = fmax(edge[0], fmax(edge[1], edge[2]));
          shortest = fmin(edge[0], fmin(edge[1], edge[2]));
          {
            const dielectric_material *material = find_material(model, cell->material_id);
            if(material == NULL)
            {
              mesh_diagnostic(diagnostic, model, box->source_line, "CELL", "material",
                              "internal unresolved material reference");
              goto failure;
            }
            cell->q = k0 * sqrt(material->epsilon_r) * longest;
          }
          cell->aspect = longest / shortest;
          box_mesh->numerical_volume += edge[0] * edge[1] * edge[2];
        }
    candidate.numerical_volume += box_mesh->numerical_volume;
  }
  if(candidate.active_volume > 0.0)
    candidate.volume_relative_error = fabs(candidate.numerical_volume - candidate.active_volume) /
                                      candidate.active_volume;
  if(candidate.volume_relative_error > 1.0e-10)
  {
    mesh_diagnostic(diagnostic, model, 0, "MESH", "volume",
                    "numerical cell volume does not conserve physical volume");
    goto failure;
  }
  if(candidate.cell_count > 0)
  {
    double *q_values = malloc(candidate.cell_count * sizeof(*q_values));
    size_t index;
    if(q_values == NULL)
    {
      mesh_diagnostic(diagnostic, model, 0, "MESH", "statistics",
                      "memory allocation failed");
      goto failure;
    }
    candidate.aspect_max = 1.0;
    for(index = 0; index < candidate.cell_count; index++)
    {
      q_values[index] = candidate.cells[index].q;
      if(candidate.cells[index].aspect > candidate.aspect_max)
        candidate.aspect_max = candidate.cells[index].aspect;
    }
    qsort(q_values, candidate.cell_count, sizeof(*q_values), compare_doubles);
    candidate.q_min = q_values[0];
    candidate.q_median = q_values[(candidate.cell_count - 1) / 2];
    candidate.q_p95 = q_values[(size_t)floor(0.95 * (double)(candidate.cell_count - 1))];
    candidate.q_max = q_values[candidate.cell_count - 1];
    free(q_values);
  }
  if(!evaluate_clearance(geometry, &candidate, diagnostic, model))
    goto failure;
  {
    size_t patch_unknowns;
    if(!checked_multiply((size_t)geometry->m, 2, &patch_unknowns) ||
       !checked_add((size_t)geometry->n, patch_unknowns,
                    &candidate.conductor_unknowns))
    {
      mesh_diagnostic(diagnostic, model, 0, "COST", "conductors",
                      "conductor equation count overflow");
      goto failure;
    }
  }
  if(!checked_multiply(candidate.cell_count, 3, &candidate.volume_unknowns) ||
     !checked_add(candidate.conductor_unknowns, candidate.volume_unknowns,
                  &candidate.total_unknowns))
  {
    mesh_diagnostic(diagnostic, model, 0, "COST", "unknowns",
                    "equation count overflow");
    goto failure;
  }
  {
    size_t squared;
    size_t matrix_values;
    size_t rhs_bytes, pivot_bytes, cell_bytes;
    if(!checked_multiply(candidate.total_unknowns, candidate.total_unknowns, &squared) ||
       !checked_multiply(squared, sizeof(complex double), &matrix_values) ||
       !checked_multiply(candidate.total_unknowns, sizeof(complex double), &rhs_bytes) ||
       !checked_multiply(candidate.total_unknowns, sizeof(int), &pivot_bytes) ||
       !checked_multiply(candidate.cell_count, sizeof(dielectric_cell), &cell_bytes) ||
       !checked_add(rhs_bytes, pivot_bytes, &candidate.workspace_bytes) ||
       !checked_add(candidate.workspace_bytes, cell_bytes, &candidate.workspace_bytes))
    {
      mesh_diagnostic(diagnostic, model, 0, "COST", "bytes",
                      "matrix or workspace byte estimate overflow");
      goto failure;
    }
    candidate.matrix_bytes = matrix_values;
  }
  candidate.built = 1;
  dielectric_mesh_free(mesh);
  *mesh = candidate;
  return 1;

failure:
  dielectric_mesh_free(&candidate);
  return 0;
}

void dielectric_print_mesh_report(FILE *stream, const dielectric_mesh_result *mesh)
{
  size_t index;
  fprintf(stream, "\n\n\n"
          "                       ---------- DIELECTRIC MESH ----------\n"
          "  MESH FREQUENCY: %.9g MHz\n"
          "  PHYSICAL BOXES: %lu\n"
          "  ACTIVE CELLS: %lu\n"
          "  VACUUM BOXES OMITTED: %lu\n",
          mesh->frequency_mhz, (unsigned long)mesh->box_count,
          (unsigned long)mesh->cell_count,
          (unsigned long)mesh->vacuum_box_count);
  for(index = 0; index < mesh->region_count; index++)
    fprintf(stream,
            "  REGION_METRIC %d MATERIAL %d DIAMETER %.9g m ELECTRICAL_DIAMETER %.9g\n",
            mesh->regions[index].region_id, mesh->regions[index].material_id,
            mesh->regions[index].diameter,
            mesh->regions[index].electrical_diameter);
  for(index = 0; index < mesh->box_count; index++)
    fprintf(stream, "  BOX_MESH %lu REGION %d MATERIAL %d NX %lu NY %lu NZ %lu CELLS %lu\n",
            (unsigned long)(index + 1), mesh->boxes[index].region_id,
            mesh->boxes[index].material_id,
            (unsigned long)mesh->boxes[index].divisions[0],
            (unsigned long)mesh->boxes[index].divisions[1],
            (unsigned long)mesh->boxes[index].divisions[2],
            (unsigned long)mesh->boxes[index].cell_count);
  fprintf(stream,
          "  DECLARED VOLUME: %.17g m^3\n"
          "  ACTIVE PHYSICAL VOLUME: %.17g m^3\n"
          "  VACUUM VOLUME OMITTED: %.17g m^3\n"
          "  NUMERICAL CELL VOLUME: %.17g m^3\n"
          "  VOLUME RELATIVE ERROR: %.6e\n"
          "  Q MIN/MEDIAN/P95/MAX: %.9g %.9g %.9g %.9g\n"
          "  MAXIMUM CELL ASPECT: %.9g\n"
          "  ONE-CELL FEATURE WARNINGS: %lu\n",
          mesh->declared_volume, mesh->active_volume, mesh->vacuum_volume,
          mesh->numerical_volume, mesh->volume_relative_error,
          mesh->q_min, mesh->q_median, mesh->q_p95, mesh->q_max,
          mesh->aspect_max, (unsigned long)mesh->feature_warning_count);
  if(mesh->accuracy_warning)
    fprintf(stream, "  ACCURACY WARNING: 0.5 < MAX_KH <= 1.0\n");
  if(mesh->has_conductors)
    fprintf(stream, "  MINIMUM CONDUCTOR GAP: %.9g m\n"
                    "  REQUIRED CONDUCTOR GAP: %.9g m\n",
            mesh->minimum_conductor_gap, mesh->required_conductor_gap);
  else
    fprintf(stream, "  CONDUCTOR CLEARANCE: NOT APPLICABLE\n");
  fprintf(stream,
          "  CONDUCTOR UNKNOWNS: %lu\n"
          "  VOLUME UNKNOWNS: %lu\n"
          "  TOTAL UNKNOWNS: %lu\n"
          "  DENSE MATRIX BYTES: %lu\n"
          "  ESTIMATED WORKSPACE BYTES: %lu\n",
          (unsigned long)mesh->conductor_unknowns,
          (unsigned long)mesh->volume_unknowns,
          (unsigned long)mesh->total_unknowns,
          (unsigned long)mesh->matrix_bytes,
          (unsigned long)mesh->workspace_bytes);
  if(mesh->cell_count > 0)
  {
    const dielectric_cell *first = &mesh->cells[0];
    const dielectric_cell *last = &mesh->cells[mesh->cell_count - 1];
    fprintf(stream,
            "  FIRST CELL: R%d M%d C(%.9g,%.9g,%.9g) D(%.9g,%.9g,%.9g)\n"
            "  LAST CELL: R%d M%d C(%.9g,%.9g,%.9g) D(%.9g,%.9g,%.9g)\n",
            first->region_id, first->material_id,
            first->center[0], first->center[1], first->center[2],
            first->size[0], first->size[1], first->size[2],
            last->region_id, last->material_id,
            last->center[0], last->center[1], last->center[2],
            last->size[0], last->size[1], last->size[2]);
  }
  fprintf(stream, "  STATUS: MESH AND COST VALIDATED; READY FOR DIELECTRIC SOLVER");
}
