#include "vie_kernel.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751
#endif

#define VIE_MAX_ORDER 64
#define VIE_MAX_SEGMENTS 5
#define VIE_GEOMETRY_TOLERANCE 64.0

typedef struct
{
  int both_vary;
  double lower;
  double upper;
  double a_min;
  double a_max;
  double b_min;
  double b_max;
} difference_segment;

typedef struct
{
  difference_segment segment[VIE_MAX_SEGMENTS];
  int count;
} segment_list;

typedef struct
{
  int dimensions;
  difference_segment segment[3];
  double fixed_squared;
  double k0;
  const double *nodes;
  const double *weights;
  int order;
  int radial_power;
} reduced_integral;

static double maximum3(double a, double b, double c)
{
  return fmax(a, fmax(b, c));
}

static int cuboid_valid(const vie_cuboid *cuboid)
{
  int axis;
  if(cuboid == NULL)
    return 0;
  for(axis = 0; axis < 3; axis++)
    if(!isfinite(cuboid->center[axis]) || !isfinite(cuboid->size[axis]) ||
       (cuboid->size[axis] <= 0.0))
      return 0;
  return 1;
}

void vie_kernel_default_options(vie_kernel_options *options)
{
  options->relative_tolerance = 2.0e-8;
  options->absolute_tolerance = 2.0e-11;
  options->maximum_order = 64;
  options->far_order = 8;
}

const char *vie_relation_name(vie_relation relation)
{
  switch(relation)
  {
    case VIE_RELATION_SELF: return "self";
    case VIE_RELATION_FACE: return "face-adjacent";
    case VIE_RELATION_EDGE: return "edge-adjacent";
    case VIE_RELATION_CORNER: return "corner-adjacent";
    case VIE_RELATION_NEAR: return "near";
    case VIE_RELATION_FAR: return "far";
  }
  return "unknown";
}

static vie_relation classify_relation(const vie_cuboid *testing,
                                      const vie_cuboid *source)
{
  int axis;
  int identical = 1;
  int separated_axes = 0;
  int touching_axes = 0;
  double distance_squared = 0.0;
  double scale = maximum3(testing->size[0], testing->size[1], testing->size[2]);
  scale = fmax(scale, maximum3(source->size[0], source->size[1], source->size[2]));
  for(axis = 0; axis < 3; axis++)
  {
    double center_distance = fabs(testing->center[axis] - source->center[axis]);
    double half_sum = 0.5 * (testing->size[axis] + source->size[axis]);
    double tolerance = VIE_GEOMETRY_TOLERANCE * DBL_EPSILON *
                       fmax(1.0, fmax(center_distance, half_sum));
    double gap = center_distance - half_sum;
    if((center_distance > tolerance) ||
       (fabs(testing->size[axis] - source->size[axis]) > tolerance))
      identical = 0;
    if(gap > tolerance)
    {
      separated_axes++;
      distance_squared += gap * gap;
    }
    else if(fabs(gap) <= tolerance)
      touching_axes++;
  }
  if(identical)
    return VIE_RELATION_SELF;
  if(separated_axes == 0)
  {
    if(touching_axes == 1) return VIE_RELATION_FACE;
    if(touching_axes == 2) return VIE_RELATION_EDGE;
    if(touching_axes == 3) return VIE_RELATION_CORNER;
    return VIE_RELATION_NEAR;
  }
  return sqrt(distance_squared) > 2.0 * scale ? VIE_RELATION_FAR : VIE_RELATION_NEAR;
}

static int interiors_overlap(const vie_cuboid *testing,
                             const vie_cuboid *source)
{
  int axis;
  int identical = 1;
  int overlap = 1;
  for(axis = 0; axis < 3; axis++)
  {
    double center_distance = fabs(testing->center[axis] - source->center[axis]);
    double half_sum = 0.5 * (testing->size[axis] + source->size[axis]);
    double tolerance = VIE_GEOMETRY_TOLERANCE * DBL_EPSILON *
                       fmax(1.0, fmax(center_distance, half_sum));
    identical &= (center_distance <= tolerance) &&
                 (fabs(testing->size[axis] - source->size[axis]) <= tolerance);
    overlap &= center_distance < half_sum - tolerance;
  }
  return overlap && !identical;
}

static int gauss_legendre(int order, double *nodes, double *weights)
{
  int root;
  int half = (order + 1) / 2;
  if((order < 2) || (order > VIE_MAX_ORDER))
    return 0;
  for(root = 0; root < half; root++)
  {
    double z = cos(PI * ((double)root + 0.75) / ((double)order + 0.5));
    double previous;
    double derivative = 0.0;
    int iteration = 0;
    do
    {
      double p0 = 1.0;
      double p1 = z;
      int degree;
      previous = z;
      for(degree = 2; degree <= order; degree++)
      {
        double p2 = (((2.0 * degree) - 1.0) * z * p1 -
                     ((double)degree - 1.0) * p0) / (double)degree;
        p0 = p1;
        p1 = p2;
      }
      derivative = (double)order * (z * p1 - p0) / (z * z - 1.0);
      z = previous - p1 / derivative;
      iteration++;
    } while((fabs(z - previous) > 4.0 * DBL_EPSILON) && (iteration < 32));
    if(iteration == 32)
      return 0;
    nodes[root] = -z;
    nodes[order - 1 - root] = z;
    weights[root] = 2.0 / ((1.0 - z * z) * derivative * derivative);
    weights[order - 1 - root] = weights[root];
  }
  return 1;
}

static double overlap_weight(const difference_segment *segment, double difference)
{
  double lower = fmax(segment->a_min, segment->b_min + difference);
  double upper = fmin(segment->a_max, segment->b_max + difference);
  return upper > lower ? upper - lower : 0.0;
}

static double segment_weight(const difference_segment *segment, double coordinate)
{
  return segment->both_vary ? overlap_weight(segment, coordinate) : 1.0;
}

static complex double green_function(double radius, double k0)
{
  return cexp(-I * k0 * radius) / (4.0 * PI * radius);
}

static complex double tensor_integral_recursive(const reduced_integral *integral,
                                                int dimension,
                                                double coordinate[3],
                                                double accumulated_weight)
{
  complex double sum = 0.0;
  const difference_segment *segment = &integral->segment[dimension];
  double midpoint = 0.5 * (segment->lower + segment->upper);
  double halfwidth = 0.5 * (segment->upper - segment->lower);
  int node;
  for(node = 0; node < integral->order; node++)
  {
    double value = midpoint + halfwidth * integral->nodes[node];
    double weight = accumulated_weight * halfwidth * integral->weights[node] *
                    segment_weight(segment, value);
    coordinate[dimension] = value;
    if(dimension + 1 == integral->dimensions)
    {
      double radius_squared = integral->fixed_squared;
      int axis;
      for(axis = 0; axis < integral->dimensions; axis++)
        radius_squared += coordinate[axis] * coordinate[axis];
      sum += weight * green_function(sqrt(radius_squared), integral->k0);
    }
    else
      sum += tensor_integral_recursive(integral, dimension + 1,
                                       coordinate, weight);
  }
  return sum;
}

static complex double duffy_integral_recursive(const reduced_integral *integral,
                                               int sector, int dimension,
                                               double normalized[3],
                                               double accumulated_weight)
{
  complex double sum = 0.0;
  int node;
  for(node = 0; node < integral->order; node++)
  {
    double u = 0.5 * (1.0 + integral->nodes[node]);
    double weight = accumulated_weight * 0.5 * integral->weights[node];
    normalized[dimension] = u;
    if(dimension + 1 == integral->dimensions)
    {
      double coordinate[3] = {0.0, 0.0, 0.0};
      double radial_parameter = normalized[sector];
      double t = pow(radial_parameter, (double)integral->radial_power);
      double radius_squared = integral->fixed_squared;
      double jacobian = 1.0;
      double physical_weight = 1.0;
      int axis;
      for(axis = 0; axis < integral->dimensions; axis++)
      {
        const difference_segment *segment = &integral->segment[axis];
        double extent = fabs(segment->upper - segment->lower);
        double sign = segment->upper == 0.0 ? -1.0 : 1.0;
        double mapped = axis == sector ? t : t * normalized[axis];
        coordinate[axis] = sign * extent * mapped;
        radius_squared += coordinate[axis] * coordinate[axis];
        jacobian *= extent;
        physical_weight *= segment_weight(segment, coordinate[axis]);
      }
      jacobian *= pow(t, (double)(integral->dimensions - 1)) *
                  (double)integral->radial_power *
                  pow(radial_parameter, (double)(integral->radial_power - 1));
      sum += weight * jacobian * physical_weight *
             green_function(sqrt(radius_squared), integral->k0);
    }
    else
      sum += duffy_integral_recursive(integral, sector, dimension + 1,
                                      normalized, weight);
  }
  return sum;
}

static int segment_has_zero_endpoint(const difference_segment *segment)
{
  return segment->lower == 0.0 || segment->upper == 0.0;
}

static complex double integrate_reduced_box(const reduced_integral *integral)
{
  double coordinate[3] = {0.0, 0.0, 0.0};
  int singular = 1;
  int axis;
  for(axis = 0; axis < integral->dimensions; axis++)
    singular &= segment_has_zero_endpoint(&integral->segment[axis]);
  if(!singular)
    return tensor_integral_recursive(integral, 0, coordinate, 1.0);
  else
  {
    complex double sum = 0.0;
    for(axis = 0; axis < integral->dimensions; axis++)
    {
      double normalized[3] = {0.0, 0.0, 0.0};
      normalized[axis] = 0.0;
      sum += duffy_integral_recursive(integral, axis, 0, normalized, 1.0);
    }
    return sum;
  }
}

static int append_break(double breaks[6], int *count, double value)
{
  int index;
  for(index = 0; index < *count; index++)
    if(value == breaks[index])
      return 1;
  if(*count >= 6)
    return 0;
  breaks[(*count)++] = value;
  return 1;
}

static void sort_breaks(double *breaks, int count)
{
  int index;
  for(index = 1; index < count; index++)
  {
    double value = breaks[index];
    int cursor = index;
    while((cursor > 0) && (breaks[cursor - 1] > value))
    {
      breaks[cursor] = breaks[cursor - 1];
      cursor--;
    }
    breaks[cursor] = value;
  }
}

static int make_segments(double a_min, double a_max, int a_varies,
                         double b_min, double b_max, int b_varies,
                         segment_list *list, double *fixed)
{
  double breaks[6];
  int break_count = 0;
  int index;
  memset(list, 0, sizeof(*list));
  if(!a_varies && !b_varies)
  {
    *fixed = a_min - b_min;
    return 1;
  }
  if(a_varies && b_varies)
  {
    if(!append_break(breaks, &break_count, a_min - b_max) ||
       !append_break(breaks, &break_count, a_min - b_min) ||
       !append_break(breaks, &break_count, a_max - b_max) ||
       !append_break(breaks, &break_count, a_max - b_min))
      return 0;
  }
  else if(a_varies)
  {
    append_break(breaks, &break_count, a_min - b_min);
    append_break(breaks, &break_count, a_max - b_min);
  }
  else
  {
    append_break(breaks, &break_count, a_min - b_max);
    append_break(breaks, &break_count, a_min - b_min);
  }
  sort_breaks(breaks, break_count);
  if((breaks[0] < 0.0) && (breaks[break_count - 1] > 0.0))
  {
    append_break(breaks, &break_count, 0.0);
    sort_breaks(breaks, break_count);
  }
  for(index = 0; index + 1 < break_count; index++)
  {
    difference_segment *segment;
    if(breaks[index + 1] == breaks[index])
      continue;
    if(list->count >= VIE_MAX_SEGMENTS)
      return 0;
    segment = &list->segment[list->count++];
    segment->both_vary = a_varies && b_varies;
    segment->lower = breaks[index];
    segment->upper = breaks[index + 1];
    segment->a_min = a_min;
    segment->a_max = a_max;
    segment->b_min = b_min;
    segment->b_max = b_max;
  }
  return 1;
}

static complex double integrate_segment_product(const segment_list lists[3],
                                                const double fixed[3],
                                                double k0, int order,
                                                const double *nodes,
                                                const double *weights)
{
  complex double sum = 0.0;
  int i, j, k;
  for(i = 0; i < (lists[0].count ? lists[0].count : 1); i++)
    for(j = 0; j < (lists[1].count ? lists[1].count : 1); j++)
      for(k = 0; k < (lists[2].count ? lists[2].count : 1); k++)
      {
        int selected[3] = {i, j, k};
        reduced_integral integral;
        int axis;
        memset(&integral, 0, sizeof(integral));
        integral.k0 = k0;
        integral.order = order;
        integral.nodes = nodes;
        integral.weights = weights;
        integral.radial_power = 1;
        for(axis = 0; axis < 3; axis++)
        {
          if(lists[axis].count)
            integral.segment[integral.dimensions++] =
                lists[axis].segment[selected[axis]];
          else
            integral.fixed_squared += fixed[axis] * fixed[axis];
        }
        if(integral.fixed_squared > 0.0)
        {
          double largest_extent = 0.0;
          double ratio;
          for(axis = 0; axis < integral.dimensions; axis++)
            largest_extent = fmax(largest_extent,
                                  fabs(integral.segment[axis].upper -
                                       integral.segment[axis].lower));
          ratio = sqrt(integral.fixed_squared) / largest_extent;
          while((ratio < 0.1) && (integral.radial_power < 8))
          {
            integral.radial_power++;
            ratio *= 10.0;
          }
        }
        sum += integrate_reduced_box(&integral);
      }
  return sum;
}

static complex double face_pair_integral(const vie_cuboid *testing,
                                         const vie_cuboid *source,
                                         int test_axis, int test_sign,
                                         int source_axis, int source_sign,
                                         double k0, int order,
                                         const double *nodes,
                                         const double *weights)
{
  segment_list lists[3];
  double fixed[3] = {0.0, 0.0, 0.0};
  int axis;
  for(axis = 0; axis < 3; axis++)
  {
    int test_varies = axis != test_axis;
    int source_varies = axis != source_axis;
    double test_min = test_varies ? -0.5 * testing->size[axis]
                                   : 0.5 * test_sign * testing->size[axis];
    double test_max = test_varies ? 0.5 * testing->size[axis] : test_min;
    double source_center = source->center[axis] - testing->center[axis];
    double source_min = source_varies ? source_center - 0.5 * source->size[axis]
                                      : source_center + 0.5 * source_sign * source->size[axis];
    double source_max = source_varies ? source_center + 0.5 * source->size[axis]
                                      : source_min;
    if(!make_segments(test_min, test_max, test_varies,
                      source_min, source_max, source_varies,
                      &lists[axis], &fixed[axis]))
      return CMPLX(NAN, NAN);
  }
  return integrate_segment_product(lists, fixed, k0, order, nodes, weights);
}

static int matrix_at_order(const vie_cuboid *testing, const vie_cuboid *source,
                           double k0, int order, complex double result[3][3])
{
  double nodes[VIE_MAX_ORDER], weights[VIE_MAX_ORDER];
  double test_volume = testing->size[0] * testing->size[1] * testing->size[2];
  int test_axis, source_axis, test_sign, source_sign;
  memset(result, 0, 9 * sizeof(complex double));
  if(!gauss_legendre(order, nodes, weights))
    return 0;
  for(test_axis = 0; test_axis < 3; test_axis++)
    for(source_axis = 0; source_axis < 3; source_axis++)
      for(test_sign = -1; test_sign <= 1; test_sign += 2)
        for(source_sign = -1; source_sign <= 1; source_sign += 2)
        {
          complex double integral = face_pair_integral(
              testing, source, test_axis, test_sign, source_axis, source_sign,
              k0, order, nodes, weights) / test_volume;
          if(!isfinite(creal(integral)) || !isfinite(cimag(integral)))
            return 0;
          if(test_axis == source_axis)
          {
            int component;
            for(component = 0; component < 3; component++)
              if(component != test_axis)
                result[component][component] +=
                    (double)(test_sign * source_sign) * integral;
          }
          else
            result[source_axis][test_axis] -=
                (double)(test_sign * source_sign) * integral;
        }
  return 1;
}

static double matrix_difference(const complex double *left,
                                const complex double *right,
                                double absolute_tolerance,
                                double relative_tolerance,
                                double *absolute_error, int *worst_index)
{
  double worst = 0.0;
  int index;
  for(index = 0; index < 9; index++)
  {
    double scale = absolute_tolerance + relative_tolerance *
                   fmax(cabs(left[index]), cabs(right[index]));
    double ratio = cabs(left[index] - right[index]) / scale;
    if(ratio > worst)
    {
      worst = ratio;
      if(absolute_error != NULL)
        *absolute_error = cabs(left[index] - right[index]);
      if(worst_index != NULL)
        *worst_index = index;
    }
  }
  return worst;
}

int vie_cuboid_n(const vie_cuboid *testing, const vie_cuboid *source,
                 double k0, const vie_kernel_options *requested,
                 complex double result[3][3], vie_kernel_report *report)
{
  static const int orders[] = {4, 6, 8, 12, 16, 24, 32, 40, 48, 56, 64};
  vie_kernel_options defaults;
  const vie_kernel_options *options = requested;
  vie_relation relation;
  complex double previous[3][3], current[3][3];
  int have_previous = 0;
  int index;
  if(!cuboid_valid(testing) || !cuboid_valid(source) || !isfinite(k0) ||
     (k0 < 0.0) || (result == NULL))
    return 0;
  if((k0 * maximum3(testing->size[0], testing->size[1], testing->size[2]) >
      1.0 + VIE_GEOMETRY_TOLERANCE * DBL_EPSILON) ||
     (k0 * maximum3(source->size[0], source->size[1], source->size[2]) >
      1.0 + VIE_GEOMETRY_TOLERANCE * DBL_EPSILON))
    return 0;
  if(interiors_overlap(testing, source))
    return 0;
  if(options == NULL)
  {
    vie_kernel_default_options(&defaults);
    options = &defaults;
  }
  if(!isfinite(options->relative_tolerance) ||
     !isfinite(options->absolute_tolerance) ||
     (options->relative_tolerance <= 0.0) ||
     (options->absolute_tolerance <= 0.0) ||
     (options->maximum_order < 4) ||
     (options->maximum_order > VIE_MAX_ORDER) ||
     (options->far_order < 2) ||
     (options->far_order > options->maximum_order))
    return 0;
  relation = classify_relation(testing, source);
  if(report != NULL)
  {
    report->relation = relation;
    report->quadrature_order = 0;
    report->estimated_error = INFINITY;
    report->estimated_absolute_error = INFINITY;
    report->worst_row = -1;
    report->worst_column = -1;
  }
  if(relation == VIE_RELATION_FAR)
  {
    if(!matrix_at_order(testing, source, k0, options->far_order, result))
      return 0;
    if(report != NULL)
    {
      report->quadrature_order = options->far_order;
      report->estimated_error = NAN;
      report->estimated_absolute_error = NAN;
    }
    return 1;
  }
  for(index = 0; index < (int)(sizeof(orders) / sizeof(orders[0])); index++)
  {
    double error;
    int order = orders[index];
    if(order > options->maximum_order)
      break;
    if(!matrix_at_order(testing, source, k0, order, current))
      return 0;
    if(have_previous)
    {
      double absolute_error = 0.0;
      int worst_index = -1;
      error = matrix_difference(&previous[0][0], &current[0][0],
                                options->absolute_tolerance,
                                options->relative_tolerance,
                                &absolute_error, &worst_index);
      if(report != NULL)
      {
        report->quadrature_order = order;
        report->estimated_error = error;
        report->estimated_absolute_error = absolute_error;
        report->worst_row = worst_index / 3;
        report->worst_column = worst_index % 3;
      }
      if(error <= 1.0)
      {
        memcpy(result, current, sizeof(current));
        return 1;
      }
    }
    memcpy(previous, current, sizeof(previous));
    have_previous = 1;
  }
  return 0;
}
