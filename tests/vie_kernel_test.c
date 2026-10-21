#include "vie_kernel.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

static void fail(const char *message)
{
  fprintf(stderr, "FAIL: %s\n", message);
  failures++;
}

static double matrix_max_difference(const complex double left[3][3],
                                    const complex double right[3][3])
{
  double maximum = 0.0;
  int row, column;
  for(row = 0; row < 3; row++)
    for(column = 0; column < 3; column++)
      maximum = fmax(maximum, cabs(left[row][column] - right[row][column]));
  return maximum;
}

static void require_kernel(const char *name, const vie_cuboid *testing,
                           const vie_cuboid *source, double k0,
                           complex double result[3][3],
                           vie_kernel_report *report)
{
  vie_kernel_options options;
  vie_kernel_default_options(&options);
  if(!vie_cuboid_n(testing, source, k0, &options, result, report))
  {
    fprintf(stderr, "FAIL: %s did not converge (order %d normalized delta %.6g)\n",
            name, report->quadrature_order, report->estimated_error);
    failures++;
  }
}

static void test_static_cube_self(void)
{
  vie_cuboid cube = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
  vie_cuboid rectangular = {{0.2, -0.4, 0.1}, {0.2, 0.5, 0.9}};
  vie_kernel_report report;
  complex double value[3][3];
  int row, column;
  require_kernel("static cube self", &cube, &cube, 0.0, value, &report);
  if(failures)
    return;
  if(report.relation != VIE_RELATION_SELF)
    fail("cube self relation classification");
  for(row = 0; row < 3; row++)
    for(column = 0; column < 3; column++)
    {
      double expected = row == column ? 2.0 / 3.0 : 0.0;
      if(cabs(value[row][column] - expected) > 2.0e-9)
      {
        fprintf(stderr, "FAIL: cube self [%d,%d] = %.17g%+.17gi\n",
                row, column, creal(value[row][column]), cimag(value[row][column]));
        failures++;
      }
    }
  require_kernel("static rectangular self", &rectangular, &rectangular,
                 0.0, value, &report);
  if(cabs(value[0][0] + value[1][1] + value[2][2] - 2.0) > 3.0e-9)
    fail("static rectangular self trace identity");
  for(row = 0; row < 3; row++)
    for(column = 0; column < 3; column++)
      if((row != column) && (cabs(value[row][column]) > 2.0e-10))
        fail("static rectangular self off-diagonal symmetry");
}

static void test_translation_and_reciprocity(void)
{
  vie_cuboid a = {{0.1, -0.2, 0.3}, {0.3, 0.2, 0.1}};
  vie_cuboid b = {{0.65, 0.15, -0.05}, {0.2, 0.25, 0.15}};
  vie_cuboid moved_a = a;
  vie_cuboid moved_b = b;
  complex double ab[3][3], ba[3][3], moved[3][3];
  vie_kernel_report report;
  double volume_a = a.size[0] * a.size[1] * a.size[2];
  double volume_b = b.size[0] * b.size[1] * b.size[2];
  int axis, row, column;
  for(axis = 0; axis < 3; axis++)
  {
    moved_a.center[axis] += 17.0 - 3.0 * axis;
    moved_b.center[axis] += 17.0 - 3.0 * axis;
  }
  require_kernel("a-b", &a, &b, 2.3, ab, &report);
  require_kernel("b-a", &b, &a, 2.3, ba, &report);
  require_kernel("translated a-b", &moved_a, &moved_b, 2.3, moved, &report);
  if(matrix_max_difference(ab, moved) > 3.0e-12)
    fail("translation invariance");
  for(row = 0; row < 3; row++)
    for(column = 0; column < 3; column++)
      if(cabs(volume_a * ab[row][column] -
               volume_b * ba[column][row]) > 2.0e-11)
        fail("volume-scaled reciprocity");
}

static void test_axis_permutation(void)
{
  const int permutation[3] = {2, 0, 1};
  vie_cuboid a = {{0.0, 0.1, 0.2}, {0.2, 0.3, 0.4}};
  vie_cuboid b = {{0.7, -0.2, 0.4}, {0.1, 0.25, 0.35}};
  vie_cuboid pa, pb;
  complex double original[3][3], permuted[3][3];
  vie_kernel_report report;
  int axis, row, column;
  for(axis = 0; axis < 3; axis++)
  {
    pa.center[axis] = a.center[permutation[axis]];
    pa.size[axis] = a.size[permutation[axis]];
    pb.center[axis] = b.center[permutation[axis]];
    pb.size[axis] = b.size[permutation[axis]];
  }
  require_kernel("axis original", &a, &b, 1.7, original, &report);
  require_kernel("axis permuted", &pa, &pb, 1.7, permuted, &report);
  for(row = 0; row < 3; row++)
    for(column = 0; column < 3; column++)
      if(cabs(permuted[row][column] -
               original[permutation[row]][permutation[column]]) > 2.0e-11)
        fail("Cartesian-axis permutation invariance");
}

static void test_static_limit(void)
{
  vie_cuboid a = {{0.0, 0.0, 0.0}, {0.2, 0.3, 0.4}};
  vie_cuboid b = {{0.2, 0.0, 0.0}, {0.2, 0.3, 0.4}};
  complex double zero[3][3], small[3][3];
  vie_kernel_report report;
  require_kernel("static face", &a, &b, 0.0, zero, &report);
  if(report.relation != VIE_RELATION_FACE)
    fail("face adjacency classification");
  require_kernel("small-k face", &a, &b, 1.0e-8, small, &report);
  if(matrix_max_difference(zero, small) > 2.0e-9)
    fail("stable static limit");
}

static void test_reference_data(void)
{
  char path[4096];
  char line[4096];
  const char *source_root = getenv("TEST_SRCDIR");
  FILE *stream;
  int cases = 0;
  if(source_root == NULL)
    source_root = ".";
  if(snprintf(path, sizeof(path), "%s/tests/vie-reference.csv", source_root) >=
     (int)sizeof(path))
  {
    fail("reference-data path is too long");
    return;
  }
  stream = fopen(path, "r");
  if(stream == NULL)
  {
    fail("cannot open reference data");
    return;
  }
  while(fgets(line, sizeof(line), stream) != NULL)
  {
    char *fields[32];
    char *token;
    int field_count = 0;
    double number[31];
    vie_cuboid testing, source;
    vie_kernel_report report;
    complex double actual[3][3];
    int index, row, column;
    if(line[0] == '#')
      continue;
    token = strtok(line, ",\r\n");
    while((token != NULL) && (field_count < 32))
    {
      fields[field_count++] = token;
      token = strtok(NULL, ",\r\n");
    }
    if(field_count != 32)
    {
      fail("malformed reference-data row");
      continue;
    }
    for(index = 1; index < 32; index++)
    {
      char *end;
      number[index - 1] = strtod(fields[index], &end);
      if((end == fields[index]) || (*end != '\0') || !isfinite(number[index - 1]))
      {
        fail("non-numeric reference-data field");
        break;
      }
    }
    if(index != 32)
      continue;
    for(index = 0; index < 3; index++)
    {
      testing.center[index] = number[1 + index];
      testing.size[index] = number[4 + index];
      source.center[index] = number[7 + index];
      source.size[index] = number[10 + index];
    }
    require_kernel(fields[0], &testing, &source, number[0], actual, &report);
    for(row = 0; row < 3; row++)
      for(column = 0; column < 3; column++)
      {
        int offset = 13 + 2 * (3 * row + column);
        complex double expected = number[offset] + I * number[offset + 1];
        double tolerance = 3.0e-10 + 3.0e-8 * cabs(expected);
        if(cabs(actual[row][column] - expected) > tolerance)
        {
          fprintf(stderr,
                  "FAIL: oracle %s [%d,%d] expected %.17g%+.17gi got %.17g%+.17gi\n",
                  fields[0], row, column, creal(expected), cimag(expected),
                  creal(actual[row][column]), cimag(actual[row][column]));
          failures++;
        }
      }
    cases++;
  }
  if(ferror(stream))
    fail("error reading reference data");
  fclose(stream);
  if(cases != 6)
    fail("reference-data case count");
}

static void test_contact_continuity_and_scaling(void)
{
  vie_cuboid a = {{0.0, 0.0, 0.0}, {0.2, 0.3, 0.4}};
  vie_cuboid contact = {{0.2, 0.0, 0.0}, {0.2, 0.3, 0.4}};
  vie_cuboid gap_large = contact;
  vie_cuboid gap_small = contact;
  vie_cuboid scaled_a = a;
  vie_cuboid scaled_contact = contact;
  complex double at_contact[3][3], at_large[3][3], at_small[3][3], scaled[3][3];
  vie_kernel_report report;
  double large_error, small_error;
  int axis;
  gap_large.center[0] += 1.0e-3;
  gap_small.center[0] += 1.0e-5;
  require_kernel("face contact", &a, &contact, 1.1, at_contact, &report);
  require_kernel("face gap 1e-3", &a, &gap_large, 1.1, at_large, &report);
  require_kernel("face gap 1e-5", &a, &gap_small, 1.1, at_small, &report);
  large_error = matrix_max_difference(at_contact, at_large);
  small_error = matrix_max_difference(at_contact, at_small);
  if(!(small_error < 0.2 * large_error))
    fail("continuity on approach to face contact");
  for(axis = 0; axis < 3; axis++)
  {
    scaled_a.center[axis] *= 7.0;
    scaled_a.size[axis] *= 7.0;
    scaled_contact.center[axis] *= 7.0;
    scaled_contact.size[axis] *= 7.0;
  }
  require_kernel("scaled dynamic face", &scaled_a, &scaled_contact,
                 1.1 / 7.0, scaled, &report);
  if(matrix_max_difference(at_contact, scaled) > 2.0e-10)
    fail("electrical-scale invariance");
}

static void test_translated_contact_roundoff(void)
{
  vie_cuboid a={{0,0,0},{.025,.025,.025}},b=a,t=a,u=a;
  complex double ref[3][3],translated[3][3];vie_kernel_report report;
  int mask,d;
  const double base[3]={.0875,-.0125,.3575};
  for(mask=1;mask<8;mask++)
  {
    for(d=0;d<3;d++) {
      b.center[d]=(mask&(1<<d)) ? .025:0;
      t.center[d]=base[d];u.center[d]=base[d]+b.center[d];
    }
    require_kernel("canonical touching cubes",&a,&b,6.2875,ref,&report);
    require_kernel("translated touching cubes",&t,&u,6.2875,translated,&report);
    if(matrix_max_difference(ref,translated)>2e-10)
      fail("translated contact quadrature consistency");
    for(d=0;d<3;d++) if(mask&(1<<d)) u.center[d]=nextafter(u.center[d],INFINITY);
    require_kernel("one-ulp touching cubes",&t,&u,6.2875,translated,&report);
    if(matrix_max_difference(ref,translated)>2e-10)
      fail("one-ulp contact quadrature consistency");
  }
}

static void test_resolved_gap_at_small_scale(void)
{
  vie_cuboid a={{0,0,0},{.025,.025,.025}},b=a,t=a,u=a;
  complex double ref[3][3],scaled[3][3];vie_kernel_report report;int d;
  const double factor=1e-10;
  b.center[0]=.025+2.5e-5;
  require_kernel("resolved positive gap",&a,&b,6.2875,ref,&report);
  for(d=0;d<3;d++) {
    t.center[d]=a.center[d]*factor;u.center[d]=b.center[d]*factor;
    t.size[d]=a.size[d]*factor;u.size[d]=b.size[d]*factor;
  }
  require_kernel("small-scale positive gap",&t,&u,6.2875/factor,scaled,&report);
  if(matrix_max_difference(ref,scaled)>2e-10)
    fail("resolved positive gap was changed by electrical scaling");
}

static void test_failure_behavior(void)
{
  vie_cuboid cube = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
  vie_cuboid overlap = {{0.25, 0.0, 0.0}, {1.0, 1.0, 1.0}};
  complex double value[3][3];
  vie_kernel_options options;
  vie_kernel_default_options(&options);
  if(vie_cuboid_n(&cube, &cube, -1.0, &options, value, NULL))
    fail("negative wavenumber was accepted");
  if(vie_cuboid_n(&cube, &cube, 1.0001, &options, value, NULL))
    fail("electrically oversized cuboid was accepted");
  if(vie_cuboid_n(&cube, &overlap, 0.5, &options, value, NULL))
    fail("overlapping distinct cuboids were accepted");
  options.maximum_order = 3;
  if(vie_cuboid_n(&cube, &cube, 1.0, &options, value, NULL))
    fail("invalid maximum quadrature order was accepted");
}

int main(void)
{
  test_static_cube_self();
  test_translation_and_reciprocity();
  test_axis_permutation();
  test_static_limit();
  test_reference_data();
  test_contact_continuity_and_scaling();
  test_failure_behavior();
  test_translated_contact_roundoff();
  test_resolved_gap_at_small_scale();
  if(failures != 0)
    return 1;
  printf("PASS: M4 cuboid kernel analytical, invariance, reciprocity, and static-limit checks\n");
  return 0;
}
