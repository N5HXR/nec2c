#ifndef DIELECTRIC_H
#define DIELECTRIC_H 1

#include "nec2c.h"

#include <stddef.h>
#include <stdio.h>

#define DIELECTRIC_RESOURCE_LINE_MAX 1024
#define DIELECTRIC_DEFAULT_MAX_KH 0.5
#define DIELECTRIC_DEFAULT_MAX_CELLS 1000
#define DIELECTRIC_DEFAULT_MAX_ASPECT 2.0

typedef struct
{
  int id;
  double epsilon_r;
  size_t source_line;
} dielectric_material;

typedef struct
{
  int id;
  int material_id;
  size_t source_line;
} dielectric_region;

typedef struct
{
  int region_id;
  double minimum[3];
  double maximum[3];
  size_t source_line;
} dielectric_box;

typedef struct
{
  int loaded;
  int format_version;
  int path_is_absolute;
  int mesh_record_present;
  char *resource_path;
  dielectric_material *materials;
  size_t material_count;
  dielectric_region *regions;
  size_t region_count;
  dielectric_box *boxes;
  size_t box_count;
  double max_kh;
  size_t max_cells;
  double max_aspect;
} dielectric_model;

typedef struct
{
  char resource[NEC2C_PATH_MAX * 4 + 1];
  char record[32];
  char value[128];
  char message[256];
  size_t line;
} dielectric_diagnostic;

typedef struct
{
  int region_id;
  int material_id;
  double center[3];
  double size[3];
  double q;
  double aspect;
} dielectric_cell;

typedef struct
{
  int region_id;
  int material_id;
  size_t divisions[3];
  size_t cell_offset;
  size_t cell_count;
  double physical_volume;
  double numerical_volume;
} dielectric_box_mesh;

typedef struct
{
  int region_id;
  int material_id;
  double diameter;
  double electrical_diameter;
} dielectric_region_mesh;

typedef struct
{
  int built;
  int accuracy_warning;
  double frequency_mhz;
  dielectric_region_mesh *regions;
  size_t region_count;
  dielectric_box_mesh *boxes;
  size_t box_count;
  dielectric_cell *cells;
  size_t cell_count;
  size_t vacuum_box_count;
  size_t feature_warning_count;
  double declared_volume;
  double active_volume;
  double vacuum_volume;
  double numerical_volume;
  double volume_relative_error;
  double q_min;
  double q_median;
  double q_p95;
  double q_max;
  double aspect_max;
  int has_conductors;
  double minimum_conductor_gap;
  double required_conductor_gap;
  size_t conductor_unknowns;
  size_t volume_unknowns;
  size_t total_unknowns;
  size_t matrix_bytes;
  size_t workspace_bytes;
} dielectric_mesh_result;

typedef struct
{
  int type;
  double theta_degrees;
  double phi_degrees;
  double eta_degrees;
  double axial_ratio;
} dielectric_plane_wave;

typedef struct
{
  int enabled;
  size_t theta_count;
  size_t phi_count;
  double theta_start_degrees;
  double phi_start_degrees;
  double theta_step_degrees;
  double phi_step_degrees;
} dielectric_pattern_request;

void dielectric_model_init(dielectric_model *model);
void dielectric_model_free(dielectric_model *model);
int dielectric_load(dielectric_model *model, const char *deck_path,
                    const char *df_field, dielectric_diagnostic *diagnostic);
void dielectric_print_diagnostic(FILE *stream,
                                 const dielectric_diagnostic *diagnostic);
void dielectric_print_summary(FILE *stream, const dielectric_model *model);
void dielectric_mesh_init(dielectric_mesh_result *mesh);
void dielectric_mesh_free(dielectric_mesh_result *mesh);
int dielectric_mesh_build(const dielectric_model *model, double frequency_mhz,
                          const data_t *geometry,
                          dielectric_mesh_result *mesh,
                          dielectric_diagnostic *diagnostic);
void dielectric_print_mesh_report(FILE *stream, const dielectric_mesh_result *mesh);
int dielectric_solve_plane_wave(const dielectric_model *model,
                                const dielectric_mesh_result *mesh,
                                double frequency_mhz,
                                const dielectric_plane_wave *wave,
                                const dielectric_pattern_request *pattern,
                                FILE *stream,
                                dielectric_diagnostic *diagnostic);
int dielectric_solve_conductor_coupled(const dielectric_model *model,
                                       const dielectric_mesh_result *mesh,
                                       double frequency_mhz, double rkh,
                                       int extended_kernel,
                                       const dielectric_plane_wave *wave,
                                       FILE *stream,
                                       dielectric_diagnostic *diagnostic);

#endif
