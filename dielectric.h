#ifndef DIELECTRIC_H
#define DIELECTRIC_H 1

#include "nec2c.h"
#include "vie_kernel.h"

#include <stddef.h>
#include <float.h>
#include <stdio.h>

#define DIELECTRIC_RESOURCE_LINE_MAX 1024
#define DIELECTRIC_DEFAULT_MAX_KH 0.5
#define DIELECTRIC_DEFAULT_MAX_CELLS 1000
#define DIELECTRIC_DEFAULT_MAX_ASPECT 2.0

typedef struct
{
  int material_id, conductivity_kind;
  double frequency_mhz, epsilon_r, loss;
  size_t source_line;
} dielectric_sample;

typedef struct
{
  int id;
  double epsilon_r;
  size_t source_line;
  double epsilon_loss; /* positive imaginary-part magnitude, exp(+i omega t) */
  double conductivity; /* S/m, nonnegative */
  int tabulated;
  const dielectric_sample *samples;
  size_t sample_count;
} dielectric_material;

/* Piecewise-linear material data, exp(+i omega t). No extrapolation. */
static inline complex double dielectric_material_epsilon(const dielectric_material *m, double k0)
{
  double er=m->epsilon_r, ei=m->epsilon_loss, sigma=m->conductivity;
  if(m->tabulated)
  {
    double f=k0*299792458.0/(2*3.14159265358979323846*1e6), t;
    size_t j=1;
    const dielectric_sample *a,*b;
    if(m->sample_count<2 || !isfinite(f)) return NAN+I*NAN;
    if(f<m->samples[0].frequency_mhz-16*DBL_EPSILON*fmax(1,fabs(f)) ||
       f>m->samples[m->sample_count-1].frequency_mhz+16*DBL_EPSILON*fmax(1,fabs(f))) return NAN+I*NAN;
    f=fmax(m->samples[0].frequency_mhz,fmin(f,m->samples[m->sample_count-1].frequency_mhz));
    while(j+1<m->sample_count && f>m->samples[j].frequency_mhz) ++j;
    a=&m->samples[j-1]; b=&m->samples[j];
    t=(f-a->frequency_mhz)/(b->frequency_mhz-a->frequency_mhz);
    er=(1-t)*a->epsilon_r+t*b->epsilon_r;
    ei=(1-t)*a->loss+t*b->loss; sigma=0;
    if(a->conductivity_kind) {sigma=ei;ei=0;}
  }
  return er-I*(ei+(sigma==0 ? 0 : sigma/(k0*299792458.0*8.8541878128e-12)));
}

/* Conservative table-wide bound: each component is bounded independently.
   For linear sigma(f), sigma/f is monotone between consecutive positive knots. */
static inline double dielectric_material_mesh_index(const dielectric_material *m,double k0)
{
  size_t j;double er=0,ei=0;
  if(!m->tabulated) return sqrt(cabs(dielectric_material_epsilon(m,k0)));
  for(j=0;j<m->sample_count;++j)
  {
    const dielectric_sample *a=&m->samples[j];
    er=fmax(er,a->epsilon_r);
    ei=fmax(ei,a->conductivity_kind ? a->loss/(2*3.14159265358979323846*a->frequency_mhz*1e6*8.8541878128e-12):a->loss);
  }
  return sqrt(hypot(er,ei));
}

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
  dielectric_sample *samples;
  size_t sample_count;
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
  int near_enabled;
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
/* Validate the supported planar-ground deck subset before any wavelength normalization. */
int dielectric_ground_validate(const dielectric_model *model,
                               const dielectric_mesh_result *mesh, int iterative,
                               const dielectric_plane_wave *wave,
                               const dielectric_pattern_request *pattern,
                               dielectric_diagnostic *diagnostic);
int dielectric_solve_plane_wave(const dielectric_model *model,
                                const dielectric_mesh_result *mesh,
                                double frequency_mhz, int iterative,
                                const dielectric_plane_wave *wave,
                                const dielectric_pattern_request *pattern,
                                FILE *stream,
                                dielectric_diagnostic *diagnostic);
int dielectric_print_near_fields(const dielectric_mesh_result *mesh,
                                  const complex double *current, double k0,
                                  int conductor, FILE *stream);
/* SI scattered far-field amplitude: E(r) = amplitude * exp(-i k r) / r.
 * Uses the currently validated NEC free-space, PEC or Sommerfeld ground background. */
void dielectric_radiation_amplitude(const dielectric_mesh_result *mesh,
                                    const complex double *current, double k0,
                                    const double direction[3],
                                    complex double amplitude[3]);
int dielectric_solve_conductor_coupled(const dielectric_model *model,
                                       const dielectric_mesh_result *mesh,
                                       double frequency_mhz, double rkh,
                                       int extended_kernel, int iterative,
                                       const dielectric_plane_wave *wave,
                                       const dielectric_pattern_request *pattern,
                                       FILE *stream,
                                       dielectric_diagnostic *diagnostic);

int dielectric_ground_power(const dielectric_mesh_result *mesh,const complex double *current,
                             double k0,int conductor,double *power,double *delta);
/* NEC background adapter, implemented in dielectric_ground.c. */
int dielectric_background_n(const vie_cuboid *testing,const vie_cuboid *source,
                           double k0,const vie_kernel_options *options,
                           complex double n[3][3]);
void dielectric_plane_wave_vectors(const dielectric_plane_wave *wave,double direction[3],complex double field[3]);
double dielectric_ground_extinction(const dielectric_mesh_result *mesh,const complex double *current,
                                    double k0,const dielectric_plane_wave *wave,int conductor);
/* Shared planar-background operations; geometry is in metres. */
complex double dielectric_ground_epsilon(double frequency_mhz);
void dielectric_ground_initialize(double frequency_mhz);
void dielectric_ground_reflect(const double downward[3], const complex double field[3],
                              double k0, complex double reflected[3]);
int dielectric_validate_sweep(const dielectric_model *model,int mode,int count,
    double first,double step,dielectric_diagnostic *diagnostic);
/* Joule/dielectric absorption from cell-averaged polarization currents. */
double dielectric_absorbed_power(const dielectric_model *model,
    const dielectric_mesh_result *mesh,const complex double *current,double k0);
#endif
