#ifndef RESULT_JSON_H
#define RESULT_JSON_H 1

#include "dielectric.h"

#include <complex.h>
#include <stddef.h>
#include <stdio.h>

typedef struct
{
  const char *requested;
  const char *used;
  int iterative_fields_defined;
  double relative_tolerance;
  size_t maximum_iterations;
  size_t iterations;
  double reported_residual;
  double true_residual;
  int converged;
  size_t accounted_storage_bytes;
  double elapsed_seconds;
} result_solver;

typedef struct
{
  const char *basis;
  int accepted_defined;
  double accepted_w;
  int air_radiated_defined;
  double air_radiated_w;
  int lower_flux_defined;
  double lower_flux_w;
  int dielectric_absorbed_defined;
  double dielectric_absorbed_w;
  int structure_loss_defined;
  double structure_loss_w;
  int network_loss_defined;
  double network_loss_w;
  int defect_defined;
  double defect_w;
  double defect_relative;
  size_t theta_count;
  size_t phi_count;
} result_power;

/* CLI/session lifecycle. All recording calls are no-ops when inactive. */
int result_json_write_capabilities(FILE *stream);
int result_json_open(const char *destination, const char *deck_path);
int result_json_active(void);
void result_json_set_validation_only(int enabled);
int result_json_validation_only(void);
void result_json_job_begin(void);
void result_json_job_end(const dielectric_model *model);
int result_json_finish(const char *status);
void result_json_abandon(void);

/* Job-local interpreted states. */
int result_json_add_mesh_state(const dielectric_model *model,
    const dielectric_mesh_result *mesh, int frequency_mode,
    int frequency_count, double first_mhz, double step);
void result_json_invalidate_mesh(void);
int result_json_current_mesh_state(void);
int result_json_ground_state(void);
void result_json_invalidate_ground(void);

/* A request index is allocated once per source control-card request. */
int result_json_request(const char *kind, const char *parameters_json);
int result_json_request_xq(int mode);
int result_json_request_rp(int mode, int theta_count, int phi_count,
    double theta_start, double phi_start, double theta_step,
    double phi_step, double range_m);
int result_json_request_near(const char *kind, int coordinate_mode,
    int x_count, int y_count, int z_count, const double start[3],
    const double step[3]);
int result_json_current_request(void);
const char *result_json_current_excitation_kind(void);

/* Execution lifecycle. A new execution finalizes the preceding successful one. */
void result_json_execution_begin(double frequency_mhz, const char *excitation_kind,
    const dielectric_plane_wave *wave);
int result_json_execution_end(int success);
void result_json_set_solver(const result_solver *solver);

void result_json_add_feed(size_t source_index, int tag, size_t segment,
    complex double voltage, complex double current, complex double impedance,
    int accepted_defined, double accepted_w);
void result_json_add_wire(int tag, size_t segment, const double center_m[3],
    double length_m, complex double current_a);
void result_json_add_patch(size_t patch_id, const double center_m[3],
    double area_m2, const double tangent1[3], const double tangent2[3],
    const complex double current_a_per_m[3]);
void result_json_add_volume(size_t cell_id, size_t box_id, int region_id,
    int material_id, const double center_m[3], const double size_m[3],
    const complex double current_a_per_m2[3]);
void result_json_add_material_evaluations(const dielectric_model *model,
    double frequency_mhz);
void result_json_set_power(const result_power *power);

void result_json_far_begin(size_t request_index, const char *kind);
void result_json_far_sample(double theta_deg, double phi_deg,
    complex double r_e_theta, complex double r_e_phi);
void result_json_near_begin(size_t request_index, const char *field,
    const char *kind, const char *units);
void result_json_near_sample(const double position_m[3],
    const complex double value[3], int valid);

void result_json_omit(int request_index, const char *kind,
    const char *code, const char *message);
void result_json_diagnostic(const char *severity, const char *code,
    const char *message, const char *context_json);
void result_json_dielectric_diagnostic(const dielectric_diagnostic *diagnostic,
    const char *severity);

#endif
