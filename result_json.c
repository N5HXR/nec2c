#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "result_json.h"
#include "shared.h"
#include "dielectric_fft.h"

#include <float.h>
#include <stdint.h>

typedef struct
{
  FILE *file;
  size_t count;
} json_spool;

typedef struct
{
  int open;
  int status_ok;
  size_t index;
  double frequency_mhz;
  int mesh_state_id;
  int ground_state_id;
  int request_index;
  char trigger_kind[8];
  char trigger_parameters[512];
  char excitation_kind[32];
  dielectric_plane_wave wave;
  int has_wave;
  int solver_set;
  result_solver solver;
  int power_set;
  result_power power;
  json_spool feeds, wires, patches, volumes, materials, far_fields, near_fields;
  int far_open, near_open;
  size_t far_samples, near_samples;
} result_execution_state;

typedef struct
{
  int open;
  size_t index;
  size_t next_execution;
  size_t next_request;
  size_t next_mesh;
  size_t next_ground;
  int current_mesh;
  int current_request;
  char current_trigger[8];
  char current_parameters[512];
  char last_ground[512];
  int last_ground_id;
  json_spool meshes, grounds, executions, omissions;
  result_execution_state execution;
} result_job_state;

typedef struct
{
  int active;
  int validation_only;
  int failed;
  char destination[NEC2C_PATH_MAX * 4 + 1];
  char temporary[NEC2C_PATH_MAX * 4 + 64];
  FILE *document;
  size_t job_count;
  json_spool diagnostics;
  result_job_state job;
} result_session;

static result_session session;

static int finite_number(double value) { return isfinite(value); }

static void json_string(FILE *f, const char *s)
{
  const unsigned char *p = (const unsigned char *)(s ? s : "");
  fputc('"', f);
  while(*p)
  {
    unsigned char c = *p++;
    switch(c)
    {
      case '"': fputs("\\\"", f); break;
      case '\\': fputs("\\\\", f); break;
      case '\b': fputs("\\b", f); break;
      case '\f': fputs("\\f", f); break;
      case '\n': fputs("\\n", f); break;
      case '\r': fputs("\\r", f); break;
      case '\t': fputs("\\t", f); break;
      default:
        if(c < 0x20) fprintf(f, "\\u%04x", (unsigned)c);
        else fputc((int)c, f);
    }
  }
  fputc('"', f);
}

/* JSON always requires '.', regardless of the process locale. Replace a
   locale decimal separator in the small local buffer without changing the
   process-global locale used by classic NEC parsing and reports. */
static int json_number(FILE *f, double value)
{
  char buffer[96], *comma;
  if(!finite_number(value))
  {
    session.failed = 1;
    fputs("null", f);
    return 0;
  }
  snprintf(buffer, sizeof(buffer), "%.17g", value);
  comma = strchr(buffer, ',');
  if(comma) *comma = '.';
  fputs(buffer, f);
  return 1;
}

static void number_text(char *buffer, size_t size, double value)
{
  char *comma;
  snprintf(buffer, size, "%.17g", value);
  comma = strchr(buffer, ',');
  if(comma) *comma = '.';
}

static void json_complex(FILE *f, complex double value)
{
  fputs("{\"re\":", f); json_number(f, creal(value));
  fputs(",\"im\":", f); json_number(f, cimag(value)); fputc('}', f);
}

static int spool_open(json_spool *spool)
{
  if(spool->file) return 1;
  spool->file = tmpfile();
  spool->count = 0;
  return spool->file != NULL;
}

static FILE *spool_row(json_spool *spool)
{
  if(!spool_open(spool)) return NULL;
  if(spool->count++) fputc(',', spool->file);
  return spool->file;
}

static void spool_copy(FILE *out, json_spool *spool)
{
  char buffer[8192];
  size_t count;
  if(!spool->file) return;
  fflush(spool->file);
  rewind(spool->file);
  while((count = fread(buffer, 1, sizeof(buffer), spool->file)) != 0)
    fwrite(buffer, 1, count, out);
}

static void spool_close(json_spool *spool)
{
  if(spool->file) fclose(spool->file);
  memset(spool, 0, sizeof(*spool));
}

static const char *base_name(const char *path)
{
  const char *a = strrchr(path ? path : "", '/');
  const char *b = strrchr(path ? path : "", '\\');
  const char *p = a;
  if(!p || (b && b > p)) p = b;
  return p ? p + 1 : (path ? path : "");
}

int result_json_write_capabilities(FILE *stream)
{
  if(!stream) return 0;
  fprintf(stream,
    "{\"program\":\"nec2c-n5hxr\",\"program_version\":");
  json_string(stream, PACKAGE_VERSION);
  fprintf(stream,
    ",\"schema\":\"nec2c-capabilities\",\"schema_version\":1,"
    "\"features\":{\"dielectric_feature_revision\":1,"
    "\"dielectric_resource_versions\":[1,2,3],\"result_json_versions\":[1],"
    "\"validate_only\":true,\"legacy_result_output\":true,"
    "\"solvers\":[\"dense\",\"gmres\",\"fft\"],"
    "\"grounds\":[\"free-space\",\"pec\",\"sommerfeld-norton\"]},"
    "\"limits\":{\"nec_card_characters\":132,\"default_max_cells\":%d,"
    "\"fft_accounted_storage_bytes\":%lu}}\n",
    DIELECTRIC_DEFAULT_MAX_CELLS, (unsigned long)DIELECTRIC_FFT_BUDGET);
  return !ferror(stream);
}

int result_json_open(const char *destination, const char *deck_path)
{
  unsigned long nonce = (unsigned long)time(NULL) ^ (unsigned long)getpid();
  if(!destination || !*destination || !strcmp(destination, "-")) return 0;
  memset(&session, 0, sizeof(session));
  if(strlen(destination) >= sizeof(session.destination)) return 0;
  strcpy(session.destination, destination);
  if(snprintf(session.temporary, sizeof(session.temporary), "%s.tmp.%lu",
      destination, nonce) >= (int)sizeof(session.temporary)) return 0;
  session.document = fopen(session.temporary, "wb");
  if(!session.document) return 0;
  if(!spool_open(&session.diagnostics)) { fclose(session.document); return 0; }
  session.active = 1;
  fputs("{\"schema\":\"nec2c-results\",\"schema_version\":1,\"producer\":{"
        "\"program\":\"nec2c-n5hxr\",\"program_version\":", session.document);
  json_string(session.document, PACKAGE_VERSION);
  fputs(",\"dielectric_feature_revision\":1},\"input\":{\"deck\":",
        session.document);
  json_string(session.document, base_name(deck_path));
  fputs("},\"jobs\":[", session.document);
  return !ferror(session.document);
}

int result_json_active(void) { return session.active; }
void result_json_set_validation_only(int enabled) { session.validation_only = enabled != 0; }
int result_json_validation_only(void) { return session.validation_only; }

static void execution_clear(result_execution_state *e)
{
  spool_close(&e->feeds); spool_close(&e->wires); spool_close(&e->patches);
  spool_close(&e->volumes); spool_close(&e->materials); spool_close(&e->far_fields);
  spool_close(&e->near_fields);
  memset(e, 0, sizeof(*e));
}

static void job_clear(result_job_state *job)
{
  execution_clear(&job->execution);
  spool_close(&job->meshes); spool_close(&job->grounds);
  spool_close(&job->executions); spool_close(&job->omissions);
  memset(job, 0, sizeof(*job));
  job->current_mesh = -1; job->current_request = -1; job->last_ground_id = -1;
}

void result_json_job_begin(void)
{
  result_job_state *job;
  if(!session.active) return;
  if(session.job.open) result_json_job_end(NULL);
  job = &session.job;
  job_clear(job);
  job->open = 1;
  job->index = session.job_count;
  spool_open(&job->meshes); spool_open(&job->grounds);
  spool_open(&job->executions); spool_open(&job->omissions);
}

static void write_material(FILE *f, const dielectric_material *m)
{
  size_t i;
  fprintf(f, "{\"id\":%d,\"law\":", m->id);
  json_string(f, m->tabulated ? "table" :
      (m->conductivity != 0.0 ? "constant-sigma" :
       (m->epsilon_loss != 0.0 ? "constant-epsilon-loss" : "constant-er")));
  if(m->tabulated)
  {
    fputs(",\"samples\":[", f);
    for(i = 0; i < m->sample_count; ++i)
    {
      const dielectric_sample *s = &m->samples[i];
      if(i) fputc(',', f);
      fputs("{\"frequency_mhz\":", f); json_number(f, s->frequency_mhz);
      fputs(",\"er\":", f); json_number(f, s->epsilon_r);
      fputs(s->conductivity_kind ? ",\"sigma_s_per_m\":" : ",\"epsilon_loss\":", f);
      json_number(f, s->loss); fputc('}', f);
    }
    fputc(']', f);
  }
  else
  {
    fputs(",\"er\":", f); json_number(f, m->epsilon_r);
    if(m->conductivity != 0.0) { fputs(",\"sigma_s_per_m\":", f); json_number(f, m->conductivity); }
    if(m->epsilon_loss != 0.0) { fputs(",\"epsilon_loss\":", f); json_number(f, m->epsilon_loss); }
  }
  fputc('}', f);
}

void result_json_job_end(const dielectric_model *model)
{
  result_job_state *job = &session.job;
  FILE *f;
  size_t i;
  if(!session.active || !job->open) return;
  if(job->execution.open) result_json_execution_end(1);
  f = session.document;
  if(session.job_count++) fputc(',', f);
  fprintf(f, "{\"job_index\":%lu,\"dielectric\":",
      (unsigned long)job->index);
  if(!model || !model->loaded) fputs("null", f);
  else
  {
    fputs("{\"resource\":", f); json_string(f, base_name(model->resource_path));
    fprintf(f, ",\"resource_version\":%d,\"materials\":[", model->format_version);
    for(i = 0; i < model->material_count; ++i) { if(i) fputc(',', f); write_material(f, &model->materials[i]); }
    fputs("],\"regions\":[", f);
    for(i = 0; i < model->region_count; ++i)
    {
      if(i) fputc(',', f);
      fprintf(f, "{\"id\":%d,\"material_id\":%d}", model->regions[i].id,
          model->regions[i].material_id);
    }
    fputs("],\"boxes\":[", f);
    for(i = 0; i < model->box_count; ++i)
    {
      const dielectric_box *b = &model->boxes[i];
      if(i) fputc(',', f);
      fprintf(f, "{\"box_id\":%lu,\"region_id\":%d,\"bounds_m\":{\"min\":[",
          (unsigned long)(i + 1), b->region_id);
      json_number(f,b->minimum[0]);fputc(',',f);json_number(f,b->minimum[1]);fputc(',',f);json_number(f,b->minimum[2]);
      fputs("],\"max\":[",f);json_number(f,b->maximum[0]);fputc(',',f);json_number(f,b->maximum[1]);fputc(',',f);json_number(f,b->maximum[2]);fputs("]}}",f);
    }
    fputs("]}", f);
  }
  if(session.validation_only)
  {
    fputs(",\"validation\":",f);
    json_string(f,model&&model->loaded?"dielectric-preflight":"classic-input-only");
  }
  fputs(",\"mesh_states\":[", f); spool_copy(f, &job->meshes);
  fputs("],\"ground_states\":[", f); spool_copy(f, &job->grounds);
  fputs("],\"executions\":[", f); spool_copy(f, &job->executions);
  fputs("],\"omitted_requests\":[", f); spool_copy(f, &job->omissions);
  fputs("]}", f);
  job_clear(job);
}

static const dielectric_box *box_for_mesh(const dielectric_model *model, size_t index)
{ return index < model->box_count ? &model->boxes[index] : NULL; }

int result_json_add_mesh_state(const dielectric_model *model,
    const dielectric_mesh_result *mesh, int mode, int count,
    double first, double step)
{
  result_job_state *job = &session.job;
  FILE *f;
  size_t i;
  int id;
  if(!session.active || !job->open || !model || !mesh || !mesh->built) return -1;
  id = (int)job->next_mesh++;
  f = spool_row(&job->meshes); if(!f) return -1;
  fprintf(f, "{\"mesh_state_id\":%d,\"realized_for\":{\"frequency_mode\":", id);
  json_string(f, mode == 1 ? "multiplicative" : "linear");
  fprintf(f, ",\"frequency_count\":%d,\"first_frequency_mhz\":", count);json_number(f,first);
  fputs(",\"frequency_step_mhz\":",f);json_number(f,step);
  fputs(",\"design_frequency_mhz\":",f);json_number(f,mesh->frequency_mhz);
  fputs("},\"boxes\":[", f);
  for(i = 0; i < mesh->box_count; ++i)
  {
    const dielectric_box_mesh *bm = &mesh->boxes[i];
    const dielectric_box *b = box_for_mesh(model, i);
    double size[3] = {0,0,0};
    int p;
    if(i) fputc(',',f);
    if(b) for(p=0;p<3;++p) if(bm->divisions[p]) size[p]=(b->maximum[p]-b->minimum[p])/(double)bm->divisions[p];
    fprintf(f,"{\"box_id\":%lu,\"region_id\":%d,\"divisions\":[%lu,%lu,%lu],\"cell_size_m\":[",
      (unsigned long)(i+1),bm->region_id,(unsigned long)bm->divisions[0],(unsigned long)bm->divisions[1],(unsigned long)bm->divisions[2]);
    json_number(f,size[0]);fputc(',',f);json_number(f,size[1]);fputc(',',f);json_number(f,size[2]);
    fprintf(f,"],\"first_cell_id\":%lu,\"cell_count\":%lu}",(unsigned long)(bm->cell_count?bm->cell_offset+1:0),(unsigned long)bm->cell_count);
  }
  fprintf(f,"],\"mesh\":{\"max_kh\":");json_number(f,model->max_kh);
  fputs(",\"max_aspect\":",f);json_number(f,model->max_aspect);
  fprintf(f,",\"max_cells\":%lu,\"active_cells\":%lu,\"volume_unknowns\":%lu,\"total_unknowns\":%lu,\"projected_dense_matrix_bytes\":%lu}}",
    (unsigned long)model->max_cells,(unsigned long)mesh->cell_count,(unsigned long)mesh->volume_unknowns,
    (unsigned long)mesh->total_unknowns,(unsigned long)mesh->matrix_bytes);
  job->current_mesh = id;
  return id;
}

void result_json_invalidate_mesh(void) { if(session.active) session.job.current_mesh = -1; }
int result_json_current_mesh_state(void) { return session.active ? session.job.current_mesh : -1; }

static void ground_description(char *buffer, size_t size)
{
  char epsilon[96], conductivity[96];
  if(gnd.ksymp == 1) snprintf(buffer,size,"{\"kind\":\"free-space\"}");
  else if(gnd.iperf == 1) snprintf(buffer,size,"{\"kind\":\"pec\"}");
  else if(gnd.iperf == 2)
  {
    number_text(epsilon,sizeof(epsilon),save.epsr);
    number_text(conductivity,sizeof(conductivity),fabs(save.sig));
    if(save.sig<0)
      snprintf(buffer,size,"{\"kind\":\"sommerfeld-norton\",\"relative_permittivity\":%s,\"epsilon_loss\":%s}",epsilon,conductivity);
    else
      snprintf(buffer,size,"{\"kind\":\"sommerfeld-norton\",\"relative_permittivity\":%s,\"conductivity_s_per_m\":%s}",epsilon,conductivity);
  }
  else snprintf(buffer,size,"{\"kind\":\"classic-unsupported\",\"gn_type\":%d,\"radial_wire_count\":%d}",gnd.iperf,gnd.nradl);
}

int result_json_ground_state(void)
{
  result_job_state *job=&session.job; char description[512]; FILE *f; int id;
  if(!session.active || !job->open) return -1;
  ground_description(description,sizeof(description));
  if(job->last_ground_id>=0 && !strcmp(description,job->last_ground)) return job->last_ground_id;
  id=(int)job->next_ground++; f=spool_row(&job->grounds); if(!f) return -1;
  fprintf(f,"{\"ground_state_id\":%d,\"ground\":%s}",id,description);
  snprintf(job->last_ground,sizeof(job->last_ground),"%s",description);job->last_ground_id=id;
  return id;
}

void result_json_invalidate_ground(void)
{
  if(session.active) session.job.last_ground_id=-1;
}

int result_json_request(const char *kind, const char *parameters_json)
{
  result_job_state *job=&session.job;
  if(!session.active || !job->open) return -1;
  job->current_request=(int)job->next_request++;
  snprintf(job->current_trigger,sizeof(job->current_trigger),"%s",kind?kind:"UNKNOWN");
  snprintf(job->current_parameters,sizeof(job->current_parameters),"%s",parameters_json?parameters_json:"{}");
  return job->current_request;
}

int result_json_request_xq(int mode)
{
  char parameters[64];snprintf(parameters,sizeof(parameters),"{\"mode\":%d}",mode);
  return result_json_request("XQ",parameters);
}

int result_json_request_rp(int mode,int nt,int np,double ts,double ps,
    double dt,double dp,double range)
{
  char a[7][32],parameters[512];
  number_text(a[0],sizeof(a[0]),ts);number_text(a[1],sizeof(a[1]),ps);
  number_text(a[2],sizeof(a[2]),dt);number_text(a[3],sizeof(a[3]),dp);
  number_text(a[4],sizeof(a[4]),range);
  snprintf(parameters,sizeof(parameters),"{\"mode\":%d,\"theta_count\":%d,\"phi_count\":%d,\"theta_start_deg\":%s,\"phi_start_deg\":%s,\"theta_step_deg\":%s,\"phi_step_deg\":%s,\"range_m\":%s}",mode,nt,np,a[0],a[1],a[2],a[3],a[4]);
  return result_json_request("RP",parameters);
}

int result_json_request_near(const char *kind,int mode,int nx,int ny,int nz,
    const double start[3],const double step[3])
{
  char n[6][32],parameters[512];int i;
  for(i=0;i<3;++i){number_text(n[i],sizeof(n[i]),start[i]);number_text(n[i+3],sizeof(n[i+3]),step[i]);}
  snprintf(parameters,sizeof(parameters),"{\"coordinate_mode\":%d,\"counts\":[%d,%d,%d],\"start\":[%s,%s,%s],\"step\":[%s,%s,%s]}",mode,nx,ny,nz,n[0],n[1],n[2],n[3],n[4],n[5]);
  return result_json_request(kind,parameters);
}

int result_json_current_request(void)
{ return session.active ? session.job.current_request : -1; }

const char *result_json_current_excitation_kind(void)
{
  return session.active && session.job.execution.open ?
      session.job.execution.excitation_kind : "unknown";
}

static void write_excitation(FILE *f, result_execution_state *e)
{
  fputs("\"excitation\":{\"kind\":",f);json_string(f,e->excitation_kind);
  fputs(",\"phasor_convention\":\"exp(+i omega t)\"",f);
  if(e->has_wave)
  {
    double d[3];complex double v[3];
    dielectric_plane_wave_vectors(&e->wave,d,v);
    fputs(",\"propagation_direction\":[",f);json_number(f,d[0]);fputc(',',f);json_number(f,d[1]);fputc(',',f);json_number(f,d[2]);
    fputs("],\"electric_field_at_origin_v_per_m\":[",f);json_complex(f,v[0]);fputc(',',f);json_complex(f,v[1]);fputc(',',f);json_complex(f,v[2]);
    fprintf(f,"],\"ex_type\":%d,\"theta_deg\":",e->wave.type);json_number(f,e->wave.theta_degrees);
    fputs(",\"phi_deg\":",f);json_number(f,e->wave.phi_degrees);fputs(",\"eta_deg\":",f);json_number(f,e->wave.eta_degrees);
    fputs(",\"axial_ratio\":",f);json_number(f,e->wave.axial_ratio);
  }
  fputc('}',f);
}

static void write_nullable(FILE *f,int defined,double value)
{ if(defined && finite_number(value)) json_number(f,value); else fputs("null",f); }

static void write_solver(FILE *f,const result_execution_state *e)
{
  const result_solver *s=&e->solver;
  if(!e->solver_set) { fputs("null",f);return; }
  fputs("{\"requested\":",f);json_string(f,s->requested);fputs(",\"used\":",f);json_string(f,s->used);
  fputs(",\"relative_tolerance\":",f);write_nullable(f,s->iterative_fields_defined,s->relative_tolerance);
  fputs(",\"maximum_iterations\":",f);if(s->iterative_fields_defined)fprintf(f,"%lu",(unsigned long)s->maximum_iterations);else fputs("null",f);
  fputs(",\"iterations\":",f);if(s->iterative_fields_defined)fprintf(f,"%lu",(unsigned long)s->iterations);else fputs("null",f);
  fputs(",\"reported_residual\":",f);write_nullable(f,s->iterative_fields_defined,s->reported_residual);
  fputs(",\"true_residual\":",f);write_nullable(f,s->iterative_fields_defined,s->true_residual);
  fputs(",\"converged\":",f);if(s->iterative_fields_defined)fputs(s->converged?"true":"false",f);else fputs("null",f);
  fprintf(f,",\"accounted_storage_bytes\":%lu,\"elapsed_seconds\":",(unsigned long)s->accounted_storage_bytes);
  write_nullable(f,finite_number(s->elapsed_seconds),s->elapsed_seconds);fputc('}',f);
}

static void write_power(FILE *f,const result_execution_state *e)
{
  const result_power *p=&e->power;
  if(!e->power_set){fputs("null",f);return;}
  fputs("{\"basis\":",f);json_string(f,p->basis);
  fputs(",\"accepted_w\":",f);write_nullable(f,p->accepted_defined,p->accepted_w);
  fputs(",\"air_radiated_w\":",f);write_nullable(f,p->air_radiated_defined,p->air_radiated_w);
  fputs(",\"lower_halfspace_flux_w\":",f);write_nullable(f,p->lower_flux_defined,p->lower_flux_w);
  fputs(",\"dielectric_absorbed_w\":",f);write_nullable(f,p->dielectric_absorbed_defined,p->dielectric_absorbed_w);
  fputs(",\"structure_loss_w\":",f);write_nullable(f,p->structure_loss_defined,p->structure_loss_w);
  fputs(",\"network_loss_w\":",f);write_nullable(f,p->network_loss_defined,p->network_loss_w);
  fputs(",\"balance_defect_w\":",f);write_nullable(f,p->defect_defined,p->defect_w);
  fputs(",\"balance_defect_relative\":",f);write_nullable(f,p->defect_defined,p->defect_relative);
  fprintf(f,",\"angular_quadrature\":{\"theta_count\":%lu,\"phi_count\":%lu}}",
    (unsigned long)p->theta_count,(unsigned long)p->phi_count);
}

void result_json_execution_begin(double frequency_mhz,const char *kind,
    const dielectric_plane_wave *wave)
{
  result_job_state *job=&session.job;result_execution_state *e;
  if(!session.active || !job->open) return;
  if(job->execution.open) result_json_execution_end(1);
  e=&job->execution;memset(e,0,sizeof(*e));e->open=1;e->status_ok=1;
  e->index=job->next_execution++;e->frequency_mhz=frequency_mhz;
  e->mesh_state_id=job->current_mesh;e->ground_state_id=result_json_ground_state();
  e->request_index=job->current_request;snprintf(e->trigger_kind,sizeof(e->trigger_kind),"%s",job->current_trigger);
  snprintf(e->trigger_parameters,sizeof(e->trigger_parameters),"%s",job->current_parameters);
  snprintf(e->excitation_kind,sizeof(e->excitation_kind),"%s",kind?kind:"unknown");
  if(wave){e->wave=*wave;e->has_wave=1;}
  spool_open(&e->feeds);spool_open(&e->wires);spool_open(&e->patches);spool_open(&e->volumes);
  spool_open(&e->materials);spool_open(&e->far_fields);spool_open(&e->near_fields);
}

int result_json_execution_end(int success)
{
  result_job_state *job=&session.job;result_execution_state *e=&job->execution;FILE *f;
  if(!session.active || !e->open)return success!=0;
  if(e->far_open){fputs("]}",e->far_fields.file);e->far_open=0;}
  if(e->near_open){fputs("]}",e->near_fields.file);e->near_open=0;}
  if(session.failed && success)
  {
    result_json_diagnostic("error","NONFINITE_RESULT",
        "a non-finite numerical result was rejected","{}");
    success=0;
  }
  f=spool_row(&job->executions);if(!f){session.failed=1;execution_clear(e);return 0;}
  fprintf(f,"{\"execution_index\":%lu,\"status\":\"%s\",\"frequency_mhz\":",
    (unsigned long)e->index,success?"ok":"failed");json_number(f,e->frequency_mhz);
  fputs(",\"mesh_state_id\":",f);if(e->mesh_state_id>=0)fprintf(f,"%d",e->mesh_state_id);else fputs("null",f);
  fprintf(f,",\"ground_state_id\":%d,",e->ground_state_id);write_excitation(f,e);
  fprintf(f,",\"trigger\":{\"kind\":");json_string(f,e->trigger_kind);
  fprintf(f,",\"request_index\":%d,\"parameters\":%s},\"solver\":",e->request_index,e->trigger_parameters);
  write_solver(f,e);
  if(success)
  {
    fputs(",\"material_evaluations\":[",f);spool_copy(f,&e->materials);
    fputs("],\"feeds\":[",f);spool_copy(f,&e->feeds);fputs("],\"wire_currents\":[",f);spool_copy(f,&e->wires);
    fputs("],\"patch_currents\":{\"encoding\":\"rows-v1\",\"units\":\"A/m\",\"columns\":[\"patch_id\",\"center_x_m\",\"center_y_m\",\"center_z_m\",\"area_m2\",\"tangent1_x\",\"tangent1_y\",\"tangent1_z\",\"tangent2_x\",\"tangent2_y\",\"tangent2_z\",\"jx_re\",\"jx_im\",\"jy_re\",\"jy_im\",\"jz_re\",\"jz_im\"],\"rows\":[",f);spool_copy(f,&e->patches);fputs("]}",f);
    if(e->volumes.count)
    {
      fputs(",\"volume_currents\":{\"encoding\":\"rows-v1\",\"columns\":[\"cell_id\",\"box_id\",\"region_id\",\"material_id\",\"center_x_m\",\"center_y_m\",\"center_z_m\",\"size_x_m\",\"size_y_m\",\"size_z_m\",\"jx_re_a_per_m2\",\"jx_im_a_per_m2\",\"jy_re_a_per_m2\",\"jy_im_a_per_m2\",\"jz_re_a_per_m2\",\"jz_im_a_per_m2\"],\"rows\":[",f);spool_copy(f,&e->volumes);fputs("]}",f);
    } else fputs(",\"volume_currents\":null",f);
    fputs(",\"power\":",f);write_power(f,e);
    fputs(",\"far_fields\":[",f);spool_copy(f,&e->far_fields);fputs("],\"near_fields\":[",f);spool_copy(f,&e->near_fields);fputc(']',f);
  }
  fputc('}',f);execution_clear(e);return success!=0;
}

void result_json_set_solver(const result_solver *solver)
{
  if(session.active&&session.job.execution.open&&solver)
  {
    session.job.execution.solver=*solver;session.job.execution.solver_set=1;
    if(!finite_number(solver->elapsed_seconds) ||
       (solver->iterative_fields_defined &&
        (!finite_number(solver->relative_tolerance) ||
         !finite_number(solver->reported_residual) ||
         !finite_number(solver->true_residual)))) session.failed=1;
  }
}

void result_json_add_feed(size_t source_index,int tag,size_t segment,complex double voltage,
    complex double current,complex double impedance,int accepted_defined,double accepted)
{
  FILE *f;if(!session.active||!session.job.execution.open)return;f=spool_row(&session.job.execution.feeds);if(!f)return;
  fprintf(f,"{\"source_index\":%lu,\"tag\":%d,\"segment\":%lu,\"voltage_v\":",
    (unsigned long)source_index,tag,(unsigned long)segment);json_complex(f,voltage);
  fputs(",\"current_a\":",f);json_complex(f,current);fputs(",\"impedance_ohm\":",f);
  if(finite_number(creal(impedance))&&finite_number(cimag(impedance)))json_complex(f,impedance);else fputs("null",f);
  fputs(",\"accepted_power_w\":",f);write_nullable(f,accepted_defined,accepted);fputc('}',f);
}

void result_json_add_wire(int tag,size_t segment,const double c[3],double length,complex double current)
{
  FILE *f;if(!session.active||!session.job.execution.open)return;f=spool_row(&session.job.execution.wires);if(!f)return;
  fprintf(f,"{\"tag\":%d,\"segment\":%lu,\"center_m\":[",tag,(unsigned long)segment);
  json_number(f,c[0]);fputc(',',f);json_number(f,c[1]);fputc(',',f);json_number(f,c[2]);
  fputs("],\"length_m\":",f);json_number(f,length);fputs(",\"current_a\":",f);json_complex(f,current);fputc('}',f);
}

void result_json_add_patch(size_t id,const double c[3],double area,const double t1[3],const double t2[3],const complex double j[3])
{
  FILE *f;int p;if(!session.active||!session.job.execution.open)return;f=spool_row(&session.job.execution.patches);if(!f)return;
  fprintf(f,"[%lu",(unsigned long)id);for(p=0;p<3;++p){fputc(',',f);json_number(f,c[p]);}fputc(',',f);json_number(f,area);
  for(p=0;p<3;++p){fputc(',',f);json_number(f,t1[p]);}for(p=0;p<3;++p){fputc(',',f);json_number(f,t2[p]);}
  for(p=0;p<3;++p){fputc(',',f);json_number(f,creal(j[p]));fputc(',',f);json_number(f,cimag(j[p]));}fputc(']',f);
}

void result_json_add_volume(size_t cell,size_t box,int region,int material,const double c[3],const double s[3],const complex double j[3])
{
  FILE *f;int p;if(!session.active||!session.job.execution.open)return;f=spool_row(&session.job.execution.volumes);if(!f)return;
  fprintf(f,"[%lu,%lu,%d,%d",(unsigned long)cell,(unsigned long)box,region,material);
  for(p=0;p<3;++p){fputc(',',f);json_number(f,c[p]);}for(p=0;p<3;++p){fputc(',',f);json_number(f,s[p]);}
  for(p=0;p<3;++p){fputc(',',f);json_number(f,creal(j[p]));fputc(',',f);json_number(f,cimag(j[p]));}fputc(']',f);
}

void result_json_add_material_evaluations(const dielectric_model *model,double frequency)
{
  size_t m;if(!session.active||!session.job.execution.open||!model)return;
  for(m=0;m<model->material_count;++m)if(model->materials[m].tabulated)
  {
    const dielectric_material *mat=&model->materials[m];size_t hi=1;double t,k0=2*3.14159265358979323846*frequency*1e6/299792458.0;complex double er;FILE *f;
    while(hi+1<mat->sample_count&&frequency>mat->samples[hi].frequency_mhz)++hi;
    t=(frequency-mat->samples[hi-1].frequency_mhz)/(mat->samples[hi].frequency_mhz-mat->samples[hi-1].frequency_mhz);er=dielectric_material_epsilon(mat,k0);
    f=spool_row(&session.job.execution.materials);if(!f)return;
    fprintf(f,"{\"material_id\":%d,\"frequency_mhz\":",mat->id);json_number(f,frequency);
    fprintf(f,",\"lower_sample_index\":%lu,\"upper_sample_index\":%lu,\"interpolation_fraction\":",(unsigned long)(hi-1),(unsigned long)hi);json_number(f,t);
    fputs(",\"er\":",f);json_number(f,(1-t)*mat->samples[hi-1].epsilon_r+t*mat->samples[hi].epsilon_r);
    fputs(",\"loss_kind\":",f);json_string(f,mat->samples[hi-1].conductivity_kind?"sigma":"epsilon-loss");
    fputs(",\"loss_value\":",f);json_number(f,(1-t)*mat->samples[hi-1].loss+t*mat->samples[hi].loss);
    fputs(",\"complex_relative_permittivity\":",f);json_complex(f,er);fputc('}',f);
  }
}

void result_json_set_power(const result_power *power)
{
  if(!session.active||!session.job.execution.open||!power)return;
  session.job.execution.power=*power;session.job.execution.power_set=1;
  if((power->accepted_defined&&!finite_number(power->accepted_w)) ||
     (power->air_radiated_defined&&!finite_number(power->air_radiated_w)) ||
     (power->lower_flux_defined&&!finite_number(power->lower_flux_w)) ||
     (power->dielectric_absorbed_defined&&!finite_number(power->dielectric_absorbed_w)) ||
     (power->structure_loss_defined&&!finite_number(power->structure_loss_w)) ||
     (power->network_loss_defined&&!finite_number(power->network_loss_w)) ||
     (power->defect_defined&&(!finite_number(power->defect_w)||!finite_number(power->defect_relative))))
    session.failed=1;
  if(power->defect_defined&&finite_number(power->defect_relative)&&power->defect_relative>0.003)
  {
    char context[192], defect[96];
    number_text(defect,sizeof(defect),power->defect_relative);
    snprintf(context,sizeof(context),"{\"job_index\":%lu,\"execution_index\":%lu,\"actual_relative_defect\":%s,\"warning_threshold\":0.003}",(unsigned long)session.job.index,(unsigned long)session.job.execution.index,defect);
    result_json_diagnostic("warning","POWER_BALANCE_EXCEEDED","relative power-balance defect exceeds 0.3 percent",context);
  }
}

void result_json_far_begin(size_t request,const char *kind)
{
  result_execution_state *e=&session.job.execution;FILE *f;if(!session.active||!e->open)return;
  if(e->far_open){fputs("]}",e->far_fields.file);e->far_open=0;}
  f=spool_row(&e->far_fields);if(!f)return;
  fprintf(f,"{\"request_index\":%lu,\"kind\":",(unsigned long)request);json_string(f,kind);
  fputs(",\"normalization\":\"range-times-electric-field\",\"units\":\"V\",\"samples\":[",f);e->far_open=1;e->far_samples=0;
}

void result_json_far_sample(double theta,double phi,complex double et,complex double ep)
{
  result_execution_state *e=&session.job.execution;FILE *f;if(!session.active||!e->open||!e->far_open)return;f=e->far_fields.file;
  if(e->far_samples++)fputc(',',f);
  fputs("{\"theta_deg\":",f);json_number(f,theta);fputs(",\"phi_deg\":",f);json_number(f,phi);
  fputs(",\"r_e_theta\":",f);json_complex(f,et);fputs(",\"r_e_phi\":",f);json_complex(f,ep);fputc('}',f);
}

void result_json_near_begin(size_t request,const char *field,const char *kind,const char *units)
{
  result_execution_state *e=&session.job.execution;FILE *f;if(!session.active||!e->open)return;
  if(e->near_open){fputs("]}",e->near_fields.file);e->near_open=0;}
  f=spool_row(&e->near_fields);if(!f)return;
  fprintf(f,"{\"request_index\":%lu,\"field\":",(unsigned long)request);json_string(f,field);fputs(",\"kind\":",f);json_string(f,kind);fputs(",\"units\":",f);json_string(f,units);fputs(",\"samples\":[",f);e->near_open=1;e->near_samples=0;
}

void result_json_near_sample(const double p[3],const complex double v[3],int valid)
{
  result_execution_state *e=&session.job.execution;FILE *f;int i;if(!session.active||!e->open||!e->near_open)return;f=e->near_fields.file;
  if(e->near_samples++)fputc(',',f);
  fputs("{\"position_m\":[",f);for(i=0;i<3;++i){if(i)fputc(',',f);json_number(f,p[i]);}
  fputs("],\"value\":{\"x\":",f);json_complex(f,v[0]);fputs(",\"y\":",f);json_complex(f,v[1]);fputs(",\"z\":",f);json_complex(f,v[2]);fputs("},\"valid\":",f);fputs(valid?"true":"false",f);fputc('}',f);
}

void result_json_omit(int request,const char *kind,const char *code,const char *message)
{
  FILE *f;if(!session.active||!session.job.open)return;f=spool_row(&session.job.omissions);if(!f)return;
  fprintf(f,"{\"request_index\":%d,\"kind\":",request);json_string(f,kind);fputs(",\"code\":",f);json_string(f,code);fputs(",\"message\":",f);json_string(f,message);fputc('}',f);
}

void result_json_diagnostic(const char *severity,const char *code,const char *message,const char *context)
{
  FILE *f;if(!session.active)return;f=spool_row(&session.diagnostics);if(!f)return;
  fputs("{\"severity\":",f);json_string(f,severity);fputs(",\"code\":",f);json_string(f,code);fputs(",\"message\":",f);json_string(f,message);fputs(",\"context\":",f);fputs(context?context:"{}",f);fputc('}',f);
}

static const char *diagnostic_code(const dielectric_diagnostic *d)
{
  if(!strcmp(d->record,"RESULT"))return "NONFINITE_RESULT";
  if(strstr(d->message,"unsupported format version"))return "DIELECTRIC_RESOURCE_VERSION";
  if(strstr(d->message,"interiors overlap"))return "DIELECTRIC_BOX_OVERLAP";
  if(strstr(d->message,"conductor clearance")||strstr(d->message,"dielectric-to-conductor"))return "DIELECTRIC_CONDUCTOR_CLEARANCE";
  if(strstr(d->message,"ground")||!strcmp(d->record,"GN"))return "DIELECTRIC_GROUND_CLEARANCE";
  if(strstr(d->message,"near-field")||strstr(d->message,"observation"))return "DIELECTRIC_OBSERVATION_CLEARANCE";
  if(strstr(d->message,"material table")||!strcmp(d->record,"TABLE"))return "DIELECTRIC_MATERIAL_RANGE";
  if(strstr(d->message,"MAX_CELLS")||strstr(d->message,"cell limit"))return "DIELECTRIC_CELL_LIMIT";
  if(!strcmp(d->record,"FFT"))return "DIELECTRIC_FFT_LATTICE";
  if(strstr(d->message,"iteration limit")||strstr(d->value,"GMRES"))return "SOLVER_NOT_CONVERGED";
  if(!strcmp(d->record,"SOLVER")||!strcmp(d->record,"COUPLED_SOLVER"))return "DIELECTRIC_SOLVER_INCOMPATIBLE";
  return "DIELECTRIC_RESOURCE_PARSE";
}

void result_json_dielectric_diagnostic(const dielectric_diagnostic *d,const char *severity)
{
  char context[192];if(!session.active||!d)return;
  snprintf(context,sizeof(context),"{\"job_index\":%lu,\"line\":%lu}",(unsigned long)session.job.index,(unsigned long)d->line);
  result_json_diagnostic(severity?severity:"error",diagnostic_code(d),d->message,context);
}

int result_json_finish(const char *status)
{
  int ok, semantic_ok;
  const char *effective=status?status:"internal-error";
  if(!session.active)return 1;
  if(session.job.open)result_json_job_end(NULL);
  fputs("],\"diagnostics\":[",session.document);spool_copy(session.document,&session.diagnostics);
  if(session.failed && !strcmp(effective,"ok")) effective="solver-error";
  semantic_ok=!session.failed;
  fputs("],\"completion\":{\"status\":",session.document);json_string(session.document,effective);fputs("}}\n",session.document);
  ok=!ferror(session.document)&&fflush(session.document)==0&&fclose(session.document)==0;
  session.document=NULL;spool_close(&session.diagnostics);
  if(ok)ok=rename(session.temporary,session.destination)==0;
  if(!ok)remove(session.temporary);
  session.active=0;return ok&&semantic_ok;
}

void result_json_abandon(void)
{
  if(session.document)fclose(session.document);
  if(session.temporary[0])remove(session.temporary);
  spool_close(&session.diagnostics);job_clear(&session.job);memset(&session,0,sizeof(session));
}
