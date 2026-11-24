#define _POSIX_C_SOURCE 200809L

#include "nec4_ground_runtime.h"
#include "nec4_ground_asymptotic.h"

#include <errno.h>
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#define NEC4_C0 299792458.0

static int runtime_cross(void *opaque,double rho,double z,double depth,
    double complex field[NEC4_GROUND_COMPONENT_COUNT])
{
  nec4_ground_runtime *runtime=opaque;
  if(runtime&&runtime->identical_media)
    return nec4_ground_matrix_homogeneous_components(
      runtime->matrix_context.frequency_hz,
      runtime->matrix_context.upper_wave_number_per_m,rho,-depth,z,field);
  return nec4_ground_runtime_evaluate(opaque,rho,z,depth,field,NULL);
}

static int runtime_cross_h(void *opaque,double rho,double z,double depth,
    double complex field[NEC4_GROUND_H_COMPONENT_COUNT])
{ return nec4_ground_runtime_evaluate_h(opaque,rho,z,depth,field,NULL); }

void nec4_ground_runtime_init(nec4_ground_runtime *runtime)
{ if(runtime) memset(runtime,0,sizeof(*runtime)); }

void nec4_ground_runtime_destroy(nec4_ground_runtime *runtime)
{
  if(!runtime) return;
  for(int index=0;index<3;++index) nec4_ground_table_destroy(&runtime->table[index]);
  nec4_ground_table_destroy(&runtime->model_grid);
  for(int index=0;index<NEC4_GROUND_ADAPTIVE_REGION_COUNT;++index)
    nec4_ground_adaptive_e_destroy(&runtime->adaptive_e[index]);
  for(int index=0;index<NEC4_GROUND_ADAPTIVE_REGION_COUNT;++index)
    nec4_ground_adaptive_h_destroy(&runtime->adaptive_h[index]);
  memset(runtime,0,sizeof(*runtime));
}

const char *nec4_ground_model_name(nec4_ground_model model)
{
  switch(model) {
    case NEC4_GROUND_MODEL_GRID:return "grid";
    case NEC4_GROUND_MODEL_DIRECT:return "direct";
    case NEC4_GROUND_MODEL_ADAPTIVE:return "adaptive";
    default:return "invalid";
  }
}

static int cache_directory_ready(const char *path)
{
  struct stat state;char copy[4096];size_t length;
  if(!path||!*path) return 1;
  if(!stat(path,&state)) return S_ISDIR(state.st_mode);
  if(errno!=ENOENT) return 0;
  length=strlen(path);
  if(length>=sizeof(copy)) return 0;
  memcpy(copy,path,length+1);
  for(char *cursor=copy+1;*cursor;++cursor) if(*cursor=='/')
  {
    *cursor='\0';
    if(stat(copy,&state))
    {
      if(errno!=ENOENT||mkdir(copy,0700)) return 0;
    }
    else if(!S_ISDIR(state.st_mode)) return 0;
    *cursor='/';
  }
  return !mkdir(copy,0700)||errno==EEXIST;
}

static int cache_path(const nec4_ground_table *table,const char *directory,
    char *path,size_t capacity)
{
  nec4_ground_cache_key key;
  uint64_t hash;
  int length;
  if(!directory||!*directory||!path||capacity==0) return 0;
  nec4_ground_table_cache_key(table,&key);
  hash=nec4_ground_cache_key_hash(&key);
  if(!hash) return 0;
  length=snprintf(path,capacity,"%s/nec4-ground-%016" PRIx64 ".bin",
      directory,hash);
  return length>0&&(size_t)length<capacity;
}

static int adaptive_cache_path(uint64_t hash,const char *directory,char *path,
    size_t capacity)
{
  int length;
  if(!hash||!directory||!*directory||!path||capacity==0) return 0;
  length=snprintf(path,capacity,"%s/nec4-adaptive-e-%016" PRIx64 ".bin",
    directory,hash);
  return length>0&&(size_t)length<capacity;
}

static int adaptive_h_cache_path(uint64_t hash,const char *directory,char *path,
    size_t capacity)
{
  int length;
  if(!hash||!directory||!*directory||!path||capacity==0) return 0;
  length=snprintf(path,capacity,"%s/nec4-adaptive-h-%016" PRIx64 ".bin",
    directory,hash);
  return length>0&&(size_t)length<capacity;
}

enum {
  ADAPTIVE_SWITCH_VOLUME=0,
  ADAPTIVE_INTERFACE_LINE=1,
  ADAPTIVE_GRAZING_VOLUME=2
};

typedef struct {
  nec4_ground_adaptive_bounds bounds;
  double tolerance;
  unsigned int maximum_depth;
  size_t maximum_cells;
  int anisotropic;
} adaptive_specification;

static adaptive_specification adaptive_spec(int region,double wavelength,
    double lg)
{
  adaptive_specification result;
  memset(&result,0,sizeof(result));result.tolerance=5.e-5;
  result.maximum_depth=30;result.anisotropic=1;
  switch(region) {
    case 0: /* route-switch volume */
      result.bounds.lower[0]=2.8*wavelength;
      result.bounds.upper[0]=3.2*wavelength;
      result.bounds.lower[1]=.015*wavelength;
      result.bounds.upper[1]=.035*wavelength;
      result.bounds.lower[2]=.002*lg;result.bounds.upper[2]=.05*lg;
      result.maximum_cells=4096;break;
    case 1: /* exact-interface line */
      result.bounds.lower[0]=.5*wavelength;
      result.bounds.upper[0]=4.*wavelength;
      result.maximum_depth=10;result.maximum_cells=2048;
      result.anisotropic=0;break;
    default: /* grazing sheet */
      result.bounds.lower[0]=3.*wavelength;
      result.bounds.upper[0]=5.*wavelength;
      result.bounds.lower[1]=.01*wavelength;
      result.bounds.upper[1]=.4*wavelength;
      result.bounds.lower[2]=result.bounds.upper[2]=.02*lg;
      result.maximum_cells=4096;break;
  }
  return result;
}

static int prepare_adaptive_e(nec4_ground_runtime *runtime,int region,
    double frequency,double complex lower_k,double complex upper_k,
    const char *cache_directory,int workers)
{
  adaptive_specification specification=adaptive_spec(region,
    runtime->wavelength0_m,runtime->lg_m);
  nec4_ground_adaptive_e *table=&runtime->adaptive_e[region];
  uint64_t hash=nec4_ground_adaptive_e_cache_key(frequency,lower_k,upper_k,
    &specification.bounds,specification.tolerance,specification.maximum_depth,
    specification.maximum_cells,specification.anisotropic,
    &runtime->direct_options);
  char path[4096];
  int have_path=adaptive_cache_path(hash,cache_directory,path,sizeof(path));
  nec4_ground_adaptive_e_init(table);
  if(have_path) {
    int loaded=specification.anisotropic?
      nec4_ground_adaptive_e_load_anisotropic(table,path,frequency,lower_k,
        upper_k,&specification.bounds,specification.tolerance,
        specification.maximum_depth,specification.maximum_cells,
        &runtime->direct_options):
      nec4_ground_adaptive_e_load(table,path,frequency,lower_k,upper_k,
        &specification.bounds,specification.tolerance,
        specification.maximum_depth,specification.maximum_cells,
        &runtime->direct_options);
    if(loaded) { ++runtime->cache_hits;runtime->adaptive_e_ready[region]=1;
      return 1; }
  }
  ++runtime->cache_misses;
  if(specification.anisotropic) {
    if(!nec4_ground_adaptive_e_build_anisotropic(table,frequency,lower_k,
        upper_k,&specification.bounds,specification.tolerance,
        specification.maximum_depth,specification.maximum_cells,workers,
        &runtime->direct_options)) return 0;
  } else if(!nec4_ground_adaptive_e_build(table,frequency,lower_k,upper_k,
      &specification.bounds,specification.tolerance,
      specification.maximum_depth,specification.maximum_cells,workers,
      &runtime->direct_options)) return 0;
  ++runtime->generated_tables;runtime->adaptive_e_ready[region]=1;
  if(have_path&&!nec4_ground_adaptive_e_store(table,path)) return 0;
  return 1;
}

static int prepare_adaptive_h(nec4_ground_runtime *runtime,int region,
    double frequency,double complex lower_k,double complex upper_k,
    const char *cache_directory,int workers)
{
  adaptive_specification specification=adaptive_spec(region,
    runtime->wavelength0_m,runtime->lg_m);
  nec4_ground_adaptive_h *table=&runtime->adaptive_h[region];
  uint64_t hash;char path[4096];int have_path;
  /* M11D qualified a two-dimensional grazing sheet at this normalized source
     depth.  Do not silently broaden it into an unqualified H volume. */
  if(region==ADAPTIVE_GRAZING_VOLUME)
    specification.bounds.lower[2]=specification.bounds.upper[2]=.02*runtime->lg_m;
  if(region==ADAPTIVE_INTERFACE_LINE) return 1;
  specification.tolerance=3.e-4;specification.maximum_cells=4096;
  specification.maximum_depth=30;specification.anisotropic=1;
  hash=nec4_ground_adaptive_h_cache_key(frequency,lower_k,upper_k,
    &specification.bounds,specification.tolerance,specification.maximum_depth,
    specification.maximum_cells,&runtime->direct_options);
  have_path=adaptive_h_cache_path(hash,cache_directory,path,sizeof(path));
  nec4_ground_adaptive_h_init(table);
  if(have_path&&nec4_ground_adaptive_h_load(table,path,frequency,lower_k,
      upper_k,&specification.bounds,specification.tolerance,
      specification.maximum_depth,specification.maximum_cells,
      &runtime->direct_options)) {
    ++runtime->cache_hits;runtime->adaptive_h_ready[region]=1;return 1;
  }
  ++runtime->cache_misses;
  if(!nec4_ground_adaptive_h_build_anisotropic(table,frequency,lower_k,upper_k,
      &specification.bounds,specification.tolerance,
      specification.maximum_depth,specification.maximum_cells,workers,
      &runtime->direct_options)) return 0;
  ++runtime->generated_tables;runtime->adaptive_h_ready[region]=1;
  if(have_path&&!nec4_ground_adaptive_h_store(table,path)) return 0;
  return 1;
}

int nec4_ground_runtime_prepare(nec4_ground_runtime *runtime,
    double frequency,double complex epsilon,double complex lower_k,
    double complex upper_k,const nec4_ground_runtime_options *options)
{
  nec4_ground_runtime_options defaults={0};
  const nec4_ground_table_direct_options direct_defaults={3,48,80.};
  const nec4_ground_table_direct_options adaptive_defaults={5,80,120.};
  double ratio,wavelength,lg,lmax;
  const nec4_ground_runtime_options *selected=options?options:&defaults;
  if(!runtime||!(frequency>0.)||cabs(lower_k)==0.||cabs(upper_k)==0.||
     !isfinite(frequency)||!isfinite(creal(epsilon))||!isfinite(cimag(epsilon)))
    return 0;
  ratio=cabs(upper_k/lower_k);wavelength=NEC4_C0/frequency;
  lg=fmin(ratio,.3)*wavelength;lmax=fmin(15.*ratio,1.)*wavelength;
  if(!(lg>0.)||!(lmax>=lg)||
     !cache_directory_ready(selected->cache_directory)) return 0;
  nec4_ground_runtime_destroy(runtime);
  runtime->wavelength0_m=wavelength;runtime->lg_m=lg;runtime->lmax_m=lmax;
  if(selected->model<NEC4_GROUND_MODEL_GRID||
     selected->model>NEC4_GROUND_MODEL_ADAPTIVE) goto failed;
  runtime->model=selected->model;runtime->workers=selected->workers;
  runtime->direct_options=selected->direct_options?*selected->direct_options:
    (runtime->model==NEC4_GROUND_MODEL_ADAPTIVE?adaptive_defaults:direct_defaults);
  if(cabs(epsilon-1.)<=64.*DBL_EPSILON&&
     cabs(lower_k-upper_k)<=64.*DBL_EPSILON*fmax(1.,cabs(upper_k)))
  {
    runtime->identical_media=1;
    runtime->matrix_context.frequency_hz=frequency;
    runtime->matrix_context.lower_epsilon_relative=epsilon;
    runtime->matrix_context.lower_wave_number_per_m=lower_k;
    runtime->matrix_context.upper_wave_number_per_m=upper_k;
    runtime->matrix_context.cross=runtime_cross;
    runtime->matrix_context.cross_h=runtime->model==NEC4_GROUND_MODEL_GRID?
      NULL:runtime_cross_h;
    runtime->matrix_context.cross_context=runtime;
    runtime->ready=1;
    return 1;
  }
  for(int index=0;index<3;++index) {
    char path[4096];
    nec4_ground_cache_status status=NEC4_CACHE_NOT_FOUND;
    int have_path;
    if(!nec4_ground_table_init(&runtime->table[index],(unsigned int)index+1,
        frequency,epsilon,wavelength,lg,lmax)) goto failed;
    runtime->payload_bytes+=runtime->table[index].value_count*
      sizeof(*runtime->table[index].values);
    have_path=cache_path(&runtime->table[index],selected->cache_directory,
      path,sizeof(path));
    if(have_path) status=nec4_ground_table_load(&runtime->table[index],path);
    if(status==NEC4_CACHE_OK) ++runtime->cache_hits;
    else
    {
      ++runtime->cache_misses;
      if(selected->sampler)
      {
        if(!nec4_ground_table_fill(&runtime->table[index],selected->sampler,
            selected->sampler_context,selected->workers)) goto failed;
      }
      else if(!nec4_ground_table_fill_direct(&runtime->table[index],lower_k,
          upper_k,selected->direct_options,selected->workers)) goto failed;
      ++runtime->generated_tables;
      if(have_path&&nec4_ground_table_store(&runtime->table[index],path)!=
          NEC4_CACHE_OK) goto failed;
    }
    runtime->matrix_tables.table[index]=&runtime->table[index];
  }
  if(runtime->model==NEC4_GROUND_MODEL_GRID||
     runtime->model==NEC4_GROUND_MODEL_ADAPTIVE) {
    char path[4096];nec4_ground_cache_status status=NEC4_CACHE_NOT_FOUND;
    int have_path;
    if(!nec4_ground_model_grid_init(&runtime->model_grid,frequency,epsilon,
                                    wavelength,lg,lmax))goto failed;
    runtime->payload_bytes+=runtime->model_grid.value_count*
                            sizeof(*runtime->model_grid.values);
    have_path=cache_path(&runtime->model_grid,selected->cache_directory,path,
                         sizeof(path));
    if(have_path)status=nec4_ground_table_load(&runtime->model_grid,path);
    if(status==NEC4_CACHE_OK)++runtime->cache_hits;
    else {
      ++runtime->cache_misses;
      if(selected->sampler) {
        if(!nec4_ground_table_fill(&runtime->model_grid,selected->sampler,
                                   selected->sampler_context,
                                   selected->workers))goto failed;
      } else if(!nec4_ground_model_grid_fill_direct(&runtime->model_grid,
          lower_k,upper_k,selected->direct_options,selected->workers))goto failed;
      ++runtime->generated_tables;
      if(have_path&&nec4_ground_table_store(&runtime->model_grid,path)!=
          NEC4_CACHE_OK)goto failed;
    }
  }
  if(runtime->model==NEC4_GROUND_MODEL_ADAPTIVE)
    for(int region=0;region<NEC4_GROUND_ADAPTIVE_REGION_COUNT;++region)
      if(!prepare_adaptive_e(runtime,region,frequency,lower_k,upper_k,
          selected->cache_directory,selected->workers)) goto failed;
  if(runtime->model==NEC4_GROUND_MODEL_ADAPTIVE)
    for(int region=0;region<NEC4_GROUND_ADAPTIVE_REGION_COUNT;++region)
      if(!prepare_adaptive_h(runtime,region,frequency,lower_k,upper_k,
          selected->cache_directory,selected->workers)) goto failed;
  runtime->matrix_tables.lower_wave_number_per_m=lower_k;
  runtime->matrix_tables.upper_wave_number_per_m=upper_k;
  runtime->matrix_context.frequency_hz=frequency;
  runtime->matrix_context.lower_epsilon_relative=epsilon;
  runtime->matrix_context.lower_wave_number_per_m=lower_k;
  runtime->matrix_context.upper_wave_number_per_m=upper_k;
  runtime->matrix_context.cross=runtime_cross;
  runtime->matrix_context.cross_h=runtime->model==NEC4_GROUND_MODEL_GRID?
    NULL:runtime_cross_h;
  runtime->matrix_context.cross_context=runtime;
  runtime->ready=1;
  return 1;
failed:
  nec4_ground_runtime_destroy(runtime);
  return 0;
}

const nec4_ground_matrix_context *nec4_ground_runtime_matrix_context(
    const nec4_ground_runtime *runtime)
{ return runtime&&runtime->ready?&runtime->matrix_context:NULL; }

int nec4_ground_runtime_evaluate(const nec4_ground_runtime *runtime,
    double rho,double z,double depth,double complex field[5],
    nec4_ground_dispatch *route_out)
{
  nec4_ground_dispatch route;
  if(!runtime||!runtime->ready||!field||
     !nec4_ground_route(rho,z,-depth,
       runtime->matrix_context.lower_wave_number_per_m,
       runtime->matrix_context.upper_wave_number_per_m,runtime->wavelength0_m,
       &route)) return 0;
  if(route.grazing_risk) {
    nec4_ground_runtime *mutable_runtime=(nec4_ground_runtime *)runtime;
#ifdef _OPENMP
#pragma omp atomic write
#endif
    mutable_runtime->grazing_seen=1;
  }
  if(route_out) *route_out=route;
  if(route.method==NEC4_GROUND_LINEAR_TABLE) {
    nec4_ground_runtime *mutable_runtime=(nec4_ground_runtime *)runtime;
#ifdef _OPENMP
#pragma omp atomic update
#endif
    ++mutable_runtime->linear_hits;
    int success=nec4_ground_matrix_table_cross((void *)&runtime->matrix_tables,
                                                rho,z,depth,field);
    if(!success) fprintf(stderr,"nec2c: NEC4ish linear table failed at "
        "rho=%.17g z=%.17g depth=%.17g lg=%.17g\n",
        rho,z,depth,runtime->lg_m);
    return success;
  }
  if(route.method==NEC4_GROUND_ASYMPTOTIC) {
    if(runtime->model==NEC4_GROUND_MODEL_ADAPTIVE&&route.grazing_risk) {
      int order[]={ADAPTIVE_GRAZING_VOLUME,ADAPTIVE_SWITCH_VOLUME};
      for(size_t index=0;index<sizeof(order)/sizeof(order[0]);++index)
        if(runtime->adaptive_e_ready[order[index]]&&
           nec4_ground_adaptive_e_evaluate(&runtime->adaptive_e[order[index]],
             rho,z,depth,field)) {
          nec4_ground_runtime *mutable_runtime=(nec4_ground_runtime *)runtime;
#ifdef _OPENMP
#pragma omp atomic update
#endif
          ++mutable_runtime->adaptive_hits;
          return 1;
        }
      /* A documented grazing-risk point is never sent to the original-NEC4
         asymptotic approximation in adaptive mode. */
      goto direct_contour;
    }
    {
      nec4_ground_runtime *mutable_runtime=(nec4_ground_runtime *)runtime;
#ifdef _OPENMP
#pragma omp atomic update
#endif
      ++mutable_runtime->asymptotic_hits;
    }
    nec4_ground_asymptotic_result result;
    if(!nec4_ground_asymptotic_second_uniform(
        runtime->matrix_context.frequency_hz,
        runtime->matrix_context.lower_wave_number_per_m,
        runtime->matrix_context.upper_wave_number_per_m,rho,z,depth,&result))
      return 0;
    memcpy(field,result.value,sizeof(result.value));
    return 1;
  }
  if(route.method==NEC4_GROUND_MODEL_FIT) {
    nec4_ground_direct_result direct;
    nec4_ground_interface_result interface;
    int success;
    if(runtime->model==NEC4_GROUND_MODEL_ADAPTIVE) {
      int order[]={ADAPTIVE_SWITCH_VOLUME,ADAPTIVE_INTERFACE_LINE,
                   ADAPTIVE_GRAZING_VOLUME};
      for(size_t index=0;index<sizeof(order)/sizeof(order[0]);++index)
        if(runtime->adaptive_e_ready[order[index]]&&
           nec4_ground_adaptive_e_evaluate(&runtime->adaptive_e[order[index]],
             rho,z,depth,field)) {
          nec4_ground_runtime *mutable_runtime=(nec4_ground_runtime *)runtime;
#ifdef _OPENMP
#pragma omp atomic update
#endif
          ++mutable_runtime->adaptive_hits;
          return 1;
        }
    }
    if((runtime->model==NEC4_GROUND_MODEL_GRID||
        runtime->model==NEC4_GROUND_MODEL_ADAPTIVE)&&
       z>=.16*runtime->wavelength0_m&&
       nec4_ground_model_grid_evaluate(&runtime->model_grid,
         runtime->matrix_context.lower_wave_number_per_m,
         runtime->matrix_context.upper_wave_number_per_m,rho,z,depth,field))
      return 1;
direct_contour:
    route.used_direct_fallback=1;if(route_out) *route_out=route;
    {
      nec4_ground_runtime *mutable_runtime=(nec4_ground_runtime *)runtime;
#ifdef _OPENMP
#pragma omp atomic write
#endif
      mutable_runtime->direct_strip_seen=1;
#ifdef _OPENMP
#pragma omp atomic update
#endif
      ++mutable_runtime->direct_fallbacks;
    }
    if(z+depth==0.) {
      success=nec4_ground_direct_cross_interface(
        runtime->matrix_context.frequency_hz,
        runtime->matrix_context.lower_wave_number_per_m,
        runtime->matrix_context.upper_wave_number_per_m,rho,
        runtime->direct_options.panel_subdivisions,
        runtime->direct_options.quadrature_order,
        runtime->direct_options.cutoff_decay,&interface);
      if(success) memcpy(field,interface.field.total,sizeof(interface.field.total));
      return success;
    }
    if(rho>0.&&z+depth<=.5*rho)
      success=nec4_ground_direct_cross_hankel(
        runtime->matrix_context.frequency_hz,
        runtime->matrix_context.lower_wave_number_per_m,
        runtime->matrix_context.upper_wave_number_per_m,rho,z,depth,
        runtime->direct_options.panel_subdivisions,
        runtime->direct_options.quadrature_order,
        runtime->direct_options.cutoff_decay,&direct);
    else {
      success=nec4_ground_direct_cross_bessel(
        runtime->matrix_context.frequency_hz,
        runtime->matrix_context.lower_wave_number_per_m,
        runtime->matrix_context.upper_wave_number_per_m,rho,z,depth,
        runtime->direct_options.quadrature_order,
        runtime->direct_options.cutoff_decay,&direct);
      /* The contour scale is derived from z+depth.  At isolated floating-
         point coincidences a quadrature node can land exactly on a spectral
         denominator zero even though the physical integral is continuous.
         Moving the observer height by one representable number selects the
         same limiting field while avoiding that artificial node collision. */
      if(!success&&z>0.)
        success=nec4_ground_direct_cross_bessel(
          runtime->matrix_context.frequency_hz,
          runtime->matrix_context.lower_wave_number_per_m,
          runtime->matrix_context.upper_wave_number_per_m,rho,
          nextafter(z,INFINITY),depth,
          runtime->direct_options.quadrature_order,
          runtime->direct_options.cutoff_decay,&direct);
      /* Near the cylindrical axis the exact radial limit is substantially
         better conditioned than evaluating Bessel derivatives at a tiny
         complex argument.  The omitted term is quadratic in rho/(z+d). */
      if(!success&&rho<.01*(z+depth))
        success=nec4_ground_direct_cross_bessel(
          runtime->matrix_context.frequency_hz,
          runtime->matrix_context.lower_wave_number_per_m,
          runtime->matrix_context.upper_wave_number_per_m,0.,z,depth,
          runtime->direct_options.quadrature_order,
          runtime->direct_options.cutoff_decay,&direct);
    }
    if(success) memcpy(field,direct.value,sizeof(direct.value));
    else fprintf(stderr,"nec2c: NEC4ish direct contour failed at "
        "rho=%.17g z=%.17g depth=%.17g\n",rho,z,depth);
    return success;
  }
  return 0;
}

int nec4_ground_runtime_grazing_seen(const nec4_ground_runtime *runtime)
{ return runtime&&runtime->grazing_seen; }

int nec4_ground_runtime_claim_grazing_warning(nec4_ground_runtime *runtime)
{
  if(!runtime||!runtime->grazing_seen||runtime->grazing_reported) return 0;
  runtime->grazing_reported=1;
  return 1;
}

int nec4_ground_runtime_direct_strip_seen(const nec4_ground_runtime *runtime)
{ return runtime&&runtime->direct_strip_seen; }

int nec4_ground_runtime_evaluate_h(const nec4_ground_runtime *runtime,
    double rho,double z,double depth,
    double complex field[NEC4_GROUND_H_COMPONENT_COUNT],
    nec4_ground_dispatch *route_out)
{
  nec4_ground_dispatch route;
  nec4_ground_direct_h_result direct;
  nec4_ground_h_interface_result interface;
  int success;
  if(!runtime||!runtime->ready||!field||
     !nec4_ground_route(rho,z,-depth,
       runtime->matrix_context.lower_wave_number_per_m,
       runtime->matrix_context.upper_wave_number_per_m,runtime->wavelength0_m,
       &route)) return 0;
  if(route_out) *route_out=route;
  if(runtime->model==NEC4_GROUND_MODEL_ADAPTIVE) {
    int order[]={ADAPTIVE_SWITCH_VOLUME,ADAPTIVE_GRAZING_VOLUME};
    for(size_t index=0;index<sizeof(order)/sizeof(order[0]);++index)
      if(runtime->adaptive_h_ready[order[index]]&&
         nec4_ground_adaptive_h_evaluate(&runtime->adaptive_h[order[index]],
           rho,z,depth,field)) {
        nec4_ground_runtime *mutable_runtime=(nec4_ground_runtime *)runtime;
#ifdef _OPENMP
#pragma omp atomic update
#endif
        ++mutable_runtime->h_adaptive_hits;
        return 1;
      }
  }
  route.used_direct_fallback=1;if(route_out) *route_out=route;
  {
    nec4_ground_runtime *mutable_runtime=(nec4_ground_runtime *)runtime;
#ifdef _OPENMP
#pragma omp atomic update
#endif
    ++mutable_runtime->h_direct_fallbacks;
  }
  if(z+depth==0.) {
    success=nec4_ground_direct_cross_h_interface(
      runtime->matrix_context.frequency_hz,
      runtime->matrix_context.lower_wave_number_per_m,
      runtime->matrix_context.upper_wave_number_per_m,rho,
      runtime->direct_options.panel_subdivisions,
      runtime->direct_options.quadrature_order,
      runtime->direct_options.cutoff_decay,&interface);
    if(success) memcpy(field,interface.field.value,sizeof(interface.field.value));
    return success;
  }
  if(rho>0.&&z+depth<=.5*rho)
    success=nec4_ground_direct_cross_h_hankel(
      runtime->matrix_context.frequency_hz,
      runtime->matrix_context.lower_wave_number_per_m,
      runtime->matrix_context.upper_wave_number_per_m,rho,z,depth,
      runtime->direct_options.panel_subdivisions,
      runtime->direct_options.quadrature_order,
      runtime->direct_options.cutoff_decay,&direct);
  else
    success=nec4_ground_direct_cross_h_bessel(
      runtime->matrix_context.frequency_hz,
      runtime->matrix_context.lower_wave_number_per_m,
      runtime->matrix_context.upper_wave_number_per_m,rho,z,depth,
      runtime->direct_options.quadrature_order,
      runtime->direct_options.cutoff_decay,&direct);
  if(success) memcpy(field,direct.value,sizeof(direct.value));
  return success;
}
