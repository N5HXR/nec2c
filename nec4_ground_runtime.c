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

void nec4_ground_runtime_init(nec4_ground_runtime *runtime)
{ if(runtime) memset(runtime,0,sizeof(*runtime)); }

void nec4_ground_runtime_destroy(nec4_ground_runtime *runtime)
{
  if(!runtime) return;
  for(int index=0;index<3;++index) nec4_ground_table_destroy(&runtime->table[index]);
  nec4_ground_table_destroy(&runtime->model_grid);
  memset(runtime,0,sizeof(*runtime));
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

int nec4_ground_runtime_prepare(nec4_ground_runtime *runtime,
    double frequency,double complex epsilon,double complex lower_k,
    double complex upper_k,const nec4_ground_runtime_options *options)
{
  nec4_ground_runtime_options defaults={0};
  const nec4_ground_table_direct_options direct_defaults={3,48,80.};
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
  runtime->direct_options=selected->direct_options?*selected->direct_options:
                          direct_defaults;
  runtime->direct_model=selected->direct_model!=0;
  if(cabs(epsilon-1.)<=64.*DBL_EPSILON&&
     cabs(lower_k-upper_k)<=64.*DBL_EPSILON*fmax(1.,cabs(upper_k)))
  {
    runtime->identical_media=1;
    runtime->matrix_context.frequency_hz=frequency;
    runtime->matrix_context.lower_epsilon_relative=epsilon;
    runtime->matrix_context.lower_wave_number_per_m=lower_k;
    runtime->matrix_context.upper_wave_number_per_m=upper_k;
    runtime->matrix_context.cross=runtime_cross;
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
  if(!runtime->direct_model) {
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
  runtime->matrix_tables.lower_wave_number_per_m=lower_k;
  runtime->matrix_tables.upper_wave_number_per_m=upper_k;
  runtime->matrix_context.frequency_hz=frequency;
  runtime->matrix_context.lower_epsilon_relative=epsilon;
  runtime->matrix_context.lower_wave_number_per_m=lower_k;
  runtime->matrix_context.upper_wave_number_per_m=upper_k;
  runtime->matrix_context.cross=runtime_cross;
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
    int success=nec4_ground_matrix_table_cross((void *)&runtime->matrix_tables,
                                                rho,z,depth,field);
    if(!success) fprintf(stderr,"nec2c: NEC4 linear table failed at "
        "rho=%.17g z=%.17g depth=%.17g lg=%.17g\n",
        rho,z,depth,runtime->lg_m);
    return success;
  }
  if(route.method==NEC4_GROUND_ASYMPTOTIC) {
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
    if(!runtime->direct_model&&z>=.16*runtime->wavelength0_m&&
       nec4_ground_model_grid_evaluate(&runtime->model_grid,
         runtime->matrix_context.lower_wave_number_per_m,
         runtime->matrix_context.upper_wave_number_per_m,rho,z,depth,field))
      return 1;
    route.used_direct_fallback=1;if(route_out) *route_out=route;
    {
      nec4_ground_runtime *mutable_runtime=(nec4_ground_runtime *)runtime;
#ifdef _OPENMP
#pragma omp atomic write
#endif
      mutable_runtime->direct_strip_seen=1;
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
    else fprintf(stderr,"nec2c: NEC4 direct contour failed at "
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
