#include "nec4_ground_table.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct { double start,step;size_t count; } table_axis;

typedef struct {
  unsigned int kind;
  double frequency_hz;
  double complex lower_k,upper_k;
  nec4_ground_table_direct_options options;
} direct_sample_context;

static int multiply_size(size_t left,size_t right,size_t *result)
{
  if(left&&right>SIZE_MAX/left) return 0;
  *result=left*right;
  return 1;
}

static int make_axis(double start,double stop,double maximum_step,
    table_axis *axis)
{
  size_t intervals;
  if(!axis||!isfinite(start)||!isfinite(stop)||!isfinite(maximum_step)||
     !(stop>start)||!(maximum_step>0.)) return 0;
  if(ceil((stop-start)/maximum_step)>(double)(SIZE_MAX-1)) return 0;
  intervals=(size_t)ceil((stop-start)/maximum_step);
  if(intervals<1) intervals=1;
  axis->start=start;axis->step=(stop-start)/intervals;
  axis->count=intervals+1;
  return axis->step>0.&&isfinite(axis->step);
}

int nec4_ground_table_init(nec4_ground_table *table,unsigned int kind,
    double frequency_hz,double complex epsilon_relative,double wavelength0_m,
    double lg_m,double lmax_m)
{
  static const double rho_steps[3]={.02,.03,.03};
  static const double observer_steps[3]={.02,.03,.03};
  static const double depth_steps[3]={.02,.015,.03};
  static const double depth_starts[3]={0.,.4,1.};
  table_axis rho,observer,depth;
  double depth_stop;
  size_t points,values;
  if(!table||kind<1||kind>3||!(frequency_hz>0.)||!(wavelength0_m>0.)||
     !(lg_m>0.)||!(lmax_m>=lg_m)||!isfinite(frequency_hz)||
     !isfinite(creal(epsilon_relative))||!isfinite(cimag(epsilon_relative))||
     !isfinite(wavelength0_m)||!isfinite(lg_m)||!isfinite(lmax_m)) return 0;
  memset(table,0,sizeof(*table));
  depth_stop=kind==1?.4:(kind==2?1.:lmax_m/lg_m);
  if(!make_axis(0.,lg_m,rho_steps[kind-1]*lg_m,&rho)||
     !make_axis(0.,lg_m,observer_steps[kind-1]*lg_m,&observer)||
     !make_axis(depth_starts[kind-1]*lg_m,depth_stop*lg_m,
                depth_steps[kind-1]*lg_m,&depth)||
     !multiply_size(rho.count,observer.count,&points)||
     !multiply_size(points,depth.count,&points)||
     !multiply_size(points,NEC4_GROUND_COMPONENT_COUNT,&values)||
     values>SIZE_MAX/sizeof(double complex)) return 0;
  table->values=calloc(values,sizeof(*table->values));
  if(!table->values) return 0;
  table->kind=kind;table->frequency_hz=frequency_hz;
  table->epsilon_relative=epsilon_relative;
  table->wavelength0_m=wavelength0_m;table->lg_m=lg_m;
  table->lmax_m=lmax_m;
  table->rho_start_m=rho.start;table->rho_step_m=rho.step;
  table->observer_start_m=observer.start;table->observer_step_m=observer.step;
  table->source_depth_start_m=depth.start;
  table->source_depth_step_m=depth.step;
  table->rho_count=rho.count;table->observer_count=observer.count;
  table->source_depth_count=depth.count;table->value_count=values;
  return 1;
}

void nec4_ground_table_destroy(nec4_ground_table *table)
{
  if(!table) return;
  free(table->values);memset(table,0,sizeof(*table));
}

int nec4_ground_table_fill(nec4_ground_table *table,
    nec4_ground_table_sampler sampler,void *context,int workers)
{
  size_t points,index;
  int failed=0,effective_workers=workers>0?workers:1;
  if(!table||!table->values||!sampler||
     !multiply_size(table->rho_count,table->observer_count,&points)||
     !multiply_size(points,table->source_depth_count,&points)) return 0;
#ifdef _OPENMP
  if(workers<=0) effective_workers=omp_get_max_threads();
#pragma omp parallel for schedule(static) num_threads(effective_workers) reduction(|:failed)
#endif
  for(index=0;index<points;++index)
  {
    size_t remaining=index;
    size_t depth_index=remaining%table->source_depth_count;
    size_t observer_index;
    size_t rho_index;
    double complex sample[NEC4_GROUND_COMPONENT_COUNT];
    remaining/=table->source_depth_count;
    observer_index=remaining%table->observer_count;
    rho_index=remaining/table->observer_count;
    if(table->kind==1&&rho_index==0&&observer_index==0&&depth_index==0)
      memset(sample,0,sizeof(sample));
    else if(!sampler(context,
        table->rho_start_m+rho_index*table->rho_step_m,
        table->observer_start_m+observer_index*table->observer_step_m,
        table->source_depth_start_m+depth_index*table->source_depth_step_m,
        sample)) failed=1;
    else memcpy(table->values+index*NEC4_GROUND_COMPONENT_COUNT,sample,
                sizeof(sample));
  }
  (void)effective_workers;
  return !failed;
}

static int storage_factor(unsigned int kind,double complex lower_k,
    double complex upper_k,double rho,double observer_z,double source_depth,
    double complex *factor)
{
  double radius;
  double complex phase;
  if(!factor||kind<1||kind>3||rho<0.||observer_z<0.||source_depth<0.)
    return 0;
  radius=hypot(rho,observer_z+source_depth);
  if(radius==0.)
  {
    if(kind!=1) return 0;
    *factor=0.;return 1;
  }
  if(kind==1) *factor=radius*cexp(I*upper_k*radius);
  else if(kind==2)
  {
    phase=upper_k*hypot(rho,observer_z)+lower_k*source_depth;
    *factor=radius*radius*cexp(I*phase);
  }
  else
  {
    phase=upper_k*observer_z+lower_k*hypot(rho,source_depth);
    *factor=radius*radius*cexp(I*phase);
  }
  return isfinite(creal(*factor))&&isfinite(cimag(*factor));
}

static int direct_sample(void *opaque,double rho,double observer_z,
    double source_depth,double complex values[NEC4_GROUND_COMPONENT_COUNT])
{
  direct_sample_context *context=opaque;
  nec4_ground_direct_result ordinary;
  nec4_ground_table1_result first;
  nec4_ground_interface_result interface;
  double complex factor;
  const double complex *physical;
  unsigned int component;
  int use_hankel;
  if(!context||!values||!storage_factor(context->kind,context->lower_k,
      context->upper_k,rho,observer_z,source_depth,&factor)) return 0;
  use_hankel=rho>0.&&observer_z+source_depth<=.5*rho;
  if(context->kind==1)
  {
    if(observer_z+source_depth==0.)
    {
      if(!nec4_ground_direct_cross_interface(context->frequency_hz,
          context->lower_k,context->upper_k,rho,
          context->options.panel_subdivisions,
          context->options.quadrature_order,context->options.cutoff_decay,
          &interface)) return 0;
      physical=interface.field.finite_remainder;
    }
    else if(use_hankel)
    {
      if(!nec4_ground_direct_cross_table1_hankel(context->frequency_hz,
          context->lower_k,context->upper_k,rho,observer_z,source_depth,
          context->options.panel_subdivisions,
          context->options.quadrature_order,context->options.cutoff_decay,
          &first)) return 0;
      physical=first.finite_remainder;
    }
    else
    {
      if(!nec4_ground_direct_cross_table1(context->frequency_hz,
          context->lower_k,context->upper_k,rho,observer_z,source_depth,
          context->options.quadrature_order,context->options.cutoff_decay,
          &first)) return 0;
      physical=first.finite_remainder;
    }
  }
  else
  {
    if(use_hankel)
    {
      if(!nec4_ground_direct_cross_hankel(context->frequency_hz,
          context->lower_k,context->upper_k,rho,observer_z,source_depth,
          context->options.panel_subdivisions,
          context->options.quadrature_order,context->options.cutoff_decay,
          &ordinary)) return 0;
    }
    else if(!nec4_ground_direct_cross_bessel(context->frequency_hz,
        context->lower_k,context->upper_k,rho,observer_z,source_depth,
        context->options.quadrature_order,context->options.cutoff_decay,
        &ordinary)) return 0;
    physical=ordinary.value;
  }
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    values[component]=physical[component]*factor;
  return 1;
}

int nec4_ground_table_fill_direct(nec4_ground_table *table,
    double complex lower_k,double complex upper_k,
    const nec4_ground_table_direct_options *options,int workers)
{
  const nec4_ground_table_direct_options defaults={3,48,80.};
  direct_sample_context context;
  double complex material_ratio;
  double material_scale;
  if(!table||!table->values||!(cabs(lower_k)>0.)||!(cabs(upper_k)>0.))
    return 0;
  material_ratio=(lower_k/upper_k)*(lower_k/upper_k);
  material_scale=fmax(1.,cabs(table->epsilon_relative));
  if(cabs(material_ratio-table->epsilon_relative)>
     1.e-10*material_scale||
     cabs(upper_k*table->wavelength0_m-
          6.28318530717958647692528676655900577)>1.e-10) return 0;
  context.kind=table->kind;context.frequency_hz=table->frequency_hz;
  context.lower_k=lower_k;context.upper_k=upper_k;
  context.options=options?*options:defaults;
  if(context.options.panel_subdivisions<1||
     context.options.quadrature_order<12||
     context.options.cutoff_decay<30.||
     !isfinite(context.options.cutoff_decay)) return 0;
  return nec4_ground_table_fill(table,direct_sample,&context,workers);
}

static int locate(double coordinate,double start,double step,size_t count,
    size_t *lower,double *fraction)
{
  double position,end,tolerance;
  if(!lower||!fraction||count<2||!(step>0.)||!isfinite(coordinate)) return 0;
  end=start+(count-1)*step;
  tolerance=64.*DBL_EPSILON*fmax(1.,fmax(fabs(start),fabs(end)));
  if(coordinate<start-tolerance||coordinate>end+tolerance) return 0;
  if(coordinate<=start) { *lower=0;*fraction=0.;return 1; }
  if(coordinate>=end) { *lower=count-2;*fraction=1.;return 1; }
  position=(coordinate-start)/step;
  *lower=(size_t)floor(position);
  if(*lower>=count-1) *lower=count-2;
  *fraction=position-*lower;
  return 1;
}

int nec4_ground_table_interpolate(const nec4_ground_table *table,double rho_m,
    double observer_z_m,double source_depth_m,
    double complex values[NEC4_GROUND_COMPONENT_COUNT])
{
  size_t lower[3],corner_rho,corner_observer,corner_depth,component;
  double fraction[3];
  if(!table||!table->values||!values||
     !locate(rho_m,table->rho_start_m,table->rho_step_m,table->rho_count,
             &lower[0],&fraction[0])||
     !locate(observer_z_m,table->observer_start_m,table->observer_step_m,
             table->observer_count,&lower[1],&fraction[1])||
     !locate(source_depth_m,table->source_depth_start_m,
             table->source_depth_step_m,table->source_depth_count,
             &lower[2],&fraction[2])) return 0;
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    values[component]=0.;
  for(corner_rho=0;corner_rho<2;++corner_rho)
    for(corner_observer=0;corner_observer<2;++corner_observer)
      for(corner_depth=0;corner_depth<2;++corner_depth)
      {
        size_t point=((lower[0]+corner_rho)*table->observer_count+
            lower[1]+corner_observer)*table->source_depth_count+
            lower[2]+corner_depth;
        double weight=(corner_rho?fraction[0]:1.-fraction[0])*
          (corner_observer?fraction[1]:1.-fraction[1])*
          (corner_depth?fraction[2]:1.-fraction[2]);
        for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
          values[component]+=weight*table->values[
              point*NEC4_GROUND_COMPONENT_COUNT+component];
      }
  return 1;
}

int nec4_ground_table_evaluate(const nec4_ground_table *table,
    double complex lower_k,double complex upper_k,double rho_m,
    double observer_z_m,double source_depth_m,
    double complex values[NEC4_GROUND_COMPONENT_COUNT])
{
  double complex stored[NEC4_GROUND_COMPONENT_COUNT],factor;
  unsigned int component;
  if(!values||!nec4_ground_table_interpolate(table,rho_m,observer_z_m,
      source_depth_m,stored)||!storage_factor(table->kind,lower_k,upper_k,
      rho_m,observer_z_m,source_depth_m,&factor)||cabs(factor)==0.) return 0;
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    stored[component]/=factor;
  if(table->kind==1)
    return nec4_ground_direct_table1_restore(table->frequency_hz,lower_k,
        upper_k,rho_m,observer_z_m,source_depth_m,stored,values);
  memcpy(values,stored,sizeof(stored));
  return 1;
}

void nec4_ground_table_cache_key(const nec4_ground_table *table,
    nec4_ground_cache_key *key)
{
  if(!key) return;
  memset(key,0,sizeof(*key));
  if(!table) return;
  key->equation_version=NEC4_GROUND_EQUATION_VERSION;
  key->precision_bits=53;key->table_kind=table->kind;
  key->component_count=NEC4_GROUND_COMPONENT_COUNT;
  key->frequency_hz=table->frequency_hz;
  key->epsilon_real=creal(table->epsilon_relative);
  key->epsilon_imag=cimag(table->epsilon_relative);
  key->wavelength0_m=table->wavelength0_m;
  key->lg_m=table->lg_m;key->lmax_m=table->lmax_m;
  key->rho_start_m=table->rho_start_m;key->rho_step_m=table->rho_step_m;
  key->observer_start_m=table->observer_start_m;
  key->observer_step_m=table->observer_step_m;
  key->source_depth_start_m=table->source_depth_start_m;
  key->source_depth_step_m=table->source_depth_step_m;
  key->rho_count=table->rho_count;key->observer_count=table->observer_count;
  key->source_depth_count=table->source_depth_count;
  key->value_count=table->value_count;
}

nec4_ground_cache_status nec4_ground_table_store(
    const nec4_ground_table *table,const char *path)
{
  nec4_ground_cache_key key;
  if(!table||!table->values) return NEC4_CACHE_IO_ERROR;
  nec4_ground_table_cache_key(table,&key);
  return nec4_ground_cache_store(path,&key,table->values);
}

nec4_ground_cache_status nec4_ground_table_load(
    nec4_ground_table *table,const char *path)
{
  nec4_ground_cache_key key;
  if(!table||!table->values) return NEC4_CACHE_IO_ERROR;
  nec4_ground_table_cache_key(table,&key);
  return nec4_ground_cache_load(path,&key,table->values);
}
