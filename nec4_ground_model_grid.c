#include "nec4_ground_model_grid.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MODEL_GRID_KIND 4u

static int product(size_t a,size_t b,size_t *result)
{ if(!result||(a!=0&&b>SIZE_MAX/a))return 0;*result=a*b;return 1; }

int nec4_ground_model_grid_init(nec4_ground_table *grid,double frequency,
    double complex epsilon,double wavelength,double lg,double lmax)
{
  size_t points,values;
  if(!grid||!(frequency>0.)||!(wavelength>0.)||!(lg>0.)||!(lmax>=lg))return 0;
  memset(grid,0,sizeof(*grid));
  grid->kind=MODEL_GRID_KIND;grid->frequency_hz=frequency;
  grid->epsilon_relative=epsilon;grid->wavelength0_m=wavelength;
  grid->lg_m=lg;grid->lmax_m=lmax;
  grid->rho_start_m=0.;grid->rho_step_m=wavelength/32.;grid->rho_count=97;
  grid->observer_start_m=5.*wavelength/32.;
  grid->observer_step_m=wavelength/32.;grid->observer_count=60;
  grid->source_depth_start_m=0.;grid->source_depth_step_m=.025*lg;
  grid->source_depth_count=33;
  if(!product(grid->rho_count,grid->observer_count,&points)||
     !product(points,grid->source_depth_count,&points)||
     !product(points,NEC4_GROUND_COMPONENT_COUNT,&values))return 0;
  grid->value_count=values;grid->values=calloc(values,sizeof(*grid->values));
  if(!grid->values){memset(grid,0,sizeof(*grid));return 0;}return 1;
}

static int f2(double complex k1,double complex k2,double rho,double z,
    double depth,double complex *factor)
{
  double radius=hypot(rho,z+depth);
  double complex phase=k2*hypot(rho,z)+k1*depth;
  if(!factor||radius==0.)return 0;
  *factor=radius*radius*cexp(I*phase);
  return isfinite(creal(*factor))&&isfinite(cimag(*factor));
}

typedef struct {
  double frequency;
  double complex k1,k2;
  nec4_ground_table_direct_options options;
} fill_context;

static int direct_sample(void *opaque,double rho,double z,double depth,
    double complex values[NEC4_GROUND_COMPONENT_COUNT])
{
  fill_context *context=opaque;nec4_ground_direct_result direct;
  double complex factor;int success;
  if(!context||!values||!f2(context->k1,context->k2,rho,z,depth,&factor))return 0;
  if(rho>0.&&z+depth<=.5*rho)
    success=nec4_ground_direct_cross_hankel(context->frequency,context->k1,
      context->k2,rho,z,depth,context->options.panel_subdivisions,
      context->options.quadrature_order,context->options.cutoff_decay,&direct);
  else success=nec4_ground_direct_cross_bessel(context->frequency,context->k1,
      context->k2,rho,z,depth,context->options.quadrature_order,
      context->options.cutoff_decay,&direct);
  if(!success)return 0;
  for(unsigned int component=0;component<NEC4_GROUND_COMPONENT_COUNT;component++)
    values[component]=direct.value[component]*factor;
  return 1;
}

int nec4_ground_model_grid_fill_direct(nec4_ground_table *grid,
    double complex k1,double complex k2,
    const nec4_ground_table_direct_options *options,int workers)
{
  const nec4_ground_table_direct_options defaults={3,48,80.};
  fill_context context;
  if(!grid||grid->kind!=MODEL_GRID_KIND||!grid->values)return 0;
  context.frequency=grid->frequency_hz;context.k1=k1;context.k2=k2;
  context.options=options?*options:defaults;
  return nec4_ground_table_fill(grid,direct_sample,&context,workers);
}

int nec4_ground_model_grid_evaluate(const nec4_ground_table *grid,
    double complex k1,double complex k2,double rho,double z,double depth,
    double complex values[NEC4_GROUND_COMPONENT_COUNT])
{
  double complex stored[NEC4_GROUND_COMPONENT_COUNT],factor;
  if(!grid||grid->kind!=MODEL_GRID_KIND||!values||
     !nec4_ground_table_interpolate(grid,rho,z,depth,stored)||
     !f2(k1,k2,rho,z,depth,&factor)||cabs(factor)==0.)return 0;
  for(unsigned int component=0;component<NEC4_GROUND_COMPONENT_COUNT;component++)
    values[component]=stored[component]/factor;
  return 1;
}
