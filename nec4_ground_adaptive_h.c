#include "nec4_ground_adaptive_h.h"

#include <math.h>
#include <string.h>

static int storage_factor(double complex lower_k,double complex upper_k,
    double rho,double z,double depth,double complex *factor)
{
  double radius=hypot(rho,z+depth);
  double complex phase=upper_k*hypot(rho,z)+lower_k*depth;
  if(!factor||radius==0.) return 0;
  *factor=radius*radius*cexp(I*phase);
  return isfinite(creal(*factor))&&isfinite(cimag(*factor));
}

static int direct_sample(void *opaque,double rho,double z,double depth,
    double complex values[NEC4_GROUND_COMPONENT_COUNT])
{
  nec4_ground_adaptive_h *table=opaque;
  nec4_ground_direct_h_result direct;
  nec4_ground_h_interface_result interface;
  const double complex *physical;
  double complex factor;
  if(!table||!values||!storage_factor(table->lower_wave_number_per_m,
      table->upper_wave_number_per_m,rho,z,depth,&factor)) return 0;
  if(z+depth==0.) {
    if(!nec4_ground_direct_cross_h_interface(table->frequency_hz,
        table->lower_wave_number_per_m,table->upper_wave_number_per_m,rho,
        table->direct_options.panel_subdivisions,
        table->direct_options.quadrature_order,
        table->direct_options.cutoff_decay,&interface)) return 0;
    physical=interface.field.value;
  } else if(rho>0.&&z+depth<=.5*rho) {
    if(!nec4_ground_direct_cross_h_hankel(table->frequency_hz,
        table->lower_wave_number_per_m,table->upper_wave_number_per_m,rho,z,
        depth,table->direct_options.panel_subdivisions,
        table->direct_options.quadrature_order,
        table->direct_options.cutoff_decay,&direct)) return 0;
    physical=direct.value;
  } else {
    if(!nec4_ground_direct_cross_h_bessel(table->frequency_hz,
        table->lower_wave_number_per_m,table->upper_wave_number_per_m,rho,z,
        depth,table->direct_options.quadrature_order,
        table->direct_options.cutoff_decay,&direct)) return 0;
    physical=direct.value;
  }
  for(unsigned int component=0;
      component<NEC4_GROUND_H_COMPONENT_COUNT;++component)
    values[component]=physical[component]*factor;
  values[NEC4_GROUND_H_COMPONENT_COUNT]=0.;
  return 1;
}

void nec4_ground_adaptive_h_init(nec4_ground_adaptive_h *table)
{
  if(!table) return;
  memset(table,0,sizeof(*table));nec4_ground_adaptive_init(&table->table);
}

void nec4_ground_adaptive_h_destroy(nec4_ground_adaptive_h *table)
{
  if(!table) return;
  nec4_ground_adaptive_destroy(&table->table);memset(table,0,sizeof(*table));
}

int nec4_ground_adaptive_h_build_anisotropic(nec4_ground_adaptive_h *table,
    double frequency,double complex lower_k,double complex upper_k,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,int workers,
    const nec4_ground_table_direct_options *options)
{
  const nec4_ground_table_direct_options defaults={3,48,80.};
  if(!table||!(frequency>0.)||cabs(lower_k)==0.||cabs(upper_k)==0.) return 0;
  nec4_ground_adaptive_h_destroy(table);
  table->frequency_hz=frequency;table->lower_wave_number_per_m=lower_k;
  table->upper_wave_number_per_m=upper_k;
  table->direct_options=options?*options:defaults;
  if(table->direct_options.panel_subdivisions<1||
     table->direct_options.quadrature_order<12||
     table->direct_options.cutoff_decay<30.) return 0;
  return nec4_ground_adaptive_build_anisotropic(&table->table,bounds,tolerance,
      maximum_depth,maximum_cells,workers,direct_sample,table);
}

int nec4_ground_adaptive_h_evaluate(const nec4_ground_adaptive_h *table,
    double rho,double z,double depth,
    double complex values[NEC4_GROUND_H_COMPONENT_COUNT])
{
  double complex stored[NEC4_GROUND_COMPONENT_COUNT],factor;
  if(!table||!values||
     !nec4_ground_adaptive_evaluate(&table->table,rho,z,depth,stored)||
     !storage_factor(table->lower_wave_number_per_m,
       table->upper_wave_number_per_m,rho,z,depth,&factor)||cabs(factor)==0.)
    return 0;
  for(unsigned int component=0;
      component<NEC4_GROUND_H_COMPONENT_COUNT;++component)
    values[component]=stored[component]/factor;
  return 1;
}
