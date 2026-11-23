#include "nec4_ground_adaptive_e.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#define ADAPTIVE_E_ALGORITHM_VERSION 0x000b0001u
#define ADAPTIVE_E_REPRESENTATION_F2 1u

typedef struct {
  uint32_t algorithm_version,representation;
  uint32_t panel_subdivisions,quadrature_order;
  double cutoff_decay,frequency_hz;
  double lower_real,lower_imag,upper_real,upper_imag;
  double bounds_lower[3],bounds_upper[3],tolerance;
  uint64_t maximum_cells;
  uint32_t maximum_depth;
} adaptive_e_identity;

static adaptive_e_identity identity(double frequency,double complex lower_k,
    double complex upper_k,const nec4_ground_table_direct_options *options,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells)
{
  adaptive_e_identity result;
  memset(&result,0,sizeof(result));
  result.algorithm_version=ADAPTIVE_E_ALGORITHM_VERSION;
  result.representation=ADAPTIVE_E_REPRESENTATION_F2;
  result.panel_subdivisions=options->panel_subdivisions;
  result.quadrature_order=options->quadrature_order;
  result.cutoff_decay=options->cutoff_decay;result.frequency_hz=frequency;
  result.lower_real=creal(lower_k);result.lower_imag=cimag(lower_k);
  result.upper_real=creal(upper_k);result.upper_imag=cimag(upper_k);
  memcpy(result.bounds_lower,bounds->lower,sizeof(result.bounds_lower));
  memcpy(result.bounds_upper,bounds->upper,sizeof(result.bounds_upper));
  result.tolerance=tolerance;result.maximum_depth=maximum_depth;
  result.maximum_cells=maximum_cells;
  return result;
}

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
    double complex values[5])
{
  nec4_ground_adaptive_e *table=opaque;
  nec4_ground_direct_result direct;
  nec4_ground_interface_result interface;
  const double complex *physical;
  double complex factor;
  if(!table||!values||!storage_factor(table->lower_wave_number_per_m,
      table->upper_wave_number_per_m,rho,z,depth,&factor)) return 0;
  if(z+depth==0.) {
    if(!nec4_ground_direct_cross_interface(table->frequency_hz,
        table->lower_wave_number_per_m,table->upper_wave_number_per_m,rho,
        table->direct_options.panel_subdivisions,
        table->direct_options.quadrature_order,
        table->direct_options.cutoff_decay,&interface)) return 0;
    physical=interface.field.total;
  } else if(rho>0.&&z+depth<=.5*rho) {
    if(!nec4_ground_direct_cross_hankel(table->frequency_hz,
        table->lower_wave_number_per_m,table->upper_wave_number_per_m,rho,z,
        depth,table->direct_options.panel_subdivisions,
        table->direct_options.quadrature_order,
        table->direct_options.cutoff_decay,&direct)) return 0;
    physical=direct.value;
  } else {
    if(!nec4_ground_direct_cross_bessel(table->frequency_hz,
        table->lower_wave_number_per_m,table->upper_wave_number_per_m,rho,z,
        depth,table->direct_options.quadrature_order,
        table->direct_options.cutoff_decay,&direct)) return 0;
    physical=direct.value;
  }
  for(unsigned int component=0;component<5;++component)
    values[component]=physical[component]*factor;
  return 1;
}

void nec4_ground_adaptive_e_init(nec4_ground_adaptive_e *table)
{
  if(!table) return;
  memset(table,0,sizeof(*table));nec4_ground_adaptive_init(&table->table);
}

void nec4_ground_adaptive_e_destroy(nec4_ground_adaptive_e *table)
{
  if(!table) return;
  nec4_ground_adaptive_destroy(&table->table);memset(table,0,sizeof(*table));
}

int nec4_ground_adaptive_e_build(nec4_ground_adaptive_e *table,
    double frequency,double complex lower_k,double complex upper_k,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,int workers,
    const nec4_ground_table_direct_options *options)
{
  const nec4_ground_table_direct_options defaults={3,48,80.};
  if(!table||!(frequency>0.)||cabs(lower_k)==0.||cabs(upper_k)==0.) return 0;
  nec4_ground_adaptive_e_destroy(table);
  table->frequency_hz=frequency;table->lower_wave_number_per_m=lower_k;
  table->upper_wave_number_per_m=upper_k;
  table->direct_options=options?*options:defaults;
  if(table->direct_options.panel_subdivisions<1||
     table->direct_options.quadrature_order<12||
     table->direct_options.cutoff_decay<30.) return 0;
  return nec4_ground_adaptive_build(&table->table,bounds,tolerance,
      maximum_depth,maximum_cells,workers,direct_sample,table);
}

int nec4_ground_adaptive_e_build_anisotropic(nec4_ground_adaptive_e *table,
    double frequency,double complex lower_k,double complex upper_k,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,int workers,
    const nec4_ground_table_direct_options *options)
{
  const nec4_ground_table_direct_options defaults={3,48,80.};
  if(!table||!(frequency>0.)||cabs(lower_k)==0.||cabs(upper_k)==0.) return 0;
  nec4_ground_adaptive_e_destroy(table);
  table->frequency_hz=frequency;table->lower_wave_number_per_m=lower_k;
  table->upper_wave_number_per_m=upper_k;
  table->direct_options=options?*options:defaults;
  if(table->direct_options.panel_subdivisions<1||
     table->direct_options.quadrature_order<12||
     table->direct_options.cutoff_decay<30.) return 0;
  return nec4_ground_adaptive_build_anisotropic(&table->table,bounds,tolerance,
      maximum_depth,maximum_cells,workers,direct_sample,table);
}

int nec4_ground_adaptive_e_evaluate(const nec4_ground_adaptive_e *table,
    double rho,double z,double depth,double complex values[5])
{
  double complex stored[5],factor;
  if(!table||!values||
     !nec4_ground_adaptive_evaluate(&table->table,rho,z,depth,stored)||
     !storage_factor(table->lower_wave_number_per_m,
       table->upper_wave_number_per_m,rho,z,depth,&factor)||cabs(factor)==0.)
    return 0;
  for(unsigned int component=0;component<5;++component)
    values[component]=stored[component]/factor;
  return 1;
}

int nec4_ground_adaptive_e_store(const nec4_ground_adaptive_e *table,
    const char *path)
{
  adaptive_e_identity key;
  if(!table) return 0;
  key=identity(table->frequency_hz,table->lower_wave_number_per_m,
               table->upper_wave_number_per_m,&table->direct_options,
               &table->table.bounds,table->table.tolerance,
               table->table.maximum_depth,table->table.maximum_cells);
  return nec4_ground_adaptive_store(&table->table,path,&key,sizeof(key));
}

int nec4_ground_adaptive_e_load(nec4_ground_adaptive_e *table,const char *path,
    double frequency,double complex lower_k,double complex upper_k,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,
    const nec4_ground_table_direct_options *options)
{
  const nec4_ground_table_direct_options defaults={3,48,80.};
  const nec4_ground_table_direct_options *selected=options?options:&defaults;
  adaptive_e_identity key;
  nec4_ground_adaptive loaded;
  if(!table||!bounds||!(frequency>0.)||!(tolerance>0.)||
     cabs(lower_k)==0.||cabs(upper_k)==0.)return 0;
  key=identity(frequency,lower_k,upper_k,selected,bounds,tolerance,
               maximum_depth,maximum_cells);
  nec4_ground_adaptive_init(&loaded);
  if(!nec4_ground_adaptive_load(&loaded,path,&key,sizeof(key)))return 0;
  nec4_ground_adaptive_e_destroy(table);table->table=loaded;
  table->frequency_hz=frequency;table->lower_wave_number_per_m=lower_k;
  table->upper_wave_number_per_m=upper_k;table->direct_options=*selected;
  return 1;
}
