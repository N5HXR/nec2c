#define _POSIX_C_SOURCE 200809L

#include "nec4_ground_table.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct { double bias;int reject_origin; } sample_context;

static int affine_sample(void *context,double rho,double observer,double depth,
    double complex values[NEC4_GROUND_COMPONENT_COUNT])
{
  unsigned int component;
  sample_context *parameters=context;
  double bias=parameters->bias;
  if(parameters->reject_origin&&rho==0.&&observer==0.&&depth==0.) return 0;
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    values[component]=(bias+component+2.*rho-3.*observer+.5*depth)+
      I*(-bias+.25*component-rho+observer+4.*depth);
  return 1;
}

int main(void)
{
  const double frequency=7.1e6,wavelength=299792458./7.1e6;
  const double lg=9.9125363429561,lmax=wavelength;
  const double complex upper_k=.14880499655856940866;
  const double complex lower_k=.58721365385418336415-
    .23866676405070434672*I;
  sample_context parameters={1.25,0};
  nec4_ground_table serial={0},parallel={0},loaded={0};
  nec4_ground_table first={0},third={0};
  nec4_ground_table production_serial={0},production_parallel={0};
  double complex interpolated[NEC4_GROUND_COMPONENT_COUNT];
  double complex expected[NEC4_GROUND_COMPONENT_COUNT];
  char path[]="/tmp/nec2c-nec4-table-XXXXXX";
  int descriptor=mkstemp(path);
  unsigned int component;
  assert(descriptor>=0);close(descriptor);unlink(path);
  assert(nec4_ground_table_init(&serial,2,frequency,13.-.1*I,
      wavelength,lg,lmax));
  assert(serial.rho_count==35&&serial.observer_count==35&&
      serial.source_depth_count==41);
  assert(serial.value_count==35u*35u*41u*NEC4_GROUND_COMPONENT_COUNT);
  assert(serial.rho_step_m<=.03*lg&&serial.observer_step_m<=.03*lg&&
      serial.source_depth_step_m<=.015*lg);
  assert(nec4_ground_table_init(&first,1,frequency,13.-.1*I,
      wavelength,lg,lmax));
  assert(nec4_ground_table_init(&third,3,frequency,13.-.1*I,
      wavelength,lg,lmax));
  assert(first.rho_count==51&&first.observer_count==51&&
      first.source_depth_count==21);
  assert(third.rho_count==35&&third.observer_count==35&&
      third.source_depth_count==110);
  assert((first.value_count+serial.value_count+third.value_count)*
      sizeof(double complex)==19167680u);
  assert(nec4_ground_table_fill(&serial,affine_sample,&parameters,1));
  assert(nec4_ground_table_init(&parallel,2,frequency,13.-.1*I,
      wavelength,lg,lmax));
  assert(nec4_ground_table_fill(&parallel,affine_sample,&parameters,4));
  assert(!memcmp(serial.values,parallel.values,
      serial.value_count*sizeof(*serial.values)));
  assert(nec4_ground_table_interpolate(&serial,.317*lg,.689*lg,.731*lg,
      interpolated));
  assert(affine_sample(&parameters,.317*lg,.689*lg,.731*lg,expected));
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    assert(cabs(interpolated[component]-expected[component])<2.e-13*
        fmax(1.,cabs(expected[component])));
  assert(!nec4_ground_table_interpolate(&serial,1.01*lg,.2*lg,.7*lg,
      interpolated));
  parameters.reject_origin=1;
  assert(nec4_ground_table_fill(&first,affine_sample,&parameters,4));
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    assert(first.values[component]==0.);
  assert(nec4_ground_table_interpolate(&first,0.,0.,0.,interpolated));
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    assert(interpolated[component]==0.);
  parameters.reject_origin=0;
  assert(nec4_ground_table_store(&serial,path)==NEC4_CACHE_OK);
  assert(nec4_ground_table_init(&loaded,2,frequency,13.-.1*I,
      wavelength,lg,lmax));
  assert(nec4_ground_table_load(&loaded,path)==NEC4_CACHE_OK);
  assert(!memcmp(serial.values,loaded.values,
      serial.value_count*sizeof(*serial.values)));
  loaded.epsilon_relative=14.-.1*I;
  assert(nec4_ground_table_load(&loaded,path)==NEC4_CACHE_KEY_MISMATCH);
  {
    nec4_ground_interface_result interface;
    nec4_ground_table1_result direct;
    nec4_ground_direct_result outer_direct;
    nec4_ground_table outer[2]={{0},{0}};
    nec4_ground_table_direct_options invalid={0,48,80.};
    double complex reconstructed[NEC4_GROUND_COMPONENT_COUNT];
    double complex factor;
    size_t interface_point;
    double maximum=0.,scale=0.,interpolation_error,outer_error=0.;
    double step=.02*lg;
    assert(nec4_ground_table_init(&production_serial,1,frequency,
        13.-12.65852365107207*I,
        wavelength,lg,lmax));
    assert(nec4_ground_table_init(&production_parallel,1,frequency,
        13.-12.65852365107207*I,
        wavelength,lg,lmax));
    production_serial.rho_count=production_serial.observer_count=
      production_serial.source_depth_count=2;
    production_parallel.rho_count=production_parallel.observer_count=
      production_parallel.source_depth_count=2;
    production_serial.rho_step_m=production_serial.observer_step_m=
      production_serial.source_depth_step_m=step;
    production_parallel.rho_step_m=production_parallel.observer_step_m=
      production_parallel.source_depth_step_m=step;
    production_serial.value_count=production_parallel.value_count=
      2u*2u*2u*NEC4_GROUND_COMPONENT_COUNT;
    assert(!nec4_ground_table_fill_direct(&production_serial,lower_k,upper_k,
        &invalid,1));
    assert(!nec4_ground_table_fill_direct(&production_serial,1.01*lower_k,
        upper_k,NULL,1));
    assert(nec4_ground_table_fill_direct(&production_serial,lower_k,upper_k,
        NULL,1));
    assert(nec4_ground_table_fill_direct(&production_parallel,lower_k,upper_k,
        NULL,4));
    assert(!memcmp(production_serial.values,production_parallel.values,
        production_serial.value_count*sizeof(*production_serial.values)));
    for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
      assert(production_serial.values[component]==0.);
    assert(nec4_ground_direct_cross_interface(frequency,lower_k,upper_k,
        step,3,48,80.,&interface));
    interface_point=((1u*2u+0u)*2u+0u)*NEC4_GROUND_COMPONENT_COUNT;
    factor=step*cexp(I*upper_k*step);
    for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
      assert(cabs(production_serial.values[interface_point+component]-
          interface.field.finite_remainder[component]*factor)<2.e-13*
          fmax(1.,cabs(interface.field.finite_remainder[component]*factor)));
    assert(nec4_ground_table_evaluate(&production_serial,lower_k,upper_k,
        step,0.,0.,reconstructed));
    assert(nec4_ground_direct_table1_restore(frequency,lower_k,upper_k,
        step,0.,0.,interface.field.finite_remainder,expected));
    for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
      assert(cabs(reconstructed[component]-expected[component])<2.e-13*
          fmax(1.,cabs(expected[component])));
    assert(nec4_ground_direct_cross_table1(frequency,lower_k,upper_k,
        step,step,step,48,80.,&direct));
    assert(nec4_ground_table_evaluate(&production_serial,lower_k,upper_k,
        step,step,step,reconstructed));
    for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    {
      scale=fmax(scale,cabs(direct.total[component]));
      maximum=fmax(maximum,cabs(reconstructed[component]-
          direct.total[component]));
    }
    assert(maximum/fmax(1.,scale)<2.e-5);
    maximum=0.;scale=0.;
    assert(nec4_ground_direct_cross_table1(frequency,lower_k,upper_k,
        step/2.,step/2.,step/2.,48,80.,&direct));
    assert(nec4_ground_table_evaluate(&production_serial,lower_k,upper_k,
        step/2.,step/2.,step/2.,reconstructed));
    for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    {
      scale=fmax(scale,cabs(direct.total[component]));
      maximum=fmax(maximum,cabs(reconstructed[component]-
          direct.total[component]));
    }
    interpolation_error=maximum/fmax(1.,scale);
    assert(interpolation_error<5.e-3);
    for(component=0;component<2;++component)
    {
      unsigned int kind=component+2;
      unsigned int field_component;
      double rho,observer,depth;
      assert(nec4_ground_table_init(&outer[component],kind,frequency,
          13.-12.65852365107207*I,wavelength,lg,lmax));
      outer[component].rho_count=outer[component].observer_count=
        outer[component].source_depth_count=2;
      outer[component].value_count=2u*2u*2u*NEC4_GROUND_COMPONENT_COUNT;
      assert(nec4_ground_table_fill_direct(&outer[component],lower_k,upper_k,
          NULL,2));
      rho=outer[component].rho_step_m;
      observer=outer[component].observer_step_m;
      depth=outer[component].source_depth_start_m+
        outer[component].source_depth_step_m;
      assert(nec4_ground_direct_cross_bessel(frequency,lower_k,upper_k,rho,
          observer,depth,48,80.,&outer_direct));
      assert(nec4_ground_table_evaluate(&outer[component],lower_k,upper_k,
          rho,observer,depth,reconstructed));
      scale=0.;maximum=0.;
      for(field_component=0;field_component<NEC4_GROUND_COMPONENT_COUNT;
          ++field_component)
      {
        scale=fmax(scale,cabs(outer_direct.value[field_component]));
        maximum=fmax(maximum,cabs(reconstructed[field_component]-
            outer_direct.value[field_component]));
      }
      outer_error=fmax(outer_error,maximum/fmax(1.,scale));
    }
    assert(outer_error<2.e-5);
    assert(!nec4_ground_table_evaluate(&production_serial,lower_k,upper_k,
        0.,0.,0.,reconstructed));
    printf("NEC4 Appendix-H table generation passed: serial/parallel exact, "
        "cell-center physical error %.3e, F2/F3 %.3e\n",
        interpolation_error,outer_error);
    nec4_ground_table_destroy(&outer[1]);
    nec4_ground_table_destroy(&outer[0]);
  }
  unlink(path);
  nec4_ground_table_destroy(&production_parallel);
  nec4_ground_table_destroy(&production_serial);
  nec4_ground_table_destroy(&third);nec4_ground_table_destroy(&first);
  nec4_ground_table_destroy(&loaded);nec4_ground_table_destroy(&parallel);
  nec4_ground_table_destroy(&serial);
  return 0;
}
