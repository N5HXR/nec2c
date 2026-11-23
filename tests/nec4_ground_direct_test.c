#include "nec4_ground_direct.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct {
  const char *name;
  double complex lower_k;
  double rho;
  double observer_z;
  double source_depth;
  double complex expected[NEC4_GROUND_COMPONENT_COUNT];
} witness;

static double maximum_error(const nec4_ground_direct_result *actual,
    const witness *reference)
{
  double scale=0.,maximum=0.;
  unsigned int component;
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    scale=fmax(scale,cabs(reference->expected[component]));
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
  {
    double denominator=fmax(cabs(reference->expected[component]),scale*1.e-12);
    maximum=fmax(maximum,cabs(actual->value[component]-
        reference->expected[component])/denominator);
  }
  return maximum;
}

static double maximum_array_error(
    const double complex actual[NEC4_GROUND_COMPONENT_COUNT],
    const double complex expected[NEC4_GROUND_COMPONENT_COUNT])
{
  double scale=0.,maximum=0.;
  unsigned int component;
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    scale=fmax(scale,cabs(expected[component]));
  for(component=0;component<NEC4_GROUND_COMPONENT_COUNT;++component)
    maximum=fmax(maximum,cabs(actual[component]-expected[component])/
        fmax(cabs(expected[component]),scale*1.e-12));
  return maximum;
}

int main(void)
{
  const double frequency=7.1e6;
  const double complex upper_k=0.14880499655856940866;
  const witness cases[]={
    {"lossless-er4",0.29760999311713881732,1.,.4,.7,
      {-.026288065897441172657-37.684652994645370107*I,
       -.45688637372782125766-17.683327630222629498*I,
       -.81142729702230278259-10.986050846749921821*I,
       .80915185810104947979-23.001470850315020365*I,
       .063478106095292008402-36.907946965375622067*I}},
    {"average-soil",.58721365385418336415-.23866676405070434672*I,
      1.,.4,.7,
      {6.3836872273565708008-8.232442339051847326*I,
       2.4744353926009315379-4.3428211443646748076*I,
       .29850819059099813487-3.3313249684161377555*I,
       5.5977242335931650525-3.8470392821995311813*I,
       6.8468407816011925787-7.1785526144168127738*I}},
    {"seawater",10.631043272109062923-10.546350267575551521*I,
      .35,.15,.25,
      {-.074037438028091987566-.041546318341716982903*I,
       -.009496663520618016771+.022118442129944226030*I,
       .020275561677981072028+.11545071378473604195*I,
       -.097193487210066490449-.19556672627807042025*I,
       -.064229316680052791582-.25957907602937218603*I}}
  };
  size_t index;
  for(index=0;index<sizeof(cases)/sizeof(cases[0]);++index)
  {
    nec4_ground_direct_result result;
    double error;
    assert(nec4_ground_direct_cross_bessel(frequency,cases[index].lower_k,
        upper_k,cases[index].rho,cases[index].observer_z,
        cases[index].source_depth,44,70.,&result));
    error=maximum_error(&result,&cases[index]);
    if(error>8.e-6)
    {
      fprintf(stderr,"%s direct cross-interface error %.12e\n",
          cases[index].name,error);
      return 1;
    }
  }
  {
    const double complex expected_total[NEC4_GROUND_COMPONENT_COUNT]={
      .77778784730762230589-1.3415792910956182931*I,
      .67312962899974722754-1.8028989234627023687*I,
      -1.0689432124699494195+.046992286802651275372*I,
      1.4922662848890811749-.67848981828712806965*I,
      1.0541553844195599744-1.0825234293645982335*I};
    const double complex expected_finite[NEC4_GROUND_COMPONENT_COUNT]={
      .000097748872008374327028-.039003562559376554750*I,
      -.17058118121800044279+.11224880850360378614*I,
      -.20434533103108121099+.9253499857319079*I,
      .220494543488437-.948509053659522*I,
      .11018197730245104-.10228295158412438*I};
    nec4_ground_table1_result result;
    double total_error,finite_error;
    assert(nec4_ground_direct_cross_table1(frequency,cases[1].lower_k,
        upper_k,1.2,1.2,1.2,44,70.,&result));
    total_error=maximum_array_error(result.total,expected_total);
    finite_error=maximum_array_error(result.finite_remainder,expected_finite);
    if(total_error>8.e-6||finite_error>8.e-6)
    {
      fprintf(stderr,"table-1 direct errors total %.12e finite %.12e\n",
          total_error,finite_error);
      return 1;
    }
  }
  {
    nec4_ground_direct_result axis,coarse,fine;
    nec4_ground_table1_result table_axis,table_coarse,table_fine;
    const double boundary_axis[2][2]={{0.,.7},{.4,0.}};
    double coarse_error,fine_error,minimum_quadratic_ratio=1.e100;
    double boundary_reconstruction_error=0.;
    unsigned int component,boundary_case;
    assert(nec4_ground_direct_cross_bessel(frequency,cases[1].lower_k,
        upper_k,0.,.4,.7,44,70.,&axis));
    assert(nec4_ground_direct_cross_table1(frequency,cases[1].lower_k,
        upper_k,0.,.4,.7,44,70.,&table_axis));
    assert(axis.value[NEC4_GROUND_RHO_VERTICAL]==0.);
    assert(axis.value[NEC4_GROUND_Z_HORIZONTAL]==0.);
    assert(cabs(axis.value[NEC4_GROUND_RHO_HORIZONTAL]+
        axis.value[NEC4_GROUND_PHI_HORIZONTAL])<2.e-14*
        cabs(axis.value[NEC4_GROUND_RHO_HORIZONTAL]));
    assert(maximum_array_error(table_axis.total,axis.value)<5.e-14);
    assert(table_axis.finite_remainder[NEC4_GROUND_RHO_VERTICAL]==0.);
    assert(table_axis.finite_remainder[NEC4_GROUND_Z_HORIZONTAL]==0.);
    assert(cabs(table_axis.finite_remainder[NEC4_GROUND_RHO_HORIZONTAL]+
        table_axis.finite_remainder[NEC4_GROUND_PHI_HORIZONTAL])<2.e-14*
        cabs(table_axis.finite_remainder[NEC4_GROUND_RHO_HORIZONTAL]));
    assert(nec4_ground_direct_cross_bessel(frequency,cases[1].lower_k,
        upper_k,1.e-3,.4,.7,44,70.,&coarse));
    assert(nec4_ground_direct_cross_bessel(frequency,cases[1].lower_k,
        upper_k,5.e-4,.4,.7,44,70.,&fine));
    assert(nec4_ground_direct_cross_table1(frequency,cases[1].lower_k,
        upper_k,1.e-3,.4,.7,44,70.,&table_coarse));
    assert(nec4_ground_direct_cross_table1(frequency,cases[1].lower_k,
        upper_k,5.e-4,.4,.7,44,70.,&table_fine));
    assert(cabs(fine.value[NEC4_GROUND_RHO_VERTICAL])/
        cabs(coarse.value[NEC4_GROUND_RHO_VERTICAL])>.4999);
    assert(cabs(fine.value[NEC4_GROUND_RHO_VERTICAL])/
        cabs(coarse.value[NEC4_GROUND_RHO_VERTICAL])<.5001);
    assert(cabs(fine.value[NEC4_GROUND_Z_HORIZONTAL])/
        cabs(coarse.value[NEC4_GROUND_Z_HORIZONTAL])>.4999);
    assert(cabs(fine.value[NEC4_GROUND_Z_HORIZONTAL])/
        cabs(coarse.value[NEC4_GROUND_Z_HORIZONTAL])<.5001);
    for(component=NEC4_GROUND_Z_VERTICAL;
        component<=NEC4_GROUND_PHI_HORIZONTAL;++component)
    {
      coarse_error=cabs(coarse.value[component]-axis.value[component]);
      fine_error=cabs(fine.value[component]-axis.value[component]);
      assert(coarse_error/fine_error>3.99);
      assert(coarse_error/fine_error<4.01);
      minimum_quadratic_ratio=fmin(minimum_quadratic_ratio,
                                  coarse_error/fine_error);
      coarse_error=cabs(table_coarse.finite_remainder[component]-
                        table_axis.finite_remainder[component]);
      fine_error=cabs(table_fine.finite_remainder[component]-
                      table_axis.finite_remainder[component]);
      assert(coarse_error/fine_error>3.95);
      assert(coarse_error/fine_error<4.05);
    }
    for(boundary_case=0;boundary_case<2;++boundary_case)
    {
      assert(nec4_ground_direct_cross_bessel(frequency,cases[1].lower_k,
          upper_k,0.,boundary_axis[boundary_case][0],
          boundary_axis[boundary_case][1],44,70.,&axis));
      assert(nec4_ground_direct_cross_table1(frequency,cases[1].lower_k,
          upper_k,0.,boundary_axis[boundary_case][0],
          boundary_axis[boundary_case][1],44,70.,&table_axis));
      assert(axis.value[NEC4_GROUND_RHO_VERTICAL]==0.);
      assert(axis.value[NEC4_GROUND_Z_HORIZONTAL]==0.);
      assert(cabs(axis.value[NEC4_GROUND_RHO_HORIZONTAL]+
          axis.value[NEC4_GROUND_PHI_HORIZONTAL])<2.e-13*
          fmax(1.,cabs(axis.value[NEC4_GROUND_RHO_HORIZONTAL])));
      coarse_error=maximum_array_error(table_axis.total,axis.value);
      boundary_reconstruction_error=fmax(boundary_reconstruction_error,
                                         coarse_error);
      assert(coarse_error<1.e-9);
    }
    assert(!nec4_ground_direct_cross_bessel(frequency,cases[1].lower_k,
        upper_k,0.,0.,0.,44,70.,&axis));
    assert(!nec4_ground_direct_cross_table1(frequency,cases[1].lower_k,
        upper_k,0.,0.,0.,44,70.,&table_axis));
    printf("NEC4 cylindrical-axis limit passed: radial 0.500x, "
        "non-radial >= %.3fx quadratic convergence, boundary %.3e\n",
        minimum_quadratic_ratio,boundary_reconstruction_error);
    assert(!nec4_ground_direct_cross_bessel(frequency,1.,upper_k,1.,.2,.3,
        8,50.,&axis));
  }
  {
    nec4_ground_direct_result serial,parallel[16];
    int item,success[16];
    assert(nec4_ground_direct_cross_bessel(frequency,cases[1].lower_k,
        upper_k,1.,.4,.7,32,50.,&serial));
#ifdef _OPENMP
#pragma omp parallel for num_threads(4) schedule(static)
#endif
    for(item=0;item<16;++item)
      success[item]=nec4_ground_direct_cross_bessel(
          frequency,cases[1].lower_k,upper_k,1.,.4,.7,32,50.,
          &parallel[item]);
    for(item=0;item<16;++item)
    {
      assert(success[item]);
      assert(!memcmp(serial.value,parallel[item].value,
          sizeof(serial.value)));
    }
  }
  {
    const double wavelength=299792458./frequency;
    const double complex interface_expected[NEC4_GROUND_COMPONENT_COUNT]={
      -.00079950159205886048-.00015938449501407253*I,
      -.0010489299281001693+.000012068155037928520*I,
      .087596274034749011+10.642573865457006*I,
      -.088170565849976623-10.642811291502829*I,
      1.692707070847741-8.0894275557541047*I};
    nec4_ground_direct_result bessel_result,hankel_result;
    nec4_ground_interface_result serial,parallel[8],point;
    nec4_ground_table1_result table1_result;
    double worst_refinement=0.,worst_path=0.,overlap_error,witness_error;
    int item,success[8];
    unsigned int material,radial_index;
    assert(!nec4_ground_direct_cross_hankel(frequency,cases[1].lower_k,
        upper_k,5.,0.,0.,3,48,80.,&hankel_result));
    assert(!nec4_ground_direct_cross_table1_hankel(frequency,
        cases[1].lower_k,upper_k,5.,0.,0.,3,48,80.,
        &table1_result));
    assert(!nec4_ground_direct_cross_interface(frequency,cases[1].lower_k,
        upper_k,0.,3,48,80.,&serial));
    for(material=0;material<sizeof(cases)/sizeof(cases[0]);++material)
    {
      double lg=fmin(cabs(upper_k/cases[material].lower_k),.3)*wavelength;
      for(radial_index=1;radial_index<=50;++radial_index)
      {
        double rho=lg*radial_index/50.;
        assert(nec4_ground_direct_cross_interface(frequency,
            cases[material].lower_k,upper_k,rho,3,48,80.,&point));
        worst_refinement=fmax(worst_refinement,point.refinement_error);
        worst_path=fmax(worst_path,point.approach_path_error);
      }
    }
    assert(worst_refinement<1.e-5&&worst_path<1.e-5);
    assert(nec4_ground_direct_cross_interface(frequency,cases[2].lower_k,
        upper_k,fmin(cabs(upper_k/cases[2].lower_k),.3)*wavelength,
        3,48,80.,&point));
    witness_error=maximum_array_error(point.field.finite_remainder,
                                      interface_expected);
    assert(witness_error<2.e-10);
    assert(nec4_ground_direct_cross_bessel(frequency,cases[1].lower_k,
        upper_k,1.,.3,.25,48,80.,&bessel_result));
    assert(nec4_ground_direct_cross_hankel(frequency,cases[1].lower_k,
        upper_k,1.,.3,.25,3,48,80.,&hankel_result));
    overlap_error=maximum_array_error(hankel_result.value,bessel_result.value);
    assert(overlap_error<8.e-5);
    assert(nec4_ground_direct_cross_interface(frequency,cases[1].lower_k,
        upper_k,5.,3,48,80.,&serial));
#ifdef _OPENMP
#pragma omp parallel for num_threads(4) schedule(static)
#endif
    for(item=0;item<8;++item)
      success[item]=nec4_ground_direct_cross_interface(frequency,
          cases[1].lower_k,upper_k,5.,3,48,80.,&parallel[item]);
    for(item=0;item<8;++item)
    {
      assert(success[item]);
      assert(!memcmp(&serial,&parallel[item],sizeof(serial)));
    }
    printf("NEC4 interface plane passed: 150 radii, refinement %.3e, "
        "approach path %.3e, overlap %.3e\n",worst_refinement,worst_path,
        overlap_error);
  }
  {
    const double complex media[]={
      cases[1].lower_k,cases[2].lower_k,9.*upper_k
    };
    double worst_recomposition=0.,worst_refinement=0.;
    for(unsigned int material=0;material<3;++material)
      for(unsigned int index=0;index<12;++index) {
        double rho=(1.5+.3*index)*299792458./frequency;
        double z=(.01+.015*(index%4))*299792458./frequency;
        double depth=(.002+.003*(index%3))*299792458./frequency;
        nec4_ground_direct_split_result coarse,fine;
        assert(nec4_ground_direct_cross_hankel_split(frequency,
            media[material],upper_k,rho,z,depth,3,48,80.,&coarse));
        assert(nec4_ground_direct_cross_hankel_split(frequency,
            media[material],upper_k,rho,z,depth,5,80,120.,&fine));
        for(unsigned int component=0;
            component<NEC4_GROUND_COMPONENT_COUNT;++component) {
          double scale=fmax(cabs(fine.value[component]),1.e-300);
          worst_recomposition=fmax(worst_recomposition,cabs(
              coarse.value[component]-coarse.lower_spectral[component]-
              coarse.upper_spectral[component])/scale);
          worst_refinement=fmax(worst_refinement,cabs(
              coarse.value[component]-fine.value[component])/scale);
        }
      }
    printf("NEC4 exact spectral split passed: recomposition %.3e, "
        "refinement %.3e\n",worst_recomposition,worst_refinement);
    assert(worst_recomposition<3.e-13);
    assert(worst_refinement<3.e-4);
  }
  return 0;
}
