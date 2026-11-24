#include "nec4_ground_matrix.h"
#include "nec4_ground_direct.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define C0 299792458.0
#define MU0 1.25663706212e-6

typedef struct { double frequency; double complex k; } cross_context;

typedef struct { double frequency; double complex lower_k,upper_k; } material_context;

static int identical_cross(void *opaque,double rho,double observer_z,
    double source_depth,double complex value[5])
{
  cross_context *context=opaque;
  nec4_ground_direct_result result;
  if(!nec4_ground_direct_cross_bessel(context->frequency,context->k,context->k,
      rho,observer_z,source_depth,64,75.,&result)) return 0;
  for(int i=0;i<5;++i) value[i]=result.value[i];
  return 1;
}

static int material_cross(void *opaque,double rho,double observer_z,
    double source_depth,double complex value[5])
{
  material_context *context=opaque;
  nec4_ground_direct_result result;
  if(!nec4_ground_direct_cross_bessel(context->frequency,context->lower_k,
      context->upper_k,rho,observer_z,source_depth,96,90.,&result)) return 0;
  for(int i=0;i<5;++i) value[i]=result.value[i];
  return 1;
}

static void homogeneous_cartesian(double frequency,double complex k,
    const double source[3],const double observer[3],const double direction[3],
    double complex field[3])
{
  double dx=observer[0]-source[0],dy=observer[1]-source[1];
  double rho=hypot(dx,dy),dz=observer[2]-source[2],r=hypot(rho,dz);
  double c=rho?dx/rho:1.,s=rho?dy/rho:0.,omega=2.*acos(-1.)*frequency;
  double complex common=-I*omega*MU0/(4.*acos(-1.)*k*k)*cexp(-I*k*r)/(r*r*r);
  double complex radial=3.+3.*I*k*r-k*k*r*r;
  double complex transverse=1.+I*k*r-k*k*r*r;
  double complex rv=common*rho*dz/(r*r)*radial;
  double complex zv=common*(dz*dz/(r*r)*radial-transverse);
  double complex rh=common*(rho*rho/(r*r)*radial-transverse);
  double complex ph=common*transverse,zh=rv;
  double complex xx=c*c*rh-s*s*ph,xy=s*c*(rh+ph),yy=s*s*rh-c*c*ph;
  field[0]=xx*direction[0]+xy*direction[1]+c*rv*direction[2];
  field[1]=xy*direction[0]+yy*direction[1]+s*rv*direction[2];
  field[2]=c*zh*direction[0]+s*zh*direction[1]+zv*direction[2];
}

static double relative(double complex a,double complex b)
{ return cabs(a-b)/fmax(1.,fmax(cabs(a),cabs(b))); }

static int synthetic_h(void *opaque,double rho,double z,double depth,
    double complex field[NEC4_GROUND_H_COMPONENT_COUNT])
{
  (void)opaque;(void)rho;(void)z;(void)depth;
  field[NEC4_GROUND_H_PHI_VERTICAL]=2.+3.*I;
  field[NEC4_GROUND_H_RHO_HORIZONTAL]=5.-I;
  field[NEC4_GROUND_H_PHI_HORIZONTAL]=-2.+4.*I;
  field[NEC4_GROUND_H_Z_HORIZONTAL]=1.-6.*I;
  return 1;
}

int main(void)
{
  const double frequency=12.e6,k=2.*acos(-1.)*frequency/C0;
  cross_context cross={frequency,k};
  nec4_ground_matrix_context context={frequency,1.,k,k,identical_cross,NULL,
    &cross};
  const double direction[3]={.36,-.48,.8};
  const struct {
    nec4_segment_medium source_medium,observer_medium;
    double source_z,observer_z;
  } pairs[]={
    {NEC4_MEDIUM_LOWER,NEC4_MEDIUM_UPPER,-.7,.4},
    {NEC4_MEDIUM_UPPER,NEC4_MEDIUM_LOWER,.7,-.4},
    {NEC4_MEDIUM_LOWER,NEC4_MEDIUM_LOWER,-.7,-.4},
    {NEC4_MEDIUM_UPPER,NEC4_MEDIUM_UPPER,.7,.4}};
  double worst=0.;
  for(size_t p=0;p<sizeof(pairs)/sizeof(pairs[0]);++p)
  {
    double source[3]={-.2,.15,pairs[p].source_z};
    double observer[3]={.31,-.27,pairs[p].observer_z};
    double complex actual[3],expected[3];
    assert(nec4_ground_matrix_point(&context,pairs[p].source_medium,
      pairs[p].observer_medium,source,observer,direction,actual));
    homogeneous_cartesian(frequency,k,source,observer,direction,expected);
    for(int component=0;component<3;++component)
      worst=fmax(worst,relative(actual[component],expected[component]));
  }
  assert(worst<2.e-5);

  /* The independently qualified cylindrical H kernels must retain their
     source-azimuth signs when assembled into Cartesian segment fields. */
  {
    nec4_ground_matrix_context hcontext=context;
    double source[3]={0.,0.,-.7},observer[3]={3.,4.,.4};
    double hdirection[3]={.6,.0,.8},cosine=.6,sine=.8;
    double source_rho=cosine*hdirection[0]+sine*hdirection[1];
    double source_phi=-sine*hdirection[0]+cosine*hdirection[1];
    double complex radial=-(5.-I)*source_phi;
    double complex azimuth=(-2.+4.*I)*source_rho+(2.+3.*I)*hdirection[2];
    double complex actual[3],expected[3]={
      cosine*radial-sine*azimuth,
      sine*radial+cosine*azimuth,
      -(1.-6.*I)*source_phi};
    hcontext.cross_h=synthetic_h;
    assert(nec4_ground_matrix_h_point_cross(&hcontext,source,observer,
      hdirection,actual));
    assert(!memcmp(actual,expected,sizeof(actual)));
  }


  /* Independent high-precision Python/M8 fixtures for epsilon_r=4 exercise
     the non-null Appendix-G transformations, not merely their vacuum limit. */
  {
    static const double complex expected[4][5]={
      {-.15408949341335565-7.117776069121258*I,
       -.8912851690329092+.4196212960317558*I,
       -1.4863304874229102-8.305192153911367*I,
       1.4488182560449283-3.8698898584893655*I,
       .31072345356017245-5.932370822839459*I},
      {-.31156543208767756+5.7156456458229545*I,
       -.7639230434130813+.5507294620424187*I,
       -1.412940100972349-8.261715902469932*I,
       1.3826263846396094-3.6924917069999315*I,
       .12090499823173552+6.876752970463969*I},
      {-.15769130055838282-5.112637709755135*I,
       -.29977069446708704+3.2670180814579712*I,
       -1.6466202781400583-13.55757212200648*I,
       1.596680560627005-5.174792773344535*I,
       .12253395039650046+.7791427212436788*I},
      {-.31577599538275125-2.564533789640512*I,
       -2.4352340313452068+15.48654202751441*I,
       -1.3217210293956012-26.638971350641445*I,
       1.299871237074812-8.531801932844497*I,
       .3203444179356354+18.26160732745499*I}};
    double k0=2.*acos(-1.)*10.e6/C0;
    material_context material={10.e6,2.*k0,k0};
    nec4_ground_matrix_context nonnull={10.e6,4.,2.*k0,k0,
      material_cross,NULL,&material};
    const double axes[3][3]={{0.,0.,1.},{1.,0.,0.},{0.,1.,0.}};
    double transformed_worst=0.;
    for(size_t p=0;p<4;++p)
    {
      double source[3]={0.,0.,pairs[p].source_z};
      double observer[3]={2.,0.,pairs[p].observer_z};
      double complex field[3];
      assert(nec4_ground_matrix_point(&nonnull,pairs[p].source_medium,
        pairs[p].observer_medium,source,observer,axes[0],field));
      transformed_worst=fmax(transformed_worst,
        cabs(field[0]-expected[p][0])/fmax(1.,cabs(expected[p][0])));
      transformed_worst=fmax(transformed_worst,
        cabs(field[2]-expected[p][1])/fmax(1.,cabs(expected[p][1])));
      assert(nec4_ground_matrix_point(&nonnull,pairs[p].source_medium,
        pairs[p].observer_medium,source,observer,axes[1],field));
      transformed_worst=fmax(transformed_worst,
        cabs(field[0]-expected[p][2])/fmax(1.,cabs(expected[p][2])));
      transformed_worst=fmax(transformed_worst,
        cabs(field[2]-expected[p][4])/fmax(1.,cabs(expected[p][4])));
      assert(nec4_ground_matrix_point(&nonnull,pairs[p].source_medium,
        pairs[p].observer_medium,source,observer,axes[2],field));
      transformed_worst=fmax(transformed_worst,
        cabs(field[1]+expected[p][3])/fmax(1.,cabs(expected[p][3])));
    }
    printf("NEC4 nonidentical 11/12/21/22 matrix fixtures passed; worst %.3e\n",
      transformed_worst);
    assert(transformed_worst<1.e-4);

    /* A whole segment and two half-segments provide an independent quadrature
       refinement while preserving the parent sine/cosine basis exactly. */
    {
      double center[3]={0.,0.,-.7},observer[3]={2.,0.,.4};
      double left_center[3],right_center[3],half=.04;
      double complex phase;
      nec4_ground_segment_field whole,left,right;
      double quadrature_worst=0.;
      for(int axis=0;axis<3;++axis)
      {
        left_center[axis]=center[axis]-.25*.08*direction[axis];
        right_center[axis]=center[axis]+.25*.08*direction[axis];
      }
      assert(nec4_ground_matrix_segment_correction(&nonnull,
        NEC4_MEDIUM_LOWER,NEC4_MEDIUM_UPPER,center,direction,.08,observer,
        .001,&whole));
      assert(nec4_ground_matrix_segment_correction(&nonnull,
        NEC4_MEDIUM_LOWER,NEC4_MEDIUM_UPPER,left_center,direction,half,
        observer,.001,&left));
      assert(nec4_ground_matrix_segment_correction(&nonnull,
        NEC4_MEDIUM_LOWER,NEC4_MEDIUM_UPPER,right_center,direction,half,
        observer,.001,&right));
      phase=nonnull.lower_wave_number_per_m*.02;
      for(int component=0;component<3;++component)
      {
        double complex refined_constant=left.constant[component]+
          right.constant[component];
        double complex refined_sine=
          left.sine[component]*ccos(phase)-left.cosine[component]*csin(phase)+
          right.sine[component]*ccos(phase)+right.cosine[component]*csin(phase);
        double complex refined_cosine=
          left.cosine[component]*ccos(phase)+left.sine[component]*csin(phase)+
          right.cosine[component]*ccos(phase)-right.sine[component]*csin(phase);
        quadrature_worst=fmax(quadrature_worst,
          cabs(whole.constant[component]-refined_constant)/
          fmax(1.,cabs(refined_constant)));
        quadrature_worst=fmax(quadrature_worst,
          cabs(whole.sine[component]-refined_sine)/fmax(1.,cabs(refined_sine)));
        quadrature_worst=fmax(quadrature_worst,
          cabs(whole.cosine[component]-refined_cosine)/
          fmax(1.,cabs(refined_cosine)));
      }
      assert(quadrature_worst<1.e-8);
      printf("NEC4 matrix segment quadrature refinement passed; worst %.3e\n",
        quadrature_worst);
    }

    /* Maxwell similarity is an absolute units gate, not a comparison of two
       evaluators sharing an adapter.  Doubling frequency and halving every
       length doubles physical E/I (ohm/m), while the NEC matrix coefficient
       in ohms remains invariant after multiplication by wavelength. */
    {
      material_context scaled_material={20.e6,4.*k0,2.*k0};
      nec4_ground_matrix_context scaled={20.e6,4.,4.*k0,2.*k0,
        material_cross,NULL,&scaled_material};
      double center[3]={0.,0.,-.7},observer[3]={2.,0.,.4};
      double scaled_center[3]={0.,0.,-.35},scaled_observer[3]={1.,0.,.2};
      nec4_ground_segment_field physical,physical_scaled,matrix_value,
        matrix_scaled;
      double physical_worst=0.,matrix_worst=0.;
      assert(nec4_ground_matrix_segment_correction(&nonnull,
        NEC4_MEDIUM_LOWER,NEC4_MEDIUM_UPPER,center,direction,.08,observer,
        .001,&physical));
      assert(nec4_ground_matrix_segment_correction(&scaled,
        NEC4_MEDIUM_LOWER,NEC4_MEDIUM_UPPER,scaled_center,direction,.04,
        scaled_observer,.0005,&physical_scaled));
      assert(nec4_ground_matrix_segment_nec_correction(&nonnull,
        NEC4_MEDIUM_LOWER,NEC4_MEDIUM_UPPER,center,direction,.08,observer,
        .001,C0/10.e6,&matrix_value));
      assert(nec4_ground_matrix_segment_nec_correction(&scaled,
        NEC4_MEDIUM_LOWER,NEC4_MEDIUM_UPPER,scaled_center,direction,.04,
        scaled_observer,.0005,C0/20.e6,&matrix_scaled));
      for(int component=0;component<3;++component)
      {
        const double complex *a[3]={physical.constant,physical.sine,
          physical.cosine};
        const double complex *b[3]={physical_scaled.constant,
          physical_scaled.sine,physical_scaled.cosine};
        const double complex *m[3]={matrix_value.constant,matrix_value.sine,
          matrix_value.cosine};
        const double complex *n[3]={matrix_scaled.constant,
          matrix_scaled.sine,matrix_scaled.cosine};
        for(int basis=0;basis<3;++basis)
        {
          physical_worst=fmax(physical_worst,
            cabs(b[basis][component]-2.*a[basis][component])/
            fmax(1.,cabs(2.*a[basis][component])));
          matrix_worst=fmax(matrix_worst,
            cabs(n[basis][component]-m[basis][component])/
            fmax(1.,cabs(m[basis][component])));
        }
      }
      assert(physical_worst<2.e-5);
      assert(matrix_worst<2.e-5);
      printf("NEC4 SI/NEC normalization similarity passed; physical %.3e, "
        "matrix %.3e\n",physical_worst,matrix_worst);
    }

    /* The interface-junction endpoint charges obey the same Maxwell scaling,
       and opposite ends must carry exactly opposite signs under exp(+jwt). */
    {
      nec4_ground_matrix_context scaled=nonnull;
      double endpoint[3]={0.,0.,-.2},observer[3]={.3,-.1,.4};
      double scaled_endpoint[3]={0.,0.,-.1},
        scaled_observer[3]={.15,-.05,.2};
      double complex positive[3],negative[3],scaled_field[3];
      double charge_worst=0.;
      scaled.frequency_hz*=2.;
      scaled.lower_wave_number_per_m*=2.;
      scaled.upper_wave_number_per_m*=2.;
      assert(nec4_ground_matrix_endpoint_charge_nec(&nonnull,
        NEC4_MEDIUM_LOWER,endpoint,observer,C0/10.e6,1,positive));
      assert(nec4_ground_matrix_endpoint_charge_nec(&nonnull,
        NEC4_MEDIUM_LOWER,endpoint,observer,C0/10.e6,-1,negative));
      assert(nec4_ground_matrix_endpoint_charge_nec(&scaled,
        NEC4_MEDIUM_LOWER,scaled_endpoint,scaled_observer,C0/20.e6,1,
        scaled_field));
      for(int component=0;component<3;++component)
      {
        charge_worst=fmax(charge_worst,relative(positive[component],
          scaled_field[component]));
        assert(cabs(positive[component]+negative[component])<1.e-13*
          fmax(1.,cabs(positive[component])));
      }
      assert(charge_worst<2.e-12);
      printf("NEC4 interface endpoint-charge sign and similarity passed; "
        "worst %.3e\n",charge_worst);
    }
  }

  /* In identical media the complete-interface minus free-space segment
     correction must vanish for every pairing and all three current bases. */
  for(size_t p=0;p<sizeof(pairs)/sizeof(pairs[0]);++p)
  {
    double source[3]={-.2,.15,pairs[p].source_z};
    double observer[3]={.31,-.27,pairs[p].observer_z};
    nec4_ground_segment_field correction;
    assert(nec4_ground_matrix_segment_correction(&context,
      pairs[p].source_medium,pairs[p].observer_medium,source,direction,.08,
      observer,.001,&correction));
    for(int component=0;component<3;++component)
    {
      assert(cabs(correction.constant[component])<2.e-5);
      assert(cabs(correction.sine[component])<2.e-5);
      assert(cabs(correction.cosine[component])<2.e-5);
    }
  }
  printf("NEC4 ground matrix pairings and segment-basis integration passed; "
      "identical-medium worst %.3e\n",worst);
  return 0;
}
