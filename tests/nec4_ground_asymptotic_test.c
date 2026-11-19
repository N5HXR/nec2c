#include "nec4_ground_asymptotic.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>

typedef struct {
  double complex k1;
  double depth;
  double complex first[5];
  double complex second[5];
  double complex pole_first[5];
  double complex pole_second[5];
} fixture;

static double relative_error(double complex actual,double complex expected)
{
  return cabs(actual-expected)/fmax(1.e-14,cabs(expected));
}

int main(void)
{
  const double frequency=7.1e6,k2=0.1488049965585694;
  const double rho=118.2280116056338,z=76.00372174647887;
  const fixture cases[]={
    {0.5952199862342776,7.917054348591549,
     {0.0004763744253449537+0.0011611185903941238*I,
      -9.103607597429004e-05-0.00022189201306039376*I,
      -0.0022245229565704895-0.00542206890674732*I,
      0.00307573709635988+0.007152785875180237*I,
      0.0034107467075806626+0.008313379557329775*I},
     {0.0006656770149792339+0.0010813099145580954*I,
      -0.000863878138601388-0.0017268547699209058*I,
      -0.0023810598134609445-0.0053664469425168905*I,
      0.0034776656114863485+0.006972892297880693*I,
      0.004247741286666199+0.007958438425393817*I},
     {0.00020371368303171656+5.418721665620934e-05*I,
      -0.000814854732126866-0.00021674886662483727*I,
      -0.000814854732126866-0.00021674886662483727*I,
      -1.269940591486782e-05+4.7742676425806783e-05*I,
      0.003259418928507465+0.0008669954664993494*I},
     {9.926124454828491e-05-0.00018238078119904916*I,
      -0.000418150759898448+0.0007174278862075225*I,
      0.0003759391964878313-0.0007416183633848708*I,
      -4.274306219144056e-05-2.3263029805188984e-05*I,
      0.0015881799127725582-0.002918092499184786*I}},
    {0.5872136538541833-0.23866676405070436*I,7.434402257217075,
     {-0.0001057613548443403+0.0001536472555787124*I,
      2.0096775600568842e-05-2.9599227234794992e-05*I,
      0.00018974904256366067-0.0009152451650157659*I,
      -0.00033514371739502613+0.0011632593650823812*I,
      -0.00028481999372874356+0.0014086854962830125*I},
     {-8.081165007480084e-05+0.00017114024470336975*I,
      0.00014648372168388825-0.0002521970918865421*I,
      0.00015846705155494217-0.000912759098408574*I,
      -0.00027060698492129486+0.0011843213929498757*I,
      -0.00014663636262255493+0.0014504633017966954*I},
     {-1.6011230511901247e-05+1.1364588301156473e-05*I,
      4.4955907474226303e-05-7.052713442646556e-05*I,
      4.4955907474226303e-05-7.052713442646556e-05*I,
      -4.039740015526783e-06-2.681321843466052e-06*I,
      -6.428708685983005e-05+0.00035041818803270287*I},
     {3.2987586112878586e-07+1.5957450893128545e-05*I,
      -2.5057085488965303e-05-6.314790668170218e-05*I,
      2.873439825758305e-05+6.173637618149448e-05*I,
      3.648109636412045e-06-1.4920272466123305e-06*I,
      0.00020628615573616325+0.0002032711202206553*I}}
  };
  double worst_first=0.,worst_second=0.,worst_pole_first=0.,worst_pole_second=0.;
  unsigned int c,n;
  for(c=0;c<sizeof(cases)/sizeof(cases[0]);c++) {
    nec4_ground_asymptotic_result first,second,pole_first,pole_second;
    assert(nec4_ground_asymptotic_first(frequency,cases[c].k1,k2,rho,z,
                                        cases[c].depth,&first));
    assert(nec4_ground_asymptotic_second_saddle(frequency,cases[c].k1,k2,
                                                rho,z,cases[c].depth,&second));
    assert(nec4_ground_asymptotic_first_pole(frequency,cases[c].k1,k2,rho,z,
                                             cases[c].depth,&pole_first));
    assert(nec4_ground_asymptotic_second_pole(frequency,cases[c].k1,k2,rho,z,
                                              cases[c].depth,&pole_second));
    for(n=0;n<5;n++) {
      worst_first=fmax(worst_first,relative_error(first.value[n],cases[c].first[n]));
      worst_second=fmax(worst_second,relative_error(second.value[n],cases[c].second[n]));
      worst_pole_first=fmax(worst_pole_first,
          relative_error(pole_first.value[n],cases[c].pole_first[n]));
      worst_pole_second=fmax(worst_pole_second,
          relative_error(pole_second.value[n],cases[c].pole_second[n]));
    }
  }
  printf("NEC4 outer terms: saddle %.3e/%.3e, pole %.3e/%.3e maximum relative error\n",
         worst_first,worst_second,worst_pole_first,worst_pole_second);
  assert(worst_first<1.e-8);
  assert(worst_second<1.e-8);
  assert(worst_pole_first<1.e-8);
  assert(worst_pole_second<1.e-8);
  for(c=0;c<sizeof(cases)/sizeof(cases[0]);c++) {
    static const double complex axis_expected[2][5]={
      {0.,-0.0009034589690876957+0.00041001173057707296*I,
       0.008884632748240206+0.02111477851209439*I,
       -0.008884632748240206-0.02111477851209439*I,0.},
      {0.,-0.00012691437432470548-8.762807076196731e-05*I,
       -0.0009347095025232755+0.0036315676793829774*I,
       0.0009347095025232755-0.0036315676793829774*I,0.}
    };
    nec4_ground_asymptotic_result axis,edge;
    double rho_a=.02*(z+cases[c].depth);
    assert(nec4_ground_asymptotic_second_saddle(frequency,cases[c].k1,k2,
                                                0.,z,cases[c].depth,&axis));
    assert(nec4_ground_asymptotic_second_saddle(frequency,cases[c].k1,k2,
                                                rho_a,z,cases[c].depth,&edge));
    for(n=0;n<5;n++) assert(relative_error(axis.value[n],axis_expected[c][n])<1.e-8);
    /* Ownership at rho_A is exactly the ordinary Appendix-J expression. */
    {
      nec4_ground_asymptotic_result outside;
      assert(nec4_ground_asymptotic_second_saddle(frequency,cases[c].k1,k2,
          nextafter(rho_a,INFINITY),z,cases[c].depth,&outside));
      for(n=0;n<5;n++) assert(relative_error(edge.value[n],outside.value[n])<1.e-12);
    }
  }
  puts("NEC4 Appendix-J rho=0 limits and rho_A interpolation boundary passed");
  for(c=0;c<sizeof(cases)/sizeof(cases[0]);c++) {
    nec4_ground_asymptotic_result saddle,pole,uniform;
    assert(nec4_ground_asymptotic_second_saddle(frequency,cases[c].k1,k2,
                                                rho,z,cases[c].depth,&saddle));
    assert(nec4_ground_asymptotic_second_pole(frequency,cases[c].k1,k2,
                                              rho,z,cases[c].depth,&pole));
    assert(nec4_ground_asymptotic_second_uniform(frequency,cases[c].k1,k2,
                                                 rho,z,cases[c].depth,&uniform));
    for(n=0;n<5;n++) assert(cabs(uniform.value[n]-saddle.value[n]-pole.value[n])<1.e-17);
  }
  puts("NEC4 complete second-order uniform outer field composition passed");
  return 0;
}
