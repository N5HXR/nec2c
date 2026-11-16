#include "junction_charge.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void require(int condition,const char *message)
{
  if(!condition) { fprintf(stderr,"junction charge test: %s\n",message);exit(1); }
}

static void solve(size_t count,const double *length,const double *radius,
                  const double direction[][3],double complex *factor)
{
  double residual,pivot;
  require(junction_charge_solve(count,length,radius,direction,factor,&residual,&pivot),
          "local solve failed");
  require(residual<=1.e-12,"local residual exceeded 1e-12");
  require(isfinite(pivot),"non-finite pivot ratio");
}

int main(void)
{
  double equal_length[2]={.0125,.0125},equal_radius[2]={.001,.001};
  double step_radius[2]={.001,.005};
  double complex factor[4],permuted[2];
  double line[2][3]={{0,0,-1},{0,0,1}};
  double reverse[2][3]={{0,0,1},{0,0,-1}};
  double t_length[3]={.02,.02,.02},t_radius[3]={.001,.001,.001};
  double tee[3][3]={{1,0,0},{-1,0,0},{0,0,1}};
  double complex tf[3];
  double star_length[4]={.02,.02,.02,.02};
  double star_radius[4]={.001,.0015,.002,.003};
  double star[4][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0}};
  double complex sf[4];

  solve(2,equal_length,equal_radius,line,factor);
  require(cabs(factor[0]-1.)<1.e-12 && cabs(factor[1]-1.)<1.e-12,
          "equal-radius null limit failed");

  solve(2,equal_length,step_radius,line,factor);
  solve(2,equal_length,(double[]){step_radius[1],step_radius[0]},reverse,
        permuted);
  require(cabs(factor[0]-permuted[1])<2.e-12 &&
          cabs(factor[1]-permuted[0])<2.e-12,
          "permutation invariance failed");

  solve(3,t_length,t_radius,tee,tf);
  require(cabs(tf[0]-tf[1])<2.e-12,"T-junction mirror symmetry failed");

  solve(4,star_length,star_radius,star,sf);
  require(creal(sf[0])<creal(sf[1]) && creal(sf[1])<creal(sf[2]) &&
          creal(sf[2])<creal(sf[3]),
          "unequal-radius star charge ordering failed");

  printf("{\"permutation_defect\":%.17g,"
         "\"step\":[[%.17g,%.17g],[%.17g,%.17g]],"
         "\"tee\":[[%.17g,%.17g],[%.17g,%.17g],[%.17g,%.17g]],"
         "\"star\":[[%.17g,%.17g],[%.17g,%.17g],[%.17g,%.17g],"
         "[%.17g,%.17g]]}\n",
         fmax(cabs(factor[0]-permuted[1]),cabs(factor[1]-permuted[0])),
         creal(factor[0]),cimag(factor[0]),creal(factor[1]),cimag(factor[1]),
         creal(tf[0]),cimag(tf[0]),creal(tf[1]),cimag(tf[1]),
         creal(tf[2]),cimag(tf[2]),creal(sf[0]),cimag(sf[0]),
         creal(sf[1]),cimag(sf[1]),creal(sf[2]),cimag(sf[2]),
         creal(sf[3]),cimag(sf[3]));
  return 0;
}
