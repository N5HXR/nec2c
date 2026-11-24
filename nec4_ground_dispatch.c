#include "nec4_ground_dispatch.h"

#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int source_region(double depth,double lg,double lmax)
{
  if(depth<=.4*lg) return 1;
  if(depth<=lg) return 2;
  if(depth<=lmax) return 3;
  return 0;
}

static int near_contains(double rho,double z,double depth,double lg,double lmax)
{
  return rho<=lg&&z<=lg&&source_region(depth,lg,lmax)!=0;
}

static int model_region(double rho,double z,double wavelength0,double lg)
{
  if(rho<=.6*wavelength0) return z>=fmin(rho,lg)?1:2;
  return z>=rho?1:3;
}

int nec4_ground_route(double rho,double z,double source_z,
    double complex k_lower,double complex k_upper,double wavelength0,
    nec4_ground_dispatch *route)
{
  double ratio,lg,lmax,depth,angle,separation;
  if(!route||!(wavelength0>0.)||rho<0.||z<0.||source_z>0.||
     !isfinite(rho)||!isfinite(z)||!isfinite(source_z)||
     !(cabs(k_lower)>0.)||!isfinite(cabs(k_lower))||
     !isfinite(cabs(k_upper))) return 0;
  ratio=cabs(k_upper/k_lower);
  lg=fmin(ratio,.3)*wavelength0;
  lmax=fmin(15.*ratio,1.)*wavelength0;
  depth=fabs(source_z);
  route->model_region=0;
  route->used_direct_fallback=0;
  if(near_contains(rho,z,depth,lg,lmax))
    route->method=NEC4_GROUND_LINEAR_TABLE;
  else if(rho<=3.*wavelength0&&z<=2.*wavelength0&&depth<=.8*lg)
  {
    route->method=NEC4_GROUND_MODEL_FIT;
    route->model_region=model_region(rho,z,wavelength0,lg);
  }
  else route->method=NEC4_GROUND_ASYMPTOTIC;
  angle=rho>0.?atan2(z+depth,rho)*180./M_PI:90.;
  separation=hypot(rho,z-source_z);
  route->grazing_angle_degrees=angle;
  route->separation_wavelengths=separation/wavelength0;
  route->grazing_risk=route->method==NEC4_GROUND_ASYMPTOTIC&&angle<=15.;
  return 1;
}

const char *nec4_ground_method_name(nec4_ground_method method)
{
  switch(method)
  {
    case NEC4_GROUND_LINEAR_TABLE:return "linear-table";
    case NEC4_GROUND_MODEL_FIT:return "model-fit";
    case NEC4_GROUND_ASYMPTOTIC:return "asymptotic";
    default:return "invalid";
  }
}
