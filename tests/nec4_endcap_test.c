#include "../nec4_endcap.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void check_close(const char *name,double complex got,
                        double complex want,double tolerance)
{
  double error=cabs(got-want)/(1.+cabs(want));
  if(error>tolerance)
  {
    fprintf(stderr,"%s error %.17g: got %.17g%+.17gj want %.17g%+.17gj\n",
      name,error,creal(got),cimag(got),creal(want),cimag(want));
    exit(1);
  }
}

int main(void)
{
  const double pi=3.141592653589793238462643383279502884;
  nec4_endcap_field f;
  double a=.01,z=.037,k=2.*pi,ra=hypot(a,z);
  double complex exact=2.*pi*cexp(-I*k*z)*
    (1.-z/ra*cexp(-I*k*(ra-z)));
  if(!nec4_endcap_disk_field(a,0.,z,k,&f) ||
     f.branch!=NEC4_ENDCAP_AXIS) return 1;
  check_close("axis",f.axial,exact,3.e-15);
  check_close("axis radial",f.radial,0.,0.);
  if(!nec4_endcap_disk_field(a,.023,.017,k,&f)) return 1;
  if(f.branch!=NEC4_ENDCAP_QUADRATURE) return 1;
  /* High-precision direct-quadrature witnesses generated independently. */
  check_close("near radial",f.radial,
    0.3013680894343398246-0.00059528019888976065*I,2.e-13);
  check_close("near axial",f.axial,
    0.2442152899578017898-0.00044007666094748538*I,2.e-13);
  if(!nec4_endcap_disk_field(a,10.,6.,k,&f) ||
     f.branch!=NEC4_ENDCAP_FAR) return 1;
  check_close("far radial",f.radial,
    -0.00012446456315795723-0.00007458677871437701*I,8.e-12);
  check_close("far axial",f.axial,
    -0.00007467816260125855-0.00004475308949834950*I,8.e-12);
  puts("NEC4 flat end-cap axis, near, and far branches passed");
  return 0;
}
