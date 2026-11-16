#include "nec4_kernel.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct
{
  const char *name;
  double length,z,rho;
  nec4_field_branch expected;
} field_case;

static void pair(double complex value)
{ printf("[%.17g,%.17g]",creal(value),cimag(value)); }

static void print_field(const char *name,double length,double z,double rho,
                        double wave_number,const nec4_segment_fields *field)
{
  printf("{\"name\":\"%s\",\"length\":%.17g,\"z\":%.17g,"
         "\"rho\":%.17g,\"wave_number\":%.17g,\"branch\":\"%s\",\"ezs\":",
         name,length,z,rho,wave_number,nec4_field_branch_name(field->branch));
  pair(field->ezs);printf(",\"ers\":");pair(field->ers);
  printf(",\"ezc_minus_ezk\":");pair(field->ezc-field->ezk);
  printf(",\"erc\":");pair(field->erc);
  printf(",\"ezk\":");pair(field->ezk);
  printf(",\"erk\":");pair(field->erk);printf("}");
}

int main(int argc,char **argv)
{
  static const field_case cases[]={
    {"exact",.1,.03,.02,NEC4_FIELD_EXACT},
    {"near-axis",.1,.08,1.e-8,NEC4_FIELD_NEAR_AXIS},
    {"small-kR",2.e-5,3.e-5,2.e-5,NEC4_FIELD_SMALL_KR},
    {"small-kR-axis",2.e-5,3.e-5,1.e-11,NEC4_FIELD_SMALL_KR_AXIS},
    {"far",1.e-7,.1,.02,NEC4_FIELD_FAR},
    {"far-small-kR",1.e-10,5.e-5,2.e-5,NEC4_FIELD_FAR_SMALL_KR},
  };
  size_t i;
  if(argc==5)
  {
    double length=strtod(argv[1],NULL),z=strtod(argv[2],NULL),
      rho=strtod(argv[3],NULL),wave_number=strtod(argv[4],NULL);
    nec4_segment_fields field;
    if(!nec4_reduced_segment_field(length,z,rho,wave_number,&field)) return 1;
    print_field("sample",length,z,rho,wave_number,&field);
    putchar('\n');
    return 0;
  }
  if(argc!=1) return 64;
  puts("[");
  for(i=0;i<sizeof(cases)/sizeof(cases[0]);++i)
  {
    nec4_segment_fields field;
    if(!nec4_reduced_segment_field(cases[i].length,cases[i].z,cases[i].rho,
                                   2.*acos(-1.),&field)) return 1;
    if(field.branch!=cases[i].expected) return 2;
    if(i) puts(",");
    print_field(cases[i].name,cases[i].length,cases[i].z,cases[i].rho,
                2.*acos(-1.),&field);
  }
  puts("\n]");
  return 0;
}
