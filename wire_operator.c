#include "nec2c.h"
#include "shared.h"
#include "wire_operator.h"
#include "junction_charge.h"

typedef struct
{
  const char *name;
  void (*sbf)(int, int, double complex *, double complex *, double complex *);
  void (*tbf)(int, int);
  void (*trio)(int);
  void (*efld)(double, double, double, double, int);
  void (*hsfld)(double, double, double, double);
} wire_operator_vtable;

static const wire_operator_vtable nec2_operator = {
  "nec2", sbf, tbf, trio, efld, hsfld
};

static const wire_operator_vtable nec4_operator = {
  "nec4", sbf, tbf, trio, efld, hsfld
};

static const wire_operator_vtable *active_operator = &nec2_operator;
static const char *requested_engine = "auto";
static const char *requested_ground_engine = "auto";
static const char *resolved_ground_engine = "nec2";
static const char *model_check_policy = "warn";

int wire_operator_set_engine(const char *name)
{
  if(!name) return 0;
  if(!strcmp(name, "auto"))
  {
    requested_engine = "auto";
    active_operator = &nec2_operator;
    return 1;
  }
  if(!strcmp(name, "nec2"))
  {
    requested_engine = "nec2";
    active_operator = &nec2_operator;
    return 1;
  }
  if(!strcmp(name, "nec4"))
  {
    requested_engine = "nec4";
    active_operator = &nec4_operator;
    return 1;
  }
  return 0;
}

int wire_operator_set_ground_engine(const char *name)
{
  if(!name) return 0;
  if(!strcmp(name, "auto"))
  {
    requested_ground_engine = "auto";
    resolved_ground_engine = "nec2";
    return 1;
  }
  if(!strcmp(name, "nec2"))
  {
    requested_ground_engine = "nec2";
    resolved_ground_engine = "nec2";
    return 1;
  }
  return 0;
}

int wire_operator_set_model_check(const char *name)
{
  if(!name) return 0;
  if(strcmp(name, "warn") && strcmp(name, "strict") && strcmp(name, "off"))
    return 0;
  model_check_policy = !strcmp(name, "strict") ? "strict" :
    (!strcmp(name, "off") ? "off" : "warn");
  return 1;
}

const char *wire_operator_requested_engine(void) { return requested_engine; }
const char *wire_operator_resolved_engine(void) { return active_operator->name; }
const char *wire_operator_requested_ground_engine(void)
{
  return requested_ground_engine;
}
const char *wire_operator_resolved_ground_engine(void)
{
  return resolved_ground_engine;
}
const char *wire_operator_model_check(void) { return model_check_policy; }
int wire_operator_nec4_basis(void) { return active_operator == &nec4_operator; }

int wire_operator_prepare_frequency(FILE *report)
{
  if(!wire_operator_nec4_basis()) return 1;
  return junction_charge_prepare(report,!strcmp(model_check_policy,"strict"));
}

double complex wire_operator_charge_factor(size_t segment,int end,
                                           double complex fallback)
{
  if(!wire_operator_nec4_basis()) return fallback;
  return junction_charge_factor(segment,end,fallback);
}

double wire_operator_cosine_shape(double angle)
{
  if(!wire_operator_nec4_basis()) return cos(angle);
  return -2.0*sin(0.5*angle)*sin(0.5*angle);
}

double complex wire_operator_center_basis(double complex constant,
                                          double complex cosine)
{
  return wire_operator_nec4_basis() ? constant : constant+cosine;
}

double wire_operator_cosine_integral(double constant_integral,
                                     double cosine_integral)
{
  return wire_operator_nec4_basis() ? cosine_integral-constant_integral :
    cosine_integral;
}

double complex wire_operator_current(double complex constant,
                                     double complex sine,
                                     double complex cosine,
                                     double angle)
{
  return constant+sine*sin(angle)+cosine*wire_operator_cosine_shape(angle);
}

void wire_operator_sbf(int i, int is, double complex *aa,
                       double complex *bb, double complex *cc)
{
  active_operator->sbf(i, is, aa, bb, cc);
}

void wire_operator_tbf(int i, int icap)
{
  active_operator->tbf(i, icap);
}

void wire_operator_trio(int j)
{
  active_operator->trio(j);
}

void wire_operator_efld(double xi, double yi, double zi, double ai, int ij)
{
  active_operator->efld(xi, yi, zi, ai, ij);
  if(wire_operator_nec4_basis())
  {
    dataj.exc -= dataj.exk;
    dataj.eyc -= dataj.eyk;
    dataj.ezc -= dataj.ezk;
  }
}

void wire_operator_hsfld(double xi, double yi, double zi, double ai)
{
  active_operator->hsfld(xi, yi, zi, ai);
  if(wire_operator_nec4_basis())
  {
    dataj.exc -= dataj.exk;
    dataj.eyc -= dataj.eyk;
    dataj.ezc -= dataj.ezk;
  }
}
