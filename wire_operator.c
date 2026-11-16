#include "nec2c.h"
#include "wire_operator.h"

typedef struct
{
  const char *name;
  void (*sbf)(int, int, double *, double *, double *);
  void (*tbf)(int, int);
  void (*trio)(int);
  void (*efld)(double, double, double, double, int);
  void (*hsfld)(double, double, double, double);
} wire_operator_vtable;

static const wire_operator_vtable nec2_operator = {
  "nec2", sbf, tbf, trio, efld, hsfld
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

void wire_operator_sbf(int i, int is, double *aa, double *bb, double *cc)
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
}

void wire_operator_hsfld(double xi, double yi, double zi, double ai)
{
  active_operator->hsfld(xi, yi, zi, ai);
}
