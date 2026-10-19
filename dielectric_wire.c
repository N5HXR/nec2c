#include "dielectric.h"
#include "shared.h"
#include "vie_kernel.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define DW_EPSILON_0 8.8541878128e-12
#define DW_C0 299792458.0
#define DW_PI 3.141592653589793238462643383279502884

static const dielectric_material *dw_material(const dielectric_model *model,
                                               int id)
{
  size_t index;
  for(index = 0; index < model->material_count; ++index)
    if(model->materials[index].id == id)
      return &model->materials[index];
  return NULL;
}

static void dw_error(dielectric_diagnostic *diagnostic,
                     const dielectric_model *model, const char *value,
                     const char *message)
{
  memset(diagnostic, 0, sizeof(*diagnostic));
  snprintf(diagnostic->resource, sizeof(diagnostic->resource), "%s",
           model->resource_path != NULL ? model->resource_path : "<dielectric>");
  snprintf(diagnostic->record, sizeof(diagnostic->record), "COUPLED_SOLVER");
  snprintf(diagnostic->value, sizeof(diagnostic->value), "%s", value);
  snprintf(diagnostic->message, sizeof(diagnostic->message), "%s", message);
}

static int dw_lu_solve(complex double *matrix, complex double *rhs, size_t n)
{
  size_t k, row, column;
  for(k = 0; k < n; ++k)
  {
    size_t pivot = k;
    double largest = cabs(matrix[k*n+k]);
    for(row = k+1; row < n; ++row)
    {
      double candidate = cabs(matrix[row*n+k]);
      if(candidate > largest) { largest = candidate; pivot = row; }
    }
    if(!(largest > 64.0*DBL_EPSILON))
      return 0;
    if(pivot != k)
    {
      for(column = k; column < n; ++column)
      {
        complex double temporary = matrix[k*n+column];
        matrix[k*n+column] = matrix[pivot*n+column];
        matrix[pivot*n+column] = temporary;
      }
      { complex double temporary = rhs[k]; rhs[k] = rhs[pivot]; rhs[pivot] = temporary; }
    }
    for(row = k+1; row < n; ++row)
    {
      complex double multiplier = matrix[row*n+k]/matrix[k*n+k];
      matrix[row*n+k] = multiplier;
      for(column = k+1; column < n; ++column)
        matrix[row*n+column] -= multiplier*matrix[k*n+column];
      rhs[row] -= multiplier*rhs[k];
    }
  }
  for(row = n; row-- > 0; )
  {
    complex double sum = rhs[row];
    for(column = row+1; column < n; ++column)
      sum -= matrix[row*n+column]*rhs[column];
    rhs[row] = sum/matrix[row*n+row];
  }
  return 1;
}

static int dw_gauss(int order, double *nodes, double *weights)
{
  int root;
  for(root = 0; root < (order+1)/2; ++root)
  {
    double z = cos(DW_PI*((double)root+0.75)/((double)order+0.5));
    double previous, derivative = 0.0;
    int iteration = 0;
    do
    {
      double p0 = 1.0, p1 = z;
      int degree;
      previous = z;
      for(degree = 2; degree <= order; ++degree)
      {
        double p2 = (((2.0*degree)-1.0)*z*p1-((double)degree-1.0)*p0)/degree;
        p0 = p1; p1 = p2;
      }
      derivative = order*(z*p1-p0)/(z*z-1.0);
      z = previous-p1/derivative;
      ++iteration;
    } while(fabs(z-previous) > 4.0*DBL_EPSILON && iteration < 32);
    if(iteration == 32) return 0;
    nodes[root] = -z; nodes[order-1-root] = z;
    weights[root] = 2.0/((1.0-z*z)*derivative*derivative);
    weights[order-1-root] = weights[root];
  }
  return 1;
}

/* Point-tested electric N operator integrated over a separated cuboid. */
static int dw_point_cuboid_n(const double observation[3],
                             const dielectric_cell *cell, double k0,
                             complex double result[3][3])
{
  enum { ORDER = 10 };
  double nodes[ORDER], weights[ORDER];
  int ix, iy, iz, p, q;
  memset(result, 0, 9*sizeof(complex double));
  if(!dw_gauss(ORDER, nodes, weights)) return 0;
  for(ix = 0; ix < ORDER; ++ix)
    for(iy = 0; iy < ORDER; ++iy)
      for(iz = 0; iz < ORDER; ++iz)
      {
        int indices[3] = {ix, iy, iz};
        double source[3], displacement[3], radius2 = 0.0;
        double weight = 1.0;
        complex double g, diagonal, radial;
        int axis;
        for(axis = 0; axis < 3; ++axis)
        {
          source[axis] = cell->center[axis]+0.5*cell->size[axis]*nodes[indices[axis]];
          displacement[axis] = observation[axis]-source[axis];
          radius2 += displacement[axis]*displacement[axis];
          weight *= 0.5*cell->size[axis]*weights[indices[axis]];
        }
        if(!(radius2 > 0.0)) return 0;
        {
          double radius = sqrt(radius2);
          g = cexp(-I*k0*radius)/(4.0*DW_PI*radius);
          diagonal = g*(k0*k0-I*k0/radius-1.0/radius2);
          radial = g*(-k0*k0+3.0*I*k0/radius+3.0/radius2);
          for(p = 0; p < 3; ++p)
            for(q = 0; q < 3; ++q)
              result[p][q] += weight*(p == q ? diagonal : 0.0)+
                  weight*radial*displacement[p]*displacement[q]/radius2;
        }
      }
  return 1;
}

/* Magnetic field at a point from each Cartesian component of a constant
   electric volume-current density.  With exp(+i omega t),
   H = integral grad(G) x J dv and G = exp(-i k r)/(4 pi r). */
static int dw_point_cuboid_h(const double observation[3],
                             const dielectric_cell *cell, double k0,
                             complex double result[3][3])
{
  enum { ORDER = 10 };
  double nodes[ORDER], weights[ORDER];
  int ix, iy, iz;
  memset(result, 0, 9*sizeof(complex double));
  if(!dw_gauss(ORDER, nodes, weights)) return 0;
  for(ix = 0; ix < ORDER; ++ix)
    for(iy = 0; iy < ORDER; ++iy)
      for(iz = 0; iz < ORDER; ++iz)
      {
        int indices[3] = {ix, iy, iz};
        double source[3], displacement[3], radius2 = 0.0, weight = 1.0;
        complex double gradient[3], factor;
        int axis;
        for(axis = 0; axis < 3; ++axis)
        {
          source[axis] = cell->center[axis]+0.5*cell->size[axis]*nodes[indices[axis]];
          displacement[axis] = observation[axis]-source[axis];
          radius2 += displacement[axis]*displacement[axis];
          weight *= 0.5*cell->size[axis]*weights[indices[axis]];
        }
        if(!(radius2 > 0.0)) return 0;
        {
          double radius = sqrt(radius2);
          complex double g = cexp(-I*k0*radius)/(4.0*DW_PI*radius);
          factor = weight*g*(-I*k0-1.0/radius)/radius;
        }
        for(axis = 0; axis < 3; ++axis)
          gradient[axis] = factor*displacement[axis];
        result[1][0] += gradient[2]; result[2][0] -= gradient[1];
        result[0][1] -= gradient[2]; result[2][1] += gradient[0];
        result[0][2] += gradient[1]; result[1][2] -= gradient[0];
      }
  return 1;
}

static void dw_restore_geometry(const double *x, const double *y, const double *z,
                                const double *length, const double *radius)
{
  int segment;
  for(segment = 0; segment < data.n; ++segment)
  {
    data.x[segment] = x[segment]; data.y[segment] = y[segment];
    data.z[segment] = z[segment]; data.si[segment] = length[segment];
    data.bi[segment] = radius[segment];
  }
}

static void dw_restore_patches(const double *x, const double *y,
                               const double *z, const double *area)
{
  int patch;
  for(patch = 0; patch < data.m; ++patch)
  {
    data.px[patch] = x[patch]; data.py[patch] = y[patch];
    data.pz[patch] = z[patch]; data.pbi[patch] = area[patch];
  }
}

static void dw_average_wire_field(const dielectric_cell *cell,
                                  complex double field[3])
{
  enum { ORDER = 6 };
  double node[ORDER], weight[ORDER];
  int ix, iy, iz;
  field[0] = field[1] = field[2] = 0.0;
  if(!dw_gauss(ORDER, node, weight))
    return;
  for(ix = 0; ix < ORDER; ++ix)
    for(iy = 0; iy < ORDER; ++iy)
      for(iz = 0; iz < ORDER; ++iz)
      {
        complex double ex, ey, ez;
        double normalized_weight = weight[ix]*weight[iy]*weight[iz]/8.0;
        double point[3] = {
          cell->center[0]+0.5*cell->size[0]*node[ix],
          cell->center[1]+0.5*cell->size[1]*node[iy],
          cell->center[2]+0.5*cell->size[2]*node[iz]
        };
        nefld(point[0]/data.wlam, point[1]/data.wlam, point[2]/data.wlam,
              &ex, &ey, &ez);
        field[0] += normalized_weight*ex;
        field[1] += normalized_weight*ey;
        field[2] += normalized_weight*ez;
      }
}

int dielectric_solve_conductor_coupled(const dielectric_model *model,
                                       const dielectric_mesh_result *mesh,
                                       double frequency_mhz, double rkh,
                                       int extended_kernel,
                                       const dielectric_plane_wave *wave,
                                       FILE *stream,
                                       dielectric_diagnostic *diagnostic)
{
  size_t nw = (size_t)data.n, np = (size_t)data.m;
  size_t conductor_unknowns = nw+2*np, nc = mesh->cell_count, nv = 3*nc;
  size_t total = conductor_unknowns+nv, row, column, cell;
  complex double *system = NULL, *original = NULL, *rhs = NULL, *rhs0 = NULL;
  complex double *conductor_matrix = NULL, *unit = NULL, *volume_current = NULL;
  complex double *conductor_current = NULL;
  double *x = NULL, *y = NULL, *z = NULL, *length = NULL, *radius = NULL;
  double *px = NULL, *py = NULL, *pz = NULL, *parea = NULL;
  double frequency_hz = frequency_mhz*1.0e6;
  double omega = 2.0*DW_PI*frequency_hz, k0 = omega/DW_C0;
  double relative_residual = 0.0, rhs_norm = 0.0;
  double wire_to_volume_norm = 0.0, patch_to_volume_norm = 0.0;
  double volume_to_wire_norm = 0.0, volume_to_patch_norm = 0.0;
  int restored = 0;
  vie_kernel_options kernel_options;
  if(conductor_unknowns == 0 || total == 0)
  {
    dw_error(diagnostic, model, "geometry", "M7 requires wire and/or conducting-patch geometry");
    return 0;
  }
  if(wave == NULL && (vsorc.nsant == 0 || vsorc.nvqd != 0))
  {
    dw_error(diagnostic, model, "EX", "M6 requires at least one ordinary voltage-source EX card and does not support discontinuity sources");
    return 0;
  }
  if(wave != NULL && (wave->type < 1 || wave->type > 3))
  {
    dw_error(diagnostic, model, "EX", "M7 plane-wave excitation must be type 1, 2, or 3");
    return 0;
  }
  if(total > SIZE_MAX/total || total*total > SIZE_MAX/sizeof(*system))
  {
    dw_error(diagnostic, model, "matrix", "coupled dense matrix size overflow");
    return 0;
  }
  system = calloc(total*total, sizeof(*system));
  original = calloc(total*total, sizeof(*original));
  rhs = calloc(total, sizeof(*rhs)); rhs0 = calloc(total, sizeof(*rhs0));
  conductor_matrix = calloc(conductor_unknowns*conductor_unknowns, sizeof(*conductor_matrix));
  unit = calloc(nw+3*np, sizeof(*unit));
  conductor_current = calloc(nw+3*np, sizeof(*conductor_current));
  volume_current = calloc(nv ? nv : 1, sizeof(*volume_current));
  x = calloc(nw ? nw : 1, sizeof(*x)); y = calloc(nw ? nw : 1, sizeof(*y));
  z = calloc(nw ? nw : 1, sizeof(*z)); length = calloc(nw ? nw : 1, sizeof(*length));
  radius = calloc(nw ? nw : 1, sizeof(*radius));
  px = calloc(np ? np : 1, sizeof(*px)); py = calloc(np ? np : 1, sizeof(*py));
  pz = calloc(np ? np : 1, sizeof(*pz)); parea = calloc(np ? np : 1, sizeof(*parea));
  if(!system || !original || !rhs || !rhs0 || !conductor_matrix || !unit ||
     !conductor_current || !volume_current || !x || !y || !z || !length ||
     !radius || !px || !py || !pz || !parea)
  {
    dw_error(diagnostic, model, "memory", "unable to allocate coupled dense system");
    goto failure;
  }
  for(row = 0; row < nw; ++row)
  {
    x[row]=data.x[row]; y[row]=data.y[row]; z[row]=data.z[row];
    length[row]=data.si[row]; radius[row]=data.bi[row];
    data.x[row] *= frequency_mhz/CVEL; data.y[row] *= frequency_mhz/CVEL;
    data.z[row] *= frequency_mhz/CVEL; data.si[row] *= frequency_mhz/CVEL;
    data.bi[row] *= frequency_mhz/CVEL;
  }
  for(row = 0; row < np; ++row)
  {
    double scale = frequency_mhz/CVEL;
    px[row]=data.px[row]; py[row]=data.py[row]; pz[row]=data.pz[row];
    parea[row]=data.pbi[row];
    data.px[row] *= scale; data.py[row] *= scale; data.pz[row] *= scale;
    data.pbi[row] *= scale*scale;
  }

  data.wlam = CVEL/frequency_mhz;
  gnd.ksymp = 1; gnd.iperf = 0;
  fblock((int)conductor_unknowns, (int)conductor_unknowns,
         (int)(conductor_unknowns*conductor_unknowns), data.ipsym);
  cmset((int)conductor_unknowns, conductor_matrix, rkh, extended_kernel);
  if(wave == NULL)
    etmns(0., 0., 0., 0., 0., 0., 0, rhs);
  else
    etmns(wave->theta_degrees*DW_PI/180.0,
          wave->phi_degrees*DW_PI/180.0,
          wave->eta_degrees*DW_PI/180.0, 0., 0., wave->axial_ratio,
          wave->type, rhs);
  for(row = 0; row < conductor_unknowns; ++row)
    for(column = 0; column < conductor_unknowns; ++column)
      /* NEC's factorizer transposes cmset() output before elimination. */
      system[row*total+column] = conductor_matrix[column+row*conductor_unknowns];

  vie_kernel_default_options(&kernel_options);
  kernel_options.relative_tolerance = 5.0e-7;
  kernel_options.absolute_tolerance = 5.0e-10;
  for(cell = 0; cell < nc; ++cell)
  {
    const dielectric_material *material = dw_material(model, mesh->cells[cell].material_id);
    double epsilon_r, chi;
    int component;
    if(material == NULL) { dw_error(diagnostic, model, "material", "unknown mesh material"); goto failure; }
    epsilon_r = material->epsilon_r; chi = epsilon_r-1.0;
    for(component = 0; component < 3; ++component)
    {
      complex double incident[3] = {0.0, 0.0, 0.0};
      if(wave != NULL)
      {
        double theta = wave->theta_degrees*DW_PI/180.0;
        double phi = wave->phi_degrees*DW_PI/180.0;
        double eta = wave->eta_degrees*DW_PI/180.0;
        double ct=cos(theta), st=sin(theta), cp=cos(phi), sp=sin(phi);
        double ce=cos(eta), se=sin(eta), direction[3], primary[3], secondary[3];
        complex double phase;
        double average=1.0, dot=0.0, ratio;
        primary[0]=ct*cp*ce-sp*se; primary[1]=ct*sp*ce+cp*se; primary[2]=-st*ce;
        direction[0]=-st*cp; direction[1]=-st*sp; direction[2]=-ct;
        secondary[0]=direction[1]*primary[2]-direction[2]*primary[1];
        secondary[1]=direction[2]*primary[0]-direction[0]*primary[2];
        secondary[2]=direction[0]*primary[1]-direction[1]*primary[0];
        for(component=0; component<3; ++component)
        {
          double argument=0.5*k0*direction[component]*mesh->cells[cell].size[component];
          average *= fabs(argument)<1.0e-8 ? 1.0 : sin(argument)/argument;
          dot += direction[component]*mesh->cells[cell].center[component];
        }
        phase=average*cexp(-I*k0*dot);
        ratio=wave->axial_ratio == 0.0 ? 1.0 : wave->axial_ratio;
        for(component=0; component<3; ++component)
        {
          incident[component]=primary[component];
          if(wave->type==2 || wave->type==3)
            incident[component] += I*(wave->type==2 ? -1.0 : 1.0)*ratio*secondary[component];
          incident[component] *= phase;
        }
      }
      for(component = 0; component < 3; ++component)
      {
        system[(conductor_unknowns+3*cell+component)*total+
               conductor_unknowns+3*cell+component] = epsilon_r;
        rhs[conductor_unknowns+3*cell+component] =
            I*omega*DW_EPSILON_0*chi*incident[component];
      }
    }
    for(column = cell; column < nc; ++column)
    {
      vie_cuboid testing, source;
      complex double block[3][3];
      double vt, vs;
      int p, q;
      memcpy(testing.center, mesh->cells[cell].center, sizeof(testing.center));
      memcpy(testing.size, mesh->cells[cell].size, sizeof(testing.size));
      memcpy(source.center, mesh->cells[column].center, sizeof(source.center));
      memcpy(source.size, mesh->cells[column].size, sizeof(source.size));
      if(!vie_cuboid_n(&testing, &source, k0, &kernel_options, block, NULL))
      { dw_error(diagnostic, model, "kernel", "volume-volume interaction failed"); goto failure; }
      vt = testing.size[0]*testing.size[1]*testing.size[2];
      vs = source.size[0]*source.size[1]*source.size[2];
      for(p = 0; p < 3; ++p)
        for(q = 0; q < 3; ++q)
        {
          system[(conductor_unknowns+3*cell+p)*total+
                 conductor_unknowns+3*column+q] -= chi*block[p][q];
          if(column != cell)
          {
            const dielectric_material *other = dw_material(model, mesh->cells[column].material_id);
            system[(conductor_unknowns+3*column+q)*total+
                   conductor_unknowns+3*cell+p] -=
                (other->epsilon_r-1.0)*(vt/vs)*block[p][q];
          }
        }
    }
  }

  /* Volume rows: averaged NEC wire field produced by each wire unknown. */
  if(nc != 0)
    for(column = 0; column < conductor_unknowns; ++column)
    {
      memset(unit, 0, (nw+3*np)*sizeof(*unit)); unit[column] = 1.0;
      cabc(unit);
      if(np != 0)
        memcpy(crnt.cur+nw, unit+nw, 3*np*sizeof(*unit));
      for(cell = 0; cell < nc; ++cell)
      {
        const dielectric_material *material = dw_material(model, mesh->cells[cell].material_id);
        double chi = material->epsilon_r-1.0;
        complex double field[3];
        int component;
        dw_average_wire_field(&mesh->cells[cell], field);
        for(component = 0; component < 3; ++component)
          system[(conductor_unknowns+3*cell+component)*total+column] =
              -I*omega*DW_EPSILON_0*chi*field[component];
      }
    }

  /* Wire rows: point-tested field of each constant cuboid current. */
  for(row = 0; row < nw; ++row)
  {
    double observation[3] = {x[row], y[row], z[row]};
    double tangent[3] = {data.cab[row], data.sab[row], data.salp[row]};
    for(cell = 0; cell < nc; ++cell)
    {
      complex double block[3][3];
      int p, q;
      if(!dw_point_cuboid_n(observation, &mesh->cells[cell], k0, block))
      { dw_error(diagnostic, model, "cross-kernel", "volume-to-wire interaction failed"); goto failure; }
      for(q = 0; q < 3; ++q)
      {
        complex double tested = 0.0;
        for(p = 0; p < 3; ++p) tested += tangent[p]*block[p][q];
        system[row*total+conductor_unknowns+3*cell+q] =
            tested/(I*omega*DW_EPSILON_0);
      }
    }
  }

  /* Patch MFIE rows: point-tested H = integral grad(G) x Jp dv.  NEC's
     first and second patch equations test -t2.H and -t1.H respectively. */
  for(row = 0; row < np; ++row)
  {
    double observation[3] = {px[row], py[row], pz[row]};
    double t1[3] = {data.t1x[row], data.t1y[row], data.t1z[row]};
    double t2[3] = {data.t2x[row], data.t2y[row], data.t2z[row]};
    for(cell = 0; cell < nc; ++cell)
    {
      complex double block[3][3];
      int p, q;
      if(!dw_point_cuboid_h(observation, &mesh->cells[cell], k0, block))
      { dw_error(diagnostic, model, "cross-kernel", "volume-to-patch interaction failed"); goto failure; }
      for(q = 0; q < 3; ++q)
      {
        complex double tested_t1=0.0, tested_t2=0.0;
        for(p = 0; p < 3; ++p)
        {
          tested_t1 += t1[p]*block[p][q];
          tested_t2 += t2[p]*block[p][q];
        }
        system[(nw+2*row)*total+conductor_unknowns+3*cell+q] =
            -data.psalp[row]*tested_t2;
        system[(nw+2*row+1)*total+conductor_unknowns+3*cell+q] =
            -data.psalp[row]*tested_t1;
      }
    }
  }

  for(row = conductor_unknowns; row < total; ++row)
  {
    for(column = 0; column < nw; ++column)
      wire_to_volume_norm += creal(system[row*total+column]*
                                   conj(system[row*total+column]));
    for(column = nw; column < conductor_unknowns; ++column)
      patch_to_volume_norm += creal(system[row*total+column]*
                                    conj(system[row*total+column]));
  }
  for(row = 0; row < nw; ++row)
    for(column = conductor_unknowns; column < total; ++column)
      volume_to_wire_norm += creal(system[row*total+column]*
                                   conj(system[row*total+column]));
  for(row = nw; row < conductor_unknowns; ++row)
    for(column = conductor_unknowns; column < total; ++column)
      volume_to_patch_norm += creal(system[row*total+column]*
                                    conj(system[row*total+column]));
  wire_to_volume_norm=sqrt(wire_to_volume_norm);
  patch_to_volume_norm=sqrt(patch_to_volume_norm);
  volume_to_wire_norm=sqrt(volume_to_wire_norm);
  volume_to_patch_norm=sqrt(volume_to_patch_norm);

  memcpy(original, system, total*total*sizeof(*system));
  memcpy(rhs0, rhs, total*sizeof(*rhs));
  if(!dw_lu_solve(system, rhs, total))
  { dw_error(diagnostic, model, "matrix", "singular coupled wire-volume system"); goto failure; }
  for(row = 0; row < total; ++row)
  {
    complex double value = -rhs0[row];
    rhs_norm += creal(rhs0[row]*conj(rhs0[row]));
    for(column = 0; column < total; ++column)
      value += original[row*total+column]*rhs[column];
    relative_residual += creal(value*conj(value));
  }
  relative_residual = sqrt(relative_residual/fmax(rhs_norm, DBL_MIN));
  memcpy(volume_current, rhs+conductor_unknowns, nv*sizeof(*volume_current));
  memcpy(conductor_current, rhs, conductor_unknowns*sizeof(*conductor_current));
  cabc(conductor_current);

  fprintf(stream, "\n\n  ---------- COUPLED CONDUCTOR-DIELECTRIC SOLUTION (M7) ----------"
                  "\n  FREQUENCY (MHZ): %12.6f"
                  "\n  WIRE UNKNOWNS:   %12zu"
                  "\n  PATCH UNKNOWNS:  %12zu"
                  "\n  VOLUME UNKNOWNS: %12zu"
                  "\n  TOTAL UNKNOWNS:  %12zu"
                  "\n  RELATIVE RESIDUAL:%11.4E"
                  "\n  WIRE-TO-VOLUME BLOCK NORM:  %11.4E"
                  "\n  VOLUME-TO-WIRE BLOCK NORM:  %11.4E"
                  "\n  PATCH-TO-VOLUME BLOCK NORM: %11.4E"
                  "\n  VOLUME-TO-PATCH BLOCK NORM: %11.4E",
          frequency_mhz, nw, 2*np, nv, total, relative_residual,
          wire_to_volume_norm, volume_to_wire_norm,
          patch_to_volume_norm, volume_to_patch_norm);
  fprintf(stream, "\n\n  COUPLED WIRE CURRENTS"
                  "\n  SEG TAG       X(M)       Y(M)       Z(M)          CURRENT(RE,IM) A");
  for(row = 0; row < nw; ++row)
  {
    complex double current = conductor_current[row]*data.wlam;
    fprintf(stream, "\n %4zu %3d %10.4E %10.4E %10.4E  %19.12E %19.12E",
            row+1, data.itag[row], x[row], y[row], z[row],
            creal(current), cimag(current));
  }
  fprintf(stream, "\n\n  COUPLED ANTENNA INPUT PARAMETERS"
                  "\n  TAG SEG       VOLTAGE(RE,IM) V          CURRENT(RE,IM) A          IMPEDANCE(RE,IM) OHM");
  for(row = 0; row < (size_t)vsorc.nsant; ++row)
  {
    size_t source = (size_t)(vsorc.isant[row]-1);
    complex double voltage = vsorc.vsant[row];
    complex double current = conductor_current[source]*data.wlam;
    complex double impedance = voltage/current;
    fprintf(stream, "\n %4d %3zu  %19.12E %19.12E  %19.12E %19.12E  %19.12E %19.12E",
            data.itag[source], source+1, creal(voltage), cimag(voltage),
            creal(current), cimag(current), creal(impedance), cimag(impedance));
  }
  if(np != 0)
  {
    fprintf(stream, "\n\n  COUPLED SURFACE PATCH CURRENTS (A/M)"
                    "\n  PATCH          JX(RE,IM)                 JY(RE,IM)                 JZ(RE,IM)");
    for(row = 0; row < np; ++row)
    {
      size_t offset = nw+3*row;
      fprintf(stream, "\n %5zu  %19.12E %19.12E  %19.12E %19.12E  %19.12E %19.12E",
              row+1, creal(conductor_current[offset]), cimag(conductor_current[offset]),
              creal(conductor_current[offset+1]), cimag(conductor_current[offset+1]),
              creal(conductor_current[offset+2]), cimag(conductor_current[offset+2]));
    }
  }
  fprintf(stream, "\n\n  COUPLED VOLUME POLARIZATION CURRENTS (A/M^2)"
                  "\n  CELL REGION MATERIAL          JX(RE,IM)                 JY(RE,IM)                 JZ(RE,IM)");
  for(cell = 0; cell < nc; ++cell)
    fprintf(stream, "\n %5zu %6d %8d  %11.4E %11.4E  %11.4E %11.4E  %11.4E %11.4E",
            cell+1, mesh->cells[cell].region_id, mesh->cells[cell].material_id,
            creal(volume_current[3*cell]), cimag(volume_current[3*cell]),
            creal(volume_current[3*cell+1]), cimag(volume_current[3*cell+1]),
            creal(volume_current[3*cell+2]), cimag(volume_current[3*cell+2]));

  dw_restore_geometry(x, y, z, length, radius);
  dw_restore_patches(px, py, pz, parea); restored = 1;
  free(parea); free(pz); free(py); free(px);
  free(radius); free(length); free(z); free(y); free(x); free(volume_current);
  free(conductor_current); free(unit); free(conductor_matrix);
  free(rhs0); free(rhs); free(original); free(system);
  return 1;

failure:
  if(!restored && x && y && z && length && radius)
    dw_restore_geometry(x, y, z, length, radius);
  if(!restored && px && py && pz && parea)
    dw_restore_patches(px, py, pz, parea);
  free(parea); free(pz); free(py); free(px);
  free(radius); free(length); free(z); free(y); free(x); free(volume_current);
  free(conductor_current); free(unit); free(conductor_matrix);
  free(rhs0); free(rhs); free(original); free(system);
  return 0;
}
