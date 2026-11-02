#include "dielectric.h"
#include "shared.h"
#include "vie_kernel.h"
#include "vie_ground.h"
#include "vie_halfspace.h"
#include "vie_iterative.h"
#include "dielectric_fft.h"
#include "result_json.h"
#include "dense_backend.h"
#include "performance.h"
#include "concurrency.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef _WIN32
#include <errno.h>
#include <signal.h>
#include <sys/mman.h>
#ifdef __linux__
#include <sys/prctl.h>
#endif
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

#define DW_EPSILON_0 8.8541878128e-12
#define DW_C0 299792458.0
#define DW_PI 3.141592653589793238462643383279502884

static size_t dw_coupling_processes = 1;
static int dw_coupling_fork_safe;
static complex double *dw_entry(complex double *dense,
    dielectric_fft_system *fft,size_t order,size_t index);

void dielectric_coupling_configure(size_t processes, int fork_safe)
{
  dw_coupling_processes = processes ? processes : 1;
  dw_coupling_fork_safe = fork_safe;
}

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
static int dw_point_cuboid_n_free(const double observation[3],
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
static int dw_point_cuboid_h_free(const double observation[3],
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

static int dw_point_background(const double observation[3],
                       const dielectric_cell *cell, double k0,
                       complex double result[3][3], int magnetic)
{
  int (*evaluate)(const double *,const dielectric_cell *,double,complex double [3][3]) =
      magnetic ? dw_point_cuboid_h_free:dw_point_cuboid_n_free;
  if(!evaluate(observation,cell,k0,result)) return 0;
  if(gnd.ksymp != 1)
  {
    dielectric_cell image=*cell;
    complex double reflected[3][3];
    int p,q;
    if(gnd.iperf==1)
    {
      image.center[2]=-image.center[2];
      if(!evaluate(observation,&image,k0,reflected)) return 0;
      for(p=0;p<3;++p) for(q=0;q<3;++q)
        result[p][q]+=(q==2 ? 1.0:-1.0)*reflected[p][q];
    }
    else
    {
      vie_cuboid source;
      complex double electric[3][3],h[3][3];
      memcpy(source.center,cell->center,sizeof(source.center));
      memcpy(source.size,cell->size,sizeof(source.size));
      if(!vie_halfspace_point_cuboid(observation,&source,k0,
           dielectric_ground_epsilon(k0*DW_C0/(2*DW_PI*1e6)),NULL,electric,h,NULL)) return 0;
      for(p=0;p<3;++p) for(q=0;q<3;++q) result[p][q]+=magnetic ? h[p][q]:electric[p][q];
    }
  }
  return 1;
}

static int dw_point_cuboid_n(const double observation[3],
                       const dielectric_cell *cell, double k0,
                       complex double result[3][3])
{ return dw_point_background(observation,cell,k0,result,0); }

static int dw_point_cuboid_h(const double observation[3],
                       const dielectric_cell *cell, double k0,
                       complex double result[3][3])
{ return dw_point_background(observation,cell,k0,result,1); }

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

/* Integrate NEC's actual sinusoidal wire current basis against the reflected
   half-space dyad. Native free-space fields retain NEC's wire-radius model. */
static int dw_wire_reflected(const dielectric_cell *cell,const double *point,
                             double k0,int magnetic,complex double field[3])
{
  enum { ORDER=8 };
  double nodes[ORDER],weights[ORDER];
  int segment,j,p,q;
  vie_cuboid testing;
  if(!dw_gauss(ORDER,nodes,weights)) return 0;
  if(cell) { memcpy(testing.center,cell->center,sizeof(testing.center)); memcpy(testing.size,cell->size,sizeof(testing.size)); }
  for(segment=0;segment<data.n;++segment)
  {
    complex double a=crnt.air[segment]+I*crnt.aii[segment];
    complex double b=crnt.bir[segment]+I*crnt.bii[segment];
    complex double c=crnt.cir[segment]+I*crnt.cii[segment];
    double tangent[3]={data.cab[segment],data.sab[segment],data.salp[segment]};
    if(a==0 && b==0 && c==0) continue;
    for(j=0;j<ORDER;++j)
    {
      double ell=.5*data.si[segment]*nodes[j],source[3];
      complex double n[3][3],h[3][3];
      complex double moment=data.wlam*data.wlam*.5*data.si[segment]*weights[j]*
                             (a+b*sin(TP*ell)+c*cos(TP*ell));
      for(p=0;p<3;++p) source[p]=data.wlam*((p==0 ? data.x[segment]:p==1 ? data.y[segment]:data.z[segment])+ell*tangent[p]);
      if(cell ? !vie_halfspace_cuboid_point(&testing,source,k0,dielectric_ground_epsilon(k0*DW_C0/(2*DW_PI*1e6)),NULL,n,h,NULL):
                !vie_halfspace_point(point,source,k0,dielectric_ground_epsilon(k0*DW_C0/(2*DW_PI*1e6)),NULL,n,h,NULL)) return 0;
      for(p=0;p<3;++p) for(q=0;q<3;++q)
        field[p]+=moment*tangent[q]*(magnetic ? h[p][q]:n[p][q]/(I*(k0*DW_C0)*DW_EPSILON_0));
    }
  }
  return 1;
}

/* Cache the three reflected NEC basis fields per segment and testing cell.
   This avoids repeating the same Sommerfeld integrations for every unknown. */
static int dw_wire_background_basis(const dielectric_cell *cell,double k0,
                                    complex double *basis)
{
  enum { ORDER=8 };
  double nodes[ORDER],weights[ORDER];
  vie_cuboid testing;
  int segment,j,p,q,kind;
  memcpy(testing.center,cell->center,sizeof(testing.center));
  memcpy(testing.size,cell->size,sizeof(testing.size));
  if(!dw_gauss(ORDER,nodes,weights)) return 0;
  for(segment=0;segment<data.n;++segment)
  {
    double tangent[3]={data.cab[segment],data.sab[segment],data.salp[segment]};
    for(j=0;j<ORDER;++j)
    {
      double ell=.5*data.si[segment]*nodes[j],source[3];
      double shape[3]={1,sin(TP*ell),cos(TP*ell)};
      complex double n[3][3],h[3][3];
      complex double weight=data.wlam*data.wlam*.5*data.si[segment]*weights[j]/
                              (I*(k0*DW_C0)*DW_EPSILON_0);
      for(p=0;p<3;++p) source[p]=data.wlam*((p==0 ? data.x[segment]:p==1 ? data.y[segment]:data.z[segment])+ell*tangent[p]);
      if(!vie_halfspace_cuboid_point(&testing,source,k0,
           dielectric_ground_epsilon(k0*DW_C0/(2*DW_PI*1e6)),NULL,n,h,NULL)) return 0;
      for(p=0;p<3;++p) for(q=0;q<3;++q) for(kind=0;kind<3;++kind)
        basis[9*segment+3*kind+p]+=weight*n[p][q]*tangent[q]*shape[kind];
    }
  }
  return 1;
}

static int dw_average_wire_field(const dielectric_cell *cell, double k0,
                                  const complex double *basis,complex double field[3])
{
  enum { ORDER = 6 };
  double node[ORDER], weight[ORDER];
  int ix, iy, iz;
  field[0] = field[1] = field[2] = 0.0;
  int background=gnd.ksymp;
  if(!dw_gauss(ORDER, node, weight)) return 0;
  if(background!=1 && gnd.iperf==2) gnd.ksymp=1;
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
  gnd.ksymp=background;
  if(background!=1 && gnd.iperf==2)
  {
    int segment,p;
    if(!basis) return dw_wire_reflected(cell,NULL,k0,0,field);
    for(segment=0;segment<data.n;++segment)
      for(p=0;p<3;++p)
        field[p]+=basis[9*segment+p]*(crnt.air[segment]+I*crnt.aii[segment])+
                  basis[9*segment+3+p]*(crnt.bir[segment]+I*crnt.bii[segment])+
                  basis[9*segment+6+p]*(crnt.cir[segment]+I*crnt.cii[segment]);
  }
  return 1;
}

/* Freeze each conductor unknown into explicit NEC basis coefficients before
   parallel work begins.  Workers install only their process-private copy. */
static int dw_prepare_conductor_sources(size_t nw,size_t np,size_t count,
                                        complex double *unit,
                                        complex double *sources)
{
  size_t column,segment,stride=3*nw+3*np;
  for(column=0;column<count;++column)
  {
    complex double *source=sources+column*stride;
    memset(unit,0,(nw+3*np)*sizeof(*unit));unit[column]=1.0;
    cabc(unit);
    for(segment=0;segment<nw;++segment)
    {
      source[3*segment]=crnt.air[segment]+I*crnt.aii[segment];
      source[3*segment+1]=crnt.bir[segment]+I*crnt.bii[segment];
      source[3*segment+2]=crnt.cir[segment]+I*crnt.cii[segment];
    }
    if(np) memcpy(source+3*nw,unit+nw,3*np*sizeof(*source));
  }
  return 1;
}

static void dw_install_conductor_source(size_t nw,size_t np,
                                         const complex double *source)
{
  size_t segment;
  for(segment=0;segment<nw;++segment)
  {
    crnt.air[segment]=creal(source[3*segment]);
    crnt.aii[segment]=cimag(source[3*segment]);
    crnt.bir[segment]=creal(source[3*segment+1]);
    crnt.bii[segment]=cimag(source[3*segment+1]);
    crnt.cir[segment]=creal(source[3*segment+2]);
    crnt.cii[segment]=cimag(source[3*segment+2]);
  }
  if(np) memcpy(crnt.cur+nw,source+3*nw,3*np*sizeof(*source));
}

static int dw_wire_to_volume_range(const dielectric_mesh_result *mesh,
    size_t nw,size_t np,size_t conductor_unknowns,const complex double *sources,
    const complex double *wire_background_basis,const complex double *volume_chi,
    double k0,double omega,size_t first,size_t last,complex double *coupling)
{
  size_t task,current_column=SIZE_MAX,stride=3*nw+3*np,nv=3*mesh->cell_count;
  for(task=first;task<last;++task)
  {
    size_t column=task/mesh->cell_count,cell=task%mesh->cell_count;
    complex double field[3];int component;
    if(column!=current_column)
    {
      dw_install_conductor_source(nw,np,sources+column*stride);
      current_column=column;
    }
    if(!dw_average_wire_field(&mesh->cells[cell],k0,
         wire_background_basis ? wire_background_basis+9*nw*cell:NULL,field))
      return 0;
    for(component=0;component<3;++component)
      coupling[column*nv+3*cell+component]=
          -I*omega*DW_EPSILON_0*volume_chi[cell]*field[component];
  }
  (void)conductor_unknowns;
  return 1;
}

static int dw_assemble_wire_to_volume(const dielectric_mesh_result *mesh,
    size_t nw,size_t np,size_t conductor_unknowns,const complex double *sources,
    const complex double *wire_background_basis,const complex double *volume_chi,
    double k0,double omega,complex double **result,int *result_mapped,
    size_t *workers_used,int *fell_back)
{
  size_t nv=3*mesh->cell_count,tasks,values,bytes,index;
  complex double *coupling=NULL;
  int mapped=0,parallel=0,attempted=0;
  *workers_used=1;
  *fell_back=0;
  if(mesh->cell_count && conductor_unknowns>SIZE_MAX/mesh->cell_count) return 0;
  tasks=conductor_unknowns*mesh->cell_count;
  if(conductor_unknowns && nv>SIZE_MAX/conductor_unknowns) return 0;
  values=conductor_unknowns*nv;
  if(values>SIZE_MAX/sizeof(*coupling)) return 0;
  bytes=values*sizeof(*coupling);
#ifndef _WIN32
  if(dw_coupling_fork_safe && dw_coupling_processes>1 && tasks>=128 && bytes)
  {
    attempted=1;
    coupling=mmap(NULL,bytes,PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANONYMOUS,-1,0);
    if(coupling!=MAP_FAILED) {mapped=1;parallel=1;}
    else coupling=NULL;
  }
#endif
  if(!coupling) coupling=calloc(values ? values:1,sizeof(*coupling));
  if(!coupling) return 0;
#ifndef _WIN32
  if(parallel)
  {
    size_t workers=dw_coupling_processes,started=0;
    pid_t *pids;
    int failed=0;
    if(workers>tasks) workers=tasks;
    pids=calloc(workers,sizeof(*pids));
    if(!pids) failed=1;
    for(index=0;!failed && index<workers;++index)
    {
      size_t quotient=tasks/workers,remainder=tasks%workers;
      size_t first=quotient*index+(index<remainder ? index:remainder);
      size_t last=first+quotient+(index<remainder ? 1:0);
      pid_t child=fork();
      if(child==0)
      {
        worker_process_enable();
#ifdef __linux__
        if(prctl(PR_SET_PDEATHSIG,SIGTERM)<0 || getppid()==1) _exit(2);
#endif
        _exit(dw_wire_to_volume_range(mesh,nw,np,conductor_unknowns,sources,
            wire_background_basis,volume_chi,k0,omega,first,last,coupling) ? 0:1);
      }
      if(child<0) {failed=1;break;}
      pids[index]=child;++started;
    }
    for(index=0;index<started;++index)
    {
      int status=0;pid_t waited;
      do waited=waitpid(pids[index],&status,0); while(waited<0 && errno==EINTR);
      if(waited<0 || !WIFEXITED(status) || WEXITSTATUS(status)!=0) failed=1;
    }
    free(pids);
    if(failed)
    {
      memset(coupling,0,bytes);
      parallel=0;
      *fell_back=1;
    }
    else *workers_used=workers;
  }
#endif
  if(attempted && !parallel) *fell_back=1;
  if(!parallel && !dw_wire_to_volume_range(mesh,nw,np,conductor_unknowns,sources,
       wire_background_basis,volume_chi,k0,omega,0,tasks,coupling))
  {
#ifndef _WIN32
    if(mapped)
      munmap(coupling,bytes);
    else
#endif
      free(coupling);
    return 0;
  }
  *result=coupling;
  *result_mapped=mapped;
  return 1;
}

static void dw_free_coupling(complex double *coupling,size_t values,int mapped)
{
  if(!coupling) return;
#ifndef _WIN32
  if(mapped)
    munmap(coupling,values*sizeof(*coupling));
  else
#else
  (void)values;(void)mapped;
#endif
    free(coupling);
}

static int dw_assemble_volume_pairs(const dielectric_mesh_result *mesh,
    double k0,const vie_kernel_options *options,const complex double *volume_chi,
    const double *volume_size,complex double *system,size_t total,
    size_t conductor_unknowns,size_t *failed_testing,size_t *failed_source)
{
  size_t nc=mesh->cell_count,cell;
  size_t *failed_columns=malloc(nc*sizeof(*failed_columns));
  ptrdiff_t cell_index;
  if(!failed_columns) return 0;
  for(cell=0;cell<nc;++cell) failed_columns[cell]=SIZE_MAX;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) if(nc >= 8)
#endif
  for(cell_index=0;cell_index<(ptrdiff_t)nc;++cell_index)
  {
    size_t testing_cell=(size_t)cell_index,source_cell;
    for(source_cell=testing_cell;source_cell<nc;++source_cell)
    {
      vie_cuboid testing,source;
      complex double block[3][3];int p,q;
      memcpy(testing.center,mesh->cells[testing_cell].center,sizeof(testing.center));
      memcpy(testing.size,mesh->cells[testing_cell].size,sizeof(testing.size));
      memcpy(source.center,mesh->cells[source_cell].center,sizeof(source.center));
      memcpy(source.size,mesh->cells[source_cell].size,sizeof(source.size));
      if(!dielectric_background_n(&testing,&source,k0,options,block))
      {failed_columns[testing_cell]=source_cell;break;}
      for(p=0;p<3;++p) for(q=0;q<3;++q)
      {
        *dw_entry(system,NULL,total,(conductor_unknowns+3*testing_cell+p)*total+
            conductor_unknowns+3*source_cell+q)-=volume_chi[testing_cell]*block[p][q];
        if(source_cell!=testing_cell)
          *dw_entry(system,NULL,total,(conductor_unknowns+3*source_cell+q)*total+
              conductor_unknowns+3*testing_cell+p)-=
              volume_chi[source_cell]*(volume_size[testing_cell]/volume_size[source_cell])*block[p][q];
      }
    }
  }
  for(cell=0;cell<nc;++cell) if(failed_columns[cell]!=SIZE_MAX)
  {
    *failed_testing=cell;*failed_source=failed_columns[cell];free(failed_columns);
    return 0;
  }
  free(failed_columns);return 1;
}

/* Exterior scattered fields only. The separated-cell rule also bounds
   observation quadrature error; internal/interface sampling needs a different
   singular integration and must never silently use this evaluator. */
int dielectric_print_near_fields(const dielectric_mesh_result *mesh,
                                  const complex double *current, double k0,
                                  int conductor, FILE *stream)
{
  int ix, iy, iz;
  double omega = k0*DW_C0;
  fprintf(stream, "\n\n  DIELECTRIC SCATTERED NEAR %s FIELD (%s)"
                  "\n  X(M) Y(M) Z(M)    FX(RE,IM) FY(RE,IM) FZ(RE,IM)",
          fpat.nfeh ? "MAGNETIC" : "ELECTRIC", fpat.nfeh ? "A/M" : "V/M");
  if(fpat.nrx <= 0 || fpat.nry <= 0 || fpat.nrz <= 0) return 0;
  result_json_near_begin((size_t)result_json_current_request(),
      fpat.nfeh ? "H" : "E", result_json_current_excitation_kind(),
      fpat.nfeh ? "A/m" : "V/m");
  for(iz = 0; iz < fpat.nrz; ++iz)
    for(iy = 0; iy < fpat.nry; ++iy)
      for(ix = 0; ix < fpat.nrx; ++ix)
      {
        double point[3] = {fpat.xnr+ix*fpat.dxnr, fpat.ynr+iy*fpat.dynr,
                           fpat.znr+iz*fpat.dznr};
        complex double field[3] = {0.0, 0.0, 0.0};
        size_t cell;
        int p, q;
        if(fpat.near != 0)
        {
          double r = point[0], phi = point[1]*DW_PI/180.0;
          double theta = point[2]*DW_PI/180.0;
          point[0] = r*sin(theta)*cos(phi);
          point[1] = r*sin(theta)*sin(phi);
          point[2] = r*cos(theta);
        }
        if(gnd.ksymp != 1 && (!isfinite(point[2]) || point[2] < 0.0)) return 0;
        for(cell = 0; cell < mesh->cell_count; ++cell)
        {
          double gap2 = 0.0, longest = 0.0;
          complex double block[3][3];
          for(p = 0; p < 3; ++p)
          {
            double gap = fmax(0.0, fabs(point[p]-mesh->cells[cell].center[p])-
                                   0.5*mesh->cells[cell].size[p]);
            gap2 += gap*gap;
            longest = fmax(longest, mesh->cells[cell].size[p]);
          }
          if(!isfinite(gap2) || gap2 < longest*longest) return 0;
          if(fpat.nfeh)
          {
            if(!dw_point_cuboid_h(point, &mesh->cells[cell], k0, block)) return 0;
          }
          else
          {
            if(!dw_point_cuboid_n(point, &mesh->cells[cell], k0, block)) return 0;
            for(p = 0; p < 3; ++p)
              for(q = 0; q < 3; ++q)
                block[p][q] /= I*omega*DW_EPSILON_0;
          }
          for(p = 0; p < 3; ++p)
            for(q = 0; q < 3; ++q)
              field[p] += block[p][q]*current[3*cell+q];
        }
        if(conductor)
        {
          complex double f[3];
          if(fpat.nfeh)
          {
            int background=gnd.ksymp;
            if(background!=1 && gnd.iperf==2) gnd.ksymp=1;
            nhfld(point[0]/data.wlam, point[1]/data.wlam, point[2]/data.wlam,
                   &f[0], &f[1], &f[2]);
            gnd.ksymp=background;
            if(background!=1 && gnd.iperf==2 && !dw_wire_reflected(NULL,point,k0,1,f)) return 0;
          }
          else
          {
            int background=gnd.ksymp;
            if(background!=1 && gnd.iperf==2) gnd.ksymp=1;
            nefld(point[0]/data.wlam, point[1]/data.wlam, point[2]/data.wlam,
                   &f[0], &f[1], &f[2]);
            gnd.ksymp=background;
            if(background!=1 && gnd.iperf==2 && !dw_wire_reflected(NULL,point,k0,0,f)) return 0;
          }
          for(p = 0; p < 3; ++p) field[p] += f[p];
        }
        fprintf(stream, "\n %19.12E %19.12E %19.12E", point[0], point[1], point[2]);
        for(p = 0; p < 3; ++p)
          fprintf(stream, " %19.12E %19.12E", creal(field[p]), cimag(field[p]));
        result_json_near_sample(point,field,1);
      }
  return 1;
}

/* ffld uses wavelength-normalized distance; convert its result to SI r E.
   Add complex fields before forming intensity to retain interference. */
static void dw_far_field(const dielectric_mesh_result *mesh,
                          const complex double *current, double k0,
                          double theta, double phi,
                          complex double *et, complex double *ep)
{
  double direction[3] = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
  double t[3] = {cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)};
  double p[3] = {-sin(phi), cos(phi), 0.0};
  complex double a[3];
  int axis;
  {
    int background=gnd.ksymp;
    if(background!=1 && gnd.iperf==2 && gnd.zrati==1.0) gnd.ksymp=1;
    ffld(theta, phi, et, ep);
    gnd.ksymp=background;
  }
  *et *= data.wlam; *ep *= data.wlam;
  dielectric_radiation_amplitude(mesh, current, k0, direction, a);
  for(axis = 0; axis < 3; ++axis)
  {
    *et += a[axis]*t[axis];
    *ep += a[axis]*p[axis];
  }
}

static double dw_radiated_power(const dielectric_mesh_result *mesh,
                                 const complex double *current, double k0,
                                 size_t count)
{
  size_t i;
  double sum = 0.0;
  for(i = 0; i < count; ++i)
  {
    double z = 1.0-(gnd.ksymp == 1 ? 2.0:1.0)*((double)i+0.5)/(double)count;
    double phi = DW_PI*(3.0-sqrt(5.0))*(double)i;
    complex double et, ep;
    dw_far_field(mesh, current, k0, acos(z), phi, &et, &ep);
    sum += creal(et*conj(et)+ep*conj(ep));
  }
  return sum*((gnd.ksymp == 1 ? 4.0:2.0)*DW_PI/(double)count)/(2.0*ETA);
}

static complex double *dw_entry(complex double *dense,dielectric_fft_system *fft,
                                  size_t order,size_t index)
{
  return fft ? dielectric_fft_entry(fft,index/order,index%order):dense+index;
}

int dielectric_solve_conductor_coupled(const dielectric_model *model,
                                       const dielectric_mesh_result *mesh,
                                       double frequency_mhz, double rkh,
                                       int extended_kernel, int iterative,
                                       const dielectric_plane_wave *wave,
                                       const dielectric_pattern_request *pattern,
                                       FILE *stream,
                                       dielectric_diagnostic *diagnostic)
{
  size_t nw = (size_t)data.n, np = (size_t)data.m;
  size_t conductor_unknowns = nw+2*np, nc = mesh->cell_count, nv = 3*nc;
  size_t total = conductor_unknowns+nv, row, column, cell, source_stride;
  complex double *system = NULL, *original = NULL, *rhs = NULL, *rhs0 = NULL;
  complex double *conductor_matrix = NULL, *unit = NULL, *volume_current = NULL;
  complex double *conductor_current = NULL, *wire_background_basis = NULL;
  complex double *volume_chi = NULL, *conductor_sources = NULL;
  complex double *wire_to_volume = NULL;
  double *volume_size = NULL;
  double *x = NULL, *y = NULL, *z = NULL, *length = NULL, *radius = NULL;
  double *px = NULL, *py = NULL, *pz = NULL, *parea = NULL;
  double frequency_hz = frequency_mhz*1.0e6;
  double omega = 2.0*DW_PI*frequency_hz, k0 = omega/DW_C0;
  double relative_residual = 0.0, rhs_norm = 0.0;
  clock_t started = clock();
  vie_gmres_report report;
  double wire_to_volume_norm = 0.0, patch_to_volume_norm = 0.0;
  double volume_to_wire_norm = 0.0, volume_to_patch_norm = 0.0;
  int restored = 0, geometry_scaled = 0, wire_to_volume_mapped = 0;
  int coupling_fell_back = 0;
  size_t coupling_workers_used = 1;
  dielectric_fft_system *accelerated = NULL;
  vie_kernel_options kernel_options;
  memset(&report,0,sizeof(report));
  if(!dielectric_ground_validate(model,mesh,iterative,wave,pattern,diagnostic))
    return 0;
  if(gnd.ksymp != 1)
    fprintf(stream,"\n\n  DIELECTRIC BACKGROUND: %s",gnd.iperf==1 ?
      "PERFECT GROUND (GN 1; DIELECTRIC Z >= 0)":"SOMMERFELD GROUND (GN 2; WIRES; POSITIVE GAP)");
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
  if(iterative == 2 &&
     (total>DIELECTRIC_FFT_BUDGET/(8*sizeof(complex double)) ||
      conductor_unknowns>DIELECTRIC_FFT_BUDGET/(sizeof(complex double)*total)))
  {
    dw_error(diagnostic,model,"memory","FFT coupling/caller storage exceeds 512 MiB budget");
    return 0;
  }
  if(iterative != 2) {
    system = calloc(total*total, sizeof(*system));
    original = calloc(total*total, sizeof(*original));
  }
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
  volume_chi = calloc(nc ? nc : 1, sizeof(*volume_chi));
  volume_size = calloc(nc ? nc : 1, sizeof(*volume_size));
  if(nw>SIZE_MAX-np || nw+np>SIZE_MAX/3)
  {
    dw_error(diagnostic,model,"size","conductor source stride overflow");
    goto failure;
  }
  source_stride=3*(nw+np);
  if(conductor_unknowns && source_stride>SIZE_MAX/conductor_unknowns)
  {
    dw_error(diagnostic,model,"size","conductor source state size overflow");
    goto failure;
  }
  conductor_sources=calloc(conductor_unknowns ?
      conductor_unknowns*source_stride:1,sizeof(*conductor_sources));
  if((iterative != 2 && (!system || !original)) || !rhs || !rhs0 || !conductor_matrix || !unit ||
     !conductor_current || !volume_current || !x || !y || !z || !length ||
     !radius || !px || !py || !pz || !parea || !volume_chi || !volume_size ||
     !conductor_sources)
  {
    dw_error(diagnostic, model, "memory", "unable to allocate coupled dense system");
    goto failure;
  }
  geometry_scaled = 1;
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
  dielectric_ground_initialize(frequency_mhz);
  fblock((int)conductor_unknowns, (int)conductor_unknowns,
         (int)(conductor_unknowns*conductor_unknowns), data.ipsym);
  {
  int background=gnd.ksymp;
  if(background!=1 && gnd.iperf==2 && gnd.zrati==1.0) gnd.ksymp=1;
  cmset((int)conductor_unknowns, conductor_matrix, rkh, extended_kernel);
  if(wave == NULL)
    etmns(0., 0., 0., 0., 0., 0., 0, rhs);
  else
    etmns(wave->theta_degrees*DW_PI/180.0,
          wave->phi_degrees*DW_PI/180.0,
          wave->eta_degrees*DW_PI/180.0, 0., 0., wave->axial_ratio,
          wave->type, rhs);
  gnd.ksymp=background;
  }
  if(system)
    for(row = 0; row < conductor_unknowns; ++row)
      for(column = 0; column < conductor_unknowns; ++column)
        /* NEC's factorizer transposes cmset() output before elimination. */
        system[row*total+column] = conductor_matrix[column+row*conductor_unknowns];

  vie_kernel_default_options(&kernel_options);
  kernel_options.relative_tolerance = 5.0e-7;
  kernel_options.absolute_tolerance = 5.0e-10;
  {
    double phase_started=performance_now();
    for(cell = 0; cell < nc; ++cell)
    {
      complex double incident[3] = {0.0, 0.0, 0.0};
      const dielectric_material *material=dw_material(model,mesh->cells[cell].material_id);
      complex double epsilon_r;
      int component;
      if(material==NULL)
      {
        performance_add("dielectric_volume_assembly",performance_now()-phase_started);
        dw_error(diagnostic,model,"material","unknown mesh material");goto failure;
      }
      epsilon_r=dielectric_material_epsilon(material,k0);
      volume_chi[cell]=epsilon_r-1.0;
      volume_size[cell]=mesh->cells[cell].size[0]*mesh->cells[cell].size[1]*
          mesh->cells[cell].size[2];
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
        }
        {
          complex double reflected[3]={0,0,0};
          if(gnd.ksymp!=1) dielectric_ground_reflect(direction,incident,k0,reflected);
          for(component=0;component<3;++component)
            incident[component]=incident[component]*phase+reflected[component]*average*
              cexp(-I*k0*(dot-2*direction[2]*mesh->cells[cell].center[2]));
        }
      }
      for(component = 0; component < 3; ++component)
      {
        if(system) system[(conductor_unknowns+3*cell+component)*total+
               conductor_unknowns+3*cell+component] = epsilon_r;
        rhs[conductor_unknowns+3*cell+component] =
            I*omega*DW_EPSILON_0*volume_chi[cell]*incident[component];
      }
    }
    performance_add("dielectric_volume_assembly",performance_now()-phase_started);
  }

  /* Bound the optional cache at 16 MiB; retain the same uncached operator if
     the cache is too large or allocation fails. Size is 9 complex/cell/wire. */
  if(nc && nw && gnd.ksymp!=1 && gnd.iperf==2 &&
     nc<=((size_t)16*1024*1024)/(9*sizeof(complex double))/nw)
  {
    wire_background_basis=calloc(nc*nw*9,sizeof(*wire_background_basis));
    if(wire_background_basis)
      for(cell=0;cell<nc;++cell)
        if(!dw_wire_background_basis(&mesh->cells[cell],k0,wire_background_basis+9*nw*cell))
        { dw_error(diagnostic,model,"cross-kernel","reflected wire basis evaluation failed"); goto failure; }
  }

  /* Volume rows: averaged NEC field produced by each explicit conductor
     source state. Processes own disjoint source/cell result blocks. */
  {
    double phase_started=performance_now();
    if(!dw_prepare_conductor_sources(nw,np,conductor_unknowns,unit,conductor_sources) ||
       !dw_assemble_wire_to_volume(mesh,nw,np,conductor_unknowns,
          conductor_sources,wire_background_basis,volume_chi,k0,omega,
          &wire_to_volume,&wire_to_volume_mapped,&coupling_workers_used,
          &coupling_fell_back))
    {
      performance_add("dielectric_wire_to_volume_assembly",
          performance_now()-phase_started);
      dw_error(diagnostic,model,"cross-kernel",
          "wire-to-volume interaction failed");
      goto failure;
    }
    performance_add("dielectric_wire_to_volume_assembly",
        performance_now()-phase_started);
  }
  free(conductor_sources);
  conductor_sources=NULL;
  fprintf(stream,"\n\n  CONDUCTOR-VOLUME COUPLING WORKERS: %zu%s",
      coupling_workers_used,coupling_fell_back ? " (SERIAL FALLBACK)":"");

  /* Fork-based coupling must finish before OpenMP/BLAS runtimes start. */
  {
    size_t numerical_threads=dw_coupling_processes;
    if(!concurrency_requested_threads() && !iterative &&
       !strcmp(dense_backend_name_for_order(total),"lapack") &&
       numerical_threads>4) numerical_threads=4;
    if(!concurrency_apply_threads(numerical_threads))
    {
      dw_error(diagnostic,model,"threads","unable to apply numerical thread budget");
      goto failure;
    }
  }

  /* FFTW may initialize a threaded runtime while creating its plans. Keep
     that work after all fork-based coupling workers have exited. */
  if(iterative == 2)
  {
    size_t reserve=sizeof(complex double)*(4*total+
        conductor_unknowns*conductor_unknowns+2*(nw+3*np))+
        sizeof(double)*(5*nw+4*np);
    size_t coupling_bytes=conductor_unknowns*nv*sizeof(complex double);
    if(reserve>SIZE_MAX-coupling_bytes)
    {
      dw_error(diagnostic,model,"memory","FFT application reserve overflow");
      goto failure;
    }
    reserve+=coupling_bytes;
    if(!(gnd.ksymp!=1 && gnd.iperf==1 ? dielectric_fft_prepare_pec :
         dielectric_fft_prepare)(&accelerated,model,mesh,conductor_unknowns,k0,
            reserve,diagnostic)) goto failure;
    performance_fft_backend(vie_fft_backend_requested_name(),
        dielectric_fft_backend(accelerated));
    fprintf(stream,"\n\n  FFT BACKEND: %s\n  FFT ACCOUNTED APPLICATION BYTES: %zu\n  FFT THREAD BUDGET: %zu",
        dielectric_fft_backend(accelerated),dielectric_fft_storage(accelerated),
        dielectric_fft_threads(accelerated));
    fflush(stream);
    for(row=0;row<conductor_unknowns;++row)
      for(column=0;column<conductor_unknowns;++column)
        *dw_entry(NULL,accelerated,total,row*total+column)=
            conductor_matrix[column+row*conductor_unknowns];
  }
  for(column=0;column<conductor_unknowns;++column)
    for(row=0;row<nv;++row)
      *dw_entry(system,accelerated,total,
          (conductor_unknowns+row)*total+column)=wire_to_volume[column*nv+row];
  dw_free_coupling(wire_to_volume,conductor_unknowns*nv,wire_to_volume_mapped);
  wire_to_volume=NULL;

  if(!accelerated)
  {
    double phase_started=performance_now();
    size_t failed_testing=SIZE_MAX,failed_source=SIZE_MAX;
    if(!dw_assemble_volume_pairs(mesh,k0,&kernel_options,volume_chi,volume_size,
         system,total,conductor_unknowns,&failed_testing,&failed_source))
    {
      char value[96];
      performance_add("dielectric_volume_assembly",performance_now()-phase_started);
      if(failed_testing==SIZE_MAX)
        dw_error(diagnostic,model,"memory","unable to allocate deterministic assembly status");
      else
      {
        snprintf(value,sizeof(value),"volume cells %zu,%zu",
            failed_testing+1,failed_source+1);
        dw_error(diagnostic,model,value,"volume-volume interaction failed");
      }
      goto failure;
    }
    performance_add("dielectric_volume_assembly",performance_now()-phase_started);
  }

  /* Wire rows: point-tested field of each constant cuboid current. */
  {
    double phase_started=performance_now();
    size_t tasks;
    if(nc && nw>SIZE_MAX/nc)
    {
      performance_add("dielectric_volume_to_wire_assembly",performance_now()-phase_started);
      dw_error(diagnostic,model,"size","wire-row task count overflow");goto failure;
    }
    tasks=nw*nc;
    unsigned char *failed=tasks ? calloc(tasks,sizeof(*failed)):NULL;
    ptrdiff_t task;
    if(tasks && failed==NULL)
    {
      performance_add("dielectric_volume_to_wire_assembly",performance_now()-phase_started);
      dw_error(diagnostic,model,"memory","unable to allocate deterministic wire-row status");
      goto failure;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) if(tasks >= 8)
#endif
    for(task=0;task<(ptrdiff_t)tasks;++task)
    {
      size_t wire=(size_t)task/nc,source_cell=(size_t)task%nc;
      double observation[3]={x[wire],y[wire],z[wire]};
      double tangent[3]={data.cab[wire],data.sab[wire],data.salp[wire]};
      complex double block[3][3];
      int p,q;
      if(!dw_point_cuboid_n(observation,&mesh->cells[source_cell],k0,block))
      {failed[task]=1;continue;}
      for(q=0;q<3;++q)
      {
        complex double tested=0.0;
        for(p=0;p<3;++p) tested+=tangent[p]*block[p][q];
        (*dw_entry(system,accelerated,total,wire*total+
            conductor_unknowns+3*source_cell+q))=
            tested/(I*omega*DW_EPSILON_0);
      }
    }
    for(task=0;task<(ptrdiff_t)tasks;++task) if(failed[task])
    {
      char value[96];size_t wire=(size_t)task/nc,source_cell=(size_t)task%nc;
      free(failed);
      performance_add("dielectric_volume_to_wire_assembly",performance_now()-phase_started);
      snprintf(value,sizeof(value),"wire %zu, volume cell %zu",wire+1,source_cell+1);
      dw_error(diagnostic,model,value,"volume-to-wire interaction failed");goto failure;
    }
    free(failed);
    performance_add("dielectric_volume_to_wire_assembly",performance_now()-phase_started);
  }

  /* Patch MFIE rows: point-tested H = integral grad(G) x Jp dv.  NEC's
     first and second patch equations test -t2.H and -t1.H respectively. */
  {
    double phase_started=performance_now();
    size_t tasks;
    if(nc && np>SIZE_MAX/nc)
    {
      performance_add("dielectric_volume_to_patch_assembly",performance_now()-phase_started);
      dw_error(diagnostic,model,"size","patch-row task count overflow");goto failure;
    }
    tasks=np*nc;
    unsigned char *failed=tasks ? calloc(tasks,sizeof(*failed)):NULL;
    ptrdiff_t task;
    if(tasks && failed==NULL)
    {
      performance_add("dielectric_volume_to_patch_assembly",performance_now()-phase_started);
      dw_error(diagnostic,model,"memory","unable to allocate deterministic patch-row status");
      goto failure;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) if(tasks >= 8)
#endif
    for(task=0;task<(ptrdiff_t)tasks;++task)
    {
      size_t patch=(size_t)task/nc,source_cell=(size_t)task%nc;
      double observation[3]={px[patch],py[patch],pz[patch]};
      double t1[3]={data.t1x[patch],data.t1y[patch],data.t1z[patch]};
      double t2[3]={data.t2x[patch],data.t2y[patch],data.t2z[patch]};
      complex double block[3][3];
      int p,q;
      if(!dw_point_cuboid_h(observation,&mesh->cells[source_cell],k0,block))
      {failed[task]=1;continue;}
      for(q=0;q<3;++q)
      {
        complex double tested_t1=0.0,tested_t2=0.0;
        for(p=0;p<3;++p)
        {
          tested_t1+=t1[p]*block[p][q];tested_t2+=t2[p]*block[p][q];
        }
        (*dw_entry(system,accelerated,total,(nw+2*patch)*total+
            conductor_unknowns+3*source_cell+q))=-data.psalp[patch]*tested_t2;
        (*dw_entry(system,accelerated,total,(nw+2*patch+1)*total+
            conductor_unknowns+3*source_cell+q))=-data.psalp[patch]*tested_t1;
      }
    }
    for(task=0;task<(ptrdiff_t)tasks;++task) if(failed[task])
    {
      char value[96];size_t patch=(size_t)task/nc,source_cell=(size_t)task%nc;
      free(failed);
      performance_add("dielectric_volume_to_patch_assembly",performance_now()-phase_started);
      snprintf(value,sizeof(value),"patch %zu, volume cell %zu",patch+1,source_cell+1);
      dw_error(diagnostic,model,value,"volume-to-patch interaction failed");goto failure;
    }
    free(failed);
    performance_add("dielectric_volume_to_patch_assembly",performance_now()-phase_started);
  }

  for(row = conductor_unknowns; row < total; ++row)
  {
    for(column = 0; column < nw; ++column)
      wire_to_volume_norm += creal((*dw_entry(system,accelerated,total,row*total+column))*
                                   conj((*dw_entry(system,accelerated,total,row*total+column))));
    for(column = nw; column < conductor_unknowns; ++column)
      patch_to_volume_norm += creal((*dw_entry(system,accelerated,total,row*total+column))*
                                    conj((*dw_entry(system,accelerated,total,row*total+column))));
  }
  for(row = 0; row < nw; ++row)
    for(column = conductor_unknowns; column < total; ++column)
      volume_to_wire_norm += creal((*dw_entry(system,accelerated,total,row*total+column))*
                                   conj((*dw_entry(system,accelerated,total,row*total+column))));
  for(row = nw; row < conductor_unknowns; ++row)
    for(column = conductor_unknowns; column < total; ++column)
      volume_to_patch_norm += creal((*dw_entry(system,accelerated,total,row*total+column))*
                                    conj((*dw_entry(system,accelerated,total,row*total+column))));
  wire_to_volume_norm=sqrt(wire_to_volume_norm);
  patch_to_volume_norm=sqrt(patch_to_volume_norm);
  volume_to_wire_norm=sqrt(volume_to_wire_norm);
  volume_to_patch_norm=sqrt(volume_to_patch_norm);

  if(!accelerated) memcpy(original, system, total*total*sizeof(*system));
  memcpy(rhs0, rhs, total*sizeof(*rhs));
  if(iterative)
  {
    memset(rhs, 0, total*sizeof(*rhs));
    int converged;
    if(accelerated) {
      vie_linear_operator op;
      dielectric_fft_get_operator(accelerated,&op);
      converged=vie_gmres(&op,rhs0,rhs,&dielectric_gmres_settings,&report);
    } else converged=vie_dense_gmres(original,total,rhs0,rhs,&dielectric_gmres_settings,&report);
    relative_residual=report.relative_residual;
    fprintf(stream, "\n\n  SOLVER: %s"
                    "\n  GMRES STATUS: %s"
                    "\n  GMRES ITERATIONS: %zu"
                    "\n  GMRES RESTARTS: %zu"
                    "\n  GMRES OPERATOR APPLICATIONS: %zu"
                    "\n  GMRES TRUE RELATIVE RESIDUAL: %.12E",
            accelerated ? "FFT-OPERATOR GMRES (UNIFORM CUBES; 512 MIB CAP)":"DENSE-OPERATOR GMRES (EXPERIMENTAL; NO ACCELERATION)",
            vie_iterative_status_name(report.status),
            report.iterations, report.restarts, report.operator_applications,
            report.relative_residual);
    if(!converged)
    {
      dw_error(diagnostic, model, "GMRES", vie_iterative_status_name(report.status));
      goto failure;
    }
  }
  else
  {
  performance_dense_backend(total);
  fprintf(stream, "\n\n  DENSE BACKEND: %s (REQUESTED %s)",
          dense_backend_name_for_order(total), dense_backend_requested_name());
  if(!dense_backend_solve_row_major(system, rhs, total))
  { dw_error(diagnostic, model, "matrix", "singular coupled wire-volume system"); goto failure; }
  }
  if(!accelerated) {
  for(row = 0; row < total; ++row)
  {
    complex double value = -rhs0[row];
    rhs_norm += creal(rhs0[row]*conj(rhs0[row]));
    for(column = 0; column < total; ++column)
      value += original[row*total+column]*rhs[column];
    relative_residual += creal(value*conj(value));
  }
  relative_residual = sqrt(relative_residual/fmax(rhs_norm, DBL_MIN));
  }
  {
    result_solver solver;
    memset(&solver,0,sizeof(solver));
    solver.requested=iterative==2?"fft":(iterative?"gmres":"dense");
    solver.used=iterative?solver.requested:dense_backend_solver_name_for_order(total);
    solver.iterative_fields_defined=1;
    solver.relative_tolerance=dielectric_gmres_settings.relative_tolerance;
    solver.maximum_iterations=dielectric_gmres_settings.maximum_iterations;
    solver.iterations=iterative?report.iterations:0;
    solver.reported_residual=relative_residual;solver.true_residual=relative_residual;
    solver.converged=1;
    solver.accounted_storage_bytes=accelerated?dielectric_fft_storage(accelerated):total*total*sizeof(complex double);
    solver.elapsed_seconds=(double)(clock()-started)/CLOCKS_PER_SEC;
    result_json_set_solver(&solver);result_json_add_material_evaluations(model,frequency_mhz);
  }
  memcpy(volume_current, rhs+conductor_unknowns, nv*sizeof(*volume_current));
  memcpy(conductor_current, rhs, conductor_unknowns*sizeof(*conductor_current));
  cabc(conductor_current);
  if(np != 0)
    memcpy(crnt.cur+nw, conductor_current+nw, 3*np*sizeof(*conductor_current));

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
    {
      double center[3]={x[row],y[row],z[row]};
      result_json_add_wire(data.itag[row],row+1,center,length[row],current);
    }
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
    result_json_add_feed(row,data.itag[source],source+1,voltage,current,
        impedance,1,0.5*creal(voltage*conj(current)));
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
      {
        double center[3]={px[row],py[row],pz[row]};
        double t1[3]={data.t1x[row],data.t1y[row],data.t1z[row]};
        double t2[3]={data.t2x[row],data.t2y[row],data.t2z[row]};
        result_json_add_patch(row+1,center,parea[row],t1,t2,conductor_current+offset);
      }
    }
  }
  fprintf(stream, "\n\n  COUPLED VOLUME POLARIZATION CURRENTS (A/M^2)"
                  "\n  CELL REGION MATERIAL          JX(RE,IM)                 JY(RE,IM)                 JZ(RE,IM)");
  for(cell = 0; cell < nc; ++cell)
  {
    fprintf(stream, "\n %5zu %6d %8d  %19.12E %19.12E  %19.12E %19.12E  %19.12E %19.12E",
            cell+1, mesh->cells[cell].region_id, mesh->cells[cell].material_id,
            creal(volume_current[3*cell]), cimag(volume_current[3*cell]),
            creal(volume_current[3*cell+1]), cimag(volume_current[3*cell+1]),
            creal(volume_current[3*cell+2]), cimag(volume_current[3*cell+2]));
    result_json_add_volume(cell+1,mesh->cells[cell].box_id,
        mesh->cells[cell].region_id,mesh->cells[cell].material_id,
        mesh->cells[cell].center,mesh->cells[cell].size,volume_current+3*cell);
  }

  if(pattern != NULL && pattern->near_enabled &&
     !dielectric_print_near_fields(mesh, volume_current, k0, 1, stream))
  {
    dw_error(diagnostic, model, "NE/NH",
        "near-field evaluation requires positive counts, at least one cell longest-edge clearance from every dielectric cell, and converged background quadrature");
    goto failure;
  }
  {
    double absorbed=dielectric_absorbed_power(model,mesh,volume_current,k0);
  fprintf(stream,"\n\n  DIELECTRIC ABSORBED POWER (W): %19.12E",absorbed);
  if(wave == NULL)
  {
    double accepted = 0.0;
    double radiated = dw_radiated_power(mesh, volume_current, k0, 3896);
    double coarse = dw_radiated_power(mesh, volume_current, k0, 974);
    for(row = 0; row < (size_t)vsorc.nsant; ++row)
      accepted += 0.5*creal(vsorc.vsant[row]*
          conj(conductor_current[vsorc.isant[row]-1]*data.wlam));
    if(gnd.ksymp!=1 && gnd.iperf==2)
    {
      double ground_power,ground_delta;
      if(!dielectric_ground_power(mesh,volume_current,k0,1,&ground_power,&ground_delta))
      { dw_error(diagnostic,model,"ground power","lower-half-space flux quadrature did not converge"); goto failure; }
      fprintf(stream,"\n\n  COUPLED REAL-GROUND POWER (W)"
                     "\n  ACCEPTED POWER: %19.12E"
                     "\n  UPPER-HEMISPHERE RADIATED POWER: %19.12E"
                     "\n  LOWER-HALF-SPACE POWER: %19.12E"
                     "\n  RELATIVE POWER DEFECT: %19.12E"
                     "\n  ANGULAR QUADRATURE CHANGE: %19.12E"
                     "\n  GROUND FLUX QUADRATURE CHANGE (W): %19.12E",
              accepted,radiated,ground_power,fabs(accepted-radiated-ground_power-absorbed)/fmax(fabs(accepted),DBL_MIN),
              fabs(radiated-coarse)/fmax(radiated,DBL_MIN),ground_delta);
      {
        result_power power;double defect=fabs(accepted-radiated-ground_power-absorbed);
        memset(&power,0,sizeof(power));power.basis="driven-accepted-power";power.accepted_defined=1;
        power.air_radiated_defined=power.lower_flux_defined=power.dielectric_absorbed_defined=power.defect_defined=1;
        power.accepted_w=accepted;power.air_radiated_w=radiated;power.lower_flux_w=ground_power;
        power.dielectric_absorbed_w=absorbed;power.defect_w=defect;
        power.defect_relative=defect/fmax(fabs(accepted),DBL_MIN);power.theta_count=3896;power.phi_count=1;
        result_json_set_power(&power);
      }
    }
    else
    {
    fprintf(stream, "\n\n  COUPLED DIELECTRIC POWER (W)"
                    "\n  ACCEPTED POWER: %19.12E"
                    "\n  RADIATED POWER: %19.12E"
                    "\n  RELATIVE POWER DEFECT: %19.12E"
                    "\n  ANGULAR QUADRATURE CHANGE: %19.12E",
            accepted, radiated, fabs(accepted-radiated-absorbed)/fmax(fabs(accepted), DBL_MIN),
            fabs(radiated-coarse)/fmax(radiated, DBL_MIN));
      {
        result_power power;double defect=fabs(accepted-radiated-absorbed);
        memset(&power,0,sizeof(power));power.basis="driven-accepted-power";power.accepted_defined=1;
        power.air_radiated_defined=power.dielectric_absorbed_defined=power.defect_defined=1;
        power.accepted_w=accepted;power.air_radiated_w=radiated;power.dielectric_absorbed_w=absorbed;
        if(gnd.ksymp!=1&&gnd.iperf==1){power.lower_flux_defined=1;power.lower_flux_w=0.0;}
        power.defect_w=defect;power.defect_relative=defect/fmax(fabs(accepted),DBL_MIN);
        power.theta_count=3896;power.phi_count=1;result_json_set_power(&power);
      }
    }
  }
  if(wave!=NULL && (gnd.ksymp==1 || gnd.iperf==1))
  {
    double extinction=dielectric_ground_extinction(mesh,volume_current,k0,wave,1);
    double radiated=dw_radiated_power(mesh,volume_current,k0,3896);
    double coarse=dw_radiated_power(mesh,volume_current,k0,974);
    fprintf(stream,"\n\n  COUPLED FREE/PEC PLANE-WAVE POWER (W)"
                   "\n  EXTINCTION POWER: %19.12E"
                   "\n  %s: %19.12E"
                   "\n  RELATIVE POWER DEFECT: %19.12E"
                   "\n  ANGULAR QUADRATURE CHANGE: %19.12E",
            extinction,gnd.ksymp==1 ? "TOTAL SCATTERED POWER":"UPPER-HEMISPHERE SCATTERED POWER",
            radiated,fabs(extinction-radiated-absorbed)/fmax(fabs(extinction),DBL_MIN),
            fabs(radiated-coarse)/fmax(radiated,DBL_MIN));
    {
      result_power power;double defect=fabs(extinction-radiated-absorbed);
      memset(&power,0,sizeof(power));power.basis="plane-wave-scattered-balance";
      power.air_radiated_defined=power.dielectric_absorbed_defined=power.defect_defined=1;
      power.air_radiated_w=radiated;power.dielectric_absorbed_w=absorbed;
      if(gnd.ksymp!=1&&gnd.iperf==1){power.lower_flux_defined=1;power.lower_flux_w=0.0;}
      power.defect_w=defect;power.defect_relative=defect/fmax(fabs(extinction),DBL_MIN);
      power.theta_count=3896;power.phi_count=1;result_json_set_power(&power);
    }
  }
  if(wave!=NULL && gnd.ksymp!=1 && gnd.iperf==2)
  {
    double extinction=dielectric_ground_extinction(mesh,volume_current,k0,wave,1);
    double radiated=dw_radiated_power(mesh,volume_current,k0,3896),ground_power,ground_delta;
    if(!dielectric_ground_power(mesh,volume_current,k0,1,&ground_power,&ground_delta))
    { dw_error(diagnostic,model,"ground power","lower-half-space flux quadrature did not converge"); goto failure; }
    fprintf(stream,"\n\n  COUPLED REAL-GROUND PLANE-WAVE POWER (W)"
                   "\n  EXTINCTION POWER: %19.12E"
                   "\n  UPPER-HEMISPHERE SCATTERED POWER: %19.12E"
                   "\n  LOWER-HALF-SPACE SCATTERED POWER: %19.12E"
                   "\n  RELATIVE POWER DEFECT: %19.12E"
                   "\n  GROUND FLUX QUADRATURE CHANGE (W): %19.12E",
            extinction,radiated,ground_power,fabs(extinction-radiated-ground_power-absorbed)/fmax(fabs(extinction),DBL_MIN),ground_delta);
    {
      result_power power;double defect=fabs(extinction-radiated-ground_power-absorbed);
      memset(&power,0,sizeof(power));power.basis="plane-wave-scattered-balance";
      power.air_radiated_defined=power.lower_flux_defined=power.dielectric_absorbed_defined=power.defect_defined=1;
      power.air_radiated_w=radiated;power.lower_flux_w=ground_power;power.dielectric_absorbed_w=absorbed;
      power.defect_w=defect;power.defect_relative=defect/fmax(fabs(extinction),DBL_MIN);
      power.theta_count=3896;power.phi_count=1;result_json_set_power(&power);
    }
  }
  } /* dielectric absorption in each power balance */
  if(pattern != NULL && pattern->enabled)
  {
    size_t it, ip;
    fprintf(stream, "\n\n  COUPLED SCATTERED FAR FIELD (R TIMES E, V)"
                    "\n  THETA      PHI          ETHETA(RE,IM)              EPHI(RE,IM)");
    result_json_far_begin((size_t)result_json_current_request(),
        wave ? "plane-wave-scattered" : "driven-total");
    for(ip = 0; ip < pattern->phi_count; ++ip)
      for(it = 0; it < pattern->theta_count; ++it)
      {
        double td = pattern->theta_start_degrees+it*pattern->theta_step_degrees;
        double pd = pattern->phi_start_degrees+ip*pattern->phi_step_degrees;
        complex double et, ep;
        dw_far_field(mesh, volume_current, k0, td*DW_PI/180.0, pd*DW_PI/180.0,
                     &et, &ep);
        fprintf(stream, "\n %8.3f %8.3f %19.12E %19.12E %19.12E %19.12E",
                td, pd, creal(et), cimag(et), creal(ep), cimag(ep));
        result_json_far_sample(td,pd,et,ep);
      }
  }

  dw_restore_geometry(x, y, z, length, radius);
  dw_restore_patches(px, py, pz, parea); restored = 1;
  free(parea); free(pz); free(py); free(px);
  free(radius); free(length); free(z); free(y); free(x); free(volume_current);
  free(volume_size); free(volume_chi);
  dw_free_coupling(wire_to_volume,conductor_unknowns*nv,wire_to_volume_mapped);
  free(conductor_sources);
  free(wire_background_basis);
  free(conductor_current); free(unit); free(conductor_matrix);
  dielectric_fft_free(accelerated);
  free(rhs0); free(rhs); free(original); free(system);
  return 1;

failure:
  if(geometry_scaled && !restored && x && y && z && length && radius)
    dw_restore_geometry(x, y, z, length, radius);
  if(geometry_scaled && !restored && px && py && pz && parea)
    dw_restore_patches(px, py, pz, parea);
  free(parea); free(pz); free(py); free(px);
  free(radius); free(length); free(z); free(y); free(x); free(volume_current);
  free(volume_size); free(volume_chi);
  dw_free_coupling(wire_to_volume,conductor_unknowns*nv,wire_to_volume_mapped);
  free(conductor_sources);
  free(wire_background_basis);
  free(conductor_current); free(unit); free(conductor_matrix);
  dielectric_fft_free(accelerated);
  free(rhs0); free(rhs); free(original); free(system);
  return 0;
}
