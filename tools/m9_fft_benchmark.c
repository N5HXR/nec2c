#define _POSIX_C_SOURCE 200809L
#include "../vie_fft.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

static double seconds(void)
{struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
typedef struct {
  vie_linear_operator original;
  size_t applications;
  double started;
} progress_context;
static int progress_apply(void *context,const complex double *x,complex double *y)
{
  progress_context *p=context;
  int ok=p->original.apply(p->original.context,x,y);
  p->applications++;
  if(p->applications%10==0) {
    fprintf(stderr,"progress: matvecs=%zu elapsed_seconds=%.3f\n",p->applications,seconds()-p->started);
    fflush(stderr);
  }
  return ok;
}
/* Integer bounds at 0.25 m, relative to the full-tree bounding box. */
static const size_t boxes[4][6]={
  {11,11,0,13,13,32},{13,11,22,24,13,24},
  {11,13,26,13,24,28},{0,0,32,24,24,48}};
static double sinc(double t) {return fabs(t)<1e-12 ? 1:sin(t)/t;}
/* Independent angular integration of radiated power, in units of
   (omega*epsilon0/2). Plane-wave extinction uses the same normalization. */
static double scattered(const size_t *xyz,size_t cells,const complex double *p,
                         double h,double k,int order)
{
  const double pi=3.14159265358979323846;
  double integral=0;int m,a,j,d;size_t c;
  for(m=0;m<order;m++) {
    double mu=cos(pi*(m+.75)/(order+.5)),old,derivative=0,weight;
    do {
      double pn=1,pm=0;
      old=mu;
      for(j=1;j<=order;j++) {double prev=pm;pm=pn;pn=((2*j-1)*mu*pm-(j-1)*prev)/j;}
      derivative=order*(mu*pn-pm)/(mu*mu-1);mu=old-pn/derivative;
    } while(fabs(mu-old)>2e-15);
    weight=2/((1-mu*mu)*derivative*derivative);
    for(a=0;a<2*order;a++) {
      double phi=pi*a/order,dir[3]={sqrt(1-mu*mu)*cos(phi),sqrt(1-mu*mu)*sin(phi),mu};
      double cell_weight=h*h*h*sinc(k*h*dir[0]/2)*sinc(k*h*dir[1]/2)*sinc(k*h*dir[2]/2),norm=0;
      complex double f[3]={0,0,0},longitudinal=0;
      for(c=0;c<cells;c++) {
        double dot=0;complex double phase;
        for(d=0;d<3;d++) dot+=dir[d]*(xyz[3*c+d]+.5);
        phase=cell_weight*cexp(I*k*h*dot);
        for(d=0;d<3;d++) f[d]+=p[3*c+d]*phase;
      }
      for(d=0;d<3;d++) {norm+=creal(f[d]*conj(f[d]));longitudinal+=dir[d]*f[d];}
      integral+=weight*pi/order*(norm-creal(longitudinal*conj(longitudinal)));
    }
  }
  return k*k*k*integral/(16*pi*pi);
}
int main(int argc,char **argv)
{
  const size_t cap=512u*1024u*1024u;
  size_t shape[3]={24,24,48},count=0,*xyz,c=0,x,y,z,n,workspace,caller;
  int full=1,nb,b,d,ok;double tol=1e-8,start,build,matvec,solve;
  const double h=.25,k=2*3.14159265358979323846*28e6/299792458.0;
  complex double *rhs,*solution,*product,sum[3]={0,0,0};
  vie_fft_plan *plan=NULL;vie_linear_operator op;progress_context progress;vie_gmres_options options;vie_gmres_report report;struct rusage usage;
  if(argc>4) return 2;
  if(argc>1 && strcmp(argv[1],"skeleton")==0) full=0;
  else if(argc>1 && strcmp(argv[1],"tree")!=0) {fprintf(stderr,"usage: %s [tree|skeleton] [rtol] [max_iterations]\n",argv[0]);return 2;}
  if(argc>2) {char *end;tol=strtod(argv[2],&end);if(*end || !isfinite(tol) || tol<=0 || tol>=1) return 2;}
  vie_gmres_default_options(&options);options.relative_tolerance=tol;
  if(argc>3) {char *end;unsigned long v=strtoul(argv[3],&end,10);if(*end || !v || v>10000) return 2;options.maximum_iterations=v;}
  nb=full ? 4:3;
  if(!full) {shape[0]=shape[1]=13;shape[2]=32;}
  for(b=0;b<nb;b++) count+=(boxes[b][3]-boxes[b][0])*(boxes[b][4]-boxes[b][1])*(boxes[b][5]-boxes[b][2]);
  n=3*count;
  workspace=16*((2*options.restart+3)*n+(options.restart+1)*options.restart+3*options.restart+1)+8*options.restart;
  caller=3*count*sizeof(size_t)+3*n*sizeof(complex double);
  if(workspace+caller>=cap) return 2;
  xyz=malloc(3*count*sizeof(*xyz));rhs=calloc(n,sizeof(*rhs));solution=calloc(n,sizeof(*solution));product=calloc(n,sizeof(*product));
  if(!xyz || !rhs || !solution || !product) return 2;
  for(b=0;b<nb;b++) for(x=boxes[b][0];x<boxes[b][3];x++)
    for(y=boxes[b][1];y<boxes[b][4];y++) for(z=boxes[b][2];z<boxes[b][5];z++) {
      xyz[3*c]=x-(full ? 0:11);xyz[3*c+1]=y-(full ? 0:11);xyz[3*c+2]=z;c++;
    }
  /* Solve for J/(i omega epsilon0) with x-polarized, +z incident plane wave.
     Exact pulse-cell average; phase origin is the grid's lower corner. */
  for(c=0;c<count;c++) rhs[3*c]=2*cexp(-I*k*h*(xyz[3*c+2]+.5))*sin(k*h/2)/(k*h/2);
  printf("model=%s cells=%zu unknowns=%zu rtol=%.3g budget_bytes=%zu\n",full ? "tree":"skeleton",count,n,tol,cap);fflush(stdout);
  start=seconds();
  if(!vie_fft_create(&plan,shape,count,xyz,h,k,3,cap-workspace-caller)) {fprintf(stderr,"FFT plan rejected or build failed\n");return 1;}
  build=seconds()-start;vie_fft_operator(plan,&op);
  progress.original=op;progress.applications=0;progress.started=seconds();
  op.apply=progress_apply;op.context=&progress;
  printf("build_seconds=%.6f plan_bytes=%zu gmres_bytes=%zu caller_bytes=%zu\n",build,vie_fft_bytes(plan),workspace,caller);fflush(stdout);
  start=seconds();op.apply(op.context,rhs,product);matvec=seconds()-start;
  start=seconds();ok=vie_gmres(&op,rhs,solution,&options,&report);solve=seconds()-start;
  for(c=0;c<count;c++) for(d=0;d<3;d++) sum[d]+=solution[3*c+d]*h*h*h;
  getrusage(RUSAGE_SELF,&usage);
  printf("matvec_seconds=%.6f solve_seconds=%.6f peak_rss_kib=%ld iterations=%zu applications=%zu residual=%.12g status=%s\n",matvec,solve,usage.ru_maxrss,report.iterations,report.operator_applications,report.relative_residual,vie_iterative_status_name(report.status));
  fflush(stdout);
  if(!ok) goto cleanup;
  for(d=0;d<3;d++) printf("integrated_normalized_current_%d=%.12g%+.12gj\n",d,creal(sum[d]),cimag(sum[d]));
  /* Complex transverse far-field sums, omitting common i*omega*eps factor. */
  for(b=0;b<3;b++) {
    complex double f[3]={0,0,0};
    for(c=0;c<count;c++) {
      complex double phase=cexp(I*k*h*(xyz[3*c+b]+.5))*h*h*h*sin(k*h/2)/(k*h/2);
      for(d=0;d<3;d++) if(d!=b) f[d]+=solution[3*c+d]*phase;
    }
    for(d=0;d<3;d++) printf("far_axis%d_component%d=%.12g%+.12gj\n",b,d,creal(f[d]),cimag(f[d]));
  }
  if(ok) {
    double extinction=0,sca16,sca32;
    for(c=0;c<count;c++) extinction-=cimag(solution[3*c]*conj(rhs[3*c]/2))*h*h*h;
    sca16=scattered(xyz,count,solution,h,k,16);sca32=scattered(xyz,count,solution,h,k,32);
    printf("extinction_normalized=%.12g scattered16_normalized=%.12g scattered32_normalized=%.12g angular_relative_change=%.12g power_relative_defect=%.12g\n",extinction,sca16,sca32,fabs(sca16-sca32)/fmax(fabs(sca32),1e-300),fabs(extinction-sca32)/fmax(fabs(extinction),1e-300));
  }
cleanup:
  vie_fft_destroy(plan);free(xyz);free(rhs);free(solution);free(product);return ok ? 0:1;
}
