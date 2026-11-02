#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "vie_fft.h"
#include "vie_kernel.h"
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_FFTW
#include <fftw3.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

struct vie_fft_plan {
  size_t cells, padded[3], grids, bytes, threads;
  int backend;
  size_t *map;
  complex double *data, *roots[3], self[3], *epsilon;
#ifdef HAVE_FFTW
  fftw_plan kernel_forward, apply_forward, apply_inverse;
  int fftw_runtime;
#endif
};
enum {FFT_BACKEND_AUTO,FFT_BACKEND_BUILTIN,FFT_BACKEND_FFTW};
#define FFTW_AUTO_MIN_GRIDS 512U
static int requested_backend=FFT_BACKEND_AUTO;
static int last_used_backend=-1;
#ifdef HAVE_FFTW
static size_t fftw_runtime_users;
static void fftw_release(vie_fft_plan *p)
{
  if(p->kernel_forward) {fftw_destroy_plan(p->kernel_forward);p->kernel_forward=NULL;}
  if(p->apply_forward) {fftw_destroy_plan(p->apply_forward);p->apply_forward=NULL;}
  if(p->apply_inverse) {fftw_destroy_plan(p->apply_inverse);p->apply_inverse=NULL;}
  if(p->fftw_runtime && fftw_runtime_users) {
    p->fftw_runtime=0;
    if(--fftw_runtime_users==0) fftw_cleanup_threads();
  }
}
#endif
static const double pi = 3.141592653589793238462643383279502884;
static int multiply(size_t a, size_t b, size_t *out)
{ if(b && a > SIZE_MAX/b) return 0; *out=a*b; return 1; }
static int add(size_t a, size_t b, size_t *out)
{ if(a > SIZE_MAX-b) return 0; *out=a+b; return 1; }
static size_t index3(const vie_fft_plan *p,size_t x,size_t y,size_t z)
{ return (x*p->padded[1]+y)*p->padded[2]+z; }
static int tensor_index(int i,int j)
{ if(i==j) return i; if(i>j) {int t=i;i=j;j=t;} return i==0 ? (j==1 ? 3:4):5; }
static size_t configured_threads(int backend)
{
#ifdef _OPENMP
  int count=omp_get_max_threads();
  (void)backend;
  return count>0 ? (size_t)count:1;
#else
  const char *text=getenv("OMP_NUM_THREADS");char *end;unsigned long count;
  if(backend!=FFT_BACKEND_FFTW) return 1;
  if(!text || !*text) return 1;
  count=strtoul(text,&end,10);
  return !*end && count>0 && count<=(unsigned long)INT_MAX ? (size_t)count:1;
#endif
}
int vie_fft_backend_set(const char *name)
{
  if(!name || !strcmp(name,"auto")) requested_backend=FFT_BACKEND_AUTO;
  else if(!strcmp(name,"builtin")) requested_backend=FFT_BACKEND_BUILTIN;
  else if(!strcmp(name,"fftw")) {
    if(!vie_fft_fftw_available()) return 0;
    requested_backend=FFT_BACKEND_FFTW;
  } else return 0;
  return 1;
}
int vie_fft_fftw_available(void)
{
#ifdef HAVE_FFTW
  return 1;
#else
  return 0;
#endif
}
const char *vie_fft_backend_requested_name(void)
{return requested_backend==FFT_BACKEND_BUILTIN ? "builtin":requested_backend==FFT_BACKEND_FFTW ? "fftw":"auto";}
const char *vie_fft_backend_name(const vie_fft_plan *p)
{return p && p->backend==FFT_BACKEND_FFTW ? "fftw":"builtin";}
const char *vie_fft_backend_last_used_name(void)
{return last_used_backend<0 ? NULL:last_used_backend==FFT_BACKEND_FFTW ? "fftw":"builtin";}
void vie_fft_backend_clear_last_used(void) {last_used_backend=-1;}
static void fft_line(complex double *a,size_t n,size_t stride,
                     const complex double *roots,int inverse)
{
  size_t i,j=0,len;
  for(i=1;i<n;i++) {
    size_t bit=n>>1;
    for(;j&bit;bit>>=1) j^=bit;
    j^=bit;
    if(i<j) {complex double t=a[i*stride];a[i*stride]=a[j*stride];a[j*stride]=t;}
  }
  for(len=2;len<=n;len*=2) {
    size_t base,k,half=len/2;
    for(base=0;base<n;base+=len) for(k=0;k<half;k++) {
      complex double w=roots[k*(n/len)],u=a[(base+k)*stride],v;
      if(inverse) w=conj(w);
      v=w*a[(base+k+half)*stride];
      a[(base+k)*stride]=u+v;a[(base+k+half)*stride]=u-v;
    }
    if(len==n) break;
  }
  if(inverse) for(i=0;i<n;i++) a[i*stride]/=(double)n;
}

static void fft3_serial(vie_fft_plan *p,complex double *a,int inverse)
{
  size_t x,y,z,nx=p->padded[0],ny=p->padded[1],nz=p->padded[2];
  for(x=0;x<nx;x++) for(y=0;y<ny;y++)
    fft_line(a+index3(p,x,y,0),nz,1,p->roots[2],inverse);
  for(x=0;x<nx;x++) for(z=0;z<nz;z++)
    fft_line(a+index3(p,x,0,z),ny,nz,p->roots[1],inverse);
  for(y=0;y<ny;y++) for(z=0;z<nz;z++)
    fft_line(a+index3(p,0,y,z),nx,ny*nz,p->roots[0],inverse);
}

/* Each transform line owns disjoint elements within a phase.  Parallelizing
   lines and tensor components therefore preserves fft_line's arithmetic order
   exactly; the implicit phase barriers preserve the serial 3-D dependency. */
static void fft3_batch(vie_fft_plan *p,complex double *a,size_t components,
                       int inverse)
{
  size_t nx=p->padded[0],ny=p->padded[1],nz=p->padded[2];
  size_t xy=nx*ny,xz=nx*nz,yz=ny*nz;
  ptrdiff_t line;size_t component;
#ifdef HAVE_FFTW
  if(p->backend==FFT_BACKEND_FFTW) {
    fftw_plan plan=components==6 ? p->kernel_forward:
        (inverse ? p->apply_inverse:p->apply_forward);
    size_t count=components*p->grids,g;
    fftw_execute_dft(plan,(fftw_complex *)a,(fftw_complex *)a);
    if(inverse) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(count >= 16384 && p->threads > 1) num_threads(p->threads)
#endif
      for(g=0;g<count;g++) a[g]/=(double)p->grids;
    }
    return;
  }
#endif
  if(p->threads==1 || p->grids<16384) {
    for(component=0;component<components;component++)
      fft3_serial(p,a+component*p->grids,inverse);
    return;
  }
#ifdef _OPENMP
#pragma omp parallel num_threads(p->threads)
#endif
  {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(line=0;line<(ptrdiff_t)(components*xy);line++) {
      size_t component=(size_t)line/xy,local=(size_t)line%xy;
      size_t x=local/ny,y=local%ny;
      fft_line(a+component*p->grids+index3(p,x,y,0),nz,1,
          p->roots[2],inverse);
    }
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(line=0;line<(ptrdiff_t)(components*xz);line++) {
      size_t component=(size_t)line/xz,local=(size_t)line%xz;
      size_t x=local/nz,z=local%nz;
      fft_line(a+component*p->grids+index3(p,x,0,z),ny,nz,
          p->roots[1],inverse);
    }
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(line=0;line<(ptrdiff_t)(components*yz);line++) {
      size_t component=(size_t)line/yz,local=(size_t)line%yz;
      size_t y=local/nz,z=local%nz;
      fft_line(a+component*p->grids+index3(p,0,y,z),nx,ny*nz,
          p->roots[0],inverse);
    }
  }
}
void vie_fft_destroy(vie_fft_plan *p)
{ if(p) {int d;
#ifdef HAVE_FFTW
  fftw_release(p);
#endif
  for(d=0;d<3;d++) free(p->roots[d]);
  free(p->map);free(p->data);free(p->epsilon);free(p);
} }
size_t vie_fft_bytes(const vie_fft_plan *p) {return p ? p->bytes:0;}
size_t vie_fft_threads(const vie_fft_plan *p) {return p ? p->threads:0;}
int vie_fft_parallel_available(void)
{
  if(requested_backend==FFT_BACKEND_BUILTIN) {
#ifdef _OPENMP
    return 1;
#else
    return 0;
#endif
  }
  if(requested_backend==FFT_BACKEND_FFTW) return vie_fft_fftw_available();
#if defined(_OPENMP) || defined(HAVE_FFTW)
  return 1;
#else
  return 0;
#endif
}
#ifdef HAVE_FFTW
static int fftw_plans(vie_fft_plan *p)
{
  int dimensions[3],d;
  for(d=0;d<3;d++) {
    if(p->padded[d]>(size_t)INT_MAX) return 0;
    dimensions[d]=(int)p->padded[d];
  }
  if(fftw_runtime_users==0 && !fftw_init_threads()) return 0;
  fftw_runtime_users++;p->fftw_runtime=1;
  fftw_plan_with_nthreads(p->threads>(size_t)INT_MAX ? INT_MAX:(int)p->threads);
  p->kernel_forward=fftw_plan_many_dft(3,dimensions,6,
      (fftw_complex *)p->data,NULL,1,(int)p->grids,
      (fftw_complex *)p->data,NULL,1,(int)p->grids,FFTW_FORWARD,FFTW_ESTIMATE);
  p->apply_forward=fftw_plan_many_dft(3,dimensions,3,
      (fftw_complex *)(p->data+6*p->grids),NULL,1,(int)p->grids,
      (fftw_complex *)(p->data+6*p->grids),NULL,1,(int)p->grids,FFTW_FORWARD,FFTW_ESTIMATE);
  p->apply_inverse=fftw_plan_many_dft(3,dimensions,3,
      (fftw_complex *)(p->data+9*p->grids),NULL,1,(int)p->grids,
      (fftw_complex *)(p->data+9*p->grids),NULL,1,(int)p->grids,FFTW_BACKWARD,FFTW_ESTIMATE);
  return p->kernel_forward && p->apply_forward && p->apply_inverse;
}
#endif
static int apply(void *context,const complex double *x,complex double *y)
{
  vie_fft_plan *p=context;
  size_t c,g;int i,j;
  complex double *src=p->data+6*p->grids,*dst=p->data+9*p->grids;
  memset(src,0,6*p->grids*sizeof(*src));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(p->cells >= 1024 && p->threads > 1) num_threads(p->threads)
#endif
  for(c=0;c<p->cells;c++) for(i=0;i<3;i++) src[i*p->grids+p->map[c]]=x[3*c+i];
  fft3_batch(p,src,3,0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(p->grids >= 16384 && p->threads > 1) num_threads(p->threads)
#endif
  for(g=0;g<3*p->grids;g++) {
    size_t component=g/p->grids,point=g%p->grids;
    for(j=0;j<3;j++) dst[g]+=p->data[tensor_index((int)component,j)*p->grids+point]*src[j*p->grids+point];
  }
  fft3_batch(p,dst,3,1);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(p->cells >= 1024 && p->threads > 1) num_threads(p->threads)
#endif
  for(c=0;c<p->cells;c++) for(i=0;i<3;i++)
    y[3*c+i]=p->epsilon[c]*x[3*c+i]-(p->epsilon[c]-1)*dst[i*p->grids+p->map[c]];
  return 1;
}
static int precondition(void *context,const complex double *x,complex double *y)
{
  vie_fft_plan *p=context;size_t c;int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(p->cells >= 1024 && p->threads > 1) num_threads(p->threads)
#endif
  for(c=0;c<p->cells;c++) for(i=0;i<3;i++) y[3*c+i]=x[3*c+i]/(p->epsilon[c]-(p->epsilon[c]-1)*p->self[i]);
  return 1;
}
void vie_fft_operator(vie_fft_plan *p,vie_linear_operator *op)
{
  memset(op,0,sizeof(*op));if(!p) return;
  op->order=3*p->cells;op->apply=apply;op->context=p;
  op->precondition=precondition;op->precondition_context=p;
}
/* Outside disjoint cubes, curl curl G = (k^2 I + grad grad) G.
   The double volume integral becomes a difference-coordinate integral with
   triangular weights, split at zero in each axis. No point-dipole sampling.
   Only far pairs use this smooth quadrature; self/near retain vie_cuboid_n. */
static void gauss(int n,double *x,double *w)
{
  int i,j;
  for(i=0;i<(n+1)/2;i++) {
    double z=cos(pi*(i+.75)/(n+.5)),old,derivative=0;
    do {
      double p=1,q=0;
      old=z;
      for(j=1;j<=n;j++) {double prev=q;q=p;p=((2*j-1)*z*q-(j-1)*prev)/j;}
      derivative=n*(z*p-q)/(z*z-1);z=old-p/derivative;
    } while(fabs(z-old)>2e-15);
    x[i]=(1-z)/2;x[n-1-i]=(1+z)/2;
    w[i]=w[n-1-i]=1/((1-z*z)*derivative*derivative);
  }
}
static void far_order(const size_t offset[3],double h,double k,int order,
                       complex double out[3][3])
{
  double x[12],w[12];int a,b,c,mask,i,j;
  gauss(order,x,w);memset(out,0,9*sizeof(complex double));
  for(a=0;a<order;a++) for(b=0;b<order;b++) for(c=0;c<order;c++) {
    double t[3]={x[a],x[b],x[c]};
    double weight=h*h*h*w[a]*w[b]*w[c]*(1-t[0])*(1-t[1])*(1-t[2]);
    for(mask=0;mask<8;mask++) {
      double rvec[3],r2=0,r;complex double g,diag,radial;
      for(i=0;i<3;i++) {rvec[i]=h*((double)offset[i]+((mask&(1<<i)) ? -t[i]:t[i]));r2+=rvec[i]*rvec[i];}
      r=sqrt(r2);g=weight*cexp(-I*k*r)/(4*pi*r);
      diag=g*(k*k-I*k/r-1/r2);
      radial=g*(-k*k+3*I*k/r+3/r2)/r2;
      for(i=0;i<3;i++) for(j=i;j<3;j++)
        out[i][j]+=(i==j ? diag:0)+radial*rvec[i]*rvec[j];
    }
  }
  for(i=0;i<3;i++) for(j=i+1;j<3;j++) out[j][i]=out[i][j];
}
static int far_kernel(const size_t v[3],double h,double k,complex double out[3][3])
{
  const int orders[]={4,6,8,12};complex double prev[3][3];int step,i,j;
  double gap2=0;
  for(i=0;i<3;i++) {double gap=v[i]>1 ? (double)v[i]-1:0;gap2+=gap*gap;}
  if(gap2<=4) return 0;
  far_order(v,h,k,orders[0],prev);
  for(step=1;step<4;step++) {
    int converged=1;far_order(v,h,k,orders[step],out);
    for(i=0;i<3;i++) for(j=0;j<3;j++)
      if(!isfinite(creal(out[i][j])) || !isfinite(cimag(out[i][j])) ||
         cabs(out[i][j]-prev[i][j])>2e-12+2e-9*cabs(out[i][j])) converged=0;
    if(converged) return 1;
    memcpy(prev,out,sizeof(prev));
  }
  return 0;
}
/* Fill signed offsets via cubic coordinate symmetries. Certified cuboid
   integration is evaluated once per sorted absolute displacement. */
static int kernels(vie_fft_plan *p,const size_t shape[3],double h,double k)
{
  static const int permutations[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
  size_t a,b,c,max=shape[0];int perm,mask,i,j;
  vie_cuboid t={{0,0,0},{h,h,h}},s=t;
  vie_kernel_options opt;
  vie_kernel_default_options(&opt);opt.relative_tolerance=5e-7;opt.absolute_tolerance=5e-10;opt.far_order=16;
  if(shape[1]>max) max=shape[1];
  if(shape[2]>max) max=shape[2];
  for(a=0;a<max;a++) for(b=0;b<=a;b++) for(c=0;c<=b;c++) {
    size_t v[3]={a,b,c};complex double n[3][3];int useful=0;
    for(perm=0;perm<6;perm++) {
      const int *q=permutations[perm];
      if(v[q[0]]<shape[0] && v[q[1]]<shape[1] && v[q[2]]<shape[2]) useful=1;
    }
    if(!useful) continue;
    for(i=0;i<3;i++) t.center[i]=h*(double)v[i];
    if(!far_kernel(v,h,k,n) && !vie_cuboid_n(&t,&s,k,&opt,n,NULL)) return 0;
    for(i=0;i<3;i++) for(j=i+1;j<3;j++) {
      complex double avg=(n[i][j]+n[j][i])*0.5;
      n[i][j]=n[j][i]=(v[i] && v[j]) ? avg:0;
    }
    if(a==0) for(i=0;i<3;i++) {
      p->self[i]=n[i][i];
    }
    for(perm=0;perm<6;perm++) {
      const int *q=permutations[perm];size_t d[3]={v[q[0]],v[q[1]],v[q[2]]};
      if(d[0]>=shape[0] || d[1]>=shape[1] || d[2]>=shape[2]) continue;
      for(mask=0;mask<8;mask++) {
        size_t pos[3],g;int sign[3];
        for(i=0;i<3;i++) {sign[i]=(mask&(1<<i)) ? -1:1;pos[i]=(sign[i]<0 && d[i]) ? p->padded[i]-d[i]:d[i];}
        g=index3(p,pos[0],pos[1],pos[2]);
        for(i=0;i<3;i++) for(j=i;j<3;j++)
          p->data[tensor_index(i,j)*p->grids+g]=sign[i]*sign[j]*n[q[i]][q[j]];
      }
    }
  }
  fft3_batch(p,p->data,6,0);
  return 1;
}
int vie_fft_create(vie_fft_plan **out,const size_t shape[3],size_t cells,
                    const size_t *xyz,double h,double k,complex double er,size_t limit)
{
  vie_fft_plan *p;size_t padded[3],g=1,bytes=sizeof(*p),tmp,root_bytes=0,c;int d;
  if(!out) return 0;
  *out=NULL;
  if(!shape || !xyz || !cells || cells>SIZE_MAX/3 || !isfinite(h) || h<=0 ||
     !isfinite(k) || k<0 || !isfinite(creal(er)) || !isfinite(cimag(er)) || creal(er)<1 || cimag(er)>0 || !isfinite(k*h) || k*h>1) return 0;
  for(d=0;d<3;d++) {
    if(!shape[d] || shape[d]>SIZE_MAX/2) return 0;
    padded[d]=1;
    while(padded[d]<2*shape[d]-1) {if(padded[d]>SIZE_MAX/2) return 0;padded[d]*=2;}
    if(!multiply(g,padded[d],&g) || !multiply(padded[d],sizeof(complex double),&tmp) || !add(root_bytes,tmp,&root_bytes)) return 0;
  }
#ifdef HAVE_FFTW
  if(g>(size_t)INT_MAX && requested_backend==FFT_BACKEND_FFTW) return 0;
#endif
  if(!multiply(g,12*sizeof(complex double),&tmp) || !add(bytes,tmp,&bytes) ||
     !multiply(cells,sizeof(size_t)+sizeof(complex double),&tmp) || !add(bytes,tmp,&bytes) || !add(bytes,root_bytes,&bytes) || bytes>limit) return 0;
  p=calloc(1,sizeof(*p));if(!p) return 0;
  p->cells=cells;p->grids=g;p->bytes=bytes;memcpy(p->padded,padded,sizeof(padded));
#ifdef HAVE_FFTW
  /* Interleaved, CPU-pinned M9 measurements put the end-to-end crossover at
     a 512-point padded grid (27 active cells in the cubic fixture).  FFTW is
     still available explicitly below that point, but its planning/runtime
     setup is counterproductive for the one- and eight-cell cases. */
  p->backend=(requested_backend==FFT_BACKEND_FFTW ||
      (requested_backend==FFT_BACKEND_AUTO && g>=(size_t)FFTW_AUTO_MIN_GRIDS &&
       g<=(size_t)INT_MAX)) ?
      FFT_BACKEND_FFTW:FFT_BACKEND_BUILTIN;
#else
  p->backend=FFT_BACKEND_BUILTIN;
#endif
  p->threads=configured_threads(p->backend);
  if(p->threads==0) p->threads=1;
  p->map=malloc(cells*sizeof(*p->map));p->data=calloc(g,12*sizeof(*p->data));
  p->epsilon=malloc(cells*sizeof(*p->epsilon));
  if(!p->map || !p->data || !p->epsilon) goto fail;
  for(c=0;c<cells;c++) p->epsilon[c]=er;
  /* Source workspace temporarily marks duplicate grid sites. */
  for(c=0;c<cells;c++) {
    size_t pos[3],idx;
    for(d=0;d<3;d++) {pos[d]=xyz[3*c+d];if(pos[d]>=shape[d]) goto fail;}
    idx=index3(p,pos[0],pos[1],pos[2]);
    if(creal(p->data[6*g+idx])!=0) goto fail;
    p->data[6*g+idx]=1;p->map[c]=idx;
  }
  memset(p->data+6*g,0,g*sizeof(*p->data));
#ifdef HAVE_FFTW
  if(p->backend==FFT_BACKEND_FFTW && !fftw_plans(p)) {
    if(requested_backend==FFT_BACKEND_FFTW) goto fail;
    fftw_release(p);p->backend=FFT_BACKEND_BUILTIN;
    p->threads=configured_threads(p->backend);
  }
#endif
  if(p->backend==FFT_BACKEND_BUILTIN) for(d=0;d<3;d++) {
    size_t j;p->roots[d]=malloc(padded[d]*sizeof(complex double));if(!p->roots[d]) goto fail;
    for(j=0;j<padded[d];j++) p->roots[d][j]=cexp(-2*pi*I*(double)j/(double)padded[d]);
  }
  if(!kernels(p,shape,h,k)) goto fail;
  for(c=0;c<cells;c++) if(!vie_fft_set_permittivity(p,c,er)) goto fail;
  last_used_backend=p->backend;*out=p;return 1;
fail:vie_fft_destroy(p);return 0;
}

/* Material multiplication is local to the testing cell; only the background
   Green tensor is convolved. No additional interface-current unknowns. */
int vie_fft_set_permittivity(vie_fft_plan *p,size_t cell,complex double er)
{
  int i;
  if(!p || cell>=p->cells || !isfinite(creal(er)) || !isfinite(cimag(er)) ||
     creal(er)<1 || cimag(er)>0) return 0;
  for(i=0;i<3;i++) {
    complex double diagonal=er-(er-1)*p->self[i];
    if(!isfinite(cabs(diagonal)) || cabs(diagonal)<1e-30) return 0;
  }
  p->epsilon[cell]=er;return 1;
}
