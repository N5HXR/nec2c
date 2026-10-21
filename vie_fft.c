#include "vie_fft.h"
#include "vie_kernel.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct vie_fft_plan {
  size_t cells, padded[3], grids, bytes;
  size_t *map;
  complex double *data, *roots[3], self[3], *epsilon;
};
static const double pi = 3.141592653589793238462643383279502884;
static int multiply(size_t a, size_t b, size_t *out)
{ if(b && a > SIZE_MAX/b) return 0; *out=a*b; return 1; }
static int add(size_t a, size_t b, size_t *out)
{ if(a > SIZE_MAX-b) return 0; *out=a+b; return 1; }
static size_t index3(const vie_fft_plan *p,size_t x,size_t y,size_t z)
{ return (x*p->padded[1]+y)*p->padded[2]+z; }
static int tensor_index(int i,int j)
{ if(i==j) return i; if(i>j) {int t=i;i=j;j=t;} return i==0 ? (j==1 ? 3:4):5; }
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
static void fft3(vie_fft_plan *p,complex double *a,int inverse)
{
  size_t x,y,z,nx=p->padded[0],ny=p->padded[1],nz=p->padded[2];
  for(x=0;x<nx;x++) for(y=0;y<ny;y++)
    fft_line(a+index3(p,x,y,0),nz,1,p->roots[2],inverse);
  for(x=0;x<nx;x++) for(z=0;z<nz;z++)
    fft_line(a+index3(p,x,0,z),ny,nz,p->roots[1],inverse);
  for(y=0;y<ny;y++) for(z=0;z<nz;z++)
    fft_line(a+index3(p,0,y,z),nx,ny*nz,p->roots[0],inverse);
}
void vie_fft_destroy(vie_fft_plan *p)
{ if(p) {int d;for(d=0;d<3;d++) free(p->roots[d]);free(p->map);free(p->data);free(p->epsilon);free(p);} }
size_t vie_fft_bytes(const vie_fft_plan *p) {return p ? p->bytes:0;}
static int apply(void *context,const complex double *x,complex double *y)
{
  vie_fft_plan *p=context;
  size_t c,g;int i,j;
  complex double *src=p->data+6*p->grids,*dst=p->data+9*p->grids;
  memset(src,0,6*p->grids*sizeof(*src));
  for(c=0;c<p->cells;c++) for(i=0;i<3;i++) src[i*p->grids+p->map[c]]=x[3*c+i];
  for(i=0;i<3;i++) fft3(p,src+i*p->grids,0);
  for(i=0;i<3;i++) {
    for(g=0;g<p->grids;g++) for(j=0;j<3;j++)
      dst[i*p->grids+g]+=p->data[tensor_index(i,j)*p->grids+g]*src[j*p->grids+g];
    fft3(p,dst+i*p->grids,1);
  }
  for(c=0;c<p->cells;c++) for(i=0;i<3;i++)
    y[3*c+i]=p->epsilon[c]*x[3*c+i]-(p->epsilon[c]-1)*dst[i*p->grids+p->map[c]];
  return 1;
}
static int precondition(void *context,const complex double *x,complex double *y)
{
  vie_fft_plan *p=context;size_t c;int i;
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
  for(i=0;i<6;i++) fft3(p,p->data+i*p->grids,0);
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
  if(!multiply(g,12*sizeof(complex double),&tmp) || !add(bytes,tmp,&bytes) ||
     !multiply(cells,sizeof(size_t)+sizeof(complex double),&tmp) || !add(bytes,tmp,&bytes) || !add(bytes,root_bytes,&bytes) || bytes>limit) return 0;
  p=calloc(1,sizeof(*p));if(!p) return 0;
  p->cells=cells;p->grids=g;p->bytes=bytes;memcpy(p->padded,padded,sizeof(padded));
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
  for(d=0;d<3;d++) {
    size_t j;p->roots[d]=malloc(padded[d]*sizeof(complex double));if(!p->roots[d]) goto fail;
    for(j=0;j<padded[d];j++) p->roots[d][j]=cexp(-2*pi*I*(double)j/(double)padded[d]);
  }
  if(!kernels(p,shape,h,k)) goto fail;
  for(c=0;c<cells;c++) if(!vie_fft_set_permittivity(p,c,er)) goto fail;
  *out=p;return 1;
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
