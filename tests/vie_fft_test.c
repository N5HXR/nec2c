#include "../vie_fft.h"
#include "../vie_kernel.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#define CHECK(x) do {if(!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;}} while(0)
static int fixture(const size_t shape[3],const size_t *xyz,size_t cells,double k,complex double er)
{
  vie_fft_plan *p=NULL;vie_linear_operator op;vie_gmres_options opts;vie_gmres_report report;
  complex double *a,*x,*y,*b,*solution;size_t n=3*cells,i,j,c,d;int u,v;size_t repeat;
  const double h=0.25;
  double maximum_error=0;
  a=calloc(n*n,sizeof(*a));x=calloc(n,sizeof(*x));y=calloc(n,sizeof(*y));b=calloc(n,sizeof(*b));solution=calloc(n,sizeof(*solution));
  CHECK(a && x && y && b && solution);
  CHECK(vie_fft_create(&p,shape,cells,xyz,h,k,er,32*1024*1024));
  vie_fft_operator(p,&op);CHECK(op.order==n);
  for(c=0;c<cells;c++) for(d=0;d<cells;d++) {
    vie_cuboid t={{0,0,0},{h,h,h}},s=t;complex double block[3][3];vie_kernel_options ko;
    vie_kernel_default_options(&ko);ko.relative_tolerance=5e-7;ko.absolute_tolerance=5e-10;ko.far_order=16;
    for(u=0;u<3;u++) {t.center[u]=h*xyz[3*c+u];s.center[u]=h*xyz[3*d+u];}
    CHECK(vie_cuboid_n(&t,&s,k,&ko,block,NULL));
    for(u=0;u<3;u++) for(v=0;v<3;v++) a[(3*c+u)*n+3*d+v]=-(er-1)*block[u][v]+((c==d && u==v) ? er:0);
  }
  for(repeat=0;repeat<n+2;repeat++) {
    for(i=0;i<n;i++) x[i]=(repeat>0 && repeat<=n) ? (i==repeat-1 ? 1+2*I:0):sin(0.7*i)+I*cos(0.3*i);
    for(i=0;i<n;i++) {b[i]=0;for(j=0;j<n;j++) b[i]+=a[i*n+j]*x[j];}
    CHECK(op.apply(op.context,x,y));
    for(i=0;i<n;i++) {
      double error=cabs(y[i]-b[i])/fmax(1,cabs(b[i]));
      if(error>maximum_error) maximum_error=error;
      CHECK(error<2e-8);
    }
  }
  vie_gmres_default_options(&opts);opts.relative_tolerance=1e-11;opts.maximum_iterations=150;
  CHECK(vie_gmres(&op,b,solution,&opts,&report));
  CHECK(report.relative_residual<=opts.relative_tolerance);
  for(i=0;i<n;i++) CHECK(cabs(solution[i]-x[i])<2e-8);
  printf("%zu cells k=%g er=(%g,%g): bytes=%zu iterations=%zu residual=%.3g product_error=%.3g\n",cells,k,creal(er),cimag(er),vie_fft_bytes(p),report.iterations,report.relative_residual,maximum_error);
  vie_fft_destroy(p);free(a);free(x);free(y);free(b);free(solution);return 0;
}
int main(void)
{
  const size_t shape[3]={5,3,2},xyz[]={0,0,0,4,2,1,1,0,0,0,1,0,0,0,1,1,1,1};
  const size_t far_shape[3]={9,5,3},far_xyz[]={0,0,0,8,4,2,4,1,0};
  const size_t one[3]={1,1,1},zero[3]={0,0,0},duplicate[]={0,0,0,0,0,0};
  const size_t overflow[3]={SIZE_MAX,2,2};vie_fft_plan *p=NULL;
  CHECK(!vie_fft_create(&p,shape,6,xyz,.25,.2,3,1));CHECK(p==NULL);
  CHECK(!vie_fft_create(&p,shape,2,duplicate,.25,.2,3,1000000));CHECK(p==NULL);
  CHECK(!vie_fft_create(&p,one,6,xyz,.25,.2,3,1000000));
  CHECK(!vie_fft_create(&p,overflow,1,zero,.25,.2,3,SIZE_MAX));
  CHECK(!vie_fft_create(&p,one,1,zero,.25,5,3,1000000));
  CHECK(!fixture(one,zero,1,.2,3));
  CHECK(!fixture(shape,xyz,6,0,3));CHECK(!fixture(shape,xyz,6,.8,3));
  CHECK(!fixture(far_shape,far_xyz,3,4,3));
  CHECK(!fixture(one,zero,1,.2,1));
  CHECK(!fixture(one,zero,1,.2,3-.4*I));
  CHECK(!fixture(shape,xyz,6,.8,3-.4*I));
  CHECK(!vie_fft_create(&p,one,1,zero,.25,.2,3+.4*I,1000000));
  CHECK(vie_fft_create(&p,one,1,zero,.25,.2,3,1000000));
  CHECK(!vie_fft_set_permittivity(p,1,2));
  CHECK(!vie_fft_set_permittivity(p,0,3+.1*I));
  CHECK(!vie_fft_set_permittivity(p,0,NAN));
  CHECK(vie_fft_set_permittivity(p,0,1-.4*I));
  vie_fft_destroy(p);
  puts("FFT operator tests passed");return 0;
}
