#include "dielectric_fft.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

vie_gmres_options dielectric_gmres_settings={30,300,1e-12,0};
struct dielectric_fft_system {
  size_t c,v,bytes;
  complex double *blocks, *image_work;
  vie_fft_plan *plan;
  vie_linear_operator volume;
};
static int add_product(size_t *sum,size_t a,size_t b)
{
  size_t n;
  if(b && a>SIZE_MAX/b) return 0;
  n=a*b;if(*sum>SIZE_MAX-n) return 0;*sum+=n;return 1;
}
static int error(dielectric_diagnostic *d,const char *message)
{
  memset(d,0,sizeof(*d));
  snprintf(d->message,sizeof(d->message),"FFT: %s",message);
  return 0;
}
void dielectric_fft_free(dielectric_fft_system *s)
{if(s) {vie_fft_destroy(s->plan);free(s->blocks);free(s->image_work);free(s);}}
size_t dielectric_fft_storage(const dielectric_fft_system *s) {return s ? s->bytes:0;}
complex double *dielectric_fft_entry(dielectric_fft_system *s,size_t row,size_t column)
{
  size_t c=s->c,v=s->v;
  if(row<c) return s->blocks+row*(c+v)+column;
  if(column<c) return s->blocks+c*(c+v)+(row-c)*c+column;
  return NULL;
}
static int apply(void *context,const complex double *x,complex double *y)
{
  dielectric_fft_system *s=context;size_t i,j,c=s->c,v=s->v;
  if(v && s->image_work) {
    complex double *input=s->image_work,*output=input+2*v;
    for(i=0;i<v;i++) {input[i]=x[c+i];input[v+i]=(i%3==2?1:-1)*x[c+i];}
    if(!s->volume.apply(s->volume.context,input,output)) return 0;
    memcpy(y+c,output,v*sizeof(*y));
  } else if(v && !s->volume.apply(s->volume.context,x+c,y+c)) return 0;
  for(i=0;i<c;i++) {
    y[i]=0;
    for(j=0;j<c+v;j++) y[i]+=s->blocks[i*(c+v)+j]*x[j];
  }
  for(i=0;i<v;i++) for(j=0;j<c;j++)
    y[c+i]+=s->blocks[c*(c+v)+i*c+j]*x[j];
  return 1;
}
static int precondition(void *context,const complex double *x,complex double *y)
{
  dielectric_fft_system *s=context;size_t i;
  if(s->v && s->image_work) {
    complex double *input=s->image_work,*output=input+2*s->v;
    memcpy(input,x+s->c,s->v*sizeof(*input));
    memset(input+s->v,0,s->v*sizeof(*input));
    if(!s->volume.precondition(s->volume.precondition_context,input,output)) return 0;
    memcpy(y+s->c,output,s->v*sizeof(*y));
  } else if(s->v && !s->volume.precondition(s->volume.precondition_context,x+s->c,y+s->c)) return 0;
  for(i=0;i<s->c;i++) {
    complex double diagonal=*dielectric_fft_entry(s,i,i);
    if(!isfinite(cabs(diagonal)) || cabs(diagonal)<DBL_MIN) return 0;
    y[i]=x[i]/diagonal;
  }
  return 1;
}
void dielectric_fft_get_operator(dielectric_fft_system *s,vie_linear_operator *op)
{
  op->order=s->c+s->v;op->apply=apply;op->context=s;
  op->precondition=precondition;op->precondition_context=s;
}
static int prepare(dielectric_fft_system **out,const dielectric_model *model,
    const dielectric_mesh_result *mesh,size_t conductors,double k0,size_t extra,
    dielectric_diagnostic *diagnostic,int pec)
{
  dielectric_fft_system *s=NULL;
  size_t cells=mesh->cell_count,n,v,bytes=sizeof(*s),blocks=0,*xyz=NULL;
  size_t shape[3]={1,1,1},i,m=dielectric_gmres_settings.restart;
  double origin[3]={DBL_MAX,DBL_MAX,DBL_MAX},h=0; complex double er=0;int d;
  *out=NULL;
  if(cells>SIZE_MAX/3 || conductors>SIZE_MAX-3*cells) return error(diagnostic,"unknown count overflow");
  v=3*cells;n=v+conductors;
  if(pec && (v>SIZE_MAX/4 || !add_product(&bytes,4*v,sizeof(complex double)) ||
      !add_product(&bytes,v,sizeof(size_t)))) return error(diagnostic,"PEC image storage overflow");
  if(!add_product(&blocks,conductors,conductors) || !add_product(&blocks,conductors,v) ||
     !add_product(&blocks,conductors,v) || !add_product(&bytes,blocks ? blocks:1,sizeof(complex double)) ||
     !add_product(&bytes,n,(2*m+3)*sizeof(complex double)) ||
     !add_product(&bytes,(m+1)*m+3*m+1,sizeof(complex double)) ||
     !add_product(&bytes,m,sizeof(double)) || !add_product(&bytes,v,sizeof(size_t)) ||
     !add_product(&bytes,cells,sizeof(dielectric_cell)) || !add_product(&bytes,1,extra) || bytes>=DIELECTRIC_FFT_BUDGET)
    return error(diagnostic,"solver storage exceeds 512 MiB budget");
  if(cells) {
    h=mesh->cells[0].size[0];
    for(i=0;i<cells;i++) {
      size_t material;complex double value=0;
      if(pec && mesh->cells[i].center[2]-.5*mesh->cells[i].size[2]<0)
        return error(diagnostic,"PEC FFT cells must lie in z >= 0");
      for(material=0;material<model->material_count;material++)
        if(model->materials[material].id==mesh->cells[i].material_id) value=dielectric_material_epsilon(&model->materials[material],k0);
      if(!value || !isfinite(creal(value)) || !isfinite(cimag(value)) || creal(value)<1 || cimag(value)>0)
        return error(diagnostic,"unknown or invalid passive material");
      if(i==0) er=value;
      for(d=0;d<3;d++) {
        if(!isfinite(h) || h<=0 || !isfinite(mesh->cells[i].size[d]) ||
           !isfinite(mesh->cells[i].center[d]) || fabs(mesh->cells[i].size[d]-h)>1e-10*h)
          return error(diagnostic,"requires equal cubic cells");
        if(mesh->cells[i].center[d]<origin[d]) origin[d]=mesh->cells[i].center[d];
      }
    }
    xyz=malloc((pec?2:1)*v*sizeof(*xyz));if(!xyz) return error(diagnostic,"lattice map allocation failed");
    for(i=0;i<cells;i++) for(d=0;d<3;d++) {
      double index=(mesh->cells[i].center[d]-origin[d])/h,nearest=round(index);
      if(!isfinite(index) || nearest<0 || nearest>=(double)(SIZE_MAX/2) || fabs(index-nearest)>1e-8) {
        free(xyz);return error(diagnostic,"cell centers do not share a common lattice");
      }
      xyz[3*i+d]=(size_t)nearest;
      if(xyz[3*i+d]>=shape[d]) shape[d]=xyz[3*i+d]+1;
    }
  }
  if(pec && cells) {
    double gap_index=2*(origin[2]-.5*h)/h,nearest=round(gap_index);
    size_t nz=shape[2],offset;
    if(!isfinite(gap_index) || nearest<0 || nearest>=(double)(SIZE_MAX/4) ||
       fabs(gap_index-nearest)>1e-8) {
      free(xyz);return error(diagnostic,"PEC FFT requires nonnegative bottom gap on the half-cell lattice");
    }
    offset=(size_t)nearest;
    if(nz>SIZE_MAX/2 || offset>SIZE_MAX-2*nz) {
      free(xyz);return error(diagnostic,"PEC image embedding overflow");
    }
    shape[2]=2*nz+offset;
    for(i=0;i<cells;i++) {
      xyz[v+3*i]=xyz[3*i];xyz[v+3*i+1]=xyz[3*i+1];
      xyz[v+3*i+2]=nz-1-xyz[3*i+2];
      xyz[3*i+2]+=nz+offset;
    }
  }
  s=calloc(1,sizeof(*s));if(!s) goto failed;
  s->c=conductors;s->v=v;
  s->blocks=calloc(blocks ? blocks:1,sizeof(*s->blocks));if(!s->blocks) goto failed;
  if(pec && v) {s->image_work=calloc(4*v,sizeof(*s->image_work));if(!s->image_work) goto failed;}
  if(cells && !vie_fft_create(&s->plan,shape,(pec?2:1)*cells,xyz,h,k0,er,DIELECTRIC_FFT_BUDGET-bytes)) goto failed;
  for(i=0;i<cells;i++) {
    size_t material;
    for(material=0;material<model->material_count;material++)
      if(model->materials[material].id==mesh->cells[i].material_id) {
        complex double value=dielectric_material_epsilon(&model->materials[material],k0);
        if(!vie_fft_set_permittivity(s->plan,i,value) ||
           (pec && !vie_fft_set_permittivity(s->plan,cells+i,value))) goto failed;
        break;
      }
  }
  if(cells) vie_fft_operator(s->plan,&s->volume);
  s->bytes=bytes+vie_fft_bytes(s->plan);free(xyz);*out=s;return 1;
failed:
  free(xyz);dielectric_fft_free(s);
  return error(diagnostic,"kernel/lattice build failed or padded FFT exceeds memory budget");
}

int dielectric_fft_prepare(dielectric_fft_system **out,const dielectric_model *model,
    const dielectric_mesh_result *mesh,size_t conductors,double k0,size_t extra,
    dielectric_diagnostic *diagnostic)
{return prepare(out,model,mesh,conductors,k0,extra,diagnostic,0);}

int dielectric_fft_prepare_pec(dielectric_fft_system **out,const dielectric_model *model,
    const dielectric_mesh_result *mesh,size_t conductors,double k0,size_t extra,
    dielectric_diagnostic *diagnostic)
{return prepare(out,model,mesh,conductors,k0,extra,diagnostic,1);}
