#include "../dielectric_fft.h"
#include "../vie_ground.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do {if(!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;}} while(0)
static int material_case(double gap,int pec)
{
  dielectric_model model={0};dielectric_mesh_result mesh={0};
  dielectric_material material[2]={{0}};dielectric_cell cells[2]={{0}};
  dielectric_diagnostic diagnostic;dielectric_fft_system *system=NULL;
  vie_linear_operator op;vie_gmres_report report;
  complex double x[6],y[6],expected[6],solution[6]={0};
  size_t i,j;int p,q;
  material[0].id=1;material[0].epsilon_r=3;material[0].epsilon_loss=.2;
  material[1].id=2;material[1].epsilon_r=2;material[1].epsilon_loss=.4;
  model.materials=material;model.material_count=2;
  mesh.cells=cells;mesh.cell_count=2;
  for(i=0;i<2;i++) {
    cells[i].material_id=i+1;
    for(j=0;j<3;j++) cells[i].size[j]=.25;
    cells[i].center[0]=.25*i;cells[i].center[1]=-.25*i;cells[i].center[2]=.125+gap+.25*i;
    for(p=0;p<3;p++) {x[3*i+p]=sin(1+3*i+p)+I*cos(.3*(3*i+p));expected[3*i+p]=dielectric_material_epsilon(&material[i],.5)*x[3*i+p];}
  }
  CHECK((pec ? dielectric_fft_prepare_pec:dielectric_fft_prepare)(&system,&model,&mesh,0,.5,0,&diagnostic));
  CHECK(dielectric_fft_storage(system)<DIELECTRIC_FFT_BUDGET);
  dielectric_fft_get_operator(system,&op);
  CHECK(op.order==6 && op.apply(op.context,x,y));
  for(i=0;i<2;i++) for(j=0;j<2;j++) {
    vie_cuboid a,b;complex double block[3][3];
    for(p=0;p<3;p++) {a.center[p]=cells[i].center[p];b.center[p]=cells[j].center[p];a.size[p]=b.size[p]=.25;}
    CHECK(pec ? vie_cuboid_pec_n(&a,&b,.5,0,NULL,block,NULL):
                vie_cuboid_n(&a,&b,.5,NULL,block,NULL));
    for(p=0;p<3;p++) for(q=0;q<3;q++) expected[3*i+p]-=(dielectric_material_epsilon(&material[i],.5)-1)*block[p][q]*x[3*j+q];
  }
  for(i=0;i<6;i++) CHECK(cabs(y[i]-expected[i])<2e-7);
  CHECK(vie_gmres(&op,y,solution,NULL,&report));
  for(i=0;i<6;i++) CHECK(cabs(solution[i]-x[i])<1e-10);
  dielectric_fft_free(system);system=NULL;
  if(!pec) return 0;
  cells[0].center[2]=.155;cells[1].center[2]=.405;
  CHECK(!dielectric_fft_prepare_pec(&system,&model,&mesh,0,.5,0,&diagnostic));
  CHECK(strstr(diagnostic.message,"half-cell lattice")!=NULL);
  return 0;
}

int main(int argc,char **argv)
{
  (void)argc;(void)argv;
  dielectric_model model;dielectric_mesh_result mesh;dielectric_material material[2];
  dielectric_cell cells[2];dielectric_diagnostic diagnostic;dielectric_fft_system *s=NULL;
  vie_linear_operator op;vie_gmres_report report;
  complex double x[8],b[8],solution[8]={0},plain[8];size_t i,j;
  memset(&model,0,sizeof(model));memset(&mesh,0,sizeof(mesh));memset(cells,0,sizeof(cells));
  memset(material,0,sizeof(material));material[0].id=1;material[0].epsilon_r=3;material[1].id=2;material[1].epsilon_r=2;
  CHECK(dielectric_material_epsilon(&material[0],0)==3);
  model.materials=material;model.material_count=2;mesh.cells=cells;mesh.cell_count=2;
  for(i=0;i<2;i++) {cells[i].material_id=1;for(j=0;j<3;j++) cells[i].size[j]=.25;}
  cells[1].center[0]=.25;
  CHECK(dielectric_fft_prepare(&s,&model,&mesh,2,.5,1024,&diagnostic));
  CHECK(dielectric_fft_storage(s)<DIELECTRIC_FFT_BUDGET);
  CHECK(dielectric_fft_entry(s,2,2)==NULL);
  for(i=0;i<8;i++) for(j=0;j<8;j++) if(i<2 || j<2)
    *dielectric_fft_entry(s,i,j)= i==j ? 4+I:.002*(1+i+2*j)+.003*I*(2+(double)i-(double)j);
  dielectric_fft_get_operator(s,&op);
  for(i=0;i<8;i++) x[i]=sin(.5*i)+I*cos(.3*i);
  CHECK(op.apply(op.context,x,b));
  CHECK(vie_gmres(&op,b,solution,NULL,&report));
  for(i=0;i<8;i++) CHECK(cabs(x[i]-solution[i])<1e-10);
  CHECK(op.precondition(op.precondition_context,b,plain));
  for(i=0;i<2;i++) CHECK(cabs(plain[i]-b[i]/(4+I))<1e-14);
  dielectric_fft_free(s);s=NULL;
  CHECK(!dielectric_fft_prepare(&s,&model,&mesh,2,.5,DIELECTRIC_FFT_BUDGET,&diagnostic));
  cells[1].material_id=2;CHECK(dielectric_fft_prepare(&s,&model,&mesh,0,.5,0,&diagnostic));
  dielectric_fft_free(s);s=NULL;cells[1].material_id=1;
  cells[1].size[2]=.125;CHECK(!dielectric_fft_prepare(&s,&model,&mesh,0,.5,0,&diagnostic));cells[1].size[2]=.25;
  cells[1].center[0]=.3;CHECK(!dielectric_fft_prepare(&s,&model,&mesh,0,.5,0,&diagnostic));
  cells[1].center[0]=0;CHECK(!dielectric_fft_prepare(&s,&model,&mesh,0,.5,0,&diagnostic));
  cells[1].center[0]=1e6;CHECK(!dielectric_fft_prepare(&s,&model,&mesh,0,.5,0,&diagnostic));
  cells[1].center[0]=NAN;CHECK(!dielectric_fft_prepare(&s,&model,&mesh,0,.5,0,&diagnostic));
  mesh.cell_count=0;CHECK(dielectric_fft_prepare(&s,&model,&mesh,2,.5,0,&diagnostic));
  *dielectric_fft_entry(s,0,0)=2;*dielectric_fft_entry(s,1,1)=3;
  dielectric_fft_get_operator(s,&op);CHECK(op.apply(op.context,x,b));
  CHECK(cabs(b[0]-2*x[0])+cabs(b[1]-3*x[1])<1e-14);
  dielectric_fft_free(s);
  CHECK(material_case(0,0)==0);CHECK(material_case(0,1)==0);CHECK(material_case(.125,1)==0);
  puts("Hybrid free-space/PEC-image FFT operator, dense agreement, lattice, and budget checks passed");return 0;
}
