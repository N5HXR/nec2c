#include "vie_ground.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
static int failures;
#define REQUIRE(c, m) do { if(!(c)) { fprintf(stderr,"FAIL: %s\n",m); ++failures; } } while(0)

static double volume(const vie_cuboid *c)
{
  return c->size[0]*c->size[1]*c->size[2];
}

/* Independent, nonsingular six-dimensional Gauss integration of the
 * closed-form electric dipole dyad. Deliberately does not call vie_kernel
 * or use its surface-reduction/difference-coordinate quadrature. */
static void oracle(const vie_cuboid *a, const vie_cuboid *b, double k,
                    complex double out[3][3])
{
  const double nodes[3] = {-0.774596669241483377,0,0.774596669241483377};
  const double weights[3] = {5.0/9.0,8.0/9.0,5.0/9.0};
  int i,j,p,q,axis;
  memset(out,0,9*sizeof(complex double));
  for(i=0;i<27;++i) for(j=0;j<27;++j)
  {
    double x[3], y[3], weight=volume(b)/64.0;
    int ti=i,sj=j;
    for(axis=0;axis<3;++axis)
    {
      int u=ti%3,v=sj%3;
      ti/=3; sj/=3;
      x[axis]=a->center[axis]+0.5*a->size[axis]*nodes[u];
      y[axis]=b->center[axis]+0.5*b->size[axis]*nodes[v];
      weight*=weights[u]*weights[v];
    }
    for(axis=0;axis<2;++axis)
    {
      double d[3]={x[0]-y[0],x[1]-y[1],x[2]-(axis ? -y[2]:y[2])};
      double r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
      double kr=k*r;
      complex double f=weight*cexp(-I*kr)/(4*PI*r*r*r);
      for(p=0;p<3;++p) for(q=0;q<3;++q)
      {
        double sign=axis && q!=2 ? -1:1;
        out[p][q]+=sign*f*((p==q ? kr*kr-1-I*kr:0)+
                      (3+3*I*kr-kr*kr)*d[p]*d[q]/(r*r));
      }
    }
  }
}

static void test_oracle_and_reciprocity(void)
{
  vie_cuboid a={{.1,-.2,.35},{.04,.05,.06}};
  vie_cuboid b={{.8,.3,.55},{.05,.03,.04}};
  complex double ab[3][3],ba[3][3],ref[3][3];
  int f,p,q;
  for(f=0;f<2;++f)
  {
    double k=f ? 2.3:0;
    REQUIRE(vie_cuboid_pec_n(&a,&b,k,.1,NULL,ab,NULL),"PEC separated pair");
    REQUIRE(vie_cuboid_pec_n(&b,&a,k,.1,NULL,ba,NULL),"PEC reverse pair");
    oracle(&a,&b,k,ref);
    for(p=0;p<3;++p) for(q=0;q<3;++q)
    {
      REQUIRE(cabs(ab[p][q]-ref[p][q])<2e-12,"independent static/dynamic dipole integration");
      REQUIRE(cabs(volume(&a)*ab[p][q]-volume(&b)*ba[q][p])<1e-16,
              "volume-weighted PEC reciprocity");
    }
  }
}

static void test_boundary(void)
{
  vie_cuboid a={{.3,-.15,1e-6},{.02,.02,1e-6}};
  vie_cuboid b={{0,0,.4},{.04,.04,.04}};
  complex double low[3][3],high[3][3];
  double tangential_low=0,tangential_high=0,normal=0;
  int p,q;
  REQUIRE(vie_cuboid_pec_n(&a,&b,1.3,0,NULL,low,NULL),"near-plane evaluation");
  a.center[2]=1e-5;
  REQUIRE(vie_cuboid_pec_n(&a,&b,1.3,0,NULL,high,NULL),"raised evaluation");
  for(p=0;p<2;++p) for(q=0;q<3;++q)
  {
    tangential_low+=cabs(low[p][q]); tangential_high+=cabs(high[p][q]);
  }
  for(q=0;q<3;++q) normal+=cabs(low[2][q]);
  REQUIRE(tangential_low<.101*tangential_high,"tangential E vanishes linearly at PEC");
  REQUIRE(tangential_low<1e-4*normal,"PEC preserves normal electric field");
}

static void test_self_and_far_ground(void)
{
  vie_cuboid a={{0,0,.2},{.04,.04,.04}},image=a,far=a;
  complex double pec[3][3],direct[3][3],reflected[3][3],remote[3][3];
  vie_pec_report report;
  int p,q;
  image.center[2]=-.2;
  REQUIRE(vie_cuboid_pec_n(&a,&a,0,.1,NULL,pec,&report),"PEC self block");
  REQUIRE(report.direct.relation==VIE_RELATION_SELF,"direct singular self retained");
  REQUIRE(vie_cuboid_n(&a,&a,0,NULL,direct,NULL),"free self");
  REQUIRE(vie_cuboid_n(&a,&image,0,NULL,reflected,NULL),"image self");
  far.center[2]=100;
  REQUIRE(vie_cuboid_pec_n(&far,&far,0,.1,NULL,remote,NULL),"remote ground self");
  for(p=0;p<3;++p) for(q=0;q<3;++q)
  {
    REQUIRE(cabs(pec[p][q]-direct[p][q]-(q==2?1:-1)*reflected[p][q])<1e-15,
            "self includes finite image correction");
    REQUIRE(cabs(remote[p][q]-direct[p][q])<2e-12,"remote ground approaches free space");
  }
}

static void test_rejection(void)
{
  vie_cuboid a={{0,0,.2},{.04,.04,.04}},bad=a;
  complex double out[3][3],saved[3][3];
  int i;
  for(i=0;i<9;++i) ((complex double *)out)[i]=7+2*I;
  memcpy(saved,out,sizeof(out));
  bad.center[2]=.019;
  REQUIRE(!vie_cuboid_pec_n(&a,&bad,1,0,NULL,out,NULL),"reject straddling source cell");
  bad.center[2]=-.1;
  REQUIRE(!vie_cuboid_pec_n(&bad,&a,1,0,NULL,out,NULL),"reject buried testing cell");
  REQUIRE(!vie_cuboid_pec_n(&a,&a,1,.19,NULL,out,NULL),"enforce clearance");
  REQUIRE(!vie_cuboid_pec_n(&a,&a,1,NAN,NULL,out,NULL),"reject NaN gap");
  REQUIRE(!vie_cuboid_pec_n(&a,&a,1,-1,NULL,out,NULL),"reject negative gap");
  REQUIRE(!vie_cuboid_pec_n(&a,&a,NAN,0,NULL,out,NULL),"reject NaN frequency");
  REQUIRE(!vie_cuboid_pec_n(NULL,&a,1,0,NULL,out,NULL),"reject null cell");
  REQUIRE(memcmp(out,saved,sizeof(out))==0,"failed calls leave result untouched");
}

static void test_pec_contact(void)
{
  /* Independent M4 high-order face oracle, rotated x->z, y->x, z->y.
     The upper cell and its image share a face, not a volume. */
  const complex double face[3]={
    -0.073160835986049028-I*0.010972943890459171,
    -0.052006267652284888-I*0.011030863498482452,
     0.18171808983515619-I*0.011131794107809491};
  vie_cuboid a={{0,0,.1},{.3,.4,.2}},raised=a,b={{.31,.07,.06},{.12,.16,.12}};
  complex double contact[3][3],direct[3][3],near[3][3],ab[3][3],ba[3][3];
  vie_pec_report report;
  const double gaps[6]={.01,.003,.001,.0001,.00001,.000001};
  double errors[6]={0,0,0,0,0,0};
  int p,q,i;
  REQUIRE(vie_cuboid_pec_n(&a,&a,2.1,0,NULL,contact,&report),"PEC face contact self");
  REQUIRE(report.direct.relation==VIE_RELATION_SELF,"contact direct self classification");
  REQUIRE(report.reflected.relation==VIE_RELATION_FACE,"contact reflected face classification");
  REQUIRE(vie_cuboid_n(&a,&a,2.1,NULL,direct,NULL),"contact direct reference");
  for(p=0;p<3;++p) for(q=0;q<3;++q)
    REQUIRE(cabs(contact[p][q]-direct[p][q]-(q==2?1:-1)*(p==q?face[p]:0))<4e-10,
            "PEC contact image agrees with independent high-order face oracle");
  for(i=0;i<6;++i)
  {
    raised.center[2]=.1+gaps[i];
    if(!vie_cuboid_pec_n(&raised,&raised,2.1,0,NULL,near,NULL))
    {
      fprintf(stderr,"gap %.9g failed\n",gaps[i]);
      REQUIRE(0,"positive-gap contact approach");
      return;
    }
    for(p=0;p<3;++p) for(q=0;q<3;++q)
      errors[i]=fmax(errors[i],cabs(near[p][q]-contact[p][q]));
  }
  for(i=1;i<6;++i)
    REQUIRE(errors[i]<.4*errors[i-1],"positive gaps approach PEC contact through micron scale");
  REQUIRE(vie_cuboid_pec_n(&a,&b,2.1,0,NULL,ab,NULL),"unequal PEC contact pair");
  REQUIRE(vie_cuboid_pec_n(&b,&a,2.1,0,NULL,ba,NULL),"reverse PEC contact pair");
  for(p=0;p<3;++p) for(q=0;q<3;++q)
    REQUIRE(cabs(volume(&a)*ab[p][q]-volume(&b)*ba[q][p])<2e-12,"contact weighted reciprocity");
  REQUIRE(!vie_cuboid_pec_n(&a,&a,2.1,1e-6,NULL,near,NULL),"positive minimum gap still excludes contact");
}

static void test_static_contact_oracle(void)
{
  const char *root=getenv("TEST_SRCDIR");
  char path[4096],line[1024];
  FILE *stream;
  int count=0;
  snprintf(path,sizeof(path),"%s/tests/m11-static-reference.csv",root?root:".");
  stream=fopen(path,"r");
  REQUIRE(stream!=NULL,"open independent closed-form PEC reference");
  if(!stream) return;
  while(fgets(line,sizeof(line),stream))
  {
    double sx,sy,sz,gap,expected[3];
    complex double actual[3][3];
    int p,q;
    vie_cuboid cell;
    if(line[0]=='#') continue;
    if(sscanf(line,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&sx,&sy,&sz,&gap,
              &expected[0],&expected[1],&expected[2])!=7)
    { REQUIRE(0,"parse independent PEC reference"); continue; }
    cell=(vie_cuboid){{0,0,sz/2+gap},{sx,sy,sz}};
    if(!vie_cuboid_pec_n(&cell,&cell,0,0,NULL,actual,NULL))
    { fprintf(stderr,"static gap %.9g failed\n",gap); REQUIRE(0,"static PEC reference convergence");continue; }
    for(p=0;p<3;++p) for(q=0;q<3;++q)
      REQUIRE(cabs(actual[p][q]-(p==q?expected[p]:0))<5e-9,
              "PEC contact/near-contact matches independent analytic Coulomb integrals");
    ++count;
  }
  fclose(stream);
  REQUIRE(count==8,"all independent static contact references checked");
}

int main(void)
{
  test_oracle_and_reciprocity(); test_boundary();
  test_self_and_far_ground(); test_rejection(); test_pec_contact(); test_static_contact_oracle();
  if(failures) return 1;
  puts("PEC volume kernel: independent dipole/contact quadrature, resolved-gap continuity, boundary, reciprocity, self, remote-ground and rejection checks passed");
  return 0;
}
