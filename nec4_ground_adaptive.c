#include "nec4_ground_adaptive.h"

#include <float.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

struct nec4_ground_adaptive_node {
  nec4_ground_adaptive_bounds bounds;
  double complex *values;
  struct nec4_ground_adaptive_node *children[8];
  double withheld_error;
  unsigned int depth;
  unsigned char active_mask,child_count,accepted;
};

typedef struct {
  nec4_ground_adaptive *table;
  nec4_ground_adaptive_sampler sampler;
  void *sampler_context;
} build_context;

static unsigned int axis_count(unsigned char mask,unsigned int axis)
{ return mask&(1u<<axis)?4u:1u; }

static unsigned char active_mask(const nec4_ground_adaptive_bounds *bounds)
{
  unsigned char mask=0;
  for(unsigned int axis=0;axis<3;++axis)
    if(bounds->upper[axis]>bounds->lower[axis]) mask|=(unsigned char)(1u<<axis);
  return mask;
}

static size_t sample_count(unsigned char mask,unsigned int per_axis)
{
  size_t count=1;
  for(unsigned int axis=0;axis<3;++axis)
    if(mask&(1u<<axis)) count*=per_axis;
  return count;
}

static void decode(size_t index,unsigned char mask,unsigned int radix,
    unsigned int coordinate[3])
{
  for(int axis=2;axis>=0;--axis) {
    if(mask&(1u<<axis)) { coordinate[axis]=(unsigned int)(index%radix);index/=radix; }
    else coordinate[axis]=0;
  }
}

static double coordinate(const nec4_ground_adaptive_bounds *bounds,
    unsigned int axis,unsigned int index,unsigned int denominator)
{
  double span=bounds->upper[axis]-bounds->lower[axis];
  return bounds->lower[axis]+span*index/denominator;
}

static void cubic_weights(double x,double weight[4])
{
  /* Lagrange basis for nodes 0, 1/3, 2/3, 1. */
  weight[0]=-.5*(3.*x-1.)*(3.*x-2.)*(x-1.);
  weight[1]= 4.5*x*(3.*x-2.)*(x-1.);
  weight[2]=-4.5*x*(3.*x-1.)*(x-1.);
  weight[3]= .5*x*(3.*x-1.)*(3.*x-2.);
}

static int interpolate_node(const nec4_ground_adaptive_node *node,
    double rho,double z,double depth,double complex values[5])
{
  const double point[3]={rho,z,depth};
  double weights[3][4]={{0}};
  unsigned int count[3],i,j,k;
  if(!node||!node->values||!values) return 0;
  for(unsigned int axis=0;axis<3;++axis) {
    double span=node->bounds.upper[axis]-node->bounds.lower[axis];
    double tolerance=64.*DBL_EPSILON*fmax(1.,
      fmax(fabs(node->bounds.lower[axis]),fabs(node->bounds.upper[axis])));
    if(point[axis]<node->bounds.lower[axis]-tolerance||
       point[axis]>node->bounds.upper[axis]+tolerance) return 0;
    count[axis]=axis_count(node->active_mask,axis);
    if(count[axis]==1) weights[axis][0]=1.;
    else {
      double fraction=(point[axis]-node->bounds.lower[axis])/span;
      fraction=fmax(0.,fmin(1.,fraction));cubic_weights(fraction,weights[axis]);
    }
  }
  for(unsigned int component=0;component<5;++component) values[component]=0.;
  for(i=0;i<count[0];++i) for(j=0;j<count[1];++j) for(k=0;k<count[2];++k) {
    size_t index=(i*count[1]+j)*count[2]+k;
    double weight=weights[0][i]*weights[1][j]*weights[2][k];
    for(unsigned int component=0;component<5;++component)
      values[component]+=weight*node->values[index*5+component];
  }
  return 1;
}

static double scaled_error(const double complex actual[5],
    const double complex expected[5],double cell_scale)
{
  double difference=0.,point_scale=0.;
  for(unsigned int component=0;component<5;++component) {
    difference=fmax(difference,cabs(actual[component]-expected[component]));
    point_scale=fmax(point_scale,cabs(expected[component]));
  }
  return difference/fmax(point_scale,1.e-12*fmax(cell_scale,1.e-300));
}

static void destroy_node(nec4_ground_adaptive_node *node)
{
  if(!node) return;
  for(unsigned int child=0;child<node->child_count;++child)
    destroy_node(node->children[child]);
  free(node->values);free(node);
}

void nec4_ground_adaptive_init(nec4_ground_adaptive *table)
{ if(table) memset(table,0,sizeof(*table)); }

void nec4_ground_adaptive_destroy(nec4_ground_adaptive *table)
{
  if(!table) return;
  destroy_node(table->root);memset(table,0,sizeof(*table));
}

static int sample_fit(build_context *context,nec4_ground_adaptive_node *node,
    double *cell_scale)
{
  size_t count=sample_count(node->active_mask,4);
  int failed=0;
  node->values=calloc(count*5,sizeof(*node->values));
  if(!node->values) return 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(context->table->workers) reduction(|:failed)
#endif
  for(size_t index=0;index<count;++index) {
    unsigned int at[3];double complex value[5];
    decode(index,node->active_mask,4,at);
    if(!context->sampler(context->sampler_context,
        coordinate(&node->bounds,0,at[0],3),
        coordinate(&node->bounds,1,at[1],3),
        coordinate(&node->bounds,2,at[2],3),value)) failed=1;
    else memcpy(node->values+index*5,value,sizeof(value));
  }
  if(failed) return 0;
  *cell_scale=0.;
  for(size_t index=0;index<count*5;++index)
    *cell_scale=fmax(*cell_scale,cabs(node->values[index]));
  context->table->stats.fitted_samples+=count;
  return 1;
}

static int measure_withheld(build_context *context,
    nec4_ground_adaptive_node *node,double cell_scale)
{
  size_t count=sample_count(node->active_mask,3);
  double worst=0.;int failed=0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(context->table->workers) reduction(|:failed) reduction(max:worst)
#endif
  for(size_t index=0;index<count;++index) {
    unsigned int at[3];double point[3];double complex expected[5],actual[5];
    decode(index,node->active_mask,3,at);
    for(unsigned int axis=0;axis<3;++axis)
      point[axis]=node->active_mask&(1u<<axis)?
        coordinate(&node->bounds,axis,2*at[axis]+1,6):node->bounds.lower[axis];
    if(!context->sampler(context->sampler_context,point[0],point[1],point[2],
                         expected)||
       !interpolate_node(node,point[0],point[1],point[2],actual)) failed=1;
    else worst=fmax(worst,scaled_error(actual,expected,cell_scale));
  }
  if(failed) return 0;
  node->withheld_error=worst;
  context->table->stats.withheld_samples+=count;
  context->table->stats.worst_encountered_scaled_error=fmax(
    context->table->stats.worst_encountered_scaled_error,worst);
  return 1;
}

static int build_node(build_context *context,nec4_ground_adaptive_node *node)
{
  double cell_scale;
  unsigned int active_axes=0;
  node->active_mask=active_mask(&node->bounds);
  if(context->table->maximum_cells&&
     context->table->stats.cells>=context->table->maximum_cells) {
    context->table->stats.leaves++;
    context->table->stats.unresolved_leaves++;return 1;
  }
  context->table->stats.cells++;
  if(node->depth>context->table->stats.deepest_level)
    context->table->stats.deepest_level=node->depth;
  if(!sample_fit(context,node,&cell_scale)||
     !measure_withheld(context,node,cell_scale)) return 0;
  if(node->withheld_error<=context->table->tolerance) {
    node->accepted=1;context->table->stats.leaves++;
    context->table->stats.accepted_leaves++;
    context->table->stats.worst_withheld_scaled_error=fmax(
      context->table->stats.worst_withheld_scaled_error,node->withheld_error);
    return 1;
  }
  if(node->depth>=context->table->maximum_depth) {
    context->table->stats.leaves++;
    context->table->stats.unresolved_leaves++;return 1;
  }
  for(unsigned int axis=0;axis<3;++axis)
    if(node->active_mask&(1u<<axis)) active_axes++;
  node->child_count=(unsigned char)(1u<<active_axes);
  for(unsigned int child=0;child<node->child_count;++child) {
    nec4_ground_adaptive_node *next=calloc(1,sizeof(*next));
    unsigned int bit=0;
    if(!next) return 0;
    next->depth=node->depth+1;
    for(unsigned int axis=0;axis<3;++axis) {
      double middle=.5*(node->bounds.lower[axis]+node->bounds.upper[axis]);
      if(node->active_mask&(1u<<axis)) {
        if(child&(1u<<bit)) { next->bounds.lower[axis]=middle;
          next->bounds.upper[axis]=node->bounds.upper[axis]; }
        else { next->bounds.lower[axis]=node->bounds.lower[axis];
          next->bounds.upper[axis]=middle; }
        bit++;
      } else next->bounds.lower[axis]=next->bounds.upper[axis]=
          node->bounds.lower[axis];
    }
    node->children[child]=next;
    if(!build_node(context,next)) return 0;
  }
  return 1;
}

int nec4_ground_adaptive_build(nec4_ground_adaptive *table,
    const nec4_ground_adaptive_bounds *bounds,double tolerance,
    unsigned int maximum_depth,size_t maximum_cells,int workers,
    nec4_ground_adaptive_sampler sampler,void *sampler_context)
{
  build_context context;
  if(!table||!bounds||!sampler||!(tolerance>0.)||!isfinite(tolerance)||
     maximum_depth>20) return 0;
  for(unsigned int axis=0;axis<3;++axis)
    if(!isfinite(bounds->lower[axis])||!isfinite(bounds->upper[axis])||
       bounds->upper[axis]<bounds->lower[axis]) return 0;
  nec4_ground_adaptive_destroy(table);
  table->bounds=*bounds;table->tolerance=tolerance;
  table->maximum_depth=maximum_depth;table->maximum_cells=maximum_cells;
  table->workers=workers>0?workers:1;
  table->root=calloc(1,sizeof(*table->root));
  if(!table->root) return 0;
  table->root->bounds=*bounds;
  memset(&context,0,sizeof(context));context.table=table;
  context.sampler=sampler;context.sampler_context=sampler_context;
  if(!build_node(&context,table->root)) {
    nec4_ground_adaptive_destroy(table);return 0;
  }
  return 1;
}

static const nec4_ground_adaptive_node *find_leaf(
    const nec4_ground_adaptive_node *node,const double point[3])
{
  if(!node||node->child_count==0) return node;
  unsigned int child=0,bit=0;
  for(unsigned int axis=0;axis<3;++axis) if(node->active_mask&(1u<<axis)) {
    double middle=.5*(node->bounds.lower[axis]+node->bounds.upper[axis]);
    if(point[axis]>=middle) child|=1u<<bit;
    bit++;
  }
  return find_leaf(node->children[child],point);
}

int nec4_ground_adaptive_evaluate(const nec4_ground_adaptive *table,
    double rho,double z,double depth,double complex values[5])
{
  const double point[3]={rho,z,depth};
  const nec4_ground_adaptive_node *leaf;
  if(!table||!table->root||!values) return 0;
  for(unsigned int axis=0;axis<3;++axis) {
    double tolerance=64.*DBL_EPSILON*fmax(1.,
      fmax(fabs(table->bounds.lower[axis]),fabs(table->bounds.upper[axis])));
    if(point[axis]<table->bounds.lower[axis]-tolerance||
       point[axis]>table->bounds.upper[axis]+tolerance) return 0;
  }
  leaf=find_leaf(table->root,point);
  return leaf&&leaf->accepted&&interpolate_node(leaf,rho,z,depth,values);
}

#define ADAPTIVE_CACHE_MAGIC "N4ADP1\0"
#define ADAPTIVE_CACHE_VERSION 1u
#define ADAPTIVE_ENDIAN_MARKER 0x01020304u
#define FNV_OFFSET UINT64_C(14695981039346656037)
#define FNV_PRIME UINT64_C(1099511628211)

static uint64_t hash_bytes(uint64_t hash,const void *data,size_t size)
{
  const unsigned char *bytes=data;
  for(size_t index=0;index<size;++index) { hash^=bytes[index];hash*=FNV_PRIME; }
  return hash;
}

static int write_hashed(FILE *stream,const void *value,size_t size,uint64_t *hash)
{
  if(fwrite(value,1,size,stream)!=size) return 0;
  *hash=hash_bytes(*hash,value,size);return 1;
}

static int read_hashed(FILE *stream,void *value,size_t size,uint64_t *hash)
{
  if(fread(value,1,size,stream)!=size) return 0;
  *hash=hash_bytes(*hash,value,size);return 1;
}

static int store_node(FILE *stream,const nec4_ground_adaptive_node *node,
    uint64_t *hash)
{
  uint32_t depth=node->depth;
  uint8_t flags[4]={node->active_mask,node->child_count,node->accepted,
                    node->values!=NULL};
  uint64_t values=node->values?sample_count(node->active_mask,4)*5:0;
  if(!write_hashed(stream,node->bounds.lower,sizeof(node->bounds.lower),hash)||
     !write_hashed(stream,node->bounds.upper,sizeof(node->bounds.upper),hash)||
     !write_hashed(stream,&node->withheld_error,sizeof(node->withheld_error),hash)||
     !write_hashed(stream,&depth,sizeof(depth),hash)||
     !write_hashed(stream,flags,sizeof(flags),hash)||
     !write_hashed(stream,&values,sizeof(values),hash)) return 0;
  for(uint64_t index=0;index<values;++index) {
    double pair[2]={creal(node->values[index]),cimag(node->values[index])};
    if(!write_hashed(stream,pair,sizeof(pair),hash)) return 0;
  }
  for(unsigned int child=0;child<node->child_count;++child)
    if(!store_node(stream,node->children[child],hash)) return 0;
  return 1;
}

static int store_stats(FILE *stream,const nec4_ground_adaptive_stats *stats,
    uint64_t *hash)
{
  uint64_t counts[6]={stats->cells,stats->leaves,stats->accepted_leaves,
    stats->unresolved_leaves,stats->fitted_samples,stats->withheld_samples};
  uint32_t depth=stats->deepest_level;
  return write_hashed(stream,counts,sizeof(counts),hash)&&
    write_hashed(stream,&depth,sizeof(depth),hash)&&
    write_hashed(stream,&stats->worst_withheld_scaled_error,sizeof(double),hash)&&
    write_hashed(stream,&stats->worst_encountered_scaled_error,sizeof(double),hash);
}

int nec4_ground_adaptive_store(const nec4_ground_adaptive *table,
    const char *path,const void *identity,size_t identity_size)
{
  FILE *stream;char *temporary;int descriptor,failed=0;size_t length;
  uint32_t version=ADAPTIVE_CACHE_VERSION,endian=ADAPTIVE_ENDIAN_MARKER;
  uint32_t double_size=sizeof(double),identity_bytes=(uint32_t)identity_size;
  uint32_t maximum_depth;uint64_t maximum_cells,hash=FNV_OFFSET;
  if(!table||!table->root||!path||!*path||(!identity&&identity_size)||
     identity_size>UINT32_MAX) return 0;
  length=strlen(path);temporary=malloc(length+12);if(!temporary)return 0;
  memcpy(temporary,path,length);memcpy(temporary+length,".tmp.XXXXXX",12);
  descriptor=mkstemp(temporary);if(descriptor<0){free(temporary);return 0;}
  stream=fdopen(descriptor,"wb");
  if(!stream){close(descriptor);unlink(temporary);free(temporary);return 0;}
  maximum_depth=table->maximum_depth;maximum_cells=table->maximum_cells;
  failed=fwrite(ADAPTIVE_CACHE_MAGIC,1,8,stream)!=8||
    !write_hashed(stream,&version,sizeof(version),&hash)||
    !write_hashed(stream,&endian,sizeof(endian),&hash)||
    !write_hashed(stream,&double_size,sizeof(double_size),&hash)||
    !write_hashed(stream,&identity_bytes,sizeof(identity_bytes),&hash)||
    (identity_size&&!write_hashed(stream,identity,identity_size,&hash))||
    !write_hashed(stream,table->bounds.lower,sizeof(table->bounds.lower),&hash)||
    !write_hashed(stream,table->bounds.upper,sizeof(table->bounds.upper),&hash)||
    !write_hashed(stream,&table->tolerance,sizeof(table->tolerance),&hash)||
    !write_hashed(stream,&maximum_depth,sizeof(maximum_depth),&hash)||
    !write_hashed(stream,&maximum_cells,sizeof(maximum_cells),&hash)||
    !store_stats(stream,&table->stats,&hash)||
    !store_node(stream,table->root,&hash)||
    fwrite(&hash,1,sizeof(hash),stream)!=sizeof(hash);
  if(!failed&&(fflush(stream)||fsync(fileno(stream)))) failed=1;
  if(fclose(stream)) failed=1;
  if(!failed&&rename(temporary,path)) failed=1;
  if(failed) unlink(temporary);
  free(temporary);return !failed;
}

static int load_stats(FILE *stream,nec4_ground_adaptive_stats *stats,
    uint64_t *hash)
{
  uint64_t counts[6];uint32_t depth;
  if(!read_hashed(stream,counts,sizeof(counts),hash)||
     !read_hashed(stream,&depth,sizeof(depth),hash)||
     !read_hashed(stream,&stats->worst_withheld_scaled_error,sizeof(double),hash)||
     !read_hashed(stream,&stats->worst_encountered_scaled_error,sizeof(double),hash))
    return 0;
  if(counts[0]>SIZE_MAX||counts[1]>SIZE_MAX||counts[2]>SIZE_MAX||
     counts[3]>SIZE_MAX||counts[4]>SIZE_MAX||counts[5]>SIZE_MAX) return 0;
  stats->cells=(size_t)counts[0];stats->leaves=(size_t)counts[1];
  stats->accepted_leaves=(size_t)counts[2];
  stats->unresolved_leaves=(size_t)counts[3];
  stats->fitted_samples=(size_t)counts[4];
  stats->withheld_samples=(size_t)counts[5];stats->deepest_level=depth;
  return 1;
}

static nec4_ground_adaptive_node *load_node(FILE *stream,uint64_t *hash,
    size_t *nodes,size_t limit,unsigned int maximum_depth)
{
  nec4_ground_adaptive_node *node;uint32_t depth;uint8_t flags[4];uint64_t values;
  if(*nodes>=limit) return NULL;
  node=calloc(1,sizeof(*node));if(!node)return NULL;(*nodes)++;
  if(!read_hashed(stream,node->bounds.lower,sizeof(node->bounds.lower),hash)||
     !read_hashed(stream,node->bounds.upper,sizeof(node->bounds.upper),hash)||
     !read_hashed(stream,&node->withheld_error,sizeof(node->withheld_error),hash)||
     !read_hashed(stream,&depth,sizeof(depth),hash)||
     !read_hashed(stream,flags,sizeof(flags),hash)||
     !read_hashed(stream,&values,sizeof(values),hash)) goto failed;
  node->depth=depth;node->active_mask=flags[0];node->child_count=flags[1];
  node->accepted=flags[2];
  if(node->depth>maximum_depth||node->active_mask!=active_mask(&node->bounds)||
     node->child_count>8||flags[3]>1||node->accepted>1||
     values!=(flags[3]?sample_count(node->active_mask,4)*5:0)) goto failed;
  if(values) {
    if(values>SIZE_MAX/sizeof(*node->values)) goto failed;
    node->values=calloc((size_t)values,sizeof(*node->values));if(!node->values)goto failed;
    for(uint64_t index=0;index<values;++index) {
      double pair[2];if(!read_hashed(stream,pair,sizeof(pair),hash))goto failed;
      node->values[index]=pair[0]+I*pair[1];
    }
  }
  for(unsigned int child=0;child<node->child_count;++child) {
    node->children[child]=load_node(stream,hash,nodes,limit,maximum_depth);
    if(!node->children[child]) goto failed;
  }
  return node;
failed:
  destroy_node(node);return NULL;
}

int nec4_ground_adaptive_load(nec4_ground_adaptive *table,const char *path,
    const void *identity,size_t identity_size)
{
  FILE *stream;char magic[8];unsigned char *actual_identity=NULL;
  uint32_t version,endian,double_size,identity_bytes,maximum_depth;
  uint64_t maximum_cells,hash=FNV_OFFSET,stored_hash;size_t nodes=0,limit;
  nec4_ground_adaptive loaded;
  int trailing;
  if(!table||!path||!*path||(!identity&&identity_size)||identity_size>UINT32_MAX)
    return 0;
  stream=fopen(path,"rb");if(!stream)return 0;
  nec4_ground_adaptive_init(&loaded);
  if(fread(magic,1,sizeof(magic),stream)!=sizeof(magic)||
     memcmp(magic,ADAPTIVE_CACHE_MAGIC,sizeof(magic))||
     !read_hashed(stream,&version,sizeof(version),&hash)||
     !read_hashed(stream,&endian,sizeof(endian),&hash)||
     !read_hashed(stream,&double_size,sizeof(double_size),&hash)||
     !read_hashed(stream,&identity_bytes,sizeof(identity_bytes),&hash)||
     version!=ADAPTIVE_CACHE_VERSION||endian!=ADAPTIVE_ENDIAN_MARKER||
     double_size!=sizeof(double)||identity_bytes!=identity_size) goto failed;
  if(identity_bytes) {
    actual_identity=malloc(identity_bytes);if(!actual_identity)goto failed;
    if(!read_hashed(stream,actual_identity,identity_bytes,&hash)||
       memcmp(actual_identity,identity,identity_bytes)) goto failed;
  }
  if(!read_hashed(stream,loaded.bounds.lower,sizeof(loaded.bounds.lower),&hash)||
     !read_hashed(stream,loaded.bounds.upper,sizeof(loaded.bounds.upper),&hash)||
     !read_hashed(stream,&loaded.tolerance,sizeof(loaded.tolerance),&hash)||
     !read_hashed(stream,&maximum_depth,sizeof(maximum_depth),&hash)||
     !read_hashed(stream,&maximum_cells,sizeof(maximum_cells),&hash)||
     !load_stats(stream,&loaded.stats,&hash)) goto failed;
  loaded.maximum_depth=maximum_depth;loaded.maximum_cells=(size_t)maximum_cells;
  if(!(loaded.tolerance>0.)||maximum_depth>20||maximum_cells>SIZE_MAX)goto failed;
  for(unsigned int axis=0;axis<3;++axis)
    if(!isfinite(loaded.bounds.lower[axis])||!isfinite(loaded.bounds.upper[axis])||
       loaded.bounds.upper[axis]<loaded.bounds.lower[axis])goto failed;
  limit=loaded.stats.cells+loaded.stats.unresolved_leaves+1;
  if(limit<loaded.stats.cells)goto failed;
  loaded.root=load_node(stream,&hash,&nodes,limit,maximum_depth);
  if(!loaded.root||fread(&stored_hash,1,sizeof(stored_hash),stream)!=sizeof(stored_hash))
    goto failed;
  trailing=fgetc(stream);if(trailing!=EOF||stored_hash!=hash)goto failed;
  free(actual_identity);fclose(stream);nec4_ground_adaptive_destroy(table);
  *table=loaded;return 1;
failed:
  free(actual_identity);fclose(stream);nec4_ground_adaptive_destroy(&loaded);
  return 0;
}
