#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dense_backend.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures;
#define CHECK(condition, message) do { if(!(condition)) { \
  fprintf(stderr, "FAIL: %s\n", message); ++failures; } } while(0)

static void fill_system(complex double matrix[9], complex double rhs[3],
                        complex double expected[3])
{
  const complex double source[9] = {
    4.0+I,  1.0-2.0*I, 0.5,
    2.0,    5.0+0.5*I, -I,
    -1.0+I, 0.25,       3.0-I
  };
  size_t row, column;
  expected[0]=1.0+2.0*I;expected[1]=-0.5+I;expected[2]=2.0-0.25*I;
  memcpy(matrix,source,sizeof(source));
  for(row=0;row<3;++row)
  {
    rhs[row]=0.0;
    for(column=0;column<3;++column)
      rhs[row]+=matrix[row*3+column]*expected[column];
  }
}

static void test_row_major(const char *backend, double tolerance)
{
  complex double matrix[9],rhs[3],expected[3];
  size_t index;
  fill_system(matrix,rhs,expected);
  CHECK(dense_backend_set(backend),"backend selection");
  CHECK(dense_backend_solve_row_major(matrix,rhs,3),"row-major solve");
  for(index=0;index<3;++index)
    CHECK(cabs(rhs[index]-expected[index])<tolerance,"row-major solution accuracy");
}

static void test_classic_lapack(void)
{
  complex double row_major[9],rhs[3],expected[3],storage[15],many[10]={0};
  complex double expected_second[3]={-I,2.0+0.5*I,-1.0};
  int pivots[3];
  size_t row,column,index;
  fill_system(row_major,rhs,expected);
  memset(storage,0,sizeof(storage));
  for(row=0;row<3;++row)
    for(column=0;column<3;++column)
      storage[row+column*5]=row_major[column*3+row];
  CHECK(dense_backend_factor_classic(3,storage,5,pivots)==1,
        "strided classic LAPACK factorization");
  CHECK(dense_backend_solve_classic(3,storage,5,pivots,rhs)==1,
        "strided classic LAPACK solve");
  for(index=0;index<3;++index)
    CHECK(cabs(rhs[index]-expected[index])<1e-12,"classic LAPACK solution accuracy");

  for(row=0;row<3;++row)
  {
    many[row]=0.0;many[5+row]=0.0;
    for(column=0;column<3;++column)
    {
      many[row]+=row_major[row*3+column]*expected[column];
      many[5+row]+=row_major[row*3+column]*expected_second[column];
    }
  }
  CHECK(dense_backend_solve_classic_many(3,storage,5,pivots,many,5,2)==1,
        "strided batched classic LAPACK solve");
  for(index=0;index<3;++index)
  {
    CHECK(cabs(many[index]-expected[index])<1e-12,
          "first batched LAPACK solution accuracy");
    CHECK(cabs(many[5+index]-expected_second[index])<1e-12,
          "second batched LAPACK solution accuracy");
  }
}

int main(void)
{
  test_row_major("builtin",1e-12);
  if(dense_backend_lapack_available())
  {
    test_row_major("lapack",1e-12);
    test_classic_lapack();
  }
  CHECK(dense_backend_set("auto"),"auto backend selection");
  CHECK(!dense_backend_set("invalid"),"invalid backend rejection");
  if(failures) return 1;
  puts("Dense builtin/LAPACK backend tests passed");
  return 0;
}
