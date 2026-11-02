#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([nec2c-n5hxr], [1.3.1], [n5hxr@thestone.zone])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([main.c])
AC_CONFIG_HEADERS([config.h])
AM_MAINTAINER_MODE

# Set CFLAGS if not set in environment
AC_MSG_CHECKING([whether configure should try to set CFLAGS])
if test "x${CFLAGS+set}" = xset; then
   enable_cflags_setting=no
else
   enable_cflags_setting=yes
fi
AC_MSG_RESULT($enable_cflags_setting)
test "x$enable_cflags_setting" = xyes && CFLAGS="-Wall -O2"

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AM_PATH_PYTHON([3.6])

# Checks for libraries.
AC_CHECK_LIB([m], [sin])
AC_ARG_WITH([lapack],
  [AS_HELP_STRING([--with-lapack@<:@=auto/yes/no@:>@],
    [enable the optional LAPACK dense backend (default: auto)])],
  [], [with_lapack=auto])
AS_IF([test "x$with_lapack" != xno], [
  save_LIBS="$LIBS"
  AC_CHECK_HEADER([lapacke.h], [have_lapacke_header=yes], [have_lapacke_header=no])
  AS_IF([test "x$have_lapacke_header" = xyes], [
  AC_SEARCH_LIBS([LAPACKE_zgetrf], [lapacke], [
    have_lapack=found
  ], [
    have_lapack=no
    LIBS="$save_LIBS"
  ])
  ], [have_lapack=no])
  AS_IF([test "x$have_lapack" = xfound], [
    AC_MSG_CHECKING([whether LAPACKE complex-double solves are compatible])
    AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <complex.h>
#include <lapacke.h>
#include <math.h>
int main(void) {
  lapack_int n=2,nrhs=1,lda=2,ldb=2,ipiv[2],info;
  double _Complex a[4]={2.0+I,1.0,3.0,4.0-I};
  double _Complex expected[2]={1.0+2.0*I,-0.5+I};
  double _Complex b[2];
  b[0]=(2.0+I)*expected[0]+1.0*expected[1];
  b[1]=3.0*expected[0]+(4.0-I)*expected[1];
  info=LAPACKE_zgetrf(LAPACK_COL_MAJOR,n,n,
      (lapack_complex_double*)a,lda,ipiv);
  if(info) return 1;
  info=LAPACKE_zgetrs(LAPACK_COL_MAJOR,'T',n,nrhs,
      (lapack_complex_double*)a,lda,ipiv,(lapack_complex_double*)b,ldb);
  return info || cabs(b[0]-expected[0])>1e-12 || cabs(b[1]-expected[1])>1e-12;
}]])], [have_lapack=yes], [have_lapack=no; LIBS="$save_LIBS"],
      [have_lapack=yes; AC_MSG_WARN([cross compiling: assuming LAPACKE complex support])])
    AC_MSG_RESULT([$have_lapack])
  ])
  AS_IF([test "x$have_lapack" = xyes], [
    AC_DEFINE([HAVE_LAPACK], [1],
      [Define when a compatible complex-double LAPACK library is available.])
    AS_IF([test "x$ac_cv_search_LAPACKE_zgetrf" != "xnone required"], [
      PACKAGE_LIBS="$ac_cv_search_LAPACKE_zgetrf $PACKAGE_LIBS"
    ])
  ])
  AS_IF([test "x$with_lapack" = xyes && test "x$have_lapack" != xyes], [
    AC_MSG_ERROR([--with-lapack=yes requested, but compatible LAPACKE headers and library were not found])
  ])
], [have_lapack=no])
AC_MSG_NOTICE([optional LAPACK dense backend: $have_lapack])
AC_SUBST(PACKAGE_CFLAGS)
AC_SUBST(PACKAGE_LIBS)

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
