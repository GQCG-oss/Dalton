#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.53)
dnl AC_INIT with arguments is too new. We use AC_INIT+AM_INIT_AUTOMAKE for
dnl now
dnl AC_INIT(qcsparsemat, 1.0.0, [emanuel@theochem.kth.se])
AC_INIT
AM_INIT_AUTOMAKE(csparsemat, 1.0.0)
AC_CONFIG_SRCDIR([Block.cpp])
dnl AM_CONFIG_HEADER([config.h])

dnl #####################################################################
dnl 2. Options
dnl #####################################################################
AC_MSG_CHECKING([which integers to use in the Fortran interface])

AC_ARG_ENABLE(64bit-int,
        AC_HELP_STRING([--enable-64bit-int(=no)],
                [Use 64-bit integers in the Fortran interface]),
	AC_DEFINE(VAR_INT64,1,[Use 64-bit integers in the Fortran interface])	
	inttype="64bit", inttype="default")
AC_MSG_RESULT($inttype)
# 3. Checks for programs.
AC_PROG_F77([xlf ifort f77 gfortran g77])
AC_PROG_CXX([xlC icpc c++ g++])
AC_PROG_CC([xlc icc cc gcc])
AC_PROG_RANLIB

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lcxa':
AC_CHECK_LIB([cxa], [main])
AC_CHECK_LIB([m], [pow])
AC_CHECK_LIB([pthread], [pthread_create])

AC_F77_LIBRARY_LDFLAGS

# work around pgf77 configuration problem needed on fimm.bccs.uib.no
if echo $FLIBS | grep -q "'"; then
   FLIBS=`echo $FLIBS | sed "s,',,g"`
fi

# FIXME: Replace `main' with a function in `-llapack':
#AC_CHECK_LIB([acml],[main],,
#  AC_CHECK_LIB([mkl],[main],,
#    AC_CHECK_LIB([essl],[main],,
#      AC_CHECK_LIB(f77blas, dswap_ ,,, -latlas)))
#  AC_CHECK_LIB([mkl_lapack],[main],,
#    AC_CHECK_LIB([lapack],[main],,
#      AC_CHECK_LIB([lapack],[dsygv],[LDFLAGS="-L/usr/local/lib $LDFLAGS"],
#      [AC_ERROR([lapack required but not found])],[-L/usr/local/lib $FLIBS]))))

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
dnl AC_HEADER_STDBOOL
AC_C_CONST
dnl code is written in C++, C inline defection just confuses...
dnl AC_C_INLINE
AC_STRUCT_TM

# Checks for library functions.
AC_CHECK_FUNCS([pow],,[AC_ERROR([pow() function needed but not found])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
