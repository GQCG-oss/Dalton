
.. _linking_to_math:

Linking to math libraries
=========================

General
-------

Dalton requires BLAS and LAPACK libraries. Typically you will want to link to
external math (BLAS and LAPACK) libraries, for instance provided by MKL or
Atlas.

By default the CMake configuration script will automatically detect these libraries::

  $ ./setup --blas=auto --lapack=auto          # this is the default

if you define MATH_ROOT, for instance::

  $ export MATH_ROOT='/opt/intel/mkl'

Do not use full path MATH_ROOT='/opt/intel/mkl/lib/ia32'. CMake will append the
correct paths depending on the processor and the default integer type.  If the
MKL libraries that you want to use reside in
/opt/intel/mkl/10.0.3.020/lib/em64t, then MATH_ROOT is defined as::

  $ export MATH_ROOT='/opt/intel/mkl/10.0.3.020'

Then::

  $ ./setup [--flags]
  $ cd build
  $ make


Intel/MKL
---------

If you compile with Intel compilers and have the MKL library available, you
should use the --mkl flag which will automatically link to the MKL libraries
(in this case you do not have to set MATH_ROOT).
You have to specify whether you want to use the sequential or parallel
(threaded) MKL version. For a parallel Dalton runs you should probably link to
the sequential MKL::

  $ ./setup --fc=mpif90 --cc=mpicc --cxx=mpicxx --mkl=sequential

For a sequential compilation you may want to link to the parallel MKL::

  $ ./setup --fc=ifort --cc=icc --cxx=icpc --mkl=parallel

The more general solution is to link to the parallel MKL and control the number
of threads using MKL environment variables.


Cray
----

Cray typically provides own optimized BLAS/LAPACK wrappers.
For this use the option --cray to disable automatic BLAS/LAPACK detection::

  $ ./setup --fc=ftn --cc=cc --cxx=CC --cray


Explicitly specifying BLAS and LAPACK libraries
-----------------------------------------------

If automatic detection of math libraries fails for whatever reason, you can
always call the libraries explicitly like here::

  $ ./setup --blas=/usr/lib/libblas.so --lapack=/usr/lib/liblapack.so

Alternatively you can use the --explicit-libs option. But in this case you should
disable BLAS/LAPACK detection::

  $ ./setup --blas=none --lapack=none --explicit-libs="-L/usr/lib -lblas -llapack"


Builtin BLAS and LAPACK implementation
--------------------------------------

If no external BLAS and LAPACK libraries are available, you can use the builtin
implementation. However note that these are not optimized and you will sacrifice
performance. This should be the last resort if nothing else is available::

  $ ./setup --blas=builtin --lapack=builtin

LSDALTON/ScaLAPACK/Intel/MKL
---------

If you compile with Intel compilers and have the MKL library available, you
can chose to compile LSDALTON using the ScaLAPACK library provided by Intel. 
In this case you should set the MATH_ROOT enviromental variable and use 
the --scalapack flag which will automatically link to the MKL libraries.

You should not use the --mkl flag for this setup::

  $ ./setup --fc=mpif90 --cc=mpicc --cxx=mpicxx --scalapack

