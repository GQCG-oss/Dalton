

=============================
Dalton |version| installation
=============================

Dalton is configured using CMake, typically via the setup script,
and subsequently compiled using make or gmake.


Basics
------

The setup script is a useful front-end to CMake. To see all options, run::

  $ ./setup --help

The setup script does nothing else than creating the directory "build" and
calling CMake with appropriate environment variables and flags::

  $ ./setup [--flags]
  $ cd build
  $ make

By default CMake builds out of source. This means that all object files and the
final binary are generated outside of the source directory. Typically the build
directory is called "build", but you can change the name of the build directory
(e.g. "build_gfortran")::

  $ ./setup [--flags] build_gfortran
  $ cd build_gfortran
  $ make

By default we compile all targets (DALTON and LSDALTON). Instead of typing
``make``, you can restrict compilation to only ``dalton.x`` or ``lsdalton.x``
like this::

  $ make dalton.x
  $ make lsdalton.x

You can compile the code on several cores::

  $ make -j4


Most typical examples
---------------------

In order to get familiar with the configuration setup, let us demonstrate
some of the most typical configuration scenarios.

Configure for parallel compilation using MPI (make sure to properly export MPI
paths)::

  $ ./setup --fc=mpif90 --cc=mpicc --cxx=mpicxx

There is a shortcut for it::

  $ ./setup --mpi

Configure for sequential compilation using ifort/icc/icpc::

  $ ./setup --fc=ifort --cc=icc --cxx=icpc

Configure for sequential compilation using gfortran/gcc/g++::

  $ ./setup --fc=gfortran --cc=gcc --cxx=g++

You get the idea. The configuration is normally good at detecting math libraries
automatically, provided you export the proper environment variable ``MATH_ROOT``
(see further below).


What to do if CMake is not available or too old?
------------------------------------------------

If it is your machine and you have an Ubuntu or Debian-based distribution::

  $ sudo apt-get install cmake

On Fedora::

  $ sudo yum install cmake

Similar mechanisms exist for other distributions or
operating systems. Please consult Google.

If it is a cluster, please ask the Administrator to install/upgrade CMake.

If it is a cluster, but you prefer to install it yourself (it's easy):

1. Download the latest precompiled tarball from http://www.cmake.org/cmake/resources/software.html
2. Extract the tarball to, say, ~/cmake-2.8.5-Linux-i386
3. Set correct PATH variable


=========================
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
performace. This should be the last resort if nothing else is available::

  $ ./setup --blas=builtin --lapack=lapack


=================
Scratch directory
=================

TODO: radovan will write this ...


===================
Basis set directory
===================

TODO: radovan will write this ...


==============
Expert options
==============


Compiling in verbose mode
-------------------------

Sometimes you want to see the actual compiler flags and definitions::

  $ make VERBOSE=1


How can I change optimization flags?
------------------------------------

You can turn optimization off (debug mode) like this::

  $ ./setup --type=debug [other flags]
  $ cd build
  $ make

You can edit compiler flags in cmake/compilers/{FortranFlags.cmake, CFlags.cmake, CXXFlags.cmake}.

Alternatively you can edit compiler flags through ccmake::

  $ cd build
  $ ccmake ..


========================
Testing the installation
========================

It is very important that you verify that your Dalton installation correctly
reproduces the reference test set before running any production calculations.

The test set driver is CTest which can be invoked with "make test".


Environment variables for testing
---------------------------------

Before testing with "make test" you should export the
following environment variables::

  $ export DALTON_TMPDIR=/scratch        # scratch space for Dalton (adapt the path of course)
  $ export CTEST_MAKE_NUM_PROCS=16       # in this case the code will be compiled with 16 processes (make -j16)
  $ export DALTON_NUM_MPI_PROCS=4        # in this case 4 processes, only relevant if you compile with MPI

Note that if you set the DALTON_NUM_MPI_PROCS to something different from 1,
the dalton script will assume you have compiled using MPI and run the mpirun
command!


Running the test set
--------------------

You can run the whole test set either using::

  $ make test

or directly through CTest::

  $ ctest

Both are equivalent ("make test" runs CTest) but running
CTest directly makes it easier to run sequential tests on several
cores::

  $ ctest -j4

You can select the subset of tests by matching test names to a regular expression::

  $ ctest -R dft

Alternatively you can select the tests with a label matching a regular expression::

  $ ctest -L rsp

The following command will give you all available labels::

  $ ctest --print-labels


Running only DALTON or only LSDALTON tests
------------------------------------------

Only DALTON tests::

  $ ctest -L dalton

Only LSDALTON tests::

  $ ctest -L linsca

All tests::

  $ ctest
