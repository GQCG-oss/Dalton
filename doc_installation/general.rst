

Basic installation
==================

General
-------

Dalton is configured using CMake, typically via the ``setup`` script,
and subsequently compiled using make or gmake.
The ``setup`` script is a useful front-end to CMake.
You need python to run ``setup``.


Basics
------

To see all options, run::

  $ ./setup --help

The ``setup`` script does nothing else than creating the directory "build" and
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


Typical examples
----------------

In order to get familiar with the configuration setup, let us demonstrate
some typical configuration scenarios.

Configure for parallel compilation using MPI (make sure to properly export MPI
paths)::

  $ ./setup --fc=mpif90 --cc=mpicc --cxx=mpicxx

There is a shortcut for it::

  $ ./setup --mpi

Configure for sequential compilation using ifort/icc/icpc and link against parallel mkl::

  $ ./setup --fc=ifort --cc=icc --cxx=icpc --mkl=parallel

Configure for sequential compilation using gfortran/gcc/g++::

  $ ./setup --fc=gfortran --cc=gcc --cxx=g++

Parallel compilation on a Cray::

  ./setup --fc=ftn --cc=cc --cxx=CC --cray --mpi

Parallel compilation on a SGI using Intel compilers and MPT::

  ./setup --fc=ifort --cc=icc --cxx=icpc --sgi-mpt

You get the idea. The configuration is usually good at detecting math libraries
automatically, provided you export the proper environment variable ``MATH_ROOT``,
see :ref:`linking_to_math`.


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

1. Download the latest pre-compiled tarball from http://www.cmake.org/cmake/resources/software.html
2. Extract the tarball to, say, ~/cmake-2.8.11-Linux-i386
3. Set correct PATH variable
