

Basic installation using the setup script
-----------------------------------------

Dalton can be configured using CMake as an alternative to configure.
The setup script is a useful front-end to CMake.
It does nothing else than creating the directory "build" and calling
CMake with appropriate environment variables and flags::

  $ ./setup [--flags]
  $ cd build
  $ make

Call setup without flags to see all available options::

  $ ./setup

You can see the CMake command using::

  $ ./setup [--flags] --show

and use it directly to call CMake without setup.


Basic installation without the setup script
-------------------------------------------

Quick start::

  $ mkdir build
  $ cd build
  $ cmake ..
  $ make

The setup script does nothing else than calling CMake with appropriate
environment variables and flags. The two following strategies are completely
equivalent:

Using CMake directly::

  $ mkdir build
  $ cd build
  $ FC=mpif90 CC=mpicc cmake -DENABLE_MPI=1 -DCMAKE_BUILD_TYPE=Release ..
  $ make

Using setup::

  $ ./setup --fc=mpif90 --cc=mpicc --mpi
  $ cd build
  $ make

If the compiler contains "mpi", then you can omit the flag --mpi, setup will set
it in this case automatically.

There is nothing special about the directory "build".
You can do this instead::

  $ mkdir /buildpath
  $ cd /buildpath
  $ cmake /sourcepath
  $ make


Linking to external math libraries
----------------------------------

Typically you will want to link to external math (BLAS and LAPACK) libraries,
for instance provided by MKL or Atlas.

The CMake configuration script will automatically find them if you define MATH_ROOT::

  $ export MATH_ROOT='/opt/intel/mkl'

Do not use full path MATH_ROOT='/opt/intel/mkl/lib/ia32'. CMake will append the
correct paths depending on the processor and the default integer type.  If the
MKL libraries that you want to use reside in
/opt/intel/mkl/10.0.3.020/lib/em64t, then MATH_ROOT is defined as::

  $ export MATH_ROOT='/opt/intel/mkl/10.0.3.020'

Then::

  $ ./setup [--flags]                 # do not need to define --math
  $ cd build
  $ make

Alternatively::

  $ cd build
  $ [FC=gfortran CC=gcc] MATH_ROOT='/opt/intel/mkl' cmake ..
  $ make

Exporting MATH_ROOT is equivalent to calling setup with --math-dir::

  $ ./setup --math-dir=/opt/intel/mkl

If automatic detection of math libraries fails for whatever reason, you can
always call the libraries explicitly like here::

  $ ./setup --math="-L/path -lfoo -lbar"


Running CMake using GUI
-----------------------

You prefer GUI? No problem. You can configure with GUI::

  $ cd build
  $ cmake ..
  $ cmake-gui ..

You may have to install cmake-gui for it, on debian/ubuntu::

  $ sudo apt-get install cmake cmake-gui


Running tests
-------------

You can run the test suite with::

  $ make test

It is HIGHLY recommended to run the test set after you have compiled
Dalton.


Compiling in verbose mode
-------------------------

Sometimes you want to see the actual compiler flags and definitions::

  $ make VERBOSE=1


Compiling on many cores
-----------------------

Yes, it works. Try::

  $ make -j4


Out of source compilation
-------------------------

By default CMake builds out of source.
You can build several binaries with the same source::

  $ cd /sourcepath
  $ ./setup --fc=gfortran --cc=gcc --build=/gfortran-buildpath
  $ cd /gfortran-builddir
  $ make
  $ cd /sourcepath
  $ ./setup --fc=ifort --cc=icc --build=/ifort-buildpath
  $ cd /ifort-buildpath
  $ make


Make install
------------

Make install is very useful to make Dalton available to other users on the same
machine::

  $ ./setup [--flags] --install=/path
  $ cd build
  $ make
  $ make install
