

=======================================================
Supported and tested platforms, compilers and libraries
=======================================================


DALTON
======


Tested on the following platforms
---------------------------------

* CentOS 6.3, Intel 13.1.3, MKL 11.0.5, OpenMPI 1.6.3
* CentOS 6.3, GNU 4.6.3, MKL 11.0.5, OpenMPI 1.6.3
* CentOS 6.3, GNU 4.4.6, builtin BLAS/LAPACK, OpenMPI 1.6.4
* CentOS 6.3, GNU 4.4.6, ATLAS 3.8.4, OpenMPI 1.6.4
* CentOS 6.3, GNU 4.4.6, system native BLAS/LAPACK 3.2.1, OpenMPI 1.6.4
* CentOS 6.4, Rocks Linux 6.0, GNU 4.4.7, MKL 11.0.0, OpenMPI 1.6.2
* CentOS 6.4, Rocks Linux 6.1, GNU 4.6.2, MKL 11.0.0, OpenMPI 1.6.4
* CentOS 6.4, Rocks Linux 6.0, GNU 4.7.2, MKL 11.0.0, OpenMPI 1.6.2
* CentOS 6.4, Rocks Linux 6.0, Intel 12.1.2, MKL 11.0.0
* CentOS 6.4, Rocks Linux 6.0, Intel 13.1.2, MKL 11.0.0
* CentOS 6.4, Rocks Linux 6.0, Intel 13.0.1, MKL 11.0.0, OpenMPI 1.6.2
* CentOS 6.4, Rocks Linux 6.1, Intel 13.0.1, MKL 11.0.0, OpenMPI 1.6.5
* CrunchBang 11, GNU 4.7.2, system native BLAS/LAPACK, OpenMPI 1.6.5
* Ubuntu 12.04, GNU 4.6.3, MKL 11.0, OpenMPI 1.6.5
* Fedora 15, GNU 4.6.0, system native BLAS/LAPACK
* Fedora 18, GNU 4.7.2, MKL 11.0.5 
* Fedora 18, Intel 13.1.3 MKL 11.0.5 
* Fedora 18, GNU 4.7.2, ATLAS 3.8.4
* Scientific Linux 6.0, GNU 4.4.4, ATLAS 64-bit version
* MAC OS X 10.9, GNU 4.8.2, Apple vecLib, MPICH 3.1rc1
* Cray XE6, Intel 12.1.5


Compiles but does not pass all tests
------------------------------------

* Fedora 15, Intel 11.4.191, OpenMPI 1.6.5 (failing tests: cc2_r12_aux, rsp_hyperpolar, rsp_dckerr, rsp_excipolar, geoopt_dckerr)


Compiles but untested on the following environments
---------------------------------------------------

* Ubuntu 12.04, Intel 14.0.0, MKL 11.0
* SGI Altix ICE 8200, Intel 11.1


Currently does not compile (but we plan to patch it)
----------------------------------------------------

* IBM AIX, XL 7.1


LSDALTON
========


Tested on the following platforms
---------------------------------

* CentOS 6.3, Intel 13.1.3, MKL 11.0.5, OpenMPI 1.6.3
* CentOS 6.3, GNU 4.6.3, MKL 11.0.5, OpenMPI 1.6.3
* CentOS 6.3, GNU 4.4.6, builtin BLAS/LAPACK, OpenMPI 1.6.4
* CentOS 6.3, GNU 4.4.6, ATLAS 3.8.4, OpenMPI 1.6.4
* CentOS 6.3, GNU 4.4.6, system native BLAS/LAPACK 3.2.1, OpenMPI 1.6.4
* CrunchBang 11, GNU 4.7.2, system native BLAS/LAPACK, OpenMPI 1.6.5
* Ubuntu 12.04, GNU 4.6.3, MKL 11.0, OpenMPI 1.6.5
* Fedora 15, GNU 4.6.0, system native BLAS/LAPACK
* Fedora 15, Intel 11.4.191, OpenMPI 1.6.5
* Fedora 18, GNU 4.7.2, MKL 11.0.5 
* Fedora 18, Intel 13.1.3 MKL 11.0.5 
* Fedora 18, GNU 4.7.2, ATLAS 3.8.4
* CentOS 6.3 Linux 64 bit, GNU 4.4.3, MKL 10.2.7
* CentOS 6.3 Linux 64 bit, GNU 4.4.3, MKL 10.2.7, OpenMPI 1.6.3
* CentOS 6.3 Linux 64 bit, PGI 12.4-0, ACML 4.1.0
* CentOS 6.3 Linux 64 bit, Intel 12.4-0, MKL 10.1.1.019
* CentOS 6.3 Linux 64 bit, Intel 11.1, MKL 10.2.7, OpenMPI 1.4.4
* CentOS 6.3 Linux 64 bit, Intel 13.0.1, MKL 11.0.1, OpenMPI 1.4.1
* MAC OS X 10.9, GNU 4.8.2, Apple vecLib, MPICH 3.1rc1


Compiles but untested on the following environments
---------------------------------------------------

* Ubuntu 12.04, Intel 14.0.0, MKL 11.0
* SGI Altix ICE 8200, Intel 11.1
* Cray XE6, Intel 12.1.5
