

=======================================================
Supported and tested platforms, compilers and libraries
=======================================================


DALTON
======


Tested on the following platforms
---------------------------------

* Horseshoe cluster (Linux x86_64), Intel 13.1.3, MKL 11.0.5, OpenMPI 1.6.3
* Horseshoe cluster (Linux x86_64), GNU 4.6.3, MKL 11.0.5, OpenMPI 1.6.3
* Horseshoe cluster (Linux x86_64), GNU 4.4.6, builtin math libraries
* CrunchBang 11, GNU 4.7.2, system native BLAS/LAPACK, OpenMPI 1.6.5
* Ubuntu 12.04, GNU 4.6.3, MKL 11.0, OpenMPI 1.6.5


Compiles but does not pass all tests
------------------------------------

* Linux x86_64, Intel 11.4.191, OpenMPI 1.6.5 (failing tests: cc2_412_aux, rsp_hyperpolar, rsp_dckerr, rsp_excipolar, geoopt_dckerr)


Compiles but untested on the following platforms
------------------------------------------------

* SGI Altix ICE 8200 (Hyperion)
* Cray XE6 (Lindgren)


Currently does not compile (but we plan to patch it)
----------------------------------------------------

* IBM AIX, XL 7.1


LSDALTON
========


Tested on the following platforms
---------------------------------

* Horseshoe cluster (Linux x86_64), Intel 13.1.3, MKL 11.0.5, OpenMPI 1.6.3
* Horseshoe cluster (Linux x86_64), GNU 4.6.3, MKL 11.0.5, OpenMPI 1.6.3
* Horseshoe cluster (Linux x86_64), GNU 4.4.6, builtin math libraries
* CrunchBang 11, GNU 4.7.2, system native BLAS/LAPACK, OpenMPI 1.6.5
* Ubuntu 12.04, GNU 4.6.3, MKL 11.0, OpenMPI 1.6.5


Compiles but untested on the following platforms
------------------------------------------------

* SGI Altix ICE 8200 (Hyperion)
* Cray XE6 (Lindgren)
