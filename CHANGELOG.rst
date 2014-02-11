

2013.2 (this will be the next patch)
====================================

COMMON
------

- Recognize CYGWIN as a LINUX and UNIX system, for proper definition of compilation flags
- Define M_PI in C-code if not already defined (problem seen with Cygwin)


DALTON
------

- Fixed a bug in printing results in CPP-QRF.
- More efficient evaluation of numerical Hessian when C1 symmetry
  (in each geometry step start wave function optimization from a
  converged wave function from a neighboring geometry rather than from scratch each time).
- Fix of error which sometimes caused a geometry optimization to stop with " *** ERROR, Wrong interval in WLKBIS".
- Fix of a bug which occasionally caused DALTON to abort a .STEX calculation


LSDALTON
--------

- Print sensible error message when running out of memory.
- Added funcitonality to search through several basis-set libraries.
- Increased max length of WRKDIR from 60 to 200.
- Fixed a bug related to improber shutdown of MPI calculation. In the case
  of wrong LSDALTON.INP for instance the calculation will issue a error 
  statement and afterward hang forever in a MPI call. 
 









-------------------------------------------------------------------------------------
-
-                         PREVIOUS PATCHES STARTS HERE
-
-------------------------------------------------------------------------------------

2013.1 (2013-12-19)
===================

DALTON
------

- Correct the printout of relativistic corrections to the shielding (thanks to M. Jaszunski).
- Compilation fix for DALTON/abacus/rma_windows.F90 (Intel 10.0.011).
- Fix of error where basis set names were changed to upper case and could not be found (reported by Yurij Rusakov).
- Each MPI slave sleeps 10 millisecond between tests for new task
  (only Intel; should enable turbomode in sequential parts of DALTON, and more efficient use of threaded MKL when combined with MPI).
- added metric scaled output of orbital response vectors in \*\*RESPONS
  (for easier interpretation of excitation operators).


LSDALTON
--------

- Fixed a bug in Jengine, related to screening for nonsymmetric density matrices.
  This may affect CCSD and some response calculation. 
- Modified the input section of the manual concerning 
  Casida-Salahub asymptotic correction CS00 (thanks to Raul Crespo).
- Changed defaults for Casida-Salahub asymptotic correction CS00 (thanks to Raul Crespo).
- Fixed errors in the MCD B terms output files (.dat files) now one file is generated
  for each B term and each A term (thanks to Raul Crespo) .
- Modified the input section of the manual concerning MCD B terms. Added desciption of MCDEXSTATES.
- Fixed a bug for lsdalton geometry optimization and dynamics related to 
  screening. The initial Cauchy-Schwartz screening matrices were incorrectly
  used in each subsequent geometry step
