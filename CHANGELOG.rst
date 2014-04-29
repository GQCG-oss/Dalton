

This is the next patch, document new changes here
=================================================

Common
------


DALTON
------

- Default DIIS space increased from 5 to 8, often resulting in 1-2 fewer SCF iterations.
- Removed the maximum of 20 excitations in summary output for second and third order transition moments.
- Bugfix for parallel calculations and some type of geometry optimizations with ano basis sets
  (this bug resulted in aborted calculations, not in wrong results)
- Print irrep names together with symmetry numbers for easier interpretation of output
- More important output with '@' in column 1 (can be obtained with 'grep @' on the output)
- Environment variable DALTON_USE_GLOBAL_SCRATCH disables copying of binaries to worker nodes.
- Environment variable DALTON_LAUNCHER introduced.
- Fixed output information about number of MPI processes and number of OpenMP threads.
- Added information in the error messages when values in maxorb.h are exceeded (which values to increase).


LSDALTON
--------

- Environment variable LSDALTON_LAUNCHER introduced.


===========================
DO NOT EDIT BELOW THIS LINE
===========================


2013.2 (2014-03-05)
===================

Common
------

- Recognize CYGWIN as a LINUX and UNIX system, for proper definition of compilation flags.
- Define M_PI in C-code if not already defined (problem seen with Cygwin).
- Added setup option --blacs to be used in combination with --scalapack; defaults to --blacs=intelmpi.


DALTON
------

- Fixed a bug in printing results in CPP-QRF.
- New CPP solver works also for non-direct calculation.
- More efficient evaluation of numerical Hessian when C1 symmetry
  (in each geometry step start wave function optimization from a
  converged wave function from a neighboring geometry rather than from scratch each time).
- Fix of error which sometimes caused a geometry optimization to stop with " *** ERROR, Wrong interval in WLKBIS".
- Fix of a bug which occasionally caused DALTON to abort a .STEX calculation.
- Print final geometry in xyz format (angstrom). File called "final_geometry.xyz" is put into the restart tarball.
- Append PID to scratch directory to avoid multiple tests running in the same directory.
- Improved manual for two-photon and non-adiabatic coupling.
- Updated/corrected g-factors for Ag, Nd, and Tl (thanks to M. Jaszunski).


LSDALTON
--------

- Print sensible error message when running out of memory.
- Added functionality to search through several basis-set libraries.
- Increased max length of WRKDIR from 60 to 200.
- Fixed a bug related to improper shutdown of MPI calculation. In the case
  of wrong LSDALTON.INP for instance the calculation will issue a error
  statement and afterward hang forever in a MPI call.
- Fixed an OpenMP bug in the calculation of how much memory there should be used during
  an exchange-correlation calculation - resulting in huge memory usage for large molecular system.


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
  for each B term and each A term (thanks to Raul Crespo).
- Modified the input section of the manual concerning MCD B terms. Added description of MCDEXSTATES.
- Fixed a bug for LSDALTON geometry optimization and dynamics related to
  screening. The initial Cauchy-Schwartz screening matrices were incorrectly
  used in each subsequent geometry step
