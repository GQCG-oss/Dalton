

2015.1 (2015-07-XX)
===================

Common
------

- Added ANO-RCC basis set.


DALTON
------

- Update PElib (v.1.2.3): Workaround for faulty system detection using Macports CMake.
- Fixed a bug with Intel Compiler 15 during initialization of Cauchy-Schwarz parameters.
- Fixed a bug for parallel build on some systems.
- Fixed a segmentation fault for approx. 3000 basis functions
  (an array was allocated in stack memory, and became too big for default size of stack memory).
- Fixed a bug in CC sumrules.
- Fixed a bug in the preoptimization, i.e. when using smaller basis sets first to geometry optimize molecule.
- Fixed some far from optimal defaults for preoptimization.
- Fixed geometry optimization for HS-ROHF and HS-RODFT with symmetry - .SINGLY input option under '\*SCF INPUT'
  (use numerical gradients as analytical gradients are only implemented without symmetry).
- Some minor corrections to the Dalton manual.
- More reasonable output for TPCD.


LSDALTON
--------

- Fixed .UNCONT for EMSL's Dalton basis set format.


2015.0 (2015-02-18)
===================

See the release notes at http://daltonprogram.org for a list of new features in
Dalton and LSDalton.


2013.4 (2014-07-10)
===================

DALTON
------

- Memory bugfix for serial PCM calculations (segmentation fault for large PCM cavities).


LSDALTON
--------

- Fixed a bug in the basis set reading. This bugfix affects almost no basis sets,
  and none of the standard basis sets, but a very few general contracted basis sets
  where the first contracted function had much smaller number of
  primitives compared to the last: Basis sets such as the pcS-1 basis set.
- Reduced the memory requirements for internal MPI buffer handling.


2013.3 (2014-06-11)
===================

Common
------

- aug-cc-pVTZ-lresc basis set added to $BASDIR.


DALTON
------

- Default DIIS space increased from 5 to 8, often resulting in 1-2 fewer SCF iterations.
- Removed the maximum of 20 excitations in summary output for second and third order transition moments.
- Warning is issued when orbitals are deleted due to linear dependencies (before SCF),
  AngPso (a 0th order LRESC diamagnetic corr) is not calculated in this case.
- Bugfix for parallel calculations and some type of geometry optimizations with ANO basis sets
  (this bug resulted in aborted calculations, not in wrong results).
- Print irrep names together with symmetry numbers for easier interpretation of output.
- More important output with '@' in column 1 (can be obtained with 'grep @' on the output).
- Environment variable DALTON_USE_GLOBAL_SCRATCH disables copying of binaries to worker nodes.
- Environment variable DALTON_LAUNCHER introduced.
- Fixed output information about number of MPI processes and number of OpenMP threads.
- Added information in the error messages when values in maxorb.h are exceeded (which values to increase).
- Increased some of the values in the common blocks:
  MXSHEL 1000 -> 1500; MXCORB 2400 -> 5000; MXPRIM 8000 -> 15000;
  MAXOCC 800 -> 1500; MXCENT 200 -> 500; MXCENT_QM 200 -> 500
  (the static size of dalton.x went from 100 MB to 165 MB).
- Do not print garbage non-zero transition moments and oscillator strengths for triplet excitations (\*EXCITA module).
- Corrected input description for transition moments between excited states (\*QUADRA with .DOUBLE RESIDUE).
- Fix for \*\*RESPONSE .EXMOM .ISPABC=1,0,1 (only half the excited state spin-orbit transition moments were calculated).
- Fix for Molden file when exponent greater than 1.0D8.
- Fix for MNF-SO (amfi) if more than 40 nuclei.
- Bugfix in quadratic response function using CPP in the tensor contraction routine of the A[2] terms.
- Added interface to ChemShell.
- Bugfix for small non-default WORK array sizes. For specific small custom values of the WORK array size
  KBLOCK was larger than MXBLCK leading to unpredictable results due to array length mismatch in DALTON/abacus/herrdn.F.


LSDALTON
--------

- Environment variable LSDALTON_LAUNCHER introduced.


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
