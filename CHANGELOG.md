# DALTON Change Log -- All notable changes to the DALTON project will be documented in this file.

## special changes for MC-srDFT

### Fixed
- Oct. 2016: state specific MC-srDFT for excited states




## [2017.alpha Unreleased]

### New features added
- New parallel 2-electron integral transformation .FCKTRA  (H. J. Aa. Jensen).
- Triplet excitation energies and polarizabilities with the AO-SOPPA code (P. A. B. Haase).
- Analytical PE-HF/DFT molecular gradients (N. H. List).
  (See e.g. test/pehf_geoopt for example of input.)
- Freezing atoms in geometry optimization (N. H. List and H. J. Aa. Jensen).
  For example for freezing capping atoms in PE geometry optimization.
  (See e.g. test/pehf_geoopt for example of input.)
- Effective external field (EEF) for one- and two-photon absorption in PE-HF/DFT calculations.
 - Reference: N. H. List, H. J. Aa. Jensen, and J. Kongsted. Local Electric Fields and Molecular Properties in Heterogeneous Environments through Polarizable Embedding.
   Phys. Chem. Chem. Phys. 18, 10070 (2016).
- Remove the most diffuse virtual orbitals after SCF or MCSCF (new .VIRTRUNC option)
- Add purely classical multipole-multipole interaction energy in PE-QM calculations

### Fixed
- open-shell doublet ROKS DFT geometry optimization
- Fix .GSPOL for parallel polarizable embedding quadratic response
- Fix text about elimination of some two-photon transitions between excited states
  because they were duplicates (text had "Third order" instead of "Second order")
- Bugfixes for CC2 in environments (incl. external fields) by Ove Christiansen

### Changed
- OK to run ECD or OECD with SOPPA
- More documentation of .STEX in manual.

### Removed


Do not make changes below this line! For the relase branch only.

## [2016.2] (2016-07-12)

### Added
- Added and documented Basis=INTGRL option for ATOMBASIS in .mol file.
- Included Be in cc-pV5Z basis set

### Fixed
- More robust code for reading exponents and contraction coefficients in Dalton-type basis set files, incl. such files from EMSL
- Work-around for Intel 15 compiler I/O problem in some response calculations
- Fix for spin-orbit coupling (SOC) between S/T excited states of same symmetry (problem reported on daltonforum.org)
- Further fixes of MCSCF in **PROPERTIES for more than 255 basis functions - hopefully it is OK now for all requests.
- Fixed an error in the manual for spin-dipole (problem reported on daltonforum.org)
- Fix of open-shell Hartree-Fock occupation output (only output, not the calculation, was wrong if ROHF was followed by MCSCF)
- Fix of Douglas-Kroll post-SCF with less than 256 contracted basis functions, but more than 255 uncontracted basis functions
- Fix of an insufficient memory error for construction of 2-el. integrals in Dirac format with more than 255 basis functions
- Removed OpenACC CMake variable (currently no OpenACC directives in Dalton).


## [2016.1] (2016-04-07)

### Added
- Possibility to read basis set files as made by the EMSL web site
  (this makes it possible to also read basis set files in basis/ based on
  emsl output as e.g. aug-pcseg-1; only LSDALTON has been able to read them so far)

### Fixed
- MCSCF in **PROPERTIES for more than 255 basis functions (fixes problem with MCSCF shielding reported on daltonforum.org)
- Make sure molecule is not moved in ADDSYM during numerical differentiation
- Fixed error in the printing of the cpu/wall time used in Sirius
- Fixed error in PBEc functional: gave NaN when rho was zero.
- Polished some format statements to reduce number of compiler warnings 
- Fixed error in memory addressing for MCSCF g-tensor calculations
- Fixed 2 errors in author list for WIRE dalton publication in dalton output
- Removed unsupported configure options.


## [2016.0] (2015-12-22)

### Changed
- Separated Dalton and LSDalton.

### Added
- New faster CC3 module (Main authors: Rolf H. Myhre and Henrik Koch)
- Parallel AO-SOPPA (Main authors: Frederik Beyer Kjær Hansen and Rasmus Faber)
- QFIT library - fitting charges and, if desired, dipoles to simulate the electrostatic potential from a QM wave function: (Main author: Casper Steinmann)
- Vibrational averaging of NMR coupling constants at SOPPA/MP2 level (Main author: Rasmus Faber)

### Fixed (since version 2015.1)
- Ahlrichs-TZV basis: Fixed error in an exponent for Boron.
- ANO-RCC basis: Fixed Carbon basis set (wrong contraction coefficients, see [MOLCAS ANO-RCC](http://www.molcas.org/ANO/).
- ANO-RCC basis: Modified the 3 Th h-functions by replacing them with the 3 Ac h-functions to Th.  
                 (A mistake was made in the original work when the 3 Th h-functions were made,
                  and this has never been redone. They are too diffuse, exponents
                  (0.3140887600, 0.1256355100, 0.0502542000) for Th, Z=90, compared to
                  (0.7947153600, 0.3149038200, 0.1259615200) for Ac, Z=89, and
                  (0.8411791300, 0.3310795400, 0.1324318200) for Pa, Z=91.  
                  We have selected to just replace the 3 Th h-functions with those from the Ac basis set,
                  because the Ac g-functions are quite close to the Th g-functions, closer than Ac g-functions,
                  and therefore differences in results compared to optimized Th h-functions should be minimal.)  
                  Thanks to Kirk Peterson for pointing out the Th problem on http://daltonforum.org.
- Fixed reading of ANO-RCC basis set library file.
- Bug fix for when more than 30 excitation energies requested (EIGENVALUES NOT PAIRED problem reported by Frank Jensen).
- Fixed some bugs for two byte packing of derivative and spin-orbit two-electron integrals.
- Fixed .NEWTRA integral transformation for 32 bit integers and exactly n\*256 orbitals and no integer overflow test
  (the first 32 bits of (n\*256)**4 are zero !!!).
- Improved performance of .NEWTRA integral transformation for response calculations.
- Do not include floating orbitals in calculation of smallest atom-atom distance.
- Enable Tamm-Dancoff approximation (.TDA) for embedding models, e.g. PE, PCM etc.
- Provide date and time stamp also for Darwin (i.e. MacOSX).
- Assume nobody uses gfortran version 4.0.2 any more (removed special test for that).


## [2015.1] (2015-07-20)

### Common
- Added ANO-RCC basis set.

### DALTON
- Fixed a bug in an LRESC correction. 
- Improved calculation of one LRESC correction.
- Update PElib (v.1.2.3): Workaround for faulty system detection using Macports CMake
- Fixed a bug with Intel Compiler 15 during initialization of Cauchy-Schwarz parameters
- Fixed a bug for parallel build on some systems
- Fixed a segmentation fault for approx. 3000 basis functions
  (an array was allocated in stack memory, and became too big for default size of stack memory).
- Fixed a bug in CC sumrules.
- Fixed a bug in the preoptimization, i.e. when using smaller basis sets first to geometry optimize molecule.
- Fixed some far from optimal defaults for preoptimization.
- Fixed geometry optimization for HS-ROHF and HS-RODFT with symmetry - .SINGLY input option under '\*SCF INPUT'
  (use numerical gradients as analytical gradients are only implemented without symmetry).
- Some minor corrections to the Dalton manual.
- More reasonable output for TPCD.

### LSDALTON
- Fixed .UNCONT for EMSLs Dalton basis set format.


## [2015.0] (2015-02-18)

### Common
- Read EMSL Dalton format
- New included basis sets: pc-seg by Frank Jensen
- Many small improvements here and there

### DALTON
- CPP (Complex Polarization Propagator): MChD, NSCD
- Extensions of the polarizable embedding model (QM/MM via PE library):
  - PE-MCSCF wave function and linear response
  - PE-CPP damped linear response
  - PE for magnetic linear response using London atomic orbitals (LAOs)
- PCM-SOPPA excitation energies
- QFIT (electrostatic potential fitted charges)
- QM/CMM approach
- DFT-D3 and DFT-D3(BJ) dispersion energy corrections

### LSDALTON
- Geometry optimizations:
  - Quasi-Newton transition state optimization
  - HOPE algorithm
- Automated Counterpoise corrected DFT, HF, DEC-MP2 and CCSD interaction energies
- Dynamics: Nose-Hoover thermostat
- DFT-D3 and DFT-D3(BJ) dispersion energy corrections
- Performance improvements:
  - Charge-constrained ADMM exchange (energy+gradients)
  - Optimized Complex Polarization Propagator (CPP) solver for LSresponse
- Quadratic Response calculation to compute full dipole moment matrices in LSDalton


## [2013.4] (2014-07-10)

### DALTON
- Memory bugfix for serial PCM calculations (segmentation fault for large PCM cavities).

### LSDALTON
- Fixed a bug in the basis set reading. This bugfix affects almost no basis sets,
  and none of the standard basis sets, but a very few general contracted basis sets
  where the first contracted function had much smaller number of
  primitives compared to the last: Basis sets such as the pcS-1 basis set.
- Reduced the memory requirements for internal MPI buffer handling.


## [2013.3] (2014-06-11)

### Common
- aug-cc-pVTZ-lresc basis set added to $BASDIR.

### DALTON
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

### LSDALTON
- Environment variable LSDALTON_LAUNCHER introduced.


## [2013.2] (2014-03-05)

### Common
- Recognize CYGWIN as a LINUX and UNIX system, for proper definition of compilation flags.
- Define M_PI in C-code if not already defined (problem seen with Cygwin).
- Added setup option --blacs to be used in combination with --scalapack; defaults to --blacs=intelmpi.

### DALTON
- Fixed a bug in printing results in CPP-QRF.
- New CPP solver works also for non-direct calculation.
- More efficient evaluation of numerical Hessian when C1 symmetry
  (in each geometry step start wave function optimization from a
  converged wave function from a neighboring geometry rather than from scratch each time).
- Fix of error which sometimes caused a geometry optimization to stop with "\*\*\* ERROR, Wrong interval in WLKBIS".
- Fix of a bug which occasionally caused DALTON to abort a .STEX calculation.
- Print final geometry in xyz format (angstrom). File called "final_geometry.xyz" is put into the restart tarball.
- Append PID to scratch directory to avoid multiple tests running in the same directory.
- Improved manual for two-photon and non-adiabatic coupling.
- Updated/corrected g-factors for Ag, Nd, and Tl (thanks to M. Jaszunski).

### LSDALTON
- Print sensible error message when running out of memory.
- Added functionality to search through several basis-set libraries.
- Increased max length of WRKDIR from 60 to 200.
- Fixed a bug related to improper shutdown of MPI calculation. In the case
  of wrong LSDALTON.INP for instance the calculation will issue a error
  statement and afterward hang forever in a MPI call.
- Fixed an OpenMP bug in the calculation of how much memory there should be used during
  an exchange-correlation calculation - resulting in huge memory usage for large molecular system.


## [2013.1] (2013-12-19)

### DALTON
- Correct the printout of relativistic corrections to the shielding (thanks to M. Jaszunski).
- Compilation fix for DALTON/abacus/rma_windows.F90 (Intel 10.0.011).
- Fix of error where basis set names were changed to upper case and could not be found (reported by Yurij Rusakov).
- Each MPI slave sleeps 10 millisecond between tests for new task
  (only Intel; should enable turbomode in sequential parts of DALTON, and more efficient use of threaded MKL when combined with MPI).
- added metric scaled output of orbital response vectors in \*\*RESPONS
  (for easier interpretation of excitation operators).

### LSDALTON
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


## [2013.0] (2013-11-11)

### DALTON
- Subsystems CC using Cholesky decomposition
- Multiscale modeling using the Polarizable Embedding (PE) library.
- Static exchange (STEX) for X-ray spectroscopy
- Damped response via Complex Polarization Propagator (CPP)
- Quadratic response for open-shell DFT
- Relativistic corrections to nuclear shielding constants
- Empirical dispersion corrections DFT-D2, DFT-D3 and DFT-D3BJ
- Various performance improvements and a few bug fixes

### LSDALTON
- Dynamics using HF and DFT
- More response properties, including some magnetic properties like MCD
- DEC-MP2 energy, density and gradient
- Local orbitals
- Improved SCF optimization routines
- Massively parallel CCSD 
- MPI parallelization of HF and DFT 
- Improved integral code including MPI parallelism
- Matrix operation parallelization using PBLAS/SCALAPACK
- ADMM exchange (energy, gradients)

