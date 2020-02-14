# DALTON Change Log -- All notable changes to the DALTON project will be documented in this file.

## [2020.alpha] (Dalton2020 alpha)

### New features added
- Added pcH-n and aug-pcH-n basis sets levels 1-4 to Dalton basis set library.
- Extended "super matrix integrals" to work with >255 basis functions. (H. J. Aa. Jensen)
- Calculate and print oscillator strengths based on .EEF dipole transition moments. (H. J. Aa. Jensen)
- Added ".TDA TR" keyword for invoking Tamm-Dancoff approximation for triplet response properties under \*\*PROPERTIES.
  Useful for avoiding (near-)triplet-instability problems, for example in DFT calculations of spin-spin coupling constants. (H. J. Aa. Jensen)
- Added ".TDA SI" keyword for invoking Tamm-Dancoff approximation for singlet response properties under \*\*PROPERTIES. (H. J. Aa. Jensen)
- Added possibility to optimize MCSCF singlet wave functions with CSFs when used for triplet properties,
  both in \*\*RESPONS and \*\*PROPERTIES. Previously .DETERMINANTS in wave function optimization was required. (H. J. Aa. Jensen)
- Added the ability in QFITLIB to fit up to and including quadrupoles (C. Steinmann)
- Added the CVS approximation for CC calculations of core-excited states (S. Coriani et al.)
- Added the possibility to calculate triplet-triplet excited state moments using the EOM-CC approximation (R. Faber)
- New features available through the Polarizable Embedding library (PElib)
  - Polarizable density embedding (PDE) model (use -DENABLE\_PDE=ON during setup to enable it [requires HDF5])
    - J. M. H. Olsen, C. Steinmann, K. Ruud, and J. Kongsted, J. Phys. Chem. A 119, 5344 (2015)
    - P. Reinholdt, J. Kongsted, and J. M. H. Olsen, J. Phys. Chem. Lett. 8, 5949 (2017)
  - PDE-CC2, PDE-CCSD, and PDE-CCSDR(3) including linear and quadratic response (also enables PE-CC through PElib)
    - D. Hrsak, J. M. H. Olsen, and J. Kongsted, J. Chem. Theory Comput. 14, 1351 (2018)
  - FixSol continuum solvation with FIXPVA2 cavity tesselation
    - M. S. Nørby, C. Steinmann, J. M. H. Olsen, H. Li, and J. Kongsted, J. Chem. Theory Comput. 12, 5050 (2016)
    - N. M. Thellamurege and H. Li, J. Chem. Phys. 137, 246101 (2012)
  - Enabled cubic response for PE-HF/DFT and PDE-HF/DFT
    - J. M. H. Olsen and J. Kongsted, Adv. Quantum Chem. 61, 107 (2011)
  - Effective external field (EEF) can now be used for all dipole properties
  - Added support for AMOEBA potential
- Added the ability to run SOPPA linear response calculations via the AOSOPPA module (R. Faber et al.)

### Added
- Added information about .MS2 input option to manual, quit if invalid value specified. (H. J. Aa. Jensen)

### Fixed
- Compilation errors with gfortran-9.2.0 and 64-bit integers
- Errors when running DFT with 64-bit integers and MPI. (P. Reinholdt and H. J. Aa. Jensen)
- Errors for Fermi-contact (FC) labels on APROPER and therefore FC properties in \*\*RESPONS when more than 99 atoms (H. J. Aa. Jensen)
- Error for \*ESR spin-dipole properties when more than 33 atoms (H. J. Aa. Jensen)
- Singlet totally-symmetric excitation energies for MCSCF in \*\*RESPONS with super-symmetry activated (.SUPSYM keyword). (H. J. Aa. Jensen)
- Make sure we include all (near-)degenerate diagonal elements for linear response excitation energies
  via \*\*PROPERTIES .EXCITA or via \*\*RESPONS \*LINEAR .SINGLE (increase .NROOTS if needed).
  Otherwise the calculation will probably exhibit spin and/or space symmetry contamination proportional
  to the convergence threshold. (H. J. Aa. Jensen)
- Never use plus combinations of determinants as start guess for singlet linear response excitation energies
  when reference wave function is not singlet (we do not want singlet states then). (H. J. Aa. Jensen)
- Dalton script: fix for using input files located in subfolders
- Fixed error from March 2015 which meant that double-hybrid DFT was not working correctly (MP2 part was ignored).
- Fixed error for MC-TDA excitation energies for RASSCF (CASSCF was OK).

### Fixes in enclosed basis set files
- Error in diffuse d-orbital exponents for Aluminum and Silicon (factor 10 too big) in aug-cc-pV(D+d)Z basis sets (H. J. Aa. Jensen)
- Error in diffuse f-orbital exponents for Aluminum and Silicon (factor 10 too big) in aug-cc-pV(Q+d)Z basis sets (H. J. Aa. Jensen)
- Error in diffuse f-orbital exponent for Aluminum (factor 10 too big) in aug-cc-pV(T+d)Z basis sets (H. J. Aa. Jensen)

### Changed
- Allow basis set(s) after BASIS in line 1 of .mol file (instead of on second line). (H. J. Aa. Jensen)


## [2018.2] (2019-03-17)

### Fixed
- Fixed error in AO-direct CC3 response calculations causing segmentation faults
- Fixed calculation of DSO contribution to spin-spin coupling for MCSCF and HSROHF when no symmetry
- Dalton script: do not set OMP\_NUM\_THREADS=1 if not MPI parallel (better performance
  for sequential calculations if threaded blas is used, e.g. MKL or openBLAS)
- Dalton script: stop if user asks for MPI run with a sequential dalton.x
- More robust .STEX input specification (old failed in some situations with gfortran 8); changed documentation accordingly

### Added
- Dalton script: -gb and -ngb options for specifying work memory in gigabytes
- Dalton script: -np as an alternative to -N for specifying number of MPI nodes

## [2018.1] (2019-01-14)

### Fixed
- Error in code for 2-el integral transformation level 4 (used in some cases for MCSCF). Error was not in Dalton2016.
- Error in export for FDE fixed.
- Compilation with PGI compilers now possible (but code with pelib or qfit using gen1int is not working).

## [2018.0] (2018-11-19)

### New features added
- New parallel 2-electron integral transformation .FCKTRA  (H. J. Aa. Jensen).
- Triplet excitation energies and polarizabilities with the AO-SOPPA code (P. A. B. Haase).
- Analytical PE-HF/DFT molecular gradients (see e.g. test/pehf\_geoopt for example of input).
  - Reference: N. H. List, M. T. B. Beerepoot, J. M. H. Olsen, B. Gao, K. Ruud, H. J. Aa. Jensen, and J. Kongsted. J. Chem. Phys. 142, 034119 (2015).
- Freezing atoms in geometry optimization (N. H. List and H. J. Aa. Jensen).
  (See e.g. test/geoopt\_freeze for example of input.)
- Effective external field (EEF) for one- and two-photon absorption in PE-HF/DFT calculations.
  - Reference: N. H. List, H. J. Aa. Jensen, and J. Kongsted. Phys. Chem. Chem. Phys. 18, 10070 (2016).
- Remove the most diffuse virtual orbitals after SCF or MCSCF (new .VIRTRUNC option).
- Add purely classical multipole-multipole interaction energy in PE-QM calculations (which can be skipped using the .SKIPMUL keyword under the \*PEQM section).
- Add basic frozen density embedding (FDE) functionality (A. Gomes, C. Jacob, L. Visscher).
- Dipole velocity complex linear polarizability with test rsp\_cpp\_veloci (N. H. List).
- Resonant-convergent (damped) cubic response at HF/DFT levels (T. Fahleson and P. Norman)
  - Reference: T. Fahleson and P. Norman. J. Chem. Phys. 147, 144109 (2017).

### Fixed
- Nuclear model keyword .NUCMOD was ignored, now the Gaussian nuclear model used in the Dirac program can be used in Dalton.
- Open-shell DFT is not implemented for many derivative properties in \*\*PROPERTIES, dalton now quits.
- Bugfix for .MNF\_SO (mean-field spin-orbit, AMFI) when basis set has big exponents (>10^9).
- Open-shell doublet ROKS DFT geometry optimization.
- Fix of .GSPOL for parallel PE-QM quadratic response.
- Corrected text about elimination of some two-photon transitions between excited states
  because they were duplicates (text had "Third order" instead of "Second order").
- Bugfixes for CC2 in environments (incl. external fields) by Ove Christiansen.
- Bugfix for .TDA for MCSCF (i.e. zero B matrix in linear response in \*\*RESPONSE).
- Bugfix for .GASCI after .HSROHF
- Fix of several library basis sets that were not read correctly for some atoms, which caused Dalton to abort.
- Now a .G-TENSOR calculation in \*\*RESPONSE:\*ESR module with a CI wave function does not abort.

### Changed
- OK to run ECD or OECD with SOPPA.
- More documentation of .STEX in manual.
- Default induced-dipole solver in polarizable embedding (through .PEQM keyword) is changed to JI/DIIS method, which improves parallel scaling performance, and default convergence threshold for induced dipoles is changed $`1.0\cdot10^{-8}>|\mu^{[k]}-\mu^{[k-1]}|`$ where $`\mu`$ is a vector containing all induced dipoles and $`k`$ is the iteration index.
- Minimum CMake version is now v3.1.

### Deprecated
- Environment variable DALTON\_NUM\_MPI\_PROCS is deprecated and will be removed in future releases, use DALTON\_LAUNCHER instead


## [2016.2] (2016-07-12)

### Added
- Added and documented Basis=INTGRL option for ATOMBASIS in .mol file.
- Included Be in cc-pV5Z basis set

### Fixed
- More robust code for reading exponents and contraction coefficients in Dalton-type basis set files, incl. such files from EMSL
- Work-around for Intel 15 compiler I/O problem in some response calculations
- Fix for spin-orbit coupling (SOC) between S/T excited states of same symmetry (problem reported on daltonforum.org)
- Further fixes of MCSCF in \*\*PROPERTIES for more than 255 basis functions - hopefully it is OK now for all requests.
- Fixed an error in the manual for spin-dipole (problem reported on daltonforum.org)
- Fix of open-shell Hartree-Fock occupation output (only output, not the calculation, was wrong if ROHF was followed by MCSCF)
- Fix of Douglas-Kroll post-SCF with less than 256 contracted basis functions, but more than 255 uncontracted basis functions
- Fix of an insufficient memory error for construction of 2-el. integrals in Dirac format with more than 255 basis functions
- Removed OpenACC CMake variable (currently no OpenACC directives in Dalton).
- Fix of Douglas-Kroll post-SCF with less than 256 contracted basis functions, but more than 255 uncontracted basis functions

## [2016.1] (2016-04-07)

### Added
- Possibility to read basis set files as made by the EMSL web site
  (this makes it possible to also read basis set files in basis/ based on
  emsl output as e.g. aug-pcseg-1; only LSDALTON has been able to read them so far)

### Fixed
- MCSCF in \*\*PROPERTIES for more than 255 basis functions (fixes problem with MCSCF shielding reported on daltonforum.org)
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
  (the first 32 bits of (n\*256)\*\*4 are zero !!!).
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
- Environment variable DALTON\_USE\_GLOBAL\_SCRATCH disables copying of binaries to worker nodes.
- Environment variable DALTON\_LAUNCHER introduced.
- Fixed output information about number of MPI processes and number of OpenMP threads.
- Added information in the error messages when values in maxorb.h are exceeded (which values to increase).
- Increased some of the values in the common blocks:
  MXSHEL 1000 -> 1500; MXCORB 2400 -> 5000; MXPRIM 8000 -> 15000;
  MAXOCC 800 -> 1500; MXCENT 200 -> 500; MXCENT\_QM 200 -> 500
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
- Environment variable LSDALTON\_LAUNCHER introduced.


## [2013.2] (2014-03-05)

### Common
- Recognize CYGWIN as a LINUX and UNIX system, for proper definition of compilation flags.
- Define M\_PI in C-code if not already defined (problem seen with Cygwin).
- Added setup option --blacs to be used in combination with --scalapack; defaults to --blacs=intelmpi.

### DALTON
- Fixed a bug in printing results in CPP-QRF.
- New CPP solver works also for non-direct calculation.
- More efficient evaluation of numerical Hessian when C1 symmetry
  (in each geometry step start wave function optimization from a
  converged wave function from a neighboring geometry rather than from scratch each time).
- Fix of error which sometimes caused a geometry optimization to stop with "\*\*\* ERROR, Wrong interval in WLKBIS".
- Fix of a bug which occasionally caused DALTON to abort a .STEX calculation.
- Print final geometry in xyz format (angstrom). File called "final\_geometry.xyz" is put into the restart tarball.
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
- Compilation fix for DALTON/abacus/rma\_windows.F90 (Intel 10.0.011).
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

