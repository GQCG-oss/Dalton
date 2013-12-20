

2013.2 (this will be the next patch)
-----------------------------------

- Document changes here.


2013.1 (2013-12-19)
-------------------

- Correct the printout of relativistic corrections to the shielding (thanks to M. Jaszunski).
- Compilation fix for DALTON/abacus/rma_windows.F90 (Intel 10.0.011).
- Fix of error where basis set names were changed to upper case and could not be found (reported by Yurij Rusakov).
- Each MPI slave sleeps 10 millisecond between tests for new task
  (only Intel; should enable turbomode in sequential parts of DALTON, and more efficient use of threaded MKL when combined with MPI).
- added metric scaled output of orbital response vectors in \*\*RESPONS
  (for easier interpretation of excitation operators).
- Fixed a bug in Jengine, related to screening for nonsymmetric density matrices.
  This may affect CCSD and some response calculation. 
- Modified the input section of the manual concerning 
  Casida-Salahub asymptotic correction CS00 (thanks to Raul Crespo).
- Changed defaults for Casida-Salahub asymptotic correction CS00 (thanks to Raul Crespo).
- Fixed errors in the MCD B terms output files (.dat files) now one file is generated
  for each B term and each A term (thanks to Raul Crespo) .
- Modified the input section of the manual concerning MCD B terms. Added desciption of MCDEXSTATES.
