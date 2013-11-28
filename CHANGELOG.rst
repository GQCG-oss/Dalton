

2013.1
------

- Correct the printout of relativistic corrections to the shielding (thanks to M. Jaszunski)
- Compilation fix for DALTON/abacus/rma_windows.F90 (Intel 10.0.011)

- Each MPI slave sleeps 1 second between tests for new task
  (should enable turbomode in sequential parts of DALTON, and
   more efficient use of threaded MKL when combined with MPI)
- added metric scaled output of orbital response vectors in **RESPONS
  (for easier interpretation of excitation operators)
