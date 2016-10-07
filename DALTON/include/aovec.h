!
!     File: aovec.h
!
!     MXAOVC = maximum number of AO input blocks for a given center
!              (in "MOLECULE.INP" input file)
!
!     IF you change MXAOVC, then rebuild with "make".
!
      INTEGER MXAOVC, MXCONT
      PARAMETER (MXAOVC = 36, MXCONT = 36)
!     hjaaj: MXCONT needs to be 35 for ANO-DK3 uncontracted for DKH Hamiltonian
!     hjaaj: MXAOVC, MXCONT needs to be 36 for some uncontracted Dyall basis sets
!  -- end of aovec.h ---
