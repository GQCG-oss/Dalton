C
C     File: aovec.h
C
C     MXAOVC = maximum number of AO input blocks for a given center
C              (in "MOLECULE.INP" input file)
C
C     IF you change MXAOVC you should do a "make depend"
C     and then rebuild the program using the command "make".
C
      INTEGER MXAOVC, MXCONT
      PARAMETER (MXAOVC = 25, MXCONT = 35)
C     hjaaj: MXCONT needs to be 35 for ANO-DK3 uncontracted for DKH Hamiltonian
