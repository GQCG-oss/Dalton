C     MXSHEL = maximum number of shells (insert shell definition here).
C     MXPRIM = maximum number of primitives.
C     MXCORB = maximum number of orbitals (possibly contracted).
C     MAXOCC = maximum number of occupied orbitals
      INTEGER MXSHEL, MXPRIM, MXCORB, MXORBT, MAXOCC
#if defined (VAR_ABASMALL)
      PARAMETER (MXSHEL = 200, MXPRIM = 800, MXCORB = 400,
     *           MXORBT = MXCORB*(MXCORB + 1)/2)
      PARAMETER ( MAXOCC = 120 )
#else
      PARAMETER (MXSHEL = 750, MXPRIM = 4000, MXCORB = 1200,
     *           MXORBT = MXCORB*(MXCORB + 1)/2)
      PARAMETER ( MAXOCC = 400 )
#endif
