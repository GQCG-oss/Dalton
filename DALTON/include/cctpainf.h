      INTEGER MXSMOP, MXSMSEL
      PARAMETER ( MXSMOP = 125 , MXSMSEL = 50 )

      LOGICAL SELSMST, HALFFR, LTPA_USE_X2, LTPA_USE_O2
      LOGICAL TPOLDW
      INTEGER ISMSEL, NSMSEL, NSMOPER
      INTEGER IASMOP, IBSMOP
      INTEGER IPRSM

#if defined (SYS_CRAY)
      REAL BSMFR
#else
      DOUBLE PRECISION BSMFR
#endif

      COMMON /CCINFTPA/  BSMFR(MXSMSEL), ISMSEL(MXSMSEL,2),
     *                 IASMOP(MXSMOP), IBSMOP(MXSMOP),
     *                 NSMSEL, IPRSM, NSMOPER,
     *                 SELSMST, HALFFR, LTPA_USE_X2, LTPA_USE_O2,
     *                 TPOLDW

