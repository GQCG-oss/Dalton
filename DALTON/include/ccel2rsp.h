      INTEGER MAXEL2LBL
      PARAMETER ( MAXEL2LBL = 20 )

      LOGICAL LEL2OPN
      INTEGER NEL2LBL, ISYOFEL2

      INTEGER ISTEL2(MAXEL2LBL)

      INTEGER ISYSEL2(MAXEL2LBL)
      INTEGER ISYOEL2(MAXEL2LBL,2)

      LOGICAL LPREL2(MAXEL2LBL)
      LOGICAL LORXEL2(MAXEL2LBL,2)

      CHARACTER*8 LBLEL2(MAXEL2LBL,2)

#if defined (SYS_CRAY)
      REAL EIGEL2(MAXEL2LBL)
      REAL FRQEL2(MAXEL2LBL,2)
#else
      DOUBLE PRECISION EIGEL2(MAXEL2LBL)
      DOUBLE PRECISION FRQEL2(MAXEL2LBL,2)
#endif

      COMMON/IEL2RSP/ ISTEL2, ISYSEL2, ISYOEL2, LPREL2, LORXEL2,
     &                NEL2LBL, ISYOFEL2(8), LEL2OPN
      COMMON/CEL2RSP/ LBLEL2
      COMMON/REL2RSP/ EIGEL2, FRQEL2
