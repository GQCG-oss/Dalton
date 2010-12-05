#if defined (SYS_CRAY)
      REAL DUMMY
#else
      DOUBLE PRECISION DUMMY
#endif
      INTEGER IDUMMY
C     make DUMMY write protected (on most computers)
      PARAMETER ( DUMMY = 1.0D20 , IDUMMY = - 9 999 999 )
