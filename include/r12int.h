      PARAMETER (MXMULB = 2)
C
C     The next parameters indicates the number of integer and
C     logical variables to be transferred in a parallel calculation.
C     Must be updated if COMR12 is changed. Note: Order of variables in
C     COMR12 must not be changed.
C
      PARAMETER (NR12I = 7, NR12L = 30)
      LOGICAL R12INT, R12TRA, U12INT, U21INT, R12SQR,
     *        R12CAL, R12NOA, R12NOP, R12NOB, R12OLD, NORXR,
     *        R12HYB, COMBSS, R12EIN, R12EOR, R12ECO, R12XXL,
     *        V12INT, NOTR12, NOTU12, NOTV12,
     *        ANTICO, DIRSCF, U12DIR, R12DIR, 
     *        V12DIR, STEPIV, DCCR12, CUSP12, DUMR12, CCPAIR
      LOGICAL LMULBS, LAUXBS
      INTEGER NMULBS, MBAS1(8), MBAS2(8), 
     *                NORB1(8), NORB2(8), MBSMAX
      CHARACTER*80 MULNAM(MXMULB)
      CHARACTER*8 LABEL
      COMMON /COMR12/ GAMAC, GAMAD, R12INT, R12TRA, U12INT, 
     *        U21INT, R12SQR, R12CAL, R12NOA, R12NOP, R12NOB, 
     *        R12OLD, NORXR, R12HYB, COMBSS, R12EIN,
     *        R12EOR, R12ECO, R12XXL, V12INT, NOTR12, NOTU12, 
     *        NOTV12, ANTICO, DIRSCF, U12DIR, R12DIR, V12DIR, 
     *        STEPIV, DCCR12, CUSP12, DUMR12,
     *        INTGAC, INTGAD, IADV12, IADR12, IADU12, IADU21, NOPP12,
     *        CCPAIR
      REAL*8 VCLTHR, SVDTHR
      COMMON /COMVCL/ VCLTHR, SVDTHR
      COMMON /CMMBAS/ LMULBS, NMULBS
      COMMON /CMMMUL/ MULNAM, NPLSH, MBAS1, MBAS2, LAUXBS,
     *                NORB1, NORB2, MBSMAX
      COMMON /CMMOLL/ LABEL
