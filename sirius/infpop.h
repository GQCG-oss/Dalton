C hjaaj sep 2005: make POPNAB etc. dynamic -- maybe remove completely infpop.h ?
      PARAMETER ( MAXTYP = 20 )
      COMMON /INFPOP/ TITLEM(24),DIPOL(3),QPOL(6),QQFAC(3),
     *                IFXYZ(3),JFXYZ(3),NAME(MXCENT),POPNUC(MXCENT),
     *                POPTYP(20,MXCENT),POPTY(20,MXCENT)
