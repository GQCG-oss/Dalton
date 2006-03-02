C File: infpop.h
C hjaaj sep 2005: make POPNAB etc. dynamic -- maybe remove completely infpop.h ?
      PARAMETER ( MAXTYP = 50 )
      COMMON /INFPOP/ TITLEM(24),DIPOL(3),QPOL(6),QQFAC(3),
     &                IFXYZ(3),JFXYZ(3),NAME(MXCENT),POPNUC(MXCENT),
     &                POPTYP(MAXTYP,MXCENT), POPTY(MAXTYP,MXCENT)
