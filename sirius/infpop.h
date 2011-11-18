C File: sirius/infpop.h
C
C Internal information for population analysis in sirius/sirpop.F  /hjaaj
C
C IPOPTYP: select type of population analysis
C
C hjaaj sep 2005: make POPNAB etc. dynamic -- maybe remove completely infpop.h ?
      PARAMETER ( MAXTYP = 50 )
      ! spherical: s-1, sp-4, spd-9, spdf-16, spdfg-25, spdfgh-36, spdfghi-48
      ! cartesian: s-1, sp-4, spd-10, spdf-20, spdfg-35
      COMMON /INFPOP/ DIPOL(3),QPOL(6),QQFAC(3),
     &                POPNUC(MXCENT),POPTYP(MAXTYP,MXCENT),
     &                POPTY(MAXTYP,MXCENT),
     &                IFXYZ(3),JFXYZ(3),IPOPTYP
      CHARACTER*14 POPANA_TYPE
      COMMON /INFPOPC/ POPANA_TYPE
C --- end of sirius/infpop.h ---
