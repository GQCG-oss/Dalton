C     File: dftcom.h
C
C     choose reasonably large MXBLLEN- so that loop unrolling gives
C     speedup but small compared with a cache size max block length
C
      INTEGER MXBLLEN
      PARAMETER (MXBLLEN=100)
C
C     HFXSET: used in input to determine if HFXFAC has been specified by user with .HFXFAC
      REAL*8  HFXFAC, HFXATT, HFXMU,
     &        DFTHR0, DFTHRL, DFTHRI, DFTELS, RADINT, WDFTMP, COPFAC
      INTEGER IPRDFT, ANGINT, ANGMIN, LEBMIN, IWINT
      LOGICAL DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC, DFTHES,
     &        DFTHRS, NOPRUN, DOVWN3, DFTEST, DOVWNI, DFTIMG
C     variables for srDFT /hjaaj
      LOGICAL DOSRX_LDA, DOSRX_GGA, DOSRBCK, DOHFEXCH, DOSRX_WIB,
     &        DOSRC_LDA, DOSRC_GGA, DOSRMULO, DOSRGGA2, DOSRLYPT,
     &        DOSRC_WIB, ISJT, DOSRX_PBEHSE, DOSRX_PBETCS, DOSRC_PBETCS,
     &        DOSRC_PBETCSJ, DOSRC_PBERI, DOSRC_PBEWI, DOSRX_PBERI,
     &        HFXSET,
     &        DOSRX_PBEGWS, DOSRC_PBEGWS, DOSRX_LDAS, DOSRC_LDAS
      COMMON /DFTCOM/ HFXFAC, HFXATT, HFXMU,
     &        DFTHR0, DFTHRL, DFTHRI, DFTELS, RADINT, WDFTMP, COPFAC,
     &        ! integer:
     &        IPRDFT, ANGINT, ANGMIN, LEBMIN, IWINT,
     &        ! logical:
     &        DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC, DFTHES,
     &        DFTHRS, NOPRUN, DOVWN3, DFTEST, DOVWNI, DFTIMG,
     &        ! srDFT (logical):
     &        DOSRX_LDA, DOSRX_GGA, DOSRBCK, DOHFEXCH, DOSRX_WIB,
     &        DOSRC_LDA, DOSRC_GGA, DOSRMULO, DOSRGGA2, DOSRLYPT,
     &        DOSRC_WIB, ISJT, DOSRX_PBEHSE, DOSRX_PBETCS, DOSRC_PBETCS,
     &        DOSRC_PBETCSJ, DOSRC_PBERI, DOSRC_PBEWI, DOSRX_PBERI,
     &        HFXSET,
     &        DOSRX_PBEGWS, DOSRC_PBEGWS, DOSRX_LDAS,  DOSRC_LDAS
C
      CHARACTER*6  DFTTYP
C     variables for srDFT /JT
      CHARACTER*10 SRXFUN, SRCFUN
      COMMON /DFTCHR/ DFTTYP,
     &                SRXFUN, SRCFUN
C -- end of dftcom.h --
