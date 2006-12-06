C
C     choose reasonably large MAXBLLEN- so that loop unrolling gives
C     speedup but small compared with a cache size max block length
C     MXBLLEN.
      INTEGER MXBLLEN
      PARAMETER (MXBLLEN=100)
C
      CHARACTER*6 DFTTYP
      REAL*8          HFXFAC, HFXATT, HFXMU, DFTHR0, DFTHRL, DFTHRI, 
     &                DFTELS, RADINT, WDFTMP
      INTEGER ANGINT, ANGMIN, LEBMIN
      LOGICAL DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC, DFTHES,
     &        DFTHRS, NOPRUN, DOVWN3, DFTEST, DOVWNI, DFTIMG
      COMMON /DFTCOM/ HFXFAC, HFXATT, HFXMU, DFTHR0, DFTHRL, DFTHRI, 
     &                DFTELS, RADINT, WDFTMP,
     &        ANGINT, ANGMIN, LEBMIN, 
     &        DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC, DFTHES, 
     &        DFTHRS, NOPRUN, DOVWN3, DFTEST, DOVWNI, DFTIMG
      COMMON /DFTCHR/ DFTTYP 


