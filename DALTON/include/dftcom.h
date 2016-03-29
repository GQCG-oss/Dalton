!     File: dftcom.h
!
!     choose reasonably large MXBLLEN - so that loop unrolling gives
!     speedup but small compared with a cache size max block length
!
      INTEGER MXBLLEN
      PARAMETER (MXBLLEN=100)
!
!     HFXSET: used in input to determine if HFXFAC has been specified by user with .HFXFAC
      REAL*8  HFXFAC, HFXATT, HFXMU,                                    &
     &        DFTHR0, DFTHRL, DFTHRI, DFTELS, RADINT, WDFTMP, COPFAC,   &
     &        XMULFAC_READIN, DSFAC, HEAVISIDE_PVALUE
      INTEGER IPRDFT, ANGINT, ANGMIN, LEBMIN, IWINT
      LOGICAL DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC, DFTHES,   &
     &        DFTHRS, NOPRUN, DOVWN3, DFTEST, DOVWNI, DFTIMG, HFXSET,   &
     &        DODFTD, GRDONE_OLD                                          ! GRDONE for REAQUA, GRDONE_OLD for REAQUA_OLD ("grid done")
!     variables for srDFT /hjaaj
      LOGICAL DOSRX_LDA, DOSRX_GGA, DOSRBCK, DOHFEXCH, DOSRX_WIB,       &
     &        DOSRC_LDA, DOSRC_GGA, DOSRC_MULOCAL,                      &
     &        DOSRGGA2, DOSRLYPT, SRCMULOFAC, DSLOCALFAC,               &
     &        DOSRC_WIB, ISJT, DOSRX_PBEHSE, DOSRX_PBETCS, DOSRC_PBETCS,&
     &        DOSRC_PBETCSJ, DOSRC_PBERI, DOSRC_PBEWI, DOSRX_PBERI,     &
     &        DOSRX_PBEGWS, DOSRC_PBEGWS, DOSRX_LDA_S, DOSRC_LDA_S,     &
     &        DOSRC_LDA_PW92, DOSRX_LDA_PW92,                           &
     &        DOSRX_PBEGWS_S, DOSRC_PBEGWS_S, DOSRC_PBEGWS_PW92,        &
     &        DOSRC_MULOC_GGA, DOSRC_MULOD_GGA, DOSRC_MULOE_GGA,        &
     &        DFT_SPINDNS,  DFT_LOCALSPIN, DOSRC_MD_LDA, DOSRC_PBELO,   &
     &        DOLAX_LDAS,DOLANSC_LDAS,DOLANSC_LDA,DOLAX_LDA,DOLASC_LDA, &
     &        DOLAX_PBEGWS,DOLANSC_PBEGWS,DOLASC_PBEGWS,DOLAX_GGABCK,   &
     &        DOLANC_GGALYP,DOLASC_GGALYP,DOSRC_LYPRI                   &
      COMMON /DFTCOM/ HFXFAC, HFXATT, HFXMU,                            &
     &        DFTHR0, DFTHRL, DFTHRI, DFTELS, RADINT, WDFTMP, COPFAC,   &
     &        XMULFAC_READIN, DSFAC, HEAVISIDE_PVALUE,                  &
! integer:
     &        IPRDFT, ANGINT, ANGMIN, LEBMIN, IWINT,                    &
! logical:
     &        DFTADD, GRDONE, DFTRUN, DFTPOT, DFTORD, DFTASC, DFTHES,   &
     &        DFTHRS, NOPRUN, DOVWN3, DFTEST, DOVWNI, DFTIMG, HFXSET,   &
     &        DODFTD, GRDONE_OLD,                                       &
! srDFT (logical):
     &        DOSRX_LDA, DOSRX_GGA, DOSRBCK, DOHFEXCH, DOSRX_WIB,       &
     &        DOSRC_LDA, DOSRC_GGA, DOSRC_MULOCAL(0:3),                 &
     &        DOSRGGA2, DOSRLYPT, SRCMULOFAC, DSLOCALFAC,               &
     &        DOSRC_WIB, ISJT, DOSRX_PBEHSE, DOSRX_PBETCS, DOSRC_PBETCS,&
     &        DOSRC_PBETCSJ, DOSRC_PBERI, DOSRC_PBEWI, DOSRX_PBERI,     &
     &        DOSRX_PBEGWS, DOSRC_PBEGWS, DOSRX_LDA_S,  DOSRC_LDA_S,    &
     &        DOSRC_LDA_PW92, DOSRX_LDA_PW92,                           &
     &        DOSRX_PBEGWS_S,  DOSRC_PBEGWS_S, DOSRC_PBEGWS_PW92,       &
     &        DOSRC_MULOC_GGA, DOSRC_MULOD_GGA, DOSRC_MULOE_GGA,        &
     &        DFT_SPINDNS,  DFT_LOCALSPIN, DOSRC_MD_LDA, DOSRC_PBELO,   &
     &        DOLAX_LDAS,DOLANSC_LDAS,DOLANSC_LDA,DOLAX_LDA,DOLASC_LDA, &
     &        DOLAX_PBEGWS,DOLANSC_PBEGWS,DOLASC_PBEGWS,DOLAX_GGABCK,   &
     &        DOLANC_GGALYP,DOLASC_GGALYP,DOSRC_LYPRI
!
      CHARACTER*6  DFTTYP
!     variables for srDFT
      CHARACTER*10 SRXFUN, SRCFUN, SRLOCALSPIN
      COMMON /DFTCHR/ DFTTYP,                                           &
     &                SRXFUN, SRCFUN, SRLOCALSPIN
! -- end of dftcom.h --
