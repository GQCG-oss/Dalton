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
     &        XMULFAC_READIN, DSFAC,  HEAVISIDE_PVALUE
      INTEGER IPRDFT, ANGINT, ANGMIN, LEBMIN, IWINT, NTOT_DFTGRID
      LOGICAL DFTADD, DFTRUN, DFTPOT, DFTORD, DFTASC, DFTHES,           &
     &        DFTHRS, NOPRUN, DOVWN3, DFTEST, DOVWNI, DFTIMG, HFXSET,   &
     &        DODFTD, DFTD_TEST,                                        &
     &        DFTGRID_DONE,   DFTGRID_DONE_OLD    
              ! DFTGRID_DONE for REAQUA, DFTGRID_DONE_OLD for REAQUA_OLD ("grid done")
!     variables for srDFT /hjaaj
      LOGICAL SRDFTRUN ! variable primarily for abacus/ routines
      LOGICAL DOSRX_LDA, DOSRX_GGA, DOSRBCK, DOHFEXCH, DOSRX_WIB,       &
     &        DOSRC_LDA, DOSRC_GGA, DOSRC_MULOCAL,                      &
     &        DOSRGGA2, DOSRLYPT, SRCMULOFAC, DSLOCALFAC,               &
     &        DOSRC_WIB, ISJT, DOSRX_PBEHSE, DOSRX_PBETCS, DOSRC_PBETCS,&
     &        DOSRC_PBETCSJ, DOSRC_PBERI, DOSRC_PBEWI, DOSRX_PBERI,     &
     &        DOSRX_PBEGWS, DOSRC_PBEGWS, DOSRX_LDA_S, DOSRC_LDA_S,     &
     &        DOSRC_LDA_PW92, DOSRC_LDA_PW92_S,                         &
     &        DOSRC_TPSS_S, DOSRC_TPSS, DOSRX_TPSS_S, DOSRX_TPSS,       &
     &        DOSRX_PBEGWS_S, DOSRC_PBEGWS_S, DOSRC_PBEGWS_VWN5,        &
     &        DOSRC_MULOC_GGA, DOSRC_MULOD_GGA, DOSRC_MULOE_GGA,        &
     &        DFT_SPINDNS,  DFT_LOCALSPIN, DOSRC_MD_LDA, DOSRC_PBELO,   &
     &        DOLAX_LDAS,DOLANSC_LDAS,DOLANSC_LDA,DOLAX_LDA,DOLASC_LDA, &
     &        DOLAX_PBEGWS,DOLANSC_PBEGWS,DOLASC_PBEGWS,DOLAX_GGABCK,   &
     &        DOLANC_GGALYP,DOLASC_GGALYP,DOSRC_LYPRI,                  &
     &        DOC_PBE_nomu_S, DOC_PBE_nomu, DOSRX_wPBE_S, DOSRX_wPBE,   &
     &        DOSRC_VWN5_S, DOSRC_VWN5, DOC_VWN5_nomu_S, DOC_VWN5_nomu, &
     &        SRLDA, SRPBEGWS, SRPBERI, SRPBE0GWS
      COMMON /DFTCOM/ HFXFAC, HFXATT, HFXMU,                            &
     &        DFTHR0, DFTHRL, DFTHRI, DFTELS, RADINT, WDFTMP, COPFAC,   &
     &        XMULFAC_READIN, DSFAC, HEAVISIDE_PVALUE,                  &
! integer:
     &        IPRDFT, ANGINT, ANGMIN, LEBMIN, IWINT, NTOT_DFTGRID,      &
! logical:
     &        DFTADD, DFTRUN, DFTPOT, DFTORD, DFTASC, DFTHES,           &
     &        DFTHRS, NOPRUN, DOVWN3, DFTEST, DOVWNI, DFTIMG, HFXSET,   &
     &        DODFTD, DFTD_TEST,                                        &
     &        DFTGRID_DONE,   DFTGRID_DONE_OLD,                         &
! srDFT (logical):
     &        SRDFTRUN,                                                 &
     &        DOSRX_LDA, DOSRX_GGA, DOSRBCK, DOHFEXCH, DOSRX_WIB,       &
     &        DOSRC_LDA, DOSRC_GGA, DOSRC_MULOCAL(0:3),                 &
     &        DOSRGGA2, DOSRLYPT, SRCMULOFAC, DSLOCALFAC,               &
     &        DOSRC_WIB, ISJT, DOSRX_PBEHSE, DOSRX_PBETCS, DOSRC_PBETCS,&
     &        DOSRC_PBETCSJ, DOSRC_PBERI, DOSRC_PBEWI, DOSRX_PBERI,     &
     &        DOSRX_PBEGWS, DOSRC_PBEGWS, DOSRX_LDA_S,  DOSRC_LDA_S,    &
     &        DOSRC_LDA_PW92, DOSRC_LDA_PW92_S,                         &
     &        DOSRC_TPSS_S, DOSRC_TPSS, DOSRX_TPSS_S, DOSRX_TPSS,       &
     &        DOSRX_PBEGWS_S,  DOSRC_PBEGWS_S, DOSRC_PBEGWS_VWN5,       &
     &        DOSRC_MULOC_GGA, DOSRC_MULOD_GGA, DOSRC_MULOE_GGA,        &
     &        DFT_SPINDNS,  DFT_LOCALSPIN, DOSRC_MD_LDA, DOSRC_PBELO,   &
     &        DOLAX_LDAS,DOLANSC_LDAS,DOLANSC_LDA,DOLAX_LDA,DOLASC_LDA, &
     &        DOLAX_PBEGWS,DOLANSC_PBEGWS,DOLASC_PBEGWS,DOLAX_GGABCK,   &
     &        DOLANC_GGALYP,DOLASC_GGALYP,DOSRC_LYPRI,                  &
     &        DOC_PBE_nomu_S, DOC_PBE_nomu, DOSRX_wPBE_S, DOSRX_wPBE,   &
     &        DOSRC_VWN5_S, DOSRC_VWN5, DOC_VWN5_nomu_S, DOC_VWN5_nomu, &
     &        SRLDA, SRPBEGWS, SRPBERI, SRPBE0GWS
      INTEGER DFTCOMlast
      COMMON /DFTCOM/ DFTCOMlast
      !  Very important !!!  Always keep DFTCOMlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
!
      CHARACTER*6  DFTTYP
!     variables for srDFT
      CHARACTER*23 SRXFUN, SRCFUN, SRLOCALSPIN
      COMMON /DFTCHR/ DFTTYP,                                           &
     &                SRXFUN, SRCFUN, SRLOCALSPIN
! -- end of dftcom.h --
