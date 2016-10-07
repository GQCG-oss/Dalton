! FILE: exeinf.h
!
!     ITRLVL_LAST : last integral transformation level (updated in sirtra.F)
!     LVLDRC_LAST : last Dirac integral transformation level (updated in sirtra.F)
!     FTRCTL : true - force new integral transformation because AOs have changed (typically new geometry)
!     NEWCMO : true - new integral transformation needed because CMO has changed
!     *XP    : variables used to control aba2tex.F(TWOEXP)
!     *TRONV : variables used to control abatro.F(TROINV)
!     (FT* true is abbreviation for FIRST call in a series)
!
      INTEGER         NASTXP, ITRLVL_LAST, LVLDRC_LAST
      LOGICAL         FTRONV, GTRONV, HTRONV, RTRONV,
     &                FTWOXP, SORTXP, PTRTXP, NPVTXP,
     &                FTRCTL, NEWCMO
      COMMON /EXEINF/ NASTXP, ITRLVL_LAST, LVLDRC_LAST,
     &                FTRONV, GTRONV, HTRONV, RTRONV,
     &                FTWOXP, SORTXP, PTRTXP, NPVTXP,
     &                FTRCTL, NEWCMO
! end of exeinf.h
