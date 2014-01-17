! File: dfterg.h
!
!     ESRDFTY : effective energy correction for SR-DFT
!     ESRDFT  : srDFT total exchange and correlation energy
!               (both ESRDFTY and ESRDFT are calculated in SIRFCK)
!
      real(8) :: EDFTX, EDFTC, EDFTY, EDFTK, ESRDFTY, ESRDFT 
      real(8) :: WDFTX, WDFTC, WDFTL, WDFTB, WDFTK
      COMMON /DFTERG/ EDFTX, EDFTC, EDFTY, EDFTK, ESRDFTY, ESRDFT,      &
     &                WDFTX, WDFTC, WDFTL, WDFTB, WDFTK
! -- end of dfterg.h --
