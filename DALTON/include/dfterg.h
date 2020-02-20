! File: dfterg.h
!
!     ESRDFTY   : effective energy correction for SR-DFT
!     ESRDFT(1) : srDFT total exchange and correlation energy = ESRDFT(2) + ESRDFT(3)
!     ESRDFT(2) : srDFT total exchange energy
!     ESRDFT(3) : srDFT total correlation energy
!               (both ESRDFTY and ESRDFT are calculated in SIRFCK)
!
      real(8) :: EDFTX, EDFTC, EDFTY, EDFTK, ESRDFTY, ESRDFT(3)
      real(8) :: WDFTX, WDFTC, WDFTL, WDFTB, WDFTK
      COMMON /DFTERG/ EDFTX, EDFTC, EDFTY, EDFTK, ESRDFTY, ESRDFT,      &
     &                WDFTX, WDFTC, WDFTL, WDFTB, WDFTK
! -- end of dfterg.h --
