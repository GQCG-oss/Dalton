! File: dfterg.h
!
!     ESRDFTY : effective energy correction for SR-DFT,
!               is calculated in SIRFCK.
      real(8) :: EDFTX, EDFTC, EDFTY, EDFTK, ESRDFTY
      real(8) :: WDFTX, WDFTC, WDFTL, WDFTB, WDFTK
      COMMON /DFTERG/ EDFTX, EDFTC, EDFTY, EDFTK, ESRDFTY,              &
     &                WDFTX, WDFTC, WDFTL, WDFTB, WDFTK
