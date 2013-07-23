!
! FILE: cbilnr.h
!
! Used to control linear response solver in ABACUS
!
      LOGICAL         ALFA, ROAA, ROAG, STATIC
      integer          MXFR,        MAXLN
      PARAMETER       (MXFR  = 300, MAXLN = 80)
      CHARACTER*8     LABALN
      COMMON /LNLBL / LABALN(MAXLN)
      real(8)         THCLNR, FRVAL
      integer         NFRVAL, IPRLNR, NABALN
      COMMON /CBILNR/ THCLNR, FRVAL(MXFR), NFRVAL, IPRLNR, NABALN,      &
     &                ALFA, ROAA, ROAG, STATIC
! -- end of cbilnr.h --
