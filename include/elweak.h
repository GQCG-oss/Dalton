C$Id: elweak.h,v 1.2 2001-10-01 13:22:18 vebjornb Exp $
C
C
C     Constants for electroweak interaction
C
      LOGICAL BGWEIL, PVFINL
      PARAMETER (WEINBG = 1D0 - 4D0 * 0.2319D0)
      PARAMETER (GFERMI = 5.73416D-17)
      COMMON /COMEPV/ ISOTPV(MXCENT), BGWEIN, PVFINN, BGWEIL, PVFINL
