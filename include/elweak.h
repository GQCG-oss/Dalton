C$Id: elweak.h,v 1.3 2002-01-17 11:48:57 alfch Exp $
C
C
C     Constants for electroweak interaction
C
      LOGICAL BGWEIL, PVFINL, ELWEAK, PVPSO, PVSO
      PARAMETER (WEINBG = 1D0 - 4D0 * 0.2319D0)
      PARAMETER (GFERMI = 5.73416D-17)
      COMMON /COMEPV/ ISOTPV(MXCENT), BGWEIN, PVFINN, BGWEIL, PVFINL, 
     &       ELWEAK
