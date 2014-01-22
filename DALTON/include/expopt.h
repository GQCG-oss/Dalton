! File : expopt.h
! needs: maxorb.h for MXPRIM
!
!     info for exponent optimization in basis functions (.EXPGRA)
!
      REAL*8  ALPGRD ! "alpha gradient", d E / d alpha_i
      LOGICAL EXPGRA
      COMMON /EXPOPT/ ALPGRD(MXPRIM),
     &                EXPGRA
! -- end of expopt.h --
