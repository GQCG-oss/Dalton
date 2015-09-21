!
!  FILE: qrinf.h
!
!  Special info for quadratic response
!
      INTEGER         MZVAR,  MZCONF, MZWOPT,                           &
     &                MZYVAR, MZYCON, MZYWOP,                           &
     &                MSYMA, MSYMB, MSYMC, MZYVMX
      LOGICAL         QRREST
      COMMON /QRINF/  MZVAR(8),  MZCONF(8), MZWOPT(8),                  &
     &                MZYVAR(8), MZYCON(8), MZYWOP(8),                  &
     &                MSYMA, MSYMB, MSYMC,  MZYVMX,                     &
     &                QRREST
! -- end of qrinf.h --

