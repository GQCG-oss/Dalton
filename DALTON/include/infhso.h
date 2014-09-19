! FILE : include/infhso.h
!
      INTEGER MXPHOS
      PARAMETER (MXPHOS = 10)

      REAL*8  PHOSMAT, PHOSMN
      REAL*8  PHOSEC, PHOTHRD, PHOFOUR, PHOFIFT, PHOSIXT
      REAL*8  PHOSMV, PHOSEV
      INTEGER IPRHSO
      LOGICAL TESTZY, DOSO1, DOSO2, X2MAT, A2MAT, X2GRAD, PHOSPH,
     &        MNFPHO, ECPHOS, PHOSPV, CPPHOL, CPPHOV, CPPHMF, CPPHEC    !MK
      COMMON /INFHSO/PHOSMAT(3,3,MXPHOS), PHOSMN(3,3,MXPHOS),                 ! real*8
     &               PHOSEC(3,3,MXPHOS),  PHOTHRD(3,3,MXPHOS),          !MK
     &               PHOFOUR(3,3,MXPHOS), PHOFIFT(3,3,MXPHOS),          !MK
     &               PHOSIXT(3,3,MXPHOS), PHOSMV(3,3,MXPHOS),           !MK
     &               PHOSEV(3,3,MXPHOS),                                !MK
     &               IPRHSO,                                                  ! integer
     &               TESTZY, DOSO1,  DOSO2,  X2MAT,  A2MAT,  X2GRAD,          ! logical
     &               PHOSPH, MNFPHO, ECPHOS, PHOSPV, CPPHOL, CPPHOV,    !MK
     &               CPPHMF, CPPHEC                                     !MK
! --  end of include/infhso.h
