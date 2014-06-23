      LOGICAL TESTZY, DOSO1, DOSO2, X2MAT, A2MAT, X2GRAD, PHOSPH,
     &        MNFPHO, ECPHOS, PHOSPV, CPPHOL, CPPHOV, CPPHMF, CPPHEC    !MK
      PARAMETER (MXPHOS = 10)
      COMMON /INFHSO/PHOSMAT(3,3,MXPHOS), PHOSMN(3,3,MXPHOS),
     &               PHOSEC(3,3,MXPHOS), PHOTHRD(3,3,MXPHOS),           !MK
     &               PHOFOUR(3,3,MXPHOS), PHOFIFT(3,3,MXPHOS),          !MK
     &               PHOSIXT(3,3,MXPHOS), PHOSMV(3,3,MXPHOS),           !MK
     &               PHOSEV(3,3,MXPHOS),                                !MK
     &               IPRHSO, TESTZY, DOSO1, DOSO2, X2MAT, A2MAT, X2GRAD,
     &               PHOSPH, MNFPHO, ECPHOS, PHOSPV, CPPHOL, CPPHOV,    !MK
     &               CPPHMF, CPPHEC                                     !MK
