C
C$Id: infmp2.h,v 1.1.1.1 2001-02-08 13:33:29 hjj Exp $
C
C     Off-set to (bj) block of kappa(ai,bj) is IADR1(B,J)
      PARAMETER (MAXVIR = MAXORB)
      LOGICAL MP2FRO
      COMMON /INFMP2/ NFRMP2(8), NPHSYM(8), NPHMAX, NP, NH,
     &       IPHORD(MAXORB), INDXPH(MAXORB), IFRMP2(MAXORB), NPHTOT(8),
     &       MP2FRO
