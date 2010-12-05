C
C  /* Deck so_onefock */
      SUBROUTINE SO_ONEFOCK(RES1E,LRES1E,RES1D,LRES1D,FOCKD,LFOCKD,
     &                      TR1E,LTR1E,TR1D,LTR1D,ISYRES,ISYMTR)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, December 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate one-particle part of the RPA Fock matrix
C              which goes into result vectors RES1E and RES1D.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1E(LRES1E), RES1D(LRES1D), TR1E(LTR1E), TR1D(LTR1D)
      DIMENSION FOCKD(LFOCKD)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_ONEFOCK')
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYRES)
C
         DO 200 I = 1,NRHF(ISYMI)
C
            DO 300 A = 1,NVIR(ISYMA)
C
               KOFF1 = IVIR(ISYMA) + A
               KOFF2 = IRHF(ISYMI) + I
               KOFF3 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA) * (I-1) + A
C
               FDIFF = FOCKD(KOFF1) - FOCKD(KOFF2)
C
               RES1E(KOFF3) = RES1E(KOFF3) + FDIFF * TR1E(KOFF3)
               RES1D(KOFF3) = RES1D(KOFF3) + FDIFF * TR1D(KOFF3)
C
  300       CONTINUE
C
  200    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_ONEFOCK')
C
      RETURN
      END
