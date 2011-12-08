C
C  /* Deck so_res_fck */
      SUBROUTINE SO_RES_FCK(RES1E,LRES1E,RES1D,LRES1D,TR1E,LTR1E,TR1D,
     &                      LTR1D,FOCKD,LFOCKD,DENSIJ,LDENSIJ,DENSAB,
     &                      LDENSAB,ISYRES,ISYMTR)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, December 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate the Fock term to the A matrix, as expressed
C              in eq. (40).
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1E(LRES1E), RES1D(LRES1D)
      DIMENSION TR1E(LTR1E),   TR1D(LTR1D)
      DIMENSION FOCKD(LFOCKD), DENSIJ(LDENSIJ), DENSAB(LDENSAB)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <soppinf.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_FCK')
C
C-------------------------------------------------------
C     Add contribution from first fock-term in eq. (40).
C-------------------------------------------------------
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMB = MULD2H(ISYMI,ISYMTR)
         ISYMA = MULD2H(ISYMI,ISYRES)
C
         DO 200 I = 1,NRHF(ISYMI)
C
            KOFF1 = IRHF(ISYMI) + I
C
            DO 300 B = 1,NVIR(ISYMB)
C
               KOFF2 = IVIR(ISYMB) + B
               KOFF3 = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I-1) + B
C
               DO 400 A = 1,NVIR(ISYMA)
C
                  KOFF4 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
                  KOFF5 = IABDEN(ISYMA,ISYMB) + NVIR(ISYMA)*(B-1) + A
                  KOFF6 = IVIR(ISYMA) + A
C
                  EIBA  = DENSAB(KOFF5) * ( FOCKD(KOFF1)
     &                  - HALF * (FOCKD(KOFF2) + FOCKD(KOFF6)) )
C
                  RES1E(KOFF4) = RES1E(KOFF4) + EIBA * TR1E(KOFF3)
                  RES1D(KOFF4) = RES1D(KOFF4) + EIBA * TR1D(KOFF3)
C
  400          CONTINUE
C
  300       CONTINUE
C
  200    CONTINUE
C
  100 CONTINUE
C
C--------------------------------------------------------
C     Add contribution from second fock-term in eq. (40).
C--------------------------------------------------------
C
      DO 1100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYRES)
         ISYMJ = MULD2H(ISYMA,ISYMTR)
C
         DO 1200 I = 1,NRHF(ISYMI)
C
            KOFF1 = IRHF(ISYMI) + I
C
            DO 1300 J = 1,NRHF(ISYMJ)
C
               KOFF2 = IRHF(ISYMJ) + J
               KOFF3 = IIJDEN(ISYMI,ISYMJ) + NRHF(ISYMJ)*(I-1) + J
C
               DO 1400 A = 1,NVIR(ISYMA)
C
                  KOFF4 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
                  KOFF5 = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J-1) + A
                  KOFF6 = IVIR(ISYMA) + A
C
                  EAJI  = DENSIJ(KOFF3) * ( FOCKD(KOFF6)
     &                  - HALF * (FOCKD(KOFF2) + FOCKD(KOFF1)) )
C
                  RES1E(KOFF4) = RES1E(KOFF4) + EAJI * TR1E(KOFF5)
                  RES1D(KOFF4) = RES1D(KOFF4) + EAJI * TR1D(KOFF5)
C
 1400          CONTINUE
C
 1300       CONTINUE
C
 1200    CONTINUE
C
 1100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_FCK')
C
      RETURN
      END
