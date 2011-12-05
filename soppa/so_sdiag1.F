C
C  /* Deck so_sdiag1 */
      SUBROUTINE SO_SDIAG1(DIAG1,LDIAG1,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                     ISYRES)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, April 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate diagonale one-particle part of the
C              SOPPA S[2] matrix.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION DIAG1(LDIAG1), DENSIJ(LDENSIJ), DENSAB(LDENSAB)
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
      CALL QENTER('SO_SDIAG1')
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMJ = ISYMI
         ISYMA = MULD2H(ISYMI,ISYRES)
         ISYMB = ISYMA
C
         KOFF1 = IIJDEN(ISYMI,ISYMJ)
         KOFF2 = IABDEN(ISYMA,ISYMB)
         KOFF3 = IT1AM(ISYMA,ISYMI)
C
         DO 200 I = 1,NRHF(ISYMI)
C
            KOFFIJ = KOFF1 + NRHF(ISYMI)*(I-1) + I
C
            DO 300 A = 1,NVIR(ISYMA)
C
               KOFFAB = KOFF2 + NVIR(ISYMA)*(A-1) + A
               KOFFAI = KOFF3 + NVIR(ISYMA)*(I-1) + A
C
               DIAG1(KOFFAI) = ONE + DENSIJ(KOFFIJ) - DENSAB(KOFFAB) 
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
      CALL QEXIT('SO_SDIAG1')
C
      RETURN
      END
