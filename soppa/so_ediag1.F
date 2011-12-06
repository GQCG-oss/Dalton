C
C  /* Deck so_ediag1 */
      SUBROUTINE SO_EDIAG1(DIAG1,LDIAG1,FOCKD,LFOCKD,ISYRES)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, April 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate diagonale one-particle part of the 
C              SOPPA E[2] matrix excluding the second order 
C              contributions.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION DIAG1(LDIAG1), FOCKD(LFOCKD)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_EDIAG1')
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
               DIAG1(KOFF3) = FOCKD(KOFF1) - FOCKD(KOFF2)
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
      CALL QEXIT('SO_EDIAG1')
C
      RETURN
      END
