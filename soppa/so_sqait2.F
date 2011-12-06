C
C  /* Deck so_sqait2 */
      SUBROUTINE SO_SQAIT2(T2SQ,LT2SQ,T2MP,LT2MP,ISYMA,ISYMI,ISYMB,
     &                     ISYMJ,J)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Set up T2-amplitudes for a given J and given
C              symmetries of A, I, and B.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION T2SQ(LT2SQ), T2MP(LT2MP)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C------------------------------
C     Statement function INDEX.
C------------------------------
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_SQAIT2')
C
      ISYMAI = MULD2H(ISYMA,ISYMI)
      ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
      IF (ISYMAI .NE. ISYMBJ) THEN
         WRITE(LUPRI,*) ' Mistake in SO_SQAIT2: ISYMAI .NE. ISYMBJ'
         CALL QUIT(' Mistake in SO_SQAIT ')
      END IF
C
      NBIA   = 0
C
      DO 100 A = 1,NVIR(ISYMA)
C
         DO 200 I = 1,NRHF(ISYMI)
C
            DO 300 B = 1,NVIR(ISYMB)
C
               NAI   = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
               NBJ   = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
               NBIA  = NBIA   + 1
               NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
C
               T2SQ(NBIA) = T2MP(NAIBJ)
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
      CALL QEXIT('SO_SQAIT2')
C
      RETURN
      END
