C  /* Deck so_res_c */
      SUBROUTINE SO_RES_C(RES1E,LRES1E,RES1D,LRES1D,SIGDA1,LSIGDA1,
     &                    SIGDA2,LSIGDA2,CMO,LCMO,ISYRES)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Calculate the forth contributions to RES1E and RES1D
C              in eqs. (34) and (35).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1E(LRES1E),   RES1D(LRES1D), SIGDA1(LSIGDA1)
      DIMENSION SIGDA2(LSIGDA2), CMO(LCMO)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_C')
C
      DO 100 ISDEL = 1, NSYM
C
	 ISYMI  = ISDEL
	 ISYMA  = MULD2H(ISYMI,ISYRES)
C
         KOFF1  = IMATAV(ISDEL,ISYMA) + 1
         KOFF2  = ILMRHF(ISYMI) + 1
         KOFF3  = IT1AM(ISYMA,ISYMI) + 1

C
	 NTOTDL = MAX(NBAS(ISDEL),1)
	 NTOTA  = MAX(NVIR(ISYMA),1)
C
C-----------------------------------------------------------------------
C        Multiply MO-coefficients C(delta,i) and Sigma(delta,a) matrices
C        to get fourth contributions to RES1E and RES1D in eqs. (34) and
C        (35).
C-----------------------------------------------------------------------
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISDEL),
     &              ONE,SIGDA1(KOFF1),NTOTDL,CMO(KOFF2),NTOTDL,
     &              ONE,RES1E(KOFF3),NTOTA)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISDEL),
     &              -ONE,SIGDA2(KOFF1),NTOTDL,CMO(KOFF2),NTOTDL,
     &              ONE,RES1D(KOFF3),NTOTA)
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_C')
C
      RETURN
      END
