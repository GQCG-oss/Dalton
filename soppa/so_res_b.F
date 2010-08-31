C  /* Deck so_res_b */
      SUBROUTINE SO_RES_B(RES1E,LRES1E,RES1D,LRES1D,SIGAI1,LSIGAI1,
     &                    SIGAI2,LSIGAI2,CMO,LCMO,ISYRES)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak and Henrik Koch, September 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Calculate the third contributions to RES1E and RES1D
C              in eqs. (34) and (35).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1E(LRES1E),   RES1D(LRES1D), SIGAI1(LSIGAI1)
      DIMENSION SIGAI2(LSIGAI2), CMO(LCMO)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_B')
C
      DO 100 ISALFA = 1, NSYM
C
	 ISYMA  = ISALFA
	 ISYMI  = MULD2H(ISYMA,ISYRES)
C
         KOFF1  = ILMVIR(ISYMA) + 1
         KOFF2  = IT1AO(ISALFA,ISYMI) + 1
         KOFF3  = IT1AM(ISYMA,ISYMI) + 1

C
	 NTOTAL = MAX(NBAS(ISALFA),1)
	 NTOTA  = MAX(NVIR(ISYMA),1)
C
C----------------------------------------------------------------------
C        Multiply MO-coefficients C(alfa,a) and Sigma(alfa,i) matrices
C        to get third contributions to RES1E and RES1D in eqs. (34) and
C        (35).
C----------------------------------------------------------------------
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISALFA),
     &              -ONE,CMO(KOFF1),NTOTAL,SIGAI1(KOFF2),NTOTAL,
     &              ONE,RES1E(KOFF3),NTOTA)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISALFA),
     &              ONE,CMO(KOFF1),NTOTAL,SIGAI2(KOFF2),NTOTAL,
     &              ONE,RES1D(KOFF3),NTOTA)
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_B')
C
      RETURN
      END
