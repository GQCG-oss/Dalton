C  /* Deck so_res_c */
      SUBROUTINE SO_RES_C(RES1,LRES1,SIGDA,LSIGDA,
     &                    CMO,LCMO,ISYRES,FACTOR)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Pi Haase 23.03.2016: Do only E or D (copied from so_res_b)
C
C     PURPOSE: Calculate the forth contributions to RES1E and RES1D
C              in eqs. (34) and (35).
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      double precision, intent(in) :: FACTOR
      DIMENSION RES1(LRES1), SIGDA(LSIGDA)
      DIMENSION CMO(LCMO)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
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
CPi-210316
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISDEL),
     &              FACTOR,SIGDA(KOFF1),NTOTDL,CMO(KOFF2),NTOTDL,
     &              ONE,RES1(KOFF3),NTOTA)
C
C         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISDEL),
C     &              -ONE,SIGDA2(KOFF1),NTOTDL,CMO(KOFF2),NTOTDL,
C     &              ONE,RES1D(KOFF3),NTOTA)
C
Cend-Pi
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
