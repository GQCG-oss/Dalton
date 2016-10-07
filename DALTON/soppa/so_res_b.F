C  /* Deck so_res_b */
      SUBROUTINE SO_RES_B(RES1,LRES1,SIGAI,LSIGAI,
     &                    CMO,LCMO,ISYRES,FACTOR)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak and Henrik Koch, September 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Rasmus Faber: 22.1.2016: Do only E or D
C              Also Sigma has been reordered as (i,alfa), though
C              the symmetry blocks are still ordered (ialpha,isymi)
C
C     PURPOSE: Calculate the third contributions to RES1E or RES1D
C              in eqs. (34) and (35).
C
C#include "implicit.h"
      implicit none
#include "priunit.h"

      double precision, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0,
     &                               ONE = 1.0D0, TWO = 2.0D0
C
      double precision, intent(inout) ::  RES1(LRES1)
      double precision, intent(in) ::  CMO(LCMO), SIGAI(LSIGAI)
      double precision, intent(in) :: FACTOR
      integer, intent(in)          :: LRES1, LSIGAI, LCMO, ISYRES
C
      integer :: ISALFA, ISYMA, ISYMI, KOFF1, KOFF2, KOFF3
      integer :: NTOTAL, NTOTA, NTOTI
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
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
         NTOTI = MAX(1,NRHF(ISYMI))
C
C----------------------------------------------------------------------
C        Multiply MO-coefficients C(alfa,a) and Sigma(alfa,i) matrices
C        to get third contributions to RES1E and RES1D in eqs. (34) and
C        (35).
CRF      Sigma now stored as (i,alfa) -> This should allow for more
C        efficient code in so_sigai.
C----------------------------------------------------------------------
C
         CALL DGEMM('T','T',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISALFA),
     &              FACTOR,CMO(KOFF1),NTOTAL,SIGAI(KOFF2),NTOTI,
     &              ONE,RES1(KOFF3),NTOTA)
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
