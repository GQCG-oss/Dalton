C  /* Deck so_res_tcb */
      SUBROUTINE SO_RES_TCB(RES1E,LRES1E,RES1D,LRES1D,
     &                      TR2EM1,LTR2E,TR2DM1,LTR2D,DSRHF,LDSRHF,
     &                      CMO,LCMO,IDEL,ISDEL,ISYDIS,ISYMTR,
     &                      WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, February 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C                        ~
C     PURPOSE: Calculate C times b contribution to 1p1h resultvectors
C              as described in eq. (72) and (73).
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1E(LRES1E), RES1D(LRES1D), DSRHF(LDSRHF)
      DIMENSION TR2EM1(LTR2E), TR2DM1(LTR2D), CMO(LCMO)
      DIMENSION WORK(LWORK)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_TCB')
C
      SQ2I  = 1.D0 / DSQRT(TWO)
C
      ISYMB = ISDEL
C
      LCDB  = NVIR(ISYMB)
      IF (LCDB.LE.0) GOTO 999
C
      KCDB    = 1
      KEND1   = KCDB  + LCDB
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_RES_TCB.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_RES_TCB.1',' ',KEND1,LWORK)
C
      KOFF1 = ILMVIR(ISDEL) + IDEL - IBAS(ISDEL)
C
C--------------------------------------------------
C     Copy delta MO-coefficients to the vector CDB.
C--------------------------------------------------
C
      CALL DCOPY(NVIR(ISYMB),CMO(KOFF1),NBAS(ISDEL),WORK(KCDB),1)
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMJ,ISYMB)
         ISALBE = MULD2H(ISYMJ,ISYDIS)
         ISYMCA = ISALBE
         ISYMKI = ISALBE
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1   = KEND1
         KEND2   = KSCR1 + LSCR1
         LWORK2  = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_TCB.2',LWORK2)
         IF (LWORK2 .LT. 0) CALL STOPIT('SO_RES_TCB.2',' ',KEND2,LWORK)
C
         DO 200 J = 1,NRHF(ISYMJ)
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE) * (J - 1) + 1
C
C-----------------------------------------------------------------------
C           Get a squared set of ( alfa beta | j delta ) for given j and
C           delta.
C-----------------------------------------------------------------------
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
            DO 300 ISYMK = 1,NSYM
C
               ISKBJ  = MULD2H(ISYMK,ISYMBJ)
               ISYMA  = MULD2H(ISKBJ,ISYMTR)
               ISYMC  = MULD2H(ISYMA,ISYMCA)
               ISYMI  = MULD2H(ISYMK,ISYMKI)
               ISYMAK = MULD2H(ISYMA,ISYMK)
               ISYMBK = MULD2H(ISYMB,ISYMK)
               ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
C-----------------------------------------------------------------
C              Generate ( a c | j delta ) and ( k i | j delta) for
C              given j and delta in KSCR2 and KSCR3, respectively.
C-----------------------------------------------------------------
C
               ISALFA = ISYMA
               ISBETA = ISYMC
C
               ISALF2 = ISYMI
               ISBET2 = ISYMK
C
               LSCR2  = NVIR(ISYMC)*NVIR(ISYMA)
               LSCR3  = NRHF(ISYMK)*NRHF(ISYMI)
C               LSCR4  = MAX(NBAS(ISALFA)*NVIR(ISYMA),
C     &                      NBAS(ISALF2)*NRHF(ISYMI))
               LSCR4  = MAX(NBAS(ISBETA)*NVIR(ISYMA),
     &                      NBAS(ISBET2)*NRHF(ISYMI))
C
               KSCR2  = KEND2
               KSCR3  = KSCR2 + LSCR2
               KSCR4  = KSCR3 + LSCR3
               KEND3  = KSCR4 + LSCR4
               LWORK3 = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_RES_TCB.3',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_TCB.3',' ',KEND3,LWORK)
C
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
C               NTOTC  = MAX(NVIR(ISYMC),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
C
               KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF3  = ILMVIR(ISYMA) + 1
               KOFF4  = ILMVIR(ISYMC) + 1
C
C              ( alpha, beta |..) * C(alpha, a) => ( beta, a | ..)
               CALL DGEMM('T','N',NBAS(ISBETA),NVIR(ISYMA),
     &                    NBAS(ISALFA),ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR4),NTOTBE)
C              ( beta, a | .. ) * C(beta, c) => ( a, c | ..)
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMC),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR4),NTOTBE,
     &                    CMO(KOFF4),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTA)


               NTOTAL = MAX(NBAS(ISALF2),1)
               NTOTBE = MAX(NBAS(ISBET2),1)
               NTOTK  = MAX(NRHF(ISYMK),1)
C
               KOFF5  = KSCR1 + IAODIS(ISALF2,ISBET2)
               KOFF6  = ILMRHF(ISYMI) + 1
               KOFF7  = ILMRHF(ISYMK) + 1
C
C              ( alpha, beta | .. ) * C(alpha,i) => ( beta, i | ..)
               CALL DGEMM('T','N',NBAS(ISBET2),NRHF(ISYMI),
     &                    NBAS(ISALF2),ONE,WORK(KOFF5),NTOTAL,
     &                    CMO(KOFF6),NTOTAL,ZERO,WORK(KSCR4),NTOTBE)
C
C              C( beta, k) * ( beta, i | ..)  = > ( k, i | j delta)
               CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),
     &                    NBAS(ISBET2),ONE,CMO(KOFF7),NTOTBE,
     &                    WORK(KSCR4),NTOTBE,ZERO,WORK(KSCR3),NTOTK)
C
               KOFFX2M1 = IT2BCD(ISYMAK,ISYMJ) + NT1AM(ISYMAK)*(J-1) +
     &                    IT1AM(ISYMA,ISYMK)+1
C
               NTOTC  = MAX(NVIR(ISYMC),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTK  = MAX(NRHF(ISYMK),1)
C
               KOFF3  = IT1AM(ISYMC,ISYMK) + 1
               KOFF4  = IT1AM(ISYMA,ISYMI) + 1
C
C               ( a, c | j, delta ) * x( a,k, j, delta ) = > x ( c, k)
               CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),NVIR(ISYMA),
     &                    SQ2I,WORK(KSCR2),NTOTA,TR2EM1(KOFFX2M1),NTOTA,
     &                    ONE,RES1E(KOFF3),NTOTC)
C              And for D vector
               CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),NVIR(ISYMA),
     &                    SQ2I,WORK(KSCR2),NTOTA,TR2DM1(KOFFX2M1),NTOTA,
     &                    ONE,RES1D(KOFF3),NTOTC)
C              ~
C              x(a,k,j,delta*(k,i | j,delta) => x( a, i)
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     &                    -SQ2I,TR2EM1(KOFFX2M1),NTOTA,WORK(KSCR3),
     &                    NTOTK,ONE,RES1E(KOFF4),NTOTA)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     &                    -SQ2I,TR2DM1(KOFFX2M1),NTOTA,WORK(KSCR3),
     &                    NTOTK,ONE,RES1D(KOFF4),NTOTA)
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
  999 CALL QEXIT('SO_RES_TCB')
C
      RETURN
      END
