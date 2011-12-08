C  /* Deck so_res_tcb */
      SUBROUTINE SO_RES_TCB(RES1E,LRES1E,RES1D,LRES1D,
     &                      TR2E,LTR2E,TR2D,LTR2D,DSRHF,LDSRHF,
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
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1E(LRES1E), RES1D(LRES1D), DSRHF(LDSRHF)
      DIMENSION TR2E(LTR2E),   TR2D(LTR2D),   CMO(LCMO)
      DIMENSION WORK(LWORK)
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
      CALL QENTER('SO_RES_TCB')
C
      SQ2I  = 1.D0 / DSQRT(TWO)
C
      ISYMB = ISDEL
C
      LCDB  = NVIR(ISYMB)
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
C              Generate ( c a | j delta ) and ( k i | j delta) for
C              given j and delta in KSCR2 and KSCR3, respectively.
C-----------------------------------------------------------------
C
               ISALFA = ISYMC
               ISBETA = ISYMA
C
               ISALF2 = ISYMK
               ISBET2 = ISYMI
C
               LSCR2  = NVIR(ISYMC)*NVIR(ISYMA)
               LSCR3  = NRHF(ISYMK)*NRHF(ISYMI)
               LSCR4  = MAX(NBAS(ISALFA)*NVIR(ISYMA),
     &                      NBAS(ISALF2)*NRHF(ISYMI))
C
               KSCR2   = KEND2
               KSCR3   = KSCR2 + LSCR2
               KSCR4   = KSCR3 + LSCR3
               KEND3   = KSCR4 + LSCR4
               LWORK3  = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_RES_TCB.3',LWORK3)
              IF (LWORK3 .LT. 0) 
     &             CALL STOPIT('SO_RES_TCB.3',' ',KEND3,LWORK)
C
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
               NTOTC  = MAX(NVIR(ISYMC),1)
C
               KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF3  = ILMVIR(ISYMA) + 1
               KOFF4  = ILMVIR(ISYMC) + 1
C
               CALL DGEMM('N','N',NBAS(ISALFA),NVIR(ISYMA),
     &                    NBAS(ISBETA),ONE,WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,ZERO,WORK(KSCR4),NTOTAL)
C
               CALL DGEMM('T','N',NVIR(ISYMC),NVIR(ISYMA),
     &                    NBAS(ISALFA),ONE,CMO(KOFF4),NTOTAL,
     &                    WORK(KSCR4),NTOTAL,ZERO,WORK(KSCR2),NTOTC)
C
               NTOTAL = MAX(NBAS(ISALF2),1)
               NTOTBE = MAX(NBAS(ISBET2),1)
               NTOTK  = MAX(NRHF(ISYMK),1)
C
               KOFF5  = KSCR1 + IAODIS(ISALF2,ISBET2)
               KOFF6  = ILMRHF(ISYMI) + 1
               KOFF7  = ILMRHF(ISYMK) + 1
C
               CALL DGEMM('N','N',NBAS(ISALF2),NRHF(ISYMI),
     &                    NBAS(ISBET2),ONE,WORK(KOFF5),NTOTAL,
     &                    CMO(KOFF6),NTOTBE,ZERO,WORK(KSCR4),NTOTAL)
C
               CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),
     &                    NBAS(ISALF2),ONE,CMO(KOFF7),NTOTAL,
     &                    WORK(KSCR4),NTOTAL,ZERO,WORK(KSCR3),NTOTK)
C
C---------------------------------------------------------------
C              Generate 2*b(a k b j) - b(b k a j) for given j in
C              KPTR2E(a k b) (and KPTR2D(a k b)) for excitations
C              (and de-excitations), respectively.
C---------------------------------------------------------------
C
               LTTR2E = NVIR(ISYMA)*NRHF(ISYMK)
               LTTR2D = NVIR(ISYMA)*NRHF(ISYMK)
               LPTR2E = NVIR(ISYMA)*NRHF(ISYMK)*NVIR(ISYMB)
               LPTR2D = NVIR(ISYMA)*NRHF(ISYMK)*NVIR(ISYMB)
C
               KTTR2E  = KSCR3  + LSCR3
               KTTR2D  = KTTR2E + LTTR2E
               KPTR2E  = KTTR2D + LTTR2D
               KPTR2D  = KPTR2E + LPTR2E
               KEND4   = KPTR2D + LPTR2D
               LWORK4  = LWORK  - KEND4
C
               CALL SO_MEMMAX ('SO_RES_TCB.4',LWORK4)
              IF (LWORK4 .LT. 0) 
     &             CALL STOPIT('SO_RES_TCB.4',' ',KEND4,LWORK)
C
               CALL SO_BEXTRACT(WORK(KPTR2E),LPTR2E,WORK(KPTR2D),LPTR2D,
     &                          TR2E,LTR2E,TR2D,LTR2D,
     &                          ISYMA,ISYMK,ISYMB,ISYMJ,J)
C
C-----------------------------------------------------------------------
C              Generate 1/SQ2 * KTTR2E(a k delta) from KPTR2E(a k b) and
C                       1/SQ2 * KTTR2D(a k delta) from KPTR2D(a k b)
C              for given delta..
C-----------------------------------------------------------------------
C
               NTOT   = MAX(NVIR(ISYMA)*NRHF(ISYMK),1)
C
               CALL DGEMV('N',NVIR(ISYMA)*NRHF(ISYMK),NVIR(ISYMB),SQ2I,
     &                    WORK(KPTR2E),NTOT,WORK(KCDB),1,ZERO,
     &                    WORK(KTTR2E),1)
               IF ( LPTR2E .EQ. 0 ) CALL DZERO(WORK(KTTR2E),LTTR2E)
C
               CALL DGEMV('N',NVIR(ISYMA)*NRHF(ISYMK),NVIR(ISYMB),SQ2I,
     &                    WORK(KPTR2D),NTOT,WORK(KCDB),1,ZERO,
     &                    WORK(KTTR2D),1)
               IF ( LPTR2D .EQ. 0 ) CALL DZERO(WORK(KTTR2D),LTTR2D)
C
C----------------------------------------------------------------------
C              Calculate and add C * b contribution to RES1E and RES1D.
C              according to eqs. (72) and (73) by multiplying
C              two-electron integrals in SCR2 and SCR3 with with
C              2p2h amplitudes in TTR2E, and TTR2D.
C----------------------------------------------------------------------
C
               NTOTC  = MAX(NVIR(ISYMC),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTK  = MAX(NRHF(ISYMK),1)
C
               KOFF3  = IT1AM(ISYMC,ISYMK) + 1
               KOFF4  = IT1AM(ISYMA,ISYMI) + 1
C
               CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMK),NVIR(ISYMA),
     &                    ONE,WORK(KSCR2),NTOTC,WORK(KTTR2E),NTOTA,
     &                    ONE,RES1E(KOFF3),NTOTC)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     &                    -ONE,WORK(KTTR2E),NTOTA,WORK(KSCR3),NTOTK,
     &                    ONE,RES1E(KOFF4),NTOTA)
C
C
               CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMK),NVIR(ISYMA),
     &                    ONE,WORK(KSCR2),NTOTC,WORK(KTTR2D),NTOTA,
     &                    ONE,RES1D(KOFF3),NTOTC)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     &                    -ONE,WORK(KTTR2D),NTOTA,WORK(KSCR3),NTOTK,
     &                    ONE,RES1D(KOFF4),NTOTA)
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
      CALL QEXIT('SO_RES_TCB')
C
      RETURN
      END
