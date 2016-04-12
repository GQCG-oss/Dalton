C
C  /* Deck so_res_tcbt */
      SUBROUTINE SO_RES_TCBT(RES1E,LRES1E,RES1D,LRES1D,
     &                       TR2E,LTR2E,TR2D,LTR2D,DSRHF,LDSRHF,
     &                       CMO,LCMO,IDEL,ISDEL,ISYDIS,ISYMTR,
     &                       WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Patricio F. Provasi & Stephan P. A. Sauer: 25.11.2010
C     based on SO_RES_TCB
C
C     Pi 29.03.16: Copied from dalton20aosoppa
C                        ~
C     PURPOSE: Calculate C times b contribution to 1p1h resultvectors
C              for the triplet case
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1E(LRES1E), RES1D(LRES1D), DSRHF(LDSRHF)
      DIMENSION TR2E(LTR2E),   TR2D(LTR2D),   CMO(LCMO)
      DIMENSION WORK(LWORK)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "soppinf.h"
#include "infsop.h"
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
      CALL QENTER('SO_RES_TCBT')
cspas
      write(LUPRI,'(A,I3)')
     &     'in start of SO_RES_TCBT: printing RES1E for IDEL=',IDEL
      write(LUPRI,'(I8,1X,F14.8)') (I,RES1E(I),I=1,LRES1E)
C
C-------------------------
C     Calculate constants.
C-------------------------
C
      SQ2I  = ONE / DSQRT(TWO)
C
      ISYMB = ISDEL
C
      LCDB  = NVIR(ISYMB)
C
      KCDB    = 1
      KEND1   = KCDB  + LCDB
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_RES_TCBT.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_RES_TCBT.1',' ',KEND1,LWORK)
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
         CALL SO_MEMMAX ('SO_RES_TCBT.2',LWORK2)
         IF (LWORK2 .LT. 0) CALL STOPIT('SO_RES_TCBT.2',' ',KEND2,LWORK)
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
               CALL SO_MEMMAX ('SO_RES_TCBT.3',LWORK3)
               IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_TCBT.3',' ',KEND3,LWORK)
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
C              Loop over the three triplet trial vectors and
C              C matrices.
C---------------------------------------------------------------
C
               DO 400 ITTRIAL = 1, 3
C
C---------------------------------------------------------------
C              Generate b(a k b j) +/- b(b k a j) for given j in
C              KPTR2EA(a k b) and KPTR2EB(a k b)
C              (and KPTR2DA(a k b) and KPTR2DB(a k b))
C              for excitations (and de-excitations), respectively.
C---------------------------------------------------------------
C
                  LTTR2E = NVIR(ISYMA)*NRHF(ISYMK)
                  LTTR2D = NVIR(ISYMA)*NRHF(ISYMK)
                  LPTR2E = NVIR(ISYMA)*NRHF(ISYMK)*NVIR(ISYMB)
                  LPTR2D = NVIR(ISYMA)*NRHF(ISYMK)*NVIR(ISYMB)
C
                  KTTR2EA = KSCR3   + LSCR3
                  KTTR2DA = KTTR2EA + LTTR2E
                  KTTR2EB = KTTR2DA + LTTR2D
                  KTTR2DB = KTTR2EB + LTTR2E
                  KPTR2EA = KTTR2DB + LTTR2D
                  KPTR2DA = KPTR2EA + LPTR2E
                  KPTR2EB = KPTR2DA + LPTR2D
                  KPTR2DB = KPTR2EB + LPTR2E
                  KEND4   = KPTR2DB + LPTR2D
                  LWORK4  = LWORK   - KEND4
C
                  CALL SO_MEMMAX ('SO_RES_TCBT.4',LWORK4)
                  IF (LWORK4 .LT. 0)
     &                CALL STOPIT('SO_RES_TCBT.4',' ',KEND4,LWORK)
C
                  CALL SO_TBEXTRACT(ITTRIAL,WORK(KPTR2EA),LPTR2E,
     &                              WORK(KPTR2EB),LPTR2E,
     &                              WORK(KPTR2DA),LPTR2D,
     &                              WORK(KPTR2DB),LPTR2D,
     &                              TR2E,LTR2E,TR2D,LTR2D,
     &                              ISYMA,ISYMK,ISYMB,ISYMJ,J)
C
C-----------------------------------------------------------------------
C              Back transform the b index of the combinations of trial
C              vectors to delta in the AO basis:
C              Generate 1/2 * KTTR2EA/B(a k delta) from KPTR2EA/B(a k b)
C                   and 1/2 * KTTR2DA/B(a k delta) from KPTR2DA/B(a k b)
C              for given delta.
C-----------------------------------------------------------------------
C
                  NTOT   = MAX(NVIR(ISYMA)*NRHF(ISYMK),1)
C
                  CALL DGEMV('N',NVIR(ISYMA)*NRHF(ISYMK),NVIR(ISYMB),
     &                       HALF,WORK(KPTR2EA),NTOT,WORK(KCDB),1,ZERO,
     &                       WORK(KTTR2EA),1)
                  IF ( LPTR2E .EQ. 0 ) CALL DZERO(WORK(KTTR2EA),LTTR2E)
C
                  CALL DGEMV('N',NVIR(ISYMA)*NRHF(ISYMK),NVIR(ISYMB),
     &                       HALF,WORK(KPTR2DA),NTOT,WORK(KCDB),1,ZERO,
     &                       WORK(KTTR2DA),1)
                  IF ( LPTR2D .EQ. 0 ) CALL DZERO(WORK(KTTR2DA),LTTR2D)
C
C
                  IF (ITTRIAL .NE. 1) THEN
                     CALL DGEMV('N',NVIR(ISYMA)*NRHF(ISYMK),NVIR(ISYMB),
     &                          HALF,WORK(KPTR2EB),NTOT,WORK(KCDB),1,
     &                          ZERO,WORK(KTTR2EB),1)
                     IF ( LPTR2E .EQ. 0 )
     &                    CALL DZERO(WORK(KTTR2EB),LTTR2E)
C
                     CALL DGEMV('N',NVIR(ISYMA)*NRHF(ISYMK),NVIR(ISYMB),
     &                          HALF,WORK(KPTR2DB),NTOT,WORK(KCDB),1,
     &                          ZERO,WORK(KTTR2DB),1)
                     IF ( LPTR2D .EQ. 0 )
     &                    CALL DZERO(WORK(KTTR2DB),LTTR2D)
                  END IF
cspas
cspas write(LUPRI,'(A,I3)')
cspas&     'in SO_RES_TCBT: printing TTR2EA for IDEL=',IDEL
cspas write(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTTR2EA-1+I),I=1,LTTR2E)
cspas write(LUPRI,'(A,I3)')
cspas&     'in SO_RES_TCBT: printing TTR2DA for IDEL=',IDEL
cspas write(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTTR2DA-1+i),I=1,LTTR2D)
C
C
C----------------------------------------------------------------------
C              Calculate and add C * b contribution to RES1E and RES1D.
C              according to
C              $(c {\widetilde{a|jb}^k}) = \sum_a (c a | j \delta)
C              ~\underline{^{E/D}\widetilde{b}}~^{2p2h}_{bk\delta j}$
C              by multiplying two-electron integrals in SCR2
C              with the triplet 2p2h amplitudes in TTR2EA and TTR2DA.
C----------------------------------------------------------------------
C
cspas
cspas write(LUPRI,'(A,I3)')
cspas&     'in 1 of SO_RES_TCBT: printing RES1E for IDEL=',IDEL
cspas write(LUPRI,'(I8,1X,F14.8)') (I,RES1E(I),I=1,LRES1E)
C
                  NTOTC  = MAX(NVIR(ISYMC),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
                  KOFF3  = IT1AM(ISYMC,ISYMK) + 1
C
                  CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMK),
     &                       NVIR(ISYMA),-ONE,WORK(KSCR2),NTOTC,
     &                       WORK(KTTR2EA),NTOTA,ONE,RES1E(KOFF3),NTOTC)
C
C
cspas
cspas write(LUPRI,'(A,I3)')
cspas&     'in 2 of SO_RES_TCBT: printing RES1E for IDEL=',IDEL
cspas write(LUPRI,'(I8,1X,F14.8)') (I,RES1E(I),I=1,LRES1E)
C
                  CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMK),
     &                       NVIR(ISYMA),-ONE,WORK(KSCR2),NTOTC,
     &                       WORK(KTTR2DA),NTOTA,ONE,RES1D(KOFF3),NTOTC)
C
C----------------------------------------------------------------------
C              Calculate and add C * b contribution to RES1E and RES1D.
C              according to
C              $(i {\widetilde{k|jb}^a}) = \sum_k (i k  | j \delta)
C              \underline{^{E/D}\widetilde{b}}~^{2p2h}_{ak\delta j}$
C              by multiplying two-electron integrals in SCR3
C              with the triplet 2p2h amplitudes in TTR2EA and TTR2DA
C              for T1 and with the triplet 2p2h amplitudes in TTR2EB
C              and TTR2DB for T2 and T3.
C----------------------------------------------------------------------
C
cspas
cspas write(LUPRI,'(A,I3)')
cspas&     'in 3 of SO_RES_TCBT: printing RES1E for IDEL=',IDEL
cspas write(LUPRI,'(I8,1X,F14.8)') (I,RES1E(I),I=1,LRES1E)
C
                  NTOTK  = MAX(NRHF(ISYMK),1)
                  KOFF4  = IT1AM(ISYMA,ISYMI) + 1
C
                  IF (ITTRIAL .EQ. 1) THEN
                     CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     &                          NRHF(ISYMK),ONE,WORK(KTTR2EA),NTOTA,
     &                          WORK(KSCR3),NTOTK,ONE,RES1E(KOFF4),
     &                          NTOTA)
C
                     CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     &                          NRHF(ISYMK),ONE,WORK(KTTR2DA),NTOTA,
     &                          WORK(KSCR3),NTOTK,ONE,RES1D(KOFF4),
     &                          NTOTA)
                  ELSE
                     CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     &                          NRHF(ISYMK),ONE,WORK(KTTR2EB),NTOTA,
     &                          WORK(KSCR3),NTOTK,ONE,RES1E(KOFF4),
     &                          NTOTA)
C
                     CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     &                          NRHF(ISYMK),ONE,WORK(KTTR2DB),NTOTA,
     &                          WORK(KSCR3),NTOTK,ONE,RES1D(KOFF4),
     &                          NTOTA)
C
                  END IF
cspas
cspas write(LUPRI,'(A,I3)')
cspas&     'in 4 of SO_RES_TCBT: printing RES1E for IDEL=',IDEL
cspas write(LUPRI,'(I8,1X,F14.8)') (I,RES1E(I),I=1,LRES1E)
C
C
  400          END DO
C
  300       END DO
C
  200    END DO
C
  100 END DO
C
C-----------------------
C     Remove from trace.
C-----------------------
C
cspas
      write(LUPRI,'(A,I3)')
     &     'in end of SO_RES_TCBT: printing RES1E for IDEL=',IDEL
      write(LUPRI,'(I8,1X,F14.8)') (I,RES1E(I),I=1,LRES1E)
C
      CALL QEXIT('SO_RES_TCBT')
C
      RETURN
      END
