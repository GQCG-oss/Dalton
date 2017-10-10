C  /* Deck so_res_td1b1 */
      SUBROUTINE SO_RES_TD1B1(RES2E,LRES2E,RES2D,LRES2D,
     &                      TR2E,LTR2E,TR2D,LTR2D,DSRHF,LDSRHF,
     &                      CMO,LCMO,IDEL,ISDEL,ISYDIS,ISYMTR,
     &                      WORK,LWORK)
C
C     This routine is part of the atomic integral direct
C     SOPPA-D1 program.
C
C     Mats Simmermacher & Stephan P. A. Sauer, October 2014
C
C                        ~
C     PURPOSE: Calculate D(1)_1 times b contribution to 2p2h
C     resultvectors.
C
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES2E(LRES2E), RES2D(LRES2D), DSRHF(LDSRHF)
      DIMENSION TR2E(LTR2E),   TR2D(LTR2D),   CMO(LCMO)
      DIMENSION WORK(LWORK)
      INTEGER NAIT, NBJT
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
C-----------------------------------------------------------------------
C     Statement function INDEX.
C-----------------------------------------------------------------------
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C-----------------------------------------------------------------------
C     Add to trace.
C-----------------------------------------------------------------------
C
      CALL QENTER('SO_RES_TD1B1')
C
      ISYMC = ISDEL
C
      LCDB  = NVIR(ISYMC)
C
      KCDB    = 1
      KEND1   = KCDB  + LCDB
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_RES_TD1B1.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_RES_TD1B1.1',' ',KEND1,LWORK)
C
      KOFF1 = ILMVIR(ISDEL) + IDEL - IBAS(ISDEL)
C
C-----------------------------------------------------------------------
C     Copy delta MO-coefficients to the vector CDB.
C-----------------------------------------------------------------------
C
      CALL DCOPY(NVIR(ISYMC),CMO(KOFF1),NBAS(ISDEL),WORK(KCDB),1)
C
C     WRITE(LUPRI,'(A)') 'Check of Contributions: New Call'
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMKC = MULD2H(ISYMK,ISYMC)
         ISALBE = MULD2H(ISYMK,ISYDIS)
         ISYMAI = ISALBE
         ISYMBJ = MULD2H(ISYMKC,ISYMTR)
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1   = KEND1
         KEND2   = KSCR1 + LSCR1
         LWORK2  = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_TD1B1.2',LWORK2)
         IF (LWORK2 .LT. 0) CALL STOPIT('SO_RES_TD1B1.2',' ',
     &                                  KEND2,LWORK)
C
         DO 200 K = 1,NRHF(ISYMK)
C
            KOFF1 = IDSRHF(ISALBE,ISYMK) + NNBST(ISALBE) * (K - 1) + 1
C
C-----------------------------------------------------------------------
C           Get a squared set of (alfa beta | k delta) for given k and
C           delta.
C-----------------------------------------------------------------------
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
            DO 300 ISYMA = 1,NSYM
C
               ISYMI  = MULD2H(ISYMA,ISYMAI)
C
C-----------------------------------------------------------------------
C              Generate (a i | k delta) for given k and delta in KSCR2.
C-----------------------------------------------------------------------
C
               LSCR2  = NVIR(ISYMA)*NRHF(ISYMI)
C
               IF (LSCR2 .NE. 0) THEN
C
               ISALFA = ISYMA
               ISBETA = ISYMI
C
               LSCR3  = NBAS(ISALFA)*NRHF(ISYMI)
C
               KSCR2   = KEND2
               KSCR3   = KSCR2 + LSCR2
               KEND3   = KSCR3 + LSCR3
               LWORK3  = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_RES_TD1B1.3',LWORK3)
               IF (LWORK3 .LT. 0) 
     &         CALL STOPIT('SO_RES_TD1B1.3',' ',KEND3,LWORK)
C
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
C
               KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF3  = ILMRHF(ISYMI) + 1
               KOFF4  = ILMVIR(ISYMA) + 1
C
               CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,ZERO,WORK(KSCR3),NTOTAL)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,CMO(KOFF4),NTOTAL,
     &                    WORK(KSCR3),NTOTAL,ZERO,WORK(KSCR2),NTOTA)
C
C-----------------------------------------------------------------------
C              Generate 2*b(b j c k) - b(c j b k) for given k in
C              KPTR2E(b j c) (and KPTR2D(b j c)) for excitations
C              (and de-excitations), respectively.
C-----------------------------------------------------------------------
C
            DO 400 ISYMB = 1,NSYM
C
               ISYMJ  = MULD2H(ISYMB,ISYMBJ)
C
               LTTR2E = NVIR(ISYMB)*NRHF(ISYMJ)
C
               IF (LTTR2E .NE. 0) THEN
C
               LTTR2D = NVIR(ISYMB)*NRHF(ISYMJ)
               LPTR2E = NVIR(ISYMB)*NRHF(ISYMJ)*NVIR(ISYMC)
               LPTR2D = NVIR(ISYMB)*NRHF(ISYMJ)*NVIR(ISYMC)
C
               KTTR2E  = KSCR2  + LSCR2
               KTTR2D  = KTTR2E + LTTR2E
               KPTR2E  = KTTR2D + LTTR2D
               KPTR2D  = KPTR2E + LPTR2E
               KEND4   = KPTR2D + LPTR2D
               LWORK4  = LWORK  - KEND4
C
               CALL SO_MEMMAX ('SO_RES_TD1B1.4',LWORK4)
               IF (LWORK4 .LT. 0) 
     &         CALL STOPIT('SO_RES_TD1B1.4',' ',KEND4,LWORK)
C
               CALL SO_BEXTRACT(WORK(KPTR2E),LPTR2E,WORK(KPTR2D),
     &                          LPTR2D,TR2E,LTR2E,TR2D,LTR2D,
     &                          ISYMB,ISYMJ,ISYMC,ISYMK,K)
C
C-----------------------------------------------------------------------
C              Generate KTTR2E(b j delta) from KPTR2E(b j c) and
C                       KTTR2D(b j delta) from KPTR2D(b j c)
C              for given delta.
C-----------------------------------------------------------------------
C
               NTOT   = MAX(NVIR(ISYMB)*NRHF(ISYMJ),1)
C
               CALL DGEMV('N',NVIR(ISYMB)*NRHF(ISYMJ),NVIR(ISYMC),ONE,
     &                    WORK(KPTR2E),NTOT,WORK(KCDB),1,ZERO,
     &                    WORK(KTTR2E),1)
               IF ( LPTR2E .EQ. 0 ) CALL DZERO(WORK(KTTR2E),LTTR2E)
C
               CALL DGEMV('N',NVIR(ISYMB)*NRHF(ISYMJ),NVIR(ISYMC),ONE,
     &                    WORK(KPTR2D),NTOT,WORK(KCDB),1,ZERO,
     &                    WORK(KTTR2D),1)
               IF ( LPTR2D .EQ. 0 ) CALL DZERO(WORK(KTTR2D),LTTR2D)
C
C-----------------------------------------------------------------------
C              Calculate and add
C              (a i | k delta) * (2*b(b j delta k) - b(delta j b k)) and
C              (b j | k delta) * (2*b(a i delta k) - b(delta i a k))
C              to 2p2h excitation and de-excitation result-vectors.
C-----------------------------------------------------------------------
C
C              WRITE(LUPRI,'(A)') 'Check of Contributions: NBJ'
C              WRITE(LUPRI,'(A)')
C    &               'IDEL,ISYMTR,ISYMK,K,ISYMA,ISYMB,ISYMAI,ISYMBJ'
C              WRITE(LUPRI,'(8(I3,2X))')
C    &               IDEL,ISYMTR,ISYMK,K,ISYMA,ISYMB,ISYMAI,ISYMBJ
C
               DO 501 NBJ = 1,NVIR(ISYMB)*NRHF(ISYMJ)
C
                  NBJT = IT1AM(ISYMB,ISYMJ) + NBJ
C
C--------------ai = bj-------------------------------------------------- 
C
               IF ( ISYMAI .EQ. ISYMBJ ) THEN
C
C              WRITE(LUPRI,'(A)') 'Part 1: ai = bj'
C
               DO 601 NAI = 1,NVIR(ISYMA)*NRHF(ISYMI)
C
                  NAIT = IT1AM(ISYMA,ISYMI) + NAI
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAIT,NBJT)
C
C                 WRITE(LUPRI,'(A)')
C    &                  'NAIBJ,NAI,NBJ,KSCR2,KTTR2E,KTTR2D'
C                 WRITE(LUPRI,'(6(I3,2X))')
C    &                  NAIBJ,NAI,NBJ,KSCR2,KTTR2E,KTTR2D
C
                  RES2E(NAIBJ) = RES2E(NAIBJ) 
     &                           + (WORK(KSCR2 + NAI - 1)
     &                           * WORK(KTTR2E + NBJ - 1))
C                 WRITE(LUPRI,'(3(E14.7,2X))') RES2E(NAIBJ),
C    &                  WORK(KSCR2 + NAI - 1),WORK(KTTR2E + NBJ - 1)
                  RES2D(NAIBJ) = RES2D(NAIBJ)
     &                           + (WORK(KSCR2 + NAI - 1)
     &                           * WORK(KTTR2D + NBJ - 1))
C                 WRITE(LUPRI,'(3(E14.7,2X))') RES2D(NAIBJ),
C    &                  WORK(KSCR2 + NAI - 1),WORK(KTTR2D + NBJ - 1)
C
  601          CONTINUE
C
C--------------ai < bj--------------------------------------------------
C
               ELSE IF ( ISYMAI .LT. ISYMBJ ) THEN
C
C              WRITE(LUPRI,'(A)') 'Part 2: ai < bj'
C               
               DO 602 NAI = 1,NVIR(ISYMA)*NRHF(ISYMI)
C
                  NAIT = IT1AM(ISYMA,ISYMI) + NAI
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)
     &                    *(NBJT - 1) + NAIT
C
C                 WRITE(LUPRI,'(A)')
C    &                  'NAIBJ,NAI,NBJ,KSCR2,KTTR2E,KTTR2D'
C                 WRITE(LUPRI,'(6(I3,2X))')
C    &                  NAIBJ,NAI,NBJ,KSCR2,KTTR2E,KTTR2D
C
                  RES2E(NAIBJ) = RES2E(NAIBJ)
     &                           + WORK(KSCR2 + NAI - 1)
     &                           * WORK(KTTR2E + NBJ - 1)
C                 WRITE(LUPRI,'(3(E14.7,2X))') RES2E(NAIBJ),
C    &                  WORK(KSCR2 + NAI - 1),WORK(KTTR2E + NBJ - 1)
                  RES2D(NAIBJ) = RES2D(NAIBJ)
     &                           + WORK(KSCR2 + NAI - 1)
     &                           * WORK(KTTR2D + NBJ - 1)
C                 WRITE(LUPRI,'(3(E14.7,2X))') RES2D(NAIBJ),
C    &                  WORK(KSCR2 + NAI - 1),WORK(KTTR2D + NBJ - 1)
C
  602          CONTINUE
C
C--------------ai > bj--------------------------------------------------
C
               ELSE IF ( ISYMAI .GT. ISYMBJ ) THEN
C
C              WRITE(LUPRI,'(A)') 'Part 3: ai > bj'
C
               DO 603 NAI = 1,NVIR(ISYMA)*NRHF(ISYMI)
C
                  NAIT = IT1AM(ISYMA,ISYMI) + NAI
                  NAIBJ = IT2AM(ISYMBJ,ISYMAI) + NT1AM(ISYMBJ)
     &                    *(NAIT - 1) + NBJT
C
C                 WRITE(LUPRI,'(A)')
C    &                  'NAIBJ,NAI,NBJ,KSCR2,KTTR2E,KTTR2D'
C                 WRITE(LUPRI,'(6(I3,2X))')
C    &                  NAIBJ,NAI,NBJ,KSCR2,KTTR2E,KTTR2D
C
                  RES2E(NAIBJ) = RES2E(NAIBJ)
     &                           + WORK(KSCR2 + NAI - 1)
     &                           * WORK(KTTR2E + NBJ - 1)
C                 WRITE(LUPRI,'(3(E14.7,2X))') RES2E(NAIBJ),
C    &                  WORK(KSCR2 + NAI - 1),WORK(KTTR2E + NBJ - 1)
                  RES2D(NAIBJ) = RES2D(NAIBJ)
     &                           + WORK(KSCR2 + NAI - 1)
     &                           * WORK(KTTR2D + NBJ - 1)
C                 WRITE(LUPRI,'(3(E14.7,2X))') RES2D(NAIBJ),
C    &                  WORK(KSCR2 + NAI - 1),WORK(KTTR2D + NBJ - 1)
C
  603          CONTINUE
C
               END IF
C
  501          CONTINUE
C
C-----------------------------------------------------------------------
C
               END IF
C
  400       CONTINUE
C
               END IF
C
  300       CONTINUE
C
  200    CONTINUE
C
  100 CONTINUE
C
C-----------------------------------------------------------------------
C     Remove from trace.
C-----------------------------------------------------------------------
C
      CALL QEXIT('SO_RES_TD1B1')
C
      RETURN
      END
