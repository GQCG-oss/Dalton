C  /* Deck so_res_a */
      SUBROUTINE SO_RES_A(RES1E,LRES1E,RES1D,LRES1D,TR1E,LTR1E,TR1D,
     &                    LTR1D,DSRHF,LDSRHF,CMO,LCMO,T2M1,LT2M1,
     &                    AIJ,LAIJ,AAB,LAAB,INEWTR,ISDEL,ISYDIS,ISYRES,
     &                    ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, October 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Transform second and third index of integral batch and
C              multiply with trial vectors and t2-amplitudes to
C              calculate the first and second contributions to
C              RES1E and RES1D in eqs. (34) and (35).
C              In addition calculate Aij and Aab as defined through
C              eqs. (42) and (43).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1E(LRES1E), RES1D(LRES1D), TR1E(LTR1E), TR1D(LTR1D)
      DIMENSION DSRHF(LDSRHF), CMO(LCMO), T2M1(LT2M1), WORK(LWORK)
      DIMENSION AIJ(LAIJ),     AAB(LAAB)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_A')
C
C-----------------------------------------------------------
C     Test if T2M1 contains any elements. If not there is no
C     contribution.
C-----------------------------------------------------------
C
      IF (LT2M1 .EQ. 0) THEN
         CALL QEXIT('SO_RES_A')
         RETURN
      END IF
C
      ISYMC  = ISDEL
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMCK = MULD2H(ISYMK,ISYMC)
         ISALBE = MULD2H(ISYMK,ISYDIS)
C
         LSCR1  = N2BST(ISALBE)
         LSCR2  = NT1AM(ISALBE)
C
         KSCR1   = 1
         KSCR2   = KSCR1 + LSCR1
         KEND1   = KSCR2 + LSCR2
         LWORK1  = LWORK - KEND1
C
         CALL SO_MEMMAX ('SO_RES_A.1',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT('SO_RES_A.1',' ',KEND1,LWORK)
C
         DO 200 K = 1,NRHF(ISYMK)
C
            KOFF1 = IDSRHF(ISALBE,ISYMK) + NNBST(ISALBE)*(K-1) + 1
C
C-------------------------------------------------------------------
C           Get a squared set of ( alfa beta | k delta ) for given k
C           and delta.
C-------------------------------------------------------------------
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
            DO 300 ISYMA = 1,NSYM
C
               ISYMAL = ISYMA
               ISYMBE = MULD2H(ISYMAL,ISALBE)
               ISYMJ  = ISYMBE
               ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
               LSCR3  = NRHF(ISYMAL)*NBAS(ISYMJ)
C
               KSCR3   = KEND1
               KEND2   = KSCR3 + LSCR3
               LWORK2  = LWORK - KEND2
C
               CALL SO_MEMMAX ('SO_RES_A.2',LWORK2)
               IF (LWORK2 .LT. 0) 
     &             CALL STOPIT('SO_RES_A.2',' ',KEND2,LWORK)
C
               KOFF2  = KSCR1 + IAODIS(ISYMAL,ISYMBE)
               KOFF3  = ILMRHF(ISYMJ) + 1
               KOFF4  = ILMVIR(ISYMA) + 1
               KOFF5  = KSCR2 + IT1AM(ISYMA,ISYMJ)
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTJ  = MAX(NRHF(ISYMJ),1)
C
C-----------------------------------------------------------------------
C              Transformation of beta index to get ( alfa j | k delta ).
C-----------------------------------------------------------------------
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),
     &                    NBAS(ISYMBE),ONE,WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,ZERO,WORK(KSCR3),
     &                    NTOTAL)
C
C-------------------------------------------------------------------
C              Transformation of alfa index to get ( a j | k delta).
C-------------------------------------------------------------------
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMJ),
     &                    NBAS(ISYMAL),ONE,CMO(KOFF4),NTOTAL,
     &                    WORK(KSCR3),NTOTAL,ZERO,WORK(KOFF5),
     &                    NTOTA)
C
               IF (INEWTR .EQ. 1) THEN
C
C-----------------------------------------------------
C              Calculate Aij one time (for INEWTR = 1).
C-----------------------------------------------------
C
               ISYMI  = ISYMJ
               ISYMAI = MULD2H(ISYMI,ISYMA)
C
               KOFF6  = IT2BCD(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1)
     &                + IT1AM(ISYMA,ISYMI) + 1
               KOFF7  = IMATIJ(ISYMI,ISYMJ) + 1
C
               NTOTI  = MAX(NRHF(ISYMI),1)
C
               CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),
     &                    NVIR(ISYMA),ONE,T2M1(KOFF6),NTOTA,
     &                    WORK(KOFF5),NTOTA,ONE,AIJ(KOFF7),NTOTI)
C
C-----------------------------------------------------
C              Calculate Aab one time (for INEWTR = 1).
C-----------------------------------------------------
C
               ISYME  = ISYMA
               ISYMEJ = MULD2H(ISYME,ISYMJ)
C
               KOFF6  = IT2BCD(ISYMEJ,ISYMK) + NT1AM(ISYMEJ)*(K-1)
     &                + IT1AM(ISYME,ISYMJ) + 1
               KOFF7  = IMATAB(ISYME,ISYMA) + 1
C
               NTOTE  = MAX(NVIR(ISYME),1)
C
               CALL DGEMM('N','T',NVIR(ISYME),NVIR(ISYMA),
     &                    NRHF(ISYMJ),ONE,T2M1(KOFF6),NTOTE,
     &                    WORK(KOFF5),NTOTA,ONE,AAB(KOFF7),NTOTE)
C
C
               END IF
C
  300       CONTINUE
C
            DO 400 ISYMA = 1,NSYM
C
               ISYMI  = MULD2H(ISYMA,ISYRES)
C
C==================================
C              1. term in Eq. (34):
C==================================
C
               ISYMD  = MULD2H(ISYMI,ISYMCK)
C
               LSCR4  = NVIR(ISYMA)*NVIR(ISYMD)
C
               KSCR4   = KEND1
               KEND3   = KSCR4 + LSCR4
               LWORK3  = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_RES_A.3',LWORK3)
               IF (LWORK3 .LT. 0) 
     &             CALL STOPIT('SO_RES_A.3',' ',KEND3,LWORK)
C
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTD  = MAX(NVIR(ISYMD),1)
C
C---------------------------------------------------------------
C              MO-integrals times trial-vector to get 1. part of
C                ~ ~
C              ( a d | k delta).
C---------------------------------------------------------------
C
               ISYMJ  = MULD2H(ISYMA,ISYMTR)
C
               KOFF6  = IT1AM(ISYMA,ISYMJ) + 1
               KOFF7  = KSCR2 + IT1AM(ISYMD,ISYMJ)
C
               CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMD),NRHF(ISYMJ),
     &                    ONE,TR1E(KOFF6),NTOTA,WORK(KOFF7),NTOTD,
     &                    ZERO,WORK(KSCR4),NTOTA)
C
C---------------------------------------------------------------
C              MO-integrals times trial-vector to get 2. part of 
C                ~ ~
C              ( a d | k delta).
C---------------------------------------------------------------
C
               ISYMJ  = MULD2H(ISYMD,ISYMTR)
C
               KOFF8  = KSCR2 + IT1AM(ISYMA,ISYMJ)
               KOFF9  = IT1AM(ISYMD,ISYMJ) + 1
C
               CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMD),NRHF(ISYMJ),
     &                    -ONE,WORK(KOFF8),NTOTA,TR1D(KOFF9),NTOTD,
     &                    ONE,WORK(KSCR4),NTOTA)
C
C--------------------------------------------------
C                         ~ ~
C              Multiply ( a d | k delta) with T2M1.
C--------------------------------------------------
C
               ISYMDI = MULD2H(ISYMD,ISYMI)
C
               KOFF10 = IT2BCD(ISYMDI,ISYMK) + NT1AM(ISYMDI)*(K-1)
     &                + IT1AM(ISYMD,ISYMI) + 1
               KOFF11 = IT1AM(ISYMA,ISYMI) + 1
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NVIR(ISYMD),-ONE,WORK(KSCR4),NTOTA,
     &                    T2M1(KOFF10),NTOTD,ONE,RES1E(KOFF11),NTOTA)
C
C==================================
C              2. term in Eq. (34):
C==================================
C
               ISYML  = MULD2H(ISYMA,ISYMCK)
C
               LSCR4  = NRHF(ISYML)*NRHF(ISYMI)
C
               KSCR4   = KEND1
               KEND4   = KSCR4 + LSCR4
               LWORK4  = LWORK - KEND4
C
               CALL SO_MEMMAX ('SO_RES_A.4',LWORK4)
               IF (LWORK4 .LT. 0) 
     &             CALL STOPIT('SO_RES_A.4',' ',KEND4,LWORK)
C
               NTOTL  = MAX(NRHF(ISYML),1)
C
C---------------------------------------------------------------
C              MO-integrals times trial-vector to get 1. part of
C                ~ ~
C              ( l i | k delta).
C---------------------------------------------------------------
C
               ISYMB  = MULD2H(ISYML,ISYMTR)
C
               NTOTB  = MAX(NVIR(ISYMB),1)
C
               KOFF12 = IT1AM(ISYMB,ISYML) + 1
               KOFF13 = KSCR2 + IT1AM(ISYMB,ISYMI)
C
               CALL DGEMM('T','N',NRHF(ISYML),NRHF(ISYMI),NVIR(ISYMB),
     &                    ONE,TR1D(KOFF12),NTOTB,WORK(KOFF13),NTOTB,
     &                    ZERO,WORK(KSCR4),NTOTL)
C
C---------------------------------------------------------------
C              MO-integrals times trial-vector to get 2. part of
C                ~ ~
C              ( l i | k delta).
C---------------------------------------------------------------
C
               ISYMB  = MULD2H(ISYMI,ISYMTR)
C
               NTOTB  = MAX(NVIR(ISYMB),1)
C
               KOFF14 = KSCR2 + IT1AM(ISYMB,ISYML)
               KOFF15 = IT1AM(ISYMB,ISYMI) + 1
C
               CALL DGEMM('T','N',NRHF(ISYML),NRHF(ISYMI),NVIR(ISYMB),
     &                    -ONE,WORK(KOFF14),NTOTB,TR1E(KOFF15),NTOTB,
     &                    ONE,WORK(KSCR4),NTOTL)
C
C--------------------------------------------------
C                                   ~ ~
C              Multiply T2M1 with ( l i | k delta).
C--------------------------------------------------
C
               ISYAL  = MULD2H(ISYMA,ISYML)
C
               KOFF16 = IT2BCD(ISYAL,ISYMK) + NT1AM(ISYAL)*(K-1)
     &                + IT1AM(ISYMA,ISYML) + 1
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NRHF(ISYML),ONE,T2M1(KOFF16),NTOTA,
     &                    WORK(KSCR4),NTOTL,ONE,RES1E(KOFF11),NTOTA)
C
C==================================
C              1. term in Eq. (35):
C==================================
C
               ISYMD  = MULD2H(ISYMI,ISYMCK)
C
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTD  = MAX(NVIR(ISYMD),1)
C
C---------------------------------------------------------------
C              MO-integrals times trial-vector to get 1. part of
C                ~ ~
C              ( d a | k delta).
C---------------------------------------------------------------
C
               ISYMJ  = MULD2H(ISYMD,ISYMTR)
C
               KOFF17 = IT1AM(ISYMD,ISYMJ) + 1
               KOFF18 = KSCR2 + IT1AM(ISYMA,ISYMJ)
C
               CALL DGEMM('N','T',NVIR(ISYMD),NVIR(ISYMA),NRHF(ISYMJ),
     &                    ONE,TR1E(KOFF17),NTOTD,WORK(KOFF18),NTOTA,
     &                    ZERO,WORK(KSCR4),NTOTD)
C
C---------------------------------------------------------------
C              MO-integrals times trial-vector to get 2. part of
C                ~ ~
C              ( d a | k delta).
C---------------------------------------------------------------
C
               ISYMJ  = MULD2H(ISYMA,ISYMTR)
C
               KOFF19 = KSCR2 + IT1AM(ISYMD,ISYMJ)
               KOFF20 = IT1AM(ISYMA,ISYMJ) + 1
C
               CALL DGEMM('N','T',NVIR(ISYMD),NVIR(ISYMA),NRHF(ISYMJ),
     &                    -ONE,WORK(KOFF19),NTOTD,TR1D(KOFF20),NTOTA,
     &                    ONE,WORK(KSCR4),NTOTD)
C
C--------------------------------------------------
C                         ~ ~
C              Multiply ( d a | k delta) with T2M1.
C--------------------------------------------------
C
               ISYMDI = MULD2H(ISYMD,ISYMI)
C
               KOFF21 = IT2BCD(ISYMDI,ISYMK) + NT1AM(ISYMDI)*(K-1)
     &                + IT1AM(ISYMD,ISYMI) + 1
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NVIR(ISYMD),ONE,WORK(KSCR4),NTOTD,
     &                    T2M1(KOFF21),NTOTD,ONE,RES1D(KOFF11),NTOTA)
C
C==================================
C              2. term in Eq. (35):
C==================================
C
               ISYML  = MULD2H(ISYMA,ISYMCK)
C
               NTOTI  = MAX(NRHF(ISYMI),1)
C
C---------------------------------------------------------------
C              MO-integrals times trial-vector to get 1. part of
C                ~ ~
C              ( i l | k delta).
C---------------------------------------------------------------
C
               ISYMB  = MULD2H(ISYMI,ISYMTR)
C
               NTOTB  = MAX(NVIR(ISYMB),1)
C
               KOFF22 = IT1AM(ISYMB,ISYMI) + 1
               KOFF23 = KSCR2 + IT1AM(ISYMB,ISYML)
C
               CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYML),NVIR(ISYMB),
     &                    ONE,TR1D(KOFF22),NTOTB,WORK(KOFF23),NTOTB,
     &                    ZERO,WORK(KSCR4),NTOTI)
C
C---------------------------------------------------------------
C              MO-integrals times trial-vector to get 2. part of
C                ~ ~
C              ( i l | k delta).
C---------------------------------------------------------------
C
               ISYMB  = MULD2H(ISYML,ISYMTR)
C
               NTOTB  = MAX(NVIR(ISYMB),1)
C
               KOFF24 = KSCR2 + IT1AM(ISYMB,ISYMI)
               KOFF25 = IT1AM(ISYMB,ISYML) + 1
C
               CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYML),NVIR(ISYMB),
     &                    -ONE,WORK(KOFF24),NTOTB,TR1E(KOFF25),NTOTB,
     &                    ONE,WORK(KSCR4),NTOTI)
C
C--------------------------------------------------
C                                   ~ ~
C              Multiply T2M1 with ( i l | k delta).
C--------------------------------------------------
C
               ISYAL  = MULD2H(ISYMA,ISYML)
C
               KOFF26 = IT2BCD(ISYAL,ISYMK) + NT1AM(ISYAL)*(K-1)
     &                + IT1AM(ISYMA,ISYML) + 1
C
               CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),
     &                    NRHF(ISYML),-ONE,T2M1(KOFF26),NTOTA,
     &                    WORK(KSCR4),NTOTI,ONE,RES1D(KOFF11),NTOTA)
C
  400       CONTINUE
C
  200    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_A')
C
      RETURN
      END
