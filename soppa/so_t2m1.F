C  /* Deck so_t2m1 */
      SUBROUTINE SO_T2M1(T2M1,LT2M1,T2MP,LT2MP,CMO,LCMO,IDEL,
     &                   ISYMD,ISYDIS,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak and Henrik Koch, September 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Calculate MP2 T2-amplitudes with one back-transformed
C              index.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION T2M1(LT2M1), T2MP(LT2MP), CMO(LCMO)
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
      CALL QENTER('SO_T2M1')
C
      ISYMB = ISYMD
C
      LCDB  = NVIR(ISYMB)
C
      KCDB    = 1
      KEND1   = KCDB  + LCDB
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_T2M1.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_T2M1.1',' ',KEND1,LWORK)
C
      KOFF1 = ILMVIR(ISYMD) + IDEL - IBAS(ISYMD)
C
C--------------------------------------------------
C     Copy delta MO-coefficients to the vector CDB.
C--------------------------------------------------
C
      CALL DCOPY(NVIR(ISYMB),CMO(KOFF1),NBAS(ISYMD),WORK(KCDB),1)
C
C------------------------------------------------------------------
C     Copy t2-amplitudes with B's having the symmetry of delta into
C     PT2MP.
C------------------------------------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMJ,ISYMB)
C
         DO 110 ISYMI = 1,NSYM
C
            ISYMAI  = ISYMBJ
            ISYMA   = MULD2H(ISYMI,ISYMAI)
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               LPT2MP  = NRHF(ISYMI)*NVIR(ISYMA)*NVIR(ISYMB)
C
               KPT2MP  = KEND1
               KEND2   = KPT2MP + LPT2MP
               LWORK2  = LWORK  - KEND2
C
               CALL SO_MEMMAX ('SO_T2M1.2',LWORK2)
               IF (LWORK2 .LT. 0) 
     &             CALL STOPIT('SO_T2M1.2',' ',KEND2,LWORK)
C
               NBJ1 = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1)
               NAIB = KPT2MP - 1
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  NBJ    = NBJ1 + B
C
                  NAIBJ1 = IT2AM(ISYMAI,ISYMBJ)
C
                  DO 140 I = 1,NRHF(ISYMI)
C
                     NAI1 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)
C
                     DO 141 A = 1,NVIR(ISYMA)
C
                        NAI   = NAI1   + A
                        NAIB  = NAIB   + 1
                        NAIBJ = NAIBJ1 + INDEX(NAI,NBJ)
C
                        WORK(NAIB) = T2MP(NAIBJ)
C
  141                CONTINUE
C
  140             CONTINUE
C
  130          CONTINUE
C
C------------------------------------------------------------
C              Calculate one-index transformed t2-amplitudes.
C------------------------------------------------------------
C
               KOFF2 = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J-1)
     &               + IT1AM(ISYMA,ISYMI) + 1
C
               NTOT  = MAX(NVIR(ISYMA)*NRHF(ISYMI),1)
C
               CALL DGEMV('N',NVIR(ISYMA)*NRHF(ISYMI),NVIR(ISYMB),ONE,
     &                    WORK(KPT2MP),NTOT,WORK(KCDB),1,ZERO,
     &                    T2M1(KOFF2),1)
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_T2M1')
C
      RETURN
      END
