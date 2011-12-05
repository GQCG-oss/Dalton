C  /* Deck so_res_sym */
      SUBROUTINE SO_RES_SYM(RES1E,LRES1E,RES1D,LRES1D,AIJ,LAIJ,AAB,
     &                      LAAB,TR1E,LTR1E,TR1D,LTR1D,ISYRES)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate the symmetry correcting term to the A matrix,
C              as expressed in eq. (44).
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1E(LRES1E), RES1D(LRES1D), AIJ(LAIJ),  AAB(LAAB)
      DIMENSION TR1E(LTR1E),   TR1D(LTR1D)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <soppinf.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_SYM')
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI = MULD2H(ISYMA,ISYRES)
C
         ISYMJ = ISYMI
C
         KOFF1 = IT1AM(ISYMA,ISYMJ)  + 1
         KOFF2 = IMATIJ(ISYMI,ISYMJ) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1
         KOFF4 = IMATIJ(ISYMJ,ISYMI) + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTI = MAX(NRHF(ISYMI),1)
         NTOTJ = MAX(NRHF(ISYMJ),1)
C
C-----------------------------------------------------------------------
C        Multiply trial-vectors with Aij elements according to eq. (44).
C-----------------------------------------------------------------------
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              HALF,TR1E(KOFF1),NTOTA,AIJ(KOFF2),NTOTI,ONE,
     &              RES1E(KOFF3),NTOTA)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              -HALF,TR1E(KOFF1),NTOTA,AIJ(KOFF4),NTOTJ,ONE,
     &              RES1E(KOFF3),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              HALF,TR1D(KOFF1),NTOTA,AIJ(KOFF2),NTOTI,ONE,
     &              RES1D(KOFF3),NTOTA)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              -HALF,TR1D(KOFF1),NTOTA,AIJ(KOFF4),NTOTJ,ONE,
     &              RES1D(KOFF3),NTOTA)
C
         ISYMB = ISYMA
C
         KOFF1 = IMATAB(ISYMA,ISYMB) + 1
         KOFF2 = IT1AM(ISYMB,ISYMI)  + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1
         KOFF4 = IMATAB(ISYMB,ISYMA) + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTB = MAX(NVIR(ISYMB),1)
C
C-----------------------------------------------------------------------
C        Multiply trial-vectors with Aab elements according to eq. (44).
C-----------------------------------------------------------------------
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &              HALF,AAB(KOFF1),NTOTA,TR1E(KOFF2),NTOTB,ONE,
     &              RES1E(KOFF3),NTOTA)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &              -HALF,AAB(KOFF4),NTOTB,TR1E(KOFF2),NTOTB,ONE,
     &              RES1E(KOFF3),NTOTA)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &              HALF,AAB(KOFF1),NTOTA,TR1D(KOFF2),NTOTB,ONE,
     &              RES1D(KOFF3),NTOTA)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &              -HALF,AAB(KOFF4),NTOTB,TR1D(KOFF2),NTOTB,ONE,
     &              RES1D(KOFF3),NTOTA)
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_SYM')
C
      RETURN
      END
