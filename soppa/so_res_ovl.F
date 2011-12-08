C
C  /* Deck so_res_ovl */
      SUBROUTINE SO_RES_OVL(RESO1E,LRESO1E,RESO1D,LRESO1D,TR1E,LTR1E,
     &                      TR1D,LTR1D,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                      ISYRES,ISYMTR)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, December 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate the overlap result vector, as expressed
C              in eq. (67).
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RESO1E(LRESO1E), RESO1D(LRESO1D)
      DIMENSION TR1E(LTR1E),     TR1D(LTR1D)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB)
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
      CALL QENTER('SO_RES_OVL')
C
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYRES)
         ISYMB = MULD2H(ISYMI,ISYMTR)
C
         KOFF1 = IABDEN(ISYMA,ISYMB) + 1
         KOFF2 = IT1AM(ISYMB,ISYMI)  + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTB = MAX(NVIR(ISYMB),1)
C
C----------------------------------------
C        Evaluate first term in eq. (68).
C----------------------------------------
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &              -ONE,DENSAB(KOFF1),NTOTA,TR1E(KOFF2),NTOTB,
     &              ONE,RESO1E(KOFF3),NTOTA)
C
C----------------------------------------
C        Evaluate first term in eq. (69).
C----------------------------------------
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &              ONE,DENSAB(KOFF1),NTOTA,TR1D(KOFF2),NTOTB,
     &              ONE,RESO1D(KOFF3),NTOTA)
C
         ISYMJ = MULD2H(ISYMA,ISYMTR)
C
         KOFF4 = IT1AM(ISYMA,ISYMJ)  + 1
         KOFF5 = IIJDEN(ISYMI,ISYMJ) + 1
C
         NTOTJ = MAX(NRHF(ISYMJ),1)
C
C-----------------------------------------
C        Evaluate second term in eq. (68).
C-----------------------------------------
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              ONE,TR1E(KOFF4),NTOTA,DENSIJ(KOFF5),NTOTJ,
     &              ONE,RESO1E(KOFF3),NTOTA)
C
C-----------------------------------------
C        Evaluate second term in eq. (69).
C-----------------------------------------
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              -ONE,TR1D(KOFF4),NTOTA,DENSIJ(KOFF5),NTOTJ,
     &              ONE,RESO1D(KOFF3),NTOTA)
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_OVL')
C
      RETURN
      END
