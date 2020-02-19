C  /* Deck so_aibj2 */
      SUBROUTINE SO_AIBJ2(XINT,XAIBJ,XLAMDP,XLAMDH,
     *                      SCR1,SCR2,IDEL,ISYMD,ISYMJ,ISYMAB,
     *                      LUFILE,FNFILE,ANTISYM)
C
C     The original CCSD_AIBJ2 is written by Henrik Koch 27-Mar-1990.
C
C     Purpose: 1) Transformation of the alpha-index to occupied.
C              2) Transformation of the beta-index  to virtual.
C              3) Adding contirubution to result vector and transforming
C                  the last index (delta) to virtual.
C
C     Written by Lilli Irene Ør Kristensen, Winter 2016
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),XAIBJ(*), SCR1(*),SCR2(*)
      DIMENSION XLAMDP(*),XLAMDH(*)
#include "priunit.h"
#include "ccinftap.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      LOGICAL ANTISYM
      CHARACTER*(*) FNFILE
C
C
      CALL QENTER('SO_AIBJ2')
C
      IF (ANTISYM) THEN
         FACDG = -ONE
      ELSE
         FACDG = ONE
      END IF
C
      DO 100 J = 1,NRHF(ISYMJ)
C
         KOFF1 = NNBST(ISYMAB)*(J-1) + 1
C
         IF (ANTISYM) THEN
            CALL CCSD_ASYMSQ(XINT(KOFF1),ISYMAB,SCR1,0,0)
         ELSE
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAB,SCR1)
         END IF
C
C--------------------------------------------------
C  1)  Transformation of the A-index to occupied.
C--------------------------------------------------
C
         KOFF3 = 1
         DO 110 ISYMI = 1,NSYM
C
            ISYMA = ISYMI
            ISYMB = MULD2H(ISYMA,ISYMAB)
C
            KOFF1 = IAODIS(ISYMA,ISYMB) + 1
            KOFF2 = ILMRHF(ISYMI) + 1
C
            NBASA = MAX(NBAS(ISYMA),1)
            NBASB = MAX(NBAS(ISYMB),1)
            CALL DGEMM('T','N',NBAS(ISYMB),NRHF(ISYMI),NBAS(ISYMA),
     *                 FACDG,SCR1(KOFF1),NBASA,XLAMDP(KOFF2),
     *                 NBASA,ZERO,SCR2(KOFF3),NBASB)
C
            KOFF3 = KOFF3 + NBAS(ISYMB)*NRHF(ISYMI)
C
  110    CONTINUE
C
C-------------------------------------------------
C    2) Transformation of the B-index to virtual.
C-------------------------------------------------
C
         KOFF2 = 1
         DO 120 ISYMI = 1,NSYM
C 
            ISYMB = MULD2H(ISYMI,ISYMAB)
            ISYMA = ISYMB
C
            KOFF1 = ILMVIR(ISYMA) + 1
            NBASB = MAX(NBAS(ISYMB),1)
C
            KOFF3 = IT1AM(ISYMA,ISYMI) + 1
            NVIRA = MAX(NVIR(ISYMA),1)
            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMB),
     *                   ONE,XLAMDH(KOFF1),NBASB,SCR2(KOFF2),
     *                   NBASB,ZERO,SCR1(KOFF3),NVIRA)
            KOFF2 = KOFF2 + NBAS(ISYMB)*NRHF(ISYMI)
C
  120    CONTINUE
C
C--------------------------------------------------
C        Add the contribution to the result vector and transform
C           delta to b using DAXPY.
C--------------------------------------------------
C
         ISYMB  = ISYMD
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = ISYMAB
C
         DO 130 B = 1, NVIR(ISYMB)
           NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
           IF (ISYMAI .EQ. ISYMBJ) THEN
               NTOTAI = NBJ
               KOFF1 = NBJ*(NBJ - 1)/2 + 1
           ELSE
               NTOTAI = NT1AM(ISYMAI)
               KOFF1  = NTOTAI*(NBJ - 1) + 1
           ENDIF
           KOFF2 = ILMVIR(ISYMB) + NBAS(ISYMD)*(B-1) + IDEL
     *             - IBAS(ISYMD)
           CALL DAXPY(NTOTAI,XLAMDH(KOFF2),SCR1,1,XAIBJ(KOFF1),1)
  130     CONTINUE
C         END IF
C
  100 CONTINUE
C
      CALL QEXIT('SO_AIBJ2')
C
      RETURN
      END
