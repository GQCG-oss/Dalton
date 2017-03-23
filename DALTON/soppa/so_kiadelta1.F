C /* Deck so_kiadelta1 */
      SUBROUTINE SO_KIADELTA1(XINT,ISYMT2,SCR1,SCR2,
     *                    XLAMDP,XLAMDH,XLAMPC,ISYMPC,
     *                    WORK,LWORK,ISYDIS,IDEL,FACTC,LUC,CFIL)
C
C     The original routine CCRHS_C1 is written by Henrik Koch 3-Jan-1994
C
C     Porpose: Calculate (ki|a delta),
C               where the delta-loop is outside in SO_2NDDEWF
C
C     List of variables:
C       XINT:   AO-integrals
C       ISYMT2: ISYMTR
C       SCR1:   Where (ki|a delta) is placed
C       SCR2:   Where we caculte (ki|alpha delta)
C       XLAMDH: 
C       XLAMPC: 
C       ISYMPC: 
C       ISYDIS: 
C       IDEL:   
C       FACTC:  HALF
C       LUC:
C       CFIL:
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "ccsdinp.h"
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,XMHALF=-0.5D0)
      PARAMETER (TWO=2.0D0)
      DIMENSION XINT(*)
      DIMENSION SCR1(*), SCR2(*),
      DIMENSION XLAMDH(*),XLAMPC(*)
      DIMENSION WORK(LWORK)
      CHARACTER CFIL*(*)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "r12int.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
      ISAIK2 = MULD2H(ISYDIS,ISYMT2)
      IF (ISYAIK .NE. ISAIK2) CALL QUIT('Symmetry mismatch in CCRHS_C')
C
C
C----------------------------------------------------------
C     Transform the two occ. indices to MO-basis: (k i | alfa delta) 
C----------------------------------------------------------
C
      DO 300 ISYMA = 1,NSYM
C
         ISYMBG = MULD2H(ISYMA,ISYDIS)
C
         KSCR10 = 1
         KEND1  = KSCR10 + N2BST(ISYMBG)
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Not enough space for allocation in CCRHS_C1')
         END IF
C
         DO 310 A = 1,NBAS(ISYMA)
C
            KOFF1 = IDSAOG(ISYMA,ISYDIS) + NNBST(ISYMBG)*(A - 1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMBG,WORK(KSCR10))
C                   XINT -> WORK(KSCR10)
C
            DO 320 ISYMG = 1,NSYM
C
               ISYMI  = ISYMG
               ISYMB  = MULD2H(ISYMG,ISYMBG)
               ISYMK  = ISYMB
               ISYMAI = MULD2H(ISYMA,ISYMI)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NBASG = MAX(NBAS(ISYMG),1)
               NRHFK = MAX(NRHF(ISYMK),1)
C
               KSCR11 = KEND1
               KSCR12 = KSCR11 + NRHF(ISYMK)*NBAS(ISYMG)
               KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
               LWRK2  = LWORK  - KEND2
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_C1')
               END IF
C
               KOFF2 = ILMRHF(ISYMK) + 1
               KOFF3 = KSCR10 + IAODIS(ISYMB,ISYMG)
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMB),
     *                    ONE,XLAMDP(KOFF2),NBASB,WORK(KOFF3),NBASB,
     *                    ZERO,WORK(KSCR11),NRHFK)
C                  (k beta | gamma delta) here (K gamma |alpha beta)
C
               KOFF5 = ILMRHF(ISYMI) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KSCR11),NRHFK,XLAMDH(KOFF5),NBASG,
     *                    ZERO,WORK(KSCR12),NRHFK)
C                   (k i| gamma delta ) here ( k i| alpha beta)
C
               DO 330 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                  KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
                  KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                  CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR2(KOFF8),1)
C                               WORK(KSCR12)->SCR2(KOFF8)
  330          CONTINUE
C
  320       CONTINUE
C
  310    CONTINUE
C
  300 CONTINUE
C
C-----------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C-----------------------------------------------
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
C
      DO 810 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYAIK)
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 820 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMAI)
            ISYMAL= MULD2H(ISYMPC,ISYMA)
            ISYALI= MULD2H(ISYMAL,ISYMI)
            NBASAL = MAX(NBAS(ISYMAL),1)
C
            DO 830 I = 1,NRHF(ISYMI)
C
               NAI  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
               MALI = IT1AO(ISYMAL,ISYMI) + NBAS(ISYMAL)*(I - 1) + 1
C
               KOFF1 = IT2BGT(ISYMK,ISYALI) + NRHF(ISYMK)*(MALI- 1) + 1
               KOFF2 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF3 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),NBAS(ISYMAL),
     *                    ONE,SCR2(KOFF1),NRHFK,XLAMPC(KOFF2),NBASAL,
     *                    FACTC,SCR1(KOFF3),NRHFK)
C
  830       CONTINUE
  820    CONTINUE
  810 CONTINUE
C
C---------------------------------------------------------
C     Dump to disk the new contribution.
C---------------------------------------------------------
C
         IOFF = IT2DEL(IDEL) + 1
C
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         CALL PUTWA2(LUC,CFIL,SCR1,IOFF,NT2BCD(ISYAIK))
      ENDIF
C
  999 CONTINUE
C
      RETURN
      END
