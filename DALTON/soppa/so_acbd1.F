C  /* Deck so_acbd1 */
      SUBROUTINE SO_ACBD1(XINT,OMEGA2,XMGD,WORK,LWORK,IDEL,ISYMD)
C
C     The original routine CCRHS_B1 is written by Henrik Koch 3-Jan-1994
C
C     Purpose: Calculate Sum_(cd) T_(ij)^(cd)*(ac|bd)
C              1) Calculate the contration of amplitude and integral
C              2) Add result to the OMEGA2 vector
C
C     List of variables:
C       XINT:   AO-integrals (alpha gamma| beta delta)
C       OMEGA2: output array
C       XMGD:   T_(ij)^(gamma delta)
C       WORK:
C       ISYMD: Symmetry of Delta
C
C     Written by Lilli Irene Ør Kristensen, Spring 2017
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),OMEGA2(*),XMGD(*)
      DIMENSION WORK(LWORK)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYDIS = MULD2H(ISYMOP,ISYMD)
C
C-----------------------------------------------------------
C      1)  Calculate the contration
C----------------------------------------------------------
C
      DO 100 ISYMB = 1,NSYM
C
         ISYMAG = MULD2H(ISYMB,ISYDIS)
C
         KSCR1 = 1
         KEND1 = KSCR1 + N2BST(ISYMAG)
         LWRK1 = LWORK - KEND1
C
         DO 110 B = 1,NBAS(ISYMB)
C
            KOFF1 = IDSAOG(ISYMB,ISYDIS) + NNBST(ISYMAG)*(B - 1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAG,WORK(KSCR1))
C
            DO 120 ISYMJ = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYMOP)
               ISYMGI = MULD2H(ISYMJ,ISYMD)
C
               KSCR2 = KEND1
               KEND2 = KSCR2 + NT1AO(ISYMAI)
               LWRK2 = LWORK - KEND2
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  CALL DZERO(WORK(KSCR2),NT1AO(ISYMAI))
C
                  DO 140 ISYMI = 1,NSYM
C
                     ISYMG  = MULD2H(ISYMI,ISYMGI)
C
                     IF (ISYMG .GT. ISYMD) GOTO 140
C
                     ISYMA  = MULD2H(ISYMG,ISYMAG)
C
                     KOFF2 = KSCR1 + IAODIS(ISYMA,ISYMG)
                     KOFF3 = IT2BGD(ISYMGI,ISYMJ)
     *                     + NT1AO(ISYMGI)*(J - 1) + 1
     *                     + IT1AO(ISYMG,ISYMI)
                     KOFF4 = KSCR2 + IT1AO(ISYMA,ISYMI)
C
                     NBASA = MAX(NBAS(ISYMA),1)
                     NBASG = MAX(NBAS(ISYMG),1)
C
                     IF (ISYMG .LT. ISYMD) THEN
                        NTOTG = NBAS(ISYMG)
                     ELSE
                        NTOTG = IDEL - IBAS(ISYMD)
                     ENDIF
C
                     CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),
     *                          NTOTG,ONE,WORK(KOFF2),NBASA,
     *                          XMGD(KOFF3),NBASG,ZERO,WORK(KOFF4),
     *                          NBASA)
C                   = Sum_(beta) Sum_(ij)(alpha gamma| beta delta)* T_(ij)^(gamma delta) 
  140             CONTINUE
C
C---------------------------------------
C                2)  Accumulate the result.
C---------------------------------------
C
                  NBJ = IT1AO(ISYMB,ISYMJ) + NBAS(ISYMB)*(J - 1) + B
C
                  IF (ISYMAI .EQ. ISYMBJ) THEN
                     WORK(KSCR2+NBJ-1) = TWO*WORK(KSCR2+NBJ-1)
                  ENDIF
C
                  DO 150 NAI = 1,NT1AO(ISYMAI)
                     KOFF5 = KSCR2 + NAI - 1
                     NAIBJ = IT2AO(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KOFF5)
  150             CONTINUE
C
  130          CONTINUE
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
C
      END
