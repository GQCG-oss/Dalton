C  /* Deck so_t2mo */
      SUBROUTINE SO_T2MO(OMEGA2,RHO2,XLAMDP,
     *                   XLAMPC,ISYMTR,WORK,LWORK)
C
C     The Origial routine CC_T2MO is written by Henrik Koch and Alfredo Sanchez 15-July-1994
C
C     Purpose: Transform the contribution in the  Omega2 vector from the AO basis to the MO
C               basis.
C
C     List of variables:
C       OMEGA2: Containing the contribution to the 2nd order DE coef.
C       RHO2:   Output array
C       XLAMDP: MO-coeffecients
C       XLAMPC: MO-coeffecients
C       ISYMTR: Symmetry of trialvector = 1
C
C     Written by Lilli Irene Ør Kristensen, Spring 2017
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION OMEGA2(*),RHO2(*),
     *          XLAMDP(*), WORK(LWORK), XLAMPC(*)
#include "ccorb.h"
#include "ccsdsym.h"
#include "symsq.h"
C#include "cclr.h"
C
C      ISYMOP = 1
C      ISYMTR = 1
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYMO2 = MULD2H(ISYMOP,ISYMTR)
C      ISYMO1 = MULD2H(ISYMO2,ISYMC2)
C
      CALL DZERO(RHO2,NT2AM(ISYMO2))
C
      DO 100 ISYMJ = 1,NSYM
         DO 110 ISYMI = 1,NSYM
C
            ISYMIJ = MULD2H(ISYMI,ISYMJ)
            ISALBE = MULD2H(ISYMIJ,ISYMOP)
            ISYMAB = MULD2H(ISYMIJ,ISYMO2)
C
            DO 120 ISYBE = 1,NSYM
C
               ISYAL  = MULD2H(ISYBE,ISALBE)
               ISYALI = MULD2H(ISYAL,ISYMI)
               ISYBEJ = MULD2H(ISYBE,ISYMJ)
C
C-----------------------------------------------
C              Dynamic allocation of work space.
C-----------------------------------------------
C
               ISYMA = MULD2H(ISYAL,ISYMTR)
               NVA = MAX(NVIR(ISYMA),NVIR(ISYAL))
               NRA = MAX(NRHF(ISYMA),NRHF(ISYAL))
               ISYMB = MULD2H(ISYBE,ISYMTR)
               NVB = MAX(NVIR(ISYMB),NVIR(ISYBE),NRHF(ISYBE))
               NRB = MAX(NRHF(ISYMB),NRHF(ISYBE))
C
               KSCR1 = 1
               KSCR2 = KSCR1 + NBAS(ISYAL)*NBAS(ISYBE)
               KSCR3 = KSCR2 + NBAS(ISYAL)*NVB
               KSCR4 = KSCR3 + NVA*NVB
               KSCR5 = KSCR4 + NBAS(ISYAL)*NRB
               KEND1 = KSCR5 + NRA*NRB
               LWRK1 = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Not enough space in SO_T2MO')
               END IF
C
               DO 130 J = 1,NRHF(ISYMJ)
                  DO 140 I = 1,NRHF(ISYMI)
C
C------------------------------------------
C                    Squareup the AB block.
C------------------------------------------
C
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NIJ = IMIJP(ISYMI,ISYMJ) + INDEX(I,J)
                        FAC1 = ONE
                        IF (I .GT. J) FAC1 = -ONE
                     ELSE IF (ISYMI .LT. ISYMJ) THEN
                        NIJ = IMIJP(ISYMI,ISYMJ)
     *                      + NRHF(ISYMI)*(J - 1) + I
                        FAC1 = ONE
                     ELSE
                        NIJ = IMIJP(ISYMI,ISYMJ)
     *                      + NRHF(ISYMJ)*(I - 1) + J
                        FAC1 = -ONE
                     ENDIF
C
                     DO 166 B = 1,NBAS(ISYBE)
                        DO 167 A = 1,NBAS(ISYAL)
C
                           IF (ISYAL .EQ. ISYBE) THEN
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + INDEX(A,B)
                              FAC2 = ONE
                              IF (A .GT. B) FAC2 = -ONE
                           ELSE IF (ISYAL .LT. ISYBE) THEN
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + NBAS(ISYAL)*(B - 1) + A
                              FAC2 = ONE
                           ELSE
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + NBAS(ISYBE)*(A - 1) + B
                              FAC2 = -ONE
                           ENDIF
C
                           NABIJP = IT2ORT(ISALBE,ISYMIJ)
     *                            + NNBST(ISALBE)*(NIJ - 1) + NABP
C
                           NABIJM = NT2ORT(ISYMOP)
     *                            + IT2ORT(ISALBE,ISYMIJ)
     *                            + NNBST(ISALBE)*(NIJ - 1) + NABP
C
                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
C
                           FAC = FAC1*FAC2
C
                           WORK(NAB) =
     *                       HALF*(OMEGA2(NABIJP) + FAC*OMEGA2(NABIJM))
C
  167                   CONTINUE
  166                CONTINUE
C
C
C------------------------------------------------------------
C                    Transform the AB block to virtual space.
C------------------------------------------------------------
C
C
                     ISYMA = MULD2H(ISYAL,ISYMTR)
                     ISYMB = ISYBE
                     ISYMAI = MULD2H(ISYMA,ISYMI)
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NBASB = MAX(NBAS(ISYBE),1)
                     NVIRA = MAX(NVIR(ISYMA),1)
C
                     KOFF1 = ILMVIR(ISYBE) + 1
C
                     CALL DGEMM('N','N',NBAS(ISYAL),NVIR(ISYMB),
     *                          NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                          XLAMDP(KOFF1),NBASB,ZERO,WORK(KSCR2),
     *                          NBASA)
C
                     KOFF2 = IGLMVI(ISYAL,ISYMA) + 1
C
                     CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                          NBAS(ISYAL),ONE,XLAMPC(KOFF2),NBASA,
     *                          WORK(KSCR2),NBASA,ZERO,WORK(KSCR3),
     *                          NVIRA)
C
C                     WRITE(LUPRI,*) 'SO_T2MO'
C                     CALL OUTPUT(WORK(KSCR3),1,NVIR(ISYMA),1,NVIR(ISYMB)
C     *               ,NVIR(ISYMA),NVIR(ISYMB),2,LUPRI)
C--------------------------------------------
C                    Store the omega2 vector.
C--------------------------------------------
C
                     DO 170 B = 1,NVIR(ISYMB)
                        NBJ   = IT1AM(ISYMB,ISYMJ)
     *                        + NVIR(ISYMB)*(J-1) + B
                        DO 180 A = 1,NVIR(ISYMA)
C
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1) + A
                           NAB   = KSCR3 + NVIR(ISYMA)*(B - 1) + A - 1
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
C
                              IF (NAI .GT. NBJ) GOTO 180
C
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSEIF (ISYMAI .LT. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           ELSEIF (ISYMAI .GT. ISYMBJ) THEN
                              GOTO 180
chjaaj: next two lines are commented because it is dead code
c                             NAIBJ = IT2AM(ISYMAI,ISYMBJ)
c    *                              + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                           ENDIF
C
                           RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
C
                           WRITE(LUPRI,*) 'NAIBJ',NAIBJ
                           WRITE(LUPRI,*) 'RHO2(NAIBJ)',RHO2(NAIBJ)
                           FLUSH(LUPRI)
  180                   CONTINUE
  170                CONTINUE
C
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
