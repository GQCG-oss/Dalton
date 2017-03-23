C  /* Deck so_t2mo */
      SUBROUTINE SO_T2MO(OMEGA2,RHO2,GAMMA,XLAMDP,
     *                   XLAMPC,ISYMTR,WORK,LWORK,ISYMOP)
C
C     The Origial routine CC_T2MO is written by Henrik Koch and Alfredo Sanchez 15-July-1994
C
C     Purpose: Transform the Omega2 vector from the AO basis to the MO
C               basis.
C
C     List of variables:
C       OMEGA2: 
C       RHO2:   Output array
C       GAMMA:  
C       XLAMDP
C       XLAMDPC
C       ISYMTR: Symmetry of trialvector = 1
C       ISYMOP: Symmetry of operator. Is 1 in AO-SOPPA
C
C     Written by Lilli Irene Ør Kristensen, Spring 2017
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION OMEGA2(*), RHO2(*), GAMMA(*),
     *          XLAMDP(*), WORK(*), XLAMPC(*)
#include "ccorb.h"
#include "ccsdsym.h"
#include "symsq.h"
#include "cclr.h"
C
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
                  CALL QUIT('Not enough space in CC_T2MO')
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
                     ENDIF
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
C                           ENDIF
C
                           RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
C
  180                   CONTINUE
  170                CONTINUE
C
C
C-------------------------------------------------------------
CDo i really need this?     Transform the AB block to occupied space.
C-------------------------------------------------------------
C
C                     NBASA = MAX(NBAS(ISYAL),1)
C                     NBASB = MAX(NBAS(ISYBE),1)
C                     NRHFA1 = MAX(NRHF(ISYAL),1)
C
C                     KOFF1 = ILMRHF(ISYBE) + 1
C
C                     CALL DGEMM('N','N',NBAS(ISYAL),NRHF(ISYBE),
C     *                          NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
C     *                          XLAMDP(KOFF1),NBASB,ZERO,WORK(KSCR4),
C     *                          NBASA)
C
C
C                     KOFF2 = ILMRHF(ISYAL) + 1
C
C                     CALL DGEMM('T','N',NRHF(ISYAL),NRHF(ISYBE),
C     *                          NBAS(ISYAL),ONE,XLAMDP(KOFF2),NBASA,
C     *                          WORK(KSCR4),NBASA,ZERO,WORK(KSCR5),
C     *                          NRHFA1)
C
C-------------------------------------------
C                    Store the gamma matrix.
C-------------------------------------------
C
C                     ISYMK = ISYAL
C                     ISYML = ISYBE
CC
C                     ISYMKI = MULD2H(ISYMK,ISYMI)
C                     ISYMLJ = MULD2H(ISYML,ISYMJ)
CC
C                     DO 190 L = 1,NRHF(ISYML)
CC
C                        NLJ = IMATIJ(ISYML,ISYMJ)
C     *                      + NRHF(ISYML)*(J - 1) + L
CC
C                        DO 200 K = 1,NRHF(ISYMK)
CC
C                           NKL = KSCR5 + NRHF(ISYMK)*(L - 1) + K - 1
CC
C                           NKI = IMATIJ(ISYMK,ISYMI)
C     *                         + NRHF(ISYMK)*(I - 1) + K
CC
C                           IF (ISYMKI .EQ. ISYMLJ) THEN
C                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
C     *                              + INDEX(NKI,NLJ)
C                              GAMMA(NKILJ) = WORK(NKL)
C                           ELSE IF (ISYMKI .LT. ISYMLJ) THEN
C                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
C     *                              + NMATIJ(ISYMKI)*(NLJ - 1) + NKI
C                              GAMMA(NKILJ) = WORK(NKL)
C                           ENDIF
CC
C  200                   CONTINUE
C  190                CONTINUE
C
C
  999                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
