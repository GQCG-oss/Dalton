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
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYMO2 = MULD2H(ISYMOP,ISYMTR)
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
               KSCR2 = KSCR1 + NBAS(ISYAL)*NBAS(ISYBE)*NRHF(ISYMI)
     *                         *NRHF(ISYMJ)
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
              CALL DCOPY(NBAS(ISYAL)*NBAS(ISYBE)*NRHF(ISYMI)*NRHF(ISYMJ)
     *                  ,OMEGA2(1),1,WORK(KSCR1),1)
              CALL DZERO(WORK(KSCR2),NBAS(ISYAL)*NVB)
              CALL DZERO(OMEGA2(1),NBAS(ISYAL)*NBAS(ISYBE)*NRHF(ISYMI)
     *                      *NRHF(ISYMJ))
              DO 130 J = 1,NRHF(ISYMJ)
                 DO 140 I = 1,NRHF(ISYMI)
C
C------------------------------------------
C                    Squareup the AB block.
C The OMEGA2 vector is already squared up, but it requires A LOT of
C memory.
C But due to the lack of understading this part of the code, and the 
C last part in so_agbd, it has been commented out. 
C
C It could be nice if someone actually would look in to this, so this
C the SO_2NDDEWF code wouldn't require as much memory.
C------------------------------------------
C
C                     IF (ISYMI .EQ. ISYMJ) THEN
C                        NIJ = IMIJP(ISYMI,ISYMJ) + INDEX(I,J)
C                        FAC1 = ONE
C                        IF (I .GT. J) FAC1 = -ONE
C                     ELSE IF (ISYMI .LT. ISYMJ) THEN
C                        NIJ = IMIJP(ISYMI,ISYMJ)
C     *                      + NRHF(ISYMI)*(J - 1) + I
C                        FAC1 = ONE
C                     ELSE
C                        NIJ = IMIJP(ISYMI,ISYMJ)
C     *                      + NRHF(ISYMJ)*(I - 1) + J
C                        FAC1 = -ONE
C                     ENDIF
CC
C                     DO 166 B = 1,NBAS(ISYBE)
C                        DO 167 A = 1,NBAS(ISYAL)
CC
C                           IF (ISYAL .EQ. ISYBE) THEN
C                              NABP = IAODPK(ISYAL,ISYBE)
C     *                             + INDEX(A,B)
C                              FAC2 = ONE
C                              IF (A .GT. B) FAC2 = -ONE
C                           ELSE IF (ISYAL .LT. ISYBE) THEN
C                              NABP = IAODPK(ISYAL,ISYBE)
C     *                             + NBAS(ISYAL)*(B - 1) + A
C                              FAC2 = ONE
C                           ELSE
C                              NABP = IAODPK(ISYAL,ISYBE)
C     *                             + NBAS(ISYBE)*(A - 1) + B
C                              FAC2 = -ONE
C                           ENDIF
CC
C                           NABIJP = IT2ORT(ISALBE,ISYMIJ)
C     *                            + NNBST(ISALBE)*(NIJ - 1) + NABP
CC
C                           NABIJM = NT2ORT(ISYMOP)
C     *                            + IT2ORT(ISALBE,ISYMIJ)
C     *                            + NNBST(ISALBE)*(NIJ - 1) + NABP
CC
C                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
CC
C                           FAC = FAC1*FAC2
CC
C                           WORK(NAB) =
C     *                       HALF*(OMEGA2(NABIJP) + FAC*OMEGA2(NABIJM))
CC
C  167                   CONTINUE
C  166                CONTINUE
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
                     KOFF = KSCR1 + (I-1)*NBAS(ISYAL)*NBAS(ISYBE)
     *                     + (J-1)*NBAS(ISYAL)*NBAS(ISYBE)*NRHF(ISYMI)
C
                     CALL DGEMM('N','N',NBAS(ISYAL),NVIR(ISYMB),
     *                          NBAS(ISYBE),ONE,WORK(KOFF),NBASA,
     *                          XLAMDP(KOFF1),NBASB,ZERO,WORK(KSCR2),
     *                          NBASA)
C
C
                     KOFF2 = IGLMVI(ISYAL,ISYMA) + 1
C
                     CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                          NBAS(ISYAL),ONE,XLAMPC(KOFF2),NBASA,
     *                          WORK(KSCR2),NBASA,ZERO,WORK(KSCR3),
     *                          NVIRA)
C
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
                           ENDIF
C
C                           RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
                            OMEGA2(NAIBJ) = OMEGA2(NAIBJ)+ WORK(NAB)
C
C                           WRITE(LUPRI,*) 'OMEGA2',NAIBJ,OMEGA2(NAIBJ)
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
