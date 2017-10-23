C  /* Deck so_agbd1 */
      SUBROUTINE SO_AGBD1(XINT,OMEGA2,XMGD,WORK,LWORK,IDEL,ISYMD)
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
C                     CALL OUTPUT(WORK(KOFF4),1,NBAS(ISYMA),1,
C     *                      NRHF(ISYMI),NBAS(ISYMA),NRHF(ISYMI),1,LUPRI)
CLilli I don't understand this if statement. And to my understanding
C      NTOTG has to be NBAS(ISYMG)
C                     IF (ISYMG .LT. ISYMD) THEN
C                        NTOTG = NBAS(ISYMG)
C                     ELSE
C                        NTOTG = IDEL - IBAS(ISYMD)
c                     ENDIF
CLEND.
                     NTOTG = NBAS(ISYMG)
C
                      CALL DGEMM('N','N',NBAS(ISYMA),
     *                         NRHF(ISYMI),NTOTG,ONE,WORK(KOFF2),NBASA,
     *                         XMGD(KOFF3),NBASG,ZERO,WORK(KOFF4),NBASA)
c                      WRITE(LUPRI,*) 'XMGD(KOFF3)'
c                      CALL OUTPUT(XMGD(KOFF3),1,NTOTG,1,NRHF(ISYMI),
c     *                      NTOTG,NRHF(ISYMI),1,LUPRI)
c                      WRITE(LUPRI,*) 'WORK(KOFF2)',KOFF2,KSCR1
c                      CALL OUTPUT(WORK(KOFF2),1,NBAS(ISYMA),1,NTOTG,
c     *  NBAS(ISYMA),NTOTG,1,LUPRI)
c                      WRITE(LUPRI,*) 'WORK(KOFF4)',KOFF4,WORK(37),
c     * OMEGA2(31)
c                      CALL OUTPUT(WORK(KOFF4),1,NBAS(ISYMA),1,
c     *                      NRHF(ISYMI),NBAS(ISYMA),NRHF(ISYMI),1,LUPRI)
C                      CALL OUTPUT(XMGD(KOFF3),1,NTOTG,1,NRHF(ISYMI),
C     *                      NTOTG,NRHF(ISYMI),1,LUPRI)
C                      WRITE(LUPRI,*) 'XMGD(KOFF3)',XMGD(KOFF3)
C                      WRITE(LUPRI,*) 'WORK(KOFF4)',WORK(KOFF4)
C
                       DO 150 I = 1,NRHF(ISYMI)
                           DO 160 A = 1,NBAS(ISYMA)
C
                              KOFF5 = KOFF4 + A-1  + (I-1)*NBAS(ISYMA)
C     *                       + (J-1)*NBAS(ISYMA)*NRHF(ISYMI)
C     *                       + (B-1)*NBAS(ISYMA)*NRHF(ISYMI)*NRHF(ISYMJ)
C
                              KOFF6 = A + (I-1)*NBAS(ISYMA)*NBAS(ISYMB)
     *                    + (J-1)*NBAS(ISYMA)*NBAS(ISYMB)*NRHF(ISYMJ)
     *                    + (B-1)*NBAS(ISYMA)
C
c                              WRITE(LUPRI,*) 'KOFF',KOFF5,KOFF6
                            OMEGA2(KOFF6) = OMEGA2(KOFF6) + WORK(KOFF5)
c                            WRITE(LUPRI,*) 'OMEGA2',OMEGA2(KOFF6)
  160                      CONTINUE
  150                   CONTINUE
C
  140             CONTINUE


C       WRITE(LUPRI,*) 'OMEGA2',OMEGA2(217),OMEGA2(218),OMEGA2(219),
C     *                  OMEGA2(220),OMEGA2(221),OMEGA2(222)

CLilli This part of the doesn't make sence
C---------------------------------------
C                2)  Accumulate the result.
C---------------------------------------
C
C                  NBJ = IT1AO(ISYMB,ISYMJ) + NBAS(ISYMB)*(J - 1) + B
C
C                  IF (ISYMAI .EQ. ISYMBJ) THEN
C                     WORK(KSCR2+NBJ-1) = TWO*WORK(KSCR2+NBJ-1)
C                     WRITE(LUPRI,*) 'ACCUMULATE'
C                     WRITE(LUPRI,*) KSCR2+NBJ-1,WORK(KSCR2+NBJ-1)
C                  ENDIF
C
C                  DO 150 NAI = 1,NT1AO(ISYMAI)
C                     KOFF5 = KSCR2 + NAI - 1
C                     NAIBJ = IT2AO(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
C
C                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KOFF5)
C                     WRITE(LUPRI,'(A)') '(ac|bd)'
C                     WRITE(LUPRI,*)'SO_AGBD',NAIBJ,OMEGA2(NAIBJ)
C                     FLUSH(LUPRI)
C  150             CONTINUE
C
  130          CONTINUE
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
C
      RETURN
      END
