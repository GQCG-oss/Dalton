C  /* Deck so_gamma1 */
      SUBROUTINE SO_GAMMA1(DSRHF,GAMMA,XLAM,
     *              XLAMDP,XLAMDH,SCR1,SCR2,SCR3,WORK,
     *              LWORK,ISYMD,ISYML,ISYMAG)
C
C     Original routine CCRHS_GAM1 is written by Henrik Koch 3-Jan-1994
C
C     Purpose: Calculate the gamma intermediate.
C
C     List of variables:
C       DSRHF:  (alpha gamma| j delta) integrals
C       GAMMA:  Output array - where gamma is placed
C       XLAM:   XLAMDH for given IDEL
C       XLAMDP:
C       XLAMDH:
C       SCR1:   DSRHF squared up is placed here
C       SCR2:   (k gamma| j alpha) is placed here
C       SCR3:   (k i | j delta) is placed here
C       WORK:   Used for storing (ki|j l/delta?)
C       ISYMD:
C       ISYML:  Loop in so_gamma (1,NSYM)
C       ISYMAG: MULD2H(ISYML,ISYDIS)
C
C     Written by Lilli Irene Ør Kristensen, Spring 2017
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION DSRHF(*),GAMMA(*),XLAM(*)
      DIMENSION SCR1(*),SCR2(*),SCR3(*),WORK(*)
      DIMENSION XLAMDP(*),XLAMDH(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      ISYMKC = ISYMAG
C
      DO 100 L = 1,NRHF(ISYML)
C
         KOFF1 = IDSRHF(ISYMAG,ISYML) + NNBST(ISYMAG)*(L - 1) + 1
C
         CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,SCR1)
C
         DO 110 ISYMG = 1,NSYM
C
            ISYMA = MULD2H(ISYMG,ISYMAG)
            ISYMK = ISYMA
            ISYMC = ISYMG
            ISYMI = ISYMG
C
            NBASA = MAX(NBAS(ISYMA),1)
            NBASG = MAX(NBAS(ISYMG),1)
            NRHFK = MAX(NRHF(ISYMK),1)
C
            KOFF2 = ILMRHF(ISYMK) + 1
            KOFF3 = IAODIS(ISYMA,ISYMG) + 1
            KOFF4 = IT1AOT(ISYMK,ISYMG) + 1

            CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMA),
     *                 ONE,XLAMDP(KOFF2),NBASA,SCR1(KOFF3),NBASA,
     *                 ZERO,SCR2(KOFF4),NRHFK)
C
            KOFF7 = ILMRHF(ISYMI) + 1
            KOFF8 = IMATIJ(ISYMK,ISYMI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                 ONE,SCR2(KOFF4),NRHFK,XLAMDH(KOFF7),NBASG,
     *                 ZERO,SCR3(KOFF8),NRHFK)

C
C            WRITE(LUPRI,*) 'CMO'
C            CALL OUTPUT(XLAMDP(KOFF2),1,NBAS(ISYMA),1,NRHF(ISYMK),
C     * NBAS(ISYMA),NRHF(ISYMK),1,LUPRI)
CC
C            WRITE(LUPRI,*) 'SCRM'
C            CALL OUTPUT(SCR1(KOFF3),1,NBAS(ISYMA),1,NBAS(ISYMG),
C     *              NBAS(ISYMA),NBAS(ISYMG),1,LUPRI)
C
C            WRITE(LUPRI,*) '(K GAMMA | j delta)'
C            CALL OUTPUT(SCR2(KOFF4),1,NRHF(ISYMK),1,NBAS(ISYMG),
C     *                      NRHF(ISYMK),NBAS(ISYMG),1,LUPRI)
C
C            WRITE(LUPRI,*) 'CMO'
C            CALL OUTPUT(XLAMDP(KOFF2),1,NBAS(ISYMA),1,NRHF(ISYMK),
C     * NBAS(ISYMA),NRHF(ISYMK),1,LUPRI)
CC
C            WRITE(LUPRI,*) 'SCRM'
C            CALL OUTPUT(SCR1(KOFF3),1,NBAS(ISYMA),1,NBAS(ISYMG),
C     *              NBAS(ISYMA),NBAS(ISYMG),1,LUPRI)
C
C            WRITE(LUPRI,*) 'CMO'
C            CALL OUTPUT(XLAMDP(KOFF7),1,NBAS(ISYMG),1,NRHF(ISYMI),
C     * NBAS(ISYMG),NRHF(ISYMI),1,LUPRI)
C
C           WRITE(LUPRI,*) '(K I | J DELTA)',ISYMG
c           CALL OUTPUT(SCR3(KOFF8),1,NRHF(ISYMK),1,NRHF(ISYMI),
c     *                      NRHF(ISYMK),NRHF(ISYMI),1,LUPRI)
C
  110    CONTINUE
C
         DO 120 ISYMJ = 1,NSYM
C           IF ISYMJ = ISYMD (line 91) the loop irrelevant
cC  When code works check wether it is needed with print statements
cC
            ISYMLJ = MULD2H(ISYML,ISYMJ)
            ISYMKI = MULD2H(ISYMLJ,ISYMOP)
            ISYMCI = MULD2H(ISYMJ,ISYMD)
C
            KSCR5 = 1
            KEND1 = KSCR5 + NMATIJ(ISYMKI)
C
            CALL DZERO(WORK(KSCR5),NMATIJ(ISYMKI))

            IF (ISYMKI .GT. ISYMLJ) GOTO 120
C
            DO 130 J = 1,NRHF(ISYMJ)
C
C
               IF (ISYMJ .EQ. ISYMD) THEN
                 CALL DAXPY(NMATIJ(ISYMKI),XLAM(J),SCR3,1,
     *                       WORK(KSCR5),1)

C                 WRITE(LUPRI,*) 'L',L,'J',J,'ISYMJ',ISYMJ
C                 WRITE(LUPRI,*) 'XLAM(J) ',XLAM(J)
C                 WRITE(LUPRI,*) 'SCR3    '
c                 CALL OUTPUT(SCR3(1),1,NMATIJ(ISYMKI),1,1,
C     *                          NMATIJ(ISYMKI),1,1,LUPRI)
C                 WRITE(LUPRI,*) 'GAMMA',KSCR5
c                CALL OUTPUT(WORK(KSCR5),1,NMATIJ(ISYMKI),1,1,
C     *              NMATIJ(ISYMKI),1,1,LUPRI)


            ENDIF
               NLJ = IMATIJ(ISYML,ISYMJ) + NRHF(ISYML)*(J - 1) + L
C
               KKILJ = IGAMMA(ISYMKI,ISYMLJ) + NLJ*(NLJ-1)/2
               DO 150 NKI = 1,NLJ
C
                  KOFF = KSCR5 + NKI - 1
                  NKILJ = KKILJ + NKI
                  GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(KOFF)
C                  WRITE(LUPRI,*) 'KOFF ',KOFF, 'NKILJ ',NKILJ
C                  WRITE(LUPRI,*) 'WORK ',WORK(KOFF)
  150          CONTINUE
C
  130       CONTINUE
  120    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
