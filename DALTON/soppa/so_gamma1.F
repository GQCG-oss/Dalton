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
       WRITE(LUPRI,*)' SO_GAMMA'
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
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
C
            CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMA),
     *                 ONE,XLAMDP(KOFF2),NBASA,SCR1(KOFF3),NBASA,
     *                 ZERO,SCR2(KOFF4),NRHFK)
C           (k gamma| j delta)=lamda(alpha k)(alpha gamme| j delta)
C            WRITE(LUPRI,*) '(K GAMMA | j delta)'
C            CALL OUTPUT(SCR2(KOFF4),1,NBAS(ISYMG),1,NBAS(ISYMA),
C     *                      NBAS(ISYMG),NBAS(ISYMA),2,LUPRI)
C
            KOFF7 = ILMRHF(ISYMI) + 1
            KOFF8 = IMATIJ(ISYMK,ISYMI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                 ONE,SCR2(KOFF4),NRHFK,XLAMDH(KOFF7),NBASG,
     *                 ZERO,SCR3(KOFF8),NRHFK)
C           (K i | j delta)=(k alpha| j delta)*Lamda(gamma i)
C            WRITE(LUPRI,*) '(K i | j delta)'
C            CALL OUTPUT(SCR3(KOFF8),1,NRHF(ISYMK),1,NRHF(ISYMI),
C     *                      NRHF(ISYMK),NRHF(ISYMI),2,LUPRI)
  110    CONTINUE
C
         DO 120 ISYMJ = 1,NSYM
C           IF ISYMJ = ISYMD (line 91) the loop irrelevant
C  When code works check wether it is needed with print statements
C
            ISYMLJ = MULD2H(ISYML,ISYMJ)
            ISYMKI = MULD2H(ISYMLJ,ISYMOP)
            ISYMCI = MULD2H(ISYMJ,ISYMD)
C
            KSCR5 = 1
            KEND1 = KSCR5 + NMATIJ(ISYMKI)
C
            CALL DZERO(WORK(KSCR5),NMATIJ(ISYMKI))
C
            IF (ISYMKI .GT. ISYMLJ) GOTO 120
C
            DO 130 J = 1,NRHF(ISYMJ)
C
C
               IF (ISYMJ .EQ. ISYMD) THEN
                 CALL DAXPY(NMATIJ(ISYMKI),XLAM(J),SCR3,1,
     *                       WORK(KSCR5),1)
               ENDIF
               NLJ = IMATIJ(ISYML,ISYMJ) + NRHF(ISYML)*(J - 1) + L
C
               KKILJ = IGAMMA(ISYMKI,ISYMLJ) + NLJ*(NLJ-1)/2
               DO 150 NKI = 1,NLJ
C
                  KOFF = KSCR5 + NKI - 1
                  NKILJ = KKILJ + NKI
                  GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(KOFF)
C
  150          CONTINUE
C
  130       CONTINUE
  120    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
