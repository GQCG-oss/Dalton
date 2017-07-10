C   /* Deck so_agbd */
      SUBROUTINE SO_AGBD(XINT,OMEGA2,XLAMDP,XLAMDH,SCRM,
     *                   WORK,LWORK,IDEL,ISYMD)
C
C     Original rutine ccrhs_b is written by Henrik Koch 3-Jan-1994
C
C     Purpose: Calculate Sum_(cd) T_(ij)^(cd)(ac|bd)
C               1) Transform the c index in T_(ij)^(c delta) to gammma.
C               2) call routine to calculate the contration of
C                   (alpha gamma| beta delta) T_(ij)^(gamma delta)
C
C     List of Variables:
C       XINT:   AO-integrals (alpha gamma | beta delta)
C       OMEGA2: ouput array
C       XLAMDP:
C       XLAMDH:
C       SCRM:   T_(ij)^(c delta)
C       IDEL:
C       ISYMD:  Symmetry of delta
C
C     Written by Lilli Irene Ør Kristensen, Spring 2017
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),OMEGA2(*),XLAMDH(*),WORK(LWORK)
      DIMENSION XLAMDP(*),SCRM(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KMGD   = 1
      KEND1  = KMGD   + NT2BGD(ISYMD)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_B')
      ENDIF
C
C-----------------------------
C     Prepare the data arrays.
C-----------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMCI = MULD2H(ISYMJ,ISYMD)
         ISYMGI = ISYMCI
C
         DO 110 ISYMI = 1,NSYM
C
            ISYMC = MULD2H(ISYMI,ISYMCI)
            ISYMG = ISYMC
C
            NVIRC = MAX(NVIR(ISYMC),1)
            NBASG = MAX(NBAS(ISYMG),1)
C
            KOFF1 = ILMVIR(ISYMC) + 1
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               KOFF2 = IT2BCD(ISYMCI,ISYMJ) + IT1AM(ISYMC,ISYMI)
     *               + NT1AM(ISYMCI)*(J - 1) + 1
               KOFF3 = IT2BGD(ISYMGI,ISYMJ) + IT1AO(ISYMG,ISYMI)
     *               + NT1AO(ISYMGI)*(J - 1) + 1
C
C------------------------------------------------
C            1) Transfrom the c index to gamma.
C-------------------------------------------------
C
               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NVIR(ISYMC),
     *                    ONE,XLAMDH(KOFF1),NBASG,SCRM(KOFF2),NVIRC,
     *                    ZERO,WORK(KOFF3),NBASG)
C                Sum_(ij) Lamda(gamma c)*T_(ij)^(c delta)
C
               LWOR3 = NBAS(ISYMG)*NRHF(ISYMI)
C               WRITE(LUPRI,*)' T_(ij)^(a delta)', KOFF2
C               CALL OUTPUT(SCRM(KOFF2),1,NBAS(ISYMG)*NRHF(ISYMI),1,1,
C     *                      NBAS(ISYMG)*NRHF(ISYMI),1,1,LUPRI)
C               WRITE(LUPRI,*)'CMO', KOFF1
C               CALL OUTPUT(XLAMDH(KOFF1),1,NRHF(ISYMI)*NVIR(ISYMC),1,1,
C     *                      NRHF(ISYMI)*NVIR(ISYMC),1,1,LUPRI)
C               WRITE(LUPRI,*) 'WORK',KOFF3
C               CALL OUTPUT(WORK(KOFF3),1,LWOR3,1,1,
C     *                      LWOR3,1,1,LUPRI)
C
C
C               IF (ISYMG .EQ. ISYMD) THEN
C                  KOFF4 = KOFF3 + IDEL - IBAS(ISYMD) - 1
C                  CALL DSCAL(NRHF(ISYMI),HALF,WORK(KOFF4),NBAS(ISYMG))
C                  CALL DCOPY(NRHF(ISYMI),WORK(KOFF3),1,WORK(KOFF4),1)
C                  WRITE(LUPRI,*)' KOFF4', KOFF4
C                  CALL OUTPUT(WORK(KOFF4),1,LWOR3,1,1,LWOR3,1,1,LUPRI)
C               ENDIF
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
C      WRITE(LUPRI,*) ' T_(ij)^(gamma delta)'
C      CALL OUTPUT(WORK(1),1,150,1,1,150,1,1,LUPRI)
C--------------------------------------------------------------------
C     2) Call routine to calculate the contraction.
C---------------------------------------------------------------------
C      WRITE(LUPRI,*) 'KMGD',(WORK(KMGD+I),I=0,NT2BGD(ISYMD)-1)
      FLUSH(LUPRI)
C
      CALL SO_AGBD1(XINT,OMEGA2,WORK(KMGD),WORK(KEND1),LWRK1,IDEL,ISYMD)
C
      RETURN
      END
