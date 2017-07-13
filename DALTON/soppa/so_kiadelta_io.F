C  /* Deck so_kiadelta_io */
      SUBROUTINE SO_KIADELTA_IO(OMEGA2,T2AM,XLAMDH,WORK,LWORK,
     *                     ISYVEC,ISYCIM,LUC,CFIL,IV,IOPT)
C
C     asm 17-aug-1994
C
C     Ove Christiansen 30-7-1995: modified to account for general
C                                 non. total symmetric vectors (ISYVEC) and
C                                 intermediates (ISYCIM). LUC and CFIL is
C                                 used to control from which file the
C                                 intermediate is obtained.
C
C                                 if iopt = 1 the C intermediate is assumed
C                                    to be as in energy calc.
C
C                                 if iopt ne. 1 we use the intermediate
C                                    on luc with address given according to
C                                    transformed vector nr iv (iv is not 1
C                                    if several vectors are transformed
C                                    simultaneously.)
C
C                                 in energy calc: iv=1,iopt=1
C
C     PURPOSE:
C             Calculate the C-term making I/O
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION OMEGA2(*),T2AM(*),XLAMDH(*)
      DIMENSION WORK(LWORK)
      CHARACTER CFIL*(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccsdio.h"
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C      IF (OMEGSQ) THEN
C         WRITE(LUPRI,*)
C     &        'I/O in C-term not implemented for square Omega2'
C         CALL QUIT('OMEGSQ = .TRUE.  in CCRHS_CIO')
C      END IF
C
      ISAIBJ = MULD2H(ISYVEC,ISYCIM)
C
      DO 100 ISYMAI = 1,NSYM
C
         IF (NT1AM(ISYMAI) .EQ. 0) GOTO 100
C
         ISYMBJ = MULD2H(ISYMAI,ISAIBJ)
         ISYMCK = MULD2H(ISYVEC,ISYMBJ)
         ISYMDK = MULD2H(ISYCIM,ISYMAI)
C
C------------------------
C        Batch structure.
C------------------------
C
         NT1AI = NT1AM(ISYMAI)
C
         LENAI  = NT1AO(ISYMDK)
         LENMIN = 2*LENAI
         IF (LENMIN .EQ. 0) GOTO 100
C
         NDISAI = LWORK / LENMIN
         IF (NDISAI .LT. 1) THEN
            CALL QUIT('Insufficient space for '//
     &           'allocation in CCRHS_CIO-1')
         END IF
         NDISAI = MIN(NDISAI,NT1AI)
C
         NBATAI = (NT1AI - 1) / NDISAI + 1
C
C--------------------------
C        Loop over batches.
C--------------------------
C
         ILSTAI = 0
         DO 110 IBATAI = 1,NBATAI
C
            IFSTAI = ILSTAI + 1
            ILSTAI = ILSTAI + NDISAI
            IF (ILSTAI .GT. NT1AI) THEN
               ILSTAI = NT1AI
               NDISAI = ILSTAI - IFSTAI + 1
            END IF
C
C-----------------------------
C           Memory allocation.
C-----------------------------
C
            KSCR1 = 1
            KSCR2 = KSCR1 + NDISAI*NT1AO(ISYMDK)
            KEND  = KSCR2 + NDISAI*NT1AO(ISYMDK)
            LWRK1 = LWORK - KEND
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient space for '//
     &              'allocation in CCRHS_CIO-2')
            END IF
C
C----------------------------------
C           Construct P(del k,#ai).
C----------------------------------
C
            KOFF1 = KSCR1
            DO 120 ISYDEL = 1,NSYM
C
               ISYMK = MULD2H(ISYDEL,ISYMDK)
C
               DO 130 IDELTA = 1,NBAS(ISYDEL)
C
                  ID = IDELTA + IBAS(ISYDEL)
C
                  IF (IOPT .EQ. 1 ) THEN
                     IOFF = IT2DEL(ID) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ELSE
                     IOFF = IT2DLR(ID,IV) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ENDIF
C
                  LEN  = NDISAI*NRHF(ISYMK)
C
                  IF (LEN .GT. 0) THEN
                     CALL GETWA2(LUC,CFIL,WORK(KOFF1),IOFF,LEN)
                  ENDIF
C
                  DO 140 NAI = IFSTAI,ILSTAI
C
                     KAI = NAI - IFSTAI + 1
C
                     KOFF2 = KOFF1 + NRHF(ISYMK)*(KAI - 1)
                     KOFF3 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                     + IT1AO(ISYDEL,ISYMK) + IDELTA - 1
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF2),1,WORK(KOFF3),
     *                          NBAS(ISYDEL))
C
  140             CONTINUE
C
                  KOFF1 = KOFF1 + LEN
C
  130          CONTINUE
  120       CONTINUE
C
C-----------------------------------------
C              Transform delta index to c.
C-----------------------------------------
C
            DO 150 NAI = IFSTAI,ILSTAI
C
               KAI = NAI - IFSTAI + 1
C
               DO 160 ISYMC = 1,NSYM
C
                  ISYDEL = ISYMC
                  ISYMK  = MULD2H(ISYMC,ISYMCK)
C
                  NBASD = MAX(NBAS(ISYDEL),1)
                  NVIRC = MAX(NVIR(ISYMC),1)
C
                  KOFF4 = ILMVIR(ISYDEL) + 1
                  KOFF5 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                  + IT1AO(ISYDEL,ISYMK)
                  KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI - 1)
     *                  + IT1AM(ISYMC,ISYMK)
C
                  CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     *                       NBAS(ISYDEL),ONE,XLAMDH(KOFF4),NBASD,
     *                       WORK(KOFF5),NBASD,ZERO,WORK(KOFF6),
     *                       NVIRC)
C                  WRITE(LUPRI,*) 'CMO'
C                  CALL OUTPUT(XLAMDH(KOFF4),1,NBAS(ISYDEL),1,NVIR(ISYMC)
C     *  ,NBAS(ISYDEL),NVIR(ISYMC),1,LUPRI)
C

C                  WRITE(LUPRI,*) '(KI|A DELTA)'
C                  CALL OUTPUT(WORK(KOFF5),1,6,1,5
C     *  ,6,5,1,LUPRI)
C
c                  WRITE(LUPRI,*) '(KI|A C)'
c                  CALL OUTPUT(WORK(KOFF6),1,NVIR(ISYMC),1,NRHF(ISYMK),
C     *                      NVIR(ISYMC),NRHF(ISYMK),1,LUPRI)
C
  160          CONTINUE
  150       CONTINUE
C
C--------------------------------------------
C           Contract P(ck,#ai) with T(bj,ck).
C--------------------------------------------
C
            NT1BJ = MAX(NT1AM(ISYMBJ),1)
            NT1CK = MAX(NT1AM(ISYMCK),1)
C
            KOFF7 = IT2SQ(ISYMBJ,ISYMCK) + 1
C
            CALL DGEMM('N','N',NT1AM(ISYMBJ),NDISAI,NT1AM(ISYMCK),
     *                 ONE,T2AM(KOFF7),NT1BJ,WORK(KSCR1),NT1CK,
     *                 ZERO,WORK(KSCR2),NT1BJ)
C
C            write(lupri,*) 't2am'
C            CALL OUTPUT(T2AM(KOFF7),1,NT1AM(ISYMBJ),1,NT1AM(ISYMCK),
C     *          NT1AM(ISYMBJ),NT1AM(ISYMCK),1,LUPRI)
C            write(lupri,*) 'kiac'
C            CALL OUTPUT(WORK(KSCR1),1,NDISAI,1,NT1AM(ISYMCK),
C     *                      NDISAI,NT1AM(ISYMCK),1,LUPRI)
C            write(lupri,*) 'kscr2'
C            CALL OUTPUT(WORK(KSCR2),1,NT1AM(ISYMBJ),1,NDISAI,
C     *                       NT1AM(ISYMBJ),NDISAI,1,LUPRI)
C------------------------------
C           Scale the diagonal.
C------------------------------
C
            IF (ISYMBJ .EQ. ISYMAI) THEN
C
               DO 170 NAI = IFSTAI,ILSTAI
                  KOFF8 = KSCR2 + NT1AM(ISYMBJ)*(NAI-IFSTAI) + NAI - 1
                  WORK(KOFF8) = TWO * WORK(KOFF8)
  170          CONTINUE
C
            END IF
C
C-----------------------------------------------
C           Add the result to the packed omega2.
C-----------------------------------------------
C
            DO 180 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 190 I = 1,NRHF(ISYMI)
C
                  DO 200 A = 1,NVIR(ISYMA)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                     IF ((NAI .LT. IFSTAI) .OR. (NAI .GT. ILSTAI))
     *                  GOTO 200
C
                     DO 210 ISYMJ = 1,NSYM
C
                        ISYMB  = MULD2H(ISYMJ,ISYMBJ)
                        ISYMAJ = MULD2H(ISYMA,ISYMJ)
                        ISYMBI = MULD2H(ISYMB,ISYMI)
C
                        DO 220 J = 1,NRHF(ISYMJ)
C
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J-1) + A
C
                           CALL CC_PUTC(WORK(KSCR2),OMEGA2,ISYMAI,
     *                                  ISYMAJ,ISYMBI,ISYMBJ,ISYMB,
     *                                  ISYMI,ISYMJ,NAI,NAJ,I,J,
     *                                  IFSTAI)
C
  220                   CONTINUE
  210                CONTINUE
  200             CONTINUE
  190          CONTINUE
  180       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
