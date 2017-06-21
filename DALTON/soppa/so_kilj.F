C  /* Deck so_kilj */
      SUBROUTINE SO_KILJ(OMEGA2,T2AM,GAMMA,WORK,LWORK,ISYGAM,ISYVEC)
C
C     Origanl routine CCRHS_A is written by Henrik Koch & Ove Christiansen 20-Jan-1994
C         Generalised to non. total sym gamma (isygam) og non. tot. sym
C         double excitation vector (isyvec) Ove Christiansen 29-7-1995
C
C     Purpose: Calculate Sum_(kl) T_(kl)^(ab)*(ki|lj)
C               Contract T_(kl)^(cd) with (ki|lj)
C
C     List of variables:
C       OMEGA2: Contains all other contributions for 2nd de WF
C       T2AM:   T2 amplitudes
C       GAMMA:  Contains the gamma intermediate (ki|lj)
C       WORK:
C       ISYGAM: ISYMTR = 1
C       ISYVEC: ISYMTR = 1
C
C     Written by Lilli Ireme Ør Kristensen, Spring 2017
C
C
#include "implicit.h"
      PARAMETER(ZERO=0.0D0, ONE=1.0D0)
      DIMENSION OMEGA2(*),GAMMA(*),T2AM(*),WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C------------------------------------------------------------
C     Calculate the contraction, and add to OMEGA2 vector
C------------------------------------------------------------
C
      ISAIBJ = MULD2H(ISYGAM,ISYVEC)
C
      DO 100 ISYMLJ = 1,NSYM
C
         ISYMKI = MULD2H(ISYMLJ,ISYGAM)
C
         KSCR1 = 1
         KEND1 = KSCR1 + NMATIJ(ISYMKI)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space for allocation in SO_KILJ')
         END IF
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYML = MULD2H(ISYMJ,ISYMLJ)
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               DO 130 L = 1,NRHF(ISYML)
C
C
                   NLJ = IMATIJ(ISYML,ISYMJ)
     *                   + NRHF(ISYML)*(J - 1) + L
C
                  DO 140 ISYMK = 1,NSYM
C
                     ISYMI = MULD2H(ISYMK,ISYMKI)
C
                     DO 150 I = 1,NRHF(ISYMI)
C
                        DO 160 K = 1,NRHF(ISYMK)
C
                           NKI = IMATIJ(ISYMK,ISYMI)
     *                            + NRHF(ISYMK)*(I - 1) + K
C
                           IF (ISYMKI .EQ. ISYMLJ) THEN
                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                              + INDEX(NKI,NLJ)
                           ELSE
                              IF (ISYMKI .LT. ISYMLJ) THEN
                                 NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                                 + NMATIJ(ISYMKI)*(NLJ - 1) + NKI
                              ELSE
                                 NKILJ = IGAMMA(ISYMLJ,ISYMKI)
     *                                 + NMATIJ(ISYMLJ)*(NKI - 1) + NLJ
                              ENDIF
                           ENDIF
C
                           NSTO = IMATIJ(ISYMK,ISYMI)
     *                          + NRHF(ISYMK)*(I - 1) + K
C
                           WORK(KSCR1 + NSTO - 1) = GAMMA(NKILJ)
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
C
                  DO 170 ISYMB = 1,NSYM
C
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
                     ISYMAI = MULD2H(ISYMBJ,ISAIBJ)
                     ISYMBL = MULD2H(ISYMB,ISYML)
                     ISYMAK = MULD2H(ISYVEC,ISYMBL)
C
                     KSCR2 = KEND1
                     KEND2 = KSCR2 + NT1AM(ISYMAI)
                     LWRK2 = LWORK - KEND2
C
                     IF (LWRK2 .LT. 0) THEN
                        CALL QUIT('Insufficient space in SO_KILJ')
                     END IF
C
                     IF (ISYMAI .GT. ISYMBJ) GOTO 170
C
                     DO 180 B = 1,NVIR(ISYMB)
C
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
                        NBL = IT1AM(ISYMB,ISYML)
     *                      + NVIR(ISYMB)*(L - 1) + B
C
                        CALL DZERO(WORK(KSCR2),NT1AM(ISYMAI))
C
                        DO 190 ISYMI = 1,NSYM
C
                           ISYMK = MULD2H(ISYMI,ISYMKI)
                           ISYMA = MULD2H(ISYMK,ISYMAK)
C
                           NVIRA = MAX(NVIR(ISYMA),1)
                           NRHFK = MAX(NRHF(ISYMK),1)
C
                           KOFF1 = IT2SQ(ISYMAK,ISYMBL)
     *                           + NT1AM(ISYMAK)*(NBL - 1)
     *                           + IT1AM(ISYMA,ISYMK) + 1
                           KOFF2 = KSCR1 + IMATIJ(ISYMK,ISYMI)
                           KOFF3 = KSCR2 + IT1AM(ISYMA,ISYMI)
C
                           CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *                                NRHF(ISYMK),ONE,T2AM(KOFF1),
     *                                NVIRA,WORK(KOFF2),NRHFK,ZERO,
     *                                WORK(KOFF3),NVIRA)
C                      T_(ij)^(ab)*(ki|lj)
C
  190                   CONTINUE
C
                        IF (ISYMAI .EQ. ISYMBJ) THEN
                           NTOT = NBJ
                        ELSE
                           NTOT = NT1AM(ISYMAI)
                        ENDIF
C
                        DO 200 NAI = 1,NTOT
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           ENDIF
C
                           OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                                   + WORK(KSCR2 + NAI - 1)
                           WRITE(LUPRI,'(A)') 'T_kl^cd (ki|lj)'
                           WRITE(LUPRI,*) OMEGA2(NAIBJ)
                           FLUSH(LUPRI)
C
  200                   CONTINUE
C
  180                CONTINUE
  170             CONTINUE
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
