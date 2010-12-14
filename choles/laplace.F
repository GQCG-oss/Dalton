C  /* Deck clap_vec */
      SUBROUTINE CLAP_VEC(FOCKD,,NPOINT,NUMIB,IB1,
     &                    ISYMB,ISYMA,ISYEMC,ISYMIJ,ISYIJC)
C
C    Alfredo Sanchez.    March 2005
C
C     Construct vector for Laplace decomposition.
C     The order of elements is:
C       1.- ijb for fixed #b from the batch over b (numib*nij)
C       2.- emc for fixed ISYEMC  (nckatr(isyemc))
C       3.- ijc for fixed ISYIJC  (nmaija(isyijc))
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRES = 1.0D-13)
C
      DIMENSION FOCKD(*),LAPVEC(*)
      DIMENSION ABSCI(NPOINT),WEI(NPOINT)
      DIMENSION IOF2(8)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      NIJ = NMATIJ(ISYMIJ)
C
C-------------------------------
C     Get abscissas and weights.
C-------------------------------
C
      CALL GAULEG(ZERO,ONE,ABSCI,WEI,NPOINT)
C
C-----------------------------------------------------
C     Calculate Laplace vectors.
C     The weight is included in virtual (second) part.
C-----------------------------------------------------
C
      IND = 0
      DO IP = 1,NPOINT
C
C     First part
C     ----------
C
      DO B = 1,NUMIB
         KOFFB = IVIR(ISYMB) + IB1 + B - 1
         DO ISYMJ = 1,NSYM
            ISYMI = MULD2H(ISYMIJ,ISYMJ)
            DO J = 1,NRHF(ISYMJ)
               KOFFJ = IRHF(ISYMJ) + J
               DO I = 1,NRHF(ISYMI)
                  KOFFI = IRHF(ISYMI) + I
                  OME = FOCKD(KOFFB) - FOCKD(KOFFI) - FOCKD(KOFFJ)
C
                  TMP = EXP(OME * ABSCI(IP)/(ABSCI(IP)-ONE)
                  CON = TMP/(ABSCI(IP)-ONE)**TWO
C
                  IND = IND + 1
                  LAPVEC(IND) = CON
C
               END DO
            END DO
         END DO
      END DO
C
C     Second part.
C     ------------
C
      DO ISYMC = 1,NSYM
         ISYMEM = MULD2H(ISYEMC,ISYMC)
         DO C = 1,NVIR(ISYMC)
            KOFFC = IVIR(ISYMC) + C
            DO ISYMM = 1,NSYM
               ISYME = MULD2H(ISYMEM,ISYMM)
               DO M = 1,NRHF(ISYMM)
                  KOFFM = IRHF(ISYMM) + M
                  DO E = 1,NVIR(ISYME)
                     KOFFE = IVIR(ISYME) + E
                     OME = FOCKD(KOFFC) + FOCKD(KOFFE) - FOCKD(KOFFM)
C
                     TMP = EXP(OME * ABSCI(IP)/(ABSCI(IP)-ONE)
                     CON = WEI(IP) * TMP/(ABSCI(IP)-ONE)**TWO
C
                     IND = IND + 1
                     LAPVEC(IND) = CON
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C     Third part.
C     -----------
C
      DO ISYMC = 1,NSYM
         JSYMIJ = MULD2H(ISYIJC,ISYMC)
         DO C = 1,NVIR(ISYMC)
            KOFFC = IVIR(ISYMC) + C
            DO ISYMJ = 1,NSYM
               ISYMI = MULD2H(JSYMIJ,ISYMJ)
               DO J = 1,NRHF(ISYMJ)
                  KOFFJ = IRHF(ISYMJ) + J
                  DO I = 1,NRHF(ISYMI)
                     KOFFI = IRHF(ISYMI) + I
                     OME = FOCKD(KOFFC) - FOCKD(KOFFI) - FOCKD(KOFFJ)
C
                     TMP = EXP(OME * ABSCI(IP)/(ABSCI(IP)-ONE)
                     CON = TMP/(ABSCI(IP)-ONE)**TWO
C
                     IND = IND + 1
                     LAPVEC(IND) = CON

C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      END DO
C
      RETURN
      END
C
