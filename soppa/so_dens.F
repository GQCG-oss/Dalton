#include <single.h>
C
C  /* Deck so_dens */
      SUBROUTINE SO_DENS(DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                   T2AM,LT2AM,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate second-order SOPPA one-electron density
C              matrices Dij and Dab from T2 amplitudes.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
C
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION T2AM(LT2AM)
      DIMENSION WORK(LWORK)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <soppinf.h>
C
C------------------------------
C     Statement function INDEX.
C------------------------------
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_DENS')
C
C---------------------------------
C     Initialize density matrices.
C---------------------------------
C
      CALL DZERO(DENSIJ,LDENSIJ)
      CALL DZERO(DENSAB,LDENSAB)
C
      DO 100 ISYMB = 1,NSYM
C
         DO 200 ISYMK = 1,NSYM
C
            ISYMBK = MULD2H(ISYMB,ISYMK)
            ISYMAJ = ISYMBK
C
            DO 300 ISYMJ = 1,NSYM
C
               ISYMA  = MULD2H(ISYMJ,ISYMAJ)
               ISYMAK = MULD2H(ISYMA,ISYMK)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMBK = MULD2H(ISYMB,ISYMK)
C
               DO 400 B = 1,NVIR(ISYMB)
C
C--------------------------------------------------------------
C                 Collect all MP2-amplitudes for the given B in
C                 PT2MP1 and PT2MP2 in different ordering and
C                 collect all T2-amplitudes for the given B in
C                 PT2AM1 and PT2AM2 in different ordering.
C--------------------------------------------------------------
C
                  LPT2MP1 = NVIR(ISYMA)*NRHF(ISYMK)*NRHF(ISYMJ)
                  LPT2MP2 = NRHF(ISYMJ)*NRHF(ISYMK)*NVIR(ISYMA)
C
                  LPT2AM1 = NRHF(ISYMJ)*NVIR(ISYMA)*NRHF(ISYMK)
                  LPT2AM2 = NVIR(ISYMA)*NRHF(ISYMJ)*NRHF(ISYMK)
C
                  KPT2MP1 = 1
                  KPT2MP2 = KPT2MP1 + LPT2MP1
                  KPT2AM1 = KPT2MP2 + LPT2MP2
                  KPT2AM2 = KPT2AM1 + LPT2AM1
                  KEND1   = KPT2AM2 + LPT2AM2
                  LWORK1  = LWORK   - KEND1
C
                  CALL SO_MEMMAX ('SO_DENS.1',LWORK1)
                  IF (LWORK1 .LT. 0) 
     &                CALL STOPIT('SO_DENS.1',' ',KEND1,LWORK)
C
                  DO 500 K = 1,NRHF(ISYMK)
C
                     DO 600 J = 1,NRHF(ISYMJ)
C
                        DO 700 A = 1,NVIR(ISYMA)
C
                           NAJ   = IT1AM(ISYMA,ISYMJ)
     &                           + NVIR(ISYMA)*(J-1) + A
                           NBK   = IT1AM(ISYMB,ISYMK)
     &                           + NVIR(ISYMB)*(K-1) + B
                           NAJBK = IT2AM(ISYMAJ,ISYMBK)
     &                           + INDEX(NAJ,NBK)
C
                           NAK   = IT1AM(ISYMA,ISYMK)
     &                           + NVIR(ISYMA)*(K-1) + A
                           NBJ   = IT1AM(ISYMB,ISYMJ)
     &                           + NVIR(ISYMB)*(J-1) + B
                           NAKBJ = IT2AM(ISYMAK,ISYMBJ)
     &                           + INDEX(NAK,NBJ)
C
                           NSQKJ = NRHF(ISYMK)*(J-1) + K
                           NSQJK = NRHF(ISYMJ)*(K-1) + J
                           NAKJ  = NVIR(ISYMA)*(NSQKJ-1) + A
                           NAJK  = NVIR(ISYMA)*(NSQJK-1) + A
C
                           NSQAK = NVIR(ISYMA)*(K-1) + A
                           NSQKA = NRHF(ISYMK)*(A-1) + K
                           NJAK  = NRHF(ISYMJ)*(NSQAK-1) + J
                           NJKA  = NRHF(ISYMJ)*(NSQKA-1) + J
C
                           T2MP  = TWO * T2AM(NAJBK) - T2AM(NAKBJ)
C
                           WORK(KPT2MP1 + NAKJ - 1) = T2MP
                           WORK(KPT2MP2 + NAJK - 1) = T2MP
C
                           WORK(KPT2AM1 + NJAK - 1) = T2AM(NAJBK)
                           WORK(KPT2AM2 + NJKA - 1) = T2AM(NAJBK)
C
  700                   CONTINUE
C
  600                CONTINUE
C
  500             CONTINUE
C
                  ISYMI = ISYMJ
                  ISYMC = ISYMA
C
C---------------------------------------------------------------------
C                 Calculate Dij by adding in contributions for each B.
C                 See eq. A.1 in J. Chem. Phys. 89, 5354 (1988)
C---------------------------------------------------------------------
C
                  NTOTJ  = MAX(NRHF(ISYMJ),1)
                  NTOTAK = MAX(NVIR(ISYMA)*NRHF(ISYMK),1)
C
                  KOFF1  = IIJDEN(ISYMJ,ISYMJ) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMJ),NRHF(ISYMJ),
     &                       NVIR(ISYMA)*NRHF(ISYMK),-ONE,WORK(KPT2AM1),
     &                       NTOTJ,WORK(KPT2MP1),NTOTAK,ONE,
     &                       DENSIJ(KOFF1),NTOTJ)
C
C---------------------------------------------------------------------
C                 Calculate Dab by adding in contributions for each B.
C                 See eq. A.2 in J. Chem. Phys. 89, 5354 (1988)
C---------------------------------------------------------------------
C
                  NTOTA  = MAX(NVIR(ISYMA),1)
                  NTOTJK = MAX(NRHF(ISYMJ)*NRHF(ISYMK),1)
C
                  KOFF2  = IABDEN(ISYMA,ISYMA) + 1
C
                  CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMA),
     &                       NRHF(ISYMJ)*NRHF(ISYMK),ONE,WORK(KPT2MP2),
     &                       NTOTA,WORK(KPT2AM2),NTOTJK,ONE,
     &                       DENSAB(KOFF2),NTOTA)
C
  400          CONTINUE
C
  300       CONTINUE
C
  200    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_DENS')
C
      RETURN
      END
