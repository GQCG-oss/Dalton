C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html 
C
C
C$Id: eri2cft.F,v 1.1.1.1 2001-02-08 13:33:13 hjj Exp $
C
#include <single.h>
C  /* Deck expcft */
      SUBROUTINE EXPCFT(ETIJ,NINNER,NPPPP,JMAX1,JMAX2,COOR1,COOR2,COORP,
     &                  IAB0,PINV,HPI,PA,PB,IELCTR,MAX1,MAX2,NCNT12,
     &                  IPRINT)
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0)
      INTEGER X
      DIMENSION ETIJ(NINNER,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3,*),
     &          COOR1(NPPPP,3), COOR2(NPPPP,3), COORP(NPPPP,3),
     &          IAB0(3), PINV(NPPPP),
     &          HPI(NINNER), PA(NINNER,3), PB(NINNER,3)
C
C     NSTRID = NPPPP/NINNER
      NSTRID = 1
      IF (IELCTR .EQ. 1) THEN
         SIGN = D1
      ELSE
         SIGN = - D1
      END IF
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('Output from EXPCFT',-1)
         WRITE (LUPRI,'(2X,A,I1,/)')
     &      'Hermite-to-Cartesian expansion coefficients for electron ',
     &       IELCTR
         WRITE (LUPRI, '(2X,A,2I5)')   'MAX1,MAX2          ',MAX1, MAX2
         WRITE (LUPRI, '(2X,A,3I5)')   'NPPPP,NINNER,NSTRID',
     &                                  NPPPP,NINNER,NSTRID
         WRITE (LUPRI, '(2X,A, I5)')   'NCNT12             ',NCNT12
         WRITE (LUPRI, '(2X,A,F12.6)') 'SIGN               ',SIGN
         WRITE (LUPRI, '(2X,A,3I5)')   'IAB0               ',
     &                                 (IAB0(I),I=1,3)
      END IF
C
C     1/(2P)
C     ======
C
      HPIFAC = SIGN*DP5
      IF (NSTRID .EQ. 1) THEN
         DO 100 I = 1, NINNER
            HPI(I) = HPIFAC*PINV(I)
  100    CONTINUE
      ELSE
         DO 110 I = 1, NINNER
            HPI(I) = HPIFAC*PINV(NSTRID*(I - 1) + 1)
  110    CONTINUE
      END IF
C
C     PA
C     ==
C
      IF (MAX1 .GT. 0) THEN
         IF (NSTRID .EQ. 1) THEN
            DO 200 I = 1, NINNER
               PA(I,1) = COORP(I,1) - COOR1(I,1)
               PA(I,2) = COORP(I,2) - COOR1(I,2)
               PA(I,3) = COORP(I,3) - COOR1(I,3)
  200       CONTINUE
         ELSE
            DO 210 I = 1, NINNER
               IJ = NSTRID*(I - 1) + 1
               PA(I,1) = COORP(IJ,1) - COOR1(IJ,1)
               PA(I,2) = COORP(IJ,2) - COOR1(IJ,2)
               PA(I,3) = COORP(IJ,3) - COOR1(IJ,3)
  210       CONTINUE
         END IF
      END IF
C
C     PB
C     ==
C
      IF (MAX2 .GT. 0) THEN
         IF (NSTRID .EQ. 1) THEN
            DO 300 I = 1, NINNER
               PB(I,1) = COORP(I,1) - COOR2(I,1)
               PB(I,2) = COORP(I,2) - COOR2(I,2)
               PB(I,3) = COORP(I,3) - COOR2(I,3)
  300       CONTINUE
         ELSE
            DO 310 I = 1, NINNER
               IJ = NSTRID*(I - 1) + 1
               PB(I,1) = COORP(IJ,1) - COOR2(IJ,1)
               PB(I,2) = COORP(IJ,2) - COOR2(IJ,2)
               PB(I,3) = COORP(IJ,3) - COOR2(IJ,3)
  310       CONTINUE
         END IF
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('HPI in EXPCFT',1)
         CALL OUTPUT(HPI,1,1,1,NINNER,1,NINNER,1,LUPRI)
         IF (MAX1 .GT. 0) THEN
            CALL HEADER('PA in EXPCFT',1)
            CALL OUTPUT(PA,1,NINNER,1,3,NINNER,3,1,LUPRI)
         END IF
         IF (MAX2 .GT. 0) THEN
            CALL HEADER('PB in EXPCFT',1)
            CALL OUTPUT(PB,1,NINNER,1,3,NINNER,3,1,LUPRI)
         END IF
      END IF
C
C     Expansion coefficients
C     ======================
C
      DO 400 X = 1, 3
         CALL ERIECF(ETIJ(1,0,0,0,X,1),NINNER,MAX1,MAX2,JMAX1,JMAX2,
     &               SIGN,IAB0(X),PA(1,X),PB(1,X),HPI,X,IPRINT)
  400 CONTINUE
C
      RETURN
      END
C  /* Deck eriecf */
      SUBROUTINE ERIECF(ETIJ,NINNER,MAXI,MAXJ,JMAXI,JMAXJ,SIGN,
     &                  IAB0,PA,PB,HPI,X,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.00 D00, D2 = 2.00 D00)
      INTEGER T, X
      CHARACTER WORD*4
      DIMENSION PA(NINNER), PB(NINNER), HPI(NINNER)
      DIMENSION ETIJ(NINNER,0:JMAXI+JMAXJ,0:JMAXI,0:JMAXJ)
#include <ibtfun.h>
#include <sdpre.h>
C
C     ****************************
C     ********** AB > 0 **********
C     ****************************
C
      IF (IAB0 .EQ. 0) THEN
C
C        Run over I (J = 0)
C        ==================
C
         DO 100 I = 0, MAXI
C
C           E(0,0)
C
            IF (I .EQ. 0) THEN
               DO 200 K = 1, NINNER
                  ETIJ(K,0,0,0) = D1
  200          CONTINUE
C
C           E(1,0)
C
            ELSE IF (I .EQ. 1) THEN
               DO 300 K = 1, NINNER
                  ETIJ(K,0,1,0) = PA(K)
                  ETIJ(K,1,1,0) = HPI(K)
  300          CONTINUE
C
C           E(2,0)
C
            ELSE IF (I .EQ. 2) THEN
               DO 400 K = 1, NINNER
                  PAK = PA(K)
                  HPK = HPI(K)
                  ETIJ(K,0,2,0) = PAK*PAK + SIGN*HPK
                  ETIJ(K,1,2,0) = D2*PAK*HPK
                  ETIJ(K,2,2,0) = HPK*HPK
  400          CONTINUE
C
C            E(I,0)
C
            ELSE
               DO 500 K = 1, NINNER
                  PAK = PA(K)
                  HPK = HPI(K)
                  ETIJ(K,  0,I,0) = PAK*ETIJ(K,  0,I-1,0)
     &                           + SIGN*ETIJ(K,  1,I-1,0)
                  ETIJ(K,I-1,I,0) = HPK*ETIJ(K,I-2,I-1,0)
     &                            + PAK*ETIJ(K,I-1,I-1,0)
                  ETIJ(K,  I,I,0) = HPK*ETIJ(K,I-1,I-1,0)

  500          CONTINUE
               DO 510 T = 1, I - 2
                  T1 = SIGN*SDPRE(T + 1)
                  DO 520 K = 1, NINNER
                     ETIJ(K,T,I,0) = HPI(K)*ETIJ(K,T-1,I-1,0)
     &                              + PA(K)*ETIJ(K,  T,I-1,0)
     &                                 + T1*ETIJ(K,T+1,I-1,0)
  520             CONTINUE
  510          CONTINUE
            END IF
C
C           Run over J
C           ==========
C
            DO 600 J = 1, MAXJ
               IJ = I + J
C
C              E(0,1)
C
               IF (IJ .EQ. 1) THEN
                  DO 700 K = 1, NINNER
                     ETIJ(K,0,0,1) = PB(K)
                     ETIJ(K,1,0,1) = HPI(K)
  700             CONTINUE
               ELSE IF (IJ .EQ. 2) THEN
C
C                 E(0,2)
C
                  IF (I .EQ. 0) THEN
                     DO 800 K = 1, NINNER
                        PBK = PB(K)
                        HPK = HPI(K)
                        ETIJ(K,0,0,2) = PBK*PBK + SIGN*HPK
                        ETIJ(K,1,0,2) = D2*PBK*HPK
                        ETIJ(K,2,0,2) = HPK*HPK
  800                CONTINUE
C
C                 E(1,1)
C
                  ELSE
                     DO 810 K = 1, NINNER
                        PAK = PA(K)
                        PBK = PB(K)
                        HPK = HPI(K)
                        ETIJ(K,0,1,1) = PAK*PBK + SIGN*HPK
                        ETIJ(K,1,1,1) = (PAK + PBK)*HPK
                        ETIJ(K,2,1,1) = HPK*HPK
  810                CONTINUE
                  END IF
C
C              E(I,J)
C
               ELSE
                  DO 900 K = 1, NINNER
                     PBK = PB(K)
                     HPK = HPI(K)
                     ETIJ(K,   0,I,J) = PBK*ETIJ(K,   0,I,J-1)
     &                               + SIGN*ETIJ(K,   1,I,J-1)
                     ETIJ(K,IJ-1,I,J) = HPK*ETIJ(K,IJ-2,I,J-1)
     &                                + PBK*ETIJ(K,IJ-1,I,J-1)
                     ETIJ(K,  IJ,I,J) = HPK*ETIJ(K,IJ-1,I,J-1)
  900             CONTINUE
                  DO 910 T = 1, IJ - 2
                     T1 = SIGN*SDPRE(T + 1)
                     DO 920 K = 1, NINNER
                        ETIJ(K,T,I,J)=HPI(K)*ETIJ(K,T-1,I,J-1)
     &                               + PB(K)*ETIJ(K,  T,I,J-1)
     &                                  + T1*ETIJ(K,T+1,I,J-1)
  920               CONTINUE
  910             CONTINUE
               END IF
  600       CONTINUE
  100    CONTINUE
C
C     ****************************
C     ********** AB = 0 **********
C     ****************************
C
      ELSE
C
C        Run over I (J = 0)
C        ==================
C
         DO 105 I = 0, MAXI
C
C           E(0,0)
C
            IF (I .EQ. 0) THEN
               DO 205 K = 1, NINNER
                  ETIJ(K,0,0,0) = D1
  205          CONTINUE
C
C           E(1,0)
C
            ELSE IF (I .EQ. 1) THEN
               DO 305 K = 1, NINNER
                  ETIJ(K,1,1,0) = HPI(K)
  305          CONTINUE
C
C           E(2,0)
C
            ELSE IF (I .EQ. 2) THEN
               DO 405 K = 1, NINNER
                  HPK = HPI(K)
                  ETIJ(K,0,2,0) = SIGN*HPK
                  ETIJ(K,2,2,0) = HPK*HPK
  405          CONTINUE
C
C            E(I,0)
C
            ELSE
               IF (IBTAND(1,I) .EQ. 0) THEN
                  DO 505 K = 1, NINNER
                     ETIJ(K,0,I,0) =   SIGN*ETIJ(K,  1,I-1,0)
                     ETIJ(K,I,I,0) = HPI(K)*ETIJ(K,I-1,I-1,0)
  505             CONTINUE
               ELSE
                  DO 506 K = 1, NINNER
                     ETIJ(K,I,I,0) = HPI(K)*ETIJ(K,I-1,I-1,0)
  506             CONTINUE
               END IF
               DO 515 T = 2 - IBTAND(1,I), I - 2, 2
                  T1 = SIGN*SDPRE(T + 1)
                  DO 525 K = 1, NINNER
                     ETIJ(K,T,I,0) = HPI(K)*ETIJ(K,T-1,I-1,0)
     &                                 + T1*ETIJ(K,T+1,I-1,0)
  525             CONTINUE
  515          CONTINUE
            END IF
C
C           Run over J
C           ==========
C
            DO 605 J = 1, MAXJ
               IJ = I + J
C
C              E(0,1)
C
               IF (IJ .EQ. 1) THEN
                  DO 705 K = 1, NINNER
                     ETIJ(K,1,0,1) = HPI(K)
  705             CONTINUE
C
C              E(1,1) AND E(0,2)
C
               ELSE IF (IJ .EQ. 2) THEN
                  DO 805 K = 1, NINNER
                     HPK = HPI(K)
                     ETIJ(K,0,I,J) = SIGN*HPK
                     ETIJ(K,2,I,J) = HPK*HPK
  805             CONTINUE
C
C              E(I,J)
C
               ELSE
                  IF (IBTAND(1,IJ) .EQ. 0) THEN
                     DO 905 K = 1, NINNER
                       ETIJ(K, 0,I,J) =   SIGN*ETIJ(K,   1,I,J-1)
                       ETIJ(K,IJ,I,J) = HPI(K)*ETIJ(K,IJ-1,I,J-1)
  905                CONTINUE
                  ELSE
                     DO 906 K = 1, NINNER
                       ETIJ(K,IJ,I,J) = HPI(K)*ETIJ(K,IJ-1,I,J-1)
  906                CONTINUE
                  END IF
                  DO 915 T = 2 - IBTAND(1,IJ), IJ - 2, 2
                     T1 = SIGN*SDPRE(T + 1)
                     DO 925 K = 1, NINNER
                        ETIJ(K,T,I,J) = HPI(K)*ETIJ(K,T-1,I,J-1)
     &                                    + T1*ETIJ(K,T+1,I,J-1)
  925                CONTINUE
  915             CONTINUE
               END IF
  605       CONTINUE
  105    CONTINUE
      END IF
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERIECF',-1)
         WRITE (LUPRI,'(2X,A,2I5)') 'JMAXI,JMAXJ ', JMAXI, JMAXJ
         WRITE (LUPRI,'(2X,A,2I5)') 'MAXI,MAXJ   ', MAXI, MAXJ
         WRITE (LUPRI,'(2X,A, I5)') 'IAB0        ', IAB0
         IF (X .EQ. 1) WORD = 'EX00'
         IF (X .EQ. 2) WORD = 'EY00'
         IF (X .EQ. 3) WORD = 'EZ00'
         IF (IPRINT .GT. 20) THEN
            DO 1000 I = 0, MAXI
            DO 1000 J = 0, MAXJ
            DO 1000 T = 0, I + J
            IF (IBTAND(I + J - T,IAB0) .EQ. 0) THEN
               WRITE (LUPRI,'(/,1X,A4,A1,I1,A1,I1,A1,I1,A1,/)')
     &              WORD, '(', I, ',', J, ', ',  T, ')'
               WRITE (LUPRI,'(1X,6F12.8)') (ETIJ(K,T,I,J),K=1,NINNER)
            END IF
 1000       CONTINUE
         END IF
      END IF
      RETURN
      END
