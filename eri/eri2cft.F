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
C$Id: eri2cft.F,v 1.3 2001-04-03 14:25:11 vebjornb Exp $
C
#include <single.h>
C  /* Deck expcft */
      SUBROUTINE EXPCFT(ETIJ,NPP12,JMAX1,JMAX2,COOR12,EXP12,
     &                  IAB0,HPI,PA,PB,IRUTIN,IELCTR,MAX1,MAX2,
     &                  NCNT12,IDL,ISCAL1,ISCAL2,IPRINT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0)
      INTEGER X
      DIMENSION ETIJ(NPP12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3,*),
     &          COOR12(NPP12,3), EXP12(NPP12,3), IAB0(3), 
     &          HPI(NPP12), PA(NPP12,3), PB(NPP12,3)
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('Output from EXPCFT',-1)
         WRITE (LUPRI,'(2X,A,I1,/)')
     &      'Hermite-to-Cartesian expansion coefficients for electron ',
     &       IELCTR
         WRITE (LUPRI, '(2X,A, I5)')   'IRUTIN    ',IRUTIN
         WRITE (LUPRI, '(2X,A, I5)')   'IELCTR    ',IELCTR
         WRITE (LUPRI, '(2X,A,2I5)')   'MAX1,MAX2 ',MAX1, MAX2
         WRITE (LUPRI, '(2X,A, I5)')   'NPP12     ',NPP12
         WRITE (LUPRI, '(2X,A, I5)')   'NCNT12    ',NCNT12
         WRITE (LUPRI, '(2X,A,3I5)')   'IAB0      ',(IAB0(I),I=1,3)
      END IF
C
C     Sign factor
C     ===========
C
      SIGN = D1
      IF (IELCTR.EQ.2) SIGN = -D1
C
C     Extract exponents and coordinates
C     =================================
C
      CALL EXPVEC(HPI,PA,PB,COOR12,EXP12,SIGN,NPP12,MAX1,MAX2,IDL,
     &            IPRINT)
C
C     Expansion coefficients
C     ======================
C
      DO 100 X = 1, 3
         CALL ERIECF(ETIJ(1,0,0,0,X,1),ETIJ(1,0,0,0,X,2),
     &               ETIJ(1,0,0,0,X,2),NPP12,MAX1,MAX2,JMAX1,JMAX2,
     &               SIGN,IAB0(X),PA(1,X),PB(1,X),EXP12,HPI,
     &               X,0,IPRINT)
  100 CONTINUE
C
C     Scale coefficients
C     ==================
C
      IF (ISCAL1 + ISCAL2 .GT. 0) THEN
         CALL ERISCL(ETIJ,EXP12,HPI,ISCAL1,ISCAL2,NPP12,JMAX1,JMAX2,
     &               MAX1,MAX2,IAB0,IPRINT)
      END IF
C
C     Derivative coefficients
C     =======================
C
      IF (IDL .GT. 0) THEN
         DO 200 X = 1, 3
            CALL ERIECF(ETIJ(1,0,0,0,X,2),ETIJ(1,0,0,0,X,1),
     &                  ETIJ(1,0,0,0,X,1),NPP12,MAX1,MAX2,JMAX1,JMAX2,
     &                  SIGN,IAB0(X),PA(1,X),PB(1,X),EXP12,HPI,
     &                  X,1,IPRINT)
  200    CONTINUE
      ENDIF
C
      IF (IDL .GT. 1) STOP 'No higher than first derivative yet'
C
      RETURN
      END
C  /* Deck EXPVEC */
      SUBROUTINE EXPVEC(HPI,PA,PB,COOR12,EXP12,SIGN,
     &                  NPP12,MAX1,MAX2,IDL,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP5 = 0.5D0)
      DIMENSION HPI(NPP12), 
     &          PA(NPP12,3), PB(NPP12,3),
     &          COOR12(NPP12,3), EXP12(NPP12,3)
C
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI,'(2X,A,I5)') 'NPP12',NPP12
         CALL HEADER('Output from EXPVEC',-1)
         WRITE (LUPRI,'(2X,A,2I5)') 'MAX1,MAX2        ',MAX1,MAX2 
         WRITE (LUPRI,'(2X,A,F12.6)') ' SIGN ', SIGN
      END IF
C
C     1/(2P)
C     ======
C
      HPIFAC = SIGN*DP5
      DO 100 I = 1, NPP12
         HPI(I) = HPIFAC*EXP12(I,3)
 100  CONTINUE
C
C     PA and PB
C     =========
C
      IF (MAX1 .GT. 0 .AND. MAX2 .GT. 0) THEN
         DO 200 I = 1, NPP12
            EXP1    =  EXP12(I,1)
            EXP2    = -EXP12(I,2)
            PA(I,1) = EXP2*COOR12(I,1)
            PA(I,2) = EXP2*COOR12(I,2)
            PA(I,3) = EXP2*COOR12(I,3)
            PB(I,1) = EXP1*COOR12(I,1)
            PB(I,2) = EXP1*COOR12(I,2)
            PB(I,3) = EXP1*COOR12(I,3)
 200     CONTINUE
      ELSE
         IF (MAX1 .GT. 0 .OR. IDL.GT.0) THEN
            DO 300 I = 1, NPP12
               EXP2    = -EXP12(I,2)
               PA(I,1) = EXP2*COOR12(I,1)
               PA(I,2) = EXP2*COOR12(I,2)
               PA(I,3) = EXP2*COOR12(I,3)
 300        CONTINUE
         END IF
         IF (MAX2 .GT. 0) THEN
            DO 400 I = 1, NPP12
               EXP1    = EXP12(I,1)
               PB(I,1) = EXP1*COOR12(I,1)
               PB(I,2) = EXP1*COOR12(I,2)
               PB(I,3) = EXP1*COOR12(I,3)
 400        CONTINUE
         END IF
      END IF
C
C     Print
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('HPI in EXPVEC',1)
         CALL OUTPUT(HPI,1,1,1,NPP12,1,NPP12,1,LUPRI)
         IF (MAX1 .GT. 0) THEN
            CALL HEADER('PA in EXPVEC',1)
            CALL OUTPUT(PA,1,NPP12,1,3,NPP12,3,1,LUPRI)
         END IF
         IF (MAX2 .GT. 0) THEN
            CALL HEADER('PB in EXPVEC',1)
            CALL OUTPUT(PB,1,NPP12,1,3,NPP12,3,1,LUPRI)
         END IF
      END IF
C
      RETURN
      END
C  /* Deck eriecf */
      SUBROUTINE ERIECF(ETIJ,ETIJM1,ETIJM2,NPP12,MAXI,MAXJ,JMAXI,JMAXJ,
     *                  SIGN,IAB0,PA,PB,EXP12,HPI,X,IDL,IPRINT)
C
C     A. Halkier & T. Helgaker Jan. 1999: Generalised Version 
C     that handles derivative coefficients, where IDL gives the order
C     of differentiation.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.00 D00, D2 = 2.00 D00)
      INTEGER T, X
      LOGICAL UNDIFF
      CHARACTER WORD*4
      DIMENSION PA(NPP12), PB(NPP12), EXP12(NPP12,3), HPI(NPP12)
      DIMENSION ETIJ(NPP12,0:JMAXI+JMAXJ,0:JMAXI,0:JMAXJ)
      DIMENSION ETIJM1(NPP12,0:JMAXI+JMAXJ,0:JMAXI,0:JMAXJ)
      DIMENSION ETIJM2(NPP12,0:JMAXI+JMAXJ,0:JMAXI,0:JMAXJ)
#include <ibtfun.h>
#include <sdpre.h>
C
      DLEVEL = SDPRE(IDL)
      UNDIFF = IDL .EQ. 0
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
               IF (IDL .EQ. 0) THEN
                  DO 200 K = 1, NPP12
                     ETIJ(K,0,0,0) = D1
  200             CONTINUE
               ELSE IF (IDL .EQ. 1) THEN
                  DO 250 K = 1,NPP12
                     PAK  = PA(K)/EXP12(K,3)
                     EXP1 = D2*EXP12(K,1)
                     ETIJ(K,0,0,0) = PAK*EXP1
  250             CONTINUE
               ELSE IF (IDL .GT. 1) THEN
                  DIDLM = D1 - DLEVEL
                  DO 260 K = 1,NPP12
                     EXPH = EXP12(K,1)/EXP12(K,3)
                     PAK1 = PA(K)*ETIJM1(K,0,0,0)
                     PAK2 = DIDLM*EXP12(K,2)*ETIJM2(K,0,0,0)
                     ETIJ(K,0,0,0) = D2*EXPH*(PAK1 + PAK2)
  260             CONTINUE
               ENDIF
C
C           E(1,0)
C
            ELSE IF (I .EQ. 1) THEN
               IF (IDL. EQ. 0) THEN
                  DO 300 K = 1, NPP12
                     ETIJ(K,0,1,0) = PA(K)
                     ETIJ(K,1,1,0) = HPI(K)
  300             CONTINUE
               ELSE 
                  DO 350 K = 1,NPP12
                     PAK  = PA(K)
                     HPK  = HPI(K)
                     EXP2 = -DLEVEL*EXP12(K,2)
                     ETIJ(K,1,1,0) =    HPK*ETIJ(K,0,0,0)
                     ETIJ(K,0,1,0) =    PAK*ETIJ(K,0,0,0)
     *                             + EXP2*ETIJM1(K,0,0,0)
  350             CONTINUE
               ENDIF
C
C           E(2,0)
C
            ELSE IF (I .EQ. 2) THEN
               IF (IDL .EQ. 0) THEN
                  DO 400 K = 1, NPP12
                     PAK = PA(K)
                     HPK = HPI(K)
                     ETIJ(K,0,2,0) = PAK*PAK + SIGN*HPK
                     ETIJ(K,1,2,0) = D2*PAK*HPK
                     ETIJ(K,2,2,0) = HPK*HPK
  400             CONTINUE
               ELSE
                  DO 450 K = 1,NPP12
                     PAK  = PA(K)
                     HPK  = HPI(K)
                     EXP2 = -DLEVEL*EXP12(K,2)
                     ETIJ(K,0,2,0) =   SIGN*ETIJ(K,1,1,0)
     *                             +    PAK*ETIJ(K,0,1,0)
     *                             + EXP2*ETIJM1(K,0,1,0)
                     ETIJ(K,1,2,0) =    HPK*ETIJ(K,0,1,0)
     *                             +    PAK*ETIJ(K,1,1,0)
     *                             + EXP2*ETIJM1(K,1,1,0)
                     ETIJ(K,2,2,0) =    HPK*ETIJ(K,1,1,0)
  450             CONTINUE
               ENDIF
C
C            E(I,0)
C
            ELSE
               DO 500 K = 1, NPP12
                  PAK = PA(K)
                  HPK = HPI(K)
                  ETIJ(K,  0,I,0) = PAK*ETIJ(K,  0,I-1,0)
     &                           + SIGN*ETIJ(K,  1,I-1,0)
                  ETIJ(K,I-1,I,0) = HPK*ETIJ(K,I-2,I-1,0)
     &                            + PAK*ETIJ(K,I-1,I-1,0)
                  ETIJ(K,  I,I,0) = HPK*ETIJ(K,I-1,I-1,0)
  500          CONTINUE
               IF (IDL .GT. 0) THEN
                  DO 550 K = 1,NPP12
                     EXP2 = -DLEVEL*EXP12(K,2)
                     ETIJ(K,  0,I,0) =        ETIJ(K,  0,  I,0)
     *                               + EXP2*ETIJM1(K,  0,I-1,0)
                     ETIJ(K,I-1,I,0) =        ETIJ(K,I-1,  I,0)
     *                               + EXP2*ETIJM1(K,I-1,I-1,0)
  550             CONTINUE
               END IF
               DO 510 T = 1, I - 2
                  T1 = SIGN*SDPRE(T + 1)
                  DO 520 K = 1, NPP12
                     ETIJ(K,T,I,0) = HPI(K)*ETIJ(K,T-1,I-1,0)
     &                              + PA(K)*ETIJ(K,  T,I-1,0)
     &                                 + T1*ETIJ(K,T+1,I-1,0)
  520             CONTINUE
  510          CONTINUE
               IF (IDL .GT. 0) THEN
                  DO 560 T = 1,I-2
                     DO 570 K = 1,NPP12
                        EXP2 = -DLEVEL*EXP12(K,2)
                        ETIJ(K,T,I,0) =        ETIJ(K,T,  I,0)
     *                                + EXP2*ETIJM1(K,T,I-1,0)
  570                CONTINUE
  560             CONTINUE
               ENDIF
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
                  IF (IDL .EQ. 0) THEN
                     DO 700 K = 1, NPP12
                        ETIJ(K,0,0,1) = PB(K)
                        ETIJ(K,1,0,1) = HPI(K)
  700                CONTINUE
                  ELSE
                     DO 750 K = 1, NPP12
                        HPK  = HPI(K)
                        PBK  = PB(K)
                        EXP1 = DLEVEL*EXP12(K,1)
                        ETIJ(K,1,0,1) =    HPK*ETIJ(K,0,0,0)
                        ETIJ(K,0,0,1) =    PBK*ETIJ(K,0,0,0)
     *                                + EXP1*ETIJM1(K,0,0,0)
  750                CONTINUE
                  ENDIF
               ELSE IF (IJ .EQ. 2) THEN
C
C                 E(0,2)
C
                  IF (IDL .EQ. 0) THEN
                     IF (I .EQ. 0) THEN
                        DO 800 K = 1, NPP12
                           PBK = PB(K)
                           HPK = HPI(K)
                           ETIJ(K,0,0,2) = PBK*PBK + SIGN*HPK
                           ETIJ(K,1,0,2) = D2*PBK*HPK
                           ETIJ(K,2,0,2) = HPK*HPK
  800                   CONTINUE
C
C                    E(1,1)
C
                     ELSE
                        DO 810 K = 1, NPP12
                           PAK = PA(K)
                           PBK = PB(K)
                           HPK = HPI(K)
                           ETIJ(K,0,1,1) = PAK*PBK + SIGN*HPK
                           ETIJ(K,1,1,1) = (PAK + PBK)*HPK
                           ETIJ(K,2,1,1) = HPK*HPK
  810                   CONTINUE
                     END IF
                  ELSE
                     IF (I .EQ. 0) THEN
                        DO 850 K = 1, NPP12
                           PBK  = PB(K)
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,0,0,2) =   SIGN*ETIJ(K,1,0,1)
     *                                   +    PBK*ETIJ(K,0,0,1)
     *                                   + EXP1*ETIJM1(K,0,0,1)
                           ETIJ(K,1,0,2) =    HPK*ETIJ(K,0,0,1)
     *                                   +    PBK*ETIJ(K,1,0,1)
     *                                   + EXP1*ETIJM1(K,1,0,1)
                           ETIJ(K,2,0,2) =    HPK*ETIJ(K,1,0,1)
  850                   CONTINUE
                     ELSE
                        DO 860 K = 1, NPP12
                           PBK  = PB(K)
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,0,1,1) =   SIGN*ETIJ(K,1,1,0)
     *                                   +    PBK*ETIJ(K,0,1,0)
     *                                   + EXP1*ETIJM1(K,0,1,0)
                           ETIJ(K,1,1,1) =    HPK*ETIJ(K,0,1,0)
     *                                   +    PBK*ETIJ(K,1,1,0)
     *                                   + EXP1*ETIJM1(K,1,1,0)
                           ETIJ(K,2,1,1) = HPK*ETIJ(K,1,1,0)
  860                   CONTINUE
                     END IF
                  ENDIF
C
C              E(I,J)
C
               ELSE
                  DO 900 K = 1, NPP12
                     PBK = PB(K)
                     HPK = HPI(K)
                     ETIJ(K,   0,I,J) = PBK*ETIJ(K,   0,I,J-1)
     &                               + SIGN*ETIJ(K,   1,I,J-1)
                     ETIJ(K,IJ-1,I,J) = HPK*ETIJ(K,IJ-2,I,J-1)
     &                                + PBK*ETIJ(K,IJ-1,I,J-1)
                     ETIJ(K,  IJ,I,J) = HPK*ETIJ(K,IJ-1,I,J-1)
  900             CONTINUE
                  IF (IDL .GT. 0) THEN
                     DO 950 K = 1, NPP12
                        EXP1 = DLEVEL*EXP12(K,1)
                        ETIJ(K,   0,I,J) =        ETIJ(K,   0,I,  J)
     *                                   + EXP1*ETIJM1(K,   0,I,J-1)
                        ETIJ(K,IJ-1,I,J) =        ETIJ(K,IJ-1,I,  J)
     *                                   + EXP1*ETIJM1(K,IJ-1,I,J-1)
  950                CONTINUE
                  END IF
                  DO 910 T = 1, IJ - 2
                     T1 = SIGN*SDPRE(T + 1)
                     DO 920 K = 1, NPP12
                        ETIJ(K,T,I,J)=HPI(K)*ETIJ(K,T-1,I,J-1)
     &                               + PB(K)*ETIJ(K,  T,I,J-1)
     &                                  + T1*ETIJ(K,T+1,I,J-1)
  920                CONTINUE
  910             CONTINUE
                  IF (IDL .GT. 0) THEN
                     DO 960 T = 1, IJ - 2
                        DO 970 K = 1, NPP12
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,T,I,J) =        ETIJ(K,T,I,  J)
     *                                   + EXP1*ETIJM1(K,T,I,J-1)
  970                   CONTINUE
  960                CONTINUE
                  ENDIF
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
               IF (IDL .EQ. 0) THEN
                  DO 205 K = 1, NPP12
                     ETIJ(K,0,0,0) = D1
  205             CONTINUE
               ELSE
                  IF (IBTAND(1,IDL) .EQ. 0) THEN
                     DIDLM = D1 - DLEVEL
                     DO 255 K = 1, NPP12
                        EXPH = EXP12(K,1)/EXP12(K,3)
                        PAKH = DIDLM*EXP12(K,2)*ETIJM2(K,0,0,0)
                        ETIJ(K,0,0,0) = D2*EXPH*PAKH
  255                CONTINUE
                  ENDIF
               ENDIF
C
C           E(1,0)
C
            ELSE IF (I .EQ. 1) THEN
               IF (IDL .EQ. 0) THEN
                  DO 305 K = 1, NPP12
                     ETIJ(K,1,1,0) = HPI(K)
  305             CONTINUE
               ELSE IF (IDL .EQ. 1) THEN
                  DO 355 K = 1, NPP12
                     ETIJ(K,0,1,0) = -DLEVEL*EXP12(K,2)
  355             CONTINUE
               ELSE IF ((IDL .GT. 1) .AND. (IBTAND(1,IDL) .EQ. 0)) THEN
                  DO 365 K = 1, NPP12
                     ETIJ(K,1,1,0) = HPI(K)*ETIJ(K,0,0,0)
  365             CONTINUE
               ELSE
                  DO 375 K = 1, NPP12
                     ETIJ(K,0,1,0) = -DLEVEL*EXP12(K,2)*ETIJM1(K,0,0,0)
  375             CONTINUE
               ENDIF
C
C           E(2,0)
C
            ELSE IF (I .EQ. 2) THEN
               IF (IDL .EQ. 0) THEN
                  DO 405 K = 1, NPP12
                     HPK = HPI(K)
                     ETIJ(K,0,2,0) = SIGN*HPK
                     ETIJ(K,2,2,0) = HPK*HPK
  405             CONTINUE
               ELSE IF (IDL .EQ. 1) THEN
                  DO 455 K = 1, NPP12
                     HPK  = HPI(K)
                     EXP2 = -D2*EXP12(K,2)
                     ETIJ(K,1,2,0) = HPK*EXP2
  455             CONTINUE
               ELSE IF ((IDL .GT. 1) .AND. (IBTAND(1,IDL) .EQ. 0)) THEN
                  DO 465 K = 1, NPP12
                     EXP2 = -DLEVEL*EXP12(K,2)
                     ETIJ(K,0,2,0) = EXP2*ETIJM1(K,0,1,0)
     *                             +   SIGN*ETIJ(K,1,1,0)
                     ETIJ(K,2,2,0) = HPI(K)*ETIJ(K,1,1,0)
  465             CONTINUE
               ELSE 
                  DO 475 K = 1, NPP12
                     HPK  = HPI(K)
                     EXP2 = -DLEVEL*EXP12(K,2)
                     ETIJ(K,1,2,0) =    HPK*ETIJ(K,0,1,0)
     *                             + EXP2*ETIJM1(K,1,1,0)
  475             CONTINUE
               ENDIF
C
C            E(I,0)
C
            ELSE
               IF ((IDL .EQ. 0) .OR.
     *             ((IDL .GT. 1) .AND. (IBTAND(1,IDL) .EQ. 0))) THEN
                  IF (IBTAND(1,I) .EQ. 0) THEN
                     DO 505 K = 1, NPP12
                        ETIJ(K,0,I,0) =   SIGN*ETIJ(K,  1,I-1,0)
                        ETIJ(K,I,I,0) = HPI(K)*ETIJ(K,I-1,I-1,0)
  505                CONTINUE
                     IF (IDL .GT. 1) THEN
                        DO 504 K = 1, NPP12
                           EXP2 = -DLEVEL*EXP12(K,2)
                           ETIJ(K,0,I,0) =        ETIJ(K,0,  I,0)
     *                                   + EXP2*ETIJM1(K,0,I-1,0)
  504                   CONTINUE
                     ENDIF
                  ELSE
                     DO 506 K = 1, NPP12
                        ETIJ(K,I,I,0) = HPI(K)*ETIJ(K,I-1,I-1,0)
  506                CONTINUE
                  END IF
                  DO 515 T = 2 - IBTAND(1,I), I - 2, 2
                     T1 = SIGN*SDPRE(T + 1)
                     DO 525 K = 1, NPP12
                        ETIJ(K,T,I,0) = HPI(K)*ETIJ(K,T-1,I-1,0)
     &                                    + T1*ETIJ(K,T+1,I-1,0)
  525                CONTINUE
                     IF (IDL .GT. 1) THEN
                        DO 524 K = 1, NPP12
                           EXP2 = -DLEVEL*EXP12(K,2)
                           ETIJ(K,T,I,0) =        ETIJ(K,T,I,0)
     *                                   + EXP2*ETIJM1(K,T,I,0)
  524                   CONTINUE
                     ENDIF
  515             CONTINUE
               ELSE 
                  IF (IBTAND(1,I) .EQ. 0) THEN
                     DO 555 K = 1,NPP12
                        HPK  = HPI(K)
                        EXP2 = -DLEVEL*EXP12(K,2)
                        ETIJ(K,I-1,I,0) =    HPK*ETIJ(K,I-2,I-1,0)
     *                                  + EXP2*ETIJM1(K,I-1,I-1,0)
  555                CONTINUE
                  ELSE
                     DO 556 K = 1,NPP12
                        HPK  = HPI(K)
                        EXP2 = -DLEVEL*EXP12(K,2)
                        ETIJ(K,I-1,I,0) =    HPK*ETIJ(K,I-2,I-1,0)
     *                                  + EXP2*ETIJM1(K,I-1,I-1,0)
                        ETIJ(K,  0,I,0) = EXP2*ETIJM1(K,  0,I-1,0)
     *                                  +   SIGN*ETIJ(K,  1,I-1,0)
  556                CONTINUE
                  ENDIF
                  DO 565 T = 1 + IBTAND(1,I), I - 3, 2
                     T1   = SIGN*SDPRE(T + 1)
                     DO 575 K = 1, NPP12
                        HPK  = HPI(K)
                        EXP2 = -DLEVEL*EXP12(K,2)
                        ETIJ(K,T,I,0) =    HPK*ETIJ(K,T-1,I-1,0)
     *                                + EXP2*ETIJM1(K,  T,I-1,0)
     *                                +     T1*ETIJ(K,T+1,I-1,0)
  575                CONTINUE
  565             CONTINUE
               ENDIF
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
                  IF (IDL .EQ. 0) THEN
                     DO 705 K = 1, NPP12
                        ETIJ(K,1,0,1) = HPI(K)
  705                CONTINUE
                  ELSE IF (IDL .EQ. 1) THEN
                     DO 755 K = 1, NPP12
                        ETIJ(K,0,0,1) = EXP12(K,1)
  755                CONTINUE
                  ELSE IF ((IDL .GT. 1) .AND.
     *                     (IBTAND(1,IDL) .EQ. 0)) THEN
                     DO 765 K = 1, NPP12
                        ETIJ(K,1,0,1) = HPI(K)*ETIJ(K,0,0,0)
  765                CONTINUE
                  ELSE
                     DO 775 K = 1, NPP12
                        ETIJ(K,0,0,1) = DLEVEL*EXP12(K,1)*ETIJ(K,0,0,0)
  775                CONTINUE
                  ENDIF
C
C              E(1,1) AND E(0,2)
C
               ELSE IF (IJ .EQ. 2) THEN
                  IF (IDL .EQ. 0) THEN
                     DO 805 K = 1, NPP12
                        HPK = HPI(K)
                        ETIJ(K,0,I,J) = SIGN*HPK
                        ETIJ(K,2,I,J) = HPK*HPK
  805                CONTINUE
                  ELSE IF (IDL .EQ. 1) THEN
                     IF (J .EQ. 1) THEN
                        DO 855 K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = EXP12(K,1)
                           EXP2 = -EXP12(K,2)
                           ETIJ(K,1,1,1) = HPK*(EXP1 + EXP2)
  855                   CONTINUE
                     ELSE
                        DO 865 K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = D2*EXP12(K,1)
                           ETIJ(K,1,0,2) = HPK*EXP1
  865                   CONTINUE
                     ENDIF
                  ELSE IF ((IDL .GT. 1) .AND.
     *                     (IBTAND(1,IDL) .EQ. 0)) THEN
                     IF (J .EQ. 1) THEN
                        DO 856 K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,2,1,1) =    HPK*ETIJ(K,1,1,0)
                           ETIJ(K,0,1,1) = EXP1*ETIJM1(K,0,1,0)
     *                                   +   SIGN*ETIJ(K,1,1,0)
  856                   CONTINUE
                     ELSE
                        DO 866 K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,2,0,2) =    HPK*ETIJ(K,1,0,1)
                           ETIJ(K,0,0,2) = EXP1*ETIJM1(K,0,0,1)
     *                                   +   SIGN*ETIJ(K,1,0,1)
  866                   CONTINUE
                     ENDIF
                  ELSE
                     IF (J .EQ. 1) THEN
                        DO 857 K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,1,1,1) =    HPK*ETIJ(K,0,1,0)
     *                                   + EXP1*ETIJM1(K,1,1,0)
  857                   CONTINUE
                     ELSE
                        DO 867 K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,1,0,2) =    HPK*ETIJ(K,0,0,1)
     *                                   + EXP1*ETIJM1(K,1,0,1)
  867                   CONTINUE
                     ENDIF
                  ENDIF
C
C              E(I,J)
C
               ELSE
                  IF ((IDL .EQ. 0) .OR.
     *               ((IDL .GT. 1) .AND. (IBTAND(1,IDL) .EQ. 0))) THEN
                     IF (IBTAND(1,IJ) .EQ. 0) THEN
                        DO 905 K = 1, NPP12
                           ETIJ(K, 0,I,J) =   SIGN*ETIJ(K,   1,I,J-1)
                           ETIJ(K,IJ,I,J) = HPI(K)*ETIJ(K,IJ-1,I,J-1)
  905                   CONTINUE
                        IF (IDL .GT. 1) THEN
                           DO 904 K = 1, NPP12
                              EXP1 = DLEVEL*EXP12(K,1)
                              ETIJ(K,0,I,J) =        ETIJ(K,0,I,  J)
     *                                      + EXP1*ETIJM1(K,0,I,J-1)
  904                      CONTINUE
                        ENDIF
                     ELSE
                        DO 906 K = 1, NPP12
                           ETIJ(K,IJ,I,J) = HPI(K)*ETIJ(K,IJ-1,I,J-1)
  906                   CONTINUE
                     END IF
                     DO 915 T = 2 - IBTAND(1,IJ), IJ - 2, 2
                        T1 = SIGN*SDPRE(T + 1)
                        DO 925 K = 1, NPP12
                           ETIJ(K,T,I,J) = HPI(K)*ETIJ(K,T-1,I,J-1)
     &                                       + T1*ETIJ(K,T+1,I,J-1)
  925                   CONTINUE
                        IF (IDL .GT. 1) THEN
                           DO 924 K = 1, NPP12
                              EXP1 = DLEVEL*EXP12(K,1)
                              ETIJ(K,T,I,J) =        ETIJ(K,T,I,  J)
     *                                      + EXP1*ETIJM1(K,T,I,J-1)
  924                      CONTINUE
                        ENDIF
  915                CONTINUE
                  ELSE
                     IF (IBTAND(1,IJ) .EQ. 0) THEN
                        DO 955 K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,IJ-1,I,J) =    HPK*ETIJ(K,IJ-2,I,J-1)
     *                                      + EXP1*ETIJM1(K,IJ-1,I,J-1)
  955                   CONTINUE
                     ELSE
                        DO 956 K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,IJ-1,I,J) =    HPK*ETIJ(K,IJ-2,I,J-1)
     *                                      + EXP1*ETIJM1(K,IJ-1,I,J-1)
                           ETIJ(K,   0,I,J) = EXP1*ETIJM1(K,   0,I,J-1)
     *                                      +   SIGN*ETIJ(K,   1,I,J-1)
  956                   CONTINUE
                     END IF
                     DO 965 T = 1 + IBTAND(1,IJ), IJ - 3, 2
                        T1 = SIGN*SDPRE(T + 1)
                        DO 975 K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,T,I,J) =    HPK*ETIJ(K,T-1,I,J-1)
     *                                   + EXP1*ETIJM1(K,  T,I,J-1)
     *                                   +     T1*ETIJ(K,T+1,I,J-1)
  975                   CONTINUE
  965                CONTINUE
                  ENDIF
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
         WRITE (LUPRI,'(2X,A,3I5)') 'IDL, IABO, X', IDL,IAB0,X
         WRITE (LUPRI,'(2X,A,2I5)') 'MAXI,MAXJ   ', MAXI, MAXJ
         IF (IPRINT .GT. 20) THEN
            IF (X .EQ. 1) WORD = 'EX00(I,J,T)'
            IF (X .EQ. 2) WORD = 'EY00(I,J,T)'
            IF (X .EQ. 3) WORD = 'EZ00(I,J,T)'
            DO 1000 I = 0, MAXI
            DO 1000 J = 0, MAXJ
            DO 1000 T = 0, I + J
            IF (IBTAND(I + J - T + IDL,IAB0) .EQ. 0) THEN
               WRITE (LUPRI,'(/,2X,A4,A1,I1,A1,I1,A1,I1,A1,/)')
     &              WORD, '(', I, ',', J, ', ',  T, ')'
               WRITE (LUPRI,'(1X,6F12.8)') (ETIJ(K,T,I,J),K=1,NPP12)
            END IF
 1000       CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck eriscl */
      SUBROUTINE ERISCL(ETIJ,EXP12,SCALE,ISCAL1,ISCAL2,NPP12,
     &                  JMAX1,JMAX2,MAX1,MAX2,IAB0,IPRINT)
#include <implicit.h>
#include <priunit.h>
      CHARACTER WORD*4
      INTEGER T, X, IAB0(3)
      DIMENSION ETIJ(NPP12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3),
     &          SCALE(NPP12), EXP12(NPP12,3)
#include <ibtfun.h>
C
      DO 100 K = 1, NPP12
         SCALE(K) = ((EXP12(K,1)/EXP12(K,3))**ISCAL1)
     &            * ((EXP12(K,2)/EXP12(K,3))**ISCAL2)
  100 CONTINUE
C
      DO 200 J = 0, MAX2
      DO 200 I = 0, MAX1
      DO 200 T = 0, I + J
      DO 200 K = 1, NPP12
         ETIJ(K,T,I,J,1) = SCALE(K)*ETIJ(K,T,I,J,1)
  200 CONTINUE
C
C     Print
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERISCL',-1)
         WRITE (LUPRI,'(/,2X,A,2I5)') ' Scale powers: ', ISCAL1,ISCAL2
         CALL HEADER('SCALE in ERISCL',1)
         CALL OUTPUT(SCALE,1,1,1,NPP12,1,NPP12,1,LUPRI)
         IF (IPRINT .GT. 20) THEN
            CALL HEADER('ETIJ in ERISCL',-1)
            DO 300 X = 1, 3
               IF (X .EQ. 1) WORD = 'EX00(I,J,T)'
               IF (X .EQ. 2) WORD = 'EY00(I,J,T)'
               IF (X .EQ. 3) WORD = 'EZ00(I,J,T)'
               DO 400 I = 0, MAX1
               DO 400 J = 0, MAX2
               DO 400 T = 0, I + J
                  IF (IBTAND(I + J - T,IAB0(X)) .EQ. 0) THEN
                    WRITE (LUPRI,'(/,2X,A4,A1,I1,A1,I1,A1,I1,A1,/)')
     &                   WORD, '(', I, ',', J, ', ',  T, ')'
                    WRITE (LUPRI,'(1X,6F12.8)')
     &                   (ETIJ(K,T,I,J,X),K=1,NPP12)
                  END IF
  400          CONTINUE
  300       CONTINUE
         END IF
      END IF
C
      RETURN
      END
