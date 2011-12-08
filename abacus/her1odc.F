!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
C  /* Deck getodc */
      SUBROUTINE GETODC(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC,KINODC,
     &                  ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,WORK,LWORK,
     &                  CORPX,CORPY,CORPZ,DONUC1,DOMOM1,ORIGIN,INTTYP)
#include <implicit.h>
#include <priunit.h>
      LOGICAL DIFODC, KINODC, ONECEN, DONUC1, DOMOM1
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          WORK(LWORK), ORIGIN(3)
C
      IF (DONUC1) THEN
         JMAXAP = JMAXA + JMAXD
         JMAXBP = JMAXB
      ELSE
         JMAXAP = JMAXA
         JMAXBP = JMAXB + JMAXD
      END IF
      IF (DOMOM1) THEN
         JMAXAP = JMAXAP + JMAXM
      ELSE
         JMAXBP = JMAXBP + JMAXM
      END IF
      ITEX = MAX(1,JMAXD + JMAXM)
      JMAXTP = JMAXT + ITEX
C
      IF ((INTTYP .EQ. 49) .OR. (INTTYP .EQ. 51)) THEN
         IF (DOMOM1) THEN
            JMAXBP = JMAXBP + 1
         ELSE
            JMAXAP = JMAXAP + 1
         END IF
         JMAXTP = JMAXTP +1
      ENDIF
C
      KODCUN = 1
      KLAST  = 1 + 3*(JMAXAP + 4)*(JMAXBP + 3)*(JMAXTP + 3)*(JMAXD + 1)
     &              *(JMAXM + 1)
      IF (KLAST .GT. LWORK) CALL STOPIT('GETODC',' ',KLAST,LWORK)
      CALL GETOD1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,JMAXAP,JMAXBP,
     &            JMAXTP,WORK(KODCUN),DIFODC,KINODC,ONECEN,EXPA,EXPB,
     &            IPRINT,SAAB13,EXPPI,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &            ORIGIN,INTTYP)
      RETURN
      END
C  /* Deck getod1 */
      SUBROUTINE GETOD1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,JMAXAP,JMAXBP,
     &                  JMAXTP,ODCUND,DIFODC,KINODC,ONECEN,EXPA,EXPB,
     &                  IPRINT,SAAB13,EXPPI,CORPX,CORPY,CORPZ,DONUC1,
     &                  DOMOM1,ORIGIN,INTTYP)
#include <implicit.h>
#include <priunit.h>
      LOGICAL DIFODC, KINODC, ONECEN, DONUC1, DOMOM1
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &      ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,0:JMAXM,3)
C
C     Clear arrays
C
      CALL DZERO(ODC,3*(JMAXA  + 1)*(JMAXB + 1)*(JMAXT + 1)*(JMAXD + 1)
     &                *(JMAXM + 1))
      CALL DZERO(ODCUND,3*(JMAXAP + 4)*(JMAXBP + 3)*(JMAXTP+3)
     &                   *(JMAXD+1)*(JMAXM + 1))
C
C     Undifferentiated expansion coefficients
C
      CALL ONEODC(ODCUND,JMAXAP,JMAXBP,JMAXTP,JMAXD,JMAXM,SAAB13,EXPPI,
     &            CORPX,CORPY,CORPZ,IPRINT)
      CALL COPODC(ODCUND,ODC,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,JMAXTP,
     &            JMAXD,JMAXM,0,0)
      IF (IPRINT .GE. 10) THEN
         CALL TITLER('Output from GETOD1','*',103)
         NROW = (JMAXA + 1)*(JMAXB + 1)
         NCOL = JMAXT + 1
         CALL AROUND('Undifferentiated ODC in GETOD1 - x component')
         CALL OUTPUT(ODC(0,0,0,0,0,1),1,NROW,1,NCOL,NROW,NCOL,1,LUPRI)
         CALL AROUND('Undifferentiated ODC in GETOD1 - y component')
         CALL OUTPUT(ODC(0,0,0,0,0,2),1,NROW,1,NCOL,NROW,NCOL,1,LUPRI)
         CALL AROUND('Undifferentiated ODC in GETOD1 - z component')
         CALL OUTPUT(ODC(0,0,0,0,0,3),1,NROW,1,NCOL,NROW,NCOL,1,LUPRI)
      END IF
C
C     Expansion coefficients for derivatives
C
      IF (DIFODC) THEN
         IF (.NOT.ONECEN) THEN
            IF (DONUC1) THEN
               CALL DODCA(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                    JMAXTP,JMAXD,JMAXM,EXPA,IPRINT)
            ELSE
               CALL DODCB(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                    JMAXTP,JMAXD,JMAXM,EXPB,IPRINT)
            END IF
            IF (INTTYP .EQ. 51) CALL DEFODC(ODC,ODCUND,JMAXA,JMAXAP,
     $           JMAXB,JMAXBP,JMAXT,JMAXTP,JMAXD,JMAXM,EXPA,EXPB,
     $           IPRINT)
         END IF
      ELSE IF (KINODC) THEN
            CALL TODC(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                JMAXTP,JMAXD,JMAXM,EXPA)
      END IF
C
C     Expansion coeffiecients for moments or electric derivatives
C
      IF (JMAXM .GT. 0) THEN
         IF ((INTTYP .EQ. 49) .OR. (INTTYP .EQ. 51)) THEN
            CALL EFODC(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                 JMAXTP,JMAXD,JMAXM,EXPA,EXPB,IPRINT)
         ELSE IF (INTTYP .EQ. 63 .OR. INTTYP .EQ. 83) THEN
            CALL DODCB2(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                  JMAXTP,JMAXD,JMAXM,EXPB,IPRINT)
         ELSE IF (DOMOM1) THEN
            CALL MODCA(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                 JMAXTP,JMAXD,JMAXM,ORIGIN,IPRINT,INTTYP)
         ELSE
            CALL MODCB(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                 JMAXTP,JMAXD,JMAXM,ORIGIN,IPRINT,INTTYP)
         END IF
      END IF
      RETURN
      END
C  /* Deck copodc */
      SUBROUTINE COPODC(ODCUND,ODC,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,
     &                JMAXTP,JMAXD,JMAXM,ICOPYD,ICOPYM)
C
C     Copy expansion coefficients from ODCUND to ODC
C
#include <implicit.h>
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,
     &                 0:JMAXD,0:JMAXM,3),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
      DO 100 IC = 1, 3
         DO 100 IM = 0, ICOPYM
         DO 100 ID = 0, ICOPYD
            DO 100 IT = 0, JMAXT
               DO 100 IB = 0, JMAXB
               DO 100 IA = 0, JMAXA
                  ODC(IA,IB,IT,ID,IM,IC) = ODCUND(IA,IB,IT,ID,IM,IC)
  100 CONTINUE
      RETURN
      END
C  /* Deck oneodc */
      SUBROUTINE ONEODC(ODCUND,JMAXAP,JMAXBP,JMAXTP,JMAXD,JMAXM,FAC,
     &                  EXPPI,CORPX,CORPY,CORPZ,IPRINT)
C
C     TUH 91
C
#include <implicit.h>
#include <priunit.h>
      INTEGER T, A, B, AB
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,
     &                 0:JMAXD,0:JMAXM,3)
#include <onecom.h>
C
      PAX = CORPX - CORAX
      PAY = CORPY - CORAY
      PAZ = CORPZ - CORAZ
      PBX = CORPX - CORBX
      PBY = CORPY - CORBY
      PBZ = CORPZ - CORBZ
      EXPPIH = EXPPI/2
      DO 100 A = 0, JMAXAP
         IF (A .EQ. 0) THEN
            ODCUND(0,0,0,0,0,1) = FAC
            ODCUND(0,0,0,0,0,2) = FAC
            ODCUND(0,0,0,0,0,3) = FAC
         ELSE
            DO 200 T = 0, A
               X_TP1 = T + 1
               ODCUND(A,0,T,0,0,1) = EXPPIH*ODCUND(A-1,0,T-1,0,0,1)
     *                                + PAX*ODCUND(A-1,0,T  ,0,0,1)
     *                              + X_TP1*ODCUND(A-1,0,T+1,0,0,1)
               ODCUND(A,0,T,0,0,2) = EXPPIH*ODCUND(A-1,0,T-1,0,0,2)
     *                                + PAY*ODCUND(A-1,0,T  ,0,0,2)
     *                              + X_TP1*ODCUND(A-1,0,T+1,0,0,2)
               ODCUND(A,0,T,0,0,3) = EXPPIH*ODCUND(A-1,0,T-1,0,0,3)
     *                                + PAZ*ODCUND(A-1,0,T  ,0,0,3)
     *                              + X_TP1*ODCUND(A-1,0,T+1,0,0,3)
  200       CONTINUE
         END IF
         DO 300 B = 1, JMAXBP
            AB = A + B
C
            DO 400 T = 0, AB
               X_TP1 = T + 1
               ODCUND(A,B,T,0,0,1) = EXPPIH*ODCUND(A,B-1,T-1,0,0,1)
     &                                + PBX*ODCUND(A,B-1,T  ,0,0,1)
     &                              + X_TP1*ODCUND(A,B-1,T+1,0,0,1)
               ODCUND(A,B,T,0,0,2) = EXPPIH*ODCUND(A,B-1,T-1,0,0,2)
     &                                + PBY*ODCUND(A,B-1,T  ,0,0,2)
     &                              + X_TP1*ODCUND(A,B-1,T+1,0,0,2)
               ODCUND(A,B,T,0,0,3) = EXPPIH*ODCUND(A,B-1,T-1,0,0,3)
     &                                + PBZ*ODCUND(A,B-1,T  ,0,0,3)
     &                              + X_TP1*ODCUND(A,B-1,T+1,0,0,3)
  400          CONTINUE
  300    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 20) THEN
         CALL TITLER('Output from ONEODC','*',103)
         NROW = (JMAXAP + 4)*(JMAXBP + 3)
         NCOL = JMAXTP + 3
C
C     NOTICE: Prints all ODCUND, including dummy elements.
C
         CALL AROUND('ODCUND in ONEODC - x component')
         CALL OUTPUT(ODCUND(-2,0,-2,0,0,1),1,NROW,1,NCOL,NROW,NCOL,
     &               1,LUPRI)
         CALL AROUND('ODCUND in ONEODC - y component')
         CALL OUTPUT(ODCUND(-2,0,-2,0,0,2),1,NROW,1,NCOL,NROW,NCOL,
     &               1,LUPRI)
         CALL AROUND('ODCUND in ONEODC - z component')
         CALL OUTPUT(ODCUND(-2,0,-2,0,0,3),1,NROW,1,NCOL,NROW,NCOL,
     &               1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck dodca */
      SUBROUTINE DODCA(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                 JMAXTP,JMAXD,JMAXM,EXPA,IPRINT)
C
C     TUH, Recursive calculation of coefficients implemented dec.91, KR
C
#include <implicit.h>
#include <priunit.h>
      INTEGER A, B, D, T, X
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,
     &                 0:JMAXD,0:JMAXM,3),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
C
      TEXPA1 = EXPA + EXPA
      DO 10 D = 1, JMAXD
         DO 20 A = 0, JMAXAP
            FAC10A = A
            DO 30 B = 0, JMAXBP
               DO 40 T = 0, A + B + D
               DO 40 X = 1, 3
                  ODCUND(A,B,T,D,0,X) = TEXPA1*ODCUND(A + 1,B,T,D-1,0,X)
     &                                - FAC10A*ODCUND(A - 1,B,T,D-1,0,X)
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      CALL COPODC(ODCUND,ODC,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,JMAXTP,
     &            JMAXD,JMAXM,JMAXD,0)
      IF (IPRINT .GE. 20) THEN
         CALL TITLER('Output from DODCA','*',103)
         JMAXAB = (JMAXA + 1)*(JMAXB + 1)
         DO 60 I = 1, JMAXD
            WRITE(LUPRI,'(1X,I2,A)') I, '. order differentiated'//
     &           ' coefficients'
            CALL AROUND('ODC in DODCA - x component')
            CALL OUTPUT(ODC(0,0,0,I,0,1),1,JMAXAB,1,JMAXT,JMAXAB,JMAXT,
     &               1,LUPRI)
            CALL AROUND('ODC in DODCA - y component')
            CALL OUTPUT(ODC(0,0,0,I,0,2),1,JMAXAB,1,JMAXT,JMAXAB,JMAXT,
     &               1,LUPRI)
            CALL AROUND('ODC in DODCA - z component')
            CALL OUTPUT(ODC(0,0,0,I,0,3),1,JMAXAB,1,JMAXT,JMAXAB,JMAXT,
     &               1,LUPRI)
 60      CONTINUE
      END IF
      RETURN
      END
C  /* Deck dodcb */
      SUBROUTINE DODCB(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                JMAXTP,JMAXD,JMAXM,EXPB,IPRINT)
C
C     TUH
C
#include <implicit.h>
#include <priunit.h>
      INTEGER A, B, D, T, X
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,
     &                 0:JMAXD,0:JMAXM,3),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
C
      TEXPB1 = EXPB + EXPB
      DO 10 D = 1, JMAXD
         DO 20 A = 0, JMAXAP
            DO 30 B = 0, JMAXBP
               FAC10B = B
               DO 40 T = 0, A + B + D
               DO 40 X = 1, 3
                  ODCUND(A,B,T,D,0,X) = TEXPB1*ODCUND(A,B + 1,T,D-1,0,X)
     &                                - FAC10B*ODCUND(A,B - 1,T,D-1,0,X)
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      CALL COPODC(ODCUND,ODC,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,JMAXTP,
     &            JMAXD,JMAXM,JMAXD,0)
      IF (IPRINT .GE. 20) THEN
         CALL TITLER('Output from DODCB','*',103)
         JMAXAB = (JMAXA + 1)*(JMAXB + 1)
         DO 60 I = 1, JMAXD
            WRITE(LUPRI,'(1X,I2,A)') I, '. order differentiated'//
     &           ' coefficients'
            CALL AROUND('ODC in DODCB - x component')
            CALL OUTPUT(ODC(0,0,0,I,0,1),1,JMAXAB,1,JMAXT,JMAXAB,JMAXT,
     &               1,LUPRI)
            CALL AROUND('ODC in DODCB - y component')
            CALL OUTPUT(ODC(0,0,0,I,0,2),1,JMAXAB,1,JMAXT,JMAXAB,JMAXT,
     &               1,LUPRI)
            CALL AROUND('ODC in DODCB - z component')
            CALL OUTPUT(ODC(0,0,0,I,0,3),1,JMAXAB,1,JMAXT,JMAXAB,JMAXT,
     &               1,LUPRI)
 60      CONTINUE
      END IF
      RETURN
      END
C  /* Deck modca */
      SUBROUTINE MODCA(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                 JMAXTP,JMAXD,JMAXM,ORIGIN,IPRINT,INTTYP)
C
C     K.Ruud Nov 1991
C
#include <implicit.h>
#include <priunit.h>
      INTEGER A, B, T, D, X
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,
     &                 0:JMAXD,0:JMAXM,3), ORIGIN(3),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include <onecom.h>
      OAX = CORAX - ORIGIN(1)
      OAY = CORAY - ORIGIN(2)
      OAZ = CORAZ - ORIGIN(3)
      DO 10 M = 1, JMAXM
         DO 20 A = 0, JMAXAP
         DO 20 B = 0, JMAXBP
         DO 20 D = 0, JMAXD
            DO 30 T = 0, A + B + D + M
               ODCUND(A,B,T,D,M,1) = ODCUND(A + 1,B,T,D,M-1,1)
     &                             + OAX*ODCUND(A,B,T,D,M-1,1)
               ODCUND(A,B,T,D,M,2) = ODCUND(A + 1,B,T,D,M-1,2)
     &                             + OAY*ODCUND(A,B,T,D,M-1,2)
               ODCUND(A,B,T,D,M,3) = ODCUND(A + 1,B,T,D,M-1,3)
     &                             + OAZ*ODCUND(A,B,T,D,M-1,3)
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      CALL COPODC(ODCUND,ODC,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,JMAXTP,
     &            JMAXD,JMAXM,JMAXD,JMAXM)
      IF (IPRINT .GE. 20) THEN
         CALL TITLER('Output from MODCA','*',103)
         JMAXAB = (JMAXA + 1)*(JMAXB + 1)
         DO 40 D = 0, JMAXD
            DO 40 M = 0, JMAXM
               WRITE(LUPRI,'(/A/)') ' '
               WRITE(LUPRI,'(10X,I2,A)')
     &              D, '.order differentiated coefficients'
               WRITE(LUPRI,'(10X,I2,A)')
     &              M, '.order moment'
               CALL AROUND('ODC in MODCA - x component')
               CALL OUTPUT(ODC(0,0,0,D,M,1),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
               CALL AROUND('ODC in MODCA - y component')
               CALL OUTPUT(ODC(0,0,0,D,M,2),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
               CALL AROUND('ODC in MODCA - z component')
               CALL OUTPUT(ODC(0,0,0,D,M,3),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
 40      CONTINUE
      END IF
      RETURN
      END
C  /* Deck modcb */
      SUBROUTINE MODCB(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                 JMAXTP,JMAXD,JMAXM,ORIGIN,IPRINT,INTTYP)
C
C     K.Ruud, Aug 1992
C
#include <implicit.h>
#include <priunit.h>
      INTEGER A, B, T, D, X
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,
     &                 0:JMAXD,0:JMAXM,3),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),ORIGIN(3)
#include <onecom.h>
C
      OBX = CORBX - ORIGIN(1)
      OBY = CORBY - ORIGIN(2)
      OBZ = CORBZ - ORIGIN(3)
      DO 10 M = 1, JMAXM
         DO 20 A = 0, JMAXAP
         DO 20 B = 0, JMAXBP
         DO 20 D = 0, JMAXD
            DO 30 T = 0, A + B + D + M
               ODCUND(A,B,T,D,M,1) = ODCUND(A,B + 1,T,D,M-1,1)
     &                             + OBX*ODCUND(A,B,T,D,M-1,1)
               ODCUND(A,B,T,D,M,2) = ODCUND(A,B + 1,T,D,M-1,2)
     &                             + OBY*ODCUND(A,B,T,D,M-1,2)
               ODCUND(A,B,T,D,M,3) = ODCUND(A,B + 1,T,D,M-1,3)
     &                             + OBZ*ODCUND(A,B,T,D,M-1,3)
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      CALL COPODC(ODCUND,ODC,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,JMAXTP,
     &            JMAXD,JMAXM,JMAXD,JMAXM)
      IF (IPRINT .GE. 20) THEN
         CALL TITLER('Output from MODCB','*',103)
         JMAXAB = (JMAXA + 1)*(JMAXB + 1)
         DO 40 D = 0, JMAXD
            DO 40 M = 0, JMAXM
               WRITE(LUPRI,'(/A/)') ' '
               WRITE(LUPRI,'(10X,I2,A)')
     &              D, '.order differentiated coefficients'
               WRITE(LUPRI,'(10X,I2,A)')
     &              M, '.order moment'
               CALL AROUND('ODC in MODCB - x component')
               CALL OUTPUT(ODC(0,0,0,D,M,1),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
               CALL AROUND('ODC in MODCB - y component')
               CALL OUTPUT(ODC(0,0,0,D,M,2),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
               CALL AROUND('ODC in MODCB - z component')
               CALL OUTPUT(ODC(0,0,0,D,M,3),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
 40      CONTINUE
      END IF
      RETURN
      END
C  /* Deck todc */
      SUBROUTINE TODC(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                JMAXTP,JMAXD,JMAXM,EXPA)
C
C     TUH
C
#include <implicit.h>
#include <priunit.h>
      INTEGER A, B, T
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,
     &                 0:JMAXM,3),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)

      TEXPA1 = EXPA + EXPA
      TEXPA2 = TEXPA1*TEXPA1
      IF (JMAXD .LT. 2) THEN
         WRITE (LUPRI,'(//,1X,A,I2,A)')
     &          ' JMAXD = ',JMAXD,' too small in TODC,',
     &          ' program cannot proceed.'
         CALL QUIT('JMAXD too small in TODC')
      END IF
      DO 100 A = 0, JMAXA
         FAC21A = (2.0D0*A + 1.0D0) *  TEXPA1
         IF (A .LT. 2) THEN
            DO 200 B = 0, JMAXB
               ODC(A,B,0,2,0,1) = TEXPA2*ODCUND(A+2,B,0,0,0,1)
     &                          - FAC21A*ODCUND(A  ,B,0,0,0,1)
               ODC(A,B,0,2,0,2) = TEXPA2*ODCUND(A+2,B,0,0,0,2)
     &                          - FAC21A*ODCUND(A  ,B,0,0,0,2)
               ODC(A,B,0,2,0,3) = TEXPA2*ODCUND(A+2,B,0,0,0,3)
     &                          - FAC21A*ODCUND(A  ,B,0,0,0,3)
  200       CONTINUE
         ELSE
            FAC20A = (A*(A - 1.0D0))
            DO 300 B = 0, JMAXB
               ODC(A,B,0,2,0,1) = TEXPA2*ODCUND(A+2,B,0,0,0,1)
     &                        - FAC21A*ODCUND(A  ,B,0,0,0,1)
     &                        + FAC20A*ODCUND(A-2,B,0,0,0,1)
               ODC(A,B,0,2,0,2) = TEXPA2*ODCUND(A+2,B,0,0,0,2)
     &                        - FAC21A*ODCUND(A  ,B,0,0,0,2)
     &                        + FAC20A*ODCUND(A-2,B,0,0,0,2)
               ODC(A,B,0,2,0,3) = TEXPA2*ODCUND(A+2,B,0,0,0,3)
     &                        - FAC21A*ODCUND(A  ,B,0,0,0,3)
     &                        + FAC20A*ODCUND(A-2,B,0,0,0,3)
  300       CONTINUE
         END IF
  100 CONTINUE
      RETURN
      END
C
C  /* Deck efodc */
      SUBROUTINE EFODC(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                JMAXTP,JMAXD,JMAXM,EXPA,EXPB,IPRINT)
C
C     Purpose: Makes expansion coefficients for electric
C     derivative of overlap. JMAXM = 1 is used for convenience
C     to store coeff. No moments involved.
C
C     H. Heiberg
C
#include <implicit.h>
#include <priunit.h>
      INTEGER A, B, T
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,
     &                 0:JMAXM,3),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include <pi.h>
      PARAMETER(D1=1.00 D00, D2=2.00 D00)
#include <onecom.h>
C
      DO 10 A = 0, JMAXA
         FA = (D2*A + 1.0D0)/(D2*EXPA)
         DO 20 B = 0, JMAXB
            FB = (D2*B + 1.0D0)/(D2*EXPB)
            DO 30 T = 0 , JMAXTP
               ODCUND(A,B,T,0,1,1) = FA*ODCUND(A + 1,B,T,0,0,1)
     &                             + FB*ODCUND(A,B + 1,T,0,0,1)
               ODCUND(A,B,T,0,1,2) = FA*ODCUND(A + 1,B,T,0,0,2)
     &                             + FB*ODCUND(A,B + 1,T,0,0,2)
               ODCUND(A,B,T,0,1,3) = FA*ODCUND(A + 1,B,T,0,0,3)
     &                             + FB*ODCUND(A,B + 1,T,0,0,3)
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      CALL COPODC(ODCUND,ODC,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,JMAXTP,
     &            JMAXD,JMAXM,JMAXD,JMAXM)
      IF (IPRINT .GE. 20) THEN
         CALL TITLER('Output from EFODC','*',103)
         JMAXAB = (JMAXA + 1)*(JMAXB + 1)
         M = 1
               CALL AROUND('ODC in EFODC - x component')
               CALL OUTPUT(ODC(0,0,0,0,M,1),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
               CALL AROUND('ODC in EFODC - y component')
               CALL OUTPUT(ODC(0,0,0,0,M,2),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
               CALL AROUND('ODC in EFODC - z component')
               CALL OUTPUT(ODC(0,0,0,0,M,3),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
      END IF
      RETURN
      END
C
C  /* Deck defodc */
      SUBROUTINE DEFODC(ODC,ODCUND,JMAXA,JMAXAP,JMAXB,JMAXBP,JMAXT,
     $           JMAXTP,JMAXD,JMAXM,EXPA,EXPB,IPRINT)
C
C     Modified differentiated coeffitions.
C
#include <implicit.h>
#include <priunit.h>
      INTEGER A, B, D
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,
     &                 0:JMAXM,3),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include <pi.h>
      PARAMETER(D1=1.00 D00, D2=2.00 D00)
#include <onecom.h>
C
      DO 10 A = 0, JMAXA
         FA = (D2*A + 1.0D0)/(D2*EXPA)
         DO 20 B = 0, JMAXB
            FB = (D2*B + 1.0D0)/(D2*EXPB)
               ODCUND(A,B,0,2,1,1) = FA*ODCUND(A + 1,B,0,2,0,1)
     &                             + FB*ODCUND(A,B + 1,0,2,0,1)
               ODCUND(A,B,0,2,1,2) = FA*ODCUND(A + 1,B,0,2,0,2)
     &                             + FB*ODCUND(A,B + 1,0,2,0,2)
               ODCUND(A,B,0,2,1,3) = FA*ODCUND(A + 1,B,0,2,0,3)
     &                             + FB*ODCUND(A,B + 1,0,2,0,3)
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      CALL COPODC(ODCUND,ODC,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,JMAXTP,
     &            JMAXD,JMAXM,JMAXD,JMAXM)
      IF (IPRINT .GE. 20) THEN
         CALL TITLER('Output from DEFODC','*',103)
         JMAXAB = (JMAXA + 1)*(JMAXB + 1)
               WRITE(LUPRI,'(/A/)') ' '
               WRITE(LUPRI,'(10X,A)')
     &              '2.order differentiated coefficients'
               CALL AROUND('ODC in DEFODC - x component')
               CALL OUTPUT(ODC(0,0,0,2,1,1),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
               CALL AROUND('ODC in DEFODC - y component')
               CALL OUTPUT(ODC(0,0,0,2,1,2),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
               CALL AROUND('ODC in DEFODC - z component')
               CALL OUTPUT(ODC(0,0,0,2,1,3),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
      END IF
      RETURN
      END
      SUBROUTINE DODCB2(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                 JMAXTP,JMAXD,JMAXM,EXPB,IPRINT)
C
C     K.Ruud, Sept 2000. Code for B-derivative integrals at the same
C     time as we have A-derivative integrals. Needed for pVp integrals
C     of the Douglas-Kroll transformation.
C
#include <implicit.h>
#include <priunit.h>
      INTEGER A, B, T, D, X, M
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,
     &                 0:JMAXD,0:JMAXM,3),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
#include <onecom.h>
C
      TEXPB1 = EXPB + EXPB
      DO 10 M = 1, JMAXM
         DO 20 A = 0, JMAXAP
            DO 30 B = 0, JMAXBP
               DO 40 D = 0, JMAXD
               FAC10B = B
               DO 40 T = 0, A + B + D + M
               DO 40 X = 1, 3
                  ODCUND(A,B,T,D,M,X) = TEXPB1*ODCUND(A,B + 1,T,D,M-1,X)
     &                                - FAC10B*ODCUND(A,B - 1,T,D,M-1,X)
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
C
      CALL COPODC(ODCUND,ODC,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,JMAXTP,
     &            JMAXD,JMAXM,JMAXD,JMAXM)
      IF (IPRINT .GE. 20) THEN
         CALL TITLER('Output from DODCB2','*',103)
         JMAXAB = (JMAXA + 1)*(JMAXB + 1)
         DO 50 D = 0, JMAXD
            DO 50 M = 0, JMAXM
               WRITE(LUPRI,'(/A/)') ' '
               WRITE(LUPRI,'(10X,I2,A)')
     &              D, '.order differentiated coefficients'
               WRITE(LUPRI,'(10X,I2,A)')
     &              M, '.order moment'
               CALL AROUND('ODC in MODCB - x component')
               CALL OUTPUT(ODC(0,0,0,D,M,1),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
               CALL AROUND('ODC in MODCB - y component')
               CALL OUTPUT(ODC(0,0,0,D,M,2),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
               CALL AROUND('ODC in MODCB - z component')
               CALL OUTPUT(ODC(0,0,0,D,M,3),1,JMAXAB,1,JMAXT,JMAXAB,
     &                     JMAXT,1,LUPRI)
 50      CONTINUE
      END IF
      RETURN
      END
