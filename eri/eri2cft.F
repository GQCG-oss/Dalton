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
C  /* Deck expcft */
      SUBROUTINE EXPCFT(ETIJ,NPP12,JMAX1,JMAX2,COOR12,EXP12,
     &                  IAB0,HPI,PA,PB,IRUTIN,IELCTR,MAX1,MAX2,
     &                  NCNT12,IDL,ISCAL1,ISCAL2,IPRINT,
     &                  TEXPA,TEXPB,AOVERP,BOVERP)
C
C     T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
#include <r12int.h>
      PARAMETER (D1 = 1.0D0)
      INTEGER X, T
      DIMENSION ETIJ(NPP12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3,*),
     &          COOR12(NPP12,3), EXP12(NPP12,3), IAB0(3), 
     &          HPI(NPP12), PA(NPP12,3), PB(NPP12,3),
     &          TEXPA(NPP12),TEXPB(NPP12),AOVERP(NPP12),BOVERP(NPP12)
C     TEXPA, TEXPB, AOVERP and BOVERP have been added; EQ30LB contains
C     labels for coefficients of subroutine WKEQ30 (WK/UniKA/04-11-2002).
      CHARACTER*4 EQ30LB(3)
      DATA EQ30LB/'EX01','EY01','EZ01'/
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
C     Compute TEXPA (= 2a) and TEXPB (= 2b) (WK/UniKA/04-11-2002).
      CALL EXPVEC(HPI,PA,PB,COOR12,EXP12,SIGN,NPP12,MAX1,MAX2,IDL,
     &            IPRINT,TEXPA,TEXPB)
C
C     Expansion coefficients
C     ======================
C
      DO X = 1, 3
         CALL ERIECF(ETIJ(1,0,0,0,X,1),ETIJ(1,0,0,0,X,2),
     &               ETIJ(1,0,0,0,X,2),NPP12,MAX1,MAX2,JMAX1,JMAX2,
     &               SIGN,IAB0(X),PA(1,X),PB(1,X),EXP12,HPI,
     &               X,0,IPRINT)
         IF (U12INT) then
C           Calculate expansion coefficients for [T1,r12] integrals (WK/UniKA/04-11-2002).
            CALL WKEQ30(1,ETIJ(1,0,0,0,X,1),ETIJ(1,0,0,0,X,2),
     &                  ETIJ(1,0,0,0,X,1),
     &                  JMAX1,JMAX2,NPP12,MAX1,MAX2,SIGN,
     &                  IAB0(X),PA(1,X),PB(1,X),TEXPA,TEXPB,AOVERP,
     &                  BOVERP,HPI,EQ30LB(X),IPRINT)
            IF(XMIADR.GT.0) THEN
C             Calculates exp. coeff. for [T1,r12] integrals derivatives.
c             C. Villani, Uni-Ka, Feb 2005
              CALL WKEQ30(2,ETIJ(1,0,0,0,X,2),ETIJ(1,0,0,0,X,3),
     &                    ETIJ(1,0,0,0,X,1),
     &                    JMAX1,JMAX2,NPP12,MAX1,MAX2,SIGN,
     &                    IAB0(X),PA(1,X),PB(1,X),TEXPA,TEXPB,AOVERP,
     &                    BOVERP,HPI,EQ30LB(X),IPRINT)
            end if
         end if
      end do
      IF ((U12INT.AND.(XMIADR.EQ.0))) THEN
C       Calculates (a-b)/(a+b) for [T1,r12] integrals (WK/UniKA/04-11-2002).
        CALL AMINBP(AOVERP,NPP12,TEXPA,TEXPB,HPI,SIGN,IPRINT)
      ELSE IF(U12INT) THEN
C       Calculates (a)/(a+b) or -(b)/(a+b) needed only for
c       [T1,r12] integral derivatives, see subroutine CMIASB
c       C. Villani, Uni-Ka, Feb 2005
        IF(AENONB) THEN
          DO I=1,NPP12
            AOVERP(I)=EXP12(I,1)
          END DO
        ELSE
          DO I=1,NPP12
            AOVERP(I)=-EXP12(I,2)
          END DO
        END IF
      END IF
      IF (MIAU12) THEN
        IF (IDL.NE.0) stop 'idl.ne.0 and miau12.eq.true.!!'
        IF (U12INT) stop 'Miau12 and u12int!'
c       for half computed [T1,r12] integrals, multiply x computed expantion
c       coefficients for the (a-b)/(a+b) factor
c       C. Villani, Uni-Ka, Feb 2005
        DO J = 0, MAX2
          DO I = 0, MAX1
            DO T = 0, I + J
              DO K = 1, NPP12
                ETIJ(K,T,I,J,1,1) = (EXP12(K,1)-EXP12(K,2))
     .                              * ETIJ(K,T,I,J,1,1)
              END DO
            END DO
          END DO
        END DO
      END IF
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
         DO X = 1, 3
            CALL ERIECF(ETIJ(1,0,0,0,X,2),ETIJ(1,0,0,0,X,1),
     &                  ETIJ(1,0,0,0,X,1),NPP12,MAX1,MAX2,JMAX1,JMAX2,
     &                  SIGN,IAB0(X),PA(1,X),PB(1,X),EXP12,HPI,
     &                  X,1,IPRINT)
         END DO
      ENDIF
C
      IF (IDL .GT. 1) STOP 'No higher than first derivative yet'
C
      RETURN
      END
C  /* Deck EXPVEC */
      SUBROUTINE EXPVEC(HPI,PA,PB,COOR12,EXP12,SIGN,
     &                  NPP12,MAX1,MAX2,IDL,IPRINT,TEXPA,TEXPB)
C     Compute TEXPA and TEXPB (WK/UniKA/04-11-2002).
#include <implicit.h>
#include <priunit.h>
#include <r12int.h>
      PARAMETER (DP5 = 0.5D0, D2 = 2.0D0)
      DIMENSION HPI(NPP12), 
     &          PA(NPP12,3), PB(NPP12,3),
     &          COOR12(NPP12,3), EXP12(NPP12,3),
     &          TEXPA(NPP12), TEXPB(NPP12)
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
      DO I = 1, NPP12
         HPI(I) = HPIFAC*EXP12(I,3)
      END DO
C
C     2A, 2B
C     ======
C
      IF (U12INT) THEN
C        2a and 2b are needed for [T1,r12] integrals (WK/UniKA/04-11-2002).
         DO I = 1, NPP12
            TEXPA(I) = D2*EXP12(I,1) / EXP12(I,3)
            TEXPB(I) = D2*EXP12(I,2) / EXP12(I,3)
         END DO
      END IF
C
C     PA and PB
C     =========
C
      IF ((MAX1 .GT. 0 .AND. MAX2 .GT. 0) .OR. U12INT) THEN
C        Are always needed for [T1,r12] integrals (WK/UniKA/04-11-2002).
         DO I = 1, NPP12
            EXP1    =  EXP12(I,1)
            EXP2    = -EXP12(I,2)
            PA(I,1) = EXP2*COOR12(I,1)
            PA(I,2) = EXP2*COOR12(I,2)
            PA(I,3) = EXP2*COOR12(I,3)
            PB(I,1) = EXP1*COOR12(I,1)
            PB(I,2) = EXP1*COOR12(I,2)
            PB(I,3) = EXP1*COOR12(I,3)
         END DO
      ELSE
         IF (MAX1 .GT. 0 .OR. IDL.GT.0) THEN
            DO I = 1, NPP12
               EXP2    = -EXP12(I,2)
               PA(I,1) = EXP2*COOR12(I,1)
               PA(I,2) = EXP2*COOR12(I,2)
               PA(I,3) = EXP2*COOR12(I,3)
            END DO
         END IF
         IF (MAX2 .GT. 0) THEN
            DO I = 1, NPP12
               EXP1    = EXP12(I,1)
               PB(I,1) = EXP1*COOR12(I,1)
               PB(I,2) = EXP1*COOR12(I,2)
               PB(I,3) = EXP1*COOR12(I,3)
            END DO
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
     &                  SIGN,IAB0,PA,PB,EXP12,HPI,X,IDL,IPRINT)
C
C     A. Halkier & T. Helgaker Jan. 1999: Generalised Version 
C     that handles derivative coefficients, where IDL gives the order
C     of differentiation.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.00 D00, D2 = 2.00 D00)
      INTEGER T, X
      CHARACTER WORD*4
      DIMENSION PA(NPP12), PB(NPP12), EXP12(NPP12,3), HPI(NPP12)
      DIMENSION ETIJ(NPP12,0:JMAXI+JMAXJ,0:JMAXI,0:JMAXJ)
      DIMENSION ETIJM1(NPP12,0:JMAXI+JMAXJ,0:JMAXI,0:JMAXJ)
      DIMENSION ETIJM2(NPP12,0:JMAXI+JMAXJ,0:JMAXI,0:JMAXJ)
#include <ibtfun.h>
C
      DLEVEL = IDL
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
         DO I = 0, MAXI
C
C           E(0,0)
C
            IF (I .EQ. 0) THEN
               IF (IDL .EQ. 0) THEN
                  DO K = 1, NPP12
                     ETIJ(K,0,0,0) = D1
                  END DO
               ELSE IF (IDL .EQ. 1) THEN
                  DO K = 1,NPP12
                     PAK  = PA(K)/EXP12(K,3)
                     EXP1 = D2*EXP12(K,1)
                     ETIJ(K,0,0,0) = PAK*EXP1
                  END DO
               ELSE IF (IDL .GT. 1) THEN
                  DIDLM = D1 - DLEVEL
                  DO K = 1,NPP12
                     EXPH = EXP12(K,1)/EXP12(K,3)
                     PAK1 = PA(K)*ETIJM1(K,0,0,0)
                     PAK2 = DIDLM*EXP12(K,2)*ETIJM2(K,0,0,0)
                     ETIJ(K,0,0,0) = D2*EXPH*(PAK1 + PAK2)
                  END DO
               ENDIF
C
C           E(1,0)
C
            ELSE IF (I .EQ. 1) THEN
               IF (IDL. EQ. 0) THEN
                  DO K = 1, NPP12
                     ETIJ(K,0,1,0) = PA(K)
                     ETIJ(K,1,1,0) = HPI(K)
                  END DO
               ELSE 
                  DO K = 1,NPP12
                     PAK  = PA(K)
                     HPK  = HPI(K)
                     EXP2 = -DLEVEL*EXP12(K,2)
                     ETIJ(K,1,1,0) =    HPK*ETIJ(K,0,0,0)
                     ETIJ(K,0,1,0) =    PAK*ETIJ(K,0,0,0)
     &                             + EXP2*ETIJM1(K,0,0,0)
                  END DO
               ENDIF
C
C           E(2,0)
C
            ELSE IF (I .EQ. 2) THEN
               IF (IDL .EQ. 0) THEN
                  DO K = 1, NPP12
                     PAK = PA(K)
                     HPK = HPI(K)
                     ETIJ(K,0,2,0) = PAK*PAK + SIGN*HPK
                     ETIJ(K,1,2,0) = D2*PAK*HPK
                     ETIJ(K,2,2,0) = HPK*HPK
                  END DO
               ELSE
                  DO K = 1,NPP12
                     PAK  = PA(K)
                     HPK  = HPI(K)
                     EXP2 = -DLEVEL*EXP12(K,2)
                     ETIJ(K,0,2,0) =   SIGN*ETIJ(K,1,1,0)
     &                             +    PAK*ETIJ(K,0,1,0)
     &                             + EXP2*ETIJM1(K,0,1,0)
                     ETIJ(K,1,2,0) =    HPK*ETIJ(K,0,1,0)
     &                             +    PAK*ETIJ(K,1,1,0)
     &                             + EXP2*ETIJM1(K,1,1,0)
                     ETIJ(K,2,2,0) =    HPK*ETIJ(K,1,1,0)
                  END DO
               ENDIF
C
C            E(I,0)
C
            ELSE
               DO K = 1, NPP12
                  PAK = PA(K)
                  HPK = HPI(K)
                  ETIJ(K,  0,I,0) = PAK*ETIJ(K,  0,I-1,0)
     &                           + SIGN*ETIJ(K,  1,I-1,0)
                  ETIJ(K,I-1,I,0) = HPK*ETIJ(K,I-2,I-1,0)
     &                            + PAK*ETIJ(K,I-1,I-1,0)
                  ETIJ(K,  I,I,0) = HPK*ETIJ(K,I-1,I-1,0)
               END DO
               IF (IDL .GT. 0) THEN
                  DO K = 1,NPP12
                     EXP2 = -DLEVEL*EXP12(K,2)
                     ETIJ(K,  0,I,0) =        ETIJ(K,  0,  I,0)
     &                               + EXP2*ETIJM1(K,  0,I-1,0)
                     ETIJ(K,I-1,I,0) =        ETIJ(K,I-1,  I,0)
     &                               + EXP2*ETIJM1(K,I-1,I-1,0)
                  END DO
               END IF
               DO T = 1, I - 2
                  T1 = SIGN*(T + 1.0D0)
                  DO K = 1, NPP12
                     ETIJ(K,T,I,0) = HPI(K)*ETIJ(K,T-1,I-1,0)
     &                              + PA(K)*ETIJ(K,  T,I-1,0)
     &                                 + T1*ETIJ(K,T+1,I-1,0)
                  END DO
               END DO
               IF (IDL .GT. 0) THEN
                  DO T = 1,I-2
                     DO K = 1,NPP12
                        EXP2 = -DLEVEL*EXP12(K,2)
                        ETIJ(K,T,I,0) =        ETIJ(K,T,  I,0)
     &                                + EXP2*ETIJM1(K,T,I-1,0)
                     END DO
                  END DO
               END IF
            END IF
C
C           Run over J
C           ==========
C
            DO J = 1, MAXJ
               IJ = I + J
C
C              E(0,1)
C
               IF (IJ .EQ. 1) THEN
                  IF (IDL .EQ. 0) THEN
                     DO K = 1, NPP12
                        ETIJ(K,0,0,1) = PB(K)
                        ETIJ(K,1,0,1) = HPI(K)
                     END DO
                  ELSE
                     DO K = 1, NPP12
                        HPK  = HPI(K)
                        PBK  = PB(K)
                        EXP1 = DLEVEL*EXP12(K,1)
                        ETIJ(K,1,0,1) =    HPK*ETIJ(K,0,0,0)
                        ETIJ(K,0,0,1) =    PBK*ETIJ(K,0,0,0)
     &                                + EXP1*ETIJM1(K,0,0,0)
                     END DO
                  ENDIF
               ELSE IF (IJ .EQ. 2) THEN
C
C                 E(0,2)
C
                  IF (IDL .EQ. 0) THEN
                     IF (I .EQ. 0) THEN
                        DO K = 1, NPP12
                           PBK = PB(K)
                           HPK = HPI(K)
                           ETIJ(K,0,0,2) = PBK*PBK + SIGN*HPK
                           ETIJ(K,1,0,2) = D2*PBK*HPK
                           ETIJ(K,2,0,2) = HPK*HPK
                        END DO
C
C                    E(1,1)
C
                     ELSE
                        DO K = 1, NPP12
                           PAK = PA(K)
                           PBK = PB(K)
                           HPK = HPI(K)
                           ETIJ(K,0,1,1) = PAK*PBK + SIGN*HPK
                           ETIJ(K,1,1,1) = (PAK + PBK)*HPK
                           ETIJ(K,2,1,1) = HPK*HPK
                        END DO
                     END IF
                  ELSE
                     IF (I .EQ. 0) THEN
                        DO K = 1, NPP12
                           PBK  = PB(K)
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,0,0,2) =   SIGN*ETIJ(K,1,0,1)
     &                                   +    PBK*ETIJ(K,0,0,1)
     &                                   + EXP1*ETIJM1(K,0,0,1)
                           ETIJ(K,1,0,2) =    HPK*ETIJ(K,0,0,1)
     &                                   +    PBK*ETIJ(K,1,0,1)
     &                                   + EXP1*ETIJM1(K,1,0,1)
                           ETIJ(K,2,0,2) =    HPK*ETIJ(K,1,0,1)
                        END DO
                     ELSE
                        DO K = 1, NPP12
                           PBK  = PB(K)
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,0,1,1) =   SIGN*ETIJ(K,1,1,0)
     &                                   +    PBK*ETIJ(K,0,1,0)
     &                                   + EXP1*ETIJM1(K,0,1,0)
                           ETIJ(K,1,1,1) =    HPK*ETIJ(K,0,1,0)
     &                                   +    PBK*ETIJ(K,1,1,0)
     &                                   + EXP1*ETIJM1(K,1,1,0)
                           ETIJ(K,2,1,1) = HPK*ETIJ(K,1,1,0)
                        END DO
                     END IF
                  ENDIF
C
C              E(I,J)
C
               ELSE
                  DO K = 1, NPP12
                     PBK = PB(K)
                     HPK = HPI(K)
                     ETIJ(K,   0,I,J) = PBK*ETIJ(K,   0,I,J-1)
     &                               + SIGN*ETIJ(K,   1,I,J-1)
                     ETIJ(K,IJ-1,I,J) = HPK*ETIJ(K,IJ-2,I,J-1)
     &                                + PBK*ETIJ(K,IJ-1,I,J-1)
                     ETIJ(K,  IJ,I,J) = HPK*ETIJ(K,IJ-1,I,J-1)
                  END DO
                  IF (IDL .GT. 0) THEN
                     DO K = 1, NPP12
                        EXP1 = DLEVEL*EXP12(K,1)
                        ETIJ(K,   0,I,J) =        ETIJ(K,   0,I,  J)
     &                                   + EXP1*ETIJM1(K,   0,I,J-1)
                        ETIJ(K,IJ-1,I,J) =        ETIJ(K,IJ-1,I,  J)
     &                                   + EXP1*ETIJM1(K,IJ-1,I,J-1)
                     END DO
                  END IF
                  DO T = 1, IJ - 2
                     T1 = SIGN*(T + 1.0D0)
                     DO K = 1, NPP12
                        ETIJ(K,T,I,J)=HPI(K)*ETIJ(K,T-1,I,J-1)
     &                               + PB(K)*ETIJ(K,  T,I,J-1)
     &                                  + T1*ETIJ(K,T+1,I,J-1)
                     END DO
                  END DO
                  IF (IDL .GT. 0) THEN
                     DO T = 1, IJ - 2
                        DO K = 1, NPP12
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,T,I,J) =        ETIJ(K,T,I,  J)
     &                                   + EXP1*ETIJM1(K,T,I,J-1)
                        END DO
                     END DO
                  ENDIF
               END IF
            END DO
         END DO
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
         DO I = 0, MAXI
C
C           E(0,0)
C
            IF (I .EQ. 0) THEN
               IF (IDL .EQ. 0) THEN
                  DO K = 1, NPP12
                     ETIJ(K,0,0,0) = D1
                  END DO
               ELSE
                  IF (IBTAND(1,IDL) .EQ. 0) THEN
                     DIDLM = D1 - DLEVEL
                     DO K = 1, NPP12
                        EXPH = EXP12(K,1)/EXP12(K,3)
                        PAKH = DIDLM*EXP12(K,2)*ETIJM2(K,0,0,0)
                        ETIJ(K,0,0,0) = D2*EXPH*PAKH
                     END DO
                  ENDIF
               ENDIF
C
C           E(1,0)
C
            ELSE IF (I .EQ. 1) THEN
               IF (IDL .EQ. 0) THEN
                  DO K = 1, NPP12
                     ETIJ(K,1,1,0) = HPI(K)
                  END DO
               ELSE IF (IDL .EQ. 1) THEN
                  DO K = 1, NPP12
                     ETIJ(K,0,1,0) = -DLEVEL*EXP12(K,2)
                  END DO
               ELSE IF ((IDL .GT. 1) .AND. (IBTAND(1,IDL) .EQ. 0)) THEN
                  DO K = 1, NPP12
                     ETIJ(K,1,1,0) = HPI(K)*ETIJ(K,0,0,0)
                  END DO
               ELSE
                  DO K = 1, NPP12
                     ETIJ(K,0,1,0) = -DLEVEL*EXP12(K,2)*ETIJM1(K,0,0,0)
                  END DO
               ENDIF
C
C           E(2,0)
C
            ELSE IF (I .EQ. 2) THEN
               IF (IDL .EQ. 0) THEN
                  DO K = 1, NPP12
                     HPK = HPI(K)
                     ETIJ(K,0,2,0) = SIGN*HPK
                     ETIJ(K,2,2,0) = HPK*HPK
                  END DO
               ELSE IF (IDL .EQ. 1) THEN
                  DO K = 1, NPP12
                     HPK  = HPI(K)
                     EXP2 = -D2*EXP12(K,2)
                     ETIJ(K,1,2,0) = HPK*EXP2
                  END DO
               ELSE IF ((IDL .GT. 1) .AND. (IBTAND(1,IDL) .EQ. 0)) THEN
                  DO K = 1, NPP12
                     EXP2 = -DLEVEL*EXP12(K,2)
                     ETIJ(K,0,2,0) = EXP2*ETIJM1(K,0,1,0)
     &                             +   SIGN*ETIJ(K,1,1,0)
                     ETIJ(K,2,2,0) = HPI(K)*ETIJ(K,1,1,0)
                  END DO
               ELSE 
                  DO K = 1, NPP12
                     HPK  = HPI(K)
                     EXP2 = -DLEVEL*EXP12(K,2)
                     ETIJ(K,1,2,0) =    HPK*ETIJ(K,0,1,0)
     &                             + EXP2*ETIJM1(K,1,1,0)
                  END DO
               ENDIF
C
C            E(I,0)
C
            ELSE
               IF ((IDL .EQ. 0) .OR.
     &             ((IDL .GT. 1) .AND. (IBTAND(1,IDL) .EQ. 0))) THEN
                  IF (IBTAND(1,I) .EQ. 0) THEN
                     DO K = 1, NPP12
                        ETIJ(K,0,I,0) =   SIGN*ETIJ(K,  1,I-1,0)
                        ETIJ(K,I,I,0) = HPI(K)*ETIJ(K,I-1,I-1,0)
                     END DO
                     IF (IDL .GT. 1) THEN
                        DO K = 1, NPP12
                           EXP2 = -DLEVEL*EXP12(K,2)
                           ETIJ(K,0,I,0) =        ETIJ(K,0,  I,0)
     &                                   + EXP2*ETIJM1(K,0,I-1,0)
                        END DO
                     ENDIF
                  ELSE
                     DO K = 1, NPP12
                        ETIJ(K,I,I,0) = HPI(K)*ETIJ(K,I-1,I-1,0)
                     END DO
                  END IF
                  DO T = 2 - IBTAND(1,I), I - 2, 2
                     T1 = SIGN*(T + 1.0D0)
                     DO K = 1, NPP12
                        ETIJ(K,T,I,0) = HPI(K)*ETIJ(K,T-1,I-1,0)
     &                                    + T1*ETIJ(K,T+1,I-1,0)
                     END DO
                     IF (IDL .GT. 1) THEN
                        DO K = 1, NPP12
                           EXP2 = -DLEVEL*EXP12(K,2)
                           ETIJ(K,T,I,0) =        ETIJ(K,T,I,0)
     &                                   + EXP2*ETIJM1(K,T,I,0)
                        END DO
                     ENDIF
                  END DO
               ELSE 
                  IF (IBTAND(1,I) .EQ. 0) THEN
                     DO K = 1,NPP12
                        HPK  = HPI(K)
                        EXP2 = -DLEVEL*EXP12(K,2)
                        ETIJ(K,I-1,I,0) =    HPK*ETIJ(K,I-2,I-1,0)
     &                                  + EXP2*ETIJM1(K,I-1,I-1,0)
                     END DO
                  ELSE
                     DO K = 1,NPP12
                        HPK  = HPI(K)
                        EXP2 = -DLEVEL*EXP12(K,2)
                        ETIJ(K,I-1,I,0) =    HPK*ETIJ(K,I-2,I-1,0)
     &                                  + EXP2*ETIJM1(K,I-1,I-1,0)
                        ETIJ(K,  0,I,0) = EXP2*ETIJM1(K,  0,I-1,0)
     &                                  +   SIGN*ETIJ(K,  1,I-1,0)
                     END DO
                  ENDIF
                  DO T = 1 + IBTAND(1,I), I - 3, 2
                     T1   = SIGN*(T + 1.0D0)
                     DO K = 1, NPP12
                        HPK  = HPI(K)
                        EXP2 = -DLEVEL*EXP12(K,2)
                        ETIJ(K,T,I,0) =    HPK*ETIJ(K,T-1,I-1,0)
     &                                + EXP2*ETIJM1(K,  T,I-1,0)
     &                                +     T1*ETIJ(K,T+1,I-1,0)
                     END DO
                  END DO
               ENDIF
            END IF
C
C           Run over J
C           ==========
C
            DO J = 1, MAXJ
               IJ = I + J
C
C              E(0,1)
C
               IF (IJ .EQ. 1) THEN
                  IF (IDL .EQ. 0) THEN
                     DO K = 1, NPP12
                        ETIJ(K,1,0,1) = HPI(K)
                     END DO
                  ELSE IF (IDL .EQ. 1) THEN
                     DO K = 1, NPP12
                        ETIJ(K,0,0,1) = EXP12(K,1)
                     END DO
                  ELSE IF ((IDL .GT. 1) .AND.
     &                     (IBTAND(1,IDL) .EQ. 0)) THEN
                     DO K = 1, NPP12
                        ETIJ(K,1,0,1) = HPI(K)*ETIJ(K,0,0,0)
                     END DO
                  ELSE
                     DO K = 1, NPP12
                        ETIJ(K,0,0,1) = DLEVEL*EXP12(K,1)*ETIJ(K,0,0,0)
                     END DO
                  ENDIF
C
C              E(1,1) AND E(0,2)
C
               ELSE IF (IJ .EQ. 2) THEN
                  IF (IDL .EQ. 0) THEN
                     DO K = 1, NPP12
                        HPK = HPI(K)
                        ETIJ(K,0,I,J) = SIGN*HPK
                        ETIJ(K,2,I,J) = HPK*HPK
                     END DO
                  ELSE IF (IDL .EQ. 1) THEN
                     IF (J .EQ. 1) THEN
                        DO K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = EXP12(K,1)
                           EXP2 = -EXP12(K,2)
                           ETIJ(K,1,1,1) = HPK*(EXP1 + EXP2)
                        END DO
                     ELSE
                        DO K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = D2*EXP12(K,1)
                           ETIJ(K,1,0,2) = HPK*EXP1
                        END DO
                     ENDIF
                  ELSE IF ((IDL .GT. 1) .AND.
     &                     (IBTAND(1,IDL) .EQ. 0)) THEN
                     IF (J .EQ. 1) THEN
                        DO K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,2,1,1) =    HPK*ETIJ(K,1,1,0)
                           ETIJ(K,0,1,1) = EXP1*ETIJM1(K,0,1,0)
     &                                   +   SIGN*ETIJ(K,1,1,0)
                        END DO
                     ELSE
                        DO K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,2,0,2) =    HPK*ETIJ(K,1,0,1)
                           ETIJ(K,0,0,2) = EXP1*ETIJM1(K,0,0,1)
     &                                   +   SIGN*ETIJ(K,1,0,1)
                        END DO
                     ENDIF
                  ELSE
                     IF (J .EQ. 1) THEN
                        DO K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,1,1,1) =    HPK*ETIJ(K,0,1,0)
     &                                   + EXP1*ETIJM1(K,1,1,0)
                        END DO
                     ELSE
                        DO K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,1,0,2) =    HPK*ETIJ(K,0,0,1)
     &                                   + EXP1*ETIJM1(K,1,0,1)
                        END DO
                     ENDIF
                  ENDIF
C
C              E(I,J)
C
               ELSE
                  IF ((IDL .EQ. 0) .OR.
     &               ((IDL .GT. 1) .AND. (IBTAND(1,IDL) .EQ. 0))) THEN
                     IF (IBTAND(1,IJ) .EQ. 0) THEN
                        DO K = 1, NPP12
                           ETIJ(K, 0,I,J) =   SIGN*ETIJ(K,   1,I,J-1)
                           ETIJ(K,IJ,I,J) = HPI(K)*ETIJ(K,IJ-1,I,J-1)
                        END DO
                        IF (IDL .GT. 1) THEN
                           DO K = 1, NPP12
                              EXP1 = DLEVEL*EXP12(K,1)
                              ETIJ(K,0,I,J) =        ETIJ(K,0,I,  J)
     &                                      + EXP1*ETIJM1(K,0,I,J-1)
                           END DO
                        ENDIF
                     ELSE
                        DO K = 1, NPP12
                           ETIJ(K,IJ,I,J) = HPI(K)*ETIJ(K,IJ-1,I,J-1)
                        END DO
                     END IF
                     DO T = 2 - IBTAND(1,IJ), IJ - 2, 2
                        T1 = SIGN*(T + 1.0D0)
                        DO K = 1, NPP12
                           ETIJ(K,T,I,J) = HPI(K)*ETIJ(K,T-1,I,J-1)
     &                                       + T1*ETIJ(K,T+1,I,J-1)
                        END DO
                        IF (IDL .GT. 1) THEN
                           DO K = 1, NPP12
                              EXP1 = DLEVEL*EXP12(K,1)
                              ETIJ(K,T,I,J) =        ETIJ(K,T,I,  J)
     &                                      + EXP1*ETIJM1(K,T,I,J-1)
                           END DO
                        ENDIF
                     END DO
                  ELSE
                     IF (IBTAND(1,IJ) .EQ. 0) THEN
                        DO K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,IJ-1,I,J) =    HPK*ETIJ(K,IJ-2,I,J-1)
     &                                      + EXP1*ETIJM1(K,IJ-1,I,J-1)
                        END DO
                     ELSE
                        DO K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,IJ-1,I,J) =    HPK*ETIJ(K,IJ-2,I,J-1)
     &                                      + EXP1*ETIJM1(K,IJ-1,I,J-1)
                           ETIJ(K,   0,I,J) = EXP1*ETIJM1(K,   0,I,J-1)
     &                                      +   SIGN*ETIJ(K,   1,I,J-1)
                        END DO
                     END IF
                     DO T = 1 + IBTAND(1,IJ), IJ - 3, 2
                        T1 = SIGN*(T + 1.0D0)
                        DO K = 1, NPP12
                           HPK  = HPI(K)
                           EXP1 = DLEVEL*EXP12(K,1)
                           ETIJ(K,T,I,J) =    HPK*ETIJ(K,T-1,I,J-1)
     &                                   + EXP1*ETIJM1(K,  T,I,J-1)
     &                                   +     T1*ETIJ(K,T+1,I,J-1)
                        END DO
                     END DO
                  ENDIF
               END IF
            END DO
         END DO
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
            IF (X .EQ. 1) WORD = 'EX00'
            IF (X .EQ. 2) WORD = 'EY00'
            IF (X .EQ. 3) WORD = 'EZ00'
            DO I = 0, MAXI
            DO J = 0, MAXJ
            DO T = 0, I + J
            IF (IBTAND(I + J - T + IDL,IAB0) .EQ. 0) THEN
               WRITE (LUPRI,'(/,2X,A4,A1,I1,A1,I1,A1,I1,A1,/)')
     &              WORD, '(', I, ',', J, ', ',  T, ')'
               WRITE (LUPRI,'(1X,6F12.8)') (ETIJ(K,T,I,J),K=1,NPP12)
            END IF
            END DO
            END DO
            END DO
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
      DO K = 1, NPP12
         SCALE(K) = ((EXP12(K,1)/EXP12(K,3))**ISCAL1)
     &            * ((EXP12(K,2)/EXP12(K,3))**ISCAL2)
      END DO
C
      DO J = 0, MAX2
      DO I = 0, MAX1
      DO T = 0, I + J
      DO K = 1, NPP12
         ETIJ(K,T,I,J,1) = SCALE(K)*ETIJ(K,T,I,J,1)
      END DO
      END DO
      END DO
      END DO
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
            DO X = 1, 3
               IF (X .EQ. 1) WORD = 'EX00'
               IF (X .EQ. 2) WORD = 'EY00'
               IF (X .EQ. 3) WORD = 'EZ00'
               DO I = 0, MAX1
               DO J = 0, MAX2
               DO T = 0, I + J
                  IF (IBTAND(I + J - T,IAB0(X)) .EQ. 0) THEN
                    WRITE (LUPRI,'(/,2X,A4,A1,I1,A1,I1,A1,I1,A1,/)')
     &                   WORD, '(', I, ',', J, ', ',  T, ')'
                    WRITE (LUPRI,'(1X,6F12.8)')
     &                   (ETIJ(K,T,I,J,X),K=1,NPP12)
                  END IF
               END DO
               END DO
               END DO
            END DO
         END IF
      END IF
C
      RETURN
      END
