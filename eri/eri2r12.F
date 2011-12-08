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
C FILE: eri/eri2r12.F
C
C  /* Deck eriifc */
      SUBROUTINE ERIIFC(NNTUV,NNPPX,NIPQ0X,NIPQ0Y,NIPQ0Z,NJMAX)
C
C     Written by Wim Klopper (University of Karlsruhe, 19 November 2002).
C
#include <implicit.h>
#include <priunit.h>
#include <ericom.h>
      NPPX = NNPPX
      NTUV = NNTUV
      IPQ0(1) = NIPQ0X
      IPQ0(2) = NIPQ0Y
      IPQ0(3) = NIPQ0Z
      JMAX = NJMAX
      RETURN
      END
C  /* Deck eriher */
      SUBROUTINE ERIHER(HRINT1,HRINT2,RJ000,COORPQ,INDHER,IODDHR,FACINT,
     &                  ALPHA,HARGE,HARGF,HALPH,HBETA,HPRE1,HPRE2,
     &                  HPRE3,RO000,HRIN1,HRIN2,HRIN3,
     &                  HEXPP,HEXPQ,WORK,LWORK,IPRINT)
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      INTEGER T, U, V, TUV
      DIMENSION HARGE(*), HARGF(*), HALPH(*), HBETA(*), HPRE1(*),
     &          HPRE2(*), HPRE3(*), RO000(*), HRIN1(*),
     &          HRIN2(*), HRIN3(*), HEXPP(*), HEXPQ(*)
      DIMENSION WORK(LWORK), HRINT1(NPPX,NTUV), HRINT2(NPPX,NTUV), 
     &          RJ000(*), COORPQ(NPPX,3), ALPHA(NPPX),
     &          INDHER(*), IODDHR(*), FACINT(NPPX)
#include <ericom.h>
#include <r12int.h>
      IF (R12EIN) THEN
        IF (SLATER) THEN
C
C        exp(r12) Hermite integrals
C        ==============================
C
         IF (INTGAC .EQ. 4 .OR. INTGAC .EQ. 6) THEN
           DO I = 1, NTOGAM
             GAMAC = GAMAX(I)
             DO J = 1, NTOGAM
              GAMAD = GAMAX(J)
              GAMFX = GAMAB(I) * GAMAB(J)
              INTGAD = 0
              CALL COEFF(HARGE,HARGF,HALPH,HBETA,HPRE1,HPRE2,HPRE3,
     &                   FACINT,HEXPP,HEXPQ,IPRINT)
              CALL ERIFGHTUV(HRIN1,HRIN2,HRIN3,COORPQ(1,1),COORPQ(1,2),
     &                       COORPQ(1,3),INDHER,HARGE,IPRINT)
              CALL WKER12(HRINT1,HRIN1,HRIN2,HRIN3,HPRE1,HPRE2,HPRE3,
     &                    GAMFX,I,J,IPRINT)
             END DO
           END DO
         ELSE
           DO I = 1, NTOGAM
             GAMAC = GAMAX(I)
             GAMAD = GAMAX(I)
             GAMFX = GAMAB(I)
             IF (INTGAC .EQ. 2) THEN
                INTGAD = 1
                CALL COEFF(HARGE,HARGF,HALPH,HBETA,HPRE1,HPRE2,HPRE3,
     &                     FACINT,HEXPP,HEXPQ,IPRINT)
                CALL ERIGOM(RO000,RJ000,HARGE,HARGF,COORPQ,
     &                      WORK,LWORK,IPRINT)
                CALL ERIOMTUV(HRIN1,HRIN2,RO000,RJ000,HARGE,HARGF,
     &                        COORPQ,INDHER,WORK,LWORK,IPRINT)
             ELSE
                INTGAD = 0
                CALL COEFF(HARGE,HARGF,HALPH,HBETA,HPRE1,HPRE2,HPRE3,
     &                     FACINT,HEXPP,HEXPQ,IPRINT)
                CALL ERIFGHTUV(HRIN1,HRIN2,HRIN3,COORPQ(1,1),
     &                         COORPQ(1,2),COORPQ(1,3),INDHER,HARGE,
     &                         IPRINT)
             END IF
             CALL WKER12(HRINT1,HRIN1,HRIN2,HRIN3,HPRE1,HPRE2,HPRE3,
     &                   GAMFX,I,1,IPRINT)
           END DO
         END IF
         IF (INTGAC .EQ. 5) CALL DCOPY(NPPX*NTUV,HRINT1,1,HRINT2,1)
        ELSE
C
C        r12*exp(r12) Hermite integrals
C        ==============================
C
         IF (INTGAC .EQ. 4 .OR. INTGAC .EQ. 6) THEN
           DO I = 1, NTOGAM
             GAMAC = GAMAX(I)
             DO J = 1, NTOGAM
              GAMAD = GAMAX(J)
              GAMFX = GAMAB(I) * GAMAB(J)
              INTGAD = 2 * MOD(INTGAC,2)
              CALL COEFF(HARGE,HARGF,HALPH,HBETA,HPRE1,HPRE2,HPRE3,
     &                   FACINT,HEXPP,HEXPQ,IPRINT)
              CALL ERIFGHTUV(HRIN1,HRIN2,HRIN3,COORPQ(1,1),COORPQ(1,2),
     &                       COORPQ(1,3),INDHER,HARGE,IPRINT)
              CALL WKER12(HRINT1,HRIN1,HRIN2,HRIN3,HPRE1,HPRE2,HPRE3,
     &                    GAMFX,I,J,IPRINT)
             END DO
           END DO
         ELSE
           DO I = 1, NTOGAM
             GAMAC = GAMAX(I)
             GAMAD = GAMAX(I)
             GAMFX = GAMAB(I)
             INTGAD = 2 * MOD(INTGAC,2)
             CALL COEFF(HARGE,HARGF,HALPH,HBETA,HPRE1,HPRE2,HPRE3,
     &                  FACINT,HEXPP,HEXPQ,IPRINT)
             CALL ERIFGHTUV(HRIN1,HRIN2,HRIN3,COORPQ(1,1),COORPQ(1,2),
     &                      COORPQ(1,3),INDHER,HARGE,IPRINT)
             IF (INTGAC .EQ. 3 .OR. INTGAC .EQ. 5) THEN
                INTGAD = 1
                CALL COEFF(HARGE,HARGF,HALPH,HBETA,HPRE1,HPRE2,HPRE3,
     &                     FACINT,HEXPP,HEXPQ,IPRINT)
                CALL ERIGOM(RO000,RJ000,HARGE,HARGF,COORPQ,
     &                      WORK,LWORK,IPRINT)
                CALL ERIOMTUV(HRIN2,HRIN3,RO000,RJ000,HARGE,HARGF,
     &                        COORPQ,INDHER,WORK,LWORK,IPRINT)
             END IF
             CALL WKER12(HRINT1,HRIN1,HRIN2,HRIN3,HPRE1,HPRE2,HPRE3,
     &                   GAMFX,I,1,IPRINT)
           END DO
         END IF
         IF (INTGAC .EQ. 5) CALL DCOPY(NPPX*NTUV,HRINT1,1,HRINT2,1)
        END IF
      ELSE
C
C        1/r12 Hermite integrals
C        =======================
C
         CALL ERITUV(HRINT1,RJ000,COORPQ,INDHER,IODDHR,
     &                     WORK,LWORK,IPRINT)
         IF (R12INT .OR. U12INT) THEN
C
C           r12 Hermite integrals
C           =====================
C
            IF (NPPX .GT. LWORK) 
     &      CALL STOPIT('ERIHER','WKEQ52',NPPX,LWORK)
            CALL WKEQ52(HRINT2,HRINT1,ALPHA,WORK,COORPQ(1,1),
     &                  COORPQ(1,2),COORPQ(1,3),INDHER,JMAX,NPPX,
     &                  NTUV,IPQ0(1),IPQ0(2),IPQ0(3),IODDHR,IPRINT)
            IF (.NOT. U12INT)
     &      CALL DCOPY(NPPX*NTUV,HRINT2,1,HRINT1,1)
         END IF
      END IF
C
C     Print section
C     =============
C
      IF (IPRINT .GE. 10) THEN
         CALL AROUND('Output from ERIHER')
         DO I=1,NTUV
            WRITE(LUPRI,'(A,I5)') ' TUV =',I
            DO J = 1, NPPX
               WRITE(LUPRI,'(I5,2F14.9)') J,HRINT1(J,I),HRINT2(J,I)
            END DO
         END DO
      END IF
      RETURN
      END
C  /* Deck wkeq52 */
      SUBROUTINE WKEQ52(Q,R,W,A,PQX,PQY,PQZ,INDHER,JMAX,NUABCD,
     &                  NTUV,IPQ0X,IPQ0Y,IPQ0Z,IODDHR,IPRINT)
C
C     Q(t,u,v) integrals are computed according to:
C
C     W. Klopper, R. Roehse, Theor. Chim. Acta (1992) 83:441, Eq. (52).
C     -----------------------------------------------------------------
C
C     Q(0,0,0) is obtained from Eq. (53).
C
C     The array R contains the integrals R(t,u,v) over 1/r12 on
C     input, and the array Q contains the integrals Q(t,u,v) over
C     r12 on output. The array A is used as scratch.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h> 
      PARAMETER (DM1 = -1.0D0, D2 = 2.0D0)
      INTEGER T, U, V, TUV
      LOGICAL PQXGT0, PQYGT0, PQZGT0
      DIMENSION Q(NUABCD,NTUV), R(NUABCD,NTUV), A(NUABCD), 
     &          W(NUABCD), IODDHR(NTUV), INDHER(0:JTOP,0:JTOP,0:JTOP), 
     &          PQX(NUABCD), PQY(NUABCD), PQZ(NUABCD)
#include <hertop.h>
C
      IPQX = IPQ0X + 1
      IPQY = IPQ0Y + 1
      IPQZ = IPQ0Z + 1
C
      PQXGT0 = IPQX .EQ. 1
      PQYGT0 = IPQY .EQ. 1
      PQZGT0 = IPQZ .EQ. 1
C
      CALL DZERO(Q,NUABCD*NTUV)
C
C     Q(0,0,0)
C     ========
C
      DO I = 1, NUABCD
         A(I) = DM1/W(I)
         Q(I,1) = D2*A(I)*R(I,1)
      END DO
      IF (PQXGT0) THEN
         DO I = 1, NUABCD
            Q(I,1) = Q(I,1) + PQX(I) * (A(I)*R(I,2) + PQX(I)*R(I,1))
         END DO
      END IF
      IF (PQYGT0) THEN
         DO I = 1, NUABCD
            Q(I,1) = Q(I,1) + PQY(I) * (A(I)*R(I,3) + PQY(I)*R(I,1))
         END DO
      END IF
      IF (PQZGT0) THEN
         DO I = 1, NUABCD
            Q(I,1) = Q(I,1) + PQZ(I) * (A(I)*R(I,4) + PQZ(I)*R(I,1))
         END DO
      END IF
C
C     Q(T,0,0)
C     ========
C
      IF (PQXGT0) THEN
         DO I = 1, NUABCD
            Q(I,2) = A(I)*R(I,2) + PQX(I)*R(I,1)
         END DO
         DO T = 2, JMAX
            TMIN1 = T - 1.0D0
            TUV   = INDHER(T  ,0,0)
            M1T   = INDHER(T-1,0,0)
            M2T   = INDHER(T-2,0,0) 
            DO I = 1, NUABCD
               Q(I,TUV) = A(I)*R(I,TUV)
     &                  + PQX(I)*R(I,M1T) + TMIN1*R(I,M2T)
            END DO
         END DO
      ELSE
         DO T = 2, JMAX, 2
            TMIN1 = T - 1.0D0
            TUV   = INDHER(T  ,0,0)
            M2T   = INDHER(T-2,0,0)
            DO I = 1, NUABCD
               Q(I,TUV) = A(I)*R(I,TUV) + TMIN1*R(I,M2T)
            END DO
         END DO
      END IF
C
C     Q(T,U,0)
C     ========
C
      IF (PQYGT0) THEN
         DO T = 0, JMAX - 1, IPQX
            TUV = INDHER(T,1,0)
            M1U = INDHER(T,0,0)
            DO I = 1, NUABCD
               Q(I,TUV) = A(I)*R(I,TUV) + PQY(I)*R(I,M1U)
            END DO
         END DO
         DO U = 2, JMAX
            UMIN1  = U - 1.0D0
            DO T = 0, JMAX - U, IPQX
               TUV = INDHER(T,U  ,0)
               M1U = INDHER(T,U-1,0)
               M2U = INDHER(T,U-2,0)
               DO I = 1, NUABCD
                  Q(I,TUV) = A(I)*R(I,TUV)
     &                     + PQY(I)*R(I,M1U) + UMIN1*R(I,M2U)
               END DO
            END DO
         END DO
      ELSE
         DO U = 2, JMAX, 2
            UMIN1  = U - 1.0D0
            DO T = 0, JMAX - U, IPQX
               TUV = INDHER(T,U  ,0)
               M2U = INDHER(T,U-2,0)
               DO I = 1, NUABCD
                  Q(I,TUV) = A(I)*R(I,TUV) + UMIN1*R(I,M2U)
               END DO
            END DO
         END DO
      END IF
C
C     Q(T,U,V)
C     ========
C
      IF (PQZGT0) THEN
         IUMAX  = JMAX - 1
         DO U = 0, IUMAX, IPQY
            DO T = 0, IUMAX - U, IPQX
               TUV = INDHER(T,U,1)
               M1V = INDHER(T,U,0)
               DO I = 1, NUABCD
                  Q(I,TUV) = A(I)*R(I,TUV) + PQZ(I)*R(I,M1V)
               END DO
            END DO
         END DO
         DO V = 2, JMAX
            VMIN1  = V - 1.0D0
            IUMAX  = JMAX - V
            DO U = 0, IUMAX, IPQY
               DO T = 0, IUMAX - U, IPQX
                  TUV = INDHER(T,U,V  )
                  M1V = INDHER(T,U,V-1)
                  M2V = INDHER(T,U,V-2)
                  DO I = 1, NUABCD 
                     Q(I,TUV) = A(I)*R(I,TUV)
     &                        + PQZ(I)*R(I,M1V) + VMIN1*R(I,M2V)
                  END DO
               END DO
            END DO
         END DO
      ELSE
         DO V = 2, JMAX, 2
            VMIN1  = V - 1.0D0
            IUMAX  = JMAX - V
            DO U = 0, IUMAX, IPQY
               DO T = 0, IUMAX - U, IPQX
                  TUV = INDHER(T,U,V  )
                  M2V = INDHER(T,U,V-2)
                  DO I = 1, NUABCD
                     Q(I,TUV) = A(I)*R(I,TUV) + VMIN1*R(I,M2V)
                  END DO
               END DO
            END DO
         END DO
      END IF
C
C     Print section
C     ============= 
C
      IF (IPRINT .GE. 10) THEN
         CALL TITLER('Output from WKEQ52','*',103)
         WRITE (LUPRI,'(2X,A,I5)') 'JMAX  ', JMAX
         WRITE (LUPRI,'(2X,A,I5)') 'NUABCD', NUABCD
         WRITE (LUPRI,'(2X,A,I5)') 'NTUV  ', NTUV
         IF (IPRINT .GE. 20) THEN
            CALL HEADER('Hermite integrals Q(t,u,v)',1)
            DO J = 0, JMAX
              DO T = J, 0, -1
                DO U = J - T, 0, -1
                  V = J - T - U
                  TUV = INDHER(T,U,V)
                  IF (IODDHR(TUV) .EQ. 0) THEN
                     WRITE (LUPRI,'(2X,3(A,I1),A,2X,5F12.8/,
     &                                                  (12X,5F12.8))')
     &              'Q(',T,',',U,',',V,')', (Q(I,TUV),I=1,NUABCD)
                     WRITE (LUPRI,'()')
                  END IF
                END DO
              END DO
            END DO
         END IF
      END IF
      RETURN   
      END
C  /* Deck wkeq30 */
      SUBROUTINE WKEQ30(N,C0,C1,CM2,JMAX1,JMAX2,
     &                  NUCAB,MAX1,MAX2,SIGN,
     &                  I120,DIFPA,DIFPB,TEXPA,TEXPB,
     &                  AOVERP,BOVERP,HEXPPI,WORD,IPRINT)
C
C     E(i,j,t;n) coefficients are computed for n=1 according to:
C    
C     W. Klopper, R. Roehse, Theor. Chim. Acta (1992) 83:441, Eq. (30).
C     -----------------------------------------------------------------
C     cf. T. Helgaker, P.R. Taylor, Theor. Chim. Acta (1992) 83:177.
C     --------------------------------------------------------------
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
c
c     Modified by Cristian Villani (Karlsruhe, Feb 2005) for n>1
C
#include <implicit.h>
#include <priunit.h> 
      INTEGER T             
      CHARACTER WORD*4
      DIMENSION DIFPA(NUCAB), DIFPB(NUCAB), HEXPPI(NUCAB)
      DIMENSION TEXPA(NUCAB), TEXPB(NUCAB)
      DIMENSION AOVERP(NUCAB), BOVERP(NUCAB)
      DIMENSION C0(NUCAB,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2)
      DIMENSION C1(NUCAB,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2)
c     if n>1, the E(k,0,0,0,n-2) coefficients are needed, for all k
      DIMENSION CM2(NUCAB)
c     if n<2, this vector is not used
c     Note that cm2(K)=1, if n=2, since the coefficients are
c     normalized to Kab
c
c     C. Villani Uni-Ka, Feb 2005
c
#include <ibtfun.h>

      DO 100 K = 1, NUCAB
         AOVERP(K) =   SIGN * TEXPA(K) * HEXPPI(K)
         BOVERP(K) = - SIGN * TEXPB(K) * HEXPPI(K)
  100 CONTINUE
C
      GO TO (1) IBTAND(I120,1)
C
      DO 1000 I = 0, MAX1
         IF (I .EQ. 0) THEN
C
C           Compute starting values from Eq.(31).
C           Note that (-2ab/p)*Rx = 2a*(Px-Ax).
C
            DO 200 K = 1, NUCAB
               C1(K,0,0,0) = TEXPA(K)*DIFPA(K)*C0(K,0,0,0)
               IF(N.GT.1) C1(K,0,0,0) = C1(K,0,0,0) + 
     .            (N-1)*TEXPA(K)* BOVERP(K)*CM2(K)
  200       CONTINUE
         ELSE IF (I .EQ. 1) THEN
            DO 300 K = 1, NUCAB
               C1(K,0,1,0) =  DIFPA(K)*C1(K,0,0,0)
     &                     + N*BOVERP(K)*C0(K,0,0,0)
               C1(K,1,1,0) = HEXPPI(K)*C1(K,0,0,0)
  300       CONTINUE
         ELSE IF (I .GE. 2) THEN
            DO 400 K = 1, NUCAB
               C1(K,0  ,I,0) =  DIFPA(K)*C1(K,0,I-1,0)
     &                       + N*BOVERP(K)*C0(K,0,I-1,0)
     &                       +      SIGN*C1(K,1,I-1,0)
               C1(K,I-1,I,0) = HEXPPI(K)*C1(K,I-2,I-1,0)
     &                       +  DIFPA(K)*C1(K,I-1,I-1,0)
     &                       + N*BOVERP(K)*C0(K,I-1,I-1,0)
               C1(K,I  ,I,0) = HEXPPI(K)*C1(K,I-1,I-1,0)
  400       CONTINUE
            IF (I .GE. 3) THEN
               DO 510 T = 1, I - 2
               TPLUS1 = SIGN*(T + 1.0D0)
               DO 520 K = 1, NUCAB
                  C1(K,T,I,0) = HEXPPI(K)*C1(K,T-1,I-1,0)
     &                        +  DIFPA(K)*C1(K,T  ,I-1,0)
     &                        + N*BOVERP(K)*C0(K,T  ,I-1,0)
     &                        +    TPLUS1*C1(K,T+1,I-1,0)
  520             CONTINUE
  510          CONTINUE
            END IF
         END IF
         IF (MAX2 .LE. 0) GO TO 1000
         DO 600 J = 1, MAX2
            IJ = I + J
            IF (IJ .EQ. 1) THEN
               DO 700 K = 1, NUCAB
                  C1(K,0,I,J) =  DIFPB(K)*C1(K,0,I,J-1)
     &                        + N*AOVERP(K)*C0(K,0,I,J-1)
                  C1(K,1,I,J) = HEXPPI(K)*C1(K,0,I,J-1)
  700          CONTINUE
            ELSE IF (IJ .GE. 2) THEN
               DO 800 K = 1, NUCAB
                  C1(K,0   ,I,J) =  DIFPB(K)*C1(K,0   ,I,J-1)
     &                           + N*AOVERP(K)*C0(K,0   ,I,J-1)
     &                           +      SIGN*C1(K,1   ,I,J-1)
                  C1(K,IJ-1,I,J) = HEXPPI(K)*C1(K,IJ-2,I,J-1)
     &                           +  DIFPB(K)*C1(K,IJ-1,I,J-1)
     &                           + N*AOVERP(K)*C0(K,IJ-1,I,J-1)
                  C1(K,IJ  ,I,J) = HEXPPI(K)*C1(K,IJ-1,I,J-1)
  800          CONTINUE
               IF (IJ .GE. 3) THEN
                  DO 910 T = 1, IJ - 2
                     TPLUS1 = SIGN*(T + 1.0D0)
                     DO 920 K = 1, NUCAB
                     C1(K,T,I,J) = HEXPPI(K)*C1(K,T-1,I,J-1)
     &                           +  DIFPB(K)*C1(K,T  ,I,J-1)
     &                           + N*AOVERP(K)*C0(K,T  ,I,J-1)
     &                           +    TPLUS1*C1(K,T+1,I,J-1)
  920                CONTINUE
  910             CONTINUE
               END IF
            END IF
  600    CONTINUE
 1000 CONTINUE
      IF (IPRINT .LT. 10) RETURN
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      WRITE (LUPRI, 1110)
      WRITE (LUPRI, 1010) MAX1, MAX2
      WRITE (LUPRI, 1030) I120
      IF (IPRINT .LT. 20) RETURN
      DO 3000 I = 0, MAX1
         DO 3100 J = 0, MAX2
            DO 3200 T = 0, I + J
               WRITE (LUPRI, 1100) WORD, I, J, T
               WRITE (LUPRI, 1130) (C1(K,T,I,J), K = 1, NUCAB)
 3200       CONTINUE
 3100    CONTINUE
 3000 CONTINUE
      RETURN
C
    1 CONTINUE
C
      DO 1001 I = 0, MAX1
         IF(I.EQ.0.AND.N.GT.1) THEN
           DO 333 K = 1, NUCAB
             C1(K,0,0,0) = (N-1)*TEXPA(K)*BOVERP(K)*CM2(K)
  333      CONTINUE
         ELSE IF (I .EQ. 1) THEN
            IF (IBTAND(N,1) .EQ. 1) THEN
               DO 301 K = 1, NUCAB
                  C1(K,0,1,0) = N*BOVERP(K)*C0(K,0,0,0)
  301          CONTINUE
            ELSE
c              case n even, C. Villani Uni-Ka, Feb 2005
               DO K = 1, NUCAB
                  C1(K,1,1,0) = HEXPPI(K)*C1(K,0,0,0)
               END DO
            END IF
         ELSE IF (I .GE. 2) THEN
            IF (IBTAND(N,1) .EQ. 1) THEN
               IF (IBTAND(I,1) .EQ. 1) THEN
                 DO 401 K = 1, NUCAB
                   C1(K,0,I,0) = N*BOVERP(K)*C0(K,0,I-1,0)
     &                         +      SIGN*C1(K,1,I-1,0)  
  401            CONTINUE
               END IF
               DO 402 K = 1, NUCAB
                  C1(K,I-1,I,0) = HEXPPI(K)*C1(K,I-2,I-1,0)
     &                          + N*BOVERP(K)*C0(K,I-1,I-1,0)
  402          CONTINUE
               IF (I .GE. 3) THEN
                  DO 511 T = 1 + IBTAND(I,1), I - 2, 2
                  TPLUS1 = SIGN*(T + 1.0D0)
                  DO 521 K = 1, NUCAB
                     C1(K,T,I,0) = HEXPPI(K)*C1(K,T-1,I-1,0)
     &                           + N*BOVERP(K)*C0(K,T  ,I-1,0)
     &                           +    TPLUS1*C1(K,T+1,I-1,0)
  521                CONTINUE
  511             CONTINUE
               END IF
            ELSE
c              case n even, C. Villani Uni-Ka, Feb 2005
               IF (IBTAND(I,1) .EQ. 0) THEN
                 DO K = 1, NUCAB
                    C1(K,0,I,0) = N*BOVERP(K)*C0(K,0,I-1,0)
     &                        +      SIGN*C1(K,1,I-1,0)  
                 END DO
               END IF
               DO K = 1, NUCAB
                  C1(K,I,I,0) = HEXPPI(K)*C1(K,I-1,I-1,0)
               END DO
               IF (I .GE. 3) THEN
                  DO T = 1 + IBTAND(I+1,1), I - 1, 2
                  TPLUS1 = SIGN*(T + 1.0D0)
                  DO K = 1, NUCAB
                     C1(K,T,I,0) = HEXPPI(K)*C1(K,T-1,I-1,0)
     &                           + N*BOVERP(K)*C0(K,T  ,I-1,0)
     &                           +    TPLUS1*C1(K,T+1,I-1,0)
                     END DO
                  END DO
               END IF
            END IF
         END IF
         IF (MAX2 .LE. 0) GO TO 1001
         DO 601 J = 1, MAX2
            IJ = I + J
            IF (IJ .EQ. 1) THEN
               IF (IBTAND(N,1) .EQ. 1) THEN
                  DO 701 K = 1, NUCAB
                     C1(K,0,I,J) = N*AOVERP(K)*C0(K,0,I,J-1)
  701             CONTINUE
               ELSE
c                 case n even, C. Villani Uni-Ka, Feb 2005
                  DO K = 1, NUCAB
                     C1(K,1,I,J) = HEXPPI(K)*C1(K,0,0,0)
                  END DO
               END IF
            ELSE IF (IJ .GE. 2) THEN
               IF (IBTAND(N,1) .EQ. 1) THEN
                  IF (IBTAND(IJ,1) .EQ. 1) THEN
                     DO 801 K = 1, NUCAB
                        C1(K,0,I,J) = N*AOVERP(K)*C0(K,0,I,J-1)
     &                              +      SIGN*C1(K,1,I,J-1)
  801                CONTINUE
                  END IF
                  DO 802 K = 1, NUCAB
                     C1(K,IJ-1,I,J) = HEXPPI(K)*C1(K,IJ-2,I,J-1)
     &                              + N*AOVERP(K)*C0(K,IJ-1,I,J-1)
  802             CONTINUE
                  IF (IJ .GE. 3) THEN
                     DO 911 T = 1 + IBTAND(IJ,1), IJ - 2, 2
                        TPLUS1 = SIGN*(T + 1.0D0)
                        DO 921 K = 1, NUCAB
                        C1(K,T,I,J) = HEXPPI(K)*C1(K,T-1,I,J-1)
     &                              + N*AOVERP(K)*C0(K,T  ,I,J-1)
     &                              +    TPLUS1*C1(K,T+1,I,J-1)
  921                   CONTINUE
  911                CONTINUE
                  END IF
               ELSE
c                 case n even, C. Villani Uni-Ka, Feb 2005
                  IF (IBTAND(IJ,1) .EQ. 0) THEN
                     DO K = 1, NUCAB
                        C1(K,0,I,J) = N*AOVERP(K)*C0(K,0,I,J-1)
     &                              +      SIGN*C1(K,1,I,J-1)
                     END DO
                  END IF
                  DO K = 1, NUCAB
                     C1(K,IJ,I,J) = HEXPPI(K)*C1(K,IJ-1,I,J-1)
                  END DO
                  IF (IJ .GE. 3) THEN
                     DO T = 1 + IBTAND(IJ+1,1), IJ - 1, 2
                        TPLUS1 = SIGN*(T + 1.0D0)
                        DO K = 1, NUCAB
                           C1(K,T,I,J) = HEXPPI(K)*C1(K,T-1,I,J-1)
     &                                 + N*AOVERP(K)*C0(K,T  ,I,J-1)
     &                                 +    TPLUS1*C1(K,T+1,I,J-1)
                        END DO
                     END DO
                  END IF
               END IF
            END IF
  601    CONTINUE
 1001 CONTINUE
      IF (IPRINT .LT. 10) RETURN
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      WRITE (LUPRI, 1110)
      WRITE (LUPRI, 1010) MAX1, MAX2
      WRITE (LUPRI, 1030) I120
      IF (IPRINT .LT. 20) RETURN
      DO 2000 I = 0, MAX1
         DO 2100 J = 0, MAX2
            DO 2200 T = 0, I + J
c              introduced n, C. Villani Uni-Ka, Feb 2005
               IF (IBTAND(I + J  - T + N,1) .EQ. 0) THEN
                  WRITE (LUPRI, 1100) WORD, I, J, T
                  WRITE (LUPRI, 1130) (C1(K,T,I,J), K = 1, NUCAB)
               END IF
 2200       CONTINUE
 2100    CONTINUE
 2000 CONTINUE
      RETURN
 1110 FORMAT (/,'  <<<<<<<<<< SUBROUTINE WKEQ30 >>>>>>>>>>',/)
 1010 FORMAT ('  MAX1/2:     ',2I7)
 1030 FORMAT ('  I120:     ',I7)
 1100 FORMAT (/,1X,A4,'(',I1,',',I1,';',I1,')',/)
 1130 FORMAT(1X,6F12.8)
      END
C  /* Deck aminbp */
      SUBROUTINE AMINBP(ABP,NUCAB,TEXPA,TEXPB,HEXPPI,SIGN,IPRINT)
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
#include <implicit.h>
#include <priunit.h> 
      DIMENSION HEXPPI(NUCAB), TEXPA(NUCAB), TEXPB(NUCAB),ABP(NUCAB)
      DO 100 I = 1, NUCAB
         ABP(I) = SIGN * (TEXPA(I) - TEXPB(I)) * HEXPPI(I)
  100 CONTINUE
      IF (IPRINT .GE. 20) THEN
         WRITE(LUPRI,'(/A/)')   ' >>> (A - B) / (A + B) <<< '
         WRITE(LUPRI,'(5F14.8)') (ABP(I), I = 1, NUCAB)
      END IF
      RETURN
      END
C  /* Deck wker12 */
      SUBROUTINE WKER12(RES,A,B,C,P1,P2,P3,GAMFX,II,JJ,IPRINT)
C     Written by Claire C.M. Samson and Wim Klopper 
C     (University of Oslo, 29 October 2001).
#include <implicit.h>
#include <priunit.h>
      DIMENSION RES(NPPX,NTUV), A(NPPX,NTUV), B(NPPX,NTUV),
     &          C(NPPX,NTUV), P1(NPPX), P2(NPPX), P3(NPPX)
#include <r12int.h>
#include <ericom.h>
      IF (II .EQ. 1 .AND. JJ .EQ. 1) THEN
         DO J = 1, NTUV
            DO I = 1, NPPX
               RES(I,J) = P1(I)*A(I,J)*GAMFX
            ENDDO      
         ENDDO  
      ELSE
         DO J = 1, NTUV
            DO I = 1, NPPX
               RES(I,J) = RES(I,J) + P1(I)*A(I,J)*GAMFX
            ENDDO      
         ENDDO  
      END IF
      IF (SLATER) THEN
         IF (INTGAC .NE. 6) GOTO 100
      ELSE
         IF (INTGAC .EQ. 2) GOTO 100
      END IF
      DO J = 1, NTUV
         DO I = 1, NPPX
            RES(I,J) = RES(I,J) + P2(I)*B(I,J)*GAMFX
         ENDDO
      ENDDO 
      IF (SLATER) GOTO 100
      IF (INTGAC .EQ. 4) GOTO 100
      DO J = 1, NTUV
         DO I = 1, NPPX
            RES(I,J) = RES(I,J) + P3(I)*C(I,J)*GAMFX
         ENDDO
      ENDDO 
  100 RETURN
      END
C  /* Deck cr1u12 */
      SUBROUTINE CR1U12(HERINT,HERR12,HCPRIM,L1,M1,N1,L2,M2,N2,
     &                  INDHER,INDHVC,IODDHH,INDHSQ,
     &                  ECOEF,EUV,ETUV,FRSTUV,
     &                  ICOMP1,ICOMP2,AOVERP,IPRINT)
C
C     This subroutine is a variant of CR1TWO for integrals over [T1,r12].
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
C
      PARAMETER (DP5 = 0.5D0)
      INTEGER T, U, V, TUV
      LOGICAL FRSTUV(NTUV34)
C
      DIMENSION HERINT(NPP12,NPRF34,NTUV),HERR12(NPP12,NPRF34,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), INDHVC(0:*),
     &          IODDHH(NRTOP), INDHSQ(NRTOP),
     &          ECOEF(NPP12,0:JMAX1+JMAX2,0:JMAX1,0:JMAX2,3,2),
     &          ETUV(NPP12), EUV(NPP12),
     &          HCPRIM(NPP12,NPRF34,NTUV34),AOVERP(NPP12)
C
#include <ericom.h>
#include <eriao.h>
#include <hertop.h>
#include <r12int.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CR1U12','*',103)
C
      INCT = I120(1) + 1
      INCU = I120(2) + 1
      INCV = I120(3) + 1
      MAXT = L1 + L2
      MAXU = M1 + M2
      MAXV = N1 + N2
      MINT = IBTAND(MAXT,I120(1))
      MINU = IBTAND(MAXU,I120(2))
      MINV = IBTAND(MAXV,I120(3))
      MINTP = IBTAND(MAXT + 1,I120(1))
      MINUP = IBTAND(MAXU + 1,I120(2))
      MINVP = IBTAND(MAXV + 1,I120(3))
C
      IF (IPRINT .GT. 25) THEN
         WRITE(LUPRI,'(/,1X,A,2I5/)')' ICOMP1, ICOMP2', ICOMP1,ICOMP2
         WRITE(LUPRI,'(1X,A,15X,3I5)')' T loop:',MINT,MAXT,INCT
         WRITE(LUPRI,'(1X,A,15X,3I5)')' U loop:',MINU,MAXU,INCU
         WRITE(LUPRI,'(1X,A,15X,3I5)')' V loop:',MINV,MAXV,INCV
      END IF
C
      DO 100 TUV = 1, NTUV34
         FRSTUV(TUV) = .TRUE.
  100 CONTINUE
C
      IF (INTGAC .EQ. 5) THEN
C
C        Expansion coefficients := (a-b)/(a+b)
C        =====================================
C
         DO 205 V = MINV, MAXV, INCV
         DO 205 U = MINU, MAXU, INCU
            DO 215 I = 1, NPP12 
               EUV(I) = ECOEF(I,V,N1,N2,3,1) * ECOEF(I,U,M1,M2,2,1)
  215       CONTINUE
            DO 305 T = MINT, MAXT, INCT
               DO 315 I = 1, NPP12
                  ETUV(I) = ECOEF(I,T,L1,L2,1,1)*EUV(I)
                  ETUV(I) = DP5*ETUV(I)*AOVERP(I)
  315          CONTINUE
C
               ITUV = INDHER(T,U,V)
               INDST = INDHSQ(INDHER(T+2,U,V))
               INDSU = INDHSQ(INDHER(T,U+2,V))
               INDSV = INDHSQ(INDHER(T,U,V+2))
               DO 405 TUV = 1, NTUV34
                  IF (IODDHH(ITUV) .EQ. IODDHH(TUV)) THEN
                     INDT = INDST + INDHSQ(TUV)
                     INDU = INDSU + INDHSQ(TUV)
                     INDV = INDSV + INDHSQ(TUV)
                     INDT = INDHVC(INDT)
                     INDU = INDHVC(INDU)
                     INDV = INDHVC(INDV)
                     IF (FRSTUV(TUV)) THEN
                        FRSTUV(TUV) = .FALSE.
                        DO 505 J = 1, NPRF34
                        DO 505 I = 1, NPP12
                           HCPRIM(I,J,TUV) = ETUV(I)*
     &                                       (HERR12(I,J,INDT)
     &                                       +HERR12(I,J,INDU)
     &                                       +HERR12(I,J,INDV))
  505                   CONTINUE
                     ELSE
                        DO 605 J = 1, NPRF34
                        DO 605 I = 1, NPP12
                           HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV) + ETUV(I)*
     &                                       (HERR12(I,J,INDT)
     &                                       +HERR12(I,J,INDU)
     &                                       +HERR12(I,J,INDV))
  605                   CONTINUE
                     END IF
                  END IF
  405          CONTINUE
C
  305       CONTINUE
  205    CONTINUE

      ELSE
C
C        Expansion coefficients := (a-b)/(a+b)
C        =====================================
C
         DO 200 V = MINV, MAXV, INCV
         DO 200 U = MINU, MAXU, INCU
            DO 210 I = 1, NPP12 
               EUV(I) = ECOEF(I,V,N1,N2,3,1) * ECOEF(I,U,M1,M2,2,1)
  210       CONTINUE
            DO 300 T = MINT, MAXT, INCT
               DO 310 I = 1, NPP12
                  ETUV(I) = ECOEF(I,T,L1,L2,1,1)*EUV(I)
                  ETUV(I) = ETUV(I)*AOVERP(I)
  310          CONTINUE
               ITUV = INDHER(T,U,V)
               INDS = INDHSQ(ITUV)
               DO 400 TUV = 1, NTUV34
                  IF (IODDHH(ITUV) .EQ. IODDHH(TUV)) THEN
#if defined (SYS_AIX)
C                    code due to AIX xlf version 2.2 bug
                     INDT = INDS + INDHSQ(TUV)
                     INDT = INDHVC(INDT)
#else
                     INDT = INDHVC(INDS + INDHSQ(TUV))
#endif
                     IF (FRSTUV(TUV)) THEN
                        FRSTUV(TUV) = .FALSE.
                        DO 500 J = 1, NPRF34
                        DO 500 I = 1, NPP12
                           HCPRIM(I,J,TUV) = ETUV(I)*HERINT(I,J,INDT)
  500                   CONTINUE
                     ELSE
                        DO 600 J = 1, NPRF34
                        DO 600 I = 1, NPP12
                           HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV)
     &                             + ETUV(I)*HERINT(I,J,INDT)
  600                   CONTINUE
                     END IF
                  END IF
  400          CONTINUE
C
  300       CONTINUE
  200    CONTINUE
      END IF
C
C     Expansion coefficients := Rx*Ey*Ez
C     ==================================
C
      DO 201 V = MINV, MAXV, INCV
      DO 201 U = MINU, MAXU, INCU
         DO 211 I = 1, NPP12 
            EUV(I) = ECOEF(I,V,N1,N2,3,1) 
     &             * ECOEF(I,U,M1,M2,2,1)
  211    CONTINUE
         DO 301 T = MINTP, MAXT, INCT
            DO 311 I = 1, NPP12
               ETUV(I) = ECOEF(I,T,L1,L2,1,2) * EUV(I)
  311       CONTINUE
C
            ITUV = INDHER(T+1,U,V)
            INDS = INDHSQ(ITUV)
            DO 401 TUV = 1, NTUV34
            IF (IODDHH(ITUV) .EQ. IODDHH(TUV)) THEN
#if defined (SYS_AIX)
C              code due to AIX xlf version 2.2 bug
               INDT = INDS + INDHSQ(TUV)
               INDT = INDHVC(INDT)
#else
               INDT = INDHVC(INDS + INDHSQ(TUV))
#endif
               IF (FRSTUV(TUV)) THEN
                  FRSTUV(TUV) = .FALSE.
                  DO 501 J = 1, NPRF34
                  DO 501 I = 1, NPP12
                     HCPRIM(I,J,TUV) = ETUV(I)*HERR12(I,J,INDT)
  501             CONTINUE
               ELSE
                  DO 601 J = 1, NPRF34
                  DO 601 I = 1, NPP12
                     HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV)
     &                       + ETUV(I)*HERR12(I,J,INDT)
  601             CONTINUE
               END IF
            END IF
  401       CONTINUE
C
  301    CONTINUE
  201 CONTINUE
C
C     Expansion coefficients := Ex*Ry*Ez
C     ==================================
C
      DO 202 V = MINV, MAXV, INCV
      DO 202 U = MINUP, MAXU, INCU
         DO 212 I = 1, NPP12 
            EUV(I) = ECOEF(I,V,N1,N2,3,1) 
     &             * ECOEF(I,U,M1,M2,2,2)
  212    CONTINUE
         DO 302 T = MINT, MAXT, INCT
            DO 312 I = 1, NPP12
               ETUV(I) = ECOEF(I,T,L1,L2,1,1) * EUV(I)
  312       CONTINUE
C
            ITUV = INDHER(T,U+1,V)
            INDS = INDHSQ(ITUV)
            DO 402 TUV = 1, NTUV34
            IF (IODDHH(ITUV) .EQ. IODDHH(TUV)) THEN
#if defined (SYS_AIX)
C              code due to AIX xlf version 2.2 bug
               INDT = INDS + INDHSQ(TUV)
               INDT = INDHVC(INDT)
#else
               INDT = INDHVC(INDS + INDHSQ(TUV))
#endif
               IF (FRSTUV(TUV)) THEN
                  FRSTUV(TUV) = .FALSE.
                  DO 502 J = 1, NPRF34
                  DO 502 I = 1, NPP12
                     HCPRIM(I,J,TUV) = ETUV(I)*HERR12(I,J,INDT)
  502             CONTINUE
               ELSE
                  DO 602 J = 1, NPRF34
                  DO 602 I = 1, NPP12
                     HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV)
     &                       + ETUV(I)*HERR12(I,J,INDT)
  602             CONTINUE
               END IF
            END IF
  402       CONTINUE
C
  302    CONTINUE
  202 CONTINUE
C
C     Expansion coefficients := Ex*Ey*Rz
C     ==================================
C
      DO 203 V = MINVP, MAXV, INCV
      DO 203 U = MINU, MAXU, INCU
         DO 213 I = 1, NPP12 
            EUV(I) = ECOEF(I,V,N1,N2,3,2) 
     &             * ECOEF(I,U,M1,M2,2,1)
  213    CONTINUE
         DO 303 T = MINT, MAXT, INCT
            DO 313 I = 1, NPP12
               ETUV(I) = ECOEF(I,T,L1,L2,1,1) * EUV(I)
  313       CONTINUE
C
            ITUV = INDHER(T,U,V+1)
            INDS = INDHSQ(ITUV)
            DO 403 TUV = 1, NTUV34
            IF (IODDHH(ITUV) .EQ. IODDHH(TUV)) THEN
#if defined (SYS_AIX)
C              code due to AIX xlf version 2.2 bug
               INDT = INDS + INDHSQ(TUV)
               INDT = INDHVC(INDT)
#else
               INDT = INDHVC(INDS + INDHSQ(TUV))
#endif
               IF (FRSTUV(TUV)) THEN
                  FRSTUV(TUV) = .FALSE.
                  DO 503 J = 1, NPRF34
                  DO 503 I = 1, NPP12
                     HCPRIM(I,J,TUV) = ETUV(I)*HERR12(I,J,INDT)
  503             CONTINUE
               ELSE
                  DO 603 J = 1, NPRF34
                  DO 603 I = 1, NPP12
                     HCPRIM(I,J,TUV) = HCPRIM(I,J,TUV)
     &                       + ETUV(I)*HERR12(I,J,INDT)
  603             CONTINUE
               END IF
            END IF
  403       CONTINUE
C
  303    CONTINUE
  203 CONTINUE
C
      RETURN
      END
C  /* Deck coeff */
      SUBROUTINE COEFF(ARGE,ARGF,ALPHA,BETA,PRE1,PRE2,
     &                 PRE3,FACINT,PP,PQ,IPRINT)
C     Written by Claire C.M. Samson (University of Oslo, 29 October 2001).
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D2 = 2.D0, D20 = 20.0D0,
     &            D3 = 3.0D0, D5 = 5.0D0, DP05 = 0.5D0, DP15 = 1.5D0, 
     &            DP25 = 2.5D0, DP375 = 3.75D0, D1M = -1.0D0, 
     &            D2M = -2.0D0, DP05M = -0.5D0, D15 = 15.0D0,
     &            D4 = 4.0D0 )
#include <pi.h>
#include <priunit.h>
#include <r12int.h>
#include <ericom.h>
C
      DIMENSION ARGE(NPPX), ARGF(NPPX),PP(NPPX),PQ(NPPX),
     &          ALPHA(NPPX), BETA(NPPX), PRE1(NPPX),PRE2(NPPX),
     &          PRE3(NPPX), FACINT(NPPX)
C
      IF (INTGAC .EQ. 2 .OR. INTGAC .EQ. 3 .OR. INTGAC .EQ. 5) THEN
         VAR = GAMAC                                         
      ELSE IF (INTGAC .EQ. 4 .OR. INTGAC .EQ. 6) THEN
         VAR = GAMAC+GAMAD
         VSQ = GAMAC*GAMAD
      ELSE
         CALL QUIT('ILLEGAL INTGAC IN COEFF')
      END IF
C
C     Calculate alpha and beta
C
      DO I = 1, NPPX
         ZZ = PP(I)*PQ(I)/(PP(I)+PQ(I))
         ALPHA(I) = ZZ
         BETA(I) = ZZ/(ZZ+VAR)
      END DO
C
C     Calculate the prefactors
C
      IF (SLATER) THEN
      IF (INTGAC .EQ. 2) THEN
         DO I = 1, NPPX
            PRE1(I) = FACINT(I)*BETA(I)*(PP(I)+PQ(I))**DP05M
            PRE2(I) = D0
            PRE3(I) = D0
         END DO
      ELSE IF (INTGAC .EQ. 3 .OR. INTGAC .EQ. 4 .OR. INTGAC .EQ. 5) THEN
C**as if INTGAC=2
         DO I = 1, NPPX
            PRE1(I) = FACINT(I)*DP05*SQRT(PI*BETA(I)**3/(PP(I)*PQ(I)))
            PRE2(I) = D0
            PRE3(I) = D0
         END DO
      ELSE IF (INTGAC .EQ. 6) THEN
C**as if INTGAC=4 ** times 4*VSQ
         DO I = 1, NPPX
            ZZ = FACINT(I)*D2*VSQ*SQRT(PI*BETA(I)**3/(PP(I)*PQ(I)))
     &         / (ALPHA(I)+VAR)
            YY = ALPHA(I)*BETA(I)
            PRE1(I) = ZZ*DP15
            PRE2(I) = ZZ*YY
            PRE3(I) = D0
         END DO
      END IF
C
C     Calculate the arguments of the exponentials and the Boys functions.
C
      IF (INTGAD .EQ. 1) THEN
         DO I = 1, NPPX
            ARGE(I) = VAR*BETA(I)
            ARGF(I) = ALPHA(I)*BETA(I) 
         END DO
      ELSE
         IF (INTGAD .NE. 0) THEN
            CALL QUIT('ILLEGAL INTGAD IN COEFF')
         END IF
         DO 501 I = 1, NPPX
            ARGE(I) = VAR*BETA(I)
  501    END DO
      END IF
      ELSE
      IF (INTGAC .EQ. 2) THEN
         DO I = 1, NPPX
            PRE1(I) = FACINT(I)*DP05*SQRT(PI*BETA(I)**3/(PP(I)*PQ(I)))
            PRE2(I) = D0
            PRE3(I) = D0  
         END DO
      ELSE IF (INTGAC .EQ. 3 .OR. INTGAC .EQ. 5) THEN
         DO I = 1, NPPX
            ZZ = FACINT(I)*DP05*(PP(I)+PQ(I))**DP05M*BETA(I)
     &           / (ALPHA(I)+VAR)
            PRE1(I) = ZZ
            PRE2(I) = ZZ
            PRE3(I) = ZZ*D2*ALPHA(I)*BETA(I)
         END DO
      ELSE IF (INTGAC .EQ. 4) THEN
         DO I = 1, NPPX
            ZZ = FACINT(I)*DP05*SQRT(PI*BETA(I)**3/(PP(I)*PQ(I)))
     &         / (ALPHA(I)+VAR)
            YY = ALPHA(I)*BETA(I)
            PRE1(I) = ZZ*DP15
            PRE2(I) = ZZ*YY
            PRE3(I) = D0
         END DO
      ELSE IF (INTGAC .EQ. 6) THEN
         DO I = 1, NPPX
            ZZ = FACINT(I)*DP05*SQRT(PI*BETA(I)**3/(PP(I)*PQ(I)))
            YY = VAR/(ALPHA(I)+VAR)
            WW = VSQ/(ALPHA(I)+VAR)**2
            XX = ALPHA(I)*BETA(I)
            PRE1(I) = ZZ*(D1-D3*YY+D15*WW)        
            PRE2(I) = ZZ*(D20*WW-D2*YY)*XX
            PRE3(I) = ZZ*D4*WW*XX**2
         END DO
      END IF
C
C     Calculate the arguments of the exponentials and the Boys functions.
C
      IF (INTGAD .EQ. 1) THEN
         DO I = 1, NPPX
            ARGE(I) = VAR*BETA(I)
            ARGF(I) = ALPHA(I)*BETA(I) 
         END DO
      ELSE
         IF (INTGAD .EQ. 2) THEN
            DO I = 1, NPPX
C**this is equal to ALPHA(I)
               ARGE(I) = (VAR+ALPHA(I))*BETA(I)
            END DO  
         ELSE
            DO 500 I = 1, NPPX
               ARGE(I) = VAR*BETA(I)
  500       END DO  
         END IF
      END IF
      END IF
C
C     Print section
C
      IF (IPRINT .GE. 10) THEN
        CALL AROUND('Output from COEFF')
        WRITE(LUPRI,'(3X,5(10X,A4))') 'PRE1','PRE2','PRE3','ARGE','ARGF'
        IF (INTGAD .EQ. 1) THEN
          DO I = 1, NPPX
            WRITE(LUPRI,'(I3,5F14.9)') I,
     &      PRE1(I),PRE2(I),PRE3(I),ARGE(I),ARGF(I)
          END DO
        ELSE
          DO I = 1, NPPX
            WRITE(LUPRI,'(I3,5F14.9)') I,
     &      PRE1(I),PRE2(I),PRE3(I),ARGE(I)
          END DO
        END IF
      END IF
      RETURN
      END
C  /* Deck tester */
      SUBROUTINE TESTER(RES,RO000,PRE1,PRE2,PRE3,
     &                  ALPHA,BETA,PQX,PQY,PQZ)
C     Written by Claire C.M. Samson (University of Oslo, 29 October 2001).
#include <implicit.h>
      PARAMETER ( D1=1.0D0, D2=2.0D0, D2M= -2.0D0 ) 
#include <pi.h>
#include <r12int.h>
#include <ericom.h>
C
      DIMENSION RES(NPPX),RO000(NPPX,0:JMAX), 
     &          PRE1(NPPX),PRE2(NPPX),PRE3(NPPX),ALPHA(NPPX),
     &          BETA(NPPX),PQX(NPPX),PQY(NPPX),PQZ(NPPX)
      DO I=1, NPPX
        R2=PQX(I)**2+PQY(I)**2+PQZ(I)**2
        DINT1=PRE1(I)*D2M*(ALPHA(I)+GAMAC)*BETA(I)*PQZ(I)*
     &       EXP(-(ALPHA(I)+GAMAC)*BETA(I)*R2)
        ARGE=GAMAC*BETA(I)
        ARGF=ALPHA(I)*BETA(I)*R2
        DINT2=PRE2(I)*(D2M*GAMAC*BETA(I)*PQZ(I)*EXP(-ARGE*R2)*
     &        RO000(I,0)
     &      +D2M*ALPHA(I)*BETA(I)*PQZ(I)*EXP(-ARGE*R2)*
     &        RO000(I,1))
        DINT3=PRE3(I)*(D2*PQZ(I)*EXP(-ARGE*R2)*RO000(I,0)
     &      +D2M*GAMAC*BETA(I)*PQZ(I)*R2*EXP(-ARGE*R2)*RO000(I,0)
     &      +D2M*ALPHA(I)*BETA(I)*PQZ(I)*R2*EXP(-ARGE*R2)*
     &       RO000(I,1))
        RES(I)=DINT1+DINT2+DINT3
        DINT0=PRE1(I)*(EXP(-(ALPHA(I)+GAMAC)*BETA(I)*R2)+
     &        (D1+D2*ALPHA(I)*BETA(I)*R2)*RO000(I,0)*
     &         EXP(-ARGE*R2)) 
      END DO
      RETURN
      END 
C  /* Deck erigom */
      SUBROUTINE ERIGOM(RO000,RM000,ARGEXP,ARGF,PQXYZ,WORK,
     &                  LWORK,IPRINT)
C     Written by Claire C.M. Samson (University of Oslo, 29 October 2001).
#include <implicit.h>
C
      DIMENSION RO000(NPPX,0:JMAX), RM000(NPPX,0:JMAX),
     &          ARGEXP(NPPX), ARGF(NPPX),PQXYZ(NPPX,3), WORK(LWORK)
C
#include <ericom.h>
C
C     Allocate ERIGOM
C
      KNDADR = 1 
      KLAST = KNDADR + NPPX
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIGOM',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      CALL ERIGO(RO000,RM000,PQXYZ,ARGEXP,ARGF,WORK(KNDADR),
     &           WORK(KLAST),LWRK,IPRINT)
C
      RETURN
      END
C  /* Deck erigo */
      SUBROUTINE ERIGO(RO000,RM000,PQXYZ,ARGEXP,ARGF,INDADR,
     &                 WORK,LWORK,IPRINT)
C     Written by Claire C.M. Samson (University of Oslo, 29 October 2001).
#include <implicit.h>
#include <iratdef.h>
      PARAMETER ( D2 = 2.D0 )
C
      DIMENSION RO000(NPPX,0:JMAX), RM000(NPPX,0:JMAX),
     &          PQXYZ(NPPX,3), ARGEXP(NPPX),ARGF(NPPX),
     &          INDADR(NPPX),WORK(LWORK)
C
#include <ericom.h>
C
      CALL DZERO(RO000(1,0),(JMAX+1)*NPPX)
      CALL DZERO(RM000(1,0),(JMAX+1)*NPPX)
C
C        Calculate gamma function with ARGF
C        ==================================
C
C        Allocate GETGAM
C
C
      KINDAD = 1
      KWVALS = KINDAD +  (3*NPPX - 1)/IRAT + 1
      KFJW   = KWVALS +  3*NPPX
      KREXPW = KFJW   +    NPPX*(JMAX + 1)
      KARGFI = KREXPW +    NPPX
      KLAST  = KARGFI +    NPPX
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIGO',' ',KLAST,LWORK)
C
      DO 100 I = 1, NPPX
         WORK(KARGFI+I-1) = ARGF(I)*(PQXYZ(I,1)**D2+PQXYZ(I,2)**D2
     &                      +PQXYZ(I,3)**D2)
      INDADR(I) = I
  100 CONTINUE
C
      CALL GETGAM(NPPX,INDADR,WORK(KARGFI),RO000,JMAX,NPPX,WORK(KFJW),
     &            WORK(KINDAD),WORK(KWVALS),WORK(KREXPW),IPRINT)
C
C        Scale gamma function
C        ====================
C
         DO 200 J = 0, JMAX
         DO 200 I = 1, NPPX
            RO000(I,J) = EXP(-ARGEXP(I)*(PQXYZ(I,1)**D2+PQXYZ(I,2)**D2
     &                      +PQXYZ(I,3)**D2))*RO000(I,J)
            RM000(I,J) = (PQXYZ(I,1)**D2+PQXYZ(I,2)**D2+PQXYZ(I,3)**D2)
     &                  *RO000(I,J)
  200    CONTINUE
C
      RETURN
      END
C  /* Deck eriomtuv */
      SUBROUTINE ERIOMTUV(HERINTO,HERINTM,RO000,RM000,ARGEXP,ARGF,
     &                    PQXYZ,INDHER,WORK,LWORK,IPRINT)
C     Written by Claire C.M. Samson (University of Oslo, 29 October 2001).
#include <implicit.h>
#include <priunit.h>
      INTEGER T, U, V, TUV
      DIMENSION HERINTO(NPPX,NTUV), HERINTM(NPPX,NTUV),
     &          RO000(NPPX,0:JMAX), RM000(NPPX,0:JMAX),
     &          ARGEXP(NPPX), ARGF(NPPX),
     &          PQXYZ(NPPX,3),INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          WORK(LWORK)
#include <ericom.h>
#include <hertop.h>
C
C     The final O(T,U,V) and M(T,U,V) integrals are arranged as follows:
C
C     O(000)
C     O(100) O(010) O(001)
C     O(200) O(110) O(101) O(020) O(011) O(002)
C     O(300) O(210) O(201) O(120) O(111) O(102) O(030) O(021) O(012)
C                                                             O(003)
C     Special case JMAX = 0
C     =====================
C
      IF (JMAX .EQ. 0) THEN
         CALL DCOPY(NPPX,RO000,1,HERINTO,1)
         CALL DCOPY(NPPX,RM000,1,HERINTM,1)
      ELSE
C
C        Allocate work space
C        ===================
C
         KHRWRKO = 1
         KHRWRKM = KHRWRKO + NTUV*NPPX
         KLAST   = KHRWRKM + NTUV*NPPX
         IF (KLAST .GT. LWORK) CALL STOPIT('HERIOM',' ',KLAST,LWORK)
C
C        Recursion loop for Hermite integrals O and M
C        ============================================
C
         IPQX = IPQ0(1) + 1
         IPQY = IPQ0(2) + 1
         IPQZ = IPQ0(3) + 1
         CALL DZERO(HERINTO,NPPX*NTUV)
         CALL DZERO(HERINTM,NPPX*NTUV)
         DO 100 JVAL = 1, JMAX
            IF (MOD(JMAX-JVAL,2).EQ.0) THEN
             CALL ERIHRCOM(HERINTO,WORK(KHRWRKO),HERINTM,WORK(KHRWRKM),
     &                     JVAL,RO000,RM000,PQXYZ(1,1),PQXYZ(1,2),
     &                     PQXYZ(1,3),INDHER,JMAX,NPPX,NTUV,
     &                     ARGEXP,ARGF,IPQX,IPQY,IPQZ)
            ELSE
             CALL ERIHRCOM(WORK(KHRWRKO),HERINTO,WORK(KHRWRKM),HERINTM,
     &                     JVAL,RO000,RM000,PQXYZ(1,1),PQXYZ(1,2),
     &                     PQXYZ(1,3),INDHER,JMAX,NPPX,NTUV,
     &                     ARGEXP,ARGF,IPQX,IPQY,IPQZ)
            END IF
  100    CONTINUE
      END IF
C
C     Print section
C     =============
C
      IF (IPRINT .GE. 10) THEN
         CALL TITLER('Output from ERIOMTUV','*',103)
         WRITE (LUPRI,'(2X,A,I10)') 'JMAX  ', JMAX
         WRITE (LUPRI,'(2X,A,I10)') 'NPPX  ', NPPX
         WRITE (LUPRI,'(2X,A,I10)') 'NTUV  ', NTUV
         IF (IPRINT .GE. 10) THEN
            CALL HEADER('Hermite integrals RO(t,u,v) and RM(t,u,v)',1)
            DO 300 J = 0, JMAX
              DO 320 T = J, 0, -1
                DO 330 U = J - T, 0, -1
                  V = J - T - U
                  TUV = INDHER(T,U,V)
C                 IF (IODDHR(TUV) .EQ. 0) THEN
                    WRITE (LUPRI,'(2X,3(A,I1),A,2X,5G12.6/,
     &                                                 (12X,5G12.6))')
     &              'RO(',T,',',U,',',V,')', (HERINTO(I,TUV),I=1,NPPX)
                    WRITE (LUPRI,'(2X,3(A,I1),A,2X,5G12.6/,
     &                                                 (12X,5G12.6))')
     &              'RM(',T,',',U,',',V,')', (HERINTM(I,TUV),I=1,NPPX)
                  WRITE (LUPRI,'()')
C                 END IF
  330           CONTINUE
  320         CONTINUE
  300      CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck erihrcom */
      SUBROUTINE ERIHRCOM(CURO,OLDO,CURM,OLDM,JVAL,RO000,RM000,
     &                   PQX,PQY,PQZ,INDHER,JMAX,NPPX,NTUV,
     &                   ARGE,ARGF,IPQX,IPQY,IPQZ)
C     Written by Claire C.M. Samson, W. Klopper and Trygve Helgaker
C     (University of Oslo, 29 October 2001).
#include <implicit.h>
      PARAMETER (D0 = 0D0)
      INTEGER T, U, V, TUV
      LOGICAL PQXGT0, PQYGT0, PQZGT0
      DIMENSION CURO(NPPX,NTUV), OLDO(NPPX,NTUV),
     &          CURM(NPPX,NTUV), OLDM(NPPX,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          PQX(NPPX), PQY(NPPX), PQZ(NPPX),
     &          RO000(NPPX,0:JMAX), RM000(NPPX,0:JMAX),
     &          ARGE(NPPX), ARGF(NPPX)
#include <hertop.h>
#include <ibtfun.h>
C
      PQXGT0 = IPQX .EQ. 1
      PQYGT0 = IPQY .EQ. 1
      PQZGT0 = IPQZ .EQ. 1
C
C     O(0,0,0) and M(0,0,0)
C
      CALL DCOPY(NPPX,RO000(1,JMAX-JVAL),1,CURO,1)
      CALL DCOPY(NPPX,RM000(1,JMAX-JVAL),1,CURM,1)
C
C     JVAL = 1
C     ========
C
      IF (JVAL .EQ. 1) THEN
         TMIN1 = D0
         UMIN1 = D0
         VMIN1 = D0
C
         CALL ERIMKO(CURO(1,2),ARGE,ARGF,TMIN1,PQX,PQXGT0,
     &               CURO(1,1),CURO(1,1),CURO(1,1),RO000(1,JMAX))
         CALL ERIMKO(CURO(1,3),ARGE,ARGF,UMIN1,PQY,PQYGT0,
     &               CURO(1,1),CURO(1,1),CURO(1,1),RO000(1,JMAX))
         CALL ERIMKO(CURO(1,4),ARGE,ARGF,VMIN1,PQZ,PQZGT0,
     &               CURO(1,1),CURO(1,1),CURO(1,1),RO000(1,JMAX))
C
         CALL ERIMKM(CURM(1,2),ARGE,ARGF,TMIN1,PQX,PQXGT0,
     &               CURO(1,1),CURO(1,1),CURM(1,1),CURM(1,1),
     &               CURM(1,1),RM000(1,JMAX))
         CALL ERIMKM(CURM(1,3),ARGE,ARGF,UMIN1,PQY,PQYGT0,
     &               CURO(1,1),CURO(1,1),CURM(1,1),CURM(1,1),
     &               CURM(1,1),RM000(1,JMAX))
         CALL ERIMKM(CURM(1,4),ARGE,ARGF,VMIN1,PQZ,PQZGT0,
     &               CURO(1,1),CURO(1,1),CURM(1,1),CURM(1,1),
     &               CURM(1,1),RM000(1,JMAX))
C
C     JVAL > 1
C     ========
C
      ELSE
C
C        O(T,0,0) and M(T,0,0) for T > 0
C
         DO T = IPQX, JVAL, IPQX
            TMIN1 = T - 1.0D0
            TUV   = INDHER(T  ,0,0)
            M1T   = INDHER(T-1,0,0)
            M2T   = INDHER(T-2,0,0)
            CALL ERIMKO(CURO(1,TUV),ARGE,ARGF,TMIN1,PQX,PQXGT0,
     &                  CURO(1,M2T),CURO(1,M1T),
     &                  OLDO(1,M2T),OLDO(1,M1T))
            CALL ERIMKM(CURM(1,TUV),ARGE,ARGF,TMIN1,PQX,PQXGT0,
     &                  CURO(1,M2T),CURO(1,M1T),
     &                  CURM(1,M2T),CURM(1,M1T),
     &                  OLDM(1,M2T),OLDM(1,M1T))
         END DO
C
C        O(T,U,0) and M(T,U,0) for U > 0 and all T
C
         DO U = IPQY, JVAL,     IPQY
         DO T = 0,    JVAL - U, IPQX
            UMIN1 = U - 1.0D0
            TUV   = INDHER(T,U  ,0)
            M1U   = INDHER(T,U-1,0)
            M2U   = INDHER(T,U-2,0)
            CALL ERIMKO(CURO(1,TUV),ARGE,ARGF,UMIN1,PQY,PQYGT0,
     &                  CURO(1,M2U),CURO(1,M1U),
     &                  OLDO(1,M2U),OLDO(1,M1U))
            CALL ERIMKM(CURM(1,TUV),ARGE,ARGF,UMIN1,PQY,PQYGT0,
     &                  CURO(1,M2U),CURO(1,M1U),
     &                  CURM(1,M2U),CURM(1,M1U),
     &                  OLDM(1,M2U),OLDM(1,M1U))
         END DO
         END DO
C
C        O(T,U,V) and M(T,U,V) for V > 0 and all (T,U)
C
         DO V = IPQZ, JVAL,         IPQZ
         DO U = 0,    JVAL - V,     IPQY
         DO T = 0,    JVAL - U - V, IPQX
            VMIN1 = V - 1.0D0
            TUV   = INDHER(T,U,V  )
            M1V   = INDHER(T,U,V-1)
            M2V   = INDHER(T,U,V-2)
            CALL ERIMKO(CURO(1,TUV),ARGE,ARGF,VMIN1,PQZ,PQZGT0,
     &                  CURO(1,M2V),CURO(1,M1V),
     &                  OLDO(1,M2V),OLDO(1,M1V))
            CALL ERIMKM(CURM(1,TUV),ARGE,ARGF,VMIN1,PQZ,PQZGT0,
     &                  CURO(1,M2V),CURO(1,M1V),
     &                  CURM(1,M2V),CURM(1,M1V),
     &                  OLDM(1,M2V),OLDM(1,M1V))
         END DO
         END DO
         END DO
      END IF
      RETURN
      END
C  /* Deck erimko */
      SUBROUTINE ERIMKO(R,A,B,T,PQX,PQXGT0,X1,X2,X3,X4)
C     Written by Wim Klopper (University of Oslo, 29 October 2001).
#include <implicit.h>
      PARAMETER (D2 = 2.D0, D2M = -D2)
      LOGICAL PQXGT0
      DIMENSION R(*),A(*),B(*),X1(*),X2(*),X3(*),X4(*),PQX(*)
      CALL ERICRX(R,D2M,A,T,X1,PQX,X2,PQXGT0)
      CALL ERICRP(R,D2M,B,T,X3,PQX,X4,PQXGT0)
      RETURN
      END
C  /* Deck erimkm */
      SUBROUTINE ERIMKM(R,A,B,T,PQX,PQXGT0,X1,X2,X3,X4,X5,X6)
C     Written by Wim Klopper (University of Oslo, 29 October 2001).
#include <implicit.h>
      PARAMETER (D2 = 2.D0, D2M = -D2)
      LOGICAL PQXGT0
      DIMENSION R(*),A(*),B(*),PQX(*),
     &          X1(*),X2(*),X3(*),X4(*),X5(*),X6(*)
      CALL ERICRI(R,D2,T,X1,PQX,X2,PQXGT0)
      CALL ERICRP(R,D2M,A,T,X3,PQX,X4,PQXGT0)
      CALL ERICRP(R,D2M,B,T,X5,PQX,X6,PQXGT0)
      RETURN
      END
C  /* Deck erimkf */
      SUBROUTINE ERIMKF(R,A,T,PQX,PQXGT0,X1,X2)
C     Written by Wim Klopper (University of Oslo, 29 October 2001).
#include <implicit.h>
      PARAMETER (D2 = 2.D0, D2M = -D2)
      LOGICAL PQXGT0
      DIMENSION R(*),A(*),X1(*),X2(*),PQX(*)
      CALL ERICRX(R,D2M,A,T,X1,PQX,X2,PQXGT0)
      RETURN
      END
C  /* Deck erimkg */
      SUBROUTINE ERIMKG(R,A,T,PQX,PQXGT0,X1,X2,X3,X4)
C     Written by Wim Klopper (University of Oslo, 29 October 2001).
#include <implicit.h>
      PARAMETER (D2 = 2.D0, D2M = -D2)
      LOGICAL PQXGT0
      DIMENSION R(*),A(*),X1(*),X2(*),X3(*),X4(*),PQX(*)
      CALL ERICRI(R,D2,T,X1,PQX,X2,PQXGT0)
      CALL ERICRP(R,D2M,A,T,X3,PQX,X4,PQXGT0)
      RETURN
      END
C  /* Deck erimkh */
      SUBROUTINE ERIMKH(R,A,T,PQX,PQXGT0,X1,X2,X3,X4)
C     Written by Wim Klopper (University of Oslo, 29 October 2001).
#include <implicit.h>
      PARAMETER (D2 = 2.D0, D2M = -D2, D4 = 4D0)
      LOGICAL PQXGT0
      DIMENSION R(*),A(*),X1(*),X2(*),X3(*),X4(*),PQX(*)
      CALL ERICRI(R,D4,T,X1,PQX,X2,PQXGT0)
      CALL ERICRP(R,D2M,A,T,X3,PQX,X4,PQXGT0)
      RETURN
      END
C  /* Deck erifghtuv */
      SUBROUTINE ERIFGHTUV(CURF,CURG,CURH,PQX,PQY,PQZ,INDHER,
     &                     ARGE,IPRINT)
C     Written by Claire C.M. Samson (University of Oslo, 29 October 2001).
#include <implicit.h>
#include <priunit.h>
      INTEGER T, U, V, TUV
      LOGICAL PQXGT0, PQYGT0, PQZGT0, DOG, DOH
      DIMENSION CURF(NPPX,NTUV), CURG(NPPX,NTUV), CURH(NPPX,NTUV),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          PQX(NPPX), PQY(NPPX), PQZ(NPPX),
     &          ARGE(NPPX)
#include <hertop.h>
#include <ericom.h>
#include <r12int.h>
#include <ibtfun.h>
C
      CALL DZERO(CURF,NPPX*NTUV)
      CALL DZERO(CURG,NPPX*NTUV)
      CALL DZERO(CURH,NPPX*NTUV)
C
      IF (SLATER) THEN
         DOH = .FALSE.
         DOG = INTGAC .EQ. 6
      ELSE
         DOH = INTGAC .EQ. 6
         DOG = INTGAC .EQ. 4 .OR. DOH
      ENDIF
C
      IPQX = IPQ0(1) + 1
      IPQY = IPQ0(2) + 1
      IPQZ = IPQ0(3) + 1
C
      PQXGT0 = IPQX .EQ. 1
      PQYGT0 = IPQY .EQ. 1
      PQZGT0 = IPQZ .EQ. 1
C
      M1T = 1 ! to avoid INDHER(-1,0,0)
      M2T = 1 ! to avoid INDHER(-2:-1,0,0)
      DO T = 0, JMAX, IPQX
         TMIN1 = T - 1.0D0
         IF (T .GT. 0) M1T = INDHER(T-1,0,0)
         IF (T .GT. 1) M2T = INDHER(T-2,0,0)
         M1U = 1 ! to avoid INDHER(T,-1,0)
         M2U = 1 ! to avoid INDHER(T,-2:-1,0)
      DO U = 0, JMAX - T, IPQY
         UMIN1 = U - 1.0D0
         IF (U .GT. 0) M1U = INDHER(T,U-1,0)
         IF (U .GT. 1) M2U = INDHER(T,U-2,0)
         M1V = 1 ! to avoid INDHER(T,U,-1)
         M2V = 1 ! to avoid INDHER(T,U,-2:-1)
      DO V = 0, JMAX - T - U , IPQZ
         VMIN1 = V - 1.0D0
         IF (V .GT. 0) M1V = INDHER(T,U,V-1)
         IF (V .GT. 1) M2V = INDHER(T,U,V-2)
C
         TUV = INDHER(T,U,V)
C
C        (0,0,0)
C
         IF (T.EQ.0 .AND. U.EQ.0 .AND. V.EQ.0) THEN
            DO I = 1, NPPX
               R2 = PQX(I)*PQX(I) + PQY(I)*PQY(I) + PQZ(I)*PQZ(I)
               CURH(I,1) = R2
               CURF(I,1) = EXP(-ARGE(I)*R2)
            END DO
            IF (DOG) THEN
               DO I = 1, NPPX
                  CURG(I,1) = CURF(I,1) * CURH(I,1) 
               END DO
            END IF
            IF (DOH) THEN
               DO I = 1, NPPX
                  CURH(I,1) = CURG(I,1) * CURH(I,1) 
               END DO
            END IF
C
C        (T,0,0)
C
         ELSE IF (U.EQ.0 .AND. V.EQ.0) THEN
            CALL ERIMKF(CURF(1,TUV),ARGE,TMIN1,PQX,PQXGT0,
     &                  CURF(1,M2T),CURF(1,M1T))
            IF (DOG) CALL ERIMKG(CURG(1,TUV),ARGE,TMIN1,PQX,PQXGT0,
     &                           CURF(1,M2T),CURF(1,M1T),
     &                           CURG(1,M2T),CURG(1,M1T))
            IF (DOH) CALL ERIMKH(CURH(1,TUV),ARGE,TMIN1,PQX,PQXGT0,
     &                           CURG(1,M2T),CURG(1,M1T),
     &                           CURH(1,M2T),CURH(1,M1T))
C
C        (T,U,0)
C
         ELSE IF (V.EQ.0) THEN
            CALL ERIMKF(CURF(1,TUV),ARGE,UMIN1,PQY,PQYGT0,
     &                  CURF(1,M2U),CURF(1,M1U))
            IF (DOG) CALL ERIMKG(CURG(1,TUV),ARGE,UMIN1,PQY,PQYGT0,
     &                           CURF(1,M2U),CURF(1,M1U),
     &                           CURG(1,M2U),CURG(1,M1U))
            IF (DOH) CALL ERIMKH(CURH(1,TUV),ARGE,UMIN1,PQY,PQYGT0,
     &                           CURG(1,M2U),CURG(1,M1U),
     &                           CURH(1,M2U),CURH(1,M1U))
C
C        (T,U,V)
C
         ELSE
            CALL ERIMKF(CURF(1,TUV),ARGE,VMIN1,PQZ,PQZGT0,
     &                  CURF(1,M2V),CURF(1,M1V))
            IF (DOG) CALL ERIMKG(CURG(1,TUV),ARGE,VMIN1,PQZ,PQZGT0,
     &                           CURF(1,M2V),CURF(1,M1V),
     &                           CURG(1,M2V),CURG(1,M1V))
            IF (DOH) CALL ERIMKH(CURH(1,TUV),ARGE,VMIN1,PQZ,PQZGT0,
     &                           CURG(1,M2V),CURG(1,M1V),
     &                           CURH(1,M2V),CURH(1,M1V))
         END IF
      ENDDO
      ENDDO
      ENDDO
C
C     Print section
C
      IF (IPRINT .GE. 10) THEN
         CALL AROUND('Output from ERIFGHTUV')
         DO ITUV = 1, NTUV
            WRITE(LUPRI,'(A5,I5)') ' TUV =', ITUV
            WRITE(LUPRI,'(3X,3(10X,A4))') 'CURF','CURG','CURH'
            DO I = 1, NPPX
               WRITE(LUPRI,'(I3,5F14.9)') I,
     &         CURF(I,ITUV),CURG(I,ITUV),CURH(I,ITUV)
            END DO
         END DO
      END IF
      RETURN
      END
C  /* Deck ericrx */
      SUBROUTINE ERICRX(R,S,X,T,A,B,C,PQGT0)
C     Written by Wim Klopper (University of Oslo, 29 October 2001).
#include <implicit.h>
      PARAMETER (D0 = 0D0)
#include <ericom.h>
      LOGICAL PQGT0
      DIMENSION A(NPPX), B(NPPX), C(NPPX), R(NPPX), X(NPPX)
C
      IF (T .GT. D0) THEN
         IF (PQGT0) THEN
            DO I = 1, NPPX
               R(I) = S * X(I) * (T*A(I) + B(I)*C(I))
            END DO
         ELSE
            DO I = 1, NPPX
               R(I) = S * X(I) * T*A(I)
            END DO
         END IF
      ELSE
         IF (PQGT0) THEN
            DO I = 1, NPPX
               R(I) = S * X(I) * B(I)*C(I)
            END DO
         ELSE
            DO I = 1, NPPX
               R(I) = D0
            END DO
         END IF
      END IF
      RETURN
      END 
C  /* Deck ericri */
      SUBROUTINE ERICRI(R,S,T,A,B,C,PQGT0)
C     Written by Wim Klopper (University of Oslo, 29 October 2001).
#include <implicit.h>
      PARAMETER (D0 = 0D0)
#include <ericom.h>
      LOGICAL PQGT0
      DIMENSION A(NPPX), B(NPPX), C(NPPX), R(NPPX)
C
      IF (T .GT. D0) THEN
         IF (PQGT0) THEN
            DO I = 1, NPPX
               R(I) = S * (T*A(I) + B(I)*C(I))
            END DO
         ELSE
            DO I = 1, NPPX
               R(I) = S * T*A(I)
            END DO
         END IF
      ELSE
         IF (PQGT0) THEN
            DO I = 1, NPPX
               R(I) = S * B(I)*C(I)
            END DO
         ELSE
            DO I = 1, NPPX
               R(I) = D0
            END DO
         END IF
      END IF
      RETURN
      END 
C  /* Deck ericrp */
      SUBROUTINE ERICRP(R,S,X,T,A,B,C,PQGT0)
C     Written by Wim Klopper (University of Oslo, 29 October 2001).
#include <implicit.h>
      PARAMETER (D0 = 0D0)
#include <ericom.h>
      LOGICAL PQGT0
      DIMENSION A(NPPX), B(NPPX), C(NPPX), R(NPPX), X(NPPX)
C
      IF (T .GT. D0) THEN
         IF (PQGT0) THEN
            DO I = 1, NPPX
               R(I) = R(I) + S * X(I) * (T*A(I) + B(I)*C(I))
            END DO
         ELSE
            DO I = 1, NPPX
               R(I) = R(I) + S * X(I) * T*A(I)
            END DO
         END IF
      ELSE
         IF (PQGT0) THEN
            DO I = 1, NPPX
               R(I) = R(I) + S * X(I) * B(I)*C(I)
            END DO
         END IF
      END IF
      RETURN
      END 
C --- end of eri2r12.F ---
