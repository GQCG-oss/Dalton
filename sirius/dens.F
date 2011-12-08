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
#ifdef UNDEF
/* Comdeck last_changes */
C modifications by jeppe sunday jan 17  '88  ( DSPDM3 ADDED )
C 880715-hjaaj: Lunar to Sirius orbital reordering.
#endif
C  /* Deck dispdm */
      SUBROUTINE DISPDM(RHOA,RHOB,RHO1,RHO2,NORB,K,L,
     *                  RHOW,ILTSOB,ISTLOB,L1,L2,IRA,IRB)
C
C INPUT:
C        RHOA(J,I) =         <L!EA(IJ)EA(KL)!R>
C                  +  SIGNLL * <L!EB(IJ)EB(KL)!R>
C                     SIGNLL = (-1)**(L1+L2)
C                   ( ALL IJ SO [IJ] .GE. [KL] )
C
C        RHOB(J,I) = <L!EB(IJ)EA(KL)!R>  ( ALL I,J )
C
C        RHO1(J,I) = <L!EA(IJ)+EB(JI)!R>
C
C        K AND L ARE FIXED AND [IJ] = (I-1)*NORB+J,
C        WHERE NORB IS NUMBER OF ORBITALS
C        /L) IS LEFT STATE AND /R) IS RIGHT STATE
C
C        RHOW(NORB,NORB) work space
C
C        ILTSOB(I) : Lunar to Sirius orbital reordering
C        ISTLOB(I) : Sirius to Lunar orbital reordering
C
C        L1 : FIRST  OPERATOR IS EA(IJ) + (-1)** L1 * EB(IJ)
C        L2 : SECOND OPERATOR IS EA(IJ) + (-1)** L2 * EB(IJ)
C OUTPUT:
C  DISTRIBUTE DENSITY MATRIX ELEMENTS OBTAINED IN DENSI1 INTO
C  TWOBODY DENSITY MATRIX (RHO2(IJKL)). RHOTYP ( TRANSFERRED IN
C  COMMON BLOCK CBRHOTYP ) DETERMINE THE TYPE OF PACKING YOU GET OF
C  RHO2
C
C  RHOTYP = 1 SYMMETRIC TWO-BODY DENSITY MATRIX FOR /L) = /R) WITH
C             TRIANGULAR PACKED DENSITIES
C             CONSTRAINTS ON INDICES [ I >= J ] >= [ K >= L ]
C
C  RHO2(IJKL)   =  [ (L/ E(IJ,KL) - DELTA(JK) E(IL) +
C                        E(IJ,LK) - DELTA(JL) E(IK) /R) ] /2
C               =  [ (L/ E(IJ,KL) - DELTA(JK) E(IL) +
C                        E(JI,KL) - DELTA(IK) E(JL) /R) ] /2
C             THE ONE ELECTRON TERMS ARE ADDED AFTER THE DENSI CALL
C
C  RHOTYP = 2 SYMMETRIC TWO-BODY TRANSITION DENSITY MATRIX FOR STATE
C             /L) AND /R) WITH SQUARE PACKED DISTRIBUTIONS
C             CONSTRAINT [IJ] >= [KL] IS REMOVED AFTER THE DENSI CALL
C
C  RHO2(I,J,K,L) =   (L/ E(IJ,KL) - DELTA(JK) E(IR)/R)
C             THE SYMMETRIZATION ARE CARRIED OUT AFTER THE DENSI CALL
C
C  RHOTYP = 3 SYMMETRIC TWO-BODY TRANSITION DENSITY MATRIX FOR STATE
C             /L) AND /R) WITH TRIANGULAR PACKED DENSITIES
C             CONSTRAINTS ON INDICES [ I >= J ] >= [ K >= L ]
C
C    RHO2 =   [ (L/ E(IJ,KL) - DELTA(JK) E(IL) +
C                   E(IJ,LK) - DELTA(JL) E(IK) /R)  +
C               (R/ E(IJ,KL) - DELTA(JK) E(IL) +
C                   E(IJ,LK) - DELTA(JL) E(IK) /L) ] / 2
C             THE ONE ELECTRON TERMS ARE ADDED AFTER THE DENSI CALL
C
C
#include <implicit.h>
C
      DIMENSION RHO1(NORB,NORB),RHOA(NORB,NORB),RHOB(NORB,NORB)
      DIMENSION RHO2(*)
      DIMENSION RHOW(NORB,NORB),ILTSOB(NORB),ISTLOB(NORB)
C
#include <cbrhotyp.h>
#include <priunit.h>
C
      GO TO (1,2,3),RHOTYP
      WRITE (LUPRI,'(/A,I8)')
     *   ' DISPDM, INCORRECT SPECIFICATION OF RHOTYP:',RHOTYP
      CALL QUIT('DISPDM,INCORRECT RHOTYP')
C
C
 1    CONTINUE
C     ... RHOTYP 1
      CALL DSPDM1(RHOA,RHOB,RHO1,RHO2,NORB,K,L,
     *            RHOW,ILTSOB,ISTLOB,L1,L2,IRA,IRB)
C
      GO TO 1000
C
C
 2    CONTINUE
C     ... RHOTYP 2
      CALL DSPDM2(RHOA,RHOB,RHO1,RHO2,NORB,K,L,
     *            RHOW,ILTSOB,ISTLOB,L1,L2,IRA,IRB)
C
      GO TO 1000
C
 3    CONTINUE
C     ... RHOTYP 3
      CALL DSPDM3(RHOA,RHOB,RHO1,RHO2,NORB,K,L,
     *            RHOW,ILTSOB,ISTLOB,L1,L2,IRA,IRB)
C
      GO TO 1000
 1000 CONTINUE
C
      RETURN
      END
C  /* Deck dspdm1 */
      SUBROUTINE DSPDM1(RHOA,RHOB,RHO1,RHO2,NORB,K,L,
     *                  RHOW,ILTSOB,ISTLOB,L1,L2,IRA,IRB)
C
C
C DISTRIBUTE DENSITY MATRIX ELEMENTS OBTAINED IN DENSI1 INTO
C SYMMETRIZED TWOBODY DENSITY MATRIX
C
C ( RHOTYP .EQ. 1 )
C
C
C... INPUT
C
C IRA = 1 :
C    RHOA(J,I) : <0!EA(IJ)EA(KL)!0> +(-1)**(L1+L2)*
C                <0!EB(IJ)EB(KL)!0>
C                FOR [IJ] .GE. [KL]
C
C IRB = 1 :
C    RHOB(J,I) : <0!EB(IJ)EA(KL)  !0>  FOR ALL IJ
C
C    RHO1(J,I) : <0!EA(IJ)+(-1)**(L1+L2)*EB(IJ)!0>
C
C    RHO2(IJ,KL): TWOBODY DENSITY MARTRIX UNDER CONSTRUCTION
C                PACKED WITH FOUR FOLD PERMUTATIONAL SYMMETRY
C                (but symmetrized eight fold)
C
C    NORB      : NUMBER OF ORBITALS
C
C    K,L       : Lunar INDICES FOR CURRENT BLOCK OF RHO2 CONSTRUCTED
C
C    RHOW(NORB,NORB) work space
C
C    ILTSOB(I) : Lunar to Sirius orbital reordering
C    ISTLOB(I) : Sirius to Lunar orbital reordering
C
C.. TASK
C
C    DISTRIBUTE ELEMENTS IN RHOA AND RHOB INTO RHO2 .
C
C    RHO2 IS DEFINED AS
C
C    RHO2(IJ,KL) = <0!(E1(IJ)(E2(KL)+E2(LK))! 0>
C
C               =
C                  (-1)**L1*     <0!(EB(IJ)+EB(JI)) EA(KL) !0>
C               +  (-1)**L2*     <0! EB(KL) (EA(IJ)+E(JI)) !0>
C
C               +                <0! EA(IJ)(EA(KL)+EA(LK)) !0>
C               +  (-1)**(L1+L2)*<0! EB(IJ)(EB(KL)+EB(LK)) !0>
C
#include <implicit.h>
      DIMENSION RHO1(NORB,NORB),RHOA(NORB,NORB),RHOB(NORB,NORB)
      DIMENSION RHO2(*)
      DIMENSION RHOW(NORB,NORB),ILTSOB(NORB),ISTLOB(NORB)
C
C Used from common blocks:
C  CBREOR : SLREOR
C
#include <cbreor.h>
C
      NNORB = (NORB*NORB + NORB)/2
      N2ORB = NORB*NORB
      KS = ILTSOB(K)
      LS = ILTSOB(L)
      KK = MAX(KS,LS)
      LL = MIN(KS,LS)
      KL = KK*(KK-1)/2 + LL
      IF (IRB .EQ. 1) THEN
C
C     Repack RHOB to Sirius orbital order
C
      IF (SLREOR) THEN
         CALL REORMT(RHOB,RHOW,NORB,NORB,ISTLOB,ISTLOB)
C        CALL REORMT(Ain ,Aout,Nrow,Ncol,Irow  ,Icol  )
      ELSE
         CALL DCOPY(N2ORB,RHOB,1,RHOW,1)
      END IF
C
C... TERM 1 : <0!(EB(IJ)+EB(JI)) EA(KL) !0>
C
      IF( KS .GE. LS) THEN
C     ... Otherwise some terms will be included twice
        SIGN1 = (-1.0D0) ** L1
        IF ( SIGN1 .NE. 1.0D0 ) CALL DSCAL(N2ORB,SIGN1,RHOW,1)
        DO 110 I = KK, NORB
          IF ( I .EQ. KK ) THEN
            JMIN = LL
          ELSE
            JMIN = 1
          END IF
          IKL = I*(I-1)/2 + NNORB*(KL-1)
          DO 105 J = JMIN, I
            RHO2(IKL+J) = RHO2(IKL+J) + RHOW(I,J) + RHOW(J,I)
  105     CONTINUE
  110   CONTINUE
      ELSE
        SIGN1 = 1.0D0
      END IF
C
C... TERM 2 :    <0! EB(IJ) (EA(KL)+E(LK)) !0> TO RHO2(KL,IJ)
C
      SIGN2 = (-1.0D0) ** L2
      FACT2 = SIGN1 * SIGN2
C     ... reverse effect of SIGN1
      IF ( K .EQ. L ) FACT2 = 2.D0*FACT2
      IF ( FACT2 .NE. 1.D0 ) CALL DSCAL(N2ORB,FACT2,RHOW,1)
      DO 185 I = 1,KK
        IF( I .EQ. KK) THEN
          JMAX = LL
        ELSE
          JMAX = I
        END IF
C
        IKL = I*(I-1)/2 + NNORB*(KL-1)
        DO 180 J = 1, JMAX
          RHO2(IKL+J) = RHO2(IKL+J) + RHOW(J,I)
  180   CONTINUE
  185 CONTINUE
      END IF
C
C... TERM 3 :  <0! ES(IJ)(ES(KL)+ES(LK)) !0>
C
      IF (IRA .EQ. 1) THEN
         KR = MAX(K,L)
         LR = MIN(K,L)
      IF (.NOT.SLREOR) THEN
         IF (K .EQ. L) CALL DSCAL(N2ORB,2.0D0,RHOA,1)
         DO 295 I = KR,NORB
            IKL = I*(I-1)/2 + NNORB*(KL-1)
            IF (I .EQ. KR) THEN
               JMIN = LR
            ELSE
               JMIN = 1
            END IF
            DO 290 J = JMIN,I
               RHO2(IKL+J) = RHO2(IKL+J) + RHOA(J,I)
  290       CONTINUE
  295    CONTINUE
      ELSE
         CALL DZERO(RHOW,N2ORB)
         DO 1195 I = KR,NORB
            IS = ILTSOB(I)
            IF (I .EQ. KR) THEN
               JMIN = LR
            ELSE
               JMIN = 1
            END IF
            DO 1190 J = JMIN,I
               JS = ILTSOB(J)
               IF (JS .LE. IS) THEN
                  RHOW(JS,IS) = RHOA(J,I)
               ELSE
                  RHOW(IS,JS) = RHOA(J,I)
               END IF
 1190       CONTINUE
 1195    CONTINUE
         IF (K .EQ. L) CALL DSCAL(N2ORB,2.0D0,RHOW,1)
C
         DO 195 I = 1,NORB
            IKL = I*(I-1)/2 + NNORB*(KL-1)
            DO 190 J = 1,I
               RHO2(IKL+J) = RHO2(IKL+J) + RHOW(J,I)
  190       CONTINUE
  195    CONTINUE
      END IF
      END IF
C
      RETURN
C     ... end of dspdm1.
      END
C  /* Deck dspdm2 */
      SUBROUTINE DSPDM2(RHOA,RHOB,RHO1,RHO2,NORB,K,L,
     *                  RHOW,ILTSOB,ISTLOB,L1,L2,IRA,IRB)
C
C DISTRIBUTE DENSITY MATRIX ELEMENTS OBTAINED IN DENSI1 INTO
C SQUARE PACKED TRANSITION TWOBODY DENSITY MATRIX
C
C ( RHOTYP .EQ. 2 )
C
C... INPUT
C
C IRA = 1 :
C    RHOA(J,I) : <0!EA(IJ)EA(KL)!0> +(-1)**(L1+L2)*
C                <0!EB(IJ)EB(KL)!0>
C                FOR [IJ] .GE. [KL]
C
C IRB = 1 :
C    RHOB(J,I) : <L!EB(IJ)EA(KL)!R>  FOR ALL IJ
C
C    RHO1(J,I) : <L!EA(IJ)+EB(IJ)!R>
C
C    RHO2(I,J,K,L): TWOBODY DENSITY MARTRIX UNDER CONSTRUCTION
C                   WITH ELEMENTS [IJ] >= [KL] CALCULATED CORRECTLY
C IF L1 + L2 = 1 THE MATRIX IS CALCULATED FOR ALL VALUES OF
C    I,J,K,L
C
C    NORB      : NUMBER OF ORBITALS
C
C    K,L       : Lunar INDICES FOR CURRENT BLOCK OF RHO2 CONSTRUCTED
C
C    RHOW(NORB,NORB) work space
C
C    ILTSOB(I) : Lunar to Sirius orbital reordering
C    ISTLOB(I) : Sirius to Lunar orbital reordering
C
C.. TASK
C
C    DISTRIBUTE ELEMENTS IN RHOA AND RHOB INTO RHO2 AND ADD ONE ELECTRON
C    MATRIX ELEMENTS
C
C    RHO2 IS DEFINED AS
C
C    RHO2(IJKL) = <L! E(IJ)E(KL)! R> - DELTA(J,K) <L! E(IL) !R>
C               =
C                  <L! EB(IJ) EA(KL) !R>
C               +  <L! EB(KL) EA(IJ) !R>
C
C               +  <L! EA(IJ) EA(KL) !R>
C               +  <L! EB(IJ) EB(KL) !R>
C
C               -   DELTA(J,K) <L! E(IL) !R>
C
C
#include <implicit.h>
C
      DIMENSION RHO1(NORB,NORB),RHOA(NORB,NORB),RHOB(NORB,NORB)
      DIMENSION RHO2(NORB,NORB,NORB,NORB)
      DIMENSION RHOW(NORB,NORB),ILTSOB(NORB),ISTLOB(NORB)
C
C
      SIGN1 = (-1.0D0) ** L1
      SIGN2 = (-1.0D0) ** L2
      KS = ILTSOB(K)
      LS = ILTSOB(L)
      IF ( L1 + L2 .NE. 1 ) THEN
C     I .eq. K:
      DO 205 J = L,NORB
         JS = ILTSOB(J)
         RHO2(KS,JS,KS,LS) = RHO2(KS,JS,KS,LS)
     *                     + RHOA(J,K) + RHOB(J,K)*SIGN1
  205 CONTINUE
C     I .gt. K:
      DO 220 I = K+1,NORB
         IS = ILTSOB(I)
         DO 210 J = 1,NORB
            JS = ILTSOB(J)
            RHO2(IS,JS,KS,LS) = RHO2(IS,JS,KS,LS)
     *                        + RHOA(J,I) + RHOB(J,I)*SIGN1
C           ... EA(ij) EA(KL) + EB(ij) EB(KL)   -- RHOA
C           ... EB(ij) EA(KL)                   -- RHOB
 210     CONTINUE
 220  CONTINUE
C
      IF (IRB .EQ. 1) THEN
C     I .eq. K:
      DO 405 J = 1,L
         JS = ILTSOB(J)
         RHO2(KS,LS,KS,JS) = RHO2(KS,LS,KS,JS) + RHOB(J,K)*SIGN2
  405 CONTINUE
C     I .lt. K:
      DO 420 I = 1,K-1
         IS = ILTSOB(I)
         DO 410 J = 1,NORB
            JS = ILTSOB(J)
            RHO2(KS,LS,IS,JS) = RHO2(KS,LS,IS,JS) + RHOB(J,I)*SIGN2
C           ... EA(KL) EB(ij)                     -- RHOB
 410     CONTINUE
 420  CONTINUE
      END IF
C
C        ... only ij .ge. kl
C
      IF (K .GE. L) THEN
         IST = K
      ELSE
         IST = K + 1
      END IF
      DO 620 I = IST,NORB
         IS = ILTSOB(I)
         RHO2(IS,KS,KS,LS) = RHO2(IS,KS,KS,LS) - RHO1(I,L)
C        ... - delta(jk) E(il)
 620  CONTINUE
      ELSE
C SPECIAL VERSION FOR ST AND TS MATRICES
C DUE TO THEIR LACK OF PARTICLE SYMMETRY
C
C
      IF (IRB .EQ. 1) THEN
      DO 1100 J = 1,NORB
        JS = ILTSOB(J)
        DO 1050 I = 1,NORB
          IS = ILTSOB(I)
          RHO2(IS,JS,KS,LS) = RHO2(IS,JS,KS,LS)
     &                      + SIGN1 * RHOB(J,I)
          RHO2(KS,LS,IS,JS) = RHO2(KS,LS,IS,JS)
     &                      + SIGN2 * RHOB(J,I)
 1050   CONTINUE
 1100 CONTINUE
      END IF
C
      DO 1200 I = K, NORB
C
        IF(I.EQ.K) THEN
          JMIN = L
        ELSE
          JMIN = 1
        END IF
        IS = ILTSOB(I)
        IF (IRA .EQ. 1) THEN
        DO 1150 J = JMIN,NORB
          JS = ILTSOB(J)
          RHO2(IS,JS,KS,LS) = RHO2(IS,JS,KS,LS)+RHOA(J,I)
C
          IF(.NOT. ( I.EQ.K.AND.J.EQ.L) ) THEN
            RHO2(KS,LS,IS,JS) = RHO2(KS,LS,IS,JS)+RHOA(J,I)
          END IF
C
 1150   CONTINUE
        END IF
        IF (JMIN .LE. K) THEN
          RHO2(IS,KS,KS,LS) = RHO2(IS,KS,KS,LS)-RHO1(I,L)
          IF(.NOT. ( I.EQ.K.AND.K.EQ.L) ) THEN
            RHO2(KS,LS,IS,KS) = RHO2(KS,LS,IS,KS)-RHO1(I,L)
          END IF
        END IF
 1200 CONTINUE
      END IF
      RETURN
C     ... end of dspdm2.
      END
C  /* Deck dspdm3 */
      SUBROUTINE DSPDM3(RHOA,RHOB,RHO1,RHO2,NORB,I1,I2,
     *                  RHOW,ILTSOB,ISTLOB,L1,L2,IRA,IRB)
C
C DISTRIBUTE DENSITY MATRIX ELEMENTS OBTAINED IN DENSI1 INTO
C SYMMETRIZED TRANSITION TWOBODY DENSITY MATRIX
C
C ( RHOTYP .EQ. 3 )
C
C
C... INPUT
C IRA = 1 :
C    RHOA(J,I) : <L!EA(IJ)EA(I1I2)!R>
C              + <L!EB(IJ)EB(I1I2)!R>  FOR [IJ] .GE. [I1I2]
C IRB = 1 :
C    RHOB(J,I) : <L!EB(IJ)EA(I1I2)  !R>  FOR ALL IJ
C
C    RHO1(J,I) : <L!EA(IJ)+EB(IJ)!R>
C
C    RHO2(IJ,KL): TWOBODY DENSITY MARTRIX UNDER CONSTRUCTION
C                PACKED WITH FOUR FOLD PERMUTATIONAL SYMMETRY
C                (but symmetrized eight fold)
C
C    NORB      : NUMBER OF ORBITALS
C
C    I1,I2     : Lunar INDICES FOR CURRENT BLOCK OF RHO2 CONSTRUCTED
C
C    RHOW(NORB,NORB) work space
C
C    ILTSOB(I) : Lunar to Sirius orbital reordering
C    ISTLOB(I) : Sirius to Lunar orbital reordering
C
C.. TASK
C
C    DISTRIBUTE ELEMENTS IN RHOA AND RHOB INTO RHO2 .
C
C    RHO2 IS DEFINED AS
C
C    RHO2(IJ,KL) = <L!E(IJ)(E(KL)+E(KL))! R>
C                + <R!E(IJ)(E(KL)+E(LK))! L>
C
C               =
C                  <L!(EB(IJ)+EB(JI)) EA(KL) !R>
C               +  <L!(EB(IJ)+EB(JI)) EA(LK) !R>
C               +  <L!(EB(KL)+EB(LK)) EA(IJ) !R>
C               +  <L!(EB(KL)+EB(LK)) EA(JI) !R>
C
C               +  <L!(ES(IJ)(ES(KL)+ES(LK)) !R>
C               +  <L!(ES(LK)+ES(KL))ES(JI)  !R>   ( S SUMS OVER A AND B)
C
#include <implicit.h>
C
      DIMENSION RHO1(NORB,NORB),RHOA(NORB,NORB),RHOB(NORB,NORB)
      DIMENSION RHO2(*)
      DIMENSION RHOW(NORB,NORB),ILTSOB(NORB),ISTLOB(NORB)
C
C Used from common blocks:
C  CBREOR : SLREOR
C
#include <cbreor.h>
C
      NNORB = (NORB*NORB + NORB)/2
      N2ORB = NORB*NORB
C
      SIGN1 = (-1.0D0) ** L1
      SIGN2 = (-1.0D0) ** L2
C     Calculate Sirius index for this distribution
C
      KS = ILTSOB(I1)
      LS = ILTSOB(I2)
      IF (KS .LT. LS) THEN
         II = KS
         KS = LS
         LS = II
      END IF
      KLS    = KS*(KS-1)/2 + LS
      IF(IRB.EQ.1) THEN
C
C
C.. TERMS FROM RHOB
C
C                  <L!(EB(IJ)+EB(JI)) EA(KL) !R>
C              OR  <L!(EB(IJ)+EB(JI)) EA(LK) !R>
C
C!
      K  = KS
      L  = LS
      KL = KLS
C
C
C     Repack RHOB to Sirius orbital order
C
      IF (.NOT. SLREOR) THEN
         CALL DCOPY(N2ORB,RHOB,1,RHOW,1)
      ELSE
         CALL REORMT(RHOB,RHOW,NORB,NORB,ISTLOB,ISTLOB)
C        CALL REORMT(Ain ,Aout,Nrow,Ncol,Irow  ,Icol  )
      END IF
C
      FACT1 = SIGN1
      IF (K .EQ. L) FACT1 = 2.0D0 * FACT1
      IF (FACT1 .NE. 1.0D0) CALL DSCAL(N2ORB,FACT1,RHOW,1)
C
      DO 110 I = K, NORB
        IF ( I .EQ. K ) THEN
          JMIN = L
        ELSE
          JMIN = 1
        END IF
        IKL = I*(I-1)/2 + (KL-1)*NNORB
        DO 105 J = JMIN, I
          RHO2(IKL+J) = RHO2(IKL+J) + RHOW(I,J) + RHOW(J,I)
  105   CONTINUE
  110 CONTINUE
C
C                  <L!(EB(KL)+EB(LK)) EA(IJ) !R>
C              OR  <L!(EB(KL)+EB(LK)) EA(JI) !R>
C
C
      I = K
      J = L
C
      IJ = I*(I-1)/2 + J
C
      FACT2 = SIGN1 * SIGN2
C     ... reverse effect of SIGN1
      IF (FACT2 .NE. 1.0D0) CALL DSCAL(N2ORB,FACT2,RHOW,1)
      DO 185 K = 1,I
        IF(K.EQ.I) THEN
         LMAX = J
        ELSE
         LMAX = K
        END IF
        IJK = K*(K-1)/2 + (IJ-1)*NNORB
        DO 180 L = 1,LMAX
          RHO2(IJK+L) = RHO2(IJK+L) + RHOW(K,L) + RHOW(L,K)
  180   CONTINUE
  185 CONTINUE
C
      END IF
C
C** TERMS FROM RHOA
C                 <L! ES(IJ)(ES(KL)+ES(LK)) !R>
C            +    <L!(ES(KL)+ES(LK))ES(JI)  !R>   ( S SUMS OVER A AND B)
C
C            =    <L!ES(IJ)ES(KL)!R>
C            +    <L!ES(IJ)ES(LK)!R>
C            +    <L!ES(LK)ES(JI)!R>
C            +    <L!ES(KL)ES(JI)!R>
C
C    SINCE THE RESTRICTION [IJ] .ge. [I1I2] IS PLACED ON THE ELEMENTS
C    IN RHO A
C    SOME GOOD OLD TYPE THINKING IS REQUIRED . AN IN DEPTH ANALYSIS
C    REVEALS THAT THE FIRST TWO TERMS NEVER CAUSES PROBLEMS
C    ([IJ] .GE. [KL] AND K.GE.L IMPLIES THAT [IJ] .GE. [LK] )
C    FOR TERM 3 AND 4 THREE CASES MUST BE CONSIDERED
C
C    CASE 1 : J .GE. K
C    CASE 2 : K .GT. J .GE. L
C    CASE 3 : L .GT. J
C
C    TERM 3 AND 4 ARE OBTAINED AS
C
C    CASE 1 :  RHOA(I,J)(I1I2 = LK ) + DELTA(J,K)*RHO1(I,L)
C                                    - DELTA(L,I)*RHO1(K,J)   ( TERM 3 )
C           +  RHOA(I,J)(I1I2 = KL ) + DELTA(J,L)*RHO1(I,K)
C                                    - DELTA(I,K)*RHO1(L,J)   ( TERM 4 )
C
C    CASE 2 :  RHOA(I,J) (I1I2 = LK )  ( TERM 3 )
C           +  RHOA(L,K) (I1I2 = JI )  ( TERM 4 )
C
C    CASE 3 :  RHOA(K,L) (I1I2 = JI )  ( TERM 3 )
C           +  RHOA(L,K) (I1I2 = JI )  ( TERM 4 )
C
C
C
      IF( I1 .GE. I2 ) THEN
C
        IF( I1 .EQ. I2 ) THEN
          FACTOR = 2.0D0
        ELSE
          FACTOR = 1.0D0
        END IF
C
C..     TERM 1
C
        IF(IRA.EQ.1) THEN
        K = I1
        L = I2
        DO 310 I = K,NORB
          IS = ILTSOB(I)
          IF( I .EQ. K ) THEN
            JMIN = L
          ELSE
            JMIN = 1
          END IF
          DO 300 J = JMIN, I
            JS = ILTSOB(J)
            IF (IS .GE. JS) THEN
               IJS = IS*(IS-1)/2 + JS
            ELSE
               IJS = JS*(JS-1)/2 + IS
            END IF
            IJKL = IJS + (KLS-1)*NNORB
            RHO2(IJKL) = RHO2(IJKL) + RHOA(J,I)*FACTOR
  300     CONTINUE
  310   CONTINUE
        END IF
C
C.. TERM 4 , CASE 1
C
        K = I1
        L = I2
        IF(IRA.EQ.1) THEN
        DO 410 I = K, NORB
          IS = ILTSOB(I)
          DO 400 J = K, I
            JS = ILTSOB(J)
            IF (IS .GE. JS) THEN
               IJS = IS*(IS-1)/2 + JS
            ELSE
               IJS = JS*(JS-1)/2 + IS
            END IF
            IJKL = IJS + (KLS-1)*NNORB
            RHO2(IJKL) = RHO2(IJKL) +  RHOA(I,J)*FACTOR
  400     CONTINUE
  410   CONTINUE
        END IF
C
        IF( K .EQ. L ) THEN
            DO 420 I = K,NORB
               IS = ILTSOB(I)
               IF (IS .GE. KS) THEN
                  IKS = IS*(IS-1)/2 + KS
               ELSE
                  IKS = KS*(KS-1)/2 + IS
               END IF
               IKKK = IKS + (KLS-1)*NNORB
               RHO2(IKKK) = RHO2(IKKK) + RHO1(K,I)*FACTOR
  420       CONTINUE
        END IF
        I1S = ILTSOB(I1)
        KKS = I1S*(I1S+1)/2
        KKKL = KKS + (KLS-1)*NNORB
        RHO2(KKKL) = RHO2(KKKL) - RHO1(K,L)*FACTOR
      END IF
C
C
      IF( I1.LT.I2 ) THEN
C
C.. TERM 2
C
        K = I2
        L = I1
        IF(IRA.EQ.1) THEN
        DO 510 I = K , NORB
          IS = ILTSOB(I)
C
          IF( I .EQ. K ) THEN
            JMIN = L
          ELSE
            JMIN = 1
          END IF
C
          DO 500 J = JMIN, I
            JS = ILTSOB(J)
            IF (IS .GE. JS) THEN
               IJS = IS*(IS-1)/2 + JS
            ELSE
               IJS = JS*(JS-1)/2 + IS
            END IF
            IJKL = IJS + (KLS-1)*NNORB
            RHO2(IJKL) = RHO2(IJKL) + RHOA(J,I)
  500     CONTINUE
  510   CONTINUE
        END IF
C
C.. CASE 1 TERM 3
C
        K = I2
        L = I1

        DO 610 I = K, NORB
          IS = ILTSOB(I)
          IF(IRA.EQ.1) THEN
          DO 600 J = K, I
            JS = ILTSOB(J)
            IF (IS .GE. JS) THEN
               IJS = IS*(IS-1)/2 + JS
            ELSE
               IJS = JS*(JS-1)/2 + IS
            END IF
            IJKL = IJS + NNORB*(KLS-1)
            RHO2(IJKL) = RHO2(IJKL) + RHOA(I,J)
  600     CONTINUE
          END IF
C
          I2S = ILTSOB(I2)
          IF (IS .GE. I2S) THEN
             IKS = IS*(IS-1)/2 + I2S
          ELSE
             IKS = I2S*(I2S-1)/2 + IS
          END IF
          IKKL = IKS + (KLS-1)*NNORB
          RHO2(IKKL) = RHO2(IKKL) + RHO1(L,I)
C
  610   CONTINUE
        IF( K.EQ.L ) THEN
           KKKK = KLS + (KLS-1)*NNORB
           RHO2(KKKK) = RHO2(KKKK) - RHO1(K,K)
        END IF
C
C.. CASE 2 , TERM 3
C
        IF(IRA.EQ.1) THEN
        K = I2
        L = I1
        DO 710 I = K, NORB
          IS = ILTSOB(I)
          DO 700 J = L, K - 1
            JS = ILTSOB(J)
            IF (IS .GE. JS) THEN
               IJS = IS*(IS-1)/2 + JS
            ELSE
               IJS = JS*(JS-1)/2 + IS
            END IF
            IJKL = IJS + (KLS-1)*NNORB
            RHO2(IJKL) = RHO2(IJKL) + RHOA(I,J)
  700     CONTINUE
  710   CONTINUE
        END IF
C
C.. CASE 2 , TERM 4
C
        IF(IRA.EQ.1) THEN
        J = I1
        I = I2
        DO 810 K = J + 1, I
          IS = ILTSOB(K)
          DO 800 L = 1 , J
            JS = ILTSOB(L)
            IF (IS .GE. JS) THEN
               IJS = IS*(IS-1)/2 + JS
            ELSE
               IJS = JS*(JS-1)/2 + IS
            END IF
            IJKL = IJS + (KLS-1)*NNORB
            RHO2(IJKL) = RHO2(IJKL) + RHOA(L,K)
  800     CONTINUE
  810   CONTINUE
        END IF
C
C.. CASE 3 TERM 3 AND 4
C
        IF(IRA.EQ.1) THEN
        J = I1
        I = I2
        DO 1010 K = J+1, I-1
           IS = ILTSOB(K)
           DO 1000 L = J+1, K
              JS = ILTSOB(L)
              IF (IS .GE. JS) THEN
                 IJS = IS*(IS-1)/2 + JS
              ELSE
                 IJS = JS*(JS-1)/2 + IS
              END IF
              IJKL = IJS + (KLS-1)*NNORB
              RHO2(IJKL) = RHO2(IJKL) + RHOA(L,K) + RHOA(K,L)
 1000      CONTINUE
 1010   CONTINUE
        END IF
C
      END IF
C
      RETURN
C     ... end of dspdm3.
      END
C  /* Deck dspdms */
#if defined (VAR_DSPDMS)
      SUBROUTINE DSPDMS(RHOA,RHOB,RHO1,RHO2,NORB,K,L,
     *                  RHOW,ILTSOB,ISTLOB,L1,L2)
C
C DISTRIBUTE DENSITY MATRIX ELEMENTS OBTAINED IN DENSI1 INTO
C SQUARE PACKED TRANSITION TWOBODY DENSITY MATRIX
C SAVED VERSION OF DSPDM2 BEFORE CHANGES FOR CORRECT TS
C
C ( RHOTYP .EQ. 2 )
C
C... INPUT
C
C    RHOA(J,I) : <L!EA(IJ)EA(KL)!R>
C              + <L!EB(IJ)EB(KL)!R>  FOR [IJ] .GE. [KL]
C
C    RHOB(J,I) : <L!EB(IJ)EA(KL)!R>  FOR ALL IJ
C
C    RHO1(J,I) : <L!EA(IJ)+EB(IJ)!R>
C
C    RHO2(I,J,K,L): TWOBODY DENSITY MARTRIX UNDER CONSTRUCTION
C                   WITH ELEMENTS [IJ] >= [KL] CALCULATED CORRECTLY
C
C IF L1 + L2 = 1 THE MATRIX IS CALCULATED FOR ALL VALUES OF
C    I,J,K,L
C    NORB      : NUMBER OF ORBITALS
C
C    K,L       : Lunar INDICES FOR CURRENT BLOCK OF RHO2 CONSTRUCTED
C
C    RHOW(NORB,NORB) work space
C
C    ILTSOB(I) : Lunar to Sirius orbital reordering
C    ISTLOB(I) : Sirius to Lunar orbital reordering
C
C.. TASK
C
C    DISTRIBUTE ELEMENTS IN RHOA AND RHOB INTO RHO2 AND ADD ONE ELECTRON
C    MATRIX ELEMENTS
C
C    RHO2 IS DEFINED AS
C
C    RHO2(IJKL) = <L! E(IJ)E(KL)! R> - DELTA(J,K) <L! E(IL) !R>
C               =
C                  <L! EB(IJ) EA(KL) !R>
C               +  <L! EB(KL) EA(IJ) !R>
C
C               +  <L! EA(IJ) EA(KL) !R>
C               +  <L! EB(IJ) EB(KL) !R>
C
C               -   DELTA(J,K) <L! E(IL) !R>
C
C
#include <implicit.h>
C
      DIMENSION RHO1(NORB,NORB),RHOA(NORB,NORB),RHOB(NORB,NORB)
      DIMENSION RHO2(NORB,NORB,NORB,NORB)
      DIMENSION RHOW(NORB,NORB),ILTSOB(NORB),ISTLOB(NORB)
C
#include <infpri.h>
C
C
      SIGN1 = (-1.0D0) ** L1
      SIGN2 = (-1.0D0) ** L2
      KS = ILTSOB(K)
      LS = ILTSOB(L)
C     I .eq. K:
      DO 205 J = L,NORB
         JS = ILTSOB(J)
         RHO2(KS,JS,KS,LS) = RHO2(KS,JS,KS,LS)
     *                     + RHOA(J,K) + RHOB(J,K)*SIGN1
  205 CONTINUE
C     I .gt. K:
      DO 220 I = K+1,NORB
         IS = ILTSOB(I)
         DO 210 J = 1,NORB
            JS = ILTSOB(J)
            RHO2(IS,JS,KS,LS) = RHO2(IS,JS,KS,LS)
     *                        + RHOA(J,I) + RHOB(J,I)*SIGN1
C           ... EA(ij) EA(KL) + EB(ij) EB(KL)   -- RHOA
C           ... EB(ij) EA(KL)                   -- RHOB
 210     CONTINUE
 220  CONTINUE
C
C     I .eq. K:
      DO 405 J = 1,L
         JS = ILTSOB(J)
         RHO2(KS,LS,KS,JS) = RHO2(KS,LS,KS,JS) + RHOB(J,K)*SIGN2
  405 CONTINUE
C     I .lt. K:
      DO 420 I = 1,K-1
         IS = ILTSOB(I)
         DO 410 J = 1,NORB
            JS = ILTSOB(J)
            RHO2(KS,LS,IS,JS) = RHO2(KS,LS,IS,JS) + RHOB(J,I)*SIGN2
C           ... EA(KL) EB(ij)                     -- RHOB
 410     CONTINUE
 420  CONTINUE
C
C        ... only ij .ge. kl
C
      IF (K .GE. L) THEN
         IST = K
      ELSE
         IST = K + 1
      END IF
      DO 620 I = IST,NORB
         IS = ILTSOB(I)
         RHO2(IS,KS,KS,LS) = RHO2(IS,KS,KS,LS) - RHO1(I,L)
C        ... - delta(jk) E(il)
 620  CONTINUE
C
      RETURN
C     ... end of dispd2.
      END
#endif
C  /* Deck r1tr21 */
      SUBROUTINE R1TR21(RHO1,RHO2,NORB,ILTSOB,ISTLOB)
C
C DISTRIBUTE TERMS FORM ONEBODY DENSITY MATRIX TO TWOBODY
C DENSITY MATRIX
C
C MODIFIED TO ALLOW FOR REORDERING OF ORBITALS
C Modified 890907-hjaaj for RHO2(NNORB,NNORB)
C    instead of packed RHO2.
C
#include <implicit.h>
      DIMENSION RHO1(NORB,NORB),RHO2(*)
      DIMENSION ILTSOB(NORB),ISTLOB(NORB)
C
      NNORB = (NORB*NORB + NORB)/2
C
C* CONTRIBUTION TO RHO2(IJ JL) IS (-  RHO1(L,I) )
C
      DO 100 I = 1, NORB
         IS = ILTSOB(I)
         DO 90 J = 1, I
            JS = ILTSOB(J)
            IF(IS .GE. JS ) THEN
               IJS = IS*(IS-1)/2 + JS
            ELSE
               IJS = JS*(JS-1)/2 + IS
            END IF
C
            DO 80 L = 1, J
               LS = ILTSOB(L)
               IF(JS.GE.LS) THEN
                 JLS = JS*(JS-1)/2+LS
               ELSE
                 JLS = LS*(LS-1)/2+JS
               END IF
C
               IJJL = JLS + (IJS-1)*NNORB
C
               RHO2( IJJL ) = RHO2( IJJL ) - RHO1(LS,IS)
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
C
C** CONTRIBUTION TO RHO2(IJ KJ ) IS  ( - RHO1(K,I))
C
      DO 200 I = 1, NORB
         IS = ILTSOB(I)
         DO 190 J = 1,I
            JS = ILTSOB(J)
            IF( IS .GE. JS ) THEN
              IJS = IS*(IS-1)/2 + JS
            ELSE
              IJS = JS*(JS-1)/2 + IS
            END IF
C
            DO 180 K = J,I
               KS = ILTSOB(K)
               IF(KS.GT.JS) THEN
                  KJS = KS*(KS-1)/2+JS
               ELSE
                  KJS = JS*(JS-1)/2+KS
               END IF
C
               KJIJ = KJS + (IJS-1)*NNORB
C
               RHO2( KJIJ ) = RHO2( KJIJ ) - RHO1(KS,IS)
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
C
C **  Symmetrize RHO2
C
      DO 400 KL = 1,NNORB
         DO 300 IJ = 1,KL-1
            IJKL = IJ + (KL-1)*NNORB
            KLIJ = KL + (IJ-1)*NNORB
            VALUE = RHO2(IJKL) + RHO2(KLIJ)
            RHO2(IJKL) = VALUE
            RHO2(KLIJ) = VALUE
  300    CONTINUE
  400 CONTINUE
C
      RETURN
      END
C  /* Deck r1tr22 */
      SUBROUTINE R1TR22(RHO1,RHO2,NORB)
C
C REMOVE CONSTRAINTS ON I,J,K,L IN RHO2
C NOTE: RHO1 was subtracted from RHO2 already in DSPDM2.
C
C INPUT:
C   RHO2(I,J,K,L) FOR [IJ] >= [KL] ; (I-1)*NORB + J >= (K-1)*NORB + L
C
C OUTPUT:
C   RHO2(I,J,K,L) FOR ALL I,J,K,L
C
#include <implicit.h>
      DIMENSION RHO1(NORB,*),RHO2(NORB,NORB,NORB,*)
C
      DO 400 K = 1,NORB
         DO 300 L = 1,NORB
            DO 200 I = 1,K-1
C
C              I<K AND J FREE
C
               DO 100 J = 1,NORB
                  VALUE = RHO2(I,J,K,L) + RHO2(K,L,I,J)
                  RHO2(I,J,K,L) = VALUE
                  RHO2(K,L,I,J) = VALUE
  100          CONTINUE
  200       CONTINUE
C
C           I=K AND J<L
C
            DO 210 J = 1,L-1
               VALUE = RHO2(K,J,K,L) + RHO2(K,L,K,J)
               RHO2(K,J,K,L) = VALUE
               RHO2(K,L,K,J) = VALUE
  210       CONTINUE
  300    CONTINUE
  400 CONTINUE
C
      RETURN
      END
C  /* Deck makdm */
      SUBROUTINE MAKDM(CREF,RHO1,RHO2,XNDXCI,WORK,LFREE)
#if defined (VAR_PACKPV)
#include <implicit.h>
      DIMENSION CREF(*),RHO1(*),RHO2(*),XNDXCI(*),WORK(*)
#include <inforb.h>
      KRHO2 = 1
      KWRK  = KRHO2 + NNASHX*NNASHX
      LWRK  = LFREE - KWRK
      CALL MAKD2(CREF,RHO1,WORK(KRHO2),XNDXCI,WORK(KWRK),LWRK)
      CALL DSITSP(NNASHX,WORK(KRHO2),RHO2)
      CALL DSCAL(NNASHY,0.5D0,RHO2,1)
C     ... old routines used PV = 0.5 * <| Eij Ekl - DELjk Eil |>
      RETURN
      END
      SUBROUTINE MAKD2(CREF,RHO1,RHO2,XNDXCI,WORK,LFREE)
#endif
C
C     CONSTRUCT ONE (RHO1) AND SYMMETRIC TWO (RHO2) PARTICLE DENSITY
C     MATRICES.
C
C  OUTPUT:
C
C  RHO1(KL)  = (0/ E(KL) /0)   K=>L
C
C
C  RHO2(IJKL): SYMMETRIC TWO-BODY DENSITY MATRIX FOR /L) = /R) PACKED
C              WITH FOUR FOLD SYMMETRY
C              CONSTRAINTS ON INDICES [ I >= J ] ; [ K >= L ]
C
C  RHO2(IJKL)= 0.5 * (0/ E(IJ,KL) - DELTA(JK) E(IL)
C                      + E(IJ,LK) - DELTA(JL) E(IK) /0)
C
#include <implicit.h>
C
      DIMENSION CREF(*),RHO1(*),RHO2(*),XNDXCI(*),WORK(*)
      PARAMETER ( D0 = 0.0 D0 )
C
C Used from common blocks:
C   INFINP : LSYM,FLAG(*)
C   INFORB : NASHT,N2ASHX,NNASHX
C   INFVAR : NCONF
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
C
      INTEGER RHOTYP
      LOGICAL USEPSM, CSFEXP, TDM, NORHO2
C
      RHOTYP = 1
      ILSYM  = LSYM
      IRSYM  = LSYM
      ISPIN1 = 0
      ISPIN2 = 0
      IF (FLAG(66)) THEN
         USEPSM = .FALSE.
      ELSE
         USEPSM = .TRUE.
      END IF
      TDM    = .FALSE.
      NORHO2 = .FALSE.
C
C
      CSFEXP = .NOT.FLAG(27)
C
      KRHO1  = 1
      KFREE  = KRHO1  + N2ASHX
      CALL SETVEC(WORK(KRHO1),D0,N2ASHX)
      CALL SETVEC(RHO2,D0,NNASHX*NNASHX)
      CALL DENSID(ILSYM,IRSYM,NCONF,NCONF,CREF,CREF,WORK(KRHO1),RHO2,
     *            RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     *            XNDXCI,WORK,KFREE,LFREE)
C     CALL DENSID(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR,RHO1,RHO2,
C    *            RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
C    *            XNDXCI,WORK,KFREE,LFREE)
C
      CALL DSITSP(NASHT,WORK(KRHO1),RHO1)
C     ... pack RHO1 as triangular matrix
C
      RETURN
C     ... end of MAKDM
      END
C  /* Deck makdv */
      SUBROUTINE MAKDV(CREF,RHO1,XNDXCI,WORK,LFREE)
C
C     CONSTRUCT ONE (RHO1) PARTICLE DENSITY MATRIX
C
C  OUTPUT:
C
C  RHO1(KL)  = (0/ E(KL) /0)   K>L
C
C
#include <implicit.h>
C
      DIMENSION CREF(*),RHO1(*),XNDXCI(*),WORK(*)
C
      PARAMETER ( D0 = 0.0 D0 )
#include <dummy.h>
C
C Used from common blocks:
C   INFINP : LSYM,FLAG(*)
C   INFORB : NASHT,N2ASHX,NNASHX
C   INFVAR : NCONF
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
C
      INTEGER RHOTYP
      LOGICAL USEPSM, CSFEXP, TDM, NORHO2
C
      RHOTYP = 1
      ISPIN1 = 0
      ISPIN2 = 0
      ILSYM  = LSYM
      IRSYM  = LSYM
      IF (FLAG(66)) THEN
         USEPSM = .FALSE.
      ELSE
         USEPSM = .TRUE.
      END IF
      TDM    = .FALSE.
      NORHO2 = .TRUE.
C
C
      CSFEXP = .NOT.FLAG(27)
C
      KRHO1  = 1
      KFREE  = KRHO1  + N2ASHX
      CALL SETVEC(WORK(KRHO1),D0,N2ASHX)
      CALL GETTIM(T1,W1)
      CALL DENSID(ILSYM,IRSYM,NCONF,NCONF,CREF,CREF,WORK(KRHO1),DUMMY,
     *            RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     *            XNDXCI,WORK,KFREE,LFREE)
C     CALL DENSID(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR,RHO1,RHO2,
C    *            RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
C    *            XNDXCI,WORK,KFREE,LFREE)
C
      CALL DSITSP(NASHT,WORK(KRHO1),RHO1)
C     ... pack RHO1 as triangular matrix
C
      RETURN
C
C     ... end of MAKDV
C
      END
C  /* Deck makdsv */
      SUBROUTINE MAKDSV(CREF,RHO1S,XNDXCI,WORK,LFREE)
C
C     CONSTRUCT (RHO1S) ONE PARTICLE SPIN DENSITY MATRIX
C
C  OUTPUT:
C
C     RHO1S(KL)  = (0/ EA(KL)-EB(KL) /0)   K>L
C
C
#include <implicit.h>
C
      DIMENSION CREF(*),RHO1S(*),XNDXCI(*),WORK(*)
C
      PARAMETER ( D0 = 0.0 D0 )
#include <dummy.h>
C
C Used from common blocks:
C   INFINP : LSYM,FLAG(*)
C   INFORB : NASHT,N2ASHX,NNASHX
C   INFVAR : NCONF
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
C
      INTEGER RHOTYP
      LOGICAL USEPSM, CSFEXP, TDM, NORHO2
C
      RHOTYP = 1
      ISPIN1 = 1
      ISPIN2 = 0
      ILSYM  = LSYM
      IRSYM  = LSYM
      IF (FLAG(66)) THEN
         USEPSM = .FALSE.
      ELSE
         USEPSM = .TRUE.
      END IF
      TDM    = .FALSE.
      NORHO2 = .TRUE.
C
C
      CSFEXP = .NOT.FLAG(27)
C
      KRHO1S  = 1
      KFREE  = KRHO1S  + N2ASHX
      CALL SETVEC(WORK(KRHO1S),D0,N2ASHX)
      CALL GETTIM(T1,W1)
      CALL DENSID(ILSYM,IRSYM,NCONF,NCONF,CREF,CREF,WORK(KRHO1S),DUMMY,
     *            RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     *            XNDXCI,WORK,KFREE,LFREE)
C     CALL DENSID(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR,RHO1S,RHO2,
C    *            RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
C    *            XNDXCI,WORK,KFREE,LFREE)
C
      CALL DSITSP(NASHT,WORK(KRHO1S),RHO1S)
C     ... pack RHO1S as triangular matrix
C
      RETURN
C
C     ... end of MAKDSV
C
      END
C  /* Deck maktdm */
      SUBROUTINE MAKTDM(NCSIM,BCVECS,CREF,RHO1,RHO2,XNDXCI,WORK,LWRK)
#if defined (VAR_PACKPV)
#include <implicit.h>
      DIMENSION BCVECS(*),CREF(*),RHO1(*)
      DIMENSION RHO2(NNASHY,*), XNDXCI(*),WORK(*)
#include <inforb.h>
      KRHO2 = 1
      KWRKD = KRHO2 + NCSIM*NNASHX*NNASHX
      LWRKD = LWRK + 1 - KWRK
      CALL MAKTD2(CREF,RHO1,WORK(KRHO2),XNDXCI,WORK(KWRKD),LWRKD)
      DO 100 I = 1,NCSIM
         JRHO2 = KRHO2 + (NCSIM-1)*NNASHX*NNASHX
         CALL DSITSP(NNASHX,WORK(JRHO2),RHO2(1,I))
  100 CONTINUE
      CALL DSCAL(NCSIM*NNASHY,0.5D0,RHO2,1)
C     ... old routines used PV = 0.5 * <| Eij Ekl - DELjk Eil |>
      RETURN
      END
      SUBROUTINE MAKTD2(NCSIM,BCVECS,CREF,RHO1,RHO2,XNDXCI,WORK,LWRK)
#endif
C
C     CONSTRUCT ONE (RHO1) AND SYMMETRIC TWO (RHO2) PARTICLE DENSITY
C     MATRICES.
C
C  OUTPUT:
C
C  RHO1(KL)   = (L/ E(KL) /R) + (R/ E(KL) /L)   K>=L
C
C  RHO2       : SYMMETRIC TWO-BODY TRANSITION DENSITY MATRIX FOR STATE
C               /L) AND /R) WITH TRIANGULAR PACKED DENSITIES
C               CONSTRAINTS ON INDICES [ I >= J ] ; [ K >= L ]
C
C  RHO2(IJKL) =   [ (L/ E(IJ,KL) - DELTA(JK) E(IL) +
C                       E(IJ,LK) - DELTA(JL) E(IK) /R) +
C                   (R/ E(IJ,KL) - DELTA(JK) E(IL) +
C                       E(IJ,LK) - DELTA(JL) E(IK) /L) ]
C
#include <implicit.h>
C
      DIMENSION BCVECS(NCONST,NCSIM),CREF(NCONRF)
      DIMENSION RHO1(NNASHX,NCSIM),RHO2(NNASHX*NNASHX,NCSIM)
      DIMENSION XNDXCI(*),WORK(*)
C
      PARAMETER ( D0 = 0.0 D0 )
C
C Used from common blocks:
C   INFINP : LSYM,FLAG(*)
C   INFORB : NASHT,N2ASHX,NNASHX
C   INFLIN : LSYMRF,LSYMST,NCONRF,NCONST
C   INFPRI : IPRDNS
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <inflin.h>
#include <infpri.h>
C
      INTEGER RHOTYP
      LOGICAL USEPSM, CSFEXP, TDM, NORHO2
C
      RHOTYP = 3
      ISPIN1 = 0
      ISPIN2 = 0
      IF (FLAG(66)) THEN
         USEPSM = .FALSE.
      ELSE
         USEPSM = .TRUE.
      END IF
      TDM    = .TRUE.
      NORHO2 = .FALSE.
C
C
      CSFEXP = .NOT.FLAG(27)
C
      KRHO1  = 1
      KWRKD  = KRHO1  + N2ASHX
      LWRKD  = LWRK + 1 - KWRKD
      DO 100 ICSIM=1,NCSIM
         KFREED = KWRKD
         CALL SETVEC(RHO2(1,ICSIM),D0,NNASHX*NNASHX)
         CALL SETVEC(WORK(KRHO1),D0,N2ASHX)
         CALL DENSID(LSYMRF,LSYMST,NCONRF,NCONST,CREF,BCVECS(1,ICSIM),
     *               WORK(KRHO1),RHO2(1,ICSIM),
     *               RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     *               XNDXCI,WORK,KFREED,LWRKD)
C        CALL DENSID(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR,RHO1,RHO2,
C    *               RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
C    *               XNDXCI,WORK,KFREE,LWRK)
C
         CALL DSITSP(NASHT,WORK(KRHO1),RHO1(1,ICSIM))
C     ... pack RHO1 as triangular matrix
         IF (IPRDNS .GE. 6) WRITE (LUPRI,'(/5X,A,2I5,2I10)')
     *      'MAKTDM : ICSIM, NCSIM, KFREE, LWRK =',
     *      ICSIM,NCSIM,KFREED,LWRK
  100 CONTINUE
C
C
        RETURN
C     ... end of MAKTDM
      END
C  /* Deck densid */
      SUBROUTINE DENSID(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR,RHO1,RHO2,
     *                  RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     *                  XNDXCI,WORK,KFRSAV,LFRSAV)
C
C   4-Sep-1989 hjaaj+pj
C   l.r. 980826-hjaaj: new CSF option for LSYM only
C
C MOTECC-90: The algorithms used in this module, DENSID,
C            are described in Chapter 8 Section D.4 of MOTECC-90
C            "Construction of Density matrices"
C
C   New interface routine for calculating density matrices
C   with determinant routines.
C   Density matrices are ADDED to RHO1 and RHO2.
C
C  RHO1 = (CR ! EA(IJ) ! CL )
C       + (-1)**(ISPIN1+ISPIN2)*(CR! EB(IJ) !CL )
C
C  RHO2 = (CR ! E1(IJ) E2(KL) ! CL ) - DELTA ( J,K ) * RHO1(I,L)
C
C  with
C       E1(IJ) = EA(IJ)+(-1)**ISPIN1 * EB(IJ)
C       E2(IJ) = EA(IJ)+(-1)**ISPIN2 * EB(IJ)
C
C
C   Input parameters:
C     RHOTYP  =1  rho2(ij,kl) density matrix
C             =2  rho2(i,j,k,l) used for response
C             =3  rho2(ij,kl) symmetrized transition density matrix
C     CSFEXP  =T  CL and/or CR in CSF basis (if of symmetry LSYM)
C             =F  CL and CR in determinant basis
C     USEPSM  =T  if MS .eq. 0 for both CL and CR, then use perm.sym.
C             =F  do not use permutation symmetry, even if MS .eq. 0
C     NORHO2  =T  only rho1, the one-particle density matrix
C             =F  rho1 and rho2, one- and two-particle density matrices
C     ISPIN1,ISPIN2 : see above
C     TDM     =T  CL(:) .ne. CR(:)
C             =F  CL(:) .eq. CR(:)
C
C   Internal:
C     PERMSM  =T  MS .eq. 0 for CL and CR
C             =F  MS .ne. 0 for CL and/or CR or do not use symmetry
C
C
#include <implicit.h>
      DIMENSION CL(NCLDIM), CR(NCRDIM), RHO1(*), RHO2(*)
      DIMENSION XNDXCI(*), WORK(*)
      LOGICAL   USEPSM,CSFEXP, PERMSM, NORHO2, TDM
      INTEGER   RHOTYP
C
#include <priunit.h>
C
C Used from common blocks:
C   INFINP : LSYM
C   CSFBAS : Pointers for CSF transformation
C   CIINFO : NDTASM
C   STRNUM : EQUAL
C   INFPRI : IPRDNS
C
#include <mxpdim.h>
#include <maxorb.h>
#include <infinp.h>
#include <csfbas.h>
#include <ciinfo.h>
#include <strnum.h>
C
#include <infpri.h>
C
C     work memory: WRK(KFRSAV:KFRSAV-1+LFRSAV)
C
      KFREE = KFRSAV
      LFREE = LFRSAV
      LWRK  = KFREE-1 + LFREE
C
C 980826-hjaaj: changed PERMSM test to the same as in CISIGD
C (including new test for IRSYM .eq. ILSYM, test today showed
C  PERMSM only works for totally symmetric density matrices)
C 990427-hjaaj: tests show PERMSM does not always work in this version,
C  we therefore never use PERMSM.
C
C old test:    IF (USEPSM) THEN
C new test:
Chj   IF (USEPSM .AND. ILSYM .EQ. IRSYM .AND.
Chj  &    ISPIN1 .EQ. 0 .AND. ISPIN2 .EQ. 0) THEN
Chj      PERMSM = EQUAL
Chj   ELSE
         PERMSM = .FALSE.
Chj   END IF
C
      IF (NORHO2) THEN
         IDMTYP = 1
      ELSE
         IDMTYP = 2
      END IF
C
C     any vectors with CSFs in this call ?
C     (only LSYM are with CSFs in this version. /980820-hjaaj)
C
      NCSF = 0
      IF (CSFEXP) THEN
         IF (ILSYM .EQ. LSYM) NCSF = NCSF + 1
         IF (IRSYM .EQ. LSYM) NCSF = NCSF + 1
      END IF
C
      IF (NCSF .GT. 0) THEN
         IERR = 0
         ILCSF = 0
         IF (ILSYM .EQ. LSYM) THEN
C           ABACUS and RESPONSE use CSF for singlet and det for triplet,
C           we find what is the case for this call by comparing
C           to NCSASM and NDTASM:
            IF (NCLDIM .EQ. NCSASM(ILSYM)) THEN
               ILCSF = 1
            ELSE IF (NCLDIM .NE. NDTASM(ILSYM)) THEN
               IERR = IERR + 1
            END IF
         ELSE
            IF (NCLDIM .NE. NDTASM(ILSYM)) IERR = IERR + 1
         END IF
         IRCSF = 0
         IF (IRSYM .EQ. LSYM) THEN
C           ABACUS and RESPONSE use CSF for singlet and det for triplet,
C           we find what is the case for this call by comparing
C           to NCSASM and NDTASM:
            IF (NCRDIM .EQ. NCSASM(IRSYM)) THEN
               IRCSF = 1
            ELSE IF (NCRDIM .NE. NDTASM(IRSYM)) THEN
               IERR = IERR + 1
            END IF
         ELSE
            IF (NCRDIM .NE. NDTASM(IRSYM)) IERR = IERR + 1
         END IF
         IF (IERR .GT. 0) THEN
            WRITE (LUPRI,*) 'CSF ERROR in DENSID, LSYM=',LSYM
            WRITE (LUPRI,*)
     &         'ILSYM, NCLDIM, NCSASM(ILSYM), NDTASM(ILSYM):',
     &          ILSYM, NCLDIM, NCSASM(ILSYM), NDTASM(ILSYM)
            WRITE (LUPRI,*)
     &         'IRSYM, NCRDIM, NCSASM(IRSYM), NDTASM(IRSYM):',
     &          IRSYM, NCRDIM, NCSASM(IRSYM), NDTASM(IRSYM)
            CALL QUIT('NCLDIM/NCRDIM disagree with NCSASM(:) in densid')
         END IF
         NLDET = NDTASM(ILSYM)
         NRDET = NDTASM(IRSYM)
         NCSF  = ILCSF + IRCSF
      ELSE
         NLDET = NCLDIM
         NRDET = NCRDIM
         IF (NLDET .NE. NDTASM(ILSYM) .OR.
     *       NRDET .NE. NDTASM(IRSYM)) THEN
            WRITE (LUPRI,*) 'ERROR in DENSID'
            WRITE (LUPRI,*) 'ILSYM, NLDET, NDTASM(ILSYM):',
     *                       ILSYM, NLDET, NDTASM(ILSYM)
            WRITE (LUPRI,*) 'IRSYM, NRDET, NDTASM(IRSYM):',
     *                       IRSYM, NRDET, NDTASM(IRSYM)
            CALL QUIT('NLDET/NRDET disagree with NDTASM(:) in densid')
         END IF
      END IF
C
C
C
      IF (NCSF .GT. 0) THEN
C
C        CSF - DET transformation of reference state
C
         NDET = MAX(NLDET,NRDET)
         CALL MEMADD(KLDET,NLDET, KFREE,2)
         IF (TDM) THEN
            CALL MEMADD(KRDET,NRDET, KFREE,2)
            ILVERV = 0
         ELSE
            KRDET = KLDET
            ILVERV = 1
         END IF
         KFREE1 = KFREE
         CALL MEMADD(KDET3,NDET, KFREE,2)
         IF (KFREE .GT. LWRK) CALL ERRWRK('DENSID CSF->DET',KFREE,LWRK)
         IF (ILCSF .EQ. 1) THEN
            CALL COPVEC(CL,WORK(KDET3),NCLDIM)
            CALL CSDTVC(WORK(KDET3),WORK(KLDET),1,XNDXCI(KDTOC),
     &                  XNDXCI(KICTS(1)),ILSYM,0,IPRDNS)
C           CSDTVC(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,IREFSM,ICOPY,NTEST)
         ELSE
            CALL COPVEC(CL,WORK(KLDET),NCLDIM)
         END IF
         IF (TDM) THEN
            IF (IRCSF .EQ. 1) THEN
               CALL COPVEC(CR,WORK(KDET3),NCRDIM)
               CALL CSDTVC(WORK(KDET3),WORK(KRDET),1,XNDXCI(KDTOC),
     &                     XNDXCI(KICTS(1)),IRSYM,0,IPRDNS)
            ELSE
               CALL COPVEC(CR,WORK(KRDET),NCRDIM)
            END IF
         END IF
         KFREE = KFREE1
         CALL DENSD2(ILSYM,IRSYM,NLDET,NRDET,WORK(KLDET),WORK(KRDET),
     *               RHO1,RHO2, RHOTYP,IDMTYP,ISPIN1,ISPIN2,PERMSM,
     *               XNDXCI,WORK,KFREE,LWRK,ILVERV)
      ELSE
C        ... determinant basis
         IF (TDM) THEN
            ILVERV = 0
         ELSE
            ILVERV = 1
         END IF
         CALL DENSD2(ILSYM,IRSYM,NLDET,NRDET,CL,CR,RHO1,RHO2,
     *               RHOTYP,IDMTYP,ISPIN1,ISPIN2,PERMSM,
     *               XNDXCI,WORK,KFREE,LWRK,ILVERV)
      END IF
C
      RETURN
C     ... end of DENSID
      END
C  /* Deck densd2 */
      SUBROUTINE DENSD2(ILSYM,IRSYM,NLDET,NRDET,CL,CR,RHO1,RHO2,
     *                  IRHOTP,IDMTYP,ISPIN1,ISPIN2,PERMSM,
     *                  XNDXCI,WORK,KFREE,LWRK,ILVERV)
C
C   Input parameters:
C     RHOTYP  =1  rho2(ij,kl) density matrix
C             =2  rho2(i,j,k,l) used for response
C             =3  rho2(ij,kl) symmetrized transition density matrix
C
C   NOTE: ILSYM = ICSYM  and KLOFF = KCOFF, KLDTAS = KCDTAS,
C                            KLCOOS = KICOOS, KLCSO = KICSO
C         IRSYM = IHCSYM and KROFF = KHCOFF, KRDTAS = KHDTAS,
C                            KRCOOS = KIHOOS, KRCSO = KIHCSO
C         This convention makes it possible to use same GETCIX call
C         for non-symmetric CISIGD and DENSID calls.  This convention
C         is the opposite of the one originally chosen by Jeppe Olsen
C         /890908-hjaaj
c
c Version with Turbo densi /oct 90 , Jeppe Olsen

C
#include <implicit.h>
      DIMENSION CL(NLDET), CR(NRDET), RHO1(*), RHO2(*)
      LOGICAL   PERMSM
      DIMENSION XNDXCI(*), WORK(LWRK)
C
      PARAMETER ( DP5 = 0.5D0 , D1 = 1.0D0 )
C
#include <priunit.h>
#include <maxash.h>
C
C Used from common blocks:
C   INFINP : ?
C   INFORB : MULD2H(8,8),NASHT,N2ASHX,NNASHX
C   INFPRI : IPRDNS
C
C   DETBAS : core allocation
C   STRNUM : EQUAL, NAEL,NBEL,NASTR,NBSTR,NAEXCI,...
C   CBRHOTYP : RHOTYP
C
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
#include <mxpdim.h>
#include <detbas.h>
#include <strnum.h>
#include <cbrhotyp.h>
#include <mxblk.h>
C
      RHOTYP = IRHOTP
C
      IF ( IDMTYP.NE.1 .AND. IDMTYP.NE.2) THEN
         WRITE (LUERR,'(//A/A,I10)')
     *      ' >>> DENSD2 fatal error, illegal density matrix code',
     *      '     IDMTYP =',IDMTYP
         CALL QUIT('DENSD2: Illegal IDMTYP')
      END IF
      IF ( IPRDNS.GT.15 ) THEN
         WRITE(LUPRI,'(/A)')' ***DENSD2 BEFORE CALLING DENSI'
         WRITE(LUPRI,'(/A,/I6,6I7,I8)')
     *     ' ILSYM  IRSYM  NLDET  NRDET RHOTYP IDMTYP  KFREE   LWRK :',
     *       ILSYM, IRSYM, NLDET, NRDET,RHOTYP,IDMTYP, KFREE,  LWRK
         WRITE(LUPRI,'(/A,/I6,I7,L7)') ' ISPIN1 ISPIN2 PERMSM :',
     *      ISPIN1,ISPIN2,PERMSM
      END IF
      IF ( IPRDNS.GT.120 ) THEN
         WRITE(LUPRI,'(/A)')' *DENSD2* LEFT REFERENCE VECTOR'
         CALL OUTPUT(CL,1,NLDET,1,1,NLDET,1,1,LUPRI)
         WRITE(LUPRI,'(/A)')' *DENSD2* RIGHT REFERENCE VECTOR'
         CALL OUTPUT(CR,1,NRDET,1,1,NRDET,1,1,LUPRI)
      END IF
      LUL = -1
      LUR = -1
c. Allocate local memory
      KFREEO = KFREE
      CALL MEMADD(KLLSOS,MAXSYM,KFREE,1)
      CALL MEMADD(KLRSOS,MAXSYM,KFREE,1)
      CALL MEMADD(KLSXST,MAXSYM ** 2, KFREE,1)
      CALL MEMADD(KLC2,MXVBLK,KFREE,2)
      CALL MEMADD(KLRHOW,NASHT ** 2, KFREE,2 )
      CALL MEMADD(KL,MXPST,KFREE,1)
      CALL MEMADD(KR,MXPST,KFREE,1)
      CALL MEMADD(KSIGNA,MXPST,KFREE,2)
      CALL MEMADD(KSGNAR,MXPST,KFREE,2)
      CALL MEMADD(KIIJSX,MXPST,KFREE,1)
      CALL MEMADD(KITSOX,MXPST,KFREE,1)
      CALL MEMADD(KVEC1,MXPST,KFREE,2)
      CALL MEMADD(KSCR1,MXPST,KFREE,2)
      CALL MEMADD(KNSXPT,MXPTP,KFREE,1)
      CALL MEMADD(KKOC,2*MXPST,KFREE,1)
      CALL MEMADD(KNFTPK,MXPTP,KFREE,1)
      CALL MEMADD(KRHO1A,NASHT**2,KFREE,2)
      CALL MEMADD(KRHO1B,NASHT**2,KFREE,2)
      CALL MEMADD(KRHO2A,NASHT**2,KFREE,2)
      CALL MEMADD(KRHO2B,NASHT**2,KFREE,2)
      IF (KFREE .GT. LWRK) CALL ERRWRK('DENSD2',KFREE,LWRK)
c. Symmetry array
      CALL ZSMOST(ILSYM,WORK(KLLSOS),MAXSYM)
      CALL ZSMOST(IRSYM,WORK(KLRSOS),MAXSYM)
      CALL SXSTSM(WORK(KLSXST),MAXSYM,MAXSYM)
c. And then , the density matrix
      CALL DENSIT(CL,CR,
     &            LUL,LUR,WORK(KLLSOS),WORK(KLRSOS),
     &            XNDXCI(KIOCOC),XNDXCI(KIOCOC),
     &            XNDXCI(KISSOA),XNDXCI(KNSSOA),
     &            XNDXCI(KISSOB),XNDXCI(KNSSOB),NOCTPA,NOCTPB,
     &            XNDXCI(KISSYM),WORK(KLSXST),MAXSYM,MAXSYM,
     &            XNDXCI(1),NASHT,RHO1,RHO2,WORK(KLC2) ,
     &            XNDXCI(KICREA),XNDXCI(KIANNI),
     &            IDMTYP,PERMSM,XNDXCI(KTPFOB),ISPIN1,ISPIN2,ILVERV,
     &            XNDXCI(KLTSOB),XNDXCI(KSTLOB),
     &            WORK(KLRHOW),WORK(KL),WORK(KR),WORK(KSIGNA),
     &            WORK(KSGNAR),WORK(KIIJSX),WORK(KITSOX),WORK(KVEC1),
     &            WORK(KSCR1),WORK(KNSXPT),WORK(KKOC),WORK(KNFTPK),
     &            WORK(KRHO1A),WORK(KRHO1B),WORK(KRHO2A),WORK(KRHO2B))
C** Symmetrize RHO1 for TDM
C
      IF (RHOTYP .EQ. 3) THEN
         CALL TRPAD(RHO1,D1,NASHT)
      END IF
c. Add contributions from one-electron density matrix
      IF ( IDMTYP .EQ. 2 ) THEN
         IF ( RHOTYP .EQ. 1 .OR. RHOTYP .EQ. 3 ) THEN
            CALL R1TR21(RHO1,RHO2,NASHT,
     &                  XNDXCI(KLTSOB),XNDXCI(KSTLOB))
            CALL DSCAL(NNASHX*NNASHX,DP5,RHO2,1)
C           ... RHO2(ij,kl) = 0.5 * (rho2(ij,kl) + rho2(ij,lk))
         ELSE IF ( RHOTYP .EQ.2 ) THEN
            IF ( (ISPIN1+ISPIN2) .NE. 1 ) THEN
               CALL R1TR22(RHO1,RHO2,NASHT)
            END IF
         ELSE
            WRITE(LUERR,*)' ERROR IN DENSD2'
            WRITE(LUERR,*)' RHOTYP out of range RHOTYP=',RHOTYP
            CALL QUIT('DENSD2 : RHOTYP out of range')
         END IF
      END IF
c
      IF ( IPRDNS.GT.15 ) THEN
         WRITE(LUPRI,'(/A)')' *** LEAVING DENSD2'
         WRITE(LUPRI,'(/A,/I6,6I7,I8)')
     *     ' ILSYM  IRSYM  NLDET  NRDET RHOTYP IDMTYP  KFREE   LWRK :',
     *       ILSYM, IRSYM, NLDET, NRDET,RHOTYP,IDMTYP, KFREE,  LWRK
         WRITE(LUPRI,'(/A,/I6,I7,L7)') ' ISPIN1 ISPIN2 PERMSM :',
     *      ISPIN1,ISPIN2,PERMSM
      END IF
      IF ( IPRDNS.GT.15 ) THEN
         CALL HEADER(
     &   'DENSD2: One-el. density matrix, active part, MO-basis',-1)
         CALL OUTPUT(RHO1,1,NASHT,1,NASHT,NASHT,NASHT,-11,LUPRI)
      ENDIF
C
      CALL MEMADD(KDUMMY,KFREEO-KFREE,KFREE,2 )
C     ... release space
C
      RETURN
C     ... end of DENSD2.
      END
C  /* Deck trpad */
      SUBROUTINE TRPAD(MAT,FACTOR,NDIM)
C
C  MAT(I,J) = MAT(I,J) + FACTOR*MAT(J,I)
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MAT(NDIM,NDIM)
C
C
      DO 100 J = 1, NDIM
        DO 90 I = J, NDIM
          MAT(I,J) =MAT(I,J) + FACTOR * MAT(J,I)
  90    CONTINUE
 100  CONTINUE
C
C
      IF( ABS(FACTOR) .NE. 1.0D0 ) THEN
        FAC2 = 1.0D0 - FACTOR**2
        DO 200 I = 1, NDIM
         DO 190 J = 1, I - 1
           MAT(J,I) = FACTOR*MAT(I,J ) + FAC2 * MAT(J,I)
 190     CONTINUE
 200    CONTINUE
      ELSE
        IF(FACTOR .EQ. 1.0D0) THEN
        DO 300 I = 1, NDIM
         DO 290 J = 1, I - 1
            MAT(J,I) = MAT(I,J )
 290     CONTINUE
 300    CONTINUE
      ELSE
        DO 400 I = 1, NDIM
         DO 390 J = 1, I - 1
            MAT(J,I) =-MAT(I,J )
 390     CONTINUE
 400    CONTINUE
      END IF
      END IF
      RETURN
      END
