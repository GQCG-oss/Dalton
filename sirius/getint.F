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
C  /* Deck getin2 */
      SUBROUTINE GETIN2( RIJKL, K1, L1, FIJKL, ILTSOB, WIJKL )
C
C  10-May-1987 hjaaj / 15-Jul-1988 hjaaj (Lunar-Sirius reorder)
C  Sep-1989 jo (distyp 4 and 5)
C  900801-hjaaj: SLREOR option
C  900914-olav vahtras: (distyp 6 and 8)
C  901109-hh and pj distyp 9
C  910416-hjaaj: transposed RIJKL; rewritten some of the code
C  920129-olav vahtras distyp 11-18
C  920228-olav vahtras distyp 19-26
C  960522-olav vahtras distyp 27-28
C
C  Get 2-electron active integrals distribution / K L ) in RIJKL
C  for Jeppe's determinant CI.
C  910416-hjaaj: NOTE that now RIJKL(I,J) = (I J / K L).
C
C  FIJKL is an auxiliary array, which may contain all or some of
C  the active integrals, or which may be used as buffer in reading
C  integrals from disk.
C
C  WIJKL is an auxiliary array, for reordering from Sirius order
C  to Lunar order.
C
C     DISTYP = 1 : SYMMETRIC TRIANGULAR PACKED FIJKL(NNASHX,NNASHX)
C                  RETURN DISTRIBUTION KL IN RIJKL
C
C            = 2 : HALF TRANSFORMED INTEGRALS IN FIJKL(N2ASHX,NNASHX)
C                  UNTRANSFORMED SECOND INDEX
C                  RETURN DISTRIBUTION KL IN RIJKL
C                  NOTE: RIJKL(I,J) = (I~1 J~2 / K~1 L~2)
C
C            = 3 : HALF TRANSFORMED INTEGRALS IN FIJKL(N2ASHX,NNASHX)
C                  UNTRANSFORMED SECOND INDEX
C                  RETURN DISTRIBUTION LK, WHICH IS TRANSPOSED INTEGRAL
C                  DISTRIBUTION FROM FIJKL, IN RIJKL
C                  NOTE: RIJKL(I,J) = (I~2 J~1 / K~2 L~1)
C
C            = 4 : HALF TRANSFORMED INTEGRALS IN FIJKL(N2ASHX,NNASHX)
C                  UNTRANSFORMED SECOND INDEX
C                  RETURN  FOR GIVEN KL IN
C                  RIJKL(I,J) = (k~ l~ | i j ) - ( i~ j~ | k l )
C
C            = 5 : HALF TRANSFORMED INTEGRALS IN FIJKL(N2ASHX,NNASHX)
C                  UNTRANSFORMED SECOND INDEX
C                  RETURN  FOR GIVEN  transposed KL IN
C                  RIJKL(I,J) = (k~ l~ | i j ) - ( i~ j~ | k l )
C
C            = 6 : SPIN-ORBIT INTEGRALS PACKED IN
C                  FIJKL(N2ASHX,N2ASHX) in the form (i j|k^l)
C                  RETURN FOR GIVEN KL IN RIJKL(I,J) =
C                  (i j|k~l) + (i~j|k l) = 3(i^j|k l) + 3(i j|k^l)
C                  (we have used (i~j|k l) = (i^j|k l) + 2(i j|k^l))
C
C            = 7 : --- not used ---
C
C            = 8 : SPIN-ORBIT INTEGRALS PACKED IN
C                  FIJKL(N2ASHX,N2ASHX) in the form (i j|k^l)
C                  RETURN FOR GIVEN KL IN RIJKL(I,J) =
C                  (i j|k~l) - (i~j|k l) = (i^j|k l) - (i j|k^l)
C                  (we have used (i~j|k l) = (i^j|k l) + 2(i j|k^l))
C
C            = 9 : SQUARE PACKED DISTRIBUTION  FIJKL OF FULLY
C                  TRANSFORMED INTEGRALS, NEEDED BY RSPEGR
C
C            = 10: SQUARE PACKED DISTRIBUTION  FIJKL OF FULLY
C                  TRANSFORMED INTEGRALS, NEEDED BY RSPEGR
C                  (TRANSPOSED)
C
C            = 11: SQUARE PACKED DISTRIBUTIONS OF DOUBLY TRANSFORMED
C                  INTEGRALS.
C                  FIJKL(1) IS ADDRESS TO (I~~J|KL)
C                  FIJKL(N2ASHX*N2ASHX+1) IS ADDRESS TO (I~J|K~L)
C                  RETURN SYMMETRIC COMBINATION
C                  W(I,J) = (I~~J|KL) + (I~J|K~L) + (K~~L|IJ) + (K~L|I~J)
C
C            = 12: SQUARE PACKED DISTRIBUTIONS OF DOUBLY TRANSFORMED
C                  INTEGRALS.
C                  FIJKL(1) IS ADDRESS TO (I~~J|KL)
C                  FIJKL(N2ASHX*N2ASHX+1) IS ADDRESS TO (I~J|K~L)
C                  RETURN SYMMETRIC COMBINATION
C                  W(I,J) = (J~~I|LK) + (J~I|L~K) + (L~~K|JI) + (L~K|J~I)
C
C            = 13: SQUARE PACKED DISTRIBUTIONS OF DOUBLY TRANSFORMED
C                  INTEGRALS.
C                  FIJKL(1) IS ADDRESS TO (I~~J|KL)
C                  FIJKL(N2ASHX*N2ASHX+1) IS ADDRESS TO (I~J|K~L)
C                  RETURN ANTISYMMETRIC COMBINATION
C                  W(I,J)=S*(I~~J|KL) + (K~~L|IJ) + S1*(I~J|K~L)  + S2*(K~L|I~J)
C
C            = 14: SQUARE PACKED DISTRIBUTIONS OF DOUBLY TRANSFORMED
C                  INTEGRALS.
C                  FIJKL(1) IS ADDRESS TO (I~~J|KL)
C                  FIJKL(N2ASHX*N2ASHX+1) IS ADDRESS TO (I~J|K~L)
C                  RETURN ANTISYMMETRIC COMBINATION
C                  W(I,J)=S*(J~~I|LK) + (L~~K|JI) + S1*(J~I|L~K)  + S2*(L~K|J~I)
C
C            = 15: SQUARE PACKED DISTRIBUTIONS OF HALF-TRANSFORMED
C                  INTEGRALS.
C                  FIJKL(1) IS ADDRESS TO (I~J|KL)
C                  RETURN SYMMETRIC COMBINATION
C                  W(I,J)=(I~J|KL) + (K~L|IJ)
C
C            = 16: SQUARE PACKED DISTRIBUTIONS OF HALF-TRANSFORMED
C                  INTEGRALS.
C                  FIJKL(1) IS ADDRESS TO (I~J|KL)
C                  RETURN SYMMETRIC COMBINATION
C                  W(I,J)=(J~I|LK) + (L~K|JI)
C
C            = 17: SQUARE PACKED DISTRIBUTIONS OF HALF-TRANSFORMED
C                  INTEGRALS.
C                  FIJKL(1) IS ADDRESS TO (I~J|KL)
C                  RETURN ANTISYMMETRIC COMBINATION
C                  W(I,J)=S1*(I~J|KL) + S2*(K~L|IJ)
C
C            = 18: SQUARE PACKED DISTRIBUTIONS OF HALF-TRANSFORMED
C                  INTEGRALS.
C                  FIJKL(1) IS ADDRESS TO (I~J|KL)
C                  RETURN ANTISYMMETRIC COMBINATION
C                  W(I,J)=S1*(J~I|LK) + S2*(L~K|JI)
C
C            = 19: SQUARE PACKED ONE-INDEX TRANSFORMED
C                  SPIN-ORBIT INTEGRALS
C                  FIJKL(1) IS ADDRESS TO (I~J|K^L)
C                  FIJKL(1+N2ASHX*N2ASHX) IS ADDRESS TO (IJ|K~^L)
C                  RETURN SYMMETRIC COMBINATION
C                  W(I,J)=3[ (I~J|K^L) + (K~L|I^J) + (IJ|K~^L) + (KL|I~^J)]
C
C            = 20: SQUARE PACKED ONE-INDEX TRANSFORMED
C                  SPIN-ORBIT INTEGRALS
C                  FIJKL(1) IS ADDRESS TO (I~J|K^L)
C                  FIJKL(1+N2ASHX*N2ASHX) IS ADDRESS TO (IJ|K~^L)
C                  RETURN SYMMETRIC COMBINATION
C                  W(I,J)=3[ (J~I|L^K) + (L~K|J^I) + (JI|L~^K) + (LK|J~^I)]
C
C            = 21: SQUARE PACKED ONE-INDEX TRANSFORMED
C                  SPIN-ORBIT INTEGRALS
C                  FIJKL(1) IS ADDRESS TO (I~J|K^L)
C                  FIJKL(1+N2ASHX*N2ASHX) IS ADDRESS TO (IJ|K~^L)
C                  RETURN ANTISYMMETRIC COMBINATION
C                  W(I,J)= (K~L|I^J) ~- (I~J|K^L) ~+ (KL|I~^J) - (IJ|K~^L)
C                  as 13 with S1=- (IOPSPI) S2=~ (ISPINV)
C
C            = 22: SQUARE PACKED ONE-INDEX TRANSFORMED
C                  SPIN-ORBIT INTEGRALS
C                  FIJKL(1) IS ADDRESS TO (I~J|K^L)
C                  FIJKL(1+N2ASHX*N2ASHX) IS ADDRESS TO (IJ|K~^L)
C                  RETURN ANTISYMMETRIC COMBINATION
C                  W(I,J)= (L~K|J^I) ~- (J~I|L^K) ~+ (LK|J~^I) - (JI|L~^K)
C                  as 14 with S1=- (IOPSPI) S2=~ (ISPINV)
C
C            = 23: see 6
C
C            = 24:  W(I,J) = 3(LK|J^I) + 3(JK|L^K)
C
C            = 25:  see 8
C
C            = 26:  W(I,J) = (LK|J^I) - (JI|L^K)
C                   S1=- S2=+
C
C            = 27: SQUARE PACKED DISTRIBUTIONS OF TRIPLY TRANSFORMED
C                  INTEGRALS.
C                  FIJKL(1) IS ADDRESS TO (I~~~J|KL)
C                  FIJKL(N2ASHX*N2ASHX+1)   IS ADDRESS TO (I~~J|K~L)
C                  FIJKL(2*N2ASHX*N2ASHX+1) IS ADDRESS TO (I~~J|K~L)
C                  FIJKL(3*N2ASHX*N2ASHX+1) IS ADDRESS TO (I~~J|K~L)
C                  RETURN SYMMETRIC COMBINATION
C                  W(I,J) = (I~~~J|KL) + (I~~J|K~L) + (I~~J|K~L) + (I~~J|K~L) 
C                         + (K~~~L|IJ) + (K~~L|I~J)  +(K~~L|I~J) + (K~~L|I~J)
C
C            = 28: SQUARE PACKED DISTRIBUTIONS OF TRIPLY TRANSFORMED
C                  INTEGRALS.
C                  FIJKL(1) IS ADDRESS TO (I~~~J|KL)
C                  FIJKL(N2ASHX*N2ASHX+1)   IS ADDRESS TO (I~~J|K~L)
C                  FIJKL(2*N2ASHX*N2ASHX+1) IS ADDRESS TO (I~~J|K~L)
C                  FIJKL(3*N2ASHX*N2ASHX+1) IS ADDRESS TO (I~~J|K~L)
C                  RETURN SYMMETRIC COMBINATION
C                  W(I,J) = (J~~~I|LK) + (J~~I|L~K) + (J~~I|L~K) + (J~~I|L~K) 
C                         + (L~~~K|JI) + (L~~K|J~I)  +(L~~K|J~I) + (L~~K|J~I)
C
#include <implicit.h>
      DIMENSION RIJKL(NASHT,NASHT), FIJKL(*), ILTSOB(*)
      DIMENSION WIJKL(NASHT,NASHT)
#include <iratdef.h>
      PARAMETER ( D1 = 1.0D0, DM1 = -1.0D0, D3 = 3.0D0 )
#include <priunit.h>
C
C Used from common blocks:
C
C  INFORB   : NASHT, N2ASHX
C  INFTAP   : LUH2AC
C  CBGETDIS : DISTYP, IADINT
C  CBREOR   : SLREOR
C  INFSPI   : ISGN1,ISGN2
C
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <cbgetdis.h>
#include <cbreor.h>
#include <infspi.h>
C
C
C     Get SIRIUS orbital numbers:
C
      K = ILTSOB(K1)
      L = ILTSOB(L1)
C
      IF (K .GE. L) THEN
         KL = K*(K-1)/2 + L
      ELSE
         KL = L*(L-1)/2 + K
      END IF
C
C
C
      IF (DISTYP .NE. 1 .AND. IADINT .GE. 0) THEN
         WRITE(LUERR,'(/A,I5)')
     &   ' GETINT ERROR: H2AC on disk is not implemented for DISTYP =',
     &   DISTYP
       CALL QUIT('GETINT: H2AC on disk not implemented for DISTYP.ne.1')
      END IF
      GO TO (1,2,3,4,5,6,99,8,9,10,11,12,13,14,15,16,17,18,
     *       19,20,13,14,23,24,25,26,27,28), DISTYP
C            19,20,21,22,23,24,25,26), DISTYP
 99      WRITE(LUERR,'(/A,I5)')
     *      ' GETINT: INCORRECT SPECIFICATION OF DISTYP, DISTYP:',DISTYP
         CALL QTRACE(LUERR)
         CALL QUIT('GETINT INCORRECT SPECIFICATION OF DISTYP')
C     ... DISTYP = 1
C     ... 2-electron integrals have full index symmetry
C         (ij/kl) = (ji/kl) = (kl/ij) = ...
C         (e.g. MC optimization or CI)
 1    CONTINUE
         IF (IADINT .GE. 0) THEN
            CALL READDX(LUH2AC,IADINT+KL,IRAT*NNASHX,FIJKL)
            KKL = 1
         ELSE
            KKL =  NNASHX*(KL-1) + 1
         END IF
         IF (SLREOR) THEN
            CALL DSPTSI(NASHT,FIJKL(KKL),WIJKL)
         ELSE
            CALL DSPTSI(NASHT,FIJKL(KKL),RIJKL)
         END IF
C        ... unpack FIJKL(**,KL) using CALL DSPTSI(N,ASP,ASI)
      GO TO 100
C
C        ... DISTYP = 2
 2    CONTINUE
         KLOFF = (L-1)*NASHT + K
         IF (SLREOR) THEN
            CALL DCOPY(NNASHX,FIJKL(KLOFF),N2ASHX,RIJKL,1)
            CALL DSPTSI(NASHT,RIJKL,WIJKL)
            CALL DAXPY(N2ASHX,D1,FIJKL(N2ASHX*(KL-1)+1),1,WIJKL,1)
         ELSE
            CALL DCOPY(NNASHX,FIJKL(KLOFF),N2ASHX,WIJKL,1)
            CALL DSPTSI(NASHT,WIJKL,RIJKL)
            CALL DAXPY(N2ASHX,D1,FIJKL(N2ASHX*(KL-1)+1),1,RIJKL,1)
         END IF
      GO TO 100
C
 3    CONTINUE
C        ... DISTYP = 3
         KLOFF = (K-1)*NASHT + L
         CALL DCOPY(NNASHX,FIJKL(KLOFF),N2ASHX,WIJKL,1)
         CALL DSPTSI(NASHT,WIJKL,RIJKL)
         CALL MTRSP(NASHT,NASHT,
     &              FIJKL(N2ASHX*(KL-1)+1),NASHT,WIJKL,NASHT)
         IF (SLREOR) THEN
            CALL DAXPY(N2ASHX,D1,RIJKL,1,WIJKL,1)
         ELSE
            CALL DAXPY(N2ASHX,D1,WIJKL,1,RIJKL,1)
         END IF
      GO TO 100
C
 4    CONTINUE
C        ... DISTYP = 4
         KLOFF = (L-1)*NASHT + K
         KKL   = N2ASHX*(KL-1) + 1
         IF (SLREOR) THEN
            CALL DCOPY(NNASHX,FIJKL(KLOFF),N2ASHX,RIJKL,1)
            CALL DSPTSI(NASHT,RIJKL,WIJKL)
            CALL DAXPY(N2ASHX,DM1,FIJKL(KKL),1,WIJKL,1)
         ELSE
            CALL DCOPY(NNASHX,FIJKL(KLOFF),N2ASHX,WIJKL,1)
            CALL DSPTSI(NASHT,WIJKL,RIJKL)
            CALL DAXPY(N2ASHX,DM1,FIJKL(KKL),1,RIJKL,1)
         END IF
      GO TO 100
C
 5    CONTINUE
C        ... DISTYP = 5
         KLOFF = (K-1)*NASHT + L
         IF (SLREOR) THEN
            CALL DCOPY(NNASHX,FIJKL(KLOFF),N2ASHX,RIJKL,1)
            CALL DSPTSI(NASHT,RIJKL,WIJKL)
            CALL MTRSP(NASHT,NASHT,
     &                 FIJKL(N2ASHX*(KL-1)+1),NASHT,RIJKL,NASHT)
            CALL DAXPY(N2ASHX,DM1,RIJKL,1,WIJKL,1)
         ELSE
            CALL DCOPY(NNASHX,FIJKL(KLOFF),N2ASHX,WIJKL,1)
            CALL DSPTSI(NASHT,WIJKL,RIJKL)
            CALL MTRSP(NASHT,NASHT,
     &                 FIJKL(N2ASHX*(KL-1)+1),NASHT,WIJKL,NASHT)
            CALL DAXPY(N2ASHX,DM1,WIJKL,1,RIJKL,1)
         END IF
      GO TO 100
C
 6    CONTINUE
C        ... DISTYP = 6
C        WIJKL(I,J) = 3(i^j|k l)+3(i j|k^l); FIJKL(I,J,K,L) = (i j|k^l)
         KL1   = 1 + (K-1 + NASHT*(L-1))*N2ASHX
         CALL DCOPY(N2ASHX,FIJKL(KL1),1,WIJKL,1)
C
         KL = K + (L-1)*NASHT
         DO 68 J=1,NASHT
            JKLOFF = ((J-1)*NASHT - 1)*N2ASHX + KL
         IF (SLREOR) THEN
            DO 66 I=1,NASHT
               WIJKL(I,J) = D3*(WIJKL(I,J) + FIJKL(JKLOFF + I*N2ASHX))
 66         CONTINUE
         ELSE
            DO 67 I=1,NASHT
               RIJKL(I,J) = D3*(WIJKL(I,J) + FIJKL(JKLOFF + I*N2ASHX))
 67         CONTINUE
         END IF
 68      CONTINUE
      GO TO 100
C
    8 CONTINUE
C        ... DISTYP = 8
C        ... DISTYP = 25
C        WIJKL(I,J) = (i^j|k l)-(i j|k^l); FIJKL(I,J,K,L) = (i j|k^l)
         KL1   = 1 + (L-1 + NASHT*(K-1))*N2ASHX
         CALL DCOPY(N2ASHX,FIJKL(KL1),1,WIJKL,1)
C        ... using -(i j|k^l) = (i j|l^k)
C
         KL = K + (L-1)*NASHT
         DO 88 J=1,NASHT
            JKLOFF = ((J-1)*NASHT - 1)*N2ASHX + KL
         IF (SLREOR) THEN
            DO 86 I=1,NASHT
               WIJKL(I,J) =  WIJKL(I,J) + FIJKL(JKLOFF + I*N2ASHX)
 86         CONTINUE
         ELSE
            DO 87 I=1,NASHT
               RIJKL(I,J) =  WIJKL(I,J) + FIJKL(JKLOFF + I*N2ASHX)
 87         CONTINUE
         END IF
 88      CONTINUE
      GO TO 100
C
C
C     ... DISTYP = 9
    9 CONTINUE
         KL = (L - 1) * NASHT + K
         KKL = N2ASHX*(KL-1) + 1
         IF (SLREOR) THEN
            CALL DCOPY(N2ASHX,FIJKL(KKL),1,WIJKL,1)
         ELSE
            CALL DCOPY(N2ASHX,FIJKL(KKL),1,RIJKL,1)
         END IF
      GO TO 100
C
C
   10 CONTINUE
C     ... DISTYP = 10
         KL = (K - 1) * NASHT + L
         KKL = N2ASHX*(KL-1) + 1
         IF (SLREOR) THEN
            CALL MTRSP(NASHT,NASHT,FIJKL(KKL),NASHT,WIJKL,NASHT)
         ELSE
            CALL MTRSP(NASHT,NASHT,FIJKL(KKL),NASHT,RIJKL,NASHT)
         END IF
      GO TO 100
C
C  W(I,J)=(I~~J|KL) + (K~~L|IJ) + (I~J|K~L)  + (K~L|I~J)
C
   11 CONTINUE
C     ... DISTYP = 11
      KLIJ = (L-1)*NASHT + K
      IJKL = (KLIJ-1)*N2ASHX + 1
      CALL DCOPY(N2ASHX,FIJKL(IJKL),1,WIJKL,1)
      CALL DAXPY(N2ASHX,D1,FIJKL(KLIJ),N2ASHX,WIJKL,1)
      CALL DAXPY(N2ASHX,D1,FIJKL(IJKL + N2ASHX*N2ASHX),1,WIJKL,1)
      CALL DAXPY(N2ASHX,D1,FIJKL(KLIJ + N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
      IF (.NOT.SLREOR) THEN
         CALL DCOPY(N2ASHX,WIJKL,1,RIJKL,1)
      END IF
      GO TO 100
C
C  W(I,J)=(J~~I|LK) + (L~~K|JI) + (J~I|L~K)  + (L~K|J~I)
C
   12 CONTINUE
C     ... DISTYP = 12
      LKJI = (K-1)*NASHT + L
      JILK = (LKJI-1)*N2ASHX + 1
      IF (SLREOR) THEN
         CALL DCOPY(N2ASHX,FIJKL(JILK),1,RIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(LKJI),N2ASHX,RIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(JILK+N2ASHX*N2ASHX),1,RIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(LKJI+N2ASHX*N2ASHX),N2ASHX,RIJKL,1)
         CALL MTRSP(NASHT,NASHT,RIJKL,NASHT,WIJKL,NASHT)
      ELSE
         CALL DCOPY(N2ASHX,FIJKL(JILK),1,WIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(LKJI),N2ASHX,WIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(JILK+N2ASHX*N2ASHX),1,WIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(LKJI+N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
         CALL MTRSP(NASHT,NASHT,WIJKL,NASHT,RIJKL,NASHT)
      END IF
      GO TO 100
C
C  W(I,J)=S*(I~~J|KL) + (K~~L|IJ) + S1*(I~J|K~L)  + S2*(K~L|I~J)
C
   13 CONTINUE
C     ... DISTYP = 13
C     ... DISTYP = 21
      KLIJ = (L-1)*NASHT + K
      IJKL = (KLIJ-1)*N2ASHX + 1
      CALL DCOPY(N2ASHX,FIJKL(KLIJ),N2ASHX,WIJKL,1)
      CALL DAXPY(N2ASHX,D1*ISGN1*ISGN2,FIJKL(IJKL),1,WIJKL,1)
      CALL DAXPY(N2ASHX,D1*ISGN1,
     *           FIJKL(IJKL + N2ASHX*N2ASHX),1,WIJKL,1)
      CALL DAXPY(N2ASHX,D1*ISGN2,
     *           FIJKL(KLIJ + N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
      IF (.NOT.SLREOR) THEN
         CALL DCOPY(N2ASHX,WIJKL,1,RIJKL,1)
      END IF
      GO TO 100
C
C  W(I,J)=S*(J~~I|LK) + (L~~K|JI) + S1*(J~I|L~K)  + S2*(L~K|J~I)
C
   14 CONTINUE
C     ... DISTYP = 14
C     ... DISTYP = 22
      LKJI = (K-1)*NASHT + L
      JILK = (LKJI-1)*N2ASHX + 1
      IF (SLREOR) THEN
         CALL DCOPY(N2ASHX,FIJKL(LKJI),N2ASHX,RIJKL,1)
         CALL DAXPY(N2ASHX,D1*ISGN1*ISGN2,FIJKL(JILK),1,RIJKL,1)
         CALL DAXPY(N2ASHX,D1*ISGN1,
     *              FIJKL(JILK + N2ASHX*N2ASHX),1,RIJKL,1)
         CALL DAXPY(N2ASHX,D1*ISGN2,
     *              FIJKL(LKJI + N2ASHX*N2ASHX),N2ASHX,RIJKL,1)
         CALL MTRSP(NASHT,NASHT,RIJKL,NASHT,WIJKL,NASHT)
      ELSE
         CALL DCOPY(N2ASHX,FIJKL(LKJI),N2ASHX,WIJKL,1)
         CALL DAXPY(N2ASHX,D1*ISGN1*ISGN2,FIJKL(JILK),1,WIJKL,1)
         CALL DAXPY(N2ASHX,D1*ISGN1,
     *              FIJKL(JILK + N2ASHX*N2ASHX),1,WIJKL,1)
         CALL DAXPY(N2ASHX,D1*ISGN2,
     *              FIJKL(LKJI + N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
         CALL MTRSP(NASHT,NASHT,WIJKL,NASHT,RIJKL,NASHT)
      END IF
      GO TO 100
C
C                  W(I,J)=(I~J|KL) + (K~L|IJ)
C
 15   CONTINUE
C     ... DISTYP = 15
         KL = (L-1)*NASHT + K
         KKL = (KL-1)*N2ASHX + 1
         IF (SLREOR) THEN
            CALL DCOPY(N2ASHX,FIJKL(KKL),1,WIJKL,1)
            CALL DAXPY(N2ASHX,D1,FIJKL(KL),N2ASHX,WIJKL,1)
         ELSE
            CALL DCOPY(N2ASHX,FIJKL(KKL),1,RIJKL,1)
            CALL DAXPY(N2ASHX,D1,FIJKL(KL),N2ASHX,RIJKL,1)
         END IF
      GO TO 100
C
C                  W(I,J)=(J~I|LK) + (L~K|JI)
C
   16 CONTINUE
C     ... DISTYP = 16
         KL = (K-1)*NASHT + L
         KKL = (KL-1)*N2ASHX + 1
         IF (SLREOR) THEN
            CALL DCOPY(N2ASHX,FIJKL(KKL),1,RIJKL,1)
            CALL DAXPY(N2ASHX,D1,FIJKL(KL),N2ASHX,RIJKL,1)
            CALL MTRSP(NASHT,NASHT,RIJKL,NASHT,WIJKL,NASHT)
         ELSE
            CALL DCOPY(N2ASHX,FIJKL(KKL),1,WIJKL,1)
            CALL DAXPY(N2ASHX,D1,FIJKL(KL),N2ASHX,WIJKL,1)
            CALL MTRSP(NASHT,NASHT,WIJKL,NASHT,RIJKL,NASHT)
         END IF
      GO TO 100
C
C                  W(J,I)=S1*(I~J|KL) + S2*(K~L|IJ)
C
   17 CONTINUE
C     ... DISTYP = 17
         KL = (L-1)*NASHT + K
         KKL = (KL-1)*N2ASHX + 1
         IF (SLREOR) THEN
            CALL DCOPY(N2ASHX,FIJKL(KKL),1,WIJKL,1)
            IF (ISGN1.EQ.-1) CALL DSCAL(N2ASHX,DM1,WIJKL,1)
            CALL DAXPY(N2ASHX,D1*ISGN2,FIJKL(KL),N2ASHX,WIJKL,1)
         ELSE
            CALL DCOPY(N2ASHX,FIJKL(KKL),1,RIJKL,1)
            IF (ISGN1.EQ.-1) CALL DSCAL(N2ASHX,DM1,RIJKL,1)
            CALL DAXPY(N2ASHX,D1*ISGN2,FIJKL(KL),N2ASHX,RIJKL,1)
         END IF
      GO TO 100
C
C                  W(I,J)=S1*(J~I|LK) + S2*(L~K|JI)
C
   18 CONTINUE
C     ... DISTYP = 18
         KL = (K-1)*NASHT + L
         KKL = (KL-1)*N2ASHX + 1
         IF (SLREOR) THEN
            CALL DCOPY(N2ASHX,FIJKL(KKL),1,RIJKL,1)
            IF (ISGN1.EQ.-1) CALL DSCAL(N2ASHX,DM1,RIJKL,1)
            CALL DAXPY(N2ASHX,D1*ISGN2,FIJKL(KL),N2ASHX,RIJKL,1)
            CALL MTRSP(NASHT,NASHT,RIJKL,NASHT,WIJKL,NASHT)
         ELSE
            CALL DCOPY(N2ASHX,FIJKL(KKL),1,WIJKL,1)
            IF (ISGN1.EQ.-1) CALL DSCAL(N2ASHX,DM1,WIJKL,1)
            CALL DAXPY(N2ASHX,D1*ISGN2,FIJKL(KL),N2ASHX,WIJKL,1)
            CALL MTRSP(NASHT,NASHT,WIJKL,NASHT,RIJKL,NASHT)
         END IF
      GO TO 100
C
C  W(I,J) = 3[ (I~J|K^L) + (K~L|I^J) + (IJ|K~^L) + (KL|I~^J)]
C
   19 CONTINUE
C     ... DISTYP = 19
      KLIJ = (L-1)*NASHT + K
      IJKL = (KLIJ-1)*N2ASHX + 1
      CALL DCOPY(N2ASHX,FIJKL(IJKL),1,WIJKL,1)
      CALL DSCAL(N2ASHX,D3,WIJKL,1)
      CALL DAXPY(N2ASHX,D3,FIJKL(KLIJ),N2ASHX,WIJKL,1)
      CALL DAXPY(N2ASHX,D3,FIJKL(IJKL + N2ASHX*N2ASHX),1,WIJKL,1)
      CALL DAXPY(N2ASHX,D3,FIJKL(KLIJ + N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
      IF (.NOT.SLREOR) THEN
         CALL DCOPY(N2ASHX,WIJKL,1,RIJKL,1)
      END IF
      GO TO 100
C
C  W(I,J)=3[ (J~I|L^K) + (L~K|J^I) + (JI|L~^K) + (LK|J~^I)]
C
   20 CONTINUE
C     ... DISTYP = 20
      LKJI = (K-1)*NASHT + L
      JILK = (LKJI-1)*N2ASHX + 1
      IF (SLREOR) THEN
         CALL DCOPY(N2ASHX,FIJKL(JILK),1,RIJKL,1)
         CALL DSCAL(N2ASHX,D3,RIJKL,1)
         CALL DAXPY(N2ASHX,D3,FIJKL(LKJI),N2ASHX,RIJKL,1)
         CALL DAXPY(N2ASHX,D3,FIJKL(JILK+N2ASHX*N2ASHX),1,RIJKL,1)
         CALL DAXPY(N2ASHX,D3,FIJKL(LKJI+N2ASHX*N2ASHX),N2ASHX,RIJKL,1)
         CALL MTRSP(NASHT,NASHT,RIJKL,NASHT,WIJKL,NASHT)
      ELSE
         CALL DCOPY(N2ASHX,FIJKL(JILK),1,WIJKL,1)
         CALL DSCAL(N2ASHX,D3,WIJKL,1)
         CALL DAXPY(N2ASHX,D3,FIJKL(LKJI),N2ASHX,WIJKL,1)
         CALL DAXPY(N2ASHX,D3,FIJKL(JILK+N2ASHX*N2ASHX),1,WIJKL,1)
         CALL DAXPY(N2ASHX,D3,FIJKL(LKJI+N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
         CALL MTRSP(NASHT,NASHT,WIJKL,NASHT,RIJKL,NASHT)
      END IF
      GO TO 100
C
C W(I,J) = 3(IJ|K^L) + 3(KL|I^J)
C
   23 CONTINUE
C     ... DISTYP = 23
         KLIJ = (L-1)*NASHT + K
         IJKL = (KLIJ-1)*N2ASHX + 1
         IF (SLREOR) THEN
            CALL DCOPY(N2ASHX,FIJKL(IJKL),1,WIJKL,1)
            CALL DSCAL(N2ASHX,D3,WIJKL,1)
            CALL DAXPY(N2ASHX,D3,FIJKL(KLIJ),N2ASHX,WIJKL,1)
         ELSE
            CALL DCOPY(N2ASHX,FIJKL(IJKL),1,RIJKL,1)
            CALL DSCAL(N2ASHX,D3,RIJKL,1)
            CALL DAXPY(N2ASHX,D3,FIJKL(KLIJ),N2ASHX,RIJKL,1)
         END IF
      GO TO 100
C
C W(I,J) = 3(JI|L^K) + 3(LK|J^I)
C
   24 CONTINUE
C     ... DISTYP = 24
         LKJI = (K-1)*NASHT + L
         JILK = (LKJI-1)*N2ASHX + 1
         IF (SLREOR) THEN
            CALL DCOPY(N2ASHX,FIJKL(JILK),1,RIJKL,1)
            CALL DSCAL(N2ASHX,D3,RIJKL,1)
            CALL DAXPY(N2ASHX,D3,FIJKL(LKJI),N2ASHX,RIJKL,1)
            CALL MTRSP(NASHT,NASHT,RIJKL,NASHT,WIJKL,NASHT)
         ELSE
            CALL DCOPY(N2ASHX,FIJKL(JILK),1,WIJKL,1)
            CALL DSCAL(N2ASHX,D3,WIJKL,1)
            CALL DAXPY(N2ASHX,D3,FIJKL(LKJI),N2ASHX,WIJKL,1)
            CALL MTRSP(NASHT,NASHT,WIJKL,NASHT,RIJKL,NASHT)
         END IF
      GO TO 100
C
C  W(I,J)= (KL|I^J) - (IJ|K^L)
C
   25 CONTINUE
C     ... DISTYP = 25
         KLIJ = (L-1)*NASHT + K
         IJKL = (KLIJ-1)*N2ASHX + 1
         IF (SLREOR) THEN
            CALL DCOPY(N2ASHX,FIJKL(IJKL),1,WIJKL,1)
            IF (ISGN1.EQ.-1) CALL DSCAL(N2ASHX,DM1,WIJKL,1)
            CALL DAXPY(N2ASHX,D1*ISGN2,FIJKL(KLIJ),N2ASHX,WIJKL,1)
         ELSE
            CALL DCOPY(N2ASHX,FIJKL(IJKL),1,RIJKL,1)
            IF (ISGN1.EQ.-1) CALL DSCAL(N2ASHX,DM1,RIJKL,1)
            CALL DAXPY(N2ASHX,D1*ISGN2,FIJKL(KLIJ),N2ASHX,RIJKL,1)
         END IF
      GO TO 100
C
C  W(I,J)= (LK|J^I) - (JI|L^K)
C
   26 CONTINUE
C     ... DISTYP = 26
         LKJI = (K-1)*NASHT + L
         JILK = (LKJI-1)*N2ASHX + 1
         IF (SLREOR) THEN
            CALL DCOPY(N2ASHX,FIJKL(JILK),1,RIJKL,1)
            IF (ISGN1.EQ.-1) CALL DSCAL(N2ASHX,DM1,RIJKL,1)
            CALL DAXPY(N2ASHX,D1*ISGN2,FIJKL(LKJI),N2ASHX,RIJKL,1)
            CALL MTRSP(NASHT,NASHT,RIJKL,NASHT,WIJKL,NASHT)
         ELSE
            CALL DCOPY(N2ASHX,FIJKL(JILK),1,WIJKL,1)
            IF (ISGN1.EQ.-1) CALL DSCAL(N2ASHX,DM1,WIJKL,1)
            CALL DAXPY(N2ASHX,D1*ISGN2,FIJKL(LKJI),N2ASHX,WIJKL,1)
            CALL MTRSP(NASHT,NASHT,WIJKL,NASHT,RIJKL,NASHT)
         END IF
      GO TO 100
C
C     W(I,J) = (I~~~J|KL) + (I~~J|K~L) + (I~~J|K~L) + (I~~J|K~L) 
C            + (K~~~L|IJ) + (K~~L|I~J)  +(K~~L|I~J) + (K~~L|I~J)
C
   27 CONTINUE
C     ... DISTYP = 27
      KLIJ = (L-1)*NASHT + K
      IJKL = (KLIJ-1)*N2ASHX + 1
      CALL DCOPY(N2ASHX,FIJKL(IJKL),1,WIJKL,1)
      CALL DAXPY(N2ASHX,D1,FIJKL(KLIJ),N2ASHX,WIJKL,1)
      CALL DAXPY(N2ASHX,D1,FIJKL(IJKL + N2ASHX*N2ASHX),1,WIJKL,1)
      CALL DAXPY(N2ASHX,D1,FIJKL(KLIJ + N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
      CALL DAXPY(N2ASHX,D1,FIJKL(IJKL + 2*N2ASHX*N2ASHX),1,WIJKL,1)
      CALL DAXPY(N2ASHX,D1,FIJKL(KLIJ + 2*N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
      CALL DAXPY(N2ASHX,D1,FIJKL(IJKL + 3*N2ASHX*N2ASHX),1,WIJKL,1)
      CALL DAXPY(N2ASHX,D1,FIJKL(KLIJ + 3*N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
      IF (.NOT.SLREOR) THEN
         CALL DCOPY(N2ASHX,WIJKL,1,RIJKL,1)
      END IF
      GO TO 100
C
C     W(I,J) = (J~~~I|LK) + (J~~I|L~K) + (J~~I|L~K) + (J~~I|L~K) 
C            + (L~~~K|JI) + (L~~K|J~I)  +(L~~K|J~I) + (L~~K|J~I)
C
   28 CONTINUE
C     ... DISTYP = 28
      LKJI = (K-1)*NASHT + L
      JILK = (LKJI-1)*N2ASHX + 1
      IF (SLREOR) THEN
         CALL DCOPY(N2ASHX,FIJKL(JILK),1,RIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(LKJI),N2ASHX,RIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(JILK+N2ASHX*N2ASHX),1,RIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(LKJI+N2ASHX*N2ASHX),N2ASHX,RIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(JILK+2*N2ASHX*N2ASHX),1,RIJKL,1)
        CALL DAXPY(N2ASHX,D1,FIJKL(LKJI+2*N2ASHX*N2ASHX),N2ASHX,RIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(JILK+3*N2ASHX*N2ASHX),1,RIJKL,1)
        CALL DAXPY(N2ASHX,D1,FIJKL(LKJI+3*N2ASHX*N2ASHX),N2ASHX,RIJKL,1)
         CALL MTRSP(NASHT,NASHT,RIJKL,NASHT,WIJKL,NASHT)
      ELSE
         CALL DCOPY(N2ASHX,FIJKL(JILK),1,WIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(LKJI),N2ASHX,WIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(JILK+N2ASHX*N2ASHX),1,WIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(LKJI+N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(JILK+2*N2ASHX*N2ASHX),1,WIJKL,1)
        CALL DAXPY(N2ASHX,D1,FIJKL(LKJI+2*N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
         CALL DAXPY(N2ASHX,D1,FIJKL(JILK+3*N2ASHX*N2ASHX),1,WIJKL,1)
        CALL DAXPY(N2ASHX,D1,FIJKL(LKJI+3*N2ASHX*N2ASHX),N2ASHX,WIJKL,1)
         CALL MTRSP(NASHT,NASHT,WIJKL,NASHT,RIJKL,NASHT)
      END IF
      GO TO 100
C
  100 CONTINUE
C
C     Reorder from Sirius order to Lunar order
C
      IF (SLREOR) THEN
         CALL REORMT(WIJKL,RIJKL,NASHT,NASHT,ILTSOB,ILTSOB)
C        CALL REORMT(AIN  ,AOUT ,NROW ,NCOL ,IROW  ,ICOL  )
      END IF
C
      RETURN
C     ... end of getint.
      END
