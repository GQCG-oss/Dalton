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
C FILE : erisph.F
C
#include <single.h>
C
C  /* Deck sphmat */
      SUBROUTINE SPHMAT(CSQ,NROWS,NCOLS,NDER,GDER,WORK,LWORK,IPRINT)
C
C     T. Helgaker
C 
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
      INTEGER C, D, E, F, G, X, Y, Z
      LOGICAL GDER
#include <nuclei.h>
#include <ccom.h>
      DIMENSION CSQ(NROWS,NCOLS), WORK(LWORK)
C
      IF (DOCART) THEN
         WRITE (LUPRI,'(1X,2A)') 
     &     ' Only spherical-harmonic basis is implemented in ERI.',
     &     ' Run TWOINT instead.' 
         CALL QUIT('ERI does not work with Cartesian basis.')
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL TITLER('Output from SPHMAT','*',103)
         WRITE (LUPRI,'(1X,A,2I5)')    ' NROWS,NCOLS: ',NROWS,NCOLS
         WRITE (LUPRI,'(1X,A,I5,L5 )') ' NDER, GDER:  ',NDER, GDER
      END IF
C
      NDIM = NROWS*NCOLS
      CALL DZERO(CSQ,NDIM)
C
C     IF ((     GDER.AND.NCOLS.NE.1) .OR.
C    &    (.NOT.GDER.AND.NCOLS.NE.NATOMS)) THEN
C        WRITE (LUPRI,'(1X,A,I5)') ' NCOLS:', NCOLS
C        CALL QUIT('Error in SPHMAT: incorrect NCOLS')
C     END IF 
C
      DO D = 0, NDER
      DO X = D, 0, -1
      DO Y = D - X, 0, -1
         Z = D - X - Y 
         DO L = 0, NHTYP - 1
         DO M = 0, 2*L 
            IF (GDER) THEN 
               CALL SPHDER(CSQ,NROWS,NCOLS,L,M,X,Y,Z,0,0,0,IPRINT)
            ELSE
               CALL SPHDER(CSQ,NROWS,NCOLS,L,M,0,0,0,X,Y,Z,IPRINT)
            END IF
         END DO
         END DO
      END DO
      END DO
      END DO
C
      IF (IPRINT.GT.10) THEN
         IF (NDIM.GT.LWORK) CALL STOPIT('SPHMAT',' ',LWORK,NDIM)
         CALL MTRSP(NROWS,NCOLS,CSQ,NROWS,WORK,NCOLS)
         CALL HEADER('CSQ in SPHMAT',-1) 
         ISTART = 0
         DO IATOM = 1, NCOLS
            IF (.NOT.GDER) WRITE (LUPRI,'(//,2X,A,I5)')
     &         ' B-field CSQ for atom ', IATOM 
            DO D = 0, NDER
            DO X = D, 0, -1
            DO Y = D - X, 0, -1
               Z = D - X - Y 
               DO L = 0, NHTYP - 1
               DO C = -D, D 
               IF (L + C .GE. 0) THEN
                  NS = 2*L+1
                  NC = (L+C+1)*(L+C+2)/2
                  IADR = ISTART + KSQADR(L,X,Y,Z,C)
                  WRITE(LUPRI,'(//,1X,A,6X,3I2,/27X,A,I2,/27X,A,I2)')
     &               ' Spherical tra. matrix for derivative:',X,Y,Z,
     &               ' angular momentum:',L, ' Cartesian level: ',C
                  CALL OUTPUT(WORK(IADR),1,NS,1,NC,NS,NC,1,LUPRI) 
               END IF
               END DO
               END DO
            END DO
            END DO
            END DO
            ISTART = ISTART + NROWS 
         END DO
      END IF
C
      RETURN
      END
C
C  /* Deck sphder */
      SUBROUTINE SPHDER(CSQ,NROWS,NCOLS,L,M1,I,J,K,E,F,G,IPRINT)
C
C     T. Helgaker 10.03.00
C
C     This routine calculates coefficient matrices for transformation
C     from an undifferentiated Cartesian basis to the differentatied
C     spherical-harmonic basis. 
C
C     Order of spherical harmonics:         L, M
C     Order of geometry derivatives:        I, J, K
C     Order of magnetic field derivatives:  E, F, G 
C       
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
      PARAMETER (DM1 = -1.0D0, D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
#include <ccom.h>
#include <maxorb.h>
#include <symmet.h>
#include <nuclei.h>
      INTEGER A, B, C, E, F, G, EE, FF, GG, P, Q, R,
     &        T, U, V, X, Y, Z, V0, ALPHA, TOT, POS, PX, PY, PZ,
     &        BX, BY, BZ
      DIMENSION CSQ(NROWS,NCOLS)
      PWD(I,J)   = FACULT(I)/FACULT(I-J)
      HCF(I,J,K) = FACULT(I)/(FACULT(J)*FACULT(K/2)*FACULT(I-J-K))
#include <ibtfun.h>
C
      M = M1 - L
      IF (L.EQ.1) THEN
        IF (M .EQ. -1) MADR =  0  
        IF (M .EQ.  0) MADR =  1 
        IF (M .EQ.  1) MADR = -1 
      ELSE
        MADR = M
      END IF
C
      NIJK = I + J + K
      NEFG = E + F + G
      MABS = ABS(M)
      V0 = 0
      IF (M .LT. 0) V0 = 1 
      FACNRM = D1
      IF (M .NE. 0) FACNRM = SQRT(D2*FACULT(L+MABS)*FACULT(L-MABS))/
     &                       (FACULT(L)*(D2**MABS))
      FACNRM = FACNRM*DM1**((NEFG-MOD(NEFG,2))/2)*D2**(-NEFG)
      FACNRM = FACNRM/SQRT(FACUL2(2*L-1))
      NDER = NIJK + NEFG 
      IOFF = NSPHER(L,I+E,J+F,K+G)
      DO 100 EE = 0, E
      DO 100 FF = 0, F
      DO 100 GG = 0, G
         BX = FF-GG+G
         BY = GG-EE+E
         BZ = EE-FF+F
         FAC1 = FACNRM*BINOM(E,EE)*BINOM(F,FF)*BINOM(G,GG)
         DO 200 II = MAX(0,I-F+FF-GG), I
         DO 200 JJ = MAX(0,J-G+GG-EE), J
         DO 200 KK = MAX(0,K-E+EE-FF), K
            FAC2 = FAC1*BINOM(I,II)*PWD(F-FF+GG,I-II)*
     &                  BINOM(J,JJ)*PWD(G-GG+EE,J-JJ)*
     &                  BINOM(K,KK)*PWD(E-EE+FF,K-KK)
            PX = F-FF+GG-I+II
            PY = G-GG+EE-J+JJ
            PZ = E-EE+FF-K+KK
            DO 300 T = 0, L - MABS, 2
            DO 300 U = 0, T, 2
            DO 300 V = V0, MABS, 2
               FAC3 = FAC2*BINOM(L,T/2)*BINOM(L-T/2,MABS+T/2)
     &                    *BINOM(T/2,U/2)*BINOM(MABS,V)
               DO 400 A = 0, MIN(II,T+MABS-U-V+BX)
               DO 400 B = 0, MIN(JJ,U+V+BY)
               DO 400 C = 0, MIN(KK,L-T-MABS+BZ)
                  FAC4 = FAC3*PWD(T+MABS-U-V+BX,A)
     &                       *PWD(U+V+BY,B)
     &                       *PWD(L-T-MABS+BZ,C)
                  DO 500 P = 0, II - A, 2
                  DO 500 Q = 0, JJ - B, 2
                  DO 500 R = 0, KK - C, 2
                     FACTOR = DM1**(A+B+C+EE+FF+GG+(T+V-V0-P-Q-R)/2)*
     &                        D2**(II+JJ+KK-A-B-C-P-Q-R-T)*FAC4*
     &                        HCF(II,A,P)*HCF(JJ,B,Q)*HCF(KK,C,R)
                     X = T+MABS-U-V+II-2*A-P+BX
                     Y = U+V+JJ-2*B-Q+BY
                     Z = L-T-MABS+KK-2*C-R+BZ
                     TOT = X + Y + Z
                     IADR = IOFF + NSPHAD(L,TOT-L,NDER) +
     &                      (2*L+1)*(NCRT(X,Y,Z)-1) + L + MADR
                     IF (NCOLS.EQ.1) THEN
                        CSQ(IADR,1) = CSQ(IADR,1) + FACTOR 
                        IF (IPRINT.GT.20) THEN
                           IATOM = 0
                           BFAC  = 0.0D0
                           ALPHA = II+JJ+KK-A-B-C-(P+Q+R)/2
                           WRITE (LUPRI,'(1X,2I3,2(2X,3I2),I3,
     &                               F10.5,I4,2X,2(3I3,2X),4I3)') 
     &                               L,M,I,J,K,E,F,G,
     &                               IATOM,BFAC,ALPHA,
     &                               PX,PY,PZ,X,Y,Z,
     &                               TOT,NCRT(X,Y,Z),IADR
                        END IF
                     ELSE
                        IA = 0
                        DO 700 IATOM = 1, NUCIND
                           AX0 = CORD(1,IATOM)
                           AY0 = CORD(2,IATOM)
                           AZ0 = CORD(3,IATOM)
                           DO 800 LA = 0, MAXOPR 
                           IF (IBTAND(ISTBNU(IATOM),LA).EQ.0) THEN
                              IA = IA + 1
                              AX = PT(IBTAND(ISYMAX(1,1),LA))*AX0
                              AY = PT(IBTAND(ISYMAX(2,1),LA))*AY0
                              AZ = PT(IBTAND(ISYMAX(3,1),LA))*AZ0
                              BFAC = FACTOR*(AX**PX)*(AY**PY)*(AZ**PZ)
                              CSQ(IADR,IA) = CSQ(IADR,IA) + BFAC
                              IF (IPRINT.GT.20) THEN
                                 ALPHA = II+JJ+KK-A-B-C-(P+Q+R)/2
                                 WRITE (LUPRI,'(1X,2I3,2(2X,3I2),I3,
     &                             F10.5,I4,2X,2(3I3,2X),3I3)') 
     &                             L,M,I,J,K,E,F,G,
     &                             IA,BFAC,ALPHA,
     &                             PX,PY,PZ,X,Y,Z,
     &                             TOT,NCRT(X,Y,Z),IADR
                              END IF
                           END IF
  800                      CONTINUE 
  700                   CONTINUE 
                     END IF
  500             CONTINUE 
  400          CONTINUE
  300       CONTINUE 
  200    CONTINUE
  100 CONTINUE 
      RETURN
      END
C  /* Deck nspher */
      INTEGER FUNCTION NSPHER(L,X,Y,Z)
#include <implicit.h>
#include <maxaqn.h>
#include <ccom.h>
      INTEGER D, L, X, Y, Z 
C
      D = X + Y + Z 
      LMAX = NHTYP - 1
      NSPHER = (24 + 2*(1 + 2*D)*L*
     &          (-2 + (3 + 2*D + 2*D**2)*L + 8*L**2 + 3*L**3)
     &       + (D*(2 + 3*D + D**2)*(1 + LMAX)*
     &         (-54 - 134*LMAX - 85*LMAX**2 - 15*LMAX**3 - 
     &          12*D**2*(1 + LMAX) + 20*D**3*(1 + LMAX) + 
     &          D*(166 + 406*LMAX + 255*LMAX**2 + 45*LMAX**3)))/30
     &       + (1 + 2*D)*(1 + LMAX)*
     &        (12 + 28*LMAX + 17*LMAX**2 + 3*LMAX**3 + 2*D*(1 + LMAX) +
     &          2*D**2*(1 + LMAX))*(y + y**2 + 2*y*z + z*(3 + z)))/24
      RETURN
      END 
C  /* Deck lspher */
      INTEGER FUNCTION LSPHER(NDER)
      INTEGER NDER
      LSPHER = NSPHER(0,NDER+1,0,0) - 1
      RETURN
      END 
C  /* Deck sphdim */
      SUBROUTINE SPHDIM
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <ericom.h>
      NCSQ1 = NSPHER(0,MAXDER+1,0,0) - 1
      IF (BDER) THEN
         NCSQ2 = NATOMS
      ELSE
         NCSQ2 = 1
      END IF 
      RETURN
      END
C  /* Deck nsphad */
      INTEGER FUNCTION NSPHAD(L,LEVEL,NDER)
#include <implicit.h>
      NSPHAD = ((LEVEL + NDER)*(1 + 2*L)*(2 + LEVEL**2 
     &        + NDER**2 + 6*L + 3*L**2 - 3*NDER*(1 + L) 
     &        + LEVEL*(3 - NDER + 3*L)))/6
      RETURN
      END
C  /* Deck ncrt */
      INTEGER FUNCTION NCRT(I,J,K)
#include <implicit.h>
      NCRT = 1 + J + 2*K + (J + K)*(J + K - 1)/2
      RETURN
      END 
C  /* Deck lmcrt */
      INTEGER FUNCTION LMCRT(L,M,I,J,K)
#include <implicit.h>
      LMCRT = (2*L+1)*(NCRT(I,J,K) - 1) + L + M + 1
      RETURN
      END 
C  /* Deck sphadr */
      INTEGER FUNCTION KSQADR(L,X,Y,Z,C)
#include <implicit.h>
#include <maxaqn.h>
#include <mxdr.h>
      LOGICAL FIRST
      INTEGER A, B, C, D, X, Y, Z
      INTEGER ACOEF(0:MXQN-1,0:MXDR), BCOEF(0:MXQN-1,0:MXDR)
#include <ccom.h>
C
      SAVE FIRST, ACOEF, BCOEF, NHOLD
      DATA FIRST,NHOLD /.TRUE.,0/
C
      IF (FIRST .OR. NHTYP.NE.NHOLD) THEN
         DO J = 0, MXQN - 1
         DO D = 0, MXDR
            CALL SPHAB(J,D,A,B)
            ACOEF(J,D) = A
            BCOEF(J,D) = B
         END DO
         END DO
         FIRST = .FALSE.
         NHOLD = NHTYP
      END IF 
      D = X + Y + Z
      IXYZ = Y + Y**2 + 2*Y*Z + Z*(3 + Z)
      KSQADR = (ACOEF(L,D) + BCOEF(L,D)*IXYZ)/24 + NSPHAD(L,C,D)
      RETURN
      END
C  /* Deck sphab */
      SUBROUTINE SPHAB(L,D,A,B)
#include <implicit.h>
#include <maxaqn.h>
      INTEGER D, A, B
#include <ccom.h>
      LX = NHTYP - 1
      A = 24 +  2*(1+2*D)*L*(-2+(3+2*D+2*D**2)*L+8*L**2+3*L**3)
     &       +  (D*(2+3*D + D**2)*(1+LX)*
     &          (-54-134*LX-85*LX**2-15*LX**3 - 
     &          12*D**2*(1+LX) + 20*D**3*(1+LX) + 
     &          D*(166+406*LX+255*LX**2+45*LX**3)))/30
      B = (1+2*D)*(1+LX)*(12+28*LX+17*LX**2+3*LX**3+2*D*(1+LX) + 
     &     2*D**2*(1+LX))
      RETURN
      END
C --- end of erisph.F ---
