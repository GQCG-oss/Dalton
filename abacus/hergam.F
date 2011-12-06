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
C FILE    : hergam.F
C
C  /* Deck gamfun */
      SUBROUTINE GAMFUN
C
C     Trygve Ulf Helgaker fall 1984
C
C     This subroutine calculates the incomplete gamma function as
C     described by McMurchie & Davidson, J. Comp. Phys. 26 (1978) 218.
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "pi.h"
      PARAMETER (D1 = 1.D0,  D2 = 2.D0, D10 = 10.D0,
     &           HALF = 0.5D0, TENTH = 0.1D0, TEN6 = 1.0D6)
      PARAMETER (SQRPIH = SQRTPI/D2)
      PARAMETER (COEF2 = HALF,  COEF3 = - D1/6.0D0,
     &           COEF4 = D1/24.D0, COEF5 = - D1/120.0D0,
     &           COEF6 = D1/720.D0)
      PARAMETER (GFAC30 =  .4999489092 D0,
     &           GFAC31 = -.2473631686 D0,
     &           GFAC32 =  .321180909  D0,
     &           GFAC33 = -.3811559346 D0,
     &           GFAC20 =  .4998436875 D0,
     &           GFAC21 = -.24249438   D0,
     &           GFAC22 =  .24642845   D0,
     &           GFAC10 =  .499093162  D0,
     &           GFAC11 = -.2152832    D0,
     &           GFAC00 =  .490        D0)
C
#include "gamcom.h"
C
      SAVE MAXJ0
      DATA MAXJ0 /-1/
C
      IPOINT = D10*MIN(WVAL,TEN6) + HALF
!     have seen problems with NINT intrinsic function here (rarely)
!     therefore the "+ HALF" before integer truncation
      IF (IPOINT .LT. 0) THEN
         write (luerr,*) 'GAMFUN warning, ipoint=',ipoint
         write (luerr,*) 'd10,wval, ten6',d10,wval,ten6
         write (luerr,*) d10*min(wval,ten6)
         call quit('Fatal error in gamfun')
      ELSE IF (IPOINT .LT. 120) THEN
         ISTART = 1 + 121*JMAX0 + IPOINT
         WDIF = WVAL - TENTH*IPOINT
         FJW(JMAX0) = (((((COEF6*TABFJW(ISTART + 726)*WDIF     ! 726 = 6*121
     &                   + COEF5*TABFJW(ISTART + 605))*WDIF
     &                    + COEF4*TABFJW(ISTART + 484))*WDIF
     &                     + COEF3*TABFJW(ISTART + 363))*WDIF
     &                      + COEF2*TABFJW(ISTART + 242))*WDIF
     &                       - TABFJW(ISTART + 121))*WDIF
     &                        + TABFJW(ISTART)
         D2WAL = D2*WVAL
         REXPW = EXP(-WVAL)
         DENOM = 2.0D0*JMAX0 + 1.0D0
         DO 100 J = JMAX0, 1, -1
            DENOM = DENOM - D2
            FJW(J - 1) = (D2WAL*FJW(J) + REXPW)/DENOM
  100    CONTINUE
      ELSE IF (IPOINT .LE. 20*JMAX0 + 360) THEN
         RWVAL = D1/WVAL
         REXPW = EXP(-WVAL)
         GVAL = GFAC30 + RWVAL*(GFAC31 + RWVAL*(GFAC32 + RWVAL*GFAC33))
         FJW(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
         FACTOR = HALF*RWVAL
         TERM = FACTOR*REXPW
         DO 200 J = 1, JMAX0
            FJW(J) = FACTOR*FJW(J - 1) - TERM
            FACTOR = RWVAL + FACTOR
  200    CONTINUE
      ELSE
         RWVAL  = D1/WVAL
         FJW(0) = SQRPIH*SQRT(RWVAL)
         FACTOR = HALF*RWVAL
         DO 300 J = 1, JMAX0
            FJW(J) = FACTOR*FJW(J-1)
            FACTOR = RWVAL + FACTOR
  300    CONTINUE
      END IF
      RETURN
C
C     ***** Tabulation of incomplete gamma function *****
C
      ENTRY GAMTAB(JMX)
C
C     For J = JMX a power series expansion is used, see for
C     example Eq.(39) given by V. Saunders in "Computational
C     Techniques in Quantum Chemistry and Molecular Physics",
C     Reidel 1975.  For J < JMX the values are calculated
C     using downward recursion in J.
C
C
      IF (JMX .GT. MAXJ) THEN
         WRITE (LUPRI,'(//A,I5,A,I3)')
     &      ' GAMTAB ERROR: JMX =',JMX,', which is greater than',MAXJ
         WRITE (LUERR,'(//A,I5,A,I3)')
     &      ' GAMTAB ERROR: JMX =',JMX,', which is greater than',MAXJ
         CALL QUIT('GAMTAB ERROR: JMX greater than limit.')
      END IF
      JMAX = JMX + 6
      MAXJ0 = JMAX
C
C     WVAL = 0.0
C
      IADR = 1
      DENOM = D1
      DO 700 J = 0,JMAX
         TABFJW(IADR) = D1/DENOM
         IADR = IADR + 121
         DENOM = DENOM + D2
  700 CONTINUE
C
C     WVAL = 0.1, 0.2, 0.3,... 12.0
C
      IADR = IADR - 121
      D2MAX1 = 2.0D0*JMAX + 1.0D0
      R2MAX1 = D1/D2MAX1
      DO 800 IPOINT = 1,120
         WVAL = TENTH*IPOINT
         D2WAL = WVAL + WVAL
         IADR = IADR + 1
         TERM = R2MAX1
         SUM = TERM
         DENOM = D2MAX1
         DO 810 IORDER = 2, 200
            DENOM = DENOM + D2
            TERM = TERM*D2WAL/DENOM
            SUM = SUM + TERM
            IF (TERM .LE. 1.0D-15) GO TO 820
  810    CONTINUE
  820    CONTINUE
         REXPW = EXP(-WVAL)
         TABFJW(IADR) = REXPW*SUM
         DENOM = D2MAX1
         JADR = IADR
         DO 830 J = 1,JMAX
            DENOM = DENOM - D2
            TABFJW(JADR - 121) = (TABFJW(JADR)*D2WAL + REXPW)/DENOM
            JADR = JADR - 121
  830    CONTINUE
  800 CONTINUE
      RETURN
      END
C  /* Deck carpow */
      SUBROUTINE CARPOW
#include "implicit.h"
C***********************************************************************
C     Calculates Cartesian powers i+j+k=L
C     tuh oct 27 90
C
C     For given L there are (L+1)(L+2)/2 distinct triplets (i,j,k).
C     They can therefore be generated by looping over the lower
C     triangle of a square matrix of dimension L+1.
C
C     From an element (a,b) of the matrix the values of the triplet
C     are given by
C
C        i = L + 1 - a; j = a - b; k = b - 1
C
C     The inverse relations are:
C
C        b = k + 1; a = j + k + 1
C
C     Commentary added by T. Saue Jun 23 2007

!     radovan: in case you are interested
!              here is a python routine that does that
!              which i use in one of my scripts
!     def lmn(k):
!         for i in range(1, k + 2):
!             for j in range(1, i + 1):
!                 l = k + 1 - i
!                 m = i - j
!                 n = j - 1
!                 print l, m, n
C
C***********************************************************************
#include "priunit.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "xyzpow.h"
      IF (MXQN .GT. MXQNM) THEN
         WRITE (LUPRI,'(//A,I5,A,I3)')
     &      ' CARPOW ERROR: MXQN =',MXQN,
     &      ', which is greater than MXQNM =',MXQNM
         WRITE (LUPRI,'(/A,I5)') ' Increase MXQNM to ',MXQN,' or more.'
         CALL QUIT('CARPOW ERROR: MXQN greater than MXQNM.')
      END IF
      IJ = 0
      DO 100 I = 1, MXQNM
         DO 200 J = 1, I
            IJ = IJ + 1
            ISTEP(IJ) = I
            MVAL(IJ) = I - J
            NVAL(IJ) = J - 1
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C -- end of hergam.F --
