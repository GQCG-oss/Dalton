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
C$Id: hergam.F,v 1.1.1.1 2001-02-08 13:33:20 hjj Exp $
C
#include <single.h>
C  /* Deck gamfun */
      SUBROUTINE GAMFUN
C
C     Trygve Ulf Helgaker fall 1984
C
C     This subroutine calculates the incomplete gamma function as
C     described by McMurchie & Davidson, J. Comp. Phys. 26 (1978) 218.
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <pi.h>
      PARAMETER (D1 = 1.D0,  D2 = 2.D0, D3 = 3.D0, D6 = 6.D0,
     &           D10 = 10.D0, D24 = 24.D0, D120 = 120.D0,
     &           D720 = 720.D0,
     &           HALF = D1/D2, TENTH = .1D0, TEN6 = 1.0D6)
      PARAMETER (COEF2 = D1/D2,  COEF3 = - D1/D6,
     &           COEF4 = D1/D24, COEF5 = - D1/D120,
     &           COEF6 = D1/D720)
      PARAMETER (GFAC30 =  .4999489092 D0,
     *           GFAC31 = -.2473631686 D0,
     *           GFAC32 =  .321180909  D0,
     *           GFAC33 = -.3811559346 D0,
     *           GFAC20 =  .4998436875 D0,
     *           GFAC21 = -.24249438   D0,
     *           GFAC22 =  .24642845   D0,
     *           GFAC10 =  .499093162  D0,
     *           GFAC11 = -.2152832    D0,
     *           GFAC00 =  .490        D0)
C
#include <gamcom.h>
C
      SAVE MAXJ0, SQRPIH
      DATA MAXJ0 /-1/
C
#include <sdpre.h>
C
      IPOINT = NINT(D10*MIN(WVAL,TEN6))
      IF (IPOINT .LT. 120) THEN
         ISTART = 1 + 121*JMAX0 + IPOINT
         WDIF = WVAL - TENTH*SDPRE(IPOINT)
         FJW(JMAX0) = (((((COEF6*TABFJW(ISTART + 726)*WDIF
     *                   + COEF5*TABFJW(ISTART + 605))*WDIF
     *                    + COEF4*TABFJW(ISTART + 484))*WDIF
     *                     + COEF3*TABFJW(ISTART + 363))*WDIF
     *                      + COEF2*TABFJW(ISTART + 242))*WDIF
     *                       - TABFJW(ISTART + 121))*WDIF
     *                        + TABFJW(ISTART)
         D2WAL = D2*WVAL
         REXPW = EXP(-WVAL)
         DENOM = SDPRE(2*JMAX0 + 1)
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
      SQRPIH = SQRT(PI)/D2
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
     *      ' GAMTAB ERROR: JMX =',JMX,', which is greater than',MAXJ
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
      D2MAX1 = SDPRE(JMAX + JMAX + 1)
      R2MAX1 = D1/D2MAX1
      DO 800 IPOINT = 1,120
         WVAL = TENTH*SDPRE(IPOINT)
         D2WAL = WVAL + WVAL
         IADR = IADR + 1
         TERM = R2MAX1
         SUM = TERM
         DENOM = D2MAX1
         DO 810 IORDER = 2,200
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
#include <implicit.h>
C
C     Calculates Cartesian powers
C     tuh oct 27 90
C
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <xyzpow.h>
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
C  /* Deck lmnval */
      SUBROUTINE LMNVAL(NHKTA,KHKTA,LVALUE,MVALUE,NVALUE)
C
C     tuh March 87
C
#include <implicit.h>
#include <maxmom.h>
#include <xyzpow.h>
      DIMENSION LVALUE(KHKTA), MVALUE(KHKTA), NVALUE(KHKTA)
#include <ibtfun.h>
C
      ICOMP = 0
      DO 100 I = 1, KHKTA
         ICOMP = ICOMP + 1
         LVALUE(ICOMP) = NHKTA - ISTEP(I)
         MVALUE(ICOMP) = MVAL(I)
         NVALUE(ICOMP) = NVAL(I)
  100 CONTINUE
      RETURN
      END
