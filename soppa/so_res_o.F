C  /* Deck so_res_o */
      SUBROUTINE SO_RES_O(RESO1E,LRESO1E,RESO1D,LRESO1D,TR1E,LTR1E,
     &                    TR1D,LTR1D,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                    ISYRES,ISYMTR)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate the 1p1h and 1h1p S[2] transformed trial
C              vectors.
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION RESO1E(LRESO1E), RESO1D(LRESO1D)
      DIMENSION TR1E(LTR1E),     TR1D(LTR1D)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB)
C
#include <soppinf.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_O')
C
C----------------------------------------------
C     Calculate the overlap matrix in eq. (67).
C----------------------------------------------
C
      CALL DZERO(RESO1E,LRESO1E)
      CALL DZERO(RESO1D,LRESO1D)
C
      DTIME      = SECOND()
      CALL SO_RES_OVL(RESO1E,LRESO1E,RESO1D,LRESO1D,TR1E,LTR1E,
     &                TR1D,LTR1D,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                ISYRES,ISYMTR)
      DTIME      = SECOND()   - DTIME
      SOTIME(22) = SOTIME(22) + DTIME
C
C--------------------------------------
C     Calculate the RPA overlap matrix.
C--------------------------------------
C
      DTIME      = SECOND()
      CALL SO_RES_OVLR(RESO1E,LRESO1E,RESO1D,LRESO1D,
     &                 TR1E,LTR1E,TR1D,LTR1D)
      DTIME      = SECOND()   - DTIME
      SOTIME(23) = SOTIME(23) + DTIME
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_O')
C
      RETURN
      END
