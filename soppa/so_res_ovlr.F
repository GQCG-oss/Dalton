C
C  /* Deck so_res_ovlr */
      SUBROUTINE SO_RES_OVLR(RESO1E,LRESO1E,RESO1D,LRESO1D,TR1E,LTR1E,
     &                      TR1D,LTR1D)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, December 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Add RPA contributions to the overlap result vectors.
C              The RPA overlap matrix is equal to the identity matrix
C              (-identity for deexcitation part) and the RPA 
C              contributions are therefore given by the
C              trial vectors.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ONE = 1.0D0)
C
      DIMENSION RESO1E(LRESO1E), RESO1D(LRESO1D)
      DIMENSION TR1E(LTR1E),     TR1D(LTR1D)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_OVLR')
C
      CALL DAXPY(LTR1E,ONE,TR1E,1,RESO1E,1)
      CALL DAXPY(LTR1D,-ONE,TR1D,1,RESO1D,1)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_OVLR')
C
      RETURN
      END
