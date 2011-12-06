C  /* Deck so_tfset */
      SUBROUTINE SO_TFSET(ARRAY,NTRIAL,LTR,ISYMTR)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, September 1996
C
C     PURPOSE: Transform a set of righthand side trial vectors to
C              a set of left hand side trial vectors which represents
C              the same vectors.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D+00, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
C---------------------------------
C     Dimensions of the arguments.
C---------------------------------
C
      DIMENSION ARRAY(1)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TFSET')
C
C---------------------------
C     Read a set of vectors.
C---------------------------
C
      DO I = 1,NTRIAL
         KOFF = 1 + LTR * (I - 1)
         CALL CCSD_TCMEPKX(ARRAY(KOFF),TWO,ISYMTR)
         CALL DSCAL(LTR,HALF,ARRAY(KOFF),1)
      END DO
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TFSET')
C
      RETURN
      END
