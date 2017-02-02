C
C  /* Deck dc_backup2 */
      SUBROUTINE DC_BACKUP2(ISYMTR, NEXCI, WORK, LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Pi A. B. Haase, January 2017
C
C     PURPOSE: Replace HRPA(D) eigenvectors with backed up HRPA 
C     eigenvectors which are needed for s-HRPA(D) calculation
C
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION WORK(LWORK)
      INTEGER   ISYMTR, NEXCI
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('DC_BACKUP2')
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LTR1E = NT1AM(ISYMTR)
      LTR1D = NT1AM(ISYMTR)
      LTR2E = N2P2HOP(ISYMTR)
      LTR2D = N2P2HOP(ISYMTR)
C
      KRES1 = 1
      KEND1  = KRES1  + LTR1E + LTR2E
C
      LWORK1 = LWORK  - KEND1
C
      CALL SO_MEMMAX ('DC_BACKUP2.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('DC_BACKUP2.1',' ',KEND1,LWORK1)
C
C----------------
C     Open files.
C----------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
C
C---------------------------
C     Loop over excitations.
C---------------------------
C
      DO 100 IEXCI = 1,NEXCI
C
C-----------------------------------------------------------------------
C        Move HRPA vectors back to the beginning of file (i.e. overwrite
C        HRPA(D) vectors)
C-----------------------------------------------------------------------
C
         CALL SO_READ(WORK(KRES1),  LTR1E,LUTR1E,FNTR1E,NEXCI+IEXCI)
C
         CALL SO_WRITE(WORK(KRES1),  LTR1E,LUTR1E,FNTR1E,IEXCI)
C
         CALL SO_READ(WORK(KRES1),  LTR1D,LUTR1D,FNTR1D,NEXCI+IEXCI)
C
         CALL SO_WRITE(WORK(KRES1),  LTR1D,LUTR1D,FNTR1D,IEXCI)
C
C-----------------------------------------------------------------------
C        Replace doubles part of HRPA vectors with zeros.
C-----------------------------------------------------------------------
C
         CALL DZERO(WORK(KRES1), LTR2E)
C
         CALL SO_WRITE(WORK(KRES1), LTR2E, LUTR2E, FNTR2E, IEXCI)
C
         CALL SO_WRITE(WORK(KRES1), LTR2D, LUTR2D, FNTR2D, IEXCI)
C
 100  CONTINUE
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
C
C------------------------------------
C     Flush the standard output file.
C------------------------------------
C
      CALL FLSHFO(LUPRI)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('DC_BACKUP2')
C
      RETURN
C
      END
