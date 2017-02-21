C
C  /* Deck dc_backup1 */
      SUBROUTINE DC_BACKUP1(ISYMTR, NEXCI, WORK, LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Pi A. B. Haase, January 2017
C
C     PURPOSE: Make backup of HRPA eigenvectors in case the HRPA(D)
C     method is used to ensure that HRPA eigenvectors are in place if
C     also the s-HRPA(D) method is requested.
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
      CALL QENTER('DC_BACKUP1')
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LTR1E = NT1AM(ISYMTR)
      LTR1D = NT1AM(ISYMTR)
C
      KRES1 = 1
      KEND1  = KRES1  + LTR1E
C
      LWORK1 = LWORK  - KEND1
C
      CALL SO_MEMMAX ('DC_BACKUP1.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('DC_BACKUP1.1',' ',KEND1,LWORK1)
C
C----------------
C     Open files.
C----------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
C
C---------------------------
C     Loop over excitations.
C---------------------------
C
      DO 100 IEXCI = 1,NEXCI
C
C-----------------------------------------------------------------------
C        Make a copy of HRPA vectors in the end of the file
C-----------------------------------------------------------------------
C
         CALL SO_READ(WORK(KRES1),  LTR1E, LUTR1E, FNTR1E, IEXCI)
C
         CALL SO_WRITE(WORK(KRES1), LTR1E, LUTR1E, FNTR1E, NEXCI+IEXCI)
C
         CALL SO_READ(WORK(KRES1),  LTR1D,LUTR1D,FNTR1D,IEXCI)
C
         CALL SO_WRITE(WORK(KRES1),  LTR1D,LUTR1D,FNTR1D,NEXCI+IEXCI)
C
 100  CONTINUE
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
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
      CALL QEXIT('DC_BACKUP1')
C
      RETURN
C
      END
