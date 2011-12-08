C
C  /* Deck so_getmo */
      SUBROUTINE SO_GETMO(CMO,LCMO,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Get MO coefficients and reorder.
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <inftap.h>
C
      DIMENSION CMO(LCMO), WORK(LWORK)
C     LOGICAL   OPENY
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_GETMO')
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LCMOT = NLAMDS
C
      KCMOT   = 1
      KEND1   = KCMOT + LCMOT
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX('SO_GETMO        ',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_GETMO',' ',KEND1,LWORK)
C
C----------------------------------------------
C     Read MO-coefficients from interface file.
C----------------------------------------------
C
      IF (LUSIFC .LE. O) CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUERR)
      READ (LUSIFC)
 
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMOT+I-1), I=1,LCMOT)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C---------------------------------------
C     Reorder the MO-coefficient matrix.
C---------------------------------------
C
      CALL CMO_REORDER(WORK(KCMOT),WORK(KEND1),LWORK1)
C
C--------------------------------------------------
C     Copy the "not frozen" MO-coefficients to CMO.
C--------------------------------------------------
C
      CALL DCOPY(LCMO,WORK(KCMOT),1,CMO,1)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_GETMO')
C
      RETURN
      END
