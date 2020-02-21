C
C  /* Deck rp_init */
      SUBROUTINE RP_INIT(WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, Januar 1998
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Stephan P. A. Sauer: 08.11.2013: Initialization of the CC common blocks
C                                      is moved to new routine CCSD_INIT0
C
C     PURPOSE: Initialize some pointers and sort AO-integrals
C     into distributions.
C
#include "implicit.h"
#include "priunit.h"
C
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "r12int.h"
C
      DIMENSION WORK(LWORK)
      CHARACTER WORD*7
      LOGICAL OPENED
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_INIT')
C
C----------------------------------------
C     Check whether SIRIFC is still open.
C----------------------------------------
C
      CALL GPINQ('SIRIFC','OPENE',OPENED)
      IF (OPENED) THEN
         INQUIRE (FILE='SIRIFC',NUMBER=LUSIFC)
         WRITE(LUPRI,'(2A,I3)') ' RP_INIT: file SIRIFC is already ',
     &                          'opened with unit number ',LUSIFC
         CALL GPCLOSE (LUSIFC,'KEEP')
         WRITE(LUPRI,'(2A,I3)') ' RP_INIT: file SIRIFC was closed again'
      ENDIF
C
CSPAS:08/11-2013: Initialization of the CC common blocks
C                 is moved to new routine CCSD_INIT0
C
C---------------------------------------------------------
C     Initialize common blocks CCSDINP, CCLR, and CCSDSYM.
C---------------------------------------------------------
C
      CALL CCSD_INIT0(WORD)
C
C-------------------------------------------------------------------------
C     Set addition variables in common blocks CCSDINP, CCSDSYM and CMMOLL.
C-------------------------------------------------------------------------
C
      ISYMOP = 1
C
      KEEPAOIN  = .TRUE.
      KEEPAOTWO = 2
CKeinSPASmehr
C
C     Use LABEL (WK/UniKA/04-11-2002).
      LABEL = 'TRCCINT '
C
      CALL CCSD_INIT1(WORK,LWORK)
C
      IF (.NOT. DIRECT) THEN
         CALL CCSD_SORTAO(WORK,LWORK)
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('RP_INIT')
C
      RETURN
      END
