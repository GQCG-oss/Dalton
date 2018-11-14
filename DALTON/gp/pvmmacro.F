!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
C-------------------------------------------------------------
C    Paal Dahle Jan-2001
C
C    Macros for PVM calls that make calls to PVM routines
C    easier and more transparent because:
C
C    a) Less arguments to PVMX... calls than to PVMF... calls
C    b) No need to specify the fpvm3.h include file in the
C       calling subroutine.
C    c) Error handling is automatically taken care of.
C-------------------------------------------------------------
C
C  /* Deck pvmxspawn */
      SUBROUTINE PVMXSPAWN(TASK,NTASK,TIDS)
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      CHARACTER TASK*6
      INTEGER TIDS
      DIMENSION TIDS(NTASK)
C
#if defined (VAR_PVM)
      CALL PVMFSPAWN(TASK,PVMDEFAULT,'*',NTASK,TIDS,NUMT)
C
      IF (NUMT.NE.NTASK) THEN
         WRITE(LUPRI,*) (TIDS(I), I=1,NTASK)
         CALL PVMXFAIL(-40)
      END IF
#endif
C
      RETURN
      END
C  /* Deck pvmxmytid */
      SUBROUTINE PVMXMYTID(TID)
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      INTEGER TID
C
#if defined (VAR_PVM)
      CALL PVMFMYTID(TID)
      IF (TID.LT.0) CALL PVMXFAIL(TID)
#endif
C
      RETURN
      END
C  /* Deck pvmxparent */
      SUBROUTINE PVMXPARENT(TID)
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      INTEGER TID
C
#if defined (VAR_PVM)
      CALL PVMFPARENT(TID)
      IF (TID.LT.0) CALL PVMXFAIL(TID)
#endif
C
      RETURN
      END
C  /* Deck pvmxexit */
      SUBROUTINE PVMXEXIT
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
C
#if defined (VAR_PVM)
      CALL PVMFEXIT(INFO)
      IF (INFO.LT.0) CALL PVMXFAIL(INFO)
#endif
C
      RETURN
      END
C  /* Deck pvmxhalt */
      SUBROUTINE PVMXHALT
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
C
#if defined (VAR_PVM)
      CALL PVMFHALT(INFO)
      IF (INFO.LT.0) CALL PVMXFAIL(INFO)
#endif
C
      RETURN
      END
C  /* Deck pvmxinitsend */
      SUBROUTINE PVMXINITSEND(ENCODE)
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      INTEGER ENCODE, BUFID
C
#if defined (VAR_PVM)
      CALL PVMFINITSEND(ENCODE,BUFID)
      IF (BUFID.LT.0) CALL PVMXFAIL(BUFID)
#endif
C
      RETURN
      END
C  /* Deck pvmxpack */
      SUBROUTINE PVMXPACK(BUFFER,COUNT,TYPE)
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      CHARACTER*6 TYPE
      INTEGER   DATATYPE, COUNT, STRIDE
      DIMENSION BUFFER(*)
C
#if defined (VAR_PVM)
      IF      (TYPE .EQ. 'INTEGE' .OR. TYPE .EQ. 'LOGICA') THEN
         DATATYPE = INTEGER4
#ifdef VAR_INT64
  ERROR: INTEGER*8 not implemented in pvmmacro.F routines
#endif
      ELSE IF (TYPE .EQ. 'DOUBLE') THEN
         DATATYPE = REAL8
      ELSE IF (TYPE .EQ. 'STRING') THEN
         DATATYPE = STRING
      ELSE
         WRITE(LUPRI,*) ' TYPE :',TYPE,' does not exist!'
         CALL QUIT('ERROR in PVMXPACK: Nonexisting DATATYPE')
      END IF
      STRIDE = 1
      CALL PVMFPACK(DATATYPE,BUFFER,COUNT,STRIDE,INFO)
      IF (INFO .LT. 0) CALL PVMXFAIL(INFO)
#endif
C
      RETURN
      END
C  /* Deck pvmxunpack */
      SUBROUTINE PVMXUNPACK(BUFFER,COUNT,TYPE)
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      CHARACTER*6 TYPE
      INTEGER DATATYPE, COUNT, STRIDE
      DIMENSION BUFFER(*)
C
#if defined (VAR_PVM)
      IF      (TYPE .EQ. 'INTEGE' .OR. TYPE .EQ. 'LOGICA') THEN
         DATATYPE = INTEGER4
#ifdef VAR_INT64
  ERROR: INTEGER*8 not implemented in pvmmacro.F routines
#endif
      ELSE IF (TYPE .EQ. 'DOUBLE') THEN
         DATATYPE = REAL8
      ELSE IF (TYPE .EQ. 'STRING') THEN
         DATATYPE = STRING
      ELSE
         WRITE(LUPRI,*) ' TYPE :',TYPE,' does not exist!'
         CALL QUIT('ERROR in PVMXUNPACK: Nonexisting DATATYPE')
      END IF
      STRIDE = 1
      CALL PVMFUNPACK(DATATYPE,BUFFER,COUNT,STRIDE,INFO)
      IF (INFO .LT. 0) CALL PVMXFAIL(INFO)
#endif
C
      RETURN
      END
C  /* Deck pvmxsend */
      SUBROUTINE PVMXSEND(DEST,TAG)
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      INTEGER DEST, TAG
C
#if defined (VAR_PVM)
      CALL PVMFSEND(DEST,TAG,INFO)
      IF (INFO .LT. 0) CALL PVMXFAIL(INFO)
#endif
C
      RETURN
      END
C  /* Deck pvmxrecv */
      SUBROUTINE PVMXRECV(SOURCE,TAG)
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      INTEGER SOURCE, TAG
C
#if defined (VAR_PVM)
      CALL PVMFRECV(SOURCE,TAG,INFO)
      IF (INFO .LT. 0) CALL PVMXFAIL(INFO)
#endif
C
      RETURN
      END
C  /* Deck pvmxmcast */
      SUBROUTINE PVMXMCAST(SLAVES,TIDS,TAG)
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      INTEGER TIDS, SLAVES, TAG
      DIMENSION TIDS(SLAVES)
C
#if defined (VAR_PVM)
      CALL PVMFMCAST(SLAVES,TIDS,TAG,INFO)
      IF (INFO .LT. 0) CALL PVMXFAIL(INFO)
#endif
C
      RETURN
      END
C  /* Deck pvmxfail */
      SUBROUTINE PVMXFAIL(NTAG)
C
#include "implicit.h"
#include "priunit.h"
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      CHARACTER FAILS(0:40)*70
C
      FAILS( 0) = 'PvmOk         (okay)'
      FAILS( 2) = 'PvmBadParam   (bad parameter)'
      FAILS( 3) = 'PvmMismatch   (barrier count mismatch)'
      FAILS( 5) = 'PvmNoData     (read past end of buffer)'
      FAILS( 6) = 'PvmNoHost     (no such host)'
      FAILS( 7) = 'PvmNoFile     (no such executable)'
      FAILS(10) = 'PvmNoMem      (cannot get memory)'
      FAILS(12) = 'PvmBadMsg     (cannot decode received message)'
      FAILS(14) = 'PvmSysErr     (pvmd not responding)'
      FAILS(15) = 'PvmNoBuf      (no current buffer)'
      FAILS(16) = 'PvmNoSuchBuf  (bad message identifier)'
      FAILS(17) = 'PvmNullGroup  (null group name is illegal)'
      FAILS(18) = 'PvmDupGroup   (already in group)'
      FAILS(19) = 'PvmNoGroup    (no group with that name)'
      FAILS(20) = 'PvmNotInGroup (not in group)'
      FAILS(21) = 'PvmNoInst     (no such instance in group)'
      FAILS(22) = 'PvmHostFail   (host failed)'
      FAILS(23) = 'PvmNoParent   (no parent task)'
      FAILS(24) = 'PvmNotImpl    (function not implemented)'
      FAILS(25) = 'PvmDSysErr    (pvmd system error)'
      FAILS(26) = 'PvmBadVersion (pvmd-pvmd protocol mismatch)'
      FAILS(27) = 'PvmOutOfRes   (out of resources)'
      FAILS(28) = 'PvmDupHost    (host already configured)'
      FAILS(29) = 'PvmCantStart  (failed to execute new slave pvmd)'
      FAILS(30) = 'PvmAlready    (already doing operation)'
      FAILS(31) = 'PvmNoTask     (no such task)'
      FAILS(32) = 'PvmNoEntry    (no such (group.instance))'
      FAILS(33) = 'PvmDupEntry   ((group.instance) already exists)'
      FAILS(40) = 'Failed to spawn process(es)'
C
      MTAG = ABS(NTAG)
      WRITE(LUPRI,'(//1X,A)') FAILS(MTAG)
C
      CALL QUIT('Error detected in PVM. Please consult dalton output!')
C
      RETURN
      END
