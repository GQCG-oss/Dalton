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
      PROGRAM NODPRG_DALTON
C
C     Node program (slave program) for PVM
C
#include "implicit.h"
#include "priunit.h"
#include "mtags.h"
!     wow thats dangerous... who has used the PVM code the last time??? sknecht - feb 2013
      PARAMETER (MWORK = 10 000 000)
      real(8), allocatable :: work(:)
!     memory-id tags for field 1 and 2+LMWORK in work (allocated in main driver routines)
      DATA WRKDLM/12345678D0/
C
C     Join PVM and return my name to master
C
      CALL PVMXMYTID(MYTID)
      CALL PVMXPARENT(MASTER)
      CALL SENDNAME(MASTER,MYTID)
C
C     Get calculation type (iprtyp)
C
      CALL PVMXRECV(MASTER,MTAG1)
      CALL PVMXUNPACK(IPRTYP,1,'INTEGER')
C
      DO WHILE(IPRTYP.GT.0)         
         CALL PVMXUNPACK(IPRINT,1,'INTEGER')
C
         IF (IPRTYP.EQ.1) THEN
            CALL HER_NODSTR(MWORK,wrkldm,IPRINT)
         ELSE IF (IPRTYP.EQ.2) THEN
           ! allocate memory using f90 utilities
           allocate(work(mwork+2),stat=i)
!          Set memory traps
           work(0)       = wrkdlm
           work(1+mwork) = wrkdlm
           CALL ERI_NODSTR(WORK(1),MWORK,IPRINT)
           deallocate(work)
C         ELSE IF (IPRTYP.EQ.3) THEN
C            CALL GEMNOD(WORK,MWORK,MASTER,MYTID,MYNUM,IPRINT)
         END IF
C
         CALL PVMXRECV(MASTER,MTAG1)
         CALL PVMXUNPACK(IPRTYP,1,'INTEGER')
      END DO
      CALL PVMXEXIT
      STOP
      END
