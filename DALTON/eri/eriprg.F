!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
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
