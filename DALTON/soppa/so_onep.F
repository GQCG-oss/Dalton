C
C  /* Deck so_onep */
      SUBROUTINE SO_ONEP (PRP1,LPRP1,LABEL,ISYMTR,RTNLBL,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997.
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Read one-electron AO property integrals from file.
C
#include "implicit.h"
#include "priunit.h"
C
#include "ccorb.h"
#include "ccsdsym.h"
#include "inftap.h"
#include "soppinf.h"
CClark:7/1/2016
#include "cbiexc.h"
CClark:end
C
      LOGICAL     FNDLB2
      DIMENSION   PRP1(LPRP1), WORK(LWORK)
      PARAMETER   (DM1 = -1.0D0, D0 = 0.0D0, D1 =1.0D0 )
      CHARACTER*8 LABEL, RTNLBL(2),
CClark:7/1/2016
     &            FNAME
CClark:end
      LOGICAL     EXTST, OPENED
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_ONEP')
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LPAOI  = N2BASX
C
      KPAOI   = 1
      KEND1   = KPAOI  + LPAOI
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_ONEP',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_ONEP',' ',KEND1,LWORK)
C
C------------------------------------------------------------
C     Open file containing one-electron AO property integrals
C     and rewind.
C------------------------------------------------------------
C
CClark:20/11/2015
      IF (LABEL(1:3) .EQ. 'COS' .OR. LABEL(1:3) .EQ. 'SIN') THEN
C
         FNAME = 'AOGOS   '
         LUFILE = LUAOGOS
C
      ELSE
         FNAME = 'AOPROPER'
         LUFILE = LUPROP
C
      ENDIF
CClark:end
C
      CALL GPINQ(FNAME,'OPENE',OPENED)
      IF (OPENED) THEN
         INQUIRE (FILE=FNAME,NUMBER=LUFILE)
C         WRITE(LUPRI,'(4A,I3)') ' SO_ONEP: file',FNAME,' is already',
C     &                          ' opened with unit number ',LUFILE
      ELSE
         IF (LUFILE .LE. 0) CALL GPOPEN(LUFILE,FNAME,'OLD',' ',
     &                                  'UNFORMATTED',IDUMMY,.FALSE.)
      ENDIF
      REWIND LUFILE
C
C-----------------------
C     Read AO integrals.
C-----------------------
C
      IF ( FNDLB2(LABEL,RTNLBL,LUFILE,LUERR)) THEN
C
         IF ( RTNLBL(2) .EQ. 'SQUARE  ' ) THEN
            WRITE(LUPRI,*) 'WARNING: The code for square property ',
     &                 'integrals is not tested'
            CALL READT(LUFILE,N2BASX,WORK(KPAOI))
         ELSE IF ( RTNLBL(2) .EQ. 'SYMMETRI' ) THEN
            CALL READT(LUFILE,NNBASX,WORK(KPAOI))
         ELSE IF ( RTNLBL(2) .EQ. 'ANTISYMM' ) THEN
            CALL READT(LUFILE,NNBASX,WORK(KPAOI))
         ELSE
            CALL QUIT('Error: No antisymmetry label on '//FNAME)
         END IF
C
      ELSE
C
         WRITE(LUPRI,'(/4A)') 'ERROR: ',LABEL,
     &                        ' integrals not on ',FNAME
         CALL QUIT('Property integrals not found on '//TRIM(FNAME)//'.')
C
      END IF
C
C----------------
C     Close file.
C----------------
C
CClark:8/1/2016
      IF (LABEL(1:3) .NE. 'COS' .AND. LABEL(1:3) .NE. 'SIN') THEN
C
         CALL GPCLOSE (LUFILE,'KEEP')
C
      ENDIF
CClark:end
C
C---------------------------------------------------------------------
C     Repack the integrals to the order used in the CC and SOPPA code.
C---------------------------------------------------------------------
C
      CALL SO_INTRP(PRP1,LRPR1,WORK(KPAOI),LPAOI,RTNLBL,ISYMTR)
C
C-------------------------------------------
C     Write repacked AO integrals to output.
C-------------------------------------------
C
      IF (IPRSOP .GE. 8) THEN
         CALL AROUND(LABEL//'repacked one electron AO integrals')
         CALL CC_PRFCKAO(PRP1,ISYMTR)
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_ONEP')
C
      RETURN
      END
