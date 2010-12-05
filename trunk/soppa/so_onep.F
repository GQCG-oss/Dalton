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
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <soppinf.h>
C
      LOGICAL     FNDLB2
      DIMENSION   PRP1(LPRP1), WORK(LWORK)
      PARAMETER   (DM1 = -1.0D0, D0 = 0.0D0, D1 =1.0D0 )
      CHARACTER*8 LABEL, RTNLBL(2)
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
      CALL GPINQ('AOPROPER','OPENE',OPENED)
      IF (OPENED) THEN
         INQUIRE (FILE='AOPROPER',NUMBER=LUPROP)
         WRITE(LUPRI,'(2A,I3)') ' SO_ONEP: file AOPROPER is already',
     &                          ' opened with unit number ',LUPROP
      ELSE
         IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER','OLD',' ',
     &                                  'UNFORMATTED',IDUMMY,.FALSE.)
      ENDIF
      REWIND LUPROP
C
C-----------------------
C     Read AO integrals.
C-----------------------
C
      IF ( FNDLB2(LABEL,RTNLBL,LUPROP,LUERR)) THEN
C
         IF ( RTNLBL(2) .EQ. 'SQUARE  ' ) THEN
            WRITE(LUPRI,*) 'WARNING: The code for square property ',
     &                 'integrals is not tested'
            CALL READT(LUPROP,N2BASX,WORK(KPAOI))
         ELSE IF ( RTNLBL(2) .EQ. 'SYMMETRI' ) THEN
            CALL READT(LUPROP,NNBASX,WORK(KPAOI))
         ELSE IF ( RTNLBL(2) .EQ. 'ANTISYMM' ) THEN
            CALL READT(LUPROP,NNBASX,WORK(KPAOI))
         ELSE
            CALL QUIT('Error: No antisymmetry label on LUPROP')
         END IF
C
      ELSE
C
         WRITE(LUPRI,'(/3A)') 'ERROR: ',LABEL,' integrals not on LUPROP'
         CALL QUIT('Property integrals not found on LUPROP.')
C
      END IF
C
C----------------
C     Close file.
C----------------
C
      CALL GPCLOSE (LUPROP,'KEEP')
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
