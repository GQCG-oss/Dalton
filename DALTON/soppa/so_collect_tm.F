C
C  /* Deck so_collect_tm */
      SUBROUTINE SO_COLLECT_TM(ISYMTR,LABEL,IEXCI,TRMOM,TRLEN,TRVEL,
CSPAS:23/5-11: second and third moment sum rules
     &                         TQLEN,TQVEL,TTLEN,
CKeinSPASmehr
     &                         TRLON,TRMAG,TRGOS,BSRLON,EXENG)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, August 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Collect transition moments in the right arrays.
C
#include "implicit.h"
#include "priunit.h"
C
CSPAS:15/03-2006: merge with Dalton-2.0
C#include "cbilrs.h"
#include "cbiexc.h"
CKeinSPASmehr
#include "ccorb.h"
#include "mxcent.h"
#include "nuclei.h"
C
      CHARACTER*8 LABEL
C
      DIMENSION TRLEN(3,NSYM,MXNEXI),        TRVEL(3,NSYM,MXNEXI)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION TQLEN(3,3,NSYM,MXNEXI), TQVEL(3,3,NSYM,MXNEXI)
      DIMENSION TTLEN(10,NSYM,MXNEXI)
CKeinSPASmehr
      DIMENSION TRLON(3,NSYM,MXNEXI),        TRMAG(3,NSYM,MXNEXI)
      DIMENSION BSRLON(3,NSYM,MXNEXI),       EXENG(NSYM,MXNEXI)
CClark:7/1/2016
      DIMENSION TRGOS(3,MXNEXI)
CClark:end
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_COLLECT_TM')
C
C-----------------------------------------
C     Check that IEXCI is less that MXNEXI.
C-----------------------------------------
C
      IF ( IEXCI .GT. MXNEXI ) THEN
         WRITE(LUPRI,*) 'WARNING: IEXCI greater than MXNEXI'
         WRITE(LUPRI,*) 'Reset MXNEXI in the file cbiexc.h and resubmit'
         RETURN
      END IF
C
C------------------------------------------------------------
C     The calculated transition moments are <f|Operator|0>.
C     We want to keep them as <O|Operator|f>. Therefore, for
C     anti-hermitian operators we multiply with a minus-sign.
C------------------------------------------------------------
C
CRF :Handle first labels starting with a coordinate axis
      IDX1 = COORD_IDX(LABEL(1:1))
      IF (IDX1 .ge. 1 .and. IDX1 .LE. 3 ) THEN
C    Dipole - length gauge
         IF (LABEL(2:7).EQ.'DIPLEN') THEN
            TRLEN(IDX1,ISYMTR,IEXCI) = TRMOM
C    Dipole - velocity gauge
         ELSE IF (LABEL(2:7).EQ.'DIPVEL') THEN
            TRVEL(IDX1,ISYMTR,IEXCI) = - TRMOM
CSPAS:23/5-11: second and third moment sum rules
         ELSE IF (LABEL(3:8).EQ.'SECMOM') THEN
            IDX2 = COORD_IDX(LABEL(2:2))
            TQLEN(IDX1,IDX2,ISYMTR,IEXCI) = TRMOM
            IF (IDX1 .NE. IDX2 ) THEN
               TQLEN(IDX2,IDX1,ISYMTR,IEXCI) = TRMOM
            ENDIF
         ELSE IF (LABEL(5:8).EQ.'3MOM') THEN
            IF (LABEL(1:3) .EQ. 'XXX') THEN
                                TTLEN(1,ISYMTR,IEXCI) = TRMOM
            ELSE IF (LABEL(1:3) .EQ. 'YYY') THEN
                                TTLEN(2,ISYMTR,IEXCI) = TRMOM
            ELSE IF (LABEL(1:3) .EQ. 'ZZZ') THEN
                                TTLEN(3,ISYMTR,IEXCI) = TRMOM
            ELSE IF (LABEL(1:3) .EQ. 'XXY') THEN
                                TTLEN(4,ISYMTR,IEXCI) = TRMOM
            ELSE IF (LABEL(1:3) .EQ. 'XXZ') THEN
                                TTLEN(5,ISYMTR,IEXCI) = TRMOM
            ELSE IF (LABEL(1:3) .EQ. 'XYY') THEN
                                TTLEN(6,ISYMTR,IEXCI) = TRMOM
            ELSE IF (LABEL(1:3) .EQ. 'XYZ') THEN
                                TTLEN(7,ISYMTR,IEXCI) = TRMOM
            ELSE IF (LABEL(1:3) .EQ. 'XZZ') THEN
                                TTLEN(8,ISYMTR,IEXCI) = TRMOM
            ELSE IF (LABEL(1:3) .EQ. 'YYZ') THEN
                                TTLEN(9,ISYMTR,IEXCI) = TRMOM
            ELSE IF (LABEL(1:3) .EQ. 'YZZ') THEN
                                TTLEN(10,ISYMTR,IEXCI) = TRMOM
            END IF
CKeinSPASmehr
C
         ELSE IF (LABEL(2:7).EQ.'ANGMOM') THEN
C
            TRMAG(IDX1,ISYMTR,IEXCI) = - TRMOM
C
         ELSE IF (LABEL(2:7).EQ.'LONMAG') THEN
C
            TRLON(IDX1,ISYMTR,IEXCI) = - TRMOM
C
         ELSE
CRF         Something is wrong
            WRITE(LUPRI,'(A,/,2A,I3)') 'ERROR in SO_COLLECT_TM:',
     &      'Label not recognized: ', LABEL, IDX1
         END IF
CRF Handle labels that does not start with X,Y OR Z
      ELSE IF (LABEL(2:6).EQ.'HBDO ') THEN
C
         IDX1 = COORD_IDX(LABEL(7:7))
         BSRLON(IDX1,ISYMTR,IEXCI) = -TRMOM *EXENG(ISYMTR,IEXCI)
C
CClark:22/12/2015 first,and changed at 17/06/2016
      ELSE IF (LABEL(1:3).EQ.'COS' .OR. LABEL(1:3).EQ.'SIN') THEN
C
         IDX1 = COORD_IDX(LABEL(6:6))
C
         TRGOS(IDX1,IEXCI) =
     &       TRGOS(IDX1,IEXCI) + TRMOM * TRMOM * EXENG(ISYMTR,IEXCI)
C
      ELSE
CClark:17/06/2016
C
         WRITE(LUPRI,'(A,/,2A)') 'ERROR in SO_COLLECT_TM:',
     &         'Label not recognized: ', LABEL
C
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_COLLECT_TM')
      RETURN

      CONTAINS
C
C        A function that converts the letters X,Y and Z to
C        the integers 1,2 and 3
         PURE FUNCTION COORD_IDX(A)
            CHARACTER(LEN=1), INTENT(IN) :: A
            INTEGER :: COORD_IDX
            COORD_IDX = ICHAR(A) - ICHAR('X') + 1
            RETURN
         END FUNCTION
      END
