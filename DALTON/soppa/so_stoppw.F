C
C  /* Deck so_stoppw */
      SUBROUTINE SO_STOPPW(STOPP,TRGOS,ISYMTR,IEXCI,EXENG,QVAL)
C
C     This routine is part of the atomic integral directSOPPA program.
C
C     The charge Z of the incoming ions is set to 1 here.
C
C     Zhiwen Shi (Clark), Stephan P. A. Sauer, January 2016
C
C     PURPOSE: Calculate Stopping Power.
C
#include "implicit.h"
#include "cbiexc.h"
#include "ccorb.h"
#include "pi.h"
C
      DIMENSION STOPP(3,LVEL,2),TRGOS(3),EXENG(NSYM,MXNEXI)
      REAL*8    QVAL,VELOC,QMAXV,QMINV
      PARAMETER ( D4 = 4.0D00 )
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_STOPPW')
C
C--------------------------
C     Loop over velocities.
C--------------------------
C
      DO IVEL=1, LVEL
C
         VELOC = VMIN+(IVEL-1)*VSTEP
C
         QMAXV = VELOC*2
         QMINV = EXENG(ISYMTR,IEXCI)/VELOC
C
         IF (QMINV .LE. QMAXV) THEN
C
            IF ((QVAL .GE. QMINV) .AND. (QVAL .LE. QMAXV)) THEN
C
               STOPP(1,IVEL,1) = STOPP(1,IVEL,1) + TRGOS(1)*QSTEP*
     &                                         PI*D4/(QVAL*VELOC*2)
               STOPP(2,IVEL,1) = STOPP(2,IVEL,1) + TRGOS(2)*QSTEP*
     &                                         PI*D4/(QVAL*VELOC*2)
               STOPP(3,IVEL,1) = STOPP(3,IVEL,1) + TRGOS(3)*QSTEP*
     &                                         PI*D4/(QVAL*VELOC*2)

            ENDIF
C
            IF ((QVAL .GE. QMINV) .AND. (QVAL .LE. QINP)) THEN
C
               STOPP(1,IVEL,2) = STOPP(1,IVEL,2) + TRGOS(1)*QSTEP*
     &                                         PI*D4/(QVAL*VELOC*2)
               STOPP(2,IVEL,2) = STOPP(2,IVEL,2) + TRGOS(2)*QSTEP*
     &                                         PI*D4/(QVAL*VELOC*2)
               STOPP(3,IVEL,2) = STOPP(3,IVEL,2) + TRGOS(3)*QSTEP*
     &                                         PI*D4/(QVAL*VELOC*2)

            ENDIF
C
         ENDIF
C
      END DO
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_STOPPW')
C
      RETURN
C
      END
