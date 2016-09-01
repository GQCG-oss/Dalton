C
C  /* Deck so_stoppw */
      SUBROUTINE SO_STOPPW(STOPP,TRGOS,ISYMTR,IEXCI,EXENG,QVAL)
C
C     This routine is part of the atomic integral directSOPPA program.
C
C     Zhiwen Shi (Clark), Stephan P. A. Sauer, January 2016
C
C     PURPOSE: Calculate Stopping Power.
C
#include "implicit.h"
#include "cbiexc.h"
#include "ccorb.h"
C
      DIMENSION STOPP(3,LVEL),TRGOS(3),EXENG(NSYM,MXNEXI)
      REAL*8    QVAL,VELOC,QMAXV,QMINV
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
               STOPP(1,IVEL) = STOPP(1,IVEL) + TRGOS(1)*QSTEP
               STOPP(2,IVEL) = STOPP(2,IVEL) + TRGOS(2)*QSTEP
               STOPP(3,IVEL) = STOPP(3,IVEL) + TRGOS(3)*QSTEP
C
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
