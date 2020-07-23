C
      SUBROUTINE SO_STOPPW(RSTOPP,TRGOS,ISYMTR,IEXCI,EXENG,QVAL)
C
C     This routine is part of the atomic integral directSOPPA program.
C     The charge Z of the incoming ions is set to 1 here.
C     Zhiwen Shi (Clark), Stephan P. A. Sauer, January 2016
C     PURPOSE: Calculate Stopping Power.
C
      implicit none
#include "cbiexc.h" ! LVEL, MXNEXI
#include "ccorb.h" ! NSYM
#include "pi.h" ! PI
C
      REAL*8, INTENT(INOUT) :: RSTOPP(3,LVEL,2)
      REAL*8, INTENT(IN) :: TRGOS(3), EXENG(NSYM,MXNEXI), QVAL
      INTEGER, INTENT(IN) :: ISYMTR, IEXCI
      REAL*8 :: VELOC, QMAXV, QMINV
      REAL*8, PARAMETER :: D4 = 4.0D0
      INTEGER :: IVEL
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
         QMAXV = VELOC*2
         QMINV = EXENG(ISYMTR,IEXCI)/VELOC
C
         IF (QMINV .LE. QMAXV) THEN
C
            IF ((QVAL .GE. QMINV) .AND. (QVAL .LE. QMAXV)) THEN
C
               RSTOPP(1,IVEL,1) = RSTOPP(1,IVEL,1) + TRGOS(1)*QSTEP*
     &                           PI*D4/(QVAL*VELOC*VELOC)
               RSTOPP(2,IVEL,1) = RSTOPP(2,IVEL,1) + TRGOS(2)*QSTEP*
     &                           PI*D4/(QVAL*VELOC*VELOC)
               RSTOPP(3,IVEL,1) = RSTOPP(3,IVEL,1) + TRGOS(3)*QSTEP*
     &                           PI*D4/(QVAL*VELOC*VELOC)

            ENDIF
C
C The next line makes sure that the split integration will be applied
C after the velocity is larger than half of Q value.
C This velocity as half of chosen Q value corresponds to
C the highest excitation energy for given basis set.
C i.e. integration can be split after this velocity.
C
            IF (VELOC .GE. QINP/2.0d0) THEN
C
               IF ((QVAL .GE. QMINV) .AND. (QVAL .LE. QINP)
     &                            .AND. (QVAL .LE. QMAXV)) THEN
C
                  RSTOPP(1,IVEL,2) = RSTOPP(1,IVEL,2) + TRGOS(1)*QSTEP*
     &                              PI*D4/(QVAL*VELOC*VELOC)
                  RSTOPP(2,IVEL,2) = RSTOPP(2,IVEL,2) + TRGOS(2)*QSTEP*
     &                              PI*D4/(QVAL*VELOC*VELOC)
                  RSTOPP(3,IVEL,2) = RSTOPP(3,IVEL,2) + TRGOS(3)*QSTEP*
     &                              PI*D4/(QVAL*VELOC*VELOC)
C
               ENDIF
C
            ELSEIF ((QVAL .GE. QMINV) .AND. (QVAL .LE. QMAXV)) THEN
C
               RSTOPP(1,IVEL,2) = RSTOPP(1,IVEL,2) + TRGOS(1)*QSTEP*
     &                           PI*D4/(QVAL*VELOC*VELOC)
               RSTOPP(2,IVEL,2) = RSTOPP(2,IVEL,2) + TRGOS(2)*QSTEP*
     &                           PI*D4/(QVAL*VELOC*VELOC)
               RSTOPP(3,IVEL,2) = RSTOPP(3,IVEL,2) + TRGOS(3)*QSTEP*
     &                           PI*D4/(QVAL*VELOC*VELOC)
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
