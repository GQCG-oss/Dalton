C
C /* extra routines only needed for standalone version of MM code */
C
      SUBROUTINE GPOPEN(LUNIT,FILEIN,STATIN,ACCEIN,FORMIN,NELEM,OLDDX)
#include <implicit.h>
      CHARACTER(*) FILEIN, STATIN, ACCEIN, FORMIN
      INTEGER LUNIT, NELEM, PREV
      LOGICAL OLDDX
      SAVE PREV
      DATA PREV/10/
      LUNIT = PREV
      PREV = PREV + 1
      OPEN(LUNIT,FILE=FILEIN,FORM=FORMIN)
      END
C
      SUBROUTINE GPCLOSE(LUNIT,STATIN)
#include <implicit.h>
      CHARACTER(*) STATIN
      CLOSE(LUNIT,STATUS='KEEP')
      END
C
      SUBROUTINE TIMTXT(TEXT,TIMUSD,LUPRIN)
C
C TIMTXT based on TIMER by TUH //900709-hjaaj
C
#include <implicit.h>
      CHARACTER(*) TEXT
      CHARACTER AHOUR*6, ASEC*8, AMIN*8
C
      ISECND = NINT(TIMUSD)
      IF (ISECND .GE. 60) THEN
         MINUTE = ISECND/60
         IHOURS = MINUTE/60
         MINUTE = MINUTE - 60*IHOURS
         ISECND = ISECND - 3600*IHOURS - 60*MINUTE
         IF (IHOURS .EQ. 1) THEN
            AHOUR = ' hour '
         ELSE
            AHOUR = ' hours'
         END IF
         IF (MINUTE .EQ. 1) THEN
            AMIN = ' minute '
         ELSE
            AMIN = ' minutes'
         END IF
         IF (ISECND .EQ. 1) THEN
            ASEC = ' second '
         ELSE
            ASEC = ' seconds'
         END IF
         IF (IHOURS .GT. 0) THEN
            WRITE(LUPRIN,100)
     *            TEXT, IHOURS, AHOUR, MINUTE, AMIN, ISECND, ASEC
         ELSE
            WRITE(LUPRIN,200) TEXT, MINUTE, AMIN, ISECND, ASEC
         END IF
      ELSE
         WRITE(LUPRIN,300) TEXT,TIMUSD
      END IF
  100 FORMAT(1X,A,I4,A,I3,A,I3,A)
  200 FORMAT(1X,A,     I3,A,I3,A)
  300 FORMAT(1X,A,F7.2,' seconds')
      RETURN
      END

      SUBROUTINE QUIT(msg)
        CHARACTER(*) msg
        WRITE(LUPRI,*) msg
        print *, msg
        STOP ">>> FATAL ERROR"
      END SUBROUTINE QUIT

