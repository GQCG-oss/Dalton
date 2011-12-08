C ----------------------------------------------------------------------
C  /* Deck so_fndpos */
      SUBROUTINE SO_FNDPOS(LABEL,FREQ,LU)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Andrea Ligabue, december 2003
C
C     Find solution or residual vector in file and position to replace
C     if exist or to append if doesn't exist
C
C     LABEL      property label
C     FREQ       frequency of the solution/residual
C     LU         file where we have to write the property
C
#include<implicit.h>
#include<priunit.h>
C
      DATA EPSI /1.0D-8/
      CHARACTER*8 LABEL,LSTARS,LB
C
      REWIND(LU)
C
    1 READ(LU,END=10) LSTARS,LB,ISYMTR,LEN,FR
C
      IF (LB .EQ. LABEL) THEN
C
         IF ((FR.GT.(FREQ-EPSI)) .AND. (FR.LT.(FREQ+EPSI))) THEN
C
            WRITE(LUPRI,'(3A,I5)')
     &         ' SO_FNDPOS: found label ',LB,' on unit ',LU
            BACKSPACE LU
            CALL FLUSH(LUPRI)
            RETURN
C
         END IF
C
      END IF
C
      READ(LU,END=10) DUMMY
      GO TO 1
C     
   10 CONTINUE
C
      RETURN
C
      END
