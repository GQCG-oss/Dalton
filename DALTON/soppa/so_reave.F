C
C  /* Deck so_reave */
      SUBROUTINE SO_REAVE(ARRAY,LARRAY,ISYMTR,LABEL,FREQ,LU)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Andrea Ligabue, december 2003
C
C     Read solution or residual vector from file
C     and all the needed information related to that vector
C
C     ARRAY      solution/residual vecotr
C     LARRAY     ARRAY length
C     ISYMTR     symmetry of the property
C     LABEL      property label
C     FREQ       frequency of the solution/residual
C     LU         file where we have to write the property
C
#include<implicit.h>
#include<priunit.h>
C
      DIMENSION ARRAY(LARRAY)
C
      CHARACTER*8 LABEL,LSTARS,LB
C
      DATA EPSI /1.0D-8/
C
      READ(LU,ERR=10,END=20) LSTARS,LB,ISYM,LEN,FR
C
      IF (LB .NE. LABEL) THEN
C
         WRITE(LUPRI,'(3A,I5,/,A,A,A,A)')
     &       ' SO_REAVE: error reading "',LABEL,'" on unit ',LU,
     &       ' found label ',LB,' instead of ',LABEL
         CALL QUIT('SO_REAVE: error reading label')
C
      ELSE IF (ISYM .NE. ISYMTR) THEN
C
         WRITE(LUPRI,'(3A,I5,/,A,I5,A,I5)')
     &       ' SO_REAVE: error reading "',LABEL,'" on unit ',LU,
     &       ' found IRREP ',ISYM,' instead of ',ISYMTR
         CALL QUIT('SO_REAVE: error reading label')
C
      ELSE IF (LEN .NE. LARRAY) THEN
C
         WRITE(LUPRI,'(3A,I5,/,A,I5,A,I5)')
     &       ' SO_REAVE: error reading "',LABEL,'" on unit ',LU,
     &       ' array length is ',LEN,' instead of ',LARRAY
         CALL QUIT('SO_REAVE: error reading label')
C
      ELSE IF ( (FR .LT. (FREQ-EPSI)) .OR. (FR .GT. (FREQ+EPSI)) ) THEN
C
         WRITE(LUPRI,'(3A,I5,/,A,F13.8,A,F13.8)')
     &       ' SO_REAVE: error reading "',LABEL,'" on unit ',LU,
     &       ' found frequency ',FR,' instead of ',FREQ
         CALL QUIT('SO_REAVE: error reading label')
C
      ELSE
C
         READ(LU) ARRAY
C
      END IF
C
      RETURN

C
 10   WRITE(LUPRI,'(3A,I5)')
     &   ' SO_REAVE: error reading "',LABEL,'" on unit ',LU
      CALL QUIT('SO_REAVE: error reading label')
C
 20   WRITE(LUPRI,'(3A,I5)')
     &   ' SO_REAVE: end of file while reading "',LABEL,'" on unit ',LU
      CALL QUIT('SO_REAVE: end of file')
C
      END
