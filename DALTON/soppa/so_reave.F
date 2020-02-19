C
C  /* Deck so_reave */
      SUBROUTINE SO_REAVE(ARRAY,LARRAY,ISYMTR,LABEL,FREQ,LU)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Andrea Ligabue, december 2003
C     Rasmus Faber, 2017 : Rewrite using linked lists
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
C#include<implicit.h>
C
      USE SO_DATA, ONLY : FILEINF
      USE SO_INFO, ONLY : SOP_DP
C
      IMPLICIT NONE
C
      REAL(SOP_DP), INTENT(OUT) :: ARRAY(LARRAY)
      REAL(SOP_DP), INTENT(IN) :: FREQ
      INTEGER, INTENT(IN) :: LARRAY, ISYMTR, LU
      CHARACTER(LEN=8) :: LABEL
C
      INTEGER :: POS
#include<priunit.h>
C
      POS = FILEINF%POS(LABEL,FREQ,ISYMTR)
C      
      IF (POS.EQ.0) GOTO 20

      READ(LU,err=10,rec=pos) ARRAY
C
      RETURN
C
 10   WRITE(LUPRI,'(3A,I5)')
     &   ' SO_REAVE: error reading "',LABEL,'" on unit ',LU
      CALL QUIT('SO_REAVE: error reading label')
C
 20   WRITE(LUPRI,'(3A,I5)')
     &   ' SO_REAVE: LABEL "',LABEL,'" not fount on unit ',LU
      CALL QUIT('SO_REAVE: end of file')
C
      END
