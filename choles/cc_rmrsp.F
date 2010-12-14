C  /* Deck cc_rmrsp */
      SUBROUTINE CC_RMRSP(LIST,IDXLST,ISYM)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Add the possibility of removing files in the
C              CC_RDRSP and CC_WRRSP functionalities.
C
C     Originally written with the purpose of deleting effective rhs
C     vector files in Cholesky CC2.
C
#include <implicit.h>
      CHARACTER*(*) LIST
#include <priunit.h>
#include <ccsdinp.h>
#include <dummy.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_RMRSP')

      CHARACTER*4  LISTI
      CHARACTER*10 FILEX

      LOGICAL LEXIST, LOPEND

C     Get file name.
C     --------------

      CALL CC_RWPRE(LIST,IDXLST,ISYM,LISTI,IDXVEC,FILEX)

C     Delete file it if it exists.
C     ----------------------------

      INQUIRE(FILE=FILEX,EXIST=LEXIST,OPENED=LOPEND,NUMBER=LUSAVE,
     &        IOSTAT=IOS,ERR=990)

      IF (LEXIST) THEN
         IF (.NOT. LOPEND) THEN
            LUSAVE = -1
            CALL GPOPEN(LUSAVE,FILEX,'OLD','SEQUENTIAL','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
         ENDIF
         CALL GPCLOSE(LUSAVE,'DELETE')
      ENDIF

      RETURN

C     Error branches.
C     ---------------

990   CONTINUE
      IF (IPRINT .GT. 0) THEN 
         WRITE (LUPRI,'(2A)')
     &   ' an error occured while inquireing file ',FILEX
      END IF
      GOTO 999

999   CONTINUE
      WRITE (LUPRI,'(A,A)')       ' Fatal I/O error in ',SECNAM
      WRITE (LUPRI,'(A,3A)')      ' LIST / LISTI   :',LIST(1:3),' / ',
     &     LISTI
      WRITE (LUPRI,'(A,2I5)')     ' IDXLST, IDXVEC :',IDXLST,IDXVEC
      WRITE (LUPRI,'(A,I5)')      ' unit number    :',LUSAVE
      WRITE (LUPRI,'(A,I5)')      ' returned IOSTAT:',IOS
      CALL QUIT ('Fatal i/o error in '//SECNAM)

      END
