*---------------------------------------------------------------------*
c/* Deck WOPEN2 */
*=====================================================================*
      SUBROUTINE WOPEN2(LUNIT,FILE,BLOCKSIZE,STATUS)
*---------------------------------------------------------------------*
* Purpose: open file and test error flag, if not 0, stop
*=====================================================================*
#include <priunit_dec.h>
#include <priunit.h>
#include <ccsdinp.h>
      INTEGER LUNIT, IADR, ILEN, IERR, BLOCKSIZE, STATUS, IUN
      INTEGER IUNTAB, J
      CHARACTER*(*) FILE
      COMMON /UNITAR/ IUNTAB(1:99)
C
      ILEN=LEN(FILE)
      DO J = LEN(FILE), 1, -1
        IF (FILE(J:J).EQ.' ') ILEN = MIN(ILEN,J)
      END DO
C
      IF (LUNIT .EQ. 0) THEN
         IUN = 0
 10      CONTINUE
         IUN = IUN + 1
         IF ((IUN .EQ. 5) .OR. (IUN .EQ. 6))
     &        GOTO 10
         IF (IUNTAB(IUN) .NE. 0) GOTO 10
         LUNIT = IUN
         IF (LUNIT .GT. 99) GOTO 9001
         IUNTAB(LUNIT) = 1
      END IF
C
      IF (DEBUG) WRITE(LUPRI,*) 'OPEN FILE "',FILE(1:ILEN),
     & '" with unit number ',LUNIT
C
#if defined (SYS_AIX)
      CALL WOPEN_(LUNIT,FILE,ILEN,BLOCKSIZE,STATUS,IERR)
#endif
#if defined (SYS_CRAY)
      CALL WOPEN(FILE,ILEN,BLOCKSIZE,STATUS,IERR)
#endif
#if !defined (SYS_AIX) && !defined (SYS_CRAY)
      CALL WOPEN(LUNIT,FILE,ILEN,BLOCKSIZE,STATUS,IERR)
#endif
C
      IF (IERR.NE.0) THEN
        WRITE (LUPRI,*) 'I/O ERROR IN WOPEN2:'
        WRITE (LUPRI,'(1x,3a)') 'FILE: >', FILE(1:ILEN),'<'
        WRITE (LUPRI,*) 'ILEN:', ILEN
        WRITE (LUPRI,*) 'UNIT:', LUNIT
        WRITE (LUPRI,*) 'IERR:', IERR
        WRITE (LUPRI,*) 'STATUS:', STATUS
        WRITE (LUPRI,*) 'BLOCKSIZE:', BLOCKSIZE
        CALL QUIT( 'I/O ERROR IN WOPEN2')
      END IF
      RETURN
C     
C error branch
C
 9001 CONTINUE
      WRITE (LUPRI,'(//A/A/A)')
     &     ' *** ERROR (WOPEN2) NO MORE AVAILABLE FILENUMBERS! ',
     &     ' *** THIS CALCULATION EITHER NEEDS TOO MANY SIMULTANEOUS '//
     &     'FILES OR', ' *** SOMEBODY HAS FORGOT TO CLOSE FILES IN '//
     &     'THE SOURCE CODE'
      WRITE (LUPRI,'(/A)') ' ### Please report the problem to '//
     &     'dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (WOPEN2) NO MORE FILE NUMBERS')
C
      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE WOPEN2
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck PUTWA2 */
*=====================================================================*
      SUBROUTINE PUTWA2(LUNIT,FILE,DATA,IADR,LEN)
*---------------------------------------------------------------------*
* Purpose: call putwa and test error flag, if not 0, stop
*=====================================================================*
#include <priunit_dec.h>
#include <priunit.h>
      INTEGER LUNIT, IADR, LEN, IERR
      CHARACTER*(*) FILE
      DIMENSION DATA(LEN)

#if defined (SYS_AIX)
      CALL PUTWA_(LUNIT,DATA(1),IADR,LEN,IERR)
#endif
#if defined (SYS_CRAY)
      CALL PUTWA(FILE,DATA(1),IADR,LEN,IERR)
#endif
#if !defined (SYS_AIX) && !defined (SYS_CRAY)
      CALL PUTWA(LUNIT,DATA(1),IADR,LEN,IERR)
#endif

      IF (IERR.NE.0) THEN
        WRITE (LUPRI,*) 'I/O ERROR IN PUTWA2:'
        WRITE (LUPRI,*) 'FILE:', FILE
        WRITE (LUPRI,*) 'UNIT:', LUNIT
        CALL QUIT('I/O ERROR IN PUTWA2.')
      END IF

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE PUTWA2
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck GETWA2 */
*=====================================================================*
      SUBROUTINE GETWA2(LUNIT,FILE,DATA,IADR,LEN)
*---------------------------------------------------------------------*
* Purpose: call GETWA and test error flag, if not 0, stop
*=====================================================================*
#include <priunit_dec.h>
#include <priunit.h>
      INTEGER LUNIT, IADR, LEN, IERR
      CHARACTER*(*) FILE
      DIMENSION DATA(LEN)

#if defined (SYS_AIX)
      CALL GETWA_(LUNIT,DATA(1),IADR,LEN,IERR)
#endif
#if defined (SYS_CRAY)
      CALL GETWA(FILE,DATA(1),IADR,LEN,IERR)
#endif
#if !defined (SYS_AIX) && !defined (SYS_CRAY)
      CALL GETWA(LUNIT,DATA(1),IADR,LEN,IERR)
#endif

      IF (IERR.NE.0) THEN
        WRITE (LUPRI,*) 'I/O ERROR IN GETWA2:'
        WRITE (LUPRI,*) 'FILE:', FILE
        WRITE (LUPRI,*) 'UNIT:', LUNIT
        WRITE (LUPRI,*) 'IADR:', IADR 
        WRITE (LUPRI,*) 'LEN :', LEN 
        WRITE (LUPRI,*) 'IERR:', IERR 
        CALL QUIT('I/O ERROR IN GETWA2.')
      END IF

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE GETWA2
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck WCLOSE2 */
*=====================================================================*
      SUBROUTINE WCLOSE2(LUNIT,FILE,MODE)
*---------------------------------------------------------------------*
* Purpose: close file, test error flag, if not 0 bump...
*          if MODE='DELETE' remove the file, else keep it
*=====================================================================*
#include <priunit_dec.h>
#include <priunit.h>
#include <ccsdinp.h>
      INTEGER LUNIT, IERR, ILEN
      CHARACTER*(*) FILE, MODE
      CHARACTER*(80) COMMAND
      COMMON /UNITAR/ IUNTAB(1:99)
C
#if defined (SYS_AIX)
      CALL WCLOSE_(LUNIT,IERR)
#endif
#if defined (SYS_CRAY)
      CALL WCLOSE(FILE,IERR)
#endif
#if !defined (SYS_AIX) && !defined (SYS_CRAY)
      CALL WCLOSE(LUNIT,IERR)
#endif
C
      IF (IERR.NE.0) THEN
        WRITE (LUPRI,*) 'I/O ERROR IN CLOSE2:'
        WRITE (LUPRI,*) 'FILE:', FILE
        WRITE (LUPRI,*) 'UNIT:', LUNIT
        WRITE (LUPRI,*) 'MODE:', MODE
        CALL QUIT('I/O ERROR IN CLOSE2')
      END IF

      ILEN = MIN(LEN(FILE),74)

      IF (DEBUG) WRITE(LUPRI,*) 'CLOSE FILE "',FILE(1:ILEN),
     & '" with unit number ',LUNIT
C
      IF (MODE(1:6).EQ.'DELETE') THEN

        WRITE(COMMAND,'(2A)') 'rm -f ', FILE(1:ILEN)

#if defined (SYS_CRAY)
        CALL INFO = ISHELL(COMMAND)
#else
        CALL SYSTEM(COMMAND)
#endif

      END IF
      IUNTAB(LUNIT) = 0
      LUNIT = 0
      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE WCLOSE2
*=====================================================================*
