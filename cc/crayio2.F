!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
*---------------------------------------------------------------------*
c/* Deck WOPEN2 */
*=====================================================================*
      SUBROUTINE WOPEN2(LUNIT,FILE,BLOCKSIZE,STATUS)
*---------------------------------------------------------------------*
* Purpose: open file and test error flag, if not 0, stop
*=====================================================================*
#include <priunit.h>
#include <ccsdinp.h>
      INTEGER LUNIT, IADR, ILEN, IERR, BLOCKSIZE, STATUS, IUN
      INTEGER IUNCRAY, J
      CHARACTER*(*) FILE
      COMMON /UNICRAYIO/ IUNCRAY(1:99)
C
      ILEN=LNBLNK(FILE)
C
      IF (LUNIT .LE. 0) THEN
         IUN = 0
 10      CONTINUE
         IUN = IUN + 1
         IF (IUN .GT. 99) GOTO 9001
         IF (IUNCRAY(IUN) .NE. 0) GOTO 10
         LUNIT = IUN
         IUNCRAY(LUNIT) = 1
      END IF
C
      IF (DEBUG) WRITE(LUPRI,*) 'WOPEN2 : open file "',FILE(1:ILEN),
     & '" with unit number ',LUNIT
C
#if defined (SYS_CRAY)
      CALL WOPEN(FILE,ILEN,BLOCKSIZE,STATUS,IERR)
#else
      CALL WOPEN(LUNIT,FILE,ILEN,BLOCKSIZE,STATUS,IERR)
#endif
C
      IF (IERR.NE.0) THEN
        WRITE (LUPRI,*) 'I/O ERROR IN WOPEN2:'
        WRITE (LUPRI,'(3a)') ' FILE: >', FILE(1:ILEN),'<'
        WRITE (LUPRI,*) 'ILEN      : ', ILEN
        WRITE (LUPRI,*) 'UNIT      : ', LUNIT
        WRITE (LUPRI,*) 'IERR      : ', IERR
        WRITE (LUPRI,*) 'STATUS    : ', STATUS
        WRITE (LUPRI,*) 'BLOCKSIZE : ', BLOCKSIZE
        CALL QUIT( 'I/O ERROR IN WOPEN2')
      END IF
      RETURN
C     
C error branch
C
 9001 CONTINUE
      WRITE (LUPRI,'(//A//A/A//A)')
     &     ' *** ERROR (WOPEN2) NO MORE FILENUMBERS AVAILABLE!',
     &     ' *** THIS CALCULATION EITHER NEEDS TOO MANY SIMULTANEOUS '//
     &     'FILES OR',
     &     ' *** SOMEBODY HAS FORGOT TO CLOSE FILES IN THE SOURCE CODE',
     &     ' ### Please report the problem to dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (WOPEN2) NO MORE FILE NUMBERS AVAILABLE')
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
#include <priunit.h>
      INTEGER LUNIT, IADR, LEN, IERR
      CHARACTER*(*) FILE
      DIMENSION DATA(LEN)

#if defined (SYS_CRAY)
      CALL PUTWA(FILE,DATA(1),IADR,LEN,IERR)
#else
      CALL PUTWA(LUNIT,DATA(1),IADR,LEN,IERR)
#endif

      IF (IERR.NE.0) THEN
        WRITE (LUPRI,*) 'I/O ERROR IN PUTWA2:'
        WRITE (LUPRI,*) 'FILE:', FILE
        WRITE (LUPRI,*) 'UNIT:', LUNIT
        WRITE (LUPRI,*) 'IADR:', IADR
        WRITE (LUPRI,*) 'LEN:',  LEN
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
#include <priunit.h>
      INTEGER LUNIT, IADR, LEN, IERR
      CHARACTER*(*) FILE
      DIMENSION DATA(LEN)

#if defined (SYS_CRAY)
      CALL GETWA(FILE,DATA(1),IADR,LEN,IERR)
#else
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
#include <priunit.h>
#include <ccsdinp.h>
      INTEGER LUNIT, IERR, ILEN
      CHARACTER*(*) FILE, MODE
      CHARACTER*(80) COMMAND
      COMMON /UNICRAYIO/ IUNCRAY(1:99)
C
#if defined (SYS_CRAY)
      CALL WCLOSE(FILE,IERR)
#else
      IF (LUNIT .LT. 1 .OR. LUNIT .GT. 99) THEN
         WRITE(LUPRI,*) 'WCLOSE2 called with illegal unit number',LUNIT
         IERR = -1
      ELSE
         CALL WCLOSE(LUNIT,IERR)
      END IF
#endif
C
      IF (IERR.NE.0) THEN
        WRITE (LUPRI,*) 'I/O ERROR IN CLOSE2:'
        WRITE (LUPRI,*) 'FILE:', FILE
        WRITE (LUPRI,*) 'UNIT:', LUNIT
        WRITE (LUPRI,*) 'MODE:', MODE
        WRITE (LUPRI,*) 'IERR:', IERR
        CALL QUIT('I/O ERROR IN WCLOSE2')
      END IF

      ILEN = LNBLNK(FILE)

      IF (DEBUG) WRITE(LUPRI,*) 'CLOSE FILE "',FILE(1:ILEN),
     & '" with unit number ',LUNIT
C
      L = MIN(6,LEN(MODE))
      IF (MODE(1:L).EQ.'DELETE') THEN

        IF (ILEN .LE. 74) THEN
           WRITE(COMMAND,'(2A)') 'rm -f ', FILE(1:ILEN)

#if defined (SYS_CRAY)
           CALL INFO = ISHELL(COMMAND)
#else
           CALL SYSTEM(COMMAND)
#endif
        ELSE
           NWARN = NWARN + 1
           WRITE(LUPRI,*) 'WARNING! "rm -f ',FILE(1:ILEN),'" failed!'
           WRITE(LUPRI,*) '(increase length of COMMAND in WCLOSE2)'
        END IF

      END IF
      IF (LUNIT .LE. 70) THEN
         IUNCRAY(LUNIT) = 0
         LUNIT = -1
      END IF
      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE WCLOSE2
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck INITWIO */
*=====================================================================*
      SUBROUTINE INITWIO()
*---------------------------------------------------------------------*
* Purpose: initialize common /UNICRAYIO/ (all files marked as unused)
*          reserve units .gt. 70 for Cholesky
*=====================================================================*
      INTEGER IUNCRAY
      COMMON /UNICRAYIO/ IUNCRAY(1:99)
 
      DO J = 1, 70
        IUNCRAY(J) = 0
      END DO
 
      DO J = 71, 99
        IUNCRAY(J) = 1
      END DO
 
      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE INITWIO
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck WCLOSEALL */
*=====================================================================*
      SUBROUTINE WCLOSEALL()
*---------------------------------------------------------------------*
* Purpose: close all files opened by WOPEN2 which have not been closed
*=====================================================================*
#include <priunit.h>
#include <ccsdinp.h>
      COMMON /UNICRAYIO/ IUNCRAY(1:99)
C
      DO LUNIT = 1, 70

         IF (IUNCRAY(LUNIT) .NE. 0) THEN

           IF (DEBUG) WRITE(LUPRI,*) 'WCLOSEALL> close unit ',LUNIT

#          if defined (SYS_CRAY)
             CONTINUE
#          else
             CALL WCLOSE(LUNIT,IERR)
#          endif
C
           IF (IERR.NE.0) THEN
             WRITE (LUPRI,*) 'I/O ERROR IN WCLOSEALL:'
c            WRITE (LUPRI,*) 'FILE:', FILE
             WRITE (LUPRI,*) 'UNIT:', LUNIT
c            WRITE (LUPRI,*) 'MODE:', MODE
             CALL QUIT('I/O ERROR IN WCLOSEALL')
           END IF
 
           IF (LUNIT .LE. 70) IUNCRAY(LUNIT) = 0

         END IF

      END DO

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE WCLOSEALL
*=====================================================================*
