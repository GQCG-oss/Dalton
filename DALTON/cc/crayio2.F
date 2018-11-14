!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
*---------------------------------------------------------------------*
c/* Deck WOPEN2 */
*=====================================================================*
      SUBROUTINE WOPEN2(LUNIT,FILE,BLOCKSIZE,STATUS)
*---------------------------------------------------------------------*
* Purpose: open file and test error flag, if not 0, stop
*=====================================================================*
#include "priunit.h"
#include "ccsdinp.h"
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
#include "priunit.h"
      INTEGER LUNIT, IADR, LEN, IERR
      INTEGER LMAX, JADR, NBATCH, LETOT, IBATCH, IOFF
      INTEGER LREAD
      CHARACTER*(*) FILE
      DOUBLE PRECISION DATA(LEN)
c     DIMENSION DATA(LEN)
C
Casm Apparently, it is not possible to read more than 2 Gb (268435456 dw)
C        
      LMAX = 250000000      ! Just to be in the (too?) safe side.
c     lmax = 80000000      ! Just to be in the (too?) safe side.
C
      LETOT = 0
      LREAD = LMAX
      JADR  = IADR
      IOFF  = 1
C
      NBATCH = (LEN - 1) / LMAX + 1
      DO IBATCH = 1, NBATCH
         LETOT = LETOT + LREAD
         IF (LETOT .GT. LEN) THEN
            LETOT = LETOT - LREAD
            LREAD = LEN   - LETOT
         END IF
C
#if defined (SYS_CRAY)
C        CALL PUTWA(FILE,DATA(1),IADR,LEN,IERR)
         CALL PUTWA(FILE,DATA(IOFF),JADR,LREAD,IERR)
#else
C        CALL PUTWA(LUNIT,DATA(1),IADR,LEN,IERR)
         CALL PUTWA(LUNIT,DATA(IOFF),JADR,LREAD,IERR)
#endif
C
         IF (IERR.NE.0) THEN
           WRITE (LUPRI,*) 'I/O ERROR IN PUTWA2:'
           WRITE (LUPRI,*) 'FILE:', FILE
           WRITE (LUPRI,*) 'UNIT:', LUNIT
           WRITE (LUPRI,*) 'IADR:', IADR
           WRITE (LUPRI,*) 'LEN:',  LEN
           CALL QUIT('I/O ERROR IN PUTWA2.')
         END IF
C
         IOFF = IOFF + LREAD
         JADR = JADR + LREAD
C
      END DO

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
#include "priunit.h"
      INTEGER LUNIT, IADR, LEN, IERR
      INTEGER LMAX, JADR, NBATCH, LETOT, IBATCH, IOFF
      INTEGER LREAD
      CHARACTER*(*) FILE
      DOUBLE PRECISION DATA(LEN)
c     DIMENSION DATA(LEN)
C
      
C
Casm Apparently, it is not possible to read more than 2 Gb (268435456 dw)
C        
      LMAX = 250000000      ! Just to be in the (too?) safe side.
c     lmax = 80000000      ! Just to be in the (too?) safe side.
C
      LETOT = 0
      LREAD = LMAX
      JADR  = IADR
      IOFF  = 1
C
      NBATCH = (LEN - 1) / LMAX + 1
      DO IBATCH = 1, NBATCH
         LETOT = LETOT + LREAD
         IF (LETOT .GT. LEN) THEN
            LETOT = LETOT - LREAD
            LREAD = LEN   - LETOT
         END IF
C
#if defined (SYS_CRAY)
C        CALL GETWA(FILE,DATA(1),IADR,LEN,IERR)
         CALL GETWA(FILE,DATA(IOFF),JADR,LREAD,IERR)
#else
C        CALL GETWA(LUNIT,DATA(1),IADR,LEN,IERR)
         CALL GETWA(LUNIT,DATA(IOFF),JADR,LREAD,IERR)
#endif
C
         IF (IERR.NE.0) THEN
           WRITE (LUPRI,*) 'I/O ERROR IN GETWA2:'
           WRITE (LUPRI,*) 'FILE:', FILE
           WRITE (LUPRI,*) 'UNIT:', LUNIT
           WRITE (LUPRI,*) 'IADR:', IADR 
           WRITE (LUPRI,*) 'LEN :', LEN 
           WRITE (LUPRI,*) 'IERR:', IERR 
           CALL QUIT('I/O ERROR IN GETWA2.')
         END IF
C
         IOFF = IOFF + LREAD
         JADR = JADR + LREAD
C
      END DO
C
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
#include "priunit.h"
#include "ccsdinp.h"
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
#include "priunit.h"
#include "ccsdinp.h"
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
