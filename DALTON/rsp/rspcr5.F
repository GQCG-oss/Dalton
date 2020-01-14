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
      SUBROUTINE WRTRSP(LU,LEN,VEC,LAB1,LAB2,FREQ1,FREQ2,
     &                  ISYM1,ISYM2,RSD,ANTSYM)
C
#include "implicit.h"
#include "priunit.h"
#include "infrsp.h"
#include "infopt.h"
#include "inforb.h"
C
      CHARACTER*8 LAB1,LAB2,BLANK,LABEL
      PARAMETER   (BLANK='        ', D0 = 0.0D0)
      DIMENSION   VEC(*)
C
      CALL QENTER('WRTRSP')

      REWIND (LU)
 100  READ(LU, END=9000) LABEL
      IF (LABEL.NE.'EOFLABEL') GOTO 100
      BACKSPACE LU
      IF (FREQ2 .EQ. D0) THEN
         WRITE(LU) LAB1,LAB2,ABS(FREQ1),ABS(FREQ2),ISYM1,ISYM2,
     &             ANTSYM,RSD,LEN,EMCSCF,NBAST,NORBT
         IF (FREQ1 .LT. D0) THEN
            CALL DSWAP(LEN/2,VEC,1,VEC(1+LEN/2),1)
            IF (ANTSYM .LT. D0) CALL DSCAL(LEN,ANTSYM,VEC,1)
         END IF
      ELSE
         WRITE(LU) LAB1,LAB2,FREQ1,FREQ2,ISYM1,ISYM2,
     &             ANTSYM,RSD,LEN,EMCSCF,NBAST,NORBT
      END IF
      IF (LEN .GT. 0) CALL WRITT(LU,LEN,VEC)
      WRITE(LU) 'EOFLABEL'
C
      CALL FLSHFO(LU)
C
      CALL QEXIT('WRTRSP')
      RETURN
 9000 CONTINUE
      CALL QUIT('FATAL ERROR:'//
     &   ' "EOFLABEL" not found on file RSPVEC as expected')
      END
      SUBROUTINE REARSP(LU,LEN,VEC,LAB1,LAB2,FREQ1,FREQ2,ISYM1,ISYM2,
     &                  RSD,FOUND,CONV,ANTSYM)
C
C Check if a linear equation needs to be solved or if a solution vector
C already exists on file. FOUND indicates the result.
C CONV indicates if the vector is converged better than RSD.
C
#include "implicit.h"
#include "thrzer.h"
#include "mxcent.h"
#include "priunit.h"
#include "infopt.h"
#include "inforb.h"
#include "gnrinf.h"
#include "optinf.h"
C
      LOGICAL FOUND, CONV, USOLVD
      CHARACTER*8 LABX1,LABX2,LAB1,LAB2,BLANK
      PARAMETER(THR = 1.0D-8, BLANK='        ', D0 = 0.0D0)
      CHARACTER*20 FILE_NAME
      INTEGER INFORB_TEST(33)

C
      CALL QENTER('REARSP')
      FOUND  = .FALSE.
      CONV   = .FALSE.
      USOLVD = .FALSE.
      REWIND (LU)
      ! rec. no. 1 = wave function info record used when opening
      READ (LU,END=200,ERR=200) INFORB_TEST
      INFORB_TEST( 1: 8) = INFORB_TEST( 1: 8) - NISH(:)
      INFORB_TEST( 9:16) = INFORB_TEST( 9:16) - NASH(:)
      INFORB_TEST(17:24) = INFORB_TEST(17:24) - NORB(:)
      INFORB_TEST(25:32) = INFORB_TEST(25:32) - NBAS(:)
      INFORB_TEST(33)    = INFORB_TEST(33)    - NSYM
      NERR = 0
      DO I = 1,33
         IF (INFORB_TEST(I) .NE. 0) NERR = NERR + 1
      END DO

      IF (NERR .NE. 0) THEN
         REWIND (LU)
         READ (LU)
         READ (LU) INFORB_TEST
         WRITE (LUPRI,'(/A/A/,(10X,8I5))')
     &      ' INFO: Orbital data found on response'//
     &      ' vector file does not match current orbital data',
     &      ' Number of inactive, active, orbitals, basis functions'//
     &      ' and number of symmetries:',
     &      INFORB_TEST(1:33)
         GOTO 198
      END IF

C
  100 READ(LU,END=200,ERR=200)
     &      LABX1,LABX2,FREQX1,FREQX2,ISYMX1,ISYMX2,ANTSYM,
     &      RSDX,LEN,EMCX,NBSX,NRBX
      IF (LABX1 .EQ. 'EOFLABEL') GOTO 200
C
      IF ((NRBX .NE. NORBT) .OR. (NBSX .NE. NBAST)) THEN
         WRITE (LUPRI,'(/A,2(/A,2I8))')
     &        ' INFO: Orbital data found on response'//
     &        ' vector file does not match current orbital data',
     &        ' Number of basis functions/orbitals on file      :',
     &        NBSX,NRBX,
     &        ' Number of basis functions/orbitals in this calc.:',
     &        NBAST,NORBT
         WRITE (LUPRI,'(/A/A,T31,A,T51,A/A,T21,2F20.5/A,T21,2I20
     &    /A,2F10.2,I10/A,F30.12/A,2I20)')
     &   'Dump of info in this record:',
     &   ' Labels',LABX1,LABX2,
     &   ' Frequencies',FREQX1,FREQX2,
     &   ' Symmetries',ISYMX1,ISYMX2,
     &   ' Antsym, residual, record length',ANTSYM,RSDX,LEN,
     &   ' Energy',EMCX,
     &   ' Nr. of basis functions and orbitals',NBSX,NRSX
         GOTO 198
      ELSE
         IF (ABS((EMCX - EMCSCF)/EMCSCF) .GT. 100D0*THRZER) THEN
Chj         IF ((OPTNEW .AND. ITRNMR .GT. 0) .OR.
Chj  &          (OPTWLK .AND. ITERNR .GT. 0)) THEN
Chj         Try to use as initial guess /hjaaj-Mar2006
               USOLVD = .TRUE.
Chj         ELSE
            IF (ABS((EMCX - EMCSCF)/EMCSCF) .GT. 1.D-3) THEN
Chj         If they don't agree to three digits, then probably
Chj         from a different calculation /hjaaj-Mar2006
               WRITE (LUPRI,'(/A,2(/A,F24.14))')
     &              ' Molecular energy found on response'//
     &              ' vector file does not match current energy',
     &              ' Energy on file            : ',EMCX,
     &              ' Energy in this calculation: ',EMCSCF
               GOTO 198
            END IF
         END IF
      END IF
      IF (LAB2 .EQ. '        ') THEN
         IF (((ABS(ABS(FREQ1)-FREQX1)).LE.THR)
     &        .AND. (LAB1.EQ.LABX1) .AND. (LAB2.EQ.LABX2) .AND.
     &        (ISYM1.EQ.ISYMX1)) THEN
            FOUND = .TRUE.
            IF ((RSDX-RSD) .LT. D0) THEN
               CONV = .TRUE.
            ELSE
               WRITE(LUPRI,'(2(/A,1P,D10.2))')
     &         ' WARNING: Solution vector on RSPVEC converged to',RSDX,
     &         ' WARNING: which is less than desired convergence',RSD

               !Write(lupri,*) 'RSPVEC: record length is',LEN
            END IF
            IF (LEN .GT. 0) CALL READT(LU,LEN,VEC)
            GOTO 199
         ELSE
            IF (LEN .GT. 0) READ (LU)
            GO TO 100
         END IF
      ELSE
         IF ((      (ABS(FREQ1-FREQX1).LE.THR)
     &        .AND. (ABS(FREQ2-FREQX2).LE.THR)
     &        .AND. LAB1.EQ.LABX1 .AND. LAB2.EQ.LABX2 .AND.
     &        ISYM1.EQ.ISYMX1 .AND. ISYM2.EQ.ISYMX2) .OR.
     &        (     (ABS(FREQ1-FREQX2).LE.THR)
     &        .AND. (ABS(FREQ2-FREQX1).LE.THR)
     &        .AND. LAB1.EQ.LABX2 .AND. LAB2.EQ.LABX1 .AND.
     &        ISYM1.EQ.ISYMX2 .AND. ISYM2.EQ.ISYMX1)) THEN
            FOUND = .TRUE.
            IF ((RSDX-RSD) .LT. D0) THEN
               CONV = .TRUE.
            ELSE
               WRITE(LUPRI,'(2(/A,1P,D10.2))')
     &         ' WARNING: Solution vector on RSPVEC converged to',RSDX,
     &         ' WARNING: which is less than desired convergence',RSD

               !write(lupri,*) 'RSPVEC: record length is',LEN
            END IF
            IF (LEN .GT. 0) CALL READT(LU,LEN,VEC)
            GOTO 199
         ELSE
            IF (LEN .GT. 0) READ (LU)
            GO TO 100
         END IF
      END IF
C
 198  CONTINUE
         INQUIRE(UNIT=LU,NAME=FILE_NAME)
         LFN = LNBLNK(FILE_NAME)
         WRITE (LUPRI,'(/3A)')
     &      ' INFO: ',FILE_NAME(1:LFN),' is reset to empty file'
         REWIND (LU)
         WRITE (LU) NISH,NASH,NORB,NBAS,NSYM
         WRITE (LU) 'EOFLABEL'
      GO TO 200
C
 199  CONTINUE
      IF (USOLVD) THEN
         CONV = .FALSE.
         WRITE (LUPRI,'(/A/10X,A,3X,A)')
     &      ' INFO: Converged vectors from'//
     &      ' previous calculation are used as starting '//
     &      'vectors for property',LABX1, LABX2
      END IF
C
 200  CALL QEXIT('REARSP')
      RETURN
      END

C
      SUBROUTINE READE3(LU,LEN,VEC,LAB1,LAB2,FREQ1,FREQ2,ISYM1,ISYM2,
     &                  RSD,ANTSYM)
C
C Modified version of REARSP with sole purpose of reading E3VEC
C
#include "implicit.h"
#include "thrzer.h"
#include "mxcent.h"
#include "priunit.h"
#include "infopt.h"
#include "inforb.h"
#include "gnrinf.h"
#include "optinf.h"
C
      CHARACTER*8 LABX1,LABX2,LAB1,LAB2,BLANK
      PARAMETER(THR = 1.0D-8, BLANK='        ', D0 = 0.0D0)
C
      CALL QENTER('READE3')

      REWIND (LU)
  100 READ(LU,END=200,ERR=200)
     &      LABX1,LABX2,FREQX1,FREQX2,ISYMX1,ISYMX2,ANTSYM,
     &      RSDX,LEN,EMCX,NBSX,NRBX
      IF (LABX1 .EQ. 'EOFLABEL') THEN
      WRITE(LUPRI,*) "Something wrong in reading E3VEC OR X2VEC &
     &    IN QRHYP"
      GOTO 200
      END IF
C
C
      IF ((LABX1 .EQ. LAB1) .AND. (LABX2 .EQ. LAB2)) THEN
      WRITE(LUPRI,*) "FOUND LABEL, writing it here:"
      WRITE(LUPRI,*) LAB1, LAB2
      ELSE 
      GOTO 100
      END IF
 200  CALL QEXIT('READE3')
      RETURN
      END

      SUBROUTINE READX2(LU,LEN,VEC,LAB1,LAB2,FREQ1,FREQ2,ISYM1,ISYM2,
     &                  RSD,ANTSYM,FIRST)
C
C Modified version of REARSP with sole purpose of reading X2VEC
C
#include "implicit.h"
#include "thrzer.h"
#include "mxcent.h"
#include "priunit.h"
#include "infopt.h"
#include "inforb.h"
#include "gnrinf.h"
#include "optinf.h"
C
      LOGICAL FIRST
      CHARACTER*8 LABX1,LABX2,LAB1,LAB2,BLANK
      PARAMETER(THR = 1.0D-8, BLANK='        ', D0 = 0.0D0)
C
      CALL QENTER('READX2')

      REWIND(LU)
  100 READ(LU,END=200,ERR=200)
     &      LABX1,LABX2,FREQX1,FREQX2,ISYMX1,ISYMX2,ANTSYM,
     &      RSDX,LEN,EMCX,NBSX,NRBX
      IF (LABX1 .EQ. 'EOFLABEL') THEN
      WRITE(LUPRI,*) "READ EOFLABEL in X2VEC, QRHYP"
      GOTO 200
      END IF
C
C
      IF ((LABX1 .EQ. LAB1) .AND. (LABX2 .EQ. LAB2) .AND. FIRST) THEN
      WRITE(LUPRI,*) "FOUND LABELS, writing it here:"
      WRITE(LUPRI,*) LAB1, LAB2
      ELSE 
      GOTO 100
      END IF
 200  CALL QEXIT('READX2')
      RETURN
      END
! -- end of rspcr5.F --
