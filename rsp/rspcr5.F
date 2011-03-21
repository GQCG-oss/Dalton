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
      SUBROUTINE WRTRSP(LU,LEN,VEC,LAB1,LAB2,FREQ1,FREQ2,
     &                  ISYM1,ISYM2,RSD,ANTSYM)
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <infrsp.h>
#include <infopt.h>
#include <inforb.h>
C
      CHARACTER*8 LAB1,LAB2,BLANK,LABEL
      PARAMETER   (BLANK='        ', D0 = 0.0D0)
      DIMENSION VEC(*)
C
      REWIND (LU)
 100  READ(LU) LABEL
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
      IF (LEN .GT. 0) CALL WRITI(LU,IRAT*LEN,VEC)
      WRITE(LU) 'EOFLABEL'
C
      CALL FLSHFO(LU)
C
      RETURN
      END
      SUBROUTINE REARSP(LU,LEN,VEC,LAB1,LAB2,FREQ1,FREQ2,ISYM1,ISYM2,
     &                  RSD,FOUND,CONV,ANTSYM)
C
C Check if a linear equation needs to be solved or if a solution vector
C already exists on file. FOUND indicates the result.
C
#include <implicit.h>
#include <iratdef.h>
#include <thrzer.h>
#include <mxcent.h>
#include <priunit.h>
#include <infopt.h>
#include <inforb.h>
#include <gnrinf.h>
#include <optinf.h>
C
      LOGICAL FOUND, CONV, USOLVD
      CHARACTER*8 LABX1,LABX2,LAB1,LAB2,BLANK
      PARAMETER(THR = 1.0D-8, BLANK='        ', D0 = 0.0D0)
C
      FOUND  = .FALSE.
      CONV   = .FALSE.
      USOLVD = .FALSE.
      REWIND (LU)
C
  100 READ(LU,END=200,ERR=200)
     &      LABX1,LABX2,FREQX1,FREQX2,ISYMX1,ISYMX2,ANTSYM,
     &      RSDX,LEN,EMCX,NBSX,NRBX
      IF (LABX1 .EQ. 'EOFLABEL') GOTO 200
C
      IF ((NRBX .NE. NORBT) .OR. (NBSX .NE. NBAST)) THEN
         WRITE (LUPRI,'(/A,2(/A,I5,2X,I5))')
     &        ' Orbital data found on response'//
     &        ' vector file does not match current orbital data',
     &        ' Number of basis functions/orbitals on file      : ',
     &        NBSX,NRBX,
     &        ' Number of basis functions/orbitals in this calc.: ',
     &        NBAST,NORBT
         WRITE (LUPRI,'(/A)') ' RSPVEC is reset to empty file'
         REWIND (LU)
         WRITE (LU) 'EOFLABEL'
         GOTO 200
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
               WRITE (LUPRI,'(/A)') ' RSPVEC is reset to empty file'
               REWIND (LU)
               WRITE (LU) 'EOFLABEL'
               GOTO 200
            END IF
         END IF
      END IF
      IF (LAB2 .EQ. '        ') THEN
         IF (((ABS(ABS(FREQ1)-FREQX1)).LE.THR)
     &        .AND. (LAB1.EQ.LABX1) .AND. (LAB2.EQ.LABX2) .AND.
     &        (ISYM1.EQ.ISYMX1)) THEN 
            FOUND = .TRUE.
            IF ((RSDX-RSD) .LT. D0) CONV = .TRUE.
            IF (LEN .GT. 0) CALL READI(LU,IRAT*LEN,VEC)
            GOTO 199
         ELSE
            IF (LEN .GT. 0) READ (LU)
            GO TO 100
         END IF
      ELSE
         IF (((ABS(FREQ1-FREQX1).LE.THR)
     &        .AND. (ABS(FREQ2-FREQX2).LE.THR)
     &        .AND. LAB1.EQ.LABX1 .AND. LAB2.EQ.LABX2 .AND.
     &        ISYM1.EQ.ISYMX1 .AND. ISYM2.EQ.ISYMX2) .OR.
     &        ((ABS(FREQ1-FREQX2).LE.THR)
     &        .AND. (ABS(FREQ2-FREQX1).LE.THR)
     &        .AND. LAB1.EQ.LABX2 .AND. LAB2.EQ.LABX1 .AND.
     &        ISYM1.EQ.ISYMX2 .AND. ISYM2.EQ.ISYMX1)) THEN 
            FOUND = .TRUE.
            IF ((RSDX-RSD) .LT. D0) CONV = .TRUE.
            IF (LEN .GT. 0) CALL READI(LU,IRAT*LEN,VEC)
            GOTO 199
         ELSE
            IF (LEN .GT. 0) READ (LU)
            GO TO 100
         END IF
      END IF
C
 199  CONTINUE
      IF (USOLVD) THEN
         CONV = .FALSE.
         WRITE (LUPRI,'(/A,/2A)') ' Converged vectors from'//
     &      ' previous calculation are used as starting '//
     &      'vectors for property',LABX1, LABX2
      END IF
 200  RETURN
      END
