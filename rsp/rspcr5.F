C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html 
C
C
C$Id: rspcr5.F,v 1.1.1.1 2001-02-08 13:33:13 hjj Exp $
C
#include <single.h>
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
      FOUND = .FALSE.
      CONV = .FALSE.
      USOLVD = .FALSE.
      REWIND (LU)
  100 READ(LU,END=200,ERR=200)
     &      LABX1,LABX2,FREQX1,FREQX2,ISYMX1,ISYMX2,ANTSYM,
     &      RSDX,LEN,EMCX,NBSX,NRBX
      IF (LABX1 .EQ. 'EOFLABEL') THEN
         RETURN
      END IF
      IF ((NRBX .NE. NORBT) .OR. (NBSX .NE. NBAST)) THEN
            WRITE (LUPRI,'(/A,2(/A,I5,2X,I5))')
     &           ' Orbital data found on response'//
     &           ' vector file does not match current orbital data',
     &           ' Number of basis functions/orbitals on file      : ',
     &           NBSX,NRBX,
     &           ' Number of basis functions/orbitals in this calc.: ',
     &           NBAST,NORBT
            WRITE (LUPRI,'(/A)') ' Calculation will be stopped'
            CALL QUIT('Inconsistent energy found on response '//
     &                'vector file')
      ELSE
         IF (ABS(EMCX - EMCSCF) .GT. THRZER) THEN
            IF ((OPTNEW .AND. ITRNMR .GT. 0) .OR. 
     &          (OPTWLK .AND. ITERNR .GT. 0)) THEN
               USOLVD = .TRUE.
               WRITE (LUPRI,'(/A)') ' Converged vectors from previous'//
     &              ' geometry used as starting vectors'
            ELSE
               WRITE (LUPRI,'(/A,2(/A,F20.10))')
     &              ' Molecular energy found on response'//
     &              ' vector file does not match current energy',
     &              ' Energy on file            : ',EMCX,
     &              ' Energy in this calculation: ',EMCSCF
               WRITE (LUPRI,'(/A)') ' Calculation will be stopped'
               CALL QUIT('Inconsistent energy found on response '//
     &                   'vector file')
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
      IF (USOLVD) CONV = .FALSE.
      RETURN
 200  END
      SUBROUTINE CR0IND(LU)
C
C Purpose: Open RSPVEC
C
#include <implicit.h>
#include <dummy.h>
#include <infcr.h>
C
      CHARACTER*8 LAB1,LAB2,BLANK
      LOGICAL     EX
C
      CALL GPINQ('RSPVEC','EXIST',EX) 
      IF (.NOT.EX) THEN 
         CALL GPOPEN(LU,'RSPVEC','NEW',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         WRITE (LU) 'EOFLABEL'
         RETURN
      ELSE
         CALL GPOPEN(LU,'RSPVEC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
      END IF
C
      RETURN
      END
