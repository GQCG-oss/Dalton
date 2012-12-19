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
C
C FILE : gpcpp.F
C
C General Purpose routines for Complex Polarization Propagator
C
      SUBROUTINE READ_VEC(LU,N,X)
      INTEGER LU,N,I
      REAL*8 X(N)
      DO I=1,N
         READ(LU) X(I)
      END DO
      END
      SUBROUTINE WRITE_VEC(LU,N,X)
      INTEGER LU,N,I
      REAL*8 X(N)
      DO I=1,N
         WRITE(LU) X(I)
      END DO
      END
      SUBROUTINE WRITE_XVEC(LU,N,VEC,LABEL,FREQ,RESIDUAL_NORM)
C
#include "abslrs.h"
C
      INTEGER N,LU
      CHARACTER*8 LABEL,LAB1
      REAL*8 FREQ, RESIDUAL_NORM, VEC(N)
C
      REWIND (LU)
 100  READ(LU) LAB1
      IF (LAB1.NE.'EOFLABEL') GOTO 100
      BACKSPACE LU
      WRITE(LU) LABEL,FREQ,ABS_DAMP,ABS_GRADSYM,RESIDUAL_NORM,N
      CALL WRITE_VEC(LU,N,VEC)
      WRITE(LU) 'EOFLABEL'
      END
      SUBROUTINE READ_XVEC(LU,N,VEC,LABEL,ISYM,FREQ,
     &                  RSD,FOUND,CONV)
C
C Check if a linear equation needs to be solved or if a solution vector
C already exists on file. FOUND indicates the result.
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "abslrs.h"
C
      INTEGER ISYM,ISYMX1
      LOGICAL FOUND,CONV
      CHARACTER*8 LABX1,BLANK,LABEL
      PARAMETER(THR = 1.0D-8, BLANK='        ', D0 = 0.0D0)
c      DIMENSION TMP(4*ABS_MKZVAR)
C
      FOUND  = .FALSE.
      CONV = .FALSE.
      REWIND (LU)
C
  100 READ(LU,END=200,ERR=200)
     &      LABX1,FREQX1,DAMP1,ISYMX1,
     &      RSDX,LEN
      IF (LABX1 .EQ. 'EOFLABEL') GOTO 200
         IF (((ABS(ABS(FREQ)-FREQX1)).LE.THR)
     &        .AND. (LABEL.EQ.LABX1) .AND. (ISYM.EQ.ISYMX1) .AND. 
     &        (DAMP1.EQ.ABS_DAMP)) THEN
c     &        .AND. (ABS_GRADSYM.EQ.ISYMX1)) THEN 
             IF (LEN .NE. N) THEN
                WRITE (LUABSPRI,*)
     &          ' Orbital data found on response'//
     &          ' vector file does not match current orbital data',
     &          ' LEN      : ',LEN,N
                GOTO 200
             ENDIF
            FOUND = .TRUE.
            IF ((RSDX-RSD) .LT. D0) CONV = .TRUE.
            IF (LEN .GT. 0) CALL READ_VEC(LU,LEN,VEC)
            GOTO 199
         ELSE
c            IF (LEN .GT. 0) READ (LU)
            IF (LEN .GT. 0) THEN
              DO I=1,LEN
               READ(LU) SCR
              END DO
            ENDIF
            GO TO 100
         END IF
C
 199  CONTINUE
 200  RETURN
      END
      SUBROUTINE ABS_CHKONFILE(LU,FOUND,LABEL,ISYM,DAMP,NFREQ_ABS,
     &                        FREQ_ABS,THD,FLAGS)
C
      LOGICAL     FOUND,FLAGS(NFREQ_ABS)
      CHARACTER*8 LABEL,LAB1
      INTEGER     LU,ISYM,NFREQ_ABS,ISYM1,LEN
      REAL*8      FREQ_ABS(NFREQ_ABS),THD,FREQ1,DAMP1,RSD,DAMP
C
      FOUND = .TRUE.
      DO I=1,NFREQ_ABS
         FLAGS(I)=.FALSE.
      END DO
C
      REWIND(LU)
 100  READ(LU,END=900,ERR=900) LAB1,FREQ1,DAMP1,ISYM1,
     &     RSD,LEN
C
      IF (LAB1.NE.LABEL .OR.
     &   ISYM1.NE.ISYM  .OR.
     &   RSD  .GT.THD   .OR. DAMP1.NE.DAMP) THEN
         IF (LEN .GT. 0) THEN
            DO I=1,LEN
              READ(LU) SCR
            END DO
         ENDIF
         GOTO 100
      ELSE
C
      DO I=1,NFREQ_ABS
         IF (FREQ1.EQ.FREQ_ABS(I)) THEN
            FLAGS(I)=.TRUE.
         END IF
      END DO
C
      DO I=1,LEN
        READ(LU,END=900,ERR=900) SCR
      END DO
C      
      GOTO 100
      ENDIF
C
 900  CONTINUE
      DO I=1,NFREQ_ABS
         FOUND = FOUND .AND. FLAGS(I)
      END DO
      RETURN
      END

