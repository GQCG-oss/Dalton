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
      SUBROUTINE READ_XVEC(LU,N,VEC,LABEL,FREQ,
     &                  RSD,FOUND,CONV)
C
C Check if a linear equation needs to be solved or if a solution vector
C already exists on file. FOUND indicates the result.
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "abslrs.h"
C
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
c      write(luabspri,*)'freq',ABS(FREQ),FREQX1,LABEL,LABX1,DAMP1,
c     &                        ABS_DAMP,ABS_GRADSYM,ISYMX1
         IF (((ABS(ABS(FREQ)-FREQX1)).LE.THR)
     &        .AND. (LABEL.EQ.LABX1) .AND. (DAMP1.EQ.ABS_DAMP)) THEN
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
