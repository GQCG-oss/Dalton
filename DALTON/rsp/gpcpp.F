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
C FILE : gpcpp.F
C
C General Purpose routines for Complex Polarization Propagator
C
      SUBROUTINE READ_VEC(LU,N,X)
      INTEGER LU,N
      REAL*8 X(N)
      READ(LU) X
      END
      SUBROUTINE WRITE_VEC(LU,N,X)
      INTEGER LU,N
      REAL*8 X(N)
      WRITE(LU) X
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
C
      SUBROUTINE WRITE_XVEC2(LU,N,VEC,LABEL1,LABEL2,FREQ1,FREQ2,
     &                       RESIDUAL_NORM)
C
#include "abslrs.h"
C
      INTEGER N,LU
      CHARACTER*8 LABEL1,LABEL2,LAB1
      REAL*8 FREQ1, FREQ2, RESIDUAL_NORM, VEC(N)
C
      REWIND (LU)
 100      READ(LU) LAB1
      IF (LAB1.NE.'EOFLABEL') GOTO 100
      BACKSPACE LU
      WRITE(LU) LABEL1, LABEL2,FREQ1, FREQ2, ABS_DAMP,ABS_GRADSYM,
     &          RESIDUAL_NORM,N
      CALL WRITE_VEC(LU,N,VEC)
      WRITE(LU) 'EOFLABEL'
C
C     End of WRITE_XVEC2
C
      END
C
C
C
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
      DIMENSION VEC(N)
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
         IF (((ABS(ABS(FREQ)-abs(FREQX1))).LE.THR)
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
            IF (LEN .GT. 0) READ (LU)
            GO TO 100
         END IF
C
 199  CONTINUE
 200  RETURN
      END
C
      SUBROUTINE READ_XVEC2(LU,N,VEC,LABEL1,LABEL2,ISYM,FREQ1,FREQ2,
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
      CHARACTER*8 LABX,LABY,BLANK,LABEL1,LABEL2
      PARAMETER(THR = 1.0D-8, BLANK='        ', D0 = 0.0D0)
      PARAMETER(THR0=1.0D-10)
      DIMENSION VEC(N)
c      DIMENSION TMP(4*ABS_MKZVAR)
C
      FOUND  = .FALSE.
      CONV = .FALSE.
      REWIND (LU)
C
 100     READ(LU,END=200,ERR=200)
     &        LABX,LABY,FREQX,FREQY,DAMP1,ISYMX1,
     &        RSDX,LEN
      IF (LABX .EQ. 'EOFLABEL') GOTO 200
c      write(luabspri,*)'freq',ABS(FREQ1),FREQX,
c     &                        DAMP1,
c     &                        ABS_DAMP,ABS_GRADSYM,ISYMX1,ISYM

      IF ((LABEL2 .EQ. BLANK) .AND. ((ABS(ABS(FREQX)-FREQ1)).LE.THR)
     &        .AND. (LABEL1.EQ.LABX) .AND. (ISYM.EQ.ISYMX1) .AND.
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
            IF ((RSDX-RSD) .LE. THR0) CONV = .TRUE.
            IF (LEN .GT. 0) CALL READ_VEC(LU,LEN,VEC)
            GOTO 199
C Removed the ABS around FREQX/Y - Tobias
        ELSEIF (((ABS(FREQX-FREQ1)).LE.THR) .AND.
     &        ((ABS(FREQY-FREQ2)).LE.THR) .AND. (LABEL1.EQ.LABX)
     &        .AND. (LABEL2.EQ.LABY) .AND. (ISYM.EQ.ISYMX1) .AND.
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
            IF ((RSDX-RSD) .LE. THR0) CONV = .TRUE.
            IF (LEN .GT. 0) CALL READ_VEC(LU,LEN,VEC)
            GOTO 199
         ELSE
            IF (LEN .GT. 0) READ (LU)
            GO TO 100
         END IF
C
 199       CONTINUE
 200         RETURN
C
C     End of READ_XVEC2
C
      END
C
C
C
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
         IF (LEN .GT. 0) READ(LU)
         GOTO 100
      ELSE
C
      DO I=1,NFREQ_ABS
         IF (FREQ1.EQ.FREQ_ABS(I)) THEN
            FLAGS(I)=.TRUE.
         END IF
      END DO
C
      IF (LEN.GT.0) READ(LU,END=900,ERR=900)
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

      SUBROUTINE ABS_REDSPACE_REBUILD(KZVAR,GD,REDE,REDS,REDGD,
     &          KNVEC,WRK,LWRK)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "abslrs.h"
C
C PURPOSE:
C 
C      rebuild reduce space from vectors on disk
C
      DIMENSION GD(KZVAR,4),REDE(ABS_MAXRM,ABS_MAXRM,2)
      DIMENSION REDS(ABS_MAXRM,ABS_MAXRM),REDGD(ABS_MAXRM,2)
      DIMENSION WRK(LWRK)
      INTEGER KNVEC(2)

      KBVEC = 1
      KEVEC = KBVEC + KZVAR
      KSVEC = KEVEC + KZVAR
      KFREE = KSVEC + KZVAR
      LFREE = LWRK - KFREE
c
      PR1=0.0d0
      REWIND(LUSB)
      DO I=1,ABS_KLRED(1)
         CALL READ_VEC(LUSB,KZVAR,WRK(KBVEC))
          REDGD(I,1)= ABS(DDOT(KZVAR,WRK(KBVEC),1,GD(1,1),1))
          PR1=MAX(PR1,REDGD(I,1))
      ENDDO
      REWIND(LUAB)
      DO I=1,ABS_KLRED(2)
         CALL READ_VEC(LUAB,KZVAR,WRK(KBVEC))
         REDGD(I,2)= ABS(DDOT(KZVAR,WRK(KBVEC),1,GD(1,2),1))
         PR1=MAX(PR1,REDGD(I,2))
      ENDDO
      IF (ABS(PR1).LE.(1.0d-8)) THEN
        REWIND(LUSB)
        REWIND(LUAB)
        REWIND(LUSS)
        REWIND(LUAS)
        CALL DZERO(REDGD(1,1),ABS_KLRED(1))
        CALL DZERO(REDGD(1,2),ABS_KLRED(2))
        ABS_KLRED(1)=0
        ABS_KLRED(2)=0
        KNVEC(1)=2*ABS_NFREQ_INTERVAL
        KNVEC(2)=2*ABS_NFREQ_INTERVAL
      ELSE
c
C        rebuild reduced spaces
c
        REWIND(LUSB)
        DO I=1,ABS_KLRED(1)
          CALL READ_VEC(LUSB,KZVAR,WRK(KBVEC))
          REWIND(LUSS)
          DO J=1,I
            CALL READ_VEC(LUSS,KZVAR,WRK(KEVEC))
            REDE(I,J,1)= 2.0d0*DDOT(KZVAR,WRK(KBVEC),1,WRK(KEVEC),1)
            IF (I .NE. J) REDE(J,I,1)=REDE(I,J,1)
          ENDDO
        REWIND(LUAB)
          DO J=1,ABS_KLRED(2)
            CALL READ_VEC(LUAB,KZVAR,WRK(KEVEC))
            CALL GETSVEC(KZVAR,1,WRK(KEVEC),WRK(KSVEC))
            REDS(I,J)=2.0d0*DDOT(KZVAR,WRK(KBVEC),1,WRK(KSVEC),1)
          ENDDO
        ENDDO
        REWIND(LUAB)
        DO I=1,ABS_KLRED(2)
          CALL READ_VEC(LUAB,KZVAR,WRK(KBVEC))
          REWIND(LUAS)
          DO J=1,I
            CALL READ_VEC(LUAS,KZVAR,WRK(KEVEC))
            REDE(I,J,2)= 2.0d0*DDOT(KZVAR,WRK(KBVEC),1,WRK(KEVEC),1)
            IF (I .NE. J) REDE(J,I,2)=REDE(I,J,2)
          ENDDO
        ENDDO
        KNVEC(1)=0
        KNVEC(2)=0
      ENDIF
c
      RETURN
      END
      SUBROUTINE ABS_READIN_RED(KZVAR,REDE,REDS,GD,
     &           REDGD,KNVEC,WRK,LWRK)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "abslrs.h"
C
C PURPOSE:
C 
C      rebuild reduce space from vectors on disk
C
      DIMENSION GD(KZVAR,4),REDE(ABS_MAXRM,ABS_MAXRM,2)
      DIMENSION REDS(ABS_MAXRM,ABS_MAXRM),REDGD(ABS_MAXRM,2)
      DIMENSION WRK(LWRK)
      INTEGER KNVEC(2)

      KBVEC = 1
      KEVEC = KBVEC + KZVAR
      KSVEC = KEVEC + KZVAR
      KFREE = KSVEC + KZVAR
      LFREE = LWRK - KFREE
c
      PR1=0.0d0
      REWIND(LUSB)
      DO I=1,ABS_KLRED(1)
         CALL READ_VEC(LUSB,KZVAR,WRK(KBVEC))
          REDGD(I,1)= ABS(DDOT(KZVAR,WRK(KBVEC),1,GD(1,1),1))
          PR1=MAX(PR1,REDGD(I,1))
      ENDDO
      REWIND(LUAB)
      DO I=1,ABS_KLRED(2)
         CALL READ_VEC(LUAB,KZVAR,WRK(KBVEC))
         REDGD(I,2)= ABS(DDOT(KZVAR,WRK(KBVEC),1,GD(1,2),1))
         PR1=MAX(PR1,REDGD(I,2))
      ENDDO
      IF (ABS(PR1).LE.(1.0d-8)) THEN
        REWIND(LUSB)
        REWIND(LUAB)
        REWIND(LUSS)
        REWIND(LUAS)
        CALL DZERO(REDGD(1,1),ABS_KLRED(1))
        CALL DZERO(REDGD(1,2),ABS_KLRED(2))
        ABS_KLRED(1)=0
        ABS_KLRED(2)=0
        KNVEC(1)=2*ABS_NFREQ_INTERVAL
        KNVEC(2)=2*ABS_NFREQ_INTERVAL
        REWIND(LUE1RED)
        REWIND(LUE2RED)
        REWIND(LUSRED)
        CALL DZERO(REDE(1,1,1),2*ABS_MAXRM*ABS_MAXRM)
        CALL DZERO(REDS(1,1),ABS_MAXRM*ABS_MAXRM)
      ELSE
c
C        rebuild reduced spaces
c
        REWIND(LUE1RED)
        REWIND(LUE2RED)
        REWIND(LUSRED)
        DO I=1,ABS_KLRED(1)
          CALL READ_VEC(LUE1RED,ABS_KLRED(1),REDE(1,I,1))
        ENDDO
        DO I=1,ABS_KLRED(2)
          CALL READ_VEC(LUE2RED,ABS_KLRED(2),REDE(1,I,2))
          CALL READ_VEC(LUSRED,ABS_KLRED(1),REDS(1,I))
        ENDDO
        KNVEC(1)=0
        KNVEC(2)=0
      ENDIF
c      write(luabspri,*)'E1 matrix',ABS_KLRED(1)
c      CALL OUTPUT(REDE(1,1,1),1,ABS_KLRED(1),1,ABS_KLRED(1),
c     &           ABS_MAXRM,ABS_MAXRM,1,LUABSPRI)
c      write(luabspri,*)'E2 matrix',ABS_KLRED(2)
c      CALL OUTPUT(REDE(1,1,2),1,ABS_KLRED(2),1,ABS_KLRED(2),
c     &           ABS_MAXRM,ABS_MAXRM,1,LUABSPRI)
c      write(luabspri,*)'S matrix'
c      CALL OUTPUT(REDS(1,1),1,ABS_KLRED(1),1,ABS_KLRED(2),
c     &           ABS_MAXRM,ABS_MAXRM,1,LUABSPRI)
c
      RETURN
      END
! end of gpcpp.F
