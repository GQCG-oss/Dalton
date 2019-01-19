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
C  /* Deck dicsym */
#if defined (VAR_DICSYM)
      SUBROUTINE DICSYM(ICSYM,CBLOCK,CTBLCK,NIASTS,NIBSTS)
C
C 28-Jun-1987 Hans Joergen Aa. Jensen
C
C Check if   CBLOCK(NIASTS,NIBSTS) is symmetric or antisymmetric
C       with CTBLCK(NIBSTS,NIASTS).
C
C On output ICSYM  = 10 if zero blocks
C                  =  1 if symmetric
C                  = -1 if antisymmetric
C                  =  0 if neither symmetric nor antisymmetric
C
#include "implicit.h"
      DIMENSION CBLOCK(NIASTS,NIBSTS), CTBLCK(NIBSTS,NIASTS)
      PARAMETER (THRSML = 1.0D-14)
      LOGICAL   SYM, ASYM
C
      SYM  = .TRUE.
      ASYM = .TRUE.
      DO 200 J = 1,NIBSTS
         DO 100 I = 1,NIASTS
            IF (ABS(CBLOCK(I,J) - CTBLCK(J,I)) .GT. THRSML) THEN
               SYM  = .FALSE.
            ELSE IF (ABS(CBLOCK(I,J) + CTBLCK(J,I)) .GT. THRSML) THEN
               ASYM = .FALSE.
            END IF
  100    CONTINUE
         IF (.NOT. SYM .AND. .NOT.ASYM) GO TO 300
  200 CONTINUE
  300 CONTINUE
      IF (SYM .AND. ASYM) THEN
C        ... zero block(s).
         ICSYM = 10
      ELSE IF (SYM) THEN
         ICSYM = 1
      ELSE IF (ASYM) THEN
         ICSYM = -1
      ELSE
         ICSYM = 0
      END IF
      RETURN
C     ... end of di_csym.
      END
#endif
C  /* Deck disngl */
#if defined (VAR_DISNGL)
      SUBROUTINE DISNGL(EQUAL,CBLOCK,CTBLCK,NIASTS,NIBSTS)
C
C  1-Jul-1987 Hans Joergen Aa. Jensen
C
C
C This is a singlet wave function, Eaa*Eaa contribution in HC,
C get Ebb*Ebb contribution by adding
C CBLOCK(NIASTS,NIBSTS) transposed to CTBLCK(NIBSTS,NIASTS)
C and vice versa.
C
#include "implicit.h"
      LOGICAL   EQUAL
      DIMENSION CBLOCK(NIASTS,NIBSTS), CTBLCK(NIBSTS,NIASTS)
C
      IF (EQUAL) THEN
         DO 1200 J = 1,NIBSTS
            DO 1100 I = 1,J
               SUM = CBLOCK(I,J) + CBLOCK(J,I)
               CBLOCK(I,J) = SUM
               CBLOCK(J,I) = SUM
 1100       CONTINUE
 1200    CONTINUE
      ELSE
         DO 2200 J = 1,NIBSTS
            DO 2100 I = 1,NIASTS
               SUM = CBLOCK(I,J) + CTBLCK(J,I)
               CBLOCK(I,J) = SUM
               CTBLCK(J,I) = SUM
 2100       CONTINUE
 2200    CONTINUE
      END IF
      RETURN
C     ... end of di_sngl.
      END
#endif
