C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
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
C$Id: snglet.F,v 1.2 2001-05-31 06:38:38 vebjornb Exp $
C
#include <single.h>
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
#include <implicit.h>
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
#include <implicit.h>
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
