C
C  /* Deck so_intrp */
      SUBROUTINE SO_INTRP(PRP1,LPRP1,PAOI,LPAOI,RTNLBL,ISYMTR)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Repack the AO property integrals to the order used
C              in the CC and SOPPA codes. The repacked integrals
C              are kept in a full square matrix.  
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER   (ONE = 1.0D0)
      DIMENSION   PAOI(LPAOI), PRP1(LPRP1)
      CHARACTER*8 RTNLBL(2)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_INTRP')
C
C------------------------------------
C     Repack square packed integrals.
C------------------------------------
C
      IF ( RTNLBL(2) .EQ. 'SQUARE  ' ) THEN
C
         NABP = 0
C
         DO 100 ISYMA = 1, NSYM
            DO 110 A = 1, NBAS(ISYMA)
               DO 120 ISYMB = 1, NSYM
                  DO 130 B = 1, NBAS(ISYMB)
                     NABP = NABP + 1
                     KOFF = IAODIS(ISYMA,ISYMB)
     &                    + NBAS(ISYMA)*(B-1) + A
                     PRP1(KOFF) = PAOI(NABP)
  130             CONTINUE
  120          CONTINUE
  110       CONTINUE
  100    CONTINUE
C
C----------------------------------------
C     Repack triangular packed integrals.
C----------------------------------------
C
      ELSE
C
         IF ( RTNLBL(2) .EQ. 'SYMMETRI' ) ANTSYM = ONE
         IF ( RTNLBL(2) .EQ. 'ANTISYMM' ) ANTSYM = -ONE
C
         NABP = 0
C
         DO 200 ISYMA = 1, NSYM
            DO 210 A = 1, NBAS(ISYMA)
               DO 220 ISYMB = 1, ISYMA
C
C--------------------------------------------------
C                 For total symmetric AO integrals.
C--------------------------------------------------
C
                  IF ( ISYMTR .EQ. 1 ) THEN
C
                     IF ( ISYMB .LT. ISYMA ) THEN
C
                        NABP = NABP + NBAS(ISYMB)
C
                     ELSE
C
                        DO 230 B = 1, A
                           NABP  = NABP + 1
                           KOFF1 = IAODIS(ISYMA,ISYMB) 
     &                           + NBAS(ISYMA)*(B-1) + A
                           KOFF2 = IAODIS(ISYMA,ISYMB) 
     &                           + NBAS(ISYMB)*(A-1) + B
                           PRP1(KOFF1) = PAOI(NABP)
                           PRP1(KOFF2) = PAOI(NABP) * ANTSYM
  230                   CONTINUE
C
                     ENDIF
C
C------------------------------------------------------
C                 For non-total symmetric AO integrals.
C------------------------------------------------------
C
                  ELSE
C
                     IF ( ISYMTR .EQ. MULD2H(ISYMA,ISYMB)) THEN
C
                        DO 240 B = 1, NBAS(ISYMB)
                           NABP  = NABP + 1
                           KOFF1 = IAODIS(ISYMA,ISYMB)
     &                           + NBAS(ISYMA)*(B-1) + A
                           KOFF2 = IAODIS(ISYMB,ISYMA)
     &                           + NBAS(ISYMB)*(A-1) + B
                           PRP1(KOFF1) = PAOI(NABP)
                           PRP1(KOFF2) = PAOI(NABP) * ANTSYM
  240                   CONTINUE
C
                     ELSE IF (ISYMB .LT. ISYMA) THEN
C
                        NABP = NABP + NBAS(ISYMB)
C
                     ELSE
C
                        DO 250 B = 1, A
                           NABP = NABP + 1
  250                   CONTINUE
C
                     ENDIF
C
                  ENDIF
C
  220          CONTINUE
  210       CONTINUE
  200    CONTINUE
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_INTRP')
C
      RETURN
      END
