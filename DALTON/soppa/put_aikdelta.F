C  /* Deck put_aikdelta */
      SUBROUTINE PUT_AIKDELTA(SCR2,OMEGA2,ISYMAI,ISYMBJ,NAI,IFSTAI)
C
C     Ove Christiansen 30-10-1995: Put in D contribution in omega vector
C                                  avoid troble on cray with optimization.
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION SCR2(*),OMEGA2(*)
C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccsdio.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF ( ISYMAI .EQ. ISYMBJ) THEN
         DO 190 NBJ = 1,NT1AM(ISYMBJ)
C
            KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + INDEX(NAI,NBJ)
C
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + HALF*SCR2(KOFF9)
C
  190    CONTINUE
C
      ENDIF
C
      IF ( ISYMAI .LT. ISYMBJ) THEN
         DO 200 NBJ = 1,NT1AM(ISYMBJ)
C
            KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + NT1AM(ISYMAI)*(NBJ - 1) + NAI
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + HALF*SCR2(KOFF9)
C
  200    CONTINUE
C
      ENDIF
C
      IF (ISYMBJ .LT. ISYMAI) THEN
         DO 210 NBJ = 1,NT1AM(ISYMBJ)
C
            KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + HALF*SCR2(KOFF9)
C
  210    CONTINUE
C
      ENDIF
C
      END
