C  /* Deck put_kiadelta */
      SUBROUTINE PUT_KIADELTA(SCR2,OMEGA2,ISYMAI,ISYMAJ,ISYMBI,ISYMBJ,
     *                   ISYMB,ISYMI,ISYMJ,NAI,NAJ,I,J,IFSTAI)
C
C     Ove Christiansen 30-10-1995: Put in C contribution in omega vector
C                                  avoid troble on cray with optimization.
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION SCR2(*),OMEGA2(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccsdio.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF ( ISYMAI .EQ. ISYMBJ ) THEN
C
         DO 100 B = 1,NVIR(ISYMB)
C
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + INDEX(NAI,NBJ)
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                      * SCR2(KOFF9)
C     *                 - HALF * SCR2(KOFF9)
C
            WRITE(LUPRI,*) 'ISYMAI=ISYMBJ'
            WRITE(LUPRI,*) 'KOFF9',KOFF9,'NAIBJ', NAIBJ
            WRITE(LUPRI,*) SCR2(KOFF9),OMEGA2(NAIBJ)
C
  100    CONTINUE
C
      ENDIF
C
      IF ( ISYMAI .LT. ISYMBJ ) THEN
C
         DO 200 B = 1,NVIR(ISYMB)
C
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + NT1AM(ISYMAI)*(NBJ - 1) + NAI
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                 - HALF * SCR2(KOFF9)
C
            WRITE(LUPRI,*) 'ISYMAI<ISYMBJ'
            WRITE(LUPRI,*) 'KOFF9',KOFF9,'NAIBJ', NAIBJ
            WRITE(LUPRI,*) SCR2(KOFF9)
  200    CONTINUE
C
      ENDIF
C
      IF ( ISYMBJ .LT. ISYMAI ) THEN
C
         DO 300 B = 1,NVIR(ISYMB)
C
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                 - HALF * SCR2(KOFF9)
C
            WRITE(LUPRI,*) 'ISYMAI>ISYMBJ'
            WRITE(LUPRI,*) 'KOFF9',KOFF9,'NAIBJ', NAIBJ
            WRITE(LUPRI,*) SCR2(KOFF9)
  300    CONTINUE
C
      ENDIF
C
      IF (ISYMAJ .EQ. ISYMBI) THEN
C
         DO 400 B = 1,NVIR(ISYMB)
C
            NBI = IT1AM(ISYMB,ISYMI)
     *          + NVIR(ISYMB)*(I-1) + B
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *            + INDEX(NAJ,NBI)
            OMEGA2(NAJBI) = OMEGA2(NAJBI) - SCR2(KOFF9)
C
            WRITE(LUPRI,*) 'ISYMAJ=ISYMBI'
            WRITE(LUPRI,*) 'KOFF9',KOFF9,'NAIBJ', NAIBJ
            WRITE(LUPRI,*) SCR2(KOFF9)
  400    CONTINUE
C
      ENDIF
C
      IF (ISYMAJ .LT. ISYMBI) THEN
C
         DO 500 B = 1,NVIR(ISYMB)
C
            NBI = IT1AM(ISYMB,ISYMI)
     *          + NVIR(ISYMB)*(I-1) + B
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *            + NT1AM(ISYMAJ)*(NBI - 1)
     *            + NAJ
            OMEGA2(NAJBI) = OMEGA2(NAJBI) - SCR2(KOFF9)
C
            WRITE(LUPRI,*) 'ISYMAJ<ISYMBI'
            WRITE(LUPRI,*) 'KOFF9',KOFF9,'NAIBJ', NAIBJ
            WRITE(LUPRI,*) SCR2(KOFF9)
  500    CONTINUE
C
      ENDIF
C
      IF (ISYMBI .LT. ISYMAJ) THEN
C
         DO 600 B = 1,NVIR(ISYMB)
C
            NBI = IT1AM(ISYMB,ISYMI)
     *          + NVIR(ISYMB)*(I-1) + B
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *            + NT1AM(ISYMBI)*(NAJ - 1)
     *            + NBI
            OMEGA2(NAJBI) = OMEGA2(NAJBI) - SCR2(KOFF9)
C
            WRITE(LUPRI,*) 'ISYMAJ>ISYMBI'
            WRITE(LUPRI,*) 'KOFF9',KOFF9,'NAIBJ', NAIBJ
            WRITE(LUPRI,*) SCR2(KOFF9)
  600    CONTINUE
C
      ENDIF
C
      END
