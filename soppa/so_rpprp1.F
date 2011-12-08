C
C  /* Deck so_rpprp1 */
      SUBROUTINE SO_RPPRP1(PRP1,LPRP1,PR1IJ,LPR1IJ,PR1AB,LPR1AB,ISYMTR)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C
C     PURPOSE: Repack the one-particle property mo-integrals (PRP1) in 
C              an occupied-occupied matirx (PR1IJ) and a virtual-virtual
C              matrix (PR1AB).
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccsdsym.h>
#include <ccorb.h>
C
      DIMENSION   PRP1(LPRP1),   PR1IJ(LPR1IJ), PR1AB(LPR1AB)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RPPRP1')
C
C--------------------------------------------
C     Make occupied-occupied property matrix.
C--------------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMI = MULD2H(ISYMJ,ISYMTR)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            KOFF1 = IFCRHF(ISYMI,ISYMJ) + NORB(ISYMI)*(J - 1) + 1
            KOFF2 = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + 1
C
            CALL DCOPY(NRHF(ISYMI),PRP1(KOFF1),1,PR1IJ(KOFF2),1)
C
  110    CONTINUE
C
  100 CONTINUE
C
C------------------------------------------
C     Make virtual-virtual property matrix.
C------------------------------------------
C
      DO 200 ISYMB = 1,NSYM
C
         ISYMA = MULD2H(ISYMB,ISYMTR)
C
         DO 210 B = 1,NVIR(ISYMB)
C
            KOFF1 = IFCVIR(ISYMA,ISYMB) + NORB(ISYMA)*(B - 1)
     &            + NRHF(ISYMA) + 1
            KOFF2 = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1) + 1
C
            CALL DCOPY(NVIR(ISYMA),PRP1(KOFF1),1,PR1AB(KOFF2),1)
C
  210    CONTINUE
C
  200 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RPPRP1')
C
      RETURN
      END
