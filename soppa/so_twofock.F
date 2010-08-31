#include <single.h>
C
C  /* Deck so_twofock */
      SUBROUTINE SO_TWOFOCK(RES1E,LRES1E,RES1D,LRES1D,FOCK,LFOCK,ISYRES)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, December 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Add the two particle part of the Fock matrix to the
C              results vectors (I.e. the two particle part of the
C              RPA problem). See eq. (19) in Chem. Phys. 172, 13 (1993).
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1E(LRES1E), RES1D(LRES1D), FOCK(LFOCK)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TWOFOCK')
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYRES)
C
         KOFFA = IFCRHF(ISYMA,ISYMI) + NRHF(ISYMA) + 1
         KOFFB = IT1AM(ISYMA,ISYMI) + 1
         KOFFC = IFCVIR(ISYMI,ISYMA)
C
         DO 200 I = 1,NRHF(ISYMI)
C
            KOFF1 = KOFFA + ( NORB(ISYMA) * (I-1) )
            KOFF2 = KOFFB + ( NVIR(ISYMA) * (I-1) )
            KOFF3 = KOFFC + I

            CALL DAXPY(NVIR(ISYMA),ONE,FOCK(KOFF1),1,RES1E(KOFF2),1)
            CALL DAXPY(NVIR(ISYMA),-ONE,FOCK(KOFF3),NORB(ISYMI),
     &                 RES1D(KOFF2),1)
C
  200    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TWOFOCK')
C
      RETURN
      END
