C
C  /* Deck rp_newtrial */
      SUBROUTINE RP_NEWTRIAL(NNEWTR,NOLDTR,TR1E,LTR1E,TR1D,LTR1D,EIVAL1,
     &                       EDIA1,LEDIA1,RESI1E,LRESI1E,RESI1D,LRESI1D)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate a new trial vector from the residual vector
C              and the diagonal parts of the Hessian and overlap
C              matrices.
C              The orthogonalization of the new trial vector is
C              not performed here but in SO_ORTH_TRN.
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (DTHRES = 1.0D-04)
C
      DIMENSION TR1E(LTR1E),     TR1D(LTR1D)
      DIMENSION EDIA1(LEDIA1)
      DIMENSION RESI1E(LRESI1E), RESI1D(LRESI1D)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_NEWTRIAL')
C
C----------------------------------------------------------------
C     Calculate 1p1h excitation part of the new raw trial vector.
C----------------------------------------------------------------
C
      DO 100 I = 1,LTR1E
C
         DIFF = EDIA1(I) - EIVAL1
C
C-----------------------------------------------------------------------
C        Ensure that the absolute denominator is greater than threshold.
C-----------------------------------------------------------------------
C
         IF ( DABS(DIFF) .GE. DTHRES ) THEN
C
            TR1E(I) = RESI1E(I) / DIFF
C
         ELSE
C
            TR1E(I) = RESI1E(I) / DSIGN(DTHRES,DIFF)
C
         END IF
C
  100 CONTINUE
C
C-------------------------------------------------------------------
C     Calculate 1p1h de-excitation part of the new raw trial vector.
C-------------------------------------------------------------------
C
      DO 101 I = 1,LTR1D
C
         SUM  = EDIA1(I) + EIVAL1
C
C-----------------------------------------------------------------------
C        Ensure that the absolute denominator is greater than threshold.
C-----------------------------------------------------------------------
C
         IF ( DABS(SUM) .GE. DTHRES ) THEN
C
            TR1D(I) = RESI1D(I) / SUM
C
         ELSE
C
            TR1D(I) = RESI1D(I) / DSIGN(DTHRES,SUM)
C
         END IF
C
  101 CONTINUE
C
C---------------------------------------------------
C     Write raw new trial vector to file and output.
C---------------------------------------------------
C
      CALL SO_WRITE(TR1E,LTR1E,LUTR1E,FNTR1E,NOLDTR+NNEWTR)
      CALL SO_WRITE(TR1D,LTR1D,LUTR1D,FNTR1D,NOLDTR+NNEWTR)
C
      IF ( IPRSOP .GE. 7 ) THEN
C
         CALL AROUND('New raw trialvector in RP_NEWTRIAL')
C
         WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &        (I,TR1E(I),TR1D(I),I=1,LTR1E)
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_NEWTRIAL')
C
      RETURN
      END
