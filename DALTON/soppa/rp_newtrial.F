C
C  /* Deck rp_newtrial */
      SUBROUTINE RP_NEWTRIAL(NNEWTR,NOLDTR,
     &                       TR1E,LTR1E,TR1D,LTR1D,
     &                       EIVAL1,EDIA1,LEDIA1,
     &                       RESI1E,LRESI1E,RESI1D,LRESI1D)
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
         use so_info, only:  sop_dp 
C
      implicit none
C#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
C
      REAL(SOP_DP), PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, 
     &                           ONE = 1.0D0, TWO = 2.0D0
C
C     Array lengths       
      INTEGER, INTENT(IN) ::  LTR1E, LTR1D, 
     &                        LEDIA1,
     &                        LRESI1E, LRESI1D

      INTEGER, INTENT(IN) ::  NNEWTR, NOLDTR

      REAL(SOP_DP), INTENT(OUT) :: TR1E(LTR1E),  TR1D(LTR1D)
      REAL(SOP_DP), INTENT(IN)  :: EDIA1(LEDIA1), EIVAL1,
     &                             RESI1E(LRESI1E), RESI1D(LRESI1D)
C      
      REAL(SOP_DP) :: DIFF, SUM
      INTEGER :: I
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
         TR1E(I) = RESI1E(I) / SAFE_DENOM(DIFF)
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
         TR1D(I) = RESI1D(I) / SAFE_DENOM(SUM)
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

      CONTAINS

         ! This function ensures that the denominator is at least
         ! "sop_dthresh"         
         pure function safe_denom(x)
            use so_info, only : sop_dthresh
            real(sop_dp), intent(in) :: x
            real(sop_dp) :: safe_denom

            safe_denom = sign(max(abs(x),sop_dthresh),x)
            return
         end function
      END
