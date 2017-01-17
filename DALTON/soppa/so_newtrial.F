C
C  /* Deck so_newtrial */
      SUBROUTINE SO_NEWTRIAL(DOUBLES,NNEWTR,NOLDTR,
     &                       TR1E,LTR1E,TR1D,LTR1D,
     &                       TR2E,LTR2E,TR2D,LTR2D,
     &                       EIVAL1,EDIA1,LEDIA1,EDIA2,
     &                       LEDIA2,SDIA1,LSDIA1,
     &                       RESI1E,LRESI1E,RESI1D,LRESI1D,
     &                       RESI2E,LRESI2E,RESI2D,LRESI2D)
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
C     INPUT:
C        DOUBLES        Include doubles contributions
C        NNEWTR         Number of the created trial-vector in new batch
C        NOLDTR         Number of trial-vectors from previous iterations
C        EIVAL1         Eigenvalue/frequecy corresponding to new vector
C                    Diagonals of elec. hessian and metric:
C        EDIA1(LEDIA1)  E[2] p-h
C        EDIA2(LEDIA2)  E[0] 2p-2h
C        SDIA1(LSDIA1)  S[2] p-h
C                    Residuals of old vectors
C        RESI1E(LRESI1E), RESI1D(LRESI1D)
C        RESI2E(LRESI2E), RESI2D(LRESI2D)
C
C     OUTPUT:
C                    New trial-vectors
C        TR1E(LTR1E), TR1D(LTR1D), TR2E(LTR2E), TR2D(LTR2D)
C
      use so_info, only:  sop_dp 
C
      implicit none
C#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
C
      REAL(sop_dp), PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, 
     &                           ONE = 1.0D0, TWO = 2.0D0
C
C     Array lengths       
      INTEGER, INTENT(IN) ::  LTR1E, LTR1D, LTR2E, LTR2D,
     &                        LEDIA1, LEDIA2, LSDIA1,
     &                        LRESI1E, LRESI1D, LRESI2E, LRESI2D
C
      INTEGER, INTENT(IN) ::  NNEWTR, NOLDTR
C
C     Output arrays. Note that these are also written to disc (is that
C     a good idea, or should it be moved to the caller?)      
      REAL(SOP_DP), INTENT(OUT) :: TR1E(LTR1E),     TR1D(LTR1D),
     &                             TR2E(LTR2E),     TR2D(LTR2D)
C
      REAL(SOP_DP), INTENT(IN)  :: 
     &            EDIA1(LEDIA1),   EDIA2(LEDIA2),   SDIA1(LSDIA1), 
     &            EIVAL1,          RESI1E(LRESI1E), RESI1D(LRESI1D),
     &                             RESI2E(LRESI2E), RESI2D(LRESI2D)
C
      REAL(SOP_DP) :: DIFF
      LOGICAL :: DOUBLES
      INTEGER :: I
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_NEWTRIAL')
C
C----------------------------------------------------------------
C     Calculate 1p1h excitation part of the new raw trial vector.
C----------------------------------------------------------------
C
      DO 100 I = 1,LTR1E
C
         DIFF = EDIA1(I) - EIVAL1 * SDIA1(I)
         TR1E(I) = RESI1E(I)/SAFE_DENOM(DIFF)
C
         DIFF = EDIA1(I) + EIVAL1 * SDIA1(I)
         TR1D(I) = RESI1D(I)/SAFE_DENOM(DIFF)
C
  100 CONTINUE
C
C---------------------------------------------------
C     Write raw new trial vector to file and output.
C---------------------------------------------------
C
      CALL SO_WRITE(TR1E,LTR1E,LUTR1E,FNTR1E,NOLDTR+NNEWTR)
      CALL SO_WRITE(TR1D,LTR1D,LUTR1D,FNTR1D,NOLDTR+NNEWTR)
C
      IF(DOUBLES)THEN
C
C----------------------------------------------------------------
C     Calculate 2p2h excitation part of the new raw trial vector.
C----------------------------------------------------------------
C
         DO 200 I = 1,LTR2E
C
            DIFF = EDIA2(I) - EIVAL1
            TR2E(I) = RESI2E(I) / SAFE_DENOM(DIFF)
C
            DIFF = EDIA2(I) + EIVAL1
            TR2D(I) = RESI2D(I) / SAFE_DENOM(DIFF)
C
  200    CONTINUE
C---------------------------------------------------
C     Write raw new trial vector to file and output.
C---------------------------------------------------
C
         CALL SO_WRITE(TR2E,LTR2E,LUTR2E,FNTR2E,NOLDTR+NNEWTR)
         CALL SO_WRITE(TR2D,LTR2D,LUTR2D,FNTR2D,NOLDTR+NNEWTR)
C
      END IF !DOUBLES
C
      IF ( IPRSOP .GE. 7 ) THEN
C
         CALL AROUND('New raw trialvector in SO_NEWTRIAL')
C
         WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &       (I,TR1E(I),TR1D(I),I=1,LTR1E)
         IF(DOUBLES)THEN
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &           (I,TR2E(I),TR2D(I),I=1,LTR2E)
         ENDIF
C
      END IF

C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_NEWTRIAL')
C
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
