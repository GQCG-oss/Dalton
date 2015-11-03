C
C  /* Deck so_residual */
      SUBROUTINE SO_RESIDUAL(DOUBLES,LTYPE,RESINM,
     &                       RESI1E,LRESI1E,RESI1D,LRESI1D,
     &                       RESI2E,LRESI2E,RESI2D,LRESI2D,
     &                       EIVAL1,
     &                       RES1E,LRES1E,RES1D,LRES1D,
     &                       RES2E,LRES2E,RES2D,LRES2D,
     &                       RESO1E,LRESO1E,RESO1D,LRESO1D,
     &                       RESO2E,LRESO2E,RESO2D,LRESO2D,
     &                       REGP1E,LREGP1E,
     &                       REGP2E,LREGP2E,
     &                       ENORM,IMAGPROP)

C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Andrea Ligabue, January 2004: linear response functions included
C     Rasmus Faber, October 2015: Merged with rp_residual
C         
C     PURPOSE: Calculate the residual vector and norm.
C
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
C
      CHARACTER*6 LTYPE
      DIMENSION RESI1E(LRESI1E), RESI1D(LRESI1D)
      DIMENSION RESI2E(LRESI2E), RESI2D(LRESI2D)
      DIMENSION RES1E(LRES1E),   RES1D(LRES1D)
      DIMENSION RES2E(LRES2E),   RES2D(LRES2D)
      DIMENSION RESO1E(LRESO1E), RESO1D(LRESO1D)
      DIMENSION RESO2E(LRESO2E), RESO2D(LRESO2D)
      DIMENSION REGP1E(LREGP1E)!, REGP1D(LREGP1D)
      DIMENSION REGP2E(LREGP2E)!, REGP2D(LREGP2D)
      LOGICAL IMAGPROP, DOUBLES
      LOGICAL STATIC
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RESIDUAL')
C
      IF ( IPRSOP .GE. 7 ) THEN
C
         CALL AROUND('E[2] linear transformed vector in SO_RESIDUAL')
C
         WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &    (I,RES1E(I),RES1D(I),I=1,LRES1E)
         WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &    (I,RES2E(I),RES2D(I),I=1,LRES2E)
C
         CALL AROUND('S[2] linear transformed vector in SO_RESIDUAL')
C
         WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &    (I,RESO1E(I),RESO1D(I),I=1,LRESO1E)
         WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &    (I,RESO2E(I),RESO2D(I),I=1,LRESO2E)
C
      END IF
C
C----------------------------------------------------------------
C     To handle static properties, ensure that only the E part is
C     non-zero.
C----------------------------------------------------------------
      STATIC = ABS(EIVAL1).LT.1.0D-8
      IF(STATIC)THEN
         IF (LTYPE.NE.'LINEAR') THEN
            WRITE(LUPRI,'(A,F7.4)') 'Warning: Very low excitation'//
     &      ' energy in SO_RESIDUAL: ', EIVAL1
         ENDIF  
         DFACTOR=-1.0D0
         IF(IMAGPROP)DFACTOR=1.0D0
C        Use explicit symmetry in static case         
         CALL DAXPY(LRES1E,DFACTOR,RES1D,1,RES1E,1)
C         IF(DOUBLES) 
C     &      CALL DAXPY(LRES2E,DFACTOR,RES2D,1,RES2E,1)
      ENDIF    
C
C----------------------------------
C     Calculate the residual vector 
C----------------------------------
C      
C     r = (E-wS)x.
C     Where Ex is on RES* and Sx is on RESO*
C
      CALL DCOPY(LRES1E,RES1E,1,RESI1E,1)
      IF(.NOT.STATIC)THEN
         CALL DCOPY(LRES1D,RES1D,1,RESI1D,1)
         CALL DAXPY(LRES1E,-EIVAL1,RESO1E,1,RESI1E,1)
         CALL DAXPY(LRES1D,-EIVAL1,RESO1D,1,RESI1D,1)
      ELSE 
         CALL DZERO(RESI1D,LRES1D)
      ENDIF
C
      IF(DOUBLES)THEN
         CALL DCOPY(LRES2E,RES2E,1,RESI2E,1)
         IF(.NOT.STATIC)THEN
            CALL DCOPY(LRES2D,RES2D,1,RESI2D,1)
            CALL DAXPY(LRES2E,-EIVAL1,RESO2E,1,RESI2E,1)
            CALL DAXPY(LRES2D,-EIVAL1,RESO2D,1,RESI2D,1)
         ELSE
            CALL DZERO(RESI2D,LRES2D)
         ENDIF
      ENDIF
C
      IF(LTYPE.EQ.'LINEAR') THEN
C
C---------------------------------------------------------
C     For linear response properties scale to the original
C     length of the solution vector and subtract the 
C     property gradient vector.
C---------------------------------------------------------
C
          ENORMINV = 1.0D0 / ENORM
          DFACTOR = 1.0D0
          IF(IMAGPROP) DFACTOR = -1.0D0
          CALL DSCAL(LRES1E,ENORMINV,RESI1E,1)
          CALL DAXPY(LRES1E,-1.0D0,REGP1E,1,RESI1E,1)
          IF(.NOT.STATIC) THEN
            CALL DSCAL(LRES1E,ENORMINV,RESI1D,1)
            CALL DAXPY(LRES1E,DFACTOR,REGP1E,1,RESI1D,1)
          ENDIF 
C
          IF(DOUBLES)THEN    
              CALL DSCAL(LRES2E,ENORMINV,RESI2E,1)
              CALL DAXPY(LRES2E,-1.0D0,REGP2E,1,RESI2E,1)
              IF(.NOT.STATIC)THEN
                  CALL DSCAL(LRES2E,ENORMINV,RESI2D,1)
                  CALL DAXPY(LRES2E,DFACTOR,REGP2E,1,RESI2D,1)
              ENDIF
          ENDIF
C
      ENDIF         
C
      IF ( IPRSOP .GE. 7 ) THEN
C
C----------------------------------------
C        Write residual vector to output.
C----------------------------------------
C
         CALL AROUND('Residual vector')
C
         WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &    (I,RESI1E(I),RESI1D(I),I=1,LRESI1E)
         IF(DOUBLES) WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &    (I,RESI2E(I),RESI2D(I),I=1,LRESI2E)
C
      END IF
C
C-----------------------------------------------
C     Calculate the norm of the residual vector.
C-----------------------------------------------
C
C      RESINM = DDOT(LRESI1E,RESI1E,1,RESI1E,1)
c     &       + DDOT(LRESI1D,RESI1D,1,RESI1D,1)
c      IF (DOUBLES) RESINM = RESINM 
C     &       + DDOT(LRESI2E,RESI2E,1,RESI2E,1)
C     &       + DDOT(LRESI2D,RESI2D,1,RESI2D,1)
      
      RESINM = DNRM2(LRESI1E,RESI1E,1)**2
      IF(.NOT.STATIC) RESINM = RESINM + DNRM2(LRESI1D,RESI1D,1)**2
      IF (DOUBLES) THEN 
         RESINM = RESINM + DNRM2(LRESI2E,RESI2E,1)**2
         IF(.NOT.STATIC) RESINM = RESINM + DNRM2(LRESI2D,RESI2D,1)**2
      ENDIF
C
CRF   For static the norm of the D part would be exactly that of the E
C     part. Should we scale RESINM apropriately?
C
      RESINM = DSQRT(RESINM)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RESIDUAL')
C
      RETURN
      END
