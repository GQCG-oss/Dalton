C
C  /* Deck rp_residual */
      SUBROUTINE RP_RESIDUAL(LTYPE,RESINM,RESI1E,LRESI1E,RESI1D,LRESI1D,
     &                       EIVAL1,RES1E,LRES1E,RES1D,LRES1D,RESO1E,
     &                       LRESO1E,RESO1D,LRESO1D,REGP1E,LREGP1E,
     &                       REGP1D,LREGP1D,ENORM)

C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Andrea Ligabue, December 2003: Linear response functions included
C
C     PURPOSE: Calculate the residual vector and norm.
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
C
      CHARACTER*6 LTYPE
      DIMENSION RESI1E(LRESI1E), RESI1D(LRESI1D)
      DIMENSION RES1E(LRES1E),   RES1D(LRES1D)
      DIMENSION RESO1E(LRESO1E), RESO1D(LRESO1D)
      DIMENSION REGP1E(LREGP1E), REGP1D(LREGP1D)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_RESIDUAL')
C
      IF ( IPRSOP .GE. 7 ) THEN
C
         CALL AROUND('E[2] linear transformed vector in RP_RESIDUAL')
C
         WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &        (I,RES1E(I),RES1D(I),I=1,LRES1E)
C
         CALL AROUND('S[2] linear transformed vector in RP_RESIDUAL')
C
         WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &        (I,RESO1E(I),RESO1D(I),I=1,LRESO1E)
C
      END IF
C
C-----------------------------------
C     Calculate the residual vector.
C-----------------------------------
C
      CALL DCOPY(LRES1E,RES1E,1,RESI1E,1)
      CALL DCOPY(LRES1D,RES1D,1,RESI1D,1)
C
      CALL DAXPY(LRES1E,-EIVAL1,RESO1E,1,RESI1E,1)
      CALL DAXPY(LRES1D,-EIVAL1,RESO1D,1,RESI1D,1)
C
C-------------------------------------------
C     For linear response properties
C     subtract the property gradient vector.
C-------------------------------------------
C
      IF (LTYPE .EQ. 'LINEAR') THEN
C
         ENORM = 1.0D0 / ENORM
         CALL DSCAL(LRES1E,ENORM,RESI1E,1)
         CALL DSCAL(LRES1D,ENORM,RESI1D,1)
C
         CALL DAXPY(LRES1E,-1.0D0,REGP1E,1,RESI1E,1)
         CALL DAXPY(LRES1D,-1.0D0,REGP1D,1,RESI1D,1)
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
     &        (I,RESI1E(I),RESI1D(I),I=1,LRESI1E)
C
      END IF
C
C-----------------------------------------------
C     Calculate the norm of the residual vector.
C-----------------------------------------------
C
      RESINM = DDOT(LRESI1E,RESI1E,1,RESI1E,1)
     &       + DDOT(LRESI1D,RESI1D,1,RESI1D,1)
C
      RESINM = DSQRT(RESINM)
C
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_RESIDUAL')
      RETURN
      END
