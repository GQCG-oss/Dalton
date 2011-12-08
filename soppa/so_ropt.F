#include <single.h>
C
C  /* Deck so_ropt */
      SUBROUTINE SO_ROPT(NOLDTR,NNEWTR,NEXCI,MXDIM,EIVEC,LEIVEC,EXVAL,
     &                   LEXVAL,OVLM,PEIV,REDE,LREDE,REDS,LREDS,LINEAR,
     &                   ENORM,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, April 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Andrea Ligabue, December 2003: LINEAR and ENORM added to the 
C                                    argument list, which are needed
C                                    for the linear equation problem
C
C     PURPOSE: Calculate an orthonormalized set of reduced vectors
C              spanning the space of reduced eigenvectors of the
C              last and the previous iterations in EIVEC.
C
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
#include <ccsdsym.h>
#include <cbiexc.h>
C
      PARAMETER   (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER   (THRNM = 1.0D-10, THRNM2 = 1.0D-12, THREX = 1.0D-5)
      PARAMETER   (T1MIN = 1.0D-8)
      CHARACTER*3 YES, NO
C
      LOGICAL   LINEAR
      DIMENSION EIVEC(LEIVEC,LEIVEC), EXVAL(LEXVAL)
      DIMENSION OVLM(LEIVEC,LEIVEC),  PEIV(LEIVEC)
      DIMENSION REDE(LREDE,LREDE), REDS(LREDS,LREDS)
      DIMENSION WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_ROPT')
C
      LEIVH = LEIVEC/2
      ENORM = ONE
C
C----------------------------------------------------
C     Determine the number of previous trial vectors.
C----------------------------------------------------
C
      NPREV = NOLDTR + NNEWTR - NEXCI
C
C------------------------------------------
C     Write reduced eigenvectors to output.
C------------------------------------------
C
      IF (IPRSOP .GE. 5) THEN
C
         CALL AROUND('Reduced normalized vec. before orthogonalization')
C
         CALL OUTPUT(EIVEC,1,LEIVEC,1,NEXCI+NPREV,LEIVEC,NEXCI+NPREV,
     &               1,LUPRI)
C
      END IF
C
C----------------------------------------------------------------------
C     Check that the NEXCI first reduced trial vectors are orthonormal.
C     They should be orthonormal at this point, small numerical
C     inaccuracies can occur and to minimize these inaccuracies, the
C     vectors are re-orthonormalized.
C----------------------------------------------------------------------
C
      DO 400 IREPET = 1,2
C
         DO 300 JVEC = 1,NEXCI
C
            EVNM = DDOT(LEIVH,EIVEC(1,JVEC),1,EIVEC(1,JVEC),1)
     &           - DDOT(LEIVH,EIVEC(LEIVH+1,JVEC),1,
     &             EIVEC(LEIVH+1,JVEC),1)
C
            IF ( (EVNM - ONE) .GE. THRNM ) THEN
C
               WRITE(LUPRI,'(/,1X,A,I3,A,/,1X,A,F18.9)') 
     &            'WARNING: Reduced trialvector ',
     &            JVEC,' is not properly normalized.',' Norm is',EVNM
C
            END IF
C
cLig        always skip that one in LR becouse nexci is always eq. 1
            DO 200 IVEC = 1,JVEC-1
C
               DOTP1 = DDOT(LEIVH,EIVEC(1,IVEC),1,EIVEC(1,JVEC),1)
     &               - DDOT(LEIVH,EIVEC(LEIVH+1,IVEC),1,
     &                      EIVEC(LEIVH+1,JVEC),1)
C
               CALL DAXPY(LEIVEC,-DOTP1,EIVEC(1,IVEC),1,
     &                    EIVEC(1,JVEC),1)
C
               IF ( (DOTP1) .GE. THRNM ) THEN
C
                  EXDIFF  = ABS( EXVAL(IVEC) - EXVAL(JVEC) )
C
                  IF (EXDIFF .GT. THREX) THEN
C
                     WRITE(LUPRI,
     &                     '(/,1X,A,I3,A,I3,A,/,1X,A,E18.9,/,1X,A)') 
     &               'WARNING: Reduced trialvectors ',IVEC,' and ',JVEC,
     &               ' are not properly orthogonal.',
     &               'Scalar product is',DOTP1,
     &               'It is not a problem if the states are degenerate.'
C
                  END IF
C
               END IF
C
               DOTP2 = DDOT(LEIVH,EIVEC(1+LEIVH,IVEC),1,
     &                      EIVEC(1,JVEC),1)
     &               - DDOT(LEIVH,EIVEC(1,IVEC),1,
     &                      EIVEC(1+LEIVH,JVEC),1)
C
                  CALL DAXPY(LEIVH,DOTP2,EIVEC(1+LEIVH,IVEC),1,
     &                       EIVEC(1,JVEC),1)
                  CALL DAXPY(LEIVH,DOTP2,EIVEC(1,IVEC),1,
     &                       EIVEC(1+LEIVH,JVEC),1)
C
               IF ( (DOTP2) .GE. THRNM ) THEN
C
                  WRITE(LUPRI,'(/,1X,A,I3,A,I3,A,/,1X,A,E18.9)') 
     &              'WARNING: Reduced trialvectors ',IVEC+LEIVH,' and ',
     &               JVEC,' are not properly orthogonal.',
     &               'Scalar product is',DOTP2
C
               END IF
C
  200       CONTINUE
cLig        I start from here again
  201       CONTINUE
C
            EVNM = DDOT(LEIVH,EIVEC(1,JVEC),1,EIVEC(1,JVEC),1)
     &           - DDOT(LEIVH,EIVEC(LEIVH+1,JVEC),1,
     &             EIVEC(LEIVH+1,JVEC),1)
C
            IF (LINEAR .AND. DABS(EVNM).LT.T1MIN) THEN
C                    
                  WRITE(LUPRI,'(/,1X,A,E18.9)') 
     &              'WARNING: Norm is ',EVNM
C
                  WRITE(LUPRI,'(/,1X,A)') 
     &              'I Try to remove the Z part of (Y Z)'
C
                  CALL DZERO(EIVEC(LEIVH+1,JVEC),LEIVH)
C
                  GOTO 201
C
            ELSE
C
              EVNM = ONE / DSQRT( DABS(EVNM) )
C
              CALL DSCAL(LEIVEC,EVNM,EIVEC(1,JVEC),1)
C
C-----------------------------------------------------
C             Keep norm for the linear equation problem.
C-----------------------------------------------------
C
              IF (LINEAR) THEN
                    ENORM = ENORM * EVNM
              ENDIF
C
            ENDIF
C
  300    CONTINUE
C
  400 CONTINUE

C
C==================================================================
C     Set up reduced vectors corresponding to optimized trial
C     vectors of previous iterations in EIVEC. Orthonormalize
C     these reduced vectors over the reduced S[2]. Notice the
C     reduced S[2] matrix is diagonal with plus one in the first
C     half an minus one the last part. The reduced S[2] is therfore
C     not explicitly present in this routine. The procedure
C     is repeated twize to ensure a VERY strict orthonormalization.
C==================================================================
C
      LEXTRA = 1
C
      DO 800 JVEC = NEXCI+1, NEXCI+NPREV
C
C----------------------------------------------------------
C        Set up reduced trial vectors in EIVEC representing
C        previously optimized trial vectors.
C----------------------------------------------------------
C
         CALL DZERO(EIVEC(1,JVEC),LEIVEC)
C
         LPOS = JVEC - NEXCI
C
         EIVEC(LPOS,JVEC) = ONE
C
  500    CONTINUE
C
         DO 700 IRPEAT = 1,3
C
C-----------------------------------------------------------------
C           Orthonormalize reduced vector against previous vectors
C           and their partners.
C-----------------------------------------------------------------
C
            DO 600 IVEC = 1,JVEC-1
C
               DOTP1 = DDOT(LEIVH,EIVEC(1,IVEC),1,EIVEC(1,JVEC),1)
     &               - DDOT(LEIVH,EIVEC(LEIVH+1,IVEC),1,
     &                      EIVEC(LEIVH+1,JVEC),1)
C
               CALL DAXPY(LEIVEC,-DOTP1,EIVEC(1,IVEC),1,
     &                    EIVEC(1,JVEC),1)
C
               DOTP2 = DDOT(LEIVH,EIVEC(1+LEIVH,IVEC),1,
     &                      EIVEC(1,JVEC),1)
     &               - DDOT(LEIVH,EIVEC(1,IVEC),1,
     &                      EIVEC(1+LEIVH,JVEC),1)
C
               CALL DAXPY(LEIVH,DOTP2,EIVEC(1+LEIVH,IVEC),1,
     &                    EIVEC(1,JVEC),1)
               CALL DAXPY(LEIVH,DOTP2,EIVEC(1,IVEC),1,
     &                    EIVEC(1+LEIVH,JVEC),1)
C
  600       CONTINUE
C
C--------------------------------------------------------------------
C           Normalization and a check to ensure linear independent
C           reduced eigenvectors. If linear dependencies are found
C           the reduced vector is reset to an older optimized vector.
C--------------------------------------------------------------------
C
            EVNM = DDOT(LEIVH,EIVEC(1,JVEC),1,EIVEC(1,JVEC),1)
     &           - DDOT(LEIVH,EIVEC(LEIVH+1,JVEC),1,
     &                  EIVEC(LEIVH+1,JVEC),1)
C
            IF ( DABS(EVNM) .LE. THRNM2) THEN
C
               IF ( IPRSOP .GE. 4 ) THEN
                  WRITE(LUPRI,'(1X,A,I3,A)') 'Reduced trialvector #',
     &                  JVEC,' is modified to avoid linear dependencies'
               END IF
C            
               CALL DZERO(EIVEC(1,JVEC),LEIVEC)
C
               LPOS = LPOS + NEXCI
C
               IF (LPOS .GT. (NOLDTR + NNEWTR)) THEN
C
                  IF (LEXTRA .GT. NEXCI ) THEN
                     CALL QUIT('ERROR in SO_ROPT 3: Linear'//
     &               ' dependencies in room of reduced'//
     &               ' trialvectors')
                  END IF 
C
                  LPOS   = NEXCI  + LEXTRA
                  LEXTRA = LEXTRA + 1
C
               END IF
C
               EIVEC(LPOS,JVEC) = ONE
C
               GO TO 500
C
            END IF
C
            IF ( EVNM .LT. ZERO ) THEN
C
C----------------------------------------------------
C              Switch X and Y part of reduced vector.
C----------------------------------------------------
C
               DO IROW = 1,LEIVH
                  TEMP                   = EIVEC(IROW,JVEC)
                  EIVEC(IROW,JVEC)       = EIVEC(IROW+LEIVH,JVEC)
                  EIVEC(IROW+LEIVH,JVEC) = TEMP
               END DO
C
            END IF
C
            EVNM = ONE / DSQRT( DABS(EVNM) )
C
            CALL DSCAL(LEIVEC,EVNM,EIVEC(1,JVEC),1)
C
  700    CONTINUE
C
  800 CONTINUE
C
C-----------------------------------------
C     Generate the paired reduced vectors.
C-----------------------------------------
C
      DO 900 IVEC = 1,NEXCI+NPREV
C
         CALL DCOPY(LEIVH,EIVEC(1,IVEC),1,EIVEC(1+LEIVH,IVEC+LEIVH),1)
         CALL DCOPY(LEIVH,EIVEC(1+LEIVH,IVEC),1,EIVEC(1,IVEC+LEIVH),1)
C
  900 CONTINUE
C
C-----------------------------------------------------------------------
C     Print orthonormalized reduced vectors and their overlap to output.
C-----------------------------------------------------------------------
C
      IF (IPRSOP .GE. 5) THEN
C
         CALL AROUND('Reduced vectors after orthonormalization')
C
         CALL OUTPUT(EIVEC,1,LEIVEC,1,NEXCI+NPREV,LEIVEC,NEXCI+NPREV,
     &               1,LUPRI)
C
         DO JVEC = 1,LEIVEC
            DO IVEC = 1,LEIVEC
               OVLM(IVEC,JVEC) = DDOT(LEIVH,EIVEC(1,IVEC),1,
     &                                EIVEC(1,JVEC),1)
     &                         - DDOT(LEIVH,EIVEC(LEIVH+1,IVEC),1,
     &                                EIVEC(LEIVH+1,JVEC),1)
            END DO
         END DO
C
         CALL AROUND('Overlap of orthonormalized reduced vectors')
C
         CALL OUTPUT(OVLM,1,LEIVEC,1,LEIVEC,LEIVEC,LEIVEC,1,LUPRI)
C
      END IF
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      LMAT1  = LEIVEC * LEIVEC
C
      KMAT1   = 1
      KEND1   = KMAT1 + LMAT1
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_ROPT',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_ROPT',' ',KEND1,LWORK)
C
C-------------------------------------------------------------------
C     Transform the reduced E[2] and S[2] matrices to the new basis.
C-------------------------------------------------------------------
C
      CALL DGEMM('T','N',LEIVEC,LEIVEC,LEIVEC,ONE,EIVEC,LEIVEC,
     &           REDE,LEIVEC,ZERO,WORK(KMAT1),LEIVEC)
C
      CALL DGEMM('N','N',LEIVEC,LEIVEC,LEIVEC,ONE,WORK(KMAT1),LEIVEC,
     &           EIVEC,LEIVEC,ZERO,REDE,LEIVEC)
C
      CALL DGEMM('T','N',LEIVEC,LEIVEC,LEIVEC,ONE,EIVEC,LEIVEC,
     &           REDS,LEIVEC,ZERO,WORK(KMAT1),LEIVEC)
C
      CALL DGEMM('N','N',LEIVEC,LEIVEC,LEIVEC,ONE,WORK(KMAT1),LEIVEC,
     &           EIVEC,LEIVEC,ZERO,REDS,LEIVEC)
C
      IF ( IPRSOP .GE. 6 ) THEN
C
C---------------------------------------------
C        Print reduced E[2] and S[2] matrices.
C---------------------------------------------
C
         CALL AROUND('Reduced E[2] Matrix in new reduced basis')
C
         CALL OUTPUT(REDE,1,LEIVEC,1,LEIVEC,LEIVEC,LEIVEC,1,LUPRI)
C
         CALL AROUND('Reduced S[2] Matrix in new reduced basis')
C
         CALL OUTPUT(REDS,1,LEIVEC,1,LEIVEC,LEIVEC,LEIVEC,1,LUPRI)
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_ROPT')
C
C
      RETURN
      END
