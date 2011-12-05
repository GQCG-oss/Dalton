#include <single.h>
C
C  /* Deck so_ordeig */
      SUBROUTINE SO_ORDEIG(NTRIAL,rede,REDS,LREDS,ALFAR,LALFAR,ALFAI,
     &                     LALFAI,BETA,LBETA,EIVEC,LEIVEC,REDST,LREDST,
     &                     WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan Sauer, June 1996
C     based on RSPORD by Jeppe Olsen, November 1984
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C        Analyse and order the eigenvalues and eigenvectors of the 
C        reduced eigenvalue problem.
C
C NTRIAL          total # of trial vectors and paired trial vectors
C REDS            reduced S[2] matrix
C LREDS           row dimension of REDS
C ALFAR           real part of the numerator of the eigenvalues
C LALFAR          dimension of ALFAR
C ALFAI           imaginary part of the numerator of the eigenvalues
C LALFAI          dimension of ALFAI
C BETA            denumerator of the eigenvalues
C LBETA           dimension of BETA
C EIVEC           eigenvectors
C LEIVEC          dimension of EIVEC
C REDST           work array for reduced S[2] matrix in eigenbasis of
C                 reduced E[2] matrix
C LREDST          dimension of REDST
C WORK            remaining work space
C LWORK           dimension of WORK
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
C
C--------------------------------
C     dimensions of the arguments
C--------------------------------
C
      DIMENSION REDS(LREDS,LREDS)
      DIMENSION ALFAR(LALFAR),ALFAI(LALFAI),BETA(LBETA)
      DIMENSION EIVEC(LEIVEC,LEIVEC)
      DIMENSION REDST(LREDST,LREDST)
      DIMENSION WORK(LWORK)
C
      PARAMETER (ZERO   = 0.0D+00, ONE    = 1.0D+00)
      PARAMETER (SMALL  = 1.0D-09, SMALL2 = 1.0D-08)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_ORDEIG')
C
C================================
C     Analyze reduced eigenvalues
C================================
C
      IF ( IPRSOP .GE. 5 ) THEN
C
C------------------------------------------------------
C        Print reduced eigenvalues and eigenvectors
C        Eigenvalue k is (ALFAR(k)+i*ALFAI(k))/BETA(k).
C------------------------------------------------------
C
         CALL AROUND('Eigenvalues of reduced E[2] before Ordering')
C
         WRITE(LUPRI,'(A,//,A,/)')
     &         '  (ALFAR + i ALFAI) / BETA are eigenvalues;',
     &         '     ALFAR           ALFAI          BETA'
         WRITE(LUPRI,'(1P,3D15.6)')
     &         (ALFAR(I),ALFAI(I),BETA(I),I=1,NTRIAL)
C
         CALL AROUND('Reduced Eigenvectors before Ordering')
C
         CALL OUTPUT(EIVEC,1,NTRIAL,1,NTRIAL,LEIVEC,LEIVEC,1,LUPRI)
C
      END IF
C
C-----------------------------------------------------------------------
C     Calculate eigenvalues of reduced E[2] and check for singularities
C     and complex eigenvalues.
C-----------------------------------------------------------------------
C
      ICMPLX = 0
      ISING  = 0
      DO I = 1, NTRIAL
         IF (ABS(BETA(I)) .GE. SMALL) THEN
            ALFAR(I) = ALFAR(I) / BETA(I)
            ALFAI(I) = ALFAI(I) / BETA(I)
         ELSE
            ISING = ISING + 1
            WRITE(LUPRI,'(/,A,/,A,I6,1P,3D13.6)') 
     &            ' *** Singularity in reduced eigenvalue problem :',
     &            '     I,ALFAR(I),ALFAI(I),BETA(I):',
     &                  I, ALFAR(I),ALFAI(I),BETA(I)
         END IF
C
         IF (ABS(ALFAI(I)) .GT. ZERO) THEN
            ICMPLX = ICMPLX + 1
            WRITE(LUPRI,'(/,A,A,/,A,I6,1P,2D13.6)') 
     &            ' *** Complex eigenvalue in reduced eigenvalue ',
     &            'problem :','     Real and imaginary part :  ',
     &            I,ALFAR(I),ALFAI(I)
C
         END IF
      END DO
      IF ((ISING .GT. 0) .OR. (ICMPLX .GT. 0)) THEN
         WRITE(LUPRI,*) 'SO_ORDEIG: ',ISING,' singular and ',ICMPLX,
     &              ' complex eigenvalues in reduced eigenvalue problem'
      ENDIF
C
C-------------------------------
C     Print reduced eigenvalues.
C-------------------------------
C
      IF ((IPRSOP .GE. 5) .OR. (ICMPLX .GE. 1) .OR. (ISING .GE. 1)) THEN
C
         DO I = 1,NTRIAL
C
            IF (ABS(BETA(I)) .GE. SMALL) THEN
C
               IF (ABS(ALFAI(I)) .GT. ZERO) THEN
                  WRITE(LUPRI,'(1X,A,I4,2P,2D14.6)') 
     &                  'I, ALFAR(I), ALFAI(I)         :',
     &                   I, ALFAR(I), ALFAI(I)
               ELSE
                  WRITE(LUPRI,'(1X,A,I4,2P,D14.6)')
     &                  'I, ALFAR(I)                   :',
     &                   I, ALFAR(I)
               END IF
C
            ELSE
C
               WRITE(LUPRI,'(1X,A,I4,2P,3D14.6)') 
     &               'I, ALFAR(I), ALFAI(I), BETA(I):',
     &                I, ALFAR(I), ALFAI(I), BETA(I)
            END IF
C
         END DO
C
      END IF
C
C=======================
C     Order eigenvectors
C=======================
C
C-----------------------------------------------------------------------
C     Transform reduced S[2] matrix to eigenbasis of reduced E[2] matrix
C-----------------------------------------------------------------------
C
      LSCR1 = NTRIAL
C
      KSCR1   = 1
      KEND1   = KSCR1  + LSCR1*LSCR1
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_ORDEIG.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_ORDEIG.1',' ',KEND1,LWORK)
C
      CALL DGEMM('N','N',NTRIAL,NTRIAL,NTRIAL,ONE,REDS,LREDS,
     &           EIVEC,LEIVEC,ZERO,WORK(KSCR1),LSCR1)
      CALL DGEMM('T','N',NTRIAL,NTRIAL,NTRIAL,ONE,EIVEC,LEIVEC,
     &           WORK(KSCR1),LSCR1,ZERO,REDST,LREDST)
C
      IF ( IPRSOP .GE. 6) THEN
C
         CALL AROUND('Reduced S[2] matrix in eigenbasis of reduced '//
     &            'E[2] matrix')
C
         CALL OUTPUT(REDST,1,NTRIAL,1,NTRIAL,LREDST,LREDST,1,LUPRI)
C
      END IF
C
C---------------------------------------------------------------
C     Count eigenvectors with positive, negative and zero metric
C     and put the ones with positive metric first
C---------------------------------------------------------------
C
      IPOS = 0
      IZER = 0
      INEG = 0
C
      DO I = 1, NTRIAL
C
         IF (BETA(I) .LE. SMALL) THEN
C
            IZER = IZER + 1
C
         ELSE IF (REDST(I,I) .GT. ZERO) THEN
C
            IPOS  = IPOS + 1
            SCALE = ONE / DSQRT(REDST(I,I))
            CALL DSCAL(NTRIAL,SCALE,EIVEC(1,I),1)
C
            IF (IPOS .NE. I) THEN
C
               CALL DSWAP(NTRIAL,EIVEC(1,I),1,EIVEC(1,IPOS),1)
C
               SAVE             = ALFAR(IPOS)
               ALFAR(IPOS)      = ALFAR(I)
               ALFAR(I)         = SAVE
C
               SAVE             = ALFAI(IPOS)
               ALFAI(IPOS)      = ALFAI(I)
               ALFAI(I)         = SAVE
C
               SAVE             = REDST(IPOS,IPOS)
               REDST(IPOS,IPOS) = REDST(I,I)
               REDST(I,I)       = SAVE
C
               SAVE             = BETA(IPOS)
               BETA(IPOS)       = BETA(I)
               BETA(I)          = SAVE
C
            END IF
C
         ELSE IF (REDST(I,I) .LT. - ZERO) THEN
C
            INEG  = INEG + 1
            SCALE = ONE / DSQRT(ABS(REDST(I,I)))
            CALL DSCAL(NTRIAL,SCALE,EIVEC(1,I),1)
C
         END IF
C
      END DO
C
      NNEG = IPOS
C
      DO I = IPOS + 1, NTRIAL
C
         IF (ABS(BETA(I)) .GT. SMALL) THEN
C
            NNEG = NNEG + 1
C
            IF (NNEG .NE. I) THEN
C
               CALL DSWAP(NTRIAL,EIVEC(1,I),1,EIVEC(1,NNEG),1)
C
               SAVE             = ALFAR(NNEG)
               ALFAR(NNEG)      = ALFAR(I)
               ALFAR(I)         = SAVE
C
               SAVE             = ALFAI(NNEG)
               ALFAI(NNEG)      = ALFAI(I)
               ALFAI(I)         = SAVE
C
               SAVE             = REDST(NNEG,NNEG)
               REDST(NNEG,NNEG) = REDST(I,I)
               REDST(I,I)       = SAVE
C
               SAVE             = BETA(NNEG)
               BETA(NNEG)       = BETA(I)
               BETA(I)          = SAVE
C
            END IF
C
         END IF
C
      END DO
C
      IF ( IPRSOP .GE. 4 ) THEN
C
         WRITE(LUPRI,'(3(/,A,I6))') 
     &          ' *** Eigenvectors with positive metric:',IPOS,
     &          '     Eigenvectors with zero metric:    ',IZER,
     &          '     Eigenvectors with negative metric:',INEG
C
      END IF
C
C------------------------------------------------------------
C     Order eigensolutions in ascending order of eigenvalues.
C------------------------------------------------------------
C
      DO I = 1, IPOS
         JMIN = I
         AMIN = ALFAR(I)
         DO J = I+1, IPOS
            IF(ALFAR(J) .LT. AMIN) THEN
               AMIN = ALFAR(J)
               JMIN = J
            ENDIF
         END DO
         IF (JMIN .NE. I) THEN
            ALFAR(JMIN) = ALFAR(I)
            ALFAR(I)    = AMIN
            CALL DSWAP(NTRIAL,EIVEC(1,I),1,EIVEC(1,JMIN),1)
         ENDIF
      END DO
C
      DO I = IPOS+1, IPOS+INEG
         JMIN = I
         AMIN = ALFAR(I)
         DO J = I+1, IPOS+INEG
            IF(ALFAR(J) .GT. AMIN) THEN
               AMIN = ALFAR(J)
               JMIN = J
            ENDIF
         END DO
         IF (JMIN .NE. I) THEN
            ALFAR(JMIN) = ALFAR(I)
            ALFAR(I)    = AMIN
            CALL DSWAP(NTRIAL,EIVEC(1,I),1,EIVEC(1,JMIN),1)
         ENDIF
      END DO
C
C------------------------------------------------
C     Print reduced eigenvalues and eigenvectors.
C------------------------------------------------
C
      IF (INEG .NE. IPOS) THEN
         WRITE(LUPRI,'(3(/A))')' ***** WARNING *********'
     &   ,' number of eigenvalues with negative metric differ from'
     &   ,' number with positive metric'
         WRITE(LUPRI,'(/A)')'   Number    Eigenvalue '
         DO I = 1, NTRIAL
            WRITE(LUPRI,'(I10,5X,1P,D15.8)') I,ALFAR(I)
         END DO
      ELSE
C
         IF ( IPRSOP .GE. 4 ) THEN
C
            CALL AROUND('Eigenvalues of reduced eigenvalue problem')
            WRITE(LUPRI,'(/A)') '    Number          Eigenvalue    '//
     &                      'Paired Eigenvalue'
C
         END IF
C
            DO I = 1, IPOS
               IF (ABS(ABS(ALFAR(I))-ABS(ALFAR(IPOS+I))).GT.SMALL2) THEN
                  WRITE(LUPRI,'(/A)') 
     &                  ' **WARNING** Eigenvalues not paired'
                  WRITE(LUPRI,'(I10,5X,1P,D15.8,6X,D15.8)')
     &                            I,ALFAR(I),ALFAR(IPOS+I)
               END IF
            END DO
C
         IF (IZER .GT. 0) THEN
            WRITE(LUPRI,'(/A)')
     &            ' Zero Metric Eigenvalue and Zero Eigenvalue'
            WRITE(LUPRI,'(/A)') 'Number    Eigenvalue'
            DO I = 1, IZER
               WRITE(LUPRI,'(I10,5X,1P,D15.8)') I,ALFAR(IPOS+INEG+I)
            END DO
         END IF
      END IF
C
      IF ( IPRSOP .GE. 5 ) THEN
C
         CALL AROUND('Reduced Eigenvectors after Ordering')
C
         CALL OUTPUT(EIVEC,1,NTRIAL,1,NTRIAL,LEIVEC,LEIVEC,1,LUPRI)
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_ORDEIG')
C
      RETURN
      END
