C
C  /* Deck so_test1 */
      SUBROUTINE SO_TEST1(TR1EI,TR1DI,TR2EI,TR2DI,TR1EJ,TR1DJ,TR2EJ,
     &                    TR2DJ,RES1EI,RES1DI,RES2EI,RES2DI,RES1EJ,
     &                    RES1DJ,RES2EJ,RES2DJ,
     &                    RESO1EI, RESO1DI, RESO1EJ, RESO1DJ,
     &                    LTR1E,LTR1D,LTR2E,LTR2D,
     &                    LRES1E,LRES1D,
     &                    LRES2E,LRES2D,LRESO1E,LRESO1D,NOLDTR,NNEWTR,
     &                    DENSIJ,  LDENSIJ,  DENSAB, LDENSAB,
     &                    T2MP,    LT2MP,   FOCKD,    LFOCKD,
     &                    ISYMTR,  WORK,    LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C 
C     Keld Bak, August 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Test of orhtonormality of trialvectors and the
C              correctnes of the linear transformed vectors 
C              which are generated as linear combinations of
C              previously calculated linear transformed vectors.
C
#include <implicit.h>
#include <priunit.h>
#include <soppinf.h>
#include <ccsdsym.h>
C
      PARAMETER (ONE = 1.0D0)
C
      DIMENSION TR1EI(LTR1E), TR1DI(LTR1E), TR2EI(LTR2E), TR2DI(LTR2D)
      DIMENSION TR1EJ(LTR1E), TR1DJ(LTR1E), TR2EJ(LTR2E), TR2DJ(LTR2D)
      DIMENSION RES1EI(LTR1E),RES1DI(LTR1D),RES2EI(LTR2E),RES2DI(LTR2D)
      DIMENSION RES1EJ(LTR1E),RES1DJ(LTR1D),RES2EJ(LTR2E),RES2DJ(LTR2D)
      DIMENSION RESO1EI(LTR1E),  RESO1DI(LTR1D)
      DIMENSION RESO1EJ(LTR1E),  RESO1DJ(LTR1D)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION T2MP(LT2MP),     FOCKD(LFOCKD)
      DIMENSION WORK(LWORK)
      LOGICAL   PROBL
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TEST1')
C
      PROBL = .FALSE.
C
C-------------------------------------------------------
C     Open files with trial and linear transformed trial
C     (solution) vectors.
C-------------------------------------------------------
C
      CALL SO_OPEN(LUTR1E, FNTR1E, LTR1E)
      CALL SO_OPEN(LUTR1D, FNTR1D, LTR1D)
      CALL SO_OPEN(LUTR2E, FNTR2E, LTR2E)
      CALL SO_OPEN(LUTR2D, FNTR2D, LTR2D)
      CALL SO_OPEN(LURS1E, FNRS1E, LRES1E)
      CALL SO_OPEN(LURS1D, FNRS1D, LRES1D)
      CALL SO_OPEN(LURS2E, FNRS2E, LRES2E)
      CALL SO_OPEN(LURS2D, FNRS2D, LRES2D)
      CALL SO_OPEN(LURO1E, FNRO1E, LRESO1E)
      CALL SO_OPEN(LURO1D, FNRO1D, LRESO1D)
C
      DO 100 ITRIAL = 1,NOLDTR+NNEWTR
C
         CALL SO_READ(TR1EI, LTR1E, LUTR1E, FNTR1E, ITRIAL)
         CALL SO_READ(TR1DI, LTR1D, LUTR1D, FNTR1D, ITRIAL)
         CALL SO_READ(TR2EI, LTR2E, LUTR2E, FNTR2E, ITRIAL)
         CALL SO_READ(TR2DI, LTR2D, LUTR2D, FNTR2D, ITRIAL)

         DO 200 JTRIAL = 1,ITRIAL
C
            CALL SO_READ(TR1EJ, LTR1E, LUTR1E, FNTR1E, JTRIAL)
            CALL SO_READ(TR1DJ, LTR1D, LUTR1D, FNTR1D, JTRIAL)
            CALL SO_READ(TR2EJ, LTR2E, LUTR2E, FNTR2E, JTRIAL)
            CALL SO_READ(TR2DJ, LTR2D, LUTR2D, FNTR2D, JTRIAL)
C
            DOTP = DDOT(LTR1E,TR1EI,1,TR1EJ,1)
     &           + DDOT(LTR1D,TR1DI,1,TR1DJ,1) 
     &           + DDOT(LTR2E,TR2EI,1,TR2EJ,1) 
     &           + DDOT(LTR2D,TR2DI,1,TR2DJ,1) 
C
            IF ( ITRIAL .EQ. JTRIAL ) THEN
C
               IF ( DABS( DOTP - ONE ) .GT. 1.D-14 ) THEN
                  WRITE(LUPRI,*) 'WARNING: Norm of trial vector',
     &                        ITRIAL,' is ',DOTP
                  PROBL = .TRUE. 
               END IF
C
            ELSE
C
               IF ( DABS( DOTP ) .GT. 1.D-14 ) THEN
                  WRITE(LUPRI,*) 'WARNING: Dotp of trial vectors',
     &                        ITRIAL,' and ',JTRIAL,' is ',DOTP
                  PROBL = .TRUE. 
               END IF
C
            END IF
C
            DOTP = DDOT(LTR1E,TR1EI,1,TR1DJ,1)
     &           + DDOT(LTR1D,TR1DI,1,TR1EJ,1) 
     &           + DDOT(LTR2E,TR2EI,1,TR2DJ,1) 
     &           + DDOT(LTR2D,TR2DI,1,TR2EJ,1) 
C
            IF ( DABS( DOTP ) .GT. 1.D-14 ) THEN
               WRITE(LUPRI,*) 'WARNING: Dotp of trial vector',
     &                     ITRIAL,' and the paired',JTRIAL,' is ',DOTP
               PROBL = .TRUE. 
            END IF
C
  200    CONTINUE
C
         IF ( .NOT. PROBL ) THEN
            WRITE(LUPRI,'(2X,A,/,11X,A)') 'SO_TEST: No problems'//
     &            ' concerning orthonormality of trial vectors',
     &            'have been found.'
         END IF
C
C-----------------------------------------------
C        Read the E[2] transformed trial vectors
C-----------------------------------------------
C
         CALL SO_READ(RES1EJ,  LRES1E, LURS1E, FNRS1E, ITRIAL)
         CALL SO_READ(RES1DJ,  LRES1D, LURS1D, FNRS1D, ITRIAL)
         CALL SO_READ(RES2EJ,  LRES2E, LURS2E, FNRS2E, ITRIAL)
         CALL SO_READ(RES2DJ,  LRES2E, LURS2D, FNRS2D, ITRIAL)
         CALL SO_READ(RESO1EJ, LRES1E, LURO1E, FNRO1E, ITRIAL)
         CALL SO_READ(RESO1DJ, LRES1D, LURO1D, FNRO1D, ITRIAL)
C
         CALL DZERO(RES1EI, LTR1E)
         CALL DZERO(RES1DI, LTR1D)
         CALL DZERO(RES2EI, LTR2E)
         CALL DZERO(RES2DI, LTR2D)
         CALL DZERO(RESO1EI,LTR1E)
         CALL DZERO(RESO1DI,LTR1D)
C
         CALL SO_RES_CP(1,       DENSIJ,  LDENSIJ,  DENSAB, LDENSAB,
     &                  T2MP,    LT2MP,   FOCKD,    LFOCKD,
     &                  ISYMTR,  1,       WORK,     LWORK,
     &                  TR1EI,   TR1DI,   TR2EI,    TR2DI,
     &                  RES1EI,  RES1DI,  RES2EI,   RES2DI,
     &                  RESO1EI, RESO1DI,
     &                  LTR1E,   LTR1D,   LTR2E,    LTR2D,
     &                  LRES1E,  LRES1D,  LRES2E,   LRES2D,
     &                  LRESO1E, LRESO1D)
C
         INR = 0
C
         DO 301 I = 1,LTR1E
            IF ( DABS(RES1EI(I) - RES1EJ(I)) .GT. 1.D-14 ) THEN
               INR = INR + 1
               WRITE(LUPRI,*) I, RES1EI(I), RES1EJ(I)
            END IF
  301    CONTINUE
         DO 302 I = 1,LTR2E
            IF ( DABS(RES2EI(I) - RES2EJ(I)) .GT. 1.D-14 ) THEN
               INR = INR + 1
               WRITE(LUPRI,*) I, RES2EI(I), RES2EJ(I)
            END IF
  302    CONTINUE
         DO 303 I = 1,LTR1D
            IF ( DABS(RES1DI(I) - RES1DJ(I)) .GT. 1.D-14 ) THEN
               INR = INR + 1
               WRITE(LUPRI,*) I, RES1DI(I), RES1DJ(I)
            END IF
  303    CONTINUE
         DO 304 I = 1,LTR2D
            IF ( DABS(RES2DI(I) - RES2DJ(I)) .GT. 1.D-14 ) THEN
               INR = INR + 1
               WRITE(LUPRI,*) I, RES2DI(I), RES2DJ(I)
            END IF
  304    CONTINUE
C
         IF (INR .GT. 0) THEN
            WRITE(LUPRI,*) 'WARNING:',INR,' of the ',ITRIAL,
     &                     ' E[2] linear'//
     &                     'transformed vector elements seems off'
         END IF
C
         INS = 0
C
         DO 401 I = 1,LTR1E
            IF ( DABS(RESO1EI(I) - RESO1EJ(I)) .GT. 1.D-14 ) THEN
               INS = INS +1
               WRITE(LUPRI,*) I, RESO1EI(I), RESO1EJ(I)
            END IF
  401    CONTINUE
         DO 402 I = 1,LTR1D
            IF ( DABS(RESO1DI(I) - RESO1DJ(I)) .GT. 1.D-14 ) THEN
               INS = INS +1
               WRITE(LUPRI,*) I, RESO1EI(I), RESO1EJ(I)
            END IF
  402    CONTINUE
C
         IF (INS .GT. 0) THEN
            WRITE(LUPRI,*) 'WARNING:',INS,' of the ',ITRIAL,
     &                     ' S[2] linear'//
     &                     'transformed vector elements seems off'
         END IF
C
         IF ( IPRSOP .GE. 7 ) THEN
C
            WRITE(LUPRI,*) ' Result-vector, calculated and determined'
            WRITE(LUPRI,9002) (I,RES1EI(I),RES1EJ(I),I=1,LTR1E)
            WRITE(LUPRI,9002) (I,RES2EI(I),RES2EJ(I),I=1,LTR2E)
            WRITE(LUPRI,9002) (I,RES1DI(I),RES1DJ(I),I=1,LTR1D)
            WRITE(LUPRI,9002) (I,RES2DI(I),RES2DJ(I),I=1,LTR2D)
            WRITE(LUPRI,*) ' Overlap-vector, calculated and determined'
            WRITE(LUPRI,9002) (I,RESO1EI(I),RESO1EJ(I),I=1,LTR1E)
            WRITE(LUPRI,9002) (I,RESO1DI(I),RESO1DJ(I),I=1,LTR1D)
C
         END IF
C
         IF ( INS. EQ. 0 .AND. INR .EQ. 0 ) THEN
            WRITE(LUPRI,'(/,2X,A,/,11X,A,/)') 'SO_TEST: No problems '//
     &            'concerning linear transformed trial vectors',
     &            'have been found.'
         END IF
C
  100 CONTINUE
C
C-------------------------------------------------------
C     Close files with trial and linear transformed trial
C     (solution) vectors.
C-------------------------------------------------------
C
      CALL SO_CLOSE(LUTR1E, FNTR1E, 'KEEP')
      CALL SO_CLOSE(LUTR1D, FNTR1D, 'KEEP')
      CALL SO_CLOSE(LUTR2E, FNTR2E, 'KEEP')
      CALL SO_CLOSE(LUTR2D, FNTR2D, 'KEEP')
      CALL SO_CLOSE(LURS1E, FNRS1E, 'KEEP')
      CALL SO_CLOSE(LURS1D, FNRS1D, 'KEEP')
      CALL SO_CLOSE(LURS2E, FNRS2E, 'KEEP')
      CALL SO_CLOSE(LURS2D, FNRS2D, 'KEEP')
      CALL SO_CLOSE(LURO1E, FNRO1E, 'KEEP')
      CALL SO_CLOSE(LURO1D, FNRO1D, 'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TEST1')
C
      RETURN
C
 9000 FORMAT(A30,I3,A5,I3,A8,F16.12)
 9001 FORMAT(A30,I3,A10,I3,A3,F16.12)
 9002 FORMAT(I3,1X,F14.9,1X,F14.9)
C
      END
