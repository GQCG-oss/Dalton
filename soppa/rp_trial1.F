C
C  /* Deck rp_trial1 */
      SUBROUTINE RP_TRIAL1(NNEWTR,POINT,LPOINT,ISYMTR,NEXCI,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Determine the initial trialvectors.
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
#include <ccsdsym.h>
C
      PARAMETER (ONE = 1.0D0, STHR = 1.0D-5)
C
      INTEGER   POINT(LPOINT)
      DIMENSION WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_TRIAL1')
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LEDIA1 = NT1AM(ISYMTR)
C
      KEDIA1  = 1
      KEND1   = KEDIA1 + LEDIA1
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_TRIAL1.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_TRIAL1.1',' ',KEND1,LWORK)
C
C---------------------------------------------------------------
C     Read diagonal E[2] elements which are taken as approximate
C     eigenvalues.
C---------------------------------------------------------------
C
      CALL GPOPEN(LUDIAG,'RP_DIAG','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUDIAG
C
      READ(LUDIAG) ( WORK(KEDIA1+I-1), I = 1,LEDIA1)
C
      CALL GPCLOSE(LUDIAG,'KEEP')
C
C---------------------------------
C     2. allocation of work space.
C---------------------------------
C
      LPARRA = LEDIA1
C
      KPARRA  = KEND1
      KEND2   = KPARRA + LPARRA
      LWORK2  = LWORK  - KEND2
C
      CALL SO_MEMMAX ('RP_TRIAL1.1',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('RP_TRIAL1.2',' ',KEND2,LWORK)
C
C--------------------------------------------------
C     Find the NEXCI lowest approximate eigenvalues
C     and leave pointers to them in POINT.
C--------------------------------------------------
C
      CALL SO_SORT(POINT,NEXCI,WORK(KEDIA1),LEDIA1,WORK(KPARRA))
C
C---------------------------------
C     3. allocation of work space.
C---------------------------------
C
      LTRIAL = LEDIA1
C
      KTRIAL  = 1
      KEND3   = KTRIAL + LTRIAL
      LWORK3  = LWORK  - KEND3
C
      CALL SO_MEMMAX ('RP_TRIAL1.1',LWORK3)
      IF (LWORK3 .LT. 0) CALL STOPIT('RP_TRIAL1.3',' ',KEND3,LWORK)
C
C-----------------------------------------------------------
C     Set the number of new trial vectors equal to number of 
C     excitations.
C-----------------------------------------------------------
C
      NNEWTR = NEXCI
C
C---------------------------------------------------------
C     Create initial new trial vectors and write to files.
C---------------------------------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LEDIA1)
      CALL SO_OPEN(LUTR1D,FNTR1D,LEDIA1)
C
      DO 200 INEWTR = 1, NNEWTR
C
         CALL DZERO(WORK(KTRIAL),LTRIAL)
C
         CALL SO_WRITE(WORK(KTRIAL),LEDIA1,LUTR1D,FNTR1D,INEWTR)
C
         WORK( KTRIAL + POINT(INEWTR) - 1 ) = ONE
C
         CALL SO_WRITE(WORK(KTRIAL),LEDIA1,LUTR1E,FNTR1E,INEWTR)
C
  200 CONTINUE
C
      IF ( IPRSOP .GE. 6 ) THEN
C
C------------------------------------------
C        Write new trial vectors to output.
C------------------------------------------
C
         DO 300 INEWTR = 1,NNEWTR
C
            WRITE(LUPRI,'(I3,A)') INEWTR,
     &         '. raw trial vector in RP_TRIAL1'
C
            CALL SO_READ(WORK(KTRIAL),LEDIA1,LUTR1E,FNTR1E,INEWTR)
            WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTRIAL+I-1),I=1,LEDIA1)
            CALL SO_READ(WORK(KTRIAL),LEDIA1,LUTR1D,FNTR1D,INEWTR)
            WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTRIAL+I-1),I=1,LEDIA1)
C
  300    CONTINUE
C
      END IF
C
C-----------------------------------------------
C     Orthogonalize new trial vectors over S[2].
C-----------------------------------------------
C
      NLINDP = 0
      NOLDTR = 0
      CALL RP_ORTH_TRN('EXCITA',NOLDTR,NNEWTR,NLINDP,ISYMTR,WORK,LWORK)
C
      IF ( IPRSOP .GE. 6 ) THEN
C
C------------------------------------------
C        Write new trial vectors to output.
C------------------------------------------
C
         DO 400 INEWTR = 1,NNEWTR
C
            WRITE(LUPRI,'(I3,A)') INEWTR,
     &         '. new trial vector from RP_TRIAL1'
C
            CALL SO_READ(WORK(KTRIAL),LEDIA1,LUTR1E,FNTR1E,INEWTR)
            WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTRIAL+I-1),I=1,LEDIA1)
            CALL SO_READ(WORK(KTRIAL),LEDIA1,LUTR1D,FNTR1D,INEWTR)
            WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTRIAL+I-1),I=1,LEDIA1)
C
  400    CONTINUE
C
      END IF
C
C------------------------------------
C     Close files with trial vectors.
C------------------------------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_TRIAL1')
C
      RETURN
C
      END
