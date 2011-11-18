#include <single.h>
C
C  /* Deck so_trial1 */
      SUBROUTINE SO_TRIAL1(MODEL,NNEWTR,POINT,LPOINT,ISYMTR,NEXCI,
     &                     DENSIJ,LDENSIJ,DENSAB,LDENSAB,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Determine the initial trialvectors.
CW
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
#include <ccsdsym.h>
C
      PARAMETER (ONE = 1.0D0, STHR = 1.0D-5)
C
      INTEGER     POINT(LPOINT)
      DIMENSION   DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION   WORK(LWORK)
      CHARACTER*5 MODEL
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TRIAL1')
C
C-----------------------------------------------------------
C     Set the number of new trial vectors equal to number of
C     excitations.
C-----------------------------------------------------------
C
      NNEWTR = NEXCI
C
C--------------------------------------------
C     Open files for storing of trialvectors.
C--------------------------------------------
C
      LTR1 = NT1AM(ISYMTR)
      LTR2 = NT2AM(ISYMTR)
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1)
      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2)
      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2)
C
C=============================================================
C     Determine how the initial trialvectors are to be chosen.
C=============================================================
C
      IF ( (MODEL .EQ. 'AOSOC') .AND. AOSOP ) THEN
C
C---------------------------------------------------------------------
C        Use SOPPA eigenvectors as first trial vectors in SOPPA(CCSD).
C        Since the SOPPA vectors are in place nothing shall be done.
C---------------------------------------------------------------------
C
         CONTINUE
C
      ELSE IF ((MODEL.EQ.'AOSOC') .AND. DCRPA .AND. (NEXCI.LE.LTR1))THEN
C
C----------------------------------------------------------------------
C        Use RPA(D) eigenvectors as first trial vectors in SOPPA(CCSD).
C        Since the RPA(D) vectors are in place nothing shall be done.
C----------------------------------------------------------------------
C
         CONTINUE
C
      ELSE IF ((MODEL.EQ.'AOSOC') .AND. AORPA .AND. (NEXCI.LE.LTR1))THEN
C
C--------------------------------------------------------------------
C        Use RPA eigenvectors as first trial vectors for SOPPA(CCSD).
C--------------------------------------------------------------------
C
         CALL SO_TR1RP(NNEWTR,LTR2,WORK,LWORK) 
C
      ELSE IF ((MODEL.EQ.'AOSOP') .AND. DCRPA .AND. (NEXCI.LE.LTR1))THEN
C
C--------------------------------------------------------------------
C        Use RPA(D) eigenvectors as first trial vectors in SOPPA.
C        Since the RPA(D) vectors are in place nothing shall be done.
C--------------------------------------------------------------------
C
         CONTINUE
C
      ELSE IF ((MODEL.EQ.'AOSOC') .AND. AORPA .AND. (NEXCI.LE.LTR1))THEN
C
C--------------------------------------------------------------
C        Use RPA eigenvectors as first trial vectors for SOPPA.
C--------------------------------------------------------------
C
         CALL SO_TR1RP(NNEWTR,LTR2,WORK,LWORK) 
C
      ELSE IF ((MODEL.EQ.'AOSOC') .OR. (MODEL.EQ.'AOSOP')) THEN
C
C---------------------------------------------------------------
C        Create brand new trial vectors for SOPPA or SOPPA(CCSD)
C        based on approximate diagonal elements of hessian.
C---------------------------------------------------------------
C
         CALL SO_TR1NEW(NNEWTR,POINT,LPOINT,ISYMTR,NEXCI,DENSIJ,
     &                  LDENSIJ,DENSAB,LDENSAB,WORK,LWORK)
C
      ELSE
C
         CALL QUIT('ERROR: Locical mistake occured in SO_TRIAL1')
C
      END IF
C
C-----------------------------------------------
C     Orthogonalize new trial vectors over S[2].
C-----------------------------------------------
C
Cekd: Same problem as for RP_ORTH_TRN, Cannot send in explicit zeros
C     as they may be assigned inside routine.
Corig  CALL SO_ORTH_TRN(0,NNEWTR,0,ISYMTR,
C     &                 DENSIJ,LDENSIJ,DENSAB,LDENSAB,WORK,LWORK)
C
      NLINDP = 0
      NOLDTR = 0
      CALL SO_ORTH_TRN('EIXVAL',NOLDTR,NNEWTR,NLINDP,ISYMTR,
     &                 DENSIJ,LDENSIJ,DENSAB,LDENSAB,WORK,LWORK)
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
     &                            '. new trial vector from SO_TRIAL1'
C
            CALL SO_READ(WORK(KTRIAL),LTR1,LUTR1E,FNTR1E,INEWTR)
            WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTRIAL+I-1),I=1,LTR1)
            CALL SO_READ(WORK(KTRIAL),LTR2,LUTR2E,FNTR2E,INEWTR)
            WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTRIAL+I-1),I=1,LTR2)
            CALL SO_READ(WORK(KTRIAL),LTR1,LUTR1D,FNTR1D,INEWTR)
            WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTRIAL+I-1),I=1,LTR1)
            CALL SO_READ(WORK(KTRIAL),LTR2,LUTR2D,FNTR2D,INEWTR)
            WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTRIAL+I-1),I=1,LTR2)
C
  300    CONTINUE
C
      END IF
C
C------------------------------------
C     Close files with trial vectors.
C------------------------------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TRIAL1')
C
      RETURN
      END
