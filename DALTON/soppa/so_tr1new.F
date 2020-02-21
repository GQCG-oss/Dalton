C
C  /* Deck so_tr1new */
      SUBROUTINE SO_TR1NEW(MODEL,NNEWTR,POINT,LPOINT,ISYMTR,NEXCI,
     &                     DENSIJ,
     &                     LDENSIJ,DENSAB,LDENSAB,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Determine new initial trialvectors.
C
C     'Atomic integral driven second order polarization propagator
C     calculations of the excitation spectra of napthalene and anthracene'
C     Keld L. Bak, Henrik Koch, Jens Oddershede, Ove Christiansen
C     and Stephan P. A Sauer
C     Journal of Chemical Physics, 112, 9, (2000)
C
C     Section C pp. 4176
C     SO_TR1NEW is called to create a brand new set of trial vectors for
C     SOPPA or SOPPA(CCSD) based on approximate diagonal elements of the
C     Hessian.
C
      use so_info, only : so_has_doubles, so_singles_second
C
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, STHR = 1.0D-5)
C
      CHARACTER*5 MODEL
      INTEGER   POINT(LPOINT)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION WORK(LWORK)
      LOGICAL   DOUBLES
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TR1NEW')
      DOUBLES = SO_HAS_DOUBLES(MODEL)
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LEDIA1 = NT1AM(ISYMTR)
      IF (DOUBLES) THEN
         LEDIA2 = N2P2HOP(ISYMTR)
      ELSE
         LEDIA2 = 0
      ENDIF
      LSDIA1 = NT1AM(ISYMTR)
C
      KEDIA1  = 1
      KEDIA2  = KEDIA1 + LEDIA1
      KSDIA1  = KEDIA2 + LEDIA2
      KEND1   = KSDIA1 + LSDIA1
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_TRIAL1.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_TRIAL1.1',' ',KEND1,LWORK)
C
C------------------------------------------
C     Read diagonal E[2] and S[2] elements.
C------------------------------------------
C
      CALL GPOPEN  (LUDIAG,'SO_DIAG','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
      REWIND LUDIAG
C
      READ(LUDIAG) ( WORK(KEDIA1+I-1), I = 1,LEDIA1)
      IF (DOUBLES) READ(LUDIAG) ( WORK(KEDIA2+I-1), I = 1,LEDIA2)
      READ(LUDIAG) ( WORK(KSDIA1+I-1), I = 1,LSDIA1)
C
      CALL GPCLOSE (LUDIAG,'KEEP')
C
C-------------------------------------------------------
C     Calculate approximate eigenvalues as diagonal E[2]
C     values devided by corresponding S[2] values. Lower
C     threshold for S[2] is STHR. Keep the results where
C     the E[2] diagonal used to be.
C-------------------------------------------------------
C
      DO 100 I = 1,LEDIA1
C
         KOFF1 = KEDIA1 + I - 1
         KOFF2 = KSDIA1 + I - 1
C
         IF ( ABS(WORK(KOFF2)) .GT. STHR ) THEN
            WORK(KOFF1) = WORK(KOFF1)/WORK(KOFF2)
         ELSE IF ( WORK(KOFF2) .GT. ZERO ) THEN
            WORK(KOFF1) = WORK(KOFF1)/STHR
         ELSE
            WORK(KOFF1) = WORK(KOFF1)/(-STHR)
         END IF
C
  100 CONTINUE

C
C---------------------------------
C     2. allocation of work space.
C---------------------------------
C
      LEDIA  = LEDIA1 + LEDIA2
      LPARRA = LEDIA
C
      KPARRA  = KEND1
      KEND2   = KPARRA + LPARRA
      LWORK2  = LWORK  - KEND2
C
      CALL SO_MEMMAX ('SO_TRIAL1.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('SO_TRIAL1.2',' ',KEND2,LWORK)
C
C--------------------------------------------------
C     Find the NEXCI lowest approximate eigenvalues
C     and leave pointers to them in POINT.
C--------------------------------------------------
C
      CALL SO_SORT(POINT,NEXCI,WORK(KEDIA1),LEDIA,WORK(KPARRA))
C
C---------------------------------
C     3. allocation of work space.
C---------------------------------
C
      LTRIAL = LEDIA
C
      KTRIAL  = 1
      KEND3   = KTRIAL + LTRIAL
      LWORK3  = LWORK  - KEND3
C
      CALL SO_MEMMAX ('SO_TRIAL1.3',LWORK3)
      IF (LWORK3 .LT. 0) CALL STOPIT('SO_TRIAL1.3',' ',KEND3,LWORK)
C
C---------------------------------------------------------
C     Create initial new trial vectors and write to files.
C---------------------------------------------------------
C
      CALL DZERO(WORK(KTRIAL),LTRIAL)
      DO 200 INEWTR = 1, NNEWTR
C
C
         CALL SO_WRITE(WORK(KTRIAL),LEDIA1,LUTR1D,FNTR1D,INEWTR)
         IF (DOUBLES) THEN
            CALL SO_WRITE(WORK(KTRIAL+LEDIA1),LEDIA2,
     &                                        LUTR2D,FNTR2D,INEWTR)
         END IF
C
         WORK( KTRIAL + POINT(INEWTR) - 1 ) = ONE
         CALL SO_WRITE(WORK(KTRIAL),LEDIA1,LUTR1E,FNTR1E,INEWTR)
         IF (DOUBLES) THEN
            CALL SO_WRITE(WORK(KTRIAL+LEDIA1),LEDIA2,
     &                                        LUTR2E,FNTR2E,INEWTR)
         END IF
         WORK( KTRIAL + POINT(INEWTR) - 1 ) = ZERO
C
  200 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TR1NEW')
C
      RETURN
      END
