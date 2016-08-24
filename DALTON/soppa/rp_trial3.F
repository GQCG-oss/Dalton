C
C  /* Deck rp_trial3 */
      SUBROUTINE RP_TRIAL3(NNEWTR,NOLDTR,ISYMTR,IMAGPROP,IFREQ,
     &                     FRVAL,NFRVAL,
     &                     NEXCI,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Andrea Ligabue, december 2002
C
C     PURPOSE: Determine the initial trialvectors for the linear
C     response equations using eq. (19) with R = (C +-C).
C
C     NNEWTR            number of new trial vectors (OUT)
C     NOLDTR            number of old trial vectors of the last frequency
C                       ... it is not used in the first calculation
C     ISYMTR            symmetry of the trial vectors
C                       (it is the symmetry of the property)
C     IFREQ             position in the frequency array we are computing
C     FRVAL(NFRVAL)     array of frequencise
C     NEXCI             number of frequencies to be computed for that property
C                       (always 1)
C
#include "implicit.h"
#include "priunit.h"
#include "soppinf.h"
#include "ccsdsym.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, EPSI = 1.0D-8)
C
      DIMENSION FRVAL(NFRVAL)
      DIMENSION WORK(LWORK)
      LOGICAL IMAGPROP
      LOGICAL     STATIC
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_TRIAL3')
C
      STATIC = ABS(FRVAL(IFREQ)).LT.EPSI
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
      CALL SO_MEMMAX ('RP_TRIAL3.1',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT('RP_TRIAL3.1',' ',KEND1,LWORK)
C
C-----------------------------------------
C     Read diagonal E[2] and S[2] elements
C-----------------------------------------
C
      CALL GPOPEN  (LUDIAG,'SO_DIAG','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &    .FALSE.)
      REWIND LUDIAG
C
      READ(LUDIAG) ( WORK(KEDIA1+I-1), I = 1,LEDIA1)
C
      CALL GPCLOSE (LUDIAG,'KEEP')
C
C---------------------------------
C     2. allocation of work space.
C---------------------------------
C
      LGPVEC = NT1AM(ISYMTR)
      LTRIAL = NT1AM(ISYMTR)
C
      KGPVEC  = KEND1
      KTRIAL  = KGPVEC + LGPVEC
      KEND2   = KTRIAL + LTRIAL
      LWORK2  = LWORK  - KEND2
C
      CALL SO_MEMMAX ('RP_TRIAL3.2',LWORK2)
      IF (LWORK2 .LT.0) CALL STOPIT('RP_TRIAL3.2',' ',KEND2,LWORK)
CRF Factor for creating D gradient from E gradient
      DFACTOR  = -ONE
      IF (IMAGPROP) DFACTOR = ONE
C
C----------------------------------
C     Open files for trial vectors.
C----------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LEDIA1)
      CALL SO_OPEN(LUTR1D,FNTR1D,LEDIA1)
C
C-------------------------------------
C     Open file for property gradient.
C-------------------------------------
C
      CALL SO_OPEN(LUGPVE,FNGPVE,LGPVEC)
C
C---------------------------------------------------------------
C     Set the number of new trial vectors equal to NEXCI
C     remember that for us NEXCI is always 1 ...
C     I'll add one trial vector each iteration
C     Loop over the new trial vectors.
C     This loop is quite stupid since NNEWTR is always 1.
C---------------------------------------------------------------
C
      NNEWTR = NEXCI
C
      DO 200 INEWTR = 1, NNEWTR
C
C----------------------------------------------------------------------
C        For the first frequency: use the GP vector as starting trial
C        vector. I need also tu use the GP vector when freq. is 0
C        Otherwise use the last trial vector of the previous frequency.
C----------------------------------------------------------------------
C
      IF ((IFREQ .EQ. 1).OR.STATIC) THEN
C
C--------------------------------------------
C           Read the GP vector from the file.
C---------------------------------------------
C
            CALL SO_READ(WORK(KGPVEC),LGPVEC,LUGPVE,FNGPVE,1)
C
C---------------------------------------------------------------
C           Calculate the excitation part of the trial vector as
C           (E[2]diag - omega)^-1 * GP
C---------------------------------------------------------------
C
            CALL DZERO(WORK(KTRIAL),LTRIAL)
C
            DO IELEM = 1, LEDIA1
C
               TMP = ONE / ( WORK(KEDIA1+IELEM-1) - FRVAL(IFREQ) )
               WORK(KTRIAL+IELEM-1) = TMP * WORK(KGPVEC + IELEM -1)
C
            END DO
C
            CALL SO_WRITE(WORK(KTRIAL),LEDIA1,LUTR1E,FNTR1E,INEWTR)
C
C------------------------------------------------------------------
C           Calculate the de-excitation part of the trial vector as
C           (E[2]diag + omega)^-1 * GP
C-------------------------------------------------------------------
C
            CALL DZERO(WORK(KTRIAL),LTRIAL)
C
            IF(.NOT.STATIC)THEN
C
               DO IELEM = 1, LEDIA1
C
                  TMP = DFACTOR / ( WORK(KEDIA1+IELEM-1) + FRVAL(IFREQ))
                  WORK(KTRIAL+IELEM-1) = TMP
     &                                 * WORK(KGPVEC + IELEM-1)
C
               END DO
            ENDIF
C
            CALL SO_WRITE(WORK(KTRIAL),LEDIA1,LUTR1D,FNTR1D,INEWTR)
C
C-----------------------------------------------------------------------
C           For the other frequencies use the last trial vector on file.
C-----------------------------------------------------------------------
C
         ELSE
C
cspas should NOLDTR not be 1 here ?
            NOLDTR = 1
C
CRF don't do a thing
C            CALL SO_READ(WORK(KTRIAL),LTRIAL,LUTR1E,FNTR1E,NOLDTR)
C            CALL SO_WRITE(WORK(KTRIAL),LTRIAL,LUTR1E,FNTR1E,INEWTR)
C
C            CALL SO_READ(WORK(KTRIAL),LTRIAL,LUTR1D,FNTR1D,NOLDTR)
C            CALL SO_WRITE(WORK(KTRIAL),LTRIAL,LUTR1D,FNTR1D,INEWTR)
C
         ENDIF
C
  200 CONTINUE
C
C--------------------------------------
C     Close file for property gradient.
C--------------------------------------
C
      CALL SO_CLOSE(LUGPVE,FNGPVE,'KEEP')
C
      IF ( IPRSOP .GE. 9 ) THEN
C
C------------------------------------------
C        Write new trial vectors to output.
C------------------------------------------
C
         DO 300 INEWTR = 1,NNEWTR
C
            WRITE(LUPRI,'(I3,A,A)') INEWTR,
     &          '. raw trial vector in RP_TRIAL3 before',
     &          ' orthonormalization'
C
            CALL SO_READ(WORK(KTRIAL),LTRIAL,LUTR1E,FNTR1E,INEWTR)
            WRITE (LUPRI,'(I8,1X,F14.8)')
     &          (I,WORK(KTRIAL+I-1),I=1,LTRIAL)
            CALL SO_READ(WORK(KTRIAL),LTRIAL,LUTR1D,FNTR1D,INEWTR)
            WRITE (LUPRI,'(I8,1X,F14.8)')
     &          (I,WORK(KTRIAL+I-1),I=1,LTRIAL)
C
  300   CONTINUE
C
      END IF
C
C-----------------------------------------------
C     Orthogonalize new trial vectors over S[2].
C-----------------------------------------------
C
      NLINDP = 0
      NOLDTR = 0
      CALL RP_ORTH_TRN('LINEAR',NOLDTR,NNEWTR,NLINDP,ISYMTR,WORK,LWORK)
C
      IF ( IPRSOP .GE. 9 ) THEN
C
C---------------------------------------------------------
C        Write new trial othonormalized vectors to output.
C---------------------------------------------------------
C
         DO 400 INEWTR = 1,NNEWTR
C
            WRITE(LUPRI,'(I3,A,A)') INEWTR,
     &          '. raw trial vector in RP_TRIAL3 after',
     &          ' orthonormalization'
C
            CALL SO_READ(WORK(KTRIAL),LTRIAL,LUTR1E,FNTR1E,INEWTR)
            WRITE (LUPRI,'(I8,1X,F14.8)')
     &          (I,WORK(KTRIAL+I-1),I=1,LTRIAL)
            CALL SO_READ(WORK(KTRIAL),LTRIAL,LUTR1D,FNTR1D,INEWTR)
            WRITE (LUPRI,'(I8,1X,F14.8)')
     &          (I,WORK(KTRIAL+I-1),I=1,LTRIAL)
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
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('RP_TRIAL3')
C
      RETURN
      END
