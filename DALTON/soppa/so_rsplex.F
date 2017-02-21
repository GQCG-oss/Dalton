C  /* Deck so_rsplex */
      SUBROUTINE SO_RSPLEX(MODEL,ISYMTR,NEXCI,EXVAL,LEXVAL,DENSIJ,
     &                     LDENSIJ,DENSAB,LDENSAB,DENSAI,LDENSAI,T2MP,
     &                     LT2MP,FOCKD,LFOCKD,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, June 1997
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Drive the calculation and analysis of SOPPA and
C              SOPPA(CCSD) excitation energies and vectors.
C
      use so_info, only: so_full_name, so_has_doubles,
     &                   so_double_correction
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION EXVAL(LEXVAL)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION T2MP(LT2MP),     FOCKD(LFOCKD)
      DIMENSION WORK(LWORK)
C
      PARAMETER ( ONE = 1.0D0, D100 = 100.0D0 )
      CHARACTER MODEL*5
      CHARACTER*8 PDENS_LABEL
C
      LOGICAL   DOUBLES, DOUB_CORR
C
#include "codata.h"
#include "wrkrsp.h"
#include "inforb.h"
#include "infpri.h"
#include "ccsdsym.h"
#include "soppinf.h"
#include "cbiexc.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RSPLEX')
C
      DOUBLES = SO_HAS_DOUBLES(MODEL)
C
      DOUB_CORR = SO_DOUBLE_CORRECTION(MODEL)
C
C----------------------------------------------------------
C     Calculate "Double corrected (H)RPA" excitation energies.
C----------------------------------------------------------
C
      IF (DOUB_CORR) THEN
C
         KEND = 1
         LWORK1 = LWORK - KEND
C
         CALL SO_MEMMAX ('SO_RSPLEX.1',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT('SO_RSPLEX.1',' ',KEND1,LWORK)
C
         CALL DC_CALC(MODEL,ISYMTR,NEXCI,EXVAL,LEXVAL,
     &                DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                DENSAI,LDENSAI,
     &                T2MP,LT2MP,FOCKD,LFOCKD,
     &                WORK(KEND),LWORK1)
C
C-----------------------------------------------
C     Solve linear response equation
C-----------------------------------------------
C
      ELSE
C
         MAXIT   = MAXITE
C
         LRESINM = NEXCI
CPi 20.06.16 Old version had this limitation on converged results
C         LCONV   = 8
         LCONV   = NEXCI
C
C         LMXRED  = (2 * NEXCI * NSAVMX)**2
CPi Adjust max size of reduced space if necessary
         IF (NEXCI*NSAVMX .GE. (NT1AM(ISYMTR)+N2P2HOP(ISYMTR))) THEN
            LMXRED = ( 2 * (NT1AM(ISYMTR)+N2P2HOP(ISYMTR)) )**2
            WRITE(LUPRI,*) 'SO_RSPLEX: Adjusting LMXRED'
         ELSE
            LMXRED  = (2 * NEXCI * NSAVMX)**2
         END IF
Cend-Pi
C
         KRESINM = 1
         KCONV   = KRESINM + LRESINM
         KREDE   = KCONV   + LCONV
         KREDS   = KREDE   + LMXRED
         KEND    = KREDS   + LMXRED
         LWORK1  = LWORK   - KEND
C
         CALL SO_MEMMAX ('SO_RSPLEX.1',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT('SO_RSPLEX.1',' ',KEND,LWORK)
C
         CALL SO_LRSOLV(MODEL,ISYMTR,NEXCI,MAXIT,EXVAL,LEXVAL,
     &                  WORK(KRESINM),LRESINM,WORK(KCONV),LCONV,DENSIJ,
     &                  LDENSIJ,DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                  T2MP,LT2MP,
     &                  FOCKD,LFOCKD,WORK(KREDE),WORK(KREDS),LMXRED,
     &                  WORK(KEND),LWORK1)
C
      END IF
C
C-----------------------------------------------------------
C     Calculate p-h and 2p-2h weight in excitation operator
C     and write to output together with excitation energies.
C-----------------------------------------------------------
C
      IF ( IPRSOP .GE. 1 .AND. NEXCI .GT. 0) THEN
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
         LTR1E   = NT1AM(ISYMTR)
         LTR1D   = NT1AM(ISYMTR)
         LRESO1E = NT1AM(ISYMTR)
         LRESO1D = NT1AM(ISYMTR)
         IF(DOUBLES) THEN
            LTR2E   = N2P2HOP(ISYMTR)
            LTR2D   = N2P2HOP(ISYMTR)
         ELSE
            LTR2E   = 0
            LTR2D   = 0
         ENDIF
C
         KTR1E   = KEND
         KTR1D   = KTR1E   + LTR1E
         KRESO1E = KTR1D   + LTR1D
         KRESO1D = KRESO1E + LRESO1E
         KTR2E   = KRESO1D + LRESO1D
         KTR2D   = KTR2E   + LTR2E
         KEND2   = KTR2D   + LTR2D
         LWORK2  = LWORK   - KEND2
C
         CALL SO_MEMMAX ('SO_RSPLEX.2',LWORK2)
         IF (LWORK2 .LT. 0) CALL STOPIT('SO_RSPLEX.2',' ',KEND2,LWORK)
C
C----------------
C     Open files.
C----------------
C
         CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
         CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
         IF ( DOUBLES ) THEN
            CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
            CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
         ENDIF
C
C-----------------------------------------------------------
C        Write excitation energies to output
C-----------------------------------------------------------
C
         DO IEXCI = 1,NEXCI
C
            WRITE(LUPRI,'(/1X,A)') '-------------------------------'//
     &                             '-------------------------------'
            IF (TRIPLET) THEN
               WRITE(LUPRI,'(/17X,A)') 'Triplet Excitation Energies'
            ELSE
               WRITE(LUPRI,'(/17X,A)') 'Singlet Excitation Energies'
            END IF
            WRITE(LUPRI,'(/1X,A)') '-------------------------------'//
     &                             '-------------------------------'
            WRITE(LUPRI,'(A,I3,A,I3)')
     &           ' Excitation symmetry',ISYMTR,', state no.',IEXCI
            WRITE(LUPRI,'(/2A,1P,G16.8,A,3(/32X,G16.8,A),/)')
     &              ADJUSTR(SO_FULL_NAME(MODEL)),
     &              ' excitation energy :',
     &              EXVAL(IEXCI),' au',
     &              EXVAL(IEXCI)*XTEV,  ' eV',
     &              EXVAL(IEXCI)*XTKAYS,' cm-1',
     &              EXVAL(IEXCI)*XTKJML,' kj / mole'
C
C-----------------------------------------------------------
C     Calculate p-h and h-p weight in excitation operator
C     and write to output.
C-----------------------------------------------------------
C
            CALL SO_READ(WORK(KTR1E), LTR1E, LUTR1E,FNTR1E,IEXCI)
            CALL SO_READ(WORK(KTR1D), LTR1D, LUTR1D,FNTR1D,IEXCI)
C
            ISYRES = MULD2H(ISYMOP,ISYMTR)
C
            IF (MODEL.EQ.'AORPA') THEN
               DO I = 0, LTR1E-1
                  WORK(KRESO1E+I) = WORK(KTR1E)
               END DO
               DO I = 0, LTR1E-1
                  WORK(KRESO1D+I) = -WORK(KTR1D)
               END DO
            ELSE
               CALL SO_RES_O(WORK(KRESO1E),LRESO1E,
     &                       WORK(KRESO1D),LRESO1D,
     &                       WORK(KTR1E),  LTR1E,  WORK(KTR1D),  LTR1D,
     &                       DENSIJ,      LDENSIJ,DENSAB,      LDENSAB,
     &                       ISYRES,      ISYMTR)
            END IF
C
            W1ENM = DDOT(LTR1E,WORK(KTR1E),1,WORK(KRESO1E),1)
            W1DNM = DDOT(LTR1D,WORK(KTR1D),1,WORK(KRESO1D),1)
C
            IF (DOUBLES) THEN
               CALL SO_READ(WORK(KTR2E), LTR2E, LUTR2E,FNTR2E,IEXCI)
               CALL DCOPY(LTR2E,WORK(KTR2E),1,WORK(KTR2D),1)
C
               IF (.NOT. TRIPLET) THEN
                  CALL SO_TFSET(WORK(KTR2D),1,LTR2E,ISYMTR)
               END IF
               W2ENM = DDOT(LTR2E,WORK(KTR2D),1,WORK(KTR2E),1)
C
               CALL SO_READ(WORK(KTR2D), LTR2D, LUTR2D,FNTR2D,IEXCI)
               CALL DCOPY(LTR2E,WORK(KTR2D),1,WORK(KTR2E),1)
               IF (.NOT. TRIPLET) THEN
                  CALL SO_TFSET(WORK(KTR2E),1,LTR2D,ISYMTR)
               END IF
C
               W2DNM = -DDOT(LTR2D,WORK(KTR2E),1,WORK(KTR2D),1)
            ELSE
               W2ENM = 0.0D0
               W2DNM = 0.0D0
            ENDIF

C
C----------------------------------------------------------
C        Normalize weight coeffictions and eigenvectors if
C        RPA(D) or HRPA(D)
C----------------------------------------------------------
C
            IF (DOUB_CORR) THEN
C
               SWNM   = ONE / (W1ENM + W1DNM + W2ENM + W2DNM)
C
               W1ENM = W1ENM * SWNM
               W1DNM = W1DNM * SWNM
               W2ENM = W2ENM * SWNM
               W2DNM = W2DNM * SWNM
C
               SQSWNM = DSQRT(SWNM)
C
               CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,IEXCI)
               CALL DSCAL(LTR1E,SQSWNM,WORK(KTR1E),1)
               CALL SO_WRITE(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,IEXCI)
C
               CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,IEXCI)
               CALL DSCAL(LTR1D,SQSWNM,WORK(KTR1D),1)
               CALL SO_WRITE(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,IEXCI)
C
               CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,IEXCI)
               CALL DSCAL(LTR2E,SQSWNM,WORK(KTR2E),1)
               CALL SO_WRITE(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,IEXCI)
C
               CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,IEXCI)
               CALL DSCAL(LTR2D,SQSWNM,WORK(KTR2D),1)
               CALL SO_WRITE(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,IEXCI)
C
            END IF
C
C----------------------------------------------------------
C           Write excitation weights to output
C----------------------------------------------------------
C
            W1ENM = (W1ENM * D100)
            W1DNM = (W1DNM * D100)
            W2ENM = (W2ENM * D100)
            W2DNM = (W2DNM * D100)
C
            WRITE(LUPRI,'(1X,A,3(F6.2,A))')
     &           '1p-1h + 1h-1p excitation weight: ',W1ENM,' +',W1DNM,
     &           '  = ',W1ENM+W1DNM,' %'
            IF (DOUBLES) THEN
               WRITE(LUPRI,'(1X,A,3(F6.2,A),/)')
     &           '2p-2h + 2h-2p excitation weight: ',W2ENM,' +',W2DNM,
     &           '  = ',W2ENM+W2DNM,' %'
            ENDIF
C
C-------------------------------------------------
C        Analyse eigenvectors and write to output.
C-------------------------------------------------
C
            THR1 = 0.1D0
            THR2 = 0.1D0
CPi TR2D is already in WORK(KTR2D)
            IF (DOUBLES) CALL SO_READ(WORK(KTR2E), LTR2E,
     &                               LUTR2E,FNTR2E,IEXCI)
C
            CALL SO_ANAL(DOUBLES,WORK(KTR1E),WORK(KTR1D),
     &                   LTR1E,WORK(KTR2E),
     &                   WORK(KTR2D),LTR2E,THR1,THR2,ISYMTR)
C
C---------------------------------------
C        Write a closing line to output.
C---------------------------------------
C
            IF ( IEXCI .EQ. NEXCI ) THEN
               WRITE(LUPRI,'(/1X,A)')'-------------------------------'//
     &                               '-------------------------------'
            END IF
C
         END DO
C
C-----------------
C     Close files.
C-----------------
C
         CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
         CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
         IF (DOUBLES) THEN
            CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
            CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
         ENDIF
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RSPLEX')
C
      RETURN
      END
