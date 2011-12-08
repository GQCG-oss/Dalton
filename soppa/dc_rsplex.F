C  /* Deck dc_rsplex */
      SUBROUTINE DC_RSPLEX(ISYMTR,NEXCI,EXVAL,LEXVAL,DENSIJ,
     &                     LDENSIJ,DENSAB,LDENSAB,DENSAI,LDENSAI,T2MP,
     &                     LT2MP,FOCKD,LFOCKD,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, June 1997.
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Controls calculation of "Double corrected RPA"
C              excitation energies.
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION EXVAL(LEXVAL)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION T2MP(LT2MP),     FOCKD(LFOCKD)
      DIMENSION WORK(LWORK)
C
      PARAMETER ( ONE = 1.0D0, D100 = 100.0D0 )
C
#include <codata.h>
#include <wrkrsp.h>
CSPAS:15/03-2006: merge with Dalton-2.0      
C#include <infpp.h>
CKeinSPASmehr
#include <inforb.h>
#include <infpri.h>
#include <ccsdsym.h>
#include <soppinf.h>
#include <cbiexc.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('DC_RSPLEX')
C
C------------------
C     Write banner.
C------------------
C 
      WRITE(LUPRI,'(/,2X,A,/,15X,A,I1,/,2X,A,/)')
     &     '*****************************************************'//
     &     '*************',
     &      'Double corrected RPA calculation, Symmetry ',
     &      ISYMTR,
     &     '*****************************************************'//
     &     '*************'
C
C----------------------------------------------------------
C     Calculate "Double corrected RPA" excitation energies.
C----------------------------------------------------------
C
      CALL DC_CALC(ISYMTR,NEXCI,EXVAL,LEXVAL,DENSIJ,LDENSIJ,
     &             DENSAB,LDENSAB,DENSAI,LDENSAI,T2MP,LT2MP,FOCKD,
     &             LFOCKD,WORK,LWORK)
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LTR1E   = NT1AM(ISYMTR)
      LTR1D   = NT1AM(ISYMTR)
      LRESO1E = NT1AM(ISYMTR)
      LRESO1D = NT1AM(ISYMTR)
      LTR2E   = NT2AM(ISYMTR)
      LTR2D   = NT2AM(ISYMTR)
C
      KTR1E   = 1
      KTR1D   = KTR1E   + LTR1E
      KRESO1E = KTR1D   + LTR1D
      KRESO1D = KRESO1E + LRESO1E
      KTR2E   = KRESO1D + LRESO1D
      KTR2D   = KTR2E   + LTR2E
      KEND2   = KTR2D   + LTR2D
      LWORK2  = LWORK   - KEND2
C
      CALL SO_MEMMAX ('DC_RSPLEX',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('DC_RSPLEX',' ',KEND2,LWORK)
C
C----------------
C     Open files.
C----------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
C
C-----------------------------------------------------------
C     Calculate p-h and 2p-2h weight in excitation operator
C     and write to output together with excitation energies.
C-----------------------------------------------------------
C
      WRITE(LUPRI,9010)
      WRITE(LUPRI,9011)
      WRITE(LUPRI,9010)
      DO 100 IEXCI = 1,NEXCI
         WRITE(LUPRI,9012) IEXCI,EXVAL(IEXCI)
 100  CONTINUE
      WRITE(LUPRI,9010)
C
      IF ( IPRSOP .GE. 1 ) THEN
C
         DO 200 IEXCI = 1,NEXCI
C     
            WRITE(LUPRI,'(/1X,A)') '-------------------------------'//
     &                             '-------------------------------'
            WRITE(LUPRI,'(A,I3,A,I3,//A,1P,G16.8,A,3(/26X,G16.8,A),/)')
     &           ' Excitation symmetry',ISYMTR,', state no.',IEXCI,
     &           ' DCRPA excitation energy :',EXVAL(IEXCI),
     &           ' au',
     &           EXVAL(IEXCI)*XTEV,  ' eV',
     &           EXVAL(IEXCI)*XTKAYS,' cm-1',
     &           EXVAL(IEXCI)*XTKJML,' kj / mole'

            CALL SO_READ(WORK(KTR1E), LTR1E, LUTR1E,FNTR1E,IEXCI)
            CALL SO_READ(WORK(KTR1D), LTR1D, LUTR1D,FNTR1D,IEXCI)
C
            ISYRES = MULD2H(ISYMOP,ISYMTR)
C
cKeld         CALL RP_RES_O(WORK(KRESO1E),LRESO1E,WORK(KRESO1D),LRESO1D,
cKeld     &                 WORK(KTR1E),  LTR1E,  WORK(KTR1D),  LTR1D,
cKeld     &                 ISYRES,      ISYMTR)
            CALL SO_RES_O(WORK(KRESO1E),LRESO1E,WORK(KRESO1D),LRESO1D,
     &                    WORK(KTR1E),  LTR1E,  WORK(KTR1D),  LTR1D,
     &                    DENSIJ,      LDENSIJ,DENSAB,      LDENSAB,
     &                    ISYRES,      ISYMTR) 
C
            W1ENM = DDOT(LTR1E,WORK(KTR1E),1,WORK(KRESO1E),1)
            W1DNM = DDOT(LTR1D,WORK(KTR1D),1,WORK(KRESO1D),1)
C
            CALL SO_READ(WORK(KTR2E), LTR2E, LUTR2E,FNTR2E,IEXCI)
            CALL DCOPY(LTR2E,WORK(KTR2E),1,WORK(KTR2D),1)
            CALL SO_TFSET(WORK(KTR2D),1,LTR2E,ISYMTR)
            W2ENM = DDOT(LTR2E,WORK(KTR2D),1,WORK(KTR2E),1)
C
            CALL SO_READ(WORK(KTR2D), LTR2D, LUTR2D,FNTR2D,IEXCI)
            CALL DCOPY(LTR2E,WORK(KTR2D),1,WORK(KTR2E),1)
            CALL SO_TFSET(WORK(KTR2E),1,LTR2D,ISYMTR)
            W2DNM = -DDOT(LTR2D,WORK(KTR2E),1,WORK(KTR2D),1)
C
C--------------------------------------
C        Normalize weight coeffictions.
C--------------------------------------
C
            SWNM   = ONE / (W1ENM + W1DNM + W2ENM + W2DNM)
C
            W1ENM = W1ENM * SWNM
            W1DNM = W1DNM * SWNM
            W2ENM = W2ENM * SWNM
            W2DNM = W2DNM * SWNM
C
            W1ENM = (W1ENM * D100)
            W1DNM = (W1DNM * D100)
            W2ENM = (W2ENM * D100)
            W2DNM = (W2DNM * D100)
C
            WRITE(LUPRI,'(1X,A,3(F6.2,A))')
     &           '1p-1h + 1h-1p excitation weight: ',W1ENM,' +',W1DNM,
     &           '  = ',W1ENM+W1DNM,' %'
            WRITE(LUPRI,'(1X,A,3(F6.2,A),/)')
     &           '2p-2h + 2h-2p excitation weight: ',W2ENM,' +',W2DNM,
     &           '  = ',W2ENM+W2DNM,' %'
C
C------------------------------------------------------------
C        Normalize RPA(D) eigenvector defined as R(0) + R(1).
C------------------------------------------------------------
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
C-------------------------------------------------
C        Analyse eigenvectors and write to output.
C-------------------------------------------------
C
            THR1 = 0.1D0
            THR2 = 0.1D0
C
            CALL SO_ANAL(WORK(KTR1E),WORK(KTR1D),LTR1E,WORK(KTR2E),
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
 200     CONTINUE
C
      END IF
C
C-----------------
C     Close files.
C-----------------
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
      CALL QEXIT('DC_RSPLEX')
C
      RETURN
C
 9010 FORMAT(25X,'------------------------')
 9011 FORMAT(25X,'Excitation   Energy (au)')
 9012 FORMAT(26X,I5,2X,F15.8)
C
      END
