#include <single.h>
C
C  /* Deck rp_rsplex */
      SUBROUTINE RP_RSPLEX(ISYMTR,NEXCI,EXVAL,LEXVAL,FOCKD,LFOCKD,
     &                     WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, June 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Drive the calculation and analysis of RPA excitation 
C              energies and vectors.
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION EXVAL(LEXVAL)
      DIMENSION FOCKD(LFOCKD)
      DIMENSION WORK(LWORK)
C
      PARAMETER ( D100 = 100.0D0 )
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
      CALL QENTER('RP_RSPLEX')
C
C---------------------------------------------
C     Solve linear response equations for RPA.
C---------------------------------------------
C
      MAXIT   = MAXITE
C
      LRESINM = NEXCI
      LCONV   = 8
      LMXRED  = (2 * NEXCI * NSAVMX)**2
C
      KRESINM = 1
      KCONV   = KRESINM + LRESINM
      KREDE   = KCONV   + LCONV
      KREDS   = KREDE   + LMXRED
      KEND    = KREDS   + LMXRED
      LWORK1  = LWORK   - KEND
C
      CALL SO_MEMMAX ('RP_RSPLEX',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT('RP_RSPLEX',' ',KEND,LWORK)
C
      CALL RP_LRSOLV(ISYMTR,NEXCI,MAXIT,EXVAL,LEXVAL,WORK(KRESINM),
     &               LRESINM,WORK(KCONV),LCONV,FOCKD,LFOCKD,
     &               WORK(KREDE),WORK(KREDS),LMXRED,WORK(KEND),LWORK1)
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LTR1E   = NT1AM(ISYMTR)
      LTR1D   = NT1AM(ISYMTR)
      LRESO1E = NT1AM(ISYMTR)
      LRESO1D = NT1AM(ISYMTR)
C
      KTR1E   = KEND
      KTR1D   = KTR1E   + LTR1E
      KRESO1E = KTR1D   + LTR1D
      KRESO1D = KRESO1E + LRESO1E
      KEND2   = KRESO1D + LRESO1D
      LWORK2  = LWORK   - KEND2
C
      CALL SO_MEMMAX ('RP_RSPLEX',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('RP_RSPLEX',' ',KEND2,LWORK)
C
C----------------
C     Open files.
C----------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
C
C-----------------------------------------------------------
C     Calculate p-h and h-p weight in excitation operator
C     and write to output together with excitation energies.
C-----------------------------------------------------------
C
      IF ( IPRSOP .GE. 1 ) THEN
C
         DO IEXCI = 1,NEXCI
C
            WRITE(LUPRI,'(/1X,A)') '-------------------------------'//
     &                             '-------------------------------'
            WRITE(LUPRI,'(A,I3,A,I3,//A,1P,G16.8,A,3(/24X,G16.8,A),/)')
     &           ' Excitation symmetry',ISYMTR,', state no.',IEXCI,
     &           ' RPA excitation energy :',EXVAL(IEXCI),
     &           ' au',
     &           EXVAL(IEXCI)*XTEV,  ' eV',
     &           EXVAL(IEXCI)*XTKAYS,' cm-1',
     &           EXVAL(IEXCI)*XTKJML,' kj / mole'

            CALL SO_READ(WORK(KTR1E), LTR1E, LUTR1E,FNTR1E,IEXCI)
            CALL SO_READ(WORK(KTR1D), LTR1D, LUTR1D,FNTR1D,IEXCI)
C
            ISYRES = MULD2H(ISYMOP,ISYMTR)
C
            CALL RP_RES_O(WORK(KRESO1E),LRESO1E,WORK(KRESO1D),LRESO1D,
     &                    WORK(KTR1E),  LTR1E,  WORK(KTR1D),  LTR1D,
     &                    ISYRES,      ISYMTR) 
C
            W1ENM = DDOT(LTR1E,WORK(KTR1E),1,WORK(KRESO1E),1)
            W1DNM = DDOT(LTR1D,WORK(KTR1D),1,WORK(KRESO1D),1)
C
            W1ENM = (W1ENM * D100)
            W1DNM = (W1DNM * D100)
C
            WRITE(LUPRI,'(1X,A,3(F6.2,A),/)')
     &           '1p-1h + 1h-1p excitation weight: ',W1ENM,' +',W1DNM,
     &           '  = ',W1ENM+W1DNM,' %'
C
C-------------------------------------------------
C        Analyse eigenvectors and write to output.
C-------------------------------------------------
C
            THR1 = 0.1D0
            THR2 = 0.1D0
C
            CALL RP_ANAL(WORK(KTR1E),WORK(KTR1D),LTR1E,THR1,ISYMTR)
C
C---------------------------------------
C        Write a closing line to output.
C---------------------------------------
C
            IF ( IEXCI .EQ. NEXCI ) THEN
               WRITE(LUPRI,'(/1X,A)') 
     &              '----------------------------------'//
     &              '----------------------------'
            END IF
C
         END DO
C
      END IF
C
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_RSPLEX')
C
      RETURN
      END
