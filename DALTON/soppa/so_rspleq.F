C
C  /* Deck so_rspleq */
      SUBROUTINE SO_RSPLEQ(MODEL,LABEL,ISYMTR,FRVAL,NFRVAL,DENSIJ,
     &                     LDENSIJ,DENSAB,LDENSAB,DENSAI,LDENSAI,T2MP,
     &                     LT2MP,FOCKD,LFOCKD,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Andrea Ligabue December 2003
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Drive the calculation and analysis of SOPPA and 
C              SOPPA(CCSD) response properties
C
C     MODEL             SOPPA or SOPPA(CCSD)
C     LABEL             property label
C     ISYMTR            property symmetry
C     FRVAL(NFRVAL)     frequencies
C
C     For us, NEXCI is always 1, since we are going to compute 
C     one frequency each run
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION FRVAL(NFRVAL)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION T2MP(LT2MP),     FOCKD(LFOCKD)
      DIMENSION WORK(LWORK)
C
      PARAMETER ( D100 = 100.0D0 )
      CHARACTER MODEL*5,LABEL*8
      LOGICAL   IMAGPROP
C
#include "codata.h"
#include "wrkrsp.h"
CSPAS:15/03-2006: merge with Dalton-2.0      
C#include "infpp.h"
CKeinSPASmehr
#include "inforb.h"
#include "infpri.h"
#include "ccsdsym.h"
#include "soppinf.h"
C cbilrs seem to be only for freq independent
C#include "cbilrs.h"
C Get ABS_MAXITER
#include "abslrs.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RSPLEQ')
C
C---------------------------------------------
C     Calculate the gradient property vectors.
C---------------------------------------------
C
      LGPVC1  =  NT1AM(ISYMTR)
CPFP   2009-Jun
      LGPVC2  =  NT2AM(ISYMTR)
C      LGPVC2  = 2 * N2P2HOP(ISYMTR)
Cend-PFP
C
      KGPVC1  = 1
      KGPVC2  = KGPVC1 + LGPVC1 
      KEND1   = KGPVC2 + LGPVC2 
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_RSPLEQ.1',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT('SO_RSPLEQ.1',' ',KEND1,LWORK)
C
      CALL SO_GETGP(WORK(KGPVC1),LGPVC1,WORK(KGPVC2),LGPVC2,
     &              LABEL,ISYMTR,IMAGPROP,MODEL,
     &              T2MP,LT2MP,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &              DENSAI,LDENSAI,WORK(KEND1),LWORK1)
C
C---------------------------------------
C     Save property gradients on a file.
C---------------------------------------
C
      CALL SO_OPEN(LUGPV1,FNGPV1,LGPVC1)
      CALL SO_OPEN(LUGPV2,FNGPV2,LGPVC2)
C
      CALL SO_WRITE(WORK(KGPVC1),LGPVC1,LUGPV1,FNGPV1,1)
      CALL SO_WRITE(WORK(KGPVC2),LGPVC2,LUGPV2,FNGPV2,1)
C
      IF(IPRSOP.GE.5) THEN
C         LGPVC1H = LGPVC1 /2
C         LGPVC2H = LGPVC2 /2
         CALL AROUND('In SO_RSPLEQ:  '//MODEL//
     &               ' gradient property vector.'//LABEL)
         IF(AOSOP) WRITE(LUPRI,'(1X,A)')
     &       'It lack the DENSAI contribution'
         CALL OUTPUT(WORK(KGPVC1),1,LGPVC1,1,1,LGPVC1,1,1,LUPRI)
         CALL OUTPUT(WORK(KGPVC2),1,LGPVC2,1,1,LGPVC2,1,1,LUPRI)
      ENDIF
C     
      CALL SO_CLOSE(LUGPV1,FNGPV1,'KEEP')
      CALL SO_CLOSE(LUGPV2,FNGPV2,'KEEP')
C      MAXIT   = MAXCLC
      MAXIT = ABS_MAXITER
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
C
CSPAS: 8/1-2004: we do not solve for the different frequencies
C                simultaneously
C      LRESINM = NEXCI
      LRESINM = 1
CKeinSPASmehr
      LCONV   = 8
      LMXRED  = (2 * LRESINM * NSAVMX)**2
C
      KRESINM = 1
      KCONV   = KRESINM + LRESINM
      KREDE   = KCONV   + LCONV
      KREDS   = KREDE   + LMXRED
      KEND    = KREDS   + LMXRED
      LWORK1  = LWORK   - KEND
C
      CALL SO_MEMMAX ('SO_RSPLEQ.2',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_RSPLEQ.2',' ',KEND1,LWORK)
C
C-----------------------------------------------
C     Solve linear response equations for SOPPA.
C-----------------------------------------------
C
      NEXCI = 1
      CALL SO_LRSOEQ(MODEL,LABEL,ISYMTR,IMAGPROP,
     &               NEXCI,MAXIT,FRVAL,NFRVAL,
     &               WORK(KRESINM),LRESINM,WORK(KCONV),LCONV,DENSIJ,
     &               LDENSIJ,DENSAB,LDENSAB,DENSAI,LDENSAI,T2MP,LT2MP,
     &               FOCKD,LFOCKD,WORK(KREDE),WORK(KREDS),LMXRED,
     &               WORK(KEND),LWORK1)
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LTR1E   = NT1AM(ISYMTR)
      LTR1D   = NT1AM(ISYMTR)
      LRESO1E = NT1AM(ISYMTR)
      LRESO1D = NT1AM(ISYMTR)
CPFP   2009-Jun
      LTR2E   = NT2AM(ISYMTR)
      LTR2D   = NT2AM(ISYMTR)
C      LTR2E   = N2P2HOP(ISYMTR)
C      LTR2D   = N2P2HOP(ISYMTR)
Cend-PFP
C      KTR1E   = KEND
C      KTR1D   = KTR1E   + LTR1E
C      KRESO1E = KTR1D   + LTR1D
c      KRESO1D = KRESO1E + LRESO1E
c      KTR2E   = KRESO1D + LRESO1D
C      KTR2D   = KTR2E   + LTR2E
C      KEND2   = KTR2D   + LTR2D
C      LWORK2  = LWORK   - KEND2
C
C      CALL SO_MEMMAX ('SO_RSPLEQ.3',LWORK2)
C      IF (LWORK2 .LT. 0) CALL STOPIT('SO_RSPLEQ.3',' ',KEND2,LWORK)
C
C----------------
C     Open files.
C----------------
C
C      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
c      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
C      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
C      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
C
C-----------------------------------------------------------
C     Calculate p-h and 2p-2h weight in excitation operator
C     and write to output together with excitation energies.
C-----------------------------------------------------------
C
C      IF ( IPRSOP .GE. 1 ) THEN
C
C         DO IEXCI = 1,NEXCI
C
C            WRITE(LUPRI,'(/1X,A)') '-------------------------------'//
C     &                             '-------------------------------'
C            WRITE(LUPRI,'(A,I3,A,I3)')
C     &           ' Excitation symmetry',ISYMTR,', state no.',IEXCI
C            IF (MODEL .EQ. 'AOSOP') THEN
C               WRITE(LUPRI,'(/A,1P,G16.8,A,3(/26X,G16.8,A),/)')
C     &              ' SOPPA excitation energy :',FRVAL(IEXCI),
C     &              ' au',
C     &              FRVAL(IEXCI)*XTEV,  ' eV',
C    &              FRVAL(IEXCI)*XTKAYS,' cm-1',
C    &              FRVAL(IEXCI)*XTKJML,' kj / mole'
C           ELSE IF (MODEL .EQ. 'AOSOC') THEN
C              WRITE(LUPRI,'(/A,1P,G16.8,A,3(/32X,G16.8,A),/)')
C    &              ' SOPPA(CCSD) excitation energy :',
C    &              FRVAL(IEXCI),' au',
C    &              FRVAL(IEXCI)*XTEV,  ' eV',
C    &              FRVAL(IEXCI)*XTKAYS,' cm-1',
C    &              FRVAL(IEXCI)*XTKJML,' kj / mole'
C           ELSE
C              WRITE(LUPRI,'(A)') 'ERROR: Unknown model in SO_RSPLEX'
C              CALL QUIT('ERROR:  Unknown model in SO_RSPLEX')
C           END IF

C           CALL SO_READ(WORK(KTR1E), LTR1E, LUTR1E,FNTR1E,IEXCI)
C           CALL SO_READ(WORK(KTR1D), LTR1D, LUTR1D,FNTR1D,IEXCI)
C
C           ISYRES = MULD2H(ISYMOP,ISYMTR)
C
C           CALL SO_RES_O(WORK(KRESO1E),LRESO1E,WORK(KRESO1D),LRESO1D,
C    &                    WORK(KTR1E),  LTR1E,  WORK(KTR1D),  LTR1D,
C    &                    DENSIJ,      LDENSIJ,DENSAB,      LDENSAB,
C    &                    ISYRES,      ISYMTR) 
C
C           W1ENM = DDOT(LTR1E,WORK(KTR1E),1,WORK(KRESO1E),1)
C           W1DNM = DDOT(LTR1D,WORK(KTR1D),1,WORK(KRESO1D),1)
C
C           CALL SO_READ(WORK(KTR2E), LTR2E, LUTR2E,FNTR2E,IEXCI)
C           CALL DCOPY(LTR2E,WORK(KTR2E),1,WORK(KTR2D),1)
C           CALL SO_TFSET(WORK(KTR2D),1,LTR2E,ISYMTR)
C           W2ENM = DDOT(LTR2E,WORK(KTR2D),1,WORK(KTR2E),1)
C
C           CALL SO_READ(WORK(KTR2D), LTR2D, LUTR2D,FNTR2D,IEXCI)
C           CALL DCOPY(LTR2E,WORK(KTR2D),1,WORK(KTR2E),1)
C           CALL SO_TFSET(WORK(KTR2E),1,LTR2D,ISYMTR)
C           W2DNM = -DDOT(LTR2D,WORK(KTR2E),1,WORK(KTR2D),1)
C
C           W1ENM = (W1ENM * D100)
C           W1DNM = (W1DNM * D100)
C           W2ENM = (W2ENM * D100)
C           W2DNM = (W2DNM * D100)
C
C           WRITE(LUPRI,'(1X,A,3(F6.2,A))')
C    &           '1p-1h + 1h-1p excitation weight: ',W1ENM,' +',W1DNM,
C    &           '  = ',W1ENM+W1DNM,' %'
C           WRITE(LUPRI,'(1X,A,3(F6.2,A),/)')
C    &           '2p-2h + 2h-2p excitation weight: ',W2ENM,' +',W2DNM,
C    &           '  = ',W2ENM+W2DNM,' %'
C
C---------------------------------------
C        Write a closing line to output.
C---------------------------------------
C
C           IF ( IEXCI .EQ. NEXCI ) THEN
C              WRITE(LUPRI,'(/1X,A)')'-------------------------------'//
C    &                               '-------------------------------'
C           END IF
C
C        END DO
C
C     END IF
C
C-----------------
C     Close files.
C-----------------
C
c      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
C      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
C      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      WRITE(LUPRI,'(A)') 'End of rspleq'
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_RSPLEQ')
C
      RETURN
      END
