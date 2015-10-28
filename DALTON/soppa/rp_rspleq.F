C
C  /* Deck rp_rspleq */
      SUBROUTINE RP_RSPLEQ(ISYMTR,LABEL,FRVAL,NFRVAL,FOCKD,LFOCKD,
     &                     WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Andrea Ligabue, December 2002
C
C     PURPOSE: Drive the calculation and analysis of RPA 
C              properties
C     ISYMTR   symmetry of the property
C
C     LABEL    label for the property to be computed
C
C     FRVAL    Frequency for the properties we are to compute
C
C     NFRVAL   number of frequencies to compute for the properties
C
C     For us, NEXCI is always 1, since we are going to compute 
C     one frequency each run
C

#include "implicit.h"
#include "priunit.h"
C
      DIMENSION FRVAL(NFRVAL)
      DIMENSION FOCKD(LFOCKD)
      DIMENSION WORK(LWORK)
C
      CHARACTER*8 LABEL
C
      PARAMETER ( D100 = 100.0D0 )
      LOGICAL IMAGPROP
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
C#include "cbilrs.h"
CGET ABS_MAXITER
#include "abslrs.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_RSPLEQ')
C
C---------------------------------------------
C     Calculate the gradient property vectors.
C---------------------------------------------
C
      LGPVC1  = NT1AM(ISYMTR)
      LGPVEC  = 2 * LGPVCH
C
      KGPVEC  = 1
      KEND1   = KGPVEC + LGPVC1 !+ LGPVEC 
      LWORK1  = LWORK  - KEND1

      LDUMMY = 1
C      KDUMMY = KEND1
C
      CALL SO_MEMMAX ('RP_RSPLEQ.1',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT('RP_RSPLEQ.1',' ',KEND1,LWORK)
C
C      CALL RP_GETGP(WORK(KGPVEC),LGPVEC,LABEL,ISYMTR,WORK(KEND1),LWORK1)
      CALL SO_GETGP(WORK(KGPVEC),LGPVC1,DUMMY,LDUMMY,
     &              LABEL,ISYMTR,IMAGPROP,'AORPA',
     &              DUMMY,LDUMMY,DUMMY,LDUMMY,DUMMY,LDUMMY,
     &              DUMMY,LDUMMY,WORK(KEND1),LWORK1)
C
C---------------------------------------
C     Save property gradients on a file.
C---------------------------------------
C
      CALL SO_OPEN(LUGPVE,FNGPVE,LGPVC1)
C
      CALL SO_WRITE(WORK(KGPVEC),LGPVC1,LUGPVE,FNGPVE,1)
C
      IF ( IPRSOP .GE. 5) THEN
C
         CALL AROUND( 'In RP_RSPLEQ:  RPA gradient property vector' //
     &                ' for : '//LABEL)
C         WRITE (LUPRI,'(2A)') 'In RP_RSPLEQ:  RPA GP vector for : ',
C     &                        LABEL
         CALL OUTPUT(WORK(KGPVEC),1,LGPVC1,1,1,LGPVC1,1,1,LUPRI)
C
      ENDIF
C
      CALL SO_CLOSE(LUGPVE,FNGPVE,'KEEP')
C
C      MAXIT   = MAXCLC
      MAXIT   = ABS_MAXITER
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      DO IFREQ=1,NFRVAL
C
        IF((DABS(FRVAL(IFREQ)).LT.EPSI).AND.NSAVMX.LT.8) 
     &      NSAVMX = MIN(NT1AM(ISYMTR),8)
C
      ENDDO
C
C
CSPAS: 8/1-2004: we do not solve for the different frequencies
C                simultaneously
C      LRESINM = NFRVAL
      LRESINM = 1
CKeinSPASmehr
      LCONV   = 8
      LMXRED  = (2 * LRESINM * NSAVMX)**2
C
      KRESINM = KEND1
      KCONV   = KRESINM + LRESINM
      KREDE   = KCONV   + LCONV
      KREDS   = KREDE   + LMXRED
      KEND2   = KREDS   + LMXRED
      LWORK2  = LWORK   - KEND2
C
      CALL SO_MEMMAX ('RP_RSPLEQ.2',LWORK2)
      IF (LWORK2 .LT.0) CALL STOPIT('RP_RSPLEQ.2',' ',KEND2,LWORK)
C
C--------------------------------------------
C     Solve linear response equations for RPA
C--------------------------------------------
C
      NEXCI = 1
      CALL RP_LRSOEQ(LABEL,ISYMTR,IMAGPROP,NEXCI,MAXIT,FRVAL,NFRVAL,
     &               WORK(KRESINM),LRESINM,WORK(KCONV),LCONV,FOCKD,
     &               LFOCKD,WORK(KREDE),WORK(KREDS),LMXRED,
     &               WORK(KEND2),LWORK2)
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
      KTR1E   = KEND2
      KTR1D   = KTR1E   + LTR1E
      KRESO1E = KTR1D   + LTR1D
      KRESO1D = KRESO1E + LRESO1E
      KEND3   = KRESO1D + LRESO1D
      LWORK3  = LWORK   - KEND3
C
      CALL SO_MEMMAX ('RP_RSPLEQ.3',LWORK3)
      IF (LWORK3 .LT.0) CALL STOPIT('RP_RSPLEQ.3',' ',KEND3,LWORK)
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
            WRITE(LUPRI,'(A,I3)') ' Property symmetry',ISYMTR
C
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

cLig I think we don't need to do that
cLigC
cLigC-------------------------------------------------
cLigC        Analyse eigenvectors and write to output.
cLigC-------------------------------------------------
cLigC
cLig         THR1 = 0.1D0
cLig         THR2 = 0.1D0
cLigC
cLig         CALL RP_ANAL(WORK(KTR1E),WORK(KTR1D),LTR1E,THR1,ISYMTR)
C
C---------------------------------------
C        Write a closing line to output.
C---------------------------------------
C
            IF ( IEXCI .EQ. NEXCI ) THEN
               WRITE(LUPRI,'(/1X,A)') 
     &             '----------------------------------'//
     &             '----------------------------'
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
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('RP_RSPLEQ')
C
      RETURN
      END
