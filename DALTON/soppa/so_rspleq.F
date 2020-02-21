C
C  /* Deck so_rspleq */
      SUBROUTINE SO_RSPLEQ(MODEL,LABEL,IMAGPROP,ISYMTR,
     &                     FRVAL,NFRVAL,DENSIJ,
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
      use so_info, only: so_has_doubles
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
      LOGICAL, INTENT(OUT) :: IMAGPROP
      LOGICAL   DOUBLES
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
      DOUBLES = SO_HAS_DOUBLES(MODEL)
C
C---------------------------------------------
C     Calculate the gradient property vectors.
C---------------------------------------------
C
      LGPVC1  =  NT1AM(ISYMTR)
CPFP   2009-Jun
      IF (DOUBLES) THEN
         LGPVC2  =  N2P2HOP(ISYMTR)
      ELSE
         LGPVC2  = 0
      ENDIF
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
      CALL SO_WRITE(WORK(KGPVC1),LGPVC1,LUGPV1,FNGPV1,1)
C
      IF (DOUBLES) THEN
         CALL SO_OPEN(LUGPV2,FNGPV2,LGPVC2)
         CALL SO_WRITE(WORK(KGPVC2),LGPVC2,LUGPV2,FNGPV2,1)
      ENDIF
C
      IF(IPRSOP.GE.5) THEN
C         LGPVC1H = LGPVC1 /2
C         LGPVC2H = LGPVC2 /2
         CALL AROUND('In SO_RSPLEQ:  '//MODEL//
     &               ' gradient property vector.'//LABEL)
         IF(MODEL.EQ.'AOSOP') WRITE(LUPRI,'(1X,A)')
     &       'It lack the DENSAI contribution'
         CALL OUTPUT(WORK(KGPVC1),1,LGPVC1,1,1,LGPVC1,1,1,LUPRI)
         IF(DOUBLES) CALL OUTPUT(WORK(KGPVC2),1,LGPVC2,
     &                           1,1,LGPVC2,1,1,LUPRI)
      ENDIF
C
      CALL SO_CLOSE(LUGPV1,FNGPV1,'KEEP')
      IF (DOUBLES)  CALL SO_CLOSE(LUGPV2,FNGPV2,'KEEP')
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
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_RSPLEQ')
C
      RETURN
      END
