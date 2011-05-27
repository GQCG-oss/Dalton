C
C  /* Deck so_getgp */
      SUBROUTINE SO_GETGP(GPVC1,LGPVC1,GPVC2,LGPVC2,LABEL,ISYMTR,
     &                    T2AM,LT2AM,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                    DENSAI,LDENSAI,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Andrea Ligabue, December 2003: linear response functions 
C                                    implemented
C
C     PURPOSE: Calculate gradient property 1p-1h and 2p-2h vectors
C              (GPVC1 and GPVC2) for the property with label LABEL
C              and symmetry ISYMTR.
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <soppinf.h>
C
      PARAMETER   (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION   GPVC1(LGPVC1), GPVC2(LGPVC2), T2AM(LT2AM)
      DIMENSION   DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION   WORK(LWORK)
      CHARACTER*8 LABEL,RTNLBL(2)
C
      LGPVC1H = LGPVC1 / 2
      LGPVC2H = LGPVC2 / 2
      SQ2 = DSQRT(TWO)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_GETGP')
C
      IF(IPRSOP.GT.10) THEN
C
         CALL AROUND('DENSIJ')
         CALL OUTPUT(DENSIJ,1,LDENSIJ,1,1,LDENSIJ,1,1,LUPRI)
         CALL AROUND('DENSAB')
         CALL OUTPUT(DENSAB,1,LDENSAB,1,1,LDENSAB,1,1,LUPRI)
         CALL AROUND('DENSAI')
         CALL OUTPUT(DENSAI,1,LDENSAI,1,1,LDENSAI,1,1,LUPRI)
         CALL AROUND('T2AM')
         CALL OUTPUT(T2AM,1,LT2AM,1,1,LT2AM,1,1,LUPRI)
C
      ENDIF
C
C----------------------------------------------------------
C     Stop if requested properties are not implemented yet.
C----------------------------------------------------------
C
      IF ( (LABEL(3:8).EQ.'SPNORB') .OR. (LABEL(5:8).EQ.'LAGR') .OR.
     &     (LABEL(2:7).EQ.'LONMAG') ) THEN
         WRITE(LUPRI,*) LABEL,' is not implemented yet.'
         CALL QUIT('ERROR, Program have to stop')
      END IF
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LPRP1  = N2BST(ISYMTR)
C
      KPRP1   = 1
      KEND1   = KPRP1 + LPRP1
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_GETGP.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_GETGP.1',' ',KEND1,LWORK)
C
C----------------------------
C     Get MO property matrix.
C----------------------------
C
      CALL SO_ONEPMO(WORK(KPRP1),LPRP1,LABEL,ISYMTR,RTNLBL,WORK(KEND1),
     &               LWORK1)
C
C-----------------------------------------------------------------
C     Initialize and determine the SOPPA gradient property vector.
C-----------------------------------------------------------------
C
      CALL DZERO(GPVC1,LGPVC1)
      CALL DZERO(GPVC2,LGPVC2)
C
C-------------------------------------------------------------
C     Calculate zeroth order contribution to 1p-1h part of the
C     gradient property vector.
C-------------------------------------------------------------
C
      CALL CCRHS_J(GPVC1,ISYMTR,WORK(KPRP1))
C
      CALL DSCAL(LGPVC1H,SQ2,GPVC1,1)
C
      IF (IPRSOP .GE. 6) THEN
         CALL AROUND( 'In SO_GETGP:  '//
     &                'RPA excitation gradient property vector.' )
         CALL CC_PRP(GPVEC,GPVEC,ISYMTR,1,0)
      ENDIF
C
C-------------------------------------------------------------
C     Calculate second order contribution to 1p-1h part of the
C     gradient property vector.
C-------------------------------------------------------------
C
      CALL SO_SECGP(GPVC1,LGPVC1H,WORK(KPRP1),LPRP1,ISYMTR,DENSIJ,
     &              LDENSIJ,DENSAB,LDENSAB,DENSAI,LDENSAI)
C
C---------------------------------
C     2. allocation of work space.
C---------------------------------
C
      LPR1IJ = NIJDEN(1)
      LPR1AB = NABDEN(1)
C
      KPR1IJ  = KEND1
      KPR1AB  = KPR1IJ + LPR1IJ
      KEND2   = KPR1AB + LPR1AB
      LWORK2  = LWORK  - KEND2
C
      CALL SO_MEMMAX ('SO_GETGP.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('SO_GETGP.2',' ',KEND2,LWORK)
C
C--------------------------------------------------------------
C     Calculate the 2p-2p part of the gradient property vector.
C--------------------------------------------------------------
C
      CALL SO_FIRGP(GPVC2,LGPVC2H,T2AM,LT2AM,WORK(KPRP1),LPRP1,
     &              WORK(KPR1IJ),LPR1IJ,WORK(KPR1AB),LPR1AB,ISYMTR,
     &              WORK(KEND2),LWORK2)
C
C----------------------------------------------
C     Print gradient property vector to output.
C----------------------------------------------
C
      IF (IPRSOP .GE. 6 ) THEN
         CALL AROUND('SO_GETGP: SOPPA '//
     &       'excitation 1p1h and 2p2h gradient property vector')
         CALL CC_PRP(GPVC1,GPVC2,ISYMTR,1,1)
      ENDIF
C
C------------------------------------------------------------
C     Calculate SOPPA de-excitation gradient property vector.
C------------------------------------------------------------
C
      IF (RTNLBL(2).EQ.'SYMMETRI') FACTOR = - ONE
      IF (RTNLBL(2).EQ.'ANTISYMM') FACTOR =   ONE
C
      CALL DAXPY(LGPVC1H,FACTOR,GPVC1,1,GPVC1(1+LGPVC1H),1)
      CALL DAXPY(LGPVC2H,FACTOR,GPVC2,1,GPVC2(1+LGPVC2H),1)
C
      IF (IPRSOP .GE. 6) THEN
         CALL AROUND('SO_GETGP: SOPPA '//
     &       'de-excitation 1p1h and 2p2h gradient property vector')
         CALL CC_PRP(GPVC1(1+LGPVC1H),GPVC2(1+LGPVC2H),ISYMTR,1,1)
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_GETGP')
C
      RETURN
      END
