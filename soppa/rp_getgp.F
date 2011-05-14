C
C  /* Deck rp_getgp */
      SUBROUTINE RP_GETGP(GPVEC,LGPVEC,LABEL,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Andrea Ligabue, December 2003: linear response functions 
C                                    implemented
C
C     PURPOSE: Calculate gradient property vector GPVEC for the
C              property with label LABEL and symmetry ISYMTR.
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccsdsym.h>
#include <soppinf.h>
C
      PARAMETER   (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION   GPVEC(LGPVEC)
      DIMENSION   WORK(LWORK)
      CHARACTER*8 LABEL,RTNLBL(2)
C
      LGPVCH = LGPVEC / 2
      SQ2    = DSQRT(TWO)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_GETGP')
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
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LPRP1  = N2BST(ISYMTR)
C
      KPRP1   = 1
      KEND1   = KPRP1  + LPRP1
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_GETGP',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_GETGP',' ',KEND1,LWORK)
C
C----------------------------
C     Get MO property matrix.
C----------------------------
C
      CALL SO_ONEPMO(WORK(KPRP1),LPRP1,LABEL,ISYMTR,RTNLBL,WORK(KEND1),
     &               LWORK1)
C
C----------------------------------------------------------------------
C     Initialize and calculate RPA excitation gradient property vector.
C----------------------------------------------------------------------
C
      CALL DZERO(GPVEC,LGPVEC)
C
      CALL CCRHS_J(GPVEC,ISYMTR,WORK(KPRP1))
C
      CALL DSCAL(LGPVCH,SQ2,GPVEC,1)
C
      IF (IPRSOP .GE. 6) THEN
         CALL AROUND( 'In RP_GETGP:  ' //
     &                'RPA excitation gradient property vector.' )
         CALL CC_PRP(GPVEC,GPVEC,ISYMTR,1,0)
      ENDIF
C
C----------------------------------------------------------
C     Calculate RPA de-excitation gradient property vector.
C----------------------------------------------------------
C
C
      IF (RTNLBL(2).EQ.'SYMMETRI') FACTOR = - ONE
      IF (RTNLBL(2).EQ.'ANTISYMM') FACTOR =   ONE
C
      CALL DAXPY(LGPVCH,FACTOR,GPVEC,1,GPVEC(1+LGPVCH),1)
C
      IF (IPRSOP .GE. 6) THEN
         CALL AROUND( 'In RP_GETGP: '//
     &                'RPA de-excitation gradient property vector.' )
         CALL CC_PRP(GPVEC(1+LGPVCH),GPVEC(1+LGPVCH),ISYMTR,1,0)
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_GETGP')
C
      RETURN
      END
