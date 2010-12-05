C  /* Deck so_onepmo */
      SUBROUTINE SO_ONEPMO(PRP1,LPRP1,LABEL,KSYMOP,RTNLBL,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1998
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Andrea Ligabue, December 2003: linear response functions 
C                                    implemented
C
C     PURPOSE: Calculate MO-property matrix for a one-electron operator
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION PRP1(LPRP1)
      DIMENSION WORK(LWORK)
C
#include <ccsdsym.h>
#include <ccorb.h>
#include <soppinf.h>
C
      CHARACTER*8 LABEL,RTNLBL(2)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_ONEPMO')
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LCMO   = NLAMDT
C
      KCMO    = 1
      KEND1   = KCMO   + LCMO
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_ONEPMO',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_ONEPMO',' ',KEND1,LWORK)
C
C-------------------------
C     Get MO coefficients.
C-------------------------
C
      CALL SO_GETMO(WORK(KCMO),LCMO,WORK(KEND1),LWORK1,KEND1)
C
C--------------------------------------------
C     Get one electron AO property integrals.
C--------------------------------------------
C
      CALL DZERO(PRP1,LPRP1)
      CALL SO_ONEP(PRP1,LPRP1,LABEL,KSYMOP,RTNLBL,WORK(KEND1),LWORK1)
C
C--------------------------------------------------------------
C     Transform AO property integrals to MO property integrals.
C--------------------------------------------------------------
C
      CALL CC_FCKMO(PRP1,WORK(KCMO),WORK(KCMO),
     &              WORK(KEND1),LWORK1,KSYMOP,1,1)
C
      IF (IPRSOP .GE. 6) THEN
         CALL AROUND( 'In SO_ONEPMO: MO property matrix' )
         CALL CC_PRFCKMO(PRP1,KSYMOP)
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_ONEPMO')
C
      RETURN
      END
