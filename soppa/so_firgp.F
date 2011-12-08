C
C  /* Deck so_firgp */
      SUBROUTINE SO_FIRGP(GPVC2,LGPVC2,T2AM,LT2AM,PRP1,LPRP1,PR1IJ,
     &                    LPR1IJ,PR1AB,LPR1AB,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate 2p-2h part of gradient property vectors.
C              There is only a first order contribution.
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccsdsym.h>
#include <ccorb.h>
#include <soppinf.h>
C
      PARAMETER   (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION   GPVC2(LGPVC2), T2AM(LT2AM)
      DIMENSION   PRP1(LPRP1),   PR1IJ(LPR1IJ), PR1AB(LPR1AB)
      DIMENSION   WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_FIRGP')
C
C---------------------------------------------------
C     Repack the one-particle property mo-integrals.
C---------------------------------------------------
C
      CALL SO_RPPRP1(PRP1,LPRP1,PR1IJ,LPR1IJ,PR1AB,LPR1AB,ISYMTR)
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LT2SQ  = NT2SQ(1)
C
      KT2SQ  = 1
      KEND1  = KT2SQ + LT2SQ
      LWORK1 = LWORK - KEND1
C
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_FIRGP',' ',KEND1,LWORK)
C
C---------------------------------
C     Square up the t2-amplitudes.
C---------------------------------
C
      CALL CC_T2SQ(T2AM,WORK(KT2SQ),1)
C
C-----------------------------------------------------------------------
C     Contract T2-amplitudes with one-particle property integrals and
C     scale with a half to calculate the 2p-2h gradient property vector.
C-----------------------------------------------------------------------
C
      CALL CCRHS_E(GPVC2,WORK(KT2SQ),PR1AB,PR1IJ,WORK(KEND1),LWORK1,
     &             1,ISYMTR)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_FIRGP')
C
      RETURN
      END
