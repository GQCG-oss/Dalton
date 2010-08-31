C  /* Deck so_tr1rp */
      SUBROUTINE SO_TR1RP(NNEWTR,LTR2,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Use RPA trial vectors as first guess for SOPPA
C              or SOPPA(CCSD). Since the 1p1h vectors are in place,
C              this routine we just zero the 2p2h vectors.
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
C
      DIMENSION WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TR1RP')
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      KTR2    = 1
      KEND    = KTR2 + LTR2
      LWORK1  = LWORK - KEND
C
      CALL SO_MEMMAX ('SO_TR1RP',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_TR1RP',' ',KEND,LWORK)
C
      CALL DZERO(WORK(KTR2),LTR2)
C
      DO 100 INEWTR = 1, NNEWTR
C
         CALL SO_WRITE(WORK(KTR2),LTR2,LUTR2E,FNTR2E,INEWTR)
         CALL SO_WRITE(WORK(KTR2),LTR2,LUTR2D,FNTR2D,INEWTR)
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TR1RP')
C
      RETURN
      END
