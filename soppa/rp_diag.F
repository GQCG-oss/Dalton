C
C  /* Deck rp_diag */
      SUBROUTINE RP_DIAG(FOCKD,LFOCKD,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate diagonale parts of the RPA E[2] and S[2]
C              matrices.
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccsdsym.h>
#include <soppinf.h>
#include <ccorb.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION FOCKD(LFOCKD)
      DIMENSION WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_DIAG')
C
C------------------------------------------------------------------
C     Determine the symmetri of the result vector from the symmetry
C     of the trial vector ISYMTR, and the opperator symmtry ISYMOP.
C------------------------------------------------------------------
C
      ISYRES  = MULD2H(ISYMOP,ISYMTR)
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LDIAG1 = NT1AM(ISYMTR)
C
      KDIAG1  = 1
      KEND1   = KDIAG1 + LDIAG1
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_DIAG',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_DIAG',' ',KEND1,LWORK)

C---------------------------------------------------------
C     Calculate the diagonal parts of the RPA E[2] matrix.
C---------------------------------------------------------
C
      CALL SO_EDIAG1(WORK(KDIAG1),LDIAG1,FOCKD,LFOCKD,ISYRES)
C
C-----------------------------------------------------------
C     Open file and write E[2] diagonals to disk and output.
C-----------------------------------------------------------
C
      LUDIAG = -1
      CALL GPOPEN(LUDIAG,'RP_DIAG','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUDIAG)
C
      CALL WRITT(LUDIAG,LDIAG1,WORK(KDIAG1))
C
      IF (IPRSOP .GE. 7) THEN
         CALL AROUND('Fock contribution to diagonal 1p1h part of E[2]')
         CALL OUTPUT(WORK(KDIAG1),1,LDIAG1,1,1,LDIAG1,1,1,LUPRI)
      END IF
C
      CALL GPCLOSE(LUDIAG,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_DIAG')
C
      RETURN
      END
