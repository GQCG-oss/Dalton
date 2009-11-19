#include <single.h>
C
C  /* Deck so_diag */
      SUBROUTINE SO_DIAG(FOCKD,LFOCKD,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                   ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate diagonale parts of the SOPPA E[2] matrix 
C              (excluding first and second order contributions in
C              the one-particle part) and of the S[2] matrix.
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
      DIMENSION FOCKD(LFOCKD), DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_DIAG')
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
      LDIAG2 = NT2AM(ISYMTR)
C
      KDIAG1  = 1
      KDIAG2  = KDIAG1 + LDIAG1
      KEND1   = KDIAG2 + LDIAG2
      LWORK1  = LWORK   - KEND1
C
      CALL SO_MEMMAX ('SO_DIAG',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_DIAG',' ',KEND1,LWORK)

C---------------------------------------------------------
C     Calculate the diagonal one and two particle parts of
C     the SOPPA E[2] matrix (In the one-particle part are
C     first and second order contributions excluded)
C---------------------------------------------------------
C
      CALL SO_EDIAG1(WORK(KDIAG1),LDIAG1,FOCKD,LFOCKD,ISYRES)
C
      CALL SO_EDIAG2(WORK(KDIAG2),LDIAG2,FOCKD,LFOCKD,ISYRES,
     &               WORK(KEND1),LWORK1)
C
C-----------------------------------------------------------
C     Open file and write E[2] diagonals to disk and output.
C-----------------------------------------------------------
C
      LUDIAG = -1
      CALL GPOPEN(LUDIAG,'SO_DIAG','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUDIAG)
C
      CALL WRITT(LUDIAG,LDIAG1,WORK(KDIAG1))
C
      CALL WRITT(LUDIAG,LDIAG2,WORK(KDIAG2))
C
      IF (IPRSOP .GE. 7) THEN
         CALL AROUND('Fock contribution to diagonal 1p1h part of E[2]')
         CALL OUTPUT(WORK(KDIAG1),1,LDIAG1,1,1,LDIAG1,1,1,LUPRI)
         CALL AROUND('Diagonal 2p2h part of E[2]')
         CALL OUTPUT(WORK(KDIAG2),1,LDIAG2,1,1,LDIAG2,1,1,LUPRI)
      END IF
C
C----------------------------------------------------
C     Calculate the diagonal one particle part of the
C     SOPPA S[2] matrix.
C----------------------------------------------------
C
      CALL SO_SDIAG1(WORK(KDIAG1),LDIAG1,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &               ISYRES)
C
C------------------------------------------------------------
C     Write S[2] diagonals to disk and output and close file.
C------------------------------------------------------------
C
      CALL WRITT(LUDIAG,LDIAG1,WORK(KDIAG1))
C
      CALL GPCLOSE (LUDIAG,'KEEP')
C
      IF (IPRSOP .GE. 7) THEN
         CALL AROUND('Diagonal 1p1h part of S[2]')
         CALL OUTPUT(WORK(KDIAG1),1,LDIAG1,1,1,LDIAG1,1,1,LUPRI)
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_DIAG')
C
      RETURN
      END
