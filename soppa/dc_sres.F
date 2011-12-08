C
C  /* Deck dc_sres */
      SUBROUTINE DC_SRES(NOLDTR, NNEWTR,  DENSIJ, LDENSIJ, 
     &                   DENSAB, LDENSAB, ISYMTR, WORK,   LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, June 1997.
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Driver routine for making a linear transformation of
C              a trialvector with the second order part of the
C              SOPPA overlap matrix S[2].
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION WORK(LWORK)
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <soppinf.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('DC_SRES')
C
C------------------------------------------------------------------
C     Determine the symmetri of the result vector from the symmetry
C     of the trial vector ISYMTR, and the opperator symmtry ISYMOP.
C------------------------------------------------------------------
C
      ISYRES  = MULD2H(ISYMOP,ISYMTR)
C
C---------------------------------
C     Work space allocation no. 1.
C---------------------------------
C
      LTR1E   = NT1AM(ISYMTR)
      LTR1D   = NT1AM(ISYMTR)
      LRESO1E = NT1AM(ISYMTR)
      LRESO1D = NT1AM(ISYMTR)
C
      KTR1E   = 1
      KTR1D   = KTR1E   + LTR1E
      KRESO1E = KTR1D   + LTR1D
      KRESO1D = KRESO1E + LRESO1E
      KEND1   = KRESO1D + LRESO1D
      LWORK1  = LWORK   - KEND1
C
      CALL SO_MEMMAX ('DC_SRES',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('DC_SRES',' ',KEND1,LWORK)
C
C----------------------------------------------
C     Open files with trial and result vectors.
C----------------------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      CALL SO_OPEN(LURO1E,FNRO1E,LRESO1E)
      CALL SO_OPEN(LURO1D,FNRO1D,LRESO1D)
C
C================================================
C     Loop over number of excitations considered.
C================================================
C
      DO 100 INEWTR = 1,NNEWTR
C
C-------------------------------------------------
C        Determine pointer to INEWTR trial vector.
C-------------------------------------------------
C
         INEW = NOLDTR + INEWTR
C
C--------------------------
C        Read trial vector.
C--------------------------
C
         CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,INEW)
         CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,INEW)
C
C-------------------------------
C        Transform trial vector.
C-------------------------------
C
         CALL DC_RES_O(WORK(KRESO1E),LRESO1E,WORK(KRESO1D),LRESO1D,
     &                 WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                 DENSIJ,LDENSIJ,DENSAB,LDENSAB,ISYRES,ISYMTR)
C
C----------------------------------------
C        Write new resultvectors to file.
C----------------------------------------
C
         CALL SO_WRITE(WORK(KRESO1E),LRESO1E,LURO1E,FNRO1E,INEW)
         CALL SO_WRITE(WORK(KRESO1D),LRESO1D,LURO1D,FNRO1D,INEW)
C
  100 CONTINUE
C
C==================================
C     End of loop over excitations.
C==================================
C
      IF ( IPRSOP. GE. 7 ) THEN
C
C------------------------------------------
C        Write new resultvectors to output.
C------------------------------------------
C
         DO 200 INEWTR = 1,NNEWTR
C
C----------------------------------------------------
C           Determine pointer to INEWTR trial vector.
C----------------------------------------------------
C
            INEW = NOLDTR + INEWTR
C
            WRITE(LUPRI,'(/,I3,A)') INEWTR,
     &          '. new S[2] linear transformed trial vector'
C
            CALL SO_READ(WORK(KRESO1E),LRESO1E,LURO1E,FNRO1E,INEW)
            CALL SO_READ(WORK(KRESO1D),LRESO1D,LURO1D,FNRO1D,INEW)
C
            WRITE(LUPRI,'(I8,1X,F14.8)')
     &           (I,WORK(KRESO1E+I-1),I=1,LRESO1E)
            WRITE(LUPRI,'(I8,1X,F14.8)') 
     &           (I,WORK(KRESO1D+I-1),I=1,LRESO1D)
C
  200    CONTINUE
C
      END IF
C
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      CALL SO_CLOSE(LURO1E,FNRO1E,'KEEP')
      CALL SO_CLOSE(LURO1D,FNRO1D,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('DC_SRES')
C
      RETURN
      END
