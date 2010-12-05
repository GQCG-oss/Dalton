C  /* Deck so_check */
      SUBROUTINE SO_CHECK(DENSIJ,LDENSIJ,DENSAB,LDENSAB,T2MP,LT2MP,
     &                    FOCKD,LFOCKD,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan P.A. Sauer.                    1-Nov-1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     Test routine that calculates the E[2] and S[2] matrices
C     explicitly by carrying out linear transformations on unit trial
C     vectors.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), T2MP(LT2MP)
      DIMENSION FOCKD(LFOCKD),   WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <soppinf.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_CHECK')
C
      CHKTIM = SECOND()
      RESTIM = ZERO
C
C-----------------------------------------------------
C     Work space allocation for E[2] and S[2] matrices
C     and trial vectors.
C-----------------------------------------------------
C
      NVAR   = 2 * (NT1AM(ISYMTR)+NT2AM(ISYMTR))
C
      LE2MAT = NVAR * NVAR
      LS2MAT = NVAR * NVAR
      LTR1E  = NT1AM(ISYMTR)
      LTR1D  = NT1AM(ISYMTR)
      LTR2E  = NT2AM(ISYMTR)
      LTR2D  = NT2AM(ISYMTR)
C
      KE2MAT  = 1
      KS2MAT  = KE2MAT + LE2MAT
      KTR1E   = KS2MAT + LS2MAT
      KTR1D   = KTR1E  + LTR1E
      KTR2E   = KTR1D  + LTR1D
      KTR2D   = KTR2E  + LTR2E
      KEND1   = KTR2D  + LTR2D
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_CHECK.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_CHECK.1',' ',KEND1,LWORK)
C
C---------------------------------------------------
C     Loop over columns of the E[2] and S[2] matrix.
C---------------------------------------------------
C
      DO 100 ICOL = 1,NVAR
C
C---------------------------------------------------
C        Work space allocation for solution vectors.
C---------------------------------------------------
C
         LE2MAT  = (ICOL - 1) * NVAR
         LS2MAT  = (ICOL - 1) * NVAR
         LRES1E  = NT1AM(ISYMTR)
         LRES1D  = NT1AM(ISYMTR)
         LRES2E  = NT2AM(ISYMTR)
         LRES2D  = NT2AM(ISYMTR)
         LRESO1E = NT1AM(ISYMTR)
         LRESO1D = NT1AM(ISYMTR)
         LRESO2E = NT2AM(ISYMTR)
         LRESO2D = NT2AM(ISYMTR)
C
         KRES1E  = KE2MAT  + LE2MAT
         KRES1D  = KRES1E  + LRES1E
         KRES2E  = KRES1D  + LRES1D
         KRES2D  = KRES2E  + LRES2E
         KRESO1E = KS2MAT  + LS2MAT
         KRESO1D = KRESO1E + LRESO1E
         KRESO2E = KRESO1D + LRESO1D
         KRESO2D = KRESO2E + LRESO2E
         KEND2   = KRESO2D + LRESO2D
         LWORK2  = LWORK   - KEND2
C
         CALL SO_MEMMAX ('SO_CHECK.2',LWORK2)
         IF (LWORK2 .LT. 0) CALL STOPIT('SO_CHECK.2',' ',KEND2,LWORK)
C
         CALL DZERO(WORK(KTR1E),LTR1E)
         CALL DZERO(WORK(KTR1D),LTR1D)
         CALL DZERO(WORK(KTR2E),LTR2E)
         CALL DZERO(WORK(KTR2D),LTR2D)
C
         IF (ICOL.LE.LTR1E) THEN
            WORK(KTR1E + ICOL - 1) = ONE
         ELSE IF (ICOL.LE.(LTR1E+LTR1D)) THEN
            WORK(KTR1D + (ICOL - LTR1E) - 1) = ONE
         ELSE IF (ICOL.LE.(LTR1E+LTR1D+LTR2E)) THEN
            WORK(KTR2E + (ICOL - LTR1E - LTR1D) - 1) = ONE
         ELSE 
            WORK(KTR2D + (ICOL - LTR1E - LTR1D - LTR2E) - 1) = ONE
         ENDIF
C
         CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
         CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
         CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
         CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
C
         CALL SO_WRITE(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,1)
         CALL SO_WRITE(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,1)
         CALL SO_WRITE(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,1)
         CALL SO_WRITE(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,1)
C
         CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
         CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
         CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
         CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
C
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
         CALL SO_ERES('AOSOP',0, 1, DENSIJ, LDENSIJ, DENSAB, LDENSAB,
     &               T2MP, LT2MP,  FOCKD,   LFOCKD, ISYMTR,
     &               WORK(KEND1), LWORK1)
         ETIME      = SECOND()   - DTIME
         SOTIME(35) = SOTIME(35) + ETIME
         DTIME      = SECOND()
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(1)  = SOWTIM(1)  + WTIMEE - WTIMES
C
         CALL GETTIM (DUMMY,WTIMES)
         CALL SO_SRES(0,1,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                ISYMTR,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(40) = SOTIME(40) + DTIME
         RESTIM = RESTIM + SECOND() - DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(2)  = SOWTIM(2)  + WTIMEE - WTIMES
C
         CALL SO_OPEN(LURS1E,FNRS1E,LRES1E)
         CALL SO_OPEN(LURS1D,FNRS1D,LRES1D)
         CALL SO_OPEN(LURS2E,FNRS2E,LRES2E)
         CALL SO_OPEN(LURS2D,FNRS2D,LRES2D)
         CALL SO_OPEN(LURO1E,FNRO1E,LRESO1E)
         CALL SO_OPEN(LURO1D,FNRO1D,LRESO1D)
C
         CALL SO_READ(WORK(KRES1E), LRES1E, LURS1E,FNRS1E,1)
         CALL SO_READ(WORK(KRES1D), LRES1D, LURS1D,FNRS1D,1)
         CALL SO_READ(WORK(KRES2E), LRES2E, LURS2E,FNRS2E,1)
         CALL SO_READ(WORK(KRES2D), LRES2D, LURS2D,FNRS2D,1)
         CALL SO_READ(WORK(KRESO1E),LRESO1E,LURO1E,FNRO1E,1)
         CALL SO_READ(WORK(KRESO1D),LRESO1D,LURO1D,FNRO1D,1)
         CALL DZERO(WORK(KRESO2E),LRESO2E)
         CALL DZERO(WORK(KRESO2D),LRESO2D)
C
         CALL SO_CLOSE(LURS1E,FNRS1E,'DELETE')
         CALL SO_CLOSE(LURS1D,FNRS1D,'DELETE')
         CALL SO_CLOSE(LURS2E,FNRS2E,'DELETE')
         CALL SO_CLOSE(LURS2D,FNRS2D,'DELETE')
         CALL SO_CLOSE(LURO1E,FNRO1E,'DELETE')
         CALL SO_CLOSE(LURO1D,FNRO1D,'DELETE')
C
  100 CONTINUE
C
      CALL AROUND('E[2] Matrix')
C
      WRITE(LUPRI,'(A,I2,A,I8)')' SYMMETRY :',ISYMTR,'  DIMENSION ',NVAR
      CALL OUTPUT(WORK(KE2MAT),1,NVAR,1,NVAR,NVAR,NVAR,1,LUPRI)
C
      CALL AROUND('S[2] Matrix')
C
      WRITE(LUPRI,'(A,I2,A,I8)')' SYMMETRY :',ISYMTR,'  DIMENSION ',NVAR
      CALL OUTPUT(WORK(KS2MAT),1,NVAR,1,NVAR,NVAR,NVAR,1,LUPRI)
C
      CHKTIM     = SECOND() - CHKTIM - RESTIM
      SOTIME(27) = SOTIME(27) + CHKTIM
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_CHECK')
C
      RETURN
      END
