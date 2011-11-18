C
C  /* Deck dc_calc */
      SUBROUTINE DC_CALC(ISYMTR,NEXCI,EXVAL,LEXVAL,
     &                   DENSIJ,LDENSIJ,DENSAB,DENSAI,LDENSAI,
     &                   LDENSAB,T2MP,LT2MP,FOCKD,LFOCKD,
     &                   WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, June 1997
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Determine "double corrected RPA" excitation energies.
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      LOGICAL   NONEWT
C
      DIMENSION EXVAL(LEXVAL)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION T2MP(LT2MP),     FOCKD(LFOCKD),   WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('DC_CALC')
C
C---------------------------------------------------------
C     Make E[2] linear transformation of RPA eigenvectors.
C---------------------------------------------------------
C
      DTIME      = SECOND()
      CALL DC_ERES(0,NEXCI,DENSIJ,LDENSIJ,DENSAB,LDENSAB,DENSAI,LDENSAI,
     &             T2MP,LT2MP,FOCKD,LFOCKD,ISYMTR,WORK,LWORK)
      DTIME      = SECOND()   - DTIME
      SOTIME(35) = SOTIME(35) + DTIME
C
C-----------------------------------------------------------
C     Make S[2] linear transformation of trialvectors giving
C     resultvectors.
C-----------------------------------------------------------
C
      DTIME      = SECOND()
      CALL DC_SRES(0,NEXCI,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &             ISYMTR,WORK,LWORK)
      DTIME      = SECOND()   - DTIME
      SOTIME(40) = SOTIME(40) + DTIME
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LTR1E = NT1AM(ISYMTR)
      LTR1D = NT1AM(ISYMTR)
      LTR2E = NT2AM(ISYMTR)
      LTR2D = NT2AM(ISYMTR)
C
      KTR1    = 1
      KRES1   = KTR1   + LTR1E
      KRESO1  = KRES1  + LTR1E
      KRES2   = KRESO1 + LTR1E
      KEND1   = KRES2  + LTR2E
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('DC_CALC',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('DC_CALC',' ',KEND1,LWORK)
C
C----------------
C     Open files.
C----------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
      CALL SO_OPEN(LURS1E,FNRS1E,LTR1E)
      CALL SO_OPEN(LURS1D,FNRS1D,LTR1D)
      CALL SO_OPEN(LURO1E,FNRO1E,LTR1E)
      CALL SO_OPEN(LURO1D,FNRO1D,LTR1D)
      CALL SO_OPEN(LURS2E,FNRS2E,LTR2E)
      CALL SO_OPEN(LURS2D,FNRS2D,LTR2D)
C
C---------------------------
C     Loop over excitations.
C---------------------------
C
      DO 100 IEXCI = 1,NEXCI
C
C-----------------------------------------------------------------------
C        Calculate 1p-1h part of double corrected RPA excitation energy.
C-----------------------------------------------------------------------
C
         CALL SO_READ(WORK(KTR1),  LTR1E,LUTR1E,FNTR1E,IEXCI)
         CALL SO_READ(WORK(KRES1), LTR1E,LURS1E,FNRS1E,IEXCI)
         CALL SO_READ(WORK(KRESO1),LTR1E,LURO1E,FNRO1E,IEXCI)
C
         CALL DAXPY(LTR1E,-EXVAL(IEXCI),WORK(KRESO1),1,WORK(KRES1),1)
C
         EX1P1H = DDOT(LTR1E,WORK(KRES1),1,WORK(KTR1),1) 
C
C-----------------------------------------------------------------------
C        Calculate 1h-1p part of double corrected RPA excitation energy.
C-----------------------------------------------------------------------
C
         CALL SO_READ(WORK(KTR1),  LTR1D,LUTR1D,FNTR1D,IEXCI)
         CALL SO_READ(WORK(KRES1), LTR1D,LURS1D,FNRS1D,IEXCI)
         CALL SO_READ(WORK(KRESO1),LTR1D,LURO1D,FNRO1D,IEXCI)
C
         CALL DAXPY(LTR1D,-EXVAL(IEXCI),WORK(KRESO1),1,WORK(KRES1),1)
C
         EX1H1P = DDOT(LTR1D,WORK(KRES1),1,WORK(KTR1),1) 
C
C-----------------------------------------------------------------------
C        Calculate 2p-2h part of double corrected RPA excitation energy.
C-----------------------------------------------------------------------
C
         CALL SO_READ(WORK(KRES2), LTR2E,LURS2E,FNRS2E,IEXCI)
C
         CALL CC_OMEC(EX2P2H,WORK(KRES2),EXVAL(IEXCI),
     &                WORK(KEND1),LWORK1,ISYMTR)
C
C--------------------------------------------------------------------
C        Calculate the first order 2p-2h eigenvector in RPA(D) theory
C        and write to output.
C--------------------------------------------------------------------
C
         CALL DC_R1VEC(WORK(KRES2), LTR2E,EXVAL(IEXCI),ISYMTR,
     &                 WORK(KEND1),LWORK1)
C
         CALL SO_WRITE(WORK(KRES2), LTR2E, LUTR2E,FNTR2E, IEXCI)
C
C-----------------------------------------------------------------------
C        Calculate 2h-2p part of double corrected RPA excitation energy.
C-----------------------------------------------------------------------
C
         CALL SO_READ(WORK(KRES2), LTR2D,LURS2D,FNRS2D,IEXCI)
C
         CALL CC_OMEC(EX2H2P,WORK(KRES2),EXVAL(IEXCI),
     &                WORK(KEND1),LWORK1,ISYMTR)
C
C--------------------------------------------------------------------
C        Calculate the first order 2h-2p eigenvector in RPA(D) theory
C        and write to output.
C--------------------------------------------------------------------
C
         CALL DC_R1VEC(WORK(KRES2), LTR2D,EXVAL(IEXCI),ISYMTR,
     &                 WORK(KEND1),LWORK1)
C
         CALL SO_WRITE(WORK(KRES2), LTR2D, LUTR2D,FNTR2D, IEXCI)
C
C-----------------------------------------------------------
C        Add contributions to give RPA(D) excitation energy.
C-----------------------------------------------------------
C
         EXVAL(IEXCI) = EX1P1H + EX1H1P + EX2P2H/TWO + EX2H2P/TWO
C
  100 CONTINUE
C
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
      CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
      CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
      CALL SO_CLOSE(LURO1E,FNRO1E,'KEEP')
      CALL SO_CLOSE(LURO1D,FNRO1D,'KEEP')
      CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
      CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
C
C------------------------------------
C     Flush the standard output file.
C------------------------------------
C
      CALL FLSHFO(LUPRI)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('DC_CALC')
C
      RETURN
      END
