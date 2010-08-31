C  /* Deck so_test */
      SUBROUTINE SO_TEST(NOLDTR,  NNEWTR, ISYMTR,  DENSIJ, 
     &                   LDENSIJ, DENSAB, LDENSAB, T2MP,   LT2MP, 
     &                   FOCKD,   LFOCKD, WORK,    LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, August 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Test of orhtonormality of trialvectors and the
C              correctnes of the linear transformed vectors 
C              which are generated as linear combinations of
C              previously calculated linear transformed vectors.
C
#include <implicit.h>
#include <priunit.h>
#include <soppinf.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO = 0.0D+00, ONE = 1.0D0)
C
C---------------------------------
C     Dimensions of the arguments.
C---------------------------------
C
      DIMENSION WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TEST')
C
C-----------------------------------------------------------------
C     Work space allocation no. 1.
C     Notice that the E[2] linear transformed trial vector and the 
C     S[2] linear transformed trial vector are of equal length and
C     that they use the same work space.
C-----------------------------------------------------------------
C
      LTR1E    = NT1AM(ISYMTR)
      LTR1D    = NT1AM(ISYMTR)
      LTR2E    = NT2AM(ISYMTR)
      LTR2D    = NT2AM(ISYMTR)
      LRES1E   = NT1AM(ISYMTR)
      LRES1D   = NT1AM(ISYMTR)
      LRES2E   = NT2AM(ISYMTR)
      LRES2D   = NT2AM(ISYMTR)
      LRESO1E  = LRES1E
      LRESO1D  = LRES1D
      LRESO2E  = LRES2E
      LRESO2D  = LRES2D
C
      KTR1EI   = 1
      KTR1DI   = KTR1EI   + LTR1E
      KTR2EI   = KTR1DI   + LTR1D
      KTR2DI   = KTR2EI   + LTR2E 
      KTR1EJ   = KTR2DI   + LTR2D
      KTR1DJ   = KTR1EJ   + LTR1E
      KTR2EJ   = KTR1DJ   + LTR1D
      KTR2DJ   = KTR2EJ   + LTR2E 
      KRES1EI  = KTR2DJ   + LTR2D
      KRES1DI  = KRES1EI  + LRES1E
      KRES2EI  = KRES1DI  + LRES1D
      KRES2DI  = KRES2EI  + LRES2E
      KRES1EJ  = KRES2DI  + LRES2D
      KRES1DJ  = KRES1EJ  + LRES1E
      KRES2EJ  = KRES1DJ  + LRES1D
      KRES2DJ  = KRES2EJ  + LRES2E
      KRESO1EI = KRES2DJ  + LRES2D
      KRESO1DI = KRESO1EI + LRES1E
      KRESO1EJ = KRESO1DI + LRES1D
      KRESO1DJ = KRESO1EJ + LRES1E
      KEND1    = KRESO1DJ + LRES1D
      LWORK1   = LWORK    - KEND1
C
      CALL SO_MEMMAX ('SO_TEST.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_TEST.1',' ',KEND1,LWORK)
C
      CALL SO_TEST1(WORK(KTR1EI),WORK(KTR1DI),WORK(KTR2EI),
     &              WORK(KTR2DI),WORK(KTR1EJ),WORK(KTR1DJ),WORK(KTR2EJ),
     &              WORK(KTR2DJ),
     &              WORK(KRES1EI),WORK(KRES1DI),WORK(KRES2EI), 
     &              WORK(KRES2DI),
     &              WORK(KRES1EJ),WORK(KRES1DJ),WORK(KRES2EJ), 
     &              WORK(KRES2DJ),
     &              WORK(KRESO1EI), WORK(KRESO1DI),
     &              WORK(KRESO1EJ), WORK(KRESO1DJ),
     &              LTR1E,LTR1D,LTR2E,LTR2D,LRES1E,LRES1D,
     &              LRES2E,LRES2D,LRESO1E,LRESO1D,NOLDTR,NNEWTR,
     &              DENSIJ, LDENSIJ, DENSAB, LDENSAB,
     &              T2MP,   LT2MP,   FOCKD,  LFOCKD,
     &              ISYMTR, WORK(KEND1),     LWORK1)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TEST')
C
      RETURN
      END
