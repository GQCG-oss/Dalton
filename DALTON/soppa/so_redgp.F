C
C  /* Deck so_redgp */
      SUBROUTINE SO_REDGP(DOUBLES,NOLDTR,NNEWTR,ISYMTR,IMAGPROP,
     &                    REDC,LREDC,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Andrea Ligabue, December 2003
C     derived by SO_INCRED Keld Bak, September 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     REDC(LREDC)       the reduced C
C     PURPOSE: Calculates the reduced eigenvalue problem.
C
C     INPUT
C        DOUBLES    Are there a doubles contribution?
C        IMAGPROP   Are there gradient imaginary?
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      PARAMETER (ZERO = 0.0D+00, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THREE = 3.0D0)
C
C---------------------------------
C     Dimensions of the arguments.
C---------------------------------
C
      DIMENSION REDC(LREDC)
      DIMENSION WORK(LWORK)
      LOGICAL IMAGPROP, DOUBLES
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_REDGP')
C
C  Factor for d-ex gradient
      DFACTOR = -1.0D0
      IF(IMAGPROP) DFACTOR = ONE
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LTR1E   = NT1AM(ISYMTR)
      IF (DOUBLES) LTR2E   = N2P2HOP(ISYMTR)
      LGPVC1  = LTR1E
      IF (DOUBLES) LGPVC2  = LTR2E
C
      NTRIAL = NNEWTR + NOLDTR
C
      NTRALL = NNEWTR +  (2 * NTRIAL)
C
      LOFF1   = NTRIAL
      LOFF2   = NTRIAL
C
      KGPVC1  = 1
      IF(DOUBLES) THEN
         KGPVC2  = KGPVC1 + LGPVC1
         KOFF1   = KGPVC2 + LGPVC2
      ELSE
         KOFF1   = KGPVC1 + LGPVC1
      ENDIF
      KOFF2   = KOFF1  + LOFF1
      KEND1   = KOFF2  + LOFF2
      LWORK1  = LWORK  - KEND1
C
C-------------------------------------------------------------------
C     Determine the length of each vector which can be held in core.
C-------------------------------------------------------------------
C
      LTR    = LWORK1 / NTRALL
C
      IF (LTR .LE. 0) CALL STOPIT('SO_REDGP.1',' ',KEND1+NTRALL,LWORK)
C
      N1READ = LTR1E / LTR
C
      L1LTR  = LTR1E - ( N1READ * LTR )
      IF (DOUBLES)THEN
         N2READ = LTR2E / LTR
C
         L2LTR  = LTR2E - ( N2READ * LTR )
      ENDIF
C
C---------------------------------
C     2. allocation of work space.
C---------------------------------
C Are these excessive?
      IF(DOUBLES) THEN
         LOFF3  = MIN((LTR * NTRIAL),(LTR2E * NTRIAL))
         LOFF4  = MIN((LTR * NTRIAL),(LTR2E * NTRIAL))
      ELSE
         LOFF3  = MIN((LTR * NTRIAL),(LTR1E * NTRIAL))
         LOFF4  = MIN((LTR * NTRIAL),(LTR1E * NTRIAL))
      ENDIF
C
      KOFF3   = KEND1
      KOFF4   = KOFF3  + LOFF3
      KEND2   = KOFF4  + LOFF4
      LWORK2  = LWORK  - KEND2
C
      CALL SO_MEMMAX ('SO_REDGP.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('SO_REDGP.2',' ',KEND2,LWORK)
C
C------------------------------------------------------
C     Open, read and close files with property vectors.
C------------------------------------------------------
C
      CALL SO_OPEN(LUGPV1,FNGPV1,LTR1E)
      CALL SO_READ(WORK(KGPVC1),LGPVC1,LUGPV1,FNGPV1,1)
      CALL SO_CLOSE(LUGPV1,FNGPV1,'KEEP')
C
      IF (DOUBLES) THEN
         CALL SO_OPEN(LUGPV2,FNGPV2,LTR2E)
         CALL SO_READ(WORK(KGPVC2),LGPVC2,LUGPV2,FNGPV2,1)
         CALL SO_CLOSE(LUGPV2,FNGPV2,'KEEP')
      ENDIF
C
      IF(IPRSOP.GT.10) THEN
C
         CALL AROUND("Right GP vector inside SO_REDGP")
         CALL OUTPUT(WORK(KGPVC1),1,LTR1E,1,1,LTR1E,1,1,LUPRI)
         IF(DOUBLES)
     &         CALL OUTPUT(WORK(KGPVC2),1,LTR2E,1,1,LTR2E,1,1,LUPRI)
C
      ENDIF
C
      CALL DZERO(WORK(KOFF1),LOFF1)
      CALL DZERO(WORK(KOFF2),LOFF2)
cLigC I think I have not to make transofrmation becouse the b are L
cLigC and so the C must be right ... and are just right on the files
cLigC
cLigC-----------------------------------
cLigC     Transform RHS GP vector to LHS
cLigC-----------------------------------
cLigC
cLig       IF(IPRSOP.GT.100) THEN
cLigC
cLig          CALL AROUND('GP Vector before changing Right to left')
cLigC
cLig          CALL OUTPUT(WORK(KGPVC2),1,LGPVC2,1,1,LGPVC2,LGPVC2,1,LUPRI)
cLigC
cLig       ENDIF
cLigC
cLig       CALL SO_TMLTR(WORK(KGPVC2),HALF,ISYMTR)
cLig       CALL SO_TMLTR(WORK(KGPVC2+LTR2E),HALF,ISYMTR)
cLigC
cLig       IF(IPRSOP.GT.10) THEN
cLigC
cLig          CALL AROUND('GP Vector after changing Right to left')
cLigC
cLig          CALL OUTPUT(WORK(KGPVC2),1,LGPVC2,1,1,LGPVC2,LGPVC2,1,LUPRI)
cLigC
cLig       ENDIF
C
C
C-------------------------------------------
C     Open files with singles trial vectors.
C-------------------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1E)
C
cLig since we create 1 NNEWTR each run we just go to compute 1 number
cLig and put it in WORK(KOFF1) (that will be 1C1E)
cLig ... and the same for WORK(KOFF2) (that will be 2C1D)
C
      IOFF = 1 - LTR
C
      DO I1READ = 1,N1READ
C
         IOFF = IOFF + LTR
C
cLig  <> read all (1 to 3) 1b1E in WORK(KOFF3) and 2b1D in WORK(KOFF4)
         CALL SO_READSET(WORK(KOFF3),LTR,NTRIAL,LUTR1E,FNTR1E,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF4),LTR,NTRIAL,LUTR1D,FNTR1D,LTR1E,
     &                   IOFF)
C
cLig  <> 1b1E * 1C1E and 2b1D * 2C1D (that is +- 1C1E) for eacch (1-3= b
cLig         CALL DGEMM('T','N',NTRIAL,1,LTR,ONE,WORK(KOFF3),LTR,
cLig     &              WORK(KGPVC1),LTR1E,ONE,WORK(KOFF1),NTRIAL)
         CALL DGEMM('T','N',NTRIAL,1,LTR,ONE,WORK(KOFF3),LTR,
     &              WORK(KGPVC1),LTR,ONE,WORK(KOFF1),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,LTR,DFACTOR,WORK(KOFF4),LTR,
     &              WORK(KGPVC1),LTR,ONE,WORK(KOFF1),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,LTR,ONE,WORK(KOFF4),LTR,
     &              WORK(KGPVC1),LTR,ONE,WORK(KOFF2),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,LTR,DFACTOR,WORK(KOFF3),LTR,
     &              WORK(KGPVC1),LTR,ONE,WORK(KOFF2),NTRIAL)
C
      END DO
C
      IF ( L1LTR .GT. 0 ) THEN
C
         IOFF = IOFF + LTR
C
cLig  <> read 1b1E in WORK(KOFF3) and 2b1D in WORK(KOFF4)
         CALL SO_READSET(WORK(KOFF3),L1LTR,NTRIAL,LUTR1E,FNTR1E,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF4),L1LTR,NTRIAL,LUTR1D,FNTR1D,LTR1E,
     &                   IOFF)
C
cLig  <> 1b1E * 1C1E and 2b1D * 2C1D (that is +- 1C1E)
cLig         CALL DGEMM('T','N',NTRIAL,1,L1LTR,ONE,WORK(KOFF3),L1LTR,
cLig     &              WORK(KGPVC1),LTR1E,ONE,WORK(KOFF1),NTRIAL)
         CALL DGEMM('T','N',NTRIAL,1,L1LTR,ONE,WORK(KOFF3),L1LTR,
     &              WORK(KGPVC1),LTR,ONE,WORK(KOFF1),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,L1LTR,DFACTOR,WORK(KOFF4),L1LTR,
     &              WORK(KGPVC1),LTR,ONE,WORK(KOFF1),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,L1LTR,ONE,WORK(KOFF4),L1LTR,
     &              WORK(KGPVC1),LTR,ONE,WORK(KOFF2),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,L1LTR,DFACTOR,WORK(KOFF3),L1LTR,
     &              WORK(KGPVC1),LTR,ONE,WORK(KOFF2),NTRIAL)
C
      END IF
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C
C---------------------------------
C     Jump to end if singles only.
C---------------------------------
C
      IF (.NOT.DOUBLES) GOTO 1010
C-----------------------------------
C     Open files with doubles trial vectors.
C-----------------------------------
C
      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2E)
C
      IOFF = 1 - LTR
C
      DO I2READ = 1,N2READ
C
         IOFF = IOFF + LTR
C
cLig  <> read 1b2E in WORK(KOFF3) and 2b2D in WORK(KOFF4)
         CALL SO_READSET(WORK(KOFF3),LTR,NTRIAL,LUTR2E,FNTR2E,LTR2E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF4),LTR,NTRIAL,LUTR2D,FNTR2D,LTR2E,
     &                   IOFF)
C
cLig  <> 1b2E * 1C2E and 2b2D * 2C2D (that is +- 1C1E)
         CALL DGEMM('T','N',NTRIAL,1,LTR,ONE,WORK(KOFF3),LTR,
     &              WORK(KGPVC2),LTR,ONE,WORK(KOFF1),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,LTR,DFACTOR,WORK(KOFF4),LTR,
     &              WORK(KGPVC2),LTR,ONE,WORK(KOFF1),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,LTR,ONE,WORK(KOFF4),LTR,
     &              WORK(KGPVC2),LTR,ONE,WORK(KOFF2),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,LTR,DFACTOR,WORK(KOFF3),LTR,
     &              WORK(KGPVC2),LTR,ONE,WORK(KOFF2),NTRIAL)
C
      END DO
C
      IF ( L2LTR .GT. 0 ) THEN
C
         IOFF = IOFF + LTR
C
cLig  <> read 1b2E in WORK(KOFF3) and 2b2D in WORK(KOFF4)
         CALL SO_READSET(WORK(KOFF3),L2LTR,NTRIAL,LUTR2E,FNTR2E,LTR2E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF4),L2LTR,NTRIAL,LUTR2D,FNTR2D,LTR2E,
     &                   IOFF)
C
cLig  <> 1b2E * 1C2E and 2b2D * 2C2D (that is +- 1C1E)
         CALL DGEMM('T','N',NTRIAL,1,L2LTR,ONE,WORK(KOFF3),L2LTR,
     &              WORK(KGPVC2),L2LTR,ONE,WORK(KOFF1),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,L2LTR,DFACTOR,WORK(KOFF4),L2LTR,
     &              WORK(KGPVC2),L2LTR,ONE,WORK(KOFF1),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,L2LTR,ONE,WORK(KOFF4),L2LTR,
     &              WORK(KGPVC2),L2LTR,ONE,WORK(KOFF2),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,L2LTR,DFACTOR,WORK(KOFF3),L2LTR,
     &              WORK(KGPVC2),L2LTR,ONE,WORK(KOFF2),NTRIAL)
C
      END IF
C
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
C  Continuation point in the end
1010  CONTINUE

C
C----------------------------------------------------
C     Move the elment in the right places of vector C
C----------------------------------------------------
C
      DO ITRVEC = 1,NTRIAL
C
         REDC(ITRVEC) = WORK(KOFF1+ITRVEC-1)
C
         REDC(ITRVEC+NTRIAL) = WORK(KOFF2+ITRVEC-1)
C
      END DO
C
      IF ( IPRSOP .GE. 6 ) THEN
C
C-----------------------
C        Print reduced C
C-----------------------
C
         CALL AROUND('New reduced C elements')
C
         CALL OUTPUT(REDC,1,LREDC,1,1,1,LREDC,1,LUPRI)
C
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_REDGP')
C
      RETURN
      END
