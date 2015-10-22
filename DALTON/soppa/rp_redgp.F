C
C  /* Deck rp_redgp */
      SUBROUTINE RP_REDGP(NOLDTR,NNEWTR,ISYMTR,REDC,LREDC,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     A.L. December 2002
C
C     PURPOSE: Calculates the reduced gradient of the properties
C
C     NOLDTR    number of trial vectors keep from the preciding iteration
C     NNEWTR    number of new trial vectors added in this iteration
C     REDC      REDC(LREDC) will be the new (C^R +- C^R) array
C     LREDC     dimension of the REDC array
C     LTR1E     number of transition moments
C     WORK
C     LWORK
C
#include "implicit.h"
#include "soppinf.h"
#include "ccsdsym.h"
#include "priunit.h"
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
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_REDGP')
C
C----------------
C     Open files.
C----------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,NT1AM(ISYMTR))
      CALL SO_OPEN(LUTR1D,FNTR1D,NT1AM(ISYMTR))
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LTR1E  = NT1AM(ISYMTR)
      LGPVEC = 2*LTR1E
      NTRIAL = NNEWTR + NOLDTR
      NTRALL = NNEWTR +  (2 * NTRIAL)
C
      LOFF1 = NTRIAL
C
      KGPVEC = 1
      KOFF1  = KGPVEC+ LGPVEC
      KOFF2  = KOFF1 + LOFF1
      KEND1  = KOFF2 + LOFF1
      LWORK1 = LWORK - KEND1
C
C-------------------------------------------------------------------
C     Determine the length of each vector which can be held in core.
C-------------------------------------------------------------------
C
      LTR    = LWORK1 / NTRALL
C
      IF (LTR .LE. 0) CALL STOPIT('RP_REDGP.1',' ',KEND1+NTRALL,LWORK)
C
      N1READ = LTR1E / LTR
C
      L1LTR  = LTR1E - ( N1READ * LTR )
C
C---------------------------------
C     2. allocation of work space.
C---------------------------------
C
      LOFF3  = MIN((LTR * NTRIAL),(LTR1E * NTRIAL))
      LOFF4  = MIN((LTR * NTRIAL),(LTR1E * NTRIAL))
      LGPVCH = LGPVEC / 2
C
      KOFF3   = KEND1
      KOFF4   = KOFF3  + LOFF3
      KEND2   = KOFF4  + LOFF4
      LWORK2  = LWORK  - KEND2
C
      CALL SO_MEMMAX ('RP_REDGP.2',LWORK2)
      IF (LWORK2 .LT.0) CALL STOPIT('RP_REDGP.2',' ',KEND2,LWORK)
C
C-------------------------------------
C     Open file for property gradient.
C-------------------------------------
C
      CALL SO_OPEN(LUGPVE,FNGPVE,LGPVEC)
C
      CALL SO_READ(WORK(KGPVEC),LGPVEC,LUGPVE,FNGPVE,1)
C
C--------------------------------------
C     Close file for property gradient.
C--------------------------------------
C
      CALL SO_CLOSE(LUGPVE,FNGPVE,'KEEP')
C
      CALL DZERO(WORK(KOFF1),LOFF1)
      CALL DZERO(WORK(KOFF2),LOFF1)
      CALL DZERO(WORK(KOFF3),LOFF3)
      CALL DZERO(WORK(KOFF4),LOFF4)
C
      IOFF = 1 - LTR
C
      DO I1READ = 1,N1READ
C
         IOFF = IOFF + LTR
C
cLig <>  that are the 1b, 2b matrix 
         CALL SO_READSET(WORK(KOFF3),LTR,NTRIAL,LUTR1E,FNTR1E,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF4),LTR,NTRIAL,LUTR1D,FNTR1D,LTR1E,
     &                   IOFF)
C
cLig  dot product between the bT C
         CALL DGEMM('T','N',NTRIAL,1,LTR,ONE,WORK(KOFF3),LTR,
     &              WORK(KGPVEC),LTR,ONE,WORK(KOFF1),NTRIAL)
         CALL DGEMM('T','N',NTRIAL,1,LTR,ONE,WORK(KOFF4),LTR,
     &              WORK(KGPVEC+LTR1E),LTR,ONE,WORK(KOFF1),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,LTR,ONE,WORK(KOFF4),LTR,
     &              WORK(KGPVEC),LTR,ONE,WORK(KOFF2),NTRIAL)
         CALL DGEMM('T','N',NTRIAL,1,LTR,ONE,WORK(KOFF3),LTR,
     &              WORK(KGPVEC+LTR1E),LTR,ONE,WORK(KOFF2),NTRIAL)
C
      END DO
C
      IF ( L1LTR .GT. 0 ) THEN
C
         IOFF = IOFF + LTR
C
cLig look the commenti in the do loppo before ...
         CALL SO_READSET(WORK(KOFF3),L1LTR,NTRIAL,LUTR1E,FNTR1E,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF4),L1LTR,NTRIAL,LUTR1D,FNTR1D,LTR1E,
     &                   IOFF)
C
         CALL DGEMM('T','N',NTRIAL,1,L1LTR,ONE,WORK(KOFF3),L1LTR,
     &              WORK(KGPVEC),L1LTR,ONE,WORK(KOFF1),NTRIAL)
         CALL DGEMM('T','N',NTRIAL,1,L1LTR,ONE,WORK(KOFF4),L1LTR,
     &              WORK(KGPVEC+LTR1E),L1LTR,ONE,WORK(KOFF1),NTRIAL)
C
         CALL DGEMM('T','N',NTRIAL,1,L1LTR,ONE,WORK(KOFF4),L1LTR,
     &              WORK(KGPVEC),L1LTR,ONE,WORK(KOFF2),NTRIAL)
         CALL DGEMM('T','N',NTRIAL,1,L1LTR,ONE,WORK(KOFF3),L1LTR,
     &              WORK(KGPVEC+LTR1E),L1LTR,ONE,WORK(KOFF2),NTRIAL)
C
      END IF
C
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C
         DO ITRVEC = 1,NTRIAL
C
           REDC(ITRVEC) = WORK(KOFF1+ITRVEC-1)
C
           REDC(ITRVEC+NTRIAL) = WORK(KOFF2+ITRVEC-1)
C
         END DO
C
C-----------------------
C        Print reduced C
C-----------------------
C
      IF( IPRSOP .GT. 5) THEN
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
      CALL QEXIT('RP_REDGP')
C
      RETURN
      END
