#include <single.h>
C
C  /* Deck so_incred */
      SUBROUTINE SO_INCRED(NOLDTR,NNEWTR,ISYMTR,REDE,LREDE,REDS,LREDS,
     &                     LREDOL,LTR1E,LTR2E,IS,IE,IN,LI,JS,JE,JN,LJ,
     &                     WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, September 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculates the reduced eigenvalue problem.
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
C
      PARAMETER (ZERO = 0.0D+00, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THREE = 3.0D0)
C
C---------------------------------
C     Dimensions of the arguments.
C---------------------------------
C
      DIMENSION REDE(LREDE,LREDE),REDS(LREDS,LREDS)
      DIMENSION IS(LI), IE(LI), IN(LI), JS(LJ), JE(LJ), JN(LJ)
      DIMENSION WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_INCRED')
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LREDEH = LREDE  / 2
      LREDOH = LREDOL / 2
C
      LOFF1  = LREDOL * LREDOL
      LOFF2  = LREDOL * LREDOL
C
      KOFF1   = 1
      KOFF2   = KOFF1  + LOFF1
      KEND1   = KOFF2  + LOFF2
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_INCRED.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_INCRED.1',' ',KEND1,LWORK)
C 
C-----------------------------------------------------------------------
C     Move matrix elements from previous reduced eigenvalue problem to
C     the right positions in the matrices of the new eigenvalue problem.
C-----------------------------------------------------------------------
C
      CALL DCOPY(LOFF1,REDE,1,WORK(KOFF1),1)
      CALL DCOPY(LOFF2,REDS,1,WORK(KOFF2),1)
C
      DO JVEC = 1,NOLDTR
C
         DO IVEC = 1,NOLDTR
C
            IOFF1 = IVEC            + ((JVEC - 1) * LREDOL)        
            IOFF2 = (IVEC + LREDOH) + ((JVEC - 1) * LREDOL)       
            IOFF3 = IVEC            + ((JVEC - 1 + LREDOH) * LREDOL)
            IOFF4 = (IVEC + LREDOH) + ((JVEC - 1 + LREDOH) * LREDOL)
C
            REDE(IVEC,JVEC)               = WORK(KOFF1+IOFF1-1)
            REDE(IVEC+LREDEH,JVEC)        = WORK(KOFF1+IOFF2-1)
            REDE(IVEC,JVEC+LREDEH)        = WORK(KOFF1+IOFF3-1)
            REDE(IVEC+LREDEH,JVEC+LREDEH) = WORK(KOFF1+IOFF4-1)
C
            REDS(IVEC,JVEC)               = WORK(KOFF2+IOFF1-1)
            REDS(IVEC+LREDEH,JVEC)        = WORK(KOFF2+IOFF2-1)
            REDS(IVEC,JVEC+LREDEH)        = WORK(KOFF2+IOFF3-1)
            REDS(IVEC+LREDEH,JVEC+LREDEH) = WORK(KOFF2+IOFF4-1)
C
         END DO
C
      END DO
C
C---------------------------------
C     2. allocation of work space.
C---------------------------------
C
      LRES1  = LTR1E
      LRES2  = LTR2E
C
      KRES1   = 1
      KRES2   = KRES1  + LRES1
      KEND2   = KRES2  + LRES2
      LWORK2  = LWORK  - KEND2
C
      CALL SO_MEMMAX ('SO_INCRED.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('SO_INCRED.2',' ',KEND2,LWORK)
C
C---------------------------------------------------------------
C     Transform right hand side result vectors to left hands and
C     save on scratch files.
C     (A.L.2004  I think that the transformation is from Left to
C     Right)
C---------------------------------------------------------------
C
      CALL SO_OPEN(LURS1E,FNRS1E,LTR1E)
      CALL SO_OPEN(LURS1D,FNRS1D,LTR1E)
      CALL SO_OPEN(LURS2E,FNRS2E,LTR2E)
      CALL SO_OPEN(LURS2D,FNRS2D,LTR2E)
C
      CALL SO_OPEN(LUSC1E,FNSC1E,LTR1E)
      CALL SO_OPEN(LUSC1D,FNSC1D,LTR1E)
      CALL SO_OPEN(LUSC2E,FNSC2E,LTR2E)
      CALL SO_OPEN(LUSC2D,FNSC2D,LTR2E)
C
      DO INEWTR = 1,NNEWTR
C
         CALL SO_READ(WORK(KRES1),LTR1E,LURS1E,FNRS1E,INEWTR+NOLDTR)
         CALL SO_WRITE(WORK(KRES1),LTR1E,LUSC1E,FNSC1E,INEWTR)
C
         CALL SO_READ(WORK(KRES2),LTR2E,LURS2E,FNRS2E,INEWTR+NOLDTR)
C
         CALL CCSD_TCMEPKX(WORK(KRES2),TWO,ISYMTR)
         CALL DSCAL(LTR2E,HALF,WORK(KRES2),1)
C
         CALL SO_WRITE(WORK(KRES2),LTR2E,LUSC2E,FNSC2E,INEWTR)
C
      END DO
C
      DO INEWTR = 1,NNEWTR
C
         CALL SO_READ(WORK(KRES1),LTR1E,LURS1D,FNRS1D,INEWTR+NOLDTR)
         CALL SO_WRITE(WORK(KRES1),LTR1E,LUSC1D,FNSC1D,INEWTR)
C
         CALL SO_READ(WORK(KRES2),LTR2E,LURS2D,FNRS2D,INEWTR+NOLDTR)
C
         CALL CCSD_TCMEPKX(WORK(KRES2),TWO,ISYMTR)
         CALL DSCAL(LTR2E,HALF,WORK(KRES2),1)
C
         CALL SO_WRITE(WORK(KRES2),LTR2E,LUSC2D,FNSC2D,INEWTR)
C
      END DO
C
C---------------------------------
C     3. allocation of work space.
C---------------------------------
C
      NTRIAL = NNEWTR + NOLDTR
C
      NTRALL = NNEWTR +  (2 * NTRIAL)
C
      LOFF1  =  NTRIAL * NNEWTR
      LOFF2  =  NTRIAL * NNEWTR
C
      KOFF1   = 1
      KOFF2   = KOFF1  + LOFF1
      KEND3   = KOFF2  + LOFF2 
      LWORK3  = LWORK  - KEND3
C
C-------------------------------------------------------------------
C     Determine the length of each vector which can be held in core.
C-------------------------------------------------------------------
C
      LTR    = LWORK3 / NTRALL
C
      IF (LTR .LE. 0) CALL STOPIT('SO_INCRED.3',' ',KEND3+NTRALL,LWORK)
C
      N1READ = LTR1E / LTR
C
      L1LTR  = LTR1E - ( N1READ * LTR )
C
      N2READ = LTR2E / LTR
C
      L2LTR  = LTR2E - ( N2READ * LTR )
C
C---------------------------------
C     4. allocation of work space.
C---------------------------------
C
      LOFF3  = MIN((LTR * NTRIAL),(LTR2E * NTRIAL))
      LOFF4  = MIN((LTR * NTRIAL),(LTR2E * NTRIAL))
      LOFF5  = MIN((LTR * NNEWTR),(LTR2E * NNEWTR))
C
      KOFF3   = KEND3
      KOFF4   = KOFF3  + LOFF3
      KOFF5   = KOFF4  + LOFF4
      KEND4   = KOFF5  + LOFF5
      LWORK4  = LWORK  - KEND4
C
      CALL SO_MEMMAX ('SO_INCRED.4',LWORK4)
      IF (LWORK4 .LT. 0) CALL STOPIT('SO_INCRED.4',' ',KEND4,LWORK)
C
      CALL DZERO(WORK(KOFF1),LOFF1)
      CALL DZERO(WORK(KOFF2),LOFF2)
C
C-----------------------------------
C     Open files with trial vectors.
C-----------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1E)
C
      IOFF = 1 - LTR
C
      DO I1READ = 1,N1READ
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF3),LTR,NTRIAL,LUTR1E,FNTR1E,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF4),LTR,NTRIAL,LUTR1D,FNTR1D,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF5),LTR,NNEWTR,LUSC1E,FNSC1E,LTR1E,
     &                   IOFF)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,LTR,ONE,WORK(KOFF5),LTR,
     &              WORK(KOFF3),LTR,ONE,WORK(KOFF1),NNEWTR)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,LTR,ONE,WORK(KOFF5),LTR,
     &              WORK(KOFF4),LTR,ONE,WORK(KOFF2),NNEWTR)
C
         CALL SO_READSET(WORK(KOFF5),LTR,NNEWTR,LUSC1D,FNSC1D,LTR1E,
     &                   IOFF)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,LTR,ONE,WORK(KOFF5),LTR,
     &              WORK(KOFF4),LTR,ONE,WORK(KOFF1),NNEWTR)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,LTR,ONE,WORK(KOFF5),LTR,
     &              WORK(KOFF3),LTR,ONE,WORK(KOFF2),NNEWTR)
C
      END DO
C
      IF ( L1LTR .GT. 0 ) THEN
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF3),L1LTR,NTRIAL,LUTR1E,FNTR1E,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF4),L1LTR,NTRIAL,LUTR1D,FNTR1D,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF5),L1LTR,NNEWTR,LUSC1E,FNSC1E,LTR1E,
     &                   IOFF)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,L1LTR,ONE,WORK(KOFF5),L1LTR,
     &              WORK(KOFF3),L1LTR,ONE,WORK(KOFF1),NNEWTR)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,L1LTR,ONE,WORK(KOFF5),L1LTR,
     &              WORK(KOFF4),L1LTR,ONE,WORK(KOFF2),NNEWTR)
C
         CALL SO_READSET(WORK(KOFF5),L1LTR,NNEWTR,LUSC1D,FNSC1D,LTR1E,
     &                   IOFF)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,L1LTR,ONE,WORK(KOFF5),L1LTR,
     &              WORK(KOFF4),L1LTR,ONE,WORK(KOFF1),NNEWTR)
C 
         CALL DGEMM('T','N',NNEWTR,NTRIAL,L1LTR,ONE,WORK(KOFF5),L1LTR,
     &              WORK(KOFF3),L1LTR,ONE,WORK(KOFF2),NNEWTR)
C 
      END IF
C
C-----------------------------------
C     Open files with trial vectors.
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
         CALL SO_READSET(WORK(KOFF3),LTR,NTRIAL,LUTR2E,FNTR2E,LTR2E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF4),LTR,NTRIAL,LUTR2D,FNTR2D,LTR2E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF5),LTR,NNEWTR,LUSC2E,FNSC2E,LTR2E,
     &                   IOFF)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,LTR,ONE,WORK(KOFF5),LTR,
     &              WORK(KOFF3),LTR,ONE,WORK(KOFF1),NNEWTR)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,LTR,ONE,WORK(KOFF5),LTR,
     &              WORK(KOFF4),LTR,ONE,WORK(KOFF2),NNEWTR)
C
         CALL SO_READSET(WORK(KOFF5),LTR,NNEWTR,LUSC2D,FNSC2D,LTR2E,
     &                   IOFF)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,LTR,ONE,WORK(KOFF5),LTR,
     &              WORK(KOFF4),LTR,ONE,WORK(KOFF1),NNEWTR)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,LTR,ONE,WORK(KOFF5),LTR,
     &              WORK(KOFF3),LTR,ONE,WORK(KOFF2),NNEWTR)
C
      END DO
C
      IF ( L2LTR .GT. 0 ) THEN
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF3),L2LTR,NTRIAL,LUTR2E,FNTR2E,LTR2E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF4),L2LTR,NTRIAL,LUTR2D,FNTR2D,LTR2E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF5),L2LTR,NNEWTR,LUSC2E,FNSC2E,LTR2E,
     &                   IOFF)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,L2LTR,ONE,WORK(KOFF5),
     &              L2LTR,WORK(KOFF3),L2LTR,ONE,WORK(KOFF1),NNEWTR)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,L2LTR,ONE,WORK(KOFF5),
     &              L2LTR,WORK(KOFF4),L2LTR,ONE,WORK(KOFF2),NNEWTR)
C
         CALL SO_READSET(WORK(KOFF5),L2LTR,NNEWTR,LUSC2D,FNSC2D,LTR2E,
     &                   IOFF)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,L2LTR,ONE,WORK(KOFF5),
     &              L2LTR,WORK(KOFF4),L2LTR,ONE,WORK(KOFF1),NNEWTR)
C
         CALL DGEMM('T','N',NNEWTR,NTRIAL,L2LTR,ONE,WORK(KOFF5),
     &              L2LTR,WORK(KOFF3),L2LTR,ONE,WORK(KOFF2),NNEWTR)
C
      END IF
C
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
      CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
      CALL SO_CLOSE(LUSC1E,FNSC1E,'DELETE')
      CALL SO_CLOSE(LUSC1D,FNSC1D,'DELETE')
      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
      CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
      CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
      CALL SO_CLOSE(LUSC2E,FNSC2E,'DELETE')
      CALL SO_CLOSE(LUSC2D,FNSC2D,'DELETE')
C
C-------------------------------------------------------
C     Move new elements of reduced E[2] matrix in place.
C-------------------------------------------------------
C
      DO JTRVEC = 1,NTRIAL
C
         DO ITRVEC = 1,NNEWTR
C
            IOFF1 = ITRVEC + ( (JTRVEC - 1) * NNEWTR ) - 1
C
            REDE(ITRVEC+NOLDTR,JTRVEC) = WORK(KOFF1+IOFF1)
            REDE(JTRVEC,ITRVEC+NOLDTR) = WORK(KOFF1+IOFF1)
C
            REDE(ITRVEC+NOLDTR+NTRIAL,JTRVEC+NTRIAL) = WORK(KOFF1+IOFF1)
            REDE(JTRVEC+NTRIAL,ITRVEC+NOLDTR+NTRIAL) = WORK(KOFF1+IOFF1)
C
            REDE(ITRVEC+NOLDTR,JTRVEC+NTRIAL) = WORK(KOFF2+IOFF1)
            REDE(JTRVEC+NTRIAL,ITRVEC+NOLDTR) = WORK(KOFF2+IOFF1)
C
            REDE(ITRVEC+NOLDTR+NTRIAL,JTRVEC) = WORK(KOFF2+IOFF1)
            REDE(JTRVEC,ITRVEC+NOLDTR+NTRIAL) = WORK(KOFF2+IOFF1)
C
         END DO
C
      END DO
C
C--------------------------------------------------------------
C     The new elements of the reduced S[2] matrix is imported
C     from the orthogonalisation in SO_ORTH_TRN and inserted in
C     the reduced S[2] here.
C--------------------------------------------------------------
C
      DO JTRVEC = 1,NTRIAL
C
         DO ITRVEC = NOLDTR+1,NTRIAL
C
            IF (ITRVEC .EQ. JTRVEC) THEN
C
               REDS(ITRVEC,JTRVEC)               =  ONE
C
               REDS(ITRVEC+NTRIAL,JTRVEC)        = ZERO
C
               REDS(ITRVEC,JTRVEC+NTRIAL)        = ZERO
C
               REDS(ITRVEC+NTRIAL,JTRVEC+NTRIAL) = -ONE
C
            ELSE
C
               REDS(ITRVEC,JTRVEC)               = ZERO
               REDS(JTRVEC,ITRVEC)               = ZERO
C
               REDS(ITRVEC+NTRIAL,JTRVEC)        = ZERO
               REDS(JTRVEC,ITRVEC+NTRIAL)        = ZERO
C
               REDS(ITRVEC,JTRVEC+NTRIAL)        = ZERO
               REDS(JTRVEC+NTRIAL,ITRVEC)        = ZERO
C
               REDS(ITRVEC+NTRIAL,JTRVEC+NTRIAL) = ZERO
               REDS(JTRVEC+NTRIAL,ITRVEC+NTRIAL) = ZERO
C
            END IF
C
         END DO
C
      END DO
C
C---------------------------------------------------------------------
      IF ( IPRSOP .GE. 6 ) THEN
C
C---------------------------------------------
C        Print reduced E[2] and S[2] matrices.
C---------------------------------------------
C
         CALL AROUND('New reduced E[2] Matrix elements')
C
         CALL OUTPUT(REDE,1,LREDE,1,LREDE,LREDE,LREDE,1,LUPRI)
C
         CALL AROUND('New reduced S[2] Matrix elements')
C
         CALL OUTPUT(REDS,1,LREDS,1,LREDS,LREDS,LREDS,1,LUPRI)
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_INCRED')
C
      RETURN
      END
