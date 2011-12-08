C  /* Deck so_optvec */
      SUBROUTINE SO_OPTVEC(NOLDTR,NNEWTR,NEXCI,EIVEC,LEIVEC,ISYMTR,
     &                     WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, April 1997
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Generate optimized trial vectors and corresponding
C              linear transformed trial vectors.
C
#include <implicit.h>
#include <priunit.h>
#include <soppinf.h>
#include <ccsdsym.h>
C
      PARAMETER   (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER   (D100 = 100.0D0)
      PARAMETER   (THRNM = 1.0D-15, OVLMIN = 1.0D-4)
      CHARACTER*3 YES, NO
C
      DIMENSION EIVEC(LEIVEC,LEIVEC)
      DIMENSION WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_OPTVEC')
C
C----------------------------------------------------------
C     Determine the number of previous trialvectors and the
C     number of trialvectors they are transformed to.
C----------------------------------------------------------
C 
      LTR1E  = NT1AM(ISYMTR)
      LTR1D  = NT1AM(ISYMTR)
      LTR2E  = NT2AM(ISYMTR)
      LTR2D  = NT2AM(ISYMTR)
C
      NTRIAL = NOLDTR + NNEWTR
      NTRMX  = MIN(NTRIAL, (NSAVMX-1)*NEXCI )
C
C-------------------------------------------------------------------
C     Determine the length of each vector which can be held in core.
C-------------------------------------------------------------------
C
      NTRALL = 2 * (NTRIAL + NTRMX)
C
      LTR    = LWORK / NTRALL
C
      IF (LTR .LE. 0) CALL STOPIT('SO_OPTVEC.0',' ',NTRALL,LWORK)
C
      LTR    = MAX ( LTR, 4000 )
C
      N1READ = LTR1E / LTR
C
      L1LTR  = LTR1E - ( N1READ * LTR )
C
      N2READ = LTR2E / LTR
C
      L2LTR  = LTR2E - ( N2READ * LTR )
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      LOFF1  = MIN((LTR * NTRIAL),(LTR2E * NTRIAL))
      LOFF2  = MIN((LTR * NTRIAL),(LTR2E * NTRIAL))
      LOFF3  = MIN((LTR * NTRMX), (LTR2E * NTRMX))
      LOFF4  = MIN((LTR * NTRMX), (LTR2E * NTRMX))
C
      KOFF1   = 1
      KOFF2   = KOFF1 + LOFF1
      KOFF3   = KOFF2 + LOFF2
      KOFF4   = KOFF3 + LOFF3
      KEND1   = KOFF4 + LOFF4 - 1
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_OPTVEC.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_OPTVEC.1',' ',KEND1,LWORK)
C
C==========================================
C     Generate the optimized trial vectors.
C==========================================
C
C-----------------------------------
C     Open files with trial vectors.
C-----------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
C
C-------------------------------------------------------------
C     Transform the 1p1h excitation and de-excitation parts of
C     the trial vectors.
C-------------------------------------------------------------
C
      IOFF = 1 - LTR
C
      DO I1READ = 1,N1READ
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF1),LTR,NTRIAL,LUTR1E,FNTR1E,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF2),LTR,NTRIAL,LUTR1D,FNTR1D,LTR1D,
     &                   IOFF)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),LTR)
C
         CALL SO_WRITESET(WORK(KOFF3),LTR,NTRMX,LUTR1E,FNTR1E,LTR1E,
     &                    IOFF)
         CALL SO_WRITESET(WORK(KOFF4),LTR,NTRMX,LUTR1D,FNTR1D,LTR1D,
     &                    IOFF)
C
      END DO
C
      IF ( L1LTR .GT. 0 ) THEN
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF1),L1LTR,NTRIAL,LUTR1E,FNTR1E,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF2),L1LTR,NTRIAL,LUTR1D,FNTR1D,LTR1D,
     &                   IOFF)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),L1LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),L1LTR)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),L1LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),L1LTR)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),L1LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),L1LTR)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),L1LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),L1LTR)
C
         CALL SO_WRITESET(WORK(KOFF3),L1LTR,NTRMX,LUTR1E,FNTR1E,LTR1E,
     &                    IOFF)
         CALL SO_WRITESET(WORK(KOFF4),L1LTR,NTRMX,LUTR1D,FNTR1D,LTR1D,
     &                    IOFF)
C
      END IF
C
C-------------------------------------------------------------
C     Transform the 2p2h excitation and de-excitation parts of
C     the trial vectors.
C-------------------------------------------------------------
C
      IOFF = 1 - LTR
C
      DO I2READ = 1,N2READ
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF1),LTR,NTRIAL,LUTR2E,FNTR2E,LTR2E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF2),LTR,NTRIAL,LUTR2D,FNTR2D,LTR2D,
     &                   IOFF)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),LTR)
C
         CALL SO_WRITESET(WORK(KOFF3),LTR,NTRMX,LUTR2E,FNTR2E,LTR2E,
     &                    IOFF)
         CALL SO_WRITESET(WORK(KOFF4),LTR,NTRMX,LUTR2D,FNTR2D,LTR2D,
     &                    IOFF)
C
      END DO
C
      IF ( L2LTR .GT. 0 ) THEN
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF1),L2LTR,NTRIAL,LUTR2E,FNTR2E,LTR2E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF2),L2LTR,NTRIAL,LUTR2D,FNTR2D,LTR2D,
     &                   IOFF)
C
         CALL DZERO(WORK(KOFF3),L2LTR*NTRMX)
         CALL DGEMM('N','N',L2LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),L2LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),L2LTR)
C
         CALL DGEMM('N','N',L2LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),L2LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),L2LTR)
C
         CALL DZERO(WORK(KOFF4),L2LTR*NTRMX)
         CALL DGEMM('N','N',L2LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),L2LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),L2LTR)
C
         CALL DGEMM('N','N',L2LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),L2LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),L2LTR)
C
         CALL SO_WRITESET(WORK(KOFF3),L2LTR,NTRMX,LUTR2E,FNTR2E,LTR2E,
     &                    IOFF)
         CALL SO_WRITESET(WORK(KOFF4),L2LTR,NTRMX,LUTR2D,FNTR2D,LTR2D,
     &                    IOFF)
C
      END IF
C
C-----------------------------------------------------------------------
C     For zero frequency linear response, we could have problem if the
C     optimized trial vector will be in the form Z = Y or Z = -Y
C
C     However, if Z = Y or Z = -Y then (Z Y) and (Y Z) are linear
C     dependent and we want instead (Z 0) and (0 Z) as trial vectors.
C     (I try to use as much as I can the name used in RSPORT)         
C-----------------------------------------------------------------------
C
      KTR1E   = 1
      KTR1D   = KTR1E  + LTR1E
      KTR2E   = KTR1D  + LTR1D
      KTR2D   = KTR2E  + LTR2E
      KTTR2   = KTR2D  + LTR2D
      KEND1   = KTTR2  + LTR2E
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_OPTVEC.2',LWORK1)
      IF (LWORK1 .LT. 0) THEN
         WRITE(LUPRI,*) 
     &      'NOTICE: Requirement is due to high print level'
         CALL STOPIT('RP_OPTVEC.2',' ',KEND1,LWORK)
      ENDIF
C
      DO ITR = 1,NTRMX
C
         CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,ITR)
         CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,ITR)
C
         T1=DDOT(LTR1E,WORK(KTR1E),1,WORK(KTR1E),1)+
     &      DDOT(LTR1D,WORK(KTR1D),1,WORK(KTR1D),1)
C         
         OVLPI=TWO*DDOT(LTR1E,WORK(KTR1E),1,WORK(KTR1D),1)
C
         CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,ITR)
         CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,ITR)
C
         CALL DCOPY(LTR2E,WORK(KTR2E),1,WORK(KTTR2),1)
         CALL CCSD_TCMEPKX(WORK(KTTR2),TWO,ISYMTR)
         CALL DSCAL(LTR2E,HALF,WORK(KTTR2),1)
C
         T1 = T1 + DDOT(LTR2E,WORK(KTR2E),1,WORK(KTTR2),1)
         OVLPI = OVLPI+DDOT(LTR2E,WORK(KTR2D),1,WORK(KTTR2),1)
C
         CALL DCOPY(LTR2D,WORK(KTR2D),1,WORK(KTTR2),1)
         CALL CCSD_TCMEPKX(WORK(KTTR2),TWO,ISYMTR)
         CALL DSCAL(LTR2D,HALF,WORK(KTTR2),1)
C
         T1 = T1 + DDOT(LTR2E,WORK(KTR2D),1,WORK(KTTR2),1)
         OVLPI = OVLPI+DDOT(LTR2E,WORK(KTR2E),1,WORK(KTTR2),1)
C
C        test if Z + Y is zero
         TZPY  = T1 + OVLPI
         IF (IPRSOP .GE. 5) WRITE (LUPRI,'(/A/,(A,1P,D15.8))')
     *      ' ** In SO_OPTVEC:  TEST OF Z = Y OR Z = -Y:',
     *      ' (Z Y) norm squared =',T1,
     *      ' Z + Y norm squared =',TZPY
         IF (TZPY .LE. OVLMIN*T1) THEN
            IF (IPRSOP .GE. 5) WRITE (LUPRI,'(/A,I4/A)')
     *         ' Z = -Y in new optimized trial vector no.',ITR,
     *         ' Y component removed.'
            CALL DZERO(WORK(KTR1D),LTR1D)
            CALL DZERO(WORK(KTR2D),LTR2D)
C
            CALL SO_WRITE(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,ITR)
            CALL SO_WRITE(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,ITR)
         ELSE IF (TZPY .GT. T1) THEN
C           test if Z - Y is zero;
C           only if TZPY = 2 * T1 can Z - Y be zero,
C           so TZPY .gt. T1 is a safe test.
            TZMY = T1 - OVLPI
            IF (IPRSOP .GE. 5) WRITE (LUPRI,'(A,1P,D15.8)')
     *         ' Z - Y norm squared =',TZMY
            IF (TZMY .LE. OVLMIN*T1) THEN
               IF (IPRSOP .GE. 5) WRITE (LUPRI,'(/A,I4/A)')
     *         ' Z = Y in new optimized trial vector no.',ITR,
     *         ' Y component removed.'
               CALL DZERO(WORK(KTR1D),LTR1D)
               CALL DZERO(WORK(KTR2D),LTR2D)
C
               CALL SO_WRITE(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,ITR)
               CALL SO_WRITE(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,ITR)
C
            END IF
C
         ENDIF
C
      ENDDO
C
      IF ( IPRSOP .GE. 8 ) THEN
C
C---------------------------------------------------
C        Write new optimized trial vector to output.
C---------------------------------------------------
C
         KTR1E   = 1
         KTR1D   = KTR1E + LTR1E
         KTR2E   = KTR1D + LTR1D
         KTR2D   = KTR2E + LTR2E
         KEND1   = KTR2D + LTR2D - 1
         LWORK1  = LWORK - KEND1
C
         DO ITR = 1,NTRMX
C
            CALL AROUND('Optimized trial vectors')
C
            CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,ITR)
            CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,ITR)
            CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,ITR)
            CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,ITR)
C
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &       (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &       (I,WORK(KTR2E+I-1),WORK(KTR2D+I-1),I=1,LTR2E)
C
         END DO
C
      END IF
C
C------------------------------------
C     Close files with trial vectors.
C------------------------------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
C
C===================================================
C     Generate the E[2] transformed optimized trial.
C===================================================
C
C----------------------------------------------------
C     Open files with E[2] transformed trial vectors.
C----------------------------------------------------
C
      CALL SO_OPEN(LURS1E,FNRS1E,LTR1E)
      CALL SO_OPEN(LURS1D,FNRS1D,LTR1D)
      CALL SO_OPEN(LURS2E,FNRS2E,LTR2E)
      CALL SO_OPEN(LURS2D,FNRS2D,LTR2D)
C
C-------------------------------------------------------------
C     Transform the 1p1h excitation and de-excitation parts of
C     the E[2] transformed trial vectors.
C-------------------------------------------------------------
C
      IOFF = 1 - LTR
C
      DO I1READ = 1,N1READ
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF1),LTR,NTRIAL,LURS1E,FNRS1E,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF2),LTR,NTRIAL,LURS1D,FNRS1D,LTR1D,
     &                   IOFF)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),LTR)
C
         CALL SO_WRITESET(WORK(KOFF3),LTR,NTRMX,LURS1E,FNRS1E,LTR1E,
     &                    IOFF)
         CALL SO_WRITESET(WORK(KOFF4),LTR,NTRMX,LURS1D,FNRS1D,LTR1D,
     &                    IOFF)
C
      END DO
C
      IF ( L1LTR .GT. 0 ) THEN
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF1),L1LTR,NTRIAL,LURS1E,FNRS1E,LTR1E,
     &                   IOFF)
         CALL SO_READSET(WORK(KOFF2),L1LTR,NTRIAL,LURS1D,FNRS1D,LTR1D,
     &                   IOFF)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),L1LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),L1LTR)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),L1LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),L1LTR)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),L1LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),L1LTR)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),L1LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),L1LTR)
C
         CALL SO_WRITESET(WORK(KOFF3),L1LTR,NTRMX,LURS1E,FNRS1E,LTR1E,
     &                    IOFF)
         CALL SO_WRITESET(WORK(KOFF4),L1LTR,NTRMX,LURS1D,FNRS1D,LTR1D,
     &                    IOFF)
C
      END IF
C
C-------------------------------------------------------------
C     Transform the 2p2h excitation and de-excitation parts of
C     the E[2] transformed trial vectors.
C-------------------------------------------------------------
C
      IOFF = 1 - LTR
C
      DO I2READ = 1,N2READ
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF1),LTR,NTRIAL,LURS2E,FNRS2E,LTR2E,
     &                        IOFF)
         CALL SO_READSET(WORK(KOFF2),LTR,NTRIAL,LURS2D,FNRS2D,LTR2D,
     &                        IOFF)
C
         CALL DZERO(WORK(KOFF3),LTR*NTRMX)
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),LTR)
C
         CALL DZERO(WORK(KOFF4),LTR*NTRMX)
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),LTR)
C
         CALL SO_WRITESET(WORK(KOFF3),LTR,NTRMX,LURS2E,FNRS2E,LTR2E,
     &                    IOFF)
         CALL SO_WRITESET(WORK(KOFF4),LTR,NTRMX,LURS2D,FNRS2D,LTR2D,
     &                    IOFF)
C
      END DO
C
      IF ( L2LTR .GT. 0 ) THEN
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF1),L2LTR,NTRIAL,LURS2E,FNRS2E,LTR2E,
     &                        IOFF)
         CALL SO_READSET(WORK(KOFF2),L2LTR,NTRIAL,LURS2D,FNRS2D,LTR2D,
     &                        IOFF)
C
         CALL DZERO(WORK(KOFF3),LTR*NTRMX)
         CALL DGEMM('N','N',L2LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),L2LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),L2LTR)
C
         CALL DGEMM('N','N',L2LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),L2LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),L2LTR)
C
         CALL DZERO(WORK(KOFF4),LTR*NTRMX)
         CALL DGEMM('N','N',L2LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),L2LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),L2LTR)
C
         CALL DGEMM('N','N',L2LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),L2LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),L2LTR)
C
         CALL SO_WRITESET(WORK(KOFF3),L2LTR,NTRMX,LURS2E,FNRS2E,LTR2E,
     &                         IOFF)
         CALL SO_WRITESET(WORK(KOFF4),L2LTR,NTRMX,LURS2D,FNRS2D,LTR2D,
     &                         IOFF)
C
      END IF
C
      IF ( IPRSOP .GE. 8 ) THEN
C
C--------------------------------------------------------------------
C        Write new optimized E[2] transformed trial vector to output.
C--------------------------------------------------------------------
C
         KTR1E   = 1
         KTR1D   = KTR1E  + LTR1E
         KTR2E   = KTR1D  + LTR1D
         KTR2D   = KTR2E  + LTR2E
         KEND1   = KTR2D  + LTR2D - 1
         LWORK1  = LWORK  - KEND1
C
         CALL SO_MEMMAX ('SO_OPTVEC.3',LWORK1)
C
         IF (LWORK1 .LT. 0) THEN
            WRITE(LUPRI,*) 
     &          'NOTICE: Requirement is due to high print level'
            CALL STOPIT('SO_OPTVEC.3',' ',KEND1,LWORK)
         ENDIF
C
         DO ITR = 1,NTRMX
C
            CALL AROUND('Optimized E[2] trial vectors')
C
            CALL SO_READ(WORK(KTR1E),LTR1E,LURS1E,FNRS1E,ITR)
            CALL SO_READ(WORK(KTR1D),LTR1D,LURS1D,FNRS1D,ITR)
            CALL SO_READ(WORK(KTR2E),LTR2E,LURS2E,FNRS2E,ITR)
            CALL SO_READ(WORK(KTR2D),LTR2D,LURS2D,FNRS2D,ITR)
C
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &         (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &         (I,WORK(KTR2E+I-1),WORK(KTR2D+I-1),I=1,LTR2E)
C
         END DO
C
      END IF
C
C-----------------------------------------------------
C     Close files with E[2] transformed trial vectors.
C-----------------------------------------------------
C
      CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
      CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
      CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
      CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
C
C===================================================
C     Generate the S[2] transformed optimized trial.
C===================================================
C
C----------------------------------------------------
C     Open files with S[2] transformed trial vectors.
C----------------------------------------------------
C
      CALL SO_OPEN(LURO1E,FNRO1E,LTR1E)
      CALL SO_OPEN(LURO1D,FNRO1D,LTR1D)
C
C-------------------------------------------------------------
C     Transform the 1p1h excitation and de-excitation parts of
C     the S[2] transformed trial vectors.
C-------------------------------------------------------------
C
      IOFF = 1 - LTR
C
      DO I1READ = 1,N1READ
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF1),LTR,NTRIAL,LURO1E,FNRO1E,LTR1E,
     &                        IOFF)
         CALL SO_READSET(WORK(KOFF2),LTR,NTRIAL,LURO1D,FNRO1D,LTR1D,
     &                        IOFF)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,-ONE,WORK(KOFF2),LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),LTR)
C
         CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,-ONE,WORK(KOFF1),LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),LTR)
C
         CALL SO_WRITESET(WORK(KOFF3),LTR,NTRMX,LURO1E,FNRO1E,LTR1E,
     &                         IOFF)
         CALL SO_WRITESET(WORK(KOFF4),LTR,NTRMX,LURO1D,FNRO1D,LTR1D,
     &                         IOFF)
C
      END DO
C
      IF ( L1LTR .GT. 0 ) THEN
C
         IOFF = IOFF + LTR
C
         CALL SO_READSET(WORK(KOFF1),L1LTR,NTRIAL,LURO1E,FNRO1E,LTR1E,
     &                        IOFF)
         CALL SO_READSET(WORK(KOFF2),L1LTR,NTRIAL,LURO1D,FNRO1D,LTR1D,
     &                        IOFF)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),L1LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),L1LTR)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,-ONE,WORK(KOFF2),L1LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),L1LTR)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),L1LTR,
     &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),L1LTR)
C
         CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,-ONE,WORK(KOFF1),L1LTR,
     &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),L1LTR)
C
         CALL SO_WRITESET(WORK(KOFF3),L1LTR,NTRMX,LURO1E,FNRO1E,LTR1E,
     &                    IOFF)
         CALL SO_WRITESET(WORK(KOFF4),L1LTR,NTRMX,LURO1D,FNRO1D,LTR1D,
     &                    IOFF)
C
      END IF
C
      IF ( IPRSOP .GE. 8 ) THEN
C
C--------------------------------------------------------------------
C        Write new optimized S[2] transformed trial vector to output.
C--------------------------------------------------------------------
C
         KTR1E   = 1
         KTR1D   = KTR1E  + LTR1E
         KEND1   = KTR1D  + LTR1D - 1
         LWORK1  = LWORK  - KEND1
C
         CALL SO_MEMMAX ('SO_OPTVEC.4',LWORK1)
C
         IF (LWORK1 .LT. 0) THEN
            WRITE(LUPRI,*) 
     &          'NOTICE: Requirement is due to high print level'
            CALL STOPIT('SO_OPTVEC.4',' ',KEND1,LWORK)
         ENDIF
C
         DO ITR = 1,NTRMX
C
            CALL SO_READ(WORK(KTR1E),LTR1E,LURO1E,FNRO1E,ITR)
            CALL SO_READ(WORK(KTR1D),LTR1D,LURO1D,FNRO1D,ITR)
C
            CALL AROUND('Optimized S[2] trial vectors')
C
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &          (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
C
         END DO
C
      END IF
C
C------------------------------------
C     Close files with trial vectors.
C------------------------------------
C
      CALL SO_CLOSE(LURO1E,FNRO1E,'KEEP')
      CALL SO_CLOSE(LURO1D,FNRO1D,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_OPTVEC')
C
      RETURN
C
      END
