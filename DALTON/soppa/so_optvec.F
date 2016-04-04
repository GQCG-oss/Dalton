C  /* Deck so_optvec */
      SUBROUTINE SO_OPTVEC(DOUBLES,NOLDTR,NNEWTR,NEXCI,
     &                     EIVEC,LEIVEC,ISYMTR,
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
#include "implicit.h"
#include "priunit.h"
#include "soppinf.h"
#include "ccsdsym.h"
C
      PARAMETER   (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER   (D100 = 100.0D0)
      PARAMETER   (THRNM = 1.0D-15, OVLMIN = 10D-20 ) ! Try (from rp) 1.0D-4) 
      CHARACTER*3 YES, NO
      LOGICAL    DOUBLES
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
      IF(DOUBLES)THEN
CPi 01.04.16
C         LTR2E  = NT2AM(ISYMTR)
C         LTR2D  = NT2AM(ISYMTR)
         LTR2E  = N2P2HOP(ISYMTR)
         LTR2D  = N2P2HOP(ISYMTR)
Cend-Pi
      ELSE
         LTR2E  = 0
         LTR2D  = 0
      ENDIF
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
      IF(DOUBLES)THEN
         N2READ = LTR2E / LTR
C
         L2LTR  = LTR2E - ( N2READ * LTR )
      ENDIF
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      IF(DOUBLES)THEN
         LOFF1  = MIN((LTR * NTRIAL),(LTR2E * NTRIAL))
         LOFF2  = MIN((LTR * NTRIAL),(LTR2E * NTRIAL))
         LOFF3  = MIN((LTR * NTRMX), (LTR2E * NTRMX))
         LOFF4  = MIN((LTR * NTRMX), (LTR2E * NTRMX))
      ELSE
         LOFF1  = MIN((LTR * NTRIAL),(LTR1E * NTRIAL))
         LOFF2  = MIN((LTR * NTRIAL),(LTR1E * NTRIAL))
         LOFF3  = MIN((LTR * NTRMX), (LTR1E * NTRMX))
         LOFF4  = MIN((LTR * NTRMX), (LTR1E * NTRMX))
      ENDIF
C
      KOFF1   = 1
      KOFF2   = KOFF1 + LOFF1
      KOFF3   = KOFF2 + LOFF2
      KOFF4   = KOFF3 + LOFF3
CPi 04.04.16 old triplet version
C      KEND1   = KOFF4 + LOFF4 - 1
      KEND1   = KOFF4 + LOFF4
Cend-Pi
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
      IF(DOUBLES)THEN
         CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
         CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
      ENDIF
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
C     Only for double excitations
      IF(DOUBLES)THEN
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
C     All methods continue at this point
      ENDIF
C
C-----------------------------------------------------------------------
C     For zero frequency linear response, we could have problem because the
C     optimized trial vector will be in the form Z = Y or Z = -Y
C
C     However, if Z = Y or Z = -Y then (Z Y) and (Y Z) are linear
C     dependent and we want instead (Z 0) and (0 Z) as trial vectors.
C     (I try to use as much as I can the name used in RSPORT)         
C-----------------------------------------------------------------------
C
C      KTR1E   = 1
C      KTR1D   = KTR1E  + LTR1E
C      KTR2E   = KTR1D  + LTR1D
C      KTR2D   = KTR2E  + LTR2E
C      KTTR2   = KTR2D  + LTR2D
C      KEND1   = KTTR2  + LTR2E
C      LWORK1  = LWORK  - KEND1
CC
C      CALL SO_MEMMAX ('SO_OPTVEC.2',LWORK1)
C      IF (LWORK1 .LT. 0) THEN
C         WRITE(LUPRI,*) 
C     &      'NOTICE: Requirement is due to high print level'
C         CALL STOPIT('SO_OPTVEC.2',' ',KEND1,LWORK)
C      ENDIF
CC
C      DO ITR = 1,NTRMX
CC
C         CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,ITR)
C         CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,ITR)
CC
C         T1=DDOT(LTR1E,WORK(KTR1E),1,WORK(KTR1E),1)+
C     &      DDOT(LTR1D,WORK(KTR1D),1,WORK(KTR1D),1)
CC         
C         OVLPI=TWO*DDOT(LTR1E,WORK(KTR1E),1,WORK(KTR1D),1)
C
C         IF(DOUBLES)THEN
CC
C            CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,ITR)
C            CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,ITR)
CC
C            CALL DCOPY(LTR2E,WORK(KTR2E),1,WORK(KTTR2),1)
C            CALL CCSD_TCMEPKX(WORK(KTTR2),TWO,ISYMTR)
C            CALL DSCAL(LTR2E,HALF,WORK(KTTR2),1)
CC
C            T1 = T1 + DDOT(LTR2E,WORK(KTR2E),1,WORK(KTTR2),1)
C            OVLPI = OVLPI+DDOT(LTR2E,WORK(KTR2D),1,WORK(KTTR2),1)
CC
C            CALL DCOPY(LTR2D,WORK(KTR2D),1,WORK(KTTR2),1)
C            CALL CCSD_TCMEPKX(WORK(KTTR2),TWO,ISYMTR)
C            CALL DSCAL(LTR2D,HALF,WORK(KTTR2),1)
CC
C            T1 = T1 + DDOT(LTR2E,WORK(KTR2D),1,WORK(KTTR2),1)
C            OVLPI = OVLPI+DDOT(LTR2E,WORK(KTR2E),1,WORK(KTTR2),1)
CC
C         ENDIF
CC        test if Z + Y is zero
C         TZPY  = T1 + OVLPI
C         IF (IPRSOP .GE. 5) WRITE (LUPRI,'(/A/,(A,1P,D15.8))')
C     *      ' ** In SO_OPTVEC:  TEST OF Z = Y OR Z = -Y:',
C     *      ' (Z Y) norm squared =',T1,
C     *      ' Z + Y norm squared =',TZPY
C         IF (TZPY .LE. OVLMIN*T1) THEN
C            IF (IPRSOP .GE. 5) WRITE (LUPRI,'(/A,I4/A)')
C     *         ' Z = -Y in new optimized trial vector no.',ITR,
C     *         ' Y component removed.'
C            CALL DZERO(WORK(KTR1D),LTR1D)
C            IF(DOUBLES) CALL DZERO(WORK(KTR2D),LTR2D)
CC
C            CALL SO_WRITE(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,ITR)
C            IF(DOUBLES) CALL SO_WRITE(WORK(KTR2D),LTR2D,LUTR2D,
C     &                                FNTR2D,ITR)
C         ELSE IF (TZPY .GT. T1) THEN
CC           test if Z - Y is zero;
CC           only if TZPY = 2 * T1 can Z - Y be zero,
CC           so TZPY .gt. T1 is a safe test.
C            TZMY = T1 - OVLPI
C            IF (IPRSOP .GE. 5) WRITE (LUPRI,'(A,1P,D15.8)')
C     *         ' Z - Y norm squared =',TZMY
C            IF (TZMY .LE. OVLMIN*T1) THEN
C               IF (IPRSOP .GE. 5) WRITE (LUPRI,'(/A,I4/A)')
C     *         ' Z = Y in new optimized trial vector no.',ITR,
C     *         ' Y component removed.'
C               CALL DZERO(WORK(KTR1D),LTR1D)
C               IF(DOUBLES) CALL DZERO(WORK(KTR2D),LTR2D)
CC
C               CALL SO_WRITE(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,ITR)
C               IF(DOUBLES) CALL SO_WRITE(WORK(KTR2D),LTR2D,LUTR2D,
C     &                                   FNTR2D,ITR)
CC
C            END IF
CC
C         ENDIF
CC
C      ENDDO
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
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &       (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
C
            IF(DOUBLES)THEN
               CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,ITR)
               CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,ITR)
               WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &         (I,WORK(KTR2E+I-1),WORK(KTR2D+I-1),I=1,LTR2E)
            ENDIF
C
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
      IF(DOUBLES)THEN
         CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
         CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
      ENDIF
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
      IF(DOUBLES)THEN
         CALL SO_OPEN(LURS2E,FNRS2E,LTR2E)
         CALL SO_OPEN(LURS2D,FNRS2D,LTR2D)
      ENDIF
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
C     Only Doubles methods do this block
      IF(DOUBLES)THEN
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
C     All methods continue here
      ENDIF
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
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &         (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
            IF(DOUBLES)THEN
               CALL SO_READ(WORK(KTR2E),LTR2E,LURS2E,FNRS2E,ITR)
               CALL SO_READ(WORK(KTR2D),LTR2D,LURS2D,FNRS2D,ITR)
C
               WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &            (I,WORK(KTR2E+I-1),WORK(KTR2D+I-1),I=1,LTR2E)
            ENDIF
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
      IF(DOUBLES)THEN
         CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
         CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
      ENDIF
C
C===================================================
C     Generate the S[2] transformed optimized trial.
C===================================================
C
C----------------------------------------------------
C     Open files with S[2] transformed trial vectors.
C----------------------------------------------------
CRF   I've stopped saving this on disk, calculate on the
C     fly if ever needed
C      
CRF     CALL SO_OPEN(LURO1E,FNRO1E,LTR1E)
CRF     CALL SO_OPEN(LURO1D,FNRO1D,LTR1D)
C
C-------------------------------------------------------------
C     Transform the 1p1h excitation and de-excitation parts of
C     the S[2] transformed trial vectors.
C-------------------------------------------------------------
C
CRF     IOFF = 1 - LTR
C
CRF     DO I1READ = 1,N1READ
C
CRF        IOFF = IOFF + LTR
C
CRF        CALL SO_READSET(WORK(KOFF1),LTR,NTRIAL,LURO1E,FNRO1E,LTR1E,
CRF    &                        IOFF)
CRF        CALL SO_READSET(WORK(KOFF2),LTR,NTRIAL,LURO1D,FNRO1D,LTR1D,
CRF    &                        IOFF)
C
CRF        CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),LTR,
CRF    &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),LTR)
C
CRF        CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,-ONE,WORK(KOFF2),LTR,
CRF    &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),LTR)
C
CRF        CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),LTR,
CRF    &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),LTR)
C
CRF        CALL DGEMM('N','N',LTR,NTRMX,NTRIAL,-ONE,WORK(KOFF1),LTR,
CRF    &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),LTR)
C
CRF        CALL SO_WRITESET(WORK(KOFF3),LTR,NTRMX,LURO1E,FNRO1E,LTR1E,
CRF    &                         IOFF)
CRF        CALL SO_WRITESET(WORK(KOFF4),LTR,NTRMX,LURO1D,FNRO1D,LTR1D,
CRF    &                         IOFF)
C
CRF     END DO
C
CRF     IF ( L1LTR .GT. 0 ) THEN
C
CRF        IOFF = IOFF + LTR
C
CRF        CALL SO_READSET(WORK(KOFF1),L1LTR,NTRIAL,LURO1E,FNRO1E,LTR1E,
CRF    &                        IOFF)
CRF        CALL SO_READSET(WORK(KOFF2),L1LTR,NTRIAL,LURO1D,FNRO1D,LTR1D,
CRF    &                        IOFF)
C
CRF        CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF1),L1LTR,
CRF    &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF3),L1LTR)
C
CRF        CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,-ONE,WORK(KOFF2),L1LTR,
CRF    &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF3),L1LTR)
C
CRF        CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,ONE,WORK(KOFF2),L1LTR,
CRF    &              EIVEC(1,1),2*NTRIAL,ZERO,WORK(KOFF4),L1LTR)
C
CRF        CALL DGEMM('N','N',L1LTR,NTRMX,NTRIAL,-ONE,WORK(KOFF1),L1LTR,
CRF    &              EIVEC(1+NTRIAL,1),2*NTRIAL,ONE,WORK(KOFF4),L1LTR)
C
CRF        CALL SO_WRITESET(WORK(KOFF3),L1LTR,NTRMX,LURO1E,FNRO1E,LTR1E,
CRF    &                    IOFF)
CRF        CALL SO_WRITESET(WORK(KOFF4),L1LTR,NTRMX,LURO1D,FNRO1D,LTR1D,
CRF    &                    IOFF)
C
CRF     END IF
C
CRF     IF ( IPRSOP .GE. 8 ) THEN
C
C--------------------------------------------------------------------
C        Write new optimized S[2] transformed trial vector to output.
C--------------------------------------------------------------------
C
CRF        KTR1E   = 1
CRF        KTR1D   = KTR1E  + LTR1E
CRF        KEND1   = KTR1D  + LTR1D - 1
CRF        LWORK1  = LWORK  - KEND1
C
CRF        CALL SO_MEMMAX ('SO_OPTVEC.4',LWORK1)
C
CRF        IF (LWORK1 .LT. 0) THEN
CRF           WRITE(LUPRI,*) 
CRF    &          'NOTICE: Requirement is due to high print level'
CRF           CALL STOPIT('SO_OPTVEC.4',' ',KEND1,LWORK)
CRF        ENDIF
C
CRF        DO ITR = 1,NTRMX
C
CRF           CALL SO_READ(WORK(KTR1E),LTR1E,LURO1E,FNRO1E,ITR)
CRF           CALL SO_READ(WORK(KTR1D),LTR1D,LURO1D,FNRO1D,ITR)
C
CRF           CALL AROUND('Optimized S[2] trial vectors')
C
CRF           WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
CRF    &          (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
C
CRF        END DO
C
CRF     END IF
C
C------------------------------------
C     Close files with trial vectors.
C------------------------------------
C
CRF     CALL SO_CLOSE(LURO1E,FNRO1E,'KEEP')
CRF     CALL SO_CLOSE(LURO1D,FNRO1D,'KEEP')
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
