C
C  /* Deck so_orth_trn */
      SUBROUTINE SO_ORTH_TRN(DOUBLES,LTYPE,NOLDTR,NNEWTR,NLINDP,ISYMTR,
     &                       DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                       WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Orthogonalize new trial vector against all previous
C              trial vectors (including the paired ones) and
C              normalize. Finally make a symmetric orthonormalization
C              of the the new trial vector and its pair trial vector.
C
#include "implicit.h"
#include "priunit.h"
C
#include "ccsdsym.h"
#include "ccorb.h"
#include "soppinf.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRLDP = 1.0D-20, THROUND = 1.0D-4, OVLMIN = 1.0D-10)
      PARAMETER (T1MIN = 1.0D-20)
C
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION WORK(LWORK)
C
      LOGICAL   DOUBLES
C
      CHARACTER*6 LTYPE
C     DOUBLES = .TRUE.
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_ORTH_TRN')
C
C-----------------------------------------------------------------
C     Initialize number of linear dependent trial vectors to zero.
C-----------------------------------------------------------------
C
      NLINDP = 0
      ISYRES = MULD2H(ISYMOP,ISYMTR)
C
C------------------------------------------------
C     Allocation of work space for trial vectors.
C------------------------------------------------
C
      LTR1E   = NT1AM(ISYMTR)
      LTR1D   = NT1AM(ISYMTR)
      IF (DOUBLES) THEN
         LTR2E   = N2P2HOP(ISYMTR)
         LTR2D   = N2P2HOP(ISYMTR)
      ELSE
         LTR2E   = 0
         LTR2D   = 0
      ENDIF
      LRSO1E  = NT1AM(ISYMTR)
      LRSO1D  = NT1AM(ISYMTR)
C
      K1TR1E  = 1
      K1TR1D  = K1TR1E + LTR1E
      K1TR2E  = K1TR1D + LTR1D
      K1TR2D  = K1TR2E + LTR2E
      KTTR2   = K1TR2D + LTR2D
      KTR1E   = KTTR2  + LTR2D
      KTR1D   = KTR1E  + LTR1E
      KTR2E   = KTR1D  + LTR1D
      KTR2D   = KTR2E  + LTR2E
      KRSO1E  = KTR2D  + LTR2D
      KRSO1D  = KRSO1E + LRSO1E
      KEND1   = KRSO1D + LRSO1D
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_ORTH_TRN',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_ORTH_TRN',' ',KEND1,LWORK)
C
C-------------------------------------
C     Loop over new raw trial vectors.
C-------------------------------------
C
      DO 100 INEWTR = 1,NNEWTR
C
C----------------------------------
C        Read new raw trial vector.
C----------------------------------
C
         CALL SO_READ(WORK(K1TR1E),LTR1E,LUTR1E,FNTR1E,NOLDTR+INEWTR)
         CALL SO_READ(WORK(K1TR1D),LTR1D,LUTR1D,FNTR1D,NOLDTR+INEWTR)
         IF (DOUBLES) THEN
            CALL SO_READ(WORK(K1TR2E),LTR2E,LUTR2E,FNTR2E,NOLDTR+INEWTR)
            CALL SO_READ(WORK(K1TR2D),LTR2D,LUTR2D,FNTR2D,NOLDTR+INEWTR)
         ENDIF
C
         IF ( IPRSOP .GE. 7 ) THEN
C
            CALL AROUND('Raw new trial vector in SO_ORTH_TRN')
C
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &       (I,WORK(K1TR1E+I-1),WORK(K1TR1D+I-1),I=1,LTR1E)
            IF (DOUBLES)
     &             WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &              (I,WORK(K1TR2E+I-1),WORK(K1TR2D+I-1),I=1,LTR2E)
C
         END IF
C
         ITURN = 0
C
  200    CONTINUE
C
         ITURN = ITURN + 1
C
C-----------------------------------------
C        Loop over previous trial vectors.
C-----------------------------------------
C
         DO 300 IPRVTR = 1,NOLDTR+(INEWTR-NLINDP)-1
C
C---------------------------------------
C           Read previous trial vectors.
C---------------------------------------
C
            CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,IPRVTR)
            CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,IPRVTR)
            IF (DOUBLES) THEN
               CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,IPRVTR)
               CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,IPRVTR)
            ENDIF
C
            IF ( (INEWTR.EQ.1) .AND. (IPRSOP .GE. 7) ) THEN
C
               CALL AROUND('Previous trial vector in SO_ORTH_TRN')
C
               WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &              (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
               IF (DOUBLES)
     &              WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &              (I,WORK(KTR2E+I-1),WORK(KTR2D+I-1),I=1,LTR2E)
C
            END IF
C
C----------------------------------------------------------------
C           Orthogonalize new trial vector against previous trial
C           vectors and their paired partners.
C           KTR1* and KTR2*  : pointers to previous trial vectors
C           K1TR1* and K1TR2*: pointers to raw new trial vectors
C----------------------------------------------------------------
C
            CALL SO_RES_O(WORK(KRSO1E),LRSO1E, WORK(KRSO1D),LRSO1D,
     &                    WORK(KTR1E), LTR1E,  WORK(KTR1D), LTR1D,
     &                    DENSIJ,      LDENSIJ,DENSAB,      LDENSAB,
     &                    ISYRES,      ISYMTR)
C
            DOTP = DDOT(LTR1E,WORK(K1TR1E),1,WORK(KRSO1E),1)
     &           + DDOT(LTR1D,WORK(K1TR1D),1,WORK(KRSO1D),1)
C
            IF (DOUBLES) THEN
C
               IF (TRIPLET) THEN
C
                  DOTP = DOTP + DDOT(LTR2E,WORK(K1TR2E),1,WORK(KTR2E),1)
                  DOTP = DOTP - DDOT(LTR2D,WORK(K1TR2D),1,WORK(KTR2D),1)
C
               ELSE
C
                  DOTP = DOTP + SO_DOUBLES_PRODUCT(WORK(K1TR2E),
     &                                             WORK(KTR2E),ISYMTR)
C
                  DOTP = DOTP - SO_DOUBLES_PRODUCT(WORK(K1TR2D),
     &                                             WORK(KTR2D),ISYMTR)
C
               END IF
C
               CALL DAXPY(LTR2E,-DOTP,WORK(KTR2E),1,WORK(K1TR2E),1)
               CALL DAXPY(LTR2D,-DOTP,WORK(KTR2D),1,WORK(K1TR2D),1)
C
            END IF
C
            CALL DAXPY(LTR1E,-DOTP,WORK(KTR1E),1,WORK(K1TR1E),1)
            CALL DAXPY(LTR1D,-DOTP,WORK(KTR1D),1,WORK(K1TR1D),1)
C
            CALL SO_RES_O(WORK(KRSO1D),LRSO1D, WORK(KRSO1E),LRSO1E,
     &                    WORK(KTR1D), LTR1D,  WORK(KTR1E), LTR1E,
     &                    DENSIJ,      LDENSIJ,DENSAB,      LDENSAB,
     &                    ISYRES,      ISYMTR)
C
            DOTP = DDOT(LTR1E,WORK(K1TR1E),1,WORK(KRSO1D),1)
     &           + DDOT(LTR1D,WORK(K1TR1D),1,WORK(KRSO1E),1)
C
            IF (DOUBLES) THEN
C
               IF (TRIPLET) THEN
C
                  DOTP = DOTP + DDOT(LTR2E,WORK(K1TR2E),1,WORK(KTR2D),1)
                  DOTP = DOTP - DDOT(LTR2D,WORK(K1TR2D),1,WORK(KTR2E),1)
C
               ELSE
C
                  DOTP = DOTP + SO_DOUBLES_PRODUCT(WORK(K1TR2E),
     &                                             WORK(KTR2D),ISYMTR)
C
                  DOTP = DOTP - SO_DOUBLES_PRODUCT(WORK(K1TR2D),
     &                                             WORK(KTR2E),ISYMTR)
C
               END IF
C
               CALL DAXPY(LTR2E,DOTP,WORK(KTR2D),1,WORK(K1TR2E),1)
               CALL DAXPY(LTR2D,DOTP,WORK(KTR2E),1,WORK(K1TR2D),1)
C
            END IF
C
            CALL DAXPY(LTR1E,DOTP,WORK(KTR1D),1,WORK(K1TR1E),1)
            CALL DAXPY(LTR1D,DOTP,WORK(KTR1E),1,WORK(K1TR1D),1)
C
  300    CONTINUE
C
C----------------------------------------------------
C        Calculate absolute norm of new trial vector.
C----------------------------------------------------
C
         CALL SO_RES_O(WORK(KRSO1E),LRSO1E, WORK(KRSO1D),LRSO1D,
     &                 WORK(K1TR1E),LTR1E,  WORK(K1TR1D),LTR1D,
     &                 DENSIJ,      LDENSIJ,DENSAB,      LDENSAB,
     &                 ISYRES,      ISYMTR)
C
         DNORM = DDOT(LTR1E,WORK(K1TR1E),1,WORK(KRSO1E),1)
     &         + DDOT(LTR1D,WORK(K1TR1D),1,WORK(KRSO1D),1)
C
         IF (DOUBLES) THEN
            IF (TRIPLET) THEN
               DNORM = DNORM + DDOT(LTR2E,WORK(K1TR2E),1,WORK(K1TR2E),1)
               DNORM = DNORM - DDOT(LTR2D,WORK(K1TR2D),1,WORK(K1TR2D),1)
            ELSE
               DNORM = DNORM + so_doubles_length(work(k1tr2e),isymtr)
               DNORM = DNORM - so_doubles_length(work(k1tr2d),isymtr)
            END IF
         ENDIF
C
C----------------------------------------------------------
C        Remove new trial vector if it is linear dependent.
C----------------------------------------------------------
C
         THNORM = THRLDP
C
         IF (DABS(DNORM) .LE. THNORM) THEN
C
            IF ( IPRSOP .GE. 1 ) WRITE(LUPRI,9002) DNORM
C
            NLINDP = NLINDP + 1
C
            GO TO 100
C
         END IF
C
         IF ( DNORM .LT. ZERO ) THEN
C
C-------------------------------------------------
C           Switch X and Y part of reduced vector.
C-------------------------------------------------
C
            DO ITR1E = 1,LTR1E
               TEMP                 = WORK(K1TR1E+ITR1E-1)
               WORK(K1TR1E+ITR1E-1) = WORK(K1TR1D+ITR1E-1)
               WORK(K1TR1D+ITR1E-1) = TEMP
            END DO
C
            DO ITR2E = 1,LTR2E
               TEMP                 = WORK(K1TR2E+ITR2E-1)
               WORK(K1TR2E+ITR2E-1) = WORK(K1TR2D+ITR2E-1)
               WORK(K1TR2D+ITR2E-1) = TEMP
            END DO
C
         END IF
C
C----------------------------------------------------
C        If Norm is little normalize the trial vector
C        a first time
C----------------------------------------------------
C
         IF((DABS(DNORM).LE.T1MIN).AND.(LTYPE.EQ.'LINEAR')) THEN
C
           IF(IPRSOP.GT.10) WRITE(LUPRI,'(1X,A,D15.8)')
     &       'I make a first normalization with norm',DNORM
C
           DNORMI = ONE  / DSQRT(DABS(DNORM))
C           print *, dnorm
CRF is the correct, don't we need to recalculate or scale KRSO1*?
           CALL DSCAL(LTR1E,DNORMI,WORK(K1TR1E),1)
           CALL DSCAL(LTR1D,DNORMI,WORK(K1TR1D),1)
C
           DNORM = DDOT(LTR1E,WORK(K1TR1E),1,WORK(KRSO1E),1)
     &           + DDOT(LTR1D,WORK(K1TR1D),1,WORK(KRSO1D),1)
C
           IF (DOUBLES) THEN
              CALL DSCAL(LTR2E,DNORMI,WORK(K1TR2E),1)
              CALL DSCAL(LTR2D,DNORMI,WORK(K1TR2D),1)
C
              IF (TRIPLET) THEN
                 DNORM = DNORM+DDOT(LTR2E,WORK(K1TR2E),1,WORK(K1TR2E),1)
                 DNORM = DNORM-DDOT(LTR2D,WORK(K1TR2D),1,WORK(K1TR2D),1)
              ELSE
                 CALL DCOPY(LTR2E,WORK(K1TR2E),1,WORK(KTTR2),1)
                 CALL CCSD_TCMEPKX(WORK(KTTR2),TWO,ISYMTR)
C
                 DNORM = DNORM +
     &                   HALF * DDOT(LTR2E,WORK(KTTR2),1,WORK(K1TR2E),1)
C
                 CALL DCOPY(LTR2E,WORK(K1TR2D),1,WORK(KTTR2),1)
                 CALL CCSD_TCMEPKX(WORK(KTTR2),TWO,ISYMTR)
C
                 DNORM = DNORM -
     &                   HALF * DDOT(LTR2D,WORK(KTTR2),1,WORK(K1TR2D),1)
              END IF
           ENDIF
C
         ENDIF
C
C-------------------------------------------------------
C        Normalize new trial vector first (second) time.
C-------------------------------------------------------
C
C         print *, dnorm
         DNORMI = ONE / DSQRT( DABS(DNORM) )
C
         CALL DSCAL(LTR1E,DNORMI,WORK(K1TR1E),1)
         CALL DSCAL(LTR1D,DNORMI,WORK(K1TR1D),1)
         IF (DOUBLES) THEN
            CALL DSCAL(LTR2E,DNORMI,WORK(K1TR2E),1)
            CALL DSCAL(LTR2D,DNORMI,WORK(K1TR2D),1)
         ENDIF
C
C--------------------------------------------------------------------
C        In case the norm of the orthogonalized new trial vector is
C        less than THROUND, the orthogonalization is repeated once by
C        looping back to line 200.
C--------------------------------------------------------------------
C
         IF (DABS(DNORM) .LT. THROUND) THEN
C
            IF (ITURN .LE. 2) THEN
               GO TO 200
            ELSE
               WRITE(LUPRI,9003)
            END IF
C
         END IF
C
C-----------------------------------------------------------------
C        Write the new orthogonalized trial vector to file (and to
C        output).
C-----------------------------------------------------------------
C
         CALL SO_WRITE(WORK(K1TR1E),LTR1E,LUTR1E,FNTR1E,
     &                 NOLDTR+INEWTR-NLINDP)
         CALL SO_WRITE(WORK(K1TR1D),LTR1D,LUTR1D,FNTR1D,
     &                 NOLDTR+INEWTR-NLINDP)
         IF (DOUBLES) THEN
            CALL SO_WRITE(WORK(K1TR2E),LTR2E,LUTR2E,FNTR2E,
     &                    NOLDTR+INEWTR-NLINDP)
            CALL SO_WRITE(WORK(K1TR2D),LTR2D,LUTR2D,FNTR2D,
     &                    NOLDTR+INEWTR-NLINDP)
         ENDIF
C
         IF ( IPRSOP .GE. 6 ) THEN
C
            CALL AROUND('Orthonormalized trialvector in SO_ORTH_TRN')
C
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &        (I,WORK(K1TR1E+I-1),WORK(K1TR1D+I-1),I=1,LTR1E)
            IF (DOUBLES) THEN
                WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &            (I,WORK(K1TR2E+I-1),WORK(K1TR2D+I-1),I=1,LTR2E)
            ENDIF
C
         END IF
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_ORTH_TRN')
C
      RETURN
C
 9001 FORMAT(/,'Square NORM of new trial vector: ',F22.18,/)
 9002 FORMAT(/,'Norm of normalized new trial vector is: ',F22.18,/
     &       'New trial vector is removed because of linear ',
     &       'dependence.')
 9003 FORMAT(/,'WARNING: Problems orthonormalizing in SO_ORTH_TRN')
C
      END
