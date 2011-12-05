C
C  /* Deck rp_orth_trn */
      SUBROUTINE RP_ORTH_TRN(LTYPE,NOLDTR,NNEWTR,NLINDP,ISYMTR,
     &                       WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Andrea Ligabue, January 2004: added the check for Z = +- Y
C       when freq=0 in linear response properties
C
C     PURPOSE: Orthogonalize new trial vector against all previous
C              trial vectors (including the paired ones) and 
C              normalize. Finally make a symmetric orthonormalization 
C              of the the new trial vector and its pair trial vector. 
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccsdsym.h>
#include <ccorb.h>
#include <soppinf.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRLDP = 1.0D-20, THROUND = 1.0D-4)
      PARAMETER (T1MIN = 1.0D-8, OVLMIN = 1.0D-20)
C
      CHARACTER*6 LTYPE
      DIMENSION WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_ORTH_TRN')
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
      LRSO1E  = NT1AM(ISYMTR)
      LRSO1D  = NT1AM(ISYMTR)
C
      K1TR1E  = 1
      K1TR1D  = K1TR1E + LTR1E
      KTR1E   = K1TR1D + LTR1D
      KTR1D   = KTR1E  + LTR1E
      KRSO1E  = KTR1D  + LTR1D
      KRSO1D  = KRSO1E + LRSO1E
      KEND1   = KRSO1D + LRSO1D
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_ORTH_TRN',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_ORTH_TRN',' ',KEND1,LWORK)
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
C
         IF ( IPRSOP .GE. 7 ) THEN
C
            CALL AROUND('Raw new trial vector in RP_ORTH_TRN')
C
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(K1TR1E+I-1),WORK(K1TR1D+I-1),I=1,LTR1E)
C
         END IF
C
         ITURN = 0
C
         IF (LTYPE .EQ. 'LINEAR') THEN
C
C-----------------------------------------------------------------------
C        For zero frequency linear response, the result will
C        be Z = Y for real perturbations and Z = -Y for imaginary
C        perturbations; trial vectors will have the same structure.
C
C
C        However, if Z = Y or Z = -Y then (Z Y) and (Y Z) are linear
C        dependent and we want instead (Z 0) and (0 Z) as trial vectors.
C        (I try to use as much as I can the name used in RSPORT)         
C-----------------------------------------------------------------------
C
         T1=DDOT(LTR1E,WORK(K1TR1E),1,WORK(K1TR1E),1)+
     &      DDOT(LTR1D,WORK(K1TR1D),1,WORK(K1TR1D),1) 
         OVLPI=TWO*DDOT(LTR1E,WORK(K1TR1E),1,WORK(K1TR1D),1)
C        test if Z + Y is zero
         TZPY  = T1 + OVLPI
         IF (IPRSOP .GE. 5) WRITE (LUPRI,'(/A/,(A,1P,D15.8))')
     *      ' ** in RP_ORTH_TRN: TEST OF Z = Y OR Z = -Y:',
     *      ' (Z Y) norm squared =',T1,
     *      ' Z + Y norm squared =',TZPY
         IF (TZPY .LE. OVLMIN*T1) THEN
            IF (IPRSOP .GE. 5) WRITE (LUPRI,'(/A,I4/A)')
     *         ' Z = -Y in new trial vector no.',INEWTR,
     *         ' Y component removed.'
            CALL DZERO(WORK(K1TR1D),LTR1D)
         ELSE IF (TZPY .GT. T1) THEN
C           test if Z - Y is zero;
C           only if TZPY = 2 * T1 can Z - Y be zero,
C           so TZPY .gt. T1 is a safe test.
            TZMY = T1 - OVLPI
            IF (IPRSOP .GE. 5) WRITE (LUPRI,'(A,1P,D15.8)')
     *         ' Z - Y norm squared =',TZMY
            IF (TZMY .LE. OVLMIN*T1) THEN
               IF (IPRSOP .GE. 5) WRITE (LUPRI,'(/A,I4/A)')
     *         ' Z = Y in new trial vector no.',INEWTR,
     *         ' Y component removed.'
               CALL DZERO(WORK(K1TR1D),LTR1D)
            END IF
C
         ENDIF
C
         ENDIF
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
C
            IF ( IPRSOP .GE. 7 ) THEN
C
               CALL AROUND('Previous trial vector in RP_ORTH_TRN')
C
               WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &              (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
C
            END IF
C
C----------------------------------------------------------------
C           Orthogonalize new trial vector against previous trial 
C           vectors and their paired partners.
C----------------------------------------------------------------
C
            CALL RP_RES_O(WORK(KRSO1E),LRSO1E, WORK(KRSO1D),LRSO1D,
     &                    WORK(KTR1E), LTR1E,  WORK(KTR1D), LTR1D,
     &                    ISYRES,      ISYMTR)
C
            DOTP = DDOT(LTR1E,WORK(K1TR1E),1,WORK(KRSO1E),1)
     &           + DDOT(LTR1D,WORK(K1TR1D),1,WORK(KRSO1D),1)
C
            CALL DAXPY(LTR1E,-DOTP,WORK(KTR1E),1,WORK(K1TR1E),1)
            CALL DAXPY(LTR1D,-DOTP,WORK(KTR1D),1,WORK(K1TR1D),1)
C
C
            CALL RP_RES_O(WORK(KRSO1D),LRSO1D, WORK(KRSO1E),LRSO1E,
     &                    WORK(KTR1D), LTR1D,  WORK(KTR1E), LTR1E,
     &                    ISYRES,      ISYMTR)
C
            DOTP = DDOT(LTR1E,WORK(K1TR1E),1,WORK(KRSO1D),1)
     &           + DDOT(LTR1D,WORK(K1TR1D),1,WORK(KRSO1E),1)
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
         CALL RP_RES_O(WORK(KRSO1E),LRSO1E, WORK(KRSO1D),LRSO1D,
     &                 WORK(K1TR1E),LTR1E,  WORK(K1TR1D),LTR1D,
     &                 ISYRES,      ISYMTR)
C
         DNORM = DDOT(LTR1E,WORK(K1TR1E),1,WORK(KRSO1E),1)
     &         + DDOT(LTR1D,WORK(K1TR1D),1,WORK(KRSO1D),1)
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
         END IF
C
C------------------------------------------------------------------
C        If Norm is little Normalize new trial vector a first time.
C------------------------------------------------------------------
C
         IF ((LTYPE .EQ. 'LINEAR') .AND. (DABS(DNORM).LE.T1MIN)) THEN
C
            DNORMI = ONE / DSQRT( DABS(DNORM) )
C
            CALL DSCAL(LTR1E,DNORMI,WORK(K1TR1E),1)
            CALL DSCAL(LTR1D,DNORMI,WORK(K1TR1D),1)
C
            DNORM = DDOT(LTR1E,WORK(K1TR1E),1,WORK(KRSO1E),1)
     &            - DDOT(LTR1D,WORK(K1TR1D),1,WORK(KRSO1D),1)
C
         ENDIF
C------------------------------------------------------------
C        Normalize new trial vector a first (or second) time.
C------------------------------------------------------------
C
         DNORMI = ONE / DSQRT( DABS(DNORM) )
C
         CALL DSCAL(LTR1E,DNORMI,WORK(K1TR1E),1)
         CALL DSCAL(LTR1D,DNORMI,WORK(K1TR1D),1)
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
C
         IF ( IPRSOP .GE. 6 ) THEN
C
            CALL AROUND('Orthonormalized trialvector in RP_ORTH_TRN')
C
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(K1TR1E+I-1),WORK(K1TR1D+I-1),I=1,LTR1E)
C
         END IF
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_ORTH_TRN')
C
      RETURN
C
 9001 FORMAT(/,'Square NORM of new trial vector: ',F22.18,/)
 9002 FORMAT(/,'Norm of normalized new trial vector is: ',F22.18,/
     &       'New trial vector is removed because of linear ',
     &       'dependence.')
 9003 FORMAT(/,'WARNING: Problems orthonormalizing in RP_ORTH_TRN')
C
      END
