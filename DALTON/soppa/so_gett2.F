C
C  /* Deck so_gett2 */
      SUBROUTINE SO_GETT2(KEY,FOCKD,LFOCKD,T2AM,LT2AM,DENSAI,LDENSAI,
     &                    DENSIJ,LDENSIJ,DENSAB,LDENSAB,WORK,LWORK)
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, December 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Eiter read from disk or calculate the appropriate T2 and possibly T1 amplitudes
C              for RPA(D), SOPPA and SOPPA(CCSD) and calculate
C              the corresponding densities.
C
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      PARAMETER (HALF = 0.5D0,ESUDIP = 64604.885D0,ESUECD = 471.44360D0)
C
      CHARACTER KEY*4
      DIMENSION FOCKD(LFOCKD),T2AM(LT2AM),DENSAI(LDENSAI)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION WORK(LWORK)
      LOGICAL DOCCSD,LEXIST
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_GETT2')
C--------------------------------------------
C     If Triplet T2AM will get double length
C--------------------------------------------
C
      LT2AMH = LT2AM
      IF (TRIPLET) LT2AMH = INT(LT2AM/2)
C
C-----------------------------------------
C     Memory allocation
C-----------------------------------------
C
      KEND1 = 1
      KT2AM = KEND1
      KEND2 = KT2AM + LT2AMH
      KT1AM = KEND2
      LWORK2 = LWORK - KEND3
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient space in GETT2')
      ENDIF
C
      IF ( KEY(1:3) .EQ. 'MP2' ) THEN
C
C----------------------------------------------------------
C        Check if T2-amplitudes has been calculated in CC:
C----------------------------------------------------------
C
         IF (DOCCSD) THEN
           INQUIRE (FILE='MP2__TAM',EXIST=LEXIST)
           IF (LEXIST) THEN
             WRITE(LUPRI, '(A)') 'The (ai|bj) integrals for the T2
     & ampliltudes from is read from the file: MP2__TAM'
C-----------------------------------------------
C          Read MP2 T2-amplitudes from disk.
C------------------------------------------------
             LUTAM = -1
             CALL GPOPEN(LUTAM,'MP2__TAM','OLD',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
             REWIND LUTAM
C
             READ(LUTAM)
             READ(LUTAM) (T2AM(I), I = 1,LT2AM)
             READ(LUTAM) (T2AM(I), I = 1,LT2AMH)
             CALL GPCLOSE(LUTAM,'KEEP')
C
           ELSE IF (.NOT. LEXIST) THEN
             WRITE(LUPRI, '(A)') 'The (ai|bj) integrals for the T2
     & amplitudes is calculated in SOPPA'
C
C----------------------------------------
C        Calculate (ia|jb) integrals .
C----------------------------------------
             CALL DZERO(WORK(KT1AM),NT1AMX)
C
             DTIME     = SECOND()
             CALL SO_IAJB(T2AM,WORK(KT1AM),WORK,LWORK2)
             DTIME     = SECOND()  - DTIME
             SOTIME(5) = SOTIME(5) + DTIME
           END IF
C
C-------------------------------------------------------
C        Calculate (ia|jb) integrals .
C-------------------------------------------------------
C 
         ELSE IF (.NOT. DOCCSD) THEN
           WRITE(LUPRI, '(A)') 'The (ai|bj) integrals for the T2
     & amplitudes is calculated in SOPPA'
           CALL DZERO(WORK(KT1AM),NT1AMX)
C
           DTIME     = SECOND()
           CALL SO_IAJB(T2AM,WORK(KT1AM),WORK,LWORK2)
           DTIME     = SECOND()  - DTIME
           SOTIME(5) = SOTIME(5) + DTIME
C
         END IF
C
C-----------------------------------------
C        Calculate MP2 T2 amplitudes
C-----------------------------------------
C
         DTIME     = SECOND()
         CALL SO_MP2T2(FOCKD,LFOCKD,T2AM,LT2AMH)
         DTIME     = SECOND()  - DTIME
         SOTIME(5) = SOTIME(5) + DTIME
C
      ELSE IF ( KEY(1:3) .EQ. 'CC2' ) THEN
C
C--------------------------------------------------
C        Read CC2 T1- and T2-amplitudes from disk.
C--------------------------------------------------
C
         LUTAM = -1
         CALL GPOPEN(LUTAM,'CC2__TAM','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUTAM
C
         READ(LUTAM) (DENSAI(I), I = 1,LDENSAI)
         CALL DSCAL(LDENSAI,HALF,DENSAI,1)
         READ(LUTAM) (T2AM(I), I = 1,LT2AMH)
         CALL GPCLOSE(LUTAM,'KEEP')
C
      ELSE IF ( KEY(1:4) .EQ. 'CCSD' ) THEN
C
C--------------------------------------------------
C        Read CCSD T1- and T2-amplitudes from disk.
C--------------------------------------------------
C
         LUTAM = -1
         CALL GPOPEN(LUTAM,'CCSD_TAM','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUTAM
C
         READ(LUTAM) (DENSAI(I), I = 1,LDENSAI)
         CALL DSCAL(LDENSAI,HALF,DENSAI,1)
         READ(LUTAM) (T2AM(I), I = 1,LT2AMH)
         CALL GPCLOSE(LUTAM,'KEEP')

      END IF
C
C----------------------------------------------------
C     Calculate the density matrices D(ij) and D(ab).
C----------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_DENS(DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &             T2AM,LT2AMH,WORK(KEND2),LWORK2)
      DTIME     = SECOND()  - DTIME
      SOTIME(4) = SOTIME(4) + DTIME
C
C------------------------------------------------------
C     Replace T2-amplitudes in T2AM with MP2-amplitudes
C     (2 * coulomb - exchange of T2-amplitudes).
C------------------------------------------------------
C
      DTIME     = SECOND()
C Using LT2AMH here and do not initialize LT2MPH in so_t2mp
      CALL SO_T2MP(T2AM,LT2AMH)
      DTIME     = SECOND()  - DTIME
      SOTIME(5) = SOTIME(5) + DTIME
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_GETT2')
C
      RETURN
      END
