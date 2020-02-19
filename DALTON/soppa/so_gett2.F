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
#ifdef VAR_MPI
      use so_parutils, only: soppa_update_amplitudes
#endif
C
      use so_info, only : sop_mp2ai_done, sop_dp
C#include "implicit.h"
      implicit none

#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      REAL(SOP_DP), PARAMETER :: HALF = 0.5D0
C
C     Arguments
      character(len=4), intent(in) :: KEY
      real(sop_dp), intent(in) :: FOCKD(LFOCKD)
      real(sop_dp), intent(out) :: T2AM(LT2AM), DENSAI(LDENSAI),
     &                             DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      real(sop_dp), intent(inout) :: WORK(LWORK)
      integer, intent(in) :: LFOCKD, LT2AM, LDENSAI, LDENSIJ, LDENSAB,
     &                       LWORK

C     Local Variables
      LOGICAL DOCCSD, LEXIST
      CHARACTER(LEN=8) :: FNAME
      LOGICAL :: FOUND_DENS
      INTEGER :: IOSTATR, IDUMMY, LUTAM
      real(sop_dp) :: dtime
C     Pointers... Do we need them?
      integer :: KT1AM, KEND1
      integer :: LWORK1
C     Time-function.
      real(sop_dp) :: second
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_GETT2')
C
C-----------------------------------------
C     Memory allocation
C-----------------------------------------
C
      KT1AM = 1
      KEND1 = KT1AM + NT1AM(1)
      LWORK1 = LWORK - KEND1
C
      IF (LWORK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in GETT2')
      ENDIF
C
      IF ( KEY .EQ. 'MP2 ' ) THEN
         FNAME = 'MP2__TAM'
      ELSE IF ( KEY .EQ. 'CC2 ' ) THEN
         FNAME = 'CC2__TAM'
      ELSE IF ( KEY .EQ. 'CCSD' ) THEN
         FNAME = 'CCSD_TAM'
      ELSE
         CALL QUIT('Invalid KEY in SO_GETT2: '//KEY)
      END IF
C
C----------------------------------------------------------
C        Check if T2-amplitudes has been calculated in CC:
C----------------------------------------------------------
C
      INQUIRE (FILE=FNAME,EXIST=LEXIST)
      IF ((.NOT.LEXIST) .AND. (KEY .EQ. 'MP2 ')) THEN
C             WRITE(LUPRI, '(2A)') 'The (ai|bj) integrals for the T2 ',
C     & 'amplitudes is calculated in SOPPA'
C
C----------------------------------------
C        Calculate (ia|jb) integrals .
C----------------------------------------
C
         DTIME     = SECOND()
         CALL SO_IAJB(T2AM,WORK(KEND1),LWORK1)
         DTIME     = SECOND()  - DTIME
         SOTIME(5) = SOTIME(5) + DTIME
C
C-----------------------------------------
C        Calculate MP2 T2 amplitudes
C-----------------------------------------
C
         DTIME     = SECOND()
         CALL SO_MP2T2(FOCKD,LFOCKD,T2AM,LT2AM,WORK(KEND1),LWORK1)
         DTIME     = SECOND()  - DTIME
         SOTIME(5) = SOTIME(5) + DTIME
         FOUND_DENS = .FALSE.
C
         LUTAM = -1
         CALL GPOPEN(LUTAM,FNAME,'NEW',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND LUTAM
         WRITE (LUTAM) T2AM
C
      ELSE IF (LEXIST) THEN
C             WRITE(LUPRI, '(2A)') 'The (ai|bj) integrals for the T2 ',
C     & 'ampliltudes from is read from the file: MP2__TAM'
C-----------------------------------------------
C          Read T1 and T2-amplitudes from disk.
C-----------------------------------------------
         LUTAM = -1
         CALL GPOPEN(LUTAM,FNAME,'OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND LUTAM
C
C        T1 amplitudes only if from CC
C
         IF ( KEY(1:2) .EQ. 'CC' ) THEN
            READ(LUTAM,END=901) (DENSAI(I), I = 1,LDENSAI)
            CALL DSCAL(LDENSAI,HALF,DENSAI,1)
         END IF
         READ(LUTAM,END=902) (T2AM(I), I = 1,LT2AM)
C
C------------------------------------------------------------
C        Attempt to read the D(ij) and D(ab) density matrices
C------------------------------------------------------------
C
         READ(LUTAM,IOSTAT=IOSTATR) DENSIJ
         IF (IOSTATR.GT.0) THEN
            CALL QUIT('Error reading file '//FNAME)
         ELSE IF (IOSTATR.LT.0) THEN
            FOUND_DENS = .FALSE.
            BACKSPACE (LUTAM)
         ELSE
            READ(LUTAM) DENSAB
C
C        And also for MP2 get D(ai) block
            if ( (KEY(1:2).NE.'CC').and.sop_mp2ai_done) then
               READ(LUTAM,END=901) (DENSAI(I), I = 1,LDENSAI)
            end if
            FOUND_DENS = .TRUE.
         END IF
C
      ELSE
         CALL QUIT('SOPPA(CC) requires running the CC module')
      END IF

      IF (.NOT.FOUND_DENS) THEN
C
101      CONTINUE ! Failed to read D(ij) above
C
C----------------------------------------------------
C     Calculate the density matrices D(ij) and D(ab).
C----------------------------------------------------
C
         DTIME     = SECOND()
         CALL SO_DENS(DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                T2AM,LT2AM,WORK(KEND1),LWORK1)
         DTIME     = SECOND()  - DTIME
         SOTIME(4) = SOTIME(4) + DTIME
C
         WRITE (LUTAM) DENSIJ
         WRITE (LUTAM) DENSAB
C
      END IF
      CALL GPCLOSE(LUTAM,'KEEP')
C
C------------------------------------------------------
C     Replace T2-amplitudes in T2AM with MP2-amplitudes
C     (2 * coulomb - exchange of T2-amplitudes).
C------------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_T2MP(T2AM,LT2AM)
      DTIME     = SECOND()  - DTIME
      SOTIME(5) = SOTIME(5) + DTIME
C
C-----------------------------------------------------
C     Send the amplitudes to the slaves
C-----------------------------------------------------
#ifdef VAR_MPI
      call soppa_update_amplitudes(t2am, lt2am)
#endif
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_GETT2')
C
      RETURN

901   CALL QUIT('End of file when reading singles amplitudes')
902   CALL QUIT('End of file when reading doubles amplitudes')


      END

      SUBROUTINE SO_DUMP_DENSAI(KEY,DENSAI,LDENSAI)
C
C  Dumps the AI density to file
C
      use so_info, only : sop_dp
      implicit none
      character(len=4), intent(in) :: KEY
      real(sop_dp), intent(in) :: DENSAI(LDENSAI)
      integer, intent(in) :: LDENSAI
C
      character(len=8) :: FNAME
      INTEGER :: LUTAM, I, IDUMMY
C
C     Set filename from key
      FNAME = KEY//'_TAM'
      if (FNAME(4:4).EQ.' ') THEN
         FNAME(4:4) = '_'
      end if
      LUTAM = -1
C  GPOPEN SEEMS TO NOT ALLOW US TO SPECIFY OPEN FOR APPENDING,
      CALL GPOPEN(LUTAM,FNAME,'OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUTAM)
      READ(LUTAM) ! T2 amplitudes
      READ(LUTAM) ! D(ij)
      READ(LUTAM) ! D(ab)
C
      WRITE(LUTAM) (DENSAI(I), I = 1,LDENSAI)
      CALL GPCLOSE(LUTAM,'KEEP')
C
      END

      SUBROUTINE SO_READ_DENS(KEY,DENSIJ,LDENSIJ,
     &                        DENSAB,LDENSAB,DENSAI,LDENSAI)
C
      use so_info, only: sop_dp
      implicit none
      character(len=4), intent(in) :: KEY
      real(sop_dp), intent(out) :: DENSIJ(LDENSIJ),
     &                            DENSAB(LDENSAB), DENSAI(LDENSAI)
      integer, intent(in) :: LDENSIJ, LDENSAB, LDENSAI
C
      real(sop_dp), parameter :: HALF = 0.5D0
      character(len=8) :: FNAME
      INTEGER :: LUTAM, I, IDUMMY
C
C     Set filename from key
      FNAME = KEY//'_TAM'
      if (FNAME(4:4).EQ.' ') THEN
         FNAME(4:4) = '_'
      end if
      LUTAM = -1
C
      CALL GPOPEN(LUTAM,FNAME,'OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUTAM)

      IF (KEY(1:2).EQ.'CC') THEN
C    CC method - T1 amplitudes are first
         READ(LUTAM) DENSAI
         CALL DSCAL(LDENSAI,HALF,DENSAI,1)
      END IF
C     Skip T2 amplitudes
      READ(LUTAM)
C     Read densij and densab
      READ(LUTAM) DENSIJ
      READ(LUTAM) DENSAB
C
C     In case of MP2: densai is last
      IF (KEY(1:3).EQ.'MP2') THEN
         READ(LUTAM) DENSAI
      END IF
      CALL GPCLOSE(LUTAM,'KEEP')

      RETURN

      END
