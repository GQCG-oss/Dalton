C  /* Deck so_lrsolv */
      SUBROUTINE SO_LRSOLV(MODEL,ISYMTR,NEXCI,MAXIT,EXVAL,LEXVAL,RESINM,
     &                     LRESINM,CONV,LCONV,DENSIJ,LDENSIJ,DENSAB,
     &                     LDENSAB,DENSAI,LDENSAI,T2MP,LT2MP,FOCKD,
     &                     LFOCKD,REDE,REDS,LMXRED,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan Sauer and Keld Bak, May 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Andrea Ligabue, December 2003: linear response functions
C                                    implemented
C
C     PURPOSE: Solve the linear response equations
C              using an AO-driven algorithm.
C
#ifdef VAR_MPI
      use so_parutils, only: parsoppa_do_eres, my_mpi_integer,
     &                       soppa_nint, sop_master, one_mpi
#endif
C
      use so_info, only: so_full_name, sop_excita, so_has_doubles,
     &                   so_model_number
#include "implicit.h"
#ifdef VAR_MPI
#include "mpif.h"
C
C MXCALL has been replaced by SOPPA_NINT for hermite compatibility
!#include "maxorb.h"
!#include "distcl.h"
C  IRAT in order to assign space for load-balancing
#include "iratdef.h"
#endif
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      LOGICAL   NONEWT
      LOGICAL   DOUBLES
C
      CHARACTER*3 CONV(LCONV)
      CHARACTER*5 MODEL
      CHARACTER*8 PDENS_LABEL
      CHARACTER*11 FULL_NAME
C
      DIMENSION EXVAL(LEXVAL),   RESINM(LRESINM)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION T2MP(LT2MP),     FOCKD(LFOCKD),   WORK(LWORK)
      DIMENSION REDE(LMXRED),    REDS(LMXRED)

      INTEGER NIT, NOLDTR, NNEWTR, IDTYPE, IMODEL
#ifdef VAR_MPI
      INTEGER   CP_ISYMTR
      INTEGER ::  MAXNUMJOBS
C     This array is only there to ensure that the four above variables
C     are allocated consecutively, so that it can be send together. Only
C     use it for this purpose.
C     The definition must match that in soppa_nodedriver
      INTEGER   INFO_ARRAY(6)
      EQUIVALENCE (info_array(1), cp_isymtr), (info_array(2),nit),
     &            (info_array(3), nnewtr),    (info_array(4),noldtr),
     &            (info_array(5), idtype), (info_array(6),imodel)
      INTEGER(MPI_INTEGER_KIND) :: ierr_mpi, numprocs_mpi
      IDTYPE = 0
#endif
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_LRSOLV')
C
      DOUBLES = SO_HAS_DOUBLES(MODEL)
C
      FULL_NAME = SO_FULL_NAME(MODEL)
      IMODEL = SO_MODEL_NUMBER(MODEL)
C
C==============================================================
C     For checking, calculate E[2] and S[2] matrices explicitly
C     by using unit vectors as trial vectors.
C==============================================================
C
      IF (SOPCHK) THEN
         CALL SO_CHECK(MODEL,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
CPi 04.04.16
C     &                 T2MP,LT2MP,
     &                 DENSAI,LDENSAI,T2MP,LT2MP,
Cend-Pi
     &                 FOCKD,LFOCKD,ISYMTR,WORK,LWORK)
      END IF
C
C-----------------------------------------------------------------
C     Calculate diagonal parts of E[2] and S[2] and write to disk.
C-----------------------------------------------------------------
C
      DTIME      = SECOND()
      CALL SO_DIAG(MODEL,FOCKD,LFOCKD,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &             ISYMTR,WORK,LWORK)
      DTIME      = SECOND()   - DTIME
      SOTIME(31) = SOTIME(31) + DTIME
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LPOINT = NEXCI
C
      KPOINT = 1
      KEND1  = KPOINT + LPOINT
      LWORK1 = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_LRSOLV.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_LRSOLV.1',' ',KEND1,LWORK)
C
C----------------------------------
C     Set up initial trial vectors.
C----------------------------------
C
      DTIME      = SECOND()
      IF (MODEL .EQ. 'AORPA') THEN
         CALL RP_TRIAL1(NNEWTR,WORK(KPOINT),LPOINT,ISYMTR,
     &                  NEXCI,WORK(KEND1),LWORK1)
      ELSE
         CALL SO_TRIAL1(MODEL,NNEWTR,WORK(KPOINT),LPOINT,ISYMTR,
     &                  NEXCI,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                  WORK(KEND1),LWORK1)
      ENDIF
      DTIME      = SECOND()   - DTIME
      SOTIME(32) = SOTIME(32) + DTIME
C
C--------------------------------------------------------------
C     Initialize iteration counter, number of old trialvectors,
C     and logical for 'no new trial vectors.'
C--------------------------------------------------------------
C
      NIT    = 0
C
      NOLDTR = 0
C
      LREDOL = 0
C
      NONEWT = .FALSE.
C
#ifdef VAR_MPI
C------------------------------------------------------------------
C     For MPI, we need some space in which to store the indices each
C     process is to work with in so_eres.
C------------------------------------------------------------------
C
      call mpi_comm_size( mpi_comm_world, numprocs_mpi, ierr_mpi)
      maxnumjobs = soppa_nint - min(soppa_nint, numprocs_mpi) + 1
      if ( numprocs_mpi .eq. 1 ) then
C Not a real parallel job, don't bother
         lAssignedIndices = 1
         kAssignedIndices = 0
      else
         lAssignedIndices = (maxnumjobs + 1) /IRAT
         kAssignedIndices = KEND1
         KEND1 = kAssignedIndices + lAssignedIndices
         LWORK1 = LWORK - KEND1
         CALL SO_MEMMAX ('SO_LRSOLV.1A',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT('SO_LRSOLV.1A',' ',KEND1,LWORK)
      endif
#endif
C------------------
C     Write banner.
C------------------
C
      WRITE(LUPRI,'(/,2X,A)') '*********************************'//
     &                        '*********************************'
         WRITE(LUPRI,'(11X,A,A,I1)') ADJUSTR(FULL_NAME),
     &        ' iterations, Excitation symmetry ',ISYMTR
      WRITE(LUPRI,'(2X,A)') '*********************************'//
     &                        '*********************************'
C
C----------------------------------------------------------------------
C     Iteration loop for solving the linear equation/eigenvalue problem
C----------------------------------------------------------------------
C
  100 CONTINUE
C
C--------------------------------------------------------------
C        Count number of iterations and write header to output.
C--------------------------------------------------------------
C
         NIT = NIT + 1
C
         IF ( IPRSOP .GE. 2 ) THEN
C
            WRITE(LUPRI,'(/,2X,A)') '================================'//
     &                              '=================================='
            WRITE(LUPRI,'(11X,I3,3A,I1)') NIT,
     &              '. ', TRIM(FULL_NAME),
     &              '  iteration, Excitation symmetry ',
     &              ISYMTR
            WRITE(LUPRI,'(2X,A,/)') '================================'//
     &                              '=================================='
C
         END IF
C
C--------------------------------------------------------------
C        Make E[2] linear transformation of trialvectors giving
C        resultvectors.
C--------------------------------------------------------------
C
#ifdef VAR_MPI
C In parallel, send slaves to so_eres
C
         call mpi_bcast( parsoppa_do_eres, one_mpi, my_mpi_integer,
     &                   sop_master,
     &                   mpi_comm_world, ierr_mpi )
C ISYMTR is a non-local parameter, we need to copy it to the info-array
         CP_ISYMTR = ISYMTR
         CALL MPI_BCAST( INFO_ARRAY, 6_mpi_integer_kind, MY_MPI_INTEGER,
     &                   sop_master, MPI_COMM_WORLD, ierr_mpi)
#endif
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
         CALL SO_ERES(MODEL,NOLDTR,NNEWTR,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                T2MP,LT2MP,FOCKD,LFOCKD,DENSAI,LDENSAI,NIT,ISYMTR,
     &                0,
#ifdef VAR_MPI
     &                WORK(kAssignedIndices),maxnumjobs,
#endif
     &                WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(35) = SOTIME(35) + DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(1)  = SOWTIM(1)  + WTIMEE - WTIMES
C
         IF ( AOTEST ) THEN
C
C--------------------------------------------------------------------
C           Test orthonormality of trial vectors and check the linear
C           transformed trial vectors.
C--------------------------------------------------------------------
C
            CALL SO_TEST(NOLDTR,NNEWTR,ISYMTR,DENSIJ,LDENSIJ,DENSAB,
     &                   LDENSAB,T2MP,LT2MP,FOCKD,LFOCKD,WORK(KEND1),
     &                   LWORK1 )
C
         END IF
C
C--------------------------------------------------------
C        Set up and solve the reduced eigenvalue problem.
C--------------------------------------------------------
C
         LREDE  = 2 * ( NOLDTR + NNEWTR )
         LREDS  = 2 * ( NOLDTR + NNEWTR )
C
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
         CALL SO_REDEV(DOUBLES,NEXCI,NOLDTR,NNEWTR,ISYMTR,
     &                 REDE,LREDE,REDS,LREDS,
     &                 LREDOL,EXVAL,LEXVAL,CONV,LCONV,
     &                 WORK(KEND1),LWORK1)
C
         DTIME      = SECOND()   - DTIME
         SOTIME(33) = SOTIME(33) + DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(3)  = SOWTIM(3)  + WTIMEE - WTIMES
C
         LREDOL = LREDE
C
C------------------------------------------
C        Reset number of old trial vectors.
C------------------------------------------
C
         NOLDTR = MIN(NOLDTR + NNEWTR, (NSAVMX - 1) * NEXCI )
C
C-------------------------------------------------------------------
C        Determine the residual from the current optimal solution
C        vectors and decide if convergence is obtained for any of
C        the vectors. For the non-converged vectors create new
C        trial-vectors. These are orthogonalized against the
C        previous optimal trial-vectors and among themself including
C        the vectors obtained from pairing.
C-------------------------------------------------------------------
C
         IF (NIT .GE. MAXIT) NONEWT = .TRUE.
C
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
         CALL SO_TRIAL2(MODEL,'EXCITA',NONEWT,NOLDTR,NNEWTR,NLINDP,
     &                  EXVAL,LEXVAL,RESINM,LRESINM,CONV,LCONV,
     &                  NCONV,ISYMTR,DUMMY,
     &                  NEXCI,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                  DUMMY,DUMMY,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(34) = SOTIME(34) + DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(4)  = SOWTIM(4)  + WTIMEE - WTIMES
C
C---------------------------------------------------------------------
C        Write calculated excitation energies and residuals to output.
C---------------------------------------------------------------------
C
         IF (IPRSOP .GE. 2) THEN
C
            WRITE(LUPRI,9010)
            WRITE(LUPRI,9011)
            WRITE(LUPRI,9010)
            DO 200 IEXCI = 1,NEXCI
               WRITE(LUPRI,9012) IEXCI,EXVAL(IEXCI),RESINM(IEXCI),
     &                           CONV(IEXCI)
 200        CONTINUE
            WRITE(LUPRI,9010)
C
         END IF
C
C---------------------------------------
C        Flush the standard output file.
C---------------------------------------
C
         CALL FLSHFO(LUPRI)
C
C---------------------------------------------------------------------
C     Go to next iteration if all eigenvalues are not converged or if
C     only complex excitations remain and if the maximum number of
C     iterations have not been reached.
CPi 10.08.16: In the case of complex excitations NNEWTR will always
C     be greater than 0. Use instead nr of converged excitations since
C     this in increased also in case of complex excitations (SO_REDEV)
C---------------------------------------------------------------------
C
C      IF ( (NNEWTR .GT. 0) .AND. (NIT .LT. MAXIT) ) GOTO 100
      IF ( (NCONV .LT. NEXCI) .AND. (NIT .LT. MAXIT) ) GOTO 100
Cend-Pi
C
C--------------------------------------------------------------
C     Calculate and save the pertubed density matrix
C-------------------------------------------------------------
C
      WRITE(PDENS_LABEL,'(A7,I1)') 'EXCITA ',ISYMTR
      CALL SO_PERTDENS(MODEL,SOP_EXCITA,NEXCI,
     &                 EXVAL,NEXCI,PDENS_LABEL,
     &                 ISYMTR,.FALSE.,1.0D0,
     &                 T2MP,LT2MP,DENSIJ,LDENSIJ,
     &                 DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                 WORK(KEND1),LWORK1)
C
C--------------------------------------------------------------
C     Write appropriate iteration info to output
C-------------------------------------------------------------
C
C--------------------------------------------------------------
C     All excitations have converged
C-------------------------------------------------------------
C
      IF ((NNEWTR .EQ. 0) .AND. (NIT.LT.MAXIT)) THEN
C
         IF ( NLINDP .EQ. 0 ) THEN
C
            WRITE(LUPRI,9001)
            WRITE(LUPRI,9002)
            WRITE(LUPRI,9003)
            WRITE(LUPRI,9006)
            WRITE(LUPRI,9008)
C
         ELSE
C
            WRITE(LUPRI,9001)
            WRITE(LUPRI,9002)
            WRITE(LUPRI,9003)
            WRITE(LUPRI,9004) NLINDP
            WRITE(LUPRI,9005) NEXCI - NLINDP
            WRITE(LUPRI,9008)
C
         END IF
CPi 10.08.16
C-------------------------------------------------------------------
C     All excitations have either converged or are imaginary.
C     Comment: Since the imaginary excitations are not converged,
C     new trial vectors are construted in the previous iteration and
C     this case will not be included in the previous if statement.
C-------------------------------------------------------------------
C
      ELSE IF ((NCONV .EQ. NEXCI) .AND. (NIT .LT. MAXIT) ) THEN
C
         WRITE(LUPRI,9001)
         WRITE(LUPRI,9002)
         WRITE(LUPRI,9003)
         WRITE(LUPRI,9009)
         WRITE(LUPRI,9008)
C
Cend-Pi
C--------------------------------------------------------------
C     Maximum number of iterations is reached
C-------------------------------------------------------------
C
      ELSE IF (NIT .EQ. MAXIT) THEN
C
         WRITE(LUPRI,9001)
         WRITE(LUPRI,9002)
         WRITE(LUPRI,9003)
         WRITE(LUPRI,9007) MAXIT
         WRITE(LUPRI,9008)
C
         NNEWTR = 0
C
      ELSE
C
         CALL QUIT('ERROR occured in SO_LRSOLV')
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_LRSOLV')
C
      RETURN
C
 9001 FORMAT(/,1X,'---------------------------------------------',
     &       '----------------')
 9002 FORMAT(1X,'Iterations stopped since: ')
 9003 FORMAT(1X,'-------------------------')
 9004 FORMAT(1X,I2,' of the excitations have linear dependent new',
     &       ' trial vectors.')
 9005 FORMAT(1X,I2,' of the excitations are converged.')
 9006 FORMAT(1X,'All of the excitations are converged.')
 9007 FORMAT(1X,'Maximum number of ',I3,' iterations is reached.')
 9009 FORMAT(1X,'All of the excitations are either converged or ',
     &       'imaginary.')
 9008 FORMAT(1X,'---------------------------------------------',
     &       '----------------',/)
 9010 FORMAT(12X,'--------------------------------------------------')
 9011 FORMAT(12X,'Excitation   Energy (au)     Residual    Converged')
 9012 FORMAT(13X,I5,2X,F15.8,1X,1P,D14.4,6X,A)
      END
