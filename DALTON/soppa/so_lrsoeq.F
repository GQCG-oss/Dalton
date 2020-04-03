C
C  /* Deck so_lrsoeq */
      SUBROUTINE SO_LRSOEQ(MODEL,LABEL,ISYMTR,IMAGPROP,
     &                     NEXCI,MAXIT,FRVAL,NFRVAL,
     &                     RESINM,LRESINM,CONV,LCONV,DENSIJ,LDENSIJ,
     &                     DENSAB,LDENSAB,DENSAI,LDENSAI,T2MP,LT2MP,
     &                     FOCKD,LFOCKD,REDE,REDS,LMXRED,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan Sauer and Keld Bak, May 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Andrea Ligabue, December 2003: linear response functions
C                                    implemented
C
C     PURPOSE: Solve the SOPPA linear response equations
C              using an AO-driven algorithm.
C
C       LABEL   label of the property to be computed
C       ISYMTR  symmetry
C       NEXCI   (number of excitation) ... for us is always 1
C       MAXIT   max # iteration
C       FRVAL   NFRVAL values of the frequencies
C       NFRVAL  number of frequencies
C       RESINM  LRESIMN values ...
C       LRESINM equal to NFRVAL
C       CONV    array of 8 elements
C       LCONV   set to 8
C       FOCKD
C       LFOCKD
C       REDE    matrix E[2]
C       REDS    matrix S[2]
C       LMXRED
C
#ifdef VAR_MPI
      use so_parutils, only: parsoppa_do_eres, my_mpi_integer,
     &                       soppa_nint, sop_master, one_mpi
#endif
      use so_info, only: sop_linres, so_full_name, so_has_doubles,
     &                   so_needs_densai, sop_dp, sop_stat_trh,
     &                   so_model_number,
     &                   sop_mp2ai_done
C
#include "implicit.h"
#ifdef VAR_MPI
#include "mpif.h"
C  IRAT in order to assign space for load-balancing
#include "iratdef.h"
#endif
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      LOGICAL   NONEWT,IMAGPROP, IMAGPROP2
      LOGICAL   DOUBLES
      LOGICAL   UPDATE_GP
C
      CHARACTER*3 CONV(LCONV)
      CHARACTER MODEL*5,LABEL*8
      CHARACTER*11 FULL_NAME
C
      DIMENSION FRVAL(NFRVAL),   RESINM(LRESINM)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION T2MP(LT2MP),     FOCKD(LFOCKD),   WORK(LWORK)
      DIMENSION REDE(LMXRED),    REDS(LMXRED)
C
      REAL(sop_dp), parameter :: TWO = 2.0D0, HALF = 0.5D0, ONE = 1.0D0
C
      INTEGER NIT, NOLDTR, NNEWTR, IDTYPE, IMODEL
#ifdef VAR_MPI
      INTEGER :: CP_ISYMTR
      INTEGER :: MAXNUMJOBS
C     This array is only there to ensure that the six above variables
C     are allocated consecutively, so that it can be send together. Only
C     use it for this purpose.
C     The definition must match that in soppa_nodedriver
      INTEGER   INFO_ARRAY(6)
      EQUIVALENCE (info_array(1), cp_isymtr), (info_array(2),nit),
     &            (info_array(3), nnewtr),    (info_array(4),noldtr),
     &            (info_array(5), idtype), (info_array(6),imodel)
      INTEGER(MPI_INTEGER_KIND) :: ierr_mpi, numprocs_mpi
#endif
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_LRSOEQ')
C
C     Basic info on what we do here
      DOUBLES = so_has_doubles(model)
      FULL_NAME = SO_FULL_NAME(MODEL)
      IMODEL = SO_MODEL_NUMBER(MODEL)
C
C     Should we update the GP vector
      UPDATE_GP = (.NOT. sop_mp2ai_done) .AND.
     &            SO_NEEDS_DENSAI(MODEL)
C
C==============================================================
C     For checking, calculate E[2] and S[2] matrices explicitly
C     by using unit vectors as trial vectors.
C==============================================================
C
C
      IF (SOPCHK) THEN
C
         CALL SO_CHECK(MODEL,DENSIJ,LDENSIJ,DENSAB,LDENSAB,T2MP,LT2MP,
     &                 FOCKD,LFOCKD,ISYMTR,WORK,LWORK)
C
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
C-------------------------
C     1. Memory allocation
C-------------------------
C
      KREDC   = 1
      KEND1   = KREDC + (2*NSAVMX)**2
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
         lAssignedIndices = (maxnumjobs + IRAT - 1) /IRAT
         kAssignedIndices = KEND1
         KEND1 = kAssignedIndices + lAssignedIndices
      endif
#endif
      LWORK1 = LWORK - KEND1
C
      CALL SO_MEMMAX('SO_LRSOEQ.1',LWORK1)
      IF(LWORK1.LT.0) CALL STOPIT('SO_LRSOEQ.1','  ',KEND1,LWORK)
C------------------
C     Write banner.
C------------------
C
      WRITE(LUPRI,'(A)') ''
C      WRITE(LUPRI,'(/,2X,A)') '*********************************'//
C     &                        '*********************************'
      WRITE(LUPRI,9001)
      WRITE(LUPRI,'(11X,3A)') ADJUSTR(FULL_NAME),
     &        ' iterations, Property ',LABEL
C      WRITE(LUPRI,'(2X,A)') '*********************************'//
C     &                        '*********************************'
      WRITE(LUPRI,9001)
      WRITE(LUPRI,9011)
      WRITE(LUPRI,9001)
C
C------------------------------
C     Loop over the frequencies
C------------------------------
C
      DO 300 IFREQ = 1,NFRVAL
C
C------------------------------------------------------------------
C     Set up initial trial vectors and write to disk using eq. (19)
C     with (C +-C) as R or the last TV of the previous frequency
C------------------------------------------------------------------
C
      DTIME      = SECOND()
      IF (MODEL.EQ.'AORPA') THEN
         CALL RP_TRIAL3(NNEWTR,NOLDTR,ISYMTR,IMAGPROP,IFREQ,
     &                  FRVAL,NFRVAL,
     &                  NEXCI,WORK(KEND1),LWORK1)
      ELSE
         CALL SO_TRIAL3(MODEL,NNEWTR,ISYMTR,IMAGPROP,IFREQ,FRVAL,NFRVAL,
     &                  NEXCI,!LABEL,
     &                  DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                  WORK(KEND1),LWORK1)
      ENDIF
      DTIME      = SECOND()   - DTIME
cLig  <> we have to define and write properly the index insotime array
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
      IF (ABS(FRVAL(IFREQ)).GT.sop_stat_trh) then
         IDTYPE = 0
      ELSE IF (IMAGPROP) THEN
         IDTYPE = 2
      ELSE
         IDTYPE = 1
      END IF
C
C----------------------------------------------------------------------
C     Iteration loop or solving the linear equation/eigenvalue problem
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
         IF ( IPRSOP .GE. 5 ) THEN
C
            WRITE(LUPRI,'(/,2X,A)') '================================'//
     &                              '=================================='
            WRITE(LUPRI,'(11X,I3,3A,I1,A,F8.5)') NIT,'. ',
     &              TRIM(FULL_NAME),
     &              '  iteration, Symmetry ', ISYMTR,
     &              ', Frequency ', FRVAL(IFREQ)
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
     &                IDTYPE,
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
         IF (UPDATE_GP) THEN
            UPDATE_GP = .FALSE.
C
C-----------------------------------------------------------------------
C          Calculate the gradient property vectors with the right DENSAI
C-----------------------------------------------------------------------
C          Only the one-partcle part need to be re-done ->
C          Maybe pass in MODEL = 'AOHRP'?
           LGPVC1  = NT1AM(ISYMTR)
           LGPVC2 = 0
           KGPVC1  = KEND1
           KGPVC2  = KGPVC1 + LGPVC1
           KEND2   = KGPVC2 + LGPVC2
           LWORK2  = LWORK  - KEND2
C
           CALL SO_MEMMAX ('SO_LRSOEQ.2',LWORK2)
           IF (LWORK2 .LT.0) CALL STOPIT('SO_LRSOEQ.2',' ',KEND2,LWORK)
C
           CALL SO_GETGP(WORK(KGPVC1),LGPVC1,WORK(KGPVC2),LGPVC2,LABEL,
     &                  ISYMTR,IMAGPROP2,'AOHRP',
     &                  T2MP,LT2MP,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                  DENSAI,LDENSAI,WORK(KEND2),LWORK2)
           IF ( .NOT.(IMAGPROP.EQV.IMAGPROP2) )
     &         WRITE(LUPRI,'(A)') 'WARNING: Problem with property type'
C
C---------------------------------------------------------
C          Save the complete property gradients on a file.
C---------------------------------------------------------
C
           CALL SO_OPEN(LUGPV1,FNGPV1,LGPVC1)
C
           CALL SO_WRITE(WORK(KGPVC1),LGPVC1,LUGPV1,FNGPV1,1)
C
           IF(IPRSOP.GE.5) THEN
              CALL AROUND('In SO_LRSOEQ:  '//MODEL//
     &                    ' gradient property vector.'//LABEL)
              CALL OUTPUT(WORK(KGPVC1),1,LGPVC1,1,1,LGPVC1,1,1,LUPRI)
           ENDIF
C
           CALL SO_CLOSE(LUGPV1,FNGPV1,'KEEP')
C
         ENDIF
C
C----------------------------------------------------
C        Set up and solve the reduced linear problem.
C----------------------------------------------------
C
         LREDE  = 2 * ( NOLDTR + NNEWTR )
         LREDS  = 2 * ( NOLDTR + NNEWTR )
         LREDC  = 2 * ( NOLDTR + NNEWTR )
C
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
         CALL SO_REDLE(DOUBLES,NEXCI,NOLDTR,NNEWTR,
     &                 LABEL,ISYMTR,IMAGPROP,
     &                 REDE,LREDE,
     &                 REDS,LREDS,WORK(KREDC),LREDC,LREDOL,FRVAL,NFRVAL,
     &                 IFREQ,ENORM,PROP,WORK(KEND1),LWORK1)
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
C        Determine the residues from the current optimal solution
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
cLig <> I have to add some parameters to that call and modify so_trial2
         CALL SO_TRIAL2(MODEL,'LINEAR',NONEWT,NOLDTR,NNEWTR,NLINDP,
     &                  FRVAL(IFREQ),1,RESINM,LRESINM,CONV,LCONV,NCONV,
     &                  ISYMTR,IMAGPROP,
     &                  NEXCI,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                  ENORM,LABEL,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(34) = SOTIME(34) + DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(4)  = SOWTIM(4)  + WTIMEE - WTIMES
C
C---------------------------------------------------------------------
C        Write residuals and "Diagonal properties"
C---------------------------------------------------------------------
C
C
C            WRITE(LUPRI,9010)
C            WRITE(LUPRI,9011)
C            WRITE(LUPRI,9010)
         WRITE(LUPRI,9012) NIT,FRVAL(IFREQ),PROP,RESINM(1),CONV(1)
C            WRITE(LUPRI,9010)
C
C
C---------------------------------------
C        Flush the standard output file.
C---------------------------------------
C
C         CALL FLSHFO(LUPRI)
C
C---------------------------------------------------------------------
C     Go to next iteration if the solution vector is not converged and
C     the maximum number of iterations have not been reached.
C---------------------------------------------------------------------
C
      IF ( (NNEWTR .GT. 0) .AND. (NIT .LT. MAXIT) ) GO TO 100
C
C--------------------------------------------------------------
C     Calculate and save the pertubed density matrix
C-------------------------------------------------------------
C
      if (triplet.AND.doubles) THEN
         CALL SO_TMLTR(T2MP,ONE,1)
         CALL DSCAL(LT2MP,HALF,T2MP,1)
      END IF
C
      ENORMINV = 1.0D0/ENORM
      CALL SO_PERTDENS(MODEL,SOP_LINRES,1,
     &                 FRVAL(IFREQ),1,LABEL,
     &                 ISYMTR,IMAGPROP,ENORMINV,
     &                 T2MP,LT2MP,DENSIJ,LDENSIJ,
     &                 DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                 WORK(KEND1),LWORK1)
C
      if (triplet.and.doubles) CALL CCSD_TCMEPKX(T2MP,ONE,1)
C
C
      IF ((NNEWTR .EQ. 0) .AND. (NIT.LT.MAXIT)) THEN
C
         IF ( NLINDP .EQ. 0 ) THEN
C
            WRITE(LUPRI,9001)
C            WRITE(LUPRI,9002)
C            WRITE(LUPRI,9003)
C            WRITE(LUPRI,9006) NIT
C            WRITE(LUPRI,9008)
C
         ELSE
C
            WRITE(LUPRI,9001)
            WRITE(LUPRI,9002)
            WRITE(LUPRI,9003)
            WRITE(LUPRI,9004)
C            WRITE(LUPRI,9005) NEXCI - NLINDP
            WRITE(LUPRI,9008)
C
         END IF
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
         CALL QUIT('ERROR occured in SO_LRSOEQ')
C
      END IF
C
C----------------------------------------
C    End of the loop over the frequencies
C----------------------------------------
C
  300 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_LRSOEQ')
C
      RETURN
C
 9001 FORMAT(1X,'---------------------------------------------',
     &       '----------------')
 9002 FORMAT(1X,'Iterations stopped since: ')
 9003 FORMAT(1X,'-------------------------')
 9004 FORMAT(1X,'The new trial vector is linearly dependent',
     &       ' on the previous one.')
!9006 FORMAT(1X,'The property converged after, ',I3,' iterations.')
 9007 FORMAT(1X,'Maximum number of ',I3,' iterations is reached.')
 9008 FORMAT(1X,'---------------------------------------------',
     &       '----------------',/)
 9010 FORMAT('--------------------------------------------------',
     &       '-----------------')
 9011 FORMAT(1X,'Iteration Frequency    Property (au)   Residual',
     &           '    Converged')
C Change format of Property to D15.8 (from F15.8) to avoid unreadable
 9012 FORMAT(2X,I5,4X,F8.6,2X,D15.8,1X,1P,D14.4,6X,A3)
      END
