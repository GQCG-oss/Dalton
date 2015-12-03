C
C  /* Deck rp_lrsoeq */
      SUBROUTINE RP_LRSOEQ(LABEL,ISYMTR,IMAGPROP,NEXCI,MAXIT,
     &                     FRVAL,NFRVAL,RESINM,
     &                     LRESINM,CONV,LCONV,FOCKD,LFOCKD,
     &                     REDE,REDS,LMXRED,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Andrea Ligabue, December 2002
C
C     PURPOSE: Solve the RPA linear response equations 
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
     &                       soppa_nint
#endif
      use so_info, only: sop_linres
#include "implicit.h"
#ifdef VAR_MPI
#include "mpif.h"
C  IRAT in order to assign space for load-balancing
#include "iratdef.h"
#endif
C      
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      PARAMETER (EPSI = 1.0D-8)
C
      LOGICAL   NONEWT, IMAGPROP
C
      CHARACTER*3 CONV(LCONV)
      CHARACTER*8 LABEL
C
      DIMENSION FRVAL(NFRVAL),  RESINM(LRESINM)
      DIMENSION REDE(LMXRED),   REDS(LMXRED)
      DIMENSION FOCKD(LFOCKD),  WORK(LWORK)
      INTEGER   NIT, NOLDTR, NNEWTR
#ifdef VAR_MPI
      INTEGER   CP_ISYMTR
C     This array is only there to ensure that the four above variables
C     are allocated consecutively, so that it can be send together. Only
C     use it for this purpose.
C     The definition must match that in soppa_nodedriver
      INTEGER   INFO_ARRAY(4)
      EQUIVALENCE (info_array(1), cp_isymtr), (info_array(2),nit), 
     &            (info_array(3), nnewtr),    (info_array(4),noldtr)
      INTEGER(MPI_INTEGER_KIND) :: IERR
#endif
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_LRSOEQ')
C
C-----------------------------------------------------------------
C     Calculate diagonal parts of E[2] and S[2] and write to disk.
C-----------------------------------------------------------------
C
      LDUMMY  = 0
      DTIME      = SECOND()
      CALL SO_DIAG('AORPA',FOCKD,LFOCKD,DUMMY,LDUMMY,DUMMY,LDUMMY,
     &              ISYMTR,WORK,LWORK)
      DTIME      = SECOND()   - DTIME
      SOTIME(31) = SOTIME(31) + DTIME
C
      KREDC   = 1
      KEND1   = KREDC  + (2*NSAVMX)**2
#ifdef VAR_MPI
C------------------------------------------------------------------
C     For MPI, we need som space in which to store the indices each
C     process is to work with in so_eres.
C------------------------------------------------------------------
C
      call mpi_comm_size( mpi_comm_world, numprocs, ierr)
      maxnumjobs = soppa_nint - min(soppa_nint, numprocs) + 1
      if ( numprocs .eq. 1 ) then
C Not a real parallel job, don't bother
         lAssignedIndices = 1
         kAssignedIndices = 0
      else  
         lAssignedIndices = (maxnumjobs + 1) /IRAT
         kAssignedIndices = KEND1
         KEND1 = kAssignedIndices + lAssignedIndices
      endif
#endif 
      LWORK1 = LWORK - KEND1
C
      CALL SO_MEMMAX ('RP_LRSOEQ.1',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT('RP_LRSOEQ.1',' ',KEND1,LWORK)
C
C------------------------------------------------------------------
C     Loop over the frequencies
C------------------------------------------------------------------
C
      DO 300 IFREQ = 1, NFRVAL
C
C------------------------------------------------------------------
C       Check tu the value of NSVMX ... 
C       rise it to 8 if FRVAL(IFREQ) eq 0
C------------------------------------------------------------------
CRF     This seems to no longer be neccesary!
C       Also changing NSAVMX here cause all kinds of trouble!
        NSAVMX = NSAVMXORIG
C
C
C        IF ( DABS(FRVAL(IFREQ)) .LE. EPSI ) THEN
C                
C            IF (NSAVMX.LT.8) THEN
C
C               NSAVMX = MIN(NT1AM(ISYMTR),8)
C
C               WRITE(LUPRI,'(/,1X,A,I4,A,/,1X,A,/,1X,A,I4,A,I4)')
C     &         'WARNING: when frequency is 0 ',NSAVMXORIG,' vectors',
C     &         ' are not enough to get a good description of the space',
C     &         ' I have temporally rised from ',NSAVMXORIG,' to',
C     &         NSAVMX,' the number of trial vector'
C
C               IF(NSAVMX.LT.8) WRITE(LUPRI,'(/,1X,A,/,A,I5)') 
C     &         'I cannot rise NSAVMX to 8 becouse',
C     &         'the numbre of 1p1h is ',NT1AM(ISYMTR)
C            ENDIF
C
C        ENDIF
C
C------------------------------------------------------------------
C     Set up initial trial vectors using eq. (19) with (C +-C) as R
C     or just using the final ones from the previous frequency
C------------------------------------------------------------------
C
      DTIME      = SECOND()
      CALL RP_TRIAL3(NNEWTR,NOLDTR,ISYMTR,IMAGPROP,IFREQ,FRVAL,NFRVAL,
     &               NEXCI,WORK(KEND1),LWORK1)
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
C-----------------------------------------------------------
C     Iteration loop for solving the linear equation problem
C-----------------------------------------------------------
C
  100 CONTINUE
C
C-----------------------------------
C        Count number of iterations.
C-----------------------------------
C
         NIT = NIT + 1
C
         WRITE(LUPRI,'(/,2X,A,/,9X,I3,A,I1,A,F8.5,/,2X,A,/)') 
     &         '====================================================='//
     &         '=============',NIT,
     &         '. RPA iteration, Property symmetry ',ISYMTR,' Freq ',
     &         FRVAL(IFREQ),'==================='//
     &         '==============================================='
         CALL FLSHFO(LUPRI)
C
C--------------------------------------------------------------
C        Make E[2] linear transformation of trialvectors giving 
C        resultvectors.
C--------------------------------------------------------------
C
#ifdef VAR_MPI
C In parallel, send slaves to so_eres
C
         call mpi_bcast( parsoppa_do_eres, 1, my_mpi_integer, 0,
     &                   mpi_comm_world, ierr )
C ISYMTR is a non-local parameter, we need to copy it to the info-array
         CP_ISYMTR = ISYMTR
         CALL MPI_BCAST( INFO_ARRAY, 4, MY_MPI_INTEGER, 0,
     &                   MPI_COMM_WORLD, IERR)
#endif
         DTIME      = SECOND()
         CALL RP_ERES(NOLDTR,NNEWTR,FOCKD,LFOCKD,NIT,ISYMTR,
#ifdef VAR_MPI
     &                WORK(kAssignedIndices),maxnumjobs,
#endif
     &                WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(35) = SOTIME(35) + DTIME
         CALL FLSHFO(LUPRI)
C
C
C--------------------------------------------------------------
C        Make S[2] linear transformation of trialvectors giving
C        resultvectors.
C--------------------------------------------------------------
CRF      This is trivial, waste of time to write it to and read it        
C        from disk
C         DTIME      = SECOND()
C         CALL RP_SRES(NOLDTR,NNEWTR,ISYMTR,WORK(KEND1),LWORK1)
c         DTIME      = SECOND()   - DTIME
C         SOTIME(40) = SOTIME(40) + DTIME
C         CALL FLSHFO(LUPRI)
C
C
C--------------------------------------------------------
C        Set up and solve the reduced eigenvalue problem.
C--------------------------------------------------------
C
         LREDE  = 2 * ( NOLDTR + NNEWTR )
         LREDS  = 2 * ( NOLDTR + NNEWTR )
         LREDC  = 2 * ( NOLDTR + NNEWTR )
C
         DTIME      = SECOND()
C
         CALL SO_REDLE(.FALSE.,NEXCI,NOLDTR,NNEWTR,LABEL,ISYMTR,
     &                 IMAGPROP,REDE,LREDE,REDS,
     &                 LREDS,WORK(KREDC),LREDC,LREDOL,FRVAL,NFRVAL,
     &                 IFREQ,ENORM,PROP,WORK(KEND1),LWORK1)
Cend-PFP
         DTIME      = SECOND()   - DTIME
         SOTIME(33) = SOTIME(33) + DTIME
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
C        vectors and decide if convergence is obtained.
C        If non-converged create new trial-vectors. 
C        These are orthogonalized against the previous optimal 
C        trial-vectors and among themself including
C        the vectors obtained from pairing.
C-------------------------------------------------------------------
C        
         IF (NIT .GE. MAXIT) NONEWT = .TRUE.
C
         LDUMMY = 1    
         DTIME      = SECOND()
         CALL SO_TRIAL2('AORPA','LINEAR',NONEWT,NOLDTR,NNEWTR,NLINDP,
     &                  FRVAL(IFREQ),1,RESINM,LRESINM,CONV,LCONV,NCONV,
     &                  ISYMTR,IMAGPROP,
     &                  NEXCI,DUMMY,LDUMMY,DUMMY,LDUMMY,
     &                  ENORM,LABEL,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(34) = SOTIME(34) + DTIME
C
C---------------------------------------------------------------------
C        Write calculated excitation energies and residuals to output.
C---------------------------------------------------------------------
C
         WRITE(LUPRI,9010)
         WRITE(LUPRI,9011)
         WRITE(LUPRI,9010)
         WRITE(LUPRI,9012) LABEL,FRVAL(IFREQ),PROP,RESINM(1),CONV(1)
         WRITE(LUPRI,9010)
C
C---------------------------------------
C        Flush the standard output file.
C---------------------------------------
C
         CALL FLSHFO(LUPRI)
C
C---------------------------------------------------------------------
C     Go to next iteration if all eigenvalues are not converged and if
C     the maximum number of iterations have not been reached.
C---------------------------------------------------------------------
C
      IF ( (NNEWTR .GT. 0) .AND. (NIT .LT. MAXIT) ) GO TO 100
C
C--------------------------------------------------------------      
C     Calculate and save the pertubed density matrix
C-------------------------------------------------------------
C 
      ENORMINV = 1.0D0/ENORM 
      CALL SO_PERTDENS('AORPA',SOP_LINRES,1,
     &                 FRVAL(IFREQ),1,LABEL,
     &                 ISYMTR,IMAGPROP,ENORMINV,
     &                 T2MP,LT2MP,DENSIJ,LDENSIJ,
     &                 DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                 WORK(KEND1),LWORK1)
C
      IF ( (NNEWTR .EQ. 0) .AND. (NIT.LT.MAXIT) ) THEN
C
         IF ( NLINDP .EQ. 0 ) THEN
C
            WRITE(LUPRI,9001)
            WRITE(LUPRI,9002)
            WRITE(LUPRI,9003)
            WRITE(LUPRI,9006) NIT
            WRITE(LUPRI,9008)
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
      ELSE ! (NNEWTR < 0 and NIT != MAXIT) or NIT > MAXIT  
C
         CALL QUIT('ERROR occured in RP_LRSOEQ')
C
      END IF
C
      CALL FLSHFO(LUPRI)
C
  300 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('RP_LRSOEQ')
C
      RETURN
C
 9001 FORMAT(/,1X,'---------------------------------------------',
     &       '----------------')
 9002 FORMAT(1X,'Iterations stopped since: ')
 9003 FORMAT(1X,'-------------------------')
 9004 FORMAT(1X,'The new trial vector is linearly dependent',
     &       ' on the previous one.')
 9006 FORMAT(1X,'The property converged after, ',I3,' iterations.')
 9007 FORMAT(1X,'Maximum number of ',I3,' iterations is reached.')
 9008 FORMAT(1X,'---------------------------------------------',
     &       '----------------',/)
 9010 FORMAT(4X,'--------------------------------------------------',
     &          '------------')
 9011 FORMAT(4X,'Label      Freq        Property (au)   Residual',
     &           '      Converged')
 9012 FORMAT(4X,A8,2X,F8.6,2X,F15.8,1X,1P,D14.4,6X,A)
      END
