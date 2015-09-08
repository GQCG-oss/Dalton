C
C  /* Deck rp_check */
      SUBROUTINE RP_CHECK(FOCKD,LFOCKD,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan P.A. Sauer.                    4-Dec-2003
C
C     Test routine that calculates the E[2] and S[2] matrices for RPA
C     explicitly by carrying out linear transformations on unit trial
C     vectors.
C
#ifdef VAR_MPI
      use so_parutils, only: parsoppa_do_eres, my_mpi_integer
#endif
#include "implicit.h"
#ifdef VAR_MPI
#include "mpif.h"
C  We need MXCALL and IRAT in order to assign space for load-balancing
C  
#include "maxorb.h"
#include "distcl.h"
#include "iratdef.h"
#endif
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION FOCKD(LFOCKD),   WORK(LWORK)
#include "ccorb.h"
#include "ccsdsym.h"
#include "soppinf.h"
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
      CALL QENTER('RP_CHECK')
C
      CHKTIM = SECOND()
      RESTIM = ZERO
C
C-----------------------------------------------------
C     Work space allocation for E[2] and S[2] matrices
C     and trial vectors.
C-----------------------------------------------------
C
      NVAR   = 2 * NT1AM(ISYMTR)
C
      LE2MAT = NVAR * NVAR
      LS2MAT = NVAR * NVAR
      LTR1E  = NT1AM(ISYMTR)
      LTR1D  = NT1AM(ISYMTR)
C
      KE2MAT  = 1
      KS2MAT  = KE2MAT + LE2MAT
      KTR1E   = KS2MAT + LS2MAT
      KTR1D   = KTR1E  + LTR1E
      KEND1   = KTR1D  + LTR1D
      LWORK1  = LWORK  - KEND1
C
#ifdef VAR_MPI
C------------------------------------------------------------------
C     For MPI, we need som space in which to store the indices each
C     process is to work with in so_eres.
C------------------------------------------------------------------
C
      call mpi_comm_size( mpi_comm_world, numprocs, ierr)
      maxnumjobs = mxcall - min(mxcall, numprocs) + 1
      if ( numprocs .eq. 1 ) then
C Not a real parallel job, don't bother
         lAssignedIndices = 1
         kAssignedIndices = 0
      else  
         lAssignedIndices = (maxnumjobs + 1) /IRAT
         kAssignedIndices = KEND1
         KEND1 = kAssignedIndices + lAssignedIndices
         LWORK1 = LWORK - KEND1
      endif
#endif 
      CALL SO_MEMMAX ('RP_CHECK.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_CHECK.1',' ',KEND1,LWORK)
C
C---------------------------------------------------
C     Loop over columns of the E[2] and S[2] matrix.
C---------------------------------------------------
C
      DO 100 ICOL = 1,NVAR
C
C---------------------------------------------------
C        Work space allocation for solution vectors.
C---------------------------------------------------
C
         LE2MAT  = (ICOL - 1) * NVAR
         LS2MAT  = (ICOL - 1) * NVAR
         LRES1E  = NT1AM(ISYMTR)
         LRES1D  = NT1AM(ISYMTR)
         LRESO1E = NT1AM(ISYMTR)
         LRESO1D = NT1AM(ISYMTR)
C
         KRES1E  = KE2MAT  + LE2MAT
         KRES1D  = KRES1E  + LRES1E
         KRESO1E = KS2MAT  + LS2MAT
         KRESO1D = KRESO1E + LRESO1E
         KEND2   = KRESO1D + LRESO1D
         LWORK2  = LWORK   - KEND2
C
         CALL SO_MEMMAX ('RP_CHECK.2',LWORK2)
         IF (LWORK2 .LT. 0) CALL STOPIT('RP_CHECK.2',' ',KEND2,LWORK)
C
         CALL DZERO(WORK(KTR1E),LTR1E)
         CALL DZERO(WORK(KTR1D),LTR1D)
C
         IF (ICOL.LE.LTR1E) THEN
            WORK(KTR1E + ICOL - 1) = ONE
         ELSE IF (ICOL.LE.(LTR1E+LTR1D)) THEN
            WORK(KTR1D + (ICOL - LTR1E) - 1) = ONE
         ENDIF
C
         CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
         CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
C
         CALL SO_WRITE(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,1)
         CALL SO_WRITE(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,1)
C
         CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
         CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
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
C Pack INFO_ARRAY variables and send them off 
         CP_ISYMTR = ISYMTR
         NOLDTR    = 0
         NNEWTR    = 1
         NIT       = ICOL 
         CALL MPI_BCAST( INFO_ARRAY, 4, MY_MPI_INTEGER, 0,
     &                   MPI_COMM_WORLD, IERR)
#endif
         DTIME      = SECOND()
         CALL RP_ERES(0,1,FOCKD,LFOCKD,ICOL,ISYMTR,
#ifdef VAR_MPI
     &                WORK(kAssignedIndices),lAssignedIndices,
#endif
     &                WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(35) = SOTIME(35) + DTIME
C
C
C--------------------------------------------------------------
C        Make S[2] linear transformation of trialvectors giving
C        resultvectors.
C--------------------------------------------------------------
C
         DTIME      = SECOND()
         CALL RP_SRES(0,1,ISYMTR,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(40) = SOTIME(40) + DTIME
C
C
         CALL SO_OPEN(LURS1E,FNRS1E,LRES1E)
         CALL SO_OPEN(LURS1D,FNRS1D,LRES1D)
         CALL SO_OPEN(LURO1E,FNRO1E,LRESO1E)
         CALL SO_OPEN(LURO1D,FNRO1D,LRESO1D)
C
         CALL SO_READ(WORK(KRES1E), LRES1E, LURS1E,FNRS1E,1)
         CALL SO_READ(WORK(KRES1D), LRES1D, LURS1D,FNRS1D,1)
         CALL SO_READ(WORK(KRESO1E),LRESO1E,LURO1E,FNRO1E,1)
         CALL SO_READ(WORK(KRESO1D),LRESO1D,LURO1D,FNRO1D,1)
C
         CALL SO_CLOSE(LURS1E,FNRS1E,'DELETE')
         CALL SO_CLOSE(LURS1D,FNRS1D,'DELETE')
         CALL SO_CLOSE(LURO1E,FNRO1E,'DELETE')
         CALL SO_CLOSE(LURO1D,FNRO1D,'DELETE')
C
  100 CONTINUE
C
      CALL AROUND('E[2] Matrix')
C
      WRITE(LUPRI,'(A,I2,A,I8)')' SYMMETRY :',ISYMTR,'  DIMENSION ',NVAR
      CALL OUTPUT(WORK(KE2MAT),1,NVAR,1,NVAR,NVAR,NVAR,1,LUPRI)
C
      CALL AROUND('S[2] Matrix')
C
      WRITE(LUPRI,'(A,I2,A,I8)')' SYMMETRY :',ISYMTR,'  DIMENSION ',NVAR
      CALL OUTPUT(WORK(KS2MAT),1,NVAR,1,NVAR,NVAR,NVAR,1,LUPRI)
C
      CHKTIM     = SECOND() - CHKTIM - RESTIM
      SOTIME(27) = SOTIME(27) + CHKTIM
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_CHECK')
C
      RETURN
      END
