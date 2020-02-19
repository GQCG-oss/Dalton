C  /* Deck so_check */
      SUBROUTINE SO_CHECK(MODEL,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
CPi 04.04.16
C     &                    T2MP,LT2MP,
     &                    DENSAI,LDENSAI,T2MP,LT2MP,
Cend-Pi
     &                    FOCKD,LFOCKD,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan P.A. Sauer.                    1-Nov-1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     Test routine that calculates the E[2] and S[2] matrices
C     explicitly by carrying out linear transformations on unit trial
C     vectors.
C
      use so_info, only: so_has_doubles
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
      CHARACTER*5 MODEL
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), T2MP(LT2MP)
      DIMENSION FOCKD(LFOCKD),   WORK(LWORK)
      LOGICAL DOUBLES
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
      INTEGER(MPI_INTEGER_KIND) :: ierr_mpi, numprocs_mpi
#endif
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_CHECK')
C
      DOUBLES = SO_HAS_DOUBLES(MODEL)
      CHKTIM = SECOND()
      RESTIM = ZERO
C
C-----------------------------------------------------
C     Work space allocation for E[2] and S[2] matrices
C     and trial vectors.
C-----------------------------------------------------
C
      IF (DOUBLES) THEN
         NVAR   = 2 * (NT1AM(ISYMTR)+N2P2HOP(ISYMTR))
      ELSE
         NVAR = 2*NT1AM(ISYMTR)
      END IF
C
      LE2MAT = NVAR * NVAR
      LS2MAT = NVAR * NVAR
      LTR1E  = NT1AM(ISYMTR)
      LTR1D  = NT1AM(ISYMTR)
CPi 04.04.16
C      LDENSAI = NT1AM(1)
Cend-Pi
      IF (DOUBLES) THEN
         LTR2E  = N2P2HOP(ISYMTR)
         LTR2D  = N2P2HOP(ISYMTR)
      ELSE
         LTR2E  = 0
         LTR2D  = 0
      END IF
C
      KE2MAT = 1
      KS2MAT = KE2MAT + LE2MAT
      KTR1E  = KS2MAT + LS2MAT
      KTR1D  = KTR1E  + LTR1E
      KTR2E  = KTR1D  + LTR1D
      KTR2D  = KTR2E  + LTR2E
      KDENSAI= KTR2D  + LTR2D
      KEND1  = KDENSAI+ LDENSAI
      LWORK1 = LWORK  - KEND1
C
#ifdef VAR_MPI
C------------------------------------------------------------------
C     For MPI, we need som space in which to store the indices each
C     process is to work with in so_eres.
C------------------------------------------------------------------
C
      call mpi_comm_size( mpi_comm_world, numprocs_mpi, ierr_mpi)
      maxnumjobs = mxcall - min(mxcall, numprocs_mpi) + 1
      if ( numprocs_mpi .eq. 1 ) then
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
C
      CALL SO_MEMMAX ('SO_CHECK.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_CHECK.1',' ',KEND1,LWORK)
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
C      Pointer for the correct column in E and S
         LE2MAT  = (ICOL - 1) * NVAR
         LS2MAT  = (ICOL - 1) * NVAR
C
         LRES1E  = NT1AM(ISYMTR)
         LRES1D  = NT1AM(ISYMTR)
         LRESO1E = NT1AM(ISYMTR)
         LRESO1D = NT1AM(ISYMTR)
         IF (DOUBLES) THEN
            LRES2E  = N2P2HOP(ISYMTR)
            LRES2D  = N2P2HOP(ISYMTR)
            LRESO2E = N2P2HOP(ISYMTR)
            LRESO2D = N2P2HOP(ISYMTR)
         ELSE
            LRES2E  = 0
            LRES2D  = 0
            LRESO2E = 0
            LRESO2D = 0
         ENDIF
C
C      Positions on various parts in E and S
         KRES1E  = KE2MAT  + LE2MAT
         KRES1D  = KRES1E  + LRES1E
         KRES2E  = KRES1D  + LRES1D
         KRES2D  = KRES2E  + LRES2E
         KRESO1E = KS2MAT  + LS2MAT
         KRESO1D = KRESO1E + LRESO1E
         KRESO2E = KRESO1D + LRESO1D
         KRESO2D = KRESO2E + LRESO2E
         KEND2   = KRESO2D + LRESO2D
         LWORK2  = LWORK   - KEND2
C
         CALL SO_MEMMAX ('SO_CHECK.2',LWORK2)
         IF (LWORK2 .LT. 0) CALL STOPIT('SO_CHECK.2',' ',KEND2,LWORK)
C
         CALL DZERO(WORK(KTR1E),LTR1E)
         CALL DZERO(WORK(KTR1D),LTR1D)
         IF (DOUBLES) THEN
            CALL DZERO(WORK(KTR2E),LTR2E)
            CALL DZERO(WORK(KTR2D),LTR2D)
            CALL DZERO(WORK(KRESO2E),LRESO2E)
            CALL DZERO(WORK(KRESO2D),LRESO2D)
         ENDIF
C
         IF (ICOL.LE.LTR1E) THEN
            WORK(KTR1E + ICOL - 1) = ONE
         ELSE IF (ICOL.LE.(LTR1E+LTR1D)) THEN
            WORK(KTR1D + (ICOL - LTR1E) - 1) = ONE
         ELSE IF (ICOL.LE.(LTR1E+LTR1D+LTR2E)) THEN ! DOUBLES ONLY
            WORK(KTR2E + (ICOL - LTR1E - LTR1D) - 1) = ONE
            WORK(KRESO2E + ICOL - LTR1E - LTR1D -1)  = ONE
         ELSE
            WORK(KTR2D + (ICOL - LTR1E - LTR1D - LTR2E) - 1) = ONE
            WORK(KRESO2D + (ICOL - LTR1E - LTR1D - LTR2E) - 1) = -ONE
         ENDIF
C
         CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
         CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
         CALL SO_WRITE(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,1)
         CALL SO_WRITE(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,1)
         CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
         CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C
         IF (DOUBLES) THEN
            CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
            CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
            CALL SO_WRITE(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,1)
            CALL SO_WRITE(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,1)
            CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
            CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
         ENDIF
C
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
     &                   mpi_comm_world, ierr_mpi )
C Pack INFO_ARRAY variables and send them off
         CP_ISYMTR = ISYMTR
         NOLDTR    = 0
         NNEWTR    = 1
         NIT       = ICOL
         CALL MPI_BCAST( INFO_ARRAY, 4, MY_MPI_INTEGER, 0,
     &                   MPI_COMM_WORLD, ierr_mpi)
#endif
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
         CALL SO_ERES(MODEL,0, 1, DENSIJ, LDENSIJ, DENSAB, LDENSAB,
     &               T2MP, LT2MP,  FOCKD,   LFOCKD,
     &               DENSAI, LDENSAI, ICOL, ISYMTR,
     &               0,
#ifdef VAR_MPI
     &               WORK(kAssignedIndices),lAssignedIndices,
#endif
     &               WORK(KEND1), LWORK1)
         ETIME      = SECOND()   - DTIME
         SOTIME(35) = SOTIME(35) + ETIME
         DTIME      = SECOND()
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(1)  = SOWTIM(1)  + WTIMEE - WTIMES
C
         CALL GETTIM (DUMMY,WTIMES)
C     Explicitly calculate S for singles if not RPA
         IF (MODEL.NE.'AORPA'.AND.(ICOL.LE.(LTR1E+LTR1D))) THEN
            CALL SO_SRES(0,1,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                ISYMTR,WORK(KEND1),LWORK1)
            CALL SO_OPEN(LURO1E,FNRO1E,LRESO1E)
            CALL SO_OPEN(LURO1D,FNRO1D,LRESO1D)
            CALL SO_READ(WORK(KRESO1E),LRESO1E,LURO1E,FNRO1E,1)
            CALL SO_READ(WORK(KRESO1D),LRESO1D,LURO1D,FNRO1D,1)
            CALL SO_CLOSE(LURO1E,FNRO1E,'DELETE')
            CALL SO_CLOSE(LURO1D,FNRO1D,'DELETE')
         ELSE
            CALL DZERO(WORK(KRESO1E),LTR1E+LTR1D)
            IF (MODEL.EQ.'AORPA') THEN
               IF (ICOL.LE.LRESO1E) THEN
                  WORK(KRESO1E-1+ICOL) = ONE
               ELSE
                  WORK(KRESO1D-1+(ICOL-LRESO1E)) = -ONE
               ENDIF
            END IF
         END IF

         DTIME      = SECOND()   - DTIME
         SOTIME(40) = SOTIME(40) + DTIME
         RESTIM = RESTIM + SECOND() - DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(2)  = SOWTIM(2)  + WTIMEE - WTIMES
C
         CALL SO_OPEN(LURS1E,FNRS1E,LRES1E)
         CALL SO_OPEN(LURS1D,FNRS1D,LRES1D)
         CALL SO_READ(WORK(KRES1E), LRES1E, LURS1E,FNRS1E,1)
         CALL SO_READ(WORK(KRES1D), LRES1D, LURS1D,FNRS1D,1)
         CALL SO_CLOSE(LURS1E,FNRS1E,'DELETE')
         CALL SO_CLOSE(LURS1D,FNRS1D,'DELETE')
C
         IF (DOUBLES) THEN
            CALL SO_OPEN(LURS2E,FNRS2E,LRES2E)
            CALL SO_OPEN(LURS2D,FNRS2D,LRES2D)
            CALL SO_READ(WORK(KRES2E), LRES2E, LURS2E,FNRS2E,1)
            CALL SO_READ(WORK(KRES2D), LRES2D, LURS2D,FNRS2D,1)
            CALL SO_CLOSE(LURS2E,FNRS2E,'DELETE')
            CALL SO_CLOSE(LURS2D,FNRS2D,'DELETE')
         ENDIF
C
C
C
  100 CONTINUE
C
      CALL AROUND('E[2] Matrix')
C
      WRITE(LUPRI,'(A,I2,A,I8)')' SYMMETRY :',ISYMTR,'  DIMENSION ',NVAR
      CALL OUTPUT(WORK(KE2MAT),1,NVAR,1,NVAR,NVAR,NVAR,1,LUPRI)
C
CPFP 2009-Aug  Write the C matrices
C
      CALL AROUND('Ce1[2] Matrix')
      WRITE(LUPRI,'(A,I2,A,5I8)')' SYMMETRY :',ISYMTR,
     &            '  DIMENSIONS ',LRESO1E,LRESO1D,LRESO2E,LRESO2D,NVAR
      CALL OUTPUT(WORK(KE2MAT),(2*LRESO1E+1),(2*LRESO1E+LRESO2E),
     &                         1,LRESO1D,NVAR,NVAR,1,LUPRI)
C
      CALL AROUND('~Cd1[2] Matrix')
      WRITE(LUPRI,'(A,I2,A,5I8)')' SYMMETRY :',ISYMTR,
     &            '  DIMENSIONS ',LRESO1E,LRESO1D,LRESO2E,LRESO2D,NVAR
      CALL OUTPUT(WORK(KE2MAT),1,LRESO1E,(2*LRESO1D+1),
     &            (2*LRESO1E+LRESO2E),NVAR,NVAR,1,LUPRI)
C
      CALL AROUND('Ce2[2] Matrix')
      WRITE(LUPRI,'(A,I2,A,5I8)')' SYMMETRY :',ISYMTR,
     &            '  DIMENSIONS ',LRESO1E,LRESO1D,LRESO2E,LRESO2D,NVAR
      CALL OUTPUT(WORK(KE2MAT),(2*LRESO1E+LRESO2E+1),NVAR,(LRESO1D+1),
     &            (2*LRESO1D),NVAR,NVAR,1,LUPRI)
C
      CALL AROUND('~Cd2[2] Matrix')
      WRITE(LUPRI,'(A,I2,A,5I8)')' SYMMETRY :',ISYMTR,
     &            '  DIMENSIONS ',LRESO1E,LRESO1D,LRESO2E,LRESO2D,NVAR
      CALL OUTPUT(WORK(KE2MAT),(LRESO1E+1),(2*LRESO1E),
     &            (2*LRESO1D+LRESO2D+1),NVAR,NVAR,NVAR,1,LUPRI)
C
CPi 02.05.16: Write D matrix
C
      CALL AROUND('D1[0] Matrix')
      WRITE(LUPRI,'(A,I2,A,5I8,/A,5I8)')' SYMMETRY :',ISYMTR,
     &            '  DIMENSIONS ',LRESO1E,LRESO1D,LRESO2E,LRESO2D,NVAR,
     &            ' DIMENSION Triplet operators',
     &            NT2AMTT(ISYMTR),NT2AMT1(ISYMTR),NT2AMT2(ISYMTR),
     &            NT2AMT3(ISYMTR)
      CALL OUTPUT(WORK(KE2MAT),(2*LRESO1E+1),(2*LRESO1E+LRESO2E),
     &            (2*LRESO1D+1),(2*LRESO1E+LRESO2E),NVAR,NVAR,1,LUPRI)
C
      CALL AROUND('D2[0] Matrix')
      WRITE(LUPRI,'(A,I2,A,5I8,/A,5I8)')' SYMMETRY :',ISYMTR,
     &            '  DIMENSIONS ',LRESO1E,LRESO1D,LRESO2E,LRESO2D,NVAR,
     &            ' DIMENSION Triplet operators',
     &            NT2AMTT(ISYMTR),NT2AMT1(ISYMTR),NT2AMT2(ISYMTR),
     &            NT2AMT3(ISYMTR)
      CALL OUTPUT(WORK(KE2MAT),(2*LRESO1E+LRESO2E+1),NVAR,
     &            (2*LRESO1D+LRESO2D+1),NVAR,NVAR,NVAR,1,LUPRI)
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
      CALL QEXIT('SO_CHECK')
C
      RETURN
      END
