C  /* Deck so_eres */
      SUBROUTINE SO_ERES(MODEL,  NOLDTR, NNEWTR,  DENSIJ,  LDENSIJ,
     &                   DENSAB, LDENSAB, T2MP,    LT2MP,   FOCKD,
     &                   LFOCKD, DENSAI,  LDENSAI, NIT,     ISYMTR,
     &                   IDTYPE,
#ifdef VAR_MPI
     &                   AssignedIndices, maxnumjobs,
#endif
     &                   WORK,   LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, October 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Frederik Beyer & Stephan P. A. Sauer: 27.08.2013: call to ERI corrected
C
C     PURPOSE: Driver routine for making a linear transformation of
C              a trialvector with the SOPPA hessian matricx E[2].
C              The trial vector consists of four parts TR1E, TR1D,
C              TR2E, and TR2D. E refers to excitations and D to
C              de-excitations. 1 refer to the one-particle part and
C              2 to the two-particle part. The linear transformed
C              trialvector is refered to as the resultvector and is
C              kept in four corresponding arrays. For the linear
C              transformation with E[2] the result vector is in RES1E,
C              RES1D, RES2E, and RES2D.
C              The linear transformation is driven over atomic orbitals,
C              and E[2] is not constructed explicitly.
C
      use so_info, only: so_singles_first, so_singles_second,
     &                   so_has_doubles, so_needs_densai,
     &                   so_double_correction,
     &                   sop_mp2ai_done

#ifdef VAR_MPI
      use so_parutils, only: soppa_comm_active, soppa_num_active,
     &                       soppa_nint, my_mpi_integer, sop_master
#ifdef USE_MPI_MOD_F90
      use mpi
#endif
#endif

#include "implicit.h"

#ifdef VAR_MPI
#ifndef USE_MPI_MOD_F90
#include "mpif.h"
#endif
#endif

#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), T2MP(LT2MP)
      DIMENSION FOCKD(LFOCKD)
      DIMENSION DENSAI(LDENSAI) !intent(inout)
      DIMENSION WORK(LWORK)
      CHARACTER MODEL*5
      integer :: inewtr, nnewtr, isymd1, ntosym
      integer :: thisinewtr
C     idtype = 0 for dynamic (handle D part explicitly),
C     idtype = 1 for real, static (Tr1D = - Tr1E)
C     idtype = 2 for imaginary, static ( Tr1D = Tr1E )
      integer, intent(in) :: idtype

#ifdef VAR_MPI
      integer :: maxnumjobs, nloopidx
      integer :: AssignedIndices(maxnumjobs)
      integer(kind=MPI_INTEGER_KIND) :: ierr_mpi, count_mpi
#endif
C
#include "ccorb.h"
#include "infind.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "eritap.h"
#include "soppinf.h"


      logical :: singles_first, singles_second, doubles, calc_densai
C
C     Logical variable which can be set to false if the dexcitation vector is
C     to avoid calculating with only zeroes. Useful for static
C     properties and first iteration of excitation energies.
      logical :: do_dex
      logical :: tr2_zero
#ifdef VAR_MPI
#include "iprtyp.h"
#include "infpar.h"
C     integer, save :: numprocs
      integer :: numprocs
      logical :: loadbal_dyn
      double precision ::  timeini, timefin
      integer(mpi_integer_kind) :: LTRTOT_mpi, lrestot, latot
      integer(mpi_integer_kind) :: my_MPI_REAL8 = MPI_REAL8
      numprocs = nodtot + 1
C
C  When to do dynamic load_balancing..?
C  The sooner the better, but is the first iteration close to the
C  average?
      loadbal_dyn = .false.
      if (nit .eq. 1) loadbal_dyn = .true.
#endif
C
      do_dex = idtype .eq. 0

      calc_densai = SO_NEEDS_DENSAI(MODEL) .AND. .NOT. SOP_MP2AI_DONE
C
C------------------
C     Add to trace.
C------------------
C
#ifdef VAR_MPI
      if (mynum .eq. 0 ) then
#endif
      CALL QENTER('SO_ERES')
#ifdef VAR_MPI
C     Slaves need to zero DENSAI if it is (re)calculated.
      elseif ( calc_densai ) THEN
         CALL DZERO(DENSAI,LDENSAI)
      endif
#endif
C
C-------------------------------------------------------------
C     Determine which terms to incude
C-------------------------------------------------------------
C
      singles_first = so_singles_first(model)
      singles_second = so_singles_second(model)
      doubles = so_has_doubles(model)
      tr2_zero = so_double_correction(model)
C
      LT2MPH = LT2MP
      IT2MPH = 0
C
C------------------------------------------------------
C     Write singlet and triplet T2 amplitudes to output
C------------------------------------------------------
C
      IF ( IPRSOP. GE. 10 ) THEN
C
         CALL AROUND('singlet T2AM in HR_ERES')
         CALL OUTPUT(T2MP,1,LT2MPH,1,1,LT2MPH,1,1,LUPRI)
         IF (TRIPLET) THEN
            CALL AROUND('triplet T2AM in HR_ERES')
            CALL OUTPUT(T2MP(LT2MPH+1),1,LT2MPH,1,1,LT2MPH,1,1,LUPRI)
         END IF
C
      END IF
C
C
C------------------------------------------------------------------
C     Determine the symmetry of the result vector from the symmetry
C     of the trial vector ISYMTR, and the opperator symmtry ISYMOP.
C------------------------------------------------------------------
C
      ISYRES  = MULD2H(ISYMOP,ISYMTR)
C
C---------------------------------
C     Work space allocation no. 1.
C---------------------------------
C
      LCMO   = NLAMDT
C
      KCMO    = 1
      KEND1   = KCMO  + LCMO
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_ERES.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_ERES.1',' ',KEND1,LWORK)
C
C-------------------------------------------------------
C     Get the matrix which contains the MO coefficients.
C-------------------------------------------------------
C
#ifdef VAR_MPI
C Only master reads...
      IF (MYNUM .EQ. 0) THEN
#endif
         DTIME      = SECOND()
         CALL SO_GETMO(WORK(KCMO),LCMO,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(1)  = SOTIME(1) + DTIME
#ifdef VAR_MPI
      ENDIF
C Should probably use non-blocking collectives, where implemented
!      IF (MPI_VERSION.GE.3) THEN
!         MPI_IBCAST(WORK(KCMO),LCMO, my_MPI_INTEGER, 0,
!     &              MPI_COMM_WORLD, ierr_mpi)
!      ELSE
      count_mpi = LCMO
      CALL MPI_BCAST(WORK(KCMO), count_mpi, my_MPI_REAL8, SOP_MASTER,
     &               SOPPA_COMM_ACTIVE, ierr_mpi)

!      ENDIF
#endif
C
C---------------------------------
C     Work space allocation no. 2.
C---------------------------------
C
      LTR1E   = NT1AM(ISYMTR)
      LBTR1E  = NT1AO(ISYMTR)
      LTR1D   = NT1AM(ISYMTR)
      LRES1E  = NT1AM(ISYMTR)
      IF (DO_DEX) THEN
         LRES1D = NT1AM(ISYMTR)
      ELSE
         LRES1D = 0
      END IF
      LFOCK   = N2BST(ISYRES)
      LDENS   = N2BST(ISYMTR)
      LBTR1D  = NT1AO(ISYMTR)
      LBTJ1E  = NMATAV(ISYMTR)
      LBTJ1D  = NMATAV(ISYMTR)

      IF(DOUBLES) THEN
         IF (TRIPLET) THEN
            LTR2E = NT2SQ(ISYMTR)
         ELSE
            LTR2E = N2P2HOP(ISYMTR)
         ENDIF
         LRES2E  = N2P2HOP(ISYMTR)
         IF (DO_DEX) THEN
            LTR2D   = LTR2E
            LRES2D  = N2P2HOP(ISYMTR)
         ELSE
            LTR2D  = 0
            LRES2D = 0
         ENDIF
      ELSE
         LTR2E   = 0
         LTR2D   = 0
         LRES2E  = 0
         LRES2D  = 0
      ENDIF

      IF (singles_second) THEN
         LAIJ    = NRHFT*NRHFT
         LAAB    = NVIRT*NVIRT
      ELSE
         LAIJ    = 0
         LAAB    = 0
      ENDIF
C
      KTR1E   = KEND1
      KTR1D   = KTR1E   + LTR1E
      KTR2E   = KTR1D   + LTR1D
      KTR2D   = KTR2E   + LTR2E

      KRES1E  = KTR2D   + LTR2D
      KRES1D  = KRES1E  + LRES1E
      KRES2E  = KRES1D  + LRES1D
      KRES2D  = KRES2E  + LRES2E

      KFOCK   = KRES2D  + LRES2D
      KDENS   = KFOCK   + LFOCK
      KBTR1E  = KDENS   + LDENS
      KBTR1D  = KBTR1E  + LBTR1E
      KBTJ1E  = KBTR1D  + LBTR1D
      KBTJ1D  = KBTJ1E  + LBTJ1E

      KAIJ    = KBTJ1D  + LBTJ1D
      KAAB    = KAIJ    + LAIJ
      KEND2   = KAAB    + LAAB
#ifdef VAR_MPI
C     MPI -- Allocate timings array
      if ( loadbal_dyn ) then
         KTIMING = KEND2
         KEND2   = KTIMING + SOPPA_NINT
         CALL DZERO(WORK(KTIMING), SOPPA_NINT)
      endif
#endif
      LWORK2  = LWORK   - KEND2
C
      CALL SO_MEMMAX ('SO_ERES.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('SO_ERES.2',' ',KEND2,LWORK)
C
C----------------------------
C     Initialize AIJ and AAB.
C----------------------------
C
      IF(SINGLES_SECOND) THEN
         CALL DZERO(WORK(KAIJ),LAIJ)
         CALL DZERO(WORK(KAAB),LAAB)
      ENDIF
#ifdef VAR_MPI
C MPI -- ONLY MASTER DOES THE READING
      IF ( MYNUM .EQ. 0 ) THEN
#endif
C
C----------------------------------------------
C     Open files with trial and result vectors.
C----------------------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LURS1E,FNRS1E,LRES1E)
      CALL SO_OPEN(LURS1D,FNRS1D,LRES1E)
      IF (DO_DEX) THEN
         CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      END IF
C
C     Note: open trial-vectors also with length LRES2E
C          singlet: LRES2E = LTR2E anyway
C          triplet: We read trial-vectors into result vector memory
C                   initially, then create intermediate in
C                   TR2 memory
      IF (DOUBLES) THEN
         CALL SO_OPEN(LUTR2E,FNTR2E,LRES2E)
         CALL SO_OPEN(LURS2E,FNRS2E,LRES2E)
            CALL SO_OPEN(LURS2D,FNRS2D,LRES2E)
         IF (DO_DEX) THEN
            CALL SO_OPEN(LUTR2D,FNTR2D,LRES2D)
         ENDIF
      ENDIF
#ifdef VAR_MPI
      ENDIF
#endif
C
      IF ( IPRSOP. GE. 7
#ifdef VAR_MPI
     &      .AND. (MYNUM .EQ.0)
#endif
     &                     ) THEN ! Only printing related.
C------------------------------------------
C        Write new trial vectors to output.
C------------------------------------------
         DO 50 INEWTR = 1,NNEWTR
C----------------------------------------------------
C           Determine pointer to INEWTR trial vector.
C----------------------------------------------------
            INEW = NOLDTR + INEWTR
C
            CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,INEW)
            IF(DO_DEX)
     &               CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,INEW)
            IF (DOUBLES .AND. .NOT. SO_DOUBLE_CORRECTION(MODEL) ) THEN
               CALL SO_READ(WORK(KRES2E),LRES2E,LUTR2E,FNTR2E,INEW)
               IF(DO_DEX)
     &              CALL SO_READ(WORK(KRES2D),LRES2D,LUTR2D,FNTR2D,INEW)
            ENDIF
C
            WRITE(LUPRI,'(/,I3,A)') INEWTR,'. new trial vector'
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &           (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
            IF (DOUBLES) THEN
               WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &            (I,WORK(KRES2E+I-1),WORK(KRES2D+I-1),I=1,LRES2E)
            ENDIF
   50    CONTINUE
      END IF
C
C================================================
C     Loop over number of excitations considered.
C================================================
C
      DO 100 INEWTR = 1,NNEWTR
C
C-------------------------------------------------
C        Determine pointer to INEWTR trial vector.
C-------------------------------------------------
C
         INEW = NOLDTR + INEWTR
C
C--------------------------------------------------------------
C        Initialize RES1E, RES1D, SIGAI1, SIGAI2,
C                   SIGDA1, SIGDA2 and FOCK
C--------------------------------------------------------------
C
         CALL DZERO(WORK(KRES1E),LRES1E)
         IF (DO_DEX) CALL DZERO(WORK(KRES1D),LRES1D)

         CALL DZERO(WORK(KFOCK),LFOCK)
C
C--------------------------
C        Read trial vector.
C--------------------------
C
#ifdef VAR_MPI
         IF ( MYNUM .EQ. 0 ) THEN
#endif
            CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,INEW)
            IF (DO_DEX) THEN
               CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,INEW)
            ELSE
C     Static case: Generate D-vector from E-vector
               CALL DCOPY(LTR1E,WORK(KTR1E),1,WORK(KTR1D),1)
               IF (IDTYPE.EQ.1) CALL DSCAL(LTR1E,-1.0D0,WORK(KTR1D),1)
            END IF
            IF (DOUBLES) THEN
C              Quick fix for RPA(D) -> Later use this variable to
C              speed up first iteration for excitation energies
               IF (TR2_ZERO) THEN
                  CALL DZERO(WORK(KTR2E),LTR2E)
                  CALL DZERO(WORK(KTR2D),LTR2D)
               ENDIF

            ENDIF
#ifdef VAR_MPI
         ENDIF
C
#endif
         IF (DOUBLES) THEN
C--------------------------------------------------------------------
C        RES2E, RES2D is initialized with the D^0*Tr2 contribution
C        --- For MPI only one process must do this
C--------------------------------------------------------------------
            IF ( (.NOT.TR2_ZERO)
#ifdef VAR_MPI
     &            .AND.(MYNUM .EQ. 0 )
#endif
     &                        ) THEN
C
               IF (TRIPLET) THEN
C-----------------------------------------------------------
C                 For triplet we read the trial-vector into
C                 solution vector memory. Then we use that
C                 to create the non-symmetric intermediate on
C                 KTR2*
C----------------------------------------------------------
                  CALL SO_READ(WORK(KRES2E),LRES2E,LUTR2E,FNTR2E,INEW)
                  CALL SO_TRANTRIP(WORK(KTR2E),WORK(KRES2E),ISYMTR)
                  IF(DO_DEX) THEN
                     CALL SO_READ(WORK(KRES2D),LRES2D,LUTR2D,
     &                     FNTR2D,INEW)
                     CALL SO_TRANTRIP(WORK(KTR2D),WORK(KRES2D),ISYMTR)
                  END IF
C
                  DTIME      = SECOND()
                  CALL SO_RES_CDT(WORK(KRES2E),LRES2E,
     &                            WORK(KRES2D),LRES2D,
     &                            FOCKD,LFOCKD,ISYRES,DO_DEX,
     &                            WORK(KEND2),LWORK2)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(30) = SOTIME(30) + DTIME
C
               ELSE
C
                  CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,INEW)
                  IF(DO_DEX)
     &               CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,INEW)
C
                  DTIME      = SECOND()
                  CALL SO_RES_CD(WORK(KRES2E),LRES2E,
     &                           WORK(KRES2D),LRES2D,
     &                           WORK(KTR2E),LTR2E,WORK(KTR2D),LTR2D,
     &                           FOCKD,LFOCKD,ISYRES,
     &                           DO_DEX,WORK(KEND2),LWORK2)

                  DTIME      = SECOND()   - DTIME
                  SOTIME(30) = SOTIME(30) + DTIME
C
C   The trial-vectors are no longer needed in the original
C   basis, transform them.
C
                  CALL CCSD_TCMEPKX(WORK(KTR2E),TWO,ISYMTR)
                  IF (DO_DEX) CALL CCSD_TCMEPKX(WORK(KTR2D),TWO,ISYMTR)
               ENDIF
C
            ELSE
               CALL DZERO(WORK(KRES2E),LRES2E)
               IF (DO_DEX) CALL DZERO(WORK(KRES2D),LRES2D)
            ENDIF
         END IF
C
#ifdef VAR_MPI
C---------------------------
C Communicate trial-vectors
C---------------------------
C Rememember that the trial-vectors are contigous in memory
         LTRTOT_mpi = LTR1E + LTR1D + LTR2E + LTR2D
         CALL MPI_BCAST(WORK(KTR1E), LTRTOT_mpi, my_MPI_REAL8,
     &                  sop_master, SOPPA_COMM_ACTIVE, ierr_mpi)
C Note for future adjustment:
C It may be worthwile using non-blocking communication here
C since the next two calls only use the singles vectors.
C
C Also, depending on the communication overhead, it may be better
C to transform on host only and then send it to the slaves..?
#endif
C
C---------------------------------------------------
C        Calculate RPA-density matrices in AO basis.
C---------------------------------------------------
C
         DTIME     = SECOND()
         CALL SO_AODENS(WORK(KDENS),LDENS,WORK(KCMO),LCMO,
     &                  WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,ISYMTR,
     &                  WORK(KEND2),LWORK2)
         DTIME     = SECOND()  - DTIME
         SOTIME(6) = SOTIME(6) + DTIME
C
C--------------------------------------------
C        Backtransformation of trial vectors.
C--------------------------------------------
C
         DTIME     = SECOND()
         CALL SO_BCKTR(WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                 WORK(KBTR1E),LBTR1E,WORK(KBTR1D),LBTR1D,
     &                 WORK(KBTJ1E),LBTJ1E,WORK(KBTJ1D),LBTJ1D,
     &                 WORK(KCMO),LCMO,ISYMTR)
         DTIME     = SECOND()  - DTIME
         SOTIME(7) = SOTIME(7) + DTIME
C
C=======================================================
C        Start the loop over distributions of integrals.
C=======================================================
C
         IF (DIRECT) THEN
            NTOSYM = 1
            DTIME  = SECOND()
            IF (HERDIR) THEN
               CALL HERDI1(WORK(KEND2),LWORK2,IPRINT)
               KINDXB = KEND2
            ELSEIF(INEWTR.EQ.1) THEN
C Can we get away with doing this once?
C Because we CAN'T move KEND2 forwards in each loop cycle,
C since that'll be leaking memory
               KCCFB1 = KEND2
               KINDXB = KCCFB1 + MXPRIM*MXCONT
               KEND2  = KINDXB + (8*MXSHEL*MXCONT)/IRAT
               LWORK2 = LWORK  - KEND2

               CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                     KODPP1,KODPP2,KRDPP1,KRDPP2,KFREE,LFREE,
     &                     KEND2,WORK(KCCFB1),WORK(KINDXB),WORK(KEND2),
     &                     LWORK2,IPRINT)

               KEND2  = KFREE
               LWORK2 = LFREE
            ENDIF
            DTIME     = SECOND()  - DTIME
            SOTIME(8) = SOTIME(8) + DTIME
         ELSE
            NTOSYM = NSYM
         ENDIF
C
         ICDEL1  = 0
C
#ifdef VAR_MPI
C-----------------------------------------------------------------
C        In MPI calculations, we need to distribute the indices.
C        On first pass, we set this using the pre-sorted approach.
C-----------------------------------------------------------------
C         IF ( (nit .ne. 1 .or. inewtr .ne. 1)
         IF ( .not. (nit .eq. 1 .and. inewtr .eq. 1)
     &                        .or.  numprocs .eq. 1 ) THEN
C           After first pass, load-balancing has been done
            CONTINUE
         ELSEIF( mynum .eq. 0 )THEN
C           Master does the balancing
            call presortloadbal_parsoppa(AssignedIndices,   maxnumjobs,
     &                     work(kindxb), work(kend2), lwork2)
         ELSE
            count_mpi = maxnumjobs
C           Slaves receives the scatter, send arguments should be ignored
            call mpi_scatter( AssignedIndices, count_mpi,my_mpi_integer,
     &                        AssignedIndices, count_mpi,my_mpi_integer,
     &                        sop_master, soppa_comm_active, ierr_mpi )
         ENDIF
#endif
         DO 210 ISYMD1 = 1,NTOSYM
C
            IF (DIRECT) THEN
               IF (HERDIR) THEN
                  NTOT = MAXSHL
               ELSE
                  NTOT = MXCALL
               ENDIF
            ELSE
               NTOT = NBAS(ISYMD1)
            ENDIF
#ifdef VAR_MPI
            IF(numprocs .gt. 1) THEN
               NLOOPIDX = maxnumjobs
            ELSE
               NLOOPIDX = NTOT
               loadbal_dyn = .false.
            ENDIF
C
#endif
C
C-------------------------------------------------
C           Main loop over integral-distributions.
C-------------------------------------------------

#ifdef VAR_MPI
C------------------------------------------------------------------
C           For_parallel calculations, we have som stuff to set up.
C------------------------------------------------------------------
            DO 220 ILLLDUMMY = 1, nloopidx
               if ( numprocs .gt. 1 ) then
                  ILLL = assignedIndices(illldummy)
C                 A zero indicates that we have no more work, exit loop
                  IF (ILLL .eq. 0) exit
                  if ( loadbal_dyn ) timeini = mpi_wtime()
               else
                  ILLL = ILLLDUMMY
               endif
#else
            DO 220 ILLL = 1,NTOT
#endif
C               print *, ILLL
C------------------------------------------------
C              If direct calculate the integrals.
C------------------------------------------------
               IF (DIRECT) THEN
C
                  DTIME  = SECOND()
                  IF (HERDIR) THEN
C
                    CALL HERDI2(WORK(KEND2),LWORK2,INDEXA,ILLL,NUMDIS,
     &                          IPRINT)
C
                  ELSE
C
                     CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                           WORK(KODCL1),WORK(KODCL2),
     &                           WORK(KODBC1),WORK(KODBC2),
     &                           WORK(KRDBC1),WORK(KRDBC2),
     &                           WORK(KODPP1),WORK(KODPP2),
     &                           WORK(KRDPP1),WORK(KRDPP2),
     &                           WORK(KCCFB1),WORK(KINDXB),
     &                           WORK(KEND2),LWORK2,IPRINT)
C
                  ENDIF
                  DTIME     = SECOND()  - DTIME
                  SOTIME(9) = SOTIME(9) + DTIME
C
                  LRECNR  = ( (NBUFX(0) -1) / IRAT ) + 1
                  KRECNR  = KEND2
                  KEND2B   = KRECNR + LRECNR
                  LWORK2B  = LWORK  - KEND2B
C
                  CALL SO_MEMMAX ('SO_ERES.2B',LWORK2B)
                  IF (LWORK2 .LT. 0)
     &                CALL STOPIT('SO_ERES.2B',' ',KEND2B,LWORK)
C
               ELSE
                  NUMDIS = 1
                  KEND2B = KEND2
               ENDIF
C
C-------------------------------------------------------------------
C   Loop over number of distributions in disk.
C   In the case of ERI there are more than one distribution and IDEL2
C   loops over them and the actual index of the delta orbital IDEL is
C   then obtained from the array INDEXA. In the case of a not direct
C   calculation there is only one distribution on the disk, which
C   implies that IDEL2 is always 1 and that IDEL is systematically
C   incremented by one each time.
C--------------------------------------------------------------------
C
               DO 230 IDEL2 = 1,NUMDIS
C
                  IF (DIRECT) THEN
                     IDEL  = INDEXA(IDEL2)
                     ISYMD = ISAO(IDEL) !keeps track of current symmetry
                  ELSE
                     IDEL  = IBAS(ISYMD1) + ILLL
                     ISYMD = ISYMD1
                  ENDIF
C
                  ISYDIS = MULD2H(ISYMD,ISYMOP)
C
                  IT2DEL(IDEL) = ICDEL1
                  ICDEL1       = ICDEL1 + NT2BCD(ISYDIS)
C
C---------------------------------------------
C                 Work space allocation no. 3.
C---------------------------------------------
C
                  LXINT  = NDISAO(ISYDIS)
C
                  KXINT   = KEND2B
                  KEND3   = KXINT + LXINT
                  LWORK3  = LWORK - KEND3
C
                  CALL SO_MEMMAX ('SO_ERES.3',LWORK3)
                  IF (LWORK3 .LT. 0)
     &                CALL STOPIT('SO_ERES.3',' ',KEND3,LWORK)
C
C--------------------------------------------
C                 Read in batch of integrals.
C--------------------------------------------
C
                  DTIME      = SECOND()
                  CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWORK3,
     &                        WORK(KRECNR),DIRECT)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(10) = SOTIME(10) + DTIME
C
C
C---------------------------------------------
C                    Work space allocation no. 4.
C---------------------------------------------
C
                  ISAIJ = MULD2H(ISYMD,1)
C
                  IF (DOUBLES.OR.SINGLES_SECOND) THEN
C
                     LT2M1 = NT2BCD(ISAIJ)
                     LX2M1 = NT2BCD(MULD2H(ISYMD,ISYMTR))
                     KT2M1  = KEND3
                     KX2EM1 = KT2M1  + LT2M1
                     KX2DM1 = KX2EM1 + LX2M1
                     KEND4  = KX2DM1 + LX2M1
                     LWORK4 = LWORK  - KEND4
C
                     CALL SO_MEMMAX ('SO_ERES.4',LWORK4)
                     IF (LWORK4 .LT. 0)
     &                   CALL STOPIT('SO_ERES.4',' ',KEND4,LWORK)
C
C
C------------------------------------------------------------
C                    Construct the partially back-transformed T2
C                    MP-amplitudes.
C------------------------------------------------------------
C
                     DTIME      = SECOND()
CPi-140316 Change LT2MP -> LT2MPH
CPi I think now we have this condition two times
                     IF (NVIR(ISYMD).GT.0) THEN
                        CALL SO_T2M1(WORK(KT2M1),LT2M1,T2MP,LT2MPH,
Cend-Pi
     &                               WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                               WORK(KEND4),LWORK4)
                     ELSE
                        CALL DZERO(WORK(KT2M1),LT2M1)
                     END IF
C
C------------------------------------------------------------
C                    Construct the partially back-transformed
C                    trial vectors.
C                    These are scaled by 1/sqrt(2)
C------------------------------------------------------------
C
                     IF(DOUBLES.AND..NOT.TRIPLET.AND.
     &                  (NVIR(ISYMD).GT.0)) THEN
                        CALL SO_X2M1(WORK(KX2EM1),LX2M1,
     &                               WORK(KTR2E),LTR2E,
     &                               WORK(KCMO),LCMO,IDEL,ISYMD,
     &                               ISYMTR,WORK(KEND4),LWORK4)
                        IF (DO_DEX) THEN
                           CALL SO_X2M1(WORK(KX2DM1),LX2M1,
     &                                  WORK(KTR2D),LTR2D,
     &                                  WORK(KCMO),LCMO,IDEL,ISYMD,
     &                                  ISYMTR,WORK(KEND4),LWORK4)
                        END IF
                     ELSE
                        CALL DZERO(WORK(KX2EM1),LX2M1)
                        IF (DO_DEX) CALL DZERO(WORK(KX2DM1),LX2M1)
                     END IF

                     IF (SINGLES_SECOND) THEN
C
C---------------------------------------------------------------
C                    Add T2^(ab)_(ij)*x1(b, delta) intermediates
C                    to the backtransformed trialvectors.
C---------------------------------------------------------------
C                    (This takes care of terms (1), and (5) of the B
C                    matrix)
                        CALL SO_T2X1(WORK(KX2EM1),LX2M1,
     &                               T2MP,LT2MPH,
     &                               WORK(KBTJ1D),LBTJ1D,IDEL,ISYMD,
     &                               ISYMTR,WORK(KEND4),LWORK4)
                        IF (DO_DEX) THEN
                           CALL SO_T2X1(WORK(KX2DM1),LX2M1,
     &                                  T2MP,LT2MPH,
     &                                  WORK(KBTJ1E),LBTJ1E,IDEL,ISYMD,
     &                                  ISYMTR,WORK(KEND4),LWORK4)
                        END IF
C
C     For triplet we need to transform the intermediates
C
                        IF (TRIPLET.AND.
     &                      NVIR(MULD2H(ISYMD,ISYMTR)).GE.1) THEN
                           CALL SO_M1SHUF(WORK(KX2EM1),LX2M1,
     &                                    ISYMD,ISYMTR)
                           IF (DO_DEX) THEN
                              CALL SO_M1SHUF(WORK(KX2DM1),LX2M1,
     &                                       ISYMD,ISYMTR)
                           END IF
                        END IF

                     END IF
C
C     In the triplet case we create the partially back-transformed
C     intermediate last.
C
                     IF(DOUBLES.AND.TRIPLET.AND.
     &                  (NVIR(ISYMD).GT.0)) THEN
                        CALL SO_X2SM1(WORK(KX2EM1),LX2M1,
     &                               WORK(KTR2E),LTR2E,
     &                               WORK(KCMO),LCMO,IDEL,ISYMD,
     &                               ISYMTR,WORK(KEND4),LWORK4)
                        IF (DO_DEX) THEN
                           CALL SO_X2SM1(WORK(KX2DM1),LX2M1,
     &                                  WORK(KTR2D),LTR2D,
     &                                  WORK(KCMO),LCMO,IDEL,ISYMD,
     &                                  ISYMTR,WORK(KEND4),LWORK4)
                        END IF
                     END IF
C
                     DTIME      = SECOND()   - DTIME
                     SOTIME(12) = SOTIME(12) + DTIME

                  ELSE
                     KT2M1 = HUGE(LWORK) ! Crash if this is accessed
                     KX2EM1 = HUGE(LWORK)
                     KX2DM1 = HUGE(LWORK)
                     LT2M1 = 0
                     LX2M1 = 0
                     KEND4 = KEND3
                     LWORK4 = LWORK3
                  ENDIF
C
C---------------------------------------------
C                 Work space allocation no. 5.
C---------------------------------------------
C
                  IF (DOUBLES.OR.SINGLES_SECOND) THEN
                     LDSRHF = NDSRHF(ISYMD)
C
                     KDSRHF = KEND4
                     KEND5  = KDSRHF + LDSRHF
                     LWORK5 = LWORK  - KEND5
C
                     CALL SO_MEMMAX ('SO_ERES.5',LWORK5)
                     IF (LWORK5 .LT. 0)
     &                  CALL STOPIT('SO_ERES.5',' ',KEND5,LWORK)
C
C----------------------------------------------------------------
C                    Transform one index in the integral batch to an
C                    occupied index.
C----------------------------------------------------------------
C
                     DTIME  = SECOND()
                     ISYMLP = 1
                     CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KCMO),
     &                        ISYMLP,WORK(KEND5),LWORK5,ISYDIS)
                     DTIME      = SECOND()   - DTIME
                     SOTIME(13) = SOTIME(13) + DTIME
                  ELSE
                     LDSRHF = 0
                     KEND5 = KEND4
                     LWORK5 = LWORK4
                     KDSRHF = HUGE(LWORK)
                  END IF
C
C-------------------------------------------------------------------
C                 Calculate part of the second order density matrix.
C-------------------------------------------------------------------
C
                  DTIME      = SECOND()
                  IF ( calc_densai .AND. (  INEWTR.EQ.1 ) ) THEN
                     CALL SO_DENSAI1(DENSAI,LDENSAI,WORK(KDSRHF),LDSRHF,
     &                               WORK(KCMO),LCMO,WORK(KT2M1),LT2M1,
     &                               ISYMD,ISYDIS,WORK(KEND5),
     &                               LWORK5)

                  END IF
                  DTIME      = SECOND()   - DTIME
                  SOTIME(41) = SOTIME(41) + DTIME
C
                  IF (TRIPLET) THEN
                     CALL LOOP_BODY_TRIPLET
                  ELSE
                     CALL LOOP_BODY_SINGLET
                  ENDIF
C
  230          CONTINUE  ! End of IDEL2 loop
C
#ifdef VAR_MPI
               IF (loadbal_dyn) then
                   timefin = mpi_wtime()
                   itimeilll = ktiming + illl - 1
                   work(itimeilll) = work(itimeilll) + timefin - timeini
               ENDIF
#endif
C
  220       CONTINUE !End of ILLL loop
C
  210    CONTINUE !end of ISYMD1 loop
C
C====================================================
C        End of loop over distributions of integrals.
C====================================================
C
#ifdef VAR_MPI
C-------------------------------------------------
C        Communicate the result vectors to master.
C-------------------------------------------------
C  Note for further development:
C  In the following, non-blocking reductions could be used,
C  since the one-particle result-vector is first needed in the
C  second call, sigai and sigda in the third and fourth call,
C  and the two-particle vector is just written to file.
C  The alternative approach, is to do all these calculations
C  on the slaves, and only reduce in the end, this would save
C  communication of sigai and sigda.
C
C  Use the fact that memory are together
         lrestot = lres1e + lres1d + lres2e + lres2d
         latot   = laij + laab
C  Master use inplace operations
         if (mynum .eq. 0 ) then
C
C           Fock-matrix
            count_mpi = lfock
            call mpi_reduce( mpi_in_place, work(kfock), count_mpi,
     &                       my_MPI_REAL8, MPI_SUM, sop_master,
     &                       soppa_comm_active, ierr_mpi)
C
C           Result-vectors
            call mpi_reduce( mpi_in_place, work(kres1e), lrestot,
     &                       my_MPI_REAL8, MPI_SUM, sop_master,
     &                       soppa_comm_active, ierr_mpi)
C
            if (singles_second) then
C
C           Aij / Aab -- only calculated in first pass
               if ( inewtr .eq. 1 ) then
C               write(lupri,*) 'Aij'
                  call mpi_reduce ( mpi_in_place, work(kaij), latot,
     &                              my_MPI_REAL8, MPI_SUM, sop_master,
     &                              soppa_comm_active, ierr_mpi)
               end if
            end if
         else
C  Slaves pass the same buffer as the recieve-buffer...
C
C           Fock-matrix
            count_mpi = lfock
            call mpi_reduce( work(kfock), work(kfock), count_mpi,
     &                       my_MPI_REAL8, MPI_SUM, sop_master,
     &                       soppa_comm_active, ierr_mpi)
C
C           Result-vectors
            call mpi_reduce( work(kres1e), work(kres1e), lrestot,
     &                       my_MPI_REAL8, MPI_SUM, sop_master,
     &                       soppa_comm_active, ierr_mpi)
C
            if (singles_second) then
C
C              Aij / Aab
               if ( inewtr .eq. 1 ) then
                  call mpi_reduce ( work(kaij), work(kaij), latot,
     &                              my_MPI_REAL8, MPI_SUM, sop_master,
     &                              soppa_comm_active, ierr_mpi)
               endif
            endif
C
C  After the reductions, the slaves are done; cycle loop
            goto 100
         endif
#endif
C
C---------------------------------------------
C        Transform AO Fock matrix to MO basis.
C---------------------------------------------
C
         DTIME      = SECOND()
         CALL TRANS_FCK(WORK(KFOCK),ISYRES)
         CALL CC_FCKMO(WORK(KFOCK),WORK(KCMO),WORK(KCMO),
     &                    WORK(KEND2),LWORK2,ISYRES,1,1)
         DTIME      = SECOND()   - DTIME
         SOTIME(24) = SOTIME(24) + DTIME
C
C------------------------------------------------------------------
C        Calculate and add the RPA two-particle parts to the result
C        vectors.
C------------------------------------------------------------------
C

         DTIME      = SECOND()
         CALL SO_TWOFOCK(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                   WORK(KFOCK),LFOCK,ISYRES,DO_DEX)
         DTIME      = SECOND()   - DTIME
         SOTIME(25) = SOTIME(25) + DTIME
C
        if (singles_second) then
C-----------------------------------------------------------------
C           Calculate and add the symmetry correcting term to A in
C           eq. (44).
C-----------------------------------------------------------------
C
            DTIME      = SECOND()

            CALL SO_RES_SYM(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                      WORK(KAIJ),LAIJ,WORK(KAAB),LAAB,WORK(KTR1E),
     &                      LTR1E,WORK(KTR1D),LTR1D,ISYRES,DO_DEX)
            DTIME      = SECOND()   - DTIME
            SOTIME(20) = SOTIME(20) + DTIME
C
C---------------------------------------------------------
C           Calculate and add the Fock-term to A in eq. (40).
C---------------------------------------------------------
C
            DTIME      = SECOND()
            CALL SO_RES_FCK(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                      WORK(KTR1E),LTR1E,WORK(KTR1D),
     &                      LTR1D,FOCKD,LFOCKD,DENSIJ,LDENSIJ,DENSAB,
     &                      LDENSAB,ISYRES,ISYMTR,DO_DEX)

            DTIME      = SECOND()   - DTIME
            SOTIME(21) = SOTIME(21) + DTIME
C
         endif
C
C------------------------------------------------------------------
C        Calculate and add the RPA one-particle parts to the result
C        vectors.
C------------------------------------------------------------------
C
         DTIME      = SECOND()

         CALL SO_ONEFOCK(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,FOCKD,
     &                   LFOCKD,WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                   ISYRES,ISYMTR,DO_DEX)
         DTIME      = SECOND()   - DTIME
         SOTIME(26) = SOTIME(26) + DTIME
C
C
#ifdef VAR_MPI
C Slaves are done
         IF (MYNUM .NE. 0) GOTO 100
#endif
C
C----------------------------------------
C        Write new result vectors to file.
C----------------------------------------
C
         CALL SO_WRITE(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,INEW)
         IF (DO_DEX)
     &      CALL SO_WRITE(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,INEW)
         IF (DOUBLES) THEN
            CALL SO_WRITE(WORK(KRES2E),LRES2E,LURS2E,FNRS2E,INEW)
            IF(DO_DEX)
     &            CALL SO_WRITE(WORK(KRES2D),LRES2D,LURS2D,FNRS2D,INEW)
         ENDIF
C
C     Write zeroes to D-solution vectors
         IF (.NOT. DO_DEX) THEN
            CALL DZERO(WORK(KRES1E),LRES1E)
            CALL SO_WRITE(WORK(KRES1E),LRES1E,LURS1D,FNRS1D,INEW)
            IF (DOUBLES) THEN
               CALL DZERO(WORK(KRES2E),LRES2E)
               CALL SO_WRITE(WORK(KRES2E),LRES2E,LURS2D,FNRS2D,INEW)
            END IF
         END IF
C
  100 CONTINUE
C
C==================================
C     End of loop over excitations.
C==================================
C
#ifdef VAR_MPI
C-----------------------------------------------------------
C     Communicate the second order density matrix if needed.
C-----------------------------------------------------------
      IF ( CALC_DENSAI ) THEN
         count_mpi = LDENSAI
         IF (MYNUM.EQ.0) THEN
            CALL MPI_REDUCE( MPI_IN_PLACE, DENSAI, count_mpi,
     &                       my_MPI_REAL8, MPI_SUM, sop_master,
     &                       SOPPA_COMM_ACTIVE, ierr_mpi)
         ELSE
            CALL MPI_REDUCE( DENSAI, DENSAI, count_mpi,
     &                       my_MPI_REAL8, MPI_SUM, sop_master,
     &                       SOPPA_COMM_ACTIVE, ierr_mpi)
         ENDIF
      ENDIF
C--------------------------------------
C     Communicate the timings if needed
C--------------------------------------
      if ( loadbal_dyn ) then
         if (mynum .eq.0) then
C  Master
C --------
C Recieve the timings
            count_mpi = soppa_nint
            call mpi_reduce( mpi_in_place, work(ktiming), count_mpi,
     &                       my_MPI_REAL8, MPI_SUM, sop_master,
     &                       SOPPA_COMM_ACTIVE, ierr_mpi)
C
C Redo the loadbalancing based on the timings and distribute them
            ksorted = kend2
            ktmp    = (soppa_num_active*maxnumjobs)
            knasjob = ksorted + ktmp/irat + mod(ktmp,irat)
            kswork  = knasjob + soppa_num_active/irat + mod(ktmp,irat)
            kendf   = kswork + soppa_num_active
C
            call dynloadbal_parsoppa( AssignedIndices, maxnumjobs,
     &                                work(ktiming), soppa_nint,
     &                                work(ksorted), work(knasjob),
     &                                work(kswork) )
            !                    add empty work-space...
            !                    currently it can start at kend2
         else
C   Slave
C  ---------
C  Send the timings
            count_mpi = soppa_nint
            call mpi_reduce( work(ktiming), work(ktiming), count_mpi,
     &                       my_MPI_REAL8, MPI_SUM, sop_master,
     &                       SOPPA_COMM_ACTIVE, ierr_mpi)
C  Recieve the new indices
            count_mpi = maxnumjobs
            call mpi_scatter( AssignedIndices,count_mpi,my_mpi_integer,
     &                        AssignedIndices,count_mpi,my_mpi_integer,
     &                        sop_master, soppa_comm_active, ierr_mpi)
         endif
      endif
C Slaves are done
      IF (MYNUM .NE. 0) THEN
         IF ( CALC_DENSAI ) SOP_MP2AI_DONE = .TRUE.
         RETURN
      END IF
#endif
C----------------------------------------------------------------
C     Calculate the last part of the second order density matrix.
C----------------------------------------------------------------
C
      DTIME      = SECOND()
      IF ( CALC_DENSAI )  THEN
         CALL SO_DENSAI2(DENSAI,LDENSAI,FOCKD,LFOCKD)
         SOP_MP2AI_DONE = .TRUE.
      END IF
C
      DTIME      = SECOND()   - DTIME
      SOTIME(41) = SOTIME(41) + DTIME
C
      IF ( IPRSOP .GE. 7 ) THEN
C------------------------------------------
C        Write new resultvectors to output.
C------------------------------------------
         DO 400 INEWTR = 1,NNEWTR
            INEW = NOLDTR + INEWTR
            WRITE(LUPRI,'(/,I3,A)') INEWTR,
     &                '. new E[2] linear transformed trial vector'
            CALL SO_READ(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,INEW)
            CALL SO_READ(WORK(KRES1D),LRES1E,LURS1D,FNRS1D,INEW)
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &           (I,WORK(KRES1E+I-1),WORK(KRES1D+I-1),I=1,LRES1E)
            IF (DOUBLES) THEN
               CALL SO_READ(WORK(KRES2E),LRES2E,LURS2E,FNRS2E,INEW)
C               IF (DO_DEX)
               CALL SO_READ(WORK(KRES2D),LRES2E,LURS2D,FNRS2D,INEW)
               WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &             (I,WORK(KRES2E+I-1),WORK(KRES2D+I-1),I=1,LRES2E)
            ENDIF
  400    CONTINUE
C
      END IF
C
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
         CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
      IF (DO_DEX) THEN
         CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      END IF
C
      IF (DOUBLES) THEN
         CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
         CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
            CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
         IF (DO_DEX) THEN
            CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
         ENDIF
      ENDIF
C
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_ERES')
C
      RETURN
C
C------------------------------------------------------------------------
C  The body of the above loop over integral distributions has been moved
C  into these internal subroutines.
C------------------------------------------------------------------------
C
      CONTAINS
C  Be aware that variables fall through from the outer scope!
C

         SUBROUTINE LOOP_BODY_SINGLET
C
         DOUBLE PRECISION DTIME
C
C
C----------------------------------------------
C        Calculate the AO-Fock matrix.
C----------------------------------------------
C
         DTIME      = SECOND()
         CALL SO_AOFOCK1(WORK(KXINT),WORK(KDENS),WORK(KFOCK),
     &                  WORK(KEND5),LWORK5,
     &                  IDEL,ISYMD,
     &                  ISYMTR)
         DTIME      = SECOND()   - DTIME
         SOTIME(11) = SOTIME(11) + DTIME
         IF (singles_second) THEN
C
C----------------------------------------------------------------------
C           Calculate part of the result vectors RES1E and RES1D,
C           specifically the first and the second term in eqs. (34,35).
C           Also calculate Aij and Aab in eqs. (43,44).
C----------------------------------------------------------------------
C
            DTIME      = SECOND()
            CALL SO_RES_A(WORK(KRES1E),LRES1E,
     &                    WORK(KRES1D),LRES1D,
     &                    WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                    WORK(KDSRHF),LDSRHF,WORK(KCMO),LCMO,
     &                    WORK(KT2M1),LT2M1,WORK(KAIJ),LAIJ,
     &                    WORK(KAAB),LAAB,INEWTR,ISYMD,ISYDIS,
     &                    ISYRES,ISYMTR,DO_DEX,
     &                    WORK(KEND5),LWORK5)
            DTIME      = SECOND()   - DTIME
            SOTIME(14) = SOTIME(14) + DTIME
         ENDIF
C
         IF (DOUBLES.OR.SINGLES_SECOND) THEN
C
C-------------------------------------------------------------------
C        Calculate the part of the result vectors RES1E and
C        RES1D which originate from the C matrices. See
C                 eqs. (72) and (73).
C-------------------------------------------------------------------
CRF      The X2M1 vectors is now an intermediate, which contain
CRF      both the 2-particle trial vector AND the T2*x1 intermediate,
CRF              ~ ~                               ~
CRF      so both C*x2 and terms (1) and (5) of B*x1 is calculated here.
C
            DTIME      = SECOND()
            CALL SO_RES_TCB(WORK(KRES1E),LRES1E,
     &                      WORK(KRES1D),LRES1D,
     &                      WORK(KX2EM1),LX2M1,
     &                      WORK(KX2DM1),LX2M1,
     &                      WORK(KDSRHF),LDSRHF,
     &                      WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                      ISYMTR,DO_DEX,WORK(KEND5),LWORK5)
            DTIME      = SECOND()   - DTIME
            SOTIME(29) = SOTIME(29) + DTIME
         END IF
C
C----------------------------------------------------------------------
C           Construct C-contribution to 2p2h result vectors RES2E
C                 and RES2D.
C----------------------------------------------------------------------
C
         IF(DOUBLES) THEN
            DTIME      = SECOND()
            CALL SO_RES_CB(WORK(KRES2E),LRES2E,
     &                    WORK(KRES2D),LRES2D,
     &                    WORK(KDSRHF),LDSRHF,
     &                    WORK(KBTR1E),LBTR1E,
     &                    WORK(KBTR1D),LBTR1D,
     &                    WORK(KBTJ1E),LBTJ1E,
     &                    WORK(KBTJ1D),LBTJ1D,WORK(KCMO),LCMO,
     &                    IDEL,ISYMD,ISYDIS,ISYMTR,DO_DEX,
     &                    WORK(KEND5),LWORK5)
            DTIME      = SECOND()   - DTIME
            SOTIME(15) = SOTIME(15) + DTIME
C
         ENDIF

         IF (SINGLES_SECOND) THEN
            DTIME   = SECOND()
            ISYDIS2 = MULD2H(ISYDIS,ISYMTR)
            KEND6 = KEND4 + NDSRHF(MULD2H(ISYMD,ISYMTR))
            LWORK6 = LWORK - KEND6
C                                     ~
C           Calculate ( alpha, beta | j delta)
            CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KBTR1D),
     &                  ISYMTR,WORK(KEND5),LWORK5,ISYDIS)
C           Calculate terms 2 and 6 of the B-matrix, by
C                     ~
C    (2)    - ( c a | j delta ) * T2M1( ci, j)
C                     ~
C    (6)      ( k i | j delta ) * T2M1( ak, j)

            CALL SO_RES_B26 ( WORK(KRES1E), LRES1E, WORK(KT2M1), LT2M1,
     &                        WORK(KDSRHF),WORK(KCMO),LCMO,IDEL,ISYMD,
     &                        ISYDIS2,ISYMTR,WORK(KEND6),LWORK6)
            IF (DO_DEX) THEN
C
C           Same for the D-part
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KBTR1E),
     &                     ISYMTR,WORK(KEND5),LWORK5,ISYDIS)
               CALL SO_RES_B26 ( WORK(KRES1D), LRES1D,
     &                           WORK(KT2M1), LT2M1,
     &                           WORK(KDSRHF),WORK(KCMO),LCMO,
     &                           IDEL,ISYMD,
     &                           ISYDIS2,ISYMTR,WORK(KEND6),LWORK6)
            END IF
            DTIME      = SECOND()   - DTIME
            SOTIME(16) = SOTIME(16) + DTIME
         ENDIF
         END SUBROUTINE
C
         SUBROUTINE LOOP_BODY_TRIPLET
C
C----------------------------------------------
C        Calculate the AO-Fock matrix.
C----------------------------------------------
C
         DTIME      = SECOND()
         CALL SO_AOFOCK3(WORK(KXINT),WORK(KDENS),
     &                   WORK(KFOCK),
     &                   IDEL,ISYMD,ISYMTR)
         DTIME      = SECOND()   - DTIME
         SOTIME(11) = SOTIME(11) + DTIME
C
         IF (singles_second) THEN
C
CRF         A^(2) Doesn't mix spins - Same for singlet and triplet
C----------------------------------------------------------------------
C           Calculate part of the result vectors RES1E and RES1D,
C           specifically the first and the second term in eqs.
C           (34,35). Also calculate Aij and Aab in eqs. (43,44).
C----------------------------------------------------------------------
C
            DTIME      = SECOND()
            CALL SO_RES_A(WORK(KRES1E),LRES1E,
     &                    WORK(KRES1D),LRES1D,
     &                    WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                    WORK(KDSRHF),LDSRHF,WORK(KCMO),LCMO,
     &                    WORK(KT2M1),LT2M1,WORK(KAIJ),LAIJ,
     &                    WORK(KAAB),LAAB,INEWTR,ISYMD,ISYDIS,
     &                    ISYRES,ISYMTR,DO_DEX,
     &                    WORK(KEND5),LWORK5)
            DTIME      = SECOND()   - DTIME
            SOTIME(14) = SOTIME(14) + DTIME
C
C--------------------------------------------------------------------
C           Transform the partially back-transformed T2 MP-amplitudes
C           Currently they are stored as
C              t2m1(ai,j) = 2t(ai,j) - t(aj,i)
C           In the following we need to have
C              t2m1(ai,j) = -t(aj,i)
C--------------------------------------------------------------------
C
            DTIME      = SECOND()
            CALL SO_M1SHUF(WORK(KT2M1),LT2M1,ISYMD,1)
            DTIME      = SECOND()   - DTIME
            SOTIME(12) = SOTIME(12) + DTIME
         END IF
         IF (DOUBLES.OR.SINGLES_SECOND) THEN
C
C-------------------------------------------------------------------
C        Calculate the part of the result vectors RES1E and
C        RES1D which originate from the C matrices. See
C                 eqs. (72) and (73).
C-------------------------------------------------------------------
CRF      The X2M1 vectors are now intermediates, which contain
CRF      the the T2*x1 intermediate,
CRF      so terms (1) and (5) of B*x1 are calculated here.
C
            DTIME      = SECOND()
            CALL SO_RES_TCB(WORK(KRES1E),LRES1E,
     &                      WORK(KRES1D),LRES1D,
     &                      WORK(KX2EM1),LX2M1,
     &                      WORK(KX2DM1),LX2M1,
     &                      WORK(KDSRHF),LDSRHF,
     &                      WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                      ISYMTR,DO_DEX,
     &                      WORK(KEND5),LWORK5)
            DTIME      = SECOND()   - DTIME
            SOTIME(29) = SOTIME(29) + DTIME

         END IF
C
         IF (DOUBLES) THEN
C
C----------------------------------------------------------------------
C           Construct C-contribution to 2p2h result vectors RES2E
C                 and RES2D.
C----------------------------------------------------------------------
C
            DTIME      = SECOND()
            CALL SO_RES_CBT(WORK(KRES2E),LRES2E,
     &                      WORK(KRES2D),LRES2D,
     &                      WORK(KDSRHF),LDSRHF,
     &                      WORK(KBTR1E),LBTR1E,
     &                      WORK(KBTR1D),LBTR1D,
     &                      WORK(KBTJ1E),LBTJ1E,
     &                      WORK(KBTJ1D),LBTJ1D,WORK(KCMO),LCMO,
     &                      IDEL,ISYMD,ISYDIS,ISYMTR,DO_DEX,
     &                      WORK(KEND5),LWORK5)
            DTIME      = SECOND()   - DTIME
            SOTIME(15) = SOTIME(15) + DTIME
C
         ENDIF
         IF(SINGLES_SECOND) THEN
            DTIME   = SECOND()
            ISYDIS2 = MULD2H(ISYDIS,ISYMTR)
            KEND6 = KEND4 + NDSRHF(MULD2H(ISYMD,ISYMTR))
            LWORK6 = LWORK - KEND6
C                                     ~
C           Calculate ( alpha, beta | j delta)
            CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KBTR1D),
     &                  ISYMTR,WORK(KEND5),LWORK5,ISYDIS)
C           Calculate terms 2 and 6 of the B-matrix, by
C                     ~
C    (2)    - ( c a | j delta ) * T2M1( ci, j)
C                     ~
C    (6)      ( k i | j delta ) * T2M1( ak, j)
C
            CALL SO_RES_B26 ( WORK(KRES1E), LRES1E, WORK(KT2M1), LT2M1,
     &                        WORK(KDSRHF),WORK(KCMO),LCMO,IDEL,ISYMD,
     &                        ISYDIS2,ISYMTR,WORK(KEND6),LWORK6)
C
            IF (DO_DEX) THEN
C           Same for the D-part
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KBTR1E),
     &                    ISYMTR,WORK(KEND5),LWORK5,ISYDIS)
               CALL SO_RES_B26 (WORK(KRES1D), LRES1D,
     &                          WORK(KT2M1), LT2M1,
     &                          WORK(KDSRHF),WORK(KCMO),LCMO,IDEL,ISYMD,
     &                          ISYDIS2,ISYMTR,WORK(KEND6),LWORK6)
            END IF
            DTIME      = SECOND()   - DTIME
            SOTIME(16) = SOTIME(16) + DTIME
         ENDIF

         END SUBROUTINE

         SUBROUTINE TRANS_FCK(FOCK,ISYMFCK)
            DOUBLE PRECISION, INTENT(INOUT) :: FOCK(*)
            INTEGER, INTENT(IN) :: ISYMFCK

            INTEGER :: IA, IB, ISIZE, NUMA, NUMB
            INTEGER :: ISYMA, ISYMB, IOFF, IOFF1, IOFF2
            DOUBLE PRECISION :: TMP

            IF ( ISYMFCK .EQ. 1) THEN
               DO ISYMA = 1, NSYM
                  IOFF = IAODIS(ISYMA,ISYMA)
                  DO IA = 2, NBAS(ISYMA)
                     DO IB = 1, IA-1
                        IPOS1 = (IA-1)*NBAS(ISYMA)+IB+IOFF
                        IPOS2 = (IB-1)*NBAS(ISYMA)+IA+IOFF
                        TMP = FOCK(IPOS1)
                        FOCK(IPOS1) = FOCK(IPOS2)
                        FOCK(IPOS2) = TMP
                     END DO
                  END DO
               END DO
            ELSE
               DO ISYMA = 1, NSYM
                  ISYMB = MULD2H(ISYMFCK,ISYMA)
                  IF (ISYMB .GT.ISYMA) CYCLE
                  NUMA = NBAS(ISYMA)
                  NUMB = NBAS(ISYMB)
                  IOFF1 = IAODIS(ISYMA,ISYMB)
                  IOFF2 = IAODIS(ISYMB,ISYMA)
                  ISIZE = NUMA*NUMB
                  DO IB = 1, NUMB
                     DO IA = 1, NBAS(ISYMA)
                        IDX1 = IOFF1 + NUMA*(IB-1) +IA
                        IDX2 = IOFF2 + NUMB*(IA-1) +IB
                        TMP = FOCK(IDX1)
                        FOCK(IDX1) = FOCK(IDX2)
                        FOCK(IDX2) = TMP
                     END DO
                  END DO
               END DO
            ENDIF

         END SUBROUTINE

         SUBROUTINE SO_AOFOCK3(XINT,DENSIT,FOCK,
     &                         IDEL,ISYMD,ISYDEN)

#include "symsq.h"

            INTEGER,INTENT(IN) ::  IDEL, ISYMD,ISYDEN
            DOUBLE PRECISION,INTENT(IN)    :: XINT(*), DENSIT(*)
            DOUBLE PRECISION,INTENT(INOUT) :: FOCK(*)

            INTEGER :: ISYMBG, ISYMAB,ISYMG,ISYDIS,ISYMB,ISYMA
            INTEGER :: IG, IA, IB
            INTEGER :: KOFFINT, KOFFINT1, KOFFDEN, KOFFOUT
            INTEGER :: NOFFB, NOFFA
            DOUBLE PRECISION :: TMP

            ISYMBG = ISYDEN
            ISYDIS = ISYMD
            ISYMA  = MULD2H(ISYDIS,ISYMBG)
            KOFFOUT = IAODIS(ISYMA,ISYMD) + NBAS(ISYMA)*
     &                (IDEL- IBAS(ISYMD) - 1)

            DO ISYMG = 1, NSYM
               ISYMB = MULD2H(ISYMBG,ISYMG)
               ISYMAB = MULD2H(ISYMA,ISYMB)
               KOFFINT1 = IDSAOG(ISYMG,ISYDIS) + IAODPK(ISYMA,ISYMB)

               IF(ISYMAB.EQ.1) THEN ! integrals stored as triangle
                  DO IG = 1, NBAS(ISYMG)
                     KOFFINT = KOFFINT1 + (IG - 1) *NNBST(ISYMAB)
                     KOFFDEN = IAODIS(ISYMB,ISYMG) + NBAS(ISYMB)*(IG-1)
                     ! alpha =< beta :
                     DO IB = 1, NBAS(ISYMB)
                        NOFFB = IB*(IB-1)/2 + KOFFINT
                        DO IA = 1, IB
                           FOCK(KOFFOUT+IA) = -XINT(NOFFB+IA)
     &                                        *DENSIT(KOFFDEN+IB)
     &                                        +FOCK(KOFFOUT+IA)
                        END DO
                     END DO

                     ! alpha > beta
                     DO IA = 2, NBAS(ISYMA)
                        NOFFA = IA*(IA-1)/2 + KOFFINT
                        TMP = 0.0D0
                        DO IB = 1, IA - 1
                           TMP = XINT(NOFFA+IB)*DENSIT(KOFFDEN+IB) + TMP
                        END DO
                        FOCK(KOFFOUT+IA) = FOCK(KOFFOUT+IA) - TMP
                     END DO
                 END DO ! LOOP IG

               ELSEIF(ISYMA .LT.ISYMB) THEN
                  ! Stored as alpha, beta
                  DO IG = 1, NBAS(ISYMG)
                     KOFFINT = KOFFINT1 + (IG - 1) *NNBST(ISYMAB)
                     KOFFDEN = IAODIS(ISYMB,ISYMG) + NBAS(ISYMB)*(IG-1)

                     ! Loop B first
                     DO IB = 1, NBAS(ISYMB)
                        NOFFB = NBAS(ISYMA)*(IB-1) + KOFFINT
                        DO IA = 1, NBAS(ISYMA)
                           FOCK(KOFFOUT+IA) = -XINT(NOFFB+IA)
     &                                        *DENSIT(KOFFDEN+IB)
     &                                        +FOCK(KOFFOUT+IA)
                        END DO
                     END DO
                  END DO
               ELSE ! ISYMA .GT. ISYMB
                  ! Stored as beta, alpha
                  DO IG = 1, NBAS(ISYMG)
                     KOFFINT = KOFFINT1 + (IG - 1) *NNBST(ISYMAB)
                     KOFFDEN = IAODIS(ISYMB,ISYMG) + NBAS(ISYMB)*(IG-1)
                     ! LOOP A first
                     DO IA = 1, NBAS(ISYMA)
                        NOFFA = NBAS(ISYMB)*(IA-1) + KOFFINT
                        TMP = 0.0D0
                        DO IB = 1, NBAS(ISYMB)
                           TMP = XINT(NOFFA+IB)*DENSIT(KOFFDEN+IB) + TMP
                        END DO
                        FOCK(KOFFOUT+IA) = FOCK(KOFFOUT+IA) - TMP
                     END DO
                  END DO ! IG
               END IF
            END DO ! ISYMG

         END SUBROUTINE
C
         SUBROUTINE SO_SYM_DENS(DENS_SQ,DENS_SYM,ISYDENS)
            !                           ~
            !  Form the symmetric array D(alpha,beta) from
            !  the array D(alpha, beta) and store it in packed form
            !  (alpha =< beta).
            !  ~
            !  D(alpha,beta) = D(alpha,beta)+ D(beta,alpha)
            !  ~
            !  D(alpha,alpha) = D(alpha,alpha)
            !
            !  D has symmetry ISYDENS
#include "symsq.h"
            DOUBLE PRECISION, INTENT(IN) :: DENS_SQ(*)
            DOUBLE PRECISION, INTENT(OUT):: DENS_SYM(*)
            INTEGER, INTENT(IN)          :: ISYDENS

            INTEGER :: ISYMA, ISYMB, IA, IB
            INTEGER :: NUMA, NUMB, IOFFA, IOFFB, IOFFPK, IOFFSQ
            INTEGER :: IOFF1, IOFF2,IOFFSQ1,IOFFSQ2

            IF (ISYDENS.EQ.1) THEN
               ! Totally symmetric case, Triangular storage
               DO ISYMA = 1, NSYM
                  IOFFPK = IAODPK(ISYMA,ISYMA)
                  IOFFSQ = IAODIS(ISYMA,ISYMA)
                  NUMA = NBAS(ISYMA)
                  NUMB = NBAS(ISYMA)
                  DO IB = 1, NUMB
                     IOFF1 = IB*(IB-1)/2 + IOFFPK
                     IOFFB = (IB-1)*NUMA + IOFFSQ
                     DO IA = 1, IB-1
                        IOFFA = (IA-1)*NUMB + IOFFSQ
                        DENS_SYM(IOFF1+IA) = DENS_SQ(IOFFB+IA) +
     &                                       DENS_SQ(IOFFA+IB)
                     END DO
                     DENS_SYM(IOFF1+IB) = DENS_SQ(IOFFB+IB)
                  END DO
               END DO
            ELSE ! Else we deal with rectangular offdiagonal blocks
               DO ISYMB = 1, NSYM
                  ISYMA = MULD2H(ISYMB,ISYDENS)
                  IF (ISYMA.GT.ISYMB) CYCLE
                  IOFFPK=IAODPK(ISYMA,ISYMB)
                  IOFFSQ1 = IAODIS(ISYMB,ISYMA)
                  IOFFSQ2 = IAODIS(ISYMA,ISYMB)
                  NUMA = NBAS(ISYMA)
                  NUMB = NBAS(ISYMB)
                  DO IB = 1, NUMB
                     IOFF1 = (IB-1)*NUMA + IOFFPK
                     IOFFB = (IB-1)*NUMA + IOFFSQ2
                     DO IA = 1, NUMA
                        IOFFA = (IA-1)*NUMB + IOFFSQ1
                        DENS_SYM(IOFF1+IA) = DENS_SQ(IOFFB+IA) +
     &                                       DENS_SQ(IOFFA+IB)
                     END DO
                  END DO
               END DO
            ENDIF

         END SUBROUTINE

         SUBROUTINE SO_AOFOCK1(XINT,DENSIT,FOCK,WORK,LWORK,
     &                         IDEL,ISYMD,ISYDEN)
            !  Calculates here
            !  2 ( alpha, beta | gamma ; delta) * D( alpha, beta) =>
            !                 F ( gamma ; delta)
            !
            !  and in call to SO_AOFOCK3 :
            !  - ( alpha, beta | gamma ; delta) * D( beta, gamma) =>
            !                 F (alpha ; delta)
            !
            DOUBLE PRECISION, PARAMETER  :: TWO = 2.0D0, ONE =1.0D0

            INTEGER,INTENT(IN) ::  IDEL, ISYMD,ISYDEN, LWORK
            DOUBLE PRECISION,INTENT(IN)    :: XINT(*), DENSIT(*)
            DOUBLE PRECISION,INTENT(INOUT) :: FOCK(*), WORK(LWORK)

            INTEGER :: ISYDIS, ISYMG
            INTEGER :: NUMG, NAB, KGAM, KOUT

            ! First do the singlet-only term
            ISYDIS = ISYMD
            ISYMG  = MULD2H(ISYDIS,ISYDEN)
            NUMG = NBAS(ISYMG)
            IF (NUMG .GE. 1) THEN
               !Symmetrice AO-DENSITY and compress to packed form (could be done outside)
               CALL SO_SYM_DENS(DENSIT,WORK,ISYDEN)

               ! Position of (alpha =< beta | gamma; delta ) block for this
               ! symmetry of gamma
               KGAM = IDSAOG(ISYMG,ISYDIS) + 1
               ! Position of the delta'th colunm in output
               KOUT = IAODIS(ISYMG,ISYMD) +
     &                NBAS(ISYMG)*(IDEL-IBAS(ISYMD)-1) + 1
               NAB  = MAX( 1, NNBST(ISYDEN) )

               CALL DGEMV('T',NNBST(ISYDEN),NUMG,TWO,XINT(KGAM),NAB,
     &                    WORK,1,ONE,FOCK(KOUT),1)
            END IF

            ! Do the term also appearing in the triplet case
            CALL SO_AOFOCK3(XINT,DENSIT,FOCK,IDEL,ISYMD,ISYDEN)


         END SUBROUTINE

      END SUBROUTINE

C
#ifdef VAR_MPI
C
C
      subroutine dynloadbal_parsoppa( localIndices, maxnumjobs,
     &                                timings, ltimings,
     &                                sortedindices, numassignjobs,
     &                                sumofwork)
C    Dynamic Load Balancing for the parallel SOPPA calculations (and parallel RPA).
C    The routine assumed that the MPI processes that does the load balancing is the master. If any other routine enters, there will be issues with the updated ILLL indices in the AssignedIndices array.  Right now, any slave that enters is immediately evicted from the routine.
C
C The routine takes an array of timings as input. The timings are creatd on the fly by every parallel process and the index corresponds to the ILLL index in the parallel calculation. The time is given in integers from calling system_clock. The timings are not given in any human-readable form. This choice was made so that one can reuse the subroutine getallocsize, which relies on an integer input array.
C The work is resorted once the actual time associated with an ILLL index is known and the indices are rebalanced once more and sent to all slaves.
C The improvement over the presorted balancing scheme is expected to be very small, but gets better the larger the basis set.
C
C F.Beyer Oct. 2014.

      use so_parutils, only : soppa_comm_active, soppa_num_active,
     &                        soppa_nint, my_mpi_integer, sop_master

      use so_info, only: sop_dp

      implicit none
#include "mpif.h"
#include "maxorb.h"
#include "distcl.h"
#include "priunit.h"

      ! Dummy parameters
      double precision, dimension(ltimings), intent(inout) :: timings
      integer, intent(inout) :: localIndices(maxnumjobs)
      integer  :: sortedindices(maxnumjobs,soppa_num_active)
      integer  :: numassignjobs(soppa_num_active)
      real(sop_dp) :: sumofwork(soppa_num_active)
      integer  :: ltimings, maxnumjobs

      ! Bookkeeping
      integer :: maxrows, maxcols, numrecipients
      integer :: getnumjobs, targetID, myid
      integer :: colindex, rowindex, assignILLL, col,i
      integer :: ntot

      integer(mpi_integer_kind) :: ierr_mpi, maxnum_mpi

      ntot = soppa_nint
      ! In case the amount of work is smaller than
      ! the number of MPI processes...
      maxcols = soppa_num_active

      ! Explicitly set the maximum number of jobs a single MPI process can be allocated.
      maxrows = maxnumjobs

      sortedindices(:,:) = 0
      numassignjobs(:) = 0
      sumofwork(:) = 0.0D0


      DO i=1, ntot
              ! Find the largest chunk of available work
              assignILLL = maxloc(timings,dim=1, mask=timings.ge.0.0D0)

              ! Find the laziest slave
              colindex = minloc(sumofwork, dim=1)

              ! Write the ILLL index to the correct row in the sortedmatrix
              numassignjobs(colindex) = numassignjobs(colindex) + 1
              rowindex = numassignjobs(colindex)
              sortedindices(rowindex, colindex) = assignILLL

              ! Update the slave's expected work/walltime and
              ! remove the work from the timings array
              sumofwork(colindex) = sumofwork(colindex)
     &                             +timings(assignILLL)
              timings(assignILLL) = -1.0D0

      ENDDO

C     PRINT HOW WORK IS DISTRIBUTED
      write (LUPRI,'(a)') 'AOSOPPA Work distribution'
      write (LUPRI,'(a)') 'NODE       Expected time '
      do i = 1, maxcols
           write (LUPRI, '(i5,f20.5)') i, sumofwork(i)
      enddo
      ! Send the info to every slave that does computation (some slaves might be stalled in the polling barrier in case there are too many MPI processes compared to the number of tasks).
      call mpi_scatter( sortedindices, maxnumjobs, my_mpi_integer,
     &                  localIndices,  maxnumjobs, my_mpi_integer,
     &                  sop_master, soppa_comm_active, ierr_mpi)

      ! Update the master's own array of ILLL indices
      ! should be done by the scatter

      return
      end subroutine






      subroutine presortloadbal_parsoppa(localIndices, maxnumjobs,
     &           indexb,
     &           work, lwork)
C     This subroutine load balances parallel SOPPA/RPA
C     calculation.
C
C     The routine makes a best guess of loadbalancing by giving
C     every MPI process an equal number of distributions to handle.
C     Testing shows this is a better first guess that giving every
C     MPI process an equal number of ILLL indices to work with
C     (since there is a different number of distributions
C     associated with every ILLL index).
C
C     F.Beyer Oct. 2014.

      use so_parutils, only: soppa_comm_active, soppa_nint, sop_master,
     &                       soppa_num_active, my_mpi_integer

      implicit none
C#include "implicit.h"
#include "priunit.h"
#include "mpif.h"
#include "maxorb.h"
#include "iratdef.h"
C fetch HERDIR
#include "ccsdinp.h"
C     Dummy parameters
      integer, intent(in) :: lwork, maxnumjobs
      ! Declare work as integer, to avoid "irats" later
      integer :: work(irat*lwork) !intent in
      integer :: indexb(*)
      integer, dimension(maxnumjobs), intent(out) :: localIndices

C     Pre-sorting load balancing variables
      integer :: getnumjobs, getindices, numprocs
      INTEGER(MPI_INTEGER_KIND) :: ierr_mpi
C      integer :: presortarray, finalsorted
      integer :: col, ntot
      integer :: kpresortarray, kend, kfinalsorted, ktmp

      ntot = soppa_nint

      ! numprocs = soppa_num_active !! not consistent with other use
!      call mpi_comm_size(soppa_comm_active, numprocs, ierr_mpi)
      numprocs = soppa_num_active


      IF (.NOT. HERDIR) THEN
C     ERI code

C     FIND THE AMOUNT OF WORK ASSOCIATED WITH EVERY AO INDEX
         kpresortarray = 1
         kend = kpresortarray + 2* ntot
         call presortaodist(ntot, indexb, work(kpresortarray) )


C     CREATE THE FINAL MATRIX OF PRE-SORTED AO INDICES
         kfinalsorted = kend
         ktmp = kfinalsorted + maxnumjobs * numprocs
         kend = ktmp + 2*numprocs
         if( kend .gt. lwork ) then
            call quit('Insufficient memory in presorted loadbalancer!')
         endif
         call partitionAOindices(ntot, maxnumjobs, numprocs,
     &              work(kpresortarray), work(kfinalsorted),
     &              work(ktmp) )
      else
C     HERDIR code
         kfinalsorted = 1
         kend = kfinalsorted + maxnumjobs * numprocs
         if( kend .gt. lwork ) then
            call quit('Insufficient memory in presorted loadbalancer!')
         endif
         call herdir_presort( work(kfinalsorted), maxnumjobs)
      endif


C     TRANSFER AO INDICES TO SLAVES
C     Use scatter
      call mpi_scatter( work(kfinalsorted), maxnumjobs, my_mpi_integer,
     &                  localIndices,       maxnumjobs, my_mpi_integer,
     &                  sop_master, soppa_comm_active, ierr_mpi )

      return
      end subroutine


      subroutine herdir_presort ( sorted, maxnumjobs )
C     This is a subroutine for doing the initial distribution of the
C     AO-indices in a parallel SOPPA calculation.
C     This subroutine works with Hermite (HERDIR) integral generator,
C     as the code by F. Beyer only supports using ERI
C
      use so_parutils, only: soppa_num_active, soppa_nint
      implicit none
      integer, intent(out) :: sorted( maxnumjobs, soppa_num_active )
      integer, intent(in)  :: maxnumjobs
      integer              :: intnum, inext, icol, inum

      sorted = 0
      inext = 0
      icol = 1
C     For now just do it the stupid way ( round robin )

C     initial implementation, is this what chokes ifort?
C      do intnum = 1, soppa_nint
C         inext = inext + 1
C         if ( inext .gt. soppa_num_active ) then
C             inext = 1
C             icol = icol + 1
C         endif
C         sorted ( icol, inext ) = intnum
C      enddo

C     Distribute first an even number
      inum = soppa_nint/soppa_num_active
      do icol = 1, soppa_num_active
         do inext = 1, inum
            sorted( inext, icol ) = inext + inum*(icol-1)
         enddo
      enddo
C     Distribute the remainder
      do icol = 1, mod(soppa_nint, soppa_num_active)
         sorted ( inum+1, soppa_num_active-icol+1) =
     &            inum*soppa_num_active + icol
      enddo
C     Debug print
      do inext = 1, soppa_num_active
         print '(10i3)', sorted(:,inext)
      enddo
      return
      end subroutine



      subroutine pollingbarrier(pollinginterval)
C     This subroutine is implemented in case there is ever a need to
C     remove certain processes from a calculation.
C
C     The processes will repeatedly poll until they receive a non-blocking
C     send from the master with a logical value equal to .true.
C
C     F.Beyer Oct. 2014

#ifdef VAR_MPI
      use so_parutils, only: my_mpi_logical
#endif

#include "implicit.h"
#include "mpif.h"

      integer, dimension(MPI_STATUS_SIZE) :: mpistatus
      integer, intent(in)                 :: pollinginterval !in milliseconds
      logical                             :: flag, exitbarrier
      volatile                            :: exitbarrier
      INTEGER(MPI_INTEGER_KIND) :: ierr_mpi, myid, request

      call mpi_comm_rank(mpi_comm_world, myid, ierr_mpi)
      print *, 'I am in the polling barrier ', myid

      ! Initiate the polling variables and ask for the non-blocking update
      exitbarrier = .false.
      call mpi_irecv(exitbarrier, 1, my_mpi_logical, 0, myid,
     &               mpi_comm_world, request, ierr_mpi)

      ! Polling barrier, the process will cycle repeatedly until released.
130   continue
#ifdef VAR_IFORT
      call sleepqq(pollinginterval)
#endif
      call mpi_test(request, flag, mpistatus, ierr_mpi)
      if (.not. flag) goto 130
      ! Warning, seems that we'll go into an infinite loop,
      ! if exitbarrier is ever sent as .false.
      if (.not.exitbarrier) then
         goto 130 ! Cycle to the top of the barrier.
      else
         return
      endif

      return
      end subroutine

C /* deck presortaodist*/
      subroutine presortaodist(Nindex, indxbt, outlist)
C A subroutine associated with the atomic integral parallel RPA/SOPPA calculations.
C Pre-calculate IDEL2 indexes before starting a parallel calculation.
C In other words, this subroutine calculates how many distributions are associated with the calculation.
C
C This routine assembles a matrix that counts the number of AOs
C associated with an ILLL distribution index. This array is used by getallocsize
C and partitionAOindices to pre-sort the integrals that need to be done
C
C The first row in outlist is the number of distributions
C The second row in outlist is the associated ILLL index
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "eridst.h"

      integer,                       intent(in) :: Nindex
      integer, dimension(*),         intent(in) :: indxbt
      integer, dimension(2, Nindex), intent(out):: outlist
      integer :: i

      do i=1, Nindex
         call getdst(i, 0, 0)
         call pickao(0)
         call eridsi(indxbt, 0)
         outlist(1, i) = ndistr
         outlist(2, i) = i
      enddo

      return

      end subroutine


C     /* deck getallocsize */
      SUBROUTINE getallocsize(ntot, originalsort, maxnumjobs)
C A subroutine associated with the atomic integral parallel RPA/SOPPA calculations.
C
C This subroutine is used to get the amount of storage that needs
C to be allocated for a parallel SO_ERES run.
C
C The subroutine calculates which process will be assigned the
C most single jobs (not the largest total amount of work) for a
C parallel RPA/SOPPA calculation. Its output is used to allocate
C the right amount of storage by the master when it starts pre-sorting
C the integrals for the parallel calculation of the E matrix.
C

#include "implicit.h"
#include "mpif.h"
      integer,                     intent(in)  :: ntot
      integer, dimension(2, ntot), intent(in)  :: originalsort
      integer,                     intent(out) :: maxnumjobs

      integer, dimension(:,:), allocatable  :: copysort
      integer, dimension(:,:), allocatable  :: sumofwork
      integer, dimension(2) :: temploc, tempwork, tempout
      integer :: allocstatus, deallocstatus, numprocs
      INTEGER(MPI_INTEGER_KIND) :: ierr_mpi, numprocs_mpi



      call mpi_comm_size(mpi_comm_world, numprocs_mpi, ierr_mpi)
      numprocs = numprocs_mpi
      allocate( copysort(2, ntot), sumofwork(2, numprocs)
     &         ,stat=allocstatus)
      if(.not.(allocstatus.eq.0) ) then
         call quit('Allocation error in GETALLOCSIZE')
      endif

      !call izero(sumofwork, (2*numprocs) )
      sumofwork = 0
      copysort = originalsort

      DO i=1, ntot
         ! Find location of largest chunk of work and the work itself
         temploc = maxloc(copysort, DIM=2, mask=copysort.gt.0)
         addwork = copysort(1, temploc(1) )
         copysort( 1,temploc(1) ) = 0

         ! Find laziest slave and simulate the workload on the slave
         tempwork = minloc(sumofwork, DIM=2)
         sumofwork(1, tempwork(1)) = sumofwork(1, tempwork(1)) + addwork! adding total number of distributions
         sumofwork(2, tempwork(1)) = sumofwork(2, tempwork(1)) + 1  !adding total number of assigned indexes
      ENDDO

      tempout = maxloc(sumofwork, dim=2)
      maxnumjobs = sumofwork( 2,tempout(2) )

      deallocate(copysort, sumofwork, stat=deallocstatus)
      if(.not.(deallocstatus.eq.0) ) then
         call quit('Deallocation error in GETALLOCSIZE')
      endif

      return

      END SUBROUTINE


C     /* deck partitionAOindices */
      SUBROUTINE  partitionAOindices(ntot, rows, cols, presortedarray,
     &                               sorted, sumofwork)
C A subroutine associated with the atomic integral parallel RPA/SOPPA calculations.
C
C The output from this routine is the array 'sorted'.
C The sorted matrix contains AO integral indexes for two-electron integrals.
C Every column contains a list of ILLL indexes that are to be assigned as work
C to a single process. The total amount of work per column is estimated based
C on the number of distributions that are related to the ILLL indexes in the column.
C
C This pre-sorting of ILLL indexes for a parallel calculation approximates an
C even distribution of total work for all processes when performing two-electron
C integrals in parallel for AOSOPPA and AORPA.
C
      !use mpi
      implicit none
#include "mpif.h"

      integer,                        intent(in)    :: ntot, rows, cols
      integer, dimension(2, ntot),    intent(inout) :: presortedarray
      integer, dimension(rows, cols), intent(out)   :: sorted

      integer, dimension(2, cols), intent(out) :: sumofwork
      integer, dimension(2) :: tempavail, templazy
      integer :: numprocs, availloc, targetrow
      integer :: allocstatus, deallocstatus
      integer :: lazyloc, aoindex, i, numdists

      numprocs = cols

      !call izero(sorted, (rows*cols) )
      !call izero(sumofwork, (2*numprocs) )
      sorted = 0
      sumofwork = 0

      DO i=1, ntot
C        FIND LARGEST CHUNK OF AVAILABLE WORK AND ITS INDEX
         tempavail = maxloc(presortedarray, DIM=2,
     &                      mask=presortedarray.gt.0)
         availloc = tempavail(1)
         numdists = presortedarray(1, availloc) ! amount of available work
         aoindex  = presortedarray(2, availloc) ! the index to be passed to process

         presortedarray(1, availloc) = 0
         presortedarray(2, availloc) = 0

C        FIND THE LAZIEST PROCESS, GIVE IT WORK & INCREMENT THE ROW COUNTER
         templazy = minloc(sumofwork, DIM=2)
         lazyloc = templazy(1) ! This is equal to: MYID+1
         sumofwork( 1,lazyloc ) = sumofwork( 1,lazyloc ) + numdists
         sumofwork( 2,lazyloc ) = sumofwork( 2,lazyloc ) + 1
         targetrow = sumofwork( 2,lazyloc )

C        ADD AO-INDEX TO FIRST AVAILABLE ROW IN THE LAZY SLAVE'S COLUMN
         sorted( targetrow, lazyloc ) = aoindex
      ENDDO

      return

      END SUBROUTINE


#endif
!VAR_MPI at the beginning of dynloadbal_parsoppa
