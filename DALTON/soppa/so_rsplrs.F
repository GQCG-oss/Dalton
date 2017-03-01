C
C  /* Deck so_rsplrs */
      SUBROUTINE SO_RSPLRS(LABELS,NUM_LABELS,LABEL_OFFSETS,
     &                     TRIP_INPUT,THRESH,FREQ,NFREQ,IMODE,
     &                     WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Rasmus Faber, Feb. 2017
C
C     PURPOSE: A RSP compatible driver for calculating linear response properties with 
C              the atomic integral direct SOPPA program.
C
      use so_info, only: sop_num_models,
     &                   sop_models, sop_mod_fullname,
     &                   so_num_active_models,
     &                   so_get_active_models,
     &                   so_has_doubles,
     &                   so_singles_second,
     &                   sop_dens_label, sop_model_rpad,
     &                   sop_model_rpa,
     &                   sop_conv_thresh,
     &                   sop_use_seller,
     &                   sop_dp,
     &                   so_needs_densai,
     &                   FN_PROP
C
#ifdef VAR_MPI
      use so_parutils, only: soppa_initialize_slaves,
     &                       soppa_release_slaves
#endif
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C  nsym
#include "inforb.h"
C for irat
#include "iratdef.h"
#include "maxaqn.h"
C for mxshel, mxprim
#include "maxorb.h"
C for mxcont
#include "aovec.h"
C REP <- Symbols of the ireducible representations.
#include "pgroup.h"
C
C  Input arguments      
      REAL(SOP_DP), INTENT(INOUT) :: WORK(LWORK) ! WORK array
      REAL(SOP_DP), INTENT(IN) :: FREQ(NFREQ) ! Array of frequencies
      REAL(SOP_DP), INTENT(IN) :: THRESH ! Convergence treshhold
      INTEGER, INTENT(IN) :: NFREQ, LWORK
      INTEGER, INTENT(IN) :: LABEL_OFFSETS(8), NUM_LABELS(8)
      CHARACTER(LEN=8),INTENT(IN) :: LABELS(*)
      LOGICAL, INTENT(IN) :: TRIP_INPUT ! Triplet flag
      INTEGER, INTENT(IN) :: IMODE ! How we are called
                                   ! IMODE = 0 : from response
                                   ! IMODE = 1, 2, 3: real, imag,
                                   ! triplet from abacus
C
      CHARACTER*5 MODEL
      CHARACTER*8 LABEL1
      CHARACTER*11 FULLNAME
      CHARACTER(len=7), parameter :: dash5 = ' ----- '
      CHARACTER(len=65), parameter :: dashl = ' -------------------'//
     &               '---------------------------------------------'
C
      character(len=9), parameter :: myname = 'SO_RSPLRS'
      integer :: num_active
      logical :: active_models(sop_num_models)
      logical :: doubles, need_t2, rpa_only
      integer :: lu_dump ! Filenumber of file property values to
#ifdef VAR_MPI
!
! This variable ensures that common blocks are only sent to the slaves
! once.
      LOGICAL update_common_blocks, get_mxcall
      LOGICAL so_get_herdir, herdir, so_get_direct

      update_common_blocks = .true.
      get_mxcall = .false.
#endif
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER(myname)
C
C     Get the active models and their number
      CALL SO_GET_ACTIVE_MODELS (ACTIVE_MODELS)
C
C     No RPA(D) properties
      ACTIVE_MODELS(SOP_MODEL_RPAD) = .FALSE.
      NUM_ACTIVE = COUNT ( ACTIVE_MODELS )
C
      RPA_ONLY = ACTIVE_MODELS(SOP_MODEL_RPA) .AND. 
     &           (NUM_ACTIVE.EQ.1) 
C
C     Set the correct convergence threshold 
      sop_conv_thresh = THRESH
C
C     Set triplet flag if given
      TRIPLET = TRIP_INPUT
C
C     Open file we will output results on
      IF (IMODE .GT. 0) THEN
         IF (IMODE .GT. 3) THEN
            WRITE(LUPRI,'(a,i3,a)') 'ERROR!! imode', imode,' is invalid'
            RETURN
         END IF
         LU_DUMP = -1
         CALL GPOPEN(LU_DUMP,FN_PROP(IMODE),'NEW',' ',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
      END IF
C
C-------------------------------------------
C     Start timing the AO-SOPPA calculation.
C-------------------------------------------
C
      CALL TIMER('START ',TIMEIN,TIMOUT)
C
      DTIME   = SECOND()
      TIMTOT  = DTIME
      CALL GETTIM (DUMMY,WTIME)
      TIMWTO  = WTIME
C
C---------------------------------------------------
C     Initialize a few pointerarrays and dimensions.
C---------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_INIT(WORK,LWORK)
      DTIME     = SECOND()  - DTIME
      SOTIME(2) = SOTIME(2) + DTIME
C
C----------------------------------
C     Set print level for AO-SOPPA.
C----------------------------------
C
      IPRSOP = 2
C
      IF (IPRLNR .GE. 0)
     &     CALL TITLER('Solving Linear Response Equations','#',118)
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LFOCKD  = NORBT
C
      KFOCKD  = 1
      KEND1   = KFOCKD  + LFOCKD
      LWORK1  = LWORK   - KEND1
C
C
      IF (.NOT.RPA_ONLY) THEN
C
         LT2AM   = NT2AMX
C        Triplet case currently requires
C        double the memory for MP2 amplitudes         
         IF(TRIPLET) LT2AM = 2*LT2AM
         LDENSIJ = NIJDEN(1)
         LDENSAB = NABDEN(1)
         LDENSAI = NAIDEN(1)
C
         KT2AM   = KFOCKD  + LFOCKD
         KDENSIJ = KT2AM   + LT2AM
         KDENSAB = KDENSIJ + LDENSIJ
         KDENSAI = KDENSAB + LDENSAB
         KEND1   = KDENSAI + LDENSAI
         LWORK1  = LWORK   - KEND1
      ELSE
         LT2AM   = 0
         LDENSIJ = 0
         LDENSAB = 0
         LDENSAI = 0
C        Force a crash if this is used
         KT2AM   = LWORK + 1
         KDENSIJ = LWORK + 1
         KDENSAB = LWORK + 1
         KDENSAI = LWORK + 1
      END IF
C
      CALL SO_MEMMAX (myname//'.1',LWORK1)
      IF (LWORK1 .LT .0) CALL STOPIT(myname//'.1',' ',KEND1,LWORK)
C
C------------------------------------------------
C     Get MO-energies (the fock matrix diagonal).
C------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_MOENERGY(WORK(KFOCKD),WORK(KEND1),LWORK1)
      DTIME     = SECOND()  - DTIME
      SOTIME(3) = SOTIME(3) + DTIME
C
C------------------------------------------------------
C     Construct property-integrals and write to LUPROP.
C------------------------------------------------------
C
      CALL SO_PRPINT('LINEAR',NLBTOT,WORK(KEND1),LWORK1)
C
C     Find some otherway of setting this?
C
      NSAVMX = 3
      WRITE(LUPRI,'(/,1X,A,I2,/)')
     &'Maximum number of trial vectors for each'//
     &' property is ',NSAVMX
C
C==========================================
C     Determine linear response properties.
C==========================================
C
C========================================
C     Loop over the posible methods
C========================================
C
      IOUT = 1
      DO IMODEL = 1, SOP_NUM_MODELS
C
C  Skip methods not treated in this calcultion
         IF (.NOT.ACTIVE_MODELS(IMODEL) ) CYCLE
C
C  Look up info on this model
C  (the model number -> model conversion is defined in so_info.F90)
C
         MODEL = SOP_MODELS(IMODEL)
         FULLNAME = SOP_MOD_FULLNAME(IMODEL)
         DOUBLES = SO_HAS_DOUBLES(IMODEL)
         NEED_T2 = DOUBLES.OR.SO_SINGLES_SECOND(MODEL)

         CALL TITLER( TRIM(FULLNAME) // ' response calculation', 
     &                '*',103)
C
C-----------------------------------------------------
C              Get T2 amplitudes and density matrices.
C-----------------------------------------------------
C
         IF (NEED_T2) THEN
             CALL SO_GETT2(sop_dens_label(imodel),
     &                     WORK(KFOCKD),LFOCKD,
     &                     WORK(KT2AM),LT2AM,
     &                     WORK(KDENSAI),LDENSAI,
     &                     WORK(KDENSIJ),LDENSIJ,
     &                     WORK(KDENSAB),LDENSAB,
     &                     WORK(KEND1),LWORK1)
         END IF
#ifdef VAR_MPI
C
C---------------------------------------------------------
C        Ready the slaves for parallel calculations.
C---------------------------------------------------------
C
         call soppa_initialize_slaves ( update_common_blocks,
     &                                   WORK(KT2AM), LT2AM, model)
         update_common_blocks = .false.
#endif
C
C-----------------------------------------
C              Initialize DENSAI if needed
C-----------------------------------------
C        
         IF ( SO_NEEDS_DENSAI(MODEL))
     &        CALL DZERO(WORK(KDENSAI),LDENSAI)
C
C--------------------------------------------------
C     Loop over symmetry of the property operators.
C--------------------------------------------------
C
         DO ISYM = 1, NSYM
C
            IF (NUM_LABELS(ISYM).LE.0) CYCLE ! Nothing to do for this sym
            ILABOFF = LABEL_OFFSETS(ISYM)
C            
            WRITE(LUPRI,'(A)') dashl
            WRITE(LUPRI,'(2A,I5,3A)') dash5,
     &       'Symmetry of excitation/property operator(s)',
     &           ISYM,'  ( ',REP(ISYM-1),')'
            WRITE(LUPRI,'(A)') dashl

            WRITE(LUPRI,'(2(/A,I5))') 
     &           ' Number of excitations of this symmetry        ',
C     &           NPPCNV(ISYM),
     &           0, ! Fix this!!!
     &           ' Number of response properties of this symmetry',
     &           NUM_LABELS(ISYM)

            WRITE(LUPRI,9000)
            WRITE(LUPRI,'(1X,A,2(/A,I8))') TRIM(FULLNAME)//':',
     &         ' Perturbation symmetry     :',ISYM,
     &         ' p-h + h-p variables       :',2*NT1AM(ISYM)
            IF (DOUBLES) THEN
               NVARPT = 2*(NT1AM(ISYM)+NT2AM(ISYM))
               WRITE(LUPRI,'(A,I8,/A,I8)')
     &         ' 2p-2h + 2h-2p variables   :',2*NT2AM(ISYM),
     &         ' Total number of variables :',NVARPT
            ENDIF
            WRITE(LUPRI,9001)

            SOP_USE_SELLER = (NUM_LABELS(ISYM).GT.1)

            DO IOPER = 1, NUM_LABELS(ISYM)
C
               LABEL1 = LABELS(ILABOFF+IOPER)
               ISYMTR = ISYM
               IOUT = 1
C
C------------------------------------------------------------------
C              Determine SOPPA linear response vectors and density
C              matrixes. Response density matrix written on file.
C------------------------------------------------------------------
C
               CALL SO_RSPLEQ(MODEL,LABEL1,ISYMTR,FREQ,NFREQ,
     &                        WORK(KDENSIJ),LDENSIJ,
     &                        WORK(KDENSAB),LDENSAB,
     &                        WORK(KDENSAI),LDENSAI,WORK(KT2AM),
     &                        LT2AM,WORK(KFOCKD),LFOCKD,
     &                        WORK(KEND1),LWORK1)
C
C               CALL FLUSH(LUPRI)
            END DO ! Loop operators
C
C----------------------------------------------
C              Calculate second order property.
C----------------------------------------------
C
            IF (NUM_LABELS(ISYM).GT.0) THEN
               KPROP = KEND1
               LPROP = NFREQ*(NUM_LABELS(ISYM)**2)
               KEND2 = KPROP + LPROP
               LWORK2 = LWORK-KEND2
               CALL SO_MEMMAX (myname//'.2     ',LWORK2)
               IF (LWORK2 .LT.0) CALL STOPIT(myname//'.2',' ',
     &               KEND2,LWORK)

               CALL SO_LRNSL(MODEL,ISYM,WORK(KPROP),
     &                       FREQ,NFREQ,LABELS(ILABOFF+1),
     &                       NUM_LABELS(ISYM),
     &                       WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                       LDENSAB,WORK(KDENSAI),LDENSAI,
     &                       WORK(KEND2),LWORK2)
C
               IF (IMODE.GT.0) THEN
                  WRITE(LU_DUMP) ISYM
                  WRITE(LU_DUMP) WORK(KPROP:KPROP+LPROP-1)
               ELSE
                  CALL SO_RSPOUT(WORK(KPROP),ISYM,FREQ,NFREQ,
     &                        LABELS(ILABOFF+1),NUM_LABELS(ISYM))
               END IF
C
            END IF
C
         END DO ! LOOP over symmetries

#ifdef VAR_MPI
C-------------------------------------------------------
C              Release slaves to the global node-driver.
C-------------------------------------------------------
         call soppa_release_slaves()
#endif
C
C        INCREMENT OUTPUT POINTER
C
         IOUT = IOUT + 1
C
      END DO ! Loop over models
C
      IF (IMODE .GT. 0) THEN
         CALL GPCLOSE(LU_DUMP,'KEEP')
      END IF
C
C---------------------------------
C     3. allocation of work space.
C---------------------------------
C
      LPARRA = LSOTIM
C
      KPARRA = KEND1
      KEND3  = KPARRA + LPARRA
      LWORK3 = LWORK  - KEND3
C
      CALL SO_MEMMAX (myname//'.3     ',LWORK3)
      IF (LWORK3 .LT.0) CALL STOPIT(myname//'.3',' ',KEND3,LWORK)
C
C---------------------------------------------------
C     Print memory statistics for SOPPA subroutines.
C---------------------------------------------------
C
      CALL SO_MEMMAX('STATISTICS      ',0)
C
C-----------------------------------------
C     Print timings for SOPPA subroutines.
C-----------------------------------------
C
      CALL SO_TIME(TIMTOT,TIMWTO,WORK(KPARRA),LPARRA)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT(myname)
C
      RETURN
C
 9000 FORMAT(/' -----------------------------------')
 9001 FORMAT(' -----------------------------------')
      END

      SUBROUTINE SO_FCSDDRV(WORK,LWORK)
C
C     Wrapper for SO_RSPLRS used for calculating the
C     FC and SD contribution to spin-spin coupling constants
C
      use so_info, only: sop_dp
C  Some of the following blocks still don't declare variables
#include "implicit.h"
C  NSYM
#include "inforb.h"
C  gdvec.h and abares.h requires MXCOOR, MXCENT
#include "mxcent.h"
C  TRPLAB <- Vector containing the required labels.
#include "abares.h"
C  NTRVEC <- The number of labels in each symmetry
#include "gdvec.h"
C  THRTRP <- Convergence criteria
#include "cbitrp.h"
      real(sop_dp), intent(inout) :: work(lwork)
      integer, intent(in) :: lwork

      real(sop_dp), dimension(1), parameter :: freq = (/ 0.0D0 /)
      integer :: label_offsets(8)
      integer :: n, ncount
C
C     Calculate label offsets
      ncount = 0
      do n = 1, nsym
         label_offsets(n) = ncount
         ncount = ncount + ntrvec(n)
      end do
C
C     We now have all info for SO_RSPLRS
      CALL SO_RSPLRS(TRPLAB(1), NTRVEC, LABEL_OFFSETS,
     &               .TRUE., THRTRP, FREQ, 1, 3,
     &               WORK, LWORK)

      RETURN
      END SUBROUTINE

      SUBROUTINE SO_READPROP(REL,NOP,IMODE,ISYM)

      use so_info, only: fn_prop, sop_dp
      implicit none
      real(sop_dp), intent(out) :: rel(nop,nop)
      integer, intent(in) :: nop, imode, isym
      integer :: lu, idummy, this_sym
      logical :: ldummy

      LU = -1
      CALL GPOPEN(LU, fn_prop(imode),'OLD',' ',
     &            'UNFORMATTED',IDUMMY,LDUMMY)
      do
         read (lu, end=900) this_sym
         if (this_sym .eq. isym) exit
         read (lu, end=900)
      end do
      read (lu) rel
C     This factor is needed for triplet properties...
      if (imode.eq.3) rel = -0.25D0*rel

      CALL GPCLOSE(lu,'KEEP')
      return

 900  CALL QUIT('ERROR reading file '//fn_prop(imode))

      END SUBROUTINE

      SUBROUTINE SO_PSODRV(WORK,LWORK)
C
C     Wrapper for SO_RSPLRS used for calculating the
C     PSO contribution to spin-spin coupling constants
C
      use so_info, only: sop_dp
C  Some of the following blocks still don't declare variables
#include "implicit.h"
C  NSYM
#include "inforb.h"
C  gdvec.h and abares.h requires MXCOOR, MXCENT
#include "mxcent.h"
C  IMGLAB <- Vector containing the required labels.
#include "abares.h"
C  NGDVEC <- The number of labels in each symmetry
#include "gdvec.h"
C  THRCLC <- Convergence criteria
#include "cbilrs.h"
      real(sop_dp), intent(inout) :: work(lwork)
      integer, intent(in) :: lwork

      real(sop_dp), dimension(1), parameter :: freq = (/ 0.0D0 /)
      integer :: label_offsets(8)
      integer :: n, ncount
C
C     Calculate label offsets
      ncount = 0
      do n = 1, nsym
         label_offsets(n) = ncount
         ncount = ncount + ngdvec(n,2)
      end do
C
C     We now have all info for SO_RSPLRS
      CALL SO_RSPLRS(IMGLAB, NGDVEC(1,2), LABEL_OFFSETS,
     &               .FALSE., THRCLC, FREQ, 1, 2,
     &               WORK, LWORK)

      RETURN
      END SUBROUTINE
