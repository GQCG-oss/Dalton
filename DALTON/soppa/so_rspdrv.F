C
C  /* Deck so_rsplr */
      SUBROUTINE SO_RSPDRV(WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Rasmus Faber, December 2015: Added MPI
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
     &                   sop_conv_thresh
C
#ifdef VAR_MPI
      use so_parutils, only: soppa_initialize_slaves,
     &                       soppa_release_slaves
#endif
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
C#include "cbiexc.h"
C From infrsp we need
C FREQ <- Array of requecies
C NFREQ <- Number of frequencies
#include "infrsp.h"
C inflr.h below requires MAXLBL from rspprp.h
#include "rspprp.h"
C We need
C LBLLR <- Label of operators
C NGPLR <- Number of operators of each symmetry
C THCLR <- Convergence threshold
#include "inflr.h"
#include "ccsdsym.h"
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
      DIMENSION WORK(LWORK)
C
      CHARACTER*5 MODEL
      CHARACTER*8 LABEL1
      CHARACTER*11 FULLNAME
      CHARACTER(len=7), parameter :: dash5 = ' ----- '
      CHARACTER(len=65), parameter :: dashl = ' -------------------'//
     &               '---------------------------------------------'
C
      integer :: num_active
      logical :: active_models(sop_num_models)
      logical :: doubles, need_t2, rpa_only
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
      CALL QENTER('SO_RSPDRV')
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
C     Set the correct convergence threshold -- THCLR from inflr.h
      sop_conv_thresh = THCLR
C
C     Set triplet flag if given -- TRPLET from infrsp.h
      IF (TRPLET) THEN 
         TRIPLET = .TRUE.
      ELSE
         TRIPLET = .FALSE.
      ENDIF
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
C-------------------------------------------
C     Initializes memory statistics routine.
C-------------------------------------------
C
      CALL SO_MEMMAX ('START',LWORK)
C
C---------------------------------------------------------------
C     If only RPA initialize some pointers and sort AO-integrals
C     into distributions.
C---------------------------------------------------------------
C     On the flip side, this cause a crash if CC has been run
      IF ( RPA_ONLY ) THEN
         CALL RP_INIT(WORK,LWORK)
#ifdef VAR_MPI
C     In pure aorpa calculations, eridi has not been initialized
         herdir = so_get_herdir()
         get_mxcall = so_get_direct()
#endif
      ENDIF

#ifdef VAR_MPI
C
C In AORPA we cannot be sure that the number of integral distributions
C is known in advance. This is needed in order to allocate memory for
C load-distribution later on
C
      if (get_mxcall) then
         if (herdir) then
            call herdi1(work(1), lwork, iprint)
         else
            kccfbt  = 1
            kindxbt = kccfbt + mxprim*mxcont
            kend3   = kindxbt + 8*mxshel*mxcont/irat
            lwork3 = lwork - kend3
            call eridi1(kodcl1, kodcl2, kodbc1, kodbc2,
     &                  krdbc1, krdbc2, kodpp1, kodpp2,
     &                  krdpp1, krdpp2, kfree,  lfree,
     &                  kend3, work(kccfbt), work(kindxbt),
     &                  work(kend3), lwork3, iprint)
            endif
         get_mxcall = .false.
      endif
#endif
C
C---------------------------------------------------
C     Initialize a few pointerarrays and dimensions.
C---------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_INIT
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
C        Force a crash of this is used
         KT2AM   = LWORK + 1
         KDENSIJ = LWORK + 1
         KDENSAB = LWORK + 1
         KDENSAI = LWORK + 1
      END IF
C
      CALL SO_MEMMAX ('SO_LNRABA.1',LWORK1)
      IF (LWORK1 .LT .0) CALL STOPIT('SO_LNRABA.1',' ',KEND1,LWORK)
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
C==========================================
C     Determine linear response properties.
C==========================================
C
C----------------------------------------------------------
C     Adjust number of trialvectors for each excitation and
C     write information to output.
C----------------------------------------------------------
C
      IF ( NSAVMX .LT. 2 ) THEN
C
         NSAVMX = 2
C
         WRITE(LUPRI,'(1X,A,/,A,I2,A)')
     &   'NOTICE: Maximum number of trial vectors for each'//
     &   ' property',' is raised to',NSAVMX, ' as this is'//
     &   ' minimal space allowed.'
C
      END IF
C
      WRITE(LUPRI,'(/,1X,A,I2,/)')
     &'Maximum number of trial vectors for each'//
     &' property is ',NSAVMX
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
         IF(NEED_T2) THEN
             CALL SO_GETT2(sop_dens_label(imodel),
     &                     WORK(KT2AM),LT2AM,WORK(KDENSAI),
     &                     LDENSAI,WORK(KDENSIJ),LDENSIJ,
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
C--------------------------------------------------
C     Loop over symmetry of the property operators.
C--------------------------------------------------
C
         DO ISYM = 1, NSYM
C
            IF (NGPLR(ISYM).LE.0) CYCLE ! Nothing to do for this sym
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
     &           NGPLR(ISYM)

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

C            DO IEXCI = 1, NPPCNV(ISYM)
C               WRITE(LUPRI,'(2A,/A)') 'Warning: Excitation energies',
C     &            'not implemented in AO-SOPPA **RESPONSE driver.',
C     &            'Use **PROPERTIES instead'
C            END DO

            DO IOPER = 1, NGPLR(ISYM)
C
               LABEL1 = LBLLR(ISYM,IOPER)
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
               CALL FLUSH(LUPRI)
            END DO ! Loop operators
C
C----------------------------------------------
C              Calculate second order property.
C----------------------------------------------
C
            IF (NGPLR(ISYM).GT.0) THEN
               CALL SO_LRNSL(MODEL,ISYM,
     &                       WORK(KT2AM),LT2AM,
     &                       WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                       LDENSAB,WORK(KDENSAI),LDENSAI,
     &                       WORK(KEND1),LWORK1)
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
      CALL SO_MEMMAX ('SO_LNRABA.3     ',LWORK3)
      IF (LWORK3 .LT.0) CALL STOPIT('SO_LNRABA.3',' ',KEND3,LWORK)
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
      CALL QEXIT('SO_RSPDRV')
C
      RETURN
C
 9000 FORMAT(/' -----------------------------------')
 9001 FORMAT(' -----------------------------------')
      END
