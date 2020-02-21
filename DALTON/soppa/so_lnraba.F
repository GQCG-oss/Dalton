C
C  /* Deck so_lnraba */
      SUBROUTINE SO_LNRABA(POLDD,POLDQ,POLDL,POLDA,WORK,LWORK,PASS)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan P. A. Sauer, Februrary 1999
C     Stephan P. A. Sauer: 3.12.2003
C     Rasmus Faber, October 2015: Added MPI
C
C     PURPOSE: Main driver routine for the calculation of frequency
C              dependent linear response properties with the atomic
C              integral direct SOPPA program.
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
     &                   so_needs_densai,
     &                   sop_mp2ai_done
C
#ifdef VAR_MPI
      use so_parutils, only: soppa_initialize_slaves,
     &                       soppa_release_slaves
#endif
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "cbilnr.h"
#include "cbiexc.h"
#include "absorp.h"
#include "ccsdsym.h"
c#include "infdim.h"
#include "inforb.h"
C for irat
#include "iratdef.h"
#include "maxaqn.h"
C for mxshel, mxprim
#include "maxorb.h"
C for mxcont
#include "aovec.h"
c#include "mxcent.h"
c#include "nuclei.h"
c#include "symmet.h"
c#include "wrkrsp.h"
C NEED DOCCSD parameter from gnring.h
#include "gnrinf.h"
c     PARAMETER (HALF = 0.5D0,ESUDIP = 64604.885D0,ESUECD = 471.44360D0)
C
      DIMENSION POLDD(2,3,3,NFRVAL,*), POLDQ(2,3,3,3,NFRVAL,*)
      DIMENSION POLDL(2,3,3,NFRVAL,*), POLDA(2,3,3,NFRVAL,*)
      DIMENSION WORK(LWORK)
C
      CHARACTER*5 MODEL
      CHARACTER*8 LABEL1
      CHARACTER*11 FULLNAME
C
      integer :: num_active
      logical :: active_models(sop_num_models)
      logical :: doubles, need_t2, rpa_only
      character(len=4) :: old_denslab, denslab
      logical :: new_densai, imagprop
#ifdef VAR_MPI
!
! This variable ensures that common blocks are only sent to the slaves
! once.
      LOGICAL update_common_blocks

      update_common_blocks = .true.
#endif
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_LNRABA')
C
C
C     Get the active models and their number
      CALL SO_GET_ACTIVE_MODELS (ACTIVE_MODELS)
      NUM_ACTIVE = COUNT ( ACTIVE_MODELS )
C
C     No RPA(D) properties
      ACTIVE_MODELS(SOP_MODEL_RPAD) = .FALSE.
C
      RPA_ONLY = ACTIVE_MODELS(SOP_MODEL_RPA) .AND. 
     &           (NUM_ACTIVE.EQ.1) 
C
C     Set the correct convergence threshold -- THCLNR from cbilnr.h
      sop_conv_thresh = THCLNR
C
C     Triplet flag not initialized...
      TRIPLET = .FALSE.
C
C     We have read no amplitudes
      OLD_DENSLAB = 'NONE'
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
      IPRSOP = IPRLNR
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
C        Force a crash if this is used
         KT2AM   = HUGE(LWORK)
         KDENSIJ = HUGE(LWORK)
         KDENSAB = HUGE(LWORK)
         KDENSAI = HUGE(LWORK)
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
C---------------------------------------------------------------------
C     Initialize arrays of transition moments and excitation energies.
C---------------------------------------------------------------------
C
C      CALL DZERO(SNDPRP,2)
      CALL DZERO(POLDD,2*9*NFRVAL*NUM_ACTIVE)
      CALL DZERO(POLDQ,2*27*NFRVAL*NUM_ACTIVE)
      CALL DZERO(POLDL,2*9*NFRVAL*NUM_ACTIVE)
      CALL DZERO(POLDA,2*9*NFRVAL*NUM_ACTIVE)
C
#ifdef VAR_MPI
C
C---------------------------------------------------------
C              Ready the slaves for parallel calculations.
C---------------------------------------------------------
C
      call soppa_initialize_slaves ( update_common_blocks, rpa_only)
      update_common_blocks = .false.
#endif
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
C--------------------------------------------------
C     Loop over symmetry of the property operators.
C--------------------------------------------------
C
      DO 200 ISYM = 1, NSYM
C
         DO 100 IOPER = 1, NOPER(ISYM)
C
            LABEL1 = LABOP(IOPER,ISYM)
C Do only for DIPLEN?
            IF (LABEL1(2:7).NE.'DIPLEN') CYCLE

CRF is there a common pattern here?
            IDIP = ICHAR(LABEL1(1:1)) - ICHAR('X') + 1

            ISYMTR = ISYM
C
C========================================
C           Loop over the posible methods
C========================================
C
            IOUT = 1
C
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
               DENSLAB = SOP_DENS_LABEL(IMODEL)
C
               IF (IPRSOP .GE. 1) THEN
                  WRITE(LUPRI,9000)
                  WRITE(LUPRI,'(1X,A,2(/A,I8))') TRIM(FULLNAME)//':',
     &            ' Perturbation symmetry     :',ISYM,
     &            ' p-h + h-p variables       :',2*NT1AM(ISYM)
                  IF (DOUBLES) THEN
                     NVARPT = 2*(NT1AM(ISYM)+NT2AM(ISYM))
                     WRITE(LUPRI,'(A,I8,/A,I8)')
     &               ' 2p-2h + 2h-2p variables   :',2*NT2AM(ISYM),
Cend-PFP
     &               ' Total number of variables :',NVARPT
                  ENDIF
                  WRITE(LUPRI,9001)
               END IF
C
C-----------------------------------------------------
C              Get T2 amplitudes and density matrices.
C-----------------------------------------------------
C
               IF ((NEED_T2).AND.
     &             (DENSLAB.NE.OLD_DENSLAB)) THEN
                  CALL SO_GETT2(DENSLAB,
     &                          WORK(KFOCKD),LFOCKD,WORK(KT2AM),LT2AM,
     &                          WORK(KDENSAI),
     &                          LDENSAI,WORK(KDENSIJ),LDENSIJ,
     &                          WORK(KDENSAB),LDENSAB,
     &                          WORK(KEND1),LWORK1)
                  OLD_DENSLAB = DENSLAB
               ENDIF
C
C-----------------------------------------
C              Initialize DENSAI if needed
C-----------------------------------------
C
               IF ( SO_NEEDS_DENSAI(MODEL) .AND.
     &              (.NOT.SOP_MP2AI_DONE) ) THEN
                  CALL DZERO(WORK(KDENSAI),LDENSAI)
                  NEW_DENSAI = .TRUE.
               ELSE
                  NEW_DENSAI = .FALSE.
               ENDIF
C
C------------------------------------------------------------------
C              Determine SOPPA excitation energies and excitation
C              vectors. The excitation vectors are written to file.
C------------------------------------------------------------------
C
               CALL SO_RSPLEQ(MODEL,LABEL1,IMAGPROP,ISYMTR,FRVAL,NFRVAL,
     &                        WORK(KDENSIJ),LDENSIJ,
     &                        WORK(KDENSAB),LDENSAB,
     &                        WORK(KDENSAI),LDENSAI,WORK(KT2AM),
     &                        LT2AM,WORK(KFOCKD),LFOCKD,
     &                        WORK(KEND1),LWORK1)
C
C              Save AI density matrix, if it has been recalculated.
               IF (NEW_DENSAI) THEN
                  CALL SO_DUMP_DENSAI(DENSLAB,WORK(KDENSAI),LDENSAI)
               END IF
C
C----------------------------------------------
C              Calculate second order property.
C----------------------------------------------
C
               CALL SO_POLAR(MODEL,ISYMTR,IDIP,LABEL1,NLBTOT,
     &                       WORK(KT2AM),LT2AM,
     &                       WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                       LDENSAB,WORK(KDENSAI),LDENSAI,
     &                       POLDD(1,1,1,1,IOUT),POLDQ(1,1,1,1,1,IOUT),
     &                       POLDL(1,1,1,1,IOUT),POLDA(1,1,1,1,IOUT),
     &                       WORK(KEND1),LWORK1)
C              INCREMENT OUTPUT POINTER
               IOUT = IOUT + 1
C
            END DO

C
  100    CONTINUE
C
  200 CONTINUE
C
#ifdef VAR_MPI
C-------------------------------------------------------
C              Release slaves to the global node-driver.
C-------------------------------------------------------
      call soppa_release_slaves()
#endif
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
      CALL QEXIT('SO_LNRABA')
C
      RETURN
C
 9000 FORMAT(/' -----------------------------------')
 9001 FORMAT(' -----------------------------------')
      END
