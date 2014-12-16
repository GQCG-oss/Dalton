!> mm_driver:
!> General module for driving multipole method calculations. 
!> 2002 (c) Mark Watson(*), Pawel Salek
! 
! Notes: traditional dalton's WORK memory is used for some of the
! larger allocations before memory allocation in dalton is "fixed".
!
#define ALLOCATION_WITHOUT_WORK 1

MODULE mm_main_driver
   use files
   use mm_mem
!ANDREAS
   USE mm_stats_mod
!
   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_get_E_via_raw_potentials,      &
             mm_get_grad_via_raw_potentials,   &
             mm_get_J_via_raw_potentials,      &
             mm_get_int_type

   REAL(REALK), POINTER, SAVE :: Vff(:,:)
   TYPE(raw_mm_data),    SAVE :: LHS_mms, RHS_mms

CONTAINS

!-------------------------------------------------------------------------------

!> \brief Initialize the mm-driver
!> \author Mark Watson
!> \date 2002
!> \param scheme Contains data for mm-calculation
!> \param work   Static memory string (if not ALLOCATION_WITHOUT_WORK)
!> \param lwork  Remaining length of work
   SUBROUTINE mm_init_driver(scheme,WORK,LWORK)

      USE density_fitting
      USE mm_interface_mod,   ONLY: mm_get_raw_data
      USE mm_aux_qlm_builder, ONLY: mm_init_qlm_builder,   &
                                    mm_get_aux_qlm

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER,      INTENT(IN) :: LWORK
      REAL(REALK),     INTENT(INOUT) :: WORK(:) !"input value of LWORK" words

      CALL mm_nullify_local_pointers
      ! get multipole data straight from interface file
      CALL mm_get_raw_data(scheme,LHS_mms,RHS_mms) 

      ! now ensure appropriate prefactors and normalisation
      CALL mm_init_qlm_builder 
      CALL mm_get_aux_qlm(scheme,LHS_mms,RHS_mms)

      ! allocate the far field potential
      CALL mm_allocate_Vff(scheme,WORK,LWORK)

   END SUBROUTINE mm_init_driver

!-------------------------------------------------------------------------------

!> \brief Frees the mm-driver
!> \author Mark Watson
!> \date 2002
   SUBROUTINE mm_free_driver

      USE mm_aux_qlm_builder, ONLY: mm_free_qlm_builder

      IMPLICIT NONE
      CALL mm_free_qlm_builder
      IF (.NOT.ASSOCIATED(Vff)) CALL LSQUIT('Vff should be allocated by now!',-1)
#ifdef ALLOCATION_WITHOUT_WORK
      call mem_dealloc_fmm(Vff)
!      DEALLOCATE(Vff)      
#endif
      CALL mm_nullify_local_pointers

   END SUBROUTINE mm_free_driver

!-------------------------------------------------------------------------------

!> \brief Nullifies local pointer
!> \author Mark Watson
!> \date 2002
   SUBROUTINE mm_nullify_local_pointers

      IMPLICIT NONE
      NULLIFY (Vff)
      NULLIFY (LHS_mms%paras, LHS_mms%dens, LHS_mms%batch_map)
      NULLIFY (LHS_mms%qlm, LHS_mms%qlm_W, LHS_mms%qlm_T)
      NULLIFY (RHS_mms%paras, RHS_mms%dens, RHS_mms%batch_map)
      NULLIFY (RHS_mms%qlm, RHS_mms%qlm_W, RHS_mms%qlm_T)
      IF(stat_do_grad) THEN
         NULLIFY (RHS_mms%qlm_der, LHS_mms%qlm_der, LHS_mms%qlm_der_T)
         NULLIFY (LHS_mms%mom2atom)
      END IF


   END SUBROUTINE mm_nullify_local_pointers

!-------------------------------------------------------------------------------

!> \brief Allocate the far-field
!> \author Mark Watson
!> \date 2002
!> \param scheme Contains data for mm-calculation
!> \param work   Static memory string (if not ALLOCATION_WITHOUT_WORK)
!> \param lwork  Remaining length of work
   SUBROUTINE mm_allocate_Vff(scheme,WORK,LWORK)

      USE mm_memory_manager_mod, ONLY: mm_allocate

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER,   INTENT(IN)    :: LWORK
      REAL(REALK),     INTENT(INOUT) :: WORK(:) ! input value of LWORK 

      INTEGER :: i, lm_dim, mms_dim, max_Lmin

      IF (.NOT.ASSOCIATED(LHS_mms%paras)) STOP 'mms ptrs not set in mm_driver!'
      IF (ASSOCIATED(Vff)) CALL LSQUIT('Vff should NOT be allocated already!',-1)

      mms_dim = SIZE(LHS_mms%paras)
      IF (scheme%dynamic_LMAX_on) THEN
         max_Lmin = 0
         DO i = 1, mms_dim
            max_Lmin = MAX(max_Lmin, LHS_mms%paras(i)%Lmin)
         END DO
         max_Lmin = MIN((max_Lmin+scheme%LEXTRA), scheme%raw_LMAX)
         lm_dim = (1+ max_Lmin)**2
      ELSE
         ! We assume we are generating a RAW potential with limited l,m
         ! ie. we shouldn't use this array for translated, BOXED potentials
         lm_dim = (1+ scheme%raw_LMAX)**2
      END IF

#ifdef ALLOCATION_WITHOUT_WORK
!      CALL mm_allocate(MEM_RAW_VFF,Vff,lm_dim,mms_dim)
      call mem_alloc_fmm(Vff,lm_dim,mms_dim)
#else
      CALL set_v(Vff,lm_dim, mms_dim, WORK, LWORK/lm_dim)
      !LWORK = LWORK - lm_dim*mms_dim
#endif
      ! must zero as Vff is build additively in general
      Vff(:,:) = zero

   CONTAINS

      SUBROUTINE set_v(v,NROWS,NCOLS,W,NWORK)

         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NROWS, NCOLS, NWORK
         REAL(REALK), POINTER      :: v(:,:)
         REAL(REALK), TARGET       :: W(NROWS,NWORK)

         IF(NWORK.LT.NCOLS) THEN
            WRITE(LUPRI,*) 'WORK space exceded for Vff!'
            WRITE(LUPRI,*) 'Require', NROWS*NCOLS, 'words'
            WRITE(LUPRI,*) 'Avaliable words =', NWORK*NROWS
            CALL LSQUIT('Work space exceded!',lupri)
         END IF
         v => w(:,1:NCOLS)

      END SUBROUTINE set_v

   END SUBROUTINE mm_allocate_Vff

!-------------------------------------------------------------------------------

!> \brief Caclulate gradient from the raw potentials
!> \author Andreas Krapp
!> \date 2010
!> \param scheme Contains data for mm-calculation
!> \param grad   The moleular gradient mm-contribution
!> \param energy The mm-energy contribution
!> \param numnuc The number of nuclei
!> \param work   Static memory string (if not ALLOCATION_WITHOUT_WORK)
!> \param lwork  Remaining length of work
   SUBROUTINE mm_get_grad_via_raw_potentials(scheme,grad,energy,numnuc,WORK,LWORK)

      USE mm_Vff_driver,    ONLY: mm_get_raw_Vff
      USE mm_Vff_processor, ONLY: mm_get_grad_from_Vff,mm_get_grad_from_pkd_Vff
      USE density_fitting

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(INOUT) :: scheme
      REAL(REALK),        INTENT(INOUT) :: WORK(:)     ! "input value of LWORK"
      REAL(REALK),        INTENT(INOUT) :: grad(numnuc*3), energy
      INTEGER,      INTENT(IN)    :: LWORK, numnuc

      ! we always have the density on both sides 
      scheme%LHS_dens = .TRUE.
      scheme%RHS_dens = .TRUE.

      ! set up T-pair generator algorithm
      IF (scheme%LHS_mm_range == scheme%RHS_mm_range) THEN
         IF (.NOT.scheme%T_searcher(DO_FQ)%all_square)  &
            scheme%T_searcher(DO_FQ)%shape = SHAPE_TRIANGULAR
         IF (.NOT.scheme%T_searcher(DO_NN)%all_square)  &
            scheme%T_searcher(DO_NN)%shape = SHAPE_TRIANGULAR
      END IF

      ! get moments and moment derivatives 
      CALL mm_init_driver(scheme,WORK,LWORK)
      ! get potential 
      CALL mm_get_raw_Vff(scheme,LHS_mms%paras,RHS_mms,Vff)
      ! get gradient
      IF (scheme%pack_LHS) THEN
         CALL mm_get_grad_from_pkd_Vff(scheme,LHS_mms,Vff,grad,energy,numnuc)
      ELSE
         CALL mm_get_grad_from_Vff(scheme,LHS_mms,Vff,grad,energy,numnuc)
      END IF

      CALL mm_free_driver

   END SUBROUTINE mm_get_grad_via_raw_potentials

!-------------------------------------------------------------------------------

!> \brief Caclulate energy from the raw potentials
!> \author Mark Watson
!> \date 2002
!> \param scheme Contains data for mm-calculation
!> \param energy The mm-energy contribution
!> \param text   The energy-contribution type being calculated
!> \param work   Static memory string (if not ALLOCATION_WITHOUT_WORK)
!> \param lwork  Remaining length of work
   SUBROUTINE mm_get_E_via_raw_potentials(scheme,energy,text,WORK,LWORK)

      USE mm_Vff_driver,    ONLY: mm_get_raw_Vff
      USE mm_Vff_processor, ONLY: mm_get_E_from_Vff, mm_get_E_from_pkd_Vff
      USE density_fitting

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(INOUT) :: scheme
      REAL(REALK),        INTENT(OUT)   :: energy
      CHARACTER(*),       INTENT(OUT)   :: text
      INTEGER,      INTENT(IN)    :: LWORK
      REAL(REALK),        INTENT(INOUT) :: WORK(:) ! "input value of LWORK"

      ! we always have the density on both sides when getting energies
      scheme%LHS_dens = .TRUE.
      scheme%RHS_dens = .TRUE.

      ! set up T-pair generator algorithm
      IF (scheme%LHS_mm_range == scheme%RHS_mm_range) THEN
         IF (.NOT.scheme%T_searcher(DO_FQ)%all_square)  & 
            scheme%T_searcher(DO_FQ)%shape = SHAPE_TRIANGULAR
         IF (.NOT.scheme%T_searcher(DO_NN)%all_square)  & 
            scheme%T_searcher(DO_NN)%shape = SHAPE_TRIANGULAR
      END IF

      ! now build energy additively
      energy = zero
      ! get moments
      CALL mm_init_driver(scheme,WORK,LWORK)
      ! get potential
      CALL mm_get_raw_Vff(scheme,LHS_mms%paras,RHS_mms,Vff)
      ! get energy
      IF (scheme%pack_LHS) THEN
         CALL mm_get_E_from_pkd_Vff(scheme,LHS_mms,Vff,energy,text)
      ELSE
         CALL mm_get_E_from_Vff(scheme,LHS_mms,Vff,energy,text)
      END IF
      CALL mm_free_driver
     
   END SUBROUTINE mm_get_E_via_raw_potentials

!-------------------------------------------------------------------------------

!> \brief Caclulate energy from the raw potentials
!> \author Mark Watson
!> \date 2002
!> \param scheme   Contains data for mm-calculation
!> \param J_matrix The mm Coulomb-matrix contribution
!> \param energy   The mm-energy contribution
!> \param text     The energy-contribution type being calculated
!> \param work     Static memory string (if not ALLOCATION_WITHOUT_WORK)
!> \param lwork    Remaining length of work
   SUBROUTINE mm_get_J_via_raw_potentials(scheme,J_matrix,energy,txt,WORK,LWRK)

      USE mm_Vff_driver,    ONLY: mm_get_raw_Vff
      USE mm_Vff_processor, ONLY: mm_get_J_from_Vff,       &
                                  mm_get_J_from_pkd_Vff,   &
                                  mm_get_E_from_pkd_Vff,   &
                                  mm_get_E_from_Vff
      USE mm_qlm_processor, ONLY: mm_factor_in_dens
      USE density_fitting
      USE LSTIMING, ONLY: LSTIMER

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(INOUT) :: scheme
      REAL(REALK),        INTENT(OUT)   :: J_matrix(:,:)
      REAL(REALK),        INTENT(OUT)   :: energy
      CHARACTER(*),       INTENT(OUT)   :: txt
      INTEGER,      INTENT(IN)    :: LWRK
      REAL(REALK),        INTENT(INOUT) :: WORK(:)

      REAL(REALK) :: TS,TE

!      CALL LSTIMER('START ',TS,TE,-1)

      ! we only have the density on RHS when getting J-matrix
      ! via far-field potential
      scheme%LHS_dens = .FALSE.
      scheme%RHS_dens = .TRUE.
      ! get moments
      CALL mm_init_driver(scheme,WORK,LWRK)


!      CALL LSTIMER('mm-ini',TS,TE,-1)

      ! get potential
      CALL mm_get_raw_Vff(scheme,LHS_mms%paras,RHS_mms,Vff)

!      CALL LSTIMER('mm-Vff',TS,TE,-1)
      ! get J-matrix
      J_matrix = zero
      energy = zero
      IF (scheme%pack_LHS) THEN
         CALL mm_get_J_from_pkd_Vff(scheme,LHS_mms,Vff,J_matrix)
         ! get energy after factoring in density to LHS
         CALL mm_factor_in_dens(LHS_mms%dens,LHS_mms%qlm_T)
         CALL mm_get_E_from_pkd_Vff(scheme,LHS_mms,Vff,energy,txt)
      ELSE
         CALL mm_get_J_from_Vff(scheme,LHS_mms,Vff,J_matrix)
         ! get energy after factoring in density to LHS
         CALL mm_factor_in_dens(LHS_mms%dens,LHS_mms%qlm_T)
         CALL mm_get_E_from_Vff(scheme,LHS_mms,Vff,energy,txt)
      END IF
!      CALL LSTIMER('mm-J-E',TS,TE,-1)

      CALL mm_free_driver

   END SUBROUTINE mm_get_J_via_raw_potentials

!-------------------------------------------------------------------------------
! When not combining moments, it is possible to save a factor of two in the
! contractions by avoiding double-counting.  This is only possible when we
! do not build the potential (which mixes RHS contributions).
! i.e. used in NN phase to obtain factor of 2 speed-up.
!
!   SUBROUTINE mm_get_J_directly(scheme,J_matrix,WORK,LWORK)
!
!      USE mm_Vff_driver,      ONLY: mm_get_raw_Vff
!      USE mm_Vff_processor,   ONLY: mm_get_J_from_Vff
!      USE mm_J_driver,        ONLY: mm_drive_NN_J_builder
!      USE mm_qlm_processor,   ONLY: mm_factor_in_dens
!
!      IMPLICIT NONE
!      TYPE(scheme_paras), INTENT(INOUT) :: scheme
!      REAL(REALK),        INTENT(OUT)   :: J_matrix(:,:)
!      INTEGER,      INTENT(IN)    :: LWORK
!      REAL(REALK),        INTENT(INOUT) :: WORK(:) ! input value of LWORK
!
!      REAL(REALK), POINTER :: dens_ptr(:)
!
!      ! set up T-pair generator algorithm
!      IF (scheme%LHS_mm_range == scheme%RHS_mm_range) THEN
!         IF (.NOT.scheme%T_searcher(DO_FQ)%all_square)  & 
!            scheme%T_searcher(DO_FQ)%shape = SHAPE_TRIANGULAR
!         IF (.NOT.scheme%T_searcher(DO_NN)%all_square)  & 
!            scheme%T_searcher(DO_NN)%shape = SHAPE_TRIANGULAR
!      END IF
!
!      J_matrix = zero
!
!      ! get raw moments for nearest neighbour phase without densities
!      scheme%LHS_dens = .FALSE.
!      scheme%RHS_dens = .FALSE.
!      CALL mm_init_driver(scheme,WORK,LWORK)
!
!      ! get NN part (or total if FQ algorithm) without going via Vff
!      IF (scheme%LHS_mm_range /= scheme%RHS_mm_range) THEN
!         CALL LSQUIT('ranges must be equal to get_J_directly')
!      END IF
!      ! set up density pointer for J-matrix build
!      dens_ptr => RHS_mms%dens(:)  ! LHS .eq. RHS
!      CALL mm_drive_NN_J_builder(scheme,LHS_mms,RHS_mms,dens_ptr,J_matrix)
!
!      ! reset moments for far-field phase
!      ! we add the density to RHS when getting J-matrix via far-field potential
!!FIMXE: add checks because this wont work if we compact RHS first
!      CALL mm_factor_in_dens(RHS_mms%dens,RHS_mms%qlm_T)
!      CALL mm_factor_in_dens(RHS_mms%dens,RHS_mms%qlm_W)
!
!      ! get potential for FF work
!      IF (scheme%algorithm /= DO_FQ) THEN
!         CALL mm_get_raw_Vff(scheme,LHS_mms%paras,RHS_mms,Vff,skip_NN=.TRUE.)
!         CALL mm_get_J_from_Vff(scheme,LHS_mms,Vff,J_matrix)
!      END IF
!
!      CALL mm_free_driver
!
!   END SUBROUTINE mm_get_J_directly
!
!-------------------------------------------------------------------------------

!> \brief Determines cacluation type based on input text-string
!> \author Mark Watson
!> \date 2002
!> \param n_el Specify type of calculation
   INTEGER FUNCTION mm_get_int_type(n_el)

      IMPLICIT NONE
      CHARACTER(6), INTENT(IN) :: n_el

      SELECT CASE(n_el)
      CASE('ONE_EL')
         mm_get_int_type = INTTYPE_ONE_EL
      CASE('TWO_EL')
         mm_get_int_type = INTTYPE_TWO_EL
      CASE('FULL_J')
         mm_get_int_type = INTTYPE_FULL_J
      CASE('NUC_AT')
         mm_get_int_type = INTTYPE_NUC_AT
      CASE('NUC_EL')
         mm_get_int_type = INTTYPE_NUC_EL
      CASE DEFAULT
         mm_get_int_type = -1
      END SELECT

   END FUNCTION mm_get_int_type

!-------------------------------------------------------------------------------

END MODULE mm_main_driver

!===============================================================================

!> \brief Caclulate energy from the raw potentials
!> \author Mark Watson
!> \date 2002
!> \param ERI_energy The mm-energy contribution
!> \param n_el       Specifying type of contribution
!> \param work       Static memory string (if not ALLOCATION_WITHOUT_WORK)
!> \param lwork      Remaining length of work
SUBROUTINE mm_get_energy(ERI_energy,n_el,WORK,LWRK)

   USE mm_global_paras_mod
   USE mm_stats_mod
   USE mm_memory_manager_mod, ONLY: mm_init_mem_man, mm_close_mem_man
   USE mm_scheme_builder, ONLY: mm_get_scheme
   USE mm_interface_mod,  ONLY: mm_connect_interface, mm_disconnect_interface
   USE mm_main_driver,    ONLY: mm_get_E_via_raw_potentials, mm_get_int_type
   use mm_mem

   IMPLICIT NONE
   REAL(REALK),   INTENT(IN)    :: ERI_energy
   CHARACTER(6),  INTENT(IN)    :: n_el
   INTEGER, INTENT(IN)    :: LWRK
   REAL(REALK),   INTENT(INOUT) :: WORK(LWRK)

   REAL(REALK)   :: energy
   CHARACTER(36) :: E_text
   LOGICAL       :: A, B
   INTEGER :: INUM, KWRK, integraltypeint
   TYPE(scheme_paras), POINTER :: scheme

   KWRK = 1

   CALL mm_init_mem_man
   CALL mm_get_scheme(scheme)
   CALL mm_init_stats(scheme,GET_ENERGY,DO_NOGRADIENT)

   CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)

   integraltypeint = mm_get_int_type(n_el)
   SELECT CASE (integraltypeint)
   CASE (INTTYPE_ONE_EL)
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = NUCLEAR_ONLY 
   CASE (INTTYPE_TWO_EL)
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = ELECTRONIC_ONLY 
   CASE (INTTYPE_NUC_AT)
      scheme%LHS_mm_range = NUCLEAR_ONLY
      scheme%RHS_mm_range = NUCLEAR_ONLY 
   CASE DEFAULT
      CALL LSQUIT ('can only do n-n e-e or e-n energy calculation!',-1)
   END SELECT

   CALL mm_get_E_via_raw_potentials(scheme,energy,E_text,              &
                                    WORK(KWRK:LWRK),LWRK-KWRK+1)

   !PRINT '(X,A," = ",F22.12)', E_text, energy
!   write(LUPRI,'(X,A," = ",F22.12)') E_text, energy

   A = (scheme%LHS_mm_range == ELECTRONIC_ONLY)
   B = (scheme%RHS_mm_range == ELECTRONIC_ONLY)
   IF (A.AND.B) THEN
     !PRINT '(" error wrt ERI e-e energy   =",E11.3)', ABS(energy-ERI_energy)
     write(LUPRI,'(" error wrt ERI e-e energy =",E11.3)')ABS(energy-ERI_energy)
   END IF

   CALL mm_print_stats
   CALL mm_disconnect_interface
   CALL mm_close_mem_man

END SUBROUTINE mm_get_energy


!===============================================================================
!> \brief Driver for calculating the mm-gradient contribution
!> \author Andreas Krapp
!> \date 2010
!> \param n_el   Specifies contribution to be calculated
!> \param grad   The moleular gradient mm-contribution
!> \param numnuc The number of nuclei
!> \param work   Static memory string (if not ALLOCATION_WITHOUT_WORK)
!> \param lwork  Remaining length of work
SUBROUTINE mm_get_grad(n_el,grad,numnuc,WORK,LWRK)

   USE mm_global_paras_mod
   USE mm_stats_mod
   USE mm_memory_manager_mod, ONLY: mm_init_mem_man, mm_close_mem_man
   USE mm_scheme_builder,     ONLY: mm_get_scheme
   USE mm_interface_mod,      ONLY: mm_connect_interface, mm_disconnect_interface
   USE mm_main_driver,        ONLY: mm_get_grad_via_raw_potentials, mm_get_int_type
   USE density_fitting
   use mm_mem

   IMPLICIT NONE
   CHARACTER(6),  INTENT(IN)    :: n_el
   REAL(REALK),   INTENT(INOUT) :: grad(3*numnuc)
   INTEGER, INTENT(IN)    :: LWRK
   REAL(REALK),   INTENT(INOUT) :: WORK(LWRK)

   TYPE(scheme_paras), POINTER  :: scheme
   REAL(REALK)                  :: energy
   LOGICAL                      :: A, B
   INTEGER                :: INUM, KWRK, integraltypeint, numnuc

   ! determine which type of calculation
   integraltypeint = mm_get_int_type(n_el)
   IF( (integraltypeint .NE. INTTYPE_TWO_EL ).AND. &
     & (integraltypeint .NE. INTTYPE_NUC_AT ).AND. &
     & (integraltypeint .NE. INTTYPE_NUC_EL )) THEN
     CALL LSQUIT ('require 1 or 2 electron gradient build!',-1)
   END IF

   ! initialisations
   CALL mm_init_mem_man
   CALL mm_get_scheme(scheme)
   CALL mm_init_stats(scheme,.false.,DO_GRADIENT)
   KWRK = 1
   energy = 0.0D0
   grad(:) = 0.0D0

   IF ( integraltypeint .EQ. INTTYPE_NUC_AT ) THEN

      call lsquit('mm gradients with option NUC_AT not tested for lsint',-1)

      IF (fit%AUX) THEN
         call lsquit('density fitting not for mm gradients with option NUC_AT',-1)
      END IF

      !
      !  get gradient contribution (rho'|nuc)
      !
      CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = NUCLEAR_ONLY
      fit%LHS_AUX = .false.
      fit%RHS_AUX = .false.
      CALL mm_get_grad_via_raw_potentials(scheme,grad,energy,numnuc,WORK(KWRK:LWRK),LWRK-KWRK+1)
      CALL mm_disconnect_interface
      !
      !  get gradient contribution (nuc'|nuc)
      !
      CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)
      scheme%LHS_mm_range = NUCLEAR_ONLY
      scheme%RHS_mm_range = NUCLEAR_ONLY
      fit%LHS_AUX = .false.
      fit%RHS_AUX = .false.
      CALL mm_get_grad_via_raw_potentials(scheme,grad,energy,numnuc,WORK(KWRK:LWRK),LWRK-KWRK+1)
      CALL mm_disconnect_interface
      !
      !  get gradient contribution (nuc'|rho)
      !
      CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)
      scheme%LHS_mm_range = NUCLEAR_ONLY
      scheme%RHS_mm_range = ELECTRONIC_ONLY
      fit%LHS_AUX = .false.
      fit%RHS_AUX = .false.
      CALL mm_get_grad_via_raw_potentials(scheme,grad,energy,numnuc,WORK(KWRK:LWRK),LWRK-KWRK+1)
      CALL mm_disconnect_interface

   ELSEIF ( integraltypeint .EQ. INTTYPE_NUC_EL) THEN
      ! 
      ! nucleus-nucleus, nucleus-electron and electron-electron gradient contributions at once
      !
      if(fit%AUX) THEN
         !
         ! gradient contribution (rho'+nuc'|rho_fit+nuc)
         !
         stat_reorder_mom = .true.
         CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)
         scheme%LHS_mm_range = ALL_MOMENTS
         scheme%RHS_mm_range = ALL_MOMENTS
         fit%LHS_AUX = .false.
         fit%RHS_AUX = .true.
         CALL mm_get_grad_via_raw_potentials(scheme,grad,energy,numnuc,WORK(KWRK:LWRK),LWRK-KWRK+1)
         CALL mm_disconnect_interface
         !
         ! gradient contribution (rho_fit'+nuc'|rho-rho_fit)
         ! 
         stat_reorder_mom = .false.
         CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)
         scheme%LHS_mm_range = ALL_MOMENTS
         scheme%RHS_mm_range = REG_MIN_AUX
         fit%LHS_AUX  = .true.
         fit%RHS_AUX  = .true.
         CALL mm_get_grad_via_raw_potentials(scheme,grad,energy,numnuc,WORK(KWRK:LWRK),LWRK-KWRK+1)
         CALL mm_disconnect_interface
      else
         !
         ! gradient contribution (rho'+nuc'|rho+nuc)
         !
         CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)
         scheme%LHS_mm_range = ALL_MOMENTS
         scheme%RHS_mm_range = ALL_MOMENTS
         CALL mm_get_grad_via_raw_potentials(scheme,grad,energy,numnuc,WORK(KWRK:LWRK),LWRK-KWRK+1)
         CALL mm_disconnect_interface
         grad(:) = grad(:)*2.0D0
      end if
   ELSEIF ( integraltypeint .EQ. INTTYPE_TWO_EL ) THEN
      !
      !  get gradient contribution (rho'|rho) 
      !  in case of density fitting get contribution (rho'|rho_fit)
      !
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = ELECTRONIC_ONLY
      CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)
      CALL mm_get_grad_via_raw_potentials(scheme,grad,energy,numnuc,WORK(KWRK:LWRK),LWRK-KWRK+1)
      CALL mm_disconnect_interface
!
      if(fit%AUX) THEN
         !
         !     get df-gradient contribution (rho_fit'|rho-rho_fit)
         !
         CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)
         scheme%LHS_mm_range = ELECTRONIC_ONLY
         scheme%RHS_mm_range = REG_MIN_AUX
         fit%LHS_AUX  = .true.
         fit%RHS_AUX  = .true.
         CALL mm_get_grad_via_raw_potentials(scheme,grad,energy,numnuc,WORK(KWRK:LWRK),LWRK-KWRK+1)
         CALL mm_disconnect_interface
      endif
   END IF

   CALL mm_print_stats
   CALL mm_close_mem_man

END SUBROUTINE mm_get_grad

!-------------------------------------------------------------------------------

!> \brief Driver for calculating the mm Coulomb-matrix contribution
!> \author Mark Watson
!> \date 2010
!> \param n_el        Specifies contribution to be calculated
!> \param Fock_matrix The matrix to add contribution to
!> \param Nfock1      The first dimension of Fock_matrix
!> \param Nfock2      The second dimension of Fock_matrix
!> \param work        Static memory string (if not ALLOCATION_WITHOUT_WORK)
!> \param lwork       Remaining length of work
SUBROUTINE mm_get_J_matrix(n_el,Fock_matrix,Nfock1,Nfock2,WORK,LWRK)

   USE mm_global_paras_mod
   USE mm_stats_mod
   USE mm_memory_manager_mod, ONLY: mm_init_mem_man, mm_close_mem_man
   USE mm_scheme_builder, ONLY: mm_get_scheme
   USE mm_interface_mod,  ONLY: mm_connect_interface, mm_disconnect_interface
   USE mm_main_driver,    ONLY: mm_get_J_via_raw_potentials, mm_get_int_type
   USE LSTIMING, ONLY: LSTIMER
   use mm_mem
#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
   use mpi
   IMPLICIT NONE
#else
   IMPLICIT NONE
   INCLUDE 'mpif.h'
#endif
#else
   IMPLICIT NONE
#endif
   CHARACTER(6),  INTENT(IN)    :: n_el
   REAL(REALK),   INTENT(INOUT) :: Fock_matrix(Nfock1,Nfock2) 
   INTEGER, INTENT(IN)    :: LWRK, Nfock1,Nfock2
   REAL(REALK),   INTENT(INOUT) :: WORK(LWRK)

!FIXME: dont want to store a separate J_matrix if we can pass Fock_matrix
   REAL(REALK), POINTER :: J_matrix(:,:)
   REAL(REALK)        :: energy
   INTEGER      :: rank, ierr, integraltypeint, KWRK
   CHARACTER(36)      :: E_text
   LOGICAL            :: branch_free_flag
   TYPE(scheme_paras), POINTER :: scheme
   DOUBLE PRECISION SECOND, T1
   REAL(REALK) :: TS,TE

   CALL LSTIMER('START ',TS,TE,-1)

   T1 = SECOND()

   CALL mm_init_stats(scheme,GET_JMAT,DO_NOGRADIENT)
!   stat_NBAST = NBAST

   CALL mm_init_mem_man
   integraltypeint = mm_get_int_type(n_el)
   stat_1el_2el_or_full = mm_get_int_type(n_el)

#if defined(VAR_MPI)&&0
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,IERR)
   CALL mm_getj_wake_up_slaves()
   IF(rank==0) CALL mm_get_scheme(scheme)
   CALL mm_getj_sync_slaves
#else
   CALL mm_get_scheme(scheme)
#endif
   KWRK = 1
!   CALL TIMER('mm-ini1',TS,TE,-1)
   CALL mm_connect_interface(scheme,WORK,KWRK,LWRK)
!   CALL TIMER('mm-conn',TS,TE,-1)

   branch_free_flag = scheme%branch_free
!   ALLOCATE(J_matrix(Nfock1,Nfock2))
   call mem_alloc_fmm(J_matrix,Nfock1,Nfock2)
   SELECT CASE (integraltypeint)
   CASE (INTTYPE_ONE_EL)
      ! always run 1-el part without "branch-free" algorithm
      scheme%branch_free = .FALSE.
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = NUCLEAR_ONLY
      CALL mm_get_J_via_raw_potentials(scheme,J_matrix,energy,E_text,  &
                                       WORK(KWRK:LWRK),LWRK-KWRK+1)
      IF ( MM_STANDALONE ) THEN
!         WRITE(LUPRI,'(X,A," = ",E20.12)') E_text, energy
      END IF
      !WRITE(LUPRI,'(A,E20.12)') ' MAX(ABS(J_matrix)) =',MAXVAL(ABS(J_matrix))
   CASE (INTTYPE_TWO_EL)
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = ELECTRONIC_ONLY
      CALL mm_get_J_via_raw_potentials(scheme,J_matrix,energy,E_text,  &
                                       WORK(KWRK:LWRK),LWRK-KWRK+1)
!      WRITE(LUPRI,'(X,A," = ",E20.12)') E_text, energy
      !WRITE(LUPRI,'(A,E20.12)') ' MAX(ABS(J_matrix)) =',MAXVAL(ABS(J_matrix))
   CASE (INTTYPE_FULL_J)
      scheme%LHS_mm_range = ELECTRONIC_ONLY
      scheme%RHS_mm_range = ALL_MOMENTS
      CALL mm_get_J_via_raw_potentials(scheme,J_matrix,energy,E_text,  &
                                       WORK(KWRK:LWRK),LWRK-KWRK+1)
!      WRITE(LUPRI,'(X,A," = ",E20.12)') E_text, energy
   CASE DEFAULT
      CALL LSQUIT ('require 1, 2, or full J_matrix build!',-1)
   END SELECT
!   CALL TIMER('mm-getJ',TS,TE,-1)
   scheme%branch_free = branch_free_flag

#if defined(VAR_MPI) &&0
   CALL mm_getj_collect(SECOND()-T1)
#endif

   Fock_matrix = Fock_matrix + J_matrix
!   DEALLOCATE(J_matrix)
   call mem_dealloc_fmm(J_matrix)
   !write(6,*) "J_matrix from MM:"
   !CALL LS_OUTPUT(J_matrix,1,Nfock1,1,Nfock2,Nfock1,Nfock2,1,6)

   CALL mm_print_stats
   CALL mm_disconnect_interface
   CALL mm_close_mem_man

!   CALL TIMER('mm-close',TS,TE,-1)
END SUBROUTINE mm_get_J_matrix

!-------------------------------------------------------------------------------

#if defined(VAR_MPI)&&0
SUBROUTINE mm_getj_sync_slaves
   ! can be split into two parts: syncing scheme, and other parts.
#ifdef USE_MPI_MOD_F90
   use mpi
   IMPLICIT NONE
#else
   IMPLICIT NONE
   INCLUDE 'mpif.h'
#endif
   INTEGER ierr
   INTEGER, PARAMETER :: master = 0

   CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNUM, IERR)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NNODES,IERR)
   CALL MPI_Bcast(integraltypeint,1,MPI_INTEGER,master,MPI_COMM_WORLD,IERR)
   CALL bcast_scheme_paras(master)
   CALL bcast_T_searcher_type(master,scheme%T_searcher)
   CALL bcast_T_contract_schm(master,scheme%T_con)
   CALL bcast_W_contract_schm(master,scheme%W_con)
END SUBROUTINE mm_getj_sync_slaves

! note that broadcasting system size is redundant, it is not an
! input PARAMETER.
SUBROUTINE bcast_scheme_paras(master)
#ifdef USE_MPI_MOD_F90
   use mpi
   IMPLICIT NONE
#else
   IMPLICIT NONE
   INCLUDE 'mpif.h'
#endif
   INTEGER, INTENT(IN) :: master
   INTEGER intEx, logEx, doubleEx, ierr, theType
   INTEGER, DIMENSION(4) :: bllen, types, disps
   CALL MPI_TYPE_EXTENT(MPI_INTEGER,intex,ierr)
   CALL MPI_TYPE_EXTENT(MPI_LOGICAL,logex,ierr)
   CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,doubleEX,ierr)
   
   types = (/ MPI_INTEGER, MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_PRECISION/)
   bllen = (/ 2,           4,          6,          4                   /)
   disps(1) = 0; 
   disps(2) = disps(1) + bllen(1)*intEx
   disps(3) = disps(2) + bllen(2)*logEx
   disps(4) = disps(3) + bllen(3)*intEx
   CALL MPI_Type_Struct(SIZE(bllen),bllen,disps,types,theType, ierr)
   CALL MPI_TYPE_COMMIT(theType,ierr)
   CALL MPI_Bcast(scheme,1,theType,master,MPI_COMM_WORLD,IERR)
   CALL MPI_TYPE_FREE(theType,ierr)
END SUBROUTINE bcast_scheme_paras

SUBROUTINE bcast_T_searcher_type(master,searcher)
#ifdef USE_MPI_MOD_F90
   use mpi
   IMPLICIT NONE
#else
   IMPLICIT NONE
   INCLUDE 'mpif.h'
#endif
   INTEGER, INTENT(IN)   :: master
   TYPE(T_searcher_type) :: searcher(5)

   INTEGER intEx, logEx, ierr, theType
   INTEGER, DIMENSION(2) :: bllen, types, disps

   CALL MPI_TYPE_EXTENT(MPI_INTEGER,intex,ierr)
   CALL MPI_TYPE_EXTENT(MPI_LOGICAL,logex,ierr)
   
   types = (/ MPI_INTEGER, MPI_LOGICAL /)
   bllen = (/ 2,           1           /)
   disps(1) = 0; 
   disps(2) = disps(1) + bllen(1)*intEx
   CALL MPI_Type_Struct(SIZE(bllen),bllen,disps,types,theType, IERR)
   CALL MPI_TYPE_COMMIT(theType,IERR)
   CALL MPI_Bcast(searcher,5,theType,master,MPI_COMM_WORLD,IERR)
   CALL MPI_TYPE_FREE(theType,ierr)
END SUBROUTINE bcast_T_searcher_type

SUBROUTINE bcast_T_contract_schm(master,contractor)
#ifdef USE_MPI_MOD_F90
   use mpi
   IMPLICIT NONE
#else
   IMPLICIT NONE
   INCLUDE 'mpif.h'
#endif
   INTEGER, INTENT(IN)                  :: master
   TYPE(T_contract_schm), INTENT(INOUT) :: contractor
   INTEGER ierr
   CALL MPI_Bcast(contractor,6,MPI_INTEGER,master,MPI_COMM_WORLD,IERR)
END SUBROUTINE bcast_T_contract_schm

SUBROUTINE bcast_W_contract_schm(master,w_con)
#ifdef USE_MPI_MOD_F90
   use mpi
   IMPLICIT NONE
#else
   IMPLICIT NONE
   INCLUDE 'mpif.h'
#endif
   INTEGER, INTENT(In)                  :: master
   TYPE(W_contract_schm), INTENT(INOUT) :: w_con
   integer ierr
   CALL MPI_Bcast(w_con,3,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
END SUBROUTINE bcast_W_contract_schm

SUBROUTINE mm_getj_collect(tm)
#ifdef USE_MPI_MOD_F90
   use mpi
   IMPLICIT NONE
#else
   IMPLICIT NONE
   INCLUDE 'mpif.h'
#endif
   INTEGER IERR, MASTER,I
   REAL(REALK) :: Work(NBAST*NBAST)
   DOUBLE PRECISION, INTENT(IN) :: tm
   DOUBLE PRECISION  :: tms(NNODES)
   
   MASTER = 0
   CALL DCOPY(NBAST*NBAST,J_matrix, 1, Work, 1)
   CALL MPI_Reduce(work,J_matrix,NBAST*NBAST,MPI_DOUBLE_PRECISION,&
        & MPI_SUM, MASTER, MPI_COMM_WORLD, IERR)
#if 1
   CALL MPI_GATHER(tm,1,MPI_DOUBLE_PRECISION,&
        &          tms,1,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,IERR)
   WRITE(LUPRI,"(I3,'Node times:',/,20(10F7.2,/))") MYNUM,&
        & (tms(i),i=1,NNODES)
#endif
END SUBROUTINE mm_getj_collect
#endif

#ifndef GFORTRAN
FUNCTION SECOND()

   USE mm_global_paras_mod
   IMPLICIT NONE
   REAL(REALK) :: SECOND
   INTEGER :: count, count_rate, count_max, tmp
   INTEGER, SAVE :: parity=0
   INTEGER, SAVE :: last_count=0

   IF (parity == 0) parity = 1  ! first call to SECOND()
   parity = -parity
   CALL SYSTEM_CLOCK(count,count_rate,count_max)
   ! must check if count has "cycled"
   tmp = count
   IF ( parity == 1 ) THEN  ! second of a paired call to SECOND()
      IF ( count < last_count ) tmp = count + (count_max - last_count)
   END IF
   SECOND = REAL(tmp,REALK)/REAL(count_rate,REALK)
   last_count = count 

 END FUNCTION SECOND
#endif

!-------------------------------------------------------------------------------
