MODULE mm_T_pair_mould

   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_init_T_pair_mould,           &
             mm_close_T_pair_mould

   INTEGER, SAVE :: LHS_LMAX, RHS_LMAX, LEXTRA
   LOGICAL,       SAVE :: mm_dynamic_LMAX_on
    ! flag to test initialisation
   CHARACTER(11), SAVE :: mm_init_mould

CONTAINS

!-------------------------------------------------------------------------------
! This routine directs the saving of functions in mm_proc_selector.c
! relevant to the formation of a T-pair entity.
! There are 4 parts to the making of a T-pair once raw data has been
! suppplied from the classical interaction search algorithm.
! These are all called consecutively from the C-code via
! mm_stored_t_pair_mould.
! The 4 functions saved are selected here and depend on whether
! pure boxed moments, or unboxed moments, are being contracted.

   SUBROUTINE mm_init_T_pair_mould(scheme,pair_type)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER,      INTENT(IN) :: pair_type
      EXTERNAL mm_store_t_pair_mould1   ! raw/box dependent variables
      EXTERNAL mm_store_t_pair_mould2   ! dynamic/non-dynamic LMAX (LHS)
      EXTERNAL mm_store_t_pair_mould3   ! dynamic/non-dynamic LMAX (RHS)
      EXTERNAL mm_store_t_pair_mould4   ! common to all T-pair builds

!$OMP MASTER
      mm_init_mould      = 'initialised' 
!$OMP END MASTER
      mm_dynamic_LMAX_on = scheme%dynamic_LMAX_on
      LEXTRA             = scheme%LEXTRA

      CALL mm_store_t_pair_mould2(mm_set_LHS_std_LMAX)
      CALL mm_store_t_pair_mould3(mm_set_RHS_std_LMAX)
      CALL mm_store_t_pair_mould4(mm_set_T_pair_basics)
      SELECT CASE (pair_type)
      CASE (LHS_raw_RHS_raw)
         LHS_LMAX = scheme%raw_LMAX   
         RHS_LMAX = scheme%raw_LMAX   
         CALL mm_store_t_pair_mould1(mm_set_RR_paras)
         IF (mm_dynamic_LMAX_on) THEN
            CALL mm_store_t_pair_mould2(mm_set_LHS_dyn_LMAX)
            CALL mm_store_t_pair_mould3(mm_set_RHS_dyn_LMAX)
         END IF
      CASE (LHS_box_RHS_box)
         LHS_LMAX = scheme%trans_LMAX   
         RHS_LMAX = scheme%trans_LMAX   
         CALL mm_store_t_pair_mould1(mm_set_BB_paras)
      CASE DEFAULT
         CALL LSQUIT ('cannot recognise T_pair type!',-1)
      END SELECT

   END SUBROUTINE mm_init_T_pair_mould

!-------------------------------------------------------------------------------

   SUBROUTINE mm_close_T_pair_mould

      IMPLICIT NONE
!$OMP MASTER
      IF (mm_init_mould /= 'initialised') STOP 'mm_T_pair_mould initialisation'
      mm_init_mould = ' '
!$OMP END MASTER
      LHS_LMAX = 0
      RHS_LMAX = 0
      LEXTRA   = 0
      mm_dynamic_LMAX_on = .FALSE.
 
   END SUBROUTINE mm_close_T_pair_mould

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_T_pair_basics(LHS,RHS,id,weight,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: LHS, RHS
      TYPE(LHS_RHS_type),  INTENT(IN)  :: id
      INTEGER,       INTENT(IN)  :: weight
      TYPE(T_pair_single), INTENT(INOUT) :: T_pair

      T_pair%N_or_T = 'N' ! not used for T contractions, so null
      ! init. normalisation scalar for other options
      T_pair%paras%ratio = one 
      ! here, T_pair%r_ab is the unnormalised vector

      ! weighting to account for double-counting in T-pair generation
      T_pair%paras%weight = weight

      T_pair%LMAX = MAX(T_pair%paras%LHS_LMAX,T_pair%paras%RHS_LMAX)
      T_pair%lm_max = (1+T_pair%LMAX)**2

   END SUBROUTINE mm_set_T_pair_basics

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_RR_paras(LHS,RHS,id,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: LHS, RHS
      TYPE(LHS_RHS_type),  INTENT(IN)  :: id
      TYPE(T_pair_single), INTENT(OUT) :: T_pair

      ! interaction vector for building T-matrix
      T_pair%r_ab = RHS%raw_paras(id%RHS)%cntr - LHS%raw_paras(id%LHS)%cntr
      ! indices to map back to actual moments in separate array
      T_pair%paras%LHS_id = LHS%raw_paras(id%LHS)%id 
      T_pair%paras%RHS_id = RHS%raw_paras(id%RHS)%id 
      ! check that paras:moments mapping was built
      IF (T_pair%paras%LHS_id == 0) CALL LSQUIT('LHS paras:moments mapping bad!',-1)
      IF (T_pair%paras%RHS_id == 0) CALL LSQUIT('RHS paras:moments mapping bad!',-1)

   END SUBROUTINE mm_set_RR_paras

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_BB_paras(LHS,RHS,id,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: LHS, RHS
      TYPE(LHS_RHS_type),  INTENT(IN)  :: id
      TYPE(T_pair_single), INTENT(OUT) :: T_pair

      ! interaction vector for building T-matrix
      T_pair%r_ab = RHS%box_paras(id%RHS)%cntr - LHS%box_paras(id%LHS)%cntr
      ! indices to map back to actual moments in separate array
      T_pair%paras%LHS_id = LHS%box_paras(id%LHS)%id 
      T_pair%paras%RHS_id = RHS%box_paras(id%RHS)%id 
      ! check that paras:moments mapping was built
      IF (T_pair%paras%LHS_id == 0) CALL LSQUIT('LHS paras:moments mapping bad!',-1)
      IF (T_pair%paras%RHS_id == 0) CALL LSQUIT('RHS paras:moments mapping bad!',-1)

   END SUBROUTINE mm_set_BB_paras

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_LHS_std_LMAX(X,Y,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: X        ! dummy variables
      TYPE(LHS_RHS_type),  INTENT(IN)  :: Y        ! dummy variables
      TYPE(T_pair_single), INTENT(OUT) :: T_pair

      T_pair%paras%LHS_LMAX = LHS_LMAX

   END SUBROUTINE mm_set_LHS_std_LMAX

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_RHS_std_LMAX(X,Y,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: X        ! dummy variables
      TYPE(LHS_RHS_type),  INTENT(IN)  :: Y        ! dummy variables
      TYPE(T_pair_single), INTENT(OUT) :: T_pair

      T_pair%paras%RHS_LMAX = RHS_LMAX

   END SUBROUTINE mm_set_RHS_std_LMAX

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_LHS_dyn_LMAX(LHS,i,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: LHS
      TYPE(LHS_RHS_type),  INTENT(IN)  :: i
      TYPE(T_pair_single), INTENT(OUT) :: T_pair

      T_pair%paras%LHS_LMAX = MIN(LHS_LMAX,LHS%raw_paras(i%LHS)%Lmin+LEXTRA)

   END SUBROUTINE mm_set_LHS_dyn_LMAX

!-------------------------------------------------------------------------------

   SUBROUTINE mm_set_RHS_dyn_LMAX(RHS,i,T_pair)

      IMPLICIT NONE
      TYPE(gen_mm_paras),  INTENT(IN)  :: RHS
      TYPE(LHS_RHS_type),  INTENT(IN)  :: i
      TYPE(T_pair_single), INTENT(OUT) :: T_pair

      T_pair%paras%RHS_LMAX = MIN(RHS_LMAX,RHS%raw_paras(i%RHS)%Lmin+LEXTRA)

   END SUBROUTINE mm_set_RHS_dyn_LMAX

!-------------------------------------------------------------------------------

END MODULE mm_T_pair_mould

!===============================================================================

MODULE mm_FQ_T_loops

   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_init_FQ_loops, &
             mm_exe_FQ_loops
   ! dimensions of FQ loops - must be initialised
   TYPE(LHS_RHS_type), SAVE :: dims

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_FQ_loops(LHS,RHS,pair_type,do_box_pretesting)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      INTEGER,      INTENT(IN) :: pair_type
      LOGICAL,            INTENT(IN) :: do_box_pretesting

      SELECT CASE (pair_type)
      CASE (LHS_raw_RHS_raw)
         IF (do_box_pretesting) THEN
            dims%LHS = SIZE(LHS%box_paras)
            dims%RHS = SIZE(RHS%box_paras)
         ELSE
            dims%LHS = SIZE(LHS%raw_paras)
            dims%RHS = SIZE(RHS%raw_paras)
         END IF
      CASE (LHS_box_RHS_box)
         dims%LHS = SIZE(LHS%box_paras)
         dims%RHS = SIZE(RHS%box_paras)
      CASE DEFAULT
         CALL LSQUIT ('cannot reconcile requested T_pair type!',-1)
      END SELECT

   END SUBROUTINE mm_init_FQ_loops

!-------------------------------------------------------------------------------
! FIXME: could get rid of triangular loop here, as we have a RHS>LHS test
! in mm_single_test instead!!!

   SUBROUTINE mm_exe_FQ_loops(LHS,RHS,scheme,phase,shape)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      INTEGER,      INTENT(IN) :: phase,shape
      TYPE(scheme_paras), INTENT(IN) :: scheme
      LOGICAL, EXTERNAL   :: mm_included_pair
      TYPE(LHS_RHS_type)  :: id
      INTEGER       :: weight, i,j,start,step
      REAL(REALK)         :: ntested
      EXTERNAL mm_test_and_buffer_T_pair

!FIXME: parallelisation broken
!      start = MYNUM+1
!      step  = NNODES

      SELECT CASE (shape)
      CASE (SHAPE_SQUARE) 
         CALL square_loop(scheme,phase)
      CASE (SHAPE_TRIANGULAR) 
         IF (dims%RHS /= dims%LHS) CALL LSQUIT('TRI_FQ_LOOP invalid!',-1)
         CALL triangular_loop
      CASE DEFAULT
         CALL LSQUIT('FQ loop type not recognised!',-1)
      END SELECT

   CONTAINS

      SUBROUTINE square_loop(scheme,phase)
      use mm_T_buffers 
      IMPLICIT NONE
      INTEGER,      INTENT(IN) :: phase
      TYPE(scheme_paras), INTENT(IN) :: scheme

         weight = 1
!         CALL mm_open_T_buffer(scheme,phase)
!$OMP    BARRIER
!$OMP    DO SCHEDULE(DYNAMIC,1)
         DO j = 1, dims%RHS 
            id%RHS = j
            DO i = 1, dims%LHS
               id%LHS = i
               CALL mm_test_and_buffer_T_pair(LHS,RHS,id,weight)
            END DO
         END DO
!$OMP    END DO
!         CALL mm_close_T_buffer(scheme,phase)
      END SUBROUTINE square_loop

      SUBROUTINE triangular_loop
!         weight = 2
!         DO j = i, dims%RHS 
!            id%RHS = j
!            DO i = start, dims%LHS, step
!               id%LHS = i
!               CALL mm_test_and_buffer_T_pair(LHS,RHS,id,weight)
!            END DO
!         END DO
      END SUBROUTINE triangular_loop

   END SUBROUTINE mm_exe_FQ_loops

!-------------------------------------------------------------------------------

END MODULE mm_FQ_T_loops

!===============================================================================

MODULE mm_T_pair_builder

   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_gen_T_pairs,          &
             mm_init_T_pair_builder,  &
             mm_close_T_pair_builder

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_T_pair_builder(scheme,phase,pair_type)

      USE mm_T_pair_mould,  ONLY: mm_init_T_pair_mould
      USE mm_T_pair_tests,  ONLY: mm_init_T_pair_tests
      USE mm_T_buffers,     ONLY: mm_open_T_buffer

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      ! note we pass "phase" in addition to scheme%algorthm because
      ! we need to distinguish the NN phase of e.g. a complete FMM run.
      INTEGER, INTENT(IN) :: phase
      INTEGER, INTENT(IN) :: pair_type

      CALL mm_init_T_pair_mould(scheme,pair_type)
      CALL mm_init_T_pair_tests(scheme,phase)
      CALL mm_open_T_buffer(scheme,phase)

   END SUBROUTINE mm_init_T_pair_builder

!-------------------------------------------------------------------------------

   SUBROUTINE mm_close_T_pair_builder(scheme,phase)

      USE mm_T_pair_mould,  ONLY: mm_close_T_pair_mould
      USE mm_T_pair_tests,  ONLY: mm_close_T_pair_tests
      USE mm_T_buffers,     ONLY: mm_close_T_buffer

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER,      INTENT(IN) :: phase
      INTEGER :: mode

      CALL mm_close_T_pair_mould
      CALL mm_close_T_pair_tests
      CALL mm_close_T_buffer(scheme,phase)

   END SUBROUTINE mm_close_T_pair_builder

!-------------------------------------------------------------------------------

   SUBROUTINE mm_gen_T_pairs(LHS,RHS,scheme,T_searcher,pair_type,phase,box_pretest_in)

      USE mm_box_to_raw_map_builder, ONLY: mm_get_box_compression_map,  &
                                           mm_free_box_compression_map
      USE mm_FQ_T_loops,             ONLY: mm_exe_FQ_loops,  &
                                           mm_init_FQ_loops
      USE mm_overlap_tree_generator, ONLY: mm_exe_overlap_tree,   &
                                           mm_init_overlap_tree
      USE mm_grid_searcher_mod,          ONLY: mm_exe_grid_search
 
      IMPLICIT NONE
      TYPE(gen_mm_paras),    INTENT(INOUT) :: LHS, RHS
      TYPE(T_searcher_type), INTENT(IN)    :: T_searcher
      INTEGER,         INTENT(IN)    :: pair_type, phase
      LOGICAL, OPTIONAL,     INTENT(IN)    :: box_pretest_in
      TYPE(scheme_paras), INTENT(IN) :: scheme

      INTEGER :: shape
      LOGICAL       :: do_box_pretest

      ! currently only the NN phase exploits box pretesting
      IF ( PRESENT(box_pretest_in) ) THEN
         do_box_pretest = box_pretest_in
      ELSE
         IF (phase == DO_NN) CALL LSQUIT('box pretest controls broken!',-1) 
         do_box_pretest = .FALSE.
      END IF

      shape = T_searcher%shape
      IF (do_box_pretest) THEN
         ! need to build explicit mapping between boxed and raw moments
         CALL mm_get_box_compression_map(LHS)
         CALL mm_get_box_compression_map(RHS)
      END IF

      SELECT CASE (T_searcher%algorithm)
      CASE (FQ_LOOPS)
         CALL mm_init_FQ_loops(LHS,RHS,pair_type,do_box_pretest)
         CALL mm_exe_FQ_loops(LHS,RHS,scheme,phase,shape)
      CASE (TREE_SEARCH)
         IF ( (phase == DO_NN) .AND. (.NOT. do_box_pretest) ) THEN
            CALL LSQUIT('must use box pre-testing with tree-search algorithm!',-1)
         END IF 
         ! currently only searches NN space boxes using raw_paras
         CALL mm_init_overlap_tree(phase)
         CALL mm_exe_overlap_tree(LHS,RHS,shape)
      CASE (GRID_SEARCH)
         IF ( phase /= DO_NN ) THEN
            CALL LSQUIT('grid search currently only implemented for NN phase',-1)
         END IF
         CALL mm_exe_grid_search(LHS,RHS,shape)
      CASE DEFAULT
         CALL LSQUIT ('cannot reconcile requested T-pair search algorithm',-1)
      END SELECT

      IF (do_box_pretest) THEN
         CALL mm_free_box_compression_map(LHS%box_map)
         CALL mm_free_box_compression_map(RHS%box_map)
      END IF

   END SUBROUTINE mm_gen_T_pairs

!-------------------------------------------------------------------------------

END MODULE mm_T_pair_builder

