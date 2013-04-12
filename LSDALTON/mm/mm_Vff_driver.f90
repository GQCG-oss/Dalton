MODULE mm_Vff_driver

   USE memory_handling, ONLY: init_threadmemvar, collect_thread_memory, &
        & mem_TurnOffThread_Memory, mem_TurnONThread_Memory
   USE mm_stats_mod
   USE mm_global_paras_mod
   USE mm_memory_manager_mod, ONLY: mm_allocate,mm_deallocate
   USE mm_T_contractors,  ONLY: mm_set_T_con_ptrs, mm_select_T_con,    &
                                mm_init_T_contractors,                 &
                                mm_free_T_contractors
   USE mm_W_contractors,  ONLY: mm_init_W_contractors,     &
                                mm_free_W_contractors
   USE mm_T_pair_builder, ONLY: mm_init_T_pair_builder,    &
                                mm_close_T_pair_builder,   &
                                mm_gen_T_pairs
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_get_raw_Vff

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_verify_data_in(RHS_in)

      IMPLICIT NONE
      TYPE(raw_mm_data),  INTENT(IN) :: RHS_in
      LOGICAL :: fail

      fail = .FALSE.
      IF (.NOT.ASSOCIATED(RHS_in%paras)) fail = .TRUE.
      IF (.NOT.ASSOCIATED(RHS_in%qlm_T)) fail = .TRUE.
      IF (.NOT.ASSOCIATED(RHS_in%qlm_W)) fail = .TRUE.
      IF (SIZE(RHS_in%paras) /= SIZE(RHS_in%qlm_T,2)) fail = .TRUE.
      IF (SIZE(RHS_in%paras) /= SIZE(RHS_in%qlm_W,2)) fail = .TRUE.
      IF (fail) CALL LSQUIT('mms pointers sent incorrectly to mm_Vff_driver',-1)

   END SUBROUTINE mm_verify_data_in

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_raw_Vff(scheme,LHS_paras,RHS_mms,Vff)

      USE mm_box_constructor, ONLY: mm_init_box_builder, mm_free_box_builder 
      USE LSTIMING, ONLY:  LSTIMER

      IMPLICIT NONE
      TYPE(raw_mm_paras),  INTENT(INOUT) :: LHS_paras(:)
      TYPE(raw_mm_data),   INTENT(INOUT) :: RHS_mms
      TYPE(scheme_paras),  INTENT(IN)    :: scheme
      REAL(REALK), TARGET, INTENT(OUT)   :: Vff(:,:)
      LOGICAL :: inc_NN

      REAL(REALK) :: TS,TE
      CALL LSTIMER('START ',TS,TE,-1)
      inc_NN = scheme%inc_NN
      CALL mm_verify_data_in(RHS_mms)

      IF ( scheme%algorithm /= DO_FQ ) THEN
         CALL mm_init_box_builder(LHS_paras,RHS_mms,scheme)
      END IF

      SELECT CASE (scheme%algorithm)
      CASE (DO_FQ)
         IF (inc_NN) CALL mm_get_FQ_Vff(scheme,LHS_paras,RHS_mms,Vff)
      CASE (DO_NN)
         IF (inc_NN) CALL mm_get_NN_Vff(scheme,LHS_paras,RHS_mms,Vff)
      CASE (DO_BQ)
         IF (inc_NN) CALL mm_get_NN_Vff(scheme,LHS_paras,RHS_mms,Vff)
         CALL mm_get_BQ_Vff(scheme,LHS_paras,Vff)
      CASE (DO_NlogN)
         IF (inc_NN) CALL mm_get_NN_Vff(scheme,LHS_paras,RHS_mms,Vff)
         CALL mm_get_NlogN_Vff(scheme,LHS_paras,Vff)
      CASE (DO_FMM)
         IF (inc_NN) CALL mm_get_NN_Vff(scheme,LHS_paras,RHS_mms,Vff)
      CALL LSTIMER('mm-nnV',TS,TE,-1)
         CALL mm_get_FMM_Vff(scheme,LHS_paras,Vff)
      CALL LSTIMER('mm-FMV',TS,TE,-1)
      CASE DEFAULT
         CALL LSQUIT('invalid algorithm requested!',-1)
      END SELECT

      IF ( scheme%algorithm /= DO_FQ ) THEN
         CALL mm_free_box_builder
      END IF

   END SUBROUTINE mm_get_raw_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_FQ_Vff(scheme,LHS_paras,RHS_mms,Vff)

      IMPLICIT NONE
      TYPE(raw_mm_paras), TARGET, INTENT(IN)    :: LHS_paras(:)
      TYPE(raw_mm_data),          INTENT(IN)    :: RHS_mms
      TYPE(scheme_paras),         INTENT(IN)    :: scheme
      REAL(REALK),        TARGET, INTENT(INOUT) :: Vff(:,:)

      INTEGER, PARAMETER :: pair_type = LHS_raw_RHS_raw
      TYPE(gen_mm_paras) :: LHS, RHS
      REAL(REALK) :: T0, SECOND, TTOTFQ

      ! this is just for diagnostic statistics (i.e. redundant to main code)
      stat_NN_not_FF = .FALSE.

      T0 = SECOND()
      NULLIFY(LHS%raw_paras,LHS%box_paras,LHS%box_map)
      NULLIFY(RHS%raw_paras,RHS%box_paras,RHS%box_map)

      LHS%raw_paras => LHS_paras(:)
      RHS%raw_paras => RHS_mms%paras(:)

      ! select the T-contractor to be stored/called via C code;
      ! note we use the (nearest neighbour ID) NN_ID here because
      ! if we run higher order methods this is the phase that is
      ! done to low order, and hence is the natural control flag
      ! for the FQ part too
      CALL mm_select_T_con(scheme%T_con%NN_ID)
      ! set T_contractor pointers   
      CALL mm_set_T_con_ptrs(Vff,RHS_mms%qlm_T)
call mem_TurnONThread_Memory()
!$OMP PARALLEL
call init_threadmemvar()
      CALL mm_init_T_contractors(scheme,DO_FQ,scheme%raw_LMAX)

      ! generate full multipole potential 
      CALL mm_init_T_pair_builder(scheme,DO_FQ,pair_type)
      CALL mm_gen_T_pairs(LHS,RHS,scheme,scheme%T_searcher(DO_FQ),pair_type,DO_FQ)
      CALL mm_close_T_pair_builder(scheme,DO_FQ)
      CALL mm_free_T_contractors
call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()

      TTOTFQ = SECOND()-T0
!      CALL TIMTXT('>>> TIME USED in mm_get_FQ_Vff   ', TTOTFQ, LUPRI)

   END SUBROUTINE mm_get_FQ_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_NN_Vff(scheme,LHS_paras,RHS_mms,Vff)

      USE mm_box_procs,   ONLY: mm_deepest_level
      USE mm_box_constructor, ONLY: mm_get_box_paras_at_level

      IMPLICIT NONE
      TYPE(raw_mm_paras), TARGET, INTENT(IN)    :: LHS_paras(:)
      TYPE(raw_mm_data),          INTENT(IN)    :: RHS_mms
      TYPE(scheme_paras),         INTENT(IN)    :: scheme
      REAL(REALK),        TARGET, INTENT(INOUT) :: Vff(:,:)

      INTEGER, PARAMETER :: pair_type = LHS_raw_RHS_raw
      TYPE(gen_mm_paras) :: LHS, RHS
      INTEGER :: level
      REAL(REALK) :: T0, SECOND, TTOTNN
      ! this is just for diagnostic statistics (i.e. redundant to main code)
      stat_NN_not_FF = .TRUE.

      T0 = SECOND()
      NULLIFY(LHS%raw_paras,LHS%box_paras,LHS%box_map)
      NULLIFY(RHS%raw_paras,RHS%box_paras,RHS%box_map)

      ! set up raw parameter pointers
      LHS%raw_paras => LHS_paras(:)
      RHS%raw_paras => RHS_mms%paras(:)
      IF ( scheme%NN_box_pretesting ) THEN
         ! using NN-box pre-testing, so need packed boxed paras too
         level = mm_deepest_level(scheme)
         CALL mm_get_box_paras_at_level(level,scheme,LHS%box_paras,'LHS')
         CALL mm_get_box_paras_at_level(level,scheme,RHS%box_paras,'RHS')
      END IF

      ! select the T-contractor to be stored/called via C code
      CALL mm_select_T_con(scheme%T_con%NN_ID)
      ! set T_contractor pointers   
      CALL mm_set_T_con_ptrs(Vff,RHS_mms%qlm_T)
call mem_TurnONThread_Memory()
!$OMP PARALLEL
call init_threadmemvar()
      CALL mm_init_T_contractors(scheme,DO_NN,scheme%raw_LMAX)

      ! generate NN Vff
      CALL mm_init_T_pair_builder(scheme,DO_NN,pair_type)
      CALL mm_gen_T_pairs(LHS,RHS,scheme,scheme%T_searcher(DO_NN),            &
                          pair_type,DO_NN,scheme%NN_box_pretesting)
      CALL mm_close_T_pair_builder(scheme,DO_NN)
      CALL mm_free_T_contractors
call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()

      TTOTNN = SECOND()-T0
!      CALL TIMTXT('>>> TIME USED in mm_get_NN_Vff   ', TTOTNN, LUPRI)

   END SUBROUTINE mm_get_NN_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_BQ_Vff(scheme,LHS_raw_paras,Vff)

      USE mm_box_constructor,ONLY: mm_get_box_qlm_at_level,            &
                                   mm_get_box_paras_at_level
      USE mm_box_procs,      ONLY: mm_deepest_level
      USE mm_W_pair_builder, ONLY: mm_get_raw_Vff_from_boxed_Vff

      IMPLICIT NONE
      TYPE(raw_mm_paras), TARGET, INTENT(IN)    :: LHS_raw_paras(:)
      TYPE(scheme_paras),         INTENT(IN)    :: scheme
      REAL(REALK),        TARGET, INTENT(INOUT) :: Vff(:,:)

      TYPE(gen_mm_paras)       :: LHS, RHS
      REAL(REALK), POINTER     :: qlm_T(:,:)
      REAL(REALK), POINTER     :: boxed_Vff(:,:)
      INTEGER            :: level, LMAX, lm_dim
      INTEGER, PARAMETER :: pair_type = LHS_box_RHS_box
      REAL(REALK) :: T0, SECOND, TTOTBQ

      ! this is just for diagnostic statistics (i.e. redundant to main code)
      stat_NN_not_FF = .FALSE.

      T0 = SECOND()
      NULLIFY(qlm_T,boxed_Vff)
      NULLIFY(LHS%raw_paras,LHS%box_paras,LHS%box_map)
      NULLIFY(RHS%raw_paras,RHS%box_paras,RHS%box_map)

      CALL mm_init_W_contractors(scheme%trans_LMAX)

      level = mm_deepest_level(scheme)
      IF (level < TOP_LEVEL) RETURN
      LMAX = scheme%trans_LMAX
      lm_dim = (1+LMAX)**2

      ! set up LHS and RHS boxed parameters
      CALL mm_get_box_paras_at_level(level,scheme,LHS%box_paras,'LHS')
      CALL mm_get_box_paras_at_level(level,scheme,RHS%box_paras,'RHS')
      CALL mm_get_box_qlm_at_level(level,scheme,qlm_T,'RHS','free')

      ! allocate temporary boxed potential
      CALL mm_allocate(MEM_BOX_VFF,boxed_Vff,lm_dim,INT(SIZE(LHS%box_paras)))
      boxed_Vff = zero

      ! select the T-contractor to be stored/called via C code
      CALL mm_select_T_con(scheme%T_con%ID)
      ! set T_contractor pointers   
      CALL mm_set_T_con_ptrs(boxed_Vff,qlm_T)
call mem_TurnONThread_Memory()
!$OMP PARALLEL
call init_threadmemvar()
      CALL mm_init_T_contractors(scheme,DO_BQ,scheme%trans_LMAX)

      ! generate FF Vff
      CALL mm_init_T_pair_builder(scheme,DO_BQ,pair_type)
      CALL mm_gen_T_pairs(LHS,RHS,scheme,scheme%T_searcher(DO_BQ),pair_type,DO_BQ)
      CALL mm_close_T_pair_builder(scheme,DO_BQ)
      CALL mm_free_T_contractors
call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()

      ! translate boxed potential to raw LHS centres
      CALL mm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,boxed_Vff,Vff)
      call mm_deallocate(boxed_Vff)
      TTOTBQ = SECOND()-T0
!      CALL TIMTXT('>>> TIME USED in mm_get_BQ_Vff   ', TTOTBQ, LUPRI)

      CALL mm_free_W_contractors

   END SUBROUTINE mm_get_BQ_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_NlogN_Vff(scheme,LHS_raw_paras,Vff)

      USE mm_box_constructor,ONLY: mm_get_box_qlm_at_level,            &
                                   mm_get_box_paras_at_level
      USE mm_box_procs,      ONLY: mm_deepest_level
      USE mm_W_pair_builder, ONLY: mm_get_raw_Vff_from_boxed_Vff

      IMPLICIT NONE
      TYPE(raw_mm_paras), TARGET, INTENT(IN)    :: LHS_raw_paras(:)
      TYPE(scheme_paras),         INTENT(IN)    :: scheme
      REAL(REALK),        TARGET, INTENT(INOUT) :: Vff(:,:)

      TYPE(gen_mm_paras)       :: LHS, RHS
      REAL(REALK), POINTER     :: qlm_T(:,:)
      REAL(REALK), POINTER     :: boxed_Vff(:,:)
      INTEGER            :: LHS_level, RHS_level, deepest_level
      INTEGER            :: LMAX, lm_dim
      INTEGER, PARAMETER :: pair_type = LHS_box_RHS_box
      REAL(REALK) :: T0, SECOND, TTOT

      ! this is just for diagnostic statistics (i.e. redundant to main code)
      stat_NN_not_FF = .FALSE.

      T0 = SECOND()
      NULLIFY(qlm_T,boxed_Vff)
      NULLIFY(LHS%raw_paras,LHS%box_paras,LHS%box_map)
      NULLIFY(RHS%raw_paras,RHS%box_paras,RHS%box_map)

      CALL mm_init_W_contractors(scheme%trans_LMAX)

      deepest_level = mm_deepest_level(scheme)
      IF (deepest_level < TOP_LEVEL) RETURN
      LMAX = scheme%trans_LMAX
      lm_dim = (1+LMAX)**2

      ! select the T-contractor to be stored/called via C code
      CALL mm_select_T_con(scheme%T_con%ID)
      ! get LHS boxed MM parameters
      LHS_level = deepest_level 
      CALL mm_get_box_paras_at_level(LHS_level,scheme,LHS%box_paras,'LHS')
      ! allocate boxed potential based on number of LHS boxes
      CALL mm_allocate(MEM_BOX_VFF,boxed_Vff,lm_dim,INT(SIZE(LHS%box_paras)))
      boxed_Vff = zero 

      ! generate far-field potential at deepest box centres
      DO RHS_level = deepest_level, TOP_LEVEL, -1
         ! get packed RHS MM paras and translated moments
         CALL mm_get_box_paras_at_level(RHS_level,scheme,RHS%box_paras,'RHS')
         CALL mm_get_box_qlm_at_level(RHS_level,scheme,qlm_T,'RHS','free')
         ! set T_contractor pointers   
         CALL mm_set_T_con_ptrs(boxed_Vff,qlm_T)
call mem_TurnONThread_Memory()
!$OMP PARALLEL
call init_threadmemvar()
         CALL mm_init_T_contractors(scheme,DO_NlogN,scheme%trans_LMAX)
         ! generate LFF potential at prescribed centres
         CALL mm_init_T_pair_builder(scheme,DO_NlogN,pair_type)
         CALL mm_gen_T_pairs(LHS,RHS,scheme,scheme%T_searcher(DO_NlogN),  &
                             pair_type,DO_NlogN)
         CALL mm_close_T_pair_builder(scheme,DO_NlogN)
         CALL mm_free_T_contractors
call collect_thread_memory()
!$OMP END PARALLEL
 call mem_TurnOffThread_Memory()
     END DO

      ! translate boxed potential to raw LHS centres
      CALL mm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,boxed_Vff,Vff)
      CALL mm_deallocate(boxed_Vff)
      TTOT = SECOND()-T0
!     CALL TIMTXT('>>> TIME USED in __NlogN_Vff', TTOT, LUPRI)

      CALL mm_free_W_contractors

   END SUBROUTINE mm_get_NlogN_Vff

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_FMM_Vff(scheme,LHS_raw_paras,Vff)

      USE mm_box_constructor,ONLY: mm_get_box_qlm_at_level,            &
                                   mm_get_box_paras_at_level
      USE mm_box_procs,      ONLY: mm_deepest_level
      USE mm_W_pair_builder, ONLY: mm_translate_parents_Vff,           &
                                   mm_get_raw_Vff_from_boxed_Vff
      USE LSTIMING, ONLY:  LSTIMER

      IMPLICIT NONE
      TYPE(raw_mm_paras), TARGET, INTENT(IN)    :: LHS_raw_paras(:)
      TYPE(scheme_paras),         INTENT(IN)    :: scheme
      REAL(REALK),                INTENT(INOUT) :: Vff(:,:)

      TYPE(gen_mm_paras)          :: LHS, RHS
      TYPE(box_mm_paras), POINTER :: p_box_paras(:)   ! parents data 
      REAL(REALK), POINTER        :: qlm_T(:,:)
      REAL(REALK), POINTER        :: Vff_tmp1(:,:)
      REAL(REALK), POINTER        :: Vff_tmp2(:,:)
      REAL(REALK), POINTER        :: Vff_ptr(:,:)
      INTEGER               :: level, l_up, deepest_level
      INTEGER               :: LMAX, lm_dim
      INTEGER               :: use_Vff_tmp
      INTEGER, PARAMETER    :: pair_type = LHS_box_RHS_box
      REAL(REALK) :: T0, SECOND, TTOT

      REAL(REALK) :: TS,TE

      CALL LSTIMER('START ',TS,TE,-1)
      ! this is just for diagnostic statistics (i.e. redundant to main code)
      stat_NN_not_FF = .FALSE.

      T0 = SECOND()
      NULLIFY(qlm_T,Vff_tmp1,Vff_tmp2,Vff_ptr)
      NULLIFY(LHS%raw_paras,LHS%box_paras,LHS%box_map)
      NULLIFY(RHS%raw_paras,RHS%box_paras,RHS%box_map)

      CALL mm_init_W_contractors(scheme%trans_LMAX) 

      use_Vff_tmp = 0
      deepest_level = mm_deepest_level(scheme)
      IF (deepest_level < TOP_LEVEL) RETURN
      LMAX = scheme%trans_LMAX
      lm_dim = (LMAX+1)**2

      ! select the T-contractor to be stored/called via C code
      CALL mm_select_T_con(scheme%T_con%ID)  
      CALL LSTIMER('**init',TS,TE,-1)
      ! generate far-field potential using whole box hierarchy
      DO level = TOP_LEVEL, deepest_level
         l_up  = level-1
         ! initialise T-contractor and allocate T-matrix
         ! get LHS boxed MM parameters
         CALL mm_get_box_paras_at_level(level,scheme,LHS%box_paras,'LHS')  
!        ! get packed RHS MM paras and translated moments
         CALL mm_get_box_paras_at_level(level,scheme,RHS%box_paras,'RHS') 
         CALL mm_get_box_qlm_at_level(level,scheme,qlm_T,'RHS','keep')    
         ! allocate temporary boxed potentials
         IF (use_Vff_tmp /= 2) THEN
            CALL mm_allocate(MEM_VFF_TMP,Vff_tmp1,lm_dim,INT(SIZE(LHS%box_paras)))
            Vff_tmp1 = zero
            Vff_ptr => Vff_tmp1(:,:) 
         ELSE
            CALL mm_allocate(MEM_VFF_TMP,Vff_tmp2,lm_dim,INT(SIZE(LHS%box_paras)))
            Vff_tmp2 = zero
            Vff_ptr => Vff_tmp2(:,:) 
         END IF
         ! set T_contractor pointers
         CALL mm_set_T_con_ptrs(Vff_ptr,qlm_T)
call mem_TurnONThread_Memory()
!$OMP PARALLEL
 call init_threadmemvar()
        CALL mm_init_T_contractors(scheme,DO_FMM,scheme%trans_LMAX) 
         ! generate LFF contribution to Vff at this level
         CALL mm_init_T_pair_builder(scheme,DO_FMM,pair_type)
         CALL mm_gen_T_pairs(LHS,RHS,scheme,scheme%T_searcher(DO_FMM),  &
                             pair_type,DO_FMM)
         CALL mm_close_T_pair_builder(scheme,DO_FMM)
         CALL mm_free_T_contractors
call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()
         ! get RFF contribution from parent level
         SELECT CASE (use_Vff_tmp)
            CASE (0)
               use_Vff_tmp = 2 
            CASE (1)
               ! RFF contribution: parents' translated FF ( parent | child )
               CALL mm_get_box_paras_at_level(l_up,scheme,p_box_paras,'LHS')  
               CALL mm_translate_parents_Vff(level,scheme,Vff_tmp2,Vff_tmp1, &  
                                             LHS%box_paras)           
               ! deallocate parent's space
               call mm_deallocate(Vff_tmp2)
               NULLIFY(Vff_tmp2)
               use_Vff_tmp = 2 
            CASE (2)
               ! RFF contribution: parents' translated FF ( parent | child )
               CALL mm_get_box_paras_at_level(l_up,scheme,p_box_paras,'LHS')
               CALL mm_translate_parents_Vff(level,scheme,Vff_tmp1,Vff_tmp2, &
                                             LHS%box_paras)
               ! deallocate parent's space
               call mm_deallocate(Vff_tmp1)
               NULLIFY(Vff_tmp1)
               use_Vff_tmp = 1 
         END SELECT 
      END DO
      CALL LSTIMER('**lvls',TS,TE,-1)

      ! translate boxed potential to raw LHS centres for final contraction
      IF (use_Vff_tmp == 1) THEN  ! Vff_tmp2 holds final boxed potential
         CALL mm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,Vff_tmp2,Vff)
      ELSE  ! Vff_tmp1 holds final boxed potential
         CALL mm_get_raw_Vff_from_boxed_Vff(LHS_raw_paras,scheme,Vff_tmp1,Vff)
      END IF
      CALL LSTIMER('**tran',TS,TE,-1)
      IF (ASSOCIATED(Vff_tmp1)) call mm_deallocate(Vff_tmp1)
      IF (ASSOCIATED(Vff_tmp2)) call mm_deallocate(Vff_tmp2)

      TTOT = SECOND()-T0
!      CALL TIMTXT('>>> TIME USED in mm_get_FMM_Vff  ', TTOT, LUPRI)

      CALL mm_free_W_contractors

   END SUBROUTINE mm_get_FMM_Vff

!-------------------------------------------------------------------------------

END MODULE mm_Vff_driver

