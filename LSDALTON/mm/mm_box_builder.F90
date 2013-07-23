MODULE mm_box_constructor

   USE mm_global_paras_mod
   USE mm_stats_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_init_box_builder,                          &
             mm_free_box_builder,                          &
             mm_get_box_paras_at_level,                    &
             mm_get_box_qlm_at_level

   ! Pointers to unboxed LHS parameters and RHS moments & parameters
   TYPE(raw_mm_data),  POINTER, SAVE :: RHS_raw_mms
   TYPE(raw_mm_paras), POINTER, SAVE :: LHS_raw_paras(:)
   ! Packed LHS & RHS parameters and RHS moments at all levels
   TYPE(box_mm_data), POINTER, SAVE :: mms_at_lev(:)
   ! Deepest level of boxes in the hierarchy 
   INTEGER, SAVE :: deepest_level

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_box_builder(LHS_paras,RHS_mms,scheme)

      USE mm_box_procs, ONLY: mm_deepest_level

      IMPLICIT NONE
      TYPE(raw_mm_paras), TARGET, INTENT(INOUT) :: LHS_paras(:)
      TYPE(raw_mm_data),  TARGET, INTENT(INOUT) :: RHS_mms
      TYPE(scheme_paras),         INTENT(IN)    :: scheme

      deepest_level = mm_deepest_level(scheme) 
      stat_deepest_level = deepest_level

      ! initialise RHS pointers to raw moments and parameters
      RHS_raw_mms => RHS_mms
      ! initialise LHS pointers to raw parameters
      LHS_raw_paras => LHS_paras(:) 
      ! allocate foundation of hierarchical moment data structure
      CALL preallocate_levels
      ! update raw_mm_paras with box and branch data
      CALL init_box_paras(LHS_paras,RHS_mms%paras,scheme) 

   END SUBROUTINE mm_init_box_builder

!-------------------------------------------------------------------------------

   SUBROUTINE init_box_paras(LHS,RHS,scheme)

      USE mm_sorting,   ONLY: mm_quicksort_wrt_branches
      USE mm_box_procs, ONLY: mm_box, mm_branch, mm_grain, mm_box_centre

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: LHS(:)
      TYPE(raw_mm_paras), INTENT(INOUT) :: RHS(:)
      TYPE(scheme_paras), INTENT(IN)    :: scheme

      INTEGER :: i, nbox_max
      REAL(REALK)   :: grain, grain_inv

      grain     = mm_grain(scheme,deepest_level)
      grain_inv = one/grain
      nbox_max  = 2**deepest_level
      
      stat_min_box_size = grain

      DO i = 1, SIZE(RHS) 
         RHS(i)%box = mm_box(RHS(i)%cntr,grain_inv,nbox_max)
         RHS(i)%bra = mm_branch(RHS(i)%ext,grain_inv,scheme)
         RHS(i)%box_cntr = mm_box_centre(RHS(i)%box,grain)
         RHS(i)%map_up = 0     ! not defined yet
      END DO
      CALL mm_quicksort_wrt_branches(RHS)

      DO i = 1, SIZE(LHS)
         LHS(i)%box = mm_box(LHS(i)%cntr,grain_inv,nbox_max)
         LHS(i)%bra = mm_branch(LHS(i)%ext,grain_inv,scheme)
         LHS(i)%box_cntr = mm_box_centre(LHS(i)%box,grain)
         LHS(i)%map_up = 0     ! not defined yet
      END DO
      CALL mm_quicksort_wrt_branches(LHS)

   END SUBROUTINE init_box_paras

!-------------------------------------------------------------------------------
!FIXME: should add whether we want RHS or LHS, or both, as an argument
! For now, automatically builds both

   SUBROUTINE build_paras_at_level(level,scheme)

      USE mm_box_packer, ONLY: mm_init_pkd_paras

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER,      INTENT(IN) :: level

      TYPE(raw_mm_paras), POINTER :: ptr(:)

      IF ((level < TOP_LEVEL) .OR. (level > deepest_level)) THEN
         CALL LSQUIT ('cannot iterate paras to this level!',-1)
      END IF
 
      ! build packed parameters at deepest level if not already done
      IF (.NOT.ASSOCIATED(mms_at_lev(deepest_level)%RHS_paras)) THEN
         ptr => RHS_raw_mms%paras(:)
         CALL mm_init_pkd_paras(deepest_level,scheme,ptr,   &
                                mms_at_lev(deepest_level)%RHS_paras)
      END IF
      IF (.NOT.ASSOCIATED(mms_at_lev(deepest_level)%LHS_paras)) THEN
!         IF (scheme%LHS_mm_range == scheme%RHS_mm_range) THEN
!print *, "*** copying RHS_paras ****"
!            mms_at_lev(deepest_level)%LHS_paras => &
!                                      mms_at_lev(deepest_level)%RHS_paras(:)
!         ELSE
            ptr => LHS_raw_paras(:)
            CALL mm_init_pkd_paras(deepest_level,scheme,ptr,  &
                                   mms_at_lev(deepest_level)%LHS_paras)
!         END IF
      END IF

      IF (level < deepest_level) THEN
         ! iterate paras through the box hierarchy to the reqd level
         CALL iterate_paras_to_level(level,scheme,'RHS')
         CALL iterate_paras_to_level(level,scheme,'LHS')
      END IF

   END SUBROUTINE build_paras_at_level

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE iterate_paras_to_level(l,scheme,side)

      USE mm_box_packer, ONLY: mm_shift_and_pack_paras

      IMPLICIT NONE
      INTEGER,      INTENT(IN) :: l         ! level in box hierarchy
      TYPE(scheme_paras), INTENT(IN) :: scheme
      CHARACTER(3),       INTENT(IN) :: side

      TYPE(box_mm_paras), POINTER :: ptr(:)
      INTEGER :: l_down
      
      l_down = l+1
      SELECT CASE (side)
      CASE ('RHS')
         IF (.NOT.ASSOCIATED(mms_at_lev(l_down)%RHS_paras)) THEN
            ! must build paras at deeper levels before higher levels
            CALL iterate_paras_to_level(l_down,scheme,'RHS')
         END IF
         ptr => mms_at_lev(l_down)%RHS_paras(:)
         ! build new paras from paras at previous level (including packing)
         CALL mm_shift_and_pack_paras(l,scheme,ptr,mms_at_lev(l)%RHS_paras)
      CASE ('LHS')
         IF (.NOT.ASSOCIATED(mms_at_lev(l_down)%LHS_paras)) THEN
            ! must build paras at deeper levels before higher levels
            CALL iterate_paras_to_level(l_down,scheme,'LHS')
         END IF
!         IF (scheme%LHS_mm_range == scheme%RHS_mm_range) THEN
!           mms_at_lev(l)%LHS_paras => mms_at_lev(l)%RHS_paras(:)
!         ELSE
           ptr => mms_at_lev(l_down)%LHS_paras(:)
           ! build new paras from paras at previous level (including packing)
           CALL mm_shift_and_pack_paras(l,scheme,ptr,mms_at_lev(l)%LHS_paras)
!         END IF
      CASE DEFAULT
         CALL LSQUIT ('must build LHS or RHS paras!',-1)
      END SELECT

   END SUBROUTINE iterate_paras_to_level

!-------------------------------------------------------------------------------

   SUBROUTINE build_qlm_at_level(level,scheme,memory)

      USE mm_W_pair_builder, ONLY: mm_translate_raw_moments

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER,      INTENT(IN) :: level
      CHARACTER(4),       INTENT(IN) :: memory

      TYPE(raw_mm_data), POINTER :: ptr1
      TYPE(box_mm_data), POINTER :: ptr2
      INTEGER :: mms_dim

      IF ((level < TOP_LEVEL) .OR. (level > deepest_level)) THEN
         CALL LSQUIT ('cannot iterate boxed moments to this level!',-1)
      END IF

      ! build boxed moments at the deepest level if not already done
      IF (.NOT.ASSOCIATED(mms_at_lev(deepest_level)%qlm_W)) THEN
         mms_dim = SIZE(mms_at_lev(deepest_level)%RHS_paras)
         CALL allocate_lm_at_level(deepest_level,mms_dim,scheme%trans_LMAX)
         IF (.NOT.ASSOCIATED(RHS_raw_mms)) &
                  CALL LSQUIT('mm_box_builder not correctly initialised!',-1)
         ptr1 => RHS_raw_mms
         ptr2 => mms_at_lev(deepest_level)
         CALL mm_translate_raw_moments(scheme,ptr1,ptr2)
      END IF

      IF (level < deepest_level) THEN
         ! iterate RHS moments through the box hierarchy to the reqd level
         CALL iterate_qlm_to_level(level,scheme,memory)
      END IF

   END SUBROUTINE build_qlm_at_level

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE iterate_qlm_to_level(level,scheme,memory)

      USE mm_W_pair_builder, ONLY: mm_translate_moments

      IMPLICIT NONE
      INTEGER,      INTENT(IN) :: level
      TYPE(scheme_paras), INTENT(IN) :: scheme
      CHARACTER(4),       INTENT(IN) :: memory

      INTEGER :: l, mms_dim, l_down
      TYPE(box_mm_data), POINTER :: ptr1, ptr2

      l_down = level+1
      ! note we always build qlm_W first, then scale to give qlm_T
      IF (.NOT.ASSOCIATED(mms_at_lev(l_down)%qlm_W)) THEN
         ! must have qlm at deeper levels before higher levels
         CALL iterate_qlm_to_level(l_down,scheme,memory)
      END IF

      ! must build boxed paras before boxed moments
      IF (.NOT.ASSOCIATED(mms_at_lev(l_down)%RHS_paras)) THEN
         CALL build_paras_at_level(l_down,scheme)
      END IF

      mms_dim = SIZE(mms_at_lev(level)%RHS_paras)
      CALL allocate_lm_at_level(level,mms_dim,scheme%trans_LMAX)
      ptr1 => mms_at_lev(l_down)
      ptr2 => mms_at_lev(level)
      CALL mm_translate_moments(scheme,ptr1,ptr2)

!FIXME: this doesn't seem to work!!!
!      ! deallocate "old" moments from unused levels
!      IF (memory == 'free') THEN
!         DO l = TOP_level, deepest_level 
!            IF (l == level) CYCLE
!            CALL deallocate_level(l,'LM_ONLY')
!         END DO
!      END IF

   END SUBROUTINE iterate_qlm_to_level

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_box_paras_at_level(l,scheme,box_paras,side)

      IMPLICIT NONE
      INTEGER,      INTENT(IN) :: l
      TYPE(scheme_paras), INTENT(IN) :: scheme
      TYPE(box_mm_paras), POINTER    :: box_paras(:)
      CHARACTER(3),       INTENT(IN) :: side 

      IF (.NOT.ASSOCIATED(mms_at_lev)) STOP 'mms_at_lev should be allocated!'

      SELECT CASE (side)
      CASE ('RHS')
         IF (.NOT.ASSOCIATED(mms_at_lev(l)%RHS_paras)) THEN
            ! RHS paras not available at this level yet. So make them!
            CALL build_paras_at_level(l,scheme)
         END IF
         box_paras => mms_at_lev(l)%RHS_paras(:)
      CASE ('LHS')
         IF (.NOT.ASSOCIATED(mms_at_lev(l)%LHS_paras)) THEN
            ! LHS paras not available at this level yet. So make them!
            CALL build_paras_at_level(l,scheme)
         END IF
         box_paras => mms_at_lev(l)%LHS_paras(:)
      CASE DEFAULT 
         CALL LSQUIT ('must select just LHS or RHS paras to use',-1)
      END SELECT

   END SUBROUTINE mm_get_box_paras_at_level

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_box_qlm_at_level(l,scheme,qlm_T,side,memory)

      IMPLICIT NONE
      INTEGER,      INTENT(IN)  :: l
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      REAL(REALK),        POINTER     :: qlm_T(:,:)
      CHARACTER(3),       INTENT(IN)  :: side 
      CHARACTER(4),       INTENT(IN)  :: memory

      IF (.NOT.ASSOCIATED(mms_at_lev)) STOP 'mms_at_lev should be allocated!'
      IF (.NOT.ASSOCIATED(mms_at_lev(l)%qlm_T)) THEN 
         ! qlm data not available at this level yet. So make them!
         CALL build_qlm_at_level(l,scheme,memory)
      END IF

      IF (side == 'LHS') CALL LSQUIT('currently no LHS boxed mms built!',-1)
      IF (side == 'RHS') THEN
         qlm_T => mms_at_lev(l)%qlm_T(:,:)
      ELSE
         CALL LSQUIT("must select LHS or RHS boxed moments!",-1)
      END IF

   END SUBROUTINE mm_get_box_qlm_at_level

!-------------------------------------------------------------------------------

   SUBROUTINE preallocate_levels

      USE mm_memory_manager_mod, ONLY: mm_allocate

      IMPLICIT NONE
      INTEGER :: i

      IF (deepest_level == 0) RETURN
      IF (ASSOCIATED(mms_at_lev)) STOP 'mms_at_lev should not be allocated!'
      IF (deepest_level < TOP_LEVEL) THEN
         CALL LSQUIT('error allocating levels in box hierarchy',-1)
      END IF
      ! We allow for possibility of all levels being used, but sub-variables
      ! are only allocated if required, so this is not a big deal
      CALL mm_allocate(MEM_BOX_QLM,mms_at_lev,deepest_level)
      DO i = LBOUND(mms_at_lev,1), UBOUND(mms_at_lev,1) 
         NULLIFY (mms_at_lev(i)%LHS_paras)
         NULLIFY (mms_at_lev(i)%RHS_paras)
         NULLIFY (mms_at_lev(i)%qlm_W)
         NULLIFY (mms_at_lev(i)%qlm_T)
      END DO

   END SUBROUTINE preallocate_levels

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_lm_at_level(l,mms_dim,LMAX)

      USE mm_memory_manager_mod, ONLY: mm_allocate

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l, mms_dim, LMAX

      INTEGER :: lm_dim
      LOGICAL :: fail

      lm_dim = (1+LMAX)**2
      IF (l > deepest_level) CALL LSQUIT('invalid level to allocate!',-1)
      IF (l < TOP_LEVEL) CALL LSQUIT('invalid level to allocate!',-1)

      ! must test if pointers already allocated (compiler may not notice)
      fail = .FALSE.
      IF (ASSOCIATED(mms_at_lev(l)%qlm_T)) fail = .TRUE.
      IF (ASSOCIATED(mms_at_lev(l)%qlm_W)) fail = .TRUE.
      IF (fail) CALL LSQUIT('box lm data already allocated!',-1)
      CALL mm_allocate(MEM_BOX_QLM,mms_at_lev(l)%qlm_T,lm_dim,mms_dim)
      CALL mm_allocate(MEM_BOX_QLM,mms_at_lev(l)%qlm_W,lm_dim,mms_dim)
      ! must zero as they are built additively
      mms_at_lev(l)%qlm_T = zero
      mms_at_lev(l)%qlm_W = zero

   END SUBROUTINE allocate_lm_at_level

!-------------------------------------------------------------------------------

   SUBROUTINE deallocate_level(l,free_terms)

      USE mm_memory_manager_mod, ONLY: mm_deallocate
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l
      CHARACTER(7),  INTENT(IN) :: free_terms

      IF (l > deepest_level) CALL LSQUIT('invalid level to deallocate!',-1)
      IF (l < TOP_LEVEL) CALL LSQUIT('invalid level to deallocate!',-1)

      SELECT CASE (free_terms)
      CASE ('PM_ONLY')
         IF (ASSOCIATED(mms_at_lev(l)%LHS_paras,mms_at_lev(l)%RHS_paras)) THEN
            ! LHS and RHS paras the same and only point the same space
            call mm_deallocate(mms_at_lev(l)%RHS_paras)
         ELSE
            IF (ASSOCIATED(mms_at_lev(l)%RHS_paras))THEN
               call mm_DEALLOCATE(mms_at_lev(l)%RHS_paras)
            ENDIF
            IF (ASSOCIATED(mms_at_lev(l)%LHS_paras))THEN
               call mm_DEALLOCATE(mms_at_lev(l)%LHS_paras)
            ENDIF
         END IF
         NULLIFY (mms_at_lev(l)%LHS_paras)
         NULLIFY (mms_at_lev(l)%RHS_paras)
      CASE ('LM_ONLY')
         IF (ASSOCIATED(mms_at_lev(l)%qlm_T))THEN
            call mm_DEALLOCATE(mms_at_lev(l)%qlm_T)
         ENDIF
         IF (ASSOCIATED(mms_at_lev(l)%qlm_W))THEN
            call mm_DEALLOCATE(mms_at_lev(l)%qlm_W)
         ENDIF
         NULLIFY (mms_at_lev(l)%qlm_T)
         NULLIFY (mms_at_lev(l)%qlm_W)
      CASE DEFAULT
         CALL LSQUIT('terms to be deallocated from mms_at_lev not recognised!',-1)
      END SELECT

   END SUBROUTINE deallocate_level

!-------------------------------------------------------------------------------

   SUBROUTINE mm_free_box_builder
     use mm_memory_manager_mod, only: mm_deallocate
      IMPLICIT NONE
      INTEGER :: l

      NULLIFY (RHS_raw_mms)
      IF (ASSOCIATED(mms_at_lev)) THEN
        DO l = LBOUND(mms_at_lev,1), UBOUND(mms_at_lev,1) 
          IF (ASSOCIATED(mms_at_lev(l)%LHS_paras) .AND. &
            & ASSOCIATED(mms_at_lev(l)%LHS_paras,mms_at_lev(l)%RHS_paras)) THEN
             ! LHS and RHS paras the same and only point the same space
             call mm_DEALLOCATE(mms_at_lev(l)%RHS_paras)
          ELSE
             IF (ASSOCIATED(mms_at_lev(l)%RHS_paras))THEN
                call mm_DEALLOCATE(mms_at_lev(l)%RHS_paras)
             ENDIF
             IF (ASSOCIATED(mms_at_lev(l)%LHS_paras))THEN
                call mm_DEALLOCATE(mms_at_lev(l)%LHS_paras)
             ENDIF
          END IF
          IF (ASSOCIATED(mms_at_lev(l)%qlm_W)) THEN
             call mm_DEALLOCATE(mms_at_lev(l)%qlm_W)
          ENDIF
          IF (ASSOCIATED(mms_at_lev(l)%qlm_T)) THEN
             call mm_DEALLOCATE(mms_at_lev(l)%qlm_T)
          ENDIF
          NULLIFY (mms_at_lev(l)%LHS_paras)
          NULLIFY (mms_at_lev(l)%RHS_paras)
          NULLIFY (mms_at_lev(l)%qlm_W)
          NULLIFY (mms_at_lev(l)%qlm_T)
        END DO
        call mm_deallocate(mms_at_lev)
      END IF
      deepest_level = 0

   END SUBROUTINE mm_free_box_builder

!-------------------------------------------------------------------------------

END MODULE mm_box_constructor

