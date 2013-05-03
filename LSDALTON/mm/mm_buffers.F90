MODULE mm_scale_T_buffer
   USE mm_stats_mod
   USE mm_global_paras_mod
   use mm_mem
   use mm_memory_manager_mod
   IMPLICIT NONE
   PRIVATE
   ! public procedures 
   PUBLIC :: mm_init_scale_T_buffer,    &
             mm_free_scale_T_buffer,    &
             mm_scale_T_buffer_add

   INTEGER, PARAMETER :: BUFFER_SIZE = 500000 
   ! module wide variables
   INTEGER,      SAVE :: ndim_max
   TYPE(T_pair_batch), SAVE :: T_pair_buffer
!$OMP THREADPRIVATE(T_pair_buffer,ndim_max)

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_scale_T_buffer
      IMPLICIT NONE

      NULLIFY (T_pair_buffer%items)
!      ALLOCATE (T_pair_buffer%items(BUFFER_SIZE))
      call mm_allocate(T_pair_buffer%items,BUFFER_SIZE)
      T_pair_buffer%ndim = 0
      stat_tpack_chunks = zero
      stat_tpack_unique = zero
      stat_tpack_total = zero

   END SUBROUTINE mm_init_scale_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE mm_free_scale_T_buffer(T_contractor)
      IMPLICIT NONE
      EXTERNAL T_contractor

      IF (.NOT.ASSOCIATED(T_pair_buffer%items))  &
         CALL LSQUIT('T_pair_buffer not alloc.',-1)
      IF ( T_pair_buffer%ndim /= 0 ) THEN
         CALL expunge_scale_buffer(T_contractor)
         T_pair_buffer%ndim = 0
      END IF
      call mm_deallocate(T_pair_buffer%items)
!      DEALLOCATE (T_pair_buffer%items)
      NULLIFY (T_pair_buffer%items)

   END SUBROUTINE mm_free_scale_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE mm_scale_T_buffer_add(T_contractor,T_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair
      EXTERNAL T_contractor
      REAL(REALK) :: ratio
      
      stat_tpack_total = stat_tpack_total + one
      T_pair_buffer%ndim = T_pair_buffer%ndim +1
      T_pair_buffer%items(T_pair_buffer%ndim) = T_pair

      ! normalise T-vectors held in buffer for sorting purposes
      ratio = SQRT(SUM(T_pair%r_ab*T_pair%r_ab))
      ! we sort wrt x-axis first, and we want opposite vectors to be together
      IF (T_pair%r_ab(1) < zero) ratio = -ratio
      T_pair_buffer%items(T_pair_buffer%ndim)%paras%ratio = ratio
      T_pair_buffer%items(T_pair_buffer%ndim)%r_ab = T_pair%r_ab/ratio

      IF ( T_pair_buffer%ndim == BUFFER_SIZE ) THEN
         ! sort the buffer and pass all T-pairs to contractor
         CALL expunge_scale_buffer(T_contractor)
      END IF

   END SUBROUTINE mm_scale_T_buffer_add

!-------------------------------------------------------------------------------

   SUBROUTINE expunge_scale_buffer(T_contractor)

      USE mm_sorting, ONLY: mm_quicksort_wrt_vector,     &
                            mm_quicksort_wrt_ratio

      IMPLICIT NONE
      EXTERNAL T_contractor

      TYPE(T_pair_batch) :: ptr, ptr2
      INTEGER :: i, lo, hi
      REAL(REALK)   :: q1,q2

      ptr%ndim = MIN(BUFFER_SIZE,T_pair_buffer%ndim)
      ptr%items => T_pair_buffer%items(1:ptr%ndim)

      ! recursively sort wrt T-vectors, starting with x-component
      CALL sort_wrt_axis(1,ptr%items)

      ! expunge in batches of the same T-vector direction
      lo = 1
      DO i = 2, ptr%ndim
         q1 = ptr%items(i)%r_ab(1)
         q2 = ptr%items(i-1)%r_ab(1)
         IF (q1 == q2) THEN
            q1 = ptr%items(i)%r_ab(2)
            q2 = ptr%items(i-1)%r_ab(2)
            IF (q1 == q2) THEN
               q1 = ptr%items(i)%r_ab(3)
               q2 = ptr%items(i-1)%r_ab(3)
               IF (q1 == q2) CYCLE 
            END IF
         END IF
         hi = i-1
         ptr2%ndim = hi-lo+1
         ptr2%items => ptr%items(lo:hi)
         stat_tpack_unique = stat_tpack_unique + one
         CALL T_contractor(ptr2)
         lo = i
      END DO

      ! finally do last batch
      hi = ptr%ndim
      ptr2%ndim = hi-lo+1
      ptr2%items => ptr%items(lo:hi)
       stat_tpack_unique = stat_tpack_unique + one
      CALL T_contractor(ptr2)

      T_pair_buffer%ndim = 0
      stat_tpack_chunks = stat_tpack_chunks + one

   CONTAINS

!-------------------------------------------------------------------------------

      RECURSIVE SUBROUTINE sort_wrt_axis(xyz,items)

         IMPLICIT NONE
         INTEGER,       INTENT(IN)    :: xyz
         TYPE(T_pair_single), INTENT(INOUT) :: items(:)

         INTEGER :: i, lo, hi
         REAL(REALK)   :: q1,q2

         IF (SIZE(items) == 1) RETURN

         ! sort only if needed
         q1 = items(1)%r_ab(xyz)
         DO i = 2, SIZE(items)
            q2 = items(i)%r_ab(xyz)
            IF ( q2 < q1 ) THEN
               CALL mm_quicksort_wrt_vector(items,xyz)
               EXIT
            END IF
            q1 = q2
         END DO

         ! sub-sort next T-vector component
         lo = 1
         DO i = 2, SIZE(items)
            q1 = items(i-1)%r_ab(xyz)
            q2 = items(i)%r_ab(xyz)
            IF ( q2 /= q1 ) THEN
               hi = i-1
               IF (xyz == 3) THEN
                  CALL mm_quicksort_wrt_ratio(items(lo:hi))
                  RETURN
               ELSE
                  CALL sort_wrt_axis(xyz+1,items(lo:hi))
               END IF
               lo = i
            END IF
         END DO

         ! do last batch
         hi = SIZE(items)
         IF (xyz == 3) THEN
            CALL mm_quicksort_wrt_ratio(items(lo:hi))
            RETURN
         ELSE
            CALL sort_wrt_axis(xyz+1,items(lo:hi))
         END IF

      END SUBROUTINE sort_wrt_axis

!-------------------------------------------------------------------------------

   END SUBROUTINE expunge_scale_buffer

!-------------------------------------------------------------------------------

END MODULE mm_scale_T_buffer


!===============================================================================

MODULE mm_multi_T_buffer

   USE mm_stats_mod
   USE mm_global_paras_mod
   use mm_mem
   use mm_memory_manager_mod
   IMPLICIT NONE
   PRIVATE
   ! public procedures 
   PUBLIC :: mm_init_multi_T_buffer,    &
             mm_free_multi_T_buffer,    &
             mm_multi_T_buffer_add

   INTEGER, PARAMETER :: BUFFER_SIZE = 1000
   ! module wide variables
   INTEGER,      SAVE :: ndim_max
   TYPE(T_pair_batch), SAVE :: T_pair_buffer
!$OMP THREADPRIVATE(T_pair_buffer,ndim_max)

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_multi_T_buffer(ndim_max_in)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ndim_max_in


!      ndim_max = ndim_max_in*2  ! we multiply by two for 'paired" algorithm
!$OMP CRITICAL
      ndim_max = ndim_max_in
!$OMP END CRITICAL
      IF (ndim_max < 1) CALL LSQUIT('invalid multiple T-matrix dimension!',-1)
      NULLIFY (T_pair_buffer%items)
!      ALLOCATE (T_pair_buffer%items(BUFFER_SIZE))
      call mm_ALLOCATE(T_pair_buffer%items,BUFFER_SIZE)
      T_pair_buffer%ndim = 0
!$OMP MASTER
      stat_tpack_chunks = zero
      stat_tpack_total = zero
!$OMP END MASTER

   END SUBROUTINE mm_init_multi_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE mm_free_multi_T_buffer(T_contractor)

      IMPLICIT NONE
      EXTERNAL T_contractor

      IF (.NOT.ASSOCIATED(T_pair_buffer%items))  &
         CALL LSQUIT('T_pair_buffer not alloc.',-1)
      IF ( T_pair_buffer%ndim /= 0 ) THEN
         CALL expunge_multi_buffer(T_contractor)
         T_pair_buffer%ndim = 0
      END IF
!      DEALLOCATE (T_pair_buffer%items)
      call mm_deallocate(T_pair_buffer%items)
      NULLIFY (T_pair_buffer%items)

   END SUBROUTINE mm_free_multi_T_buffer

!-------------------------------------------------------------------------------
!
!   SUBROUTINE mm_multi_T_buffer_add(T_contractor,T_pair)
!
!      IMPLICIT NONE
!      TYPE(T_pair_single), INTENT(IN) :: T_pair
!      EXTERNAL T_contractor
!      
!      INTEGER, SAVE :: iRHS_last = 0
!      INTEGER :: iRHS
!
!      iRHS = T_pair%paras%RHS_id
!
!!!!!      IF ( BTEST(T_pair_buffer%ndim+1,0) ) THEN
!         ! number of buffer entries is even; try to expunge
!         IF ( (T_pair_buffer%ndim == ndim_max)   &
!              .OR. ((iRHS /= iRHS_last) .AND. (iRHS_last /= 0)) ) THEN
!            ! expunge buffer and build all the T-matrices at once
!            CALL T_contractor(T_pair_buffer)
!            T_pair_buffer%ndim = 0
!         END IF
!         iRHS_last = T_pair%paras%RHS_id
!!!!!      END IF
!
!      T_pair_buffer%ndim = T_pair_buffer%ndim +1
!      T_pair_buffer%items(T_pair_buffer%ndim) = T_pair
!
!   END SUBROUTINE mm_multi_T_buffer_add
!
!-------------------------------------------------------------------------------

   SUBROUTINE mm_multi_T_buffer_add(T_contractor,T_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair
      EXTERNAL T_contractor
      
      IF ( T_pair_buffer%ndim == BUFFER_SIZE ) THEN
         ! expunge buffer and build all the T-matrices at once
         CALL expunge_multi_buffer(T_contractor)
      END IF

      stat_tpack_total = stat_tpack_total + one
      T_pair_buffer%ndim = T_pair_buffer%ndim +1
      T_pair_buffer%items(T_pair_buffer%ndim) = T_pair

   END SUBROUTINE mm_multi_T_buffer_add

!-------------------------------------------------------------------------------

   SUBROUTINE expunge_multi_buffer(T_contractor)

      USE mm_sorting, ONLY: mm_quicksort_wrt_RHS

      IMPLICIT NONE
      EXTERNAL T_contractor

      INTEGER :: i, lo,hi
      INTEGER :: iRHS, iRHS_next, item_max
      TYPE(T_pair_batch) :: ptr

      lo = 1  
      item_max = MIN((BUFFER_SIZE-1),(T_pair_buffer%ndim-1))

      ! sort only if needed
      iRHS = T_pair_buffer%items(1)%paras%RHS_id
      DO i = 2, item_max
         iRHS_next = T_pair_buffer%items(i)%paras%RHS_id
         IF ( iRHS_next < iRHS ) THEN
            CALL mm_quicksort_wrt_RHS(T_pair_buffer%items(1:item_max))
            EXIT
         END IF
         iRHS = iRHS_next
      END DO

      DO i = 1, item_max
         iRHS = T_pair_buffer%items(i)%paras%RHS_id
         iRHS_next = T_pair_buffer%items(i+1)%paras%RHS_id
         ptr%ndim = i-lo+1
         IF ((iRHS /= iRHS_next) .OR. (ptr%ndim == ndim_max)) THEN
            ptr%items => T_pair_buffer%items(lo:i)
            CALL T_contractor(ptr)
            lo = i+1  
         END IF
      END DO

      ptr%ndim = (item_max+1)-lo+1
      ptr%items => T_pair_buffer%items(lo:(item_max+1))
      CALL T_contractor(ptr)

      T_pair_buffer%ndim = 0
      stat_tpack_chunks = stat_tpack_chunks + one

   END SUBROUTINE expunge_multi_buffer

!-------------------------------------------------------------------------------

END MODULE mm_multi_T_buffer


!===============================================================================

MODULE mm_T_buffers
   USE mm_global_paras_mod
   USE mm_stats_mod
   use mm_mem
   use mm_memory_manager_mod
   IMPLICIT NONE
   PRIVATE
   ! public procedures
   PUBLIC :: mm_add_to_T_buffer,   &
             mm_open_T_buffer,     &
             mm_close_T_buffer

   ! diagnostic flag
   CHARACTER(4), SAVE :: T_buffer_stat

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_add_to_T_buffer(T_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair
      EXTERNAL mm_selected_t_buffer
      INTERFACE
         SUBROUTINE mm_selected_t_contractor
         END SUBROUTINE mm_selected_t_contractor
      END INTERFACE

      CALL mm_selected_t_buffer(mm_selected_t_contractor,T_pair)

   END SUBROUTINE mm_add_to_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE mm_open_T_buffer(scheme,phase) 

      USE mm_T_contractors,  ONLY: mm_lock_T_con
      USE mm_tree_T_buffer_mod,  ONLY: mm_tree_T_buffer_init,      &
                                   mm_tree_T_buffer_add
      USE mm_multi_T_buffer, ONLY: mm_init_multi_T_buffer,     &
                                   mm_multi_T_buffer_add
      USE mm_scale_T_buffer, ONLY: mm_init_scale_T_buffer,     &
                                   mm_scale_T_buffer_add

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER,      INTENT(IN) :: phase
      INTEGER :: mode
      EXTERNAL mm_store_t_buffer

!$OMP MASTER
      CALL mm_init_buffer_stats('T')
      IF (T_buffer_stat == 'OPEN') CALL LSQUIT ('cannot reopen T_buffer',-1)
!$OMP END MASTER

      IF ( phase == DO_FQ  .OR.  phase == DO_NN ) THEN
         mode = scheme%T_con%NN_T_buffer
      ELSE
         mode = scheme%T_con%T_buffer
      END IF

      SELECT CASE (mode)
         CASE (SKIP_T_BUFFER)
            ! all T-contractions will be skipped by this choice of buffer
            CALL mm_store_t_buffer(mm_skip_T_buffer)
         CASE (NULL_T_BUFFER)
            CALL mm_store_t_buffer(mm_null_T_buffer)
         CASE (TREE_T_BUFFER)
             ! use tree-based sorting/evaluating module
            CALL mm_store_t_buffer(mm_tree_T_buffer_add)
            CALL mm_tree_T_buffer_init(TREE_LENGTH,scheme%T_con%sort_para)
         CASE (SCALE_T_BUFFER)
             ! use tree-based sorting/evaluating module
            CALL mm_store_t_buffer(mm_scale_T_buffer_add)
            CALL mm_init_scale_T_buffer
         CASE (MULTI_T_BUFFER)
             ! use buffer to drive multiple T matrix simultaneous build
            CALL mm_store_t_buffer(mm_multi_T_buffer_add)
            CALL mm_init_multi_T_buffer(TMATM_DF)
         CASE DEFAULT
            CALL LSQUIT ('cannot reconcile list type in mm_open_T_buffer',-1)
      END SELECT

!$OMP MASTER
      T_buffer_stat = 'OPEN'
      mm_lock_T_con = .TRUE.
!$OMP END MASTER

   END SUBROUTINE mm_open_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE mm_close_T_buffer(scheme,phase) 

      USE mm_T_contractors,  ONLY: mm_lock_T_con
      USE mm_tree_T_buffer_mod,  ONLY: mm_tree_T_buffer_finish
      USE mm_multi_T_buffer, ONLY: mm_free_multi_T_buffer
      USE mm_scale_T_buffer, ONLY: mm_free_scale_T_buffer

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER,      INTENT(IN) :: phase
      INTEGER :: mode
      INTERFACE
         SUBROUTINE mm_selected_t_contractor
         END SUBROUTINE mm_selected_t_contractor
      END INTERFACE


!$OMP MASTER
      IF (T_buffer_stat /= 'OPEN') CALL LSQUIT ('T_buffer already closed!',-1)
!$OMP END MASTER

      IF ( phase == DO_FQ  .OR.  phase == DO_NN ) THEN
         mode = scheme%T_con%NN_T_buffer
      ELSE
         mode = scheme%T_con%T_buffer
      END IF

      SELECT CASE (mode)
         CASE (SKIP_T_BUFFER)
            ! do nothing
         CASE (NULL_T_BUFFER)
            ! do nothing
         CASE (TREE_T_BUFFER)
            CALL mm_tree_T_buffer_finish(mm_selected_t_contractor)
         CASE (MULTI_T_BUFFER)
            CALL mm_free_multi_T_buffer(mm_selected_t_contractor)
         CASE (SCALE_T_BUFFER)
            CALL mm_free_scale_T_buffer(mm_selected_t_contractor)
         CASE DEFAULT
            CALL LSQUIT ('cannot reconcile list type in mm_close_T_buffer',-1)
      END SELECT

!$OMP MASTER
      T_buffer_stat = 'FREE'
!$OMP END MASTER
      mm_lock_T_con = .FALSE.

   END SUBROUTINE mm_close_T_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE mm_null_T_buffer(T_contractor,T_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair
      EXTERNAL T_contractor

      stat_tpack_total = stat_tpack_total + one
      CALL T_contractor(T_pair)

   END SUBROUTINE mm_null_T_buffer

!-------------------------------------------------------------------------------
! for diagnostic use only

   SUBROUTINE mm_skip_T_buffer(T_contractor,T_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair
      EXTERNAL T_contractor

      stat_tvect_builds = stat_tvect_builds + 1
      RETURN

   END SUBROUTINE mm_skip_T_buffer

!-------------------------------------------------------------------------------

END MODULE mm_T_buffers

!===============================================================================

MODULE mm_W_buffers

   use mm_memory_manager_mod
   USE mm_global_paras_mod
   USE mm_stats_mod
   use mm_mem
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_add_to_W_buffer,   &
             mm_open_W_buffer,     &
             mm_close_W_buffer

   ! diagnostic flag
   CHARACTER(4), SAVE :: W_buffer_stat

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_add_to_W_buffer(W_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: W_pair
      EXTERNAL mm_selected_w_buffer
      INTERFACE
         SUBROUTINE mm_selected_w_contractor
         END SUBROUTINE mm_selected_w_contractor
      END INTERFACE

      CALL mm_selected_w_buffer(mm_selected_w_contractor,W_pair)

   END SUBROUTINE mm_add_to_W_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE mm_open_W_buffer(scheme) 

      USE mm_W_contractors, ONLY: mm_lock_W_con
      USE mm_tree_T_buffer_mod, ONLY: mm_tree_T_buffer_init,      &
                                  mm_tree_T_buffer_add

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      EXTERNAL mm_store_w_buffer

      CALL mm_init_buffer_stats('W')
      IF (W_buffer_stat == 'OPEN') CALL LSQUIT ('cannot reopen W_buffer',-1)

      SELECT CASE (scheme%W_con%W_buffer)
         CASE (SKIP_W_BUFFER)
            ! all W-contractions will be skipped by this choice of buffer
            CALL mm_store_w_buffer(mm_skip_W_buffer)
         CASE (NULL_W_BUFFER)
            CALL mm_store_w_buffer(mm_null_W_buffer)
         CASE (TREE_W_BUFFER)
             ! use tree-based sorting/evaluating module
            CALL mm_store_w_buffer(mm_tree_T_buffer_add)
            CALL mm_tree_T_buffer_init(TREE_LENGTH,scheme%W_con%sort_para)
         CASE DEFAULT
         CALL LSQUIT ('cannot reconcile list type in mm_open_W_buffer',-1)
      END SELECT

      W_buffer_stat = 'OPEN'
      mm_lock_W_con = .TRUE.

   END SUBROUTINE mm_open_W_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE mm_close_W_buffer(scheme) 

      USE mm_W_contractors, ONLY: mm_lock_W_con
      USE mm_tree_T_buffer_mod, ONLY: mm_tree_T_buffer_finish

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTERFACE
         SUBROUTINE mm_selected_w_contractor
         END SUBROUTINE mm_selected_w_contractor
      END INTERFACE


      IF (W_buffer_stat /= 'OPEN') CALL LSQUIT ('W_buffer already closed!',-1)
      SELECT CASE (scheme%W_con%W_buffer)
         CASE (SKIP_W_BUFFER)
            ! do nothing
         CASE (NULL_W_BUFFER)
            ! do nothing
         CASE (TREE_W_BUFFER)
            CALL mm_tree_T_buffer_finish(mm_selected_w_contractor)
         CASE DEFAULT
            CALL LSQUIT ('cannot reconcile list type in mm_close_W_buffer',-1)
      END SELECT
      W_buffer_stat = 'FREE'
      mm_lock_W_con = .FALSE.

   END SUBROUTINE mm_close_W_buffer

!-------------------------------------------------------------------------------

   SUBROUTINE mm_null_W_buffer(W_contractor,W_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: W_pair
      EXTERNAL W_contractor

      CALL W_contractor(W_pair)

   END SUBROUTINE mm_null_W_buffer

!-------------------------------------------------------------------------------
! for diagnostic use only

   SUBROUTINE mm_skip_W_buffer(W_contractor,W_pair)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: W_pair
      EXTERNAL W_contractor

      RETURN

   END SUBROUTINE mm_skip_W_buffer

!-------------------------------------------------------------------------------

END MODULE mm_W_buffers
