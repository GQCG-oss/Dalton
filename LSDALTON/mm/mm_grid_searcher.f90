MODULE mm_grid_searcher_mod

! Mark A. Watson.  Dec. 2002.  Oslo.
! Module to rapidly generate interaction pairs avoiding quadratic loops.
! Algorithm pre-computes a look-up table (in 3D) of occupied boxes for each
! branch of boxes.
! Look-up can then be done directly over a given box space (e.g. NN space)
! for a given branch, with sub-quadratic cost.

!FIXME: for now only implemented for NN clean-up phase with ext-test
! i.e. we only get NN space
! not generalised for different test, dynamic LMAX etc... 

   USE mm_global_paras_mod
   USE mm_grid_types
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_exe_grid_search
   REAL(REALK) :: NN_boxes, empty_boxes, full_boxes
   REAL(REALK) :: NN_overs, cls_overs, ncl_overs

CONTAINS

!===============================================================================

   SUBROUTINE mm_exe_grid_search(LHS,RHS,shape)

      USE mm_grid_map_builder, ONLY: mm_get_grid, mm_free_grid

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      INTEGER,      INTENT(IN) :: shape
      TYPE(grid_type),    POINTER    :: grid_head 
      INTEGER :: i, box(3), bra
      TYPE(hi_lo)   :: xtbox

      ! generate look-up table of RHS occupied boxes with O(n) cost
      CALL mm_get_grid(RHS%raw_paras,grid_head,xtbox)

      ! generate all interaction pairs and pass into evaluation buffer
      CALL gen_pairs(LHS,RHS,shape,grid_head,xtbox)

      ! deallocate grid index (data still held in RHS%raw_paras)
      CALL mm_free_grid(grid_head,xtbox)

   END SUBROUTINE mm_exe_grid_search

!-------------------------------------------------------------------------------
! Routine to identify Nearest Neighbour boxes for two interacting branches
! around a reference box

   SUBROUTINE get_NN_space(box,bra1,bra2,xtbox,NN_space)

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: box(3), bra1, bra2
      TYPE(hi_lo),   INTENT(IN)  :: xtbox
      TYPE(hi_lo),   INTENT(OUT) :: NN_space
      INTEGER :: WSnn

      WSnn = (bra1+bra2)/2
      IF ( (((bra1+bra2)/2)*2) /= (bra1+bra2) ) &
         CALL LSQUIT('Error in generating NN_space',-1)

      NN_space%hi(:) = MIN( (box(:) + WSnn), xtbox%hi(:) )
      NN_space%lo(:) = MAX( (box(:) - WSnn), xtbox%lo(:) )

   END SUBROUTINE get_NN_space

!-------------------------------------------------------------------------------

   SUBROUTINE gen_pairs(LHS,RHS,shape,grid_head,xtbox)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      INTEGER,      INTENT(IN) :: shape
      INTEGER,      EXTERNAL   :: LMAX_proc
      TYPE(grid_type),    TARGET     :: grid_head
      TYPE(hi_lo),        INTENT(IN) :: xtbox
      TYPE(grid_type), POINTER :: grid 
      TYPE(hi_lo)   :: NN_space
      INTEGER :: p,i,j,k, box(3), bra, wt

      NN_boxes = zero
      empty_boxes = zero
      full_boxes = zero
      NN_overs = zero
      ncl_overs = zero
      cls_overs = zero

      wt = 1
      IF (shape == SHAPE_TRIANGULAR) wt = 2
      DO p = 1, SIZE(LHS%raw_paras)
         box = LHS%raw_paras(p)%box
         bra = LHS%raw_paras(p)%bra
         grid => grid_head
         DO ! loop over grids for all branches
            CALL get_NN_space(box,bra,grid%bra,xtbox,NN_space)
            ! loop over NN space and get box contents
            DO k = NN_space%lo(3), NN_space%hi(3) 
               DO j = NN_space%lo(2), NN_space%hi(2) 
                  DO i = NN_space%lo(1), NN_space%hi(1) 
                     NN_boxes = NN_boxes +one
                     ! if box is occupied, test against list of occupants
                     IF (grid%box(i,j,k)%occ == 0) THEN
                        empty_boxes = empty_boxes +one
                        CYCLE
                     END IF
                     full_boxes = full_boxes +1
                     CALL test_pairs(p,LHS,RHS,grid%box(i,j,k)%head,shape,wt)
                  END DO
               END DO
            END DO
            grid => grid%next
            IF (.NOT.ASSOCIATED(grid)) EXIT
         END DO
      END DO

      IF (mm_stats_printed) THEN
         WRITE(lupri,'(/,A,I9)') "   LHS moments       =", SIZE(LHS%raw_paras) 
         WRITE(lupri,'(A,E20.12)') "   NN_boxes tested   =", NN_boxes
         WRITE(lupri,'(A,E20.12)')   "   empty_boxes       =", empty_boxes
         WRITE(lupri,'(A,E20.12)')   "   full_boxes        =", full_boxes
         WRITE(lupri,'(A,E20.12)')   "   NN_overs tested   =", NN_overs
         WRITE(lupri,'(A,E20.12)')   "   ncl_overs         =", ncl_overs
         WRITE(lupri,'(A,E20.12)')   "   cls_overs         =", cls_overs
         WRITE(lupri,'(A,I3,/)') "   shape =", shape
      END IF

   END SUBROUTINE gen_pairs

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE test_pairs(LHS_id,LHS,RHS,RHS_id,shape,weight)

      IMPLICIT NONE
      INTEGER,      INTENT(IN) :: LHS_id, shape, weight
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(occ_list),     INTENT(IN) :: RHS_id

      TYPE(T_pair_single) :: T_pair
      TYPE(LHS_RHS_type)  :: id
      LOGICAL :: A,B
      EXTERNAL mm_test_and_buffer_T_pair

      NN_overs = NN_overs +one
      A = (shape /= SHAPE_TRIANGULAR)
!radovan: id not set!
      call lsquit('id not set in test_pairs',-1)
!      B = (id%LHS < id%RHS)
!      IF ( A .OR. B ) THEN
!         id%LHS = LHS_id
!         id%RHS = RHS_id%id
!         CALL mm_test_and_buffer_T_pair(LHS,RHS,id,weight)
!      END IF
!
!      ! if there are more moments in the box, test them against LHS too
!      IF (ASSOCIATED(RHS_id%next)) THEN
!         CALL test_pairs(LHS_id,LHS,RHS,RHS_id%next,shape,weight) 
!      END IF
         
   END SUBROUTINE test_pairs

!-------------------------------------------------------------------------------

END MODULE mm_grid_searcher_mod

