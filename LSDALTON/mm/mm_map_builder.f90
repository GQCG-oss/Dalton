MODULE mm_box_to_raw_map_builder

! Module to build a map between one set of boxed parameters and all the 
! the sets of raw parameters that contribute to that box and branch.
! i.e. for given boxed moment, a linked list is generated of all the
! raw (unboxed) moments contributing to it.  This is very similar to
! the grid mapping below, except only occupied boxes are described, and
! there is no spacial structure in the list.
! Note this mapping is not completely flexible in that it will break
! if the raw_paras array is subsequently rearranged or sorted.

   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_get_box_compression_map,      &
             mm_free_box_compression_map

CONTAINS

!-------------------------------------------------------------------------------
! Routine to collect together raw_paras contributing to the same
! boxed moment (given by "map_up").

   SUBROUTINE mm_get_box_compression_map(data)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(INOUT) :: data
      INTEGER :: i, id

!$OMP MASTER
      IF (ASSOCIATED(data%box_map)) CALL LSQUIT ('map should not be built now!',-1)
      IF (.NOT.ASSOCIATED(data%raw_paras)) CALL LSQUIT ('need raw paras to map!',-1)
      ALLOCATE(data%box_map(SIZE(data%box_paras)))

      DO i = 1, SIZE(data%box_paras) 
         NULLIFY(data%box_map(i)%head)
         data%box_map(i)%occ = 0
      END DO
 
      DO i = 1, SIZE(data%raw_paras) 
         ! use "map_up" to identify common box and branch elements
         id = data%raw_paras(i)%map_up
         IF (id == 0) CALL LSQUIT ('must build packed paras first!',-1)
         ! add new element to list of occupants
         IF (data%box_map(id)%occ == 0) THEN
            ! linked list for this box is empty so start one
            data%box_map(id)%occ = 1
            ALLOCATE(data%box_map(id)%head)
            data%box_map(id)%head%id = i
            NULLIFY(data%box_map(id)%head%next)   ! list is empty
         ELSE
            ! make new entry and pre-pend to existing list
            CALL add_item(data%box_map(id),i)
         END IF
      END DO
!$OMP END MASTER

!   call print_map(data%box_map,data%raw_paras)

   CONTAINS

      SUBROUTINE add_item(box_map,raw_id)
   
         IMPLICIT NONE
         TYPE(id_list), INTENT(INOUT) :: box_map
         INTEGER, INTENT(IN) :: raw_id
         TYPE(id_node), POINTER :: new_node
   
         box_map%occ = box_map%occ + 1
         ALLOCATE(new_node)
         new_node%id = raw_id
         IF (ASSOCIATED(box_map%head%next)) THEN
            ! more than one entry in list (including head)
            ! so point new_node to old second entry
            new_node%next => box_map%head%next
            ! point head to new_node
            NULLIFY(box_map%head%next)
            box_map%head%next => new_node
         ELSE
            ! only head so far; make new_node our second entry
            box_map%head%next => new_node
            NULLIFY(new_node%next)   ! end of list
         END IF
   
      END SUBROUTINE add_item

   END SUBROUTINE mm_get_box_compression_map

!-------------------------------------------------------------------------------

   SUBROUTINE mm_free_box_compression_map(box_map)

      IMPLICIT NONE
      TYPE(id_list), POINTER :: box_map(:)
      INTEGER :: i

!$OMP MASTER
      IF (.NOT.ASSOCIATED(box_map)) THEN
!         RETURN  ! map not in fact built
      ELSE
         DO i = 1, SIZE(box_map)
            IF (box_map(i)%occ == 0) STOP 'should not have 0 occupancy'
            CALL free_linked_list(box_map(i)%head)
         END DO
         DEALLOCATE(box_map)
         NULLIFY(box_map)
      ENDIF
!$OMP END MASTER

   END SUBROUTINE mm_free_box_compression_map

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE free_linked_list(occ_node)

      TYPE(id_node), POINTER :: occ_node

      ! traverse linked-list freeing from bottom up recursively
      IF (ASSOCIATED(occ_node%next)) THEN
         CALL free_linked_list(occ_node%next)
         IF (ASSOCIATED(occ_node)) THEN
            DEALLOCATE(occ_node)
         ENDIF
         NULLIFY(occ_node)
      END IF
      IF (ASSOCIATED(occ_node)) THEN
         DEALLOCATE(occ_node)
      ENDIF
      NULLIFY(occ_node)   

   END SUBROUTINE free_linked_list

!-------------------------------------------------------------------------------

   subroutine print_map(map,raw_paras)

      implicit none
      type(raw_mm_paras), intent(in) :: raw_paras(:) 
      type(id_list), intent(in) :: map(:)
      type(id_node), pointer :: ptr
      integer :: i

      write (lupri,*) "printing map:"
      do i = 1, size(map)
         write(lupri,*)
         ptr => map(i)%head
         do 
            write(lupri,'(2I5,5X,I5,3X,3I5)') i, ptr%id, &
                         raw_paras(ptr%id)%bra, raw_paras(ptr%id)%box
            if (.not.associated(ptr%next)) exit 
            ptr => ptr%next
         end do
      end do

   end subroutine print_map

!-------------------------------------------------------------------------------

END MODULE mm_box_to_raw_map_builder

!===============================================================================

MODULE mm_grid_types

   USE mm_global_paras_mod
   IMPLICIT NONE
   PUBLIC

   TYPE hi_lo
      INTEGER :: hi(3), lo(3)
   END TYPE  hi_lo

   TYPE occ_list
      INTEGER :: id   ! ID of box member (position in raw input array)
      TYPE(occ_list), POINTER :: next   ! pointer to next box member ID 
   END TYPE occ_list

   TYPE box_type
      INTEGER :: occ
      TYPE(occ_list), POINTER :: head
   END TYPE box_type

   TYPE grid_type
      INTEGER :: bra      ! branch (one grid for each branch)
      INTEGER :: level    ! level FIXME (should be constant for now)
      TYPE(box_type),  POINTER :: box(:,:,:)
      TYPE(grid_type), POINTER :: next
   END TYPE grid_type

END MODULE mm_grid_types

!===============================================================================
! Mark A. Watson.  Dec. 2002.  Oslo.
!
! Module to rapidly generate interaction pairs avoiding quadratic loops.
! Algorithm pre-computes a look-up table (in 3D) of occupied boxes for each
! branch of boxes.
! Look-up can then be done directly over a given box space (e.g. NN space)
! for a given branch, with sub-quadratic cost.
!
!===============================================================================

MODULE mm_grid_map_builder

   USE mm_global_paras_mod
   USE mm_grid_types
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mm_get_grid,     &
             mm_free_grid

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_get_grid(overs,grid_head,xtbox)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN) :: overs(:)    ! overlap distributions
      TYPE(grid_type),    POINTER    :: grid_head
      TYPE(hi_lo)   :: xtbox
      INTEGER :: i

      ! get extremum box indices to avoid over-allocating grid
      CALL get_xtbox(overs,xtbox)
      ! initialise first grid in linked list over branches
      ALLOCATE(grid_head)
!      CALL init_new_grid(overs(1)%bra,overs(1)%level,grid_head,xtbox)
      CALL init_new_grid(overs(1)%bra,0,grid_head,xtbox)
      ! add remaining overlap distributions to this grid, or make/add to next 
      DO i = 1, SIZE(overs)
         CALL add_to_grid(i,overs(i),grid_head,xtbox)
      END DO
!      CALL print_grid(grid_head,xtbox)

   END SUBROUTINE mm_get_grid

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE mm_free_grid(grid,xtbox)

      IMPLICIT NONE
      TYPE(grid_type), POINTER    :: grid
      TYPE(hi_lo),     INTENT(IN) :: xtbox

      CALL free_boxes(grid%box,xtbox)
      IF (ASSOCIATED(grid%next)) THEN
         CALL mm_free_grid(grid%next,xtbox)
         IF (ASSOCIATED(grid))THEN
            DEALLOCATE(grid)
         ENDIF
         NULLIFY(grid)
      END IF 
      ! deallocation complete
      IF (ASSOCIATED(grid))THEN
         DEALLOCATE(grid)
      ENDIF
      NULLIFY(grid)

   END SUBROUTINE mm_free_grid

!-------------------------------------------------------------------------------

   SUBROUTINE free_boxes(boxes,xtbox)

      IMPLICIT NONE
      TYPE(box_type), POINTER    :: boxes(:,:,:)
      TYPE(hi_lo),    INTENT(IN) :: xtbox
      INTEGER :: i,j,k

      DO k = xtbox%lo(3), xtbox%hi(3) 
         DO j = xtbox%lo(2), xtbox%hi(2) 
            DO i = xtbox%lo(1), xtbox%hi(1) 
               IF (boxes(i,j,k)%occ == 0) CYCLE
               CALL free_occ(boxes(i,j,k)%head)
            END DO
         END DO
      END DO
      IF (ASSOCIATED(boxes)) THEN
         DEALLOCATE(boxes)
      ENDIF
      NULLIFY(boxes)

   END SUBROUTINE free_boxes

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE free_occ(occ_node)

      TYPE(occ_list), POINTER :: occ_node

      ! traverse linked-list freeing as we go
      IF (ASSOCIATED(occ_node%next)) THEN
         CALL free_occ(occ_node%next)
         IF (ASSOCIATED(occ_node))THEN
            DEALLOCATE(occ_node)
         ENDIF
         NULLIFY(occ_node)
      END IF
      ! deallocation of single box in grid complete
      IF (ASSOCIATED(occ_node))THEN
         DEALLOCATE(occ_node)
      ENDIF
      NULLIFY(occ_node)   

   END SUBROUTINE free_occ

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE add_to_grid(indx,over,grid,xtbox)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN)    :: over
      INTEGER,      INTENT(IN)    :: indx
      TYPE(hi_lo),        INTENT(IN)    :: xtbox
      TYPE(grid_type),    INTENT(INOUT) :: grid
      TYPE(grid_type),    POINTER       :: new_grid

      IF (over%bra == grid%bra) THEN
!         IF (over%level /= grid%level) CALL LSQUIT('can only do one level!')
         ! grid located, so now add to correct box occupancy list
         CALL add_to_box(indx,over%box,grid)
      ELSE
         IF (.NOT.ASSOCIATED(grid%next)) THEN
            ! no grid for this branch, so make one
            ALLOCATE(new_grid)
!            CALL init_new_grid(over%bra,over%level,new_grid,xtbox)
            CALL init_new_grid(over%bra,0,new_grid,xtbox)
            grid%next => new_grid
         END IF
         ! try adding again to next grid in list
         CALL add_to_grid(indx,over,grid%next,xtbox)
      END IF

   END SUBROUTINE add_to_grid

!-------------------------------------------------------------------------------
! simple routine to delimit the boced space that needs to be allocated
! in the look-up table.  i,e, we find the extremum box indices
!FIXME: could make more elaborate by having a separate set of parameters
!       for each grid (ie. also test against branch below) but not easy
!       with current algorithm

   SUBROUTINE get_xtbox(overs,xtbox)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(IN)  :: overs(:) 
      TYPE(hi_lo),        INTENT(OUT) :: xtbox
      INTEGER :: i

      xtbox%lo(:) = overs(1)%box(:)
      xtbox%hi(:) = overs(1)%box(:)
      DO i = 2, SIZE(overs)
         xtbox%lo(:) = MIN(xtbox%lo(:), overs(i)%box(:))
         xtbox%hi(:) = MAX(xtbox%hi(:), overs(i)%box(:))
      END DO 
!      print *
!      print '("xtbox%lo ",3I6)', xtbox%lo(:)
!      print '("xtbox%hi ",3I6)', xtbox%hi(:)

   END SUBROUTINE get_xtbox

!-------------------------------------------------------------------------------

   SUBROUTINE init_new_grid(branch,level,grid,xtbox)

      IMPLICIT NONE
      TYPE(grid_type), INTENT(OUT) :: grid
      INTEGER,   INTENT(IN)  :: branch, level
      TYPE(hi_lo),     INTENT(IN)  :: xtbox
      INTEGER :: i,j,k

      grid%bra = branch
      grid%level = level
      ! only allocate that part of the grid beyond which all boxes are empty
      ALLOCATE(grid%box( xtbox%lo(1) : xtbox%hi(1),              &
                         xtbox%lo(2) : xtbox%hi(2),              &
                         xtbox%lo(3) : xtbox%hi(3) ))
      ! list of grids for different branches is empty!
      NULLIFY(grid%next)
      ! initialise empty boxes
      DO k = xtbox%lo(3), xtbox%hi(3) 
         DO j = xtbox%lo(2), xtbox%hi(2) 
            DO i = xtbox%lo(1), xtbox%hi(1) 
               ! list of box occupants is empty!
               grid%box(i,j,k)%occ = 0
               NULLIFY(grid%box(i,j,k)%head)
            END DO
         END DO
      END DO

   END SUBROUTINE init_new_grid

!-------------------------------------------------------------------------------

   SUBROUTINE add_to_box(indx,box,grid)

      IMPLICIT NONE
      TYPE(grid_type), INTENT(INOUT) :: grid
      INTEGER,   INTENT(IN)    :: indx, box(3) 
      TYPE(occ_list),  POINTER :: new_node 
      INTEGER :: i,j,k

      i = box(1); j = box(2); k = box(3)
      ! add new element to list of occupants
      IF (grid%box(i,j,k)%occ == 0) THEN
         ! linked list of occupants is empty so start one
         grid%box(i,j,k)%occ = 1
         ALLOCATE(grid%box(i,j,k)%head)
         grid%box(i,j,k)%head%id = indx
         NULLIFY(grid%box(i,j,k)%head%next)   ! list is empty
      ELSE
         ! make new entry and pre-pend to list
         grid%box(i,j,k)%occ = grid%box(i,j,k)%occ + 1
         ALLOCATE(new_node)
         new_node%id = indx
         IF (ASSOCIATED(grid%box(i,j,k)%head%next)) THEN
            ! more than one entry in list (including head)
            ! so point new_node to old second entry
            new_node%next => grid%box(i,j,k)%head%next
            ! point head to new_node
            NULLIFY(grid%box(i,j,k)%head%next)
            grid%box(i,j,k)%head%next => new_node
         ELSE
            ! only head so far; make new_node our second entry
            grid%box(i,j,k)%head%next => new_node
            NULLIFY(new_node%next)   ! end of list
         END IF
      END IF

   END SUBROUTINE add_to_box

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE print_grid(grid,xtbox)

      IMPLICIT NONE
      TYPE(grid_type), INTENT(IN) :: grid
      TYPE(hi_lo),     INTENT(IN) :: xtbox
      INTEGER :: i,j,k
      REAL(REALK) :: boxes, full_boxes, empty_boxes
      REAL(REALK) :: max_occ, tot_occ, n_per_box

      write(lupri,'(/,A)') "analysis of box occupancies:"
      write(lupri,'(A,I4)') "branch", grid%bra
!      print *, "level ", grid%level
!      print *, "box entries:"

      max_occ = 0 
      tot_occ = 0 
      boxes = 0
      full_boxes = 0
      empty_boxes = 0

      DO k = LBOUND(grid%box,3), UBOUND(grid%box,3)
         DO j = LBOUND(grid%box,2), UBOUND(grid%box,2)
            DO i = LBOUND(grid%box,1), UBOUND(grid%box,1)
               boxes = boxes + 1
               tot_occ = tot_occ + grid%box(i,j,k)%occ
               max_occ = MAX(max_occ,REAL(grid%box(i,j,k)%occ,REALK))
               IF (grid%box(i,j,k)%occ == 0) THEN
                  empty_boxes = empty_boxes + 1
                  CYCLE
               ELSE
                  full_boxes = full_boxes + 1
!                  print '(" box:",3I5,7X,"occ:",I5,7X,"id:",I5)', i,j,k,   &
!                                           grid%box(i,j,k)%occ, &
!                                           grid%box(i,j,k)%head%id
!                  IF (ASSOCIATED(grid%box(i,j,k)%head%next)) & 
!                                    CALL print_box(grid%box(i,j,k)%head%next)
               END IF
            END DO
         END DO
      END DO

      write(lupri,'(A,E20.10)') "boxes           =", boxes
      write(lupri,'(A,E20.10)') "full boxes      =", full_boxes
      write(lupri,'(A,E20.10)') "empty boxes     =", empty_boxes
      write(lupri,'(A,E20.10)') "total moments   =", tot_occ
      write(lupri,'(A,E20.10)') "max occupancy   =", max_occ
      write(lupri,'(A,E20.10)') "mean occupancy  =", tot_occ/full_boxes
      write(lupri,'(A,E20.10)') "mean density    =", tot_occ/boxes


      IF (.NOT.ASSOCIATED(grid%next)) RETURN
      CALL print_grid(grid%next,xtbox)

   CONTAINS

      RECURSIVE SUBROUTINE print_box(box)
         IMPLICIT NONE
         TYPE(occ_list), INTENT(IN) :: box 
         print '(43X,"id:",I5)', box%id
         IF (.NOT.ASSOCIATED(box%next)) RETURN
         CALL print_box(box%next)
      END SUBROUTINE print_box

   END SUBROUTINE print_grid

!-------------------------------------------------------------------------------

END MODULE mm_grid_map_builder

