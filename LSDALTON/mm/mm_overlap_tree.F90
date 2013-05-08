MODULE mm_overlap_tree_generator

! Mark A. Watson, Dec-Jan 2002-2003. Oslo/Cambridge.
!
! Original overview
!--------------------
! Module to generate all classical interaction pairs between 
! nearest neighbour (NN) boxes.  This is done using a binary
! space partitioning algorithm to search for all occupied
! NN boxes in sub-quadratic time.
! Each node in the deepest tree level will then contain a
! linked list of all overlaps found in mutual NN boxes.
! This O(N) number of overlaps can then be searched rapidly
! for classically interacting pairs to obtain the final NN classical
! (non-overlapping) contributions.

! Multipole branches must be looped over quadratically, and for each
! branch-pair a separate NN space is defined and overlap tree generated. 
! For a single MM branch system (all overlaps similar extents) the
! mean cost to find all interactions grows as O(NlogN).  For more
! than one branch, this is then on average an O((N/A)*log(N/A)) cost
! for all pairs of the A branches.

! Extension - Feb 2003 (Mark Watson, Cambridge)
!-----------------------------------------------
! Generalisation of local space definition to include either NN boxes
! only or LFF boxes + NN boxes (for use in FMM and NlogN hierarchy).
! i.e. there is now a WS_para ("well-separatedness" parameter) which
! defines the size of the local space.

   USE mm_global_paras_mod
   use mm_mem
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_exe_overlap_tree,   &
             mm_init_overlap_tree

!-------------------------------------------------------------------------------

   TYPE plane_type
      INTEGER :: n     ! normal vector to the plane (along x,y, or z)
      INTEGER :: d     ! box position along axis where plane is drawn
   END TYPE plane_type

   TYPE vol_type
      ! these are integers since we currently hard-wire for NN box search
      INTEGER :: lo(3), hi(3)
   END TYPE vol_type

   TYPE item_type
      INTEGER :: id       ! overlap ID to recall from main LHS/RHS arrays
      INTEGER :: box(3)   ! box in which moment is held
   END TYPE item_type

   TYPE node_type
      TYPE(node_type), POINTER :: left, right   ! sub-trees
      INTEGER,   POINTER :: items(:)      ! id to main RHS array
      LOGICAL          :: leaf      ! end of tree (no sub-trees)
      INTEGER    :: nitems    ! if a leaf, number of items in it
      TYPE(plane_type) :: plane     ! split plane data
      TYPE(vol_type)   :: vol       ! volume of node (i.e. min/max xyz)
   END TYPE node_type

   TYPE batches_by_bra
      TYPE(batches_by_bra), POINTER :: next
      TYPE(mm_range) :: id
   END TYPE batches_by_bra

   ! Max no. of items per leaf
   INTEGER, PARAMETER :: LEAF_MAX = 10 
   ! Named constants for this module
   INTEGER, PARAMETER :: LEFT_OF_PLANE = 1
   INTEGER, PARAMETER :: RIGHT_OF_PLANE = 2
   INTEGER, PARAMETER :: CROSS_PLANE = 3
   ! Max. number of allowed tree recursions
   INTEGER, PARAMETER :: MAX_DIV = 100000000
   ! pointer to box and branch mm parameters
   TYPE(box_mm_paras), POINTER, SAVE :: RHS_paras(:)
   ! diagnostic variable to see how effectively tests are done
   REAL(REALK), SAVE :: npairs, ntested, nleaves

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_overlap_tree(phase)

      USE mm_T_pair_tests, ONLY: mm_def_WS_NN, mm_def_WS_RFF

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: phase
      EXTERNAL mm_store_ws_para

      SELECT CASE(phase)
         CASE (DO_FQ)
            CALL LSQUIT ('inappropriate use of TREE_SEARCH for non-local tests',-1)
         CASE (DO_NN)
            CALL mm_store_ws_para(mm_def_WS_NN)
         CASE (DO_BQ)
            CALL LSQUIT ('inappropriate use of TREE_SEARCH for non-local tests',-1)
         CASE (DO_NlogN)
            CALL mm_store_ws_para(mm_def_WS_RFF)
         CASE (DO_FMM)
            CALL mm_store_ws_para(mm_def_WS_RFF)
      END SELECT

   END SUBROUTINE mm_init_overlap_tree

!-------------------------------------------------------------------------------
! 
   SUBROUTINE mm_exe_overlap_tree(LHS,RHS,shape)

      USE mm_T_pair_tests, ONLY: FF_overs, NN_overs, cls_overs, ncl_overs

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS 
      INTEGER,      INTENT(IN) :: shape

      TYPE(node_type), POINTER :: tree 
      TYPE(item_type) :: ref_item
      TYPE(batches_by_bra), POINTER :: LHS_by_bra, RHS_by_bra
      TYPE(batches_by_bra), POINTER :: LHS_batch, RHS_batch
      TYPE(LHS_RHS_type) :: id
      INTEGER :: WS_para, weight, NOFF
      INTEGER :: i, common_level
      REAL(REALK) :: rawpairs
      LOGICAL, pointer :: tested(:)
      EXTERNAL mm_get_ws_para

      ! diagnostic variables
      FF_overs = zero
      NN_overs = zero
      ncl_overs = zero
      cls_overs = zero
      ntested = zero
      nleaves = zero

      IF ( (.NOT.ASSOCIATED(LHS%box_paras))   .OR.   &
           (.NOT.ASSOCIATED(RHS%box_paras)) ) THEN
         CALL LSQUIT('must build box_paras for binary tree search!',-1)
      END IF

      SELECT CASE (shape)
      CASE (SHAPE_SQUARE)
         npairs = SIZE(RHS%box_paras)*SIZE(LHS%box_paras)
         IF (ASSOCIATED(RHS%raw_paras)) THEN
            rawpairs = SIZE(RHS%raw_paras)*SIZE(LHS%raw_paras)
         END IF
         weight = 1
      CASE (SHAPE_TRIANGULAR)
         IF (SIZE(RHS%box_paras) /= SIZE(LHS%box_paras)) STOP 'bad shape!'
         npairs = half*SIZE(RHS%box_paras)*(1+SIZE(LHS%box_paras))
         IF (ASSOCIATED(RHS%raw_paras)) THEN
            rawpairs = half*SIZE(RHS%raw_paras)*(1+SIZE(LHS%raw_paras))
         END IF
         weight = 2
      CASE DEFAULT
         CALL LSQUIT('loop shape not recognised!',-1)
      END SELECT

      ! get linked-list of ID ranges sorted by branches
      CALL get_paras_by_bra(LHS%box_paras,LHS_by_bra)
      CALL get_paras_by_bra(RHS%box_paras,RHS_by_bra)

      NULLIFY(RHS_paras)
      RHS_paras => RHS%box_paras(:)
      ! loop over all branch combinations
      NOFF = MYNUM
      LHS_batch => LHS_by_bra
      LHS_branches: DO  ! until pointer disassociated
         RHS_batch => RHS_by_bra
         RHS_branches: DO  ! until pointer disassociated

            ! set "well-separatedness" parameter for local space;
            ! note that for the NlogN scheme, LHS and RHS are
            ! at different levels, and mm_def_WS_RFF returns WS_para
            ! for the largest box size; hence things will only work if
            ! the RHS data corresponds to the larger boxes!
            id%LHS = LHS_batch%id%lo
            id%RHS = RHS_batch%id%lo
            CALL mm_get_ws_para(LHS,RHS,id,WS_para)

            ! build the tree with primary RHS data only (for LHS to compare to)
            NULLIFY(tree)
            CALL build_tree(RHS_batch%id,tree,WS_para) 
            ! print tree for diagnostics
            CALL print_tree(tree)
      
            ! walk over tree and generate overlapping pairs
            ! "tested" is to avoid over-testing when particles cross leaves
            call mem_alloc_fmm(tested,RHS_batch%id%hi)  ! note LBOUND is conservative!
            DO i = LHS_batch%id%lo+NOFF, LHS_batch%id%hi, NNODES
               ref_item%id  = i
               common_level = RHS%box_paras(id%RHS)%level
               ref_item%box = translated_box(LHS%box_paras(i),common_level)
               tested = .FALSE.
               CALL search_tree(tree,ref_item,WS_para,weight,LHS,RHS,tested)
            END DO
            NOFF = MOD(NOFF+LHS_batch%id%hi-LHS_batch%id%lo+1,INT(NNODES))
            call mem_dealloc_fmm(tested)
      
            ! Deallocate tree memory for this branch pair
            CALL free_tree(tree) 

            IF (.NOT.ASSOCIATED (RHS_batch%next)) EXIT RHS_branches
            RHS_batch => RHS_batch%next
         END DO RHS_branches
         IF (.NOT.ASSOCIATED (LHS_batch%next)) EXIT LHS_branches
         LHS_batch => LHS_batch%next
      END DO LHS_branches

      NULLIFY(RHS_paras)
      ! deallocate linked-list of pointers to branches
      CALL free_paras_by_bra(LHS_by_bra)
      CALL free_paras_by_bra(RHS_by_bra)

      IF(mm_stats_printed) THEN
         WRITE(lupri,'(/,A,E20.12)') "   number of pairs    =", npairs
         WRITE(lupri,'(A,E20.12)')   "   number of leaves   =", nleaves
         WRITE(lupri,'(A,E20.12)')   "   box pairs tested   =", ntested
         WRITE(lupri,'(A,F7.2)') "   percentage tested  =", ntested/npairs*100
         WRITE(lupri,'(A,E20.12)')   "   FF boxes skipped   =", FF_overs
         WRITE(lupri,'(A,E20.12)')   "   NN boxes tested    =", NN_overs
         WRITE(lupri,'(A,E20.12)')   "   LFF boxes tested   =", &
              ntested-NN_overs-FF_overs
         IF (ASSOCIATED(RHS%raw_paras)) THEN
            WRITE(lupri,'(A,E20.12)')"   NCL overs skipped  =", ncl_overs
            WRITE(lupri,'(A,E20.12)')"   CLS overs included =", cls_overs
            ntested = ncl_overs + cls_overs
            cls_overs = 100*cls_overs
            ntested = 100*ntested
            WRITE(lupri,'(A,F7.2)')"   percentage tested  =", ntested/rawpairs
            WRITE(lupri,'(A,F7.2)')"   CLS percentage     =",cls_overs/rawpairs
         END IF
         WRITE(lupri,'(A,I3,/)') "   weight =", weight
      END IF
   END SUBROUTINE mm_exe_overlap_tree

!-------------------------------------------------------------------------------
! function to express a box at a deeper level of the hierarchy in terms
! of a box higher in the hierarchy (i.e. it's parent's parent's parent's...)
! so that separations can be expressed in terms of a common grid

   FUNCTION translated_box(box_paras,common_level)

      USE mm_box_procs, ONLY: mm_parent_box

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(IN) :: box_paras
      INTEGER,      INTENT(IN) :: common_level
      INTEGER :: level, tmp_box(3), translated_box(3)

      IF (box_paras%level < common_level) THEN
         ! must have secondary data at the deeper level in search tree
         CALL LSQUIT ('can only translate boxes UP to common level!',-1)
      END IF
      level = box_paras%level
      tmp_box = box_paras%box
      DO WHILE (level > common_level) 
         tmp_box = mm_parent_box(tmp_box)
         level = level -1
      END DO
      translated_box = tmp_box

   END FUNCTION translated_box

!-------------------------------------------------------------------------------
! Subroutine to build a linked list of indices partitioning the original
! array into chunks - each corresponding to a given branch.
! Note this only makes sense if the original array was sorted wrt branches
! in the first place!  This code still works if it's not, but we lose
! the power of the whole approach.

   SUBROUTINE get_paras_by_bra(box_paras,paras_by_bra)

      IMPLICIT NONE
      TYPE(box_mm_paras),   INTENT(IN) :: box_paras(:)
      TYPE(batches_by_bra), POINTER    :: paras_by_bra, tail
      INTEGER :: i, hi

      hi = SIZE(box_paras)
      NULLIFY(paras_by_bra)
      ALLOCATE(paras_by_bra)
      tail => paras_by_bra
      NULLIFY(tail%next)
      tail%id%lo = 1
      ! build list of indices
      DO i = 2, SIZE(box_paras)
         IF (box_paras(i)%bra == box_paras(i-1)%bra) CYCLE
         ! include check that data _is_ sorted here
         IF (box_paras(i)%bra < box_paras(i-1)%bra) STOP 'should sort data!!'
         tail%id%hi = i-1
         ALLOCATE(tail%next)
         tail => tail%next
         NULLIFY(tail%next)
         tail%id%lo = i
         tail%id%hi = i
      END DO
      tail%id%hi = hi

   END SUBROUTINE get_paras_by_bra

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE free_paras_by_bra(node)

      IMPLICIT NONE
      TYPE(batches_by_bra), POINTER :: node

      IF (ASSOCIATED(node%next)) CALL free_paras_by_bra(node%next)           
      IF (ASSOCIATED(node))THEN
         DEALLOCATE(node) 
      ENDIF
      NULLIFY(node)

   END SUBROUTINE free_paras_by_bra

!-------------------------------------------------------------------------------

   SUBROUTINE build_tree(range,root,WS_para)
   
      IMPLICIT NONE
      TYPE(mm_range),  INTENT(IN) :: range
      TYPE(node_type), POINTER    :: root 
      INTEGER,   INTENT(IN) :: WS_para
      INTEGER :: i,j, nitems

      ! initialise tree with all items in root
      ALLOCATE(root)      
      CALL initialise_node(root)
      nitems = 1 + range%hi - range%lo 
      ALLOCATE(root%items(nitems))
      root%nitems = nitems
      DO i = 1, nitems
         j = i-1 + range%lo
         root%items(i) = j     ! local array referring to global array
      END DO
      ! recursively sub-divide space and generate sub-trees
      CALL iterate_tree(root,WS_para)

   END SUBROUTINE build_tree

!-------------------------------------------------------------------------------

   SUBROUTINE initialise_node(node)

      IMPLICIT NONE
      TYPE(node_type), INTENT(OUT) :: node 

      node%leaf    = .TRUE.       ! assume a leaf until proven otherwise
      node%nitems  = -1           ! not filled yet 
      node%plane%n = 0            ! non-existent plane
      node%plane%d = 0 
      node%vol%hi  = 0            ! not filled, so zero volume
      node%vol%lo  = 0
      NULLIFY(node%items)               ! no contents yet in leaf
      NULLIFY(node%left, node%right)    ! no sub-trees yet

   END SUBROUTINE initialise_node

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE free_tree(node)

      IMPLICIT NONE
      TYPE(node_type), POINTER :: node

      IF (ASSOCIATED(node%items))THEN
         DEALLOCATE(node%items)
      ENDIF
      NULLIFY(node%items)
      IF (ASSOCIATED(node%left)) THEN
         CALL free_tree(node%left)
         NULLIFY(node%left)
      END IF
      IF (ASSOCIATED(node%right)) THEN
         CALL free_tree(node%right)
         NULLIFY(node%right)
      END IF
      IF (ASSOCIATED(node))THEN
         DEALLOCATE(node)
      ENDIF
      NULLIFY(node)

   END SUBROUTINE free_tree

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE iterate_tree(node,WS_para)

      IMPLICIT NONE
      TYPE(node_type), INTENT(INOUT) :: node 
      INTEGER,   INTENT(IN)    :: WS_para
      INTEGER :: i, nitems

      TYPE(node_type), POINTER :: left, right  ! sub-trees  
      INTEGER :: divisions = 1
      LOGICAL :: done_tree

      IF (.NOT.subdivide_criteria(node)) RETURN
      divisions = divisions +1
      IF (divisions > MAX_DIV) CALL LSQUIT('binary search tree too big!',-1)
      ! This should not be a leaf so we subdivide further
      node%leaf  = .FALSE.
      ! Initialise sub-trees
      ALLOCATE(left,right)
      node%left  => left
      node%right => right
      CALL initialise_node(left)
      CALL initialise_node(right)
      ! Define split-plane and redistribute items in sub-trees
      CALL get_min_max_xyz(node%items,node%vol)
      CALL get_split_plane(node%vol,node%plane)
      CALL share_out_items(node%items,node%plane,WS_para,left,right,done_tree)
      
      IF (done_tree) THEN
         ! both sub-trees identical so branch finished afterall
         node%leaf = .TRUE.
         NULLIFY(node%left, node%right)
         DEALLOCATE(left,right)
         NULLIFY(left,right)
         RETURN
      ELSE
         DEALLOCATE(node%items)   ! all info now in sub-trees
         NULLIFY(node%items)
         ! expand sub-trees if different from parent
         IF (left%nitems /= node%nitems) CALL iterate_tree(left,WS_para) 
         IF (right%nitems /= node%nitems) CALL iterate_tree(right,WS_para) 
      END IF

   END SUBROUTINE iterate_tree

!-------------------------------------------------------------------------------
! Routine to allocate the children nodes and fill with items
! according to split-plane criteria 
! "done_tree" flags if the orginal node should in fact not be split afterall

   SUBROUTINE share_out_items(items,plane,WS_para,left,right,done_tree)

      IMPLICIT NONE
      INTEGER,    INTENT(IN)  :: items(:)
      TYPE(plane_type), INTENT(IN)  :: plane
      INTEGER,    INTENT(IN)  :: WS_para
      TYPE(node_type),  INTENT(OUT) :: left, right
      LOGICAL,          INTENT(OUT) :: done_tree

      INTEGER  :: sub_tree(SIZE(items),2)  ! temporary
      TYPE(vol_type) :: vol(2)                   ! tmp volumes (L and R)
      INTEGER  :: box(3), i, nleft, nright
      INTEGER, PARAMETER :: iLEFT=1, iRIGHT=2
      
      done_tree = .TRUE.
      nleft = 0
      nright = 0
      DO i = 1, SIZE(items)
         box = RHS_paras(items(i))%box
         SELECT CASE (side_of_plane(plane,box,WS_para))
         CASE (LEFT_OF_PLANE)
            done_tree = .FALSE.    ! flag branches as different
            CALL add_to_subtree(iLEFT,nleft,i)
         CASE (RIGHT_OF_PLANE)
            done_tree = .FALSE.    ! flag branches as different
            CALL add_to_subtree(iRIGHT,nright,i)
         CASE (CROSS_PLANE)
            CALL add_to_subtree(iLEFT,nleft,i)
            CALL add_to_subtree(iRIGHT,nright,i)
         END SELECT
      END DO

      IF (done_tree) THEN
         ! left and right sub-trees identical so quit recursion
         RETURN
      ELSE
         ! allocate and build actual sub-nodes
         ! (deallocated in iterate_tree on next pass)
         ALLOCATE(left%items(nleft), right%items(nright))
         left%nitems  = nleft
         right%nitems = nright
         left%items(:)  = sub_tree(1:nleft,iLEFT) 
         right%items(:) = sub_tree(1:nright,iRIGHT) 
      END IF

   CONTAINS

      SUBROUTINE add_to_subtree(side,nitems,ip)
         IMPLICIT NONE
         INTEGER,       INTENT(IN)    :: side
         INTEGER, INTENT(IN)    :: ip
         INTEGER, INTENT(INOUT) :: nitems
         nitems = nitems+1
         sub_tree(nitems,side) = items(ip) 
      END SUBROUTINE add_to_subtree

   END SUBROUTINE share_out_items

!-------------------------------------------------------------------------------
! routine to return the extremum co-ordinates of a group of items
! FIXME: could use only for the root and update on-the-fly for deeper nodes

   SUBROUTINE get_min_max_xyz(node_items,vol)

      IMPLICIT NONE
      INTEGER,  INTENT(IN)  :: node_items(:)
      TYPE(vol_type), INTENT(OUT) :: vol
      INTEGER :: i, id
    
      IF (.NOT.ASSOCIATED(RHS_paras)) STOP 'RHS_paras not initialised!'
      id = node_items(1)
      vol%hi(:) = RHS_paras(id)%box(:)
      vol%lo(:) = RHS_paras(id)%box(:)
      DO i = 1, SIZE(node_items)
         id = node_items(i)
         vol%hi(:) = MAX( vol%hi(:), RHS_paras(id)%box(:))
         vol%lo(:) = MIN( vol%lo(:), RHS_paras(id)%box(:))
      END DO

   END SUBROUTINE get_min_max_xyz

!-------------------------------------------------------------------------------
!  Determine if a given group of items should be subdivided or not
!  For now just a simple test of number in group

   FUNCTION subdivide_criteria(node)

      IMPLICIT NONE
      TYPE(node_type), INTENT(IN) :: node 
      LOGICAL :: subdivide_criteria

      subdivide_criteria = (node%nitems > LEAF_MAX)
!      subdivide_criteria = .FALSE.
!      IF (node%nitems > LEAF_MAX) subdivide_criteria = .TRUE.
!      IF (node%nitems <= LEAF_MAX) THEN
!         print *, "few items in node = leaf"
!         print *, "node%nitems", node%nitems
!         print *, "plane", node%plane%n, node%plane%d
!         print *, "leaf", node%leaf
!         print *, "vol", node%vol%hi, node%vol%lo
!      END IF

   END FUNCTION subdivide_criteria

!-------------------------------------------------------------------------------

   SUBROUTINE get_split_plane(vol,plane)

      IMPLICIT NONE
      TYPE(vol_type),   INTENT(IN)  :: vol
      TYPE(plane_type), INTENT(OUT) :: plane
      INTEGER :: n(1)

      ! return longest dimension (n=1 impies x) as array of size 1
      n = MAXLOC(vol%hi(:) - vol%lo(:)) 
      ! define plane
      plane%n = n(1)
      plane%d = (vol%hi(n(1)) + vol%lo(n(1)))/2

   END SUBROUTINE get_split_plane

!-------------------------------------------------------------------------------
! note the definition of splitting is such that plane%d = 4 implies
! the plane splits between boxes 4 and 5.

   FUNCTION side_of_plane(plane,box,WS_para)

      IMPLICIT NONE
      TYPE(plane_type), INTENT(IN) :: plane
      INTEGER,    INTENT(IN) :: box(3)
      INTEGER,    INTENT(IN) :: WS_para
      INTEGER :: side_of_plane
      INTEGER :: s

      s = box(plane%n) - plane%d
      IF ( (s - WS_para) > 0) THEN
         side_of_plane = RIGHT_OF_PLANE
      ELSE IF ( (s + WS_para) <= 0) THEN
         side_of_plane = LEFT_OF_PLANE
      ELSE
         side_of_plane = CROSS_PLANE
      END IF

   END FUNCTION side_of_plane

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE search_tree(node,ref_item,WS_para,wt,LHS,RHS,tested)

      IMPLICIT NONE
      TYPE(node_type),    INTENT(IN)    :: node 
      TYPE(item_type),    INTENT(IN)    :: ref_item 
      INTEGER,      INTENT(IN)    :: WS_para
      INTEGER,      INTENT(IN)    :: wt
      TYPE(gen_mm_paras), INTENT(IN)    :: LHS, RHS
      LOGICAL,            INTENT(INOUT) :: tested(:)

      IF (node%leaf) THEN
         ! check against all items in the leaf
         CALL search_leaf(ref_item,node%items,LHS,RHS,wt,tested)
         RETURN
      ELSE
         ! else walk down tree to next node
         SELECT CASE (side_of_plane(node%plane,ref_item%box,WS_para))
         CASE (LEFT_OF_PLANE)
            CALL search_tree(node%left,ref_item,WS_para,wt,LHS,RHS,tested)
         CASE (RIGHT_OF_PLANE)
            CALL search_tree(node%right,ref_item,WS_para,wt,LHS,RHS,tested)
         CASE DEFAULT
            CALL search_tree(node%left,ref_item,WS_para,wt,LHS,RHS,tested)
            CALL search_tree(node%right,ref_item,WS_para,wt,LHS,RHS,tested)
         END SELECT
      END IF  

   END SUBROUTINE search_tree

!-------------------------------------------------------------------------------

   SUBROUTINE search_leaf(ref_item,leaf_items,LHS,RHS,wt,tested)

      IMPLICIT NONE
      TYPE(item_type),    INTENT(IN)    :: ref_item
      INTEGER,      INTENT(IN)    :: leaf_items(:)
      TYPE(gen_mm_paras), INTENT(IN)    :: LHS, RHS
      INTEGER,      INTENT(IN)    :: wt
      LOGICAL,            INTENT(INOUT) :: tested(:)
      EXTERNAL mm_test_and_buffer_T_pair

      INTEGER :: i
      TYPE(LHS_RHS_type) :: id

      DO i = 1, SIZE(leaf_items)
         IF (tested(leaf_items(i))) CYCLE  ! already tested against
         tested(leaf_items(i)) = .TRUE.
         id%LHS = ref_item%id
         id%RHS = leaf_items(i)
         ! test to skip half if triangular loop
         IF (id%RHS >= id%LHS .OR. wt ==1) THEN
            ntested = ntested +1
            CALL mm_test_and_buffer_T_pair(LHS,RHS,id,wt)
         END IF
      END DO

   END SUBROUTINE search_leaf

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE print_tree(node)

      IMPLICIT NONE
      TYPE(node_type), INTENT(IN) :: node 
      INTEGER :: i

      IF (node%leaf) THEN
         nleaves = nleaves +1
         ! print all items in the leaf
!         write(lupri,'(/,A)') "leaf items:"
!         write(lupri,'("nitems",I6)') node%nitems
!         write(lupri,'("plane",2I5)') node%plane%n, node%plane%d
!         DO i = 1, SIZE(node%items)
!            write(lupri,'("id",I5,5X,"box",3I5)') node%items(i)
!         END DO
      ELSE
         ! else walk down tree to next node
         CALL print_tree(node%left) 
         CALL print_tree(node%right) 
      END IF

   END SUBROUTINE print_tree

!-------------------------------------------------------------------------------

END MODULE mm_overlap_tree_generator
