MODULE mm_tree_T_buffer_mod
   ! Dalton code, copyright, revised approx. 2002.05.13
   ! Pawel Salek
   ! 
   ! mm_tree_T_buffer packs together interaction pairs having same R vector
   ! to reduce the cost of generating interactions matrices T - the matrix
   ! is constructed only once for entire set.
   ! This module needs a routine that will evaluate entire set corresponding
   ! to same R vector.
   !
   ! One limitation cast on the implementation is to use possibly
   ! large but limited memory resources. This requires an ability to evaluate
   ! and expunge occasionally the interaction sets as they are created.
   !
   ! The pack_inter_tree is used to find out unique translation vectors.
   ! This is done by sorting all translation vectors first into
   ! surfaces, then into lines within surfaces, and finally into
   ! unique points.
   !
   ! Mark Watson : Sep. 2002.
   ! Generalised to allow for translation pairs too.
   ! i.e. can store and evaluate on the fly W-pairs 
   ! Also included field for dynamic LMAX

   USE mm_global_paras_mod
   USE mm_stats_mod

   ! Public functions:
   PUBLIC :: mm_tree_T_buffer_init,     &
             mm_tree_T_buffer_finish,   &
             mm_tree_T_buffer_add
   
   ! Public types:
   PUBLIC :: PointNode

   PRIVATE
   INTEGER, SAVE :: pack_sort_order

   ! PointNode: private structure used 
   ! for on-the-fly packing of the interactions.
   !
   TYPE PointNode
      INTEGER            :: level !for debugging:1-plane, 2-line, etc
      REAL(REALK)              :: coord ! normalized, the vector scale==1
      TYPE(PointNode), POINTER :: left, right
      ! the one below used by non-leaves only
      TYPE(PointNode), POINTER :: this ! for level==1: this plane, 
      !     level==2: this line
      !     level==3: self (for consistency)
      ! the three below used by leaves only (i.e. ones with level==3)
      TYPE(T_paras), POINTER :: abl1l2(:)
      ! largest LMAX in each abl1l2 set stored here
      INTEGER          :: entries_used
      INTEGER          :: LHS_LMAX, RHS_LMAX, LMAX
     ! For translations must know if qlm or Vff translation
      CHARACTER(1)           :: N_or_T 
   END TYPE PointNode

   ! PRIVATE DATA
   ! pack_inter_tree is a pack-on-the-fly, evaluate-when-needed structure for
   ! efficient, memory-adapted, R-vector packed evaluation of the
   ! interaction pairs.
   TYPE(PointNode), ALLOCATABLE,TARGET,SAVE :: pack_inter_tree(:)
   INTEGER, SAVE                      :: pack_inter_tree_used
   INTEGER,SAVE                       :: pack_total_in_current_chunk

CONTAINS   
   !---------------------------------------------------------------------------
   ! mm_tree_T_buffer_init:
   ! initialize data structures and counters when starting packing the tree.
   ! interaction sorting depends on whether the contraction scheme is
   ! T-matrix driven or U-matrix driven.  Here we initialise a 
   ! module-wide logical to determine this for all schemes
   !
   SUBROUTINE mm_tree_T_buffer_init(t_size, sort_order)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: t_size
      INTEGER, INTENT(IN) :: sort_order

      ALLOCATE(pack_inter_tree(t_size))
      pack_inter_tree_used = 0
      stat_tpack_total     = 0
      stat_tpack_unique    = 0
      stat_tpack_chunks    = one
      pack_sort_order      = sort_order
      pack_total_in_current_chunk = 0
   END SUBROUTINE mm_tree_T_buffer_init

   !---------------------------------------------------------------------------

   SUBROUTINE mm_tree_T_buffer_finish(pack_ev)
      IMPLICIT NONE
      EXTERNAL                    pack_ev
      TYPE(PointNode), POINTER :: root

      IF(pack_inter_tree_used>0) THEN
         root => pack_inter_tree(1)
         CALL mm_tpack_process(root, pack_ev)
         pack_inter_tree_used = 0
         pack_total_in_current_chunk = 0
      END IF
      IF (ALLOCATED(pack_inter_tree)) DEALLOCATE(pack_inter_tree)
   END SUBROUTINE mm_tree_T_buffer_finish

   !---------------------------------------------------------------------------

   SUBROUTINE mm_tpack_process(root, pack_ev)
      IMPLICIT NONE
      TYPE(PointNode), POINTER :: root
      EXTERNAL pack_ev 

      INTEGER :: i

      CALL node_evaluator(root, 1, pack_ev)

      DO i=1, pack_inter_tree_used
         IF(ASSOCIATED(pack_inter_tree(i)%abl1l2))&
              & DEALLOCATE(pack_inter_tree(i)%abl1l2)
         NULLIFY(pack_inter_tree(i)%abl1l2) ! not strictly neccesary
      END DO
      pack_inter_tree_used = 0    
      pack_total_in_current_chunk = 0
      NULLIFY(root)
   END SUBROUTINE mm_tpack_process

   ! ===================================================================
   ! point_node_new:
   ! NOTE: we choose for conveninence that:
   ! a). non-leaves have abl1l2 not associated.
   ! b). leaves have leaf%this => leaf
   RECURSIVE SUBROUTINE point_node_new(node,level,r,ll,lr,lm,NT,T_pair)
      IMPLICIT NONE
      TYPE(PointNode), POINTER    :: node
      TYPE(T_paras),   INTENT(IN) :: T_pair
      REAL(REALK),     INTENT(IN) :: r(3)
      INTEGER,   INTENT(IN) :: level
      ! put these variables in a TYPE ??
      INTEGER,   INTENT(IN) :: ll,lr,lm
     ! For translations must know if qlm or Vff translation
      CHARACTER(1),    INTENT(IN) :: NT 

      pack_inter_tree_used = pack_inter_tree_used + 1
      node => pack_inter_tree(pack_inter_tree_used)
      node%level = level
      node%coord = r(level)
      NULLIFY(node%left); NULLIFY(node%right)
      IF(level<3) THEN
         NULLIFY(node%abl1l2)
         CALL point_node_new(node%this, level+1, r,ll,lr,lm,NT, T_pair)
      ELSE 
         stat_tpack_unique = stat_tpack_unique + 1
         node%this => node
         node%entries_used = 1
         node%LMAX     = lm
         node%LHS_LMAX = ll
         node%RHS_LMAX = lr
         node%N_or_T   = NT
         ALLOCATE(node%abl1l2(START_LEN))
         node%abl1l2(1) = T_pair
      END IF
   END SUBROUTINE point_node_new

   !---------------------------------------------------------------------------
   ! node_evaluator:
   ! walks over the tree and evaluates the interaction pairs.
   ! 
   RECURSIVE SUBROUTINE node_evaluator(node, level, pack_ev)
      IMPLICIT NONE
      TYPE(PointNode), POINTER  :: node
      INTEGER, INTENT(IN) :: level
      EXTERNAL                     pack_ev   ! pair evaluator

      REAL(REALK), SAVE :: r(3)

      IF(.NOT.ASSOCIATED(node)) RETURN
      IF(ASSOCIATED(node%left)) CALL node_evaluator(node%left, level,pack_ev)
      IF(ASSOCIATED(node%right))CALL node_evaluator(node%right,level,pack_ev)
      r(level) = node%coord
      IF(level<3) THEN
         CALL node_evaluator(node%this, level+1, pack_ev)
      ELSE
         ! evaluate current point using provided evaluator.
         ! We sort node%ratio node%abl1l2 wrt to abl1l2(:)%RHS_id first.
         IF(node%entries_used>1) THEN
            CALL momentsort(node%abl1l2, node%entries_used) 
         ENDIF
         CALL mm_interface_T_pair_out(pack_ev,r,node)
      END IF

   END SUBROUTINE node_evaluator

   !---------------------------------------------------------------------------

   SUBROUTINE momentsort(abl1l2, N)
      IMPLICIT NONE
      TYPE(T_paras), INTENT(INOUT) :: abl1l2(:)
      INTEGER, INTENT(IN)    :: N
      TYPE(T_paras) :: idxs
      INTEGER :: i

      DO i = N/2, 1, -1
         CALL downheap(abl1l2, i, N)
      END DO

      ! abl1l2[1..N] is a heap now

      DO i=N,1,-1
         idxs      = abl1l2(i)
         abl1l2(i) = abl1l2(1)
         abl1l2(1) = idxs
         CALL downheap(abl1l2,1,i-1) ! restore a[1..i-1] heap
      END DO
   END SUBROUTINE momentsort

   !---------------------------------------------------------------------------

   SUBROUTINE downheap(abl1l2, lo, hi)
      IMPLICIT NONE
      TYPE(T_paras), INTENT(INOUT) :: abl1l2(:)
      INTEGER, INTENT(IN)    :: lo, hi

      !  PRE: a[lo+1..hi] is a heap 
      ! POST:  a[lo..hi]  is a heap 
      TYPE(T_paras) :: idxs
      INTEGER :: idx, child

      idxs  = abl1l2(lo)

      idx = lo
      makeheap: DO WHILE(idx <= hi/2) ! while k has child(s)
         child = 2*idx
         ! pick larger child...
         IF(child<hi) THEN ! exists right branch
            IF (pack_sort_order==SORT_BY_RHS_MMS) THEN
               IF(abl1l2(child)%RHS_id<abl1l2(child+1)%RHS_id .OR. &
                    & ((abl1l2(child)%RHS_id==abl1l2(child+1)%RHS_id &
                    & .AND. abl1l2(child)%ratio<abl1l2(child+1)%ratio)))&
                    &          child = child+1
            ELSE
               IF(abl1l2(child)%ratio<abl1l2(child+1)%ratio .OR. &
                 (ABS(abl1l2(child)%ratio-abl1l2(child+1)%ratio)&
                      <DISTINCT_T_TOL .AND. &
                  abl1l2(child)%RHS_id<abl1l2(child+1)%RHS_id)) child = child+1
            END IF
         END IF
         IF (pack_sort_order==SORT_BY_RHS_MMS) THEN
            IF(idxs%RHS_id >= abl1l2(child)%RHS_id) EXIT makeheap
         ELSE
            IF(idxs%ratio >= abl1l2(child)%ratio) EXIT makeheap
         END IF
         abl1l2(idx) = abl1l2(child)
         idx = child
      END DO makeheap

      abl1l2(idx) = idxs
   END SUBROUTINE downheap

   !---------------------------------------------------------------------------
   ! mm_tree_T_buffer_add:
   ! creates packed list of interactions corresponding to the same R.
   ! packing is done as-you-go, without creation of the intermediate
   ! unpacked list. Tree structure is used for sorting.
   ! Appends given data to specified tree.
   ! direction is a normalized r_pq with positive x component.
   SUBROUTINE mm_tree_T_buffer_add(pack_ev,T_pair_in)
      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(IN) :: T_pair_in 
      EXTERNAL pack_ev ! packed pair evaluator

      ! Local variables.
      TYPE(PointNode), POINTER :: node
      TYPE(T_paras),   POINTER :: newblock(:)
      TYPE(T_paras)            :: new_T_pair
      INTEGER            :: level
      REAL(REALK)              :: r_pq(3)
      REAL(REALK)              :: direction(3), r_pq_scale
      ! orders of multipole expanions per node
      INTEGER            :: ll,lr,lm
      ! for translations must know if qlm or Vff translation per node
      CHARACTER(1)             :: NT 

      r_pq = T_pair_in%r_ab
      new_T_pair = T_pair_in%paras
      NT = T_pair_in%N_or_T 
      ll = T_pair_in%paras%LHS_LMAX
      lr = T_pair_in%paras%RHS_LMAX
      lm = T_pair_in%LMAX

      node => pack_inter_tree(1) ! root is the first element of pack_inter_tree
      IF(pack_inter_tree_used+3 > SIZE(pack_inter_tree).OR.&
         pack_total_in_current_chunk>MAX_AVG_PER_NODE*SIZE(pack_inter_tree))&
         THEN
         stat_tpack_chunks = stat_tpack_chunks + one
         CALL mm_tpack_process(node, pack_ev)
      END IF

      ! do the job
      stat_tpack_total = stat_tpack_total + 1
      pack_total_in_current_chunk = pack_total_in_current_chunk+1
      r_pq_scale = SQRT(SUM(r_pq*r_pq))
      IF(r_pq(1).LT.0.0) r_pq_scale = -r_pq_scale
      direction = r_pq/r_pq_scale
      ! update ratio (r_ab is now the normalised vector)
      new_T_pair%ratio = r_pq_scale

      IF(pack_inter_tree_used==0) THEN
         CALL point_node_new(node,1,direction,ll,lr,lm,NT,new_T_pair)
         RETURN
      END IF

      ! level = 1 - searching for plane, 2 - searching for line, 3 - for point
      DO level = 1, 3
         nodesearch: DO WHILE(ABS(direction(level)-node%coord)>DISTINCT_T_TOL)
            IF(direction(level)<node%coord) THEN
               IF(.NOT.ASSOCIATED(node%left)) THEN
                  CALL point_node_new(node%left,level,direction,&
                                      ll,lr,lm,NT,new_T_pair)
                  RETURN
               ELSE
                  node => node%left
               ENDIF
            ELSE
               IF(.NOT.ASSOCIATED(node%right)) THEN
                  CALL point_node_new(node%right,level,direction,&
                                      ll,lr,lm,NT,new_T_pair)
                  RETURN
               ELSE
                  node => node%right
               ENDIF
            END IF
         END DO nodesearch
         node => node%this ! node%this points to node for leaves
      END DO

      ! append to the current point. It is guaranteed to be the final point
      ! expand first data block if needed.
      IF(node%entries_used==SIZE(node%abl1l2)) THEN
         ALLOCATE(newblock(node%entries_used*2))
         !node%abl1l2(1:node%entries_used): index removed to avoid crash on IBM
         newblock(1:node%entries_used) = node%abl1l2
         DEALLOCATE(node%abl1l2)
         node%abl1l2 => newblock
      END IF
      node%entries_used = node%entries_used +1
      IF (node%N_or_T /= NT) CALL LSQUIT('inconsistent data in one buffer node!',-1)
      node%N_or_T = NT
      node%LHS_LMAX = MAX(node%LHS_LMAX, ll)
      node%RHS_LMAX = MAX(node%RHS_LMAX, lr)
      node%LMAX = MAX(node%LMAX, lm)
      node%abl1l2(node%entries_used) = new_T_pair 

   END SUBROUTINE mm_tree_T_buffer_add

   !---------------------------------------------------------------------------
   ! Here we transform the data of this packer (old) into the new format

   SUBROUTINE mm_interface_T_pair_out(pack_ev,r,node)

      IMPLICIT NONE
      TYPE(PointNode), INTENT(IN) :: node
      REAL(REALK),     INTENT(IN) :: r(3)
      EXTERNAL pack_ev

      TYPE(T_pair_list) :: T_pairs_out
      INTEGER :: i

      T_pairs_out%N_or_T = node%N_or_T
      T_pairs_out%r_ab = r  ! normalised vector
      T_pairs_out%paras => node%abl1l2(1:node%entries_used)

      T_pairs_out%LMAX = node%LMAX
      T_pairs_out%lm_max = (1+node%LMAX)**2
      T_pairs_out%LHS_LMAX = node%LHS_LMAX
      T_pairs_out%RHS_LMAX = node%RHS_LMAX

      CALL pack_ev(T_pairs_out)

   END SUBROUTINE mm_interface_T_pair_out

   !---------------------------------------------------------------------------

END MODULE mm_tree_T_buffer_mod
