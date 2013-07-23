MODULE mm_T_pair_tests

   USE mm_global_paras_mod
   USE mm_stats_mod

   IMPLICIT NONE
   PRIVATE
   ! public procedures
   PUBLIC :: mm_init_T_pair_tests,         &
             mm_close_T_pair_tests,        &
             mm_def_WS_NN,                 &
             mm_def_WS_RFF
   PUBLIC :: FF_overs, NN_overs,           &
             cls_overs, ncl_overs
 
   REAL(REALK),   SAVE :: cls_radius
   INTEGER, SAVE :: shape       ! i.e. triangular or square loops
    ! flag to test initialisation
   CHARACTER(11), SAVE :: init_tests

   REAL(REALK), SAVE :: FF_overs, NN_overs
   REAL(REALK), SAVE :: cls_overs, ncl_overs

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_T_pair_tests(scheme,phase)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER,      INTENT(IN) :: phase
      EXTERNAL mm_store_t_tester
      EXTERNAL mm_store_test

      ! for diagnostic purposes only (i.e. redundant call)
      CALL mm_init_test_stats(phase)

      cls_radius = scheme%cls_radius
      shape = scheme%T_searcher(phase)%shape

      SELECT CASE (phase)
         CASE (DO_FQ)
            CALL mm_store_t_tester(mm_single_T_pair_test)
            IF (scheme%branch_free) THEN
               CALL mm_store_test(mm_test_r_ij)
            ELSE
               CALL mm_store_test(mm_test_ext)
            END IF
         CASE (DO_NN)
            IF (scheme%NN_box_pretesting) THEN
               ! use box_paras pretesting scheme to first find
               ! NN boxes, then raw_paras.
               CALL mm_store_t_tester(mm_multiple_T_pair_test)
               IF (scheme%branch_free) THEN
                  CALL mm_store_test(mm_test_r_ij)
               ELSE
                  CALL mm_store_test(mm_test_ext)
               END IF
            ELSE
               ! test over ALL raw parameters
               CALL mm_store_t_tester(mm_single_T_pair_test)
               IF (scheme%branch_free) THEN
                  CALL mm_store_test(mm_test_NN_r_ij)
               ELSE
                  CALL mm_store_test(mm_test_NN_ext)
               END IF
            END IF
         CASE (DO_BQ)
            CALL mm_store_t_tester(mm_single_T_pair_test)
            CALL mm_store_test(mm_test_FF)
         CASE (DO_NlogN)
            CALL mm_store_t_tester(mm_single_T_pair_test)
            CALL mm_store_test(mm_test_LFF)
         CASE (DO_FMM)
            CALL mm_store_t_tester(mm_single_T_pair_test)
            CALL mm_store_test(mm_test_LFF)
         CASE DEFAULT
            CALL LSQUIT ('unable to initialise T_pair_tests',-1)
      END SELECT 
!$OMP MASTER
      init_tests = 'initialised' 
!$OMP END MASTER

   END SUBROUTINE mm_init_T_pair_tests

!-------------------------------------------------------------------------------

   SUBROUTINE mm_close_T_pair_tests
      IMPLICIT NONE

!$OMP MASTER
      IF (init_tests /= 'initialised') STOP 'must initialise pair_tests!'
      init_tests = ' '
!$OMP END MASTER
 
   END SUBROUTINE mm_close_T_pair_tests

!-------------------------------------------------------------------------------
! routine to perform T-pair test on two raw or boxed moments, and if they
! pass the pair is sent for evaluation.

   SUBROUTINE mm_single_T_pair_test(LHS,RHS,id,weight)

      USE mm_T_buffers, ONLY: mm_add_to_T_buffer

      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      INTEGER,      INTENT(IN) :: weight
      TYPE(T_pair_single) :: T_pair
      LOGICAL, EXTERNAL   :: mm_included_pair
      EXTERNAL mm_stored_t_pair_mould


       IF (mm_included_pair(LHS,RHS,id)) THEN  
         ! pour all relevant info into the T-pair mould
          CALL mm_stored_t_pair_mould(LHS,RHS,id,weight,T_pair)
         ! pass single T_pair entity to contractors via buffer
          CALL mm_add_to_T_buffer(T_pair)
       END IF

   END SUBROUTINE mm_single_T_pair_test

!-------------------------------------------------------------------------------
! routine to perform T-pair tests on multiple raw moments, corresponding to
! the occupants of a single pair of boxes.  If each raw pair passes the test
! it is sent for evaluation.  Designed for efficient NN clean-up phase,
! therefore the "outer" box test is currently only "NN_boxes".

   SUBROUTINE mm_multiple_T_pair_test(LHS,RHS,box_id,weight)

      USE mm_T_buffers, ONLY: mm_add_to_T_buffer

      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: box_id
      INTEGER,      INTENT(IN) :: weight
      TYPE(id_node), POINTER :: LHS_box_map, RHS_box_map
      TYPE(LHS_RHS_type) :: map_id, raw_id
      LOGICAL :: included
      TYPE(T_pair_single) :: T_pair
      LOGICAL, EXTERNAL   :: mm_included_pair
      EXTERNAL mm_stored_t_pair_mould

      stat_multi_tests = stat_multi_tests + 1

      NULLIFY(LHS_box_map, RHS_box_map)
      included = NN_boxes(LHS%box_paras(box_id%LHS),RHS%box_paras(box_id%RHS))
      IF (.NOT.included) THEN
!         FF_overs = FF_overs + 1
         RETURN
      ELSE
!         NN_overs = NN_overs + 1
         map_id%LHS = LHS%box_paras(box_id%LHS)%id
         map_id%RHS = RHS%box_paras(box_id%RHS)%id

         ! loop over all raw moments in this pair of NN boxes

         RHS_box_map => RHS%box_map(map_id%RHS)%head
         RHS_raw_loop: DO ! until pointer disassociated
            raw_id%RHS = RHS_box_map%id

            LHS_box_map => LHS%box_map(map_id%LHS)%head
            LHS_raw_loop: DO ! until pointer disassociated
               raw_id%LHS = LHS_box_map%id

!       stat_single_tests = stat_single_tests + 1
       IF (mm_included_pair(LHS,RHS,raw_id)) THEN  
         ! pour all relevant info into the T-pair mould
          CALL mm_stored_t_pair_mould(LHS,RHS,raw_id,weight,T_pair)
         ! pass single T_pair entity to contractors via buffer
          CALL mm_add_to_T_buffer(T_pair)
       END IF
!                CALL mm_single_T_pair_test(LHS,RHS,raw_id,weight)

               ! only do next raw item in the LHS box if it exists 
               IF (.NOT.ASSOCIATED(LHS_box_map%next)) EXIT LHS_raw_loop
               LHS_box_map => LHS_box_map%next
            END DO LHS_raw_loop

            ! only do next raw item in the RHS box if it exists 
            IF (.NOT.ASSOCIATED(RHS_box_map%next)) EXIT RHS_raw_loop
            RHS_box_map => RHS_box_map%next
         END DO RHS_raw_loop

      END IF

   END SUBROUTINE mm_multiple_T_pair_test

!-------------------------------------------------------------------------------
! logical function to test if two gaussian overlaps are separated by more
! than a given cut-off radius.
! mm_test_r_ij is TRUE if overlaps are sufficiently well-separated.

   FUNCTION mm_test_r_ij(LHS,RHS,id)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL :: mm_test_r_ij
      REAL(REALK) :: r_ij(3), r_mod, r_zero

      mm_test_r_ij = .TRUE.
      r_ij   = RHS%raw_paras(id%RHS)%cntr - LHS%raw_paras(id%LHS)%cntr
      r_mod  = r_ij(1)*r_ij(1) + r_ij(2)*r_ij(2) + r_ij(3)*r_ij(3)
      r_zero = ZERO_DIST_TOL*ZERO_DIST_TOL

      IF ((r_mod - cls_radius*cls_radius) > r_zero) RETURN
      ! we must now also do the extent test because FCK3 doesn't
      ! include purely classical integrals even within the cut-off radius
      mm_test_r_ij = mm_test_ext(LHS,RHS,id)

   END FUNCTION mm_test_r_ij

!-------------------------------------------------------------------------------
! logical function to test if two gaussian overlaps have negligible 
! non-classical interaction energy based on separation of centres wrt extents.
! mm_test_ext is TRUE if interaction is "classical" in nature.

   FUNCTION mm_test_ext(LHS,RHS,id)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL :: mm_test_ext
      REAL(REALK) :: r_ij(3), ext_ij, r_mod, r_zero

      ext_ij = RHS%raw_paras(id%RHS)%ext + LHS%raw_paras(id%LHS)%ext
      r_ij   = RHS%raw_paras(id%RHS)%cntr - LHS%raw_paras(id%LHS)%cntr
      r_mod  = r_ij(1)*r_ij(1) + r_ij(2)*r_ij(2) + r_ij(3)*r_ij(3)
      ext_ij = ext_ij*ext_ij
      r_zero = ZERO_DIST_TOL*ZERO_DIST_TOL
      mm_test_ext = ((r_mod - ext_ij) > r_zero)

   END FUNCTION mm_test_ext

!-------------------------------------------------------------------------------
! interaction pair test routine for NN clean-up phase in branch-free scheme.
! an interaction should be treated by multipoles if in NN boxes and not too
! close.  i.e. if NN, test is purely against separation of centres and
! independent of extent of gaussian overlaps.

   FUNCTION mm_test_NN_r_ij(LHS,RHS,id)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL :: mm_test_NN_r_ij
      REAL(REALK)   :: r_ij(3), ext_ij, r_mod, r_zero
      INTEGER :: i,j, WS_NN
      INTEGER :: space

      i = id%LHS
      j = id%RHS
      WS_NN = (LHS%raw_paras(i)%bra + RHS%raw_paras(j)%bra)/2
      mm_test_NN_r_ij = .FALSE.

      space = ABS(LHS%raw_paras(i)%box(3) - RHS%raw_paras(j)%box(3))
      IF (space > WS_NN) THEN 
!         FF_overs = FF_overs + 1
         RETURN
      ELSE
         space = ABS(LHS%raw_paras(i)%box(2) - RHS%raw_paras(j)%box(2))
         IF (space > WS_NN) THEN 
!            FF_overs = FF_overs + 1
            RETURN
         ELSE
            space = ABS(LHS%raw_paras(i)%box(1) - RHS%raw_paras(j)%box(1))
            IF (space > WS_NN) THEN 
!               FF_overs = FF_overs + 1
               RETURN
            END IF
         END IF
      END IF

      ! pair are NN's, but is it a classical overlap?
!      NN_overs = NN_overs + 1
      r_ij   = RHS%raw_paras(j)%cntr - LHS%raw_paras(i)%cntr
      r_mod  = r_ij(1)*r_ij(1) + r_ij(2)*r_ij(2) + r_ij(3)*r_ij(3)
      r_zero = ZERO_DIST_TOL*ZERO_DIST_TOL

      IF ( (r_mod-cls_radius*cls_radius) < r_zero) THEN
         ! overlap is too close to be treated by MM
!         ncl_overs = ncl_overs +1
         RETURN
      END IF

!      cls_overs = cls_overs +1
      mm_test_NN_r_ij = .TRUE.

   END FUNCTION mm_test_NN_r_ij

!-------------------------------------------------------------------------------
! interaction pair test routine for NN clean-up phase in standard scheme.
! an interaction should be treated by multipoles if in NN boxes and without
! significant non-classical overlap. i.e. if NN, test is against separation of
! centres compared to extents of gaussian overlaps.

   FUNCTION mm_test_NN_ext(LHS,RHS,id)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL :: mm_test_NN_ext
      REAL(REALK)   :: r_ij(3), ext_ij, r_mod, r_zero
      INTEGER :: i,j, WS_NN
      INTEGER :: space

      i = id%LHS
      j = id%RHS
      WS_NN = (LHS%raw_paras(i)%bra + RHS%raw_paras(j)%bra)/2
      mm_test_NN_ext = .FALSE.

      space = ABS(LHS%raw_paras(i)%box(3) - RHS%raw_paras(j)%box(3))
      IF (space > WS_NN) THEN 
!         FF_overs = FF_overs + 1
         RETURN
      ELSE
         space = ABS(LHS%raw_paras(i)%box(2) - RHS%raw_paras(j)%box(2))
         IF (space > WS_NN) THEN 
!            FF_overs = FF_overs + 1
            RETURN
         ELSE
            space = ABS(LHS%raw_paras(i)%box(1) - RHS%raw_paras(j)%box(1))
            IF (space > WS_NN) THEN 
!               FF_overs = FF_overs + 1
               RETURN
            END IF
         END IF
      END IF

      ! pair are NN's, but is it a classical overlap?
!      NN_overs = NN_overs + 1
      ext_ij = RHS%raw_paras(j)%ext + LHS%raw_paras(i)%ext
      r_ij   = RHS%raw_paras(j)%cntr - LHS%raw_paras(i)%cntr
      r_mod  = DOT_PRODUCT(r_ij,r_ij)
      ext_ij = ext_ij*ext_ij
      r_zero = ZERO_DIST_TOL*ZERO_DIST_TOL

      IF ( (r_mod-ext_ij) < r_zero ) THEN
         ! overlap is not purely classical
!         ncl_overs = ncl_overs +1
         RETURN
      END IF

!      cls_overs = cls_overs +1
      mm_test_NN_ext = .TRUE.

   END FUNCTION mm_test_NN_ext

!-------------------------------------------------------------------------------
! logical function to test if two boxes are in each others "Far Field"
! based on separation.  Used by Boxed Quadratic (BQ) algorithm. 

   FUNCTION mm_test_FF(LHS,RHS,id)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL :: mm_test_FF

      mm_test_FF = .NOT.(NN_boxes(LHS%box_paras(id%LHS),RHS%box_paras(id%RHS))) 
      IF(mm_test_FF) THEN
!         FF_overs = FF_overs +1
      ELSE
!         NN_overs = NN_overs +1
      END IF

   END FUNCTION mm_test_FF

!-------------------------------------------------------------------------------
! logical function to test if two boxes are in each others "Local Far Field"
! based on separation;  used by FMM and NlogN algorithms;
! note that in the NlogN scheme boxes are at different levels, but the
! definition of LFF is such that we can then translate the box at the
! deeper level to the higher level and apply the same criteria as for
! two boxes at the deeper level (assuming branches join at higher levels) 
! FIXME: how about branch scheme where branches don't join at higher levels??

   FUNCTION mm_test_LFF(LHS_in,RHS_in,id)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN) :: LHS_in, RHS_in
      TYPE(LHS_RHS_type), INTENT(IN) :: id
      LOGICAL :: mm_test_LFF

      ! introduce these so we can translate to common grid if needed
      TYPE(box_mm_paras) :: LHS, RHS

      LHS = LHS_in%box_paras(id%LHS)
      RHS = RHS_in%box_paras(id%RHS)
      IF (LHS%level /= RHS%level) THEN
         ! All NN, LFF and RFF definitions are based on boxes in ONE grid.
         ! We must therefore translate the boxes into a common grid first.
         ! Of course, this can be avoided by passing the transformed paras.
         CALL translate_to_common_grid(LHS,RHS)
      END IF

      mm_test_LFF = .FALSE.
      IF (NN_boxes(LHS,RHS)) THEN
!         NN_overs = NN_overs +1
         RETURN 
      ELSE IF (RFF_boxes(LHS,RHS)) THEN
!         FF_overs = FF_overs +1
         RETURN
      END IF
      mm_test_LFF = .TRUE.

   END FUNCTION mm_test_LFF

!-------------------------------------------------------------------------------
! routine to define the separation parameter for NN interactions
! i.e. boxes closer than this separation are nearest neighbours
! only designed for cases where boxes are at the same level (same size)

   SUBROUTINE mm_def_WS_NN(LHS,RHS,id,WS_para)

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN)  :: LHS, RHS
      TYPE(LHS_RHS_type), INTENT(IN)  :: id
      INTEGER,      INTENT(OUT) :: WS_para

      WS_para = (LHS%box_paras(id%LHS)%bra + RHS%box_paras(id%RHS)%bra) /2

   END SUBROUTINE mm_def_WS_NN

!-------------------------------------------------------------------------------
! routine to define local space that contains both LFF and NN space;
! i.e. boxes further apart than this separation are RFF boxes;
! note that this routine is generalised for cases where the boxes are
! at different levels (of different sizes) as is needed for the NlogN
! contraction scheme; in this case, we must choose which box size to express
! WS_para in terms of; by default we take the box size at highest level

   SUBROUTINE mm_def_WS_RFF(LHS_in,RHS_in,id,WS_para)

      USE mm_box_procs, ONLY: mm_parent_bra

      IMPLICIT NONE
      TYPE(gen_mm_paras), INTENT(IN)  :: LHS_in, RHS_in
      TYPE(LHS_RHS_type), INTENT(IN)  :: id
      INTEGER,      INTENT(OUT) :: WS_para

      TYPE(box_mm_paras) :: LHS, RHS
      LHS = LHS_in%box_paras(id%LHS)
      RHS = RHS_in%box_paras(id%RHS)

      ! we define WS_para in terms of the highest common grid
      CALL translate_to_common_grid(LHS,RHS)

      ! LFF is related to the parent's NN space
      LHS%bra = mm_parent_bra(LHS%bra)
      RHS%bra = mm_parent_bra(RHS%bra)

      ! hence NN-space for _parents_ of common grid
      WS_para = (LHS%bra + RHS%bra) /2
      ! hence LFF-space in common grid itself is
      WS_para = 1 + 2*WS_para
      ! +1 because of asymmetry in LFF space, and we are conservative.

   END SUBROUTINE mm_def_WS_RFF

!-------------------------------------------------------------------------------
! Are 2 boxes within their Nearest Neighbour space (hence too close to be
!  treated by boxed multipole expansions exactly)?  
! If number of boxes between A and B < ((Br(A)+Br(B))/2) they are NN

   FUNCTION NN_boxes(LHS,RHS)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(IN) :: LHS, RHS
      LOGICAL :: NN_boxes
      INTEGER :: WS_NN
      INTEGER :: space

      ! test only valid if box grids are at the same level depth 
      IF (LHS%level /= RHS%level) CALL LSQUIT('levels not equal in NN_boxes',-1)

      WS_NN = (LHS%bra + RHS%bra)/2

      NN_boxes = .FALSE.
      space = ABS(LHS%box(3) - RHS%box(3))
      IF (space > WS_NN) THEN
         RETURN
      ELSE
         space = ABS(LHS%box(2) - RHS%box(2))
         IF (space > WS_NN) THEN 
            RETURN
         ELSE
            space = ABS(LHS%box(1) - RHS%box(1))
            IF (space > WS_NN) THEN 
               RETURN
            END IF
         END IF
      END IF
      ! pair are NN's
      NN_boxes = .TRUE.

   END FUNCTION NN_boxes

!-------------------------------------------------------------------------------
! logical function to test if two boxes are in each others "Remote Far Field"
! based on separation.  Used by FMM and NlogN algorithms. 

   FUNCTION RFF_boxes(LHS,RHS)

      USE mm_box_procs, ONLY: mm_parent_box, mm_parent_bra

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(IN) :: LHS, RHS
      LOGICAL :: RFF_boxes

      TYPE(box_mm_paras) :: LHS_up, RHS_up
      INTEGER :: WS_RFF, br1,br2, brp,brq
      INTEGER :: box1(3), box2(3), space(3) 

      ! test only valid if box grids are at the same level depth 
      IF (LHS%level /= RHS%level) CALL LSQUIT('levels not equal in RFF_boxes',-1)

      LHS_up = LHS
      RHS_up = RHS
      LHS_up%box = mm_parent_box(LHS%box)
      LHS_up%bra = mm_parent_bra(LHS%bra)
      RHS_up%box = mm_parent_box(RHS%box)
      RHS_up%bra = mm_parent_bra(RHS%bra)

      RFF_boxes = .NOT.(NN_boxes(LHS_up,RHS_up))

   END FUNCTION RFF_boxes

!-------------------------------------------------------------------------------

   FUNCTION same_box(LHS,RHS)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(IN) :: LHS, RHS
      LOGICAL :: same_box 

      IF (LHS%level /= RHS%level) CALL LSQUIT('levels not equal in same_box',-1)
      same_box = .FALSE.
      IF (LHS%box(1) == RHS%box(1)  .AND.   &
          LHS%box(2) == RHS%box(2)  .AND.   &
          LHS%box(3) == RHS%box(3))         &
      THEN
         same_box = .TRUE.
      END IF

   END FUNCTION same_box

!-------------------------------------------------------------------------------
! routine to generate a temporary set of boxes parameters at a common grid
! level.  parameters at the deeper level are translated to their parent's
! grid until the levels match.

   SUBROUTINE translate_to_common_grid(LHS,RHS)

      USE mm_box_procs, ONLY: mm_parent_box, mm_parent_bra

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: LHS, RHS

      IF (LHS%level == RHS%level) RETURN 
      IF (LHS%level > RHS%level) THEN
         DO WHILE (LHS%level > RHS%level)
            LHS%box = mm_parent_box(LHS%box)
            LHS%bra = mm_parent_bra(LHS%bra) 
            LHS%level = LHS%level -1
         END DO
      ELSE
         DO WHILE (RHS%level > LHS%level)
            RHS%box = mm_parent_box(RHS%box)
            RHS%bra = mm_parent_bra(RHS%bra) 
            RHS%level = RHS%level -1
         END DO
      END IF

   END SUBROUTINE translate_to_common_grid

!-------------------------------------------------------------------------------

END MODULE mm_T_pair_tests

