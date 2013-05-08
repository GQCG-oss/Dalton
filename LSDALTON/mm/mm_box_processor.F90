MODULE mm_box_procs

   USE mm_global_paras_mod
   USE mm_stats_mod
   use mm_mem
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_box,                          &
             mm_box_centre,                   &
             mm_branch,                       &
             mm_grain,                        &
             mm_parent_box,                   &
             mm_parent_bra,                   & 
             mm_deepest_level

   ! if the input grain is not a power of 2 we have a choice:
   ! (a) let the largest box be bigger than the longest system dimension
   ! (b) round the grain off to the nearest "good" value
   LOGICAL, PARAMETER :: USE_SYS_SIZE_AS_BIGGEST_BOX = .FALSE.

CONTAINS

!-------------------------------------------------------------------------------

   FUNCTION mm_parent_box(box)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: box(3)
      INTEGER :: mm_parent_box(3)
     
      mm_parent_box = 1+((box-1)/2)    ! get largest integer after /2

   END FUNCTION mm_parent_box

!-------------------------------------------------------------------------------

   FUNCTION mm_parent_bra(branch)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: branch
      INTEGER :: mm_parent_bra

      IF (JOIN_BRANCHES) THEN
         mm_parent_bra = (((branch-1)/2)/2) *2 +2
      ELSE
         mm_parent_bra = branch
      END IF

   END FUNCTION mm_parent_bra

!-------------------------------------------------------------------------------

   FUNCTION mm_box(centre,grain_inv,nbox_max)

      IMPLICIT NONE
      REAL(REALK),   INTENT(IN) :: centre(3)   
      REAL(REALK),   INTENT(IN) :: grain_inv   
      INTEGER, INTENT(IN) :: nbox_max 
      INTEGER :: mm_box(3)

      mm_box(:) = 1+ INT(centre(:)*grain_inv)
      IF (MAXVAL(mm_box) > nbox_max) CALL LSQUIT('bad box index',-1)
      IF (MINVAL(mm_box) < 1) CALL LSQUIT('negative/zero box indices invalid!',-1)

   END FUNCTION mm_box

!-------------------------------------------------------------------------------

   FUNCTION mm_box_centre(box,grain)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: box(3)
      REAL(REALK),   INTENT(IN) :: grain
      REAL(REALK) :: mm_box_centre(3)   

      mm_box_centre(:) = grain*(box(:)-half)

   END FUNCTION mm_box_centre

!-------------------------------------------------------------------------------
! If "Branch Free" then a single branch must be chosen so that
!  nothing is treated via boxes inside the classical radius

   FUNCTION mm_branch(extent,grain_inv,scheme)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      REAL(REALK),        INTENT(IN) :: extent, grain_inv   
      INTEGER :: mm_branch 

      IF (scheme%branch_free) THEN 
         mm_branch = MAX(WS_MIN,INT(CEILING(scheme%cls_radius*grain_inv)))
         ! we ensure branch is EVEN (for consistency with normal code)
         mm_branch = 2*((mm_branch+1)/2)
         stat_common_branch = mm_branch
      ELSE
         mm_branch = MAX(WS_MIN,INT(2*CEILING(extent*grain_inv)))
      END IF

   END FUNCTION mm_branch

!-------------------------------------------------------------------------------
! Generates deepest level required to produce input grain (smallest box size)
!  If used so that largest box = system size then
!  "grain_input"/2 < unscaled "grain" <= "grain_input"
!  else grain = input_grain throughout

   FUNCTION mm_deepest_level(scheme)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER :: mm_deepest_level
      REAL(REALK) :: x

      IF (scheme%dynamic_levels) THEN
         x = scheme%system_size/scheme%grain_input
         ! Note we round UP the integer
         mm_deepest_level = MAX(TOP_LEVEL,(1+INT(LOG(x)/LOG(two))))
      ELSE
         mm_deepest_level = scheme%NLEVEL
      END IF
      IF (mm_deepest_level > MAX_LEVEL) CALL LSQUIT('level depth too great!',-1)

   END FUNCTION mm_deepest_level

!-------------------------------------------------------------------------------

   FUNCTION mm_grain(scheme,level)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      INTEGER,      INTENT(IN) :: level
      REAL(REALK) :: mm_grain
      INTEGER :: deepest_level

      IF (USE_SYS_SIZE_AS_BIGGEST_BOX) THEN
         mm_grain = GRAIN_BUFFER*scheme%system_size*(half**level)
      ELSE IF (scheme%dynamic_levels) THEN
         ! maintain grain as input value
         deepest_level = mm_deepest_level(scheme)
         mm_grain = (scheme%grain_input)*(2**(deepest_level-level))
      ELSE
         mm_grain = GRAIN_BUFFER*scheme%system_size*(half**level)
      END IF

   END FUNCTION mm_grain

!-------------------------------------------------------------------------------

END MODULE mm_box_procs

!===============================================================================

MODULE mm_box_packer

   USE mm_global_paras_mod
   use mm_mem
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_init_pkd_paras,              &
             mm_shift_and_pack_paras

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_pkd_paras(deepest_level,scheme,raw_paras,pkd_paras)

      USE mm_memory_manager_mod, ONLY: mm_allocate, mm_deallocate
      USE mm_box_procs,      ONLY: mm_box_centre, mm_parent_box, mm_grain

      IMPLICIT NONE
      INTEGER,      INTENT(IN)    :: deepest_level
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(raw_mm_paras), INTENT(INOUT) :: raw_paras(:)
      TYPE(box_mm_paras), POINTER       :: pkd_paras(:)

      TYPE(box_mm_paras),pointer :: tmp_paras(:)
!      TYPE(box_mm_paras) :: tmp_paras(SIZE(raw_paras))
      INTEGER      :: i, l_up, box(3)
      INTEGER,pointer :: tmp_map(:)
!      INTEGER      :: tmp_map(SIZE(raw_paras))
      REAL(REALK)        :: grain

!      ALLOCATE(tmp_paras(SIZE(raw_paras)))
!      ALLOCATE(tmp_map(SIZE(raw_paras)))
      call mm_allocate(1,tmp_paras,SIZE(raw_paras))
      call mem_alloc_fmm(tmp_map,SIZE(raw_paras))

      l_up = deepest_level-1
      grain = mm_grain(scheme,l_up)
      DO i = 1, SIZE(raw_paras)
         tmp_paras(i)%box(:) = raw_paras(i)%box(:)
         tmp_paras(i)%cntr(:) = raw_paras(i)%box_cntr(:)
         tmp_paras(i)%bra = raw_paras(i)%bra
         tmp_paras(i)%level = deepest_level
         tmp_paras(i)%id = i
         box(:) = mm_parent_box(tmp_paras(i)%box(:))
         tmp_paras(i)%cntr_up(:) = mm_box_centre(box,grain) 
         tmp_paras(i)%map_up = 0   ! not defined yet
      END DO

      IF (ASSOCIATED(pkd_paras)) CALL LSQUIT ('paras should be nullified!',-1)
      IF (.FALSE.) THEN
         ! just build unpacked paras
         CALL mm_allocate(MEM_BOX_QLM,pkd_paras,INT(SIZE(raw_paras)))
         pkd_paras(:) = tmp_paras(:)
         tmp_map = (/ (i,i=1,SIZE(raw_paras)) /)
      ELSE
         ! combine paras in same box and branch (i.e. packed form)
         CALL pack_boxed_paras(tmp_paras,pkd_paras,tmp_map)
      END IF
      ! store map of packing indices (raw:boxed)
      raw_paras(:)%map_up = tmp_map(:) 

      call mem_dealloc_fmm(tmp_map)
      call mm_deallocate(tmp_paras)

   END SUBROUTINE mm_init_pkd_paras

!-------------------------------------------------------------------------------

   SUBROUTINE mm_shift_and_pack_paras(level,scheme,paras_in,paras_out)

      USE mm_memory_manager_mod, ONLY: mm_allocate
      USE mm_box_procs,      ONLY: mm_box_centre, mm_grain,                &
                                   mm_parent_box, mm_parent_bra

      IMPLICIT NONE
      INTEGER,      INTENT(IN)    :: level
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      TYPE(box_mm_paras), INTENT(INOUT) :: paras_in(:)
      TYPE(box_mm_paras), POINTER       :: paras_out(:)

      TYPE(box_mm_paras) :: tmp_paras(SIZE(paras_in))
      INTEGER      :: i, l_down, l_up, box(3)
      INTEGER      :: tmp_map(SIZE(paras_in))
      REAL(REALK)        :: grain, grain_up

      ! build tmp array for unpacked paras at next level up
      l_down = level+1
      l_up = level-1
      grain = mm_grain(scheme,level)
      grain_up = mm_grain(scheme,l_up)
      DO i = 1, SIZE(paras_in)
         tmp_paras(i)%box = mm_parent_box(paras_in(i)%box)
         tmp_paras(i)%cntr(:) = mm_box_centre(tmp_paras(i)%box,grain)
         tmp_paras(i)%bra = mm_parent_bra(paras_in(i)%bra)
         tmp_paras(i)%level = level
         tmp_paras(i)%id = i
         box(:) = mm_parent_box(tmp_paras(i)%box(:))
         tmp_paras(i)%cntr_up(:) = mm_box_centre(box,grain_up) 
         tmp_paras(i)%map_up = 0   ! not defined yet
      END DO

      IF (.FALSE.) THEN
         ! just build unpacked paras
         CALL mm_allocate(MEM_BOX_QLM,paras_out,INT(SIZE(paras_in)))
         paras_out(:) = tmp_paras(:)
         tmp_map = (/ (i,i=1,SIZE(paras_in)) /)
      ELSE
         ! combine paras in same box and branch (i.e. packed form)
         CALL pack_boxed_paras(tmp_paras,paras_out,tmp_map)
      END IF
      ! store map of packing indices (raw:boxed)
      paras_in(:)%map_up = tmp_map(:) 

   END SUBROUTINE mm_shift_and_pack_paras

!-------------------------------------------------------------------------------

   SUBROUTINE pack_boxed_paras(paras_in,paras_out,map) 

      USE mm_memory_manager_mod, ONLY: mm_allocate, mm_deallocate 

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(INOUT) :: paras_in(:)
      TYPE(box_mm_paras), POINTER       :: paras_out(:)
      INTEGER,      INTENT(OUT)   :: map(:)

      TYPE(box_mm_paras),pointer :: tmp_paras(:)
      INTEGER :: i, lo, k, test
      LOGICAL :: new_entry

      call mm_ALLOCATE(1,tmp_paras,SIZE(paras_in))

      map = 0
      tmp_paras(1) = paras_in(1)
      tmp_paras(1)%id = 1 
      map(paras_in(1)%id) = 1
      lo = 1
      k = 1

      DO i = 2, SIZE(paras_in)
         new_entry = .FALSE.
         ! we've sorted wrt branches, so we can pack each batch separately
         IF (paras_in(i)%bra /= paras_in(i-1)%bra) THEN
            IF (paras_in(i)%bra < paras_in(i-1)%bra) STOP 'must sort data!' 
            new_entry = .TRUE.
            lo = k+1 
         ELSE
            ! same branch of parameters, so test if also same box
            test = test_new_parameters(paras_in(i),tmp_paras,lo,k)
            IF (test /= 0) THEN
               ! same box/branch, so pack and store mapping at previous level
               map(paras_in(i)%id) = test  !id for moment at upper level
               CYCLE
            END IF
            new_entry = .TRUE.
         END IF
         IF (new_entry) THEN
            k = k+1 
            tmp_paras(k) = paras_in(i)
            tmp_paras(k)%id = k
            map(paras_in(i)%id) = k 
         END IF
      END DO

      ! store packed parameters in exactly allocated array
      CALL mm_allocate(MEM_BOX_QLM,paras_out,k)
      paras_out = tmp_paras(1:k)
      call mm_DEALLOCATE(tmp_paras)

   END SUBROUTINE pack_boxed_paras

!-------------------------------------------------------------------------------

   FUNCTION test_new_parameters(test_item,ref_items,lo,hi)

      IMPLICIT NONE
      TYPE(box_mm_paras), INTENT(IN) :: test_item
      TYPE(box_mm_paras), INTENT(IN) :: ref_items(:)
      INTEGER,      INTENT(IN) :: lo, hi
      INTEGER :: i, test_new_parameters

      test_new_parameters = 0
      DO i = lo, hi
         IF (test_item%box(1) /= ref_items(i)%box(1)) CYCLE
         IF (test_item%box(2) /= ref_items(i)%box(2)) CYCLE
         IF (test_item%box(3) /= ref_items(i)%box(3)) CYCLE
         ! boxes identical
         test_new_parameters = i
         RETURN
      END DO

   END FUNCTION test_new_parameters

!-------------------------------------------------------------------------------

END MODULE mm_box_packer

!===============================================================================

MODULE mm_level_optimiser

   USE mm_global_paras_mod
   use mm_mem
   IMPLICIT NONE
   PRIVATE

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE a

   END SUBROUTINE a

!-------------------------------------------------------------------------------

END MODULE mm_level_optimiser

