MODULE mm_sorting

! Sorting module
! 1. sort raw_mm_paras wrt %bra
!    with a very simple O(N^2) insertion sort algorithm
!    or using the average O(NlogN) Quicksort algorithm
! 2. sort array of interaction T-pairs using Quicksort wrt
!     (a) RHS index
!     (b) T-vector component
!     (c) T-vector modulus (ratio)

   USE mm_global_paras_mod
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_quicksort_wrt_branches,      &
             mm_quicksort_wrt_RHS,           &
             mm_quicksort_wrt_vector,        &
             mm_quicksort_wrt_ratio

   ! criteria for switching from Quicksort to Insertion-Sort algorithms
   INTEGER, PARAMETER :: quicksort_CUTOFF = 10 

CONTAINS

!-------------------------------------------------------------------------------
! Insertion Sort
!----------------

   SUBROUTINE insertion_sort_wrt_branches(arr)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: arr(:)
      TYPE(raw_mm_paras) :: tmp
      INTEGER :: i, j

      DO i = 1, SIZE(arr) 
         tmp = arr(i)
         DO j = (i-1), 0, -1
            IF (j>0) THEN 
               IF (arr(j)%bra > tmp%bra) THEN
                 arr(j+1) = arr(j) 
               ELSE
                 arr(j+1) = tmp
                 EXIT
               ENDIF 
            ELSE
               arr(j+1) = tmp
               EXIT
            END IF
         END DO
      END DO

   END SUBROUTINE insertion_sort_wrt_branches

!-------------------------------------------------------------------------------
! Quicksort
!-----------

   SUBROUTINE bra_swap_elements(arr,i,j)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER,      INTENT(IN)    :: i, j
      TYPE(raw_mm_paras) :: tmp

      tmp = arr(i)
      arr(i) = arr(j)
      arr(j) = tmp

   END SUBROUTINE bra_swap_elements

!-------------------------------------------------------------------------------
! subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

   FUNCTION bra_median_of_three(arr,left,right)

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: arr(:)
      INTEGER,      INTENT(IN)    :: left, right
      INTEGER :: bra_median_of_three
      INTEGER :: cntr

      cntr = (left + right)/2
      IF(arr(left)%bra > arr(cntr)%bra)  CALL bra_swap_elements(arr,left,cntr)
      IF(arr(left)%bra > arr(right)%bra) CALL bra_swap_elements(arr,left,right)
      IF(arr(cntr)%bra > arr(right)%bra) CALL bra_swap_elements(arr,cntr,right)
      CALL bra_swap_elements(arr,cntr,right-1)
      bra_median_of_three = arr(right-1)%bra

   END FUNCTION bra_median_of_three

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE mm_quicksort_wrt_branches(arr) 

      IMPLICIT NONE
      TYPE(raw_mm_paras), INTENT(INOUT) :: arr(:)

      INTEGER :: left, right
      INTEGER :: i,j, pivot

      left = 1
      right = SIZE(arr)

      IF (left + quicksort_CUTOFF > right) THEN
         ! array too small for Quicksort to be efficient, so don't use it
         CALL insertion_sort_wrt_branches(arr)
         RETURN
      END IF

      pivot = bra_median_of_three(arr,left,right)
      ! partition array elements about chosen pivot
      i = left
      j = right-2
      DO ! until left, right pointers cross
         ! increment left until large element found
         DO WHILE ( arr(i)%bra < pivot )
            i = i+1
         END DO
         ! increment right until small element found
         DO WHILE ( arr(j)%bra > pivot )
            j = j-1
         END DO
         IF (i < j) THEN
            ! swap items (i,j)
            CALL bra_swap_elements(arr,i,j)
            ! now increment (i,j) to avoid infinite loop if both = pivot
            i = i+1
            j = j-1
         ELSE
            EXIT
         END IF
      END DO
      ! i is now the position where the pivot should be 
      ! swap pivot back to middle of array
      CALL bra_swap_elements(arr,i,right-1)

      ! now sort sub-arrays either side of pivot
      CALL mm_quicksort_wrt_branches( arr( left:(i-1) ) )
      CALL mm_quicksort_wrt_branches( arr( (i+1):right ) )

   END SUBROUTINE mm_quicksort_wrt_branches

!-------------------------------------------------------------------------------
! Insertion sort for sorting of T-pair batch wrt RHS moments
!-----------------------------------------------------------

   SUBROUTINE insertion_sort_wrt_RHS(arr)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      TYPE(T_pair_single) :: tmp
      INTEGER :: i, j

      DO i = 1, SIZE(arr) 
         tmp = arr(i)
         DO j = (i-1), 0, -1
            IF (j>0) THEN
               IF  (arr(j)%paras%RHS_id > tmp%paras%RHS_id) THEN
                 arr(j+1) = arr(j) 
               ELSE
                 arr(j+1) = tmp
                 EXIT
               END IF
            ELSE
               arr(j+1) = tmp
               EXIT
            END IF
         END DO
      END DO

   END SUBROUTINE insertion_sort_wrt_RHS

!-------------------------------------------------------------------------------
! Quicksort for sorting of T-pair batch wrt RHS moments
!------------------------------------------------------

   SUBROUTINE RHS_swap_elements(arr,i,j)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER,       INTENT(IN)    :: i, j
      TYPE(T_pair_single) :: tmp

      tmp = arr(i)
      arr(i) = arr(j)
      arr(j) = tmp

   END SUBROUTINE RHS_swap_elements

!-------------------------------------------------------------------------------
! subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

   FUNCTION RHS_median_of_three(arr,left,right)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER,       INTENT(IN)    :: left, right
      INTEGER :: RHS_median_of_three
      INTEGER :: cntr, id_left, id_right, id_cntr

      cntr = (left + right)/2

      IF (arr(left)%paras%RHS_id > arr(cntr)%paras%RHS_id)  &
         CALL RHS_swap_elements(arr,left,cntr)
      IF (arr(left)%paras%RHS_id > arr(right)%paras%RHS_id)  &
         CALL RHS_swap_elements(arr,left,right)
      IF (arr(cntr)%paras%RHS_id > arr(right)%paras%RHS_id)  &
         CALL RHS_swap_elements(arr,cntr,right)

      CALL RHS_swap_elements(arr,cntr,right-1)
      RHS_median_of_three = arr(right-1)%paras%RHS_id

   END FUNCTION RHS_median_of_three

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE mm_quicksort_wrt_RHS(arr) 

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)

      INTEGER :: left, right
      INTEGER :: i,j, pivot

      left = 1
      right = SIZE(arr)

      IF (left + quicksort_CUTOFF > right) THEN
         ! array too small for Quicksort to be efficient, so don't use it
         CALL insertion_sort_wrt_RHS(arr)
         RETURN
      END IF

      pivot = RHS_median_of_three(arr,left,right)
      ! partition array elements about chosen pivot
      i = left
      j = right-2
      DO ! until left, right pointers cross
         ! increment left until large element found
         DO WHILE ( arr(i)%paras%RHS_id < pivot )
            i = i+1
         END DO
         ! increment right until small element found
         DO WHILE ( arr(j)%paras%RHS_id > pivot )
            j = j-1
         END DO
         IF (i < j) THEN
            ! swap items (i,j)
            CALL RHS_swap_elements(arr,i,j)
            ! now increment (i,j) to avoid infinite loop if both = pivot
            i = i+1
            j = j-1
         ELSE
            EXIT
         END IF
      END DO
      ! i is now the position where the pivot should be 
      ! swap pivot back to middle of array
      CALL RHS_swap_elements(arr,i,right-1)

      ! now sort sub-arrays either side of pivot
      CALL mm_quicksort_wrt_RHS( arr( left:(i-1) ) )
      CALL mm_quicksort_wrt_RHS( arr( (i+1):right ) )

   END SUBROUTINE mm_quicksort_wrt_RHS

!-------------------------------------------------------------------------------
! Insertion sort for sorting of T-pair batch wrt 1 component of T-vector 
!-----------------------------------------------------------------------

   SUBROUTINE insertion_sort_wrt_vector(arr,xyz)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER,       INTENT(IN)    :: xyz

      TYPE(T_pair_single) :: tmp
      INTEGER :: i, j

      DO i = 1, SIZE(arr) 
         tmp = arr(i)
         DO j = (i-1), 0, -1
            IF (j>0) THEN
               IF(arr(j)%r_ab(xyz) > tmp%r_ab(xyz)) THEN
                  arr(j+1) = arr(j) 
               ELSE
                  arr(j+1) = tmp
                  EXIT
               ENDIF
            ELSE
               arr(j+1) = tmp
               EXIT
            END IF
         END DO
      END DO

   END SUBROUTINE insertion_sort_wrt_vector

!-------------------------------------------------------------------------------
! Quicksort for sorting of T-pair batch wrt 1 component of T-vector
!------------------------------------------------------------------

   SUBROUTINE vector_swap_elements(arr,i,j)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER,       INTENT(IN)    :: i, j
      TYPE(T_pair_single) :: tmp

      tmp = arr(i)
      arr(i) = arr(j)
      arr(j) = tmp

   END SUBROUTINE vector_swap_elements

!-------------------------------------------------------------------------------
! subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

   FUNCTION vector_median_of_three(arr,left,right,xyz)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER,       INTENT(IN)    :: left, right
      INTEGER,       INTENT(IN)    :: xyz

      REAL(REALK)   :: vector_median_of_three
      INTEGER :: cntr

      cntr = (left + right)/2

      IF (arr(left)%r_ab(xyz) > arr(cntr)%r_ab(xyz))  &
         CALL vector_swap_elements(arr,left,cntr)
      IF (arr(left)%r_ab(xyz) > arr(right)%r_ab(xyz))  &
         CALL vector_swap_elements(arr,left,right)
      IF (arr(cntr)%r_ab(xyz) > arr(right)%r_ab(xyz))  &
         CALL vector_swap_elements(arr,cntr,right)

      CALL vector_swap_elements(arr,cntr,right-1)
      vector_median_of_three = arr(right-1)%r_ab(xyz)

   END FUNCTION vector_median_of_three

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE mm_quicksort_wrt_vector(arr,xyz) 

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER,       INTENT(IN)    :: xyz

      INTEGER :: left, right
      INTEGER :: i,j
      REAL(REALK)   :: pivot

      left = 1
      right = SIZE(arr)

      IF (left + quicksort_CUTOFF > right) THEN
         ! array too small for Quicksort to be efficient, so don't use it
         CALL insertion_sort_wrt_vector(arr,xyz)
         RETURN
      END IF

      pivot = vector_median_of_three(arr,left,right,xyz)
      ! partition array elements about chosen pivot
      i = left
      j = right-2
      DO ! until left, right pointers cross
         ! increment left until large element found
         DO WHILE ( arr(i)%r_ab(xyz) < pivot )
            i = i+1
         END DO
         ! increment right until small element found
         DO WHILE ( arr(j)%r_ab(xyz) > pivot )
            j = j-1
         END DO
         IF (i < j) THEN
            ! swap items (i,j)
            CALL vector_swap_elements(arr,i,j)
            ! now increment (i,j) to avoid infinite loop if both = pivot
            i = i+1
            j = j-1
         ELSE
            EXIT
         END IF
      END DO
      ! i is now the position where the pivot should be 
      ! swap pivot back to middle of array
      CALL vector_swap_elements(arr,i,right-1)

      ! now sort sub-arrays either side of pivot
      CALL mm_quicksort_wrt_vector( arr( left:(i-1) ),xyz )
      CALL mm_quicksort_wrt_vector( arr( (i+1):right ),xyz )

   END SUBROUTINE mm_quicksort_wrt_vector

!-------------------------------------------------------------------------------
! Insertion sort for sorting of T-pair batch wrt T-vector modulus (ratio)
!------------------------------------------------------------------------

   SUBROUTINE insertion_sort_wrt_ratio(arr)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)

      TYPE(T_pair_single) :: tmp
      INTEGER :: i, j

      DO i = 1, SIZE(arr) 
         tmp = arr(i)
         DO j = (i-1), 0, -1
            IF (j>0) THEN
               IF (arr(j)%paras%ratio > tmp%paras%ratio) THEN
                  arr(j+1) = arr(j)
               ELSE
                  arr(j+1) = tmp
                  EXIT
               END IF
            ELSE
               arr(j+1) = tmp
               EXIT
            END IF
         END DO
      END DO

   END SUBROUTINE insertion_sort_wrt_ratio

!-------------------------------------------------------------------------------
! Quicksort for sorting of T-pair batch wrt T-vector modulus (ratio)
!-------------------------------------------------------------------

   SUBROUTINE ratio_swap_elements(arr,i,j)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER,       INTENT(IN)    :: i, j
      TYPE(T_pair_single) :: tmp

      tmp = arr(i)
      arr(i) = arr(j)
      arr(j) = tmp

   END SUBROUTINE ratio_swap_elements

!-------------------------------------------------------------------------------
! subroutine to return the median of first, last, and centre array elements.
! order these and hide the pivot (median) as the second to last array element
! (acts as a sentinel whilst the last is the largest and put in correct place.)

   FUNCTION ratio_median_of_three(arr,left,right)

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)
      INTEGER,       INTENT(IN)    :: left, right

      REAL(REALK)   :: ratio_median_of_three
      REAL(REALK)   :: xyz_cntr, xyz_left, xyz_right
      INTEGER :: cntr

      cntr = (left + right)/2

      IF (arr(left)%paras%ratio > arr(cntr)%paras%ratio)  &
         CALL ratio_swap_elements(arr,left,cntr)
      IF (arr(left)%paras%ratio > arr(right)%paras%ratio)  &
         CALL ratio_swap_elements(arr,left,right)
      IF (arr(cntr)%paras%ratio > arr(right)%paras%ratio)  &
         CALL ratio_swap_elements(arr,cntr,right)

      CALL ratio_swap_elements(arr,cntr,right-1)
      ratio_median_of_three = arr(right-1)%paras%ratio

   END FUNCTION ratio_median_of_three

!-------------------------------------------------------------------------------

   RECURSIVE SUBROUTINE mm_quicksort_wrt_ratio(arr) 

      IMPLICIT NONE
      TYPE(T_pair_single), INTENT(INOUT) :: arr(:)

      INTEGER :: left, right
      INTEGER :: i,j
      REAL(REALK)   :: pivot

      left = 1
      right = SIZE(arr)

      IF (left + quicksort_CUTOFF > right) THEN
         ! array too small for Quicksort to be efficient, so don't use it
         CALL insertion_sort_wrt_ratio(arr)
         RETURN
      END IF

      pivot = ratio_median_of_three(arr,left,right)
      ! partition array elements about chosen pivot
      i = left
      j = right-2
      DO ! until left, right pointers cross
         ! increment left until large element found
         DO WHILE ( arr(i)%paras%ratio < pivot )
            i = i+1
         END DO
         ! increment right until small element found
         DO WHILE ( arr(j)%paras%ratio > pivot )
            j = j-1
         END DO
         IF (i < j) THEN
            ! swap items (i,j)
            CALL ratio_swap_elements(arr,i,j)
            ! now increment (i,j) to avoid infinite loop if both = pivot
            i = i+1
            j = j-1
         ELSE
            EXIT
         END IF
      END DO
      ! i is now the position where the pivot should be 
      ! swap pivot back to middle of array
      CALL ratio_swap_elements(arr,i,right-1)

      ! now sort sub-arrays either side of pivot
      CALL mm_quicksort_wrt_ratio( arr( left:(i-1) ) )
      CALL mm_quicksort_wrt_ratio( arr( (i+1):right ) )

   END SUBROUTINE mm_quicksort_wrt_ratio

!-------------------------------------------------------------------------------

END MODULE mm_sorting


