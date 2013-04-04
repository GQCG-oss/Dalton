MODULE mm_memory_manager_mod

! Module to record memory allocations.
! Note it only records total allocations, not
! peak allocation.

   USE mm_global_paras_mod
   use mm_mem, only: mem_alloc_fmm,mem_dealloc_fmm
   use memory_handling, only: mem_realsize,mem_intsize,mem_complexsize, mem_allocated_mem_fmm, mem_deallocated_mem_fmm
   use precision

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: mm_allocate,        &
             mm_deallocate,      &
             mm_init_mem_man,    &
             mm_close_mem_man

   INTERFACE mm_allocate
      MODULE PROCEDURE allocate_1
      MODULE PROCEDURE allocate_2
      MODULE PROCEDURE allocate_3
      MODULE PROCEDURE allocate_4
      MODULE PROCEDURE allocate_5
      MODULE PROCEDURE allocate_6
      MODULE PROCEDURE allocate_7
      MODULE PROCEDURE allocate_8
   END INTERFACE

   INTERFACE mm_deallocate
      MODULE PROCEDURE deallocate_1
      MODULE PROCEDURE deallocate_2
      MODULE PROCEDURE deallocate_3
      MODULE PROCEDURE deallocate_4
      MODULE PROCEDURE deallocate_5
      MODULE PROCEDURE deallocate_6
      MODULE PROCEDURE deallocate_7
      MODULE PROCEDURE deallocate_8
   END INTERFACE

   INTEGER, PARAMETER :: Frealk = 8
   INTEGER, PARAMETER :: Fraw_mm_paras = 92
   INTEGER, PARAMETER :: Fbox_mm_data = 4
   INTEGER, PARAMETER :: Fbox_mm_paras = 76

   INTEGER, PARAMETER :: MBytes = 1048576
   REAL(REALK), SAVE :: mem_use(NMEMDIVS)

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE mm_close_mem_man

      IMPLICIT NONE
      INTEGER :: i
      LOGICAL, PARAMETER :: PRINT_DEBUG_MEMORY=.FALSE.

      IF (PRINT_DEBUG_MEMORY) THEN 
         mem_use(:) = mem_use(:)/Mbytes
         
         WRITE(LUPRI,'(/,A)')      "  Cumulative memory allocation / MB"
         WRITE(LUPRI,'(A,/)')      "  ---------------------------------"
         DO i = 1, NMEMDIVS
            WRITE(LUPRI,'(3X,A,"  =",F11.2)') NSPACE(i), mem_use(i)
         END DO
         WRITE(LUPRI,*)
      END IF

   END SUBROUTINE mm_close_mem_man

!-------------------------------------------------------------------------------

   SUBROUTINE mm_init_mem_man

      IMPLICIT NONE
      mem_use = zero

   END SUBROUTINE mm_init_mem_man

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_1(name,array,a)
      IMPLICIT NONE
      REAL(REALK),   POINTER    :: array(:) 
      INTEGER, INTENT(IN) :: a, name
      call mem_alloc_fmm(array,a)
   END SUBROUTINE allocate_1

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_2(name,array,a,b)
      IMPLICIT NONE
      REAL(REALK),   POINTER    :: array(:,:) 
      INTEGER, INTENT(IN) :: a,b, name
      call mem_alloc_fmm(array,a,b)
   END SUBROUTINE allocate_2

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_3(name,array,a)
      IMPLICIT NONE
      TYPE(raw_mm_paras), POINTER    :: array(:) 
      INTEGER,      INTENT(IN) :: a, name
      INTEGER :: error
      integer(kind=long) :: nsize
      ALLOCATE(array(a),STAT=error)
      IF (error /= 0) CALL LSQUIT('allocation error in allocate_3',-1)
      mem_use(name) = mem_use(name) + a*Fraw_mm_paras
      
      nsize = a*Fraw_mm_paras
      call mem_allocated_mem_FMM(nsize,.TRUE.)

   END SUBROUTINE allocate_3

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_4(name,array,a)
      IMPLICIT NONE
      TYPE(box_mm_data), POINTER    :: array(:) 
      INTEGER,     INTENT(IN) :: a, name
      INTEGER :: error
      integer(KIND=long) :: nsize

      ALLOCATE(array(a),STAT=error)
      IF (error /= 0) CALL LSQUIT('allocation error in allocate_4',-1)
      mem_use(name) = mem_use(name) + a*Fbox_mm_data

      nsize = a*Fbox_mm_data
      call mem_allocated_mem_FMM(nsize,.TRUE.)

   END SUBROUTINE allocate_4

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_5(name,array,a)
      IMPLICIT NONE
      TYPE(box_mm_paras), POINTER    :: array(:) 
      INTEGER,      INTENT(IN) :: a, name
      INTEGER :: error
      integer(KIND=long) :: nsize

      ALLOCATE(array(a),STAT=error)
      IF (error /= 0) CALL LSQUIT('allocation error in allocate_5',-1)
      mem_use(name) = mem_use(name) + a*Fbox_mm_paras

      nsize = a*Fbox_mm_paras
      call mem_allocated_mem_FMM(nsize,.TRUE.)
   END SUBROUTINE allocate_5

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_6(name,array,a,b,c)
      IMPLICIT NONE
      REAL(REALK),   POINTER    :: array(:,:,:)
      INTEGER, INTENT(IN) :: a,b,c, name

      call mem_alloc_fmm(array,a,b,c)

   END SUBROUTINE allocate_6

!-------------------------------------------------------------------------------

   SUBROUTINE allocate_7(name,array,a,b)
      IMPLICIT NONE
      INTEGER, POINTER :: array(:,:)
      INTEGER, INTENT(IN) :: a,b, name
      call mem_alloc_fmm(array,a,b)

   END SUBROUTINE allocate_7

   SUBROUTINE allocate_8(array,a)
      IMPLICIT NONE
      TYPE(T_pair_single), POINTER :: array(:)
      INTEGER, INTENT(IN) :: a
      INTEGER :: error
      integer(KIND=long) :: nsize,size2

      ALLOCATE(array(a),STAT=error)
      IF (error /= 0) CALL LSQUIT('allocation error in allocate_8',-1)
      !size of T_pair_single = 7int,4realk,1Character
      size2 = mem_realsize*4+7*mem_intsize+1*mem_complexsize

      nsize = size(array)*size2
      call mem_allocated_mem_FMM(nsize,.TRUE.)

    END SUBROUTINE allocate_8
 !-------------------------------------------------------------------------------

   SUBROUTINE deallocate_1(array)
      IMPLICIT NONE
      REAL(REALK),   POINTER    :: array(:) 
      call mem_dealloc_fmm(array)
    END SUBROUTINE deallocate_1

!-------------------------------------------------------------------------------

   SUBROUTINE deallocate_2(array)
      IMPLICIT NONE
      REAL(REALK),   POINTER    :: array(:,:) 
      call mem_dealloc_fmm(array)
   END SUBROUTINE deallocate_2

!-------------------------------------------------------------------------------

   SUBROUTINE deallocate_3(array)
      IMPLICIT NONE
      TYPE(raw_mm_paras), POINTER    :: array(:) 
      INTEGER :: error
      integer(KIND=long) :: nsize

      nsize = size(array)*Fraw_mm_paras
      call mem_deallocated_mem_FMM(nsize,.TRUE.)
      deallocate(array)
   END SUBROUTINE deallocate_3

!-------------------------------------------------------------------------------

   SUBROUTINE deallocate_4(array)
      IMPLICIT NONE
      TYPE(box_mm_data), POINTER    :: array(:) 
      integer(KIND=long) :: nsize

      nsize = size(array)*Fbox_mm_data
      call mem_deallocated_mem_FMM(nsize,.TRUE.)
      deallocate(array)
   END SUBROUTINE deallocate_4

!-------------------------------------------------------------------------------

   SUBROUTINE deallocate_5(array)
      IMPLICIT NONE
      TYPE(box_mm_paras), POINTER    :: array(:) 
      integer(KIND=long) :: nsize

      nsize = size(array)*Fbox_mm_paras
      call mem_deallocated_mem_FMM(nsize,.TRUE.)
      deallocate(array)
   END SUBROUTINE deallocate_5

!-------------------------------------------------------------------------------

   SUBROUTINE deallocate_6(array)
      IMPLICIT NONE
      REAL(REALK),   POINTER    :: array(:,:,:)
      call mem_dealloc_fmm(array)
   END SUBROUTINE deallocate_6

!-------------------------------------------------------------------------------

   SUBROUTINE deallocate_7(array)
      IMPLICIT NONE
      INTEGER, POINTER :: array(:,:)
      call mem_dealloc_fmm(array)
   END SUBROUTINE deallocate_7
 !-------------------------------------------------------------------------------
   SUBROUTINE deallocate_8(array)
      IMPLICIT NONE
      TYPE(T_pair_single), POINTER :: array(:)
      integer(KIND=long) :: size2,nsize
      size2 = mem_realsize*4+7*mem_intsize+1*mem_complexsize
      nsize = size(array)*size2
      call mem_deallocated_mem_FMM(nsize,.TRUE.)
      deallocate(array)
   END SUBROUTINE deallocate_8
 !------------------------------------------------------------------------------
END MODULE mm_memory_manager_mod
