!Monitor memory NON THREAD SAFE
MODULE IchorMemory
  use IchorCommonModule
  use IchorprecisionModule
   private
   public set_ichor_memvar
   public retrieve_ichor_memvar
   public stats_ichor_mem
   public mem_ichor_alloc
   public mem_ichor_dealloc
   public mem_ichor_alloc_dryrun
   public mem_ichor_dealloc_dryrun
   public ADD_OMP_MEM
   public REMOVE_OMP_MEM
!GLOBAL VARIABLES
   integer(KIND=long),save :: mem_allocated_ichor, max_mem_used_ichor  !Count all memory
   integer(KIND=long),save :: maxMemLimit_ichor !maximum allowed memory usage
   !Count 'real' memory, integral code
   integer(KIND=long),save :: mem_allocated_real, max_mem_used_real
   !Count 'integer' memory, integral code
   integer(KIND=long),save :: mem_allocated_integer, max_mem_used_integer
   !Count 'logical' memory, integral code
   integer(KIND=long),save :: mem_allocated_logical, max_mem_used_logical
! sizes
   integer(KIND=long),parameter :: mem_realsize=8
   integer(KIND=long),parameter :: mem_logicalsize=4
#if VAR_INT64
   integer(KIND=long),parameter :: mem_intsize=8_long
#else
   integer(KIND=long),parameter :: mem_intsize=4_long
#endif
!Interfaces for counting memory 
INTERFACE mem_ichor_alloc_dryrun
  MODULE PROCEDURE alloc_dryrun_1dim, alloc_dryrun_2dim, &
       & alloc_dryrun_3dim, alloc_dryrun_4dim

END INTERFACE
!
INTERFACE mem_ichor_dealloc_dryrun
  MODULE PROCEDURE dealloc_dryrun_1dim, dealloc_dryrun_2dim, &
       & dealloc_dryrun_3dim, dealloc_dryrun_4dim
END INTERFACE

INTERFACE mem_ichor_alloc
  MODULE PROCEDURE real_allocate_1dim, real_allocate_2dim, &
       & real_allocate_3dim, real_allocate_4dim, real_allocate_5dim, &
       & int_allocate_1dim, int_allocate_2dim, int_allocate_3dim, &
       & logic_allocate_1dim, logic_allocate_2dim

END INTERFACE
!
INTERFACE mem_ichor_dealloc
  MODULE PROCEDURE real_deallocate_1dim, real_deallocate_2dim,  &
     & real_deallocate_3dim, real_deallocate_4dim, real_deallocate_5dim, &
     & int_deallocate_1dim, int_deallocate_2dim, int_deallocate_3dim,&
     & logic_deallocate_1dim, logic_deallocate_2dim
END INTERFACE
!
CONTAINS
subroutine set_ichor_memvar(MaxMemAllocated,MemAllocated,MaxMem)
implicit none
integer(KIND=long),intent(in) :: MaxMemAllocated,MemAllocated,MaxMem
max_mem_used_ichor = MaxMemAllocated
mem_allocated_ichor = MemAllocated
IF(MaxMem.EQ.0)THEN
   maxMemLimit_ichor = HUGE(MaxMem)
ELSE
   maxMemLimit_ichor = MaxMem
ENDIF
mem_allocated_real = 0
max_mem_used_real = 0
mem_allocated_integer = 0
max_mem_used_integer = 0
mem_allocated_logical = 0
max_mem_used_logical = 0
end subroutine set_ichor_memvar

subroutine retrieve_ichor_memvar(MaxMemAllocated,MemAllocated)
implicit none
integer(KIND=long),intent(inout) :: MaxMemAllocated,MemAllocated
MemAllocated = mem_allocated_ichor
MaxMemAllocated = max_mem_used_ichor
end subroutine retrieve_ichor_memvar

subroutine stats_ichor_mem(lupri)
  implicit none
  !> Logical unit number for output file.
  integer,intent(in) :: lupri
  WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
  WRITE(LUPRI,'("            Ichor Memory statistics          ")')
  WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
  WRITE(LUPRI,'("  Allocated memory (TOTAL):           ",i9," byte  &
       &- Should be zero - otherwise a leakage is present")') mem_allocated_ichor
  WRITE(LUPRI,'("  Allocated memory (real):            ",i9," byte  &
       &- Should be zero - otherwise a leakage is present")') mem_allocated_real
  WRITE(LUPRI,'("  Allocated memory (integer):         ",i9," byte  &
       &- Should be zero - otherwise a leakage is present")') mem_allocated_integer
  WRITE(LUPRI,'("  Allocated memory (logical):         ",i9," byte  &
       &- Should be zero - otherwise a leakage is present")') mem_allocated_logical
  call print_ichor_maxmem(lupri,max_mem_used_ichor,'TOTAL')
  CALL print_ichor_maxmem(lupri,max_mem_used_real,'real(realk)')
  CALL print_ichor_maxmem(lupri,max_mem_used_integer,'integer')
  CALL print_ichor_maxmem(lupri,max_mem_used_logical,'logical')
  WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
  WRITE(LUPRI,*)
end subroutine stats_ichor_mem

!> \brief Print routine for memory statistics (max. allocated memory).
!> \author S. Host
!> \date 2009
subroutine print_ichor_maxmem(lupri,max_mem_used,STRING)
  implicit none
  !> Logical unit number for output file
  integer :: lupri
  !> Amount of memory used
  integer(KIND=long) :: max_mem_used
  !> Name of data type for which memory usage is printed
  Character*(*)    ::  STRING
  character(len=27) :: trimstring
  !Character(len=3) ::  FORMATSTRING1
  !Character(len=23) ::  FORMATSTRING2
  trimstring = adjustl(string)
  if (max_mem_used < 100) then !Divide by 1 to typecast real
     WRITE(LUPRI,'("  Max allocated memory, ", a27, f9.3, " Byte")') trimstring, max_mem_used/1.0E0_realk
  else if (max_mem_used < 1000000) then
     WRITE(LUPRI,'("  Max allocated memory, ", a27, f9.3, " kB")') trimstring, max_mem_used/1000.0E0_realk
  else if (max_mem_used < 1000000000) then
     WRITE(LUPRI,'("  Max allocated memory, ", a27, f9.3, " MB")') trimstring, max_mem_used/(1000.0E0_realk*1000.0E0_realk)
  else
     WRITE(LUPRI,'("  Max allocated memory, ", a27, f9.3, " GB")') trimstring, &
          & max_mem_used/(1000.0E0_realk*1000.0E0_realk*1000.0E0_realk)
  endif
end subroutine print_ichor_maxmem

SUBROUTINE alloc_dryrun_1dim(n1)
  implicit none
  integer :: n1
  integer (kind=long) :: nsize
  nsize = n1*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE ALLOC_DRYRUN_1DIM

SUBROUTINE alloc_dryrun_2dim(n1,n2)
  implicit none
  integer :: n1,n2
  integer (kind=long) :: nsize
  nsize = n1*n2*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE ALLOC_DRYRUN_2DIM

SUBROUTINE alloc_dryrun_3dim(n1,n2,n3)
  implicit none
  integer :: n1,n2,n3
  integer (kind=long) :: nsize
  nsize = n1*n2*n3*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE ALLOC_DRYRUN_3DIM

SUBROUTINE alloc_dryrun_4dim(n1,n2,n3,n4)
  implicit none
  integer :: n1,n2,n3,n4
  integer (kind=long) :: nsize
  nsize = n1*n2*n3*n4*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE ALLOC_DRYRUN_4DIM

SUBROUTINE dealloc_dryrun_1dim(n1)
  implicit none
  integer :: n1
  integer (kind=long) :: nsize
  nsize = n1*mem_realsize
  call mem_deallocated_ichor_mem_real(nsize)
END SUBROUTINE DEALLOC_DRYRUN_1DIM

SUBROUTINE dealloc_dryrun_2dim(n1,n2)
  implicit none
  integer :: n1,n2
  integer (kind=long) :: nsize
  nsize = n1*n2*mem_realsize
  call mem_deallocated_ichor_mem_real(nsize)
END SUBROUTINE DEALLOC_DRYRUN_2DIM

SUBROUTINE dealloc_dryrun_3dim(n1,n2,n3)
  implicit none
  integer :: n1,n2,n3
  integer (kind=long) :: nsize
  nsize = n1*n2*n3*mem_realsize
  call mem_deallocated_ichor_mem_real(nsize)
END SUBROUTINE DEALLOC_DRYRUN_3DIM

SUBROUTINE dealloc_dryrun_4dim(n1,n2,n3,n4)
  implicit none
  integer :: n1,n2,n3,n4
  integer (kind=long) :: nsize
  nsize = n1*n2*n3*n4*mem_realsize
  call mem_deallocated_ichor_mem_real(nsize)
END SUBROUTINE DEALLOC_DRYRUN_4DIM

SUBROUTINE real_allocate_1dim(A)
  implicit none
  REAL(REALK) :: A(:)
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE real_allocate_1dim

SUBROUTINE real_allocate_2dim(A)
  implicit none
  REAL(REALK) :: A(:,:)
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE real_allocate_2dim

SUBROUTINE ADD_OMP_MEM(nsize1)
  implicit none
  integer :: nsize1
  integer (kind=long) :: nsize
  nsize = nsize1
  call mem_allocated_ichor_mem_real2(nsize)
END SUBROUTINE ADD_OMP_MEM

SUBROUTINE REMOVE_OMP_MEM(nsize1)
  implicit none
  integer :: nsize1
  integer (kind=long) :: nsize
  nsize = nsize1
  call mem_deallocated_ichor_mem_real2(nsize)
END SUBROUTINE REMOVE_OMP_MEM

SUBROUTINE real_allocate_3dim(A)
  implicit none
  REAL(REALK) :: A(:,:,:)
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE real_allocate_3dim

SUBROUTINE real_allocate_4dim(A)
  implicit none
  REAL(REALK) :: A(:,:,:,:)
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE real_allocate_4dim

SUBROUTINE real_allocate_5dim(A)
  implicit none
  REAL(REALK) :: A(:,:,:,:,:)
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE real_allocate_5dim

SUBROUTINE real_deallocate_1dim(A)
  implicit none
  REAL(REALK) :: A(:)
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_deallocated_ichor_mem_real(nsize)
END SUBROUTINE real_deallocate_1dim

SUBROUTINE real_deallocate_2dim(A)
implicit none
REAL(REALK) :: A(:,:)
integer (kind=long) :: nsize
nsize = size(A,KIND=long)*mem_realsize
call mem_deallocated_ichor_mem_real(nsize)
END SUBROUTINE real_deallocate_2dim

SUBROUTINE real_deallocate_3dim(A)
implicit none
REAL(REALK) :: A(:,:,:)
integer (kind=long) :: nsize
nsize = size(A,KIND=long)*mem_realsize
call mem_deallocated_ichor_mem_real(nsize)
END SUBROUTINE real_deallocate_3dim

SUBROUTINE real_deallocate_4dim(A)
implicit none
REAL(REALK) :: A(:,:,:,:)
integer (kind=long) :: nsize
nsize = size(A,KIND=long)*mem_realsize
call mem_deallocated_ichor_mem_real(nsize)
END SUBROUTINE real_deallocate_4dim

SUBROUTINE real_deallocate_5dim(A)
implicit none
REAL(REALK) :: A(:,:,:,:,:)
integer (kind=long) :: nsize
nsize = size(A,KIND=long)*mem_realsize
call mem_deallocated_ichor_mem_real(nsize)
END SUBROUTINE real_deallocate_5dim

SUBROUTINE int_allocate_1dim(I)
implicit none
INTEGER             :: I(:)
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_ichor_mem_integer(nsize)
END SUBROUTINE int_allocate_1dim

SUBROUTINE int_allocate_2dim(I)
implicit none
INTEGER    :: I(:,:)
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_ichor_mem_integer(nsize)
END SUBROUTINE int_allocate_2dim

SUBROUTINE int_allocate_3dim(I)
implicit none
INTEGER :: I(:,:,:)
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_ichor_mem_integer(nsize)
END SUBROUTINE int_allocate_3dim

SUBROUTINE int_allocate_4dim(I)
implicit none
INTEGER :: I(:,:,:,:)
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_ichor_mem_integer(nsize)
END SUBROUTINE int_allocate_4dim

SUBROUTINE int_deallocate_1dim(I)
implicit none
INTEGER :: I(:)
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_deallocated_ichor_mem_integer(nsize)
END SUBROUTINE int_deallocate_1dim

SUBROUTINE int_deallocate_2dim(I)
implicit none
INTEGER :: I(:,:)
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_deallocated_ichor_mem_integer(nsize)
END SUBROUTINE int_deallocate_2dim

SUBROUTINE int_deallocate_3dim(I)
implicit none
INTEGER :: I(:,:,:)
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_deallocated_ichor_mem_integer(nsize)
END SUBROUTINE int_deallocate_3dim

SUBROUTINE logic_allocate_1dim(L)
implicit none
logical     :: L(:)
integer (kind=long) :: nsize
nsize = size(L,KIND=long)*mem_logicalsize
call mem_allocated_ichor_mem_logical(nsize)
END SUBROUTINE logic_allocate_1dim

SUBROUTINE logic_allocate_2dim(L)
implicit none
logical    :: L(:,:)
integer (kind=long) :: nsize
nsize = size(L,KIND=long)*mem_logicalsize
call mem_allocated_ichor_mem_logical(nsize)
END SUBROUTINE logic_allocate_2dim

SUBROUTINE logic_deallocate_1dim(L)
implicit none
logical :: L(:)
integer (kind=long) :: nsize
nsize = size(L,KIND=long)*mem_logicalsize
call mem_deallocated_ichor_mem_logical(nsize)
END SUBROUTINE logic_deallocate_1dim

SUBROUTINE logic_deallocate_2dim(L)
implicit none
LOGICAL :: L(:,:)
integer (kind=long) :: nsize
nsize = size(L,KIND=long)*mem_logicalsize
call mem_deallocated_ichor_mem_logical(nsize)
END SUBROUTINE logic_deallocate_2dim

!----- MEMORY HANDLING -----!

subroutine mem_allocated_ichor_mem_real(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  mem_allocated_real = mem_allocated_real + nsize
  if (mem_allocated_real < 0) then
     write(*,*) 'Real memory negative! mem_allocated_real =', mem_allocated_real
     write(*,*) 'Real memory negative! nsize =', nsize
     call IchorQuit('Error in mem_allocated_ichor_mem_real - probably integer overflow!',-1)
  endif
  max_mem_used_real = MAX(max_mem_used_real,mem_allocated_real)
  !Count also the total memory:
  mem_allocated_ichor = mem_allocated_ichor  + nsize
  if (mem_allocated_ichor < 0) then
     write(*,*) 'Total memory negative! mem_allocated_ichor =', mem_allocated_ichor
     call IchorQuit('Error in mem_allocated_ichor_mem_real - probably integer overflow!',-1)
  endif
  max_mem_used_ichor = MAX(max_mem_used_ichor,mem_allocated_ichor)
  IF(max_mem_used_ichor.GT.maxMemLimit_ichor)THEN
     print*,'Maximum Memory Limit',maxMemLimit_ichor
     print*,'Memory Usage        ',max_mem_used_ichor
     call IchorQuit('Error in mem_allocated_ichor_mem_real - Maximum Memory Limit Exceeded!',-1)     
  ENDIF
end subroutine mem_allocated_ichor_mem_real

subroutine mem_deallocated_ichor_mem_real(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  
  mem_allocated_real = mem_allocated_real - nsize
  if (mem_allocated_real < 0) then
     write(*,*) 'Real memory negative! mem_allocated_real =', mem_allocated_real
     call IchorQuit('Error in mem_deallocated_ichor_mem_real - something wrong with deallocation!',-1)
  endif
  !Count also the total memory:
  mem_allocated_ichor = mem_allocated_ichor - nsize
  if (mem_allocated_ichor < 0) then
     write(*,*) 'Total memory negative! mem_allocated_ichor =', mem_allocated_ichor
     call IchorQuit('Error in mem_deallocated_ichor_mem_real - something wrong with deallocation!',-1)
  endif
end subroutine mem_deallocated_ichor_mem_real

subroutine mem_allocated_ichor_mem_real2(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  mem_allocated_real = mem_allocated_real + nsize
  max_mem_used_real = MAX(max_mem_used_real,mem_allocated_real)
  mem_allocated_ichor = mem_allocated_ichor  + nsize
  max_mem_used_ichor = MAX(max_mem_used_ichor,mem_allocated_ichor)
  IF(max_mem_used_ichor.GT.maxMemLimit_ichor)THEN
     print*,'Maximum Memory Limit',maxMemLimit_ichor
     print*,'Memory Usage        ',max_mem_used_ichor
     call IchorQuit('Error in mem_allocated_ichor_mem_real2 - Maximum Memory Limit Exceeded!',-1)     
  ENDIF
end subroutine mem_allocated_ichor_mem_real2

subroutine mem_deallocated_ichor_mem_real2(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  mem_allocated_real = mem_allocated_real - nsize
  mem_allocated_ichor = mem_allocated_ichor - nsize
end subroutine mem_deallocated_ichor_mem_real2

subroutine mem_allocated_ichor_mem_integer(nsize)
  implicit none
  integer(kind=long), intent(in) :: nsize
  mem_allocated_integer = mem_allocated_integer + nsize
  max_mem_used_integer = MAX(max_mem_used_integer,mem_allocated_integer)
  !Count also the total memory:
  mem_allocated_ichor = mem_allocated_ichor  + nsize
  max_mem_used_ichor = MAX(max_mem_used_ichor,mem_allocated_ichor)
  IF(max_mem_used_ichor.GT.maxMemLimit_ichor)THEN
     print*,'Maximum Memory Limit',maxMemLimit_ichor
     print*,'Memory Usage        ',max_mem_used_ichor
     call IchorQuit('Error in mem_allocated_ichor_mem_integer - Maximum Memory Limit Exceeded!',-1)     
  ENDIF
end subroutine mem_allocated_ichor_mem_integer

subroutine mem_deallocated_ichor_mem_integer(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  mem_allocated_integer = mem_allocated_integer - nsize
  if (mem_allocated_integer < 0) then
     call IchorQuit('Error in mem_deallocated_ichor_mem_integer1 - probably integer overflow!',-1)
  endif
  !Count also the total memory:
  mem_allocated_ichor = mem_allocated_ichor - nsize
  if (mem_allocated_ichor < 0) then
     call IchorQuit('Error in mem_deallocated_ichor_mem_integer2 - probably integer overflow!',-1)
  endif
end subroutine mem_deallocated_ichor_mem_integer

subroutine mem_allocated_ichor_mem_logical(nsize)
  implicit none
  integer(kind=long), intent(in) :: nsize
  mem_allocated_logical = mem_allocated_logical + nsize
  max_mem_used_logical = MAX(max_mem_used_logical,mem_allocated_logical)
  !Count also the total memory:
  mem_allocated_ichor = mem_allocated_ichor  + nsize
  max_mem_used_ichor = MAX(max_mem_used_ichor,mem_allocated_ichor)
  IF(max_mem_used_ichor.GT.maxMemLimit_ichor)THEN
     print*,'Maximum Memory Limit',maxMemLimit_ichor
     print*,'Memory Usage        ',max_mem_used_ichor
     call IchorQuit('Error in mem_allocated_ichor_mem_logical - Maximum Memory Limit Exceeded!',-1)     
  ENDIF
end subroutine mem_allocated_ichor_mem_logical

subroutine mem_deallocated_ichor_mem_logical(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  mem_allocated_logical = mem_allocated_logical - nsize
  if (mem_allocated_logical < 0) then
     call IchorQuit('Error in mem_deallocated_ichor_mem_logical1 - probably integer overflow!',-1)
  endif
  !Count also the total memory:
  mem_allocated_ichor = mem_allocated_ichor - nsize
  if (mem_allocated_ichor < 0) then
     call IchorQuit('Error in mem_deallocated_ichor_mem_logical2 - probably integer overflow!',-1)
  endif
end subroutine mem_deallocated_ichor_mem_logical

END MODULE IchorMemory

