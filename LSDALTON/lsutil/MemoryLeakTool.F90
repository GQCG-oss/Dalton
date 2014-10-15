!Monitor memory NON THREAD SAFE
MODULE MemoryLeakToolMod
  use precision
  use memory_handling
!GLOBAL VARIABLES
   integer(KIND=long),save :: mem_allocated_LeakTool, max_mem_used_LeakTool       !Count all memory
   integer(KIND=long),save :: mem_allocated_real, max_mem_used_real       !Count 'real' memory   
   integer(KIND=long),save :: mem_allocated_integer, max_mem_used_integer !Count 'integer' memory
   integer(KIND=long),save :: mem_allocated_logical, max_mem_used_logical !Count 'logical' memory
   integer(KIND=long),save :: maxMemLimit_LeakTool !maximum allowed memory usage
!sizes
   integer,parameter :: nAllocTags = 10
   integer,parameter :: LT_RIMP2_AlphaBeta_inv = 1
   integer,parameter :: LT_RIMP2_TMPAlphaBeta_inv = 2
   integer,parameter :: LT_TMP = 3
   integer,parameter :: LT_Calpha = 4
   integer,parameter :: LT_AlphaCD5 = 5
   integer,parameter :: LT_AlphaCD6 = 6
   integer,parameter :: LT_Eps = 7 
   integer,parameter :: LT_AlphaCD2 = 8
   integer,parameter :: LT_Calpha2 = 9 
   integer,parameter :: LT_AlphaCD = 10 
   integer :: nAlloced(nAllocTags)
   

INTERFACE mem_LeakTool_alloc
  MODULE PROCEDURE real_allocate_1dim, real_allocate_2dim, &
       & real_allocate_3dim, real_allocate_4dim, real_allocate_5dim, &
       & int_allocate_1dim, int_allocate_2dim, int_allocate_3dim, &
       & logic_allocate_1dim, logic_allocate_2dim

END INTERFACE
!
INTERFACE mem_LeakTool_dealloc
  MODULE PROCEDURE real_deallocate_1dim, real_deallocate_2dim,  &
     & real_deallocate_3dim, real_deallocate_4dim, real_deallocate_5dim, &
     & int_deallocate_1dim, int_deallocate_2dim, int_deallocate_3dim,&
     & logic_deallocate_1dim, logic_deallocate_2dim
END INTERFACE
!
CONTAINS
subroutine set_LeakTool_memvar(MaxMemAllocated,MemAllocated,MaxMem)
implicit none
integer(KIND=long),intent(in) :: MaxMemAllocated,MemAllocated
integer(KIND=long),optional :: MaxMem
Integer :: I
max_mem_used_LeakTool = MaxMemAllocated
mem_allocated_LeakTool = MemAllocated
IF(present(MaxMem))THEN
   maxMemLimit_LeakTool = MaxMem
ELSE
   maxMemLimit_LeakTool = HUGE(MaxMemAllocated)
ENDIF
mem_allocated_real = 0
max_mem_used_real = 0
mem_allocated_integer = 0
max_mem_used_integer = 0
mem_allocated_logical = 0
max_mem_used_logical = 0
DO I = 1,nAllocTags
   nAlloced(I) = 0 
ENDDO
end subroutine set_LeakTool_memvar

subroutine retrieve_LeakTool_memvar(MaxMemAllocated,MemAllocated)
implicit none
integer(KIND=long),intent(inout) :: MaxMemAllocated,MemAllocated
MemAllocated = mem_allocated_LeakTool
MaxMemAllocated = max_mem_used_LeakTool
end subroutine retrieve_LeakTool_memvar

subroutine LeakTools_stat_mem(lupri)
  implicit none
  !> Logical unit number for output file.
  integer,intent(in) :: lupri
  !
  integer :: I
  logical :: FoundLeak
  WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
  WRITE(LUPRI,'("            LeakTool Memory statistics          ")')
  WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
  WRITE(LUPRI,'("  Allocated memory (TOTAL):           ",i9," byte  &
       &- Should be zero - otherwise a leakage is present")') mem_allocated_LeakTool
  WRITE(LUPRI,'("  Allocated memory (real):            ",i9," byte  &
       &- Should be zero - otherwise a leakage is present")') mem_allocated_real
  WRITE(LUPRI,'("  Allocated memory (integer):         ",i9," byte  &
       &- Should be zero - otherwise a leakage is present")') mem_allocated_integer
  WRITE(LUPRI,'("  Allocated memory (logical):         ",i9," byte  &
       &- Should be zero - otherwise a leakage is present")') mem_allocated_logical
  call print_LeakTool_maxmem(lupri,max_mem_used_LeakTool,'TOTAL')
  CALL print_LeakTool_maxmem(lupri,max_mem_used_real,'real(realk)')
  CALL print_LeakTool_maxmem(lupri,max_mem_used_integer,'integer')
  CALL print_LeakTool_maxmem(lupri,max_mem_used_logical,'logical')
  WRITE(LUPRI,'("*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*")')
  WRITE(LUPRI,*)
  FoundLeak = .FALSE.
  DO I = 1,nAllocTags
     IF(nAlloced(I).NE.0)THEN
        WRITE(LUPRI,*)'Memory leak for Identifyier: ',I        
        FoundLeak = .TRUE.
     ENDIF
  ENDDO
  IF(FoundLeak)THEN
     WRITE(lupri,*)'Find Identifyier in the list in statmemory.F90'
     WRITE(lupri,*)'integer,parameter :: LT_RIMP2_AlphaBeta_inv = 1'
     WRITE(lupri,*)'integer,parameter :: LT_RIMP2_TMPAlphaBeta_inv = 2'
     WRITE(lupri,*)'integer,parameter :: LT_TMP = 3'
     WRITE(lupri,*)' ... '
  ENDIF
end subroutine LeakTools_stat_mem

!> \brief Print routine for memory LeakToolistics (max. allocated memory).
!> \author S. Host
!> \date 2009
subroutine print_LeakTool_maxmem(lupri,max_mem_used,STRING)
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
end subroutine print_LeakTool_maxmem

SUBROUTINE real_allocate_1dim(A,II)
  implicit none
  REAL(REALK) :: A(:)
  integer,optional :: II
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_LeakTool_mem_real(nsize)
  IF(present(II)) nAlloced(II) = nAlloced(II) + 1
END SUBROUTINE real_allocate_1dim

SUBROUTINE real_allocate_2dim(A,II)
  implicit none
  REAL(REALK) :: A(:,:)
  integer,optional :: II
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_LeakTool_mem_real(nsize)
  IF(present(II)) nAlloced(II) = nAlloced(II) + 1
END SUBROUTINE real_allocate_2dim

SUBROUTINE ADD_OMP_MEM(nsize1)
  implicit none
  integer :: nsize1
  integer (kind=long) :: nsize
  nsize = nsize1
  call mem_allocated_LeakTool_mem_real2(nsize)
END SUBROUTINE ADD_OMP_MEM

SUBROUTINE REMOVE_OMP_MEM(nsize1)
  implicit none
  integer :: nsize1
  integer (kind=long) :: nsize
  nsize = nsize1
  call mem_deallocated_LeakTool_mem_real2(nsize)
END SUBROUTINE REMOVE_OMP_MEM

SUBROUTINE real_allocate_3dim(A,II)
  implicit none
  REAL(REALK) :: A(:,:,:)
  integer,optional :: II
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_LeakTool_mem_real(nsize)
  IF(present(II)) nAlloced(II) = nAlloced(II) + 1
END SUBROUTINE real_allocate_3dim

SUBROUTINE real_allocate_4dim(A,II)
  implicit none
  REAL(REALK) :: A(:,:,:,:)
  integer,optional :: II
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_LeakTool_mem_real(nsize)
  IF(present(II)) nAlloced(II) = nAlloced(II) + 1
END SUBROUTINE real_allocate_4dim

SUBROUTINE real_allocate_5dim(A,II)
  implicit none
  REAL(REALK) :: A(:,:,:,:,:)
  integer,optional :: II
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_LeakTool_mem_real(nsize)
  IF(present(II)) nAlloced(II) = nAlloced(II) + 1
END SUBROUTINE real_allocate_5dim

SUBROUTINE real_deallocate_1dim(A,II)
  implicit none
  REAL(REALK) :: A(:)
  integer,optional :: II
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_deallocated_LeakTool_mem_real(nsize)
  IF(present(II)) nAlloced(II) = nAlloced(II) - 1
END SUBROUTINE real_deallocate_1dim

SUBROUTINE real_deallocate_2dim(A,II)
implicit none
REAL(REALK) :: A(:,:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(A,KIND=long)*mem_realsize
call mem_deallocated_LeakTool_mem_real(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) - 1
END SUBROUTINE real_deallocate_2dim

SUBROUTINE real_deallocate_3dim(A,II)
implicit none
REAL(REALK) :: A(:,:,:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(A,KIND=long)*mem_realsize
call mem_deallocated_LeakTool_mem_real(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) - 1
END SUBROUTINE real_deallocate_3dim

SUBROUTINE real_deallocate_4dim(A,II)
implicit none
REAL(REALK) :: A(:,:,:,:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(A,KIND=long)*mem_realsize
call mem_deallocated_LeakTool_mem_real(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) - 1
END SUBROUTINE real_deallocate_4dim

SUBROUTINE real_deallocate_5dim(A,II)
implicit none
REAL(REALK) :: A(:,:,:,:,:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(A,KIND=long)*mem_realsize
call mem_deallocated_LeakTool_mem_real(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) - 1
END SUBROUTINE real_deallocate_5dim

SUBROUTINE int_allocate_1dim(I,II)
implicit none
INTEGER             :: I(:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_LeakTool_mem_integer(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) + 1
END SUBROUTINE int_allocate_1dim

SUBROUTINE int_allocate_2dim(I,II)
implicit none
INTEGER    :: I(:,:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_LeakTool_mem_integer(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) + 1
END SUBROUTINE int_allocate_2dim

SUBROUTINE int_allocate_3dim(I,II)
implicit none
INTEGER :: I(:,:,:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_LeakTool_mem_integer(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) + 1
END SUBROUTINE int_allocate_3dim

SUBROUTINE int_allocate_4dim(I,II)
implicit none
INTEGER :: I(:,:,:,:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_LeakTool_mem_integer(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) + 1
END SUBROUTINE int_allocate_4dim

SUBROUTINE int_deallocate_1dim(I,II)
implicit none
INTEGER :: I(:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_deallocated_LeakTool_mem_integer(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) - 1
END SUBROUTINE int_deallocate_1dim

SUBROUTINE int_deallocate_2dim(I,II)
implicit none
INTEGER :: I(:,:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_deallocated_LeakTool_mem_integer(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) - 1
END SUBROUTINE int_deallocate_2dim

SUBROUTINE int_deallocate_3dim(I,II)
implicit none
INTEGER :: I(:,:,:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_deallocated_LeakTool_mem_integer(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) - 1
END SUBROUTINE int_deallocate_3dim

SUBROUTINE logic_allocate_1dim(L,II)
implicit none
logical     :: L(:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(L,KIND=long)*mem_logicalsize
call mem_allocated_LeakTool_mem_logical(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) + 1
END SUBROUTINE logic_allocate_1dim

SUBROUTINE logic_allocate_2dim(L,II)
implicit none
logical    :: L(:,:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(L,KIND=long)*mem_logicalsize
call mem_allocated_LeakTool_mem_logical(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) + 1
END SUBROUTINE logic_allocate_2dim

SUBROUTINE logic_deallocate_1dim(L,II)
implicit none
logical :: L(:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(L,KIND=long)*mem_logicalsize
call mem_deallocated_LeakTool_mem_logical(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) - 1
END SUBROUTINE logic_deallocate_1dim

SUBROUTINE logic_deallocate_2dim(L,II)
implicit none
LOGICAL :: L(:,:)
integer,optional :: II
integer (kind=long) :: nsize
nsize = size(L,KIND=long)*mem_logicalsize
call mem_deallocated_LeakTool_mem_logical(nsize)
IF(present(II)) nAlloced(II) = nAlloced(II) - 1
END SUBROUTINE logic_deallocate_2dim

!----- MEMORY HANDLING -----!

subroutine mem_allocated_LeakTool_mem_real(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  mem_allocated_real = mem_allocated_real + nsize
  if (mem_allocated_real < 0) then
     write(*,*) 'Real memory negative! mem_allocated_real =', mem_allocated_real
     write(*,*) 'Real memory negative! nsize =', nsize
     call Lsquit('Error in mem_allocated_LeakTool_mem_real - probably integer overflow!',-1)
  endif
  max_mem_used_real = MAX(max_mem_used_real,mem_allocated_real)
  !Count also the total memory:
  mem_allocated_LeakTool = mem_allocated_LeakTool  + nsize
  if (mem_allocated_LeakTool < 0) then
     write(*,*) 'Total memory negative! mem_allocated_LeakTool =', mem_allocated_LeakTool
     call Lsquit('Error in mem_allocated_LeakTool_mem_real - probably integer overflow!',-1)
  endif
  max_mem_used_LeakTool = MAX(max_mem_used_LeakTool,mem_allocated_LeakTool)
  IF(max_mem_used_LeakTool.GT.maxMemLimit_LeakTool)THEN
     print*,'Maximum Memory Limit',maxMemLimit_LeakTool
     print*,'Memory Usage        ',max_mem_used_LeakTool
     call Lsquit('Error in mem_allocated_LeakTool_mem_real - Maximum Memory Limit Exceeded!',-1)     
  ENDIF
end subroutine mem_allocated_LeakTool_mem_real

subroutine mem_deallocated_LeakTool_mem_real(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  
  mem_allocated_real = mem_allocated_real - nsize
  if (mem_allocated_real < 0) then
     write(*,*) 'Real memory negative! mem_allocated_real =', mem_allocated_real
     call Lsquit('Error in mem_deallocated_LeakTool_mem_real - something wrong with deallocation!',-1)
  endif
  !Count also the total memory:
  mem_allocated_LeakTool = mem_allocated_LeakTool - nsize
  if (mem_allocated_LeakTool < 0) then
     write(*,*) 'Total memory negative! mem_allocated_LeakTool =', mem_allocated_LeakTool
     call Lsquit('Error in mem_deallocated_LeakTool_mem_real - something wrong with deallocation!',-1)
  endif
end subroutine mem_deallocated_LeakTool_mem_real

subroutine mem_allocated_LeakTool_mem_real2(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  mem_allocated_real = mem_allocated_real + nsize
  max_mem_used_real = MAX(max_mem_used_real,mem_allocated_real)
  mem_allocated_LeakTool = mem_allocated_LeakTool  + nsize
  max_mem_used_LeakTool = MAX(max_mem_used_LeakTool,mem_allocated_LeakTool)
  IF(max_mem_used_LeakTool.GT.maxMemLimit_LeakTool)THEN
     print*,'Maximum Memory Limit',maxMemLimit_LeakTool
     print*,'Memory Usage        ',max_mem_used_LeakTool
     call Lsquit('Error in mem_allocated_LeakTool_mem_real2 - Maximum Memory Limit Exceeded!',-1)     
  ENDIF
end subroutine mem_allocated_LeakTool_mem_real2

subroutine mem_deallocated_LeakTool_mem_real2(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  mem_allocated_real = mem_allocated_real - nsize
  mem_allocated_LeakTool = mem_allocated_LeakTool - nsize
end subroutine mem_deallocated_LeakTool_mem_real2

subroutine mem_allocated_LeakTool_mem_integer(nsize)
  implicit none
  integer(kind=long), intent(in) :: nsize
  mem_allocated_integer = mem_allocated_integer + nsize
  max_mem_used_integer = MAX(max_mem_used_integer,mem_allocated_integer)
  !Count also the total memory:
  mem_allocated_LeakTool = mem_allocated_LeakTool  + nsize
  max_mem_used_LeakTool = MAX(max_mem_used_LeakTool,mem_allocated_LeakTool)
  IF(max_mem_used_LeakTool.GT.maxMemLimit_LeakTool)THEN
     print*,'Maximum Memory Limit',maxMemLimit_LeakTool
     print*,'Memory Usage        ',max_mem_used_LeakTool
     call Lsquit('Error in mem_allocated_LeakTool_mem_integer - Maximum Memory Limit Exceeded!',-1)     
  ENDIF
end subroutine mem_allocated_LeakTool_mem_integer

subroutine mem_deallocated_LeakTool_mem_integer(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  mem_allocated_integer = mem_allocated_integer - nsize
  if (mem_allocated_integer < 0) then
     call Lsquit('Error in mem_deallocated_LeakTool_mem_integer1 - probably integer overflow!',-1)
  endif
  !Count also the total memory:
  mem_allocated_LeakTool = mem_allocated_LeakTool - nsize
  if (mem_allocated_LeakTool < 0) then
     call Lsquit('Error in mem_deallocated_LeakTool_mem_integer2 - probably integer overflow!',-1)
  endif
end subroutine mem_deallocated_LeakTool_mem_integer

subroutine mem_allocated_LeakTool_mem_logical(nsize)
  implicit none
  integer(kind=long), intent(in) :: nsize
  mem_allocated_logical = mem_allocated_logical + nsize
  max_mem_used_logical = MAX(max_mem_used_logical,mem_allocated_logical)
  !Count also the total memory:
  mem_allocated_LeakTool = mem_allocated_LeakTool  + nsize
  max_mem_used_LeakTool = MAX(max_mem_used_LeakTool,mem_allocated_LeakTool)
  IF(max_mem_used_LeakTool.GT.maxMemLimit_LeakTool)THEN
     print*,'Maximum Memory Limit',maxMemLimit_LeakTool
     print*,'Memory Usage        ',max_mem_used_LeakTool
     call Lsquit('Error in mem_allocated_LeakTool_mem_logical - Maximum Memory Limit Exceeded!',-1)     
  ENDIF
end subroutine mem_allocated_LeakTool_mem_logical

subroutine mem_deallocated_LeakTool_mem_logical(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  mem_allocated_logical = mem_allocated_logical - nsize
  if (mem_allocated_logical < 0) then
     call Lsquit('Error in mem_deallocated_LeakTool_mem_logical1 - probably integer overflow!',-1)
  endif
  !Count also the total memory:
  mem_allocated_LeakTool = mem_allocated_LeakTool - nsize
  if (mem_allocated_LeakTool < 0) then
     call Lsquit('Error in mem_deallocated_LeakTool_mem_logical2 - probably integer overflow!',-1)
  endif
end subroutine mem_deallocated_LeakTool_mem_logical

END MODULE MemoryLeakToolMod

