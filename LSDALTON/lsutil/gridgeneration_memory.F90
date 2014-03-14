MODULE grid_memory_handling
   use precision
   use memory_handling
   private
   public mem_grid_alloc
   public mem_grid_dealloc
   public init_gridmemvar
   public mem_grid_TurnOffThread_Memory
   public mem_grid_TurnONThread_Memory
   public collect_thread_grid_memory
   public init_grid_threadmemvar
   public stats_grid_mem
   integer(KIND=long),save :: mem_allocated_grid, max_mem_used_grid
   logical,save :: mem_grid_InsideOMPsection
!THREADPRIVATE LOCAL VARIABLES 
integer(KIND=long),save :: mem_tp_allocated_grid, max_mem_tp_used_grid
!$OMP THREADPRIVATE(mem_tp_allocated_grid,max_mem_tp_used_grid)

!Interfaces for allocating/deallocating pointers
INTERFACE mem_grid_alloc
  MODULE PROCEDURE real_grid_allocate_1dim, real_grid_allocate_2dim, &
     &             real_grid_allocate_3dim, int_grid_allocate_1dim,  &
     &             int_grid_allocate_2dim,  int_grid_allocate_3dim, &
#ifndef VAR_INT64
     &             intlong_grid_allocate_1dim, & !in case of 64bit this is the same as int_allocate_1dim
#endif
     &             logic_grid_allocate_1dim, logic_grid_allocate_2dim
END INTERFACE
!
INTERFACE mem_grid_dealloc
  MODULE PROCEDURE real_grid_deallocate_1dim, real_grid_deallocate_2dim, &
     &             real_grid_deallocate_3dim, int_grid_deallocate_1dim, &
     &             int_grid_deallocate_2dim, int_grid_deallocate_3dim,  &
#ifndef VAR_INT64
     &             intlong_grid_deallocate_1dim, & !in case of 64bit this is the same as int_allocate_1dim
#endif
     &             logic_grid_deallocate_1dim, logic_grid_deallocate_2dim
END INTERFACE

CONTAINS
subroutine init_gridmemvar()
implicit none
mem_allocated_grid = 0
max_mem_used_grid = 0 
mem_grid_InsideOMPsection=.FALSE.
call init_grid_threadmemvar()
end subroutine init_gridmemvar

subroutine init_grid_threadmemvar()
implicit none
mem_tp_allocated_grid = 0
max_mem_tp_used_grid = 0    
end subroutine init_grid_threadmemvar

subroutine mem_grid_TurnOffThread_Memory()
implicit none
mem_grid_InsideOMPsection = .FALSE.
call mem_TurnOffThread_Memory()
end subroutine mem_grid_TurnOffThread_Memory

subroutine mem_grid_TurnONThread_Memory()
implicit none
mem_grid_InsideOMPsection = .TRUE.
call mem_TurnONThread_Memory()
end subroutine mem_grid_TurnONThread_Memory

subroutine collect_thread_grid_memory()
  implicit none
!$OMP CRITICAL
mem_allocated_grid = mem_allocated_grid+mem_tp_allocated_grid
max_mem_used_grid = max_mem_used_grid+max_mem_tp_used_grid
!$OMP END CRITICAL
call collect_thread_memory()
end subroutine collect_thread_grid_memory

!> \brief Print current and max. amount of memory allocated
!> \author T. Kjaergaard
!> \date 2011
subroutine stats_grid_mem(lupri)
  implicit none
  !> Logical unit number for output file.
  integer,intent(in) :: lupri
  WRITE(LUPRI,*)
  IF(mem_allocated_grid.GT. 0)THEN
     WRITE(LUPRI,'(" WARNING Allocated memory in Grid:",i9," byte  &
       &- Should be zero - A leakage is present")') mem_allocated_grid
  ENDIF
  call print_grid_maxmem(lupri,max_mem_used_grid,'Grid')
  WRITE(LUPRI,*)
end subroutine stats_grid_mem

!> \brief Print routine for memory statistics (max. allocated memory). 
!> \author S. Host 
!> \date 2009
subroutine print_grid_maxmem(lupri,max_mem_used,STRING)
  implicit none
  !> Logical unit number for output file
  integer :: lupri
  !> Amount of memory used
  integer(KIND=long) :: max_mem_used
  !> Name of data type for which memory usage is printed
  Character*(*)    ::  STRING
  character(len=20) :: trimstring
  trimstring = adjustl(string)
  if (max_mem_used < 100) then !Divide by 1 to typecast real
     WRITE(LUPRI,'("  Max allocated memory, ", a20, f9.3, " Byte")') trimstring, max_mem_used/1.0E0_realk
  else if (max_mem_used < 1000000) then
     WRITE(LUPRI,'("  Max allocated memory, ", a20, f9.3, " kB")') trimstring, max_mem_used/1024.0E0_realk
  else if (max_mem_used < 1000000000) then
     WRITE(LUPRI,'("  Max allocated memory, ", a20, f9.3, " MB")') trimstring, max_mem_used/(1024.0E0_realk*1024.0E0_realk)
  else
     WRITE(LUPRI,'("  Max allocated memory, ", a20, f9.3, " GB")') &
      & trimstring, max_mem_used/(1024.0E0_realk*1024.0E0_realk*1024.0E0_realk)
  endif
end subroutine print_grid_maxmem

SUBROUTINE real_grid_allocate_1dim(A,n)
implicit none
integer,intent(in)  :: n
REAL(REALK),pointer :: A(:)
integer (kind=long) :: nsize
call mem_alloc(A,n)
nsize = size(A)*mem_realsize
call mem_allocated_mem_grid(nsize)
END SUBROUTINE real_grid_allocate_1dim

SUBROUTINE real_grid_allocate_2dim(A,n1,n2)
implicit none
integer,intent(in)  :: n1, n2
REAL(REALK),pointer :: A(:,:)
integer (kind=long) :: nsize
call mem_alloc(A,n1,n2)
nsize = size(A)*mem_realsize
call mem_allocated_mem_grid(nsize)
END SUBROUTINE real_grid_allocate_2dim

SUBROUTINE real_grid_allocate_3dim(A,n1,n2,n3)
implicit none
integer,intent(in)  :: n1, n2, n3
REAL(REALK),pointer :: A(:,:,:)
integer (kind=long) :: nsize
call mem_alloc(A,n1,n2,n3)
nsize = size(A)*mem_realsize
call mem_allocated_mem_grid(nsize)
END SUBROUTINE real_grid_allocate_3dim

!----- DEALLOCATE REAL POINTERS -----!

SUBROUTINE real_grid_deallocate_1dim(A)
implicit none
REAL(REALK),pointer :: A(:)
integer (kind=long) :: nsize
nsize = size(A)*mem_realsize
call mem_deallocated_mem_grid(nsize)
call mem_dealloc(A)
END SUBROUTINE real_grid_deallocate_1dim

SUBROUTINE real_grid_deallocate_2dim(A)
implicit none
REAL(REALK),pointer :: A(:,:)
integer (kind=long) :: nsize
nsize = size(A)*mem_realsize
call mem_deallocated_mem_grid(nsize)
call mem_dealloc(A)
END SUBROUTINE real_grid_deallocate_2dim

SUBROUTINE real_grid_deallocate_3dim(A)
implicit none
REAL(REALK),pointer :: A(:,:,:)
integer (kind=long) :: nsize
nsize = size(A)*mem_realsize
call mem_deallocated_mem_grid(nsize)
call mem_dealloc(A)
END SUBROUTINE real_grid_deallocate_3dim

!----- ALLOCATE INTEGER POINTERS -----!

SUBROUTINE int_grid_allocate_1dim(I,n)
implicit none
integer,intent(in)  :: n
INTEGER,pointer     :: I(:)
integer (kind=long) :: nsize
call mem_alloc(I,n)
nsize = size(I)*mem_intsize
call mem_allocated_mem_grid(nsize)
END SUBROUTINE int_grid_allocate_1dim

#ifndef VAR_INT64
SUBROUTINE intlong_grid_allocate_1dim(I,n)
implicit none
integer,intent(in)  :: n
INTEGER(kind=long),pointer     :: I(:)
integer (kind=long) :: nsize
call mem_alloc(I,n)
nsize = size(I)*mem_intsize*2
call mem_allocated_mem_grid(nsize)
END SUBROUTINE intlong_grid_allocate_1dim
#endif

SUBROUTINE int_grid_allocate_2dim(I,n1,n2)
implicit none
integer,intent(in) :: n1,n2
INTEGER,pointer    :: I(:,:)
integer (kind=long) :: nsize
call mem_alloc(I,n1,n2)
nsize = size(I)*mem_intsize
call mem_allocated_mem_grid(nsize)
END SUBROUTINE int_grid_allocate_2dim

SUBROUTINE int_grid_allocate_3dim(I,n1,n2,n3)
implicit none
integer,intent(in) :: n1,n2,n3
INTEGER,pointer    :: I(:,:,:)
integer (kind=long) :: nsize
call mem_alloc(I,n1,n2,n3)
nsize = size(I)*mem_intsize
call mem_allocated_mem_grid(nsize)
END SUBROUTINE int_grid_allocate_3dim

!----- DEALLOCATE INTEGER POINTERS -----!

SUBROUTINE int_grid_deallocate_1dim(I)
implicit none
INTEGER,pointer :: I(:)
integer (kind=long) :: nsize
nsize = size(I)*mem_intsize
call mem_deallocated_mem_grid(nsize)
call mem_dealloc(I)
END SUBROUTINE int_grid_deallocate_1dim

#ifndef VAR_INT64
SUBROUTINE intlong_grid_deallocate_1dim(I)
implicit none
INTEGER(kind=long),pointer :: I(:)
integer (kind=long) :: nsize
nsize = size(I)*mem_intsize*2
call mem_deallocated_mem_grid(nsize)
call mem_dealloc(I)
 END SUBROUTINE intlong_grid_deallocate_1dim
#endif

SUBROUTINE int_grid_deallocate_2dim(I)
implicit none
INTEGER,pointer :: I(:,:)
integer (kind=long) :: nsize
nsize = size(I)*mem_intsize
call mem_deallocated_mem_grid(nsize)
call mem_dealloc(I)
END SUBROUTINE int_grid_deallocate_2dim

SUBROUTINE int_grid_deallocate_3dim(I)
implicit none
INTEGER,pointer :: I(:,:,:)
integer (kind=long) :: nsize
nsize = size(I)*mem_intsize
call mem_deallocated_mem_grid(nsize)
call mem_dealloc(I)
END SUBROUTINE int_grid_deallocate_3dim

!----- ALLOCATE LOGICAL POINTERS -----!

SUBROUTINE logic_grid_allocate_1dim(L,n)
implicit none
integer,intent(in) :: n
LOGICAL,pointer    :: L(:)
integer (kind=long) :: nsize
call mem_alloc(L,n)
nsize = size(L)*mem_logicalsize 
call mem_allocated_mem_grid(nsize)
END SUBROUTINE logic_grid_allocate_1dim

SUBROUTINE logic_grid_allocate_2dim(L,n1,n2)
implicit none
integer,intent(in) :: n1,n2
LOGICAL,pointer    :: L(:,:)
integer (kind=long) :: nsize
call mem_alloc(L,n1,n2)
nsize = size(L)*mem_logicalsize 
call mem_allocated_mem_grid(nsize)
END SUBROUTINE logic_grid_allocate_2dim

!----- DEALLOCATE LOGICAL POINTERS -----!

SUBROUTINE logic_grid_deallocate_1dim(L)
implicit none
LOGICAL,pointer :: L(:)
integer :: IERR
integer (kind=long) :: nsize
nsize = size(L)*mem_logicalsize
call mem_deallocated_mem_grid(nsize)
call mem_dealloc(L)
END SUBROUTINE logic_grid_deallocate_1dim

SUBROUTINE logic_grid_deallocate_2dim(L)
implicit none
LOGICAL,pointer :: L(:,:)
integer :: IERR
integer (kind=long) :: nsize
nsize = size(L)*mem_logicalsize
call mem_deallocated_mem_grid(nsize)
call mem_dealloc(L)
END SUBROUTINE logic_grid_deallocate_2dim

!----- MEMORY HANDLING -----!

subroutine mem_allocated_mem_grid(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  IF(mem_grid_InsideOMPsection)THEN!we add to thread private variables
     mem_tp_allocated_grid = mem_tp_allocated_grid + nsize
     if (mem_tp_allocated_grid < 0) then
        write(*,*) 'Real memory negative! mem_tp_allocated_grid =', mem_tp_allocated_grid
        write(*,*) 'Real memory negative! nsize =', nsize
        call lsQUIT('Error in mem_tp_allocated_mem_grid-probably integer overflow!',-1)
     endif
     max_mem_tp_used_grid = MAX(max_mem_tp_used_grid,mem_tp_allocated_grid)
  ELSE
     mem_allocated_grid = mem_allocated_grid + nsize
     if (mem_allocated_grid < 0) then
        write(*,*) 'Real memory negative! mem_allocated_grid =', mem_allocated_grid
        write(*,*) 'Real memory negative! nsize =', nsize
        call lsQUIT('Error in mem_allocated_mem_grid - probably integer overflow!',-1)
     endif
     max_mem_used_grid = MAX(max_mem_used_grid,mem_allocated_grid)
  ENDIF
end subroutine mem_allocated_mem_grid

subroutine mem_deallocated_mem_grid(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  IF(mem_grid_InsideOMPsection)THEN!we add to thread private variables
     mem_tp_allocated_grid = mem_tp_allocated_grid - nsize
     if (mem_tp_allocated_grid < 0) then
        write(*,*) 'Real memory negative! mem_tp_allocated_grid =', mem_tp_allocated_grid
        call lsQUIT('Error in mem_tp_deallocated_mem_grid - something wrong with deallocation!',-1)
     endif
  else
     mem_allocated_grid = mem_allocated_grid - nsize
     if (mem_allocated_grid < 0) then
        write(*,*) 'Real memory negative! mem_allocated_grid =', mem_allocated_grid
        call lsQUIT('Error in mem_deallocated_mem_grid - something wrong with deallocation!',-1)
     endif
  endif
end subroutine mem_deallocated_mem_grid

END MODULE grid_memory_handling

