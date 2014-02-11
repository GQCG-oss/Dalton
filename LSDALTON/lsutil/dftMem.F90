MODULE dft_memory_handling
   use precision
   use memory_handling
   private
   public mem_dft_alloc
   public mem_dft_dealloc
   public init_dftmemvar
   public mem_dft_TurnOffThread_Memory
   public mem_dft_TurnONThread_Memory
   public collect_thread_dft_memory
   public init_dft_threadmemvar
   public stats_dft_mem
   public stats_dft_mem_debug
   integer(KIND=long),save :: mem_allocated_dft, max_mem_used_dft
   logical,save :: mem_dft_InsideOMPsection
!THREADPRIVATE LOCAL VARIABLES 
integer(KIND=long),save :: mem_tp_allocated_dft, max_mem_tp_used_dft
!$OMP THREADPRIVATE(mem_tp_allocated_dft,max_mem_tp_used_dft)

!Interfaces for allocating/deallocating pointers
INTERFACE mem_dft_alloc
  MODULE PROCEDURE real_dft_allocate_1dim, real_dft_allocate_2dim, &
     &             real_dft_allocate_3dim, real_dft_allocate_4dim, &
     &             real_dft_allocate_5dim, int_dft_allocate_1dim,  &
     &             int_dft_allocate_2dim,  int_dft_allocate_3dim, &
#ifndef VAR_INT64
     &             intlong_dft_allocate_1dim, & !in case of 64bit this is the same as int_allocate_1dim
#endif
     &             char_dft_allocate_1dim, &
     &             logic_dft_allocate_1dim, logic_dft_allocate_2dim
END INTERFACE
!
INTERFACE mem_dft_dealloc
  MODULE PROCEDURE real_dft_deallocate_1dim, real_dft_deallocate_2dim, &
     &             real_dft_deallocate_3dim, real_dft_deallocate_4dim, &
     &             real_dft_deallocate_5dim, int_dft_deallocate_1dim, &
     &             int_dft_deallocate_2dim, int_dft_deallocate_3dim,  &
#ifndef VAR_INT64
     &             intlong_dft_deallocate_1dim, & !in case of 64bit this is the same as int_allocate_1dim
#endif
     &             char_dft_deallocate_1dim, &
     &             logic_dft_deallocate_1dim, logic_dft_deallocate_2dim
END INTERFACE

CONTAINS
subroutine init_dftmemvar()
implicit none
mem_allocated_dft = 0
max_mem_used_dft = 0 
call set_mem_XCcalc(mem_allocated_dft,max_mem_used_dft)
mem_dft_InsideOMPsection=.FALSE.
call init_dft_threadmemvar()
end subroutine init_dftmemvar

subroutine init_dft_threadmemvar()
implicit none
mem_tp_allocated_dft = 0
max_mem_tp_used_dft = 0
call init_threadmemvar()
end subroutine init_dft_threadmemvar

subroutine mem_dft_TurnOffThread_Memory()
implicit none
mem_dft_InsideOMPsection = .FALSE.
call mem_TurnOffThread_Memory()
end subroutine mem_dft_TurnOffThread_Memory

subroutine mem_dft_TurnONThread_Memory()
implicit none
mem_dft_InsideOMPsection = .TRUE.
call mem_TurnONThread_Memory()
end subroutine mem_dft_TurnONThread_Memory

subroutine collect_thread_dft_memory()
  implicit none
!$OMP CRITICAL
mem_allocated_dft = mem_allocated_dft+mem_tp_allocated_dft
max_mem_used_dft = max_mem_used_dft+max_mem_tp_used_dft
!$OMP END CRITICAL
call collect_thread_memory()
end subroutine collect_thread_dft_memory

!> \brief Print current and max. amount of memory allocated
!> \author T. Kjaergaard
!> \date 2011
subroutine stats_dft_mem(lupri)
  implicit none
  !> Logical unit number for output file.
  integer,intent(in) :: lupri
  IF(mem_allocated_dft.GT. 0)THEN
     WRITE(LUPRI,'(" WARNING Allocated memory in XC:",i9," byte  &
       &- Should be zero - A leakage is present")') mem_allocated_dft
  ENDIF
  call print_dft_maxmem(lupri,max_mem_used_dft,'XC')
  call set_mem_XCcalc(mem_allocated_dft,max_mem_used_dft)
end subroutine stats_dft_mem

!> \brief Print current and max. amount of memory allocated, for debuggin
!> \author T. Kjaergaard
!> \date 2011
subroutine stats_dft_mem_debug(lupri)
  implicit none
  !> Logical unit number for output file.
  integer,intent(in) :: lupri

  if (mem_allocated_dft < 100) then !Divide by 1 to typecast real
     WRITE(LUPRI,'("  Memory allocated ", f9.3, " Byte")') mem_allocated_dft/1.0E0_realk
  else if (mem_allocated_dft < 1000000) then
     WRITE(LUPRI,'("  Memory allocated ", f9.3, " kB")') mem_allocated_dft/1000.0E0_realk
  else if (mem_allocated_dft < 1000000000) then
     WRITE(LUPRI,'("  Memory allocated ", f9.3, " MB")') mem_allocated_dft/(1000.0E0_realk*1000.0E0_realk)
  else
     WRITE(LUPRI,'("  Memory allocated ", f9.3, " GB")') &
      & mem_allocated_dft/(1000.0E0_realk*1000.0E0_realk*1000.0E0_realk)
  endif

  if (max_mem_used_dft < 100) then !Divide by 1 to typecast real
     WRITE(LUPRI,'("  Max Memory used  ", f9.3, " Byte")') max_mem_used_dft/1.0E0_realk
  else if (max_mem_used_dft < 1000000) then
     WRITE(LUPRI,'("  Max Memory used  ", f9.3, " kB")') max_mem_used_dft/1000.0E0_realk
  else if (max_mem_used_dft < 1000000000) then
     WRITE(LUPRI,'("  Max Memory used  ", f9.3, " MB")') max_mem_used_dft/(1000.0E0_realk*1000.0E0_realk)
  else
     WRITE(LUPRI,'("  Max Memory used  ", f9.3, " GB")') &
      & max_mem_used_dft/(1000.0E0_realk*1000.0E0_realk*1000.0E0_realk)
  endif

  if (mem_tp_allocated_dft < 100) then !Divide by 1 to typecast real
     WRITE(LUPRI,'("  TP Memory allocated ", f9.3, " Byte")') mem_tp_allocated_dft/1.0E0_realk
  else if (mem_tp_allocated_dft < 1000000) then
     WRITE(LUPRI,'("  TP Memory allocated ", f9.3, " kB")') mem_tp_allocated_dft/1000.0E0_realk
  else if (mem_tp_allocated_dft < 1000000000) then
     WRITE(LUPRI,'("  TP Memory allocated ", f9.3, " MB")') mem_tp_allocated_dft/(1000.0E0_realk*1000.0E0_realk)
  else
     WRITE(LUPRI,'("  TP Memory allocated ", f9.3, " GB")') &
      & mem_tp_allocated_dft/(1000.0E0_realk*1000.0E0_realk*1000.0E0_realk)
  endif

  if (max_mem_tp_used_dft < 100) then !Divide by 1 to typecast real
     WRITE(LUPRI,'("  Max TP Memory used  ", f9.3, " Byte")') max_mem_tp_used_dft/1.0E0_realk
  else if (max_mem_tp_used_dft < 1000000) then
     WRITE(LUPRI,'("  Max TP Memory used  ", f9.3, " kB")') max_mem_tp_used_dft/1000.0E0_realk
  else if (max_mem_tp_used_dft < 1000000000) then
     WRITE(LUPRI,'("  Max TP Memory used  ", f9.3, " MB")') max_mem_tp_used_dft/(1000.0E0_realk*1000.0E0_realk)
  else
     WRITE(LUPRI,'("  Max TP Memory used  ", f9.3, " GB")') &
      & max_mem_tp_used_dft/(1000.0E0_realk*1000.0E0_realk*1000.0E0_realk)
  endif

end subroutine stats_dft_mem_debug

!> \brief Print routine for memory statistics (max. allocated memory). 
!> \author S. Host 
!> \date 2009
subroutine print_dft_maxmem(lupri,max_mem_used,STRING)
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
     WRITE(LUPRI,'("  Max allocated memory, ", a20, f9.3, " kB")') trimstring, max_mem_used/1000.0E0_realk
  else if (max_mem_used < 1000000000) then
     WRITE(LUPRI,'("  Max allocated memory, ", a20, f9.3, " MB")') trimstring, max_mem_used/(1000.0E0_realk*1000.0E0_realk)
  else
     WRITE(LUPRI,'("  Max allocated memory, ", a20, f9.3, " GB")') &
      & trimstring, max_mem_used/(1000.0E0_realk*1000.0E0_realk*1000.0E0_realk)
  endif
end subroutine print_dft_maxmem

SUBROUTINE real_dft_allocate_1dim(A,n)
implicit none
integer,intent(in)  :: n
REAL(REALK),pointer :: A(:)
integer (kind=long) :: nsize
call mem_alloc(A,n)
nsize = size(A)*mem_realsize
call mem_allocated_mem_dft(nsize)
END SUBROUTINE real_dft_allocate_1dim

SUBROUTINE real_dft_allocate_2dim(A,n1,n2)
implicit none
integer,intent(in)  :: n1, n2
REAL(REALK),pointer :: A(:,:)
integer (kind=long) :: nsize
call mem_alloc(A,n1,n2)
nsize = size(A)*mem_realsize
call mem_allocated_mem_dft(nsize)
END SUBROUTINE real_dft_allocate_2dim

SUBROUTINE real_dft_allocate_3dim(A,n1,n2,n3)
implicit none
integer,intent(in)  :: n1, n2, n3
REAL(REALK),pointer :: A(:,:,:)
integer (kind=long) :: nsize
call mem_alloc(A,n1,n2,n3)
nsize = size(A)*mem_realsize
call mem_allocated_mem_dft(nsize)
END SUBROUTINE real_dft_allocate_3dim

SUBROUTINE real_dft_allocate_4dim(A,n1,n2,n3,n4)
implicit none
integer,intent(in)  :: n1,n2,n3,n4
REAL(REALK),pointer :: A(:,:,:,:)
integer (kind=long) :: nsize
call mem_alloc(A,n1,n2,n3,n4)
nsize = size(A)*mem_realsize
call mem_allocated_mem_dft(nsize)
END SUBROUTINE real_dft_allocate_4dim

SUBROUTINE real_dft_allocate_5dim(A,n1,n2,n3,n4,n5)
implicit none
integer,intent(in)  :: n1,n2,n3,n4,n5
REAL(REALK),pointer :: A(:,:,:,:,:)
integer (kind=long) :: nsize
call mem_alloc(A,n1,n2,n3,n4,n5)
nsize = size(A)*mem_realsize
call mem_allocated_mem_real(nsize)
END SUBROUTINE real_dft_allocate_5dim

!----- DEALLOCATE REAL POINTERS -----!

SUBROUTINE real_dft_deallocate_1dim(A)
implicit none
REAL(REALK),pointer :: A(:)
integer (kind=long) :: nsize
nsize = size(A)*mem_realsize
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(A)
END SUBROUTINE real_dft_deallocate_1dim

SUBROUTINE real_dft_deallocate_2dim(A)
implicit none
REAL(REALK),pointer :: A(:,:)
integer (kind=long) :: nsize
nsize = size(A)*mem_realsize
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(A)
END SUBROUTINE real_dft_deallocate_2dim

SUBROUTINE real_dft_deallocate_3dim(A)
implicit none
REAL(REALK),pointer :: A(:,:,:)
integer (kind=long) :: nsize
nsize = size(A)*mem_realsize
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(A)
END SUBROUTINE real_dft_deallocate_3dim

SUBROUTINE real_dft_deallocate_4dim(A)
implicit none
REAL(REALK),pointer :: A(:,:,:,:)
integer (kind=long) :: nsize
nsize = size(A)*mem_realsize
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(A)
END SUBROUTINE real_dft_deallocate_4dim

SUBROUTINE real_dft_deallocate_5dim(A)
implicit none
REAL(REALK),pointer :: A(:,:,:,:,:)
integer (kind=long) :: nsize
nsize = size(A)*mem_realsize
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(A)
END SUBROUTINE real_dft_deallocate_5dim

!----- ALLOCATE INTEGER POINTERS -----!

SUBROUTINE int_dft_allocate_1dim(I,n)
implicit none
integer,intent(in)  :: n
INTEGER,pointer     :: I(:)
integer (kind=long) :: nsize
call mem_alloc(I,n)
nsize = size(I)*mem_intsize
call mem_allocated_mem_dft(nsize)
END SUBROUTINE int_dft_allocate_1dim

#ifndef VAR_INT64
SUBROUTINE intlong_dft_allocate_1dim(I,n)
implicit none
integer,intent(in)  :: n
INTEGER(kind=long),pointer     :: I(:)
integer (kind=long) :: nsize
call mem_alloc(I,n)
nsize = size(I)*mem_intsize*2
call mem_allocated_mem_dft(nsize)
END SUBROUTINE intlong_dft_allocate_1dim
#endif

SUBROUTINE int_dft_allocate_2dim(I,n1,n2)
implicit none
integer,intent(in) :: n1,n2
INTEGER,pointer    :: I(:,:)
integer (kind=long) :: nsize
call mem_alloc(I,n1,n2)
nsize = size(I)*mem_intsize
call mem_allocated_mem_dft(nsize)
END SUBROUTINE int_dft_allocate_2dim

SUBROUTINE int_dft_allocate_3dim(I,n1,n2,n3)
implicit none
integer,intent(in) :: n1,n2,n3
INTEGER,pointer    :: I(:,:,:)
integer (kind=long) :: nsize
call mem_alloc(I,n1,n2,n3)
nsize = size(I)*mem_intsize
call mem_allocated_mem_dft(nsize)
END SUBROUTINE int_dft_allocate_3dim

!----- DEALLOCATE INTEGER POINTERS -----!

SUBROUTINE int_dft_deallocate_1dim(I)
implicit none
INTEGER,pointer :: I(:)
integer (kind=long) :: nsize
nsize = size(I)*mem_intsize
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(I)
END SUBROUTINE int_dft_deallocate_1dim

#ifndef VAR_INT64
SUBROUTINE intlong_dft_deallocate_1dim(I)
implicit none
INTEGER(kind=long),pointer :: I(:)
integer (kind=long) :: nsize
nsize = size(I)*mem_intsize*2
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(I)
 END SUBROUTINE intlong_dft_deallocate_1dim
#endif

SUBROUTINE int_dft_deallocate_2dim(I)
implicit none
INTEGER,pointer :: I(:,:)
integer (kind=long) :: nsize
nsize = size(I)*mem_intsize
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(I)
END SUBROUTINE int_dft_deallocate_2dim

SUBROUTINE int_dft_deallocate_3dim(I)
implicit none
INTEGER,pointer :: I(:,:,:)
integer (kind=long) :: nsize
nsize = size(I)*mem_intsize
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(I)
END SUBROUTINE int_dft_deallocate_3dim

!----- ALLOCATE CHARACTER POINTERS -----!

SUBROUTINE char_dft_allocate_1dim(C,n)
implicit none
integer,intent(in)         :: n
CHARACTER(LEN=*),pointer :: C(:)
integer (kind=long) :: nsize
call mem_alloc(C,n)
nsize = mem_complexsize*size(C)
call mem_allocated_mem_dft(nsize)
END SUBROUTINE char_dft_allocate_1dim

!----- DEALLOCATE CHARACTER POINTERS -----!

SUBROUTINE char_dft_deallocate_1dim(C)
implicit none
CHARACTER(LEN=*),pointer :: C(:)
integer :: IERR
integer (kind=long) :: nsize
nsize = mem_complexsize*size(C)
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(C)
END SUBROUTINE char_dft_deallocate_1dim

!----- ALLOCATE LOGICAL POINTERS -----!

SUBROUTINE logic_dft_allocate_1dim(L,n)
implicit none
integer,intent(in) :: n
LOGICAL,pointer    :: L(:)
integer (kind=long) :: nsize
call mem_alloc(L,n)
nsize = size(L)*mem_logicalsize 
call mem_allocated_mem_dft(nsize)
END SUBROUTINE logic_dft_allocate_1dim

SUBROUTINE logic_dft_allocate_2dim(L,n1,n2)
implicit none
integer,intent(in) :: n1,n2
LOGICAL,pointer    :: L(:,:)
integer (kind=long) :: nsize
call mem_alloc(L,n1,n2)
nsize = size(L)*mem_logicalsize 
call mem_allocated_mem_dft(nsize)
END SUBROUTINE logic_dft_allocate_2dim

!----- DEALLOCATE LOGICAL POINTERS -----!

SUBROUTINE logic_dft_deallocate_1dim(L)
implicit none
LOGICAL,pointer :: L(:)
integer :: IERR
integer (kind=long) :: nsize
nsize = size(L)*mem_logicalsize
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(L)
END SUBROUTINE logic_dft_deallocate_1dim

SUBROUTINE logic_dft_deallocate_2dim(L)
implicit none
LOGICAL,pointer :: L(:,:)
integer :: IERR
integer (kind=long) :: nsize
nsize = size(L)*mem_logicalsize
call mem_deallocated_mem_dft(nsize)
call mem_dealloc(L)
END SUBROUTINE logic_dft_deallocate_2dim

!----- MEMORY HANDLING -----!

subroutine mem_allocated_mem_dft(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  IF(mem_dft_InsideOMPsection)THEN!we add to thread private variables
     mem_tp_allocated_dft = mem_tp_allocated_dft + nsize
     if (mem_tp_allocated_dft < 0) then
        write(*,*) 'XC memory negative! mem_tp_allocated_dft =', mem_tp_allocated_dft
        write(*,*) 'XC memory negative! nsize =', nsize
        call lsQUIT('Error in mem_tp_allocated_mem_dft-probably integer overflow!',-1)
     endif
     max_mem_tp_used_dft = MAX(max_mem_tp_used_dft,mem_tp_allocated_dft)
  ELSE
     mem_allocated_dft = mem_allocated_dft + nsize
     if (mem_allocated_dft < 0) then
        write(*,*) 'XC memory negative! mem_allocated_dft =', mem_allocated_dft
        write(*,*) 'XC memory negative! nsize =', nsize
        call lsQUIT('Error in mem_allocated_mem_dft - probably integer overflow!',-1)
     endif
     max_mem_used_dft = MAX(max_mem_used_dft,mem_allocated_dft)
  ENDIF
end subroutine mem_allocated_mem_dft

subroutine mem_deallocated_mem_dft(nsize)
  implicit none
  integer (kind=long), intent(in) :: nsize
  IF(mem_dft_InsideOMPsection)THEN!we add to thread private variables
     mem_tp_allocated_dft = mem_tp_allocated_dft - nsize
     if (mem_tp_allocated_dft < 0) then
        write(*,*) 'XC memory negative! mem_tp_allocated_dft =', mem_tp_allocated_dft
        call lsQUIT('Error in mem_tp_deallocated_mem_dft - something wrong with deallocation!',-1)
     endif
  else
     mem_allocated_dft = mem_allocated_dft - nsize
     if (mem_allocated_dft < 0) then
        write(*,*) 'XC memory negative! mem_allocated_dft =', mem_allocated_dft
        call lsQUIT('Error in mem_deallocated_mem_dft - something wrong with deallocation!',-1)
     endif
  endif
end subroutine mem_deallocated_mem_dft

END MODULE dft_memory_handling

