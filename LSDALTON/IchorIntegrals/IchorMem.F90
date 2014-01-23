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
   public ADD_OMP_MEM
   public REMOVE_OMP_MEM
!GLOBAL VARIABLES
   integer(KIND=long),save :: mem_allocated_ichor, max_mem_used_ichor  !Count all memory
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
!Interfaces for allocating/deallocating pointers
INTERFACE mem_ichor_alloc
  MODULE PROCEDURE real_allocate_1dim, real_allocate_2dim, &
       & real_allocate_2dim_zero, real_allocate_3dim, &
       & real_allocate_3dim_zero, real_allocate_4dim, &
       & real_allocate_5dim, &
       & int_allocate_1dim,int_allocate_2dim,int_allocate_3dim, &
       & int_allocate_1dim_zero, int_allocate_2dim_zero,  &
       & int_allocate_3dim_zero,&
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
subroutine set_ichor_memvar(MaxMemAllocated,MemAllocated)
implicit none
integer(KIND=long),intent(in) :: MaxMemAllocated,MemAllocated
mem_allocated_ichor = MemAllocated
max_mem_used_ichor = MaxMemAllocated
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

!----- ALLOCATE REAL POINTERS -----!

SUBROUTINE real_allocate_1dim(A,n)
  implicit none
  integer,intent(in)  :: n
  REAL(REALK),pointer :: A(:)
  integer :: IERR
  integer (kind=long) :: nsize
  nullify(A)
  ALLOCATE(A(n),STAT = IERR)
  IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_allocate_1dim',IERR,n
     CALL ichorquit('Error in real_allocate_1dim',-1)
  ENDIF
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE real_allocate_1dim

SUBROUTINE real_allocate_2dim(A,n1,n2)
  implicit none
  integer,intent(in)  :: n1, n2
  REAL(REALK),pointer :: A(:,:)
  integer :: IERR
  integer (kind=long) :: nsize
  nullify(A)
  ALLOCATE(A(n1,n2),STAT = IERR)
  IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_allocate_2dim',IERR,n1,n2
     call IchorQuit('Error in real_allocate_2dim',-1)
  ENDIF
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

SUBROUTINE real_allocate_2dim_zero(A,n1,n2,z1,z2)
  ! Allocates 2d arrays starting from zero index
  ! for first,second or both dimensions
  implicit none
  integer,intent(in)  :: n1, n2
  logical,intent(in)  :: z1,z2
  REAL(REALK),pointer :: A(:,:)
  integer :: IERR
  integer (kind=long) :: nsize
  integer             :: i1,i2
  !
  nullify(A)
  i1=1
  IF(z1)i1=0
  i2=1
  IF(z2)i2=0
  ALLOCATE(A(i1:n1,i2:n2),STAT = IERR)
  IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_allocate_2dim',IERR,n1,n2
     CALL ICHORQUIT('Error in real_allocate_2dim',-1)
  ENDIF
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE real_allocate_2dim_zero

SUBROUTINE real_allocate_3dim(A,n1,n2,n3)
  implicit none
  integer,intent(in)  :: n1, n2, n3
  REAL(REALK),pointer :: A(:,:,:)
  integer :: IERR
  integer (kind=long) :: nsize
  nullify(A)
  ALLOCATE(A(n1,n2,n3),STAT = IERR)
  IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_allocate_3dim',IERR,n1,n2,n3
     CALL ICHORQUIT('Error in real_allocate_3dim',-1)
  ENDIF
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE real_allocate_3dim

SUBROUTINE real_allocate_3dim_zero(A,n1,n2,n3,z1,z2,z3)
  implicit none
  integer,intent(in)  :: n1, n2, n3
  logical,intent(in)  :: z1,z2,z3
  REAL(REALK),pointer :: A(:,:,:)
  integer :: IERR
  integer (kind=long) :: nsize
  integer             :: i1,i2,i3
  nullify(A)
  i1=1
  IF(z1)i1=0
  i2=1
  IF(z2)i2=0
  i3=1
  IF(z3)i3=0
  ALLOCATE(A(i1:n1,i2:n2,i3:n3),STAT = IERR)
  IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_allocate_3dim',IERR,n1,n2,n3
     CALL ICHORQUIT('Error in real_allocate_3dim',-1)
  ENDIF
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE real_allocate_3dim_zero

SUBROUTINE real_allocate_4dim(A,n1,n2,n3,n4)
  implicit none
  integer,intent(in)  :: n1,n2,n3,n4
  REAL(REALK),pointer :: A(:,:,:,:)
  integer :: IERR
  integer (kind=long) :: nsize
  nullify(A)
  ALLOCATE(A(n1,n2,n3,n4),STAT = IERR)
  IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_allocate_4dim',IERR,n1,n2,n3,n4
     CALL ICHORQUIT('Error in real_allocate_4dim',-1)
  ENDIF
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE real_allocate_4dim

SUBROUTINE real_allocate_5dim(A,n1,n2,n3,n4,n5)
  implicit none
  integer,intent(in)  :: n1,n2,n3,n4,n5
  REAL(REALK),pointer :: A(:,:,:,:,:)
  integer :: IERR
  integer (kind=long) :: nsize
  nullify(A)
  ALLOCATE(A(n1,n2,n3,n4,n5),STAT = IERR)
  IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_allocate_5dim',IERR,n1,n2,n3,n4,n5
     CALL ICHORQUIT('Error in real_allocate_5dim',-1)
  ENDIF
  nsize = size(A,KIND=long)*mem_realsize
  call mem_allocated_ichor_mem_real(nsize)
END SUBROUTINE real_allocate_5dim

!----- DEALLOCATE REAL POINTERS -----!

SUBROUTINE real_deallocate_1dim(A)
  implicit none
  REAL(REALK),pointer :: A(:)
  integer :: IERR
  integer (kind=long) :: nsize
  nsize = size(A,KIND=long)*mem_realsize
  call mem_deallocated_ichor_mem_real(nsize)
  if (.not.ASSOCIATED(A)) then
     print *,'Memory previously released!!'
     call IchorQuit('Error in real_deallocate_1dim - memory previously released',-1)
  endif
  DEALLOCATE(A,STAT = IERR)
  IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_deallocate_1dim',IERR
     CALL ICHORQUIT('Error in real_deallocate_1dim',-1)
  ENDIF
  nullify(A)
END SUBROUTINE real_deallocate_1dim

SUBROUTINE real_deallocate_2dim(A)
implicit none
REAL(REALK),pointer :: A(:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_realsize
   call mem_deallocated_ichor_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call IchorQuit('Error in real_deallocate_2dim - memory previously released',-1)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_deallocate_2dim',IERR
     CALL ICHORQUIT('Error in real_deallocate_2dim',-1)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_2dim

SUBROUTINE real_deallocate_3dim(A)
implicit none
REAL(REALK),pointer :: A(:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_realsize
   call mem_deallocated_ichor_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call IchorQuit('Error in real_deallocate_3dim - memory previously released',-1)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_deallocate_3dim',IERR
     CALL ICHORQUIT('Error in real_deallocate_3dim',-1)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_3dim

SUBROUTINE real_deallocate_4dim(A)
implicit none
REAL(REALK),pointer :: A(:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_realsize
   call mem_deallocated_ichor_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call IchorQuit('Error in real_deallocate_4dim - memory previously released',-1)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_deallocate_4dim',IERR
     CALL ICHORQUIT('Error in real_deallocate_4dim',-1)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_4dim

SUBROUTINE real_deallocate_5dim(A)
implicit none
REAL(REALK),pointer :: A(:,:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A,KIND=long)*mem_realsize
   call mem_deallocated_ichor_mem_real(nsize)
   if (.not.ASSOCIATED(A)) then
      print *,'Memory previously released!!'
      call IchorQuit('Error in real_deallocate_5dim - memory previously released',-1)
   endif
   DEALLOCATE(A,STAT = IERR)
   IF (IERR.NE. 0) THEN
     write(*,*) 'Error in real_deallocate_5dim',IERR
     CALL ICHORQUIT('Error in real_deallocate_5dim',-1)
   ENDIF
   nullify(A)
END SUBROUTINE real_deallocate_5dim

!----- ALLOCATE INTEGER POINTERS -----!

SUBROUTINE int_allocate_1dim(I,n)
implicit none
integer,intent(in)  :: n
INTEGER,pointer     :: I(:)
integer :: IERR
integer (kind=long) :: nsize
nullify(I)
ALLOCATE(I(n),STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in int_allocate_1dim',IERR,n
   call IchorQuit('Error in int_allocate_1dim',-1)
ENDIF
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_ichor_mem_integer(nsize)
END SUBROUTINE int_allocate_1dim

SUBROUTINE int_allocate_2dim(I,n1,n2)
implicit none
integer,intent(in) :: n1,n2
INTEGER,pointer    :: I(:,:)
integer :: IERR
integer (kind=long) :: nsize
nullify(I)
ALLOCATE(I(n1,n2),STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in int_allocate_2dim',IERR,n1,n2
   call IchorQuit('Error in int_allocate_2dim',-1)
ENDIF
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_ichor_mem_integer(nsize)
END SUBROUTINE int_allocate_2dim

SUBROUTINE int_allocate_3dim(I,n1,n2,n3)
implicit none
integer,intent(in) :: n1,n2,n3
INTEGER,pointer    :: I(:,:,:)
integer :: IERR
integer (kind=long) :: nsize
nullify(I)
ALLOCATE(I(n1,n2,n3),STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in int_allocate_3dim',IERR,n1,n2,n3
   call IchorQuit('Error in int_allocate_3dim',-1)
ENDIF
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_ichor_mem_integer(nsize)
END SUBROUTINE int_allocate_3dim

SUBROUTINE int_allocate_4dim(I,n1,n2,n3,n4)
implicit none
integer,intent(in) :: n1,n2,n3,n4
INTEGER,pointer    :: I(:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
nullify(I)
ALLOCATE(I(n1,n2,n3,n4),STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in int_allocate_4dim',IERR,n1,n2,n3,n4
   call IchorQuit('Error in int_allocate_4dim',-1)
ENDIF
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_ichor_mem_integer(nsize)
END SUBROUTINE int_allocate_4dim

SUBROUTINE int_allocate_1dim_zero(I,n,z1)
implicit none
integer,intent(in)  :: n
INTEGER,pointer     :: I(:)
logical,intent(in) :: z1
integer :: IERR,i1
integer (kind=long) :: nsize
i1=1
IF(z1)i1=0
nullify(I)
ALLOCATE(I(i1:n),STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in int_allocate_1dim',IERR,n
   call IchorQuit('Error in int_allocate_1dim',-1)
ENDIF
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_ichor_mem_integer(nsize)
END SUBROUTINE int_allocate_1dim_zero

SUBROUTINE int_allocate_2dim_zero(I,n1,n2,z1,z2)
implicit none
integer,intent(in) :: n1,n2
INTEGER,pointer    :: I(:,:)
logical,intent(in) :: z1,z2
integer :: IERR,i1,i2
integer (kind=long) :: nsize
i1=1
IF(z1)i1=0
i2=1
IF(z2)i2=0
nullify(I)
ALLOCATE(I(i1:n1,i2:n2),STAT = IERR)
   IF (IERR.NE. 0) THEN
     write(*,*) 'Error in int_allocate_2dim',IERR,n1,n2
     call IchorQuit('Error in int_allocate_2dim',-1)
   ENDIF
   nsize = size(I,KIND=long)*mem_intsize
   call mem_allocated_ichor_mem_integer(nsize)
END SUBROUTINE int_allocate_2dim_zero

SUBROUTINE int_allocate_3dim_zero(I,n1,n2,n3,z1,z2,z3)
implicit none
integer,intent(in) :: n1,n2,n3
logical,intent(in) :: z1,z2,z3
INTEGER,pointer    :: I(:,:,:)
integer :: IERR,i1,i2,i3
integer (kind=long) :: nsize
i1=1
IF(z1)i1=0
i2=1
IF(z2)i2=0
i3=1
IF(z3)i3=0
nullify(I)
ALLOCATE(I(i1:n1,i2:n2,i3:n3),STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in int_allocate_3dim',IERR,n1,n2,n3
   call IchorQuit('Error in int_allocate_3dim',-1)
ENDIF
nsize = size(I,KIND=long)*mem_intsize
call mem_allocated_ichor_mem_integer(nsize)
END SUBROUTINE int_allocate_3dim_zero

!----- DEALLOCATE INTEGER POINTERS -----!
SUBROUTINE int_deallocate_1dim(I)
implicit none
INTEGER,pointer :: I(:)
integer :: IERR
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_deallocated_ichor_mem_integer(nsize)
if (.not.ASSOCIATED(I)) then
   print *,'Memory previously released!!'
   call IchorQuit('Error in int_deallocate_1dim - memory previously released',-1)
endif
DEALLOCATE(I,STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in int_deallocate_1dim',IERR
   CALL ICHORQUIT('Error in int_deallocate_1dim',-1)
ENDIF
nullify(I)
END SUBROUTINE int_deallocate_1dim

SUBROUTINE int_deallocate_2dim(I)
implicit none
INTEGER,pointer :: I(:,:)
integer :: IERR
integer (kind=long) :: nsize
nsize = size(I,KIND=long)*mem_intsize
call mem_deallocated_ichor_mem_integer(nsize)
if (.not.ASSOCIATED(I)) then
   print *,'Memory previously released!!'
   call IchorQuit('Error in int_deallocate_2dim - memory previously released',-1)
endif
DEALLOCATE(I,STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in int_deallocate_2dim',IERR
   CALL ICHORQUIT('Error in int_deallocate_2dim',-1)
ENDIF
nullify(I)
END SUBROUTINE int_deallocate_2dim

SUBROUTINE int_deallocate_3dim(I)
implicit none
INTEGER,pointer :: I(:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(I,KIND=long)*mem_intsize
   call mem_deallocated_ichor_mem_integer(nsize)
   if (.not.ASSOCIATED(I)) then
      print *,'Memory previously released!!'
      call IchorQuit('Error in int_deallocate_3dim - memory previously released',-1)
   endif
   DEALLOCATE(I,STAT = IERR)
   IF (IERR.NE. 0) THEN
      write(*,*) 'Error in int_deallocate_3dim',IERR
      CALL ICHORQUIT('Error in int_deallocate_3dim',-1)
   ENDIF
   nullify(I)
 END SUBROUTINE int_deallocate_3dim

SUBROUTINE logic_allocate_1dim(L,n)
implicit none
integer,intent(in)  :: n
logical,pointer     :: L(:)
integer :: IERR
integer (kind=long) :: nsize
nullify(L)
ALLOCATE(L(n),STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in logic_allocate_1dim',IERR,n
   call IchorQuit('Error in logic_allocate_1dim',-1)
ENDIF
nsize = size(L,KIND=long)*mem_logicalsize
call mem_allocated_ichor_mem_logical(nsize)
END SUBROUTINE logic_allocate_1dim

SUBROUTINE logic_allocate_2dim(L,n1,n2)
implicit none
integer,intent(in) :: n1,n2
logical,pointer    :: L(:,:)
integer :: IERR
integer (kind=long) :: nsize
nullify(L)
ALLOCATE(L(n1,n2),STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in logic_allocate_2dim',IERR,n1,n2
   call IchorQuit('Error in logic_allocate_2dim',-1)
ENDIF
nsize = size(L,KIND=long)*mem_logicalsize
call mem_allocated_ichor_mem_logical(nsize)
END SUBROUTINE logic_allocate_2dim

!----- DEALLOCATE LOGICAL POINTERS -----!
SUBROUTINE logic_deallocate_1dim(L)
implicit none
logical,pointer :: L(:)
integer :: IERR
integer (kind=long) :: nsize
nsize = size(L,KIND=long)*mem_logicalsize
call mem_deallocated_ichor_mem_logical(nsize)
if (.not.ASSOCIATED(L)) then
   print *,'Memory previously released!!'
   call IchorQuit('Error in logic_deallocate_1dim - memory previously released',-1)
endif
DEALLOCATE(L,STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in logic_deallocate_1dim',IERR
   CALL ICHORQUIT('Error in logic_deallocate_1dim',-1)
ENDIF
nullify(L)
END SUBROUTINE logic_deallocate_1dim

SUBROUTINE logic_deallocate_2dim(L)
implicit none
LOGICAL,pointer :: L(:,:)
integer :: IERR
integer (kind=long) :: nsize
nsize = size(L,KIND=long)*mem_logicalsize
call mem_deallocated_ichor_mem_logical(nsize)
if (.not.ASSOCIATED(L)) then
   print *,'Memory previously released!!'
   call IchorQuit('Error in logic_deallocate_2dim - memory previously released',-1)
endif
DEALLOCATE(L,STAT = IERR)
IF (IERR.NE. 0) THEN
   write(*,*) 'Error in logic_deallocate_2dim',IERR
   CALL ICHORQUIT('Error in logic_deallocate_2dim',-1)
ENDIF
nullify(L)
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

