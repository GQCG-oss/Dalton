MODULE mm_mem
  USE memory_handling
  use precision
INTERFACE mem_alloc_fmm
  MODULE PROCEDURE real_allocate_1dim_fmm, real_allocate_2dim_fmm, &
     &             real_allocate_2dim_zero_fmm,                & 
     &             real_allocate_3dim_fmm, real_allocate_4dim_fmm, &
     &             real_allocate_5dim_fmm, real_allocate_5dim_zero_fmm, &
     &             int_allocate_1dim_fmm,  int_allocate_2dim_fmm,  &
     &             char_allocate_1dim_fmm,                     &
     &             logic_allocate_1dim_fmm
END INTERFACE

INTERFACE mem_dealloc_fmm
  MODULE PROCEDURE real_deallocate_1dim_fmm, real_deallocate_2dim_fmm, &
     &             real_deallocate_3dim_fmm, real_deallocate_4dim_fmm, &
     &             real_deallocate_5dim_fmm,                       &
     &             int_deallocate_1dim_fmm,  int_deallocate_2dim_fmm,  &
     &             char_deallocate_1dim_fmm,                       &
     &             logic_deallocate_1dim_fmm
END INTERFACE

PUBLIC:: mem_alloc_fmm, mem_dealloc_fmm, add_fmm_mem, remove_fmm_mem 

PRIVATE

CONTAINS
!----- ALLOCATE REAL POINTERS -----! 

SUBROUTINE real_allocate_1dim_fmm(A,n)
implicit none
integer,intent(in)  :: n
REAL(REALK),pointer :: A(:)
integer (kind=long) :: nsize
call mem_alloc(A,n)
nsize = size(A)*mem_realsize
call add_fmm_mem(nsize,1)
END SUBROUTINE real_allocate_1dim_fmm

SUBROUTINE real_allocate_2dim_fmm(A,n1,n2)
implicit none
integer,intent(in)  :: n1, n2
REAL(REALK),pointer :: A(:,:)
integer (kind=long) :: nsize
call mem_alloc(A,n1,n2)
nsize = size(A)*mem_realsize
call add_fmm_mem(nsize,2)
END SUBROUTINE real_allocate_2dim_fmm

SUBROUTINE real_allocate_2dim_zero_fmm(A,n1,n2,First,Second)
implicit none
integer,intent(in)  :: n1, n2
REAL(REALK),pointer :: A(:,:)
integer (kind=long) :: nsize
Logical :: First, Second
call mem_alloc(A,n1,n2,First,Second)
nsize = size(A)*mem_realsize
call add_fmm_mem(nsize,2)
END SUBROUTINE real_allocate_2dim_zero_fmm

SUBROUTINE real_allocate_3dim_fmm(A,n1,n2,n3)
implicit none
integer,intent(in)  :: n1, n2, n3
REAL(REALK),pointer :: A(:,:,:)
integer (kind=long) :: nsize
call mem_alloc(A,n1,n2,n3)
nsize = size(A)*mem_realsize
call add_fmm_mem(nsize,6)
END SUBROUTINE real_allocate_3dim_fmm

SUBROUTINE real_allocate_4dim_fmm(A,n1,n2,n3,n4)
implicit none
integer,intent(in)  :: n1,n2,n3,n4
REAL(REALK),pointer :: A(:,:,:,:)
integer (kind=long) :: nsize
call mem_alloc(A,n1,n2,n3,n4)
nsize = size(A)*mem_realsize
call add_fmm_mem(nsize,0)
END SUBROUTINE real_allocate_4dim_fmm

SUBROUTINE real_allocate_5dim_fmm(A,n1,n2,n3,n4,n5)
implicit none
integer,intent(in)  :: n1,n2,n3,n4,n5
REAL(REALK),pointer :: A(:,:,:,:,:)
integer (kind=long) :: itest, nsize
call mem_alloc(A,n1,n2,n3,n4,n5)
nsize = size(A)*mem_realsize
call add_fmm_mem(nsize,0)
END SUBROUTINE real_allocate_5dim_fmm

SUBROUTINE real_allocate_5dim_zero_fmm(A,n1,n2,n3,n4,n5,z1,z2,z3,z4,z5)
! Allocates 5d arrays starting from zero index
! for some or all dimensions
implicit none
integer,intent(in)  :: n1,n2,n3,n4,n5
logical,intent(in)  :: z1,z2,z3,z4,z5
REAL(REALK),pointer :: A(:,:,:,:,:)
integer (kind=long) :: itest, nsize
call mem_alloc(A,n1,n2,n3,n4,n5,z1,z2,z3,z4,z5)
nsize = size(A)*mem_realsize
call add_fmm_mem(nsize,0)
END SUBROUTINE real_allocate_5dim_zero_fmm

!----- DEALLOCATE REAL POINTERS -----!

SUBROUTINE real_deallocate_1dim_fmm(A)
implicit none
REAL(REALK),pointer :: A(:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call remove_fmm_mem(nsize,1)
   call mem_dealloc(A)
END SUBROUTINE real_deallocate_1dim_fmm

SUBROUTINE real_deallocate_2dim_fmm(A)
implicit none
REAL(REALK),pointer :: A(:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call remove_fmm_mem(nsize,2)
   call mem_dealloc(A)
END SUBROUTINE real_deallocate_2dim_fmm

SUBROUTINE real_deallocate_3dim_fmm(A)
implicit none
REAL(REALK),pointer :: A(:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call remove_fmm_mem(nsize,6)
   call mem_dealloc(A)
END SUBROUTINE real_deallocate_3dim_fmm

SUBROUTINE real_deallocate_4dim_fmm(A)
implicit none
REAL(REALK),pointer :: A(:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call remove_fmm_mem(nsize,0)
   call mem_dealloc(A)
END SUBROUTINE real_deallocate_4dim_fmm

SUBROUTINE real_deallocate_5dim_fmm(A)
implicit none
REAL(REALK),pointer :: A(:,:,:,:,:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(A)*mem_realsize
   call remove_fmm_mem(nsize,0)
   call mem_dealloc(A)
END SUBROUTINE real_deallocate_5dim_fmm

!----- ALLOCATE INTEGER POINTERS -----!

SUBROUTINE int_allocate_1dim_fmm(I,n)
implicit none
integer,intent(in)  :: n
INTEGER,pointer     :: I(:)
integer (kind=long) :: nsize
   call mem_alloc(I,n)
   nsize = size(I)*mem_intsize
   call add_fmm_mem(nsize,7)
END SUBROUTINE int_allocate_1dim_fmm

SUBROUTINE int_allocate_2dim_fmm(I,n1,n2)
implicit none
integer,intent(in) :: n1,n2
INTEGER,pointer    :: I(:,:)
integer (kind=long) :: nsize
   call mem_alloc(I,n1,n2)
   nsize = size(I)*mem_intsize
   call add_fmm_mem(nsize,7)
END SUBROUTINE int_allocate_2dim_fmm

!----- DEALLOCATE INTEGER POINTERS -----!

SUBROUTINE int_deallocate_1dim_fmm(I)
implicit none
INTEGER,pointer :: I(:)
integer (kind=long) :: nsize
   nsize = size(I)*mem_intsize
   call remove_fmm_mem(nsize,7)
   call mem_dealloc(I)
END SUBROUTINE int_deallocate_1dim_fmm

SUBROUTINE int_deallocate_2dim_fmm(I)
implicit none
INTEGER,pointer :: I(:,:)
integer (kind=long) :: nsize
   nsize = size(I)*mem_intsize
   call remove_fmm_mem(nsize,7)
   call mem_dealloc(I)
END SUBROUTINE int_deallocate_2dim_fmm

!----- ALLOCATE CHARACTER POINTERS -----!

SUBROUTINE char_allocate_1dim_fmm(C,n)
implicit none
integer,intent(in)         :: n
CHARACTER(LEN=*),pointer :: C(:)
integer :: IERR
integer (kind=long) :: nsize
   call mem_alloc(C,n)
   nsize = mem_complexsize*size(C)
   call add_fmm_mem(nsize,0)
END SUBROUTINE char_allocate_1dim_fmm

!----- DEALLOCATE CHARACTER POINTERS -----!

SUBROUTINE char_deallocate_1dim_fmm(C)
implicit none
CHARACTER(LEN=*),pointer :: C(:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = mem_complexsize*size(C)
   call remove_fmm_mem(nsize,0)
   call mem_dealloc(C)
END SUBROUTINE char_deallocate_1dim_fmm

!----- ALLOCATE LOGICAL POINTERS -----!

SUBROUTINE logic_allocate_1dim_fmm(L,n)
implicit none
integer,intent(in) :: n
LOGICAL,pointer    :: L(:)
integer (kind=long) :: nsize
   call mem_alloc(L,n)
   nsize = mem_logicalsize*size(L)
   call add_fmm_mem(nsize,0)
END SUBROUTINE logic_allocate_1dim_fmm

!----- DEALLOCATE LOGICAL POINTERS -----!

SUBROUTINE logic_deallocate_1dim_fmm(L)
implicit none
LOGICAL,pointer :: L(:)
integer :: IERR
integer (kind=long) :: nsize
   nsize = size(L)*mem_logicalsize
   call remove_fmm_mem(nsize,0)
   call mem_dealloc(L)
END SUBROUTINE logic_deallocate_1dim_fmm

SUBROUTINE add_fmm_mem(nsize,opt)
implicit none
integer (kind=long) :: nsize
integer :: opt

call mem_allocated_mem_FMM(nsize,.false.)

END SUBROUTINE add_fmm_mem

SUBROUTINE remove_fmm_mem(nsize,opt)
implicit none
integer (kind=long) :: nsize
integer :: opt

call mem_deallocated_mem_FMM(nsize,.false.)

END SUBROUTINE remove_fmm_mem

END MODULE mm_mem
