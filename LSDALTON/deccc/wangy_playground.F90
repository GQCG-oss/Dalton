!> @file
!> Author: Yang M. Wang
!> Date: April 2013
!> Have fun and live life

module wangy_playground_module

!!$  use dec_fragment_utils       !Routines: dec_simple_dgemm
!!$  use array4_simple_operations !Routines: mat_transpose 
!!$
!!$  !NB! Need this to make the routine available for f12_integrals.F90
!!$  public :: timings_for_dgemm 
!!$  private
!!$
!!$contains
!!$
!!$  !> Brief: Compares the timining from Kaspers dgemm wrapper with the lapack dgemm.
!!$  !> We initiate random matrices for the test.
!!$  !> Author: Yang M. Wang
!!$  !> Date: April 2013
!!$  subroutine timings_for_dgemm
!!$    implicit none
!!$    !real(realk), pointer :: A(:,:)  !F90 way of doing dynamical memory allocation
!!$    !call mem_alloc(A,2,2)
!!$    real(realk) :: A(2,2), tA(2,2), B(2,2), C(2,2), D(2,3), E(3,1), F(2,1), H(2,3), G(2,2) 
!!$    real(realk) :: m,n,k
!!$
!!$    A = 0E0_realk         ! Zeroing out to be on the safe side
!!$    print *, 
!!$    print *, "YOLO!!! --- Welcome to wangy playground ---" 
!!$    print *, 
!!$    print *, " --- Testing matrix multiplication --- "
!!$    print *,
!!$    print *, "Matrix A:"
!!$    A(1,1) = 1
!!$    A(2,2) = A(1,1)
!!$    A(1,2) = 2
!!$    A(2,1) = A(1,2) 
!!$    call matrix_print(A,2,2)
!!$    print *, "Matrix B:"
!!$    B(1,1) = 2
!!$    B(2,2) = B(1,1)
!!$    B(1,2) = 2
!!$    B(2,1) = B(1,2) 
!!$    call matrix_print(B,2,2)
!!$    call dec_simple_dgemm(2,2,2,A,B,C,'n','n')  
!!$    print *, "Matrix C:"
!!$    call matrix_print(C,2,2)
!!$
!!$    D(1,1) = 1
!!$    D(1,2) = 2
!!$    D(1,3) = 3
!!$    D(2,1) = D(1,1) 
!!$    D(2,2) = D(1,2) 
!!$    D(2,3) = D(1,3) 
!!$
!!$    E(1,1) = 1
!!$    E(2,1) = 1
!!$    E(3,1) = 1
!!$
!!$    H(1,1) = -1
!!$    H(1,2) = 2
!!$    H(1,3) = 1
!!$
!!$    H(2,1) = 1
!!$    H(2,2) = 0
!!$    H(2,3) = 2
!!$
!!$    print *, "Matrix D:"
!!$    call matrix_print(D,2,3)
!!$
!!$    print *, "Matrix E:"
!!$    call matrix_print(E,3,1)
!!$
!!$    call dec_simple_dgemm(2,3,1,D,E,F,'n','n')
!!$    print *, "Matrix F:"
!!$    call matrix_print(F,2,1)
!!$
!!$    print *, " --- Testing matrix multiplication rectangular --- "
!!$    print *, "Matrix D:"
!!$    call matrix_print(D,2,3)
!!$
!!$    print *, "Matrix H:"
!!$    call matrix_print(H,2,3)
!!$
!!$    call dec_simple_dgemm(2,3,2,D,H,G,'n','n')
!!$    print *, "Matrix G:"
!!$    call matrix_print(F,2,2)
!!$
!!$    call dec_simple_dgemm(2,3,2,D,H,G,'n','t')
!!$    print *, "Matrix G:"
!!$    call matrix_print(F,2,2)
!!$  
!!$    m = 2;
!!$    n = 2;
!!$    k = 2;
!!$    !call DGEMM('N','N',m,n,k,1.0E0_realk,A,2,B,2,0.0E0_realk,C,2)   
!!$    print *,
!!$    print *, " --- Testing matrix transpose --- "
!!$    print *,
!!$           
!!$    A(1,1) = 1
!!$    A(2,2) = 0
!!$    A(1,2) = 3
!!$    A(2,1) = 2 
!!$
!!$    print *, "Matrix A:"
!!$    call matrix_print(A,2,2)
!!$    call mat_transpose(A,2,2,tA)
!!$    print *, "Matrix A^T:"
!!$    call matrix_print(tA,2,2)
!!$
!!$    !REAL :: r(5,5)
!!$    !call init_random_seed()         ! see example of RANDOM_SEED
!!$    !call random_number(r)
!!$ 
!!$  end subroutine timings_for_dgemm
!!$
!!$  subroutine matrix_print(A,m,n)
!!$    implicit none
!!$    integer :: i, j 
!!$    integer :: m, n
!!$    real(realk) :: A(m,n) 
!!$
!!$    !do, i=1,m
!!$    !   write(*,5) (A(i,j), j=1, n)
!!$    !enddo
!!$    
!!$    do i=1,m
!!$       do j=1,n
!!$          write(*,"(F5.2)",advance='no') A(i,j)
!!$       enddo
!!$       write(*,*)
!!$    enddo
!!$    write(*,*)
!!$    
!!$  end subroutine matrix_print
!!$ 
end module wangy_playground_module
