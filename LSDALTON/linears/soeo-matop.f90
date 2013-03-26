module soeo_matop

use precision
use matrix_module
use Matrix_Operations
use Matrix_Operations_aux
!use matrix_util
!use dal_interface
!use decompMod

implicit none

!Contains:
!  soeo_diag_elma
!  soeo_diag_elmb
!  soeo_dotproduct
!  soeo_norm
!  soeo_normalize
!  soeo_scal
!  soeo_orthonormalize
!  soeo_orthogonalize
!  soeo_daxpy
!!  determinant (should not really be here, only used for debugging)
!!  invert (should not really be here, only used for debugging)
!!  matmult (should not really be here, only used for debugging)

contains

!> \brief Gives the diagonal element (i,i) of a matrix A (alpha part)
!> \author C. Nygaard
!> \date 2010
!> \param A The matrix
!> \param i The place of the diagonal element
!=======================================================================
function soeo_diag_elma (A, i)

implicit none

real(realk) soeo_diag_elma
type(matrix), intent(in) :: A
integer, intent(in)      :: i
real(realk)              :: tmp(2)

call mat_get_ab_elms (A, i, i, tmp)
soeo_diag_elma = tmp(1)

end function soeo_diag_elma
!=======================================================================

!> \brief Gives the diagonal element (i,i) of a matrix A (beta part)
!> \author C. Nygaard
!> \date 2010
!> \param A The matrix
!> \param i The place of the diagonal element
!=======================================================================
function soeo_diag_elmb (A, i)

implicit none

real(realk) soeo_diag_elmb
type(matrix), intent(in) :: A
integer, intent(in)      :: i
real(realk)              :: tmp(2)

call mat_get_ab_elms (A, i, i, tmp)
soeo_diag_elmb = tmp(2)

end function soeo_diag_elmb
!=======================================================================

!> \brief Calculates the dot-product of a matrix-vector
!> \author C. Nygaard
!> \date 2010
!> \param a_m Matrix part of the first matrix-vector
!> \param a_v Vector part of the first matrix-vector
!> \param b_m Matrix part of the second matrix-vector
!> \param b_v Vector part of the second matrix-vector
!   a_m, a_v dot b_m, b_v
!>  soeo_dotproduct = 1/2*dot(a_m, b_m) + dot(a_v, b_v)
!=======================================================================
function soeo_dotproduct (a_m, a_v, b_m, b_v)

implicit none

!I/O
real(realk) soeo_dotproduct
type(matrix), intent(in) :: a_m, b_m, a_v, b_v
!Other
real(realk)              :: dot_m, dot_v

!Works because the m-part is allways antisymmetric, so that it contains
! half as many independent variables as nonzero elements

!dot_m = mat_dotproduct (a_m, b_m)
dot_m = 0.5E0_realk*mat_dotproduct (a_m, b_m)
dot_v = mat_dotproduct (a_v, b_v)
soeo_dotproduct = dot_m + dot_v

end function soeo_dotproduct
!=======================================================================

!> \brief Calculates the norm of a matrix-vector
!> \author C. Nygaard
!> \date 2010
!> \param x_m Matrix part of the matrix-vector
!> \param x_v Vector part of the matrix-vector
!=======================================================================
function soeo_norm (x_m, x_v)

implicit none

!I/O
real(realk) soeo_norm
type(matrix), intent(in) :: x_m, x_v
!Other
real(realk)              :: sqnorm

sqnorm = soeo_dotproduct (x_m, x_v, x_m, x_v)
soeo_norm = sqrt(sqnorm)

end function soeo_norm
!=======================================================================

!> \brief Normalizes a matrix-vector
!> \author C. Nygaard
!> \date 2010
!> \param a_m Matrix part of the matrix-vector
!> \param a_v Vector part of the matrix-vector
!=======================================================================
subroutine soeo_normalize (a_m, a_v)

implicit none

!I/O
type(matrix), intent(inout) :: a_m, a_v
!Other
real(realk)                 :: anorm

anorm = soeo_norm (a_m, a_v)
call soeo_scal (1.0E0_realk/anorm, a_m, a_v)

end subroutine soeo_normalize
!=======================================================================

!> \brief Multiplies a matrix-vector with a scalar alpha
!> \author C. Nygaard
!> \date 2010
!> \param alpha The scalar
!> \param a_m Matrix part of the matrix-vector
!> \param a_v Vector part of the matrix-vector
!=======================================================================
subroutine soeo_scal (alpha, a_m, a_v)

implicit none

!I/O
real(realk), intent(in)        :: alpha
type(matrix), intent(inout)    :: a_m, a_v
!Other
integer                        :: N

call mat_scal (alpha, a_m)
call mat_scal (alpha, a_v)

end subroutine soeo_scal
!=======================================================================

!> \brief Orthonormalizes matrix-vector b to matrix-vector a
!> \author C. Nygaard
!> \date 2010
!> \param a_m Matrix part of the first matrix-vector (in)
!> \param a_v Vector part of the first matrix-vector (in)
!> \param b_m Matrix part of the second matrix-vector (inout)
!> \param b_v Vector part of the second matrix-vector (inout)
!> \param err Is true if norm of b after orthogonalization is too small
!>
!> 1/2a_m*b_m + a_v*b_v = 0
!> soeo_norm(b_m,b_v) = 1/2norm(b_m) + norm(b_v) = 1
!=======================================================================
subroutine soeo_orthonormalize (a_m, a_v, b_m, b_v, err)

implicit none

!I/O
type(matrix), intent(in)    :: a_m, a_v
type(matrix), intent(inout) :: b_m, b_v
logical, intent(out)        :: err
!Other
real(realk)                 :: num, denom, proj, normb

err = .false.

num = soeo_dotproduct (a_m, a_v, b_m, b_v)
denom = soeo_dotproduct (a_m, a_v, a_m, a_v)
if (denom > 1.0E-8_realk) then
  proj = num / denom
else
  proj = 0.0E0_realk
endif
call soeo_daxpy (-proj, a_m, a_v, b_m, b_v) !b = b - bproj

normb = soeo_norm (b_m, b_v)

if (normb < 1.0E-8_realk) then
  err = .true.
else
  call soeo_normalize (b_m, b_v)
endif

end subroutine soeo_orthonormalize
!=======================================================================

!!> \brief Orthonormalizes matrix-vector b to matrix-vector a
!!> \author C. Nygaard
!!> \date Mar 21 2012
!!> \param a_m Matrix part of the first matrix-vector (in)
!!> \param a_v Vector part of the first matrix-vector (in)
!!> \param b_m Matrix part of the second matrix-vector (inout)
!!> \param b_v Vector part of the second matrix-vector (inout)
!!> \param err Is true if norm of b after orthogonalization is too small
!!>
!!> a_m*b_m = 0
!!> a_v*b_v = 0
!!> soeo_norm(b_m,b_v) = 1/2norm(b_m) + norm(b_v) = 1
!!=======================================================================
!subroutine soeo_strict_orthonormalize (a_m, a_v, b_m, b_v, err)
!
!implicit none
!
!!I/O
!type(matrix), intent(in)    :: a_m, a_v
!type(matrix), intent(inout) :: b_m, b_v
!logical, intent(out)        :: err
!!Other
!real(realk)                 :: proj_m, proj_v, normb
!
!err = .false.
!
!proj_m = mat_dotproduct (a_m, b_m)/mat_sqnorm2(a_m)
!proj_v = mat_dotproduct (a_v, b_v)/mat_sqnorm2(a_v)
!!b = b - bproj
!call mat_daxpy (-proj_m, a_m, b_m)
!call mat_daxpy (-proj_v, a_v, b_v)
!
!normb = soeo_norm (b_m, b_v)
!
!if (normb < 1.0E-8_realk) then
!  err = .true.
!else
!  call soeo_normalize (b_m, b_v)
!endif
!
!end subroutine soeo_strict_orthonormalize
!!=======================================================================
!
!!> \brief Makes matrix-vector b orthogonal to matrix-vector a
!!> \author C. Nygaard
!!> \date 2010
!!> \param a_m Matrix part of the first matrix-vector (in)
!!> \param a_v Vector part of the first matrix-vector (in)
!!> \param b_m Matrix part of the second matrix-vector (inout)
!!> \param b_v Vector part of the second matrix-vector (inout)
!!> \param err Is true if norm of b after orthogonalization is too small
!!=======================================================================
!subroutine soeo_orthogonalize (a_m, a_v, b_m, b_v, err)
!
!implicit none
!
!!I/O
!type(matrix), intent(in)    :: a_m, a_v
!type(matrix), intent(inout) :: b_m, b_v
!logical, intent(out)        :: err
!!Other
!real(realk)                 :: num, denom, proj, normb
!
!err = .false.
!
!num = soeo_dotproduct (a_m, a_v, b_m, b_v)
!denom = soeo_dotproduct (a_m, a_v, a_m, a_v)
!proj = num / denom
!call soeo_daxpy (-proj, a_m, a_v, b_m, b_v) !b = b - bproj
!
!normb = soeo_norm (b_m, b_v)
!if (normb < 1.0E-8_realk) then
!  err = .true.
!endif
!
!end subroutine soeo_orthogonalize
!!=======================================================================

!> \brief Calculates y = alpha*x + y for matrix-vectors
!> \author C. Nygaard
!> \date 2010
!> \param alpha The scalar
!> \param x_m Matrix part of the first matrix-vector (in)
!> \param x_v Vector part of the first matrix-vector (in)
!> \param y_m Matrix part of the second matrix-vector (inout)
!> \param y_v Vector part of the second matrix-vector (inout)
!=======================================================================
subroutine soeo_daxpy (alpha, x_m, x_v, y_m, y_v)

implicit none

!I/O
real(realk), intent(in)     :: alpha
type(matrix), intent(in)    :: x_m, x_v
type(matrix), intent(inout) :: y_m, y_v

if (x_m%nrow == y_m%nrow .and. x_m%ncol == y_m%ncol &
                       & .and. x_v%nrow == y_v%nrow) then

  call mat_daxpy (alpha, x_m, y_m)
  call mat_daxpy (alpha, x_v, y_v)

else
  call lsquit ("Wrong dimensions in soeo_daxpy",-1)
endif

end subroutine soeo_daxpy
!=======================================================================

!function determinant (Ain, N)
!
!implicit none
!
!real(realk) determinant
!integer, intent(in)     :: N
!real(realk), intent(in) :: Ain(:,:)
!
!real(realk)             :: A(N,N)
!integer                 :: i, j, IERR
!integer                 :: IV1(N)
!real(realk)             :: er(N), ei(N), X(N,N), FV1(N), FV2(N)
!logical                 :: symmetric
!
!A = Ain(1:N,1:N)
!
!symmetric = .true.
!do i=1,N
!  do j=1,N
!    if (A(i,j)-A(j,i) > 1.0E-5_realk .or. A(i,j)-A(j,i) < -1.0E-5_realk) symmetric = .false.
!  enddo
!enddo
!
!if (symmetric) then
!  call RS (N, N, A, er, 1, X, FV1, FV2, IERR)
!  determinant=1.0E0_realk
!  do i=1,N
!    determinant = determinant * er(i)
!  enddo
!else
!  call RG (N, N, A, er, ei, 1, X, IV1, FV1, IERR)
!  do i=1,N
!    if (ei(i)<-1.0E-5_realk .or. ei(i)>1.0E-5_realk) then
!      call lsquit ('A has complex eigenvalues, soeo_redspace/determinant',-1)
!    endif
!  enddo
!  determinant=1.0E0_realk
!  do i=1,N
!    determinant = determinant * er(i)
!  enddo
!endif
!
!end function determinant
!
!subroutine invert (A, Am1)
!
!implicit none
!
!real(realk), intent(in) :: A(:,:)
!real(realk), intent(out) :: Am1(:,:)
!
!integer :: N, IERR, i, j
!integer, allocatable :: IV1(:)
!real(realk), allocatable :: er(:), X(:,:), FV1(:), FV2(:), D(:,:), tmp(:,:)
!logical :: sym
!
!N = size(A(1,:))
!allocate (IV1(N), er(N), X(N,N), FV1(N), FV2(N), D(N,N), tmp(N,N))
!
!sym = .true.
!do i=1,N
!  do j=i,N
!    if (A(i,j)-A(j,i) > 1.0E-5_realk .or. A(i,j)-A(j,i) < -1.0E-5_realk) sym = .false.
!  enddo
!enddo
!
!if (sym) then
!  call RS (N, N, A, er, 1, X, FV1, FV2, IERR)
!else
!  call lsquit ('A is not symmetric (soeo_matop/invert)',-1)
!endif
!
!D = 0.0E0_realk
!do i=1,N
!  D(i,i) = 1.0E0_realk / er(i)
!enddo
!call matmul22 (X, .false., D, .false., tmp)
!call matmul22 (tmp, .false., X, .true., Am1)
!
!deallocate (IV1, er, X, FV1, FV2, D, tmp)
!
!end subroutine invert
!
!subroutine matmul21 (A, B, C)
!
!implicit none
!
!!A*B = C
!
!real(realk) :: A(:,:), B(:), C(:)
!integer :: N, i, j
!
!N = size(B)
!
!do i=1,N
!  C(i) = 0.0E0_realk
!  do j=1,N
!    C(i) = C(i) + A(i,j)*B(j)
!  enddo
!enddo
!
!end subroutine matmul21
!
!subroutine matmul22 (A, ta, B, tb, C)
!
!implicit none
!
!!A*B = C
!
!real(realk) :: A(:,:), B(:,:), C(:,:)
!logical :: ta, tb
!integer :: N, i, j, k
!real(realk) :: aval, bval
!
!N = size(A(1,:))
!
!do i=1,N
!  do j=1,N
!    C(i,j) = 0.0E0_realk
!    do k=1,N
!      if (ta) then
!        aval = A(k,i)
!      else
!        aval = A(i,k)
!      endif
!      if (tb) then
!        bval = B(j,k)
!      else
!        bval = B(k,j)
!      endif
!      C(i,j) = C(i,j) + aval*bval
!    enddo
!  enddo
!enddo
!
!end subroutine matmul22

end module soeo_matop
