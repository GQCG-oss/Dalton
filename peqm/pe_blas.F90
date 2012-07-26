module blas_f90

    use double_precision

    implicit none

contains

!------------------------------------------------------------------------------

function nrm2(x)

    real(dp), external :: dnrm2

    real(dp) :: nrm2
    real(dp), dimension(:), intent(in) :: x

    integer :: n, incx

    incx = 1

    n = size(x)

    nrm2 = dnrm2(n, x, incx)

end function nrm2

!------------------------------------------------------------------------------

function dot(x,y)

    real(dp), external :: ddot

    real(dp) :: dot
    real(dp), dimension(:), intent(in) :: x, y

    integer :: n, incx, incy

    incx = 1
    incy = 1

    n = size(x)

    dot = ddot(n, x, incx, y, incy)

end function dot

!------------------------------------------------------------------------------

subroutine axpy(x, y, a)

    external :: daxpy

    real(dp), intent(in), optional :: a
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(inout) :: y

    real(dp) :: o_a
    integer :: n, incx, incy

    if (present(a)) then
        o_a = a
    else
        o_a = 1.0d0
    end if

    incx = 1
    incy = 1

    n = size(x)

    call daxpy(n, o_a, x, incx, y, incy)

end subroutine axpy

!------------------------------------------------------------------------------

subroutine gemm(a, b, c, transa, transb, alpha, beta)

    external :: dgemm

    real(dp), intent(in), optional :: alpha, beta
    character(len=1), intent(in), optional :: transa, transb
    real(dp), dimension(:,:), intent(in) :: a, b
    real(dp), dimension(:,:) , intent(inout) :: c

    integer :: m, n, k, lda, ldb, ldc
    character(len=1) :: o_transa, o_transb
    real(dp) :: o_alpha, o_beta

    if (present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1.0d0
    end if

    if (present(beta)) then
        o_beta = beta
    else
        o_beta = 0.0d0
    end if

    if (present(transa)) then
        o_transa = transa
    else
        o_transa = 'N'
    end if

    if (present(transb)) then
        o_transb = transb
    else
        o_transb = 'N'
    end if

    if (o_transa == 'N') then
        k = size(a, 2)
    else
        k = size(a, 1)
    end if

    m = size(c, 1)
    n = size(c, 2)
    lda = max(1, size(a, 1))
    ldb = max(1, size(b, 1))
    ldc = max(1, size(c, 1))

    call dgemm(o_transa, o_transb, m, n, k, o_alpha, a, lda, b, ldb, o_beta, c, ldc)

end subroutine gemm

!------------------------------------------------------------------------------

subroutine spmv(ap, x, y, uplo, alpha, beta)

    external :: dspmv

    real(dp), dimension(:), intent(in) :: ap
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(inout) :: y
    character(len=1), intent(in), optional :: uplo
    real(dp), intent(in), optional :: alpha, beta

    integer :: n, incx, incy
    real(dp) :: o_alpha, o_beta
    character(len=1) :: o_uplo

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    if (present(alpha)) then
        o_alpha = alpha
    else
        o_alpha = 1.0d0
    end if

    if (present(beta)) then
        o_beta = beta
    else
        o_beta = 0.0d0
    end if

    incx = 1
    incy = 1
    n = size(x)

    call dspmv(o_uplo, n, o_alpha, ap, x, incx, o_beta, y, incy)

end subroutine spmv

!------------------------------------------------------------------------------

end module blas_f90
