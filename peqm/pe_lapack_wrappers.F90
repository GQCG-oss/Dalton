module pe_lapack_wrappers

    use pe_precision

    implicit none

contains

!------------------------------------------------------------------------------

subroutine sptrf(ap, uplo, ipiv, info)

    external :: dsptrf

    real(dp), dimension(:), intent(inout) :: ap
    character(len=1), optional :: uplo
    integer, dimension(:), optional, target :: ipiv
    integer, intent(out), optional :: info

    character(len=1) :: o_uplo
    integer :: o_info
    integer :: n, nn
    integer, dimension(:), pointer :: o_ipiv

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1 + sqrt(1.0d0 + 8.0d0 * real(nn, dp))) * 0.5d0)
    end if

    if (present(ipiv)) then
        o_ipiv => ipiv
    else
        allocate(o_ipiv(n))
    end if

    call dsptrf(o_uplo, n, ap, o_ipiv, o_info)

    if (.not. present(ipiv)) then
        deallocate(o_ipiv)
    end if

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('SPTRF', - o_info)
    end if

end subroutine sptrf

!------------------------------------------------------------------------------

subroutine sptri(ap, ipiv, uplo, info)

    external :: dsptri

    real(dp), dimension(:), intent(inout) :: ap
    integer, dimension(:), intent(in) :: ipiv
    character(len=1), intent(in), optional :: uplo
    integer, intent(out), optional :: info

    character(len=1) :: o_uplo
    integer :: o_info
    integer :: n, nn
    real(dp), dimension(:), allocatable :: work

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1 + sqrt(1.0d0 + 8.0d0 * real(nn, dp))) * 0.5d0)
    end if

    allocate(work(n))

    call dsptri(o_uplo, n, ap, ipiv, work, o_info)

    deallocate(work)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('SPTRI', - o_info)
    end if

end subroutine sptri

!------------------------------------------------------------------------------

subroutine sptrs(ap, b, ipiv, uplo, info)

    external :: dsptrs

    real(dp), dimension(:), intent(in) :: ap
    real(dp), dimension(:,:), intent(inout) :: b
    integer, dimension(:), intent(in) :: ipiv
    character(len=1), optional :: uplo
    integer, optional :: info

    character(len=1) :: o_uplo
    integer :: o_info
    integer :: n, nn, nrhs, ldb

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)
    nrhs = size(b, 2)
    ldb = max(1, size(b, 1))

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1 + sqrt(1.0d0 + 8.0d0 * real(nn, dp))) * 0.5d0)
    end if

    call dsptrs(o_uplo, n, nrhs, ap, ipiv, b, ldb, o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('SPTRS', - o_info)
    end if

end subroutine sptrs

!------------------------------------------------------------------------------

subroutine pptrf(ap, uplo, info)

    external :: dpptrf

    real(dp), dimension(:), intent(inout) :: ap
    character(len=1), optional :: uplo
    integer, optional :: info

    character(len=1) :: o_uplo
    integer :: o_info
    integer :: n, nn

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1 + sqrt(1.0d0 + 8.0d0 * real(nn, dp))) * 0.5d0)
    end if

    call dpptrf(o_uplo, n, ap, o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('PPTRF', - o_info)
    end if

end subroutine pptrf

!------------------------------------------------------------------------------

subroutine pptri(ap, uplo, info)

    external :: dpptri

    real(dp), dimension(:), intent(inout) :: ap
    character(len=1), optional :: uplo
    integer, optional :: info

    character(len=1) :: o_uplo
    integer :: o_info
    integer :: n, nn

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn 
    else
        n = int((- 1 + sqrt(1.0d0 + 8.0d0 * real(nn, dp))) * 0.5d0)
    end if

    call dpptri(o_uplo, n, ap, o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('PPTRI', - o_info)
    end if

end subroutine pptri

!------------------------------------------------------------------------------

subroutine pptrs(ap, b, uplo, info)

    external :: dpptrs

    real(dp), dimension(:), intent(inout) :: ap
    real(dp), dimension(:,:), intent(inout) :: b
    character(len=1), optional :: uplo
    integer, optional :: info

    character(len=1) :: o_uplo
    integer :: o_info
    integer :: n, nn, nrhs, ldb

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)
    nrhs = size(b, 2)
    ldb = max(1, size(b, 1))

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1 + sqrt(1.0d0 + 8.0d0 * real(nn, dp))) * 0.5d0)
    end if

    call dpptrs(o_uplo, n, nrhs, ap, b, ldb, o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('PPTRS', - o_info)
    end if

end subroutine pptrs

!------------------------------------------------------------------------------

subroutine spcon(ap, ipiv, anorm, rcond, uplo, info)

    external :: dspcon

    real(dp), intent(in) :: anorm
    real(dp), intent(out) :: rcond
    real(dp), dimension(:), intent(in) :: ap
    integer, dimension(:), intent(in) :: ipiv
    character(len=1), optional :: uplo
    integer, intent(out), optional :: info

    character(len=1) :: o_uplo
    integer :: o_info
    integer :: n, nn
    real(dp), dimension(:), allocatable :: iwork, work

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1 + sqrt(1.0d0 + 8.0d0 * real(nn, dp))) * 0.5d0)
    end if

    allocate(iwork(n))
    allocate(work(2*n))

    call dspcon(o_uplo, n, ap, ipiv, anorm, rcond, work, iwork, o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('DSPCON', - o_info)
    end if

end subroutine spcon

!------------------------------------------------------------------------------

subroutine ppcon(ap, anorm, rcond, uplo, info)

    external :: dppcon

    real(dp), intent(in) :: anorm
    real(dp), intent(out) :: rcond
    real(dp), dimension(:), intent(in) :: ap
    character(len=1), optional :: uplo
    integer, intent(out), optional :: info

    character(len=1) :: o_uplo
    integer :: o_info
    integer :: n, nn
    real(dp), dimension(:), allocatable :: iwork, work

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1 + sqrt(1.0d0 + 8.0d0 * real(nn, dp))) * 0.5d0)
    end if

    allocate(iwork(n))
    allocate(work(3*n))

    call dppcon(o_uplo, n, ap, anorm, rcond, work, iwork, o_info)

    if (present(info)) then
        info = o_info
    else if (o_info <= - 1000) then
        call xerbla('DPPCON', - o_info)
    end if

end subroutine ppcon

!------------------------------------------------------------------------------

function lansp(norm, ap, uplo)

    real(dp), external :: dlansp

    real(dp) :: lansp
    real(dp), dimension(:), intent(in) :: ap
    character(len=1), intent(in) :: norm
    character(len=1), optional :: uplo

    character(len=1) :: o_uplo
    integer :: n, nn
    real(dp), dimension(:), allocatable :: work

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

    nn = size(ap)

    if (nn <= 0) then
        n = nn
    else
        n = int((- 1 + sqrt(1.0d0 + 8.0d0 * real(nn, dp))) * 0.5d0)
    end if

    allocate(work(n))

    lansp = dlansp(norm, o_uplo, n, ap, work)
    lansp = 0.0d0

    deallocate(work)

end function lansp

!------------------------------------------------------------------------------

end module pe_lapack_wrappers
