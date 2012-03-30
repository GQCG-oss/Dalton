#ifdef BUILD_GEN1INT
subroutine Tk_integrals(Tk_ints, nnbas, ncomps, R, work, nwrk)

    use dalton_shell

    implicit none

    integer, intent(in) :: nnbas, ncomps, nwrk
    real(8), dimension(3,1), intent(in) :: R
    real(8), dimension(nnbas,ncomps), intent(out) :: Tk_ints
    real(8), dimension(nwrk), intent(inout) :: work

    type(one_prop_t) :: prop_operator
    integer :: num_ao
    integer :: num_prop
    integer :: kind_prop
    logical :: triangular
    logical :: symmetric
    type(matrix), dimension(:), allocatable :: intmats
    integer :: io
    integer :: printlvl = 0
    integer :: ierr

    integer :: nbas
    integer :: i, j, k, x, y, z
    integer, dimension(3,ncomps) :: row2col
    real(8), dimension(1) :: charge

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * ncomps) - 1.0d0)) - 1

    nbas = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * nnbas) - 1.0d0))

    if (mod(k,2) == 0) then
        charge = -1.0d0
    else if (mod(k,2) /= 0) then
        charge = 1.0d0
    end if

    call OnePropCreate(prop_name=INT_POT_ENERGY,&
                       one_prop=prop_operator,  &
                       info_prop=ierr,          &
                       num_prop=num_prop,       &
                       kind_prop=kind_prop,     &
                       idx_nuclei=(/-1/),       &
                       coord_nuclei=R,          &
                       charge_nuclei=charge,    &
                       order_geo_pot=k)
    if (ierr /= 0) stop 'Failed to create property operator.'
    if (num_prop /= ncomps) stop 'Wrong number of components.'

    call DaltonShellGetNumAO(num_ao=num_ao)
    if (num_ao /= nbas) stop 'Array size inconsistency.'

    select case(kind_prop)
        case(SYMM_INT_MAT, ANTI_INT_MAT)
            triangular = .true.
            symmetric = (kind_prop == SYMM_INT_MAT)
        case default
            triangular = .false.
            symmetric = .false.
            stop 'Integral matrices not symmetric!'
    end select

    allocate(intmats(num_prop), stat=ierr)
    if (ierr /= 0) stop 'Failed to allocate matrices.'

    i = 1
    do x = k, 0, -1
        do y = k, 0, -1
            do z = k, 0, -1
                if (x+y+z > k .or. x+y+z < k) cycle
                row2col(:,i) = (/ x, y, z /)
                i = i + 1
            end do
        end do
    end do

    i = 1
    do z = 0, k
        do y = 0, k
            do x = 0, k
                if (x+y+z > k .or. x+y+z < k) cycle
                do j = 1, ncomps
                    if (x == row2col(1,j) .and.&
                        y == row2col(2,j) .and.&
                        z == row2col(3,j)) then
                        call MatAssociate(work_alpha=Tk_ints(:,j),  &
                                          num_row=nbas,             &
                                          A=intmats(i),             &
                                          info_mat=ierr,            &
                                          triangular=triangular,    &
                                          symmetric=symmetric)
                        if (ierr /= 0) stop 'Failed to associate matrices.'
                    end if
                end do
                i = i + 1
            end do
        end do
    end do

    call DaltonShellIntegral(one_prop=prop_operator,&
                             num_ints=num_prop,     &
                             val_ints=intmats,      &
                             num_dens=1,            &
                             io_viewer=io,          &
                             level_print=printlvl)

    call OnePropDestroy(one_prop=prop_operator)

    do i = 1, num_prop
        call MatNullify(A=intmats(i))
    end do

    deallocate(intmats)

end subroutine Tk_integrals
#else
subroutine Tk_integrals(Tk_ints, nints, ncomps, coord, work, nwrk)

#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qm3.h"
#include "orgcom.h"

    character(len=8), dimension(9*mxcent) :: labint
    integer, dimension(9*mxcent) :: intrep, intadr

    integer, intent(in) :: nints, ncomps, nwrk
    real(8), dimension(nints*ncomps), intent(inout) :: Tk_ints
    real(8), dimension(3), intent(in) :: coord
    real(8), dimension(nwrk), intent(inout) :: work

    logical :: trimat
    integer :: i, j, k, l, m, n
    character(len=7) :: inttype
    real(8), dimension(3) :: backup

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * real(ncomps)) - 1.0d0)) - 1

    backup = diporg
    diporg = coord
    runqm3 = .true.

    trimat = .true.

    if (k == 0) then
        inttype = 'NPETES '
    else if (k == 1) then
        inttype = 'NEFIELD'
    else if (k == 2) then
        inttype = 'ELFGRDC'
    else if (k == 3) then
        if (nwrk < 24 * nints) then
            print *, 'Not enough work space for T^(3) integrals!'
        end if
        inttype = 'ELFGRDC'
    end if

    Tk_ints = 0.0d0

    if (k <= 2) then
        n = ncomps
        call get1in(Tk_ints(1), inttype, n, work(1), nwrk, labint, intrep,  &
                    intadr, 0, .false., 0, trimat, dummy, .false., dummy, 1)
!    else if (k == 2) then
!        n = 3
!        m = 3 * nints
!        ! xx, xy, xz
!        diporg(1) = coord(1) + 0.01d0
!        call get1in(work(1:m), 'NEFIELD', n, work(m+1), nwrk,               &
!                    labint, intrep, intadr, 0, .false., 0, trimat, dummy,   &
!                    .false., dummy, 1)
!        diporg(1) = coord(1) - 0.01d0
!        call get1in(work(1+m:2*m), 'NEFIELD', n, work(2*m+1), nwrk,         &
!                    labint, intrep, intadr, 0, .false., 0, trimat,          &
!                    dummy, .false., dummy, 1)
!        diporg = coord
!        Tk_ints(1:3*nints) = - (work(1:m) - work(1+m:2*m))                  &
!                                / (2.0d0 * 0.01d0)
! 
!        ! yy, yz
!        diporg(2) = diporg(2) + 0.01d0
!        call get1in(work(1:m), 'NEFIELD', n, work(m+1), nwrk,               &
!                    labint, intrep, intadr, 0, .false., 0, trimat,          &
!                    dummy, .false., dummy, 1)
!        diporg(2) = diporg(2) - 2.0d0 * 0.01d0
!        call get1in(work(1+m:2*m), 'NEFIELD', n, work(2*m+1), nwrk,         &
!                    labint, intrep, intadr, 0, .false., 0, trimat,          &
!                    dummy, .false., dummy, 1)
!        diporg = coord
!        Tk_ints(1+3*nints:5*nints) = - (work(1+nints:m)                     &
!                                        - work(1+4*nints:2*m))              &
!                                        / (2.0d0 * 0.01d0)
! 
!        ! zz
!        diporg(3) = diporg(3) + 0.01d0
!        call get1in(work(1:m), 'NEFIELD', n, work(m+1), nwrk,               &
!                    labint, intrep, intadr, 0, .false., 0, trimat,          &
!                    dummy, .false., dummy, 1)
!        diporg(3) = diporg(3) - 2.0d0 * 0.01d0
!        call get1in(work(1+m:2*m), 'NEFIELD', n, work(2*m+1),nwrk,          &
!                    labint, intrep, intadr, 0, .false., 0, trimat,          &
!                    dummy, .false., dummy, 1)
!        diporg = coord
!        Tk_ints(1+5*nints:6*nints) = - (work(1+2*nints:m)                   &
!                                        - work(1+5*nints:2*m))              &
!                                        / (2.0d0 * 0.01d0)
    else if (k == 3) then
        n = 6
        m = nints
        i = n * nints
        l = 0
        do j = 1, 3
            diporg(j) = diporg(j) + 0.01d0
            call get1in(work(l*i+1), inttype, n, work(j*i+1), nwrk,     &
                        labint, intrep, intadr, 0, .false., 0, trimat,  &
                        dummy, .false., dummy, 1)
            diporg(j) = diporg(j) - 2.0d0 * 0.01d0
            call get1in(work(j*i+1), inttype, n, work((j+1)*i+1), nwrk, &
                        labint, intrep, intadr, 0, .false., 0, trimat,  &
                        dummy, .false., dummy, 1)
            diporg(j) = coord(j)
            work(l*i+1:j*i) = (work(l*i+1:j*i) - work(j*i+1:(j+1)*i))   &
                              / (2.0d0 * 0.01d0)
            l = l + 1
        end do
        Tk_ints(1:m) = work(1:m)
        Tk_ints(m+1:2*m) = work(i+1:i+1+m)
        Tk_ints(2*m+1:3*m) = work(2*i+1:2*i+1+m)
        Tk_ints(3*m+1:4*m) = work(i+1+m:i+1+2*m)
        Tk_ints(4*m+1:5*m) = work(2*i+1+m:2*i+1+2*m)
        Tk_ints(5*m+1:6*m) = work(2*i+1+2*m:2*i+1+3*m)
        Tk_ints(6*m+1:7*m) = work(i+1+3*m:i+1+4*m)
        Tk_ints(7*m+1:8*m) = work(2*i+1+3*m:2*i+1+4*m)
        Tk_ints(8*m+1:9*m) = work(2*i+1+4*m:2*i+1+5*m)
        Tk_ints(9*m+1:10*m) = work(2*i+1+5*m:2*i+1+6*m)
        Tk_ints = -1.0d0 * Tk_ints
    end if

    runqm3 = .false.
    diporg = backup

end subroutine Tk_integrals
#endif
