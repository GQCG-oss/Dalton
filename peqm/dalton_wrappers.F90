#ifdef BUILD_GEN1INT
subroutine Tk_integrals(Tk_ints, nnbas, ncomps, R, work, nwrk)

    use dalton_shell

    implicit none

    intrinsic :: real

    integer, intent(in) :: nnbas, ncomps, nwrk
    real(8), dimension(3,1), intent(in) :: R
    real(8), dimension(nnbas,ncomps), intent(out) :: Tk_ints
    real(8), dimension(nwrk), intent(inout) :: work

    type(one_prop_t) :: prop_operator
    integer :: num_ao
    integer :: size_ints
    integer :: num_prop
    integer :: kind_prop
    logical :: triangular
    logical :: symmetric
    type(matrix), dimension(:), allocatable :: Int_Matrix
    integer :: io_viewer
    integer :: level_print = 0
    integer :: imat
    integer :: ierr

    integer :: i, j, k
    real(8), dimension(1) :: charge

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * real(ncomps,8)) - 1.0d0)) - 1

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

    call DaltonShellGetNumAO(num_ao=num_ao)

    select case(kind_prop)
        case(SYMM_INT_MAT, ANTI_INT_MAT)
            triangular = .true.
            symmetric = (kind_prop == SYMM_INT_MAT)
            size_ints = num_ao * (num_ao + 1) / 2
        case default
            triangular = .false.
            symmetric = .false.
            size_ints = num_ao * num_ao
    end select

    if (size_ints /= size(Tk_ints,1)) stop 'Integral matrix is wrong size.'

    allocate(Int_Matrix(num_prop), stat=ierr)
    if (ierr /= 0) stop 'Failed to allocate matrices.'

    if (k == 0) then
        call MatAssociate(Tk_ints(1:size_ints,1), num_ao, Int_Matrix(1), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices for k = 0.'
    else if (k == 1) then
        do i = 1, num_prop
            call MatAssociate(Tk_ints(1:size_ints,i), num_ao, Int_Matrix(i), ierr,&
                              triangular, symmetric)
            if (ierr /= 0) stop 'Failed to associate integral matrices.'
        end do
    else if (k == 2) then
        call MatAssociate(Tk_ints(1:size_ints,1), num_ao, Int_Matrix(1), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,2), num_ao, Int_Matrix(2), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,4), num_ao, Int_Matrix(3), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,3), num_ao, Int_Matrix(4), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,5), num_ao, Int_Matrix(5), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,6), num_ao, Int_Matrix(6), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
    else if (k == 3) then
        call MatAssociate(Tk_ints(1:size_ints,1), num_ao, Int_Matrix(1), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,2), num_ao, Int_Matrix(2), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,4), num_ao, Int_Matrix(3), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,7), num_ao, Int_Matrix(4), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,3), num_ao, Int_Matrix(5), ierr,&
                         triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,5), num_ao, Int_Matrix(6), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,8), num_ao, Int_Matrix(7), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,6), num_ao, Int_Matrix(8), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,9), num_ao, Int_Matrix(9), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
        call MatAssociate(Tk_ints(1:size_ints,10), num_ao, Int_Matrix(10), ierr,&
                          triangular, symmetric)
        if (ierr /= 0) stop 'Failed to associate integral matrices.'
    end if

    call DaltonShellIntegral(one_prop=prop_operator,&
                             num_ints=num_prop,     &
                             val_ints=Int_Matrix,     &
                             num_dens=1,            &
                             io_viewer=io_viewer,   &
                             level_print=level_print)

    call OnePropDestroy(one_prop=prop_operator)

    do i = 1, num_prop
        call MatNullify(A=Int_Matrix(i))
    end do

    deallocate(Int_Matrix)

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
        call get1in(Tk_ints(1), inttype, ncomps, work(1), nwrk,         &
                    labint, intrep, intadr, 0, .false., 0, trimat,      &
                    dummy, .false., dummy, 1)
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
            work(l*i+1:j*i) = (work(l*i+1:j*i) - work(j*i+1:(j+1)*i)) &
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
