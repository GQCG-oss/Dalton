!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
!------------------------------------------------------------------------------

subroutine Tk_integrals(inttype, Tk_ints, nnbas, ncomps, coord)

    ! Gen1Int API
    use gen1int_api
    use pe_precision

    implicit none

    intrinsic :: present

    integer, intent(in) :: nnbas, ncomps
    real(dp), dimension(3,1), intent(in) :: coord
    real(dp), dimension(nnbas,ncomps), intent(out) :: Tk_ints
    character(*), intent(in) :: inttype

    integer :: num_ao
    integer :: num_prop
    integer :: prop_sym
    integer :: io
    integer :: printlvl = 0
    integer :: ierr
    logical :: triangular
    logical :: symmetric
    type(one_prop_t) :: prop_operator
    type(matrix), dimension(:), allocatable :: intmats

    integer :: nbas
    integer :: i, j, k, x, y, z
    integer, dimension(3,ncomps) :: row2col
    real(dp), dimension(1) :: charge

    ! non-zero components for the operator, the first dimension is for bra and
    ! ket sub-shells, the last is the number of non-zero components, which should
    ! be 1 for non-relativistic calcualtions
    integer nnz_comp(2,1)

    Tk_ints = 0.0d0

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * ncomps) - 1.0d0)) - 1

    nbas = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * nnbas) - 1.0d0))

    if (mod(k,2) == 0) then
        charge = -1.0d0
    else if (mod(k,2) /= 0) then
        charge = 1.0d0
    end if

    if (inttype == 'es') then
        call OnePropCreate(prop_name=INT_POT_ENERGY,&
                           one_prop=prop_operator,  &
                           info_prop=ierr,          &
                           idx_nuclei=(/-1/),       &
                           coord_nuclei=coord,      &
                           charge_nuclei=charge,    &
                           order_geo_pot=k)
    else
        stop 'ERROR: unknown integral type'
    end if
    if (ierr /= 0) stop 'Failed to create property operator.'
    ! gets the number of property integrals and their symmetry
    call OnePropGetNumProp(one_prop=prop_operator, &
                           num_prop=num_prop)
    call OnePropGetSymmetry(one_prop=prop_operator, &
                            prop_sym=prop_sym)
    if (num_prop /= ncomps) stop 'Wrong number of components.'

    call Gen1IntAPIGetNumAO(num_ao=num_ao)
    if (num_ao /= nbas) stop 'Array size inconsistency.'

    select case(prop_sym)
        case(SYMM_INT_MAT, ANTI_INT_MAT)
            triangular = .true.
            symmetric = (prop_sym == SYMM_INT_MAT)
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
                if (x+y+z /= k) cycle
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

    ! sets the non-zero components for the one-electron operator, here we only
    ! consider the (large,large) part
    nnz_comp(1,1) = 1
    nnz_comp(2,1) = 1

    ! calculates the integrals, please refer to the comments in subroutine
    ! \fn(Gen1IntOnePropGetIntExpt) in gen1int_api.F90
    call Gen1IntOnePropGetIntExpt(nnz_comp=nnz_comp,      &
                                  one_prop=prop_operator, &
                                  num_ints=num_prop,      &
                                  val_ints=intmats,       &
                                  num_dens=1,             &
                                  io_viewer=io,           &
                                  level_print=printlvl)

    call OnePropDestroy(one_prop=prop_operator)

    do i = 1, num_prop
        call MatNullify(A=intmats(i))
    end do

    deallocate(intmats)

end subroutine Tk_integrals
