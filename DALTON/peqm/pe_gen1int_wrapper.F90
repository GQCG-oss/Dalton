subroutine Tk_integrals(inttype, Tk_ints, nnbas, ncomps, coord)

    ! Gen1Int API
    use gen1int_api
    use gen1int_matrix
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
    !N-ary tree for partial geometric derivatives on bra center
    type(nary_tree_t) :: nary_tree_bra
    !N-ary tree for partial geometric derivatives on ket center
    type(nary_tree_t) :: nary_tree_ket
    !N-ary tree for total geometric derivatives
    type(nary_tree_t) :: nary_tree_total
    !number of partial geometric derivatives on bra center
    integer :: num_geo_bra
    !number of partial geometric derivatives on ket center
    integer :: num_geo_ket
    !number of total geometric derivatives
    integer :: num_geo_total
    type(matrix), dimension(:), allocatable :: intmats

    integer :: nbas
    integer :: i, j, k, x, y, z
    integer, dimension(3,ncomps) :: row2col
    real(dp), dimension(1) :: charge
    real(dp), dimension(:,:), allocatable :: temp

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
    else if (inttype == 'gaussian') then
!        call OnePropCreate(prop_name=INT_GAUSSIAN_POT,&
!                           one_prop=prop_operator,    &
!                           info_prop=ierr,            &
!                           idx_gauorg=(/-1/),         &
!                           gaupot_origin=coord,       &
!                           gaupot_charge=charge,      &
!                           gaupot_expt=gauexp,        &
!                           order_geo_pot=k)
    else if (inttype == 'molgrad') then
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

    ! creates N-ary tree for geometric derivatives
    if (inttype == 'molgrad') then  !for geometry opt
    call Gen1IntAPINaryTreeCreate(max_num_cent=1,      &
                                  order_geo=1,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeCreate(max_num_cent=1,      &
                                  order_geo=1,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_total)


    else !if not geometry opt
    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_total)

    end if !NANNA: this end if should be moved further down!


    ! gets the number of property integrals and their symmetry
    call OnePropGetNumProp(one_prop=prop_operator, &
                           num_prop=num_prop)
    call OnePropGetSymmetry(one_prop=prop_operator, &
                            prop_sym=prop_sym)

    ! gets the number of geometric derivatives
    call NaryTreeGetNumGeo(nary_tree=nary_tree_bra, num_unique_geo=num_geo_bra)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_ket, num_unique_geo=num_geo_ket)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_total, num_unique_geo=num_geo_total)

    ! updates the number and symmetry of property integrals
    num_prop = num_prop*num_geo_bra*num_geo_ket*num_geo_total
    ! FIXME: if there are partial geometric derivatives, please set prop_sym = SQUARE_INT_MAT

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

#if defined(BUILD_OPENRSP)
    do i = 1, num_prop
            call MatCreate(A=intmats(i), num_row=nbas, info_mat=ierr,&
                          & triangular=triangular, symmetric=symmetric)
    end do
#else
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
#endif

    ! sets the non-zero components for the one-electron operator, here we only
    ! consider the (large,large) part
    nnz_comp(1,1) = 1
    nnz_comp(2,1) = 1

    ! calculates the integrals, please refer to the comments in subroutine
    ! \fn(Gen1IntOnePropGetIntExpt) in gen1int_api.F90
    call Gen1IntOnePropGetIntExpt(nnz_comp=nnz_comp,               &
                                  one_prop=prop_operator,          &
                                  nary_tree_bra=nary_tree_bra,     &
                                  nary_tree_ket=nary_tree_ket,     &
                                  nary_tree_total=nary_tree_total, &
                                  num_ints=num_prop,               &
                                  val_ints=intmats,                &
                                  num_dens=1,                      &
                                  io_viewer=io,                    &
                                  level_print=printlvl)

    call OnePropDestroy(one_prop=prop_operator)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_total)

#if defined(BUILD_OPENRSP)
    allocate(temp(nbas,nbas))
    temp = 0.0d0
    i = 1
    do z = 0, k
        do y = 0, k
            do x = 0, k
                if (x+y+z /= k) cycle
                do j = 1, ncomps
                    if (x == row2col(1,j) .and.&
                        y == row2col(2,j) .and.&
                        z == row2col(3,j)) then
                        call MatGetValues(intmats(i), 1, nbas, 1, nbas, temp)
                        call dgetsp(nbas, temp, Tk_ints(:,j))
                    end if
                end do
                i = i + 1
            end do
        end do
    end do
    do i = 1, num_prop
        call MatDestroy(A=intmats(i))
    end do
#else
    do i = 1, num_prop
        call MatNullify(A=intmats(i))
    end do
#endif

    deallocate(intmats)

end subroutine Tk_integrals


subroutine Tkb_integrals(Tk_ints, nnbas, ncomps, coord)

    ! Gen1Int API
    use gen1int_api
    use pe_precision

    implicit none

    intrinsic :: present

    integer, intent(in) :: nnbas, ncomps
    real(dp), dimension(3,1), intent(in) :: coord
    real(dp), dimension(nnbas,ncomps), intent(out) :: Tk_ints
!    logical, intent(in) :: gauss
!    real(dp), dimension(1), intent(in) :: gauexp

    integer :: num_ao
    integer :: num_prop
    integer :: prop_sym
    integer :: io
    integer :: printlvl = 0
    integer :: ierr
    logical :: triangular
    logical :: symmetric
    type(one_prop_t) :: prop_operator
    type(nary_tree_t) :: nary_tree_bra
    type(nary_tree_t) :: nary_tree_ket
    type(nary_tree_t) :: nary_tree_total
    integer :: num_geo_bra
    integer :: num_geo_ket
    integer :: num_geo_total
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

    ! for magnetic properties, the true number of components are
    ! multiplies by 3 (Bx, By, Bz)

    !nbas = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * nnbas) - 1.0d0))
    nbas = int(sqrt(1.0d0 * nnbas))

    if (mod(k,2) == 0) then
        charge = -1.0d0
    else if (mod(k,2) /= 0) then
        charge = 1.0d0
    end if

    call OnePropCreate(prop_name=INT_POT_ENERGY,&
                       one_prop=prop_operator,  &
                       info_prop=ierr,          &
                       idx_nuclei=(/-1/),       &
                       coord_nuclei=coord,      &
                       charge_nuclei=charge,    &
                       order_geo_pot=0)

    if (ierr /= 0) stop 'Failed to create property operator.'

    ! setup different properties
    call OnePropSetMag(one_prop=prop_operator, &
                       order_mag=1)
    call OnePropSetGTO(one_prop=prop_operator, &
                       gto_type=LONDON,        &
                       info_prop=ierr)
    if (ierr /= 0) stop 'Failed to set london orbitals.'

    ! creates N-ary tree for geometric derivatives on bra center
    ! here is zeroth order geometric derivatives
    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_bra)
    ! creates N-ary tree for geometric derivatives on ket center
    ! here is zeroth order geometric derivatives
    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_ket)
    ! creates N-ary tree for total geometric derivatives
    ! here is zeroth order geometric derivatives
    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_total)
    ! gets the number of geometric derivatives, since we are asking zeroth order geometric
    ! derivatives, the number is simply 1; otherwise, we need to
    ! call NaryTreeGetNumGeo(nary_tree=nary_tree_bra, num_unique_geo=num_geo_bra)
    ! call NaryTreeGetNumGeo(nary_tree=nary_tree_ket, num_unique_geo=num_geo_ket)
    ! call NaryTreeGetNumGeo(nary_tree=nary_tree_total, num_unique_geo=num_geo_total)
    num_geo_bra = 1
    num_geo_ket = 1
    num_geo_total = 1

    ! gets the number of property integrals and their symmetry
    call OnePropGetNumProp(one_prop=prop_operator, &
                           num_prop=num_prop)
    num_prop = num_prop*num_geo_bra*num_geo_ket*num_geo_total
    call OnePropGetSymmetry(one_prop=prop_operator, &
                            prop_sym=prop_sym)
    call OnePropView(one_prop=prop_operator, io_viewer=6)
    write(*,*) "CSS:NUMPROP", num_prop
    ! Gao: not sure what ncomps mean? and what order of geometric derivatives you may need?
    if (num_prop /= 3*ncomps) stop 'Wrong number of components.'

    call Gen1IntAPIGetNumAO(num_ao=num_ao)
    write(6,'(A5,3I5)') 'CSS:B', num_ao, nbas, nnbas
    if (num_ao /= nbas) stop 'Array size inconsistency.'

    ! if you have partial gemetric derivatives on bra and/or ket center, please
    ! reset prop_sym = SQUARE_INT_MAT
    write(6,'(A5,I5)') 'CSS:P', prop_sym
    select case(prop_sym)
        case(SYMM_INT_MAT, ANTI_INT_MAT)
            triangular = .true.
            symmetric = (prop_sym == SYMM_INT_MAT)
        case default
            triangular = .false.
            symmetric = .false.
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
                do j = 1, 3*ncomps
                    if (x == row2col(1,j) .and.&
                        y == row2col(2,j) .and.&
                        z == row2col(3,j)) then
                        call MatAssociate(work_alpha=Tk_ints(:,j),  &
                                          num_row=nbas,             &
                                          A=intmats(i),             &
                                          info_mat=ierr,            &
                                          triangular=triangular,    &
                                          symmetric=symmetric)
                        if (ierr /= 0) write(6,*) "nerr:",ierr,nbas,size(Tk_ints(:,j))
                        if (ierr /= 0) write(6,*) "nerr:",triangular, symmetric
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
    call Gen1IntOnePropGetIntExpt(nnz_comp=nnz_comp,               &
                                  one_prop=prop_operator,          &
                                  nary_tree_bra=nary_tree_bra,     &
                                  nary_tree_ket=nary_tree_ket,     &
                                  nary_tree_total=nary_tree_total, &
                                  num_ints=num_prop,               &
                                  val_ints=intmats,                &
                                  num_dens=1,                      &
                                  io_viewer=io,                    &
                                  level_print=printlvl)

    call OnePropDestroy(one_prop=prop_operator)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_total)

    do i = 1, num_prop
        call MatNullify(A=intmats(i))
    end do

    deallocate(intmats)

end subroutine Tkb_integrals
