module polarizable_embedding

    implicit none

    external get_Tk_integrals

    private

    public :: pe_read_input, pe_main

    integer, save :: ncents, lenexlst
    integer, save :: mulorder, l_isopols, polorder

    ! thresholds
    real(8), parameter :: zero_thr = 1.0d-10

    ! elements, coordinates and exclusion list
    real(8), dimension(:), allocatable, save :: elems
    real(8), dimension(:,:), allocatable, save :: coords
    real(8), dimension(:,:), allocatable, save :: exlist

    ! multipole moments
    real(8), dimension(:), allocatable, save :: monopoles
    real(8), dimension(:,:), allocatable, save :: dipoles
    real(8), dimension(:,:), allocatable, save :: quadrupoles
    real(8), dimension(:,:), allocatable, save :: octopoles
    real(8), dimension(:,:), allocatable, save :: hexadecapoles

    ! (hyper)polarizabilities
    real(8), dimension(:), allocatable, save :: isopols
    real(8), dimension(:,:), allocatable, save :: dipols
    real(8), dimension(:,:), allocatable, save :: dihypols

! TODO:
! memory management
! add error catching
! use optional?

contains

subroutine pe_read_input(work, nwrk)

    ! input parameters could be options given in dalton input
    ! so that cutoffs etc. are handled in here.

    integer :: nwrk
    real(8), dimension(nwrk) :: work

    integer :: i, j, k, l
    integer :: lupot, nlines

    print *, size(work), nwrk

    call openfile('POTENTIAL.INP', lupot, 'old', 'formatted')

    read(lupot,*) mulorder, l_isopols, polorder

    read(lupot,*) ncents

    allocate(elems(ncents), coords(ncents,3))
    elems = 0.0d0; coords = 0.0d0

    do i = 1, ncents
        read(lupot,*) elems(i), (coords(i,j), j = 1, 3)
    end do

    if (mulorder >= 0) then
        read(lupot,*) nlines
        allocate(monopoles(ncents))
        monopoles = 0.0d0
        do i = 1, nlines
            read(lupot,*) k
            read(lupot,*) monopoles(k)
        end do
    end if

    if (mulorder >= 1) then
        read(lupot,*) nlines
        allocate(dipoles(ncents,3))
        dipoles = 0.0d0
        do i = 1, nlines
            read(lupot,*) k
            read(lupot,*) (dipoles(k,j), j = 1, 3)
        end do
    end if

    if (mulorder >= 2) then
        read(lupot,*) nlines
        allocate(quadrupoles(ncents,6))
        quadrupoles = 0.0d0
        do i = 1, nlines
            read(lupot,*) k
            read(lupot,*) (quadrupoles(k,j), j = 1, 6)
        end do
    end if

    if (mulorder >= 3) then
        read(lupot,*) nlines
        allocate(octopoles(ncents,10))
        octopoles = 0.0d0
        do i = 1, nlines
            read(lupot,*) k
            read(lupot,*) (octopoles(k,j), j = 1, 10)
        end do
    end if

    if (mulorder >= 4) then
        read(lupot,*) nlines
        allocate(hexadecapoles(ncents,15))
        hexadecapoles = 0.0d0
        do i = 1, nlines
            read(lupot,*) k
            read(lupot,*) (hexadecapoles(k,j), j = 1, 15)
        end do
    end if

    if (l_isopols == 1) then
        read(lupot,*) nlines
        allocate(isopols(ncents))
        isopols = 0.0d0
        do i = 1, nlines
            read(lupot,*) k
            read(lupot,*) isopols(k)
        end do
    end if

    if (polorder >= 1) then
        read(lupot,*) nlines
        allocate(dipols(ncents,6))
        dipols = 0.0d0
        do i = 1, nlines
            read(lupot,*) k
            read(lupot,*) (dipols(k,j), j = 1, 6)
        end do
    end if

    if (polorder >= 2) then
        read(lupot,*) nlines
        allocate(dihypols(ncents,10))
        dihypols = 0.0d0
        do i = 1, nlines
            read(lupot,*) k
            read(lupot,*) (dihypols(k,j), j = 1, 10)
        end do
    end if

    if (polorder >= 1 .or. l_isopols == 1) then
        read(lupot,*) nlines, lenexlst
        allocate(exlist(ncents,lenexlst))
        exlist = 0.0d0
        do i = 1, nlines
            read(lupot,*) k
            read(lupot,*) (exlist(k,j), j = 1, lenexlst)
        end do
    end if

    close(lupot)

end subroutine pe_read_input

!##############################################################################

subroutine pe_main(density, fock, E_pe, qm_coords, qm_charges, work)

    real(8), intent(out) :: E_pe
    real(8), dimension(:,:), intent(in) :: qm_coords
    real(8), dimension(:), intent(in) :: qm_charges
    real(8), dimension(:), intent(in) :: density
    real(8), dimension(:), intent(out) :: fock
    real(8), dimension(:), intent(inout) :: work

    real(8) :: E_es, E_ind

    fock = 0.0d0

    call pe_electrostatic(density, fock, E_es, qm_coords, qm_charges, work)

    E_pe = E_pe + E_es

end subroutine pe_main

!##############################################################################

subroutine pe_electrostatic(density, fock, E_es, qm_coords, qm_charges, work)

    real(8), intent(out) :: E_es
    real(8), dimension(:,:), intent(in) :: qm_coords
    real(8), dimension(:), intent(in) :: qm_charges
    real(8), dimension(:), intent(in) :: density
    real(8), dimension(:), intent(inout) :: fock
    real(8), dimension(:), intent(inout) :: work

    real(8) :: E_el, E_nuc

    E_es = 0.0d0

    if (mulorder >= 0) then
        call pe_monopoles(density, fock, E_el, E_nuc, qm_coords, qm_charges, work)
        E_es = E_es + E_el + E_nuc
    end if

    if (mulorder >= 1) then
        call pe_dipoles(density, fock, E_el, E_nuc, qm_coords, qm_charges, work)
        E_es = E_es + E_el + E_nuc
    end if

    if (mulorder >= 2) then
        call pe_quadrupoles(density, fock, E_el, E_nuc, qm_coords, qm_charges, work)
        E_es = E_es + E_el + E_nuc
    end if

    if (mulorder >= 3) then
        call pe_octopoles(density, fock, E_el, E_nuc, qm_coords, qm_charges, work)
        E_es = E_es + E_el + E_nuc
    end if

end subroutine pe_electrostatic

!##############################################################################

subroutine pe_monopoles(density, fock, E_el, E_nuc, qm_coords, qm_charges, work)

    real(8), dimension(:), intent(in) :: density
    real(8), dimension(:), intent(inout) :: fock
    real(8), dimension(:,:), intent(in) :: qm_coords
    real(8), dimension(:), intent(in) :: qm_charges
    real(8), dimension(:), intent(inout) :: work
    real(8), intent(out) :: E_el, E_nuc

    integer :: ndim
    integer :: i, j, k, l
    real(8), dimension(3) :: Rij
    real(8), dimension(:), allocatable :: Tk_ints

    ndim = size(fock)
    allocate(Tk_ints(ndim))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        ! Skip if charge is "zero"
        if (abs(monopoles(i)) < zero_thr) cycle

        ! get T^(0) integrals
        call get_Tk_integrals(Tk_ints, ndim, 0, 'packed', coords(i,:),&
                              work, size(work))

        ! electron - monopole interaction
        Tk_ints = monopoles(i) * Tk_ints
        E_el = E_el + dot_product(density, Tk_ints)
        fock = fock + Tk_ints

        ! nuclei - monopole interaction
        do j = 1, size(qm_charges)
            Rij = 0.0d0
            Rij = qm_coords(j,:) - coords(i,:)
            E_nuc = E_nuc + (monopoles(i) * qm_charges(j)) / length(Rij)
        end do
    end do

    deallocate(Tk_ints)

    print *, 'monopole'
    print *, E_el, E_nuc, E_el + E_nuc

end subroutine pe_monopoles

!##############################################################################

subroutine pe_dipoles(density, fock, E_el, E_nuc, qm_coords, qm_charges, work)

    real(8), dimension(:), intent(in) :: density
    real(8), dimension(:), intent(inout) :: fock
    real(8), dimension(:,:), intent(in) :: qm_coords
    real(8), dimension(:), intent(in) :: qm_charges
    real(8), dimension(:), intent(inout) :: work
    real(8), intent(out) :: E_el, E_nuc

    integer :: ndim
    integer :: i, j, k, l
    real(8), dimension(3) :: Rij
    real(8), dimension(:), allocatable :: Tk_ints

    ndim = size(fock)
    allocate(Tk_ints(3*ndim))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        ! Skip if dipole length is "zero"
        if (abs(length(dipoles(i,:))) < zero_thr) cycle

        ! get T^(1) integrals
        call get_Tk_integrals(Tk_ints, 3*ndim, 1, 'packed', coords(i,:),&
                              work, size(work))

        ! electron - dipole interaction
        ! x
        Tk_ints(1:ndim) = dipoles(i,1) * Tk_ints(1:ndim)
        E_el = E_el - dot_product(density, Tk_ints(1:ndim))
        fock = fock - Tk_ints(1:ndim)
        ! y
        Tk_ints(ndim+1:2*ndim) = dipoles(i,2) * Tk_ints(ndim+1:2*ndim)
        E_el = E_el - dot_product(density, Tk_ints(ndim+1:2*ndim))
        fock = fock - Tk_ints(ndim+1:2*ndim)
        ! z
        Tk_ints(2*ndim+1:3*ndim) = dipoles(i,3) * Tk_ints(2*ndim+1:3*ndim)
        E_el = E_el - dot_product(density, Tk_ints(2*ndim+1:3*ndim))
        fock = fock - Tk_ints(2*ndim+1:3*ndim)

        ! nuclei - dipole interaction energy
        do j = 1, size(qm_charges)
            Rij = 0.0d0
            Rij = qm_coords(j,:) - coords(i,:)
            do k = 1, 3
                E_nuc = E_nuc + qm_charges(j) * dipoles(i,k) * Rij(k) / length(Rij)**3
            end do
        end do
    end do

    deallocate(Tk_ints)

    print *, 'dipole'
    print *, E_el, E_nuc, E_el + E_nuc

end subroutine pe_dipoles

!##############################################################################

subroutine pe_quadrupoles(density, fock, E_el, E_nuc, qm_coords, qm_charges, work)

    real(8), dimension(:), intent(in) :: density
    real(8), dimension(:), intent(inout) :: fock
    real(8), dimension(:,:), intent(in) :: qm_coords
    real(8), dimension(:), intent(in) :: qm_charges
    real(8), dimension(:), intent(inout) :: work
    real(8), intent(out) :: E_el, E_nuc

    integer :: ndim
    integer :: i, j, k, l
    real(8), dimension(3) :: Rij
    real(8), dimension(6) :: T
    real(8), dimension(:), allocatable :: Tk_ints

    ndim = size(fock)
    allocate(Tk_ints(6*ndim))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        ! Skip if norm is "zero"
        if (abs(length(quadrupoles(i,:))) < zero_thr) cycle

        ! get T^(2) integrals
        call get_Tk_integrals(Tk_ints, 6*ndim, 2, 'packed', coords(i,:),&
                              work, size(work))

        ! electron - quadrupole interaction
        !
        ! the off-diagonal tensor components are multiplied by two due to
        ! symmetry, e.g. xy == yx etc. and all elements by factor 1/2 from
        ! Taylor expansion
        !
        ! xx
        Tk_ints(0*ndim+1:1*ndim) = 0.5d0 * quadrupoles(i,1) * Tk_ints(0*ndim+1:1*ndim)
        ! xy
        Tk_ints(1*ndim+1:2*ndim) = 1.0d0 * quadrupoles(i,2) * Tk_ints(1*ndim+1:2*ndim)
        ! xz
        Tk_ints(2*ndim+1:3*ndim) = 1.0d0 * quadrupoles(i,3) * Tk_ints(2*ndim+1:3*ndim)
        ! yy
        Tk_ints(3*ndim+1:4*ndim) = 0.5d0 * quadrupoles(i,4) * Tk_ints(3*ndim+1:4*ndim)
        ! yz
        Tk_ints(4*ndim+1:5*ndim) = 1.0d0 * quadrupoles(i,5) * Tk_ints(4*ndim+1:5*ndim)
        ! zz
        Tk_ints(5*ndim+1:6*ndim) = 0.5d0 * quadrupoles(i,6) * Tk_ints(5*ndim+1:6*ndim)

        ! electron - quadrupole interaction energy
        l = 0
        do k = 1, 6
            fock = fock + Tk_ints(l*ndim+1:k*ndim)
            E_el = E_el + dot_product(density, Tk_ints(l*ndim+1:k*ndim))
            l = l + 1
        end do

        ! nuclei - quadrupole interaction energy
        do j = 1, size(qm_charges)

            ! distance vector betweem nucleus and quadrupole
            Rij = qm_coords(j,:) - coords(i,:)

            ! the T^(2) interaction tensor
            ! the off-diagonal tensor components are multiplied by two due to
            ! symmetry, e.g. xy == yx etc.
            ! Txx
            T(1) = (3.0d0 * Rij(1)**2 - length(Rij)**2) / length(Rij)**5
            ! Txy
            T(2) = 2.0d0 * (3.0d0 * Rij(1) * Rij(2)) / length(Rij)**5
            ! Txz
            T(3) = 2.0d0 * (3.0d0 * Rij(1) * Rij(3)) / length(Rij)**5
            ! Tyy
            T(4) = (3.0d0 * Rij(2)**2 - length(Rij)**2) / length(Rij)**5
            ! Tyz
            T(5) = 2.0d0 * (3.0d0 * Rij(2) * Rij(3)) / length(Rij)**5
            ! Tzz
            T(6) = (3.0d0 * Rij(3)**2 - length(Rij)**2) / length(Rij)**5
    
            do k = 1, 6
                E_nuc = E_nuc + 0.5d0 * qm_charges(j) * quadrupoles(i,k) * T(k)
            end do

        end do
    
    end do

    print *, 'quadrupole'
    print *, E_el, E_nuc, E_el + E_nuc

    deallocate(Tk_ints)

end subroutine pe_quadrupoles


!##############################################################################

subroutine pe_octopoles(density, fock, E_el, E_nuc, qm_coords, qm_charges, work)

    real(8), dimension(:), intent(in) :: density
    real(8), dimension(:), intent(inout) :: fock
    real(8), dimension(:,:), intent(in) :: qm_coords
    real(8), dimension(:), intent(in) :: qm_charges
    real(8), dimension(:), intent(inout) :: work
    real(8), intent(out) :: E_el, E_nuc

    integer :: ndim
    integer :: i, j, k, l
    real(8), dimension(3) :: Rij
    real(8), dimension(10) :: T
    real(8), dimension(:), allocatable :: Tk_ints

    ndim = size(fock)
    allocate(Tk_ints(10*ndim))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        ! Skip if norm is "zero"
        if (abs(length(octopoles(i,:))) < zero_thr) cycle

        ! get T^(3) integrals
        call get_Tk_integrals(Tk_ints, 10*ndim, 3, 'packed', coords(i,:),&
                              work, size(work))

        ! electron - octopole interaction
        !
        ! the off-diagonal tensor components are multiplied by a factor due to
        ! symmetry, e.g. xxy == xyx == yxx etc. and all elements by factor 1/6 from
        ! Taylor expansion
        !
        ! xxx
        Tk_ints(0*ndim+1:1*ndim) =         octopoles(i,1) * Tk_ints(0*ndim+1:1*ndim)
        ! xxy
        Tk_ints(1*ndim+1:2*ndim) = 3.0d0 * octopoles(i,2) * Tk_ints(1*ndim+1:2*ndim)
        ! xxz
        Tk_ints(2*ndim+1:3*ndim) = 3.0d0 * octopoles(i,3) * Tk_ints(2*ndim+1:3*ndim)
        ! xyy
        Tk_ints(3*ndim+1:4*ndim) = 3.0d0 * octopoles(i,4) * Tk_ints(3*ndim+1:4*ndim)
        ! xyz
        Tk_ints(4*ndim+1:5*ndim) = 6.0d0 * octopoles(i,5) * Tk_ints(4*ndim+1:5*ndim)
        ! xzz
        Tk_ints(5*ndim+1:6*ndim) = 3.0d0 * octopoles(i,6) * Tk_ints(5*ndim+1:6*ndim)
        ! yyy
        Tk_ints(6*ndim+1:7*ndim) =         octopoles(i,7) * Tk_ints(6*ndim+1:7*ndim)
        ! yyz
        Tk_ints(7*ndim+1:8*ndim) = 3.0d0 * octopoles(i,8) * Tk_ints(7*ndim+1:8*ndim)
        ! yzz
        Tk_ints(8*ndim+1:9*ndim) = 3.0d0 * octopoles(i,9) * Tk_ints(8*ndim+1:9*ndim)
        ! zzz
        Tk_ints(9*ndim+1:10*ndim) =        octopoles(i,10) * Tk_ints(9*ndim+1:10*ndim)

        Tk_ints = (1.0d0 / 6.0d0) * Tk_ints

        ! electron - octopole interaction energy
        l = 0
        do k = 1, 10
            fock = fock + Tk_ints(l*ndim+1:k*ndim)
            E_el = E_el + dot_product(density, Tk_ints(l*ndim+1:k*ndim))
            l = l + 1
        end do

        ! nuclei - octopole interaction energy
        do j = 1, size(qm_charges)

            ! distance vector betweem nucleus and octopole
            Rij = qm_coords(j,:) - coords(i,:)

            ! the T^(3) interaction tensor
            ! the off-diagonal tensor components are multiplied by two due to
            ! symmetry, e.g. xy == yx etc.
            ! Txxx
            T(1) = (15.0d0 * Rij(1)**3 - 3.0d0 * length(Rij)**2 * (3.0d0 * Rij(1))) / length(Rij)**7
            ! Txxy
            T(2) = 3.0d0 * (15.0d0 * Rij(1)**2 * Rij(2) - 3.0d0 * length(Rij)**2 * (Rij(2))) / length(Rij)**7
            ! Txxz
            T(3) = 3.0d0 * (15.0d0 * Rij(1)**2 * Rij(3) - 3.0d0 * length(Rij)**2 * (Rij(3))) / length(Rij)**7
            ! Txyy
            T(4) = 3.0d0 * (15.0d0 * Rij(1) * Rij(2)**2 - 3.0d0 * length(Rij)**2 * (Rij(1))) / length(Rij)**7
            ! Txyz
            T(5) = 6.0d0 * (15.0d0 * Rij(1) * Rij(2) * Rij(3)) / length(Rij)**7
            ! Txzz
            T(6) = 3.0d0 * (15.0d0 * Rij(1) * Rij(3)**2 - 3.0d0 * length(Rij)**2 * (Rij(1))) / length(Rij)**7
            ! Tyyy
            T(7) = (15.0d0 * Rij(2)**3 - 3.0d0 * length(Rij)**2 * (3.0d0 * Rij(2))) / length(Rij)**7
            ! Tyyz
            T(8) = 3.0d0 * (15.0d0 * Rij(2)**2 * Rij(3) - 3.0d0 * length(Rij)**2 * (Rij(3))) / length(Rij)**7
            ! Tyzz
            T(9) = 3.0d0 * (15.0d0 * Rij(2) * Rij(3)**2 - 3.0d0 * length(Rij)**2 * (Rij(2))) / length(Rij)**7
            ! Tzzz
            T(10) = (15.0d0 * Rij(3)**3 - 3.0d0 * length(Rij)**2 * (3.0d0 * Rij(3))) / length(Rij)**7

            do k = 1, 10
                E_nuc = E_nuc + (1.0d0 / 6.0d0) * qm_charges(j) * octopoles(i,k) * T(k)
            end do

        end do

    end do

    print *, 'octopole'
    print *, E_el, E_nuc, E_el + E_nuc

    deallocate(Tk_ints)

end subroutine pe_octopoles

!##############################################################################

!##############################################################################

subroutine openfile(filename, lunit, stat, frmt)

    character(*), intent(in) :: filename, stat, frmt
    integer, intent(out) :: lunit
    integer :: i
    logical :: lexist, lopen

    if (stat == 'old') then
      inquire(file=filename, exist=lexist)

      if (.not. lexist) then
        print *, filename, ' not found!'
        stop
      endif
    endif

    do i = 10, 99
      inquire(unit=i, opened=lopen)
      if (lopen) then
        cycle
      else
        lunit = i
        open(unit=lunit, file=filename, status=stat, form=frmt)
        exit
      endif
    enddo

    return

end subroutine openfile

!##############################################################################

function length(vector)

    integer :: i, j
    real(8) :: length
    real(8), dimension(:), intent(in) :: vector

    length = 0.0d0
    do i = 1, size(vector)
        length = length + vector(i)**2
    end do

    length = sqrt(length)

end function length

end module polarizable_embedding
