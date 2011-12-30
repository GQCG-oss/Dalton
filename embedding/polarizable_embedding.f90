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
    real(8), dimension(:,:), allocatable, save :: monopoles
    real(8), dimension(:,:), allocatable, save :: dipoles
    real(8), dimension(:,:), allocatable, save :: quadrupoles
    real(8), dimension(:,:), allocatable, save :: octopoles
    real(8), dimension(:,:), allocatable, save :: hexadecapoles

    ! (hyper)polarizabilities
    real(8), dimension(:,:), allocatable, save :: isopols
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
        allocate(monopoles(ncents,1))
        monopoles = 0.0d0
        do i = 1, nlines
            read(lupot,*) k
            read(lupot,*) monopoles(k,1)
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
        allocate(isopols(ncents,1))
        isopols = 0.0d0
        do i = 1, nlines
            read(lupot,*) k
            read(lupot,*) isopols(k,1)
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

    logical :: lexist
    integer :: i, j
    integer :: lutemp
    real(8) :: E_el, E_nuc
    real(8), dimension(0:3) :: E_Qk_el, E_Qk_nuc

    E_es = 0.0d0
    E_Qk_el = 0.0d0; E_Qk_nuc = 0.0d0

    inquire(file='pe_temp.bin', exist=lexist)
    if (lexist) then
        call openfile('pe_temp.bin', lutemp, 'old', 'unformatted')
        rewind(lutemp)
        read(lutemp) E_nuc, fock
        close(lutemp)
        E_el = dot_product(density, fock)
    else
        if (mulorder >= 0) then
            call pe_monopoles(density, fock, E_Qk_el(0), E_Qk_nuc(0),&
                              qm_coords, qm_charges, work)
        end if

        if (mulorder >= 1) then
            call pe_dipoles(density, fock, E_Qk_el(1), E_Qk_nuc(1),&
                            qm_coords, qm_charges, work)
        end if

        if (mulorder >= 2) then
            call pe_quadrupoles(density, fock, E_Qk_el(2), E_Qk_nuc(2),&
                                qm_coords, qm_charges, work)
        end if

        if (mulorder >= 3) then
            call pe_octopoles(density, fock, E_Qk_el(3), E_Qk_nuc(3),&
                              qm_coords, qm_charges, work)
        end if

        E_el = sum(E_Qk_el)
        E_nuc = sum(E_Qk_nuc)

        call openfile('pe_temp.bin', lutemp, 'new', 'unformatted')
        rewind(lutemp)
        write(lutemp) E_nuc, fock
        close(lutemp)
    end if

    E_es = E_nuc + E_el

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
    integer :: i, j, k
    real(8), dimension(3) :: Rij
    real(8), dimension(1) :: Tk
    real(8), dimension(:,:), allocatable :: Qk_ints

    ndim = size(fock)
    allocate(Qk_ints(ndim,1))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        ! Skip if charge is "zero"
        if (abs(monopoles(i,1)) < zero_thr) cycle

        ! get Q^(0) integrals
        call get_Qk_integrals(Qk_ints, 0, coords(i,:), monopoles(i,:), work)

        ! nuclei - monopole interaction
        do j = 1, size(qm_charges)

            ! distance vector between nucleus and monopole
            Rij = qm_coords(j,:) - coords(i,:)

            ! get the T^(0) tensor (which is a scalar)
            call get_Tk_tensor(Tk, 0, Rij)

            E_nuc = E_nuc + monopoles(i,1) * qm_charges(j) * Tk(1)

        end do

        ! electron - monopole interaction
        E_el = E_el + dot_product(density, Qk_ints(:,1))
        fock = fock + Qk_ints(:,1)

    end do

    deallocate(Qk_ints)

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
    integer :: i, j, k
    real(8), dimension(3) :: Rij, Tk
    real(8), dimension(:,:), allocatable :: Qk_ints

    ndim = size(fock)
    allocate(Qk_ints(ndim,3))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        ! Skip if dipole length is "zero"
        if (abs(length(dipoles(i,:))) < zero_thr) cycle

        ! get Q^(1) integrals
        call get_Qk_integrals(Qk_ints, 1, coords(i,:), dipoles(i,:), work)

        ! nuclei - dipole interaction energy
        do j = 1, size(qm_charges)

            ! distance vector between nucleus and dipole
            Rij = qm_coords(j,:) - coords(i,:)

            ! get the T^(1) tensor
            call get_Tk_tensor(Tk, 1, Rij)

            do k = 1, 3
                E_nuc = E_nuc + qm_charges(j) * dipoles(i,k) * Tk(k)
            end do
        end do

        ! electron - dipole interaction
        do j = 1, 3
            E_el = E_el + dot_product(density, Qk_ints(:,j))
            fock = fock + Qk_ints(:,j)
        end do

    end do

    deallocate(Qk_ints)

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
    integer :: i, j, k
    real(8), dimension(3) :: Rij
    real(8), dimension(6) :: Tk
    real(8), dimension(:,:), allocatable :: Qk_ints

    ndim = size(fock)
    allocate(Qk_ints(ndim,6))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        ! Skip if norm is "zero"
        if (abs(length(quadrupoles(i,:))) < zero_thr) cycle

        ! get Q^(2) integrals
        call get_Qk_integrals(Qk_ints, 2, coords(i,:), quadrupoles(i,:), work)

        ! nuclei - quadrupole interaction energy
        do j = 1, size(qm_charges)

            ! distance vector between nucleus and quadrupole
            Rij = qm_coords(j,:) - coords(i,:)

            ! get the T^(2) tensor
            call get_Tk_tensor(Tk, 2, Rij)

            do k = 1, 6
                E_nuc = E_nuc + 0.5d0 * qm_charges(j) * quadrupoles(i,k) * Tk(k)
            end do

        end do

        ! electron - quadrupole interaction energy
        do j = 1, 6
            fock = fock + Qk_ints(:,j)
            E_el = E_el + dot_product(density, Qk_ints(:,j))
        end do

    end do

    deallocate(Qk_ints)

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
    integer :: i, j, k
    real(8), dimension(3) :: Rij
    real(8), dimension(10) :: Tk
    real(8), dimension(:,:), allocatable :: Qk_ints

    ndim = size(fock)
    allocate(Qk_ints(ndim,10))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        ! Skip if norm is "zero"
        if (abs(length(octopoles(i,:))) < zero_thr) cycle

        ! get Q^(3) integrals
        call get_Qk_integrals(Qk_ints, 3, coords(i,:), octopoles(i,:), work)

        ! nuclei - octopole interaction energy
        do j = 1, size(qm_charges)

            ! distance vector between nucleus and octopole
            Rij = qm_coords(j,:) - coords(i,:)

            ! get the T^(3) tensor
            call get_Tk_tensor(Tk, 3, Rij)

            do k = 1, 10
                E_nuc = E_nuc + 6.0d0**(-1) * qm_charges(j) * octopoles(i,k) * Tk(k)
            end do

        end do

        ! electron - octopole interaction energy
        do j = 1, 10
            fock = fock + Qk_ints(:,j)
            E_el = E_el + dot_product(density, Qk_ints(:,j))
        end do

    end do

    deallocate(Qk_ints)

end subroutine pe_octopoles

!##############################################################################

subroutine get_Tk_tensor(Tk, k, Rij)

    integer, intent(in) :: k
    real(8), dimension(:), intent(out) :: Tk
    real(8), dimension(3), intent(in) :: Rij

    real(8) :: R, R3i, R5i, R7i

    R = length(Rij)

    if (k == 0) then
        Tk(1) = R**(-1)
    else if (k == 1) then
        R3i = R**(-3)
        ! Tx
        Tk(1) = Rij(1) * R3i
        ! Ty
        Tk(2) = Rij(2) * R3i
        ! Tz
        Tk(3) = Rij(3) * R3i
    else if (k == 2) then
        R3i = R**(-3)
        R5i = R**(-5)
        ! Txx
        Tk(1) = 3.0d0 * Rij(1)**2 * R5i - R3i
        ! Txy
        Tk(2) = 6.0d0 * Rij(1) * Rij(2) * R5i
        ! Txz
        Tk(3) = 6.0d0 * Rij(1) * Rij(3) * R5i
        ! Tyy
        Tk(4) = 3.0d0 * Rij(2)**2 * R5i - R3i
        ! Tyz
        Tk(5) = 6.0d0 * Rij(2) * Rij(3) * R5i
        ! Tzz
        Tk(6) = 3.0d0 * Rij(3)**2 * R5i - R3i
    else if (k == 3) then
        R7i = R**(-7)
        R5i = R**(-5)
        ! Txxx
        Tk(1) = 15.0d0 * Rij(1)**3 * R7i - 9.0d0 * Rij(1) * R5i
        ! Txxy
        Tk(2) = 45.0d0 * Rij(1)**2 * Rij(2) * R7i - 9.0d0 * Rij(2) * R5i
        ! Txxz
        Tk(3) = 45.0d0 * Rij(1)**2 * Rij(3) * R7i - 9.0d0 * Rij(3) * R5i
        ! Txyy
        Tk(4) = 45.0d0 * Rij(1) * Rij(2)**2 * R7i - 9.0d0 * Rij(1) * R5i
        ! Txyz
        Tk(5) = 90.0d0 * Rij(1) * Rij(2) * Rij(3) * R7i
        ! Txzz
        Tk(6) = 45.0d0 * Rij(1) * Rij(3)**2 * R7i - 9.0d0 * Rij(1) * R5i
        ! Tyyy
        Tk(7) = 15.0d0 * Rij(2)**3 * R7i - 9.0d0 * Rij(2) * R5i
        ! Tyyz
        Tk(8) = 45.0d0 * Rij(2)**2 * Rij(3) * R7i - 9.0d0 * Rij(3) * R5i
        ! Tyzz
        Tk(9) = 45.0d0 * Rij(2) * Rij(3)**2 * R7i - 9.0d0 * Rij(2) * R5i
        ! Tzzz
        Tk(10) = 45.0d0 * Rij(3)**3 * R7i - 9.0d0 * Rij(3) * R5i
    end if


end subroutine get_Tk_tensor

!##############################################################################

subroutine get_Qk_integrals(Qk_ints, k, coord, multipole, work)

    ! order of the multipole, i.e. Q^(k)
    integer, intent(in) :: k
    ! the Q^(k) integral matrix
    real(8), dimension(:,:), intent(out) :: Qk_ints
    ! the Q^(k) multipole moment
    real(8), dimension(:), intent(in) :: multipole
    ! xyz position of the multipole moment
    real(8), dimension(3), intent(in) :: coord
    ! work array
    real(8), dimension(:) :: work
    ! symmetry factors
    real(8), dimension(0:3,10) :: factors
    ! number of unique components of the multipole moment
    integer :: ncomps
    ! dimension of Q^(k) integral matrix
    integer :: ndim
    ! loop indices
    integer :: i, j

    ! symmetry factors and Taylor expansion factors
    factors = 0.0d0
    ! monopole
    factors(0,1) = 1.0d0
    ! dipole
    factors(1,1) = -1.0d0
    factors(1,2) = -1.0d0
    factors(1,3) = -1.0d0
    ! quadrupole (incl. 1/2 from Taylor expansion)
    factors(2,1) = 0.5d0
    factors(2,2) = 1.0d0
    factors(2,3) = 1.0d0
    factors(2,4) = 0.5d0
    factors(2,5) = 1.0d0
    factors(2,6) = 0.5d0
    ! octopole (incl. 1/6 from Taylor expansion)
    factors(3,1) = -1.0d0/6.0d0
    factors(3,2) = -0.5d0
    factors(3,3) = -0.5d0
    factors(3,4) = -0.5d0
    factors(3,5) = -1.0d0
    factors(3,6) = -0.5d0
    factors(3,7) = -1.0d0/6.0d0
    factors(3,8) = -0.5d0
    factors(3,9) = -0.5d0
    factors(3,10) = -1.0d0/6.0d0
    ! hexadecapole
!    factors(4,1) = 1.0d0
!    factors(4,2) = 4.0d0
!    factors(4,3) = 4.0d0
!    factors(4,4) = 6.0d0
!    factors(4,5) = 12.0d0
!    factors(4,6) = 6.0d0
!    factors(4,7) = 4.0d0
!    factors(4,8) = 12.0d0
!    factors(4,9) = 12.0d0
!    factors(4,10) = 4.0d0
!    factors(4,11) = 1.0d0
!    factors(4,12) = 4.0d0
!    factors(4,13) = 6.0d0
!    factors(4,14) = 4.0d0
!    factors(4,15) = 1.0d0

    ndim = size(Qk_ints, 1)
    ncomps = size(Qk_ints, 2)

    ! get T^(k) integrals (incl. negative sign)
    call get_Tk_integrals(Qk_ints, ncomps*ndim, k, 'packed', coord,&
                          work, size(work))

    ! dot T^(k) integrals with multipole to get Q^(k) integrals
    do i = 1, ncomps
        Qk_ints(:,i) = factors(k,i) * multipole(i) * Qk_ints(:,i)
    end do

end subroutine get_Qk_integrals

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
