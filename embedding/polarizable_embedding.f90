module polarizable_embedding

    implicit none

    private

    public :: pe_read_input, pe_main

    integer, parameter :: r8 = selected_real_kind(15)
    integer, parameter :: i8 = selected_int_kind(18)

    ! number of centers and length of exclusion list
    integer, save :: ncents, nexlist
    ! order of multipole moment expansion
    integer, save :: mulorder
    ! type of polarization, i.e. 0: isotropic alpha, 1: alpha, 2: A, 3: C
    integer, save :: poltype
    ! nonlinear polarizabilities
    integer, save :: hypoltype

    ! thresholds
    real(r8), parameter :: zero = 1.0d-8

    ! elements, coordinates and exclusion list
    real(r8), dimension(:), allocatable, save :: elems
    real(r8), dimension(:,:), allocatable, save :: Rs
    real(r8), dimension(:,:), allocatable, save :: exlist

    ! multipole moments
    ! monopoles
    real(r8), dimension(:,:), allocatable, save :: Q0
    ! dipoles
    real(r8), dimension(:,:), allocatable, save :: Q1
    ! quadrupoles
    real(r8), dimension(:,:), allocatable, save :: Q2
    ! octopoles
    real(r8), dimension(:,:), allocatable, save :: Q3
    ! hexadecapoles
    real(r8), dimension(:,:), allocatable, save :: Q4

    ! (hyper)polarizabilities
    ! dipole-dipole polarizability
    real(r8), dimension(:,:), allocatable, save :: alpha
    ! dipole-dipole-dipole polarizability / 1st hyperpolarizability
    real(r8), dimension(:,:), allocatable, save :: beta
    ! dipole-quadrupole polarizability
    real(r8), dimension(:,:), allocatable, save :: A
    ! quadrupole-quadrupole polarizability
    real(r8), dimension(:,:), allocatable, save :: C

    ! QM info: number of nuclei, nuclear charges and nuclear coordinates
    integer, save :: qmnucs
    real(r8), dimension(:), allocatable, save :: Zm
    real(r8), dimension(:,:), allocatable, save :: Rm

! TODO:
! cutoffs and damping
! memory management
! add error catching
! use optional?
! ddot or dot_product and dgemm or matmul
! parallelization (openMP, MPI, CUDA/openCL)

contains

subroutine pe_read_input(cord, charge, natoms, work, nwrk)

    ! input parameters could be options given in dalton input
    ! so that cutoffs etc. are handled in here.

    integer :: natoms, nwrk
    real(r8), dimension(natoms) :: charge
    real(r8), dimension(3,natoms) :: cord
    real(r8), dimension(nwrk) :: work

    integer :: i, j, s
    integer :: lupot, nlines

    qmnucs = natoms
    allocate(Rm(qmnucs,3), Zm(qmnucs))
    do i = 1, natoms
        Zm(i) = charge(i)
        do j = 1, 3
            Rm(i,j) = cord(j,i)
        end do
    end do

    call openfile('POTENTIAL.INP', lupot, 'old', 'formatted')

    read(lupot,*) mulorder, poltype, hypoltype
    read(lupot,*) ncents

    allocate(elems(ncents), Rs(ncents,3))
    elems = 0.0d0; Rs = 0.0d0

    do i = 1, ncents
        read(lupot,*) elems(i), (Rs(i,j), j = 1, 3)
    end do

    if (mulorder >= 0) then
        read(lupot,*) nlines
        allocate(Q0(ncents,1))
        Q0 = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) Q0(s,1)
        end do
    end if

    if (mulorder >= 1) then
        read(lupot,*) nlines
        allocate(Q1(ncents,3))
        Q1 = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (Q1(s,j), j = 1, 3)
        end do
    end if

    if (mulorder >= 2) then
        read(lupot,*) nlines
        allocate(Q2(ncents,6))
        Q2 = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (Q2(s,j), j = 1, 6)
        end do
    end if

    if (mulorder >= 3) then
        read(lupot,*) nlines
        allocate(Q3(ncents,10))
        Q3 = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (Q3(s,j), j = 1, 10)
        end do
    end if

    if (mulorder >= 4) then
        read(lupot,*) nlines
        allocate(Q4(ncents,15))
        Q4 = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (Q4(s,j), j = 1, 15)
        end do
    end if

    if (poltype >= 0) then
        allocate(alpha(ncents,6))
        alpha = 0.0d0
        read(lupot,*) nlines
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) alpha(s,1)
            alpha(s,4) = alpha(s,1)
            alpha(s,6) = alpha(s,1)
        end do
    end if

    if (poltype >= 1) then
        if (.not. allocated(alpha)) then
            allocate(alpha(ncents,6))
            alpha = 0.0d0
        end if
        read(lupot,*) nlines
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (alpha(s,j), j = 1, 6)
        end do
    end if

    if (poltype >= 2) then
        allocate(A(ncents,10))
        A = 0.0d0
        read(lupot,*) nlines
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (A(s,j), j = 1, 10)
        end do
    end if

    if (poltype >= 3) then
        allocate(C(ncents,15))
        C = 0.0d0
        read(lupot,*) nlines
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (C(s,j), j = 1, 15)
        end do
    end if

    if (poltype >= 1) then
        allocate(beta(ncents,10))
        beta = 0.0d0
        read(lupot,*) nlines
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (beta(s,j), j = 1, 10)
        end do
    end if

    if (poltype >= 1) then
        read(lupot,*) nlines, nexlist
        allocate(exlist(ncents,nexlist))
        exlist = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (exlist(s,j), j = 1, nexlist)
        end do
    end if

    close(lupot)

    ! check memory requirements?

end subroutine pe_read_input

!------------------------------------------------------------------------------

subroutine pe_main(density, fock, E_pe, work)

    real(r8), intent(out) :: E_pe
    real(r8), dimension(:), intent(in) :: density
    real(r8), dimension(:), intent(out) :: fock
    real(r8), dimension(:), intent(inout) :: work

    real(r8) :: E_es, E_ind

    fock = 0.0d0; E_pe = 0.0d0

    if (mulorder >= 0) then
        call pe_electrostatic(density, fock, E_es, work)
        E_pe = E_pe + E_es
    end if

    if (poltype >= 0) then
        call pe_polarization(density, fock, E_ind, work)
        E_pe = E_pe + E_ind
    end if

end subroutine pe_main

!------------------------------------------------------------------------------

subroutine pe_electrostatic(density, fock, E_es, work)

    real(r8), intent(out) :: E_es
    real(r8), dimension(:), intent(in) :: density
    real(r8), dimension(:), intent(inout) :: fock, work

    logical :: lexist
    integer :: lutemp
    real(r8), dimension(0:3) :: E_el, E_nuc

    E_es = 0.0d0; E_el = 0.0d0; E_nuc = 0.0d0

    inquire(file='pe_electrostatics.bin', exist=lexist)

    if (lexist) then
        call openfile('pe_electrostatics.bin', lutemp, 'old', 'unformatted')
        rewind(lutemp)
        read(lutemp) E_es, fock
        close(lutemp)
        E_es = E_es + dot(density, fock)
    else
        print *, 'q0'
        if (mulorder >= 0) then
            call es_monopoles(density, fock, E_el(0), E_nuc(0), work)
        end if
        print *, 'q1'
        if (mulorder >= 1) then
            call es_dipoles(density, fock, E_el(1), E_nuc(1), work)
        end if
        print *, 'q2'
        if (mulorder >= 2) then
            call es_quadrupoles(density, fock, E_el(2), E_nuc(2), work)
        end if
        print *, 'q3'
        if (mulorder >= 3) then
            call es_octopoles(density, fock, E_el(3), E_nuc(3), work)
        end if

        E_es = sum(E_el) + sum(E_nuc)

        call openfile('pe_electrostatics.bin', lutemp, 'new', 'unformatted')
        rewind(lutemp)
        write(lutemp) sum(E_nuc), fock
        close(lutemp)
    end if

end subroutine pe_electrostatic

!------------------------------------------------------------------------------

subroutine es_monopoles(density, fock, E_el, E_nuc, work)

    real(r8), dimension(:), intent(in) :: density
    real(r8), dimension(:), intent(inout) :: fock, work
    real(r8), intent(out) :: E_el, E_nuc

    integer :: nbas
    integer :: i, j
    integer, parameter :: k = 0
    real(r8), dimension(3) :: Rsm
    real(r8), dimension(1) :: Tsm
    real(r8), dimension(:,:), allocatable :: Q0_ints

    nbas = size(fock)

    allocate(Q0_ints(nbas,1))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        if (abs(Q0(i,1)) < zero) cycle

        call get_Qk_integrals(Q0_ints, k, Rs(i,:), Q0(i,:), work)

        ! nuclei - monopole interaction
        do j = 1, qmnucs
            Rsm = Rm(j,:) - Rs(i,:)
            call get_Tk_tensor(Tsm, k, Rsm)
            E_nuc = E_nuc + Q0(i,1) * Zm(j) * Tsm(1)
        end do

        ! electron - monopole interaction
        E_el = E_el + dot(density, Q0_ints(:,1))
        fock = fock + Q0_ints(:,1)

    end do

    deallocate(Q0_ints)

end subroutine es_monopoles

!------------------------------------------------------------------------------

subroutine es_dipoles(density, fock, E_el, E_nuc, work)

    real(r8), dimension(:), intent(in) :: density
    real(r8), dimension(:), intent(inout) :: fock, work
    real(r8), intent(out) :: E_el, E_nuc

    integer :: nbas
    integer :: i, j, l
    integer, parameter :: k = 1
    real(r8), dimension(3) :: Rsm, Tsm
    real(r8), dimension(:,:), allocatable :: Q1_ints

    nbas = size(fock)

    allocate(Q1_ints(nbas,3))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        if (abs(maxval(Q1(i,:))) < zero) cycle

        call get_Qk_integrals(Q1_ints, k, Rs(i,:), Q1(i,:), work)

        ! nuclei - dipole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(j,:) - Rs(i,:)
            call get_Tk_tensor(Tsm, k, Rsm)
            do l = 1, 3
                E_nuc = E_nuc - Zm(j) * Q1(i,l) * Tsm(l)
            end do
        end do

        ! electron - dipole interaction
        do j = 1, 3
            E_el = E_el + dot(density, Q1_ints(:,j))
            fock = fock + Q1_ints(:,j)
        end do

    end do

    deallocate(Q1_ints)

end subroutine es_dipoles

!------------------------------------------------------------------------------

subroutine es_quadrupoles(density, fock, E_el, E_nuc, work)

    real(r8), dimension(:), intent(in) :: density
    real(r8), dimension(:), intent(inout) :: fock, work
    real(r8), intent(out) :: E_el, E_nuc

    integer :: nbas
    integer :: i, j, l
    integer, parameter :: k = 2
    real(r8), dimension(3) :: Rsm
    real(r8), dimension(6) :: Tsm, factors
    real(r8), dimension(:,:), allocatable :: Q2_ints

    nbas = size(fock)

    allocate(Q2_ints(nbas,6))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        if (abs(maxval(Q2(i,:))) < zero) cycle

        call get_Qk_integrals(Q2_ints, k, Rs(i,:), Q2(i,:), work)

        ! nuclei - quadrupole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(j,:) - Rs(i,:)
            call get_Tk_tensor(Tsm, k, Rsm)
            call get_symmetry_factors(factors, k)
            do l = 1, 6
                E_nuc = E_nuc + 0.5d0 * factors(l) * Zm(j) * Q2(i,l) * Tsm(l)
            end do

        end do

        ! electron - quadrupole interaction energy
        do j = 1, 6
            fock = fock + Q2_ints(:,j)
            E_el = E_el + dot(density, Q2_ints(:,j))
        end do

    end do

    deallocate(Q2_ints)

end subroutine es_quadrupoles

!------------------------------------------------------------------------------

subroutine es_octopoles(density, fock, E_el, E_nuc, work)

    real(r8), dimension(:), intent(in) :: density
    real(r8), dimension(:), intent(inout) :: fock, work
    real(r8), intent(out) :: E_el, E_nuc

    integer :: nbas
    integer :: i, j, l
    integer, parameter :: k = 3
    real(r8), dimension(3) :: Rsm
    real(r8), dimension(10) :: Tsm, factors
    real(r8), dimension(:,:), allocatable :: Q3_ints

    nbas = size(fock)

    allocate(Q3_ints(nbas,10))

    E_el = 0.0d0; E_nuc = 0.0d0

    do i = 1, ncents

        if (abs(maxval(Q3(i,:))) < zero) cycle

        call get_Qk_integrals(Q3_ints, k, Rs(i,:), Q3(i,:), work)

        ! nuclei - octopole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(j,:) - Rs(i,:)
            call get_Tk_tensor(Tsm, k, Rsm)
            call get_symmetry_factors(factors, k)
            do l = 1, 10
                E_nuc = E_nuc - factors(l) * Zm(j) * Q3(i,l) * Tsm(l) / 6.0d0
            end do
        end do

        ! electron - octopole interaction energy
        do j = 1, 10
            fock = fock + Q3_ints(:,j)
            E_el = E_el + dot(density, Q3_ints(:,j))
        end do

    end do

    deallocate(Q3_ints)

end subroutine es_octopoles

!------------------------------------------------------------------------------

subroutine pe_polarization(density, fock, E_ind, work)

    real(r8), intent(out) :: E_ind
    real(r8), dimension(:), intent(in) :: density
    real(r8), dimension(:), intent(inout) :: fock, work

    integer :: npol
    integer :: i, j
    real(r8), dimension(:,:), allocatable :: Mu

    E_ind = 0.0d0

    if (poltype >= 0) then

        npol = 0
        do i = 1, ncents
            if (abs(maxval(alpha(i,:))) >= zero) then
                npol = npol + 1
            end if
        end do

        allocate(Mu(npol,3))

        Mu = 0.0d0

        call get_induced_dipoles(Mu, work)

        deallocate(Mu)

    end if

end subroutine pe_polarization

!------------------------------------------------------------------------------

subroutine get_induced_dipoles(Mu, work)

    real(r8), dimension(:,:), intent(out) :: Mu
    real(r8), dimension(:), intent(inout) :: work

    integer :: npol
    real(r8), dimension(:,:), allocatable :: Ftot

    npol = size(Mu, 1)

    allocate(Ftot(npol,3))

    Ftot = 0.0d0
    
    call get_electric_fields(Ftot, work)

    deallocate(Ftot)

end subroutine get_induced_dipoles

!------------------------------------------------------------------------------

subroutine get_electric_fields(Ftot, work)

    real(r8), dimension(:,:), intent(out) :: Ftot
    real(r8), dimension(:), intent(inout) :: work

    integer :: i, j
    integer :: npol
    real(r8), dimension(:,:), allocatable :: Fel, Fnuc, Fmul

    npol = size(Ftot, 1)

    allocate(Fel(npol,3), Fnuc(npol,3), Fmul(npol,3))

    Fel = 0.0d0; Fnuc = 0.0d0; Fmul = 0.0d0

    call get_electron_fields(Fel, work)

    call get_nuclear_fields(Fnuc, work)

    call get_multipole_fields(Fmul, work)

end subroutine get_electric_fields

!------------------------------------------------------------------------------

subroutine get_electron_fields(F, work)

    real(r8), dimension(:,:), intent(out) :: F
    real(r8), dimension(:), intent(inout) :: work

    integer :: i, j

    do i = 1, ncents

        print *, 'temp'

    end do

end subroutine get_electron_fields

!------------------------------------------------------------------------------

subroutine get_nuclear_fields(F, work)

    real(r8), dimension(:,:), intent(out) :: F
    real(r8), dimension(:), intent(inout) :: work

    integer :: i, j

    do i = 1, ncents

        print *, 'temp'

    end do

end subroutine get_nuclear_fields

!------------------------------------------------------------------------------

subroutine get_multipole_fields(Fmul, work)

    real(r8), dimension(:,:), intent(out) :: Fmul
    real(r8), dimension(:), intent(inout) :: work

    logical :: skip
    integer :: i, j, k
    real(r8) :: Rss
    real(r8), dimension(3) :: Fs

    do i = 1, ncents
        do j = 1, ncents

            ! check if j is allowed to polarize i
            if (i == j) cycle
            skip = .false.
            do k = 1, nexlist
                if (exlist(i,1) == exlist(j,k)) then
                    skip = .true.
                    exit
                end if
            end do
            if (skip) cycle

            ! cutoff
!            Rij = Rs(j,:) - Rs(i,:)

            ! get electric field at i due to monopole at j
            if (mulorder >= 0) then
                ! skip if monopole is 'zero'
                if (abs(maxval(Q0(j,:))) >= zero) then
                    call get_monopole_field(Fs, Rs(i,:), Rs(j,:), Q0(j,:))
                    Fmul(i,:) = Fmul(i,:) + Fs
                end if
            end if

            ! get electric field at i due to dipole at j
            if (mulorder >= 1) then
                ! skip if dipole is 'zero'
                if (abs(maxval(Q1(j,:))) >= zero) then
                    call get_dipole_field(Fs, Rs(i,:), Rs(j,:), Q1(j,:))
                    Fmul(i,:) = Fmul(i,:) + Fs
                end if
            end if

            ! get electric field at i due to quadrupole at j
!            if (mulorder >= 2) then
!                ! skip if quadrupole is 'zero'
!                if (abs(maxval(Q2(j,:))) >= zero) then
!                    call get_quadrupole_field(F)
!                    Fmul(i,:) = Fmul(i,:) + F
!                end if
!            end if

            ! get electric field at i due to octopole at j
!            if (mulorder >= 3) then
!                ! skip if octopole is 'zero'
!                if (abs(maxval(Q3(j,:))) >= zero) then
!                    call get_octopole_field(F)
!                    Fmul(i,:) = Fmul(i,:) + F
!                end if
!            end if
        end do
    end do

end subroutine get_multipole_fields

!------------------------------------------------------------------------------

subroutine get_monopole_field(Fa, Ra, Rb, Q0b)

    real(r8), dimension(3), intent(out) :: Fa
    real(r8), dimension(3), intent(in) :: Ra, Rb
    real(r8), dimension(1), intent(in) :: Q0b

    integer :: i
    integer, parameter :: k = 1
    real(r8), dimension(3) :: Rab
    real(r8), dimension(3) :: Tab

    Rab = Rb - Ra

    call get_Tk_tensor(Tab, k, Rab)

    Fa = 0.0d0

    do i = 1, 3
        Fa(i) = - Q0b(1) * Tab(i)
    end do

end subroutine get_monopole_field

!------------------------------------------------------------------------------

subroutine get_dipole_field(Fa, Ra, Rb, Q1b)

    real(r8), dimension(3), intent(out) :: Fa
    real(r8), dimension(3), intent(in) :: Ra, Rb, Q1b

    integer :: i, j
    integer, parameter :: k = 2
    real(r8), dimension(3) :: Rab
    real(r8), dimension(6) :: Tab
    real(r8), dimension(3,3) :: Tf

    Rab = Rb - Ra

    call get_Tk_tensor(Tab, k, Rab)

    call get_full_2nd_tensor(Tf, Tab)

    Fa = 0.0d0

    do i = 1, 3
        do j = 1, 3
            Fa(i) = Fa(i) + Tf(i,j) * Q1b(j)
        end do
    end do

!    Fa = matmul(Tf, Q1b)

!    Fa(1) = Tab(1) * Q1b(1) + Tab(2) * Q1b(2) + Tab(3) * Q1b(3)
!    Fa(2) = Tab(2) * Q1b(1) + Tab(4) * Q1b(2) + Tab(5) * Q1b(3)
!    Fa(3) = Tab(3) * Q1b(1) + Tab(5) * Q1b(2) + Tab(6) * Q1b(3)

end subroutine get_dipole_field

!------------------------------------------------------------------------------

subroutine get_quadrupole_field(Fa, Ra, Rb, Q2b)

    real(r8), dimension(3), intent(out) :: Fa
    real(r8), dimension(3), intent(in) :: Ra, Rb, Q2b

    integer :: i, j, l
    integer, parameter :: k = 3
    real(r8), dimension(3) :: Rab
    real(r8), dimension(10) :: Tab
    real(r8), dimension(3,3) :: Q2f
    real(r8), dimension(3,3,3) :: Tf

    Rab = Rb - Ra

    call get_Tk_tensor(Tab, k, Rab)

    call get_full_2nd_tensor(Q2f, Q2b)
    call get_full_3rd_tensor(Tf, Tab)

    Fa = 0.0d0

    do i = 1, 3
        do j = 1, 3
            do l = 1,3
                Fa(i) = Fa(i) + 0.5d0 * Tf(i,j,l) * Q2f(j,l)
            end do
        end do
    end do

end subroutine get_quadrupole_field

!------------------------------------------------------------------------------

subroutine get_octopole_field(Fa, Ra, Rb, Q3b)

    real(r8), dimension(3), intent(out) :: Fa
    real(r8), dimension(3), intent(in) :: Ra, Rb, Q3b

    integer :: i, j, l, m
    integer, parameter :: k = 4
    real(r8), dimension(3) :: Rab
    real(r8), dimension(15) :: Tab
    real(r8), dimension(3,3,3) :: Q3f
    real(r8), dimension(3,3,3,3) :: Tf

    Rab = Rb - Ra

    call get_Tk_tensor(Tab, k, Rab)

    call get_full_3rd_tensor(Q3f, Q3b)
    call get_full_4th_tensor(Tf, Tab)

    Fa = 0.0d0

    do i = 1, 3
        do j = 1, 3
            do l = 1,3
                do m = 1, 3
                    Fa(i) = Fa(i) + Tf(i,j,l,m) * Q3f(j,l,m) / 6.0d0
                end do
            end do
        end do
    end do

end subroutine get_octopole_field

!------------------------------------------------------------------------------

!ubroutine create_response_matrix(A, work)

! TODO: Damping schemes
!       Cutoff radius
!       Check if A is allocated

!   real(r8), dimension(:), intent(out) :: A
!   real(r8), dimension(:), intent(inout) :: work

!   integer :: i, j, k, l, m, n
!   integer :: info, stat_alloc
!   real(r8), dimension(6) :: alpha
!   real(r8) :: Rij, T2
!   logical :: exclude

!   

! A = 0.0d0
! alpha = 0.0d0
! T2 = 0.0d0
! exclude = .false.

! m = 0

! do i = 1, ncents

!   if (zeropol(i) == -1) cycle

!   ! Inverse of diagonal isotropic polarizability tensor
!   if (poltype == 0) then

!     alpha(1) = 1.0d0/isoalpha(i)
!     alpha(4) = 1.0d0/isoalpha(i)
!     alpha(6) = 1.0d0/isoalpha(i)

!   ! Inverse of symmetric and positive-definite polarizability tensor
!   else if (poltype == 1) then

!     alpha = dipols(i,:)

!     ! First we factorize
!     call dpptrf('L', 3, alpha, info)
!     if (info > 0) then
!       print *, 'WARNING'
!       print *, 'Non-positive-definite polarizability found.'
!       print *, 'Offending polarizabilty is located at site no.: ', i
!       print *, 'Polarizability (xx, xy, xz, yy, yz, zz) :'
!       do k = 1, 6
!         print *, '       ', dipdippols(i,k)
!       end do
!       print *, 'Skipping the polarizability at this site'
!       cycle
!     else if (info < 0) then
!       print *, 'ERROR'
!       print *, 'Illegal value in polarizability at site no.: ', i
!       print *, 'Polarizability (xx, xy, xz, yy, yz, zz) :'
!       do k = 1, 6
!         print *, '       ', dipdippols(i,k)
!       end do
!       stop 'Quitting...'
!     end if

!     ! and then we invert
!     call dpptri('L', 3, alpha, info)
!     if (info > 0 .or. info < 0) then
!       print *, 'ERROR'
!       print *, 'Inversion of polarizability tensor failed.'
!       print *, 'Offending polarizabilty tensor is located at site no.: ', i
!       print *, 'Polarizability (xx, xy, xz, yy, yz, zz) :'
!       do k = 1, 6
!         print *, '       ', dipdippols(i,k)
!       end do
!       stop 'Quitting...'
!     end if

!   end if

!   do l = 3, 1, -1

!     do j = i, nsites

!       if (zeropol(j) == -1) cycle

!       if (j == i) then

!         if (l == 3) then
!           do k = 1, l
!             A(m+k) = alpha(k)
!           end do
!         else if (l == 2) then
!           do k = 1, l
!             A(m+k) = alpha(3+k)
!           end do
!         else if (l == 1) then
!             A(m+1) = alpha(5+1)
!         end if
!         m = m + l

!       else

!         R = 0.0d0
!         do k = 1, 3
!           R = R + (Rs(i,k) - Rs(j,k))**2
!         end do
!         R = sqrt(R)

!         if (R > cutoff) then
!           m = m + 3
!           cycle
!         end if

!         exclude = .false.
!         do k = 1, nexclist
!           if (exclist(i,k) == exclist(j,1)) exclude = .true.
!         end do

!         if (exclude) then
!           m = m + 3
!           cycle
!         end if

!         if (l == 3) then
!           do k = 1, 3
!             T = (3.0d0*(Rs(i,1) - Rs(j,1))*(Rs(i,k) -
!oords(j,k)))/R**5
!             if (k == 1) T = T - 1.0d0/R**3
!             A(m+k) = - T
!           end do
!         else if (l == 2) then
!           do k = 1, 3
!             T = (3.0d0*(Rs(i,2) - Rs(j,2))*(Rs(i,k) -
!oords(j,k)))/R**5
!             if (k == 2) T = T - 1.0d0/R**3
!             A(m+k) = - T
!           end do
!         else if (l == 1) then
!           do k = 1, 3
!             T = (3.0d0*(Rs(i,3) - Rs(j,3))*(Rs(i,k) -
!oords(j,k)))/R**5
!             if (k == 3) T = T - 1.0d0/R**3
!             A(m+k) = - T
!           end do
!         end if
!         m = m + 3
!       end if
!     end do
!   end do
! end do

!nd subroutine create_response_matrix

!------------------------------------------------------------------------------

subroutine get_Tk_tensor(Tk, k, Rij)

    integer, intent(in) :: k
    real(r8), dimension(:), intent(out) :: Tk
    real(r8), dimension(3), intent(in) :: Rij

    real(r8) :: R, R2, R4

    R = nrm2(Rij)

    if (k == 0) then
        Tk(1) = 1.0d0 / R
    else if (k == 1) then
        Tk(1) = Rij(1)
        Tk(2) = Rij(2)
        Tk(3) = Rij(3)
        Tk = - Tk / R**3
    else if (k == 2) then
        R2 = R**2
        Tk(1) = 3.0d0 * Rij(1)**2 - R2
        Tk(2) = 3.0d0 * Rij(1) * Rij(2)
        Tk(3) = 3.0d0 * Rij(1) * Rij(3)
        Tk(4) = 3.0d0 * Rij(2)**2 - R2
        Tk(5) = 3.0d0 * Rij(2) * Rij(3)
        Tk(6) = 3.0d0 * Rij(3)**2 - R2
        Tk = Tk / R**5
    else if (k == 3) then
        R2 = R**2
        Tk(1) = 15.0d0 * Rij(1)**3 - 9.0d0 * R2 * Rij(1)
        Tk(2) = 15.0d0 * Rij(1)**2 * Rij(2) - 3.0d0 * R2 * Rij(2)
        Tk(3) = 15.0d0 * Rij(1)**2 * Rij(3) - 3.0d0 * R2 * Rij(3)
        Tk(4) = 15.0d0 * Rij(1) * Rij(2)**2 - 3.0d0 * R2 * Rij(1)
        Tk(5) = 15.0d0 * Rij(1) * Rij(2) * Rij(3)
        Tk(6) = 15.0d0 * Rij(1) * Rij(3)**2 - 3.0d0 * R2 * Rij(1)
        Tk(7) = 15.0d0 * Rij(2)**3 - 9.0d0 * R2 * Rij(2)
        Tk(8) = 15.0d0 * Rij(2)**2 * Rij(3) - 3.0d0 * R2 * Rij(3)
        Tk(9) = 15.0d0 * Rij(2) * Rij(3)**2 - 3.0d0 * R2 * Rij(2)
        Tk(10) = 15.0d0 * Rij(3)**3 - 9.0d0 * R2 * Rij(3)
        Tk = - Tk / R**7
    else if (k == 4) then
        R2 = R**2
        R4 = R**4
        Tk(1) = 105.0d0 * Rij(1)**4 - 90.0d0 * R2 * Rij(1)**2 + 9.0d0 * R4
        Tk(2) = 105.0d0 * Rij(1)**3 * Rij(2) - 45.0d0 * R2 * Rij(1) * Rij(2)
        Tk(3) = 105.0d0 * Rij(1)**3 * Rij(3) - 45.0d0 * R2 * Rij(1) * Rij(3)
        Tk(4) = 105.0d0 * Rij(1)**2 * Rij(2)**2 &
                - 30.0d0 * R2 * (Rij(1)**2 + Rij(2)**2) + 3.0d0 * R4
        Tk(5) = 105.0d0 * Rij(1)**2 * Rij(2) * Rij(3) &
                - 15.0d0 * R2 * Rij(2) * Rij(3)
        Tk(6) = 105.0d0 * Rij(1)**2 * Rij(3)**2 &
                - 30.0d0 * R2 * (Rij(1)**2 + Rij(3)**2) + 3.0d0 * R4
        Tk(7) = 105.0d0 * Rij(1) * Rij(2)**3 - 45.0d0 * R2 * Rij(1) * Rij(2)
        Tk(8) = 105.0d0 * Rij(1) * Rij(2)**2 * Rij(3) &
                - 15.0d0 * R2 * Rij(1) * Rij(3)
        Tk(9) = 105.0d0 * Rij(1) * Rij(2) * Rij(3)**2 &
                - 15.0d0 * R2 * Rij(1) * Rij(2)
        Tk(10) = 105.0d0 * Rij(1) * Rij(3)**3 - 45.0d0 * R2 * Rij(1) * Rij(3)
        Tk(11) = 105.0d0 * Rij(2)**4 - 90.0d0 * R2 * Rij(2)**2 + 9.0d0 * R4
        Tk(12) = 105.0d0 * Rij(2)**3 * Rij(3) - 45.0d0 * R2 * Rij(2) * Rij(3)
        Tk(13) = 105.0d0 * Rij(2)**2 * Rij(3)**2 &
                - 30.0d0 * R2 * (Rij(2)**2 + Rij(3)**2) + 3.0d0 * R4
        Tk(14) = 105.0d0 * Rij(2) * Rij(3)**3 - 45.0d0 * R2 * Rij(2) * Rij(3)
        Tk(15) = 105.0d0 * Rij(3)**4 - 90.0d0 * R2 * Rij(3)**2 + 9.0d0 * R4
        Tk = Tk / R**9
    end if

end subroutine get_Tk_tensor

!------------------------------------------------------------------------------

subroutine get_full_2nd_tensor(Tf, Ts)

    real(r8), dimension(:), intent(in) :: Ts
    real(r8), dimension(3,3), intent(out) :: Tf

    Tf = 0.0d0

    Tf(1,1) = Ts(1); Tf(1,2) = Ts(2); Tf(1,3) = Ts(3)
    Tf(2,1) = Ts(2); Tf(2,2) = Ts(4); Tf(2,3) = Ts(5)
    Tf(3,1) = Ts(3); Tf(3,2) = Ts(5); Tf(3,3) = Ts(6)

end subroutine get_full_2nd_tensor

!------------------------------------------------------------------------------

subroutine get_full_3rd_tensor(Tf, Ts)

    real(r8), dimension(:), intent(in) :: Ts
    real(r8), dimension(3,3,3), intent(out) :: Tf

    Tf = 0.0d0

    Tf(1,1,1) = Ts(1); Tf(1,1,2) = Ts(2); Tf(1,1,3) = Ts(3)
    Tf(1,2,1) = Ts(2); Tf(1,2,2) = Ts(4); Tf(1,2,3) = Ts(5)
    Tf(1,3,1) = Ts(3); Tf(1,3,2) = Ts(5); Tf(1,3,3) = Ts(6)
    Tf(2,1,1) = Ts(2); Tf(2,1,2) = Ts(4); Tf(2,1,3) = Ts(5)
    Tf(2,2,1) = Ts(4); Tf(2,2,2) = Ts(7); Tf(2,2,3) = Ts(8)
    Tf(2,3,1) = Ts(5); Tf(2,3,2) = Ts(8); Tf(2,3,3) = Ts(9)
    Tf(3,1,1) = Ts(3); Tf(3,1,2) = Ts(5); Tf(3,1,3) = Ts(6)
    Tf(3,2,1) = Ts(5); Tf(3,2,2) = Ts(8); Tf(3,2,3) = Ts(9)
    Tf(3,3,1) = Ts(6); Tf(3,3,2) = Ts(9); Tf(3,3,3) = Ts(10)

end subroutine get_full_3rd_tensor

!------------------------------------------------------------------------------

subroutine get_full_4th_tensor(Tf, Ts)

    real(r8), dimension(:), intent(in) :: Ts
    real(r8), dimension(3,3,3,3), intent(out) :: Tf

    Tf = 0.0d0

    Tf(1,1,1,1) = Ts(1); Tf(1,1,1,2) = Ts(2); Tf(1,1,1,3) = Ts(3)
    Tf(1,1,2,1) = Ts(2); Tf(1,1,2,2) = Ts(4); Tf(1,1,2,3) = Ts(5)
    Tf(1,1,3,1) = Ts(3); Tf(1,1,3,2) = Ts(5); Tf(1,1,3,3) = Ts(6)
    Tf(1,2,1,1) = Ts(2); Tf(1,2,1,2) = Ts(4); Tf(1,2,1,3) = Ts(5)
    Tf(1,2,2,1) = Ts(4); Tf(1,2,2,2) = Ts(7); Tf(1,2,2,3) = Ts(8)
    Tf(1,2,3,1) = Ts(5); Tf(1,2,3,2) = Ts(8); Tf(1,2,3,3) = Ts(9)
    Tf(1,3,1,1) = Ts(3); Tf(1,3,1,2) = Ts(5); Tf(1,3,1,3) = Ts(6)
    Tf(1,3,2,1) = Ts(5); Tf(1,3,2,2) = Ts(8); Tf(1,3,2,3) = Ts(9)
    Tf(1,3,3,1) = Ts(6); Tf(1,3,3,2) = Ts(9); Tf(1,3,3,3) = Ts(10)
    Tf(2,1,1,1) = Ts(2); Tf(2,1,1,2) = Ts(4); Tf(2,1,1,3) = Ts(5)
    Tf(2,1,2,1) = Ts(4); Tf(2,1,2,2) = Ts(7); Tf(2,1,2,3) = Ts(8)
    Tf(2,1,3,1) = Ts(5); Tf(2,1,3,2) = Ts(8); Tf(2,1,3,3) = Ts(9)
    Tf(2,2,1,1) = Ts(4); Tf(2,2,1,2) = Ts(7); Tf(2,2,1,3) = Ts(8)
    Tf(2,2,2,1) = Ts(7); Tf(2,2,2,2) = Ts(11); Tf(2,2,2,3) = Ts(12)
    Tf(2,2,3,1) = Ts(8); Tf(2,2,3,2) = Ts(12); Tf(2,2,3,3) = Ts(13)
    Tf(2,3,1,1) = Ts(5); Tf(2,3,1,2) = Ts(8); Tf(2,3,1,3) = Ts(9)
    Tf(2,3,2,1) = Ts(8); Tf(2,3,2,2) = Ts(12); Tf(2,3,2,3) = Ts(13)
    Tf(2,3,3,1) = Ts(9); Tf(2,3,3,2) = Ts(13); Tf(2,3,3,3) = Ts(14)
    Tf(3,1,1,1) = Ts(3); Tf(3,1,1,2) = Ts(5); Tf(3,1,1,3) = Ts(6)
    Tf(3,1,2,1) = Ts(5); Tf(3,1,2,2) = Ts(8); Tf(3,1,2,3) = Ts(9)
    Tf(3,1,3,1) = Ts(6); Tf(3,1,3,2) = Ts(9); Tf(3,1,3,3) = Ts(10)
    Tf(3,2,1,1) = Ts(5); Tf(3,2,1,2) = Ts(8); Tf(3,2,1,3) = Ts(9)
    Tf(3,2,2,1) = Ts(8); Tf(3,2,2,2) = Ts(12); Tf(3,2,2,3) = Ts(13)
    Tf(3,2,3,1) = Ts(9); Tf(3,2,3,2) = Ts(13); Tf(3,2,3,3) = Ts(14)
    Tf(3,3,1,1) = Ts(6); Tf(3,3,1,2) = Ts(9); Tf(3,3,1,3) = Ts(10)
    Tf(3,3,2,1) = Ts(9); Tf(3,3,2,2) = Ts(13); Tf(3,3,2,3) = Ts(14)
    Tf(3,3,3,1) = Ts(10); Tf(3,3,3,2) = Ts(14); Tf(3,3,3,3) = Ts(15)

end subroutine get_full_4th_tensor

!------------------------------------------------------------------------------

subroutine get_Qk_integrals(Qk_ints, k, coord, Qk, work)

    external get_Tk_integrals

    integer, intent(in) :: k
    real(r8), dimension(:,:), intent(out) :: Qk_ints
    real(r8), dimension(:), intent(in) :: Qk
    real(r8), dimension(3), intent(in) :: coord
    real(r8), dimension(:), intent(inout) :: work

    integer :: i
    integer :: ncomps, nbas
    real(r8) :: taylor
    real(r8), dimension(:), allocatable :: factors

    if (k == 0) then
        taylor = 1.0d0
    else if (k == 1) then
        taylor = - 1.0d0
    else if (k == 2) then
        taylor = 0.5d0
    else if (k == 3) then
        taylor = - 1.0d0 / 6.0d0
    else if (k == 4) then
        taylor = 1.0d0 / 24.0d0
    end if

    nbas = size(Qk_ints, 1)
    ncomps = size(Qk_ints, 2)

    ! get T^(k) integrals (incl. negative sign from electron density)
    call get_Tk_integrals(Qk_ints, ncomps*nbas, k, 'packed', coord,&
                          work, size(work))

    ! get symmetry factors
    allocate(factors(ncomps))
    call get_symmetry_factors(factors, k)

    ! dot T^(k) integrals with multipole to get Q^(k) integrals
    do i = 1, ncomps
        Qk_ints(:,i) = taylor * factors(i) * Qk(i) * Qk_ints(:,i)
    end do

    deallocate(factors)

end subroutine get_Qk_integrals

!------------------------------------------------------------------------------

subroutine get_symmetry_factors(factors, k)

    integer, intent(in) :: k
    real(r8), dimension(:), intent(out) :: factors

    factors = 0.0d0

    if (k == 0) then
        factors(1) = 1.0d0
    else if (k == 1) then
        factors(1) = 1.0d0
        factors(2) = 1.0d0
        factors(3) = 1.0d0
    else if (k== 2) then
        factors(1) = 1.0d0
        factors(2) = 2.0d0
        factors(3) = 2.0d0
        factors(4) = 1.0d0
        factors(5) = 2.0d0
        factors(6) = 1.0d0
    else if (k == 3) then
        factors(1) = 1.0d0
        factors(2) = 3.0d0
        factors(3) = 3.0d0
        factors(4) = 3.0d0
        factors(5) = 6.0d0
        factors(6) = 3.0d0
        factors(7) = 1.0d0
        factors(8) = 3.0d0
        factors(9) = 3.0d0
        factors(10) = 1.0d0
    else if (k == 4) then
        factors(1) = 1.0d0
        factors(2) = 4.0d0
        factors(3) = 4.0d0
        factors(4) = 6.0d0
        factors(5) = 12.0d0
        factors(6) = 6.0d0
        factors(7) = 4.0d0
        factors(8) = 12.0d0
        factors(9) = 12.0d0
        factors(10) = 4.0d0
        factors(11) = 1.0d0
        factors(12) = 4.0d0
        factors(13) = 6.0d0
        factors(14) = 4.0d0
        factors(15) = 1.0d0
    end if

end subroutine get_symmetry_factors

!------------------------------------------------------------------------------

function nrm2(x)

    real(r8), external :: dnrm2
    real(r8) :: nrm2
    real(r8), dimension(:), intent(in) :: x

    nrm2 = dnrm2(size(x), x, 1)

end function nrm2

!------------------------------------------------------------------------------

function dot(x,y)

    real(r8), external :: ddot
    real(r8) :: dot
    real(r8), dimension(:), intent(in) :: x, y

    dot = ddot(size(x), x, 1, y, 1)

end function dot

!------------------------------------------------------------------------------

!subroutine symm()

!    real(r8), external :: dsymm
!    real(r8) :: symm
!    real(r8), dimension(:) :: 

!end subroutine symm

!------------------------------------------------------------------------------

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

!------------------------------------------------------------------------------

end module polarizable_embedding
