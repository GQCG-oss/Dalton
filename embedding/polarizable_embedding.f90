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
    real(r8), dimension(:,:), allocatable, save :: alphas
    ! dipole-dipole-dipole polarizability / 1st hyperpolarizability
    real(r8), dimension(:,:), allocatable, save :: betas
    ! dipole-quadrupole polarizability
    real(r8), dimension(:,:), allocatable, save :: As
    ! quadrupole-quadrupole polarizability
    real(r8), dimension(:,:), allocatable, save :: Cs

    ! QM info: number of nuclei, nuclear charges and nuclear coordinates
    integer, save :: qmnucs
    real(r8), dimension(:), allocatable, save :: Zm
    real(r8), dimension(:,:), allocatable, save :: Rm

! TODO:
! AA and AU
! initialize when allocate
! transpose to follow array element ordering
! save individual one-electron integrals and reuse
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

    integer, intent(in) :: natoms, nwrk
    real(r8), dimension(natoms), intent(in) :: charge
    real(r8), dimension(3,natoms), intent(in) :: cord
    real(r8), dimension(nwrk), intent(inout) :: work

    integer :: i, j, s
    integer :: lupot, nlines

    qmnucs = natoms

    allocate(Rm(3,qmnucs), Zm(qmnucs))

    Rm = cord
    Zm = charge

    call openfile('POTENTIAL.INP', lupot, 'old', 'formatted')

    read(lupot,*) mulorder, poltype, hypoltype
    read(lupot,*) ncents

    allocate(elems(ncents), Rs(3,ncents))
    elems = 0.0d0; Rs = 0.0d0

    do i = 1, ncents
        read(lupot,*) elems(i), (Rs(j,i), j = 1, 3)
    end do

    if (mulorder >= 0) then
        read(lupot,*) nlines
        allocate(Q0(1,ncents))
        Q0 = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) Q0(1,s)
        end do
    end if

    if (mulorder >= 1) then
        read(lupot,*) nlines
        allocate(Q1(3,ncents))
        Q1 = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (Q1(j,s), j = 1, 3)
        end do
    end if

    if (mulorder >= 2) then
        read(lupot,*) nlines
        allocate(Q2(6,ncents))
        Q2 = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (Q2(j,s), j = 1, 6)
        end do
    end if

    if (mulorder >= 3) then
        read(lupot,*) nlines
        allocate(Q3(10,ncents))
        Q3 = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (Q3(j,s), j = 1, 10)
        end do
    end if

    if (mulorder >= 4) then
        read(lupot,*) nlines
        allocate(Q4(15,ncents))
        Q4 = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (Q4(j,s), j = 1, 15)
        end do
    end if

!    if (poltype >= 0) then
!        allocate(alphas(6,ncents))
!        alphas = 0.0d0
!        read(lupot,*) nlines
!        do i = 1, nlines
!            read(lupot,*) s
!            read(lupot,*) alphas(1,s)
!            alphas(4,s) = alphas(1,s)
!            alphas(6,s) = alphas(1,s)
!        end do
!    end if

    if (poltype >= 1) then
        if (.not. allocated(alphas)) then
            allocate(alphas(6,ncents))
            alphas = 0.0d0
        end if
        read(lupot,*) nlines
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (alphas(j,s), j = 1, 6)
        end do
    end if

    if (poltype >= 2) then
        allocate(As(10,ncents))
        As = 0.0d0
        read(lupot,*) nlines
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (As(j,s), j = 1, 10)
        end do
    end if

    if (poltype >= 3) then
        allocate(Cs(15,ncents))
        Cs = 0.0d0
        read(lupot,*) nlines
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (Cs(j,s), j = 1, 15)
        end do
    end if

    if (hypoltype >= 1) then
        allocate(betas(10,ncents))
        betas = 0.0d0
        read(lupot,*) nlines
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (betas(j,s), j = 1, 10)
        end do
    end if

    if (poltype >= 1) then
        read(lupot,*) nlines, nexlist
        allocate(exlist(nexlist,ncents))
        exlist = 0.0d0
        do i = 1, nlines
            read(lupot,*) s
            read(lupot,*) (exlist(j,s), j = 1, nexlist)
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
        if (mulorder >= 0) then
            call es_monopoles(density, fock, E_el(0), E_nuc(0), work)
        end if

        if (mulorder >= 1) then
            call es_dipoles(density, fock, E_el(1), E_nuc(1), work)
        end if

        if (mulorder >= 2) then
            call es_quadrupoles(density, fock, E_el(2), E_nuc(2), work)
        end if

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

        if (abs(Q0(1,i)) < zero) cycle

        call get_Qk_integrals(Q0_ints, k, Rs(:,i), Q0(:,i), work)

        ! nuclei - monopole interaction
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call get_Tk_tensor(Tsm, k, Rsm)
            E_nuc = E_nuc + Q0(1,i) * Zm(j) * Tsm(1)
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

        if (abs(maxval(Q1(:,i))) < zero) cycle

        call get_Qk_integrals(Q1_ints, k, Rs(:,i), Q1(:,i), work)

        ! nuclei - dipole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call get_Tk_tensor(Tsm, k, Rsm)
            do l = 1, 3
                E_nuc = E_nuc - Zm(j) * Q1(l,i) * Tsm(l)
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

        if (abs(maxval(Q2(:,i))) < zero) cycle

        call get_Qk_integrals(Q2_ints, k, Rs(:,i), Q2(:,i), work)

        ! nuclei - quadrupole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call get_Tk_tensor(Tsm, k, Rsm)
            call get_symmetry_factors(factors, k)
            do l = 1, 6
                E_nuc = E_nuc + 0.5d0 * factors(l) * Zm(j) * Q2(l,i) * Tsm(l)
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

        if (abs(maxval(Q3(:,i))) < zero) cycle

        call get_Qk_integrals(Q3_ints, k, Rs(:,i), Q3(:,i), work)

        ! nuclei - octopole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call get_Tk_tensor(Tsm, k, Rsm)
            call get_symmetry_factors(factors, k)
            do l = 1, 10
                E_nuc = E_nuc - factors(l) * Zm(j) * Q3(l,i) * Tsm(l) / 6.0d0
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

    integer :: npol, nbas
    integer :: i, j, l
    integer, parameter :: k = 1
    real(r8) :: E_el, E_nuc, E_mul
    real(r8), dimension(0:3) :: E_Qk
    real(r8), dimension(:), allocatable :: Mu, Fel, Fnuc, Fmul
    real(r8), dimension(:,:), allocatable :: Mu_ints

    E_ind = 0.0d0; E_el = 0.0d0; E_nuc = 0.0d0; E_mul = 0.0d0; E_Qk = 0.0d0

    if (poltype >= 1) then

        npol = 0
        do i = 1, ncents
            if (abs(maxval(alphas(:,i))) >= zero) then
                npol = npol + 1
            end if
        end do

        allocate(Mu_ints(size(fock),3))

        Mu_ints = 0.0d0

        allocate(Mu(3*npol), Fel(3*npol), Fnuc(3*npol), Fmul(3*npol))

        Mu = 0.0d0; Fel = 0.0d0; Fnuc = 0.0d0; Fmul = 0.0d0

        call get_electron_fields(Fel, density, work)

        call get_nuclear_fields(Fnuc, work)

! TODO: check octopoles
        call get_multipole_fields(Fmul, work)

!        call get_electric_fields(Ftot, density, work)

        call get_induced_dipoles(Mu, Fel + Fnuc + Fmul, density, work)

        l = 1

        do i = 1, ncents

            if (abs(maxval(alphas(:,i))) <= zero) cycle

            call get_Qk_integrals(Mu_ints, k, Rs(:,i), Mu(l:l+2), work)

            do j = 1, 3
                fock = fock + Mu_ints(:,j)
            end do

            l = l + 3

        end do

!        E_el = E_el - 0.5d0 * dot(Mu, Fel)
!        E_nuc = E_nuc - 0.5d0 * dot(Mu, Fnuc)
!        E_mul = E_mul - 0.5d0 * dot(Mu, Fmul)

!        E_ind = E_el + E_nuc + E_mul
        E_ind = - 0.5d0 * dot(Mu, Fel + Fnuc + Fmul)

        deallocate(Mu_ints)
        deallocate(Mu, Fel, Fnuc, Fmul)

    end if

end subroutine pe_polarization

!------------------------------------------------------------------------------

subroutine get_induced_dipoles(Mu, F, density, work)

    real(r8), dimension(:), intent(out) :: Mu
    real(r8), dimension(:), intent(in) :: F
    real(r8), dimension(:), intent(in) :: density
    real(r8), dimension(:), intent(inout) :: work

    integer :: i, j, n
    real(r8), dimension(:), allocatable :: A

    n = size(Mu)

    allocate(A(int(n*(n+1)*0.5d0)))

    A = 0.0d0
    
    call get_response_matrix(A, .false., work)

    Mu = F

    call solve(A, Mu)

    deallocate(A)

end subroutine get_induced_dipoles

!------------------------------------------------------------------------------

subroutine get_electric_fields(Ftot, density, work)

    real(r8), dimension(:), intent(out) :: Ftot
    real(r8), dimension(:), intent(in) :: density
    real(r8), dimension(:), intent(inout) :: work

    integer :: i, j
    integer :: n
    real(r8), dimension(:), allocatable :: Fel, Fnuc, Fmul

    n = size(Ftot)

    allocate(Fel(n), Fnuc(n), Fmul(n))

    Fel = 0.0d0; Fnuc = 0.0d0; Fmul = 0.0d0

    call get_electron_fields(Fel, density, work)

    call get_nuclear_fields(Fnuc, work)

    call get_multipole_fields(Fmul, work)

    Ftot = Fel + Fnuc + Fmul

end subroutine get_electric_fields

!------------------------------------------------------------------------------

subroutine get_electron_fields(Fel, density, work)

    external get_Tk_integrals

    real(r8), dimension(:), intent(out) :: Fel
    real(r8), dimension(:), intent(in) :: density
    real(r8), dimension(:), intent(inout) :: work

    integer :: nbas
    integer :: i, j, l
    integer, parameter :: k = 1
    real(r8), dimension(:,:), allocatable :: Fel_ints

    nbas = size(density)

    allocate(Fel_ints(nbas,3))

    l = 0

    do i = 1, ncents

        if (abs(maxval(alphas(:,i))) <= zero) cycle

        call get_Tk_integrals(Fel_ints, 3*nbas, k, 'packed', Rs(:,i),&
                              work, size(work))

        do j = 1, 3
            Fel(l+j) = dot(density, Fel_ints(:,j))
        end do

        l = l + 3

    end do

    deallocate(Fel_ints)

end subroutine get_electron_fields

!------------------------------------------------------------------------------

subroutine get_nuclear_fields(Fnuc, work)

    real(r8), dimension(:), intent(out) :: Fnuc
    real(r8), dimension(:), intent(inout) :: work

    logical :: exclude, lexist
    integer :: lutemp
    integer :: i, j, l, m
    integer, parameter :: k = 1
    real(r8), dimension(3) :: Rms, Tms

    inquire(file='pe_nuclear_field.bin', exist=lexist)

    if (lexist) then
        call openfile('pe_nuclear_field.bin', lutemp, 'old', 'unformatted')
        rewind(lutemp)
        read(lutemp) Fnuc
        close(lutemp)
    else

        l = 0

        do i = 1, ncents

            if (abs(maxval(alphas(:,i))) <= zero) cycle

            do j = 1, qmnucs

                Rms = Rs(:,i) - Rm(:,j)

                call get_Tk_tensor(Tms, k, Rms)

                do m = 1, 3
                    Fnuc(l+m) = Fnuc(l+m) - Zm(j) * Tms(m)
                end do

            end do

            l = l + 3

        end do

        call openfile('pe_nuclear_field.bin', lutemp, 'new', 'unformatted')
        rewind(lutemp)
        write(lutemp) Fnuc
        close(lutemp)

    end if

end subroutine get_nuclear_fields

!------------------------------------------------------------------------------

subroutine get_multipole_fields(Fmul, work)

    real(r8), dimension(:), intent(out) :: Fmul
    real(r8), dimension(:), intent(inout) :: work

    logical :: exclude, lexist
    integer :: lutemp
    integer :: i, j, k, l
    real(r8) :: Rij
    real(r8), dimension(3) :: Fs

    inquire(file='pe_multipole_field.bin', exist=lexist)

    if (lexist) then
        call openfile('pe_multipole_field.bin', lutemp, 'old', 'unformatted')
        rewind(lutemp)
        read(lutemp) Fmul
        close(lutemp)
    else

        l = 1

        do i = 1, ncents

            if (abs(maxval(alphas(:,i))) <= zero) cycle

            do j = 1, ncents

                ! check if j is allowed to polarize i
                if (i == j) cycle

                exclude = .false.

                do k = 1, nexlist
                    if (exlist(k,i) == exlist(1,j)) then
                        exclude = .true.
                        exit
                    end if
                end do

                if (exclude) cycle

                ! cutoff
!                Rij = Rs(:,j) - Rs(:,j)

                ! get electric field at i due to monopole at j
                if (mulorder >= 0) then
                    ! skip if monopole is 'zero'
                    if (abs(maxval(Q0(:,j))) >= zero) then
                        call get_monopole_field(Fs, Rs(:,i), Rs(:,j), Q0(:,j))
                        Fmul(l:l+2) = Fmul(l:l+2) + Fs
                    end if
                end if

                ! get electric field at i due to dipole at j
                if (mulorder >= 1) then
                    ! skip if dipole is 'zero'
                    if (abs(maxval(Q1(:,j))) >= zero) then
                        call get_dipole_field(Fs, Rs(:,i), Rs(:,j), Q1(:,j))
                        Fmul(l:l+2) = Fmul(l:l+2) + Fs
                    end if
                end if

                ! get electric field at i due to quadrupole at j
                if (mulorder >= 2) then
                    ! skip if quadrupole is 'zero'
                    if (abs(maxval(Q2(:,j))) >= zero) then
                        call get_quadrupole_field(Fs, Rs(:,i), Rs(:,j), Q2(:,j))
                        Fmul(l:l+2) = Fmul(l:l+2) + Fs
                    end if
                end if

                ! get electric field at i due to octopole at j
                if (mulorder >= 3) then
                    ! skip if octopole is 'zero'
                    if (abs(maxval(Q3(:,j))) >= zero) then
                        call get_octopole_field(Fs, Rs(:,i), Rs(:,j), Q3(:,j))
                        Fmul(l:l+2) = Fmul(l:l+2) + Fs
                    end if
                end if

            end do

            l = l + 3

        end do

        call openfile('pe_multipole_field.bin', lutemp, 'new', 'unformatted')
        rewind(lutemp)
        write(lutemp) Fmul
        close(lutemp)

    end if

end subroutine get_multipole_fields

!------------------------------------------------------------------------------

subroutine get_monopole_field(Fi, Ri, Rj, Q0j)

    real(r8), dimension(3), intent(out) :: Fi
    real(r8), dimension(3), intent(in) :: Ri, Rj
    real(r8), dimension(1), intent(in) :: Q0j

    integer :: a
    integer, parameter :: k = 1
    real(r8), dimension(3) :: Rij, Tij

    Rij = Rj - Ri

    call get_Tk_tensor(Tij, k, Rij)

    Fi = 0.0d0

    do a = 1, 3
        Fi(a) = Fi(a) + Tij(a) *  Q0j(1)
    end do

end subroutine get_monopole_field

!------------------------------------------------------------------------------

subroutine get_dipole_field(Fi, Ri, Rj, Q1j)

    real(r8), dimension(3), intent(out) :: Fi
    real(r8), dimension(3), intent(in) :: Ri, Rj, Q1j

    integer :: a, b
    integer, parameter :: k = 2
    real(r8), dimension(3) :: Rij
    real(r8), dimension(6) :: Tij
    real(r8), dimension(3,3) :: Tf

    Rij = Rj - Ri

    call get_Tk_tensor(Tij, k, Rij)

    call get_full_2nd_tensor(Tf, Tij)

    Fi = 0.0d0

    do a = 1, 3
        do b = 1, 3
            Fi(a) = Fi(a) + Tf(a,b) * Q1j(b)
        end do
    end do

!    Fa = matmul(Tf, Q1b)

!    Fa(1) = Tab(1) * Q1b(1) + Tab(2) * Q1b(2) + Tab(3) * Q1b(3)
!    Fa(2) = Tab(2) * Q1b(1) + Tab(4) * Q1b(2) + Tab(5) * Q1b(3)
!    Fa(3) = Tab(3) * Q1b(1) + Tab(5) * Q1b(2) + Tab(6) * Q1b(3)

end subroutine get_dipole_field

!------------------------------------------------------------------------------

subroutine get_quadrupole_field(Fi, Ri, Rj, Q2j)

    real(r8), dimension(3), intent(out) :: Fi
    real(r8), dimension(3), intent(in) :: Ri, Rj, Q2j

    integer :: a, b, g
    integer, parameter :: k = 3
    real(r8), dimension(3) :: Rij
    real(r8), dimension(10) :: Tij
    real(r8), dimension(3,3) :: Q2f
    real(r8), dimension(3,3,3) :: Tf

    Rij = Rj - Ri

    call get_Tk_tensor(Tij, k, Rij)

    call get_full_2nd_tensor(Q2f, Q2j)
    call get_full_3rd_tensor(Tf, Tij)

    Fi = 0.0d0

    do a = 1, 3
        do b = 1, 3
            do g = 1, 3
                Fi(a) = Fi(a) + 0.5d0 * Tf(a,b,g) * Q2f(b,g)
            end do
        end do
    end do

end subroutine get_quadrupole_field

!------------------------------------------------------------------------------

subroutine get_octopole_field(Fi, Ri, Rj, Q3j)

    real(r8), dimension(3), intent(out) :: Fi
    real(r8), dimension(3), intent(in) :: Ri, Rj, Q3j

    integer :: a, b, g, d
    integer, parameter :: k = 4
    real(r8), dimension(3) :: Rij
    real(r8), dimension(15) :: Tij
    real(r8), dimension(3,3,3) :: Q3f
    real(r8), dimension(3,3,3,3) :: Tf

    Rij = Rj - Ri

    call get_Tk_tensor(Tij, k, Rij)

    call get_full_3rd_tensor(Q3f, Q3j)
    call get_full_4th_tensor(Tf, Tij)

    Fi = 0.0d0

    do a = 1, 3
        do b = 1, 3
            do g = 1,3
                do d = 1, 3
                    Fi(a) = Fi(a) + Tf(a,b,g,d) * Q3f(b,g,d) / 6.0d0
                end do
            end do
        end do
    end do

end subroutine get_octopole_field

!------------------------------------------------------------------------------

subroutine get_response_matrix(A, invert, work)

! TODO: Damping schemes
!       Cutoff radius
!       Check if A is allocated

    real(r8), dimension(:), intent(out) :: A
    logical, intent(in) :: invert
    real(r8), dimension(:), intent(inout) :: work

    logical :: exclude, lexist
    integer :: info, lutemp
    integer :: i, j, k, l, m, n
    real(r8) :: R, R3, R5, T
    real(r8), dimension(3) :: Rij
    real(r8), dimension(6) :: alphainv

    A = 0.0d0; R = 0.0d0; alphainv = 0.0d0

    inquire(file='pe_response_matrix.bin', exist=lexist)

    if (lexist) then
        call openfile('pe_response_matrix.bin', lutemp, 'old', 'unformatted')
        rewind(lutemp)
        read(lutemp) A
        close(lutemp)
    else

        m = 0

        do i = 1, ncents

            if (abs(maxval(alphas(:,i))) <= zero) cycle

            alphainv = alphas(:,i)

            call invert_packed_matrix(alphainv, 'p', work)

            do l = 3, 1, -1

                do j = i, ncents

                    if (abs(maxval(alphas(:,j))) <= zero) cycle

                    if (j == i) then

                        if (l == 3) then
                            do k = 1, l
                                A(m+k) = alphainv(k)
                            end do
                        else if (l == 2) then
                            do k = 1, l
                                A(m+k) = alphainv(3+k)
                            end do
                        else if (l == 1) then
                            do k = 1, l
                                A(m+k) = alphainv(5+k)
                            end do
                        end if

                        m = m + l

                    else

! Remove manybody polarization
!                        if (?) then
!                            m = m + 3
!                            cycle
!                        end if

                        exclude = .false.
                        do k = 1, nexlist
                            if (exlist(k,i) == exlist(1,j)) then
                                exclude = .true.
                                exit
                            end if
                        end do

                        if (exclude) then
                            m = m + 3
                            cycle
                        end if

                        Rij = Rs(:,j) - Rs(:,i)
                        R = nrm2(Rij)
                        R3 = R**3
                        R5 = R**5

! Implement cutoff
!                        if (R > cutoff) then
!                            m = m + 3
!                            cycle
!                        end if

                        if (l == 3) then
                            do k = 1, 3
                                T = 3.0d0 * Rij(1) * Rij(k) / R5
                                if (k == 1) T = T - 1.0d0/R3
                                A(m+k) = - T
                            end do
                        else if (l == 2) then
                            do k = 1, 3
                                T = 3.0d0 * Rij(2) * Rij(k) / R5
                                if (k == 2) T = T - 1.0d0/R3
                                A(m+k) = - T
                            end do
                        else if (l == 1) then
                            do k = 1, 3
                                T = 3.0d0 * Rij(3) * Rij(k) / R5
                                if (k == 3) T = T - 1.0d0/R3
                                A(m+k) = - T
                            end do
                        end if

                        m = m + 3

                    end if
                end do
            end do
        end do

        call openfile('pe_response_matrix.bin', lutemp, 'new', 'unformatted')
        rewind(lutemp)
        write(lutemp) A
        close(lutemp)
    end if

    if (invert) then
        call invert_packed_matrix(A, 's', work)
    end if

end subroutine get_response_matrix

!------------------------------------------------------------------------------

subroutine get_Tk_tensor(Tk, k, Rij)

    integer, intent(in) :: k
    real(r8), dimension(:), intent(out) :: Tk
    real(r8), dimension(3), intent(in) :: Rij

    integer :: i, a, b, g, d
    real(r8) :: R, R2, R4

    R = nrm2(Rij)

    if (k == 0) then
        Tk(1) = 1.0d0 / R
    else if (k == 1) then
        do a = 1, 3
            Tk(a) = Rij(a)
        end do
        Tk = - Tk / R**3
    else if (k == 2) then
        R2 = R**2
        i = 1
        do a = 1, 3
            do b = a, 3
                Tk(i) = 3.0d0 * Rij(a) * Rij(b)
                if (a == b) then
                    Tk(i) = Tk(i) - R2
                end if
                i = i + 1
            end do
        end do
!        Tk(1) = 3.0d0 * Rij(1)**2 - R2
!        Tk(2) = 3.0d0 * Rij(1) * Rij(2)
!        Tk(3) = 3.0d0 * Rij(1) * Rij(3)
!        Tk(4) = 3.0d0 * Rij(2)**2 - R2
!        Tk(5) = 3.0d0 * Rij(2) * Rij(3)
!        Tk(6) = 3.0d0 * Rij(3)**2 - R2
        Tk = Tk / R**5
    else if (k == 3) then
        R2 = R**2
        i = 1
        do a = 1, 3
            do b = a, 3
                do g = b, 3
                    Tk(i) = 15.0d0 * Rij(a) * Rij(b) * Rij(g)
                    if (b == g) then
                        Tk(i) = Tk(i) - 3.0d0 * R2 * Rij(a)
                    end if
                    if (a == g) then
                        Tk(i) = Tk(i) - 3.0d0 * R2 * Rij(b)
                    end if
                    if (a == b) then
                        Tk(i) = Tk(i) - 3.0d0 * R2 * Rij(g)
                    end if
                    i = i + 1
                end do
            end do
        end do
!        Tk(1) = 15.0d0 * Rij(1)**3 - 9.0d0 * R2 * Rij(1)
!        Tk(2) = 15.0d0 * Rij(1)**2 * Rij(2) - 3.0d0 * R2 * Rij(2)
!        Tk(3) = 15.0d0 * Rij(1)**2 * Rij(3) - 3.0d0 * R2 * Rij(3)
!        Tk(4) = 15.0d0 * Rij(1) * Rij(2)**2 - 3.0d0 * R2 * Rij(1)
!        Tk(5) = 15.0d0 * Rij(1) * Rij(2) * Rij(3)
!        Tk(6) = 15.0d0 * Rij(1) * Rij(3)**2 - 3.0d0 * R2 * Rij(1)
!        Tk(7) = 15.0d0 * Rij(2)**3 - 9.0d0 * R2 * Rij(2)
!        Tk(8) = 15.0d0 * Rij(2)**2 * Rij(3) - 3.0d0 * R2 * Rij(3)
!        Tk(9) = 15.0d0 * Rij(2) * Rij(3)**2 - 3.0d0 * R2 * Rij(2)
!        Tk(10) = 15.0d0 * Rij(3)**3 - 9.0d0 * R2 * Rij(3)
        Tk = - Tk / R**7
    else if (k == 4) then
        R2 = R**2
        R4 = R**4
        i = 1
        do a = 1, 3
            do b = a, 3
                do g = b, 3
                    do d = g, 3
                        Tk(i) = 105.0d0 * Rij(a) * Rij(b) * Rij(g) * Rij(d)
                        if (g == d) then
                            Tk(i) = Tk(i) - 15.0d0 * R2 * Rij(a) * Rij(b)
                        end if
                        if (b == d) then
                            Tk(i) = Tk(i) - 15.0d0 * R2 * Rij(a) * Rij(g)
                        end if
                        if (b == g) then
                            Tk(i) = Tk(i) - 15.0d0 * R2 * Rij(a) * Rij(d)
                        end if
                        if (a == d) then
                            Tk(i) = Tk(i) - 15.0d0 * R2 * Rij(b) * Rij(g)
                        end if
                        if (a == g) then
                            Tk(i) = Tk(i) - 15.0d0 * R2 * Rij(b) * Rij(d)
                        end if
                        if (a == b) then
                            Tk(i) = Tk(i) - 15.0d0 * R2 * Rij(g) * Rij(d)
                        end if
                        if (a == b .and. g == d) then
                            Tk(i) = Tk(i) + 3.0d0 * R4
                        end if
                        if (a == g .and. b == d) then
                            Tk(i) = Tk(i) + 3.0d0 * R4
                        end if
                        if (a == d .and. b == g) then
                            Tk(i) = Tk(i) + 3.0d0 * R4
                        end if
                        i = i + 1
                    end do
                end do
            end do
        end do
!        Tk(1) = 105.0d0 * Rij(1)**4 - 90.0d0 * R2 * Rij(1)**2 + 9.0d0 * R4
!        Tk(2) = 105.0d0 * Rij(1)**3 * Rij(2) - 45.0d0 * R2 * Rij(1) * Rij(2)
!        Tk(3) = 105.0d0 * Rij(1)**3 * Rij(3) - 45.0d0 * R2 * Rij(1) * Rij(3)
!        Tk(4) = 105.0d0 * Rij(1)**2 * Rij(2)**2 &
!                - 30.0d0 * R2 * (Rij(1)**2 + Rij(2)**2) + 3.0d0 * R4
!        Tk(5) = 105.0d0 * Rij(1)**2 * Rij(2) * Rij(3) &
!                - 15.0d0 * R2 * Rij(2) * Rij(3)
!        Tk(6) = 105.0d0 * Rij(1)**2 * Rij(3)**2 &
!                - 30.0d0 * R2 * (Rij(1)**2 + Rij(3)**2) + 3.0d0 * R4
!        Tk(7) = 105.0d0 * Rij(1) * Rij(2)**3 - 45.0d0 * R2 * Rij(1) * Rij(2)
!        Tk(8) = 105.0d0 * Rij(1) * Rij(2)**2 * Rij(3) &
!                - 15.0d0 * R2 * Rij(1) * Rij(3)
!        Tk(9) = 105.0d0 * Rij(1) * Rij(2) * Rij(3)**2 &
!                - 15.0d0 * R2 * Rij(1) * Rij(2)
!        Tk(10) = 105.0d0 * Rij(1) * Rij(3)**3 - 45.0d0 * R2 * Rij(1) * Rij(3)
!        Tk(11) = 105.0d0 * Rij(2)**4 - 90.0d0 * R2 * Rij(2)**2 + 9.0d0 * R4
!        Tk(12) = 105.0d0 * Rij(2)**3 * Rij(3) - 45.0d0 * R2 * Rij(2) * Rij(3)
!        Tk(13) = 105.0d0 * Rij(2)**2 * Rij(3)**2 &
!                - 30.0d0 * R2 * (Rij(2)**2 + Rij(3)**2) + 3.0d0 * R4
!        Tk(14) = 105.0d0 * Rij(2) * Rij(3)**3 - 45.0d0 * R2 * Rij(2) * Rij(3)
!        Tk(15) = 105.0d0 * Rij(3)**4 - 90.0d0 * R2 * Rij(3)**2 + 9.0d0 * R4
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

    Tf(1,1,1,1) = Ts(1);  Tf(1,1,1,2) = Ts(2);  Tf(1,1,1,3) = Ts(3)
    Tf(1,1,2,1) = Ts(2);  Tf(1,1,2,2) = Ts(4);  Tf(1,1,2,3) = Ts(5)
    Tf(1,1,3,1) = Ts(3);  Tf(1,1,3,2) = Ts(5);  Tf(1,1,3,3) = Ts(6)
    Tf(1,2,1,1) = Ts(2);  Tf(1,2,1,2) = Ts(4);  Tf(1,2,1,3) = Ts(5)
    Tf(1,2,2,1) = Ts(4);  Tf(1,2,2,2) = Ts(7);  Tf(1,2,2,3) = Ts(8)
    Tf(1,2,3,1) = Ts(5);  Tf(1,2,3,2) = Ts(8);  Tf(1,2,3,3) = Ts(9)
    Tf(1,3,1,1) = Ts(3);  Tf(1,3,1,2) = Ts(5);  Tf(1,3,1,3) = Ts(6)
    Tf(1,3,2,1) = Ts(5);  Tf(1,3,2,2) = Ts(8);  Tf(1,3,2,3) = Ts(9)
    Tf(1,3,3,1) = Ts(6);  Tf(1,3,3,2) = Ts(9);  Tf(1,3,3,3) = Ts(10)

    Tf(2,1,1,1) = Ts(2);  Tf(2,1,1,2) = Ts(4);  Tf(2,1,1,3) = Ts(5)
    Tf(2,1,2,1) = Ts(4);  Tf(2,1,2,2) = Ts(7);  Tf(2,1,2,3) = Ts(8)
    Tf(2,1,3,1) = Ts(5);  Tf(2,1,3,2) = Ts(8);  Tf(2,1,3,3) = Ts(9)
    Tf(2,2,1,1) = Ts(4);  Tf(2,2,1,2) = Ts(7);  Tf(2,2,1,3) = Ts(8)
    Tf(2,2,2,1) = Ts(7);  Tf(2,2,2,2) = Ts(11); Tf(2,2,2,3) = Ts(12)
    Tf(2,2,3,1) = Ts(8);  Tf(2,2,3,2) = Ts(12); Tf(2,2,3,3) = Ts(13)
    Tf(2,3,1,1) = Ts(5);  Tf(2,3,1,2) = Ts(8);  Tf(2,3,1,3) = Ts(9)
    Tf(2,3,2,1) = Ts(8);  Tf(2,3,2,2) = Ts(12); Tf(2,3,2,3) = Ts(13)
    Tf(2,3,3,1) = Ts(9);  Tf(2,3,3,2) = Ts(13); Tf(2,3,3,3) = Ts(14)

    Tf(3,1,1,1) = Ts(3);  Tf(3,1,1,2) = Ts(5);  Tf(3,1,1,3) = Ts(6)
    Tf(3,1,2,1) = Ts(5);  Tf(3,1,2,2) = Ts(8);  Tf(3,1,2,3) = Ts(9)
    Tf(3,1,3,1) = Ts(6);  Tf(3,1,3,2) = Ts(9);  Tf(3,1,3,3) = Ts(10)
    Tf(3,2,1,1) = Ts(5);  Tf(3,2,1,2) = Ts(8);  Tf(3,2,1,3) = Ts(9)
    Tf(3,2,2,1) = Ts(8);  Tf(3,2,2,2) = Ts(12); Tf(3,2,2,3) = Ts(13)
    Tf(3,2,3,1) = Ts(9);  Tf(3,2,3,2) = Ts(13); Tf(3,2,3,3) = Ts(14)
    Tf(3,3,1,1) = Ts(6);  Tf(3,3,1,2) = Ts(9);  Tf(3,3,1,3) = Ts(10)
    Tf(3,3,2,1) = Ts(9);  Tf(3,3,2,2) = Ts(13); Tf(3,3,2,3) = Ts(14)
    Tf(3,3,3,1) = Ts(10); Tf(3,3,3,2) = Ts(14); Tf(3,3,3,3) = Ts(15)

end subroutine get_full_4th_tensor

!------------------------------------------------------------------------------

subroutine get_Qk_integrals(Qk_ints, k, Rij, Qk, work)

    external get_Tk_integrals

    integer, intent(in) :: k
    real(r8), dimension(:,:), intent(out) :: Qk_ints
    real(r8), dimension(:), intent(in) :: Qk
    real(r8), dimension(3), intent(in) :: Rij
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
    call get_Tk_integrals(Qk_ints, ncomps*nbas, k, 'packed', Rij,&
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

subroutine invert_packed_matrix(ap, sp, work)

    external dpptrf, dpptri, dsptrf, dsptri, xerbla

    character(*), optional :: sp
    real(r8), dimension(:), intent(inout) :: ap, work

    integer :: n, info
    integer, dimension(:), allocatable :: ipiv

    n = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * real(size(ap), r8)) - 1.0d0))

    if (.not. present(sp)) then
        sp = 'p'
    end if

    if (sp == 'p') then
        call dpptrf('L', n, ap, info)
        if (info /= 0) call xerbla('pptrf', info)
        call dpptri('L', n, ap, info)
        if (info /= 0) call xerbla('pptri', info)
    else if (sp == 's') then
        if (size(work) < n) then
            stop('Not enough work memory to invert matrix.')
        end if
        allocate(ipiv(n))
        call dsptrf('L', n, ap, ipiv, info)
        if (info /= 0) call xerbla('sptrf', info)
        call dsptri('L', n, ap, ipiv, work, info)
        if (info /= 0) call xerbla('sptri', info)
        deallocate(ipiv)
    end if

end subroutine invert_packed_matrix

!------------------------------------------------------------------------------

subroutine solve(ap, b)

    external dspsv

    real(r8), dimension(:), intent(inout) :: ap, b

    integer :: n, info
    integer, dimension(:), allocatable :: ipiv

    n = size(b)

    allocate(ipiv(n))

    info = 0

    call dspsv('L', n, 1, ap, ipiv, b, n, info)

    if (info /= 0) call xerbla('spsv', info)

end subroutine solve

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
