module polarizable_embedding

    implicit none

    private

    public :: pe_dalton_input, pe_read_potential, pe_fock, pe_energy
    public :: pe_frozen_density, pe_intermolecular, pe_repulsion
    public :: pe_response


    ! logicals
    logical, save :: energy = .false.

    ! options
    logical, public, save :: peqm = .false.
    logical, public, save :: pe_twoint = .false.
    logical, public, save :: pe_repuls = .false.
    logical, public, save :: pe_savden = .false.
    logical, public, save :: pe_qmes = .false.
    
    ! logical unit from dalton
    integer, save :: luout = 0

    ! precision
    integer, parameter :: dp = selected_real_kind(15, 200)

    ! constants (codata 2002)
    ! 1 bohr = 0.5291772108 Aa
    real(dp), parameter :: au2aa = 1.8897261249935897d0

    ! thresholds
    real(dp), parameter :: zero = 1.0d-6

    ! polarizable embedding potential info
    ! ------------------------------------

    ! number of sites
    integer, public, save :: nsites = 0
    ! number of polarizable sites
    integer, save :: npols = 0
    ! exclusion list length
    integer, save :: lexlst
    
    ! specifies what type of parameters are present
    ! lmul(0): monopoles, lmul(1): dipoles etc.
    logical, dimension(0:4), public, save :: lmul = .false.
    ! lpol(0): isotropic dipole-dipole polarizabilities
    ! lpol(1): anisotropic dipole-dipole polarizabilities
    logical, dimension(0:1), public, save :: lpol = .false.
    ! lhypol(1): dipole-dipole-dipole polarizabilities/1st hyperpolarizability
!    logical, dimension(1), public, save :: lhypol

    ! elements, coordinates and exclusion lists
    ! nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zs
    ! coordinates
    real(dp), dimension(:,:), allocatable, save :: Rs
    ! exclusion list
    integer, dimension(:,:), allocatable, save :: exlists

    ! energy contributions
    ! electrostatic (last spot reserved for frozen densities)
    real(dp), dimension(0:4), public, save :: Ees = 0.0d0
    ! polarization (1: electron, 2: nuclear, 3: multipole, 4: frozen density)
    real(dp), dimension(4), public, save :: Epol = 0.0d0

    ! multipole moments
    ! monopoles
    real(dp), dimension(:,:), allocatable, save :: Q0
    ! dipoles
    real(dp), dimension(:,:), allocatable, save :: Q1
    ! quadrupoles
    real(dp), dimension(:,:), allocatable, save :: Q2
    ! octopoles
    real(dp), dimension(:,:), allocatable, save :: Q3
    ! hexadecapoles
!    real(dp), dimension(:,:), allocatable, save :: Q4

    ! (hyper)polarizabilities
    ! dipole-dipole polarizabilities
    real(dp), dimension(:,:), allocatable, save :: alphas
    ! .true. if alpha > 0 else .false.
    logical, dimension(:), allocatable, save :: zeroalphas
    ! dipole-dipole-dipole polarizabilities / 1st hyperpolarizabilities
!    real(dp), dimension(:,:), allocatable, save :: betas
    ! dipole-quadrupole polarizabilities
!    real(dp), dimension(:,:), allocatable, save :: As
    ! quadrupole-quadrupole polarizabilities
!    real(dp), dimension(:,:), allocatable, save :: Cs

    ! QM core fragment info
    ! ---------------------

    ! number of density matrices
!    integer :: ndens
    ! number of basis functions and orbitals
!    integer, save :: nbas, norb
    ! number nuclei in qm region
    integer, save :: qmnucs
    ! nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zm
    ! nuclear coordinates
    real(dp), dimension(:,:), allocatable, save :: Rm

    ! frozen density fragment info
    ! ----------------------------

    ! number of frozen densities
    integer, public, save :: nfds
    ! number of nuclei in current frozen density
    integer, public, save :: fdnucs
    ! nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zfd
    ! nuclear coordinates
    real(dp), dimension(:,:), allocatable, save :: Rfd

! TODO:
! find better solution for electric field calculation from frozen densities
! hexadecapoles and higher order polarizabilities
! write list of publications which should be cited
! write output related to QMES
! avoid dimensions as input
! remove double zeroing and unecessary zeroing?
! quit if symmetry or implement symmetry
! Environment properties in herrdn.F
! write to output
! memory checks and handle memory better
! response properties (incl. magnetic)
! AA and AU
! save individual one-electron integrals and reuse
! cutoffs and damping
! memory management
! add error catching
! use optional?
! ddot or dot_product and dgemm or matmul
! parallelization (openMP, MPI, CUDA/openCL)

contains

!------------------------------------------------------------------------------

subroutine pe_dalton_input(word, luinp, lupri)

    external :: upcase

    character(7), intent(inout) :: word
    integer, intent(in) :: luinp
    integer, intent(in) :: lupri

    character(7) :: option

    luout = lupri

    do
        read(luinp,'(a7)') option
        call upcase(option)

        ! do a Polarizable Embedding calculation
        if (trim(option(2:)) == 'PEQM') then
            peqm = .true.
        ! calculate intermolecular two-electron integrals
        else if (trim(option(2:)) == 'TWOINT') then
            read(luinp,*) fdnucs
            pe_twoint = .true.
        ! save frozen density
        else if (trim(option(2:)) == 'SAVDEN') then
            pe_savden = .true.
        ! get fock matrix for repulsion potential
        else if (trim(option(2:)) == 'REPULS') then
            pe_repuls = .true.
        else if (trim(option(2:)) == 'QMES') then
            ! number of frozen densities
            read(luinp,*) nfds
            pe_qmes = .true.
        else if (option(1:1) == '*') then
            word = option
            exit
        else if (option(1:1) == '!' .or. option(1:1) == '#') then
            cycle
        end if
    end do

end subroutine pe_dalton_input

!------------------------------------------------------------------------------

subroutine pe_read_potential(work, coords, charges)

    ! input parameters could be options given in dalton input
    ! so that cutoffs etc. are handled in here.

    real(dp), dimension(:), intent(in), optional :: charges
    real(dp), dimension(:,:), intent(in), optional :: coords
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, k, s
    integer :: lupot, nlines
    real(dp), dimension(15) :: temp
    character(2) :: auoraa
    character(80) :: word

    if (present(coords) .and. present(charges)) then
        qmnucs = size(charges)
        allocate(Rm(3,qmnucs), Zm(1,qmnucs))
        Rm = coords
        Zm(1,:) = charges
    else if (present(coords) .and. .not. present(charges)) then
        print *, 'ERROR: nuclear charges of the QM system are missing.'
        stop
    else if (.not. present(coords) .and. present(charges)) then
        print *, 'ERROR: nuclear coordinates of the QM system are missing.'
    end if

    call openfile('POTENTIAL.INP', lupot, 'old', 'formatted')

    do
        read(lupot,*,end=100) word

        if (trim(word) == 'coordinates') then
            read(lupot,*) nsites
            read(lupot,*) auoraa
            allocate(Zs(1,nsites), Rs(3,nsites))
            do i = 1, nsites
                read(lupot,*) Zs(1,i), (Rs(j,i), j = 1, 3)
            end do
        else if (trim(word) == 'monopoles') then
            lmul(0) = .true.
            allocate(Q0(1,nsites)); Q0 = 0.0d0
            temp = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, temp(1)
                Q0(1,s) = temp(1)
            end do
        else if (trim(word) == 'dipoles') then
            lmul(1) = .true.
            allocate(Q1(3,nsites)); Q1 = 0.0d0
            temp = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 3)
                Q1(:,s) = temp(1:3)
            end do
        else if (trim(word) == 'quadrupoles') then
            lmul(2) = .true.
            allocate(Q2(6,nsites)); Q2 = 0.0d0
            temp = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 6)
                Q2(:,s) = temp(1:6)
            end do
        else if (trim(word) == 'octopoles') then
            lmul(3) = .true.
            allocate(Q3(10,nsites)); Q3 = 0.0d0
            temp = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 10)
                Q3(:,s) = temp(1:10)
            end do
!        else if (trim(word) == 'hexadecapoles') then
!            lmul(4) = .true.
!            allocate(Q4(15,nsites)); Q4 = 0.0d0
!            temp = 0.0d0
!            read(lupot,*) nlines
!            do i = 1, nlines
!                read(lupot,*) s, (temp(j), j = 1, 15)
!                Q4(:,s) = temp(1:15)
!            end do
        else if (trim(word) == 'isoalphas') then
            lpol(0) = .true.
            allocate(alphas(6,nsites)); alphas = 0.0d0
            temp = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, temp(1)
                alphas(1,s) = temp(1)
                alphas(4,s) = temp(1)
                alphas(6,s) = temp(1)
            end do
        else if (trim(word) == 'alphas') then
            lpol(1) = .true.
            if (.not. allocated(alphas)) then
                allocate(alphas(6,nsites)); alphas = 0.0d0
            end if
            temp = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 6)
                alphas(:,s) = temp(1:6)
            end do
!        else if (trim(word) == 'As') then
!            lpol(2) = .true.
!            allocate(As(10,nsites)); As = 0.0d0
!            temp = 0.0d0
!            read(lupot,*) nlines
!            do i = 1, nlines
!                read(lupot,*) s, (temp(j), j = 1, 10)
!                As(:,s) = temp(1:10)
!            end do
!        else if (trim(word) == 'Cs') then
!            lpol(3) = .true.
!            allocate(Cs(15,nsites)); Cs = 0.0d0
!            temp = 0.0d0
!            read(lupot,*) nlines
!            do i = 1, nlines
!                read(lupot,*) s, (temp(j), j = 1, 15)
!                Cs(:,s) = temp(1:15)
!            end do
!        else if (trim(word) == 'betas') then
!            lhypol(1) = .true.
!            allocate(betas(10,nsites)); betas = 0.0d0
!            temp = 0.0d0
!            read(lupot,*) nlines
!            do i = 1, nlines
!                read(lupot,*) s, (temp(j), j = 1, 10)
!                betas(:,s) = temp(1:10)
!            end do
        else if (trim(word) == 'exlists') then
            read(lupot,*) lexlst
            allocate(exlists(lexlst,nsites)); exlists = 0
            do i = 1, nsites
                read(lupot,*) (exlists(j,i), j = 1, lexlst)
            end do
        else if (word(1:1) == '!' .or. word(1:1) == '#') then
            cycle
        end if
    end do

100 continue

    ! if coordinates are in AA then convert to AU
    if (auoraa == 'AA') then
        Rs = Rs * au2aa
    end if

    ! default exclusion list (everything polarizes everything)
    if (.not. allocated(exlists)) then
        allocate(exlists(1,nsites))
        do i = 1, nsites
            exlists(1,i) = i
        end do
    end if

    ! number of polarizabilities different from zero
    npols = 0
    allocate(zeroalphas(nsites))
    do i = 1, nsites
        if (abs(maxval(alphas(:,i))) >= zero) then
            zeroalphas(i) = .false.
            npols = npols + 1
        else
            zeroalphas = .true.
        end if
    end do

    close(lupot)

    ! check memory requirements?

end subroutine pe_read_potential

!------------------------------------------------------------------------------

subroutine pe_fock(density, fock, Epe, work)

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(out) :: fock
    real(dp), intent(out) :: Epe
    real(dp), dimension(:), intent(inout) :: work

    Epe = 0.0d0

    fock = 0.0d0

    call pe_electrostatic(density, fock, work)

    call pe_polarization(density, fock, work)

    Epe = Epe + sum(Ees) + sum(Epol)

end subroutine pe_fock

!------------------------------------------------------------------------------

subroutine pe_energy(density, work)

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(inout) :: work

    energy = .true.

    call pe_electrostatic(density, work(1:size(density)), work(size(density)+1:))

    call pe_polarization(density, work(1:size(density)), work(size(density)+1:))

    energy = .false.

end subroutine pe_energy

!------------------------------------------------------------------------------

subroutine pe_response(density, fock, ndens, work)

    external :: get_Tk_integrals

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(out) :: fock
    integer, intent(in) :: ndens
    real(dp), dimension(:), intent(inout) :: work

    logical :: skip
    integer :: i, j, l, m
    integer :: nnbas
    integer, parameter :: k = 1
    real(dp), dimension(:,:), allocatable :: Fel, Mu
    real(dp), dimension(:,:), allocatable :: Fel_ints, Mu_ints

    nnbas = size(density) / ndens

    fock = 0.0d0

    allocate(Fel(3*npols,ndens)); Fel = 0.0d0
    allocate(Fel_ints(nnbas,3))

    l = 0

    do i = 1, nsites

        if (zeroalphas(i)) cycle

        call get_Tk_integrals(Fel_ints, 3*nnbas, k, Rs(:,i), work, size(work))

        do m = 1, ndens
            do j = 1, 3
                Fel(l+j,m) = dot(density((m-1)*nnbas+1:m*nnbas), Fel_ints(:,j))
            end do
        end do

        l = l + 3

    end do

    deallocate(Fel_ints)

    allocate(Mu(3*npols,ndens))
    allocate(Mu_ints(nnbas,3))

    do m = 1, ndens
        call get_induced_dipoles(Mu(:,m), Fel(:,m), work)
    end do

    deallocate(Fel)

    l = 1
    do i = 1, nsites
        if (zeroalphas(i)) cycle
        do m = 1, ndens
            call get_Qk_integrals(Mu_ints, k, Rs(:,i), Mu(l:l+2,m), work)
            do j = 1, 3
                fock((m-1)*nnbas+1:m*nnbas) = fock((m-1)*nnbas+1:m*nnbas) + Mu_ints(:,j)
            end do
         end do
         l = l + 3
    end do

    deallocate(Mu, Mu_ints)

end subroutine pe_response

!------------------------------------------------------------------------------

subroutine pe_electrostatic(density, fock, work)

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(inout) :: fock
    real(dp), dimension(:), intent(inout) :: work

    logical :: lexist
    integer :: lutemp
    real(dp) :: Eel, Enuc, Esave

    Ees = 0.0d0

    inquire(file='pe_electrostatics.bin', exist=lexist)

    if (lexist .and. .not. energy) then

        call openfile('pe_electrostatics.bin', lutemp, 'old', 'unformatted')
        rewind(lutemp)
        read(lutemp) Esave, fock
        close(lutemp)
        Ees(0) = Esave
        Ees(0) = Ees(0) + dot(density, fock)

    else

        Esave = 0.0d0

        if (lmul(0)) then
            call es_monopoles(density, fock, Eel, Enuc, work)
            Ees(0) = Ees(0) + Eel + Enuc
            Esave = Esave + Enuc
        end if

        if (lmul(1)) then
            call es_dipoles(density, fock, Eel, Enuc, work)
            Ees(1) = Ees(1) + Eel + Enuc
            Esave = Esave + Enuc
        end if

        if (lmul(2)) then
            call es_quadrupoles(density, fock, Eel, Enuc, work)
            Ees(2) = Ees(2) + Eel + Enuc
            Esave = Esave + Enuc
        end if

        if (lmul(3)) then
            call es_octopoles(density, fock, Eel, Enuc, work)
            Ees(3) = Ees(3) + Eel + Enuc
            Esave = Esave + Enuc
        end if

        if (pe_qmes) then
            call es_frozen_densities(density, fock, Eel, Enuc, work)
            Ees(4) = Ees(4) + Eel + Enuc
            Esave = Esave + Enuc
        end if

        if (.not. energy) then
            call openfile('pe_electrostatics.bin', lutemp, 'new', 'unformatted')
            rewind(lutemp)
            write(lutemp) Esave, fock
            close(lutemp)
        end if

    end if

end subroutine pe_electrostatic

!------------------------------------------------------------------------------

subroutine es_frozen_densities(density, fock, Eel, Enuc, work)

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(inout) :: fock
    real(dp), intent(out) :: Eel, Enuc
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l
    integer, parameter :: k = 0
    integer :: lufck, lexist, lutemp
    real(dp) :: Ene, Een, Eee, Enn
    real(dp), dimension(1) :: Tfm
    real(dp), dimension(3) :: Rfm
    real(dp), dimension(:), allocatable :: fd_fock
    real(dp), dimension(:,:), allocatable :: Zfd_ints
    character(99) :: ci
    character(99) :: filename

    Eel = 0.0d0; Enuc = 0.0d0
    Ene = 0.0d0; Een = 0.0d0; Eee = 0.0d0; Enn = 0.0d0

    allocate(fd_fock(size(density)), Zfd_ints(size(density),1))

    do i = 1, nfds

        write(ci,*) i
        ci = adjustl(ci)
        filename = 'pe_fock_'//trim(ci)//'.bin'
        call openfile(trim(filename), lufck, 'old', 'unformatted')
        rewind(lufck)
        read(lufck) work(1:3*npols)
        read(lufck) Ene
        read(lufck) fd_fock
        read(lufck) fdnucs
        allocate(Rfd(3,fdnucs), Zfd(1,fdnucs))
        read(lufck) Rfd, Zfd
        close(lufck)

        fock = fock + fd_fock

        Eee = dot(density, fd_fock)

        do j = 1, fdnucs
            call get_Qk_integrals(Zfd_ints, k, Rfd(:,j), Zfd(:,j), work)
            Een = Een + dot(density, Zfd_ints(:,1))
            fock = fock + Zfd_ints(:,1)
        end do

        do j = 1, fdnucs
            do l = 1, qmnucs
                Rfm = Rm(:,l) - Rfd(:,j)
                call get_Tk_tensor(Tfm, k, Rfm)
                Enn = Enn + Zm(1,l) * Zfd(1,j) * Tfm(1)
            end do
        end do

        deallocate(Rfd, Zfd)

        Enuc = Enuc + Ene + Enn
        Eel = Eel + Een + Eee

    end do

    deallocate(fd_fock, Zfd_ints)

end subroutine es_frozen_densities

!------------------------------------------------------------------------------

subroutine es_monopoles(density, fock, Eel, Enuc, work)

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(inout) :: fock
    real(dp), intent(out) :: Eel, Enuc
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j
    integer, parameter :: k = 0
    real(dp), dimension(3) :: Rsm
    real(dp), dimension(1) :: Tsm
    real(dp), dimension(:,:), allocatable :: Q0_ints

    allocate(Q0_ints(size(density),1)); Q0_ints = 0.0d0

    Eel = 0.0d0; Enuc = 0.0d0

    do i = 1, nsites

        if (abs(Q0(1,i)) < zero) cycle

        call get_Qk_integrals(Q0_ints, k, Rs(:,i), Q0(:,i), work)

        ! nuclei - monopole interaction
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call get_Tk_tensor(Tsm, k, Rsm)
            Enuc = Enuc + Q0(1,i) * Zm(1,j) * Tsm(1)
        end do

        ! electron - monopole interaction
        Eel = Eel + dot(density, Q0_ints(:,1))
        if (.not. energy) then
            fock = fock + Q0_ints(:,1)
        end if

    end do

    deallocate(Q0_ints)

end subroutine es_monopoles

!------------------------------------------------------------------------------

subroutine es_dipoles(density, fock, Eel, Enuc, work)

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(inout) :: fock
    real(dp), intent(out) :: Eel, Enuc
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l
    integer, parameter :: k = 1
    real(dp), dimension(3) :: Rsm, Tsm
    real(dp), dimension(:,:), allocatable :: Q1_ints

    allocate(Q1_ints(size(density),3)); Q1_ints = 0.0d0

    Eel = 0.0d0; Enuc = 0.0d0

    do i = 1, nsites

        if (abs(maxval(Q1(:,i))) < zero) cycle

        call get_Qk_integrals(Q1_ints, k, Rs(:,i), Q1(:,i), work)

        ! nuclei - dipole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call get_Tk_tensor(Tsm, k, Rsm)
            do l = 1, 3
                Enuc = Enuc - Zm(1,j) * Q1(l,i) * Tsm(l)
            end do
        end do

        ! electron - dipole interaction
        do j = 1, 3
            Eel = Eel + dot(density, Q1_ints(:,j))
            if (.not. energy) then
                fock = fock + Q1_ints(:,j)
            end if
        end do

    end do

    deallocate(Q1_ints)

end subroutine es_dipoles

!------------------------------------------------------------------------------

subroutine es_quadrupoles(density, fock, Eel, Enuc, work)

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(inout) :: fock
    real(dp), intent(out) :: Eel, Enuc
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l
    integer, parameter :: k = 2
    real(dp), dimension(3) :: Rsm
    real(dp), dimension(6) :: Tsm, factors
    real(dp), dimension(:,:), allocatable :: Q2_ints

    allocate(Q2_ints(size(density),6)); Q2_ints = 0.0d0

    Eel = 0.0d0; Enuc = 0.0d0

    do i = 1, nsites

        if (abs(maxval(Q2(:,i))) < zero) cycle

        call get_Qk_integrals(Q2_ints, k, Rs(:,i), Q2(:,i), work)

        ! nuclei - quadrupole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call get_Tk_tensor(Tsm, k, Rsm)
            call get_symmetry_factors(factors, k)
            do l = 1, 6
                Enuc = Enuc + 0.5d0 * factors(l) * Zm(1,j) * Q2(l,i) * Tsm(l)
            end do

        end do

        ! electron - quadrupole interaction energy
        do j = 1, 6
            Eel = Eel + dot(density, Q2_ints(:,j))
            if (.not. energy) then
                fock = fock + Q2_ints(:,j)
            end if
        end do

    end do

    deallocate(Q2_ints)

end subroutine es_quadrupoles

!------------------------------------------------------------------------------

subroutine es_octopoles(density, fock, Eel, Enuc, work)

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(inout) :: fock
    real(dp), intent(out) :: Eel, Enuc
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l
    integer, parameter :: k = 3
    real(dp), dimension(3) :: Rsm
    real(dp), dimension(10) :: Tsm, factors
    real(dp), dimension(:,:), allocatable :: Q3_ints

    allocate(Q3_ints(size(density),10)); Q3_ints = 0.0d0

    Eel = 0.0d0; Enuc = 0.0d0

    do i = 1, nsites

        if (abs(maxval(Q3(:,i))) < zero) cycle

        call get_Qk_integrals(Q3_ints, k, Rs(:,i), Q3(:,i), work)

        ! nuclei - octopole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call get_Tk_tensor(Tsm, k, Rsm)
            call get_symmetry_factors(factors, k)
            do l = 1, 10
                Enuc = Enuc - factors(l) * Zm(1,j) * Q3(l,i) * Tsm(l) / 6.0d0
            end do
        end do

        ! electron - octopole interaction energy
        do j = 1, 10
            Eel = Eel + dot(density, Q3_ints(:,j))
            if (.not. energy) then
                fock = fock + Q3_ints(:,j)
            end if
        end do

    end do

    deallocate(Q3_ints)

end subroutine es_octopoles

!------------------------------------------------------------------------------

subroutine pe_polarization(density, fock, work)

    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(inout) :: fock
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l
    integer, parameter :: k = 1
    real(dp), dimension(:), allocatable :: Mu, Fel, Fnuc, Fmul, Ffd
    real(dp), dimension(:,:), allocatable :: Mu_ints

    Epol = 0.0d0

    if (lpol(0) .or. lpol(1)) then

        allocate(Mu(3*npols)); Mu = 0.0d0
        allocate(Fel(3*npols)); Fel = 0.0d0
        allocate(Fnuc(3*npols)); Fnuc = 0.0d0
        allocate(Fmul(3*npols)); Fmul = 0.0d0
        allocate(Ffd(3*npols)); Ffd = 0.0d0

        if (energy) then

            call get_electron_fields(Fel, density, work)
            call get_nuclear_fields(Fnuc, work)
            call get_multipole_fields(Fmul, work)
            if (pe_qmes) call get_frozen_density_field(Ffd, work)

            if (pe_qmes) then
                call get_induced_dipoles(Mu, Fel + Fnuc + Fmul + Ffd, work)
            else
                call get_induced_dipoles(Mu, Fel + Fnuc + Fmul, work)
            end if

            Epol(1) = - 0.5d0 * dot(Mu, Fel)
            Epol(2) = - 0.5d0 * dot(Mu, Fnuc)
            Epol(3) = - 0.5d0 * dot(Mu, Fmul)
            if (pe_qmes) Epol(4) = - 0.5d0 * dot(Mu, Ffd)

        else

            allocate(Mu_ints(size(fock),3))

            call get_electron_fields(Fel, density, work)
            call get_nuclear_fields(Fnuc, work)
            call get_multipole_fields(Fmul, work)
            if (pe_qmes) call get_frozen_density_field(Ffd, work)

            if (pe_qmes) then
                call get_induced_dipoles(Mu, Fel + Fnuc + Fmul + Ffd, work)
            else
                call get_induced_dipoles(Mu, Fel + Fnuc + Fmul, work)
            end if

            l = 1
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                call get_Qk_integrals(Mu_ints, k, Rs(:,i), Mu(l:l+2), work)
                do j = 1, 3
                    fock = fock + Mu_ints(:,j)
                end do
                l = l + 3
            end do

            Epol(1) = - 0.5d0 * dot(Mu, Fel)
            Epol(2) = - 0.5d0 * dot(Mu, Fnuc)
            Epol(3) = - 0.5d0 * dot(Mu, Fmul)
            if (pe_qmes) Epol(4) = - 0.5d0 * dot(Mu, Ffd)

            deallocate(Mu_ints)

        end if

        deallocate(Mu, Fel, Fnuc, Fmul, Ffd)

    end if

end subroutine pe_polarization

!------------------------------------------------------------------------------

subroutine get_induced_dipoles(Mu, F, work)

    real(dp), dimension(:), intent(out) :: Mu
    real(dp), dimension(:), intent(in) :: F
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, n
    real(dp), dimension(:), allocatable :: A

    n = size(Mu)

    allocate(A(n*(n+1)/2))
    
    call get_response_matrix(A, .false., work)

    Mu = F

    call solve(A, Mu)

    deallocate(A)

end subroutine get_induced_dipoles

!------------------------------------------------------------------------------

subroutine get_electric_fields(Ftot, density, work)

    real(dp), dimension(:), intent(out) :: Ftot
    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j
    integer :: n
    real(dp), dimension(:), allocatable :: Fel, Fnuc, Fmul

! TODO: electric field from frozen densities

    n = size(Ftot)

    allocate(Fel(n), Fnuc(n), Fmul(n)); Fel = 0.0d0; Fnuc = 0.0d0; Fmul = 0.0d0

    call get_electron_fields(Fel, density, work)

    call get_nuclear_fields(Fnuc, work)

    call get_multipole_fields(Fmul, work)

    Ftot = Fel + Fnuc + Fmul

end subroutine get_electric_fields

!------------------------------------------------------------------------------

subroutine get_electron_fields(Fel, density, work)

    external :: get_Tk_integrals

    real(dp), dimension(:), intent(out) :: Fel
    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(inout) :: work

    logical :: skip
    integer :: i, j, l
    integer :: nnbas
    integer, parameter :: k = 1
    real(dp), dimension(:,:), allocatable :: Fel_ints

    Fel = 0.0d0

    nnbas = size(density)

    allocate(Fel_ints(nnbas,3))

    l = 0

    do i = 1, nsites

        if (pe_savden) then
! TODO: finde better threshold or better solution
            skip = .false.
            do j = 1, qmnucs
                if (nrm2(Rs(:,i) - Rm(:,j)) <= 1.0d0) skip = .true.
            end do
            if (skip) cycle
        end if

        if (zeroalphas(i)) cycle

        call get_Tk_integrals(Fel_ints, 3*nnbas, k, Rs(:,i), work, size(work))

        do j = 1, 3
            Fel(l+j) = dot(density, Fel_ints(:,j))
        end do

        l = l + 3

    end do

    deallocate(Fel_ints)

end subroutine get_electron_fields

!------------------------------------------------------------------------------

subroutine get_nuclear_fields(Fnuc, work)

    real(dp), dimension(:), intent(out) :: Fnuc
    real(dp), dimension(:), intent(inout) :: work

    logical :: exclude, lexist, skip
    integer :: lutemp
    integer :: i, j, l, m
    integer, parameter :: k = 1
    real(dp), dimension(3) :: Rms, Tms

    Fnuc = 0.0d0

    inquire(file='pe_nuclear_field.bin', exist=lexist)

    if (lexist) then

        call openfile('pe_nuclear_field.bin', lutemp, 'old', 'unformatted')
        rewind(lutemp)
        read(lutemp) Fnuc
        close(lutemp)

    else

        l = 0
        do i = 1, nsites

            if (pe_savden) then
! TODO: finde better threshold or better solution
                skip = .false.
                do j = 1, qmnucs
                    if (nrm2(Rs(:,i) - Rm(:,j)) <= 1.0d0) skip = .true.
                end do
                if (skip) cycle
            end if

            if (zeroalphas(i)) cycle
            do j = 1, qmnucs
                Rms = Rs(:,i) - Rm(:,j)
                call get_Tk_tensor(Tms, k, Rms)
                do m = 1, 3
                    Fnuc(l+m) = Fnuc(l+m) - Zm(1,j) * Tms(m)
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

subroutine get_frozen_density_field(Ffd, work)

    real(dp), dimension(:), intent(out) :: Ffd
    real(dp), dimension(:), intent(inout) :: work

    integer :: i
    integer :: lutemp
    character(99) :: ci
    character(80) :: filename
    real(dp), dimension(3*npols) :: Ftmp

    Ffd = 0.0d0

! TODO: nuclear field???
    do i = 1, nfds

        Ftmp = 0.0d0
        write(ci,*) i
        ci = adjustl(ci)
        filename = 'pe_fock_'//trim(ci)//'.bin'
        call openfile(trim(filename), lutemp, 'old', 'unformatted')
        rewind(lutemp)
        read(lutemp) Ftmp
        close(lutemp)

        Ffd = Ffd + Ftmp

    end do

end subroutine get_frozen_density_field

!------------------------------------------------------------------------------

subroutine get_multipole_fields(Fmul, work)

    real(dp), dimension(:), intent(out) :: Fmul
    real(dp), dimension(:), intent(inout) :: work

    logical :: exclude, lexist
    integer :: lutemp
    integer :: i, j, k, l
    real(dp) :: Rij
    real(dp), dimension(3) :: Fs

    Fmul = 0.0d0

    inquire(file='pe_multipole_field.bin', exist=lexist)

    if (lexist) then

        call openfile('pe_multipole_field.bin', lutemp, 'old', 'unformatted')
        rewind(lutemp)
        read(lutemp) Fmul
        close(lutemp)

    else

        l = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            do j = 1, nsites
                ! check if j is allowed to polarize i
                if (i == j) cycle
                exclude = .false.
                do k = 1, lexlst
                    if (exlists(k,i) == exlists(1,j)) then
                        exclude = .true.
                        exit
                    end if
                end do
                if (exclude) cycle

                ! cutoff
!                Rij = Rs(:,j) - Rs(:,j)

                ! get electric field at i due to monopole at j
                if (lmul(0)) then
                    ! skip if monopole is 'zero'
                    if (abs(maxval(Q0(:,j))) >= zero) then
                        call get_monopole_field(Fs, Rs(:,i), Rs(:,j), Q0(:,j))
                        Fmul(l:l+2) = Fmul(l:l+2) + Fs
                    end if
                end if

                ! get electric field at i due to dipole at j
                if (lmul(1)) then
                    ! skip if dipole is 'zero'
                    if (abs(maxval(Q1(:,j))) >= zero) then
                        call get_dipole_field(Fs, Rs(:,i), Rs(:,j), Q1(:,j))
                        Fmul(l:l+2) = Fmul(l:l+2) + Fs
                    end if
                end if

                ! get electric field at i due to quadrupole at j
                if (lmul(2)) then
                    ! skip if quadrupole is 'zero'
                    if (abs(maxval(Q2(:,j))) >= zero) then
                        call get_quadrupole_field(Fs, Rs(:,i), Rs(:,j), Q2(:,j))
                        Fmul(l:l+2) = Fmul(l:l+2) + Fs
                    end if
                end if

                ! get electric field at i due to octopole at j
                if (lmul(3)) then
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

    real(dp), dimension(3), intent(out) :: Fi
    real(dp), dimension(3), intent(in) :: Ri, Rj
    real(dp), dimension(1), intent(in) :: Q0j

    integer :: a
    integer, parameter :: k = 1
    real(dp), dimension(3) :: Rji, Tji

    Rji = Ri - Rj

    call get_Tk_tensor(Tji, k, Rji)

    Fi = 0.0d0

    do a = 1, 3
        Fi(a) = Fi(a) - Tji(a) *  Q0j(1)
    end do

end subroutine get_monopole_field

!------------------------------------------------------------------------------

subroutine get_dipole_field(Fi, Ri, Rj, Q1j)

    real(dp), dimension(3), intent(out) :: Fi
    real(dp), dimension(3), intent(in) :: Ri, Rj, Q1j

    integer :: a, b
    integer, parameter :: k = 2
    real(dp), dimension(3) :: Rji
    real(dp), dimension(6) :: Tji
    real(dp), dimension(3,3) :: Tf

    Rji = Ri - Rj

    call get_Tk_tensor(Tji, k, Rji)

    call get_full_2nd_tensor(Tf, Tji)

    Fi = 0.0d0

    do a = 1, 3
        do b = 1, 3
            Fi(a) = Fi(a) + Tf(a,b) * Q1j(b)
        end do
    end do

end subroutine get_dipole_field

!------------------------------------------------------------------------------

subroutine get_quadrupole_field(Fi, Ri, Rj, Q2j)

    real(dp), dimension(3), intent(out) :: Fi
    real(dp), dimension(3), intent(in) :: Ri, Rj, Q2j

    integer :: a, b, g
    integer, parameter :: k = 3
    real(dp), dimension(3) :: Rji
    real(dp), dimension(10) :: Tji
    real(dp), dimension(3,3) :: Q2f
    real(dp), dimension(3,3,3) :: Tf

    Rji = Ri - Rj

    call get_Tk_tensor(Tji, k, Rji)

    call get_full_2nd_tensor(Q2f, Q2j)
    call get_full_3rd_tensor(Tf, Tji)

    Fi = 0.0d0

    do a = 1, 3
        do b = 1, 3
            do g = 1, 3
                Fi(a) = Fi(a) - 0.5d0 * Tf(a,b,g) * Q2f(b,g)
            end do
        end do
    end do

end subroutine get_quadrupole_field

!------------------------------------------------------------------------------

subroutine get_octopole_field(Fi, Ri, Rj, Q3j)

    real(dp), dimension(3), intent(out) :: Fi
    real(dp), dimension(3), intent(in) :: Ri, Rj, Q3j

    integer :: a, b, g, d
    integer, parameter :: k = 4
    real(dp), dimension(3) :: Rji
    real(dp), dimension(15) :: Tji
    real(dp), dimension(3,3,3) :: Q3f
    real(dp), dimension(3,3,3,3) :: Tf

    Rji = Ri - Rj

    call get_Tk_tensor(Tji, k, Rji)

    call get_full_3rd_tensor(Q3f, Q3j)
    call get_full_4th_tensor(Tf, Tji)

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

    real(dp), dimension(:), intent(out) :: A
    logical, intent(in) :: invert
    real(dp), dimension(:), intent(inout) :: work

    logical :: exclude, lexist
    integer :: info, lutemp
    integer :: i, j, k, l, m, n
    real(dp) :: R, R3, R5, T
    real(dp), dimension(3) :: Rij
    real(dp), dimension(6) :: alphainv

    A = 0.0d0; R = 0.0d0; alphainv = 0.0d0

    inquire(file='pe_response_matrix.bin', exist=lexist)

    if (lexist) then
        call openfile('pe_response_matrix.bin', lutemp, 'old', 'unformatted')
        rewind(lutemp)
        read(lutemp) A
        close(lutemp)
    else

        m = 0

        do i = 1, nsites

            if (zeroalphas(i)) cycle

            alphainv = alphas(:,i)

            call invert_packed_matrix(alphainv, 'p', work)

            do l = 3, 1, -1

                do j = i, nsites

                    if (zeroalphas(j)) cycle

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
                        do k = 1, lexlst
                            if (exlists(k,i) == exlists(1,j)) then
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
    real(dp), dimension(:), intent(out) :: Tk
    real(dp), dimension(3), intent(in) :: Rij

    integer :: i, a, b, g, d
    real(dp) :: R, R2, R4

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
        Tk = Tk / R**9
    end if

end subroutine get_Tk_tensor

!------------------------------------------------------------------------------

subroutine get_full_2nd_tensor(Tf, Ts)

    real(dp), dimension(:), intent(in) :: Ts
    real(dp), dimension(3,3), intent(out) :: Tf

    Tf = 0.0d0

    Tf(1,1) = Ts(1); Tf(1,2) = Ts(2); Tf(1,3) = Ts(3)
    Tf(2,1) = Ts(2); Tf(2,2) = Ts(4); Tf(2,3) = Ts(5)
    Tf(3,1) = Ts(3); Tf(3,2) = Ts(5); Tf(3,3) = Ts(6)

end subroutine get_full_2nd_tensor

!------------------------------------------------------------------------------

subroutine get_full_3rd_tensor(Tf, Ts)

    real(dp), dimension(:), intent(in) :: Ts
    real(dp), dimension(3,3,3), intent(out) :: Tf

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

    real(dp), dimension(:), intent(in) :: Ts
    real(dp), dimension(3,3,3,3), intent(out) :: Tf

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

    external :: get_Tk_integrals

    integer, intent(in) :: k
    real(dp), dimension(:,:), intent(out) :: Qk_ints
    real(dp), dimension(:), intent(in) :: Qk
    real(dp), dimension(3), intent(in) :: Rij
    real(dp), dimension(:), intent(inout) :: work

    integer :: i
    integer :: ncomps, nnbas
    real(dp) :: taylor
    real(dp), dimension(:), allocatable :: factors

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

    nnbas = size(Qk_ints, 1)
    ncomps = size(Qk_ints, 2)

    ! get T^(k) integrals (incl. negative sign from electron density)
    call get_Tk_integrals(Qk_ints, ncomps*nnbas, k, Rij, work, size(work))

    ! get symmetry factors
    allocate(factors(ncomps)); factors = 0.0d0
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
    real(dp), dimension(:), intent(out) :: factors

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

subroutine pe_frozen_density(density, nbas, coords, charges, work)

    external :: get_Tk_integrals

    integer, intent(in) :: nbas
    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(in) :: charges
    real(dp), dimension(:,:), intent(in) :: coords
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l
    integer :: nnbas, corenucs
    integer, parameter :: k = 0
    integer :: lucore, luden
    character(2) :: auoraa
    real(dp) :: Ene
    real(dp), dimension(:,:), allocatable :: Rc, Zc
    real(dp), dimension(:,:), allocatable :: full_density
    real(dp), dimension(:), allocatable :: T0_ints
    real(dp), dimension(:), allocatable :: Ffd, Ftmp

    ! frozen density nuclear charges and coordinates
    qmnucs = size(charges)
    allocate(Rm(3,qmnucs), Zm(1,qmnucs))
    Rm = coords
    Zm(1,:) = charges
    
    ! read in information about qm core
    call openfile('core.dat', lucore, 'old', 'formatted')
    rewind(lucore)
    read(lucore,*) auoraa
    read(lucore,*) corenucs
    allocate(Zc(1,corenucs), Rc(3,corenucs))
    do i = 1, corenucs
        read(lucore,*) Zc(1,i), (Rc(j,i), j = 1, 3)
    end do
    close(lucore)
    if (auoraa == 'AA') then
        Rc = Rc * au2aa
    end if

    ! unfold density matrix
    allocate(full_density(nbas,nbas)); full_density = 0.0d0
!    call dunfld(nbas, density, full_density)
    l = 1
    do i = 1, nbas
        do j = 1, i
            if (j == i) then
                full_density(i,j) = density(l)
            else
                full_density(i,j) = 0.5d0 * density(l)
                full_density(j,i) = 0.5d0 * density(l)
            end if
            l = l + 1
        end do
    end do

    ! get electric field from frozen density at polarizable sites
    allocate(Ftmp(3*npols), Ffd(3*npols)); Ftmp = 0.0d0
    call get_electron_fields(Ftmp, density, work)
    Ffd = Ftmp

    l = 1
    do i = 1, npols
        print *, i, Ftmp(l:l+2)
        l = l + 3
    end do
    call get_nuclear_fields(Ftmp, work)
    Ffd = Ffd + Ftmp
    l = 1
    do i = 1, npols
        print *, i, Ftmp(l:l+2)
        l = l + 3
    end do


    ! calculate nuclear - electron energy contribution
    nnbas = nbas*(nbas+1)/2
    allocate(T0_ints(nnbas)); T0_ints = 0.0d0
    Ene = 0.0d0
    do i = 1, corenucs
        call get_Tk_integrals(T0_ints, nnbas, k, Rc(:,i), work, size(work))
        T0_ints = Zc(1,i) * T0_ints
        Ene = Ene + dot(density, T0_ints)
    end do
    deallocate(T0_ints)

    ! save density and energy for subsequent calculations
    call openfile('pe_density.bin', luden, 'new', 'unformatted')
    rewind(luden)
    write(luden) Ene
    write(luden) qmnucs
    write(luden) Rm, Zm
    write(luden) npols
    write(luden) Ffd
    write(luden) nbas
    write(luden) full_density
    close(luden)

    deallocate(full_density, Ffd, Ftmp)

end subroutine pe_frozen_density

!------------------------------------------------------------------------------

subroutine pe_intermolecular(nbas, work)

    external :: sirfck

    integer, intent(in) :: nbas
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, k, l
    integer :: fbas, cbas
    integer :: luden, lufck
    integer, dimension(1) :: isymdm, ifctyp
    real(dp) :: Ene
    real(dp), dimension(:), allocatable :: core_fock, Ffd
    real(dp), dimension(:,:), allocatable :: frozen_density, full_density
    real(dp), dimension(:,:), allocatable :: full_fock


    call openfile('pe_density.bin', luden, 'old', 'unformatted')
    rewind(luden)
    read(luden) Ene
    read(luden) fdnucs
    allocate(Rfd(3,fdnucs), Zfd(1,fdnucs))
    read(luden) Rfd, Zfd
    read(luden) npols
    allocate(Ffd(3*npols))
    read(luden) Ffd
    read(luden) fbas
    allocate(frozen_density(fbas, fbas))
    read(luden) frozen_density
    close(luden)

    cbas = nbas - fbas

    ! resize density matrix to full size and fill with frozen density in first
    ! block
    allocate(full_density(nbas,nbas)); full_density = 0.0d0
    full_density(1:fbas,1:fbas) = frozen_density

    ! get two-electron part of Fock matrix using resized density matrix
    allocate(full_fock(nbas,nbas)); full_fock = 0.0d0
!     IFCTYP = +/-XY
!       X indicates symmetry about diagonal
!         X = 0 No symmetry
!         X = 1 Symmetric
!         X = 2 Anti-symmetric
!       Y indicates contributions
!         Y = 0 no contribution !
!         Y = 1 Coulomb
!         Y = 2 Exchange
!         Y = 3 Coulomb + Exchange
!       + sign: alpha + beta matrix (singlet)
!       - sign: alpha - beta matrix (triplet)
!     sirfck(fock, density, ?, isymdm, ifctyp, direct, work, nwrk)
    isymdm = 1
    ifctyp = 11
    call sirfck(full_fock, full_density, 1, isymdm, ifctyp, .false.,&
                work, size(work))

    deallocate(full_density)

    ! extract upper triangle part of full Fock matrix corresponding to
    ! core fragment
    allocate(core_fock(cbas*(cbas+1)/2)); core_fock = 0.0d0
    l = 1
    do i = fbas + 1, nbas
        do j = fbas + 1, i
            core_fock(l) = full_fock(j,i)
            l = l + 1
        end do
    end do

    deallocate(full_fock)

    ! save core Fock matrix
    call openfile('pe_fock.bin', lufck, 'new', 'unformatted')
    rewind(lufck)
    write(lufck) Ffd
    write(lufck) Ene
    write(lufck) core_fock
    write(lufck) fdnucs
    write(lufck) Rfd, Zfd
    close(lufck)

    deallocate(core_fock, Rfd, Zfd, Ffd)

end subroutine pe_intermolecular

!------------------------------------------------------------------------------

subroutine pe_repulsion()

end subroutine pe_repulsion

!------------------------------------------------------------------------------

function nrm2(x)

    real(dp), external :: dnrm2

    real(dp) :: nrm2
    real(dp), dimension(:), intent(in) :: x

    nrm2 = dnrm2(size(x), x, 1)

end function nrm2

!------------------------------------------------------------------------------

function dot(x,y)

    real(dp), external :: ddot

    real(dp) :: dot
    real(dp), dimension(:), intent(in) :: x, y

    dot = ddot(size(x), x, 1, y, 1)

end function dot

!------------------------------------------------------------------------------

subroutine axpy(x, y, a)

    external :: daxpy

    real(dp), optional :: a
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(inout) :: y

    if (.not. present(a)) a = 1.0d0

    call daxpy(size(x), a, x, 1, y, 1)

end subroutine axpy

!------------------------------------------------------------------------------

subroutine gemm(a, b, c, transa, transb, alpha, beta)

    external :: dgemm

    real(dp), intent(in), optional :: alpha, beta
    character(1), intent(in), optional :: transa, transb
    real(dp), dimension(:,:), intent(in) :: a, b
    real(dp), dimension(:,:) , intent(inout) :: c

    integer :: m, n, k, lda, ldb, ldc
    character(1) :: o_transa, o_transb
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

    if (transa == 'N') then
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

subroutine invert_packed_matrix(ap, sp, work)

    external :: dpptrf, dpptri, dsptrf, dsptri, xerbla

    character(*), optional :: sp
    real(dp), dimension(:), intent(inout) :: ap, work

    integer :: n, info
    integer, dimension(:), allocatable :: ipiv

    n = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * real(size(ap), dp)) - 1.0d0))

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
        allocate(ipiv(n)); ipiv = 0.0d0
        call dsptrf('L', n, ap, ipiv, info)
        if (info /= 0) call xerbla('sptrf', info)
        call dsptri('L', n, ap, ipiv, work, info)
        if (info /= 0) call xerbla('sptri', info)
        deallocate(ipiv)
    end if

end subroutine invert_packed_matrix

!------------------------------------------------------------------------------

subroutine solve(ap, b)

    external :: dspsv

    real(dp), dimension(:), intent(inout) :: ap, b

    integer :: n, info
    integer, dimension(:), allocatable :: ipiv

    n = size(b)

    allocate(ipiv(n)); ipiv = 0.0d0

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

    do i = 21, 99
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
