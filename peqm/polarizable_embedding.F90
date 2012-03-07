module polarizable_embedding

#ifdef VAR_MPI
    use mpi
#endif

    implicit none

    private

    intrinsic :: present, size, cpu_time

    public :: pe_dalton_input, pe_read_potential, pe_master
    public :: pe_save_density, pe_intmol_twoints, pe_repulsion

    ! options
    logical, public, save :: peqm = .false.
    logical, public, save :: pe_damp = .false.
    logical, public, save :: pe_gspol = .false.
    logical, public, save :: pe_nomb = .false.
    logical, public, save :: pe_twoint = .false.
    logical, public, save :: pe_repuls = .false.
    logical, public, save :: pe_savden = .false.
    logical, public, save :: pe_fd = .false.
    logical, public, save :: pe_timing = .false.

    ! calculation type
    logical :: fock = .false.
    logical :: energy = .false.
    logical :: response = .false.

#ifdef VAR_MPI
    ! MPI stuff
    public :: pe_mpi
    logical, save :: initialized = .false.
    integer, dimension(:), save, allocatable :: ndists, displs
#endif

    ! logical unit from dalton
    integer, save :: luout = 0

    ! precision
    integer, parameter :: dp = selected_real_kind(15, 200)

    ! constants (codata 2002)
    ! 1 bohr = 0.5291772108 Aa
    real(dp), parameter :: aa2au = 1.8897261249935897d0

    ! thresholds
    real(dp), parameter :: zero = 1.0d-6

    ! damping parameter
    real(dp) :: damp = 2.1304d0

    ! variables used for timings
    real(dp) :: t1, t2

    ! polarizable embedding potential info
    ! ------------------------------------

    ! number of sites
    integer, public, save :: nsites = 0
    ! number of polarizable sites
    integer, save :: npols = 0
    ! loop variable (nloop = nsites if serial)
    integer, save :: nloop = 0
    ! exclusion list length
    integer, save :: lexlst = 0
    
    ! specifies what type of parameters are present
    ! lmul(0): monopoles, lmul(1): dipoles etc.
    logical, dimension(0:3), public, save :: lmul = .false.
    ! lpol(0): isotropic dipole-dipole polarizabilities
    ! lpol(1): anisotropic dipole-dipole polarizabilities
    logical, dimension(0:1), public, save :: lpol = .false.
    ! lhypol(1): dipole-dipole-dipole polarizabilities/1st hyperpolarizability
!    logical, dimension(1), public, save :: lhypol

    ! nuclear charges, coordinates, elements and exclusion lists
    ! nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zs
    ! coordinates
    real(dp), dimension(:,:), allocatable, save :: Rs
    ! elements
    character(len=2), dimension(:,:), allocatable, save :: elems
    ! exclusion list
    integer, dimension(:,:), allocatable, save :: exlists

    ! energy contributions
    ! electrostatic
    real(dp), dimension(:,:), allocatable, public, save :: Ees
    ! polarization
    real(dp), dimension(:,:), allocatable, public, save :: Epol

    ! multipole moments
    ! monopoles
    real(dp), dimension(:,:), allocatable, save :: Q0s
    ! dipoles
    real(dp), dimension(:,:), allocatable, save :: Q1s
    ! quadrupoles
    real(dp), dimension(:,:), allocatable, save :: Q2s
    ! octopoles
    real(dp), dimension(:,:), allocatable, save :: Q3s
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
    integer :: ndens = 0
    ! size of packed matrix in AO basis
    integer, save :: nnbas = 0
    ! number of nuclei in qm region
    integer, save :: qmnucs = 0
    ! nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zm
    ! nuclear coordinates
    real(dp), dimension(:,:), allocatable, save :: Rm

    ! frozen density fragment info
    ! ----------------------------

    ! number of frozen densities
    integer, public, save :: nfds = 0
    ! number of nuclei in current frozen density
    integer, public, save :: fdnucs = 0
    ! nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zfd
    ! nuclear coordinates
    real(dp), dimension(:,:), allocatable, save :: Rfd

! TODO:
! use allocate/deallocate where possible?
! insert quit if symmetry
! insert quits inside dalton if QM3, QMMM etc.
! consistent implementation where several densities are taken as input
! find better solution for electric field calculation from frozen densities
! hexadecapoles and higher order polarizabilities
! write list of publications which should be cited
! write output related to QMES
! avoid dimensions as input
! remove double zeroing and unecessary zeroing
! write to output
! memory checks and handle memory better
! nonlinear response properties
! magnetic properties
! AA and AU
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

    character(len=7), intent(inout) :: word
    integer, intent(in) :: luinp
    integer, intent(in) :: lupri

    character(len=7) :: option

    luout = lupri

    do
        read(luinp,'(a7)') option
        call upcase(option)

        ! do a Polarizable Embedding calculation
        if (trim(option(2:)) == 'PEQM') then
            peqm = .true.
        ! induced dipole - induced dipole damping
        else if (trim(option(2:)) == 'DAMP') then
            read(luinp,*) option
            backspace(luinp)
            if (option(1:1) /= '.' .and. option(1:1) /= '*'&
               &.and. option(1:1) /= '!') then
                read(luinp,*) damp
            end if
            pe_damp = .true.
        ! neglect dynamic response from environment
        else if (trim(option(2:)) == 'GSPOL') then
            pe_gspol = .true.
        else if (trim(option(2:)) == 'NOMB') then
            pe_nomb = .true.
        ! calculate intermolecular two-electron integrals
        else if (trim(option(2:)) == 'TWOINT') then
            read(luinp,*) fdnucs
            pe_twoint = .true.
        ! save density matrix
        else if (trim(option(2:)) == 'SAVDEN') then
            pe_savden = .true.
        ! get fock matrix for repulsion potential
        else if (trim(option(2:)) == 'REPULS') then
            pe_repuls = .true.
        ! electrostatics from frozen densities
        else if (trim(option(2:)) == 'FD') then
            ! number of frozen densities
            read(luinp,*) nfds
            pe_fd = .true.
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

    real(dp), dimension(:), intent(in), optional :: charges
    real(dp), dimension(:,:), intent(in), optional :: coords
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, k, s
    integer :: lupot, nlines
    integer, dimension(:), allocatable :: idxs
    real(dp) :: r
    real(dp), dimension(15) :: temp
    character(len=2) :: auoraa
    character(len=80) :: word
    logical :: lexist

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

    inquire(file='POTENTIAL.INP', exist=lexist)
    if (lexist) then
        call openfile('POTENTIAL.INP', lupot, 'old', 'formatted')
    else
        if (pe_savden) then
            return
        else if (pe_fd) then
            goto 101
        else
            stop('POTENTIAL.INP not found!')
        end if
    end if

    do
        read(lupot,*,end=100) word

        if (trim(word) == 'coordinates') then
            read(lupot,*) nsites
            read(lupot,*) auoraa
            allocate(elems(1,nsites), Zs(1,nsites), Rs(3,nsites))
            do i = 1, nsites
                read(lupot,*) elems(1,i), (Rs(j,i), j = 1, 3)
                Zs(1,i) = elem2charge(elems(1,i))
            end do
        else if (trim(word) == 'monopoles') then
            lmul(0) = .true.
            allocate(Q0s(1,nsites)); Q0s = 0.0d0
            temp = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, temp(1)
                Q0s(1,s) = temp(1)
            end do
        else if (trim(word) == 'dipoles') then
            lmul(1) = .true.
            allocate(Q1s(3,nsites)); Q1s = 0.0d0
            temp = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 3)
                Q1s(:,s) = temp(1:3)
            end do
        else if (trim(word) == 'quadrupoles') then
            lmul(2) = .true.
            allocate(Q2s(6,nsites)); Q2s = 0.0d0
            temp = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 6)
                Q2s(:,s) = temp(1:6)
            end do
        else if (trim(word) == 'octopoles') then
            lmul(3) = .true.
            allocate(Q3s(10,nsites)); Q3s = 0.0d0
            temp = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 10)
                Q3s(:,s) = temp(1:10)
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
            allocate(exlists(lexlst,nsites))
            do i = 1, nsites
                read(lupot,*) (exlists(j,i), j = 1, lexlst)
            end do
        else if (word(1:1) == '!' .or. word(1:1) == '#') then
            cycle
        end if
    end do

100 continue

    close(lupot)

    ! initialize nloop
    nloop = nsites

    ! if coordinates are in AA then convert to AU
    if (auoraa == 'AA') then
        Rs = Rs * aa2au
    end if

    ! default exclusion list (everything polarizes everything)
    if (.not. allocated(exlists)) then
        lexlst = 1
        allocate(exlists(lexlst,nsites))
        do i = 1, nsites
            exlists(1,i) = i
        end do
    end if

!   ! Handle quantum-classical boundary
    mindist = 2.0d0
    allocate(idxs(nsites))
!   r = 1.0d10
!   do i = 1, qmnucs
!       do j = 1, nsites
!           if (nrm2(Rm(:,i) - Rs(:,j)) < r) then
!               r = nrm2(Rm(:,i) - Rs(:,j))
!               idx = j
!           end if
!       end do
!   end do

    ! number of polarizabilities different from zero
    if (lpol(0) .or. lpol(1)) then
        allocate(zeroalphas(nsites))
        do i = 1, nsites
            if (abs(maxval(alphas(:,i))) >= zero) then
                zeroalphas(i) = .false.
                npols = npols + 1
            else
                zeroalphas(i) = .true.
            end if
        end do
    end if

    ! write to Dalton output file
101 write(luout,'(//2x,a)') 'Polarizable Embedding potential'
    write(luout,'(2x,a)')   '-------------------------------'
    write(luout,'(/4x,a,i6)') 'Number of classical sites: ', nsites
    if (lmul(3)) then
        write(luout,'(4x,a)') 'Multipole moments upto 3rd order.'
    else if (lmul(2)) then
        write(luout,'(4x,a)') 'Multipole moments upto 2nd order.'
    else if (lmul(1)) then
        write(luout,'(4x,a)') 'Multipole moments upto 1st order.'
    else if (lmul(0)) then
        write(luout,'(4x,a)') 'Multipole moments upto 0th order.'
    end if
    if (lpol(0) .and. lpol(1)) then
        write(luout,'(4x,a)') 'Isotropic and anisotropic'//&
                              ' dipole-dipole polarizabilities.'
    else if (lpol(0) .and. .not. lpol(1)) then
        write(luout,'(4x,a)') 'Isotropic'//&
                              ' dipole-dipole polarizabilities.'
    else if (.not. lpol(0) .and. lpol(1)) then
        write(luout,'(4x,a)') 'Anisotropic'//&
                              ' dipole-dipole polarizabilities.'
    end if
    if (pe_fd) then
        write(luout,'(4x,a,i4)') 'Number of frozen densities:', nfds
    end if

end subroutine pe_read_potential

!------------------------------------------------------------------------------

subroutine pe_master(runtype, denmats, fckmats, nmats, Epe, work)

    character(*), intent(in) :: runtype
    integer, intent(in) :: nmats
    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out), optional :: fckmats
    real(dp), dimension(:), intent(out), optional :: Epe
    real(dp), dimension(:), intent(inout) :: work

#ifdef VAR_MPI
    integer :: myid, ncores, ierr

    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, ncores, ierr)
#endif

    ! determine what to calculate and do consistency check
    if (runtype == 'fock') then
        fock = .true.
        energy = .false.
        response = .false.
        if (.not. present(fckmats)) then
            stop('Output matrices are missing from input!')
        else if (.not. present(Epe)) then
            stop('The energy variable is missing from input!')
        end if
    else if (runtype == 'energy') then
        energy = .true.
        fock = .false.
        response = .false.
    else if (runtype == 'response') then
        if (pe_gspol) return
        if (.not. lpol(0) .and. .not. lpol(1)) return
        response = .true.
        fock = .false.
        energy = .false.
        if (.not. present(fckmats)) then
            stop('Output matrices are missing from input!')
        end if
    else
        stop('Could not determine calculation type.')
    end if

    ndens = nmats
    nnbas = size(denmats) / ndens

#ifdef VAR_MPI
    if (myid == 0 .and. ncores > 1) then
        call mpi_bcast(44, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        if (fock) then
            call mpi_bcast(1, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        else if (energy) then
            call mpi_bcast(2, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        else if (response) then
            call mpi_bcast(3, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        end if

        call mpi_bcast(nnbas, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        call mpi_bcast(ndens, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        call mpi_bcast(denmats, nnbas*ndens, MPI_REAL8, myid, MPI_COMM_WORLD, ierr)

        if (.not. initialized) then
            call pe_sync(work)
        end if
    end if
#endif

    if (fock) then
        call pe_fock(denmats, fckmats, Epe, work)
    else if (energy) then
        call pe_fock(denmats, work=work)
    else if (response) then
        call pe_response(denmats, fckmats, work)
!        call pe_polarization(denmats, fckmats, work)
    end if

#ifdef VAR_MPI
    if (myid == 0 .and. ncores > 1) then
        if (fock .or. response) then
            call mpi_reduce(MPI_IN_PLACE, fckmats, ndens*nnbas, MPI_REAL8,&
                           &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        end if
    end if
#endif

end subroutine pe_master

!------------------------------------------------------------------------------

#ifdef VAR_MPI
subroutine pe_mpi(work, runtype)

    real(dp), dimension(:), intent(inout) :: work
    integer :: runtype

    integer :: i
    integer :: nwrk
    integer :: myid, ncores, ierr

    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, ncores, ierr)

    nwrk = size(work)

    if (runtype == 1) then
        fock = .true.
        energy = .false.
        response = .false.
    else if (runtype == 2) then
        energy = .true.
        fock = .false.
        response = .false.
    else if (runtype == 3) then
        response = .true.
        fock = .false.
        energy = .false.
    end if

    call mpi_bcast(nnbas, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ndens, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(work(1), nnbas*ndens, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    if (.not. initialized) then
        call pe_sync(work(nnbas*ndens+1:nwrk))
    end if

    if (fock) then
        call pe_fock(work(1:ndens*nnbas), work(ndens*nnbas+1:2*ndens*nnbas),&
                    &work(2*ndens*nnbas+1:2*ndens*nnbas+ndens),&
                    &work(2*ndens*nnbas+ndens+1:nwrk))
    else if (energy) then
        call pe_fock(work(1:ndens*nnbas), work=work(ndens*nnbas+1:nwrk))
    else if (response) then
        call pe_response(work(1:ndens*nnbas),&
                        &work(ndens*nnbas+1:2*ndens*nnbas),&
                        &work(2*ndens*nnbas+1:nwrk))
    end if

    if (fock .or. response) then
        call mpi_reduce(work(ndens*nnbas+1), 0, ndens*nnbas, MPI_REAL8,&
                       &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if

end subroutine pe_mpi

!------------------------------------------------------------------------------

subroutine pe_sync(work)

    real(dp), dimension(:) :: work

    integer :: i
    integer :: ndist, nrest
    integer :: myid, ncores, ierr

    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, ncores, ierr)

    if (myid == 0) then
        allocate(ndists(0:ncores-1), displs(0:ncores-1))
        ndist = nsites / ncores
        ndists = ndist
        if (ncores * ndist < nsites) then
            nrest = nsites - ncores * ndist
            do i = 0, nrest-1
                ndists(i) = ndists(i) + 1
            end do
        end if
        nloop = ndists(0)
        call mpi_scatter(ndists, 1, MPI_INTEGER,&
                        &MPI_IN_PLACE, 1, MPI_INTEGER,&
                        &0, MPI_COMM_WORLD, ierr)
    else
        call mpi_scatter(0, 0, MPI_INTEGER,&
                        &nloop, 1, MPI_INTEGER,&
                        &0, MPI_COMM_WORLD, ierr)
    end if

    if (myid == 0) then
        displs(0) = 0
        do i = 1, ncores-1
            displs(i) = displs(i-1) + 3 * ndists(i-1)
        end do
        call mpi_scatterv(Rs, 3*ndists, displs, MPI_REAL8,&
                         &MPI_IN_PLACE, 0, MPI_REAL8,&
                         &0, MPI_COMM_WORLD, ierr)
    else
        allocate(Rs(3,nloop))
        call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                         &Rs, 3*nloop, MPI_REAL8,&
                         &0, MPI_COMM_WORLD, ierr)
    end if

    call mpi_bcast(lpol, 2, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    if (lpol(0) .or. lpol(1)) then
        call mpi_bcast(npols, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (myid == 0) then
            displs(0) = 0
            do i = 1, ncores-1
                displs(i) = displs(i-1) + ndists(i-1)
            end do
            call mpi_scatterv(zeroalphas, ndists, displs, MPI_LOGICAL,&
                             &MPI_IN_PLACE, 0, MPI_LOGICAL,&
                             &0, MPI_COMM_WORLD, ierr)
        else
            allocate(zeroalphas(nloop))
            call mpi_scatterv(0, 0, 0, MPI_LOGICAL,&
                             &zeroalphas, nloop, MPI_LOGICAL,&
                             &0, MPI_COMM_WORLD, ierr)
        end if
    end if

    call mpi_bcast(qmnucs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if (myid /= 0) allocate(Zm(1,qmnucs))
    call mpi_bcast(Zm, qmnucs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    if (myid /= 0) allocate(Rm(3,qmnucs))
    call mpi_bcast(Rm, 3*qmnucs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    call mpi_bcast(lmul, 4, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    if (lmul(0)) then
        if (myid == 0) then
            displs(0) = 0
            do i = 1, ncores-1
                displs(i) = displs(i-1) + ndists(i-1)
            end do
            call mpi_scatterv(Q0s, ndists, displs, MPI_REAL8,&
                             &MPI_IN_PLACE, 0, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        else
            allocate(Q0s(1,nloop))
            call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                             &Q0s, nloop, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        end if
    end if

    if (lmul(1)) then
        if (myid == 0) then
            displs(0) = 0
            do i = 1, ncores-1
                displs(i) = displs(i-1) + 3 * ndists(i-1)
            end do
            call mpi_scatterv(Q1s, 3*ndists, displs, MPI_REAL8,&
                             &MPI_IN_PLACE, 0, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        else
            allocate(Q1s(3,nloop))
            call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                             &Q1s, 3*nloop, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        end if
    end if

    if (lmul(2)) then
        if (myid == 0) then
            displs(0) = 0
            do i = 1, ncores-1
                displs(i) = displs(i-1) + 6 * ndists(i-1)
            end do
            call mpi_scatterv(Q2s, 6*ndists, displs, MPI_REAL8,&
                             &MPI_IN_PLACE, 0, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        else
            allocate(Q2s(6,nloop))
            call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                             &Q2s, 6*nloop, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        end if
    end if

    if (lmul(3)) then
        if (myid == 0) then
            displs(0) = 0
            do i = 1, ncores-1
                displs(i) = displs(i-1) + 10 * ndists(i-1)
            end do
            call mpi_scatterv(Q3s, 10*ndists, displs, MPI_REAL8,&
                             &MPI_IN_PLACE, 0, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        else
            allocate(Q3s(10,nloop))
            call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                             &Q3s, 10*nloop, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        end if
    end if

    initialized = .true.

end subroutine pe_sync
#endif

!------------------------------------------------------------------------------

subroutine pe_fock(denmats, fckmats, Epe, work)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out), optional :: fckmats
    real(dp), dimension(:), intent(out), optional :: Epe
    real(dp), dimension(:), intent(inout) :: work

    integer :: i
    logical :: es = .false.
    logical :: pol = .false.

    do i = 0, 3
        if (lmul(i)) es = .true.
    end do
    if (pe_fd) es = .true.
    do i = 0, 1
        if (lpol(i)) pol = .true.
    end do

    if (allocated(Ees)) deallocate(Ees)
    if (allocated(Epol)) deallocate(Epol)
    allocate(Ees(0:4,ndens), Epol(4,ndens))
    Ees = 0.0d0; Epol = 0.0d0

    if (fock) fckmats = 0.0d0

    if (fock) then
        if (es) call pe_electrostatic(denmats, fckmats, work)
        if (pol) call pe_polarization(denmats, fckmats, work)
    else if (energy) then
        if (es) call pe_electrostatic(denmats, work=work)
        if (pol) call pe_polarization(denmats, work=work)
    end if

    if (fock) then
        Epe = 0.0d0
        do i = 1, ndens
            Epe(i) = sum(Ees(:,i)) + sum(Epol(:,i))
        end do
    end if

end subroutine pe_fock

!------------------------------------------------------------------------------

subroutine pe_response(denmats, fckmats, work)

    external :: Tk_integrals

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out) :: fckmats
    real(dp), dimension(:), intent(inout) :: work

    logical :: skip
    integer :: i, j, l, m, n, o
    integer, parameter :: k = 1
    real(dp), dimension(3*npols,ndens) :: Fel, Mu
    real(dp), dimension(nnbas,3) :: Fel_ints

#ifdef VAR_MPI
    integer :: ndist, nrest
    integer :: myid, ncores, ierr

    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, ncores, ierr)
#endif

    Fel = 0.0d0; fckmats = 0.0d0

    l = 0
    do i = 1, nloop
        if (zeroalphas(i)) cycle
        call Tk_integrals(Fel_ints, 3*nnbas, k, Rs(:,i), work, size(work))
        do j = 1, 3
            do m = 1, ndens
                n = (m - 1) * nnbas + 1
                o = m * nnbas
                Fel(l+j,m) = dot(denmats(n:o), Fel_ints(:,j))
            end do
        end do
        l = l + 3
    end do

#ifdef VAR_MPI
    if (myid == 0 .and. ncores > 1) then
        ndist = npols / ncores
        ndists = ndist
        if (ncores * ndist < npols) then
            nrest = npols - ncores * ndist
            do i = 0, nrest-1
                ndists(i) = ndists(i) + 1
            end do
        end if
        call mpi_scatter(ndists, 1, MPI_INTEGER,&
                        &MPI_IN_PLACE, 1, MPI_INTEGER,&
                        &0, MPI_COMM_WORLD, ierr)
    else if (myid /= 0) then
            call mpi_scatter(0, 0, MPI_INTEGER,&
                            &ndist, 1, MPI_INTEGER,&
                            &0, MPI_COMM_WORLD, ierr)
    end if

    if (myid == 0 .and. ncores > 1) then
        displs(0) = 0
        do i = 1, ncores-1
            displs(i) = displs(i-1) + 3 * ndists(i-1) 
        end do
        do i = 1, ndens
            call mpi_gatherv(MPI_IN_PLACE, 0, MPI_REAL8,&
                            &Fel(:,i), 3*ndists, displs, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        end do
    else if (myid /= 0) then
        do i = 1, ndens
            call mpi_gatherv(Fel(:,i), 3*ndist, MPI_REAL8,&
                            &0, 0, 0, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        end do
    end if

    if (myid == 0) then
#endif
        call induced_dipoles(Mu, Fel)
#ifdef VAR_MPI
    end if

    if (myid == 0 .and. ncores > 1) then
        displs(0) = 0
        do i = 1, ncores-1
            displs(i) = displs(i-1) + 3 * ndists(i-1)
        end do
        do i = 1, ndens
            call mpi_scatterv(Mu(:,i), 3*ndists, displs, MPI_REAL8,&
                             &MPI_IN_PLACE, 0, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        end do
    else if (myid /= 0) then
        do i = 1, ndens
            call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                             &Mu(:,i), 3*ndist, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        end do
    end if
#endif

    l = 0
    do i = 1, nloop
        if (zeroalphas(i)) cycle
        call Tk_integrals(Fel_ints, 3*nnbas, k, Rs(:,i), work, size(work))
        do j = 1, 3
            do m = 1, ndens
                n = (m - 1) * nnbas + 1
                o = m * nnbas
                fckmats(n:o) = fckmats(n:o) - Mu(l+j,m) * Fel_ints(:,j)
            end do
        end do
        l = l + 3
    end do

end subroutine pe_response

!------------------------------------------------------------------------------

subroutine pe_electrostatic(denmats, fckmats, work)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, k
    logical :: lexist
    integer :: lutemp
    real(dp) :: Enuc, Esave
    real(dp), dimension(ndens) :: Eel

#ifdef VAR_MPI
    real(dp), dimension(:), allocatable :: tmpfcks
    integer :: myid, ncores, ierr

    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, ncores, ierr)

    if (myid == 0) then
#endif
        inquire(file='pe_electrostatics.bin', exist=lexist)
#ifdef VAR_MPI
    end if

    if (ncores > 1) then
        call mpi_bcast(lexist, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    end if
#endif

    if (lexist .and. fock) then
#ifdef VAR_MPI
        if (myid == 0) then
#endif
            call openfile('pe_electrostatics.bin', lutemp, 'old', 'unformatted')
            rewind(lutemp)
            read(lutemp) Esave, fckmats
            close(lutemp)
            do i = 1, ndens
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                Ees(0,i) = Ees(0,i) + dot(denmats(j:k), fckmats(j:k)) + Esave
            end do
#ifdef VAR_MPI
        end if
#endif
    else
        Esave = 0.0d0
        if (lmul(0)) then
            if (fock) then
                call es_monopoles(denmats, Eel, Enuc, fckmats, work)
            else if (energy) then
                call es_monopoles(denmats, Eel, Enuc, work=work)
            end if
            do i = 1, ndens
                Ees(0,i) = Ees(0,i) + Eel(i) + Enuc
            end do
            Esave = Esave + Enuc
        end if
        if (lmul(1)) then
            if (fock) then
                call es_dipoles(denmats, Eel, Enuc, fckmats, work)
            else if (energy) then
                call es_dipoles(denmats, Eel, Enuc, work=work)
            end if
            do i = 1, ndens
                Ees(1,i) = Ees(1,i) + Eel(i) + Enuc
            end do
            Esave = Esave + Enuc
        end if
        if (lmul(2)) then
            if (fock) then
                call es_quadrupoles(denmats, Eel, Enuc, fckmats, work)
            else if (energy) then
                call es_quadrupoles(denmats, Eel, Enuc, work=work)
            end if
            do i = 1, ndens
                Ees(2,i) = Ees(2,i) + Eel(i) + Enuc
            end do
            Esave = Esave + Enuc
        end if
        if (lmul(3)) then
            if (fock) then
                call es_octopoles(denmats, Eel, Enuc, fckmats, work)
            else if (energy) then
                call es_octopoles(denmats, Eel, Enuc, work=work)
            end if
            do i = 1, ndens
                Ees(3,i) = Ees(3,i) + Eel(i) + Enuc
            end do
            Esave = Esave + Enuc
        end if
#ifdef VAR_MPI
        if (myid == 0) then
#endif
            if (pe_fd) then
                if (fock) then
                    call es_frozen_densities(denmats, Eel, Enuc, fckmats, work)
                else if (energy) then
                    call es_frozen_densities(denmats, Eel, Enuc, work=work)
                end if
                do i = 1, ndens
                    Ees(4,i) = Ees(4,i) + Eel(i) + Enuc
                end do
                Esave = Esave + Enuc
            end if
#ifdef VAR_MPI
        end if
#endif
        if (fock) then
#ifdef VAR_MPI
            if (myid == 0 .and. ncores > 1) then
                allocate(tmpfcks(ndens*nnbas))
                tmpfcks = fckmats
                call mpi_reduce(MPI_IN_PLACE, fckmats, ndens*nnbas, MPI_REAL8,&
                               &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(MPI_IN_PLACE, Esave, 1, MPI_REAL8, MPI_SUM,&
                               &0, MPI_COMM_WORLD, ierr)
            else if (myid /= 0) then
                call mpi_reduce(fckmats, 0, ndens*nnbas, MPI_REAL8, MPI_SUM,&
                               &0, MPI_COMM_WORLD, ierr)
                call mpi_reduce(Esave, 0, 1, MPI_REAL8, MPI_SUM,&
                               &0, MPI_COMM_WORLD, ierr)
            end if
            if (myid == 0) then
#endif
                call openfile('pe_electrostatics.bin', lutemp, 'new', 'unformatted')
                rewind(lutemp)
                write(lutemp) Esave, fckmats
                close(lutemp)
#ifdef VAR_MPI
            end if

            if (myid == 0 .and. ncores > 1) then
                fckmats = tmpfcks
                deallocate(tmpfcks)
            end if
#endif
        end if
#ifdef VAR_MPI
        if (myid == 0 .and. ncores > 1) then
            call mpi_reduce(MPI_IN_PLACE, Ees, 5*ndens, MPI_REAL8, MPI_SUM,&
                           &0, MPI_COMM_WORLD, ierr)
        else if (myid /= 0) then
            call mpi_reduce(Ees, 0, 5*ndens, MPI_REAL8, MPI_SUM,&
                           &0, MPI_COMM_WORLD, ierr)
        end if
#endif
    end if

end subroutine pe_electrostatic

!------------------------------------------------------------------------------

subroutine es_frozen_densities(denmats, Eel, Enuc, fckmats, work)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out) :: Eel
    real(dp), intent(out) :: Enuc
    real(dp), dimension(:), intent(inout), optional :: fckmats
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l, m, n, o
    integer, parameter :: k = 0
    integer :: lufck, lexist, lutemp
    real(dp) :: Ene, Enn
    real(dp), dimension(ndens) :: Een, Eee
    real(dp), dimension(1) :: Tfm
    real(dp), dimension(3) :: Rfm
    real(dp), dimension(3*npols) :: temp
    real(dp), dimension(nnbas) :: fd_fock
    real(dp), dimension(nnbas,1) :: Zfd_ints
    character(len=99) :: ci
    character(len=99) :: filename

    Eel = 0.0d0; Enuc = 0.0d0

    do i = 1, nfds
        Eee = 0.0d0; Een = 0.0d0; Ene = 0.0d0; Enn = 0.0d0
        write(ci,*) i
        ci = adjustl(ci)
        filename = 'pe_fock_'//trim(ci)//'.bin'
        call openfile(trim(filename), lufck, 'old', 'unformatted')
        rewind(lufck)
        read(lufck) temp
        read(lufck) Ene
        read(lufck) fd_fock
        read(lufck) fdnucs
        allocate(Rfd(3,fdnucs), Zfd(1,fdnucs))
        read(lufck) Rfd, Zfd
        close(lufck)

        do j = 1, ndens
            l = (j - 1) * nnbas + 1
            m = j * nnbas
            if (fock) fckmats(l:m) = fckmats(l:m) + fd_fock
            Eee(j) = dot(denmats(l:m), fd_fock)
        end do

        do j = 1, fdnucs
            call Qk_integrals(Zfd_ints, k, Rfd(:,j), Zfd(:,j), work)
            do m = 1, ndens
                n = (m - 1) * nnbas + 1
                o = m * nnbas
                Een(m) = Een(m) + dot(denmats(n:o), Zfd_ints(:,1))
                if (fock) fckmats(n:o) = fckmats(n:o) + Zfd_ints(:,1)
            end do
            do l = 1, qmnucs
                Rfm = Rm(:,l) - Rfd(:,j)
                call Tk_tensor(Tfm, k, Rfm)
                Enn = Enn + Zm(1,l) * Zfd(1,j) * Tfm(1)
            end do
        end do

        deallocate(Rfd, Zfd)

        Enuc = Enuc + Ene + Enn
        do j = 1, ndens
            Eel(j) = Eel(j) + Een(j) + Eee(j)
        end do
    end do

end subroutine es_frozen_densities

!------------------------------------------------------------------------------

subroutine es_monopoles(denmats, Eel, Enuc, fckmats, work)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out) :: Eel
    real(dp), intent(out) :: Enuc
    real(dp), dimension(:), intent(inout), optional :: fckmats
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l, m
    integer, parameter :: k = 0
    real(dp), dimension(3) :: Rsm
    real(dp), dimension(1) :: Tsm
    real(dp), dimension(nnbas,1) :: Q0_ints

    Eel = 0.0d0; Enuc = 0.0d0

    do i = 1, nloop
        if (abs(Q0s(1,i)) < zero) cycle

        ! nuclei - monopole interaction
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call Tk_tensor(Tsm, k, Rsm)
            Enuc = Enuc + Q0s(1,i) * Zm(1,j) * Tsm(1)
        end do

        ! electron - monopole interaction
        call Qk_integrals(Q0_ints, k, Rs(:,i), Q0s(:,i), work)
        do j = 1, ndens
            l = (j - 1) * nnbas + 1
            m = j * nnbas
            Eel(j) = Eel(j) + dot(denmats(l:m), Q0_ints(:,1))
            if (fock) fckmats(l:m) = fckmats(l:m) + Q0_ints(:,1)
        end do
    end do

end subroutine es_monopoles

!------------------------------------------------------------------------------

subroutine es_dipoles(denmats, Eel, Enuc, fckmats, work)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out) :: Eel
    real(dp), intent(out) :: Enuc
    real(dp), dimension(:), intent(inout), optional :: fckmats
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l, m, n
    integer, parameter :: k = 1
    real(dp), dimension(3) :: Rsm, Tsm
    real(dp), dimension(nnbas,3) :: Q1_ints

    Eel = 0.0d0; Enuc = 0.0d0

    do i = 1, nloop
        if (abs(maxval(Q1s(:,i))) < zero) cycle

        ! nuclei - dipole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call Tk_tensor(Tsm, k, Rsm)
            do l = 1, 3
                Enuc = Enuc - Zm(1,j) * Q1s(l,i) * Tsm(l)
            end do
        end do

        ! electron - dipole interaction
        call Qk_integrals(Q1_ints, k, Rs(:,i), Q1s(:,i), work)
        do j = 1, 3
            do l = 1, ndens
                m = (l - 1) * nnbas + 1
                n = l * nnbas
                Eel(l) = Eel(l) + dot(denmats(m:n), Q1_ints(:,j))
                if (fock) fckmats(m:n) = fckmats(m:n) + Q1_ints(:,j)
            end do
        end do
    end do

end subroutine es_dipoles

!------------------------------------------------------------------------------

subroutine es_quadrupoles(denmats, Eel, Enuc, fckmats, work)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out) :: Eel
    real(dp), intent(out) :: Enuc
    real(dp), dimension(:), intent(inout), optional :: fckmats
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l, m, n
    integer, parameter :: k = 2
    real(dp), dimension(3) :: Rsm
    real(dp), dimension(6) :: Tsm, factors
    real(dp), dimension(nnbas,6) :: Q2_ints

    Eel = 0.0d0; Enuc = 0.0d0

    do i = 1, nloop
        if (abs(maxval(Q2s(:,i))) < zero) cycle

        ! nuclei - quadrupole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call Tk_tensor(Tsm, k, Rsm)
            call symmetry_factors(factors, k)
            do l = 1, 6
                Enuc = Enuc + 0.5d0 * factors(l) * Zm(1,j) * Q2s(l,i) * Tsm(l)
            end do
        end do

        ! electron - quadrupole interaction energy
        call Qk_integrals(Q2_ints, k, Rs(:,i), Q2s(:,i), work)
        do j = 1, 6
            do l = 1, ndens
                m = (l - 1) * nnbas + 1
                n = l * nnbas
                Eel(l) = Eel(l) + dot(denmats(m:n), Q2_ints(:,j))
                if (fock) fckmats(m:n) = fckmats(m:n) + Q2_ints(:,j)
            end do
        end do
    end do

end subroutine es_quadrupoles

!------------------------------------------------------------------------------

subroutine es_octopoles(denmats, Eel, Enuc, fckmats, work)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out) :: Eel
    real(dp), intent(out) :: Enuc
    real(dp), dimension(:), intent(inout), optional :: fckmats
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l, m, n
    integer, parameter :: k = 3
    real(dp), dimension(3) :: Rsm
    real(dp), dimension(10) :: Tsm, factors
    real(dp), dimension(nnbas,10) :: Q3_ints

    Eel = 0.0d0; Enuc = 0.0d0

    do i = 1, nloop
        if (abs(maxval(Q3s(:,i))) < zero) cycle

        ! nuclei - octopole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,i)
            call Tk_tensor(Tsm, k, Rsm)
            call symmetry_factors(factors, k)
            do l = 1, 10
                Enuc = Enuc - factors(l) * Zm(1,j) * Q3s(l,i) * Tsm(l) / 6.0d0
            end do
        end do

        ! electron - octopole interaction energy
        call Qk_integrals(Q3_ints, k, Rs(:,i), Q3s(:,i), work)
        do j = 1, 10
            do l = 1, ndens
                m = (l - 1) * nnbas + 1
                n = l * nnbas
                Eel(l) = Eel(l) + dot(denmats(m:n), Q3_ints(:,j))
                if (fock) fckmats(m:n) = fckmats(m:n) + Q3_ints(:,j)
            end do
        end do
    end do

end subroutine es_octopoles

!------------------------------------------------------------------------------

subroutine pe_polarization(denmats, fckmats, work)

    external :: Tk_integrals

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l, m, n, o
    integer, parameter :: k = 1
    logical :: skip
    real(dp), dimension(3*npols) :: Fnuc, Fmul, Ffd
    real(dp), dimension(3*npols,ndens) :: Mu, Fel, Ftot
    real(dp), dimension(nnbas,3) :: Fel_ints

#ifdef VAR_MPI
    integer :: ndist, nrest
    integer :: myid, ncores, ierr

    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, ncores, ierr)
#endif

    Fel = 0.0d0

    l = 0
    do i = 1, nloop
        if (zeroalphas(i)) cycle
        if (pe_savden) then
! TODO: find better threshold or better solution
            skip = .false.
            do j = 1, qmnucs
                if (nrm2(Rs(:,i) - Rm(:,j)) <= 1.0d0) skip = .true.
            end do
            if (skip) cycle
        end if
        call Tk_integrals(Fel_ints, 3*nnbas, k, Rs(:,i), work, size(work))
        do j = 1, 3
            do m = 1, ndens
                n = (m - 1) * nnbas + 1
                o = m * nnbas
                Fel(l+j,m) = dot(denmats(n:o), Fel_ints(:,j))
            end do
        end do
        l = l + 3
    end do

#ifdef VAR_MPI
    if (myid == 0 .and. ncores > 1) then
        ndist = npols / ncores
        ndists = ndist
        if (ncores * ndist < npols) then
            nrest = npols - ncores * ndist
            do i = 0, nrest-1
                ndists(i) = ndists(i) + 1
            end do
        end if
        call mpi_scatter(ndists, 1, MPI_INTEGER,&
                        &MPI_IN_PLACE, 1, MPI_INTEGER,&
                        &0, MPI_COMM_WORLD, ierr)
    else if (myid /= 0) then
        call mpi_scatter(0, 0, MPI_INTEGER,&
                        &ndist, 1, MPI_INTEGER,&
                        &0, MPI_COMM_WORLD, ierr)
    end if

    if (myid == 0 .and. ncores > 1) then
        displs(0) = 0
        do i = 1, ncores-1
            displs(i) = displs(i-1) + 3 * ndists(i-1)
        end do
        do i = 1, ndens
            call mpi_gatherv(MPI_IN_PLACE, 0, MPI_REAL8,&
                            &Fel(:,i), 3*ndists, displs, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        end do
    else if (myid /= 0) then
        do i = 1, ndens
            call mpi_gatherv(Fel(:,i), 3*ndist, MPI_REAL8,&
                            &0, 0, 0, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        end do
    end if

    if (myid == 0) then
#endif
    call nuclear_fields(Fnuc)
    call multipole_fields(Fmul)
    if (pe_fd) then
        call frozen_density_field(Ffd)
    else
        Ffd = 0.0d0
    end if

    do i = 1, ndens
        Ftot(:,i) = Fel(:,i) + Fnuc + Fmul + Ffd
    end do

    call induced_dipoles(Mu, Ftot)

    do i = 1, ndens
        Epol(1,i) = - 0.5d0 * dot(Mu(:,i), Fel(:,i))
        Epol(2,i) = - 0.5d0 * dot(Mu(:,i), Fnuc)
        Epol(3,i) = - 0.5d0 * dot(Mu(:,i), Fmul)
        if (pe_fd) Epol(4,i) = - 0.5d0 * dot(Mu(:,i), Ffd)
    end do

#ifdef VAR_MPI
    end if

    if (myid == 0 .and. ncores > 1) then
        displs(0) = 0
        do i = 1, ncores-1
            displs(i) = displs(i-1) + 3 * ndists(i-1)
        end do
        do i = 1, ndens
            call mpi_scatterv(Mu(:,i), 3*ndists, displs, MPI_REAL8,&
                             &MPI_IN_PLACE, 0, MPI_REAL8,&
                             &myid, MPI_COMM_WORLD, ierr)
        end do
    else if (myid /= 0) then
        do i = 1, ndens
            call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                             &Mu(:,i), 3*ndist, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        end do
    end if
#endif

    if (fock) then
        l = 0
        do i = 1, nloop
            if (zeroalphas(i)) cycle
            call Tk_integrals(Fel_ints, 3*nnbas, k, Rs(:,i), work, size(work))
            do j = 1, 3
                do m = 1, ndens
                    n = (m - 1) * nnbas + 1
                    o = m * nnbas
                    fckmats(n:o) = fckmats(n:o) - Mu(j+l,m) * Fel_ints(:,j)
                end do
            end do
            l = l + 3
        end do
    end if

end subroutine pe_polarization

!------------------------------------------------------------------------------

subroutine induced_dipoles(Mu, F)

    real(dp), dimension(:,:), intent(out) :: Mu
    real(dp), dimension(:,:), intent(in) :: F

    integer :: i, j, l
    real(dp), dimension(:), allocatable :: B

    allocate(B(3*npols*(3*npols+1)/2))

    call response_matrix(B)

    do i = 1, ndens
        call spmv(B, F(:,i), Mu(:,i), 'L')
    end do

    do i = 1, ndens
        l = 1
        do j = 1, npols
            if (nrm2(Mu(l:l+2,i)) > 1.0d0) then
                print *, 'large induced dipole encountered'
                print *, Mu(l:l+2,i)
            end if
            l = l + 3
        end do
    end do

    deallocate(B)

end subroutine induced_dipoles

!------------------------------------------------------------------------------

subroutine electron_fields(Fel, denmats, work)

    external :: Tk_integrals

    real(dp), dimension(:,:), intent(out) :: Fel
    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(inout) :: work

    logical :: skip
    integer :: i, j, l, m, n, o
    integer, parameter :: k = 1
    real(dp), dimension(nnbas,3) :: Fel_ints

    Fel = 0.0d0

    l = 0
    do i = 1, nloop
        if (zeroalphas(i)) cycle
        if (pe_savden) then
! TODO: find better threshold or better solution
            skip = .false.
            do j = 1, qmnucs
                if (nrm2(Rs(:,i) - Rm(:,j)) <= 1.0d0) skip = .true.
            end do
            if (skip) cycle
        end if
        call Tk_integrals(Fel_ints, 3*nnbas, k, Rs(:,i), work, size(work))
        do j = 1, 3
            do m = 1, ndens
                n = (m - 1) * nnbas + 1
                o = m * nnbas
                Fel(l+j,m) = dot(denmats(n:o), Fel_ints(:,j))
            end do
        end do
        l = l + 3
    end do

end subroutine electron_fields

!------------------------------------------------------------------------------

subroutine nuclear_fields(Fnuc)

    real(dp), dimension(:), intent(out) :: Fnuc

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
            if (zeroalphas(i)) cycle
            if (pe_savden) then
! TODO: finde better threshold or better solution
                skip = .false.
                do j = 1, qmnucs
                    if (nrm2(Rs(:,i) - Rm(:,j)) <= 1.0d0) skip = .true.
                end do
                if (skip) cycle
            end if
            do j = 1, qmnucs
                Rms = Rs(:,i) - Rm(:,j)
                call Tk_tensor(Tms, k, Rms)
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

end subroutine nuclear_fields

!------------------------------------------------------------------------------

subroutine frozen_density_field(Ffd)

    real(dp), dimension(:), intent(out) :: Ffd

    integer :: i
    integer :: lutemp
    character(len=99) :: ci
    character(len=80) :: filename
    real(dp), dimension(3*npols) :: Ftmp

    Ffd = 0.0d0

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

end subroutine frozen_density_field

!------------------------------------------------------------------------------

subroutine multipole_fields(Fmul)

    real(dp), dimension(:), intent(out) :: Fmul

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
                if (i == j) cycle
                exclude = .false.
                do k = 1, lexlst
                    if (exlists(k,i) == exlists(1,j)) then
                        exclude = .true.
                        exit
                    end if
                end do
                if (exclude) cycle
! TODO: cutoff???
!                Rij = Rs(:,j) - Rs(:,j)
!
                ! get electric field at i due to monopole at j
                if (lmul(0)) then
                    if (abs(maxval(Q0s(:,j))) >= zero) then
                        call monopole_field(Fs, Rs(:,i), Rs(:,j), Q0s(:,j))
                        Fmul(l:l+2) = Fmul(l:l+2) + Fs
                    end if
                end if
                ! get electric field at i due to dipole at j
                if (lmul(1)) then
                    if (abs(maxval(Q1s(:,j))) >= zero) then
                        call dipole_field(Fs, Rs(:,i), Rs(:,j), Q1s(:,j))
                        Fmul(l:l+2) = Fmul(l:l+2) + Fs
                    end if
                end if
                ! get electric field at i due to quadrupole at j
                if (lmul(2)) then
                    if (abs(maxval(Q2s(:,j))) >= zero) then
                        call quadrupole_field(Fs, Rs(:,i), Rs(:,j), Q2s(:,j))
                        Fmul(l:l+2) = Fmul(l:l+2) + Fs
                    end if
                end if
                ! get electric field at i due to octopole at j
                if (lmul(3)) then
                    if (abs(maxval(Q3s(:,j))) >= zero) then
                        call octopole_field(Fs, Rs(:,i), Rs(:,j), Q3s(:,j))
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

end subroutine multipole_fields

!------------------------------------------------------------------------------

subroutine monopole_field(Fi, Ri, Rj, Q0j)

    real(dp), dimension(3), intent(out) :: Fi
    real(dp), dimension(3), intent(in) :: Ri, Rj
    real(dp), dimension(1), intent(in) :: Q0j

    integer :: a
    integer, parameter :: k = 1
    real(dp), dimension(3) :: Rji, Tji

    Rji = Ri - Rj

    call Tk_tensor(Tji, k, Rji)

    Fi = 0.0d0

    do a = 1, 3
        Fi(a) = Fi(a) - Tji(a) *  Q0j(1)
    end do

end subroutine monopole_field

!------------------------------------------------------------------------------

subroutine dipole_field(Fi, Ri, Rj, Q1j)

    real(dp), dimension(3), intent(out) :: Fi
    real(dp), dimension(3), intent(in) :: Ri, Rj
    real(dp), dimension(3), intent(in) :: Q1j

    integer :: a, b
    integer, parameter :: k = 2
    real(dp), dimension(3) :: Rji
    real(dp), dimension(6) :: Tji
    real(dp), dimension(3,3) :: Tf

    Rji = Ri - Rj

    call Tk_tensor(Tji, k, Rji)

    call unpack_tensor(Tf, Tji)

    Fi = 0.0d0

    do a = 1, 3
        do b = 1, 3
            Fi(a) = Fi(a) + Tf(a,b) * Q1j(b)
        end do
    end do

end subroutine dipole_field

!------------------------------------------------------------------------------

subroutine quadrupole_field(Fi, Ri, Rj, Q2j)

    real(dp), dimension(3), intent(out) :: Fi
    real(dp), dimension(3), intent(in) :: Ri, Rj
    real(dp), dimension(6), intent(in) :: Q2j

    integer :: a, b, g
    integer, parameter :: k = 3
    real(dp), dimension(3) :: Rji
    real(dp), dimension(10) :: Tji
    real(dp), dimension(3,3) :: Q2f
    real(dp), dimension(3,3,3) :: Tf

    Rji = Ri - Rj

    call Tk_tensor(Tji, k, Rji)

    call unpack_tensor(Q2f, Q2j)
    call unpack_tensor(Tf, Tji)

    Fi = 0.0d0

    do a = 1, 3
        do b = 1, 3
            do g = 1, 3
                Fi(a) = Fi(a) - 0.5d0 * Tf(a,b,g) * Q2f(b,g)
            end do
        end do
    end do

end subroutine quadrupole_field

!------------------------------------------------------------------------------

subroutine octopole_field(Fi, Ri, Rj, Q3j)

    real(dp), dimension(3), intent(out) :: Fi
    real(dp), dimension(3), intent(in) :: Ri, Rj
    real(dp), dimension(10), intent(in) :: Q3j

    integer :: a, b, g, d
    integer, parameter :: k = 4
    real(dp), dimension(3) :: Rji
    real(dp), dimension(15) :: Tji
    real(dp), dimension(3,3,3) :: Q3f
    real(dp), dimension(3,3,3,3) :: Tf

    Rji = Ri - Rj

    call Tk_tensor(Tji, k, Rji)

    call unpack_tensor(Q3f, Q3j)
    call unpack_tensor(Tf, Tji)

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

end subroutine octopole_field

!------------------------------------------------------------------------------

subroutine response_matrix(B, invert, wrt2file)

! TODO: Damping schemes
!       Cutoff radius

    real(dp), dimension(:), intent(out) :: B
    logical, intent(in), optional :: invert, wrt2file

    logical :: exclude, lexist, inv, wrt
    integer :: info, lutemp
    integer :: i, j, k, l, m, n
    real(dp) :: damp, d5, d3
    real(dp) :: R, R2, R3, R5, T
    real(dp), dimension(3) :: Rij
    real(dp), dimension(6) :: alphainv

    if (present(invert)) then
        inv = invert
    else
        inv = .true.
    end if

    if (present(wrt2file)) then
        wrt = wrt2file
    else
        wrt = .true.
    end if

    B = 0.0d0

    inquire(file='pe_response_matrix.bin', exist=lexist)

    if (lexist) then
        call openfile('pe_response_matrix.bin', lutemp, 'old', 'unformatted')
        rewind(lutemp)
        read(lutemp) B
        close(lutemp)
    else
        m = 0
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            alphainv = alphas(:,i)
            call invert_packed_matrix(alphainv, 's')
            do l = 3, 1, -1
                do j = i, nsites
                    if (zeroalphas(j)) cycle
                    if (j == i) then
                        if (l == 3) then
                            do k = 1, l
                                B(m+k) = alphainv(k)
                            end do
                        else if (l == 2) then
                            do k = 1, l
                                B(m+k) = alphainv(3+k)
                            end do
                        else if (l == 1) then
                            do k = 1, l
                                B(m+k) = alphainv(5+k)
                            end do
                        end if
                        m = m + l
                    else
                        if (pe_nomb) then
                            m = m + 3
                            cycle
                        end if
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
                        R2 = R**2
                        R3 = R**3
                        R5 = R**5
! TODO: cutoff radius
!                        if (R > cutoff) then
!                            m = m + 3
!                            cycle
!                        end if
                        ! damping parameters
                        ! JPC A 102 (1998) 2399 & Mol. Sim. 32 (2006) 471
                        ! a = 2.1304 = damp
                        ! d3 = 1-(a²r²/2+ar+1)*exp(-ar)
                        ! d5 = 1-(a³r³/6+a²r²/2+ar+1)*exp(-ar)
                        if (pe_damp) then
                            d3 = damp**2 * R2 / 2.0d0
                            d3 = d3 + damp * R
                            d3 = d3 + 1.0d0
                            d3 = 1.0d0 - d3 * exp(-damp * R)
                            d5 = damp**3 * R3 / 6.0d0
                            d5 = d5 + damp**2 * R2 / 2.0d0
                            d5 = d5 + damp * R
                            d5 = d5 + 1.0d0
                            d5 = 1.0d0 - d5 * exp(-damp * R)
                            R3 = R3 * d3
                            R5 = R5 * d5
                        end if
                        if (l == 3) then
                            do k = 1, 3
                                T = 3.0d0 * Rij(1) * Rij(k) / R5
                                if (k == 1) T = T - 1.0d0 / R3
                                B(m+k) = - T
                            end do
                        else if (l == 2) then
                            do k = 1, 3
                                T = 3.0d0 * Rij(2) * Rij(k) / R5
                                if (k == 2) T = T - 1.0d0 / R3
                                B(m+k) = - T
                            end do
                        else if (l == 1) then
                            do k = 1, 3
                                T = 3.0d0 * Rij(3) * Rij(k) / R5
                                if (k == 3) T = T - 1.0d0 / R3
                                B(m+k) = - T
                            end do
                        end if
                        m = m + 3
                    end if
                end do
            end do
        end do
        if (inv) then
            call invert_packed_matrix(B, 's')
        end if
        if (wrt) then
            call openfile('pe_response_matrix.bin', lutemp, 'new', 'unformatted')
            rewind(lutemp)
            write(lutemp) B
            close(lutemp)
        end if
    end if

end subroutine response_matrix

!------------------------------------------------------------------------------

subroutine Tk_tensor(Tk, k, Rij)

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

end subroutine Tk_tensor

!------------------------------------------------------------------------------

subroutine unpack_tensor(Tf, Ts)

    real(dp), dimension(:), intent(in) :: Ts
    real(dp), dimension(*), intent(out) :: Tf

    if (size(Ts) == 1) then
        stop('Error in unpack_tensor: no unpacking necessary')
    else if (size(Ts) == 3) then
        stop('Error in unpack_tensor: no unpacking necessary')
    else if (size(Ts) == 6) then
        call unpack_2nd_order(Tf, Ts)
    else if (size(Ts) == 10) then
        call unpack_3rd_order(Tf, Ts)
    else if (size(Ts) == 15) then
        call unpack_4th_order(Tf, Ts)
    else if (size(Ts) > 15) then
        stop('Error in unpack_tensor: wrong size or not implemented')
    else
        stop('Error in unpack_tensor: packed tensor is wrong size')
    end if

end subroutine unpack_tensor

!------------------------------------------------------------------------------

subroutine unpack_2nd_order(Tf, Ts)

    real(dp), dimension(:), intent(in) :: Ts
    real(dp), dimension(3,3), intent(out) :: Tf

    Tf = 0.0d0

    Tf(1,1) = Ts(1); Tf(1,2) = Ts(2); Tf(1,3) = Ts(3)
    Tf(2,1) = Ts(2); Tf(2,2) = Ts(4); Tf(2,3) = Ts(5)
    Tf(3,1) = Ts(3); Tf(3,2) = Ts(5); Tf(3,3) = Ts(6)

end subroutine unpack_2nd_order

!------------------------------------------------------------------------------

subroutine unpack_3rd_order(Tf, Ts)

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

end subroutine unpack_3rd_order

!------------------------------------------------------------------------------

subroutine unpack_4th_order(Tf, Ts)

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

end subroutine unpack_4th_order

!------------------------------------------------------------------------------

subroutine Qk_integrals(Qk_ints, k, Rij, Qk, work)

    external :: Tk_integrals

    integer, intent(in) :: k
    real(dp), dimension(:,:), intent(out) :: Qk_ints
    real(dp), dimension(:), intent(in) :: Qk
    real(dp), dimension(3), intent(in) :: Rij
    real(dp), dimension(:), intent(inout) :: work

    integer :: i
    integer :: ncomps
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

    ncomps = size(Qk_ints, 2)

    ! get T^(k) integrals (incl. negative sign from electron density)
    call Tk_integrals(Qk_ints, ncomps*nnbas, k, Rij, work, size(work))

    ! get symmetry factors
    allocate(factors(ncomps)); factors = 0.0d0
    call symmetry_factors(factors, k)

    ! dot T^(k) integrals with multipole to get Q^(k) integrals
    do i = 1, ncomps
        Qk_ints(:,i) = taylor * factors(i) * Qk(i) * Qk_ints(:,i)
    end do

    deallocate(factors)

end subroutine Qk_integrals

!------------------------------------------------------------------------------

subroutine symmetry_factors(factors, k)

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

end subroutine symmetry_factors

!------------------------------------------------------------------------------

subroutine pe_save_density(density, nbas, coords, charges, work)

    external :: Tk_integrals

    integer, intent(in) :: nbas
    real(dp), dimension(:), intent(in) :: density
    real(dp), dimension(:), intent(in) :: charges
    real(dp), dimension(:,:), intent(in) :: coords
    real(dp), dimension(:), intent(inout) :: work

    integer :: i, j, l
    integer :: corenucs
    integer, parameter :: k = 0
    integer :: lucore, luden
    character(len=2) :: auoraa
    real(dp) :: Ene
    real(dp), dimension(:,:), allocatable :: Rc, Zc
    real(dp), dimension(:,:), allocatable :: full_density
    real(dp), dimension(:), allocatable :: T0_ints
    real(dp), dimension(:), allocatable :: Ffd
    real(dp), dimension(:,:), allocatable :: Ftmp

    ndens = 1
    nnbas = nbas * (nbas + 1) / 2

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
        Rc = Rc * aa2au
    end if

    ! unfold density matrix
    allocate(full_density(nbas,nbas)); full_density = 0.0d0
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
    allocate(Ftmp(3*npols,1), Ffd(3*npols)); Ftmp = 0.0d0
    call electron_fields(Ftmp, density, work)
    Ffd = Ftmp(:,1)
    call nuclear_fields(Ftmp(:,1))
    Ffd = Ffd + Ftmp(:,1)

    ! calculate nuclear - electron energy contribution
    allocate(T0_ints(nnbas)); T0_ints = 0.0d0
    Ene = 0.0d0
    do i = 1, corenucs
        call Tk_integrals(T0_ints, nnbas, k, Rc(:,i), work, size(work))
        T0_ints = Zc(1,i) * T0_ints
        Ene = Ene + dot(density, T0_ints)
    end do
    deallocate(T0_ints)

    ! save density, energy and field for subsequent calculations
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

end subroutine pe_save_density

!------------------------------------------------------------------------------

subroutine pe_intmol_twoints(nbas, work)

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

end subroutine pe_intmol_twoints

!------------------------------------------------------------------------------

subroutine pe_repulsion()

end subroutine pe_repulsion

!------------------------------------------------------------------------------

function nrm2(x)

    real(dp), external :: dnrm2

    real(dp) :: nrm2
    real(dp), dimension(:), intent(in) :: x

    integer :: n
    integer :: incx

    incx = 1

    n = size(x)

    nrm2 = dnrm2(n, x, incx)

end function nrm2

!------------------------------------------------------------------------------

function dot(x,y)

    real(dp), external :: ddot

    real(dp) :: dot
    real(dp), dimension(:), intent(in) :: x, y

    integer :: n
    integer :: incx, incy

    incx = 1
    incy = 1

    n = size(x)

    dot = ddot(n, x, incx, y, incy)

end function dot

!------------------------------------------------------------------------------

subroutine axpy(x, y, a)

    external :: daxpy

    real(dp), intent(in), optional :: a
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(inout) :: y

    real(dp) :: o_a
    integer :: n
    integer :: incx, incy

    if (present(a)) then
        o_a = a
    else
        o_a = 1.0d0
    end if

    incx = 1
    incy = 1

    n = size(x)

    call daxpy(n, o_a, x, incx, y, incy)

end subroutine axpy

!------------------------------------------------------------------------------

subroutine gemm(a, b, c, transa, transb, alpha, beta)

    external :: dgemm

    real(dp), intent(in), optional :: alpha, beta
    character(len=1), intent(in), optional :: transa, transb
    real(dp), dimension(:,:), intent(in) :: a, b
    real(dp), dimension(:,:) , intent(inout) :: c

    integer :: m, n, k, lda, ldb, ldc
    character(len=1) :: o_transa, o_transb
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

subroutine spmv(ap, x, y, uplo, alpha, beta)

    external :: dspmv
    intrinsic :: present, size

    real(dp), dimension(:), intent(in) :: ap
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(inout) :: y
    character(len=1), optional :: uplo
    real(dp), intent(in), optional :: alpha, beta

    integer :: n, incx, incy
    real(dp) :: o_alpha, o_beta
    character(len=1) :: o_uplo

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'U'
    end if

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

    incx = 1
    incy = 1
    n = size(x)

    call dspmv(o_uplo, n, o_alpha, ap, x, incx, o_beta, y, incy)

end subroutine spmv

!------------------------------------------------------------------------------

subroutine tpttr(ap, a, uplo)

    external :: dtpttr

    real(dp), dimension(:), intent(in) :: ap
    real(dp), dimension(:,:), intent(out) :: a
    character(len=1), intent(in), optional :: uplo

    integer :: n, lda, info
    character(len=1) :: o_uplo

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'L'
    end if

    n = size(a, 2)
    lda = max(1, size(a, 1))

    call dtpttr(o_uplo, n, ap, a, lda, info)

    if (info /= 0) call xerbla('tpttr', info)

end subroutine tpttr

!------------------------------------------------------------------------------

subroutine invert_packed_matrix(ap, sp)

    external :: dpptrf, dpptri, dsptrf, dsptri, xerbla

    character(*), optional :: sp
    real(dp), dimension(:), intent(inout) :: ap

    integer :: n, info
    integer, dimension(:), allocatable :: ipiv
    real(dp), dimension(:), allocatable :: work

    n = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * real(size(ap), dp)) - 1.0d0))

    if (.not. present(sp)) then
        sp = 's'
    end if

    if (sp == 'p') then
        call dpptrf('L', n, ap, info)
        if (info /= 0) call xerbla('pptrf', info)
        call dpptri('L', n, ap, info)
        if (info /= 0) call xerbla('pptri', info)
    else if (sp == 's') then
        allocate(ipiv(n), work(n))
        call dsptrf('L', n, ap, ipiv, info)
        if (info /= 0) call xerbla('sptrf', info)
        call dsptri('L', n, ap, ipiv, work, info)
        if (info /= 0) call xerbla('sptri', info)
        deallocate(ipiv, work)
    end if

end subroutine invert_packed_matrix

!------------------------------------------------------------------------------

subroutine solve(ap, b)

    external :: dspsv

    real(dp), dimension(:), intent(inout) :: ap
    real(dp), dimension(:,:), intent(inout) :: b

    integer :: n, nrhs, info
    integer, dimension(:), allocatable :: ipiv

    n = size(b, 1)
    nrhs = size(b, 2)

    allocate(ipiv(n))

    info = 0

    call dspsv('L', n, nrhs, ap, ipiv, b, n, info)

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
      end if
    end if

    do i = 21, 99
      inquire(unit=i, opened=lopen)
      if (lopen) then
        cycle
      else
        lunit = i
        open(unit=lunit, file=filename, status=stat, form=frmt)
        exit
      end if
    end do

    return

end subroutine openfile

!------------------------------------------------------------------------------

function elem2charge(elem)

    character(*) :: elem
    real(dp) :: elem2charge

    integer :: i
    real(dp), dimension(89) :: charges
    character(len=2), dimension(89) :: elements

    elements = (/ 'H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne',&
                & 'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' , 'Ca',&
                & 'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',&
                & 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',&
                & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',&
                & 'Sb', 'Te', 'I' , 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',&
                & 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',&
                & 'Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',&
                & 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'X'/)

    do i = 1, 88
        charges(i) = real(i, dp)
    end do
    charges(89) = 0.0d0

    do i = 1, 88
        if (elem == elements(i)) then
            elem2charge = charges(i)
            exit
        end if
    end do

end function elem2charge

!------------------------------------------------------------------------------

end module polarizable_embedding
