module polarizable_embedding

    use pe_precision
    use pe_blas_wrappers
    use pe_lapack_wrappers
    use pe_variables

#if defined(VAR_MPI)
#if defined(VAR_USE_MPIF)
    implicit none
#include "mpif.h"
#else
    use mpi
    implicit none
#endif
#else
    implicit none
#endif

    private

    intrinsic :: allocated, present, min, minval, max, maxval, size, cpu_time

    ! public subroutines/functions
    public :: pe_init, pe_master, pe_dalton_input
    public :: pe_save_density, pe_twoints
#if defined(VAR_MPI)
    public :: pe_mpi
#endif

! TODO:
! handle interface better, e.g. scale or remove higher order moments and pols
! damping of electric field from QM system?
! find better solution for electric field calculation from fragment densities
! higher order polarizabilities
! write list of publications which should be cited
! remove double zeroing and unecessary zeroing
! cutoffs and damping
! memory management
! add error catching

contains

!------------------------------------------------------------------------------

subroutine pe_init(lupri, coords, charges, dalwrk)

    ! Initialization routine for the PE module.
    integer :: lupri
    real(dp), dimension(:), intent(in), optional :: charges
    real(dp), dimension(:,:), intent(in), optional :: coords
    real(dp), dimension(:), target, intent(inout) :: dalwrk

    integer :: i, j, k, l
    integer :: idx, jdx, kdx, nidx
    integer, dimension(:), allocatable :: idxs
    logical, dimension(:), allocatable :: redists
    logical :: lexist
    real(dp) :: rclose

    ! Assume geometry optimization
    if (allocated(Rm) .and. allocated(Zm)) then
        Rm(:,:) = coords
        Zm(1,:) = charges
        return
    end if

    luout = lupri

    if (present(coords) .and. present(charges)) then
        qmnucs = size(charges)
        allocate(Rm(3,qmnucs), Zm(1,qmnucs))
        Rm(:,:) = coords
        Zm(1,:) = charges
    else if (present(coords) .and. .not. present(charges)) then
        stop 'ERROR in pe_init: coords present but charges missing'
    else if (.not. present(coords) .and. present(charges)) then
        stop 'ERROR in pe_init: charges present but coords missing'
    end if

    work => dalwrk

    ! setting up grid for MEP calculation
    if (pe_mep) then
        origin(1) = minval(Rm(1,:)) - xsize
        origin(2) = minval(Rm(2,:)) - ysize
        origin(3) = minval(Rm(3,:)) - zsize
        step(1) = 1.0 / xgrid
        step(2) = 1.0 / ygrid
        step(3) = 1.0 / zgrid
        xsteps = int((maxval(Rm(1,:)) + xsize - origin(1)) / step(1))
        ysteps = int((maxval(Rm(2,:)) + ysize - origin(2)) / step(2))
        zsteps = int((maxval(Rm(3,:)) + zsize - origin(3)) / step(3))
        npoints = 0
        npoints = xsteps * ysteps * zsteps
        allocate(mepgrid(3,npoints))
        l = 1
        do i = 1, xsteps
            do j = 1, ysteps
                do k = 1, zsteps
                    mepgrid(1,l) = origin(1) + (i - 1) * step(1)
                    mepgrid(2,l) = origin(2) + (j - 1) * step(2)
                    mepgrid(3,l) = origin(3) + (k - 1) * step(3)
                    l = l + 1
                end do
            end do
        end do
    end if

    call read_potential(trim(potfile))

    if (pe_sol) then
           call setup_solvent()
           call setup_cavity()
           if (.not. fixpva) call read_surface(trim(surfile))
    end if

    write(luout,'(//2x,a)') 'Polarizable embedding information'
    write(luout,'(2x,a)')   '---------------------------------'
    if (nsites > 0) then
        write(luout,'(/4x,a,i6)') 'Number of classical sites: ', nsites
    end if
    if (mulorder == 5) then
        write(luout,'(/4x,a)') 'Multipole moments upto 5th order.'
    else if (mulorder == 4) then
        write(luout,'(/4x,a)') 'Multipole moments upto 4th order.'
    else if (mulorder == 3) then
        write(luout,'(/4x,a)') 'Multipole moments upto 3rd order.'
    else if (mulorder == 2) then
        write(luout,'(/4x,a)') 'Multipole moments upto 2nd order.'
    else if (mulorder == 1) then
        write(luout,'(/4x,a)') 'Multipole moments upto 1st order.'
    else if (mulorder == 0) then
        write(luout,'(/4x,a)') 'Multipole moments upto 0th order.'
    end if
    if (pe_polar) then
        if (lpol(1)) write(luout,'(/4x,a)') 'Dipole-dipole polarizabilities.'
        if (pe_gspol) then
            write(luout,'(/4x,a)') 'Dynamic response from environment will be&
                                   & neglected during response calculation.'
        end if
        if (pe_nomb) then
            write(luout,'(/4x,a)') 'Many-body interactions will be neglected.'
        end if
        if (pe_iter) then
            write(luout,'(/4x,a)') 'Iterative solver for induced moments will&
                                   & be used'
            write(luout,'(4x,a,es7.1)') 'with convergence threshold: ', thriter
            if (pe_diis) then
                write(luout,'(4x,a)') 'and DIIS accelerator.'
            end if
            if (pe_redthr) then
                write(luout,'(/4x,a)') 'Using reduced threshold in first few&
                                       & SCF iterations.'
            end if
        else
            write(luout,'(/4x,a)') 'Direct solver for induced moments will be&
                                  & used.'
        end if
        if (pe_damp) then
            write(luout,'(/4x,a)') 'Interactions between inducible moments&
                                   & will be damped'
            write(luout,'(4x,a,f8.4)') 'using damping coefficient:', damp
        end if
    end if
    if (pe_fd) then
        write(luout,'(/4x,a,i4)') 'Number of fragment densities: ', nfds
        if (pe_repuls) then
            write(luout,'(/4x,a,f5.3)') 'Repulsion operator will be used for&
                                        & fragment densities'
        end if
    end if
    if (pe_sol) then
        write(luout,'(/4x,3a)') 'Continuum solvation in ', trim(solvent),&
                                & ' solvent.'
        write(luout,'(/4x,a,i5)') 'Number of surface points:', nsurp
    end if
    if (pe_restart) then
         write(luout,'(/4x,a)') 'Existing files will be used to restart if&
                                & possible.'
    end if

    ! handling sites near border
    ! -----------------------------------------------
    if (pe_border) then
        ! first locate all sites within given threshold of QM nuclei
        allocate(idxs(nsites))
        idxs = 0; nidx = 0
        do i = 1, qmnucs
            do j = 1, nsites
                lexist = .false.
                do k = 1, nidx
                    if (j == idxs(k)) then
                        lexist = .true.
                        exit
                    end if
                end do
                if (lexist) cycle
                if (nrm2(Rm(:,i) - Rs(:,j)) <= Rmin) then
                    nidx = nidx + 1
                    idxs(nidx) = j
                end if
            end do
        end do

        if (border_type == 'REMOVE') then
            write(luout,*) ''
            do i = 1, nidx
                write(luout,'(4x,a,i6,2x,a)') 'Removing parameters on site:',&
                                              & idxs(i) !, elems(idxs(i))
                if (lmul(0)) then
                    M0s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(1)) then
                    M1s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(2)) then
                    M2s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(3)) then
                    M3s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(4)) then
                    M4s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(5)) then
                    M5s(:,idxs(i)) = 0.0d0
                endif
                if (lpol(1)) then
                    P1s(:,idxs(i)) = 0.0d0
                end if
            end do
        else if (border_type == 'REDIST') then
            allocate(redists(nsites))
            redists = .false.
            do i = 1, nidx
                rclose = 1.0d10
                do j = 1, nsites
                    lexist = .false.
                    do k = 1, nidx
                        if (j == idxs(k)) then
                            lexist = .true.
                            exit
                        end if
                    end do
                    if (lexist) cycle
                    if (nrm2(Rs(:,idxs(i)) - Rs(:,j)) <= rclose) then
                        rclose = nrm2(Rs(:,idxs(i)) - Rs(:,j))
                        idx = j
                    end if
                end do
                if (lmul(0)) then
                    M0s(:,idx) = M0s(:,idx) + M0s(:,idxs(i)) / 3.0d0
                endif
!                if (lmul(1)) then
!                    M1s(:,idx) = M1s(:,idx) + M1s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(2)) then
!                    M2s(:,idx) = M2s(:,idx) + M2s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(3)) then
!                    M3s(:,idx) = M3s(:,idx) + M3s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(4)) then
!                    M4s(:,idx) = M4s(:,idx) + M4s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(5)) then
!                    M5s(:,idx) = M5s(:,idx) + M5s(:,idxs(i)) / 3.0d0
!                endif
!                if (lpol(1)) then
!                    P1s(:,idx) = P1s(:,idx) + P1s(:,idxs(i)) / 3.0d0
!                end if
                rclose = 1.0d10
                do j = 1, nsites
                    if (j == idx) cycle
                    lexist = .false.
                    do k = 1, nidx
                        if (j == idxs(k)) then
                            lexist = .true.
                            exit
                        end if
                    end do
                    if (lexist) cycle
                    if (nrm2(Rs(:,idxs(i)) - Rs(:,j)) <= rclose) then
                        rclose = nrm2(Rs(:,idxs(i)) - Rs(:,j))
                        jdx = j
                    end if
                end do
                if (lmul(0)) then
                    M0s(:,jdx) = M0s(:,jdx) + M0s(:,idxs(i)) / 3.0d0
                endif
!                if (lmul(1)) then
!                    M1s(:,jdx) = M1s(:,jdx) + M1s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(2)) then
!                    M2s(:,jdx) = M2s(:,jdx) + M2s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(3)) then
!                    M3s(:,jdx) = M3s(:,jdx) + M3s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(4)) then
!                    M4s(:,jdx) = M4s(:,jdx) + M4s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(5)) then
!                    M5s(:,jdx) = M5s(:,jdx) + M5s(:,idxs(i)) / 3.0d0
!                endif
!                if (lpol(1)) then
!                    P1s(:,jdx) = P1s(:,jdx) + P1s(:,idxs(i)) / 3.0d0
!                end if
                rclose = 1.0d10
                do j = 1, nsites
                    if (j == idx .or. j == jdx) cycle
                    lexist = .false.
                    do k = 1, nidx
                        if (j == idxs(k)) then
                            lexist = .true.
                            exit
                        end if
                    end do
                    if (lexist) cycle
                    if (nrm2(Rs(:,idxs(i)) - Rs(:,j)) <= rclose) then
                        rclose = nrm2(Rs(:,idxs(i)) - Rs(:,j))
                        kdx = j
                    end if
                end do
                if (lmul(0)) then
                    M0s(:,kdx) = M0s(:,kdx) + M0s(:,idxs(i)) / 3.0d0
                endif
!                if (lmul(1)) then
!                    M1s(:,kdx) = M1s(:,kdx) + M1s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(2)) then
!                    M2s(:,kdx) = M2s(:,kdx) + M2s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(3)) then
!                    M3s(:,kdx) = M3s(:,kdx) + M3s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(4)) then
!                    M4s(:,kdx) = M4s(:,kdx) + M4s(:,idxs(i)) / 3.0d0
!                endif
!                if (lmul(5)) then
!                    M5s(:,kdx) = M5s(:,kdx) + M5s(:,idxs(i)) / 3.0d0
!                endif
!                if (lpol(1)) then
!                    P1s(:,kdx) = P1s(:,kdx) + P1s(:,idxs(i)) / 3.0d0
!                end if
                if (lmul(0)) then
                    M0s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(1)) then
                    M1s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(2)) then
                    M2s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(3)) then
                    M3s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(4)) then
                    M4s(:,idxs(i)) = 0.0d0
                endif
                if (lmul(5)) then
                    M5s(:,idxs(i)) = 0.0d0
                endif
                if (lpol(1)) then
                    P1s(:,idxs(i)) = 0.0d0
                end if
                write(luout,'(/4x,a,i6)') 'Redistributing charges from site:',&
                                          & idxs(i)
                write(luout,'(4x,a,3i6)') 'to neighbouring sites:', idx, jdx,&
                                          & kdx
                write(luout,'(4x,a)') 'and removing all other parameters.'
                redists(idx) = .true.
                redists(jdx) = .true.
                redists(kdx) = .true.
            end do
            if (lmul(0)) then
                write(luout,'(/6x,a)') ' Resulting monopoles: '
                write(luout,'(6x,a)') '----------------------'
                do i = 1, nsites
                    if (redists(i)) then
                        write(luout,'(7x,a,1x,i6,2x,f9.4)') elems(i), i,&
                                                            & M0s(:,i)
                    end if
                end do
            end if
            deallocate(redists)
        end if
        deallocate(idxs)
    end if

    ! number of polarizabilities different from zero
    if (lpol(1)) then
        allocate(zeroalphas(nsites))
        do i = 1, nsites
            if (abs(maxval(P1s(:,i))) <= zero) then
                zeroalphas(i) = .true.
            else
                zeroalphas(i) = .false.
                npols = npols + 1
            end if
        end do
    end if

    nullify(work)

end subroutine pe_init

!------------------------------------------------------------------------------

subroutine pe_dalton_input(word, luinp, lupri)

    character(len=7), intent(inout) :: word
    integer, intent(in) :: luinp
    integer, intent(in) :: lupri

    integer :: i, j
    character(len=7) :: option
    character(len=2) :: auoraa

    luout = lupri

    do
        read(luinp,'(a7)') option
        call chcase(option)

        ! do a Polarizable Embedding calculation
        if (trim(option(2:)) == 'PEQM') then
            peqm = .true.
        else if (trim(option(2:)) == 'POTENT') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) potfile
            end if
        ! direct solver for induced moments
        else if (trim(option(2:)) == 'DIRECT') then
            pe_iter = .false.
        ! iterative solver for induced moments (default)
        else if (trim(option(2:)) == 'ITERAT') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) thriter
            end if
            pe_iter = .true.
        else if (trim(option(2:)) == 'DIIS') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) fixtol  
            end if
            pe_iter = .true.
            pe_diis = .true.
        ! mixed solver for induced moments (work in progress)
        else if (trim(option(2:)) == 'MIXED') then
            write(luout,*) 'WARNING: mixed solver is work in progress'
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) thriter
            end if
            pe_iter = .true.
            pe_mixed = .true.
        ! use reduced threshold in iterative induced moments solver
        else if (trim(option(2:)) == 'REDTHR') then
            pe_redthr = .true.
        ! handling sites near quantum-classical border
        else if (trim(option(2:)) == 'BORDER') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) border_type, Rmin, auoraa
                call chcase(border_type)
                if ((trim(border_type) /= 'REMOVE') .and.&
                   & (trim(border_type) /= 'REDIST')) then
                    stop 'ERROR: unknown handling of border sites!'
                end if
                call chcase(auoraa)
                if (trim(auoraa) == 'AA') Rmin = Rmin * aa2au
            end if
            pe_border = .true.
        ! damping interactions between inducible moments
        else if (trim(option(2:)) == 'DAMP') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) damp
            end if
            pe_damp = .true.
        ! neglect dynamic response from environment
        else if (trim(option(2:)) == 'GSPOL') then
            pe_gspol = .true.
        ! neglect many-body interactions
        else if (trim(option(2:)) == 'NOMB') then
            pe_nomb = .true.
!        ! use Gaussian broadened multipoles
!        else if (trim(option(2:)) == 'GAUSS') then
!            read(luinp,*) option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp,*) gauss_factor
!            end if
!            pe_gauss = .true.
        ! Use existing files for restart
        else if (trim(option(2:)) == 'RESTAR') then
            pe_restart = .true.
        ! calculate intermolecular two-electron integrals
        else if (trim(option(2:)) == 'TWOINT') then
            read(luinp,*) fdnucs
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) repfacs
            else
                repfacs(1) = 1.0d0
                repfacs(2) = 0.0d0
                repfacs(3) = 0.0d0
            end if
            pe_twoint = .true.
        ! save density matrix
        else if (trim(option(2:)) == 'SAVDEN') then
            pe_savden = .true.
        ! get fock matrix for repulsion potential
        else if (trim(option(2:)) == 'REPULS') then
            pe_repuls = .true.
        ! electrostatics from fragment densities
        else if (trim(option(2:)) == 'FD') then
            ! number of fragment densities
            read(luinp,*) nfds
            pe_fd = .true.
        ! skip electrostatics from fragment densities
        else if (trim(option(2:)) == 'NOFDES') then
            pe_fdes = .false.
        ! skip QM calculations, i.e. go directly into PE module
        else if (trim(option(2:)) == 'SKIPQM') then
            pe_skipqm = .true.
        ! calculate internal field
        else if (trim(option(2:)) == 'INFLD') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) ncrds
                allocate(crds(3,ncrds))
                do i = 1, ncrds 
                    read(luinp,*) (crds(j,i), j = 1, 3) 
                end do
            end if 
            pe_infld = .true.
        ! evaluate molecular electrostatic potential
        else if (trim(option(2:)) == 'MEP') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                do
                    read(luinp,*) option
                    call chcase(option)
                    if (trim(option(1:)) == 'GRID') then
                        read(luinp,*) xsize, xgrid, ysize, ygrid, zsize, zgrid
                    else if (trim(option(1:)) == 'FIELD') then
                        mep_field = .true.
                    else if (trim(option(1:)) == 'FLDNRM') then
                        mep_fldnrm = .true.
                    else if (trim(option(1:)) == 'EXTFLD') then
                        read(luinp,*) (extfld(i), i = 1, 3)
                        mep_extfld = .true.
                    else if (trim(option(1:)) == 'SKIPQM') then
                        mep_qmcube = .false.
                    else if (trim(option(1:)) == 'SKIPMUL') then
                        mep_mulcube = .false.
                    else if (option(1:1) == '.' .or. option(1:1) == '*') then
                        backspace(luinp)
                        exit
                    else if (option(1:1) == '!' .or. option(1:1) == '#') then
                        cycle
                    else
                        stop 'ERROR: unknown option present in .MEP section.'
                    end if
                end do
            end if
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) 
            end if
            pe_mep = .true.
        ! continuum solvation (COSMO) calculation 
        else if (trim(option(2:)) == 'SOLVAT') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) solvent
            else
                solvent = 'H2O'
            end if
            pe_sol = .true.
            pe_polar = .true.
        ! continuum solvation (COSMO) calculation 
        else if (trim(option(2:)) == 'EQSOL') then
            pe_noneq = .false.
        ! specify surface file
        else if (trim(option(2:)) == 'SURFAC') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) surfile
            end if
        ! FIXPVA
        else if (trim(option(2:)) == 'FIXPVA') then
             fixpva = .true.
        ! Read in maximum number of tessera
        else if (trim(option(2:)) == 'MAXTES') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) MXFFTS
            end if
        ! Read in number of tessera per atom
        else if (trim(option(2:)) == 'NTES') then
            read(luinp,*) option
            backspace(luinp)
            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.& 
               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
                read(luinp,*) NTSATM
                if (NTSATM .ne. 60 .and. ntsatm .ne. 240 .and. ntsatm .ne. 960) then
                   STOP 'NTSATM must be equal to 60, 240 or 960'
                end if
            end if
        ! verbose output
        else if (trim(option(2:)) == 'VERBOS') then
            pe_verbose = .true.
        ! debug output
        else if (trim(option(2:)) == 'DEBUG') then
            pe_debug = .true.
            pe_verbose = .true.
        else if (option(1:1) == '*') then
            word = option
            exit
        else if (option(1:1) == '!' .or. option(1:1) == '#') then
            cycle
        else
            write(luout,*) 'Unknown option:', option
        end if
    end do

! compatibility checks
    if (pe_sol .and. ndens .gt. 1) then
        stop 'Solvation model only implemented for ndens = 1'
    end if

end subroutine pe_dalton_input

!------------------------------------------------------------------------------

subroutine read_potential(filename)

    character(len=*) :: filename

    integer :: i, j, s
    integer :: nlines
    integer :: lupot
    real(dp) :: trace
    real(dp), dimension(21) :: temp
    character(len=2) :: auoraa
    character(len=80) :: word
    logical :: lexist

    inquire(file=filename, exist=lexist)
    if (lexist) then
        call openfile(filename, lupot, 'old', 'formatted')
    else
        return
!        write(luout,*) 'ERROR: input potential not found: ', filename
!        stop 'ERROR: input potential not found'
    end if

    do
        read(lupot,*, end=100) word

        if (trim(word) == 'coordinates') then
            read(lupot,*) nsites
            read(lupot,*) auoraa
            allocate(elems(nsites), Zs(1,nsites), Rs(3,nsites))
            do i = 1, nsites
                read(lupot,*) elems(i), (Rs(j,i), j = 1, 3)
                Zs(1,i) = elem2charge(elems(i))
            end do
            call chcase(auoraa)
            if (auoraa == 'AA') then
                Rs = Rs * aa2au
            end if
        else if (trim(word) == 'monopoles') then
            lmul(0) = .true.
            if (mulorder < 0) mulorder = 0
            allocate(M0s(1,nsites))
            M0s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, temp(1)
                M0s(1,s) = temp(1)
            end do
        else if (trim(word) == 'dipoles') then
            lmul(1) = .true.
            if (mulorder < 1) mulorder = 1
            allocate(M1s(3,nsites))
            M1s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 3)
                M1s(:,s) = temp(1:3)
            end do
        else if (trim(word) == 'quadrupoles') then
            lmul(2) = .true.
            if (mulorder < 2) mulorder = 2
            allocate(M2s(6,nsites))
            M2s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 6)
                ! remove trace
                trace = (temp(1) + temp(4) + temp(6)) / 3.0d0
                temp(1) = temp(1) - trace
                temp(4) = temp(4) - trace
                temp(6) = temp(6) - trace
                M2s(:,s) = temp(1:6)
            end do
        else if (trim(word) == 'octopoles') then
            lmul(3) = .true.
            if (mulorder < 3) mulorder = 3
            allocate(M3s(10,nsites))
            M3s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 10)
                trace = (temp(1) + temp(4) + temp(6)) / 5.0d0
                temp(1) = temp(1) - 3.0d0 * trace
                temp(4) = temp(4) - trace
                temp(6) = temp(6) - trace
                trace = (temp(2) + temp(7) + temp(9)) / 5.0d0
                temp(2) = temp(2) - trace
                temp(7) = temp(7) - 3.0d0 * trace
                temp(9) = temp(9) - trace
                trace = (temp(3) + temp(8) + temp(10)) / 5.0d0
                temp(3) = temp(3) - trace
                temp(8) = temp(8) - trace
                temp(10) = temp(10) - 3.0d0 * trace
                M3s(:,s) = temp(1:10)
            end do
        else if (trim(word) == 'hexadecapoles') then
            stop 'Hexadecapoles not supported currently'
            lmul(4) = .true.
            if (mulorder < 4) mulorder = 4
            allocate(M4s(15,nsites))
            M4s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 15)
                M4s(:,s) = temp(1:15)
            end do
        else if (trim(word) == 'ditriacontapoles') then
            stop 'Ditriacontapoles not supported currently'
            lmul(5) = .true.
            if (mulorder < 5) mulorder = 5
            allocate(M5s(21,nsites))
            M5s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 21)
                M5s(:,s) = temp(1:21)
            end do
        else if (trim(word) == 'isoalphas') then
            lpol(1) = .true.
            polorder = 1
            pe_polar = .true.
            if (.not. allocated(P1s)) then
                allocate(P1s(6,nsites))
                P1s = 0.0d0
            end if
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, temp(1)
                P1s(1,s) = temp(1)
                P1s(4,s) = temp(1)
                P1s(6,s) = temp(1)
            end do
        else if (trim(word) == 'alphas') then
            lpol(1) = .true.
            polorder = 1
            pe_polar = .true.
            if (.not. allocated(P1s)) then
                allocate(P1s(6,nsites))
                P1s = 0.0d0
            end if
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 6)
                P1s(:,s) = temp(1:6)
            end do
        else if (trim(word) == 'exclists' .or. trim(word) == 'exlists') then
            read(lupot,*) lexlst
            allocate(exclists(lexlst,nsites))
            do i = 1, nsites
                read(lupot,*) (exclists(j,i), j = 1, lexlst)
            end do
        else if (word(1:1) == '!' .or. word(1:1) == '#') then
            cycle
        end if
    end do

100 continue

    close(lupot)

   ! default exclusion list (everything polarizes everything)
    if (.not. allocated(exclists)) then
        lexlst = 1
        allocate(exclists(lexlst,nsites))
        do i = 1, nsites
            exclists(1,i) = i
        end do
    end if

end subroutine read_potential

!------------------------------------------------------------------------------

subroutine read_surface(filename)

    character(len=*) :: filename

    logical :: lexist
    integer :: i, j, lusurf
    character(len=2) :: auoraa
    real(dp), dimension(3) :: Rji
    
    inquire(file=filename, exist=lexist)
    if (lexist) then
        call openfile(filename, lusurf, 'old', 'formatted')
    else
        write(luout,*) 'ERROR: surface file not found: ', filename
        stop 'ERROR: surface file not found'
    end if

    read(lusurf,*) nsurp
    read(lusurf,*) auoraa

    allocate(Sp(3,nsurp))
    allocate(Sa(nsurp))

    do i = 1, nsurp
        read(lusurf,*) (Sp(j,i), j = 1, 3), Sa(i)
    end do

    close(lusurf)

    call chcase(auoraa)
    if (auoraa == 'AA') then
        Sp = Sp * aa2au
        Sa = Sa * aa2au2
    end if

    do i = 1, nsurp
       do j = 1, nsites
          Rji = Rs(:,j) - Sp(:,i)
          if (nrm2(Rji) < 1.2d0 ) then
              write(luout,'(a,f12.8)') 'WARNING: Cavity to close to classical&
                                       & site:', nrm2(Rji)
              write(luout,'(a,f12.8)') 'Surface point:', Sp(:,i)
              write(luout,'(a,f12.8)') 'Classical site:', Rs(:,j)
          end if
       end do
    end do

!    if (pe_debug) then
!       write(luout,*) 'Sp in read_surface, number of surface points:',nsurp
!       do i=1,nsurp
!              write (luout,*) i, Sp(:,i)
!       end do
!       write(luout,*) 'Sa in read_surface'
!       do i=1,nsurp
!          write (luout,*) i, Sa(i)
!       end do
!       write(luout,*) 'Sp in read_surface in AU'
!       do i=1,nsurp
!              write (luout,*) Sp(:,i)
!       end do
!       write(luout,*) 'Sa in read_surface in AU'
!       do i=1,nsurp
!          write (luout,*) Sa(i)
!       end do
!    end if
    
end subroutine read_surface

!------------------------------------------------------------------------------

subroutine pe_master(runtype, denmats, fckmats, nmats, energies, dalwrk)

    character(*), intent(in) :: runtype
    integer, intent(in) :: nmats
    real(dp), dimension(:), intent(in), optional :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats
    real(dp), dimension(:), intent(out), optional :: energies
    real(dp), dimension(:), target, intent(inout) :: dalwrk

#if defined(VAR_MPI)
    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)
#else
    myid = 0
    nprocs = 1
#endif

    work => dalwrk

    ! determine what to calculate and do consistency check
    if (runtype == 'fock') then
        fock = .true.
        energy = .false.
        response = .false.
        mep = .false.
        noneq = .false.
        london = .false.
        scfcycle = scfcycle + 1
        if (.not. present(fckmats)) then
            stop 'Output matrices are missing from input'
        else if (.not. present(energies)) then
            stop 'The energy variable is missing from input'
        end if
    else if (runtype == 'energy') then
        fock = .false.
        energy = .true.
        response = .false.
        mep = .false.
        noneq = .false.
        london = .false.
    else if (runtype == 'response') then
        if (pe_gspol) return
        fock = .false.
        energy = .false.
        response = .true.
        mep = .false.
        noneq = .false.
        london = .false.
        if (.not. present(fckmats)) then
            stop 'Output matrices are missing from input'
        end if
    else if (runtype == 'mep') then
        fock = .false.
        energy = .false.
        response = .false.
        mep = .true.
        noneq = .false.
        london = .false.
    else if (runtype == 'noneq') then
        fock = .false.
        energy = .false.
        response = .false.
        mep = .false.
!        noneq = .false.
        noneq = .true.
        london = .false.
        if (.not. present(fckmats)) then
            stop 'Output matrices are missing from input'
        end if
    else if (runtype == 'london') then
        fock = .false.
        energy = .false.
        response = .false.
        mep = .false.
        noneq = .false.
        london = .true.
        if (.not. present(fckmats)) then
            stop 'Output matrix missing from input.'
        end if
    else
        stop 'Could not determine calculation type.'
    end if

    ndens = nmats
    if (london) then
        if (ndens > 1) stop 'WARNING: LAO contribution for nmats > 1 requested'
        nbas = int(sqrt(size(fckmats) / 3.0d0))
    else
        nnbas = size(denmats) / ndens
        nbas = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * nnbas) - 1.0d0))
    endif

    if (.not. allocated(Epe)) allocate(Epe(ndens))
    Epe = 0.0d0
    if (.not. allocated(Ees)) allocate(Ees(0:5,ndens))
    Ees = 0.0d0
    if (.not. allocated(Efd)) allocate(Efd(3,ndens))
    Efd = 0.0d0
    if (.not. allocated(Epol)) allocate(Epol(3,ndens))
    Epol = 0.0d0
    if (.not. allocated(Esol)) allocate(Esol(3,ndens))
    Esol = 0.0d0

#if defined(VAR_MPI)
    if (myid == 0 .and. nprocs > 1) then
        call mpi_bcast(44, 1, impi, 0, comm, ierr)
        if (fock) then
            call mpi_bcast(1, 1, impi, 0, comm, ierr)
        else if (energy) then
            call mpi_bcast(2, 1, impi, 0, comm, ierr)
        else if (response) then
            call mpi_bcast(3, 1, impi, 0, comm, ierr)
        else if (mep) then
            call mpi_bcast(4, 1, impi, 0, comm, ierr)
        else if (london) then
            call mpi_bcast(5, 1, impi, 0, comm, ierr)
        end if

        call mpi_bcast(nbas, 1, impi, 0, comm, ierr)
        call mpi_bcast(ndens, 1, impi, 0, comm, ierr)
        if (.not. london) then
            call mpi_bcast(nnbas, 1, impi, 0, comm, ierr)
            call mpi_bcast(denmats, nnbas * ndens, rmpi, 0, comm, ierr)
        end if

        if (.not. synced) then
            call pe_sync()
        end if
    end if
#else
    site_start = 1
    site_finish = nsites
    surp_start = 1
    surp_finish = nsurp
    mep_start = 1
    mep_finish = npoints
#endif

    if (fock) then
        call pe_fock(denmats, fckmats, energies)
    else if (energy) then
        call pe_fock(denmats)
        write(luout,'(/6x,a)') 'Polarizable embedding energy contributions:'
        write(luout,'(5x,a)') '---------------------------------------------'
        if (mulorder >= 0) write(luout,'(/7x,a)') 'Electrostatic contributions:'
        if (lmul(0)) write(luout,'(9x,a16,5x,f20.12)') 'Monopoles       ',&
                                                       & Ees(0,1)
        if (lmul(1)) write(luout,'(9x,a16,5x,f20.12)') 'Dipoles         ',&
                                                       & Ees(1,1)
        if (lmul(2)) write(luout,'(9x,a16,5x,f20.12)') 'Quadrupoles     ',&
                                                       & Ees(2,1)
        if (lmul(3)) write(luout,'(9x,a16,5x,f20.12)') 'Octopoles       ',&
                                                       & Ees(3,1)
        if (lmul(4)) write(luout,'(9x,a16,5x,f20.12)') 'Hexadecapoles   ',&
                                                       & Ees(4,1)
        if (lmul(5)) write(luout,'(9x,a16,5x,f20.12)') 'Ditriacontapoles',&
                                                       & Ees(5,1)
        if (lpol(1)) then
            write(luout,'(/7x,a)') 'Polarization contributions:'
            write(luout,'(9x,a16,5x,f20.12)') 'Electronic      ', Epol(1,1)
            write(luout,'(9x,a16,5x,f20.12)') 'Nuclear         ', Epol(2,1)
            if (mulorder >= 0) then
                write(luout,'(9x,a16,5x,f20.12)') 'Multipole       ', Epol(3,1)
            end if
        end if
        if (pe_sol) then
            write(luout,'(/7x,a)') 'Continuum solvation contributions:'
            write(luout,'(9x,a16,5x,f20.12)') 'Electronic      ', Esol(1,1)
            write(luout,'(9x,a16,5x,f20.12)') 'Nuclear         ', Esol(2,1)
            if (mulorder >= 0) then
                write(luout,'(9x,a16,5x,f20.12)') 'Multipole       ', Esol(3,1)
            end if
        end if
        if (pe_fd) then
            write(luout,'(/7x,a)') 'Fragment density contributions:'
            write(luout,'(9x,a16,5x,f20.12)') 'Electrostatic   ', Efd(1,1)
            write(luout,'(9x,a16,5x,f20.12)') 'Polarization    ', Efd(2,1)
            if (pe_repuls) then
                write(luout,'(9x,a16,5x,f20.12)') 'Repulsion       ', Efd(3,1)
            end if
        end if
        write(luout,'(/3x,a18,9x,f20.12)') 'Total PE energy: ', Epe(1)
        if (pe_infld .and. nprocs > 1) then
            stop 'infld not parallelized yet.'
        else if (pe_infld) then
            if (allocated(crds)) then
                call pe_mappot2points(crds)
            else 
                call pe_mappot2points()
            end if
        endif
    else if (response) then
        call pe_fock(denmats, fckmats)
    else if (noneq) then
        call pe_fock(denmats, fckmats)
! TODO read in fckmats from file and subtrack from new fckmats 
    else if (mep) then
        if (ndens > 1) stop 'Not implemented for more than 1 density matrix'
        call pe_compute_mep(denmats)
    else if (london) then
        call pe_compute_london(fckmats)
    end if

!    deallocate(Epe, Ees, Epol, Esol, Efd)
    nullify(work)

end subroutine pe_master

!------------------------------------------------------------------------------

#if defined(VAR_MPI)
subroutine pe_mpi(dalwrk, runtype)

    real(dp), dimension(:), target, intent(inout) :: dalwrk
    integer :: runtype

    integer :: i, nwrk

    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)

    work => dalwrk

    nwrk = size(work)

    if (runtype == 1) then
        fock = .true.
        energy = .false.
        response = .false.
        mep = .false.
        scfcycle = scfcycle + 1
    else if (runtype == 2) then
        fock = .false.
        energy = .true.
        response = .false.
        mep = .false.
    else if (runtype == 3) then
        fock = .false.
        energy = .false.
        response = .true.
        mep = .false.
    else if (runtype == 4) then
        fock = .false.
        energy = .false.
        response = .false.
        mep = .true.
    else if (runtype == 5) then
        fock = .false.
        energy = .false.
        response = .false.
        mep = .false.
        london = .true.
    end if

    call mpi_bcast(nbas, 1, impi, 0, comm, ierr)
    call mpi_bcast(ndens, 1, impi, 0, comm, ierr)
    if (.not. london) then
        call mpi_bcast(nnbas, 1, impi, 0, comm, ierr)
        call mpi_bcast(work(1), nnbas * ndens, rmpi, 0, comm, ierr)
    end if

    if (.not. allocated(Epe)) allocate(Epe(ndens))
    Epe = 0.0d0
    if (.not. allocated(Ees)) allocate(Ees(0:5,ndens))
    Ees = 0.0d0
    if (.not. allocated(Efd)) allocate(Efd(3,ndens))
    Efd = 0.0d0
    if (.not. allocated(Epol)) allocate(Epol(3,ndens))
    Epol = 0.0d0
    if (.not. allocated(Esol)) allocate(Esol(3,ndens))
    Esol = 0.0d0

    if (.not. synced) then
        call pe_sync()
    end if

    if (fock) then
        call pe_fock(work(1:ndens*nnbas), work(ndens*nnbas+1:2*ndens*nnbas),&
                    & work(2*ndens*nnbas+1:2*ndens*nnbas+ndens))
    else if (energy) then
        call pe_fock(work(1:ndens*nnbas))
    else if (response) then
        call pe_fock(work(1:ndens*nnbas), work(ndens*nnbas+1:2*ndens*nnbas))
    else if (mep) then
        call pe_compute_mep(work(1:ndens*nnbas))
    else if (london) then
        call pe_compute_london(work(1:ndens*nbas*nbas))
    end if

!    deallocate(Epe, Ees, Epol, Esol, Efd)
    nullify(work)

end subroutine pe_mpi

!------------------------------------------------------------------------------

subroutine pe_sync()

    integer :: i, j, k
    integer :: quotient, remainder

    allocate(siteloops(0:nprocs))
    allocate(sitedists(0:nprocs-1))
    if (myid == 0) then
        quotient = nsites / nprocs
        sitedists = quotient
        if (nprocs * quotient < nsites) then
            remainder = nsites - nprocs * quotient
            do i = 1, remainder
                sitedists(i-1) = sitedists(i-1) + 1
            end do
        end if
        siteloops(0) = 0
        do i = 1, nprocs
            siteloops(i) = sum(sitedists(0:i-1))
        end do
        allocate(displs(0:nprocs))
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + 3 * sitedists(i-1)
        end do
    end if

    call mpi_bcast(nsites, 1, impi, 0, comm, ierr)
    call mpi_bcast(sitedists, nprocs, impi, 0, comm, ierr)
    call mpi_bcast(siteloops, nprocs + 1, impi, 0, comm, ierr)

    site_start = siteloops(myid) + 1
    site_finish = siteloops(myid+1)

    call mpi_bcast(qmnucs, 1, impi, 0, comm, ierr)

    if (myid /= 0) allocate(Zm(1,qmnucs))
    call mpi_bcast(Zm, qmnucs, rmpi, 0, comm, ierr)

    if (myid /= 0) allocate(Rm(3,qmnucs))
    call mpi_bcast(Rm, 3 * qmnucs, rmpi, 0, comm, ierr)

    if (myid /= 0) allocate(Rs(3,nsites))
    call mpi_bcast(Rs, 3 * nsites, rmpi, 0, comm, ierr)

    call mpi_bcast(pe_polar, 1, lmpi, 0, comm, ierr)
    call mpi_bcast(lpol, 1, lmpi, 0, comm, ierr)
    call mpi_bcast(pe_sol, 1, lmpi, 0, comm, ierr)

    if (pe_polar) then
        if (lpol(1)) then
            allocate(poldists(0:nprocs-1))
            if (myid == 0) then
                poldists = 0
                do i = 1, nprocs
                    do j = siteloops(i-1) + 1, siteloops(i)
                        if (zeroalphas(j)) then
                            continue
                        else
                            poldists(i-1) = poldists(i-1) + 1
                        end if
                    end do
                end do
            end if
            call mpi_bcast(poldists, nprocs, impi, 0, comm, ierr)
            call mpi_bcast(npols, 1, impi, 0, comm, ierr)
            call mpi_bcast(lexlst, 1, impi, 0, comm, ierr)
            if (myid /= 0) allocate(exclists(lexlst,nsites))
            call mpi_bcast(exclists, lexlst * nsites, impi, 0, comm, ierr)
            if (myid /= 0) allocate(zeroalphas(nsites))
            call mpi_bcast(zeroalphas, nsites, lmpi, 0, comm, ierr)
        end if
        if (pe_sol) then
            allocate(surploops(0:nprocs))
            allocate(surpdists(0:nprocs-1))
            if (myid == 0) then
                quotient = nsurp / nprocs
                surpdists = quotient
                if (nprocs * quotient < nsurp) then
                    remainder = nsurp - nprocs * quotient
                    do i = 1, remainder
                        surpdists(i-1) = surpdists(i-1) + 1
                    end do
                end if
                surploops(0) = 0
                do i = 1, nprocs
                    surploops(i) = sum(surpdists(0:i-1))
                end do
            end if
            call mpi_bcast(nsurp, 1, impi, 0, comm, ierr)
            call mpi_bcast(surpdists, nprocs, impi, 0, comm, ierr)
            call mpi_bcast(surploops, nprocs + 1, impi, 0, comm, ierr)
            surp_start = surploops(myid) + 1
            surp_finish = surploops(myid+1)
            if (myid /= 0) allocate(Sp(3,nsurp))
            call mpi_bcast(Sp, 3 * nsurp, rmpi, 0, comm, ierr)
            if (myid /= 0) allocate(Sa(nsurp))
            call mpi_bcast(Sa, nsurp, rmpi, 0, comm, ierr)
            call mpi_bcast(eps, 1, rmpi, 0, comm, ierr)
            call mpi_bcast(epsinf, 1, rmpi, 0, comm, ierr)
        end if
        call mpi_bcast(pe_iter, 1, lmpi, 0, comm, ierr)
        call mpi_bcast(pe_mixed, 1, lmpi, 0, comm, ierr)
        if (pe_iter .or. pe_mixed) then
            if (myid /= 0) allocate(P1s(6,nsites))
            call mpi_bcast(P1s, 6 * nsites, rmpi, 0, comm, ierr)
        end if
        call mpi_bcast(pe_nomb, 1, lmpi, 0, comm, ierr)
        call mpi_bcast(pe_damp, 1, lmpi, 0, comm, ierr)
        call mpi_bcast(damp, 1, rmpi, 0, comm, ierr)
    end if

    call mpi_bcast(pe_fd, 1, lmpi, 0, comm, ierr)
    call mpi_bcast(pe_restart, 1, lmpi, 0, comm, ierr)

    call mpi_bcast(mulorder, 1, impi, 0, comm, ierr)
    call mpi_bcast(lmul, 6, lmpi, 0, comm, ierr)

    if (lmul(0)) then
        if (myid /= 0) allocate(M0s(1,nsites))
        call mpi_bcast(M0s, nsites, rmpi, 0, comm, ierr)
    end if
    if (lmul(1)) then
        if (myid /= 0) allocate(M1s(3,nsites))
        call mpi_bcast(M1s, 3 * nsites, rmpi, 0, comm, ierr)
    end if
    if (lmul(2)) then
        if (myid /= 0) allocate(M2s(6,nsites))
        call mpi_bcast(M2s, 6 * nsites, rmpi, 0, comm, ierr)
    end if
    if (lmul(3)) then
        if (myid /= 0) allocate(M3s(10,nsites))
        call mpi_bcast(M3s, 10 * nsites, rmpi, 0, comm, ierr)
    end if
    if (lmul(4)) then
        if (myid /= 0) allocate(M4s(15,nsites))
        call mpi_bcast(M4s, 15 * nsites, rmpi, 0, comm, ierr)
    end if
    if (lmul(5)) then
        if (myid /= 0) allocate(M5s(21,nsites))
        call mpi_bcast(M5s, 21 * nsites, rmpi, 0, comm, ierr)
    end if

    if (mep) then
        call mpi_bcast(mep_field, 1, lmpi, 0, comm, ierr)
        call mpi_bcast(mep_fldnrm, 1, lmpi, 0, comm, ierr)
        call mpi_bcast(mep_extfld, 1, lmpi, 0, comm, ierr)
        call mpi_bcast(mep_qmcube, 1, lmpi, 0, comm, ierr)
        call mpi_bcast(mep_mulcube, 1, lmpi, 0, comm, ierr)
        allocate(meploops(0:nprocs))
        allocate(mepdists(0:nprocs-1))
        if (myid == 0) then
            quotient = npoints / nprocs
            mepdists = quotient
            if (nprocs * quotient < npoints) then
                remainder = npoints - nprocs * quotient
                do i = 1, remainder
                    mepdists(i-1) = mepdists(i-1) + 1
                end do
            end if
            meploops(0) = 0
            do i = 1, nprocs
                meploops(i) = sum(mepdists(0:i-1))
            end do
        end if
        call mpi_bcast(mepdists, nprocs, impi, 0, comm, ierr)
        call mpi_bcast(npoints, nprocs + 1, impi, 0, comm, ierr)
        if (myid == 0) then
            displs(0) = 0
            do i = 1, nprocs
             displs(i) = displs(i-1) + 3 * mepdists(i-1)
            end do
            call mpi_scatterv(mepgrid, 3 * mepdists, displs, rmpi,&
                             & mpi_in_place, 0, rmpi, 0, comm, ierr)
        else if (myid /= 0) then
            allocate(mepgrid(3,mepdists(myid)))
            call mpi_scatterv(0, 0, 0, rmpi, mepgrid, 3 * mepdists(myid),&
                             & rmpi, 0, comm, ierr)
        end if
    end if

    synced = .true.

end subroutine pe_sync
#endif

!------------------------------------------------------------------------------

subroutine pe_fock(denmats, fckmats, energies)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out), optional :: fckmats
    real(dp), dimension(:), intent(out), optional :: energies

    integer :: i
    logical :: es = .false.
    logical :: pol = .false.

    if (((mulorder >= 0) .or. pe_fd) .and. .not. response) es = .true.
    if (pe_polar) pol = .true.

    if (fock .or. response) fckmats = 0.0d0

    if (fock) then
        if (es) call pe_electrostatic(denmats, fckmats)
        if (pol) call pe_polarization(denmats, fckmats)
    else if (energy) then
        if (es) call pe_electrostatic(denmats)
        if (pol) call pe_polarization(denmats)
    else if (response .or. noneq) then
        if (pol) call pe_polarization(denmats, fckmats)
    end if

    if (fock .or. energy) then
        if (myid == 0) then
            do i = 1, ndens
                Epe(i) = sum(Ees(:,i)) + sum(Efd(:,i)) + sum(Epol(:,i)) +&
                       & sum(Esol(:,i))
            end do
            if (fock) energies = Epe
        end if
    end if

#if defined(VAR_MPI)
    if (fock .or. response) then
        if (myid == 0 .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, fckmats, ndens * nnbas, rmpi,&
                           & mpi_sum, 0, comm, ierr)
        else
            call mpi_reduce(fckmats, 0, ndens * nnbas, rmpi, mpi_sum, 0,&
                           & comm, ierr)
        end if
    end if
#endif

end subroutine pe_fock

!------------------------------------------------------------------------------

subroutine pe_electrostatic(denmats, fckmats)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats

    logical :: lexist
    integer :: lu
    integer :: i, j, k
    real(dp) :: Enuc, Etmp
    real(dp), dimension(ndens) :: Eel
    real(dp), dimension(:), allocatable :: tmpfcks

    if (myid == 0) then
        inquire(file='pe_electrostatics.bin', exist=lexist)
    end if
#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif
    if (lexist .and. fock) then
        if (myid == 0) then
            call openfile('pe_electrostatics.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Etmp, fckmats
            close(lu)
            do i = 1, ndens
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                Ees(0,i) = Ees(0,i) + dot(denmats(j:k), fckmats(j:k)) + Etmp
            end do
        end if
    else
        Etmp = 0.0d0
        if (lmul(0)) then
            if (fock) then
                call es_multipoles(M0s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M0s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(0,i) = Ees(0,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (lmul(1)) then
            if (fock) then
                call es_multipoles(M1s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M1s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(1,i) = Ees(1,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (lmul(2)) then
            if (fock) then
                call es_multipoles(M2s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M2s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(2,i) = Ees(2,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (lmul(3)) then
            if (fock) then
                call es_multipoles(M3s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M3s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(3,i) = Ees(3,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (lmul(4)) then
            if (fock) then
                call es_multipoles(M4s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M4s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(4,i) = Ees(4,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (lmul(5)) then
            if (fock) then
                call es_multipoles(M5s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M5s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(5,i) = Ees(5,i) + Eel(i) + Enuc
            end do
            Etmp = Etmp + Enuc
        end if
        if (myid == 0) then
            if (pe_fd) then
                if (fock) then
                    call es_fragment_densities(denmats, Eel, Enuc, fckmats)
                else if (energy) then
                    call es_fragment_densities(denmats, Eel, Enuc)
                end if
                do i = 1, ndens
                    Efd(1,i) = Efd(1,i) + Eel(i) + Enuc
                end do
                Etmp = Etmp + Enuc
            end if
        end if
        if (fock) then
#if defined(VAR_MPI)
            if (myid == 0 .and. nprocs > 1) then
                allocate(tmpfcks(ndens * nnbas))
                tmpfcks = fckmats
                call mpi_reduce(mpi_in_place, fckmats, ndens * nnbas, rmpi,&
                               & mpi_sum, 0, comm, ierr)
                call mpi_reduce(mpi_in_place, Etmp, 1, rmpi, mpi_sum, 0,&
                               & comm, ierr)
            else if (myid /= 0) then
                call mpi_reduce(fckmats, 0, ndens * nnbas, rmpi, mpi_sum, 0,&
                               & comm, ierr)
                call mpi_reduce(Etmp, 0, 1, rmpi, mpi_sum, 0, comm, ierr)
            end if
#endif
            if (myid == 0) then
                call openfile('pe_electrostatics.bin', lu, 'new', 'unformatted')
                rewind(lu)
                write(lu) Etmp, fckmats
                close(lu)
            end if
            if (myid == 0 .and. nprocs > 1) then
                fckmats = tmpfcks
                deallocate(tmpfcks)
            end if
        end if
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, Ees, 6 * ndens, rmpi, mpi_sum, 0,&
                           & comm, ierr)
        else if (myid /= 0) then
            call mpi_reduce(Ees, 0, 6 * ndens, rmpi, mpi_sum, 0, comm, ierr)
        end if
#endif
    end if

end subroutine pe_electrostatic

!------------------------------------------------------------------------------

subroutine es_fragment_densities(denmats, Eel, Enuc, fckmats)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out) :: Eel
    real(dp), intent(out) :: Enuc
    real(dp), dimension(:), intent(inout), optional :: fckmats

    integer :: i, j, k, l, m, n, o
    integer :: lufck, lexist, lu
    real(dp) :: Ene, Enn
    real(dp), dimension(ndens) :: Een, Eee
    real(dp), dimension(1) :: Tfm
    real(dp), dimension(3) :: Rfm
    real(dp), dimension(3*npols) :: temp
    real(dp), dimension(nnbas) :: fd_fckmat, fd_repmat
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
        read(lufck) fd_fckmat
        read(lufck) fd_repmat
        read(lufck) fdnucs
        allocate(Rfd(3,fdnucs), Zfd(1,fdnucs))
        read(lufck) Rfd, Zfd
        close(lufck)

        do j = 1, ndens
            l = (j - 1) * nnbas + 1
            m = j * nnbas
            if (pe_fdes) then
                if (fock) fckmats(l:m) = fckmats(l:m) + fd_fckmat
                Eee(j) = dot(denmats(l:m), fd_fckmat)
            end if
            if (pe_repuls) then
                if (fock) fckmats(l:m) = fckmats(l:m) + fd_repmat
                Efd(3,j) = Efd(3,j) + dot(denmats(l:m), fd_repmat)
            end if
        end do

        if (pe_fdes) then
            do j = 1, fdnucs
    !            gauss = (8.0d0  * gauss_factor) /&
    !                    ((P1s(1,j) + P1s(4,j) + P1s(6,j)) / 3.0d0)**(2.0d0/3.0d0)
                do k = 1, qmnucs
                    Rfm = Rm(:,k) - Rfd(:,j)
                    call Tk_tensor(Tfm, Rfm)
                    Enn = Enn + Zm(1,k) * Zfd(1,j) * Tfm(1)
                end do
    !            call Tk_integrals('es', Zfd_ints, nnbas, 1, Rfd(:,j)) 
    !            Zfd_ints = Zfd(1,j) * Zfd_ints
                call Mk_integrals(Zfd_ints, Rfd(:,j), Zfd(:,j))
                do m = 1, ndens
                    n = (m - 1) * nnbas + 1
                    o = m * nnbas
                    Een(m) = Een(m) + dot(denmats(n:o), Zfd_ints(:,1))
                    if (fock) fckmats(n:o) = fckmats(n:o) + Zfd_ints(:,1)
                end do
            end do
        end if

        deallocate(Rfd, Zfd)

        if (pe_fdes) then
            Enuc = Enuc + Ene + Enn
            do j = 1, ndens
                Eel(j) = Eel(j) + Een(j) + Eee(j)
            end do
        end if
    end do

end subroutine es_fragment_densities

!------------------------------------------------------------------------------

subroutine es_multipoles(Mks, denmats, Eel, Enuc, fckmats)

    real(dp), dimension(:,:), intent(in) :: Mks
    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats
    real(dp), dimension(:), intent(out) :: Eel
    real(dp), intent(out) :: Enuc

    integer :: site, ncomps
    integer :: i, j, k, l, m, n
    real(dp) :: taylor
    real(dp), dimension(3) :: Rsm
    real(dp), dimension(:), allocatable :: Tsm, symfacs
    real(dp), dimension(:,:), allocatable :: Mk_ints

    ncomps = size(Mks,1)

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * ncomps) - 1.0d0)) - 1

    if (mod(k,2) == 0) then
        taylor = 1.0d0 / factorial(k)
    else if (mod(k,2) /= 0) then
        taylor = - 1.0d0 / factorial(k)
    end if

    allocate(Tsm(ncomps))
    allocate(symfacs(ncomps))
    allocate(Mk_ints(nnbas,ncomps))

    Eel = 0.0d0; Enuc = 0.0d0

    i = 1
    do site = site_start, site_finish
        if (abs(maxval(Mks(:,site))) < zero) then
            i = i + 1
            cycle
        end if

        ! nuclei - multipole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,site)
            call Tk_tensor(Tsm, Rsm)
            call symmetry_factors(symfacs)
            do l = 1, ncomps
                Enuc = Enuc + taylor * symfacs(l) * Zm(1,j) * Mks(l,site) *&
                     & Tsm(l)
            end do
        end do

        ! electron - multipole interaction energy
        call Mk_integrals(Mk_ints, Rs(:,site), Mks(:,site))
        do l = 1, ndens
            m = (l - 1) * nnbas + 1
            n = l * nnbas
            Eel(l) = Eel(l) + dot(denmats(m:n), sum(Mk_ints, 2))
            if (fock) fckmats(m:n) = fckmats(m:n) + sum(Mk_ints, 2)
        end do
        i = i + 1
    end do

    deallocate(Tsm, symfacs, Mk_ints)

end subroutine es_multipoles

!------------------------------------------------------------------------------

subroutine pe_polarization(denmats, fckmats)

    external :: Tk_integrals

    real(dp), dimension(:), intent(in), optional :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats

    integer :: site, ndist, nrest
    integer :: i, j, k, l, m, n
    integer :: lunoneq
    logical :: skip
    real(dp) :: eps_fac
    real(dp), dimension(3) :: indtot
    real(dp), dimension(:,:), allocatable :: Vels
    real(dp), dimension(:), allocatable :: Vnucs, Vmuls, Vfds
    real(dp), dimension(:,:), allocatable :: Fels
    real(dp), dimension(:), allocatable :: Fnucs, Fmuls, Ffds
    real(dp), dimension(:,:), allocatable :: Mkinds, Fktots
    real(dp), dimension(:,:), allocatable :: Fel_ints, Vel_ints
    
    allocate(Mkinds(3*npols+nsurp,ndens), Fktots(3*npols+nsurp,ndens))
    Mkinds = 0.0d0; Fktots = 0.0d0
    if (lpol(1)) then
        allocate(Fels(3*npols,ndens))
        allocate(Fnucs(3*npols), Fmuls(3*npols), Ffds(3*npols))
        Ffds = 0.0d0
    end if
    if (pe_sol) then
        allocate(Vels(nsurp,ndens))
        allocate(Vnucs(nsurp), Vmuls(nsurp), Vfds(nsurp))
        if (noneq) then
!            eps_fac = (eps - epsinf) / ((eps - epsinf) - 1.0d0)
            eps_fac = ((eps - epsinf) - 1.0d0) / (eps - epsinf)
        else if (response) then
!            eps_fac = epsinf / (epsinf - 1.0d0)
            eps_fac = (epsinf - 1.0d0) / epsinf
        else if (fock .or. energy) then
            eps_fac = (eps - 1.0d0) / eps
        end if
    end if

    if (response .or. noneq) then
        if (lpol(1)) then
            call electron_fields(Fels, denmats)
            do i = 1, ndens
                Fktots(:3*npols,i) = Fels(:,i)
            end do
        end if
        if (pe_sol) then
            call electron_potentials(Vels, denmats)
            do i = 1, ndens
                if (pe_iter) then
                  ! Fktots(3*npols+1:,i) = - Vels(:,i)
                   Fktots(3*npols+1:,i) =  Vels(:,i)
                else
                   Fktots(3*npols+1:,i) = - eps_fac * Vels(:,i)
                end if
            end do 
        end if
        if (pe_debug) then
            write(luout,'(a)') 'Total electric field at polarizable sites:'
            do n = 1, ndens
                write(luout,'(a,i2)') 'From density matrix:', n
                i = 1
                do j = 1, npols
                    write(luout,'(i3,3f12.5)') j, Fktots(i,n), Fktots(i+1,n),&
                                             &Fktots(i+2,n)
                    i = i + 3
                end do
            end do
        end if
        call induced_moments(Mkinds, Fktots)
    else
        if (lpol(1)) then
            call electron_fields(Fels, denmats)
            call nuclear_fields(Fnucs)
            call multipole_fields(Fmuls)
        end if
        if (pe_sol) then 
            call electron_potentials(Vels, denmats)
            call nuclear_potentials(Vnucs)
            call multipole_potentials(Vmuls)
        end if 
        if (myid == 0) then
            if (pe_fd) then
! TODO frozen density potential
                call fragment_density_field(Ffds)
            end if
            do i = 1, ndens
                if (lpol(1)) then
                    Fktots(:3*npols,i) = Fels(:,i) + Fnucs + Fmuls + Ffds
                end if
                if (pe_sol) then
                    if (pe_iter) then
                       ! minus and eps_fac are included in the iterativ solver routine
                       Fktots(3*npols+1:,i) =  Vels(:,i) + Vnucs + Vmuls
                    else
                       Fktots(3*npols+1:,i) = - eps_fac * ( Vels(:,i) + Vnucs + Vmuls )
                    end if
                end if
            end do
        end if
        call induced_moments(Mkinds, Fktots)
        if (myid == 0) then
            do i = 1, ndens
                if (lpol(1)) then
                    Epol(1,i) = - 0.5d0 * dot(Mkinds(:3*npols,i), Fels(:,i))
                    Epol(2,i) = - 0.5d0 * dot(Mkinds(:3*npols,i), Fnucs)
                    Epol(3,i) = - 0.5d0 * dot(Mkinds(:3*npols,i), Fmuls)
                end if
                if (pe_sol) then
                    Esol(1,i) = 0.5d0 * dot(Mkinds(3*npols+1:,1), Vels(:,i))
                    Esol(2,i) = 0.5d0 * dot(Mkinds(3*npols+1:,1), Vnucs)
                    Esol(3,i) = 0.5d0 * dot(Mkinds(3*npols+1:,1), Vmuls)
                end if
                if (pe_fd) then
                    Efd(2,i) = - 0.5d0 * dot(Mkinds(:3*npols,i), Ffds)
                end if
            end do
        end if
    end if
    if (myid == 0) then
        if (pe_verbose .or. pe_debug) then
            write(luout,'(/2x,a)') ' Polarizable embedding information'
            write(luout,'(2x,a/)') '-----------------------------------'
        end if
        do i = 1, ndens
            if (pe_verbose .or. pe_debug) then
                write(luout,'(4x,a,i3)') 'Input density no.: ', i
            end if
            if (lpol(1) .and. pe_debug) then
                write(luout,'(/25x,a)') 'Induced dipole moments'
                write(luout,'(6x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
                k = 1
                do j = 1, npols
                    write(luout,'(4x,i6,3f15.8)') j, Mkinds(k:k+2,i)
                    k = k + 3
                end do
            end if
            if (lpol(1) .and. pe_verbose) then
                indtot = 0.0d0
                k = 1
                do j = 1, npols
                    indtot(1) = indtot(1) + Mkinds(k,i)
                    indtot(2) = indtot(2) + Mkinds(k+1,i)
                    indtot(3) = indtot(3) + Mkinds(k+2,i)
                end do
                write(luout,'(/22x,a)') 'Total induced dipole moment'
                write(luout,'(21x,a,14x,a,14x,a)') 'X', 'Y', 'Z'
                write(luout,'(10x,3f15.8/)') indtot
            end if
            if (pe_sol .and. pe_debug) then
                write(luout,'(28x,a)') 'Induced charges'
                do j = 3 * npols + 1, 3 * npols + nsurp
                    write(luout,'(25x,f15.8)') Mkinds(j,i)
                end do
            end if
            if (pe_sol .and. pe_verbose) then
                write(luout,'(4x,a,f15.8/)') 'Sum of induced charges: ',&
                                       & sum(Mkinds(3*npols+1:3*npols+nsurp,i))
            end if
        end do
    end if
#if defined(VAR_MPI)
    if (myid == 0 .and. nprocs > 1) then
        do i = 1, ndens
            if (lpol(1)) then
                displs(0) = 0
                do j = 1, nprocs
                    displs(j) = displs(j-1) + 3 * poldists(j-1)
                end do
                call mpi_scatterv(Mkinds(1:3*npols,i), 3 * poldists, displs,&
                                 & rmpi, mpi_in_place, 0, rmpi, 0, comm, ierr)
            end if
            if (pe_sol) then
                displs(0) = 0
                do j = 1, nprocs
                    displs(j) = displs(j-1) + surpdists(j-1)
                end do
                call mpi_scatterv(Mkinds(3*npols+1:3*npols+nsurp,i),&
                                 & surpdists, displs, rmpi, mpi_in_place, 0,&
                                 & rmpi, 0, comm, ierr)
            end if
        end do
    else if (myid /= 0) then
        do i = 1, ndens
            if (lpol(1)) then
                call mpi_scatterv(0, 0, 0, rmpi, Mkinds(1:3*npols,i),&
                                 & 3 * poldists(myid), rmpi, 0, comm, ierr)
            end if
            if (pe_sol) then
                call mpi_scatterv(0, 0, 0, rmpi,&
                                 & Mkinds(3*npols+1:3*npols+nsurp,i),& 
                                 & surpdists(myid), rmpi, 0, comm, ierr)
            end if
        end do
    end if
#endif

    if (fock .or. response) then
        if (lpol(1)) then
            allocate(Fel_ints(nnbas,3))
            i = 0
            do site = site_start, site_finish
                if (zeroalphas(site)) cycle
                call Tk_integrals('es', Fel_ints, nnbas, 3, Rs(:,site))
                do j = 1, 3
                    do k = 1, ndens
                        l = (k - 1) * nnbas + 1
                        m = k * nnbas
                        fckmats(l:m) = fckmats(l:m) - Mkinds(i+j,k) *&
                                     & Fel_ints(:,j)
                    end do
                end do
                i = i + 3
            end do
            deallocate(Fel_ints)
        end if
        if (pe_sol) then
            i = 1
            allocate(Vel_ints(nnbas,1))
            do site = surp_start, surp_finish
                call Tk_integrals('es', Vel_ints, nnbas, 1, Sp(:,site))
                do k = 1, ndens
                    l = (k - 1) * nnbas + 1
                    m = k * nnbas
                    fckmats(l:m) = fckmats(l:m) + Mkinds(3*npols+i,k) *&
                                 & Vel_ints(:,1)
                end do
                if (pe_noneq .and. fock) then
                    call openfile('pe_noneq.bin', lunoneq, 'unknown', 'unformatted')
                    rewind(lunoneq)
                    do k = 1, ndens
                        l = (k - 1) * nnbas + 1
                        m = k * nnbas
                        write(lunoneq) Mkinds(3*npols+i,k) * Vel_ints(:,1)
                    end do
                    close(lunoneq)
                end if
                i = i + 1
            end do
            deallocate(Vel_ints)
        end if
    else if (noneq) then
        if (pe_sol) then
            i = 1
            allocate(Vel_ints(nnbas,1))
            do site = surp_start, surp_finish
                call openfile('pe_noneq.bin', lunoneq, 'old', 'unformatted')
                rewind(lunoneq)
                do k = 1, ndens
                    l = (k - 1) * nnbas + 1
                    m = k * nnbas
                    write(lunoneq) fckmats(l:m)
                end do
                close(lunoneq)
                call Tk_integrals('es', Vel_ints, nnbas, 1, Sp(:,site))
                do k = 1, ndens
                    l = (k - 1) * nnbas + 1
                    m = k * nnbas
                    fckmats(l:m) = Mkinds(3*npols+i,k) *&
                                 & Vel_ints(:,1) - fckmats(l:m)
                end do
                i = i + 1
            end do
            deallocate(Vel_ints)
        end if
    end if

    deallocate(Mkinds, Fktots)
    if (lpol(1)) then
        deallocate(Fels)
        deallocate(Fnucs, Fmuls, Ffds)
    end if
    if (pe_sol) then
        deallocate(Vels)
        deallocate(Vnucs, Vmuls, Vfds)
    end if

end subroutine pe_polarization

!------------------------------------------------------------------------------

subroutine induced_moments(Mkinds, Fs)

    real(dp), dimension(:,:), intent(out) :: Mkinds
    real(dp), dimension(:,:), intent(in) :: Fs

    integer :: i, j, k

    if (pe_iter) then
        if (pe_diis) then
!            call pe_diis_solver(Mkinds, Fs)
            call pe_diis_solver_charges(Mkinds,Fs)
        else if (pe_mixed) then
            call mixed_solver(Mkinds, Fs)
        else
            call iterative_solver(Mkinds, Fs)
        end if
    else
        if (myid == 0) then
            call direct_solver(Mkinds, Fs)
        end if
    end if
!    do i = 1, 3*npols+nsurp
!        write(luout,*) 'Induced dipole i', Mkinds(i,1), i
!    end do

    ! check induced dipoles
    if (myid == 0) then
        do i = 1, ndens
            k = 1
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                if (nrm2(Mkinds(k:k+2,i)) > 1.0d0) then
                    write(luout,'(4x,a,i6)') 'Large induced dipole encountered&
                                             & at site:', j
                    write(luout,'(f10.4)') nrm2(Mkinds(k:k+2,i))
                end if
                k = k + 3
            end do
        end do
    end if

end subroutine induced_moments

!------------------------------------------------------------------------------

subroutine iterative_solver(Mkinds, Fs)

    real(dp), dimension(:,:), intent(out) :: Mkinds
    real(dp), dimension(:,:), intent(in) :: Fs

    integer :: lu, iter
    integer :: i, j, k, l, m, n, o, p, q
    logical :: exclude, lexist
    logical :: converged = .false.
    real(dp) :: fe = 1.0d0
    real(dp) :: ft = 1.0d0
    real(dp) :: R, R3, R5, Rd, ai, aj, norm, redthr
    real(dp), parameter :: d3i = 1.0d0 / 3.0d0
    real(dp), parameter :: d6i = 1.0d0 / 6.0d0
    real(dp), dimension(:), allocatable :: T, Rij, Ftmp, M1tmp

    if (myid == 0) then
        if (fock .and. scfcycle <= - nint(log10(thriter)) .and. .not.&
           & pe_restart .and. pe_redthr) then
            redthr = 10**(- log10(thriter) - scfcycle)
            write(luout,'(a)') 'INFO: using reduced threshold to determine&
                               & induced dipole moments.'
        else
            redthr = 1.0d0
        end if
    end if

    if (myid == 0) then
        inquire(file='pe_induced_dipoles.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif

    if (lexist .and. (fock .or. energy)) then
        if (myid == 0) then
            call openfile('pe_induced_dipoles.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Mkinds
            close(lu)
        end if
    end if

    allocate(T(6), Rij(3), Ftmp(3), M1tmp(3))
    do n = 1, ndens
        if (.not.lexist .or. response) then
#if defined(VAR_MPI)
            if (myid == 0 .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + 3 * poldists(i-1)
                end do
                call mpi_scatterv(Fs(:,n), 3 * poldists, displs, rmpi,&
                                 &mpi_in_place, 0, rmpi, 0, comm, ierr)
            else if (myid /= 0) then
                call mpi_scatterv(0, 0, 0, rmpi, Fs(:,n), 3 * poldists(myid),&
                                 & rmpi, 0, comm, ierr)
            end if
#endif

            l = 1
            do i = site_start, site_finish
                if (zeroalphas(i)) cycle
                call spmv(P1s(:,i), Fs(l:l+2,n), Mkinds(l:l+2,n), 'L')
                l = l + 3
            end do

#if defined(VAR_MPI)
            if (myid == 0 .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + 3 * poldists(i-1)
                end do
                call mpi_gatherv(mpi_in_place, 0, rmpi, Mkinds(:,n),&
                                & 3 * poldists, displs, rmpi, 0, comm, ierr)
            else if (myid /= 0) then
                call mpi_gatherv(Mkinds(:,n), 3 * poldists(myid), rmpi, 0, 0,&
                                & 0, rmpi, 0, comm, ierr)
            end if
#endif
        end if

        if (pe_nomb) cycle
! TODO (maybe)
! use only static field first two iterations
!        if (fock .and. (scfcycle <= 2) .and. .not. pe_restart) then
!            write(luout,'(a)') 'INFO: using static electric fields only to&
!                               & converge induced dipole moments.'
!            cycle
!        end if

#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + 3 * poldists(i-1)
            end do
            call mpi_scatterv(Mkinds(:,n), 3 * poldists, displs, rmpi,&
                             & mpi_in_place, 0, rmpi, 0, comm, ierr)
        else if (myid /= 0) then
            call mpi_scatterv(0, 0, 0, rmpi, Mkinds(:,n), 3 * poldists(myid),&
                             & rmpi, 0, comm, ierr)
        end if
#endif

        iter = 1
        do
            norm = 0.0d0
            l = 1
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                if (pe_damp) then
                    ai = (P1s(1,i) + P1s(4,i) + P1s(6,i)) * d3i
                end if
                m = 1
                Ftmp = 0.0d0
                do j = site_start, site_finish
                    if (zeroalphas(j)) cycle
                    exclude = .false.
                    do k = 1, lexlst
                        if (exclists(k,i) == exclists(1,j)) then
                            exclude = .true.
                            exit
                        end if
                    end do
                    if (i == j .or. exclude) then
                        m = m + 3
                        cycle
                    end if
                    Rij = Rs(:,j) - Rs(:,i)
                    R = nrm2(Rij)
                    R3 = R**3
                    R5 = R**5
                    ! damping parameters
                    ! JPC A 102 (1998) 2399 & Mol. Sim. 32 (2006) 471
                    ! a = 2.1304 = damp
                    ! u = R / (alpha_i * alpha_j)**(1/6)
                    ! fe = 1-(a²u²/2+au+1)*exp(-au)
                    ! ft = 1-(a³u³/6+a²u²/2+au+1)*exp(-au)
                    if (pe_damp) then
                        aj = (P1s(1,j) + P1s(4,j) + P1s(6,j)) * d3i
                        Rd = damp * R / (ai * aj)**d6i
                        fe = 1.0d0 - (0.5d0 * Rd**2 + Rd + 1.0d0) * exp(-Rd)
                        ft = fe - d6i * Rd**3 * exp(-Rd)
                    end if
                    q = 1
                    do o = 1, 3
                        do p = o, 3
                            T(q) = 3.0d0 * Rij(o) * Rij(p) * ft / R5
                            if (o == p) then
                                T(q) = T(q) - fe / R3
                            end if
                            q = q + 1
                        end do
                    end do
                    call spmv(T, Mkinds(m:m+2,n), Ftmp, 'L', 1.0d0, 1.0d0)
                    m = m + 3
                end do

#if defined(VAR_MPI)
                if (myid == 0 .and. nprocs > 1) then
                    call mpi_reduce(mpi_in_place, Ftmp, 3, rmpi, mpi_sum, 0,&
                                   & comm, ierr)
                else if (myid /= 0) then
                    call mpi_reduce(Ftmp, 0, 3, rmpi, mpi_sum, 0, comm, ierr)
                end if
#endif

                if (myid == 0) then
                    M1tmp = Mkinds(l:l+2,n)
                    Ftmp = Ftmp + Fs(l:l+2,n)
                    call spmv(P1s(:,i), Ftmp, Mkinds(l:l+2,n), 'L')
                    norm = norm + nrm2(Mkinds(l:l+2,n) - M1tmp)
                end if

#if defined(VAR_MPI)
                if (myid == 0 .and. nprocs > 1) then
                    displs(0) = 0
                    do j = 1, nprocs
                        displs(j) = displs(j-1) + 3 * poldists(j-1)
                    end do
                    call mpi_scatterv(Mkinds(:,n), 3*poldists, displs, rmpi,&
                                     & mpi_in_place, 0, rmpi, 0, comm, ierr)
                else if (myid /= 0) then
                    call mpi_scatterv(0, 0, 0, rmpi, Mkinds(:,n),&
                                     & 3*poldists(myid), rmpi, 0, comm, ierr)
                end if
#endif
                l = l + 3
            end do

            if (myid == 0) then
                if (norm < redthr * thriter) then
                    if (pe_verbose) then
                        write (luout,'(4x,a,i2,a)') 'Induced dipole moments&
                                                    & converged in ', iter,&
                                                    & ' iterations.'
                    end if
                    converged = .true.
                else if (iter > 50) then
                    write(luout,*) 'ERROR: could not converge induced dipole&
                                   & moments.'
                    stop 'ERROR: could not converge induced dipole moments.'
                else
                    converged = .false.
                    iter = iter + 1
                end if
            end if

#if defined(VAR_MPI)
            if (nprocs > 1) then
                call mpi_bcast(converged, 1, lmpi, 0, comm, ierr)
            end if
#endif
            if (converged) exit
        end do
    end do

    if (fock) then
        if (myid == 0) then
            call openfile('pe_induced_dipoles.bin', lu, 'unknown',&
                         & 'unformatted')
            rewind(lu)
            write(lu) Mkinds
            close(lu)
        end if
    end if

    deallocate(T, Rij, Ftmp, M1tmp)

end subroutine iterative_solver

!------------------------------------------------------------------------------

subroutine direct_solver(Mkinds, Fs)

    real(dp), dimension(:,:), intent(out) :: Mkinds
    real(dp), dimension(:,:), intent(in) :: Fs

    logical :: lexist
    integer :: lu, info, i
    integer, dimension(:), allocatable :: ipiv
    real(dp) :: anorm, rcond, prefac, eps_fac
    real(dp), dimension(:), allocatable :: B

    allocate(B((3*npols+nsurp)*(3*npols+nsurp+1)/2))

    inquire(file='pe_response_matrix.bin', exist=lexist)
    if (lexist .and. .not. noneq .and. .not. (response .and. rsp_first)) then
       
        call openfile('pe_response_matrix.bin', lu, 'old', 'unformatted')
        rewind(lu)
        if (chol) then
            read(lu) B
        else
            allocate(ipiv(3*npols+nsurp))
            read(lu) B, ipiv
        end if
        close(lu)

    else

        if (response .and. rsp_first) rsp_first = .false.

        call response_matrix(B)

        if (pe_debug) then
            anorm = lansp('1', B, 'L')
            write(luout,'(/4x,a,f15.8)') '1-norm of response matrix B: ', anorm
        end if

        if (chol) then

            call pptrf(B, 'L', info)
            if (info /= 0) then
                print *, 'Cholesky factorization failed. Trying regular...'
                if (pe_sol) then
                    allocate(ipiv(3*npols+nsurp))
                else
                    allocate(ipiv(3*npols))
                end if
                call sptrf(B, 'L', ipiv, info)
                if (info /= 0) then
                    stop 'ERROR: cannot create response matrix.'
                else
                    chol = .false.
                end if
            end if

        else

            if (pe_sol) then
                allocate(ipiv(3*npols+nsurp))
            else
                allocate(ipiv(3*npols))
            end if
            call sptrf(B, 'L', ipiv, info)
            if (info /= 0) then
                stop 'ERROR: cannot create response matrix.'
            end if

        end if

        if (pe_debug) then
            if (chol) then
                call ppcon(B, anorm, rcond, 'L')
                write(luout,'(4x,a,f15.8/)') 'Condition number of response&
                                             & matrix B: ', 1.0d0 / rcond
            else
                call spcon(B, ipiv, anorm, rcond, 'L')
                write(luout,'(4x,a,f15.8/)') 'Condition number of response&
                                             & matrix B: ', 1.0d0 / rcond
            end if
        end if

        call openfile('pe_response_matrix.bin', lu, 'unknown', 'unformatted')
        rewind(lu)
        if (chol) then
            write(lu) B
        else
            write(lu) B, ipiv
        end if
        close(lu)

    end if

    Mkinds = Fs
    if (chol) then
        call pptrs(B, Mkinds, 'L')
        deallocate(B)
    else
        call sptrs(B, Mkinds, ipiv, 'L')
        deallocate(B, ipiv)
    end if

end subroutine direct_solver

!------------------------------------------------------------------------------

subroutine mixed_solver(Mkinds, Fs)

    real(dp), dimension(:,:), intent(out) :: Mkinds
    real(dp), dimension(:,:), intent(in) :: Fs

    integer :: lu, iter, info, bsize
    integer :: i, j, k, l, m, n, o, p, q
    logical :: exclude, lexist
    logical :: converged = .false.
    real(dp) :: fe = 1.0d0
    real(dp) :: ft = 1.0d0
    real(dp) :: R, R3, R5, Rd, ai, aj, norm, redthr
    real(dp), parameter :: d3i = 1.0d0 / 3.0d0
    real(dp), parameter :: d6i = 1.0d0 / 6.0d0
    real(dp), dimension(:), allocatable :: T, Rij, Ftmp, M1tmp

    character(len=2) :: tno
    character(len=99) :: no
    integer, dimension(:), allocatable :: ipiv
    real(dp) :: anorm, rcond
    real(dp), dimension(:), allocatable :: B

    if (myid == 0) then
        if (fock .and. scfcycle <= - nint(log10(thriter)) .and. .not.&
           & pe_restart .and. pe_redthr) then
            redthr = 10**(- log10(thriter) - scfcycle)
            write(luout,'(a)') 'INFO: using reduced threshold to determine&
                               & induced dipole moments.'
        else
            redthr = 1.0d0
        end if
    end if

    if (myid == 0) then
        inquire(file='pe_induced_dipoles.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif

    if (lexist .and. (fock .or. energy)) then
        if (myid == 0) then
            call openfile('pe_induced_dipoles.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Mkinds
            close(lu)
        end if
    end if

    if (.not.lexist .or. response) then
#if defined(VAR_MPI)
        do n = 1, ndens
            if (myid == 0 .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + 3 * poldists(i-1)
                end do
                call mpi_scatterv(Fs(:,n), 3 * poldists, displs, rmpi,&
                                 &mpi_in_place, 0, rmpi, 0, comm, ierr)
            else if (myid /= 0) then
                call mpi_scatterv(0, 0, 0, rmpi, Fs(:,n), 3 * poldists(myid),&
                                 & rmpi, 0, comm, ierr)
            end if
        end do
#endif
        write(no,*) myid
        tno = trim(adjustl(no))
        bsize = 3 * site_finish - 3 * (site_start - 1)
        allocate(B(bsize*(bsize+1)/2), ipiv(bsize))
        inquire(file='pe_response_matrix_block'//tno//'.bin', exist=lexist)
        if (lexist) then
            call openfile('pe_response_matrix_block'//tno//'.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) B, ipiv
            close(lu)
        else
            call response_matrix_block(B, site_start, site_finish)
            call sptrf(B, 'L', ipiv, info)
            if (info /= 0) then
                stop 'ERROR: cannot create response matrix.'
            end if
            call openfile('pe_response_matrix_block'//tno//'.bin', lu, 'new', 'unformatted')
            rewind(lu)
            write(lu) B, ipiv
            close(lu)
        end if
        Mkinds(3*(site_start-1)+1:3*site_finish,:) = Fs(3*(site_start-1)+1:3*site_finish,:)
        call sptrs(B, Mkinds(3*(site_start-1)+1:3*site_finish,:), ipiv, 'L')
        deallocate(B, ipiv)
#if defined(VAR_MPI)
        do n = 1, ndens
            if (myid == 0 .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + 3 * poldists(i-1)
                end do
                call mpi_gatherv(mpi_in_place, 0, rmpi, Mkinds(:,n),&
                                & 3 * poldists, displs, rmpi, 0, comm, ierr)
            else if (myid /= 0) then
                call mpi_gatherv(Mkinds(:,n), 3 * poldists(myid), rmpi, 0, 0,&
                                & 0, rmpi, 0, comm, ierr)
            end if
        end do
#endif
    end if

    allocate(T(6), Rij(3), Ftmp(3), M1tmp(3))
    do n = 1, ndens
        if (pe_nomb) cycle
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + 3 * poldists(i-1)
            end do
            call mpi_scatterv(Mkinds(:,n), 3 * poldists, displs, rmpi,&
                             & mpi_in_place, 0, rmpi, 0, comm, ierr)
        else if (myid /= 0) then
            call mpi_scatterv(0, 0, 0, rmpi, Mkinds(:,n), 3 * poldists(myid),&
                             & rmpi, 0, comm, ierr)
        end if
#endif

        iter = 1
        do
            norm = 0.0d0
            l = 1
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                if (pe_damp) then
                    ai = (P1s(1,i) + P1s(4,i) + P1s(6,i)) * d3i
                end if
                m = 1
                Ftmp = 0.0d0
                do j = site_start, site_finish
                    if (zeroalphas(j)) cycle
                    exclude = .false.
                    do k = 1, lexlst
                        if (exclists(k,i) == exclists(1,j)) then
                            exclude = .true.
                            exit
                        end if
                    end do
                    if (i == j .or. exclude) then
                        m = m + 3
                        cycle
                    end if
                    Rij = Rs(:,j) - Rs(:,i)
                    R = nrm2(Rij)
                    R3 = R**3
                    R5 = R**5
                    ! damping parameters
                    ! JPC A 102 (1998) 2399 & Mol. Sim. 32 (2006) 471
                    ! a = 2.1304 = damp
                    ! u = R / (alpha_i * alpha_j)**(1/6)
                    ! fe = 1-(a²u²/2+au+1)*exp(-au)
                    ! ft = 1-(a³u³/6+a²u²/2+au+1)*exp(-au)
                    if (pe_damp) then
                        aj = (P1s(1,j) + P1s(4,j) + P1s(6,j)) * d3i
                        Rd = damp * R / (ai * aj)**d6i
                        fe = 1.0d0 - (0.5d0 * Rd**2 + Rd + 1.0d0) * exp(-Rd)
                        ft = fe - d6i * Rd**3 * exp(-Rd)
                    end if
                    q = 1
                    do o = 1, 3
                        do p = o, 3
                            T(q) = 3.0d0 * Rij(o) * Rij(p) * ft / R5
                            if (o == p) then
                                T(q) = T(q) - fe / R3
                            end if
                            q = q + 1
                        end do
                    end do
                    call spmv(T, Mkinds(m:m+2,n), Ftmp, 'L', 1.0d0, 1.0d0)
                    m = m + 3
                end do

#if defined(VAR_MPI)
                if (myid == 0 .and. nprocs > 1) then
                    call mpi_reduce(mpi_in_place, Ftmp, 3, rmpi, mpi_sum, 0,&
                                   & comm, ierr)
                else if (myid /= 0) then
                    call mpi_reduce(Ftmp, 0, 3, rmpi, mpi_sum, 0, comm, ierr)
                end if
#endif

                if (myid == 0) then
                    M1tmp = Mkinds(l:l+2,n)
                    Ftmp = Ftmp + Fs(l:l+2,n)
                    call spmv(P1s(:,i), Ftmp, Mkinds(l:l+2,n), 'L')
                    norm = norm + nrm2(Mkinds(l:l+2,n) - M1tmp)
                end if

#if defined(VAR_MPI)
                if (myid == 0 .and. nprocs > 1) then
                    displs(0) = 0
                    do j = 1, nprocs
                        displs(j) = displs(j-1) + 3 * poldists(j-1)
                    end do
                    call mpi_scatterv(Mkinds(:,n), 3*poldists, displs, rmpi,&
                                     & mpi_in_place, 0, rmpi, 0, comm, ierr)
                else if (myid /= 0) then
                    call mpi_scatterv(0, 0, 0, rmpi, Mkinds(:,n),&
                                     & 3*poldists(myid), rmpi, 0, comm, ierr)
                end if
#endif
                l = l + 3
            end do

            if (myid == 0) then
                if (norm < redthr * thriter) then
                    if (pe_verbose) then
                        write (luout,'(4x,a,i2,a)') 'Induced dipole moments&
                                                    & converged in ', iter,&
                                                    & ' iterations.'
                    end if
                    converged = .true.
                else if (iter > 50) then
                    write(luout,*) 'ERROR: could not converge induced dipole&
                                   & moments.'
                    stop 'ERROR: could not converge induced dipole moments.'
                else
                    converged = .false.
                    iter = iter + 1
                end if
            end if

#if defined(VAR_MPI)
            if (nprocs > 1) then
                call mpi_bcast(converged, 1, lmpi, 0, comm, ierr)
            end if
#endif
            if (converged) exit
        end do
    end do

    if (fock) then
        if (myid == 0) then
            call openfile('pe_induced_dipoles.bin', lu, 'unknown',&
                         & 'unformatted')
            rewind(lu)
            write(lu) Mkinds
            close(lu)
        end if
    end if

    deallocate(T, Rij, Ftmp, M1tmp)

end subroutine mixed_solver

!------------------------------------------------------------------------------

subroutine electron_potentials(Vels, denmats)

    external :: Tk_integrals

    real(dp), dimension(:,:), intent(out) :: Vels
    real(dp), dimension(:), intent(in) :: denmats

    logical :: skip
    integer :: site
    integer :: i, j, k, l, m
    real(dp), dimension(nnbas,1) :: Vel_ints

    Vels = 0.0d0

    i = 1
    do site = surp_start, surp_finish 
        call Tk_integrals('es', Vel_ints, nnbas, 1, Sp(:,site))
        do k = 1, ndens
            l = (k - 1) * nnbas + 1
            m = k * nnbas
            Vels(i,k) = dot(denmats(l:m), Vel_ints(:,1))
        end do
        i = i + 1
    end do

#if defined(VAR_MPI)
    if (myid == 0 .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + surpdists(i-1)
        end do
        do i = 1, ndens
            call mpi_gatherv(mpi_in_place, 0, rmpi, Vels(:,i), surpdists,&
                            & displs, rmpi, 0, comm, ierr)
        end do
    else if (myid /= 0) then
        do i = 1, ndens
            call mpi_gatherv(Vels(:,i), surpdists(myid), rmpi, 0, 0, 0, rmpi,&
                            & 0, comm, ierr)
        end do
    end if
#endif

    if (pe_debug) then
         do i = 1, nsurp 
             write (luout,*) 'i, Vels(i)' , i, Vels(i,:)
         end do
    end if

end subroutine electron_potentials

!------------------------------------------------------------------------------

subroutine electron_fields(Fels, denmats)

    external :: Tk_integrals

    real(dp), dimension(:,:), intent(inout) :: Fels
    real(dp), dimension(:), intent(in) :: denmats

    logical :: skip
    integer :: site
    integer :: i, j, k, l, m
    real(dp), dimension(nnbas,3) :: Fel_ints

    Fels = 0.0d0

    i = 0
    do site = site_start, site_finish
        if (zeroalphas(site)) cycle
        if (pe_savden) then
            skip = .false.
            do j = 1, qmnucs
                if (nrm2(Rs(:,site) - Rm(:,j)) <= 1.0d0) skip = .true.
            end do
            if (skip) then
                i = i + 3
                cycle
            end if
        end if
        call Tk_integrals('es', Fel_ints, nnbas, 3, Rs(:,site))
        do j = 1, 3
            do k = 1, ndens
                l = (k - 1) * nnbas + 1
                m = k * nnbas
                Fels(i+j,k) = dot(denmats(l:m), Fel_ints(:,j))
            end do
        end do
        i = i + 3
    end do

#if defined(VAR_MPI)
    if (myid == 0 .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + 3 * poldists(i-1)
        end do
        do i = 1, ndens
            call mpi_gatherv(mpi_in_place, 0, rmpi, Fels(:,i), 3 * poldists,&
                            & displs, rmpi, 0, comm, ierr)
        end do
    else if (myid /= 0) then
        do i = 1, ndens
            call mpi_gatherv(Fels(:,i), 3 * poldists(myid), rmpi, 0, 0, 0,&
                            & rmpi, 0, comm, ierr)
        end do
    end if
#endif

end subroutine electron_fields

!------------------------------------------------------------------------------

subroutine nuclear_potentials(Vnucs)

    real(dp), dimension(:), intent(out) :: Vnucs

    logical :: lexist, skip
    integer :: lu, site
    integer :: i, j
    real(dp), dimension(3) :: Rmsp
    real(dp), dimension(1) :: Tmsp

    if (myid == 0) then
        inquire(file='pe_nuclear_potential.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif

    if (lexist) then
        if (myid == 0) then
            call openfile('pe_nuclear_potential.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Vnucs
            close(lu)
        end if
    else
        Vnucs = 0.0d0
        i = 1
        do site = surp_start, surp_finish
            do j = 1, qmnucs
                Rmsp = Sp(:,site) - Rm(:,j)
                call Tk_tensor(Tmsp, Rmsp)
                Vnucs(i) = Vnucs(i) + Zm(1,j) * Tmsp(1)
            end do
            i = i + 1
        end do
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + surpdists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, 0, rmpi, Vnucs, surpdists, displs,&
                            & rmpi, 0, comm, ierr)
        else if (myid /= 0) then
            call mpi_gatherv(Vnucs, surpdists(myid), rmpi, 0, 0, 0, rmpi, 0,&
                            & comm, ierr)
        end if
#endif
        if (myid == 0) then
            call openfile('pe_nuclear_potential.bin', lu, 'new', 'unformatted')
            rewind(lu)
            write(lu) Vnucs
            close(lu)
        end if
    end if

    if (pe_debug) then
        do i=1, nsurp 
            write (luout,*) 'Vnucs(i)' ,i, Vnucs(i)
        end do
    end if

end subroutine nuclear_potentials

!------------------------------------------------------------------------------

subroutine nuclear_fields(Fnucs)

    real(dp), dimension(:), intent(out) :: Fnucs

    logical :: lexist, skip
    integer :: lu, site
    integer :: i, j, k
    real(dp), dimension(3) :: Rms, Tms

    if (myid == 0) then
        inquire(file='pe_nuclear_field.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif

    if (lexist) then
        if (myid == 0) then
            call openfile('pe_nuclear_field.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Fnucs
            close(lu)
        end if
    else
        Fnucs = 0.0d0
        i = 0
        do site = site_start, site_finish
            if (zeroalphas(site)) cycle
            if (pe_savden) then
                skip = .false.
                do j = 1, qmnucs
                    if (nrm2(Rs(:,site) - Rm(:,j)) <= 1.0d0) skip = .true.
                end do
                if (skip) then
                    i = i + 3
                    cycle
                end if
            end if
            do j = 1, qmnucs
                Rms = Rs(:,site) - Rm(:,j)
                call Tk_tensor(Tms, Rms)
                do k = 1, 3
                    Fnucs(i+k) = Fnucs(i+k) - Zm(1,j) * Tms(k)
                end do
            end do
            i = i + 3
        end do
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + 3 * poldists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, 0, rmpi, Fnucs, 3*poldists, displs,&
                            & rmpi, 0, comm, ierr)
        else if (myid /= 0) then
            call mpi_gatherv(Fnucs, 3*poldists(myid), rmpi, 0, 0, 0, rmpi, 0,&
                            & comm, ierr)
        end if
#endif
        if (myid == 0) then
            call openfile('pe_nuclear_field.bin', lu, 'new', 'unformatted')
            rewind(lu)
            write(lu) Fnucs
            close(lu)
        end if
    end if

end subroutine nuclear_fields

!------------------------------------------------------------------------------

subroutine fragment_density_field(Ffd)

    real(dp), dimension(:), intent(out) :: Ffd

    integer :: i
    integer :: lu
    character(len=99) :: ci
    character(len=80) :: filename
    real(dp), dimension(3*npols) :: Ftmp

    Ffd = 0.0d0

    do i = 1, nfds
        Ftmp = 0.0d0
        write(ci,*) i
        ci = adjustl(ci)
        filename = 'pe_fock_'//trim(ci)//'.bin'
        call openfile(trim(filename), lu, 'old', 'unformatted')
        rewind(lu)
        read(lu) Ftmp
        close(lu)
        Ffd = Ffd + Ftmp
    end do

end subroutine fragment_density_field

!------------------------------------------------------------------------------

subroutine multipole_potentials(Vmuls)

    real(dp), dimension(:), intent(out) :: Vmuls

    logical :: exclude, lexist
    integer :: lu
    integer :: i, j, k
    real(dp), dimension(3) :: Rji

    if (myid == 0) then
        inquire(file='pe_multipole_potential.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif

    if (lexist) then
        if (myid == 0) then
            call openfile('pe_multipole_potential.bin', lu, 'old',&
                         & 'unformatted')
            rewind(lu)
            read(lu) Vmuls
            close(lu)
        end if
    else
        Vmuls = 0.0d0
        k = 1
        do i = 1, nsurp
            do j = site_start, site_finish
                Rji = Sp(:,i) - Rs(:,j)
                if (lmul(0)) then
                    if (abs(maxval(M0s(:,j))) >= zero) then
                        call multipole_potential(Vmuls(k), Rji, M0s(:,j))
                    end if
                end if
                if (lmul(1)) then
                    if (abs(maxval(M1s(:,j))) >= zero) then
                        call multipole_potential(Vmuls(k), Rji, M1s(:,j))
                    end if
                end if
                if (lmul(2)) then
                    if (abs(maxval(M2s(:,j))) >= zero) then
                        call multipole_potential(Vmuls(k), Rji, M2s(:,j))
                    end if
                end if
                if (lmul(3)) then
                    if (abs(maxval(M3s(:,j))) >= zero) then
                        call multipole_potential(Vmuls(k), Rji, M3s(:,j))
                    end if
                end if
                if (lmul(4)) then
                    if (abs(maxval(M4s(:,j))) >= zero) then
                        call multipole_potential(Vmuls(k), Rji, M4s(:,j))
                    end if
                end if
                if (lmul(5)) then
                    if (abs(maxval(M5s(:,j))) >= zero) then
                        call multipole_potential(Vmuls(k), Rji, M5s(:,j))
                    end if
                end if
            end do
            k = k + 1
        end do
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, Vmuls, nsurp, rmpi, mpi_sum, 0,&
                           & comm, ierr)
        else if (myid /= 0) then
            call mpi_reduce(Vmuls, 0, nsurp, rmpi, mpi_sum, 0, comm, ierr)
        end if
#endif
        if (myid == 0) then
            call openfile('pe_multipole_potential.bin', lu, 'new',&
                         & 'unformatted')
            rewind(lu)
            write(lu) Vmuls
            close(lu)
        end if
     end if

!    if (pe_debug) then
!        do i=1, nsurp 
!            write (luout,*) 'Vmuls(i)' ,i, Vmuls(i)
!        end do
!    end if

end subroutine multipole_potentials

!------------------------------------------------------------------------------

subroutine multipole_fields(Fmuls)

    real(dp), dimension(:), intent(out) :: Fmuls

    logical :: exclude, lexist
    integer :: lu
    integer :: i, j, k, l
    real(dp), dimension(3) :: Rji

    if (myid == 0) then
        inquire(file='pe_multipole_field.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif

    if (lexist) then
        if (myid == 0) then
            call openfile('pe_multipole_field.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Fmuls
            close(lu)
        end if
    else
        Fmuls = 0.0d0
        l = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            do j = site_start, site_finish
                if (i == j) then
                    cycle
                end if
                exclude = .false.
                do k = 1, lexlst
                    if (exclists(k,i) == exclists(1,j)) then
                        exclude = .true.
                        exit
                    end if
                end do
                if (exclude) cycle
! TODO: cutoff???
                Rji = Rs(:,i) - Rs(:,j)
                if (lmul(0)) then
                    if (abs(maxval(M0s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M0s(:,j))
                    end if
                end if
                if (lmul(1)) then
                    if (abs(maxval(M1s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M1s(:,j))
                    end if
                end if
                if (lmul(2)) then
                    if (abs(maxval(M2s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M2s(:,j))
                    end if
                end if
                if (lmul(3)) then
                    if (abs(maxval(M3s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M3s(:,j))
                    end if
                end if
                if (lmul(4)) then
                    if (abs(maxval(M4s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M4s(:,j))
                    end if
                end if
                if (lmul(5)) then
                    if (abs(maxval(M5s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M5s(:,j))
                    end if
                end if
            end do
            l = l + 3
        end do
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, Fmuls, 3*npols, rmpi, mpi_sum, 0,&
                           & comm, ierr)
        else if (myid /= 0) then
            call mpi_reduce(Fmuls, 0, 3*npols, rmpi, mpi_sum, 0, comm, ierr)
        end if
#endif
        if (myid == 0) then
            call openfile('pe_multipole_field.bin', lu, 'new', 'unformatted')
            rewind(lu)
            write(lu) Fmuls
            close(lu)
        end if
     end if

end subroutine multipole_fields

!------------------------------------------------------------------------------

subroutine multipole_potential(Vi, Rji, Mkj)

    real(dp), intent(inout) :: Vi
    real(dp), dimension(3), intent(in) :: Rji
    real(dp), dimension(:), intent(in) :: Mkj

    integer :: k
    integer :: a, x, y, z
    real(dp) :: taylor

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * size(Mkj)) - 1.0d0)) - 1

    if (mod(k,2) == 0) then
        taylor = 1.0d0 / factorial(k)
    else if (mod(k,2) /= 0) then
        taylor = - 1.0d0 / factorial(k)
    end if

    a = 1
    do x = k, 0, -1
        do y = k, 0, -1
            do z = k, 0, -1
                if (x+y+z /= k) cycle
                Vi = Vi + taylor * symfac(x,y,z) * T(Rji,x,y,z) * Mkj(a)
                a = a + 1
            end do
        end do
     end do

end subroutine multipole_potential

!------------------------------------------------------------------------------

subroutine multipole_field(Fi, Rji, Mkj)

    real(dp), dimension(3), intent(inout) :: Fi
    real(dp), dimension(3), intent(in) :: Rji
    real(dp), dimension(:), intent(in) :: Mkj

    integer :: k
    integer :: a, b, c, x, y, z
    real(dp) :: taylor

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * size(Mkj)) - 1.0d0))

    if (mod(k,2) == 0) then
        taylor = 1.0d0 / factorial(k-1)
    else if (mod(k,2) /= 0) then
        taylor = - 1.0d0 / factorial(k-1)
    end if

    a = 1; b = 1; c = 1
    do x = k, 0, -1
        do y = k, 0, -1
            do z = k, 0, -1
                if (x+y+z /= k) cycle
                if (x /= 0) then
                    Fi(1) = Fi(1) + taylor * symfac(x-1,y,z) * T(Rji,x,y,z) *&
                          & Mkj(a)
                    a = a + 1
                end if
                if (y /= 0) then
                    Fi(2) = Fi(2) + taylor * symfac(x,y-1,z) * T(Rji,x,y,z) *&
                          & Mkj(b)
                    b = b + 1
                end if
                if (z /= 0) then
                    Fi(3) = Fi(3) + taylor * symfac(x,y,z-1) * T(Rji,x,y,z) *&
                          & Mkj(c)
                    c = c + 1
                end if
            end do
        end do
     end do

end subroutine multipole_field

!------------------------------------------------------------------------------

subroutine response_matrix(B)

! TODO: Cutoff radius

    real(dp), dimension(:), intent(out) :: B

    logical :: exclude
    integer :: info
    integer :: i, j, k, l, m, n, o
    integer, dimension(3) :: ipiv
    real(dp), parameter :: fourpi = 4.0d0 * pi
    real(dp), parameter :: d3i = 1.0d0 / 3.0d0
    real(dp), parameter :: d6i = 1.0d0 / 6.0d0
    real(dp), parameter :: dism0 = 0.4d0 * aa2au
    real(dp) :: fe = 1.0d0
    real(dp) :: ft = 1.0d0
    real(dp) :: Rd, ai, aj
    real(dp) :: R, R3, R5, T
    real(dp) :: prefac, eps_fac
    real(dp), dimension(3) :: Rij
    real(dp), dimension(6) :: P1inv

    B = 0.0d0

!   if (pe_sol) then
!       if (noneq) then
!           eps_fac = (eps - epsinf) / ((eps - epsinf) - 1.0d0)
!       else if (response) then
!           eps_fac = epsinf / (epsinf - 1.0d0)
!       else if (fock) then
!           eps_fac = eps / (eps - 1.0d0)
!       end if
!       prefac = 1.07d0 * eps_fac
!   end if

    m = 0
    do i = 1, nsites
        if (zeroalphas(i)) cycle
        P1inv = P1s(:,i)
        call pptrf(P1inv, 'L', info)
        if (info /= 0) then
            P1inv = P1s(:,i)
            call sptrf(P1inv, 'L', ipiv, info)
            if (info /= 0) then
                stop 'ERROR: could not factorize polarizability.'
            else if (chol) then
                chol = .false.
            end if
            call sptri(P1inv, ipiv, 'L')
        else
            call pptri(P1inv, 'L')
        end if
        if (pe_damp) then
            ai = (P1s(1,i) + P1s(4,i) + P1s(6,i)) * d3i
        end if
        do l = 3, 1, -1
            do j = i, nsites
                if (zeroalphas(j)) cycle
                if (j == i) then
                    if (l == 3) then
                        do k = 1, l
                            B(m+k) = P1inv(k)
                        end do
                    else if (l == 2) then
                        do k = 1, l
                            B(m+k) = P1inv(3+k)
                        end do
                    else if (l == 1) then
                        do k = 1, l
                            B(m+k) = P1inv(5+k)
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
                        if (exclists(k,i) == exclists(1,j)) then
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
! TODO: cutoff radius
!                        if (R > cutoff) then
!                            m = m + 3
!                            cycle
!                        end if
                    ! damping parameters
                    ! JPC A 102 (1998) 2399 & Mol. Sim. 32 (2006) 471
                    ! a = 2.1304 = damp
                    ! u = R / (alpha_i * alpha_j)**(1/6)
                    ! fe = 1-(a²u²/2+au+1)*exp(-au)
                    ! ft = 1-(a³u³/6+a²u²/2+au+1)*exp(-au)
                    if (pe_damp) then
                        aj = (P1s(1,j) + P1s(4,j) + P1s(6,j)) * d3i
                        Rd = damp * R / (ai * aj)**d6i
                        fe = 1.0d0 - (0.5d0 * Rd**2 + Rd + 1.0d0) * exp(-Rd)
                        ft = fe - d6i * Rd**3 * exp(-Rd)
                    end if
                    if (l == 3) then
                        do k = 1, 3
                            T = 3.0d0 * Rij(1) * Rij(k) * ft / R5
                            if (k == 1) T = T - fe / R3
                            B(m+k) = - T
                        end do
                    else if (l == 2) then
                        do k = 1, 3
                            T = 3.0d0 * Rij(2) * Rij(k) * ft / R5
                            if (k == 2) T = T - fe / R3
                            B(m+k) = - T
                        end do
                    else if (l == 1) then
                        do k = 1, 3
                            T = 3.0d0 * Rij(3) * Rij(k) * ft / R5
                            if (k == 3) T = T - fe / R3
                            B(m+k) = - T
                        end do
                    end if
                    m = m + 3
                end if ! i /= j
            end do  ! do j = i, nsites
            if (pe_sol) then
                do j = 1, nsurp
                    Rij = Sp(:,j) - Rs(:,i)
                    R3 = nrm2(Rij)**3
                    if (l == 3) then
                        B(m+1) = - Rij(1) / R3
                    else if (l == 2) then
                        B(m+1) = - Rij(2) / R3
                    else if (l == 1) then
                        B(m+1) = - Rij(3) / R3
                    end if
                    m = m + 1
                end do ! j = 1, nsurp
            end if ! pe_sol
        end do  ! do l = 3, 1, -1
    end do  ! do i = 1, nsites
    if (pe_sol) then
        do i = 1, nsurp
            B(m+1) = 1.07d0 * sqrt(fourpi / Sa(i))
            m = m + 1
            do j = i + 1, nsurp
                Rij = Sp(:,j) - Sp(:,i)
                R = nrm2(Rij)
                if (fixpva) then
                   if ( R .le. dism0 ) then
                      B(m+1) = 0.5d0*(3.0d0/dism0 - R**2/dism0**3) 
                   else
                      B(m+1) = 1.0d0 / R
                   end if
                   m = m + 1
                else
                    B(m+1) = 1.0d0 / R
                    m = m + 1
                end if
            end do ! j = i + 1, nsurp
        end do ! i = 1, nsurp
    end if

    if (pe_debug) then
        do i = 1, (3 * npols + nsurp) * (3 * npols + nsurp + 1) / 2
            write (luout,*) 'Response matrix(i)',i, B(i)
        end do
    end if

end subroutine response_matrix

!------------------------------------------------------------------------------

subroutine response_matrix_block(B, start, finish)

! TODO: Cutoff radius

    integer :: start, finish
    real(dp), dimension(:), intent(out) :: B

    logical :: exclude
    integer :: info
    integer :: i, j, k, l, m, n, o
    integer, dimension(3) :: ipiv
    real(dp), parameter :: d3i = 1.0d0 / 3.0d0
    real(dp), parameter :: d6i = 1.0d0 / 6.0d0
    real(dp) :: fe = 1.0d0
    real(dp) :: ft = 1.0d0
    real(dp) :: Rd, ai, aj
    real(dp) :: R, R3, R5, T, eps_fac
    real(dp), dimension(3) :: Rij
    real(dp), dimension(6) :: P1inv

    B = 0.0d0

!    if (pe_sol) then
!        if (noneq) then
!            eps_fac = (eps - epsinf) / ((eps - epsinf) - 1.0d0)
!        else if (response) then
!            eps_fac = epsinf / (epsinf - 1.0d0)
!        else if (fock) then
!            eps_fac = eps / (eps - 1.0d0)
!        end if
!    end if

    m = 0
    do i = start, finish
        if (zeroalphas(i)) cycle
        P1inv = P1s(:,i)
        call sptrf(P1inv, 'L', ipiv, info)
        if (info /= 0) then
            stop 'ERROR: could not factorize polarizability.'
        end if
        call sptri(P1inv, ipiv, 'L')
        if (pe_damp) then
            ai = (P1s(1,i) + P1s(4,i) + P1s(6,i)) * d3i
        end if
        do l = 3, 1, -1
            do j = i, finish
                if (zeroalphas(j)) cycle
                if (j == i) then
                    if (l == 3) then
                        do k = 1, l
                            B(m+k) = P1inv(k)
                        end do
                    else if (l == 2) then
                        do k = 1, l
                            B(m+k) = P1inv(3+k)
                        end do
                    else if (l == 1) then
                        do k = 1, l
                            B(m+k) = P1inv(5+k)
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
                        if (exclists(k,i) == exclists(1,j)) then
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
! TODO: cutoff radius
!                        if (R > cutoff) then
!                            m = m + 3
!                            cycle
!                        end if
                    ! damping parameters
                    ! JPC A 102 (1998) 2399 & Mol. Sim. 32 (2006) 471
                    ! a = 2.1304 = damp
                    ! u = R / (alpha_i * alpha_j)**(1/6)
                    ! fe = 1-(a²u²/2+au+1)*exp(-au)
                    ! ft = 1-(a³u³/6+a²u²/2+au+1)*exp(-au)
                    if (pe_damp) then
                        aj = (P1s(1,j) + P1s(4,j) + P1s(6,j)) * d3i
                        Rd = damp * R / (ai * aj)**d6i
                        fe = 1.0d0 - (0.5d0 * Rd**2 + Rd + 1.0d0) * exp(-Rd)
                        ft = fe - d6i * Rd**3 * exp(-Rd)
                    end if
                    if (l == 3) then
                        do k = 1, 3
                            T = 3.0d0 * Rij(1) * Rij(k) * ft / R5
                            if (k == 1) T = T - fe / R3
                            B(m+k) = - T
                        end do
                    else if (l == 2) then
                        do k = 1, 3
                            T = 3.0d0 * Rij(2) * Rij(k) * ft / R5
                            if (k == 2) T = T - fe / R3
                            B(m+k) = - T
                        end do
                    else if (l == 1) then
                        do k = 1, 3
                            T = 3.0d0 * Rij(3) * Rij(k) * ft / R5
                            if (k == 3) T = T - fe / R3
                            B(m+k) = - T
                        end do
                    end if
                    m = m + 3
                end if ! i /= j
            end do  ! do j = i, nsites
!            if (pe_sol) then
!                do j = 1, nsurp
!                    Rij = Sp(:,j) - Rs(:,i)
!                    R3 = nrm2(Rij)**3
!                    if (l == 3) then
!                        B(m+1) = - Rij(1) / R3
!                    else if (l == 2) then
!                        B(m+1) = - Rij(2) / R3
!                    else if (l == 1) then
!                        B(m+1) = - Rij(3) / R3
!                    end if
!                    m = m + 1
!                end do ! j = 1, nsurp
!            end if ! pe_sol
        end do  ! do l = 3, 1, -1
    end do  ! do i = 1, nsites
!    if (pe_sol) then
!        do i = 1, nsurp
!            B(m+1) = 1.07d0 * eps_fac * sqrt((4.0d0 * pi) / Sa(i))
!            m = m + 1
!            do j = i + 1, nsurp
!                Rij = Sp(:,j) - Sp(:,i)
!                R = nrm2(Rij)
!                B(m+1) = eps_fac / R
!                m = m + 1
!            end do ! j = i + 1, nsurp
!        end do ! i = 1, nsurp
!    end if

!    if (pe_debug) then
!        do i = 1, 3 * npols + nsurp * (3 * npols + nsurp + 1) / 2
!            write (luout,*) 'Response matrix(i)',i, B(i)
!        end do
!    end if

end subroutine response_matrix_block

!------------------------------------------------------------------------------

subroutine Tk_coefficients

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    integer :: i, j, k, l, m, n

! TODO
!    i = max(mulorder, polorder)
!    allocate(Cnij(2*i+3,0:i+1,0:i+1))
    allocate(Cnij(2*5+3,0:5+1,0:5+1))

    Cnij = 0
    Cnij(:,0,0) = 1
    do n = 1, 2*5+3
        if (mod(n,2) == 0) cycle
        do i = 1, 5+1
            if (mod(i,2) /= 0) then
                k = i - 1
            else if (mod(i,2) == 0) then
                k = i
            end if
            do j = 0, i
                if (mod(i+j,2) /= 0) cycle
                if (j == 0) then
                    Cnij(n,i,j) = Cnij(n,i-1,j+1)
                else if (j /= i) then
                    Cnij(n,i,j) = (j + 1) * Cnij(n,i-1,j+1)
                    Cnij(n,i,j) = Cnij(n,i,j) - (n + k) * Cnij(n,i-1,j-1)
                    k = k + 2
                else if (j == i) then
                    Cnij(n,i,j) = - (n + k) * Cnij(n,i-1,j-1)
                end if
            end do
        end do
    end do

end subroutine Tk_coefficients

!------------------------------------------------------------------------------

function T(Rij, x, y, z)

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    integer, intent(in) :: x, y, z
    real(dp), dimension(3), intent(in) :: Rij

    integer :: l, m, n
    real(dp) :: T
    real(dp) :: R, Cx, Cy, Cz

    if (.not. allocated(Cnij)) call Tk_coefficients

    R = nrm2(Rij)

    T = 0.0d0

    do l = 0, x
        Cx = Cnij(1,x,l)*(Rij(1)/R)**l
        do m = 0, y
            Cy = Cx * Cnij(l+x+1,y,m)*(Rij(2)/R)**m
            do n = 0, z
                Cz = Cy * Cnij(l+x+m+y+1,z,n)*(Rij(3)/R)**n
                T = T + Cz
            end do
        end do
    end do
    T = T / R**(x+y+z+1)

end function T

!------------------------------------------------------------------------------

subroutine Tk_tensor(Tk, Rij)

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    real(dp), dimension(:), intent(out) :: Tk
    real(dp), dimension(3), intent(in) :: Rij

    integer :: k, i
    integer :: x, y, z

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * size(Tk)) - 1.0d0)) - 1

    i = 1
    do x = k, 0, -1
        do y = k, 0, -1
            do z = k, 0, -1
                if (x+y+z /= k) cycle
                Tk(i) = T(Rij, x, y, z)
                i = i + 1
            end do
        end do
    end do

end subroutine Tk_tensor

!------------------------------------------------------------------------------

subroutine Mk_integrals(Mk_ints, Rij, Mk)

    external :: Tk_integrals

    real(dp), dimension(:,:), intent(out) :: Mk_ints
    real(dp), dimension(:), intent(in) :: Mk
    real(dp), dimension(3), intent(in) :: Rij

    integer :: i, k
    integer :: ncomps
    real(dp) :: taylor
    real(dp), dimension(:), allocatable :: factors

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * size(Mk)) - 1.0d0)) - 1

    if (mod(k,2) == 0) then
        taylor = 1.0d0 / factorial(k)
    else if (mod(k,2) /= 0) then
        taylor = - 1.0d0 / factorial(k)
    end if

    ncomps = size(Mk_ints, 2)

    call Tk_integrals('es', Mk_ints, nnbas, ncomps, Rij)

    ! get symmetry factors
    allocate(factors(ncomps))
    call symmetry_factors(factors)

    ! dot T^(k) integrals with multipole to get M^(k) integrals
    do i = 1, ncomps
        Mk_ints(:,i) = taylor * factors(i) * Mk(i) * Mk_ints(:,i)
    end do

    deallocate(factors)

end subroutine Mk_integrals

!------------------------------------------------------------------------------

subroutine symmetry_factors(factors)

    real(dp), dimension(:), intent(out) :: factors

    integer :: idx, x, y, z, k

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * size(factors)) - 1.0d0)) - 1

    idx = 1
    do x = k, 0, -1
        do y = k, 0, -1
            do z = k, 0, -1
                if (x+y+z /= k) cycle
                factors(idx) = symfac(x, y, z)
                idx = idx + 1
            end do
        end do
     end do

end subroutine symmetry_factors

!------------------------------------------------------------------------------

function symfac(i, j, k)

    ! trinomial coefficient

    integer, intent(in) :: i, j, k

    integer :: symfac

    symfac = factorial(i+j+k) / (factorial(i) * factorial(j) * factorial(k))

end function symfac

!------------------------------------------------------------------------------

recursive function factorial(n) result(nfact)

    ! Clive Page, http://www.star.le.ac.uk/~cgp/f90course/f90.html

    integer, intent(in) :: n

    integer :: nfact

    if (n > 0) then
        nfact = n * factorial(n-1)
    else
        nfact = 1
    end if

end function factorial

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

function elem2charge(elem) result(charge)

    character(len=*), intent(in) :: elem

    integer :: i
    real(dp) :: charge
    character(len=2), dimension(112) :: elements

    elements = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',&
                & 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',&
                & 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',&
                & 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',&
                & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',&
                & 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',&
                & 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',&
                & 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',&
                & 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',&
                & 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',&
                & 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',&
                & 'Rg', 'Cn' /)

    if (elem == 'X') then
        charge = 0.0d0
        return
    end if

    do i = 1, 112
        if (elem == trim(elements(i))) then
            charge = real(i, dp)
            exit
        else
            charge = 0.0d0
        end if
    end do

end function elem2charge

!------------------------------------------------------------------------------

function charge2vdw(charge) result(vdw)

    real(dp), intent(in) :: charge

    integer :: i
    real(dp) :: vdw
    real(dp), dimension(19) :: radii

    radii = (/ 1.20, 1.40, 2.20, 1.90, 1.80, 1.70, 1.60, 1.55, 1.50, 1.54,&
            &  2.40, 2.20, 2.10, 2.10, 1.95, 1.80, 1.80, 1.88, 1.90 /)

    i = nint(charge)
    if (i > 19) stop 'vdw radius not defined for Z > 19'
    vdw = radii(i) * aa2au

end function charge2vdw

!------------------------------------------------------------------------------

subroutine setup_solvent()

    integer :: i
    logical :: notfound = .true.
    character(len=6), dimension(3) :: solvents
    real(dp), dimension(3) :: epslist, epsinflist

    solvents = (/ 'H2O   ', 'CH3OH ', 'C2H5OH' /)

    epslist = (/ 78.39, 32.63, 24.55 /)

    epsinflist = (/ 1.776, 1.758, 1.847 /)

    do i = 1, 3
        if (trim(solvent) == trim(solvents(i))) then
            eps = epslist(i)
            epsinf = epsinflist(i)
            notfound = .false.
        end if
    end do

    if (notfound) stop 'ERROR: unknown solvent'

end subroutine setup_solvent

!------------------------------------------------------------------------------

subroutine setup_cavity()

     integer :: i, nz
     real(dp), dimension(:,:), allocatable :: all_coords
     real(dp), dimension(:), allocatable :: all_charges

     allocate(all_coords(3,qmnucs+nsites))
     allocate(all_charges(qmnucs+nsites))
     all_coords(:,1:qmnucs) = Rm
     all_charges(1:qmnucs) = Zm(1,:)
     nz = 0
     do i = 1, nsites
!         if (Zs(1,i) <= zero) then
!            nz = nz + 1
!            cycle
!         end if
         all_coords(:,qmnucs+i-nz) = Rs(:,i)
         all_charges(qmnucs+i-nz) = Zs(1,i)
     end do

    if (fixpva) then
       call fixtes(all_coords, all_charges)
    else
       call dalton_cavity(all_coords, all_charges, qmnucs + nsites - nz, work, size(work)) 
    end if

end subroutine setup_cavity

!------------------------------------------------------------------------------

subroutine chcase(string, uplo)

    character(len=*), intent(inout) :: string
    character(len=*), intent(in), optional :: uplo

    integer :: i, gap
    character(len=1) :: a, z, o_uplo

    if (present(uplo)) then
        o_uplo = uplo
    else
        o_uplo = 'u'
    end if

    gap = iachar('a') - iachar('A')

    if (o_uplo == 'u' .or. o_uplo == 'U') then
        a = 'a'
        z = 'z'
    else if (o_uplo == 'l' .or. o_uplo == 'L') then
        a = 'A'
        z = 'Z'
        gap = - gap
    else
        stop 'Unknown case specified'
    end if

    do i = 1, len_trim(string)
        if (lge(string(i:i), a) .and. lle(string(i:i), z)) then
            string(i:i) = achar(iachar(string(i:i)) - gap)
        end if
    end do

end subroutine chcase

!------------------------------------------------------------------------------

subroutine pe_save_density(denmat, mofckmat, cmo, nbas, nocc, norb, coords,&
                          & charges, dalwrk)

    external :: Tk_integrals

    integer, intent(in) :: nbas, nocc, norb
    real(dp), dimension(:), intent(in) :: denmat
    real(dp), dimension(:), intent(in) :: mofckmat
    real(dp), dimension(nbas,norb), intent(in) :: cmo
    real(dp), dimension(:), intent(in) :: charges
    real(dp), dimension(:,:), intent(in) :: coords
    real(dp), dimension(:), target, intent(inout) :: dalwrk

    integer :: i, j, l
    integer :: corenucs
    integer, parameter :: k = 0
    integer :: lucore, luden
    character(len=2) :: auoraa
    real(dp) :: Ene
    real(dp), dimension(:,:), allocatable :: Rc, Zc
    real(dp), dimension(:,:), allocatable :: full_denmat
    real(dp), dimension(:,:), allocatable :: T0_ints
    real(dp), dimension(:), allocatable :: Ffd
    real(dp), dimension(:,:), allocatable :: Ftmp
    real(dp), dimension(:), allocatable :: mo_energies

    work => dalwrk

    site_start = 1
    site_finish = nsites
    surp_start = 1
    surp_finish = nsurp

    ndens = 1
    nnbas = nbas * (nbas + 1) / 2

    ! fragment density nuclear charges and coordinates
!    qmnucs = size(charges)
!    allocate(Rm(3,qmnucs), Zm(1,qmnucs))
!    Rm = coords
!    Zm(1,:) = charges

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
    allocate(full_denmat(nbas,nbas)); full_denmat = 0.0d0
    l = 1
    do i = 1, nbas
        do j = 1, i
            if (j == i) then
                full_denmat(i,j) = denmat(l)
            else
                full_denmat(i,j) = 0.5d0 * denmat(l)
                full_denmat(j,i) = 0.5d0 * denmat(l)
            end if
            l = l + 1
        end do
    end do

    ! get electric field from fragment density at polarizable sites
    ! TODO: better solution for neglecting polarization
    ! TODO: potential from fragment density at surface points
    if (pe_polar) then
        allocate(Ftmp(3*npols,1), Ffd(3*npols)); Ftmp = 0.0d0
        call electron_fields(Ftmp, denmat)
        Ffd = Ftmp(:,1)
        call nuclear_fields(Ftmp(:,1))
        Ffd = Ffd + Ftmp(:,1)
        deallocate(Ftmp)
    else
        allocate(Ffd(3*npols))
        Ffd = 0.0d0
    end if

    ! calculate nuclear - electron energy contribution
    allocate(T0_ints(nnbas,1)); T0_ints = 0.0d0
    Ene = 0.0d0
    do i = 1, corenucs
        call Tk_integrals('es', T0_ints, nnbas, 1, Rc(:,i))
        T0_ints = Zc(1,i) * T0_ints
        Ene = Ene + dot(denmat, T0_ints(:,1))
    end do
    deallocate(T0_ints)

    allocate(mo_energies(nocc))
    do i = 1, nocc
        mo_energies(i) = mofckmat(i*(i+1)/2)
    end do

    ! save density, energy and field for subsequent calculations
    call openfile('pe_density.bin', luden, 'new', 'unformatted')
    rewind(luden)
    write(luden) Ene
    write(luden) qmnucs
    write(luden) Rm, Zm
    write(luden) npols
    write(luden) Ffd
    write(luden) nbas, nocc
    write(luden) full_denmat
    write(luden) cmo(:,1:nocc)
    write(luden) mo_energies
    close(luden)

    deallocate(full_denmat, Ffd, mo_energies)

end subroutine pe_save_density

!------------------------------------------------------------------------------

subroutine pe_twoints(nbas, nocc, norb, dalwrk)

    external :: sirfck, rdonel, dsptge

    integer, intent(in) :: nbas, nocc, norb
    real(dp), dimension(:), target, intent(inout) :: dalwrk

    integer :: i, j, k, l, m
    integer :: fbas, focc, cbas, cocc
    integer :: luden, lufck
    integer, dimension(1) :: isymdm, ifctyp
    real(dp) :: Ene
    real(dp), dimension(:), allocatable :: core_fckmat, Ffd
    real(dp), dimension(:,:), allocatable :: frag_denmat, full_denmat
    real(dp), dimension(:,:), allocatable :: full_fckmat
    real(dp), dimension(:), allocatable :: overlap, repmat
    real(dp), dimension(:,:), allocatable :: full_overlap
    real(dp), dimension(:,:), allocatable :: full_repmat
    real(dp), dimension(:), allocatable :: mo_energies
    real(dp), dimension(:,:), allocatable :: cmo, ecmo
    real(dp), dimension(:,:), allocatable :: weighted_denmat

    work => dalwrk

    site_start = 1
    site_finish = nsites
    surp_start = 1
    surp_finish = nsurp

    call openfile('pe_density.bin', luden, 'old', 'unformatted')
    rewind(luden)
    read(luden) Ene
    read(luden) fdnucs
    allocate(Rfd(3,fdnucs), Zfd(1,fdnucs))
    read(luden) Rfd, Zfd
    read(luden) npols
    allocate(Ffd(3*npols))
    read(luden) Ffd
    read(luden) fbas, focc
    allocate(frag_denmat(fbas, fbas))
    read(luden) frag_denmat
    allocate(cmo(fbas,focc))
    read(luden) cmo
    allocate(mo_energies(focc))
    read(luden) mo_energies
    close(luden)

    cbas = nbas - fbas
    cocc = nocc - focc

    ! full density matrix with fragment density in first block
    allocate(full_denmat(nbas,nbas))
    full_denmat = 0.0d0
    full_denmat(1:fbas,1:fbas) = frag_denmat
    deallocate(frag_denmat)

    ! get two-electron part of Fock matrix using resized density matrix
    allocate(full_fckmat(nbas,nbas))
    full_fckmat = 0.0d0
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
!     sirfck(fckmat, denmat, ?, isymdm, ifctyp, direct, work, nwrk)
    isymdm = 1
    ifctyp = 11
    call sirfck(full_fckmat, full_denmat, 1, isymdm, ifctyp, .true., work,&
               & size(work))
    deallocate(full_denmat)

    ! extract upper triangle part of full Fock matrix corresponding to
    ! core fragment
    allocate(core_fckmat(cbas*(cbas+1)/2))
    l = 1
    do j = fbas + 1, nbas
        do i = fbas + 1, j
            core_fckmat(l) = full_fckmat(i,j)
            l = l + 1
        end do
    end do

    deallocate(full_fckmat)

    ! Repulsion stuff from here
    allocate(overlap(nbas*(nbas+1)/2))
    allocate(full_overlap(nbas,nbas))
    call rdonel('OVERLAP', .true., overlap, nbas*(nbas+1)/2)
    call dsptge(nbas, overlap, full_overlap)
    deallocate(overlap)
!    call gemm(full_overlap(fbas+1:nbas,1:fbas),&
!             &full_overlap(1:fbas,fbas+1:nbas),&
!             &intmol_overlap)

    allocate(weighted_denmat(fbas,fbas))
    allocate(ecmo(fbas,focc))
    do i = 1, focc
        ecmo(:,i) = mo_energies(i) * cmo(:,i)
    end do
    weighted_denmat = - matmul(cmo, transpose(ecmo))
    deallocate(mo_energies, ecmo, cmo)

    allocate(full_repmat(cbas,cbas))
    full_repmat = matmul(matmul(full_overlap(fbas+1:nbas,1:fbas), weighted_denmat),&
                        & full_overlap(1:fbas,fbas+1:nbas))
    full_repmat = repfacs(1) * full_repmat + repfacs(2) * full_repmat**2&
                & + repfacs(3) * full_repmat**4
    deallocate(full_overlap, weighted_denmat)

    allocate(repmat(cbas*(cbas+1)/2))
    l = 1
    do j = 1, cbas
        do i = 1, j
            repmat(l) = full_repmat(i,j)
            l = l + 1
        end do
    end do
    deallocate(full_repmat)

    ! save core Fock matrix
    call openfile('pe_fock.bin', lufck, 'new', 'unformatted')
    rewind(lufck)
    write(lufck) Ffd
    write(lufck) Ene
    write(lufck) core_fckmat
    write(lufck) repmat
    write(lufck) fdnucs
    write(lufck) Rfd, Zfd
    close(lufck)

    deallocate(core_fckmat, repmat, Rfd, Zfd, Ffd)

end subroutine pe_twoints

!------------------------------------------------------------------------------

subroutine pe_mappot2points(o_coords)
    ! Calculates the electric potential and field at specific points.
    ! Default points are the positions of the qm nuclei.

    real(dp), dimension(:,:), intent(in), optional, target :: o_coords    

    character(len=1) :: tcmul
    character(len=99) :: cmul
    integer :: ncoords
    integer :: i, j, k, l
    integer :: lu
    real(dp) :: taylor, t_Vind, t_Vind_qmconv
    real(dp), dimension(3) :: Fs, Fs_qmconv, Rsp
    real(dp), dimension(:), allocatable :: t_Vpe, t_Fpe
    real(dp), dimension(:,:), allocatable :: Vpe, Vind, Vind_qmconv
    real(dp), dimension(:,:), allocatable :: Vtot, Ftot
    real(dp), dimension(:,:), allocatable :: Vtot_qmconv, Ftot_qmconv
    real(dp), dimension(:,:), allocatable :: Fpe, Fnrm, Fnrm_qmconv 
    real(dp), dimension(:,:), allocatable :: Fmuls, Find, Find_qmconv
    real(dp), dimension(:,:), allocatable :: M1inds, M1inds_qmconv
    real(dp), dimension(:), allocatable :: factors, Tsp
    real(dp), dimension(:,:), pointer :: coords
    
    if (present(o_coords)) then
        coords => o_coords
        ncoords = size(coords)/3
    else
        allocate(coords(3,qmnucs))
        ncoords = qmnucs         
        coords = Rm
    end if

    if (mulorder >= 0) then  
        allocate(Vpe(0:mulorder,1:ncoords))
        allocate(t_Vpe(0:mulorder))
        Vpe=0.0d0
        write(luout,*)
        write(luout,*) 'Potential from static multipoles'
        do i=1,ncoords   !positions on which pot and field are calculated.
            write(luout,'(a,2x,i4)') 'QM site no: ', i 
            do j = 1, nsites !MM sites 
                Rsp = coords(:,i) - Rs(:,j)
                t_Vpe = 0.0d0
                if (lmul(0)) then
                    allocate(Tsp(1),factors(1))
                    call symmetry_factors(factors)
                    taylor = 1.0d0 / factorial(0)
                    call Tk_tensor(Tsp, Rsp)
                    t_Vpe(0) = t_Vpe(0) + taylor * factors(1) * Tsp(1) *&
                             & M0s(1,j)
                    Vpe(0,i) = Vpe(0,i) + t_Vpe(0)
                    deallocate(Tsp, factors)
                end if
                if (lmul(1)) then
                    allocate(Tsp(3), factors(3))
                    call symmetry_factors(factors)
                    taylor = - 1.0d0 / factorial(1)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 3
                        t_Vpe(1) = t_Vpe(1) + taylor * factors(k) * Tsp(k) *&
                                 & M1s(k,j)
                    end do
                    Vpe(1,i) = Vpe(1,i) + t_Vpe(1)
                    deallocate(Tsp, factors)
                end if
                if (lmul(2)) then
                    allocate(Tsp(6), factors(6))
                    call symmetry_factors(factors)
                    taylor = 1.0d0 / factorial(2)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 6
                        t_Vpe(2) = t_Vpe(2) + taylor * factors(k) * Tsp(k) *&
                                 & M2s(k,j)
                    end do
                    Vpe(2,i) = Vpe(2,i)+t_Vpe(2)
                    deallocate(Tsp, factors)
                end if
                if (lmul(3)) then
                    allocate(Tsp(10), factors(10))
                    call symmetry_factors(factors)
                    taylor = - 1.0d0 / factorial(3)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 10
                        t_Vpe(3) = t_Vpe(3) + taylor * factors(k) * Tsp(k) *&
                                 & M3s(k,j)
                    end do
                    Vpe(3,i) = Vpe(3,i) + t_Vpe(3)
                    deallocate(Tsp, factors)
                end if
                if (lmul(4)) then
                    allocate(Tsp(15), factors(15))
                    call symmetry_factors(factors)
                    taylor = 1.0d0 / factorial(4)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 15
                        t_Vpe(4) = t_Vpe(4) + taylor * factors(k) * Tsp(k) *&
                                 & M4s(k,j)
                    end do
                    Vpe(4,i) = Vpe(4,i) + t_Vpe(4)
                    deallocate(Tsp, factors)
                end if
                if (lmul(5)) then
                    allocate(Tsp(21), factors(21))
                    call symmetry_factors(factors)
                    taylor = - 1.0d0 / factorial(5)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 21
                        t_Vpe(5) = t_Vpe(5) + taylor * factors(k) * Tsp(k) *&
                                 & M5s(k,j)
                    end do
                    Vpe(5,i) = Vpe(5,i) + t_Vpe(5)
                    deallocate(Tsp, factors)
                end if
                write(luout,'(a4,2x,i6,(e13.5))') 'Site',j, sum(t_Vpe)
            end do
        end do
        deallocate(t_Vpe)
    end if

    if (lpol(1)) then 
        !induced dipoles when QM present
        allocate(M1inds_qmconv(3*npols,1))
        call openfile('pe_induced_dipoles.bin', lu, 'old', 'unformatted')
        rewind(lu)
        read(lu) M1inds_qmconv 
!        write(luout,*) 'Now qmconv induced dipoles have been read'
        close(lu)

        !calculate induced dipoles when QM absent
        allocate(Fmuls(3*npols,1))
        allocate(M1inds(3*npols,1))
        call multipole_fields(Fmuls(:,1))
        call induced_moments(M1inds, Fmuls)
 !       write(luout,*) 'After call induced_dipoles'
        deallocate(Fmuls)         

        allocate(Vind(1,ncoords))
        allocate(Vind_qmconv(1,ncoords))
        allocate(Tsp(3), factors(3))
        Vind = 0.0d0
        Vind_qmconv = 0.0d0
        call symmetry_factors(factors)
        taylor = -1.0d0 / factorial(1)

        write(luout,*)
        write(luout,*) 'Induced potential from qm_conv (au)' 
        do i=1, ncoords
            write(luout,*) 'QM site no: ', i
            l = 0
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                Rsp = coords(:,i) - Rs(:,j)
                Tsp = 0.0d0
                t_Vind = 0.0d0
                t_Vind_qmconv = 0.0d0
                call Tk_tensor(Tsp, Rsp)
              !  write(luout,*) 'After call Tk_tensor'
                do k = 1, 3
                    t_Vind = t_Vind + taylor * factors(k) * Tsp(k) *&
                           & M1inds(l+k,1)
                    t_Vind_qmconv = t_Vind_qmconv + taylor * factors(k) *&
                                  & Tsp(k) * M1inds_qmconv(l+k,1)
                !write(luout,*) t_Vind_qmconv
                end do 
                Vind(1,i) = Vind(1,i) + t_Vind 
                Vind_qmconv(1,i) = Vind_qmconv(1,i) + t_Vind_qmconv
                write(luout,'(a4,2x,i6,(e13.5))') 'Site',j, t_Vind_qmconv
                l = l + 3
            end do 
        end do
        deallocate(Tsp, factors)
        !deallocate(M1inds, M1inds_qmconv)
!        write(luout,*) 'We have now computed Vind'
    end if

    if (mulorder >= 0) then
        write(luout,*) 
        write(luout,*) 'Field from static multipoles (x, y, z)'
        allocate(Fpe(3*(mulorder+1),ncoords)) 
        Fpe = 0.0d0
        do i = 1, ncoords
            write(luout,*) 'QM site no: ', i
            do j = 1, nsites
                Rsp = coords(:,i) - Rs(:,j)
                allocate(t_Fpe(3))
                t_Fpe = 0.0d0
                if (lmul(0)) then
                    Fs = 0.0d0
                    call multipole_field(Fs, Rsp, M0s(:,j))
                    Fpe(1:3,i) = Fpe(1:3,i) + Fs
                    t_Fpe(1:3) = Fs
                end if
                if (lmul(1)) then
                    Fs = 0.0d0
                    call multipole_field(Fs, Rsp, M1s(:,j))
                    Fpe(4:6,i) = Fpe(4:6,i) + Fs
                    t_Fpe(1:3) = t_Fpe(1:3) + Fs
                end if
                if (lmul(2)) then
                    Fs = 0.0d0
                    call multipole_field(Fs, Rsp, M2s(:,j))
                    Fpe(7:9,i) = Fpe(7:9,i) + Fs
                    t_Fpe(1:3) = t_Fpe(1:3) + Fs
                end if
                if (lmul(3)) then
                    Fs = 0.0d0
                    call multipole_field(Fs, Rsp, M3s(:,j))
                    Fpe(10:12,i) = Fpe(10:12,i) + Fs
                    t_Fpe(1:3) = t_Fpe(1:3) + Fs
                end if
                if (lmul(4)) then
                    Fs = 0.0d0
                    call multipole_field(Fs, Rsp, M4s(:,j))
                    Fpe(13:15,i) = Fpe(13:15,i) + Fs
                    t_Fpe(1:3) = t_Fpe(1:3) + Fs
                    end if
                if (lmul(5)) then
                    Fs = 0.0d0
                    call multipole_field(Fs, Rsp, M5s(:,j))
                    Fpe(16:18,i) = Fpe(16:18,i) + Fs
                    t_Fpe(1:3) = t_Fpe(1:3) + Fs
                end if
!                l = 1
!                do k = 1,mulorder+1
!                    t_Fpe(1) = t_Fpe(1) + Fpe(l,i)
!                    t_Fpe(2) = t_Fpe(2) + Fpe(l+1,i)
!                    t_Fpe(3) = t_Fpe(3) + Fpe(l+2,i)
!                    l = l + 3
!                end do  
                write(luout,'(a4,2x,i6,(3e13.5))') 'Site',j, t_Fpe
                deallocate(t_Fpe)
            end do
        end do
    end if

    if (lpol(1)) then 
        write(luout,*)
        write(luout,*) 'Induced field from qm_conv (x, y, z)'
        allocate(Find(3,ncoords))
        allocate(Find_qmconv(3,ncoords))
        Find = 0.0d0
        Find_qmconv = 0.0d0

        do i = 1, ncoords
            write(luout,*) 'QM site no: ', i
            l = 0
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                Rsp = coords(:,i) - Rs(:,j)
                Fs = 0.0d0     
                Fs_qmconv = 0.0d0
                call multipole_field(Fs, Rsp, M1inds(l:l+2,1))
                call multipole_field(Fs_qmconv, Rsp, M1inds_qmconv(l:l+2,1))
                Find(:,i) = Find(:,i) + Fs 
                Find_qmconv(:,i) = Find_qmconv(:,i) + Fs_qmconv
                write(luout,'(a4,2x,i6,(3e13.5))') 'Site',j, Fs_qmconv
                l = l + 3
            end do
        end do
    end if

! Calculate total multipole potential and field
    if (mulorder >= 0) then
        allocate(Vtot(1,ncoords))
        allocate(Ftot(3,ncoords))
        allocate(Fnrm(1,ncoords))
        Vtot = 0.0d0
        Ftot = 0.0d0
        Fnrm = 0.0d0
!        if (lpol(1)) then
!            allocate(Vtot_qmconv(1,ncoords))
!            allocate(Ftot_qmconv(3,ncoords))
!            allocate(Fnrm_qmconv(1,ncoords))
!            Vtot_qmconv = 0.0d0
!            Ftot_qmconv = 0.0d0
!            Fnrm_qmconv = 0.0d0
!        end if
        do i = 1, ncoords
            if (mulorder == 0) then
                Vtot(1,i) = Vpe(0,i)
                Ftot(1:3,i) = Fpe(1:3,i)

            else if (mulorder == 1) then
                Vtot(1,i) = Vpe(0,i) + Vpe(1,i)
                Ftot(1:3,i) = Fpe(1:3,i) + Fpe(4:6,i)

            else if (mulorder == 2) then
                Vtot(1,i) = Vpe(0,i) + Vpe(1,i) + Vpe(2,i)
                Ftot(1:3,i) = Fpe(1:3,i) + Fpe(4:6,i) + Fpe(7:9,i)

            else if (mulorder == 3) then
                Vtot(1,i) = Vpe(0,i) + Vpe(1,i) + Vpe(2,i) + Vpe(3,i)
                Ftot(1:3,i) = Fpe(1:3,i) + Fpe(4:6,i) + Fpe(7:9,i) +&
                            & Fpe(10:12,i)

            else if (mulorder == 4) then
                Vtot(1,i) = Vpe(0,i) + Vpe(1,i) + Vpe(2,i) + Vpe(3,i) +&
                          & Vpe(4,i)
                Ftot(1:3,i) = Fpe(1:3,i) + Fpe(4:6,i) + Fpe(7:9,i) +&
                            & Fpe(10:12,i) + Fpe(13:15,i)

            else if (mulorder == 5) then
                Vtot(1,i) = Vpe(0,1) + Vpe(1,i) + Vpe(2,i) + Vpe(3,i) +&
                          & Vpe(4,i) + Vpe(5,i)
                Ftot(1:3,i) = Fpe(1:3,i) + Fpe(4:6,i) + Fpe(7:9,i) +&
                            & Fpe(10:12,i) + Fpe(13:15,i) + Fpe(16:18,i)
            end if

            ! add contribution from induced dipoles to potential and field
            ! if present
!            if (lpol(1)) then 
!                ! the order of E/Vtot_qmconv and E/Vtot IS important because
!                ! E/Vtot is changed
!                Vtot_qmconv(1,i) = Vtot(1,i) + Vind_qmconv(1,i)
!                Vtot(1,i) = Vtot(1,i) + Vind(1,i)
!                Ftot_qmconv(1:3,i) = Ftot(1:3,i) + Find_qmconv(1:3,i)
!                Ftot(1:3,i) = Ftot(1:3,i) + Find(1:3,i)

            !calculate norm of Efield in a given site i
!                Fnrm_qmconv(1,i) = nrm2(Ftot_qmconv(1:3,i))
!                !calculate only if Fnrm_qmconv present
            Fnrm(1,i) = nrm2(Ftot(1:3,i))    
!            end if
        end do
    end if

    if (mulorder>=0) then
        write(luout,'(//a)') repeat('*', 54)
        write(luout,'(a)')  '*** Internal electric potential and field&
                            & analysis ***'
        write(luout,'(a//)') repeat('*', 54)

        if (lpol(1)) then
            write(luout,'(a)') repeat('-', 72)
            write(luout,'(20x,a,24x)') 'QM present', 'QM absent'
            write(luout,'(x,a,6x,a,9x,a,11x,a,8x,a,9x,a)')  'Site', 'Vtot',&
                                                            & 'Vpe', 'Vind',&
                                                            & 'Vtot', 'Vind'
            write(luout,'(a)') repeat('=', 72)
            do i = 1, ncoords
                write(luout,'(i4,2x,(5e13.5))') i, (Vtot(1,i) +&
                                                & Vind_qmconv(1,i)),&
                                                & Vtot(1,i),&
                                                & Vind_qmconv(1,i),&
                                                & (Vtot(1,i) + Vind(1,i)),&
                                                & Vind(1,i)       
            end do
            write(luout,'(a)') repeat('=', 72)
            write(luout,'(/a)') repeat('-', 72)
            write(luout,'(8x,a)') 'QM present'
            write(luout,'(x,a,6x,a,23x,a,34x,a)') 'Site', 'Enrm', 'Epe',&
                                                  & 'Eind'
            write(luout,'(25x,a,12x,a,12x,a,12x,a,12x,a,12x,a)') 'x', 'y',&
                                                                 & 'z', 'x',&
                                                                 & 'y', 'z'
            write(luout,'(a)') repeat('=', 72)
            do i = 1, ncoords
                write(luout,'(i4,2x,(7e13.5))') i, (nrm2(Ftot(1:3,i) +&
                                                & Find_qmconv(1:3,i))),&
                                                & Ftot(1:3,i),&
                                                & Find_qmconv(1:3,i) 
            end do
            write(luout,'(a/)') repeat('=', 72)
            write(luout,'(a)') repeat('-', 72)
            write(luout,'(8x,a)') 'QM absent'
            write(luout,'(x,a,6x,a,23x,a,34x,a)') 'Site', 'Enrm', 'Epe',&
                                                  & 'Eind'
            write(luout,'(25x,a,12x,a,12x,a,12x,a,12x,a,12x,a)') 'x', 'y',& 
                                                                 & 'z', 'x',& 
                                                                 & 'y', 'z'
            write(luout,'(a)') repeat('=', 72)
            do i = 1, ncoords
                write(luout,'(i4,2x,(7e13.5))') i, (nrm2(Ftot(1:3,i) +&
                                                & Find(1:3,i))), Ftot(1:3,i),&
                                                & Find(1:3,i)
            end do
            write(luout,'(a)') repeat('=', 72)
        else 
            write(luout,'(a)') repeat('-', 72)
            write(luout,'(x,a,6x,a,23x,a,21x,a)') 'Site', 'Vtot', 'Etot', 'Enrm'
            write(luout,'(25x,a,14x,a,14x,a)')  'x', 'y', 'z'
            write(luout,'(a)') repeat('=', 72)
            do i = 1, ncoords
                write(luout,'(i4,2x,(5e13.5))') i, Vtot(1,i), Ftot(1:3,i),&
                                                & nrm2(Ftot(1:3,i))
            end do
            write(luout,'(a/)') repeat('=', 72)
        end if

        write(luout,'(/a/)') 'Contributions from static multipoles MXs'

        do j = 1, mulorder + 1
            write(cmul,*) j-1
            tcmul = trim(adjustl(cmul))
            write(luout,'(2a)') 'M'//tcmul
            write(luout,'(a)') repeat('-', 72)
            write(luout,'(x,a,6x,a,25x)')  'Site', 'Vpe', 'E' 
            write(luout,'(25x,a,14x,a,14x,a)')  'x', 'y', 'z'
            write(luout,'(a)') repeat('=', 72)
            do i = 1, ncoords
                write(luout,'(i4,2x,(4e13.5))') i, Vpe(j-1,i), Fpe(j:j+2,i)
            end do

            write(luout,'(a/)') repeat('=', 72)
        end do    
    end if

end subroutine pe_mappot2points

!------------------------------------------------------------------------------

subroutine pe_compute_mep(denmats)

    real(dp), dimension(:), intent(in) :: denmats

    character(len=1) :: tcl
    character(len=99) :: cl
    integer :: point
    integer :: i, j, k, l
    integer :: ndist, nrest
    integer :: lu, lum0, lum1, lum2, lum3, lum4, lum5
    logical :: exclude
    real(dp) :: taylor
    real(dp), dimension(3) :: Tm, Rsp, Rji, Fs
    real(dp), dimension(:,:), allocatable :: Vqm, Vpe, Vind
    real(dp), dimension(:,:), allocatable :: Fqm, Find
    real(dp), dimension(:,:,:), allocatable :: Fpe
    real(dp), dimension(:,:), allocatable :: Fmuls, M1inds
    real(dp), dimension(:,:), allocatable :: Tk_ints
    real(dp), dimension(:), allocatable :: factors, Tsp

    if (mep_qmcube) then
        allocate(Vqm(1,npoints))
        allocate(Tk_ints(nnbas,1))
        i = 1
        do point = site_start, site_finish
            call Tk_integrals('es', Tk_ints(:,1), nnbas, 1, mepgrid(:,i))
            Vqm(1,i) = dot(denmats, Tk_ints(:,1))
            do j = 1, qmnucs
                call Tk_tensor(Tm(1:1), mepgrid(:,i) - Rm(:,j))
                Vqm(1,i) = Vqm(1,i) + Zm(1,j) * Tm(1)
            end do
            i = i + 1
        end do
        deallocate(Tk_ints)
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + mepdists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, 0, rmpi,&
                            &Vqm, mepdists, displs, rmpi,&
                            &0, comm, ierr)
        else if (myid /= 0) then
            call mpi_gatherv(Vqm, mepdists(myid), rmpi,&
                            &0, 0, 0, rmpi,&
                            &0, comm, ierr)
        end if
#endif
        if (myid == 0) then
            call openfile('qm_mep.cube', lu, 'new', 'formatted')
            write(lu,'(a)') 'QM MEP'
            write(lu,'(a)') 'Generated by the polarizable embedding (PE) module'
            write(lu,'(i5,3f12.6)') qmnucs, origin
            write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0d0, 0.0d0
            write(lu,'(i5,3f12.6)') ysteps, 0.0d0, step(2), 0.0d0
            write(lu,'(i5,3f12.6)') zsteps, 0.0d0, 0.0d0, step(3)
            do j = 1, qmnucs
                write(lu,'(i5,4f12.6)') int(Zm(1,j)), Zm(1,j), Rm(:,j)
            end do
            do i = 1, xsteps * ysteps
                j = (i - 1) * zsteps + 1
                k = j - 1 + zsteps
                write(lu,'(6e13.5)') Vqm(1,j:k)
            end do
            close(lu)
        end if
        deallocate(Vqm)
    end if

    if (mulorder >= 0 .and. mep_mulcube) then
        allocate(Vpe(0:mulorder,npoints))
        Vpe = 0.0d0
        i = 1
        do point = site_start, site_finish
            do j = 1, nsites
                Rsp = mepgrid(:,i) - Rs(:,j)
                if (lmul(0)) then
                    allocate(Tsp(1), factors(1))
                    call symmetry_factors(factors)
                    taylor = 1.0d0 / factorial(0)
                    call Tk_tensor(Tsp, Rsp)
                    Vpe(0,i) = Vpe(0,i) + taylor * factors(1) * Tsp(1) *&
                             & M0s(1,j)
                    deallocate(Tsp, factors)
                end if
                if (lmul(1)) then
                    allocate(Tsp(3), factors(3))
                    call symmetry_factors(factors)
                    taylor = - 1.0d0 / factorial(1)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 3
                        Vpe(1,i) = Vpe(1,i) + taylor * factors(k) * Tsp(k) *&
                                 & M1s(k,j)
                    end do
                    deallocate(Tsp, factors)
                end if
                if (lmul(2)) then
                    allocate(Tsp(6), factors(6))
                    call symmetry_factors(factors)
                    taylor = 1.0d0 / factorial(2)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 6
                        Vpe(2,i) = Vpe(2,i) + taylor * factors(k) * Tsp(k) *&
                                 & M2s(k,j)
                    end do
                    deallocate(Tsp, factors)
                end if
                if (lmul(3)) then
                    allocate(Tsp(10), factors(10))
                    call symmetry_factors(factors)
                    taylor = - 1.0d0 / factorial(3)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 10
                        Vpe(3,i) = Vpe(3,i) + taylor * factors(k) * Tsp(k) *&
                                 & M3s(k,j)
                    end do
                    deallocate(Tsp, factors)
                end if
                if (lmul(4)) then
                    allocate(Tsp(15), factors(15))
                    call symmetry_factors(factors)
                    taylor = 1.0d0 / factorial(4)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 15
                        Vpe(4,i) = Vpe(4,i) + taylor * factors(k) * Tsp(k) *&
                                 & M4s(k,j)
                    end do
                    deallocate(Tsp, factors)
                end if
                if (lmul(5)) then
                    allocate(Tsp(21), factors(21))
                    call symmetry_factors(factors)
                    taylor = - 1.0d0 / factorial(5)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 21
                        Vpe(5,i) = Vpe(5,i) + taylor * factors(k) * Tsp(k) *&
                                 & M5s(k,j)
                    end do
                    deallocate(Tsp, factors)
                end if
            end do
            i = i + 1
        end do
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + (mulorder + 1) * mepdists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, 0, rmpi,&
                            &Vpe, (mulorder+1)*mepdists, displs, rmpi,&
                            &0, comm, ierr)
        else if (myid /= 0) then
            call mpi_gatherv(Vpe, (mulorder+1)*mepdists(myid), rmpi,&
                            &0, 0, 0, rmpi,&
                            &0, comm, ierr)
        end if
#endif
        if (myid == 0) then
            if (mulorder >= 0) then
                call openfile('m0_mep.cube', lum0, 'new', 'formatted')
            end if
            if (mulorder >= 1) then
                call openfile('m1_mep.cube', lum1, 'new', 'formatted')
            end if
            if (mulorder >= 2) then
                call openfile('m2_mep.cube', lum2, 'new', 'formatted')
            end if
            if (mulorder >= 3) then
                call openfile('m3_mep.cube', lum3, 'new', 'formatted')
            end if
            if (mulorder >= 4) then
                call openfile('m4_mep.cube', lum4, 'new', 'formatted')
            end if
            if (mulorder >= 5) then
                call openfile('m5_mep.cube', lum5, 'new', 'formatted')
            end if
            do i = 1, mulorder + 1
                if (i == 1) then
                    lu = lum0
                else if (i == 2) then
                    lu = lum1
                else if (i == 3) then
                    lu = lum2
                else if (i == 4) then
                    lu = lum3
                else if (i == 5) then
                    lu = lum4
                else if (i == 6) then
                    lu = lum5
                end if
                write(lu,'(a)') 'PE electrostatic potential'
                write(lu,'(a)') 'Generated by the polarizable embedding (PE) module'
                write(lu,'(i5,3f12.6)') qmnucs, origin
                write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0d0, 0.0d0
                write(lu,'(i5,3f12.6)') ysteps, 0.0d0, step(2), 0.0d0
                write(lu,'(i5,3f12.6)') zsteps, 0.0d0, 0.0d0, step(3)
                do j = 1, qmnucs
                    write(lu,'(i5,4f12.6)') int(Zm(1,j)), Zm(1,j), Rm(:,j)
                end do
            end do
            do i = 1, xsteps * ysteps
                j = (i - 1) * zsteps + 1
                k = j - 1 + zsteps
                if (mulorder >= 0) then
                    write(lum0,'(6e13.5)') Vpe(0,j:k)
                end if
                if (mulorder >= 1) then
                    write(lum1,'(6e13.5)') (Vpe(0,j:k) + Vpe(1,j:k))
                end if
                if (mulorder >= 2) then
                    write(lum2,'(6e13.5)') (Vpe(0,j:k) + Vpe(1,j:k) +&
                                           & Vpe(2,j:k))
                end if
                if (mulorder >= 3) then
                    write(lum3,'(6e13.5)') (Vpe(0,j:k) + Vpe(1,j:k) +&
                                           & Vpe(2,j:k) + Vpe(3,j:k))
                end if
                if (mulorder >= 4) then
                    write(lum4,'(6e13.5)') (Vpe(0,j:k) + Vpe(1,j:k) +&
                                           & Vpe(2,j:k) + Vpe(3,j:k) +&
                                           & Vpe(4,j:k))
                end if
                if (mulorder >= 5) then
                    write(lum5,'(6e13.5)') (Vpe(0,j:k) + Vpe(1,j:k) +&
                                           & Vpe(2,j:k) + Vpe(3,j:k) +&
                                           & Vpe(4,j:k) + Vpe(5,j:k))
                end if
            end do
            if (mulorder >= 0) then
                close(lum0)
            end if
            if (mulorder >= 1) then
                close(lum1)
            end if
            if (mulorder >= 2) then
                close(lum2)
            end if
            if (mulorder >= 3) then
                close(lum3)
            end if
            if (mulorder >= 4) then
                close(lum4)
            end if
            if (mulorder >= 5) then
                close(lum5)
            end if
        endif
        deallocate(Vpe)
    end if

    if (lpol(1)) then
        allocate(Fmuls(3*npols,1))
        Fmuls = 0.0d0
        l = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            do j = 1, nsites
                if (i == j) then
                    cycle
                end if
                exclude = .false.
                do k = 1, lexlst
                    if (exclists(k,i) == exclists(1,j)) then
                        exclude = .true.
                        exit
                    end if
                end do
                if (exclude) cycle
                Rji = Rs(:,i) - Rs(:,j)
                if (lmul(0)) then
                    if (abs(maxval(M0s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2,1), Rji, M0s(:,j))
                    end if
                end if
                if (lmul(1)) then
                    if (abs(maxval(M1s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2,1), Rji, M1s(:,j))
                    end if
                end if
                if (lmul(2)) then
                    if (abs(maxval(M2s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2,1), Rji, M2s(:,j))
                    end if
                end if
                if (lmul(3)) then
                    if (abs(maxval(M3s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2,1), Rji, M3s(:,j))
                    end if
                end if
                if (lmul(4)) then
                    if (abs(maxval(M4s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2,1), Rji, M4s(:,j))
                    end if
                end if
                if (lmul(5)) then
                    if (abs(maxval(M5s(:,j))) >= zero) then
                        call multipole_field(Fmuls(l:l+2,1), Rji, M5s(:,j))
                    end if
                end if
            end do
            l = l + 3
        end do
        allocate(M1inds(3*npols,1))
!        call multipole_fields(Fmuls(:,1))
        if (myid == 0) then
            write(luout,*) 'Electric fields from static multipole moments:'
            write(luout,'(3f12.6)') Fmuls
        end if
        if (mep_extfld) then
            j = 1
            do i = 1, npols
                Fmuls(j:j+2,1) = Fmuls(j:j+2,1) + extfld
                j = j + 3
            end do
        end if
        call induced_moments(M1inds, Fmuls)
        if (myid == 0) then
            write(luout,*) 'Induced dipole moments:'
            write(luout,'(3f12.6)') M1inds
        end if
        deallocate(Fmuls)
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(M1inds, 3*npols, rmpi,&
                          &0, comm, ierr)
        end if
#endif
        allocate(Vind(1,npoints))
        Vind = 0.0d0
        allocate(Tsp(3), factors(3))
        call symmetry_factors(factors)
        taylor = - 1.0d0 / factorial(1)
        i = 1
        do point = site_start, site_finish
            l = 0
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                Rsp = mepgrid(:,i) - Rs(:,j)
                Tsp = 0.0d0
                call Tk_tensor(Tsp, Rsp)
                do k = 1, 3
                    Vind(1,i) = Vind(1,i) + taylor * factors(k) * Tsp(k) *&
                              & M1inds(l+k,1)
                end do
                l = l + 3
            end do
            i = i + 1
        end do
        deallocate(Tsp, factors)
        deallocate(M1inds)
#if defined(VAR_MPI)
        if (myid == 0 .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + mepdists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, 0, rmpi,&
                            &Vind, mepdists, displs, rmpi,&
                            &0, comm, ierr)
        else if (myid /= 0) then
            call mpi_gatherv(Vind, mepdists(myid), rmpi,&
                            &0, 0, 0, rmpi,&
                            &0, comm, ierr)
        end if
#endif
        if (myid == 0) then
            call openfile('ind_mep.cube', lu, 'new', 'formatted')
            write(lu,'(a)') 'PE induced potential'
            write(lu,'(a)') 'Generated by the polarizable embedding (PE) module'
            write(lu,'(i5,3f12.6)') qmnucs, origin
            write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0d0, 0.0d0
            write(lu,'(i5,3f12.6)') ysteps, 0.0d0, step(2), 0.0d0
            write(lu,'(i5,3f12.6)') zsteps, 0.0d0, 0.0d0, step(3)
            do j = 1, qmnucs
                write(lu,'(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j), Rm(:,j)
            end do
            do i = 1, xsteps * ysteps
                j = (i - 1) * zsteps + 1
                k = j - 1 + zsteps
                write(lu,'(6e13.5)') Vind(1,j:k)
            end do
            close(lu)
        end if
        deallocate(Vind)
    end if

    if (mep_field) then
        if (mep_qmcube) then
            allocate(Fqm(3,npoints))
            allocate(Tk_ints(nnbas,3))
            i = 1
            do point = site_start, site_finish
                call Tk_integrals('es', Tk_ints, nnbas, 3, mepgrid(:,i))
                do j = 1, 3
                    Fqm(j,i) = dot(denmats, Tk_ints(:,j))
                end do
                do j = 1, qmnucs
                    call Tk_tensor(Tm, mepgrid(:,i) - Rm(:,j))
                    do k = 1, 3
                        Fqm(k,i) = Fqm(k,i) - Zm(1,j) * Tm(k)
                    end do
                end do
                i = i + 1
            end do
            deallocate(Tk_ints)
#if defined(VAR_MPI)
            if (myid == 0 .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + 3 * mepdists(i-1)
                end do
                call mpi_gatherv(mpi_in_place, 0, rmpi,&
                                &Fqm, 3*mepdists, displs, rmpi,&
                                &0, comm, ierr)
            else if (myid /= 0) then
                call mpi_gatherv(Fqm, 3*mepdists(myid), rmpi,&
                                &0, 0, 0, rmpi,&
                                &0, comm, ierr)
            end if
#endif
            if (myid == 0) then
                if (mep_fldnrm) then
                    call openfile('qm_field.cube', lu, 'new', 'formatted')
                    write(lu,'(a)') 'QM electric field norm'
                    write(lu,'(a)') 'Generated by the polarizable embedding&
                                    & (PE) module'
                    write(lu,'(i5,3f12.6)') qmnucs, origin
                    write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0d0, 0.0d0
                    write(lu,'(i5,3f12.6)') ysteps, 0.0d0, step(2), 0.0d0
                    write(lu,'(i5,3f12.6)') zsteps, 0.0d0, 0.0d0, step(3)
                    do j = 1, qmnucs
                        write(lu,'(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j), Rm(:,j)
                    end do
                    do i = 1, xsteps * ysteps
                        j = (i - 1) * zsteps + 1
                        k = j - 1 + zsteps
                        write(lu,'(6e13.5)') (nrm2(Fqm(:,l)), l = j, k)
                    end do
                    close(lu)
                else
                    do l = 1, 3
                        write(cl,*) l
                        tcl = trim(adjustl(cl))
                        call openfile('qm_field_'//tcl//'.cube', lu, 'new',&
                                     & 'formatted')
                        write(lu,'(a)') 'QM electric field component '//tcl
                        write(lu,'(a)') 'Generated by the Polarizable&
                                        & Embedding module'
                        write(lu,'(i5,3f12.6)') qmnucs, origin
                        write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0d0, 0.0d0
                        write(lu,'(i5,3f12.6)') ysteps, 0.0d0, step(2), 0.0d0
                        write(lu,'(i5,3f12.6)') zsteps, 0.0d0, 0.0d0, step(3)
                        do j = 1, qmnucs
                            write(lu,'(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j),&
                                                    & Rm(:,j)
                        end do
                        do i = 1, xsteps * ysteps
                            j = (i - 1) * zsteps + 1
                            k = j - 1 + zsteps
                            write(lu,'(6e13.5)') Fqm(l,j:k)
                        end do
                        close(lu)
                    end do
                end if
            end if
            deallocate(Fqm)
        end if

        if (mulorder >= 0 .and. mep_mulcube) then
            allocate(Fpe(3,0:mulorder,npoints))
            Fpe = 0.0d0
            i = 1
            do point = site_start, site_finish
                do j = 1, nsites
                    Rsp = mepgrid(:,i) - Rs(:,j)
                    if (lmul(0)) then
                        Fs = 0.0d0
                        call multipole_field(Fs, Rsp, M0s(:,j))
                        Fpe(:,0,i) = Fpe(:,0,i) + Fs
                    end if
                    if (lmul(1)) then
                        Fs = 0.0d0
                        call multipole_field(Fs, Rsp, M1s(:,j))
                        Fpe(:,1,i) = Fpe(:,1,i) + Fs
                    end if
                    if (lmul(2)) then
                        Fs = 0.0d0
                        call multipole_field(Fs, Rsp, M2s(:,j))
                        Fpe(:,2,i) = Fpe(:,2,i) + Fs
                    end if
                    if (lmul(3)) then
                        Fs = 0.0d0
                        call multipole_field(Fs, Rsp, M3s(:,j))
                        Fpe(:,3,i) = Fpe(:,3,i) + Fs
                    end if
                    if (lmul(4)) then
                        Fs = 0.0d0
                        call multipole_field(Fs, Rsp, M4s(:,j))
                        Fpe(:,4,i) = Fpe(:,4,i) + Fs
                    end if
                    if (lmul(5)) then
                        Fs = 0.0d0
                        call multipole_field(Fs, Rsp, M5s(:,j))
                        Fpe(:,5,i) = Fpe(:,5,i) + Fs
                    end if
                end do
                i = i + 1
            end do
#if defined(VAR_MPI)
            if (myid == 0 .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + 3 * (mulorder + 1) * mepdists(i-1)
                end do
                call mpi_gatherv(mpi_in_place, 0, rmpi,&
                                &Fpe, 3*(mulorder+1)*mepdists, displs, rmpi,&
                                &0, comm, ierr)
            else if (myid /= 0) then
                call mpi_gatherv(Fpe, 3*(mulorder+1)*mepdists(myid), rmpi,&
                                &0, 0, 0, rmpi,&
                                &0, comm, ierr)
            end if
#endif
            if (myid == 0) then
                if (mep_fldnrm) then
                    if (mulorder >= 0) then
                        call openfile('m0_field.cube', lum0, 'new', 'formatted')
                    end if
                    if (mulorder >= 1) then
                        call openfile('m1_field.cube', lum1, 'new', 'formatted')
                    end if
                    if (mulorder >= 2) then
                        call openfile('m2_field.cube', lum2, 'new', 'formatted')
                    end if
                    if (mulorder >= 3) then
                        call openfile('m3_field.cube', lum3, 'new', 'formatted')
                    end if
                    if (mulorder >= 4) then
                        call openfile('m4_field.cube', lum4, 'new', 'formatted')
                    end if
                    if (mulorder >= 5) then
                        call openfile('m5_field.cube', lum5, 'new', 'formatted')
                    end if
                    do i = 1, mulorder + 1
                        if (i == 1) then
                            lu = lum0
                        else if (i == 2) then
                            lu = lum1
                        else if (i == 3) then
                            lu = lum2
                        else if (i == 4) then
                            lu = lum3
                        else if (i == 5) then
                            lu = lum4
                        else if (i == 6) then
                            lu = lum5
                        end if
                        write(lu,'(a)') 'PE electric field norm'
                        write(lu,'(a)') 'Generated by the Polarizable&
                                        & Embedding module'
                        write(lu,'(i5,3f12.6)') qmnucs, origin
                        write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0d0, 0.0d0
                        write(lu,'(i5,3f12.6)') ysteps, 0.0d0, step(2), 0.0d0
                        write(lu,'(i5,3f12.6)') zsteps, 0.0d0, 0.0d0, step(3)
                        do j = 1, qmnucs
                            write(lu,'(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j),&
                                                    & Rm(:,j)
                        end do
                    end do
                    do i = 1, xsteps * ysteps
                        j = (i - 1) * zsteps + 1
                        k = j - 1 + zsteps
                        if (mulorder >= 0) then
                            write(lum0,'(6e13.5)') (nrm2(Fpe(:,0,l)), l = j, k)
                        end if
                        if (mulorder >= 1) then
                            write(lum1,'(6e13.5)') (nrm2(Fpe(:,0,l) +&
                                                        & Fpe(:,1,l)), l = j, k)
                        end if
                        if (mulorder >= 2) then
                            write(lum2,'(6e13.5)') (nrm2(Fpe(:,0,l) +&
                                                        & Fpe(:,1,l) +&
                                                        & Fpe(:,2,l)), l = j, k)
                        end if
                        if (mulorder >= 3) then
                            write(lum3,'(6e13.5)') (nrm2(Fpe(:,0,l) +&
                                                        & Fpe(:,1,l) +&
                                                        & Fpe(:,2,l) +&
                                                        & Fpe(:,3,l)), l = j, k)
                        end if
                        if (mulorder >= 4) then
                            write(lum4,'(6e13.5)') (nrm2(Fpe(:,0,l) +&
                                                        & Fpe(:,1,l) +&
                                                        & Fpe(:,2,l) +&
                                                        & Fpe(:,3,l) +&
                                                        & Fpe(:,4,l)), l = j, k)
                        end if
                        if (mulorder >= 5) then
                            write(lum5,'(6e13.5)') (nrm2(Fpe(:,0,l) +&
                                                        & Fpe(:,1,l) +&
                                                        & Fpe(:,2,l) +&
                                                        & Fpe(:,3,l) +&
                                                        & Fpe(:,4,l) +&
                                                        & Fpe(:,5,l)), l = j, k)
                        end if
                    end do
                    if (mulorder >= 0) then
                        close(lum0)
                    end if
                    if (mulorder >= 1) then
                        close(lum1)
                    end if
                    if (mulorder >= 2) then
                        close(lum2)
                    end if
                    if (mulorder >= 3) then
                        close(lum3)
                    end if
                    if (mulorder >= 4) then
                        close(lum4)
                    end if
                    if (mulorder >= 5) then
                        close(lum5)
                    end if
                else
                    do l = 1, 3
                        write(cl,*) l
                        tcl = trim(adjustl(cl))
                        if (mulorder >= 0) then
                            call openfile('m0_field_'//tcl//'.cube', lum0,&
                                          'new', 'formatted')
                        end if
                        if (mulorder >= 1) then
                            call openfile('m1_field_'//tcl//'.cube', lum1,&
                                          'new', 'formatted')
                        end if
                        if (mulorder >= 2) then
                            call openfile('m2_field_'//tcl//'.cube', lum2,&
                                          'new', 'formatted')
                        end if
                        if (mulorder >= 3) then
                            call openfile('m3_field_'//tcl//'.cube', lum3,&
                                          'new', 'formatted')
                        end if
                        if (mulorder >= 4) then
                            call openfile('m4_field_'//tcl//'.cube', lum4,&
                                          'new', 'formatted')
                        end if
                        if (mulorder >= 5) then
                            call openfile('m5_field_'//tcl//'.cube', lum5,&
                                          'new', 'formatted')
                        end if
                        do i = 1, mulorder + 1
                            if (i == 1) then
                                lu = lum0
                            else if (i == 2) then
                                lu = lum1
                            else if (i == 3) then
                                lu = lum2
                            else if (i == 4) then
                                lu = lum3
                            else if (i == 5) then
                                lu = lum4
                            else if (i == 6) then
                                lu = lum5
                            end if
                            write(lu,'(a)') 'PE electric field component '//tcl
                            write(lu,'(a)') 'Generated by the Polarizable&
                                            & Embedding module'
                            write(lu,'(i5,3f12.6)') qmnucs, origin
                            write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0d0,&
                                                    & 0.0d0
                            write(lu,'(i5,3f12.6)') ysteps, 0.0d0, step(2),&
                                                    & 0.0d0
                            write(lu,'(i5,3f12.6)') zsteps, 0.0d0, 0.0d0,&
                                                    & step(3)
                            do j = 1, qmnucs
                                write(lu,'(i5,4f12.6)') nint(Zm(1,j)),&
                                                        & Zm(1,j), Rm(:,j)
                            end do
                        end do
                        do i = 1, xsteps * ysteps
                            j = (i - 1) * zsteps + 1
                            k = j - 1 + zsteps
                            if (mulorder >= 0) then
                                write(lum0,'(6e13.5)') Fpe(l,0,j:k)
                            end if
                            if (mulorder >= 1) then
                                write(lum1,'(6e13.5)') (Fpe(l,0,j:k) +&
                                                       & Fpe(l,1,j:k))
                            end if
                            if (mulorder >= 2) then
                                write(lum2,'(6e13.5)') (Fpe(l,0,j:k) +&
                                                       & Fpe(l,1,j:k) +&
                                                       & Fpe(l,2,j:k))
                            end if
                            if (mulorder >= 3) then
                                write(lum3,'(6e13.5)') (Fpe(l,0,j:k) +&
                                                       & Fpe(l,1,j:k) +&
                                                       & Fpe(l,2,j:k) +&
                                                       & Fpe(l,3,j:k))
                            end if
                            if (mulorder >= 4) then
                                write(lum4,'(6e13.5)') (Fpe(l,0,j:k) +&
                                                       & Fpe(l,1,j:k) +&
                                                       & Fpe(l,2,j:k) +&
                                                       & Fpe(l,3,j:k) +&
                                                       & Fpe(l,4,j:k))
                            end if
                            if (mulorder >= 5) then
                                write(lum5,'(6e13.5)') (Fpe(l,0,j:k) +&
                                                       & Fpe(l,1,j:k) +&
                                                       & Fpe(l,2,j:k) +&
                                                       & Fpe(l,3,j:k) +&
                                                       & Fpe(l,4,j:k) +&
                                                       & Fpe(l,5,j:k))
                            end if
                        end do
                        if (mulorder >= 0) then
                            close(lum0)
                        end if
                        if (mulorder >= 1) then
                            close(lum1)
                        end if
                        if (mulorder >= 2) then
                            close(lum2)
                        end if
                        if (mulorder >= 3) then
                            close(lum3)
                        end if
                        if (mulorder >= 4) then
                            close(lum4)
                        end if
                        if (mulorder >= 5) then
                            close(lum5)
                        end if
                    end do
                end if
            end if
            if (.not. lpol(1)) deallocate(Fpe)
        end if

        if (lpol(1)) then
            allocate(Fmuls(3*npols,1))
            Fmuls = 0.0d0
            l = 1
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                do j = 1, nsites
                    if (i == j) then
                        cycle
                    end if
                    exclude = .false.
                    do k = 1, lexlst
                        if (exclists(k,i) == exclists(1,j)) then
                            exclude = .true.
                            exit
                        end if
                    end do
                    if (exclude) cycle
                    Rji = Rs(:,i) - Rs(:,j)
                    if (lmul(0)) then
                        if (abs(maxval(M0s(:,j))) >= zero) then
                            call multipole_field(Fmuls(l:l+2,1), Rji, M0s(:,j))
                        end if
                    end if
                    if (lmul(1)) then
                        if (abs(maxval(M1s(:,j))) >= zero) then
                            call multipole_field(Fmuls(l:l+2,1), Rji, M1s(:,j))
                        end if
                    end if
                    if (lmul(2)) then
                        if (abs(maxval(M2s(:,j))) >= zero) then
                            call multipole_field(Fmuls(l:l+2,1), Rji, M2s(:,j))
                        end if
                    end if
                    if (lmul(3)) then
                        if (abs(maxval(M3s(:,j))) >= zero) then
                            call multipole_field(Fmuls(l:l+2,1), Rji, M3s(:,j))
                        end if
                    end if
                    if (lmul(4)) then
                        if (abs(maxval(M4s(:,j))) >= zero) then
                            call multipole_field(Fmuls(l:l+2,1), Rji, M4s(:,j))
                        end if
                    end if
                    if (lmul(5)) then
                        if (abs(maxval(M5s(:,j))) >= zero) then
                            call multipole_field(Fmuls(l:l+2,1), Rji, M5s(:,j))
                        end if
                    end if
                end do
                l = l + 3
            end do
            allocate(M1inds(3*npols,1))
!            call multipole_fields(Fmuls(:,1))
            if (mep_extfld) then
                j = 1
                do i = 1, npols
                    Fmuls(j:j+2,1) = Fmuls(j:j+2,1) + extfld
                    j = j + 3
                end do
            end if
            call induced_moments(M1inds, Fmuls)
            deallocate(Fmuls)
#if defined(VAR_MPI)
            if (nprocs > 1) then
                call mpi_bcast(M1inds, 3*npols, rmpi,&
                              &0, comm, ierr)
            end if
#endif
            allocate(Find(3,npoints))
            Find = 0.0d0
            i = 1
            do point = site_start, site_finish
                l = 1
                do j = 1, nsites
                    if (zeroalphas(j)) cycle
                    Rsp = mepgrid(:,i) - Rs(:,j)
                    Fs = 0.0d0
                    call multipole_field(Fs, Rsp, M1inds(l:l+2,1))
                    Find(:,i) = Find(:,i) + Fs
                    l = l + 3
                end do
                i = i + 1
            end do
            deallocate(M1inds)
#if defined(VAR_MPI)
            if (myid == 0 .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + 3 * mepdists(i-1)
                end do
                call mpi_gatherv(mpi_in_place, 0, rmpi,&
                                &Find, 3*mepdists, displs, rmpi,&
                                &0, comm, ierr)
            else if (myid /= 0) then
                call mpi_gatherv(Find, 3*mepdists(myid), rmpi,&
                                &0, 0, 0, rmpi,&
                                &0, comm, ierr)
            end if
#endif
            if (myid == 0) then
                if (mep_fldnrm) then
                    call openfile('ind_field.cube', lu, 'new', 'formatted')
                    write(lu,'(a)') 'PE induced electric field norm'
                    write(lu,'(a)') 'Generated by the polarizable embedding&
                                    & (PE) module'
                    write(lu,'(i5,3f12.6)') qmnucs, origin
                    write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0d0, 0.0d0
                    write(lu,'(i5,3f12.6)') ysteps, 0.0d0, step(2), 0.0d0
                    write(lu,'(i5,3f12.6)') zsteps, 0.0d0, 0.0d0, step(3)
                    do j = 1, qmnucs
                        write(lu,'(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j), Rm(:,j)
                    end do
                    do i = 1, xsteps * ysteps
                        j = (i - 1) * zsteps + 1
                        k = j - 1 + zsteps
                        write(lu,'(6e13.5)') (nrm2(Find(:,l)), l = j, k)
                    end do
                    close(lu)
                else
                    do l = 1, 3
                        write(cl,*) l
                        tcl = trim(adjustl(cl))
                        call openfile('ind_field_'//tcl//'.cube', lu,&
                                      'new', 'formatted')
                        write(lu,'(a)') 'PE induced electric field'
                        write(lu,'(a)') 'Generated by the Polarizable&
                                        & Embedding module'
                        write(lu,'(i5,3f12.6)') qmnucs, origin
                        write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0d0, 0.0d0
                        write(lu,'(i5,3f12.6)') ysteps, 0.0d0, step(2), 0.0d0
                        write(lu,'(i5,3f12.6)') zsteps, 0.0d0, 0.0d0, step(3)
                        do j = 1, qmnucs
                            write(lu,'(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j),&
                                                    & Rm(:,j)
                        end do
                        do i = 1, xsteps * ysteps
                            j = (i - 1) * zsteps + 1
                            k = j - 1 + zsteps
                            write(lu,'(6e13.5)') Find(l,j:k)
                        end do
                        close(lu)
                    end do
                end if
            end if
            deallocate(Find)
        end if
    end if

end subroutine pe_compute_mep

!------------------------------------------------------------------------------

subroutine pe_diis_solver(Mkinds,Fs)

    real(dp), dimension(:,:), intent(out) :: Mkinds
    real(dp), dimension(:,:), intent(in) :: Fs

    integer :: lu, itdiis, info, ndiis
    integer :: i, j, k, l, m, n, o, p, q
    logical :: exclude, lexist
    logical :: converged = .false.
    real(dp) :: fe = 1.0d0
    real(dp) :: ft = 1.0d0
    real(dp) :: R, R3, R5, Rd, ai, aj, norm, redthr, eps_fac, eps_inf, error, chk_sum
    real(dp), parameter :: d3i = 1.0d0 / 3.0d0
    real(dp), parameter :: d6i = 1.0d0 / 6.0d0
    integer, parameter :: mxdiis = 50
    real(dp), dimension(:), allocatable :: diis_vec
    integer, dimension(:), allocatable :: ipiv
    real(dp), dimension(:,:), allocatable :: Mkinds_diis, diis_mat, Mkinds_tmp
    real(dp), dimension(:), allocatable :: T, Rij, Ftmp, M1tmp,temp1,temp2

    allocate(Mkinds_tmp(3*npols,mxdiis))
    allocate(Mkinds_diis(3*npols,mxdiis))
    allocate(T(6), Rij(3), Ftmp(3), M1tmp(3))
! Start guess is in Mkinds_diis(:,1)
    do n = 1, ndens
        l = 1
        do i = site_start, site_finish
            if (zeroalphas(i)) cycle
            call spmv(P1s(:,i), Fs(l:l+2,n), Mkinds_diis(l:l+2,1), 'L')
            write(luout,*) Mkinds_diis(l:l+2,1), i
            l = l + 3
        end do
        do itdiis = 1, mxdiis
            if (itdiis == 1) then
               Mkinds_tmp(:,1) = Mkinds_diis(:,1)
            else
               Mkinds_tmp(:,itdiis) = 0.0d0
               do j = 1, itdiis
                   Mkinds_tmp(:,itdiis) = Mkinds_tmp(:,itdiis) + diis_vec(j+1) * Mkinds_diis(:,j)
               end do
               deallocate(diis_vec)
            end if

!             Start Jacobi:
!             Mkinds_diis(:,itdiis+1) = M0*( F - M1* Mkinds_temp)       
     
            l = 1
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                if (pe_damp) then
                    ai = (P1s(1,i) + P1s(4,i) + P1s(6,i)) * d3i
                end if
                m = 1
                Ftmp = 0.0d0
                do j = site_start, site_finish
                    if (zeroalphas(j)) cycle
                    exclude = .false.
                    do k = 1, lexlst
                        if (exclists(k,i) == exclists(1,j)) then
                            exclude = .true.
                            exit
                        end if
                    end do
                    if (i == j .or. exclude) then
                        m = m + 3
                        cycle
                    end if
                    Rij = Rs(:,j) - Rs(:,i)
                    R = nrm2(Rij)
                    R3 = R**3
                    R5 = R**5
                    ! damping parameters
                    ! JPC A 102 (1998) 2399 & Mol. Sim. 32 (2006) 471
                    ! a = 2.1304 = damp
                    ! u = R / (alpha_i * alpha_j)**(1/6)
                    ! fe = 1-(a²u²/2+au+1)*exp(-au)
                    ! ft = 1-(a³u³/6+a²u²/2+au+1)*exp(-au)
                    if (pe_damp) then
                        aj = (P1s(1,j) + P1s(4,j) + P1s(6,j)) * d3i
                        Rd = damp * R / (ai * aj)**d6i
                        fe = 1.0d0 - (0.5d0 * Rd**2 + Rd + 1.0d0) * exp(-Rd)
                        ft = fe - d6i * Rd**3 * exp(-Rd)
                    end if
                    q = 1
                    do o = 1, 3
                        do p = o, 3
                            T(q) = 3.0d0 * Rij(o) * Rij(p) * ft / R5
                            if (o == p) then
                                T(q) = T(q) - fe / R3
                            end if
                            q = q + 1
                        end do
                    end do
                    call spmv(T, Mkinds_tmp(m:m+2,itdiis), Ftmp, 'L', 1.0d0, 1.0d0)
                    m = m + 3
                end do
     
!    #if defined(VAR_MPI)
!                if (myid == 0 .and. nprocs > 1) then
!                    call mpi_reduce(mpi_in_place, Ftmp, 3, rmpi, mpi_sum, 0,&
!                                   & comm, ierr)
!                else if (myid /= 0) then
!                    call mpi_reduce(Ftmp, 0, 3, rmpi, mpi_sum, 0, comm, ierr)
!                end if
!    #endif
     
                if (myid == 0) then
                    Ftmp = Ftmp + Fs(l:l+2,n)
                    call spmv(P1s(:,i), Ftmp, Mkinds_diis(l:l+2,itdiis), 'L')
                    write(luout,*) 'Mkinds_diis, itdiis, i', Mkinds_diis(l:l+2,itdiis), itdiis, i
                end if
     
!    if defined(VAR_MPI)
!               if (myid == 0 .and. nprocs > 1) then
!                   displs(0) = 0
!                   do j = 1, nprocs
!                       displs(j) = displs(j-1) + 3 * poldists(j-1)
!                   end do
!                   call mpi_scatterv(Mkinds(:,n), 3*poldists, displs, rmpi,&
!                                    & mpi_in_place, 0, rmpi, 0, comm, ierr)
!               else if (myid /= 0) then
!                   call mpi_scatterv(0, 0, 0, rmpi, Mkinds(:,n),&
!                                    & 3*poldists(myid), rmpi, 0, comm, ierr)
!               end if
!    endif
                    l = l + 3
            end do ! i = 1, nsites
               
!      Check convergance
            error = sqrt(dot((Mkinds_diis(:,itdiis)-Mkinds_tmp(:,itdiis)),(Mkinds_diis(:,itdiis)-Mkinds_tmp(:,itdiis)))/(3*npols))
            if (myid == 0) then
                if (error < thriter) then
                    if (pe_verbose) then
                        write (luout,'(4x,a,i2,a)') 'Induced dipole moments&
                                                    & converged in ', itdiis,&
                                                    & ' iterations.'
                    end if
                    converged = .true.
                else if (itdiis == mxdiis) then
                    write(luout,*) 'ERROR: could not converge induced dipole&
                                   & moments.'
                    stop 'ERROR: could not converge induced dipole moments.'
                else
                ! solve DIIS equations
                    converged = .false.
                    ndiis = itdiis + 1
       
                    allocate(diis_mat(ndiis,ndiis))
                    allocate(diis_vec(ndiis))
                    allocate(temp1(ndiis))
                    allocate(temp2(ndiis))
          
                    diis_mat = 0.0d0
                    diis_vec = 0.0d0
                    diis_mat(1,:) = -1.0d0
                    diis_mat(:,1) = -1.0d0
                    diis_mat(1,1) = 0.0d0
                    diis_vec(1) = -1.0d0
       
                    do i = 2, ndiis
                        do j = 2, ndiis 
                            temp1 = Mkinds_diis(:,i-1) - Mkinds_tmp(:,i-1)
                            temp2 = Mkinds_diis(:,j-1) - Mkinds_tmp(:,j-1)
                            diis_mat(i,j) = dot(temp1,temp2)
                        end do
                    end do
                    deallocate(temp1,temp2) 
                    if (pe_verbose) then
                        do i = 1, ndiis 
                            do j = 1, ndiis
                                write(luout,*) 'diis_mat, i, j',diis_mat(i,j), i, j
                            end do
                        end do
                    end if
                    allocate(ipiv(ndiis)) 
                    call dgetrf(ndiis,ndiis,diis_mat,ndiis,ipiv,info)
                    call dgetrs('N', ndiis, 1, diis_mat, ndiis, ipiv, diis_vec, ndiis, info)
                    chk_sum = 0.0d0
                    if (pe_verbose) then
                        write(luout,*) 'diis_vec',diis_vec(1)
                        do i = 2, ndiis
                            chk_sum = chk_sum + diis_vec(i)
                            write(luout,*) 'diis_vec',diis_vec(i)
                        end do
                        write(luout,*) 'Sum of weights in diis_vec', chk_sum
                    end if
                    deallocate(diis_mat,ipiv)
     
                end if
            end if
!     #if defined(VAR_MPI)
!                 if (nprocs > 1) then
!                     call mpi_bcast(converged, 1, lmpi, 0, comm, ierr)
!                 end if
!     #endif
            if (converged) then
                Mkinds(:,n) = Mkinds_diis(:,itdiis)
                exit
            end if
     
        end do !itdiis 
    end do ! n= 1, ndens

end subroutine pe_diis_solver

!------------------------------------------------------------------------------

subroutine pe_diis_solver_charges(Mkinds,Fs)

    real(dp), dimension(:,:), intent(out) :: Mkinds
    real(dp), dimension(:,:), intent(in) :: Fs

    integer :: lu, itdiis, info, ndiis
    integer :: i, j, n, m
    logical :: exclude, lexist
    logical :: converged = .false.
    integer, parameter :: mxdiis = 100
    real(dp) :: FACTOR, DSCALE
    real(dp) :: error, X, Y, Z, R2, DISM0, DUM, temp, bla, R
    real(dp) :: oner, oner2, oner3, xi, yi, zi, xj, yj, zj
    real(dp) :: dipjx, dipjy, dipjz , qj
    real(dp), dimension(nsurp) :: QFIX, QNEW, VFIX2
    real(dp), dimension(:,:), allocatable :: tmpmat, dimat, field2
    real(dp), dimension(:), allocatable :: tmp, ipvt
    real(dp), dimension(:,:,:), allocatable :: qrep

    allocate(field2(3*npols, ndens))
    allocate(tmpmat(mxdiis+1,mxdiis+1))
    allocate(dimat(mxdiis+1,mxdiis+1))
    allocate(tmp(mxdiis+1))
    allocate(qrep(nsurp,mxdiis,2))
    allocate(ipvt(mxdiis+1)) 

    if (myid == 0) then
        inquire(file='pe_induced_charges.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
    end if
#endif

    if (lexist .and. (fock .or. energy)) then
        if (myid == 0) then
            call openfile('pe_induced_charges.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) QFIX 
            close(lu)
        end if
    else
        QFIX = 0.0d0
    end if
    IF(NTSATM.EQ. 60) DISM0 = 0.4000D+00*aa2au 
    IF(NTSATM.EQ.240) DISM0 = 0.2000D+00*aa2au 
    IF(NTSATM.EQ.960) DISM0 = 0.1000D+00*aa2au 
    FACTOR  = 1.0d0/SQRT(4.0d0*pi)/1.07D+00
    DSCALE   = (eps - 1.0d0)/eps 
    print *, dscale
    
    do n = 1, ndens
      
        do itdiis = 1, mxdiis

!             Start Jacobi:
!             Mkinds_diis(:,itdiis+1) = M0*( F - M1* Mkinds_temp)       
!     -- VFIX2: POTENTIAL AT NFFTS DUE TO SURFACE CHARGES --
           vfix2 = 0.0d0
           do i = 1, nsurp
              if(Sa(i) ==  0.0d0) cycle   
              do j = i+1, nsurp
                 if(Sa(j) ==  0.0d0) cycle   
                 X   = Sp(1,i) - Sp(1,j)
                 Y   = Sp(2,i) - Sp(2,j)
                 Z   = Sp(3,i) - Sp(3,j)
                 R2  = X*X+Y*Y+Z*Z
                 R   = SQRT(R2)
                 IF(R.GT.DISM0) THEN
                    VFIX2(I)=VFIX2(I) + QFIX(J)/R
                    VFIX2(J)=VFIX2(J) + QFIX(I)/R
                 ELSE
                    DUM = 1.5D+00/DISM0 - 0.5d0*R2/DISM0**3
                    VFIX2(I)=VFIX2(I) + QFIX(J)*DUM
                    VFIX2(J)=VFIX2(J) + QFIX(I)*DUM
                 END IF
              end do   
           end do   
  
           if (lpol(1)) then
           ! calculate field at induced dipoles due to surface charges
               field2 = 0.0d0
               m = 1
               do i = 1, nsites 
                   xi = Rs(1,i)
                   yi = Rs(2,i)
                   zi = Rs(3,i)
                   do j = 1, nsurp
                       xj = Sp(1,j)
                       yj = Sp(2,j)
                       zj = Sp(3,j)
                       qj = qfix(j)*dscale
                       x = xi - xj
                       y = yi - yj
                       z = zi - zj
                       r2 = x*x + y*y +z*z
                       oner2 = 1.0d0/r2
                       oner = sqrt(oner2)
                       oner3 = oner*oner2
                       dum = qj*oner3
                       field2(m,n)   = field2(m,n)   - dum*x
                       field2(m+1,n) = field2(m+1,n) - dum*y
                       field2(m+2,n) = field2(m+2,n) - dum*z
                   end do
                   m = m + 3
               end do

! Get induced dipols
               call iterative_solver(Mkinds(1:3*npols,n:n), (Fs(1:3*npols,n:n)+field2(:,n:n)))

! Potential at surface charge due to mm pol's
               do i = 1, nsurp
                    xi = Sp(1,i)
                    yi = Sp(2,i)
                    zi = Sp(3,i)
                    m = 1
                    do j = 1, nsites
                        xj = Rs(1,j)
                        yj = Rs(2,j)
                        zj = Rs(3,j)
                        dipjx = Mkinds( m, n)
                        dipjy = Mkinds( m + 1, n)
                        dipjz = Mkinds( m + 2, n)
                        x = xi - xj
                        y = yi - yj
                        z = zi - zj
                        r2 = x*x + y*y + z*z
                        oner2 = 1.0d0/r2
                        oner  = sqrt(oner2)
                        oner3 = oner2*oner
                        vfix2(i) = vfix2(i) - (dipjx*x + dipjy*y + dipjz*z)*oner3
                        m = m + 3
                    end do
               end do
           end if

           do i = 1, nsurp 
              QNEW(i) = -FACTOR*SQRT(Sa(i))*(Fs(3*npols+i,n) + VFIX2(I))
           end do

!      Check convergance
            error = 0.0d0
            do i = 1, nsurp
                error = error + abs(QFIX(i) - QNEW(i))
            end do 
            if (myid == 0) then
                if (error < nsurp*fixtol) then
                    if (pe_verbose) then
                        write (luout,'(4x,a,i2,a)') 'Induced surface charges&
                                                    & converged in ', itdiis,&
                                                    & ' iterations.'
                    end if
                    converged = .true.
                else if (itdiis == mxdiis) then
                    write(luout,*) 'ERROR: could not converge induced surface&
                                   & charges.'
                    stop 'ERROR: could not converge induced surface charges.'
                else
                ! solve DIIS equations
                    converged = .false.
                    call fixdiis(nsurp,itdiis,mxdiis,QNEW,QFIX,DIMAT&
                               &,QREP,TMP,TMPMAT,IPVT,nsurp) 
                    CALL DCOPY(nsurp,QNEW,1,QFIX,1)
                end if
            end if
!     #if defined(VAR_MPI)
!                 if (nprocs > 1) then
!                     call mpi_bcast(converged, 1, lmpi, 0, comm, ierr)
!                 end if
!     #endif
            if (converged) then
                CALL DSCAL(nsurp,DSCALE,QFIX,1)
                CALL DCOPY(nsurp,QFIX,1,Mkinds(3*npols+1:,n),1)
                exit
            end if
     
        end do !itdiis 
    end do ! n= 1, ndens

    if (fock) then
        if (myid == 0) then
            call openfile('pe_induced_charges.bin', lu, 'unknown',&
                         & 'unformatted')
            rewind(lu)
            write(lu) Mkinds(3*npols+1:,1)
            close(lu)
        end if
    end if

end subroutine pe_diis_solver_charges


!------------------------------------------------------------------------------
!
! The following code has been implemented into the gamess program by Hui Li 
! and adabted for use in dalton by MNP.
!
! The tesselation is FIXPVA2
!     Nandun M. Thellamurege and Hui Li THE JOURNAL OF CHEMICAL PHYSICS 137, 246101 (2012)
!
!------------------------------------------------------------------------------
subroutine fixtes(all_centers, all_z)

    real(dp), dimension(960) :: AST
    real(dp), dimension(360) :: IDUM
    real(dp), dimension(24) :: THEV
    real(dp), dimension(24) :: FIV
    real(dp), dimension(3,960) :: CDTST
    real(dp), dimension(3,20) :: VT20
    real(dp), dimension(6,60) :: JVT1
    real(dp), dimension(122,3) :: CV
    real(dp), dimension(3,nsites+qmnucs) :: all_centers
    real(dp), dimension(3,40,960) :: DAIT
    integer, dimension(41,960) :: IDTMP
    real(dp), dimension(nsites+qmnucs) :: all_z
    real(dp), dimension(:), allocatable :: XTS, YTS, ZTS
    real(dp), dimension(:), allocatable :: AFIX, RFIX
    real(dp), dimension(:,:,:), allocatable :: DAI
    integer, dimension(:), allocatable :: IDATOM
    integer, dimension(:,:), allocatable :: IDDAI
    real(dp), dimension(:,:), allocatable :: TMPTS
    integer :: II, I, J, INUC, NFFTS, mxtemp, IFFAT
    integer :: KFFTS, ITS, JJJ, IFFTS, surf
    real(dp) :: TH, FI, CTH, STH, FIR, FIXA
    real(dp) :: GOLD, ONEGOLD, SQRT13
    PARAMETER (GOLD=1.618033988749895D+00, ONEGOLD=1.0D+00/GOLD, SQRT13=0.577350269189626D+00)
    
    EQUIVALENCE (IDUM(1),JVT1(1,1))
    MXTEMP = qmnucs + nsites
    IF(MXFFTS .eq. 1) MXFFTS=MXTEMP*NTSATM/4
    allocate(XTS(mxffts))
    allocate(YTS(mxffts))
    allocate(ZTS(mxffts))
    allocate(AFIX(mxffts))
    allocate(RFIX(mxffts))
    allocate(IDATOM(mxffts))
    allocate(IDDAI(41,mxffts))
    allocate(DAI(3,40,mxffts))
    allocate(TMPTS(3,mxffts))
    

    THEV  = (/0.6523581398D+00,1.107148718D+00,1.382085796D+00,  &
            &  1.759506858D+00,2.034443936D+00,2.489234514D+00,  &
            &                 0.3261790699D+00,0.5535743589D+00, &
            &  0.8559571251D+00,0.8559571251D+00,1.017221968D+00,&
            &  1.229116717D+00,1.229116717D+00,1.433327788D+00,  &
            &  1.570796327D+00,1.570796327D+00,1.708264866D+00,  &
            &  1.912475937D+00,1.912475937D+00,2.124370686D+00,  &
            &  2.285635528D+00,2.285635528D+00,2.588018295D+00,  &
            &  2.815413584D+00/)

    FIV   = (/ 0.6283185307D+00,0.0000000000D+00,                 & 
            &  0.6283185307D+00,0.0000000000D+00,0.6283185307D+00,&
            &  0.0000000000D+00,0.6283185307D+00,0.0000000000D+00,&
            &  0.2520539002D+00,1.004583161D+00,0.6283185307D+00, &
            &  0.3293628477D+00,0.9272742138D+00,0.0000000000D+00,&
            &  0.3141592654D+00,0.9424777961D+00,0.6283185307D+00,&
            &  0.2989556830D+00,0.9576813784D+00,0.0000000000D+00,&
            &  0.3762646305D+00,0.8803724309D+00,0.6283188307D+00,&
            &  0.0000000000D+00/)

    FIR   = 1.256637061D+00

    IDUM = (/1, 6, 2, 32, 36, 37, 1, 2, 3, 33, 32, 38, 1, 3, 4, 34,  &
     &   33, 39, 1, 4, 5, 35, 34, 40, 1, 5, 6, 36, 35, 41, 7, 2, 6, 51,&
     &   42, 37, 8, 3, 2, 47, 43, 38, 9, 4, 3, 48, 44, 39, 10, 5, 4,   &
     &   49, 45, 40, 11, 6, 5, 50, 46, 41, 8, 2, 12, 62, 47, 52, 9,    &
     &   3, 13, 63, 48, 53, 10, 4, 14, 64, 49, 54, 11, 5, 15, 65, 50,  &
     &   55, 7, 6, 16, 66, 51, 56, 7, 12, 2, 42, 57, 52, 8, 13, 3,     &
     &   43, 58, 53, 9, 14, 4, 44, 59, 54, 10, 15, 5, 45, 60, 55, 11,  &
     &   16, 6, 46, 61, 56, 8, 12, 18, 68, 62, 77, 9, 13, 19, 69, 63,  &
     &   78, 10, 14, 20, 70, 64, 79, 11, 15, 21, 71, 65, 80, 7, 16,    &
     &   17, 67, 66, 81, 7, 17, 12, 57, 67, 72, 8, 18, 13, 58, 68, 73, &
     &   9, 19, 14, 59, 69, 74, 10, 20, 15, 60, 70, 75, 11, 21, 16,    &
     &   61, 71, 76, 22, 12, 17, 87, 82, 72, 23, 13, 18, 88, 83, 73,   &
     &   24, 14, 19, 89, 84, 74, 25, 15, 20, 90, 85, 75, 26, 16, 21,   &
     &   91, 86, 76, 22, 18, 12, 82, 92, 77, 23, 19, 13, 83, 93, 78,   &
     &   24, 20, 14, 84, 94, 79, 25, 21, 15, 85, 95, 80, 26, 17, 16,   &
     &   86, 96, 81, 22, 17, 27, 102, 87, 97, 23, 18, 28, 103, 88, 98, &
     &   24, 19, 29, 104, 89, 99, 25, 20, 30, 105, 90, 100, 26, 21,    &
     &   31, 106, 91, 101, 22, 28, 18, 92, 107, 98, 23, 29, 19, 93,    &
     &   108, 99, 24, 30, 20, 94, 109, 100, 25, 31, 21, 95, 110, 101,  &
     &   26, 27, 17, 96, 111, 97, 22, 27, 28, 107, 102, 112, 23, 28,   &
     &   29, 108, 103, 113, 24, 29, 30, 109, 104, 114, 25, 30, 31,     &
     &   110, 105, 115, 26, 31, 27, 111, 106, 116, 122, 28, 27, 117,   &
     &   118, 112, 122, 29, 28, 118, 119, 113, 122, 30, 29, 119, 120,  &
     &   114, 122, 31, 30, 120, 121, 115, 122, 27, 31, 121, 117, 116 /)

!
!     ADAPTED FROM PCM CODE (IN GAMESS)
!     HUI LI, NOV 26, 2011
!
!     COORDINATES OF VERTICES OF TESSERAE IN A SPHERE WITH UNIT RADIUS.
!
!                                    1
!
!                                 4     5
!
!                              3     6     2
!
      CV(  1,1) =  0.0D+00
      CV(  1,2) =  0.0D+00
      CV(  1,3) =  1.0D+00
      CV(122,1) =  0.0D+00
      CV(122,2) =  0.0D+00
      CV(122,3) = -1.0D+00
      II=1
      DO I=1,24
         TH=THEV(I)
         FI=FIV(I)
         CTH=COS(TH)
         STH=SIN(TH)
         DO J=1,5
            FI=FI+FIR
            IF(J.EQ.1) FI=FIV(I)
            II=II+1
            CV(II,1)=STH*COS(FI)
            CV(II,2)=STH*SIN(FI)
            CV(II,3)=CTH
         END DO
      END DO
!
      IF(NTSATM.EQ.4) THEN
         VT20(1,1)  = 1.0D+00
         VT20(2,1)  = 1.0D+00
         VT20(3,1)  = 1.0D+00
         VT20(1,2)  =-1.0D+00
         VT20(2,2)  =-1.0D+00
         VT20(3,2)  = 1.0D+00
         VT20(1,3)  =-1.0D+00
         VT20(2,3)  = 1.0D+00
         VT20(3,3)  =-1.0D+00
         VT20(1,4)  = 1.0D+00
         VT20(2,4)  =-1.0D+00
         VT20(3,4)  =-1.0D+00
         CALL DSCAL(12,SQRT13,VT20,1)
      END IF
!
      IF(NTSATM.EQ.6) THEN
         VT20(1,1)  = 1.0D+00
         VT20(2,1)  = 0.0D+00
         VT20(3,1)  = 0.0D+00
         VT20(1,2)  =-1.0D+00
         VT20(2,2)  = 0.0D+00
         VT20(3,2)  = 0.0D+00
         VT20(1,3)  = 0.0D+00
         VT20(2,3)  = 1.0D+00
         VT20(3,3)  = 0.0D+00
         VT20(1,4)  = 0.0D+00
         VT20(2,4)  =-1.0D+00
         VT20(3,4)  = 0.0D+00
         VT20(1,5)  = 0.0D+00
         VT20(2,5)  = 0.0D+00
         VT20(3,5)  = 1.0D+00
         VT20(1,6)  = 0.0D+00
         VT20(2,6)  = 0.0D+00
         VT20(3,6)  =-1.0D+00
      END IF
!
      IF(NTSATM.EQ.8) THEN
         VT20(1,1)  = 1.0D+00
         VT20(2,1)  = 1.0D+00
         VT20(3,1)  = 1.0D+00
         VT20(1,2)  =-1.0D+00
         VT20(2,2)  = 1.0D+00
         VT20(3,2)  = 1.0D+00
         VT20(1,3)  = 1.0D+00
         VT20(2,3)  =-1.0D+00
         VT20(3,3)  = 1.0D+00
         VT20(1,4)  = 1.0D+00
         VT20(2,4)  = 1.0D+00
         VT20(3,4)  =-1.0D+00
         VT20(1,5)  =-1.0D+00
         VT20(2,5)  =-1.0D+00
         VT20(3,5)  = 1.0D+00
         VT20(1,6)  =-1.0D+00
         VT20(2,6)  = 1.0D+00
         VT20(3,6)  =-1.0D+00
         VT20(1,7)  = 1.0D+00
         VT20(2,7)  =-1.0D+00
         VT20(3,7)  =-1.0D+00
         VT20(1,8)  =-1.0D+00
         VT20(2,8)  =-1.0D+00
         VT20(3,8)  =-1.0D+00
         CALL DSCAL(24,SQRT13,VT20,1)
      END IF
!
      IF(NTSATM.EQ.12) THEN
         VT20(1,1)  = 0.0D+00
         VT20(2,1)  = 1.0D+00
         VT20(3,1)  =    GOLD
         VT20(1,2)  = 0.0D+00
         VT20(2,2)  = 1.0D+00
         VT20(3,2)  =   -GOLD
         VT20(1,3)  = 0.0D+00
         VT20(2,3)  =-1.0D+00
         VT20(3,3)  =    GOLD
         VT20(1,4)  = 0.0D+00
         VT20(2,4)  =-1.0D+00
         VT20(3,4)  =   -GOLD
         VT20(3,5)  = 0.0D+00
         VT20(1,5)  = 1.0D+00
         VT20(2,5)  =    GOLD
         VT20(3,6)  = 0.0D+00
         VT20(1,6)  = 1.0D+00
         VT20(2,6)  =   -GOLD
         VT20(3,7)  = 0.0D+00
         VT20(1,7)  =-1.0D+00
         VT20(2,7)  =    GOLD
         VT20(3,8)  = 0.0D+00
         VT20(1,8)  =-1.0D+00
         VT20(2,8)  =   -GOLD
         VT20(2,9)  = 0.0D+00
         VT20(3,9)  = 1.0D+00
         VT20(1,9)  =    GOLD
         VT20(2,10) = 0.0D+00
         VT20(3,10) = 1.0D+00
         VT20(1,10) =   -GOLD
         VT20(2,11) = 0.0D+00
         VT20(3,11) =-1.0D+00
         VT20(1,11) =    GOLD
         VT20(2,12) = 0.0D+00
         VT20(3,12) =-1.0D+00
         VT20(1,12) =   -GOLD
         CALL DSCAL(36,0.525731112119134D+00,VT20,1)
      END IF
!
      IF(NTSATM.EQ.20) THEN
         VT20(1,1)  = 1.0D+00
         VT20(2,1)  = 1.0D+00
         VT20(3,1)  = 1.0D+00
         VT20(1,2)  = 1.0D+00
         VT20(2,2)  = 1.0D+00
         VT20(3,2)  =-1.0D+00
         VT20(1,3)  = 1.0D+00
         VT20(2,3)  =-1.0D+00
         VT20(3,3)  = 1.0D+00
         VT20(1,4)  =-1.0D+00
         VT20(2,4)  = 1.0D+00
         VT20(3,4)  = 1.0D+00
         VT20(1,5)  = 1.0D+00
         VT20(2,5)  =-1.0D+00
         VT20(3,5)  =-1.0D+00
         VT20(1,6)  =-1.0D+00
         VT20(2,6)  = 1.0D+00
         VT20(3,6)  =-1.0D+00
         VT20(1,7)  =-1.0D+00
         VT20(2,7)  =-1.0D+00
         VT20(3,7)  = 1.0D+00
         VT20(1,8)  =-1.0D+00
         VT20(2,8)  =-1.0D+00
         VT20(3,8)  =-1.0D+00
         VT20(1,9)  = 0.0D+00
         VT20(2,9)  = ONEGOLD
         VT20(3,9)  =    GOLD
         VT20(1,10) = 0.0D+00
         VT20(2,10) = ONEGOLD
         VT20(3,10) =   -GOLD
         VT20(1,11) = 0.0D+00
         VT20(2,11) =-ONEGOLD
         VT20(3,11) =    GOLD
         VT20(1,12) = 0.0D+00
         VT20(2,12) =-ONEGOLD
         VT20(3,12) =   -GOLD
         VT20(3,13) = 0.0D+00
         VT20(1,13) = ONEGOLD
         VT20(2,13) =    GOLD
         VT20(3,14) = 0.0D+00
         VT20(1,14) = ONEGOLD
         VT20(2,14) =   -GOLD
         VT20(3,15) = 0.0D+00
         VT20(1,15) =-ONEGOLD
         VT20(2,15) =    GOLD
         VT20(3,16) = 0.0D+00
         VT20(1,16) =-ONEGOLD
         VT20(2,16) =   -GOLD
         VT20(2,17) = 0.0D+00
         VT20(3,17) = ONEGOLD
         VT20(1,17) =    GOLD
         VT20(2,18) = 0.0D+00
         VT20(3,18) = ONEGOLD
         VT20(1,18) =   -GOLD
         VT20(2,19) = 0.0D+00
         VT20(3,19) =-ONEGOLD
         VT20(1,19) =    GOLD
         VT20(2,20) = 0.0D+00
         VT20(3,20) =-ONEGOLD
         VT20(1,20) =   -GOLD
         CALL DSCAL(60,SQRT13,VT20,1)
      END IF
!
      DO i = 1, nsites + qmnucs
         inuc = int(all_z(i))
         RFIX(I) = 2.400D+00*aa2au 

         IF(INUC.EQ. 0) RFIX(I) = 0.001D+00*aa2au 
         IF(INUC.EQ. 1) RFIX(I) = 1.400D+00*aa2au 
         IF(INUC.EQ. 3) RFIX(I) = 1.400D+00*aa2au 
         IF(INUC.EQ. 4) RFIX(I) = 1.400D+00*aa2au 
         IF(INUC.EQ. 5) RFIX(I) = 1.400D+00*aa2au 
         IF(INUC.EQ. 6) RFIX(I) = 2.100D+00*aa2au 
         IF(INUC.EQ. 7) RFIX(I) = 2.000D+00*aa2au 
         IF(INUC.EQ. 8) RFIX(I) = 1.900D+00*aa2au 
         IF(INUC.EQ. 9) RFIX(I) = 1.800D+00*aa2au 
         IF(INUC.EQ.10) RFIX(I) = 1.800D+00*aa2au 
         IF(INUC.EQ.11) RFIX(I) = 1.800D+00*aa2au 
         IF(INUC.EQ.12) RFIX(I) = 1.800D+00*aa2au 
         IF(INUC.EQ.13) RFIX(I) = 1.800D+00*aa2au 
         IF(INUC.EQ.14) RFIX(I) = 2.000D+00*aa2au 
         IF(INUC.EQ.15) RFIX(I) = 2.200D+00*aa2au 
         IF(INUC.EQ.16) RFIX(I) = 2.400D+00*aa2au 
         IF(INUC.EQ.17) RFIX(I) = 2.760D+00*aa2au 
         IF(INUC.EQ.18) RFIX(I) = 3.000D+00*aa2au 
! MNP LISTQM ????
!         IF(I.GT.NFFAT.AND.LISTQM(I).GT.0) RFIX(I) = 0.001D+00*TOBOHR
      END DO
!
!     -- NOTE: 'ME' STARTS AT 0 --
!
      xts = 0.0d0
      yts = 0.0d0
      zts = 0.0d0
      afix = 0.0d0
      idatom = 0
      dai = 0.0d0
      iddai = 0
!      NFFTS   = ME*(MXFFTS/NPROC - 1)  
!      IPCOUNT = ME - 1
      DO IFFAT = 1, nsites + qmnucs 
         IF(RFIX(IFFAT).LE.0.1D+00) cycle
!         IF(GOPARR) THEN
!           IPCOUNT = IPCOUNT + 1
!           IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 500
!         END IF
         CALL FIXPVA2(IFFAT,all_centers,JVT1,CV,CDTST,AST,  &
                    & XTS,YTS,ZTS,AFIX,RFIX,IDATOM,  &
                    & DAI,IDDAI,DAIT,IDTMP,TMPTS,VT20)
      end do
!     Parallel code not available in dalton yet
!     IF(GOPARR) THEN
!        CALL DDI_GSUMF(2418,XTS   ,      MXFFTS)
!        CALL DDI_GSUMF(2419,YTS   ,      MXFFTS)
!        CALL DDI_GSUMF(2420,ZTS   ,      MXFFTS)
!        CALL DDI_GSUMF(2421,AFIX  ,      MXFFTS)
!        CALL DDI_GSUMI(2422,IDATOM,      MXFFTS)
!        CALL DDI_GSUMF(2423,DAI   , 3*40*MXFFTS)
!        CALL DDI_GSUMI(2424,IDDAI ,   41*MXFFTS)
!     END IF
!
         
      KFFTS = 0
      DO ITS=1,MXFFTS
         IF(AFIX(ITS).LT.1.0D-04) cycle  ! 1.0D-04
         KFFTS = KFFTS + 1
         IF(KFFTS.GT.MXFFTS) THEN
            STOP 'FIXTES ERROR: PLEASE INCREASE MXFFTS.'
         END IF
         XTS(KFFTS)      = XTS(ITS)
         YTS(KFFTS)      = YTS(ITS)
         ZTS(KFFTS)      = ZTS(ITS)
         AFIX(KFFTS)     = AFIX(ITS)
         IDATOM(KFFTS)   = IDATOM(ITS)
         IDDAI(41,KFFTS) = IDDAI(41,ITS)
         DO JJJ = 1, IDDAI(41,ITS)
            IDDAI(JJJ,KFFTS) = IDDAI(JJJ,ITS)
            DAI(1,JJJ,KFFTS) = DAI(1,JJJ,ITS)
            DAI(2,JJJ,KFFTS) = DAI(2,JJJ,ITS)
            DAI(3,JJJ,KFFTS) = DAI(3,JJJ,ITS)
         END DO
      end do
      NFFTS = KFFTS   ! NFFTS IS GLOBAL
!
      FIXA = 0.0D+00
      DO IFFTS=1,NFFTS
         FIXA = FIXA + AFIX(IFFTS)
      END DO
!      write(luout,*) 'Surface point coordinates'
!      do iffts= 1, nffts
!         write(luout,*) XTS(IFFTS), YTS(IFFTS), ZTS(IFFTS)
!      end do
      FIXA = FIXA*TOANGS**2
!
      write(luout,*) 'Surface area FIXA', FIXA, '/A**2'
      write(luout,*) 'Number of surface charges', NFFTS
!     call openfile(surfile, surf, 'new', 'formatted')
!     rewind(lu)
!     write(surf,'(i6)') nffts
!     write(surf,'(a)') 'AA'
!     do i = 1, nffts
!         write(surf,'(4f12.6)') XTS(I), YTS(I), ZTS(I), AFIX(I)
!     end do
      allocate(Sp(3,nffts))
      allocate(Sa(nffts))
      nsurp = nffts
      do i = 1, nffts
         Sa(i) = AFIX(i)
         Sp(1,i) = XTS(i)
         Sp(2,i) = YTS(i)
         Sp(3,i) = ZTS(i)
      end do
 
!     close(lu)
!      RETURN
end subroutine fixtes

!-----------------------------------------------------------------------------

subroutine fixpva2(iffat,cord,jvt1,cv,cdtst,ast,&
                & xts,yts,zts,afix,rfix,idatom,&
                & dai, iddai, dait,idtmp,tmpts,vt20)

   real(dp), dimension(:), intent(out) :: XTS, YTS, ZTS, AFIX, RFIX

   real(dp), dimension(3,nsites+qmnucs) :: cord
   real(dp), dimension(960) :: AST
   real(dp), dimension(3,960) :: CDTST
   real(dp), dimension(3,40,960) :: DAIT
   integer, dimension(41,960) :: IDTMP
   integer, dimension(41) :: IDTMPTS
   real(dp), dimension(6,60) :: JVT1
   real(dp), dimension(122,3) :: CV
   real(dp), dimension(3,20) :: VT20
   real(dp), dimension(3,mxffts) :: TMPTS
   real(dp) :: XII, YII, ZII, RII, AREA0, DUMMY
   real(dp) :: P1X, P1Y, P1Z, P2X, P2Y, P2Z, P3X, P3Y, P3Z, P4X, P4Y, P4Z
   real(dp) :: P5X, P5Y, P5Z, P6X, P6Y, P6Z, P7X, P7Y, P7Z
   real(dp) :: PTS11, PTS21, PTS31, PTS12, PTS22, PTS32, PTS13, PTS23, PTS33
   real(dp) :: DNORM4, DNORM5, DNORM6, DNORM7, FACTOR, SCALE4, SCALE5, SCALE6, SCALE7
   integer :: IFFAT, KFFAT, ITS, KKK, MFFAT, III
   integer :: N1, N2, N3
   integer :: JJJ, KFFMX, JTS, KTS, LTS
   real(dp) :: FOURPI = 4.0d0*pi
   real(dp) :: ONE3RD=1.0D+00/3.0D+00
   real(dp), dimension(3,40,mxffts) :: DAI
   integer, dimension(nsites+qmnucs) :: IDATOM
   integer, dimension(41,mxffts) :: IDDAI


!     PARTITION OF THE CAVITY SURFACE INTO TESSERAE
!     HUI LI, NOV 26, 2011

     XII    = CORD(1,IFFAT)
     YII    = CORD(2,IFFAT)
     ZII    = CORD(3,IFFAT)
     RII    = RFIX(IFFAT)
     AREA0  = FOURPI*RII*RII/NTSATM
     cdtst  = 0.0d0
     ast    = 0.0d0
     dait   = 0.0d0
     idtmp  = 0



!     --- USE 20 TESSERAE FOR EACH SPHERE ---
     IF(NTSATM.EQ.20) THEN
     DO ITS = 1, NTSATM

        CDTST(1,ITS) = XII + VT20(1,ITS)*RII
        CDTST(2,ITS) = YII + VT20(2,ITS)*RII
        CDTST(3,ITS) = ZII + VT20(3,ITS)*RII

        tmpts = 0.0d0
        idtmpts = 0
        KFFAT=IFFAT
        AST(ITS) = AREA0
        CALL FIXPVASWF(KFFAT,CORD,CDTST,ITS,AST(ITS),RFIX,&
                     & TMPTS,IDTMPTS)
        IF(AST(ITS).GT.0.9D-04) THEN
           DO KKK = 1, IDTMPTS(41)
              MFFAT = IDTMPTS(KKK)
              DUMMY=ABS(TMPTS(1,MFFAT))+&
                  & ABS(TMPTS(2,MFFAT))+&
                  & ABS(TMPTS(3,MFFAT))
              IF(DUMMY.GT.1.0D-05) THEN
                 IDTMP(41,ITS) = IDTMP(41,ITS) + 1
                 III = IDTMP(41,ITS)
                 IDTMP(III,ITS)  = MFFAT
                 DAIT(1,III,ITS) = TMPTS(1,MFFAT)
                 DAIT(2,III,ITS) = TMPTS(2,MFFAT)
                 DAIT(3,III,ITS) = TMPTS(3,MFFAT)
              END IF
           ENDDO
        END IF
     end do
     END IF



!     --- USE 60 TESSERAE FOR EACH SPHERE ---
     IF(NTSATM.EQ.60) THEN
     DO ITS = 1, NTSATM

        N1 = JVT1(1,ITS)
        N2 = JVT1(2,ITS)
        N3 = JVT1(3,ITS)
        P4X = (CV(N1,1)+CV(N2,1)+CV(N3,1))*ONE3RD
        P4Y = (CV(N1,2)+CV(N2,2)+CV(N3,2))*ONE3RD
        P4Z = (CV(N1,3)+CV(N2,3)+CV(N3,3))*ONE3RD
        DNORM4 = P4X**2+P4Y**2+P4Z**2
        SCALE4 = 1.0D+00/SQRT(DNORM4)
        CDTST(1,ITS) = XII + P4X*SCALE4*RII
        CDTST(2,ITS) = YII + P4Y*SCALE4*RII
        CDTST(3,ITS) = ZII + P4Z*SCALE4*RII

        tmpts = 0.0d0
        idtmpts = 0
        KFFAT=IFFAT
        AST(ITS) = AREA0
        CALL FIXPVASWF(KFFAT,CORD,CDTST,ITS,AST(ITS),RFIX,&
                     & TMPTS,IDTMPTS)
        IF(AST(ITS).GT.0.9D-04) THEN
           DO KKK = 1, IDTMPTS(41)
              MFFAT = IDTMPTS(KKK)
              DUMMY=ABS(TMPTS(1,MFFAT))+&
                  & ABS(TMPTS(2,MFFAT))+&
                  & ABS(TMPTS(3,MFFAT))
              IF(DUMMY.GT.1.0D-05) THEN
                 IDTMP(41,ITS) = IDTMP(41,ITS) + 1
                 III = IDTMP(41,ITS)
                 IDTMP(III,ITS)  = MFFAT
                 DAIT(1,III,ITS) = TMPTS(1,MFFAT)
                 DAIT(2,III,ITS) = TMPTS(2,MFFAT)
                 DAIT(3,III,ITS) = TMPTS(3,MFFAT)
              END IF
           ENDDO
        END IF
     end do
     END IF



!     --- USE 240 TESSERAE FOR EACH SPHERE ---
     IF(NTSATM.EQ.240) THEN
     DO KTS = 1, 60

     DO JTS = 1, 4
        ITS = (KTS-1)*4 + JTS
        IF(JTS.EQ.1) THEN
           N1 = JVT1(1,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(5,KTS)
        END IF
        IF(JTS.EQ.2) THEN
           N1 = JVT1(6,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(5,KTS)
        END IF
        IF(JTS.EQ.3) THEN
           N1 = JVT1(3,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(6,KTS)
        END IF
        IF(JTS.EQ.4) THEN
           N1 = JVT1(2,KTS)
           N2 = JVT1(6,KTS)
           N3 = JVT1(5,KTS)
        END IF
        P4X = (CV(N1,1)+CV(N2,1)+CV(N3,1))*ONE3RD
        P4Y = (CV(N1,2)+CV(N2,2)+CV(N3,2))*ONE3RD
        P4Z = (CV(N1,3)+CV(N2,3)+CV(N3,3))*ONE3RD
        DNORM4 = P4X**2+P4Y**2+P4Z**2
        SCALE4 = 1.0D+00/SQRT(DNORM4)
        CDTST(1,ITS) = XII + P4X*SCALE4*RII
        CDTST(2,ITS) = YII + P4Y*SCALE4*RII
        CDTST(3,ITS) = ZII + P4Z*SCALE4*RII
        KKK = MOD(ITS,4)
        IF(KKK.EQ.1) FACTOR = 0.97527227808D+00
        IF(KKK.EQ.2) FACTOR = 1.04680294088D+00
        IF(KKK.EQ.3) FACTOR = 0.98896281888D+00
        IF(KKK.EQ.0) FACTOR = 0.98896281888D+00

        tmpts = 0.0d0
        idtmpts = 0
        KFFAT=IFFAT
        AST(ITS) = AREA0*FACTOR
        CALL FIXPVASWF(KFFAT,CORD,CDTST,ITS,AST(ITS),RFIX,&
                     & TMPTS,IDTMPTS)
        IF(AST(ITS).GT.0.9D-04) THEN
           DO KKK = 1, IDTMPTS(41)
              MFFAT = IDTMPTS(KKK)
              DUMMY=ABS(TMPTS(1,MFFAT))+&
                  & ABS(TMPTS(2,MFFAT))+&
                  & ABS(TMPTS(3,MFFAT))
              IF(DUMMY.GT.1.0D-05) THEN
                 IDTMP(41,ITS) = IDTMP(41,ITS) + 1
                 III = IDTMP(41,ITS)
                 IDTMP(III,ITS)  = MFFAT
                 DAIT(1,III,ITS) = TMPTS(1,MFFAT)
                 DAIT(2,III,ITS) = TMPTS(2,MFFAT)
                 DAIT(3,III,ITS) = TMPTS(3,MFFAT)
              END IF
           ENDDO
        END IF
     end do   
     end do  
     END IF



!     --- USE 960 TESSERAE FOR EACH SPHERE ---
     IF(NTSATM.EQ.960) THEN
     DO KTS = 1, 60

     DO JTS = 1, 4
        IF(JTS.EQ.1) THEN
           N1 = JVT1(1,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(5,KTS)
        END IF
        IF(JTS.EQ.2) THEN
           N1 = JVT1(6,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(5,KTS)
        END IF
        IF(JTS.EQ.3) THEN
           N1 = JVT1(3,KTS)
           N2 = JVT1(4,KTS)
           N3 = JVT1(6,KTS)
        END IF
        IF(JTS.EQ.4) THEN
           N1 = JVT1(2,KTS)
           N2 = JVT1(6,KTS)
           N3 = JVT1(5,KTS)
        END IF
        P1X = CV(N1,1)
        P1Y = CV(N1,2)
        P1Z = CV(N1,3)
        P2X = CV(N2,1)
        P2Y = CV(N2,2)
        P2Z = CV(N2,3)
        P3X = CV(N3,1)
        P3Y = CV(N3,2)
        P3Z = CV(N3,3)
!          COMPUTE THE COORDINATES OF POINTS 4, 5, 6
!
!                         1
!
!                      4     5
!
!                   3     6     2


        P4X = (P1X+P3X)/2.0D+00
        P4Y = (P1Y+P3Y)/2.0D+00
        P4Z = (P1Z+P3Z)/2.0D+00
        P5X = (P1X+P2X)/2.0D+00
        P5Y = (P1Y+P2Y)/2.0D+00
        P5Z = (P1Z+P2Z)/2.0D+00
        P6X = (P2X+P3X)/2.0D+00
        P6Y = (P2Y+P3Y)/2.0D+00
        P6Z = (P2Z+P3Z)/2.0D+00

        DNORM4 = P4X**2+P4Y**2+P4Z**2
        DNORM5 = P5X**2+P5Y**2+P5Z**2
        DNORM6 = P6X**2+P6Y**2+P6Z**2
        SCALE4 = 1.0D+00/SQRT(DNORM4)
        SCALE5 = 1.0D+00/SQRT(DNORM5)
        SCALE6 = 1.0D+00/SQRT(DNORM6)
        P4X = P4X*SCALE4
        P4Y = P4Y*SCALE4
        P4Z = P4Z*SCALE4
        P5X = P5X*SCALE5
        P5Y = P5Y*SCALE5
        P5Z = P5Z*SCALE5
        P6X = P6X*SCALE6
        P6Y = P6Y*SCALE6
        P6Z = P6Z*SCALE6

    DO  LTS = 1, 4
        ITS = ((KTS-1)*4 + JTS-1)*4 + LTS
        IF(LTS.EQ.1) THEN
          PTS11=P1X
          PTS21=P1Y
          PTS31=P1Z
          PTS12=P4X
          PTS22=P4Y
          PTS32=P4Z
          PTS13=P5X
          PTS23=P5Y
          PTS33=P5Z
        ELSE IF(LTS.EQ.2) THEN
          PTS11=P6X
          PTS21=P6Y
          PTS31=P6Z
          PTS12=P4X
          PTS22=P4Y
          PTS32=P4Z
          PTS13=P5X
          PTS23=P5Y
          PTS33=P5Z
        ELSE IF(LTS.EQ.3) THEN
          PTS11=P3X
          PTS21=P3Y
          PTS31=P3Z
          PTS12=P4X
          PTS22=P4Y
          PTS32=P4Z
          PTS13=P6X
          PTS23=P6Y
          PTS33=P6Z
        ELSE IF(LTS.EQ.4) THEN
          PTS11=P2X
          PTS21=P2Y
          PTS31=P2Z
          PTS12=P6X
          PTS22=P6Y
          PTS32=P6Z
          PTS13=P5X
          PTS23=P5Y
          PTS33=P5Z
        END IF

        P7X = (PTS11+PTS12+PTS13)*ONE3RD
        P7Y = (PTS21+PTS22+PTS23)*ONE3RD
        P7Z = (PTS31+PTS32+PTS33)*ONE3RD
        DNORM7 = P7X**2+P7Y**2+P7Z**2
        SCALE7 = 1.0D+00/SQRT(DNORM7)
        CDTST(1,ITS) = XII + P7X*SCALE7*RII
        CDTST(2,ITS) = YII + P7Y*SCALE7*RII
        CDTST(3,ITS) = ZII + P7Z*SCALE7*RII
        KKK = MOD(ITS,16)
        IF(KKK.EQ. 1) FACTOR = 0.96853843384D+00
        IF(KKK.EQ. 2) FACTOR = 0.98634758299D+00
        IF(KKK.EQ. 3) FACTOR = 0.97310155328D+00
        IF(KKK.EQ. 4) FACTOR = 0.97310155328D+00
        IF(KKK.EQ. 5) FACTOR = 1.04026074020D+00
        IF(KKK.EQ. 6) FACTOR = 1.05941468448D+00
        IF(KKK.EQ. 7) FACTOR = 1.04376817374D+00
        IF(KKK.EQ. 8) FACTOR = 1.04376817369D+00
        IF(KKK.EQ. 9) FACTOR = 0.98544774054D+00
        IF(KKK.EQ.10) FACTOR = 1.00020007354D+00
        IF(KKK.EQ.11) FACTOR = 0.98676349755D+00
        IF(KKK.EQ.12) FACTOR = 0.98343824680D+00
        IF(KKK.EQ.13) FACTOR = 0.98544774307D+00
        IF(KKK.EQ.14) FACTOR = 1.00020007615D+00
        IF(KKK.EQ.15) FACTOR = 0.98343824928D+00
        IF(KKK.EQ. 0) FACTOR = 0.98676350013D+00

        tmpts = 0.0d0
        idtmpts = 0
        KFFAT=IFFAT
        AST(ITS) = AREA0*FACTOR
        CALL FIXPVASWF(KFFAT,CORD,CDTST,ITS,AST(ITS),RFIX,&
                     & TMPTS,IDTMPTS)
       IF(AST(ITS).GT.0.9D-04) THEN
           DO KKK = 1, IDTMPTS(41)
              MFFAT = IDTMPTS(KKK)
              DUMMY=ABS(TMPTS(1,MFFAT))+&
                  & ABS(TMPTS(2,MFFAT))+&
                  & ABS(TMPTS(3,MFFAT))
              IF(DUMMY.GT.1.0D-05) THEN
                 IDTMP(41,ITS) = IDTMP(41,ITS) + 1
                 III = IDTMP(41,ITS)
                 IDTMP(III,ITS)  = MFFAT
                 DAIT(1,III,ITS) = TMPTS(1,MFFAT)
                 DAIT(2,III,ITS) = TMPTS(2,MFFAT)
                 DAIT(3,III,ITS) = TMPTS(3,MFFAT)
              END IF
           ENDDO
        END IF
     end do  ! lts = 1, 4
     end do  ! jts = 1, 4
     end do  ! kts = 1, 60
     END IF


     DO ITS=1,NTSATM
        IF(AST(ITS).LT.0.95D-04) cycle      !  0.95D-04
        NFFTS = NFFTS + 1
        if (NFFTS .GT. MXFFTS) then
           STOP 'FIXPVA2 ERROR: PLEASE INCREASE MXFFTS.'
        end if
!        IF(NFFTS.GT.(ME+1)*(MXFFTS/NPROC-1)) THEN
!           STOP 'FIXPVA2 ERROR: PLEASE INCREASE MXFFTS.'
!        END IF
        XTS(NFFTS)      = CDTST(1,ITS)
        YTS(NFFTS)      = CDTST(2,ITS)
        ZTS(NFFTS)      = CDTST(3,ITS)
        AFIX(NFFTS)     = AST(ITS)
        IDATOM(NFFTS)   = IFFAT
        IDDAI(41,NFFTS) = IDTMP(41,ITS)
        DO JJJ = 1, IDTMP(41,ITS)
           IDDAI(JJJ,NFFTS) = IDTMP(JJJ,ITS)
           DAI(1,JJJ,NFFTS) = DAIT(1,JJJ,ITS)
           DAI(2,JJJ,NFFTS) = DAIT(2,JJJ,ITS)
           DAI(3,JJJ,NFFTS) = DAIT(3,JJJ,ITS)
        ENDDO
     end do  

     RETURN

end subroutine fixpva2

!-----------------------------------------------------------------------------

subroutine FIXPVASWF(IFFAT,CORD,CDTST,ITS,AREA,RFIX,TMP,IDTMPTS)

   real(dp), intent(out) :: AREA
   integer, intent(in) :: iffat
   real(dp), dimension(3,nsites+qmnucs) :: cord, tmp
   real(dp), dimension(3,960) :: CDTST
   real(dp), dimension(mxffts) :: RFIX
   integer, dimension(41) :: IDTMPTS
   real(dp) :: ZERO=0.0d0
   real(dp) :: ONE=1.0d0
   real(dp) :: TWO=2.0d0
   real(dp) :: PT5=0.5d0
   real(dp) :: TEN=10.0d0
   real(dp) :: DISM0, ONEDISM0, DISN1, DISN2, RA, RB
   real(dp) :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X5, Y5, Z5
   real(dp) :: DISC2, ONEDISC2, DISC, ONEDISC, DISC3, DIS132, DIS13, DUM
   real(dp) :: DISC1, DISN, SWF2, DSWF2X, DSWF2Y, DSWF2Z
   real(dp) :: DUW, DX5DX3, DY5DX3, DZ5DX3, DX5DY3, DY5DY3, DZ5DY3
   real(dp) :: DX5DZ3, DY5DZ3, DZ5DZ3, DDISN2X, DDISN2Y, DDISN2Z
   real(dp) :: SWF1, DSWF1X, DSWF1Y, DSWF1Z, DISD
   real(dp) :: FB, FA, R1, R2, R3, R4, R5, R6, R7, R8, FAB, ONEFAB, DISM
   real(dp) :: ONEDISM, ONEDISD, DDDIS
   real(dp) :: CX3, CY3, CZ3, DX3, DY3, DZ3, R1X, R1Y, R1Z
   real(dp) :: R2X, R2Y, R2Z, R3X, R3Y, R3Z, R4X, R4Y, R4Z, R5X, R5Y, R5Z
   real(dp) :: R6X, R6Y, R6Z, R7X, R7Y, R7Z, R8X, R8Y, R8Z
   real(dp) :: FABX, FABY, FABZ, TEMP, DDISM2X, DDISM2Y, DDISM2Z, DUMMY
   integer :: JFFAT, KKK, ITS, KFFAT, i

!
!     DETERMINE THE AREA AND AREA DERIVATIVES FOR A TESSERA
!     USING FIXPVA2
!     HUI LI, NOV 27, 2011, LINCOLN
!
!     P1 = X1, Y1, Z1 (THE CURRENT TESSERA)
!     P2 = X2, Y2, Z2 (THE CENTER OF THE SPHERE OF THE TESSERA)
!     P3 = X3, Y3, Z3 (THE CENTER OF A NEIGHBOR SPHERE)
!
!
      DISM0 = 0.30D+00*aa2au 
      ONEDISM0 = ONE/DISM0
      DISN1 = 0.70D+00*aa2au 
      DISN2 = 1.20D+00*aa2au    ! DISN2 MUST BE < RA
!
      X1 = CDTST(1,ITS)
      Y1 = CDTST(2,ITS)
      Z1 = CDTST(3,ITS)
      X2 = CORD(1,IFFAT)
      Y2 = CORD(2,IFFAT)
      Z2 = CORD(3,IFFAT)
      RA = RFIX(IFFAT)
!
      DO JFFAT=1,nsites+qmnucs
!
!        IF AREA=0, RETURN
!        SKIP DISTANT SPHERES
!
         IF(JFFAT.EQ.IFFAT .OR. RFIX(JFFAT).LT.0.10D+00) cycle
         X3 = CORD(1,JFFAT)
         Y3 = CORD(2,JFFAT)
         Z3 = CORD(3,JFFAT)
         RB = RFIX(JFFAT)
         IF(ABS(X2-X3).GT.(RA+RB+DISN2)) cycle
         IF(ABS(Y2-Y3).GT.(RA+RB+DISN2)) cycle
         IF(ABS(Z2-Z3).GT.(RA+RB+DISN2)) cycle
         DISC2 = (X2-X3)**2 + (Y2-Y3)**2 + (Z2-Z3)**2
         ONEDISC2 = ONE/DISC2
         IF(DISC2.GE.(RA+RB+DISN2)**2) cycle
         IF(DISC2.LE.(RA-RB)**2 .AND. RB.LT.RA) cycle
         IF(DISC2.LE.(RA-RB)**2 .AND. RB.GT.RA) THEN
            AREA = ZERO
            RETURN
         END IF
         DISC   = SQRT(DISC2)
         ONEDISC= ONE/DISC
         DISC3  = DISC*DISC2
         DIS132 = (X1-X3)**2+(Y1-Y3)**2+(Z1-Z3)**2
         DIS13  = SQRT(DIS132)
!
         DUM = RB*ONEDISC
         X5  = X3 + (X2-X3)*DUM
         Y5  = Y3 + (Y2-Y3)*DUM
         Z5  = Z3 + (Z2-Z3)*DUM
         DISN = SQRT((X1-X5)**2+(Y1-Y5)**2+(Z1-Z5)**2)
         IF(DISN.GE.DISN2 .OR. DISC.LE.RB) THEN
            SWF2   = ONE
            DSWF2X = ZERO
            DSWF2Y = ZERO
            DSWF2Z = ZERO
         ELSE IF (DISN.LE.DISN1) THEN
            AREA   = ZERO
            RETURN
         ELSE
            DUW  = (DISN**2 - DISN1**2)/(DISN2**2-DISN1**2)
            SWF2 = 10.0D+00*DUW**3 - 15.0D+00*DUW**4 + 6.0D+00*DUW**5
            DUM    = RB/DISC3
            DX5DX3 = ONE - RB*ONEDISC + DUM*(X2-X3)*(X2-X3)
            DY5DX3 =                    DUM*(Y2-Y3)*(X2-X3)
            DZ5DX3 =                    DUM*(Z2-Z3)*(X2-X3)
            DX5DY3 =                    DUM*(X2-X3)*(Y2-Y3)
            DY5DY3 = ONE - RB*ONEDISC + DUM*(Y2-Y3)*(Y2-Y3)
            DZ5DY3 =                    DUM*(Z2-Z3)*(Y2-Y3)
            DX5DZ3 =                    DUM*(X2-X3)*(Z2-Z3)
            DY5DZ3 =                    DUM*(Y2-Y3)*(Z2-Z3)
            DZ5DZ3 = ONE - RB*ONEDISC + DUM*(Z2-Z3)*(Z2-Z3)
            DDISN2X=-TWO*((X1-X5)*DX5DX3+(Y1-Y5)*DY5DX3+(Z1-Z5)*DZ5DX3)
            DDISN2Y=-TWO*((X1-X5)*DX5DY3+(Y1-Y5)*DY5DY3+(Z1-Z5)*DZ5DY3)
            DDISN2Z=-TWO*((X1-X5)*DX5DZ3+(Y1-Y5)*DY5DZ3+(Z1-Z5)*DZ5DZ3)
            DUM = (30.0D+00*DUW**2-60.0D+00*DUW**3+30.0D+00*DUW**4)/&
               &  (DISN2**2-DISN1**2)
            DSWF2X = DUM*DDISN2X
            DSWF2Y = DUM*DDISN2Y
            DSWF2Z = DUM*DDISN2Z
         END IF
!
!
         IF(DISC.GE.(RA+RB)) THEN
            SWF1   = ONE
            DSWF1X = ZERO
            DSWF1Y = ZERO
            DSWF1Z = ZERO
         ELSE
            DISD = DIS13
            FB =  RA**2 + DISC2 - DISD**2
            FA =  RA**2 + DISC2 - RB**2
            R1 =  RA    + DISC  + DISD
            R2 =  RA    + DISC  - DISD
            R3 =  RA    - DISC  + DISD
            R4 = -RA    + DISC  + DISD
            R5 =  RA    + DISC  + RB
            R6 =  RA    + DISC  - RB
            R7 =  RA    - DISC  + RB
            R8 = -RA    + DISC  + RB
            FAB = SQRT(ABS(R1*R2*R3*R4*R5*R6*R7*R8))
            ONEFAB = ONE/FAB
            DISM = SQRT(ABS(TWO*RA*RA - (FA*FB + FAB)*PT5*ONEDISC2))
            IF(DISM.GE.DISM0) THEN
               IF(DIS13.LT.RB) THEN
                  AREA   = ZERO
                  RETURN
               END IF
               SWF1   = ONE
               DSWF1X = ZERO
               DSWF1Y = ZERO
               DSWF1Z = ZERO
            ELSE
               ONEDISM  = 1.0D+04
               IF(DISM.GT.1.0D-04) ONEDISM  = ONE/DISM
               ONEDISD  = ONE/DISD
               DDDIS    = (DISM-DISM0)*ONEDISM0
               SWF1     = PT5*DDDIS*DDDIS
               IF(DIS13.GE.RB) SWF1 = ONE - SWF1
               CX3 =  (X3-X2)*ONEDISC
               CY3 =  (Y3-Y2)*ONEDISC
               CZ3 =  (Z3-Z2)*ONEDISC
               DX3 =  (X3-X1)*ONEDISD
               DY3 =  (Y3-Y1)*ONEDISD
               DZ3 =  (Z3-Z1)*ONEDISD
               R1X =  CX3 + DX3
               R1Y =  CY3 + DY3
               R1Z =  CZ3 + DZ3
               R2X =  CX3 - DX3
               R2Y =  CY3 - DY3
               R2Z =  CZ3 - DZ3
               R3X = -CX3 + DX3
               R3Y = -CY3 + DY3
               R3Z = -CZ3 + DZ3
               R4X =  CX3 + DX3
               R4Y =  CY3 + DY3
               R4Z =  CZ3 + DZ3
               R5X =  CX3
               R5Y =  CY3
               R5Z =  CZ3
               R6X =  CX3
               R6Y =  CY3
               R6Z =  CZ3
               R7X = -CX3
               R7Y = -CY3
               R7Z = -CZ3
               R8X =  CX3
               R8Y =  CY3
               R8Z =  CZ3
               FABX = (R1X*R2*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2X*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3X*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4X*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5X*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6X*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7X*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7*R8X)*ONEFAB*PT5
               FABY = (R1Y*R2*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2Y*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3Y*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4Y*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5Y*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6Y*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7Y*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7*R8Y)*ONEFAB*PT5
               FABZ = (R1Z*R2*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2Z*R3*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3Z*R4*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4Z*R5*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5Z*R6*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6Z*R7*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7Z*R8+&
                    &  R1*R2*R3*R4*R5*R6*R7*R8Z)*ONEFAB*PT5
               TEMP = (FA*FB + FAB)*ONEDISC2*ONEDISC2
               DDISM2X= (X3-X2)*TEMP&
                    &  -((X1-X2)*FA+(X3-X2)*FB)*ONEDISC2&
                    &  -FABX*PT5*ONEDISC2
               DDISM2Y= (Y3-Y2)*TEMP&
                    &  -((Y1-Y2)*FA+(Y3-Y2)*FB)*ONEDISC2&
                    &  -FABY*PT5*ONEDISC2
               DDISM2Z= (Z3-Z2)*TEMP&
                    &  -((Z1-Z2)*FA+(Z3-Z2)*FB)*ONEDISC2&
                    &  -FABZ*PT5*ONEDISC2
               DUM =    (DISM-DISM0)*ONEDISM0*ONEDISM0&
                    &   *PT5*ONEDISM
               IF(DIS13.GE.RB) DUM = -DUM
               DSWF1X = DUM*DDISM2X
               DSWF1Y = DUM*DDISM2Y
               DSWF1Z = DUM*DDISM2Z
               IF(DSWF1X.GT. TEN) DSWF1X =  TEN  ! AVOID SINGULARITY WHEN DISM=0
               IF(DSWF1X.LT.-TEN) DSWF1X = -TEN
               IF(DSWF1Y.GT. TEN) DSWF1Y =  TEN
               IF(DSWF1Y.LT.-TEN) DSWF1Y = -TEN
               IF(DSWF1Z.GT. TEN) DSWF1Z =  TEN
               IF(DSWF1Z.LT.-TEN) DSWF1Z = -TEN
            END IF
         END IF
!
!        - NOTE: LOOP 150 MEANS MULTI-SPHERE SCALING
!                AREA MUST BE SCALED AFTER AREA DERIVATIVES
!
         TMP(1,JFFAT) = AREA*(SWF1*DSWF2X + SWF2*DSWF1X)
         TMP(2,JFFAT) = AREA*(SWF1*DSWF2Y + SWF2*DSWF1Y)
         TMP(3,JFFAT) = AREA*(SWF1*DSWF2Z + SWF2*DSWF1Z)
         IF(JFFAT.GE.2) THEN
            DO KFFAT = 1, JFFAT-1
               TMP(1,KFFAT) = TMP(1,KFFAT)*SWF1*SWF2
               TMP(2,KFFAT) = TMP(2,KFFAT)*SWF1*SWF2
               TMP(3,KFFAT) = TMP(3,KFFAT)*SWF1*SWF2
            ENDDO
         END IF
         AREA = AREA*SWF1*SWF2
         IF(AREA.LT.0.9D-04) THEN   ! 0.9D-04
            AREA = ZERO
            RETURN
         END IF
         DUMMY=ABS(TMP(1,JFFAT))+ABS(TMP(2,JFFAT))+ABS(TMP(3,JFFAT))
         IF(DUMMY.LT.0.9D-05) cycle
         IDTMPTS(41)  = IDTMPTS(41) + 1
         KKK          = IDTMPTS(41)
         IF(KKK.EQ.39) THEN
            STOP 'FIXPVASWF ERROR: IDTMPTS EXCEEDED 40.'
         END IF
         IDTMPTS(KKK) = JFFAT
      end do ! jffat =1,nsites+qmnucs  
!
!
      IDTMPTS(41)  = IDTMPTS(41) + 1
      KKK          = IDTMPTS(41)
      IDTMPTS(KKK) = IFFAT
      DO KFFAT = 1, nsites+qmnucs
         IF(KFFAT.NE.IFFAT) THEN
            TMP(1,IFFAT) = TMP(1,IFFAT) - TMP(1,KFFAT)
            TMP(2,IFFAT) = TMP(2,IFFAT) - TMP(2,KFFAT)
            TMP(3,IFFAT) = TMP(3,IFFAT) - TMP(3,KFFAT)
         END IF
      ENDDO
!
      RETURN
end subroutine FIXPVASWF

!-----------------------------------------------------------------------------

subroutine FIXDIIS(NFFPAR, NIT, MXDIIS, QOUT, QIN, DIMAT,&
                 & QREP, TMP, TMPMAT, IPVT, NSIZE)

    real(dp), dimension(:) :: qin, qout
    real(dp), dimension(:,:) :: DIMAT, tmpmat
    real(dp), dimension(:,:,:) :: qrep
    real(dp), dimension(:) :: ipvt, tmp
    real(dp), dimension(:), allocatable :: bla, bla2
    integer :: nitmax, nit0, I0, nffme, info
! ME used for parallelisation in GAMES here we just set it to zero
    integer :: ME, NIT, i, j, nffpar, NSIZE, MXDIIS, IMAX, IMIN

      nitmax = min(nit, mxdiis)
      nit0 = mod(nit-1,mxdiis) + 1

      ME = 0
      IMIN=ME*NFFPAR+1
      IMAX=MIN((ME+1)*NFFPAR,NSIZE)
      NFFME=IMAX-IMIN+1
!
!     -- STORE THE CHARGES
!        QREP(,,1)=QOUT
!        QREP(,,2)=QOUT-QIN
!
      CALL DCOPY(NFFME,QOUT(IMIN),1,QREP(1,NIT0,1),1)
      CALL DCOPY(NFFME,QOUT(IMIN),1,QREP(1,NIT0,2),1)
      CALL DAXPY(NFFME,-1.0d0,QIN(IMIN),1,QREP(1,NIT0,2),1)
!
!     -- UPGRADE THE INTERPOLATION MATRIX
!
      IF(NIT.GT.MXDIIS) THEN
         DO I=1,MXDIIS-1
            DO J=1,MXDIIS-1
               DIMAT(I+1,J+1)=DIMAT(I+2,J+2)
            ENDDO
         ENDDO
      END IF
!
      I0=NIT0
      allocate(bla(nffme))
      allocate(bla2(nffme))
      bla = qrep(:,nit0,2)
      DO I=NITMAX,1,-1
         bla2 = qrep(:,i0,2)
         TMP(I)=DOT(bla,bla2)
         I0=I0-1
         IF(I0.EQ.0) I0=I0+MXDIIS
      ENDDO
      deallocate(bla,bla2)
!
      DIMAT(NITMAX+1,1)=-1.0D+00
      DIMAT(1,NITMAX+1)=-1.0D+00
      DO I=NITMAX,1,-1
         DIMAT(NITMAX+1,I+1)=TMP(I)
         DIMAT(I+1,NITMAX+1)=TMP(I)
      ENDDO
!
!     -- AT THE FIRST ITERATION ONLY MATRIX INITIALIZATION
!
      IF (NIT.EQ.1) THEN
         DIMAT(1,1)=0.0D+00
         RETURN
      END IF
!
!     -- VECTOR INITIALIZATION
!
      tmp = 0.0d0
      TMP(1)=-1.0D+00
!
!     -- COPY THE MATRIX (SHOULD BE DESTROYED)
!
      DO I=1,NITMAX+1
         DO J=1,NITMAX+1
            TMPMAT(I,J)=DIMAT(I,J)
         ENDDO
      ENDDO
!
!     -- SOLVE THE LINEAR SYSTEM
!
!      call dgetrf(mxdiis+1,mxdiis+1,TMPMAT,mxdiis+1,ipvt,info)
      CALL DGEFA(TMPMAT,MXDIIS+1,NITMAX+1,IPVT,INFO)
      IF (INFO.NE.0) THEN
         WRITE(*,*) 'ERROR: SINGULAR MATRIX IN FIXDIIS.'
         STOP      
      END IF
!      call dgetrs('N', mxdiis+1, 1, TMPMAT, mxdiis+1, ipvt, TMP, mxdiis+1, info)
      CALL DGESL(TMPMAT,MXDIIS+1,NITMAX+1,IPVT,TMP,0)
!
!     -- INTERPOLATE
!
      qout = 0.0d0
      I0=NIT0
      DO I=NITMAX,1,-1
         CALL DAXPY(NFFME,TMP(I+1),QREP(1,I0,1),1,QOUT(IMIN),1)
         I0=I0-1
         IF(I0.EQ.0) I0=I0+MXDIIS
      ENDDO
      RETURN
!
end subroutine fixdiis

!------------------------------------------------------------------------------

subroutine pe_compute_london(fckmats)

    real(dp), dimension(:), intent(out) :: fckmats

    integer :: k, lu
    real(dp), dimension(:,:), allocatable :: Mkinds
    logical :: lexist

    fckmats = 0.0d0

    if (lmul(0)) call lao_multipoles(M0s, fckmats)
    if (lmul(1)) call lao_multipoles(M1s, fckmats)
    if (lmul(2)) call lao_multipoles(M2s, fckmats)

    if (lpol(1)) then
        allocate(Mkinds(3,npols))
        if (myid == 0) then
            inquire(file='pe_induced_dipoles.bin', exist=lexist)
        end if
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(lexist, 1, lmpi, 0, comm, ierr)
        end if
#endif
        if (lexist) then
            if (myid == 0) then
                call openfile('pe_induced_dipoles.bin', lu, 'old', 'unformatted')
                rewind(lu)
                read(lu) Mkinds
                close(lu)
            end if
        else
            stop 'ERROR: pe_induced_dipoles.bin does not exist'
        end if
        call lao_multipoles(Mkinds, fckmats, linduced=.true.)
        deallocate(Mkinds)
    endif

#if defined(VAR_MPI)
    if (myid == 0 .and. nprocs > 1) then
        call mpi_reduce(mpi_in_place, fckmats, ndens * n2bas, rmpi, mpi_sum,&
                       & 0, comm, ierr)
    else
        call mpi_reduce(fckmats, 0, ndens * n2bas, rmpi, mpi_sum, 0,&
                       & comm, ierr)
    end if
#endif

end subroutine pe_compute_london

!------------------------------------------------------------------------------

subroutine lao_multipoles(Mks, fckmats, linduced)

    real(dp), dimension(:,:), intent(in) :: Mks
    real(dp), dimension(:), intent(inout) :: fckmats
    logical, intent(in), optional :: linduced

    integer :: i, j, ifrom, ito
    integer :: site, ncomps, n2bas
    real(dp), dimension(:), allocatable :: symfacs
    real(dp), dimension(:,:,:), allocatable :: Mk_ints
    logical :: pol

    pol = .false.
    if (present(linduced)) pol = linduced

    ncomps = size(Mks, 1)
    n2bas = nbas * nbas

    ! notice change to n2bas instead of nnbas here
    ! also, the integrals now have a Bx, By and Bz factor on them
    allocate(Mk_ints(n2bas,3,ncomps))

    do site = site_start, site_finish
        if (abs(maxval(Mks(:,site))) < zero) then
            cycle
        end if
        if (pol .and. zeroalphas(site)) cycle

        call Mk_lao_integrals(Mk_ints, Rs(:,site), Mks(:,site))

        ! update for all directions of the magnetic field
        do j = 1, 3
            ifrom = (j - 1) * n2bas + 1
            ito   = j * n2bas
            fckmats(ifrom:ito) = fckmats(ifrom:ito) + sum(Mk_ints(:,j,:), 2)
        end do
    end do

    deallocate(Mk_ints)
    
end subroutine lao_multipoles

!-----------------------------------------------------------------------------

subroutine Mk_lao_integrals(Mk_ints, Rij, Mk)

    external :: Tk_lao_integrals

    real(dp), dimension(:,:,:), intent(out) :: Mk_ints
    real(dp), dimension(:), intent(in) :: Mk
    real(dp), dimension(3), intent(in) :: Rij

    integer :: i, j, k
    integer :: ncomps, nints, nwrk
    real(dp) :: taylor
    real(dp), dimension(:), allocatable :: factors

    k = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * size(Mk)) - 1.0d0)) - 1

    nwrk = size(work)

    if (mod(k,2) == 0) then
        taylor = 1.0d0 / factorial(k)
    else if (mod(k,2) /= 0) then
        taylor = - 1.0d0 / factorial(k)
    end if

    ncomps = size(Mk_ints, 3)
    nints = size(Mk_ints, 1)

    call Tk_lao_integrals(Mk_ints, nints, ncomps, Rij, work, nwrk)

    allocate(factors(ncomps))
    call symmetry_factors(factors)

    ! dot T^(k) integrals with multipole to get M^(k) integrals
    do i = 1, ncomps
        do j = 1, 3
            Mk_ints(:,j,i) = taylor * factors(i) * Mk(i) * Mk_ints(:,j,i)
        end do
    enddo

    deallocate(factors)

end subroutine Mk_lao_integrals

end module polarizable_embedding
