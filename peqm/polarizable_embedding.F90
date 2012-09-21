module polarizable_embedding

    use double_precision
    use blas_f90
    use lapack_f90

#if defined(VAR_MPI)
    use mpi
#endif

    implicit none

    private

    intrinsic :: allocated, present, min, minval, max, maxval, size, cpu_time

    public :: pe_dalton_input, pe_read_potential, pe_master
    public :: pe_save_density, pe_twoints
#if defined(VAR_MPI)
    public :: pe_mpi
#endif

    ! options
    logical, public, save :: peqm = .false.
    logical, save :: pe_iter = .true.
    logical, save :: pe_border = .false.
    logical, save :: pe_damp = .false.
    logical, save :: pe_gspol = .false.
    logical, save :: pe_nomb = .false.
    logical, save :: pe_gauss = .false.
    logical, public, save :: pe_polar = .false.
    logical, public, save :: pe_mep = .false.
    logical, public, save :: pe_skipqm = .false.
    logical, public, save :: pe_twoint = .false.
    logical, public, save :: pe_repuls = .false.
    logical, public, save :: pe_savden = .false.
    logical, public, save :: pe_pd = .false.
    logical, save :: pe_timing = .false.
    logical, save :: pe_infld = .false.

    ! calculation type
    logical, save :: fock = .false.
    logical, save :: energy = .false.
    logical, save :: response = .false.
    logical, save :: mep = .false.

    ! temporary solution for work array thing
    real(dp), dimension(:), pointer :: work

    ! MPI stuff
    integer, save :: myid, ncores, ierr
    logical, save :: initialized = .false.
    integer, dimension(:), save, allocatable :: npoldists, ndists, displs

    ! logical unit from dalton
    integer, save :: luout = 0

    ! constants, thresholds and stuff
    ! 1 bohr = 0.5291772108 Aa (codata 2002)
    real(dp), parameter :: aa2au = 1.8897261249935897d0
    real(dp), parameter :: pi = 3.141592653589793d0
    real(dp), parameter :: zero = 1.0d-6
    integer, save :: scfcycle = 0
    real(dp), save :: thriter = 1.0d-8
    real(dp), save :: damp = 2.1304d0
    real(dp), save :: gauss_factor = 1.0d0
    real(dp), save :: rep_factor = 1.0d0
    real(dp), save :: Rmin = 2.2d0
    character(len=6), save :: border_type = 'REDIST'
    logical, save :: chol = .true.

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476
    ! C^(n)_ij coefficients for calculating T(k) tensor elements
    integer, dimension(:,:,:), allocatable, save :: Cnij

    ! variables used for timings
    real(dp) :: t1, t2

    ! polarizable embedding potential info
    ! ------------------------------------

    ! number of sites
    integer, dimension(:), allocatable, public, save :: nsites
    ! number of polarizable sites
    integer, save :: npols = 0
    ! exclusion list length
    integer, save :: lexlst = 0

    ! specifies what type of parameters are present
    ! lmul(0): monopoles, lmul(1): dipoles etc.
    logical, dimension(0:5), public, save :: lmul = .false.
    ! lpol(1): (an)isotropic dipole-dipole polarizabilities
    logical, dimension(1), public, save :: lpol = .false.
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
    integer, dimension(:,:), allocatable, save :: exclists

    ! energy contributions
    ! electrostatic
    real(dp), dimension(:,:), allocatable, public, save :: Ees
    ! polarization
    real(dp), dimension(:,:), allocatable, public, save :: Epol

    ! multipole moments
    ! order of the highest order multipole moment
    integer, save :: mulorder = -1
    ! monopoles
    real(dp), dimension(:,:), allocatable, save :: M0s
    ! dipoles
    real(dp), dimension(:,:), allocatable, save :: M1s
    ! quadrupoles
    real(dp), dimension(:,:), allocatable, save :: M2s
    ! octopoles
    real(dp), dimension(:,:), allocatable, save :: M3s
    ! hexadecapoles
    real(dp), dimension(:,:), allocatable, save :: M4s
    ! ditriacontapoles
    real(dp), dimension(:,:), allocatable, save :: M5s

    ! (hyper)polarizabilities
    ! order of highest order polarizability
    integer, save :: polorder = -1
    ! dipole-dipole polarizabilities
    real(dp), dimension(:,:), allocatable, save :: P1s
    ! .true. if P1 > 0 else .false.
    logical, dimension(:), allocatable, save :: zeroalphas


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


    ! polarizable density fragment info
    ! ----------------------------

    ! number of polarizable densities
    integer, public, save :: npds = 0
    ! number of nuclei in current polarizable density
    integer, public, save :: pdnucs = 0
    ! nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zpd
    ! nuclear coordinates
    real(dp), dimension(:,:), allocatable, save :: Rpd
    ! energy contributions
    real(dp), dimension(:,:), allocatable, public, save :: Epd


    ! MEP stuff
    ! ---------
    ! create QM cubes
    logical, save :: mep_qmcube = .true.
    ! create multipole cubes
    logical, save :: mep_mulcube = .true.
    ! external electric field
    logical, save :: mep_extfld = .true.
    real(dp), dimension(3), save :: extfld = 0.0d0
    ! calculate electric field
    logical, save :: mep_field = .false.
    logical, save :: mep_fldnrm = .false.
    ! number of grid points
    integer, dimension(:), allocatable, save :: npoints
    ! point distribution
    integer, dimension(:), save, allocatable :: nmepdists
    ! grid points
    real(dp), dimension(:,:), allocatable, save :: mepgrid
    ! CUBE file origin and step sizes
    real(dp), dimension(3), save :: origin, step
    ! grid density in x, y and z direction
    integer, save :: xgrid = 20
    integer, save :: ygrid = 20
    integer, save :: zgrid = 20
    ! numberof steps in x, y and z direction
    integer, save :: xsteps
    integer, save :: ysteps
    integer, save :: zsteps
    ! box size relative to molecule size
    real(dp), save :: xsize = 5.0d0
    real(dp), save :: ysize = 5.0d0
    real(dp), save :: zsize = 5.0d0

    ! Internal field stuff
    ! --------------------
    ! Coordinates on which potential and field are calculated
    real(dp), dimension(:,:), allocatable, save :: crds
    ! Number of coordinates (length of crds)/3
    integer, save :: ncrds

! TODO:
! make only two routines accesible from outside module: pe_init and pe_master
! handle interface better, e.g. scale or remove higher order moments and pols
! write better output (e.g. in abadrv,F)
! damping of electric field from QM system?
! check for positive definiteness of the response matrix?
! write results after redistributing parameters
! better solution for lmul and lpol
! use allocate/deallocate where possible?
! insert quit if symmetry or QM3, QMMM etc.
! find better solution for electric field calculation from polarizable densities
! higher order polarizabilities
! write list of publications which should be cited
! write output related to FDs
! remove double zeroing and unecessary zeroing
! nonlinear response properties
! magnetic properties
! cutoffs and damping
! memory management (dalton work, pointer, allocate)
! add error catching
! parallelization (openMP, MPI, CUDA/openCL?)

contains

!------------------------------------------------------------------------------

subroutine pe_dalton_input(word, luinp, lupri)

    character(len=7), intent(inout) :: word
    integer, intent(in) :: luinp
    integer, intent(in) :: lupri

    character(len=7) :: option
    character(len=2) :: aaorau
    integer :: i
    luout = lupri

    do
        read(luinp,'(a7)') option
        call chcase(option)

        ! do a Polarizable Embedding calculation
        if (trim(option(2:)) == 'PEQM') then
            peqm = .true.
        ! direct solver for induced dipoles
        else if (trim(option(2:)) == 'DIRECT') then
            pe_iter = .false.
        ! iterative solver for induced dipoles (defaults to true)
        else if (trim(option(2:)) == 'ITERAT') then
            read(luinp,*) option
            backspace(luinp)
            if (option(1:1) /= '.' .and. option(1:1) /= '*'&
               &.and. option(1:1) /= '!' .and. option(1:1) /= '#') then
                read(luinp,*) thriter
            end if
            pe_iter = .true.
        ! handling sites near quantum-classical border
        else if (trim(option(2:)) == 'BORDER') then
            read(luinp,*) option
            backspace(luinp)
            if (option(1:1) /= '.' .and. option(1:1) /= '*' .and.&
            & option(1:1) /= '!') then
                read(luinp,*) border_type, Rmin, aaorau
                call chcase(border_type)
                if (trim(border_type) /= 'REMOVE' .and.&
                & trim(border_type) /= 'REDIST') then
                    stop 'ERROR: unknown handling of border sites!'
                end if
                call chcase(aaorau)
                if (trim(aaorau) == 'AA') Rmin = Rmin * aa2au
            end if
            pe_border = .true.
        ! induced dipole - induced dipole damping
        else if (trim(option(2:)) == 'DAMP') then
            read(luinp,*) option
            backspace(luinp)
            if (option(1:1) /= '.' .and. option(1:1) /= '*'&
               &.and. option(1:1) /= '!' .and. option(1:1) /= '#') then
                read(luinp,*) damp
            end if
            pe_damp = .true.
        ! neglect dynamic response from environment
        else if (trim(option(2:)) == 'GSPOL') then
            pe_gspol = .true.
        ! neglect many-body interactions
        else if (trim(option(2:)) == 'NOMB') then
            pe_nomb = .true.
        ! use Gaussian broadened multipoles and FD nuclear charges
        else if (trim(option(2:)) == 'GAUSS') then
            read(luinp,*) option
            backspace(luinp)
            if (option(1:1) /= '.' .and. option(1:1) /= '*'&
               &.and. option(1:1) /= '!' .and. option(1:1) /= '#') then
                read(luinp,*) gauss_factor
            end if
            pe_gauss = .true.
        ! calculate intermolecular two-electron integrals
        else if (trim(option(2:)) == 'TWOINT') then
            read(luinp,*) pdnucs
            pe_twoint = .true.
        ! save density matrix
        else if (trim(option(2:)) == 'SAVDEN') then
            pe_savden = .true.
        ! get fock matrix for repulsion potential
        else if (trim(option(2:)) == 'REPULS') then
            read(luinp,*) option
            backspace(luinp)
            if (option(1:1) /= '.' .and. option(1:1) /= '*'&
               &.and. option(1:1) /= '!' .and. option(1:1) /= '#') then
                read(luinp,*) rep_factor
            end if
            pe_repuls = .true.
        ! electrostatics from polarizable densities
        else if (trim(option(2:)) == 'PD') then
            ! number of polarizable densities
            read(luinp,*) npds
            pe_pd = .true.
        ! skip QM calculations, i.e. go directly into PE module
        else if (trim(option(2:)) == 'SKIPQM') then
            pe_skipqm = .true.
        ! calculate internal field
        else if (trim(option(2:)) == 'INFLD') then
            read(luinp,*) option
            backspace(luinp)
            if (option(1:1) /= '.' .and. option(1:1) /= '*'&
               &.and. option(1:1) /= '!' .and. option(1:1) /= '#') then
                read(luinp,*) ncrds
                allocate(crds(3,ncrds))
                do i=1, ncrds 
                    read(luinp,*) crds(1,i), crds(2,i), crds(3,i) 
                end do
            end if 
            pe_infld = .true.
        ! evaluate molecular electrostatic potential
        else if (trim(option(2:)) == 'MEP') then
            read(luinp,*) option
            backspace(luinp)
            if (option(1:1) /= '.' .and. option(1:1) /= '*'&
               &.and. option(1:1) /= '!' .and. option(1:1) /= '#') then
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
                        read(luinp,*) extfld(1), extfld(2), extfld(3)
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
            if (option(1:1) /= '.' .and. option(1:1) /= '*'&
               &.and. option(1:1) /= '!' .and. option(1:1) /= '#') then
                read(luinp,*) 
            end if
            pe_mep = .true.
        else if (option(1:1) == '*') then
            word = option
            exit
        else if (option(1:1) == '!' .or. option(1:1) == '#') then
            cycle
        end if
    end do

    if (pe_mep .and. peqm) stop 'PEQM and MEP are not compatible'
    if (pe_nomb .and. pe_iter) stop 'NOMB and ITERATIVE are not compatible'
    if (peqm .and. pe_savden) stop 'PEQM and SAVDEN are not compatible'
    if (peqm .and. pe_twoint) stop 'PEQM and TWOINT are not compatible'

end subroutine pe_dalton_input

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
            do j=1, nsites(ncores-1) !MM sites 
                Rsp = coords(:,i) - Rs(:,j)
                t_Vpe = 0.0d0
                if (lmul(0)) then
                    allocate(Tsp(1),factors(1))
                    call symmetry_factors(factors)
                    taylor = 1.0d0 / factorial(0)
                    call Tk_tensor(Tsp, Rsp)
                    t_Vpe(0) = t_Vpe(0) + taylor * factors(1) * Tsp(1) * M0s(1,j)
                    Vpe(0,i) = Vpe(0,i) + t_Vpe(0)
                    deallocate(Tsp, factors)
                end if
                if (lmul(1)) then
                    allocate(Tsp(3), factors(3))
                    call symmetry_factors(factors)
                    taylor = - 1.0d0 / factorial(1)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 3
                        t_Vpe(1) = t_Vpe(1) + taylor * factors(k) * Tsp(k) * M1s(k,j)
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
                        !Vpe(2,i) = Vpe(2,i) + taylor * factors(k) * Tsp(k) * M2s(k,j)
                        t_Vpe(2) = t_Vpe(2) + taylor * factors(k) * Tsp(k) * M2s(k,j)
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
                        t_Vpe(3) = t_Vpe(3) + taylor * factors(k) * Tsp(k) * M3s(k,j)
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
                        t_Vpe(4) = t_Vpe(4) + taylor * factors(k) * Tsp(k) * M4s(k,j)
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
                        t_Vpe(5) = t_Vpe(5) + taylor * factors(k) * Tsp(k) * M5s(k,j)
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
        call induced_dipoles(M1inds, Fmuls)
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
            do j=1, nsites(ncores-1)
                if (zeroalphas(j)) cycle
                Rsp = coords(:,i) - Rs(:,j)
                Tsp = 0.0d0
                t_Vind = 0.0d0
                t_Vind_qmconv = 0.0d0
                call Tk_tensor(Tsp, Rsp)
              !  write(luout,*) 'After call Tk_tensor'
                do k = 1, 3
                    t_Vind = t_Vind + taylor * factors(k) * Tsp(k) * M1inds(l+k,1)
                    t_Vind_qmconv = t_Vind_qmconv &
                                    + taylor * factors(k) * Tsp(k) * M1inds_qmconv(l+k,1)
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
            do j = 1, nsites(ncores-1)
                Rsp = coords(:,i) - Rs(:,j)
                allocate(t_Fpe(3))
                t_Fpe = 0.0d0
                if (lmul(0)) then
                    Fs = 0.0d0 !field from M0s
                    call multipole_field(Fs, Rsp, M0s(:,j))
                    Fpe(1:3,i) = Fpe(1:3,i) + Fs
                    t_Fpe(1:3) = Fs
                    !write(luout,*) 'Fpe from M0', Fpe(1,i), Fpe(2,i), Fpe(3,i)
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
            do j = 1, nsites(ncores-1)
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
                Vtot(1,i) = Vpe(0,i)+Vpe(1,i)
                Ftot(1:3,i) = Fpe(1:3,i)+Fpe(4:6,i)

            else if (mulorder == 2) then
                Vtot(1,i) = Vpe(0,i)+Vpe(1,i)+Vpe(2,i)
                Ftot(1:3,i) = Fpe(1:3,i)+Fpe(4:6,i)+Fpe(7:9,i)

            else if (mulorder == 3) then
                Vtot(1,i) = Vpe(0,i)+Vpe(1,i)+Vpe(2,i)+Vpe(3,i)
                Ftot(1:3,i) = Fpe(1:3,i)+Fpe(4:6,i)+Fpe(7:9,i)+Fpe(10:12,i)

            else if (mulorder == 4) then
                Vtot(1,i) = Vpe(0,i)+Vpe(1,i)+Vpe(2,i)+Vpe(3,i)+Vpe(4,i)
                Ftot(1:3,i) = Fpe(1:3,i)+Fpe(4:6,i)+Fpe(7:9,i)+Fpe(10:12,i)&
                            +Fpe(13:15,i)

            else if (mulorder == 5) then
                Vtot(1,i) = Vpe(0,1)+Vpe(1,i)+Vpe(2,i)+Vpe(3,i)+Vpe(5,i)&
                                +Vpe(5,i)
                Ftot(1:3,i) = Fpe(1:3,i)+Fpe(4:6,i)+Fpe(7:9,i)+Fpe(10:12,i)&
                                +Fpe(13:15,i)+Fpe(16:18,i)
            end if
            

            ! add contribution from induced dipoles to potential and field if present
!            if (lpol(1)) then 
                ! the order of E/Vtot_qmconv and E/Vtot IS important because E/Vtot is changed
!                Vtot_qmconv(1,i) = Vtot(1,i) + Vind_qmconv(1,i)
!                Vtot(1,i) = Vtot(1,i) + Vind(1,i)
!                Ftot_qmconv(1:3,i) = Ftot(1:3,i) + Find_qmconv(1:3,i)
!                Ftot(1:3,i) = Ftot(1:3,i) + Find(1:3,i)

            !calculate norm of Efield in a given site i
!                Fnrm_qmconv(1,i) = nrm2(Ftot_qmconv(1:3,i)) !calculate only if Fnrm_qmconv present
            Fnrm(1,i) = nrm2(Ftot(1:3,i))    
!            end if
        end do
    end if

    if (mulorder>=0) then
        write(luout,*) 'mulorder', mulorder
        write(luout,*)
        write(luout,*)
        write(luout,'(a)') '******************************************************' 
        write(luout,'(a)') '*** Internal electric potential and field analysis ***'
        write(luout,'(a)') '******************************************************' 
        write(luout,*)
        write(luout,*)

        if (lpol(1)) then
            write(luout,*)'-------------------------------------------------------&
                          &-----------------'
            write(luout,'(a)') '                    QM present                 &
                           &       QM absent             '
            write(luout,'(a)')  ' Site      Vtot         Vpe           Vind   &
                        &     Vtot         Vind    ' 
            write(luout,'(a)') '=======================================================&
                             &================='     
            do i = 1, ncoords
                write(luout,'(i4,2x,(5e13.5))') i, (Vtot(1,i) + Vind_qmconv(1,i)),&
                         Vtot(1,i), Vind_qmconv(1,i), (Vtot(1,i)+Vind(1,i)), Vind(1,i)       
            end do

            write(luout,'(a)') '=======================================================&
                          &================='
            write(luout,*)
    
            write(luout,*)'-------------------------------------------------------&
                         &------------------------------------------'
            write(luout,'(a)') '        QM present       &
                             &                                              '                                         
            write(luout,'(a)')  ' Site      Enrm                       Epe   &
                           &                              Eind          '
            write(luout,'(a)')  '                         x             y   &
                           &         z            x            y  &
                           &        z '
            write(luout,'(a)') '=======================================================&
                                &=========================================='
            do i = 1, ncoords
                write(luout,'(i4,2x,(7e13.5))') i, (nrm2(Ftot(1:3,i)+Find_qmconv(1:3,i))),&
                                    &Ftot(1:3,i),Find_qmconv(1:3,i) 
            end do
            write(luout,'(a)') '=======================================================&
                                &=========================================='

            write(luout,*)
    
            write(luout,*)'-------------------------------------------------------&
                            &------------------------------------------'
            write(luout,'(a)') '        QM absent       &                            
                            &                                              '                                        
            write(luout,'(a)')  ' Site      Enrm                       Epe   &
                           &                                Eind        '
            write(luout,'(a)')  '                         x             y   &
                           &         z            x            y  &
                           &        z '
            write(luout,'(a)') '=======================================================&
                                &=========================================='
            do i = 1, ncoords
                write(luout,'(i4,2x,(7e13.5))') i, (nrm2(Ftot(1:3,i)+Find(1:3,i))),&
                                        &Ftot(1:3,i),Find(1:3,i) 
            end do
            write(luout,'(a)') '=======================================================&
                                &=========================================='



        else 
            write(luout,*)'-------------------------------------------------------&
                          &-----------------'
            write(luout,'(a)')  ' Site      Vtot                       Etot   &
                       &                  Enrm          '
            write(luout,'(a)')  '                         x             y   &
                           &           z                      '
            write(luout,'(a)') '=======================================================&
                                &================='
            do i = 1, ncoords
                write(luout,'(i4,2x,(5e13.5))') i, Vtot(1,i),&
                         Ftot(1:3,i), (nrm2(Ftot(1:3,i)))
            end do
            write(luout,'(a)') '=======================================================&
                            &================='
            write(luout,*)
        end if

        write(luout,*)
        write(luout,'(a)')'Contributions from static multipoles MXs'
        write(luout,*) 

        do j = 1, mulorder + 1
            write(cmul,*) j-1
            tcmul = trim(adjustl(cmul))
            write(luout,*) 'M'//tcmul//'                                     &
                       &                    '
            write(luout,*)'-------------------------------------------------------&
                          &---'
            write(luout,'(a)')  ' Site      Vpe                         E   &
                           &         ' 
            write(luout,'(a)')  '                         x             y   &
                           &        z                      '
            write(luout,'(a)') '=======================================================&
                                &==='     
            do i = 1, ncoords
                write(luout,'(i4,2x,(4e13.5))') i, Vpe(j-1,i), Fpe(j:j+2,i)     
            end do

            write(luout,'(a)') '=======================================================&
                        &==='
            write(luout,*)
        end do    
    end if

end subroutine pe_mappot2points
!------------------------------------------------------------------------------

subroutine pe_read_potential(coords, charges)

    real(dp), dimension(:), intent(in), optional :: charges
    real(dp), dimension(:,:), intent(in), optional :: coords

    integer :: i, j, k, l, s
    integer :: lupot, lumep, nlines
    integer :: nidx, idx, jdx, kdx, ldx
    integer, dimension(:), allocatable :: idxs
    real(dp) :: rclose
    real(dp), dimension(21) :: temp
    character(len=2) :: auoraa
    character(len=80) :: word
    logical :: lexist

#if defined(VAR_MPI)
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, ncores, ierr)
#else
    myid = 0
    ncores = 1
#endif

    allocate(nsites(-1:ncores-1))
    nsites = 0

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
        allocate(npoints(-1:ncores-1))
        npoints = 0
        npoints(0) = xsteps * ysteps * zsteps
        allocate(mepgrid(3,npoints(0)))
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

    inquire(file='POTENTIAL.INP', exist=lexist)
    if (lexist) then
        call openfile('POTENTIAL.INP', lupot, 'old', 'formatted')
    else
        if (pe_savden) then
            return
        else if (pe_pd) then
            goto 101
        else if (pe_mep) then
            return
        else
            stop 'POTENTIAL.INP not found!'
        end if
    end if

    do
        read(lupot,*,end=100) word

        if (trim(word) == 'coordinates') then
            read(lupot,*) nsites(0)
            read(lupot,*) auoraa
            allocate(elems(1,nsites(0)), Zs(1,nsites(0)), Rs(3,nsites(0)))
            do i = 1, nsites(0)
                read(lupot,*) elems(1,i), (Rs(j,i), j = 1, 3)
                Zs(1,i) = elem2charge(elems(1,i))
            end do
        else if (trim(word) == 'monopoles') then
            lmul(0) = .true.
            if (mulorder < 0) mulorder = 0
            allocate(M0s(1,nsites(0)))
            M0s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, temp(1)
                M0s(1,s) = temp(1)
            end do
        else if (trim(word) == 'dipoles') then
            lmul(1) = .true.
            if (mulorder < 1) mulorder = 1
            allocate(M1s(3,nsites(0)))
            M1s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 3)
                M1s(:,s) = temp(1:3)
            end do
        else if (trim(word) == 'quadrupoles') then
            lmul(2) = .true.
            if (mulorder < 2) mulorder = 2
            allocate(M2s(6,nsites(0)))
            M2s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 6)
                M2s(:,s) = temp(1:6)
            end do
        else if (trim(word) == 'octopoles') then
            lmul(3) = .true.
            if (mulorder < 3) mulorder = 3
            allocate(M3s(10,nsites(0)))
            M3s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 10)
                M3s(:,s) = temp(1:10)
            end do
        else if (trim(word) == 'hexadecapoles') then
            lmul(4) = .true.
            if (mulorder < 4) mulorder = 4
            allocate(M4s(15,nsites(0)))
            M4s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 15)
                M4s(:,s) = temp(1:15)
            end do
        else if (trim(word) == 'ditriacontapoles') then
            lmul(5) = .true.
            if (mulorder < 5) mulorder = 5
            allocate(M5s(21,nsites(0)))
            M5s = 0.0d0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 21)
                M5s(:,s) = temp(1:21)
            end do
        else if (trim(word) == 'isoalphas') then
            lpol(1) = .true.
            pe_polar = .true.
            if (.not. allocated(P1s)) then
                allocate(P1s(6,nsites(0)))
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
            pe_polar = .true.
            if (.not. allocated(P1s)) then
                allocate(P1s(6,nsites(0)))
                P1s = 0.0d0
            end if
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 6)
                P1s(:,s) = temp(1:6)
            end do
        else if (trim(word) == 'exclists' .or. trim(word) == 'exlists') then
            read(lupot,*) lexlst
            allocate(exclists(lexlst,nsites(0)))
            do i = 1, nsites(0)
                read(lupot,*) (exclists(j,i), j = 1, lexlst)
            end do
        else if (word(1:1) == '!' .or. word(1:1) == '#') then
            cycle
        end if
    end do

100 continue

    close(lupot)

    ! if coordinates are in AA then convert to AU
    if (auoraa == 'AA') then
        Rs = Rs * aa2au
    end if

101 write(luout,'(//2x,a)') 'Polarizable Embedding potential'
    write(luout,'(2x,a)')   '-------------------------------'
    if (nsites(0) > 0) then
        write(luout,'(/4x,a,i6)') 'Number of classical sites: ', nsites(0)
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
    if (lpol(1)) then
        write(luout,'(/4x,a)') '(An)isotropic dipole-dipole polarizabilities.'
        if (pe_damp) then
            write(luout,'(/4x,a)') 'Induced dipole-induced dipole&
                                   & interactions will be damped'
            write(luout,'(4x,a,f8.4)') 'using damping coefficient:', damp
        end if
        if (pe_gspol) then
            write(luout,'(/4x,a)') 'Dynamic response from environment will be&
                                   & neglected during response calculation.'
        end if
        if (pe_nomb) then
            write(luout,'(/4x,a)') 'Many-body interactions will be neglected.'
        end if
        if (pe_iter) then
            write(luout,'(/4x,a)') 'Iterative solver for induced dipoles will&
                                   & be used'
            write(luout,'(4x,a,es7.1)') 'with convergence threshold: ', thriter
        else
            write(luout,'(/4x,a)') 'Direct solver for induced dipoles will be&
                                  & used.'
        end if
    end if
    if (pe_pd) then
        write(luout,'(/4x,a,i4)') 'Number of polarizable densities: ', npds
        if (pe_repuls) then
            write(luout,'(/4x,a,f5.3)') 'Repulsion operator will be used for&
                                        & PDs using the scaling factor: ',&
                                        & rep_factor
        end if
        if (pe_gauss) then
            write(luout,'(/4x,a,f5.3)') 'Gaussian nuclear charges will be&
                                        & used for PDs using a scaling&
                                        & factor: ', gauss_factor
        end if
    end if

   ! default exclusion list (everything polarizes everything)
    if (.not. allocated(exclists)) then
        lexlst = 1
        allocate(exclists(lexlst,nsites(0)))
        do i = 1, nsites(0)
            exclists(1,i) = i
        end do
    end if

    ! handling sites near quantum-classical border
    ! -----------------------------------------------
    if (pe_border) then
        ! first locate all sites within given threshold of QM nuclei
        allocate(idxs(nsites(0)))
        idxs = 0; nidx = 0
        do i = 1, qmnucs
            do j = 1, nsites(0)
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
            do i = 1, nidx
                write(luout,'(/4x,a,i6,2x,a)') 'Removing parameters on site:',&
                                               & idxs(i), elems(1,idxs(i))
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
            do i = 1, nidx
                rclose = 1.0d10
                do j = 1, nsites(0)
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
                if (lmul(1)) then
                    M1s(:,idx) = M1s(:,idx) + M1s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(2)) then
                    M2s(:,idx) = M2s(:,idx) + M2s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(3)) then
                    M3s(:,idx) = M3s(:,idx) + M3s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(4)) then
                    M4s(:,idx) = M4s(:,idx) + M4s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(5)) then
                    M5s(:,idx) = M5s(:,idx) + M5s(:,idxs(i)) / 3.0d0
                endif
                if (lpol(1)) then
                    P1s(:,idx) = P1s(:,idx) + P1s(:,idxs(i)) / 3.0d0
                end if

                rclose = 1.0d10
                do j = 1, nsites(0)
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
                if (lmul(1)) then
                    M1s(:,jdx) = M1s(:,jdx) + M1s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(2)) then
                    M2s(:,jdx) = M2s(:,jdx) + M2s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(3)) then
                    M3s(:,jdx) = M3s(:,jdx) + M3s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(4)) then
                    M4s(:,jdx) = M4s(:,jdx) + M4s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(5)) then
                    M5s(:,jdx) = M5s(:,jdx) + M5s(:,idxs(i)) / 3.0d0
                endif
                if (lpol(1)) then
                    P1s(:,jdx) = P1s(:,jdx) + P1s(:,idxs(i)) / 3.0d0
                end if

                rclose = 1.0d10
                do j = 1, nsites(0)
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
                if (lmul(1)) then
                    M1s(:,kdx) = M1s(:,kdx) + M1s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(2)) then
                    M2s(:,kdx) = M2s(:,kdx) + M2s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(3)) then
                    M3s(:,kdx) = M3s(:,kdx) + M3s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(4)) then
                    M4s(:,kdx) = M4s(:,kdx) + M4s(:,idxs(i)) / 3.0d0
                endif
                if (lmul(5)) then
                    M5s(:,kdx) = M5s(:,kdx) + M5s(:,idxs(i)) / 3.0d0
                endif
                if (lpol(1)) then
                    P1s(:,kdx) = P1s(:,kdx) + P1s(:,idxs(i)) / 3.0d0
                end if

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

                write(luout,'(/4x,a,i6)') 'Redistributing parameters on&
                                          & site:', idxs(i)
                write(luout,'(4x,a,3i6)') 'to neighbouring sites:', idx, jdx, kdx

                if (lmul(0)) then
                    write(luout,'(/4x,a)') 'Resulting monopoles:'
                    write(luout,'(4x,a)') '--------------------'
                    write(luout,'(4x,i6,2x,f9.4)') idx, M0s(:,idx)
                    write(luout,'(4x,i6,2x,f9.4)') jdx, M0s(:,jdx)
                    write(luout,'(4x,i6,2x,f9.4)') kdx, M0s(:,kdx)
                end if
                if (lmul(1)) then
                    write(luout,'(/4x,a)') 'Resulting dipoles:'
                    write(luout,'(4x,a)') '------------------'
                    write(luout,'(4x,i6,2x,3f9.4)') idx, M1s(:,idx)
                    write(luout,'(4x,i6,2x,3f9.4)') jdx, M1s(:,jdx)
                    write(luout,'(4x,i6,2x,3f9.4)') kdx, M1s(:,kdx)
                end if
                if (lmul(2)) then
                    write(luout,'(/4x,a)') 'Resulting quadrupoles:'
                    write(luout,'(4x,a)') '----------------------'
                    write(luout,'(4x,i6,2x,6f9.4)') idx, M2s(:,idx)
                    write(luout,'(4x,i6,2x,6f9.4)') jdx, M2s(:,jdx)
                    write(luout,'(4x,i6,2x,6f9.4)') kdx, M2s(:,kdx)
                end if
                if (lmul(3)) then
                    write(luout,'(/4x,a)') 'Resulting octopoles:'
                    write(luout,'(4x,a)') '--------------------'
                    write(luout,'(4x,i6,2x,10f9.4)') idx, M3s(:,idx)
                    write(luout,'(4x,i6,2x,10f9.4)') jdx, M3s(:,jdx)
                    write(luout,'(4x,i6,2x,10f9.4)') kdx, M3s(:,kdx)
                end if
                if (lmul(4)) then
                    write(luout,'(/4x,a)') 'Resulting hexadecapoles:'
                    write(luout,'(4x,a)') '------------------------'
                    write(luout,'(4x,i6,2x,15f9.4)') idx, M4s(:,idx)
                    write(luout,'(4x,i6,2x,15f9.4)') jdx, M4s(:,jdx)
                    write(luout,'(4x,i6,2x,15f9.4)') kdx, M4s(:,kdx)
                end if
                if (lmul(5)) then
                    write(luout,'(/4x,a)') 'Resulting ditriacontapoles:'
                    write(luout,'(4x,a)') '---------------------------'
                    write(luout,'(4x,i6,2x,21f9.4)') idx, M5s(:,idx)
                    write(luout,'(4x,i6,2x,21f9.4)') jdx, M5s(:,jdx)
                    write(luout,'(4x,i6,2x,21f9.4)') kdx, M5s(:,kdx)
                end if
                if (lpol(1)) then
                    write(luout,'(/4x,a)') 'Resulting polarizabilities:'
                    write(luout,'(4x,a)') '---------------------------'
                    write(luout,'(4x,i6,2x,6f9.4)') idx, P1s(:,idx)
                    write(luout,'(4x,i6,2x,6f9.4)') jdx, P1s(:,jdx)
                    write(luout,'(4x,i6,2x,6f9.4)') kdx, P1s(:,kdx)
                end if
            end do
        end if
    end if

    ! number of polarizabilities different from zero
    if (lpol(1)) then
        allocate(zeroalphas(nsites(0)))
        do i = 1, nsites(0)
            if (abs(maxval(P1s(:,i))) <= zero) then
                zeroalphas(i) = .true.
            else
                zeroalphas(i) = .false.
                npols = npols + 1
            end if
        end do
    end if

end subroutine pe_read_potential

!------------------------------------------------------------------------------

subroutine pe_master(runtype, denmats, fckmats, nmats, Epe, dalwrk)

    character(*), intent(in) :: runtype
    integer, intent(in) :: nmats
    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out), optional :: fckmats
    real(dp), dimension(:), intent(out), optional :: Epe
    real(dp), dimension(:), target, intent(inout) :: dalwrk

#if defined(VAR_MPI)
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, ncores, ierr)
#else
    myid = 0
    ncores = 1
#endif

    work => dalwrk

    ! determine what to calculate and do consistency check
    if (runtype == 'fock') then
        fock = .true.
        energy = .false.
        response = .false.
        mep = .false.
        scfcycle = scfcycle + 1
        if (.not. present(fckmats)) then
            stop 'Output matrices are missing from input!'
        else if (.not. present(Epe)) then
            stop 'The energy variable is missing from input!'
        end if
    else if (runtype == 'energy') then
        fock = .false.
        energy = .true.
        response = .false.
        mep = .false.
    else if (runtype == 'response') then
        if (pe_gspol) return
        if (npols < 1) return
        if (.not. lpol(1)) return
        fock = .false.
        energy = .false.
        response = .true.
        mep = .false.
        if (.not. present(fckmats)) then
            stop 'Output matrices are missing from input!'
        end if
    else if (runtype == 'mep') then
        fock = .false.
        energy = .false.
        response = .false.
        mep = .true.
    else
        stop 'Could not determine calculation type.'
    end if

    ndens = nmats
    nnbas = size(denmats) / ndens

#if defined(VAR_MPI)
    if (myid == 0 .and. ncores > 1) then
        call mpi_bcast(44, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        if (fock) then
            call mpi_bcast(1, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        else if (energy) then
            call mpi_bcast(2, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        else if (response) then
            call mpi_bcast(3, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        else if (mep) then
            call mpi_bcast(4, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        end if

        call mpi_bcast(nnbas, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        call mpi_bcast(ndens, 1, MPI_INTEGER, myid, MPI_COMM_WORLD, ierr)
        call mpi_bcast(denmats, nnbas*ndens, MPI_REAL8,&
                      &myid, MPI_COMM_WORLD, ierr)

        if (.not. initialized) then
            call pe_sync()
        end if
    end if
#endif

    if (fock) then
        call pe_fock(denmats, fckmats, Epe)
    else if (energy) then
        call pe_fock(denmats)
        if (pe_infld) then
            if (allocated(crds)) then
                call pe_mappot2points(crds)
            else 
                call pe_mappot2points()
            end if
        endif
! call pe_pepot_point
    else if (response) then
        call pe_polarization(denmats, fckmats)
    else if (mep) then
        if (ndens > 1) stop 'Not implemented for more than 1 density matrix'
        call pe_compute_mep(denmats)
    end if

#if defined(VAR_MPI)
    if (myid == 0 .and. ncores > 1) then
        if (fock .or. response) then
            call mpi_reduce(MPI_IN_PLACE, fckmats, ndens*nnbas, MPI_REAL8,&
                           &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        end if
    end if
#endif

    nullify(work)

end subroutine pe_master

!------------------------------------------------------------------------------

#if defined(VAR_MPI)
subroutine pe_mpi(dalwrk, runtype)

    real(dp), dimension(:), target, intent(inout) :: dalwrk
    integer :: runtype

    integer :: i
    integer :: nwrk

    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, ncores, ierr)

    work => dalwrk

    nwrk = size(work)

    if (runtype == 1) then
        fock = .true.
        energy = .false.
        response = .false.
        mep = .false.
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
    end if

    call mpi_bcast(nnbas, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ndens, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(work(1), nnbas*ndens, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    if (.not. initialized) then
        call pe_sync()
    end if

    if (fock) then
        call pe_fock(work(1:ndens*nnbas), work(ndens*nnbas+1:2*ndens*nnbas),&
                    &work(2*ndens*nnbas+1:2*ndens*nnbas+ndens))
    else if (energy) then
        call pe_fock(work(1:ndens*nnbas))
    else if (response) then
        call pe_polarization(work(1:ndens*nnbas),&
                            &work(ndens*nnbas+1:2*ndens*nnbas))
    else if (mep) then
        call pe_compute_mep(work(1:ndens*nnbas))
    end if

    if (fock .or. response) then
        call mpi_reduce(work(ndens*nnbas+1), 0, ndens*nnbas, MPI_REAL8,&
                       &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if

    nullify(work)

end subroutine pe_mpi

!------------------------------------------------------------------------------

subroutine pe_sync()

    integer :: i, j, k
    integer :: ndist, nrest

    allocate(ndists(0:ncores-1))

    if (myid == 0) then
        ndist = nsites(0) / ncores
        ndists = ndist
        if (ncores * ndist < nsites(0)) then
            nrest = nsites(0) - ncores * ndist
            do i = 0, nrest-1
                ndists(i) = ndists(i) + 1
            end do
        end if
        do i = 1, ncores-1
            nsites(i) = sum(ndists(0:i))
        end do
        nsites(0) = ndists(0)
    else if (myid /= 0) then
        allocate(nsites(-1:ncores-1))
    end if

    call mpi_bcast(ndists, ncores, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nsites, ncores+1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    call mpi_bcast(qmnucs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if (myid /= 0) allocate(Zm(1,qmnucs))
    call mpi_bcast(Zm, qmnucs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    if (myid /= 0) allocate(Rm(3,qmnucs))
    call mpi_bcast(Rm, 3*qmnucs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    if (myid == 0) then
        allocate(displs(0:ncores-1))
        displs(0) = 0
        do i = 1, ncores-1
            displs(i) = displs(i-1) + 3 * ndists(i-1)
        end do
    else if (myid /= 0) then
        allocate(Rs(3,sum(ndists)))
    end if

    call mpi_bcast(Rs, 3*nsites(ncores-1), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    call mpi_bcast(lpol, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    if (lpol(1)) then
        allocate(npoldists(0:ncores-1))
        if (myid == 0) then
            npoldists = 0
            do i = 0, ncores-1
                do j = nsites(i-1)+1, nsites(i)
                    if (zeroalphas(j)) then
                        continue
                    else
                        npoldists(i) = npoldists(i) + 1
                    end if
                end do
            end do
        end if
        call mpi_bcast(npoldists, ncores, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(npols, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(lexlst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (myid /= 0) allocate(exclists(lexlst,nsites(ncores-1)))
        call mpi_bcast(exclists, lexlst*nsites(ncores-1), MPI_INTEGER,&
                      &0, MPI_COMM_WORLD, ierr)
        if (myid /= 0) allocate(zeroalphas(nsites(ncores-1)))
        call mpi_bcast(zeroalphas, nsites(ncores-1), MPI_LOGICAL,&
                      &0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(pe_pd, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(pe_nomb, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(pe_iter, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(pe_damp, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(damp, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        if (pe_iter) then
            if (myid /= 0) allocate(P1s(6,nsites(ncores-1)))
            call mpi_bcast(P1s, 6*nsites(ncores-1), MPI_REAL8,&
                          &0, MPI_COMM_WORLD, ierr)
        end if
    end if

    call mpi_bcast(mulorder, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(lmul, 6, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    if (mep) then
        if (lmul(0)) then
            if (myid /= 0) allocate(M0s(1,nsites(ncores-1)))
            call mpi_bcast(M0s, nsites(ncores-1), MPI_REAL8,&
                           0, MPI_COMM_WORLD, ierr)
        end if
        if (lmul(1)) then
            if (myid /= 0) allocate(M1s(3,nsites(ncores-1)))
            call mpi_bcast(M1s, 3*nsites(ncores-1), MPI_REAL8,&
                           0, MPI_COMM_WORLD, ierr)
        end if
        if (lmul(2)) then
            if (myid /= 0) allocate(M2s(6,nsites(ncores-1)))
            call mpi_bcast(M2s, 6*nsites(ncores-1), MPI_REAL8,&
                           0, MPI_COMM_WORLD, ierr)
        end if
        if (lmul(3)) then
            if (myid /= 0) allocate(M3s(10,nsites(ncores-1)))
            call mpi_bcast(M3s, 10*nsites(ncores-1), MPI_REAL8,&
                           0, MPI_COMM_WORLD, ierr)
        end if
        if (lmul(4)) then
            if (myid /= 0) allocate(M4s(15,nsites(ncores-1)))
            call mpi_bcast(M4s, 15*nsites(ncores-1), MPI_REAL8,&
                           0, MPI_COMM_WORLD, ierr)
        end if
        if (lmul(5)) then
            if (myid /= 0) allocate(M5s(21,nsites(ncores-1)))
            call mpi_bcast(M5s, 21*nsites(ncores-1), MPI_REAL8,&
                           0, MPI_COMM_WORLD, ierr)
        end if
    else
        if (lmul(0)) then
            if (myid == 0) then
                displs(0) = 0
                do i = 1, ncores-1
                    displs(i) = displs(i-1) + ndists(i-1)
                end do
                call mpi_scatterv(M0s, ndists, displs, MPI_REAL8,&
                                 &MPI_IN_PLACE, 0, MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            else if (myid /= 0) then
                allocate(M0s(1,ndists(myid)))
                call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                                 &M0s, ndists(myid), MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            end if
        end if
    
        if (lmul(1)) then
            if (myid == 0) then
                displs(0) = 0
                do i = 1, ncores-1
                    displs(i) = displs(i-1) + 3 * ndists(i-1)
                end do
                call mpi_scatterv(M1s, 3*ndists, displs, MPI_REAL8,&
                                 &MPI_IN_PLACE, 0, MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            else if (myid /= 0) then
                allocate(M1s(3,ndists(myid)))
                call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                                 &M1s, 3*ndists(myid), MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            end if
        end if
    
        if (lmul(2)) then
            if (myid == 0) then
                displs(0) = 0
                do i = 1, ncores-1
                    displs(i) = displs(i-1) + 6 * ndists(i-1)
                end do
                call mpi_scatterv(M2s, 6*ndists, displs, MPI_REAL8,&
                                 &MPI_IN_PLACE, 0, MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            else if (myid /= 0) then
                allocate(M2s(6,ndists(myid)))
                call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                                 &M2s, 6*ndists(myid), MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            end if
        end if
    
        if (lmul(3)) then
            if (myid == 0) then
                displs(0) = 0
                do i = 1, ncores-1
                    displs(i) = displs(i-1) + 10 * ndists(i-1)
                end do
                call mpi_scatterv(M3s, 10*ndists, displs, MPI_REAL8,&
                                 &MPI_IN_PLACE, 0, MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            else if (myid /= 0) then
                allocate(M3s(10,ndists(myid)))
                call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                                 &M3s, 10*ndists(myid), MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            end if
        end if
    
        if (lmul(4)) then
            if (myid == 0) then
                displs(0) = 0
                do i = 1, ncores-1
                    displs(i) = displs(i-1) + 15 * ndists(i-1)
                end do
                call mpi_scatterv(M4s, 15*ndists, displs, MPI_REAL8,&
                                 &MPI_IN_PLACE, 0, MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            else if (myid /= 0) then
                allocate(M4s(15,ndists(myid)))
                call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                                 &M4s, 15*ndists(myid), MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            end if
        end if
    
        if (lmul(5)) then
            if (myid == 0) then
                displs(0) = 0
                do i = 1, ncores-1
                    displs(i) = displs(i-1) + 21 * ndists(i-1)
                end do
                call mpi_scatterv(M5s, 21*ndists, displs, MPI_REAL8,&
                                 &MPI_IN_PLACE, 0, MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            else if (myid /= 0) then
                allocate(M5s(21,ndists(myid)))
                call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                                 &M5s, 21*ndists(myid), MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            end if
        end if
    end if

    if (mep) then
        call mpi_bcast(mep_field, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(mep_fldnrm, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(mep_extfld, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(mep_qmcube, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(mep_mulcube, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        allocate(nmepdists(0:ncores-1))
        if (myid == 0) then
            ndist = npoints(0) / ncores
            nmepdists = ndist
            if (ncores * ndist < npoints(0)) then
                nrest = npoints(0) - ncores * ndist
                do i = 0, nrest-1
                    nmepdists(i) = nmepdists(i) + 1
                end do
            end if
            do i = 1, ncores-1
                npoints(i) = sum(nmepdists(0:i))
            end do
            npoints(0) = nmepdists(0)
        else if (myid /= 0) then
            allocate(npoints(-1:ncores-1))
        end if
        call mpi_bcast(nmepdists, ncores, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(npoints, ncores+1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (myid == 0) then
            displs(0) = 0
            do i = 1, ncores-1
                displs(i) = displs(i-1) + 3 * nmepdists(i-1)
            end do
            call mpi_scatterv(mepgrid, 3*nmepdists, displs, MPI_REAL8,&
                             &MPI_IN_PLACE, 0, MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        else if (myid /= 0) then
            allocate(mepgrid(3,nmepdists(myid)))
            call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                             &mepgrid, 3*nmepdists(myid), MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        end if
    end if

    initialized = .true.

end subroutine pe_sync
#endif

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
        allocate(Vqm(1,npoints(ncores-1)))
        allocate(Tk_ints(nnbas,1))
        i = 1
        do point = npoints(myid-1)+1, npoints(myid)
            call Tk_integrals(Tk_ints(:,1), nnbas, 1, mepgrid(:,i), .false., 0.0d0)
            Vqm(1,i) = dot(denmats, Tk_ints(:,1))
            do j = 1, qmnucs
                call Tk_tensor(Tm(1:1), mepgrid(:,i) - Rm(:,j))
                Vqm(1,i) = Vqm(1,i) + Zm(1,j) * Tm(1)
            end do
            i = i + 1
        end do
        deallocate(Tk_ints)
#if defined(VAR_MPI)
        if (myid == 0 .and. ncores > 1) then
            displs(0) = 0
            do i = 1, ncores-1
                displs(i) = displs(i-1) + nmepdists(i-1)
            end do
            call mpi_gatherv(MPI_IN_PLACE, 0, MPI_REAL8,&
                            &Vqm, nmepdists, displs, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        else if (myid /= 0) then
            call mpi_gatherv(Vqm, nmepdists(myid), MPI_REAL8,&
                            &0, 0, 0, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        end if
#endif
        if (myid == 0) then
            call openfile('qm_mep.cube', lu, 'new', 'formatted')
            write(lu,'(a)') 'QM MEP'
            write(lu,'(a)') 'Generated by the Polarizable Embedding module'
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
        allocate(Vpe(0:mulorder,npoints(ncores-1)))
        Vpe = 0.0d0
        i = 1
        do point = npoints(myid-1)+1, npoints(myid)
            do j = 1, nsites(ncores-1)
                Rsp = mepgrid(:,i) - Rs(:,j)
                if (lmul(0)) then
                    allocate(Tsp(1), factors(1))
                    call symmetry_factors(factors)
                    taylor = 1.0d0 / factorial(0)
                    call Tk_tensor(Tsp, Rsp)
                    Vpe(0,i) = Vpe(0,i) + taylor * factors(1) * Tsp(1) * M0s(1,j)
                    deallocate(Tsp, factors)
                end if
                if (lmul(1)) then
                    allocate(Tsp(3), factors(3))
                    call symmetry_factors(factors)
                    taylor = - 1.0d0 / factorial(1)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 3
                        Vpe(1,i) = Vpe(1,i) + taylor * factors(k) * Tsp(k) * M1s(k,j)
                    end do
                    deallocate(Tsp, factors)
                end if
                if (lmul(2)) then
                    allocate(Tsp(6), factors(6))
                    call symmetry_factors(factors)
                    taylor = 1.0d0 / factorial(2)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 6
                        Vpe(2,i) = Vpe(2,i) + taylor * factors(k) * Tsp(k) * M2s(k,j)
                    end do
                    deallocate(Tsp, factors)
                end if
                if (lmul(3)) then
                    allocate(Tsp(10), factors(10))
                    call symmetry_factors(factors)
                    taylor = - 1.0d0 / factorial(3)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 10
                        Vpe(3,i) = Vpe(3,i) + taylor * factors(k) * Tsp(k) * M3s(k,j)
                    end do
                    deallocate(Tsp, factors)
                end if
                if (lmul(4)) then
                    allocate(Tsp(15), factors(15))
                    call symmetry_factors(factors)
                    taylor = 1.0d0 / factorial(4)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 15
                        Vpe(4,i) = Vpe(4,i) + taylor * factors(k) * Tsp(k) * M4s(k,j)
                    end do
                    deallocate(Tsp, factors)
                end if
                if (lmul(5)) then
                    allocate(Tsp(21), factors(21))
                    call symmetry_factors(factors)
                    taylor = - 1.0d0 / factorial(5)
                    call Tk_tensor(Tsp, Rsp)
                    do k = 1, 21
                        Vpe(5,i) = Vpe(5,i) + taylor * factors(k) * Tsp(k) * M5s(k,j)
                    end do
                    deallocate(Tsp, factors)
                end if
            end do
            i = i + 1
        end do
#if defined(VAR_MPI)
        if (myid == 0 .and. ncores > 1) then
            displs(0) = 0
            do i = 1, ncores-1
                displs(i) = displs(i-1) + (mulorder + 1) * nmepdists(i-1)
            end do
            call mpi_gatherv(MPI_IN_PLACE, 0, MPI_REAL8,&
                            &Vpe, (mulorder+1)*nmepdists, displs, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        else if (myid /= 0) then
            call mpi_gatherv(Vpe, (mulorder+1)*nmepdists(myid), MPI_REAL8,&
                            &0, 0, 0, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
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
                write(lu,'(a)') 'Generated by the Polarizable Embedding module'
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
                    write(lum2,'(6e13.5)') (Vpe(0,j:k) + Vpe(1,j:k) + Vpe(2,j:k))
                end if
                if (mulorder >= 3) then
                    write(lum3,'(6e13.5)') (Vpe(0,j:k) + Vpe(1,j:k) + Vpe(2,j:k)&
                                            + Vpe(3,j:k))
                end if
                if (mulorder >= 4) then
                    write(lum4,'(6e13.5)') (Vpe(0,j:k) + Vpe(1,j:k) + Vpe(2,j:k)&
                                            + Vpe(3,j:k) + Vpe(4,j:k))
                end if
                if (mulorder >= 5) then
                    write(lum5,'(6e13.5)') (Vpe(0,j:k) + Vpe(1,j:k) + Vpe(2,j:k)&
                                            + Vpe(3,j:k) + Vpe(4,j:k) + Vpe(5,j:k))
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
        do i = 1, nsites(ncores-1)
            if (zeroalphas(i)) cycle
            do j = 1, nsites(ncores-1)
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
        call induced_dipoles(M1inds, Fmuls)
        if (myid == 0) then
            write(luout,*) 'Induced dipole moments:'
            write(luout,'(3f12.6)') M1inds
        end if
        deallocate(Fmuls)
#if defined(VAR_MPI)
        if (ncores > 1) then
            call mpi_bcast(M1inds, 3*npols, MPI_REAL8,&
                          &0, MPI_COMM_WORLD, ierr)
        end if
#endif
        allocate(Vind(1,npoints(ncores-1)))
        Vind = 0.0d0
        allocate(Tsp(3), factors(3))
        call symmetry_factors(factors)
        taylor = - 1.0d0 / factorial(1)
        i = 1
        do point = npoints(myid-1)+1, npoints(myid)
            l = 0
            do j = 1, nsites(ncores-1)
                if (zeroalphas(j)) cycle
                Rsp = mepgrid(:,i) - Rs(:,j)
                Tsp = 0.0d0
                call Tk_tensor(Tsp, Rsp)
                do k = 1, 3
                    Vind(1,i) = Vind(1,i) +&
                                taylor * factors(k) * Tsp(k) * M1inds(l+k,1)
                end do
                l = l + 3
            end do
            i = i + 1
        end do
        deallocate(Tsp, factors)
        deallocate(M1inds)
#if defined(VAR_MPI)
        if (myid == 0 .and. ncores > 1) then
            displs(0) = 0
            do i = 1, ncores-1
                displs(i) = displs(i-1) + nmepdists(i-1)
            end do
            call mpi_gatherv(MPI_IN_PLACE, 0, MPI_REAL8,&
                            &Vind, nmepdists, displs, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        else if (myid /= 0) then
            call mpi_gatherv(Vind, nmepdists(myid), MPI_REAL8,&
                            &0, 0, 0, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        end if
#endif
        if (myid == 0) then
            call openfile('ind_mep.cube', lu, 'new', 'formatted')
            write(lu,'(a)') 'PE induced potential'
            write(lu,'(a)') 'Generated by the Polarizable Embedding module'
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
                write(lu,'(6e13.5)') Vind(1,j:k)
            end do
            close(lu)
        end if
        deallocate(Vind)
    end if

    if (mep_field) then
        if (mep_qmcube) then
            allocate(Fqm(3,npoints(ncores-1)))
            allocate(Tk_ints(nnbas,3))
            i = 1
            do point = npoints(myid-1)+1, npoints(myid)
                call Tk_integrals(Tk_ints, nnbas, 3, mepgrid(:,i), .false., 0.0d0)
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
            if (myid == 0 .and. ncores > 1) then
                displs(0) = 0
                do i = 1, ncores-1
                    displs(i) = displs(i-1) + 3 * nmepdists(i-1)
                end do
                call mpi_gatherv(MPI_IN_PLACE, 0, MPI_REAL8,&
                                &Fqm, 3*nmepdists, displs, MPI_REAL8,&
                                &0, MPI_COMM_WORLD, ierr)
            else if (myid /= 0) then
                call mpi_gatherv(Fqm, 3*nmepdists(myid), MPI_REAL8,&
                                &0, 0, 0, MPI_REAL8,&
                                &0, MPI_COMM_WORLD, ierr)
            end if
#endif
            if (myid == 0) then
                if (mep_fldnrm) then
                    call openfile('qm_field.cube', lu, 'new', 'formatted')
                    write(lu,'(a)') 'QM electric field norm'
                    write(lu,'(a)') 'Generated by the Polarizable Embedding module'
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
                        write(lu,'(6e13.5)') (nrm2(Fqm(:,l)), l = j, k)
                    end do
                    close(lu)
                else
                    do l = 1, 3
                        write(cl,*) l
                        tcl = trim(adjustl(cl))
                        call openfile('qm_field_'//tcl//'.cube', lu, 'new', 'formatted')
                        write(lu,'(a)') 'QM electric field component '//tcl
                        write(lu,'(a)') 'Generated by the Polarizable Embedding module'
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
                            write(lu,'(6e13.5)') Fqm(l,j:k)
                        end do
                        close(lu)
                    end do
                end if
            end if
            deallocate(Fqm)
        end if

        if (mulorder >= 0 .and. mep_mulcube) then
            allocate(Fpe(3,0:mulorder,npoints(ncores-1)))
            Fpe = 0.0d0
            i = 1
            do point = npoints(myid-1)+1, npoints(myid)
                do j = 1, nsites(ncores-1)
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
            if (myid == 0 .and. ncores > 1) then
                displs(0) = 0
                do i = 1, ncores-1
                    displs(i) = displs(i-1) + 3 * (mulorder + 1) * nmepdists(i-1)
                end do
                call mpi_gatherv(MPI_IN_PLACE, 0, MPI_REAL8,&
                                &Fpe, 3*(mulorder+1)*nmepdists, displs, MPI_REAL8,&
                                &0, MPI_COMM_WORLD, ierr)
            else if (myid /= 0) then
                call mpi_gatherv(Fpe, 3*(mulorder+1)*nmepdists(myid), MPI_REAL8,&
                                &0, 0, 0, MPI_REAL8,&
                                &0, MPI_COMM_WORLD, ierr)
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
                        write(lu,'(a)') 'Generated by the Polarizable Embedding module'
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
                            write(lum0,'(6e13.5)') (nrm2(Fpe(:,0,l)), l = j, k)
                        end if
                        if (mulorder >= 1) then
                            write(lum1,'(6e13.5)') (nrm2(Fpe(:,0,l) +&
                                                         Fpe(:,1,l)), l = j, k)
                        end if
                        if (mulorder >= 2) then
                            write(lum2,'(6e13.5)') (nrm2(Fpe(:,0,l) +&
                                                         Fpe(:,1,l) +&
                                                         Fpe(:,2,l)), l = j, k)
                        end if
                        if (mulorder >= 3) then
                            write(lum3,'(6e13.5)') (nrm2(Fpe(:,0,l) +&
                                                         Fpe(:,1,l) +&
                                                         Fpe(:,2,l) +&
                                                         Fpe(:,3,l)), l = j, k)
                        end if
                        if (mulorder >= 4) then
                            write(lum4,'(6e13.5)') (nrm2(Fpe(:,0,l) +&
                                                         Fpe(:,1,l) +&
                                                         Fpe(:,2,l) +&
                                                         Fpe(:,3,l) +&
                                                         Fpe(:,4,l)), l = j, k)
                        end if
                        if (mulorder >= 5) then
                            write(lum5,'(6e13.5)') (nrm2(Fpe(:,0,l) +&
                                                         Fpe(:,1,l) +&
                                                         Fpe(:,2,l) +&
                                                         Fpe(:,3,l) +&
                                                         Fpe(:,4,l) +&
                                                         Fpe(:,5,l)), l = j, k)
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
                                write(lum0,'(6e13.5)') Fpe(l,0,j:k)
                            end if
                            if (mulorder >= 1) then
                                write(lum1,'(6e13.5)') (Fpe(l,0,j:k) + Fpe(l,1,j:k))
                            end if
                            if (mulorder >= 2) then
                                write(lum2,'(6e13.5)') (Fpe(l,0,j:k) + Fpe(l,1,j:k)&
                                                        + Fpe(l,2,j:k))
                            end if
                            if (mulorder >= 3) then
                                write(lum3,'(6e13.5)') (Fpe(l,0,j:k) + Fpe(l,1,j:k)&
                                                        + Fpe(l,2,j:k) + Fpe(l,3,j:k))
                            end if
                            if (mulorder >= 4) then
                                write(lum4,'(6e13.5)') (Fpe(l,0,j:k) + Fpe(l,1,j:k)&
                                                        + Fpe(l,2,j:k) + Fpe(l,3,j:k)&
                                                        + Fpe(l,4,j:k))
                            end if
                            if (mulorder >= 5) then
                                write(lum5,'(6e13.5)') (Fpe(l,0,j:k) + Fpe(l,1,j:k)&
                                                        + Fpe(l,2,j:k) + Fpe(l,3,j:k)&
                                                        + Fpe(l,4,j:k) + Fpe(l,5,j:k))
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
            do i = 1, nsites(ncores-1)
                if (zeroalphas(i)) cycle
                do j = 1, nsites(ncores-1)
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
            call induced_dipoles(M1inds, Fmuls)
            deallocate(Fmuls)
#if defined(VAR_MPI)
            if (ncores > 1) then
                call mpi_bcast(M1inds, 3*npols, MPI_REAL8,&
                              &0, MPI_COMM_WORLD, ierr)
            end if
#endif
            allocate(Find(3,npoints(ncores-1)))
            Find = 0.0d0
            i = 1
            do point = npoints(myid-1)+1, npoints(myid)
                l = 1
                do j = 1, nsites(ncores-1)
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
            if (myid == 0 .and. ncores > 1) then
                displs(0) = 0
                do i = 1, ncores-1
                    displs(i) = displs(i-1) + 3 * nmepdists(i-1)
                end do
                call mpi_gatherv(MPI_IN_PLACE, 0, MPI_REAL8,&
                                &Find, 3*nmepdists, displs, MPI_REAL8,&
                                &0, MPI_COMM_WORLD, ierr)
            else if (myid /= 0) then
                call mpi_gatherv(Find, 3*nmepdists(myid), MPI_REAL8,&
                                &0, 0, 0, MPI_REAL8,&
                                &0, MPI_COMM_WORLD, ierr)
            end if
#endif
            if (myid == 0) then
                if (mep_fldnrm) then
                    call openfile('ind_field.cube', lu, 'new', 'formatted')
                    write(lu,'(a)') 'PE induced electric field norm'
                    write(lu,'(a)') 'Generated by the Polarizable Embedding module'
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
                        write(lu,'(a)') 'Generated by the Polarizable Embedding module'
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
                            write(lu,'(6e13.5)') Find(l,j:k)
                        end do
                        close(lu)
                    end do
                end if
            end if
            deallocate(Find)
        end if
    end if

#if defined(VAR_MPI)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

end subroutine pe_compute_mep

!------------------------------------------------------------------------------

subroutine pe_fock(denmats, fckmats, Epe)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out), optional :: fckmats
    real(dp), dimension(:), intent(out), optional :: Epe

    integer :: i
    logical :: es = .false.
    logical :: pol = .false.

    if ((mulorder >= 0) .or. pe_pd) es = .true.
    if (lpol(1)) pol = .true.

    if (allocated(Ees)) deallocate(Ees)
    if (allocated(Epol)) deallocate(Epol)
    if (allocated(Epd)) deallocate(Epd)
    allocate(Ees(0:5,ndens))
    allocate(Epol(3,ndens))
    allocate(Epd(3,ndens))
    Ees = 0.0d0
    Epol = 0.0d0
    Epd = 0.0d0

    if (fock) fckmats = 0.0d0

    if (fock) then
        if (es) call pe_electrostatic(denmats, fckmats)
        if (pol) call pe_polarization(denmats, fckmats)
    else if (energy) then
        if (es) call pe_electrostatic(denmats)
        if (pol) call pe_polarization(denmats)
    end if

    if (fock) then
        Epe = 0.0d0
        do i = 1, ndens
            Epe(i) = sum(Ees(:,i)) + sum(Epol(:,i)) + sum(Epd(:,i))
        end do
    end if

end subroutine pe_fock

!------------------------------------------------------------------------------

subroutine pe_electrostatic(denmats, fckmats)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats

    logical :: lexist
    integer :: lu
    integer :: i, j, k
    real(dp) :: Enuc, Esave
    real(dp), dimension(ndens) :: Eel
    real(dp), dimension(:), allocatable :: tmpfcks

    if (myid == 0) then
        inquire(file='pe_electrostatics.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (ncores > 1) then
        call mpi_bcast(lexist, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    end if
#endif

    if (lexist .and. fock) then
        if (myid == 0) then
            call openfile('pe_electrostatics.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Esave, fckmats
            close(lu)
            do i = 1, ndens
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                Ees(0,i) = Ees(0,i) + dot(denmats(j:k), fckmats(j:k)) + Esave
            end do
        end if
    else
        Esave = 0.0d0
        if (lmul(0)) then
            if (fock) then
                call es_multipoles(M0s, denmats, Eel, Enuc, fckmats)
            else if (energy) then
                call es_multipoles(M0s, denmats, Eel, Enuc)
            end if
            do i = 1, ndens
                Ees(0,i) = Ees(0,i) + Eel(i) + Enuc
            end do
            Esave = Esave + Enuc
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
            Esave = Esave + Enuc
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
            Esave = Esave + Enuc
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
            Esave = Esave + Enuc
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
            Esave = Esave + Enuc
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
            Esave = Esave + Enuc
        end if
        if (myid == 0) then
            if (pe_pd) then
                if (fock) then
                    call es_polarizable_densities(denmats, Eel, Enuc, fckmats)
                else if (energy) then
                    call es_polarizable_densities(denmats, Eel, Enuc)
                end if
                do i = 1, ndens
                    Epd(1,i) = Epd(1,i) + Eel(i) + Enuc
                end do
                Esave = Esave + Enuc
            end if
        end if
        if (fock) then
#if defined(VAR_MPI)
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
#endif
            if (myid == 0) then
                call openfile('pe_electrostatics.bin', lu, 'new', 'unformatted')
                rewind(lu)
                write(lu) Esave, fckmats
                close(lu)
            end if
            if (myid == 0 .and. ncores > 1) then
                fckmats = tmpfcks
                deallocate(tmpfcks)
            end if
        end if
#if defined(VAR_MPI)
        if (myid == 0 .and. ncores > 1) then
            call mpi_reduce(MPI_IN_PLACE, Ees, 7*ndens, MPI_REAL8, MPI_SUM,&
                           &0, MPI_COMM_WORLD, ierr)
        else if (myid /= 0) then
            call mpi_reduce(Ees, 0, 7*ndens, MPI_REAL8, MPI_SUM,&
                           &0, MPI_COMM_WORLD, ierr)
        end if
#endif
    end if

end subroutine pe_electrostatic

!------------------------------------------------------------------------------

subroutine es_polarizable_densities(denmats, Eel, Enuc, fckmats)

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(out) :: Eel
    real(dp), intent(out) :: Enuc
    real(dp), dimension(:), intent(inout), optional :: fckmats

    integer :: i, j, k, l, m, n, o
    integer :: lufck, lexist, lu
    real(dp) :: gauss = 0.0d0
    real(dp) :: Ene, Enn
    real(dp), dimension(ndens) :: Een, Eee
    real(dp), dimension(1) :: Tfm
    real(dp), dimension(3) :: Rfm
    real(dp), dimension(3*npols) :: temp
    real(dp), dimension(nnbas) :: pd_fckmat, pd_repmat
    real(dp), dimension(nnbas,1) :: Zpd_ints
    character(len=99) :: ci
    character(len=99) :: filename

    Eel = 0.0d0; Enuc = 0.0d0

    do i = 1, npds
        Eee = 0.0d0; Een = 0.0d0; Ene = 0.0d0; Enn = 0.0d0
        write(ci,*) i
        ci = adjustl(ci)
        filename = 'pe_fock_'//trim(ci)//'.bin'
        call openfile(trim(filename), lufck, 'old', 'unformatted')
        rewind(lufck)
        read(lufck) temp
        read(lufck) Ene
        read(lufck) pd_fckmat
        read(lufck) pd_repmat
        read(lufck) pdnucs
        allocate(Rpd(3,pdnucs), Zpd(1,pdnucs))
        read(lufck) Rpd, Zpd
        close(lufck)

        do j = 1, ndens
            l = (j - 1) * nnbas + 1
            m = j * nnbas
            if (fock) fckmats(l:m) = fckmats(l:m) + pd_fckmat
            Eee(j) = dot(denmats(l:m), pd_fckmat)
            if (pe_repuls) then
                if (fock) fckmats(l:m) = fckmats(l:m) - rep_factor * pd_repmat
                Epd(3,j) = Epd(3,j) + dot(denmats(l:m), - rep_factor * pd_repmat)
            end if
        end do

        do j = 1, pdnucs
!            gauss = (8.0d0  * gauss_factor) /&
!                    ((P1s(1,j) + P1s(4,j) + P1s(6,j)) / 3.0d0)**(2.0d0/3.0d0)
            do k = 1, qmnucs
                Rfm = Rm(:,k) - Rpd(:,j)
                call Tk_tensor(Tfm, Rfm)
                Enn = Enn + Zm(1,k) * Zpd(1,j) * Tfm(1)
            end do
            call Tk_integrals(Zpd_ints, nnbas, 1, Rpd(:,j), .false., gauss) 
            Zpd_ints = Zpd(1,j) * Zpd_ints
!            call Mk_integrals(Zpd_ints, Rpd(:,j), Zpd(:,j))
            do m = 1, ndens
                n = (m - 1) * nnbas + 1
                o = m * nnbas
                Een(m) = Een(m) + dot(denmats(n:o), Zpd_ints(:,1))
                if (fock) fckmats(n:o) = fckmats(n:o) + Zpd_ints(:,1)
            end do
        end do

        deallocate(Rpd, Zpd)

        Enuc = Enuc + Ene + Enn
        do j = 1, ndens
            Eel(j) = Eel(j) + Een(j) + Eee(j)
        end do
    end do

end subroutine es_polarizable_densities

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
    do site = nsites(myid-1)+1, nsites(myid)
        if (abs(maxval(Mks(:,i))) < zero) then
            i = i + 1
            cycle
        end if

        ! nuclei - multipole interaction energy
        do j = 1, qmnucs
            Rsm = Rm(:,j) - Rs(:,site)
            call Tk_tensor(Tsm, Rsm)
            call symmetry_factors(symfacs)
            do l = 1, ncomps
                Enuc = Enuc + taylor * symfacs(l) * Zm(1,j) * Mks(l,i) * Tsm(l)
            end do
        end do

        ! electron - multipole interaction energy
        call Mk_integrals(Mk_ints, Rs(:,site), Mks(:,i))
        do j = 1, ncomps
            do l = 1, ndens
                m = (l - 1) * nnbas + 1
                n = l * nnbas
                Eel(l) = Eel(l) + dot(denmats(m:n), Mk_ints(:,j))
                if (fock) fckmats(m:n) = fckmats(m:n) + Mk_ints(:,j)
            end do
        end do
        i = i + 1
    end do

    deallocate(Tsm, symfacs, Mk_ints)

end subroutine es_multipoles

!------------------------------------------------------------------------------

subroutine pe_polarization(denmats, fckmats)

    external :: Tk_integrals

    real(dp), dimension(:), intent(in) :: denmats
    real(dp), dimension(:), intent(inout), optional :: fckmats

    integer :: site, ndist, nrest
    integer :: i, j, k, l, m
    logical :: skip
    real(dp), dimension(3*npols) :: Fnucs, Fmuls, Fpd
    real(dp), dimension(3*npols,ndens) :: M1inds, Fels, Ftots
    real(dp), dimension(:,:), allocatable :: Fel_ints

    if (response) fckmats = 0.0d0

    if (response) then
        call electron_fields(Fels, denmats)
        call induced_dipoles(M1inds, Fels)
    else
        call electron_fields(Fels, denmats)
        call nuclear_fields(Fnucs)
        call multipole_fields(Fmuls)
        if (myid == 0) then
            if (pe_pd) then
                call polarizable_density_field(Fpd)
            else
                Fpd = 0.0d0
            end if
            do i = 1, ndens
                Ftots(:,i) = Fels(:,i) + Fnucs + Fmuls + Fpd
            end do
        end if
        call induced_dipoles(M1inds, Ftots)
        if (myid == 0) then
            do i = 1, ndens
                Epol(1,i) = - 0.5d0 * dot(M1inds(:,i), Fels(:,i))
                Epol(2,i) = - 0.5d0 * dot(M1inds(:,i), Fnucs)
                Epol(3,i) = - 0.5d0 * dot(M1inds(:,i), Fmuls)
                if (pe_pd) Epd(2,i) = - 0.5d0 * dot(M1inds(:,i), Fpd)
            end do
        end if
    end if

#if defined(VAR_MPI)
    if (myid == 0 .and. ncores > 1) then
        do i = 1, ndens
            displs(0) = 0
            do j = 1, ncores-1
                displs(j) = displs(j-1) + 3 * npoldists(j-1)
            end do
            call mpi_scatterv(M1inds(:,i), 3*npoldists, displs, MPI_REAL8,&
                             &MPI_IN_PLACE, 0, MPI_REAL8,&
                             &myid, MPI_COMM_WORLD, ierr)
        end do
    else if (myid /= 0) then
        do i = 1, ndens
            call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                             &M1inds(:,i), 3*npoldists(myid), MPI_REAL8,&
                             &0, MPI_COMM_WORLD, ierr)
        end do
    end if
#endif

    if (fock .or. response) then
        allocate(Fel_ints(nnbas,3))
        i = 0
        do site = nsites(myid-1)+1, nsites(myid)
            if (zeroalphas(site)) cycle
            call Tk_integrals(Fel_ints, nnbas, 3, Rs(:,site), .false., 0.0d0)
            do j = 1, 3
                do k = 1, ndens
                    l = (k - 1) * nnbas + 1
                    m = k * nnbas
                    fckmats(l:m) = fckmats(l:m) - M1inds(i+j,k) * Fel_ints(:,j)
                end do
            end do
            i = i + 3
        end do
    end if

end subroutine pe_polarization

!------------------------------------------------------------------------------

subroutine induced_dipoles(M1inds, Fs)

    real(dp), dimension(:,:), intent(out) :: M1inds
    real(dp), dimension(:,:), intent(in) :: Fs

    integer :: lu, iter, info
    integer :: i, j, k, l, m, n, o, p, q
    integer, dimension(:), allocatable :: ipiv
    logical :: exclude, lexist
    logical :: converged = .false.
    real(dp) :: fe = 1.0d0
    real(dp) :: ft = 1.0d0
    real(dp) :: R, R3, R5, Rd, ai, aj, norm, redthr
    real(dp), parameter :: d3i = 1.0d0 / 3.0d0
    real(dp), parameter :: d6i = 1.0d0 / 6.0d0
    real(dp), dimension(:), allocatable :: B, T, Rij, Ftmp, M1tmp

    if (pe_iter) then
        if (myid == 0) then
            if (fock .and. scfcycle <= 5) then
                redthr = 10**(5-scfcycle)
            else
                redthr = 1.0d0
            end if
        end if

        if (myid == 0) then
            inquire(file='pe_induced_dipoles.bin', exist=lexist)
        end if

#if defined(VAR_MPI)
        if (ncores > 1) then
            call mpi_bcast(lexist, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        end if
#endif

        if (lexist .and. (fock .or. energy)) then
            if (myid == 0) then
                call openfile('pe_induced_dipoles.bin', lu, 'old', 'unformatted')
                rewind(lu)
                read(lu) M1inds
                close(lu)
            end if
        end if

        allocate(T(6), Rij(3), Ftmp(3), M1tmp(3))
        do n = 1, ndens
            if (.not.lexist .or. response) then
#if defined(VAR_MPI)
                if (myid == 0 .and. ncores > 1) then
                    displs(0) = 0
                    do i = 1, ncores-1
                        displs(i) = displs(i-1) + 3 * npoldists(i-1)
                    end do
                    call mpi_scatterv(Fs(:,n), 3*npoldists, displs, MPI_REAL8,&
                                     &MPI_IN_PLACE, 0, MPI_REAL8,&
                                     &0, MPI_COMM_WORLD, ierr)
                else if (myid /= 0) then
                    call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                                     &Fs(:,n), 3*npoldists(myid), MPI_REAL8,&
                                     &0, MPI_COMM_WORLD, ierr)
                end if
#endif

                l = 1
                do i = nsites(myid-1)+1, nsites(myid)
                    if (zeroalphas(i)) cycle
                    call spmv(P1s(:,i), Fs(l:l+2,n), M1inds(l:l+2,n), 'L')
                    l = l + 3
                end do

#if defined(VAR_MPI)
                if (myid == 0 .and. ncores > 1) then
                    displs(0) = 0
                    do i = 1, ncores-1
                        displs(i) = displs(i-1) + 3 * npoldists(i-1)
                    end do
                    call mpi_gatherv(MPI_IN_PLACE, 0, MPI_REAL8,&
                                    &M1inds(:,n), 3*npoldists, displs, MPI_REAL8,&
                                    &0, MPI_COMM_WORLD, ierr)
                else if (myid /= 0) then
                    call mpi_gatherv(M1inds(:,n), 3*npoldists(myid), MPI_REAL8,&
                                    &0, 0, 0, MPI_REAL8,&
                                    &0, MPI_COMM_WORLD, ierr)
                end if
#endif
            end if

            if (pe_nomb) cycle

#if defined(VAR_MPI)
            if (myid == 0 .and. ncores > 1) then
                displs(0) = 0
                do i = 1, ncores-1
                    displs(i) = displs(i-1) + 3 * npoldists(i-1)
                end do
                call mpi_scatterv(M1inds(:,n), 3*npoldists, displs, MPI_REAL8,&
                                 &MPI_IN_PLACE, 0, MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            else if (myid /= 0) then
                call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                                 &M1inds(:,n), 3*npoldists(myid), MPI_REAL8,&
                                 &0, MPI_COMM_WORLD, ierr)
            end if
#endif

            iter = 1
            do
                norm = 0.0d0
                l = 1
                do i = 1, nsites(ncores-1)
                    if (zeroalphas(i)) cycle
                    if (pe_damp) then
                        ai = (P1s(1,i) + P1s(4,i) + P1s(6,i)) * d3i
                    end if
                    m = 1
                    Ftmp = 0.0d0
                    do j = nsites(myid-1)+1, nsites(myid)
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
                        call spmv(T, M1inds(m:m+2,n), Ftmp, 'L', 1.0d0, 1.0d0)
                        m = m + 3
                    end do

#if defined(VAR_MPI)
                    if (myid == 0 .and. ncores > 1) then
                        call mpi_reduce(MPI_IN_PLACE, Ftmp, 3, MPI_REAL8,&
                                       &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
                    else if (myid /= 0) then
                        call mpi_reduce(Ftmp, 0, 3, MPI_REAL8,&
                                       &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
                    end if
#endif

                    if (myid == 0) then
                        M1tmp = M1inds(l:l+2,n)
                        Ftmp = Ftmp + Fs(l:l+2,n)
                        call spmv(P1s(:,i), Ftmp, M1inds(l:l+2,n), 'L')
                        norm = norm + nrm2(M1inds(l:l+2,n) - M1tmp)
                    end if

#if defined(VAR_MPI)
                    if (myid == 0 .and. ncores > 1) then
                        displs(0) = 0
                        do j = 1, ncores-1
                            displs(j) = displs(j-1) + 3 * npoldists(j-1)
                        end do
                        call mpi_scatterv(M1inds(:,n), 3*npoldists, displs,&
                                         &MPI_REAL8, MPI_IN_PLACE, 0,&
                                         &MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
                    else if (myid /= 0) then
                        call mpi_scatterv(0, 0, 0, MPI_REAL8,&
                                         &M1inds(:,n), 3*npoldists(myid),&
                                         &MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
                    end if
#endif
                    l = l + 3
                end do

                if (myid == 0) then
                    if (norm < thriter) then
                        write (luout,'(4x,a,i2,a)') 'Induced dipole moments&
                                                    & converged in ',&
                                                    & iter, ' iterations.'
                        converged = .true.
                    else if (iter > 50) then
                        write(luout,*) 'ERROR: could not converge induced&
                                       & dipole moments.'
                        stop 'ERROR: could not converge induced dipole&
                             & moments.'
                    else
                        converged = .false.
                        iter = iter + 1
                    end if
                end if

#if defined(VAR_MPI)
                if (ncores > 1) then
                    call mpi_bcast(converged, 1, MPI_LOGICAL, 0,&
                                  &MPI_COMM_WORLD, ierr)
                end if
#endif
                if (converged) exit
            end do
        end do
        if (fock) then
            if (myid == 0) then
                call openfile('pe_induced_dipoles.bin', lu, 'unknown', 'unformatted')
                rewind(lu)
                write(lu) M1inds
                close(lu)
            end if
        end if
        deallocate(T, Rij, Ftmp, M1tmp)
    else
        if (myid == 0) then
            allocate(B(3*npols*(3*npols+1)/2))
            inquire(file='pe_response_matrix.bin', exist=lexist)
            if (lexist) then
                call openfile('pe_response_matrix.bin', lu, 'old', 'unformatted')
                rewind(lu)
                if (chol) then
                    read(lu) B
                else
                    allocate(ipiv(3*npols))
                    read(lu) B, ipiv
                end if
                close(lu)
            else
                call response_matrix(B)
                if (chol) then
                    call pptrf(B, 'L', info)
                    if (info /= 0) then
                        write(luout,*) 'Cholesky factorization failed.&
                                       & Trying regular...'
                        allocate(ipiv(3*npols))
                        call sptrf(B, 'L', ipiv, info)
                        if (info /= 0) then
                            stop 'ERROR: cannot create response matrix.'
                        else
                            chol = .false.
                        end if
                    end if
                else
                    allocate(ipiv(3*npols))
                    call sptrf(B, 'L', ipiv, info)
                    if (info /= 0) then
                        stop 'ERROR: cannot create response matrix.'
                    end if
                end if
                call openfile('pe_response_matrix.bin', lu, 'new', 'unformatted')
                rewind(lu)
                if (chol) then
                    write(lu) B
                else
                    write(lu) B, ipiv
                end if
                close(lu)
            end if
            M1inds = Fs
            if (chol) then
                call pptrs(B, M1inds, 'L')
                deallocate(B)
            else
                call sptrs(B, M1inds, ipiv, 'L')
                deallocate(B, ipiv)
            end if
        end if
    end if

    ! check induced dipoles
    if (myid == 0) then
        do n = 1, ndens
            l = 1
            do i = 1, nsites(ncores-1)
                if (zeroalphas(i)) cycle
                if (nrm2(M1inds(l:l+2,n)) > 1.0d0) then
                    write(luout,'(4x,a,i6)') 'Large induced dipole encountered&
                                             & at site:', i
                    write(luout,'(f10.4)') nrm2(M1inds(l:l+2,n))
                end if
                l = l + 3
            end do
        end do
    end if

end subroutine induced_dipoles

!------------------------------------------------------------------------------

subroutine electron_fields(Fels, denmats)

    external :: Tk_integrals

    real(dp), dimension(:,:), intent(out) :: Fels
    real(dp), dimension(:), intent(in) :: denmats

    logical :: skip
    integer :: site
    integer :: i, j, k, l, m
    real(dp), dimension(nnbas,3) :: Fel_ints

    Fels = 0.0d0

    i = 0
    do site = nsites(myid-1)+1, nsites(myid)
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
        call Tk_integrals(Fel_ints, nnbas, 3, Rs(:,site), .false., 0.0d0)
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
    if (myid == 0 .and. ncores > 1) then
        displs(0) = 0
        do i = 1, ncores-1
            displs(i) = displs(i-1) + 3 * npoldists(i-1)
        end do
        do i = 1, ndens
            call mpi_gatherv(MPI_IN_PLACE, 0, MPI_REAL8,&
                            &Fels(:,i), 3*npoldists, displs, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        end do
    else if (myid /= 0) then
        do i = 1, ndens
            call mpi_gatherv(Fels(:,i), 3*npoldists(myid), MPI_REAL8,&
                            &0, 0, 0, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        end do
    end if
#endif

end subroutine electron_fields

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
    if (ncores > 1) then
        call mpi_bcast(lexist, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
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
        do site = nsites(myid-1)+1, nsites(myid)
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
        if (myid == 0 .and. ncores > 1) then
            displs(0) = 0
            do i = 1, ncores-1
                displs(i) = displs(i-1) + 3 * npoldists(i-1)
            end do
            call mpi_gatherv(MPI_IN_PLACE, 0, MPI_REAL8,&
                            &Fnucs, 3*npoldists, displs, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
        else if (myid /= 0) then
            call mpi_gatherv(Fnucs, 3*npoldists(myid), MPI_REAL8,&
                            &0, 0, 0, MPI_REAL8,&
                            &0, MPI_COMM_WORLD, ierr)
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

subroutine polarizable_density_field(Fpd)

    real(dp), dimension(:), intent(out) :: Fpd

    integer :: i
    integer :: lu
    character(len=99) :: ci
    character(len=80) :: filename
    real(dp), dimension(3*npols) :: Ftmp

    Fpd = 0.0d0

    do i = 1, npds
        Ftmp = 0.0d0
        write(ci,*) i
        ci = adjustl(ci)
        filename = 'pe_fock_'//trim(ci)//'.bin'
        call openfile(trim(filename), lu, 'old', 'unformatted')
        rewind(lu)
        read(lu) Ftmp
        close(lu)
        Fpd = Fpd + Ftmp
    end do

end subroutine polarizable_density_field

!------------------------------------------------------------------------------

subroutine multipole_fields(Fmuls)

    real(dp), dimension(:), intent(out) :: Fmuls

    logical :: exclude, lexist
    integer :: lu
    integer :: i, j, k, l, m
    real(dp), dimension(3) :: Rji

    if (myid == 0) then
        inquire(file='pe_multipole_field.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (ncores > 1) then
        call mpi_bcast(lexist, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
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
        do i = 1, nsites(ncores-1)
            if (zeroalphas(i)) cycle
            k = 1
            do j = nsites(myid-1)+1, nsites(myid)
                if (i == j) then
                    k = k + 1
                    cycle
                end if
                exclude = .false.
                do m = 1, lexlst
                    if (exclists(m,i) == exclists(1,j)) then
                        exclude = .true.
                        k = k + 1
                        exit
                    end if
                end do
                if (exclude) cycle
! TODO: cutoff???
                Rji = Rs(:,i) - Rs(:,j)
                if (lmul(0)) then
                    if (abs(maxval(M0s(:,k))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M0s(:,k))
                    end if
                end if
                if (lmul(1)) then
                    if (abs(maxval(M1s(:,k))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M1s(:,k))
                    end if
                end if
                if (lmul(2)) then
                    if (abs(maxval(M2s(:,k))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M2s(:,k))
                    end if
                end if
                if (lmul(3)) then
                    if (abs(maxval(M3s(:,k))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M3s(:,k))
                    end if
                end if
                if (lmul(4)) then
                    if (abs(maxval(M4s(:,k))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M4s(:,k))
                    end if
                end if
                if (lmul(5)) then
                    if (abs(maxval(M5s(:,k))) >= zero) then
                        call multipole_field(Fmuls(l:l+2), Rji, M5s(:,k))
                    end if
                end if
                k = k + 1
            end do
            l = l + 3
        end do
#if defined(VAR_MPI)
        if (myid == 0 .and. ncores > 1) then
            call mpi_reduce(MPI_IN_PLACE, Fmuls, 3*npols, MPI_REAL8,&
                           &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        else if (myid /= 0) then
            call mpi_reduce(Fmuls, 0, 3*npols, MPI_REAL8,&
                           &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
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
                    Fi(1) = Fi(1) + taylor * symfac(x-1,y,z) * T(Rji,x,y,z) * Mkj(a)
                    a = a + 1
                end if
                if (y /= 0) then
                    Fi(2) = Fi(2) + taylor * symfac(x,y-1,z) * T(Rji,x,y,z) * Mkj(b)
                    b = b + 1
                end if
                if (z /= 0) then
                    Fi(3) = Fi(3) + taylor * symfac(x,y,z-1) * T(Rji,x,y,z) * Mkj(c)
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
    integer :: i, j, k, l, m, n
    integer, dimension(3) :: ipiv
    real(dp), parameter :: d3i = 1.0d0 / 3.0d0
    real(dp), parameter :: d6i = 1.0d0 / 6.0d0
    real(dp) :: fe = 1.0d0
    real(dp) :: ft = 1.0d0
    real(dp) :: Rd, ai, aj
    real(dp) :: R, R3, R5, T
    real(dp), dimension(3) :: Rij
    real(dp), dimension(6) :: P1inv

    B = 0.0d0
    m = 0
    do i = 1, nsites(ncores-1)
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
            do j = i, nsites(ncores-1)
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
                end if
            end do
        end do
    end do

end subroutine response_matrix

!------------------------------------------------------------------------------

subroutine Tk_coefficients

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    integer :: i, j, k, l, m, n

! TODO: if mulorder less than polorder then it will be too small?
!    allocate(Cnij(2*mulorder+3,0:mulorder+1,0:mulorder+1))
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

    call Tk_integrals(Mk_ints, nnbas, ncomps, Rij, .false., 0.0d0)

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

subroutine pe_save_density(denmat, mofckmat, cmo, nbas, nocc, norb,&
                           coords, charges, dalwrk)

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
    real(dp), dimension(:), allocatable :: Fpd
    real(dp), dimension(:,:), allocatable :: Ftmp
    real(dp), dimension(:), allocatable :: Emo

    work => dalwrk

    ndens = 1
    nnbas = nbas * (nbas + 1) / 2

    ! polarizable density nuclear charges and coordinates
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
    if (lpol(1)) then
        allocate(Ftmp(3*npols,1), Fpd(3*npols)); Ftmp = 0.0d0
        call electron_fields(Ftmp, denmat)
        Fpd = Ftmp(:,1)
        call nuclear_fields(Ftmp(:,1))
        Fpd = Fpd + Ftmp(:,1)
        deallocate(Ftmp)
    else
        allocate(Fpd(3*npols))
        Fpd = 0.0d0
    end if

    ! calculate nuclear - electron energy contribution
    allocate(T0_ints(nnbas,1)); T0_ints = 0.0d0
    Ene = 0.0d0
    do i = 1, corenucs
        call Tk_integrals(T0_ints, nnbas, 1, Rc(:,i), .false., 0.0d0)
        T0_ints = Zc(1,i) * T0_ints
        Ene = Ene + dot(denmat, T0_ints(:,1))
    end do
    deallocate(T0_ints)

    allocate(Emo(nocc))
    do i = 1, nocc
        Emo(i) = mofckmat(i*(i+1)/2)
    end do

    ! save density, energy and field for subsequent calculations
    call openfile('pe_density.bin', luden, 'new', 'unformatted')
    rewind(luden)
    write(luden) Ene
    write(luden) qmnucs
    write(luden) Rm, Zm
    write(luden) npols
    write(luden) Fpd
    write(luden) nbas, nocc
    write(luden) full_denmat
    write(luden) cmo(:,1:nocc)
    write(luden) Emo
    close(luden)

    deallocate(full_denmat, Fpd, Emo)

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
    real(dp), dimension(:), allocatable :: core_fckmat, Fpd
    real(dp), dimension(:,:), allocatable :: frag_denmat, full_denmat
    real(dp), dimension(:,:), allocatable :: full_fckmat
    real(dp), dimension(:), allocatable :: overlap, repmat
    real(dp), dimension(:,:), allocatable :: full_overlap
    real(dp), dimension(:,:), allocatable :: full_rep
    real(dp), dimension(:), allocatable :: Emo
    real(dp), dimension(:,:), allocatable :: cmo

    work => dalwrk

    call openfile('pe_density.bin', luden, 'old', 'unformatted')
    rewind(luden)
    read(luden) Ene
    read(luden) pdnucs
    allocate(Rpd(3,pdnucs), Zpd(1,pdnucs))
    read(luden) Rpd, Zpd
    read(luden) npols
    allocate(Fpd(3*npols))
    read(luden) Fpd
    read(luden) fbas, focc
    allocate(frag_denmat(fbas, fbas))
    read(luden) frag_denmat
    allocate(cmo(fbas,focc))
    read(luden) cmo
    allocate(Emo(focc))
    read(luden) Emo
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
    call sirfck(full_fckmat, full_denmat, 1, isymdm, ifctyp, .true.,&
                work, size(work))
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

    do i = 1, focc
        cmo(:,i) = Emo(i) * cmo(:,i)
    end do
    write(lupri,*) Emo
    allocate(full_rep(cbas,cbas))
    full_rep = 0.0d0
    full_rep = - 1.0d0 * matmul(matmul(full_overlap(fbas+1:nbas,1:fbas), cmo),&
               matmul(transpose(cmo), full_overlap(1:fbas,fbas+1:nbas)))

    deallocate(full_overlap, cmo)

    allocate(repmat(cbas*(cbas+1)/2))
    l = 1
    do j = 1, cbas
        do i = 1, j
            repmat(l) = full_rep(i,j)
            l = l + 1
        end do
    end do

    deallocate(full_rep)

    ! save core Fock matrix
    call openfile('pe_fock.bin', lufck, 'new', 'unformatted')
    rewind(lufck)
    write(lufck) Fpd
    write(lufck) Ene
    write(lufck) core_fckmat
    write(lufck) repmat
    write(lufck) pdnucs
    write(lufck) Rpd, Zpd
    close(lufck)

    deallocate(core_fckmat, repmat, Rpd, Zpd, Fpd)

end subroutine pe_twoints

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
        elem2charge = 0.0d0
        return
    end if

    do i = 1, 112
        if (elem == trim(elements(i))) then
            elem2charge = real(i, dp)
            exit
        else
            elem2charge = 0.0d0
        end if
    end do

end function elem2charge

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

end module polarizable_embedding
