module pe_variables

    use pe_precision

#if defined(VAR_MPI)
#if defined(VAR_USE_MPIF)
    implicit none
#include "mpif.h"
#else
    use mpi
    implicit none
#endif
#endif

    ! options
    logical, save :: peqm = .false.
    logical, save :: pe_pot = .false.
    logical, save :: pe_iter = .true.
    logical, save :: pe_diis = .false.
    logical, save :: pe_mixed = .false.
    logical, save :: pe_redthr = .false.
    logical, save :: pe_border = .false.
    logical, save :: pe_damp = .false.
    logical, save :: pe_gspol = .false.
    logical, save :: pe_nomb = .false.
    logical, save :: pe_gauss = .false.
    logical, save :: pe_polar = .false.
    logical, save :: pe_mep = .false.
    logical, save :: pe_skipqm = .false.
    logical, save :: pe_twoint = .false.
    logical, save :: pe_repuls = .false.
    logical, save :: pe_savden = .false.
    logical, save :: pe_fd = .false.
    logical, save :: pe_fdes = .true.
    logical, save :: lvdw  = .false.
    logical, save :: pe_sol = .false.
    logical, save :: pe_noneq = .true.
    logical, save :: pe_infld = .false.
    logical, save :: pe_restart = .false.
    logical, save :: pe_verbose = .false.
    logical, save :: pe_debug = .false.
    logical, save :: rsp_first = .true.
    logical, save :: fixpva = .false.

    ! calculation type
    logical, save :: fock = .false.
    logical, save :: energy = .false.
    logical, save :: response = .false.
    logical, save :: molgrad = .false.
    logical, save :: mep = .false.
    logical, save :: noneq = .false.
    logical, save :: london = .false.

    ! temporary solution for work array thing
    real(dp), dimension(:), pointer :: work

    ! filenames
    character(len=80) :: potfile = 'POTENTIAL.INP'
    character(len=80) :: surfile = 'SURFACE.INP'

    ! MPI stuff
#if defined(VAR_MPI)
    integer, parameter :: comm = MPI_COMM_WORLD
    integer, parameter :: impi = MPI_INTEGER
    integer, parameter :: rmpi = MPI_REAL8
    integer, parameter :: lmpi = MPI_LOGICAL
#endif
    integer, save :: myid, nprocs, ierr
    integer, save :: site_start, site_finish
    integer, save :: surp_start, surp_finish
    integer, save :: mep_start, mep_finish
    logical, save :: synced = .false.
    integer, dimension(:), save, allocatable :: siteloops, surploops, meploops
    integer, dimension(:), save, allocatable :: poldists, sitedists, surpdists
    integer, dimension(:), save, allocatable :: mepdists
    integer, dimension(:), save, allocatable :: displs

    ! logical unit for output file (default is stdout)
    integer, save :: luout = 6

    ! constants, thresholds and stuff
    ! 1 bohr = 0.5291772109217 Aa (codata 2010)
    real(dp), parameter :: aa2au = 1.0 / 0.5291772109217
    real(dp), parameter :: aa2au2 = aa2au * aa2au
    real(dp), parameter :: pi = 3.141592653589793
    real(dp), parameter :: zero = 1.0d-8
    integer, save :: scfcycle = 0
    real(dp), save :: thriter = 1.0d-5
    real(dp), save :: damp = 2.1304
    real(dp), save :: gauss_factor = 1.0
    real(dp), dimension(3), save :: repfacs = 1.0
    real(dp), save :: Rmin = 2.2
    character(len=6), save :: border_type = 'REDIST'
    ! use Cholesky factorization of classical response matrix
    logical, save :: chol = .true.
    ! solvent and dielectric constant (defaults to water)
    character(len=80) :: solvent
    real(dp), save :: eps = 0.0
    real(dp), save :: epsinf = 0.0


    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476
    ! C^(n)_ij coefficients for calculating T(k) tensor elements
    integer, dimension(:,:,:), allocatable, save :: Cnij

    ! polarizable embedding potential info
    ! ------------------------------------
    ! total number of classical sites
    integer, save :: nsites = 0
    ! number of polarizable sites
    integer, save :: npols = 0
    ! number of surface points
    integer, save :: nsurp = 0
    ! number fragment densities
    integer, save :: nfds = 0
    ! number of nuclei in fragment density
    integer, save :: fdnucs = 0
    ! exclusion list length
    integer, save :: lexlst = 0
    ! number of density matrices
    integer :: ndens = 0
    ! number of basis functions in core fragment
    integer, save :: nbas
    ! size of packed matrices
    integer, save :: nnbas
    ! number of nuclei in core region
    integer, save :: qmnucs = 0

    ! specifies what type of parameters are present
    ! lmul(0): monopoles, lmul(1): dipoles etc.
    logical, dimension(0:5), save :: lmul = .false.
    ! lpol(1): (an)isotropic dipole-dipole polarizabilities
    logical, dimension(1), save :: lpol = .false.
    ! lhypol(1): dipole-dipole-dipole polarizabilities/1st hyperpolarizability
!    logical, dimension(1), save :: lhypol
!    ! lvdw: LJ parameters
!    logical, dimension(1), save :: lvdw = .false.

    ! charges, areas, coordinates, elements and exclusion lists
    ! site nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zs
    ! fragment density nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zfd
    ! core fragment nuclear charges
    real(dp), dimension(:,:), allocatable, save :: Zm
    ! surface point areas
    real(dp), dimension(:), allocatable, save :: Sa
    ! site coordinates
    real(dp), dimension(:,:), allocatable, save :: Rs
    ! fragment density nuclear coordinates
    real(dp), dimension(:,:), allocatable, save :: Rfd
    ! core fragment nuclear coordinates
    real(dp), dimension(:,:), allocatable, save :: Rm
    ! surface point coordinates
    real(dp), dimension(:,:), allocatable, save :: Sp
    ! site elements
    character(len=2), dimension(:), allocatable, save :: elems
    ! exclusion list
    integer, dimension(:,:), allocatable, save :: exclists

    ! energy contributions
    ! total
    real(dp), dimension(:), allocatable, save :: Epe
    ! electrostatic
    real(dp), dimension(:,:), allocatable, save :: Ees
    ! polarization
    real(dp), dimension(:,:), allocatable, save :: Epol
    ! continuum solvation
    real(dp), dimension(:,:), allocatable, save :: Esol
    ! fragment density
    real(dp), dimension(:,:), allocatable, save :: Efd
    ! LJ 
    real(dp), dimension(:), allocatable, save :: Elj

    ! multipole moments
    ! order of the highest order multipole moment
    integer, save :: mulorder = -1
    ! monopoles, dipoles, quadrupoles, octopoles, etc.
    real(dp), dimension(:,:), allocatable, save :: M0s, M1s, M2s, M3s, M4s, M5s
    ! (hyper)polarizabilities
    ! order of highest order polarizability
    integer, save :: polorder = -1
    ! dipole-dipole polarizabilities
    real(dp), dimension(:,:), allocatable, save :: P1s
    ! .true. if P1 > 0 else .false.
    logical, dimension(:), allocatable, save :: zeroalphas

    ! LJ parameters for MM region - from pot file
    real(dp), dimension(:,:), allocatable, save :: LJs
    ! LJ paramters for QM region - from dal file
    real(dp), dimension(:,:), allocatable, save :: qmLJs

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
    integer, save :: npoints
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
    real(dp), save :: xsize = 5.0
    real(dp), save :: ysize = 5.0
    real(dp), save :: zsize = 5.0

    ! Internal field stuff
    ! --------------------
    ! Coordinates on which potential and field are calculated
    real(dp), dimension(:,:), allocatable, save :: crds
    ! Number of coordinates (length of crds)/3
    integer, save :: ncrds

    ! FIXPVA2 stuff       
    ! Maximum number of tessera
    integer, save :: MXFFTS = 1
    ! Number of tessera per atom
    integer, save :: NTSATM = 60 
    ! Number of tessera
    integer, save :: NFFTS = 1
    !
    real(dp), save :: TOANGS = 0.52917724924D+00
    !
    real(dp), save :: fixtol = 1.0d-10
    ! --------------------

end module pe_variables
