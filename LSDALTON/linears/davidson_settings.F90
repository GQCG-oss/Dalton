module davidson_settings
use kurtosis, only: PFMitem !KurtosisItem
use precision
use matrix_module!, only: matrix
use loc_types!, only: OrbitalLoc
use decompMod!, only: orbspread_data, decompitem
use queue_module!, only: modFIFO
use arhDensity!, only: solveritem
use ChargePrecMod !charge_precond, pmlowdin_lintra, pmmull_lintra, cllineartrans
TYPE RedSpaceItem

! *********************
! * General info      *
! *********************
! file
integer :: lupri

! *********************
! * General settings  *
! *********************
! Stepsize constraint
real(realk)  :: stepsize
! Maximum stepsize 
real(realk)  :: max_stepsize
! if using max elements
logical      :: use_max_element
! Maximum number of iterations
integer      :: max_it
! exit if approaching lin.dep.
logical      :: dep_exit
! use linesearch in localization
logical      :: linesearch
! T if using preconditioning
logical      :: precond
! True if only want locality info to be printed
logical      :: PRINT_INFO
! True if no localization on L2 is wanted
logical :: NOL2OPT
! True if we want to read in orbitals and localize with no HF calc
logical :: OnlyLocalize
! True if you want to print info useful for debug
logical :: debug_info
! Number of starting vectors
integer :: start_it
! TRUE if b3 coupling to b1 and b2 is too weak
logical :: singularity
! if using quadratic fit for line search
logical :: lines_fit
! *************************
! * convergence settings  *
! *************************
! Threshold for when residual is converged
real(realk)  :: conv_thresh
! thresh in local region
real(realk)  :: local_conv_thresh
! thresh in global region
real(realk)  :: global_conv_thresh
! Thresh for converging macro iterations
real(realk)  :: macro_thresh
! Maximum number of macro iterations 
integer :: max_macroit


!*******************
!* Used in solver  *
!*******************
!Ared : Reduced augmented hessian
real(realk), pointer :: Ared(:,:)
!Gred : reduced gradient ,vector
real(realk), pointer:: Gred(:)
!xred : vector of expansion coeff 
real(realk), pointer :: xred(:)
!all trial vectors and linear trans. are stored in mem
type(matrix), pointer :: Allb(:),AllSigma(:)
! Pointer to vector of diagonal elements regardless of function
type(matrix),pointer :: P
!previous levelshift.
real(realk) :: old_mu
! Preconditioner in realk
real(realk),pointer :: Pmat(:,:)


!**********************************
!* Orbital localization settings  *
!**********************************
! Info used when localizing orbitals using PSM
 type(orbspread_data) :: orbspread_inp
! Info used when localizing orbitals with Pipek/Lowdin
 type(PMitem)  :: PM_input
! Info used when localizing orbitals using PFM
 type(PFMitem) :: PFM_input
! true if using charge localiztion scheme
logical :: PM
! true if orbital variance localization scheme
logical :: orbspread
! true if using PFM localization scheme
logical :: PFM
! enable debug print outs for orbital localization
logical :: orb_debug
! print locality measures for all orbitals
logical :: all_orb_locality
! Specifies which orbital should be plotted
character*(10) :: plt_orbital
! set to true if orbital plt file should be made
logical   :: make_orb_plot
! quit after 10 iterations withoug significant changes in the localization measure
logical :: quit_after_10it
! which orbital indices reference the least/most local occ/virt
integer :: mostl_occ
integer :: mostl_virt
integer :: leastl_occ
integer :: leastl_virt
! offset that keeps track of where in CMO we are loclizing  
! offset = 1 for core, offset = ncore+1 for valence and offset = nocc+1 for virt
integer :: offset
! write CMOs to restart file each orbital_save_interval-iteration
integer :: orbital_save_interval
!**************************
!* ARH specific settings  *
!**************************
! contains D,G from prev. iterations
type(modFIFO),pointer :: fifoqueue
!> Contains solver info (ARH/TrFD)
type(solverItem),pointer :: arh
!> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
type(decompItem),pointer :: decomp
!> integer; 2 if antisymmetric
integer :: symm
!> if using arh linear transformations
logical :: arh_lintrans
!> if using arh precond
logical :: arh_precond
!> use davidson solver for arh
logical :: arh_davidson
!> Gradient norm, used for conver
real(realk) :: arh_gradnorm
!> linesearch on or off at given point in calc
logical :: arh_linesearch
!> linesearch on or off
logical :: arh_inp_linesearch
!> enable debug printouts
logical :: arh_davidson_debug
!> turned on if HOMO-LUMP gap is needed in reduced space Hessian
logical :: arh_extravec
logical :: arh_inp_extravec
!> if true, must perturb extravec
logical :: arh_dodft


!> debug 
logical :: arh_debug_linesearch
!> Energy computed in linesearch
real(realk) :: arh_linesE
!> difference in energies found via linesearch
real(realk) :: MaxLineSearchEnergyDiff
!> difference between actual energy from Fock matrix and energy found via linesearch
real(realk) :: ActualEnergyDiff
!> modification of thresholds for the line search energy compared to Fock Matrix
real(realk) :: LSmodthresh
!> determines if the ActualEnergyDiff have been set
logical :: EnergyDiffset

!*****************
!* Useful output *
!*****************
! Number of micro iterations
integer :: it
! levelshift
real(realk) :: mu
! Norm of solution vector X
real(realk) :: NormX
! lowest eigval of S
real(realk) :: eigval_S
! denominator of r : Q-En=GX+0.5xHx
real(realk) :: r_denom
! If step should be accepted
logical     :: step_accepted


ENDTYPE RedSpaceItem

CONTAINS


!> \brief contains default settings for use of davidson solver
!> both for orbital localization and arh.
subroutine davidson_default(CFG)
implicit none
type(RedSpaceItem) :: CFG

call davidson_reset(CFG)

nullify(CFG%decomp)
nullify(CFG%arh)
nullify(CFG%fifoqueue)
nullify(CFG%P)

! ARH related keywords
CFG%debug_info     = .false.
CFG%arh_davidson   = .false.
CFG%arh_precond    = .true.
CFG%arh_lintrans   = .false.
CFG%arh_linesearch = .false.
CFG%arh_extravec  = .false.
CFG%arh_inp_extravec  = .false.
CFG%arh_inp_linesearch   = .false.
CFG%arh_debug_linesearch = .false.
CFG%arh_davidson_debug   = .false.
CFG%arh_dodft            = .false.

! Orbital loc. related keywords
CFG%precond     =.true.
CFG%PRINT_INFO  = .false.
CFG%PM          = .false.
CFG%orbspread   = .false.
CFG%PFM         = .false.
CFG%NOL2OPT     = .false.
CFG%OnlyLocalize       = .false.
CFG%PFM_input%TESTCASE = .false.
CFG%all_orb_locality   = .false.
CFG%orb_debug          = .false.
CFG%make_orb_plot      = .false.
CFG%linesearch         =.true.
CFG%PM_input%PipekMezeyLowdin  = .false.
CFG%PM_input%linesearch        = .false.
CFG%PM_input%precond           = .true.
CFG%PM_input%orb_debug         = .false.
CFG%leastl_occ=0
CFG%leastl_virt=0
CFG%mostl_occ=0
CFG%mostl_virt=0
! Convergence thresh for micro
CFG%conv_thresh = 0.01_realk
! Global convergence threshold for micro
CFG%global_conv_thresh = 0.01_realk
! Local convergence threshold for micro
CFG%local_conv_thresh = 0.005_realk
CFG%lines_fit =.true.
CFG%quit_after_10it = .true.
CFG%orbital_save_interval = 20



!General solver related keywords
CFG%max_it = 25
CFG%max_macroit = 200
end subroutine davidson_default

subroutine davidson_default_SCF(CFG)
implicit none
type(RedSpaceItem) :: CFG



! line search spec. keywords
CFG%EnergyDiffset=.FALSE.
CFG%arh_linesE = 0.0E0_realk
CFG%ActualEnergyDiff = 0.0E0_realk
CFG%MaxLineSearchEnergyDiff = 0.0E0_realk
CFG%LSmodthresh = 1000.0E0_realk !correspond to 10^-5 on the energy on 1. iteration

end subroutine davidson_default_SCF


subroutine davidson_reset(CFG)
implicit none
type(RedSpaceItem) :: CFG



CFG%singularity=.false.
CFG%start_it = 0
!**********************************
!*       Solver settings          *
!**********************************
! convergence threshold for macro
CFG%macro_thresh = 0.0001_realk

!*********************************
!* Stepsize specific settings    *
!*********************************
! stepsize
CFG%stepsize =0.5_realk
! maximum stepsize
CFG%max_stepsize =0.5_realk
! Initialize mu
CFG%mu = 0.0_realk
! Initialize denominator
CFG%r_denom = 1.0_realk
! Initialize micro it
CFG%it = 0

end subroutine davidson_reset




end module davidson_settings










