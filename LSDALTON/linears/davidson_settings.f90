module davidson_settings
use kurtosis!, only: KurtosisItem
use precision
use matrix_module!, only: matrix
use loc_types!, only: OrbitalLoc
use decompMod!, only: orbspread_data, decompitem
use queue_module!, only: modFIFO
use arhDensity!, only: solveritem

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
 type(orbspread_data),pointer :: orbspread_input
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
!> linesearch on or off
logical :: arh_linesearch

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
! true if we want some extra trial vecs for ARH 
logical :: arh_extravecs

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

subroutine davidson_default_SCF(CFG)
implicit none
type(RedSpaceItem) :: CFG

call davidson_reset(CFG)

nullify(CFG%decomp)
nullify(CFG%arh)
nullify(CFG%fifoqueue)
nullify(CFG%P)

CFG%EnergyDiffset=.FALSE.
CFG%arh_linesE = 0.0E0_realk
CFG%ActualEnergyDiff = 0.0E0_realk
CFG%MaxLineSearchEnergyDiff = 0.0E0_realk
CFG%LSmodthresh = 1000.0E0_realk !correspond to 10^-5 on the energy on 1. iteration
CFG%max_it = 20

CFG%debug_info     = .false.
CFG%arh_davidson   = .false.
CFG%arh_precond    = .true.
CFG%arh_lintrans   = .false.
CFG%arh_linesearch = .false.
CFG%arh_extravecs  = .false.
CFG%arh_debug_linesearch = .false.
CFG%precond=.true.
CFG%PRINT_INFO = .false.
CFG%PM = .false.
CFG%orbspread = .false.
CFG%PFM  = .false.
CFG%NOL2OPT   = .false.
CFG%OnlyLocalize   = .false.
CFG%PFM_input%TESTCASE = .false.
CFG%PM_input%ChargeLocMulliken = .false.
CFG%PM_input%ChargeLocLowdin   = .false.
CFG%PM_input%PipekMezeyLowdin  = .false.
CFG%PM_input%linesearch        = .false.

end subroutine davidson_default_SCF

subroutine davidson_default_OrbLoc(CFG)
implicit none
type(RedSpaceItem) :: CFG

call davidson_reset(CFG)

CFG%precond=.true.
CFG%PRINT_INFO = .false.
CFG%PM = .false.
CFG%orbspread = .false.
CFG%PFM  = .false.
CFG%NOL2OPT   = .false.
CFG%OnlyLocalize   = .false.
CFG%PFM_input%TESTCASE = .false.
CFG%PM_input%ChargeLocMulliken = .false.
CFG%PM_input%ChargeLocLowdin   = .false.
CFG%PM_input%PipekMezeyLowdin  = .false.
CFG%PM_input%linesearch        = .false.
CFG%PM_input%precond           = .true.
CFG%arh_davidson   = .false.
CFG%arh_precond    = .true.
CFG%arh_lintrans   = .false.
CFG%arh_linesearch = .false.
CFG%arh_extravecs  = .false.
CFG%arh_debug_linesearch = .false.
CFG%debug_info = .false.
CFG%max_it = 25

end subroutine davidson_default_OrbLoc

subroutine davidson_reset(CFG)
implicit none
type(RedSpaceItem) :: CFG

!**********************************
!*       Solver settings          *
!**********************************
! Convergence thresh for micro
CFG%conv_thresh = 0.05_realk
! Global convergence threshold for micro
CFG%global_conv_thresh = 0.05_realk
! Local convergence threshold for micro
CFG%local_conv_thresh = 0.005_realk
! convergence threshold for macro
CFG%macro_thresh = 0.0001_realk

!*********************************
!* Stepsize specific settings    *
!*********************************
! stepsize
CFG%stepsize =0.75_realk
! maximum stepsize
CFG%max_stepsize =0.75_realk
! Initialize mu
CFG%mu = 0d0
end subroutine davidson_reset



end module davidson_settings










