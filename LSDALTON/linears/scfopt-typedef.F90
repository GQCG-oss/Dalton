!> \author S. Host
!> \date March 2010
!> \brief Contains structure with info about SCF optimization type.
!> 
!> **********************************
!> ** SCF optimization parameters. **
!> **********************************
!> 
module opttype
use precision

type OptItem
   !CONFIGURATION:
   !==============    
      !> Logical unit number for LSDALTON.OUT
      integer :: lupri
      !> Logical unit number for LSDALTON.ERR
      integer :: luerr
      !> Nuclear repulsion
      real(realk)          :: potnuc
      !> Which type of starting guess?
      character(len=16) :: cfg_start_guess
      !> Should the H1DIAG guess be asymmetric? (for unrestricted)
      logical :: cfg_asym
      !> True is unrestricted calculation
      logical :: cfg_unres
      !> Defines what type of density optimization should be used
      integer     :: cfg_density_method
      !
      !** Which F_(i) => D_(i+1) scheme should be used? **
      !
      !The potential-to-density phase configuration. A number of options
      !are available.
      !> Standard diagonalization - the Roothaan way
      integer :: cfg_f2d_roothaan     ! = 1  
      !> Direct density optimization. Minimize E = 2TrFD(X) 
      integer :: cfg_F2D_direct_dens  ! = 2 
      !> Purification scheme - NO LONGER SUPPORTED! /Stinne 16-08-2010
      !integer :: cfg_F2D_purification ! = 3
      !> Augmented Roothaan-Hall optimization
      integer :: cfg_F2D_arh          ! = 4               
      !> Augmented Roothaan-Hall optimization
      logical :: cfg_oao_gradnrm      ! = .TRUE.               
      !> Which purification method? (1=TCP, 2=TRS, 3=MCW, 4=PM)
      !integer :: cfg_purification_method - NO LONGER SUPPORTED! /Stinne 16-08-2010
      !> Convergence threshold on the gradient norm
      real(realk) :: cfg_convergence_threshold
      !> Multiply convergence threshold for level 2 with this factor (when doing trilevel)
      real(realk) :: cfg_level2_convfactor
      !> Maximum number of SCF iterations
      integer     :: cfg_max_linscf_iterations
      !> For large calculations, it is advantageous to use dynamic threshold instead of gradient norm
      logical     :: cfg_convdyn
      !> If dynamic threshold is used, should this be sloppy, standard, or tight? 
      character*5 :: cfg_convdyn_type
      !> Some parameters depend on whether we do HF or DFT, so keep track of this here:
      integer     :: calctype
      !> Calctype can be set to HF or dft:
      integer     :: hfcalc
      !> Calctype can be set to HF or dft:
      integer     :: dftcalc
      !> Use Block-Sparse Matrices:
      logical     :: cfg_prefer_CSR
      !> Use SCALAPACK Matrices:
      logical     :: cfg_prefer_SCALAPACK
      !> Use PDMM Matrices:
      logical     :: cfg_prefer_PDMM
      !> Should we crash the calculation - for debugging purposes
      logical     :: crashcalc
      !> Should incremental scheme be used for integrals?
      logical     :: cfg_incremental
      !> Should we save old F0 and D0 related to incremental and line search
      logical     :: cfg_saveF0andD0
      !> Dump previous density/Fock matrices to disk when constructing new Fock matrix (saves memory)
      logical     :: cfg_queue_on_disk
      !> If true, calculate only energy (i.e., no optimization) and hes eigenval for constructed by weighted
      !> average of densities for two minima. Use with files D1 and D2 (densities). 
      logical     :: cfg_hesonly
      !> Same as above, but find hessian eigenval by setting up explicit hessian
      !> and diagonalozing, instead of using iterative scheme. Expensive, of course!
      logical     :: cfg_diaghesonly
      !> Read from input. Weight between densities for the two minima
      real(realk) :: cfg_weight_param
      !> True if trust-region scheme is used for SCF optimization
      logical     :: do_trustregion      
      !> Do Second Order Ensemble Optimization
      logical     :: cfg_soeo
      !> Dump F, D matrices to disk in every SCF iteration (for investigating sparsity)
      logical     :: dumpmatrices
      !> Indicates which level of optimization is currently being carried out
      integer     :: optlevel  !1=atoms, 2=valence, 3=full 
   !SETTINGS:
   !=========
      !> Current setting for convergence threshold on the gradient norm (might be changed during opt.)
      real(realk) :: set_convergence_threshold
      !> !Should scaling of virtual orbitals be used?
      logical     :: cfg_scale_virt
   !INFO OPTIONS:
   !=============
      logical     :: info_matop
      !> Should the final MO be printed to output and file? (cmo.out)
      logical     :: print_final_cmo
   !DEBUG OPTIONS:
   !==============
      !> Convert files between formatted and unformatted form to avoid little/big-endian problems
      logical     :: debug_convert
      !> How should file be converted? (see debug.f90)
      integer     :: cfg_which_conversion
      !> Diagonalize full second order Hessian after SCF optimization
      logical     :: debug_diag_hessian
      !whteher calling lsquit if convergence is not reached
      logical     :: opt_quit
      !if Atoms start density should be disregarded for trilevel
      logical     :: add_atoms_start
      !Perform McWeeny purification on the non idempotent Atoms Density
      logical     :: MWPURIFYATOMSTART
      !Dense Matrix Type in Level 2 of Trilevel (requires .START = TRILEVEL)
      logical     :: DENSELEVEL2
end type OptItem

contains

!> \brief Set default configuration for SCF optimization.
!> \author S. Host
!> \date March 2010
subroutine opt_set_default_config(opt)
implicit none
   !> Used to store info about SCF optimization
   type(OptItem), intent(inout) :: opt

   !CONFIGURATION:
   !==============
   opt%cfg_f2d_roothaan     = 1  
   opt%cfg_F2D_direct_dens  = 2 
   !opt%cfg_F2D_purification = 3
   opt%cfg_F2D_arh          = 4

   opt%cfg_oao_gradnrm           =.true.
   opt%cfg_start_guess           ='ATOMS'
   opt%cfg_asym                  =.false.
   opt%cfg_unres                 =.false.

   opt%cfg_density_method        = opt%cfg_F2D_arh  !Default is Augmented Roothaan-Hall
   opt%cfg_convergence_threshold = 1e-4
   opt%cfg_level2_convfactor     = 10E0_realk !Level 2 usually doesn't need to be converged as hard as level 3
                                         !Default level 2 threshold is 100*level 3 threshold (change with .L2THR)
   opt%cfg_max_linscf_iterations = 100
   opt%cfg_convdyn               = .false.
   !opt%cfg_purification_method   = 1 !purification method (1=TCP, 2=TRS, 3=MCW, 4=PM)
   opt%hfcalc                    = 1
   opt%dftcalc                   = 2
   opt%calctype                  = opt%hfcalc !Default is HF
   opt%cfg_prefer_CSR            = .false.
   opt%cfg_prefer_SCALAPACK      = .false.
   opt%cfg_prefer_PDMM           = .false.
   opt%crashcalc                 = .false.
   opt%cfg_incremental           = .false.
   opt%cfg_saveF0andD0           = .false.
   opt%cfg_queue_on_disk         = .false.
   opt%cfg_hesonly               = .false.
   opt%cfg_diaghesonly           = .false.
   opt%cfg_weight_param          = 0.5E0_realk
   opt%do_trustregion            = .false.
   opt%cfg_soeo                  = .false.
   opt%dumpmatrices              = .false.
   opt%optlevel                  = 3
   !SETTINGS:
   !=========
   opt%set_convergence_threshold = 1e-4
   opt%cfg_scale_virt            = .false.
   !INFO OPTIONS:
   !=============
   opt%info_matop                = .FALSE.
   opt%print_final_cmo           = .false.
   !DEBUG OPTIONS:
   !==============
   opt%debug_convert             = .false.
   opt%debug_diag_hessian        = .false.
   opt%opt_quit                  = .true.
   opt%add_atoms_start           = .true.
   opt%MWPURIFYATOMSTART         = .false.
   opt%DENSELEVEL2               = .false.

end subroutine opt_set_default_config

end module opttype
