!> @file
!> Contains ARH optimization routines.

!> \brief Contains routines used by the ARH/Direct Dens solver.
!> \author S. Host
!> \date 2007
!>
!> ARH: S. Høst, B. Jansík, J. Olsen et al. PCCP 10, 5344 (2008)   \n
!>      S. Høst, J. Olsen, B. Jansík et al. JCP 129, 124106 (2008) \n
!>                                                                 \n
!> Direct density: P. Sałek, S. Høst, L. Thøgersen et al. JCP 126, 114110 (2007)
!>
module arhDensity
use precision
use matrix_module
use decompMod
use queue_ops
use queue_module
use lsdalton_fock_module
use memory_handling
use dal_interface
use matrix_operations_aux
use matrix_operations
use files
use matrix_util
use typedeftype
use II_XC_interfaceModule
private
public :: arh_symmetric, arh_antisymmetric, SolverItem, arh_set_default_config,&
     & arh_test_convergence, fifo_inverse_metric,&
     & arh_get_weights, arh_get_M, fifo_inv_metric_times_vector,&
     & arh_precond, arh_PCG, arh_lintrans, arh_xdep_matrices,&
     & arh_crop_x_and_res,  arh_crop_intermed_sub,&
     & arh_crop_setup_redsp, arh_crop_optimal_x, arh_crop_extra_dim,&
     & arh_get_TR_denom!, epred, debug_get_hessian
!> Used to pass info about symmetry of trial vectors/matrices
integer, parameter :: arh_symmetric = 1
!> Used to pass info about symmetry of trial vectors/matrices
integer, parameter :: arh_antisymmetric = 2

!> \author S. Host
!> \date March 2010
!> \brief Contains settings for ARH/Direct Dens solver
!>
!> The idea with 'configurations' and 'settings' is that the configuration is
!> never changed, after it has been set (from default or input). In this way,
!> the 'setting', which can be changed during the optimization, may always be
!> reset to whatever it originally was.
!>
type SolverItem
   !CONFIGURATIONS:
      !> Logical unit number for LSDALTON.OUT
      integer      :: lupri
      !> True if DFT calculation
      logical      :: do_dft
      !> MO coefficients to be saved if local Link
      type(matrix) :: C
!      !> Integral settings (for 2nd order opt)
!      type(lssetting), pointer :: integr_settings 

      !Solver configuration:
      !---------------------
      !> If true, do ARH. If false, do direct density optimization.
      logical     :: cfg_arhterms
      !> If true, use Conjugate Residual Optimal Vectors scheme
      logical     :: cfg_arh_crop
      !> If true, do an extra preconditioning in each micro it. Preliminary evidence suggests this is necessary
      logical     :: cfg_arh_crop_safe
      !> If true, use Conjugate Residual Optimal Vectors scheme and save only few vectors
      logical     :: cfg_arh_truncate     
      !> Number of vectors to keep when using cfg_arh_truncate
      integer     :: cfg_arh_microvecs
      !> If true, do not use diagonal preconditioning when solving equations.
      logical     :: cfg_noprec
      !> Do 2nd order optimization in the whole calculation
      logical     :: cfg_2nd_order_all
      !> Switch to 2nd order optimization in the local region
      logical     :: cfg_2nd_order_local
      !> Second order optimization, i.e. construct Fock matrix in each micro it. Expensive!!!
      logical     :: cfg_do_2nd_order
      !> True if trial and sigma vectors should be kept on disk instead of in core
      logical     :: cfg_arh_disk_micro
      !> True if density and Fock matrices should be kept on disk instead of in core
      !logical     :: cfg_arh_disk_macro !Not active - get_from_modFIFO_disk won't work!
      !> True if the solver is used to localize orbitals
      logical     :: cfg_orbspread
      !> Check if residual is going down for every nits_check iterations
      integer     :: nits_check
      !> We require that residual is reduced by (1-error_decrease)% for every nits_check iteration
      real(realk) :: error_decrease
      !> If error has not decreased enough, contract trust radius by this factor
      real(realk) :: contract_factor

      type(orbspread_data), pointer :: orbspread_input

      !Trust-region update configuration:
      !----------------------------------
      !> Expand trust radius if ratio is larger than cfg_arh_contract_crit
      REAL(REALK) :: cfg_arh_contract_crit
      !> Contract trust radius if ratio is smaller than cfg_arh_expand_crit
      REAL(REALK) :: cfg_arh_expand_crit
      !> On expansion, trust radius is expanded by the factor cfg_arh_expand 
      REAL(REALK) :: cfg_arh_expand
      !> On contraction, trust radius is contracted by the factor cfg_arh_contract 
      REAL(REALK) :: cfg_arh_contract

      !Levelshift parameters:
      !----------------------
      !> Level shift scheme which should be better suited for CROP. Not fully tested.
      logical     :: cfg_arh_newdamp 
      !> Do no level shifting
      logical     :: cfg_nodamp 
      !> A minimum level shift can be set
      real(realk) :: cfg_min_lshift
      !> Run with fixed level shift
      logical     :: cfg_fixed_shift
      !> Value of fixed level shift
      real(realk) :: cfg_fixed_shift_param
      !> Absolute max allowed Frobenius norm of X, regardless of trust radius
      real(realk) :: cfg_max_step
      !> Absolute max allowed element of X, regardless of trust radius
      real(realk) :: cfg_max_element
      !> Level shift by calculating HOMO-LUMO gap in each SCF iteration
      logical     :: lshift_by_hlgap

      !Thresholds:
      !-----------
      !> !Convergence threshold for micro iterations (PCG or reduced space)
      real(realk) :: cfg_micro_thresh
   !SETTINGS:
   !=========
      !> If true, do ARH. If false, do direct density optimization.
      logical     :: set_arhterms
      !> Second order optimization, i.e. construct Fock matrix in each micro it. Expensive!!!
      logical     :: set_do_2nd_order
      !> Current max allowed Frobenius norm of X (i.e. trust radius)
      real(realk) :: set_max_step
      !> Current max allowed element of X (i.e. trust radius)
      real(realk) :: set_max_element
      !> Set true when we decide that we are local (if xnorm > 0.2 and SCF it > 4)
      logical     :: set_local
      !> When true, optimize max element instead of Frobenius norm (initial SCF iterations). This is size extensive.  
      logical     :: set_optxelm

   !INFO VARIABLES:
   !===============
      !> Print detailed info from CROP solver
      logical     :: info_crop
      !> Print detailed info about ARH weights
      logical     :: info_arh_weights
      !> Print info about convergence of linear equations
      logical     :: info_lineq
      !> Print info about level shift
      logical     :: info_levelshift

   !DEBUG VARIABLES:
   !================
      !> Print detailed info from ARH linear transformation
      logical     :: debug_arh_lintra
      !> Do PCG instead of diagonal preconditioning
      logical     :: debug_arh_precond
      !> Current SCF iteration. For printing debug info
      integer     :: debug_arh_scfit
      !> Detailed debug info from ARH/DD calculation (expensive!)
      logical     :: debug_dd 
      !> Diagonalize FUP and FUQ in each SCF iteration to get all orbital energies
      logical     :: debug_dd_homolumo
      !> Print the RH percentage of total ARH linear transformation
      logical     :: debug_dd_lintra
      !> Currently lowest eigenvalue in reduced space in each iteration
      logical     :: debug_diag_redspace
      !> Set up and diagonalize ARH Hessian in each SCF iteration (expensive!)
      logical     :: debug_hessian
      !> As above, but full 2nd order Hessian instead of ARH Hessian (very, very expensive!). 
      logical     :: debug_hessian_exact
      !> Very detailed debug info can be printed after this number of SCF its. Currently disabled.
      integer     :: cfg_nits_debug
   !DATA:
   !=====
      !> Energy from previous SCF iteration
      real(realk)          :: old_energy 
      !> Contains reduced ARH Hessian
      real(realk), pointer :: Ared(:,:)
      !> Contains reduced gradient
      real(realk), pointer :: Gred(:)
      !> Contains overlap of trial vectors (for debugging. Diagonalize to test for linear dependencies).
      real(realk), pointer :: Sred(:,:)
      !> Intermediate matrix used with CROP scheme
      real(realk), pointer :: CROPmat(:,:)
      !> ARH metric: Metr_ij = Tr(Di - D0)(Dj - D0)
      real(realk), pointer :: fifometric(:,:)
      !> Inverted ARH metric
      real(realk), pointer :: inv_fifometric(:,:)
      !> ndens*ndens matrix needed for ARH: M_ij = 1/2 * Tr((D_i - D0)*(F(D_j)-F(D0)) + (D_j - D0)*(F(D_i)-F(D0)) )
      real(realk), pointer :: fifoM(:,:)
      !> Number of rejection steps in current SCF iteration
      integer                  :: Nrejections
      !> Part of new density that can be expanded in subspace of densities
      real(realk)              :: D_para_D_tot
      !> Final dimension of reduced space after solution of micro equations
      real(realk)              :: final_redspacedim
      !> Denominator of trust radius
      real(realk)              :: denom
      !> True if step is accepted, i.e. if SCF energy is decresed
      logical                  :: step_accepted
      !> Max element of final X
      real(realk)              :: maxelm
      !> Frobenius norm of final X
      real(realk)              :: xnorm
      !> Final mu in solution of linear equations (save here for later printing)
      real(realk)              :: current_mu
      !> True if trust radius had to decreased during solution of linear equations to obtain convergence 
      logical                  :: trustradius_decreased
      !> Gradient norm in OAO basis
      real(realk)              :: OAO_gradnrm
      !> Whether OAO_gradnrm is used or not
      logical                  :: OAO_gradnrm_exist

      !> Starting guess available
      logical                  :: starting_guess_defined
      !> Starting guess matrix pointer
      type(Matrix), pointer    :: starting_guess

      integer                  :: newton_it
      !> Matrix containg diagonal Hessian elements for preconditioning
      type(matrix),pointer :: P

      !DEBUG DATA:
      !============
      !> SCF iteration
      integer     :: scfit
      !> Final lowest eigenvalue in reduced space
      real(realk) :: final_redspace_eival
      !> HOMO-LUMO gap obtained by diagonalization
      real(realk) :: diag_hlgap
      !> HOMO-LUMO gap obtained iteratively
      real(realk) :: iter_hlgap
      !> Lowest full 2nd order Hessian eigenvalue
      real(realk) :: heseival
      !> Lowest ARH Hessian eigenvalue
      real(realk) :: arheival
end type SolverItem

contains

   !> \brief Set default configuration for Augmented Roothaan-Hall.
   !> \author S. Host
   !> \date March 2010
   subroutine arh_set_default_config(arh)
   implicit none
      !> Contains solver info (ARH/TrFD)
      type(SolverItem) :: arh

      !DFT calculation?
      arh%do_dft               = .false.
      !nullify(arh%integr_settings)
      !Solver configuration:
      !---------------------
      arh%cfg_arhterms         = .true.
      arh%cfg_arh_crop         = .true.
      arh%cfg_arh_crop_safe    = .true.
      arh%cfg_arh_truncate     = .true.
      arh%cfg_arh_microvecs    = 2
      arh%cfg_noprec           = .false.
      arh%cfg_2nd_order_all    = .false.
      arh%cfg_2nd_order_local  = .false.
      arh%cfg_do_2nd_order     = .false.
      arh%cfg_arh_disk_micro   = .false.
      !arh%cfg_arh_disk_macro   = .false.
      arh%cfg_orbspread        = .false.

      !Trust-region update configuration:
      !----------------------------------
      arh%cfg_arh_expand        = 1.2E0_realk
      arh%cfg_arh_contract      = 0.7E0_realk
      arh%cfg_arh_contract_crit = 0.25E0_realk
      arh%cfg_arh_expand_crit   = 0.75E0_realk
      !Levelshift parameters:
      !----------------------
      arh%cfg_arh_newdamp       = .false.
      arh%cfg_nodamp            = .false.
      arh%cfg_min_lshift        = 0.0E0_realk
      arh%cfg_fixed_shift       = .false.
      arh%cfg_fixed_shift_param = 0.0E0_realk
      arh%cfg_max_step          = 0.6E0_realk
      arh%cfg_max_element       = 0.35E0_realk
      arh%lshift_by_hlgap       = .true.
      arh%nits_check            = 15
      arh%error_decrease        = 0.85E0_realk
      arh%contract_factor       = 0.9E0_realk
      !Thresholds:
      !-----------
      arh%cfg_micro_thresh  = 1.0E-2_realk
   !INFO VARIABLES:
   !=========================
      arh%info_crop           = .false.
      arh%info_arh_weights    = .false.
      arh%info_lineq          = .true.
      arh%info_levelshift     = .false.
   !DEBUG VARIABLES:
   !================
      arh%debug_arh_lintra    = .false.
      arh%debug_arh_precond   = .false.
      arh%debug_dd            = .false.
      arh%debug_dd_homolumo   = .false.
      arh%debug_dd_lintra     = .false.
      arh%debug_diag_redspace = .false.
      arh%debug_hessian       = .false.
      arh%debug_hessian_exact = .false.
   !SETTINGS:
   !=========
      arh%set_arhterms     = .true.
      arh%set_do_2nd_order = .false.
      arh%set_max_step    = 0.6E0_realk
      arh%set_max_element = 0.35E0_realk
      arh%set_local       = .false.
      arh%set_optxelm     = .false.
      arh%current_mu      = 0.0E0_realk
   !DATA:
   !=====
      arh%Nrejections     = 0
      arh%step_accepted   = .true.
      arh%starting_guess_defined = .false.
      nullify(arh%Ared)
      nullify(arh%Gred)
      nullify(arh%Sred)
      nullify(arh%CROPmat)
      nullify(arh%fifometric)
      nullify(arh%inv_fifometric)
      nullify(arh%fifoM)
      arh%OAO_gradnrm = 1E10_realk
      arh%OAO_gradnrm_exist= .false.

   end subroutine arh_set_default_config

   !> \brief Test if linear equations have converged.
   !> \author S. Host
   !> \date March 2010
   subroutine arh_test_convergence(arh,err,errsave,x,mu,mumax,thresh,it,j,k,&
   &                               xsave_lu,redspacedim_save,maxvec,do_LS,done)
   implicit none
        !> Contains solver info (ARH/TrFD)
        type(solverItem),intent(inout) :: arh
        !> Frobenius norm of current residual
        real(realk), intent(in)        :: err
        !> Save residual with certain interval so we can check if residual is decreasing
        real(realk), intent(inout)     :: errsave
        !> Current trial vector
        type(matrix), intent(inout)    :: x
        !> Current level shift
        real(realk), intent(inout)     :: mu
        !> Max level shift. Obtained when linear equations have converged using Frobenius norm as trust radius
        real(realk), intent(inout)     :: mumax
        !> Threshold for convergence of linear equations
        real(realk), intent(in)        :: thresh
        !> Number of current micro iteration
        integer, intent(inout)         :: it
        !> j = 0 when using Frobenius norm as trust radius, j = 1 when using max element as trust radius
        integer, intent(inout)         :: j
        !> Keeps track of how often we should check if residual is decreasing
        integer, intent(inout)         :: k
        !> Logical unit number for saving current trial vector to file, when equations have converged using Frobenius norm as trust radius
        integer, intent(inout)         :: xsave_lu
        !> Size of reduced space when equations have converged using Frobenius norm as trust radius
        integer, intent(inout)         :: redspacedim_save
        !> Max number of vectors to be kept in subspace for CROP solver
        integer, intent(in)            :: maxvec
        !> Set to true if level shift should be called on exit
        logical, intent(out)           :: do_LS
        !> Set to true if linear equations have converged and we are done
        logical, intent(out)           :: done
        real(realk)                    :: xnorm, maxelm
        type(matrix)                   :: xsave
        integer                        :: m, n
        logical                        :: OnMaster

   do_LS = .false.
   done = .false.

   if (it == 1) then
      errsave = err
   endif

   if (arh%info_lineq) then
      xnorm = sqrt(mat_sqnorm2(x))
      call mat_max_elm(x, maxelm)
      write (arh%lupri, '("  Error of red space iteration", i3, ":  ", F18.10, " (mu = ", F10.2, ", xmax = ", &
            & E14.6, ", xnorm = ", E14.6, ", THR = ", E14.6, ")")') it, err, mu, maxelm, xnorm, thresh
   endif

   if (err < thresh .or. err < 1.0E-9_realk) then
      call mat_max_elm(x, maxelm)
      if (j == 0 .AND. ABS(mu) > 1.0E-2_realk .and. err > 1.0E-9_realk .and. .not. arh%set_local &
          & .and. .not. arh%trustradius_decreased) then 
         arh%set_optxelm = .true.
         j = 1
         k = it
         mumax = mu
         xsave_lu = -1
         CALL LSOPEN(xsave_lu,'xsave','UNKNOWN','UNFORMATTED')
         OnMaster = .TRUE.
         call mat_write_to_disk(xsave_lu,x,OnMaster)
         if (arh%cfg_arh_truncate) then
            redspacedim_save = maxvec
         else
            redspacedim_save = it
         endif
         if (arh%info_lineq) write(arh%lupri,*) 'mumax =', mumax
      else if (maxelm > 1.1E0_realk*arh%set_max_element) then
         if (arh%info_lineq) write (arh%lupri,*) 'xmax', maxelm
         if (arh%lshift_by_hlgap) then
            write (arh%lupri,'("Newton equations converged in", i3, " iterations!")') it; arh%newton_it=it
            write (arh%lupri,'("- but scale solution since it is too large...")')
            call mat_scal(arh%set_max_element/maxelm, x)
            !Test that scaling is correct:
            call mat_max_elm(x, maxelm)
            if (arh%info_lineq) write (arh%lupri,*) 'xmax after scaling:', maxelm
            done = .true.
         else if (arh%set_optxelm) then 
            write (arh%lupri,*) 'Convergence, but xmax too large. No way to recover!'
            if (arh%cfg_arh_truncate) then
               write (arh%lupri,*) 'I recommend increasing the number of microvectors'
               write (arh%lupri,*) 'using:'
               write (arh%lupri,*) '.MICROVECS'
               write (arh%lupri,*) '<Number of desired microvectors>'
               write (arh%lupri,*) 'Currently set to:', arh%cfg_arh_microvecs
               call lsquit('Increase number of microvectors',arh%lupri)
            else
               write (arh%lupri,*) 'Either code is buggy or your example is really weird!'
               call lsquit('Augmented Roothaan-Hall equations failed to converge!',arh%lupri)
            endif
         else 
            if (arh%info_lineq) write(arh%lupri,*) 'mu found =', mumax
            write (arh%lupri,*)  "mu = 0 but xmax too large"
            call mat_max_elm(x, maxelm)
            arh%set_max_element = maxelm
            if (arh%cfg_arh_truncate) then
               arh%final_redspacedim = maxvec
            else
               arh%final_redspacedim = it
            endif
            done = .true.
         endif
      else 
         if (arh%info_lineq) then
            if (j==0 .and. .not. arh%set_local .and. maxelm/arh%set_max_element > 0.85E0_realk .and. &
            & abs(mu) > 1.0E-2_realk) then
               write (arh%lupri,*) 'No search in xmax attempted - size of X already at TR limit...'
            endif
            write (arh%lupri,'("Newton equations converged in", i3, " iterations!")') it; arh%newton_it=it
            if (.not. arh%set_local) arh%set_optxelm = .true. !For printing meaningful output
         endif
         if (arh%cfg_arh_truncate) then
            arh%final_redspacedim = maxvec
         else
            arh%final_redspacedim = it
         endif
         done = .true.
      endif
   endif 
   if (it == 1) then
      errsave = err
   endif
   call mat_max_elm(x, maxelm)
   if (it == k + arh%nits_check .and. .not. done) then
      k = it
      if (arh%info_lineq) then
         write (arh%lupri, '("Error went down by", F8.2, "%")') 100.0E0_realk - (err/errsave)*100E0_realk
      endif
      if (err > arh%error_decrease*errsave .or. (maxelm > 1.1E0_realk*arh%set_max_element .and. arh%lshift_by_hlgap)) then
          if (arh%lshift_by_hlgap) then
             if (err > arh%error_decrease*errsave) then
                write (arh%lupri,*) 'Residual is not going down, increase level shift'
             else
                write (arh%lupri,*) 'Solution vector too large, increase level shift'
             endif
             write (arh%lupri,*) 'Current level shift is', mu
             arh%cfg_fixed_shift = .true.
             !FIXME: make the 2E0_realk a variable
             if (abs(mu) > 1.0E-2_realk) then
                arh%cfg_fixed_shift_param = abs(mu*2E0_realk)
             else
                if (abs(maxelm) > 0.1) then
                   arh%cfg_fixed_shift_param = 0.5E0_realk
                else
                   arh%cfg_fixed_shift_param = 0.1E0_realk
                endif
             endif
             write (arh%lupri,*) 'New level shift set to:', -arh%cfg_fixed_shift_param
          else
             if (arh%set_optxelm) then
                write (arh%lupri,*) 'Residual is not going down, decrease trust radius (xmax)'
                call mat_max_elm(x, maxelm)
                arh%set_max_element = maxelm*arh%contract_factor
                write (arh%lupri,*) 'New trust radius (max element) =', arh%set_max_element
                arh%trustradius_decreased = .true.
             else
                write (arh%lupri,*) 'Residual is not going down, decrease trust radius (xnorm)'
                arh%trustradius_decreased = .true.
                xnorm = sqrt(mat_sqnorm2(x))
                arh%set_max_step = xnorm*arh%contract_factor
                write (arh%lupri,*) 'New trust radius (max step)=', arh%set_max_step
                call mat_max_elm(x, maxelm)
                arh%set_max_element = maxelm*arh%contract_factor
                write (arh%lupri,*) 'New trust radius (max element)=', arh%set_max_element
             endif
          endif
          do_LS = .true.
      endif
   endif
   if (it == k + 1) then
      errsave = err
   endif
   end subroutine arh_test_convergence

   !> \brief Set up ARH metric: Metr_ij = Tr(Di - D0)(Dj - D0)
   !> \author S. Host
   !> \date 2007
   subroutine fifo_inverse_metric(arh,fifoqueue)
   implicit none
        !> Contains solver info (ARH/TrFD)
        type(solverItem), intent(inout)         :: arh
        !> Contains Fock/KS and density matrices from previous SCF iterations
        TYPE(modFIFO),intent(inout)             :: fifoqueue
        integer                                 :: i, j, ndens, ndim, k
        real(realk),pointer                 :: scr1(:), scr2(:), vec(:), testmetric(:,:)
        type(matrix), pointer                   :: Di, Dj, D0, dummy
        type(matrix)                            :: idiff, jdiff
        !real(realk),allocatable :: metric(:,:), unitmat(:,:), test(:,:) for debug
        interface
           subroutine dsyevx_interface(A,ndim,print_eivecs,lupri,string)
             use precision
             implicit none
             character(6), intent(in), optional :: string
             integer, intent(in)      :: ndim, lupri
             real(realk), intent(in)  :: A(ndim,ndim)
             logical, intent(in)      :: print_eivecs 
           end subroutine dsyevx_interface
        end interface

   if (.not. arh%set_arhterms) then
      !nothing
   else
      ndens = fifoqueue%offset
      ndim = fifoqueue%D_exp%nrow
      call mat_init(idiff, ndim, ndim)
      call mat_init(jdiff, ndim, ndim)

      !DEBUG:
      !allocate(metric(ndens,ndens))
      !allocate(unitmat(ndens,ndens))
      !allocate(test(ndens,ndens))
      !END DEBUG 
      call mem_alloc(scr1,ndens*ndens)
      call mem_alloc(scr2,ndens*ndens)
      call mem_alloc(vec,ndens)

      !Set up modFIFO metric: Metr_ij = Tr(Di - D0)(Dj - D0)
      call get_from_modFIFO(fifoqueue, 0, dummy, D0)
      do i = 1, ndens
         call get_from_modFIFO(fifoqueue, i, dummy, Di)
         call MAT_ADD(1.0E0_realk, Di, -1.0E0_realk, D0, idiff) !Di = Di - D0
         do j = 1, ndens
            call get_from_modFIFO(fifoqueue, j, dummy, Dj)
            call MAT_ADD(1.0E0_realk, Dj, -1.0E0_realk, D0, jdiff) !Dj = Dj - D0
            !arh%fifometric(i,j) = mat_TrAB(idiff,jdiff)
            !We exploit symmetry to write
            arh%fifometric(i,j) = mat_dotproduct(idiff,jdiff)
         enddo
      enddo
      !write (arh%lupri,*) 'Metric, fifoqueue:'
      !call LS_OUTPUT(arh%fifometric, 1, ndens, 1, ndens, ndens, ndens, 1, arh%lupri)
      arh%inv_fifometric = arh%fifometric
      call INVERT_BY_DIAG2(arh%inv_fifometric,scr1,scr2,vec,ndens)
      !write (arh%lupri,*) 'Inverse metric, fifoqueue:'
      !call LS_OUTPUT(arh%inv_fifometric, 1, ndens, 1, ndens, ndens, ndens, 1, arh%lupri)

      if (.false.) then
         write(arh%lupri,*) 'Debug: test for linear dependencies'
         call mem_alloc(testmetric,ndens,ndens)
         do i = 1, ndens
            do j = 1, ndens
               testmetric(i,j) = (1.0E0_realk/sqrt(arh%fifometric(i,i))) &
                   & * arh%fifometric(i,j) * (1.0E0_realk/sqrt(arh%fifometric(j,j)))
            enddo
         enddo
         k = ndens
         !write(arh%lupri,*) 'Scaled metric:'
         !CALL LS_OUTPUT(testmetric, 1, k-1, 1, k-1, k-1, k-1, 1,arh%lupri)
         call dsyevx_interface(testmetric,k,.false.,arh%lupri)
         call mem_dealloc(testmetric)
      endif

      call mem_dealloc(scr1)
      call mem_dealloc(scr2)
      call mem_dealloc(vec)
      call mat_free(idiff)
      call mat_free(jdiff)
   endif
   end subroutine fifo_inverse_metric  

   !> \brief Get the weights of previous Fock/densities in the new D
   !> \author S. Host
   !> \date 2007
   subroutine arh_get_weights(arh,Dnew,x,fifoqueue,weights)
   implicit none
        !> Contains solver info (ARH/TrFD)
        type(solverItem), intent(in)        :: arh
        !> New density matrix D(X)
        type(Matrix), intent(in)            :: Dnew
        !> X from which new density is constructed, D(X)
        type(Matrix), intent(in)            :: x
        !> Contains Fock/KS and density matrices from previous SCF iterations
        TYPE(modFIFO),intent(inout)         :: fifoqueue
        !> The weights of previous Fock/densities in the new density
        real(realk), intent(inout)          :: weights(:)
        type(Matrix)                        :: wrk1, wrk2, com
        type(Matrix),pointer                :: D0, Dj, dummy 
        integer                             :: j, ndens, ndim
        real(realk)                         :: last_weight
        real(realk), pointer                :: vec(:)
        real(realk), pointer                :: T(:,:), Thalf(:,:), Tminushalf(:,:), scr(:)

   if (.not. arh%set_arhterms) then
      !nothing
   else
      ndens = fifoqueue%offset
      ndim = Dnew%nrow
      call mat_init(wrk1,ndim,ndim)
      call mat_init(wrk2,ndim,ndim)
      call mem_alloc(vec,ndens)

      if (arh%info_arh_weights) then
         call mat_init(com,ndim,ndim)
         !Construct weights using [D0,X] instead of total D(X) - for debugging.
         call mat_mul(fifoqueue%D_exp,x,'n','n',1E0_realk,0E0_realk,wrk1)
         call mat_trans(wrk1,wrk2)
         call MAT_ADD(1.0E0_realk, wrk1, 1.0E0_realk, wrk2, com) !com = DX - XD 

         call get_from_modFIFO(fifoqueue, 0, dummy, D0)
         do j = 1, ndens
            call get_from_modFIFO(fifoqueue, j, dummy, Dj)
            call MAT_ADD(1.0E0_realk, Dj, -1.0E0_realk, D0, wrk2)
!            vec(j) = mat_TrAB(com,wrk2)
            !We exploit that wrk2 is symmetric to do
            vec(j) = mat_dotproduct(com,wrk2)
         enddo

         last_weight = 0.0E0_realk
         weights = matmul(arh%inv_fifometric,vec)
         !write (arh%lupri,*) 'Weights using [D0,X]:'
         do j = 1, ndens
            !write (arh%lupri,'(F12.6)') weights(j)
            !Find 'last weight' : 1 - sum_i c_i
            last_weight = last_weight + weights(j)
         enddo
         last_weight = 1 - last_weight
         !write (arh%lupri,'(F12.6, "   (last weight)")') last_weight
         !write (arh%lupri,*)
         call mat_free(com)
         if (ndens > 0) then !Multiply weights by T^-1/2
            call mem_alloc(scr,2*ndens**2 + ndens*(ndens+1)/2)
            call mem_alloc(T,ndens,ndens)
            call mem_alloc(Thalf,ndens,ndens)
            call mem_alloc(Tminushalf,ndens,ndens)
            !T = queue%metric(1:ndens,1:ndens)
            T = arh%fifometric

            call LS_SQRTMT(T,ndens,2,Thalf,Tminushalf,scr)

            vec = matmul(Thalf,weights)
            !write (arh%lupri,*) 'Weights using [D0,X] multiplied by T^1/2:'
            last_weight = 0.0E0_realk
            do j = 1, ndens
               !write (arh%lupri,'(F12.6)') vec(j)
               !Find 'last weight' : 1 - sum_i c_i
               last_weight = last_weight + vec(j)
            enddo
            last_weight = 1 - last_weight
            !write (arh%lupri,'(F12.6, "   (last weight)")') last_weight
            !write (arh%lupri,*)
            !Do reverse operation to test:
            !weights = matmul(Tminushalf,vec)
            !write (arh%lupri,*) 'Weights using [D0,X] multiplied by T^-1/2*T^1/2: '
            !last_weight = 0.0E0_realk
            !do j = 1, ndens
            !   write (arh%lupri,'(F12.6)') weights(j)
            !   !Find 'last weight' : 1 - sum_i c_i
            !   last_weight = last_weight + weights(j)
            !enddo
            !last_weight = 1 - last_weight
            !write (arh%lupri,'(F12.6, "   (last weight)")') last_weight
            !write (arh%lupri,*)

            call mem_dealloc(scr)
            call mem_dealloc(T)
            call mem_dealloc(Thalf)
            call mem_dealloc(Tminushalf)
         endif
      endif

      call get_from_modFIFO(fifoqueue, 0, dummy, D0)
      do j = 1, ndens
         call MAT_ADD(1.0E0_realk, Dnew, -1.0E0_realk, D0, wrk1)
         call get_from_modFIFO(fifoqueue, j, dummy, Dj)
         call MAT_ADD(1.0E0_realk, Dj, -1.0E0_realk, D0, wrk2)
         !vec(j) = mat_TrAB(wrk1,wrk2)
         !We exploit that wrk1 and wrk2 is symmetric to do
         vec(j) = mat_dotproduct(wrk1,wrk2)         
      enddo
      weights = matmul(arh%inv_fifometric,vec)
      if (arh%info_arh_weights) write (arh%lupri,*) 'FIFO Weights using full D(X):'
      last_weight = 0.0E0_realk
      do j = 1, fifoqueue%offset
         last_weight = last_weight + weights(j)
         if (arh%info_arh_weights) write (arh%lupri,'(i5, F12.6)') j, weights(j)
      enddo
      last_weight = 1 - last_weight
      if (arh%info_arh_weights) write (arh%lupri,'(i5, F12.6, "   (FIFO last weight)")') fifoqueue%offset+1, last_weight
      if (arh%info_arh_weights) write (arh%lupri,*)
   
      call mat_free(wrk1)
      call mat_free(wrk2)
      call mem_dealloc(vec)
   endif
   end subroutine arh_get_weights   

   !> \brief Set up ndens*ndens matrix needed for ARH: M_ij = 0.5*Tr((D_i - D0)*(F(D_j)-F(D0)) + (D_j - D0)*(F(D_i)-F(D0)))
   !> \author S. Host
   !> \date 2007
   subroutine arh_get_M(arh,fifoqueue)
   implicit none
      !> Contains solver info (ARH/TrFD)
      type(SolverItem), intent(inout) :: arh
      !> Contains Fock/KS and density matrices from previous SCF iterations
      TYPE(modFIFO),intent(inout) :: fifoqueue
      type(Matrix)         :: FiminusF0, FjminusF0, DiminusD0, DjminusD0, Yi, Yj, scr1, scr2
      type(Matrix),pointer :: D0, F0, Di, Dj, Fi, Fj
      real(realk)          :: tr1, tr2
      integer              :: i, j, ndens, ndim
 
      !Extra Augmented Roothaan-Hall transformed matrices:
      !    M_ij = 1/2 * Tr((D_i - D0)*(F(D_j)-F(D0)) + (D_j - D0)*(F(D_i)-F(D0)) )

      ndim = fifoqueue%D_exp%nrow
      ndens = fifoqueue%offset

      call mat_init(FiminusF0,ndim,ndim)
      call mat_init(FjminusF0,ndim,ndim)
      call mat_init(DiminusD0,ndim,ndim)
      call mat_init(DjminusD0,ndim,ndim)

     call get_from_modFIFO(fifoqueue, 0, F0, D0)

     do i = 1, ndens
        call get_from_modFIFO(fifoqueue, i, Fi, Di)
        call MAT_ADD(1.0E0_realk, Di, -1.0E0_realk, D0, DiminusD0) !Di = Di - D0
        call MAT_ADD(1.0E0_realk, Fi, -1.0E0_realk, F0, FiminusF0) !Fi = Fi - F0
        do j = 1, ndens
           call get_from_modFIFO(fifoqueue, j, Fj, Dj)
           call MAT_ADD(1.0E0_realk, Fj, -1.0E0_realk, F0, FjminusF0) !Dj = Dj - D0
           call MAT_ADD(1.0E0_realk, Dj, -1.0E0_realk, D0, DjminusD0) !Fj = Fj - D0
!           tr1 = mat_TrAB(DiminusD0,FjminusF0)
!           tr2 = mat_TrAB(DjminusD0,FiminusF0)
           !we exploit symmetry of matrices to do 
           tr1 = mat_dotproduct(DiminusD0,FjminusF0)
           tr2 = mat_dotproduct(DjminusD0,FiminusF0)
           arh%fifoM(i,j) = tr1 + tr2
        enddo
     enddo

     arh%fifoM = 0.5E0_realk*arh%fifoM

     !write (lupri,*) 'FIFO M:'
     !call LS_OUTPUT(arh%fifoM, 1, ndens, 1, ndens, ndens, ndens, 1, arh%lupri)

      call mat_free(FiminusF0)
      call mat_free(FjminusF0)
      call mat_free(DiminusD0)
      call mat_free(DjminusD0)
   end subroutine arh_get_M

   !> \brief Multiply inverse ARH metric by vector
   !> \author S. Host
   !> \date 2007
   subroutine fifo_inv_metric_times_vector(arh,ndens,Vin, Vout)
   implicit none
        !> Contains solver info (ARH/TrFD)
        type(solverItem),intent(inout) :: arh
        !> Number of previous Fock/densities
        integer, intent(in)      :: ndens
        !> Input vector to be multiplied with inverse ARH metric 
        real(realk), intent(in)  :: Vin(:)
        !> Ouput vector = inverse ARH metric * input vector 
        real(realk), intent(inout) :: Vout(:)
        real(realk), pointer     :: A(:,:)
        integer, pointer         :: IPIV(:)
        integer                  :: INFO
   INFO=0

   if (.true.) then
      call mem_alloc(A,ndens,ndens)
      call mem_alloc(IPIV,ndens)
      A = arh%fifometric
      call DGESV(ndens, 1, A, ndens, IPIV, Vin, ndens, INFO)
      if (INFO /= 0) then
         write(arh%lupri,*) 'Error in DGESV =', INFO
         STOP 'Error in DGESV'
      endif
      Vout = Vin
      call mem_dealloc(A)
      call mem_dealloc(IPIV)
   else
      Vout = matmul(arh%inv_fifometric,Vin)
   endif
   end subroutine fifo_inv_metric_times_vector

   !> \brief Preconditioning of residual in linear equations
   !> \author S. Host
   !> \date 2007
   subroutine arh_precond(arh,decomp,res,symm,mu,res_prec)
   implicit none
         !> Contains solver info (ARH/TrFD)
         type(solverItem)         :: arh
         !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
         type(decompItem)         :: decomp
         !> Residual to be preconditioned
         TYPE(matrix), intent(in) :: res 
         !> Symmetry of residual: 1 = symmetric, 2 = antisymmetric
         integer, intent(in)      :: symm
         !> Level shift
         real(realk), intent(in)  :: mu               
         !> Preconditioned residual - output
         TYPE(matrix), intent(inout) :: res_prec 
         TYPE(matrix)                :: scr, scr2

   if (arh%cfg_noprec) then
      call mat_assign(res_prec,res)
   else if (arh%debug_arh_precond .and. .not. arh%set_local) then !PCG preconditioning in local area may force a saddlepoint solution
      call mat_init(scr, res%nrow, res%ncol)
      call mat_init(scr2, res%nrow, res%ncol)

      arh%debug_arh_precond = .false. !Otherwise this routine is called recursively!
      arh%set_arhterms = .false.

      call project_oao_basis(decomp, res, symm, scr)
      call arh_PCG(arh, decomp, scr, scr2, mu, symm)
      call project_oao_basis(decomp, scr2, symm, res_prec)

      arh%debug_arh_precond = .true.
      arh%set_arhterms = arh%cfg_arhterms

      call mat_free(scr)
      call mat_free(scr2)
   else
      call mat_init(scr, res%nrow, res%ncol)

      call mat_assign(scr,res)
      call mat_ao_precond(symm,mu,decomp%FUP,decomp%FUQ,decomp%DU,scr)
      call project_oao_basis(decomp, scr, symm, res_prec)

      call mat_free(scr)
   endif
   end subroutine arh_precond

   !> \brief Preconditioned Conjugate Gradient solver for linear equations (A-omega*I)x = b
   !> \author S. Host
   !> \date 2005
   subroutine arh_PCG(arh, decomp, b, x, omega, symm, fifoqueue)
   implicit none
         !> Contains solver info (ARH/TrFD)
         type(solverItem)         :: arh
         !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
         type(decompItem)         :: decomp
         !> Right hand side in linear equations (gradient)
         TYPE(matrix), intent(in) :: b
         !> Solution to system of linear equations (output)
         TYPE(matrix),intent(inout) :: x 
         !> Level shift
         real(realk), intent(in)  :: omega               
         !> Symmetry of matrices, 1 = symmetric, 2 = antisymmetric
         integer, intent(in)      :: symm
         !> Contains Fock/KS and density matrices from previous SCF iterations 
         TYPE(modFIFO),intent(inout),optional    :: fifoqueue !Not used if we do preconditioning
         integer                  :: N, j, i, k
         real(realk)              :: alfa, beta, err, alfa_num, alfa_den, test, &
                                   & beta_num, beta_den, gradnorm, thresh, norm, first_err
         TYPE(matrix)             :: p, Hr, Hr_new, res, scr
         integer                  :: ndim, maxit
         logical                  :: queue

   if (present(fifoqueue)) then
      queue = .true.
   else
      queue = .false.
   endif

   if (.not. queue .and. arh%set_arhterms) then
      WRITE(arh%LUPRI,'(/A)') &
      &     'ARH terms requested in arh_PCG but queue not present!'
      CALL lsQUIT('Queue not present in subroutine arh_PCG',arh%lupri)
   endif

   if (arh%set_arhterms) then
      maxit = 4
   else
      maxit = 100
   endif

   ndim = b%nrow
   gradnorm = SQRT(mat_sqnorm2(b))
   call MAT_INIT(p,ndim,ndim)
   call MAT_INIT(Hr,ndim,ndim)
   call MAT_INIT(Hr_new,ndim,ndim)
   call MAT_INIT(scr,ndim,ndim)
   call MAT_INIT(res,ndim,ndim)
   !Threshold of convergence:
   !Initialize iteration number and x:
   N=0
   call mat_assign(x,b)
   thresh = gradnorm*arh%cfg_micro_thresh
   !thresh = 1.0E-8_realk
   !Solving:
   do
      if (N == 10 + 1 .and. decomp%pcg_preconditioning) then
         write(decomp%lupri,'("                    Preconditioning: Force exit from PCG after", i4, " iterations")') N
         exit
      else if (N == maxit) then
         write(decomp%lupri,'("                    Force exit from PCG after", i4, " iterations")') N
         exit
      endif
      !If first iteration, initialize:
      if (N == 0) then
         if (queue) then
            call arh_lintrans(arh,decomp,x,symm,omega,scr,fifoqueue)
         else
            call arh_lintrans(arh,decomp,x,symm,omega,scr)
         endif
         norm = SQRT(mat_sqnorm2(scr))
         call MAT_ADD(1.0E0_realk, b, -1.0E0_realk, scr, res)    !res = b - A*x
         !first_err = SQRT(mat_sqnorm2(res))
         call arh_precond(arh,decomp,res,symm,omega,p)
         call mat_assign(Hr,p)
      else
         beta_den = mat_dotproduct(res, Hr)
         if (queue) then
            call arh_lintrans(arh,decomp,p,symm,omega,scr,fifoqueue)
         else
            call arh_lintrans(arh,decomp,p,symm,omega,scr)
         endif
         alfa_num = beta_den
         alfa_den = mat_dotproduct(p, scr)
         if (ABS(alfa_den) < 1.0E-15_realk*gradnorm) then
            write (decomp%lupri,*) "Warning: Exit from PCG loops due to denominator = 0!"
            if (N == 1) then  !fix made 1/11-05
               !call mat_copy(1.0E0_realk,b,x)
               call mat_assign(x,b) !Stinne 14/9-2010
            endif
            EXIT
         endif
         alfa = alfa_num/alfa_den
         !x = x + alfa*p
         call mat_daxpy(alfa, p, x)
         !res = res - alfa*VEC     !r(n+1) = r(n) - alfa*A*p(n)
         call mat_daxpy(-alfa, scr, res)
         !Convergence of PCG?
         err = SQRT(mat_sqnorm2(res))
         if (N == 1) first_err = err
         if (arh%info_lineq) then
            write (arh%LUPRI, '("                    Error of PCG iteration", i3, ":  ", F16.10, " (mu = ", F6.2, ")")') &
            & N, err, omega
         endif
         if (err < thresh) then
            write (arh%LUPRI, '("                    PCG converged in ", i3, " iterations!")') N
            EXIT
         endif
         call arh_precond(arh,decomp,res,symm,omega,Hr_new)
         beta_num = mat_dotproduct(res, Hr_new)
         beta = beta_num/beta_den
         !p = Hr_new + beta*p
         call mat_scal(beta,p) !Scale p by beta.
         call mat_daxpy(1.0E0_realk, Hr_new, p) !p = Hr_new + beta*p
         call mat_assign(Hr,Hr_new)
      endif
      N = N + 1 
   enddo

   !Test how much residual went down - if less than a factor of 10, damping is too small
   test = (first_err - err)/first_err
   write (arh%LUPRI, '("                    PCG: residual reduced by ", F9.2, "%  (mu = ", F6.2, ")")') &
                  & test*100.0E0_realk, omega 

   call MAT_FREE(p)
   call MAT_FREE(Hr)
   call MAT_FREE(Hr_new)
   call MAT_FREE(scr)
   call MAT_FREE(res)
   end subroutine arh_PCG

   !> @callergraph
   !> \brief Linear transformation for TrFD/Second order opt/ARH
   !> \author S. Host
   !> \date 2005
   !>
   !> \callgraph
   !>
   subroutine arh_lintrans(arh,decomp,x_in,symm,mu,AX,fifoqueue)
   implicit none  
         !> Contains solver info (ARH/TrFD)
         type(solverItem)          :: arh
         !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
         type(decompItem)          :: decomp
         !> Input matrix to be linearly transformed
         TYPE(matrix)              :: x_in
         !> Symmetry of input matrix: 1 = symmetric, 2 = antisymmetric
         integer, intent(in)       :: symm 
         !> Level shift
         real(realk), intent(in)   :: mu
         !> Linearly transformed vector (output)
         TYPE(matrix)              :: AX !Intent(out)
         !> Contains Fock/KS and density matrices from previous SCF iterations
         TYPE(modFIFO),intent(inout),optional  :: fifoqueue
         integer                   :: i, j, k, ndens, nbast
         TYPE(matrix)              :: scr1, scr2, D_AO, Gxc(1)
         real(realk)               :: err, testnorm, norm1, norm2
         TYPE(matrix)              :: FX, XF(1), G_xc, Dmat, scr1_mat(1)

      if (.not. present(fifoqueue) .and. arh%set_arhterms .and. &
        & .not. arh%set_do_2nd_order) then
         WRITE(arh%LUPRI,'(/A)') &
         &     'ARH terms requested in arh_lintrans but queue not present!'
         CALL lsQUIT('Queue not present in subroutine arh_lintrans',arh%lupri)
      endif
      !write (decomp%lupri,*) 'x, linear transformation, arh:'
      !call MAT_PRINT(x_in, 1, AX%nrow, 1, AX%ncol, LUPRI)

      call MAT_INIT(scr1,X_in%nrow,X_in%ncol)
      
      !TEST: 1st test proj
      !call project_oao_basis(x_in, symm, scr1)
      !write (decomp%lupri,*) 'scr1, linear transformation, arh, after proj:'
      !call MAT_PRINT(scr1, 1, AX%nrow, 1, AX%ncol, LUPRI)
      !x_in = scr1
      !END TEST

      call MAT_ADD(1.0E0_realk,decomp%FUQ,-1.0E0_realk,decomp%FUP,scr1) !scr1= FUQ - FUP

      call MAT_MUL(scr1,x_in,'n','n',1.0E0_realk,0.0E0_realk,AX) !AX = (FUQ-FUP)*X
      call mat_trans(AX,scr1)
      if (symm == arh_antisymmetric) then
         call mat_scal(-1.0E0_realk,scr1)
      endif                                          !scr1 = X*(FUQ-FUP)

      call mat_daxpy(1.0E0_realk,scr1,AX)    !AX = (FUQ-FUP)*X + X*(FUQ-FUP)
      call mat_daxpy(-mu,x_in,AX)      !AX = (FUQ-FUP)*X + X*(FUQ-FUP) - mu*X
      !testnorm = sqrt(mat_sqnorm2(AX))
      !write(decomp%lupri,*) 'Norm of TrFD contrib to lintrans:', testnorm 
      if (arh%debug_dd_lintra) then
         !write (decomp%lupri,*) 'Output from linear transformation, arh, TrFD contribution:'
         !call MAT_PRINT(AX, 1, AX%nrow, 1, AX%ncol, LUPRI)
         norm1 = sqrt(mat_sqnorm2(AX))
         write (decomp%lupri,*) 'Linear transf: norm of RH contribution is ', norm1
      endif

      if (arh%set_arhterms .and. .not. arh%set_do_2nd_order) then
         ndens = fifoqueue%offset
         if (ndens > 0) then
            call arh_xdep_matrices(arh,decomp,x_in,symm,fifoqueue,scr1) !scr is the
                                                                !arh part of the linear transformation
            !WRITE(decomp%lupri,*) 'arh contrib to lintrans:' 
            !call mat_print(scr1,1,FUP%nrow,1,FUP%ncol,decomp%lupri)
            call mat_daxpy(1.0E0_realk,scr1,AX) 
         endif
      endif

      if (arh%set_do_2nd_order) then !Add G(D) type contribution
         call MAT_INIT(FX,x_in%nrow,x_in%ncol)
         call MAT_INIT(XF(1),x_in%nrow,x_in%ncol)
         call MAT_INIT(scr2,X_in%nrow,X_in%ncol)
         call MAT_MUL(x_in,decomp%DU,'n','n',1.0E0_realk,0.0E0_realk,XF(1)) !XF is XD here!
         call mat_trans(XF(1),FX)
         if (symm == 2) then
            call mat_scal(-1.0E0_realk,FX)
         endif
         call MAT_ADD(1.0E0_realk,XF(1),-1.0E0_realk,FX,scr1)
         call x_from_oao_basis(decomp, scr1, XF(1)) !XF is XD_AO
         if (decomp%nocca+decomp%noccb /= 1) then !if only 1 electron -> 2 el part = 0
         
            !call di_GET_GbDs(decomp%lupri,decomp%lupri,XF(1),scr1) !scr1 is G_AO without DFT
            nbast = scr1%nrow
            !call mat_init(G_xc,nbast,nbast)
            !call mat_init(Dmat,nbast,nbast)
            call di_GET_GbDs_and_XC_linrsp(scr1, G_xc, decomp%lupri, decomp%lupri, XF(1), nbast, Dmat, .false.)
            !call mat_free(G_xc)
            !call mat_free(Dmat)
            !-------------------------------------------
			
         else
            call mat_zero(scr1)
         endif
         if(lsint_fock_data%ls%setting%do_dft) then
            !Add extra G XC contributions
            call MAT_INIT(Gxc(1),X_in%nrow,X_in%ncol)
            call mat_zero(Gxc(1))
            call MAT_INIT(D_AO,X_in%nrow,X_in%ncol)  
            call x_from_oao_basis(decomp,decomp%DU, D_AO)  ! density in AO basis
            call II_get_xc_linrsp(lsint_fock_data%ls%lupri,lsint_fock_data%ls%luerr,&
                    & lsint_fock_data%ls%setting,XF(1)%nrow,XF,D_AO,Gxc,1)  ! Get XC contribution to G
            call mat_daxpy(1E0_realk,Gxc(1),scr1)  ! Add to existing G (stored in scr1)
            call mat_free(D_AO)
            call mat_free(Gxc(1))
         end if
         call res_to_oao_basis(decomp, scr1, scr2) !scr2 is G in chol basis
         !Add contribution to AX, 1.0*(GU*DU -DU*GU)
         call MAT_MUL(scr2,decomp%DU,'n','n',1.0E0_realk,0.0E0_realk,XF(1)) !XF is GU*DU here!
         call mat_trans(XF(1),FX) !FX is DU*GU here (both D and G are always symm)
         call MAT_ADD(1.0E0_realk,XF(1),-1.0E0_realk,FX,scr1)
         call mat_daxpy(1.0E0_realk,scr1,AX)
         call project_oao_basis(decomp, AX, symm, scr1) 
         call mat_assign(AX,scr1)
         call MAT_FREE(scr2)
         call MAT_FREE(XF(1))
         call MAT_FREE(FX)
      endif

      if (arh%debug_dd_lintra) then
         !write (decomp%lupri,*) 'Total linear transformation, arh:'
         !call MAT_PRINT(AX, 1, AX%nrow, 1, AX%ncol, LUPRI)
         norm2 = sqrt(mat_sqnorm2(AX))
         write (decomp%lupri,*) 'Linear transf: total norm is              ', norm2
         write (decomp%lupri,*) 'Percentage of RH contribution in lintrans:', 100*norm1/norm2
      endif

      call MAT_FREE(scr1)
   end subroutine arh_lintrans

   !> \brief ARH part of linear transformation
   !> \author S. Host
   !> \date 2005
   subroutine arh_xdep_matrices(arh,decomp,x,symm,fifoqueue,Ax)
   implicit none
        !> Contains solver info (ARH/TrFD)
        type(solverItem),intent(inout) :: arh
        !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
        type(decompItem),intent(inout) :: decomp
        !> Input matrix to be linearly transformed
        type(Matrix), intent(in) :: x 
        !> Symmetry of input matrix: 1 = symmetric, 2 = antisymmetric
        integer, intent(in)      :: symm 
        !> Contains Fock/KS and density matrices from previous SCF iterations
        TYPE(modFIFO),intent(inout) :: fifoqueue
        !> Linearly transformed matrix (output)
        type(matrix), intent(inout) :: Ax
        integer                  :: i, ndens
        type(Matrix)             :: wrk1, wrk2
        type(Matrix),pointer     :: F0, D0, Fi, Di
        !real(realk), allocatable, dimension(:) :: V1, V2, scrvec
        !real(realk), allocatable, dimension(:) :: vec1, vec2, vec3
        real(realk), pointer    :: vec2(:), V2(:)
        real(realk) :: test, norm1, norm2, norm3, maxelm1, maxelm2, maxelm3 
        integer                  :: ndim

    ndim = x%nrow

    ndens = fifoqueue%offset

    !allocate(scrvec(ndens))
    !allocate(vec1(ndens), vec2(ndens), vec3(ndens)) 
    !allocate(V1(ndens),V2(ndens))
    call mem_alloc(vec2,ndens)
    call mem_alloc(V2,ndens)

    call mat_init(wrk1,ndim,ndim)
    call mat_init(wrk2,ndim,ndim)

    !Get projected input matrix = DX - XD:
    call MAT_MUL(decomp%DU,x,'n','n',1.0E0_realk,0.0E0_realk,wrk1)
    call mat_trans(wrk1,wrk2)
    if (symm == arh_antisymmetric) then
       call mat_scal(-1.0E0_realk,wrk2)
    endif
    call MAT_ADD(1.0E0_realk, wrk1, -1.0E0_realk, wrk2, Ax) !Ax = x projected, DX - XD

    ! Get the following
    !    V1_i = Tr((DX - XD)*(Fi-F0)) not used anymore!
    !    V2_i = Tr((DX - XD)*(Di-D0))

    call get_from_modFIFO(fifoqueue, 0, F0, D0)
    do i = 1, ndens
       call get_from_modFIFO(fifoqueue, i, Fi, Di)
       call MAT_ADD(1.0E0_realk, Di, -1.0E0_realk, D0, wrk1) !wrk1 = Di - D0
       !call MAT_ADD(1.0E0_realk, Fi, -1.0E0_realk, F0, wrk2) !wrk2 = Fi - F0
       !V1(i) = mat_trAB(wrk2,Ax)
!       V2(i) = mat_TrAB(wrk1,Ax)
       !we exploit symmetry to write
       V2(i) = mat_dotproduct(wrk1,Ax)
!       print*,'mat_get_isym(wrk1,1.0E-10_realk)',mat_get_isym(wrk1,1.0E-10_realk)
!       print*,'mat_get_isym(Ax,1.0E-10_realk)',mat_get_isym(Ax,1.0E-10_realk)
!       if(mat_get_isym(wrk1,1.0E-10_realk).NE.1)call lsquit('wrk1 C',-1)
!       if(mat_get_isym(Ax,1.0E-10_realk).NE.1)call lsquit('Ax C',-1)
    enddo

    if (arh%debug_arh_lintra) then
       if (ndens > 0) then 
          !write (lupri,*) 'V1:'
          !call LS_OUTPUT(V1, 1, ndens, 1, 1, ndens, 1, 1, arh%lupri)
          write (arh%lupri,*) 'V2:'
          call LS_OUTPUT(V2, 1, ndens, 1, 1, ndens, 1, 1, arh%lupri)
       endif
    endif

    ! Get the following
    !    vec1 = [T]^-1 * V1 not used anymore!
    !    vec2 = [T]^-1 * V2
    !    vec3 = [T]^-1 * M * [T]^-1 * V2 not used anymore!

    call fifo_inv_metric_times_vector(arh,ndens,V2,vec2)

    !scrvec = matmul(arh%fifoM,vec2)
    !call fifo_inv_metric_times_vector(arh,ndens,scrvec,vec3)

    if (arh%debug_arh_lintra) then
       if (ndens > 0) then 
          !write (arh%lupri,*) 'vec1 = [T]^-1 * V1'
          !call LS_OUTPUT(vec1, 1, ndens, 1, 1, ndens, 1, 1, arh%lupri)
          write (arh%lupri,*) 'vec2 = [T]^-1 * V2'
          call LS_OUTPUT(vec2, 1, ndens, 1, 1, ndens, 1, 1, arh%lupri)
          !write (arh%lupri,*) 'vec3 = [T]^-1 * M * [T]^-1 * V2'
          !call LS_OUTPUT(vec3, 1, ndens, 1, 1, ndens, 1, 1, arh%lupri)
       endif
    endif

    call MAT_ZERO(Ax)

    do i = 1, ndens
       call get_from_modFIFO(fifoqueue, i, Fi, Di) 
       call MAT_ADD(1.0E0_realk, Fi, -1.0E0_realk, F0, wrk1) !wrk1 = Fi - F0
       call mat_daxpy(vec2(i),wrk1,Ax)           !Ax = SUM_i (Fi -F0)*[T^-1*V2]_i
    enddo

    !Project:
    call MAT_MUL(decomp%DU,Ax,'n','n',1.0E0_realk,0.0E0_realk,wrk1)   !wrk1 = DU*Ax
    call MAT_MUL(wrk1,decomp%QU,'n','n',1.0E0_realk,0.0E0_realk,wrk2) !wrk2 = DU*Ax*QU
    call mat_trans(wrk2,wrk1)
    call MAT_ADD(1.0E0_realk, wrk2, -1.0E0_realk, wrk1, Ax)           !Ax = DU*Ax*QU - QU*Ax* DU 

    !deallocate(scrvec)
    !deallocate(V1,V2)
    !deallocate(vec1,vec2,vec3)
    call mem_dealloc(vec2)
    call mem_dealloc(V2)
    call mat_free(wrk1)
    call mat_free(wrk2)
   end subroutine arh_xdep_matrices

!###### Debug section ###############################
!!$
!!$   !> \brief Set up ARH or exact Hessian by linear transformation
!!$   !> \author S. Host
!!$   !> \date 2005
!!$   subroutine debug_get_hessian(arh,decomp,fifoqueue,hes)
!!$   implicit none
!!$           !> Contains solver info (ARH/TrFD)
!!$           type(SolverItem),intent(inout)    :: arh
!!$           !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
!!$           type(decompItem),intent(in)       :: decomp
!!$           !> Contains Fock/KS and density matrices from previous SCF iterations
!!$           TYPE(modFIFO),intent(inout)          :: fifoqueue
!!$           !> Hessian. Dimension should be nbas*(nbas+1)/2 - nbas (output)
!!$           TYPE(matrix),intent(inout)        :: hes 
!!$           integer                 :: m, l, vecdim, matdim
!!$           TYPE(matrix)            :: x_trial, x_trial_mat, column, scr
!!$
!!$      matdim = fifoqueue%D_exp%nrow
!!$      vecdim = matdim*(matdim+1)/2 - matdim
!!$
!!$      call MAT_INIT(scr,matdim,matdim)
!!$      call MAT_INIT(x_trial,vecdim,1)
!!$      call MAT_INIT(x_trial_mat,matdim,matdim)
!!$      call MAT_INIT(column,vecdim,1)
!!$
!!$      if (arh%debug_hessian_exact) then
!!$         arh%set_do_2nd_order = .true.
!!$         arh%set_arhterms     = .false.
!!$      endif
!!$
!!$      do m = 1, vecdim
!!$         call MAT_ZERO(x_trial)
!!$         call MAT_ZERO(x_trial_mat)
!!$         call mat_create_elm(m, 1, 1.0E0_realk, x_trial)
!!$         call mat_VEC_TO_MAT('a', x_trial, scr)
!!$      !write (LUPRI,*) "xmat:"
!!$      !call MAT_PRINT(scr, 1, scr%nrow, 1, scr%ncol, LUPRI)
!!$         call project_oao_basis(decomp, scr, arh_antisymmetric, x_trial_mat)  
!!$      !write (LUPRI,*) "xmat projected:"
!!$      !call MAT_PRINT(x_trial_mat, 1, x_trial_mat%nrow, 1, x_trial_mat%ncol, LUPRI)
!!$         call project_oao_basis(decomp, x_trial_mat, arh_antisymmetric, scr)
!!$      !write (LUPRI,*) "xmat projected again:"
!!$      !call MAT_PRINT(scr, 1, x_trial_mat%nrow, 1, x_trial_mat%ncol, LUPRI)
!!$         !column = m'th column of hessian:
!!$      !write(lupri,*) 'antisymmetric:', antisymmetric
!!$         call arh_lintrans(arh,decomp,x_trial_mat,arh_antisymmetric,0.0E0_realk,scr,fifoqueue)
!!$         call MAT_TO_VEC('a', scr, column)
!!$         do l = 1, vecdim !Put elements of column in m'th column of hessian
!!$            call mat_create_elm(l,m,column%elms(l),hes)
!!$         enddo         
!!$      enddo
!!$
!!$      call mat_scal(4.0E0_realk, hes)
!!$
!!$      if (arh%debug_hessian_exact) then
!!$         arh%set_do_2nd_order = arh%cfg_do_2nd_order
!!$         arh%set_arhterms     = arh%cfg_arhterms
!!$      endif
!!$
!!$      !write (LUPRI,*) "Hessian:"
!!$      !call MAT_PRINT(hes, 1, hes%nrow, 1, hes%ncol, LUPRI)
!!$
!!$      call MAT_FREE(x_trial)
!!$      call MAT_FREE(x_trial_mat)
!!$      call MAT_FREE(column)
!!$      call MAT_FREE(scr)
!!$   end subroutine debug_get_hessian

!####################################################################################
!  The following section contains routines for using the Conjugate Residual Optimal
!  Vectors (CROP) scheme (imposed with cfg_arh_crop) for the solution of the arh micro equations.
!  CROP:  M. Ziolkowski, V. Weijo, P. Jorgensen et al. JCP 128, 204105
!####################################################################################

   !> \brief Get x and residual in CROP scheme
   !> \author S. Host
   !> \date 2007
   subroutine arh_crop_x_and_res(arh,decomp,fifoqueue,G,symm,mu,Xn,sigma_n,res_new,res)
   implicit none
         !> Contains solver info (ARH/TrFD)
         type(solverItem), intent(inout)  :: arh
         !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
         type(decompItem), intent(in)  :: decomp
         !> Contains Fock/KS and density matrices from previous SCF iterations
         TYPE(modFIFO),intent(inout)    :: fifoqueue
         !> Gradient
         TYPE(matrix), intent(in)    :: G
         !> Symmetry of matrices. 1 = symmetric, 2 = antisymmetric
         integer, intent(in)         :: symm
         !> Level shift
         real(realk), intent(in)     :: mu
         !>  Input: Current optimal solution vector x_n par. Output: n+1'th solution vector
         TYPE(matrix), intent(inout) :: Xn     
         !>  Input: Current optimal linear transf. sigma_n par. Output: n+1'th linear transf. 
         TYPE(matrix), intent(inout) :: sigma_n
         !> Preconditioned current optimal residual P(r_n par) (may be either input or constructed here)
         TYPE(matrix), intent(inout) :: res_new 
         !> n+1'th residual (output)
         TYPE(matrix), intent(inout) :: res
         integer                     :: ndim

   if (arh%cfg_arh_crop_safe) then !Otherwise use incoming residual
      ndim = G%nrow
      !1. Preconditioning + projection of r_n par = G - sigma_n par
      call mat_add(-1E0_realk,sigma_n,-1E0_realk, G, res) !scr = r_n par
      call mat_daxpy(mu,Xn,res) 
      !Preconditioning + projection of r_n par
      call arh_precond(arh,decomp,res,symm,mu,res_new) !res_new = P(r_n par)
   endif

   !2. Calculate x_n+1 = x_n par - res_n par 
   call mat_daxpy(1.0E0_realk,res_new,Xn) !Xn is now x_n+1

   !3. Linear transformation of new solution vector:
   call arh_lintrans(arh,decomp,Xn,symm,0.0E0_realk,sigma_n,fifoqueue) !sigma_n is now sigma_n+1

   !4. Calculate residual r_n+1 = G - (H(Xn) + mu*Xn) 
   call mat_add(-1E0_realk,G,-1E0_realk, sigma_n, res)
   call mat_daxpy(mu,Xn,res)

   end subroutine arh_crop_x_and_res

   !> \brief Set up intermediate subspace in CROP scheme
   !> \author S. Host
   !> \date 2007
   subroutine arh_crop_intermed_sub(arh,decomp,lub,lusigma,vectorsubspace,n,symm,G,mu,r_new,x,sigma,xF,sigmaF)
   implicit none
         !> Contains solver info (ARH/TrFD)
         type(solverItem),intent(inout) :: arh
         !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
         type(decompItem),intent(in) :: decomp
         !> Logical unit number of file containing trial vectors. Not referenced if trial vectors are kept in core
         integer, intent(in)         :: lub
         !> Logical unit number of file containing sigma vectors. Not referenced if sigma vectors are kept in core
         integer, intent(in)         :: lusigma
         !> Subspace of previous trial and sigma vectors
         TYPE(modFIFO),intent(inout) :: vectorsubspace
         !> Size of subspace
         integer, intent(in)         :: n 
         !> Symmetry of matrices. 1 = symmetric, 2 = antisymmetric
         integer, intent(in)         :: symm
         !> Gradient
         TYPE(matrix), intent(in)    :: G
         !> Level shift
         real(realk), intent(in)     :: mu
         !> Newest residual r_n+1
         TYPE(matrix), intent(in)    :: r_new
         !> Input: n+1'th solution vector. Output: n+1'th optimal solution vector x_n+1 par
         TYPE(matrix), intent(inout) :: x
         !> Input:  n+1'th linear transf.  Output: n+1'th optimal linear transf. sigma_n+1 par 
         TYPE(matrix), intent(inout) :: sigma
         !> For alternative CROP level shift scheme (see description in arh_crop_optimal_x)
         TYPE(matrix), intent(inout) :: xF
         !> For alternative CROP level shift scheme (see description in arh_crop_optimal_x)
         TYPE(matrix), intent(inout) :: sigmaF
         real(realk), pointer        :: ResRed(:,:), RHS(:)
         TYPE(matrix)                :: scr, scr2, r_newP, xvec
         type(matrix), pointer       :: xpointer, sigmapointer
         integer, pointer            :: IPIV(:)
         integer                     :: i, rowdim, coldim, IERR, redsize, qsize, printsize, Asize, newsize
         real(realk), pointer    :: tempA(:,:), tempS(:,:)
         logical                     :: space_expanded,OnMaster
        interface
           subroutine dsyevx_interface(A,ndim,print_eivecs,lupri,string)
             use precision
             implicit none
             character(6), intent(in), optional :: string
             integer, intent(in)      :: ndim, lupri
             real(realk), intent(in)  :: A(ndim,ndim)
             logical, intent(in)      :: print_eivecs 
           end subroutine dsyevx_interface
        end interface
   IERR=0
   !ndim = G%nrow
   rowdim = G%nrow
   coldim = G%ncol
   if (arh%cfg_arh_truncate) then
      redsize = vectorsubspace%offset
   else
      redsize = n
   endif

   call mem_alloc(ResRed,redsize+2,redsize+2)
   !Set up intermediate 'residual subspace' with proper dimension:

   ResRed(1:redsize,1:redsize) = arh%CROPmat(1:redsize,1:redsize)

   if (arh%info_crop) then
      write (arh%lupri,*) 'Intermediate subspace before augmentation:'
      call LS_OUTPUT(ResRed, 1, redsize, 1, redsize, redsize+2, redsize+2, 1, arh%lupri)
      !call ls_flshfo(arh%lupri)
   endif

   !Preconditioning of r_n+1:
   call mat_init(r_newP, rowdim,coldim)
   call arh_precond(arh,decomp,r_new,symm,mu,r_newP)

   if (.not. arh%cfg_arh_truncate) then
!      call mat_init(xvec, rowdim,coldim)
      call mat_init(scr2, rowdim,coldim)
   endif
   call mat_init(scr, rowdim,coldim)

   !Diagonal elements:
   ResRed(redsize+1,redsize+1) = mat_dotproduct(r_new,r_newP)
   ResRed(redsize+2,redsize+2) = 0.0E0_realk
   if (.not. arh%cfg_arh_truncate) then
      rewind(lusigma) ; rewind(lub)
   endif
   do i = 1, redsize  !Setup n+1'th dimension
      if (arh%cfg_arh_truncate) then
         call get_from_modFIFO(vectorsubspace, i, xpointer, sigmapointer)
         !call get_sigma(lusigma,sigma_n)
         call mat_add(-1E0_realk,G,-1E0_realk, sigmapointer, scr) !scr = i'th residual
         call mat_daxpy(mu,xpointer,scr)
      else
         OnMaster = .TRUE.
         call mat_read_from_disk(lub,scr2,OnMaster)
         call mat_add(-1E0_realk,G,mu, scr2, scr) !scr = i'th residual 
         call mat_read_from_disk(lusigma,scr2,OnMaster)
         call mat_daxpy(-1E0_realk,scr2,scr)
         !call get_sigma(lusigma,sigma_n)
!         call mat_read_from_disk(lub,xvec,OnMaster)
!         call mat_read_from_disk(lusigma,scr2,OnMaster)
!         call mat_add(-1E0_realk,G,-1E0_realk, scr2, scr) !scr = i'th residual 
!         call mat_daxpy(mu,xvec,scr)
      endif

      ResRed(redsize+1,i) = mat_dotproduct(r_newP,scr)
      ResRed(i,redsize+1) = ResRed(redsize+1,i)
   enddo 
!   if (.not. arh%cfg_arh_truncate) then
!      call mat_free(xvec)
!   endif
   call mat_free(r_newP)
   do i = 1, redsize+1  !'Augmented dimension':
      ResRed(redsize+2,i) = -1.0E0_realk
      ResRed(i,redsize+2) = ResRed(redsize+2,i)
   enddo

   if (arh%info_crop) then
      write (arh%lupri,*) 'Intermediate subspace after augmentation:'
      call LS_OUTPUT(ResRed, 1, redsize+2, 1, redsize+2, redsize+2, redsize+2, 1, arh%lupri)
      !call ls_flshfo(arh%lupri)
   endif

   call mem_alloc(RHS,redsize+2)
   call mem_alloc(IPIV,redsize+2)
   RHS = 0.0E0_realk ; RHS(redsize+2) = -1.0E0_realk

   !Solve set of linear equations CROPmat*c = RHS:
   call DGESV(redsize+2, 1, ResRed, redsize+2, IPIV, RHS, redsize+2, IERR) !Solution vector is found in RHS.

   call mem_dealloc(ResRed)
   call mem_dealloc(IPIV)
   if (IERR /= 0) then
      WRITE(arh%LUPRI,'(/A, i4)') &
      &     'Problem in DGESV, IERR = ', IERR
      CALL lsQUIT(' Problem in DGESV',arh%lupri)
   endif

   if (arh%info_crop) then
      write (arh%lupri,*) 'Solution vector:'
      call LS_OUTPUT(RHS, 1, redsize+2, 1, 1, redsize+2, 1, 1, arh%lupri) ; call ls_flshfo(arh%lupri)
   endif

   !Construct x_n+1 par and sigma_n+1 par:
   !x_n+1 par = c_n+1* x_n+1 + sum(i=1,n) c_i * x_i par
   !sigma_n+1 par = c_n+1* sigma_n+1 + sum(i=1,n) c_i * sigma_i par
   call mat_scal(RHS(redsize+1),x) ; call mat_scal(RHS(redsize+1),sigma)
   if (.not. arh%cfg_arh_truncate) then
      rewind(lusigma) ; rewind(lub)
   endif
   do i = 1, redsize
      if (arh%cfg_arh_truncate) then
         call get_from_modFIFO(vectorsubspace, i, xpointer, sigmapointer)
         call mat_daxpy(RHS(i),xpointer,x)
         call mat_daxpy(RHS(i),sigmapointer,sigma)         
      else
         OnMaster=.TRUE.
         call mat_read_from_disk(lub,scr,OnMaster)
         call mat_read_from_disk(lusigma,scr2,OnMaster)
         call mat_daxpy(RHS(i),scr,x)
         call mat_daxpy(RHS(i),scr2,sigma)
      endif
   enddo
   call mem_dealloc(RHS)
   if (.not. arh%cfg_arh_truncate) then
      call mat_free(scr2)
   endif
   call mat_free(scr)

   if (arh%cfg_arh_truncate .and. arh%cfg_arh_newdamp) then
      qsize = vectorsubspace%queuesize-1
      !write(arh%lupri,*) 'qsize, redsize:', qsize, redsize
      !write(arh%lupri,*) 'vectorsubspace%offset, n', vectorsubspace%offset, n
      if (vectorsubspace%offset == qsize .and. n > qsize-1) then
         !space_expanded = .true.
         Asize = qsize+1
         !if (vectorsubspace%offset == qsize .and. n > qsize-1) space_expanded = .false.
         if (n == qsize) then
            !Special case: Subspace is full for the first time
            newsize = redsize
         else
            newsize = redsize+1
         endif
         if (arh%info_crop) then
            !if (space_expanded) then
            !   printsize = redsize-1
            !else
            !   printsize = redsize
            !endif
            write (arh%lupri,*) 'Ared subspace, get extra mats:'
            call LS_OUTPUT(arh%Ared, 1, newsize, 1, newsize, Asize, Asize, 1, arh%lupri)
            write (arh%lupri,*) 'Sred subspace, get extra mats:'
            call LS_OUTPUT(arh%Sred, 1, newsize, 1, newsize, Asize, Asize, 1, arh%lupri)
            write (arh%lupri,*) 'Gred subspace, get extra mats:'
            call LS_OUTPUT(arh%Gred, 1, newsize, 1, 1, Asize, 1, 1, arh%lupri)
         endif
         ! if n = qsize-1, it is the first time that offset = qsize and we do not truncate yet
         !Determine xF and sigmaF by diagonalizing current
         !redspace:
         call mem_alloc(tempA,newsize,newsize)
	 call mem_alloc(tempS,newsize,newsize)
         tempA(1:newsize,1:newsize) = arh%Ared(1:newsize,1:newsize)
         tempS(1:newsize,1:newsize) = arh%Sred(1:newsize,1:newsize)
         if (.false.) then
            call dsyevx_interface(tempS,newsize,.false.,arh%lupri,'Metric')
            tempS(1:newsize,1:newsize) = arh%Sred(1:newsize,1:newsize)
         endif
         call arh_crop_optimal_x(arh%lupri,tempA,tempS,newsize,vectorsubspace,xF,sigmaF)
         call mem_dealloc(tempA)
	 call mem_dealloc(tempS)
      endif
   endif

   end subroutine arh_crop_intermed_sub

   !> \brief Set up reduced space in CROP scheme
   !> \author S. Host
   !> \date 2007
   subroutine arh_crop_setup_redsp(arh,decomp,lub,lusigma,vectorsubspace,n,symm,G,mu,x,sigma,resP,scr,xF,sigmaF)
   implicit none
         !> Contains solver info (ARH/TrFD)
         type(solverItem),intent(inout) :: arh
         !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
         type(decompItem),intent(in) :: decomp
         !> Logical unit number of file containing trial vectors. Not referenced if trial vectors are kept in core
         integer, intent(in)         :: lub
         !> Logical unit number of file containing sigma vectors. Not referenced if sigma vectors are kept in core
         integer, intent(in)         :: lusigma
         !> Subspace of previous trial and sigma vectors
         TYPE(modFIFO),intent(inout) :: vectorsubspace
         !> Size of subspace
         integer, intent(in)         :: n
         !> Symmetry of matrices. 1 = symmetric, 2 = antisymmetric
         integer, intent(in)         :: symm
         !> Gradient
         TYPE(matrix), intent(in)    :: G
         !> Level shift
         real(realk), intent(in)     :: mu
         !>  n+1'th optimal solution vector, x_n+1 par
         TYPE(matrix), intent(in)    :: x  
         !>  n+1'th optimal sigma vector, sigma_n+1 par
         TYPE(matrix), intent(in)    :: sigma
         !>  New preconditioned optimal residual P(r_n+1 par) (output)
         TYPE(matrix), intent(inout) :: resP
         !> Work space
         TYPE(matrix), intent(inout) :: scr
         !> For alternative CROP level shift scheme (see description in arh_crop_optimal_x)
         TYPE(matrix), intent(inout) :: xF
         !> For alternative CROP level shift scheme (see description in arh_crop_optimal_x)
         TYPE(matrix), intent(inout) :: sigmaF
         TYPE(matrix)                :: scr2, xvec
         integer                     :: i, rowdim, coldim, redsize, qsize, newsize, Asize, printsize, min_eival_index(1), nsize
         type(matrix), pointer       :: xpointer, sigmapointer
         real(realk), pointer    :: tempvec(:), tempmat(:,:), tempA(:,:), tempS(:,:)
         real(realk)                 :: min_eival
         logical                     :: new_scheme, space_expanded,OnMaster !Temporary!!
        interface
           subroutine dsyevx_interface(A,ndim,print_eivecs,lupri,string)
             use precision
             implicit none
             character(6), intent(in), optional :: string
             integer, intent(in)      :: ndim, lupri
             real(realk), intent(in)  :: A(ndim,ndim)
             logical, intent(in)      :: print_eivecs 
           end subroutine dsyevx_interface
        end interface

         OnMaster = .TRUE.
   space_expanded = .true.
   !ndim = G%nrow
   rowdim = G%nrow
   coldim = G%ncol
   qsize = vectorsubspace%queuesize-1
   if (arh%cfg_arh_newdamp) then
      Asize = qsize+1
   else
      Asize = qsize
   endif
   if (arh%cfg_arh_truncate) then
      redsize = vectorsubspace%offset
      !write(arh%lupri,*) 'qsize, redsize:', qsize, redsize
      if (vectorsubspace%offset == qsize .and. n > qsize-1) space_expanded = .false.
      if (arh%info_crop) then
         if (space_expanded) then
            printsize = redsize-1
         else
            printsize = redsize
         endif
         write (arh%lupri,*) 'Residual subspace initial:'
         call LS_OUTPUT(arh%CROPmat, 1, printsize, 1, printsize, qsize, qsize, 1, arh%lupri)
         write (arh%lupri,*) 'Ared subspace initial:'
         call LS_OUTPUT(arh%Ared, 1, printsize, 1, printsize, Asize, Asize, 1, arh%lupri)
         write (arh%lupri,*) 'Sred subspace initial:'
         call LS_OUTPUT(arh%Sred, 1, printsize, 1, printsize, Asize, Asize, 1, arh%lupri)
         write (arh%lupri,*) 'Gred subspace initial:'
         call LS_OUTPUT(arh%Gred, 1, printsize, 1, 1, Asize, 1, 1, arh%lupri)
         !call ls_flshfo(arh%lupri)
      endif
      if (vectorsubspace%offset == qsize .and. n > qsize-1) then
         ! if n = qsize-1, it is the first time that offset = qsize and we do not truncate yet
         !Subspace is full and we start to overwrite old vectors
         !if (cfg_arh_newdamp) then
         !   !Start by determining xF and sigmaF by diagonalizing current
         !   !redspace:
         !   if (n == qsize) then
         !      !Special case: Subspace is full for the first time
         !      newsize = redsize
         !   else
         !      newsize = redsize+1
         !   endif
         !   allocate(tempA(newsize,newsize),tempS(newsize,newsize))
         !   tempA(1:newsize,1:newsize) = Ared(1:newsize,1:newsize)
         !   tempS(1:newsize,1:newsize) = Sred(1:newsize,1:newsize)
   !!Diagonalize Sred to test for linear dependencies:
   !if (.true.) then
   !   call dsyevx_interface(tempS,newsize,.false.,arh%lupri,'Metric')
   !   tempS(1:newsize,1:newsize) = Sred(1:newsize,1:newsize)
   !endif
         !   call arh_crop_optimal_x(tempA,tempS,newsize,vectorsubspace,xF,sigmaF)
         !   deallocate(tempA,tempS)
         !endif

         call mem_alloc(tempvec,redsize)
	 call mem_alloc(tempmat,redsize,redsize)
         tempvec(1:redsize) = arh%Gred(1:redsize)
         arh%Gred = 0.0E0_realk !For debug only
         arh%Gred(1:redsize-1) = tempvec(2:redsize)

         tempmat(1:redsize,1:redsize) = arh%Ared(1:redsize,1:redsize)
         arh%Ared = 0.0E0_realk !For debug only
         arh%Ared(1:redsize-1,1:redsize-1) = tempmat(2:redsize,2:redsize) 

         tempmat(1:redsize,1:redsize) = arh%Sred(1:redsize,1:redsize)
         arh%Sred = 0.0E0_realk !For debug only
         arh%Sred(1:redsize-1,1:redsize-1) = tempmat(2:redsize,2:redsize)

         tempmat(1:redsize,1:redsize) = arh%CROPmat(1:redsize,1:redsize)
         arh%CROPmat = 0.0E0_realk !For debug only
         arh%CROPmat(1:redsize-1,1:redsize-1) = tempmat(2:redsize,2:redsize)

         if (arh%info_crop) then
            printsize = printsize+1
            write (arh%lupri,*) 'Residual subspace after truncate:'
            call LS_OUTPUT(arh%CROPmat, 1, qsize, 1, qsize, qsize, qsize, 1, arh%lupri)
            write (arh%lupri,*) 'Ared after truncate:'
            call LS_OUTPUT(arh%Ared, 1, printsize, 1, printsize, Asize, Asize, 1, arh%lupri)
            write (arh%lupri,*) 'Sred after truncate:'
            call LS_OUTPUT(arh%Sred, 1, printsize, 1, printsize, Asize, Asize, 1, arh%lupri)
            write (arh%lupri,*) 'Gred after truncate:'
            call LS_OUTPUT(arh%Gred, 1, printsize, 1, 1, Asize, 1, 1, arh%lupri)
         endif

         call mem_dealloc(tempvec)
	 call mem_dealloc(tempmat)
         !call ls_flshfo(arh%lupri)
      endif
      arh%Gred(redsize) = mat_dotproduct(G,x)
   else
      if (arh%info_crop) then
         write (arh%lupri,*) 'Residual subspace before augmentation:'
         call LS_OUTPUT(arh%CROPmat, 1, n, 1, n, 200, 200, 1, arh%lupri)
         write (arh%lupri,*) 'Ared before augmentation:'
         call LS_OUTPUT(arh%Ared, 1, n, 1, n, 200, 200, 1, arh%lupri)
         write (arh%lupri,*) 'Sred before augmentation:'
         call LS_OUTPUT(arh%Sred, 1, n, 1, n, 200, 200, 1, arh%lupri)
         write (arh%lupri,*) 'Gred before augmentation:'
         call LS_OUTPUT(arh%Gred, 1, n, 1, 1, 200, 1, 1, arh%lupri)
      endif

      redsize = n
      arh%Gred(redsize+1) = mat_dotproduct(G,x)
   endif

   !Set up reduced Hessian Ared, overlap Sred:

   !Diagonal elements - not really necessary:
   if (.not. arh%cfg_arh_truncate) then
      arh%Ared(redsize+1,redsize+1) = mat_dotproduct(x,sigma)
      arh%Sred(redsize+1,redsize+1) = mat_dotproduct(x,x)
   endif
   call mat_add(-1E0_realk,G,-1E0_realk, sigma, scr) !scr = n+1'th optimal residual, r_n+1 par 
   call mat_daxpy(mu,x,scr)
   !Preconditioning:
   call arh_precond(arh,decomp,scr,symm,mu,resP)
   if (.not. arh%cfg_arh_truncate) then
      arh%CROPmat(redsize+1,redsize+1) = mat_dotproduct(scr,resP)
      rewind(lusigma) ; rewind(lub)
   endif
   if (.not. arh%cfg_arh_truncate) then
      call mat_init(xvec,rowdim,coldim)
      call mat_init(scr2,rowdim,coldim)
   endif
   do i = 1, redsize  
      if (arh%cfg_arh_truncate) then
         call get_from_modFIFO(vectorsubspace, i, xpointer, sigmapointer)
         arh%Ared(redsize,i) = mat_dotproduct(x,sigmapointer)
         call mat_add(-1E0_realk,G,-1E0_realk, sigmapointer, scr) !scr = i'th optimal residual, r_i par
         call mat_daxpy(mu,xpointer,scr) 
         arh%CROPmat(redsize,i) = mat_dotproduct(resP,scr)
         arh%CROPmat(i,redsize) = arh%CROPmat(redsize,i)
      else
         call mat_read_from_disk(lusigma,scr2,OnMaster)
         call mat_read_from_disk(lub,xvec,OnMaster)
         arh%Ared(redsize+1,i) = mat_dotproduct(x,scr2)
         call mat_add(-1E0_realk,G,-1E0_realk, scr2, scr) !scr = i'th optimal residual, r_i par
         call mat_daxpy(mu,xvec,scr) 
         arh%CROPmat(redsize+1,i) = mat_dotproduct(resP,scr)
         arh%CROPmat(i,redsize+1) = arh%CROPmat(redsize+1,i)
      endif
   enddo 
   if (.not. arh%cfg_arh_truncate) then
      call mat_free(xvec)
      call mat_free(scr2)
   endif
   !Explicitly calculate upper half: 
   if (.not. arh%cfg_arh_truncate) then
      rewind(lub)
   endif
   do i = 1, redsize
      if (arh%cfg_arh_truncate) then
         call get_from_modFIFO(vectorsubspace, i, xpointer, sigmapointer)
         arh%Ared(i,redsize) = mat_dotproduct(xpointer,sigma)
         arh%Sred(i,redsize) = mat_dotproduct(xpointer,x)
         arh%Sred(redsize,i) = arh%Sred(i,redsize)
      else
         call mat_read_from_disk(lub,scr,OnMaster)
         arh%Ared(i,redsize+1) = mat_dotproduct(scr,sigma)
         arh%Sred(i,redsize+1) = mat_dotproduct(scr,x)
         arh%Sred(redsize+1,i) = arh%Sred(i,redsize+1)
      endif
   enddo

   if (arh%cfg_arh_truncate) then
      if (arh%info_crop) then
         write (arh%lupri,*) 'Residual subspace after augmentation:'
         call LS_OUTPUT(arh%CROPmat, 1, redsize, 1, redsize, qsize, qsize, 1, arh%lupri)
         write (arh%lupri,*) 'Ared after augmentation:'
         call LS_OUTPUT(arh%Ared, 1, redsize, 1, redsize, Asize, Asize, 1, arh%lupri)
         write (arh%lupri,*) 'Sred after augmentation:'
         call LS_OUTPUT(arh%Sred, 1, redsize, 1, redsize, Asize, Asize, 1, arh%lupri)
         write (arh%lupri,*) 'Gred after augmentation:'
         call LS_OUTPUT(arh%Gred, 1, redsize, 1, 1, Asize, 1, 1, arh%lupri)
      endif
   else
      if (arh%info_crop) then
         write (arh%lupri,*) 'Residual subspace after augmentation:'
         call LS_OUTPUT(arh%CROPmat, 1, n+1, 1, n+1, 200, 200, 1, arh%lupri)
         write (arh%lupri,*) 'Ared after augmentation:'
         call LS_OUTPUT(arh%Ared, 1, n+1, 1, n+1, 200, 200, 1, arh%lupri)
         write (arh%lupri,*) 'Sred after augmentation:'
         call LS_OUTPUT(arh%Sred, 1, n+1, 1, n+1, 200, 200, 1, arh%lupri)
         write (arh%lupri,*) 'Gred after augmentation:'
         call LS_OUTPUT(arh%Gred, 1, n+1, 1, 1, 200, 1, 1, arh%lupri)
      endif
   endif

   if (arh%debug_diag_redspace .or. arh%debug_dd) then
      if (arh%cfg_arh_truncate) then
         nsize = redsize
      else
         nsize = n+1
      endif
      call mem_alloc(tempA,nsize,nsize)
      call mem_alloc(tempS,nsize,nsize)
      call mem_alloc(tempmat,nsize,nsize)
      call mem_alloc(tempvec,nsize)
      tempA(1:nsize,1:nsize) = arh%Ared(1:nsize,1:nsize)
      tempS(1:nsize,1:nsize) = arh%Sred(1:nsize,1:nsize)
      call rgg_interface(tempA,tempS,nsize,arh%lupri,tempvec,tempmat)
      call find_min_eival(nsize,tempvec,min_eival,min_eival_index)
      write (arh%lupri,'("Lowest eigenvalue in reduced space: ",F12.6, "      ")') min_eival
      arh%final_redspace_eival = min_eival
      call mem_dealloc(tempA)
      call mem_dealloc(tempS)
   endif

   if (arh%cfg_arh_newdamp .and. n > qsize-1) then
      !We only add the extra dimension after truncation has begun
      call arh_crop_extra_dim(arh,xF,sigmaF,vectorsubspace,G)
   endif

   !Diagonalize Sred to test for linear dependencies:
   if (.false.) then
      call mem_alloc(tempmat,redsize,redsize)
      tempmat(1:redsize,1:redsize) = arh%Sred(1:redsize,1:redsize)
      call dsyevx_interface(tempmat,redsize,.false.,arh%lupri,'Metric')
      call mem_dealloc(tempmat)
   endif

   !call mat_free(res)
   !call mat_free(scr)

   !if (cfg_arh_newdamp) then
   !   call mat_free(xF)
   !   call mat_free(sigmaF)
   !endif
   end subroutine arh_crop_setup_redsp

   !> \brief Construct x and sigma corresponding to the lowest eigenvalue in reduced space.
   !> \author S. Host
   !> \date 2007
   !> 
   !> When we truncate the number of vectors we keep, the levelshift becomes
   !> unhappy! Here we try to save the vector corresponding to the lowest
   !> eigenvalue instead of the standard vector. This will hopefully force the
   !> lowest eigenvalue in the reduced space to be similar to the one obtained
   !> when the full space is kept. \n
   !> 
   !> The 'standard' reduced space has been augmented. \n
   !> 1) Diagonalize reduced Hessian                   \n
   !> 2) Construct x and sigma corresponding to lowest eigenvalue in full space \n
   !>    as      xF = sum_i c_i x_i \n
   !>    and sigmaF = sum_i c_i sigma_i where c_i are elements of eigenvector
   !>    corresponding to lowest eigenvalue, and x_i and sigma_i are the stored
   !>    vectors. 
   !>  - NOTE THAT THIS CODE HAS NOT BEEN FULLY TESTED!
   !>
   subroutine arh_crop_optimal_x(lupri,A,S,reddim,vectorsubspace,xF,sigmaF)
   implicit none
         !> Logical unit number for output file
         integer, intent(in)         :: lupri
         !> Size of reduced space
         integer, intent(in)         :: reddim
         !> Reduced Hessian
         real(realk)                 :: A(reddim,reddim)
         !> Overlap of trial vectors
         real(realk)                 :: S(reddim,reddim)
         !> Subspace of previous trial and sigma vectors
         TYPE(modFIFO),intent(inout) :: vectorsubspace
         !> The trial vector corresponding to the lowest eigenvalue
         TYPE(matrix),intent(inout)  :: xF
         !> The sigma vector corresponding to the lowest eigenvalue
         TYPE(matrix),intent(inout)  :: sigmaF
         integer                     :: i, ndim, min_index(1), nvecs
         type(matrix), pointer       :: xpointer, sigmapointer
         real(realk), pointer        :: eival(:), eivec(:,:)
         real(realk)                 :: norm

   call mem_alloc(eival,reddim)
   call mem_alloc(eivec,reddim,reddim)

   nvecs = vectorsubspace%offset

   !1) Diagonalize reduced Hessian:
   call rgg_interface(A,S,reddim,lupri,eival,eivec)
   !write (arh%lupri,*) 'TEST!!! All eivecs:'
   !call LS_OUTPUT(eivec, 1, reddim, 1, reddim, reddim, reddim, 1, arh%lupri)

   !Locate lowest eigenvalue:
   min_index = minloc(eival)
   !write(arh%lupri,*) 'Lowest eigenval is number', min_index(1)

   !write (arh%lupri,*) 'TEST!!! Lowest eivec:'
   !call LS_OUTPUT(eivec, 1, reddim, min_index(1), min_index(1), reddim, reddim, 1, arh%lupri)

   !3) Construct x and sigma corresponding to lowest eigenvalue in full space
   !   as      xF = sum_i c_i x_i
   !   and sigmaF = sum_i c_i sigma_i where c_i are elements of eigenvector 
   !   corresponding to lowest eigenvalue, and x_i and sigma_i are the stored
   !   vectors.

   !Include the old xF and sigmaF in the new one...
   if (nvecs < reddim) then
      call mat_scal(eivec(reddim,min_index(1)),xF)
      call mat_scal(eivec(reddim,min_index(1)),sigmaF)
   else
      call mat_zero(xF) ; call mat_zero(sigmaF)
   endif
   do i = 1, nvecs
      call get_from_modFIFO(vectorsubspace, i, xpointer, sigmapointer)
      call mat_daxpy(eivec(i,min_index(1)),xpointer,xF)
      call mat_daxpy(eivec(i,min_index(1)),sigmapointer,sigmaF)         
   enddo
   !Normalize xF and sigma F:
   norm = sqrt(mat_sqnorm2(xF))
   call mat_scal(1.0E0_realk/norm,xF)
   call mat_scal(1.0E0_realk/norm,sigmaF)

   call mem_dealloc(eival)
   call mem_dealloc(eivec)
   end subroutine arh_crop_optimal_x

   !> \brief Construct x and sigma corresponding to the lowest eigenvalue in reduced space.
   !> \author S. Host
   !> \date 2007
   !> 
   !> When we truncate the number of vectors we keep, the levelshift becomes
   !> unhappy! Here we try to save the vector corresponding to the lowest
   !> eigenvalue instead of the standard vector. This will hopefully force the
   !> lowest eigenvalue in the reduced space to be similar to the one obtained
   !> when the full space is kept. \n
   !> The 'standard' reduced space has been augmented. Now
   !> construct the extra dimension using the incoming xF and sigmaF.
   !>
   subroutine arh_crop_extra_dim(arh,xF,sigmaF,vectorsubspace,G)
   implicit none
         !> Contains solver info (ARH/TrFD)
         type(solverItem)            :: arh
         !> The trial vector corresponding to the lowest eigenvalue
         TYPE(matrix),intent(in)     :: xF
         !> The sigma vector corresponding to the lowest eigenvalue
         TYPE(matrix),intent(in)     :: sigmaF
         !> Subspace of previous trial and sigma vectors
         TYPE(modFIFO),intent(inout) :: vectorsubspace
         !> Gradient
         TYPE(matrix), intent(in)    :: G
         integer                     :: i, Asize, nvecs
         type(matrix), pointer       :: xpointer, sigmapointer
         real(realk), pointer        :: tempmat(:,:)
        interface
           subroutine dsyevx_interface(A,ndim,print_eivecs,lupri,string)
             use precision
             implicit none
             character(6), intent(in), optional :: string
             integer, intent(in)      :: ndim, lupri
             real(realk), intent(in)  :: A(ndim,ndim)
             logical, intent(in)      :: print_eivecs 
           end subroutine dsyevx_interface
        end interface

   Asize = vectorsubspace%queuesize
   nvecs = vectorsubspace%offset
   !write(arh%lupri,*) 'Asize, nvecs:', Asize, nvecs

   if (arh%info_crop) then
      write (arh%lupri,*) 'Ared subspace initial, extra_dim:'
      call LS_OUTPUT(arh%Ared, 1, Asize-1, 1, Asize-1, Asize, Asize, 1, arh%lupri)
      write (arh%lupri,*) 'Sred subspace initial, extra_dim:'
      call LS_OUTPUT(arh%Sred, 1, Asize-1, 1, Asize-1, Asize, Asize, 1, arh%lupri)
      write (arh%lupri,*) 'Gred subspace initial, extra_dim:'
      call LS_OUTPUT(arh%Gred, 1, Asize-1, 1, 1, Asize, 1, 1, arh%lupri)
      !call ls_flshfo(arh%lupri)
   endif

   arh%Gred(Asize) = mat_dotproduct(G,xF)

   !Set up reduced Hessian Ared, overlap Sred:

   !Diagonal elements:
   arh%Ared(Asize,Asize) = mat_dotproduct(xF,sigmaF)
   arh%Sred(Asize,Asize) = mat_dotproduct(xF,xF)
   do i = 1, Asize-1
      call get_from_modFIFO(vectorsubspace, i, xpointer, sigmapointer)
      arh%Ared(Asize,i) = mat_dotproduct(xF,sigmapointer)
   enddo 
   !Explicitly calculate upper half: 
   do i = 1, Asize-1
      call get_from_modFIFO(vectorsubspace, i, xpointer, sigmapointer)
      arh%Ared(i,Asize) = mat_dotproduct(xpointer,sigmaF)
      arh%Sred(i,Asize) = mat_dotproduct(xpointer,xF)
      arh%Sred(Asize,i) = arh%Sred(i,Asize)
   enddo

   if (arh%info_crop) then
      write (arh%lupri,*) 'Ared after augmentation, extra_dim:'
      call LS_OUTPUT(arh%Ared, 1, Asize, 1, Asize, Asize, Asize, 1, arh%lupri)
      write (arh%lupri,*) 'Sred after augmentation, extra_dim:'
      call LS_OUTPUT(arh%Sred, 1, Asize, 1, Asize, Asize, Asize, 1, arh%lupri)
      write (arh%lupri,*) 'Gred after augmentation, extra_dim:'
      call LS_OUTPUT(arh%Gred, 1, Asize, 1, 1, Asize, 1, 1, arh%lupri)
   endif

   !Diagonalize Sred to test for linear dependencies:
   if (.false.) then
      call mem_alloc(tempmat,Asize,Asize)
      tempmat(1:Asize,1:Asize) = arh%Sred(1:Asize,1:Asize)
      call dsyevx_interface(tempmat,Asize,.false.,arh%lupri,'Metric')
      call mem_dealloc(tempmat)
   endif
   end subroutine arh_crop_extra_dim

   !> \brief Construct denominator of ratio used in trust-region update.
   !> \author S. Host
   !> \date 2010
   !> 
   !> See Molecular Electronic Structure Theory p. 615
   !> The update of trust-radius is called from get_fock in linscf.f90
   !>
   subroutine arh_get_TR_denom(arh, G, x, lintra, unres, mu)
   implicit none
      !> Contains solver info (ARH/TrFD)
      type(SolverItem),intent(inout) :: arh
      !> Gradient
      type(matrix), intent(in)       :: G
      !> Final solution vector
      type(matrix), intent(in)       :: x 
      !> Linear transformation of final solution vector
      type(matrix), intent(in)       :: lintra
      !> True if unrestricted calculation
      logical, intent(in)            :: unres
      !> level shift
      real(realk), intent(in)        :: mu

   call mat_max_elm(x, arh%maxelm)
   arh%xnorm = SQRT(mat_sqnorm2(x))

   !The predicted SCF energy change is
   !Q = E(SCF)_n + (E1_red)T * c + 1/2*cT * (E2_red) * c
   !and the ratio between actual and predicted energy change is
   !    E(SCF)_n+1 - E(SCF)_n
   !r = _____________________
   !        Q - E(SCF)_n
 
   !Calculate ratio denominator
   !Q - E(SCF)_n = (E1_red)T * c + 1/2*cT * (E2_red) * c  
   
   arh%denom = epred(x, lintra, G, mu)
   if (.not. unres) then !For closed shell, the matrices have been scaled by
                         !1/2, so we have to multiply denominator by two.
      arh%denom = 2.0E0_realk*arh%denom
   endif

   if (arh%denom > 0.0E0_realk) then
      write(arh%lupri,*) 'WARNING!!!! Predicted energy is positive!'
   endif

   end subroutine arh_get_TR_denom

   !> \brief Get predicted energy Epred = 1/2 x*A*x - x*b
   !> \author S. Host
   !> \date 2005
   !>
   !> Original formula is Epred = 1/2 x*A*x + x*b
   !> but because of inconsistencies with signs, the 
   !> following should be correct: -1/2 x*A*x - x*b - 1/2 x*x
   !> This is because linear transformation Ax has been carried out
   !> with solution x, which has afterwards been scaled to -x.
   !>
   real(realk) function epred(x, Ax, b, mu)
   implicit none
        type(matrix), intent(in)   :: x, Ax, b
        !> level shift
        real(realk), intent(in)    :: mu

   !Stinne, Brano: New (hopefully correct) Epredicted, October 2010
   epred = - mat_dotproduct(x,b) - 0.5E0_realk*mat_dotproduct(x,Ax) 
   end function epred

end module arhDensity

