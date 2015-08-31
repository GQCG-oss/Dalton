!> @file
!> Contains property wrappers to the response_driver_module in response_driver.f90

!> Calculates A and B terms of Magnetic Circular Dichroisme (MCD)
!> based on J. Chem. Theory Comput. 2009, 5, 1997-2020  but reformulated using JCP, 129, 214108 (2008).
!> \author Thomas Kjaergaard
!> \date 2010-03
module response_wrapper_type_module
  use precision

  type mcdinputitem
     !> Number of excited states for MCD calculation.
     Integer :: nexci
     !> Number of excited states for MCD calculation for which A and B terms
     !> should be calculated 
     Integer :: nMCDexci
     logical :: specific_states_in_input
     logical :: degeneratestates
     logical :: london
     logical :: nolondon
     logical :: simulate
     logical :: lorentz
     logical :: useinputgamma
     real(realk) :: gamma
     Integer :: nsteps
     Integer :: nVecForPeak
     logical :: doAterms
     logical :: doBterms
     logical :: dampedMCD
     Integer :: nXcoor
     !> Which specific excited states are requested.
     integer, pointer :: ExStates(:)
     real(realk),pointer :: Xcoor(:)
  end  type mcdinputitem

  !> Information for polarizability (alpha) calculation.
  type ALPHAinputitem
     !> Number of real and imaginary frequencies, respectively.
     integer :: nfreq, nimfreq
     !> Are any real frequencies defined in the input?
     logical :: real_frequencies_in_input
     !> Are any imaginary frequencies defined in the input?
     logical :: imag_frequencies_in_input
     real(realk), pointer :: BFREQ(:),IMBFREQ(:)
  end type ALPHAinputitem

  !> Information for 1st polarizability (beta) calculation.
  type BETAinputitem
     !> Number of real and imaginary frequencies for B and C operators.
     integer :: nbfreq, nimbfreq, ncfreq, nimcfreq
     !> Are any real frequencies for operator B defined in the input?
     logical :: real_bfrequencies_in_input
     !> Are any imag frequencies for operator B defined in the input?
     logical :: imag_bfrequencies_in_input
     !> Are any real frequencies for operator C defined in the input?
     logical :: real_cfrequencies_in_input
     !> Are any imag frequencies for operator C defined in the input?
     logical :: imag_cfrequencies_in_input
     real(realk), pointer :: BFREQ(:),IMBFREQ(:), CFREQ(:),IMCFREQ(:)
  end type BETAinputitem

  !> Information for 2nd polarizability (gamma) calculation.
  type GAMMAinputitem
     !> Number of real and imaginary frequencies for operator B,C, and D.
     integer :: nbfreq, nimbfreq, ncfreq, nimcfreq, ndfreq, nimdfreq
     !> Are any real frequencies for operator B defined in the input?
     logical :: real_bfrequencies_in_input
     !> Are any imag frequencies for operator B defined in the input?
     logical :: imag_bfrequencies_in_input
     !> Are any real frequencies for operator C defined in the input?
     logical :: real_cfrequencies_in_input
     !> Are any imag frequencies for operator C defined in the input?
     logical :: imag_cfrequencies_in_input
     !> Are any real frequencies for operator D defined in the input?
     logical :: real_dfrequencies_in_input
     !> Are any imag frequencies for operator D defined in the input?
     logical :: imag_dfrequencies_in_input
     real(realk), pointer :: BFREQ(:),IMBFREQ(:), &
          & CFREQ(:),IMCFREQ(:), DFREQ(:),IMDFREQ(:)
  end type GAMMAinputitem


  !> Information for two-photon absorption calculation.
  !> Note: The maximum excited state is determined by the .NEXCIT keyword 
  !> under **RESPONSE.
  !> Also note that we assume the two individual photons have the same
  !> frequency = 1/2 * excitation energy.
  type TPAinputitem
     !> Specific excited states defined in input
     logical :: specific_states_in_input
     !> Number of excited states for TPA calculation.
     !> Note that the position of the maximum excited states must be
     !> defined by the .NEXCIT ketword under **RESPONSE.
     integer :: tpa_nexci
     !> Which specific excited states are requested.
     integer, pointer :: ExStates(:)
  end type TPAinputitem



  !> Information for DAMPED two-photon absorption calculation.
  !> We assume the two individual photons have the same
  !> frequency = 1/2 * excitation energy.
  type DTPAinputitem
     !> Number of different frequencies.
     integer :: nfreq
     !> Damping parameter gamma
     real(realk) :: gamma
     !> Gamma specified in input?
     logical :: gamma_specified
     !> Which frequencies
     real(realk), pointer :: FREQ(:)
  end type DTPAinputitem


  !> Information for excited state gradient (ESG).
  !> Note: The maximum excited state is determined by the .NEXCIT keyword 
  !> under **RESPONSE.
  type ESGinputitem
     !> Specific excited states defined in input
     logical :: specific_states_in_input
     !> Number of excited states for ESG calculation.
     !> Note that the position of the maximum excited states must be
     !> defined by the .NEXCIT ketword under **RESPONSE.
     integer :: esg_nexci
     !> Which specific excited states are requested.
     integer, pointer :: ExStates(:)
  end type ESGinputitem


  !> Information for excited state dipole moment (ESD).
  !> Note: The maximum excited state is determined by the .NEXCIT keyword 
  !> under **RESPONSE.
  type ESDinputitem
     !> Specific excited states defined in input
     logical :: specific_states_in_input
     !> Number of excited states for ESD calculation.
     !> Note that the position of the maximum excited states must be
     !> defined by the .NEXCIT ketword under **RESPONSE.
     integer :: esd_nexci
     !> Which specific excited states are requested.
     integer, pointer :: ExStates(:)
  end type ESDinputitem

  !> General information for rsp solver.
  type RSPSOLVERinputitem
     !> damped response calculation requested, by default solver with symmetrized trial vectors is used
     logical :: rsp_complex
     !> damping parameter gamma
     real(realk) :: rsp_gamma
     !> standard response solver with symmetriazed trial vectors
     logical :: rsp_stdnew
     !> complex response solver with symmetrixed trial vectors
     logical :: rsp_cmplxnew
     !cpp solver
     logical :: rsp_cpp
     !> dynamic convergence requested
     logical :: rsp_convdyn
     !> Olsen algorithm requested for the eigenvalue equation
     logical :: rsp_olsen
     !> sloppy convergence threshold requested
     logical :: rsp_single_norm
     !> threshold for rsp calculation
     real(realk) :: rsp_thresh
     !> threshold for dynamic convergence
     real(realk):: rsp_dyn_thresh
     !> max number of iterations
     integer :: rsp_maxit
     !> max size of reduced space
     integer :: rsp_maxred
     !> no prints during calculation
     logical :: rsp_quiet
     !> dynamic convergence type
     character*5 :: rsp_convdyn_type
     !> dynamic convergence factor
     real(realk) :: rsp_conv_factor
     !> MO preconditioning
     logical :: rsp_mo_precond
     !> 
     logical :: rsp_mostart
     !> NO preconditioning
     logical :: rsp_no_precond
     !> maximum number of trial vectors
     integer :: rsp_maxvec
     !> maximum number of gradients
     integer :: rsp_maxgd
     !> number of linear transformation
     INTEGER     :: rsp_nlintra
     !> 
     LOGICAL     :: rsp_restart_exci
     !> 
     INTEGER     :: rsp_restart_nexci
     !> convergence threshold of residuals
     REAL(REALK) :: rsp_conv_thr
     !> (threshold found in orthonormalize)
     REAL(REALK) :: rsp_ovlmin
     !> (threshold found in orthonormalize)
     REAL(REALK) :: rsp_thr_lin_depend
     !> (threshold found in normalize)
     REAL(REALK) :: rsp_thr_round
     !> (threshold found in orthonormalize)
     REAL(REALK) :: rsp_t1min
     !> tolerance in NEX og start
     REAL(REALK) :: rsp_tolerance
     !>
     LOGICAL     :: rsp_startvectors
     !>
     INTEGER     :: rsp_no_of_startvectors
     !> unrestricted calculation?
     LOGICAL     :: cfg_unres
     !> print info about the soulution of the response calc
     LOGICAL     :: info_rsp
     !> print info about the soulution of the response calc reduced space
     LOGICAL     :: info_rsp_redspace
     !> print info about the sparsity of the trial vectors
     LOGICAL     :: info_rsp_sparsity
     !>
     LOGICAL     :: rsp_damp_2start
     !> Degenerate yes if degenerate states are possible
     logical :: degeneratestates
     !> Threshold for when excited states are considered degenerate
     real(realk) :: degenerateTHR
     !> Should the excitation vectors be used for linear equations
     logical :: UseExcitationVecs     
     !>
     INTEGER     :: rsp_eigenvecs
     !> Do a SVD decomposition of the Residual to improve convergence
     !> by generating several excitation vectors.
     logical :: DoSVD
  end type RSPSOLVERinputitem

end module response_wrapper_type_module
