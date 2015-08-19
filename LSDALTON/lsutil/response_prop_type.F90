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

  public mcdinputitem_set_default_config, mcdinputitem, &
       & free_MCDinputitem, & 
       & ALPHAinputitem, ALPHAinputitem_set_default_config, &
       & BETAinputitem, BETAinputitem_set_default_config, &
       & GAMMAinputitem, GAMMAinputitem_set_default_config, &
       & TPAinputitem, TPAinputitem_set_default_config, &
       & DTPAinputitem, DTPAinputitem_set_default_config, &
       & ESGinputitem, ESGinputitem_set_default_config, &
       & ESDinputitem, ESDinputitem_set_default_config, &
       & RSPSOLVERinputitem, RSPSOLVERiputitem_set_default_config
  private

Contains
  !> \brief free memory used in the MCD input handling.
  !> \author T. Kjaergaard
  !> \date October 2010  
  !> \param MCDinput the MCDinputItem structure to be freed 
  subroutine free_MCDinputitem(MCDinput)
    implicit none
    type(MCDinputItem),intent(inout)       :: MCDinput
    
    IF(MCDinput%nXcoor .NE. 0)THEN
       call mem_dealloc(MCDinput%Xcoor)
       MCDinput%nXcoor = 0
    ENDIF
    
  end subroutine free_mcdinputitem
  
  !> \brief Sets default parameters for MCD calculation.
  !> \author T. Kjaergaard
  !> \date October 2010
  !> \param MCDinput the MCDinputItem structure to be set 
  subroutine mcdinputitem_set_default_config(MCDinput)
    implicit none
    type(MCDinputItem),intent(inout)       :: MCDinput

    MCDinput%nexci = 0
    MCDinput%nMCDexci = 0
    MCDinput%london = .TRUE.
    MCDinput%nolondon = .TRUE.
    MCDinput%simulate = .TRUE.
    MCDinput%lorentz = .TRUE.
    !gamma is not used unless specifeid by input
    !because it is 0.005E0_realk for lorentz
    !and 0.0070851079363103793E0_realk for gaussian
    MCDinput%useinputgamma=.FALSE.   
    MCDinput%specific_states_in_input = .FALSE.
    MCDinput%gamma = 0E0_realk
    MCDinput%nsteps = 5000
    MCDinput%nVecForPeak = 10
    MCDinput%doAterms = .TRUE.
    MCDinput%doBterms = .TRUE.
    MCDinput%dampedMCD = .TRUE.
    MCDinput%nXcoor = 0
    nullify(MCDinput%Xcoor)
    nullify(MCDinput%EXSTATES)
  end subroutine mcdinputitem_set_default_config

  !> \brief Sets default parameters for polarizability calculation.
  !> \author K. Kristensen
  !> \date August 2010  
  subroutine ALPHAinputitem_set_default_config(alphainput)
    implicit none
    type(alphainputItem),intent(inout)       :: alphainput

    alphainput%nfreq = 1
    alphainput%nimfreq = 1
    alphainput%real_frequencies_in_input=.false.
    alphainput%imag_frequencies_in_input=.false.
    nullify(alphainput%BFREQ)
    nullify(alphainput%IMBFREQ)

  end subroutine ALPHAinputitem_set_default_config



  !> \brief Sets default parameters for 1st hyperpolarizability calculation.
  !> \author K. Kristensen
  !> \date August 2010  
  subroutine BETAinputitem_set_default_config(betainput)
    implicit none
    type(betainputItem),intent(inout)       :: betainput

    betainput%nbfreq=1
    betainput%nimbfreq=1
    betainput%ncfreq=1
    betainput%nimcfreq=1
    betainput%real_bfrequencies_in_input=.false.
    betainput%imag_bfrequencies_in_input=.false.
    betainput%real_cfrequencies_in_input=.false.
    betainput%imag_cfrequencies_in_input=.false.
    nullify(betainput%BFREQ)
    nullify(betainput%IMBFREQ)
    nullify(betainput%CFREQ)
    nullify(betainput%IMCFREQ)

  end subroutine BETAinputitem_set_default_config


  !> \brief Sets default parameters for 2nd hyperpolarizability calculation.
  !> \author K. Kristensen
  !> \date August 2010  
  subroutine GAMMAinputitem_set_default_config(gammainput)
    implicit none
    type(gammainputItem),intent(inout)       :: gammainput

    gammainput%nbfreq=1
    gammainput%nimbfreq=1
    gammainput%ncfreq=1
    gammainput%nimcfreq=1
    gammainput%ndfreq=1
    gammainput%nimdfreq=1
    gammainput%real_bfrequencies_in_input=.false.
    gammainput%imag_bfrequencies_in_input=.false.
    gammainput%real_cfrequencies_in_input=.false.
    gammainput%imag_cfrequencies_in_input=.false.
    gammainput%real_dfrequencies_in_input=.false.
    gammainput%imag_dfrequencies_in_input=.false.
    nullify(gammainput%BFREQ)
    nullify(gammainput%IMBFREQ)
    nullify(gammainput%CFREQ)
    nullify(gammainput%IMCFREQ)
    nullify(gammainput%DFREQ)
    nullify(gammainput%IMDFREQ)

  end subroutine GAMMAinputitem_set_default_config



  !> \brief Sets default parameters for TPA calculation.
  !> \author K. Kristensen
  !> \date August 2010  
  subroutine TPAinputitem_set_default_config(tpainput)
    implicit none
    type(tpainputItem),intent(inout)       :: tpainput

    tpainput%tpa_nexci = 1
    tpainput%specific_states_in_input = .false.
    nullify(tpainput%ExStates)

  end subroutine TPAinputitem_set_default_config


  !> \brief Sets default parameters for damped TPA calculation.
  !> \author K. Kristensen
  !> \date August 2010  
  subroutine DTPAinputitem_set_default_config(dtpainput)
    implicit none
    type(dtpainputItem),intent(inout)       :: dtpainput


    dtpainput%nfreq=0
    dtpainput%gamma=0.005E0_realk
    dtpainput%gamma_specified=.false.
    nullify(dtpainput%FREQ)

  end subroutine DTPAinputitem_set_default_config



  !> \brief Sets default parameters for excited state gradient calculation.
  !> \author K. Kristensen
  !> \date August 2010  
  subroutine ESGinputitem_set_default_config(esginput)
    implicit none
    type(esginputItem),intent(inout)       :: esginput

    esginput%esg_nexci = 1
    esginput%specific_states_in_input = .false.
    nullify(esginput%ExStates)

  end subroutine ESGinputitem_set_default_config


  !> \brief Sets default parameters for excited state dipole calculation.
  !> \author K. Kristensen
  !> \date August 2010  
  subroutine ESDinputitem_set_default_config(esdinput)
    implicit none
    type(esdinputItem),intent(inout)       :: esdinput

    esdinput%esd_nexci = 1
    esdinput%specific_states_in_input = .false.
    nullify(esdinput%ExStates)

  end subroutine ESDinputitem_set_default_config


  !> \brief Sets default parameters for the solver.
  !> \author J. Kauczor
  !> \date August 2010
  subroutine RSPSOLVERiputitem_set_default_config(rspsolverinput)
    implicit none
    type(rspsolverinputItem),intent(inout)  :: rspsolverinput

    rspsolverinput%rsp_complex = .false.
    rspsolverinput%rsp_gamma=0E0_realk
    rspsolverinput%rsp_cpp = .true.
    rspsolverinput%rsp_cmplxnew = .false.
    rspsolverinput%rsp_stdnew = .false.
    rspsolverinput%rsp_convdyn = .false.   
    rspsolverinput%rsp_olsen = .false.
    rspsolverinput%rsp_single_norm = .false.
    rspsolverinput%rsp_thresh = 1E-4_realk
    rspsolverinput%rsp_maxit = 100
    rspsolverinput%rsp_maxred = 200
    rspsolverinput%rsp_quiet=.false.
    rspsolverinput%rsp_convdyn_type = 'TIGHT'
    rspsolverinput%rsp_conv_factor = 1.0E-3_realk
    rspsolverinput%rsp_dyn_thresh = 1E-4_realk
    ! KK, do MO preconditioning by default
    rspsolverinput%rsp_mo_precond = .true.
    rspsolverinput%rsp_mostart = .true.
    rspsolverinput%rsp_no_precond = .false.
    rspsolverinput%rsp_maxvec = 120
    rspsolverinput%rsp_maxgd  = 100
    rspsolverinput%rsp_nlintra = 0
    rspsolverinput%rsp_restart_exci = .false.
    rspsolverinput%rsp_restart_nexci = 1
    rspsolverinput%rsp_conv_thr = 1.0E-8_realk
    rspsolverinput%rsp_ovlmin = 1E-5_realk
    rspsolverinput%rsp_thr_lin_depend = 1E-20_realk
    rspsolverinput%rsp_thr_round = 1E-10_realk
    rspsolverinput%rsp_t1min = 1E-8_realk
    rspsolverinput%rsp_tolerance = 1.0E-7_realk 
    !rspsolverinput%no_of_startvectors
    rspsolverinput%rsp_startvectors = .false.
    rspsolverinput%rsp_no_of_startvectors = 0
    rspsolverinput%cfg_unres = .false.
    rspsolverinput%info_rsp = .false.
    rspsolverinput%info_rsp_redspace = .false.
    rspsolverinput%info_rsp_sparsity = .false.
    rspsolverinput%rsp_damp_2start = .false.
    rspsolverinput%degeneratestates = .FALSE.
    rspsolverinput%degenerateTHR = 1.d-5
    rspsolverinput%UseExcitationVecs = .FALSE.     
    rspsolverinput%rsp_eigenvecs = 0
    rspsolverinput%DoSVD = .FALSE.
  endsubroutine RSPSOLVERiputitem_set_default_config

end module response_wrapper_type_module
