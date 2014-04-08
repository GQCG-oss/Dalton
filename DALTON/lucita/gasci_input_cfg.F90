!dalton_copyright_start
!
!
!dalton_copyright_end

module gasci_input_cfg

  implicit none

  save

! parameters

  integer, parameter, public :: gasci_input_max_number_of_gas_spaces = 6
  integer, parameter, public :: gasci_input_max_number_of_ptg_irreps = 8
! Mickael : parameters for CT control 
!  integer, parameter, public :: gasci_input_ct_parameter_1           = 10
!  integer, parameter, public :: gasci_input_ct_parameter_2           = 30
!  integer, parameter, public :: gasci_input_max_number_of_occlses    = 5000 
   
! character block

  character (len = 72), public :: gasci_input_run_title        = 'no title specified'
  character (len = 72), public :: gasci_input_ini_wavef        = 'none'
  character (len = 72), public :: gasci_input_ci_type          = 'none'
  character (len = 72), public :: gasci_input_calculation_size = 'NORMAL'
! Mickael : character block for CT control 
!  character (len =  9), public :: CT(gasci_input_max_number_of_ct_types)     = 'none'
!  character (len =  9), public :: CTAB(gasci_input_max_number_of_ctab_types) = 'none'
!  character (len = 72), public :: GASAB                                      = 'none'

! logical block
  logical, public :: gasci_input_inactive_shell_set   =  .false.
  logical, public :: gasci_input_nr_gas_space_set     =  .false.
  logical, public :: gasci_input_minmax_occ_gas_set   =  .false.
  logical, public :: gasci_input_ras1_set             =  .false.
  logical, public :: gasci_input_ras2_set             =  .false.
  logical, public :: gasci_input_ras3_set             =  .false.
  logical, public :: gasci_input_analyze_cvec         =  .false.
  logical, public :: gasci_input_timing_par           =  .false.
  logical, public :: gasci_input_natural_orb_occ_nr   =  .true.
  logical, public :: gasci_input_skip_4index_trafo    =  .false.
  logical, public :: gasci_input_fci_dump             =  .false.
! Mickael : logical block for CT control 
!  logical, public :: gasci_input_charge_transfer      =  .false.
!  logical, public :: gasci_input_occlse_manual_control=  .false.
!  logical, public :: gasci_input_gas_conf_for_ct      =  .false.
  logical, public :: gasci_input_plus_combi           =  .false.

! double precision block

  real(8), public :: gasci_input_accepted_truncation = 1.0d-10
  real(8), public :: gasci_input_convergence_thr     = 2.0d-04
  real(8), public :: gasci_input_aux_convergence_thr = 2.0d-04
  real(8), public :: gasci_input_core_energy         = 0.0d0

! integer block

  integer, public :: gasci_input_nr_roots             =  1
  integer, public :: gasci_input_ptg_symmetry         =  1
  integer, public :: gasci_input_nr_active_e          = -1
  integer, public :: gasci_input_is_spin_multiplett   = -1
  integer, public :: gasci_input_global_print_lvl     =  0
  integer, public :: gasci_input_local_print_lvl      =  0
  integer, public :: gasci_input_nr_gas_spaces        =  0
  integer, public :: gasci_input_nr_ptg_irreps        =  0
  integer, public :: gasci_input_min_e_ras1           =  0
  integer, public :: gasci_input_max_e_ras1           =  0
  integer, public :: gasci_input_min_e_ras3           =  0
  integer, public :: gasci_input_max_e_ras3           =  0
  integer, public :: gasci_input_density_calc_lvl     =  0
  integer, public :: gasci_input_spindensity_calc_lvl =  0
  integer, public :: gasci_input_restart_ci           =  0
  integer, public :: gasci_input_max_dav_subspace_dim =  0
  integer, public :: gasci_input_max_nr_dav_ci_iter   =  30
  integer, public :: gasci_input_max_batch_size       =  64000000
  integer, public :: gasci_input_init_wave_f_type     =  0
  integer, public :: gasci_input_init_input_type      =  0
  integer, public :: gasci_input_ttss_dist_strategy   =  2
  integer, public :: gasci_input_mem_reduction_multp  =  3
! Mickael : integer block for CT control 
!  integer, public :: INCTAB                           =  0 
!  integer, public :: CTLIST()                         =  0 
!  integer, public :: CTLISTAB()                       =  0 
!  integer, public :: NCTLIST()                        =  0 
!  integer, public :: CTABN                            =  0 
!  integer, public :: CTN()                            =  0 
!  integer, public :: NCTLISTAB()                      =  0 
!  integer, public :: GASLIST()                        =  0 
!  integer, public :: NGASLIST                         =  0 
!  integer, public :: GASCTN                           =  0 
!  integer, public :: GAS_CT()                         =  0 
!  integer, public :: NSYS                             =  0 
!  integer, public :: NCT                              =  0 
!  integer, public :: NCLSOCC                          =  0 
!  integer, public :: NCLSOCC                          =  0 


! input orbital spaces
  integer, public :: nish_gasci_input(gasci_input_max_number_of_ptg_irreps)
  integer, public :: naos_gasci_input(gasci_input_max_number_of_ptg_irreps)
  integer, public :: nmos_gasci_input(gasci_input_max_number_of_ptg_irreps)
  integer, public :: nas1_gasci_input(gasci_input_max_number_of_ptg_irreps)
  integer, public :: nas2_gasci_input(gasci_input_max_number_of_ptg_irreps)
  integer, public :: nas3_gasci_input(gasci_input_max_number_of_ptg_irreps)
  integer, public :: ngsh_gasci_input(gasci_input_max_number_of_gas_spaces,gasci_input_max_number_of_ptg_irreps)
  integer, public :: ngso_gasci_input(gasci_input_max_number_of_gas_spaces,2)

! calculated orbital spaces based on input orbital spaces
  integer, public :: nfro_gasci_input(gasci_input_max_number_of_ptg_irreps) ! might become input in future
  integer, public :: nash_gasci_input(gasci_input_max_number_of_ptg_irreps)
  integer, public :: nocc_gasci_input(gasci_input_max_number_of_ptg_irreps)

end module
