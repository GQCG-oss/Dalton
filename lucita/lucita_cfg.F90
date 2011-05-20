!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_cfg

! stefan: please note that if you modify this module (e.g. add new variables)
!         do not forget to update sync_coworkers_cfg (in module sync_coworkers
!         in file par_sync_cw.F90)

  implicit none

  save

! parameters

  integer, parameter, public :: max_number_of_gas_spaces     =  6
  integer, parameter, public :: max_number_of_ptg_irreps     =  8
  integer, parameter, public :: lucita_cfg_nr_calc_sequences =  1 ! always one in this version, therefore a "parameter"

! character block

  character (len = 72), public :: lucita_cfg_run_title        = 'no title specified'
  character (len = 72), public :: lucita_cfg_ini_wavef        = 'none'
  character (len = 72), public :: lucita_cfg_ci_type          = 'none'
  character (len = 72), public :: lucita_cfg_calculation_size = 'NORMAL'

! logical block
  logical, public :: lucita_cfg_inactive_shell_set   =  .false.
  logical, public :: lucita_cfg_nr_gas_space_set     =  .false.
  logical, public :: lucita_cfg_minmax_occ_gas_set   =  .false.
  logical, public :: lucita_cfg_ras1_set             =  .false.
  logical, public :: lucita_cfg_ras2_set             =  .false.
  logical, public :: lucita_cfg_ras3_set             =  .false.
  logical, public :: lucita_cfg_analyze_cvec         =  .false.
  logical, public :: lucita_cfg_timing_par           =  .false.
  logical, public :: lucita_cfg_natural_orb_occ_nr   =  .false.


! double precision block

  real(8), public :: lucita_cfg_accepted_truncation = 1.0d-10

! integer block

  integer, public :: lucita_cfg_nr_roots             =  1
  integer, public :: lucita_cfg_ptg_symmetry         =  1
  integer, public :: lucita_cfg_nr_active_e          = -1
  integer, public :: lucita_cfg_is_spin_multiplett   = -1
  integer, public :: lucita_cfg_global_print_lvl     =  0
  integer, public :: lucita_cfg_local_print_lvl      =  0
  integer, public :: lucita_cfg_nr_gas_spaces        =  0
  integer, public :: lucita_cfg_nr_ptg_irreps        =  0
 !integer, public :: lucita_cfg_nr_calc_sequences    =  1 ! number of different GAS specifications
  integer, public :: lucita_cfg_max_holes_ras1       =  0
  integer, public :: lucita_cfg_max_e_ras3           =  0
  integer, public :: lucita_cfg_density_calc_lvl     =  0
  integer, public :: lucita_cfg_restart_ci           =  0
  integer, public :: lucita_cfg_max_dav_subspace_dim =  0
  integer, public :: lucita_cfg_max_nr_dav_ci_iter   =  30
  integer, public :: lucita_cfg_max_batch_size       =  100000000
  integer, public :: lucita_cfg_init_wave_f_type     =  0
  integer, public :: lucita_cfg_init_input_type      =  0
  integer, public :: lucipar_cfg_ttss_dist_strategy  =  2
  integer, public :: lucipar_cfg_mem_reduction_multp =  3

! input orbital spaces
  integer, public :: nish_lucita(max_number_of_ptg_irreps)
  integer, public :: nas1_lucita(max_number_of_ptg_irreps)
  integer, public :: nas2_lucita(max_number_of_ptg_irreps)
  integer, public :: nas3_lucita(max_number_of_ptg_irreps)
  integer, public :: ngsh_lucita(max_number_of_gas_spaces,max_number_of_ptg_irreps)
  integer, public :: ngso_lucita(max_number_of_gas_spaces,2)

! calculated orbital spaces based on input orbital spaces
  integer, public :: nfro_lucita(max_number_of_ptg_irreps) ! might become input in future
  integer, public :: nash_lucita(max_number_of_ptg_irreps)
  integer, public :: nocc_lucita(max_number_of_ptg_irreps)

end module
