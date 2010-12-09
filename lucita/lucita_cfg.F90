!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_cfg

! stefan: please note that if you modify this module (e.g. add new variables)
!         do not forget to update sync_coworkers_cfg (in module sync_coworkers)

  implicit none

  public lucita_orb_spaces_setup_delete_total

  save

! parameters

  integer, parameter, public :: max_number_of_gas_spaces = 6
  integer, parameter, public :: max_number_of_ptg_irreps = 8

! character block

  character (len = 72), public :: lucita_cfg_run_title = 'no title specified'
  character (len = 72), public :: lucita_cfg_ini_wavef = 'none'
  character (len = 72), public :: lucita_cfg_ci_type   = 'none'
  character (len = 72), public :: lucita_cfg_calculation_size = 'NORMAL'

! logical block
  logical, public :: lucita_cfg_inactive_shell_set   =  .false.
  logical, public :: lucita_cfg_nr_gas_space_set     =  .false.
  logical, public :: lucita_cfg_minmax_occ_gas_set   =  .false.
  logical, public :: lucita_cfg_ras1_set             =  .false.
  logical, public :: lucita_cfg_ras2_set             =  .false.
  logical, public :: lucita_cfg_ras3_set             =  .false.

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
  integer, public :: lucita_cfg_nr_calc_sequences    =  0
  integer, public :: lucita_cfg_max_holes_ras1       =  0
  integer, public :: lucita_cfg_max_e_ras3           =  0
  integer, public :: lucita_cfg_density_calc_lvl     =  0
  integer, public :: lucita_cfg_restart_ci           =  0
  integer, public :: lucita_cfg_max_dav_subspace_dim =  0
  integer, public :: lucita_cfg_max_nr_dav_ci_iter   =  30
  integer, public :: lucita_cfg_max_batch_size       =  100000000
  integer, public :: lucita_cfg_init_wave_f_type     =  0
  integer, public :: lucipar_cfg_ttss_dist_strategy  =  2
  integer, public :: lucipar_cfg_mem_reduction_multp =  3

  integer, allocatable, public :: nish_lucita(:)
  integer, allocatable, public :: nash_lucita(:)
  integer, allocatable, public :: nas1_lucita(:)
  integer, allocatable, public :: nas2_lucita(:)
  integer, allocatable, public :: nas3_lucita(:)
  integer, allocatable, public :: nocc_lucita(:)
  integer, allocatable, public :: nssh_lucita(:)
  integer, allocatable, public :: ngsh_lucita(:,:)
  integer, allocatable, public :: ngso_lucita(:,:)

  private

contains

 subroutine lucita_orb_spaces_setup_delete_total
!*******************************************************************************
!
!    purpose: deallocate all CI orbital space arrays for LUCITA.
!
!*******************************************************************************
!-------------------------------------------------------------------------------

      if(allocated(nish_lucita))deallocate(nish_lucita)
      if(allocated(nash_lucita))deallocate(nash_lucita)
      if(allocated(nas1_lucita))deallocate(nas1_lucita)
      if(allocated(nas2_lucita))deallocate(nas2_lucita)
      if(allocated(nas3_lucita))deallocate(nas3_lucita)
      if(allocated(nocc_lucita))deallocate(nocc_lucita)
      if(allocated(nssh_lucita))deallocate(nssh_lucita)
      if(allocated(ngsh_lucita))deallocate(ngsh_lucita)
      if(allocated(ngso_lucita))deallocate(ngso_lucita)

  end subroutine lucita_orb_spaces_setup_delete_total
!*******************************************************************************

end module
