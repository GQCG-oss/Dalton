!dalton_copyright_start
!
!
!dalton_copyright_end

#ifdef VAR_MPI
module sync_coworkers

! stefan: - this module provides all necessary functionality
!           to synchronize the co-workers in parallel mcscf/ci 
!           calculations.
!
!           written by sknecht for DALTON, december 2010.
  use dalton_mpi
  use lucita_cfg
  use lucita_mcscf_ci_cfg
  use vector_xc_file_type
  use file_type_module, only : file_type, file_info
#ifdef USE_MPI_MOD_F90
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif

#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

  public sync_coworkers_cfg
  public set_sync_default

  private

  save

  integer(kind=MPI_INTEGER_KIND), private   :: istat(MPI_STATUS_SIZE)
  integer(kind=MPI_INTEGER_KIND), private   :: ierr
  integer                       , parameter :: my_comm_world = MPI_COMM_WORLD

  integer, public, parameter :: dim_sync_ctrl_array                    = 6

  logical, public            :: sync_ctrl_array(1:dim_sync_ctrl_array) = .false.

#ifdef PRG_DIRAC
#include "dgroup.h"
#endif

contains 
 
  subroutine sync_coworkers_cfg(xarray1,xarray2)
!******************************************************************************
!
!    purpose:  synchronize (if necessary) co-workers for parallel CI/MCSCF
!
!*******************************************************************************
    real(8), optional, intent(inout) :: xarray1(*)
    real(8), optional, intent(inout) :: xarray2(*)
!-------------------------------------------------------------------------------
    integer                          :: select_sync
!-------------------------------------------------------------------------------

      select_sync = 0     

      do ! loop over synchronization processes

        select_sync = select_sync + 1

        if(select_sync > dim_sync_ctrl_array) exit

        select case(sync_ctrl_array(select_sync))
          case(.true.)
!           select synchronization module
#ifdef LUCI_DEBUG
            print *,' *** select_sync = ***',select_sync
#endif
            select case(select_sync)
              case(1) ! CI default settings and MCSCF/CI input variables
                call sync_coworkers_ci_cfg()
              case(2) ! MCSCF environment settings
                call sync_coworkers_mc_cfg()
              case(3) ! synchronize 1-el (i|j) and 2-el (ab|cd) integrals
                if(present(xarray1) .and. present(xarray2))then
                  call sync_coworkers_ij_abcd(xarray1,xarray2)
                else
                  call quit('*** sync_coworkers_cfg: sync of 1-/2-el integrals requires'//          &
                            ' the output arrays in the argument list.***')
                end if
              case(4) ! distribute previous solution vector (i.e., we restart a CI run)
              case(5) ! update information for vector exchange
                call sync_coworkers_mc_vector_xc()
              case(6) ! set current CI task
                call sync_coworkers_citask()
              case default ! no option found, quit.
                print *,' *** sync_coworkers_cfg: synchronization module does not exist, program exits. ***'
                call quit('*** sync_coworkers_cfg: unknown synchronization module.***')
            end select
          case default ! .false.
#ifdef LUCI_DEBUG
!           issue a warning
            print *,' *** sync_coworkers_cfg: synchronization requested but not needed, program continues. ***'
#endif
        end select

!       reset control array
        sync_ctrl_array(select_sync) = .false.

      end do

  end subroutine sync_coworkers_cfg

  subroutine sync_coworkers_ci_cfg()
!*******************************************************************************
!
!    purpose:  provide co-workers with basic common block/orbital space
!              knowledge in parallel CI/MCSCF runs.
!
!*******************************************************************************
!-------------------------------------------------------------------------------

!     part 1 - the dynamic variables
!     ------------------------------

!     logical block
      call dalton_mpi_bcast(lucita_cfg_analyze_cvec,        0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_natural_orb_occ_nr,  0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_transition_densm,    0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_initialize_cb,       0, my_comm_world)

!     real(8) block
      call dalton_mpi_bcast(lucita_cfg_accepted_truncation, 0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_convergence_c      , 0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_auxiliary_conv_c   , 0, my_comm_world)

!     integer block
      call dalton_mpi_bcast(lucita_cfg_csym,                0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_hcsym,               0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_icstate,             0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_spin1_2el_op,        0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_spin2_2el_op,        0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_ptg_symmetry,        0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_el_operator_level,   0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_nr_roots,            0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_global_print_lvl,    0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_local_print_lvl,     0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_density_calc_lvl,    0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_spindensity_calc_lvl,0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_is_spin_multiplett,  0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_max_dav_subspace_dim,0, my_comm_world)
      call dalton_mpi_bcast(lucita_cfg_max_nr_dav_ci_iter,  0, my_comm_world)

!     part 2 - the static variables (update only if necessary, i.e. SETCI has been called)
!     ------------------------------------------------------------------------------------

      if(lucita_cfg_initialize_cb)then

#ifdef PRG_DIRAC
!       transfer information from Dirac cb to LUCITA
        lucita_cfg_nr_ptg_irreps = nbsym
#endif

!       character
        call dalton_mpi_bcast(lucita_cfg_run_title,           0, my_comm_world)
        call dalton_mpi_bcast(lucita_cfg_ci_type,             0, my_comm_world)
!       real(8)
        call dalton_mpi_bcast(lucita_cfg_core_energy        , 0, my_comm_world)
!       logical
        call dalton_mpi_bcast(lucita_cfg_inactive_shell_set,  0, my_comm_world)
        call dalton_mpi_bcast(lucita_cfg_minmax_occ_gas_set,  0, my_comm_world)
        call dalton_mpi_bcast(lucita_cfg_timing_par,          0, my_comm_world)
!       integer
!       a. general definitions/settings
        call dalton_mpi_bcast(lucita_cfg_nr_active_e,         0, my_comm_world)
        call dalton_mpi_bcast(lucita_cfg_restart_ci,          0, my_comm_world)
        call dalton_mpi_bcast(lucita_cfg_max_batch_size,      0, my_comm_world)
        call dalton_mpi_bcast(lucipar_cfg_ttss_dist_strategy, 0, my_comm_world)
        call dalton_mpi_bcast(lucipar_cfg_mem_reduction_multp,0, my_comm_world)
!       call dalton_mpi_bcast(lucita_cfg_nr_calc_sequences,   0, my_comm_world)

!       b. symmetry, spin and orbital spaces
        call dalton_mpi_bcast(lucita_cfg_nr_ptg_irreps,       0, my_comm_world)
        call dalton_mpi_bcast(lucita_cfg_nr_gas_spaces,       0, my_comm_world)
        call dalton_mpi_bcast(lucita_cfg_init_wave_f_type,    0, my_comm_world)
        call dalton_mpi_bcast(nish_lucita,                    0, my_comm_world)
        call dalton_mpi_bcast(naos_lucita,                    0, my_comm_world)
        call dalton_mpi_bcast(nmos_lucita,                    0, my_comm_world)
        call dalton_mpi_bcast(nash_lucita,                    0, my_comm_world)
        call dalton_mpi_bcast(nfro_lucita,                    0, my_comm_world)
        call dalton_mpi_bcast(nocc_lucita,                    0, my_comm_world)

        call dalton_mpi_bcast(ngsh_lucita,                    0, my_comm_world)
        call dalton_mpi_bcast(ngso_lucita,                    0, my_comm_world)

        select case(lucita_cfg_init_wave_f_type)
          case(1) ! GAS settings
!           nothing particular here
          case(2) ! RAS settings
!           logical
!           call dalton_mpi_bcast(lucita_cfg_ras1_set,        0, my_comm_world)
!           call dalton_mpi_bcast(lucita_cfg_ras2_set,        0, my_comm_world)
!           call dalton_mpi_bcast(lucita_cfg_ras3_set,        0, my_comm_world)
!           integer 
!           call dalton_mpi_bcast(nas1_lucita,                0, my_comm_world)
!           call dalton_mpi_bcast(nas2_lucita,                0, my_comm_world)
!           call dalton_mpi_bcast(nas3_lucita,                0, my_comm_world)
            call dalton_mpi_bcast(lucita_cfg_min_e_ras1,      0, my_comm_world)
            call dalton_mpi_bcast(lucita_cfg_max_e_ras1,      0, my_comm_world)
            call dalton_mpi_bcast(lucita_cfg_min_e_ras3,      0, my_comm_world)
            call dalton_mpi_bcast(lucita_cfg_max_e_ras3,      0, my_comm_world)
        end select
      end if

  end subroutine sync_coworkers_ci_cfg
!******************************************************************************

  subroutine sync_coworkers_mc_cfg()
!*******************************************************************************
!
!    purpose:  provide co-workers with the basic settings from the mcscf
!              environment.
!
!*******************************************************************************
#ifdef MOD_SRDFT
      use lucita_mcscf_srdftci_cfg
#endif
!-------------------------------------------------------------------------------

!     logical
        
      call dalton_mpi_bcast(docisrdft_mc2lu,                  0,my_comm_world)
      srdft_ci_with_lucita = docisrdft_mc2lu
      call dalton_mpi_bcast(integrals_from_mcscf_env,         0,my_comm_world)
      call dalton_mpi_bcast(mcscf_ci_update_ijkl,             0,my_comm_world)
      call dalton_mpi_bcast(mcscf_orbital_trial_vector,       0,my_comm_world)
      call dalton_mpi_bcast(mcscf_ci_trial_vector,            0,my_comm_world)
      call dalton_mpi_bcast(io2io_vector_exchange_mc2lu_lu2mc,0,my_comm_world)
      call dalton_mpi_bcast(vector_update_mc2lu_lu2mc,        0,my_comm_world)

!     real(8)
!     call dalton_mpi_bcast(einact_mc2lu,                     0, my_comm_world) ! no need to sync - we sync ecore/ecore_orig in lucita internally

!     integer
      call dalton_mpi_bcast(len_cref_mc2lu,                   0, my_comm_world)
      call dalton_mpi_bcast(len_hc_mc2lu,                     0, my_comm_world)
      call dalton_mpi_bcast(len_resolution_mat_mc2lu,         0, my_comm_world)
      call dalton_mpi_bcast(len_int1_or_rho1_mc2lu,           0, my_comm_world)
      call dalton_mpi_bcast(len_int2_or_rho2_mc2lu,           0, my_comm_world)
      call dalton_mpi_bcast(vector_exchange_type1,            0, my_comm_world)
      call dalton_mpi_bcast(vector_exchange_type2,            0, my_comm_world)

      call dalton_mpi_bcast(file_info%current_file_nr_active1,0, my_comm_world)
      call dalton_mpi_bcast(file_info%current_file_nr_active2,0, my_comm_world)

  end subroutine sync_coworkers_mc_cfg
!******************************************************************************

  subroutine sync_coworkers_mc_vector_xc()
!*******************************************************************************
!
!    purpose:  provide co-workers with the basic information to initiate
!              the vector exchange process.
!
!*******************************************************************************

!     logical
      call dalton_mpi_bcast(exchange_f_info%exchange_file_io2io, 0,my_comm_world)

!     integer
      call dalton_mpi_bcast(exchange_f_info%present_sym_irrep,   0, my_comm_world)
      call dalton_mpi_bcast(exchange_f_info%push_pull_switch,    0, my_comm_world)
      call dalton_mpi_bcast(exchange_f_info%total_nr_vectors,    0, my_comm_world)

      call dalton_mpi_bcast(file_info%current_file_nr_active1,   0, my_comm_world)
      call dalton_mpi_bcast(file_info%current_file_nr_active2,   0, my_comm_world)

  end subroutine sync_coworkers_mc_vector_xc
!******************************************************************************

  subroutine sync_coworkers_ij_abcd(xarray1,xarray2)
!******************************************************************************
!
!    purpose:  provide co-workers with 1-/2-electron integrals in parallel 
!              CI/MCSCF runs.
!
!******************************************************************************
    use lucita_energy_types
    real(8), intent(inout) :: xarray1(*)
    real(8), intent(inout) :: xarray2(*)
!-------------------------------------------------------------------------------
!   common block from lucita (note: the one below does not have its own include file)
    integer    :: i12s,i34s,i1234s,nint1,nint2,nbint1,nbint2
    COMMON/CINTFO/i12s,i34s,i1234s,nint1,nint2,nbint1,nbint2
!#include "parluci.h"
!-------------------------------------------------------------------------------

!     total # 1-/2-electron integrals
      call dalton_mpi_bcast(nint1,      0, my_comm_world)
      call dalton_mpi_bcast(nint2,      0, my_comm_world)
!     core energy + inactive energy and core energy
      call dalton_mpi_bcast(ecore,      0, my_comm_world)
      call dalton_mpi_bcast(ecore_orig, 0, my_comm_world)

!     synchronize the 1-electron integrals (use the generic bcast because
!     the 1-/2-electron integral arrays have too large dimensions in the calling
!     lucita routines (for historic reasons/other purposes, i do not know. stefan dec 2010)

!     call mpi_bcast(xarray1, nint1, mpi_real8, 0, my_comm_world, ierr)
!     if(.not.mcscf_ci_trial_vector) call mpi_bcast(xarray2, nint2, mpi_real8, 0, my_comm_world, ierr)
!     31-jan-2020 hjaaj: try again with dalton_mpi_bcast ...
      call dalton_mpi_bcast(xarray1(1:nint1), 0, my_comm_world)
      if(.not.mcscf_ci_trial_vector) call dalton_mpi_bcast(xarray2(1:nint2), 0, my_comm_world)

!     set sync_ctrl option
      sync_ctrl_array(3) = .false.

  end subroutine sync_coworkers_ij_abcd
!******************************************************************************

  subroutine sync_coworkers_citask()
!*******************************************************************************
!
!    purpose:  provide co-workers with knowledge about the parallel CI run task.
!
!*******************************************************************************
!-------------------------------------------------------------------------------

!     character
      call dalton_mpi_bcast(lucita_ci_run_id, 0, my_comm_world)

  end subroutine sync_coworkers_citask
!******************************************************************************

  subroutine set_sync_default(set_value_sync_ctrl_array,sync_ctrl_array_pos)
!*******************************************************************************
!
!    purpose:  set/reset the co-workers synchronization control array.
!
!*******************************************************************************
    logical, intent(in) :: set_value_sync_ctrl_array
    integer, intent(in) :: sync_ctrl_array_pos
!-------------------------------------------------------------------------------
    logical             :: set_sync_default_val
!-------------------------------------------------------------------------------

      set_sync_default_val = set_value_sync_ctrl_array

      call dalton_mpi_bcast(set_sync_default_val, 0, my_comm_world)
     
!     insert value in ctrl-array
      sync_ctrl_array(sync_ctrl_array_pos) = set_sync_default_val

  end subroutine set_sync_default
!******************************************************************************
  
end module
#else 
subroutine sync_coworkers
! dummy routine for non-mpi compilation
end
#endif
