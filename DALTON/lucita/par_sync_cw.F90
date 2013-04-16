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
  use lucita_mcscf_ci_task
#ifndef VAR_USE_MPIF
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

  integer, private :: istat(MPI_STATUS_SIZE)
  integer, private :: ierr

  integer, public, parameter :: dim_sync_ctrl_array                    = 6

  logical, public            :: sync_ctrl_array(1:dim_sync_ctrl_array) = .true.

#ifdef PRG_DIRAC
#include "dgroup.h"
#else
#include "inforb.h"
#endif

contains 
 
  subroutine sync_coworkers_cfg(select_sync,xarray1,xarray2)
!******************************************************************************
!
!    purpose:  synchronize (if necessary) co-workers for parallel CI/MCSCF
!
!*******************************************************************************
    integer, intent(in)              :: select_sync
    real(8), optional, intent(inout) :: xarray1(*)
    real(8), optional, intent(inout) :: xarray2(*)
!-------------------------------------------------------------------------------
    logical             :: need_sync
!-------------------------------------------------------------------------------
     
!     check for need of master-co-workers synchronization
      need_sync = sync_ctrl_array(select_sync)

      select case(need_sync)
        case(.true.)
!         select synchronization module
#ifdef LUCI_DEBUG
          print *,' *** select_sync = ***',select_sync
#endif
          select case(select_sync)
            case(1) ! CI default settings and input variables
              call sync_coworkers_ci_cfg()
            case(2) ! MCSCF default settings and input variables
            case(3) ! synchronize 1-el (i|j) and 2-el (ab|cd) integrals
              if(present(xarray1) .and. present(xarray2))then
                call sync_coworkers_ij_abcd(xarray1,xarray2)
              else
                call quit('*** sync_coworkers_cfg: sync of 1-/2-el integrals requires& 
 the output arrays in the argument list.***')
              end if
            case(4) ! distribute previous solution vector (i.e., we restart a CI run)
            case(5) ! synchronize with current expansion point (CEP) vector (MCSCF run)
            case(6) ! set current CI task
              call sync_coworkers_citask()
            case default ! no option found, quit.
              print *,' *** sync_coworkers_cfg: synchronization module does not exist, program exits. ***'
              call quit('*** sync_coworkers_cfg: unknown synchronization module.***')
          end select
        case default ! .false.
#ifdef LUCI_DEBUG
!         issue a warning
          print *,' *** sync_coworkers_cfg: synchronization requested but not needed, program continues. ***'
#endif
      end select
#ifdef LUCI_DEBUG
          print *,' *** sync_coworkerscfg: synchronization done. ***'
#endif
  end subroutine sync_coworkers_cfg

  subroutine sync_coworkers_ci_cfg()
!*******************************************************************************
!
!    purpose:  provide co-workers with basic common block/orbital space
!              knowledge in parallel CI/MCSCF runs.
!
!*******************************************************************************
!-------------------------------------------------------------------------------

!     character
      call dalton_mpi_bcast(lucita_cfg_run_title,           0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_ini_wavef,           0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_ci_type,             0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_calculation_size,    0, mpi_comm_world)
!     real(8)
      call dalton_mpi_bcast(lucita_cfg_accepted_truncation, 0, mpi_comm_world)
!     logical
      call dalton_mpi_bcast(lucita_cfg_inactive_shell_set,  0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_minmax_occ_gas_set,  0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_analyze_cvec,        0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_timing_par,          0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_natural_orb_occ_nr,  0, mpi_comm_world)
!     integer
#ifdef PRG_DIRAC
      lucita_cfg_nr_ptg_irreps = nbsym
#else
      lucita_cfg_nr_ptg_irreps = nsym
#endif
!     general definitions/settings
      call dalton_mpi_bcast(lucita_cfg_nr_roots,            0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_nr_active_e,         0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_global_print_lvl,    0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_local_print_lvl,     0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_density_calc_lvl,    0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_restart_ci,          0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_max_dav_subspace_dim,0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_max_nr_dav_ci_iter,  0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_max_batch_size,      0, mpi_comm_world)
      call dalton_mpi_bcast(lucipar_cfg_ttss_dist_strategy, 0, mpi_comm_world)
      call dalton_mpi_bcast(lucipar_cfg_mem_reduction_multp,0, mpi_comm_world)
!     call dalton_mpi_bcast(lucita_cfg_nr_calc_sequences,   0, mpi_comm_world)

!     symmetry, spin and orbital spaces
      call dalton_mpi_bcast(lucita_cfg_ptg_symmetry,        0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_is_spin_multiplett,  0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_nr_ptg_irreps,       0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_nr_gas_spaces,       0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_init_input_type,     0, mpi_comm_world)
      call dalton_mpi_bcast(lucita_cfg_init_wave_f_type,    0, mpi_comm_world)
      call dalton_mpi_bcast(nish_lucita,                    0, mpi_comm_world)

      select case(lucita_cfg_init_wave_f_type)
        case(1) ! GAS settings
          call dalton_mpi_bcast(ngsh_lucita,                0, mpi_comm_world)
          call dalton_mpi_bcast(ngso_lucita,                0, mpi_comm_world)
        case(2) ! RAS settings
!         logical
          call dalton_mpi_bcast(lucita_cfg_ras1_set,        0, mpi_comm_world)
          call dalton_mpi_bcast(lucita_cfg_ras2_set,        0, mpi_comm_world)
          call dalton_mpi_bcast(lucita_cfg_ras3_set,        0, mpi_comm_world)
!         integer 
          call dalton_mpi_bcast(nas1_lucita,                0, mpi_comm_world)
          call dalton_mpi_bcast(nas2_lucita,                0, mpi_comm_world)
          call dalton_mpi_bcast(nas3_lucita,                0, mpi_comm_world)
          call dalton_mpi_bcast(lucita_cfg_max_holes_ras1,  0, mpi_comm_world)
          call dalton_mpi_bcast(lucita_cfg_max_e_ras3,      0, mpi_comm_world)
      end select

  end subroutine sync_coworkers_ci_cfg
!******************************************************************************

  subroutine sync_coworkers_ij_abcd(xarray1,xarray2)
!******************************************************************************
!
!    purpose:  provide co-workers with 1-/2-electron integrals in parallel 
!              CI/MCSCF runs.
!
!******************************************************************************
    real(8), intent(inout) :: xarray1(*)
    real(8), intent(inout) :: xarray2(*)
!-------------------------------------------------------------------------------
!   common blocks from lucita (note: they do not have their own include file)
    real(8)    :: ecore,ecore_orig,ecore_h,ecore_hex
    COMMON/CECORE/ecore,ecore_orig,ecore_h,ecore_hex
    integer    :: i12s,i34s,i1234s,nint1,nint2,nbint1,nbint2
    COMMON/CINTFO/i12s,i34s,i1234s,nint1,nint2,nbint1,nbint2
!#include "parluci.h"
!-------------------------------------------------------------------------------

!     total # 1-/2-electron integrals
      call dalton_mpi_bcast(nint1,      0, mpi_comm_world)
      call dalton_mpi_bcast(nint2,      0, mpi_comm_world)
!     core energy + inactive energy and core energy
      call dalton_mpi_bcast(ecore,      0, mpi_comm_world)
      call dalton_mpi_bcast(ecore_orig, 0, mpi_comm_world)

!     synchronize the 1-electron integrals (use the generic bcast because
!     the 1-/2-electron integral arrays have too large dimensions in the calling
!     lucita routines (for historic reasons/other purposes, i do not know. stefan dec 2010)

      call mpi_bcast(xarray1, nint1, mpi_real8, 0, mpi_comm_world, ierr)
      call mpi_bcast(xarray2, nint2, mpi_real8, 0, mpi_comm_world, ierr)

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
      call dalton_mpi_bcast(lucita_citask_id, 0, mpi_comm_world)

!     more to come once we actually have an mcscf-lucita interface

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

      call dalton_mpi_bcast(set_sync_default_val, 0, mpi_comm_world)
     
!     insert value in ctrl-array
      sync_ctrl_array(sync_ctrl_array_pos) = set_sync_default_val

  end subroutine set_sync_default
!******************************************************************************
  
#else 
module dummy_sync_coworkers
#endif
end module
