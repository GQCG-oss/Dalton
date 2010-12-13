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
  use mpi
  use dalton_mpi
  use lucita_cfg

#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

  implicit none

  public sync_coworkers_cfg

  private

  save

  integer, private :: istat(MPI_STATUS_SIZE)
  integer, private :: ierr

  logical, public  :: sync_cw_ci_cfg = .true.
  logical, public  :: sync_cw_mc_cfg = .true.
  logical, public  :: sync_cw_abcd   = .true.
  logical, public  :: sync_cw_cvec   = .true.
  logical, public  :: sync_cw_refvec = .true.

contains 
 
  subroutine sync_coworkers_cfg
!******************************************************************************
!
!    purpose:  synchronize (if necessary) co-workers for parallel CI/MCSCF
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!     case select
  end subroutine sync_coworkers_cfg

  subroutine sync_coworkers_ci_cfg
!*******************************************************************************
!
!    purpose:  provide co-workers with basic common block/orbital space
!              knowledge in parallel CI/MCSCF runs.
!
!*******************************************************************************
!-------------------------------------------------------------------------------

!     real(8)
      call dalton_mpi_bcast(lucita_cfg_accepted_truncation, 0, mpi_comm_world)
!     logical
      call dalton_mpi_bcast(lucita_cfg_ras1_set,            0, mpi_comm_world)
!     integer
      call dalton_mpi_bcast(ngsh_lucita,                    0, mpi_comm_world)

  end subroutine sync_coworkers_ci_cfg
!******************************************************************************
  subroutine sync_coworkers_abcd
!******************************************************************************
!
!    purpose:  provide co-workers with 1-/2-electron integrals in parallel 
!              CI/MCSCF runs.
!
!******************************************************************************
  end subroutine sync_coworkers_abcd
!******************************************************************************
  
#else 
module dummy_sync_coworkers
#endif
end module
