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

#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

  implicit none

  private

  save

  integer, private                       :: istat(MPI_STATUS_SIZE)
  integer, private                       :: ierr

contains 
 
  subroutine sync_coworkers_cfg
!******************************************************************************
!
!    purpose:  synchronize (if necessary) co-workers for parallel CI/MCSCF
!
!*******************************************************************************
!-------------------------------------------------------------------------------
  end subroutine sync_coworkers_cfg

  subroutine sync_coworkers_ijkl
!******************************************************************************
!
!    purpose:  provide co-workers with 1-/2-electron integrals in parallel 
!              CI/MCSCF runs.
!
!******************************************************************************
  end subroutine sync_coworkers_ijkl
!******************************************************************************
  
#else 
module dummy_sync_coworkers
#endif
end module
