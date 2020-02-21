!dalton_copyright_start
!
!
!dalton_copyright_end

! module for creating, checking and deleting basic parallel communication models beyond MPI_COMM_WORLD 
! which will be created at the start of the host program and could be used throughout the various 
! modules (if wanted)
!
! so far, we will use it in Dalton:
! - lucita (after merge with hjaaj-srdft)
! - mcscf with lucita (after merge with hjaaj-srdft)
! - hermit
!
! written by sknecht - linkoeping jan 2014
!
! note for merge with hjaaj-srdft in Dalton: routines are modified copies of the existing routines for lucita in this branch. 
! (whoever is doing that one day)            merging could be a pain in the a** but i guess i will do that anyway. (for now i
!                                            protect everything "old" with NEW_LUCITA)

module parallel_models_mpi

  use parallel_communication_models_mpi
  use parallel_models_lucita
! use parallel_file_io_models_mpi

#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif
#else
  implicit none
  integer, parameter :: mpi_comm_world = -1
#endif


  public parallel_models_initialize_mpi
  public parallel_models_finalize_mpi

  logical, public :: lucita_models_enabled = .false.

  private

  save

contains

  subroutine parallel_models_initialize_mpi()
!******************************************************************************
!
!    purpose: initialize parallel communication/file-I/O models
!             (call this module also in sequential runs to allocate temporary 
!             arrays)
!
!*******************************************************************************
!-------------------------------------------------------------------------------
  integer                        :: nr_of_process_glb,my_process_id_glb
  integer                        :: communication_glb_world  ! internal global communicator handle

#ifdef VAR_MPI
  integer(kind=MPI_INTEGER_KIND) :: nr_procs_glb_mpi,my_proc_id_glb_mpi
  integer(kind=MPI_INTEGER_KIND) :: ierr
  logical(kind=MPI_INTEGER_KIND) :: is_init
!-------------------------------------------------------------------------------

      call mpi_initialized(is_init, ierr)
      if(.not.is_init)then 
        print *, ' warning: initialization of parallel communication models requires MPI_init to be called first.'
        print *, ' warning: assuming DALTON was compiled with MPI support but run in serial mode... continuing...'
        my_process_id_glb = 0
        nr_of_process_glb = 1
      else
        call mpi_comm_rank(mpi_comm_world, my_proc_id_glb_mpi, ierr) 
        call mpi_comm_size(mpi_comm_world, nr_procs_glb_mpi, ierr) 
        my_process_id_glb = my_proc_id_glb_mpi
        nr_of_process_glb = nr_procs_glb_mpi
      end if
#else
      my_process_id_glb = 0
      nr_of_process_glb = 1
#endif

!     initialize parallel communication models
      communication_glb_world = MPI_COMM_WORLD
      call communication_init_mpi(communication_info_mpi, &
                                  my_process_id_glb,      &
                                  nr_of_process_glb,      &
                                  communication_glb_world)

  end subroutine parallel_models_initialize_mpi
!*******************************************************************************

  subroutine parallel_models_finalize_mpi(nr_of_process_glb)
!******************************************************************************
!
!    purpose: finalize parallel communication/file-I/O models
!             (call this module also in sequential runs to deallocate temporary 
!             arrays)
!
!*******************************************************************************
  integer, intent(in) :: nr_of_process_glb
!-------------------------------------------------------------------------------

!     finalize parallel communication models
      call communication_free_mpi(communication_info_mpi,  &
                                  nr_of_process_glb)
      if(lucita_models_enabled) call check_parallel_models_mpi('lucita')

  end subroutine parallel_models_finalize_mpi
!*******************************************************************************

  subroutine check_parallel_models_mpi(models_check)
!******************************************************************************
!
!    purpose: check for existing/enabled parallel communication/file-I/O models
!             and initialize their shutdown (call this module also in sequential 
!             runs to deallocate temporary arrays)
!
!*******************************************************************************
  character*(*), intent(in) :: models_check
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      select case(models_check)
        case('all')
!         call check_all_models
          call quit('parallel model check for all modules not implemented yet')
        case('lucita')
          call check_lucita_models()
        case default 
          call quit('parallel model check for all modules not implemented yet')
      end select

  end subroutine check_parallel_models_mpi
!*******************************************************************************

end module parallel_models_mpi
