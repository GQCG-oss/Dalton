!dalton_copyright_start
!
!
!dalton_copyright_end

! module for creating, checking and deleting basic parallel communication models beyond MPI_COMM_WORLD 
! which will be created at the start of a DALTON calculation and could be used throughout the various 
! modules (if wanted)
!
! so far, we will use it in:
! - lucita (after merge with hjaaj-srdft)
! - mcscf with lucita (after merge with hjaaj-srdft)
! - hermit
!
! written by sknecht - linkoeping jan 2013
!
! note for merge with hjaaj-srdft: routines are modified copies of the existing routines for lucita in this branch. 
! (whoever is doing that one day)  merging could be a pain in the a** but i guess i will do that anyway. (for now i
!                                  protect everything "old" with NEW_LUCITA)

module parallel_models

  implicit none

  public parallel_models_initialize
  public parallel_models_finalize
  public check_parallel_models

  private

  save

#ifdef NEW_LUCITA
  logical, public :: lucita_models_enabled = .false.
#endif

contains

  subroutine parallel_models_initialize()
!******************************************************************************
!
!    purpose: initialize parallel communication/file-I/O models
!             (call this module also in sequential runs to allocate temporary 
!             arrays)
!
!*******************************************************************************
#include "maxorb.h"
#include "infpar.h"
!*******************************************************************************
!-------------------------------------------------------------------------------

  end subroutine parallel_models_initialize
!*******************************************************************************

  subroutine parallel_models_finalize()
!******************************************************************************
!
!    purpose: finalize parallel communication/file-I/O models
!             (call this module also in sequential runs to deallocate temporary 
!             arrays)
!
!*******************************************************************************
#include "maxorb.h"
#include "infpar.h"
!*******************************************************************************
!-------------------------------------------------------------------------------

  end subroutine parallel_models_finalize
!*******************************************************************************

  subroutine check_parallel_models()
!******************************************************************************
!
!    purpose: check for existing/enabled parallel communication/file-I/O models
!             and initialize their shutdown (call this module also in sequential 
!             runs to deallocate temporary arrays)
!
!*******************************************************************************
#ifdef NEW_LUCITA
  use ttss_block_module
  use parallel_setup
  use parallel_task_distribution_type_module
#endif
#include "maxorb.h"
#include "infpar.h"
#include "priunit.h"
!*******************************************************************************
    logical :: file_open
!-------------------------------------------------------------------------------

#ifdef NEW_LUCITA
  if(lucita_models_enabled)then
!   possibly free ttss-type used in LUCITA/MCSCF-LUCITA
    call ttss_free(ttss_info)
!   possibly free communicator + file type used in LUCITA/MCSCF-LUCITA
    call lucita_close_parallel_model(.true.)
!   possibly free parallel task list distribution in LUCITA/MCSCF-LUCITA
    call parallel_task_distribution_free_lucipar(ptask_distribution)
    if(mytid > 0)then
      inquire(unit=lupri,opened=file_open)
      if(file_open) close(lupri,status="keep")
    end if
  end if
#endif

  end subroutine check_parallel_models

!*******************************************************************************
end module parallel_models
