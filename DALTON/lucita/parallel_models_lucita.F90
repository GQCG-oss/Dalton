!dalton_copyright_start
!
!
!dalton_copyright_end

module parallel_models_lucita

  implicit none

  public check_lucita_models

  private

  save

contains

  subroutine check_lucita_models()
!******************************************************************************
!
!    purpose: check for existing/enabled parallel communication/file-I/O models
!             and initialize their shutdown (call this module also in sequential 
!             runs to deallocate temporary arrays)
!
!*******************************************************************************
  use ttss_block_module
  use parallel_setup
  use parallel_task_distribution_type_module
#include "maxorb.h"
#include "infpar.h"
#include "priunit.h"
!*******************************************************************************
    logical :: file_open
!-------------------------------------------------------------------------------

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

  end subroutine check_lucita_models

!*******************************************************************************
end module parallel_models_lucita
