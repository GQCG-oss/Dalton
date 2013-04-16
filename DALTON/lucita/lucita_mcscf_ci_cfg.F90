!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_mcscf_ci_task

! stefan: please note that if you modify this module (e.g. add new variables)
!         do not forget to update sync_coworkers_citask (in module sync_coworkers
!         in file par_sync_cw.F90)

  implicit none

  save

! character block

  character (len = 72), public :: lucita_citask_id            = 'undefined'

end module
