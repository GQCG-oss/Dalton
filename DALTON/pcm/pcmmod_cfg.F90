module pcmmod_cfg

  implicit none

  save

! logical block
  logical, public :: pcmmod_is_pcm_calculation = .false.
  logical, public :: pcmmod_old_integration    = .false.
  logical, public :: pcmmod_separate           = .false.
! printlevel
  integer, public :: pcmmod_print              = 0

end module
