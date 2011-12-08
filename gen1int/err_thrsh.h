  ! threshold of error
  real(REALK), parameter :: ERR_THRSH = 10.0_REALK**(-8)
  ! threshold of ratio to the referenced result
  real(REALK), parameter :: RATIO_THRSH(2) = (/ &
    1.0_REALK-10.0_REALK**(-6), 1.0_REALK+10.0_REALK**(-6)/)
