!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_mcscf_srdftci_cfg

!
! add some explanation(s)...
!

  implicit none

  save

! parameters

! logical block
  logical,              public :: srdft_ci_1pdens_cref_restore = .false. 

! double precision block

  real(8), allocatable, public :: srdft_srac_lucita(:)
  real(8), allocatable, public :: srdft_cmo_lucita(:)
  real(8),              public :: emydft_mc2lu                 =  0.0d0
  real(8),              public :: ejcsr_mc2lu                  =  0.0d0
  real(8),              public :: ejvsr_mc2lu                  =  0.0d0
  real(8),              public :: edsr_mc2lu                   =  0.0d0
  real(8),              public :: edft_mc2lu                   =  0.0d0

end module
