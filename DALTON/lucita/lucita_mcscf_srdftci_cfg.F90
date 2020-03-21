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

! integer   
  integer, parameter,   public :: nr_of_weights_max            =  100

! input parameters
  real(8),              public :: weights(1:nr_of_weights_max) =  0.0d0
  integer,              public :: nr_of_weights                = -1

  integer,              public :: srdft_restart_ci             =  0

! parameters

! logical block
  logical,              public :: srdft_ci_1pdens_cref_restore = .false. 
  logical,              public :: srdft_ci_with_lucita         = .false. 
  logical,              public :: do_sc_ensemble_dft           = .false. 

! double precision block

  real(8), allocatable, public :: srdft_srac_lucita(:)
  real(8), allocatable, public :: srdft_cmo_lucita(:)
  real(8),              public :: emydft_mc2lu                 =  0.0d0
  real(8),              public :: emydftaux_mc2lu              =  0.0d0
  real(8),              public :: ejcsr_mc2lu                  =  0.0d0
  real(8),              public :: ejvsr_mc2lu                  =  0.0d0
  real(8),              public :: edsr_mc2lu                   =  0.0d0
  real(8),              public :: edft_mc2lu(3)                =  0.0d0


end module
