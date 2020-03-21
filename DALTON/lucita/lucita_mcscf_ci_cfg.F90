!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_mcscf_ci_cfg

! stefan: please note that if you modify this module (e.g. add new variables)
!         do not forget to update sync_coworkers_citask (in module sync_coworkers
!         in file par_sync_cw.F90)

  implicit none

  save

! parameter list
  integer, parameter,   public :: vector_exchange_types      =  4

! logical block

  logical,              public :: docisrdft_mc2lu                                      = .false.
  logical,              public :: integrals_from_mcscf_env                             = .false.
  logical,              public :: mcscf_ci_update_ijkl                                 = .true.
  logical,              public :: mcscf_orbital_trial_vector                           = .false.
  logical,              public :: mcscf_ci_trial_vector                                = .false.
  logical,              public :: io2io_vector_exchange_mc2lu_lu2mc                    = .false.
  logical,              public :: cref_is_active_bvec_for_sigma                        = .false.
  logical,              public :: vector_exchange_mc2lu_lu2mc_active                   = .false.
  logical,              public :: vector_update_mc2lu_lu2mc(1:2*vector_exchange_types) = .false.

! character block

! real(8) block

  real(8),              public :: einact_mc2lu               =  0.0d0

! integer block

  integer,              public :: len_cref_mc2lu             =  0
  integer,              public :: len_hc_mc2lu               =  0
  integer,              public :: len_resolution_mat_mc2lu   =  0
  integer,              public :: len_int1_or_rho1_mc2lu     =  0
  integer,              public :: len_int2_or_rho2_mc2lu     =  0

  integer,              public :: vector_exchange_type1      = -1
  integer,              public :: vector_exchange_type2      = -1
end module
