!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_integral_density_interface

  implicit none

#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

  public lucita_pointer_integral_driver

contains

!**********************************************************************

  subroutine lucita_pointer_integral_driver(citask_id,                      &
                                            int1_or_rho1,                   &
                                            int2_or_rho2,                   &
                                            update_ijkl_from_env,           &
                                            orbital_trial_vector,           &
                                            ci_trial_vector,                &
                                            mcscf_provides_integrals,       &
                                            print_level)
!*******************************************************************************
!
!    purpose:  interface routine to the LUCITA setup routine which 
!              will create the mandatory integral indices work space pointers 
!              as well as read in integrals (if applicable).
!
!*******************************************************************************
    character (len=12), intent(in)  :: citask_id
    real(8), intent(inout)          :: int1_or_rho1(*)
    real(8), intent(inout)          :: int2_or_rho2(*)
    integer, intent(in)             :: print_level
    logical, intent(inout)          :: update_ijkl_from_env
    logical, intent(in)             :: orbital_trial_vector
    logical, intent(in)             :: ci_trial_vector
    logical, intent(in)             :: mcscf_provides_integrals
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!     setup integral array pointers / read in integrals (if enabled)
      call intim(citask_id,                &
                 int1_or_rho1,             &
                 int2_or_rho2,             &
                 update_ijkl_from_env,     &
                 orbital_trial_vector,     &
                 ci_trial_vector,          &
                 mcscf_provides_integrals)

  end subroutine lucita_pointer_integral_driver
!*******************************************************************************

end module lucita_integral_density_interface
