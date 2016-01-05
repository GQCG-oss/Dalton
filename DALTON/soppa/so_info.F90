module so_info
!
!  This module contains configuration information for the AOSOPPA code.
!  It also contains a few routines convenience functions. 
!  
!
   implicit none

   !
   real(8), parameter :: sop_stat_trh = 1.0D-8 ! Threshold for considering a
                                               ! frequency as zero in linear
                                               ! response

   integer, parameter :: sop_num_models = 5

   ! Parameter for calculation types
   integer, parameter :: sop_linres = 1, & ! linear response
                         sop_excita = 2    ! excitation energy

   ! Array of the allowed model labels
   character(len=5), dimension(sop_num_models), parameter :: sop_models = &
                  (/ 'AORPA','DCRPA','AOHRP','AOSOP','AOSOC' /) 

   ! Array of method full model names names 
   character(len=11), dimension(sop_num_models), parameter :: sop_mod_fullname = &
      (/'RPA        ','RPA(D)     ','Higher RPA ','SOPPA      ','SOPPA(CCSD)'/)

   ! Additional SOPPA filenames (Added here instead of soppinf.h)
   character(len=11), parameter :: FN_RDENS  = 'soppa_densp', &
                                   FN_RDENSE = 'soppa_dense', &
                                   FN_RDENSD = 'soppa_densd'
      
   
contains

   pure function so_full_name(model)
      !
      !  Returns the name of the model corresponding to
      !  the label "model"
      !
      character(len=5),intent(in) :: model
      character(len=11)  :: so_full_name
      integer :: i

      so_full_name = '!!INVALID!!'
      do i = 1, sop_num_models
         if ( model .eq. sop_models(i) ) then
            so_full_name = sop_mod_fullname(i)
            exit
         end if
      end do

      return
   end function

   pure function so_has_doubles(model)
      ! 
      !  Return true if model includes double
      !  excitations
      !
      character(len=5), intent(in) :: model
      logical :: so_has_doubles

      so_has_doubles = (model.eq.'DCRPA').or.(model.eq.'AOSOP').or.(model.eq.'AOSOC')
      return
   end function

   pure function so_singles_second(model)
      ! 
      !  Return true if model includes second order singles
      !  contributions
      !
      character(len=5), intent(in) :: model
      logical :: so_singles_second

      so_singles_second = (model.eq.'DCRPA').or.(model.eq.'AOSOP').or.&
                          (model.eq.'AOSOC').or.(model.eq.'AOHRP')
      return
   end function

   pure function so_singles_first(model) 
      !
      !  Returns true if model includes second order singles.
      !  Needed if we calculate perturbation correction.
      character(len=5), intent(in) :: model
      logical :: so_singles_first
      
      so_singles_first = .not.(model.eq.'DCRPA')
   end function

end module so_info
