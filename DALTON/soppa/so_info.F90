module so_info
!
!  This module contains configuration information for the AOSOPPA code.
!  It also contains a few convenience functions.
!
!
   implicit none

!#include"soppinf.h"
   !
   !  implicit.h currently uses real*8 ?? 
   integer, parameter :: sop_dp = kind(1.0D-8)
   !
   real(sop_dp), parameter :: sop_stat_trh = 1.0D-8 ! Threshold for considering a
                                                    ! frequency as zero in linear
                                                    ! response

   real(sop_dp), parameter :: sop_dthresh = 1.0D-4  ! Smallest denominator used when 
                                                    ! preconditioning trial-vectors.

   integer, parameter :: sop_num_models = 8 ! Number of currently defined methods

   ! Parameter for calculation types
   integer, parameter :: sop_linres = 1, & ! linear response
                         sop_excita = 2    ! excitation energy

   ! Integer constants giving the position of each method in the following
   ! arrays. The methods will be used in the order listed here, so order them
   ! such that reuse of results from previous methods makes sense.
   ! Otherwise code outside of this module should not make any assumptions on the
   ! actual values of these constats!
   integer, parameter :: sop_model_rpa     = 1, &
                         sop_model_rpad    = 2, &
                         sop_model_hrpa    = 3, &
                         sop_model_hrpad   = 4, &
                         sop_model_shrpad  = 5, &
                         sop_model_soppa   = 6, &
                         sop_model_sopcc2  = 7, &
                         sop_model_sopccsd = 8

   ! Array of the allowed model labels
   character(len=5), dimension(sop_num_models), parameter :: sop_models = &
                  (/ 'AORPA','DCRPA','AOHRP','DCHRP','SDCHR', &
                     'AOSOP','AOCC2','AOSOC' /)

   ! Array of method full model names names
   character(len=11), dimension(sop_num_models), parameter :: sop_mod_fullname = &
      (/'RPA        ','RPA(D)     ','Higher RPA ','HRPA(D)    ', &
        's-HRPA(D)  ','SOPPA      ','SOPPA(CC2) ','SOPPA(CCSD)'/)

   ! Arrays of arguments the method needs to pass to GET_DENS
   character(len=4), dimension(sop_num_models), parameter :: sop_dens_label = &
      (/'NONE','MP2 ','MP2 ','MP2 ','MP2 ','MP2 ','CC2 ','CCSD'/)

   ! Additional SOPPA filenames (Added here instead of soppinf.h)
   character(len=11), parameter :: FN_RDENS  = 'soppa_densp', &
                                   FN_RDENSE = 'soppa_dense', &
                                   FN_RDENSD = 'soppa_densd'
   character(len=10), dimension(3), parameter :: FN_PROP =    &
      (/ 'soppa_real', 'soppa_imag', 'soppa_trip' /)
   ! Flags stating which models are active
   logical :: AORPA = .false., AOHRP = .false., DCRPA = .false., &
              DCHRP = .false., SDCHR = .false., &
              AOSOP = .false., AOSOC = .false., AOCC2 = .false.
   ! Flag for ensuring that second order, first order density is written
   ! only once
   logical :: sop_mp2ai_done = .false.

   ! Convergence threshold
   real(sop_dp) :: sop_conv_thresh = 1.D-04

   ! Whether to save vectors for Seller's algorithm
   logical :: sop_use_seller = .false.

   interface so_has_doubles
      module procedure so_has_doubles_name, so_has_doubles_num
   end interface

contains

   pure function so_full_name(model)
      !  Returns the name of the model corresponding to
      !  the label "model"
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

   pure function so_needs_densai(model)
      !  Returns true if the "model" needs recalculation of
      !  the O-V part of the density matrix
      character(len=5),intent(in) :: model
      logical :: so_needs_densai
      integer :: i
      do i = 1, sop_num_models
         if ( model .eq. sop_models(i) ) then
            so_needs_densai = (sop_dens_label(i) .eq. 'MP2 ')
            exit
         end if
      end do
      return
   end function

   pure function so_has_doubles_num(nmodel)
      !  Return true if model includes double
      !  excitations, interger argument
      integer, intent(in) :: nmodel
      logical :: so_has_doubles_num
      so_has_doubles_num = any ( nmodel .EQ. (/ sop_model_rpad,   &
                                                sop_model_hrpad,  &
                                                sop_model_shrpad, &
                                                sop_model_soppa,  &
                                                sop_model_sopcc2, &
                                                sop_model_sopccsd /) )
      return
   end function

   pure subroutine so_get_active_models ( list )
      !  Returns list of active models in execution order.
      !  (Work around to the fact that the models are controlled by
      !   individual logical variables)
      logical, intent(out) :: list(sop_num_models)
      list = (/ AORPA, DCRPA, AOHRP, DCHRP, SDCHR, AOSOP, AOCC2, AOSOC /)
      return
   end subroutine

   pure function so_num_active_models ()
      ! Returns number of models treated in this AOSOPPA run
      integer :: so_num_active_models
      logical :: active_models (sop_num_models)
      call so_get_active_models( active_models )
      so_num_active_models = count( active_models )
      return
   end function

   pure function so_any_active_models ()
      logical :: so_any_active_models
      logical :: active_models (sop_num_models)
      call so_get_active_models(active_models)
      so_any_active_models = any(active_models)
      return
   end function

   pure function so_has_doubles_name(model)
      !  Return true if model includes double
      !  excitations
      character(len=5), intent(in) :: model
      logical :: so_has_doubles_name

      so_has_doubles_name = (model.eq.'DCRPA').or.(model.eq.'AOSOP').or. &
                            (model.eq.'AOCC2').or.(model.eq.'AOSOC').or. &
                            (model.eq.'DCHRP').or.(model.eq.'SDCHR')
      return
   end function

   pure function so_double_correction(model)
      !  Return true if model includes double
      !  correction (RPA(D) and HRPA(D)
      character(len=5), intent(in) :: model
      logical :: so_double_correction

      so_double_correction = (model.eq.'DCRPA').or.(model.eq.'DCHRP').or. &
                             (model.eq.'SDCHR')
      return
   end function

   pure function so_singles_second(model)
      !  Return true if model includes second order singles
      !  contributions
      character(len=5), intent(in) :: model
      logical :: so_singles_second

      so_singles_second = (model.eq.'DCRPA').or.(model.eq.'AOSOP').or.&
                          (model.eq.'AOSOC').or.(model.eq.'AOHRP').or.&
                          (model.eq.'AOCC2').or.(model.eq.'DCHRP').or.&
                          (model.eq.'SDCHR')
      return
   end function

   pure function so_singles_first(model)
      !  Returns true if model includes second order singles.
      !  Needed if we calculate perturbation correction.
      character(len=5), intent(in) :: model
      logical :: so_singles_first
      so_singles_first = (model.ne.'DCRPA').or.(model.ne.'DCHRP').or.&
                         (model.ne.'SDCHR')
   end function

   pure function so_model_number(model)
      !  Reverse lookup function to get AOSOPPA
      !  model number from old model labels
      character(len=5), intent(in) :: model
      integer :: i
      integer :: so_model_number
      so_model_number = 0 !! invalid number !
      do i = 1, sop_num_models
         if (model .eq. sop_models(i) ) so_model_number = i
      end do
      return
   end function

end module so_info
