module so_info

   implicit none

   !
   real(8), parameter :: sop_stat_trh = 1.0D-8 ! Threshold for considering a
                                               ! frequency as zero in linear
                                               ! response

   ! Parameter for calculation types
   integer, parameter :: sop_linres = 1, & ! linear response
                         sop_excita = 2    ! excitation energy

   ! Array of the allowed model labels
   character(len=5), dimension(5), parameter :: sop_models = &
                  (/ 'AORPA','DCRPA','AOHRP','AOSOP','AOSOC' /) 

   ! Array of method full model names names 
   character(len=11), dimension(5), parameter :: sot_mod_fullname = &
      (/'RPA        ','RPA(D)     ','Higher RPA ','SOPPA      ','SOPPA(CCSD)'/)

   ! Additional SOPPA filenames (Added here instead of soppinf.h)
   character(len=11), parameter :: FN_RDENS  = 'soppa_densp', &
                                   FN_RDENSE = 'soppa_dense', &
                                   FN_RDENSD = 'soppa_densd'


end module so_info
