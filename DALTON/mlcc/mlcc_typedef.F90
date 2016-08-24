module mlcc_typedef
!
!
!  mlcc types
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: define the types required for mlcc
!
   implicit none
!
!
   integer, parameter                  :: dp = selected_real_kind(15,307)
!
   type pointer_list
      integer                          :: length, int_length
      real(dp), pointer, dimension(:)  :: point
      integer,  pointer, dimension(:)  :: int_point
      type(pointer_list), pointer      :: next
      type(pointer_list), pointer      :: previous
      logical                          :: in_use
   end type pointer_list
!   
   real(dp), parameter                 :: zero = 0.0D0, one = 1.0D0, two = 2.0D0, half = 0.5D0
   real(dp), parameter                 :: three = 3.0D0, four = 4.0D0, five = 5.0D0, six = 6.0D0
!
!
!  The mlc_dummy subroutine below is an ugly hack to quench these 3 output lines every time rebuilding dalton.x:
!    /opt/local/bin/ranlib: file: lib/libdalton.a(mlcc_typedef.F90.o) has no symbols
!    /opt/local/bin/ranlib: file: lib/libdalton.a(mlcc_typedef.F90.o) has no symbols
!    /opt/local/bin/ranlib: file: lib/libdalton.a(mlcc_typedef.F90.o) has no symbols
!  If anyone has a prettier solution, nice! /hjaaj Aug. 2016
contains
   subroutine mlcc_dummy
   end subroutine mlcc_dummy
end module mlcc_typedef
