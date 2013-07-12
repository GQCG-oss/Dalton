module unittestDriver
use unittest_util
use unittest_mat

contains

subroutine TestDriver(test)
implicit none
   type(TestItem), intent(inout) :: test

   !Branch out here for what tests to run!
   if (test%run_all) then
      call MatrixTest
      write(*,*) "Only matrix tests implemented"
   else if (test%run_matrix) then
      call MatrixTest
   else
      write(*,*) "No tests chosen"
   endif
end subroutine TestDriver

end module unittestDriver

!> \brief Beginning of a unit testing framework for LSDALTON.
!> \author S. Host
!> \date August 2010 
program unittest
use unittestDriver
use unittest_util
implicit none

   type(TestItem) :: test
   integer        :: choice

   call set_default_TestItem(test)
   
   ! 1. What to test? Should be chosen here. Answer is put in a 
   !    structure, because there could be many questions here, and probably a
   !    lot of info will be carried around.
   write(*,*) "Welcome to the unit testing framework for LSDALTON" 
   write(*,*) "Your options are:"
   write(*,*) "1. Run all unit tests"
   write(*,*) "2. Run unit tests related to matrix framework"
   write(*,*) "Please enter your choice:"
   read(*,*)   choice
   write(*,*)

   call set_TestItem(choice,test)

   call TestDriver(test)
end program unittest


