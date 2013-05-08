module unittest_util

   integer, save :: TotalTests
   integer, save :: FailedTests
   integer, save :: SuccesfullTests

type TestItem
   logical :: run_all
   logical :: run_matrix

end type TestItem

interface AssertEquals
   module procedure AssertEqInt
   module procedure AssertEqReal
   module procedure AssertEqLogical
   module procedure AssertEqMatrix
end interface AssertEquals

contains 
 
subroutine reset_counters
implicit none

   TotalTests = 0
   FailedTests = 0
   SuccesfullTests = 0 

end subroutine reset_counters

subroutine set_default_testItem(test)
implicit none
   type(TestItem), intent(inout) :: test

   test%run_all = .false.
   test%run_matrix = .false.

end subroutine set_default_testItem

subroutine set_TestItem(choice,test)
implicit none
   integer, intent(in)           :: choice
   type(TestItem), intent(inout) :: test

   if (choice == 1) then
      test%run_all = .true.
   else if (choice == 2) then
      test%run_matrix = .true.
   else
      STOP 'Invalid option chosen!'
   endif

end subroutine set_TestItem

subroutine print_status(string)
implicit none
   character(len=*), intent(in) :: string

   write(*,*) "   =========================================================== "
   write(*,*) "      Status for unit testing of ", adjustr(string)
   write(*,*) "   =========================================================== "
   write(*, '(4x,i3, " tests carried out")'), TotalTests
   write(*, '(4x,"- of these", i3, " failed, and", i3, " were succesfull.")') FailedTests, SuccesfullTests
   write(*,*)
end subroutine print_status

subroutine AssertEqInt(int1, int2, string)
implicit none 
   integer, intent(in) :: int1
   integer, intent(in) :: int2
   character(len=*), intent(in) :: string

   if (int1 == int2) then
      SuccesfullTests = SuccesfullTests + 1
   else
      write(*,*) "TEST FAILED:"
      write(*,*) string
      write(*,*)
      FailedTests = FailedTests + 1
   endif
   TotalTests = TotalTests + 1

end subroutine AssertEqInt 

subroutine AssertEqReal(real1, real2, tol, string)
use precision
implicit none 
   !> First variable to test
   real(realk), intent(in) :: real1
   !> Second variable to test
   real(realk), intent(in) :: real2
   !> Tolerance for equality
   real(realk), intent(in) :: tol
   !> String to print if test fails
   character(len=*), intent(in) :: string

   if (abs(real1 - real2) <= tol) then
      SuccesfullTests = SuccesfullTests + 1
   else
      write(*,*) "TEST FAILED:"
      write(*,*) string
      write(*,*)
      FailedTests = FailedTests + 1
   endif
   TotalTests = TotalTests + 1

end subroutine AssertEqReal

subroutine AssertEqLogical(logical1, logical2, string)
use precision
implicit none 
   !> First variable to test
   logical, intent(in) :: logical1
   !> Second variable to test
   logical, intent(in) :: logical2
   !> String to print if test fails
   character(len=*), intent(in) :: string

   if (logical1 == logical2) then
      SuccesfullTests = SuccesfullTests + 1
   else
      write(*,*) "TEST FAILED:"
      write(*,*) string
      write(*,*)
      FailedTests = FailedTests + 1
   endif
   TotalTests = TotalTests + 1

end subroutine AssertEqLogical

subroutine AssertEqMatrix(mat1, mat2, string)
use matrix_module
implicit none 
   !> First variable to test
   type(matrix), intent(in) :: mat1
   !> Second variable to test
   type(matrix), intent(in) :: mat2
   !> String to print if test fails
   character(len=*), intent(in) :: string
   real(realk) :: sumElms
   integer     :: i

   if (mat1%nrow /= mat2%nrow .or. mat1%ncol /= mat2%ncol) STOP 'Wrong matrix dimensions in AssertEqMatrix'

   sumElms = 0.0d0
   do i = 1, mat1%nrow*mat1%ncol
      sumElms = sumElms + mat1%elms(i) - mat2%elms(i)
   enddo

   if (abs(sumElms) < 1.0d-15) then
      SuccesfullTests = SuccesfullTests + 1
   else
      write(*,*) "TEST FAILED:"
      write(*,*) string
      write(*,*)
      FailedTests = FailedTests + 1
   endif
   TotalTests = TotalTests + 1

end subroutine AssertEqMatrix

end module unittest_util
