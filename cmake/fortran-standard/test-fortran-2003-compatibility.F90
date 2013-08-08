   !SIMPLE TEST PROGRAM FOR A FEW FORTRAN 2003 STANDARD TESTS
  
   module test_module
     use,intrinsic :: iso_c_binding
     implicit none

     interface add
       module procedure add_int,add_real
     end interface add
     interface mult
       module procedure mult_int,mult_real
     end interface mult

     abstract interface
       subroutine testroutine(a,b,c)
         real,intent(in)  :: a,b
         real,intent(out) :: c
       end subroutine
     end interface

     procedure(mult_real),pointer :: mp => null()
     procedure(add_real),pointer  :: ap => null()


     contains

     subroutine my_test
       implicit none
       type(c_ptr)  :: cptr
       real,pointer :: fptr1(:)
       real,pointer :: fptr2(:)
       logical      :: success
       procedure(testroutine), pointer :: funcptr => null()
      
       allocate(fptr1(8))
      
       cptr = c_loc(fptr1(1))
      
       if(c_associated(cptr))then
         call c_f_pointer(cptr,fptr2,[8])
       endif
      
       cptr = c_null_ptr
      
       if(associated(fptr2))then
         call random_number(fptr2)
         fptr2(3) = fptr2(1) + fptr2(2)
         fptr2(6) = fptr2(4) * fptr2(5)
       endif
       mp => mult_real
       ap => add_real
      
       success = .true.
     
       funcptr => ap
      
       call funcptr(fptr1(1),fptr1(2),fptr1(7))
      
       funcptr => mp
       
       call funcptr(fptr1(4),fptr1(5),fptr1(8))
      
       if(fptr1(3)/=fptr1(7).or.fptr1(6)/=fptr1(8)) success = .false.
       if(c_associated(cptr)) success = .false.
       deallocate(fptr1)
       fptr2 => null()
       print *,success
     end subroutine my_test

     subroutine mult_real(a,b,c)
       real, intent(in)  :: a,b
       real, intent(out) :: c
       c = a * b
     end subroutine mult_real
     subroutine mult_int(a,b,c)
       integer, intent(in)  :: a,b
       integer, intent(out) :: c
       c = a * b
     end subroutine mult_int

     subroutine add_real(a,b,c)
       real, intent(in)  :: a,b
       real, intent(out) :: c
       c = a + b
     end subroutine add_real
     subroutine add_int(a,b,c)
       integer, intent(in)  :: a,b
       integer, intent(out) :: c
       c = a + b
     end subroutine add_int

   end module test_module

   program test
     use test_module     
     implicit none

     call my_test
     
   end program

