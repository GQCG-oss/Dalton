   !SIMPLE TEST PROGRAM FOR A FEW FORTRAN 2003 STANDARD TESTS
   program test
     use, intrinsic :: iso_c_binding
     implicit none
     abstract interface
       subroutine testroutine(a,b,c)
         real,intent(in)  :: a,b
         real,intent(out) :: c
       end subroutine
     end interface
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

     success = .true.
    
     funcptr => add

     call funcptr(fptr1(1),fptr1(2),fptr1(7))

     funcptr => mult
     
     call funcptr(fptr1(4),fptr1(5),fptr1(8))

     if(fptr1(3)/=fptr1(7).or.fptr1(6)/=fptr1(8)) success = .false.
     if(c_associated(cptr)) success = .false.
     deallocate(fptr1)
     fptr2 => null()
     

     contains
       subroutine mult(a,b,c)
         real, intent(in)  :: a,b
         real, intent(out) :: c
         c = a * b
       end subroutine mult
       subroutine add(a,b,c)
         real, intent(in)  :: a,b
         real, intent(out) :: c
         c = a + b
       end subroutine add

   end program
