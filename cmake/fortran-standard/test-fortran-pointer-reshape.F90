!SIMPLE TEST PROGRAM FOR TESTING POINTER RESHAPES

program test
   implicit none
   real,pointer :: fptr1(:)
   real,pointer :: fptr2(:)
   real,pointer,contiguous :: fptr3(:,:,:)
   logical      :: success
   integer      :: a,b,c
   
   allocate(fptr1(12))
   call random_number(fptr1)
   
   !Test pointer reshape II
   fptr3(1:2,1:2,1:2) => fptr1(4:)
   success = .true.

   do a=1,2
      do b=1,2
         do c=1,2
            if(fptr3(a,b,c) /= fptr1(a+(b-1)*2+(c-1)*4+3))then
               success = .false.
            endif
         enddo
      enddo
   enddo

   deallocate(fptr1)

   fptr1 => null()
   fptr3 => null()

   print *,success

   contains

   !PGI can do the upper one but not the subroutine
   subroutine test_in_subroutine(finp,n1,n2)
      implicit none
      integer, intent(in) :: n1,n2
      real, contiguous, target :: finp(n1,n2)
      real, pointer :: f(:)
      f(1:n1*n2) => finp
   end subroutine test_in_subroutine
  
end program

