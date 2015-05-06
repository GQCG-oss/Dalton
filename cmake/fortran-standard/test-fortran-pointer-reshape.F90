!SIMPLE TEST PROGRAM FOR TESTING POINTER RESHAPES                                                                  
module typedefmod
type vec
   real,pointer,contiguous :: t(:)
end type vec

type tile
   real,pointer,contiguous :: t(:) => null()
end type tile

end module typedefmod

program test
  use typedefmod
   implicit none
   real,pointer :: fptr1(:)
   real,pointer :: fptr2(:)
   real,pointer,contiguous :: fptr3(:,:,:)
   logical      :: success
   integer      :: a,b,c,n1,n2
   real,pointer :: matV(:)
   type(vec)    :: vect
   type(tile),pointer   :: tt(:)
   real,pointer,contiguous :: tpm(:,:)
   n1=3
   n2=5
   allocate(tt(3))
   tpm(1:n1,1:n2) => tt(1)%t

   matV(1:5) => vect%t

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
      real, target :: finp(n1,n2)
      real, pointer :: f(:)
      f(1:n1*n2) => finp
   end subroutine test_in_subroutine
  
end program

