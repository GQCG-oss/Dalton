module mlcc3_various
!
!
!  mlcc3 routines to calculate fock matrix and other intermediates
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
   use mlcc_typedef
   use mlcc3_data
   use mlcc_block_import
   use mlcc_work
!
!   
contains
!      
   subroutine mlcc3_square_packed(packed,unpacked,dime)
!
!  Square a packed matrix
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: take in a packed square matrix and return the unpacked one
!
      implicit none
!      
      integer,intent(in)                                 :: dime
      integer                                            :: a, b, inda, indb, indp
      real(dp), dimension(dime*(dime+1)/2), intent(in)   :: packed
      real(dp), dimension(dime**2), intent(inout)        :: unpacked
!
!     loop over a and b .le. a
!
!$omp parallel do private(b,a,inda,indb,indp)
      do a = 1,dime
         do b = 1,a
!
            inda = (a-1)*dime + b
            indb = (b-1)*dime + a
            indp = max(a,b)*(max(a,b)-3)/2 + a + b
!            
            unpacked(inda) = packed(indp)
!            
            if (b .ne. a) then
               unpacked(indb) = packed(indp)
            end if
!            
         end do
      end do
!$omp end parallel do
!            
   end subroutine mlcc3_square_packed
!
   function new_unit()
!      
!     Return a new unit
!     Authors Henrik Koch and Rolf H. Myhre
!     January 2015
!
!     Return a new I/O unit. Adapted from the fortran wiki: http://fortranwiki.org/fortran/show/newunit
!
      implicit none
!      
      integer              :: new_unit
      integer, parameter   :: lun_min=10, lun_max=1000
      logical              :: opened
      integer :: lun
!      
      new_unit=-1
!      
      do lun=lun_min,lun_max
         inquire(unit=lun,opened=opened)
         if (.not. opened) then
            new_unit=lun
            exit
         end if
      end do
!
      if (new_unit .eq. -1) then
         call quit('No free print units found in new_unit')
      end if
!      
   end function new_unit
!   
end module mlcc3_various

