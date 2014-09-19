!
! asm september 2014
!
module dyn_iadrpk
!
  integer, dimension(:), allocatable :: iadrpk
  integer :: iadrpk_dim
  logical :: iadrpk_exists
!
  contains
!
     integer function index90(i,j)
       implicit none
       integer :: i,j
       index90 = max(i,j) * (max(i,j)-3) / 2 + i + j
     end function index90
!
end module dyn_iadrpk
!
!
!
subroutine get_iadrpk(luout,nsym,muld2h,nbas,nbast,i2bst,iaodis,iaodpk)
!
   use dyn_iadrpk
!
   implicit none
!
   interface
     subroutine alloc_iadrpk(luout)
       integer :: luout
     end subroutine alloc_iadrpk
   end interface
!
   interface
     subroutine deall_iadrpk(luout)
       integer :: luout
     end subroutine deall_iadrpk
   end interface
!
   integer :: luout, nsym, nbast
   integer, dimension(8) :: nbas, i2bst
   integer, dimension(8,8) :: muld2h, iaodis, iaodpk
!
   integer :: isymab, isyma, isymb, a, b, nab
   integer :: nabsq, nabpk
!
!
   iadrpk_exists = allocated(iadrpk)
   if (iadrpk_exists) call deall_iadrpk(luout)
!
   iadrpk_dim = nbast*nbast
   call alloc_iadrpk(luout)
!
   do isymab = 1,nsym
      do isymb = 1,nsym
         isyma = muld2h(isymab,isymb)
!
         do b = 1,nbas(isymb)
            do a = 1,nbas(isyma)
!
               nabsq = i2bst(isymab) + iaodis(isyma,isymb) &
                     + nbas(isyma)*(b-1) + a
!
               if (isyma .eq. isymb) then
                  nab = index90(a,b)
               else if (isymb .gt. isyma) then
                  nab = nbas(isyma) * (b-1) + a
               else
                  nab = nbas(isymb) * (a-1) + b
               end if
               nabpk = iaodpk(isyma,isymb) + nab
!
               iadrpk(nabsq) = nabpk
!
            end do
         end do
      end do
   end do
!
   return
   end subroutine get_iadrpk
!
!
!
subroutine deall_iadrpk(luout)
!
   use dyn_iadrpk
!
   implicit none
!
   integer :: luout, iadrpk_ok
!
!
   deallocate(iadrpk,stat=iadrpk_ok)
   if (iadrpk_ok .ne. 0) then
      write(luout,*)'Error deallocating existing IADRPK in GET_IADRPK' 
      call quit('Error deallocating existing IADRPK in GET_IADRPK')
!  else
!     write(luout,*) 'Existing IADRPK deallocated by GET_IADRPK'
   end if
!
end subroutine deall_iadrpk
!
!
!
subroutine alloc_iadrpk(luout)
!
   use dyn_iadrpk
!
   implicit none
!
   integer :: luout, iadrpk_ok
!
!
   allocate(iadrpk(iadrpk_dim),stat=iadrpk_ok)
   if (iadrpk_ok .ne. 0) then
      write(luout,*) 'Error allocating IADRPK in GET_IADRPK'
      call quit('Error allocating IADRPK in GET_IADRPK')
!  else
!     write(luout,*) 'IADRPK allocated by GET_IADRPK'
   end if
!
end subroutine alloc_iadrpk
