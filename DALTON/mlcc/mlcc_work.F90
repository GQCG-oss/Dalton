module mlcc_work
!
!
!  mlcc work definitions
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: define pointers needed for memory management
!
   use mlcc_typedef
   use mlcc_block_import
!   
   implicit none
!
   real(dp), pointer, private             :: work_point_start(:) => null()
   real(dp), pointer, private             :: work_point_end(:) => null()
   integer, pointer, private              :: work_point_start_int(:) => null()
   integer, pointer, private              :: work_point_end_int(:) => null()
   integer, private                       :: work_length = 0, n_arrays=0
   integer, private                       :: work_remains = 0
   integer, private                       :: work_used = 0
   type(pointer_list), private, pointer   :: list_head,list_tail
   integer, private                       :: real_to_int
!
!
contains
!
   subroutine work_setup(work,int_work,lwork)
!
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: set up pointers to work and end of work
!
      implicit none
!   
      real(dp), target, dimension(:)   :: work
      integer, target, dimension(:)    :: int_work
      integer, intent(in)              :: lwork
      real(dp)                         :: x
!   
      real_to_int = sizeof(x)/sizeof(lwork)
!      
      n_arrays = 0
!      
      work_point_start     => work
      work_point_start_int => int_work
      work_point_end       => work
      work_point_end_int   => int_work
!
      work_length = lwork
      work_remains = lwork
!      
      allocate(list_head)
      list_head%length     = 0
      list_head%int_length = 0
      list_head%point      => null()
      list_head%int_point  => null()
      list_head%next       => null()
      list_head%previous   => null()
      list_head%in_use     = .true.
      list_tail            => list_head
!      
   end subroutine work_setup
!   
   subroutine work_allocator(point,length)
!      
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: allocate arrays of length in work and return pointer
!
      implicit none
!   
      real(dp), dimension(:), pointer  :: point
      integer, intent(in)              :: length
      integer                          :: int_length
      type(pointer_list), pointer      :: current,new_point
!   
      int_length = real_to_int*length
!      
      if (work_remains - length .gt. 0) then
!      
         work_remains = work_remains - length
         work_used    = work_used + length
!      
         point        => work_point_end(1:length)
!
!      
         current => list_head
!
         do while (associated(current%next))
            current => current%next
         end do
!
         allocate(current%next)
!
         new_point            => current%next
         new_point%previous   => current
         new_point%next       => null()
         new_point%length     = length
         new_point%int_length = int_length
         new_point%point      => point
         new_point%int_point  => work_point_end_int(1:int_length)
         new_point%in_use     = .true.
!
         list_tail => new_point
!         
         work_point_end       => work_point_end(length+1:work_remains)
         work_point_end_int   => work_point_end_int(int_length+1:work_remains)
!
      else
         write(lupri,*) 'Not enough memory in work_allocator'
         write(lupri,*) 'Available: ', work_remains
         write(lupri,*) 'Needs:     ', length
         call quit('Not enough work space.')
      end if
!   
   end subroutine work_allocator
!   
!   
   subroutine work_int_allocator(point,length)
!      
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: allocate integer arrays in work. Uses sizeof() to figure out the size
!           of an integer compared to a real(dp). This is Fortran 2008, but at least it's
!           Fortran standard
!
      implicit none
!   
      integer, dimension(:), pointer   :: point
      integer, intent(in)              :: length
      type(pointer_list), pointer      :: current,new_point
      integer                          :: real_length, int_pad, int_length
!   
!     Integer and real storage requirements
!
      int_pad     = length - (length/real_to_int)*real_to_int
      int_length  = length + int_pad
      real_length = int_length/real_to_int
!      
      if (work_remains - real_length .gt. 0) then
!      
         work_remains   = work_remains - real_length
         work_used      = work_used + real_length
!      
         point          => work_point_end_int(1:int_length)
!
         current => list_head
!
         do while (associated(current%next))
            current => current%next
         end do
!
         allocate(current%next)
!
         new_point            => current%next
         new_point%previous   => current
         new_point%next       => null()
         new_point%length     = real_length
         new_point%int_length = int_length
         new_point%point      => work_point_end(1:real_length)
         new_point%int_point  => point
         new_point%in_use     = .true.
!
         list_tail => new_point
!
         work_point_end       => work_point_end(length+1:work_remains)
         work_point_end_int   => work_point_end_int(int_length+1:work_remains)
!      
      else
         write(lupri,*) 'Not enough memory in work_int_allocator'
         write(lupri,*) 'Available: ', work_remains
         write(lupri,*) 'Needs:     ', real_length
         call quit('Not enough work space.')
      end if
!   
   end subroutine work_int_allocator
!   
!
   subroutine work_deallocator(point)
!      
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: deallocate pointers allocated in the allocator
!
      implicit none
!   
      real(dp), pointer, dimension(:)  :: point
!
      if(associated(list_head,target=list_tail)) then
         call quit('No pointers to deallocate in work_deallocater')
      end if
!
      if(associated(point,target=list_tail%point)) then
         work_point_end       => point
         work_point_end_int   => list_tail%int_point
         work_remains         = work_remains + list_tail%length
         work_used            = work_used - list_tail%length
         point                => null()
         list_tail            => list_tail%previous
         deallocate(list_tail%next)
         list_tail%next       => null()
      else
         call quit('point in work_deallocator is not last in list')
      end if

!
   end subroutine work_deallocator
!   
!
   subroutine work_int_deallocator(point)
!      
!  Authors Henrik Koch and Rolf H. Myhre
!  January 2015
!
!  Purpose: deallocate integer pointers allocated in the allocator
!
      implicit none
!   
      integer, pointer, dimension(:)   :: point
!   
!     Integer and real storage requirements
!
      if(associated(list_head,target=list_tail)) then
         call quit('No pointers to deallocate in work_int_deallocater')
      end if
!
      if(associated(point,target=list_tail%int_point)) then
!  
         work_point_end       => list_tail%point
         work_point_end_int   => point
         work_remains         = work_remains + list_tail%length
         work_used            = work_used - list_tail%length
         point                => null()
         list_tail            => list_tail%previous
         deallocate(list_tail%next)
         list_tail%next       => null()
      else
         call quit('point in work_int_deallocator is not last in list')
      end if

!
   end subroutine work_int_deallocator
!
!
   subroutine work_info
!      
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: print info about contents in work_list
!
      implicit none
!      
      integer                       :: counter
      type(pointer_list), pointer   :: current
!      
      counter = 0
      current => list_head
!
      write(lupri,*)
      write(lupri,*) 'free space: ', work_remains
      write(lupri,*) 'work used:  ', work_used
      write(lupri,*)
      do while (associated(current%next))
!         
         counter = counter + 1
         current => current%next
!         
         write(lupri,*) 'counter: ', counter
         write(lupri,*) 'length:  ', current%length
         write(lupri,*)
!         
      end do
!      
      call flshfo(lupri)
!      
   end subroutine work_info
!  
   function work_end()
!      
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: return pointer to first free element in work
!
      implicit none
!      
      real(dp), dimension(:), pointer :: work_end
! 
      work_end => work_point_end
!      
   end function work_end
!   
   function work_free()
!      
!  Authors Henrik Koch and Rolf H. Myhre
!  December 2014
!
!  Purpose: resturn size of free work
!
      implicit none
!      
      integer :: work_free
! 
      work_free = work_remains
!      
   end function work_free
!   
!  
end module mlcc_work

