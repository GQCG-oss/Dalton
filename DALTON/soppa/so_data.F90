module so_data
   
   use so_info

   implicit none
   private
   public :: link_root, list_element

   type list_element
      character(len=8) :: label
      real(sop_dp) :: freq
      integer :: sym
      type(list_element), pointer :: next => null()
   end type

   type link_root
      integer :: num_elements = 0
      type(list_element), pointer :: first => null(), &
                                     last => null()
      contains
         procedure :: new => new_element
         procedure :: pos => element_position
         procedure :: get => get_element
         procedure :: empty => empty_list
   end type

   type(link_root),public :: fileinf

contains

   function element_position(self,label,freq,sym)
      !
      !  Returns the position of the matching element in
      !  the list
      !
      class(link_root), intent(in) :: self
      character(len=8), intent(in) :: label
      real(sop_dp), intent(in) :: freq
      integer, intent(in) :: sym
      integer :: element_position
      type(list_element), pointer :: current
      integer :: n
      !
      !  Check for any elements
      if (self%num_elements .le. 0) then
         element_position = 0
         return
      end if
      current => self%first
      n = 1
      !
      !  Loop over the elements 
      do
         ! Check current element
         if ( (current%label.eq.label).and. &
              (current%freq.eq.freq) ) then
            element_position = n
            return
         !  Check if we've run out of elements
         elseif ( n .ge. self%num_elements) then
            element_position = 0
            return
         else
            n = n + 1
            current => current%next
         end if
      end do
   end function

   function new_element(self,label,freq,sym)
      !
      !  Adds a new element to the list, if it is not
      !  present.
      !  Returns the index of the element
      !
      !  Arguments
      class(link_root) :: self
      character(len=8), intent(in) :: label
      real(sop_dp), intent(in) :: freq
      integer, intent(in) :: sym
      integer :: new_element
      !
      ! Locals
      integer :: n
      !
      ! Check if the element excists, in that case return that
      n = self%pos(label,freq,sym)
      if ( n .gt. 0 ) then
         new_element = n
         return
      end if
      !
      ! If the list is empty, make the first element
      if (self%num_elements .eq. 0 ) then
         allocate(self%first)
         self%last => self%first
      else
         ! Add the new element
         allocate(self%last%next)
         self%last => self%last%next
      end if
      ! Fill in the info
      self%last%label = label
      self%last%freq = freq
      self%last%sym = sym
      self%num_elements = self%num_elements + 1
      new_element = self%num_elements
      
      return
   end function

   function get_element(self,idx)
      !
      !  Returns the idx'th element of the list
      !
      !  Arguments
      integer, intent(in) :: idx
      class(link_root), intent(in) :: self
      type(list_element), pointer :: get_element
      !
      !  Local variables
      integer :: n
      type(list_element), pointer :: current
      !
      !  Check sanity of input
      if ((idx.le.0).or.(idx.gt.self%num_elements)) then
         get_element => null()
         return
      end if
      !
      !  Loop through the list 
      current => self%first
      do n = 2, idx
         current => current%next
      end do
      get_element => current
      return
   end function

   subroutine empty_list(self)
      !
      !  Deallocate all list elements
      !
      class(link_root),intent(inout) :: self
      !
      type(list_element), pointer :: current, previous
      integer :: n
      self%num_elements = 0
      current => self%first
      ! Loop through the list, deallocating as we go
      do 
         if (.not.associated(current%next)) exit
         previous => current
         current => current%next
         deallocate(previous)
      end do
      deallocate(current)
      self%first => null()
      return
   end subroutine 
end module


