!\author Patrick Ettenhuber
! The purpose of this module is to provide a background allcated space of
! memory, which the application may use. This can prevent alloc-dealloc overhead
! and prevent the memory from getting fragmented
module tensor_bg_buf_module
  use,intrinsic :: iso_c_binding,only:c_ptr,c_null_ptr,c_loc,c_associated

  use tensor_error_handler
  use tensor_parameters_and_counters

  public tensor_bg_init, tensor_bg_n, tensor_bg_free
  public tensor_bg_buf_dp_type
#ifdef TENSORS_IN_LSDALTON
  public tensor_max_bg_pointers
  public buf_tensor_dp
  public tensor_initialize_bg_buf_from_lsdalton_bg_buf
  public tensor_free_bg_buf
#endif


  private

  type tensor_bg_buf_dp_type
     logical, pointer     :: init
     integer, pointer     :: offset
     integer, pointer     :: nmax
     integer, pointer     :: max_usage
     real(tensor_dp), pointer :: p(:)
     type(c_ptr), pointer :: c
     integer, pointer     :: n
     integer, pointer     :: f_addr(:)
     type(c_ptr), pointer :: c_addr(:)
     !all deletion handling counters
     type(c_ptr), pointer :: c_mdel(:)
     integer, pointer     :: e_mdel(:)
     integer, pointer     :: n_mdel
     integer, pointer     :: n_prev
     logical, pointer     :: l_mdel
     type(c_ptr), pointer :: f_mdel
     integer :: void
     contains
     procedure :: clear_md => tensor_clear_mdel_bg_buf
  end type tensor_bg_buf_dp_type

  save

  type(tensor_bg_buf_dp_type) :: buf_tensor_dp
  integer, pointer            :: tensor_max_bg_pointers

  contains

  subroutine tensor_initialize_bg_buf_from_lsdalton_bg_buf(max_pointers,init,offset,nmax,max_usage,p,&
        &c,n,f_addr,c_addr,c_mdel,e_mdel,n_mdel,n_prev,l_mdel,f_mdel)
     integer, target     :: max_pointers
     logical, target     :: init
     integer, target     :: offset
     integer, target     :: nmax
     integer, target     :: max_usage
     real(tensor_dp), target :: p(:)
     type(c_ptr), target :: c
     integer, target     :: n
     integer, target     :: f_addr(:)
     type(c_ptr), target :: c_addr(:)
     type(c_ptr), target :: c_mdel(:)
     integer, target     :: e_mdel(:)
     integer, target     :: n_mdel
     integer, target     :: n_prev
     logical, target     :: l_mdel
     type(c_ptr), target :: f_mdel
     tensor_max_bg_pointers   =>   max_pointers
     buf_tensor_dp%init       =>   init         
     buf_tensor_dp%offset     =>   offset
     buf_tensor_dp%nmax       =>   nmax
     buf_tensor_dp%max_usage  =>   max_usage
     buf_tensor_dp%p          =>   p(:)
     buf_tensor_dp%c          =>   c
     buf_tensor_dp%n          =>   n
     buf_tensor_dp%f_addr     =>   f_addr 
     buf_tensor_dp%c_addr     =>   c_addr
     buf_tensor_dp%c_mdel     =>   c_mdel
     buf_tensor_dp%e_mdel     =>   e_mdel
     buf_tensor_dp%n_mdel     =>   n_mdel
     buf_tensor_dp%n_prev     =>   n_prev
     buf_tensor_dp%l_mdel     =>   l_mdel
     buf_tensor_dp%f_mdel     =>   f_mdel
  end subroutine tensor_initialize_bg_buf_from_lsdalton_bg_buf
  subroutine tensor_free_bg_buf()
     implicit none
     tensor_max_bg_pointers   =>  null()
     buf_tensor_dp%init       =>  null() 
     buf_tensor_dp%offset     =>  null() 
     buf_tensor_dp%nmax       =>  null() 
     buf_tensor_dp%max_usage  =>  null() 
     buf_tensor_dp%p          =>  null() 
     buf_tensor_dp%c          =>  null() 
     buf_tensor_dp%n          =>  null() 
     buf_tensor_dp%f_addr     =>  null() 
     buf_tensor_dp%c_addr     =>  null() 
     buf_tensor_dp%c_mdel     =>  null() 
     buf_tensor_dp%e_mdel     =>  null() 
     buf_tensor_dp%n_mdel     =>  null() 
     buf_tensor_dp%n_prev     =>  null() 
     buf_tensor_dp%l_mdel     =>  null() 
     buf_tensor_dp%f_mdel     =>  null() 
  end subroutine tensor_free_bg_buf
  subroutine tensor_clear_mdel_bg_buf(this)
     implicit none
     class(tensor_bg_buf_dp_type), intent(inout) :: this
     integer :: i
     if (this%l_mdel) then
        do i=1,this%n_mdel

           if(.not.c_associated(this%c_mdel(i),c_loc(this%p(this%f_addr(this%n-1)))))then
              call tensor_status_quit("ERROR(tensor_clear_mdel_bg_buf): wrong&
                 & sequence of&
                 & deallocating, make sure you dealloc in the opposite seqence as&
                 & allocating, also when using mark_deleted",-1)
           endif

           this%n = this%n - 1
           this%offset = this%offset - this%e_mdel(i)

        enddo

        this%c_mdel = c_null_ptr
        this%f_mdel = c_null_ptr
        this%e_mdel = 0
        this%n_mdel = 0
        this%l_mdel = .false.
     endif
  end subroutine tensor_clear_mdel_bg_buf

  function tensor_bg_init() result(init)
     implicit none
     logical :: init
     if(associated(buf_tensor_dp%init))then
        init = buf_tensor_dp%init
     else
        init = .false.
     endif
  end function tensor_bg_init

  function tensor_bg_n() result(n)
     implicit none
     integer(kind=tensor_long_int) :: n
     if(associated(buf_tensor_dp%nmax))then
        n = buf_tensor_dp%nmax
     else
        n = 0
     endif
  end function tensor_bg_n

  function tensor_bg_free() result(n)
     implicit none
     integer(kind=tensor_long_int) :: n
     if(associated(buf_tensor_dp%nmax).and.associated(buf_tensor_dp%offset))then
        n = buf_tensor_dp%nmax-buf_tensor_dp%offset
     else
        n = 0
     endif
  end function tensor_bg_free


end module tensor_bg_buf_module
