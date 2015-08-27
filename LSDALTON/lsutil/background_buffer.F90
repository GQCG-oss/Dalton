!\author Patrick Ettenhuber
! The purpose of this module is to provide a background allcated space of
! memory, which the application may use. This can prevent alloc-dealloc overhead
! and prevent the memory from getting fragmented
module background_buffer_module
  use,intrinsic :: iso_c_binding,only:c_ptr,c_null_ptr
  use precision

  private

  public :: max_n_pointers,buf_realk,&
       & mem_is_background_buf_init,mem_get_bg_buf_n,mem_get_bg_buf_free

  ! maximum number of pointer to be associated ot bg buf, -> allocated on stack,
  ! may be changed if necessary
  integer, parameter   :: max_n_pointers = 500

  type bg_buf_realk
     logical              :: init = .false.
     integer              :: offset
     integer              :: nmax
     integer              :: max_usage
     real(realk), pointer :: p(:)
     type(c_ptr)          :: c = c_null_ptr
     integer              :: n
     integer              :: f_addr(max_n_pointers)
     type(c_ptr)          :: c_addr(max_n_pointers)
     !all deletion handling counters
     type(c_ptr)          :: c_mdel(max_n_pointers)
     integer              :: e_mdel(max_n_pointers)
     integer              :: n_mdel
     integer              :: n_prev
     logical              :: l_mdel
     type(c_ptr)          :: f_mdel
  end type bg_buf_realk

  save

  type(bg_buf_realk) :: buf_realk

  contains
  function mem_is_background_buf_init() result(init)
     implicit none
     logical :: init
     init = buf_realk%init
  end function mem_is_background_buf_init

  function mem_get_bg_buf_n() result(n)
     implicit none
     integer(kind=8) :: n
     n = buf_realk%nmax
  end function mem_get_bg_buf_n

  function mem_get_bg_buf_free() result(n)
     implicit none
     integer(kind=8) :: n
     n = buf_realk%nmax-buf_realk%offset
  end function mem_get_bg_buf_free


end module background_buffer_module
