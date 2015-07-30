module tensor_allocator
   use, intrinsic :: iso_c_binding

#ifdef TENSORS_IN_LSDALTON
   use background_buffer_module, only: mem_is_background_buf_init
   use memory_handling, only: mem_pseudo_alloc, mem_pseudo_dealloc
#endif

   use tensor_error_handler
   use tensor_parameters_and_counters
   use tensor_bg_buf_module
   use tensor_mpi_binding_module
   use tensor_type_def_module

   !Routines for the allocation and deallocation within the tensor module
   public :: tensor_alloc_mem
   public :: tensor_free_mem
   ! Memory inquiry funcions
   public :: tensor_get_total_mem
   public :: tensor_get_total_heap_mem
   public :: tensor_get_total_bg_mem

   private

   interface tensor_alloc_mem
      module procedure tensor_allocate_tensor_dp_1d,&
                      &tensor_allocate_tensor_dp_1d_std,&
                      &tensor_allocate_tensor_dp_2d_ll,&
                      &tensor_allocate_tensor_dp_2d_ls,&
                      &tensor_allocate_tensor_dp_2d_sl,&
                      &tensor_allocate_tensor_dp_2d_ss,&
                      &tensor_allocate_tensor_dp_1d_mpi_std,&
                      &tensor_allocate_tensor_dp_1d_mpi,&
                      &tensor_allocate_tensor_long_int_1d,&
                      &tensor_allocate_tensor_long_int_1d_std,&
                      &tensor_allocate_tensor_long_int_2d_ll,&
                      &tensor_allocate_tensor_long_int_2d_ls,&
                      &tensor_allocate_tensor_long_int_2d_sl,&
                      &tensor_allocate_tensor_long_int_2d_ss,&
                      &tensor_allocate_tensor_standard_int_1d,&
                      &tensor_allocate_tensor_standard_int_1d_std,&
                      &tensor_allocate_tensor_standard_int_2d_ll,&
                      &tensor_allocate_tensor_standard_int_2d_ls,&
                      &tensor_allocate_tensor_standard_int_2d_sl,&
                      &tensor_allocate_tensor_standard_int_2d_ss,&
                      &tensor_allocate_tensor_long_log_1d,&
                      &tensor_allocate_tensor_long_log_1d_std,&
                      &tensor_allocate_tensor_standard_log_1d,&
                      &tensor_allocate_tensor_standard_log_1d_std,&
                      &tensor_allocate_character_1d,&
                      &tensor_allocate_character_1d_std,&
                      &tensor_allocate_tile_1d,&
                      &tensor_allocate_tile_1d_std, &
                      &tensor_allocate_tensor_1d,&
                      &tensor_allocate_tensor_1d_std
   end interface tensor_alloc_mem

   interface tensor_free_mem
      module procedure tensor_free_tensor_dp_1d,&
                      &tensor_free_tensor_dp_2d,&
                      &tensor_free_tensor_dp_1d_mpi,&
                      &tensor_free_tensor_long_int_1d,&
                      &tensor_free_tensor_long_int_2d,&
                      &tensor_free_tensor_standard_int_1d,&
                      &tensor_free_tensor_standard_int_2d,&
                      &tensor_free_tensor_long_log_1d,&
                      &tensor_free_tensor_standard_log_1d,&
                      &tensor_free_character_1d,&
                      &tensor_free_tile_1d, &
                      &tensor_free_tensor_1d
   end interface tensor_free_mem


   abstract interface
      subroutine tensor_aif_long_int(p,idx,stat)
         use tensor_parameters_and_counters
         import
         implicit none
         integer(kind=tensor_long_int),pointer, intent(inout) :: p(:)
         integer(kind=tensor_standard_int), intent(in) :: idx
         integer, intent(out), optional                :: stat
      end subroutine tensor_aif_long_int
      subroutine tensor_aif_standard_int(p,idx,stat)
         use tensor_parameters_and_counters
         import
         implicit none
         integer(kind=tensor_standard_int),pointer, intent(inout) :: p(:)
         integer(kind=tensor_standard_int), intent(in) :: idx
         integer, intent(out), optional                :: stat
      end subroutine tensor_aif_standard_int
   end interface

   contains


   !get memory currents
   function tensor_get_total_mem() result(bytes)
      implicit none
      integer(kind=tensor_long_int) :: bytes
      integer :: b_heap, b_bg

      b_heap = tensor_get_total_heap_mem()
      b_bg   = tensor_get_total_bg_mem()

      bytes = b_heap + b_bg
   end function tensor_get_total_mem
   function tensor_get_total_heap_mem() result(bytes)
      implicit none
      integer(kind=tensor_long_int) :: bytes
      integer :: i
      bytes = 0
      do i=1,tensor_nmem_idx
         bytes = bytes + counters(i)%curr_
      enddo
   end function tensor_get_total_heap_mem
   function tensor_get_total_bg_mem() result(bytes)
      implicit none
      integer(kind=tensor_long_int) :: bytes
      integer :: i
      bytes = 0
      do i=1,tensor_nmem_idx
         bytes = bytes + counters_bg(i)%curr_
      enddo
   end function tensor_get_total_bg_mem



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!              ATOMIC TYPE: DOUPLE PRECISION: tensor_dp                     !!!!!!!!!!!!!!!!!!!!!!!!                                      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !ALLOCATION INTERFACES
   subroutine tensor_allocate_tensor_dp_1d_std(p,n1,bg,stat)
      implicit none
      real(tensor_dp), pointer, intent(inout) :: p(:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      logical, intent(in),  optional         :: bg
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      logical :: bg_
      bg_=.false.
      if(present(bg)) bg_ = bg
      
      n = n1
      if(bg_)then
         call tensor_allocate_tensor_dp_basic_bg(p,n,tensor_mem_idx_tensor_dp,buf_tensor_dp,stat=stat)
      else
         call tensor_allocate_tensor_dp_basic(p,n,tensor_mem_idx_tensor_dp,stat=stat)
      endif

   end subroutine tensor_allocate_tensor_dp_1d_std
   subroutine tensor_allocate_tensor_dp_1d(p,n1,bg,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p(:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      logical, intent(in),  optional         :: bg
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      logical :: bg_
      bg_=.false.
      if(present(bg)) bg_ = bg
      
      n = n1
      if(bg_)then
         call tensor_allocate_tensor_dp_basic_bg(p,n,tensor_mem_idx_tensor_dp,buf_tensor_dp,stat=stat)
      else
         call tensor_allocate_tensor_dp_basic(p,n,tensor_mem_idx_tensor_dp,stat=stat)
      endif

   end subroutine tensor_allocate_tensor_dp_1d
   subroutine tensor_allocate_tensor_dp_1d_mpi_std(p,c,n1,bg,stat)
      implicit none
      real(tensor_dp), pointer, intent(inout) :: p(:)
      type(c_ptr), intent(inout)              :: c
      integer(kind=tensor_standard_int), intent(in) :: n1
      logical, intent(in),  optional         :: bg
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      logical :: bg_
      bg_=.false.
      if(present(bg)) bg_ = bg
      
      n = n1
      if(bg_)then
         call tensor_allocate_tensor_dp_basic_bg(p,n,tensor_mem_idx_tensor_dp_mpi,buf_tensor_dp,stat=stat)
         c = c_loc(p(1))
      else
         call tensor_allocate_tensor_dp_basic_mpi(p,c,n,tensor_mem_idx_tensor_dp_mpi,stat=stat)
      endif

   end subroutine tensor_allocate_tensor_dp_1d_mpi_std
   subroutine tensor_allocate_tensor_dp_1d_mpi(p,c,n1,bg,stat)
      implicit none
      real(tensor_dp), pointer, intent(inout) :: p(:)
      type(c_ptr), intent(inout)              :: c
      integer(kind=tensor_long_int), intent(in) :: n1
      logical, intent(in),  optional         :: bg
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      logical :: bg_
      bg_=.false.
      if(present(bg)) bg_ = bg
      
      n = n1
      if(bg_)then
         call tensor_allocate_tensor_dp_basic_bg(p,n,tensor_mem_idx_tensor_dp_mpi,buf_tensor_dp,stat=stat)
         c = c_loc(p(1))
      else
         call tensor_allocate_tensor_dp_basic_mpi(p,c,n,tensor_mem_idx_tensor_dp_mpi,stat=stat)
      endif

   end subroutine tensor_allocate_tensor_dp_1d_mpi

   !2D allocations
   subroutine tensor_allocate_tensor_dp_2d_ll(p2,n1,n2,bg,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer(kind=tensor_long_int), intent(in)   :: n2
      logical, intent(in),  optional         :: bg
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n, n1l, n2l
      logical :: bg_
      bg_=.false.
      if(present(bg)) bg_ = bg
      
      n1l = n1
      n2l = n2

      if(bg_)then
         call tensor_status_quit("ERROR(2d alloc): not implemented with bg bufffer",3434)
      else
         call tensor_allocate_tensor_dp_2d_basic(p2,n1l,n2l,tensor_mem_idx_tensor_dp,stat=stat)
      endif

   end subroutine tensor_allocate_tensor_dp_2d_ll
   subroutine tensor_allocate_tensor_dp_2d_sl(p2,n1,n2,bg,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer(kind=tensor_long_int), intent(in)   :: n2
      logical, intent(in),  optional         :: bg
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n, n1l,n2l
      logical :: bg_
      bg_=.false.
      if(present(bg)) bg_ = bg
      
      n1l = n1
      n2l = n2

      if(bg_)then
         call tensor_status_quit("ERROR(2d alloc): not implemented with bg bufffer",3434)
      else
         call tensor_allocate_tensor_dp_2d_basic(p2,n1l,n2l,tensor_mem_idx_tensor_dp,stat=stat)
      endif

   end subroutine tensor_allocate_tensor_dp_2d_sl
   subroutine tensor_allocate_tensor_dp_2d_ls(p2,n1,n2,bg,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer(kind=tensor_standard_int), intent(in)   :: n2
      logical, intent(in),  optional         :: bg
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n,n1l,n2l
      logical :: bg_
      bg_=.false.
      if(present(bg)) bg_ = bg
      
      n1l = n1
      n2l = n2
      
      if(bg_)then
         call tensor_status_quit("ERROR(2d alloc): not implemented with bg bufffer",3434)
      else
         call tensor_allocate_tensor_dp_2d_basic(p2,n1l,n2l,tensor_mem_idx_tensor_dp,stat=stat)
      endif

   end subroutine tensor_allocate_tensor_dp_2d_ls
   subroutine tensor_allocate_tensor_dp_2d_ss(p2,n1,n2,bg,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer(kind=tensor_standard_int), intent(in)   :: n2
      logical, intent(in),  optional         :: bg
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n,n1l,n2l
      logical :: bg_
      bg_=.false.
      if(present(bg)) bg_ = bg
      
      n1l = n1
      n2l = n2
      
      if(bg_)then
         call tensor_status_quit("ERROR(2d alloc): not implemented with bg bufffer",3434)
      else
         call tensor_allocate_tensor_dp_2d_basic(p2,n1l,n2l,tensor_mem_idx_tensor_dp,stat=stat)
      endif

   end subroutine tensor_allocate_tensor_dp_2d_ss

   !DEALLOCATION
   subroutine tensor_free_tensor_dp_1d(p,bg,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p(:)
      logical, intent(in),  optional         :: bg
      integer, intent(out), optional         :: stat
      logical :: bg_
      bg_=.false.
      if(present(bg)) bg_ = bg
      
      if(bg_)then
         call tensor_free_tensor_dp_basic_bg(p,tensor_mem_idx_tensor_dp,buf_tensor_dp,stat=stat)
      else
         call tensor_free_tensor_dp_basic(p,tensor_mem_idx_tensor_dp,stat=stat)
      endif

   end subroutine tensor_free_tensor_dp_1d
   subroutine tensor_free_tensor_dp_1d_mpi(p,c,bg,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p(:)
      type(c_ptr), intent(inout)             :: c
      logical, intent(in),  optional         :: bg
      integer, intent(out), optional         :: stat
      logical :: bg_
      bg_=.false.
      if(present(bg)) bg_ = bg
      
      if(bg_)then
         !if( .not. c_associated( c, c_loc(p(1))) )then
         !   call tensor_status_quit("ERROR(tensor_free_tensor_dp_1d_mpi): invalid c/p combination",23)
         !endif
         c = c_null_ptr
         call tensor_free_tensor_dp_basic_bg(p,tensor_mem_idx_tensor_dp_mpi,buf_tensor_dp,stat=stat)
      else
         call tensor_free_tensor_dp_basic_mpi(p,c,tensor_mem_idx_tensor_dp_mpi,stat=stat)
      endif

   end subroutine tensor_free_tensor_dp_1d_mpi

   subroutine tensor_free_tensor_dp_2d(p2,bg,stat)
      implicit none
#ifdef VAR_PTR_RESHAPE
      real(tensor_dp),pointer, contiguous, intent(inout) :: p2(:,:)
#else
      real(tensor_dp),pointer, intent(inout) :: p2(:,:)
#endif
      logical, intent(in),  optional         :: bg
      integer, intent(out), optional         :: stat

      call tensor_free_tensor_dp_2d_basic(p2,tensor_mem_idx_tensor_dp,stat=stat)

   end subroutine tensor_free_tensor_dp_2d



   !BASIC ALLOCATOR NORMAL
   subroutine tensor_allocate_tensor_dp_2d_basic(p,n1,n2,idx,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p(:,:)

      include "standard_2d_allocation.inc"

   end subroutine tensor_allocate_tensor_dp_2d_basic
   subroutine tensor_allocate_tensor_dp_basic(p,n,idx,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p(:)

      include "standard_allocation.inc"

   end subroutine tensor_allocate_tensor_dp_basic
   !BASIC ALLOCATOR BACKGROUND BUFFER
   subroutine tensor_allocate_tensor_dp_basic_bg(p,n,idx,buf,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout)        :: p(:)
      integer(kind=tensor_standard_int), intent(in) :: idx
      integer(kind=tensor_long_int), intent(in)     :: n
      type(tensor_bg_buf_dp_type), intent(inout)    :: buf
      integer, intent(out), optional                :: stat

#ifdef TENSORS_IN_LSDALTON
      include "bg_allocation.inc"
#else
      call tensor_status_quit("ERROR(tensor_allocate_tensor_dp_basic_bg):&
      & currently only available for LSDALTON",-1)
#endif

   end subroutine tensor_allocate_tensor_dp_basic_bg
   !BASIC ALLOCATOR MPI
   subroutine tensor_allocate_tensor_dp_basic_mpi(p,c,n,idx,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p(:)
      include "mpi_alloc_input.inc"

#ifdef VAR_MPI
      include "mpi_allocation.inc"
#else
      call tensor_status_quit("ERROR(tensor_free_tensor_dp_basic_mpi): called without mpi",223)
#endif

   end subroutine tensor_allocate_tensor_dp_basic_mpi

   !BASIC DEALLOCATOR
   subroutine tensor_free_tensor_dp_2d_basic(p,idx,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p(:,:)

      include "standard_deallocation.inc"

   end subroutine tensor_free_tensor_dp_2d_basic
   subroutine tensor_free_tensor_dp_basic(p,idx,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p(:)

      include "standard_deallocation.inc"

   end subroutine tensor_free_tensor_dp_basic
   !BASIC DEALLOCATOR BACKGROUND BUFFER
   subroutine tensor_free_tensor_dp_basic_bg(p,idx,buf,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p(:)
      integer(kind=tensor_standard_int), intent(in) :: idx
      type(tensor_bg_buf_dp_type), intent(inout)    :: buf
      integer, intent(out), optional                :: stat

#ifdef TENSORS_IN_LSDALTON
      include "bg_deallocation.inc"
#else
      call tensor_status_quit("ERROR(tensor_free_tensor_dp_basic_bg):&
      & currently only available for LSDALTON",-1)
#endif

   end subroutine tensor_free_tensor_dp_basic_bg
   !BASIC DEALLOCATOR MPI
   subroutine tensor_free_tensor_dp_basic_mpi(p,c,idx,stat)
      implicit none
      real(tensor_dp),pointer, intent(inout) :: p(:)
      include "mpi_dealloc_input.inc"
#ifdef VAR_MPI
      include "mpi_deallocation.inc"
#else
      call tensor_status_quit("ERROR(tensor_free_tensor_dp_basic_mpi): called without mpi",223)
#endif

   end subroutine tensor_free_tensor_dp_basic_mpi


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!         ATOMIC TYPE: 64bit INTEGER: tensor_long_int             !!!!!!!!!!!!!!!!!!!!!!!!!!                                      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !ALLOCATION INTERFACES
   subroutine tensor_allocate_tensor_long_int_1d_std(p,n1,stat)
      implicit none
      integer(kind=tensor_long_int),pointer, intent(inout) :: p(:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tensor_long_int_basic(p,n,tensor_mem_idx_tensor_long_int,stat=stat)

   end subroutine tensor_allocate_tensor_long_int_1d_std
   subroutine tensor_allocate_tensor_long_int_1d(p,n1,stat)
      implicit none
      integer(kind=tensor_long_int),pointer, intent(inout) :: p(:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tensor_long_int_basic(p,n,tensor_mem_idx_tensor_long_int,stat=stat)

   end subroutine tensor_allocate_tensor_long_int_1d
   subroutine tensor_allocate_tensor_long_int_2d_ll(p2,n1,n2,stat)
      implicit none
      integer(kind=tensor_long_int),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer(kind=tensor_long_int), intent(in)   :: n2
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n,n1l,n2l
      integer(kind=tensor_long_int),pointer :: p(:)
      
      n1l = n1
      n2l = n2
      
      n = n1l*n2l
      call tensor_allocate_tensor_long_int_basic(p,n,tensor_mem_idx_tensor_long_int,stat=stat)
#ifdef VAR_PTR_RESHAPE
      p2(1:n1l,1:n2l) => p(1:n)
#else
      call c_f_pointer(c_loc(p(1)),p2,[n1l,n2l])
#endif
      p => null()

   end subroutine tensor_allocate_tensor_long_int_2d_ll
   subroutine tensor_allocate_tensor_long_int_2d_sl(p2,n1,n2,stat)
      implicit none
      integer(kind=tensor_long_int),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer(kind=tensor_long_int), intent(in)   :: n2
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n,n1l,n2l
      integer(kind=tensor_long_int),pointer :: p(:)
      
      n1l = n1
      n2l = n2
      
      n = n1l*n2l
      call tensor_allocate_tensor_long_int_basic(p,n,tensor_mem_idx_tensor_long_int,stat=stat)
#ifdef VAR_PTR_RESHAPE
      p2(1:n1l,1:n2l) => p(1:n)
#else
      call c_f_pointer(c_loc(p(1)),p2,[n1l,n2l])
#endif
      p => null()

   end subroutine tensor_allocate_tensor_long_int_2d_sl
   subroutine tensor_allocate_tensor_long_int_2d_ls(p2,n1,n2,stat)
      implicit none
      integer(kind=tensor_long_int),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer(kind=tensor_standard_int), intent(in)   :: n2
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n,n1l,n2l
      integer(kind=tensor_long_int),pointer :: p(:)
      
      n1l = n1
      n2l = n2
      
      n = n1l*n2l
      call tensor_allocate_tensor_long_int_basic(p,n,tensor_mem_idx_tensor_long_int,stat=stat)
#ifdef VAR_PTR_RESHAPE
      p2(1:n1l,1:n2l) => p(1:n)
#else
      call c_f_pointer(c_loc(p(1)),p2,[n1l,n2l])
#endif
      p => null()

   end subroutine tensor_allocate_tensor_long_int_2d_ls
   subroutine tensor_allocate_tensor_long_int_2d_ss(p2,n1,n2,stat)
      implicit none
      integer(kind=tensor_long_int),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer(kind=tensor_standard_int), intent(in)   :: n2
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n,n1l,n2l
      integer(kind=tensor_long_int),pointer :: p(:)
      
      n1l = n1
      n2l = n2
      
      n = n1l*n2l
      call tensor_allocate_tensor_long_int_basic(p,n,tensor_mem_idx_tensor_long_int,stat=stat)
#ifdef VAR_PTR_RESHAPE
      p2(1:n1l,1:n2l) => p(1:n)
#else
      call c_f_pointer(c_loc(p(1)),p2,[n1l,n2l])
#endif
      p => null()

   end subroutine tensor_allocate_tensor_long_int_2d_ss


   !DEALLOCATION
   subroutine tensor_free_tensor_long_int_1d(p,stat)
      implicit none
      integer(kind=tensor_long_int),pointer, intent(inout) :: p(:)
      integer, intent(out), optional         :: stat

      call tensor_free_tensor_long_int_basic(p,tensor_mem_idx_tensor_long_int,stat=stat)

   end subroutine tensor_free_tensor_long_int_1d
   subroutine tensor_free_tensor_long_int_2d(p2,stat)
      implicit none
#ifdef VAR_PTR_RESHAPE
      integer(kind=tensor_long_int),pointer,contiguous, intent(inout) :: p2(:,:)
#else
      integer(kind=tensor_long_int),pointer, intent(inout) :: p2(:,:)
#endif
      integer(kind=tensor_long_int),pointer  :: p(:)
      integer(kind=tensor_standard_int)       :: idx  = tensor_mem_idx_tensor_long_int
      procedure(tensor_aif_long_int), pointer :: free

      include "standard_dealloc2d_var.inc"

      free => tensor_free_tensor_long_int_basic

#ifdef VAR_PTR_RESHAPE
      include "standard_dealloc2d.inc"
#else
      include "standard_dealloc2d_cptr.inc"
#endif

   end subroutine tensor_free_tensor_long_int_2d

   !BASIC ALLOCATOR
   subroutine tensor_allocate_tensor_long_int_basic(p,n,idx,stat)
      implicit none
      integer(kind=tensor_long_int),pointer, intent(inout) :: p(:)

      include "standard_allocation.inc"

   end subroutine tensor_allocate_tensor_long_int_basic
   !BASIC DEALLOCATOR
   subroutine tensor_free_tensor_long_int_basic(p,idx,stat)
      implicit none
      integer(kind=tensor_long_int),pointer, intent(inout) :: p(:)

      include "standard_deallocation.inc"

   end subroutine tensor_free_tensor_long_int_basic



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!      ATOMIC TYPE: 32bit INTEGER: tensor_standard_int          !!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !ALLOCATION INTERFACES
   subroutine tensor_allocate_tensor_standard_int_1d_std(p,n1,stat)
      implicit none
      integer(kind=tensor_standard_int),pointer, intent(inout) :: p(:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tensor_standard_int_basic(p,n,tensor_mem_idx_tensor_standard_int,stat=stat)

   end subroutine tensor_allocate_tensor_standard_int_1d_std
   subroutine tensor_allocate_tensor_standard_int_1d(p,n1,stat)
      implicit none
      integer(kind=tensor_standard_int),pointer, intent(inout) :: p(:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tensor_standard_int_basic(p,n,tensor_mem_idx_tensor_standard_int,stat=stat)

   end subroutine tensor_allocate_tensor_standard_int_1d
   subroutine tensor_allocate_tensor_standard_int_2d_ll(p2,n1,n2,stat)
      implicit none
      integer(kind=tensor_standard_int),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer(kind=tensor_long_int), intent(in)   :: n2
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n,n1l,n2l
      integer(kind=tensor_standard_int),pointer :: p(:)
      n1l = n1
      n2l = n2
      
      n = n1l*n2l
      call tensor_allocate_tensor_standard_int_basic(p,n,tensor_mem_idx_tensor_standard_int,stat=stat)
#ifdef VAR_PTR_RESHAPE
      p2(1:n1l,1:n2l) => p(1:n)
#else
      call c_f_pointer(c_loc(p(1)),p2,[n1l,n2l])
#endif
      p => null()

   end subroutine tensor_allocate_tensor_standard_int_2d_ll
   subroutine tensor_allocate_tensor_standard_int_2d_ls(p2,n1,n2,stat)
      implicit none
      integer(kind=tensor_standard_int),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer(kind=tensor_standard_int), intent(in)   :: n2
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n,n1l,n2l
      integer(kind=tensor_standard_int),pointer :: p(:)
      
      n1l = n1
      n2l = n2
      
      n = n1l*n2l
      call tensor_allocate_tensor_standard_int_basic(p,n,tensor_mem_idx_tensor_standard_int,stat=stat)
#ifdef VAR_PTR_RESHAPE
      p2(1:n1l,1:n2l) => p(1:n)
#else
      call c_f_pointer(c_loc(p(1)),p2,[n1l,n2l])
#endif
      p => null()

   end subroutine tensor_allocate_tensor_standard_int_2d_ls
   subroutine tensor_allocate_tensor_standard_int_2d_sl(p2,n1,n2,stat)
      implicit none
      integer(kind=tensor_standard_int),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer(kind=tensor_long_int), intent(in)   :: n2
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n,n1l,n2l
      integer(kind=tensor_standard_int),pointer :: p(:)
      
      n1l = n1
      n2l = n2
      
      n = n1l*n2l
      call tensor_allocate_tensor_standard_int_basic(p,n,tensor_mem_idx_tensor_standard_int,stat=stat)
#ifdef VAR_PTR_RESHAPE
      p2(1:n1l,1:n2l) => p(1:n)
#else
      call c_f_pointer(c_loc(p(1)),p2,[n1l,n2l])
#endif
      p => null()

   end subroutine tensor_allocate_tensor_standard_int_2d_sl
   subroutine tensor_allocate_tensor_standard_int_2d_ss(p2,n1,n2,stat)
      implicit none
      integer(kind=tensor_standard_int),pointer, intent(inout) :: p2(:,:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer(kind=tensor_standard_int), intent(in)   :: n2
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n, n1l, n2l
      integer(kind=tensor_standard_int),pointer :: p(:)
     
      n1l = n1
      n2l = n2
      
      n = n1l*n2l
      call tensor_allocate_tensor_standard_int_basic(p,n,tensor_mem_idx_tensor_standard_int,stat=stat)
#ifdef VAR_PTR_RESHAPE
      p2(1:n1l,1:n2l) => p(1:n)
#else
      call c_f_pointer(c_loc(p(1)),p2,[n1l,n2l])
#endif
      p => null()

   end subroutine tensor_allocate_tensor_standard_int_2d_ss


   !DEALLOCATION
   subroutine tensor_free_tensor_standard_int_1d(p,stat)
      implicit none
      integer(kind=tensor_standard_int),pointer, intent(inout) :: p(:)
      integer, intent(out), optional         :: stat

      call tensor_free_tensor_standard_int_basic(p,tensor_mem_idx_tensor_standard_int,stat=stat)

   end subroutine tensor_free_tensor_standard_int_1d
   subroutine tensor_free_tensor_standard_int_2d(p2,stat)
      implicit none
#ifdef VAR_PTR_RESHAPE
      integer(kind=tensor_standard_int), pointer, intent(inout), contiguous :: p2(:,:)
#else
      integer(kind=tensor_standard_int), pointer, intent(inout) :: p2(:,:)
#endif
      integer(kind=tensor_standard_int),pointer :: p(:)
      integer(kind=tensor_standard_int)           :: idx  = tensor_mem_idx_tensor_standard_int
      procedure(tensor_aif_standard_int), pointer :: free 
      include "standard_dealloc2d_var.inc"
      free => tensor_free_tensor_standard_int_basic
#ifdef VAR_PTR_RESHAPE
      include "standard_dealloc2d.inc"
#else
      include "standard_dealloc2d_cptr.inc"
#endif

   end subroutine tensor_free_tensor_standard_int_2d

   !BASIC ALLOCATOR
   subroutine tensor_allocate_tensor_standard_int_basic(p,n,idx,stat)
      implicit none
      integer(kind=tensor_standard_int),pointer, intent(inout) :: p(:)

      include "standard_allocation.inc"

   end subroutine tensor_allocate_tensor_standard_int_basic
   !BASIC DEALLOCATOR
   subroutine tensor_free_tensor_standard_int_basic(p,idx,stat)
      implicit none
      integer(kind=tensor_standard_int),pointer, intent(inout) :: p(:)

      include "standard_deallocation.inc"

   end subroutine tensor_free_tensor_standard_int_basic


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!         ATOMIC TYPE: 64bit LOGICAL: tensor_long_log             !!!!!!!!!!!!!!!!!!!!!!!!!!                                      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !ALLOCATION INTERFACES
   subroutine tensor_allocate_tensor_long_log_1d_std(p,n1,stat)
      implicit none
      logical(kind=tensor_long_log),pointer, intent(inout) :: p(:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tensor_long_log_basic(p,n,tensor_mem_idx_tensor_long_log,stat=stat)

   end subroutine tensor_allocate_tensor_long_log_1d_std
   subroutine tensor_allocate_tensor_long_log_1d(p,n1,stat)
      implicit none
      logical(kind=tensor_long_log),pointer, intent(inout) :: p(:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tensor_long_log_basic(p,n,tensor_mem_idx_tensor_long_log,stat=stat)

   end subroutine tensor_allocate_tensor_long_log_1d


   !DEALLOCATION
   subroutine tensor_free_tensor_long_log_1d(p,stat)
      implicit none
      logical(kind=tensor_long_log),pointer, intent(inout) :: p(:)
      integer, intent(out), optional         :: stat

      call tensor_free_tensor_long_log_basic(p,tensor_mem_idx_tensor_long_log,stat=stat)

   end subroutine tensor_free_tensor_long_log_1d

   !BASIC ALLOCATOR
   subroutine tensor_allocate_tensor_long_log_basic(p,n,idx,stat)
      implicit none
      logical(kind=tensor_long_log),pointer, intent(inout) :: p(:)

      include "standard_allocation.inc"

   end subroutine tensor_allocate_tensor_long_log_basic
   !BASIC DEALLOCATOR
   subroutine tensor_free_tensor_long_log_basic(p,idx,stat)
      implicit none
      logical(kind=tensor_long_log),pointer, intent(inout) :: p(:)

      include "standard_deallocation.inc"

   end subroutine tensor_free_tensor_long_log_basic



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!      ATOMIC TYPE: 32bit LOGICAL: tensor_standard_log          !!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !ALLOCATION INTERFACES
   subroutine tensor_allocate_tensor_standard_log_1d_std(p,n1,stat)
      implicit none
      logical(kind=tensor_standard_log),pointer, intent(inout) :: p(:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tensor_standard_log_basic(p,n,tensor_mem_idx_tensor_standard_log,stat=stat)

   end subroutine tensor_allocate_tensor_standard_log_1d_std
   subroutine tensor_allocate_tensor_standard_log_1d(p,n1,stat)
      implicit none
      logical(kind=tensor_standard_log),pointer, intent(inout) :: p(:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tensor_standard_log_basic(p,n,tensor_mem_idx_tensor_standard_log,stat=stat)

   end subroutine tensor_allocate_tensor_standard_log_1d


   !DEALLOCATION
   subroutine tensor_free_tensor_standard_log_1d(p,stat)
      implicit none
      logical(kind=tensor_standard_log),pointer, intent(inout) :: p(:)
      integer, intent(out), optional         :: stat

      call tensor_free_tensor_standard_log_basic(p,tensor_mem_idx_tensor_standard_log,stat=stat)

   end subroutine tensor_free_tensor_standard_log_1d

   !BASIC ALLOCATOR
   subroutine tensor_allocate_tensor_standard_log_basic(p,n,idx,stat)
      implicit none
      logical(kind=tensor_standard_log),pointer, intent(inout) :: p(:)

      include "standard_allocation.inc"

   end subroutine tensor_allocate_tensor_standard_log_basic
   !BASIC DEALLOCATOR
   subroutine tensor_free_tensor_standard_log_basic(p,idx,stat)
      implicit none
      logical(kind=tensor_standard_log),pointer, intent(inout) :: p(:)

      include "standard_deallocation.inc"

   end subroutine tensor_free_tensor_standard_log_basic

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!          ATOMIC TYPE: character: tensor_character             !!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !ALLOCATION INTERFACES
   subroutine tensor_allocate_character_1d_std(p,n1,stat)
      implicit none
      character, pointer, intent(inout)             :: p(:)
      integer(kind=tensor_standard_int), intent(in) :: n1
      integer, intent(out), optional                :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_character_basic(p,n,tensor_mem_idx_character,stat=stat)

   end subroutine tensor_allocate_character_1d_std
   subroutine tensor_allocate_character_1d(p,n1,stat)
      implicit none
      character,pointer, intent(inout)          :: p(:)
      integer(kind=tensor_long_int), intent(in) :: n1
      integer, intent(out), optional            :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_character_basic(p,n,tensor_mem_idx_character,stat=stat)

   end subroutine tensor_allocate_character_1d


   !DEALLOCATION
   subroutine tensor_free_character_1d(p,stat)
      implicit none
      character,pointer, intent(inout) :: p(:)
      integer, intent(out), optional   :: stat

      call tensor_free_character_basic(p,tensor_mem_idx_character,stat=stat)

   end subroutine tensor_free_character_1d

   !BASIC ALLOCATOR
   subroutine tensor_allocate_character_basic(p,n,idx,stat)
      implicit none
      character,pointer, intent(inout) :: p(:)

      include "standard_allocation.inc"

   end subroutine tensor_allocate_character_basic
   !BASIC DEALLOCATOR
   subroutine tensor_free_character_basic(p,idx,stat)
      implicit none
      character,pointer, intent(inout) :: p(:)

      include "standard_deallocation.inc"

   end subroutine tensor_free_character_basic

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             DERIVED TYPE: tile               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !ALLOCATION INTERFACES
   subroutine tensor_allocate_tile_1d_std(p,n1,stat)
      implicit none
      type(tile), pointer, intent(inout) :: p(:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tile_basic(p,n,tensor_mem_idx_tile,stat=stat)

   end subroutine tensor_allocate_tile_1d_std
   subroutine tensor_allocate_tile_1d(p,n1,stat)
      implicit none
      type(tile),pointer, intent(inout) :: p(:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tile_basic(p,n,tensor_mem_idx_tile,stat=stat)

   end subroutine tensor_allocate_tile_1d


   !DEALLOCATION
   subroutine tensor_free_tile_1d(p,stat)
      implicit none
      type(tile),pointer, intent(inout) :: p(:)
      integer, intent(out), optional         :: stat

      call tensor_free_tile_basic(p,tensor_mem_idx_tile,stat=stat)

   end subroutine tensor_free_tile_1d

   !BASIC ALLOCATOR
   subroutine tensor_allocate_tile_basic(p,n,idx,stat)
      implicit none
      type(tile),pointer, intent(inout) :: p(:)

      include "standard_allocation.inc"

   end subroutine tensor_allocate_tile_basic
   !BASIC DEALLOCATOR
   subroutine tensor_free_tile_basic(p,idx,stat)
      implicit none
      type(tile),pointer, intent(inout) :: p(:)

      include "standard_deallocation.inc"

   end subroutine tensor_free_tile_basic

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             DERIVED TYPE: tensor               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !ALLOCATION INTERFACES
   subroutine tensor_allocate_tensor_1d_std(p,n1,stat)
      implicit none
      type(tensor), pointer, intent(inout) :: p(:)
      integer(kind=tensor_standard_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tensor_basic(p,n,tensor_mem_idx_tensor,stat=stat)

   end subroutine tensor_allocate_tensor_1d_std
   subroutine tensor_allocate_tensor_1d(p,n1,stat)
      implicit none
      type(tensor),pointer, intent(inout) :: p(:)
      integer(kind=tensor_long_int), intent(in)   :: n1
      integer, intent(out), optional         :: stat
      integer(kind=tensor_long_int ) :: n
      
      n = n1
      call tensor_allocate_tensor_basic(p,n,tensor_mem_idx_tensor,stat=stat)

   end subroutine tensor_allocate_tensor_1d


   !DEALLOCATION
   subroutine tensor_free_tensor_1d(p,stat)
      implicit none
      type(tensor),pointer, intent(inout) :: p(:)
      integer, intent(out), optional         :: stat

      call tensor_free_tensor_basic(p,tensor_mem_idx_tensor,stat=stat)

   end subroutine tensor_free_tensor_1d

   !BASIC ALLOCATOR
   subroutine tensor_allocate_tensor_basic(p,n,idx,stat)
      implicit none
      type(tensor),pointer, intent(inout) :: p(:)

      include "standard_allocation.inc"

   end subroutine tensor_allocate_tensor_basic
   !BASIC DEALLOCATOR
   subroutine tensor_free_tensor_basic(p,idx,stat)
      implicit none
      type(tensor),pointer, intent(inout) :: p(:)

      include "standard_deallocation.inc"

   end subroutine tensor_free_tensor_basic

end module tensor_allocator
