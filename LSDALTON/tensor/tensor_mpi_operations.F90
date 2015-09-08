module tensor_mpi_operations_module

   use tensor_error_handler
   use tensor_parameters_and_counters
   use tensor_mpi_binding_module
   use tensor_mpi_interface_module
   use tensor_allocator, only: tensor_alloc_mem, tensor_free_mem

   public :: tensor_buffer
   public :: tensor_get_rank
   public :: tensor_get_size
   private

#ifdef VAR_MPI
   type tensor_mpi_buffer_type
      logical                       :: is_initialized = .false.
      integer(kind=tensor_long_int) :: std_size       = 10000000
      integer(kind=tensor_long_int) :: all_size       = 0
      integer(kind=tensor_long_int) :: nmsg           = 0
      integer(kind=tensor_mpi_kind) :: root           = -1
      integer(kind=tensor_mpi_kind) :: comm           = -1
      integer(kind=tensor_long_int) :: offset         = 0
      integer(kind=tensor_long_int) :: nbytes         = 0
      character, pointer            :: buffer(:)      => null()
      character, pointer            :: tmp(:)         => null()
      contains
      procedure :: init_ => tensor_initialize_mpi_buffer
      procedure :: incr_ => tensor_mpi_buffer_increase_buffer
      procedure :: free_ => tensor_free_mpi_buffer
   end type tensor_mpi_buffer_type

   type(tensor_mpi_buffer_type) :: tensor_mpi_buffer
#endif
   
   interface tensor_buffer
#ifdef VAR_MPI
      module procedure tensor_buffer_dp_l, &
                     & tensor_buffer_dp_s, &
                     & tensor_buffer_dp, &
                     & tensor_buffer_long_int_l, &
                     & tensor_buffer_long_int_s, &
                     & tensor_buffer_long_int, &
                     & tensor_buffer_standard_int_l, &
                     & tensor_buffer_standard_int_s, &
                     & tensor_buffer_standard_int, &
                     & tensor_buffer_log_l, &
                     & tensor_buffer_log_s, &
                     & tensor_buffer_log, &
                     & tensor_buffer_char_l, &
                     & tensor_buffer_char_s, &
                     & tensor_buffer_char
#else
      module procedure tensor_buffer_dummy
#endif
   end interface tensor_buffer

   contains

   subroutine tensor_get_rank(rank)
      implicit none
      integer(kind=tensor_mpi_kind) :: rank
#ifdef VAR_MPI
      if( tensor_work_comm /= tensor_comm_null)then
         call tensor_get_rank_for_comm(tensor_work_comm,rank)
      else
         rank=0
      endif
#else
      rank=0
#endif
   end subroutine tensor_get_rank
   subroutine tensor_get_size(size_)
      implicit none
      integer(kind=tensor_mpi_kind) :: size_
#ifdef VAR_MPI
      if( tensor_work_comm /= tensor_comm_null)then
         call tensor_get_size_for_comm(tensor_work_comm,size_)
      else
         size_=1
      endif
#else
      size_=1
#endif
   end subroutine tensor_get_size


   subroutine tensor_buffer_dummy()
      implicit none
      call tensor_status_quit("ERROR(tensor_buffer_dummy): you should not have been able to compile this",666)
   end subroutine tensor_buffer_dummy

#ifdef VAR_MPI

   subroutine tensor_initialize_mpi_buffer(this,sze,root,comm)
      implicit none
      class(tensor_mpi_buffer_type), intent(inout) :: this
      integer(kind=tensor_long_int), intent(inout) :: sze
      integer(kind=tensor_mpi_kind), intent(in)    :: root,comm
      integer(kind=tensor_mpi_kind) :: rank

      call tensor_get_rank_for_comm(comm,rank)

      !RECIVER RECIEVE AT INITIALIZE
      if( rank /= root )then
         call tensor_mpi_bcast(sze,root,comm)
      endif

      call tensor_alloc_mem(this%buffer,sze)
      this%nbytes         = sze
      this%root           = root
      this%comm           = comm
      this%is_initialized = .true.
      this%offset         = 1

      if( rank /= root )then
         call tensor_mpi_bcast(tensor_mpi_buffer%buffer,sze,this%root,this%comm)
      endif

   end subroutine tensor_initialize_mpi_buffer
   subroutine tensor_free_mpi_buffer(this)
      implicit none
      class(tensor_mpi_buffer_type), intent(inout) :: this
      integer(kind=tensor_mpi_kind) :: rank
      integer(kind=tensor_long_int) :: sze

      call tensor_get_rank_for_comm(this%comm,rank)

      !SENDER SEND AT FINALIZE
      if( rank == this%root )then
         sze = tensor_mpi_buffer%offset - 1
         call tensor_mpi_bcast(sze,this%root,this%comm)
         call tensor_mpi_bcast(tensor_mpi_buffer%buffer,sze,this%root,this%comm)
      endif

      !calculate new standard message size
      this%nmsg     = this%nmsg + 1
      this%all_size = this%all_size + this%offset
      this%std_size = ceiling(float(this%all_size)/float(this%nmsg))

      !fee and set standard values
      call tensor_free_mem(this%buffer)
      this%nbytes         = 0
      this%root           = -1
      this%comm           = -1
      this%is_initialized = .false.
      this%offset         = 0
   end subroutine tensor_free_mpi_buffer
   subroutine tensor_mpi_buffer_increase_buffer(this,inc)
      implicit none
      class(tensor_mpi_buffer_type), intent(inout) :: this
      integer(kind=tensor_long_int), intent(in)    :: inc
      integer(kind=tensor_long_int) :: oldbytes
#ifdef VAR_LSDEBUG
      print *,"WARNING(tensor_mpi_buffer_increase_buffer): increasing MPI buffer from:",&
         &this%nbytes," to ",this%nbytes + max(this%std_size,inc)
#endif

      oldbytes     = this%nbytes

      call tensor_alloc_mem(this%tmp,this%nbytes)

      this%tmp     = this%buffer
      this%nbytes  = this%nbytes + max(this%std_size,inc)

      call tensor_free_mem(this%buffer)
      call tensor_alloc_mem(this%buffer,this%nbytes)

      this%buffer(1:oldbytes) = this%tmp

      call tensor_free_mem(this%tmp)
   end subroutine tensor_mpi_buffer_increase_buffer

   subroutine tensor_buffer_dp(buffer,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_long_int), parameter :: datatype = tensor_dp
      real(tensor_dp) :: buffer
      real(tensor_dp) :: buf(1)
      integer, intent(in), optional :: init_size
      logical, optional :: finalize

#ifdef USE_MPI_MOD_F08
      type(MPI_Comm),intent(in),optional :: comm
      integer(kind=tensor_mpi_kind),intent(in), optional :: root
#else
      integer(kind=tensor_mpi_kind),intent(in), optional :: comm
      integer(kind=tensor_mpi_kind),intent(in), optional :: root
#endif

      integer(kind=tensor_long_int) :: n
      integer(kind=tensor_mpi_kind) :: i, j, c

      if(present(comm))then
         c=comm
      else if(tensor_mpi_buffer%is_initialized)then
         c=tensor_mpi_buffer%comm
      else
         call tensor_status_quit("ERROR(tensor_buffer_dp): buffer needs to be initialized or comm given at first call",22)
      endif

      call tensor_get_rank_for_comm(c,i)
      j = tensor_mpi_buffer%root

      n=1
      buf(1) = buffer
      call tensor_buffer_dp_l(buf,n,init_size=init_size,root=root,comm=comm,finalize=finalize)
      if(i /= j) buffer = buf(1)
   end subroutine tensor_buffer_dp
   subroutine tensor_buffer_dp_l(buffer,n1,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      integer(kind=tensor_long_int), parameter :: datatype = tensor_dp
      real(tensor_dp) :: buffer(n1)
      include "tensor_buffer_body.inc"
   end subroutine tensor_buffer_dp_l
   subroutine tensor_buffer_dp_s(buffer,n1,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      integer(kind=tensor_long_int), parameter :: datatype = tensor_dp
      real(tensor_dp) :: buffer(n1)
      include "tensor_buffer_body.inc"
   end subroutine tensor_buffer_dp_s
   subroutine tensor_buffer_long_int(buffer,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_long_int), parameter :: datatype = tensor_long_int
      integer(kind=tensor_long_int) :: buffer
      integer(kind=tensor_long_int) :: buf(1)
      integer, intent(in), optional :: init_size
      logical, optional :: finalize

#ifdef USE_MPI_MOD_F08
      type(MPI_Comm),intent(in),optional :: comm
      integer(kind=tensor_mpi_kind),intent(in), optional :: root
#else
      integer(kind=tensor_mpi_kind),intent(in), optional :: comm
      integer(kind=tensor_mpi_kind),intent(in), optional :: root
#endif

      integer(kind=tensor_long_int) :: n
      integer(kind=tensor_mpi_kind) :: i, j, c

      if(present(comm))then
         c=comm
      else if(tensor_mpi_buffer%is_initialized)then
         c=tensor_mpi_buffer%comm
      else
         call tensor_status_quit("ERROR(tensor_buffer_long_int): buffer needs to be initialized or comm given at first call",22)
      endif

      call tensor_get_rank_for_comm(c,i)
      j = tensor_mpi_buffer%root

      n=1
      buf(1) = buffer
      call tensor_buffer_long_int_l(buf,n,init_size=init_size,root=root,comm=comm,finalize=finalize)
      if(i /= j) buffer = buf(1)
   end subroutine tensor_buffer_long_int
   subroutine tensor_buffer_long_int_l(buffer,n1,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      integer(kind=tensor_long_int), parameter :: datatype = tensor_long_int
      integer(kind=tensor_long_int) :: buffer(n1)
      include "tensor_buffer_body.inc"
   end subroutine tensor_buffer_long_int_l
   subroutine tensor_buffer_long_int_s(buffer,n1,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      integer(kind=tensor_long_int), parameter :: datatype = tensor_long_int
      integer(kind=tensor_long_int) :: buffer(n1)
      include "tensor_buffer_body.inc"
   end subroutine tensor_buffer_long_int_s
   subroutine tensor_buffer_standard_int(buffer,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_long_int), parameter :: datatype = tensor_standard_int
      integer(kind=tensor_standard_int) :: buffer
      integer(kind=tensor_standard_int) :: buf(1)
      integer, intent(in), optional :: init_size
      logical, optional :: finalize

#ifdef USE_MPI_MOD_F08
      type(MPI_Comm),intent(in),optional :: comm
      integer(kind=tensor_mpi_kind),intent(in), optional :: root
#else
      integer(kind=tensor_mpi_kind),intent(in), optional :: comm
      integer(kind=tensor_mpi_kind),intent(in), optional :: root
#endif

      integer(kind=tensor_long_int) :: n
      integer(kind=tensor_mpi_kind) :: i, j, c

      if(present(comm))then
         c=comm
      else if(tensor_mpi_buffer%is_initialized)then
         c=tensor_mpi_buffer%comm
      else
         call tensor_status_quit("ERROR(tensor_buffer_standard_int): buffer needs to be initialized or comm given at first call",22)
      endif

      call tensor_get_rank_for_comm(c,i)
      j = tensor_mpi_buffer%root

      n=1
      buf(1) = buffer
      call tensor_buffer_standard_int_l(buf,n,init_size=init_size,root=root,comm=comm,finalize=finalize)
      if(i/= j) buffer = buf(1)
   end subroutine tensor_buffer_standard_int
   subroutine tensor_buffer_standard_int_l(buffer,n1,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      integer(kind=tensor_long_int), parameter :: datatype = tensor_standard_int
      integer(kind=tensor_standard_int) :: buffer(n1)
      include "tensor_buffer_body.inc"
   end subroutine tensor_buffer_standard_int_l
   subroutine tensor_buffer_standard_int_s(buffer,n1,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      integer(kind=tensor_long_int), parameter :: datatype = tensor_standard_int
      integer(kind=tensor_standard_int) :: buffer(n1)
      include "tensor_buffer_body.inc"
   end subroutine tensor_buffer_standard_int_s
   subroutine tensor_buffer_log(buffer,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_long_int), parameter :: datatype = tensor_log
      logical :: buffer
      logical :: buf(1)
      integer, intent(in), optional :: init_size
      logical, optional :: finalize

#ifdef USE_MPI_MOD_F08
      type(MPI_Comm),intent(in),optional :: comm
      integer(kind=tensor_mpi_kind),intent(in), optional :: root
#else
      integer(kind=tensor_mpi_kind),intent(in), optional :: comm
      integer(kind=tensor_mpi_kind),intent(in), optional :: root
#endif

      integer(kind=tensor_long_int) :: n
      integer(kind=tensor_mpi_kind) :: i, j, c

      if(present(comm))then
         c=comm
      else if(tensor_mpi_buffer%is_initialized)then
         c=tensor_mpi_buffer%comm
      else
         call tensor_status_quit("ERROR(tensor_buffer_log): buffer needs to be initialized or comm given at first call",22)
      endif

      call tensor_get_rank_for_comm(c,i)
      j = tensor_mpi_buffer%root

      n=1
      buf(1) = buffer
      call tensor_buffer_log_l(buf,n,init_size=init_size,root=root,comm=comm,finalize=finalize)
      if(i/= j) buffer = buf(1)
   end subroutine tensor_buffer_log
   subroutine tensor_buffer_log_l(buffer,n1,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      integer(kind=tensor_long_int), parameter :: datatype = tensor_log
      logical :: buffer(n1)
      include "tensor_buffer_body.inc"
   end subroutine tensor_buffer_log_l
   subroutine tensor_buffer_log_s(buffer,n1,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      integer(kind=tensor_long_int), parameter :: datatype = tensor_log
      logical :: buffer(n1)
      include "tensor_buffer_body.inc"
   end subroutine tensor_buffer_log_s
   subroutine tensor_buffer_char(buffer,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_long_int), parameter :: datatype = tensor_char_size
      character*(*) :: buffer
      character :: buf(len(buffer))
      integer, intent(in), optional :: init_size
      logical, optional :: finalize

#ifdef USE_MPI_MOD_F08
      type(MPI_Comm),intent(in),optional :: comm
      integer(kind=tensor_mpi_kind),intent(in), optional :: root
#else
      integer(kind=tensor_mpi_kind),intent(in), optional :: comm
      integer(kind=tensor_mpi_kind),intent(in), optional :: root
#endif

      integer(kind=tensor_long_int) :: n
      integer(kind=tensor_mpi_kind) :: i, j, c

      if(present(comm))then
         c=comm
      else if(tensor_mpi_buffer%is_initialized)then
         c=tensor_mpi_buffer%comm
      else
         call tensor_status_quit("ERROR(tensor_buffer_log): buffer needs to be initialized or comm given at first call",22)
      endif

      call tensor_get_rank_for_comm(c,i)
      j = tensor_mpi_buffer%root

      n   = len(buffer,kind=tensor_long_int)
      buf = transfer(buffer,buf)
      call tensor_buffer_char_l(buf,n,init_size=init_size,root=root,comm=comm,finalize=finalize)
      if(i/= j) buffer = transfer(buf,buffer)
   end subroutine tensor_buffer_char
   subroutine tensor_buffer_char_l(buffer,n1,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      integer(kind=tensor_long_int), parameter :: datatype = tensor_char_size
      character :: buffer(n1)
      include "tensor_buffer_body.inc"
   end subroutine tensor_buffer_char_l
   subroutine tensor_buffer_char_s(buffer,n1,init_size,root,comm,finalize)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      integer(kind=tensor_long_int), parameter :: datatype = tensor_char_size
      character :: buffer(n1)
      include "tensor_buffer_body.inc"
   end subroutine tensor_buffer_char_s

!ENDIF VAR_MPI
#endif

end module tensor_mpi_operations_module
