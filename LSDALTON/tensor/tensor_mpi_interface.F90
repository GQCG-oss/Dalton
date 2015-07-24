module tensor_mpi_interface_module
#ifdef VAR_MPI
   use lsmpi_module
   use infpar_module, only: infpar
#endif

   !use tensor_mpi_binding_module
   use tensor_parameters_and_counters
   use tensor_type_def_module

   !NOTE THAT ONLY MPI_IN_PLACE HAS BEEN IMPLEMENTED FOR THE MPI COLLECTIVES
   public :: tensor_mpi_bcast
   public :: tensor_mpi_reduce
   public :: tensor_mpi_allreduce
   public :: tensor_mpi_barrier
   public :: tensor_get_rank_for_comm
   public :: tensor_get_size_for_comm
   private

   interface tensor_mpi_bcast
#ifdef VAR_MPI
      module procedure tensor_mpi_bcast_std_int_s,&
                     & tensor_mpi_bcast_std_int_l, &
                     & tensor_mpi_bcast_std_int, &
                     & tensor_mpi_bcast_long_int_s,&
                     & tensor_mpi_bcast_long_int_l, &
                     & tensor_mpi_bcast_long_int, &
                     & tensor_mpi_bcast_tensor_dp_s,&
                     & tensor_mpi_bcast_tensor_dp_l, &
                     & tensor_mpi_bcast_tensor_dp, &
                     & tensor_mpi_bcast_log_s,&
                     & tensor_mpi_bcast_log_l, &
                     & tensor_mpi_bcast_log, &
                     & tensor_mpi_bcast_char_s,&
                     & tensor_mpi_bcast_char_l, &
                     & tensor_mpi_bcast_char
#else
      module procedure tensor_mpi_dummy
#endif
   end interface tensor_mpi_bcast

   interface tensor_mpi_reduce
#ifdef VAR_MPI
      module procedure tensor_mpi_reduce_tensor_dp_s,&
                     & tensor_mpi_reduce_tensor_dp_l, &
                     & tensor_mpi_reduce_tensor_dp
#else
      module procedure tensor_mpi_dummy
#endif
   end interface tensor_mpi_reduce

   interface tensor_mpi_allreduce
#ifdef VAR_MPI
      module procedure tensor_mpi_dummy
#else
      module procedure tensor_mpi_dummy
#endif
   end interface tensor_mpi_allreduce

   contains

    subroutine tensor_mpi_dummy()
       implicit none
       call tensor_status_quit("ERROR(tensor_mpi_no_mpi): you should not have been able to compile this!!",666)
    end subroutine tensor_mpi_dummy

#ifdef VAR_MPI

    subroutine tensor_mpi_barrier(comm)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Comm),intent(in)                :: comm
#else
       integer(kind=tensor_mpi_kind),intent(in) :: comm
#endif
       integer(kind=tensor_mpi_kind) :: ierr = 0

       call mpi_barrier(comm, ierr)

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_barrier): barrier failed",220)
    end subroutine tensor_mpi_barrier

    subroutine tensor_get_rank_for_comm(comm,rank)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Comm),intent(in)                   :: comm
#else
       integer(kind=tensor_mpi_kind),intent(in)    :: comm
#endif
       integer(kind=tensor_mpi_kind),intent(inout) :: rank
       integer(kind=tensor_mpi_kind) :: ierr = 0

       call mpi_comm_rank(comm,rank,ierr)

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_get_rank_for_comm): failed",220)
    end subroutine tensor_get_rank_for_comm

    subroutine tensor_get_size_for_comm(comm,nodtot)
#ifdef USE_MPI_MOD_F08
       type(MPI_Comm),intent(in)                   :: comm
#else
       integer(kind=tensor_mpi_kind),intent(in)    :: comm
#endif
       integer(kind=tensor_mpi_kind),intent(inout) :: nodtot
       integer(kind=tensor_mpi_kind) :: ierr = 0

       call mpi_comm_size(comm,nodtot,ierr)

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_get_size_for_comm): failed",220)
    end subroutine tensor_get_size_for_comm

    !STD INTEGER BCAST
    subroutine tensor_mpi_bcast_std_int(b,root,comm)
      implicit none
      integer(kind=tensor_standard_int) :: b
      include "mpi_collective_vars.inc"
      integer(kind=tensor_standard_int) :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_std_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_standard_int))
      if( infpar%lg_mynum /= root) b = buffer(1)
    end subroutine tensor_mpi_bcast_std_int
    subroutine tensor_mpi_bcast_std_int_s(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in)    :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_std_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_standard_int))
    end subroutine tensor_mpi_bcast_std_int_s
    subroutine tensor_mpi_bcast_std_int_l(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_std_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_standard_int))
    end subroutine tensor_mpi_bcast_std_int_l

    subroutine tensor_mpi_bcast_std_int_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)        :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)        :: stats
      include "mpi_collective_vars.inc"
      include "mpi_bcast_std.inc"
    end subroutine tensor_mpi_bcast_std_int_basic

    !LONG INTEGER BCAST
    subroutine tensor_mpi_bcast_long_int(b,root,comm)
      implicit none
      integer(kind=tensor_long_int) :: b
      include "mpi_collective_vars.inc"
      integer(kind=tensor_long_int) :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_long_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_long_int))
      if( infpar%lg_mynum /= root) b = buffer(1)
    end subroutine tensor_mpi_bcast_long_int
    subroutine tensor_mpi_bcast_long_int_s(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      integer(kind=tensor_long_int), intent(inout)  :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_long_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_long_int))
    end subroutine tensor_mpi_bcast_long_int_s
    subroutine tensor_mpi_bcast_long_int_l(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      integer(kind=tensor_long_int), intent(inout) :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_long_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_long_int))
    end subroutine tensor_mpi_bcast_long_int_l

    subroutine tensor_mpi_bcast_long_int_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      integer(kind=tensor_long_int), intent(inout) :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_collective_vars.inc"
      include "mpi_bcast_std.inc"
    end subroutine tensor_mpi_bcast_long_int_basic

    !LOGICAL BCAST
    subroutine tensor_mpi_bcast_log(b,root,comm)
      implicit none
      logical, intent(inout) :: b
      include "mpi_collective_vars.inc"
      logical :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_log_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_log))
      if( infpar%lg_mynum /= root) b = buffer(1)
    end subroutine tensor_mpi_bcast_log
    subroutine tensor_mpi_bcast_log_s(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in)    :: n1
      logical, intent(inout) :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_log_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_log))
    end subroutine tensor_mpi_bcast_log_s
    subroutine tensor_mpi_bcast_log_l(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      logical, intent(inout) :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_log_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_log))
    end subroutine tensor_mpi_bcast_log_l

    subroutine tensor_mpi_bcast_log_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      logical, intent(inout)                       :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_collective_vars.inc"
      include "mpi_bcast_std.inc"
    end subroutine tensor_mpi_bcast_log_basic

    !DOUBLE PRECISION BCAST
    subroutine tensor_mpi_bcast_tensor_dp(b,root,comm)
      implicit none
      real(tensor_dp), intent(inout) :: b
      include "mpi_collective_vars.inc"
      real(tensor_dp) :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_dp_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_dp))
      if( infpar%lg_mynum /= root) b = buffer(1)
    end subroutine tensor_mpi_bcast_tensor_dp
    subroutine tensor_mpi_bcast_tensor_dp_s(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      real(tensor_dp), intent(inout)                :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_dp_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_dp))
    end subroutine tensor_mpi_bcast_tensor_dp_s
    subroutine tensor_mpi_bcast_tensor_dp_l(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      real(tensor_dp), intent(inout)            :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_dp_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_dp))
    end subroutine tensor_mpi_bcast_tensor_dp_l
    subroutine tensor_mpi_bcast_dp_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      real(tensor_dp), intent(inout)               :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_collective_vars.inc"
      include "mpi_bcast_std.inc"
    end subroutine tensor_mpi_bcast_dp_basic

    !CHARACTER BCAST
    subroutine tensor_mpi_bcast_char(buffer,root,comm)
      implicit none
      character*(*), intent(inout) :: buffer
      include "mpi_collective_vars.inc"
      n=1
      call tensor_mpi_bcast_char_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_char))
    end subroutine tensor_mpi_bcast_char
    subroutine tensor_mpi_bcast_char_s(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      character*(*), intent(inout)                  :: buffer
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_char_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_char))
    end subroutine tensor_mpi_bcast_char_s
    subroutine tensor_mpi_bcast_char_l(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      character*(*), intent(inout)              :: buffer
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_char_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_char))
    end subroutine tensor_mpi_bcast_char_l
    subroutine tensor_mpi_bcast_char_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      character*(*), intent(inout)              :: buffer
      type(tensor_mpi_stats_type),intent(inout) :: stats
      include "mpi_collective_vars.inc"
      include "mpi_bcast_std.inc"
    end subroutine tensor_mpi_bcast_char_basic


    !DOUBLE PRECISION REDUCE
    subroutine tensor_mpi_reduce_tensor_dp(b,root,comm)
      implicit none
      real(tensor_dp), intent(inout) :: b
      include "mpi_collective_vars.inc"
      real(tensor_dp) :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_reduce_dp_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_reduce,tensor_mpi_idx_tensor_dp))
      if( infpar%lg_mynum /= root) b = buffer(1)
    end subroutine tensor_mpi_reduce_tensor_dp
    subroutine tensor_mpi_reduce_tensor_dp_s(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      real(tensor_dp), intent(inout)                :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_reduce_dp_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_reduce,tensor_mpi_idx_tensor_dp))
    end subroutine tensor_mpi_reduce_tensor_dp_s
    subroutine tensor_mpi_reduce_tensor_dp_l(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      real(tensor_dp), intent(inout)            :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_reduce_dp_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_reduce,tensor_mpi_idx_tensor_dp))
    end subroutine tensor_mpi_reduce_tensor_dp_l
    subroutine tensor_mpi_reduce_dp_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      real(tensor_dp), intent(inout)               :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_collective_vars.inc"
      include "mpi_reduce_std.inc"
    end subroutine tensor_mpi_reduce_dp_basic
#endif

end module tensor_mpi_interface_module
