module tensor_mpi_interface_module

   use tensor_mpi_binding_module
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
   public :: tensor_mpi_wait
   public :: tensor_mpi_win_lock
   public :: tensor_mpi_win_unlock
   public :: tensor_mpi_win_lock_all
   public :: tensor_mpi_win_unlock_all
   public :: tensor_mpi_win_flush
   public :: tensor_mpi_probe
   public :: tensor_mpi_get_count
   public :: tensor_mpi_sendrecv
   public :: tensor_mpi_win_create
   public :: tensor_mpi_win_free
   public :: tensor_mpi_put
   public :: tensor_mpi_get
   public :: tensor_mpi_acc
   public :: tensor_mpi_win_fence

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
                     & tensor_mpi_bcast_char, &
                     & tensor_mpi_bcast_char_vec_l, &
                     & tensor_mpi_bcast_char_vec_s
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
      module procedure tensor_mpi_allreduce_tensor_dp_s,&
                     & tensor_mpi_allreduce_tensor_dp_l, &
                     & tensor_mpi_allreduce_tensor_dp, &
                     & tensor_mpi_allreduce_long_int_s,&
                     & tensor_mpi_allreduce_long_int_l, &
                     & tensor_mpi_allreduce_long_int, &
                     & tensor_mpi_allreduce_std_int_s,&
                     & tensor_mpi_allreduce_std_int_l, &
                     & tensor_mpi_allreduce_std_int
#else
      module procedure tensor_mpi_dummy
#endif
   end interface tensor_mpi_allreduce

   interface tensor_mpi_sendrecv
#ifdef VAR_MPI
      module procedure tensor_mpi_sendrecv_tensor_dp_s, &
                     & tensor_mpi_sendrecv_tensor_dp_l, &
                     & tensor_mpi_sendrecv_tensor_dp, &
                     & tensor_mpi_sendrecv_tensor_standard_int_s, &
                     & tensor_mpi_sendrecv_tensor_standard_int_l, &
                     & tensor_mpi_sendrecv_tensor_standard_int, &
                     & tensor_mpi_sendrecv_tensor_long_int_s, &
                     & tensor_mpi_sendrecv_tensor_long_int_l, &
                     & tensor_mpi_sendrecv_tensor_long_int
#else
      module procedure tensor_mpi_dummy
#endif
   end interface tensor_mpi_sendrecv

   interface tensor_mpi_win_create
#ifdef VAR_MPI
      module procedure tensor_mpi_win_create_dp_s,&
                     & tensor_mpi_win_create_dp_l
#else
      module procedure tensor_mpi_dummy
#endif
   end interface tensor_mpi_win_create

   interface tensor_mpi_win_fence
#ifdef VAR_MPI
      module procedure tensor_mpi_win_fence_simple, tensor_mpi_win_fence_special
#else
      module procedure tensor_mpi_dummy
#endif
   end interface tensor_mpi_win_fence

   interface tensor_mpi_put
#ifdef VAR_MPI
      module procedure tensor_mpi_put_dp_s, &
                     & tensor_mpi_put_dp_l, &
                     & tensor_mpi_rput_dp_s, &
                     & tensor_mpi_rput_dp_l
#else
      module procedure tensor_mpi_dummy
#endif
   end interface tensor_mpi_put

   interface tensor_mpi_get
#ifdef VAR_MPI
      module procedure tensor_mpi_get_dp_s, &
                     & tensor_mpi_get_dp_l, &
                     & tensor_mpi_rget_dp_s, &
                     & tensor_mpi_rget_dp_l
#else
      module procedure tensor_mpi_dummy
#endif
   end interface tensor_mpi_get

   interface tensor_mpi_acc
#ifdef VAR_MPI
      module procedure tensor_mpi_acc_dp_s, &
                     & tensor_mpi_acc_dp_l, &
                     & tensor_mpi_racc_dp_s, &
                     & tensor_mpi_racc_dp_l
#else
      module procedure tensor_mpi_dummy
#endif
   end interface tensor_mpi_acc

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

    subroutine tensor_mpi_wait(request_handle)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Request),intent(inout) :: request_handle
#else
       integer(kind=tensor_mpi_kind),intent(inout) :: request_handle
#endif
       integer(kind=tensor_mpi_kind) :: ierr = 0

       call mpi_wait(request_handle,MPI_STATUS_IGNORE,ierr)

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_wait): failed",220)
    end subroutine tensor_mpi_wait

    subroutine tensor_mpi_win_lock(dest,win,typeoflock,ass)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Win),intent(in) :: win
#else
       integer(kind=tensor_mpi_kind),intent(in) :: win
#endif
       integer(kind=tensor_mpi_kind),intent(in) :: dest
       integer(kind=tensor_mpi_kind),intent(in),optional :: ass
       character, intent(in) :: typeoflock
       integer(kind=tensor_mpi_kind) :: assert = 0
       integer(kind=tensor_mpi_kind) :: ierr   = 0

       if(present(ass))assert=ass
       if(typeoflock=='e')then
          CALL mpi_win_lock(MPI_LOCK_EXCLUSIVE,dest,assert,win,ierr)
       else if(typeoflock=='s')then
          CALL mpi_win_lock(MPI_LOCK_SHARED,dest,assert,win,ierr)
       else
          call tensor_status_quit("ERROR(tensor_mpi_win_lock): no valid lock type selected",636)
       endif

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_win_lock): failed",220)
    end subroutine tensor_mpi_win_lock

    subroutine tensor_mpi_win_unlock(dest,win)
       implicit none
       integer(kind=tensor_mpi_kind),intent(in) :: dest
#ifdef USE_MPI_MOD_F08
       type(MPI_Win), intent(in) :: win
#else
       integer(kind=tensor_mpi_kind), intent(in) :: win
#endif
       integer(kind=tensor_mpi_kind) :: ierr = 0

       call mpi_win_unlock(dest,win,ierr)     

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_win_unlock): failed",220)
    end subroutine tensor_mpi_win_unlock

    subroutine tensor_mpi_win_flush(win,rank,local)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Win), intent(in) :: win
#else
       integer(kind=tensor_mpi_kind), intent(in) :: win
#endif
       integer(kind=tensor_mpi_kind), optional :: rank
       logical, optional :: local
       integer(kind=tensor_mpi_kind) :: ierr=0
       logical :: loc
       loc  = .false.

       if(present(local))loc = local

#ifdef VAR_HAVE_MPI3
       if(loc)then
          if(present(rank))then
             call mpi_win_flush_local(rank,win,ierr)
          else
             call mpi_win_flush_local_all(win,ierr)
          endif
       else
          if(present(rank))then
             call mpi_win_flush(rank,win,ierr)
          else
             call mpi_win_flush_all(win,ierr)
          endif
       endif

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_win_flush): failed",220)
#else
       print *,"WARNING(tensor_mpi_win_flush)called without MPI3, unlock should force&
          & the sync"
#endif
    end subroutine tensor_mpi_win_flush

    subroutine tensor_mpi_win_lock_all(win,ass)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Win), intent(in) :: win
#else
       integer(kind=tensor_mpi_kind), intent(in) :: win
#endif
       integer(kind=tensor_mpi_kind),intent(in),optional :: ass
       integer(kind=tensor_mpi_kind) :: assert = 0
       integer(kind=tensor_mpi_kind) :: ierr   = 0

       if(present(ass))assert=ass

#ifdef VAR_HAVE_MPI3
       CALL mpi_win_lock_all(assert,win,ierr)
#else
       call tensor_status_quit("ERROR(tensor_mpi_win_lock_all): this routine is MPI 3 only ",-1)
#endif

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_win_lock_all): failed",220)
    end subroutine tensor_mpi_win_lock_all

    subroutine tensor_mpi_win_unlock_all(win)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Win), intent(in) :: win
#else
       integer(kind=tensor_mpi_kind), intent(in) :: win
#endif
       integer(kind=tensor_mpi_kind) :: ierr   = 0

#ifdef VAR_HAVE_MPI3
       CALL mpi_win_unlock_all(win,ierr)
#else
       call tensor_status_quit("ERROR(tensor_mpi_win_unlock_all): this routine is MPI 3 only ",636)
#endif

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_win_unlock_all): failed",220)
    end subroutine tensor_mpi_win_unlock_all

    subroutine tensor_mpi_win_free(win)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Win), intent(inout) :: win
#else
       integer(kind=tensor_mpi_kind), intent(inout) :: win
#endif
       integer(kind=tensor_mpi_kind) :: ierr = 0

       call mpi_win_free(win,ierr)

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_win_free): failed",220)
    end subroutine tensor_mpi_win_free

    subroutine tensor_mpi_win_fence_simple(win,assert)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Win), intent(inout) :: win
#else
       integer(kind=tensor_mpi_kind), intent(inout) :: win
#endif
       integer(kind=tensor_mpi_kind),intent(in),optional :: assert
       integer(kind=tensor_mpi_kind) :: as   = 0
       integer(kind=tensor_mpi_kind) :: ierr = 0

       if(present(assert)) as = assert

       call mpi_win_fence(as,win,ierr)

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_win_fence): failed",220)
    end subroutine tensor_mpi_win_fence_simple

    subroutine tensor_mpi_win_fence_special(win,openwin)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Win), intent(inout) :: win
#else
       integer(kind=tensor_mpi_kind), intent(inout) :: win
#endif
       logical,intent(in) :: openwin
       integer(kind=tensor_mpi_kind) :: ierr = 0
       if(openwin)then
          call mpi_win_fence(MPI_MODE_NOPRECEDE,win,ierr)
          if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_win_fence_noprecede): failed",220)
       else
          call mpi_win_fence(MPI_MODE_NOSUCCEED,win,ierr)
          if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_win_fence_noprecede): failed",220)
       endif
    end subroutine tensor_mpi_win_fence_special

    subroutine tensor_mpi_win_create_dp_s(darr,win,n1,comm)
       implicit none
       integer(kind=tensor_standard_int), intent(in) :: n1
       real(tensor_dp),intent(in) :: darr(n1)
#ifdef USE_MPI_MOD_F08
       include "mpi_win_create_vars_mpif08.inc"
#else
       include "mpi_win_create_vars_std.inc"
#endif
       datatype = MPI_DOUBLE_PRECISION
       n        = n1
       call tensor_mpi_win_create_dp_basic(darr,win,n,comm)
    end subroutine tensor_mpi_win_create_dp_s
    subroutine tensor_mpi_win_create_dp_l(darr,win,n1,comm)
       implicit none
       integer(kind=tensor_long_int), intent(in) :: n1
       real(tensor_dp),intent(in) :: darr(n1)
#ifdef USE_MPI_MOD_F08
       include "mpi_win_create_vars_mpif08.inc"
#else
       include "mpi_win_create_vars_std.inc"
#endif
       datatype = MPI_DOUBLE_PRECISION
       n        = n1
       call tensor_mpi_win_create_dp_basic(darr,win,n,comm)
    end subroutine tensor_mpi_win_create_dp_l
    subroutine tensor_mpi_win_create_dp_basic(darr,win,n1,comm)
       implicit none
       integer(kind=tensor_long_int), intent(in) :: n1
       real(tensor_dp),intent(in) :: darr(n1)
#ifdef USE_MPI_MOD_F08
       include "mpi_win_create_vars_mpif08.inc"
#else
       include "mpi_win_create_vars_std.inc"
#endif
       datatype = MPI_DOUBLE_PRECISION

       include "mpi_win_create_generic.inc"
    end subroutine tensor_mpi_win_create_dp_basic

    subroutine tensor_mpi_probe(stat,comm,source,tag)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Comm),intent(in)                   :: comm
       type(MPI_Status),intent(inout)              :: stat
#else
       integer(kind=tensor_mpi_kind),intent(in)    :: comm
       integer(kind=tensor_mpi_kind),intent(inout) :: stat(MPI_STATUS_SIZE)
#endif
       integer(kind=tensor_mpi_kind),intent(in), optional :: source, tag
       integer(kind=tensor_mpi_kind) :: ierr = 0
       integer(kind=tensor_mpi_kind) :: s = MPI_ANY_SOURCE
       integer(kind=tensor_mpi_kind) :: t = MPI_ANY_TAG
       
       if(present(source))s = source
       if(present(tag))   t = tag

       call mpi_probe(s, t, comm, stat, ierr)

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_probe): probe failed",220)
    end subroutine tensor_mpi_probe

    subroutine tensor_mpi_get_count(stat,dat,nelms)
       implicit none
#ifdef USE_MPI_MOD_F08
       type(MPI_Datatype),intent(in)                :: dat
       type(MPI_Status),intent(inout)               :: stat
#else
       integer(kind=tensor_mpi_kind),intent(in)    :: dat
       integer(kind=tensor_mpi_kind),intent(inout) :: stat(MPI_STATUS_SIZE)
#endif
       integer(kind=tensor_mpi_kind),intent(out) :: nelms
       integer(kind=tensor_mpi_kind) :: ierr = 0

       call mpi_get_count(stat, dat, nelms, ierr)

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_get_count): failed",220)
    end subroutine tensor_mpi_get_count

    !STD INTEGER BCAST
    subroutine tensor_mpi_bcast_std_int(b,root,comm)
      implicit none
      integer(kind=tensor_standard_int) :: b
      integer(kind=tensor_standard_int) :: buffer(1)
      include "mpi_collective_vars.inc"
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_std_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_standard_int))
      if( rank /= root ) b = buffer(1)
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
      integer(kind=tensor_long_int) :: buffer(1)
      include "mpi_collective_vars.inc"
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_long_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_long_int))
      if( rank /= root ) b = buffer(1)
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
      logical :: buffer(1)
      include "mpi_collective_vars.inc"
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_log_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_log))
      if( rank /= root ) b = buffer(1)
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
      real(tensor_dp) :: buffer(1)
      include "mpi_collective_vars.inc"
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_dp_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_dp))
      if( rank /= root ) b = buffer(1)
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
    subroutine tensor_mpi_bcast_char_vec_s(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      character, intent(inout)                      :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_char_vec_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_char))
    end subroutine tensor_mpi_bcast_char_vec_s
    subroutine tensor_mpi_bcast_char_vec_l(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      character, intent(inout)                  :: buffer(n1)
      include "mpi_collective_vars.inc"
      n=n1
      call tensor_mpi_bcast_char_vec_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_char))
    end subroutine tensor_mpi_bcast_char_vec_l
    subroutine tensor_mpi_bcast_char_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      character*(*), intent(inout)              :: buffer
      type(tensor_mpi_stats_type),intent(inout) :: stats
      include "mpi_collective_vars.inc"
      include "mpi_bcast_std.inc"
    end subroutine tensor_mpi_bcast_char_basic
    subroutine tensor_mpi_bcast_char_vec_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      character, intent(inout)                  :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout) :: stats
      include "mpi_collective_vars.inc"
      include "mpi_bcast_std.inc"
    end subroutine tensor_mpi_bcast_char_vec_basic


    !DOUBLE PRECISION REDUCE
    subroutine tensor_mpi_reduce_tensor_dp(b,root,comm)
      implicit none
      real(tensor_dp), intent(inout) :: b
      real(tensor_dp) :: buffer(1)
      include "mpi_collective_vars.inc"
      n=1
      buffer(1) = b
      call tensor_mpi_reduce_dp_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_reduce,tensor_mpi_idx_tensor_dp))
      b = buffer(1)
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
      real(tensor_dp), intent(inout)               :: buffer(:)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      !CHANGE THIS ACCORDING TO THE DATATYPE OF BUFFER
      real(tensor_dp) :: noelm = 0.0E0_tensor_dp
      include "mpi_collective_vars.inc"
      include "mpi_reduce_std.inc"
    end subroutine tensor_mpi_reduce_dp_basic

    !STANDARD INTEGER ALLREDUCE
    subroutine tensor_mpi_allreduce_std_int(b,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(inout) :: b
      include "mpi_collective_vars_noroot.inc"
      integer(kind=tensor_standard_int) :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_allreduce_std_int_basic(buffer,n,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_allreduce,tensor_mpi_idx_tensor_standard_int))
      b = buffer(1)
    end subroutine tensor_mpi_allreduce_std_int
    subroutine tensor_mpi_allreduce_std_int_s(buffer,n1,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in)    :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(n1)
      include "mpi_collective_vars_noroot.inc"
      n=n1
      call tensor_mpi_allreduce_std_int_basic(buffer,n,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_allreduce,tensor_mpi_idx_tensor_standard_int))
    end subroutine tensor_mpi_allreduce_std_int_s
    subroutine tensor_mpi_allreduce_std_int_l(buffer,n1,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in)        :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(n1)
      include "mpi_collective_vars_noroot.inc"
      n=n1
      call tensor_mpi_allreduce_std_int_basic(buffer,n,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_allreduce,tensor_mpi_idx_tensor_standard_int))
    end subroutine tensor_mpi_allreduce_std_int_l
    subroutine tensor_mpi_allreduce_std_int_basic(buffer,n1,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(:)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_collective_vars_noroot.inc"
      include "mpi_allreduce_std.inc"
    end subroutine tensor_mpi_allreduce_std_int_basic

    !LONG INTEGER ALLREDUCE
    subroutine tensor_mpi_allreduce_long_int(b,comm)
      implicit none
      integer(kind=tensor_long_int), intent(inout) :: b
      include "mpi_collective_vars_noroot.inc"
      integer(kind=tensor_long_int) :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_allreduce_long_int_basic(buffer,n,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_allreduce,tensor_mpi_idx_tensor_long_int))
      b = buffer(1)
    end subroutine tensor_mpi_allreduce_long_int
    subroutine tensor_mpi_allreduce_long_int_s(buffer,n1,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      integer(kind=tensor_long_int), intent(inout)  :: buffer(n1)
      include "mpi_collective_vars_noroot.inc"
      n=n1
      call tensor_mpi_allreduce_long_int_basic(buffer,n,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_allreduce,tensor_mpi_idx_tensor_long_int))
    end subroutine tensor_mpi_allreduce_long_int_s
    subroutine tensor_mpi_allreduce_long_int_l(buffer,n1,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      integer(kind=tensor_long_int), intent(inout) :: buffer(n1)
      include "mpi_collective_vars_noroot.inc"
      n=n1
      call tensor_mpi_allreduce_long_int_basic(buffer,n,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_allreduce,tensor_mpi_idx_tensor_long_int))
    end subroutine tensor_mpi_allreduce_long_int_l
    subroutine tensor_mpi_allreduce_long_int_basic(buffer,n1,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      integer(kind=tensor_long_int), intent(inout) :: buffer(:)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_collective_vars_noroot.inc"
      include "mpi_allreduce_std.inc"
    end subroutine tensor_mpi_allreduce_long_int_basic

    !DOUBLE PRECISION ALLREDUCE
    subroutine tensor_mpi_allreduce_tensor_dp(b,comm)
      implicit none
      real(tensor_dp), intent(inout) :: b
      include "mpi_collective_vars_noroot.inc"
      real(tensor_dp) :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_allreduce_dp_basic(buffer,n,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_allreduce,tensor_mpi_idx_tensor_dp))
      b = buffer(1)
    end subroutine tensor_mpi_allreduce_tensor_dp
    subroutine tensor_mpi_allreduce_tensor_dp_s(buffer,n1,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      real(tensor_dp), intent(inout)                :: buffer(n1)
      include "mpi_collective_vars_noroot.inc"
      n=n1
      call tensor_mpi_allreduce_dp_basic(buffer,n,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_allreduce,tensor_mpi_idx_tensor_dp))
    end subroutine tensor_mpi_allreduce_tensor_dp_s
    subroutine tensor_mpi_allreduce_tensor_dp_l(buffer,n1,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      real(tensor_dp), intent(inout)            :: buffer(n1)
      include "mpi_collective_vars_noroot.inc"
      n=n1
      call tensor_mpi_allreduce_dp_basic(buffer,n,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_allreduce,tensor_mpi_idx_tensor_dp))
    end subroutine tensor_mpi_allreduce_tensor_dp_l
    subroutine tensor_mpi_allreduce_dp_basic(buffer,n1,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      real(tensor_dp), intent(inout)               :: buffer(:)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_collective_vars_noroot.inc"
      include "mpi_allreduce_std.inc"
    end subroutine tensor_mpi_allreduce_dp_basic

    !DOUBLE PRECISION SENDRCV
    subroutine tensor_mpi_sendrecv_tensor_dp(b,comm,sender,receiver)
      implicit none
      real(tensor_dp) :: b
      real(tensor_dp) :: buffer(1)
      include "mpi_sendrecv_vars.inc"
      n=1
      buffer(1) = b
      call tensor_mpi_sendrecv_tensor_dp_basic(buffer,n,comm,sender,receiver,&
         & tensor_mpi_stats(tensor_mpi_idx_sendrecv,tensor_mpi_idx_tensor_dp))
      if( rank == receiver) b = buffer(1)
    end subroutine tensor_mpi_sendrecv_tensor_dp
    subroutine tensor_mpi_sendrecv_tensor_dp_s(buffer,n1,comm,sender,receiver)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      real(tensor_dp), intent(inout)            :: buffer(n1)
      include "mpi_sendrecv_vars.inc"
      n=n1
      call tensor_mpi_sendrecv_tensor_dp_basic(buffer,n,comm,sender,receiver,&
         & tensor_mpi_stats(tensor_mpi_idx_sendrecv,tensor_mpi_idx_tensor_dp))
    end subroutine tensor_mpi_sendrecv_tensor_dp_s
    subroutine tensor_mpi_sendrecv_tensor_dp_l(buffer,n1,comm,sender,receiver)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      real(tensor_dp), intent(inout)            :: buffer(n1)
      include "mpi_sendrecv_vars.inc"
      n=n1
      call tensor_mpi_sendrecv_tensor_dp_basic(buffer,n,comm,sender,receiver,&
         & tensor_mpi_stats(tensor_mpi_idx_sendrecv,tensor_mpi_idx_tensor_dp))
    end subroutine tensor_mpi_sendrecv_tensor_dp_l
    subroutine tensor_mpi_sendrecv_tensor_dp_basic(buffer,n1,comm,sender,receiver,stats)
      implicit none
      integer(kind=tensor_long_int)  :: n1
      real(tensor_dp), intent(inout) :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_sendrecv_vars.inc"
      include "mpi_sendrecv_std.inc"
    end subroutine tensor_mpi_sendrecv_tensor_dp_basic

    !STANDARD INTEGER SENDRCV
    subroutine tensor_mpi_sendrecv_tensor_standard_int(b,comm,sender,receiver)
      implicit none
      integer(kind=tensor_standard_int) :: b
      integer(kind=tensor_standard_int) :: buffer(1)
      include "mpi_sendrecv_vars.inc"
      n=1
      buffer(1) = b
      call tensor_mpi_sendrecv_tensor_standard_int_basic(buffer,n,comm,sender,receiver,&
         & tensor_mpi_stats(tensor_mpi_idx_sendrecv,tensor_mpi_idx_tensor_standard_int))
      if( rank == receiver ) b = buffer(1)
    end subroutine tensor_mpi_sendrecv_tensor_standard_int
    subroutine tensor_mpi_sendrecv_tensor_standard_int_s(buffer,n1,comm,sender,receiver)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(n1)
      include "mpi_sendrecv_vars.inc"
      n=n1
      call tensor_mpi_sendrecv_tensor_standard_int_basic(buffer,n,comm,sender,receiver,&
         & tensor_mpi_stats(tensor_mpi_idx_sendrecv,tensor_mpi_idx_tensor_standard_int))
    end subroutine tensor_mpi_sendrecv_tensor_standard_int_s
    subroutine tensor_mpi_sendrecv_tensor_standard_int_l(buffer,n1,comm,sender,receiver)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(n1)
      include "mpi_sendrecv_vars.inc"
      n=n1
      call tensor_mpi_sendrecv_tensor_standard_int_basic(buffer,n,comm,sender,receiver,&
         & tensor_mpi_stats(tensor_mpi_idx_sendrecv,tensor_mpi_idx_tensor_standard_int))
    end subroutine tensor_mpi_sendrecv_tensor_standard_int_l
    subroutine tensor_mpi_sendrecv_tensor_standard_int_basic(buffer,n1,comm,sender,receiver,stats)
      implicit none
      integer(kind=tensor_long_int)  :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_sendrecv_vars.inc"
      include "mpi_sendrecv_std.inc"
    end subroutine tensor_mpi_sendrecv_tensor_standard_int_basic

    !LONG INTEGER SENDRCV
    subroutine tensor_mpi_sendrecv_tensor_long_int(b,comm,sender,receiver)
      implicit none
      integer(kind=tensor_long_int) :: b
      integer(kind=tensor_long_int) :: buffer(1)
      include "mpi_sendrecv_vars.inc"
      n=1
      buffer(1) = b
      call tensor_mpi_sendrecv_tensor_long_int_basic(buffer,n,comm,sender,receiver,&
         & tensor_mpi_stats(tensor_mpi_idx_sendrecv,tensor_mpi_idx_tensor_long_int))
      if( rank == receiver ) b = buffer(1)
    end subroutine tensor_mpi_sendrecv_tensor_long_int
    subroutine tensor_mpi_sendrecv_tensor_long_int_s(buffer,n1,comm,sender,receiver)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      integer(kind=tensor_long_int), intent(inout)            :: buffer(n1)
      include "mpi_sendrecv_vars.inc"
      n=n1
      call tensor_mpi_sendrecv_tensor_long_int_basic(buffer,n,comm,sender,receiver,&
         & tensor_mpi_stats(tensor_mpi_idx_sendrecv,tensor_mpi_idx_tensor_long_int))
    end subroutine tensor_mpi_sendrecv_tensor_long_int_s
    subroutine tensor_mpi_sendrecv_tensor_long_int_l(buffer,n1,comm,sender,receiver)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      integer(kind=tensor_long_int), intent(inout)            :: buffer(n1)
      include "mpi_sendrecv_vars.inc"
      n=n1
      call tensor_mpi_sendrecv_tensor_long_int_basic(buffer,n,comm,sender,receiver,&
         & tensor_mpi_stats(tensor_mpi_idx_sendrecv,tensor_mpi_idx_tensor_long_int))
    end subroutine tensor_mpi_sendrecv_tensor_long_int_l
    subroutine tensor_mpi_sendrecv_tensor_long_int_basic(buffer,n1,comm,sender,receiver,stats)
      implicit none
      integer(kind=tensor_long_int)  :: n1
      integer(kind=tensor_long_int), intent(inout) :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_sendrecv_vars.inc"
      include "mpi_sendrecv_std.inc"
    end subroutine tensor_mpi_sendrecv_tensor_long_int_basic

    !DOUBLE PRECISION PUT/RPUT
    subroutine tensor_mpi_put_dp_s(buffer,n1,pos,dest,win,fl)
       implicit none
       integer(kind=tensor_standard_int)         :: n1
       real(tensor_dp),intent(in)                :: buffer(n1)
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_put_dp_basic(buffer,n,pos,dest,win, &
          & tensor_mpi_stats(tensor_mpi_idx_put,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_put_dp_s
    subroutine tensor_mpi_put_dp_l(buffer,n1,pos,dest,win,fl)
       implicit none
       integer(kind=tensor_long_int)             :: n1
       real(tensor_dp),intent(in)                :: buffer(n1)
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_put_dp_basic(buffer,n,pos,dest,win,&
          & tensor_mpi_stats(tensor_mpi_idx_put,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_put_dp_l
    subroutine tensor_mpi_put_dp_basic(buffer,n1,pos,dest,win,stats,fl)
       implicit none
       integer(kind=tensor_long_int)             :: n1
       real(tensor_dp),intent(in)                :: buffer(n1)
       type(tensor_mpi_stats_type),intent(inout) :: stats
       include "mpi_one_sided_vars.inc"
       include "mpi_put_std.inc"
    end subroutine tensor_mpi_put_dp_basic

    subroutine tensor_mpi_rput_dp_s(buffer,n1,pos,dest,win,req,fl)
       implicit none
       integer(kind=tensor_standard_int)            :: n1
       real(tensor_dp),intent(in)                   :: buffer(n1)
       integer(kind=tensor_mpi_kind), intent(inout) :: req
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_rput_dp_basic(buffer,n,pos,dest,win,req,&
          & tensor_mpi_stats(tensor_mpi_idx_put,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_rput_dp_s
    subroutine tensor_mpi_rput_dp_l(buffer,n1,pos,dest,win,req,fl)
       implicit none
       integer(kind=tensor_long_int)                :: n1
       real(tensor_dp),intent(in)                   :: buffer(n1)
       integer(kind=tensor_mpi_kind), intent(inout) :: req
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_rput_dp_basic(buffer,n,pos,dest,win,req,&
          & tensor_mpi_stats(tensor_mpi_idx_put,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_rput_dp_l
    subroutine tensor_mpi_rput_dp_basic(buffer,n1,pos,dest,win,req,stats,fl)
       implicit none
       integer(kind=tensor_long_int)                :: n1
       real(tensor_dp),intent(in)                   :: buffer(n1)
       integer(kind=tensor_mpi_kind), intent(inout) :: req
       type(tensor_mpi_stats_type),intent(inout)    :: stats
       include "mpi_one_sided_vars.inc"
#ifdef VAR_HAVE_MPI3
       include "mpi_rput_std.inc"
#else
       call tensor_status_quit("ERROR(tensor_mpi_rput_dp_basic): MPI_RPUT can only be used with MPI 3.0+ ",-1)
#endif
    end subroutine tensor_mpi_rput_dp_basic

    !DOUBLE PRECISION GET/RGET
    subroutine tensor_mpi_get_dp_s(buffer,n1,pos,dest,win,fl)
       implicit none
       integer(kind=tensor_standard_int)         :: n1
       real(tensor_dp),intent(in)                :: buffer(n1)
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_get_dp_basic(buffer,n,pos,dest,win, &
          & tensor_mpi_stats(tensor_mpi_idx_get,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_get_dp_s
    subroutine tensor_mpi_get_dp_l(buffer,n1,pos,dest,win,fl)
       implicit none
       integer(kind=tensor_long_int)             :: n1
       real(tensor_dp),intent(in)                :: buffer(n1)
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_get_dp_basic(buffer,n,pos,dest,win,&
          & tensor_mpi_stats(tensor_mpi_idx_get,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_get_dp_l
    subroutine tensor_mpi_get_dp_basic(buffer,n1,pos,dest,win,stats,fl)
       implicit none
       integer(kind=tensor_long_int)             :: n1
       real(tensor_dp),intent(in)                :: buffer(n1)
       type(tensor_mpi_stats_type),intent(inout) :: stats
       include "mpi_one_sided_vars.inc"
       include "mpi_get_std.inc"
    end subroutine tensor_mpi_get_dp_basic

    subroutine tensor_mpi_rget_dp_s(buffer,n1,pos,dest,win,req,fl)
       implicit none
       integer(kind=tensor_standard_int)            :: n1
       real(tensor_dp),intent(in)                   :: buffer(n1)
       integer(kind=tensor_mpi_kind), intent(inout) :: req
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_rget_dp_basic(buffer,n,pos,dest,win,req,&
          & tensor_mpi_stats(tensor_mpi_idx_get,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_rget_dp_s
    subroutine tensor_mpi_rget_dp_l(buffer,n1,pos,dest,win,req,fl)
       implicit none
       integer(kind=tensor_long_int)                :: n1
       real(tensor_dp),intent(in)                   :: buffer(n1)
       integer(kind=tensor_mpi_kind), intent(inout) :: req
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_rget_dp_basic(buffer,n,pos,dest,win,req,&
          & tensor_mpi_stats(tensor_mpi_idx_get,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_rget_dp_l
    subroutine tensor_mpi_rget_dp_basic(buffer,n1,pos,dest,win,req,stats,fl)
       implicit none
       integer(kind=tensor_long_int)                :: n1
       real(tensor_dp),intent(in)                   :: buffer(n1)
       integer(kind=tensor_mpi_kind), intent(inout) :: req
       type(tensor_mpi_stats_type),intent(inout)    :: stats
       include "mpi_one_sided_vars.inc"
#ifdef VAR_HAVE_MPI3
       include "mpi_rget_std.inc"
#else
       call tensor_status_quit("ERROR(tensor_mpi_rget_dp_basic): MPI_RPUT can only be used with MPI 3.0+ ",-1)
#endif
    end subroutine tensor_mpi_rget_dp_basic

    !DOUBLE PRECISION PUT/RPUT
    subroutine tensor_mpi_acc_dp_s(buffer,n1,pos,dest,win,fl)
       implicit none
       integer(kind=tensor_standard_int)         :: n1
       real(tensor_dp),intent(in)                :: buffer(n1)
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_acc_dp_basic(buffer,n,pos,dest,win, &
          & tensor_mpi_stats(tensor_mpi_idx_acc,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_acc_dp_s
    subroutine tensor_mpi_acc_dp_l(buffer,n1,pos,dest,win,fl)
       implicit none
       integer(kind=tensor_long_int)             :: n1
       real(tensor_dp),intent(in)                :: buffer(n1)
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_acc_dp_basic(buffer,n,pos,dest,win,&
          & tensor_mpi_stats(tensor_mpi_idx_acc,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_acc_dp_l
    subroutine tensor_mpi_acc_dp_basic(buffer,n1,pos,dest,win,stats,fl)
       implicit none
       integer(kind=tensor_long_int)             :: n1
       real(tensor_dp),intent(in)                :: buffer(n1)
       type(tensor_mpi_stats_type),intent(inout) :: stats
       include "mpi_one_sided_vars.inc"
       include "mpi_acc_std.inc"
    end subroutine tensor_mpi_acc_dp_basic

    subroutine tensor_mpi_racc_dp_s(buffer,n1,pos,dest,win,req,fl)
       implicit none
       integer(kind=tensor_standard_int)            :: n1
       real(tensor_dp),intent(in)                   :: buffer(n1)
       integer(kind=tensor_mpi_kind), intent(inout) :: req
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_racc_dp_basic(buffer,n,pos,dest,win,req,&
          & tensor_mpi_stats(tensor_mpi_idx_acc,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_racc_dp_s
    subroutine tensor_mpi_racc_dp_l(buffer,n1,pos,dest,win,req,fl)
       implicit none
       integer(kind=tensor_long_int)                :: n1
       real(tensor_dp),intent(in)                   :: buffer(n1)
       integer(kind=tensor_mpi_kind), intent(inout) :: req
       include "mpi_one_sided_vars.inc"
       n=n1
       call tensor_mpi_racc_dp_basic(buffer,n,pos,dest,win,req,&
          & tensor_mpi_stats(tensor_mpi_idx_acc,tensor_mpi_idx_tensor_dp),fl=fl)
    end subroutine tensor_mpi_racc_dp_l
    subroutine tensor_mpi_racc_dp_basic(buffer,n1,pos,dest,win,req,stats,fl)
       implicit none
       integer(kind=tensor_long_int)                :: n1
       real(tensor_dp),intent(in)                   :: buffer(n1)
       integer(kind=tensor_mpi_kind), intent(inout) :: req
       type(tensor_mpi_stats_type),intent(inout)    :: stats
       include "mpi_one_sided_vars.inc"
#ifdef VAR_HAVE_MPI3
       include "mpi_racc_std.inc"
#else
       call tensor_status_quit("ERROR(tensor_mpi_racc_dp_basic): MPI_RPUT can only be used with MPI 3.0+ ",-1)
#endif
    end subroutine tensor_mpi_racc_dp_basic


#endif
end module tensor_mpi_interface_module
