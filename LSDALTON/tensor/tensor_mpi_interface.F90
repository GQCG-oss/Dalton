module tensor_mpi_interface_module
   use tensor_parameters_and_counters
   use tensor_mpi_binding_module
   use tensor_type_def_module

   !NOTE THAT ONLY MPI_IN_PLACE HAS BEEN IMPLEMENTED FOR THE MPI COLLECTIVES
   public :: tensor_mpi_bcast
   public :: tensor_mpi_reduce
   public :: tensor_mpi_allreduce
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
                     & tensor_mpi_bcast_log
#else
      module procedure tensor_mpi_dummy
#endif
   end interface tensor_mpi_bcast

   interface tensor_mpi_reduce
#ifdef VAR_MPI
      module procedure tensor_mpi_dummy
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
    !STD INTEGER BCAST
    subroutine tensor_mpi_bcast_std_int(b,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(inout) :: b
      include "mpi_bcast_vars.inc"
      integer(kind=tensor_standard_int) :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_std_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_standard_int))
      b = buffer(1)
    end subroutine tensor_mpi_bcast_std_int
    subroutine tensor_mpi_bcast_std_int_s(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in)    :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(n1)
      include "mpi_bcast_vars.inc"
      n=n1
      call tensor_mpi_bcast_std_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_standard_int))
    end subroutine tensor_mpi_bcast_std_int_s
    subroutine tensor_mpi_bcast_std_int_l(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(n1)
      include "mpi_bcast_vars.inc"
      n=n1
      call tensor_mpi_bcast_std_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_standard_int))
    end subroutine tensor_mpi_bcast_std_int_l

    subroutine tensor_mpi_bcast_std_int_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)        :: n1
      integer(kind=tensor_standard_int), intent(inout) :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)        :: stats
      include "mpi_bcast_vars.inc"
      include "mpi_bcast_std.inc"
    end subroutine tensor_mpi_bcast_std_int_basic

    !LONG INTEGER BCAST
    subroutine tensor_mpi_bcast_long_int(b,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(inout)  :: b
      include "mpi_bcast_vars.inc"
      integer(kind=tensor_long_int) :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_long_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_long_int))
      b = buffer(1)
    end subroutine tensor_mpi_bcast_long_int
    subroutine tensor_mpi_bcast_long_int_s(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      integer(kind=tensor_long_int), intent(inout)  :: buffer(n1)
      include "mpi_bcast_vars.inc"
      n=n1
      call tensor_mpi_bcast_long_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_long_int))
    end subroutine tensor_mpi_bcast_long_int_s
    subroutine tensor_mpi_bcast_long_int_l(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      integer(kind=tensor_long_int), intent(inout) :: buffer(n1)
      include "mpi_bcast_vars.inc"
      n=n1
      call tensor_mpi_bcast_long_int_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_long_int))
    end subroutine tensor_mpi_bcast_long_int_l

    subroutine tensor_mpi_bcast_long_int_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      integer(kind=tensor_long_int), intent(inout) :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_bcast_vars.inc"
      include "mpi_bcast_std.inc"
    end subroutine tensor_mpi_bcast_long_int_basic

    !LOGICAL BCAST
    subroutine tensor_mpi_bcast_log(b,root,comm)
      implicit none
      logical, intent(inout) :: b
      include "mpi_bcast_vars.inc"
      logical :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_log_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_log))
      b = buffer(1)
    end subroutine tensor_mpi_bcast_log
    subroutine tensor_mpi_bcast_log_s(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in)    :: n1
      logical, intent(inout) :: buffer(n1)
      include "mpi_bcast_vars.inc"
      n=n1
      call tensor_mpi_bcast_log_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_log))
    end subroutine tensor_mpi_bcast_log_s
    subroutine tensor_mpi_bcast_log_l(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      logical, intent(inout) :: buffer(n1)
      include "mpi_bcast_vars.inc"
      n=n1
      call tensor_mpi_bcast_log_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_log))
    end subroutine tensor_mpi_bcast_log_l

    subroutine tensor_mpi_bcast_log_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      logical, intent(inout)                       :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_bcast_vars.inc"
      include "mpi_bcast_std.inc"
    end subroutine tensor_mpi_bcast_log_basic

    !DOUBLE PRECISION BCAST
    subroutine tensor_mpi_bcast_tensor_dp(b,root,comm)
      implicit none
      real(tensor_dp), intent(inout) :: b
      include "mpi_bcast_vars.inc"
      real(tensor_dp) :: buffer(1)
      n=1
      buffer(1) = b
      call tensor_mpi_bcast_dp_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_dp))
      b = buffer(1)
    end subroutine tensor_mpi_bcast_tensor_dp
    subroutine tensor_mpi_bcast_tensor_dp_s(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_standard_int), intent(in) :: n1
      real(tensor_dp), intent(inout)                :: buffer(n1)
      include "mpi_bcast_vars.inc"
      n=n1
      call tensor_mpi_bcast_dp_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_dp))
    end subroutine tensor_mpi_bcast_tensor_dp_s
    subroutine tensor_mpi_bcast_tensor_dp_l(buffer,n1,root,comm)
      implicit none
      integer(kind=tensor_long_int), intent(in) :: n1
      real(tensor_dp), intent(inout)            :: buffer(n1)
      include "mpi_bcast_vars.inc"
      n=n1
      call tensor_mpi_bcast_dp_basic(buffer,n,root,comm,&
         & tensor_mpi_stats(tensor_mpi_idx_bcast,tensor_mpi_idx_tensor_dp))
    end subroutine tensor_mpi_bcast_tensor_dp_l
    subroutine tensor_mpi_bcast_dp_basic(buffer,n1,root,comm,stats)
      implicit none
      integer(kind=tensor_long_int), intent(in)    :: n1
      real(tensor_dp), intent(inout)               :: buffer(n1)
      type(tensor_mpi_stats_type),intent(inout)    :: stats
      include "mpi_bcast_vars.inc"
      include "mpi_bcast_std.inc"
    end subroutine tensor_mpi_bcast_dp_basic
#endif

end module tensor_mpi_interface_module
