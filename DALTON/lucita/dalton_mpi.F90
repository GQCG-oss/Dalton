!dalton_copyright_start
!
!
!dalton_copyright_end

#ifdef VAR_MPI
module dalton_mpi

! stefan: - this module provides generic interfaces to major MPI routines
!           in parallel calculations.
!
!           written by radovan bast, for DIRAC.
!           adapted and extended for DALTON by sknecht, december 2010.

#ifdef USE_MPI_MOD_F90
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif

  public dalton_mpi_bcast
  public dalton_mpi_reduce
  public dalton_mpi_allreduce

  private

  save 


#if defined (VAR_INT64)
  integer(kind=MPI_INTEGER_KIND), private   :: my_MPI_INTEGER = MPI_INTEGER8
  integer(kind=MPI_INTEGER_KIND), private   :: my_MPI_LOGICAL = MPI_INTEGER8
#else
  integer(kind=MPI_INTEGER_KIND), private   :: my_MPI_INTEGER = MPI_INTEGER4
  integer(kind=MPI_INTEGER_KIND), private   :: my_MPI_LOGICAL = MPI_INTEGER4
#endif
  integer(kind=MPI_INTEGER_KIND), private   :: ierr_mpi, size_x, size_y
  integer(kind=MPI_INTEGER_KIND), private   :: root_proc_mpi, comm_mpi
  integer(kind=MPI_INTEGER_KIND), parameter :: one_mpi = 1

  interface dalton_mpi_bcast
    module procedure dalton_mpi_bcast_i0
    module procedure dalton_mpi_bcast_i1
    module procedure dalton_mpi_bcast_i2
    module procedure dalton_mpi_bcast_i3
    module procedure dalton_mpi_bcast_i4
    module procedure dalton_mpi_bcast_r0
    module procedure dalton_mpi_bcast_r1
    module procedure dalton_mpi_bcast_r2
    module procedure dalton_mpi_bcast_r3
    module procedure dalton_mpi_bcast_r4
    module procedure dalton_mpi_bcast_l0
    module procedure dalton_mpi_bcast_l1
    module procedure dalton_mpi_bcast_l2
    module procedure dalton_mpi_bcast_c0
  end interface

  interface dalton_mpi_reduce
    module procedure dalton_mpi_reduce_i0
    module procedure dalton_mpi_reduce_i1
    module procedure dalton_mpi_reduce_i2
    module procedure dalton_mpi_reduce_i3
    module procedure dalton_mpi_reduce_i4
    module procedure dalton_mpi_reduce_r0
    module procedure dalton_mpi_reduce_r1
    module procedure dalton_mpi_reduce_r2
    module procedure dalton_mpi_reduce_r3
    module procedure dalton_mpi_reduce_r4
  end interface

  interface dalton_mpi_allreduce
    module procedure dalton_mpi_allreduce_i0
    module procedure dalton_mpi_allreduce_i1
    module procedure dalton_mpi_allreduce_i2
    module procedure dalton_mpi_allreduce_i3
    module procedure dalton_mpi_allreduce_i4
    module procedure dalton_mpi_allreduce_r0
    module procedure dalton_mpi_allreduce_r1
    module procedure dalton_mpi_allreduce_r2
    module procedure dalton_mpi_allreduce_r3
    module procedure dalton_mpi_allreduce_r4
  end interface

contains

  subroutine dalton_mpi_bcast_i0(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, one_mpi, my_MPI_INTEGER, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_i1(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, size_x, my_MPI_INTEGER, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_i2(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, size_x, my_MPI_INTEGER, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_i3(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, size_x, my_MPI_INTEGER, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_i4(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, size_x, my_MPI_INTEGER, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_r0(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, one_mpi, mpi_double_precision, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_r1(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, size_x, mpi_double_precision, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_r2(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, size_x, mpi_double_precision, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_r3(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, size_x, mpi_double_precision, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_r4(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, size_x, mpi_double_precision, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_l0(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    logical :: x
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, one_mpi, my_MPI_LOGICAL, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_l1(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    logical :: x(:)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, size_x, my_MPI_LOGICAL, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_l2(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    logical :: x(:, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, size_x, my_MPI_LOGICAL, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_bcast_c0(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer(kind=MPI_INTEGER_KIND), parameter :: len_x_mpi = 72
    character (len=len_x_mpi) :: x
    integer            :: root_proc
    integer            :: communicator
!   ----------------------------------------------------------------------------

    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_bcast(x, len_x_mpi, mpi_character, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

!   ************************* reduce interface *********************************
  subroutine dalton_mpi_reduce_i0(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x
    integer :: y
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_reduce(x, y, one_mpi, my_MPI_INTEGER, mpi_sum, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_reduce_i1(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:)
    integer :: y(:)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_reduce(x, y, size_x, my_MPI_INTEGER, mpi_sum, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_reduce_i2(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :)
    integer :: y(:, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_reduce(x, y, size_x, my_MPI_INTEGER, mpi_sum, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_reduce_i3(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :)
    integer :: y(:, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_reduce(x, y, size_x, my_MPI_INTEGER, mpi_sum, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_reduce_i4(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :, :)
    integer :: y(:, :, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_reduce(x, y, size_x, my_MPI_INTEGER, mpi_sum, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_reduce_r0(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x
    real(8) :: y
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_reduce(x, y, one_mpi, mpi_double_precision, mpi_sum, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_reduce_r1(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:)
    real(8) :: y(:)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_reduce(x, y, size_x, mpi_double_precision, mpi_sum, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_reduce_r2(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :)
    real(8) :: y(:, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_reduce(x, y, size_x, mpi_double_precision, mpi_sum, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_reduce_r3(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :)
    real(8) :: y(:, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_reduce(x, y, size_x, mpi_double_precision, mpi_sum, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_reduce_r4(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :, :)
    real(8) :: y(:, :, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_x = size(x)
    root_proc_mpi = root_proc
    comm_mpi      = communicator
    call mpi_reduce(x, y, size_x, mpi_double_precision, mpi_sum, root_proc_mpi, comm_mpi, ierr_mpi)

  end subroutine

!   ************************* allreduce interface ******************************
  subroutine dalton_mpi_allreduce_i0(x,y,communicator)

!   ----------------------------------------------------------------------------
    integer :: x
    integer :: y
    integer :: communicator
!   ----------------------------------------------------------------------------

    comm_mpi      = communicator
    call mpi_allreduce(x, y, one_mpi, my_MPI_INTEGER, mpi_sum, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_allreduce_i1(x,y,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:)
    integer :: y(:)
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_y = size(y)
    comm_mpi      = communicator
    call mpi_allreduce(x, y, size_y, my_MPI_INTEGER, mpi_sum, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_allreduce_i2(x,y,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :)
    integer :: y(:, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_y = size(y)
    comm_mpi      = communicator
    call mpi_allreduce(x, y, size_y, my_MPI_INTEGER, mpi_sum, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_allreduce_i3(x,y,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :)
    integer :: y(:, :, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_y = size(y)
    comm_mpi      = communicator
    call mpi_allreduce(x, y, size_y, my_MPI_INTEGER, mpi_sum, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_allreduce_i4(x,y,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :, :)
    integer :: y(:, :, :, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_y = size(y)
    comm_mpi      = communicator
    call mpi_allreduce(x, y, size_y, my_MPI_INTEGER, mpi_sum, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_allreduce_r0(x,y,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x
    real(8) :: y
    integer :: communicator
!   ----------------------------------------------------------------------------

    comm_mpi      = communicator
    call mpi_allreduce(x, y, one_mpi, mpi_double_precision, mpi_sum, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_allreduce_r1(x,y,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:)
    real(8) :: y(:)
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_y = size(y)
    comm_mpi      = communicator
    call mpi_allreduce(x, y, size_y, mpi_double_precision, mpi_sum, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_allreduce_r2(x,y,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :)
    real(8) :: y(:, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_y = size(y)
    comm_mpi      = communicator
    call mpi_allreduce(x, y, size_y, mpi_double_precision, mpi_sum, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_allreduce_r3(x,y,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :)
    real(8) :: y(:, :, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_y = size(y)
    comm_mpi      = communicator
    call mpi_allreduce(x, y, size_y, mpi_double_precision, mpi_sum, comm_mpi, ierr_mpi)

  end subroutine

  subroutine dalton_mpi_allreduce_r4(x,y,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :, :)
    real(8) :: y(:, :, :, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    size_y = size(y)
    comm_mpi      = communicator
    call mpi_allreduce(x, y, size_y, mpi_double_precision, mpi_sum, comm_mpi, ierr_mpi)

  end subroutine
end module
#else /* ifdef VAR_MPI */
subroutine dalton_mpi
   call quit('ERROR: dummy dalton_mpi called')
end
#endif
