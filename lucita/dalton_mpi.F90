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

#ifndef VAR_USE_MPIF
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif

#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

  public dalton_mpi_bcast
  public dalton_mpi_reduce
  public dalton_mpi_allreduce

  private

  save 

  integer, private                       :: istat(MPI_STATUS_SIZE)
  integer, private                       :: ierr

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

    call mpi_bcast(x, 1, my_MPI_INTEGER, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_i1(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), my_MPI_INTEGER, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_i2(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), my_MPI_INTEGER, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_i3(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), my_MPI_INTEGER, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_i4(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), my_MPI_INTEGER, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_r0(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, 1, mpi_double_precision, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_r1(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_double_precision, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_r2(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_double_precision, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_r3(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_double_precision, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_r4(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_double_precision, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_l0(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    logical :: x
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, 1, mpi_logical, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_l1(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    logical :: x(:)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_logical, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_l2(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    logical :: x(:, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, size(x), mpi_logical, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_bcast_c0(x,root_proc,communicator)

!   ----------------------------------------------------------------------------
    character (len=72) :: x
    integer            :: root_proc
    integer            :: communicator
!   ----------------------------------------------------------------------------

    call mpi_bcast(x, 72, mpi_character, root_proc, communicator, ierr)

  end subroutine

!   ************************* reduce interface *********************************
  subroutine dalton_mpi_reduce_i0(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x
    integer :: y
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_reduce(x, y, 1, my_MPI_INTEGER, mpi_sum, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_reduce_i1(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:)
    integer :: y(:)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_reduce(x, y, size(x), my_MPI_INTEGER, mpi_sum, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_reduce_i2(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :)
    integer :: y(:, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_reduce(x, y, size(x), my_MPI_INTEGER, mpi_sum, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_reduce_i3(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :)
    integer :: y(:, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_reduce(x, y, size(x), my_MPI_INTEGER, mpi_sum, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_reduce_i4(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :, :)
    integer :: y(:, :, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_reduce(x, y, size(x), my_MPI_INTEGER, mpi_sum, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_reduce_r0(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x
    real(8) :: y
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_reduce(x, y, 1, mpi_double_precision, mpi_sum, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_reduce_r1(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:)
    real(8) :: y(:)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_reduce(x, y, size(x), mpi_double_precision, mpi_sum, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_reduce_r2(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :)
    real(8) :: y(:, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_reduce(x, y, size(x), mpi_double_precision, mpi_sum, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_reduce_r3(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :)
    real(8) :: y(:, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_reduce(x, y, size(x), mpi_double_precision, mpi_sum, root_proc, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_reduce_r4(x,y,root_proc,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :, :)
    real(8) :: y(:, :, :, :)
    integer :: root_proc
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_reduce(x, y, size(x), mpi_double_precision, mpi_sum, root_proc, communicator, ierr)

  end subroutine

!   ************************* allreduce interface ******************************
  subroutine dalton_mpi_allreduce_i0(x,y,communicator)

!   ----------------------------------------------------------------------------
    integer :: x
    integer :: y
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_allreduce(x, y, 1, my_MPI_INTEGER, mpi_sum, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_allreduce_i1(x,y,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:)
    integer :: y(:)
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_allreduce(x, y, size(y), my_MPI_INTEGER, mpi_sum, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_allreduce_i2(x,y,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :)
    integer :: y(:, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_allreduce(x, y, size(y), my_MPI_INTEGER, mpi_sum, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_allreduce_i3(x,y,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :)
    integer :: y(:, :, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_allreduce(x, y, size(y), my_MPI_INTEGER, mpi_sum, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_allreduce_i4(x,y,communicator)

!   ----------------------------------------------------------------------------
    integer :: x(:, :, :, :)
    integer :: y(:, :, :, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_allreduce(x, y, size(y), my_MPI_INTEGER, mpi_sum, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_allreduce_r0(x,y,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x
    real(8) :: y
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_allreduce(x, y, 1, mpi_double_precision, mpi_sum, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_allreduce_r1(x,y,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:)
    real(8) :: y(:)
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_allreduce(x, y, size(y), mpi_double_precision, mpi_sum, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_allreduce_r2(x,y,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :)
    real(8) :: y(:, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_allreduce(x, y, size(y), mpi_double_precision, mpi_sum, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_allreduce_r3(x,y,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :)
    real(8) :: y(:, :, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_allreduce(x, y, size(y), mpi_double_precision, mpi_sum, communicator, ierr)

  end subroutine

  subroutine dalton_mpi_allreduce_r4(x,y,communicator)

!   ----------------------------------------------------------------------------
    real(8) :: x(:, :, :, :)
    real(8) :: y(:, :, :, :)
    integer :: communicator
!   ----------------------------------------------------------------------------

    call mpi_allreduce(x, y, size(y), mpi_double_precision, mpi_sum, communicator, ierr)

  end subroutine
#else /* ifdef VAR_MPI */
module dummy_dalton_mpi
#endif
end module
