module memory_parallel

  implicit none

#ifdef VAR_MPI
#include "mpif.h"

  public memallocmpi
  public memfreempi

  private

  integer(kind=MPI_INTEGER_KIND) :: ierr

contains

  subroutine memallocmpi(nelement, ptr)
  integer(kind=mpi_address_kind), intent(inout) :: ptr
  integer(kind=mpi_address_kind), intent(in)    :: nelement
  call mpi_alloc_mem(nelement, mpi_info_null, ptr, ierr)
  end subroutine memallocmpi

  subroutine memfreempi(buf)
  real(8), intent(inout) :: buf(*)
  call mpi_free_mem(buf, ierr)
  end subroutine memfreempi
#endif

end module
