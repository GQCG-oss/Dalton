module memory_parallel
implicit none
integer :: blubb, ierr

#ifdef VAR_MPI
public memallocmpi
public memfreempi

contains
subroutine memallocmpi(nelement, ptr)
#include "mpif.h"
integer(kind=mpi_address_kind), intent(inout) :: ptr
integer(kind=mpi_address_kind), intent(in)    :: nelement
call mpi_alloc_mem(nelement, mpi_info_null, ptr, ierr)
end subroutine memallocmpi

subroutine memfreempi(buf)
#include "mpif.h"
real(8), intent(inout) :: buf(*)
call mpi_free_mem(buf, ierr)
end subroutine memfreempi
#endif
end module
