program test_mpi

   use mpi

   implicit none
   integer(kind=4) :: ierr
   integer(kind=mpi_address_kind) :: lb, extent
   integer(kind=4) :: stat(MPI_STATUS_SIZE),comm
   integer(kind=4) :: s = MPI_ANY_SOURCE
   integer(kind=4) :: t = MPI_ANY_TAG

   call mpi_init(ierr)
   call mpi_type_get_extent(mpi_double_precision,lb,extent,ierr)
   call mpi_probe(s, t, comm, stat, ierr)
   call mpi_finalize(ierr)

end program test_mpi
