program test
   use mpi_f08
   TYPE(MPI_Comm) :: comm1
   integer :: ierr
   call MPI_INIT(ierr)
   call MPI_COMM_DUP(MPI_COMM_WORLD,comm1,ierr)
   call MPI_FINALIZE(ierr)
end program test
