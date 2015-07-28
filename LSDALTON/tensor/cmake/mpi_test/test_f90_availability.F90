program test
   use mpi
   integer :: comm1,ierr
   call MPI_INIT(ierr)
   call MPI_COMM_DUP(MPI_COMM_WORLD,comm1,ierr)
   call MPI_FINALIZE(ierr)
end program test
