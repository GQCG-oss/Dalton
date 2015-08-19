program test

#ifdef USE_MPI_MOD_F08
   use mpi_f08
#else if defined(USE_MPI_MOD_F90)
   use mpi
#else
   include 'mpif.h'
#endif

   integer :: ierr

   call MPI_INIT(ierr)
   call MPI_FINALIZE(ierr)

end program test
