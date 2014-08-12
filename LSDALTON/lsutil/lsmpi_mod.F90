module lsmpi_module

#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
  use mpi
#else
  include 'mpif.h'
#endif
#endif

  contains

  subroutine void_mpi_sub
  end subroutine void_mpi_sub

end module lsmpi_module
