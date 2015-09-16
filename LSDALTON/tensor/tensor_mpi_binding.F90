
module tensor_mpi_binding_module

#ifdef VAR_MPI

#ifdef USE_MPI_MOD_F08
  use mpi_f08
#elif defined(USE_MPI_MOD_F90)
  use mpi
#else
  include 'mpif.h'
#endif

#endif

  contains

  subroutine tensor_void_mpi()
     implicit none
  end subroutine tensor_void_mpi

end module tensor_mpi_binding_module
