
module tensor_mpi_binding_module

#ifdef TENSORS_IN_LSDALTON
   use lsmpi_module
#else
#ifdef VAR_MPI

#ifdef USE_MPI_MOD_F08
  use mpi_f08
#elif defined(USE_MPI_MOD_F90)
  use mpi
#else
  include 'mpif.h'
#endif
!ENDIF VAR_MPI
#endif
!ENDIF TENSORS_IN_LSDALTON
#endif

  contains

  subroutine tensor_void_mpi()
     implicit none
  end subroutine tensor_void_mpi

end module tensor_mpi_binding_module
