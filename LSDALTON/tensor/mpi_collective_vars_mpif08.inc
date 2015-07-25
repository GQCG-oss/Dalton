      type(MPI_Comm),intent(in)                :: comm
      integer(kind=tensor_mpi_kind),intent(in) :: root
      type(MPI_Datatype)            :: datatype
      integer(kind=tensor_mpi_kind) :: rank, ierr, nMPI
      integer(kind=tensor_long_int) :: n, chunk, first_el, datatype_size
