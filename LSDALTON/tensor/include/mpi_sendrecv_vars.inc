      integer(kind=tensor_mpi_kind),intent(in) :: comm
      integer(kind=tensor_mpi_kind),intent(in) :: sender, receiver
      integer(kind=tensor_mpi_kind) :: ierr, datatype, nMPI, tag, stat(MPI_STATUS_SIZE), rank
      integer(kind=tensor_long_int) :: n, chunk, first_el, datatype_size
