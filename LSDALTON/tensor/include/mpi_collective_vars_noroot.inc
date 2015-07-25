      integer(kind=tensor_mpi_kind),intent(in) :: comm
      integer(kind=tensor_mpi_kind) :: rank, ierr, datatype, nMPI
      integer(kind=tensor_long_int) :: n, chunk, first_el, datatype_size
