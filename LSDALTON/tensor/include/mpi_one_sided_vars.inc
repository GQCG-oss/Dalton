       integer(kind=tensor_long_int), intent(in) :: pos
       integer(kind=tensor_mpi_kind), intent(in) :: dest
       integer(kind=tensor_mpi_kind), intent(in) :: win
       logical, intent(in), optional             :: fl
       integer(kind=tensor_mpi_kind)  :: ierr, nMPI, datatype
       integer(kind=tensor_long_int)  :: n, chunk, first_el
       integer(kind=MPI_ADDRESS_KIND) :: offset
       logical :: f
