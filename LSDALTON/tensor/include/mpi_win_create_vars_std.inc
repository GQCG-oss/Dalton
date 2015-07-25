
       integer(kind=tensor_mpi_kind),intent(inout) :: win,comm
       integer(kind=tensor_mpi_kind) :: datatype
       integer(kind=tensor_mpi_kind)  :: ierr=0
       integer(kind=tensor_mpi_kind)  :: info,s
       integer(kind=MPI_ADDRESS_KIND) :: mpi_size,lb,b
