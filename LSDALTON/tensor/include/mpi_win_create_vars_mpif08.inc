
       type(MPI_Win), intent(inout)  :: win
       type(MPI_Comm), intent(inout) :: comm
       type(MPI_Datatype) :: datatype
       integer(kind=tensor_mpi_kind)  :: ierr=0
       integer(kind=tensor_mpi_kind)  :: info,s
       integer(kind=MPI_ADDRESS_KIND) :: mpi_size,lb,b
       integer(kind=tensor_long_int)  :: n
