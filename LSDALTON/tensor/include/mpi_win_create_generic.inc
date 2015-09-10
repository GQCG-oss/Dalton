
       n=n1

       info = MPI_INFO_NULL
       call MPI_TYPE_GET_EXTENT(datatype,lb,mpi_size,ierr)

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_win_create_generic): failed",220)

       b = n*mpi_size
       s = int(mpi_size,kind=tensor_mpi_kind)

       call mpi_win_create(darr,b,s,info,comm,win,ierr)

       if(ierr /= 0) call tensor_status_quit("ERROR(tensor_mpi_win_create_generic): failed",220)
