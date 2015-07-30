       ierr        = 0
       chunk       = TENSOR_MPI_MSG_LEN
       n           = n1
       f           = .false.
       if(present(fl)) f = fl

       datatype    = stats%d_mpi

       stats%bytes = stats%bytes + stats%size_ * n
       stats%time_ = stats%time_ - mpi_wtime()

       do first_el=1,n,chunk

          if( n-first_el+1 < chunk )then
             nMPI = mod(n,chunk)
          else
             nMPI = chunk
          endif

          offset = int(first_el-1+pos-1,kind=MPI_ADDRESS_KIND)

          call mpi_accumulate(buffer(first_el:first_el+nMPI-1),nMPI,datatype,dest,offset,nMPI,datatype,MPI_SUM,win,ierr)

          if (ierr/= 0) call tensor_status_quit("ERROR(tensor_general_put): mpi returns ierr=",ierr)

          if( f ) call tensor_mpi_win_flush(win,local=.true.)
       enddo

       stats%time_ = stats%time_ + mpi_wtime()
