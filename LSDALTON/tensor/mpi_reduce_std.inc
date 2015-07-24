      ierr        = 0
      chunk       = TENSOR_MPI_MSG_LEN
      n           = n1
      !call tensor_get_rank_for_comm(comm,rank)

      datatype    = stats%d_mpi

      stats%bytes = stats%bytes + stats%size_ * n
      stats%time_ = stats%time_ - mpi_wtime()

      do first_el=1,n,chunk

         if( n-first_el+1 < chunk )then
            nMPI = mod(n,chunk)
         else
            nMPI = chunk
         endif

         stats%ncall = stats%ncall + 1

         !if( rank == root)then
            call mpi_reduce(buffer(first_el:first_el+nMPI-1),MPI_IN_PLACE,nMPI,datatype,MPI_SUM,root,comm,ierr)
         !else
         !   call mpi_reduce(buffer(first_el:first_el+nMPI-1),null(),nMPI,datatype,MPI_SUM,root,comm,ierr)
         !endif

         if (ierr/= 0) call tensor_status_quit("ERROR(tensor_general_bacst): mpi returns ierr=",ierr)

      enddo

      stats%time_ = stats%time_ + mpi_wtime()
