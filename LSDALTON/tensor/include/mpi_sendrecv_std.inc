      ierr        = 0
      chunk       = TENSOR_MPI_MSG_LEN
      n           = n1
      call tensor_get_rank_for_comm(comm,rank)
      tag         = 33

      datatype    = stats%d_mpi

      stats%bytes = stats%bytes + stats%size_ * n
      stats%time_ = stats%time_ - mpi_wtime()

      if( n <= chunk )then

         nMPI = n
         stats%ncall = stats%ncall + 1
         if( rank == sender)then
            call mpi_send(buffer,nMPI,datatype,receiver,tag, comm, ierr)
         else if( rank == receiver )then
            call mpi_recv(buffer,nMPI,datatype,sender, tag, comm,stat,ierr)
         else
            call tensor_status_quit("ERROR(mpi_sendrecv_generic): wrong MPI rank in sendrecv",233)
         endif

      else

         do first_el=1,n,chunk

            if( n-first_el+1 < chunk )then
               nMPI = mod(n,chunk)
            else
               nMPI = chunk
            endif

            stats%ncall = stats%ncall + 1

            !call mpi_sendrecv(buffer(first_el:first_el+nMPI-1), nMPI, datatype, receiver, tag,&
            !                & buffer(first_el:first_el+nMPI-1), nMPI, datatype, sender  , tag,&
            !                & comm, stat, ierr)

            if( rank == sender)then
               call mpi_send(buffer(first_el:first_el+nMPI-1),nMPI,datatype,receiver,tag, comm, ierr)
            else if( rank == receiver )then
               call mpi_recv(buffer(first_el:first_el+nMPI-1),nMPI,datatype,sender, tag, comm,stat,ierr)
            else
               call tensor_status_quit("ERROR(mpi_sendrecv_generic): wrong MPI rank in sendrecv",233)
            endif



         enddo
      endif
      if (ierr/= 0) call tensor_status_quit("ERROR(tensor_general_reduce): mpi returns ierr=",ierr)

      stats%time_ = stats%time_ + mpi_wtime()
