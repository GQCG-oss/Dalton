      ierr        = 0
      chunk       = TENSOR_MPI_MSG_LEN
      n           = n1

      datatype    = stats%d_mpi

      stats%bytes = stats%bytes + stats%size_ * n
      stats%time_ = stats%time_ - mpi_wtime()

      if( n <= chunk )then

         nMPI = n
         stats%ncall = stats%ncall + 1
         call mpi_bcast(buffer,nMPI,datatype,root,comm,ierr)

      else

         do first_el=1,n,chunk

            if( n-first_el+1 < chunk )then
               nMPI = mod(n,chunk)
            else
               nMPI = chunk
            endif

            stats%ncall = stats%ncall + 1

            call mpi_bcast(buffer(first_el:first_el+nMPI-1),nMPI,datatype,root,comm,ierr)

         enddo
      endif

      if (ierr/= 0) call tensor_status_quit("ERROR(tensor_general_bacst): mpi returns ierr=",ierr)

      stats%time_ = stats%time_ + mpi_wtime()