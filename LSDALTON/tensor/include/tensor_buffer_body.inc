      integer, intent(in), optional :: init_size
      logical, optional :: finalize

!#ifdef USE_MPI_MOD_F08
!      type(MPI_Comm),intent(in),optional :: comm
!      integer(kind=tensor_mpi_kind),intent(in), optional :: root
!#else
      integer(kind=tensor_mpi_kind),intent(in), optional :: comm
      integer(kind=tensor_mpi_kind),intent(in), optional :: root
!#endif

      integer(kind=tensor_long_int) :: isize, tsize
      integer(kind=tensor_long_int) :: b,e,n
      integer(kind=tensor_mpi_kind) :: r, c, rank
      logical :: fin

      n = n1


      if(tensor_mpi_buffer%is_initialized)then
         r = tensor_mpi_buffer%root
         c = tensor_mpi_buffer%comm
      else if(present(root).and.present(comm))then
         r     = root
         c     = comm
         isize = tensor_mpi_buffer%std_size
      else
         call tensor_status_quit("ERROR(tensor_buffer_add): root and comm need to/can be specified in the FIRST call, only",33)
      endif

      call tensor_get_rank_for_comm(c,rank)

      fin = .false.
      if(present(finalize))  fin = finalize
      if(present(init_size)) isize = init_size


      if(.not.tensor_mpi_buffer%is_initialized)then
         call tensor_mpi_buffer%init_(isize,r,c)
      else if(present(init_size))then
         call tensor_status_quit("ERROR(tensor_buffer_add): buffer is already initialized, the buffer size is unnecessary",33)
      endif

      tsize = n * datatype

      b = tensor_mpi_buffer%offset
      e = b + tsize - 1

      if( rank == r )then
         !SENDER WRITE DATA TO BUFFER, IF BUFFER TOO SMALL, INCREASE
         if( e > tensor_mpi_buffer%nbytes ) then
            call tensor_mpi_buffer%incr_(tsize)
         endif

         tensor_mpi_buffer%buffer(b:e) = transfer(buffer,tensor_mpi_buffer%buffer(b:e))

      else
         !RECEIVER EXTRACT DATA FROM BUFFER, BUFFER CANNOT BE TOO SMALL
         if( e > tensor_mpi_buffer%nbytes )  call tensor_status_quit("ERROR(tensor_buffer_add): error in reading the message, &
            &possibly wrong order of arguments for the receiver",33)

         buffer = transfer(tensor_mpi_buffer%buffer(b:e),buffer)

      endif

      tensor_mpi_buffer%offset = e + 1

      if(fin) then
         call tensor_mpi_buffer%free_()
      endif
