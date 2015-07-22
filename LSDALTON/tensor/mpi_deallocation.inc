      integer(kind=tensor_mpi_kind)  :: alloc_stat 
      integer(kind=MPI_ADDRESS_KIND) :: b,lb,sze
      integer(kind=tensor_mpi_kind)  :: info
      integer(tensor_long_int) :: n

      if(  .not. c_associated( c, c_loc(p) )  )then
         call tensor_status_quit("ERROR(tensor_free_generic_mpi): wrong c/p combination",440)
      endif

      info = MPI_INFO_NULL

      alloc_stat = 0
      n          = size(p)

      !bytes are the number of elements times the size
      b = int(n * counters(idx)%size_,kind=MPI_ADDRESS_KIND)

      call MPI_FREE_MEM(p,alloc_stat)

      !$OMP CRITICAL
      counters(idx)%curr_ = counters(idx)%curr_ - b
      !$OMP END CRITICAL
         

      if( alloc_stat /= 0 ) then
         if(.not. present(stat))then
            call tensor_status_quit("ERROR(tensor_free_generic_mpi): allocation failed",alloc_stat)
         else
            stat = alloc_stat
         endif
      endif

      p => null()
      c = c_null_ptr
