      type(c_ptr), intent(inout) :: c
      integer(kind=tensor_standard_int),intent(in) :: idx
      integer(kind=tensor_long_int), intent(in)    :: n
      integer, intent(out), optional               :: stat
      integer(kind=tensor_mpi_kind)  :: alloc_stat 
      integer(kind=MPI_ADDRESS_KIND) :: b
      integer(kind=tensor_mpi_kind)  :: info

      info = MPI_INFO_NULL

      alloc_stat = 0

      !bytes are the number of elements times the size
      b = int(n * counters(idx)%size_,kind=MPI_ADDRESS_KIND)

      call MPI_ALLOC_MEM(b,info,c,alloc_stat)

      !$OMP CRITICAL
      counters(idx)%curr_ = counters(idx)%curr_ + b
      counters(idx)%high_ = max(counters(idx)%high_,counters(idx)%curr_)
      tensor_counter_max_hp_mem = max(tensor_counter_max_hp_mem,tensor_get_total_heap_mem())
      !$OMP END CRITICAL
         

      if( alloc_stat /= 0 ) then
         if(.not. present(stat))then
            call tensor_status_quit("ERROR(tensor_allocate_generic_mpi): allocation failed",alloc_stat)
         else
            stat = alloc_stat
         endif
      endif

      call c_f_pointer(c,p,[n])
