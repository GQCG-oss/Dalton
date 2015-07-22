      integer(kind=tensor_standard_int), intent(in) :: idx
      integer(kind=tensor_long_int), intent(in)     :: n
      integer, intent(out), optional                :: stat
      integer :: alloc_stat 
      integer(kind=tensor_long_int) :: b

      alloc_stat = 0

      !bytes are the number of elements times the size
      b = n * counters(idx)%size_

      allocate(p(n),stat=alloc_stat)
      !$OMP CRITICAL
      counters(idx)%curr_ = counters(idx)%curr_ + b
      counters(idx)%high_ = max(counters(idx)%high_,counters(idx)%curr_)
      tensor_counter_max_hp_mem = max(tensor_counter_max_hp_mem,tensor_get_total_heap_mem())
      !$OMP END CRITICAL
         

      if( alloc_stat /= 0 ) then
         if(.not. present(stat))then
            call tensor_status_quit("ERROR(tensor_allocate_generic): allocation failed",alloc_stat)
         else
            stat = alloc_stat
         endif
      endif
