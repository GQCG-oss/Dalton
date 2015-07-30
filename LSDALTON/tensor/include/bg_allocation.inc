      integer :: alloc_stat 
      integer(kind=tensor_long_int) :: b

      alloc_stat = 0

      !bytes are the number of elements times the size
      b = n * counters_bg(idx)%size_

      call mem_pseudo_alloc(p,n)

      !$OMP CRITICAL
      counters_bg(idx)%curr_ = counters_bg(idx)%curr_ + b
      counters_bg(idx)%high_ = max(counters_bg(idx)%high_,counters_bg(idx)%curr_)
      tensor_counter_max_bg_mem = max(tensor_counter_max_bg_mem,tensor_get_total_bg_mem())
      !$OMP END CRITICAL
         

      if( alloc_stat /= 0 ) then
         if(.not. present(stat))then
            call tensor_status_quit("ERROR(tensor_allocate_generic_bg): allocation failed",alloc_stat)
         else
            stat = alloc_stat
         endif
      endif
