      integer :: alloc_stat 
      integer(kind=tensor_long_int) :: b

      alloc_stat = 0

      !bytes are the number of elements times the size
      b = n * counters_bg(idx)%size_

      !call mem_pseudo_alloc(p,n)
      if (buf%offset+n > buf%nmax)then
         print *,"Buffer Space (#elements):",buf%nmax," Used:",buf%offset," Requested:",n
         call tensor_status_quit("ERROR(tensor_bg_alloc_generic): more requested than available",3773)
      endif

      p => buf%p(buf%offset+1:buf%offset+n)

      buf%f_addr(buf%n) = buf%offset+1
      buf%c_addr(buf%n) = c_loc(p(1))

      buf%offset = buf%offset+n
      buf%max_usage = max(buf%max_usage,buf%offset) 


      buf%n = buf%n + 1

      if(buf%n > tensor_max_bg_pointers)then
         call tensor_status_quit("ERROR(tensor_bg_alloc_generic): more pointers associated then currently supported, &
            &please change max_n_pointers in LSDALTON/lsutil/background_buffer.F90",3773)
      endif

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
