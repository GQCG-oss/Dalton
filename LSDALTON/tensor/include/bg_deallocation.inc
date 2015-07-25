      integer(kind=tensor_standard_int), intent(in) :: idx
      integer, intent(out), optional                :: stat
      integer :: alloc_stat 
      integer(tensor_long_int) :: b,n

      alloc_stat = 0
      n          = size(p)

      !bytes are the number of elements times the size
      b = n * counters_bg(idx)%size_

      call mem_pseudo_dealloc(p,mark_deleted=.true.)

      !$OMP CRITICAL
      counters_bg(idx)%curr_ = counters_bg(idx)%curr_ - b
      !$OMP END CRITICAL
         

      if( alloc_stat /= 0 ) then
         if(.not. present(stat))then
            call tensor_status_quit("ERROR(tensor_free_generic_bg): deallocation failed",alloc_stat)
         else
            stat = alloc_stat
         endif
      endif
