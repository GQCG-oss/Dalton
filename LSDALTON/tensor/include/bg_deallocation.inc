      integer :: alloc_stat 
      integer(tensor_long_int) :: b,n
      logical :: md, last_assoc, del_assoc
      integer :: i
      type(c_ptr) :: loc1,loc2

      alloc_stat = 0
      n          = size(p)

      !bytes are the number of elements times the size
      b = n * counters_bg(idx)%size_

      !call mem_pseudo_dealloc(p,mark_deleted=.true.)

      md=.true.
      !no need to mark_deleted if it is the correct last element in the buffer
      last_assoc = c_associated(buf%c_addr(buf%n-1),c_loc(p(1)))
      del_assoc  = c_associated(buf%f_mdel,         c_loc(p(1)))
      if(last_assoc) md = .false.

      if(md)then
         if(.not.buf%l_mdel)then
            FindPos:do i=1,buf%n
               loc1 = c_loc(p(1))
               loc2 = c_loc(buf%p(buf%f_addr(i)))
               if( c_associated(loc1,loc2))then
                  buf%f_mdel = c_loc(buf%p(buf%f_addr(i+1)))
                  exit FindPos
               endif
            enddo FindPos
         endif
         buf%l_mdel = .true.
         buf%n_mdel = buf%n_mdel + 1
      else
         buf%n = buf%n - 1
      endif

      if(buf%n < 0)then
         call tensor_status_quit("ERROR(pseudo_dealloc_generic): programming error, more&
         & pointers freed than associated",-1)
      endif

      if(.not.md.and..not.last_assoc)then
         call tensor_status_quit("ERROR(pseudo_dealloc_generic): wrong sequence of &
         &deallocating, make sure you dealloc in the opposite seqence as allocating, &
         &otherwise you will corrupt your data",3737)
      endif

      n = size(p,kind=tensor_long_int)

      if (buf%offset-n < 0)then
         call  tensor_status_quit("ERROR(pseudo_dealloc_generic): more freed than allocated",33883)
      endif


      if(md)then
         buf%e_mdel(buf%n_mdel) = n
         buf%c_mdel(buf%n_mdel) = c_loc(p(1))
      else
         buf%offset = buf%offset-n
      end if

      p => null()
      buf%c_addr(buf%n) = c_null_ptr

      if(del_assoc)then
         call buf%clear_md()
      endif


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
