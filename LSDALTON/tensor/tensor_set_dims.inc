      integer(kind=tensor_long_int) :: vector_size
      integer :: i
      if (present(mode))then
        if((arr%mode/=0 .and. arr%mode/=mode).or.arr%mode==0)then
          print *,"mode",mode,"arr%mode",arr%mode      
          call lsquit("wrong use of tensor_set_dims",lspdm_errout)
        else
          arr%mode=mode
        endif
      endif
      if (associated(arr%dims))then
!$OMP CRITICAL
        vector_size = int(size(arr%dims)*tensor_standard_int,kind=tensor_long_int)
        tensor_counter_aux_f_mem = tensor_counter_aux_f_mem + vector_size
        tensor_counter_memory_in_use    = tensor_counter_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%dims)
      endif
      if(.not.associated(arr%dims))then
        call mem_alloc(arr%dims,arr%mode)
!$OMP CRITICAL
        vector_size = int(size(arr%dims)*tensor_standard_int,kind=tensor_long_int)
        tensor_counter_aux_a_mem = tensor_counter_aux_a_mem + vector_size
        tensor_counter_memory_in_use  = tensor_counter_memory_in_use  + vector_size
!$OMP END CRITICAL
      endif
      do i=1,arr%mode
        arr%dims(i)=dims(i)
      enddo