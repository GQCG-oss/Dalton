      integer(kind=tensor_long_int) :: vector_size
      integer :: i
      if (present(mode))then
        if((arr%mode/=0 .and. arr%mode/=mode).or.arr%mode==0)then
          print *,"mode",mode,"arr%mode",arr%mode      
          call lsquit("wrong use of tensor_set_tdim",lspdm_errout)
        else
          arr%mode=mode
        endif
      endif
      vector_size = int(arr%mode*tensor_standard_int,kind=tensor_long_int)
      if (associated(arr%tdim))then
!$OMP CRITICAL
        tensor_counter_aux_f_mem = tensor_counter_aux_f_mem + vector_size
        tensor_counter_memory_in_use    = tensor_counter_memory_in_use  - vector_size
!$OMP END CRITICAL
        call tensor_free_mem(arr%tdim)
      endif
      if(.not.associated(arr%tdim))then
        call tensor_alloc_mem(arr%tdim,arr%mode)
!$OMP CRITICAL
        tensor_counter_aux_a_mem = tensor_counter_aux_a_mem + vector_size
        tensor_counter_memory_in_use  = tensor_counter_memory_in_use  + vector_size
!$OMP END CRITICAL
      endif
      do i=1,arr%mode
        arr%tdim(i)=tdims(i)
      enddo
