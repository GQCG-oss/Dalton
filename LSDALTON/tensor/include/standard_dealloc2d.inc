      n1 = size(p2,1,kind=tensor_long_int)
      n2 = size(p2,2,kind=tensor_long_int)

      n = n1 * n2

!#ifdef VAR_PTR_RESHAPE
      p(1:n) => p2(:,:)
!#else
!      call c_f_pointer(c_loc(p2(1,1)),p,[n])
!#endif
      call free(p,idx,stat=stat)
      p2 => null()
