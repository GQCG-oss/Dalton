      n1 = size(p2,1,kind=tensor_long_int)
      n2 = size(p2,2,kind=tensor_long_int)

      n = n1 * n2

      call c_f_pointer(c_loc(p2(1,1)),p,[n])

      call free(p,idx,stat=stat)
      p2 => null()
