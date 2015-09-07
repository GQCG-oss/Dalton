
     integer(kind=tensor_long_int) :: rank8
     integer(kind=tensor_long_int) :: globaltilenumber8
     integer(kind=tensor_long_int) :: pos8
     integer(kind=tensor_long_int) :: idx8
     integer(kind=tensor_long_int) :: widx8
     integer :: nnod

     globaltilenumber8 = globaltilenumber

     call get_residence_of_tile_basic(arr,globaltilenumber8,rank8,pos8,idx8,widx8)

     rank = rank8
     pos  = pos8
     idx  = idx8
     widx = widx8
