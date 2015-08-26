    nels=1
    do j=1, arr%mode
      if(((arr%dims(j)-(orig_addr(j)-1)*arr%tdim(j))/arr%tdim(j))>=1)then
        nels=nels*arr%tdim(j)
      else
        nels=nels*mod(arr%dims(j),arr%tdim(j))
      endif
    enddo
