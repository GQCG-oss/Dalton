integer(kind=tensor_long_int) :: i,cind,ndim
cind=a
do i=1,modes
   ndim = dims(i)
   inds(i)=mod(cind-1,ndim)+1
   cind=(cind-inds(i))/dims(i) + 1
enddo
