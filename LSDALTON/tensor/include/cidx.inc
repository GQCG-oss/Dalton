integer :: a,b,i,j,cdim,cd(modes-2)
select case(modes)
case(2)
  a=inds(1)+(inds(2)-1)*dims(1)
case(3)
  cd(1) = dims(1)*dims(2)
  a=inds(1)+(inds(2)-1)*dims(1)+(inds(3)-1)*cd(1)
case(4)
  cd(1) = dims(1)*dims(2)
  cd(2) = dims(1)*dims(2)*dims(3)
  a=inds(1)+(inds(2)-1)*dims(1)+(inds(3)-1)*cd(1)+(inds(4)-1)*cd(2)
case default
  a=1
  do i=1,modes
    cdim=1
    do j=i-1,1,-1
      cdim=cdim*dims(j)
    enddo
    a=a+(inds(i)-1)*cdim
  enddo
end select
