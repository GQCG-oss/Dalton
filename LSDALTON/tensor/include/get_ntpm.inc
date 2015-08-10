     integer :: j
     if(present(ntiles)) ntiles = 1
     do j=1,mode
        ntpm(j)=dims(j)/tdim(j)
        if(mod(dims(j),tdim(j))>0)ntpm(j)=ntpm(j)+1
        if(present(ntiles)) ntiles = ntiles * ntpm(j)
     enddo
