      subroutine daxpint(from,to,fact,ndim1,ndim2,ndim3,ndim4)
      implicit double precision (a-h,o-z)
cbs   subroutine similar to daxpy with interchange of two indices 
cbs   change from physicists notation to chemists notaion 
cbs   to(i,j,k,l)=to(i,j,k,l)+fact*from(i,k,j,l) 
      dimension from(ndim1,ndim2,ndim3,ndim4),
     *to(ndim1,ndim3,ndim2,ndim4)
      if (fact.eq.0d0) return 
      do irun4=1,ndim4
      do irun3=1,ndim3
      do irun2=1,ndim2
      do irun1=1,ndim1
      to(irun1,irun3,irun2,irun4)=to(irun1,irun3,irun2,irun4)+
     *fact*from(irun1,irun2,irun3,irun4)  
      enddo
      enddo
      enddo
      enddo
      return   
      end 
