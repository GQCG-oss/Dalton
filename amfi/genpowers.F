      subroutine genpowers(Lhigh)   
      implicit double precision (a-h,o-z) 
      include 'param.h' 
      include 'dofuc.h'
cbs   set some often used powers of exponents 
      do L2=0,Lhigh 
      do L1=0,L2    
      do irun1=1,nprimit(L1)
      do irun2=1,nprimit(L2)
      powexp(irun1,irun2,L1,L2,0)=1d0
      enddo
      enddo
      enddo
      enddo
      do L2=0,Lhigh 
      do L1=0,L2    
      do Lrun=1,(L1+L2+5)
      do irun2=1,nprimit(L2)
      do irun1=1,nprimit(L1)
      fact=dsqrt(0.5d0*(exponents(irun1,L1)+exponents(irun2,L2)))
      powexp(irun1,irun2,L1,L2,Lrun)= powexp(irun1,irun2,L1,L2,Lrun-1)*
     *fact
      enddo
      enddo
      enddo
      enddo
      enddo
cbs   generate coulovlp = overlap for normalized functions, but sometimes
cbs   with shifted l-values
      do l2=0,lhigh
      do incl2=-1,1
         if (l2+incl2.ge.0) then  ! do not lower l for s-functions 
         n2=l2+incl2+1
         df2=1d0/dsqrt(df(n2+n2-1))
         do l1=0,l2
         do incl1=-1,1
         if (l1+incl1.ge.0) then ! do not lower l for s-functions 
         n1=l1+incl1+1
         df1=1d0/dsqrt(df(n1+n1-1))
         df12=df(n1+n2-1)
         do iprim2=1,nprimit(l2)
         fact2=dsqrt(powexp(iprim2,iprim2,l2,l2,n2+n2+1))
         factor=fact2*df1*df2*df12
         do iprim1=1,nprimit(l1)
         fact1=dsqrt(powexp(iprim1,iprim1,l1,l1,n1+n1+1))
         coulovlp(iprim1,iprim2,incl1,incl2,l1,l2)= 
     *   fact1*factor/powexp(iprim1,iprim2,l1,l2,n1+n2+1) 
c        write(6,*) 'l1,l2,incl1,incl2,iprim1,iprim2,fact1,factor',
c    *   ',powexp(iprim1,iprim2,l1,l2,n1+n2+1)',
c    *   l1,l2,incl1,incl2,iprim1,iprim2,fact1,factor,
c    *   powexp(iprim1,iprim2,l1,l2,n1+n2+1)
         enddo
         enddo
         endif
         enddo
         enddo
         endif
      enddo
      enddo
      return 
      end 
 
    
