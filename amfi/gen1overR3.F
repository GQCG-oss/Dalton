      subroutine gen1overR3(Lhigh)
      implicit double precision (a-h,o-z) 
cbs   generates the radial integrals  for the one electron spin orbit integrals   
cbs   taken the 1/r**3 formula from the documentation and included additional 
cbs   factors for normalization   
      include 'param.h' 
      include 'dofuc.h'
      common /pifact/ pi,root8ovpi 
      do L=1,Lhigh    
      icount=0
      do iprim2=1,nprimit(L)
      alpha2=exponents(iprim2,L)
      do iprim1=1,iprim2    
      alpha1=exponents(iprim1,L)
      icount=icount+1
      oneoverR3(icount,L)=dsqrt(2d0/pi)*
     *(df(L+L-2)*2**(L+3)*
     *(alpha1*alpha2)**(0.25d0*
     *(L+L+3)))/((alpha1+alpha2)**L*df(L+L+1))
      enddo
      enddo
      enddo 
      return 
      end   
