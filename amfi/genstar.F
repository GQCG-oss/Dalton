      subroutine genstar(Lhigh) 
      implicit double precision (a-h,o-z)
cbs   purpose: generate start adresses of contraction coeffs on 
cbs   contrarray for the different L-Blocks 
      include 'param.h' 
      istart=1
      do L=0,Lhigh
      inc=nprimit(L)*ncontrac(L)
      iaddori(L)=istart
      istart=istart+inc  
      iaddtyp1(L)=istart 
      istart=istart+inc  
      iaddtyp2(L)=istart 
      istart=istart+inc  
      iaddtyp3(L)=istart 
      istart=istart+inc  
      iaddtyp4(L)=istart 
      istart=istart+inc  
      enddo 
      return 
      end  
