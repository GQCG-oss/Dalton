      subroutine incrpow(max,new,quot,quotpow,
     *nprim1,nprim2,nprim3,nprim4)
cbs   generates some powers of for the prefactors of cfunct(X)
cbs   look out for details there and in initfrac   
      implicit double precision (a-h,o-z)
      include 'para.h' 
      dimension quotpow(nprim1,nprim2,
     *nprim3,nprim4,*),
     *quot(nprim1,nprim2,nprim3,nprim4) 
c     write(6,*) 'incrpow was called ',new
      if ((new*nprim1*nprim2*nprim3*nprim4)
     *.gt.mxfrac) then 
      write(6,*) 'error in incrpow: '
      write(6,*) 'increase mxfrac to at least ',
     *new*nprim1*nprim2*nprim3*nprim4
      stop 
      endif 
      do irun=max+1,new
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      quotpow(irun1,irun2,irun3,irun4,irun)=
     *quotpow(irun1,irun2,irun3,irun4,irun-1)*
     *quot(irun1,irun2,irun3,irun4)
      enddo
      enddo
      enddo
      enddo
      enddo
      max=new
      return 
      end 
