      subroutine initfrac(nprimit1,nprimit2,
     *nprimit3,nprimit4,
     *quot1,quotpow1,quot2,quotpow2,expo1,expo2,
     *expo3,expo4)
cbs   initialize some arrays with factors  needed for cfunct(x) 
      implicit double precision(a-h,o-z) 
      dimension expo1(*),expo2(*),expo3(*),expo4(*), 
     *quot1(nprimit1,nprimit2,nprimit3,nprimit4),
     *quot2(nprimit1,nprimit2,nprimit3,nprimit4),
     *quotpow1(nprimit1,nprimit2,nprimit3,nprimit4,*),
     *quotpow2(nprimit1,nprimit2,nprimit3,nprimit4,*)
      do irun4=1,nprimit4 
      do irun3=1,nprimit3 
      do irun2=1,nprimit2 
      	sum24=expo2(irun2)+expo4(irun4)
      		do irun1=1,nprimit1 
      		quot1(irun1,irun2,irun3,irun4)= 
     *		1d0/(1d0+(expo1(irun1)+expo3(irun3))/
     *		sum24)                             
      		enddo
      enddo
      enddo
      enddo
      do irun4=1,nprimit4 
      do irun3=1,nprimit3 
      do irun2=1,nprimit2 
      	sum24=expo2(irun2)+expo4(irun4)
      		do irun1=1,nprimit1 
      		quot2(irun1,irun2,irun3,irun4)= 
     *		1d0/(1d0+sum24/
     *		(expo1(irun1)+expo3(irun3)))
      		enddo
      enddo
      enddo
      enddo
      do irun4=1,nprimit4 
      do irun3=1,nprimit3 
      do irun2=1,nprimit2 
      do irun1=1,nprimit1 
      	quotpow1(irun1,irun2,irun3,irun4,1)= 
     *	dsqrt(quot1(irun1,irun2,irun3,irun4))
      enddo
      enddo
      enddo
      enddo
      do irun4=1,nprimit4 
      do irun3=1,nprimit3 
      do irun2=1,nprimit2 
      do irun1=1,nprimit1 
      	quotpow2(irun1,irun2,irun3,irun4,1)= 
     *	dsqrt(quot2(irun1,irun2,irun3,irun4))
      enddo
      enddo
      enddo
      enddo
      return
      end
