      subroutine contone(L,oneoverr3,onecontr,Lmax,
     *contcoeff,nprim,ncont,MxcontL,dummy,
     *onecartx,onecartY,onecartZ,charge)
cbs   contracts one-electron integrals and multiplies with l,m-dependent 
cbs   factors for L-,L0,L+ 
      implicit double precision (a-h,o-z)
      dimension oneoverR3(*),
     *onecontr(MxcontL,MxcontL,-Lmax:Lmax,3),
     *contcoeff(nprim,ncont),dummy(ncont,ncont),
     *onecartx(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1)),  
     *onecarty(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1)),  
     *onecartz(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1))   
      ipnt(I,J)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
c     write(6,*) 'one-electron integrals L= ',L
cbs   first of all cleaning dummy and onecontr
      do jrun=1,ncont
      do irun=1,ncont
      dummy(irun,jrun)=0d0   
      enddo 
      enddo 
      iprod=3*(Lmax+lmax+1)*MxcontL*MxcontL
      call dzero(iprod,onecontr)
cbs   contract onto dummy 
      do icont2=1,ncont
      do icont1=1,ncont 
      do iprim2=1,nprim
      do iprim1=1,nprim
      dummy(icont1,icont2)=dummy(icont1,icont2)+
     *contcoeff(iprim1,icont1)*contcoeff(iprim2,icont2)*
     *oneoverR3(ipnt(iprim1,iprim2))
c     write(6,*) 'the uncontracted integral ', L,
c    *oneoverR3(ipnt(iprim1,iprim2)),contcoeff(iprim1,icont1),
c    *contcoeff(iprim2,icont2)
      enddo
      enddo
      enddo 
      enddo 
      do icont2=1,ncont
      do icont1=1,ncont 
      dummy(icont1,icont2)=dummy(icont1,icont2)*charge 
      enddo
      enddo
c     write(6,*) 
c     write(6,'(A,I3)') 'the contracted Z/r**3 integrals for L=',L
c     write(6,*) 
c     do icont2=1,ncont
c     write(6,'(4E20.14)') (dummy(icont1,icont2),icont1=1,icont2) 
c     write(6,*) 
c     enddo 
cbs   start to add l,m dependent factors    
      do M=-L,L 
      factormin=dsqrt(dfloat(L*L-M*M+L+M))
      factor0=dfloat(M)
      factorplus=dsqrt(dfloat(L*L-M*M+L-M))
      do irun=1,ncont 
      do jrun=1,ncont  
      onecontr(irun,jrun,M,1)=dummy(jrun,irun)*factormin  ! L-minus 
      enddo
      enddo
      do irun=1,ncont 
      do jrun=1,ncont  
      onecontr(irun,jrun,M,2)=dummy(jrun,irun)*factor0    ! L-0     
      enddo
      enddo
      do irun=1,ncont 
      do jrun=1,ncont  
      onecontr(irun,jrun,M,3)=dummy(jrun,irun)*factorplus ! L-plus  
      enddo
      enddo
      enddo
c     write(6,*) 'integrals in spherical basis(angular factor only)'
c     write(6,*) 'L-value ',L 
c     www=1d0/dummy(1,1)
c     do M=-L,L
c     write(6,'(I5,3E20.4)') M,onecontr(1,1,M,1)*www,
c    *onecontr(1,1,M,2)*www,onecontr(1,1,M,3)*www
c     enddo
cbs   make the final cartesian integrals 
      call cartoneX(L,Lmax,onecontr,ncont,
     *MxcontL,onecartX(1,1,1))
      call cartoneY(L,Lmax,onecontr,ncont,
     *MxcontL,onecartY(1,1,1))
      call cartoneZ(L,Lmax,onecontr,ncont,
     *MxcontL,onecartZ(1,1,1))
      return 
      end
