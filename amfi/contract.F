      subroutine contract(
     *coeffs1, !(nprim(1),ncont(1)) modified contraction coefficients
     *coeffs2, !(nprim(2),ncont(2)) modified contraction coefficients
     *coeffs3, !(nprim(3),ncont(3)) modified contraction coefficients
     *coeffs4, !(nprim(4),ncont(4)) modified contraction coefficients
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
cbs  array one contains at the beginning the uncontracted integrals 
     *arr1,  ! array of size (nprim(1)*nprim(2)*nprim(3)*nprim(4))
     *arr2   ! array of size (nprim(1)*nprim(2)*nprim(3)*nprim(4))
     *)
      implicit double precision (a-h,o-z)
      dimension coeffs1(*),coeffs2(*),coeffs3(*),coeffs4(*),   
     *arr1(*),arr2(*),ncont(4),nprim(4),nolds(4),nnew(4)
C
cbs   makes four indextransformations in a row....
cbs   try to find out, which indices should be transformed first...
c
c     write(6,*) 'begin contract ' 
      ratio1=dfloat(nprim(1))/dfloat(ncont(1))
      ratio2=dfloat(nprim(2))/dfloat(ncont(2))
      ratio3=dfloat(nprim(3))/dfloat(ncont(3))
      ratio4=dfloat(nprim(4))/dfloat(ncont(4))
      do IBM=1,4
      nolds(IBM)=nprim(IBM)  
      nnew(IBM)=nprim(IBM)  
      enddo
c    test
c     nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
c     write(6,*) nprod,'integrals before the contractions' 
c     write(6,'(6E12.4)') (arr1(I),I=1,nprod) 
c     write(6,*) ' '
cbs   determine first, second,third and last index
cbs   determine the first
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      ifirst=1
      ratio1=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
      else if (xmax.eq.ratio2) then 
      ifirst=2
      ratio2=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio3) then 
      ifirst=3
      ratio3=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio4) then 
      ifirst=4
      ratio4=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      endif
      nolds(ifirst)=nnew(ifirst)
c     write(6,*) 'ifirst = ',ifirst   
cbs   determine the second   
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      isec=1
      ratio1=0
      nnew(isec)=ncont(isec)
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio2) then 
      isec=2
      ratio2=0
      nnew(isec)=ncont(isec)
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio3) then 
      isec=3
      ratio3=0
      nnew(isec)=ncont(isec)
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio4) then 
      isec=4
      ratio4=0
      nnew(isec)=ncont(isec)
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      endif
      nolds(isec)=nnew(isec)
c     write(6,*) 'isecond = ',isec   
cbs   determine the third    
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      ithird=1
      ratio1=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio2) then 
      ithird=2
      ratio2=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio3) then 
      ithird=3
      ratio3=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio4) then 
      ithird=4
      ratio4=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      endif
      nolds(ithird)=nnew(ithird)
c     write(6,*) 'ithird = ',ithird    
cbs   determine the last    
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      ifourth=1
      ratio1=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio2) then 
      ifourth=2
      ratio2=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio3) then 
      ifourth=3
      ratio3=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio4) then 
      ifourth=4
      ratio4=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      endif
c     write(6,*) 'ifourth= ',ifourth   
cbs   contracted integrals are now on 
cbs   arr1(ncont1,ncont2,ncont3,ncont4)   
c     write(6,*) 'end contract ' 
      return 
      end
