      subroutine contcasb2SO(l1,l2,l3,l4,nstart,writetwo)
cbs   contraction for powers (0)  with alpha3
cbs   this is one of the cases b in the documentation 
      implicit double precision (a-h,o-z)
      include 'param.h' 
      logical*4 writetwo
      dimension ncont(4),nprim(4) 
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp4(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      call dcopy(nprod,scratch1(1),1,cont4SO(nstart),1)
      if (writetwo) then
      write(6,*) 'contracted radial integrals type B2 same-orbit'
      call writecont2(cont4SO(nstart),ncont(1),ncont(2),
     *ncont(3),ncont(4))
      endif  
      return 
      end
