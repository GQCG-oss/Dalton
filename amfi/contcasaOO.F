      subroutine contcasaOO(l1,l2,l3,l4,nstart,writetwo)
cbs   contraction for powers (+2)  with alpha1*alpha3
cbs   other-orbit term 
cbs   use averaged integrals by interchanging kinematic factors   
cbs   this is case a in the documentation 
      implicit double precision (a-h,o-z)
      include 'param.h' 
      logical*4 writetwo
      dimension ncont(4),nprim(4) 
c     write(6,*) 'begin contcasa '
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp2(l1)), !A *alpha                                       
     *contrarray(iaddtyp3(l2)), !A/E+m                                         
     *contrarray(iaddtyp4(l3)), !A/E+m *alpha     
     *contrarray(iaddtyp1(l4)), !A                                       
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=0.25d0*scratch1(irun)
      enddo 
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)), 
     *contrarray(iaddtyp3(l2)), 
     *contrarray(iaddtyp2(l3)), 
     *contrarray(iaddtyp1(l4)), 
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+0.25d0*
     *scratch1(irun)
      enddo 
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp2(l1)),
     *contrarray(iaddtyp1(l2)),
     *contrarray(iaddtyp4(l3)),
     *contrarray(iaddtyp3(l4)),
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+0.25d0*
     *scratch1(irun)
      enddo 
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)), 
     *contrarray(iaddtyp1(l2)),
     *contrarray(iaddtyp2(l3)),
     *contrarray(iaddtyp3(l4)),
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+0.25d0*
     *scratch1(irun)
      enddo 
C
C
C
      if (writetwo) then 
      write(6,*) 'contracted radial integrals type A other-orbit'
       call writecont2(cont4OO(nstart),ncont(1),ncont(2),
     *ncont(3),ncont(4))   
      endif
c     write(6,*) 'contcasaOO ',scratch1(1)
c     write(6,*) 'end contcasaOO '
      return 
      end
