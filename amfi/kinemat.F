      Subroutine kinemat(L,ndim,evtkin,type1,type2,Energy,HERMIT)
      implicit double precision (a-h,o-z) 
      parameter (fine=7.29735308D-03) !TO_BE_CHECKED    
cbs   at least it's identical with Odd's valuE
      parameter (speed=1d0/fine)     
      parameter (speed2=speed*speed) 
      parameter (speed4=speed2*speed2) 
cbs   this routine generates the kinematic A-factors=dsqrt((E+mc^2)/(2E))  
cbs   (type1) and   c*A/(E+mc^2) (type2)
cbs   The c in the second kinematic factor comes from Jan Almloef and 
cbs   Odd Gropen in Rev in Comp.Chem. 8(1996)
      dimension evtkin(*),type1(*),type2(*),Energy(*)  
      logical*4 HERMIT  
c     E= dsqrt(p**2 c**2 + m**2 c**4) 
c     p**2= 2*m*TKIN    
c     with m = 1 
      do Irun=1,ndim
      if (evtkin(Irun).lt.0) stop 'strange kinetic energy '
      Energy(Irun)=(evtkin(Irun)+evtkin(Irun))*speed2+speed4 
      enddo
      do Irun=1,ndim
      Energy(Irun)=dsqrt(energy(irun))
      enddo
      do Irun=1,ndim
!     dsqrt((E+mc^2)/(2E)):
      type1(Irun)=dsqrt(0.5d0*(1d0+speed2/Energy(Irun)))
      enddo
!      c*A/(E+mc^2) 
      do Irun=1,ndim
      type2(Irun)=speed*type1(Irun)/(Energy(Irun)+speed2)
      enddo
      if (HERMIT) then  ! to get a final factor 4 c**2 
	      do Irun=1,ndim
  	      type2(Irun)=2*speed*type2(Irun)
   	      enddo
      endif 
c     write(6,*) ' ' 
c     write(6,'(A,I3)') 'Kinematic factors for L=',L
c     write(6,*) ' ' 
c     write(6,*) 
c    *'Kinetic energy | Energy |  Kinematic factors A and A/(E+m):' 
c     do Irun=1,ndim
c     write(6,'(2(E20.14,X),F18.14,X,F18.14)')  evtkin(Irun),
c    *energy(irun)-speed2,type1(Irun),type2(irun)  
c     enddo
C     write(6,*) ' ' 
      return 
      end
