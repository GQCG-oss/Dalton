      subroutine genovlp(Lhigh)
      implicit double precision (a-h,o-z) 
cbs   generates overlap of normalized  primitives. 
      include 'param.h' 
      dimension evecinv(MxprimL,MxprimL)
      do L=0,Lhigh 
      	do Jrun=1,nprimit(L)
      	do Irun=1,nprimit(L)
        normovlp(Irun,Jrun,L)=coulovlp(irun,jrun,0,0,
     *  L,L)   
      	enddo 
      	enddo 
cbs   invert the matrix, not very elegant, but sufficient
      ipnt=0
      do jrun=1,nprimit(L)
      do irun=1,jrun
      ipnt=ipnt+1
      scratchinv(ipnt)=normovlp(irun,jrun,L)
      enddo
      enddo
c     write(6,*)  'scratchinv '
c     write(6,*) 
c     ie=0
c     do jrun=1,nprimit(L)
c     ia=ie+1
c     ie=ie+jrun
c     write(6,'(4(E18.12,X))') 
c    *(scratchinv(iii),iii=ia,ie)
c     write(6,*)
c     enddo
      do Jrun=1,nprimit(L)
      do Irun=1,nprimit(L)
      evecinv(Irun,Jrun)=0d0                                
      enddo 
      enddo 
      do Jrun=1,nprimit(L)
      evecinv(jrun,jrun)=1d0
      enddo 
      call jacobi(scratchinv,evecinv,nprimit(L),MxprimL) 
c     write(6,*) 'eigenvalues of overlap '
c     write(6,'(4E20.14)') (scratchinv((irun*irun+irun)/2),
c    *irun=1,nprimit(L))
      do irun=1,nprimit(L)
      eval(irun)=dsqrt(scratchinv((irun*irun+irun)/2))
      enddo
cbs   ensure normalization of the vectors.
      do IRUN=1,nprimit(L)
      fact=0d0
      do JRUN=1,nprimit(L)
      fact=fact+evecinv(JRUN,IRUN)*evecinv(JRUN,IRUN)
      enddo
      fact=1d0/dsqrt(fact)
      do JRUN=1,nprimit(L)
      evecinv(JRUN,IRUN)=fact*evecinv(JRUN,IRUN)
      enddo
      enddo  
cbs   now generate rootOVLP           
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      rootOVLP(irun,jrun,l)=0d0
      enddo  
      enddo  
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      rootOVLP(irun,jrun,L)=rootOVLP(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo  
      enddo  
      enddo  
cbs   now generate rootOVLPinv           
      do irun=1,nprimit(L)
      eval(irun)=1d0/eval(irun)                            
      enddo
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      rootOVLPinv(irun,jrun,l)=0d0
      enddo  
      enddo  
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      rootOVLPinv(irun,jrun,L)=rootOVLPinv(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo  
      enddo  
      enddo  
cbs   now generate OVLPinv           
      do irun=1,nprimit(L)
      eval(irun)=eval(irun)*eval(irun)
      enddo
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      OVLPinv(irun,jrun,l)=0d0
      enddo  
      enddo  
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      OVLPinv(irun,jrun,L)=OVLPinv(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo  
      enddo  
      enddo  
      enddo 
      return 
      end   
