      subroutine kindiag(TKIN,TKINTRIA,ndim,evec,eval,breit)
      implicit double precision (a-h,o-z)
cbs   determines eigenvectors and -values of TKIN  
      dimension tkin(ndim,ndim),
     *TKINTRIA((ndim*ndim+ndim)/2),eval(ndim),evec(ndim,ndim)
      logical*4 breit
cbs   move symmetric matrix to triangular matrix 
      itria=1
      do irun2=1,ndim
      do irun1=1,irun2 
      TKINTRIA(itria)=TKIN(irun1,irun2)
      itria=itria+1 
      enddo
      enddo
      do irun2=1,ndim
      do irun1=1,ndim
      evec(irun1,irun2)=0d0
      enddo
      enddo
      do irun1=1,ndim
      evec(irun1,irun1)=1d0
      enddo
cbs   now diagonalize  
c     write(6,*) 'calling eigen dimension ',ndim
c     ie=0
c     write(6,*) 'the matrix '
c     do irun=1,ndim
c     ia=ie+1
c     ie=ia+irun-1
c     write(6,'(6E14.6)') (TKINTRIA(III),III=ia,ie)
c     enddo
c           CALL EIGEN(TKINTRIA,evec,ndim,0,1)     
            CALL jacobi(TKINTRIA,evec,ndim,ndim)    
c     call dsyev('V','U',ndim,evec,ndim,eval,
c    *work,ndim*ndim*ndim+10,info)
c     if (info.ne.0) then 
c     write(6,*) 'Error at diagonalization in kindiag. error code ',info 
c     stop 
c     endif  
cbs   get the eigenvalues   
      do irun=1,ndim
      eval(irun)=TKINTRIA((irun*irun+irun)/2)
      enddo
c     write(6,*) 'the eigenvalues' 
c     write(6,'(4E20.14)') (eval(I),I=1,ndim)
      if (breit) then
      do irun=1,ndim
      eval(irun)=0d0 
      enddo
      endif 
cbs   ensure normalization of the vectors. 
      do IRUN=1,ndim
      fact=0d0
      do JRUN=1,ndim 
      fact=fact+evec(JRUN,IRUN)*evec(JRUN,IRUN) 
      enddo
      fact=1d0/dsqrt(fact)
      do JRUN=1,ndim
      evec(JRUN,IRUN)=fact*evec(JRUN,IRUN)
      enddo 
      enddo
      return   
      end   
