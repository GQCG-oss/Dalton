      subroutine inidf
cbs   initializes the df on common block  with double facultatives
      implicit double precision (a-h,o-z)
      include 'param.h'
      include 'dofuc.h'
c     write(6,*) 'begin inidf'
      df(0)=1.d0
      df(1)=1.d0
      do irun=2,ndfmx
      df(irun)=dfloat(irun)*df(irun-2)
      enddo
      do jbm=0,ndfmx-1  
      do ibm=jbm,ndfmx   
      dffrac(ibm,jbm)=df(ibm)/df(jbm) 
      enddo 
      enddo 
      do jbm=1,ndfmx  
      do ibm=0,jbm-1       
      dffrac(ibm,jbm)=1d0/dffrac(jbm,ibm)   
      enddo 
      enddo 
c     do ibm=0,ndfmx
c     write(6,*) ibm,df(ibm)
c     enddo
c     write(6,*) 'end inidf'
      return 
      end
