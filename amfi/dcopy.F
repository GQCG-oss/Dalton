      subroutine dcopy(N,X,INCX,Y,INCY)
cbs   BLAS-substitute  
      Double precision X(*),Y(*)
      if (N.le.0) return
      IX=1
      IY=1
      IF(INCX.lt.0) IX=1-(N-1)*INCX
      IF(INCY.lt.0) IY=1-(N-1)*INCY
      do I=1,N
      Y(IY)=X(IX)
      IX=IX+INCX
      IY=IY+INCY
      enddo
      return 
      end
