      subroutine dzero(ianz,vec)
      implicit double precision (a-h,o-z)
cbs   cleans up the vector vec
      dimension vec(*)
      do irun=1,ianz
      vec(irun)=0d0
      enddo 
      return
      end
