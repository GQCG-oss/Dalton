      subroutine writecont2(contintgrls,ncont1,ncont2,ncont3,ncont4)
      implicit double precision (a-h,o-z) 
      dimension contintgrls(ncont1,ncont2,ncont3,ncont4)
      do irun1=1,ncont1
      do irun2=1,ncont2
      do irun3=1,ncont3
      do irun4=1,ncont4
      write(6,'(4I3,E20.14)') irun1,irun2,irun3,irun4,
     *contintgrls(irun1,irun2,irun3,irun4)
      enddo
      enddo
      enddo
      enddo
      return
      end
