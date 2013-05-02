program ppsphmat
use print_moorb_grid_mod
use grid_utilities_module

real(4), allocatable           :: sphmat(:)
integer :: J,angmom,P,PN,spSIZE  

   J=5

   angmom = J-1
   P = (angmom+1)*(angmom+2)/2
   PN = 2*angmom+1
   spSIZE=PN*P
   allocate(SPHMAT(spSIZE))
   SPHMAT = 0_4
   CALL BUILD_SPHMAT(angmom,spSIZE,SPHMAT)


   write(6,'(5E16.8)') (SPHMAT(i),i=1,spSIZE)

   deallocate(sphmat)
endprogram
