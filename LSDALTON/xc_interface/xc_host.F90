!> @file
!> Module containing main exchange-correlation integral driver, and routines to evaluate AOs and electron-densities
!> \brief 
!> \author T. Kjaergaard
!> \date 2009 
MODULE XCHOST
!use gridgenerationmodule
!use memory_handling
!use LSparameters
use precision
use TYPEDEF
use dft_type
use BUILDAOBATCH
! Notes 
! OLD NOTATION !   NEW NOTATION
!==========================================================
! NUCO         !   SHELLNPRIM   !nprimitives for the shell
! NHKT         !   SHELLANGMOM  !angmom for the shell
! NCENT        !   SHELL2ATOM   !Atom for the shell
! KMAX         !   MAXNSHELL    !maximum number of shells
! NHTYP        !   maxAngmom    !maximum angmom+1
! JSTRT        !   PRIEXPSTART  !index to start in PRIEXP for shell
!========================================================== 
CONTAINS
!> \brief wrapper exchange-correlation integral routine that build basinf.
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE XC_HOST_interface(LUPRI,IPRINT,SETTING,DMAT,NBAST,NDMAT,NGEODRV,DOLND)
IMPLICIT NONE
INTEGER,intent(in)     :: LUPRI,IPRINT,NBAST,NDMAT,NGEODRV
REAL(REALK),intent(in) :: DMAT(NBAST,NBAST,NDMAT)
TYPE(LSSETTING) :: SETTING
LOGICAL         :: DOLND !do london
!
TYPE(BASINF)  :: BAS
integer :: GRDONE,MAXREP,nCC,maxprim,iCC,iprim,i
logical :: DoSpherical,dopriexp
integer,pointer :: nPrim(:)
real(realk),pointer :: Coord(:,:),Charge(:),CC(:,:),priexp(:,:)
GRDONE=0
dopriexp = .TRUE.
CALL BUILD_BASINF(LUPRI,IPRINT,BAS,SETTING,GRDONE,dopriexp)

maxrep = 0
nCC = BAS%ushells
maxprim = 0
DO iCC = 1, nCC
   maxprim = MAX(maxprim,BAS%CC(iCC)%nrow)
ENDDO
call mem_alloc(CC,maxprim,nCC)
call mem_alloc(priexp,maxprim,nCC)
call mem_alloc(nPrim,nCC)
DO iCC = 1, nCC
   nPrim(iCC) = BAS%CC(iCC)%nrow
   do iprim = 1, BAS%CC(iCC)%nrow
      CC(iprim,iCC) = BAS%CC(iCC)%elms(iprim)
   enddo
   do iprim = 1, BAS%CC(iCC)%nrow
      priexp(iprim,iCC) = BAS%priexpM(iCC)%elms(iprim)
   enddo
ENDDO
call mem_alloc(Coord,3,BAS%nAtoms)
do I=1,BAS%nAtoms
   Coord(1,I) = BAS%X(I)
   Coord(2,I) = BAS%Y(I)
   Coord(3,I) = BAS%Z(I)
enddo
call mem_alloc(Charge,BAS%nAtoms)
do I=1,BAS%nAtoms
   charge(I) = BAS%Charge(I)
enddo
DoSpherical=.TRUE.
call interface_ao_write_general(maxrep,BAS%MAXNSHELL,BAS%natoms,maxprim,nCC,&
     & BAS%SHELLANGMOM,BAS%nAtoms,BAS%CCINDEX,nPrim,Charge,Coord,priexp,CC,&
     & DoSpherical)

call mem_dealloc(Charge)
call mem_dealloc(Coord)
call mem_dealloc(CC)
call mem_dealloc(priexp)
call mem_dealloc(nPrim)
CALL FREE_BASINF(BAS)

END SUBROUTINE XC_HOST_INTERFACE

END MODULE XCHOST
