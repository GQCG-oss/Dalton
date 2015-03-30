!> @file
!> Contains the general McM driver 

!> \brief General McMurchie-Davidson Integral scheme
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorEriCoulombintegralCPUMcMspecL2EcoeffMod
use IchorPrecisionMod
use IchorEriCoulombintegralCPUMcMGeneralWTUVMod

CONTAINS
subroutine ICHOR_Ecoeffn_maxAngP4_maxAngA3_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,0:4,0:3,0:1,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass, 35, 30)
!  
  integer :: i,ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,1) = ETIJ(i,ii,0,3,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,1) = ETIJ(i,ii,1,3,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,1) = ETIJ(i,ii,2,3,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,11,1) = ETIJ(i,ii,3,3,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,21,1) = ETIJ(i,ii,4,3,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,2) = ETIJ(i,ii,0,2,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,2) = ETIJ(i,ii,0,2,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,2) = ETIJ(i,ii,1,2,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,2) = ETIJ(i,ii,1,2,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,2) = ETIJ(i,ii,2,2,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,12,2) = ETIJ(i,ii,2,2,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,11,2) = ETIJ(i,ii,3,2,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,22,2) = ETIJ(i,ii,3,2,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,3) = ETIJ(i,ii,0,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,3) = ETIJ(i,ii,0,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,2,3) = ETIJ(i,ii,1,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,3) = ETIJ(i,ii,1,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,5,3) = ETIJ(i,ii,2,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,13,3) = ETIJ(i,ii,2,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,11,3) = ETIJ(i,ii,3,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,23,3) = ETIJ(i,ii,3,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,4) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,4) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,4) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,4) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,4) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,14,4) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,5,4) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,12,4) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,24,4) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,5) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,5) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,5) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,5) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,5) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,5) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,6,5) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,15,5) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,5,5) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,13,5) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,12,5) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,25,5) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,6) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,4,6) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,10,6) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,2,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,6) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,7,6) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,16,6) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,2,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,5,6) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,13,6) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,26,6) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,2,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,7) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,3,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,7) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,3,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,7) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,2,3,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,17,7) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,3,3,0,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,7) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,3,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,7) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,3,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,14,7) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,2,3,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,27,7) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,3,3,0,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,8) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,8) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,8) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,8) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,8,8) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,18,8) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,1,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,8) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,8) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,6,8) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,15,8) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,14,8) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,28,8) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,1,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,9) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,4,9) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,10,9) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,3,9) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,9,9) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,19,9) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,2,2,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,9) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,7,9) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,16,9) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,6,9) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,15,9) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,29,9) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,2,2,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,10) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,3,0,3)
     EcoeffN(i,ii,4,10) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,3,0,3)
     EcoeffN(i,ii,10,10) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,3,0,3)
     EcoeffN(i,ii,20,10) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,3,3,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,10) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,3,0,3)
     EcoeffN(i,ii,7,10) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,3,0,3)
     EcoeffN(i,ii,16,10) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,3,0,3)
     EcoeffN(i,ii,30,10) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,3,3,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,11) = ETIJ(i,ii,0,3,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,11) = ETIJ(i,ii,0,3,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,11) = ETIJ(i,ii,1,3,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,11) = ETIJ(i,ii,1,3,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,11) = ETIJ(i,ii,2,3,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,12,11) = ETIJ(i,ii,2,3,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,11,11) = ETIJ(i,ii,3,3,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,22,11) = ETIJ(i,ii,3,3,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,12) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,12) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,12) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,12) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,12) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,14,12) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,12) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,12,12) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,24,12) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,13) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,13) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,13) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,13) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,13) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,13) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,6,13) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,15,13) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,5,13) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,13,13) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,12,13) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,25,13) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,14) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,2,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,14) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,2,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,14) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,2,2,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,17,14) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,3,2,1,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,14) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,2,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,14) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,2,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,14,14) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,2,2,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,27,14) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,3,2,1,2)*ETIJ(i,ii,0,0,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,15) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,15) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,15) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,15) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,8,15) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,18,15) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,1,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,15) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,15) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,6,15) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,15,15) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,14,15) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,28,15) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,1,1,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,16) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,4,16) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,10,16) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,3,16) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,9,16) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,19,16) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,2,2,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,16) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,7,16) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,16,16) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,6,16) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,15,16) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,29,16) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,2,2,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,3,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,3,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,3,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,17,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,3,3,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,31,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,4,3,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,2,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,2,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,2,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,2,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,8,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,2,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,18,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,2,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,17,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,3,2,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,32,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,3,2,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,19) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,4,19) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,10,19) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,3,19) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,9,19) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,19,19) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,8,19) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,18,19) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,33,19) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,1,20) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,3,0,3)
     EcoeffN(i,ii,4,20) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,3,0,3)
     EcoeffN(i,ii,10,20) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,2,3,0,3)
     EcoeffN(i,ii,20,20) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,3,3,0,3)
     EcoeffN(i,ii,3,20) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,3,0,3)
     EcoeffN(i,ii,9,20) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,3,0,3)
     EcoeffN(i,ii,19,20) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,2,3,0,3)
     EcoeffN(i,ii,34,20) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,3,3,0,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,21) = ETIJ(i,ii,0,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,21) = ETIJ(i,ii,0,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,2,21) = ETIJ(i,ii,1,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,21) = ETIJ(i,ii,1,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,5,21) = ETIJ(i,ii,2,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,13,21) = ETIJ(i,ii,2,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,11,21) = ETIJ(i,ii,3,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,23,21) = ETIJ(i,ii,3,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,22) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,22) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,22) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,22) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,0,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,22) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,22) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,6,22) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,15,22) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,0,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,5,22) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,13,22) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,12,22) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,25,22) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,0,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,23) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,4,23) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,10,23) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,1,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,23) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,7,23) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,16,23) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,1,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,5,23) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,13,23) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,26,23) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,1,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,24) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,24) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,24) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,24) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,8,24) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,18,24) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,1,0,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,24) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,24) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,6,24) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,15,24) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,14,24) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,28,24) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,1,0,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,25) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,4,25) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,10,25) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,3,25) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,9,25) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,19,25) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,2,1,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,25) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,7,25) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,16,25) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,6,25) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,15,25) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,29,25) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,2,1,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,26) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,2,1,3)
     EcoeffN(i,ii,4,26) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,2,1,3)
     EcoeffN(i,ii,10,26) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,2,1,3)
     EcoeffN(i,ii,20,26) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,3,2,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,2,26) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,2,1,3)
     EcoeffN(i,ii,7,26) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,2,1,3)
     EcoeffN(i,ii,16,26) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,2,1,3)
     EcoeffN(i,ii,30,26) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,3,2,1,3)
   ENDDO
  ENDDO
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,27) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,3,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,27) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,3,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,27) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,3,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,27) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,3,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,8,27) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,3,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,18,27) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,3,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,17,27) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,3,3,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,32,27) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,3,3,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,28) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,4,28) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,10,28) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,3,28) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,9,28) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,19,28) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,8,28) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,18,28) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,33,28) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,1,29) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,2,1,3)
     EcoeffN(i,ii,4,29) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,2,1,3)
     EcoeffN(i,ii,10,29) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,2,2,1,3)
     EcoeffN(i,ii,20,29) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,3,2,1,3)
     EcoeffN(i,ii,3,29) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,2,1,3)
     EcoeffN(i,ii,9,29) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,2,1,3)
     EcoeffN(i,ii,19,29) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,2,2,1,3)
     EcoeffN(i,ii,34,29) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,3,2,1,3)
     EcoeffN(i,ii,1,30) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,3,1,3)
     EcoeffN(i,ii,4,30) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,3,1,3)
     EcoeffN(i,ii,10,30) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,3,1,3)
     EcoeffN(i,ii,20,30) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,3,3,1,3)
     EcoeffN(i,ii,35,30) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,4,3,1,3)
   ENDDO
  ENDDO
!!$OMP END DO
!!$OMP DO COLLAPSE(2) PRIVATE(I,II,iAtomA,iAtomB)
  DO ii = 1, nPass
   DO i = 1, nPrimP
    iAtomA = IatomAPass(ii)
    iAtomB = IatomBPass(ii)
     EcoeffN(i,ii,1,1) = EcoeffN(i,ii,1,1)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,1) = EcoeffN(i,ii,2,1)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,1) = EcoeffN(i,ii,5,1)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,11,1) = EcoeffN(i,ii,11,1)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,21,1) = EcoeffN(i,ii,21,1)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,2) = EcoeffN(i,ii,1,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,2) = EcoeffN(i,ii,3,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,2) = EcoeffN(i,ii,2,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,2) = EcoeffN(i,ii,6,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,2) = EcoeffN(i,ii,5,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,2) = EcoeffN(i,ii,12,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,11,2) = EcoeffN(i,ii,11,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,22,2) = EcoeffN(i,ii,22,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,3) = EcoeffN(i,ii,1,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,3) = EcoeffN(i,ii,4,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,3) = EcoeffN(i,ii,2,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,3) = EcoeffN(i,ii,7,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,3) = EcoeffN(i,ii,5,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,3) = EcoeffN(i,ii,13,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,11,3) = EcoeffN(i,ii,11,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,23,3) = EcoeffN(i,ii,23,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,4) = EcoeffN(i,ii,1,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,4) = EcoeffN(i,ii,3,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,4) = EcoeffN(i,ii,8,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,4) = EcoeffN(i,ii,2,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,4) = EcoeffN(i,ii,6,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,4) = EcoeffN(i,ii,14,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,4) = EcoeffN(i,ii,5,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,4) = EcoeffN(i,ii,12,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,24,4) = EcoeffN(i,ii,24,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,5) = EcoeffN(i,ii,1,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,5) = EcoeffN(i,ii,4,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,5) = EcoeffN(i,ii,3,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,5) = EcoeffN(i,ii,9,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,5) = EcoeffN(i,ii,2,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,5) = EcoeffN(i,ii,7,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,5) = EcoeffN(i,ii,6,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,5) = EcoeffN(i,ii,15,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,5) = EcoeffN(i,ii,5,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,5) = EcoeffN(i,ii,13,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,5) = EcoeffN(i,ii,12,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,25,5) = EcoeffN(i,ii,25,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,6) = EcoeffN(i,ii,1,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,6) = EcoeffN(i,ii,4,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,6) = EcoeffN(i,ii,10,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,6) = EcoeffN(i,ii,2,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,6) = EcoeffN(i,ii,7,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,6) = EcoeffN(i,ii,16,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,6) = EcoeffN(i,ii,5,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,6) = EcoeffN(i,ii,13,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,26,6) = EcoeffN(i,ii,26,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,7) = EcoeffN(i,ii,1,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,7) = EcoeffN(i,ii,3,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,7) = EcoeffN(i,ii,8,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,17,7) = EcoeffN(i,ii,17,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,7) = EcoeffN(i,ii,2,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,7) = EcoeffN(i,ii,6,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,7) = EcoeffN(i,ii,14,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,27,7) = EcoeffN(i,ii,27,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,8) = EcoeffN(i,ii,1,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,8) = EcoeffN(i,ii,4,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,8) = EcoeffN(i,ii,3,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,8) = EcoeffN(i,ii,9,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,8) = EcoeffN(i,ii,8,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,8) = EcoeffN(i,ii,18,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,8) = EcoeffN(i,ii,2,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,8) = EcoeffN(i,ii,7,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,8) = EcoeffN(i,ii,6,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,8) = EcoeffN(i,ii,15,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,8) = EcoeffN(i,ii,14,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,28,8) = EcoeffN(i,ii,28,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,9) = EcoeffN(i,ii,1,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,9) = EcoeffN(i,ii,4,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,9) = EcoeffN(i,ii,10,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,9) = EcoeffN(i,ii,3,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,9) = EcoeffN(i,ii,9,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,9) = EcoeffN(i,ii,19,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,9) = EcoeffN(i,ii,2,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,9) = EcoeffN(i,ii,7,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,9) = EcoeffN(i,ii,16,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,9) = EcoeffN(i,ii,6,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,9) = EcoeffN(i,ii,15,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,29,9) = EcoeffN(i,ii,29,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,10) = EcoeffN(i,ii,1,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,10) = EcoeffN(i,ii,4,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,10) = EcoeffN(i,ii,10,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,20,10) = EcoeffN(i,ii,20,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,10) = EcoeffN(i,ii,2,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,10) = EcoeffN(i,ii,7,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,10) = EcoeffN(i,ii,16,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,30,10) = EcoeffN(i,ii,30,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,11) = EcoeffN(i,ii,1,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,11) = EcoeffN(i,ii,3,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,11) = EcoeffN(i,ii,2,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,11) = EcoeffN(i,ii,6,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,11) = EcoeffN(i,ii,5,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,11) = EcoeffN(i,ii,12,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,11,11) = EcoeffN(i,ii,11,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,22,11) = EcoeffN(i,ii,22,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,12) = EcoeffN(i,ii,1,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,12) = EcoeffN(i,ii,3,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,12) = EcoeffN(i,ii,8,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,12) = EcoeffN(i,ii,2,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,12) = EcoeffN(i,ii,6,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,12) = EcoeffN(i,ii,14,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,12) = EcoeffN(i,ii,5,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,12) = EcoeffN(i,ii,12,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,24,12) = EcoeffN(i,ii,24,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,13) = EcoeffN(i,ii,1,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,13) = EcoeffN(i,ii,4,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,13) = EcoeffN(i,ii,3,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,13) = EcoeffN(i,ii,9,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,13) = EcoeffN(i,ii,2,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,13) = EcoeffN(i,ii,7,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,13) = EcoeffN(i,ii,6,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,13) = EcoeffN(i,ii,15,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,13) = EcoeffN(i,ii,5,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,13) = EcoeffN(i,ii,13,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,13) = EcoeffN(i,ii,12,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,25,13) = EcoeffN(i,ii,25,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,14) = EcoeffN(i,ii,1,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,14) = EcoeffN(i,ii,3,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,14) = EcoeffN(i,ii,8,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,17,14) = EcoeffN(i,ii,17,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,14) = EcoeffN(i,ii,2,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,14) = EcoeffN(i,ii,6,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,14) = EcoeffN(i,ii,14,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,27,14) = EcoeffN(i,ii,27,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,15) = EcoeffN(i,ii,1,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,15) = EcoeffN(i,ii,4,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,15) = EcoeffN(i,ii,3,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,15) = EcoeffN(i,ii,9,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,15) = EcoeffN(i,ii,8,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,15) = EcoeffN(i,ii,18,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,15) = EcoeffN(i,ii,2,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,15) = EcoeffN(i,ii,7,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,15) = EcoeffN(i,ii,6,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,15) = EcoeffN(i,ii,15,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,15) = EcoeffN(i,ii,14,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,28,15) = EcoeffN(i,ii,28,15)*PreExpFac(i,iAtomA,iAtomB)
    ENDDO
   ENDDO
   DO ii = 1, nPass
    DO i = 1, nPrimP
     iAtomA = IatomAPass(ii)
     iAtomB = IatomBPass(ii)
     EcoeffN(i,ii,1,16) = EcoeffN(i,ii,1,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,16) = EcoeffN(i,ii,4,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,16) = EcoeffN(i,ii,10,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,16) = EcoeffN(i,ii,3,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,16) = EcoeffN(i,ii,9,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,16) = EcoeffN(i,ii,19,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,16) = EcoeffN(i,ii,2,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,16) = EcoeffN(i,ii,7,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,16) = EcoeffN(i,ii,16,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,16) = EcoeffN(i,ii,6,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,16) = EcoeffN(i,ii,15,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,29,16) = EcoeffN(i,ii,29,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,17) = EcoeffN(i,ii,1,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,17) = EcoeffN(i,ii,3,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,17) = EcoeffN(i,ii,8,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,17,17) = EcoeffN(i,ii,17,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,31,17) = EcoeffN(i,ii,31,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,18) = EcoeffN(i,ii,1,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,18) = EcoeffN(i,ii,4,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,18) = EcoeffN(i,ii,3,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,18) = EcoeffN(i,ii,9,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,18) = EcoeffN(i,ii,8,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,18) = EcoeffN(i,ii,18,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,17,18) = EcoeffN(i,ii,17,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,32,18) = EcoeffN(i,ii,32,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,19) = EcoeffN(i,ii,1,19)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,19) = EcoeffN(i,ii,4,19)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,19) = EcoeffN(i,ii,10,19)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,19) = EcoeffN(i,ii,3,19)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,19) = EcoeffN(i,ii,9,19)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,19) = EcoeffN(i,ii,19,19)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,19) = EcoeffN(i,ii,8,19)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,19) = EcoeffN(i,ii,18,19)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,33,19) = EcoeffN(i,ii,33,19)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,20) = EcoeffN(i,ii,1,20)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,20) = EcoeffN(i,ii,4,20)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,20) = EcoeffN(i,ii,10,20)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,20,20) = EcoeffN(i,ii,20,20)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,20) = EcoeffN(i,ii,3,20)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,20) = EcoeffN(i,ii,9,20)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,20) = EcoeffN(i,ii,19,20)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,34,20) = EcoeffN(i,ii,34,20)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,21) = EcoeffN(i,ii,1,21)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,21) = EcoeffN(i,ii,4,21)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,21) = EcoeffN(i,ii,2,21)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,21) = EcoeffN(i,ii,7,21)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,21) = EcoeffN(i,ii,5,21)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,21) = EcoeffN(i,ii,13,21)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,11,21) = EcoeffN(i,ii,11,21)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,23,21) = EcoeffN(i,ii,23,21)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,22) = EcoeffN(i,ii,1,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,22) = EcoeffN(i,ii,4,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,22) = EcoeffN(i,ii,3,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,22) = EcoeffN(i,ii,9,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,22) = EcoeffN(i,ii,2,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,22) = EcoeffN(i,ii,7,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,22) = EcoeffN(i,ii,6,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,22) = EcoeffN(i,ii,15,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,22) = EcoeffN(i,ii,5,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,22) = EcoeffN(i,ii,13,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,22) = EcoeffN(i,ii,12,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,25,22) = EcoeffN(i,ii,25,22)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,23) = EcoeffN(i,ii,1,23)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,23) = EcoeffN(i,ii,4,23)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,23) = EcoeffN(i,ii,10,23)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,23) = EcoeffN(i,ii,2,23)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,23) = EcoeffN(i,ii,7,23)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,23) = EcoeffN(i,ii,16,23)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,23) = EcoeffN(i,ii,5,23)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,23) = EcoeffN(i,ii,13,23)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,26,23) = EcoeffN(i,ii,26,23)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,24) = EcoeffN(i,ii,1,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,24) = EcoeffN(i,ii,4,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,24) = EcoeffN(i,ii,3,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,24) = EcoeffN(i,ii,9,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,24) = EcoeffN(i,ii,8,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,24) = EcoeffN(i,ii,18,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,24) = EcoeffN(i,ii,2,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,24) = EcoeffN(i,ii,7,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,24) = EcoeffN(i,ii,6,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,24) = EcoeffN(i,ii,15,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,24) = EcoeffN(i,ii,14,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,28,24) = EcoeffN(i,ii,28,24)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,25) = EcoeffN(i,ii,1,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,25) = EcoeffN(i,ii,4,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,25) = EcoeffN(i,ii,10,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,25) = EcoeffN(i,ii,3,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,25) = EcoeffN(i,ii,9,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,25) = EcoeffN(i,ii,19,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,25) = EcoeffN(i,ii,2,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,25) = EcoeffN(i,ii,7,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,25) = EcoeffN(i,ii,16,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,25) = EcoeffN(i,ii,6,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,25) = EcoeffN(i,ii,15,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,29,25) = EcoeffN(i,ii,29,25)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,26) = EcoeffN(i,ii,1,26)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,26) = EcoeffN(i,ii,4,26)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,26) = EcoeffN(i,ii,10,26)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,20,26) = EcoeffN(i,ii,20,26)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,26) = EcoeffN(i,ii,2,26)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,26) = EcoeffN(i,ii,7,26)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,26) = EcoeffN(i,ii,16,26)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,30,26) = EcoeffN(i,ii,30,26)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,27) = EcoeffN(i,ii,1,27)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,27) = EcoeffN(i,ii,4,27)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,27) = EcoeffN(i,ii,3,27)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,27) = EcoeffN(i,ii,9,27)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,27) = EcoeffN(i,ii,8,27)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,27) = EcoeffN(i,ii,18,27)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,17,27) = EcoeffN(i,ii,17,27)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,32,27) = EcoeffN(i,ii,32,27)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,28) = EcoeffN(i,ii,1,28)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,28) = EcoeffN(i,ii,4,28)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,28) = EcoeffN(i,ii,10,28)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,28) = EcoeffN(i,ii,3,28)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,28) = EcoeffN(i,ii,9,28)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,28) = EcoeffN(i,ii,19,28)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,28) = EcoeffN(i,ii,8,28)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,28) = EcoeffN(i,ii,18,28)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,33,28) = EcoeffN(i,ii,33,28)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,29) = EcoeffN(i,ii,1,29)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,29) = EcoeffN(i,ii,4,29)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,29) = EcoeffN(i,ii,10,29)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,20,29) = EcoeffN(i,ii,20,29)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,29) = EcoeffN(i,ii,3,29)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,29) = EcoeffN(i,ii,9,29)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,29) = EcoeffN(i,ii,19,29)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,34,29) = EcoeffN(i,ii,34,29)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,30) = EcoeffN(i,ii,1,30)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,30) = EcoeffN(i,ii,4,30)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,30) = EcoeffN(i,ii,10,30)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,20,30) = EcoeffN(i,ii,20,30)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,35,30) = EcoeffN(i,ii,35,30)*PreExpFac(i,iAtomA,iAtomB)
   ENDDO
  ENDDO
!!$OMP END DO
!$OMP END MASTER
!$OMP BARRIER 
End subroutine ICHOR_Ecoeffn_maxAngP4_maxAngA3_LHS  
 
end MODULE IchorEriCoulombintegralCPUMcMspecL2EcoeffMod
