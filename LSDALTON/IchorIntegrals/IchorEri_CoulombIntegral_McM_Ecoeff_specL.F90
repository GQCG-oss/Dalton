!> @file
!> Contains the general McM driver 

!> \brief General McMurchie-Davidson Integral scheme
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorEriCoulombintegralCPUMcMspecLEcoeffMod
use IchorPrecisionMod
use IchorEriCoulombintegralCPUMcMGeneralWTUVMod

CONTAINS
subroutine ICHOR_Ecoeffn_maxAngP0_maxAngA0_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,3),PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass,1,  1)
!  
  integer :: i,ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
     DO i = 1, nPrimP
        EcoeffN(i,ii,1,1) = ETIJ(i,ii,1)*ETIJ(i,ii,2)*ETIJ(i,ii,3)
     ENDDO
  ENDDO
!!$OMP END DO 
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
     DO i = 1, nPrimP
        iAtomA = IatomAPass(ii)
        iAtomB = IatomBPass(ii)
        EcoeffN(i,ii,1,1) = EcoeffN(i,ii,1,1)*PreExpFac(i,iAtomA,iAtomB)
     ENDDO
  ENDDO
!!$OMP END DO 
!$OMP END MASTER
!$OMP BARRIER 
End subroutine ICHOR_Ecoeffn_maxAngP0_maxAngA0_LHS
  
subroutine ICHOR_Ecoeffn_maxAngP1_maxAngA1_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,0:1,0:1,0:0,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass,  4,  3)
!  
  integer :: i,ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,1) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,1) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,2) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,2) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,3) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,3) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
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
     EcoeffN(i,ii,1,2) = EcoeffN(i,ii,1,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,2) = EcoeffN(i,ii,3,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,3) = EcoeffN(i,ii,1,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,3) = EcoeffN(i,ii,4,3)*PreExpFac(i,iAtomA,iAtomB)
   ENDDO
  ENDDO
!!$OMP END DO
!$OMP END MASTER
!$OMP BARRIER 
End subroutine ICHOR_Ecoeffn_maxAngP1_maxAngA1_LHS
  
  
subroutine ICHOR_Ecoeffn_maxAngP1_maxAngA0_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,0:1,0:0,0:1,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass,  4,  3)
!  
  integer :: i,ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,1) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,1) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,2) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,2) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,3) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,3) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
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
     EcoeffN(i,ii,1,2) = EcoeffN(i,ii,1,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,2) = EcoeffN(i,ii,3,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,3) = EcoeffN(i,ii,1,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,3) = EcoeffN(i,ii,4,3)*PreExpFac(i,iAtomA,iAtomB)
   ENDDO
  ENDDO
!!$OMP END DO
!$OMP END MASTER
!$OMP BARRIER 
End subroutine ICHOR_Ecoeffn_maxAngP1_maxAngA0_LHS
  
  
subroutine ICHOR_Ecoeffn_maxAngP2_maxAngA2_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,0:2,0:2,0:0,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass, 10,  6)
!  
  integer :: i,ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,1) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,1) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,1) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,2) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,2) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,2) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,2) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,3) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,3) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,2,3) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,3) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,4) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,4) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,4) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,5) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,5) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,5) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,5) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,6) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,4,6) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,10,6) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,2,0,3)
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
     EcoeffN(i,ii,1,2) = EcoeffN(i,ii,1,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,2) = EcoeffN(i,ii,3,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,2) = EcoeffN(i,ii,2,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,2) = EcoeffN(i,ii,6,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,3) = EcoeffN(i,ii,1,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,3) = EcoeffN(i,ii,4,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,3) = EcoeffN(i,ii,2,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,3) = EcoeffN(i,ii,7,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,4) = EcoeffN(i,ii,1,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,4) = EcoeffN(i,ii,3,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,4) = EcoeffN(i,ii,8,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,5) = EcoeffN(i,ii,1,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,5) = EcoeffN(i,ii,4,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,5) = EcoeffN(i,ii,3,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,5) = EcoeffN(i,ii,9,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,6) = EcoeffN(i,ii,1,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,6) = EcoeffN(i,ii,4,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,6) = EcoeffN(i,ii,10,6)*PreExpFac(i,iAtomA,iAtomB)
   ENDDO
  ENDDO
!!$OMP END DO
!$OMP END MASTER
!$OMP BARRIER 
End subroutine ICHOR_Ecoeffn_maxAngP2_maxAngA2_LHS
  
  
subroutine ICHOR_Ecoeffn_maxAngP2_maxAngA1_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,0:2,0:1,0:1,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass, 10,  9)
!  
  integer :: i,ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,1) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,1) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,1) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,2) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,2) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,2) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,2) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,3) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,3) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,2,3) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,3) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,4) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,4) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,4) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,4) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,5) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,5) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,5) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,6) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,6) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,6) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,6) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,7) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,7) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,2,7) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,7) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,4,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,10,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,1,1,3)
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
     EcoeffN(i,ii,1,2) = EcoeffN(i,ii,1,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,2) = EcoeffN(i,ii,3,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,2) = EcoeffN(i,ii,2,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,2) = EcoeffN(i,ii,6,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,3) = EcoeffN(i,ii,1,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,3) = EcoeffN(i,ii,4,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,3) = EcoeffN(i,ii,2,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,3) = EcoeffN(i,ii,7,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,4) = EcoeffN(i,ii,1,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,4) = EcoeffN(i,ii,3,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,4) = EcoeffN(i,ii,2,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,4) = EcoeffN(i,ii,6,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,5) = EcoeffN(i,ii,1,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,5) = EcoeffN(i,ii,3,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,5) = EcoeffN(i,ii,8,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,6) = EcoeffN(i,ii,1,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,6) = EcoeffN(i,ii,4,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,6) = EcoeffN(i,ii,3,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,6) = EcoeffN(i,ii,9,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,7) = EcoeffN(i,ii,1,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,7) = EcoeffN(i,ii,4,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,7) = EcoeffN(i,ii,2,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,7) = EcoeffN(i,ii,7,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,8) = EcoeffN(i,ii,1,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,8) = EcoeffN(i,ii,4,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,8) = EcoeffN(i,ii,3,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,8) = EcoeffN(i,ii,9,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,9) = EcoeffN(i,ii,1,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,9) = EcoeffN(i,ii,4,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,9) = EcoeffN(i,ii,10,9)*PreExpFac(i,iAtomA,iAtomB)
   ENDDO
  ENDDO
!!$OMP END DO
!$OMP END MASTER
!$OMP BARRIER 
End subroutine ICHOR_Ecoeffn_maxAngP2_maxAngA1_LHS

  
subroutine ICHOR_Ecoeffn_maxAngP2_maxAngA0_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,0:2,0:0,0:2,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass, 10,  6)
!  
  integer :: i,ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,1) = ETIJ(i,ii,0,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,1) = ETIJ(i,ii,1,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,1) = ETIJ(i,ii,2,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,2) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,2) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,2) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,2) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,3) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,3) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,2,3) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,3) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,4) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,4) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,4) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,5) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,5) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,5) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,5) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,6) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,2,3)
     EcoeffN(i,ii,4,6) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,2,3)
     EcoeffN(i,ii,10,6) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,0,2,3)
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
     EcoeffN(i,ii,1,2) = EcoeffN(i,ii,1,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,2) = EcoeffN(i,ii,3,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,2) = EcoeffN(i,ii,2,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,2) = EcoeffN(i,ii,6,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,3) = EcoeffN(i,ii,1,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,3) = EcoeffN(i,ii,4,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,3) = EcoeffN(i,ii,2,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,3) = EcoeffN(i,ii,7,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,4) = EcoeffN(i,ii,1,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,4) = EcoeffN(i,ii,3,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,4) = EcoeffN(i,ii,8,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,5) = EcoeffN(i,ii,1,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,5) = EcoeffN(i,ii,4,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,5) = EcoeffN(i,ii,3,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,5) = EcoeffN(i,ii,9,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,6) = EcoeffN(i,ii,1,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,6) = EcoeffN(i,ii,4,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,6) = EcoeffN(i,ii,10,6)*PreExpFac(i,iAtomA,iAtomB)
   ENDDO
  ENDDO
!!$OMP END DO
!$OMP END MASTER
!$OMP BARRIER 
End subroutine ICHOR_Ecoeffn_maxAngP2_maxAngA0_LHS
  
  
subroutine ICHOR_Ecoeffn_maxAngP3_maxAngA3_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,0:3,0:3,0:0,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass, 20, 10)
!  
  integer :: i,ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,1) = ETIJ(i,ii,0,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,1) = ETIJ(i,ii,1,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,1) = ETIJ(i,ii,2,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,11,1) = ETIJ(i,ii,3,3,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,2) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,2) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,2) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,2) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,2) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,12,2) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,3) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,3) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,2,3) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,3) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,5,3) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,13,3) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,4) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,4) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,4) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,4) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,4) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,14,4) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,5) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,5) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,5) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,5) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,2,5) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,5) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,6,5) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,15,5) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,6) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,4,6) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,10,6) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,2,6) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,7,6) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,16,6) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,1,7) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,3,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,7) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,3,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,7) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,3,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,17,7) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,3,3,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,8,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,18,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,4,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,10,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,3,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,9,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,19,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,1,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,3,0,3)
     EcoeffN(i,ii,4,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,3,0,3)
     EcoeffN(i,ii,10,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,3,0,3)
     EcoeffN(i,ii,20,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,3,3,0,3)
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
     EcoeffN(i,ii,1,2) = EcoeffN(i,ii,1,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,2) = EcoeffN(i,ii,3,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,2) = EcoeffN(i,ii,2,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,2) = EcoeffN(i,ii,6,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,2) = EcoeffN(i,ii,5,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,2) = EcoeffN(i,ii,12,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,3) = EcoeffN(i,ii,1,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,3) = EcoeffN(i,ii,4,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,3) = EcoeffN(i,ii,2,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,3) = EcoeffN(i,ii,7,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,3) = EcoeffN(i,ii,5,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,3) = EcoeffN(i,ii,13,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,4) = EcoeffN(i,ii,1,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,4) = EcoeffN(i,ii,3,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,4) = EcoeffN(i,ii,8,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,4) = EcoeffN(i,ii,2,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,4) = EcoeffN(i,ii,6,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,4) = EcoeffN(i,ii,14,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,5) = EcoeffN(i,ii,1,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,5) = EcoeffN(i,ii,4,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,5) = EcoeffN(i,ii,3,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,5) = EcoeffN(i,ii,9,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,5) = EcoeffN(i,ii,2,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,5) = EcoeffN(i,ii,7,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,5) = EcoeffN(i,ii,6,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,5) = EcoeffN(i,ii,15,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,6) = EcoeffN(i,ii,1,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,6) = EcoeffN(i,ii,4,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,6) = EcoeffN(i,ii,10,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,6) = EcoeffN(i,ii,2,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,6) = EcoeffN(i,ii,7,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,6) = EcoeffN(i,ii,16,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,7) = EcoeffN(i,ii,1,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,7) = EcoeffN(i,ii,3,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,7) = EcoeffN(i,ii,8,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,17,7) = EcoeffN(i,ii,17,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,8) = EcoeffN(i,ii,1,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,8) = EcoeffN(i,ii,4,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,8) = EcoeffN(i,ii,3,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,8) = EcoeffN(i,ii,9,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,8) = EcoeffN(i,ii,8,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,8) = EcoeffN(i,ii,18,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,9) = EcoeffN(i,ii,1,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,9) = EcoeffN(i,ii,4,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,9) = EcoeffN(i,ii,10,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,9) = EcoeffN(i,ii,3,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,9) = EcoeffN(i,ii,9,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,9) = EcoeffN(i,ii,19,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,10) = EcoeffN(i,ii,1,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,10) = EcoeffN(i,ii,4,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,10) = EcoeffN(i,ii,10,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,20,10) = EcoeffN(i,ii,20,10)*PreExpFac(i,iAtomA,iAtomB)
   ENDDO
  ENDDO
!!$OMP END DO
!$OMP END MASTER
!$OMP BARRIER 
End subroutine ICHOR_Ecoeffn_maxAngP3_maxAngA3_LHS
  
  
subroutine ICHOR_Ecoeffn_maxAngP3_maxAngA2_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,0:3,0:2,0:1,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass, 20, 18)
!  
  integer :: i,ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,1) = ETIJ(i,ii,0,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,1) = ETIJ(i,ii,1,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,1) = ETIJ(i,ii,2,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,11,1) = ETIJ(i,ii,3,2,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,2) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,2) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,2) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,2) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,2) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,12,2) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,3) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,3) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,2,3) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,3) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,5,3) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,13,3) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,4) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,4) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,4) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,4) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,4) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,14,4) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,5) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,5) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,5) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,5) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,2,5) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,5) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,6,5) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,15,5) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,6) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,4,6) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,10,6) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,2,6) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,7,6) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,16,6) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,1,7) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,7) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,7) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,7) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,7) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,12,7) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,8) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,8) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,8) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,8) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,8) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,14,8) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,9) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,9) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,9) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,9) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,2,9) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,9) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,6,9) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,15,9) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,2,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,2,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,2,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,17,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,3,2,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,11) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,11) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,11) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,11) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,8,11) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,18,11) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,4,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,10,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,3,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,2,0,3)
     EcoeffN(i,ii,9,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,2,0,3)
     EcoeffN(i,ii,19,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,2,2,0,3)
     EcoeffN(i,ii,1,13) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,13) = ETIJ(i,ii,0,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,2,13) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,13) = ETIJ(i,ii,1,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,5,13) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,13,13) = ETIJ(i,ii,2,2,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,14) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,14) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,14) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,14) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,2,14) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,14) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,6,14) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,15,14) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,15) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,4,15) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,10,15) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,2,15) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,7,15) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,16,15) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,1,16) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,16) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,2,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,16) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,16) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,2,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,8,16) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,18,16) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,2,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,4,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,10,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,3,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,9,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,19,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,1,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,2,1,3)
     EcoeffN(i,ii,4,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,2,1,3)
     EcoeffN(i,ii,10,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,2,1,3)
     EcoeffN(i,ii,20,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,3,2,1,3)
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
     EcoeffN(i,ii,1,2) = EcoeffN(i,ii,1,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,2) = EcoeffN(i,ii,3,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,2) = EcoeffN(i,ii,2,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,2) = EcoeffN(i,ii,6,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,2) = EcoeffN(i,ii,5,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,2) = EcoeffN(i,ii,12,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,3) = EcoeffN(i,ii,1,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,3) = EcoeffN(i,ii,4,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,3) = EcoeffN(i,ii,2,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,3) = EcoeffN(i,ii,7,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,3) = EcoeffN(i,ii,5,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,3) = EcoeffN(i,ii,13,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,4) = EcoeffN(i,ii,1,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,4) = EcoeffN(i,ii,3,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,4) = EcoeffN(i,ii,8,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,4) = EcoeffN(i,ii,2,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,4) = EcoeffN(i,ii,6,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,4) = EcoeffN(i,ii,14,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,5) = EcoeffN(i,ii,1,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,5) = EcoeffN(i,ii,4,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,5) = EcoeffN(i,ii,3,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,5) = EcoeffN(i,ii,9,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,5) = EcoeffN(i,ii,2,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,5) = EcoeffN(i,ii,7,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,5) = EcoeffN(i,ii,6,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,5) = EcoeffN(i,ii,15,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,6) = EcoeffN(i,ii,1,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,6) = EcoeffN(i,ii,4,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,6) = EcoeffN(i,ii,10,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,6) = EcoeffN(i,ii,2,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,6) = EcoeffN(i,ii,7,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,6) = EcoeffN(i,ii,16,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,7) = EcoeffN(i,ii,1,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,7) = EcoeffN(i,ii,3,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,7) = EcoeffN(i,ii,2,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,7) = EcoeffN(i,ii,6,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,7) = EcoeffN(i,ii,5,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,7) = EcoeffN(i,ii,12,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,8) = EcoeffN(i,ii,1,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,8) = EcoeffN(i,ii,3,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,8) = EcoeffN(i,ii,8,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,8) = EcoeffN(i,ii,2,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,8) = EcoeffN(i,ii,6,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,8) = EcoeffN(i,ii,14,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,9) = EcoeffN(i,ii,1,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,9) = EcoeffN(i,ii,4,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,9) = EcoeffN(i,ii,3,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,9) = EcoeffN(i,ii,9,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,9) = EcoeffN(i,ii,2,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,9) = EcoeffN(i,ii,7,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,9) = EcoeffN(i,ii,6,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,9) = EcoeffN(i,ii,15,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,10) = EcoeffN(i,ii,1,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,10) = EcoeffN(i,ii,3,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,10) = EcoeffN(i,ii,8,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,17,10) = EcoeffN(i,ii,17,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,11) = EcoeffN(i,ii,1,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,11) = EcoeffN(i,ii,4,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,11) = EcoeffN(i,ii,3,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,11) = EcoeffN(i,ii,9,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,11) = EcoeffN(i,ii,8,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,11) = EcoeffN(i,ii,18,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,12) = EcoeffN(i,ii,1,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,12) = EcoeffN(i,ii,4,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,12) = EcoeffN(i,ii,10,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,12) = EcoeffN(i,ii,3,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,12) = EcoeffN(i,ii,9,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,12) = EcoeffN(i,ii,19,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,13) = EcoeffN(i,ii,1,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,13) = EcoeffN(i,ii,4,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,13) = EcoeffN(i,ii,2,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,13) = EcoeffN(i,ii,7,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,13) = EcoeffN(i,ii,5,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,13) = EcoeffN(i,ii,13,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,14) = EcoeffN(i,ii,1,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,14) = EcoeffN(i,ii,4,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,14) = EcoeffN(i,ii,3,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,14) = EcoeffN(i,ii,9,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,14) = EcoeffN(i,ii,2,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,14) = EcoeffN(i,ii,7,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,14) = EcoeffN(i,ii,6,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,14) = EcoeffN(i,ii,15,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,15) = EcoeffN(i,ii,1,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,15) = EcoeffN(i,ii,4,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,15) = EcoeffN(i,ii,10,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,15) = EcoeffN(i,ii,2,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,15) = EcoeffN(i,ii,7,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,15) = EcoeffN(i,ii,16,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,16) = EcoeffN(i,ii,1,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,16) = EcoeffN(i,ii,4,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,16) = EcoeffN(i,ii,3,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,16) = EcoeffN(i,ii,9,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,16) = EcoeffN(i,ii,8,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,16) = EcoeffN(i,ii,18,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,17) = EcoeffN(i,ii,1,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,17) = EcoeffN(i,ii,4,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,17) = EcoeffN(i,ii,10,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,17) = EcoeffN(i,ii,3,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,17) = EcoeffN(i,ii,9,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,17) = EcoeffN(i,ii,19,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,18) = EcoeffN(i,ii,1,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,18) = EcoeffN(i,ii,4,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,18) = EcoeffN(i,ii,10,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,20,18) = EcoeffN(i,ii,20,18)*PreExpFac(i,iAtomA,iAtomB)
   ENDDO
  ENDDO
!!$OMP END DO
!$OMP END MASTER
!$OMP BARRIER 
End subroutine ICHOR_Ecoeffn_maxAngP3_maxAngA2_LHS
  
  
subroutine ICHOR_Ecoeffn_maxAngP3_maxAngA1_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,0:3,0:1,0:2,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass, 20, 18)
!  
  integer :: i,ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,1) = ETIJ(i,ii,0,1,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,1) = ETIJ(i,ii,1,1,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,1) = ETIJ(i,ii,2,1,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,11,1) = ETIJ(i,ii,3,1,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,2) = ETIJ(i,ii,0,0,2,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,2) = ETIJ(i,ii,0,0,2,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,2) = ETIJ(i,ii,1,0,2,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,2) = ETIJ(i,ii,1,0,2,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,2) = ETIJ(i,ii,2,0,2,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,12,2) = ETIJ(i,ii,2,0,2,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,3) = ETIJ(i,ii,0,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,3) = ETIJ(i,ii,0,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,2,3) = ETIJ(i,ii,1,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,3) = ETIJ(i,ii,1,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,5,3) = ETIJ(i,ii,2,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,13,3) = ETIJ(i,ii,2,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,4) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,4) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,4) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,4) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,4) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,12,4) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,5) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,5) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,5) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,5) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,5) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,14,5) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,6) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,6) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,6) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,6) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,2,6) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,7,6) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,6,6) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,15,6) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,7) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,7) = ETIJ(i,ii,0,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,2,7) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,7) = ETIJ(i,ii,1,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,5,7) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,13,7) = ETIJ(i,ii,2,1,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,8) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,8) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,8) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,8) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,2,8) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,8) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,6,8) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,15,8) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,9) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,4,9) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,10,9) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,2,9) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,7,9) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,16,9) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,1,10) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,10) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,10) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,2,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,10) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,10) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,14,10) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,2,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,11) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,11) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,11) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,1,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,17,11) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,3,1,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,2,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,4,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,2,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,3,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,2,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,9,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,2,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,8,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,0,2,2)*ETIJ(i,ii,0,1,0,3)
     EcoeffN(i,ii,18,12) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,0,2,2)*ETIJ(i,ii,1,1,0,3)
     EcoeffN(i,ii,1,13) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,13) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,13) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,13) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,2,13) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,13) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,6,13) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,15,13) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,14) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,14) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,14) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,14) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,8,14) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,18,14) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,1,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,15) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,4,15) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,10,15) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,3,15) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,1,1,3)
     EcoeffN(i,ii,9,15) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,1,1,3)
     EcoeffN(i,ii,19,15) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,2,1,1,3)
     EcoeffN(i,ii,1,16) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,2,3)
     EcoeffN(i,ii,4,16) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,2,3)
     EcoeffN(i,ii,10,16) = ETIJ(i,ii,0,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,0,2,3)
     EcoeffN(i,ii,2,16) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,2,3)
     EcoeffN(i,ii,7,16) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,2,3)
     EcoeffN(i,ii,16,16) = ETIJ(i,ii,1,1,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,0,2,3)
     EcoeffN(i,ii,1,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,0,0,2,3)
     EcoeffN(i,ii,4,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,1,0,2,3)
     EcoeffN(i,ii,10,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,1,0,2)*ETIJ(i,ii,2,0,2,3)
     EcoeffN(i,ii,3,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,0,0,2,3)
     EcoeffN(i,ii,9,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,1,0,2,3)
     EcoeffN(i,ii,19,17) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,1,0,2)*ETIJ(i,ii,2,0,2,3)
     EcoeffN(i,ii,1,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,1,2,3)
     EcoeffN(i,ii,4,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,1,2,3)
     EcoeffN(i,ii,10,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,1,2,3)
     EcoeffN(i,ii,20,18) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,3,1,2,3)
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
     EcoeffN(i,ii,1,2) = EcoeffN(i,ii,1,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,2) = EcoeffN(i,ii,3,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,2) = EcoeffN(i,ii,2,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,2) = EcoeffN(i,ii,6,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,2) = EcoeffN(i,ii,5,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,2) = EcoeffN(i,ii,12,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,3) = EcoeffN(i,ii,1,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,3) = EcoeffN(i,ii,4,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,3) = EcoeffN(i,ii,2,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,3) = EcoeffN(i,ii,7,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,3) = EcoeffN(i,ii,5,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,3) = EcoeffN(i,ii,13,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,4) = EcoeffN(i,ii,1,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,4) = EcoeffN(i,ii,3,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,4) = EcoeffN(i,ii,2,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,4) = EcoeffN(i,ii,6,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,4) = EcoeffN(i,ii,5,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,4) = EcoeffN(i,ii,12,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,5) = EcoeffN(i,ii,1,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,5) = EcoeffN(i,ii,3,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,5) = EcoeffN(i,ii,8,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,5) = EcoeffN(i,ii,2,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,5) = EcoeffN(i,ii,6,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,5) = EcoeffN(i,ii,14,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,6) = EcoeffN(i,ii,1,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,6) = EcoeffN(i,ii,4,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,6) = EcoeffN(i,ii,3,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,6) = EcoeffN(i,ii,9,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,6) = EcoeffN(i,ii,2,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,6) = EcoeffN(i,ii,7,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,6) = EcoeffN(i,ii,6,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,6) = EcoeffN(i,ii,15,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,7) = EcoeffN(i,ii,1,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,7) = EcoeffN(i,ii,4,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,7) = EcoeffN(i,ii,2,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,7) = EcoeffN(i,ii,7,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,7) = EcoeffN(i,ii,5,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,7) = EcoeffN(i,ii,13,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,8) = EcoeffN(i,ii,1,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,8) = EcoeffN(i,ii,4,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,8) = EcoeffN(i,ii,3,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,8) = EcoeffN(i,ii,9,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,8) = EcoeffN(i,ii,2,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,8) = EcoeffN(i,ii,7,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,8) = EcoeffN(i,ii,6,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,8) = EcoeffN(i,ii,15,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,9) = EcoeffN(i,ii,1,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,9) = EcoeffN(i,ii,4,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,9) = EcoeffN(i,ii,10,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,9) = EcoeffN(i,ii,2,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,9) = EcoeffN(i,ii,7,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,9) = EcoeffN(i,ii,16,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,10) = EcoeffN(i,ii,1,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,10) = EcoeffN(i,ii,3,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,10) = EcoeffN(i,ii,8,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,10) = EcoeffN(i,ii,2,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,10) = EcoeffN(i,ii,6,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,10) = EcoeffN(i,ii,14,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,11) = EcoeffN(i,ii,1,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,11) = EcoeffN(i,ii,3,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,11) = EcoeffN(i,ii,8,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,17,11) = EcoeffN(i,ii,17,11)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,12) = EcoeffN(i,ii,1,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,12) = EcoeffN(i,ii,4,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,12) = EcoeffN(i,ii,3,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,12) = EcoeffN(i,ii,9,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,12) = EcoeffN(i,ii,8,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,12) = EcoeffN(i,ii,18,12)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,13) = EcoeffN(i,ii,1,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,13) = EcoeffN(i,ii,4,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,13) = EcoeffN(i,ii,3,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,13) = EcoeffN(i,ii,9,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,13) = EcoeffN(i,ii,2,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,13) = EcoeffN(i,ii,7,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,13) = EcoeffN(i,ii,6,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,13) = EcoeffN(i,ii,15,13)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,14) = EcoeffN(i,ii,1,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,14) = EcoeffN(i,ii,4,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,14) = EcoeffN(i,ii,3,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,14) = EcoeffN(i,ii,9,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,14) = EcoeffN(i,ii,8,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,14) = EcoeffN(i,ii,18,14)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,15) = EcoeffN(i,ii,1,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,15) = EcoeffN(i,ii,4,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,15) = EcoeffN(i,ii,10,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,15) = EcoeffN(i,ii,3,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,15) = EcoeffN(i,ii,9,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,15) = EcoeffN(i,ii,19,15)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,16) = EcoeffN(i,ii,1,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,16) = EcoeffN(i,ii,4,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,16) = EcoeffN(i,ii,10,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,16) = EcoeffN(i,ii,2,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,16) = EcoeffN(i,ii,7,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,16) = EcoeffN(i,ii,16,16)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,17) = EcoeffN(i,ii,1,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,17) = EcoeffN(i,ii,4,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,17) = EcoeffN(i,ii,10,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,17) = EcoeffN(i,ii,3,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,17) = EcoeffN(i,ii,9,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,17) = EcoeffN(i,ii,19,17)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,18) = EcoeffN(i,ii,1,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,18) = EcoeffN(i,ii,4,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,18) = EcoeffN(i,ii,10,18)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,20,18) = EcoeffN(i,ii,20,18)*PreExpFac(i,iAtomA,iAtomB)
   ENDDO
  ENDDO
!!$OMP END DO
!$OMP END MASTER
!$OMP BARRIER 
End subroutine ICHOR_Ecoeffn_maxAngP3_maxAngA1_LHS
  
  
subroutine ICHOR_Ecoeffn_maxAngP3_maxAngA0_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,0:3,0:0,0:3,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass, 20, 10)
!  
  integer :: i,ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO COLLAPSE(2) PRIVATE(I,II)
  DO ii = 1, nPass
   DO i = 1, nPrimP
     EcoeffN(i,ii,1,1) = ETIJ(i,ii,0,0,3,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,1) = ETIJ(i,ii,1,0,3,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,1) = ETIJ(i,ii,2,0,3,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,11,1) = ETIJ(i,ii,3,0,3,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,2) = ETIJ(i,ii,0,0,2,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,2) = ETIJ(i,ii,0,0,2,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,2) = ETIJ(i,ii,1,0,2,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,2) = ETIJ(i,ii,1,0,2,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,5,2) = ETIJ(i,ii,2,0,2,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,12,2) = ETIJ(i,ii,2,0,2,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,3) = ETIJ(i,ii,0,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,3) = ETIJ(i,ii,0,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,2,3) = ETIJ(i,ii,1,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,3) = ETIJ(i,ii,1,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,5,3) = ETIJ(i,ii,2,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,13,3) = ETIJ(i,ii,2,0,2,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,4) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,4) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,4) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,2,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,2,4) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,6,4) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,14,4) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,2,0,2,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,5) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,5) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,5) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,5) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,2,5) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,7,5) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,6,5) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,15,5) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,6) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,2,3)
     EcoeffN(i,ii,4,6) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,2,3)
     EcoeffN(i,ii,10,6) = ETIJ(i,ii,0,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,0,2,3)
     EcoeffN(i,ii,2,6) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,2,3)
     EcoeffN(i,ii,7,6) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,2,3)
     EcoeffN(i,ii,16,6) = ETIJ(i,ii,1,0,1,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,0,2,3)
     EcoeffN(i,ii,1,7) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,3,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,3,7) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,3,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,8,7) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,0,3,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,17,7) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,3,0,3,2)*ETIJ(i,ii,0,0,0,3)
     EcoeffN(i,ii,1,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,2,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,4,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,2,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,3,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,2,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,9,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,2,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,8,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,0,2,2)*ETIJ(i,ii,0,0,1,3)
     EcoeffN(i,ii,18,8) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,2,0,2,2)*ETIJ(i,ii,1,0,1,3)
     EcoeffN(i,ii,1,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,0,0,2,3)
     EcoeffN(i,ii,4,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,1,0,2,3)
     EcoeffN(i,ii,10,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,1,2)*ETIJ(i,ii,2,0,2,3)
     EcoeffN(i,ii,3,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,0,0,2,3)
     EcoeffN(i,ii,9,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,1,0,2,3)
     EcoeffN(i,ii,19,9) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,1,0,1,2)*ETIJ(i,ii,2,0,2,3)
     EcoeffN(i,ii,1,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,0,0,3,3)
     EcoeffN(i,ii,4,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,1,0,3,3)
     EcoeffN(i,ii,10,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,2,0,3,3)
     EcoeffN(i,ii,20,10) = ETIJ(i,ii,0,0,0,1)*ETIJ(i,ii,0,0,0,2)*ETIJ(i,ii,3,0,3,3)
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
     EcoeffN(i,ii,1,2) = EcoeffN(i,ii,1,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,2) = EcoeffN(i,ii,3,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,2) = EcoeffN(i,ii,2,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,2) = EcoeffN(i,ii,6,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,2) = EcoeffN(i,ii,5,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,12,2) = EcoeffN(i,ii,12,2)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,3) = EcoeffN(i,ii,1,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,3) = EcoeffN(i,ii,4,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,3) = EcoeffN(i,ii,2,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,3) = EcoeffN(i,ii,7,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,5,3) = EcoeffN(i,ii,5,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,13,3) = EcoeffN(i,ii,13,3)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,4) = EcoeffN(i,ii,1,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,4) = EcoeffN(i,ii,3,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,4) = EcoeffN(i,ii,8,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,4) = EcoeffN(i,ii,2,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,4) = EcoeffN(i,ii,6,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,14,4) = EcoeffN(i,ii,14,4)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,5) = EcoeffN(i,ii,1,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,5) = EcoeffN(i,ii,4,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,5) = EcoeffN(i,ii,3,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,5) = EcoeffN(i,ii,9,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,5) = EcoeffN(i,ii,2,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,5) = EcoeffN(i,ii,7,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,6,5) = EcoeffN(i,ii,6,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,15,5) = EcoeffN(i,ii,15,5)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,6) = EcoeffN(i,ii,1,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,6) = EcoeffN(i,ii,4,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,6) = EcoeffN(i,ii,10,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,2,6) = EcoeffN(i,ii,2,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,7,6) = EcoeffN(i,ii,7,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,16,6) = EcoeffN(i,ii,16,6)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,7) = EcoeffN(i,ii,1,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,7) = EcoeffN(i,ii,3,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,7) = EcoeffN(i,ii,8,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,17,7) = EcoeffN(i,ii,17,7)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,8) = EcoeffN(i,ii,1,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,8) = EcoeffN(i,ii,4,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,8) = EcoeffN(i,ii,3,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,8) = EcoeffN(i,ii,9,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,8,8) = EcoeffN(i,ii,8,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,18,8) = EcoeffN(i,ii,18,8)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,9) = EcoeffN(i,ii,1,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,9) = EcoeffN(i,ii,4,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,9) = EcoeffN(i,ii,10,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,3,9) = EcoeffN(i,ii,3,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,9,9) = EcoeffN(i,ii,9,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,19,9) = EcoeffN(i,ii,19,9)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,1,10) = EcoeffN(i,ii,1,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,4,10) = EcoeffN(i,ii,4,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,10,10) = EcoeffN(i,ii,10,10)*PreExpFac(i,iAtomA,iAtomB)
     EcoeffN(i,ii,20,10) = EcoeffN(i,ii,20,10)*PreExpFac(i,iAtomA,iAtomB)
   ENDDO
  ENDDO
!!$OMP END DO
!$OMP END MASTER
!$OMP BARRIER 
End subroutine ICHOR_Ecoeffn_maxAngP3_maxAngA0_LHS
  
subroutine ICHOR_Ecoeffn_general_LHS(nPrimP,nPass,nTUV,ijk,JMAX,l1,l2,EcoeffN,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,nTUV,ijk,JMAX,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nPrimP,nPass,0:JMAX,0:l1,0:l2,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: EcoeffN(nPrimP,nPass,nTUV,ijk)
  !
  integer :: l2,l1,iP2,JP2,kp2,jp1,ip1,kp1,tp,up,vp,ituvp,J,ijkQ,i
  integer :: ii,iAtomA,iAtomB
!$OMP MASTER
!!$OMP DO PRIVATE(ii,i,iAtomA,iAtomB,ijkQ,iP2,jP2,kP2,iP1,&
!!$OMP            jP1,kP1,tp,up,vp,ituvP)
  DO ii = 1, nPass
   DO i = 1, nPrimP
    iAtomA = IatomAPass(ii)
    iAtomB = IatomBPass(ii)
    ijkQ=0
    DO iP2=l2,0,-1
     DO jP2=l2-iP2,0,-1
      kP2=l2-iP2-jP2
      DO iP1=l1,0,-1
       DO jP1=l1-iP1,0,-1
        kP1=l1-iP1-jP1
        ijkQ=ijkQ+1
        DO tP=0,iP1+iP2
         DO uP=0,jP1+jP2
          DO vP=0,kP1+kP2
             ituvP = IchorTUVindexFuncFull(tp,up,vp)
             EcoeffN(i,ii,ituvP,ijkQ) = ETIJ(i,ii,tp,iP1,iP2,1) &
                  & *ETIJ(i,ii,up,jP1,jP2,2)*ETIJ(i,ii,vp,kP1,kP2,3)*PreExpFac(i,iAtomA,iAtomB)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!!$OMP END DO
!$OMP END MASTER
!$OMP BARRIER 
end subroutine ICHOR_Ecoeffn_general_LHS

end MODULE IchorEriCoulombintegralCPUMcMspecLEcoeffMod
