module SPAGC_CPU_OBS_HorizontalRecurrenceLHSModBtoA
 use IchorPrecisionMod
  
 CONTAINS

subroutine SPHorizontalRR_CPU_LHS_P1A0B1BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(reals),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(reals),intent(in) :: AuxCont(    4,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(inout) :: ThetaP(    1:    1,    2:    4,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iTUVB) 
  DO iP = 1,nTUVQ*nContQP*nPasses
     DO iTUVB=  2,  4
        ThetaP(1,iTUVB,iP) = AuxCont(iTUVB,iP)
     ENDDO
  ENDDO
!$OMP END DO
end subroutine SPHorizontalRR_CPU_LHS_P1A0B1BtoA

subroutine SPHorizontalRR_CPU_LHS_P2A0B2BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(reals),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(reals),intent(in) :: AuxCont(   10,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(inout) :: ThetaP(    1:    1,    5:   10,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iTUVB) 
  DO iP = 1,nTUVQ*nContQP*nPasses
     DO iTUVB=  5, 10
        ThetaP(1,iTUVB,iP) = AuxCont(iTUVB,iP)
     ENDDO
  ENDDO
!$OMP END DO
end subroutine SPHorizontalRR_CPU_LHS_P2A0B2BtoA

subroutine SPHorizontalRR_CPU_LHS_P3A1B2BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(reals),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(reals),intent(in) :: AuxCont(   20,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(inout) :: ThetaP(    2:    4,    5:   10,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
  real(reals) :: Xab,Yab,Zab
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iPassP,iTUVB,iAtomA,iAtomB,Xab,Yab,Zab) 
  DO iP = 1,nTUVQ*nContQP*nPasses
   iPassP = (iP-1)/(nTUVQ*nContQP)+1
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = -Pdistance12(1,iAtomA,iAtomB)
   Yab = -Pdistance12(2,iAtomA,iAtomB)
   Zab = -Pdistance12(3,iAtomA,iAtomB)
     ThetaP(2,5,IP) = AuxCont(11,iP)+ Xab*AuxCont(5,ip) 
     ThetaP(2,6,IP) = AuxCont(12,iP)+ Xab*AuxCont(6,ip) 
     ThetaP(2,7,IP) = AuxCont(13,iP)+ Xab*AuxCont(7,ip) 
     ThetaP(2,8,IP) = AuxCont(14,iP)+ Xab*AuxCont(8,ip) 
     ThetaP(2,9,IP) = AuxCont(15,iP)+ Xab*AuxCont(9,ip) 
     ThetaP(2,10,IP) = AuxCont(16,iP)+ Xab*AuxCont(10,ip) 
     ThetaP(3,5,IP) = AuxCont(12,iP)+ Yab*AuxCont(5,ip) 
     ThetaP(3,6,IP) = AuxCont(14,iP)+ Yab*AuxCont(6,ip) 
     ThetaP(3,7,IP) = AuxCont(15,iP)+ Yab*AuxCont(7,ip) 
     ThetaP(3,8,IP) = AuxCont(17,iP)+ Yab*AuxCont(8,ip) 
     ThetaP(3,9,IP) = AuxCont(18,iP)+ Yab*AuxCont(9,ip) 
     ThetaP(3,10,IP) = AuxCont(19,iP)+ Yab*AuxCont(10,ip) 
     ThetaP(4,5,IP) = AuxCont(13,iP)+ Zab*AuxCont(5,ip) 
     ThetaP(4,6,IP) = AuxCont(15,iP)+ Zab*AuxCont(6,ip) 
     ThetaP(4,7,IP) = AuxCont(16,iP)+ Zab*AuxCont(7,ip) 
     ThetaP(4,8,IP) = AuxCont(18,iP)+ Zab*AuxCont(8,ip) 
     ThetaP(4,9,IP) = AuxCont(19,iP)+ Zab*AuxCont(9,ip) 
     ThetaP(4,10,IP) = AuxCont(20,iP)+ Zab*AuxCont(10,ip) 
  ENDDO
!$OMP END DO
end subroutine SPHorizontalRR_CPU_LHS_P3A1B2BtoA
end module
