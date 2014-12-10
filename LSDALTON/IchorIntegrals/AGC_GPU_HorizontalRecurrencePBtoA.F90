MODULE AGC_GPU_OBS_HorizontalRecurrenceLHSModBtoA
 use IchorPrecisionMod
  
 CONTAINS

subroutine HorizontalRR_GPU_LHS_P1A0B1BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri,iASync)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(nContQP*nPasses,    4,nTUVQ)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(nContQP*nPasses,    1:    1,    2:    4,nTUVQ)
  integer(kind=acckind) :: iASync
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVQ,&
!$ACC         iTUVB) &
!$ACC PRESENT(nPasses,&
!$ACC         AuxCont,ThetaP) ASYNC(iASync)
  DO iP = 1,nContQP*nPasses
   DO iTUVQ = 1,nTUVQ
     DO iTUVB=  2,  4
        ThetaP(iP,1,iTUVB,iTUVQ) = AuxCont(iP,iTUVB,iTUVQ)
     ENDDO
   ENDDO
  ENDDO
end subroutine HorizontalRR_GPU_LHS_P1A0B1BtoA

subroutine HorizontalRR_GPU_LHS_P2A0B2BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri,iASync)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(nContQP*nPasses,   10,nTUVQ)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(nContQP*nPasses,    1:    1,    5:   10,nTUVQ)
  integer(kind=acckind) :: iASync
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVQ,&
!$ACC         iTUVB) &
!$ACC PRESENT(nPasses,&
!$ACC         AuxCont,ThetaP) ASYNC(iASync)
  DO iP = 1,nContQP*nPasses
   DO iTUVQ = 1,nTUVQ
     DO iTUVB=  5, 10
        ThetaP(iP,1,iTUVB,iTUVQ) = AuxCont(iP,iTUVB,iTUVQ)
     ENDDO
   ENDDO
  ENDDO
end subroutine HorizontalRR_GPU_LHS_P2A0B2BtoA

subroutine HorizontalRR_GPU_LHS_P3A1B2BtoA(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri,iASync)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(nContQP*nPasses,   20,nTUVQ)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(nContQP*nPasses,    2:    4,    5:   10,nTUVQ)
  integer(kind=acckind) :: iASync
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVB,iAtomA,iAtomB
  real(realk) :: Xab,Yab,Zab
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVQ,&
!$ACC         iPassP,iTUVB,iAtomA,iAtomB,Xab,Yab,Zab) &
!$ACC PRESENT(nPasses,&
!$ACC         iAtomApass,iAtomBpass,Pdistance12,AuxCont,ThetaP) ASYNC(iASync)
  DO iP = 1,nContQP*nPasses
   DO iTUVQ = 1,nTUVQ
   iPassP = (iP-1)/(nContQP)+1
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = -Pdistance12(1,iAtomA,iAtomB)
   Yab = -Pdistance12(2,iAtomA,iAtomB)
   Zab = -Pdistance12(3,iAtomA,iAtomB)
     ThetaP(iP,2,5,iTUVQ) = AuxCont(iP,11,iTUVQ)+ Xab*AuxCont(iP,5,iTUVQ) 
     ThetaP(iP,2,6,iTUVQ) = AuxCont(iP,12,iTUVQ)+ Xab*AuxCont(iP,6,iTUVQ) 
     ThetaP(iP,2,7,iTUVQ) = AuxCont(iP,13,iTUVQ)+ Xab*AuxCont(iP,7,iTUVQ) 
     ThetaP(iP,2,8,iTUVQ) = AuxCont(iP,14,iTUVQ)+ Xab*AuxCont(iP,8,iTUVQ) 
     ThetaP(iP,2,9,iTUVQ) = AuxCont(iP,15,iTUVQ)+ Xab*AuxCont(iP,9,iTUVQ) 
     ThetaP(iP,2,10,iTUVQ) = AuxCont(iP,16,iTUVQ)+ Xab*AuxCont(iP,10,iTUVQ) 
     ThetaP(iP,3,5,iTUVQ) = AuxCont(iP,12,iTUVQ)+ Yab*AuxCont(iP,5,iTUVQ) 
     ThetaP(iP,3,6,iTUVQ) = AuxCont(iP,14,iTUVQ)+ Yab*AuxCont(iP,6,iTUVQ) 
     ThetaP(iP,3,7,iTUVQ) = AuxCont(iP,15,iTUVQ)+ Yab*AuxCont(iP,7,iTUVQ) 
     ThetaP(iP,3,8,iTUVQ) = AuxCont(iP,17,iTUVQ)+ Yab*AuxCont(iP,8,iTUVQ) 
     ThetaP(iP,3,9,iTUVQ) = AuxCont(iP,18,iTUVQ)+ Yab*AuxCont(iP,9,iTUVQ) 
     ThetaP(iP,3,10,iTUVQ) = AuxCont(iP,19,iTUVQ)+ Yab*AuxCont(iP,10,iTUVQ) 
     ThetaP(iP,4,5,iTUVQ) = AuxCont(iP,13,iTUVQ)+ Zab*AuxCont(iP,5,iTUVQ) 
     ThetaP(iP,4,6,iTUVQ) = AuxCont(iP,15,iTUVQ)+ Zab*AuxCont(iP,6,iTUVQ) 
     ThetaP(iP,4,7,iTUVQ) = AuxCont(iP,16,iTUVQ)+ Zab*AuxCont(iP,7,iTUVQ) 
     ThetaP(iP,4,8,iTUVQ) = AuxCont(iP,18,iTUVQ)+ Zab*AuxCont(iP,8,iTUVQ) 
     ThetaP(iP,4,9,iTUVQ) = AuxCont(iP,19,iTUVQ)+ Zab*AuxCont(iP,9,iTUVQ) 
     ThetaP(iP,4,10,iTUVQ) = AuxCont(iP,20,iTUVQ)+ Zab*AuxCont(iP,10,iTUVQ) 
   ENDDO
  ENDDO
end subroutine HorizontalRR_GPU_LHS_P3A1B2BtoA
end module
