MODULE AGC_GPU_OBS_HorizontalRecurrenceLHSModAtoB
 use IchorPrecisionModule
  
 CONTAINS

subroutine HorizontalRR_GPU_LHS_P1A1B0AtoB(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(nContQP*nPasses,    4,nTUVQ)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(nContQP*nPasses,    2:    4,1,nTUVQ)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVQ,&
!$ACC         iTUVA) &
!$ACC PRESENT(nPasses,&
!$ACC         AuxCont,ThetaP)
  DO iP = 1,nContQP*nPasses
   DO iTUVQ = 1,nTUVQ
     DO iTUVA=  2,  4
        ThetaP(IP,iTUVA,1,iTUVQ) = AuxCont(IP,iTUVA,iTUVQ)
     ENDDO
  ENDDO
   ENDDO
end subroutine HorizontalRR_GPU_LHS_P1A1B0AtoB

subroutine HorizontalRR_GPU_LHS_P2A1B1AtoB(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(nContQP*nPasses,   10,nTUVQ)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(nContQP*nPasses,    2:    4,    2:    4,nTUVQ)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB
  real(realk) :: Xab,Yab,Zab
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVQ,&
!$ACC         iPassP,iTUVA,iAtomA,iAtomB,Xab,Yab,Zab) &
!$ACC PRESENT(nPasses,&
!$ACC         iAtomApass,iAtomBpass,Pdistance12,AuxCont,ThetaP)
  DO iP = 1,nContQP*nPasses
   DO iTUVQ = 1,nTUVQ
    iPassP = (IP-1)/(nContQP) + 1
    iAtomA = iAtomApass(iPassP)
    iAtomB = iAtomBpass(iPassP)
    Xab = Pdistance12(1,iAtomA,iAtomB)
    Yab = Pdistance12(2,iAtomA,iAtomB)
    Zab = Pdistance12(3,iAtomA,iAtomB)
     ThetaP(iP,2,2,iTUVQ) = AuxCont(iP,5,iTUVQ) + Xab*AuxCont(iP,2,iTUVQ) 
     ThetaP(iP,3,2,iTUVQ) = AuxCont(iP,6,iTUVQ) + Xab*AuxCont(iP,3,iTUVQ) 
     ThetaP(iP,4,2,iTUVQ) = AuxCont(iP,7,iTUVQ) + Xab*AuxCont(iP,4,iTUVQ) 
     ThetaP(iP,2,3,iTUVQ) = AuxCont(iP,6,iTUVQ) + Yab*AuxCont(iP,2,iTUVQ) 
     ThetaP(iP,3,3,iTUVQ) = AuxCont(iP,8,iTUVQ) + Yab*AuxCont(iP,3,iTUVQ) 
     ThetaP(iP,4,3,iTUVQ) = AuxCont(iP,9,iTUVQ) + Yab*AuxCont(iP,4,iTUVQ) 
     ThetaP(iP,2,4,iTUVQ) = AuxCont(iP,7,iTUVQ) + Zab*AuxCont(iP,2,iTUVQ) 
     ThetaP(iP,3,4,iTUVQ) = AuxCont(iP,9,iTUVQ) + Zab*AuxCont(iP,3,iTUVQ) 
     ThetaP(iP,4,4,iTUVQ) = AuxCont(iP,10,iTUVQ) + Zab*AuxCont(iP,4,iTUVQ) 
  ENDDO
   ENDDO
end subroutine HorizontalRR_GPU_LHS_P2A1B1AtoB

subroutine HorizontalRR_GPU_LHS_P2A2B0AtoB(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(nContQP*nPasses,   10,nTUVQ)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(nContQP*nPasses,    5:   10,1,nTUVQ)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVQ,&
!$ACC         iTUVA) &
!$ACC PRESENT(nPasses,&
!$ACC         AuxCont,ThetaP)
  DO iP = 1,nContQP*nPasses
   DO iTUVQ = 1,nTUVQ
     DO iTUVA=  5, 10
        ThetaP(IP,iTUVA,1,iTUVQ) = AuxCont(IP,iTUVA,iTUVQ)
     ENDDO
  ENDDO
   ENDDO
end subroutine HorizontalRR_GPU_LHS_P2A2B0AtoB

subroutine HorizontalRR_GPU_LHS_P3A2B1AtoB(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(nContQP*nPasses,   20,nTUVQ)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(nContQP*nPasses,    5:   10,    2:    4,nTUVQ)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB
  real(realk) :: Xab,Yab,Zab
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVQ,&
!$ACC         iPassP,iTUVA,iAtomA,iAtomB,Xab,Yab,Zab) &
!$ACC PRESENT(nPasses,&
!$ACC         iAtomApass,iAtomBpass,Pdistance12,AuxCont,ThetaP)
  DO iP = 1,nContQP*nPasses
   DO iTUVQ = 1,nTUVQ
    iPassP = (IP-1)/(nContQP) + 1
    iAtomA = iAtomApass(iPassP)
    iAtomB = iAtomBpass(iPassP)
    Xab = Pdistance12(1,iAtomA,iAtomB)
    Yab = Pdistance12(2,iAtomA,iAtomB)
    Zab = Pdistance12(3,iAtomA,iAtomB)
     ThetaP(iP,5,2,iTUVQ) = AuxCont(iP,11,iTUVQ) + Xab*AuxCont(iP,5,iTUVQ) 
     ThetaP(iP,6,2,iTUVQ) = AuxCont(iP,12,iTUVQ) + Xab*AuxCont(iP,6,iTUVQ) 
     ThetaP(iP,7,2,iTUVQ) = AuxCont(iP,13,iTUVQ) + Xab*AuxCont(iP,7,iTUVQ) 
     ThetaP(iP,8,2,iTUVQ) = AuxCont(iP,14,iTUVQ) + Xab*AuxCont(iP,8,iTUVQ) 
     ThetaP(iP,9,2,iTUVQ) = AuxCont(iP,15,iTUVQ) + Xab*AuxCont(iP,9,iTUVQ) 
     ThetaP(iP,10,2,iTUVQ) = AuxCont(iP,16,iTUVQ) + Xab*AuxCont(iP,10,iTUVQ) 
     ThetaP(iP,5,3,iTUVQ) = AuxCont(iP,12,iTUVQ) + Yab*AuxCont(iP,5,iTUVQ) 
     ThetaP(iP,6,3,iTUVQ) = AuxCont(iP,14,iTUVQ) + Yab*AuxCont(iP,6,iTUVQ) 
     ThetaP(iP,7,3,iTUVQ) = AuxCont(iP,15,iTUVQ) + Yab*AuxCont(iP,7,iTUVQ) 
     ThetaP(iP,8,3,iTUVQ) = AuxCont(iP,17,iTUVQ) + Yab*AuxCont(iP,8,iTUVQ) 
     ThetaP(iP,9,3,iTUVQ) = AuxCont(iP,18,iTUVQ) + Yab*AuxCont(iP,9,iTUVQ) 
     ThetaP(iP,10,3,iTUVQ) = AuxCont(iP,19,iTUVQ) + Yab*AuxCont(iP,10,iTUVQ) 
     ThetaP(iP,5,4,iTUVQ) = AuxCont(iP,13,iTUVQ) + Zab*AuxCont(iP,5,iTUVQ) 
     ThetaP(iP,6,4,iTUVQ) = AuxCont(iP,15,iTUVQ) + Zab*AuxCont(iP,6,iTUVQ) 
     ThetaP(iP,7,4,iTUVQ) = AuxCont(iP,16,iTUVQ) + Zab*AuxCont(iP,7,iTUVQ) 
     ThetaP(iP,8,4,iTUVQ) = AuxCont(iP,18,iTUVQ) + Zab*AuxCont(iP,8,iTUVQ) 
     ThetaP(iP,9,4,iTUVQ) = AuxCont(iP,19,iTUVQ) + Zab*AuxCont(iP,9,iTUVQ) 
     ThetaP(iP,10,4,iTUVQ) = AuxCont(iP,20,iTUVQ) + Zab*AuxCont(iP,10,iTUVQ) 
  ENDDO
   ENDDO
end subroutine HorizontalRR_GPU_LHS_P3A2B1AtoB

subroutine HorizontalRR_GPU_LHS_P4A2B2AtoB(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(nContQP*nPasses,   35,nTUVQ)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(nContQP*nPasses,    5:   10,    5:   10,nTUVQ)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB
  real(realk) :: Xab,Yab,Zab
  real(realk) :: Tmp1(  5: 20,  2:  4)
!  real(realk) :: Tmp(nTUVA,nTUVB) ordering
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iP,iTUVQ,&
!$ACC         Tmp1,&
!$ACC         iPassP,iTUVA,iAtomA,iAtomB,Xab,Yab,Zab) &
!$ACC PRESENT(nPasses,&
!$ACC         iAtomApass,iAtomBpass,Pdistance12,AuxCont,ThetaP)
  DO iP = 1,nContQP*nPasses
   DO iTUVQ = 1,nTUVQ
    iPassP = (IP-1)/(nContQP) + 1
    iAtomA = iAtomApass(iPassP)
    iAtomB = iAtomBpass(iPassP)
    Xab = Pdistance12(1,iAtomA,iAtomB)
    Yab = Pdistance12(2,iAtomA,iAtomB)
    Zab = Pdistance12(3,iAtomA,iAtomB)
     Tmp1(5,2) = AuxCont(iP,11,iTUVQ) + Xab*AuxCont(iP,5,iTUVQ) 
     Tmp1(6,2) = AuxCont(iP,12,iTUVQ) + Xab*AuxCont(iP,6,iTUVQ) 
     Tmp1(7,2) = AuxCont(iP,13,iTUVQ) + Xab*AuxCont(iP,7,iTUVQ) 
     Tmp1(8,2) = AuxCont(iP,14,iTUVQ) + Xab*AuxCont(iP,8,iTUVQ) 
     Tmp1(9,2) = AuxCont(iP,15,iTUVQ) + Xab*AuxCont(iP,9,iTUVQ) 
     Tmp1(10,2) = AuxCont(iP,16,iTUVQ) + Xab*AuxCont(iP,10,iTUVQ) 
     Tmp1(11,2) = AuxCont(iP,21,iTUVQ) + Xab*AuxCont(iP,11,iTUVQ) 
     Tmp1(12,2) = AuxCont(iP,22,iTUVQ) + Xab*AuxCont(iP,12,iTUVQ) 
     Tmp1(13,2) = AuxCont(iP,23,iTUVQ) + Xab*AuxCont(iP,13,iTUVQ) 
     Tmp1(14,2) = AuxCont(iP,24,iTUVQ) + Xab*AuxCont(iP,14,iTUVQ) 
     Tmp1(15,2) = AuxCont(iP,25,iTUVQ) + Xab*AuxCont(iP,15,iTUVQ) 
     Tmp1(16,2) = AuxCont(iP,26,iTUVQ) + Xab*AuxCont(iP,16,iTUVQ) 
     Tmp1(17,2) = AuxCont(iP,27,iTUVQ) + Xab*AuxCont(iP,17,iTUVQ) 
     Tmp1(18,2) = AuxCont(iP,28,iTUVQ) + Xab*AuxCont(iP,18,iTUVQ) 
     Tmp1(19,2) = AuxCont(iP,29,iTUVQ) + Xab*AuxCont(iP,19,iTUVQ) 
     Tmp1(20,2) = AuxCont(iP,30,iTUVQ) + Xab*AuxCont(iP,20,iTUVQ) 
     Tmp1(5,3) = AuxCont(iP,12,iTUVQ) + Yab*AuxCont(iP,5,iTUVQ) 
     Tmp1(6,3) = AuxCont(iP,14,iTUVQ) + Yab*AuxCont(iP,6,iTUVQ) 
     Tmp1(7,3) = AuxCont(iP,15,iTUVQ) + Yab*AuxCont(iP,7,iTUVQ) 
     Tmp1(8,3) = AuxCont(iP,17,iTUVQ) + Yab*AuxCont(iP,8,iTUVQ) 
     Tmp1(9,3) = AuxCont(iP,18,iTUVQ) + Yab*AuxCont(iP,9,iTUVQ) 
     Tmp1(10,3) = AuxCont(iP,19,iTUVQ) + Yab*AuxCont(iP,10,iTUVQ) 
     Tmp1(11,3) = AuxCont(iP,22,iTUVQ) + Yab*AuxCont(iP,11,iTUVQ) 
     Tmp1(12,3) = AuxCont(iP,24,iTUVQ) + Yab*AuxCont(iP,12,iTUVQ) 
     Tmp1(13,3) = AuxCont(iP,25,iTUVQ) + Yab*AuxCont(iP,13,iTUVQ) 
     Tmp1(14,3) = AuxCont(iP,27,iTUVQ) + Yab*AuxCont(iP,14,iTUVQ) 
     Tmp1(15,3) = AuxCont(iP,28,iTUVQ) + Yab*AuxCont(iP,15,iTUVQ) 
     Tmp1(16,3) = AuxCont(iP,29,iTUVQ) + Yab*AuxCont(iP,16,iTUVQ) 
     Tmp1(17,3) = AuxCont(iP,31,iTUVQ) + Yab*AuxCont(iP,17,iTUVQ) 
     Tmp1(18,3) = AuxCont(iP,32,iTUVQ) + Yab*AuxCont(iP,18,iTUVQ) 
     Tmp1(19,3) = AuxCont(iP,33,iTUVQ) + Yab*AuxCont(iP,19,iTUVQ) 
     Tmp1(20,3) = AuxCont(iP,34,iTUVQ) + Yab*AuxCont(iP,20,iTUVQ) 
     Tmp1(5,4) = AuxCont(iP,13,iTUVQ) + Zab*AuxCont(iP,5,iTUVQ) 
     Tmp1(6,4) = AuxCont(iP,15,iTUVQ) + Zab*AuxCont(iP,6,iTUVQ) 
     Tmp1(7,4) = AuxCont(iP,16,iTUVQ) + Zab*AuxCont(iP,7,iTUVQ) 
     Tmp1(8,4) = AuxCont(iP,18,iTUVQ) + Zab*AuxCont(iP,8,iTUVQ) 
     Tmp1(9,4) = AuxCont(iP,19,iTUVQ) + Zab*AuxCont(iP,9,iTUVQ) 
     Tmp1(10,4) = AuxCont(iP,20,iTUVQ) + Zab*AuxCont(iP,10,iTUVQ) 
     Tmp1(11,4) = AuxCont(iP,23,iTUVQ) + Zab*AuxCont(iP,11,iTUVQ) 
     Tmp1(12,4) = AuxCont(iP,25,iTUVQ) + Zab*AuxCont(iP,12,iTUVQ) 
     Tmp1(13,4) = AuxCont(iP,26,iTUVQ) + Zab*AuxCont(iP,13,iTUVQ) 
     Tmp1(14,4) = AuxCont(iP,28,iTUVQ) + Zab*AuxCont(iP,14,iTUVQ) 
     Tmp1(15,4) = AuxCont(iP,29,iTUVQ) + Zab*AuxCont(iP,15,iTUVQ) 
     Tmp1(16,4) = AuxCont(iP,30,iTUVQ) + Zab*AuxCont(iP,16,iTUVQ) 
     Tmp1(17,4) = AuxCont(iP,32,iTUVQ) + Zab*AuxCont(iP,17,iTUVQ) 
     Tmp1(18,4) = AuxCont(iP,33,iTUVQ) + Zab*AuxCont(iP,18,iTUVQ) 
     Tmp1(19,4) = AuxCont(iP,34,iTUVQ) + Zab*AuxCont(iP,19,iTUVQ) 
     Tmp1(20,4) = AuxCont(iP,35,iTUVQ) + Zab*AuxCont(iP,20,iTUVQ) 
     ThetaP(iP,5,5,iTUVQ) = Tmp1(11,2) + Xab*Tmp1(5,2) 
     ThetaP(iP,6,5,iTUVQ) = Tmp1(12,2) + Xab*Tmp1(6,2) 
     ThetaP(iP,7,5,iTUVQ) = Tmp1(13,2) + Xab*Tmp1(7,2) 
     ThetaP(iP,8,5,iTUVQ) = Tmp1(14,2) + Xab*Tmp1(8,2) 
     ThetaP(iP,9,5,iTUVQ) = Tmp1(15,2) + Xab*Tmp1(9,2) 
     ThetaP(iP,10,5,iTUVQ) = Tmp1(16,2) + Xab*Tmp1(10,2) 
     ThetaP(iP,5,6,iTUVQ) = Tmp1(11,3) + Xab*Tmp1(5,3) 
     ThetaP(iP,6,6,iTUVQ) = Tmp1(12,3) + Xab*Tmp1(6,3) 
     ThetaP(iP,7,6,iTUVQ) = Tmp1(13,3) + Xab*Tmp1(7,3) 
     ThetaP(iP,8,6,iTUVQ) = Tmp1(14,3) + Xab*Tmp1(8,3) 
     ThetaP(iP,9,6,iTUVQ) = Tmp1(15,3) + Xab*Tmp1(9,3) 
     ThetaP(iP,10,6,iTUVQ) = Tmp1(16,3) + Xab*Tmp1(10,3) 
     ThetaP(iP,5,7,iTUVQ) = Tmp1(11,4) + Xab*Tmp1(5,4) 
     ThetaP(iP,6,7,iTUVQ) = Tmp1(12,4) + Xab*Tmp1(6,4) 
     ThetaP(iP,7,7,iTUVQ) = Tmp1(13,4) + Xab*Tmp1(7,4) 
     ThetaP(iP,8,7,iTUVQ) = Tmp1(14,4) + Xab*Tmp1(8,4) 
     ThetaP(iP,9,7,iTUVQ) = Tmp1(15,4) + Xab*Tmp1(9,4) 
     ThetaP(iP,10,7,iTUVQ) = Tmp1(16,4) + Xab*Tmp1(10,4) 
     ThetaP(iP,5,8,iTUVQ) = Tmp1(12,3) + Yab*Tmp1(5,3) 
     ThetaP(iP,6,8,iTUVQ) = Tmp1(14,3) + Yab*Tmp1(6,3) 
     ThetaP(iP,7,8,iTUVQ) = Tmp1(15,3) + Yab*Tmp1(7,3) 
     ThetaP(iP,8,8,iTUVQ) = Tmp1(17,3) + Yab*Tmp1(8,3) 
     ThetaP(iP,9,8,iTUVQ) = Tmp1(18,3) + Yab*Tmp1(9,3) 
     ThetaP(iP,10,8,iTUVQ) = Tmp1(19,3) + Yab*Tmp1(10,3) 
     ThetaP(iP,5,9,iTUVQ) = Tmp1(12,4) + Yab*Tmp1(5,4) 
     ThetaP(iP,6,9,iTUVQ) = Tmp1(14,4) + Yab*Tmp1(6,4) 
     ThetaP(iP,7,9,iTUVQ) = Tmp1(15,4) + Yab*Tmp1(7,4) 
     ThetaP(iP,8,9,iTUVQ) = Tmp1(17,4) + Yab*Tmp1(8,4) 
     ThetaP(iP,9,9,iTUVQ) = Tmp1(18,4) + Yab*Tmp1(9,4) 
     ThetaP(iP,10,9,iTUVQ) = Tmp1(19,4) + Yab*Tmp1(10,4) 
     ThetaP(iP,5,10,iTUVQ) = Tmp1(13,4) + Zab*Tmp1(5,4) 
     ThetaP(iP,6,10,iTUVQ) = Tmp1(15,4) + Zab*Tmp1(6,4) 
     ThetaP(iP,7,10,iTUVQ) = Tmp1(16,4) + Zab*Tmp1(7,4) 
     ThetaP(iP,8,10,iTUVQ) = Tmp1(18,4) + Zab*Tmp1(8,4) 
     ThetaP(iP,9,10,iTUVQ) = Tmp1(19,4) + Zab*Tmp1(9,4) 
     ThetaP(iP,10,10,iTUVQ) = Tmp1(20,4) + Zab*Tmp1(10,4) 
  ENDDO
   ENDDO
end subroutine HorizontalRR_GPU_LHS_P4A2B2AtoB
end module
