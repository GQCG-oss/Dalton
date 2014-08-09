MODULE AGC_CPU_OBS_HorizontalRecurrenceLHSModAtoB
 use IchorPrecisionModule
  
 CONTAINS

subroutine HorizontalRR_CPU_LHS_P1A1B0AtoB(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(    4,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(    2:    4,1,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iTUVA) 
  DO iP = 1,nTUVQ*nContQP*nPasses
     DO iTUVA=  2,  4
        ThetaP(iTUVA,1,IP) = AuxCont(iTUVA,IP)
     ENDDO
  ENDDO
!$OMP END DO
end subroutine HorizontalRR_CPU_LHS_P1A1B0AtoB

subroutine HorizontalRR_CPU_LHS_P2A1B1AtoB(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(   10,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(    2:    4,    2:    4,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB
  real(realk) :: Xab,Yab,Zab
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iPassP,iTUVA,iAtomA,iAtomB,Xab,Yab,Zab) 
  DO iP = 1,nTUVQ*nContQP*nPasses
!    iTUVQ = mod(IP-1,nTUVQ)+1
!    iContQP = mod((IP-(mod(IP-1,nTUVQ)+1))/nTUVQ,nContQP)+1
    iPassP = (IP-1)/(nTUVQ*nContQP) + 1
    iAtomA = iAtomApass(iPassP)
    iAtomB = iAtomBpass(iPassP)
    Xab = Pdistance12(1,iAtomA,iAtomB)
    Yab = Pdistance12(2,iAtomA,iAtomB)
    Zab = Pdistance12(3,iAtomA,iAtomB)
     ThetaP( 2, 2      ,IP) = AuxCont( 5,      IP) + Xab*AuxCont( 2,      IP) 
     ThetaP( 3, 2      ,IP) = AuxCont( 6,      IP) + Xab*AuxCont( 3,      IP) 
     ThetaP( 4, 2      ,IP) = AuxCont( 7,      IP) + Xab*AuxCont( 4,      IP) 
     ThetaP( 2, 3      ,IP) = AuxCont( 6,      IP) + Yab*AuxCont( 2,      IP) 
     ThetaP( 3, 3      ,IP) = AuxCont( 8,      IP) + Yab*AuxCont( 3,      IP) 
     ThetaP( 4, 3      ,IP) = AuxCont( 9,      IP) + Yab*AuxCont( 4,      IP) 
     ThetaP( 2, 4      ,IP) = AuxCont( 7,      IP) + Zab*AuxCont( 2,      IP) 
     ThetaP( 3, 4      ,IP) = AuxCont( 9,      IP) + Zab*AuxCont( 3,      IP) 
     ThetaP( 4, 4      ,IP) = AuxCont(10,      IP) + Zab*AuxCont( 4,      IP) 
  ENDDO
!$OMP END DO
end subroutine HorizontalRR_CPU_LHS_P2A1B1AtoB

subroutine HorizontalRR_CPU_LHS_P2A2B0AtoB(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(   10,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(    5:   10,1,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iTUVA) 
  DO iP = 1,nTUVQ*nContQP*nPasses
     DO iTUVA=  5, 10
        ThetaP(iTUVA,1,IP) = AuxCont(iTUVA,IP)
     ENDDO
  ENDDO
!$OMP END DO
end subroutine HorizontalRR_CPU_LHS_P2A2B0AtoB

subroutine HorizontalRR_CPU_LHS_P3A2B1AtoB(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(   20,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(    5:   10,    2:    4,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB
  real(realk) :: Xab,Yab,Zab
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         iPassP,iTUVA,iAtomA,iAtomB,Xab,Yab,Zab) 
  DO iP = 1,nTUVQ*nContQP*nPasses
!    iTUVQ = mod(IP-1,nTUVQ)+1
!    iContQP = mod((IP-(mod(IP-1,nTUVQ)+1))/nTUVQ,nContQP)+1
    iPassP = (IP-1)/(nTUVQ*nContQP) + 1
    iAtomA = iAtomApass(iPassP)
    iAtomB = iAtomBpass(iPassP)
    Xab = Pdistance12(1,iAtomA,iAtomB)
    Yab = Pdistance12(2,iAtomA,iAtomB)
    Zab = Pdistance12(3,iAtomA,iAtomB)
     ThetaP( 5, 2      ,IP) = AuxCont(11,      IP) + Xab*AuxCont( 5,      IP) 
     ThetaP( 6, 2      ,IP) = AuxCont(12,      IP) + Xab*AuxCont( 6,      IP) 
     ThetaP( 7, 2      ,IP) = AuxCont(13,      IP) + Xab*AuxCont( 7,      IP) 
     ThetaP( 8, 2      ,IP) = AuxCont(14,      IP) + Xab*AuxCont( 8,      IP) 
     ThetaP( 9, 2      ,IP) = AuxCont(15,      IP) + Xab*AuxCont( 9,      IP) 
     ThetaP(10, 2      ,IP) = AuxCont(16,      IP) + Xab*AuxCont(10,      IP) 
     ThetaP( 5, 3      ,IP) = AuxCont(12,      IP) + Yab*AuxCont( 5,      IP) 
     ThetaP( 6, 3      ,IP) = AuxCont(14,      IP) + Yab*AuxCont( 6,      IP) 
     ThetaP( 7, 3      ,IP) = AuxCont(15,      IP) + Yab*AuxCont( 7,      IP) 
     ThetaP( 8, 3      ,IP) = AuxCont(17,      IP) + Yab*AuxCont( 8,      IP) 
     ThetaP( 9, 3      ,IP) = AuxCont(18,      IP) + Yab*AuxCont( 9,      IP) 
     ThetaP(10, 3      ,IP) = AuxCont(19,      IP) + Yab*AuxCont(10,      IP) 
     ThetaP( 5, 4      ,IP) = AuxCont(13,      IP) + Zab*AuxCont( 5,      IP) 
     ThetaP( 6, 4      ,IP) = AuxCont(15,      IP) + Zab*AuxCont( 6,      IP) 
     ThetaP( 7, 4      ,IP) = AuxCont(16,      IP) + Zab*AuxCont( 7,      IP) 
     ThetaP( 8, 4      ,IP) = AuxCont(18,      IP) + Zab*AuxCont( 8,      IP) 
     ThetaP( 9, 4      ,IP) = AuxCont(19,      IP) + Zab*AuxCont( 9,      IP) 
     ThetaP(10, 4      ,IP) = AuxCont(20,      IP) + Zab*AuxCont(10,      IP) 
  ENDDO
!$OMP END DO
end subroutine HorizontalRR_CPU_LHS_P3A2B1AtoB

subroutine HorizontalRR_CPU_LHS_P4A2B2AtoB(nContQP,nPasses,nTUVQ,&
         & Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,AuxCont,ThetaP,lupri)
  implicit none
  integer,intent(in) :: nContQP,nPasses,nTUVQ,lupri,MaxPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in) :: AuxCont(   35,nTUVQ*nContQP*nPasses)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(inout) :: ThetaP(    5:   10,    5:   10,nTUVQ*nContQP*nPasses)
  !Local variables
  integer :: iPassP,iP,iTUVQ,iTUVA,iAtomA,iAtomB
  real(realk) :: Xab,Yab,Zab
  real(realk) :: Tmp1(  5: 20,  2:  4)
!  real(realk) :: Tmp(nTUVA,nTUVB) ordering
!$OMP DO &
!$OMP PRIVATE(iP,&
!$OMP         Tmp1,&
!$OMP         iPassP,iTUVA,iAtomA,iAtomB,Xab,Yab,Zab) 
  DO iP = 1,nTUVQ*nContQP*nPasses
!    iTUVQ = mod(IP-1,nTUVQ)+1
!    iContQP = mod((IP-(mod(IP-1,nTUVQ)+1))/nTUVQ,nContQP)+1
    iPassP = (IP-1)/(nTUVQ*nContQP) + 1
    iAtomA = iAtomApass(iPassP)
    iAtomB = iAtomBpass(iPassP)
    Xab = Pdistance12(1,iAtomA,iAtomB)
    Yab = Pdistance12(2,iAtomA,iAtomB)
    Zab = Pdistance12(3,iAtomA,iAtomB)
     Tmp1( 5, 2) = AuxCont(11,      IP) + Xab*AuxCont( 5,      IP) 
     Tmp1( 6, 2) = AuxCont(12,      IP) + Xab*AuxCont( 6,      IP) 
     Tmp1( 7, 2) = AuxCont(13,      IP) + Xab*AuxCont( 7,      IP) 
     Tmp1( 8, 2) = AuxCont(14,      IP) + Xab*AuxCont( 8,      IP) 
     Tmp1( 9, 2) = AuxCont(15,      IP) + Xab*AuxCont( 9,      IP) 
     Tmp1(10, 2) = AuxCont(16,      IP) + Xab*AuxCont(10,      IP) 
     Tmp1(11, 2) = AuxCont(21,      IP) + Xab*AuxCont(11,      IP) 
     Tmp1(12, 2) = AuxCont(22,      IP) + Xab*AuxCont(12,      IP) 
     Tmp1(13, 2) = AuxCont(23,      IP) + Xab*AuxCont(13,      IP) 
     Tmp1(14, 2) = AuxCont(24,      IP) + Xab*AuxCont(14,      IP) 
     Tmp1(15, 2) = AuxCont(25,      IP) + Xab*AuxCont(15,      IP) 
     Tmp1(16, 2) = AuxCont(26,      IP) + Xab*AuxCont(16,      IP) 
     Tmp1(17, 2) = AuxCont(27,      IP) + Xab*AuxCont(17,      IP) 
     Tmp1(18, 2) = AuxCont(28,      IP) + Xab*AuxCont(18,      IP) 
     Tmp1(19, 2) = AuxCont(29,      IP) + Xab*AuxCont(19,      IP) 
     Tmp1(20, 2) = AuxCont(30,      IP) + Xab*AuxCont(20,      IP) 
     Tmp1( 5, 3) = AuxCont(12,      IP) + Yab*AuxCont( 5,      IP) 
     Tmp1( 6, 3) = AuxCont(14,      IP) + Yab*AuxCont( 6,      IP) 
     Tmp1( 7, 3) = AuxCont(15,      IP) + Yab*AuxCont( 7,      IP) 
     Tmp1( 8, 3) = AuxCont(17,      IP) + Yab*AuxCont( 8,      IP) 
     Tmp1( 9, 3) = AuxCont(18,      IP) + Yab*AuxCont( 9,      IP) 
     Tmp1(10, 3) = AuxCont(19,      IP) + Yab*AuxCont(10,      IP) 
     Tmp1(11, 3) = AuxCont(22,      IP) + Yab*AuxCont(11,      IP) 
     Tmp1(12, 3) = AuxCont(24,      IP) + Yab*AuxCont(12,      IP) 
     Tmp1(13, 3) = AuxCont(25,      IP) + Yab*AuxCont(13,      IP) 
     Tmp1(14, 3) = AuxCont(27,      IP) + Yab*AuxCont(14,      IP) 
     Tmp1(15, 3) = AuxCont(28,      IP) + Yab*AuxCont(15,      IP) 
     Tmp1(16, 3) = AuxCont(29,      IP) + Yab*AuxCont(16,      IP) 
     Tmp1(17, 3) = AuxCont(31,      IP) + Yab*AuxCont(17,      IP) 
     Tmp1(18, 3) = AuxCont(32,      IP) + Yab*AuxCont(18,      IP) 
     Tmp1(19, 3) = AuxCont(33,      IP) + Yab*AuxCont(19,      IP) 
     Tmp1(20, 3) = AuxCont(34,      IP) + Yab*AuxCont(20,      IP) 
     Tmp1( 5, 4) = AuxCont(13,      IP) + Zab*AuxCont( 5,      IP) 
     Tmp1( 6, 4) = AuxCont(15,      IP) + Zab*AuxCont( 6,      IP) 
     Tmp1( 7, 4) = AuxCont(16,      IP) + Zab*AuxCont( 7,      IP) 
     Tmp1( 8, 4) = AuxCont(18,      IP) + Zab*AuxCont( 8,      IP) 
     Tmp1( 9, 4) = AuxCont(19,      IP) + Zab*AuxCont( 9,      IP) 
     Tmp1(10, 4) = AuxCont(20,      IP) + Zab*AuxCont(10,      IP) 
     Tmp1(11, 4) = AuxCont(23,      IP) + Zab*AuxCont(11,      IP) 
     Tmp1(12, 4) = AuxCont(25,      IP) + Zab*AuxCont(12,      IP) 
     Tmp1(13, 4) = AuxCont(26,      IP) + Zab*AuxCont(13,      IP) 
     Tmp1(14, 4) = AuxCont(28,      IP) + Zab*AuxCont(14,      IP) 
     Tmp1(15, 4) = AuxCont(29,      IP) + Zab*AuxCont(15,      IP) 
     Tmp1(16, 4) = AuxCont(30,      IP) + Zab*AuxCont(16,      IP) 
     Tmp1(17, 4) = AuxCont(32,      IP) + Zab*AuxCont(17,      IP) 
     Tmp1(18, 4) = AuxCont(33,      IP) + Zab*AuxCont(18,      IP) 
     Tmp1(19, 4) = AuxCont(34,      IP) + Zab*AuxCont(19,      IP) 
     Tmp1(20, 4) = AuxCont(35,      IP) + Zab*AuxCont(20,      IP) 
     ThetaP( 5, 5      ,IP) = Tmp1(11, 2) + Xab*Tmp1( 5, 2) 
     ThetaP( 6, 5      ,IP) = Tmp1(12, 2) + Xab*Tmp1( 6, 2) 
     ThetaP( 7, 5      ,IP) = Tmp1(13, 2) + Xab*Tmp1( 7, 2) 
     ThetaP( 8, 5      ,IP) = Tmp1(14, 2) + Xab*Tmp1( 8, 2) 
     ThetaP( 9, 5      ,IP) = Tmp1(15, 2) + Xab*Tmp1( 9, 2) 
     ThetaP(10, 5      ,IP) = Tmp1(16, 2) + Xab*Tmp1(10, 2) 
     ThetaP( 5, 6      ,IP) = Tmp1(11, 3) + Xab*Tmp1( 5, 3) 
     ThetaP( 6, 6      ,IP) = Tmp1(12, 3) + Xab*Tmp1( 6, 3) 
     ThetaP( 7, 6      ,IP) = Tmp1(13, 3) + Xab*Tmp1( 7, 3) 
     ThetaP( 8, 6      ,IP) = Tmp1(14, 3) + Xab*Tmp1( 8, 3) 
     ThetaP( 9, 6      ,IP) = Tmp1(15, 3) + Xab*Tmp1( 9, 3) 
     ThetaP(10, 6      ,IP) = Tmp1(16, 3) + Xab*Tmp1(10, 3) 
     ThetaP( 5, 7      ,IP) = Tmp1(11, 4) + Xab*Tmp1( 5, 4) 
     ThetaP( 6, 7      ,IP) = Tmp1(12, 4) + Xab*Tmp1( 6, 4) 
     ThetaP( 7, 7      ,IP) = Tmp1(13, 4) + Xab*Tmp1( 7, 4) 
     ThetaP( 8, 7      ,IP) = Tmp1(14, 4) + Xab*Tmp1( 8, 4) 
     ThetaP( 9, 7      ,IP) = Tmp1(15, 4) + Xab*Tmp1( 9, 4) 
     ThetaP(10, 7      ,IP) = Tmp1(16, 4) + Xab*Tmp1(10, 4) 
     ThetaP( 5, 8      ,IP) = Tmp1(12, 3) + Yab*Tmp1( 5, 3) 
     ThetaP( 6, 8      ,IP) = Tmp1(14, 3) + Yab*Tmp1( 6, 3) 
     ThetaP( 7, 8      ,IP) = Tmp1(15, 3) + Yab*Tmp1( 7, 3) 
     ThetaP( 8, 8      ,IP) = Tmp1(17, 3) + Yab*Tmp1( 8, 3) 
     ThetaP( 9, 8      ,IP) = Tmp1(18, 3) + Yab*Tmp1( 9, 3) 
     ThetaP(10, 8      ,IP) = Tmp1(19, 3) + Yab*Tmp1(10, 3) 
     ThetaP( 5, 9      ,IP) = Tmp1(12, 4) + Yab*Tmp1( 5, 4) 
     ThetaP( 6, 9      ,IP) = Tmp1(14, 4) + Yab*Tmp1( 6, 4) 
     ThetaP( 7, 9      ,IP) = Tmp1(15, 4) + Yab*Tmp1( 7, 4) 
     ThetaP( 8, 9      ,IP) = Tmp1(17, 4) + Yab*Tmp1( 8, 4) 
     ThetaP( 9, 9      ,IP) = Tmp1(18, 4) + Yab*Tmp1( 9, 4) 
     ThetaP(10, 9      ,IP) = Tmp1(19, 4) + Yab*Tmp1(10, 4) 
     ThetaP( 5,10      ,IP) = Tmp1(13, 4) + Zab*Tmp1( 5, 4) 
     ThetaP( 6,10      ,IP) = Tmp1(15, 4) + Zab*Tmp1( 6, 4) 
     ThetaP( 7,10      ,IP) = Tmp1(16, 4) + Zab*Tmp1( 7, 4) 
     ThetaP( 8,10      ,IP) = Tmp1(18, 4) + Zab*Tmp1( 8, 4) 
     ThetaP( 9,10      ,IP) = Tmp1(19, 4) + Zab*Tmp1( 9, 4) 
     ThetaP(10,10      ,IP) = Tmp1(20, 4) + Zab*Tmp1(10, 4) 
  ENDDO
!$OMP END DO
end subroutine HorizontalRR_CPU_LHS_P4A2B2AtoB
end module
