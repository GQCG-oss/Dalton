MODULE AGC_CPU_OBS_TRMODAtoDSegP1
 use IchorPrecisionMod
  
 CONTAINS
 subroutine TransferRecurrenceCPUP1Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(   10,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(    4,    4,nPrimQ*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0(  4,  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimB,iPrimC,&
!$OMP         Tmp0,&
!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1,  4
    DO iTUVP=1,  4
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = Pdistance12(1,iAtomA,iAtomB)
   Yab = Pdistance12(2,iAtomA,iAtomB)
   Zab = Pdistance12(3,iAtomA,iAtomB)
    expP = Pexp(iPrimP)
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpQ
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpQ
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
     DO iTUVP=1,  4
      Tmp0(iTUVP,1) = Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jq = 1
     do iTUVP = 1,  4
      Tmp0(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,  4
      Tmp0(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,  4
      Tmp0(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(2,iPrimQ,iPrimP,iPassP)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(5,iPrimQ,iPrimP,iPassP)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(3,iPrimQ,iPrimP,iPassP)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(8,iPrimQ,iPrimP,iPassP)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(4,iPrimQ,iPrimP,iPassP)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(10,iPrimQ,iPrimP,iPassP)
     DO iTUVQ=1,  4
      DO iTUVP=1,  4
        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP1Q1AtoDSegP
 subroutine TransferRecurrenceCPUP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(   20,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   10,    4,nPrimQ*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 10,  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimB,iPrimC,&
!$OMP         Tmp0,&
!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1,  4
    DO iTUVP=1, 10
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = Pdistance12(1,iAtomA,iAtomB)
   Yab = Pdistance12(2,iAtomA,iAtomB)
   Zab = Pdistance12(3,iAtomA,iAtomB)
    expP = Pexp(iPrimP)
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpQ
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpQ
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
     DO iTUVP=1, 10
      Tmp0(iTUVP,1) = Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jq = 1
     do iTUVP = 1, 10
      Tmp0(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 10
      Tmp0(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 10
      Tmp0(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expQ*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expQ*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expQ*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expQ*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expQ*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expQ*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expQ*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expQ*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expQ*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(2,iPrimQ,iPrimP,iPassP)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(5,iPrimQ,iPrimP,iPassP)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(5,2) = Tmp0(5,2) + pinvq*Aux(11,iPrimQ,iPrimP,iPassP)
     Tmp0(6,2) = Tmp0(6,2) + pinvq*Aux(12,iPrimQ,iPrimP,iPassP)
     Tmp0(7,2) = Tmp0(7,2) + pinvq*Aux(13,iPrimQ,iPrimP,iPassP)
     Tmp0(8,2) = Tmp0(8,2) + pinvq*Aux(14,iPrimQ,iPrimP,iPassP)
     Tmp0(9,2) = Tmp0(9,2) + pinvq*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(10,2) = Tmp0(10,2) + pinvq*Aux(16,iPrimQ,iPrimP,iPassP)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(3,iPrimQ,iPrimP,iPassP)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(8,iPrimQ,iPrimP,iPassP)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(5,3) = Tmp0(5,3) + pinvq*Aux(12,iPrimQ,iPrimP,iPassP)
     Tmp0(6,3) = Tmp0(6,3) + pinvq*Aux(14,iPrimQ,iPrimP,iPassP)
     Tmp0(7,3) = Tmp0(7,3) + pinvq*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(8,3) = Tmp0(8,3) + pinvq*Aux(17,iPrimQ,iPrimP,iPassP)
     Tmp0(9,3) = Tmp0(9,3) + pinvq*Aux(18,iPrimQ,iPrimP,iPassP)
     Tmp0(10,3) = Tmp0(10,3) + pinvq*Aux(19,iPrimQ,iPrimP,iPassP)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(4,iPrimQ,iPrimP,iPassP)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(10,iPrimQ,iPrimP,iPassP)
     Tmp0(5,4) = Tmp0(5,4) + pinvq*Aux(13,iPrimQ,iPrimP,iPassP)
     Tmp0(6,4) = Tmp0(6,4) + pinvq*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(7,4) = Tmp0(7,4) + pinvq*Aux(16,iPrimQ,iPrimP,iPassP)
     Tmp0(8,4) = Tmp0(8,4) + pinvq*Aux(18,iPrimQ,iPrimP,iPassP)
     Tmp0(9,4) = Tmp0(9,4) + pinvq*Aux(19,iPrimQ,iPrimP,iPassP)
     Tmp0(10,4) = Tmp0(10,4) + pinvq*Aux(20,iPrimQ,iPrimP,iPassP)
     DO iTUVQ=1,  4
      DO iTUVP=1, 10
        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP2Q1AtoDSegP
 subroutine TransferRecurrenceCPUP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(   35,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   10,   10,nPrimQ*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 10, 10)
  real(realk) :: Tmp1( 11: 20,  2:  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimB,iPrimC,&
!$OMP         Tmp0,&
!$OMP         Tmp1,&
!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1, 10
    DO iTUVP=1, 10
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = Pdistance12(1,iAtomA,iAtomB)
   Yab = Pdistance12(2,iAtomA,iAtomB)
   Zab = Pdistance12(3,iAtomA,iAtomB)
    expP = Pexp(iPrimP)
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpQ
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpQ
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
     DO iTUVP=1, 10
      Tmp0(iTUVP,1) = Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jq = 1
     do iTUVP = 1, 10
      Tmp0(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 10
      Tmp0(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 10
      Tmp0(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  11, 20
      Tmp1(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  11, 20
      Tmp1(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  11, 20
      Tmp1(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expQ*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expQ*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expQ*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expQ*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expQ*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expQ*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expQ*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expQ*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expQ*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp1(11,2) = Tmp1(11,2) + 3*inv2expQ*Aux(5,iPrimQ,iPrimP,iPassP) 
     Tmp1(12,2) = Tmp1(12,2) + 2*inv2expQ*Aux(6,iPrimQ,iPrimP,iPassP) 
     Tmp1(13,2) = Tmp1(13,2) + 2*inv2expQ*Aux(7,iPrimQ,iPrimP,iPassP) 
     Tmp1(14,2) = Tmp1(14,2) + inv2expQ*Aux(8,iPrimQ,iPrimP,iPassP) 
     Tmp1(15,2) = Tmp1(15,2) + inv2expQ*Aux(9,iPrimQ,iPrimP,iPassP) 
     Tmp1(16,2) = Tmp1(16,2) + inv2expQ*Aux(10,iPrimQ,iPrimP,iPassP) 
     Tmp1(12,3) = Tmp1(12,3) + inv2expQ*Aux(5,iPrimQ,iPrimP,iPassP) 
     Tmp1(14,3) = Tmp1(14,3) + 2*inv2expQ*Aux(6,iPrimQ,iPrimP,iPassP) 
     Tmp1(15,3) = Tmp1(15,3) + inv2expQ*Aux(7,iPrimQ,iPrimP,iPassP) 
     Tmp1(17,3) = Tmp1(17,3) + 3*inv2expQ*Aux(8,iPrimQ,iPrimP,iPassP) 
     Tmp1(18,3) = Tmp1(18,3) + 2*inv2expQ*Aux(9,iPrimQ,iPrimP,iPassP) 
     Tmp1(19,3) = Tmp1(19,3) + inv2expQ*Aux(10,iPrimQ,iPrimP,iPassP) 
     Tmp1(13,4) = Tmp1(13,4) + inv2expQ*Aux(5,iPrimQ,iPrimP,iPassP) 
     Tmp1(15,4) = Tmp1(15,4) + inv2expQ*Aux(6,iPrimQ,iPrimP,iPassP) 
     Tmp1(16,4) = Tmp1(16,4) + 2*inv2expQ*Aux(7,iPrimQ,iPrimP,iPassP) 
     Tmp1(18,4) = Tmp1(18,4) + inv2expQ*Aux(8,iPrimQ,iPrimP,iPassP) 
     Tmp1(19,4) = Tmp1(19,4) + 2*inv2expQ*Aux(9,iPrimQ,iPrimP,iPassP) 
     Tmp1(20,4) = Tmp1(20,4) + 3*inv2expQ*Aux(10,iPrimQ,iPrimP,iPassP) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(2,iPrimQ,iPrimP,iPassP)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(5,iPrimQ,iPrimP,iPassP)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(5,2) = Tmp0(5,2) + pinvq*Aux(11,iPrimQ,iPrimP,iPassP)
     Tmp0(6,2) = Tmp0(6,2) + pinvq*Aux(12,iPrimQ,iPrimP,iPassP)
     Tmp0(7,2) = Tmp0(7,2) + pinvq*Aux(13,iPrimQ,iPrimP,iPassP)
     Tmp0(8,2) = Tmp0(8,2) + pinvq*Aux(14,iPrimQ,iPrimP,iPassP)
     Tmp0(9,2) = Tmp0(9,2) + pinvq*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(10,2) = Tmp0(10,2) + pinvq*Aux(16,iPrimQ,iPrimP,iPassP)
     Tmp1(11,2) = Tmp1(11,2) + pinvq*Aux(21,iPrimQ,iPrimP,iPassP)
     Tmp1(12,2) = Tmp1(12,2) + pinvq*Aux(22,iPrimQ,iPrimP,iPassP)
     Tmp1(13,2) = Tmp1(13,2) + pinvq*Aux(23,iPrimQ,iPrimP,iPassP)
     Tmp1(14,2) = Tmp1(14,2) + pinvq*Aux(24,iPrimQ,iPrimP,iPassP)
     Tmp1(15,2) = Tmp1(15,2) + pinvq*Aux(25,iPrimQ,iPrimP,iPassP)
     Tmp1(16,2) = Tmp1(16,2) + pinvq*Aux(26,iPrimQ,iPrimP,iPassP)
     Tmp1(17,2) = Tmp1(17,2) + pinvq*Aux(27,iPrimQ,iPrimP,iPassP)
     Tmp1(18,2) = Tmp1(18,2) + pinvq*Aux(28,iPrimQ,iPrimP,iPassP)
     Tmp1(19,2) = Tmp1(19,2) + pinvq*Aux(29,iPrimQ,iPrimP,iPassP)
     Tmp1(20,2) = Tmp1(20,2) + pinvq*Aux(30,iPrimQ,iPrimP,iPassP)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(3,iPrimQ,iPrimP,iPassP)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(8,iPrimQ,iPrimP,iPassP)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(5,3) = Tmp0(5,3) + pinvq*Aux(12,iPrimQ,iPrimP,iPassP)
     Tmp0(6,3) = Tmp0(6,3) + pinvq*Aux(14,iPrimQ,iPrimP,iPassP)
     Tmp0(7,3) = Tmp0(7,3) + pinvq*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(8,3) = Tmp0(8,3) + pinvq*Aux(17,iPrimQ,iPrimP,iPassP)
     Tmp0(9,3) = Tmp0(9,3) + pinvq*Aux(18,iPrimQ,iPrimP,iPassP)
     Tmp0(10,3) = Tmp0(10,3) + pinvq*Aux(19,iPrimQ,iPrimP,iPassP)
     Tmp1(11,3) = Tmp1(11,3) + pinvq*Aux(22,iPrimQ,iPrimP,iPassP)
     Tmp1(12,3) = Tmp1(12,3) + pinvq*Aux(24,iPrimQ,iPrimP,iPassP)
     Tmp1(13,3) = Tmp1(13,3) + pinvq*Aux(25,iPrimQ,iPrimP,iPassP)
     Tmp1(14,3) = Tmp1(14,3) + pinvq*Aux(27,iPrimQ,iPrimP,iPassP)
     Tmp1(15,3) = Tmp1(15,3) + pinvq*Aux(28,iPrimQ,iPrimP,iPassP)
     Tmp1(16,3) = Tmp1(16,3) + pinvq*Aux(29,iPrimQ,iPrimP,iPassP)
     Tmp1(17,3) = Tmp1(17,3) + pinvq*Aux(31,iPrimQ,iPrimP,iPassP)
     Tmp1(18,3) = Tmp1(18,3) + pinvq*Aux(32,iPrimQ,iPrimP,iPassP)
     Tmp1(19,3) = Tmp1(19,3) + pinvq*Aux(33,iPrimQ,iPrimP,iPassP)
     Tmp1(20,3) = Tmp1(20,3) + pinvq*Aux(34,iPrimQ,iPrimP,iPassP)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(4,iPrimQ,iPrimP,iPassP)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(10,iPrimQ,iPrimP,iPassP)
     Tmp0(5,4) = Tmp0(5,4) + pinvq*Aux(13,iPrimQ,iPrimP,iPassP)
     Tmp0(6,4) = Tmp0(6,4) + pinvq*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(7,4) = Tmp0(7,4) + pinvq*Aux(16,iPrimQ,iPrimP,iPassP)
     Tmp0(8,4) = Tmp0(8,4) + pinvq*Aux(18,iPrimQ,iPrimP,iPassP)
     Tmp0(9,4) = Tmp0(9,4) + pinvq*Aux(19,iPrimQ,iPrimP,iPassP)
     Tmp0(10,4) = Tmp0(10,4) + pinvq*Aux(20,iPrimQ,iPrimP,iPassP)
     Tmp1(11,4) = Tmp1(11,4) + pinvq*Aux(23,iPrimQ,iPrimP,iPassP)
     Tmp1(12,4) = Tmp1(12,4) + pinvq*Aux(25,iPrimQ,iPrimP,iPassP)
     Tmp1(13,4) = Tmp1(13,4) + pinvq*Aux(26,iPrimQ,iPrimP,iPassP)
     Tmp1(14,4) = Tmp1(14,4) + pinvq*Aux(28,iPrimQ,iPrimP,iPassP)
     Tmp1(15,4) = Tmp1(15,4) + pinvq*Aux(29,iPrimQ,iPrimP,iPassP)
     Tmp1(16,4) = Tmp1(16,4) + pinvq*Aux(30,iPrimQ,iPrimP,iPassP)
     Tmp1(17,4) = Tmp1(17,4) + pinvq*Aux(32,iPrimQ,iPrimP,iPassP)
     Tmp1(18,4) = Tmp1(18,4) + pinvq*Aux(33,iPrimQ,iPrimP,iPassP)
     Tmp1(19,4) = Tmp1(19,4) + pinvq*Aux(34,iPrimQ,iPrimP,iPassP)
     Tmp1(20,4) = Tmp1(20,4) + pinvq*Aux(35,iPrimQ,iPrimP,iPassP)
 ! Building for Angular momentum Jq = 2
     do iTUVP = 1, 10
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 10
      Tmp0(iTUVP,6) = facX*Tmp0(iTUVP,3)
     enddo
     do iTUVP = 1, 10
      Tmp0(iTUVP,7) = facX*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 10
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 10
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 10
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     Tmp0(2,5) = Tmp0(2,5) + inv2expQ*Tmp0(1,2) 
     Tmp0(5,5) = Tmp0(5,5) + 2*inv2expQ*Tmp0(2,2) 
     Tmp0(6,5) = Tmp0(6,5) + inv2expQ*Tmp0(3,2) 
     Tmp0(7,5) = Tmp0(7,5) + inv2expQ*Tmp0(4,2) 
     Tmp0(2,6) = Tmp0(2,6) + inv2expQ*Tmp0(1,3) 
     Tmp0(5,6) = Tmp0(5,6) + 2*inv2expQ*Tmp0(2,3) 
     Tmp0(6,6) = Tmp0(6,6) + inv2expQ*Tmp0(3,3) 
     Tmp0(7,6) = Tmp0(7,6) + inv2expQ*Tmp0(4,3) 
     Tmp0(2,7) = Tmp0(2,7) + inv2expQ*Tmp0(1,4) 
     Tmp0(5,7) = Tmp0(5,7) + 2*inv2expQ*Tmp0(2,4) 
     Tmp0(6,7) = Tmp0(6,7) + inv2expQ*Tmp0(3,4) 
     Tmp0(7,7) = Tmp0(7,7) + inv2expQ*Tmp0(4,4) 
     Tmp0(3,8) = Tmp0(3,8) + inv2expQ*Tmp0(1,3) 
     Tmp0(6,8) = Tmp0(6,8) + inv2expQ*Tmp0(2,3) 
     Tmp0(8,8) = Tmp0(8,8) + 2*inv2expQ*Tmp0(3,3) 
     Tmp0(9,8) = Tmp0(9,8) + inv2expQ*Tmp0(4,3) 
     Tmp0(3,9) = Tmp0(3,9) + inv2expQ*Tmp0(1,4) 
     Tmp0(6,9) = Tmp0(6,9) + inv2expQ*Tmp0(2,4) 
     Tmp0(8,9) = Tmp0(8,9) + 2*inv2expQ*Tmp0(3,4) 
     Tmp0(9,9) = Tmp0(9,9) + inv2expQ*Tmp0(4,4) 
     Tmp0(4,10) = Tmp0(4,10) + inv2expQ*Tmp0(1,4) 
     Tmp0(7,10) = Tmp0(7,10) + inv2expQ*Tmp0(2,4) 
     Tmp0(9,10) = Tmp0(9,10) + inv2expQ*Tmp0(3,4) 
     Tmp0(10,10) = Tmp0(10,10) + 2*inv2expQ*Tmp0(4,4) 
     Tmp0(1,5) = Tmp0(1,5) + pinvq*Tmp0(2,2)
     Tmp0(2,5) = Tmp0(2,5) + pinvq*Tmp0(5,2)
     Tmp0(3,5) = Tmp0(3,5) + pinvq*Tmp0(6,2)
     Tmp0(4,5) = Tmp0(4,5) + pinvq*Tmp0(7,2)
     Tmp0(5,5) = Tmp0(5,5) + pinvq*Tmp1(11,2)
     Tmp0(6,5) = Tmp0(6,5) + pinvq*Tmp1(12,2)
     Tmp0(7,5) = Tmp0(7,5) + pinvq*Tmp1(13,2)
     Tmp0(8,5) = Tmp0(8,5) + pinvq*Tmp1(14,2)
     Tmp0(9,5) = Tmp0(9,5) + pinvq*Tmp1(15,2)
     Tmp0(10,5) = Tmp0(10,5) + pinvq*Tmp1(16,2)
     Tmp0(1,6) = Tmp0(1,6) + pinvq*Tmp0(2,3)
     Tmp0(2,6) = Tmp0(2,6) + pinvq*Tmp0(5,3)
     Tmp0(3,6) = Tmp0(3,6) + pinvq*Tmp0(6,3)
     Tmp0(4,6) = Tmp0(4,6) + pinvq*Tmp0(7,3)
     Tmp0(5,6) = Tmp0(5,6) + pinvq*Tmp1(11,3)
     Tmp0(6,6) = Tmp0(6,6) + pinvq*Tmp1(12,3)
     Tmp0(7,6) = Tmp0(7,6) + pinvq*Tmp1(13,3)
     Tmp0(8,6) = Tmp0(8,6) + pinvq*Tmp1(14,3)
     Tmp0(9,6) = Tmp0(9,6) + pinvq*Tmp1(15,3)
     Tmp0(10,6) = Tmp0(10,6) + pinvq*Tmp1(16,3)
     Tmp0(1,7) = Tmp0(1,7) + pinvq*Tmp0(2,4)
     Tmp0(2,7) = Tmp0(2,7) + pinvq*Tmp0(5,4)
     Tmp0(3,7) = Tmp0(3,7) + pinvq*Tmp0(6,4)
     Tmp0(4,7) = Tmp0(4,7) + pinvq*Tmp0(7,4)
     Tmp0(5,7) = Tmp0(5,7) + pinvq*Tmp1(11,4)
     Tmp0(6,7) = Tmp0(6,7) + pinvq*Tmp1(12,4)
     Tmp0(7,7) = Tmp0(7,7) + pinvq*Tmp1(13,4)
     Tmp0(8,7) = Tmp0(8,7) + pinvq*Tmp1(14,4)
     Tmp0(9,7) = Tmp0(9,7) + pinvq*Tmp1(15,4)
     Tmp0(10,7) = Tmp0(10,7) + pinvq*Tmp1(16,4)
     Tmp0(1,8) = Tmp0(1,8) + pinvq*Tmp0(3,3)
     Tmp0(2,8) = Tmp0(2,8) + pinvq*Tmp0(6,3)
     Tmp0(3,8) = Tmp0(3,8) + pinvq*Tmp0(8,3)
     Tmp0(4,8) = Tmp0(4,8) + pinvq*Tmp0(9,3)
     Tmp0(5,8) = Tmp0(5,8) + pinvq*Tmp1(12,3)
     Tmp0(6,8) = Tmp0(6,8) + pinvq*Tmp1(14,3)
     Tmp0(7,8) = Tmp0(7,8) + pinvq*Tmp1(15,3)
     Tmp0(8,8) = Tmp0(8,8) + pinvq*Tmp1(17,3)
     Tmp0(9,8) = Tmp0(9,8) + pinvq*Tmp1(18,3)
     Tmp0(10,8) = Tmp0(10,8) + pinvq*Tmp1(19,3)
     Tmp0(1,9) = Tmp0(1,9) + pinvq*Tmp0(3,4)
     Tmp0(2,9) = Tmp0(2,9) + pinvq*Tmp0(6,4)
     Tmp0(3,9) = Tmp0(3,9) + pinvq*Tmp0(8,4)
     Tmp0(4,9) = Tmp0(4,9) + pinvq*Tmp0(9,4)
     Tmp0(5,9) = Tmp0(5,9) + pinvq*Tmp1(12,4)
     Tmp0(6,9) = Tmp0(6,9) + pinvq*Tmp1(14,4)
     Tmp0(7,9) = Tmp0(7,9) + pinvq*Tmp1(15,4)
     Tmp0(8,9) = Tmp0(8,9) + pinvq*Tmp1(17,4)
     Tmp0(9,9) = Tmp0(9,9) + pinvq*Tmp1(18,4)
     Tmp0(10,9) = Tmp0(10,9) + pinvq*Tmp1(19,4)
     Tmp0(1,10) = Tmp0(1,10) + pinvq*Tmp0(4,4)
     Tmp0(2,10) = Tmp0(2,10) + pinvq*Tmp0(7,4)
     Tmp0(3,10) = Tmp0(3,10) + pinvq*Tmp0(9,4)
     Tmp0(4,10) = Tmp0(4,10) + pinvq*Tmp0(10,4)
     Tmp0(5,10) = Tmp0(5,10) + pinvq*Tmp1(13,4)
     Tmp0(6,10) = Tmp0(6,10) + pinvq*Tmp1(15,4)
     Tmp0(7,10) = Tmp0(7,10) + pinvq*Tmp1(16,4)
     Tmp0(8,10) = Tmp0(8,10) + pinvq*Tmp1(18,4)
     Tmp0(9,10) = Tmp0(9,10) + pinvq*Tmp1(19,4)
     Tmp0(10,10) = Tmp0(10,10) + pinvq*Tmp1(20,4)
     DO iTUVQ=1, 10
      DO iTUVP=1, 10
        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP2Q2AtoDSegP
 subroutine TransferRecurrenceCPUP3Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(   35,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   20,    4,nPrimQ*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 20,  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimB,iPrimC,&
!$OMP         Tmp0,&
!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1,  4
    DO iTUVP=1, 20
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = Pdistance12(1,iAtomA,iAtomB)
   Yab = Pdistance12(2,iAtomA,iAtomB)
   Zab = Pdistance12(3,iAtomA,iAtomB)
    expP = Pexp(iPrimP)
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpQ
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpQ
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
     DO iTUVP=1, 20
      Tmp0(iTUVP,1) = Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jq = 1
     do iTUVP = 1, 20
      Tmp0(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expQ*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expQ*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expQ*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(11,2) = Tmp0(11,2) + 3*inv2expQ*Aux(5,iPrimQ,iPrimP,iPassP) 
     Tmp0(12,2) = Tmp0(12,2) + 2*inv2expQ*Aux(6,iPrimQ,iPrimP,iPassP) 
     Tmp0(13,2) = Tmp0(13,2) + 2*inv2expQ*Aux(7,iPrimQ,iPrimP,iPassP) 
     Tmp0(14,2) = Tmp0(14,2) + inv2expQ*Aux(8,iPrimQ,iPrimP,iPassP) 
     Tmp0(15,2) = Tmp0(15,2) + inv2expQ*Aux(9,iPrimQ,iPrimP,iPassP) 
     Tmp0(16,2) = Tmp0(16,2) + inv2expQ*Aux(10,iPrimQ,iPrimP,iPassP) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expQ*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expQ*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expQ*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(12,3) = Tmp0(12,3) + inv2expQ*Aux(5,iPrimQ,iPrimP,iPassP) 
     Tmp0(14,3) = Tmp0(14,3) + 2*inv2expQ*Aux(6,iPrimQ,iPrimP,iPassP) 
     Tmp0(15,3) = Tmp0(15,3) + inv2expQ*Aux(7,iPrimQ,iPrimP,iPassP) 
     Tmp0(17,3) = Tmp0(17,3) + 3*inv2expQ*Aux(8,iPrimQ,iPrimP,iPassP) 
     Tmp0(18,3) = Tmp0(18,3) + 2*inv2expQ*Aux(9,iPrimQ,iPrimP,iPassP) 
     Tmp0(19,3) = Tmp0(19,3) + inv2expQ*Aux(10,iPrimQ,iPrimP,iPassP) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expQ*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expQ*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expQ*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(13,4) = Tmp0(13,4) + inv2expQ*Aux(5,iPrimQ,iPrimP,iPassP) 
     Tmp0(15,4) = Tmp0(15,4) + inv2expQ*Aux(6,iPrimQ,iPrimP,iPassP) 
     Tmp0(16,4) = Tmp0(16,4) + 2*inv2expQ*Aux(7,iPrimQ,iPrimP,iPassP) 
     Tmp0(18,4) = Tmp0(18,4) + inv2expQ*Aux(8,iPrimQ,iPrimP,iPassP) 
     Tmp0(19,4) = Tmp0(19,4) + 2*inv2expQ*Aux(9,iPrimQ,iPrimP,iPassP) 
     Tmp0(20,4) = Tmp0(20,4) + 3*inv2expQ*Aux(10,iPrimQ,iPrimP,iPassP) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(2,iPrimQ,iPrimP,iPassP)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(5,iPrimQ,iPrimP,iPassP)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(5,2) = Tmp0(5,2) + pinvq*Aux(11,iPrimQ,iPrimP,iPassP)
     Tmp0(6,2) = Tmp0(6,2) + pinvq*Aux(12,iPrimQ,iPrimP,iPassP)
     Tmp0(7,2) = Tmp0(7,2) + pinvq*Aux(13,iPrimQ,iPrimP,iPassP)
     Tmp0(8,2) = Tmp0(8,2) + pinvq*Aux(14,iPrimQ,iPrimP,iPassP)
     Tmp0(9,2) = Tmp0(9,2) + pinvq*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(10,2) = Tmp0(10,2) + pinvq*Aux(16,iPrimQ,iPrimP,iPassP)
     Tmp0(11,2) = Tmp0(11,2) + pinvq*Aux(21,iPrimQ,iPrimP,iPassP)
     Tmp0(12,2) = Tmp0(12,2) + pinvq*Aux(22,iPrimQ,iPrimP,iPassP)
     Tmp0(13,2) = Tmp0(13,2) + pinvq*Aux(23,iPrimQ,iPrimP,iPassP)
     Tmp0(14,2) = Tmp0(14,2) + pinvq*Aux(24,iPrimQ,iPrimP,iPassP)
     Tmp0(15,2) = Tmp0(15,2) + pinvq*Aux(25,iPrimQ,iPrimP,iPassP)
     Tmp0(16,2) = Tmp0(16,2) + pinvq*Aux(26,iPrimQ,iPrimP,iPassP)
     Tmp0(17,2) = Tmp0(17,2) + pinvq*Aux(27,iPrimQ,iPrimP,iPassP)
     Tmp0(18,2) = Tmp0(18,2) + pinvq*Aux(28,iPrimQ,iPrimP,iPassP)
     Tmp0(19,2) = Tmp0(19,2) + pinvq*Aux(29,iPrimQ,iPrimP,iPassP)
     Tmp0(20,2) = Tmp0(20,2) + pinvq*Aux(30,iPrimQ,iPrimP,iPassP)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(3,iPrimQ,iPrimP,iPassP)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(8,iPrimQ,iPrimP,iPassP)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(5,3) = Tmp0(5,3) + pinvq*Aux(12,iPrimQ,iPrimP,iPassP)
     Tmp0(6,3) = Tmp0(6,3) + pinvq*Aux(14,iPrimQ,iPrimP,iPassP)
     Tmp0(7,3) = Tmp0(7,3) + pinvq*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(8,3) = Tmp0(8,3) + pinvq*Aux(17,iPrimQ,iPrimP,iPassP)
     Tmp0(9,3) = Tmp0(9,3) + pinvq*Aux(18,iPrimQ,iPrimP,iPassP)
     Tmp0(10,3) = Tmp0(10,3) + pinvq*Aux(19,iPrimQ,iPrimP,iPassP)
     Tmp0(11,3) = Tmp0(11,3) + pinvq*Aux(22,iPrimQ,iPrimP,iPassP)
     Tmp0(12,3) = Tmp0(12,3) + pinvq*Aux(24,iPrimQ,iPrimP,iPassP)
     Tmp0(13,3) = Tmp0(13,3) + pinvq*Aux(25,iPrimQ,iPrimP,iPassP)
     Tmp0(14,3) = Tmp0(14,3) + pinvq*Aux(27,iPrimQ,iPrimP,iPassP)
     Tmp0(15,3) = Tmp0(15,3) + pinvq*Aux(28,iPrimQ,iPrimP,iPassP)
     Tmp0(16,3) = Tmp0(16,3) + pinvq*Aux(29,iPrimQ,iPrimP,iPassP)
     Tmp0(17,3) = Tmp0(17,3) + pinvq*Aux(31,iPrimQ,iPrimP,iPassP)
     Tmp0(18,3) = Tmp0(18,3) + pinvq*Aux(32,iPrimQ,iPrimP,iPassP)
     Tmp0(19,3) = Tmp0(19,3) + pinvq*Aux(33,iPrimQ,iPrimP,iPassP)
     Tmp0(20,3) = Tmp0(20,3) + pinvq*Aux(34,iPrimQ,iPrimP,iPassP)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(4,iPrimQ,iPrimP,iPassP)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(10,iPrimQ,iPrimP,iPassP)
     Tmp0(5,4) = Tmp0(5,4) + pinvq*Aux(13,iPrimQ,iPrimP,iPassP)
     Tmp0(6,4) = Tmp0(6,4) + pinvq*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(7,4) = Tmp0(7,4) + pinvq*Aux(16,iPrimQ,iPrimP,iPassP)
     Tmp0(8,4) = Tmp0(8,4) + pinvq*Aux(18,iPrimQ,iPrimP,iPassP)
     Tmp0(9,4) = Tmp0(9,4) + pinvq*Aux(19,iPrimQ,iPrimP,iPassP)
     Tmp0(10,4) = Tmp0(10,4) + pinvq*Aux(20,iPrimQ,iPrimP,iPassP)
     Tmp0(11,4) = Tmp0(11,4) + pinvq*Aux(23,iPrimQ,iPrimP,iPassP)
     Tmp0(12,4) = Tmp0(12,4) + pinvq*Aux(25,iPrimQ,iPrimP,iPassP)
     Tmp0(13,4) = Tmp0(13,4) + pinvq*Aux(26,iPrimQ,iPrimP,iPassP)
     Tmp0(14,4) = Tmp0(14,4) + pinvq*Aux(28,iPrimQ,iPrimP,iPassP)
     Tmp0(15,4) = Tmp0(15,4) + pinvq*Aux(29,iPrimQ,iPrimP,iPassP)
     Tmp0(16,4) = Tmp0(16,4) + pinvq*Aux(30,iPrimQ,iPrimP,iPassP)
     Tmp0(17,4) = Tmp0(17,4) + pinvq*Aux(32,iPrimQ,iPrimP,iPassP)
     Tmp0(18,4) = Tmp0(18,4) + pinvq*Aux(33,iPrimQ,iPrimP,iPassP)
     Tmp0(19,4) = Tmp0(19,4) + pinvq*Aux(34,iPrimQ,iPrimP,iPassP)
     Tmp0(20,4) = Tmp0(20,4) + pinvq*Aux(35,iPrimQ,iPrimP,iPassP)
     DO iTUVQ=1,  4
      DO iTUVP=1, 20
        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP3Q1AtoDSegP
 subroutine TransferRecurrenceCPUP3Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  use AGC_CPU_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(   56,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   20,   10,nPrimQ*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 20, 10)
  real(realk) :: Tmp1( 21: 35,  2:  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimB,iPrimC,&
!$OMP         Tmp0,&
!$OMP         Tmp1,&
!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1, 10
    DO iTUVP=1, 20
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = Pdistance12(1,iAtomA,iAtomB)
   Yab = Pdistance12(2,iAtomA,iAtomB)
   Zab = Pdistance12(3,iAtomA,iAtomB)
    expP = Pexp(iPrimP)
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpQ
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpQ
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
     DO iTUVP=1, 20
      Tmp0(iTUVP,1) = Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jq = 1
     do iTUVP = 1, 20
      Tmp0(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  21, 35
      Tmp1(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  21, 35
      Tmp1(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  21, 35
      Tmp1(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_35(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1_20(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2_35(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2_20(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3_35(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3_20(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX1_35(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1_20(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX2_35(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2_20(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX3_35(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3_20(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_35(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_35(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_35(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_35(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_35(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_35(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
 ! Building for Angular momentum Jq = 2
     do iTUVP = 1, 20
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,6) = facX*Tmp0(iTUVP,3)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,7) = facX*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_35(ituvpminus1)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + IfacX1_20(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_35(ituvpminus1)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + IfacX1_20(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_35(ituvpminus1)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + IfacX1_20(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2_35(ituvpminus1)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + IfacX2_20(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2_35(ituvpminus1)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + IfacX2_20(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3_35(ituvpminus1)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + IfacX3_20(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1_35(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp0(iTUVplus1,2)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1_35(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1_35(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp0(iTUVplus1,3)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1_35(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1_35(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1_35(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2_35(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp0(iTUVplus1,3)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2_35(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2_35(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2_35(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX3_35(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX3_35(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     DO iTUVQ=1, 10
      DO iTUVP=1, 20
        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP3Q2AtoDSegP
 subroutine TransferRecurrenceCPUP4Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  use AGC_CPU_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(   56,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   35,    4,nPrimQ*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 35,  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimB,iPrimC,&
!$OMP         Tmp0,&
!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1,  4
    DO iTUVP=1, 35
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = Pdistance12(1,iAtomA,iAtomB)
   Yab = Pdistance12(2,iAtomA,iAtomB)
   Zab = Pdistance12(3,iAtomA,iAtomB)
    expP = Pexp(iPrimP)
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpQ
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpQ
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
     DO iTUVP=1, 35
      Tmp0(iTUVP,1) = Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jq = 1
     do iTUVP = 1, 35
      Tmp0(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_35(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1_20(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_35(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2_20(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_35(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3_20(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1_35(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2_35(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3_35(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     DO iTUVQ=1,  4
      DO iTUVP=1, 35
        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP4Q1AtoDSegP
 subroutine TransferRecurrenceCPUP3Q3AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  use AGC_CPU_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(   84,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   20,   20,nPrimQ*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 20, 20)
  real(realk) :: Tmp1( 21: 56,  2:  4)
  real(realk) :: Tmp2( 21: 35,  5: 10)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimB,iPrimC,&
!$OMP         Tmp0,&
!$OMP         Tmp1,&
!$OMP         Tmp2,&
!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1, 20
    DO iTUVP=1, 20
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = Pdistance12(1,iAtomA,iAtomB)
   Yab = Pdistance12(2,iAtomA,iAtomB)
   Zab = Pdistance12(3,iAtomA,iAtomB)
    expP = Pexp(iPrimP)
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpQ
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpQ
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
     DO iTUVP=1, 20
      Tmp0(iTUVP,1) = Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jq = 1
     do iTUVP = 1, 20
      Tmp0(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  21, 56
      Tmp1(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  21, 56
      Tmp1(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  21, 56
      Tmp1(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 11,35
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 11,35
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 11,35
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 21,56
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 21,56
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 21,56
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
 ! Building for Angular momentum Jq = 2
     do iTUVP = 1, 20
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,6) = facX*Tmp0(iTUVP,3)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,7) = facX*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  21, 35
      Tmp2(iTUVP,5) = facX*Tmp1(iTUVP,2)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  21, 35
      Tmp2(iTUVP,6) = facX*Tmp1(iTUVP,3)
     enddo
     do iTUVP =  21, 35
      Tmp2(iTUVP,7) = facX*Tmp1(iTUVP,4)
     enddo
     do iTUVP =  21, 35
      Tmp2(iTUVP,8) = facY*Tmp1(iTUVP,3)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  21, 35
      Tmp2(iTUVP,9) = facY*Tmp1(iTUVP,4)
     enddo
     do iTUVP =  21, 35
      Tmp2(iTUVP,10) = facZ*Tmp1(iTUVP,4)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + IfacX2_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + IfacX2_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + IfacX3_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + IfacX2_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + IfacX2_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + IfacX3_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp0(iTUVplus1,2)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp0(iTUVplus1,3)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp0(iTUVplus1,3)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
 ! Building for Angular momentum Jq = 3
     do iTUVP = 1, 20
      Tmp0(iTUVP,11) = facX*Tmp0(iTUVP,5)+2*inv2expQ*Tmp0(iTUVP,2)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,12) = facY*Tmp0(iTUVP,5)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,13) = facZ*Tmp0(iTUVP,5)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,14) = facX*Tmp0(iTUVP,8)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,15) = facX*Tmp0(iTUVP,9)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,16) = facX*Tmp0(iTUVP,10)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,17) = facY*Tmp0(iTUVP,8)+2*inv2expQ*Tmp0(iTUVP,3)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,18) = facZ*Tmp0(iTUVP,8)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,19) = facY*Tmp0(iTUVP,10)
     enddo
     do iTUVP = 1, 20
      Tmp0(iTUVP,20) = facZ*Tmp0(iTUVP,10)+2*inv2expQ*Tmp0(iTUVP,4)
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + IfacX2_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + IfacX3_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,9) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + IfacX2_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + IfacX3_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + IfacX2_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + IfacX3_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp0(iTUVplus1,5)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp0(iTUVplus1,5)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp0(iTUVplus1,5)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp0(iTUVplus1,8)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp0(iTUVplus1,9)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp2(iTUVplus1,9)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp0(iTUVplus1,10)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp2(iTUVplus1,10)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp0(iTUVplus1,8)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp0(iTUVplus1,8)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp0(iTUVplus1,10)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp2(iTUVplus1,10)
     enddo
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp0(iTUVplus1,10)
     enddo
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp2(iTUVplus1,10)
     enddo
     DO iTUVQ=1, 20
      DO iTUVP=1, 20
        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP3Q3AtoDSegP
END MODULE AGC_CPU_OBS_TRMODAtoDSegP1
