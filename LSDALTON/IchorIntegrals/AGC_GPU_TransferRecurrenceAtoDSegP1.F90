MODULE AGC_GPU_OBS_TRMODAtoDSegP1
 use IchorPrecisionMod
  
 CONTAINS
 subroutine TransferRecurrenceGPUP1Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   10)
  real(realk),intent(inout) :: Aux2(nPrimQ*nPasses,    4,    4)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0(  4,  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,nPrimP,Aux2) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1,  4
    DO iTUVP=1,  4
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         iPrimB,iPrimC,&
!$ACC         Tmp0,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
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
!$ACC LOOP SEQ
     DO iTUVP=1,  4
      Tmp0(iTUVP,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1,  4
      Tmp0(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,  4
      Tmp0(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,  4
      Tmp0(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,2)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,5)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,3)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,8)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,4)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,10)
!$ACC LOOP SEQ
     DO iTUVQ=1,  4
!$ACC LOOP SEQ
      DO iTUVP=1,  4
        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
 end subroutine TransferRecurrenceGPUP1Q1AtoDSegP
 subroutine TransferRecurrenceGPUP2Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   20)
  real(realk),intent(inout) :: Aux2(nPrimQ*nPasses,   10,    4)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 10,  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,nPrimP,Aux2) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1,  4
    DO iTUVP=1, 10
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         iPrimB,iPrimC,&
!$ACC         Tmp0,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
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
!$ACC LOOP SEQ
     DO iTUVP=1, 10
      Tmp0(iTUVP,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,2)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,5)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(5,2) = Tmp0(5,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,11)
     Tmp0(6,2) = Tmp0(6,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(7,2) = Tmp0(7,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(8,2) = Tmp0(8,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(9,2) = Tmp0(9,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(10,2) = Tmp0(10,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,3)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,8)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(5,3) = Tmp0(5,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(6,3) = Tmp0(6,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(7,3) = Tmp0(7,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(8,3) = Tmp0(8,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,17)
     Tmp0(9,3) = Tmp0(9,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(10,3) = Tmp0(10,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,4)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,10)
     Tmp0(5,4) = Tmp0(5,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(6,4) = Tmp0(6,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(7,4) = Tmp0(7,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp0(8,4) = Tmp0(8,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(9,4) = Tmp0(9,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp0(10,4) = Tmp0(10,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,20)
!$ACC LOOP SEQ
     DO iTUVQ=1,  4
!$ACC LOOP SEQ
      DO iTUVP=1, 10
        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
 end subroutine TransferRecurrenceGPUP2Q1AtoDSegP
 subroutine TransferRecurrenceGPUP2Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   35)
  real(realk),intent(inout) :: Aux2(nPrimQ*nPasses,   10,   10)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
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
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,nPrimP,Aux2) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1, 10
    DO iTUVP=1, 10
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         iPrimB,iPrimC,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
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
!$ACC LOOP SEQ
     DO iTUVP=1, 10
      Tmp0(iTUVP,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  11, 20
      Tmp1(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  11, 20
      Tmp1(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  11, 20
      Tmp1(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp1(11,2) = Tmp1(11,2) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp1(12,2) = Tmp1(12,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp1(13,2) = Tmp1(13,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp1(14,2) = Tmp1(14,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp1(15,2) = Tmp1(15,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp1(16,2) = Tmp1(16,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp1(12,3) = Tmp1(12,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp1(14,3) = Tmp1(14,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp1(15,3) = Tmp1(15,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp1(17,3) = Tmp1(17,3) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp1(18,3) = Tmp1(18,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp1(19,3) = Tmp1(19,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp1(13,4) = Tmp1(13,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp1(15,4) = Tmp1(15,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp1(16,4) = Tmp1(16,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp1(18,4) = Tmp1(18,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp1(19,4) = Tmp1(19,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp1(20,4) = Tmp1(20,4) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,2)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,5)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(5,2) = Tmp0(5,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,11)
     Tmp0(6,2) = Tmp0(6,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(7,2) = Tmp0(7,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(8,2) = Tmp0(8,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(9,2) = Tmp0(9,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(10,2) = Tmp0(10,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp1(11,2) = Tmp1(11,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,21)
     Tmp1(12,2) = Tmp1(12,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,22)
     Tmp1(13,2) = Tmp1(13,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,23)
     Tmp1(14,2) = Tmp1(14,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,24)
     Tmp1(15,2) = Tmp1(15,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp1(16,2) = Tmp1(16,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,26)
     Tmp1(17,2) = Tmp1(17,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,27)
     Tmp1(18,2) = Tmp1(18,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp1(19,2) = Tmp1(19,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp1(20,2) = Tmp1(20,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,30)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,3)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,8)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(5,3) = Tmp0(5,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(6,3) = Tmp0(6,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(7,3) = Tmp0(7,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(8,3) = Tmp0(8,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,17)
     Tmp0(9,3) = Tmp0(9,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(10,3) = Tmp0(10,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp1(11,3) = Tmp1(11,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,22)
     Tmp1(12,3) = Tmp1(12,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,24)
     Tmp1(13,3) = Tmp1(13,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp1(14,3) = Tmp1(14,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,27)
     Tmp1(15,3) = Tmp1(15,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp1(16,3) = Tmp1(16,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp1(17,3) = Tmp1(17,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,31)
     Tmp1(18,3) = Tmp1(18,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,32)
     Tmp1(19,3) = Tmp1(19,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,33)
     Tmp1(20,3) = Tmp1(20,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,34)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,4)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,10)
     Tmp0(5,4) = Tmp0(5,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(6,4) = Tmp0(6,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(7,4) = Tmp0(7,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp0(8,4) = Tmp0(8,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(9,4) = Tmp0(9,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp0(10,4) = Tmp0(10,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,20)
     Tmp1(11,4) = Tmp1(11,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,23)
     Tmp1(12,4) = Tmp1(12,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp1(13,4) = Tmp1(13,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,26)
     Tmp1(14,4) = Tmp1(14,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp1(15,4) = Tmp1(15,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp1(16,4) = Tmp1(16,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,30)
     Tmp1(17,4) = Tmp1(17,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,32)
     Tmp1(18,4) = Tmp1(18,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,33)
     Tmp1(19,4) = Tmp1(19,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,34)
     Tmp1(20,4) = Tmp1(20,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,35)
 ! Building for Angular momentum Jq = 2
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,6) = facX*Tmp0(iTUVP,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,7) = facX*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
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
!$ACC LOOP SEQ
     DO iTUVQ=1, 10
!$ACC LOOP SEQ
      DO iTUVP=1, 10
        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
 end subroutine TransferRecurrenceGPUP2Q2AtoDSegP
 subroutine TransferRecurrenceGPUP3Q1AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   35)
  real(realk),intent(inout) :: Aux2(nPrimQ*nPasses,   20,    4)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 20,  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,nPrimP,Aux2) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1,  4
    DO iTUVP=1, 20
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         iPrimB,iPrimC,&
!$ACC         Tmp0,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
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
!$ACC LOOP SEQ
     DO iTUVP=1, 20
      Tmp0(iTUVP,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(11,2) = Tmp0(11,2) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp0(12,2) = Tmp0(12,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp0(13,2) = Tmp0(13,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp0(14,2) = Tmp0(14,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp0(15,2) = Tmp0(15,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp0(16,2) = Tmp0(16,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(12,3) = Tmp0(12,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp0(14,3) = Tmp0(14,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp0(15,3) = Tmp0(15,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp0(17,3) = Tmp0(17,3) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp0(18,3) = Tmp0(18,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp0(19,3) = Tmp0(19,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(13,4) = Tmp0(13,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp0(15,4) = Tmp0(15,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp0(16,4) = Tmp0(16,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp0(18,4) = Tmp0(18,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp0(19,4) = Tmp0(19,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp0(20,4) = Tmp0(20,4) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,2)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,5)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(5,2) = Tmp0(5,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,11)
     Tmp0(6,2) = Tmp0(6,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(7,2) = Tmp0(7,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(8,2) = Tmp0(8,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(9,2) = Tmp0(9,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(10,2) = Tmp0(10,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp0(11,2) = Tmp0(11,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,21)
     Tmp0(12,2) = Tmp0(12,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,22)
     Tmp0(13,2) = Tmp0(13,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,23)
     Tmp0(14,2) = Tmp0(14,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,24)
     Tmp0(15,2) = Tmp0(15,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp0(16,2) = Tmp0(16,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,26)
     Tmp0(17,2) = Tmp0(17,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,27)
     Tmp0(18,2) = Tmp0(18,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp0(19,2) = Tmp0(19,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp0(20,2) = Tmp0(20,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,30)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,3)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,8)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(5,3) = Tmp0(5,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(6,3) = Tmp0(6,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(7,3) = Tmp0(7,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(8,3) = Tmp0(8,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,17)
     Tmp0(9,3) = Tmp0(9,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(10,3) = Tmp0(10,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp0(11,3) = Tmp0(11,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,22)
     Tmp0(12,3) = Tmp0(12,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,24)
     Tmp0(13,3) = Tmp0(13,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp0(14,3) = Tmp0(14,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,27)
     Tmp0(15,3) = Tmp0(15,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp0(16,3) = Tmp0(16,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp0(17,3) = Tmp0(17,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,31)
     Tmp0(18,3) = Tmp0(18,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,32)
     Tmp0(19,3) = Tmp0(19,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,33)
     Tmp0(20,3) = Tmp0(20,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,34)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,4)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,10)
     Tmp0(5,4) = Tmp0(5,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(6,4) = Tmp0(6,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(7,4) = Tmp0(7,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp0(8,4) = Tmp0(8,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(9,4) = Tmp0(9,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp0(10,4) = Tmp0(10,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,20)
     Tmp0(11,4) = Tmp0(11,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,23)
     Tmp0(12,4) = Tmp0(12,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp0(13,4) = Tmp0(13,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,26)
     Tmp0(14,4) = Tmp0(14,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp0(15,4) = Tmp0(15,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp0(16,4) = Tmp0(16,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,30)
     Tmp0(17,4) = Tmp0(17,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,32)
     Tmp0(18,4) = Tmp0(18,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,33)
     Tmp0(19,4) = Tmp0(19,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,34)
     Tmp0(20,4) = Tmp0(20,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,35)
!$ACC LOOP SEQ
     DO iTUVQ=1,  4
!$ACC LOOP SEQ
      DO iTUVP=1, 20
        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
 end subroutine TransferRecurrenceGPUP3Q1AtoDSegP
 subroutine TransferRecurrenceGPUP3Q2AtoDSegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   56)
  real(realk),intent(inout) :: Aux2(nPrimQ*nPasses,   20,   10)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
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
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,nPrimP,Aux2) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1, 10
    DO iTUVP=1, 20
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         iPrimB,iPrimC,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
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
!$ACC LOOP SEQ
     DO iTUVP=1, 20
      Tmp0(iTUVP,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp1(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp1(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp1(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(11,2) = Tmp0(11,2) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp0(12,2) = Tmp0(12,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp0(13,2) = Tmp0(13,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp0(14,2) = Tmp0(14,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp0(15,2) = Tmp0(15,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp0(16,2) = Tmp0(16,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(12,3) = Tmp0(12,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp0(14,3) = Tmp0(14,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp0(15,3) = Tmp0(15,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp0(17,3) = Tmp0(17,3) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp0(18,3) = Tmp0(18,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp0(19,3) = Tmp0(19,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(13,4) = Tmp0(13,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp0(15,4) = Tmp0(15,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp0(16,4) = Tmp0(16,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp0(18,4) = Tmp0(18,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp0(19,4) = Tmp0(19,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp0(20,4) = Tmp0(20,4) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp1(21,2) = Tmp1(21,2) + 4*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,11) 
     Tmp1(22,2) = Tmp1(22,2) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,12) 
     Tmp1(23,2) = Tmp1(23,2) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,13) 
     Tmp1(24,2) = Tmp1(24,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,14) 
     Tmp1(25,2) = Tmp1(25,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,15) 
     Tmp1(26,2) = Tmp1(26,2) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,16) 
     Tmp1(27,2) = Tmp1(27,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,17) 
     Tmp1(28,2) = Tmp1(28,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,18) 
     Tmp1(29,2) = Tmp1(29,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,19) 
     Tmp1(30,2) = Tmp1(30,2) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,20) 
     Tmp1(22,3) = Tmp1(22,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,11) 
     Tmp1(24,3) = Tmp1(24,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,12) 
     Tmp1(25,3) = Tmp1(25,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,13) 
     Tmp1(27,3) = Tmp1(27,3) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,14) 
     Tmp1(28,3) = Tmp1(28,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,15) 
     Tmp1(29,3) = Tmp1(29,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,16) 
     Tmp1(31,3) = Tmp1(31,3) + 4*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,17) 
     Tmp1(32,3) = Tmp1(32,3) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,18) 
     Tmp1(33,3) = Tmp1(33,3) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,19) 
     Tmp1(34,3) = Tmp1(34,3) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,20) 
     Tmp1(23,4) = Tmp1(23,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,11) 
     Tmp1(25,4) = Tmp1(25,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,12) 
     Tmp1(26,4) = Tmp1(26,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,13) 
     Tmp1(28,4) = Tmp1(28,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,14) 
     Tmp1(29,4) = Tmp1(29,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,15) 
     Tmp1(30,4) = Tmp1(30,4) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,16) 
     Tmp1(32,4) = Tmp1(32,4) + inv2expQ*Aux(iPrimQ,iPrimP,iPassP,17) 
     Tmp1(33,4) = Tmp1(33,4) + 2*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,18) 
     Tmp1(34,4) = Tmp1(34,4) + 3*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,19) 
     Tmp1(35,4) = Tmp1(35,4) + 4*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,20) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,2)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,5)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(5,2) = Tmp0(5,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,11)
     Tmp0(6,2) = Tmp0(6,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(7,2) = Tmp0(7,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(8,2) = Tmp0(8,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(9,2) = Tmp0(9,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(10,2) = Tmp0(10,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp0(11,2) = Tmp0(11,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,21)
     Tmp0(12,2) = Tmp0(12,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,22)
     Tmp0(13,2) = Tmp0(13,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,23)
     Tmp0(14,2) = Tmp0(14,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,24)
     Tmp0(15,2) = Tmp0(15,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp0(16,2) = Tmp0(16,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,26)
     Tmp0(17,2) = Tmp0(17,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,27)
     Tmp0(18,2) = Tmp0(18,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp0(19,2) = Tmp0(19,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp0(20,2) = Tmp0(20,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,30)
     Tmp1(21,2) = Tmp1(21,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,36)
     Tmp1(22,2) = Tmp1(22,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,37)
     Tmp1(23,2) = Tmp1(23,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,38)
     Tmp1(24,2) = Tmp1(24,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,39)
     Tmp1(25,2) = Tmp1(25,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,40)
     Tmp1(26,2) = Tmp1(26,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,41)
     Tmp1(27,2) = Tmp1(27,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,42)
     Tmp1(28,2) = Tmp1(28,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,43)
     Tmp1(29,2) = Tmp1(29,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,44)
     Tmp1(30,2) = Tmp1(30,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,45)
     Tmp1(31,2) = Tmp1(31,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,46)
     Tmp1(32,2) = Tmp1(32,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,47)
     Tmp1(33,2) = Tmp1(33,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,48)
     Tmp1(34,2) = Tmp1(34,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,49)
     Tmp1(35,2) = Tmp1(35,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,50)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,3)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,8)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(5,3) = Tmp0(5,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(6,3) = Tmp0(6,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(7,3) = Tmp0(7,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(8,3) = Tmp0(8,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,17)
     Tmp0(9,3) = Tmp0(9,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(10,3) = Tmp0(10,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp0(11,3) = Tmp0(11,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,22)
     Tmp0(12,3) = Tmp0(12,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,24)
     Tmp0(13,3) = Tmp0(13,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp0(14,3) = Tmp0(14,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,27)
     Tmp0(15,3) = Tmp0(15,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp0(16,3) = Tmp0(16,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp0(17,3) = Tmp0(17,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,31)
     Tmp0(18,3) = Tmp0(18,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,32)
     Tmp0(19,3) = Tmp0(19,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,33)
     Tmp0(20,3) = Tmp0(20,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,34)
     Tmp1(21,3) = Tmp1(21,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,37)
     Tmp1(22,3) = Tmp1(22,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,39)
     Tmp1(23,3) = Tmp1(23,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,40)
     Tmp1(24,3) = Tmp1(24,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,42)
     Tmp1(25,3) = Tmp1(25,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,43)
     Tmp1(26,3) = Tmp1(26,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,44)
     Tmp1(27,3) = Tmp1(27,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,46)
     Tmp1(28,3) = Tmp1(28,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,47)
     Tmp1(29,3) = Tmp1(29,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,48)
     Tmp1(30,3) = Tmp1(30,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,49)
     Tmp1(31,3) = Tmp1(31,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,51)
     Tmp1(32,3) = Tmp1(32,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,52)
     Tmp1(33,3) = Tmp1(33,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,53)
     Tmp1(34,3) = Tmp1(34,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,54)
     Tmp1(35,3) = Tmp1(35,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,55)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,4)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,10)
     Tmp0(5,4) = Tmp0(5,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(6,4) = Tmp0(6,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(7,4) = Tmp0(7,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp0(8,4) = Tmp0(8,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(9,4) = Tmp0(9,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp0(10,4) = Tmp0(10,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,20)
     Tmp0(11,4) = Tmp0(11,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,23)
     Tmp0(12,4) = Tmp0(12,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp0(13,4) = Tmp0(13,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,26)
     Tmp0(14,4) = Tmp0(14,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp0(15,4) = Tmp0(15,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp0(16,4) = Tmp0(16,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,30)
     Tmp0(17,4) = Tmp0(17,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,32)
     Tmp0(18,4) = Tmp0(18,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,33)
     Tmp0(19,4) = Tmp0(19,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,34)
     Tmp0(20,4) = Tmp0(20,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,35)
     Tmp1(21,4) = Tmp1(21,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,38)
     Tmp1(22,4) = Tmp1(22,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,40)
     Tmp1(23,4) = Tmp1(23,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,41)
     Tmp1(24,4) = Tmp1(24,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,43)
     Tmp1(25,4) = Tmp1(25,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,44)
     Tmp1(26,4) = Tmp1(26,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,45)
     Tmp1(27,4) = Tmp1(27,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,47)
     Tmp1(28,4) = Tmp1(28,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,48)
     Tmp1(29,4) = Tmp1(29,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,49)
     Tmp1(30,4) = Tmp1(30,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,50)
     Tmp1(31,4) = Tmp1(31,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,52)
     Tmp1(32,4) = Tmp1(32,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,53)
     Tmp1(33,4) = Tmp1(33,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,54)
     Tmp1(34,4) = Tmp1(34,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,55)
     Tmp1(35,4) = Tmp1(35,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,56)
 ! Building for Angular momentum Jq = 2
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,6) = facX*Tmp0(iTUVP,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,7) = facX*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
     Tmp0(2,5) = Tmp0(2,5) + inv2expQ*Tmp0(1,2) 
     Tmp0(5,5) = Tmp0(5,5) + 2*inv2expQ*Tmp0(2,2) 
     Tmp0(6,5) = Tmp0(6,5) + inv2expQ*Tmp0(3,2) 
     Tmp0(7,5) = Tmp0(7,5) + inv2expQ*Tmp0(4,2) 
     Tmp0(11,5) = Tmp0(11,5) + 3*inv2expQ*Tmp0(5,2) 
     Tmp0(12,5) = Tmp0(12,5) + 2*inv2expQ*Tmp0(6,2) 
     Tmp0(13,5) = Tmp0(13,5) + 2*inv2expQ*Tmp0(7,2) 
     Tmp0(14,5) = Tmp0(14,5) + inv2expQ*Tmp0(8,2) 
     Tmp0(15,5) = Tmp0(15,5) + inv2expQ*Tmp0(9,2) 
     Tmp0(16,5) = Tmp0(16,5) + inv2expQ*Tmp0(10,2) 
     Tmp0(2,6) = Tmp0(2,6) + inv2expQ*Tmp0(1,3) 
     Tmp0(5,6) = Tmp0(5,6) + 2*inv2expQ*Tmp0(2,3) 
     Tmp0(6,6) = Tmp0(6,6) + inv2expQ*Tmp0(3,3) 
     Tmp0(7,6) = Tmp0(7,6) + inv2expQ*Tmp0(4,3) 
     Tmp0(11,6) = Tmp0(11,6) + 3*inv2expQ*Tmp0(5,3) 
     Tmp0(12,6) = Tmp0(12,6) + 2*inv2expQ*Tmp0(6,3) 
     Tmp0(13,6) = Tmp0(13,6) + 2*inv2expQ*Tmp0(7,3) 
     Tmp0(14,6) = Tmp0(14,6) + inv2expQ*Tmp0(8,3) 
     Tmp0(15,6) = Tmp0(15,6) + inv2expQ*Tmp0(9,3) 
     Tmp0(16,6) = Tmp0(16,6) + inv2expQ*Tmp0(10,3) 
     Tmp0(2,7) = Tmp0(2,7) + inv2expQ*Tmp0(1,4) 
     Tmp0(5,7) = Tmp0(5,7) + 2*inv2expQ*Tmp0(2,4) 
     Tmp0(6,7) = Tmp0(6,7) + inv2expQ*Tmp0(3,4) 
     Tmp0(7,7) = Tmp0(7,7) + inv2expQ*Tmp0(4,4) 
     Tmp0(11,7) = Tmp0(11,7) + 3*inv2expQ*Tmp0(5,4) 
     Tmp0(12,7) = Tmp0(12,7) + 2*inv2expQ*Tmp0(6,4) 
     Tmp0(13,7) = Tmp0(13,7) + 2*inv2expQ*Tmp0(7,4) 
     Tmp0(14,7) = Tmp0(14,7) + inv2expQ*Tmp0(8,4) 
     Tmp0(15,7) = Tmp0(15,7) + inv2expQ*Tmp0(9,4) 
     Tmp0(16,7) = Tmp0(16,7) + inv2expQ*Tmp0(10,4) 
     Tmp0(3,8) = Tmp0(3,8) + inv2expQ*Tmp0(1,3) 
     Tmp0(6,8) = Tmp0(6,8) + inv2expQ*Tmp0(2,3) 
     Tmp0(8,8) = Tmp0(8,8) + 2*inv2expQ*Tmp0(3,3) 
     Tmp0(9,8) = Tmp0(9,8) + inv2expQ*Tmp0(4,3) 
     Tmp0(12,8) = Tmp0(12,8) + inv2expQ*Tmp0(5,3) 
     Tmp0(14,8) = Tmp0(14,8) + 2*inv2expQ*Tmp0(6,3) 
     Tmp0(15,8) = Tmp0(15,8) + inv2expQ*Tmp0(7,3) 
     Tmp0(17,8) = Tmp0(17,8) + 3*inv2expQ*Tmp0(8,3) 
     Tmp0(18,8) = Tmp0(18,8) + 2*inv2expQ*Tmp0(9,3) 
     Tmp0(19,8) = Tmp0(19,8) + inv2expQ*Tmp0(10,3) 
     Tmp0(3,9) = Tmp0(3,9) + inv2expQ*Tmp0(1,4) 
     Tmp0(6,9) = Tmp0(6,9) + inv2expQ*Tmp0(2,4) 
     Tmp0(8,9) = Tmp0(8,9) + 2*inv2expQ*Tmp0(3,4) 
     Tmp0(9,9) = Tmp0(9,9) + inv2expQ*Tmp0(4,4) 
     Tmp0(12,9) = Tmp0(12,9) + inv2expQ*Tmp0(5,4) 
     Tmp0(14,9) = Tmp0(14,9) + 2*inv2expQ*Tmp0(6,4) 
     Tmp0(15,9) = Tmp0(15,9) + inv2expQ*Tmp0(7,4) 
     Tmp0(17,9) = Tmp0(17,9) + 3*inv2expQ*Tmp0(8,4) 
     Tmp0(18,9) = Tmp0(18,9) + 2*inv2expQ*Tmp0(9,4) 
     Tmp0(19,9) = Tmp0(19,9) + inv2expQ*Tmp0(10,4) 
     Tmp0(4,10) = Tmp0(4,10) + inv2expQ*Tmp0(1,4) 
     Tmp0(7,10) = Tmp0(7,10) + inv2expQ*Tmp0(2,4) 
     Tmp0(9,10) = Tmp0(9,10) + inv2expQ*Tmp0(3,4) 
     Tmp0(10,10) = Tmp0(10,10) + 2*inv2expQ*Tmp0(4,4) 
     Tmp0(13,10) = Tmp0(13,10) + inv2expQ*Tmp0(5,4) 
     Tmp0(15,10) = Tmp0(15,10) + inv2expQ*Tmp0(6,4) 
     Tmp0(16,10) = Tmp0(16,10) + 2*inv2expQ*Tmp0(7,4) 
     Tmp0(18,10) = Tmp0(18,10) + inv2expQ*Tmp0(8,4) 
     Tmp0(19,10) = Tmp0(19,10) + 2*inv2expQ*Tmp0(9,4) 
     Tmp0(20,10) = Tmp0(20,10) + 3*inv2expQ*Tmp0(10,4) 
     Tmp0(1,5) = Tmp0(1,5) + pinvq*Tmp0(2,2)
     Tmp0(2,5) = Tmp0(2,5) + pinvq*Tmp0(5,2)
     Tmp0(3,5) = Tmp0(3,5) + pinvq*Tmp0(6,2)
     Tmp0(4,5) = Tmp0(4,5) + pinvq*Tmp0(7,2)
     Tmp0(5,5) = Tmp0(5,5) + pinvq*Tmp0(11,2)
     Tmp0(6,5) = Tmp0(6,5) + pinvq*Tmp0(12,2)
     Tmp0(7,5) = Tmp0(7,5) + pinvq*Tmp0(13,2)
     Tmp0(8,5) = Tmp0(8,5) + pinvq*Tmp0(14,2)
     Tmp0(9,5) = Tmp0(9,5) + pinvq*Tmp0(15,2)
     Tmp0(10,5) = Tmp0(10,5) + pinvq*Tmp0(16,2)
     Tmp0(11,5) = Tmp0(11,5) + pinvq*Tmp1(21,2)
     Tmp0(12,5) = Tmp0(12,5) + pinvq*Tmp1(22,2)
     Tmp0(13,5) = Tmp0(13,5) + pinvq*Tmp1(23,2)
     Tmp0(14,5) = Tmp0(14,5) + pinvq*Tmp1(24,2)
     Tmp0(15,5) = Tmp0(15,5) + pinvq*Tmp1(25,2)
     Tmp0(16,5) = Tmp0(16,5) + pinvq*Tmp1(26,2)
     Tmp0(17,5) = Tmp0(17,5) + pinvq*Tmp1(27,2)
     Tmp0(18,5) = Tmp0(18,5) + pinvq*Tmp1(28,2)
     Tmp0(19,5) = Tmp0(19,5) + pinvq*Tmp1(29,2)
     Tmp0(20,5) = Tmp0(20,5) + pinvq*Tmp1(30,2)
     Tmp0(1,6) = Tmp0(1,6) + pinvq*Tmp0(2,3)
     Tmp0(2,6) = Tmp0(2,6) + pinvq*Tmp0(5,3)
     Tmp0(3,6) = Tmp0(3,6) + pinvq*Tmp0(6,3)
     Tmp0(4,6) = Tmp0(4,6) + pinvq*Tmp0(7,3)
     Tmp0(5,6) = Tmp0(5,6) + pinvq*Tmp0(11,3)
     Tmp0(6,6) = Tmp0(6,6) + pinvq*Tmp0(12,3)
     Tmp0(7,6) = Tmp0(7,6) + pinvq*Tmp0(13,3)
     Tmp0(8,6) = Tmp0(8,6) + pinvq*Tmp0(14,3)
     Tmp0(9,6) = Tmp0(9,6) + pinvq*Tmp0(15,3)
     Tmp0(10,6) = Tmp0(10,6) + pinvq*Tmp0(16,3)
     Tmp0(11,6) = Tmp0(11,6) + pinvq*Tmp1(21,3)
     Tmp0(12,6) = Tmp0(12,6) + pinvq*Tmp1(22,3)
     Tmp0(13,6) = Tmp0(13,6) + pinvq*Tmp1(23,3)
     Tmp0(14,6) = Tmp0(14,6) + pinvq*Tmp1(24,3)
     Tmp0(15,6) = Tmp0(15,6) + pinvq*Tmp1(25,3)
     Tmp0(16,6) = Tmp0(16,6) + pinvq*Tmp1(26,3)
     Tmp0(17,6) = Tmp0(17,6) + pinvq*Tmp1(27,3)
     Tmp0(18,6) = Tmp0(18,6) + pinvq*Tmp1(28,3)
     Tmp0(19,6) = Tmp0(19,6) + pinvq*Tmp1(29,3)
     Tmp0(20,6) = Tmp0(20,6) + pinvq*Tmp1(30,3)
     Tmp0(1,7) = Tmp0(1,7) + pinvq*Tmp0(2,4)
     Tmp0(2,7) = Tmp0(2,7) + pinvq*Tmp0(5,4)
     Tmp0(3,7) = Tmp0(3,7) + pinvq*Tmp0(6,4)
     Tmp0(4,7) = Tmp0(4,7) + pinvq*Tmp0(7,4)
     Tmp0(5,7) = Tmp0(5,7) + pinvq*Tmp0(11,4)
     Tmp0(6,7) = Tmp0(6,7) + pinvq*Tmp0(12,4)
     Tmp0(7,7) = Tmp0(7,7) + pinvq*Tmp0(13,4)
     Tmp0(8,7) = Tmp0(8,7) + pinvq*Tmp0(14,4)
     Tmp0(9,7) = Tmp0(9,7) + pinvq*Tmp0(15,4)
     Tmp0(10,7) = Tmp0(10,7) + pinvq*Tmp0(16,4)
     Tmp0(11,7) = Tmp0(11,7) + pinvq*Tmp1(21,4)
     Tmp0(12,7) = Tmp0(12,7) + pinvq*Tmp1(22,4)
     Tmp0(13,7) = Tmp0(13,7) + pinvq*Tmp1(23,4)
     Tmp0(14,7) = Tmp0(14,7) + pinvq*Tmp1(24,4)
     Tmp0(15,7) = Tmp0(15,7) + pinvq*Tmp1(25,4)
     Tmp0(16,7) = Tmp0(16,7) + pinvq*Tmp1(26,4)
     Tmp0(17,7) = Tmp0(17,7) + pinvq*Tmp1(27,4)
     Tmp0(18,7) = Tmp0(18,7) + pinvq*Tmp1(28,4)
     Tmp0(19,7) = Tmp0(19,7) + pinvq*Tmp1(29,4)
     Tmp0(20,7) = Tmp0(20,7) + pinvq*Tmp1(30,4)
     Tmp0(1,8) = Tmp0(1,8) + pinvq*Tmp0(3,3)
     Tmp0(2,8) = Tmp0(2,8) + pinvq*Tmp0(6,3)
     Tmp0(3,8) = Tmp0(3,8) + pinvq*Tmp0(8,3)
     Tmp0(4,8) = Tmp0(4,8) + pinvq*Tmp0(9,3)
     Tmp0(5,8) = Tmp0(5,8) + pinvq*Tmp0(12,3)
     Tmp0(6,8) = Tmp0(6,8) + pinvq*Tmp0(14,3)
     Tmp0(7,8) = Tmp0(7,8) + pinvq*Tmp0(15,3)
     Tmp0(8,8) = Tmp0(8,8) + pinvq*Tmp0(17,3)
     Tmp0(9,8) = Tmp0(9,8) + pinvq*Tmp0(18,3)
     Tmp0(10,8) = Tmp0(10,8) + pinvq*Tmp0(19,3)
     Tmp0(11,8) = Tmp0(11,8) + pinvq*Tmp1(22,3)
     Tmp0(12,8) = Tmp0(12,8) + pinvq*Tmp1(24,3)
     Tmp0(13,8) = Tmp0(13,8) + pinvq*Tmp1(25,3)
     Tmp0(14,8) = Tmp0(14,8) + pinvq*Tmp1(27,3)
     Tmp0(15,8) = Tmp0(15,8) + pinvq*Tmp1(28,3)
     Tmp0(16,8) = Tmp0(16,8) + pinvq*Tmp1(29,3)
     Tmp0(17,8) = Tmp0(17,8) + pinvq*Tmp1(31,3)
     Tmp0(18,8) = Tmp0(18,8) + pinvq*Tmp1(32,3)
     Tmp0(19,8) = Tmp0(19,8) + pinvq*Tmp1(33,3)
     Tmp0(20,8) = Tmp0(20,8) + pinvq*Tmp1(34,3)
     Tmp0(1,9) = Tmp0(1,9) + pinvq*Tmp0(3,4)
     Tmp0(2,9) = Tmp0(2,9) + pinvq*Tmp0(6,4)
     Tmp0(3,9) = Tmp0(3,9) + pinvq*Tmp0(8,4)
     Tmp0(4,9) = Tmp0(4,9) + pinvq*Tmp0(9,4)
     Tmp0(5,9) = Tmp0(5,9) + pinvq*Tmp0(12,4)
     Tmp0(6,9) = Tmp0(6,9) + pinvq*Tmp0(14,4)
     Tmp0(7,9) = Tmp0(7,9) + pinvq*Tmp0(15,4)
     Tmp0(8,9) = Tmp0(8,9) + pinvq*Tmp0(17,4)
     Tmp0(9,9) = Tmp0(9,9) + pinvq*Tmp0(18,4)
     Tmp0(10,9) = Tmp0(10,9) + pinvq*Tmp0(19,4)
     Tmp0(11,9) = Tmp0(11,9) + pinvq*Tmp1(22,4)
     Tmp0(12,9) = Tmp0(12,9) + pinvq*Tmp1(24,4)
     Tmp0(13,9) = Tmp0(13,9) + pinvq*Tmp1(25,4)
     Tmp0(14,9) = Tmp0(14,9) + pinvq*Tmp1(27,4)
     Tmp0(15,9) = Tmp0(15,9) + pinvq*Tmp1(28,4)
     Tmp0(16,9) = Tmp0(16,9) + pinvq*Tmp1(29,4)
     Tmp0(17,9) = Tmp0(17,9) + pinvq*Tmp1(31,4)
     Tmp0(18,9) = Tmp0(18,9) + pinvq*Tmp1(32,4)
     Tmp0(19,9) = Tmp0(19,9) + pinvq*Tmp1(33,4)
     Tmp0(20,9) = Tmp0(20,9) + pinvq*Tmp1(34,4)
     Tmp0(1,10) = Tmp0(1,10) + pinvq*Tmp0(4,4)
     Tmp0(2,10) = Tmp0(2,10) + pinvq*Tmp0(7,4)
     Tmp0(3,10) = Tmp0(3,10) + pinvq*Tmp0(9,4)
     Tmp0(4,10) = Tmp0(4,10) + pinvq*Tmp0(10,4)
     Tmp0(5,10) = Tmp0(5,10) + pinvq*Tmp0(13,4)
     Tmp0(6,10) = Tmp0(6,10) + pinvq*Tmp0(15,4)
     Tmp0(7,10) = Tmp0(7,10) + pinvq*Tmp0(16,4)
     Tmp0(8,10) = Tmp0(8,10) + pinvq*Tmp0(18,4)
     Tmp0(9,10) = Tmp0(9,10) + pinvq*Tmp0(19,4)
     Tmp0(10,10) = Tmp0(10,10) + pinvq*Tmp0(20,4)
     Tmp0(11,10) = Tmp0(11,10) + pinvq*Tmp1(23,4)
     Tmp0(12,10) = Tmp0(12,10) + pinvq*Tmp1(25,4)
     Tmp0(13,10) = Tmp0(13,10) + pinvq*Tmp1(26,4)
     Tmp0(14,10) = Tmp0(14,10) + pinvq*Tmp1(28,4)
     Tmp0(15,10) = Tmp0(15,10) + pinvq*Tmp1(29,4)
     Tmp0(16,10) = Tmp0(16,10) + pinvq*Tmp1(30,4)
     Tmp0(17,10) = Tmp0(17,10) + pinvq*Tmp1(32,4)
     Tmp0(18,10) = Tmp0(18,10) + pinvq*Tmp1(33,4)
     Tmp0(19,10) = Tmp0(19,10) + pinvq*Tmp1(34,4)
     Tmp0(20,10) = Tmp0(20,10) + pinvq*Tmp1(35,4)
!$ACC LOOP SEQ
     DO iTUVQ=1, 10
!$ACC LOOP SEQ
      DO iTUVP=1, 20
        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
 end subroutine TransferRecurrenceGPUP3Q2AtoDSegP
END MODULE AGC_GPU_OBS_TRMODAtoDSegP1
