MODULE AGC_GPU_OBS_TRMODAtoDSeg2
 use IchorPrecisionMod
  
 CONTAINS
 subroutine TransferRecurrenceGPUP4Q2AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  use AGC_GPU_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   84)
  real(realk),intent(inout) :: Aux2(nPasses,   35,   10)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 35, 10)
  real(realk) :: Tmp1( 36: 56,  2:  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         iPrimB,iPrimC,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        TUVindexX1_56,TUVindexX2_56,TUVindexX3_56, &
!$ACC        IfacX1_35,IfacX2_35,IfacX3_35, &
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPasses
   DO iTUVQ=1, 10
    DO iTUVP=1, 35
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
   DO iPrimQP=1,nPrimQ*nPrimP
    iPrimQ = iPrimQP - ((iPrimQP-1)/nPrimQ)*nPrimQ
    iPrimP = (iPrimQP-1)/nPrimQ + 1
    iPassP = iP
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
     DO iTUVP=1, 35
      Tmp0(iTUVP,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp1(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp1(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp1(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1_35(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2_35(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3_35(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1_35(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2_35(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3_35(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jq = 2
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,6) = facX*Tmp0(iTUVP,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,7) = facX*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + IfacX2_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + IfacX2_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + IfacX3_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp0(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     DO iTUVQ=1, 10
!$ACC LOOP SEQ
      DO iTUVP=1, 35
        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP4Q2AtoDSeg
 subroutine TransferRecurrenceGPUP4Q3AtoDSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  use AGC_GPU_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,  120)
  real(realk),intent(inout) :: Aux2(nPasses,   35,   20)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 35, 20)
  real(realk) :: Tmp1( 36: 84,  2:  4)
  real(realk) :: Tmp2( 36: 56,  5: 10)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimC
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         iPrimB,iPrimC,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         Tmp2,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        TUVindexX1_84,TUVindexX2_84,TUVindexX3_84, &
!$ACC        IfacX1_56,IfacX2_56,IfacX3_56, &
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPasses
   DO iTUVQ=1, 20
    DO iTUVP=1, 35
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
   DO iPrimQP=1,nPrimQ*nPrimP
    iPrimQ = iPrimQP - ((iPrimQP-1)/nPrimQ)*nPrimQ
    iPrimP = (iPrimQP-1)/nPrimQ + 1
    iPassP = iP
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
     DO iTUVP=1, 35
      Tmp0(iTUVP,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp1(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp1(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp1(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1_56(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2_56(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3_56(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1_56(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2_56(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3_56(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jq = 2
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,6) = facX*Tmp0(iTUVP,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,7) = facX*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp2(iTUVP,5) = facX*Tmp1(iTUVP,2)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp2(iTUVP,6) = facX*Tmp1(iTUVP,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp2(iTUVP,7) = facX*Tmp1(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp2(iTUVP,8) = facY*Tmp1(iTUVP,3)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp2(iTUVP,9) = facY*Tmp1(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp2(iTUVP,10) = facZ*Tmp1(iTUVP,4)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + IfacX3_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + IfacX3_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp0(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
 ! Building for Angular momentum Jq = 3
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,11) = facX*Tmp0(iTUVP,5)+2*inv2expQ*Tmp0(iTUVP,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,12) = facY*Tmp0(iTUVP,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,13) = facZ*Tmp0(iTUVP,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,14) = facX*Tmp0(iTUVP,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,15) = facX*Tmp0(iTUVP,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,16) = facX*Tmp0(iTUVP,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,17) = facY*Tmp0(iTUVP,8)+2*inv2expQ*Tmp0(iTUVP,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,18) = facZ*Tmp0(iTUVP,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,19) = facY*Tmp0(iTUVP,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,20) = facZ*Tmp0(iTUVP,10)+2*inv2expQ*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + IfacX3_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,9) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + IfacX3_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + IfacX3_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp0(iTUVplus1,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp2(iTUVplus1,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     DO iTUVQ=1, 20
!$ACC LOOP SEQ
      DO iTUVP=1, 35
        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP4Q3AtoDSeg
end module
