MODULE AGC_GPU_OBS_TRMODAtoCSeg1Prim2
 use IchorPrecisionMod
  
 CONTAINS
 subroutine TransferRecurrenceGPUP4Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  use AGC_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(1),Dexp(1)
  real(realk),intent(in) :: Aux(nPasses,   84)
  real(realk),intent(inout) :: Aux2(nPasses,   35,   10)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 35, 10)
  real(realk) :: Tmp1( 36: 56,  2:  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        TUVindexX1_56,TUVindexX2_56,TUVindexX3_56, &
!$ACC        IfacX1_35,IfacX2_35,IfacX3_35, &
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Dexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPasses
   iPassP = iP
   iPrimP=1
   iPrimQ=1
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = Pdistance12(1,iAtomA,iAtomB)
   Yab = Pdistance12(2,iAtomA,iAtomB)
   Zab = Pdistance12(3,iAtomA,iAtomB)
    expP = Pexp(iPrimP)
    expBX = Bexp(1)*Xab
    expBY = Bexp(1)*Yab
    expBZ = Bexp(1)*Zab
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     facX = -(expBX+Dexp(1)*Xcd)*invexpQ
     facY = -(expBY+Dexp(1)*Ycd)*invexpQ
     facZ = -(expBZ+Dexp(1)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
!$ACC LOOP SEQ
     DO iTUVP=1, 35
      Tmp0(iTUVP,1) = Aux(IP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp1(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp1(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp1(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1_35(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2_35(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3_35(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1_35(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2_35(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3_35(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jq = 2
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(IP,iTUVP)
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
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(IP,iTUVP)
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
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP4Q2AtoCSeg1Prim
 subroutine TransferRecurrenceGPUP4Q3AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  use AGC_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(1),Dexp(1)
  real(realk),intent(in) :: Aux(nPasses,  120)
  real(realk),intent(inout) :: Aux2(nPasses,   35,   20)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 35, 20)
  real(realk) :: Tmp1( 36: 84,  2:  4)
  real(realk) :: Tmp2( 36: 56,  5: 10)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         Tmp2,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        TUVindexX1_84,TUVindexX2_84,TUVindexX3_84, &
!$ACC        IfacX1_56,IfacX2_56,IfacX3_56, &
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Dexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPasses
   iPassP = iP
   iPrimP=1
   iPrimQ=1
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = Pdistance12(1,iAtomA,iAtomB)
   Yab = Pdistance12(2,iAtomA,iAtomB)
   Zab = Pdistance12(3,iAtomA,iAtomB)
    expP = Pexp(iPrimP)
    expBX = Bexp(1)*Xab
    expBY = Bexp(1)*Yab
    expBZ = Bexp(1)*Zab
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     facX = -(expBX+Dexp(1)*Xcd)*invexpQ
     facY = -(expBY+Dexp(1)*Ycd)*invexpQ
     facZ = -(expBZ+Dexp(1)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
!$ACC LOOP SEQ
     DO iTUVP=1, 35
      Tmp0(iTUVP,1) = Aux(IP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp1(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp1(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp1(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1_56(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2_56(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3_56(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1_56(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2_56(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3_56(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jq = 2
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(IP,iTUVP)
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
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp2(iTUVP,5) = facX*Tmp1(iTUVP,2)+ inv2expQ*Aux(IP,iTUVP)
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
      Tmp2(iTUVP,8) = facY*Tmp1(iTUVP,3)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp2(iTUVP,9) = facY*Tmp1(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp2(iTUVP,10) = facZ*Tmp1(iTUVP,4)+ inv2expQ*Aux(IP,iTUVP)
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
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP4Q3AtoCSeg1Prim
 subroutine TransferRecurrenceGPUP4Q4AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  use AGC_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(1),Dexp(1)
  real(realk),intent(in) :: Aux(nPasses,  165)
  real(realk),intent(inout) :: Aux2(nPasses,   35,   35)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 35, 35)
  real(realk) :: Tmp1( 36:120,  2:  4)
  real(realk) :: Tmp2( 36: 84,  5: 10)
  real(realk) :: Tmp3( 36: 56, 11: 20)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         Tmp2,&
!$ACC         Tmp3,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        TUVindexX1_120,TUVindexX2_120,TUVindexX3_120, &
!$ACC        IfacX1_84,IfacX2_84,IfacX3_84, &
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Dexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPasses
   iPassP = iP
   iPrimP=1
   iPrimQ=1
   Xcd = Qdistance12(1)
   Ycd = Qdistance12(2)
   Zcd = Qdistance12(3)
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Xab = Pdistance12(1,iAtomA,iAtomB)
   Yab = Pdistance12(2,iAtomA,iAtomB)
   Zab = Pdistance12(3,iAtomA,iAtomB)
    expP = Pexp(iPrimP)
    expBX = Bexp(1)*Xab
    expBY = Bexp(1)*Yab
    expBZ = Bexp(1)*Zab
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     facX = -(expBX+Dexp(1)*Xcd)*invexpQ
     facY = -(expBY+Dexp(1)*Ycd)*invexpQ
     facZ = -(expBZ+Dexp(1)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
!$ACC LOOP SEQ
     DO iTUVP=1, 35
      Tmp0(iTUVP,1) = Aux(IP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36,120
      Tmp1(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36,120
      Tmp1(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36,120
      Tmp1(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1_84(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2_84(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3_84(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,84
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1_84(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,84
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2_84(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,84
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3_84(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,120
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,120
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,120
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jq = 2
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(IP,iTUVP)
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
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp2(iTUVP,5) = facX*Tmp1(iTUVP,2)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp2(iTUVP,6) = facX*Tmp1(iTUVP,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp2(iTUVP,7) = facX*Tmp1(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp2(iTUVP,8) = facY*Tmp1(iTUVP,3)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp2(iTUVP,9) = facY*Tmp1(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 84
      Tmp2(iTUVP,10) = facZ*Tmp1(iTUVP,4)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp0(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX3_120(iTUVP)
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
     do iTUVP =  36, 56
      Tmp3(iTUVP,11) = facX*Tmp2(iTUVP,5)+2*inv2expQ*Tmp1(iTUVP,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp3(iTUVP,12) = facY*Tmp2(iTUVP,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp3(iTUVP,13) = facZ*Tmp2(iTUVP,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp3(iTUVP,14) = facX*Tmp2(iTUVP,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp3(iTUVP,15) = facX*Tmp2(iTUVP,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp3(iTUVP,16) = facX*Tmp2(iTUVP,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp3(iTUVP,17) = facY*Tmp2(iTUVP,8)+2*inv2expQ*Tmp1(iTUVP,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp3(iTUVP,18) = facZ*Tmp2(iTUVP,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp3(iTUVP,19) = facY*Tmp2(iTUVP,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  36, 56
      Tmp3(iTUVP,20) = facZ*Tmp2(iTUVP,10)+2*inv2expQ*Tmp1(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,9) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp3(iTUVP,11) = Tmp3(iTUVP,11) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp3(iTUVP,12) = Tmp3(iTUVP,12) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp3(iTUVP,13) = Tmp3(iTUVP,13) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp3(iTUVP,14) = Tmp3(iTUVP,14) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp3(iTUVP,15) = Tmp3(iTUVP,15) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,9) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp3(iTUVP,16) = Tmp3(iTUVP,16) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp3(iTUVP,17) = Tmp3(iTUVP,17) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp3(iTUVP,18) = Tmp3(iTUVP,18) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp3(iTUVP,19) = Tmp3(iTUVP,19) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp3(iTUVP,20) = Tmp3(iTUVP,20) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp3(iTUVP,11) = Tmp3(iTUVP,11) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp3(iTUVP,12) = Tmp3(iTUVP,12) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp3(iTUVP,13) = Tmp3(iTUVP,13) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp3(iTUVP,14) = Tmp3(iTUVP,14) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp0(iTUVplus1,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp2(iTUVplus1,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp3(iTUVP,15) = Tmp3(iTUVP,15) + pinvq*Tmp2(iTUVplus1,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp3(iTUVP,16) = Tmp3(iTUVP,16) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp3(iTUVP,17) = Tmp3(iTUVP,17) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp3(iTUVP,18) = Tmp3(iTUVP,18) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp3(iTUVP,19) = Tmp3(iTUVP,19) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp3(iTUVP,20) = Tmp3(iTUVP,20) + pinvq*Tmp2(iTUVplus1,10)
     enddo
 ! Building for Angular momentum Jq = 4
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,21) = facX*Tmp0(iTUVP,11)+3*inv2expQ*Tmp0(iTUVP,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,22) = facY*Tmp0(iTUVP,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,23) = facZ*Tmp0(iTUVP,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,24) = facX*Tmp0(iTUVP,14)+ inv2expQ*Tmp0(iTUVP,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,25) = facY*Tmp0(iTUVP,13)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,26) = facX*Tmp0(iTUVP,16)+ inv2expQ*Tmp0(iTUVP,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,27) = facX*Tmp0(iTUVP,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,28) = facX*Tmp0(iTUVP,18)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,29) = facX*Tmp0(iTUVP,19)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,30) = facX*Tmp0(iTUVP,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,31) = facY*Tmp0(iTUVP,17)+3*inv2expQ*Tmp0(iTUVP,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,32) = facZ*Tmp0(iTUVP,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,33) = facY*Tmp0(iTUVP,19)+ inv2expQ*Tmp0(iTUVP,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,34) = facY*Tmp0(iTUVP,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 35
      Tmp0(iTUVP,35) = facZ*Tmp0(iTUVP,20)+3*inv2expQ*Tmp0(iTUVP,10)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,21) = Tmp0(iTUVP,21) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,11) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,22) = Tmp0(iTUVP,22) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,11) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,23) = Tmp0(iTUVP,23) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,11) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,24) = Tmp0(iTUVP,24) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,14) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,25) = Tmp0(iTUVP,25) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,13) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,26) = Tmp0(iTUVP,26) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,16) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,27) = Tmp0(iTUVP,27) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,17) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,28) = Tmp0(iTUVP,28) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,18) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,29) = Tmp0(iTUVP,29) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,19) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,30) = Tmp0(iTUVP,30) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,20) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,31) = Tmp0(iTUVP,31) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,17) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,32) = Tmp0(iTUVP,32) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,17) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,33) = Tmp0(iTUVP,33) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,19) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,34) = Tmp0(iTUVP,34) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,20) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,35) = Tmp0(iTUVP,35) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,20) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,21) = Tmp0(iTUVP,21) + pinvq*Tmp0(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,21) = Tmp0(iTUVP,21) + pinvq*Tmp3(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,22) = Tmp0(iTUVP,22) + pinvq*Tmp0(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,22) = Tmp0(iTUVP,22) + pinvq*Tmp3(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,23) = Tmp0(iTUVP,23) + pinvq*Tmp0(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,23) = Tmp0(iTUVP,23) + pinvq*Tmp3(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,24) = Tmp0(iTUVP,24) + pinvq*Tmp0(iTUVplus1,14)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,24) = Tmp0(iTUVP,24) + pinvq*Tmp3(iTUVplus1,14)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,25) = Tmp0(iTUVP,25) + pinvq*Tmp0(iTUVplus1,13)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,25) = Tmp0(iTUVP,25) + pinvq*Tmp3(iTUVplus1,13)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,26) = Tmp0(iTUVP,26) + pinvq*Tmp0(iTUVplus1,16)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,26) = Tmp0(iTUVP,26) + pinvq*Tmp3(iTUVplus1,16)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,27) = Tmp0(iTUVP,27) + pinvq*Tmp0(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,27) = Tmp0(iTUVP,27) + pinvq*Tmp3(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,28) = Tmp0(iTUVP,28) + pinvq*Tmp0(iTUVplus1,18)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,28) = Tmp0(iTUVP,28) + pinvq*Tmp3(iTUVplus1,18)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,29) = Tmp0(iTUVP,29) + pinvq*Tmp0(iTUVplus1,19)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,29) = Tmp0(iTUVP,29) + pinvq*Tmp3(iTUVplus1,19)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,30) = Tmp0(iTUVP,30) + pinvq*Tmp0(iTUVplus1,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,30) = Tmp0(iTUVP,30) + pinvq*Tmp3(iTUVplus1,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,31) = Tmp0(iTUVP,31) + pinvq*Tmp0(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,31) = Tmp0(iTUVP,31) + pinvq*Tmp3(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,32) = Tmp0(iTUVP,32) + pinvq*Tmp0(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,32) = Tmp0(iTUVP,32) + pinvq*Tmp3(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,33) = Tmp0(iTUVP,33) + pinvq*Tmp0(iTUVplus1,19)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,33) = Tmp0(iTUVP,33) + pinvq*Tmp3(iTUVplus1,19)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,34) = Tmp0(iTUVP,34) + pinvq*Tmp0(iTUVplus1,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,34) = Tmp0(iTUVP,34) + pinvq*Tmp3(iTUVplus1,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,35) = Tmp0(iTUVP,35) + pinvq*Tmp0(iTUVplus1,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,35) = Tmp0(iTUVP,35) + pinvq*Tmp3(iTUVplus1,20)
     enddo
!$ACC LOOP SEQ
     DO iTUVQ=1, 35
!$ACC LOOP SEQ
      DO iTUVP=1, 35
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP4Q4AtoCSeg1Prim
end module
