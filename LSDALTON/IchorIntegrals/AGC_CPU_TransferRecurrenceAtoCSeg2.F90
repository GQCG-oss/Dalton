MODULE AGC_CPU_OBS_TRMODAtoCSeg2
 use IchorPrecisionMod
  
 CONTAINS
 subroutine TransferRecurrenceCPUP4Q2AtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  use AGC_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Dexp(nPrimD)
  real(realk),intent(in) :: Aux(   84,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   35,   10,nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 35, 10)
  real(realk) :: Tmp1( 36: 56,  2:  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimD
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$OMP DO COLLAPSE(3) PRIVATE(iP,iTUVP,iTUVQ)
  DO iP = 1,nPasses
   DO iTUVQ=1, 10
    DO iTUVP=1, 35
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimB,iPrimD,&
!$OMP         Tmp0,&
!$OMP         Tmp1,&
!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1)
  DO iP = 1,nPasses
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
     iPrimD = (iPrimQ-1)/nPrimC+1                
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpQ
     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpQ
     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpQ
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
     do iTUVP =  36, 56
      Tmp1(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36, 56
      Tmp1(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36, 56
      Tmp1(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3_35(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
 ! Building for Angular momentum Jq = 2
     do iTUVP = 1, 35
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,6) = facX*Tmp0(iTUVP,3)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,7) = facX*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_56(ituvpminus1)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + IfacX1_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + IfacX2_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_56(ituvpminus1)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + IfacX2_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_56(ituvpminus1)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + IfacX3_35(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp0(iTUVplus1,2)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp0(iTUVplus1,3)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp0(iTUVplus1,3)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_56(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_56(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     DO iTUVQ=1, 10
      DO iTUVP=1, 35
        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP4Q2AtoCSeg
 subroutine TransferRecurrenceCPUP4Q3AtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  use AGC_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Dexp(nPrimD)
  real(realk),intent(in) :: Aux(  120,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   35,   20,nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 35, 20)
  real(realk) :: Tmp1( 36: 84,  2:  4)
  real(realk) :: Tmp2( 36: 56,  5: 10)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimD
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$OMP DO COLLAPSE(3) PRIVATE(iP,iTUVP,iTUVQ)
  DO iP = 1,nPasses
   DO iTUVQ=1, 20
    DO iTUVP=1, 35
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimB,iPrimD,&
!$OMP         Tmp0,&
!$OMP         Tmp1,&
!$OMP         Tmp2,&
!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1)
  DO iP = 1,nPasses
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
     iPrimD = (iPrimQ-1)/nPrimC+1                
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpQ
     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpQ
     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpQ
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
     do iTUVP =  36, 84
      Tmp1(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36, 84
      Tmp1(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36, 84
      Tmp1(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1_56(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2_56(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3_56(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1_56(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2_56(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3_56(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
 ! Building for Angular momentum Jq = 2
     do iTUVP = 1, 35
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,6) = facX*Tmp0(iTUVP,3)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,7) = facX*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36, 56
      Tmp2(iTUVP,5) = facX*Tmp1(iTUVP,2)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36, 56
      Tmp2(iTUVP,6) = facX*Tmp1(iTUVP,3)
     enddo
     do iTUVP =  36, 56
      Tmp2(iTUVP,7) = facX*Tmp1(iTUVP,4)
     enddo
     do iTUVP =  36, 56
      Tmp2(iTUVP,8) = facY*Tmp1(iTUVP,3)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36, 56
      Tmp2(iTUVP,9) = facY*Tmp1(iTUVP,4)
     enddo
     do iTUVP =  36, 56
      Tmp2(iTUVP,10) = facZ*Tmp1(iTUVP,4)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + IfacX3_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + IfacX3_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp0(iTUVplus1,2)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp0(iTUVplus1,3)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp0(iTUVplus1,3)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
 ! Building for Angular momentum Jq = 3
     do iTUVP = 1, 35
      Tmp0(iTUVP,11) = facX*Tmp0(iTUVP,5)+2*inv2expQ*Tmp0(iTUVP,2)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,12) = facY*Tmp0(iTUVP,5)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,13) = facZ*Tmp0(iTUVP,5)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,14) = facX*Tmp0(iTUVP,8)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,15) = facX*Tmp0(iTUVP,9)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,16) = facX*Tmp0(iTUVP,10)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,17) = facY*Tmp0(iTUVP,8)+2*inv2expQ*Tmp0(iTUVP,3)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,18) = facZ*Tmp0(iTUVP,8)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,19) = facY*Tmp0(iTUVP,10)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,20) = facZ*Tmp0(iTUVP,10)+2*inv2expQ*Tmp0(iTUVP,4)
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + IfacX3_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,9) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_84(ituvpminus1)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + IfacX1_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + IfacX3_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_84(ituvpminus1)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + IfacX2_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_84(ituvpminus1)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + IfacX3_56(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp0(iTUVplus1,5)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp0(iTUVplus1,5)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp0(iTUVplus1,5)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp0(iTUVplus1,8)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp0(iTUVplus1,9)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp2(iTUVplus1,9)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp0(iTUVplus1,10)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp2(iTUVplus1,10)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp0(iTUVplus1,8)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp0(iTUVplus1,8)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp0(iTUVplus1,10)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp2(iTUVplus1,10)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp0(iTUVplus1,10)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp2(iTUVplus1,10)
     enddo
     DO iTUVQ=1, 20
      DO iTUVP=1, 35
        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP4Q3AtoCSeg
 subroutine TransferRecurrenceCPUP4Q4AtoCSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  use AGC_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Dexp(nPrimD)
  real(realk),intent(in) :: Aux(  165,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   35,   35,nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 35, 35)
  real(realk) :: Tmp1( 36:120,  2:  4)
  real(realk) :: Tmp2( 36: 84,  5: 10)
  real(realk) :: Tmp3( 36: 56, 11: 20)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  integer :: iPrimB,iPrimD
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$OMP DO COLLAPSE(3) PRIVATE(iP,iTUVP,iTUVQ)
  DO iP = 1,nPasses
   DO iTUVQ=1, 35
    DO iTUVP=1, 35
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimB,iPrimD,&
!$OMP         Tmp0,&
!$OMP         Tmp1,&
!$OMP         Tmp2,&
!$OMP         Tmp3,&
!$OMP         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1)
  DO iP = 1,nPasses
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
     iPrimD = (iPrimQ-1)/nPrimC+1                
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpQ
     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpQ
     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpQ
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
     do iTUVP =  36,120
      Tmp1(iTUVP,2) = facX*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36,120
      Tmp1(iTUVP,3) = facY*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36,120
      Tmp1(iTUVP,4) = facZ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1_84(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2_84(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3_84(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 21,84
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1_84(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 21,84
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2_84(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvpminus1 = 21,84
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3_84(ituvpminus1)*inv2expQ*Aux(ituvpminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 36,120
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 36,120
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 36,120
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
 ! Building for Angular momentum Jq = 2
     do iTUVP = 1, 35
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,6) = facX*Tmp0(iTUVP,3)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,7) = facX*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36, 84
      Tmp2(iTUVP,5) = facX*Tmp1(iTUVP,2)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36, 84
      Tmp2(iTUVP,6) = facX*Tmp1(iTUVP,3)
     enddo
     do iTUVP =  36, 84
      Tmp2(iTUVP,7) = facX*Tmp1(iTUVP,4)
     enddo
     do iTUVP =  36, 84
      Tmp2(iTUVP,8) = facY*Tmp1(iTUVP,3)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVP =  36, 84
      Tmp2(iTUVP,9) = facY*Tmp1(iTUVP,4)
     enddo
     do iTUVP =  36, 84
      Tmp2(iTUVP,10) = facZ*Tmp1(iTUVP,4)+ inv2expQ*Aux(iTUVP,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,2) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,3) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,4) 
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp0(iTUVplus1,2)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp0(iTUVplus1,3)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp0(iTUVplus1,3)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp0(iTUVplus1,4)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
 ! Building for Angular momentum Jq = 3
     do iTUVP = 1, 35
      Tmp0(iTUVP,11) = facX*Tmp0(iTUVP,5)+2*inv2expQ*Tmp0(iTUVP,2)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,12) = facY*Tmp0(iTUVP,5)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,13) = facZ*Tmp0(iTUVP,5)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,14) = facX*Tmp0(iTUVP,8)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,15) = facX*Tmp0(iTUVP,9)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,16) = facX*Tmp0(iTUVP,10)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,17) = facY*Tmp0(iTUVP,8)+2*inv2expQ*Tmp0(iTUVP,3)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,18) = facZ*Tmp0(iTUVP,8)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,19) = facY*Tmp0(iTUVP,10)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,20) = facZ*Tmp0(iTUVP,10)+2*inv2expQ*Tmp0(iTUVP,4)
     enddo
     do iTUVP =  36, 56
      Tmp3(iTUVP,11) = facX*Tmp2(iTUVP,5)+2*inv2expQ*Tmp1(iTUVP,2)
     enddo
     do iTUVP =  36, 56
      Tmp3(iTUVP,12) = facY*Tmp2(iTUVP,5)
     enddo
     do iTUVP =  36, 56
      Tmp3(iTUVP,13) = facZ*Tmp2(iTUVP,5)
     enddo
     do iTUVP =  36, 56
      Tmp3(iTUVP,14) = facX*Tmp2(iTUVP,8)
     enddo
     do iTUVP =  36, 56
      Tmp3(iTUVP,15) = facX*Tmp2(iTUVP,9)
     enddo
     do iTUVP =  36, 56
      Tmp3(iTUVP,16) = facX*Tmp2(iTUVP,10)
     enddo
     do iTUVP =  36, 56
      Tmp3(iTUVP,17) = facY*Tmp2(iTUVP,8)+2*inv2expQ*Tmp1(iTUVP,3)
     enddo
     do iTUVP =  36, 56
      Tmp3(iTUVP,18) = facZ*Tmp2(iTUVP,8)
     enddo
     do iTUVP =  36, 56
      Tmp3(iTUVP,19) = facY*Tmp2(iTUVP,10)
     enddo
     do iTUVP =  36, 56
      Tmp3(iTUVP,20) = facZ*Tmp2(iTUVP,10)+2*inv2expQ*Tmp1(iTUVP,4)
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,9) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp3(iTUVP,11) = Tmp3(iTUVP,11) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp3(iTUVP,12) = Tmp3(iTUVP,12) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp3(iTUVP,13) = Tmp3(iTUVP,13) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp3(iTUVP,14) = Tmp3(iTUVP,14) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp3(iTUVP,15) = Tmp3(iTUVP,15) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,9) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp3(iTUVP,16) = Tmp3(iTUVP,16) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp3(iTUVP,17) = Tmp3(iTUVP,17) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp3(iTUVP,18) = Tmp3(iTUVP,18) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp3(iTUVP,19) = Tmp3(iTUVP,19) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp3(iTUVP,20) = Tmp3(iTUVP,20) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp0(iTUVplus1,5)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp3(iTUVP,11) = Tmp3(iTUVP,11) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp0(iTUVplus1,5)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp3(iTUVP,12) = Tmp3(iTUVP,12) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp0(iTUVplus1,5)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp3(iTUVP,13) = Tmp3(iTUVP,13) + pinvq*Tmp2(iTUVplus1,5)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp0(iTUVplus1,8)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp3(iTUVP,14) = Tmp3(iTUVP,14) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp0(iTUVplus1,9)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp2(iTUVplus1,9)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp3(iTUVP,15) = Tmp3(iTUVP,15) + pinvq*Tmp2(iTUVplus1,9)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp0(iTUVplus1,10)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp2(iTUVplus1,10)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp3(iTUVP,16) = Tmp3(iTUVP,16) + pinvq*Tmp2(iTUVplus1,10)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp0(iTUVplus1,8)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp3(iTUVP,17) = Tmp3(iTUVP,17) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp0(iTUVplus1,8)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp3(iTUVP,18) = Tmp3(iTUVP,18) + pinvq*Tmp2(iTUVplus1,8)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp0(iTUVplus1,10)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp2(iTUVplus1,10)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp3(iTUVP,19) = Tmp3(iTUVP,19) + pinvq*Tmp2(iTUVplus1,10)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp0(iTUVplus1,10)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp2(iTUVplus1,10)
     enddo
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp3(iTUVP,20) = Tmp3(iTUVP,20) + pinvq*Tmp2(iTUVplus1,10)
     enddo
 ! Building for Angular momentum Jq = 4
     do iTUVP = 1, 35
      Tmp0(iTUVP,21) = facX*Tmp0(iTUVP,11)+3*inv2expQ*Tmp0(iTUVP,5)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,22) = facY*Tmp0(iTUVP,11)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,23) = facZ*Tmp0(iTUVP,11)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,24) = facX*Tmp0(iTUVP,14)+ inv2expQ*Tmp0(iTUVP,8)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,25) = facY*Tmp0(iTUVP,13)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,26) = facX*Tmp0(iTUVP,16)+ inv2expQ*Tmp0(iTUVP,10)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,27) = facX*Tmp0(iTUVP,17)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,28) = facX*Tmp0(iTUVP,18)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,29) = facX*Tmp0(iTUVP,19)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,30) = facX*Tmp0(iTUVP,20)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,31) = facY*Tmp0(iTUVP,17)+3*inv2expQ*Tmp0(iTUVP,8)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,32) = facZ*Tmp0(iTUVP,17)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,33) = facY*Tmp0(iTUVP,19)+ inv2expQ*Tmp0(iTUVP,10)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,34) = facY*Tmp0(iTUVP,20)
     enddo
     do iTUVP = 1, 35
      Tmp0(iTUVP,35) = facZ*Tmp0(iTUVP,20)+3*inv2expQ*Tmp0(iTUVP,10)
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,21) = Tmp0(iTUVP,21) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,11) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,22) = Tmp0(iTUVP,22) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,11) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,23) = Tmp0(iTUVP,23) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,11) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,24) = Tmp0(iTUVP,24) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,14) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,25) = Tmp0(iTUVP,25) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,13) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,26) = Tmp0(iTUVP,26) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,16) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,27) = Tmp0(iTUVP,27) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,17) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,28) = Tmp0(iTUVP,28) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,18) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,29) = Tmp0(iTUVP,29) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,19) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1_120(ituvpminus1)
      Tmp0(iTUVP,30) = Tmp0(iTUVP,30) + IfacX1_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,20) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,31) = Tmp0(iTUVP,31) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,17) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,32) = Tmp0(iTUVP,32) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,17) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,33) = Tmp0(iTUVP,33) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,19) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2_120(ituvpminus1)
      Tmp0(iTUVP,34) = Tmp0(iTUVP,34) + IfacX2_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,20) 
     enddo
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3_120(ituvpminus1)
      Tmp0(iTUVP,35) = Tmp0(iTUVP,35) + IfacX3_84(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,20) 
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,21) = Tmp0(iTUVP,21) + pinvq*Tmp0(iTUVplus1,11)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,21) = Tmp0(iTUVP,21) + pinvq*Tmp3(iTUVplus1,11)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,22) = Tmp0(iTUVP,22) + pinvq*Tmp0(iTUVplus1,11)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,22) = Tmp0(iTUVP,22) + pinvq*Tmp3(iTUVplus1,11)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,23) = Tmp0(iTUVP,23) + pinvq*Tmp0(iTUVplus1,11)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,23) = Tmp0(iTUVP,23) + pinvq*Tmp3(iTUVplus1,11)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,24) = Tmp0(iTUVP,24) + pinvq*Tmp0(iTUVplus1,14)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,24) = Tmp0(iTUVP,24) + pinvq*Tmp3(iTUVplus1,14)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,25) = Tmp0(iTUVP,25) + pinvq*Tmp0(iTUVplus1,13)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,25) = Tmp0(iTUVP,25) + pinvq*Tmp3(iTUVplus1,13)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,26) = Tmp0(iTUVP,26) + pinvq*Tmp0(iTUVplus1,16)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,26) = Tmp0(iTUVP,26) + pinvq*Tmp3(iTUVplus1,16)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,27) = Tmp0(iTUVP,27) + pinvq*Tmp0(iTUVplus1,17)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,27) = Tmp0(iTUVP,27) + pinvq*Tmp3(iTUVplus1,17)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,28) = Tmp0(iTUVP,28) + pinvq*Tmp0(iTUVplus1,18)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,28) = Tmp0(iTUVP,28) + pinvq*Tmp3(iTUVplus1,18)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,29) = Tmp0(iTUVP,29) + pinvq*Tmp0(iTUVplus1,19)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,29) = Tmp0(iTUVP,29) + pinvq*Tmp3(iTUVplus1,19)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,30) = Tmp0(iTUVP,30) + pinvq*Tmp0(iTUVplus1,20)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1_120(iTUVP)
      Tmp0(iTUVP,30) = Tmp0(iTUVP,30) + pinvq*Tmp3(iTUVplus1,20)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,31) = Tmp0(iTUVP,31) + pinvq*Tmp0(iTUVplus1,17)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,31) = Tmp0(iTUVP,31) + pinvq*Tmp3(iTUVplus1,17)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,32) = Tmp0(iTUVP,32) + pinvq*Tmp0(iTUVplus1,17)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,32) = Tmp0(iTUVP,32) + pinvq*Tmp3(iTUVplus1,17)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,33) = Tmp0(iTUVP,33) + pinvq*Tmp0(iTUVplus1,19)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,33) = Tmp0(iTUVP,33) + pinvq*Tmp3(iTUVplus1,19)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,34) = Tmp0(iTUVP,34) + pinvq*Tmp0(iTUVplus1,20)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2_120(iTUVP)
      Tmp0(iTUVP,34) = Tmp0(iTUVP,34) + pinvq*Tmp3(iTUVplus1,20)
     enddo
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,35) = Tmp0(iTUVP,35) + pinvq*Tmp0(iTUVplus1,20)
     enddo
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3_120(iTUVP)
      Tmp0(iTUVP,35) = Tmp0(iTUVP,35) + pinvq*Tmp3(iTUVplus1,20)
     enddo
     DO iTUVQ=1, 35
      DO iTUVP=1, 35
        Aux2(iTUVP,iTUVQ,IP) = Aux2(iTUVP,iTUVQ,IP) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP4Q4AtoCSeg
end module
