MODULE AGC_GPU_OBS_TRMODAtoDSegQ2
 use IchorPrecisionModule
  
 CONTAINS
 subroutine TransferRecurrenceGPUP3Q3AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   84)
  real(realk),intent(inout) :: Aux2(nPrimP*nPasses,   20,   20)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
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
  !CARTDIR = 1
  integer,parameter, dimension(56) :: TUVindexX1 = (/ 2,5,6,7,11,12,13,&
          & 14,15,16,21,22,23,24,25,26,27,28,29,30,36,37,38,39,&
          & 40,41,42,43,44,45,46,47,48,49,50,57,58,59,60,61,62,&
          & 63,64,65,66,67,68,69,70,71,72,73,74,75,76,77 /)
  !CARTDIR = 2
  integer,parameter, dimension(56) :: TUVindexX2 = (/ 3,6,8,9,12,14,15,&
          & 17,18,19,22,24,25,27,28,29,31,32,33,34,37,39,40,42,&
          & 43,44,46,47,48,49,51,52,53,54,55,58,60,61,63,64,65,&
          & 67,68,69,70,72,73,74,75,76,78,79,80,81,82,83 /)
  !CARTDIR = 3
  integer,parameter, dimension(56) :: TUVindexX3 = (/ 4,7,9,10,13,15,16,&
          & 18,19,20,23,25,26,28,29,30,32,33,34,35,38,40,41,43,&
          & 44,45,47,48,49,50,52,53,54,55,56,59,61,62,64,65,66,&
          & 68,69,70,71,73,74,75,76,77,79,80,81,82,83,84 /)
  !CARTDIR = 1
  integer,parameter, dimension(35) :: IfacX1 = (/ 1,2,1,1,3,2,2,&
          & 1,1,1,4,3,3,2,2,2,1,1,1,1,5,4,4,3,&
          & 3,3,2,2,2,2,1,1,1,1,1 /)
  !CARTDIR = 2
  integer,parameter, dimension(35) :: IfacX2 = (/ 1,1,2,1,1,2,1,&
          & 3,2,1,1,2,1,3,2,1,4,3,2,1,1,2,1,3,&
          & 2,1,4,3,2,1,5,4,3,2,1 /)
  !CARTDIR = 3
  integer,parameter, dimension(35) :: IfacX3 = (/ 1,1,1,2,1,1,2,&
          & 1,2,3,1,1,2,1,2,3,1,2,3,4,1,1,2,1,&
          & 2,3,1,2,3,4,1,2,3,4,5 /)
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,nPrimP,Aux2) ASYNC(iASync)
  DO iP = 1,nPrimP*nPasses
   DO iTUVQ=1, 20
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
!$ACC         Tmp2,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPrimP*nPasses
   DO iPrimQ=1, nPrimQ
    iPrimP = iP - ((iP-1)/nPrimP)*nPrimP
    iPassP = (iP-1)/nPrimP + 1
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
     do iTUVP =  21, 56
      Tmp1(iTUVP,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 56
      Tmp1(iTUVP,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 56
      Tmp1(iTUVP,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,35
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,35
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,56
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,56
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,56
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
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
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp2(iTUVP,5) = facX*Tmp1(iTUVP,2)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp2(iTUVP,6) = facX*Tmp1(iTUVP,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp2(iTUVP,7) = facX*Tmp1(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp2(iTUVP,8) = facY*Tmp1(iTUVP,3)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp2(iTUVP,9) = facY*Tmp1(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp2(iTUVP,10) = facZ*Tmp1(iTUVP,4)+ inv2expQ*Aux(iPrimQ,iPrimP,iPassP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp0(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
 ! Building for Angular momentum Jq = 3
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,11) = facX*Tmp0(iTUVP,5)+2*inv2expQ*Tmp0(iTUVP,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,12) = facY*Tmp0(iTUVP,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,13) = facZ*Tmp0(iTUVP,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,14) = facX*Tmp0(iTUVP,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,15) = facX*Tmp0(iTUVP,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,16) = facX*Tmp0(iTUVP,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,17) = facY*Tmp0(iTUVP,8)+2*inv2expQ*Tmp0(iTUVP,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,18) = facZ*Tmp0(iTUVP,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,19) = facY*Tmp0(iTUVP,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,20) = facZ*Tmp0(iTUVP,10)+2*inv2expQ*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,9) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp0(iTUVplus1,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp2(iTUVplus1,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,10
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 11,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     DO iTUVQ=1, 20
!$ACC LOOP SEQ
      DO iTUVP=1, 20
        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimP=1, nPrimP
  ENDDO !iP = 1,nPrimQ*nPasses
 end subroutine TransferRecurrenceGPUP3Q3AtoDSegQ
 subroutine TransferRecurrenceGPUP4Q2AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   84)
  real(realk),intent(inout) :: Aux2(nPrimP*nPasses,   35,   10)
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
  !CARTDIR = 1
  integer,parameter, dimension(56) :: TUVindexX1 = (/ 2,5,6,7,11,12,13,&
          & 14,15,16,21,22,23,24,25,26,27,28,29,30,36,37,38,39,&
          & 40,41,42,43,44,45,46,47,48,49,50,57,58,59,60,61,62,&
          & 63,64,65,66,67,68,69,70,71,72,73,74,75,76,77 /)
  !CARTDIR = 2
  integer,parameter, dimension(56) :: TUVindexX2 = (/ 3,6,8,9,12,14,15,&
          & 17,18,19,22,24,25,27,28,29,31,32,33,34,37,39,40,42,&
          & 43,44,46,47,48,49,51,52,53,54,55,58,60,61,63,64,65,&
          & 67,68,69,70,72,73,74,75,76,78,79,80,81,82,83 /)
  !CARTDIR = 3
  integer,parameter, dimension(56) :: TUVindexX3 = (/ 4,7,9,10,13,15,16,&
          & 18,19,20,23,25,26,28,29,30,32,33,34,35,38,40,41,43,&
          & 44,45,47,48,49,50,52,53,54,55,56,59,61,62,64,65,66,&
          & 68,69,70,71,73,74,75,76,77,79,80,81,82,83,84 /)
  !CARTDIR = 1
  integer,parameter, dimension(35) :: IfacX1 = (/ 1,2,1,1,3,2,2,&
          & 1,1,1,4,3,3,2,2,2,1,1,1,1,5,4,4,3,&
          & 3,3,2,2,2,2,1,1,1,1,1 /)
  !CARTDIR = 2
  integer,parameter, dimension(35) :: IfacX2 = (/ 1,1,2,1,1,2,1,&
          & 3,2,1,1,2,1,3,2,1,4,3,2,1,1,2,1,3,&
          & 2,1,4,3,2,1,5,4,3,2,1 /)
  !CARTDIR = 3
  integer,parameter, dimension(35) :: IfacX3 = (/ 1,1,1,2,1,1,2,&
          & 1,2,3,1,1,2,1,2,3,1,2,3,4,1,1,2,1,&
          & 2,3,1,2,3,4,1,2,3,4,5 /)
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,nPrimP,Aux2) ASYNC(iASync)
  DO iP = 1,nPrimP*nPasses
   DO iTUVQ=1, 10
    DO iTUVP=1, 35
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
  DO iP = 1,nPrimP*nPasses
   DO iPrimQ=1, nPrimQ
    iPrimP = iP - ((iP-1)/nPrimP)*nPrimP
    iPassP = (iP-1)/nPrimP + 1
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
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3(iTUVP)
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
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp0(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     DO iTUVQ=1, 10
!$ACC LOOP SEQ
      DO iTUVP=1, 35
        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimP=1, nPrimP
  ENDDO !iP = 1,nPrimQ*nPasses
 end subroutine TransferRecurrenceGPUP4Q2AtoDSegQ
 subroutine TransferRecurrenceGPUP4Q3AtoDSegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Bexp(nPrimB),Cexp(nPrimC)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,  120)
  real(realk),intent(inout) :: Aux2(nPrimP*nPasses,   35,   20)
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
  !CARTDIR = 1
  integer,parameter, dimension(84) :: TUVindexX1 = (/ 2,5,6,7,11,12,13,&
          & 14,15,16,21,22,23,24,25,26,27,28,29,30,36,37,38,39,&
          & 40,41,42,43,44,45,46,47,48,49,50,57,58,59,60,61,62,&
          & 63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,85,86,&
          & 87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,&
          & 104,105,106,107,108,109,110,111,112 /)
  !CARTDIR = 2
  integer,parameter, dimension(84) :: TUVindexX2 = (/ 3,6,8,9,12,14,15,&
          & 17,18,19,22,24,25,27,28,29,31,32,33,34,37,39,40,42,&
          & 43,44,46,47,48,49,51,52,53,54,55,58,60,61,63,64,65,&
          & 67,68,69,70,72,73,74,75,76,78,79,80,81,82,83,86,88,&
          & 89,91,92,93,95,96,97,98,100,101,102,103,104,106,107,108,109,&
          & 110,111,113,114,115,116,117,118,119 /)
  !CARTDIR = 3
  integer,parameter, dimension(84) :: TUVindexX3 = (/ 4,7,9,10,13,15,16,&
          & 18,19,20,23,25,26,28,29,30,32,33,34,35,38,40,41,43,&
          & 44,45,47,48,49,50,52,53,54,55,56,59,61,62,64,65,66,&
          & 68,69,70,71,73,74,75,76,77,79,80,81,82,83,84,87,89,&
          & 90,92,93,94,96,97,98,99,101,102,103,104,105,107,108,109,110,&
          & 111,112,114,115,116,117,118,119,120 /)
  !CARTDIR = 1
  integer,parameter, dimension(56) :: IfacX1 = (/ 1,2,1,1,3,2,2,&
          & 1,1,1,4,3,3,2,2,2,1,1,1,1,5,4,4,3,&
          & 3,3,2,2,2,2,1,1,1,1,1,6,5,5,4,4,4,&
          & 3,3,3,3,2,2,2,2,2,1,1,1,1,1,1 /)
  !CARTDIR = 2
  integer,parameter, dimension(56) :: IfacX2 = (/ 1,1,2,1,1,2,1,&
          & 3,2,1,1,2,1,3,2,1,4,3,2,1,1,2,1,3,&
          & 2,1,4,3,2,1,5,4,3,2,1,1,2,1,3,2,1,&
          & 4,3,2,1,5,4,3,2,1,6,5,4,3,2,1 /)
  !CARTDIR = 3
  integer,parameter, dimension(56) :: IfacX3 = (/ 1,1,1,2,1,1,2,&
          & 1,2,3,1,1,2,1,2,3,1,2,3,4,1,1,2,1,&
          & 2,3,1,2,3,4,1,2,3,4,5,1,1,2,1,2,3,&
          & 1,2,3,4,1,2,3,4,5,1,2,3,4,5,6 /)
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,nPrimP,Aux2) ASYNC(iASync)
  DO iP = 1,nPrimP*nPasses
   DO iTUVQ=1, 20
    DO iTUVP=1, 35
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
!$ACC         Tmp2,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Bexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPrimP*nPasses
   DO iPrimQ=1, nPrimQ
    iPrimP = iP - ((iP-1)/nPrimP)*nPrimP
    iPassP = (iP-1)/nPrimP + 1
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
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(iPrimQ,iPrimP,iPassP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX3(iTUVP)
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
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp0(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,5) = Tmp0(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + pinvq*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,6) = Tmp0(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,7) = Tmp0(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,8) = Tmp0(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + pinvq*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,9) = Tmp0(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,10) = Tmp0(iTUVP,10) + pinvq*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3(iTUVP)
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
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,9) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,11) = Tmp0(iTUVP,11) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,12) = Tmp0(iTUVP,12) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,13) = Tmp0(iTUVP,13) + pinvq*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,14) = Tmp0(iTUVP,14) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp0(iTUVplus1,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,15) = Tmp0(iTUVP,15) + pinvq*Tmp2(iTUVplus1,9)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,16) = Tmp0(iTUVP,16) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,17) = Tmp0(iTUVP,17) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,18) = Tmp0(iTUVP,18) + pinvq*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,19) = Tmp0(iTUVP,19) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,20) = Tmp0(iTUVP,20) + pinvq*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     DO iTUVQ=1, 20
!$ACC LOOP SEQ
      DO iTUVP=1, 35
        Aux2(IP,iTUVP,iTUVQ) = Aux2(IP,iTUVP,iTUVQ) + Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
   ENDDO !iPrimP=1, nPrimP
  ENDDO !iP = 1,nPrimQ*nPasses
 end subroutine TransferRecurrenceGPUP4Q3AtoDSegQ
end module
