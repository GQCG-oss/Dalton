MODULE AGC_GPU_OBS_TRMODAtoCSeg1Prim2
 use IchorPrecisionModule
  
 CONTAINS
 subroutine TransferRecurrenceGPUP4Q2AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
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
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
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
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3(iTUVP)
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
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP4Q2AtoCSeg1Prim
 subroutine TransferRecurrenceGPUP4Q3AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
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
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         Tmp2,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
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
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,56
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,84
      iTUVplus1 = TUVindexX3(iTUVP)
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
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP4Q3AtoCSeg1Prim
 subroutine TransferRecurrenceGPUP4Q4AtoCSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
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
  !CARTDIR = 1
  integer,parameter, dimension(120) :: TUVindexX1 = (/ 2,5,6,7,11,12,13,&
          & 14,15,16,21,22,23,24,25,26,27,28,29,30,36,37,38,39,&
          & 40,41,42,43,44,45,46,47,48,49,50,57,58,59,60,61,62,&
          & 63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,85,86,&
          & 87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,&
          & 104,105,106,107,108,109,110,111,112,121,122,123,124,125,126,127,128,&
          & 129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,&
          & 146,147,148,149,150,151,152,153,154,155,156 /)
  !CARTDIR = 2
  integer,parameter, dimension(120) :: TUVindexX2 = (/ 3,6,8,9,12,14,15,&
          & 17,18,19,22,24,25,27,28,29,31,32,33,34,37,39,40,42,&
          & 43,44,46,47,48,49,51,52,53,54,55,58,60,61,63,64,65,&
          & 67,68,69,70,72,73,74,75,76,78,79,80,81,82,83,86,88,&
          & 89,91,92,93,95,96,97,98,100,101,102,103,104,106,107,108,109,&
          & 110,111,113,114,115,116,117,118,119,122,124,125,127,128,129,131,132,&
          & 133,134,136,137,138,139,140,142,143,144,145,146,147,149,150,151,152,&
          & 153,154,155,157,158,159,160,161,162,163,164 /)
  !CARTDIR = 3
  integer,parameter, dimension(120) :: TUVindexX3 = (/ 4,7,9,10,13,15,16,&
          & 18,19,20,23,25,26,28,29,30,32,33,34,35,38,40,41,43,&
          & 44,45,47,48,49,50,52,53,54,55,56,59,61,62,64,65,66,&
          & 68,69,70,71,73,74,75,76,77,79,80,81,82,83,84,87,89,&
          & 90,92,93,94,96,97,98,99,101,102,103,104,105,107,108,109,110,&
          & 111,112,114,115,116,117,118,119,120,123,125,126,128,129,130,132,133,&
          & 134,135,137,138,139,140,141,143,144,145,146,147,148,150,151,152,153,&
          & 154,155,156,158,159,160,161,162,163,164,165 /)
  !CARTDIR = 1
  integer,parameter, dimension(84) :: IfacX1 = (/ 1,2,1,1,3,2,2,&
          & 1,1,1,4,3,3,2,2,2,1,1,1,1,5,4,4,3,&
          & 3,3,2,2,2,2,1,1,1,1,1,6,5,5,4,4,4,&
          & 3,3,3,3,2,2,2,2,2,1,1,1,1,1,1,7,6,&
          & 6,5,5,5,4,4,4,4,3,3,3,3,3,2,2,2,2,&
          & 2,2,1,1,1,1,1,1,1 /)
  !CARTDIR = 2
  integer,parameter, dimension(84) :: IfacX2 = (/ 1,1,2,1,1,2,1,&
          & 3,2,1,1,2,1,3,2,1,4,3,2,1,1,2,1,3,&
          & 2,1,4,3,2,1,5,4,3,2,1,1,2,1,3,2,1,&
          & 4,3,2,1,5,4,3,2,1,6,5,4,3,2,1,1,2,&
          & 1,3,2,1,4,3,2,1,5,4,3,2,1,6,5,4,3,&
          & 2,1,7,6,5,4,3,2,1 /)
  !CARTDIR = 3
  integer,parameter, dimension(84) :: IfacX3 = (/ 1,1,1,2,1,1,2,&
          & 1,2,3,1,1,2,1,2,3,1,2,3,4,1,1,2,1,&
          & 2,3,1,2,3,4,1,2,3,4,5,1,1,2,1,2,3,&
          & 1,2,3,4,1,2,3,4,5,1,2,3,4,5,6,1,1,&
          & 2,1,2,3,1,2,3,4,1,2,3,4,5,1,2,3,4,&
          & 5,6,1,2,3,4,5,6,7 /)
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
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,84
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,84
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,84
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,120
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,120
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 36,120
      iTUVplus1 = TUVindexX3(iTUVP)
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
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp2(iTUVP,5) = Tmp2(iTUVP,5) + IfacX1(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp2(iTUVP,6) = Tmp2(iTUVP,6) + IfacX1(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp2(iTUVP,7) = Tmp2(iTUVP,7) + IfacX1(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp2(iTUVP,8) = Tmp2(iTUVP,8) + IfacX2(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp2(iTUVP,9) = Tmp2(iTUVP,9) + IfacX2(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 36,56
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp2(iTUVP,10) = Tmp2(iTUVP,10) + IfacX3(ituvpminus1)*inv2expQ*Tmp1(ituvpminus1,4) 
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
     do iTUVP = 36,84
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
     do iTUVP = 36,84
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
     do iTUVP = 36,84
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
     do iTUVP = 36,84
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
     do iTUVP = 36,84
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
     do iTUVP = 36,84
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
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp3(iTUVP,11) = Tmp3(iTUVP,11) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp3(iTUVP,12) = Tmp3(iTUVP,12) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp3(iTUVP,13) = Tmp3(iTUVP,13) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp3(iTUVP,14) = Tmp3(iTUVP,14) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp3(iTUVP,15) = Tmp3(iTUVP,15) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,9) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp3(iTUVP,16) = Tmp3(iTUVP,16) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp3(iTUVP,17) = Tmp3(iTUVP,17) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp3(iTUVP,18) = Tmp3(iTUVP,18) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp3(iTUVP,19) = Tmp3(iTUVP,19) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 21,35
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp3(iTUVP,20) = Tmp3(iTUVP,20) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,10) 
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
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp3(iTUVP,11) = Tmp3(iTUVP,11) + pinvq*Tmp2(iTUVplus1,5)
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
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp3(iTUVP,12) = Tmp3(iTUVP,12) + pinvq*Tmp2(iTUVplus1,5)
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
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp3(iTUVP,13) = Tmp3(iTUVP,13) + pinvq*Tmp2(iTUVplus1,5)
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
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp3(iTUVP,14) = Tmp3(iTUVP,14) + pinvq*Tmp2(iTUVplus1,8)
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
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp3(iTUVP,15) = Tmp3(iTUVP,15) + pinvq*Tmp2(iTUVplus1,9)
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
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp3(iTUVP,16) = Tmp3(iTUVP,16) + pinvq*Tmp2(iTUVplus1,10)
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
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp3(iTUVP,17) = Tmp3(iTUVP,17) + pinvq*Tmp2(iTUVplus1,8)
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
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp3(iTUVP,18) = Tmp3(iTUVP,18) + pinvq*Tmp2(iTUVplus1,8)
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
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp3(iTUVP,19) = Tmp3(iTUVP,19) + pinvq*Tmp2(iTUVplus1,10)
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
     do iTUVP = 36,56
      iTUVplus1 = TUVindexX3(iTUVP)
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
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,21) = Tmp0(iTUVP,21) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,11) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,22) = Tmp0(iTUVP,22) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,11) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,23) = Tmp0(iTUVP,23) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,11) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,24) = Tmp0(iTUVP,24) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,14) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,25) = Tmp0(iTUVP,25) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,13) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,26) = Tmp0(iTUVP,26) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,16) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,27) = Tmp0(iTUVP,27) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,17) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,28) = Tmp0(iTUVP,28) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,18) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,29) = Tmp0(iTUVP,29) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,19) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,30) = Tmp0(iTUVP,30) + IfacX1(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,20) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,31) = Tmp0(iTUVP,31) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,17) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,32) = Tmp0(iTUVP,32) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,17) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,33) = Tmp0(iTUVP,33) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,19) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,34) = Tmp0(iTUVP,34) + IfacX2(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,20) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,35) = Tmp0(iTUVP,35) + IfacX3(ituvpminus1)*inv2expQ*Tmp0(ituvpminus1,20) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,21) = Tmp0(iTUVP,21) + pinvq*Tmp0(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,21) = Tmp0(iTUVP,21) + pinvq*Tmp3(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,22) = Tmp0(iTUVP,22) + pinvq*Tmp0(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,22) = Tmp0(iTUVP,22) + pinvq*Tmp3(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,23) = Tmp0(iTUVP,23) + pinvq*Tmp0(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,23) = Tmp0(iTUVP,23) + pinvq*Tmp3(iTUVplus1,11)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,24) = Tmp0(iTUVP,24) + pinvq*Tmp0(iTUVplus1,14)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,24) = Tmp0(iTUVP,24) + pinvq*Tmp3(iTUVplus1,14)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,25) = Tmp0(iTUVP,25) + pinvq*Tmp0(iTUVplus1,13)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,25) = Tmp0(iTUVP,25) + pinvq*Tmp3(iTUVplus1,13)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,26) = Tmp0(iTUVP,26) + pinvq*Tmp0(iTUVplus1,16)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,26) = Tmp0(iTUVP,26) + pinvq*Tmp3(iTUVplus1,16)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,27) = Tmp0(iTUVP,27) + pinvq*Tmp0(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,27) = Tmp0(iTUVP,27) + pinvq*Tmp3(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,28) = Tmp0(iTUVP,28) + pinvq*Tmp0(iTUVplus1,18)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,28) = Tmp0(iTUVP,28) + pinvq*Tmp3(iTUVplus1,18)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,29) = Tmp0(iTUVP,29) + pinvq*Tmp0(iTUVplus1,19)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,29) = Tmp0(iTUVP,29) + pinvq*Tmp3(iTUVplus1,19)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,30) = Tmp0(iTUVP,30) + pinvq*Tmp0(iTUVplus1,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,30) = Tmp0(iTUVP,30) + pinvq*Tmp3(iTUVplus1,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,31) = Tmp0(iTUVP,31) + pinvq*Tmp0(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,31) = Tmp0(iTUVP,31) + pinvq*Tmp3(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,32) = Tmp0(iTUVP,32) + pinvq*Tmp0(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,32) = Tmp0(iTUVP,32) + pinvq*Tmp3(iTUVplus1,17)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,33) = Tmp0(iTUVP,33) + pinvq*Tmp0(iTUVplus1,19)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,33) = Tmp0(iTUVP,33) + pinvq*Tmp3(iTUVplus1,19)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,34) = Tmp0(iTUVP,34) + pinvq*Tmp0(iTUVplus1,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,34) = Tmp0(iTUVP,34) + pinvq*Tmp3(iTUVplus1,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,35) = Tmp0(iTUVP,35) + pinvq*Tmp0(iTUVplus1,20)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3(iTUVP)
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
