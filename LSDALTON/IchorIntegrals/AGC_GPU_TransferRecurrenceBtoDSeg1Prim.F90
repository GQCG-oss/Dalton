MODULE AGC_GPU_OBS_TRMODBtoDSeg1Prim
 use IchorPrecisionModule
  
 CONTAINS
 subroutine TransferRecurrenceGPUP1Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Aexp(1),Cexp(1)
  real(realk),intent(in) :: Aux(nPasses,   10)
  real(realk),intent(inout) :: Aux2(nPasses,    4,    4)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0(  4,  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expAX,expAY,expAZ,&
!$ACC         Tmp0,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Aexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux)
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
    expAX = -Aexp(1)*Xab
    expAY = -Aexp(1)*Yab
    expAZ = -Aexp(1)*Zab
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     facX = -(expAX-Cexp(1)*Xcd)*invexpQ
     facY = -(expAY-Cexp(1)*Ycd)*invexpQ
     facZ = -(expAZ-Cexp(1)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
!$ACC LOOP SEQ
     DO iTUVP=1,  4
      Tmp0(iTUVP,1) = Aux(IP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1,  4
      Tmp0(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,  4
      Tmp0(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,  4
      Tmp0(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(IP,1) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(IP,1) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(IP,1) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(IP,2)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(IP,5)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(IP,6)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(IP,7)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(IP,3)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(IP,6)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(IP,8)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(IP,9)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(IP,4)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(IP,7)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(IP,9)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(IP,10)
!$ACC LOOP SEQ
     DO iTUVQ=1,  4
!$ACC LOOP SEQ
      DO iTUVP=1,  4
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP1Q1BtoDSeg1Prim
 subroutine TransferRecurrenceGPUP2Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Aexp(1),Cexp(1)
  real(realk),intent(in) :: Aux(nPasses,   20)
  real(realk),intent(inout) :: Aux2(nPasses,   10,    4)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0( 10,  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expAX,expAY,expAZ,&
!$ACC         Tmp0,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Aexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux)
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
    expAX = -Aexp(1)*Xab
    expAY = -Aexp(1)*Yab
    expAZ = -Aexp(1)*Zab
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     facX = -(expAX-Cexp(1)*Xcd)*invexpQ
     facY = -(expAY-Cexp(1)*Ycd)*invexpQ
     facZ = -(expAZ-Cexp(1)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
!$ACC LOOP SEQ
     DO iTUVP=1, 10
      Tmp0(iTUVP,1) = Aux(IP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(IP,1) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expQ*Aux(IP,2) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expQ*Aux(IP,3) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expQ*Aux(IP,4) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(IP,1) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expQ*Aux(IP,2) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expQ*Aux(IP,3) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expQ*Aux(IP,4) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(IP,1) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expQ*Aux(IP,2) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expQ*Aux(IP,3) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expQ*Aux(IP,4) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(IP,2)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(IP,5)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(IP,6)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(IP,7)
     Tmp0(5,2) = Tmp0(5,2) + pinvq*Aux(IP,11)
     Tmp0(6,2) = Tmp0(6,2) + pinvq*Aux(IP,12)
     Tmp0(7,2) = Tmp0(7,2) + pinvq*Aux(IP,13)
     Tmp0(8,2) = Tmp0(8,2) + pinvq*Aux(IP,14)
     Tmp0(9,2) = Tmp0(9,2) + pinvq*Aux(IP,15)
     Tmp0(10,2) = Tmp0(10,2) + pinvq*Aux(IP,16)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(IP,3)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(IP,6)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(IP,8)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(IP,9)
     Tmp0(5,3) = Tmp0(5,3) + pinvq*Aux(IP,12)
     Tmp0(6,3) = Tmp0(6,3) + pinvq*Aux(IP,14)
     Tmp0(7,3) = Tmp0(7,3) + pinvq*Aux(IP,15)
     Tmp0(8,3) = Tmp0(8,3) + pinvq*Aux(IP,17)
     Tmp0(9,3) = Tmp0(9,3) + pinvq*Aux(IP,18)
     Tmp0(10,3) = Tmp0(10,3) + pinvq*Aux(IP,19)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(IP,4)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(IP,7)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(IP,9)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(IP,10)
     Tmp0(5,4) = Tmp0(5,4) + pinvq*Aux(IP,13)
     Tmp0(6,4) = Tmp0(6,4) + pinvq*Aux(IP,15)
     Tmp0(7,4) = Tmp0(7,4) + pinvq*Aux(IP,16)
     Tmp0(8,4) = Tmp0(8,4) + pinvq*Aux(IP,18)
     Tmp0(9,4) = Tmp0(9,4) + pinvq*Aux(IP,19)
     Tmp0(10,4) = Tmp0(10,4) + pinvq*Aux(IP,20)
!$ACC LOOP SEQ
     DO iTUVQ=1,  4
!$ACC LOOP SEQ
      DO iTUVP=1, 10
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP2Q1BtoDSeg1Prim
 subroutine TransferRecurrenceGPUP2Q2BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Aexp(1),Cexp(1)
  real(realk),intent(in) :: Aux(nPasses,   35)
  real(realk),intent(inout) :: Aux2(nPasses,   10,   10)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0( 10, 10)
  real(realk) :: Tmp1( 11: 20,  2:  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expAX,expAY,expAZ,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Aexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux)
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
    expAX = -Aexp(1)*Xab
    expAY = -Aexp(1)*Yab
    expAZ = -Aexp(1)*Zab
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     facX = -(expAX-Cexp(1)*Xcd)*invexpQ
     facY = -(expAY-Cexp(1)*Ycd)*invexpQ
     facZ = -(expAZ-Cexp(1)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
!$ACC LOOP SEQ
     DO iTUVP=1, 10
      Tmp0(iTUVP,1) = Aux(IP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  11, 20
      Tmp1(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  11, 20
      Tmp1(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  11, 20
      Tmp1(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(IP,1) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expQ*Aux(IP,2) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expQ*Aux(IP,3) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expQ*Aux(IP,4) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(IP,1) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expQ*Aux(IP,2) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expQ*Aux(IP,3) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expQ*Aux(IP,4) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(IP,1) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expQ*Aux(IP,2) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expQ*Aux(IP,3) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expQ*Aux(IP,4) 
     Tmp1(11,2) = Tmp1(11,2) + 3*inv2expQ*Aux(IP,5) 
     Tmp1(12,2) = Tmp1(12,2) + 2*inv2expQ*Aux(IP,6) 
     Tmp1(13,2) = Tmp1(13,2) + 2*inv2expQ*Aux(IP,7) 
     Tmp1(14,2) = Tmp1(14,2) + inv2expQ*Aux(IP,8) 
     Tmp1(15,2) = Tmp1(15,2) + inv2expQ*Aux(IP,9) 
     Tmp1(16,2) = Tmp1(16,2) + inv2expQ*Aux(IP,10) 
     Tmp1(12,3) = Tmp1(12,3) + inv2expQ*Aux(IP,5) 
     Tmp1(14,3) = Tmp1(14,3) + 2*inv2expQ*Aux(IP,6) 
     Tmp1(15,3) = Tmp1(15,3) + inv2expQ*Aux(IP,7) 
     Tmp1(17,3) = Tmp1(17,3) + 3*inv2expQ*Aux(IP,8) 
     Tmp1(18,3) = Tmp1(18,3) + 2*inv2expQ*Aux(IP,9) 
     Tmp1(19,3) = Tmp1(19,3) + inv2expQ*Aux(IP,10) 
     Tmp1(13,4) = Tmp1(13,4) + inv2expQ*Aux(IP,5) 
     Tmp1(15,4) = Tmp1(15,4) + inv2expQ*Aux(IP,6) 
     Tmp1(16,4) = Tmp1(16,4) + 2*inv2expQ*Aux(IP,7) 
     Tmp1(18,4) = Tmp1(18,4) + inv2expQ*Aux(IP,8) 
     Tmp1(19,4) = Tmp1(19,4) + 2*inv2expQ*Aux(IP,9) 
     Tmp1(20,4) = Tmp1(20,4) + 3*inv2expQ*Aux(IP,10) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(IP,2)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(IP,5)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(IP,6)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(IP,7)
     Tmp0(5,2) = Tmp0(5,2) + pinvq*Aux(IP,11)
     Tmp0(6,2) = Tmp0(6,2) + pinvq*Aux(IP,12)
     Tmp0(7,2) = Tmp0(7,2) + pinvq*Aux(IP,13)
     Tmp0(8,2) = Tmp0(8,2) + pinvq*Aux(IP,14)
     Tmp0(9,2) = Tmp0(9,2) + pinvq*Aux(IP,15)
     Tmp0(10,2) = Tmp0(10,2) + pinvq*Aux(IP,16)
     Tmp1(11,2) = Tmp1(11,2) + pinvq*Aux(IP,21)
     Tmp1(12,2) = Tmp1(12,2) + pinvq*Aux(IP,22)
     Tmp1(13,2) = Tmp1(13,2) + pinvq*Aux(IP,23)
     Tmp1(14,2) = Tmp1(14,2) + pinvq*Aux(IP,24)
     Tmp1(15,2) = Tmp1(15,2) + pinvq*Aux(IP,25)
     Tmp1(16,2) = Tmp1(16,2) + pinvq*Aux(IP,26)
     Tmp1(17,2) = Tmp1(17,2) + pinvq*Aux(IP,27)
     Tmp1(18,2) = Tmp1(18,2) + pinvq*Aux(IP,28)
     Tmp1(19,2) = Tmp1(19,2) + pinvq*Aux(IP,29)
     Tmp1(20,2) = Tmp1(20,2) + pinvq*Aux(IP,30)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(IP,3)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(IP,6)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(IP,8)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(IP,9)
     Tmp0(5,3) = Tmp0(5,3) + pinvq*Aux(IP,12)
     Tmp0(6,3) = Tmp0(6,3) + pinvq*Aux(IP,14)
     Tmp0(7,3) = Tmp0(7,3) + pinvq*Aux(IP,15)
     Tmp0(8,3) = Tmp0(8,3) + pinvq*Aux(IP,17)
     Tmp0(9,3) = Tmp0(9,3) + pinvq*Aux(IP,18)
     Tmp0(10,3) = Tmp0(10,3) + pinvq*Aux(IP,19)
     Tmp1(11,3) = Tmp1(11,3) + pinvq*Aux(IP,22)
     Tmp1(12,3) = Tmp1(12,3) + pinvq*Aux(IP,24)
     Tmp1(13,3) = Tmp1(13,3) + pinvq*Aux(IP,25)
     Tmp1(14,3) = Tmp1(14,3) + pinvq*Aux(IP,27)
     Tmp1(15,3) = Tmp1(15,3) + pinvq*Aux(IP,28)
     Tmp1(16,3) = Tmp1(16,3) + pinvq*Aux(IP,29)
     Tmp1(17,3) = Tmp1(17,3) + pinvq*Aux(IP,31)
     Tmp1(18,3) = Tmp1(18,3) + pinvq*Aux(IP,32)
     Tmp1(19,3) = Tmp1(19,3) + pinvq*Aux(IP,33)
     Tmp1(20,3) = Tmp1(20,3) + pinvq*Aux(IP,34)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(IP,4)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(IP,7)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(IP,9)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(IP,10)
     Tmp0(5,4) = Tmp0(5,4) + pinvq*Aux(IP,13)
     Tmp0(6,4) = Tmp0(6,4) + pinvq*Aux(IP,15)
     Tmp0(7,4) = Tmp0(7,4) + pinvq*Aux(IP,16)
     Tmp0(8,4) = Tmp0(8,4) + pinvq*Aux(IP,18)
     Tmp0(9,4) = Tmp0(9,4) + pinvq*Aux(IP,19)
     Tmp0(10,4) = Tmp0(10,4) + pinvq*Aux(IP,20)
     Tmp1(11,4) = Tmp1(11,4) + pinvq*Aux(IP,23)
     Tmp1(12,4) = Tmp1(12,4) + pinvq*Aux(IP,25)
     Tmp1(13,4) = Tmp1(13,4) + pinvq*Aux(IP,26)
     Tmp1(14,4) = Tmp1(14,4) + pinvq*Aux(IP,28)
     Tmp1(15,4) = Tmp1(15,4) + pinvq*Aux(IP,29)
     Tmp1(16,4) = Tmp1(16,4) + pinvq*Aux(IP,30)
     Tmp1(17,4) = Tmp1(17,4) + pinvq*Aux(IP,32)
     Tmp1(18,4) = Tmp1(18,4) + pinvq*Aux(IP,33)
     Tmp1(19,4) = Tmp1(19,4) + pinvq*Aux(IP,34)
     Tmp1(20,4) = Tmp1(20,4) + pinvq*Aux(IP,35)
 ! Building for Angular momentum Jq = 2
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(IP,iTUVP)
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
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 10
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(IP,iTUVP)
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
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP2Q2BtoDSeg1Prim
 subroutine TransferRecurrenceGPUP3Q1BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Aexp(1),Cexp(1)
  real(realk),intent(in) :: Aux(nPasses,   35)
  real(realk),intent(inout) :: Aux2(nPasses,   20,    4)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0( 20,  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expAX,expAY,expAZ,&
!$ACC         Tmp0,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Aexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux)
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
    expAX = -Aexp(1)*Xab
    expAY = -Aexp(1)*Yab
    expAZ = -Aexp(1)*Zab
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     facX = -(expAX-Cexp(1)*Xcd)*invexpQ
     facY = -(expAY-Cexp(1)*Ycd)*invexpQ
     facZ = -(expAZ-Cexp(1)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
!$ACC LOOP SEQ
     DO iTUVP=1, 20
      Tmp0(iTUVP,1) = Aux(IP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expQ*Aux(IP,1) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expQ*Aux(IP,2) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expQ*Aux(IP,3) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expQ*Aux(IP,4) 
     Tmp0(11,2) = Tmp0(11,2) + 3*inv2expQ*Aux(IP,5) 
     Tmp0(12,2) = Tmp0(12,2) + 2*inv2expQ*Aux(IP,6) 
     Tmp0(13,2) = Tmp0(13,2) + 2*inv2expQ*Aux(IP,7) 
     Tmp0(14,2) = Tmp0(14,2) + inv2expQ*Aux(IP,8) 
     Tmp0(15,2) = Tmp0(15,2) + inv2expQ*Aux(IP,9) 
     Tmp0(16,2) = Tmp0(16,2) + inv2expQ*Aux(IP,10) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expQ*Aux(IP,1) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expQ*Aux(IP,2) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expQ*Aux(IP,3) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expQ*Aux(IP,4) 
     Tmp0(12,3) = Tmp0(12,3) + inv2expQ*Aux(IP,5) 
     Tmp0(14,3) = Tmp0(14,3) + 2*inv2expQ*Aux(IP,6) 
     Tmp0(15,3) = Tmp0(15,3) + inv2expQ*Aux(IP,7) 
     Tmp0(17,3) = Tmp0(17,3) + 3*inv2expQ*Aux(IP,8) 
     Tmp0(18,3) = Tmp0(18,3) + 2*inv2expQ*Aux(IP,9) 
     Tmp0(19,3) = Tmp0(19,3) + inv2expQ*Aux(IP,10) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expQ*Aux(IP,1) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expQ*Aux(IP,2) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expQ*Aux(IP,3) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expQ*Aux(IP,4) 
     Tmp0(13,4) = Tmp0(13,4) + inv2expQ*Aux(IP,5) 
     Tmp0(15,4) = Tmp0(15,4) + inv2expQ*Aux(IP,6) 
     Tmp0(16,4) = Tmp0(16,4) + 2*inv2expQ*Aux(IP,7) 
     Tmp0(18,4) = Tmp0(18,4) + inv2expQ*Aux(IP,8) 
     Tmp0(19,4) = Tmp0(19,4) + 2*inv2expQ*Aux(IP,9) 
     Tmp0(20,4) = Tmp0(20,4) + 3*inv2expQ*Aux(IP,10) 
     Tmp0(1,2) = Tmp0(1,2) + pinvq*Aux(IP,2)
     Tmp0(2,2) = Tmp0(2,2) + pinvq*Aux(IP,5)
     Tmp0(3,2) = Tmp0(3,2) + pinvq*Aux(IP,6)
     Tmp0(4,2) = Tmp0(4,2) + pinvq*Aux(IP,7)
     Tmp0(5,2) = Tmp0(5,2) + pinvq*Aux(IP,11)
     Tmp0(6,2) = Tmp0(6,2) + pinvq*Aux(IP,12)
     Tmp0(7,2) = Tmp0(7,2) + pinvq*Aux(IP,13)
     Tmp0(8,2) = Tmp0(8,2) + pinvq*Aux(IP,14)
     Tmp0(9,2) = Tmp0(9,2) + pinvq*Aux(IP,15)
     Tmp0(10,2) = Tmp0(10,2) + pinvq*Aux(IP,16)
     Tmp0(11,2) = Tmp0(11,2) + pinvq*Aux(IP,21)
     Tmp0(12,2) = Tmp0(12,2) + pinvq*Aux(IP,22)
     Tmp0(13,2) = Tmp0(13,2) + pinvq*Aux(IP,23)
     Tmp0(14,2) = Tmp0(14,2) + pinvq*Aux(IP,24)
     Tmp0(15,2) = Tmp0(15,2) + pinvq*Aux(IP,25)
     Tmp0(16,2) = Tmp0(16,2) + pinvq*Aux(IP,26)
     Tmp0(17,2) = Tmp0(17,2) + pinvq*Aux(IP,27)
     Tmp0(18,2) = Tmp0(18,2) + pinvq*Aux(IP,28)
     Tmp0(19,2) = Tmp0(19,2) + pinvq*Aux(IP,29)
     Tmp0(20,2) = Tmp0(20,2) + pinvq*Aux(IP,30)
     Tmp0(1,3) = Tmp0(1,3) + pinvq*Aux(IP,3)
     Tmp0(2,3) = Tmp0(2,3) + pinvq*Aux(IP,6)
     Tmp0(3,3) = Tmp0(3,3) + pinvq*Aux(IP,8)
     Tmp0(4,3) = Tmp0(4,3) + pinvq*Aux(IP,9)
     Tmp0(5,3) = Tmp0(5,3) + pinvq*Aux(IP,12)
     Tmp0(6,3) = Tmp0(6,3) + pinvq*Aux(IP,14)
     Tmp0(7,3) = Tmp0(7,3) + pinvq*Aux(IP,15)
     Tmp0(8,3) = Tmp0(8,3) + pinvq*Aux(IP,17)
     Tmp0(9,3) = Tmp0(9,3) + pinvq*Aux(IP,18)
     Tmp0(10,3) = Tmp0(10,3) + pinvq*Aux(IP,19)
     Tmp0(11,3) = Tmp0(11,3) + pinvq*Aux(IP,22)
     Tmp0(12,3) = Tmp0(12,3) + pinvq*Aux(IP,24)
     Tmp0(13,3) = Tmp0(13,3) + pinvq*Aux(IP,25)
     Tmp0(14,3) = Tmp0(14,3) + pinvq*Aux(IP,27)
     Tmp0(15,3) = Tmp0(15,3) + pinvq*Aux(IP,28)
     Tmp0(16,3) = Tmp0(16,3) + pinvq*Aux(IP,29)
     Tmp0(17,3) = Tmp0(17,3) + pinvq*Aux(IP,31)
     Tmp0(18,3) = Tmp0(18,3) + pinvq*Aux(IP,32)
     Tmp0(19,3) = Tmp0(19,3) + pinvq*Aux(IP,33)
     Tmp0(20,3) = Tmp0(20,3) + pinvq*Aux(IP,34)
     Tmp0(1,4) = Tmp0(1,4) + pinvq*Aux(IP,4)
     Tmp0(2,4) = Tmp0(2,4) + pinvq*Aux(IP,7)
     Tmp0(3,4) = Tmp0(3,4) + pinvq*Aux(IP,9)
     Tmp0(4,4) = Tmp0(4,4) + pinvq*Aux(IP,10)
     Tmp0(5,4) = Tmp0(5,4) + pinvq*Aux(IP,13)
     Tmp0(6,4) = Tmp0(6,4) + pinvq*Aux(IP,15)
     Tmp0(7,4) = Tmp0(7,4) + pinvq*Aux(IP,16)
     Tmp0(8,4) = Tmp0(8,4) + pinvq*Aux(IP,18)
     Tmp0(9,4) = Tmp0(9,4) + pinvq*Aux(IP,19)
     Tmp0(10,4) = Tmp0(10,4) + pinvq*Aux(IP,20)
     Tmp0(11,4) = Tmp0(11,4) + pinvq*Aux(IP,23)
     Tmp0(12,4) = Tmp0(12,4) + pinvq*Aux(IP,25)
     Tmp0(13,4) = Tmp0(13,4) + pinvq*Aux(IP,26)
     Tmp0(14,4) = Tmp0(14,4) + pinvq*Aux(IP,28)
     Tmp0(15,4) = Tmp0(15,4) + pinvq*Aux(IP,29)
     Tmp0(16,4) = Tmp0(16,4) + pinvq*Aux(IP,30)
     Tmp0(17,4) = Tmp0(17,4) + pinvq*Aux(IP,32)
     Tmp0(18,4) = Tmp0(18,4) + pinvq*Aux(IP,33)
     Tmp0(19,4) = Tmp0(19,4) + pinvq*Aux(IP,34)
     Tmp0(20,4) = Tmp0(20,4) + pinvq*Aux(IP,35)
!$ACC LOOP SEQ
     DO iTUVQ=1,  4
!$ACC LOOP SEQ
      DO iTUVP=1, 20
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP3Q1BtoDSeg1Prim
 subroutine TransferRecurrenceGPUP3Q2BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Aexp(1),Cexp(1)
  real(realk),intent(in) :: Aux(nPasses,   56)
  real(realk),intent(inout) :: Aux2(nPasses,   20,   10)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0( 20, 10)
  real(realk) :: Tmp1( 21: 35,  2:  4)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
  real(realk) :: invexpQ,inv2expQ,facX,facY,facZ,pinvq
  !CARTDIR = 1
  integer,parameter, dimension(35) :: TUVindexX1 = (/ 2,5,6,7,11,12,13,&
          & 14,15,16,21,22,23,24,25,26,27,28,29,30,36,37,38,39,&
          & 40,41,42,43,44,45,46,47,48,49,50 /)
  !CARTDIR = 2
  integer,parameter, dimension(35) :: TUVindexX2 = (/ 3,6,8,9,12,14,15,&
          & 17,18,19,22,24,25,27,28,29,31,32,33,34,37,39,40,42,&
          & 43,44,46,47,48,49,51,52,53,54,55 /)
  !CARTDIR = 3
  integer,parameter, dimension(35) :: TUVindexX3 = (/ 4,7,9,10,13,15,16,&
          & 18,19,20,23,25,26,28,29,30,32,33,34,35,38,40,41,43,&
          & 44,45,47,48,49,50,52,53,54,55,56 /)
  !CARTDIR = 1
  integer,parameter, dimension(20) :: IfacX1 = (/ 1,2,1,1,3,2,2,&
          & 1,1,1,4,3,3,2,2,2,1,1,1,1 /)
  !CARTDIR = 2
  integer,parameter, dimension(20) :: IfacX2 = (/ 1,1,2,1,1,2,1,&
          & 3,2,1,1,2,1,3,2,1,4,3,2,1 /)
  !CARTDIR = 3
  integer,parameter, dimension(20) :: IfacX3 = (/ 1,1,1,2,1,1,2,&
          & 1,2,3,1,1,2,1,2,3,1,2,3,4 /)
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expAX,expAY,expAZ,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Aexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux)
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
    expAX = -Aexp(1)*Xab
    expAY = -Aexp(1)*Yab
    expAZ = -Aexp(1)*Zab
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     facX = -(expAX-Cexp(1)*Xcd)*invexpQ
     facY = -(expAY-Cexp(1)*Ycd)*invexpQ
     facZ = -(expAZ-Cexp(1)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
!$ACC LOOP SEQ
     DO iTUVP=1, 20
      Tmp0(iTUVP,1) = Aux(IP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp1(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp1(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp1(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,20
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,35
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jq = 2
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(IP,iTUVP)
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
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(IP,iTUVP)
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
     DO iTUVQ=1, 10
!$ACC LOOP SEQ
      DO iTUVP=1, 20
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP3Q2BtoDSeg1Prim
 subroutine TransferRecurrenceGPUP3Q3BtoDSeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Aexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Aexp(1),Cexp(1)
  real(realk),intent(in) :: Aux(nPasses,   84)
  real(realk),intent(inout) :: Aux2(nPasses,   20,   20)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0( 20, 20)
  real(realk) :: Tmp1( 21: 56,  2:  4)
  real(realk) :: Tmp2( 21: 35,  5: 10)
!  real(realk) :: Tmp(nTUVP,nTUVQ) ordering
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvpminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
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
!$ACC         expAX,expAY,expAZ,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         Tmp2,&
!$ACC         invexpQ,inv2expQ,facX,facY,facZ,pinvq,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,&
!$ACC        reducedExponents,Pexp,Qexp,Pdistance12,Qdistance12,&
!$ACC       Aexp,Cexp,&
!$ACC        IatomApass,IatomBpass,Aux2,Aux)
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
    expAX = -Aexp(1)*Xab
    expAY = -Aexp(1)*Yab
    expAZ = -Aexp(1)*Zab
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     facX = -(expAX-Cexp(1)*Xcd)*invexpQ
     facY = -(expAY-Cexp(1)*Ycd)*invexpQ
     facZ = -(expAZ-Cexp(1)*Zcd)*invexpQ
     pinvq = -expP*invexpQ
 ! Building for Angular momentum Jq = 0
!$ACC LOOP SEQ
     DO iTUVP=1, 20
      Tmp0(iTUVP,1) = Aux(IP,iTUVP)
     ENDDO
 ! Building for Angular momentum Jq = 1
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 56
      Tmp1(iTUVP,2) = facX*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 56
      Tmp1(iTUVP,3) = facY*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 56
      Tmp1(iTUVP,4) = facZ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 1,10
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,35
      iTUVP = TUVindexX1(ituvpminus1)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + IfacX1(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,35
      iTUVP = TUVindexX2(ituvpminus1)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + IfacX2(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvpminus1 = 11,35
      iTUVP = TUVindexX3(ituvpminus1)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + IfacX3(ituvpminus1)*inv2expQ*Aux(IP,ituvpminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp0(iTUVP,2) = Tmp0(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,56
      iTUVplus1 = TUVindexX1(iTUVP)
      Tmp1(iTUVP,2) = Tmp1(iTUVP,2) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp0(iTUVP,3) = Tmp0(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,56
      iTUVplus1 = TUVindexX2(iTUVP)
      Tmp1(iTUVP,3) = Tmp1(iTUVP,3) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1,20
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp0(iTUVP,4) = Tmp0(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 21,56
      iTUVplus1 = TUVindexX3(iTUVP)
      Tmp1(iTUVP,4) = Tmp1(iTUVP,4) + pinvq*Aux(IP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jq = 2
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,5) = facX*Tmp0(iTUVP,2)+ inv2expQ*Aux(IP,iTUVP)
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
      Tmp0(iTUVP,8) = facY*Tmp0(iTUVP,3)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,9) = facY*Tmp0(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP = 1, 20
      Tmp0(iTUVP,10) = facZ*Tmp0(iTUVP,4)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp2(iTUVP,5) = facX*Tmp1(iTUVP,2)+ inv2expQ*Aux(IP,iTUVP)
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
      Tmp2(iTUVP,8) = facY*Tmp1(iTUVP,3)+ inv2expQ*Aux(IP,iTUVP)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp2(iTUVP,9) = facY*Tmp1(iTUVP,4)
     enddo
!$ACC LOOP SEQ
     do iTUVP =  21, 35
      Tmp2(iTUVP,10) = facZ*Tmp1(iTUVP,4)+ inv2expQ*Aux(IP,iTUVP)
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
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVP,iTUVQ)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP3Q3BtoDSeg1Prim
end module
