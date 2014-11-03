MODULE AGC_GPU_OBS_TRMODCtoASeg1Prim
 use IchorPrecisionMod
  
 CONTAINS
 subroutine TransferRecurrenceGPUP1Q2CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(1),Bexp(1)
  real(realk),intent(in) :: Aux(nPasses,   20)
  real(realk),intent(inout) :: Aux2(nPasses,    4,   10)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 10,  4)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         Tmp0,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        Bexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
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
    expP = Pexp(1)
    invexpP = D1/Pexp(1)
    inv2expP = D05*invexpP
    expBX = Bexp(1)*Xab
    expBY = Bexp(1)*Yab
    expBZ = Bexp(1)*Zab
     facX = -(expBX+Dexp(1)*Xcd)*invexpP
     facY = -(expBY+Dexp(1)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(1)*Zcd)*invexpP
     qinvp = -Qexp(1)*invexpP
 ! Building for Angular momentum Jp = 0
!$ACC LOOP SEQ
     DO iTUVQ=1, 10
      Tmp0(iTUVQ,1) = Aux(iP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
!$ACC LOOP SEQ
     do iTUVQ = 1, 10
      Tmp0(iTUVQ,2) = facX*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 10
      Tmp0(iTUVQ,3) = facY*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 10
      Tmp0(iTUVQ,4) = facZ*Aux(iP,iTUVQ)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expP*Aux(iP,1) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expP*Aux(iP,2) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expP*Aux(iP,3) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expP*Aux(iP,4) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expP*Aux(iP,1) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expP*Aux(iP,2) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expP*Aux(iP,3) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expP*Aux(iP,4) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expP*Aux(iP,1) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expP*Aux(iP,2) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expP*Aux(iP,3) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expP*Aux(iP,4) 
     Tmp0(1,2) = Tmp0(1,2) + qinvp*Aux(iP,2)
     Tmp0(2,2) = Tmp0(2,2) + qinvp*Aux(iP,5)
     Tmp0(3,2) = Tmp0(3,2) + qinvp*Aux(iP,6)
     Tmp0(4,2) = Tmp0(4,2) + qinvp*Aux(iP,7)
     Tmp0(5,2) = Tmp0(5,2) + qinvp*Aux(iP,11)
     Tmp0(6,2) = Tmp0(6,2) + qinvp*Aux(iP,12)
     Tmp0(7,2) = Tmp0(7,2) + qinvp*Aux(iP,13)
     Tmp0(8,2) = Tmp0(8,2) + qinvp*Aux(iP,14)
     Tmp0(9,2) = Tmp0(9,2) + qinvp*Aux(iP,15)
     Tmp0(10,2) = Tmp0(10,2) + qinvp*Aux(iP,16)
     Tmp0(1,3) = Tmp0(1,3) + qinvp*Aux(iP,3)
     Tmp0(2,3) = Tmp0(2,3) + qinvp*Aux(iP,6)
     Tmp0(3,3) = Tmp0(3,3) + qinvp*Aux(iP,8)
     Tmp0(4,3) = Tmp0(4,3) + qinvp*Aux(iP,9)
     Tmp0(5,3) = Tmp0(5,3) + qinvp*Aux(iP,12)
     Tmp0(6,3) = Tmp0(6,3) + qinvp*Aux(iP,14)
     Tmp0(7,3) = Tmp0(7,3) + qinvp*Aux(iP,15)
     Tmp0(8,3) = Tmp0(8,3) + qinvp*Aux(iP,17)
     Tmp0(9,3) = Tmp0(9,3) + qinvp*Aux(iP,18)
     Tmp0(10,3) = Tmp0(10,3) + qinvp*Aux(iP,19)
     Tmp0(1,4) = Tmp0(1,4) + qinvp*Aux(iP,4)
     Tmp0(2,4) = Tmp0(2,4) + qinvp*Aux(iP,7)
     Tmp0(3,4) = Tmp0(3,4) + qinvp*Aux(iP,9)
     Tmp0(4,4) = Tmp0(4,4) + qinvp*Aux(iP,10)
     Tmp0(5,4) = Tmp0(5,4) + qinvp*Aux(iP,13)
     Tmp0(6,4) = Tmp0(6,4) + qinvp*Aux(iP,15)
     Tmp0(7,4) = Tmp0(7,4) + qinvp*Aux(iP,16)
     Tmp0(8,4) = Tmp0(8,4) + qinvp*Aux(iP,18)
     Tmp0(9,4) = Tmp0(9,4) + qinvp*Aux(iP,19)
     Tmp0(10,4) = Tmp0(10,4) + qinvp*Aux(iP,20)
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
!$ACC LOOP SEQ
     DO iTUVQ=1, 10
!$ACC LOOP SEQ
      DO iTUVP=1,  4
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP1Q2CtoASeg1Prim
 subroutine TransferRecurrenceGPUP1Q3CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(1),Bexp(1)
  real(realk),intent(in) :: Aux(nPasses,   35)
  real(realk),intent(inout) :: Aux2(nPasses,    4,   20)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 20,  4)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         Tmp0,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        Bexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
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
    expP = Pexp(1)
    invexpP = D1/Pexp(1)
    inv2expP = D05*invexpP
    expBX = Bexp(1)*Xab
    expBY = Bexp(1)*Yab
    expBZ = Bexp(1)*Zab
     facX = -(expBX+Dexp(1)*Xcd)*invexpP
     facY = -(expBY+Dexp(1)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(1)*Zcd)*invexpP
     qinvp = -Qexp(1)*invexpP
 ! Building for Angular momentum Jp = 0
!$ACC LOOP SEQ
     DO iTUVQ=1, 20
      Tmp0(iTUVQ,1) = Aux(iP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,2) = facX*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,3) = facY*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,4) = facZ*Aux(iP,iTUVQ)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expP*Aux(iP,1) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expP*Aux(iP,2) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expP*Aux(iP,3) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expP*Aux(iP,4) 
     Tmp0(11,2) = Tmp0(11,2) + 3*inv2expP*Aux(iP,5) 
     Tmp0(12,2) = Tmp0(12,2) + 2*inv2expP*Aux(iP,6) 
     Tmp0(13,2) = Tmp0(13,2) + 2*inv2expP*Aux(iP,7) 
     Tmp0(14,2) = Tmp0(14,2) + inv2expP*Aux(iP,8) 
     Tmp0(15,2) = Tmp0(15,2) + inv2expP*Aux(iP,9) 
     Tmp0(16,2) = Tmp0(16,2) + inv2expP*Aux(iP,10) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expP*Aux(iP,1) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expP*Aux(iP,2) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expP*Aux(iP,3) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expP*Aux(iP,4) 
     Tmp0(12,3) = Tmp0(12,3) + inv2expP*Aux(iP,5) 
     Tmp0(14,3) = Tmp0(14,3) + 2*inv2expP*Aux(iP,6) 
     Tmp0(15,3) = Tmp0(15,3) + inv2expP*Aux(iP,7) 
     Tmp0(17,3) = Tmp0(17,3) + 3*inv2expP*Aux(iP,8) 
     Tmp0(18,3) = Tmp0(18,3) + 2*inv2expP*Aux(iP,9) 
     Tmp0(19,3) = Tmp0(19,3) + inv2expP*Aux(iP,10) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expP*Aux(iP,1) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expP*Aux(iP,2) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expP*Aux(iP,3) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expP*Aux(iP,4) 
     Tmp0(13,4) = Tmp0(13,4) + inv2expP*Aux(iP,5) 
     Tmp0(15,4) = Tmp0(15,4) + inv2expP*Aux(iP,6) 
     Tmp0(16,4) = Tmp0(16,4) + 2*inv2expP*Aux(iP,7) 
     Tmp0(18,4) = Tmp0(18,4) + inv2expP*Aux(iP,8) 
     Tmp0(19,4) = Tmp0(19,4) + 2*inv2expP*Aux(iP,9) 
     Tmp0(20,4) = Tmp0(20,4) + 3*inv2expP*Aux(iP,10) 
     Tmp0(1,2) = Tmp0(1,2) + qinvp*Aux(iP,2)
     Tmp0(2,2) = Tmp0(2,2) + qinvp*Aux(iP,5)
     Tmp0(3,2) = Tmp0(3,2) + qinvp*Aux(iP,6)
     Tmp0(4,2) = Tmp0(4,2) + qinvp*Aux(iP,7)
     Tmp0(5,2) = Tmp0(5,2) + qinvp*Aux(iP,11)
     Tmp0(6,2) = Tmp0(6,2) + qinvp*Aux(iP,12)
     Tmp0(7,2) = Tmp0(7,2) + qinvp*Aux(iP,13)
     Tmp0(8,2) = Tmp0(8,2) + qinvp*Aux(iP,14)
     Tmp0(9,2) = Tmp0(9,2) + qinvp*Aux(iP,15)
     Tmp0(10,2) = Tmp0(10,2) + qinvp*Aux(iP,16)
     Tmp0(11,2) = Tmp0(11,2) + qinvp*Aux(iP,21)
     Tmp0(12,2) = Tmp0(12,2) + qinvp*Aux(iP,22)
     Tmp0(13,2) = Tmp0(13,2) + qinvp*Aux(iP,23)
     Tmp0(14,2) = Tmp0(14,2) + qinvp*Aux(iP,24)
     Tmp0(15,2) = Tmp0(15,2) + qinvp*Aux(iP,25)
     Tmp0(16,2) = Tmp0(16,2) + qinvp*Aux(iP,26)
     Tmp0(17,2) = Tmp0(17,2) + qinvp*Aux(iP,27)
     Tmp0(18,2) = Tmp0(18,2) + qinvp*Aux(iP,28)
     Tmp0(19,2) = Tmp0(19,2) + qinvp*Aux(iP,29)
     Tmp0(20,2) = Tmp0(20,2) + qinvp*Aux(iP,30)
     Tmp0(1,3) = Tmp0(1,3) + qinvp*Aux(iP,3)
     Tmp0(2,3) = Tmp0(2,3) + qinvp*Aux(iP,6)
     Tmp0(3,3) = Tmp0(3,3) + qinvp*Aux(iP,8)
     Tmp0(4,3) = Tmp0(4,3) + qinvp*Aux(iP,9)
     Tmp0(5,3) = Tmp0(5,3) + qinvp*Aux(iP,12)
     Tmp0(6,3) = Tmp0(6,3) + qinvp*Aux(iP,14)
     Tmp0(7,3) = Tmp0(7,3) + qinvp*Aux(iP,15)
     Tmp0(8,3) = Tmp0(8,3) + qinvp*Aux(iP,17)
     Tmp0(9,3) = Tmp0(9,3) + qinvp*Aux(iP,18)
     Tmp0(10,3) = Tmp0(10,3) + qinvp*Aux(iP,19)
     Tmp0(11,3) = Tmp0(11,3) + qinvp*Aux(iP,22)
     Tmp0(12,3) = Tmp0(12,3) + qinvp*Aux(iP,24)
     Tmp0(13,3) = Tmp0(13,3) + qinvp*Aux(iP,25)
     Tmp0(14,3) = Tmp0(14,3) + qinvp*Aux(iP,27)
     Tmp0(15,3) = Tmp0(15,3) + qinvp*Aux(iP,28)
     Tmp0(16,3) = Tmp0(16,3) + qinvp*Aux(iP,29)
     Tmp0(17,3) = Tmp0(17,3) + qinvp*Aux(iP,31)
     Tmp0(18,3) = Tmp0(18,3) + qinvp*Aux(iP,32)
     Tmp0(19,3) = Tmp0(19,3) + qinvp*Aux(iP,33)
     Tmp0(20,3) = Tmp0(20,3) + qinvp*Aux(iP,34)
     Tmp0(1,4) = Tmp0(1,4) + qinvp*Aux(iP,4)
     Tmp0(2,4) = Tmp0(2,4) + qinvp*Aux(iP,7)
     Tmp0(3,4) = Tmp0(3,4) + qinvp*Aux(iP,9)
     Tmp0(4,4) = Tmp0(4,4) + qinvp*Aux(iP,10)
     Tmp0(5,4) = Tmp0(5,4) + qinvp*Aux(iP,13)
     Tmp0(6,4) = Tmp0(6,4) + qinvp*Aux(iP,15)
     Tmp0(7,4) = Tmp0(7,4) + qinvp*Aux(iP,16)
     Tmp0(8,4) = Tmp0(8,4) + qinvp*Aux(iP,18)
     Tmp0(9,4) = Tmp0(9,4) + qinvp*Aux(iP,19)
     Tmp0(10,4) = Tmp0(10,4) + qinvp*Aux(iP,20)
     Tmp0(11,4) = Tmp0(11,4) + qinvp*Aux(iP,23)
     Tmp0(12,4) = Tmp0(12,4) + qinvp*Aux(iP,25)
     Tmp0(13,4) = Tmp0(13,4) + qinvp*Aux(iP,26)
     Tmp0(14,4) = Tmp0(14,4) + qinvp*Aux(iP,28)
     Tmp0(15,4) = Tmp0(15,4) + qinvp*Aux(iP,29)
     Tmp0(16,4) = Tmp0(16,4) + qinvp*Aux(iP,30)
     Tmp0(17,4) = Tmp0(17,4) + qinvp*Aux(iP,32)
     Tmp0(18,4) = Tmp0(18,4) + qinvp*Aux(iP,33)
     Tmp0(19,4) = Tmp0(19,4) + qinvp*Aux(iP,34)
     Tmp0(20,4) = Tmp0(20,4) + qinvp*Aux(iP,35)
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
!$ACC LOOP SEQ
     DO iTUVQ=1, 20
!$ACC LOOP SEQ
      DO iTUVP=1,  4
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP1Q3CtoASeg1Prim
 subroutine TransferRecurrenceGPUP1Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  use AGC_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(1),Bexp(1)
  real(realk),intent(in) :: Aux(nPasses,   56)
  real(realk),intent(inout) :: Aux2(nPasses,    4,   35)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 35,  4)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         Tmp0,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        TUVindexX1_35,TUVindexX2_35,TUVindexX3_35, &
!$ACC        IfacX1_20,IfacX2_20,IfacX3_20, &
!$ACC        Bexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
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
    expP = Pexp(1)
    invexpP = D1/Pexp(1)
    inv2expP = D05*invexpP
    expBX = Bexp(1)*Xab
    expBY = Bexp(1)*Yab
    expBZ = Bexp(1)*Zab
     facX = -(expBX+Dexp(1)*Xcd)*invexpP
     facY = -(expBY+Dexp(1)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(1)*Zcd)*invexpP
     qinvp = -Qexp(1)*invexpP
 ! Building for Angular momentum Jp = 0
!$ACC LOOP SEQ
     DO iTUVQ=1, 35
      Tmp0(iTUVQ,1) = Aux(iP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,2) = facX*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,3) = facY*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,4) = facZ*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1_20(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2_20(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_35(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3_20(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iP,iTUVplus1)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
!$ACC LOOP SEQ
     DO iTUVQ=1, 35
!$ACC LOOP SEQ
      DO iTUVP=1,  4
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP1Q4CtoASeg1Prim
 subroutine TransferRecurrenceGPUP2Q3CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  use AGC_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(1),Bexp(1)
  real(realk),intent(in) :: Aux(nPasses,   56)
  real(realk),intent(inout) :: Aux2(nPasses,   10,   20)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 20, 10)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
  real(realk) :: Tmp1( 21: 35,  2:  4)
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        TUVindexX1_35,TUVindexX2_35,TUVindexX3_35, &
!$ACC        IfacX1_20,IfacX2_20,IfacX3_20, &
!$ACC        Bexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
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
    expP = Pexp(1)
    invexpP = D1/Pexp(1)
    inv2expP = D05*invexpP
    expBX = Bexp(1)*Xab
    expBY = Bexp(1)*Yab
    expBZ = Bexp(1)*Zab
     facX = -(expBX+Dexp(1)*Xcd)*invexpP
     facY = -(expBY+Dexp(1)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(1)*Zcd)*invexpP
     qinvp = -Qexp(1)*invexpP
 ! Building for Angular momentum Jp = 0
!$ACC LOOP SEQ
     DO iTUVQ=1, 20
      Tmp0(iTUVQ,1) = Aux(iP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,2) = facX*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,3) = facY*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,4) = facZ*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,2) = facX*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,3) = facY*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,4) = facZ*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1_20(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2_20(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX3_35(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3_20(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + IfacX1_20(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + IfacX2_20(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX3_35(ituvqminus1)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + IfacX3_20(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + qinvp*Aux(iP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jp = 2
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,5) = facX*Tmp0(iTUVQ,2)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,6) = facX*Tmp0(iTUVQ,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,7) = facX*Tmp0(iTUVQ,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,8) = facY*Tmp0(iTUVQ,3)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,9) = facY*Tmp0(iTUVQ,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,10) = facZ*Tmp0(iTUVQ,4)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + IfacX1_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + IfacX1_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + IfacX1_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + IfacX2_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + IfacX2_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX3_35(ituvqminus1)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + IfacX3_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp0(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
!$ACC LOOP SEQ
     DO iTUVQ=1, 20
!$ACC LOOP SEQ
      DO iTUVP=1, 10
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP2Q3CtoASeg1Prim
 subroutine TransferRecurrenceGPUP2Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  use AGC_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(1),Bexp(1)
  real(realk),intent(in) :: Aux(nPasses,   84)
  real(realk),intent(inout) :: Aux2(nPasses,   10,   35)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 35, 10)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
  real(realk) :: Tmp1( 36: 56,  2:  4)
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        TUVindexX1_56,TUVindexX2_56,TUVindexX3_56, &
!$ACC        IfacX1_35,IfacX2_35,IfacX3_35, &
!$ACC        Bexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
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
    expP = Pexp(1)
    invexpP = D1/Pexp(1)
    inv2expP = D05*invexpP
    expBX = Bexp(1)*Xab
    expBY = Bexp(1)*Yab
    expBZ = Bexp(1)*Zab
     facX = -(expBX+Dexp(1)*Xcd)*invexpP
     facY = -(expBY+Dexp(1)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(1)*Zcd)*invexpP
     qinvp = -Qexp(1)*invexpP
 ! Building for Angular momentum Jp = 0
!$ACC LOOP SEQ
     DO iTUVQ=1, 35
      Tmp0(iTUVQ,1) = Aux(iP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,2) = facX*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,3) = facY*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,4) = facZ*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 56
      Tmp1(iTUVQ,2) = facX*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 56
      Tmp1(iTUVQ,3) = facY*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 56
      Tmp1(iTUVQ,4) = facZ*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_56(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1_35(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_56(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2_35(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_56(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3_35(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1_56(ituvqminus1)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + IfacX1_35(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX2_56(ituvqminus1)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + IfacX2_35(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX3_56(ituvqminus1)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + IfacX3_35(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX3_56(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX3_56(iTUVQ)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + qinvp*Aux(iP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jp = 2
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,5) = facX*Tmp0(iTUVQ,2)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,6) = facX*Tmp0(iTUVQ,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,7) = facX*Tmp0(iTUVQ,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,8) = facY*Tmp0(iTUVQ,3)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,9) = facY*Tmp0(iTUVQ,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,10) = facZ*Tmp0(iTUVQ,4)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_56(ituvqminus1)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + IfacX1_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_56(ituvqminus1)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + IfacX1_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_56(ituvqminus1)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + IfacX1_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_56(ituvqminus1)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + IfacX2_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_56(ituvqminus1)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + IfacX2_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_56(ituvqminus1)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + IfacX3_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp0(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_56(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_56(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
!$ACC LOOP SEQ
     DO iTUVQ=1, 35
!$ACC LOOP SEQ
      DO iTUVP=1, 10
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP2Q4CtoASeg1Prim
 subroutine TransferRecurrenceGPUP3Q4CtoASeg1Prim(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  use AGC_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(1,1),Pexp(1),Qexp(1)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(1),Bexp(1)
  real(realk),intent(in) :: Aux(nPasses,  120)
  real(realk),intent(inout) :: Aux2(nPasses,   20,   35)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 35, 20)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
  real(realk) :: Tmp1( 36: 84,  2:  4)
  real(realk) :: Tmp2( 36: 56,  5: 10)
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         Tmp2,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        TUVindexX1_84,TUVindexX2_84,TUVindexX3_84, &
!$ACC        IfacX1_56,IfacX2_56,IfacX3_56, &
!$ACC        Bexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
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
    expP = Pexp(1)
    invexpP = D1/Pexp(1)
    inv2expP = D05*invexpP
    expBX = Bexp(1)*Xab
    expBY = Bexp(1)*Yab
    expBZ = Bexp(1)*Zab
     facX = -(expBX+Dexp(1)*Xcd)*invexpP
     facY = -(expBY+Dexp(1)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(1)*Zcd)*invexpP
     qinvp = -Qexp(1)*invexpP
 ! Building for Angular momentum Jp = 0
!$ACC LOOP SEQ
     DO iTUVQ=1, 35
      Tmp0(iTUVQ,1) = Aux(iP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,2) = facX*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,3) = facY*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,4) = facZ*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 84
      Tmp1(iTUVQ,2) = facX*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 84
      Tmp1(iTUVQ,3) = facY*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 84
      Tmp1(iTUVQ,4) = facZ*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1_56(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2_56(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3_56(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,56
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + IfacX1_56(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,56
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + IfacX2_56(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,56
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + IfacX3_56(ituvqminus1)*inv2expP*Aux(iP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,84
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,84
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,84
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + qinvp*Aux(iP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jp = 2
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,5) = facX*Tmp0(iTUVQ,2)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,6) = facX*Tmp0(iTUVQ,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,7) = facX*Tmp0(iTUVQ,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,8) = facY*Tmp0(iTUVQ,3)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,9) = facY*Tmp0(iTUVQ,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,10) = facZ*Tmp0(iTUVQ,4)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,5) = facX*Tmp1(iTUVQ,2)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,6) = facX*Tmp1(iTUVQ,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,7) = facX*Tmp1(iTUVQ,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,8) = facY*Tmp1(iTUVQ,3)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,9) = facY*Tmp1(iTUVQ,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,10) = facZ*Tmp1(iTUVQ,4)+ inv2expP*Aux(iP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + IfacX3_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp2(iTUVQ,5) = Tmp2(iTUVQ,5) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp2(iTUVQ,6) = Tmp2(iTUVQ,6) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp2(iTUVQ,7) = Tmp2(iTUVQ,7) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp2(iTUVQ,8) = Tmp2(iTUVQ,8) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp2(iTUVQ,9) = Tmp2(iTUVQ,9) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp2(iTUVQ,10) = Tmp2(iTUVQ,10) + IfacX3_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp0(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp2(iTUVQ,5) = Tmp2(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp2(iTUVQ,6) = Tmp2(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp2(iTUVQ,7) = Tmp2(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp0(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp2(iTUVQ,8) = Tmp2(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp2(iTUVQ,9) = Tmp2(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp0(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp2(iTUVQ,10) = Tmp2(iTUVQ,10) + qinvp*Tmp1(iTUVplus1,4)
     enddo
 ! Building for Angular momentum Jp = 3
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,11) = facX*Tmp0(iTUVQ,5)+2*inv2expP*Tmp0(iTUVQ,2)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,12) = facY*Tmp0(iTUVQ,5)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,13) = facZ*Tmp0(iTUVQ,5)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,14) = facX*Tmp0(iTUVQ,8)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,15) = facX*Tmp0(iTUVQ,9)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,16) = facX*Tmp0(iTUVQ,10)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,17) = facY*Tmp0(iTUVQ,8)+2*inv2expP*Tmp0(iTUVQ,3)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,18) = facZ*Tmp0(iTUVQ,8)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,19) = facY*Tmp0(iTUVQ,10)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,20) = facZ*Tmp0(iTUVQ,10)+2*inv2expP*Tmp0(iTUVQ,4)
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,11) = Tmp0(iTUVQ,11) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,12) = Tmp0(iTUVQ,12) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp0(iTUVQ,13) = Tmp0(iTUVQ,13) + IfacX3_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,5) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,14) = Tmp0(iTUVQ,14) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,15) = Tmp0(iTUVQ,15) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,9) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,16) = Tmp0(iTUVQ,16) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,17) = Tmp0(iTUVQ,17) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp0(iTUVQ,18) = Tmp0(iTUVQ,18) + IfacX3_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,8) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,19) = Tmp0(iTUVQ,19) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,10) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp0(iTUVQ,20) = Tmp0(iTUVQ,20) + IfacX3_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,10) 
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,11) = Tmp0(iTUVQ,11) + qinvp*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,11) = Tmp0(iTUVQ,11) + qinvp*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,12) = Tmp0(iTUVQ,12) + qinvp*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,12) = Tmp0(iTUVQ,12) + qinvp*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,13) = Tmp0(iTUVQ,13) + qinvp*Tmp0(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,13) = Tmp0(iTUVQ,13) + qinvp*Tmp2(iTUVplus1,5)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,14) = Tmp0(iTUVQ,14) + qinvp*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,14) = Tmp0(iTUVQ,14) + qinvp*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,15) = Tmp0(iTUVQ,15) + qinvp*Tmp0(iTUVplus1,9)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,15) = Tmp0(iTUVQ,15) + qinvp*Tmp2(iTUVplus1,9)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,16) = Tmp0(iTUVQ,16) + qinvp*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,16) = Tmp0(iTUVQ,16) + qinvp*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,17) = Tmp0(iTUVQ,17) + qinvp*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,17) = Tmp0(iTUVQ,17) + qinvp*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,18) = Tmp0(iTUVQ,18) + qinvp*Tmp0(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,18) = Tmp0(iTUVQ,18) + qinvp*Tmp2(iTUVplus1,8)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,19) = Tmp0(iTUVQ,19) + qinvp*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,19) = Tmp0(iTUVQ,19) + qinvp*Tmp2(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,20) = Tmp0(iTUVQ,20) + qinvp*Tmp0(iTUVplus1,10)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,20) = Tmp0(iTUVQ,20) + qinvp*Tmp2(iTUVplus1,10)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
!$ACC LOOP SEQ
     DO iTUVQ=1, 35
!$ACC LOOP SEQ
      DO iTUVP=1, 20
        Aux2(IP,iTUVP,iTUVQ) = Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP3Q4CtoASeg1Prim
end module
