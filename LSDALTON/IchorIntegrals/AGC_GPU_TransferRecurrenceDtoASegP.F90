MODULE AGC_GPU_OBS_TRMODDtoASegP
 use IchorPrecisionMod
  
 CONTAINS
 subroutine TransferRecurrenceGPUP1Q2DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Cexp(nPrimC),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   20)
  real(realk),intent(inout) :: Aux2(nPrimQ*nPasses,    4,   10)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 10,  4)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimC,iPrimB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         iPrimC,iPrimB,&
!$ACC         Tmp0,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        Bexp,Cexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
!$ACC LOOP SEQ
   DO iTUVQ=1, 10
!$ACC LOOP SEQ
    DO iTUVP=1,  4
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpP
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpP
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
!$ACC LOOP SEQ
     DO iTUVQ=1, 10
      Tmp0(iTUVQ,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
!$ACC LOOP SEQ
     do iTUVQ = 1, 10
      Tmp0(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 10
      Tmp0(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 10
      Tmp0(iTUVQ,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(1,2) = Tmp0(1,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,2)
     Tmp0(2,2) = Tmp0(2,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,5)
     Tmp0(3,2) = Tmp0(3,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(4,2) = Tmp0(4,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(5,2) = Tmp0(5,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,11)
     Tmp0(6,2) = Tmp0(6,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(7,2) = Tmp0(7,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(8,2) = Tmp0(8,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(9,2) = Tmp0(9,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(10,2) = Tmp0(10,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp0(1,3) = Tmp0(1,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,3)
     Tmp0(2,3) = Tmp0(2,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(3,3) = Tmp0(3,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,8)
     Tmp0(4,3) = Tmp0(4,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(5,3) = Tmp0(5,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(6,3) = Tmp0(6,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(7,3) = Tmp0(7,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(8,3) = Tmp0(8,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,17)
     Tmp0(9,3) = Tmp0(9,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(10,3) = Tmp0(10,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp0(1,4) = Tmp0(1,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,4)
     Tmp0(2,4) = Tmp0(2,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(3,4) = Tmp0(3,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(4,4) = Tmp0(4,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,10)
     Tmp0(5,4) = Tmp0(5,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(6,4) = Tmp0(6,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(7,4) = Tmp0(7,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp0(8,4) = Tmp0(8,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(9,4) = Tmp0(9,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp0(10,4) = Tmp0(10,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,20)
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
!$ACC LOOP SEQ
     DO iTUVQ=1, 10
!$ACC LOOP SEQ
      DO iTUVP=1,  4
        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
 end subroutine TransferRecurrenceGPUP1Q2DtoASegP
 subroutine TransferRecurrenceGPUP1Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Cexp(nPrimC),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   35)
  real(realk),intent(inout) :: Aux2(nPrimQ*nPasses,    4,   20)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 20,  4)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimC,iPrimB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         iPrimC,iPrimB,&
!$ACC         Tmp0,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        Bexp,Cexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
!$ACC LOOP SEQ
   DO iTUVQ=1, 20
!$ACC LOOP SEQ
    DO iTUVP=1,  4
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpP
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpP
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
!$ACC LOOP SEQ
     DO iTUVQ=1, 20
      Tmp0(iTUVQ,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(11,2) = Tmp0(11,2) + 3*inv2expP*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp0(12,2) = Tmp0(12,2) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp0(13,2) = Tmp0(13,2) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp0(14,2) = Tmp0(14,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp0(15,2) = Tmp0(15,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp0(16,2) = Tmp0(16,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(12,3) = Tmp0(12,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp0(14,3) = Tmp0(14,3) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp0(15,3) = Tmp0(15,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp0(17,3) = Tmp0(17,3) + 3*inv2expP*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp0(18,3) = Tmp0(18,3) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp0(19,3) = Tmp0(19,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,1) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,2) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,3) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,4) 
     Tmp0(13,4) = Tmp0(13,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,5) 
     Tmp0(15,4) = Tmp0(15,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,6) 
     Tmp0(16,4) = Tmp0(16,4) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,7) 
     Tmp0(18,4) = Tmp0(18,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,8) 
     Tmp0(19,4) = Tmp0(19,4) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,9) 
     Tmp0(20,4) = Tmp0(20,4) + 3*inv2expP*Aux(iPrimQ,iPrimP,iPassP,10) 
     Tmp0(1,2) = Tmp0(1,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,2)
     Tmp0(2,2) = Tmp0(2,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,5)
     Tmp0(3,2) = Tmp0(3,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(4,2) = Tmp0(4,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(5,2) = Tmp0(5,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,11)
     Tmp0(6,2) = Tmp0(6,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(7,2) = Tmp0(7,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(8,2) = Tmp0(8,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(9,2) = Tmp0(9,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(10,2) = Tmp0(10,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp0(11,2) = Tmp0(11,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,21)
     Tmp0(12,2) = Tmp0(12,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,22)
     Tmp0(13,2) = Tmp0(13,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,23)
     Tmp0(14,2) = Tmp0(14,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,24)
     Tmp0(15,2) = Tmp0(15,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp0(16,2) = Tmp0(16,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,26)
     Tmp0(17,2) = Tmp0(17,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,27)
     Tmp0(18,2) = Tmp0(18,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp0(19,2) = Tmp0(19,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp0(20,2) = Tmp0(20,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,30)
     Tmp0(1,3) = Tmp0(1,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,3)
     Tmp0(2,3) = Tmp0(2,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,6)
     Tmp0(3,3) = Tmp0(3,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,8)
     Tmp0(4,3) = Tmp0(4,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(5,3) = Tmp0(5,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,12)
     Tmp0(6,3) = Tmp0(6,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,14)
     Tmp0(7,3) = Tmp0(7,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(8,3) = Tmp0(8,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,17)
     Tmp0(9,3) = Tmp0(9,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(10,3) = Tmp0(10,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp0(11,3) = Tmp0(11,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,22)
     Tmp0(12,3) = Tmp0(12,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,24)
     Tmp0(13,3) = Tmp0(13,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp0(14,3) = Tmp0(14,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,27)
     Tmp0(15,3) = Tmp0(15,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp0(16,3) = Tmp0(16,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp0(17,3) = Tmp0(17,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,31)
     Tmp0(18,3) = Tmp0(18,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,32)
     Tmp0(19,3) = Tmp0(19,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,33)
     Tmp0(20,3) = Tmp0(20,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,34)
     Tmp0(1,4) = Tmp0(1,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,4)
     Tmp0(2,4) = Tmp0(2,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,7)
     Tmp0(3,4) = Tmp0(3,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,9)
     Tmp0(4,4) = Tmp0(4,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,10)
     Tmp0(5,4) = Tmp0(5,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,13)
     Tmp0(6,4) = Tmp0(6,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,15)
     Tmp0(7,4) = Tmp0(7,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,16)
     Tmp0(8,4) = Tmp0(8,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,18)
     Tmp0(9,4) = Tmp0(9,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,19)
     Tmp0(10,4) = Tmp0(10,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,20)
     Tmp0(11,4) = Tmp0(11,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,23)
     Tmp0(12,4) = Tmp0(12,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,25)
     Tmp0(13,4) = Tmp0(13,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,26)
     Tmp0(14,4) = Tmp0(14,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,28)
     Tmp0(15,4) = Tmp0(15,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,29)
     Tmp0(16,4) = Tmp0(16,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,30)
     Tmp0(17,4) = Tmp0(17,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,32)
     Tmp0(18,4) = Tmp0(18,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,33)
     Tmp0(19,4) = Tmp0(19,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,34)
     Tmp0(20,4) = Tmp0(20,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,35)
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
!$ACC LOOP SEQ
     DO iTUVQ=1, 20
!$ACC LOOP SEQ
      DO iTUVP=1,  4
        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
 end subroutine TransferRecurrenceGPUP1Q3DtoASegP
 subroutine TransferRecurrenceGPUP2Q3DtoASegP(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2,iASync)
  use AGC_GPU_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Cexp(nPrimC),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   56)
  real(realk),intent(inout) :: Aux2(nPrimQ*nPasses,   10,   20)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  integer(kind=acckind),intent(in) :: iASync
  !Local variables
  real(realk) :: Tmp0( 20, 10)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
  real(realk) :: Tmp1( 21: 35,  2:  4)
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimC,iPrimB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         iPrimC,iPrimB,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        TUVindexX1_35,TUVindexX2_35,TUVindexX3_35, &
!$ACC        IfacX1_20,IfacX2_20,IfacX3_20, &
!$ACC        Bexp,Cexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
!$ACC LOOP SEQ
   DO iTUVQ=1, 20
!$ACC LOOP SEQ
    DO iTUVP=1, 10
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpP
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpP
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
!$ACC LOOP SEQ
     DO iTUVQ=1, 20
      Tmp0(iTUVQ,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1_20(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2_20(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX3_35(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3_20(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + IfacX1_20(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + IfacX2_20(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX3_35(ituvqminus1)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + IfacX3_20(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jp = 2
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,5) = facX*Tmp0(iTUVQ,2)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
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
      Tmp0(iTUVQ,8) = facY*Tmp0(iTUVQ,3)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,9) = facY*Tmp0(iTUVQ,4)
     enddo
!$ACC LOOP SEQ
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,10) = facZ*Tmp0(iTUVQ,4)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
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
        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPasses
 end subroutine TransferRecurrenceGPUP2Q3DtoASegP
end module
