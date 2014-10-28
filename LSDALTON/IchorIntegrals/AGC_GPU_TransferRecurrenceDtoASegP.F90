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
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,Aux2) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1, 10
    DO iTUVP=1,  4
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
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
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,Aux2) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1, 20
    DO iTUVP=1,  4
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
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
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPrimQ,nPasses,Aux2) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPasses
   DO iTUVQ=1, 20
    DO iTUVP=1, 10
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPassP,&
!$ACC         expBX,expBY,expBZ,&
!$ACC         iPrimC,iPrimB,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        Bexp,Cexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux) ASYNC(iASync)
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
     Tmp1(21,2) = Tmp1(21,2) + 4*inv2expP*Aux(iPrimQ,iPrimP,iPassP,11) 
     Tmp1(22,2) = Tmp1(22,2) + 3*inv2expP*Aux(iPrimQ,iPrimP,iPassP,12) 
     Tmp1(23,2) = Tmp1(23,2) + 3*inv2expP*Aux(iPrimQ,iPrimP,iPassP,13) 
     Tmp1(24,2) = Tmp1(24,2) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,14) 
     Tmp1(25,2) = Tmp1(25,2) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,15) 
     Tmp1(26,2) = Tmp1(26,2) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,16) 
     Tmp1(27,2) = Tmp1(27,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,17) 
     Tmp1(28,2) = Tmp1(28,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,18) 
     Tmp1(29,2) = Tmp1(29,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,19) 
     Tmp1(30,2) = Tmp1(30,2) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,20) 
     Tmp1(22,3) = Tmp1(22,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,11) 
     Tmp1(24,3) = Tmp1(24,3) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,12) 
     Tmp1(25,3) = Tmp1(25,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,13) 
     Tmp1(27,3) = Tmp1(27,3) + 3*inv2expP*Aux(iPrimQ,iPrimP,iPassP,14) 
     Tmp1(28,3) = Tmp1(28,3) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,15) 
     Tmp1(29,3) = Tmp1(29,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,16) 
     Tmp1(31,3) = Tmp1(31,3) + 4*inv2expP*Aux(iPrimQ,iPrimP,iPassP,17) 
     Tmp1(32,3) = Tmp1(32,3) + 3*inv2expP*Aux(iPrimQ,iPrimP,iPassP,18) 
     Tmp1(33,3) = Tmp1(33,3) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,19) 
     Tmp1(34,3) = Tmp1(34,3) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,20) 
     Tmp1(23,4) = Tmp1(23,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,11) 
     Tmp1(25,4) = Tmp1(25,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,12) 
     Tmp1(26,4) = Tmp1(26,4) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,13) 
     Tmp1(28,4) = Tmp1(28,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,14) 
     Tmp1(29,4) = Tmp1(29,4) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,15) 
     Tmp1(30,4) = Tmp1(30,4) + 3*inv2expP*Aux(iPrimQ,iPrimP,iPassP,16) 
     Tmp1(32,4) = Tmp1(32,4) + inv2expP*Aux(iPrimQ,iPrimP,iPassP,17) 
     Tmp1(33,4) = Tmp1(33,4) + 2*inv2expP*Aux(iPrimQ,iPrimP,iPassP,18) 
     Tmp1(34,4) = Tmp1(34,4) + 3*inv2expP*Aux(iPrimQ,iPrimP,iPassP,19) 
     Tmp1(35,4) = Tmp1(35,4) + 4*inv2expP*Aux(iPrimQ,iPrimP,iPassP,20) 
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
     Tmp1(21,2) = Tmp1(21,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,36)
     Tmp1(22,2) = Tmp1(22,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,37)
     Tmp1(23,2) = Tmp1(23,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,38)
     Tmp1(24,2) = Tmp1(24,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,39)
     Tmp1(25,2) = Tmp1(25,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,40)
     Tmp1(26,2) = Tmp1(26,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,41)
     Tmp1(27,2) = Tmp1(27,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,42)
     Tmp1(28,2) = Tmp1(28,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,43)
     Tmp1(29,2) = Tmp1(29,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,44)
     Tmp1(30,2) = Tmp1(30,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,45)
     Tmp1(31,2) = Tmp1(31,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,46)
     Tmp1(32,2) = Tmp1(32,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,47)
     Tmp1(33,2) = Tmp1(33,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,48)
     Tmp1(34,2) = Tmp1(34,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,49)
     Tmp1(35,2) = Tmp1(35,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,50)
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
     Tmp1(21,3) = Tmp1(21,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,37)
     Tmp1(22,3) = Tmp1(22,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,39)
     Tmp1(23,3) = Tmp1(23,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,40)
     Tmp1(24,3) = Tmp1(24,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,42)
     Tmp1(25,3) = Tmp1(25,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,43)
     Tmp1(26,3) = Tmp1(26,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,44)
     Tmp1(27,3) = Tmp1(27,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,46)
     Tmp1(28,3) = Tmp1(28,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,47)
     Tmp1(29,3) = Tmp1(29,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,48)
     Tmp1(30,3) = Tmp1(30,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,49)
     Tmp1(31,3) = Tmp1(31,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,51)
     Tmp1(32,3) = Tmp1(32,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,52)
     Tmp1(33,3) = Tmp1(33,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,53)
     Tmp1(34,3) = Tmp1(34,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,54)
     Tmp1(35,3) = Tmp1(35,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,55)
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
     Tmp1(21,4) = Tmp1(21,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,38)
     Tmp1(22,4) = Tmp1(22,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,40)
     Tmp1(23,4) = Tmp1(23,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,41)
     Tmp1(24,4) = Tmp1(24,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,43)
     Tmp1(25,4) = Tmp1(25,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,44)
     Tmp1(26,4) = Tmp1(26,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,45)
     Tmp1(27,4) = Tmp1(27,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,47)
     Tmp1(28,4) = Tmp1(28,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,48)
     Tmp1(29,4) = Tmp1(29,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,49)
     Tmp1(30,4) = Tmp1(30,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,50)
     Tmp1(31,4) = Tmp1(31,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,52)
     Tmp1(32,4) = Tmp1(32,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,53)
     Tmp1(33,4) = Tmp1(33,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,54)
     Tmp1(34,4) = Tmp1(34,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,55)
     Tmp1(35,4) = Tmp1(35,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,56)
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
     Tmp0(2,5) = Tmp0(2,5) + inv2expP*Tmp0(1,2) 
     Tmp0(5,5) = Tmp0(5,5) + 2*inv2expP*Tmp0(2,2) 
     Tmp0(6,5) = Tmp0(6,5) + inv2expP*Tmp0(3,2) 
     Tmp0(7,5) = Tmp0(7,5) + inv2expP*Tmp0(4,2) 
     Tmp0(11,5) = Tmp0(11,5) + 3*inv2expP*Tmp0(5,2) 
     Tmp0(12,5) = Tmp0(12,5) + 2*inv2expP*Tmp0(6,2) 
     Tmp0(13,5) = Tmp0(13,5) + 2*inv2expP*Tmp0(7,2) 
     Tmp0(14,5) = Tmp0(14,5) + inv2expP*Tmp0(8,2) 
     Tmp0(15,5) = Tmp0(15,5) + inv2expP*Tmp0(9,2) 
     Tmp0(16,5) = Tmp0(16,5) + inv2expP*Tmp0(10,2) 
     Tmp0(2,6) = Tmp0(2,6) + inv2expP*Tmp0(1,3) 
     Tmp0(5,6) = Tmp0(5,6) + 2*inv2expP*Tmp0(2,3) 
     Tmp0(6,6) = Tmp0(6,6) + inv2expP*Tmp0(3,3) 
     Tmp0(7,6) = Tmp0(7,6) + inv2expP*Tmp0(4,3) 
     Tmp0(11,6) = Tmp0(11,6) + 3*inv2expP*Tmp0(5,3) 
     Tmp0(12,6) = Tmp0(12,6) + 2*inv2expP*Tmp0(6,3) 
     Tmp0(13,6) = Tmp0(13,6) + 2*inv2expP*Tmp0(7,3) 
     Tmp0(14,6) = Tmp0(14,6) + inv2expP*Tmp0(8,3) 
     Tmp0(15,6) = Tmp0(15,6) + inv2expP*Tmp0(9,3) 
     Tmp0(16,6) = Tmp0(16,6) + inv2expP*Tmp0(10,3) 
     Tmp0(2,7) = Tmp0(2,7) + inv2expP*Tmp0(1,4) 
     Tmp0(5,7) = Tmp0(5,7) + 2*inv2expP*Tmp0(2,4) 
     Tmp0(6,7) = Tmp0(6,7) + inv2expP*Tmp0(3,4) 
     Tmp0(7,7) = Tmp0(7,7) + inv2expP*Tmp0(4,4) 
     Tmp0(11,7) = Tmp0(11,7) + 3*inv2expP*Tmp0(5,4) 
     Tmp0(12,7) = Tmp0(12,7) + 2*inv2expP*Tmp0(6,4) 
     Tmp0(13,7) = Tmp0(13,7) + 2*inv2expP*Tmp0(7,4) 
     Tmp0(14,7) = Tmp0(14,7) + inv2expP*Tmp0(8,4) 
     Tmp0(15,7) = Tmp0(15,7) + inv2expP*Tmp0(9,4) 
     Tmp0(16,7) = Tmp0(16,7) + inv2expP*Tmp0(10,4) 
     Tmp0(3,8) = Tmp0(3,8) + inv2expP*Tmp0(1,3) 
     Tmp0(6,8) = Tmp0(6,8) + inv2expP*Tmp0(2,3) 
     Tmp0(8,8) = Tmp0(8,8) + 2*inv2expP*Tmp0(3,3) 
     Tmp0(9,8) = Tmp0(9,8) + inv2expP*Tmp0(4,3) 
     Tmp0(12,8) = Tmp0(12,8) + inv2expP*Tmp0(5,3) 
     Tmp0(14,8) = Tmp0(14,8) + 2*inv2expP*Tmp0(6,3) 
     Tmp0(15,8) = Tmp0(15,8) + inv2expP*Tmp0(7,3) 
     Tmp0(17,8) = Tmp0(17,8) + 3*inv2expP*Tmp0(8,3) 
     Tmp0(18,8) = Tmp0(18,8) + 2*inv2expP*Tmp0(9,3) 
     Tmp0(19,8) = Tmp0(19,8) + inv2expP*Tmp0(10,3) 
     Tmp0(3,9) = Tmp0(3,9) + inv2expP*Tmp0(1,4) 
     Tmp0(6,9) = Tmp0(6,9) + inv2expP*Tmp0(2,4) 
     Tmp0(8,9) = Tmp0(8,9) + 2*inv2expP*Tmp0(3,4) 
     Tmp0(9,9) = Tmp0(9,9) + inv2expP*Tmp0(4,4) 
     Tmp0(12,9) = Tmp0(12,9) + inv2expP*Tmp0(5,4) 
     Tmp0(14,9) = Tmp0(14,9) + 2*inv2expP*Tmp0(6,4) 
     Tmp0(15,9) = Tmp0(15,9) + inv2expP*Tmp0(7,4) 
     Tmp0(17,9) = Tmp0(17,9) + 3*inv2expP*Tmp0(8,4) 
     Tmp0(18,9) = Tmp0(18,9) + 2*inv2expP*Tmp0(9,4) 
     Tmp0(19,9) = Tmp0(19,9) + inv2expP*Tmp0(10,4) 
     Tmp0(4,10) = Tmp0(4,10) + inv2expP*Tmp0(1,4) 
     Tmp0(7,10) = Tmp0(7,10) + inv2expP*Tmp0(2,4) 
     Tmp0(9,10) = Tmp0(9,10) + inv2expP*Tmp0(3,4) 
     Tmp0(10,10) = Tmp0(10,10) + 2*inv2expP*Tmp0(4,4) 
     Tmp0(13,10) = Tmp0(13,10) + inv2expP*Tmp0(5,4) 
     Tmp0(15,10) = Tmp0(15,10) + inv2expP*Tmp0(6,4) 
     Tmp0(16,10) = Tmp0(16,10) + 2*inv2expP*Tmp0(7,4) 
     Tmp0(18,10) = Tmp0(18,10) + inv2expP*Tmp0(8,4) 
     Tmp0(19,10) = Tmp0(19,10) + 2*inv2expP*Tmp0(9,4) 
     Tmp0(20,10) = Tmp0(20,10) + 3*inv2expP*Tmp0(10,4) 
     Tmp0(1,5) = Tmp0(1,5) + qinvp*Tmp0(2,2)
     Tmp0(2,5) = Tmp0(2,5) + qinvp*Tmp0(5,2)
     Tmp0(3,5) = Tmp0(3,5) + qinvp*Tmp0(6,2)
     Tmp0(4,5) = Tmp0(4,5) + qinvp*Tmp0(7,2)
     Tmp0(5,5) = Tmp0(5,5) + qinvp*Tmp0(11,2)
     Tmp0(6,5) = Tmp0(6,5) + qinvp*Tmp0(12,2)
     Tmp0(7,5) = Tmp0(7,5) + qinvp*Tmp0(13,2)
     Tmp0(8,5) = Tmp0(8,5) + qinvp*Tmp0(14,2)
     Tmp0(9,5) = Tmp0(9,5) + qinvp*Tmp0(15,2)
     Tmp0(10,5) = Tmp0(10,5) + qinvp*Tmp0(16,2)
     Tmp0(11,5) = Tmp0(11,5) + qinvp*Tmp1(21,2)
     Tmp0(12,5) = Tmp0(12,5) + qinvp*Tmp1(22,2)
     Tmp0(13,5) = Tmp0(13,5) + qinvp*Tmp1(23,2)
     Tmp0(14,5) = Tmp0(14,5) + qinvp*Tmp1(24,2)
     Tmp0(15,5) = Tmp0(15,5) + qinvp*Tmp1(25,2)
     Tmp0(16,5) = Tmp0(16,5) + qinvp*Tmp1(26,2)
     Tmp0(17,5) = Tmp0(17,5) + qinvp*Tmp1(27,2)
     Tmp0(18,5) = Tmp0(18,5) + qinvp*Tmp1(28,2)
     Tmp0(19,5) = Tmp0(19,5) + qinvp*Tmp1(29,2)
     Tmp0(20,5) = Tmp0(20,5) + qinvp*Tmp1(30,2)
     Tmp0(1,6) = Tmp0(1,6) + qinvp*Tmp0(2,3)
     Tmp0(2,6) = Tmp0(2,6) + qinvp*Tmp0(5,3)
     Tmp0(3,6) = Tmp0(3,6) + qinvp*Tmp0(6,3)
     Tmp0(4,6) = Tmp0(4,6) + qinvp*Tmp0(7,3)
     Tmp0(5,6) = Tmp0(5,6) + qinvp*Tmp0(11,3)
     Tmp0(6,6) = Tmp0(6,6) + qinvp*Tmp0(12,3)
     Tmp0(7,6) = Tmp0(7,6) + qinvp*Tmp0(13,3)
     Tmp0(8,6) = Tmp0(8,6) + qinvp*Tmp0(14,3)
     Tmp0(9,6) = Tmp0(9,6) + qinvp*Tmp0(15,3)
     Tmp0(10,6) = Tmp0(10,6) + qinvp*Tmp0(16,3)
     Tmp0(11,6) = Tmp0(11,6) + qinvp*Tmp1(21,3)
     Tmp0(12,6) = Tmp0(12,6) + qinvp*Tmp1(22,3)
     Tmp0(13,6) = Tmp0(13,6) + qinvp*Tmp1(23,3)
     Tmp0(14,6) = Tmp0(14,6) + qinvp*Tmp1(24,3)
     Tmp0(15,6) = Tmp0(15,6) + qinvp*Tmp1(25,3)
     Tmp0(16,6) = Tmp0(16,6) + qinvp*Tmp1(26,3)
     Tmp0(17,6) = Tmp0(17,6) + qinvp*Tmp1(27,3)
     Tmp0(18,6) = Tmp0(18,6) + qinvp*Tmp1(28,3)
     Tmp0(19,6) = Tmp0(19,6) + qinvp*Tmp1(29,3)
     Tmp0(20,6) = Tmp0(20,6) + qinvp*Tmp1(30,3)
     Tmp0(1,7) = Tmp0(1,7) + qinvp*Tmp0(2,4)
     Tmp0(2,7) = Tmp0(2,7) + qinvp*Tmp0(5,4)
     Tmp0(3,7) = Tmp0(3,7) + qinvp*Tmp0(6,4)
     Tmp0(4,7) = Tmp0(4,7) + qinvp*Tmp0(7,4)
     Tmp0(5,7) = Tmp0(5,7) + qinvp*Tmp0(11,4)
     Tmp0(6,7) = Tmp0(6,7) + qinvp*Tmp0(12,4)
     Tmp0(7,7) = Tmp0(7,7) + qinvp*Tmp0(13,4)
     Tmp0(8,7) = Tmp0(8,7) + qinvp*Tmp0(14,4)
     Tmp0(9,7) = Tmp0(9,7) + qinvp*Tmp0(15,4)
     Tmp0(10,7) = Tmp0(10,7) + qinvp*Tmp0(16,4)
     Tmp0(11,7) = Tmp0(11,7) + qinvp*Tmp1(21,4)
     Tmp0(12,7) = Tmp0(12,7) + qinvp*Tmp1(22,4)
     Tmp0(13,7) = Tmp0(13,7) + qinvp*Tmp1(23,4)
     Tmp0(14,7) = Tmp0(14,7) + qinvp*Tmp1(24,4)
     Tmp0(15,7) = Tmp0(15,7) + qinvp*Tmp1(25,4)
     Tmp0(16,7) = Tmp0(16,7) + qinvp*Tmp1(26,4)
     Tmp0(17,7) = Tmp0(17,7) + qinvp*Tmp1(27,4)
     Tmp0(18,7) = Tmp0(18,7) + qinvp*Tmp1(28,4)
     Tmp0(19,7) = Tmp0(19,7) + qinvp*Tmp1(29,4)
     Tmp0(20,7) = Tmp0(20,7) + qinvp*Tmp1(30,4)
     Tmp0(1,8) = Tmp0(1,8) + qinvp*Tmp0(3,3)
     Tmp0(2,8) = Tmp0(2,8) + qinvp*Tmp0(6,3)
     Tmp0(3,8) = Tmp0(3,8) + qinvp*Tmp0(8,3)
     Tmp0(4,8) = Tmp0(4,8) + qinvp*Tmp0(9,3)
     Tmp0(5,8) = Tmp0(5,8) + qinvp*Tmp0(12,3)
     Tmp0(6,8) = Tmp0(6,8) + qinvp*Tmp0(14,3)
     Tmp0(7,8) = Tmp0(7,8) + qinvp*Tmp0(15,3)
     Tmp0(8,8) = Tmp0(8,8) + qinvp*Tmp0(17,3)
     Tmp0(9,8) = Tmp0(9,8) + qinvp*Tmp0(18,3)
     Tmp0(10,8) = Tmp0(10,8) + qinvp*Tmp0(19,3)
     Tmp0(11,8) = Tmp0(11,8) + qinvp*Tmp1(22,3)
     Tmp0(12,8) = Tmp0(12,8) + qinvp*Tmp1(24,3)
     Tmp0(13,8) = Tmp0(13,8) + qinvp*Tmp1(25,3)
     Tmp0(14,8) = Tmp0(14,8) + qinvp*Tmp1(27,3)
     Tmp0(15,8) = Tmp0(15,8) + qinvp*Tmp1(28,3)
     Tmp0(16,8) = Tmp0(16,8) + qinvp*Tmp1(29,3)
     Tmp0(17,8) = Tmp0(17,8) + qinvp*Tmp1(31,3)
     Tmp0(18,8) = Tmp0(18,8) + qinvp*Tmp1(32,3)
     Tmp0(19,8) = Tmp0(19,8) + qinvp*Tmp1(33,3)
     Tmp0(20,8) = Tmp0(20,8) + qinvp*Tmp1(34,3)
     Tmp0(1,9) = Tmp0(1,9) + qinvp*Tmp0(3,4)
     Tmp0(2,9) = Tmp0(2,9) + qinvp*Tmp0(6,4)
     Tmp0(3,9) = Tmp0(3,9) + qinvp*Tmp0(8,4)
     Tmp0(4,9) = Tmp0(4,9) + qinvp*Tmp0(9,4)
     Tmp0(5,9) = Tmp0(5,9) + qinvp*Tmp0(12,4)
     Tmp0(6,9) = Tmp0(6,9) + qinvp*Tmp0(14,4)
     Tmp0(7,9) = Tmp0(7,9) + qinvp*Tmp0(15,4)
     Tmp0(8,9) = Tmp0(8,9) + qinvp*Tmp0(17,4)
     Tmp0(9,9) = Tmp0(9,9) + qinvp*Tmp0(18,4)
     Tmp0(10,9) = Tmp0(10,9) + qinvp*Tmp0(19,4)
     Tmp0(11,9) = Tmp0(11,9) + qinvp*Tmp1(22,4)
     Tmp0(12,9) = Tmp0(12,9) + qinvp*Tmp1(24,4)
     Tmp0(13,9) = Tmp0(13,9) + qinvp*Tmp1(25,4)
     Tmp0(14,9) = Tmp0(14,9) + qinvp*Tmp1(27,4)
     Tmp0(15,9) = Tmp0(15,9) + qinvp*Tmp1(28,4)
     Tmp0(16,9) = Tmp0(16,9) + qinvp*Tmp1(29,4)
     Tmp0(17,9) = Tmp0(17,9) + qinvp*Tmp1(31,4)
     Tmp0(18,9) = Tmp0(18,9) + qinvp*Tmp1(32,4)
     Tmp0(19,9) = Tmp0(19,9) + qinvp*Tmp1(33,4)
     Tmp0(20,9) = Tmp0(20,9) + qinvp*Tmp1(34,4)
     Tmp0(1,10) = Tmp0(1,10) + qinvp*Tmp0(4,4)
     Tmp0(2,10) = Tmp0(2,10) + qinvp*Tmp0(7,4)
     Tmp0(3,10) = Tmp0(3,10) + qinvp*Tmp0(9,4)
     Tmp0(4,10) = Tmp0(4,10) + qinvp*Tmp0(10,4)
     Tmp0(5,10) = Tmp0(5,10) + qinvp*Tmp0(13,4)
     Tmp0(6,10) = Tmp0(6,10) + qinvp*Tmp0(15,4)
     Tmp0(7,10) = Tmp0(7,10) + qinvp*Tmp0(16,4)
     Tmp0(8,10) = Tmp0(8,10) + qinvp*Tmp0(18,4)
     Tmp0(9,10) = Tmp0(9,10) + qinvp*Tmp0(19,4)
     Tmp0(10,10) = Tmp0(10,10) + qinvp*Tmp0(20,4)
     Tmp0(11,10) = Tmp0(11,10) + qinvp*Tmp1(23,4)
     Tmp0(12,10) = Tmp0(12,10) + qinvp*Tmp1(25,4)
     Tmp0(13,10) = Tmp0(13,10) + qinvp*Tmp1(26,4)
     Tmp0(14,10) = Tmp0(14,10) + qinvp*Tmp1(28,4)
     Tmp0(15,10) = Tmp0(15,10) + qinvp*Tmp1(29,4)
     Tmp0(16,10) = Tmp0(16,10) + qinvp*Tmp1(30,4)
     Tmp0(17,10) = Tmp0(17,10) + qinvp*Tmp1(32,4)
     Tmp0(18,10) = Tmp0(18,10) + qinvp*Tmp1(33,4)
     Tmp0(19,10) = Tmp0(19,10) + qinvp*Tmp1(34,4)
     Tmp0(20,10) = Tmp0(20,10) + qinvp*Tmp1(35,4)
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
