MODULE AGC_GPU_OBS_TRMODCtoBSeg
 use IchorPrecisionModule
  
 CONTAINS
 subroutine TransferRecurrenceGPUP1Q2CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Aexp(nPrimA)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   20)
  real(realk),intent(inout) :: Aux2(nPasses,    4,   10)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0( 10,  4)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimA
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPasses,Aux2)
  DO iP = 1,nPasses
   DO iTUVQ=1, 10
    DO iTUVP=1,  4
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&
!$ACC         expAX,expAY,expAZ,&
!$ACC         iPrimD,iPrimA,&
!$ACC         Tmp0,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        Aexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux)
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimA = iPrimP - ((iPrimP-1)/nPrimA)*nPrimA
    expAX = -Aexp(iPrimA)*Xab
    expAY = -Aexp(iPrimA)*Yab
    expAZ = -Aexp(iPrimA)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expAX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expAY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expAZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 10
      Tmp0(iTUVQ,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 10
      Tmp0(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 10
      Tmp0(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
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
     DO iTUVQ=1, 10
      DO iTUVP=1,  4
        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP1Q2CtoBSeg
 subroutine TransferRecurrenceGPUP1Q3CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Aexp(nPrimA)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   35)
  real(realk),intent(inout) :: Aux2(nPasses,    4,   20)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0( 20,  4)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimA
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPasses,Aux2)
  DO iP = 1,nPasses
   DO iTUVQ=1, 20
    DO iTUVP=1,  4
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&
!$ACC         expAX,expAY,expAZ,&
!$ACC         iPrimD,iPrimA,&
!$ACC         Tmp0,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        Aexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux)
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimA = iPrimP - ((iPrimP-1)/nPrimA)*nPrimA
    expAX = -Aexp(iPrimA)*Xab
    expAY = -Aexp(iPrimA)*Yab
    expAZ = -Aexp(iPrimA)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expAX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expAY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expAZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 20
      Tmp0(iTUVQ,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
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
     DO iTUVQ=1, 20
      DO iTUVP=1,  4
        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP1Q3CtoBSeg
 subroutine TransferRecurrenceGPUP1Q4CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Aexp(nPrimA)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   56)
  real(realk),intent(inout) :: Aux2(nPasses,    4,   35)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0( 35,  4)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimA
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
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
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPasses,Aux2)
  DO iP = 1,nPasses
   DO iTUVQ=1, 35
    DO iTUVP=1,  4
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&
!$ACC         expAX,expAY,expAZ,&
!$ACC         iPrimD,iPrimA,&
!$ACC         Tmp0,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        Aexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux)
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimA = iPrimP - ((iPrimP-1)/nPrimA)*nPrimA
    expAX = -Aexp(iPrimA)*Xab
    expAY = -Aexp(iPrimA)*Yab
    expAZ = -Aexp(iPrimA)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expAX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expAY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expAZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 35
      Tmp0(iTUVQ,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
     DO iTUVQ=1, 35
      DO iTUVP=1,  4
        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP1Q4CtoBSeg
 subroutine TransferRecurrenceGPUP2Q3CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Aexp(nPrimA)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   56)
  real(realk),intent(inout) :: Aux2(nPasses,   10,   20)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0( 20, 10)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
  real(realk) :: Tmp1( 21: 35,  2:  4)
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimA
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
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
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPasses,Aux2)
  DO iP = 1,nPasses
   DO iTUVQ=1, 20
    DO iTUVP=1, 10
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&
!$ACC         expAX,expAY,expAZ,&
!$ACC         iPrimD,iPrimA,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        Aexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux)
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimA = iPrimP - ((iPrimP-1)/nPrimA)*nPrimA
    expAX = -Aexp(iPrimA)*Xab
    expAY = -Aexp(iPrimA)*Yab
    expAZ = -Aexp(iPrimA)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expAX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expAY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expAZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 20
      Tmp0(iTUVQ,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + IfacX1(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + IfacX2(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + IfacX3(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jp = 2
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,5) = facX*Tmp0(iTUVQ,2)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,6) = facX*Tmp0(iTUVQ,3)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,7) = facX*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,8) = facY*Tmp0(iTUVQ,3)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,9) = facY*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,10) = facZ*Tmp0(iTUVQ,4)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + IfacX2(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + IfacX2(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + IfacX3(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp0(iTUVplus1,2)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
     DO iTUVQ=1, 20
      DO iTUVP=1, 10
        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP2Q3CtoBSeg
 subroutine TransferRecurrenceGPUP2Q4CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Aexp(nPrimA)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,   84)
  real(realk),intent(inout) :: Aux2(nPasses,   10,   35)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0( 35, 10)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
  real(realk) :: Tmp1( 36: 56,  2:  4)
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimA
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
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
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPasses,Aux2)
  DO iP = 1,nPasses
   DO iTUVQ=1, 35
    DO iTUVP=1, 10
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&
!$ACC         expAX,expAY,expAZ,&
!$ACC         iPrimD,iPrimA,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        Aexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux)
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimA = iPrimP - ((iPrimP-1)/nPrimA)*nPrimA
    expAX = -Aexp(iPrimA)*Xab
    expAY = -Aexp(iPrimA)*Yab
    expAZ = -Aexp(iPrimA)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expAX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expAY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expAZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 35
      Tmp0(iTUVQ,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  36, 56
      Tmp1(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  36, 56
      Tmp1(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  36, 56
      Tmp1(iTUVQ,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + IfacX1(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + IfacX2(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + IfacX3(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jp = 2
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,5) = facX*Tmp0(iTUVQ,2)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,6) = facX*Tmp0(iTUVQ,3)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,7) = facX*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,8) = facY*Tmp0(iTUVQ,3)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,9) = facY*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,10) = facZ*Tmp0(iTUVQ,4)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + IfacX2(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + IfacX2(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + IfacX3(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp0(iTUVplus1,2)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
     DO iTUVQ=1, 35
      DO iTUVP=1, 10
        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP2Q4CtoBSeg
 subroutine TransferRecurrenceGPUP3Q4CtoBSeg(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Aexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Aexp(nPrimA)
  real(realk),intent(in) :: Aux(nPrimQ,nPrimP,nPasses,  120)
  real(realk),intent(inout) :: Aux2(nPasses,   20,   35)
!  real(realk),intent(inout) :: Aux2(nPrimQ,nPrimP,nPasses,nTUVP,nTUVQ)
  !Local variables
  real(realk) :: Tmp0( 35, 20)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
  real(realk) :: Tmp1( 36: 84,  2:  4)
  real(realk) :: Tmp2( 36: 56,  5: 10)
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimA
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expAX,expAY,expAZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
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
!$ACC PARALLEL LOOP PRIVATE(iP,iTUVP,iTUVQ) PRESENT(nPasses,Aux2)
  DO iP = 1,nPasses
   DO iTUVQ=1, 35
    DO iTUVP=1, 20
     Aux2(iP,iTUVP,iTUVQ) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP,&
!$ACC         expAX,expAY,expAZ,&
!$ACC         iPrimD,iPrimA,&
!$ACC         Tmp0,&
!$ACC         Tmp1,&
!$ACC         Tmp2,&
!$ACC         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) &
!$ACC PRESENT(nPasses,nPrimP,nPrimQ,nPrimA,nPrimC,reducedExponents,Pexp,Qexp,&
!$ACC        Aexp,Dexp,&
!$ACC        Pdistance12,Qdistance12,IatomApass,IatomBpass,Aux2,Aux)
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimA = iPrimP - ((iPrimP-1)/nPrimA)*nPrimA
    expAX = -Aexp(iPrimA)*Xab
    expAY = -Aexp(iPrimA)*Yab
    expAZ = -Aexp(iPrimA)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expAX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expAY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expAZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 35
      Tmp0(iTUVQ,1) = Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  36, 84
      Tmp1(iTUVQ,2) = facX*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  36, 84
      Tmp1(iTUVQ,3) = facY*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  36, 84
      Tmp1(iTUVQ,4) = facZ*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 21,56
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + IfacX1(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 21,56
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + IfacX2(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do ituvqminus1 = 21,56
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + IfacX3(ituvqminus1)*inv2expP*Aux(iPrimQ,iPrimP,iPassP,ituvqminus1) 
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 36,84
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 36,84
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
     do iTUVQ = 36,84
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + qinvp*Aux(iPrimQ,iPrimP,iPassP,iTUVplus1)
     enddo
 ! Building for Angular momentum Jp = 2
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,5) = facX*Tmp0(iTUVQ,2)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,6) = facX*Tmp0(iTUVQ,3)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,7) = facX*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,8) = facY*Tmp0(iTUVQ,3)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,9) = facY*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,10) = facZ*Tmp0(iTUVQ,4)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,5) = facX*Tmp1(iTUVQ,2)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,6) = facX*Tmp1(iTUVQ,3)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,7) = facX*Tmp1(iTUVQ,4)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,8) = facY*Tmp1(iTUVQ,3)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,9) = facY*Tmp1(iTUVQ,4)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,10) = facZ*Tmp1(iTUVQ,4)+ inv2expP*Aux(iPrimQ,iPrimP,iPassP,iTUVQ)
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + IfacX2(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + IfacX2(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + IfacX3(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp2(iTUVQ,5) = Tmp2(iTUVQ,5) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp2(iTUVQ,6) = Tmp2(iTUVQ,6) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp2(iTUVQ,7) = Tmp2(iTUVQ,7) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp2(iTUVQ,8) = Tmp2(iTUVQ,8) + IfacX2(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp2(iTUVQ,9) = Tmp2(iTUVQ,9) + IfacX2(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp2(iTUVQ,10) = Tmp2(iTUVQ,10) + IfacX3(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp0(iTUVplus1,2)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp2(iTUVQ,5) = Tmp2(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp2(iTUVQ,6) = Tmp2(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp2(iTUVQ,7) = Tmp2(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp2(iTUVQ,8) = Tmp2(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp2(iTUVQ,9) = Tmp2(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp2(iTUVQ,10) = Tmp2(iTUVQ,10) + qinvp*Tmp1(iTUVplus1,4)
     enddo
 ! Building for Angular momentum Jp = 3
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,11) = facX*Tmp0(iTUVQ,5)+2*inv2expP*Tmp0(iTUVQ,2)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,12) = facY*Tmp0(iTUVQ,5)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,13) = facZ*Tmp0(iTUVQ,5)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,14) = facX*Tmp0(iTUVQ,8)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,15) = facX*Tmp0(iTUVQ,9)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,16) = facX*Tmp0(iTUVQ,10)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,17) = facY*Tmp0(iTUVQ,8)+2*inv2expP*Tmp0(iTUVQ,3)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,18) = facZ*Tmp0(iTUVQ,8)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,19) = facY*Tmp0(iTUVQ,10)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,20) = facZ*Tmp0(iTUVQ,10)+2*inv2expP*Tmp0(iTUVQ,4)
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,11) = Tmp0(iTUVQ,11) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,5) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,12) = Tmp0(iTUVQ,12) + IfacX2(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,5) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp0(iTUVQ,13) = Tmp0(iTUVQ,13) + IfacX3(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,5) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,14) = Tmp0(iTUVQ,14) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,8) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,15) = Tmp0(iTUVQ,15) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,9) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,16) = Tmp0(iTUVQ,16) + IfacX1(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,10) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,17) = Tmp0(iTUVQ,17) + IfacX2(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,8) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp0(iTUVQ,18) = Tmp0(iTUVQ,18) + IfacX3(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,8) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,19) = Tmp0(iTUVQ,19) + IfacX2(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,10) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp0(iTUVQ,20) = Tmp0(iTUVQ,20) + IfacX3(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,10) 
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,11) = Tmp0(iTUVQ,11) + qinvp*Tmp0(iTUVplus1,5)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,11) = Tmp0(iTUVQ,11) + qinvp*Tmp2(iTUVplus1,5)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,12) = Tmp0(iTUVQ,12) + qinvp*Tmp0(iTUVplus1,5)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,12) = Tmp0(iTUVQ,12) + qinvp*Tmp2(iTUVplus1,5)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,13) = Tmp0(iTUVQ,13) + qinvp*Tmp0(iTUVplus1,5)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,13) = Tmp0(iTUVQ,13) + qinvp*Tmp2(iTUVplus1,5)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,14) = Tmp0(iTUVQ,14) + qinvp*Tmp0(iTUVplus1,8)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,14) = Tmp0(iTUVQ,14) + qinvp*Tmp2(iTUVplus1,8)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,15) = Tmp0(iTUVQ,15) + qinvp*Tmp0(iTUVplus1,9)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,15) = Tmp0(iTUVQ,15) + qinvp*Tmp2(iTUVplus1,9)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,16) = Tmp0(iTUVQ,16) + qinvp*Tmp0(iTUVplus1,10)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,16) = Tmp0(iTUVQ,16) + qinvp*Tmp2(iTUVplus1,10)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,17) = Tmp0(iTUVQ,17) + qinvp*Tmp0(iTUVplus1,8)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,17) = Tmp0(iTUVQ,17) + qinvp*Tmp2(iTUVplus1,8)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,18) = Tmp0(iTUVQ,18) + qinvp*Tmp0(iTUVplus1,8)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,18) = Tmp0(iTUVQ,18) + qinvp*Tmp2(iTUVplus1,8)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,19) = Tmp0(iTUVQ,19) + qinvp*Tmp0(iTUVplus1,10)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,19) = Tmp0(iTUVQ,19) + qinvp*Tmp2(iTUVplus1,10)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,20) = Tmp0(iTUVQ,20) + qinvp*Tmp0(iTUVplus1,10)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,20) = Tmp0(iTUVQ,20) + qinvp*Tmp2(iTUVplus1,10)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
     DO iTUVQ=1, 35
      DO iTUVP=1, 20
        Aux2(IP,iTUVP,iTUVQ) = Aux2(iP,iTUVP,iTUVQ) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPasses
 end subroutine TransferRecurrenceGPUP3Q4CtoBSeg
end module
