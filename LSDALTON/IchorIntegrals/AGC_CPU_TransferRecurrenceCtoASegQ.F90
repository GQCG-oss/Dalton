MODULE AGC_CPU_OBS_TRMODCtoASegQ
 use IchorPrecisionMod
  
 CONTAINS
 subroutine TransferRecurrenceCPUP1Q2CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(   20,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(    4,   10,nPrimP*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 10,  4)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$OMP DO COLLAPSE(3) PRIVATE(iP,iTUVP,iTUVQ)
  DO iP = 1,nPrimP*nPasses
   DO iTUVQ=1, 10
    DO iTUVP=1,  4
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimD,iPrimB,&
!$OMP         Tmp0,&
!$OMP         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) 
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 10
      Tmp0(iTUVQ,1) = Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 10
      Tmp0(iTUVQ,2) = facX*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 10
      Tmp0(iTUVQ,3) = facY*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 10
      Tmp0(iTUVQ,4) = facZ*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expP*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expP*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expP*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expP*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expP*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expP*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expP*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expP*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expP*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expP*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expP*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expP*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(1,2) = Tmp0(1,2) + qinvp*Aux(2,iPrimQ,iPrimP,iPassP)
     Tmp0(2,2) = Tmp0(2,2) + qinvp*Aux(5,iPrimQ,iPrimP,iPassP)
     Tmp0(3,2) = Tmp0(3,2) + qinvp*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(4,2) = Tmp0(4,2) + qinvp*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(5,2) = Tmp0(5,2) + qinvp*Aux(11,iPrimQ,iPrimP,iPassP)
     Tmp0(6,2) = Tmp0(6,2) + qinvp*Aux(12,iPrimQ,iPrimP,iPassP)
     Tmp0(7,2) = Tmp0(7,2) + qinvp*Aux(13,iPrimQ,iPrimP,iPassP)
     Tmp0(8,2) = Tmp0(8,2) + qinvp*Aux(14,iPrimQ,iPrimP,iPassP)
     Tmp0(9,2) = Tmp0(9,2) + qinvp*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(10,2) = Tmp0(10,2) + qinvp*Aux(16,iPrimQ,iPrimP,iPassP)
     Tmp0(1,3) = Tmp0(1,3) + qinvp*Aux(3,iPrimQ,iPrimP,iPassP)
     Tmp0(2,3) = Tmp0(2,3) + qinvp*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(3,3) = Tmp0(3,3) + qinvp*Aux(8,iPrimQ,iPrimP,iPassP)
     Tmp0(4,3) = Tmp0(4,3) + qinvp*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(5,3) = Tmp0(5,3) + qinvp*Aux(12,iPrimQ,iPrimP,iPassP)
     Tmp0(6,3) = Tmp0(6,3) + qinvp*Aux(14,iPrimQ,iPrimP,iPassP)
     Tmp0(7,3) = Tmp0(7,3) + qinvp*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(8,3) = Tmp0(8,3) + qinvp*Aux(17,iPrimQ,iPrimP,iPassP)
     Tmp0(9,3) = Tmp0(9,3) + qinvp*Aux(18,iPrimQ,iPrimP,iPassP)
     Tmp0(10,3) = Tmp0(10,3) + qinvp*Aux(19,iPrimQ,iPrimP,iPassP)
     Tmp0(1,4) = Tmp0(1,4) + qinvp*Aux(4,iPrimQ,iPrimP,iPassP)
     Tmp0(2,4) = Tmp0(2,4) + qinvp*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(3,4) = Tmp0(3,4) + qinvp*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(4,4) = Tmp0(4,4) + qinvp*Aux(10,iPrimQ,iPrimP,iPassP)
     Tmp0(5,4) = Tmp0(5,4) + qinvp*Aux(13,iPrimQ,iPrimP,iPassP)
     Tmp0(6,4) = Tmp0(6,4) + qinvp*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(7,4) = Tmp0(7,4) + qinvp*Aux(16,iPrimQ,iPrimP,iPassP)
     Tmp0(8,4) = Tmp0(8,4) + qinvp*Aux(18,iPrimQ,iPrimP,iPassP)
     Tmp0(9,4) = Tmp0(9,4) + qinvp*Aux(19,iPrimQ,iPrimP,iPassP)
     Tmp0(10,4) = Tmp0(10,4) + qinvp*Aux(20,iPrimQ,iPrimP,iPassP)
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
     DO iTUVQ=1, 10
      DO iTUVP=1,  4
        Aux2(iTUVP,iTUVQ,iP) = Aux2(iTUVP,iTUVQ,iP) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimP=1, nPrimP
  ENDDO !iP = 1,nPrimQ*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP1Q2CtoASegQ
 subroutine TransferRecurrenceCPUP1Q3CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(   35,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(    4,   20,nPrimP*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 20,  4)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$OMP DO COLLAPSE(3) PRIVATE(iP,iTUVP,iTUVQ)
  DO iP = 1,nPrimP*nPasses
   DO iTUVQ=1, 20
    DO iTUVP=1,  4
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimD,iPrimB,&
!$OMP         Tmp0,&
!$OMP         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) 
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 20
      Tmp0(iTUVQ,1) = Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,2) = facX*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,3) = facY*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,4) = facZ*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     Tmp0(2,2) = Tmp0(2,2) + inv2expP*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(5,2) = Tmp0(5,2) + 2*inv2expP*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(6,2) = Tmp0(6,2) + inv2expP*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(7,2) = Tmp0(7,2) + inv2expP*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(11,2) = Tmp0(11,2) + 3*inv2expP*Aux(5,iPrimQ,iPrimP,iPassP) 
     Tmp0(12,2) = Tmp0(12,2) + 2*inv2expP*Aux(6,iPrimQ,iPrimP,iPassP) 
     Tmp0(13,2) = Tmp0(13,2) + 2*inv2expP*Aux(7,iPrimQ,iPrimP,iPassP) 
     Tmp0(14,2) = Tmp0(14,2) + inv2expP*Aux(8,iPrimQ,iPrimP,iPassP) 
     Tmp0(15,2) = Tmp0(15,2) + inv2expP*Aux(9,iPrimQ,iPrimP,iPassP) 
     Tmp0(16,2) = Tmp0(16,2) + inv2expP*Aux(10,iPrimQ,iPrimP,iPassP) 
     Tmp0(3,3) = Tmp0(3,3) + inv2expP*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(6,3) = Tmp0(6,3) + inv2expP*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(8,3) = Tmp0(8,3) + 2*inv2expP*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(9,3) = Tmp0(9,3) + inv2expP*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(12,3) = Tmp0(12,3) + inv2expP*Aux(5,iPrimQ,iPrimP,iPassP) 
     Tmp0(14,3) = Tmp0(14,3) + 2*inv2expP*Aux(6,iPrimQ,iPrimP,iPassP) 
     Tmp0(15,3) = Tmp0(15,3) + inv2expP*Aux(7,iPrimQ,iPrimP,iPassP) 
     Tmp0(17,3) = Tmp0(17,3) + 3*inv2expP*Aux(8,iPrimQ,iPrimP,iPassP) 
     Tmp0(18,3) = Tmp0(18,3) + 2*inv2expP*Aux(9,iPrimQ,iPrimP,iPassP) 
     Tmp0(19,3) = Tmp0(19,3) + inv2expP*Aux(10,iPrimQ,iPrimP,iPassP) 
     Tmp0(4,4) = Tmp0(4,4) + inv2expP*Aux(1,iPrimQ,iPrimP,iPassP) 
     Tmp0(7,4) = Tmp0(7,4) + inv2expP*Aux(2,iPrimQ,iPrimP,iPassP) 
     Tmp0(9,4) = Tmp0(9,4) + inv2expP*Aux(3,iPrimQ,iPrimP,iPassP) 
     Tmp0(10,4) = Tmp0(10,4) + 2*inv2expP*Aux(4,iPrimQ,iPrimP,iPassP) 
     Tmp0(13,4) = Tmp0(13,4) + inv2expP*Aux(5,iPrimQ,iPrimP,iPassP) 
     Tmp0(15,4) = Tmp0(15,4) + inv2expP*Aux(6,iPrimQ,iPrimP,iPassP) 
     Tmp0(16,4) = Tmp0(16,4) + 2*inv2expP*Aux(7,iPrimQ,iPrimP,iPassP) 
     Tmp0(18,4) = Tmp0(18,4) + inv2expP*Aux(8,iPrimQ,iPrimP,iPassP) 
     Tmp0(19,4) = Tmp0(19,4) + 2*inv2expP*Aux(9,iPrimQ,iPrimP,iPassP) 
     Tmp0(20,4) = Tmp0(20,4) + 3*inv2expP*Aux(10,iPrimQ,iPrimP,iPassP) 
     Tmp0(1,2) = Tmp0(1,2) + qinvp*Aux(2,iPrimQ,iPrimP,iPassP)
     Tmp0(2,2) = Tmp0(2,2) + qinvp*Aux(5,iPrimQ,iPrimP,iPassP)
     Tmp0(3,2) = Tmp0(3,2) + qinvp*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(4,2) = Tmp0(4,2) + qinvp*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(5,2) = Tmp0(5,2) + qinvp*Aux(11,iPrimQ,iPrimP,iPassP)
     Tmp0(6,2) = Tmp0(6,2) + qinvp*Aux(12,iPrimQ,iPrimP,iPassP)
     Tmp0(7,2) = Tmp0(7,2) + qinvp*Aux(13,iPrimQ,iPrimP,iPassP)
     Tmp0(8,2) = Tmp0(8,2) + qinvp*Aux(14,iPrimQ,iPrimP,iPassP)
     Tmp0(9,2) = Tmp0(9,2) + qinvp*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(10,2) = Tmp0(10,2) + qinvp*Aux(16,iPrimQ,iPrimP,iPassP)
     Tmp0(11,2) = Tmp0(11,2) + qinvp*Aux(21,iPrimQ,iPrimP,iPassP)
     Tmp0(12,2) = Tmp0(12,2) + qinvp*Aux(22,iPrimQ,iPrimP,iPassP)
     Tmp0(13,2) = Tmp0(13,2) + qinvp*Aux(23,iPrimQ,iPrimP,iPassP)
     Tmp0(14,2) = Tmp0(14,2) + qinvp*Aux(24,iPrimQ,iPrimP,iPassP)
     Tmp0(15,2) = Tmp0(15,2) + qinvp*Aux(25,iPrimQ,iPrimP,iPassP)
     Tmp0(16,2) = Tmp0(16,2) + qinvp*Aux(26,iPrimQ,iPrimP,iPassP)
     Tmp0(17,2) = Tmp0(17,2) + qinvp*Aux(27,iPrimQ,iPrimP,iPassP)
     Tmp0(18,2) = Tmp0(18,2) + qinvp*Aux(28,iPrimQ,iPrimP,iPassP)
     Tmp0(19,2) = Tmp0(19,2) + qinvp*Aux(29,iPrimQ,iPrimP,iPassP)
     Tmp0(20,2) = Tmp0(20,2) + qinvp*Aux(30,iPrimQ,iPrimP,iPassP)
     Tmp0(1,3) = Tmp0(1,3) + qinvp*Aux(3,iPrimQ,iPrimP,iPassP)
     Tmp0(2,3) = Tmp0(2,3) + qinvp*Aux(6,iPrimQ,iPrimP,iPassP)
     Tmp0(3,3) = Tmp0(3,3) + qinvp*Aux(8,iPrimQ,iPrimP,iPassP)
     Tmp0(4,3) = Tmp0(4,3) + qinvp*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(5,3) = Tmp0(5,3) + qinvp*Aux(12,iPrimQ,iPrimP,iPassP)
     Tmp0(6,3) = Tmp0(6,3) + qinvp*Aux(14,iPrimQ,iPrimP,iPassP)
     Tmp0(7,3) = Tmp0(7,3) + qinvp*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(8,3) = Tmp0(8,3) + qinvp*Aux(17,iPrimQ,iPrimP,iPassP)
     Tmp0(9,3) = Tmp0(9,3) + qinvp*Aux(18,iPrimQ,iPrimP,iPassP)
     Tmp0(10,3) = Tmp0(10,3) + qinvp*Aux(19,iPrimQ,iPrimP,iPassP)
     Tmp0(11,3) = Tmp0(11,3) + qinvp*Aux(22,iPrimQ,iPrimP,iPassP)
     Tmp0(12,3) = Tmp0(12,3) + qinvp*Aux(24,iPrimQ,iPrimP,iPassP)
     Tmp0(13,3) = Tmp0(13,3) + qinvp*Aux(25,iPrimQ,iPrimP,iPassP)
     Tmp0(14,3) = Tmp0(14,3) + qinvp*Aux(27,iPrimQ,iPrimP,iPassP)
     Tmp0(15,3) = Tmp0(15,3) + qinvp*Aux(28,iPrimQ,iPrimP,iPassP)
     Tmp0(16,3) = Tmp0(16,3) + qinvp*Aux(29,iPrimQ,iPrimP,iPassP)
     Tmp0(17,3) = Tmp0(17,3) + qinvp*Aux(31,iPrimQ,iPrimP,iPassP)
     Tmp0(18,3) = Tmp0(18,3) + qinvp*Aux(32,iPrimQ,iPrimP,iPassP)
     Tmp0(19,3) = Tmp0(19,3) + qinvp*Aux(33,iPrimQ,iPrimP,iPassP)
     Tmp0(20,3) = Tmp0(20,3) + qinvp*Aux(34,iPrimQ,iPrimP,iPassP)
     Tmp0(1,4) = Tmp0(1,4) + qinvp*Aux(4,iPrimQ,iPrimP,iPassP)
     Tmp0(2,4) = Tmp0(2,4) + qinvp*Aux(7,iPrimQ,iPrimP,iPassP)
     Tmp0(3,4) = Tmp0(3,4) + qinvp*Aux(9,iPrimQ,iPrimP,iPassP)
     Tmp0(4,4) = Tmp0(4,4) + qinvp*Aux(10,iPrimQ,iPrimP,iPassP)
     Tmp0(5,4) = Tmp0(5,4) + qinvp*Aux(13,iPrimQ,iPrimP,iPassP)
     Tmp0(6,4) = Tmp0(6,4) + qinvp*Aux(15,iPrimQ,iPrimP,iPassP)
     Tmp0(7,4) = Tmp0(7,4) + qinvp*Aux(16,iPrimQ,iPrimP,iPassP)
     Tmp0(8,4) = Tmp0(8,4) + qinvp*Aux(18,iPrimQ,iPrimP,iPassP)
     Tmp0(9,4) = Tmp0(9,4) + qinvp*Aux(19,iPrimQ,iPrimP,iPassP)
     Tmp0(10,4) = Tmp0(10,4) + qinvp*Aux(20,iPrimQ,iPrimP,iPassP)
     Tmp0(11,4) = Tmp0(11,4) + qinvp*Aux(23,iPrimQ,iPrimP,iPassP)
     Tmp0(12,4) = Tmp0(12,4) + qinvp*Aux(25,iPrimQ,iPrimP,iPassP)
     Tmp0(13,4) = Tmp0(13,4) + qinvp*Aux(26,iPrimQ,iPrimP,iPassP)
     Tmp0(14,4) = Tmp0(14,4) + qinvp*Aux(28,iPrimQ,iPrimP,iPassP)
     Tmp0(15,4) = Tmp0(15,4) + qinvp*Aux(29,iPrimQ,iPrimP,iPassP)
     Tmp0(16,4) = Tmp0(16,4) + qinvp*Aux(30,iPrimQ,iPrimP,iPassP)
     Tmp0(17,4) = Tmp0(17,4) + qinvp*Aux(32,iPrimQ,iPrimP,iPassP)
     Tmp0(18,4) = Tmp0(18,4) + qinvp*Aux(33,iPrimQ,iPrimP,iPassP)
     Tmp0(19,4) = Tmp0(19,4) + qinvp*Aux(34,iPrimQ,iPrimP,iPassP)
     Tmp0(20,4) = Tmp0(20,4) + qinvp*Aux(35,iPrimQ,iPrimP,iPassP)
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
     DO iTUVQ=1, 20
      DO iTUVP=1,  4
        Aux2(iTUVP,iTUVQ,iP) = Aux2(iTUVP,iTUVQ,iP) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimP=1, nPrimP
  ENDDO !iP = 1,nPrimQ*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP1Q3CtoASegQ
 subroutine TransferRecurrenceCPUP1Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  use AGC_CPU_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(   56,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(    4,   35,nPrimP*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 35,  4)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$OMP DO COLLAPSE(3) PRIVATE(iP,iTUVP,iTUVQ)
  DO iP = 1,nPrimP*nPasses
   DO iTUVQ=1, 35
    DO iTUVP=1,  4
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimD,iPrimB,&
!$OMP         Tmp0,&
!$OMP         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) 
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 35
      Tmp0(iTUVQ,1) = Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,2) = facX*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,3) = facY*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,4) = facZ*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1_20(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2_20(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_35(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3_20(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
     DO iTUVQ=1, 35
      DO iTUVP=1,  4
        Aux2(iTUVP,iTUVQ,iP) = Aux2(iTUVP,iTUVQ,iP) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimP=1, nPrimP
  ENDDO !iP = 1,nPrimQ*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP1Q4CtoASegQ
 subroutine TransferRecurrenceCPUP2Q3CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  use AGC_CPU_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(   56,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   10,   20,nPrimP*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 20, 10)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
  real(realk) :: Tmp1( 21: 35,  2:  4)
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$OMP DO COLLAPSE(3) PRIVATE(iP,iTUVP,iTUVQ)
  DO iP = 1,nPrimP*nPasses
   DO iTUVQ=1, 20
    DO iTUVP=1, 10
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimD,iPrimB,&
!$OMP         Tmp0,&
!$OMP         Tmp1,&
!$OMP         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) 
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 20
      Tmp0(iTUVQ,1) = Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,2) = facX*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,3) = facY*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,4) = facZ*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,2) = facX*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,3) = facY*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  21, 35
      Tmp1(iTUVQ,4) = facZ*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1_20(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2_20(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX3_35(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3_20(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + IfacX1_20(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + IfacX2_20(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX3_35(ituvqminus1)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + IfacX3_20(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
 ! Building for Angular momentum Jp = 2
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,5) = facX*Tmp0(iTUVQ,2)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,6) = facX*Tmp0(iTUVQ,3)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,7) = facX*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,8) = facY*Tmp0(iTUVQ,3)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,9) = facY*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 20
      Tmp0(iTUVQ,10) = facZ*Tmp0(iTUVQ,4)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + IfacX1_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + IfacX1_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX1_35(ituvqminus1)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + IfacX1_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + IfacX2_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX2_35(ituvqminus1)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + IfacX2_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX3_35(ituvqminus1)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + IfacX3_20(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp0(iTUVplus1,2)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX1_35(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX2_35(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,10
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 11,20
      iTUVplus1 = TUVindexX3_35(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
     DO iTUVQ=1, 20
      DO iTUVP=1, 10
        Aux2(iTUVP,iTUVQ,iP) = Aux2(iTUVP,iTUVQ,iP) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimP=1, nPrimP
  ENDDO !iP = 1,nPrimQ*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP2Q3CtoASegQ
 subroutine TransferRecurrenceCPUP2Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  use AGC_CPU_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(   84,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   10,   35,nPrimP*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 35, 10)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
  real(realk) :: Tmp1( 36: 56,  2:  4)
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$OMP DO COLLAPSE(3) PRIVATE(iP,iTUVP,iTUVQ)
  DO iP = 1,nPrimP*nPasses
   DO iTUVQ=1, 35
    DO iTUVP=1, 10
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimD,iPrimB,&
!$OMP         Tmp0,&
!$OMP         Tmp1,&
!$OMP         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) 
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 35
      Tmp0(iTUVQ,1) = Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,2) = facX*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,3) = facY*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,4) = facZ*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  36, 56
      Tmp1(iTUVQ,2) = facX*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  36, 56
      Tmp1(iTUVQ,3) = facY*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  36, 56
      Tmp1(iTUVQ,4) = facZ*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_56(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1_35(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_56(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2_35(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_56(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3_35(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1_56(ituvqminus1)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + IfacX1_35(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX2_56(ituvqminus1)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + IfacX2_35(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX3_56(ituvqminus1)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + IfacX3_35(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX3_56(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX3_56(iTUVQ)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
 ! Building for Angular momentum Jp = 2
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,5) = facX*Tmp0(iTUVQ,2)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,6) = facX*Tmp0(iTUVQ,3)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,7) = facX*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,8) = facY*Tmp0(iTUVQ,3)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,9) = facY*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,10) = facZ*Tmp0(iTUVQ,4)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_56(ituvqminus1)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + IfacX1_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_56(ituvqminus1)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + IfacX1_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_56(ituvqminus1)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + IfacX1_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_56(ituvqminus1)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + IfacX2_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_56(ituvqminus1)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + IfacX2_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_56(ituvqminus1)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + IfacX3_35(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp0(iTUVplus1,2)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_56(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_56(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_56(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_56(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp1(iTUVplus1,4)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
     DO iTUVQ=1, 35
      DO iTUVP=1, 10
        Aux2(iTUVP,iTUVQ,iP) = Aux2(iTUVP,iTUVQ,iP) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimP=1, nPrimP
  ENDDO !iP = 1,nPrimQ*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP2Q4CtoASegQ
 subroutine TransferRecurrenceCPUP3Q4CtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Dexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  use AGC_CPU_OBS_TRParamMod
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Dexp(nPrimD),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(  120,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   20,   35,nPrimP*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
  !Local variables
  real(realk) :: Tmp0( 35, 20)
! Note that Tmp0 have the opposite order Tmp0(nTUVQ,nTUVP), than the Aux2
  real(realk) :: Tmp1( 36: 84,  2:  4)
  real(realk) :: Tmp2( 36: 56,  5: 10)
!  Note Tmp(nTUVQ,nTUVP) ordering different from Aux2 and the PtoQ routines 
  integer :: iPassP,iPrimP,iPrimQ,iPrimQP,IP,iTUVP,iTUVQ,iTUVplus1,ituvqminus1,iAtomA,iAtomB
  integer :: iPrimD,iPrimB
  real(realk),parameter :: D1=1.0E0_realk,D05=0.5E0_realk
  real(realk) :: Xab,Yab,Zab,Xcd,Ycd,Zcd,expP
  real(realk) :: expBX,expBY,expBZ
  real(realk) :: invexpP,inv2expP,facX,facY,facZ,qinvp
!$OMP DO COLLAPSE(3) PRIVATE(iP,iTUVP,iTUVQ)
  DO iP = 1,nPrimP*nPasses
   DO iTUVQ=1, 35
    DO iTUVP=1, 20
     Aux2(iTUVP,iTUVQ,iP) = 0.0E0_realk
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xab,Yab,Zab,Xcd,Ycd,Zcd,expP,&
!$OMP         iP,iPrimQ,iPrimP,iPassP,&
!$OMP         expBX,expBY,expBZ,&
!$OMP         iPrimD,iPrimB,&
!$OMP         Tmp0,&
!$OMP         Tmp1,&
!$OMP         Tmp2,&
!$OMP         invexpP,inv2expP,facX,facY,facZ,qinvp,iTUVQ,iTUVP,iTUVplus1) 
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
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    iPrimB = (iPrimP-1)/nPrimA+1                
    expBX = Bexp(iPrimB)*Xab
    expBY = Bexp(iPrimB)*Yab
    expBZ = Bexp(iPrimB)*Zab
     iPrimD = (iPrimQ-1)/nPrimC+1                
     facX = -(expBX+Dexp(iPrimD)*Xcd)*invexpP
     facY = -(expBY+Dexp(iPrimD)*Ycd)*invexpP
     facZ = -(expBZ+Dexp(iPrimD)*Zcd)*invexpP
     qinvp = -Qexp(iPrimQ)*invexpP
 ! Building for Angular momentum Jp = 0
     DO iTUVQ=1, 35
      Tmp0(iTUVQ,1) = Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     ENDDO
 ! Building for Angular momentum Jp = 1
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,2) = facX*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,3) = facY*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,4) = facZ*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  36, 84
      Tmp1(iTUVQ,2) = facX*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  36, 84
      Tmp1(iTUVQ,3) = facY*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  36, 84
      Tmp1(iTUVQ,4) = facZ*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1_56(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2_56(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3_56(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 21,56
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + IfacX1_56(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 21,56
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + IfacX2_56(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 21,56
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + IfacX3_56(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 36,84
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 36,84
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1,35
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 36,84
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
 ! Building for Angular momentum Jp = 2
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,5) = facX*Tmp0(iTUVQ,2)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,6) = facX*Tmp0(iTUVQ,3)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,7) = facX*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,8) = facY*Tmp0(iTUVQ,3)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,9) = facY*Tmp0(iTUVQ,4)
     enddo
     do iTUVQ = 1, 35
      Tmp0(iTUVQ,10) = facZ*Tmp0(iTUVQ,4)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,5) = facX*Tmp1(iTUVQ,2)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,6) = facX*Tmp1(iTUVQ,3)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,7) = facX*Tmp1(iTUVQ,4)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,8) = facY*Tmp1(iTUVQ,3)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,9) = facY*Tmp1(iTUVQ,4)
     enddo
     do iTUVQ =  36, 56
      Tmp2(iTUVQ,10) = facZ*Tmp1(iTUVQ,4)+ inv2expP*Aux(iTUVQ,iPrimQ,iPrimP,iPassP)
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + IfacX3_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp2(iTUVQ,5) = Tmp2(iTUVQ,5) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,2) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp2(iTUVQ,6) = Tmp2(iTUVQ,6) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp2(iTUVQ,7) = Tmp2(iTUVQ,7) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp2(iTUVQ,8) = Tmp2(iTUVQ,8) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,3) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp2(iTUVQ,9) = Tmp2(iTUVQ,9) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do ituvqminus1 = 21,35
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp2(iTUVQ,10) = Tmp2(iTUVQ,10) + IfacX3_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,4) 
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp0(iTUVplus1,2)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,5) = Tmp0(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp2(iTUVQ,5) = Tmp2(iTUVQ,5) + qinvp*Tmp1(iTUVplus1,2)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,6) = Tmp0(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp2(iTUVQ,6) = Tmp2(iTUVQ,6) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,7) = Tmp0(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp2(iTUVQ,7) = Tmp2(iTUVQ,7) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp0(iTUVplus1,3)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,8) = Tmp0(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp2(iTUVQ,8) = Tmp2(iTUVQ,8) + qinvp*Tmp1(iTUVplus1,3)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,9) = Tmp0(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp2(iTUVQ,9) = Tmp2(iTUVQ,9) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp0(iTUVplus1,4)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,10) = Tmp0(iTUVQ,10) + qinvp*Tmp1(iTUVplus1,4)
     enddo
     do iTUVQ = 36,56
      iTUVplus1 = TUVindexX3_84(iTUVQ)
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
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,11) = Tmp0(iTUVQ,11) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,5) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,12) = Tmp0(iTUVQ,12) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,5) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp0(iTUVQ,13) = Tmp0(iTUVQ,13) + IfacX3_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,5) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,14) = Tmp0(iTUVQ,14) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,8) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,15) = Tmp0(iTUVQ,15) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,9) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX1_84(ituvqminus1)
      Tmp0(iTUVQ,16) = Tmp0(iTUVQ,16) + IfacX1_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,10) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,17) = Tmp0(iTUVQ,17) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,8) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp0(iTUVQ,18) = Tmp0(iTUVQ,18) + IfacX3_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,8) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX2_84(ituvqminus1)
      Tmp0(iTUVQ,19) = Tmp0(iTUVQ,19) + IfacX2_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,10) 
     enddo
     do ituvqminus1 = 1,20
      iTUVQ = TUVindexX3_84(ituvqminus1)
      Tmp0(iTUVQ,20) = Tmp0(iTUVQ,20) + IfacX3_56(ituvqminus1)*inv2expP*Tmp0(ituvqminus1,10) 
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,11) = Tmp0(iTUVQ,11) + qinvp*Tmp0(iTUVplus1,5)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,11) = Tmp0(iTUVQ,11) + qinvp*Tmp2(iTUVplus1,5)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,12) = Tmp0(iTUVQ,12) + qinvp*Tmp0(iTUVplus1,5)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,12) = Tmp0(iTUVQ,12) + qinvp*Tmp2(iTUVplus1,5)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,13) = Tmp0(iTUVQ,13) + qinvp*Tmp0(iTUVplus1,5)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,13) = Tmp0(iTUVQ,13) + qinvp*Tmp2(iTUVplus1,5)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,14) = Tmp0(iTUVQ,14) + qinvp*Tmp0(iTUVplus1,8)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,14) = Tmp0(iTUVQ,14) + qinvp*Tmp2(iTUVplus1,8)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,15) = Tmp0(iTUVQ,15) + qinvp*Tmp0(iTUVplus1,9)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,15) = Tmp0(iTUVQ,15) + qinvp*Tmp2(iTUVplus1,9)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,16) = Tmp0(iTUVQ,16) + qinvp*Tmp0(iTUVplus1,10)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1_84(iTUVQ)
      Tmp0(iTUVQ,16) = Tmp0(iTUVQ,16) + qinvp*Tmp2(iTUVplus1,10)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,17) = Tmp0(iTUVQ,17) + qinvp*Tmp0(iTUVplus1,8)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,17) = Tmp0(iTUVQ,17) + qinvp*Tmp2(iTUVplus1,8)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,18) = Tmp0(iTUVQ,18) + qinvp*Tmp0(iTUVplus1,8)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,18) = Tmp0(iTUVQ,18) + qinvp*Tmp2(iTUVplus1,8)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,19) = Tmp0(iTUVQ,19) + qinvp*Tmp0(iTUVplus1,10)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2_84(iTUVQ)
      Tmp0(iTUVQ,19) = Tmp0(iTUVQ,19) + qinvp*Tmp2(iTUVplus1,10)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,20) = Tmp0(iTUVQ,20) + qinvp*Tmp0(iTUVplus1,10)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3_84(iTUVQ)
      Tmp0(iTUVQ,20) = Tmp0(iTUVQ,20) + qinvp*Tmp2(iTUVplus1,10)
     enddo
!    Warning Note Tmp0 have the opposite ordering so this is not that efficient. 
!    Hopefully Tmp0 is small enough that it can be in cache. 
     DO iTUVQ=1, 35
      DO iTUVP=1, 20
        Aux2(iTUVP,iTUVQ,iP) = Aux2(iTUVP,iTUVQ,iP) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimP=1, nPrimP
  ENDDO !iP = 1,nPrimQ*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP3Q4CtoASegQ
end module
