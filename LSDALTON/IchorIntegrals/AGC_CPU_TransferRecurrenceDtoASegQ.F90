MODULE AGC_CPU_OBS_TRMODDtoASegQ
 use IchorPrecisionMod
  
 CONTAINS
 subroutine TransferRecurrenceCPUP1Q2DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Cexp(nPrimC),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(   20,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(    4,   10,nPrimP*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
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
!$OMP         iPrimC,iPrimB,&
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
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpP
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpP
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpP
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
 end subroutine TransferRecurrenceCPUP1Q2DtoASegQ
 subroutine TransferRecurrenceCPUP1Q3DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Cexp(nPrimC),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(   35,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(    4,   20,nPrimP*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
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
!$OMP         iPrimC,iPrimB,&
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
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpP
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpP
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpP
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
 end subroutine TransferRecurrenceCPUP1Q3DtoASegQ
 subroutine TransferRecurrenceCPUP2Q3DtoASegQ(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & Pexp,Qexp,Pdistance12,Qdistance12,Cexp,Bexp,nPrimA,nPrimB,nPrimC,nPrimD,&
         & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,Aux,Aux2)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ,nPrimA,nPrimB,nPrimC,nPrimD,nAtomsA,nAtomsB,MaxPasses
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pdistance12(3,nAtomsA,nAtomsB),Qdistance12(3)
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: Cexp(nPrimC),Bexp(nPrimB)
  real(realk),intent(in) :: Aux(   56,nPrimQ,nPrimP,nPasses)
  real(realk),intent(inout) :: Aux2(   10,   20,nPrimP*nPasses)
!  real(realk),intent(inout) :: Aux2(nTUVP,nTUVQ,nPrimQ,nPrimP,nPasses)
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
!$OMP         iPrimC,iPrimB,&
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
     iPrimC = iPrimQ - ((iPrimQ-1)/nPrimC)*nPrimC
     facX = -(expBX-Cexp(iPrimC)*Xcd)*invexpP
     facY = -(expBY-Cexp(iPrimC)*Ycd)*invexpP
     facZ = -(expBZ-Cexp(iPrimC)*Zcd)*invexpP
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
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + IfacX1(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + IfacX2(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 1,10
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + IfacX3(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX1(ituvqminus1)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + IfacX1(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX2(ituvqminus1)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + IfacX2(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do ituvqminus1 = 11,20
      iTUVQ = TUVindexX3(ituvqminus1)
      Tmp1(iTUVQ,4) = Tmp1(iTUVQ,4) + IfacX3(ituvqminus1)*inv2expP*Aux(ituvqminus1,iPrimQ,iPrimP,iPassP) 
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp0(iTUVQ,2) = Tmp0(iTUVQ,2) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX1(iTUVQ)
      Tmp1(iTUVQ,2) = Tmp1(iTUVQ,2) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp0(iTUVQ,3) = Tmp0(iTUVQ,3) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX2(iTUVQ)
      Tmp1(iTUVQ,3) = Tmp1(iTUVQ,3) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 1,20
      iTUVplus1 = TUVindexX3(iTUVQ)
      Tmp0(iTUVQ,4) = Tmp0(iTUVQ,4) + qinvp*Aux(iTUVplus1,iPrimQ,iPrimP,iPassP)
     enddo
     do iTUVQ = 21,35
      iTUVplus1 = TUVindexX3(iTUVQ)
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
        Aux2(iTUVP,iTUVQ,iP) = Aux2(iTUVP,iTUVQ,iP) + Tmp0(iTUVQ,iTUVP)
      ENDDO
     ENDDO
   ENDDO !iPrimP=1, nPrimP
  ENDDO !iP = 1,nPrimQ*nPasses
!$OMP END DO
 end subroutine TransferRecurrenceCPUP2Q3DtoASegQ
end module
