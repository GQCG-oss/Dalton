MODULE AGC_GPU_OBS_VERTICALRECURRENCEMODBSegP
 use IchorPrecisionMod
  
 CONTAINS


subroutine VerticalRecurrenceGPUSegP1B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: TABFJW(0:4,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPassP,4)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  Integer :: iP,iPrimQ,iPrimP,iPassP,ipnt,iAtomA,iAtomB
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,alphaP,RJ000(0:1)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(realk),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D2=2.0E0_realk
  real(realk),PARAMETER :: D2JP36=  3.8000000000000000E+01_realk
  real(realk),parameter :: D1=1.0E0_realk,D03333=1.0E0_realk/3.0E0_realk
  real(realk),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk
  real(realk),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
  real(realk),PARAMETER :: SMALL = 1E-15_realk,D12 = 12.0E0_realk
  real(realk), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
  real(realk), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
  real(realk), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
  real(realk), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
  real(realk), parameter :: PI=3.14159265358979323846E0_realk
  real(realk), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
  real(realk), PARAMETER :: SQRPIH = SQRTPI/D2
  real(realk), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
  !ThetaAux(n,1,0,0) = Xpb*ThetaAux(n,0,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,0,0,0)
  !i = 0 last 2 term vanish
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP PRIVATE(iP) PRESENT(AUXarray)
  DO iP = 1,nPrimQ*nPassP
    AUXarray(iP,1)=0.0E0_realk
    AUXarray(iP,2)=0.0E0_realk
    AUXarray(iP,3)=0.0E0_realk
    AUXarray(iP,4)=0.0E0_realk
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         RWVAL,GVAL,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpP,&
!$ACC         PREF,&
!$ACC         TMP1,TMP2,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,TABFJW,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Pexp,Bcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPassP
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Bx = -Bcenter(1,iAtomB)
   By = -Bcenter(2,iAtomB)
   Bz = -Bcenter(3,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    Xpb = Pcent(1,iPrimP,iAtomA,iAtomB) + Bx
    Ypb = Pcent(2,iPrimP,iAtomA,iAtomB) + By
    Zpb = Pcent(3,iPrimP,iAtomA,iAtomB) + Bz
   Xpq = mPX + Qcent(1,iPrimQ)
   Ypq = mPY + Qcent(2,iPrimQ)
   Zpq = mPZ + Qcent(3,iPrimQ)
     alphaP = reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = alphaP*Xpq
     alphaYpq = alphaP*Ypq
     alphaZpq = alphaP*Zpq
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000(0)=TABFJW(0,IPNT)-TABFJW(1,IPNT)*WDIFF+TABFJW(2,IPNT)*W2+TABFJW(3,IPNT)*W3
      RJ000(1)=TABFJW(1,IPNT)-TABFJW(2,IPNT)*WDIFF+TABFJW(3,IPNT)*W2+TABFJW(4,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000(1) = RWVAL*(D05*RJ000(0)-REXPW)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RJ000(1) = RWVAL*PID4I*D05*RJ000(0)
     ENDIF
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)
     TMP1 = PREF*RJ000(0)
     TMP2 = PREF*RJ000(1)
     AUXarray(iP,1) = AUXarray(iP,1) + TMP1
     AUXarray(iP,2) = AUXarray(iP,2) + Xpb*TMP1 + alphaXpq*TMP2
     AUXarray(iP,3) = AUXarray(iP,3) + Ypb*TMP1 + alphaYpq*TMP2
     AUXarray(iP,4) = AUXarray(iP,4) + Zpb*TMP1 + alphaZpq*TMP2
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPassP
end subroutine VerticalRecurrenceGPUSegP1B

subroutine VerticalRecurrenceGPUSegP2B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 2)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPassP,   10)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: TMPAUXarray(    4)
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(   1)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:3)
  real(realk) :: TMParray2(  2:  4,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP PRIVATE(iP,iTUV) PRESENT(AUXarray)
  DO iP = 1,nPrimQ*nPassP
    DO iTUV=1,   10
     AUXarray(iP,iTUV)=0.0E0_realk
    ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpP,inv2expP,&
!$ACC         PREF,&
!$ACC         TMPAUXarray,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Pexp,Bcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPassP
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
   Bx = -Bcenter(1,iAtomB)
   By = -Bcenter(2,iAtomB)
   Bz = -Bcenter(3,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpb = Pcent(1,iPrimP,iAtomA,iAtomB) + Bx
    Ypb = Pcent(2,iPrimP,iAtomA,iAtomB) + By
    Zpb = Pcent(3,iPrimP,iAtomA,iAtomB) + Bz
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)
     TMPAuxarray(1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     TMPAuxArray(2) = Xpb*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypb*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpb*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(1) + alphaP*TmpArray1(1,2))
!$ACC LOOP SEQ
     do iTUV = 1,    4
      AuxArray(iP,iTUV) = AuxArray(iP,iTUV) + TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,5) = AuxArray(iP,5) + Xpb*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = AuxArray(iP,6) + Xpb*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = AuxArray(iP,7) + Xpb*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = AuxArray(iP,8) + Ypb*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = AuxArray(iP,9) + Ypb*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = AuxArray(iP,10) + Zpb*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUSegP3B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 3)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPassP,   20)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: TMPAUXarray(   10)
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(   3)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:4)
  real(realk) :: TMParray2(  2:  4,2:3)
  real(realk) :: TMParray3(  5: 10,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP PRIVATE(iP,iTUV) PRESENT(AUXarray)
  DO iP = 1,nPrimQ*nPassP
    DO iTUV=1,   20
     AUXarray(iP,iTUV)=0.0E0_realk
    ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpP,inv2expP,&
!$ACC         PREF,&
!$ACC         TMPAUXarray,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Pexp,Bcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPassP
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
   Bx = -Bcenter(1,iAtomB)
   By = -Bcenter(2,iAtomB)
   Bz = -Bcenter(3,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpb = Pcent(1,iPrimP,iAtomA,iAtomB) + Bx
    Ypb = Pcent(2,iPrimP,iAtomA,iAtomB) + By
    Zpb = Pcent(3,iPrimP,iAtomA,iAtomB) + Bz
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)
     TMPAuxarray(1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 3)
     TMPAuxArray(2) = Xpb*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypb*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpb*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpb*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypb*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpb*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expP*(TMPAuxArray(1) + alphaP*TmpArray1(1,2))
     TMPAuxArray(5) = Xpb*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     TMPAuxArray(6) = Xpb*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     TMPAuxArray(7) = Xpb*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     TMPAuxArray(8) = Ypb*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     TMPAuxArray(9) = Ypb*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     TMPAuxArray(10) = Zpb*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,2) + alphaP*TmpArray1(1,3))
     tmpArray3(5,2) = Xpb*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xpb*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xpb*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Ypb*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Ypb*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zpb*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TMPAuxArray(2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(4) + alphaP*TmpArray2(4,2))
!$ACC LOOP SEQ
     do iTUV = 1,   10
      AuxArray(iP,iTUV) = AuxArray(iP,iTUV) + TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,11) = AuxArray(iP,11) + Xpb*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(iP,12) = AuxArray(iP,12) + Ypb*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     AuxArray(iP,13) = AuxArray(iP,13) + Zpb*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     AuxArray(iP,14) = AuxArray(iP,14) + Xpb*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     AuxArray(iP,15) = AuxArray(iP,15) + Xpb*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     AuxArray(iP,16) = AuxArray(iP,16) + Xpb*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     AuxArray(iP,17) = AuxArray(iP,17) + Ypb*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(iP,18) = AuxArray(iP,18) + Zpb*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     AuxArray(iP,19) = AuxArray(iP,19) + Ypb*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     AuxArray(iP,20) = AuxArray(iP,20) + Zpb*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPassP
 end subroutine
end module
