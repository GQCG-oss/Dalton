MODULE AGC_GPU_OBS_VERTICALRECURRENCEMODBGen
 use IchorPrecisionMod
  
 CONTAINS


subroutine VerticalRecurrenceGPUGen1B(nPassP,nPrimP,nPrimQ,&
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
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPrimP*nPassP,4)
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
  DO iP = 1,nPrimQ*nPrimP*nPassP
   iPrimQ = mod(IP-1,nPrimQ)+1
   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
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
     AUXarray(iP,1) = TMP1
     AUXarray(iP,2) = Xpb*TMP1 + alphaXpq*TMP2
     AUXarray(iP,3) = Ypb*TMP1 + alphaYpq*TMP2
     AUXarray(iP,4) = Zpb*TMP1 + alphaZpq*TMP2
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
end subroutine VerticalRecurrenceGPUGen1B

subroutine VerticalRecurrenceGPUGen2B(nPassP,nPrimP,nPrimQ,&
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
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPrimP*nPassP,   10)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
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
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpP,inv2expP,&
!$ACC         PREF,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Pexp,Bcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPrimP*nPassP
   iPrimQ = mod(IP-1,nPrimQ)+1
   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
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
     Auxarray(IP,1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     AuxArray(iP,2) = Xpb*AuxArray(iP,1) + alphaXpq*TmpArray1(1,2)
     AuxArray(iP,3) = Ypb*AuxArray(iP,1) + alphaYpq*TmpArray1(1,2)
     AuxArray(iP,4) = Zpb*AuxArray(iP,1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expP*(AuxArray(IP,1) + alphaP*TmpArray1(1,2))
     AuxArray(iP,5) = Xpb*AuxArray(iP,2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = Xpb*AuxArray(iP,3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = Xpb*AuxArray(iP,4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = Ypb*AuxArray(iP,3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = Ypb*AuxArray(iP,4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = Zpb*AuxArray(iP,4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen3B(nPassP,nPrimP,nPrimQ,&
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
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPrimP*nPassP,   20)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
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
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpP,inv2expP,&
!$ACC         PREF,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Pexp,Bcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPrimP*nPassP
   iPrimQ = mod(IP-1,nPrimQ)+1
   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
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
     Auxarray(IP,1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 3)
     AuxArray(iP,2) = Xpb*AuxArray(iP,1) + alphaXpq*TmpArray1(1,2)
     AuxArray(iP,3) = Ypb*AuxArray(iP,1) + alphaYpq*TmpArray1(1,2)
     AuxArray(iP,4) = Zpb*AuxArray(iP,1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpb*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypb*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpb*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expP*(AuxArray(IP,1) + alphaP*TmpArray1(1,2))
     AuxArray(iP,5) = Xpb*AuxArray(iP,2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = Xpb*AuxArray(iP,3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = Xpb*AuxArray(iP,4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = Ypb*AuxArray(iP,3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = Ypb*AuxArray(iP,4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = Zpb*AuxArray(iP,4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,2) + alphaP*TmpArray1(1,3))
     tmpArray3(5,2) = Xpb*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xpb*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xpb*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Ypb*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Ypb*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zpb*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(AuxArray(IP,2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(IP,3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(IP,4) + alphaP*TmpArray2(4,2))
     AuxArray(iP,11) = Xpb*AuxArray(iP,5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(iP,12) = Ypb*AuxArray(iP,5) + alphaYpq*TmpArray3(5,2)
     AuxArray(iP,13) = Zpb*AuxArray(iP,5) + alphaZpq*TmpArray3(5,2)
     AuxArray(iP,14) = Xpb*AuxArray(iP,8) + alphaXpq*TmpArray3(8,2)
     AuxArray(iP,15) = Xpb*AuxArray(iP,9) + alphaXpq*TmpArray3(9,2)
     AuxArray(iP,16) = Xpb*AuxArray(iP,10) + alphaXpq*TmpArray3(10,2)
     AuxArray(iP,17) = Ypb*AuxArray(iP,8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(iP,18) = Zpb*AuxArray(iP,8) + alphaZpq*TmpArray3(8,2)
     AuxArray(iP,19) = Ypb*AuxArray(iP,10) + alphaYpq*TmpArray3(10,2)
     AuxArray(iP,20) = Zpb*AuxArray(iP,10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen4B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 4)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPrimP*nPassP,   35)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(   6)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:5)
  real(realk) :: TMParray2(  2:  4,2:4)
  real(realk) :: TMParray3(  5: 10,2:3)
  real(realk) :: TMParray4( 11: 20,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpP,inv2expP,&
!$ACC         PREF,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TMParray4,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Pexp,Bcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPrimP*nPassP
   iPrimQ = mod(IP-1,nPrimQ)+1
   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
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
     Auxarray(IP,1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 3)
     TMParray1(1, 5) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 4)
     AuxArray(iP,2) = Xpb*AuxArray(iP,1) + alphaXpq*TmpArray1(1,2)
     AuxArray(iP,3) = Ypb*AuxArray(iP,1) + alphaYpq*TmpArray1(1,2)
     AuxArray(iP,4) = Zpb*AuxArray(iP,1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpb*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypb*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpb*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xpb*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Ypb*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zpb*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     TwoTerms(1) = inv2expP*(AuxArray(IP,1) + alphaP*TmpArray1(1,2))
     AuxArray(iP,5) = Xpb*AuxArray(iP,2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = Xpb*AuxArray(iP,3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = Xpb*AuxArray(iP,4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = Ypb*AuxArray(iP,3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = Ypb*AuxArray(iP,4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = Zpb*AuxArray(iP,4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,2) + alphaP*TmpArray1(1,3))
     tmpArray3(5,2) = Xpb*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xpb*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xpb*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Ypb*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Ypb*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zpb*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,3) + alphaP*TmpArray1(1,4))
     tmpArray3(5,3) = Xpb*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xpb*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xpb*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Ypb*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Ypb*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zpb*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(AuxArray(IP,2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(IP,3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(IP,4) + alphaP*TmpArray2(4,2))
     AuxArray(iP,11) = Xpb*AuxArray(iP,5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(iP,12) = Ypb*AuxArray(iP,5) + alphaYpq*TmpArray3(5,2)
     AuxArray(iP,13) = Zpb*AuxArray(iP,5) + alphaZpq*TmpArray3(5,2)
     AuxArray(iP,14) = Xpb*AuxArray(iP,8) + alphaXpq*TmpArray3(8,2)
     AuxArray(iP,15) = Xpb*AuxArray(iP,9) + alphaXpq*TmpArray3(9,2)
     AuxArray(iP,16) = Xpb*AuxArray(iP,10) + alphaXpq*TmpArray3(10,2)
     AuxArray(iP,17) = Ypb*AuxArray(iP,8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(iP,18) = Zpb*AuxArray(iP,8) + alphaZpq*TmpArray3(8,2)
     AuxArray(iP,19) = Ypb*AuxArray(iP,10) + alphaYpq*TmpArray3(10,2)
     AuxArray(iP,20) = Zpb*AuxArray(iP,10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,2) + alphaP*TmpArray2(2,3))
     TwoTerms(2) = inv2expP*(TmpArray2(3,2) + alphaP*TmpArray2(3,3))
     TwoTerms(3) = inv2expP*(TmpArray2(4,2) + alphaP*TmpArray2(4,3))
     tmpArray4(11,2) = Xpb*tmpArray3(5,2) + alphaXpq*TmpArray3(5,3) + 2*TwoTerms(1)
     tmpArray4(12,2) = Ypb*tmpArray3(5,2) + alphaYpq*TmpArray3(5,3)
     tmpArray4(13,2) = Zpb*tmpArray3(5,2) + alphaZpq*TmpArray3(5,3)
     tmpArray4(14,2) = Xpb*tmpArray3(8,2) + alphaXpq*TmpArray3(8,3)
     tmpArray4(15,2) = Xpb*tmpArray3(9,2) + alphaXpq*TmpArray3(9,3)
     tmpArray4(16,2) = Xpb*tmpArray3(10,2) + alphaXpq*TmpArray3(10,3)
     tmpArray4(17,2) = Ypb*tmpArray3(8,2) + alphaYpq*TmpArray3(8,3) + 2*TwoTerms(2)
     tmpArray4(18,2) = Zpb*tmpArray3(8,2) + alphaZpq*TmpArray3(8,3)
     tmpArray4(19,2) = Ypb*tmpArray3(10,2) + alphaYpq*TmpArray3(10,3)
     tmpArray4(20,2) = Zpb*tmpArray3(10,2) + alphaZpq*TmpArray3(10,3) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(AuxArray(IP,5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(AuxArray(IP,8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(AuxArray(IP,10) + alphaP*TmpArray3(10,2))
     AuxArray(iP,21) = Xpb*AuxArray(iP,11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(iP,22) = Ypb*AuxArray(iP,11) + alphaYpq*TmpArray4(11,2)
     AuxArray(iP,23) = Zpb*AuxArray(iP,11) + alphaZpq*TmpArray4(11,2)
     AuxArray(iP,24) = Xpb*AuxArray(iP,14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(iP,25) = Ypb*AuxArray(iP,13) + alphaYpq*TmpArray4(13,2)
     AuxArray(iP,26) = Xpb*AuxArray(iP,16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(iP,27) = Xpb*AuxArray(iP,17) + alphaXpq*TmpArray4(17,2)
     AuxArray(iP,28) = Xpb*AuxArray(iP,18) + alphaXpq*TmpArray4(18,2)
     AuxArray(iP,29) = Xpb*AuxArray(iP,19) + alphaXpq*TmpArray4(19,2)
     AuxArray(iP,30) = Xpb*AuxArray(iP,20) + alphaXpq*TmpArray4(20,2)
     AuxArray(iP,31) = Ypb*AuxArray(iP,17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(iP,32) = Zpb*AuxArray(iP,17) + alphaZpq*TmpArray4(17,2)
     AuxArray(iP,33) = Ypb*AuxArray(iP,19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(iP,34) = Ypb*AuxArray(iP,20) + alphaYpq*TmpArray4(20,2)
     AuxArray(iP,35) = Zpb*AuxArray(iP,20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen5B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 5)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPrimP*nPassP,   56)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(  10)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:6)
  real(realk) :: TMParray2(  2:  4,2:5)
  real(realk) :: TMParray3(  5: 10,2:4)
  real(realk) :: TMParray4( 11: 20,2:3)
  real(realk) :: TMParray5( 21: 35,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpP,inv2expP,&
!$ACC         PREF,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TMParray4,&
!$ACC         TMParray5,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Pexp,Bcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPrimP*nPassP
   iPrimQ = mod(IP-1,nPrimQ)+1
   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
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
     Auxarray(IP,1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 3)
     TMParray1(1, 5) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 4)
     TMParray1(1, 6) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 5)
     AuxArray(iP,2) = Xpb*AuxArray(iP,1) + alphaXpq*TmpArray1(1,2)
     AuxArray(iP,3) = Ypb*AuxArray(iP,1) + alphaYpq*TmpArray1(1,2)
     AuxArray(iP,4) = Zpb*AuxArray(iP,1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpb*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypb*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpb*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xpb*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Ypb*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zpb*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     tmpArray2(2,5) = Xpb*tmpArray1(1,5) + alphaXpq*TmpArray1(1,6)
     tmpArray2(3,5) = Ypb*tmpArray1(1,5) + alphaYpq*TmpArray1(1,6)
     tmpArray2(4,5) = Zpb*tmpArray1(1,5) + alphaZpq*TmpArray1(1,6)
     TwoTerms(1) = inv2expP*(AuxArray(IP,1) + alphaP*TmpArray1(1,2))
     AuxArray(iP,5) = Xpb*AuxArray(iP,2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = Xpb*AuxArray(iP,3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = Xpb*AuxArray(iP,4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = Ypb*AuxArray(iP,3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = Ypb*AuxArray(iP,4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = Zpb*AuxArray(iP,4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,2) + alphaP*TmpArray1(1,3))
     tmpArray3(5,2) = Xpb*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xpb*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xpb*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Ypb*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Ypb*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zpb*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,3) + alphaP*TmpArray1(1,4))
     tmpArray3(5,3) = Xpb*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xpb*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xpb*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Ypb*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Ypb*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zpb*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,4) + alphaP*TmpArray1(1,5))
     tmpArray3(5,4) = Xpb*tmpArray2(2,4) + alphaXpq*TmpArray2(2,5) + TwoTerms(1)
     tmpArray3(6,4) = Xpb*tmpArray2(3,4) + alphaXpq*TmpArray2(3,5)
     tmpArray3(7,4) = Xpb*tmpArray2(4,4) + alphaXpq*TmpArray2(4,5)
     tmpArray3(8,4) = Ypb*tmpArray2(3,4) + alphaYpq*TmpArray2(3,5) + TwoTerms(1)
     tmpArray3(9,4) = Ypb*tmpArray2(4,4) + alphaYpq*TmpArray2(4,5)
     tmpArray3(10,4) = Zpb*tmpArray2(4,4) + alphaZpq*TmpArray2(4,5) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(AuxArray(IP,2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(IP,3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(IP,4) + alphaP*TmpArray2(4,2))
     AuxArray(iP,11) = Xpb*AuxArray(iP,5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(iP,12) = Ypb*AuxArray(iP,5) + alphaYpq*TmpArray3(5,2)
     AuxArray(iP,13) = Zpb*AuxArray(iP,5) + alphaZpq*TmpArray3(5,2)
     AuxArray(iP,14) = Xpb*AuxArray(iP,8) + alphaXpq*TmpArray3(8,2)
     AuxArray(iP,15) = Xpb*AuxArray(iP,9) + alphaXpq*TmpArray3(9,2)
     AuxArray(iP,16) = Xpb*AuxArray(iP,10) + alphaXpq*TmpArray3(10,2)
     AuxArray(iP,17) = Ypb*AuxArray(iP,8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(iP,18) = Zpb*AuxArray(iP,8) + alphaZpq*TmpArray3(8,2)
     AuxArray(iP,19) = Ypb*AuxArray(iP,10) + alphaYpq*TmpArray3(10,2)
     AuxArray(iP,20) = Zpb*AuxArray(iP,10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,2) + alphaP*TmpArray2(2,3))
     TwoTerms(2) = inv2expP*(TmpArray2(3,2) + alphaP*TmpArray2(3,3))
     TwoTerms(3) = inv2expP*(TmpArray2(4,2) + alphaP*TmpArray2(4,3))
     tmpArray4(11,2) = Xpb*tmpArray3(5,2) + alphaXpq*TmpArray3(5,3) + 2*TwoTerms(1)
     tmpArray4(12,2) = Ypb*tmpArray3(5,2) + alphaYpq*TmpArray3(5,3)
     tmpArray4(13,2) = Zpb*tmpArray3(5,2) + alphaZpq*TmpArray3(5,3)
     tmpArray4(14,2) = Xpb*tmpArray3(8,2) + alphaXpq*TmpArray3(8,3)
     tmpArray4(15,2) = Xpb*tmpArray3(9,2) + alphaXpq*TmpArray3(9,3)
     tmpArray4(16,2) = Xpb*tmpArray3(10,2) + alphaXpq*TmpArray3(10,3)
     tmpArray4(17,2) = Ypb*tmpArray3(8,2) + alphaYpq*TmpArray3(8,3) + 2*TwoTerms(2)
     tmpArray4(18,2) = Zpb*tmpArray3(8,2) + alphaZpq*TmpArray3(8,3)
     tmpArray4(19,2) = Ypb*tmpArray3(10,2) + alphaYpq*TmpArray3(10,3)
     tmpArray4(20,2) = Zpb*tmpArray3(10,2) + alphaZpq*TmpArray3(10,3) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,3) + alphaP*TmpArray2(2,4))
     TwoTerms(2) = inv2expP*(TmpArray2(3,3) + alphaP*TmpArray2(3,4))
     TwoTerms(3) = inv2expP*(TmpArray2(4,3) + alphaP*TmpArray2(4,4))
     tmpArray4(11,3) = Xpb*tmpArray3(5,3) + alphaXpq*TmpArray3(5,4) + 2*TwoTerms(1)
     tmpArray4(12,3) = Ypb*tmpArray3(5,3) + alphaYpq*TmpArray3(5,4)
     tmpArray4(13,3) = Zpb*tmpArray3(5,3) + alphaZpq*TmpArray3(5,4)
     tmpArray4(14,3) = Xpb*tmpArray3(8,3) + alphaXpq*TmpArray3(8,4)
     tmpArray4(15,3) = Xpb*tmpArray3(9,3) + alphaXpq*TmpArray3(9,4)
     tmpArray4(16,3) = Xpb*tmpArray3(10,3) + alphaXpq*TmpArray3(10,4)
     tmpArray4(17,3) = Ypb*tmpArray3(8,3) + alphaYpq*TmpArray3(8,4) + 2*TwoTerms(2)
     tmpArray4(18,3) = Zpb*tmpArray3(8,3) + alphaZpq*TmpArray3(8,4)
     tmpArray4(19,3) = Ypb*tmpArray3(10,3) + alphaYpq*TmpArray3(10,4)
     tmpArray4(20,3) = Zpb*tmpArray3(10,3) + alphaZpq*TmpArray3(10,4) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(AuxArray(IP,5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(AuxArray(IP,8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(AuxArray(IP,10) + alphaP*TmpArray3(10,2))
     AuxArray(iP,21) = Xpb*AuxArray(iP,11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(iP,22) = Ypb*AuxArray(iP,11) + alphaYpq*TmpArray4(11,2)
     AuxArray(iP,23) = Zpb*AuxArray(iP,11) + alphaZpq*TmpArray4(11,2)
     AuxArray(iP,24) = Xpb*AuxArray(iP,14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(iP,25) = Ypb*AuxArray(iP,13) + alphaYpq*TmpArray4(13,2)
     AuxArray(iP,26) = Xpb*AuxArray(iP,16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(iP,27) = Xpb*AuxArray(iP,17) + alphaXpq*TmpArray4(17,2)
     AuxArray(iP,28) = Xpb*AuxArray(iP,18) + alphaXpq*TmpArray4(18,2)
     AuxArray(iP,29) = Xpb*AuxArray(iP,19) + alphaXpq*TmpArray4(19,2)
     AuxArray(iP,30) = Xpb*AuxArray(iP,20) + alphaXpq*TmpArray4(20,2)
     AuxArray(iP,31) = Ypb*AuxArray(iP,17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(iP,32) = Zpb*AuxArray(iP,17) + alphaZpq*TmpArray4(17,2)
     AuxArray(iP,33) = Ypb*AuxArray(iP,19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(iP,34) = Ypb*AuxArray(iP,20) + alphaYpq*TmpArray4(20,2)
     AuxArray(iP,35) = Zpb*AuxArray(iP,20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,2) + alphaP*TmpArray3(5,3))
     TwoTerms(2) = inv2expP*(TmpArray3(8,2) + alphaP*TmpArray3(8,3))
     TwoTerms(3) = inv2expP*(TmpArray3(10,2) + alphaP*TmpArray3(10,3))
     tmpArray5(21,2) = Xpb*tmpArray4(11,2) + alphaXpq*TmpArray4(11,3) + 3*TwoTerms(1)
     tmpArray5(22,2) = Ypb*tmpArray4(11,2) + alphaYpq*TmpArray4(11,3)
     tmpArray5(23,2) = Zpb*tmpArray4(11,2) + alphaZpq*TmpArray4(11,3)
     tmpArray5(24,2) = Xpb*tmpArray4(14,2) + alphaXpq*TmpArray4(14,3) + TwoTerms(2)
     tmpArray5(25,2) = Ypb*tmpArray4(13,2) + alphaYpq*TmpArray4(13,3)
     tmpArray5(26,2) = Xpb*tmpArray4(16,2) + alphaXpq*TmpArray4(16,3) + TwoTerms(3)
     tmpArray5(27,2) = Xpb*tmpArray4(17,2) + alphaXpq*TmpArray4(17,3)
     tmpArray5(28,2) = Xpb*tmpArray4(18,2) + alphaXpq*TmpArray4(18,3)
     tmpArray5(29,2) = Xpb*tmpArray4(19,2) + alphaXpq*TmpArray4(19,3)
     tmpArray5(30,2) = Xpb*tmpArray4(20,2) + alphaXpq*TmpArray4(20,3)
     tmpArray5(31,2) = Ypb*tmpArray4(17,2) + alphaYpq*TmpArray4(17,3) + 3*TwoTerms(2)
     tmpArray5(32,2) = Zpb*tmpArray4(17,2) + alphaZpq*TmpArray4(17,3)
     tmpArray5(33,2) = Ypb*tmpArray4(19,2) + alphaYpq*TmpArray4(19,3) + TwoTerms(3)
     tmpArray5(34,2) = Ypb*tmpArray4(20,2) + alphaYpq*TmpArray4(20,3)
     tmpArray5(35,2) = Zpb*tmpArray4(20,2) + alphaZpq*TmpArray4(20,3) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(AuxArray(IP,11) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(AuxArray(IP,14) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(AuxArray(IP,16) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(AuxArray(IP,17) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(AuxArray(IP,19) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(AuxArray(IP,20) + alphaP*TmpArray4(20,2))
     AuxArray(iP,36) = Xpb*AuxArray(iP,21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(iP,37) = Ypb*AuxArray(iP,21) + alphaYpq*TmpArray5(21,2)
     AuxArray(iP,38) = Zpb*AuxArray(iP,21) + alphaZpq*TmpArray5(21,2)
     AuxArray(iP,39) = Xpb*AuxArray(iP,24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(iP,40) = Ypb*AuxArray(iP,23) + alphaYpq*TmpArray5(23,2)
     AuxArray(iP,41) = Xpb*AuxArray(iP,26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(iP,42) = Xpb*AuxArray(iP,27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(iP,43) = Zpb*AuxArray(iP,24) + alphaZpq*TmpArray5(24,2)
     AuxArray(iP,44) = Ypb*AuxArray(iP,26) + alphaYpq*TmpArray5(26,2)
     AuxArray(iP,45) = Xpb*AuxArray(iP,30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(iP,46) = Xpb*AuxArray(iP,31) + alphaXpq*TmpArray5(31,2)
     AuxArray(iP,47) = Xpb*AuxArray(iP,32) + alphaXpq*TmpArray5(32,2)
     AuxArray(iP,48) = Xpb*AuxArray(iP,33) + alphaXpq*TmpArray5(33,2)
     AuxArray(iP,49) = Xpb*AuxArray(iP,34) + alphaXpq*TmpArray5(34,2)
     AuxArray(iP,50) = Xpb*AuxArray(iP,35) + alphaXpq*TmpArray5(35,2)
     AuxArray(iP,51) = Ypb*AuxArray(iP,31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(iP,52) = Zpb*AuxArray(iP,31) + alphaZpq*TmpArray5(31,2)
     AuxArray(iP,53) = Ypb*AuxArray(iP,33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(iP,54) = Ypb*AuxArray(iP,34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(iP,55) = Ypb*AuxArray(iP,35) + alphaYpq*TmpArray5(35,2)
     AuxArray(iP,56) = Zpb*AuxArray(iP,35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen6B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 6)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPrimP*nPassP,   84)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(  15)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:7)
  real(realk) :: TMParray2(  2:  4,2:6)
  real(realk) :: TMParray3(  5: 10,2:5)
  real(realk) :: TMParray4( 11: 20,2:4)
  real(realk) :: TMParray5( 21: 35,2:3)
  real(realk) :: TMParray6( 36: 56,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpP,inv2expP,&
!$ACC         PREF,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TMParray4,&
!$ACC         TMParray5,&
!$ACC         TMParray6,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Pexp,Bcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPrimP*nPassP
   iPrimQ = mod(IP-1,nPrimQ)+1
   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
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
     Auxarray(IP,1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 3)
     TMParray1(1, 5) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 4)
     TMParray1(1, 6) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 5)
     TMParray1(1, 7) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 6)
     AuxArray(iP,2) = Xpb*AuxArray(iP,1) + alphaXpq*TmpArray1(1,2)
     AuxArray(iP,3) = Ypb*AuxArray(iP,1) + alphaYpq*TmpArray1(1,2)
     AuxArray(iP,4) = Zpb*AuxArray(iP,1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpb*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypb*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpb*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xpb*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Ypb*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zpb*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     tmpArray2(2,5) = Xpb*tmpArray1(1,5) + alphaXpq*TmpArray1(1,6)
     tmpArray2(3,5) = Ypb*tmpArray1(1,5) + alphaYpq*TmpArray1(1,6)
     tmpArray2(4,5) = Zpb*tmpArray1(1,5) + alphaZpq*TmpArray1(1,6)
     tmpArray2(2,6) = Xpb*tmpArray1(1,6) + alphaXpq*TmpArray1(1,7)
     tmpArray2(3,6) = Ypb*tmpArray1(1,6) + alphaYpq*TmpArray1(1,7)
     tmpArray2(4,6) = Zpb*tmpArray1(1,6) + alphaZpq*TmpArray1(1,7)
     TwoTerms(1) = inv2expP*(AuxArray(IP,1) + alphaP*TmpArray1(1,2))
     AuxArray(iP,5) = Xpb*AuxArray(iP,2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = Xpb*AuxArray(iP,3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = Xpb*AuxArray(iP,4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = Ypb*AuxArray(iP,3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = Ypb*AuxArray(iP,4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = Zpb*AuxArray(iP,4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,2) + alphaP*TmpArray1(1,3))
     tmpArray3(5,2) = Xpb*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xpb*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xpb*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Ypb*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Ypb*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zpb*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,3) + alphaP*TmpArray1(1,4))
     tmpArray3(5,3) = Xpb*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xpb*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xpb*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Ypb*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Ypb*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zpb*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,4) + alphaP*TmpArray1(1,5))
     tmpArray3(5,4) = Xpb*tmpArray2(2,4) + alphaXpq*TmpArray2(2,5) + TwoTerms(1)
     tmpArray3(6,4) = Xpb*tmpArray2(3,4) + alphaXpq*TmpArray2(3,5)
     tmpArray3(7,4) = Xpb*tmpArray2(4,4) + alphaXpq*TmpArray2(4,5)
     tmpArray3(8,4) = Ypb*tmpArray2(3,4) + alphaYpq*TmpArray2(3,5) + TwoTerms(1)
     tmpArray3(9,4) = Ypb*tmpArray2(4,4) + alphaYpq*TmpArray2(4,5)
     tmpArray3(10,4) = Zpb*tmpArray2(4,4) + alphaZpq*TmpArray2(4,5) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,5) + alphaP*TmpArray1(1,6))
     tmpArray3(5,5) = Xpb*tmpArray2(2,5) + alphaXpq*TmpArray2(2,6) + TwoTerms(1)
     tmpArray3(6,5) = Xpb*tmpArray2(3,5) + alphaXpq*TmpArray2(3,6)
     tmpArray3(7,5) = Xpb*tmpArray2(4,5) + alphaXpq*TmpArray2(4,6)
     tmpArray3(8,5) = Ypb*tmpArray2(3,5) + alphaYpq*TmpArray2(3,6) + TwoTerms(1)
     tmpArray3(9,5) = Ypb*tmpArray2(4,5) + alphaYpq*TmpArray2(4,6)
     tmpArray3(10,5) = Zpb*tmpArray2(4,5) + alphaZpq*TmpArray2(4,6) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(AuxArray(IP,2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(IP,3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(IP,4) + alphaP*TmpArray2(4,2))
     AuxArray(iP,11) = Xpb*AuxArray(iP,5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(iP,12) = Ypb*AuxArray(iP,5) + alphaYpq*TmpArray3(5,2)
     AuxArray(iP,13) = Zpb*AuxArray(iP,5) + alphaZpq*TmpArray3(5,2)
     AuxArray(iP,14) = Xpb*AuxArray(iP,8) + alphaXpq*TmpArray3(8,2)
     AuxArray(iP,15) = Xpb*AuxArray(iP,9) + alphaXpq*TmpArray3(9,2)
     AuxArray(iP,16) = Xpb*AuxArray(iP,10) + alphaXpq*TmpArray3(10,2)
     AuxArray(iP,17) = Ypb*AuxArray(iP,8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(iP,18) = Zpb*AuxArray(iP,8) + alphaZpq*TmpArray3(8,2)
     AuxArray(iP,19) = Ypb*AuxArray(iP,10) + alphaYpq*TmpArray3(10,2)
     AuxArray(iP,20) = Zpb*AuxArray(iP,10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,2) + alphaP*TmpArray2(2,3))
     TwoTerms(2) = inv2expP*(TmpArray2(3,2) + alphaP*TmpArray2(3,3))
     TwoTerms(3) = inv2expP*(TmpArray2(4,2) + alphaP*TmpArray2(4,3))
     tmpArray4(11,2) = Xpb*tmpArray3(5,2) + alphaXpq*TmpArray3(5,3) + 2*TwoTerms(1)
     tmpArray4(12,2) = Ypb*tmpArray3(5,2) + alphaYpq*TmpArray3(5,3)
     tmpArray4(13,2) = Zpb*tmpArray3(5,2) + alphaZpq*TmpArray3(5,3)
     tmpArray4(14,2) = Xpb*tmpArray3(8,2) + alphaXpq*TmpArray3(8,3)
     tmpArray4(15,2) = Xpb*tmpArray3(9,2) + alphaXpq*TmpArray3(9,3)
     tmpArray4(16,2) = Xpb*tmpArray3(10,2) + alphaXpq*TmpArray3(10,3)
     tmpArray4(17,2) = Ypb*tmpArray3(8,2) + alphaYpq*TmpArray3(8,3) + 2*TwoTerms(2)
     tmpArray4(18,2) = Zpb*tmpArray3(8,2) + alphaZpq*TmpArray3(8,3)
     tmpArray4(19,2) = Ypb*tmpArray3(10,2) + alphaYpq*TmpArray3(10,3)
     tmpArray4(20,2) = Zpb*tmpArray3(10,2) + alphaZpq*TmpArray3(10,3) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,3) + alphaP*TmpArray2(2,4))
     TwoTerms(2) = inv2expP*(TmpArray2(3,3) + alphaP*TmpArray2(3,4))
     TwoTerms(3) = inv2expP*(TmpArray2(4,3) + alphaP*TmpArray2(4,4))
     tmpArray4(11,3) = Xpb*tmpArray3(5,3) + alphaXpq*TmpArray3(5,4) + 2*TwoTerms(1)
     tmpArray4(12,3) = Ypb*tmpArray3(5,3) + alphaYpq*TmpArray3(5,4)
     tmpArray4(13,3) = Zpb*tmpArray3(5,3) + alphaZpq*TmpArray3(5,4)
     tmpArray4(14,3) = Xpb*tmpArray3(8,3) + alphaXpq*TmpArray3(8,4)
     tmpArray4(15,3) = Xpb*tmpArray3(9,3) + alphaXpq*TmpArray3(9,4)
     tmpArray4(16,3) = Xpb*tmpArray3(10,3) + alphaXpq*TmpArray3(10,4)
     tmpArray4(17,3) = Ypb*tmpArray3(8,3) + alphaYpq*TmpArray3(8,4) + 2*TwoTerms(2)
     tmpArray4(18,3) = Zpb*tmpArray3(8,3) + alphaZpq*TmpArray3(8,4)
     tmpArray4(19,3) = Ypb*tmpArray3(10,3) + alphaYpq*TmpArray3(10,4)
     tmpArray4(20,3) = Zpb*tmpArray3(10,3) + alphaZpq*TmpArray3(10,4) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,4) + alphaP*TmpArray2(2,5))
     TwoTerms(2) = inv2expP*(TmpArray2(3,4) + alphaP*TmpArray2(3,5))
     TwoTerms(3) = inv2expP*(TmpArray2(4,4) + alphaP*TmpArray2(4,5))
     tmpArray4(11,4) = Xpb*tmpArray3(5,4) + alphaXpq*TmpArray3(5,5) + 2*TwoTerms(1)
     tmpArray4(12,4) = Ypb*tmpArray3(5,4) + alphaYpq*TmpArray3(5,5)
     tmpArray4(13,4) = Zpb*tmpArray3(5,4) + alphaZpq*TmpArray3(5,5)
     tmpArray4(14,4) = Xpb*tmpArray3(8,4) + alphaXpq*TmpArray3(8,5)
     tmpArray4(15,4) = Xpb*tmpArray3(9,4) + alphaXpq*TmpArray3(9,5)
     tmpArray4(16,4) = Xpb*tmpArray3(10,4) + alphaXpq*TmpArray3(10,5)
     tmpArray4(17,4) = Ypb*tmpArray3(8,4) + alphaYpq*TmpArray3(8,5) + 2*TwoTerms(2)
     tmpArray4(18,4) = Zpb*tmpArray3(8,4) + alphaZpq*TmpArray3(8,5)
     tmpArray4(19,4) = Ypb*tmpArray3(10,4) + alphaYpq*TmpArray3(10,5)
     tmpArray4(20,4) = Zpb*tmpArray3(10,4) + alphaZpq*TmpArray3(10,5) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(AuxArray(IP,5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(AuxArray(IP,8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(AuxArray(IP,10) + alphaP*TmpArray3(10,2))
     AuxArray(iP,21) = Xpb*AuxArray(iP,11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(iP,22) = Ypb*AuxArray(iP,11) + alphaYpq*TmpArray4(11,2)
     AuxArray(iP,23) = Zpb*AuxArray(iP,11) + alphaZpq*TmpArray4(11,2)
     AuxArray(iP,24) = Xpb*AuxArray(iP,14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(iP,25) = Ypb*AuxArray(iP,13) + alphaYpq*TmpArray4(13,2)
     AuxArray(iP,26) = Xpb*AuxArray(iP,16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(iP,27) = Xpb*AuxArray(iP,17) + alphaXpq*TmpArray4(17,2)
     AuxArray(iP,28) = Xpb*AuxArray(iP,18) + alphaXpq*TmpArray4(18,2)
     AuxArray(iP,29) = Xpb*AuxArray(iP,19) + alphaXpq*TmpArray4(19,2)
     AuxArray(iP,30) = Xpb*AuxArray(iP,20) + alphaXpq*TmpArray4(20,2)
     AuxArray(iP,31) = Ypb*AuxArray(iP,17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(iP,32) = Zpb*AuxArray(iP,17) + alphaZpq*TmpArray4(17,2)
     AuxArray(iP,33) = Ypb*AuxArray(iP,19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(iP,34) = Ypb*AuxArray(iP,20) + alphaYpq*TmpArray4(20,2)
     AuxArray(iP,35) = Zpb*AuxArray(iP,20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,2) + alphaP*TmpArray3(5,3))
     TwoTerms(2) = inv2expP*(TmpArray3(8,2) + alphaP*TmpArray3(8,3))
     TwoTerms(3) = inv2expP*(TmpArray3(10,2) + alphaP*TmpArray3(10,3))
     tmpArray5(21,2) = Xpb*tmpArray4(11,2) + alphaXpq*TmpArray4(11,3) + 3*TwoTerms(1)
     tmpArray5(22,2) = Ypb*tmpArray4(11,2) + alphaYpq*TmpArray4(11,3)
     tmpArray5(23,2) = Zpb*tmpArray4(11,2) + alphaZpq*TmpArray4(11,3)
     tmpArray5(24,2) = Xpb*tmpArray4(14,2) + alphaXpq*TmpArray4(14,3) + TwoTerms(2)
     tmpArray5(25,2) = Ypb*tmpArray4(13,2) + alphaYpq*TmpArray4(13,3)
     tmpArray5(26,2) = Xpb*tmpArray4(16,2) + alphaXpq*TmpArray4(16,3) + TwoTerms(3)
     tmpArray5(27,2) = Xpb*tmpArray4(17,2) + alphaXpq*TmpArray4(17,3)
     tmpArray5(28,2) = Xpb*tmpArray4(18,2) + alphaXpq*TmpArray4(18,3)
     tmpArray5(29,2) = Xpb*tmpArray4(19,2) + alphaXpq*TmpArray4(19,3)
     tmpArray5(30,2) = Xpb*tmpArray4(20,2) + alphaXpq*TmpArray4(20,3)
     tmpArray5(31,2) = Ypb*tmpArray4(17,2) + alphaYpq*TmpArray4(17,3) + 3*TwoTerms(2)
     tmpArray5(32,2) = Zpb*tmpArray4(17,2) + alphaZpq*TmpArray4(17,3)
     tmpArray5(33,2) = Ypb*tmpArray4(19,2) + alphaYpq*TmpArray4(19,3) + TwoTerms(3)
     tmpArray5(34,2) = Ypb*tmpArray4(20,2) + alphaYpq*TmpArray4(20,3)
     tmpArray5(35,2) = Zpb*tmpArray4(20,2) + alphaZpq*TmpArray4(20,3) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,3) + alphaP*TmpArray3(5,4))
     TwoTerms(2) = inv2expP*(TmpArray3(8,3) + alphaP*TmpArray3(8,4))
     TwoTerms(3) = inv2expP*(TmpArray3(10,3) + alphaP*TmpArray3(10,4))
     tmpArray5(21,3) = Xpb*tmpArray4(11,3) + alphaXpq*TmpArray4(11,4) + 3*TwoTerms(1)
     tmpArray5(22,3) = Ypb*tmpArray4(11,3) + alphaYpq*TmpArray4(11,4)
     tmpArray5(23,3) = Zpb*tmpArray4(11,3) + alphaZpq*TmpArray4(11,4)
     tmpArray5(24,3) = Xpb*tmpArray4(14,3) + alphaXpq*TmpArray4(14,4) + TwoTerms(2)
     tmpArray5(25,3) = Ypb*tmpArray4(13,3) + alphaYpq*TmpArray4(13,4)
     tmpArray5(26,3) = Xpb*tmpArray4(16,3) + alphaXpq*TmpArray4(16,4) + TwoTerms(3)
     tmpArray5(27,3) = Xpb*tmpArray4(17,3) + alphaXpq*TmpArray4(17,4)
     tmpArray5(28,3) = Xpb*tmpArray4(18,3) + alphaXpq*TmpArray4(18,4)
     tmpArray5(29,3) = Xpb*tmpArray4(19,3) + alphaXpq*TmpArray4(19,4)
     tmpArray5(30,3) = Xpb*tmpArray4(20,3) + alphaXpq*TmpArray4(20,4)
     tmpArray5(31,3) = Ypb*tmpArray4(17,3) + alphaYpq*TmpArray4(17,4) + 3*TwoTerms(2)
     tmpArray5(32,3) = Zpb*tmpArray4(17,3) + alphaZpq*TmpArray4(17,4)
     tmpArray5(33,3) = Ypb*tmpArray4(19,3) + alphaYpq*TmpArray4(19,4) + TwoTerms(3)
     tmpArray5(34,3) = Ypb*tmpArray4(20,3) + alphaYpq*TmpArray4(20,4)
     tmpArray5(35,3) = Zpb*tmpArray4(20,3) + alphaZpq*TmpArray4(20,4) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(AuxArray(IP,11) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(AuxArray(IP,14) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(AuxArray(IP,16) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(AuxArray(IP,17) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(AuxArray(IP,19) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(AuxArray(IP,20) + alphaP*TmpArray4(20,2))
     AuxArray(iP,36) = Xpb*AuxArray(iP,21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(iP,37) = Ypb*AuxArray(iP,21) + alphaYpq*TmpArray5(21,2)
     AuxArray(iP,38) = Zpb*AuxArray(iP,21) + alphaZpq*TmpArray5(21,2)
     AuxArray(iP,39) = Xpb*AuxArray(iP,24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(iP,40) = Ypb*AuxArray(iP,23) + alphaYpq*TmpArray5(23,2)
     AuxArray(iP,41) = Xpb*AuxArray(iP,26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(iP,42) = Xpb*AuxArray(iP,27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(iP,43) = Zpb*AuxArray(iP,24) + alphaZpq*TmpArray5(24,2)
     AuxArray(iP,44) = Ypb*AuxArray(iP,26) + alphaYpq*TmpArray5(26,2)
     AuxArray(iP,45) = Xpb*AuxArray(iP,30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(iP,46) = Xpb*AuxArray(iP,31) + alphaXpq*TmpArray5(31,2)
     AuxArray(iP,47) = Xpb*AuxArray(iP,32) + alphaXpq*TmpArray5(32,2)
     AuxArray(iP,48) = Xpb*AuxArray(iP,33) + alphaXpq*TmpArray5(33,2)
     AuxArray(iP,49) = Xpb*AuxArray(iP,34) + alphaXpq*TmpArray5(34,2)
     AuxArray(iP,50) = Xpb*AuxArray(iP,35) + alphaXpq*TmpArray5(35,2)
     AuxArray(iP,51) = Ypb*AuxArray(iP,31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(iP,52) = Zpb*AuxArray(iP,31) + alphaZpq*TmpArray5(31,2)
     AuxArray(iP,53) = Ypb*AuxArray(iP,33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(iP,54) = Ypb*AuxArray(iP,34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(iP,55) = Ypb*AuxArray(iP,35) + alphaYpq*TmpArray5(35,2)
     AuxArray(iP,56) = Zpb*AuxArray(iP,35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TmpArray4(11,2) + alphaP*TmpArray4(11,3))
     TwoTerms(2) = inv2expP*(TmpArray4(14,2) + alphaP*TmpArray4(14,3))
     TwoTerms(3) = inv2expP*(TmpArray4(16,2) + alphaP*TmpArray4(16,3))
     TwoTerms(4) = inv2expP*(TmpArray4(17,2) + alphaP*TmpArray4(17,3))
     TwoTerms(5) = inv2expP*(TmpArray4(19,2) + alphaP*TmpArray4(19,3))
     TwoTerms(6) = inv2expP*(TmpArray4(20,2) + alphaP*TmpArray4(20,3))
     tmpArray6(36,2) = Xpb*tmpArray5(21,2) + alphaXpq*TmpArray5(21,3) + 4*TwoTerms(1)
     tmpArray6(37,2) = Ypb*tmpArray5(21,2) + alphaYpq*TmpArray5(21,3)
     tmpArray6(38,2) = Zpb*tmpArray5(21,2) + alphaZpq*TmpArray5(21,3)
     tmpArray6(39,2) = Xpb*tmpArray5(24,2) + alphaXpq*TmpArray5(24,3) + 2*TwoTerms(2)
     tmpArray6(40,2) = Ypb*tmpArray5(23,2) + alphaYpq*TmpArray5(23,3)
     tmpArray6(41,2) = Xpb*tmpArray5(26,2) + alphaXpq*TmpArray5(26,3) + 2*TwoTerms(3)
     tmpArray6(42,2) = Xpb*tmpArray5(27,2) + alphaXpq*TmpArray5(27,3) + TwoTerms(4)
     tmpArray6(43,2) = Zpb*tmpArray5(24,2) + alphaZpq*TmpArray5(24,3)
     tmpArray6(44,2) = Ypb*tmpArray5(26,2) + alphaYpq*TmpArray5(26,3)
     tmpArray6(45,2) = Xpb*tmpArray5(30,2) + alphaXpq*TmpArray5(30,3) + TwoTerms(6)
     tmpArray6(46,2) = Xpb*tmpArray5(31,2) + alphaXpq*TmpArray5(31,3)
     tmpArray6(47,2) = Xpb*tmpArray5(32,2) + alphaXpq*TmpArray5(32,3)
     tmpArray6(48,2) = Xpb*tmpArray5(33,2) + alphaXpq*TmpArray5(33,3)
     tmpArray6(49,2) = Xpb*tmpArray5(34,2) + alphaXpq*TmpArray5(34,3)
     tmpArray6(50,2) = Xpb*tmpArray5(35,2) + alphaXpq*TmpArray5(35,3)
     tmpArray6(51,2) = Ypb*tmpArray5(31,2) + alphaYpq*TmpArray5(31,3) + 4*TwoTerms(4)
     tmpArray6(52,2) = Zpb*tmpArray5(31,2) + alphaZpq*TmpArray5(31,3)
     tmpArray6(53,2) = Ypb*tmpArray5(33,2) + alphaYpq*TmpArray5(33,3) + 2*TwoTerms(5)
     tmpArray6(54,2) = Ypb*tmpArray5(34,2) + alphaYpq*TmpArray5(34,3) + TwoTerms(6)
     tmpArray6(55,2) = Ypb*tmpArray5(35,2) + alphaYpq*TmpArray5(35,3)
     tmpArray6(56,2) = Zpb*tmpArray5(35,2) + alphaZpq*TmpArray5(35,3) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(AuxArray(IP,21) + alphaP*TmpArray5(21,2))
     TwoTerms(2) = inv2expP*(AuxArray(IP,24) + alphaP*TmpArray5(24,2))
     TwoTerms(3) = inv2expP*(AuxArray(IP,26) + alphaP*TmpArray5(26,2))
     TwoTerms(4) = inv2expP*(AuxArray(IP,27) + alphaP*TmpArray5(27,2))
     TwoTerms(5) = inv2expP*(AuxArray(IP,30) + alphaP*TmpArray5(30,2))
     TwoTerms(6) = inv2expP*(AuxArray(IP,31) + alphaP*TmpArray5(31,2))
     TwoTerms(7) = inv2expP*(AuxArray(IP,33) + alphaP*TmpArray5(33,2))
     TwoTerms(8) = inv2expP*(AuxArray(IP,34) + alphaP*TmpArray5(34,2))
     TwoTerms(9) = inv2expP*(AuxArray(IP,35) + alphaP*TmpArray5(35,2))
     AuxArray(iP,57) = Xpb*AuxArray(iP,36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(iP,58) = Ypb*AuxArray(iP,36) + alphaYpq*TmpArray6(36,2)
     AuxArray(iP,59) = Zpb*AuxArray(iP,36) + alphaZpq*TmpArray6(36,2)
     AuxArray(iP,60) = Xpb*AuxArray(iP,39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(iP,61) = Ypb*AuxArray(iP,38) + alphaYpq*TmpArray6(38,2)
     AuxArray(iP,62) = Xpb*AuxArray(iP,41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(iP,63) = Xpb*AuxArray(iP,42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(iP,64) = Zpb*AuxArray(iP,39) + alphaZpq*TmpArray6(39,2)
     AuxArray(iP,65) = Ypb*AuxArray(iP,41) + alphaYpq*TmpArray6(41,2)
     AuxArray(iP,66) = Xpb*AuxArray(iP,45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(iP,67) = Xpb*AuxArray(iP,46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(iP,68) = Zpb*AuxArray(iP,42) + alphaZpq*TmpArray6(42,2)
     AuxArray(iP,69) = Xpb*AuxArray(iP,48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(iP,70) = Ypb*AuxArray(iP,45) + alphaYpq*TmpArray6(45,2)
     AuxArray(iP,71) = Xpb*AuxArray(iP,50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(iP,72) = Xpb*AuxArray(iP,51) + alphaXpq*TmpArray6(51,2)
     AuxArray(iP,73) = Xpb*AuxArray(iP,52) + alphaXpq*TmpArray6(52,2)
     AuxArray(iP,74) = Xpb*AuxArray(iP,53) + alphaXpq*TmpArray6(53,2)
     AuxArray(iP,75) = Xpb*AuxArray(iP,54) + alphaXpq*TmpArray6(54,2)
     AuxArray(iP,76) = Xpb*AuxArray(iP,55) + alphaXpq*TmpArray6(55,2)
     AuxArray(iP,77) = Xpb*AuxArray(iP,56) + alphaXpq*TmpArray6(56,2)
     AuxArray(iP,78) = Ypb*AuxArray(iP,51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(iP,79) = Zpb*AuxArray(iP,51) + alphaZpq*TmpArray6(51,2)
     AuxArray(iP,80) = Ypb*AuxArray(iP,53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(iP,81) = Ypb*AuxArray(iP,54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(iP,82) = Ypb*AuxArray(iP,55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(iP,83) = Ypb*AuxArray(iP,56) + alphaYpq*TmpArray6(56,2)
     AuxArray(iP,84) = Zpb*AuxArray(iP,56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine
end module
