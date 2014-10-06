MODULE AGC_GPU_OBS_VERTICALRECURRENCEMODDGen
 use IchorPrecisionModule
  
 CONTAINS


subroutine VerticalRecurrenceGPUGen1D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: TABFJW(0:4,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPrimP*nPassP,4)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  Integer :: iP,iPrimQ,iPrimP,iPassP,ipnt,iAtomA,iAtomB
  real(realk) :: Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ,alphaQ,RJ000(0:1)
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
  !ThetaAux(n,1,0,0) = Xqd*ThetaAux(n,0,0,0) + (-alpha/q*Xpq)*ThetaAux(n+1,0,0,0)
  !i = 0 last 2 term vanish
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         RWVAL,GVAL,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqd,Yqd,Zqd,alphaQ,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpQ,&
!$ACC         PREF,&
!$ACC         TMP1,TMP2,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,TABFJW,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Dcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPrimP*nPassP
   iPrimQ = mod(IP-1,nPrimQ)+1
   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpQ = D1/Qexp(iPrimQ)
   Xpq = mPX + Qcent(1,iPrimQ)
   Ypq = mPY + Qcent(2,iPrimQ)
   Zpq = mPZ + Qcent(3,iPrimQ)
     Xqd = Qcent(1,iPrimQ) - Dcenter(1)
     Yqd = Qcent(2,iPrimQ) - Dcenter(2)
     Zqd = Qcent(3,iPrimQ) - Dcenter(3)
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
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
     AUXarray(iP,2) = Xqd*TMP1 + alphaXpq*TMP2
     AUXarray(iP,3) = Yqd*TMP1 + alphaYpq*TMP2
     AUXarray(iP,4) = Zqd*TMP1 + alphaZpq*TMP2
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
end subroutine VerticalRecurrenceGPUGen1D

subroutine VerticalRecurrenceGPUGen2D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 2)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPrimP*nPassP,   10)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ
  real(realk) :: TwoTerms(   1)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:3)
  real(realk) :: TMParray2(  2:  4,2:2)
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqd,Yqd,Zqd,alphaQ,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpQ,inv2expQ,&
!$ACC         PREF,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Dcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPrimP*nPassP
   iPrimQ = mod(IP-1,nPrimQ)+1
   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     Xqd = Qcent(1,iPrimQ) - Dcenter(1)
     Yqd = Qcent(2,iPrimQ) - Dcenter(2)
     Zqd = Qcent(3,iPrimQ) - Dcenter(3)
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)
     Auxarray(IP,1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     AuxArray(iP,2) = Xqd*AuxArray(iP,1) + alphaXpq*TmpArray1(1,2)
     AuxArray(iP,3) = Yqd*AuxArray(iP,1) + alphaYpq*TmpArray1(1,2)
     AuxArray(iP,4) = Zqd*AuxArray(iP,1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expQ*(AuxArray(IP,1) + alphaQ*TmpArray1(1,2))
     AuxArray(iP,5) = Xqd*AuxArray(iP,2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = Xqd*AuxArray(iP,3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = Xqd*AuxArray(iP,4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = Yqd*AuxArray(iP,3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = Yqd*AuxArray(iP,4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = Zqd*AuxArray(iP,4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen3D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 3)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPrimP*nPassP,   20)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ
  real(realk) :: TwoTerms(   3)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:4)
  real(realk) :: TMParray2(  2:  4,2:3)
  real(realk) :: TMParray3(  5: 10,2:2)
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqd,Yqd,Zqd,alphaQ,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpQ,inv2expQ,&
!$ACC         PREF,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Dcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPrimP*nPassP
   iPrimQ = mod(IP-1,nPrimQ)+1
   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     Xqd = Qcent(1,iPrimQ) - Dcenter(1)
     Yqd = Qcent(2,iPrimQ) - Dcenter(2)
     Zqd = Qcent(3,iPrimQ) - Dcenter(3)
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)
     Auxarray(IP,1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 3)
     AuxArray(iP,2) = Xqd*AuxArray(iP,1) + alphaXpq*TmpArray1(1,2)
     AuxArray(iP,3) = Yqd*AuxArray(iP,1) + alphaYpq*TmpArray1(1,2)
     AuxArray(iP,4) = Zqd*AuxArray(iP,1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqd*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqd*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqd*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expQ*(AuxArray(IP,1) + alphaQ*TmpArray1(1,2))
     AuxArray(iP,5) = Xqd*AuxArray(iP,2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = Xqd*AuxArray(iP,3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = Xqd*AuxArray(iP,4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = Yqd*AuxArray(iP,3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = Yqd*AuxArray(iP,4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = Zqd*AuxArray(iP,4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,2) + alphaQ*TmpArray1(1,3))
     tmpArray3(5,2) = Xqd*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xqd*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xqd*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Yqd*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Yqd*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zqd*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(AuxArray(IP,2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(IP,3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(IP,4) + alphaQ*TmpArray2(4,2))
     AuxArray(iP,11) = Xqd*AuxArray(iP,5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(iP,12) = Yqd*AuxArray(iP,5) + alphaYpq*TmpArray3(5,2)
     AuxArray(iP,13) = Zqd*AuxArray(iP,5) + alphaZpq*TmpArray3(5,2)
     AuxArray(iP,14) = Xqd*AuxArray(iP,8) + alphaXpq*TmpArray3(8,2)
     AuxArray(iP,15) = Xqd*AuxArray(iP,9) + alphaXpq*TmpArray3(9,2)
     AuxArray(iP,16) = Xqd*AuxArray(iP,10) + alphaXpq*TmpArray3(10,2)
     AuxArray(iP,17) = Yqd*AuxArray(iP,8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(iP,18) = Zqd*AuxArray(iP,8) + alphaZpq*TmpArray3(8,2)
     AuxArray(iP,19) = Yqd*AuxArray(iP,10) + alphaYpq*TmpArray3(10,2)
     AuxArray(iP,20) = Zqd*AuxArray(iP,10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen4D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 4)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPrimP*nPassP,   35)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ
  real(realk) :: TwoTerms(   6)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:5)
  real(realk) :: TMParray2(  2:  4,2:4)
  real(realk) :: TMParray3(  5: 10,2:3)
  real(realk) :: TMParray4( 11: 20,2:2)
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqd,Yqd,Zqd,alphaQ,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpQ,inv2expQ,&
!$ACC         PREF,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TMParray4,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Dcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPrimP*nPassP
   iPrimQ = mod(IP-1,nPrimQ)+1
   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     Xqd = Qcent(1,iPrimQ) - Dcenter(1)
     Yqd = Qcent(2,iPrimQ) - Dcenter(2)
     Zqd = Qcent(3,iPrimQ) - Dcenter(3)
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)
     Auxarray(IP,1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 3)
     TMParray1(1, 5) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 4)
     AuxArray(iP,2) = Xqd*AuxArray(iP,1) + alphaXpq*TmpArray1(1,2)
     AuxArray(iP,3) = Yqd*AuxArray(iP,1) + alphaYpq*TmpArray1(1,2)
     AuxArray(iP,4) = Zqd*AuxArray(iP,1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqd*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqd*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqd*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xqd*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Yqd*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zqd*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     TwoTerms(1) = inv2expQ*(AuxArray(IP,1) + alphaQ*TmpArray1(1,2))
     AuxArray(iP,5) = Xqd*AuxArray(iP,2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = Xqd*AuxArray(iP,3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = Xqd*AuxArray(iP,4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = Yqd*AuxArray(iP,3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = Yqd*AuxArray(iP,4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = Zqd*AuxArray(iP,4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,2) + alphaQ*TmpArray1(1,3))
     tmpArray3(5,2) = Xqd*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xqd*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xqd*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Yqd*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Yqd*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zqd*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,3) + alphaQ*TmpArray1(1,4))
     tmpArray3(5,3) = Xqd*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xqd*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xqd*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Yqd*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Yqd*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zqd*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(AuxArray(IP,2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(IP,3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(IP,4) + alphaQ*TmpArray2(4,2))
     AuxArray(iP,11) = Xqd*AuxArray(iP,5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(iP,12) = Yqd*AuxArray(iP,5) + alphaYpq*TmpArray3(5,2)
     AuxArray(iP,13) = Zqd*AuxArray(iP,5) + alphaZpq*TmpArray3(5,2)
     AuxArray(iP,14) = Xqd*AuxArray(iP,8) + alphaXpq*TmpArray3(8,2)
     AuxArray(iP,15) = Xqd*AuxArray(iP,9) + alphaXpq*TmpArray3(9,2)
     AuxArray(iP,16) = Xqd*AuxArray(iP,10) + alphaXpq*TmpArray3(10,2)
     AuxArray(iP,17) = Yqd*AuxArray(iP,8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(iP,18) = Zqd*AuxArray(iP,8) + alphaZpq*TmpArray3(8,2)
     AuxArray(iP,19) = Yqd*AuxArray(iP,10) + alphaYpq*TmpArray3(10,2)
     AuxArray(iP,20) = Zqd*AuxArray(iP,10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,2) + alphaQ*TmpArray2(2,3))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,2) + alphaQ*TmpArray2(3,3))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,2) + alphaQ*TmpArray2(4,3))
     tmpArray4(11,2) = Xqd*tmpArray3(5,2) + alphaXpq*TmpArray3(5,3) + 2*TwoTerms(1)
     tmpArray4(12,2) = Yqd*tmpArray3(5,2) + alphaYpq*TmpArray3(5,3)
     tmpArray4(13,2) = Zqd*tmpArray3(5,2) + alphaZpq*TmpArray3(5,3)
     tmpArray4(14,2) = Xqd*tmpArray3(8,2) + alphaXpq*TmpArray3(8,3)
     tmpArray4(15,2) = Xqd*tmpArray3(9,2) + alphaXpq*TmpArray3(9,3)
     tmpArray4(16,2) = Xqd*tmpArray3(10,2) + alphaXpq*TmpArray3(10,3)
     tmpArray4(17,2) = Yqd*tmpArray3(8,2) + alphaYpq*TmpArray3(8,3) + 2*TwoTerms(2)
     tmpArray4(18,2) = Zqd*tmpArray3(8,2) + alphaZpq*TmpArray3(8,3)
     tmpArray4(19,2) = Yqd*tmpArray3(10,2) + alphaYpq*TmpArray3(10,3)
     tmpArray4(20,2) = Zqd*tmpArray3(10,2) + alphaZpq*TmpArray3(10,3) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(AuxArray(IP,5) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(AuxArray(IP,8) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(AuxArray(IP,10) + alphaQ*TmpArray3(10,2))
     AuxArray(iP,21) = Xqd*AuxArray(iP,11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(iP,22) = Yqd*AuxArray(iP,11) + alphaYpq*TmpArray4(11,2)
     AuxArray(iP,23) = Zqd*AuxArray(iP,11) + alphaZpq*TmpArray4(11,2)
     AuxArray(iP,24) = Xqd*AuxArray(iP,14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(iP,25) = Yqd*AuxArray(iP,13) + alphaYpq*TmpArray4(13,2)
     AuxArray(iP,26) = Xqd*AuxArray(iP,16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(iP,27) = Xqd*AuxArray(iP,17) + alphaXpq*TmpArray4(17,2)
     AuxArray(iP,28) = Xqd*AuxArray(iP,18) + alphaXpq*TmpArray4(18,2)
     AuxArray(iP,29) = Xqd*AuxArray(iP,19) + alphaXpq*TmpArray4(19,2)
     AuxArray(iP,30) = Xqd*AuxArray(iP,20) + alphaXpq*TmpArray4(20,2)
     AuxArray(iP,31) = Yqd*AuxArray(iP,17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(iP,32) = Zqd*AuxArray(iP,17) + alphaZpq*TmpArray4(17,2)
     AuxArray(iP,33) = Yqd*AuxArray(iP,19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(iP,34) = Yqd*AuxArray(iP,20) + alphaYpq*TmpArray4(20,2)
     AuxArray(iP,35) = Zqd*AuxArray(iP,20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen5D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 5)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(nPrimQ*nPrimP*nPassP,   56)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ
  real(realk) :: TwoTerms(  10)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:6)
  real(realk) :: TMParray2(  2:  4,2:5)
  real(realk) :: TMParray3(  5: 10,2:4)
  real(realk) :: TMParray4( 11: 20,2:3)
  real(realk) :: TMParray5( 21: 35,2:2)
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqd,Yqd,Zqd,alphaQ,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpQ,inv2expQ,&
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
!$ACC        Qexp,Dcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPrimP*nPassP
   iPrimQ = mod(IP-1,nPrimQ)+1
   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1
   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
     invexpQ = D1/Qexp(iPrimQ)
     inv2expQ = D05*invexpQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     Xqd = Qcent(1,iPrimQ) - Dcenter(1)
     Yqd = Qcent(2,iPrimQ) - Dcenter(2)
     Zqd = Qcent(3,iPrimQ) - Dcenter(3)
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)
     Auxarray(IP,1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 3)
     TMParray1(1, 5) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 4)
     TMParray1(1, 6) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 5)
     AuxArray(iP,2) = Xqd*AuxArray(iP,1) + alphaXpq*TmpArray1(1,2)
     AuxArray(iP,3) = Yqd*AuxArray(iP,1) + alphaYpq*TmpArray1(1,2)
     AuxArray(iP,4) = Zqd*AuxArray(iP,1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqd*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqd*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqd*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xqd*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Yqd*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zqd*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     tmpArray2(2,5) = Xqd*tmpArray1(1,5) + alphaXpq*TmpArray1(1,6)
     tmpArray2(3,5) = Yqd*tmpArray1(1,5) + alphaYpq*TmpArray1(1,6)
     tmpArray2(4,5) = Zqd*tmpArray1(1,5) + alphaZpq*TmpArray1(1,6)
     TwoTerms(1) = inv2expQ*(AuxArray(IP,1) + alphaQ*TmpArray1(1,2))
     AuxArray(iP,5) = Xqd*AuxArray(iP,2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = Xqd*AuxArray(iP,3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = Xqd*AuxArray(iP,4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = Yqd*AuxArray(iP,3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = Yqd*AuxArray(iP,4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = Zqd*AuxArray(iP,4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,2) + alphaQ*TmpArray1(1,3))
     tmpArray3(5,2) = Xqd*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xqd*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xqd*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Yqd*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Yqd*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zqd*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,3) + alphaQ*TmpArray1(1,4))
     tmpArray3(5,3) = Xqd*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xqd*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xqd*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Yqd*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Yqd*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zqd*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,4) + alphaQ*TmpArray1(1,5))
     tmpArray3(5,4) = Xqd*tmpArray2(2,4) + alphaXpq*TmpArray2(2,5) + TwoTerms(1)
     tmpArray3(6,4) = Xqd*tmpArray2(3,4) + alphaXpq*TmpArray2(3,5)
     tmpArray3(7,4) = Xqd*tmpArray2(4,4) + alphaXpq*TmpArray2(4,5)
     tmpArray3(8,4) = Yqd*tmpArray2(3,4) + alphaYpq*TmpArray2(3,5) + TwoTerms(1)
     tmpArray3(9,4) = Yqd*tmpArray2(4,4) + alphaYpq*TmpArray2(4,5)
     tmpArray3(10,4) = Zqd*tmpArray2(4,4) + alphaZpq*TmpArray2(4,5) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(AuxArray(IP,2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(IP,3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(IP,4) + alphaQ*TmpArray2(4,2))
     AuxArray(iP,11) = Xqd*AuxArray(iP,5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(iP,12) = Yqd*AuxArray(iP,5) + alphaYpq*TmpArray3(5,2)
     AuxArray(iP,13) = Zqd*AuxArray(iP,5) + alphaZpq*TmpArray3(5,2)
     AuxArray(iP,14) = Xqd*AuxArray(iP,8) + alphaXpq*TmpArray3(8,2)
     AuxArray(iP,15) = Xqd*AuxArray(iP,9) + alphaXpq*TmpArray3(9,2)
     AuxArray(iP,16) = Xqd*AuxArray(iP,10) + alphaXpq*TmpArray3(10,2)
     AuxArray(iP,17) = Yqd*AuxArray(iP,8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(iP,18) = Zqd*AuxArray(iP,8) + alphaZpq*TmpArray3(8,2)
     AuxArray(iP,19) = Yqd*AuxArray(iP,10) + alphaYpq*TmpArray3(10,2)
     AuxArray(iP,20) = Zqd*AuxArray(iP,10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,2) + alphaQ*TmpArray2(2,3))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,2) + alphaQ*TmpArray2(3,3))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,2) + alphaQ*TmpArray2(4,3))
     tmpArray4(11,2) = Xqd*tmpArray3(5,2) + alphaXpq*TmpArray3(5,3) + 2*TwoTerms(1)
     tmpArray4(12,2) = Yqd*tmpArray3(5,2) + alphaYpq*TmpArray3(5,3)
     tmpArray4(13,2) = Zqd*tmpArray3(5,2) + alphaZpq*TmpArray3(5,3)
     tmpArray4(14,2) = Xqd*tmpArray3(8,2) + alphaXpq*TmpArray3(8,3)
     tmpArray4(15,2) = Xqd*tmpArray3(9,2) + alphaXpq*TmpArray3(9,3)
     tmpArray4(16,2) = Xqd*tmpArray3(10,2) + alphaXpq*TmpArray3(10,3)
     tmpArray4(17,2) = Yqd*tmpArray3(8,2) + alphaYpq*TmpArray3(8,3) + 2*TwoTerms(2)
     tmpArray4(18,2) = Zqd*tmpArray3(8,2) + alphaZpq*TmpArray3(8,3)
     tmpArray4(19,2) = Yqd*tmpArray3(10,2) + alphaYpq*TmpArray3(10,3)
     tmpArray4(20,2) = Zqd*tmpArray3(10,2) + alphaZpq*TmpArray3(10,3) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,3) + alphaQ*TmpArray2(2,4))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,3) + alphaQ*TmpArray2(3,4))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,3) + alphaQ*TmpArray2(4,4))
     tmpArray4(11,3) = Xqd*tmpArray3(5,3) + alphaXpq*TmpArray3(5,4) + 2*TwoTerms(1)
     tmpArray4(12,3) = Yqd*tmpArray3(5,3) + alphaYpq*TmpArray3(5,4)
     tmpArray4(13,3) = Zqd*tmpArray3(5,3) + alphaZpq*TmpArray3(5,4)
     tmpArray4(14,3) = Xqd*tmpArray3(8,3) + alphaXpq*TmpArray3(8,4)
     tmpArray4(15,3) = Xqd*tmpArray3(9,3) + alphaXpq*TmpArray3(9,4)
     tmpArray4(16,3) = Xqd*tmpArray3(10,3) + alphaXpq*TmpArray3(10,4)
     tmpArray4(17,3) = Yqd*tmpArray3(8,3) + alphaYpq*TmpArray3(8,4) + 2*TwoTerms(2)
     tmpArray4(18,3) = Zqd*tmpArray3(8,3) + alphaZpq*TmpArray3(8,4)
     tmpArray4(19,3) = Yqd*tmpArray3(10,3) + alphaYpq*TmpArray3(10,4)
     tmpArray4(20,3) = Zqd*tmpArray3(10,3) + alphaZpq*TmpArray3(10,4) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(AuxArray(IP,5) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(AuxArray(IP,8) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(AuxArray(IP,10) + alphaQ*TmpArray3(10,2))
     AuxArray(iP,21) = Xqd*AuxArray(iP,11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(iP,22) = Yqd*AuxArray(iP,11) + alphaYpq*TmpArray4(11,2)
     AuxArray(iP,23) = Zqd*AuxArray(iP,11) + alphaZpq*TmpArray4(11,2)
     AuxArray(iP,24) = Xqd*AuxArray(iP,14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(iP,25) = Yqd*AuxArray(iP,13) + alphaYpq*TmpArray4(13,2)
     AuxArray(iP,26) = Xqd*AuxArray(iP,16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(iP,27) = Xqd*AuxArray(iP,17) + alphaXpq*TmpArray4(17,2)
     AuxArray(iP,28) = Xqd*AuxArray(iP,18) + alphaXpq*TmpArray4(18,2)
     AuxArray(iP,29) = Xqd*AuxArray(iP,19) + alphaXpq*TmpArray4(19,2)
     AuxArray(iP,30) = Xqd*AuxArray(iP,20) + alphaXpq*TmpArray4(20,2)
     AuxArray(iP,31) = Yqd*AuxArray(iP,17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(iP,32) = Zqd*AuxArray(iP,17) + alphaZpq*TmpArray4(17,2)
     AuxArray(iP,33) = Yqd*AuxArray(iP,19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(iP,34) = Yqd*AuxArray(iP,20) + alphaYpq*TmpArray4(20,2)
     AuxArray(iP,35) = Zqd*AuxArray(iP,20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray3(5,2) + alphaQ*TmpArray3(5,3))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,2) + alphaQ*TmpArray3(8,3))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,2) + alphaQ*TmpArray3(10,3))
     tmpArray5(21,2) = Xqd*tmpArray4(11,2) + alphaXpq*TmpArray4(11,3) + 3*TwoTerms(1)
     tmpArray5(22,2) = Yqd*tmpArray4(11,2) + alphaYpq*TmpArray4(11,3)
     tmpArray5(23,2) = Zqd*tmpArray4(11,2) + alphaZpq*TmpArray4(11,3)
     tmpArray5(24,2) = Xqd*tmpArray4(14,2) + alphaXpq*TmpArray4(14,3) + TwoTerms(2)
     tmpArray5(25,2) = Yqd*tmpArray4(13,2) + alphaYpq*TmpArray4(13,3)
     tmpArray5(26,2) = Xqd*tmpArray4(16,2) + alphaXpq*TmpArray4(16,3) + TwoTerms(3)
     tmpArray5(27,2) = Xqd*tmpArray4(17,2) + alphaXpq*TmpArray4(17,3)
     tmpArray5(28,2) = Xqd*tmpArray4(18,2) + alphaXpq*TmpArray4(18,3)
     tmpArray5(29,2) = Xqd*tmpArray4(19,2) + alphaXpq*TmpArray4(19,3)
     tmpArray5(30,2) = Xqd*tmpArray4(20,2) + alphaXpq*TmpArray4(20,3)
     tmpArray5(31,2) = Yqd*tmpArray4(17,2) + alphaYpq*TmpArray4(17,3) + 3*TwoTerms(2)
     tmpArray5(32,2) = Zqd*tmpArray4(17,2) + alphaZpq*TmpArray4(17,3)
     tmpArray5(33,2) = Yqd*tmpArray4(19,2) + alphaYpq*TmpArray4(19,3) + TwoTerms(3)
     tmpArray5(34,2) = Yqd*tmpArray4(20,2) + alphaYpq*TmpArray4(20,3)
     tmpArray5(35,2) = Zqd*tmpArray4(20,2) + alphaZpq*TmpArray4(20,3) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(AuxArray(IP,11) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(AuxArray(IP,14) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(AuxArray(IP,16) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(AuxArray(IP,17) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(AuxArray(IP,19) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(AuxArray(IP,20) + alphaQ*TmpArray4(20,2))
     AuxArray(iP,36) = Xqd*AuxArray(iP,21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(iP,37) = Yqd*AuxArray(iP,21) + alphaYpq*TmpArray5(21,2)
     AuxArray(iP,38) = Zqd*AuxArray(iP,21) + alphaZpq*TmpArray5(21,2)
     AuxArray(iP,39) = Xqd*AuxArray(iP,24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(iP,40) = Yqd*AuxArray(iP,23) + alphaYpq*TmpArray5(23,2)
     AuxArray(iP,41) = Xqd*AuxArray(iP,26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(iP,42) = Xqd*AuxArray(iP,27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(iP,43) = Zqd*AuxArray(iP,24) + alphaZpq*TmpArray5(24,2)
     AuxArray(iP,44) = Yqd*AuxArray(iP,26) + alphaYpq*TmpArray5(26,2)
     AuxArray(iP,45) = Xqd*AuxArray(iP,30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(iP,46) = Xqd*AuxArray(iP,31) + alphaXpq*TmpArray5(31,2)
     AuxArray(iP,47) = Xqd*AuxArray(iP,32) + alphaXpq*TmpArray5(32,2)
     AuxArray(iP,48) = Xqd*AuxArray(iP,33) + alphaXpq*TmpArray5(33,2)
     AuxArray(iP,49) = Xqd*AuxArray(iP,34) + alphaXpq*TmpArray5(34,2)
     AuxArray(iP,50) = Xqd*AuxArray(iP,35) + alphaXpq*TmpArray5(35,2)
     AuxArray(iP,51) = Yqd*AuxArray(iP,31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(iP,52) = Zqd*AuxArray(iP,31) + alphaZpq*TmpArray5(31,2)
     AuxArray(iP,53) = Yqd*AuxArray(iP,33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(iP,54) = Yqd*AuxArray(iP,34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(iP,55) = Yqd*AuxArray(iP,35) + alphaYpq*TmpArray5(35,2)
     AuxArray(iP,56) = Zqd*AuxArray(iP,35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine
end module
