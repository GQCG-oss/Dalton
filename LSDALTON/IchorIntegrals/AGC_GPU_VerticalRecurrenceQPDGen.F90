MODULE AGC_GPU_OBS_VERTICALRECURRENCEMODDGen
 use IchorPrecisionModule
  
 CONTAINS


subroutine VerticalRecurrenceGPUGen1D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: TABFJW(0:4,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(4,nPrimQ*nPrimP*nPassP)
  !local variables
  Integer :: iP,iPrimQ,iPrimP,iPassP,ipnt,iAtomA,iAtomB
  real(realk) :: Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ,alphaQ,RJ000(0:1)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D2=2.0E0_realk
  REAL(REALK),PARAMETER :: D2JP36=  3.8000000000000000E+01_realk
  real(realk),parameter :: D1=1.0E0_realk,D03333=1.0E0_realk/3.0E0_realk
  REAL(REALK),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk
  REAL(REALK),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
  REAL(REALK),PARAMETER :: SMALL = 1E-15_realk,D12 = 12.0E0_realk
  REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
  REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
  REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
  REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
  Real(realk), parameter :: PI=3.14159265358979323846E0_realk
  REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
  REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
  REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
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
!$ACC        nPrimP,nPrimQ,nPassP)
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
     AUXarray(1,iP) = TMP1
     AUXarray(2,iP) = Xqd*TMP1 + alphaXpq*TMP2
     AUXarray(3,iP) = Yqd*TMP1 + alphaYpq*TMP2
     AUXarray(4,iP) = Zqd*TMP1 + alphaZpq*TMP2
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
end subroutine VerticalRecurrenceGPUGen1D

subroutine VerticalRecurrenceGPUGen2D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 2,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   10,nPrimQ*nPrimP*nPassP)
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
!$ACC        nPrimP,nPrimQ,nPassP)
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
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xqd*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Yqd*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zqd*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,iP) = Xqd*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xqd*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xqd*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Yqd*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Yqd*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zqd*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen3D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 3,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   20,nPrimQ*nPrimP*nPassP)
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
!$ACC        nPrimP,nPrimQ,nPassP)
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
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xqd*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Yqd*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zqd*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqd*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqd*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqd*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,iP) = Xqd*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xqd*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xqd*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Yqd*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Yqd*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zqd*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,2) + alphaQ*TmpArray1(1,3))
     tmpArray3(5,2) = Xqd*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xqd*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xqd*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Yqd*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Yqd*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zqd*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(AuxArray(2,IP) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(3,IP) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(4,IP) + alphaQ*TmpArray2(4,2))
     AuxArray(11,iP) = Xqd*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Yqd*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zqd*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xqd*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xqd*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xqd*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Yqd*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zqd*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Yqd*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zqd*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen4D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 4,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   35,nPrimQ*nPrimP*nPassP)
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
!$ACC        nPrimP,nPrimQ,nPassP)
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
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xqd*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Yqd*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zqd*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqd*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqd*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqd*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xqd*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Yqd*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zqd*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,iP) = Xqd*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xqd*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xqd*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Yqd*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Yqd*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zqd*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expQ*(AuxArray(2,IP) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(3,IP) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(4,IP) + alphaQ*TmpArray2(4,2))
     AuxArray(11,iP) = Xqd*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Yqd*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zqd*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xqd*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xqd*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xqd*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Yqd*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zqd*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Yqd*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zqd*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(AuxArray(5,IP) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(AuxArray(8,IP) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(AuxArray(10,IP) + alphaQ*TmpArray3(10,2))
     AuxArray(21,iP) = Xqd*AuxArray(11,iP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,iP) = Yqd*AuxArray(11,iP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,iP) = Zqd*AuxArray(11,iP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,iP) = Xqd*AuxArray(14,iP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,iP) = Yqd*AuxArray(13,iP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,iP) = Xqd*AuxArray(16,iP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,iP) = Xqd*AuxArray(17,iP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,iP) = Xqd*AuxArray(18,iP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,iP) = Xqd*AuxArray(19,iP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,iP) = Xqd*AuxArray(20,iP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,iP) = Yqd*AuxArray(17,iP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,iP) = Zqd*AuxArray(17,iP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,iP) = Yqd*AuxArray(19,iP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,iP) = Yqd*AuxArray(20,iP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,iP) = Zqd*AuxArray(20,iP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen5D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 5,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   56,nPrimQ*nPrimP*nPassP)
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
!$ACC        nPrimP,nPrimQ,nPassP)
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
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 6) = PREF*RJ000Array( 5,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xqd*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Yqd*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zqd*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,iP) = Xqd*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xqd*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xqd*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Yqd*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Yqd*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zqd*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expQ*(AuxArray(2,IP) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(3,IP) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(4,IP) + alphaQ*TmpArray2(4,2))
     AuxArray(11,iP) = Xqd*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Yqd*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zqd*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xqd*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xqd*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xqd*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Yqd*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zqd*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Yqd*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zqd*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(AuxArray(5,IP) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(AuxArray(8,IP) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(AuxArray(10,IP) + alphaQ*TmpArray3(10,2))
     AuxArray(21,iP) = Xqd*AuxArray(11,iP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,iP) = Yqd*AuxArray(11,iP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,iP) = Zqd*AuxArray(11,iP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,iP) = Xqd*AuxArray(14,iP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,iP) = Yqd*AuxArray(13,iP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,iP) = Xqd*AuxArray(16,iP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,iP) = Xqd*AuxArray(17,iP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,iP) = Xqd*AuxArray(18,iP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,iP) = Xqd*AuxArray(19,iP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,iP) = Xqd*AuxArray(20,iP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,iP) = Yqd*AuxArray(17,iP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,iP) = Zqd*AuxArray(17,iP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,iP) = Yqd*AuxArray(19,iP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,iP) = Yqd*AuxArray(20,iP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,iP) = Zqd*AuxArray(20,iP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(AuxArray(11,IP) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(AuxArray(14,IP) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(AuxArray(16,IP) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(AuxArray(17,IP) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(AuxArray(19,IP) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(AuxArray(20,IP) + alphaQ*TmpArray4(20,2))
     AuxArray(36,iP) = Xqd*AuxArray(21,iP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,iP) = Yqd*AuxArray(21,iP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,iP) = Zqd*AuxArray(21,iP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,iP) = Xqd*AuxArray(24,iP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,iP) = Yqd*AuxArray(23,iP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,iP) = Xqd*AuxArray(26,iP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,iP) = Xqd*AuxArray(27,iP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,iP) = Zqd*AuxArray(24,iP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,iP) = Yqd*AuxArray(26,iP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,iP) = Xqd*AuxArray(30,iP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,iP) = Xqd*AuxArray(31,iP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,iP) = Xqd*AuxArray(32,iP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,iP) = Xqd*AuxArray(33,iP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,iP) = Xqd*AuxArray(34,iP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,iP) = Xqd*AuxArray(35,iP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,iP) = Yqd*AuxArray(31,iP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,iP) = Zqd*AuxArray(31,iP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,iP) = Yqd*AuxArray(33,iP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,iP) = Yqd*AuxArray(34,iP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,iP) = Yqd*AuxArray(35,iP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,iP) = Zqd*AuxArray(35,iP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine
end module
