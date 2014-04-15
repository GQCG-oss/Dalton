MODULE AGC_GPU_OBS_VERTICALRECURRENCEMODCGen
 use IchorPrecisionModule
  
 CONTAINS


subroutine VerticalRecurrenceGPUGen1C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
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
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(4,nPrimQ*nPrimP*nPassP)
  !local variables
  Integer :: iP,iPrimQ,iPrimP,iPassP,ipnt,iAtomA,iAtomB
  real(realk) :: Xqc,Yqc,Zqc
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
  !ThetaAux(n,1,0,0) = Xqc*ThetaAux(n,0,0,0) + (-alpha/q*Xpq)*ThetaAux(n+1,0,0,0)
  !i = 0 last 2 term vanish
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         RWVAL,GVAL,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqc,Yqc,Zqc,alphaQ,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpQ,&
!$ACC         PREF,&
!$ACC         TMP1,TMP2,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,TABFJW,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Ccenter, &
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
     Xqc = Qcent(1,iPrimQ) - Ccenter(1)
     Yqc = Qcent(2,iPrimQ) - Ccenter(2)
     Zqc = Qcent(3,iPrimQ) - Ccenter(3)
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
     AUXarray(2,iP) = Xqc*TMP1 + alphaXpq*TMP2
     AUXarray(3,iP) = Yqc*TMP1 + alphaYpq*TMP2
     AUXarray(4,iP) = Zqc*TMP1 + alphaZpq*TMP2
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
end subroutine VerticalRecurrenceGPUGen1C

subroutine VerticalRecurrenceGPUGen2C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 2,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(   10,nPrimQ*nPrimP*nPassP)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ
  real(realk) :: TwoTerms(   1)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:3)
  real(realk) :: TMParray2(  2:  4,2:2)
  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqc,Yqc,Zqc,alphaQ,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpQ,inv2expQ,&
!$ACC         PREF,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Ccenter, &
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
     Xqc = Qcent(1,iPrimQ) - Ccenter(1)
     Yqc = Qcent(2,iPrimQ) - Ccenter(2)
     Zqc = Qcent(3,iPrimQ) - Ccenter(3)
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xqc*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Yqc*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zqc*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqc*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqc*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqc*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,iP) = Xqc*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xqc*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xqc*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Yqc*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Yqc*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zqc*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen3C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 3,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(   20,nPrimQ*nPrimP*nPassP)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ
  real(realk) :: TwoTerms(   3)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:4)
  real(realk) :: TMParray2(  2:  4,2:3)
  real(realk) :: TMParray3(  5: 10,2:2)
  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqc,Yqc,Zqc,alphaQ,&
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
!$ACC        Qexp,Ccenter, &
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
     Xqc = Qcent(1,iPrimQ) - Ccenter(1)
     Yqc = Qcent(2,iPrimQ) - Ccenter(2)
     Zqc = Qcent(3,iPrimQ) - Ccenter(3)
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xqc*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Yqc*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zqc*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqc*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqc*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqc*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqc*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqc*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqc*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,iP) = Xqc*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xqc*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xqc*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Yqc*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Yqc*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zqc*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,2) + alphaQ*TmpArray1(1,3))
     tmpArray3(5,2) = Xqc*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xqc*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xqc*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Yqc*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Yqc*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zqc*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(AuxArray(2,IP) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(3,IP) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(4,IP) + alphaQ*TmpArray2(4,2))
     AuxArray(11,iP) = Xqc*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Yqc*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zqc*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xqc*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xqc*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xqc*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Yqc*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zqc*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Yqc*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zqc*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen4C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 4,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(   35,nPrimQ*nPrimP*nPassP)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ
  real(realk) :: TwoTerms(   6)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:5)
  real(realk) :: TMParray2(  2:  4,2:4)
  real(realk) :: TMParray3(  5: 10,2:3)
  real(realk) :: TMParray4( 11: 20,2:2)
  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqc,Yqc,Zqc,alphaQ,&
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
!$ACC        Qexp,Ccenter, &
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
     Xqc = Qcent(1,iPrimQ) - Ccenter(1)
     Yqc = Qcent(2,iPrimQ) - Ccenter(2)
     Zqc = Qcent(3,iPrimQ) - Ccenter(3)
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
     AuxArray(2,iP) = Xqc*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Yqc*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zqc*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqc*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqc*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqc*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqc*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqc*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqc*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xqc*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Yqc*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zqc*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,iP) = Xqc*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xqc*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xqc*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Yqc*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Yqc*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zqc*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,2) + alphaQ*TmpArray1(1,3))
     tmpArray3(5,2) = Xqc*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xqc*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xqc*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Yqc*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Yqc*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zqc*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,3) + alphaQ*TmpArray1(1,4))
     tmpArray3(5,3) = Xqc*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xqc*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xqc*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Yqc*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Yqc*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zqc*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(AuxArray(2,IP) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(3,IP) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(4,IP) + alphaQ*TmpArray2(4,2))
     AuxArray(11,iP) = Xqc*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Yqc*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zqc*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xqc*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xqc*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xqc*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Yqc*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zqc*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Yqc*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zqc*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,2) + alphaQ*TmpArray2(2,3))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,2) + alphaQ*TmpArray2(3,3))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,2) + alphaQ*TmpArray2(4,3))
     tmpArray4(11,2) = Xqc*tmpArray3(5,2) + alphaXpq*TmpArray3(5,3) + 2*TwoTerms(1)
     tmpArray4(12,2) = Yqc*tmpArray3(5,2) + alphaYpq*TmpArray3(5,3)
     tmpArray4(13,2) = Zqc*tmpArray3(5,2) + alphaZpq*TmpArray3(5,3)
     tmpArray4(14,2) = Xqc*tmpArray3(8,2) + alphaXpq*TmpArray3(8,3)
     tmpArray4(15,2) = Xqc*tmpArray3(9,2) + alphaXpq*TmpArray3(9,3)
     tmpArray4(16,2) = Xqc*tmpArray3(10,2) + alphaXpq*TmpArray3(10,3)
     tmpArray4(17,2) = Yqc*tmpArray3(8,2) + alphaYpq*TmpArray3(8,3) + 2*TwoTerms(2)
     tmpArray4(18,2) = Zqc*tmpArray3(8,2) + alphaZpq*TmpArray3(8,3)
     tmpArray4(19,2) = Yqc*tmpArray3(10,2) + alphaYpq*TmpArray3(10,3)
     tmpArray4(20,2) = Zqc*tmpArray3(10,2) + alphaZpq*TmpArray3(10,3) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(AuxArray(5,IP) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(AuxArray(8,IP) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(AuxArray(10,IP) + alphaQ*TmpArray3(10,2))
     AuxArray(21,iP) = Xqc*AuxArray(11,iP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,iP) = Yqc*AuxArray(11,iP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,iP) = Zqc*AuxArray(11,iP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,iP) = Xqc*AuxArray(14,iP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,iP) = Yqc*AuxArray(13,iP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,iP) = Xqc*AuxArray(16,iP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,iP) = Xqc*AuxArray(17,iP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,iP) = Xqc*AuxArray(18,iP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,iP) = Xqc*AuxArray(19,iP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,iP) = Xqc*AuxArray(20,iP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,iP) = Yqc*AuxArray(17,iP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,iP) = Zqc*AuxArray(17,iP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,iP) = Yqc*AuxArray(19,iP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,iP) = Yqc*AuxArray(20,iP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,iP) = Zqc*AuxArray(20,iP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen5C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 5,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(   56,nPrimQ*nPrimP*nPassP)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Xqc,Yqc,Zqc
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
  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqc,Yqc,Zqc,alphaQ,&
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
!$ACC        Qexp,Ccenter, &
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
     Xqc = Qcent(1,iPrimQ) - Ccenter(1)
     Yqc = Qcent(2,iPrimQ) - Ccenter(2)
     Zqc = Qcent(3,iPrimQ) - Ccenter(3)
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
     AuxArray(2,iP) = Xqc*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Yqc*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zqc*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqc*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqc*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqc*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqc*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqc*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqc*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xqc*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Yqc*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zqc*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     tmpArray2(2,5) = Xqc*tmpArray1(1,5) + alphaXpq*TmpArray1(1,6)
     tmpArray2(3,5) = Yqc*tmpArray1(1,5) + alphaYpq*TmpArray1(1,6)
     tmpArray2(4,5) = Zqc*tmpArray1(1,5) + alphaZpq*TmpArray1(1,6)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,iP) = Xqc*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xqc*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xqc*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Yqc*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Yqc*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zqc*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,2) + alphaQ*TmpArray1(1,3))
     tmpArray3(5,2) = Xqc*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xqc*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xqc*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Yqc*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Yqc*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zqc*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,3) + alphaQ*TmpArray1(1,4))
     tmpArray3(5,3) = Xqc*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xqc*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xqc*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Yqc*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Yqc*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zqc*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,4) + alphaQ*TmpArray1(1,5))
     tmpArray3(5,4) = Xqc*tmpArray2(2,4) + alphaXpq*TmpArray2(2,5) + TwoTerms(1)
     tmpArray3(6,4) = Xqc*tmpArray2(3,4) + alphaXpq*TmpArray2(3,5)
     tmpArray3(7,4) = Xqc*tmpArray2(4,4) + alphaXpq*TmpArray2(4,5)
     tmpArray3(8,4) = Yqc*tmpArray2(3,4) + alphaYpq*TmpArray2(3,5) + TwoTerms(1)
     tmpArray3(9,4) = Yqc*tmpArray2(4,4) + alphaYpq*TmpArray2(4,5)
     tmpArray3(10,4) = Zqc*tmpArray2(4,4) + alphaZpq*TmpArray2(4,5) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(AuxArray(2,IP) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(3,IP) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(4,IP) + alphaQ*TmpArray2(4,2))
     AuxArray(11,iP) = Xqc*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Yqc*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zqc*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xqc*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xqc*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xqc*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Yqc*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zqc*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Yqc*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zqc*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,2) + alphaQ*TmpArray2(2,3))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,2) + alphaQ*TmpArray2(3,3))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,2) + alphaQ*TmpArray2(4,3))
     tmpArray4(11,2) = Xqc*tmpArray3(5,2) + alphaXpq*TmpArray3(5,3) + 2*TwoTerms(1)
     tmpArray4(12,2) = Yqc*tmpArray3(5,2) + alphaYpq*TmpArray3(5,3)
     tmpArray4(13,2) = Zqc*tmpArray3(5,2) + alphaZpq*TmpArray3(5,3)
     tmpArray4(14,2) = Xqc*tmpArray3(8,2) + alphaXpq*TmpArray3(8,3)
     tmpArray4(15,2) = Xqc*tmpArray3(9,2) + alphaXpq*TmpArray3(9,3)
     tmpArray4(16,2) = Xqc*tmpArray3(10,2) + alphaXpq*TmpArray3(10,3)
     tmpArray4(17,2) = Yqc*tmpArray3(8,2) + alphaYpq*TmpArray3(8,3) + 2*TwoTerms(2)
     tmpArray4(18,2) = Zqc*tmpArray3(8,2) + alphaZpq*TmpArray3(8,3)
     tmpArray4(19,2) = Yqc*tmpArray3(10,2) + alphaYpq*TmpArray3(10,3)
     tmpArray4(20,2) = Zqc*tmpArray3(10,2) + alphaZpq*TmpArray3(10,3) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,3) + alphaQ*TmpArray2(2,4))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,3) + alphaQ*TmpArray2(3,4))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,3) + alphaQ*TmpArray2(4,4))
     tmpArray4(11,3) = Xqc*tmpArray3(5,3) + alphaXpq*TmpArray3(5,4) + 2*TwoTerms(1)
     tmpArray4(12,3) = Yqc*tmpArray3(5,3) + alphaYpq*TmpArray3(5,4)
     tmpArray4(13,3) = Zqc*tmpArray3(5,3) + alphaZpq*TmpArray3(5,4)
     tmpArray4(14,3) = Xqc*tmpArray3(8,3) + alphaXpq*TmpArray3(8,4)
     tmpArray4(15,3) = Xqc*tmpArray3(9,3) + alphaXpq*TmpArray3(9,4)
     tmpArray4(16,3) = Xqc*tmpArray3(10,3) + alphaXpq*TmpArray3(10,4)
     tmpArray4(17,3) = Yqc*tmpArray3(8,3) + alphaYpq*TmpArray3(8,4) + 2*TwoTerms(2)
     tmpArray4(18,3) = Zqc*tmpArray3(8,3) + alphaZpq*TmpArray3(8,4)
     tmpArray4(19,3) = Yqc*tmpArray3(10,3) + alphaYpq*TmpArray3(10,4)
     tmpArray4(20,3) = Zqc*tmpArray3(10,3) + alphaZpq*TmpArray3(10,4) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(AuxArray(5,IP) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(AuxArray(8,IP) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(AuxArray(10,IP) + alphaQ*TmpArray3(10,2))
     AuxArray(21,iP) = Xqc*AuxArray(11,iP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,iP) = Yqc*AuxArray(11,iP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,iP) = Zqc*AuxArray(11,iP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,iP) = Xqc*AuxArray(14,iP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,iP) = Yqc*AuxArray(13,iP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,iP) = Xqc*AuxArray(16,iP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,iP) = Xqc*AuxArray(17,iP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,iP) = Xqc*AuxArray(18,iP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,iP) = Xqc*AuxArray(19,iP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,iP) = Xqc*AuxArray(20,iP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,iP) = Yqc*AuxArray(17,iP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,iP) = Zqc*AuxArray(17,iP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,iP) = Yqc*AuxArray(19,iP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,iP) = Yqc*AuxArray(20,iP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,iP) = Zqc*AuxArray(20,iP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray3(5,2) + alphaQ*TmpArray3(5,3))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,2) + alphaQ*TmpArray3(8,3))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,2) + alphaQ*TmpArray3(10,3))
     tmpArray5(21,2) = Xqc*tmpArray4(11,2) + alphaXpq*TmpArray4(11,3) + 3*TwoTerms(1)
     tmpArray5(22,2) = Yqc*tmpArray4(11,2) + alphaYpq*TmpArray4(11,3)
     tmpArray5(23,2) = Zqc*tmpArray4(11,2) + alphaZpq*TmpArray4(11,3)
     tmpArray5(24,2) = Xqc*tmpArray4(14,2) + alphaXpq*TmpArray4(14,3) + TwoTerms(2)
     tmpArray5(25,2) = Yqc*tmpArray4(13,2) + alphaYpq*TmpArray4(13,3)
     tmpArray5(26,2) = Xqc*tmpArray4(16,2) + alphaXpq*TmpArray4(16,3) + TwoTerms(3)
     tmpArray5(27,2) = Xqc*tmpArray4(17,2) + alphaXpq*TmpArray4(17,3)
     tmpArray5(28,2) = Xqc*tmpArray4(18,2) + alphaXpq*TmpArray4(18,3)
     tmpArray5(29,2) = Xqc*tmpArray4(19,2) + alphaXpq*TmpArray4(19,3)
     tmpArray5(30,2) = Xqc*tmpArray4(20,2) + alphaXpq*TmpArray4(20,3)
     tmpArray5(31,2) = Yqc*tmpArray4(17,2) + alphaYpq*TmpArray4(17,3) + 3*TwoTerms(2)
     tmpArray5(32,2) = Zqc*tmpArray4(17,2) + alphaZpq*TmpArray4(17,3)
     tmpArray5(33,2) = Yqc*tmpArray4(19,2) + alphaYpq*TmpArray4(19,3) + TwoTerms(3)
     tmpArray5(34,2) = Yqc*tmpArray4(20,2) + alphaYpq*TmpArray4(20,3)
     tmpArray5(35,2) = Zqc*tmpArray4(20,2) + alphaZpq*TmpArray4(20,3) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(AuxArray(11,IP) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(AuxArray(14,IP) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(AuxArray(16,IP) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(AuxArray(17,IP) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(AuxArray(19,IP) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(AuxArray(20,IP) + alphaQ*TmpArray4(20,2))
     AuxArray(36,iP) = Xqc*AuxArray(21,iP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,iP) = Yqc*AuxArray(21,iP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,iP) = Zqc*AuxArray(21,iP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,iP) = Xqc*AuxArray(24,iP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,iP) = Yqc*AuxArray(23,iP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,iP) = Xqc*AuxArray(26,iP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,iP) = Xqc*AuxArray(27,iP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,iP) = Zqc*AuxArray(24,iP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,iP) = Yqc*AuxArray(26,iP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,iP) = Xqc*AuxArray(30,iP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,iP) = Xqc*AuxArray(31,iP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,iP) = Xqc*AuxArray(32,iP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,iP) = Xqc*AuxArray(33,iP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,iP) = Xqc*AuxArray(34,iP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,iP) = Xqc*AuxArray(35,iP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,iP) = Yqc*AuxArray(31,iP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,iP) = Zqc*AuxArray(31,iP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,iP) = Yqc*AuxArray(33,iP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,iP) = Yqc*AuxArray(34,iP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,iP) = Yqc*AuxArray(35,iP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,iP) = Zqc*AuxArray(35,iP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen6C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 6,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(   84,nPrimQ*nPrimP*nPassP)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ
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
  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqc,Yqc,Zqc,alphaQ,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpQ,inv2expQ,&
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
!$ACC        Qexp,Ccenter, &
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
     Xqc = Qcent(1,iPrimQ) - Ccenter(1)
     Yqc = Qcent(2,iPrimQ) - Ccenter(2)
     Zqc = Qcent(3,iPrimQ) - Ccenter(3)
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
     TMParray1(1, 7) = PREF*RJ000Array( 6,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xqc*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Yqc*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zqc*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqc*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqc*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqc*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqc*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqc*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqc*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xqc*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Yqc*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zqc*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     tmpArray2(2,5) = Xqc*tmpArray1(1,5) + alphaXpq*TmpArray1(1,6)
     tmpArray2(3,5) = Yqc*tmpArray1(1,5) + alphaYpq*TmpArray1(1,6)
     tmpArray2(4,5) = Zqc*tmpArray1(1,5) + alphaZpq*TmpArray1(1,6)
     tmpArray2(2,6) = Xqc*tmpArray1(1,6) + alphaXpq*TmpArray1(1,7)
     tmpArray2(3,6) = Yqc*tmpArray1(1,6) + alphaYpq*TmpArray1(1,7)
     tmpArray2(4,6) = Zqc*tmpArray1(1,6) + alphaZpq*TmpArray1(1,7)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,iP) = Xqc*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xqc*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xqc*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Yqc*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Yqc*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zqc*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,2) + alphaQ*TmpArray1(1,3))
     tmpArray3(5,2) = Xqc*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xqc*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xqc*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Yqc*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Yqc*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zqc*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,3) + alphaQ*TmpArray1(1,4))
     tmpArray3(5,3) = Xqc*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xqc*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xqc*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Yqc*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Yqc*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zqc*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,4) + alphaQ*TmpArray1(1,5))
     tmpArray3(5,4) = Xqc*tmpArray2(2,4) + alphaXpq*TmpArray2(2,5) + TwoTerms(1)
     tmpArray3(6,4) = Xqc*tmpArray2(3,4) + alphaXpq*TmpArray2(3,5)
     tmpArray3(7,4) = Xqc*tmpArray2(4,4) + alphaXpq*TmpArray2(4,5)
     tmpArray3(8,4) = Yqc*tmpArray2(3,4) + alphaYpq*TmpArray2(3,5) + TwoTerms(1)
     tmpArray3(9,4) = Yqc*tmpArray2(4,4) + alphaYpq*TmpArray2(4,5)
     tmpArray3(10,4) = Zqc*tmpArray2(4,4) + alphaZpq*TmpArray2(4,5) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,5) + alphaQ*TmpArray1(1,6))
     tmpArray3(5,5) = Xqc*tmpArray2(2,5) + alphaXpq*TmpArray2(2,6) + TwoTerms(1)
     tmpArray3(6,5) = Xqc*tmpArray2(3,5) + alphaXpq*TmpArray2(3,6)
     tmpArray3(7,5) = Xqc*tmpArray2(4,5) + alphaXpq*TmpArray2(4,6)
     tmpArray3(8,5) = Yqc*tmpArray2(3,5) + alphaYpq*TmpArray2(3,6) + TwoTerms(1)
     tmpArray3(9,5) = Yqc*tmpArray2(4,5) + alphaYpq*TmpArray2(4,6)
     tmpArray3(10,5) = Zqc*tmpArray2(4,5) + alphaZpq*TmpArray2(4,6) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(AuxArray(2,IP) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(3,IP) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(4,IP) + alphaQ*TmpArray2(4,2))
     AuxArray(11,iP) = Xqc*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Yqc*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zqc*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xqc*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xqc*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xqc*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Yqc*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zqc*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Yqc*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zqc*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,2) + alphaQ*TmpArray2(2,3))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,2) + alphaQ*TmpArray2(3,3))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,2) + alphaQ*TmpArray2(4,3))
     tmpArray4(11,2) = Xqc*tmpArray3(5,2) + alphaXpq*TmpArray3(5,3) + 2*TwoTerms(1)
     tmpArray4(12,2) = Yqc*tmpArray3(5,2) + alphaYpq*TmpArray3(5,3)
     tmpArray4(13,2) = Zqc*tmpArray3(5,2) + alphaZpq*TmpArray3(5,3)
     tmpArray4(14,2) = Xqc*tmpArray3(8,2) + alphaXpq*TmpArray3(8,3)
     tmpArray4(15,2) = Xqc*tmpArray3(9,2) + alphaXpq*TmpArray3(9,3)
     tmpArray4(16,2) = Xqc*tmpArray3(10,2) + alphaXpq*TmpArray3(10,3)
     tmpArray4(17,2) = Yqc*tmpArray3(8,2) + alphaYpq*TmpArray3(8,3) + 2*TwoTerms(2)
     tmpArray4(18,2) = Zqc*tmpArray3(8,2) + alphaZpq*TmpArray3(8,3)
     tmpArray4(19,2) = Yqc*tmpArray3(10,2) + alphaYpq*TmpArray3(10,3)
     tmpArray4(20,2) = Zqc*tmpArray3(10,2) + alphaZpq*TmpArray3(10,3) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,3) + alphaQ*TmpArray2(2,4))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,3) + alphaQ*TmpArray2(3,4))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,3) + alphaQ*TmpArray2(4,4))
     tmpArray4(11,3) = Xqc*tmpArray3(5,3) + alphaXpq*TmpArray3(5,4) + 2*TwoTerms(1)
     tmpArray4(12,3) = Yqc*tmpArray3(5,3) + alphaYpq*TmpArray3(5,4)
     tmpArray4(13,3) = Zqc*tmpArray3(5,3) + alphaZpq*TmpArray3(5,4)
     tmpArray4(14,3) = Xqc*tmpArray3(8,3) + alphaXpq*TmpArray3(8,4)
     tmpArray4(15,3) = Xqc*tmpArray3(9,3) + alphaXpq*TmpArray3(9,4)
     tmpArray4(16,3) = Xqc*tmpArray3(10,3) + alphaXpq*TmpArray3(10,4)
     tmpArray4(17,3) = Yqc*tmpArray3(8,3) + alphaYpq*TmpArray3(8,4) + 2*TwoTerms(2)
     tmpArray4(18,3) = Zqc*tmpArray3(8,3) + alphaZpq*TmpArray3(8,4)
     tmpArray4(19,3) = Yqc*tmpArray3(10,3) + alphaYpq*TmpArray3(10,4)
     tmpArray4(20,3) = Zqc*tmpArray3(10,3) + alphaZpq*TmpArray3(10,4) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,4) + alphaQ*TmpArray2(2,5))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,4) + alphaQ*TmpArray2(3,5))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,4) + alphaQ*TmpArray2(4,5))
     tmpArray4(11,4) = Xqc*tmpArray3(5,4) + alphaXpq*TmpArray3(5,5) + 2*TwoTerms(1)
     tmpArray4(12,4) = Yqc*tmpArray3(5,4) + alphaYpq*TmpArray3(5,5)
     tmpArray4(13,4) = Zqc*tmpArray3(5,4) + alphaZpq*TmpArray3(5,5)
     tmpArray4(14,4) = Xqc*tmpArray3(8,4) + alphaXpq*TmpArray3(8,5)
     tmpArray4(15,4) = Xqc*tmpArray3(9,4) + alphaXpq*TmpArray3(9,5)
     tmpArray4(16,4) = Xqc*tmpArray3(10,4) + alphaXpq*TmpArray3(10,5)
     tmpArray4(17,4) = Yqc*tmpArray3(8,4) + alphaYpq*TmpArray3(8,5) + 2*TwoTerms(2)
     tmpArray4(18,4) = Zqc*tmpArray3(8,4) + alphaZpq*TmpArray3(8,5)
     tmpArray4(19,4) = Yqc*tmpArray3(10,4) + alphaYpq*TmpArray3(10,5)
     tmpArray4(20,4) = Zqc*tmpArray3(10,4) + alphaZpq*TmpArray3(10,5) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(AuxArray(5,IP) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(AuxArray(8,IP) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(AuxArray(10,IP) + alphaQ*TmpArray3(10,2))
     AuxArray(21,iP) = Xqc*AuxArray(11,iP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,iP) = Yqc*AuxArray(11,iP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,iP) = Zqc*AuxArray(11,iP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,iP) = Xqc*AuxArray(14,iP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,iP) = Yqc*AuxArray(13,iP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,iP) = Xqc*AuxArray(16,iP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,iP) = Xqc*AuxArray(17,iP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,iP) = Xqc*AuxArray(18,iP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,iP) = Xqc*AuxArray(19,iP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,iP) = Xqc*AuxArray(20,iP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,iP) = Yqc*AuxArray(17,iP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,iP) = Zqc*AuxArray(17,iP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,iP) = Yqc*AuxArray(19,iP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,iP) = Yqc*AuxArray(20,iP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,iP) = Zqc*AuxArray(20,iP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray3(5,2) + alphaQ*TmpArray3(5,3))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,2) + alphaQ*TmpArray3(8,3))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,2) + alphaQ*TmpArray3(10,3))
     tmpArray5(21,2) = Xqc*tmpArray4(11,2) + alphaXpq*TmpArray4(11,3) + 3*TwoTerms(1)
     tmpArray5(22,2) = Yqc*tmpArray4(11,2) + alphaYpq*TmpArray4(11,3)
     tmpArray5(23,2) = Zqc*tmpArray4(11,2) + alphaZpq*TmpArray4(11,3)
     tmpArray5(24,2) = Xqc*tmpArray4(14,2) + alphaXpq*TmpArray4(14,3) + TwoTerms(2)
     tmpArray5(25,2) = Yqc*tmpArray4(13,2) + alphaYpq*TmpArray4(13,3)
     tmpArray5(26,2) = Xqc*tmpArray4(16,2) + alphaXpq*TmpArray4(16,3) + TwoTerms(3)
     tmpArray5(27,2) = Xqc*tmpArray4(17,2) + alphaXpq*TmpArray4(17,3)
     tmpArray5(28,2) = Xqc*tmpArray4(18,2) + alphaXpq*TmpArray4(18,3)
     tmpArray5(29,2) = Xqc*tmpArray4(19,2) + alphaXpq*TmpArray4(19,3)
     tmpArray5(30,2) = Xqc*tmpArray4(20,2) + alphaXpq*TmpArray4(20,3)
     tmpArray5(31,2) = Yqc*tmpArray4(17,2) + alphaYpq*TmpArray4(17,3) + 3*TwoTerms(2)
     tmpArray5(32,2) = Zqc*tmpArray4(17,2) + alphaZpq*TmpArray4(17,3)
     tmpArray5(33,2) = Yqc*tmpArray4(19,2) + alphaYpq*TmpArray4(19,3) + TwoTerms(3)
     tmpArray5(34,2) = Yqc*tmpArray4(20,2) + alphaYpq*TmpArray4(20,3)
     tmpArray5(35,2) = Zqc*tmpArray4(20,2) + alphaZpq*TmpArray4(20,3) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray3(5,3) + alphaQ*TmpArray3(5,4))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,3) + alphaQ*TmpArray3(8,4))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,3) + alphaQ*TmpArray3(10,4))
     tmpArray5(21,3) = Xqc*tmpArray4(11,3) + alphaXpq*TmpArray4(11,4) + 3*TwoTerms(1)
     tmpArray5(22,3) = Yqc*tmpArray4(11,3) + alphaYpq*TmpArray4(11,4)
     tmpArray5(23,3) = Zqc*tmpArray4(11,3) + alphaZpq*TmpArray4(11,4)
     tmpArray5(24,3) = Xqc*tmpArray4(14,3) + alphaXpq*TmpArray4(14,4) + TwoTerms(2)
     tmpArray5(25,3) = Yqc*tmpArray4(13,3) + alphaYpq*TmpArray4(13,4)
     tmpArray5(26,3) = Xqc*tmpArray4(16,3) + alphaXpq*TmpArray4(16,4) + TwoTerms(3)
     tmpArray5(27,3) = Xqc*tmpArray4(17,3) + alphaXpq*TmpArray4(17,4)
     tmpArray5(28,3) = Xqc*tmpArray4(18,3) + alphaXpq*TmpArray4(18,4)
     tmpArray5(29,3) = Xqc*tmpArray4(19,3) + alphaXpq*TmpArray4(19,4)
     tmpArray5(30,3) = Xqc*tmpArray4(20,3) + alphaXpq*TmpArray4(20,4)
     tmpArray5(31,3) = Yqc*tmpArray4(17,3) + alphaYpq*TmpArray4(17,4) + 3*TwoTerms(2)
     tmpArray5(32,3) = Zqc*tmpArray4(17,3) + alphaZpq*TmpArray4(17,4)
     tmpArray5(33,3) = Yqc*tmpArray4(19,3) + alphaYpq*TmpArray4(19,4) + TwoTerms(3)
     tmpArray5(34,3) = Yqc*tmpArray4(20,3) + alphaYpq*TmpArray4(20,4)
     tmpArray5(35,3) = Zqc*tmpArray4(20,3) + alphaZpq*TmpArray4(20,4) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(AuxArray(11,IP) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(AuxArray(14,IP) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(AuxArray(16,IP) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(AuxArray(17,IP) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(AuxArray(19,IP) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(AuxArray(20,IP) + alphaQ*TmpArray4(20,2))
     AuxArray(36,iP) = Xqc*AuxArray(21,iP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,iP) = Yqc*AuxArray(21,iP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,iP) = Zqc*AuxArray(21,iP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,iP) = Xqc*AuxArray(24,iP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,iP) = Yqc*AuxArray(23,iP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,iP) = Xqc*AuxArray(26,iP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,iP) = Xqc*AuxArray(27,iP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,iP) = Zqc*AuxArray(24,iP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,iP) = Yqc*AuxArray(26,iP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,iP) = Xqc*AuxArray(30,iP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,iP) = Xqc*AuxArray(31,iP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,iP) = Xqc*AuxArray(32,iP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,iP) = Xqc*AuxArray(33,iP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,iP) = Xqc*AuxArray(34,iP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,iP) = Xqc*AuxArray(35,iP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,iP) = Yqc*AuxArray(31,iP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,iP) = Zqc*AuxArray(31,iP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,iP) = Yqc*AuxArray(33,iP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,iP) = Yqc*AuxArray(34,iP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,iP) = Yqc*AuxArray(35,iP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,iP) = Zqc*AuxArray(35,iP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TmpArray4(11,2) + alphaQ*TmpArray4(11,3))
     TwoTerms(2) = inv2expQ*(TmpArray4(14,2) + alphaQ*TmpArray4(14,3))
     TwoTerms(3) = inv2expQ*(TmpArray4(16,2) + alphaQ*TmpArray4(16,3))
     TwoTerms(4) = inv2expQ*(TmpArray4(17,2) + alphaQ*TmpArray4(17,3))
     TwoTerms(5) = inv2expQ*(TmpArray4(19,2) + alphaQ*TmpArray4(19,3))
     TwoTerms(6) = inv2expQ*(TmpArray4(20,2) + alphaQ*TmpArray4(20,3))
     tmpArray6(36,2) = Xqc*tmpArray5(21,2) + alphaXpq*TmpArray5(21,3) + 4*TwoTerms(1)
     tmpArray6(37,2) = Yqc*tmpArray5(21,2) + alphaYpq*TmpArray5(21,3)
     tmpArray6(38,2) = Zqc*tmpArray5(21,2) + alphaZpq*TmpArray5(21,3)
     tmpArray6(39,2) = Xqc*tmpArray5(24,2) + alphaXpq*TmpArray5(24,3) + 2*TwoTerms(2)
     tmpArray6(40,2) = Yqc*tmpArray5(23,2) + alphaYpq*TmpArray5(23,3)
     tmpArray6(41,2) = Xqc*tmpArray5(26,2) + alphaXpq*TmpArray5(26,3) + 2*TwoTerms(3)
     tmpArray6(42,2) = Xqc*tmpArray5(27,2) + alphaXpq*TmpArray5(27,3) + TwoTerms(4)
     tmpArray6(43,2) = Zqc*tmpArray5(24,2) + alphaZpq*TmpArray5(24,3)
     tmpArray6(44,2) = Yqc*tmpArray5(26,2) + alphaYpq*TmpArray5(26,3)
     tmpArray6(45,2) = Xqc*tmpArray5(30,2) + alphaXpq*TmpArray5(30,3) + TwoTerms(6)
     tmpArray6(46,2) = Xqc*tmpArray5(31,2) + alphaXpq*TmpArray5(31,3)
     tmpArray6(47,2) = Xqc*tmpArray5(32,2) + alphaXpq*TmpArray5(32,3)
     tmpArray6(48,2) = Xqc*tmpArray5(33,2) + alphaXpq*TmpArray5(33,3)
     tmpArray6(49,2) = Xqc*tmpArray5(34,2) + alphaXpq*TmpArray5(34,3)
     tmpArray6(50,2) = Xqc*tmpArray5(35,2) + alphaXpq*TmpArray5(35,3)
     tmpArray6(51,2) = Yqc*tmpArray5(31,2) + alphaYpq*TmpArray5(31,3) + 4*TwoTerms(4)
     tmpArray6(52,2) = Zqc*tmpArray5(31,2) + alphaZpq*TmpArray5(31,3)
     tmpArray6(53,2) = Yqc*tmpArray5(33,2) + alphaYpq*TmpArray5(33,3) + 2*TwoTerms(5)
     tmpArray6(54,2) = Yqc*tmpArray5(34,2) + alphaYpq*TmpArray5(34,3) + TwoTerms(6)
     tmpArray6(55,2) = Yqc*tmpArray5(35,2) + alphaYpq*TmpArray5(35,3)
     tmpArray6(56,2) = Zqc*tmpArray5(35,2) + alphaZpq*TmpArray5(35,3) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(AuxArray(21,IP) + alphaQ*TmpArray5(21,2))
     TwoTerms(2) = inv2expQ*(AuxArray(24,IP) + alphaQ*TmpArray5(24,2))
     TwoTerms(3) = inv2expQ*(AuxArray(26,IP) + alphaQ*TmpArray5(26,2))
     TwoTerms(4) = inv2expQ*(AuxArray(27,IP) + alphaQ*TmpArray5(27,2))
     TwoTerms(5) = inv2expQ*(AuxArray(30,IP) + alphaQ*TmpArray5(30,2))
     TwoTerms(6) = inv2expQ*(AuxArray(31,IP) + alphaQ*TmpArray5(31,2))
     TwoTerms(7) = inv2expQ*(AuxArray(33,IP) + alphaQ*TmpArray5(33,2))
     TwoTerms(8) = inv2expQ*(AuxArray(34,IP) + alphaQ*TmpArray5(34,2))
     TwoTerms(9) = inv2expQ*(AuxArray(35,IP) + alphaQ*TmpArray5(35,2))
     AuxArray(57,iP) = Xqc*AuxArray(36,iP) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,iP) = Yqc*AuxArray(36,iP) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,iP) = Zqc*AuxArray(36,iP) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,iP) = Xqc*AuxArray(39,iP) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,iP) = Yqc*AuxArray(38,iP) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,iP) = Xqc*AuxArray(41,iP) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,iP) = Xqc*AuxArray(42,iP) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,iP) = Zqc*AuxArray(39,iP) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,iP) = Yqc*AuxArray(41,iP) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,iP) = Xqc*AuxArray(45,iP) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,iP) = Xqc*AuxArray(46,iP) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,iP) = Zqc*AuxArray(42,iP) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,iP) = Xqc*AuxArray(48,iP) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,iP) = Yqc*AuxArray(45,iP) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,iP) = Xqc*AuxArray(50,iP) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,iP) = Xqc*AuxArray(51,iP) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,iP) = Xqc*AuxArray(52,iP) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,iP) = Xqc*AuxArray(53,iP) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,iP) = Xqc*AuxArray(54,iP) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,iP) = Xqc*AuxArray(55,iP) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,iP) = Xqc*AuxArray(56,iP) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,iP) = Yqc*AuxArray(51,iP) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,iP) = Zqc*AuxArray(51,iP) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,iP) = Yqc*AuxArray(53,iP) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,iP) = Yqc*AuxArray(54,iP) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,iP) = Yqc*AuxArray(55,iP) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,iP) = Yqc*AuxArray(56,iP) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,iP) = Zqc*AuxArray(56,iP) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine

subroutine VerticalRecurrenceGPUGen7C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 7,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(  120,nPrimQ*nPrimP*nPassP)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ
  real(realk) :: TwoTerms(  21)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:8)
  real(realk) :: TMParray2(  2:  4,2:7)
  real(realk) :: TMParray3(  5: 10,2:6)
  real(realk) :: TMParray4( 11: 20,2:5)
  real(realk) :: TMParray5( 21: 35,2:4)
  real(realk) :: TMParray6( 36: 56,2:3)
  real(realk) :: TMParray7( 57: 84,2:2)
  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqc,Yqc,Zqc,alphaQ,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpQ,inv2expQ,&
!$ACC         PREF,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TMParray4,&
!$ACC         TMParray5,&
!$ACC         TMParray6,&
!$ACC         TMParray7,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Ccenter, &
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
     Xqc = Qcent(1,iPrimQ) - Ccenter(1)
     Yqc = Qcent(2,iPrimQ) - Ccenter(2)
     Zqc = Qcent(3,iPrimQ) - Ccenter(3)
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
     TMParray1(1, 7) = PREF*RJ000Array( 6,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 8) = PREF*RJ000Array( 7,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xqc*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Yqc*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zqc*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqc*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqc*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqc*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqc*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqc*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqc*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xqc*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Yqc*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zqc*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     tmpArray2(2,5) = Xqc*tmpArray1(1,5) + alphaXpq*TmpArray1(1,6)
     tmpArray2(3,5) = Yqc*tmpArray1(1,5) + alphaYpq*TmpArray1(1,6)
     tmpArray2(4,5) = Zqc*tmpArray1(1,5) + alphaZpq*TmpArray1(1,6)
     tmpArray2(2,6) = Xqc*tmpArray1(1,6) + alphaXpq*TmpArray1(1,7)
     tmpArray2(3,6) = Yqc*tmpArray1(1,6) + alphaYpq*TmpArray1(1,7)
     tmpArray2(4,6) = Zqc*tmpArray1(1,6) + alphaZpq*TmpArray1(1,7)
     tmpArray2(2,7) = Xqc*tmpArray1(1,7) + alphaXpq*TmpArray1(1,8)
     tmpArray2(3,7) = Yqc*tmpArray1(1,7) + alphaYpq*TmpArray1(1,8)
     tmpArray2(4,7) = Zqc*tmpArray1(1,7) + alphaZpq*TmpArray1(1,8)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,iP) = Xqc*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xqc*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xqc*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Yqc*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Yqc*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zqc*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,2) + alphaQ*TmpArray1(1,3))
     tmpArray3(5,2) = Xqc*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xqc*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xqc*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Yqc*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Yqc*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zqc*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,3) + alphaQ*TmpArray1(1,4))
     tmpArray3(5,3) = Xqc*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xqc*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xqc*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Yqc*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Yqc*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zqc*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,4) + alphaQ*TmpArray1(1,5))
     tmpArray3(5,4) = Xqc*tmpArray2(2,4) + alphaXpq*TmpArray2(2,5) + TwoTerms(1)
     tmpArray3(6,4) = Xqc*tmpArray2(3,4) + alphaXpq*TmpArray2(3,5)
     tmpArray3(7,4) = Xqc*tmpArray2(4,4) + alphaXpq*TmpArray2(4,5)
     tmpArray3(8,4) = Yqc*tmpArray2(3,4) + alphaYpq*TmpArray2(3,5) + TwoTerms(1)
     tmpArray3(9,4) = Yqc*tmpArray2(4,4) + alphaYpq*TmpArray2(4,5)
     tmpArray3(10,4) = Zqc*tmpArray2(4,4) + alphaZpq*TmpArray2(4,5) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,5) + alphaQ*TmpArray1(1,6))
     tmpArray3(5,5) = Xqc*tmpArray2(2,5) + alphaXpq*TmpArray2(2,6) + TwoTerms(1)
     tmpArray3(6,5) = Xqc*tmpArray2(3,5) + alphaXpq*TmpArray2(3,6)
     tmpArray3(7,5) = Xqc*tmpArray2(4,5) + alphaXpq*TmpArray2(4,6)
     tmpArray3(8,5) = Yqc*tmpArray2(3,5) + alphaYpq*TmpArray2(3,6) + TwoTerms(1)
     tmpArray3(9,5) = Yqc*tmpArray2(4,5) + alphaYpq*TmpArray2(4,6)
     tmpArray3(10,5) = Zqc*tmpArray2(4,5) + alphaZpq*TmpArray2(4,6) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,6) + alphaQ*TmpArray1(1,7))
     tmpArray3(5,6) = Xqc*tmpArray2(2,6) + alphaXpq*TmpArray2(2,7) + TwoTerms(1)
     tmpArray3(6,6) = Xqc*tmpArray2(3,6) + alphaXpq*TmpArray2(3,7)
     tmpArray3(7,6) = Xqc*tmpArray2(4,6) + alphaXpq*TmpArray2(4,7)
     tmpArray3(8,6) = Yqc*tmpArray2(3,6) + alphaYpq*TmpArray2(3,7) + TwoTerms(1)
     tmpArray3(9,6) = Yqc*tmpArray2(4,6) + alphaYpq*TmpArray2(4,7)
     tmpArray3(10,6) = Zqc*tmpArray2(4,6) + alphaZpq*TmpArray2(4,7) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(AuxArray(2,IP) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(3,IP) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(4,IP) + alphaQ*TmpArray2(4,2))
     AuxArray(11,iP) = Xqc*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Yqc*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zqc*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xqc*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xqc*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xqc*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Yqc*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zqc*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Yqc*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zqc*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,2) + alphaQ*TmpArray2(2,3))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,2) + alphaQ*TmpArray2(3,3))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,2) + alphaQ*TmpArray2(4,3))
     tmpArray4(11,2) = Xqc*tmpArray3(5,2) + alphaXpq*TmpArray3(5,3) + 2*TwoTerms(1)
     tmpArray4(12,2) = Yqc*tmpArray3(5,2) + alphaYpq*TmpArray3(5,3)
     tmpArray4(13,2) = Zqc*tmpArray3(5,2) + alphaZpq*TmpArray3(5,3)
     tmpArray4(14,2) = Xqc*tmpArray3(8,2) + alphaXpq*TmpArray3(8,3)
     tmpArray4(15,2) = Xqc*tmpArray3(9,2) + alphaXpq*TmpArray3(9,3)
     tmpArray4(16,2) = Xqc*tmpArray3(10,2) + alphaXpq*TmpArray3(10,3)
     tmpArray4(17,2) = Yqc*tmpArray3(8,2) + alphaYpq*TmpArray3(8,3) + 2*TwoTerms(2)
     tmpArray4(18,2) = Zqc*tmpArray3(8,2) + alphaZpq*TmpArray3(8,3)
     tmpArray4(19,2) = Yqc*tmpArray3(10,2) + alphaYpq*TmpArray3(10,3)
     tmpArray4(20,2) = Zqc*tmpArray3(10,2) + alphaZpq*TmpArray3(10,3) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,3) + alphaQ*TmpArray2(2,4))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,3) + alphaQ*TmpArray2(3,4))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,3) + alphaQ*TmpArray2(4,4))
     tmpArray4(11,3) = Xqc*tmpArray3(5,3) + alphaXpq*TmpArray3(5,4) + 2*TwoTerms(1)
     tmpArray4(12,3) = Yqc*tmpArray3(5,3) + alphaYpq*TmpArray3(5,4)
     tmpArray4(13,3) = Zqc*tmpArray3(5,3) + alphaZpq*TmpArray3(5,4)
     tmpArray4(14,3) = Xqc*tmpArray3(8,3) + alphaXpq*TmpArray3(8,4)
     tmpArray4(15,3) = Xqc*tmpArray3(9,3) + alphaXpq*TmpArray3(9,4)
     tmpArray4(16,3) = Xqc*tmpArray3(10,3) + alphaXpq*TmpArray3(10,4)
     tmpArray4(17,3) = Yqc*tmpArray3(8,3) + alphaYpq*TmpArray3(8,4) + 2*TwoTerms(2)
     tmpArray4(18,3) = Zqc*tmpArray3(8,3) + alphaZpq*TmpArray3(8,4)
     tmpArray4(19,3) = Yqc*tmpArray3(10,3) + alphaYpq*TmpArray3(10,4)
     tmpArray4(20,3) = Zqc*tmpArray3(10,3) + alphaZpq*TmpArray3(10,4) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,4) + alphaQ*TmpArray2(2,5))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,4) + alphaQ*TmpArray2(3,5))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,4) + alphaQ*TmpArray2(4,5))
     tmpArray4(11,4) = Xqc*tmpArray3(5,4) + alphaXpq*TmpArray3(5,5) + 2*TwoTerms(1)
     tmpArray4(12,4) = Yqc*tmpArray3(5,4) + alphaYpq*TmpArray3(5,5)
     tmpArray4(13,4) = Zqc*tmpArray3(5,4) + alphaZpq*TmpArray3(5,5)
     tmpArray4(14,4) = Xqc*tmpArray3(8,4) + alphaXpq*TmpArray3(8,5)
     tmpArray4(15,4) = Xqc*tmpArray3(9,4) + alphaXpq*TmpArray3(9,5)
     tmpArray4(16,4) = Xqc*tmpArray3(10,4) + alphaXpq*TmpArray3(10,5)
     tmpArray4(17,4) = Yqc*tmpArray3(8,4) + alphaYpq*TmpArray3(8,5) + 2*TwoTerms(2)
     tmpArray4(18,4) = Zqc*tmpArray3(8,4) + alphaZpq*TmpArray3(8,5)
     tmpArray4(19,4) = Yqc*tmpArray3(10,4) + alphaYpq*TmpArray3(10,5)
     tmpArray4(20,4) = Zqc*tmpArray3(10,4) + alphaZpq*TmpArray3(10,5) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,5) + alphaQ*TmpArray2(2,6))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,5) + alphaQ*TmpArray2(3,6))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,5) + alphaQ*TmpArray2(4,6))
     tmpArray4(11,5) = Xqc*tmpArray3(5,5) + alphaXpq*TmpArray3(5,6) + 2*TwoTerms(1)
     tmpArray4(12,5) = Yqc*tmpArray3(5,5) + alphaYpq*TmpArray3(5,6)
     tmpArray4(13,5) = Zqc*tmpArray3(5,5) + alphaZpq*TmpArray3(5,6)
     tmpArray4(14,5) = Xqc*tmpArray3(8,5) + alphaXpq*TmpArray3(8,6)
     tmpArray4(15,5) = Xqc*tmpArray3(9,5) + alphaXpq*TmpArray3(9,6)
     tmpArray4(16,5) = Xqc*tmpArray3(10,5) + alphaXpq*TmpArray3(10,6)
     tmpArray4(17,5) = Yqc*tmpArray3(8,5) + alphaYpq*TmpArray3(8,6) + 2*TwoTerms(2)
     tmpArray4(18,5) = Zqc*tmpArray3(8,5) + alphaZpq*TmpArray3(8,6)
     tmpArray4(19,5) = Yqc*tmpArray3(10,5) + alphaYpq*TmpArray3(10,6)
     tmpArray4(20,5) = Zqc*tmpArray3(10,5) + alphaZpq*TmpArray3(10,6) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(AuxArray(5,IP) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(AuxArray(8,IP) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(AuxArray(10,IP) + alphaQ*TmpArray3(10,2))
     AuxArray(21,iP) = Xqc*AuxArray(11,iP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,iP) = Yqc*AuxArray(11,iP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,iP) = Zqc*AuxArray(11,iP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,iP) = Xqc*AuxArray(14,iP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,iP) = Yqc*AuxArray(13,iP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,iP) = Xqc*AuxArray(16,iP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,iP) = Xqc*AuxArray(17,iP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,iP) = Xqc*AuxArray(18,iP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,iP) = Xqc*AuxArray(19,iP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,iP) = Xqc*AuxArray(20,iP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,iP) = Yqc*AuxArray(17,iP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,iP) = Zqc*AuxArray(17,iP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,iP) = Yqc*AuxArray(19,iP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,iP) = Yqc*AuxArray(20,iP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,iP) = Zqc*AuxArray(20,iP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray3(5,2) + alphaQ*TmpArray3(5,3))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,2) + alphaQ*TmpArray3(8,3))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,2) + alphaQ*TmpArray3(10,3))
     tmpArray5(21,2) = Xqc*tmpArray4(11,2) + alphaXpq*TmpArray4(11,3) + 3*TwoTerms(1)
     tmpArray5(22,2) = Yqc*tmpArray4(11,2) + alphaYpq*TmpArray4(11,3)
     tmpArray5(23,2) = Zqc*tmpArray4(11,2) + alphaZpq*TmpArray4(11,3)
     tmpArray5(24,2) = Xqc*tmpArray4(14,2) + alphaXpq*TmpArray4(14,3) + TwoTerms(2)
     tmpArray5(25,2) = Yqc*tmpArray4(13,2) + alphaYpq*TmpArray4(13,3)
     tmpArray5(26,2) = Xqc*tmpArray4(16,2) + alphaXpq*TmpArray4(16,3) + TwoTerms(3)
     tmpArray5(27,2) = Xqc*tmpArray4(17,2) + alphaXpq*TmpArray4(17,3)
     tmpArray5(28,2) = Xqc*tmpArray4(18,2) + alphaXpq*TmpArray4(18,3)
     tmpArray5(29,2) = Xqc*tmpArray4(19,2) + alphaXpq*TmpArray4(19,3)
     tmpArray5(30,2) = Xqc*tmpArray4(20,2) + alphaXpq*TmpArray4(20,3)
     tmpArray5(31,2) = Yqc*tmpArray4(17,2) + alphaYpq*TmpArray4(17,3) + 3*TwoTerms(2)
     tmpArray5(32,2) = Zqc*tmpArray4(17,2) + alphaZpq*TmpArray4(17,3)
     tmpArray5(33,2) = Yqc*tmpArray4(19,2) + alphaYpq*TmpArray4(19,3) + TwoTerms(3)
     tmpArray5(34,2) = Yqc*tmpArray4(20,2) + alphaYpq*TmpArray4(20,3)
     tmpArray5(35,2) = Zqc*tmpArray4(20,2) + alphaZpq*TmpArray4(20,3) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray3(5,3) + alphaQ*TmpArray3(5,4))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,3) + alphaQ*TmpArray3(8,4))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,3) + alphaQ*TmpArray3(10,4))
     tmpArray5(21,3) = Xqc*tmpArray4(11,3) + alphaXpq*TmpArray4(11,4) + 3*TwoTerms(1)
     tmpArray5(22,3) = Yqc*tmpArray4(11,3) + alphaYpq*TmpArray4(11,4)
     tmpArray5(23,3) = Zqc*tmpArray4(11,3) + alphaZpq*TmpArray4(11,4)
     tmpArray5(24,3) = Xqc*tmpArray4(14,3) + alphaXpq*TmpArray4(14,4) + TwoTerms(2)
     tmpArray5(25,3) = Yqc*tmpArray4(13,3) + alphaYpq*TmpArray4(13,4)
     tmpArray5(26,3) = Xqc*tmpArray4(16,3) + alphaXpq*TmpArray4(16,4) + TwoTerms(3)
     tmpArray5(27,3) = Xqc*tmpArray4(17,3) + alphaXpq*TmpArray4(17,4)
     tmpArray5(28,3) = Xqc*tmpArray4(18,3) + alphaXpq*TmpArray4(18,4)
     tmpArray5(29,3) = Xqc*tmpArray4(19,3) + alphaXpq*TmpArray4(19,4)
     tmpArray5(30,3) = Xqc*tmpArray4(20,3) + alphaXpq*TmpArray4(20,4)
     tmpArray5(31,3) = Yqc*tmpArray4(17,3) + alphaYpq*TmpArray4(17,4) + 3*TwoTerms(2)
     tmpArray5(32,3) = Zqc*tmpArray4(17,3) + alphaZpq*TmpArray4(17,4)
     tmpArray5(33,3) = Yqc*tmpArray4(19,3) + alphaYpq*TmpArray4(19,4) + TwoTerms(3)
     tmpArray5(34,3) = Yqc*tmpArray4(20,3) + alphaYpq*TmpArray4(20,4)
     tmpArray5(35,3) = Zqc*tmpArray4(20,3) + alphaZpq*TmpArray4(20,4) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray3(5,4) + alphaQ*TmpArray3(5,5))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,4) + alphaQ*TmpArray3(8,5))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,4) + alphaQ*TmpArray3(10,5))
     tmpArray5(21,4) = Xqc*tmpArray4(11,4) + alphaXpq*TmpArray4(11,5) + 3*TwoTerms(1)
     tmpArray5(22,4) = Yqc*tmpArray4(11,4) + alphaYpq*TmpArray4(11,5)
     tmpArray5(23,4) = Zqc*tmpArray4(11,4) + alphaZpq*TmpArray4(11,5)
     tmpArray5(24,4) = Xqc*tmpArray4(14,4) + alphaXpq*TmpArray4(14,5) + TwoTerms(2)
     tmpArray5(25,4) = Yqc*tmpArray4(13,4) + alphaYpq*TmpArray4(13,5)
     tmpArray5(26,4) = Xqc*tmpArray4(16,4) + alphaXpq*TmpArray4(16,5) + TwoTerms(3)
     tmpArray5(27,4) = Xqc*tmpArray4(17,4) + alphaXpq*TmpArray4(17,5)
     tmpArray5(28,4) = Xqc*tmpArray4(18,4) + alphaXpq*TmpArray4(18,5)
     tmpArray5(29,4) = Xqc*tmpArray4(19,4) + alphaXpq*TmpArray4(19,5)
     tmpArray5(30,4) = Xqc*tmpArray4(20,4) + alphaXpq*TmpArray4(20,5)
     tmpArray5(31,4) = Yqc*tmpArray4(17,4) + alphaYpq*TmpArray4(17,5) + 3*TwoTerms(2)
     tmpArray5(32,4) = Zqc*tmpArray4(17,4) + alphaZpq*TmpArray4(17,5)
     tmpArray5(33,4) = Yqc*tmpArray4(19,4) + alphaYpq*TmpArray4(19,5) + TwoTerms(3)
     tmpArray5(34,4) = Yqc*tmpArray4(20,4) + alphaYpq*TmpArray4(20,5)
     tmpArray5(35,4) = Zqc*tmpArray4(20,4) + alphaZpq*TmpArray4(20,5) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(AuxArray(11,IP) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(AuxArray(14,IP) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(AuxArray(16,IP) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(AuxArray(17,IP) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(AuxArray(19,IP) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(AuxArray(20,IP) + alphaQ*TmpArray4(20,2))
     AuxArray(36,iP) = Xqc*AuxArray(21,iP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,iP) = Yqc*AuxArray(21,iP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,iP) = Zqc*AuxArray(21,iP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,iP) = Xqc*AuxArray(24,iP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,iP) = Yqc*AuxArray(23,iP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,iP) = Xqc*AuxArray(26,iP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,iP) = Xqc*AuxArray(27,iP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,iP) = Zqc*AuxArray(24,iP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,iP) = Yqc*AuxArray(26,iP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,iP) = Xqc*AuxArray(30,iP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,iP) = Xqc*AuxArray(31,iP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,iP) = Xqc*AuxArray(32,iP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,iP) = Xqc*AuxArray(33,iP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,iP) = Xqc*AuxArray(34,iP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,iP) = Xqc*AuxArray(35,iP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,iP) = Yqc*AuxArray(31,iP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,iP) = Zqc*AuxArray(31,iP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,iP) = Yqc*AuxArray(33,iP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,iP) = Yqc*AuxArray(34,iP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,iP) = Yqc*AuxArray(35,iP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,iP) = Zqc*AuxArray(35,iP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TmpArray4(11,2) + alphaQ*TmpArray4(11,3))
     TwoTerms(2) = inv2expQ*(TmpArray4(14,2) + alphaQ*TmpArray4(14,3))
     TwoTerms(3) = inv2expQ*(TmpArray4(16,2) + alphaQ*TmpArray4(16,3))
     TwoTerms(4) = inv2expQ*(TmpArray4(17,2) + alphaQ*TmpArray4(17,3))
     TwoTerms(5) = inv2expQ*(TmpArray4(19,2) + alphaQ*TmpArray4(19,3))
     TwoTerms(6) = inv2expQ*(TmpArray4(20,2) + alphaQ*TmpArray4(20,3))
     tmpArray6(36,2) = Xqc*tmpArray5(21,2) + alphaXpq*TmpArray5(21,3) + 4*TwoTerms(1)
     tmpArray6(37,2) = Yqc*tmpArray5(21,2) + alphaYpq*TmpArray5(21,3)
     tmpArray6(38,2) = Zqc*tmpArray5(21,2) + alphaZpq*TmpArray5(21,3)
     tmpArray6(39,2) = Xqc*tmpArray5(24,2) + alphaXpq*TmpArray5(24,3) + 2*TwoTerms(2)
     tmpArray6(40,2) = Yqc*tmpArray5(23,2) + alphaYpq*TmpArray5(23,3)
     tmpArray6(41,2) = Xqc*tmpArray5(26,2) + alphaXpq*TmpArray5(26,3) + 2*TwoTerms(3)
     tmpArray6(42,2) = Xqc*tmpArray5(27,2) + alphaXpq*TmpArray5(27,3) + TwoTerms(4)
     tmpArray6(43,2) = Zqc*tmpArray5(24,2) + alphaZpq*TmpArray5(24,3)
     tmpArray6(44,2) = Yqc*tmpArray5(26,2) + alphaYpq*TmpArray5(26,3)
     tmpArray6(45,2) = Xqc*tmpArray5(30,2) + alphaXpq*TmpArray5(30,3) + TwoTerms(6)
     tmpArray6(46,2) = Xqc*tmpArray5(31,2) + alphaXpq*TmpArray5(31,3)
     tmpArray6(47,2) = Xqc*tmpArray5(32,2) + alphaXpq*TmpArray5(32,3)
     tmpArray6(48,2) = Xqc*tmpArray5(33,2) + alphaXpq*TmpArray5(33,3)
     tmpArray6(49,2) = Xqc*tmpArray5(34,2) + alphaXpq*TmpArray5(34,3)
     tmpArray6(50,2) = Xqc*tmpArray5(35,2) + alphaXpq*TmpArray5(35,3)
     tmpArray6(51,2) = Yqc*tmpArray5(31,2) + alphaYpq*TmpArray5(31,3) + 4*TwoTerms(4)
     tmpArray6(52,2) = Zqc*tmpArray5(31,2) + alphaZpq*TmpArray5(31,3)
     tmpArray6(53,2) = Yqc*tmpArray5(33,2) + alphaYpq*TmpArray5(33,3) + 2*TwoTerms(5)
     tmpArray6(54,2) = Yqc*tmpArray5(34,2) + alphaYpq*TmpArray5(34,3) + TwoTerms(6)
     tmpArray6(55,2) = Yqc*tmpArray5(35,2) + alphaYpq*TmpArray5(35,3)
     tmpArray6(56,2) = Zqc*tmpArray5(35,2) + alphaZpq*TmpArray5(35,3) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TmpArray4(11,3) + alphaQ*TmpArray4(11,4))
     TwoTerms(2) = inv2expQ*(TmpArray4(14,3) + alphaQ*TmpArray4(14,4))
     TwoTerms(3) = inv2expQ*(TmpArray4(16,3) + alphaQ*TmpArray4(16,4))
     TwoTerms(4) = inv2expQ*(TmpArray4(17,3) + alphaQ*TmpArray4(17,4))
     TwoTerms(5) = inv2expQ*(TmpArray4(19,3) + alphaQ*TmpArray4(19,4))
     TwoTerms(6) = inv2expQ*(TmpArray4(20,3) + alphaQ*TmpArray4(20,4))
     tmpArray6(36,3) = Xqc*tmpArray5(21,3) + alphaXpq*TmpArray5(21,4) + 4*TwoTerms(1)
     tmpArray6(37,3) = Yqc*tmpArray5(21,3) + alphaYpq*TmpArray5(21,4)
     tmpArray6(38,3) = Zqc*tmpArray5(21,3) + alphaZpq*TmpArray5(21,4)
     tmpArray6(39,3) = Xqc*tmpArray5(24,3) + alphaXpq*TmpArray5(24,4) + 2*TwoTerms(2)
     tmpArray6(40,3) = Yqc*tmpArray5(23,3) + alphaYpq*TmpArray5(23,4)
     tmpArray6(41,3) = Xqc*tmpArray5(26,3) + alphaXpq*TmpArray5(26,4) + 2*TwoTerms(3)
     tmpArray6(42,3) = Xqc*tmpArray5(27,3) + alphaXpq*TmpArray5(27,4) + TwoTerms(4)
     tmpArray6(43,3) = Zqc*tmpArray5(24,3) + alphaZpq*TmpArray5(24,4)
     tmpArray6(44,3) = Yqc*tmpArray5(26,3) + alphaYpq*TmpArray5(26,4)
     tmpArray6(45,3) = Xqc*tmpArray5(30,3) + alphaXpq*TmpArray5(30,4) + TwoTerms(6)
     tmpArray6(46,3) = Xqc*tmpArray5(31,3) + alphaXpq*TmpArray5(31,4)
     tmpArray6(47,3) = Xqc*tmpArray5(32,3) + alphaXpq*TmpArray5(32,4)
     tmpArray6(48,3) = Xqc*tmpArray5(33,3) + alphaXpq*TmpArray5(33,4)
     tmpArray6(49,3) = Xqc*tmpArray5(34,3) + alphaXpq*TmpArray5(34,4)
     tmpArray6(50,3) = Xqc*tmpArray5(35,3) + alphaXpq*TmpArray5(35,4)
     tmpArray6(51,3) = Yqc*tmpArray5(31,3) + alphaYpq*TmpArray5(31,4) + 4*TwoTerms(4)
     tmpArray6(52,3) = Zqc*tmpArray5(31,3) + alphaZpq*TmpArray5(31,4)
     tmpArray6(53,3) = Yqc*tmpArray5(33,3) + alphaYpq*TmpArray5(33,4) + 2*TwoTerms(5)
     tmpArray6(54,3) = Yqc*tmpArray5(34,3) + alphaYpq*TmpArray5(34,4) + TwoTerms(6)
     tmpArray6(55,3) = Yqc*tmpArray5(35,3) + alphaYpq*TmpArray5(35,4)
     tmpArray6(56,3) = Zqc*tmpArray5(35,3) + alphaZpq*TmpArray5(35,4) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(AuxArray(21,IP) + alphaQ*TmpArray5(21,2))
     TwoTerms(2) = inv2expQ*(AuxArray(24,IP) + alphaQ*TmpArray5(24,2))
     TwoTerms(3) = inv2expQ*(AuxArray(26,IP) + alphaQ*TmpArray5(26,2))
     TwoTerms(4) = inv2expQ*(AuxArray(27,IP) + alphaQ*TmpArray5(27,2))
     TwoTerms(5) = inv2expQ*(AuxArray(30,IP) + alphaQ*TmpArray5(30,2))
     TwoTerms(6) = inv2expQ*(AuxArray(31,IP) + alphaQ*TmpArray5(31,2))
     TwoTerms(7) = inv2expQ*(AuxArray(33,IP) + alphaQ*TmpArray5(33,2))
     TwoTerms(8) = inv2expQ*(AuxArray(34,IP) + alphaQ*TmpArray5(34,2))
     TwoTerms(9) = inv2expQ*(AuxArray(35,IP) + alphaQ*TmpArray5(35,2))
     AuxArray(57,iP) = Xqc*AuxArray(36,iP) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,iP) = Yqc*AuxArray(36,iP) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,iP) = Zqc*AuxArray(36,iP) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,iP) = Xqc*AuxArray(39,iP) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,iP) = Yqc*AuxArray(38,iP) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,iP) = Xqc*AuxArray(41,iP) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,iP) = Xqc*AuxArray(42,iP) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,iP) = Zqc*AuxArray(39,iP) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,iP) = Yqc*AuxArray(41,iP) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,iP) = Xqc*AuxArray(45,iP) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,iP) = Xqc*AuxArray(46,iP) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,iP) = Zqc*AuxArray(42,iP) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,iP) = Xqc*AuxArray(48,iP) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,iP) = Yqc*AuxArray(45,iP) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,iP) = Xqc*AuxArray(50,iP) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,iP) = Xqc*AuxArray(51,iP) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,iP) = Xqc*AuxArray(52,iP) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,iP) = Xqc*AuxArray(53,iP) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,iP) = Xqc*AuxArray(54,iP) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,iP) = Xqc*AuxArray(55,iP) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,iP) = Xqc*AuxArray(56,iP) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,iP) = Yqc*AuxArray(51,iP) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,iP) = Zqc*AuxArray(51,iP) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,iP) = Yqc*AuxArray(53,iP) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,iP) = Yqc*AuxArray(54,iP) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,iP) = Yqc*AuxArray(55,iP) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,iP) = Yqc*AuxArray(56,iP) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,iP) = Zqc*AuxArray(56,iP) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expQ*(TmpArray5(21,2) + alphaQ*TmpArray5(21,3))
     TwoTerms(2) = inv2expQ*(TmpArray5(24,2) + alphaQ*TmpArray5(24,3))
     TwoTerms(3) = inv2expQ*(TmpArray5(26,2) + alphaQ*TmpArray5(26,3))
     TwoTerms(4) = inv2expQ*(TmpArray5(27,2) + alphaQ*TmpArray5(27,3))
     TwoTerms(5) = inv2expQ*(TmpArray5(30,2) + alphaQ*TmpArray5(30,3))
     TwoTerms(6) = inv2expQ*(TmpArray5(31,2) + alphaQ*TmpArray5(31,3))
     TwoTerms(7) = inv2expQ*(TmpArray5(33,2) + alphaQ*TmpArray5(33,3))
     TwoTerms(8) = inv2expQ*(TmpArray5(34,2) + alphaQ*TmpArray5(34,3))
     TwoTerms(9) = inv2expQ*(TmpArray5(35,2) + alphaQ*TmpArray5(35,3))
     tmpArray7(57,2) = Xqc*tmpArray6(36,2) + alphaXpq*TmpArray6(36,3) + 5*TwoTerms(1)
     tmpArray7(58,2) = Yqc*tmpArray6(36,2) + alphaYpq*TmpArray6(36,3)
     tmpArray7(59,2) = Zqc*tmpArray6(36,2) + alphaZpq*TmpArray6(36,3)
     tmpArray7(60,2) = Xqc*tmpArray6(39,2) + alphaXpq*TmpArray6(39,3) + 3*TwoTerms(2)
     tmpArray7(61,2) = Yqc*tmpArray6(38,2) + alphaYpq*TmpArray6(38,3)
     tmpArray7(62,2) = Xqc*tmpArray6(41,2) + alphaXpq*TmpArray6(41,3) + 3*TwoTerms(3)
     tmpArray7(63,2) = Xqc*tmpArray6(42,2) + alphaXpq*TmpArray6(42,3) + 2*TwoTerms(4)
     tmpArray7(64,2) = Zqc*tmpArray6(39,2) + alphaZpq*TmpArray6(39,3)
     tmpArray7(65,2) = Yqc*tmpArray6(41,2) + alphaYpq*TmpArray6(41,3)
     tmpArray7(66,2) = Xqc*tmpArray6(45,2) + alphaXpq*TmpArray6(45,3) + 2*TwoTerms(5)
     tmpArray7(67,2) = Xqc*tmpArray6(46,2) + alphaXpq*TmpArray6(46,3) + TwoTerms(6)
     tmpArray7(68,2) = Zqc*tmpArray6(42,2) + alphaZpq*TmpArray6(42,3)
     tmpArray7(69,2) = Xqc*tmpArray6(48,2) + alphaXpq*TmpArray6(48,3) + TwoTerms(7)
     tmpArray7(70,2) = Yqc*tmpArray6(45,2) + alphaYpq*TmpArray6(45,3)
     tmpArray7(71,2) = Xqc*tmpArray6(50,2) + alphaXpq*TmpArray6(50,3) + TwoTerms(9)
     tmpArray7(72,2) = Xqc*tmpArray6(51,2) + alphaXpq*TmpArray6(51,3)
     tmpArray7(73,2) = Xqc*tmpArray6(52,2) + alphaXpq*TmpArray6(52,3)
     tmpArray7(74,2) = Xqc*tmpArray6(53,2) + alphaXpq*TmpArray6(53,3)
     tmpArray7(75,2) = Xqc*tmpArray6(54,2) + alphaXpq*TmpArray6(54,3)
     tmpArray7(76,2) = Xqc*tmpArray6(55,2) + alphaXpq*TmpArray6(55,3)
     tmpArray7(77,2) = Xqc*tmpArray6(56,2) + alphaXpq*TmpArray6(56,3)
     tmpArray7(78,2) = Yqc*tmpArray6(51,2) + alphaYpq*TmpArray6(51,3) + 5*TwoTerms(6)
     tmpArray7(79,2) = Zqc*tmpArray6(51,2) + alphaZpq*TmpArray6(51,3)
     tmpArray7(80,2) = Yqc*tmpArray6(53,2) + alphaYpq*TmpArray6(53,3) + 3*TwoTerms(7)
     tmpArray7(81,2) = Yqc*tmpArray6(54,2) + alphaYpq*TmpArray6(54,3) + 2*TwoTerms(8)
     tmpArray7(82,2) = Yqc*tmpArray6(55,2) + alphaYpq*TmpArray6(55,3) + TwoTerms(9)
     tmpArray7(83,2) = Yqc*tmpArray6(56,2) + alphaYpq*TmpArray6(56,3)
     tmpArray7(84,2) = Zqc*tmpArray6(56,2) + alphaZpq*TmpArray6(56,3) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expQ*(AuxArray(36,IP) + alphaQ*TmpArray6(36,2))
     TwoTerms(2) = inv2expQ*(AuxArray(39,IP) + alphaQ*TmpArray6(39,2))
     TwoTerms(3) = inv2expQ*(AuxArray(41,IP) + alphaQ*TmpArray6(41,2))
     TwoTerms(4) = inv2expQ*(AuxArray(42,IP) + alphaQ*TmpArray6(42,2))
     TwoTerms(5) = inv2expQ*(AuxArray(45,IP) + alphaQ*TmpArray6(45,2))
     TwoTerms(6) = inv2expQ*(AuxArray(46,IP) + alphaQ*TmpArray6(46,2))
     TwoTerms(7) = inv2expQ*(AuxArray(48,IP) + alphaQ*TmpArray6(48,2))
     TwoTerms(8) = inv2expQ*(AuxArray(50,IP) + alphaQ*TmpArray6(50,2))
     TwoTerms(9) = inv2expQ*(AuxArray(51,IP) + alphaQ*TmpArray6(51,2))
     TwoTerms(10) = inv2expQ*(AuxArray(53,IP) + alphaQ*TmpArray6(53,2))
     TwoTerms(11) = inv2expQ*(AuxArray(54,IP) + alphaQ*TmpArray6(54,2))
     TwoTerms(12) = inv2expQ*(AuxArray(55,IP) + alphaQ*TmpArray6(55,2))
     TwoTerms(13) = inv2expQ*(AuxArray(56,IP) + alphaQ*TmpArray6(56,2))
     AuxArray(85,iP) = Xqc*AuxArray(57,iP) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     AuxArray(86,iP) = Yqc*AuxArray(57,iP) + alphaYpq*TmpArray7(57,2)
     AuxArray(87,iP) = Zqc*AuxArray(57,iP) + alphaZpq*TmpArray7(57,2)
     AuxArray(88,iP) = Xqc*AuxArray(60,iP) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     AuxArray(89,iP) = Yqc*AuxArray(59,iP) + alphaYpq*TmpArray7(59,2)
     AuxArray(90,iP) = Xqc*AuxArray(62,iP) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     AuxArray(91,iP) = Xqc*AuxArray(63,iP) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     AuxArray(92,iP) = Zqc*AuxArray(60,iP) + alphaZpq*TmpArray7(60,2)
     AuxArray(93,iP) = Yqc*AuxArray(62,iP) + alphaYpq*TmpArray7(62,2)
     AuxArray(94,iP) = Xqc*AuxArray(66,iP) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     AuxArray(95,iP) = Xqc*AuxArray(67,iP) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     AuxArray(96,iP) = Zqc*AuxArray(63,iP) + alphaZpq*TmpArray7(63,2)
     AuxArray(97,iP) = Xqc*AuxArray(69,iP) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     AuxArray(98,iP) = Yqc*AuxArray(66,iP) + alphaYpq*TmpArray7(66,2)
     AuxArray(99,iP) = Xqc*AuxArray(71,iP) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     AuxArray(100,iP) = Xqc*AuxArray(72,iP) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     AuxArray(101,iP) = Zqc*AuxArray(67,iP) + alphaZpq*TmpArray7(67,2)
     AuxArray(102,iP) = Xqc*AuxArray(74,iP) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     AuxArray(103,iP) = Xqc*AuxArray(75,iP) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     AuxArray(104,iP) = Yqc*AuxArray(71,iP) + alphaYpq*TmpArray7(71,2)
     AuxArray(105,iP) = Xqc*AuxArray(77,iP) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     AuxArray(106,iP) = Xqc*AuxArray(78,iP) + alphaXpq*TmpArray7(78,2)
     AuxArray(107,iP) = Xqc*AuxArray(79,iP) + alphaXpq*TmpArray7(79,2)
     AuxArray(108,iP) = Xqc*AuxArray(80,iP) + alphaXpq*TmpArray7(80,2)
     AuxArray(109,iP) = Xqc*AuxArray(81,iP) + alphaXpq*TmpArray7(81,2)
     AuxArray(110,iP) = Xqc*AuxArray(82,iP) + alphaXpq*TmpArray7(82,2)
     AuxArray(111,iP) = Xqc*AuxArray(83,iP) + alphaXpq*TmpArray7(83,2)
     AuxArray(112,iP) = Xqc*AuxArray(84,iP) + alphaXpq*TmpArray7(84,2)
     AuxArray(113,iP) = Yqc*AuxArray(78,iP) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     AuxArray(114,iP) = Zqc*AuxArray(78,iP) + alphaZpq*TmpArray7(78,2)
     AuxArray(115,iP) = Yqc*AuxArray(80,iP) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     AuxArray(116,iP) = Yqc*AuxArray(81,iP) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     AuxArray(117,iP) = Yqc*AuxArray(82,iP) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     AuxArray(118,iP) = Yqc*AuxArray(83,iP) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     AuxArray(119,iP) = Yqc*AuxArray(84,iP) + alphaYpq*TmpArray7(84,2)
     AuxArray(120,iP) = Zqc*AuxArray(84,iP) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
 end subroutine
end module
