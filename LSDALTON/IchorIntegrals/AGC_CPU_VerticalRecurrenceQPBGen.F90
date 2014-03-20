MODULE AGC_CPU_OBS_VERTICALRECURRENCEMODBGen
 use IchorPrecisionModule
  
 CONTAINS


subroutine VerticalRecurrenceCPUGen1B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: TABFJW(0:4,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(4,nPrimQ*nPrimP*nPassP)
  !local variables
  Integer :: iP,iPrimQ,iPrimP,iPassP,ipnt,iAtomA,iAtomB
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,alphaP,RJ000(0:1)
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
  !ThetaAux(n,1,0,0) = Xpb*ThetaAux(n,0,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,0,0,0)
  !i = 0 last 2 term vanish
  !We include scaling of RJ000 
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$OMP         RWVAL,GVAL,&
!$OMP         mPx,mPy,mPz,&
!$OMP         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,&
!$OMP         invexpP,&
!$OMP         PREF,&
!$OMP         TMP1,TMP2,&
!$OMP         iP,iPrimQ,iPrimP,iPassP)
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
     AUXarray(1,iP) = TMP1
     AUXarray(2,iP) = Xpb*TMP1 + alphaXpq*TMP2
     AUXarray(3,iP) = Ypb*TMP1 + alphaYpq*TMP2
     AUXarray(4,iP) = Zpb*TMP1 + alphaZpq*TMP2
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
!$OMP END DO
end subroutine VerticalRecurrenceCPUGen1B

subroutine VerticalRecurrenceCPUGen2B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 2,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(   10,nPrimQ*nPrimP*nPassP)
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
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         mPx,mPy,mPz,&
!$OMP         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,&
!$OMP         invexpP,inv2expP,&
!$OMP         PREF,&
!$OMP         TMParray1,&
!$OMP         TMParray2,&
!$OMP         TwoTerms,&
!$OMP         iP,iPrimQ,iPrimP,iPassP)
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
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xpb*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Ypb*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zpb*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,iP) = Xpb*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xpb*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xpb*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Ypb*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Ypb*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zpb*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPUGen3B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 3,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(   20,nPrimQ*nPrimP*nPassP)
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
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         mPx,mPy,mPz,&
!$OMP         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,&
!$OMP         invexpP,inv2expP,&
!$OMP         PREF,&
!$OMP         TMParray1,&
!$OMP         TMParray2,&
!$OMP         TMParray3,&
!$OMP         TwoTerms,&
!$OMP         iP,iPrimQ,iPrimP,iPassP)
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
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xpb*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Ypb*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zpb*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpb*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypb*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpb*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,iP) = Xpb*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xpb*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xpb*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Ypb*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Ypb*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zpb*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,2) + alphaP*TmpArray1(1,3))
     tmpArray3(5,2) = Xpb*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xpb*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xpb*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Ypb*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Ypb*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zpb*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(AuxArray(2,IP) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(3,IP) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(4,IP) + alphaP*TmpArray2(4,2))
     AuxArray(11,iP) = Xpb*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Ypb*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zpb*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xpb*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xpb*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xpb*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Ypb*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zpb*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Ypb*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zpb*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPUGen4B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 4,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(   35,nPrimQ*nPrimP*nPassP)
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
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         mPx,mPy,mPz,&
!$OMP         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,&
!$OMP         invexpP,inv2expP,&
!$OMP         PREF,&
!$OMP         TMParray1,&
!$OMP         TMParray2,&
!$OMP         TMParray3,&
!$OMP         TMParray4,&
!$OMP         TwoTerms,&
!$OMP         iP,iPrimQ,iPrimP,iPassP)
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
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xpb*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Ypb*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zpb*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpb*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypb*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpb*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xpb*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Ypb*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zpb*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,iP) = Xpb*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xpb*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xpb*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Ypb*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Ypb*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zpb*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(AuxArray(2,IP) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(3,IP) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(4,IP) + alphaP*TmpArray2(4,2))
     AuxArray(11,iP) = Xpb*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Ypb*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zpb*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xpb*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xpb*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xpb*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Ypb*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zpb*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Ypb*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zpb*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(5,IP) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(AuxArray(8,IP) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(AuxArray(10,IP) + alphaP*TmpArray3(10,2))
     AuxArray(21,iP) = Xpb*AuxArray(11,iP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,iP) = Ypb*AuxArray(11,iP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,iP) = Zpb*AuxArray(11,iP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,iP) = Xpb*AuxArray(14,iP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,iP) = Ypb*AuxArray(13,iP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,iP) = Xpb*AuxArray(16,iP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,iP) = Xpb*AuxArray(17,iP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,iP) = Xpb*AuxArray(18,iP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,iP) = Xpb*AuxArray(19,iP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,iP) = Xpb*AuxArray(20,iP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,iP) = Ypb*AuxArray(17,iP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,iP) = Zpb*AuxArray(17,iP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,iP) = Ypb*AuxArray(19,iP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,iP) = Ypb*AuxArray(20,iP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,iP) = Zpb*AuxArray(20,iP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPUGen5B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 5,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(   56,nPrimQ*nPrimP*nPassP)
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
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         mPx,mPy,mPz,&
!$OMP         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,&
!$OMP         invexpP,inv2expP,&
!$OMP         PREF,&
!$OMP         TMParray1,&
!$OMP         TMParray2,&
!$OMP         TMParray3,&
!$OMP         TMParray4,&
!$OMP         TMParray5,&
!$OMP         TwoTerms,&
!$OMP         iP,iPrimQ,iPrimP,iPassP)
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
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 6) = PREF*RJ000Array( 5,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xpb*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Ypb*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zpb*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,iP) = Xpb*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xpb*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xpb*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Ypb*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Ypb*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zpb*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(AuxArray(2,IP) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(3,IP) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(4,IP) + alphaP*TmpArray2(4,2))
     AuxArray(11,iP) = Xpb*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Ypb*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zpb*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xpb*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xpb*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xpb*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Ypb*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zpb*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Ypb*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zpb*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(5,IP) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(AuxArray(8,IP) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(AuxArray(10,IP) + alphaP*TmpArray3(10,2))
     AuxArray(21,iP) = Xpb*AuxArray(11,iP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,iP) = Ypb*AuxArray(11,iP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,iP) = Zpb*AuxArray(11,iP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,iP) = Xpb*AuxArray(14,iP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,iP) = Ypb*AuxArray(13,iP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,iP) = Xpb*AuxArray(16,iP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,iP) = Xpb*AuxArray(17,iP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,iP) = Xpb*AuxArray(18,iP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,iP) = Xpb*AuxArray(19,iP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,iP) = Xpb*AuxArray(20,iP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,iP) = Ypb*AuxArray(17,iP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,iP) = Zpb*AuxArray(17,iP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,iP) = Ypb*AuxArray(19,iP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,iP) = Ypb*AuxArray(20,iP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,iP) = Zpb*AuxArray(20,iP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(11,IP) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(AuxArray(14,IP) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(AuxArray(16,IP) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(AuxArray(17,IP) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(AuxArray(19,IP) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(AuxArray(20,IP) + alphaP*TmpArray4(20,2))
     AuxArray(36,iP) = Xpb*AuxArray(21,iP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,iP) = Ypb*AuxArray(21,iP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,iP) = Zpb*AuxArray(21,iP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,iP) = Xpb*AuxArray(24,iP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,iP) = Ypb*AuxArray(23,iP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,iP) = Xpb*AuxArray(26,iP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,iP) = Xpb*AuxArray(27,iP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,iP) = Zpb*AuxArray(24,iP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,iP) = Ypb*AuxArray(26,iP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,iP) = Xpb*AuxArray(30,iP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,iP) = Xpb*AuxArray(31,iP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,iP) = Xpb*AuxArray(32,iP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,iP) = Xpb*AuxArray(33,iP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,iP) = Xpb*AuxArray(34,iP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,iP) = Xpb*AuxArray(35,iP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,iP) = Ypb*AuxArray(31,iP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,iP) = Zpb*AuxArray(31,iP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,iP) = Ypb*AuxArray(33,iP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,iP) = Ypb*AuxArray(34,iP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,iP) = Ypb*AuxArray(35,iP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,iP) = Zpb*AuxArray(35,iP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPUGen6B(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 6,nPrimQ,nPrimP,nPassP)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Bcenter(3,nAtomsB)
  real(realk),intent(inout) :: AUXarray(   84,nPrimQ*nPrimP*nPassP)
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
!$OMP DO &
!$OMP PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         mPx,mPy,mPz,&
!$OMP         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,&
!$OMP         invexpP,inv2expP,&
!$OMP         PREF,&
!$OMP         TMParray1,&
!$OMP         TMParray2,&
!$OMP         TMParray3,&
!$OMP         TMParray4,&
!$OMP         TMParray5,&
!$OMP         TMParray6,&
!$OMP         TwoTerms,&
!$OMP         iP,iPrimQ,iPrimP,iPassP)
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
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 6) = PREF*RJ000Array( 5,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 7) = PREF*RJ000Array( 6,iPrimQ,iPrimP,iPassP)
     AuxArray(2,iP) = Xpb*AuxArray(1,iP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,iP) = Ypb*AuxArray(1,iP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,iP) = Zpb*AuxArray(1,iP) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,iP) = Xpb*AuxArray(2,iP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,iP) = Xpb*AuxArray(3,iP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,iP) = Xpb*AuxArray(4,iP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,iP) = Ypb*AuxArray(3,iP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,iP) = Ypb*AuxArray(4,iP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,iP) = Zpb*AuxArray(4,iP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(AuxArray(2,IP) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(3,IP) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(4,IP) + alphaP*TmpArray2(4,2))
     AuxArray(11,iP) = Xpb*AuxArray(5,iP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,iP) = Ypb*AuxArray(5,iP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,iP) = Zpb*AuxArray(5,iP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,iP) = Xpb*AuxArray(8,iP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,iP) = Xpb*AuxArray(9,iP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,iP) = Xpb*AuxArray(10,iP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,iP) = Ypb*AuxArray(8,iP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,iP) = Zpb*AuxArray(8,iP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,iP) = Ypb*AuxArray(10,iP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,iP) = Zpb*AuxArray(10,iP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(5,IP) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(AuxArray(8,IP) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(AuxArray(10,IP) + alphaP*TmpArray3(10,2))
     AuxArray(21,iP) = Xpb*AuxArray(11,iP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,iP) = Ypb*AuxArray(11,iP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,iP) = Zpb*AuxArray(11,iP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,iP) = Xpb*AuxArray(14,iP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,iP) = Ypb*AuxArray(13,iP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,iP) = Xpb*AuxArray(16,iP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,iP) = Xpb*AuxArray(17,iP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,iP) = Xpb*AuxArray(18,iP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,iP) = Xpb*AuxArray(19,iP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,iP) = Xpb*AuxArray(20,iP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,iP) = Ypb*AuxArray(17,iP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,iP) = Zpb*AuxArray(17,iP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,iP) = Ypb*AuxArray(19,iP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,iP) = Ypb*AuxArray(20,iP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,iP) = Zpb*AuxArray(20,iP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(11,IP) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(AuxArray(14,IP) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(AuxArray(16,IP) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(AuxArray(17,IP) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(AuxArray(19,IP) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(AuxArray(20,IP) + alphaP*TmpArray4(20,2))
     AuxArray(36,iP) = Xpb*AuxArray(21,iP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,iP) = Ypb*AuxArray(21,iP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,iP) = Zpb*AuxArray(21,iP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,iP) = Xpb*AuxArray(24,iP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,iP) = Ypb*AuxArray(23,iP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,iP) = Xpb*AuxArray(26,iP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,iP) = Xpb*AuxArray(27,iP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,iP) = Zpb*AuxArray(24,iP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,iP) = Ypb*AuxArray(26,iP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,iP) = Xpb*AuxArray(30,iP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,iP) = Xpb*AuxArray(31,iP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,iP) = Xpb*AuxArray(32,iP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,iP) = Xpb*AuxArray(33,iP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,iP) = Xpb*AuxArray(34,iP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,iP) = Xpb*AuxArray(35,iP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,iP) = Ypb*AuxArray(31,iP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,iP) = Zpb*AuxArray(31,iP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,iP) = Ypb*AuxArray(33,iP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,iP) = Ypb*AuxArray(34,iP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,iP) = Ypb*AuxArray(35,iP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,iP) = Zpb*AuxArray(35,iP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     TwoTerms(1) = inv2expP*(AuxArray(21,IP) + alphaP*TmpArray5(21,2))
     TwoTerms(2) = inv2expP*(AuxArray(24,IP) + alphaP*TmpArray5(24,2))
     TwoTerms(3) = inv2expP*(AuxArray(26,IP) + alphaP*TmpArray5(26,2))
     TwoTerms(4) = inv2expP*(AuxArray(27,IP) + alphaP*TmpArray5(27,2))
     TwoTerms(5) = inv2expP*(AuxArray(30,IP) + alphaP*TmpArray5(30,2))
     TwoTerms(6) = inv2expP*(AuxArray(31,IP) + alphaP*TmpArray5(31,2))
     TwoTerms(7) = inv2expP*(AuxArray(33,IP) + alphaP*TmpArray5(33,2))
     TwoTerms(8) = inv2expP*(AuxArray(34,IP) + alphaP*TmpArray5(34,2))
     TwoTerms(9) = inv2expP*(AuxArray(35,IP) + alphaP*TmpArray5(35,2))
     AuxArray(57,iP) = Xpb*AuxArray(36,iP) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,iP) = Ypb*AuxArray(36,iP) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,iP) = Zpb*AuxArray(36,iP) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,iP) = Xpb*AuxArray(39,iP) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,iP) = Ypb*AuxArray(38,iP) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,iP) = Xpb*AuxArray(41,iP) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,iP) = Xpb*AuxArray(42,iP) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,iP) = Zpb*AuxArray(39,iP) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,iP) = Ypb*AuxArray(41,iP) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,iP) = Xpb*AuxArray(45,iP) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,iP) = Xpb*AuxArray(46,iP) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,iP) = Zpb*AuxArray(42,iP) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,iP) = Xpb*AuxArray(48,iP) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,iP) = Ypb*AuxArray(45,iP) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,iP) = Xpb*AuxArray(50,iP) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,iP) = Xpb*AuxArray(51,iP) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,iP) = Xpb*AuxArray(52,iP) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,iP) = Xpb*AuxArray(53,iP) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,iP) = Xpb*AuxArray(54,iP) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,iP) = Xpb*AuxArray(55,iP) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,iP) = Xpb*AuxArray(56,iP) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,iP) = Ypb*AuxArray(51,iP) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,iP) = Zpb*AuxArray(51,iP) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,iP) = Ypb*AuxArray(53,iP) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,iP) = Ypb*AuxArray(54,iP) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,iP) = Ypb*AuxArray(55,iP) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,iP) = Ypb*AuxArray(56,iP) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,iP) = Zpb*AuxArray(56,iP) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP
!$OMP END DO
 end subroutine
end module
