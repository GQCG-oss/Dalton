module SPAGC_GPU_OBS_VERTICALRECURRENCEMODDSeg
 use IchorPrecisionMod
  
 CONTAINS


subroutine SPVerticalRecurrenceGPUSeg1D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: TABFJW(0:4,0:1200)
  real(reals),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(reals),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(reals),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(reals),intent(in) :: Dcenter(3)
  real(reals),intent(inout) :: AUXarray(nPassP,4)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  Integer :: iP,iPrimQ,iPrimP,iPassP,ipnt,iAtomA,iAtomB
  Integer :: iPrimQP
  real(reals) :: Xqd,Yqd,Zqd
  real(reals) :: mPX,mPY,mPZ,invexpQ,alphaQ,RJ000(0:1)
  real(reals) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(reals) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(reals),PARAMETER :: TENTH = 0.01E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D2=2.0E0_reals
  real(reals),PARAMETER :: D2JP36=  3.8000000000000000E+01_reals
  real(reals),parameter :: D1=1.0E0_reals,D03333=1.0E0_reals/3.0E0_reals
  real(reals),PARAMETER :: D4 = 4E0_reals, D100=100E0_reals
  real(reals),PARAMETER :: COEF3 = - D1/6E0_reals, COEF4 = D1/24E0_reals
  real(reals),PARAMETER :: SMALL = 1E-15_reals,D12 = 12.0E0_reals
  real(reals), PARAMETER :: GFAC0 =  D2*0.4999489092E0_reals
  real(reals), PARAMETER :: GFAC1 = -D2*0.2473631686E0_reals
  real(reals), PARAMETER :: GFAC2 =  D2*0.321180909E0_reals
  real(reals), PARAMETER :: GFAC3 = -D2*0.3811559346E0_reals
  real(reals), parameter :: PI=3.14159265358979323846E0_reals
  real(reals), PARAMETER :: SQRTPI = 1.77245385090551602730E00_reals
  real(reals), PARAMETER :: SQRPIH = SQRTPI/D2
  real(reals), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
  !ThetaAux(n,1,0,0) = Xqd*ThetaAux(n,0,0,0) + (-alpha/q*Xpq)*ThetaAux(n+1,0,0,0)
  !i = 0 last 2 term vanish
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP PRIVATE(iPassP) PRESENT(AUXarray)
  DO iPassP = 1,nPassP
   AUXarray(iPassP,1)=0.0E0_reals
   AUXarray(iPassP,2)=0.0E0_reals
   AUXarray(iPassP,3)=0.0E0_reals
   AUXarray(iPassP,4)=0.0E0_reals
  ENDDO
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
!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,TABFJW,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Dcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPassP
   DO iPrimQP=1,nPrimQ*nPrimP
    iPrimQ = iPrimQP - ((iPrimQP-1)/nPrimQ)*nPrimQ
    iPrimP = (iPrimQP-1)/nPrimQ + 1
    iPassP = iP
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
     AUXarray(iP,1) = AUXarray(iP,1) + TMP1
     AUXarray(iP,2) = AUXarray(iP,2) + Xqd*TMP1 + alphaXpq*TMP2
     AUXarray(iP,3) = AUXarray(iP,3) + Yqd*TMP1 + alphaYpq*TMP2
     AUXarray(iP,4) = AUXarray(iP,4) + Zqd*TMP1 + alphaZpq*TMP2
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPassP
end subroutine SPVerticalRecurrenceGPUSeg1D

subroutine SPVerticalRecurrenceGPUSeg2D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 2)
  real(reals),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(reals),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(reals),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(reals),intent(in) :: Dcenter(3)
  real(reals),intent(inout) :: AUXarray(nPassP,   10)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  Integer :: iPrimQP
  real(reals) :: TMPAUXarray(    4)
  real(reals) :: Xqd,Yqd,Zqd
  real(reals) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ
  real(reals) :: TwoTerms(   1)
  real(reals) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(reals),parameter :: D2=2.0E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D1=1.0E0_reals
  real(reals) :: TMParray1(  1:  1,2:3)
  real(reals) :: TMParray2(  2:  4,2:2)
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP PRIVATE(iPassP,iTUV) PRESENT(AUXarray)
  DO iPassP = 1,nPassP
   DO iTUV=1,   10
    AUXarray(iPassP,iTUV)=0.0E0_reals
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqd,Yqd,Zqd,alphaQ,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpQ,inv2expQ,&
!$ACC         PREF,&
!$ACC         TMPAUXarray,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Dcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPassP
   DO iPrimQP=1,nPrimQ*nPrimP
    iPrimQ = iPrimQP - ((iPrimQP-1)/nPrimQ)*nPrimQ
    iPrimP = (iPrimQP-1)/nPrimQ + 1
    iPassP = iP
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
     TMPAuxarray(1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     TMPAuxArray(2) = Xqd*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqd*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqd*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(1) + alphaQ*TmpArray1(1,2))
!$ACC LOOP SEQ
     do iTUV = 1,    4
      AuxArray(iP,iTUV) = AuxArray(iP,iTUV) + TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,5) = AuxArray(iP,5) + Xqd*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = AuxArray(iP,6) + Xqd*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = AuxArray(iP,7) + Xqd*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = AuxArray(iP,8) + Yqd*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = AuxArray(iP,9) + Yqd*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = AuxArray(iP,10) + Zqd*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPassP
 end subroutine

subroutine SPVerticalRecurrenceGPUSeg3D(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 3)
  real(reals),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(reals),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(reals),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(reals),intent(in) :: Dcenter(3)
  real(reals),intent(inout) :: AUXarray(nPassP,   20)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  Integer :: iPrimQP
  real(reals) :: TMPAUXarray(   10)
  real(reals) :: Xqd,Yqd,Zqd
  real(reals) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ
  real(reals) :: TwoTerms(   3)
  real(reals) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(reals),parameter :: D2=2.0E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D1=1.0E0_reals
  real(reals) :: TMParray1(  1:  1,2:4)
  real(reals) :: TMParray2(  2:  4,2:3)
  real(reals) :: TMParray3(  5: 10,2:2)
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP PRIVATE(iPassP,iTUV) PRESENT(AUXarray)
  DO iPassP = 1,nPassP
   DO iTUV=1,   20
    AUXarray(iPassP,iTUV)=0.0E0_reals
   ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Xqd,Yqd,Zqd,alphaQ,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpQ,inv2expQ,&
!$ACC         PREF,&
!$ACC         TMPAUXarray,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPrimQP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Dcenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPassP
   DO iPrimQP=1,nPrimQ*nPrimP
    iPrimQ = iPrimQP - ((iPrimQP-1)/nPrimQ)*nPrimQ
    iPrimP = (iPrimQP-1)/nPrimQ + 1
    iPassP = iP
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
     TMPAuxarray(1) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 3)
     TMPAuxArray(2) = Xqd*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqd*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqd*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqd*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqd*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqd*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(1) + alphaQ*TmpArray1(1,2))
     TMPAuxArray(5) = Xqd*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     TMPAuxArray(6) = Xqd*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     TMPAuxArray(7) = Xqd*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     TMPAuxArray(8) = Yqd*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     TMPAuxArray(9) = Yqd*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     TMPAuxArray(10) = Zqd*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,2) + alphaQ*TmpArray1(1,3))
     tmpArray3(5,2) = Xqd*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xqd*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xqd*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Yqd*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Yqd*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zqd*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(4) + alphaQ*TmpArray2(4,2))
!$ACC LOOP SEQ
     do iTUV = 1,   10
      AuxArray(iP,iTUV) = AuxArray(iP,iTUV) + TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,11) = AuxArray(iP,11) + Xqd*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(iP,12) = AuxArray(iP,12) + Yqd*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     AuxArray(iP,13) = AuxArray(iP,13) + Zqd*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     AuxArray(iP,14) = AuxArray(iP,14) + Xqd*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     AuxArray(iP,15) = AuxArray(iP,15) + Xqd*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     AuxArray(iP,16) = AuxArray(iP,16) + Xqd*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     AuxArray(iP,17) = AuxArray(iP,17) + Yqd*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(iP,18) = AuxArray(iP,18) + Zqd*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     AuxArray(iP,19) = AuxArray(iP,19) + Yqd*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     AuxArray(iP,20) = AuxArray(iP,20) + Zqd*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
   ENDDO !iPrimQP = 1,nPrimQ*nPrimP
  ENDDO !iP = 1,nPassP
 end subroutine
end module
