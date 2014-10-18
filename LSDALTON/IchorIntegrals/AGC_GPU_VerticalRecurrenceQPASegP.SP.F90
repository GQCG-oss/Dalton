module SPAGC_GPU_OBS_VERTICALRECURRENCEMODASegP
 use IchorPrecisionMod
  
 CONTAINS

subroutine SPVerticalRecurrenceGPUSegP0(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,&
         & Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
         & AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: TABFJW(0:3,0:1200)
  real(reals),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(reals),intent(in) :: integralPrefactor(nprimQ,nPrimP)
  real(reals),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(reals),intent(in) :: QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(reals),intent(inout) :: AUXarray(nPrimQ*nPassP)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  real(reals),PARAMETER :: D2JP36=  3.6000000000000000E+01_reals
  real(reals),parameter :: D2=2.0E0_reals
  real(reals),PARAMETER :: D05 =0.5E0_reals,D1=1E0_reals
  real(reals),PARAMETER :: D4 = 4E0_reals, D100=100E0_reals
  real(reals),parameter :: D12 = 12E0_reals, TENTH = 0.01E0_reals
  real(reals),PARAMETER :: COEF3 = - D1/6E0_reals, COEF4 = D1/24E0_reals
  real(reals),PARAMETER :: COEF5 = - D1/120E0_reals, COEF6 = D1/720E0_reals
  real(reals),PARAMETER :: GFAC0 =  D2*0.4999489092E0_reals
  real(reals),PARAMETER :: GFAC1 = -D2*0.2473631686E0_reals
  real(reals),PARAMETER :: GFAC2 =  D2*0.321180909E0_reals
  real(reals),PARAMETER :: GFAC3 = -D2*0.3811559346E0_reals
  real(reals),parameter :: PI=3.14159265358979323846E0_reals
  real(reals),PARAMETER :: SQRTPI = 1.77245385090551602730E00_reals
  real(reals),PARAMETER :: SQRPIH = SQRTPI/D2
  real(reals),PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
!  real(reals),PARAMETER :: SMALL = 1E-15_reals
  real(reals) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA,WVAL
  real(reals) :: W2,W3,PX,PY,PZ,XPQ,YPQ,ZPQ,squaredDistance,RJ000
  Integer :: IPNT,iAtomA,iAtomB
  Integer :: iP,iPrimQ,iPrimP,iPassP
!$ACC PARALLEL LOOP PRIVATE(iP) PRESENT(AUXarray)
  DO iP = 1,nPrimQ*nPassP
    AUXarray(iP)=0.0E0_reals
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         RWVAL,GVAL,&
!$ACC         Px,Py,Pz,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,TABFJW,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPassP
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
    px = Pcent(1,iPrimP,iAtomA,iAtomB)
    py = Pcent(2,iPrimP,iAtomA,iAtomB)
    pz = Pcent(3,iPrimP,iAtomA,iAtomB)
     Xpq = px - Qcent(1,iPrimQ)
     Ypq = py - Qcent(2,iPrimQ)
     Zpq = pz - Qcent(3,iPrimQ)
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance
     !  0 < WVAL < 0.000001
!     IF (ABS(WVAL) .LT. SMALL) THEN
!      RJ000 = D1
!     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000 = TABFJW(0,IPNT)-TABFJW(1,IPNT)*WDIFF+TABFJW(2,IPNT)*W2+TABFJW(3,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000 = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
     !  (2J+36) < WVAL 
     ELSE
      RJ000 = SQRT(PID4/WVAL)
     ENDIF
     AUXarray(iP)=AUXarray(iP) + integralPrefactor(iPrimQ,iPrimP)*&
          & QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)*RJ000
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPassP
end subroutine SPVerticalRecurrenceGPUSegP0

subroutine SPVerticalRecurrenceGPUSegP1A(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: TABFJW(0:4,0:1200)
  real(reals),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(reals),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(reals),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(reals),intent(in) :: Acenter(3,nAtomsA)
  real(reals),intent(inout) :: AUXarray(nPrimQ*nPassP,4)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  Integer :: iP,iPrimQ,iPrimP,iPassP,ipnt,iAtomA,iAtomB
  real(reals) :: Ax,Ay,Az,Xpa,Ypa,Zpa
  real(reals) :: mPX,mPY,mPZ,invexpP,alphaP,RJ000(0:1)
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
  !ThetaAux(n,1,0,0) = Xpa*ThetaAux(n,0,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,0,0,0)
  !i = 0 last 2 term vanish
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP PRIVATE(iP) PRESENT(AUXarray)
  DO iP = 1,nPrimQ*nPassP
    AUXarray(iP,1)=0.0E0_reals
    AUXarray(iP,2)=0.0E0_reals
    AUXarray(iP,3)=0.0E0_reals
    AUXarray(iP,4)=0.0E0_reals
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         RWVAL,GVAL,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpP,&
!$ACC         PREF,&
!$ACC         TMP1,TMP2,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,TABFJW,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Pexp,Acenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPassP
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
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
     AUXarray(iP,2) = AUXarray(iP,2) + Xpa*TMP1 + alphaXpq*TMP2
     AUXarray(iP,3) = AUXarray(iP,3) + Ypa*TMP1 + alphaYpq*TMP2
     AUXarray(iP,4) = AUXarray(iP,4) + Zpa*TMP1 + alphaZpq*TMP2
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPassP
end subroutine SPVerticalRecurrenceGPUSegP1A

subroutine SPVerticalRecurrenceGPUSegP2A(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 2)
  real(reals),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(reals),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(reals),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(reals),intent(in) :: Acenter(3,nAtomsA)
  real(reals),intent(inout) :: AUXarray(nPrimQ*nPassP,   10)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(reals) :: TMPAUXarray(    4)
  real(reals) :: Ax,Ay,Az,Xpa,Ypa,Zpa
  real(reals) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(reals) :: TwoTerms(   1)
  real(reals) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(reals),parameter :: D2=2.0E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D1=1.0E0_reals
  real(reals) :: TMParray1(  1:  1,2:3)
  real(reals) :: TMParray2(  2:  4,2:2)
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP PRIVATE(iP,iTUV) PRESENT(AUXarray)
  DO iP = 1,nPrimQ*nPassP
    DO iTUV=1,   10
     AUXarray(iP,iTUV)=0.0E0_reals
    ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,&
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
!$ACC        Pexp,Acenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPassP
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
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
     TMPAuxArray(2) = Xpa*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypa*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpa*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpa*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypa*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpa*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(1) + alphaP*TmpArray1(1,2))
!$ACC LOOP SEQ
     do iTUV = 1,    4
      AuxArray(iP,iTUV) = AuxArray(iP,iTUV) + TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,5) = AuxArray(iP,5) + Xpa*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = AuxArray(iP,6) + Xpa*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = AuxArray(iP,7) + Xpa*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = AuxArray(iP,8) + Ypa*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = AuxArray(iP,9) + Ypa*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = AuxArray(iP,10) + Zpa*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPassP
 end subroutine

subroutine SPVerticalRecurrenceGPUSegP3A(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 3)
  real(reals),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(reals),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(reals),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(reals),intent(in) :: Acenter(3,nAtomsA)
  real(reals),intent(inout) :: AUXarray(nPrimQ*nPassP,   20)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(reals) :: TMPAUXarray(   10)
  real(reals) :: Ax,Ay,Az,Xpa,Ypa,Zpa
  real(reals) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(reals) :: TwoTerms(   3)
  real(reals) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(reals),parameter :: D2=2.0E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D1=1.0E0_reals
  real(reals) :: TMParray1(  1:  1,2:4)
  real(reals) :: TMParray2(  2:  4,2:3)
  real(reals) :: TMParray3(  5: 10,2:2)
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP PRIVATE(iP,iTUV) PRESENT(AUXarray)
  DO iP = 1,nPrimQ*nPassP
    DO iTUV=1,   20
     AUXarray(iP,iTUV)=0.0E0_reals
    ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,&
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
!$ACC        Pexp,Acenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPassP
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
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
     TMPAuxArray(2) = Xpa*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypa*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpa*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpa*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypa*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpa*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpa*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypa*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpa*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expP*(TMPAuxArray(1) + alphaP*TmpArray1(1,2))
     TMPAuxArray(5) = Xpa*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     TMPAuxArray(6) = Xpa*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     TMPAuxArray(7) = Xpa*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     TMPAuxArray(8) = Ypa*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     TMPAuxArray(9) = Ypa*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     TMPAuxArray(10) = Zpa*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,2) + alphaP*TmpArray1(1,3))
     tmpArray3(5,2) = Xpa*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xpa*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xpa*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Ypa*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Ypa*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zpa*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TMPAuxArray(2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(4) + alphaP*TmpArray2(4,2))
!$ACC LOOP SEQ
     do iTUV = 1,   10
      AuxArray(iP,iTUV) = AuxArray(iP,iTUV) + TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,11) = AuxArray(iP,11) + Xpa*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(iP,12) = AuxArray(iP,12) + Ypa*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     AuxArray(iP,13) = AuxArray(iP,13) + Zpa*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     AuxArray(iP,14) = AuxArray(iP,14) + Xpa*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     AuxArray(iP,15) = AuxArray(iP,15) + Xpa*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     AuxArray(iP,16) = AuxArray(iP,16) + Xpa*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     AuxArray(iP,17) = AuxArray(iP,17) + Ypa*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(iP,18) = AuxArray(iP,18) + Zpa*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     AuxArray(iP,19) = AuxArray(iP,19) + Ypa*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     AuxArray(iP,20) = AuxArray(iP,20) + Zpa*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPassP
 end subroutine

subroutine SPVerticalRecurrenceGPUSegP4A(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: RJ000Array(nPrimQ,nPrimP,nPassP,0: 4)
  real(reals),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(reals),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(reals),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(reals),intent(in) :: Acenter(3,nAtomsA)
  real(reals),intent(inout) :: AUXarray(nPrimQ*nPassP,   35)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(reals) :: TMPAUXarray(   20)
  real(reals) :: Ax,Ay,Az,Xpa,Ypa,Zpa
  real(reals) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(reals) :: TwoTerms(   6)
  real(reals) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(reals),parameter :: D2=2.0E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D1=1.0E0_reals
  real(reals) :: TMParray1(  1:  1,2:5)
  real(reals) :: TMParray2(  2:  4,2:4)
  real(reals) :: TMParray3(  5: 10,2:3)
  real(reals) :: TMParray4( 11: 20,2:2)
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$ACC PARALLEL LOOP PRIVATE(iP,iTUV) PRESENT(AUXarray)
  DO iP = 1,nPrimQ*nPassP
    DO iTUV=1,   35
     AUXarray(iP,iTUV)=0.0E0_reals
    ENDDO
  ENDDO
!$ACC PARALLEL LOOP &
!$ACC PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         mPx,mPy,mPz,&
!$ACC         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,&
!$ACC         alphaXpq,alphaYpq,alphaZpq,&
!$ACC         invexpP,inv2expP,&
!$ACC         PREF,&
!$ACC         TMPAUXarray,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TMParray4,&
!$ACC         TwoTerms,&
!$ACC         iP,iPrimQ,iPrimP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Pexp,Acenter, &
!$ACC        nPrimP,nPrimQ,nPassP) ASYNC(iASync)
  DO iP = 1,nPrimQ*nPassP
   DO iPrimP=1, nPrimP
    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ
    iPassP = (iP-1)/nPrimQ + 1
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
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
     TMParray1(1, 5) = PREF*RJ000Array(iPrimQ,iPrimP,iPassP, 4)
     TMPAuxArray(2) = Xpa*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypa*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpa*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpa*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypa*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpa*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpa*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypa*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpa*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xpa*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Ypa*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zpa*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     TwoTerms(1) = inv2expP*(TMPAuxArray(1) + alphaP*TmpArray1(1,2))
     TMPAuxArray(5) = Xpa*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     TMPAuxArray(6) = Xpa*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     TMPAuxArray(7) = Xpa*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     TMPAuxArray(8) = Ypa*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     TMPAuxArray(9) = Ypa*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     TMPAuxArray(10) = Zpa*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,2) + alphaP*TmpArray1(1,3))
     tmpArray3(5,2) = Xpa*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xpa*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xpa*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Ypa*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Ypa*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zpa*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,3) + alphaP*TmpArray1(1,4))
     tmpArray3(5,3) = Xpa*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xpa*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xpa*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Ypa*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Ypa*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zpa*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TMPAuxArray(2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(4) + alphaP*TmpArray2(4,2))
     TMPAuxArray(11) = Xpa*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Ypa*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zpa*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xpa*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xpa*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xpa*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Ypa*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zpa*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Ypa*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zpa*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,2) + alphaP*TmpArray2(2,3))
     TwoTerms(2) = inv2expP*(TmpArray2(3,2) + alphaP*TmpArray2(3,3))
     TwoTerms(3) = inv2expP*(TmpArray2(4,2) + alphaP*TmpArray2(4,3))
     tmpArray4(11,2) = Xpa*tmpArray3(5,2) + alphaXpq*TmpArray3(5,3) + 2*TwoTerms(1)
     tmpArray4(12,2) = Ypa*tmpArray3(5,2) + alphaYpq*TmpArray3(5,3)
     tmpArray4(13,2) = Zpa*tmpArray3(5,2) + alphaZpq*TmpArray3(5,3)
     tmpArray4(14,2) = Xpa*tmpArray3(8,2) + alphaXpq*TmpArray3(8,3)
     tmpArray4(15,2) = Xpa*tmpArray3(9,2) + alphaXpq*TmpArray3(9,3)
     tmpArray4(16,2) = Xpa*tmpArray3(10,2) + alphaXpq*TmpArray3(10,3)
     tmpArray4(17,2) = Ypa*tmpArray3(8,2) + alphaYpq*TmpArray3(8,3) + 2*TwoTerms(2)
     tmpArray4(18,2) = Zpa*tmpArray3(8,2) + alphaZpq*TmpArray3(8,3)
     tmpArray4(19,2) = Ypa*tmpArray3(10,2) + alphaYpq*TmpArray3(10,3)
     tmpArray4(20,2) = Zpa*tmpArray3(10,2) + alphaZpq*TmpArray3(10,3) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(10) + alphaP*TmpArray3(10,2))
!$ACC LOOP SEQ
     do iTUV = 1,   20
      AuxArray(iP,iTUV) = AuxArray(iP,iTUV) + TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,21) = AuxArray(iP,21) + Xpa*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(iP,22) = AuxArray(iP,22) + Ypa*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     AuxArray(iP,23) = AuxArray(iP,23) + Zpa*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     AuxArray(iP,24) = AuxArray(iP,24) + Xpa*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(iP,25) = AuxArray(iP,25) + Ypa*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     AuxArray(iP,26) = AuxArray(iP,26) + Xpa*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(iP,27) = AuxArray(iP,27) + Xpa*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     AuxArray(iP,28) = AuxArray(iP,28) + Xpa*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     AuxArray(iP,29) = AuxArray(iP,29) + Xpa*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     AuxArray(iP,30) = AuxArray(iP,30) + Xpa*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     AuxArray(iP,31) = AuxArray(iP,31) + Ypa*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(iP,32) = AuxArray(iP,32) + Zpa*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     AuxArray(iP,33) = AuxArray(iP,33) + Ypa*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(iP,34) = AuxArray(iP,34) + Ypa*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     AuxArray(iP,35) = AuxArray(iP,35) + Zpa*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
   ENDDO !iPrimQ=1, nPrimQ
  ENDDO !iP = 1,nPrimP*nPassP
 end subroutine
end module
