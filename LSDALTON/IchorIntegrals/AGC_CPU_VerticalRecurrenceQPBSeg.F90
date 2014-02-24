MODULE AGC_CPU_OBS_VERTICALRECURRENCEMODBSeg
 use IchorPrecisionModule
  
 CONTAINS

subroutine VerticalRecurrenceCPUSeg1B(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,&
         & QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:4,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Bcenter(3)
  real(realk),intent(inout) :: AUXarray(4,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,alphaP,RJ000(0:1)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
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
  Bx = -Bcenter(1)
  By = -Bcenter(2)
  Bz = -Bcenter(3)
  DO iPassQ = 1,nPasses
   AUXarray(1,iPassQ)=0.0E0_realk
   AUXarray(2,iPassQ)=0.0E0_realk
   AUXarray(3,iPassQ)=0.0E0_realk
   AUXarray(4,iPassQ)=0.0E0_realk
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    Xpb = Pcent(1,iPrimP) + Bx
    Ypb = Pcent(2,iPrimP) + By
    Zpb = Pcent(3,iPrimP) + Bz
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
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
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac
     TMP1 = PREF*RJ000(0)
     TMP2 = PREF*RJ000(1)
     AUXarray(1,iPassQ) = AUXarray(1,iPassQ) + TMP1
     AUXarray(2,iPassQ) = AUXarray(2,iPassQ) + Xpb*TMP1 + alphaXpq*TMP2
     AUXarray(3,iPassQ) = AUXarray(3,iPassQ) + Ypb*TMP1 + alphaYpq*TMP2
     AUXarray(4,iPassQ) = AUXarray(4,iPassQ) + Zpb*TMP1 + alphaZpq*TMP2
    enddo
   enddo
  enddo
end subroutine

subroutine VerticalRecurrenceCPUSeg2B(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 5,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Bcenter(3)
  real(realk),intent(inout) :: AUXarray(   10,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(    4)
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP,RJ000(0: 2)
  real(realk) :: TwoTerms(   1)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D2=2.0E0_realk
  REAL(REALK),PARAMETER :: D2JP36=  4.0000000000000000E+01_realk
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
  real(realk) :: TMParray1(  1:  1,2:3)
  real(realk) :: TMParray2(  2:  4,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Bx = -Bcenter(1)
  By = -Bcenter(2)
  Bz = -Bcenter(3)
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,   10
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpb = Pcent(1,iPrimP) + Bx
    Ypb = Pcent(2,iPrimP) + By
    Zpb = Pcent(3,iPrimP) + Bz
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000( 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000( 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000( 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000( 1) = RWVAL*(( 1 - D05)*RJ000( 0)-REXPW)
      RJ000( 2) = RWVAL*(( 2 - D05)*RJ000( 1)-REXPW)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
     ENDIF
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac
     TMPAuxarray(1) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMPAuxArray(2) = Xpb*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypb*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpb*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(1) + alphaP*TmpArray1(1,2))
     do iTUV = 1,    4
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(5,IPassQ) = AuxArray(5,IPassQ) + Xpb*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IPassQ) = AuxArray(6,IPassQ) + Xpb*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IPassQ) = AuxArray(7,IPassQ) + Xpb*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IPassQ) = AuxArray(8,IPassQ) + Ypb*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IPassQ) = AuxArray(9,IPassQ) + Ypb*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IPassQ) = AuxArray(10,IPassQ) + Zpb*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg3B(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 6,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Bcenter(3)
  real(realk),intent(inout) :: AUXarray(   20,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   10)
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP,RJ000(0: 3)
  real(realk) :: TwoTerms(   3)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D2=2.0E0_realk
  REAL(REALK),PARAMETER :: D2JP36=  4.2000000000000000E+01_realk
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
  real(realk) :: TMParray1(  1:  1,2:4)
  real(realk) :: TMParray2(  2:  4,2:3)
  real(realk) :: TMParray3(  5: 10,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Bx = -Bcenter(1)
  By = -Bcenter(2)
  Bz = -Bcenter(3)
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,   20
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpb = Pcent(1,iPrimP) + Bx
    Ypb = Pcent(2,iPrimP) + By
    Zpb = Pcent(3,iPrimP) + Bz
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000( 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000( 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000( 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000( 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000( 1) = RWVAL*(( 1 - D05)*RJ000( 0)-REXPW)
      RJ000( 2) = RWVAL*(( 2 - D05)*RJ000( 1)-REXPW)
      RJ000( 3) = RWVAL*(( 3 - D05)*RJ000( 2)-REXPW)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000( 3) = RWVAL*( 3 - D05)*RJ000( 2)
     ENDIF
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac
     TMPAuxarray(1) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
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
     do iTUV = 1,   10
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(11,IPassQ) = AuxArray(11,IPassQ) + Xpb*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IPassQ) = AuxArray(12,IPassQ) + Ypb*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IPassQ) = AuxArray(13,IPassQ) + Zpb*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IPassQ) = AuxArray(14,IPassQ) + Xpb*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IPassQ) = AuxArray(15,IPassQ) + Xpb*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IPassQ) = AuxArray(16,IPassQ) + Xpb*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IPassQ) = AuxArray(17,IPassQ) + Ypb*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IPassQ) = AuxArray(18,IPassQ) + Zpb*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IPassQ) = AuxArray(19,IPassQ) + Ypb*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IPassQ) = AuxArray(20,IPassQ) + Zpb*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg4B(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 7,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Bcenter(3)
  real(realk),intent(inout) :: AUXarray(   35,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   20)
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP,RJ000(0: 4)
  real(realk) :: TwoTerms(   6)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D2=2.0E0_realk
  REAL(REALK),PARAMETER :: D2JP36=  4.4000000000000000E+01_realk
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
  real(realk) :: TMParray1(  1:  1,2:5)
  real(realk) :: TMParray2(  2:  4,2:4)
  real(realk) :: TMParray3(  5: 10,2:3)
  real(realk) :: TMParray4( 11: 20,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Bx = -Bcenter(1)
  By = -Bcenter(2)
  Bz = -Bcenter(3)
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,   35
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpb = Pcent(1,iPrimP) + Bx
    Ypb = Pcent(2,iPrimP) + By
    Zpb = Pcent(3,iPrimP) + Bz
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000( 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000( 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000( 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000( 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000( 4) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000( 1) = RWVAL*(( 1 - D05)*RJ000( 0)-REXPW)
      RJ000( 2) = RWVAL*(( 2 - D05)*RJ000( 1)-REXPW)
      RJ000( 3) = RWVAL*(( 3 - D05)*RJ000( 2)-REXPW)
      RJ000( 4) = RWVAL*(( 4 - D05)*RJ000( 3)-REXPW)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000( 3) = RWVAL*( 3 - D05)*RJ000( 2)
      RJ000( 4) = RWVAL*( 4 - D05)*RJ000( 3)
     ENDIF
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac
     TMPAuxarray(1) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
     TMParray1(1, 5) = PREF*RJ000( 4)
     TMPAuxArray(2) = Xpb*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypb*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpb*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpb*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypb*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpb*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpb*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypb*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpb*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xpb*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Ypb*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zpb*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
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
     TwoTerms(1) = inv2expP*(TmpArray1(1,3) + alphaP*TmpArray1(1,4))
     tmpArray3(5,3) = Xpb*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xpb*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xpb*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Ypb*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Ypb*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zpb*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TMPAuxArray(2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(4) + alphaP*TmpArray2(4,2))
     TMPAuxArray(11) = Xpb*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Ypb*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zpb*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xpb*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xpb*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xpb*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Ypb*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zpb*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Ypb*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zpb*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(TMPAuxArray(5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(10) + alphaP*TmpArray3(10,2))
     do iTUV = 1,   20
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(21,IPassQ) = AuxArray(21,IPassQ) + Xpb*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IPassQ) = AuxArray(22,IPassQ) + Ypb*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IPassQ) = AuxArray(23,IPassQ) + Zpb*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IPassQ) = AuxArray(24,IPassQ) + Xpb*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IPassQ) = AuxArray(25,IPassQ) + Ypb*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IPassQ) = AuxArray(26,IPassQ) + Xpb*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IPassQ) = AuxArray(27,IPassQ) + Xpb*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IPassQ) = AuxArray(28,IPassQ) + Xpb*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IPassQ) = AuxArray(29,IPassQ) + Xpb*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IPassQ) = AuxArray(30,IPassQ) + Xpb*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IPassQ) = AuxArray(31,IPassQ) + Ypb*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IPassQ) = AuxArray(32,IPassQ) + Zpb*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IPassQ) = AuxArray(33,IPassQ) + Ypb*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IPassQ) = AuxArray(34,IPassQ) + Ypb*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IPassQ) = AuxArray(35,IPassQ) + Zpb*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg5B(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 8,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Bcenter(3)
  real(realk),intent(inout) :: AUXarray(   56,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   35)
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP,RJ000(0: 5)
  real(realk) :: TwoTerms(  10)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D2=2.0E0_realk
  REAL(REALK),PARAMETER :: D2JP36=  4.6000000000000000E+01_realk
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
  real(realk) :: TMParray1(  1:  1,2:6)
  real(realk) :: TMParray2(  2:  4,2:5)
  real(realk) :: TMParray3(  5: 10,2:4)
  real(realk) :: TMParray4( 11: 20,2:3)
  real(realk) :: TMParray5( 21: 35,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Bx = -Bcenter(1)
  By = -Bcenter(2)
  Bz = -Bcenter(3)
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,   56
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpb = Pcent(1,iPrimP) + Bx
    Ypb = Pcent(2,iPrimP) + By
    Zpb = Pcent(3,iPrimP) + Bz
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000( 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000( 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000( 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000( 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000( 4) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000( 5) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000( 1) = RWVAL*(( 1 - D05)*RJ000( 0)-REXPW)
      RJ000( 2) = RWVAL*(( 2 - D05)*RJ000( 1)-REXPW)
      RJ000( 3) = RWVAL*(( 3 - D05)*RJ000( 2)-REXPW)
      RJ000( 4) = RWVAL*(( 4 - D05)*RJ000( 3)-REXPW)
      RJ000( 5) = RWVAL*(( 5 - D05)*RJ000( 4)-REXPW)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000( 3) = RWVAL*( 3 - D05)*RJ000( 2)
      RJ000( 4) = RWVAL*( 4 - D05)*RJ000( 3)
      RJ000( 5) = RWVAL*( 5 - D05)*RJ000( 4)
     ENDIF
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac
     TMPAuxarray(1) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
     TMParray1(1, 5) = PREF*RJ000( 4)
     TMParray1(1, 6) = PREF*RJ000( 5)
     TMPAuxArray(2) = Xpb*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypb*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpb*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expP*(TMPAuxArray(2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(4) + alphaP*TmpArray2(4,2))
     TMPAuxArray(11) = Xpb*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Ypb*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zpb*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xpb*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xpb*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xpb*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Ypb*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zpb*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Ypb*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zpb*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(TMPAuxArray(5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(10) + alphaP*TmpArray3(10,2))
     TMPAuxArray(21) = Xpb*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Ypb*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zpb*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xpb*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Ypb*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xpb*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xpb*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xpb*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xpb*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xpb*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Ypb*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zpb*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Ypb*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Ypb*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zpb*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(TMPAuxArray(11) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(14) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(16) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(17) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(19) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(20) + alphaP*TmpArray4(20,2))
     do iTUV = 1,   35
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(36,IPassQ) = AuxArray(36,IPassQ) + Xpb*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IPassQ) = AuxArray(37,IPassQ) + Ypb*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IPassQ) = AuxArray(38,IPassQ) + Zpb*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IPassQ) = AuxArray(39,IPassQ) + Xpb*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IPassQ) = AuxArray(40,IPassQ) + Ypb*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IPassQ) = AuxArray(41,IPassQ) + Xpb*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IPassQ) = AuxArray(42,IPassQ) + Xpb*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IPassQ) = AuxArray(43,IPassQ) + Zpb*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IPassQ) = AuxArray(44,IPassQ) + Ypb*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IPassQ) = AuxArray(45,IPassQ) + Xpb*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IPassQ) = AuxArray(46,IPassQ) + Xpb*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IPassQ) = AuxArray(47,IPassQ) + Xpb*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IPassQ) = AuxArray(48,IPassQ) + Xpb*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IPassQ) = AuxArray(49,IPassQ) + Xpb*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IPassQ) = AuxArray(50,IPassQ) + Xpb*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IPassQ) = AuxArray(51,IPassQ) + Ypb*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IPassQ) = AuxArray(52,IPassQ) + Zpb*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IPassQ) = AuxArray(53,IPassQ) + Ypb*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IPassQ) = AuxArray(54,IPassQ) + Ypb*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IPassQ) = AuxArray(55,IPassQ) + Ypb*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IPassQ) = AuxArray(56,IPassQ) + Zpb*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg6B(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 9,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Bcenter(3)
  real(realk),intent(inout) :: AUXarray(   84,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   56)
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP,RJ000(0: 6)
  real(realk) :: TwoTerms(  15)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D2=2.0E0_realk
  REAL(REALK),PARAMETER :: D2JP36=  4.8000000000000000E+01_realk
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
  real(realk) :: TMParray1(  1:  1,2:7)
  real(realk) :: TMParray2(  2:  4,2:6)
  real(realk) :: TMParray3(  5: 10,2:5)
  real(realk) :: TMParray4( 11: 20,2:4)
  real(realk) :: TMParray5( 21: 35,2:3)
  real(realk) :: TMParray6( 36: 56,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Bx = -Bcenter(1)
  By = -Bcenter(2)
  Bz = -Bcenter(3)
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,   84
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpb = Pcent(1,iPrimP) + Bx
    Ypb = Pcent(2,iPrimP) + By
    Zpb = Pcent(3,iPrimP) + Bz
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000( 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000( 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000( 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000( 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000( 4) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000( 5) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
      RJ000( 6) = TABFJW( 6,IPNT)-TABFJW( 7,IPNT)*WDIFF+TABFJW( 8,IPNT)*W2+TABFJW( 9,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000( 1) = RWVAL*(( 1 - D05)*RJ000( 0)-REXPW)
      RJ000( 2) = RWVAL*(( 2 - D05)*RJ000( 1)-REXPW)
      RJ000( 3) = RWVAL*(( 3 - D05)*RJ000( 2)-REXPW)
      RJ000( 4) = RWVAL*(( 4 - D05)*RJ000( 3)-REXPW)
      RJ000( 5) = RWVAL*(( 5 - D05)*RJ000( 4)-REXPW)
      RJ000( 6) = RWVAL*(( 6 - D05)*RJ000( 5)-REXPW)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000( 3) = RWVAL*( 3 - D05)*RJ000( 2)
      RJ000( 4) = RWVAL*( 4 - D05)*RJ000( 3)
      RJ000( 5) = RWVAL*( 5 - D05)*RJ000( 4)
      RJ000( 6) = RWVAL*( 6 - D05)*RJ000( 5)
     ENDIF
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac
     TMPAuxarray(1) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
     TMParray1(1, 5) = PREF*RJ000( 4)
     TMParray1(1, 6) = PREF*RJ000( 5)
     TMParray1(1, 7) = PREF*RJ000( 6)
     TMPAuxArray(2) = Xpb*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypb*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpb*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expP*(TMPAuxArray(2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(4) + alphaP*TmpArray2(4,2))
     TMPAuxArray(11) = Xpb*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Ypb*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zpb*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xpb*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xpb*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xpb*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Ypb*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zpb*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Ypb*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zpb*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(TMPAuxArray(5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(10) + alphaP*TmpArray3(10,2))
     TMPAuxArray(21) = Xpb*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Ypb*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zpb*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xpb*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Ypb*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xpb*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xpb*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xpb*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xpb*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xpb*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Ypb*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zpb*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Ypb*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Ypb*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zpb*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(TMPAuxArray(11) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(14) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(16) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(17) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(19) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(20) + alphaP*TmpArray4(20,2))
     TMPAuxArray(36) = Xpb*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     TMPAuxArray(37) = Ypb*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     TMPAuxArray(38) = Zpb*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     TMPAuxArray(39) = Xpb*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     TMPAuxArray(40) = Ypb*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     TMPAuxArray(41) = Xpb*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     TMPAuxArray(42) = Xpb*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     TMPAuxArray(43) = Zpb*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     TMPAuxArray(44) = Ypb*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     TMPAuxArray(45) = Xpb*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     TMPAuxArray(46) = Xpb*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     TMPAuxArray(47) = Xpb*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     TMPAuxArray(48) = Xpb*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     TMPAuxArray(49) = Xpb*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     TMPAuxArray(50) = Xpb*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     TMPAuxArray(51) = Ypb*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     TMPAuxArray(52) = Zpb*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     TMPAuxArray(53) = Ypb*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     TMPAuxArray(54) = Ypb*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     TMPAuxArray(55) = Ypb*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     TMPAuxArray(56) = Zpb*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     TwoTerms(1) = inv2expP*(TMPAuxArray(21) + alphaP*TmpArray5(21,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(24) + alphaP*TmpArray5(24,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(26) + alphaP*TmpArray5(26,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(27) + alphaP*TmpArray5(27,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(30) + alphaP*TmpArray5(30,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(31) + alphaP*TmpArray5(31,2))
     TwoTerms(7) = inv2expP*(TMPAuxArray(33) + alphaP*TmpArray5(33,2))
     TwoTerms(8) = inv2expP*(TMPAuxArray(34) + alphaP*TmpArray5(34,2))
     TwoTerms(9) = inv2expP*(TMPAuxArray(35) + alphaP*TmpArray5(35,2))
     do iTUV = 1,   56
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(57,IPassQ) = AuxArray(57,IPassQ) + Xpb*TMPAuxArray(36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,IPassQ) = AuxArray(58,IPassQ) + Ypb*TMPAuxArray(36) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,IPassQ) = AuxArray(59,IPassQ) + Zpb*TMPAuxArray(36) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,IPassQ) = AuxArray(60,IPassQ) + Xpb*TMPAuxArray(39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,IPassQ) = AuxArray(61,IPassQ) + Ypb*TMPAuxArray(38) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,IPassQ) = AuxArray(62,IPassQ) + Xpb*TMPAuxArray(41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,IPassQ) = AuxArray(63,IPassQ) + Xpb*TMPAuxArray(42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,IPassQ) = AuxArray(64,IPassQ) + Zpb*TMPAuxArray(39) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,IPassQ) = AuxArray(65,IPassQ) + Ypb*TMPAuxArray(41) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,IPassQ) = AuxArray(66,IPassQ) + Xpb*TMPAuxArray(45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,IPassQ) = AuxArray(67,IPassQ) + Xpb*TMPAuxArray(46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,IPassQ) = AuxArray(68,IPassQ) + Zpb*TMPAuxArray(42) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,IPassQ) = AuxArray(69,IPassQ) + Xpb*TMPAuxArray(48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,IPassQ) = AuxArray(70,IPassQ) + Ypb*TMPAuxArray(45) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,IPassQ) = AuxArray(71,IPassQ) + Xpb*TMPAuxArray(50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,IPassQ) = AuxArray(72,IPassQ) + Xpb*TMPAuxArray(51) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,IPassQ) = AuxArray(73,IPassQ) + Xpb*TMPAuxArray(52) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,IPassQ) = AuxArray(74,IPassQ) + Xpb*TMPAuxArray(53) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,IPassQ) = AuxArray(75,IPassQ) + Xpb*TMPAuxArray(54) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,IPassQ) = AuxArray(76,IPassQ) + Xpb*TMPAuxArray(55) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,IPassQ) = AuxArray(77,IPassQ) + Xpb*TMPAuxArray(56) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,IPassQ) = AuxArray(78,IPassQ) + Ypb*TMPAuxArray(51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,IPassQ) = AuxArray(79,IPassQ) + Zpb*TMPAuxArray(51) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,IPassQ) = AuxArray(80,IPassQ) + Ypb*TMPAuxArray(53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,IPassQ) = AuxArray(81,IPassQ) + Ypb*TMPAuxArray(54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,IPassQ) = AuxArray(82,IPassQ) + Ypb*TMPAuxArray(55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,IPassQ) = AuxArray(83,IPassQ) + Ypb*TMPAuxArray(56) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,IPassQ) = AuxArray(84,IPassQ) + Zpb*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg7B(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:10,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Bcenter(3)
  real(realk),intent(inout) :: AUXarray(  120,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   84)
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP,RJ000(0: 7)
  real(realk) :: TwoTerms(  21)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D2=2.0E0_realk
  REAL(REALK),PARAMETER :: D2JP36=  5.0000000000000000E+01_realk
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
  real(realk) :: TMParray1(  1:  1,2:8)
  real(realk) :: TMParray2(  2:  4,2:7)
  real(realk) :: TMParray3(  5: 10,2:6)
  real(realk) :: TMParray4( 11: 20,2:5)
  real(realk) :: TMParray5( 21: 35,2:4)
  real(realk) :: TMParray6( 36: 56,2:3)
  real(realk) :: TMParray7( 57: 84,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Bx = -Bcenter(1)
  By = -Bcenter(2)
  Bz = -Bcenter(3)
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,  120
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpb = Pcent(1,iPrimP) + Bx
    Ypb = Pcent(2,iPrimP) + By
    Zpb = Pcent(3,iPrimP) + Bz
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000( 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000( 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000( 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000( 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000( 4) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000( 5) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
      RJ000( 6) = TABFJW( 6,IPNT)-TABFJW( 7,IPNT)*WDIFF+TABFJW( 8,IPNT)*W2+TABFJW( 9,IPNT)*W3
      RJ000( 7) = TABFJW( 7,IPNT)-TABFJW( 8,IPNT)*WDIFF+TABFJW( 9,IPNT)*W2+TABFJW(10,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000( 1) = RWVAL*(( 1 - D05)*RJ000( 0)-REXPW)
      RJ000( 2) = RWVAL*(( 2 - D05)*RJ000( 1)-REXPW)
      RJ000( 3) = RWVAL*(( 3 - D05)*RJ000( 2)-REXPW)
      RJ000( 4) = RWVAL*(( 4 - D05)*RJ000( 3)-REXPW)
      RJ000( 5) = RWVAL*(( 5 - D05)*RJ000( 4)-REXPW)
      RJ000( 6) = RWVAL*(( 6 - D05)*RJ000( 5)-REXPW)
      RJ000( 7) = RWVAL*(( 7 - D05)*RJ000( 6)-REXPW)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000( 3) = RWVAL*( 3 - D05)*RJ000( 2)
      RJ000( 4) = RWVAL*( 4 - D05)*RJ000( 3)
      RJ000( 5) = RWVAL*( 5 - D05)*RJ000( 4)
      RJ000( 6) = RWVAL*( 6 - D05)*RJ000( 5)
      RJ000( 7) = RWVAL*( 7 - D05)*RJ000( 6)
     ENDIF
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac
     TMPAuxarray(1) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
     TMParray1(1, 5) = PREF*RJ000( 4)
     TMParray1(1, 6) = PREF*RJ000( 5)
     TMParray1(1, 7) = PREF*RJ000( 6)
     TMParray1(1, 8) = PREF*RJ000( 7)
     TMPAuxArray(2) = Xpb*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypb*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpb*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
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
     tmpArray2(2,7) = Xpb*tmpArray1(1,7) + alphaXpq*TmpArray1(1,8)
     tmpArray2(3,7) = Ypb*tmpArray1(1,7) + alphaYpq*TmpArray1(1,8)
     tmpArray2(4,7) = Zpb*tmpArray1(1,7) + alphaZpq*TmpArray1(1,8)
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
     TwoTerms(1) = inv2expP*(TmpArray1(1,6) + alphaP*TmpArray1(1,7))
     tmpArray3(5,6) = Xpb*tmpArray2(2,6) + alphaXpq*TmpArray2(2,7) + TwoTerms(1)
     tmpArray3(6,6) = Xpb*tmpArray2(3,6) + alphaXpq*TmpArray2(3,7)
     tmpArray3(7,6) = Xpb*tmpArray2(4,6) + alphaXpq*TmpArray2(4,7)
     tmpArray3(8,6) = Ypb*tmpArray2(3,6) + alphaYpq*TmpArray2(3,7) + TwoTerms(1)
     tmpArray3(9,6) = Ypb*tmpArray2(4,6) + alphaYpq*TmpArray2(4,7)
     tmpArray3(10,6) = Zpb*tmpArray2(4,6) + alphaZpq*TmpArray2(4,7) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TMPAuxArray(2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(4) + alphaP*TmpArray2(4,2))
     TMPAuxArray(11) = Xpb*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Ypb*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zpb*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xpb*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xpb*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xpb*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Ypb*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zpb*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Ypb*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zpb*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(TmpArray2(2,5) + alphaP*TmpArray2(2,6))
     TwoTerms(2) = inv2expP*(TmpArray2(3,5) + alphaP*TmpArray2(3,6))
     TwoTerms(3) = inv2expP*(TmpArray2(4,5) + alphaP*TmpArray2(4,6))
     tmpArray4(11,5) = Xpb*tmpArray3(5,5) + alphaXpq*TmpArray3(5,6) + 2*TwoTerms(1)
     tmpArray4(12,5) = Ypb*tmpArray3(5,5) + alphaYpq*TmpArray3(5,6)
     tmpArray4(13,5) = Zpb*tmpArray3(5,5) + alphaZpq*TmpArray3(5,6)
     tmpArray4(14,5) = Xpb*tmpArray3(8,5) + alphaXpq*TmpArray3(8,6)
     tmpArray4(15,5) = Xpb*tmpArray3(9,5) + alphaXpq*TmpArray3(9,6)
     tmpArray4(16,5) = Xpb*tmpArray3(10,5) + alphaXpq*TmpArray3(10,6)
     tmpArray4(17,5) = Ypb*tmpArray3(8,5) + alphaYpq*TmpArray3(8,6) + 2*TwoTerms(2)
     tmpArray4(18,5) = Zpb*tmpArray3(8,5) + alphaZpq*TmpArray3(8,6)
     tmpArray4(19,5) = Ypb*tmpArray3(10,5) + alphaYpq*TmpArray3(10,6)
     tmpArray4(20,5) = Zpb*tmpArray3(10,5) + alphaZpq*TmpArray3(10,6) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(10) + alphaP*TmpArray3(10,2))
     TMPAuxArray(21) = Xpb*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Ypb*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zpb*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xpb*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Ypb*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xpb*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xpb*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xpb*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xpb*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xpb*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Ypb*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zpb*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Ypb*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Ypb*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zpb*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(TmpArray3(5,4) + alphaP*TmpArray3(5,5))
     TwoTerms(2) = inv2expP*(TmpArray3(8,4) + alphaP*TmpArray3(8,5))
     TwoTerms(3) = inv2expP*(TmpArray3(10,4) + alphaP*TmpArray3(10,5))
     tmpArray5(21,4) = Xpb*tmpArray4(11,4) + alphaXpq*TmpArray4(11,5) + 3*TwoTerms(1)
     tmpArray5(22,4) = Ypb*tmpArray4(11,4) + alphaYpq*TmpArray4(11,5)
     tmpArray5(23,4) = Zpb*tmpArray4(11,4) + alphaZpq*TmpArray4(11,5)
     tmpArray5(24,4) = Xpb*tmpArray4(14,4) + alphaXpq*TmpArray4(14,5) + TwoTerms(2)
     tmpArray5(25,4) = Ypb*tmpArray4(13,4) + alphaYpq*TmpArray4(13,5)
     tmpArray5(26,4) = Xpb*tmpArray4(16,4) + alphaXpq*TmpArray4(16,5) + TwoTerms(3)
     tmpArray5(27,4) = Xpb*tmpArray4(17,4) + alphaXpq*TmpArray4(17,5)
     tmpArray5(28,4) = Xpb*tmpArray4(18,4) + alphaXpq*TmpArray4(18,5)
     tmpArray5(29,4) = Xpb*tmpArray4(19,4) + alphaXpq*TmpArray4(19,5)
     tmpArray5(30,4) = Xpb*tmpArray4(20,4) + alphaXpq*TmpArray4(20,5)
     tmpArray5(31,4) = Ypb*tmpArray4(17,4) + alphaYpq*TmpArray4(17,5) + 3*TwoTerms(2)
     tmpArray5(32,4) = Zpb*tmpArray4(17,4) + alphaZpq*TmpArray4(17,5)
     tmpArray5(33,4) = Ypb*tmpArray4(19,4) + alphaYpq*TmpArray4(19,5) + TwoTerms(3)
     tmpArray5(34,4) = Ypb*tmpArray4(20,4) + alphaYpq*TmpArray4(20,5)
     tmpArray5(35,4) = Zpb*tmpArray4(20,4) + alphaZpq*TmpArray4(20,5) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(11) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(14) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(16) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(17) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(19) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(20) + alphaP*TmpArray4(20,2))
     TMPAuxArray(36) = Xpb*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     TMPAuxArray(37) = Ypb*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     TMPAuxArray(38) = Zpb*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     TMPAuxArray(39) = Xpb*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     TMPAuxArray(40) = Ypb*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     TMPAuxArray(41) = Xpb*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     TMPAuxArray(42) = Xpb*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     TMPAuxArray(43) = Zpb*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     TMPAuxArray(44) = Ypb*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     TMPAuxArray(45) = Xpb*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     TMPAuxArray(46) = Xpb*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     TMPAuxArray(47) = Xpb*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     TMPAuxArray(48) = Xpb*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     TMPAuxArray(49) = Xpb*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     TMPAuxArray(50) = Xpb*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     TMPAuxArray(51) = Ypb*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     TMPAuxArray(52) = Zpb*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     TMPAuxArray(53) = Ypb*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     TMPAuxArray(54) = Ypb*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     TMPAuxArray(55) = Ypb*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     TMPAuxArray(56) = Zpb*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     TwoTerms(1) = inv2expP*(TmpArray4(11,3) + alphaP*TmpArray4(11,4))
     TwoTerms(2) = inv2expP*(TmpArray4(14,3) + alphaP*TmpArray4(14,4))
     TwoTerms(3) = inv2expP*(TmpArray4(16,3) + alphaP*TmpArray4(16,4))
     TwoTerms(4) = inv2expP*(TmpArray4(17,3) + alphaP*TmpArray4(17,4))
     TwoTerms(5) = inv2expP*(TmpArray4(19,3) + alphaP*TmpArray4(19,4))
     TwoTerms(6) = inv2expP*(TmpArray4(20,3) + alphaP*TmpArray4(20,4))
     tmpArray6(36,3) = Xpb*tmpArray5(21,3) + alphaXpq*TmpArray5(21,4) + 4*TwoTerms(1)
     tmpArray6(37,3) = Ypb*tmpArray5(21,3) + alphaYpq*TmpArray5(21,4)
     tmpArray6(38,3) = Zpb*tmpArray5(21,3) + alphaZpq*TmpArray5(21,4)
     tmpArray6(39,3) = Xpb*tmpArray5(24,3) + alphaXpq*TmpArray5(24,4) + 2*TwoTerms(2)
     tmpArray6(40,3) = Ypb*tmpArray5(23,3) + alphaYpq*TmpArray5(23,4)
     tmpArray6(41,3) = Xpb*tmpArray5(26,3) + alphaXpq*TmpArray5(26,4) + 2*TwoTerms(3)
     tmpArray6(42,3) = Xpb*tmpArray5(27,3) + alphaXpq*TmpArray5(27,4) + TwoTerms(4)
     tmpArray6(43,3) = Zpb*tmpArray5(24,3) + alphaZpq*TmpArray5(24,4)
     tmpArray6(44,3) = Ypb*tmpArray5(26,3) + alphaYpq*TmpArray5(26,4)
     tmpArray6(45,3) = Xpb*tmpArray5(30,3) + alphaXpq*TmpArray5(30,4) + TwoTerms(6)
     tmpArray6(46,3) = Xpb*tmpArray5(31,3) + alphaXpq*TmpArray5(31,4)
     tmpArray6(47,3) = Xpb*tmpArray5(32,3) + alphaXpq*TmpArray5(32,4)
     tmpArray6(48,3) = Xpb*tmpArray5(33,3) + alphaXpq*TmpArray5(33,4)
     tmpArray6(49,3) = Xpb*tmpArray5(34,3) + alphaXpq*TmpArray5(34,4)
     tmpArray6(50,3) = Xpb*tmpArray5(35,3) + alphaXpq*TmpArray5(35,4)
     tmpArray6(51,3) = Ypb*tmpArray5(31,3) + alphaYpq*TmpArray5(31,4) + 4*TwoTerms(4)
     tmpArray6(52,3) = Zpb*tmpArray5(31,3) + alphaZpq*TmpArray5(31,4)
     tmpArray6(53,3) = Ypb*tmpArray5(33,3) + alphaYpq*TmpArray5(33,4) + 2*TwoTerms(5)
     tmpArray6(54,3) = Ypb*tmpArray5(34,3) + alphaYpq*TmpArray5(34,4) + TwoTerms(6)
     tmpArray6(55,3) = Ypb*tmpArray5(35,3) + alphaYpq*TmpArray5(35,4)
     tmpArray6(56,3) = Zpb*tmpArray5(35,3) + alphaZpq*TmpArray5(35,4) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TMPAuxArray(21) + alphaP*TmpArray5(21,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(24) + alphaP*TmpArray5(24,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(26) + alphaP*TmpArray5(26,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(27) + alphaP*TmpArray5(27,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(30) + alphaP*TmpArray5(30,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(31) + alphaP*TmpArray5(31,2))
     TwoTerms(7) = inv2expP*(TMPAuxArray(33) + alphaP*TmpArray5(33,2))
     TwoTerms(8) = inv2expP*(TMPAuxArray(34) + alphaP*TmpArray5(34,2))
     TwoTerms(9) = inv2expP*(TMPAuxArray(35) + alphaP*TmpArray5(35,2))
     TMPAuxArray(57) = Xpb*TMPAuxArray(36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     TMPAuxArray(58) = Ypb*TMPAuxArray(36) + alphaYpq*TmpArray6(36,2)
     TMPAuxArray(59) = Zpb*TMPAuxArray(36) + alphaZpq*TmpArray6(36,2)
     TMPAuxArray(60) = Xpb*TMPAuxArray(39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     TMPAuxArray(61) = Ypb*TMPAuxArray(38) + alphaYpq*TmpArray6(38,2)
     TMPAuxArray(62) = Xpb*TMPAuxArray(41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     TMPAuxArray(63) = Xpb*TMPAuxArray(42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     TMPAuxArray(64) = Zpb*TMPAuxArray(39) + alphaZpq*TmpArray6(39,2)
     TMPAuxArray(65) = Ypb*TMPAuxArray(41) + alphaYpq*TmpArray6(41,2)
     TMPAuxArray(66) = Xpb*TMPAuxArray(45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     TMPAuxArray(67) = Xpb*TMPAuxArray(46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     TMPAuxArray(68) = Zpb*TMPAuxArray(42) + alphaZpq*TmpArray6(42,2)
     TMPAuxArray(69) = Xpb*TMPAuxArray(48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     TMPAuxArray(70) = Ypb*TMPAuxArray(45) + alphaYpq*TmpArray6(45,2)
     TMPAuxArray(71) = Xpb*TMPAuxArray(50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     TMPAuxArray(72) = Xpb*TMPAuxArray(51) + alphaXpq*TmpArray6(51,2)
     TMPAuxArray(73) = Xpb*TMPAuxArray(52) + alphaXpq*TmpArray6(52,2)
     TMPAuxArray(74) = Xpb*TMPAuxArray(53) + alphaXpq*TmpArray6(53,2)
     TMPAuxArray(75) = Xpb*TMPAuxArray(54) + alphaXpq*TmpArray6(54,2)
     TMPAuxArray(76) = Xpb*TMPAuxArray(55) + alphaXpq*TmpArray6(55,2)
     TMPAuxArray(77) = Xpb*TMPAuxArray(56) + alphaXpq*TmpArray6(56,2)
     TMPAuxArray(78) = Ypb*TMPAuxArray(51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     TMPAuxArray(79) = Zpb*TMPAuxArray(51) + alphaZpq*TmpArray6(51,2)
     TMPAuxArray(80) = Ypb*TMPAuxArray(53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     TMPAuxArray(81) = Ypb*TMPAuxArray(54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     TMPAuxArray(82) = Ypb*TMPAuxArray(55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     TMPAuxArray(83) = Ypb*TMPAuxArray(56) + alphaYpq*TmpArray6(56,2)
     TMPAuxArray(84) = Zpb*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expP*(TmpArray5(21,2) + alphaP*TmpArray5(21,3))
     TwoTerms(2) = inv2expP*(TmpArray5(24,2) + alphaP*TmpArray5(24,3))
     TwoTerms(3) = inv2expP*(TmpArray5(26,2) + alphaP*TmpArray5(26,3))
     TwoTerms(4) = inv2expP*(TmpArray5(27,2) + alphaP*TmpArray5(27,3))
     TwoTerms(5) = inv2expP*(TmpArray5(30,2) + alphaP*TmpArray5(30,3))
     TwoTerms(6) = inv2expP*(TmpArray5(31,2) + alphaP*TmpArray5(31,3))
     TwoTerms(7) = inv2expP*(TmpArray5(33,2) + alphaP*TmpArray5(33,3))
     TwoTerms(8) = inv2expP*(TmpArray5(34,2) + alphaP*TmpArray5(34,3))
     TwoTerms(9) = inv2expP*(TmpArray5(35,2) + alphaP*TmpArray5(35,3))
     tmpArray7(57,2) = Xpb*tmpArray6(36,2) + alphaXpq*TmpArray6(36,3) + 5*TwoTerms(1)
     tmpArray7(58,2) = Ypb*tmpArray6(36,2) + alphaYpq*TmpArray6(36,3)
     tmpArray7(59,2) = Zpb*tmpArray6(36,2) + alphaZpq*TmpArray6(36,3)
     tmpArray7(60,2) = Xpb*tmpArray6(39,2) + alphaXpq*TmpArray6(39,3) + 3*TwoTerms(2)
     tmpArray7(61,2) = Ypb*tmpArray6(38,2) + alphaYpq*TmpArray6(38,3)
     tmpArray7(62,2) = Xpb*tmpArray6(41,2) + alphaXpq*TmpArray6(41,3) + 3*TwoTerms(3)
     tmpArray7(63,2) = Xpb*tmpArray6(42,2) + alphaXpq*TmpArray6(42,3) + 2*TwoTerms(4)
     tmpArray7(64,2) = Zpb*tmpArray6(39,2) + alphaZpq*TmpArray6(39,3)
     tmpArray7(65,2) = Ypb*tmpArray6(41,2) + alphaYpq*TmpArray6(41,3)
     tmpArray7(66,2) = Xpb*tmpArray6(45,2) + alphaXpq*TmpArray6(45,3) + 2*TwoTerms(5)
     tmpArray7(67,2) = Xpb*tmpArray6(46,2) + alphaXpq*TmpArray6(46,3) + TwoTerms(6)
     tmpArray7(68,2) = Zpb*tmpArray6(42,2) + alphaZpq*TmpArray6(42,3)
     tmpArray7(69,2) = Xpb*tmpArray6(48,2) + alphaXpq*TmpArray6(48,3) + TwoTerms(7)
     tmpArray7(70,2) = Ypb*tmpArray6(45,2) + alphaYpq*TmpArray6(45,3)
     tmpArray7(71,2) = Xpb*tmpArray6(50,2) + alphaXpq*TmpArray6(50,3) + TwoTerms(9)
     tmpArray7(72,2) = Xpb*tmpArray6(51,2) + alphaXpq*TmpArray6(51,3)
     tmpArray7(73,2) = Xpb*tmpArray6(52,2) + alphaXpq*TmpArray6(52,3)
     tmpArray7(74,2) = Xpb*tmpArray6(53,2) + alphaXpq*TmpArray6(53,3)
     tmpArray7(75,2) = Xpb*tmpArray6(54,2) + alphaXpq*TmpArray6(54,3)
     tmpArray7(76,2) = Xpb*tmpArray6(55,2) + alphaXpq*TmpArray6(55,3)
     tmpArray7(77,2) = Xpb*tmpArray6(56,2) + alphaXpq*TmpArray6(56,3)
     tmpArray7(78,2) = Ypb*tmpArray6(51,2) + alphaYpq*TmpArray6(51,3) + 5*TwoTerms(6)
     tmpArray7(79,2) = Zpb*tmpArray6(51,2) + alphaZpq*TmpArray6(51,3)
     tmpArray7(80,2) = Ypb*tmpArray6(53,2) + alphaYpq*TmpArray6(53,3) + 3*TwoTerms(7)
     tmpArray7(81,2) = Ypb*tmpArray6(54,2) + alphaYpq*TmpArray6(54,3) + 2*TwoTerms(8)
     tmpArray7(82,2) = Ypb*tmpArray6(55,2) + alphaYpq*TmpArray6(55,3) + TwoTerms(9)
     tmpArray7(83,2) = Ypb*tmpArray6(56,2) + alphaYpq*TmpArray6(56,3)
     tmpArray7(84,2) = Zpb*tmpArray6(56,2) + alphaZpq*TmpArray6(56,3) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expP*(TMPAuxArray(36) + alphaP*TmpArray6(36,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(39) + alphaP*TmpArray6(39,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(41) + alphaP*TmpArray6(41,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(42) + alphaP*TmpArray6(42,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(45) + alphaP*TmpArray6(45,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(46) + alphaP*TmpArray6(46,2))
     TwoTerms(7) = inv2expP*(TMPAuxArray(48) + alphaP*TmpArray6(48,2))
     TwoTerms(8) = inv2expP*(TMPAuxArray(50) + alphaP*TmpArray6(50,2))
     TwoTerms(9) = inv2expP*(TMPAuxArray(51) + alphaP*TmpArray6(51,2))
     TwoTerms(10) = inv2expP*(TMPAuxArray(53) + alphaP*TmpArray6(53,2))
     TwoTerms(11) = inv2expP*(TMPAuxArray(54) + alphaP*TmpArray6(54,2))
     TwoTerms(12) = inv2expP*(TMPAuxArray(55) + alphaP*TmpArray6(55,2))
     TwoTerms(13) = inv2expP*(TMPAuxArray(56) + alphaP*TmpArray6(56,2))
     do iTUV = 1,   84
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(85,IPassQ) = AuxArray(85,IPassQ) + Xpb*TMPAuxArray(57) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     AuxArray(86,IPassQ) = AuxArray(86,IPassQ) + Ypb*TMPAuxArray(57) + alphaYpq*TmpArray7(57,2)
     AuxArray(87,IPassQ) = AuxArray(87,IPassQ) + Zpb*TMPAuxArray(57) + alphaZpq*TmpArray7(57,2)
     AuxArray(88,IPassQ) = AuxArray(88,IPassQ) + Xpb*TMPAuxArray(60) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     AuxArray(89,IPassQ) = AuxArray(89,IPassQ) + Ypb*TMPAuxArray(59) + alphaYpq*TmpArray7(59,2)
     AuxArray(90,IPassQ) = AuxArray(90,IPassQ) + Xpb*TMPAuxArray(62) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     AuxArray(91,IPassQ) = AuxArray(91,IPassQ) + Xpb*TMPAuxArray(63) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     AuxArray(92,IPassQ) = AuxArray(92,IPassQ) + Zpb*TMPAuxArray(60) + alphaZpq*TmpArray7(60,2)
     AuxArray(93,IPassQ) = AuxArray(93,IPassQ) + Ypb*TMPAuxArray(62) + alphaYpq*TmpArray7(62,2)
     AuxArray(94,IPassQ) = AuxArray(94,IPassQ) + Xpb*TMPAuxArray(66) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     AuxArray(95,IPassQ) = AuxArray(95,IPassQ) + Xpb*TMPAuxArray(67) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     AuxArray(96,IPassQ) = AuxArray(96,IPassQ) + Zpb*TMPAuxArray(63) + alphaZpq*TmpArray7(63,2)
     AuxArray(97,IPassQ) = AuxArray(97,IPassQ) + Xpb*TMPAuxArray(69) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     AuxArray(98,IPassQ) = AuxArray(98,IPassQ) + Ypb*TMPAuxArray(66) + alphaYpq*TmpArray7(66,2)
     AuxArray(99,IPassQ) = AuxArray(99,IPassQ) + Xpb*TMPAuxArray(71) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     AuxArray(100,IPassQ) = AuxArray(100,IPassQ) + Xpb*TMPAuxArray(72) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     AuxArray(101,IPassQ) = AuxArray(101,IPassQ) + Zpb*TMPAuxArray(67) + alphaZpq*TmpArray7(67,2)
     AuxArray(102,IPassQ) = AuxArray(102,IPassQ) + Xpb*TMPAuxArray(74) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     AuxArray(103,IPassQ) = AuxArray(103,IPassQ) + Xpb*TMPAuxArray(75) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     AuxArray(104,IPassQ) = AuxArray(104,IPassQ) + Ypb*TMPAuxArray(71) + alphaYpq*TmpArray7(71,2)
     AuxArray(105,IPassQ) = AuxArray(105,IPassQ) + Xpb*TMPAuxArray(77) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     AuxArray(106,IPassQ) = AuxArray(106,IPassQ) + Xpb*TMPAuxArray(78) + alphaXpq*TmpArray7(78,2)
     AuxArray(107,IPassQ) = AuxArray(107,IPassQ) + Xpb*TMPAuxArray(79) + alphaXpq*TmpArray7(79,2)
     AuxArray(108,IPassQ) = AuxArray(108,IPassQ) + Xpb*TMPAuxArray(80) + alphaXpq*TmpArray7(80,2)
     AuxArray(109,IPassQ) = AuxArray(109,IPassQ) + Xpb*TMPAuxArray(81) + alphaXpq*TmpArray7(81,2)
     AuxArray(110,IPassQ) = AuxArray(110,IPassQ) + Xpb*TMPAuxArray(82) + alphaXpq*TmpArray7(82,2)
     AuxArray(111,IPassQ) = AuxArray(111,IPassQ) + Xpb*TMPAuxArray(83) + alphaXpq*TmpArray7(83,2)
     AuxArray(112,IPassQ) = AuxArray(112,IPassQ) + Xpb*TMPAuxArray(84) + alphaXpq*TmpArray7(84,2)
     AuxArray(113,IPassQ) = AuxArray(113,IPassQ) + Ypb*TMPAuxArray(78) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     AuxArray(114,IPassQ) = AuxArray(114,IPassQ) + Zpb*TMPAuxArray(78) + alphaZpq*TmpArray7(78,2)
     AuxArray(115,IPassQ) = AuxArray(115,IPassQ) + Ypb*TMPAuxArray(80) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     AuxArray(116,IPassQ) = AuxArray(116,IPassQ) + Ypb*TMPAuxArray(81) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     AuxArray(117,IPassQ) = AuxArray(117,IPassQ) + Ypb*TMPAuxArray(82) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     AuxArray(118,IPassQ) = AuxArray(118,IPassQ) + Ypb*TMPAuxArray(83) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     AuxArray(119,IPassQ) = AuxArray(119,IPassQ) + Ypb*TMPAuxArray(84) + alphaYpq*TmpArray7(84,2)
     AuxArray(120,IPassQ) = AuxArray(120,IPassQ) + Zpb*TMPAuxArray(84) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg8B(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:11,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Bcenter(3)
  real(realk),intent(inout) :: AUXarray(  165,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(  120)
  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP,RJ000(0: 8)
  real(realk) :: TwoTerms(  28)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D2=2.0E0_realk
  REAL(REALK),PARAMETER :: D2JP36=  5.2000000000000000E+01_realk
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
  real(realk) :: TMParray1(  1:  1,2:9)
  real(realk) :: TMParray2(  2:  4,2:8)
  real(realk) :: TMParray3(  5: 10,2:7)
  real(realk) :: TMParray4( 11: 20,2:6)
  real(realk) :: TMParray5( 21: 35,2:5)
  real(realk) :: TMParray6( 36: 56,2:4)
  real(realk) :: TMParray7( 57: 84,2:3)
  real(realk) :: TMParray8( 85:120,2:2)
  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Bx = -Bcenter(1)
  By = -Bcenter(2)
  Bz = -Bcenter(3)
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,  165
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpb = Pcent(1,iPrimP) + Bx
    Ypb = Pcent(2,iPrimP) + By
    Zpb = Pcent(3,iPrimP) + Bz
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000( 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000( 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000( 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000( 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000( 4) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000( 5) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
      RJ000( 6) = TABFJW( 6,IPNT)-TABFJW( 7,IPNT)*WDIFF+TABFJW( 8,IPNT)*W2+TABFJW( 9,IPNT)*W3
      RJ000( 7) = TABFJW( 7,IPNT)-TABFJW( 8,IPNT)*WDIFF+TABFJW( 9,IPNT)*W2+TABFJW(10,IPNT)*W3
      RJ000( 8) = TABFJW( 8,IPNT)-TABFJW( 9,IPNT)*WDIFF+TABFJW(10,IPNT)*W2+TABFJW(11,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000( 1) = RWVAL*(( 1 - D05)*RJ000( 0)-REXPW)
      RJ000( 2) = RWVAL*(( 2 - D05)*RJ000( 1)-REXPW)
      RJ000( 3) = RWVAL*(( 3 - D05)*RJ000( 2)-REXPW)
      RJ000( 4) = RWVAL*(( 4 - D05)*RJ000( 3)-REXPW)
      RJ000( 5) = RWVAL*(( 5 - D05)*RJ000( 4)-REXPW)
      RJ000( 6) = RWVAL*(( 6 - D05)*RJ000( 5)-REXPW)
      RJ000( 7) = RWVAL*(( 7 - D05)*RJ000( 6)-REXPW)
      RJ000( 8) = RWVAL*(( 8 - D05)*RJ000( 7)-REXPW)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000( 3) = RWVAL*( 3 - D05)*RJ000( 2)
      RJ000( 4) = RWVAL*( 4 - D05)*RJ000( 3)
      RJ000( 5) = RWVAL*( 5 - D05)*RJ000( 4)
      RJ000( 6) = RWVAL*( 6 - D05)*RJ000( 5)
      RJ000( 7) = RWVAL*( 7 - D05)*RJ000( 6)
      RJ000( 8) = RWVAL*( 8 - D05)*RJ000( 7)
     ENDIF
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac
     TMPAuxarray(1) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
     TMParray1(1, 5) = PREF*RJ000( 4)
     TMParray1(1, 6) = PREF*RJ000( 5)
     TMParray1(1, 7) = PREF*RJ000( 6)
     TMParray1(1, 8) = PREF*RJ000( 7)
     TMParray1(1, 9) = PREF*RJ000( 8)
     TMPAuxArray(2) = Xpb*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypb*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpb*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
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
     tmpArray2(2,7) = Xpb*tmpArray1(1,7) + alphaXpq*TmpArray1(1,8)
     tmpArray2(3,7) = Ypb*tmpArray1(1,7) + alphaYpq*TmpArray1(1,8)
     tmpArray2(4,7) = Zpb*tmpArray1(1,7) + alphaZpq*TmpArray1(1,8)
     tmpArray2(2,8) = Xpb*tmpArray1(1,8) + alphaXpq*TmpArray1(1,9)
     tmpArray2(3,8) = Ypb*tmpArray1(1,8) + alphaYpq*TmpArray1(1,9)
     tmpArray2(4,8) = Zpb*tmpArray1(1,8) + alphaZpq*TmpArray1(1,9)
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
     TwoTerms(1) = inv2expP*(TmpArray1(1,6) + alphaP*TmpArray1(1,7))
     tmpArray3(5,6) = Xpb*tmpArray2(2,6) + alphaXpq*TmpArray2(2,7) + TwoTerms(1)
     tmpArray3(6,6) = Xpb*tmpArray2(3,6) + alphaXpq*TmpArray2(3,7)
     tmpArray3(7,6) = Xpb*tmpArray2(4,6) + alphaXpq*TmpArray2(4,7)
     tmpArray3(8,6) = Ypb*tmpArray2(3,6) + alphaYpq*TmpArray2(3,7) + TwoTerms(1)
     tmpArray3(9,6) = Ypb*tmpArray2(4,6) + alphaYpq*TmpArray2(4,7)
     tmpArray3(10,6) = Zpb*tmpArray2(4,6) + alphaZpq*TmpArray2(4,7) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,7) + alphaP*TmpArray1(1,8))
     tmpArray3(5,7) = Xpb*tmpArray2(2,7) + alphaXpq*TmpArray2(2,8) + TwoTerms(1)
     tmpArray3(6,7) = Xpb*tmpArray2(3,7) + alphaXpq*TmpArray2(3,8)
     tmpArray3(7,7) = Xpb*tmpArray2(4,7) + alphaXpq*TmpArray2(4,8)
     tmpArray3(8,7) = Ypb*tmpArray2(3,7) + alphaYpq*TmpArray2(3,8) + TwoTerms(1)
     tmpArray3(9,7) = Ypb*tmpArray2(4,7) + alphaYpq*TmpArray2(4,8)
     tmpArray3(10,7) = Zpb*tmpArray2(4,7) + alphaZpq*TmpArray2(4,8) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TMPAuxArray(2) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(3) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(4) + alphaP*TmpArray2(4,2))
     TMPAuxArray(11) = Xpb*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Ypb*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zpb*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xpb*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xpb*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xpb*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Ypb*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zpb*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Ypb*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zpb*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(TmpArray2(2,5) + alphaP*TmpArray2(2,6))
     TwoTerms(2) = inv2expP*(TmpArray2(3,5) + alphaP*TmpArray2(3,6))
     TwoTerms(3) = inv2expP*(TmpArray2(4,5) + alphaP*TmpArray2(4,6))
     tmpArray4(11,5) = Xpb*tmpArray3(5,5) + alphaXpq*TmpArray3(5,6) + 2*TwoTerms(1)
     tmpArray4(12,5) = Ypb*tmpArray3(5,5) + alphaYpq*TmpArray3(5,6)
     tmpArray4(13,5) = Zpb*tmpArray3(5,5) + alphaZpq*TmpArray3(5,6)
     tmpArray4(14,5) = Xpb*tmpArray3(8,5) + alphaXpq*TmpArray3(8,6)
     tmpArray4(15,5) = Xpb*tmpArray3(9,5) + alphaXpq*TmpArray3(9,6)
     tmpArray4(16,5) = Xpb*tmpArray3(10,5) + alphaXpq*TmpArray3(10,6)
     tmpArray4(17,5) = Ypb*tmpArray3(8,5) + alphaYpq*TmpArray3(8,6) + 2*TwoTerms(2)
     tmpArray4(18,5) = Zpb*tmpArray3(8,5) + alphaZpq*TmpArray3(8,6)
     tmpArray4(19,5) = Ypb*tmpArray3(10,5) + alphaYpq*TmpArray3(10,6)
     tmpArray4(20,5) = Zpb*tmpArray3(10,5) + alphaZpq*TmpArray3(10,6) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,6) + alphaP*TmpArray2(2,7))
     TwoTerms(2) = inv2expP*(TmpArray2(3,6) + alphaP*TmpArray2(3,7))
     TwoTerms(3) = inv2expP*(TmpArray2(4,6) + alphaP*TmpArray2(4,7))
     tmpArray4(11,6) = Xpb*tmpArray3(5,6) + alphaXpq*TmpArray3(5,7) + 2*TwoTerms(1)
     tmpArray4(12,6) = Ypb*tmpArray3(5,6) + alphaYpq*TmpArray3(5,7)
     tmpArray4(13,6) = Zpb*tmpArray3(5,6) + alphaZpq*TmpArray3(5,7)
     tmpArray4(14,6) = Xpb*tmpArray3(8,6) + alphaXpq*TmpArray3(8,7)
     tmpArray4(15,6) = Xpb*tmpArray3(9,6) + alphaXpq*TmpArray3(9,7)
     tmpArray4(16,6) = Xpb*tmpArray3(10,6) + alphaXpq*TmpArray3(10,7)
     tmpArray4(17,6) = Ypb*tmpArray3(8,6) + alphaYpq*TmpArray3(8,7) + 2*TwoTerms(2)
     tmpArray4(18,6) = Zpb*tmpArray3(8,6) + alphaZpq*TmpArray3(8,7)
     tmpArray4(19,6) = Ypb*tmpArray3(10,6) + alphaYpq*TmpArray3(10,7)
     tmpArray4(20,6) = Zpb*tmpArray3(10,6) + alphaZpq*TmpArray3(10,7) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(10) + alphaP*TmpArray3(10,2))
     TMPAuxArray(21) = Xpb*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Ypb*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zpb*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xpb*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Ypb*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xpb*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xpb*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xpb*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xpb*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xpb*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Ypb*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zpb*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Ypb*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Ypb*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zpb*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(TmpArray3(5,4) + alphaP*TmpArray3(5,5))
     TwoTerms(2) = inv2expP*(TmpArray3(8,4) + alphaP*TmpArray3(8,5))
     TwoTerms(3) = inv2expP*(TmpArray3(10,4) + alphaP*TmpArray3(10,5))
     tmpArray5(21,4) = Xpb*tmpArray4(11,4) + alphaXpq*TmpArray4(11,5) + 3*TwoTerms(1)
     tmpArray5(22,4) = Ypb*tmpArray4(11,4) + alphaYpq*TmpArray4(11,5)
     tmpArray5(23,4) = Zpb*tmpArray4(11,4) + alphaZpq*TmpArray4(11,5)
     tmpArray5(24,4) = Xpb*tmpArray4(14,4) + alphaXpq*TmpArray4(14,5) + TwoTerms(2)
     tmpArray5(25,4) = Ypb*tmpArray4(13,4) + alphaYpq*TmpArray4(13,5)
     tmpArray5(26,4) = Xpb*tmpArray4(16,4) + alphaXpq*TmpArray4(16,5) + TwoTerms(3)
     tmpArray5(27,4) = Xpb*tmpArray4(17,4) + alphaXpq*TmpArray4(17,5)
     tmpArray5(28,4) = Xpb*tmpArray4(18,4) + alphaXpq*TmpArray4(18,5)
     tmpArray5(29,4) = Xpb*tmpArray4(19,4) + alphaXpq*TmpArray4(19,5)
     tmpArray5(30,4) = Xpb*tmpArray4(20,4) + alphaXpq*TmpArray4(20,5)
     tmpArray5(31,4) = Ypb*tmpArray4(17,4) + alphaYpq*TmpArray4(17,5) + 3*TwoTerms(2)
     tmpArray5(32,4) = Zpb*tmpArray4(17,4) + alphaZpq*TmpArray4(17,5)
     tmpArray5(33,4) = Ypb*tmpArray4(19,4) + alphaYpq*TmpArray4(19,5) + TwoTerms(3)
     tmpArray5(34,4) = Ypb*tmpArray4(20,4) + alphaYpq*TmpArray4(20,5)
     tmpArray5(35,4) = Zpb*tmpArray4(20,4) + alphaZpq*TmpArray4(20,5) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,5) + alphaP*TmpArray3(5,6))
     TwoTerms(2) = inv2expP*(TmpArray3(8,5) + alphaP*TmpArray3(8,6))
     TwoTerms(3) = inv2expP*(TmpArray3(10,5) + alphaP*TmpArray3(10,6))
     tmpArray5(21,5) = Xpb*tmpArray4(11,5) + alphaXpq*TmpArray4(11,6) + 3*TwoTerms(1)
     tmpArray5(22,5) = Ypb*tmpArray4(11,5) + alphaYpq*TmpArray4(11,6)
     tmpArray5(23,5) = Zpb*tmpArray4(11,5) + alphaZpq*TmpArray4(11,6)
     tmpArray5(24,5) = Xpb*tmpArray4(14,5) + alphaXpq*TmpArray4(14,6) + TwoTerms(2)
     tmpArray5(25,5) = Ypb*tmpArray4(13,5) + alphaYpq*TmpArray4(13,6)
     tmpArray5(26,5) = Xpb*tmpArray4(16,5) + alphaXpq*TmpArray4(16,6) + TwoTerms(3)
     tmpArray5(27,5) = Xpb*tmpArray4(17,5) + alphaXpq*TmpArray4(17,6)
     tmpArray5(28,5) = Xpb*tmpArray4(18,5) + alphaXpq*TmpArray4(18,6)
     tmpArray5(29,5) = Xpb*tmpArray4(19,5) + alphaXpq*TmpArray4(19,6)
     tmpArray5(30,5) = Xpb*tmpArray4(20,5) + alphaXpq*TmpArray4(20,6)
     tmpArray5(31,5) = Ypb*tmpArray4(17,5) + alphaYpq*TmpArray4(17,6) + 3*TwoTerms(2)
     tmpArray5(32,5) = Zpb*tmpArray4(17,5) + alphaZpq*TmpArray4(17,6)
     tmpArray5(33,5) = Ypb*tmpArray4(19,5) + alphaYpq*TmpArray4(19,6) + TwoTerms(3)
     tmpArray5(34,5) = Ypb*tmpArray4(20,5) + alphaYpq*TmpArray4(20,6)
     tmpArray5(35,5) = Zpb*tmpArray4(20,5) + alphaZpq*TmpArray4(20,6) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(11) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(14) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(16) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(17) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(19) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(20) + alphaP*TmpArray4(20,2))
     TMPAuxArray(36) = Xpb*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     TMPAuxArray(37) = Ypb*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     TMPAuxArray(38) = Zpb*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     TMPAuxArray(39) = Xpb*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     TMPAuxArray(40) = Ypb*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     TMPAuxArray(41) = Xpb*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     TMPAuxArray(42) = Xpb*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     TMPAuxArray(43) = Zpb*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     TMPAuxArray(44) = Ypb*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     TMPAuxArray(45) = Xpb*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     TMPAuxArray(46) = Xpb*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     TMPAuxArray(47) = Xpb*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     TMPAuxArray(48) = Xpb*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     TMPAuxArray(49) = Xpb*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     TMPAuxArray(50) = Xpb*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     TMPAuxArray(51) = Ypb*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     TMPAuxArray(52) = Zpb*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     TMPAuxArray(53) = Ypb*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     TMPAuxArray(54) = Ypb*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     TMPAuxArray(55) = Ypb*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     TMPAuxArray(56) = Zpb*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     TwoTerms(1) = inv2expP*(TmpArray4(11,3) + alphaP*TmpArray4(11,4))
     TwoTerms(2) = inv2expP*(TmpArray4(14,3) + alphaP*TmpArray4(14,4))
     TwoTerms(3) = inv2expP*(TmpArray4(16,3) + alphaP*TmpArray4(16,4))
     TwoTerms(4) = inv2expP*(TmpArray4(17,3) + alphaP*TmpArray4(17,4))
     TwoTerms(5) = inv2expP*(TmpArray4(19,3) + alphaP*TmpArray4(19,4))
     TwoTerms(6) = inv2expP*(TmpArray4(20,3) + alphaP*TmpArray4(20,4))
     tmpArray6(36,3) = Xpb*tmpArray5(21,3) + alphaXpq*TmpArray5(21,4) + 4*TwoTerms(1)
     tmpArray6(37,3) = Ypb*tmpArray5(21,3) + alphaYpq*TmpArray5(21,4)
     tmpArray6(38,3) = Zpb*tmpArray5(21,3) + alphaZpq*TmpArray5(21,4)
     tmpArray6(39,3) = Xpb*tmpArray5(24,3) + alphaXpq*TmpArray5(24,4) + 2*TwoTerms(2)
     tmpArray6(40,3) = Ypb*tmpArray5(23,3) + alphaYpq*TmpArray5(23,4)
     tmpArray6(41,3) = Xpb*tmpArray5(26,3) + alphaXpq*TmpArray5(26,4) + 2*TwoTerms(3)
     tmpArray6(42,3) = Xpb*tmpArray5(27,3) + alphaXpq*TmpArray5(27,4) + TwoTerms(4)
     tmpArray6(43,3) = Zpb*tmpArray5(24,3) + alphaZpq*TmpArray5(24,4)
     tmpArray6(44,3) = Ypb*tmpArray5(26,3) + alphaYpq*TmpArray5(26,4)
     tmpArray6(45,3) = Xpb*tmpArray5(30,3) + alphaXpq*TmpArray5(30,4) + TwoTerms(6)
     tmpArray6(46,3) = Xpb*tmpArray5(31,3) + alphaXpq*TmpArray5(31,4)
     tmpArray6(47,3) = Xpb*tmpArray5(32,3) + alphaXpq*TmpArray5(32,4)
     tmpArray6(48,3) = Xpb*tmpArray5(33,3) + alphaXpq*TmpArray5(33,4)
     tmpArray6(49,3) = Xpb*tmpArray5(34,3) + alphaXpq*TmpArray5(34,4)
     tmpArray6(50,3) = Xpb*tmpArray5(35,3) + alphaXpq*TmpArray5(35,4)
     tmpArray6(51,3) = Ypb*tmpArray5(31,3) + alphaYpq*TmpArray5(31,4) + 4*TwoTerms(4)
     tmpArray6(52,3) = Zpb*tmpArray5(31,3) + alphaZpq*TmpArray5(31,4)
     tmpArray6(53,3) = Ypb*tmpArray5(33,3) + alphaYpq*TmpArray5(33,4) + 2*TwoTerms(5)
     tmpArray6(54,3) = Ypb*tmpArray5(34,3) + alphaYpq*TmpArray5(34,4) + TwoTerms(6)
     tmpArray6(55,3) = Ypb*tmpArray5(35,3) + alphaYpq*TmpArray5(35,4)
     tmpArray6(56,3) = Zpb*tmpArray5(35,3) + alphaZpq*TmpArray5(35,4) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TmpArray4(11,4) + alphaP*TmpArray4(11,5))
     TwoTerms(2) = inv2expP*(TmpArray4(14,4) + alphaP*TmpArray4(14,5))
     TwoTerms(3) = inv2expP*(TmpArray4(16,4) + alphaP*TmpArray4(16,5))
     TwoTerms(4) = inv2expP*(TmpArray4(17,4) + alphaP*TmpArray4(17,5))
     TwoTerms(5) = inv2expP*(TmpArray4(19,4) + alphaP*TmpArray4(19,5))
     TwoTerms(6) = inv2expP*(TmpArray4(20,4) + alphaP*TmpArray4(20,5))
     tmpArray6(36,4) = Xpb*tmpArray5(21,4) + alphaXpq*TmpArray5(21,5) + 4*TwoTerms(1)
     tmpArray6(37,4) = Ypb*tmpArray5(21,4) + alphaYpq*TmpArray5(21,5)
     tmpArray6(38,4) = Zpb*tmpArray5(21,4) + alphaZpq*TmpArray5(21,5)
     tmpArray6(39,4) = Xpb*tmpArray5(24,4) + alphaXpq*TmpArray5(24,5) + 2*TwoTerms(2)
     tmpArray6(40,4) = Ypb*tmpArray5(23,4) + alphaYpq*TmpArray5(23,5)
     tmpArray6(41,4) = Xpb*tmpArray5(26,4) + alphaXpq*TmpArray5(26,5) + 2*TwoTerms(3)
     tmpArray6(42,4) = Xpb*tmpArray5(27,4) + alphaXpq*TmpArray5(27,5) + TwoTerms(4)
     tmpArray6(43,4) = Zpb*tmpArray5(24,4) + alphaZpq*TmpArray5(24,5)
     tmpArray6(44,4) = Ypb*tmpArray5(26,4) + alphaYpq*TmpArray5(26,5)
     tmpArray6(45,4) = Xpb*tmpArray5(30,4) + alphaXpq*TmpArray5(30,5) + TwoTerms(6)
     tmpArray6(46,4) = Xpb*tmpArray5(31,4) + alphaXpq*TmpArray5(31,5)
     tmpArray6(47,4) = Xpb*tmpArray5(32,4) + alphaXpq*TmpArray5(32,5)
     tmpArray6(48,4) = Xpb*tmpArray5(33,4) + alphaXpq*TmpArray5(33,5)
     tmpArray6(49,4) = Xpb*tmpArray5(34,4) + alphaXpq*TmpArray5(34,5)
     tmpArray6(50,4) = Xpb*tmpArray5(35,4) + alphaXpq*TmpArray5(35,5)
     tmpArray6(51,4) = Ypb*tmpArray5(31,4) + alphaYpq*TmpArray5(31,5) + 4*TwoTerms(4)
     tmpArray6(52,4) = Zpb*tmpArray5(31,4) + alphaZpq*TmpArray5(31,5)
     tmpArray6(53,4) = Ypb*tmpArray5(33,4) + alphaYpq*TmpArray5(33,5) + 2*TwoTerms(5)
     tmpArray6(54,4) = Ypb*tmpArray5(34,4) + alphaYpq*TmpArray5(34,5) + TwoTerms(6)
     tmpArray6(55,4) = Ypb*tmpArray5(35,4) + alphaYpq*TmpArray5(35,5)
     tmpArray6(56,4) = Zpb*tmpArray5(35,4) + alphaZpq*TmpArray5(35,5) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TMPAuxArray(21) + alphaP*TmpArray5(21,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(24) + alphaP*TmpArray5(24,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(26) + alphaP*TmpArray5(26,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(27) + alphaP*TmpArray5(27,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(30) + alphaP*TmpArray5(30,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(31) + alphaP*TmpArray5(31,2))
     TwoTerms(7) = inv2expP*(TMPAuxArray(33) + alphaP*TmpArray5(33,2))
     TwoTerms(8) = inv2expP*(TMPAuxArray(34) + alphaP*TmpArray5(34,2))
     TwoTerms(9) = inv2expP*(TMPAuxArray(35) + alphaP*TmpArray5(35,2))
     TMPAuxArray(57) = Xpb*TMPAuxArray(36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     TMPAuxArray(58) = Ypb*TMPAuxArray(36) + alphaYpq*TmpArray6(36,2)
     TMPAuxArray(59) = Zpb*TMPAuxArray(36) + alphaZpq*TmpArray6(36,2)
     TMPAuxArray(60) = Xpb*TMPAuxArray(39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     TMPAuxArray(61) = Ypb*TMPAuxArray(38) + alphaYpq*TmpArray6(38,2)
     TMPAuxArray(62) = Xpb*TMPAuxArray(41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     TMPAuxArray(63) = Xpb*TMPAuxArray(42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     TMPAuxArray(64) = Zpb*TMPAuxArray(39) + alphaZpq*TmpArray6(39,2)
     TMPAuxArray(65) = Ypb*TMPAuxArray(41) + alphaYpq*TmpArray6(41,2)
     TMPAuxArray(66) = Xpb*TMPAuxArray(45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     TMPAuxArray(67) = Xpb*TMPAuxArray(46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     TMPAuxArray(68) = Zpb*TMPAuxArray(42) + alphaZpq*TmpArray6(42,2)
     TMPAuxArray(69) = Xpb*TMPAuxArray(48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     TMPAuxArray(70) = Ypb*TMPAuxArray(45) + alphaYpq*TmpArray6(45,2)
     TMPAuxArray(71) = Xpb*TMPAuxArray(50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     TMPAuxArray(72) = Xpb*TMPAuxArray(51) + alphaXpq*TmpArray6(51,2)
     TMPAuxArray(73) = Xpb*TMPAuxArray(52) + alphaXpq*TmpArray6(52,2)
     TMPAuxArray(74) = Xpb*TMPAuxArray(53) + alphaXpq*TmpArray6(53,2)
     TMPAuxArray(75) = Xpb*TMPAuxArray(54) + alphaXpq*TmpArray6(54,2)
     TMPAuxArray(76) = Xpb*TMPAuxArray(55) + alphaXpq*TmpArray6(55,2)
     TMPAuxArray(77) = Xpb*TMPAuxArray(56) + alphaXpq*TmpArray6(56,2)
     TMPAuxArray(78) = Ypb*TMPAuxArray(51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     TMPAuxArray(79) = Zpb*TMPAuxArray(51) + alphaZpq*TmpArray6(51,2)
     TMPAuxArray(80) = Ypb*TMPAuxArray(53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     TMPAuxArray(81) = Ypb*TMPAuxArray(54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     TMPAuxArray(82) = Ypb*TMPAuxArray(55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     TMPAuxArray(83) = Ypb*TMPAuxArray(56) + alphaYpq*TmpArray6(56,2)
     TMPAuxArray(84) = Zpb*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expP*(TmpArray5(21,2) + alphaP*TmpArray5(21,3))
     TwoTerms(2) = inv2expP*(TmpArray5(24,2) + alphaP*TmpArray5(24,3))
     TwoTerms(3) = inv2expP*(TmpArray5(26,2) + alphaP*TmpArray5(26,3))
     TwoTerms(4) = inv2expP*(TmpArray5(27,2) + alphaP*TmpArray5(27,3))
     TwoTerms(5) = inv2expP*(TmpArray5(30,2) + alphaP*TmpArray5(30,3))
     TwoTerms(6) = inv2expP*(TmpArray5(31,2) + alphaP*TmpArray5(31,3))
     TwoTerms(7) = inv2expP*(TmpArray5(33,2) + alphaP*TmpArray5(33,3))
     TwoTerms(8) = inv2expP*(TmpArray5(34,2) + alphaP*TmpArray5(34,3))
     TwoTerms(9) = inv2expP*(TmpArray5(35,2) + alphaP*TmpArray5(35,3))
     tmpArray7(57,2) = Xpb*tmpArray6(36,2) + alphaXpq*TmpArray6(36,3) + 5*TwoTerms(1)
     tmpArray7(58,2) = Ypb*tmpArray6(36,2) + alphaYpq*TmpArray6(36,3)
     tmpArray7(59,2) = Zpb*tmpArray6(36,2) + alphaZpq*TmpArray6(36,3)
     tmpArray7(60,2) = Xpb*tmpArray6(39,2) + alphaXpq*TmpArray6(39,3) + 3*TwoTerms(2)
     tmpArray7(61,2) = Ypb*tmpArray6(38,2) + alphaYpq*TmpArray6(38,3)
     tmpArray7(62,2) = Xpb*tmpArray6(41,2) + alphaXpq*TmpArray6(41,3) + 3*TwoTerms(3)
     tmpArray7(63,2) = Xpb*tmpArray6(42,2) + alphaXpq*TmpArray6(42,3) + 2*TwoTerms(4)
     tmpArray7(64,2) = Zpb*tmpArray6(39,2) + alphaZpq*TmpArray6(39,3)
     tmpArray7(65,2) = Ypb*tmpArray6(41,2) + alphaYpq*TmpArray6(41,3)
     tmpArray7(66,2) = Xpb*tmpArray6(45,2) + alphaXpq*TmpArray6(45,3) + 2*TwoTerms(5)
     tmpArray7(67,2) = Xpb*tmpArray6(46,2) + alphaXpq*TmpArray6(46,3) + TwoTerms(6)
     tmpArray7(68,2) = Zpb*tmpArray6(42,2) + alphaZpq*TmpArray6(42,3)
     tmpArray7(69,2) = Xpb*tmpArray6(48,2) + alphaXpq*TmpArray6(48,3) + TwoTerms(7)
     tmpArray7(70,2) = Ypb*tmpArray6(45,2) + alphaYpq*TmpArray6(45,3)
     tmpArray7(71,2) = Xpb*tmpArray6(50,2) + alphaXpq*TmpArray6(50,3) + TwoTerms(9)
     tmpArray7(72,2) = Xpb*tmpArray6(51,2) + alphaXpq*TmpArray6(51,3)
     tmpArray7(73,2) = Xpb*tmpArray6(52,2) + alphaXpq*TmpArray6(52,3)
     tmpArray7(74,2) = Xpb*tmpArray6(53,2) + alphaXpq*TmpArray6(53,3)
     tmpArray7(75,2) = Xpb*tmpArray6(54,2) + alphaXpq*TmpArray6(54,3)
     tmpArray7(76,2) = Xpb*tmpArray6(55,2) + alphaXpq*TmpArray6(55,3)
     tmpArray7(77,2) = Xpb*tmpArray6(56,2) + alphaXpq*TmpArray6(56,3)
     tmpArray7(78,2) = Ypb*tmpArray6(51,2) + alphaYpq*TmpArray6(51,3) + 5*TwoTerms(6)
     tmpArray7(79,2) = Zpb*tmpArray6(51,2) + alphaZpq*TmpArray6(51,3)
     tmpArray7(80,2) = Ypb*tmpArray6(53,2) + alphaYpq*TmpArray6(53,3) + 3*TwoTerms(7)
     tmpArray7(81,2) = Ypb*tmpArray6(54,2) + alphaYpq*TmpArray6(54,3) + 2*TwoTerms(8)
     tmpArray7(82,2) = Ypb*tmpArray6(55,2) + alphaYpq*TmpArray6(55,3) + TwoTerms(9)
     tmpArray7(83,2) = Ypb*tmpArray6(56,2) + alphaYpq*TmpArray6(56,3)
     tmpArray7(84,2) = Zpb*tmpArray6(56,2) + alphaZpq*TmpArray6(56,3) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expP*(TmpArray5(21,3) + alphaP*TmpArray5(21,4))
     TwoTerms(2) = inv2expP*(TmpArray5(24,3) + alphaP*TmpArray5(24,4))
     TwoTerms(3) = inv2expP*(TmpArray5(26,3) + alphaP*TmpArray5(26,4))
     TwoTerms(4) = inv2expP*(TmpArray5(27,3) + alphaP*TmpArray5(27,4))
     TwoTerms(5) = inv2expP*(TmpArray5(30,3) + alphaP*TmpArray5(30,4))
     TwoTerms(6) = inv2expP*(TmpArray5(31,3) + alphaP*TmpArray5(31,4))
     TwoTerms(7) = inv2expP*(TmpArray5(33,3) + alphaP*TmpArray5(33,4))
     TwoTerms(8) = inv2expP*(TmpArray5(34,3) + alphaP*TmpArray5(34,4))
     TwoTerms(9) = inv2expP*(TmpArray5(35,3) + alphaP*TmpArray5(35,4))
     tmpArray7(57,3) = Xpb*tmpArray6(36,3) + alphaXpq*TmpArray6(36,4) + 5*TwoTerms(1)
     tmpArray7(58,3) = Ypb*tmpArray6(36,3) + alphaYpq*TmpArray6(36,4)
     tmpArray7(59,3) = Zpb*tmpArray6(36,3) + alphaZpq*TmpArray6(36,4)
     tmpArray7(60,3) = Xpb*tmpArray6(39,3) + alphaXpq*TmpArray6(39,4) + 3*TwoTerms(2)
     tmpArray7(61,3) = Ypb*tmpArray6(38,3) + alphaYpq*TmpArray6(38,4)
     tmpArray7(62,3) = Xpb*tmpArray6(41,3) + alphaXpq*TmpArray6(41,4) + 3*TwoTerms(3)
     tmpArray7(63,3) = Xpb*tmpArray6(42,3) + alphaXpq*TmpArray6(42,4) + 2*TwoTerms(4)
     tmpArray7(64,3) = Zpb*tmpArray6(39,3) + alphaZpq*TmpArray6(39,4)
     tmpArray7(65,3) = Ypb*tmpArray6(41,3) + alphaYpq*TmpArray6(41,4)
     tmpArray7(66,3) = Xpb*tmpArray6(45,3) + alphaXpq*TmpArray6(45,4) + 2*TwoTerms(5)
     tmpArray7(67,3) = Xpb*tmpArray6(46,3) + alphaXpq*TmpArray6(46,4) + TwoTerms(6)
     tmpArray7(68,3) = Zpb*tmpArray6(42,3) + alphaZpq*TmpArray6(42,4)
     tmpArray7(69,3) = Xpb*tmpArray6(48,3) + alphaXpq*TmpArray6(48,4) + TwoTerms(7)
     tmpArray7(70,3) = Ypb*tmpArray6(45,3) + alphaYpq*TmpArray6(45,4)
     tmpArray7(71,3) = Xpb*tmpArray6(50,3) + alphaXpq*TmpArray6(50,4) + TwoTerms(9)
     tmpArray7(72,3) = Xpb*tmpArray6(51,3) + alphaXpq*TmpArray6(51,4)
     tmpArray7(73,3) = Xpb*tmpArray6(52,3) + alphaXpq*TmpArray6(52,4)
     tmpArray7(74,3) = Xpb*tmpArray6(53,3) + alphaXpq*TmpArray6(53,4)
     tmpArray7(75,3) = Xpb*tmpArray6(54,3) + alphaXpq*TmpArray6(54,4)
     tmpArray7(76,3) = Xpb*tmpArray6(55,3) + alphaXpq*TmpArray6(55,4)
     tmpArray7(77,3) = Xpb*tmpArray6(56,3) + alphaXpq*TmpArray6(56,4)
     tmpArray7(78,3) = Ypb*tmpArray6(51,3) + alphaYpq*TmpArray6(51,4) + 5*TwoTerms(6)
     tmpArray7(79,3) = Zpb*tmpArray6(51,3) + alphaZpq*TmpArray6(51,4)
     tmpArray7(80,3) = Ypb*tmpArray6(53,3) + alphaYpq*TmpArray6(53,4) + 3*TwoTerms(7)
     tmpArray7(81,3) = Ypb*tmpArray6(54,3) + alphaYpq*TmpArray6(54,4) + 2*TwoTerms(8)
     tmpArray7(82,3) = Ypb*tmpArray6(55,3) + alphaYpq*TmpArray6(55,4) + TwoTerms(9)
     tmpArray7(83,3) = Ypb*tmpArray6(56,3) + alphaYpq*TmpArray6(56,4)
     tmpArray7(84,3) = Zpb*tmpArray6(56,3) + alphaZpq*TmpArray6(56,4) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expP*(TMPAuxArray(36) + alphaP*TmpArray6(36,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(39) + alphaP*TmpArray6(39,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(41) + alphaP*TmpArray6(41,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(42) + alphaP*TmpArray6(42,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(45) + alphaP*TmpArray6(45,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(46) + alphaP*TmpArray6(46,2))
     TwoTerms(7) = inv2expP*(TMPAuxArray(48) + alphaP*TmpArray6(48,2))
     TwoTerms(8) = inv2expP*(TMPAuxArray(50) + alphaP*TmpArray6(50,2))
     TwoTerms(9) = inv2expP*(TMPAuxArray(51) + alphaP*TmpArray6(51,2))
     TwoTerms(10) = inv2expP*(TMPAuxArray(53) + alphaP*TmpArray6(53,2))
     TwoTerms(11) = inv2expP*(TMPAuxArray(54) + alphaP*TmpArray6(54,2))
     TwoTerms(12) = inv2expP*(TMPAuxArray(55) + alphaP*TmpArray6(55,2))
     TwoTerms(13) = inv2expP*(TMPAuxArray(56) + alphaP*TmpArray6(56,2))
     TMPAuxArray(85) = Xpb*TMPAuxArray(57) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     TMPAuxArray(86) = Ypb*TMPAuxArray(57) + alphaYpq*TmpArray7(57,2)
     TMPAuxArray(87) = Zpb*TMPAuxArray(57) + alphaZpq*TmpArray7(57,2)
     TMPAuxArray(88) = Xpb*TMPAuxArray(60) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     TMPAuxArray(89) = Ypb*TMPAuxArray(59) + alphaYpq*TmpArray7(59,2)
     TMPAuxArray(90) = Xpb*TMPAuxArray(62) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     TMPAuxArray(91) = Xpb*TMPAuxArray(63) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     TMPAuxArray(92) = Zpb*TMPAuxArray(60) + alphaZpq*TmpArray7(60,2)
     TMPAuxArray(93) = Ypb*TMPAuxArray(62) + alphaYpq*TmpArray7(62,2)
     TMPAuxArray(94) = Xpb*TMPAuxArray(66) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     TMPAuxArray(95) = Xpb*TMPAuxArray(67) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     TMPAuxArray(96) = Zpb*TMPAuxArray(63) + alphaZpq*TmpArray7(63,2)
     TMPAuxArray(97) = Xpb*TMPAuxArray(69) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     TMPAuxArray(98) = Ypb*TMPAuxArray(66) + alphaYpq*TmpArray7(66,2)
     TMPAuxArray(99) = Xpb*TMPAuxArray(71) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     TMPAuxArray(100) = Xpb*TMPAuxArray(72) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     TMPAuxArray(101) = Zpb*TMPAuxArray(67) + alphaZpq*TmpArray7(67,2)
     TMPAuxArray(102) = Xpb*TMPAuxArray(74) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     TMPAuxArray(103) = Xpb*TMPAuxArray(75) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     TMPAuxArray(104) = Ypb*TMPAuxArray(71) + alphaYpq*TmpArray7(71,2)
     TMPAuxArray(105) = Xpb*TMPAuxArray(77) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     TMPAuxArray(106) = Xpb*TMPAuxArray(78) + alphaXpq*TmpArray7(78,2)
     TMPAuxArray(107) = Xpb*TMPAuxArray(79) + alphaXpq*TmpArray7(79,2)
     TMPAuxArray(108) = Xpb*TMPAuxArray(80) + alphaXpq*TmpArray7(80,2)
     TMPAuxArray(109) = Xpb*TMPAuxArray(81) + alphaXpq*TmpArray7(81,2)
     TMPAuxArray(110) = Xpb*TMPAuxArray(82) + alphaXpq*TmpArray7(82,2)
     TMPAuxArray(111) = Xpb*TMPAuxArray(83) + alphaXpq*TmpArray7(83,2)
     TMPAuxArray(112) = Xpb*TMPAuxArray(84) + alphaXpq*TmpArray7(84,2)
     TMPAuxArray(113) = Ypb*TMPAuxArray(78) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     TMPAuxArray(114) = Zpb*TMPAuxArray(78) + alphaZpq*TmpArray7(78,2)
     TMPAuxArray(115) = Ypb*TMPAuxArray(80) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     TMPAuxArray(116) = Ypb*TMPAuxArray(81) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     TMPAuxArray(117) = Ypb*TMPAuxArray(82) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     TMPAuxArray(118) = Ypb*TMPAuxArray(83) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     TMPAuxArray(119) = Ypb*TMPAuxArray(84) + alphaYpq*TmpArray7(84,2)
     TMPAuxArray(120) = Zpb*TMPAuxArray(84) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
     TwoTerms(1) = inv2expP*(TmpArray6(36,2) + alphaP*TmpArray6(36,3))
     TwoTerms(2) = inv2expP*(TmpArray6(39,2) + alphaP*TmpArray6(39,3))
     TwoTerms(3) = inv2expP*(TmpArray6(41,2) + alphaP*TmpArray6(41,3))
     TwoTerms(4) = inv2expP*(TmpArray6(42,2) + alphaP*TmpArray6(42,3))
     TwoTerms(5) = inv2expP*(TmpArray6(45,2) + alphaP*TmpArray6(45,3))
     TwoTerms(6) = inv2expP*(TmpArray6(46,2) + alphaP*TmpArray6(46,3))
     TwoTerms(7) = inv2expP*(TmpArray6(48,2) + alphaP*TmpArray6(48,3))
     TwoTerms(8) = inv2expP*(TmpArray6(50,2) + alphaP*TmpArray6(50,3))
     TwoTerms(9) = inv2expP*(TmpArray6(51,2) + alphaP*TmpArray6(51,3))
     TwoTerms(10) = inv2expP*(TmpArray6(53,2) + alphaP*TmpArray6(53,3))
     TwoTerms(11) = inv2expP*(TmpArray6(54,2) + alphaP*TmpArray6(54,3))
     TwoTerms(12) = inv2expP*(TmpArray6(55,2) + alphaP*TmpArray6(55,3))
     TwoTerms(13) = inv2expP*(TmpArray6(56,2) + alphaP*TmpArray6(56,3))
     tmpArray8(85,2) = Xpb*tmpArray7(57,2) + alphaXpq*TmpArray7(57,3) + 6*TwoTerms(1)
     tmpArray8(86,2) = Ypb*tmpArray7(57,2) + alphaYpq*TmpArray7(57,3)
     tmpArray8(87,2) = Zpb*tmpArray7(57,2) + alphaZpq*TmpArray7(57,3)
     tmpArray8(88,2) = Xpb*tmpArray7(60,2) + alphaXpq*TmpArray7(60,3) + 4*TwoTerms(2)
     tmpArray8(89,2) = Ypb*tmpArray7(59,2) + alphaYpq*TmpArray7(59,3)
     tmpArray8(90,2) = Xpb*tmpArray7(62,2) + alphaXpq*TmpArray7(62,3) + 4*TwoTerms(3)
     tmpArray8(91,2) = Xpb*tmpArray7(63,2) + alphaXpq*TmpArray7(63,3) + 3*TwoTerms(4)
     tmpArray8(92,2) = Zpb*tmpArray7(60,2) + alphaZpq*TmpArray7(60,3)
     tmpArray8(93,2) = Ypb*tmpArray7(62,2) + alphaYpq*TmpArray7(62,3)
     tmpArray8(94,2) = Xpb*tmpArray7(66,2) + alphaXpq*TmpArray7(66,3) + 3*TwoTerms(5)
     tmpArray8(95,2) = Xpb*tmpArray7(67,2) + alphaXpq*TmpArray7(67,3) + 2*TwoTerms(6)
     tmpArray8(96,2) = Zpb*tmpArray7(63,2) + alphaZpq*TmpArray7(63,3)
     tmpArray8(97,2) = Xpb*tmpArray7(69,2) + alphaXpq*TmpArray7(69,3) + 2*TwoTerms(7)
     tmpArray8(98,2) = Ypb*tmpArray7(66,2) + alphaYpq*TmpArray7(66,3)
     tmpArray8(99,2) = Xpb*tmpArray7(71,2) + alphaXpq*TmpArray7(71,3) + 2*TwoTerms(8)
     tmpArray8(100,2) = Xpb*tmpArray7(72,2) + alphaXpq*TmpArray7(72,3) + TwoTerms(9)
     tmpArray8(101,2) = Zpb*tmpArray7(67,2) + alphaZpq*TmpArray7(67,3)
     tmpArray8(102,2) = Xpb*tmpArray7(74,2) + alphaXpq*TmpArray7(74,3) + TwoTerms(10)
     tmpArray8(103,2) = Xpb*tmpArray7(75,2) + alphaXpq*TmpArray7(75,3) + TwoTerms(11)
     tmpArray8(104,2) = Ypb*tmpArray7(71,2) + alphaYpq*TmpArray7(71,3)
     tmpArray8(105,2) = Xpb*tmpArray7(77,2) + alphaXpq*TmpArray7(77,3) + TwoTerms(13)
     tmpArray8(106,2) = Xpb*tmpArray7(78,2) + alphaXpq*TmpArray7(78,3)
     tmpArray8(107,2) = Xpb*tmpArray7(79,2) + alphaXpq*TmpArray7(79,3)
     tmpArray8(108,2) = Xpb*tmpArray7(80,2) + alphaXpq*TmpArray7(80,3)
     tmpArray8(109,2) = Xpb*tmpArray7(81,2) + alphaXpq*TmpArray7(81,3)
     tmpArray8(110,2) = Xpb*tmpArray7(82,2) + alphaXpq*TmpArray7(82,3)
     tmpArray8(111,2) = Xpb*tmpArray7(83,2) + alphaXpq*TmpArray7(83,3)
     tmpArray8(112,2) = Xpb*tmpArray7(84,2) + alphaXpq*TmpArray7(84,3)
     tmpArray8(113,2) = Ypb*tmpArray7(78,2) + alphaYpq*TmpArray7(78,3) + 6*TwoTerms(9)
     tmpArray8(114,2) = Zpb*tmpArray7(78,2) + alphaZpq*TmpArray7(78,3)
     tmpArray8(115,2) = Ypb*tmpArray7(80,2) + alphaYpq*TmpArray7(80,3) + 4*TwoTerms(10)
     tmpArray8(116,2) = Ypb*tmpArray7(81,2) + alphaYpq*TmpArray7(81,3) + 3*TwoTerms(11)
     tmpArray8(117,2) = Ypb*tmpArray7(82,2) + alphaYpq*TmpArray7(82,3) + 2*TwoTerms(12)
     tmpArray8(118,2) = Ypb*tmpArray7(83,2) + alphaYpq*TmpArray7(83,3) + TwoTerms(13)
     tmpArray8(119,2) = Ypb*tmpArray7(84,2) + alphaYpq*TmpArray7(84,3)
     tmpArray8(120,2) = Zpb*tmpArray7(84,2) + alphaZpq*TmpArray7(84,3) + 6*TwoTerms(13)
     TwoTerms(1) = inv2expP*(TMPAuxArray(57) + alphaP*TmpArray7(57,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(60) + alphaP*TmpArray7(60,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(62) + alphaP*TmpArray7(62,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(63) + alphaP*TmpArray7(63,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(66) + alphaP*TmpArray7(66,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(67) + alphaP*TmpArray7(67,2))
     TwoTerms(7) = inv2expP*(TMPAuxArray(69) + alphaP*TmpArray7(69,2))
     TwoTerms(8) = inv2expP*(TMPAuxArray(71) + alphaP*TmpArray7(71,2))
     TwoTerms(9) = inv2expP*(TMPAuxArray(72) + alphaP*TmpArray7(72,2))
     TwoTerms(10) = inv2expP*(TMPAuxArray(74) + alphaP*TmpArray7(74,2))
     TwoTerms(11) = inv2expP*(TMPAuxArray(75) + alphaP*TmpArray7(75,2))
     TwoTerms(12) = inv2expP*(TMPAuxArray(77) + alphaP*TmpArray7(77,2))
     TwoTerms(13) = inv2expP*(TMPAuxArray(78) + alphaP*TmpArray7(78,2))
     TwoTerms(14) = inv2expP*(TMPAuxArray(80) + alphaP*TmpArray7(80,2))
     TwoTerms(15) = inv2expP*(TMPAuxArray(81) + alphaP*TmpArray7(81,2))
     TwoTerms(16) = inv2expP*(TMPAuxArray(82) + alphaP*TmpArray7(82,2))
     TwoTerms(17) = inv2expP*(TMPAuxArray(83) + alphaP*TmpArray7(83,2))
     TwoTerms(18) = inv2expP*(TMPAuxArray(84) + alphaP*TmpArray7(84,2))
     do iTUV = 1,  120
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(121,IPassQ) = AuxArray(121,IPassQ) + Xpb*TMPAuxArray(85) + alphaXpq*TmpArray8(85,2) + 7*TwoTerms(1)
     AuxArray(122,IPassQ) = AuxArray(122,IPassQ) + Ypb*TMPAuxArray(85) + alphaYpq*TmpArray8(85,2)
     AuxArray(123,IPassQ) = AuxArray(123,IPassQ) + Zpb*TMPAuxArray(85) + alphaZpq*TmpArray8(85,2)
     AuxArray(124,IPassQ) = AuxArray(124,IPassQ) + Xpb*TMPAuxArray(88) + alphaXpq*TmpArray8(88,2) + 5*TwoTerms(2)
     AuxArray(125,IPassQ) = AuxArray(125,IPassQ) + Ypb*TMPAuxArray(87) + alphaYpq*TmpArray8(87,2)
     AuxArray(126,IPassQ) = AuxArray(126,IPassQ) + Xpb*TMPAuxArray(90) + alphaXpq*TmpArray8(90,2) + 5*TwoTerms(3)
     AuxArray(127,IPassQ) = AuxArray(127,IPassQ) + Xpb*TMPAuxArray(91) + alphaXpq*TmpArray8(91,2) + 4*TwoTerms(4)
     AuxArray(128,IPassQ) = AuxArray(128,IPassQ) + Zpb*TMPAuxArray(88) + alphaZpq*TmpArray8(88,2)
     AuxArray(129,IPassQ) = AuxArray(129,IPassQ) + Ypb*TMPAuxArray(90) + alphaYpq*TmpArray8(90,2)
     AuxArray(130,IPassQ) = AuxArray(130,IPassQ) + Xpb*TMPAuxArray(94) + alphaXpq*TmpArray8(94,2) + 4*TwoTerms(5)
     AuxArray(131,IPassQ) = AuxArray(131,IPassQ) + Xpb*TMPAuxArray(95) + alphaXpq*TmpArray8(95,2) + 3*TwoTerms(6)
     AuxArray(132,IPassQ) = AuxArray(132,IPassQ) + Zpb*TMPAuxArray(91) + alphaZpq*TmpArray8(91,2)
     AuxArray(133,IPassQ) = AuxArray(133,IPassQ) + Xpb*TMPAuxArray(97) + alphaXpq*TmpArray8(97,2) + 3*TwoTerms(7)
     AuxArray(134,IPassQ) = AuxArray(134,IPassQ) + Ypb*TMPAuxArray(94) + alphaYpq*TmpArray8(94,2)
     AuxArray(135,IPassQ) = AuxArray(135,IPassQ) + Xpb*TMPAuxArray(99) + alphaXpq*TmpArray8(99,2) + 3*TwoTerms(8)
     AuxArray(136,IPassQ) = AuxArray(136,IPassQ) + Xpb*TMPAuxArray(100) + alphaXpq*TmpArray8(100,2) + 2*TwoTerms(9)
     AuxArray(137,IPassQ) = AuxArray(137,IPassQ) + Zpb*TMPAuxArray(95) + alphaZpq*TmpArray8(95,2)
     AuxArray(138,IPassQ) = AuxArray(138,IPassQ) + Xpb*TMPAuxArray(102) + alphaXpq*TmpArray8(102,2) + 2*TwoTerms(10)
     AuxArray(139,IPassQ) = AuxArray(139,IPassQ) + Xpb*TMPAuxArray(103) + alphaXpq*TmpArray8(103,2) + 2*TwoTerms(11)
     AuxArray(140,IPassQ) = AuxArray(140,IPassQ) + Ypb*TMPAuxArray(99) + alphaYpq*TmpArray8(99,2)
     AuxArray(141,IPassQ) = AuxArray(141,IPassQ) + Xpb*TMPAuxArray(105) + alphaXpq*TmpArray8(105,2) + 2*TwoTerms(12)
     AuxArray(142,IPassQ) = AuxArray(142,IPassQ) + Xpb*TMPAuxArray(106) + alphaXpq*TmpArray8(106,2) + TwoTerms(13)
     AuxArray(143,IPassQ) = AuxArray(143,IPassQ) + Zpb*TMPAuxArray(100) + alphaZpq*TmpArray8(100,2)
     AuxArray(144,IPassQ) = AuxArray(144,IPassQ) + Xpb*TMPAuxArray(108) + alphaXpq*TmpArray8(108,2) + TwoTerms(14)
     AuxArray(145,IPassQ) = AuxArray(145,IPassQ) + Xpb*TMPAuxArray(109) + alphaXpq*TmpArray8(109,2) + TwoTerms(15)
     AuxArray(146,IPassQ) = AuxArray(146,IPassQ) + Xpb*TMPAuxArray(110) + alphaXpq*TmpArray8(110,2) + TwoTerms(16)
     AuxArray(147,IPassQ) = AuxArray(147,IPassQ) + Ypb*TMPAuxArray(105) + alphaYpq*TmpArray8(105,2)
     AuxArray(148,IPassQ) = AuxArray(148,IPassQ) + Xpb*TMPAuxArray(112) + alphaXpq*TmpArray8(112,2) + TwoTerms(18)
     AuxArray(149,IPassQ) = AuxArray(149,IPassQ) + Xpb*TMPAuxArray(113) + alphaXpq*TmpArray8(113,2)
     AuxArray(150,IPassQ) = AuxArray(150,IPassQ) + Xpb*TMPAuxArray(114) + alphaXpq*TmpArray8(114,2)
     AuxArray(151,IPassQ) = AuxArray(151,IPassQ) + Xpb*TMPAuxArray(115) + alphaXpq*TmpArray8(115,2)
     AuxArray(152,IPassQ) = AuxArray(152,IPassQ) + Xpb*TMPAuxArray(116) + alphaXpq*TmpArray8(116,2)
     AuxArray(153,IPassQ) = AuxArray(153,IPassQ) + Xpb*TMPAuxArray(117) + alphaXpq*TmpArray8(117,2)
     AuxArray(154,IPassQ) = AuxArray(154,IPassQ) + Xpb*TMPAuxArray(118) + alphaXpq*TmpArray8(118,2)
     AuxArray(155,IPassQ) = AuxArray(155,IPassQ) + Xpb*TMPAuxArray(119) + alphaXpq*TmpArray8(119,2)
     AuxArray(156,IPassQ) = AuxArray(156,IPassQ) + Xpb*TMPAuxArray(120) + alphaXpq*TmpArray8(120,2)
     AuxArray(157,IPassQ) = AuxArray(157,IPassQ) + Ypb*TMPAuxArray(113) + alphaYpq*TmpArray8(113,2) + 7*TwoTerms(13)
     AuxArray(158,IPassQ) = AuxArray(158,IPassQ) + Zpb*TMPAuxArray(113) + alphaZpq*TmpArray8(113,2)
     AuxArray(159,IPassQ) = AuxArray(159,IPassQ) + Ypb*TMPAuxArray(115) + alphaYpq*TmpArray8(115,2) + 5*TwoTerms(14)
     AuxArray(160,IPassQ) = AuxArray(160,IPassQ) + Ypb*TMPAuxArray(116) + alphaYpq*TmpArray8(116,2) + 4*TwoTerms(15)
     AuxArray(161,IPassQ) = AuxArray(161,IPassQ) + Ypb*TMPAuxArray(117) + alphaYpq*TmpArray8(117,2) + 3*TwoTerms(16)
     AuxArray(162,IPassQ) = AuxArray(162,IPassQ) + Ypb*TMPAuxArray(118) + alphaYpq*TmpArray8(118,2) + 2*TwoTerms(17)
     AuxArray(163,IPassQ) = AuxArray(163,IPassQ) + Ypb*TMPAuxArray(119) + alphaYpq*TmpArray8(119,2) + TwoTerms(18)
     AuxArray(164,IPassQ) = AuxArray(164,IPassQ) + Ypb*TMPAuxArray(120) + alphaYpq*TmpArray8(120,2)
     AuxArray(165,IPassQ) = AuxArray(165,IPassQ) + Zpb*TMPAuxArray(120) + alphaZpq*TmpArray8(120,2) + 7*TwoTerms(18)
    ENDDO
   ENDDO
  ENDDO
 end subroutine
end module
