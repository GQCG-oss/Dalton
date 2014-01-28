MODULE AGC_CPU_OBS_VERTICALRECURRENCEMODDSeg
 use IchorPrecisionModule
  
 CONTAINS

subroutine VerticalRecurrenceCPUSeg1D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,&
         & integralPrefactor,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:4,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(4,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt
  real(realk) :: Dx,Dy,Dz,Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),alphaQ,RJ000(0:1)
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
  !ThetaAux(n,1,0,0) = Xqd*ThetaAux(n,0,0,0) + (-alpha/q*Xpq)*ThetaAux(n+1,0,0,0)
  !i = 0 last 2 term vanish
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
  DO iPassQ = 1,nPasses
   AUXarray(1,iPassQ)=0.0E0_realk
   AUXarray(2,iPassQ)=0.0E0_realk
   AUXarray(3,iPassQ)=0.0E0_realk
   AUXarray(4,iPassQ)=0.0E0_realk
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     Xqd = Qcent(1,iPrimQ,iPassQ) + Dx
     Yqd = Qcent(2,iPrimQ,iPassQ) + Dy
     Zqd = Qcent(3,iPrimQ,iPassQ) + Dz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
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
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac
     TMP1 = PREF*RJ000(0)
     TMP2 = PREF*RJ000(1)
     AUXarray(1,iPassQ) = AUXarray(1,iPassQ) + TMP1
     AUXarray(2,iPassQ) = AUXarray(2,iPassQ) + Xqd*TMP1 + alphaXpq*TMP2
     AUXarray(3,iPassQ) = AUXarray(3,iPassQ) + Yqd*TMP1 + alphaYpq*TMP2
     AUXarray(4,iPassQ) = AUXarray(4,iPassQ) + Zqd*TMP1 + alphaZpq*TMP2
    enddo
   enddo
  enddo
end subroutine

subroutine VerticalRecurrenceCPUSeg2D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 5,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   10,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(    4)
  real(realk) :: Dx,Dy,Dz,Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ,RJ000(0: 2)
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
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,   10
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     Xqd = Qcent(1,iPrimQ,iPassQ) + Dx
     Yqd = Qcent(2,iPrimQ,iPassQ) + Dy
     Zqd = Qcent(3,iPrimQ,iPassQ) + Dz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
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
     TMPAuxArray(2) = Xqd*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqd*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqd*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(1) + alphaQ*TmpArray1(1,2))
     do iTUV = 1,    4
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(5,IPassQ) = AuxArray(5,IPassQ) + Xqd*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IPassQ) = AuxArray(6,IPassQ) + Xqd*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IPassQ) = AuxArray(7,IPassQ) + Xqd*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IPassQ) = AuxArray(8,IPassQ) + Yqd*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IPassQ) = AuxArray(9,IPassQ) + Yqd*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IPassQ) = AuxArray(10,IPassQ) + Zqd*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg3D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 6,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   20,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   10)
  real(realk) :: Dx,Dy,Dz,Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ,RJ000(0: 3)
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
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,   20
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     Xqd = Qcent(1,iPrimQ,iPassQ) + Dx
     Yqd = Qcent(2,iPrimQ,iPassQ) + Dy
     Zqd = Qcent(3,iPrimQ,iPassQ) + Dz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
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
     do iTUV = 1,   10
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(11,IPassQ) = AuxArray(11,IPassQ) + Xqd*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IPassQ) = AuxArray(12,IPassQ) + Yqd*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IPassQ) = AuxArray(13,IPassQ) + Zqd*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IPassQ) = AuxArray(14,IPassQ) + Xqd*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IPassQ) = AuxArray(15,IPassQ) + Xqd*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IPassQ) = AuxArray(16,IPassQ) + Xqd*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IPassQ) = AuxArray(17,IPassQ) + Yqd*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IPassQ) = AuxArray(18,IPassQ) + Zqd*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IPassQ) = AuxArray(19,IPassQ) + Yqd*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IPassQ) = AuxArray(20,IPassQ) + Zqd*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg4D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 7,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   35,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   20)
  real(realk) :: Dx,Dy,Dz,Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ,RJ000(0: 4)
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
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,   35
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     Xqd = Qcent(1,iPrimQ,iPassQ) + Dx
     Yqd = Qcent(2,iPrimQ,iPassQ) + Dy
     Zqd = Qcent(3,iPrimQ,iPassQ) + Dz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
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
     TMPAuxArray(2) = Xqd*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqd*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqd*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqd*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqd*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqd*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xqd*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Yqd*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zqd*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
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
     TwoTerms(1) = inv2expQ*(TmpArray1(1,3) + alphaQ*TmpArray1(1,4))
     tmpArray3(5,3) = Xqd*tmpArray2(2,3) + alphaXpq*TmpArray2(2,4) + TwoTerms(1)
     tmpArray3(6,3) = Xqd*tmpArray2(3,3) + alphaXpq*TmpArray2(3,4)
     tmpArray3(7,3) = Xqd*tmpArray2(4,3) + alphaXpq*TmpArray2(4,4)
     tmpArray3(8,3) = Yqd*tmpArray2(3,3) + alphaYpq*TmpArray2(3,4) + TwoTerms(1)
     tmpArray3(9,3) = Yqd*tmpArray2(4,3) + alphaYpq*TmpArray2(4,4)
     tmpArray3(10,3) = Zqd*tmpArray2(4,3) + alphaZpq*TmpArray2(4,4) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(4) + alphaQ*TmpArray2(4,2))
     TMPAuxArray(11) = Xqd*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Yqd*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zqd*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xqd*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xqd*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xqd*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Yqd*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zqd*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Yqd*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zqd*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(5) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(8) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(10) + alphaQ*TmpArray3(10,2))
     do iTUV = 1,   20
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(21,IPassQ) = AuxArray(21,IPassQ) + Xqd*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IPassQ) = AuxArray(22,IPassQ) + Yqd*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IPassQ) = AuxArray(23,IPassQ) + Zqd*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IPassQ) = AuxArray(24,IPassQ) + Xqd*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IPassQ) = AuxArray(25,IPassQ) + Yqd*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IPassQ) = AuxArray(26,IPassQ) + Xqd*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IPassQ) = AuxArray(27,IPassQ) + Xqd*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IPassQ) = AuxArray(28,IPassQ) + Xqd*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IPassQ) = AuxArray(29,IPassQ) + Xqd*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IPassQ) = AuxArray(30,IPassQ) + Xqd*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IPassQ) = AuxArray(31,IPassQ) + Yqd*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IPassQ) = AuxArray(32,IPassQ) + Zqd*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IPassQ) = AuxArray(33,IPassQ) + Yqd*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IPassQ) = AuxArray(34,IPassQ) + Yqd*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IPassQ) = AuxArray(35,IPassQ) + Zqd*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg5D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 8,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   56,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   35)
  real(realk) :: Dx,Dy,Dz,Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ,RJ000(0: 5)
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
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,   56
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     Xqd = Qcent(1,iPrimQ,iPassQ) + Dx
     Yqd = Qcent(2,iPrimQ,iPassQ) + Dy
     Zqd = Qcent(3,iPrimQ,iPassQ) + Dz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
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
     TMPAuxArray(2) = Xqd*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqd*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqd*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(4) + alphaQ*TmpArray2(4,2))
     TMPAuxArray(11) = Xqd*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Yqd*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zqd*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xqd*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xqd*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xqd*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Yqd*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zqd*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Yqd*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zqd*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(5) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(8) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(10) + alphaQ*TmpArray3(10,2))
     TMPAuxArray(21) = Xqd*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Yqd*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zqd*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xqd*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Yqd*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xqd*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xqd*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xqd*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xqd*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xqd*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Yqd*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zqd*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Yqd*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Yqd*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zqd*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(11) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(14) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(16) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(17) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(19) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(20) + alphaQ*TmpArray4(20,2))
     do iTUV = 1,   35
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(36,IPassQ) = AuxArray(36,IPassQ) + Xqd*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IPassQ) = AuxArray(37,IPassQ) + Yqd*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IPassQ) = AuxArray(38,IPassQ) + Zqd*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IPassQ) = AuxArray(39,IPassQ) + Xqd*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IPassQ) = AuxArray(40,IPassQ) + Yqd*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IPassQ) = AuxArray(41,IPassQ) + Xqd*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IPassQ) = AuxArray(42,IPassQ) + Xqd*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IPassQ) = AuxArray(43,IPassQ) + Zqd*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IPassQ) = AuxArray(44,IPassQ) + Yqd*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IPassQ) = AuxArray(45,IPassQ) + Xqd*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IPassQ) = AuxArray(46,IPassQ) + Xqd*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IPassQ) = AuxArray(47,IPassQ) + Xqd*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IPassQ) = AuxArray(48,IPassQ) + Xqd*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IPassQ) = AuxArray(49,IPassQ) + Xqd*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IPassQ) = AuxArray(50,IPassQ) + Xqd*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IPassQ) = AuxArray(51,IPassQ) + Yqd*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IPassQ) = AuxArray(52,IPassQ) + Zqd*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IPassQ) = AuxArray(53,IPassQ) + Yqd*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IPassQ) = AuxArray(54,IPassQ) + Yqd*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IPassQ) = AuxArray(55,IPassQ) + Yqd*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IPassQ) = AuxArray(56,IPassQ) + Zqd*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg6D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 9,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   84,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   56)
  real(realk) :: Dx,Dy,Dz,Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ,RJ000(0: 6)
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
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,   84
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     Xqd = Qcent(1,iPrimQ,iPassQ) + Dx
     Yqd = Qcent(2,iPrimQ,iPassQ) + Dy
     Zqd = Qcent(3,iPrimQ,iPassQ) + Dz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
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
     TMPAuxArray(2) = Xqd*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqd*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqd*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
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
     tmpArray2(2,6) = Xqd*tmpArray1(1,6) + alphaXpq*TmpArray1(1,7)
     tmpArray2(3,6) = Yqd*tmpArray1(1,6) + alphaYpq*TmpArray1(1,7)
     tmpArray2(4,6) = Zqd*tmpArray1(1,6) + alphaZpq*TmpArray1(1,7)
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
     TwoTerms(1) = inv2expQ*(TmpArray1(1,5) + alphaQ*TmpArray1(1,6))
     tmpArray3(5,5) = Xqd*tmpArray2(2,5) + alphaXpq*TmpArray2(2,6) + TwoTerms(1)
     tmpArray3(6,5) = Xqd*tmpArray2(3,5) + alphaXpq*TmpArray2(3,6)
     tmpArray3(7,5) = Xqd*tmpArray2(4,5) + alphaXpq*TmpArray2(4,6)
     tmpArray3(8,5) = Yqd*tmpArray2(3,5) + alphaYpq*TmpArray2(3,6) + TwoTerms(1)
     tmpArray3(9,5) = Yqd*tmpArray2(4,5) + alphaYpq*TmpArray2(4,6)
     tmpArray3(10,5) = Zqd*tmpArray2(4,5) + alphaZpq*TmpArray2(4,6) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(4) + alphaQ*TmpArray2(4,2))
     TMPAuxArray(11) = Xqd*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Yqd*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zqd*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xqd*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xqd*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xqd*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Yqd*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zqd*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Yqd*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zqd*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TmpArray2(2,4) + alphaQ*TmpArray2(2,5))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,4) + alphaQ*TmpArray2(3,5))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,4) + alphaQ*TmpArray2(4,5))
     tmpArray4(11,4) = Xqd*tmpArray3(5,4) + alphaXpq*TmpArray3(5,5) + 2*TwoTerms(1)
     tmpArray4(12,4) = Yqd*tmpArray3(5,4) + alphaYpq*TmpArray3(5,5)
     tmpArray4(13,4) = Zqd*tmpArray3(5,4) + alphaZpq*TmpArray3(5,5)
     tmpArray4(14,4) = Xqd*tmpArray3(8,4) + alphaXpq*TmpArray3(8,5)
     tmpArray4(15,4) = Xqd*tmpArray3(9,4) + alphaXpq*TmpArray3(9,5)
     tmpArray4(16,4) = Xqd*tmpArray3(10,4) + alphaXpq*TmpArray3(10,5)
     tmpArray4(17,4) = Yqd*tmpArray3(8,4) + alphaYpq*TmpArray3(8,5) + 2*TwoTerms(2)
     tmpArray4(18,4) = Zqd*tmpArray3(8,4) + alphaZpq*TmpArray3(8,5)
     tmpArray4(19,4) = Yqd*tmpArray3(10,4) + alphaYpq*TmpArray3(10,5)
     tmpArray4(20,4) = Zqd*tmpArray3(10,4) + alphaZpq*TmpArray3(10,5) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(5) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(8) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(10) + alphaQ*TmpArray3(10,2))
     TMPAuxArray(21) = Xqd*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Yqd*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zqd*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xqd*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Yqd*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xqd*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xqd*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xqd*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xqd*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xqd*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Yqd*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zqd*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Yqd*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Yqd*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zqd*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TmpArray3(5,3) + alphaQ*TmpArray3(5,4))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,3) + alphaQ*TmpArray3(8,4))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,3) + alphaQ*TmpArray3(10,4))
     tmpArray5(21,3) = Xqd*tmpArray4(11,3) + alphaXpq*TmpArray4(11,4) + 3*TwoTerms(1)
     tmpArray5(22,3) = Yqd*tmpArray4(11,3) + alphaYpq*TmpArray4(11,4)
     tmpArray5(23,3) = Zqd*tmpArray4(11,3) + alphaZpq*TmpArray4(11,4)
     tmpArray5(24,3) = Xqd*tmpArray4(14,3) + alphaXpq*TmpArray4(14,4) + TwoTerms(2)
     tmpArray5(25,3) = Yqd*tmpArray4(13,3) + alphaYpq*TmpArray4(13,4)
     tmpArray5(26,3) = Xqd*tmpArray4(16,3) + alphaXpq*TmpArray4(16,4) + TwoTerms(3)
     tmpArray5(27,3) = Xqd*tmpArray4(17,3) + alphaXpq*TmpArray4(17,4)
     tmpArray5(28,3) = Xqd*tmpArray4(18,3) + alphaXpq*TmpArray4(18,4)
     tmpArray5(29,3) = Xqd*tmpArray4(19,3) + alphaXpq*TmpArray4(19,4)
     tmpArray5(30,3) = Xqd*tmpArray4(20,3) + alphaXpq*TmpArray4(20,4)
     tmpArray5(31,3) = Yqd*tmpArray4(17,3) + alphaYpq*TmpArray4(17,4) + 3*TwoTerms(2)
     tmpArray5(32,3) = Zqd*tmpArray4(17,3) + alphaZpq*TmpArray4(17,4)
     tmpArray5(33,3) = Yqd*tmpArray4(19,3) + alphaYpq*TmpArray4(19,4) + TwoTerms(3)
     tmpArray5(34,3) = Yqd*tmpArray4(20,3) + alphaYpq*TmpArray4(20,4)
     tmpArray5(35,3) = Zqd*tmpArray4(20,3) + alphaZpq*TmpArray4(20,4) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(11) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(14) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(16) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(17) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(19) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(20) + alphaQ*TmpArray4(20,2))
     TMPAuxArray(36) = Xqd*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     TMPAuxArray(37) = Yqd*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     TMPAuxArray(38) = Zqd*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     TMPAuxArray(39) = Xqd*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     TMPAuxArray(40) = Yqd*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     TMPAuxArray(41) = Xqd*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     TMPAuxArray(42) = Xqd*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     TMPAuxArray(43) = Zqd*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     TMPAuxArray(44) = Yqd*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     TMPAuxArray(45) = Xqd*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     TMPAuxArray(46) = Xqd*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     TMPAuxArray(47) = Xqd*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     TMPAuxArray(48) = Xqd*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     TMPAuxArray(49) = Xqd*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     TMPAuxArray(50) = Xqd*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     TMPAuxArray(51) = Yqd*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     TMPAuxArray(52) = Zqd*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     TMPAuxArray(53) = Yqd*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     TMPAuxArray(54) = Yqd*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     TMPAuxArray(55) = Yqd*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     TMPAuxArray(56) = Zqd*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TmpArray4(11,2) + alphaQ*TmpArray4(11,3))
     TwoTerms(2) = inv2expQ*(TmpArray4(14,2) + alphaQ*TmpArray4(14,3))
     TwoTerms(3) = inv2expQ*(TmpArray4(16,2) + alphaQ*TmpArray4(16,3))
     TwoTerms(4) = inv2expQ*(TmpArray4(17,2) + alphaQ*TmpArray4(17,3))
     TwoTerms(5) = inv2expQ*(TmpArray4(19,2) + alphaQ*TmpArray4(19,3))
     TwoTerms(6) = inv2expQ*(TmpArray4(20,2) + alphaQ*TmpArray4(20,3))
     tmpArray6(36,2) = Xqd*tmpArray5(21,2) + alphaXpq*TmpArray5(21,3) + 4*TwoTerms(1)
     tmpArray6(37,2) = Yqd*tmpArray5(21,2) + alphaYpq*TmpArray5(21,3)
     tmpArray6(38,2) = Zqd*tmpArray5(21,2) + alphaZpq*TmpArray5(21,3)
     tmpArray6(39,2) = Xqd*tmpArray5(24,2) + alphaXpq*TmpArray5(24,3) + 2*TwoTerms(2)
     tmpArray6(40,2) = Yqd*tmpArray5(23,2) + alphaYpq*TmpArray5(23,3)
     tmpArray6(41,2) = Xqd*tmpArray5(26,2) + alphaXpq*TmpArray5(26,3) + 2*TwoTerms(3)
     tmpArray6(42,2) = Xqd*tmpArray5(27,2) + alphaXpq*TmpArray5(27,3) + TwoTerms(4)
     tmpArray6(43,2) = Zqd*tmpArray5(24,2) + alphaZpq*TmpArray5(24,3)
     tmpArray6(44,2) = Yqd*tmpArray5(26,2) + alphaYpq*TmpArray5(26,3)
     tmpArray6(45,2) = Xqd*tmpArray5(30,2) + alphaXpq*TmpArray5(30,3) + TwoTerms(6)
     tmpArray6(46,2) = Xqd*tmpArray5(31,2) + alphaXpq*TmpArray5(31,3)
     tmpArray6(47,2) = Xqd*tmpArray5(32,2) + alphaXpq*TmpArray5(32,3)
     tmpArray6(48,2) = Xqd*tmpArray5(33,2) + alphaXpq*TmpArray5(33,3)
     tmpArray6(49,2) = Xqd*tmpArray5(34,2) + alphaXpq*TmpArray5(34,3)
     tmpArray6(50,2) = Xqd*tmpArray5(35,2) + alphaXpq*TmpArray5(35,3)
     tmpArray6(51,2) = Yqd*tmpArray5(31,2) + alphaYpq*TmpArray5(31,3) + 4*TwoTerms(4)
     tmpArray6(52,2) = Zqd*tmpArray5(31,2) + alphaZpq*TmpArray5(31,3)
     tmpArray6(53,2) = Yqd*tmpArray5(33,2) + alphaYpq*TmpArray5(33,3) + 2*TwoTerms(5)
     tmpArray6(54,2) = Yqd*tmpArray5(34,2) + alphaYpq*TmpArray5(34,3) + TwoTerms(6)
     tmpArray6(55,2) = Yqd*tmpArray5(35,2) + alphaYpq*TmpArray5(35,3)
     tmpArray6(56,2) = Zqd*tmpArray5(35,2) + alphaZpq*TmpArray5(35,3) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(21) + alphaQ*TmpArray5(21,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(24) + alphaQ*TmpArray5(24,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(26) + alphaQ*TmpArray5(26,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(27) + alphaQ*TmpArray5(27,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(30) + alphaQ*TmpArray5(30,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(31) + alphaQ*TmpArray5(31,2))
     TwoTerms(7) = inv2expQ*(TMPAuxArray(33) + alphaQ*TmpArray5(33,2))
     TwoTerms(8) = inv2expQ*(TMPAuxArray(34) + alphaQ*TmpArray5(34,2))
     TwoTerms(9) = inv2expQ*(TMPAuxArray(35) + alphaQ*TmpArray5(35,2))
     do iTUV = 1,   56
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(57,IPassQ) = AuxArray(57,IPassQ) + Xqd*TMPAuxArray(36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,IPassQ) = AuxArray(58,IPassQ) + Yqd*TMPAuxArray(36) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,IPassQ) = AuxArray(59,IPassQ) + Zqd*TMPAuxArray(36) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,IPassQ) = AuxArray(60,IPassQ) + Xqd*TMPAuxArray(39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,IPassQ) = AuxArray(61,IPassQ) + Yqd*TMPAuxArray(38) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,IPassQ) = AuxArray(62,IPassQ) + Xqd*TMPAuxArray(41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,IPassQ) = AuxArray(63,IPassQ) + Xqd*TMPAuxArray(42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,IPassQ) = AuxArray(64,IPassQ) + Zqd*TMPAuxArray(39) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,IPassQ) = AuxArray(65,IPassQ) + Yqd*TMPAuxArray(41) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,IPassQ) = AuxArray(66,IPassQ) + Xqd*TMPAuxArray(45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,IPassQ) = AuxArray(67,IPassQ) + Xqd*TMPAuxArray(46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,IPassQ) = AuxArray(68,IPassQ) + Zqd*TMPAuxArray(42) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,IPassQ) = AuxArray(69,IPassQ) + Xqd*TMPAuxArray(48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,IPassQ) = AuxArray(70,IPassQ) + Yqd*TMPAuxArray(45) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,IPassQ) = AuxArray(71,IPassQ) + Xqd*TMPAuxArray(50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,IPassQ) = AuxArray(72,IPassQ) + Xqd*TMPAuxArray(51) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,IPassQ) = AuxArray(73,IPassQ) + Xqd*TMPAuxArray(52) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,IPassQ) = AuxArray(74,IPassQ) + Xqd*TMPAuxArray(53) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,IPassQ) = AuxArray(75,IPassQ) + Xqd*TMPAuxArray(54) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,IPassQ) = AuxArray(76,IPassQ) + Xqd*TMPAuxArray(55) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,IPassQ) = AuxArray(77,IPassQ) + Xqd*TMPAuxArray(56) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,IPassQ) = AuxArray(78,IPassQ) + Yqd*TMPAuxArray(51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,IPassQ) = AuxArray(79,IPassQ) + Zqd*TMPAuxArray(51) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,IPassQ) = AuxArray(80,IPassQ) + Yqd*TMPAuxArray(53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,IPassQ) = AuxArray(81,IPassQ) + Yqd*TMPAuxArray(54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,IPassQ) = AuxArray(82,IPassQ) + Yqd*TMPAuxArray(55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,IPassQ) = AuxArray(83,IPassQ) + Yqd*TMPAuxArray(56) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,IPassQ) = AuxArray(84,IPassQ) + Zqd*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg7D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:10,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(  120,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   84)
  real(realk) :: Dx,Dy,Dz,Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ,RJ000(0: 7)
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
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,  120
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     Xqd = Qcent(1,iPrimQ,iPassQ) + Dx
     Yqd = Qcent(2,iPrimQ,iPassQ) + Dy
     Zqd = Qcent(3,iPrimQ,iPassQ) + Dz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
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
     TMPAuxArray(2) = Xqd*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqd*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqd*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
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
     tmpArray2(2,6) = Xqd*tmpArray1(1,6) + alphaXpq*TmpArray1(1,7)
     tmpArray2(3,6) = Yqd*tmpArray1(1,6) + alphaYpq*TmpArray1(1,7)
     tmpArray2(4,6) = Zqd*tmpArray1(1,6) + alphaZpq*TmpArray1(1,7)
     tmpArray2(2,7) = Xqd*tmpArray1(1,7) + alphaXpq*TmpArray1(1,8)
     tmpArray2(3,7) = Yqd*tmpArray1(1,7) + alphaYpq*TmpArray1(1,8)
     tmpArray2(4,7) = Zqd*tmpArray1(1,7) + alphaZpq*TmpArray1(1,8)
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
     TwoTerms(1) = inv2expQ*(TmpArray1(1,5) + alphaQ*TmpArray1(1,6))
     tmpArray3(5,5) = Xqd*tmpArray2(2,5) + alphaXpq*TmpArray2(2,6) + TwoTerms(1)
     tmpArray3(6,5) = Xqd*tmpArray2(3,5) + alphaXpq*TmpArray2(3,6)
     tmpArray3(7,5) = Xqd*tmpArray2(4,5) + alphaXpq*TmpArray2(4,6)
     tmpArray3(8,5) = Yqd*tmpArray2(3,5) + alphaYpq*TmpArray2(3,6) + TwoTerms(1)
     tmpArray3(9,5) = Yqd*tmpArray2(4,5) + alphaYpq*TmpArray2(4,6)
     tmpArray3(10,5) = Zqd*tmpArray2(4,5) + alphaZpq*TmpArray2(4,6) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,6) + alphaQ*TmpArray1(1,7))
     tmpArray3(5,6) = Xqd*tmpArray2(2,6) + alphaXpq*TmpArray2(2,7) + TwoTerms(1)
     tmpArray3(6,6) = Xqd*tmpArray2(3,6) + alphaXpq*TmpArray2(3,7)
     tmpArray3(7,6) = Xqd*tmpArray2(4,6) + alphaXpq*TmpArray2(4,7)
     tmpArray3(8,6) = Yqd*tmpArray2(3,6) + alphaYpq*TmpArray2(3,7) + TwoTerms(1)
     tmpArray3(9,6) = Yqd*tmpArray2(4,6) + alphaYpq*TmpArray2(4,7)
     tmpArray3(10,6) = Zqd*tmpArray2(4,6) + alphaZpq*TmpArray2(4,7) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(4) + alphaQ*TmpArray2(4,2))
     TMPAuxArray(11) = Xqd*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Yqd*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zqd*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xqd*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xqd*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xqd*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Yqd*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zqd*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Yqd*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zqd*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TmpArray2(2,4) + alphaQ*TmpArray2(2,5))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,4) + alphaQ*TmpArray2(3,5))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,4) + alphaQ*TmpArray2(4,5))
     tmpArray4(11,4) = Xqd*tmpArray3(5,4) + alphaXpq*TmpArray3(5,5) + 2*TwoTerms(1)
     tmpArray4(12,4) = Yqd*tmpArray3(5,4) + alphaYpq*TmpArray3(5,5)
     tmpArray4(13,4) = Zqd*tmpArray3(5,4) + alphaZpq*TmpArray3(5,5)
     tmpArray4(14,4) = Xqd*tmpArray3(8,4) + alphaXpq*TmpArray3(8,5)
     tmpArray4(15,4) = Xqd*tmpArray3(9,4) + alphaXpq*TmpArray3(9,5)
     tmpArray4(16,4) = Xqd*tmpArray3(10,4) + alphaXpq*TmpArray3(10,5)
     tmpArray4(17,4) = Yqd*tmpArray3(8,4) + alphaYpq*TmpArray3(8,5) + 2*TwoTerms(2)
     tmpArray4(18,4) = Zqd*tmpArray3(8,4) + alphaZpq*TmpArray3(8,5)
     tmpArray4(19,4) = Yqd*tmpArray3(10,4) + alphaYpq*TmpArray3(10,5)
     tmpArray4(20,4) = Zqd*tmpArray3(10,4) + alphaZpq*TmpArray3(10,5) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,5) + alphaQ*TmpArray2(2,6))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,5) + alphaQ*TmpArray2(3,6))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,5) + alphaQ*TmpArray2(4,6))
     tmpArray4(11,5) = Xqd*tmpArray3(5,5) + alphaXpq*TmpArray3(5,6) + 2*TwoTerms(1)
     tmpArray4(12,5) = Yqd*tmpArray3(5,5) + alphaYpq*TmpArray3(5,6)
     tmpArray4(13,5) = Zqd*tmpArray3(5,5) + alphaZpq*TmpArray3(5,6)
     tmpArray4(14,5) = Xqd*tmpArray3(8,5) + alphaXpq*TmpArray3(8,6)
     tmpArray4(15,5) = Xqd*tmpArray3(9,5) + alphaXpq*TmpArray3(9,6)
     tmpArray4(16,5) = Xqd*tmpArray3(10,5) + alphaXpq*TmpArray3(10,6)
     tmpArray4(17,5) = Yqd*tmpArray3(8,5) + alphaYpq*TmpArray3(8,6) + 2*TwoTerms(2)
     tmpArray4(18,5) = Zqd*tmpArray3(8,5) + alphaZpq*TmpArray3(8,6)
     tmpArray4(19,5) = Yqd*tmpArray3(10,5) + alphaYpq*TmpArray3(10,6)
     tmpArray4(20,5) = Zqd*tmpArray3(10,5) + alphaZpq*TmpArray3(10,6) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(5) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(8) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(10) + alphaQ*TmpArray3(10,2))
     TMPAuxArray(21) = Xqd*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Yqd*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zqd*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xqd*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Yqd*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xqd*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xqd*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xqd*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xqd*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xqd*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Yqd*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zqd*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Yqd*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Yqd*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zqd*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TmpArray3(5,3) + alphaQ*TmpArray3(5,4))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,3) + alphaQ*TmpArray3(8,4))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,3) + alphaQ*TmpArray3(10,4))
     tmpArray5(21,3) = Xqd*tmpArray4(11,3) + alphaXpq*TmpArray4(11,4) + 3*TwoTerms(1)
     tmpArray5(22,3) = Yqd*tmpArray4(11,3) + alphaYpq*TmpArray4(11,4)
     tmpArray5(23,3) = Zqd*tmpArray4(11,3) + alphaZpq*TmpArray4(11,4)
     tmpArray5(24,3) = Xqd*tmpArray4(14,3) + alphaXpq*TmpArray4(14,4) + TwoTerms(2)
     tmpArray5(25,3) = Yqd*tmpArray4(13,3) + alphaYpq*TmpArray4(13,4)
     tmpArray5(26,3) = Xqd*tmpArray4(16,3) + alphaXpq*TmpArray4(16,4) + TwoTerms(3)
     tmpArray5(27,3) = Xqd*tmpArray4(17,3) + alphaXpq*TmpArray4(17,4)
     tmpArray5(28,3) = Xqd*tmpArray4(18,3) + alphaXpq*TmpArray4(18,4)
     tmpArray5(29,3) = Xqd*tmpArray4(19,3) + alphaXpq*TmpArray4(19,4)
     tmpArray5(30,3) = Xqd*tmpArray4(20,3) + alphaXpq*TmpArray4(20,4)
     tmpArray5(31,3) = Yqd*tmpArray4(17,3) + alphaYpq*TmpArray4(17,4) + 3*TwoTerms(2)
     tmpArray5(32,3) = Zqd*tmpArray4(17,3) + alphaZpq*TmpArray4(17,4)
     tmpArray5(33,3) = Yqd*tmpArray4(19,3) + alphaYpq*TmpArray4(19,4) + TwoTerms(3)
     tmpArray5(34,3) = Yqd*tmpArray4(20,3) + alphaYpq*TmpArray4(20,4)
     tmpArray5(35,3) = Zqd*tmpArray4(20,3) + alphaZpq*TmpArray4(20,4) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray3(5,4) + alphaQ*TmpArray3(5,5))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,4) + alphaQ*TmpArray3(8,5))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,4) + alphaQ*TmpArray3(10,5))
     tmpArray5(21,4) = Xqd*tmpArray4(11,4) + alphaXpq*TmpArray4(11,5) + 3*TwoTerms(1)
     tmpArray5(22,4) = Yqd*tmpArray4(11,4) + alphaYpq*TmpArray4(11,5)
     tmpArray5(23,4) = Zqd*tmpArray4(11,4) + alphaZpq*TmpArray4(11,5)
     tmpArray5(24,4) = Xqd*tmpArray4(14,4) + alphaXpq*TmpArray4(14,5) + TwoTerms(2)
     tmpArray5(25,4) = Yqd*tmpArray4(13,4) + alphaYpq*TmpArray4(13,5)
     tmpArray5(26,4) = Xqd*tmpArray4(16,4) + alphaXpq*TmpArray4(16,5) + TwoTerms(3)
     tmpArray5(27,4) = Xqd*tmpArray4(17,4) + alphaXpq*TmpArray4(17,5)
     tmpArray5(28,4) = Xqd*tmpArray4(18,4) + alphaXpq*TmpArray4(18,5)
     tmpArray5(29,4) = Xqd*tmpArray4(19,4) + alphaXpq*TmpArray4(19,5)
     tmpArray5(30,4) = Xqd*tmpArray4(20,4) + alphaXpq*TmpArray4(20,5)
     tmpArray5(31,4) = Yqd*tmpArray4(17,4) + alphaYpq*TmpArray4(17,5) + 3*TwoTerms(2)
     tmpArray5(32,4) = Zqd*tmpArray4(17,4) + alphaZpq*TmpArray4(17,5)
     tmpArray5(33,4) = Yqd*tmpArray4(19,4) + alphaYpq*TmpArray4(19,5) + TwoTerms(3)
     tmpArray5(34,4) = Yqd*tmpArray4(20,4) + alphaYpq*TmpArray4(20,5)
     tmpArray5(35,4) = Zqd*tmpArray4(20,4) + alphaZpq*TmpArray4(20,5) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(11) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(14) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(16) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(17) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(19) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(20) + alphaQ*TmpArray4(20,2))
     TMPAuxArray(36) = Xqd*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     TMPAuxArray(37) = Yqd*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     TMPAuxArray(38) = Zqd*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     TMPAuxArray(39) = Xqd*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     TMPAuxArray(40) = Yqd*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     TMPAuxArray(41) = Xqd*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     TMPAuxArray(42) = Xqd*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     TMPAuxArray(43) = Zqd*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     TMPAuxArray(44) = Yqd*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     TMPAuxArray(45) = Xqd*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     TMPAuxArray(46) = Xqd*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     TMPAuxArray(47) = Xqd*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     TMPAuxArray(48) = Xqd*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     TMPAuxArray(49) = Xqd*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     TMPAuxArray(50) = Xqd*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     TMPAuxArray(51) = Yqd*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     TMPAuxArray(52) = Zqd*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     TMPAuxArray(53) = Yqd*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     TMPAuxArray(54) = Yqd*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     TMPAuxArray(55) = Yqd*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     TMPAuxArray(56) = Zqd*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TmpArray4(11,2) + alphaQ*TmpArray4(11,3))
     TwoTerms(2) = inv2expQ*(TmpArray4(14,2) + alphaQ*TmpArray4(14,3))
     TwoTerms(3) = inv2expQ*(TmpArray4(16,2) + alphaQ*TmpArray4(16,3))
     TwoTerms(4) = inv2expQ*(TmpArray4(17,2) + alphaQ*TmpArray4(17,3))
     TwoTerms(5) = inv2expQ*(TmpArray4(19,2) + alphaQ*TmpArray4(19,3))
     TwoTerms(6) = inv2expQ*(TmpArray4(20,2) + alphaQ*TmpArray4(20,3))
     tmpArray6(36,2) = Xqd*tmpArray5(21,2) + alphaXpq*TmpArray5(21,3) + 4*TwoTerms(1)
     tmpArray6(37,2) = Yqd*tmpArray5(21,2) + alphaYpq*TmpArray5(21,3)
     tmpArray6(38,2) = Zqd*tmpArray5(21,2) + alphaZpq*TmpArray5(21,3)
     tmpArray6(39,2) = Xqd*tmpArray5(24,2) + alphaXpq*TmpArray5(24,3) + 2*TwoTerms(2)
     tmpArray6(40,2) = Yqd*tmpArray5(23,2) + alphaYpq*TmpArray5(23,3)
     tmpArray6(41,2) = Xqd*tmpArray5(26,2) + alphaXpq*TmpArray5(26,3) + 2*TwoTerms(3)
     tmpArray6(42,2) = Xqd*tmpArray5(27,2) + alphaXpq*TmpArray5(27,3) + TwoTerms(4)
     tmpArray6(43,2) = Zqd*tmpArray5(24,2) + alphaZpq*TmpArray5(24,3)
     tmpArray6(44,2) = Yqd*tmpArray5(26,2) + alphaYpq*TmpArray5(26,3)
     tmpArray6(45,2) = Xqd*tmpArray5(30,2) + alphaXpq*TmpArray5(30,3) + TwoTerms(6)
     tmpArray6(46,2) = Xqd*tmpArray5(31,2) + alphaXpq*TmpArray5(31,3)
     tmpArray6(47,2) = Xqd*tmpArray5(32,2) + alphaXpq*TmpArray5(32,3)
     tmpArray6(48,2) = Xqd*tmpArray5(33,2) + alphaXpq*TmpArray5(33,3)
     tmpArray6(49,2) = Xqd*tmpArray5(34,2) + alphaXpq*TmpArray5(34,3)
     tmpArray6(50,2) = Xqd*tmpArray5(35,2) + alphaXpq*TmpArray5(35,3)
     tmpArray6(51,2) = Yqd*tmpArray5(31,2) + alphaYpq*TmpArray5(31,3) + 4*TwoTerms(4)
     tmpArray6(52,2) = Zqd*tmpArray5(31,2) + alphaZpq*TmpArray5(31,3)
     tmpArray6(53,2) = Yqd*tmpArray5(33,2) + alphaYpq*TmpArray5(33,3) + 2*TwoTerms(5)
     tmpArray6(54,2) = Yqd*tmpArray5(34,2) + alphaYpq*TmpArray5(34,3) + TwoTerms(6)
     tmpArray6(55,2) = Yqd*tmpArray5(35,2) + alphaYpq*TmpArray5(35,3)
     tmpArray6(56,2) = Zqd*tmpArray5(35,2) + alphaZpq*TmpArray5(35,3) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TmpArray4(11,3) + alphaQ*TmpArray4(11,4))
     TwoTerms(2) = inv2expQ*(TmpArray4(14,3) + alphaQ*TmpArray4(14,4))
     TwoTerms(3) = inv2expQ*(TmpArray4(16,3) + alphaQ*TmpArray4(16,4))
     TwoTerms(4) = inv2expQ*(TmpArray4(17,3) + alphaQ*TmpArray4(17,4))
     TwoTerms(5) = inv2expQ*(TmpArray4(19,3) + alphaQ*TmpArray4(19,4))
     TwoTerms(6) = inv2expQ*(TmpArray4(20,3) + alphaQ*TmpArray4(20,4))
     tmpArray6(36,3) = Xqd*tmpArray5(21,3) + alphaXpq*TmpArray5(21,4) + 4*TwoTerms(1)
     tmpArray6(37,3) = Yqd*tmpArray5(21,3) + alphaYpq*TmpArray5(21,4)
     tmpArray6(38,3) = Zqd*tmpArray5(21,3) + alphaZpq*TmpArray5(21,4)
     tmpArray6(39,3) = Xqd*tmpArray5(24,3) + alphaXpq*TmpArray5(24,4) + 2*TwoTerms(2)
     tmpArray6(40,3) = Yqd*tmpArray5(23,3) + alphaYpq*TmpArray5(23,4)
     tmpArray6(41,3) = Xqd*tmpArray5(26,3) + alphaXpq*TmpArray5(26,4) + 2*TwoTerms(3)
     tmpArray6(42,3) = Xqd*tmpArray5(27,3) + alphaXpq*TmpArray5(27,4) + TwoTerms(4)
     tmpArray6(43,3) = Zqd*tmpArray5(24,3) + alphaZpq*TmpArray5(24,4)
     tmpArray6(44,3) = Yqd*tmpArray5(26,3) + alphaYpq*TmpArray5(26,4)
     tmpArray6(45,3) = Xqd*tmpArray5(30,3) + alphaXpq*TmpArray5(30,4) + TwoTerms(6)
     tmpArray6(46,3) = Xqd*tmpArray5(31,3) + alphaXpq*TmpArray5(31,4)
     tmpArray6(47,3) = Xqd*tmpArray5(32,3) + alphaXpq*TmpArray5(32,4)
     tmpArray6(48,3) = Xqd*tmpArray5(33,3) + alphaXpq*TmpArray5(33,4)
     tmpArray6(49,3) = Xqd*tmpArray5(34,3) + alphaXpq*TmpArray5(34,4)
     tmpArray6(50,3) = Xqd*tmpArray5(35,3) + alphaXpq*TmpArray5(35,4)
     tmpArray6(51,3) = Yqd*tmpArray5(31,3) + alphaYpq*TmpArray5(31,4) + 4*TwoTerms(4)
     tmpArray6(52,3) = Zqd*tmpArray5(31,3) + alphaZpq*TmpArray5(31,4)
     tmpArray6(53,3) = Yqd*tmpArray5(33,3) + alphaYpq*TmpArray5(33,4) + 2*TwoTerms(5)
     tmpArray6(54,3) = Yqd*tmpArray5(34,3) + alphaYpq*TmpArray5(34,4) + TwoTerms(6)
     tmpArray6(55,3) = Yqd*tmpArray5(35,3) + alphaYpq*TmpArray5(35,4)
     tmpArray6(56,3) = Zqd*tmpArray5(35,3) + alphaZpq*TmpArray5(35,4) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(21) + alphaQ*TmpArray5(21,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(24) + alphaQ*TmpArray5(24,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(26) + alphaQ*TmpArray5(26,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(27) + alphaQ*TmpArray5(27,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(30) + alphaQ*TmpArray5(30,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(31) + alphaQ*TmpArray5(31,2))
     TwoTerms(7) = inv2expQ*(TMPAuxArray(33) + alphaQ*TmpArray5(33,2))
     TwoTerms(8) = inv2expQ*(TMPAuxArray(34) + alphaQ*TmpArray5(34,2))
     TwoTerms(9) = inv2expQ*(TMPAuxArray(35) + alphaQ*TmpArray5(35,2))
     TMPAuxArray(57) = Xqd*TMPAuxArray(36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     TMPAuxArray(58) = Yqd*TMPAuxArray(36) + alphaYpq*TmpArray6(36,2)
     TMPAuxArray(59) = Zqd*TMPAuxArray(36) + alphaZpq*TmpArray6(36,2)
     TMPAuxArray(60) = Xqd*TMPAuxArray(39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     TMPAuxArray(61) = Yqd*TMPAuxArray(38) + alphaYpq*TmpArray6(38,2)
     TMPAuxArray(62) = Xqd*TMPAuxArray(41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     TMPAuxArray(63) = Xqd*TMPAuxArray(42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     TMPAuxArray(64) = Zqd*TMPAuxArray(39) + alphaZpq*TmpArray6(39,2)
     TMPAuxArray(65) = Yqd*TMPAuxArray(41) + alphaYpq*TmpArray6(41,2)
     TMPAuxArray(66) = Xqd*TMPAuxArray(45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     TMPAuxArray(67) = Xqd*TMPAuxArray(46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     TMPAuxArray(68) = Zqd*TMPAuxArray(42) + alphaZpq*TmpArray6(42,2)
     TMPAuxArray(69) = Xqd*TMPAuxArray(48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     TMPAuxArray(70) = Yqd*TMPAuxArray(45) + alphaYpq*TmpArray6(45,2)
     TMPAuxArray(71) = Xqd*TMPAuxArray(50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     TMPAuxArray(72) = Xqd*TMPAuxArray(51) + alphaXpq*TmpArray6(51,2)
     TMPAuxArray(73) = Xqd*TMPAuxArray(52) + alphaXpq*TmpArray6(52,2)
     TMPAuxArray(74) = Xqd*TMPAuxArray(53) + alphaXpq*TmpArray6(53,2)
     TMPAuxArray(75) = Xqd*TMPAuxArray(54) + alphaXpq*TmpArray6(54,2)
     TMPAuxArray(76) = Xqd*TMPAuxArray(55) + alphaXpq*TmpArray6(55,2)
     TMPAuxArray(77) = Xqd*TMPAuxArray(56) + alphaXpq*TmpArray6(56,2)
     TMPAuxArray(78) = Yqd*TMPAuxArray(51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     TMPAuxArray(79) = Zqd*TMPAuxArray(51) + alphaZpq*TmpArray6(51,2)
     TMPAuxArray(80) = Yqd*TMPAuxArray(53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     TMPAuxArray(81) = Yqd*TMPAuxArray(54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     TMPAuxArray(82) = Yqd*TMPAuxArray(55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     TMPAuxArray(83) = Yqd*TMPAuxArray(56) + alphaYpq*TmpArray6(56,2)
     TMPAuxArray(84) = Zqd*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expQ*(TmpArray5(21,2) + alphaQ*TmpArray5(21,3))
     TwoTerms(2) = inv2expQ*(TmpArray5(24,2) + alphaQ*TmpArray5(24,3))
     TwoTerms(3) = inv2expQ*(TmpArray5(26,2) + alphaQ*TmpArray5(26,3))
     TwoTerms(4) = inv2expQ*(TmpArray5(27,2) + alphaQ*TmpArray5(27,3))
     TwoTerms(5) = inv2expQ*(TmpArray5(30,2) + alphaQ*TmpArray5(30,3))
     TwoTerms(6) = inv2expQ*(TmpArray5(31,2) + alphaQ*TmpArray5(31,3))
     TwoTerms(7) = inv2expQ*(TmpArray5(33,2) + alphaQ*TmpArray5(33,3))
     TwoTerms(8) = inv2expQ*(TmpArray5(34,2) + alphaQ*TmpArray5(34,3))
     TwoTerms(9) = inv2expQ*(TmpArray5(35,2) + alphaQ*TmpArray5(35,3))
     tmpArray7(57,2) = Xqd*tmpArray6(36,2) + alphaXpq*TmpArray6(36,3) + 5*TwoTerms(1)
     tmpArray7(58,2) = Yqd*tmpArray6(36,2) + alphaYpq*TmpArray6(36,3)
     tmpArray7(59,2) = Zqd*tmpArray6(36,2) + alphaZpq*TmpArray6(36,3)
     tmpArray7(60,2) = Xqd*tmpArray6(39,2) + alphaXpq*TmpArray6(39,3) + 3*TwoTerms(2)
     tmpArray7(61,2) = Yqd*tmpArray6(38,2) + alphaYpq*TmpArray6(38,3)
     tmpArray7(62,2) = Xqd*tmpArray6(41,2) + alphaXpq*TmpArray6(41,3) + 3*TwoTerms(3)
     tmpArray7(63,2) = Xqd*tmpArray6(42,2) + alphaXpq*TmpArray6(42,3) + 2*TwoTerms(4)
     tmpArray7(64,2) = Zqd*tmpArray6(39,2) + alphaZpq*TmpArray6(39,3)
     tmpArray7(65,2) = Yqd*tmpArray6(41,2) + alphaYpq*TmpArray6(41,3)
     tmpArray7(66,2) = Xqd*tmpArray6(45,2) + alphaXpq*TmpArray6(45,3) + 2*TwoTerms(5)
     tmpArray7(67,2) = Xqd*tmpArray6(46,2) + alphaXpq*TmpArray6(46,3) + TwoTerms(6)
     tmpArray7(68,2) = Zqd*tmpArray6(42,2) + alphaZpq*TmpArray6(42,3)
     tmpArray7(69,2) = Xqd*tmpArray6(48,2) + alphaXpq*TmpArray6(48,3) + TwoTerms(7)
     tmpArray7(70,2) = Yqd*tmpArray6(45,2) + alphaYpq*TmpArray6(45,3)
     tmpArray7(71,2) = Xqd*tmpArray6(50,2) + alphaXpq*TmpArray6(50,3) + TwoTerms(9)
     tmpArray7(72,2) = Xqd*tmpArray6(51,2) + alphaXpq*TmpArray6(51,3)
     tmpArray7(73,2) = Xqd*tmpArray6(52,2) + alphaXpq*TmpArray6(52,3)
     tmpArray7(74,2) = Xqd*tmpArray6(53,2) + alphaXpq*TmpArray6(53,3)
     tmpArray7(75,2) = Xqd*tmpArray6(54,2) + alphaXpq*TmpArray6(54,3)
     tmpArray7(76,2) = Xqd*tmpArray6(55,2) + alphaXpq*TmpArray6(55,3)
     tmpArray7(77,2) = Xqd*tmpArray6(56,2) + alphaXpq*TmpArray6(56,3)
     tmpArray7(78,2) = Yqd*tmpArray6(51,2) + alphaYpq*TmpArray6(51,3) + 5*TwoTerms(6)
     tmpArray7(79,2) = Zqd*tmpArray6(51,2) + alphaZpq*TmpArray6(51,3)
     tmpArray7(80,2) = Yqd*tmpArray6(53,2) + alphaYpq*TmpArray6(53,3) + 3*TwoTerms(7)
     tmpArray7(81,2) = Yqd*tmpArray6(54,2) + alphaYpq*TmpArray6(54,3) + 2*TwoTerms(8)
     tmpArray7(82,2) = Yqd*tmpArray6(55,2) + alphaYpq*TmpArray6(55,3) + TwoTerms(9)
     tmpArray7(83,2) = Yqd*tmpArray6(56,2) + alphaYpq*TmpArray6(56,3)
     tmpArray7(84,2) = Zqd*tmpArray6(56,2) + alphaZpq*TmpArray6(56,3) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(36) + alphaQ*TmpArray6(36,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(39) + alphaQ*TmpArray6(39,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(41) + alphaQ*TmpArray6(41,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(42) + alphaQ*TmpArray6(42,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(45) + alphaQ*TmpArray6(45,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(46) + alphaQ*TmpArray6(46,2))
     TwoTerms(7) = inv2expQ*(TMPAuxArray(48) + alphaQ*TmpArray6(48,2))
     TwoTerms(8) = inv2expQ*(TMPAuxArray(50) + alphaQ*TmpArray6(50,2))
     TwoTerms(9) = inv2expQ*(TMPAuxArray(51) + alphaQ*TmpArray6(51,2))
     TwoTerms(10) = inv2expQ*(TMPAuxArray(53) + alphaQ*TmpArray6(53,2))
     TwoTerms(11) = inv2expQ*(TMPAuxArray(54) + alphaQ*TmpArray6(54,2))
     TwoTerms(12) = inv2expQ*(TMPAuxArray(55) + alphaQ*TmpArray6(55,2))
     TwoTerms(13) = inv2expQ*(TMPAuxArray(56) + alphaQ*TmpArray6(56,2))
     do iTUV = 1,   84
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(85,IPassQ) = AuxArray(85,IPassQ) + Xqd*TMPAuxArray(57) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     AuxArray(86,IPassQ) = AuxArray(86,IPassQ) + Yqd*TMPAuxArray(57) + alphaYpq*TmpArray7(57,2)
     AuxArray(87,IPassQ) = AuxArray(87,IPassQ) + Zqd*TMPAuxArray(57) + alphaZpq*TmpArray7(57,2)
     AuxArray(88,IPassQ) = AuxArray(88,IPassQ) + Xqd*TMPAuxArray(60) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     AuxArray(89,IPassQ) = AuxArray(89,IPassQ) + Yqd*TMPAuxArray(59) + alphaYpq*TmpArray7(59,2)
     AuxArray(90,IPassQ) = AuxArray(90,IPassQ) + Xqd*TMPAuxArray(62) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     AuxArray(91,IPassQ) = AuxArray(91,IPassQ) + Xqd*TMPAuxArray(63) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     AuxArray(92,IPassQ) = AuxArray(92,IPassQ) + Zqd*TMPAuxArray(60) + alphaZpq*TmpArray7(60,2)
     AuxArray(93,IPassQ) = AuxArray(93,IPassQ) + Yqd*TMPAuxArray(62) + alphaYpq*TmpArray7(62,2)
     AuxArray(94,IPassQ) = AuxArray(94,IPassQ) + Xqd*TMPAuxArray(66) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     AuxArray(95,IPassQ) = AuxArray(95,IPassQ) + Xqd*TMPAuxArray(67) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     AuxArray(96,IPassQ) = AuxArray(96,IPassQ) + Zqd*TMPAuxArray(63) + alphaZpq*TmpArray7(63,2)
     AuxArray(97,IPassQ) = AuxArray(97,IPassQ) + Xqd*TMPAuxArray(69) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     AuxArray(98,IPassQ) = AuxArray(98,IPassQ) + Yqd*TMPAuxArray(66) + alphaYpq*TmpArray7(66,2)
     AuxArray(99,IPassQ) = AuxArray(99,IPassQ) + Xqd*TMPAuxArray(71) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     AuxArray(100,IPassQ) = AuxArray(100,IPassQ) + Xqd*TMPAuxArray(72) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     AuxArray(101,IPassQ) = AuxArray(101,IPassQ) + Zqd*TMPAuxArray(67) + alphaZpq*TmpArray7(67,2)
     AuxArray(102,IPassQ) = AuxArray(102,IPassQ) + Xqd*TMPAuxArray(74) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     AuxArray(103,IPassQ) = AuxArray(103,IPassQ) + Xqd*TMPAuxArray(75) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     AuxArray(104,IPassQ) = AuxArray(104,IPassQ) + Yqd*TMPAuxArray(71) + alphaYpq*TmpArray7(71,2)
     AuxArray(105,IPassQ) = AuxArray(105,IPassQ) + Xqd*TMPAuxArray(77) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     AuxArray(106,IPassQ) = AuxArray(106,IPassQ) + Xqd*TMPAuxArray(78) + alphaXpq*TmpArray7(78,2)
     AuxArray(107,IPassQ) = AuxArray(107,IPassQ) + Xqd*TMPAuxArray(79) + alphaXpq*TmpArray7(79,2)
     AuxArray(108,IPassQ) = AuxArray(108,IPassQ) + Xqd*TMPAuxArray(80) + alphaXpq*TmpArray7(80,2)
     AuxArray(109,IPassQ) = AuxArray(109,IPassQ) + Xqd*TMPAuxArray(81) + alphaXpq*TmpArray7(81,2)
     AuxArray(110,IPassQ) = AuxArray(110,IPassQ) + Xqd*TMPAuxArray(82) + alphaXpq*TmpArray7(82,2)
     AuxArray(111,IPassQ) = AuxArray(111,IPassQ) + Xqd*TMPAuxArray(83) + alphaXpq*TmpArray7(83,2)
     AuxArray(112,IPassQ) = AuxArray(112,IPassQ) + Xqd*TMPAuxArray(84) + alphaXpq*TmpArray7(84,2)
     AuxArray(113,IPassQ) = AuxArray(113,IPassQ) + Yqd*TMPAuxArray(78) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     AuxArray(114,IPassQ) = AuxArray(114,IPassQ) + Zqd*TMPAuxArray(78) + alphaZpq*TmpArray7(78,2)
     AuxArray(115,IPassQ) = AuxArray(115,IPassQ) + Yqd*TMPAuxArray(80) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     AuxArray(116,IPassQ) = AuxArray(116,IPassQ) + Yqd*TMPAuxArray(81) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     AuxArray(117,IPassQ) = AuxArray(117,IPassQ) + Yqd*TMPAuxArray(82) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     AuxArray(118,IPassQ) = AuxArray(118,IPassQ) + Yqd*TMPAuxArray(83) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     AuxArray(119,IPassQ) = AuxArray(119,IPassQ) + Yqd*TMPAuxArray(84) + alphaYpq*TmpArray7(84,2)
     AuxArray(120,IPassQ) = AuxArray(120,IPassQ) + Zqd*TMPAuxArray(84) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSeg8D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:11,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(  165,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(  120)
  real(realk) :: Dx,Dy,Dz,Xqd,Yqd,Zqd
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ,RJ000(0: 8)
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
  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
  DO iPassQ = 1,nPasses
   iP = iPassQ
   DO iTUV=1,  165
    AUXarray(iTUV,iPassQ)=0.0E0_realk
   ENDDO
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)
     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)
     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)
     Xqd = Qcent(1,iPrimQ,iPassQ) + Dx
     Yqd = Qcent(2,iPrimQ,iPassQ) + Dy
     Zqd = Qcent(3,iPrimQ,iPassQ) + Dz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
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
     TMPAuxArray(2) = Xqd*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqd*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqd*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
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
     tmpArray2(2,6) = Xqd*tmpArray1(1,6) + alphaXpq*TmpArray1(1,7)
     tmpArray2(3,6) = Yqd*tmpArray1(1,6) + alphaYpq*TmpArray1(1,7)
     tmpArray2(4,6) = Zqd*tmpArray1(1,6) + alphaZpq*TmpArray1(1,7)
     tmpArray2(2,7) = Xqd*tmpArray1(1,7) + alphaXpq*TmpArray1(1,8)
     tmpArray2(3,7) = Yqd*tmpArray1(1,7) + alphaYpq*TmpArray1(1,8)
     tmpArray2(4,7) = Zqd*tmpArray1(1,7) + alphaZpq*TmpArray1(1,8)
     tmpArray2(2,8) = Xqd*tmpArray1(1,8) + alphaXpq*TmpArray1(1,9)
     tmpArray2(3,8) = Yqd*tmpArray1(1,8) + alphaYpq*TmpArray1(1,9)
     tmpArray2(4,8) = Zqd*tmpArray1(1,8) + alphaZpq*TmpArray1(1,9)
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
     TwoTerms(1) = inv2expQ*(TmpArray1(1,5) + alphaQ*TmpArray1(1,6))
     tmpArray3(5,5) = Xqd*tmpArray2(2,5) + alphaXpq*TmpArray2(2,6) + TwoTerms(1)
     tmpArray3(6,5) = Xqd*tmpArray2(3,5) + alphaXpq*TmpArray2(3,6)
     tmpArray3(7,5) = Xqd*tmpArray2(4,5) + alphaXpq*TmpArray2(4,6)
     tmpArray3(8,5) = Yqd*tmpArray2(3,5) + alphaYpq*TmpArray2(3,6) + TwoTerms(1)
     tmpArray3(9,5) = Yqd*tmpArray2(4,5) + alphaYpq*TmpArray2(4,6)
     tmpArray3(10,5) = Zqd*tmpArray2(4,5) + alphaZpq*TmpArray2(4,6) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,6) + alphaQ*TmpArray1(1,7))
     tmpArray3(5,6) = Xqd*tmpArray2(2,6) + alphaXpq*TmpArray2(2,7) + TwoTerms(1)
     tmpArray3(6,6) = Xqd*tmpArray2(3,6) + alphaXpq*TmpArray2(3,7)
     tmpArray3(7,6) = Xqd*tmpArray2(4,6) + alphaXpq*TmpArray2(4,7)
     tmpArray3(8,6) = Yqd*tmpArray2(3,6) + alphaYpq*TmpArray2(3,7) + TwoTerms(1)
     tmpArray3(9,6) = Yqd*tmpArray2(4,6) + alphaYpq*TmpArray2(4,7)
     tmpArray3(10,6) = Zqd*tmpArray2(4,6) + alphaZpq*TmpArray2(4,7) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,7) + alphaQ*TmpArray1(1,8))
     tmpArray3(5,7) = Xqd*tmpArray2(2,7) + alphaXpq*TmpArray2(2,8) + TwoTerms(1)
     tmpArray3(6,7) = Xqd*tmpArray2(3,7) + alphaXpq*TmpArray2(3,8)
     tmpArray3(7,7) = Xqd*tmpArray2(4,7) + alphaXpq*TmpArray2(4,8)
     tmpArray3(8,7) = Yqd*tmpArray2(3,7) + alphaYpq*TmpArray2(3,8) + TwoTerms(1)
     tmpArray3(9,7) = Yqd*tmpArray2(4,7) + alphaYpq*TmpArray2(4,8)
     tmpArray3(10,7) = Zqd*tmpArray2(4,7) + alphaZpq*TmpArray2(4,8) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(4) + alphaQ*TmpArray2(4,2))
     TMPAuxArray(11) = Xqd*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Yqd*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zqd*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xqd*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xqd*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xqd*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Yqd*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zqd*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Yqd*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zqd*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TmpArray2(2,4) + alphaQ*TmpArray2(2,5))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,4) + alphaQ*TmpArray2(3,5))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,4) + alphaQ*TmpArray2(4,5))
     tmpArray4(11,4) = Xqd*tmpArray3(5,4) + alphaXpq*TmpArray3(5,5) + 2*TwoTerms(1)
     tmpArray4(12,4) = Yqd*tmpArray3(5,4) + alphaYpq*TmpArray3(5,5)
     tmpArray4(13,4) = Zqd*tmpArray3(5,4) + alphaZpq*TmpArray3(5,5)
     tmpArray4(14,4) = Xqd*tmpArray3(8,4) + alphaXpq*TmpArray3(8,5)
     tmpArray4(15,4) = Xqd*tmpArray3(9,4) + alphaXpq*TmpArray3(9,5)
     tmpArray4(16,4) = Xqd*tmpArray3(10,4) + alphaXpq*TmpArray3(10,5)
     tmpArray4(17,4) = Yqd*tmpArray3(8,4) + alphaYpq*TmpArray3(8,5) + 2*TwoTerms(2)
     tmpArray4(18,4) = Zqd*tmpArray3(8,4) + alphaZpq*TmpArray3(8,5)
     tmpArray4(19,4) = Yqd*tmpArray3(10,4) + alphaYpq*TmpArray3(10,5)
     tmpArray4(20,4) = Zqd*tmpArray3(10,4) + alphaZpq*TmpArray3(10,5) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,5) + alphaQ*TmpArray2(2,6))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,5) + alphaQ*TmpArray2(3,6))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,5) + alphaQ*TmpArray2(4,6))
     tmpArray4(11,5) = Xqd*tmpArray3(5,5) + alphaXpq*TmpArray3(5,6) + 2*TwoTerms(1)
     tmpArray4(12,5) = Yqd*tmpArray3(5,5) + alphaYpq*TmpArray3(5,6)
     tmpArray4(13,5) = Zqd*tmpArray3(5,5) + alphaZpq*TmpArray3(5,6)
     tmpArray4(14,5) = Xqd*tmpArray3(8,5) + alphaXpq*TmpArray3(8,6)
     tmpArray4(15,5) = Xqd*tmpArray3(9,5) + alphaXpq*TmpArray3(9,6)
     tmpArray4(16,5) = Xqd*tmpArray3(10,5) + alphaXpq*TmpArray3(10,6)
     tmpArray4(17,5) = Yqd*tmpArray3(8,5) + alphaYpq*TmpArray3(8,6) + 2*TwoTerms(2)
     tmpArray4(18,5) = Zqd*tmpArray3(8,5) + alphaZpq*TmpArray3(8,6)
     tmpArray4(19,5) = Yqd*tmpArray3(10,5) + alphaYpq*TmpArray3(10,6)
     tmpArray4(20,5) = Zqd*tmpArray3(10,5) + alphaZpq*TmpArray3(10,6) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray2(2,6) + alphaQ*TmpArray2(2,7))
     TwoTerms(2) = inv2expQ*(TmpArray2(3,6) + alphaQ*TmpArray2(3,7))
     TwoTerms(3) = inv2expQ*(TmpArray2(4,6) + alphaQ*TmpArray2(4,7))
     tmpArray4(11,6) = Xqd*tmpArray3(5,6) + alphaXpq*TmpArray3(5,7) + 2*TwoTerms(1)
     tmpArray4(12,6) = Yqd*tmpArray3(5,6) + alphaYpq*TmpArray3(5,7)
     tmpArray4(13,6) = Zqd*tmpArray3(5,6) + alphaZpq*TmpArray3(5,7)
     tmpArray4(14,6) = Xqd*tmpArray3(8,6) + alphaXpq*TmpArray3(8,7)
     tmpArray4(15,6) = Xqd*tmpArray3(9,6) + alphaXpq*TmpArray3(9,7)
     tmpArray4(16,6) = Xqd*tmpArray3(10,6) + alphaXpq*TmpArray3(10,7)
     tmpArray4(17,6) = Yqd*tmpArray3(8,6) + alphaYpq*TmpArray3(8,7) + 2*TwoTerms(2)
     tmpArray4(18,6) = Zqd*tmpArray3(8,6) + alphaZpq*TmpArray3(8,7)
     tmpArray4(19,6) = Yqd*tmpArray3(10,6) + alphaYpq*TmpArray3(10,7)
     tmpArray4(20,6) = Zqd*tmpArray3(10,6) + alphaZpq*TmpArray3(10,7) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(5) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(8) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(10) + alphaQ*TmpArray3(10,2))
     TMPAuxArray(21) = Xqd*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Yqd*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zqd*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xqd*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Yqd*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xqd*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xqd*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xqd*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xqd*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xqd*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Yqd*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zqd*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Yqd*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Yqd*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zqd*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TmpArray3(5,3) + alphaQ*TmpArray3(5,4))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,3) + alphaQ*TmpArray3(8,4))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,3) + alphaQ*TmpArray3(10,4))
     tmpArray5(21,3) = Xqd*tmpArray4(11,3) + alphaXpq*TmpArray4(11,4) + 3*TwoTerms(1)
     tmpArray5(22,3) = Yqd*tmpArray4(11,3) + alphaYpq*TmpArray4(11,4)
     tmpArray5(23,3) = Zqd*tmpArray4(11,3) + alphaZpq*TmpArray4(11,4)
     tmpArray5(24,3) = Xqd*tmpArray4(14,3) + alphaXpq*TmpArray4(14,4) + TwoTerms(2)
     tmpArray5(25,3) = Yqd*tmpArray4(13,3) + alphaYpq*TmpArray4(13,4)
     tmpArray5(26,3) = Xqd*tmpArray4(16,3) + alphaXpq*TmpArray4(16,4) + TwoTerms(3)
     tmpArray5(27,3) = Xqd*tmpArray4(17,3) + alphaXpq*TmpArray4(17,4)
     tmpArray5(28,3) = Xqd*tmpArray4(18,3) + alphaXpq*TmpArray4(18,4)
     tmpArray5(29,3) = Xqd*tmpArray4(19,3) + alphaXpq*TmpArray4(19,4)
     tmpArray5(30,3) = Xqd*tmpArray4(20,3) + alphaXpq*TmpArray4(20,4)
     tmpArray5(31,3) = Yqd*tmpArray4(17,3) + alphaYpq*TmpArray4(17,4) + 3*TwoTerms(2)
     tmpArray5(32,3) = Zqd*tmpArray4(17,3) + alphaZpq*TmpArray4(17,4)
     tmpArray5(33,3) = Yqd*tmpArray4(19,3) + alphaYpq*TmpArray4(19,4) + TwoTerms(3)
     tmpArray5(34,3) = Yqd*tmpArray4(20,3) + alphaYpq*TmpArray4(20,4)
     tmpArray5(35,3) = Zqd*tmpArray4(20,3) + alphaZpq*TmpArray4(20,4) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray3(5,4) + alphaQ*TmpArray3(5,5))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,4) + alphaQ*TmpArray3(8,5))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,4) + alphaQ*TmpArray3(10,5))
     tmpArray5(21,4) = Xqd*tmpArray4(11,4) + alphaXpq*TmpArray4(11,5) + 3*TwoTerms(1)
     tmpArray5(22,4) = Yqd*tmpArray4(11,4) + alphaYpq*TmpArray4(11,5)
     tmpArray5(23,4) = Zqd*tmpArray4(11,4) + alphaZpq*TmpArray4(11,5)
     tmpArray5(24,4) = Xqd*tmpArray4(14,4) + alphaXpq*TmpArray4(14,5) + TwoTerms(2)
     tmpArray5(25,4) = Yqd*tmpArray4(13,4) + alphaYpq*TmpArray4(13,5)
     tmpArray5(26,4) = Xqd*tmpArray4(16,4) + alphaXpq*TmpArray4(16,5) + TwoTerms(3)
     tmpArray5(27,4) = Xqd*tmpArray4(17,4) + alphaXpq*TmpArray4(17,5)
     tmpArray5(28,4) = Xqd*tmpArray4(18,4) + alphaXpq*TmpArray4(18,5)
     tmpArray5(29,4) = Xqd*tmpArray4(19,4) + alphaXpq*TmpArray4(19,5)
     tmpArray5(30,4) = Xqd*tmpArray4(20,4) + alphaXpq*TmpArray4(20,5)
     tmpArray5(31,4) = Yqd*tmpArray4(17,4) + alphaYpq*TmpArray4(17,5) + 3*TwoTerms(2)
     tmpArray5(32,4) = Zqd*tmpArray4(17,4) + alphaZpq*TmpArray4(17,5)
     tmpArray5(33,4) = Yqd*tmpArray4(19,4) + alphaYpq*TmpArray4(19,5) + TwoTerms(3)
     tmpArray5(34,4) = Yqd*tmpArray4(20,4) + alphaYpq*TmpArray4(20,5)
     tmpArray5(35,4) = Zqd*tmpArray4(20,4) + alphaZpq*TmpArray4(20,5) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TmpArray3(5,5) + alphaQ*TmpArray3(5,6))
     TwoTerms(2) = inv2expQ*(TmpArray3(8,5) + alphaQ*TmpArray3(8,6))
     TwoTerms(3) = inv2expQ*(TmpArray3(10,5) + alphaQ*TmpArray3(10,6))
     tmpArray5(21,5) = Xqd*tmpArray4(11,5) + alphaXpq*TmpArray4(11,6) + 3*TwoTerms(1)
     tmpArray5(22,5) = Yqd*tmpArray4(11,5) + alphaYpq*TmpArray4(11,6)
     tmpArray5(23,5) = Zqd*tmpArray4(11,5) + alphaZpq*TmpArray4(11,6)
     tmpArray5(24,5) = Xqd*tmpArray4(14,5) + alphaXpq*TmpArray4(14,6) + TwoTerms(2)
     tmpArray5(25,5) = Yqd*tmpArray4(13,5) + alphaYpq*TmpArray4(13,6)
     tmpArray5(26,5) = Xqd*tmpArray4(16,5) + alphaXpq*TmpArray4(16,6) + TwoTerms(3)
     tmpArray5(27,5) = Xqd*tmpArray4(17,5) + alphaXpq*TmpArray4(17,6)
     tmpArray5(28,5) = Xqd*tmpArray4(18,5) + alphaXpq*TmpArray4(18,6)
     tmpArray5(29,5) = Xqd*tmpArray4(19,5) + alphaXpq*TmpArray4(19,6)
     tmpArray5(30,5) = Xqd*tmpArray4(20,5) + alphaXpq*TmpArray4(20,6)
     tmpArray5(31,5) = Yqd*tmpArray4(17,5) + alphaYpq*TmpArray4(17,6) + 3*TwoTerms(2)
     tmpArray5(32,5) = Zqd*tmpArray4(17,5) + alphaZpq*TmpArray4(17,6)
     tmpArray5(33,5) = Yqd*tmpArray4(19,5) + alphaYpq*TmpArray4(19,6) + TwoTerms(3)
     tmpArray5(34,5) = Yqd*tmpArray4(20,5) + alphaYpq*TmpArray4(20,6)
     tmpArray5(35,5) = Zqd*tmpArray4(20,5) + alphaZpq*TmpArray4(20,6) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(11) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(14) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(16) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(17) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(19) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(20) + alphaQ*TmpArray4(20,2))
     TMPAuxArray(36) = Xqd*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     TMPAuxArray(37) = Yqd*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     TMPAuxArray(38) = Zqd*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     TMPAuxArray(39) = Xqd*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     TMPAuxArray(40) = Yqd*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     TMPAuxArray(41) = Xqd*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     TMPAuxArray(42) = Xqd*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     TMPAuxArray(43) = Zqd*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     TMPAuxArray(44) = Yqd*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     TMPAuxArray(45) = Xqd*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     TMPAuxArray(46) = Xqd*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     TMPAuxArray(47) = Xqd*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     TMPAuxArray(48) = Xqd*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     TMPAuxArray(49) = Xqd*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     TMPAuxArray(50) = Xqd*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     TMPAuxArray(51) = Yqd*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     TMPAuxArray(52) = Zqd*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     TMPAuxArray(53) = Yqd*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     TMPAuxArray(54) = Yqd*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     TMPAuxArray(55) = Yqd*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     TMPAuxArray(56) = Zqd*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TmpArray4(11,2) + alphaQ*TmpArray4(11,3))
     TwoTerms(2) = inv2expQ*(TmpArray4(14,2) + alphaQ*TmpArray4(14,3))
     TwoTerms(3) = inv2expQ*(TmpArray4(16,2) + alphaQ*TmpArray4(16,3))
     TwoTerms(4) = inv2expQ*(TmpArray4(17,2) + alphaQ*TmpArray4(17,3))
     TwoTerms(5) = inv2expQ*(TmpArray4(19,2) + alphaQ*TmpArray4(19,3))
     TwoTerms(6) = inv2expQ*(TmpArray4(20,2) + alphaQ*TmpArray4(20,3))
     tmpArray6(36,2) = Xqd*tmpArray5(21,2) + alphaXpq*TmpArray5(21,3) + 4*TwoTerms(1)
     tmpArray6(37,2) = Yqd*tmpArray5(21,2) + alphaYpq*TmpArray5(21,3)
     tmpArray6(38,2) = Zqd*tmpArray5(21,2) + alphaZpq*TmpArray5(21,3)
     tmpArray6(39,2) = Xqd*tmpArray5(24,2) + alphaXpq*TmpArray5(24,3) + 2*TwoTerms(2)
     tmpArray6(40,2) = Yqd*tmpArray5(23,2) + alphaYpq*TmpArray5(23,3)
     tmpArray6(41,2) = Xqd*tmpArray5(26,2) + alphaXpq*TmpArray5(26,3) + 2*TwoTerms(3)
     tmpArray6(42,2) = Xqd*tmpArray5(27,2) + alphaXpq*TmpArray5(27,3) + TwoTerms(4)
     tmpArray6(43,2) = Zqd*tmpArray5(24,2) + alphaZpq*TmpArray5(24,3)
     tmpArray6(44,2) = Yqd*tmpArray5(26,2) + alphaYpq*TmpArray5(26,3)
     tmpArray6(45,2) = Xqd*tmpArray5(30,2) + alphaXpq*TmpArray5(30,3) + TwoTerms(6)
     tmpArray6(46,2) = Xqd*tmpArray5(31,2) + alphaXpq*TmpArray5(31,3)
     tmpArray6(47,2) = Xqd*tmpArray5(32,2) + alphaXpq*TmpArray5(32,3)
     tmpArray6(48,2) = Xqd*tmpArray5(33,2) + alphaXpq*TmpArray5(33,3)
     tmpArray6(49,2) = Xqd*tmpArray5(34,2) + alphaXpq*TmpArray5(34,3)
     tmpArray6(50,2) = Xqd*tmpArray5(35,2) + alphaXpq*TmpArray5(35,3)
     tmpArray6(51,2) = Yqd*tmpArray5(31,2) + alphaYpq*TmpArray5(31,3) + 4*TwoTerms(4)
     tmpArray6(52,2) = Zqd*tmpArray5(31,2) + alphaZpq*TmpArray5(31,3)
     tmpArray6(53,2) = Yqd*tmpArray5(33,2) + alphaYpq*TmpArray5(33,3) + 2*TwoTerms(5)
     tmpArray6(54,2) = Yqd*tmpArray5(34,2) + alphaYpq*TmpArray5(34,3) + TwoTerms(6)
     tmpArray6(55,2) = Yqd*tmpArray5(35,2) + alphaYpq*TmpArray5(35,3)
     tmpArray6(56,2) = Zqd*tmpArray5(35,2) + alphaZpq*TmpArray5(35,3) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TmpArray4(11,3) + alphaQ*TmpArray4(11,4))
     TwoTerms(2) = inv2expQ*(TmpArray4(14,3) + alphaQ*TmpArray4(14,4))
     TwoTerms(3) = inv2expQ*(TmpArray4(16,3) + alphaQ*TmpArray4(16,4))
     TwoTerms(4) = inv2expQ*(TmpArray4(17,3) + alphaQ*TmpArray4(17,4))
     TwoTerms(5) = inv2expQ*(TmpArray4(19,3) + alphaQ*TmpArray4(19,4))
     TwoTerms(6) = inv2expQ*(TmpArray4(20,3) + alphaQ*TmpArray4(20,4))
     tmpArray6(36,3) = Xqd*tmpArray5(21,3) + alphaXpq*TmpArray5(21,4) + 4*TwoTerms(1)
     tmpArray6(37,3) = Yqd*tmpArray5(21,3) + alphaYpq*TmpArray5(21,4)
     tmpArray6(38,3) = Zqd*tmpArray5(21,3) + alphaZpq*TmpArray5(21,4)
     tmpArray6(39,3) = Xqd*tmpArray5(24,3) + alphaXpq*TmpArray5(24,4) + 2*TwoTerms(2)
     tmpArray6(40,3) = Yqd*tmpArray5(23,3) + alphaYpq*TmpArray5(23,4)
     tmpArray6(41,3) = Xqd*tmpArray5(26,3) + alphaXpq*TmpArray5(26,4) + 2*TwoTerms(3)
     tmpArray6(42,3) = Xqd*tmpArray5(27,3) + alphaXpq*TmpArray5(27,4) + TwoTerms(4)
     tmpArray6(43,3) = Zqd*tmpArray5(24,3) + alphaZpq*TmpArray5(24,4)
     tmpArray6(44,3) = Yqd*tmpArray5(26,3) + alphaYpq*TmpArray5(26,4)
     tmpArray6(45,3) = Xqd*tmpArray5(30,3) + alphaXpq*TmpArray5(30,4) + TwoTerms(6)
     tmpArray6(46,3) = Xqd*tmpArray5(31,3) + alphaXpq*TmpArray5(31,4)
     tmpArray6(47,3) = Xqd*tmpArray5(32,3) + alphaXpq*TmpArray5(32,4)
     tmpArray6(48,3) = Xqd*tmpArray5(33,3) + alphaXpq*TmpArray5(33,4)
     tmpArray6(49,3) = Xqd*tmpArray5(34,3) + alphaXpq*TmpArray5(34,4)
     tmpArray6(50,3) = Xqd*tmpArray5(35,3) + alphaXpq*TmpArray5(35,4)
     tmpArray6(51,3) = Yqd*tmpArray5(31,3) + alphaYpq*TmpArray5(31,4) + 4*TwoTerms(4)
     tmpArray6(52,3) = Zqd*tmpArray5(31,3) + alphaZpq*TmpArray5(31,4)
     tmpArray6(53,3) = Yqd*tmpArray5(33,3) + alphaYpq*TmpArray5(33,4) + 2*TwoTerms(5)
     tmpArray6(54,3) = Yqd*tmpArray5(34,3) + alphaYpq*TmpArray5(34,4) + TwoTerms(6)
     tmpArray6(55,3) = Yqd*tmpArray5(35,3) + alphaYpq*TmpArray5(35,4)
     tmpArray6(56,3) = Zqd*tmpArray5(35,3) + alphaZpq*TmpArray5(35,4) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TmpArray4(11,4) + alphaQ*TmpArray4(11,5))
     TwoTerms(2) = inv2expQ*(TmpArray4(14,4) + alphaQ*TmpArray4(14,5))
     TwoTerms(3) = inv2expQ*(TmpArray4(16,4) + alphaQ*TmpArray4(16,5))
     TwoTerms(4) = inv2expQ*(TmpArray4(17,4) + alphaQ*TmpArray4(17,5))
     TwoTerms(5) = inv2expQ*(TmpArray4(19,4) + alphaQ*TmpArray4(19,5))
     TwoTerms(6) = inv2expQ*(TmpArray4(20,4) + alphaQ*TmpArray4(20,5))
     tmpArray6(36,4) = Xqd*tmpArray5(21,4) + alphaXpq*TmpArray5(21,5) + 4*TwoTerms(1)
     tmpArray6(37,4) = Yqd*tmpArray5(21,4) + alphaYpq*TmpArray5(21,5)
     tmpArray6(38,4) = Zqd*tmpArray5(21,4) + alphaZpq*TmpArray5(21,5)
     tmpArray6(39,4) = Xqd*tmpArray5(24,4) + alphaXpq*TmpArray5(24,5) + 2*TwoTerms(2)
     tmpArray6(40,4) = Yqd*tmpArray5(23,4) + alphaYpq*TmpArray5(23,5)
     tmpArray6(41,4) = Xqd*tmpArray5(26,4) + alphaXpq*TmpArray5(26,5) + 2*TwoTerms(3)
     tmpArray6(42,4) = Xqd*tmpArray5(27,4) + alphaXpq*TmpArray5(27,5) + TwoTerms(4)
     tmpArray6(43,4) = Zqd*tmpArray5(24,4) + alphaZpq*TmpArray5(24,5)
     tmpArray6(44,4) = Yqd*tmpArray5(26,4) + alphaYpq*TmpArray5(26,5)
     tmpArray6(45,4) = Xqd*tmpArray5(30,4) + alphaXpq*TmpArray5(30,5) + TwoTerms(6)
     tmpArray6(46,4) = Xqd*tmpArray5(31,4) + alphaXpq*TmpArray5(31,5)
     tmpArray6(47,4) = Xqd*tmpArray5(32,4) + alphaXpq*TmpArray5(32,5)
     tmpArray6(48,4) = Xqd*tmpArray5(33,4) + alphaXpq*TmpArray5(33,5)
     tmpArray6(49,4) = Xqd*tmpArray5(34,4) + alphaXpq*TmpArray5(34,5)
     tmpArray6(50,4) = Xqd*tmpArray5(35,4) + alphaXpq*TmpArray5(35,5)
     tmpArray6(51,4) = Yqd*tmpArray5(31,4) + alphaYpq*TmpArray5(31,5) + 4*TwoTerms(4)
     tmpArray6(52,4) = Zqd*tmpArray5(31,4) + alphaZpq*TmpArray5(31,5)
     tmpArray6(53,4) = Yqd*tmpArray5(33,4) + alphaYpq*TmpArray5(33,5) + 2*TwoTerms(5)
     tmpArray6(54,4) = Yqd*tmpArray5(34,4) + alphaYpq*TmpArray5(34,5) + TwoTerms(6)
     tmpArray6(55,4) = Yqd*tmpArray5(35,4) + alphaYpq*TmpArray5(35,5)
     tmpArray6(56,4) = Zqd*tmpArray5(35,4) + alphaZpq*TmpArray5(35,5) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(21) + alphaQ*TmpArray5(21,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(24) + alphaQ*TmpArray5(24,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(26) + alphaQ*TmpArray5(26,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(27) + alphaQ*TmpArray5(27,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(30) + alphaQ*TmpArray5(30,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(31) + alphaQ*TmpArray5(31,2))
     TwoTerms(7) = inv2expQ*(TMPAuxArray(33) + alphaQ*TmpArray5(33,2))
     TwoTerms(8) = inv2expQ*(TMPAuxArray(34) + alphaQ*TmpArray5(34,2))
     TwoTerms(9) = inv2expQ*(TMPAuxArray(35) + alphaQ*TmpArray5(35,2))
     TMPAuxArray(57) = Xqd*TMPAuxArray(36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     TMPAuxArray(58) = Yqd*TMPAuxArray(36) + alphaYpq*TmpArray6(36,2)
     TMPAuxArray(59) = Zqd*TMPAuxArray(36) + alphaZpq*TmpArray6(36,2)
     TMPAuxArray(60) = Xqd*TMPAuxArray(39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     TMPAuxArray(61) = Yqd*TMPAuxArray(38) + alphaYpq*TmpArray6(38,2)
     TMPAuxArray(62) = Xqd*TMPAuxArray(41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     TMPAuxArray(63) = Xqd*TMPAuxArray(42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     TMPAuxArray(64) = Zqd*TMPAuxArray(39) + alphaZpq*TmpArray6(39,2)
     TMPAuxArray(65) = Yqd*TMPAuxArray(41) + alphaYpq*TmpArray6(41,2)
     TMPAuxArray(66) = Xqd*TMPAuxArray(45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     TMPAuxArray(67) = Xqd*TMPAuxArray(46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     TMPAuxArray(68) = Zqd*TMPAuxArray(42) + alphaZpq*TmpArray6(42,2)
     TMPAuxArray(69) = Xqd*TMPAuxArray(48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     TMPAuxArray(70) = Yqd*TMPAuxArray(45) + alphaYpq*TmpArray6(45,2)
     TMPAuxArray(71) = Xqd*TMPAuxArray(50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     TMPAuxArray(72) = Xqd*TMPAuxArray(51) + alphaXpq*TmpArray6(51,2)
     TMPAuxArray(73) = Xqd*TMPAuxArray(52) + alphaXpq*TmpArray6(52,2)
     TMPAuxArray(74) = Xqd*TMPAuxArray(53) + alphaXpq*TmpArray6(53,2)
     TMPAuxArray(75) = Xqd*TMPAuxArray(54) + alphaXpq*TmpArray6(54,2)
     TMPAuxArray(76) = Xqd*TMPAuxArray(55) + alphaXpq*TmpArray6(55,2)
     TMPAuxArray(77) = Xqd*TMPAuxArray(56) + alphaXpq*TmpArray6(56,2)
     TMPAuxArray(78) = Yqd*TMPAuxArray(51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     TMPAuxArray(79) = Zqd*TMPAuxArray(51) + alphaZpq*TmpArray6(51,2)
     TMPAuxArray(80) = Yqd*TMPAuxArray(53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     TMPAuxArray(81) = Yqd*TMPAuxArray(54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     TMPAuxArray(82) = Yqd*TMPAuxArray(55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     TMPAuxArray(83) = Yqd*TMPAuxArray(56) + alphaYpq*TmpArray6(56,2)
     TMPAuxArray(84) = Zqd*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expQ*(TmpArray5(21,2) + alphaQ*TmpArray5(21,3))
     TwoTerms(2) = inv2expQ*(TmpArray5(24,2) + alphaQ*TmpArray5(24,3))
     TwoTerms(3) = inv2expQ*(TmpArray5(26,2) + alphaQ*TmpArray5(26,3))
     TwoTerms(4) = inv2expQ*(TmpArray5(27,2) + alphaQ*TmpArray5(27,3))
     TwoTerms(5) = inv2expQ*(TmpArray5(30,2) + alphaQ*TmpArray5(30,3))
     TwoTerms(6) = inv2expQ*(TmpArray5(31,2) + alphaQ*TmpArray5(31,3))
     TwoTerms(7) = inv2expQ*(TmpArray5(33,2) + alphaQ*TmpArray5(33,3))
     TwoTerms(8) = inv2expQ*(TmpArray5(34,2) + alphaQ*TmpArray5(34,3))
     TwoTerms(9) = inv2expQ*(TmpArray5(35,2) + alphaQ*TmpArray5(35,3))
     tmpArray7(57,2) = Xqd*tmpArray6(36,2) + alphaXpq*TmpArray6(36,3) + 5*TwoTerms(1)
     tmpArray7(58,2) = Yqd*tmpArray6(36,2) + alphaYpq*TmpArray6(36,3)
     tmpArray7(59,2) = Zqd*tmpArray6(36,2) + alphaZpq*TmpArray6(36,3)
     tmpArray7(60,2) = Xqd*tmpArray6(39,2) + alphaXpq*TmpArray6(39,3) + 3*TwoTerms(2)
     tmpArray7(61,2) = Yqd*tmpArray6(38,2) + alphaYpq*TmpArray6(38,3)
     tmpArray7(62,2) = Xqd*tmpArray6(41,2) + alphaXpq*TmpArray6(41,3) + 3*TwoTerms(3)
     tmpArray7(63,2) = Xqd*tmpArray6(42,2) + alphaXpq*TmpArray6(42,3) + 2*TwoTerms(4)
     tmpArray7(64,2) = Zqd*tmpArray6(39,2) + alphaZpq*TmpArray6(39,3)
     tmpArray7(65,2) = Yqd*tmpArray6(41,2) + alphaYpq*TmpArray6(41,3)
     tmpArray7(66,2) = Xqd*tmpArray6(45,2) + alphaXpq*TmpArray6(45,3) + 2*TwoTerms(5)
     tmpArray7(67,2) = Xqd*tmpArray6(46,2) + alphaXpq*TmpArray6(46,3) + TwoTerms(6)
     tmpArray7(68,2) = Zqd*tmpArray6(42,2) + alphaZpq*TmpArray6(42,3)
     tmpArray7(69,2) = Xqd*tmpArray6(48,2) + alphaXpq*TmpArray6(48,3) + TwoTerms(7)
     tmpArray7(70,2) = Yqd*tmpArray6(45,2) + alphaYpq*TmpArray6(45,3)
     tmpArray7(71,2) = Xqd*tmpArray6(50,2) + alphaXpq*TmpArray6(50,3) + TwoTerms(9)
     tmpArray7(72,2) = Xqd*tmpArray6(51,2) + alphaXpq*TmpArray6(51,3)
     tmpArray7(73,2) = Xqd*tmpArray6(52,2) + alphaXpq*TmpArray6(52,3)
     tmpArray7(74,2) = Xqd*tmpArray6(53,2) + alphaXpq*TmpArray6(53,3)
     tmpArray7(75,2) = Xqd*tmpArray6(54,2) + alphaXpq*TmpArray6(54,3)
     tmpArray7(76,2) = Xqd*tmpArray6(55,2) + alphaXpq*TmpArray6(55,3)
     tmpArray7(77,2) = Xqd*tmpArray6(56,2) + alphaXpq*TmpArray6(56,3)
     tmpArray7(78,2) = Yqd*tmpArray6(51,2) + alphaYpq*TmpArray6(51,3) + 5*TwoTerms(6)
     tmpArray7(79,2) = Zqd*tmpArray6(51,2) + alphaZpq*TmpArray6(51,3)
     tmpArray7(80,2) = Yqd*tmpArray6(53,2) + alphaYpq*TmpArray6(53,3) + 3*TwoTerms(7)
     tmpArray7(81,2) = Yqd*tmpArray6(54,2) + alphaYpq*TmpArray6(54,3) + 2*TwoTerms(8)
     tmpArray7(82,2) = Yqd*tmpArray6(55,2) + alphaYpq*TmpArray6(55,3) + TwoTerms(9)
     tmpArray7(83,2) = Yqd*tmpArray6(56,2) + alphaYpq*TmpArray6(56,3)
     tmpArray7(84,2) = Zqd*tmpArray6(56,2) + alphaZpq*TmpArray6(56,3) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expQ*(TmpArray5(21,3) + alphaQ*TmpArray5(21,4))
     TwoTerms(2) = inv2expQ*(TmpArray5(24,3) + alphaQ*TmpArray5(24,4))
     TwoTerms(3) = inv2expQ*(TmpArray5(26,3) + alphaQ*TmpArray5(26,4))
     TwoTerms(4) = inv2expQ*(TmpArray5(27,3) + alphaQ*TmpArray5(27,4))
     TwoTerms(5) = inv2expQ*(TmpArray5(30,3) + alphaQ*TmpArray5(30,4))
     TwoTerms(6) = inv2expQ*(TmpArray5(31,3) + alphaQ*TmpArray5(31,4))
     TwoTerms(7) = inv2expQ*(TmpArray5(33,3) + alphaQ*TmpArray5(33,4))
     TwoTerms(8) = inv2expQ*(TmpArray5(34,3) + alphaQ*TmpArray5(34,4))
     TwoTerms(9) = inv2expQ*(TmpArray5(35,3) + alphaQ*TmpArray5(35,4))
     tmpArray7(57,3) = Xqd*tmpArray6(36,3) + alphaXpq*TmpArray6(36,4) + 5*TwoTerms(1)
     tmpArray7(58,3) = Yqd*tmpArray6(36,3) + alphaYpq*TmpArray6(36,4)
     tmpArray7(59,3) = Zqd*tmpArray6(36,3) + alphaZpq*TmpArray6(36,4)
     tmpArray7(60,3) = Xqd*tmpArray6(39,3) + alphaXpq*TmpArray6(39,4) + 3*TwoTerms(2)
     tmpArray7(61,3) = Yqd*tmpArray6(38,3) + alphaYpq*TmpArray6(38,4)
     tmpArray7(62,3) = Xqd*tmpArray6(41,3) + alphaXpq*TmpArray6(41,4) + 3*TwoTerms(3)
     tmpArray7(63,3) = Xqd*tmpArray6(42,3) + alphaXpq*TmpArray6(42,4) + 2*TwoTerms(4)
     tmpArray7(64,3) = Zqd*tmpArray6(39,3) + alphaZpq*TmpArray6(39,4)
     tmpArray7(65,3) = Yqd*tmpArray6(41,3) + alphaYpq*TmpArray6(41,4)
     tmpArray7(66,3) = Xqd*tmpArray6(45,3) + alphaXpq*TmpArray6(45,4) + 2*TwoTerms(5)
     tmpArray7(67,3) = Xqd*tmpArray6(46,3) + alphaXpq*TmpArray6(46,4) + TwoTerms(6)
     tmpArray7(68,3) = Zqd*tmpArray6(42,3) + alphaZpq*TmpArray6(42,4)
     tmpArray7(69,3) = Xqd*tmpArray6(48,3) + alphaXpq*TmpArray6(48,4) + TwoTerms(7)
     tmpArray7(70,3) = Yqd*tmpArray6(45,3) + alphaYpq*TmpArray6(45,4)
     tmpArray7(71,3) = Xqd*tmpArray6(50,3) + alphaXpq*TmpArray6(50,4) + TwoTerms(9)
     tmpArray7(72,3) = Xqd*tmpArray6(51,3) + alphaXpq*TmpArray6(51,4)
     tmpArray7(73,3) = Xqd*tmpArray6(52,3) + alphaXpq*TmpArray6(52,4)
     tmpArray7(74,3) = Xqd*tmpArray6(53,3) + alphaXpq*TmpArray6(53,4)
     tmpArray7(75,3) = Xqd*tmpArray6(54,3) + alphaXpq*TmpArray6(54,4)
     tmpArray7(76,3) = Xqd*tmpArray6(55,3) + alphaXpq*TmpArray6(55,4)
     tmpArray7(77,3) = Xqd*tmpArray6(56,3) + alphaXpq*TmpArray6(56,4)
     tmpArray7(78,3) = Yqd*tmpArray6(51,3) + alphaYpq*TmpArray6(51,4) + 5*TwoTerms(6)
     tmpArray7(79,3) = Zqd*tmpArray6(51,3) + alphaZpq*TmpArray6(51,4)
     tmpArray7(80,3) = Yqd*tmpArray6(53,3) + alphaYpq*TmpArray6(53,4) + 3*TwoTerms(7)
     tmpArray7(81,3) = Yqd*tmpArray6(54,3) + alphaYpq*TmpArray6(54,4) + 2*TwoTerms(8)
     tmpArray7(82,3) = Yqd*tmpArray6(55,3) + alphaYpq*TmpArray6(55,4) + TwoTerms(9)
     tmpArray7(83,3) = Yqd*tmpArray6(56,3) + alphaYpq*TmpArray6(56,4)
     tmpArray7(84,3) = Zqd*tmpArray6(56,3) + alphaZpq*TmpArray6(56,4) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(36) + alphaQ*TmpArray6(36,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(39) + alphaQ*TmpArray6(39,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(41) + alphaQ*TmpArray6(41,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(42) + alphaQ*TmpArray6(42,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(45) + alphaQ*TmpArray6(45,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(46) + alphaQ*TmpArray6(46,2))
     TwoTerms(7) = inv2expQ*(TMPAuxArray(48) + alphaQ*TmpArray6(48,2))
     TwoTerms(8) = inv2expQ*(TMPAuxArray(50) + alphaQ*TmpArray6(50,2))
     TwoTerms(9) = inv2expQ*(TMPAuxArray(51) + alphaQ*TmpArray6(51,2))
     TwoTerms(10) = inv2expQ*(TMPAuxArray(53) + alphaQ*TmpArray6(53,2))
     TwoTerms(11) = inv2expQ*(TMPAuxArray(54) + alphaQ*TmpArray6(54,2))
     TwoTerms(12) = inv2expQ*(TMPAuxArray(55) + alphaQ*TmpArray6(55,2))
     TwoTerms(13) = inv2expQ*(TMPAuxArray(56) + alphaQ*TmpArray6(56,2))
     TMPAuxArray(85) = Xqd*TMPAuxArray(57) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     TMPAuxArray(86) = Yqd*TMPAuxArray(57) + alphaYpq*TmpArray7(57,2)
     TMPAuxArray(87) = Zqd*TMPAuxArray(57) + alphaZpq*TmpArray7(57,2)
     TMPAuxArray(88) = Xqd*TMPAuxArray(60) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     TMPAuxArray(89) = Yqd*TMPAuxArray(59) + alphaYpq*TmpArray7(59,2)
     TMPAuxArray(90) = Xqd*TMPAuxArray(62) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     TMPAuxArray(91) = Xqd*TMPAuxArray(63) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     TMPAuxArray(92) = Zqd*TMPAuxArray(60) + alphaZpq*TmpArray7(60,2)
     TMPAuxArray(93) = Yqd*TMPAuxArray(62) + alphaYpq*TmpArray7(62,2)
     TMPAuxArray(94) = Xqd*TMPAuxArray(66) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     TMPAuxArray(95) = Xqd*TMPAuxArray(67) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     TMPAuxArray(96) = Zqd*TMPAuxArray(63) + alphaZpq*TmpArray7(63,2)
     TMPAuxArray(97) = Xqd*TMPAuxArray(69) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     TMPAuxArray(98) = Yqd*TMPAuxArray(66) + alphaYpq*TmpArray7(66,2)
     TMPAuxArray(99) = Xqd*TMPAuxArray(71) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     TMPAuxArray(100) = Xqd*TMPAuxArray(72) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     TMPAuxArray(101) = Zqd*TMPAuxArray(67) + alphaZpq*TmpArray7(67,2)
     TMPAuxArray(102) = Xqd*TMPAuxArray(74) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     TMPAuxArray(103) = Xqd*TMPAuxArray(75) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     TMPAuxArray(104) = Yqd*TMPAuxArray(71) + alphaYpq*TmpArray7(71,2)
     TMPAuxArray(105) = Xqd*TMPAuxArray(77) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     TMPAuxArray(106) = Xqd*TMPAuxArray(78) + alphaXpq*TmpArray7(78,2)
     TMPAuxArray(107) = Xqd*TMPAuxArray(79) + alphaXpq*TmpArray7(79,2)
     TMPAuxArray(108) = Xqd*TMPAuxArray(80) + alphaXpq*TmpArray7(80,2)
     TMPAuxArray(109) = Xqd*TMPAuxArray(81) + alphaXpq*TmpArray7(81,2)
     TMPAuxArray(110) = Xqd*TMPAuxArray(82) + alphaXpq*TmpArray7(82,2)
     TMPAuxArray(111) = Xqd*TMPAuxArray(83) + alphaXpq*TmpArray7(83,2)
     TMPAuxArray(112) = Xqd*TMPAuxArray(84) + alphaXpq*TmpArray7(84,2)
     TMPAuxArray(113) = Yqd*TMPAuxArray(78) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     TMPAuxArray(114) = Zqd*TMPAuxArray(78) + alphaZpq*TmpArray7(78,2)
     TMPAuxArray(115) = Yqd*TMPAuxArray(80) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     TMPAuxArray(116) = Yqd*TMPAuxArray(81) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     TMPAuxArray(117) = Yqd*TMPAuxArray(82) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     TMPAuxArray(118) = Yqd*TMPAuxArray(83) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     TMPAuxArray(119) = Yqd*TMPAuxArray(84) + alphaYpq*TmpArray7(84,2)
     TMPAuxArray(120) = Zqd*TMPAuxArray(84) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
     TwoTerms(1) = inv2expQ*(TmpArray6(36,2) + alphaQ*TmpArray6(36,3))
     TwoTerms(2) = inv2expQ*(TmpArray6(39,2) + alphaQ*TmpArray6(39,3))
     TwoTerms(3) = inv2expQ*(TmpArray6(41,2) + alphaQ*TmpArray6(41,3))
     TwoTerms(4) = inv2expQ*(TmpArray6(42,2) + alphaQ*TmpArray6(42,3))
     TwoTerms(5) = inv2expQ*(TmpArray6(45,2) + alphaQ*TmpArray6(45,3))
     TwoTerms(6) = inv2expQ*(TmpArray6(46,2) + alphaQ*TmpArray6(46,3))
     TwoTerms(7) = inv2expQ*(TmpArray6(48,2) + alphaQ*TmpArray6(48,3))
     TwoTerms(8) = inv2expQ*(TmpArray6(50,2) + alphaQ*TmpArray6(50,3))
     TwoTerms(9) = inv2expQ*(TmpArray6(51,2) + alphaQ*TmpArray6(51,3))
     TwoTerms(10) = inv2expQ*(TmpArray6(53,2) + alphaQ*TmpArray6(53,3))
     TwoTerms(11) = inv2expQ*(TmpArray6(54,2) + alphaQ*TmpArray6(54,3))
     TwoTerms(12) = inv2expQ*(TmpArray6(55,2) + alphaQ*TmpArray6(55,3))
     TwoTerms(13) = inv2expQ*(TmpArray6(56,2) + alphaQ*TmpArray6(56,3))
     tmpArray8(85,2) = Xqd*tmpArray7(57,2) + alphaXpq*TmpArray7(57,3) + 6*TwoTerms(1)
     tmpArray8(86,2) = Yqd*tmpArray7(57,2) + alphaYpq*TmpArray7(57,3)
     tmpArray8(87,2) = Zqd*tmpArray7(57,2) + alphaZpq*TmpArray7(57,3)
     tmpArray8(88,2) = Xqd*tmpArray7(60,2) + alphaXpq*TmpArray7(60,3) + 4*TwoTerms(2)
     tmpArray8(89,2) = Yqd*tmpArray7(59,2) + alphaYpq*TmpArray7(59,3)
     tmpArray8(90,2) = Xqd*tmpArray7(62,2) + alphaXpq*TmpArray7(62,3) + 4*TwoTerms(3)
     tmpArray8(91,2) = Xqd*tmpArray7(63,2) + alphaXpq*TmpArray7(63,3) + 3*TwoTerms(4)
     tmpArray8(92,2) = Zqd*tmpArray7(60,2) + alphaZpq*TmpArray7(60,3)
     tmpArray8(93,2) = Yqd*tmpArray7(62,2) + alphaYpq*TmpArray7(62,3)
     tmpArray8(94,2) = Xqd*tmpArray7(66,2) + alphaXpq*TmpArray7(66,3) + 3*TwoTerms(5)
     tmpArray8(95,2) = Xqd*tmpArray7(67,2) + alphaXpq*TmpArray7(67,3) + 2*TwoTerms(6)
     tmpArray8(96,2) = Zqd*tmpArray7(63,2) + alphaZpq*TmpArray7(63,3)
     tmpArray8(97,2) = Xqd*tmpArray7(69,2) + alphaXpq*TmpArray7(69,3) + 2*TwoTerms(7)
     tmpArray8(98,2) = Yqd*tmpArray7(66,2) + alphaYpq*TmpArray7(66,3)
     tmpArray8(99,2) = Xqd*tmpArray7(71,2) + alphaXpq*TmpArray7(71,3) + 2*TwoTerms(8)
     tmpArray8(100,2) = Xqd*tmpArray7(72,2) + alphaXpq*TmpArray7(72,3) + TwoTerms(9)
     tmpArray8(101,2) = Zqd*tmpArray7(67,2) + alphaZpq*TmpArray7(67,3)
     tmpArray8(102,2) = Xqd*tmpArray7(74,2) + alphaXpq*TmpArray7(74,3) + TwoTerms(10)
     tmpArray8(103,2) = Xqd*tmpArray7(75,2) + alphaXpq*TmpArray7(75,3) + TwoTerms(11)
     tmpArray8(104,2) = Yqd*tmpArray7(71,2) + alphaYpq*TmpArray7(71,3)
     tmpArray8(105,2) = Xqd*tmpArray7(77,2) + alphaXpq*TmpArray7(77,3) + TwoTerms(13)
     tmpArray8(106,2) = Xqd*tmpArray7(78,2) + alphaXpq*TmpArray7(78,3)
     tmpArray8(107,2) = Xqd*tmpArray7(79,2) + alphaXpq*TmpArray7(79,3)
     tmpArray8(108,2) = Xqd*tmpArray7(80,2) + alphaXpq*TmpArray7(80,3)
     tmpArray8(109,2) = Xqd*tmpArray7(81,2) + alphaXpq*TmpArray7(81,3)
     tmpArray8(110,2) = Xqd*tmpArray7(82,2) + alphaXpq*TmpArray7(82,3)
     tmpArray8(111,2) = Xqd*tmpArray7(83,2) + alphaXpq*TmpArray7(83,3)
     tmpArray8(112,2) = Xqd*tmpArray7(84,2) + alphaXpq*TmpArray7(84,3)
     tmpArray8(113,2) = Yqd*tmpArray7(78,2) + alphaYpq*TmpArray7(78,3) + 6*TwoTerms(9)
     tmpArray8(114,2) = Zqd*tmpArray7(78,2) + alphaZpq*TmpArray7(78,3)
     tmpArray8(115,2) = Yqd*tmpArray7(80,2) + alphaYpq*TmpArray7(80,3) + 4*TwoTerms(10)
     tmpArray8(116,2) = Yqd*tmpArray7(81,2) + alphaYpq*TmpArray7(81,3) + 3*TwoTerms(11)
     tmpArray8(117,2) = Yqd*tmpArray7(82,2) + alphaYpq*TmpArray7(82,3) + 2*TwoTerms(12)
     tmpArray8(118,2) = Yqd*tmpArray7(83,2) + alphaYpq*TmpArray7(83,3) + TwoTerms(13)
     tmpArray8(119,2) = Yqd*tmpArray7(84,2) + alphaYpq*TmpArray7(84,3)
     tmpArray8(120,2) = Zqd*tmpArray7(84,2) + alphaZpq*TmpArray7(84,3) + 6*TwoTerms(13)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(57) + alphaQ*TmpArray7(57,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(60) + alphaQ*TmpArray7(60,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(62) + alphaQ*TmpArray7(62,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(63) + alphaQ*TmpArray7(63,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(66) + alphaQ*TmpArray7(66,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(67) + alphaQ*TmpArray7(67,2))
     TwoTerms(7) = inv2expQ*(TMPAuxArray(69) + alphaQ*TmpArray7(69,2))
     TwoTerms(8) = inv2expQ*(TMPAuxArray(71) + alphaQ*TmpArray7(71,2))
     TwoTerms(9) = inv2expQ*(TMPAuxArray(72) + alphaQ*TmpArray7(72,2))
     TwoTerms(10) = inv2expQ*(TMPAuxArray(74) + alphaQ*TmpArray7(74,2))
     TwoTerms(11) = inv2expQ*(TMPAuxArray(75) + alphaQ*TmpArray7(75,2))
     TwoTerms(12) = inv2expQ*(TMPAuxArray(77) + alphaQ*TmpArray7(77,2))
     TwoTerms(13) = inv2expQ*(TMPAuxArray(78) + alphaQ*TmpArray7(78,2))
     TwoTerms(14) = inv2expQ*(TMPAuxArray(80) + alphaQ*TmpArray7(80,2))
     TwoTerms(15) = inv2expQ*(TMPAuxArray(81) + alphaQ*TmpArray7(81,2))
     TwoTerms(16) = inv2expQ*(TMPAuxArray(82) + alphaQ*TmpArray7(82,2))
     TwoTerms(17) = inv2expQ*(TMPAuxArray(83) + alphaQ*TmpArray7(83,2))
     TwoTerms(18) = inv2expQ*(TMPAuxArray(84) + alphaQ*TmpArray7(84,2))
     do iTUV = 1,  120
      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)
     enddo
     AuxArray(121,IPassQ) = AuxArray(121,IPassQ) + Xqd*TMPAuxArray(85) + alphaXpq*TmpArray8(85,2) + 7*TwoTerms(1)
     AuxArray(122,IPassQ) = AuxArray(122,IPassQ) + Yqd*TMPAuxArray(85) + alphaYpq*TmpArray8(85,2)
     AuxArray(123,IPassQ) = AuxArray(123,IPassQ) + Zqd*TMPAuxArray(85) + alphaZpq*TmpArray8(85,2)
     AuxArray(124,IPassQ) = AuxArray(124,IPassQ) + Xqd*TMPAuxArray(88) + alphaXpq*TmpArray8(88,2) + 5*TwoTerms(2)
     AuxArray(125,IPassQ) = AuxArray(125,IPassQ) + Yqd*TMPAuxArray(87) + alphaYpq*TmpArray8(87,2)
     AuxArray(126,IPassQ) = AuxArray(126,IPassQ) + Xqd*TMPAuxArray(90) + alphaXpq*TmpArray8(90,2) + 5*TwoTerms(3)
     AuxArray(127,IPassQ) = AuxArray(127,IPassQ) + Xqd*TMPAuxArray(91) + alphaXpq*TmpArray8(91,2) + 4*TwoTerms(4)
     AuxArray(128,IPassQ) = AuxArray(128,IPassQ) + Zqd*TMPAuxArray(88) + alphaZpq*TmpArray8(88,2)
     AuxArray(129,IPassQ) = AuxArray(129,IPassQ) + Yqd*TMPAuxArray(90) + alphaYpq*TmpArray8(90,2)
     AuxArray(130,IPassQ) = AuxArray(130,IPassQ) + Xqd*TMPAuxArray(94) + alphaXpq*TmpArray8(94,2) + 4*TwoTerms(5)
     AuxArray(131,IPassQ) = AuxArray(131,IPassQ) + Xqd*TMPAuxArray(95) + alphaXpq*TmpArray8(95,2) + 3*TwoTerms(6)
     AuxArray(132,IPassQ) = AuxArray(132,IPassQ) + Zqd*TMPAuxArray(91) + alphaZpq*TmpArray8(91,2)
     AuxArray(133,IPassQ) = AuxArray(133,IPassQ) + Xqd*TMPAuxArray(97) + alphaXpq*TmpArray8(97,2) + 3*TwoTerms(7)
     AuxArray(134,IPassQ) = AuxArray(134,IPassQ) + Yqd*TMPAuxArray(94) + alphaYpq*TmpArray8(94,2)
     AuxArray(135,IPassQ) = AuxArray(135,IPassQ) + Xqd*TMPAuxArray(99) + alphaXpq*TmpArray8(99,2) + 3*TwoTerms(8)
     AuxArray(136,IPassQ) = AuxArray(136,IPassQ) + Xqd*TMPAuxArray(100) + alphaXpq*TmpArray8(100,2) + 2*TwoTerms(9)
     AuxArray(137,IPassQ) = AuxArray(137,IPassQ) + Zqd*TMPAuxArray(95) + alphaZpq*TmpArray8(95,2)
     AuxArray(138,IPassQ) = AuxArray(138,IPassQ) + Xqd*TMPAuxArray(102) + alphaXpq*TmpArray8(102,2) + 2*TwoTerms(10)
     AuxArray(139,IPassQ) = AuxArray(139,IPassQ) + Xqd*TMPAuxArray(103) + alphaXpq*TmpArray8(103,2) + 2*TwoTerms(11)
     AuxArray(140,IPassQ) = AuxArray(140,IPassQ) + Yqd*TMPAuxArray(99) + alphaYpq*TmpArray8(99,2)
     AuxArray(141,IPassQ) = AuxArray(141,IPassQ) + Xqd*TMPAuxArray(105) + alphaXpq*TmpArray8(105,2) + 2*TwoTerms(12)
     AuxArray(142,IPassQ) = AuxArray(142,IPassQ) + Xqd*TMPAuxArray(106) + alphaXpq*TmpArray8(106,2) + TwoTerms(13)
     AuxArray(143,IPassQ) = AuxArray(143,IPassQ) + Zqd*TMPAuxArray(100) + alphaZpq*TmpArray8(100,2)
     AuxArray(144,IPassQ) = AuxArray(144,IPassQ) + Xqd*TMPAuxArray(108) + alphaXpq*TmpArray8(108,2) + TwoTerms(14)
     AuxArray(145,IPassQ) = AuxArray(145,IPassQ) + Xqd*TMPAuxArray(109) + alphaXpq*TmpArray8(109,2) + TwoTerms(15)
     AuxArray(146,IPassQ) = AuxArray(146,IPassQ) + Xqd*TMPAuxArray(110) + alphaXpq*TmpArray8(110,2) + TwoTerms(16)
     AuxArray(147,IPassQ) = AuxArray(147,IPassQ) + Yqd*TMPAuxArray(105) + alphaYpq*TmpArray8(105,2)
     AuxArray(148,IPassQ) = AuxArray(148,IPassQ) + Xqd*TMPAuxArray(112) + alphaXpq*TmpArray8(112,2) + TwoTerms(18)
     AuxArray(149,IPassQ) = AuxArray(149,IPassQ) + Xqd*TMPAuxArray(113) + alphaXpq*TmpArray8(113,2)
     AuxArray(150,IPassQ) = AuxArray(150,IPassQ) + Xqd*TMPAuxArray(114) + alphaXpq*TmpArray8(114,2)
     AuxArray(151,IPassQ) = AuxArray(151,IPassQ) + Xqd*TMPAuxArray(115) + alphaXpq*TmpArray8(115,2)
     AuxArray(152,IPassQ) = AuxArray(152,IPassQ) + Xqd*TMPAuxArray(116) + alphaXpq*TmpArray8(116,2)
     AuxArray(153,IPassQ) = AuxArray(153,IPassQ) + Xqd*TMPAuxArray(117) + alphaXpq*TmpArray8(117,2)
     AuxArray(154,IPassQ) = AuxArray(154,IPassQ) + Xqd*TMPAuxArray(118) + alphaXpq*TmpArray8(118,2)
     AuxArray(155,IPassQ) = AuxArray(155,IPassQ) + Xqd*TMPAuxArray(119) + alphaXpq*TmpArray8(119,2)
     AuxArray(156,IPassQ) = AuxArray(156,IPassQ) + Xqd*TMPAuxArray(120) + alphaXpq*TmpArray8(120,2)
     AuxArray(157,IPassQ) = AuxArray(157,IPassQ) + Yqd*TMPAuxArray(113) + alphaYpq*TmpArray8(113,2) + 7*TwoTerms(13)
     AuxArray(158,IPassQ) = AuxArray(158,IPassQ) + Zqd*TMPAuxArray(113) + alphaZpq*TmpArray8(113,2)
     AuxArray(159,IPassQ) = AuxArray(159,IPassQ) + Yqd*TMPAuxArray(115) + alphaYpq*TmpArray8(115,2) + 5*TwoTerms(14)
     AuxArray(160,IPassQ) = AuxArray(160,IPassQ) + Yqd*TMPAuxArray(116) + alphaYpq*TmpArray8(116,2) + 4*TwoTerms(15)
     AuxArray(161,IPassQ) = AuxArray(161,IPassQ) + Yqd*TMPAuxArray(117) + alphaYpq*TmpArray8(117,2) + 3*TwoTerms(16)
     AuxArray(162,IPassQ) = AuxArray(162,IPassQ) + Yqd*TMPAuxArray(118) + alphaYpq*TmpArray8(118,2) + 2*TwoTerms(17)
     AuxArray(163,IPassQ) = AuxArray(163,IPassQ) + Yqd*TMPAuxArray(119) + alphaYpq*TmpArray8(119,2) + TwoTerms(18)
     AuxArray(164,IPassQ) = AuxArray(164,IPassQ) + Yqd*TMPAuxArray(120) + alphaYpq*TmpArray8(120,2)
     AuxArray(165,IPassQ) = AuxArray(165,IPassQ) + Zqd*TMPAuxArray(120) + alphaZpq*TmpArray8(120,2) + 7*TwoTerms(18)
    ENDDO
   ENDDO
  ENDDO
 end subroutine
end module
