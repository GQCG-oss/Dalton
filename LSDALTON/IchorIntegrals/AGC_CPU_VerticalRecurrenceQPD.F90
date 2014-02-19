MODULE AGC_CPU_OBS_VERTICALRECURRENCEMODD
 use IchorPrecisionModule
  
 CONTAINS

subroutine VerticalRecurrenceCPU1D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,&
         & integralPrefactor,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:4,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(4,nPrimQ,nPrimP,nPasses)
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
     AUXarray(1,iPrimQ,iPrimP,iPassQ) = TMP1
     AUXarray(2,iPrimQ,iPrimP,iPassQ) = Xqd*TMP1 + alphaXpq*TMP2
     AUXarray(3,iPrimQ,iPrimP,iPassQ) = Yqd*TMP1 + alphaYpq*TMP2
     AUXarray(4,iPrimQ,iPrimP,iPassQ) = Zqd*TMP1 + alphaZpq*TMP2
    enddo
   enddo
  enddo
end subroutine

subroutine VerticalRecurrenceCPU2D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 5,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   10,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
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
   iP = (iPassQ-1)*nPrimQ*nPrimP
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
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
     Auxarray(1,IP) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     AuxArray(2,IP) = Xqd*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqd*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqd*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,IP) = Xqd*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqd*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqd*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqd*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqd*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqd*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPU3D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 6,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   20,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
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
   iP = (iPassQ-1)*nPrimQ*nPrimP
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
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
     Auxarray(1,IP) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
     AuxArray(2,IP) = Xqd*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqd*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqd*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqd*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqd*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqd*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqd*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqd*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqd*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,IP) = Xqd*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqd*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqd*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqd*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqd*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqd*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     AuxArray(11,IP) = Xqd*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Yqd*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zqd*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xqd*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xqd*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xqd*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Yqd*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zqd*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Yqd*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zqd*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPU4D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 7,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   35,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
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
   iP = (iPassQ-1)*nPrimQ*nPrimP
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
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
     Auxarray(1,IP) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
     TMParray1(1, 5) = PREF*RJ000( 4)
     AuxArray(2,IP) = Xqd*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqd*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqd*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     AuxArray(5,IP) = Xqd*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqd*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqd*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqd*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqd*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqd*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     AuxArray(11,IP) = Xqd*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Yqd*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zqd*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xqd*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xqd*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xqd*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Yqd*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zqd*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Yqd*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zqd*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     AuxArray(21,IP) = Xqd*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Yqd*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zqd*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xqd*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Yqd*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xqd*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xqd*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xqd*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xqd*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xqd*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Yqd*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zqd*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Yqd*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Yqd*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zqd*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPU5D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 8,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   56,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
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
   iP = (iPassQ-1)*nPrimQ*nPrimP
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
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
     Auxarray(1,IP) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
     TMParray1(1, 5) = PREF*RJ000( 4)
     TMParray1(1, 6) = PREF*RJ000( 5)
     AuxArray(2,IP) = Xqd*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqd*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqd*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     AuxArray(5,IP) = Xqd*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqd*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqd*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqd*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqd*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqd*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     AuxArray(11,IP) = Xqd*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Yqd*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zqd*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xqd*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xqd*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xqd*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Yqd*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zqd*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Yqd*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zqd*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     AuxArray(21,IP) = Xqd*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Yqd*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zqd*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xqd*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Yqd*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xqd*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xqd*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xqd*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xqd*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xqd*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Yqd*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zqd*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Yqd*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Yqd*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zqd*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     AuxArray(36,IP) = Xqd*AuxArray(21,IP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IP) = Yqd*AuxArray(21,IP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IP) = Zqd*AuxArray(21,IP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IP) = Xqd*AuxArray(24,IP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IP) = Yqd*AuxArray(23,IP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IP) = Xqd*AuxArray(26,IP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IP) = Xqd*AuxArray(27,IP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IP) = Zqd*AuxArray(24,IP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IP) = Yqd*AuxArray(26,IP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IP) = Xqd*AuxArray(30,IP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IP) = Xqd*AuxArray(31,IP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IP) = Xqd*AuxArray(32,IP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IP) = Xqd*AuxArray(33,IP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IP) = Xqd*AuxArray(34,IP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IP) = Xqd*AuxArray(35,IP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IP) = Yqd*AuxArray(31,IP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IP) = Zqd*AuxArray(31,IP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IP) = Yqd*AuxArray(33,IP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IP) = Yqd*AuxArray(34,IP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IP) = Yqd*AuxArray(35,IP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IP) = Zqd*AuxArray(35,IP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPU6D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 9,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(   84,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
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
   iP = (iPassQ-1)*nPrimQ*nPrimP
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
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
     Auxarray(1,IP) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
     TMParray1(1, 5) = PREF*RJ000( 4)
     TMParray1(1, 6) = PREF*RJ000( 5)
     TMParray1(1, 7) = PREF*RJ000( 6)
     AuxArray(2,IP) = Xqd*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqd*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqd*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,IP) = Xqd*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqd*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqd*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqd*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqd*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqd*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expQ*(AuxArray(2,IP) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(3,IP) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(4,IP) + alphaQ*TmpArray2(4,2))
     AuxArray(11,IP) = Xqd*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Yqd*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zqd*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xqd*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xqd*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xqd*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Yqd*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zqd*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Yqd*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zqd*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(AuxArray(5,IP) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(AuxArray(8,IP) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(AuxArray(10,IP) + alphaQ*TmpArray3(10,2))
     AuxArray(21,IP) = Xqd*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Yqd*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zqd*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xqd*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Yqd*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xqd*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xqd*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xqd*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xqd*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xqd*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Yqd*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zqd*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Yqd*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Yqd*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zqd*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(AuxArray(11,IP) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(AuxArray(14,IP) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(AuxArray(16,IP) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(AuxArray(17,IP) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(AuxArray(19,IP) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(AuxArray(20,IP) + alphaQ*TmpArray4(20,2))
     AuxArray(36,IP) = Xqd*AuxArray(21,IP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IP) = Yqd*AuxArray(21,IP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IP) = Zqd*AuxArray(21,IP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IP) = Xqd*AuxArray(24,IP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IP) = Yqd*AuxArray(23,IP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IP) = Xqd*AuxArray(26,IP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IP) = Xqd*AuxArray(27,IP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IP) = Zqd*AuxArray(24,IP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IP) = Yqd*AuxArray(26,IP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IP) = Xqd*AuxArray(30,IP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IP) = Xqd*AuxArray(31,IP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IP) = Xqd*AuxArray(32,IP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IP) = Xqd*AuxArray(33,IP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IP) = Xqd*AuxArray(34,IP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IP) = Xqd*AuxArray(35,IP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IP) = Yqd*AuxArray(31,IP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IP) = Zqd*AuxArray(31,IP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IP) = Yqd*AuxArray(33,IP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IP) = Yqd*AuxArray(34,IP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IP) = Yqd*AuxArray(35,IP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IP) = Zqd*AuxArray(35,IP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     TwoTerms(1) = inv2expQ*(AuxArray(21,IP) + alphaQ*TmpArray5(21,2))
     TwoTerms(2) = inv2expQ*(AuxArray(24,IP) + alphaQ*TmpArray5(24,2))
     TwoTerms(3) = inv2expQ*(AuxArray(26,IP) + alphaQ*TmpArray5(26,2))
     TwoTerms(4) = inv2expQ*(AuxArray(27,IP) + alphaQ*TmpArray5(27,2))
     TwoTerms(5) = inv2expQ*(AuxArray(30,IP) + alphaQ*TmpArray5(30,2))
     TwoTerms(6) = inv2expQ*(AuxArray(31,IP) + alphaQ*TmpArray5(31,2))
     TwoTerms(7) = inv2expQ*(AuxArray(33,IP) + alphaQ*TmpArray5(33,2))
     TwoTerms(8) = inv2expQ*(AuxArray(34,IP) + alphaQ*TmpArray5(34,2))
     TwoTerms(9) = inv2expQ*(AuxArray(35,IP) + alphaQ*TmpArray5(35,2))
     AuxArray(57,IP) = Xqd*AuxArray(36,IP) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,IP) = Yqd*AuxArray(36,IP) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,IP) = Zqd*AuxArray(36,IP) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,IP) = Xqd*AuxArray(39,IP) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,IP) = Yqd*AuxArray(38,IP) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,IP) = Xqd*AuxArray(41,IP) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,IP) = Xqd*AuxArray(42,IP) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,IP) = Zqd*AuxArray(39,IP) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,IP) = Yqd*AuxArray(41,IP) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,IP) = Xqd*AuxArray(45,IP) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,IP) = Xqd*AuxArray(46,IP) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,IP) = Zqd*AuxArray(42,IP) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,IP) = Xqd*AuxArray(48,IP) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,IP) = Yqd*AuxArray(45,IP) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,IP) = Xqd*AuxArray(50,IP) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,IP) = Xqd*AuxArray(51,IP) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,IP) = Xqd*AuxArray(52,IP) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,IP) = Xqd*AuxArray(53,IP) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,IP) = Xqd*AuxArray(54,IP) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,IP) = Xqd*AuxArray(55,IP) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,IP) = Xqd*AuxArray(56,IP) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,IP) = Yqd*AuxArray(51,IP) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,IP) = Zqd*AuxArray(51,IP) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,IP) = Yqd*AuxArray(53,IP) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,IP) = Yqd*AuxArray(54,IP) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,IP) = Yqd*AuxArray(55,IP) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,IP) = Yqd*AuxArray(56,IP) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,IP) = Zqd*AuxArray(56,IP) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPU7D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:10,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(  120,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
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
   iP = (iPassQ-1)*nPrimQ*nPrimP
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
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
     Auxarray(1,IP) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
     TMParray1(1, 5) = PREF*RJ000( 4)
     TMParray1(1, 6) = PREF*RJ000( 5)
     TMParray1(1, 7) = PREF*RJ000( 6)
     TMParray1(1, 8) = PREF*RJ000( 7)
     AuxArray(2,IP) = Xqd*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqd*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqd*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,IP) = Xqd*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqd*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqd*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqd*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqd*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqd*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expQ*(AuxArray(2,IP) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(3,IP) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(4,IP) + alphaQ*TmpArray2(4,2))
     AuxArray(11,IP) = Xqd*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Yqd*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zqd*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xqd*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xqd*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xqd*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Yqd*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zqd*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Yqd*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zqd*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(AuxArray(5,IP) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(AuxArray(8,IP) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(AuxArray(10,IP) + alphaQ*TmpArray3(10,2))
     AuxArray(21,IP) = Xqd*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Yqd*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zqd*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xqd*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Yqd*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xqd*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xqd*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xqd*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xqd*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xqd*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Yqd*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zqd*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Yqd*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Yqd*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zqd*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(AuxArray(11,IP) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(AuxArray(14,IP) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(AuxArray(16,IP) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(AuxArray(17,IP) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(AuxArray(19,IP) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(AuxArray(20,IP) + alphaQ*TmpArray4(20,2))
     AuxArray(36,IP) = Xqd*AuxArray(21,IP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IP) = Yqd*AuxArray(21,IP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IP) = Zqd*AuxArray(21,IP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IP) = Xqd*AuxArray(24,IP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IP) = Yqd*AuxArray(23,IP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IP) = Xqd*AuxArray(26,IP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IP) = Xqd*AuxArray(27,IP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IP) = Zqd*AuxArray(24,IP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IP) = Yqd*AuxArray(26,IP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IP) = Xqd*AuxArray(30,IP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IP) = Xqd*AuxArray(31,IP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IP) = Xqd*AuxArray(32,IP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IP) = Xqd*AuxArray(33,IP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IP) = Xqd*AuxArray(34,IP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IP) = Xqd*AuxArray(35,IP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IP) = Yqd*AuxArray(31,IP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IP) = Zqd*AuxArray(31,IP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IP) = Yqd*AuxArray(33,IP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IP) = Yqd*AuxArray(34,IP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IP) = Yqd*AuxArray(35,IP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IP) = Zqd*AuxArray(35,IP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     TwoTerms(1) = inv2expQ*(AuxArray(21,IP) + alphaQ*TmpArray5(21,2))
     TwoTerms(2) = inv2expQ*(AuxArray(24,IP) + alphaQ*TmpArray5(24,2))
     TwoTerms(3) = inv2expQ*(AuxArray(26,IP) + alphaQ*TmpArray5(26,2))
     TwoTerms(4) = inv2expQ*(AuxArray(27,IP) + alphaQ*TmpArray5(27,2))
     TwoTerms(5) = inv2expQ*(AuxArray(30,IP) + alphaQ*TmpArray5(30,2))
     TwoTerms(6) = inv2expQ*(AuxArray(31,IP) + alphaQ*TmpArray5(31,2))
     TwoTerms(7) = inv2expQ*(AuxArray(33,IP) + alphaQ*TmpArray5(33,2))
     TwoTerms(8) = inv2expQ*(AuxArray(34,IP) + alphaQ*TmpArray5(34,2))
     TwoTerms(9) = inv2expQ*(AuxArray(35,IP) + alphaQ*TmpArray5(35,2))
     AuxArray(57,IP) = Xqd*AuxArray(36,IP) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,IP) = Yqd*AuxArray(36,IP) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,IP) = Zqd*AuxArray(36,IP) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,IP) = Xqd*AuxArray(39,IP) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,IP) = Yqd*AuxArray(38,IP) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,IP) = Xqd*AuxArray(41,IP) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,IP) = Xqd*AuxArray(42,IP) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,IP) = Zqd*AuxArray(39,IP) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,IP) = Yqd*AuxArray(41,IP) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,IP) = Xqd*AuxArray(45,IP) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,IP) = Xqd*AuxArray(46,IP) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,IP) = Zqd*AuxArray(42,IP) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,IP) = Xqd*AuxArray(48,IP) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,IP) = Yqd*AuxArray(45,IP) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,IP) = Xqd*AuxArray(50,IP) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,IP) = Xqd*AuxArray(51,IP) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,IP) = Xqd*AuxArray(52,IP) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,IP) = Xqd*AuxArray(53,IP) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,IP) = Xqd*AuxArray(54,IP) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,IP) = Xqd*AuxArray(55,IP) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,IP) = Xqd*AuxArray(56,IP) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,IP) = Yqd*AuxArray(51,IP) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,IP) = Zqd*AuxArray(51,IP) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,IP) = Yqd*AuxArray(53,IP) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,IP) = Yqd*AuxArray(54,IP) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,IP) = Yqd*AuxArray(55,IP) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,IP) = Yqd*AuxArray(56,IP) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,IP) = Zqd*AuxArray(56,IP) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
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
     AuxArray(85,IP) = Xqd*AuxArray(57,IP) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     AuxArray(86,IP) = Yqd*AuxArray(57,IP) + alphaYpq*TmpArray7(57,2)
     AuxArray(87,IP) = Zqd*AuxArray(57,IP) + alphaZpq*TmpArray7(57,2)
     AuxArray(88,IP) = Xqd*AuxArray(60,IP) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     AuxArray(89,IP) = Yqd*AuxArray(59,IP) + alphaYpq*TmpArray7(59,2)
     AuxArray(90,IP) = Xqd*AuxArray(62,IP) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     AuxArray(91,IP) = Xqd*AuxArray(63,IP) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     AuxArray(92,IP) = Zqd*AuxArray(60,IP) + alphaZpq*TmpArray7(60,2)
     AuxArray(93,IP) = Yqd*AuxArray(62,IP) + alphaYpq*TmpArray7(62,2)
     AuxArray(94,IP) = Xqd*AuxArray(66,IP) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     AuxArray(95,IP) = Xqd*AuxArray(67,IP) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     AuxArray(96,IP) = Zqd*AuxArray(63,IP) + alphaZpq*TmpArray7(63,2)
     AuxArray(97,IP) = Xqd*AuxArray(69,IP) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     AuxArray(98,IP) = Yqd*AuxArray(66,IP) + alphaYpq*TmpArray7(66,2)
     AuxArray(99,IP) = Xqd*AuxArray(71,IP) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     AuxArray(100,IP) = Xqd*AuxArray(72,IP) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     AuxArray(101,IP) = Zqd*AuxArray(67,IP) + alphaZpq*TmpArray7(67,2)
     AuxArray(102,IP) = Xqd*AuxArray(74,IP) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     AuxArray(103,IP) = Xqd*AuxArray(75,IP) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     AuxArray(104,IP) = Yqd*AuxArray(71,IP) + alphaYpq*TmpArray7(71,2)
     AuxArray(105,IP) = Xqd*AuxArray(77,IP) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     AuxArray(106,IP) = Xqd*AuxArray(78,IP) + alphaXpq*TmpArray7(78,2)
     AuxArray(107,IP) = Xqd*AuxArray(79,IP) + alphaXpq*TmpArray7(79,2)
     AuxArray(108,IP) = Xqd*AuxArray(80,IP) + alphaXpq*TmpArray7(80,2)
     AuxArray(109,IP) = Xqd*AuxArray(81,IP) + alphaXpq*TmpArray7(81,2)
     AuxArray(110,IP) = Xqd*AuxArray(82,IP) + alphaXpq*TmpArray7(82,2)
     AuxArray(111,IP) = Xqd*AuxArray(83,IP) + alphaXpq*TmpArray7(83,2)
     AuxArray(112,IP) = Xqd*AuxArray(84,IP) + alphaXpq*TmpArray7(84,2)
     AuxArray(113,IP) = Yqd*AuxArray(78,IP) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     AuxArray(114,IP) = Zqd*AuxArray(78,IP) + alphaZpq*TmpArray7(78,2)
     AuxArray(115,IP) = Yqd*AuxArray(80,IP) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     AuxArray(116,IP) = Yqd*AuxArray(81,IP) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     AuxArray(117,IP) = Yqd*AuxArray(82,IP) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     AuxArray(118,IP) = Yqd*AuxArray(83,IP) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     AuxArray(119,IP) = Yqd*AuxArray(84,IP) + alphaYpq*TmpArray7(84,2)
     AuxArray(120,IP) = Zqd*AuxArray(84,IP) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPU8D(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:11,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Dcenter(3)
  real(realk),intent(inout) :: AUXarray(  165,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
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
   iP = (iPassQ-1)*nPrimQ*nPrimP
   Dx = -Dcenter(1)
   Dy = -Dcenter(2)
   Dz = -Dcenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
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
     Auxarray(1,IP) = PREF*RJ000(0)
     TMParray1(1, 2) = PREF*RJ000( 1)
     TMParray1(1, 3) = PREF*RJ000( 2)
     TMParray1(1, 4) = PREF*RJ000( 3)
     TMParray1(1, 5) = PREF*RJ000( 4)
     TMParray1(1, 6) = PREF*RJ000( 5)
     TMParray1(1, 7) = PREF*RJ000( 6)
     TMParray1(1, 8) = PREF*RJ000( 7)
     TMParray1(1, 9) = PREF*RJ000( 8)
     AuxArray(2,IP) = Xqd*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqd*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqd*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,IP) = Xqd*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqd*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqd*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqd*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqd*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqd*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expQ*(AuxArray(2,IP) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(AuxArray(3,IP) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(AuxArray(4,IP) + alphaQ*TmpArray2(4,2))
     AuxArray(11,IP) = Xqd*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Yqd*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zqd*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xqd*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xqd*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xqd*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Yqd*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zqd*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Yqd*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zqd*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(AuxArray(5,IP) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(AuxArray(8,IP) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(AuxArray(10,IP) + alphaQ*TmpArray3(10,2))
     AuxArray(21,IP) = Xqd*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Yqd*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zqd*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xqd*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Yqd*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xqd*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xqd*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xqd*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xqd*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xqd*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Yqd*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zqd*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Yqd*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Yqd*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zqd*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(AuxArray(11,IP) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(AuxArray(14,IP) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(AuxArray(16,IP) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(AuxArray(17,IP) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(AuxArray(19,IP) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(AuxArray(20,IP) + alphaQ*TmpArray4(20,2))
     AuxArray(36,IP) = Xqd*AuxArray(21,IP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IP) = Yqd*AuxArray(21,IP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IP) = Zqd*AuxArray(21,IP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IP) = Xqd*AuxArray(24,IP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IP) = Yqd*AuxArray(23,IP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IP) = Xqd*AuxArray(26,IP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IP) = Xqd*AuxArray(27,IP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IP) = Zqd*AuxArray(24,IP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IP) = Yqd*AuxArray(26,IP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IP) = Xqd*AuxArray(30,IP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IP) = Xqd*AuxArray(31,IP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IP) = Xqd*AuxArray(32,IP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IP) = Xqd*AuxArray(33,IP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IP) = Xqd*AuxArray(34,IP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IP) = Xqd*AuxArray(35,IP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IP) = Yqd*AuxArray(31,IP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IP) = Zqd*AuxArray(31,IP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IP) = Yqd*AuxArray(33,IP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IP) = Yqd*AuxArray(34,IP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IP) = Yqd*AuxArray(35,IP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IP) = Zqd*AuxArray(35,IP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     TwoTerms(1) = inv2expQ*(AuxArray(21,IP) + alphaQ*TmpArray5(21,2))
     TwoTerms(2) = inv2expQ*(AuxArray(24,IP) + alphaQ*TmpArray5(24,2))
     TwoTerms(3) = inv2expQ*(AuxArray(26,IP) + alphaQ*TmpArray5(26,2))
     TwoTerms(4) = inv2expQ*(AuxArray(27,IP) + alphaQ*TmpArray5(27,2))
     TwoTerms(5) = inv2expQ*(AuxArray(30,IP) + alphaQ*TmpArray5(30,2))
     TwoTerms(6) = inv2expQ*(AuxArray(31,IP) + alphaQ*TmpArray5(31,2))
     TwoTerms(7) = inv2expQ*(AuxArray(33,IP) + alphaQ*TmpArray5(33,2))
     TwoTerms(8) = inv2expQ*(AuxArray(34,IP) + alphaQ*TmpArray5(34,2))
     TwoTerms(9) = inv2expQ*(AuxArray(35,IP) + alphaQ*TmpArray5(35,2))
     AuxArray(57,IP) = Xqd*AuxArray(36,IP) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,IP) = Yqd*AuxArray(36,IP) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,IP) = Zqd*AuxArray(36,IP) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,IP) = Xqd*AuxArray(39,IP) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,IP) = Yqd*AuxArray(38,IP) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,IP) = Xqd*AuxArray(41,IP) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,IP) = Xqd*AuxArray(42,IP) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,IP) = Zqd*AuxArray(39,IP) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,IP) = Yqd*AuxArray(41,IP) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,IP) = Xqd*AuxArray(45,IP) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,IP) = Xqd*AuxArray(46,IP) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,IP) = Zqd*AuxArray(42,IP) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,IP) = Xqd*AuxArray(48,IP) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,IP) = Yqd*AuxArray(45,IP) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,IP) = Xqd*AuxArray(50,IP) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,IP) = Xqd*AuxArray(51,IP) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,IP) = Xqd*AuxArray(52,IP) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,IP) = Xqd*AuxArray(53,IP) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,IP) = Xqd*AuxArray(54,IP) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,IP) = Xqd*AuxArray(55,IP) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,IP) = Xqd*AuxArray(56,IP) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,IP) = Yqd*AuxArray(51,IP) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,IP) = Zqd*AuxArray(51,IP) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,IP) = Yqd*AuxArray(53,IP) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,IP) = Yqd*AuxArray(54,IP) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,IP) = Yqd*AuxArray(55,IP) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,IP) = Yqd*AuxArray(56,IP) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,IP) = Zqd*AuxArray(56,IP) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
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
     AuxArray(85,IP) = Xqd*AuxArray(57,IP) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     AuxArray(86,IP) = Yqd*AuxArray(57,IP) + alphaYpq*TmpArray7(57,2)
     AuxArray(87,IP) = Zqd*AuxArray(57,IP) + alphaZpq*TmpArray7(57,2)
     AuxArray(88,IP) = Xqd*AuxArray(60,IP) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     AuxArray(89,IP) = Yqd*AuxArray(59,IP) + alphaYpq*TmpArray7(59,2)
     AuxArray(90,IP) = Xqd*AuxArray(62,IP) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     AuxArray(91,IP) = Xqd*AuxArray(63,IP) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     AuxArray(92,IP) = Zqd*AuxArray(60,IP) + alphaZpq*TmpArray7(60,2)
     AuxArray(93,IP) = Yqd*AuxArray(62,IP) + alphaYpq*TmpArray7(62,2)
     AuxArray(94,IP) = Xqd*AuxArray(66,IP) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     AuxArray(95,IP) = Xqd*AuxArray(67,IP) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     AuxArray(96,IP) = Zqd*AuxArray(63,IP) + alphaZpq*TmpArray7(63,2)
     AuxArray(97,IP) = Xqd*AuxArray(69,IP) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     AuxArray(98,IP) = Yqd*AuxArray(66,IP) + alphaYpq*TmpArray7(66,2)
     AuxArray(99,IP) = Xqd*AuxArray(71,IP) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     AuxArray(100,IP) = Xqd*AuxArray(72,IP) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     AuxArray(101,IP) = Zqd*AuxArray(67,IP) + alphaZpq*TmpArray7(67,2)
     AuxArray(102,IP) = Xqd*AuxArray(74,IP) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     AuxArray(103,IP) = Xqd*AuxArray(75,IP) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     AuxArray(104,IP) = Yqd*AuxArray(71,IP) + alphaYpq*TmpArray7(71,2)
     AuxArray(105,IP) = Xqd*AuxArray(77,IP) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     AuxArray(106,IP) = Xqd*AuxArray(78,IP) + alphaXpq*TmpArray7(78,2)
     AuxArray(107,IP) = Xqd*AuxArray(79,IP) + alphaXpq*TmpArray7(79,2)
     AuxArray(108,IP) = Xqd*AuxArray(80,IP) + alphaXpq*TmpArray7(80,2)
     AuxArray(109,IP) = Xqd*AuxArray(81,IP) + alphaXpq*TmpArray7(81,2)
     AuxArray(110,IP) = Xqd*AuxArray(82,IP) + alphaXpq*TmpArray7(82,2)
     AuxArray(111,IP) = Xqd*AuxArray(83,IP) + alphaXpq*TmpArray7(83,2)
     AuxArray(112,IP) = Xqd*AuxArray(84,IP) + alphaXpq*TmpArray7(84,2)
     AuxArray(113,IP) = Yqd*AuxArray(78,IP) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     AuxArray(114,IP) = Zqd*AuxArray(78,IP) + alphaZpq*TmpArray7(78,2)
     AuxArray(115,IP) = Yqd*AuxArray(80,IP) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     AuxArray(116,IP) = Yqd*AuxArray(81,IP) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     AuxArray(117,IP) = Yqd*AuxArray(82,IP) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     AuxArray(118,IP) = Yqd*AuxArray(83,IP) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     AuxArray(119,IP) = Yqd*AuxArray(84,IP) + alphaYpq*TmpArray7(84,2)
     AuxArray(120,IP) = Zqd*AuxArray(84,IP) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
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
     TwoTerms(1) = inv2expQ*(AuxArray(57,IP) + alphaQ*TmpArray7(57,2))
     TwoTerms(2) = inv2expQ*(AuxArray(60,IP) + alphaQ*TmpArray7(60,2))
     TwoTerms(3) = inv2expQ*(AuxArray(62,IP) + alphaQ*TmpArray7(62,2))
     TwoTerms(4) = inv2expQ*(AuxArray(63,IP) + alphaQ*TmpArray7(63,2))
     TwoTerms(5) = inv2expQ*(AuxArray(66,IP) + alphaQ*TmpArray7(66,2))
     TwoTerms(6) = inv2expQ*(AuxArray(67,IP) + alphaQ*TmpArray7(67,2))
     TwoTerms(7) = inv2expQ*(AuxArray(69,IP) + alphaQ*TmpArray7(69,2))
     TwoTerms(8) = inv2expQ*(AuxArray(71,IP) + alphaQ*TmpArray7(71,2))
     TwoTerms(9) = inv2expQ*(AuxArray(72,IP) + alphaQ*TmpArray7(72,2))
     TwoTerms(10) = inv2expQ*(AuxArray(74,IP) + alphaQ*TmpArray7(74,2))
     TwoTerms(11) = inv2expQ*(AuxArray(75,IP) + alphaQ*TmpArray7(75,2))
     TwoTerms(12) = inv2expQ*(AuxArray(77,IP) + alphaQ*TmpArray7(77,2))
     TwoTerms(13) = inv2expQ*(AuxArray(78,IP) + alphaQ*TmpArray7(78,2))
     TwoTerms(14) = inv2expQ*(AuxArray(80,IP) + alphaQ*TmpArray7(80,2))
     TwoTerms(15) = inv2expQ*(AuxArray(81,IP) + alphaQ*TmpArray7(81,2))
     TwoTerms(16) = inv2expQ*(AuxArray(82,IP) + alphaQ*TmpArray7(82,2))
     TwoTerms(17) = inv2expQ*(AuxArray(83,IP) + alphaQ*TmpArray7(83,2))
     TwoTerms(18) = inv2expQ*(AuxArray(84,IP) + alphaQ*TmpArray7(84,2))
     AuxArray(121,IP) = Xqd*AuxArray(85,IP) + alphaXpq*TmpArray8(85,2) + 7*TwoTerms(1)
     AuxArray(122,IP) = Yqd*AuxArray(85,IP) + alphaYpq*TmpArray8(85,2)
     AuxArray(123,IP) = Zqd*AuxArray(85,IP) + alphaZpq*TmpArray8(85,2)
     AuxArray(124,IP) = Xqd*AuxArray(88,IP) + alphaXpq*TmpArray8(88,2) + 5*TwoTerms(2)
     AuxArray(125,IP) = Yqd*AuxArray(87,IP) + alphaYpq*TmpArray8(87,2)
     AuxArray(126,IP) = Xqd*AuxArray(90,IP) + alphaXpq*TmpArray8(90,2) + 5*TwoTerms(3)
     AuxArray(127,IP) = Xqd*AuxArray(91,IP) + alphaXpq*TmpArray8(91,2) + 4*TwoTerms(4)
     AuxArray(128,IP) = Zqd*AuxArray(88,IP) + alphaZpq*TmpArray8(88,2)
     AuxArray(129,IP) = Yqd*AuxArray(90,IP) + alphaYpq*TmpArray8(90,2)
     AuxArray(130,IP) = Xqd*AuxArray(94,IP) + alphaXpq*TmpArray8(94,2) + 4*TwoTerms(5)
     AuxArray(131,IP) = Xqd*AuxArray(95,IP) + alphaXpq*TmpArray8(95,2) + 3*TwoTerms(6)
     AuxArray(132,IP) = Zqd*AuxArray(91,IP) + alphaZpq*TmpArray8(91,2)
     AuxArray(133,IP) = Xqd*AuxArray(97,IP) + alphaXpq*TmpArray8(97,2) + 3*TwoTerms(7)
     AuxArray(134,IP) = Yqd*AuxArray(94,IP) + alphaYpq*TmpArray8(94,2)
     AuxArray(135,IP) = Xqd*AuxArray(99,IP) + alphaXpq*TmpArray8(99,2) + 3*TwoTerms(8)
     AuxArray(136,IP) = Xqd*AuxArray(100,IP) + alphaXpq*TmpArray8(100,2) + 2*TwoTerms(9)
     AuxArray(137,IP) = Zqd*AuxArray(95,IP) + alphaZpq*TmpArray8(95,2)
     AuxArray(138,IP) = Xqd*AuxArray(102,IP) + alphaXpq*TmpArray8(102,2) + 2*TwoTerms(10)
     AuxArray(139,IP) = Xqd*AuxArray(103,IP) + alphaXpq*TmpArray8(103,2) + 2*TwoTerms(11)
     AuxArray(140,IP) = Yqd*AuxArray(99,IP) + alphaYpq*TmpArray8(99,2)
     AuxArray(141,IP) = Xqd*AuxArray(105,IP) + alphaXpq*TmpArray8(105,2) + 2*TwoTerms(12)
     AuxArray(142,IP) = Xqd*AuxArray(106,IP) + alphaXpq*TmpArray8(106,2) + TwoTerms(13)
     AuxArray(143,IP) = Zqd*AuxArray(100,IP) + alphaZpq*TmpArray8(100,2)
     AuxArray(144,IP) = Xqd*AuxArray(108,IP) + alphaXpq*TmpArray8(108,2) + TwoTerms(14)
     AuxArray(145,IP) = Xqd*AuxArray(109,IP) + alphaXpq*TmpArray8(109,2) + TwoTerms(15)
     AuxArray(146,IP) = Xqd*AuxArray(110,IP) + alphaXpq*TmpArray8(110,2) + TwoTerms(16)
     AuxArray(147,IP) = Yqd*AuxArray(105,IP) + alphaYpq*TmpArray8(105,2)
     AuxArray(148,IP) = Xqd*AuxArray(112,IP) + alphaXpq*TmpArray8(112,2) + TwoTerms(18)
     AuxArray(149,IP) = Xqd*AuxArray(113,IP) + alphaXpq*TmpArray8(113,2)
     AuxArray(150,IP) = Xqd*AuxArray(114,IP) + alphaXpq*TmpArray8(114,2)
     AuxArray(151,IP) = Xqd*AuxArray(115,IP) + alphaXpq*TmpArray8(115,2)
     AuxArray(152,IP) = Xqd*AuxArray(116,IP) + alphaXpq*TmpArray8(116,2)
     AuxArray(153,IP) = Xqd*AuxArray(117,IP) + alphaXpq*TmpArray8(117,2)
     AuxArray(154,IP) = Xqd*AuxArray(118,IP) + alphaXpq*TmpArray8(118,2)
     AuxArray(155,IP) = Xqd*AuxArray(119,IP) + alphaXpq*TmpArray8(119,2)
     AuxArray(156,IP) = Xqd*AuxArray(120,IP) + alphaXpq*TmpArray8(120,2)
     AuxArray(157,IP) = Yqd*AuxArray(113,IP) + alphaYpq*TmpArray8(113,2) + 7*TwoTerms(13)
     AuxArray(158,IP) = Zqd*AuxArray(113,IP) + alphaZpq*TmpArray8(113,2)
     AuxArray(159,IP) = Yqd*AuxArray(115,IP) + alphaYpq*TmpArray8(115,2) + 5*TwoTerms(14)
     AuxArray(160,IP) = Yqd*AuxArray(116,IP) + alphaYpq*TmpArray8(116,2) + 4*TwoTerms(15)
     AuxArray(161,IP) = Yqd*AuxArray(117,IP) + alphaYpq*TmpArray8(117,2) + 3*TwoTerms(16)
     AuxArray(162,IP) = Yqd*AuxArray(118,IP) + alphaYpq*TmpArray8(118,2) + 2*TwoTerms(17)
     AuxArray(163,IP) = Yqd*AuxArray(119,IP) + alphaYpq*TmpArray8(119,2) + TwoTerms(18)
     AuxArray(164,IP) = Yqd*AuxArray(120,IP) + alphaYpq*TmpArray8(120,2)
     AuxArray(165,IP) = Zqd*AuxArray(120,IP) + alphaZpq*TmpArray8(120,2) + 7*TwoTerms(18)
    ENDDO
   ENDDO
  ENDDO
 end subroutine
end module
