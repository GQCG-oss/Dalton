MODULE AGC_CPU_OBS_VERTICALRECURRENCEMODASegQ
 use IchorPrecisionModule
  
 CONTAINS

subroutine VerticalRecurrenceCPUSegQ0(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,&
         & Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  REAL(REALK),intent(in) :: integralPrefactor(nprimQ,nPrimP)
  REAL(REALK),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  REAL(REALK),intent(in) :: TABFJW(0:3,0:1200)
  REAL(REALK),intent(in) :: QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(inout) :: AUXarray(nPrimP,nPasses)
  !local variables
  REAL(REALK),PARAMETER :: D2JP36=  3.6000000000000000E+01_realk
  real(realk),parameter :: D2=2.0E0_realk
  REAL(REALK),PARAMETER :: D05 =0.5E0_realk,D1=1E0_realk
  REAL(REALK),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk
  Real(realk),parameter :: D12 = 12E0_realk, TENTH = 0.01E0_realk
  REAL(REALK),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
  REAL(REALK),PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk
  REAL(REALK),PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
  REAL(REALK),PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
  REAL(REALK),PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
  REAL(REALK),PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
  Real(realk),parameter :: PI=3.14159265358979323846E0_realk
  REAL(REALK),PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
  REAL(REALK),PARAMETER :: SQRPIH = SQRTPI/D2
  REAL(REALK),PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
!  REAL(REALK),PARAMETER :: SMALL = 1E-15_realk
  Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA,WVAL,Pexpfac
  Real(realk) :: W2,W3,PX,PY,PZ,PQX,PQY,PQZ,squaredDistance,RJ000
  Integer :: IPNT,iPassQ,iPrimP,iPrimQ,iPQ
  DO iPassQ = 1,nPasses
   DO iPrimP=1, nPrimP
    AUXarray(iPrimP,iPassQ)=0.0E0_realk
    Pexpfac = PpreExpFac(iPrimP)
    px = Pcent(1,iPrimP)
    py = Pcent(2,iPrimP)
    pz = Pcent(3,iPrimP)
    DO iPrimQ=1, nPrimQ
     pqx = px - Qcent(1,iPrimQ,iPassQ)
     pqy = py - Qcent(2,iPrimQ,iPassQ)
     pqz = pz - Qcent(3,iPrimQ,iPassQ)
     squaredDistance = pqx*pqx+pqy*pqy+pqz*pqz
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
     AUXarray(iPrimP,iPassQ)=AUXarray(iPrimP,iPassQ) + integralPrefactor(iPrimQ,iPrimP)*&
          & QpreExpFac(iPrimQ,iPassQ)*Pexpfac*RJ000
    enddo
   enddo
  enddo
end subroutine VerticalRecurrenceCPUSegQ0

subroutine VerticalRecurrenceCPUSegQ1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,&
         & QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:4,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Acenter(3)
  real(realk),intent(inout) :: AUXarray(4,nPrimP,nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
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
  !ThetaAux(n,1,0,0) = Xpa*ThetaAux(n,0,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,0,0,0)
  !i = 0 last 2 term vanish
  !We include scaling of RJ000 
  Ax = -Acenter(1)
  Ay = -Acenter(2)
  Az = -Acenter(3)
  DO iPassQ = 1,nPasses
   DO iPrimP=1, nPrimP
    AUXarray(1,iPrimP,iPassQ)=0.0E0_realk
    AUXarray(2,iPrimP,iPassQ)=0.0E0_realk
    AUXarray(3,iPrimP,iPassQ)=0.0E0_realk
    AUXarray(4,iPrimP,iPassQ)=0.0E0_realk
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    Xpa = Pcent(1,iPrimP) + Ax
    Ypa = Pcent(2,iPrimP) + Ay
    Zpa = Pcent(3,iPrimP) + Az
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
     AUXarray(1,iPrimP,iPassQ) = AUXarray(1,iPrimP,iPassQ) + TMP1
     AUXarray(2,iPrimP,iPassQ) = AUXarray(2,iPrimP,iPassQ) + Xpa*TMP1 + alphaXpq*TMP2
     AUXarray(3,iPrimP,iPassQ) = AUXarray(3,iPrimP,iPassQ) + Ypa*TMP1 + alphaYpq*TMP2
     AUXarray(4,iPrimP,iPassQ) = AUXarray(4,iPrimP,iPassQ) + Zpa*TMP1 + alphaZpq*TMP2
    enddo
   enddo
  enddo
end subroutine

subroutine VerticalRecurrenceCPUSegQ2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 5,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Acenter(3)
  real(realk),intent(inout) :: AUXarray(   10,nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(    4)
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
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
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Ax = -Acenter(1)
  Ay = -Acenter(2)
  Az = -Acenter(3)
  DO iPassQ = 1,nPasses
   iP = (iPassQ-1)*nPrimP
   DO iPrimP=1, nPrimP
    iP = iP + 1
    DO iTUV=1,   10
     AUXarray(iTUV,iP)=0.0E0_realk
    ENDDO
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP) + Ax
    Ypa = Pcent(2,iPrimP) + Ay
    Zpa = Pcent(3,iPrimP) + Az
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
     TMPAuxArray(2) = Xpa*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Ypa*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zpa*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpa*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypa*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpa*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(1) + alphaP*TmpArray1(1,2))
     do iTUV = 1,    4
      AuxArray(iTUV,IP) = AuxArray(iTUV,IP) + TMPAuxarray(iTUV)
     enddo
     AuxArray(5,IP) = AuxArray(5,IP) + Xpa*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = AuxArray(6,IP) + Xpa*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = AuxArray(7,IP) + Xpa*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = AuxArray(8,IP) + Ypa*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = AuxArray(9,IP) + Ypa*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = AuxArray(10,IP) + Zpa*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSegQ3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 6,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Acenter(3)
  real(realk),intent(inout) :: AUXarray(   20,nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   10)
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
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
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Ax = -Acenter(1)
  Ay = -Acenter(2)
  Az = -Acenter(3)
  DO iPassQ = 1,nPasses
   iP = (iPassQ-1)*nPrimP
   DO iPrimP=1, nPrimP
    iP = iP + 1
    DO iTUV=1,   20
     AUXarray(iTUV,iP)=0.0E0_realk
    ENDDO
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP) + Ax
    Ypa = Pcent(2,iPrimP) + Ay
    Zpa = Pcent(3,iPrimP) + Az
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
     do iTUV = 1,   10
      AuxArray(iTUV,IP) = AuxArray(iTUV,IP) + TMPAuxarray(iTUV)
     enddo
     AuxArray(11,IP) = AuxArray(11,IP) + Xpa*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = AuxArray(12,IP) + Ypa*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = AuxArray(13,IP) + Zpa*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = AuxArray(14,IP) + Xpa*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = AuxArray(15,IP) + Xpa*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = AuxArray(16,IP) + Xpa*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = AuxArray(17,IP) + Ypa*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = AuxArray(18,IP) + Zpa*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = AuxArray(19,IP) + Ypa*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = AuxArray(20,IP) + Zpa*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSegQ4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 7,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Acenter(3)
  real(realk),intent(inout) :: AUXarray(   35,nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   20)
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
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
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Ax = -Acenter(1)
  Ay = -Acenter(2)
  Az = -Acenter(3)
  DO iPassQ = 1,nPasses
   iP = (iPassQ-1)*nPrimP
   DO iPrimP=1, nPrimP
    iP = iP + 1
    DO iTUV=1,   35
     AUXarray(iTUV,iP)=0.0E0_realk
    ENDDO
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP) + Ax
    Ypa = Pcent(2,iPrimP) + Ay
    Zpa = Pcent(3,iPrimP) + Az
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
     do iTUV = 1,   20
      AuxArray(iTUV,IP) = AuxArray(iTUV,IP) + TMPAuxarray(iTUV)
     enddo
     AuxArray(21,IP) = AuxArray(21,IP) + Xpa*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = AuxArray(22,IP) + Ypa*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = AuxArray(23,IP) + Zpa*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = AuxArray(24,IP) + Xpa*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = AuxArray(25,IP) + Ypa*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = AuxArray(26,IP) + Xpa*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = AuxArray(27,IP) + Xpa*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = AuxArray(28,IP) + Xpa*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = AuxArray(29,IP) + Xpa*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = AuxArray(30,IP) + Xpa*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = AuxArray(31,IP) + Ypa*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = AuxArray(32,IP) + Zpa*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = AuxArray(33,IP) + Ypa*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = AuxArray(34,IP) + Ypa*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = AuxArray(35,IP) + Zpa*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSegQ5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 8,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Acenter(3)
  real(realk),intent(inout) :: AUXarray(   56,nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   35)
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
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
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Ax = -Acenter(1)
  Ay = -Acenter(2)
  Az = -Acenter(3)
  DO iPassQ = 1,nPasses
   iP = (iPassQ-1)*nPrimP
   DO iPrimP=1, nPrimP
    iP = iP + 1
    DO iTUV=1,   56
     AUXarray(iTUV,iP)=0.0E0_realk
    ENDDO
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP) + Ax
    Ypa = Pcent(2,iPrimP) + Ay
    Zpa = Pcent(3,iPrimP) + Az
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
     tmpArray2(2,5) = Xpa*tmpArray1(1,5) + alphaXpq*TmpArray1(1,6)
     tmpArray2(3,5) = Ypa*tmpArray1(1,5) + alphaYpq*TmpArray1(1,6)
     tmpArray2(4,5) = Zpa*tmpArray1(1,5) + alphaZpq*TmpArray1(1,6)
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
     TwoTerms(1) = inv2expP*(TmpArray1(1,4) + alphaP*TmpArray1(1,5))
     tmpArray3(5,4) = Xpa*tmpArray2(2,4) + alphaXpq*TmpArray2(2,5) + TwoTerms(1)
     tmpArray3(6,4) = Xpa*tmpArray2(3,4) + alphaXpq*TmpArray2(3,5)
     tmpArray3(7,4) = Xpa*tmpArray2(4,4) + alphaXpq*TmpArray2(4,5)
     tmpArray3(8,4) = Ypa*tmpArray2(3,4) + alphaYpq*TmpArray2(3,5) + TwoTerms(1)
     tmpArray3(9,4) = Ypa*tmpArray2(4,4) + alphaYpq*TmpArray2(4,5)
     tmpArray3(10,4) = Zpa*tmpArray2(4,4) + alphaZpq*TmpArray2(4,5) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(TmpArray2(2,3) + alphaP*TmpArray2(2,4))
     TwoTerms(2) = inv2expP*(TmpArray2(3,3) + alphaP*TmpArray2(3,4))
     TwoTerms(3) = inv2expP*(TmpArray2(4,3) + alphaP*TmpArray2(4,4))
     tmpArray4(11,3) = Xpa*tmpArray3(5,3) + alphaXpq*TmpArray3(5,4) + 2*TwoTerms(1)
     tmpArray4(12,3) = Ypa*tmpArray3(5,3) + alphaYpq*TmpArray3(5,4)
     tmpArray4(13,3) = Zpa*tmpArray3(5,3) + alphaZpq*TmpArray3(5,4)
     tmpArray4(14,3) = Xpa*tmpArray3(8,3) + alphaXpq*TmpArray3(8,4)
     tmpArray4(15,3) = Xpa*tmpArray3(9,3) + alphaXpq*TmpArray3(9,4)
     tmpArray4(16,3) = Xpa*tmpArray3(10,3) + alphaXpq*TmpArray3(10,4)
     tmpArray4(17,3) = Ypa*tmpArray3(8,3) + alphaYpq*TmpArray3(8,4) + 2*TwoTerms(2)
     tmpArray4(18,3) = Zpa*tmpArray3(8,3) + alphaZpq*TmpArray3(8,4)
     tmpArray4(19,3) = Ypa*tmpArray3(10,3) + alphaYpq*TmpArray3(10,4)
     tmpArray4(20,3) = Zpa*tmpArray3(10,3) + alphaZpq*TmpArray3(10,4) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(10) + alphaP*TmpArray3(10,2))
     TMPAuxArray(21) = Xpa*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Ypa*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zpa*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xpa*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Ypa*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xpa*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xpa*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xpa*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xpa*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xpa*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Ypa*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zpa*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Ypa*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Ypa*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zpa*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,2) + alphaP*TmpArray3(5,3))
     TwoTerms(2) = inv2expP*(TmpArray3(8,2) + alphaP*TmpArray3(8,3))
     TwoTerms(3) = inv2expP*(TmpArray3(10,2) + alphaP*TmpArray3(10,3))
     tmpArray5(21,2) = Xpa*tmpArray4(11,2) + alphaXpq*TmpArray4(11,3) + 3*TwoTerms(1)
     tmpArray5(22,2) = Ypa*tmpArray4(11,2) + alphaYpq*TmpArray4(11,3)
     tmpArray5(23,2) = Zpa*tmpArray4(11,2) + alphaZpq*TmpArray4(11,3)
     tmpArray5(24,2) = Xpa*tmpArray4(14,2) + alphaXpq*TmpArray4(14,3) + TwoTerms(2)
     tmpArray5(25,2) = Ypa*tmpArray4(13,2) + alphaYpq*TmpArray4(13,3)
     tmpArray5(26,2) = Xpa*tmpArray4(16,2) + alphaXpq*TmpArray4(16,3) + TwoTerms(3)
     tmpArray5(27,2) = Xpa*tmpArray4(17,2) + alphaXpq*TmpArray4(17,3)
     tmpArray5(28,2) = Xpa*tmpArray4(18,2) + alphaXpq*TmpArray4(18,3)
     tmpArray5(29,2) = Xpa*tmpArray4(19,2) + alphaXpq*TmpArray4(19,3)
     tmpArray5(30,2) = Xpa*tmpArray4(20,2) + alphaXpq*TmpArray4(20,3)
     tmpArray5(31,2) = Ypa*tmpArray4(17,2) + alphaYpq*TmpArray4(17,3) + 3*TwoTerms(2)
     tmpArray5(32,2) = Zpa*tmpArray4(17,2) + alphaZpq*TmpArray4(17,3)
     tmpArray5(33,2) = Ypa*tmpArray4(19,2) + alphaYpq*TmpArray4(19,3) + TwoTerms(3)
     tmpArray5(34,2) = Ypa*tmpArray4(20,2) + alphaYpq*TmpArray4(20,3)
     tmpArray5(35,2) = Zpa*tmpArray4(20,2) + alphaZpq*TmpArray4(20,3) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(11) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(14) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(16) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(17) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(19) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(20) + alphaP*TmpArray4(20,2))
     do iTUV = 1,   35
      AuxArray(iTUV,IP) = AuxArray(iTUV,IP) + TMPAuxarray(iTUV)
     enddo
     AuxArray(36,IP) = AuxArray(36,IP) + Xpa*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IP) = AuxArray(37,IP) + Ypa*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IP) = AuxArray(38,IP) + Zpa*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IP) = AuxArray(39,IP) + Xpa*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IP) = AuxArray(40,IP) + Ypa*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IP) = AuxArray(41,IP) + Xpa*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IP) = AuxArray(42,IP) + Xpa*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IP) = AuxArray(43,IP) + Zpa*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IP) = AuxArray(44,IP) + Ypa*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IP) = AuxArray(45,IP) + Xpa*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IP) = AuxArray(46,IP) + Xpa*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IP) = AuxArray(47,IP) + Xpa*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IP) = AuxArray(48,IP) + Xpa*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IP) = AuxArray(49,IP) + Xpa*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IP) = AuxArray(50,IP) + Xpa*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IP) = AuxArray(51,IP) + Ypa*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IP) = AuxArray(52,IP) + Zpa*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IP) = AuxArray(53,IP) + Ypa*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IP) = AuxArray(54,IP) + Ypa*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IP) = AuxArray(55,IP) + Ypa*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IP) = AuxArray(56,IP) + Zpa*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSegQ6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0: 9,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Acenter(3)
  real(realk),intent(inout) :: AUXarray(   84,nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   56)
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
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
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Ax = -Acenter(1)
  Ay = -Acenter(2)
  Az = -Acenter(3)
  DO iPassQ = 1,nPasses
   iP = (iPassQ-1)*nPrimP
   DO iPrimP=1, nPrimP
    iP = iP + 1
    DO iTUV=1,   84
     AUXarray(iTUV,iP)=0.0E0_realk
    ENDDO
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP) + Ax
    Ypa = Pcent(2,iPrimP) + Ay
    Zpa = Pcent(3,iPrimP) + Az
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
     tmpArray2(2,5) = Xpa*tmpArray1(1,5) + alphaXpq*TmpArray1(1,6)
     tmpArray2(3,5) = Ypa*tmpArray1(1,5) + alphaYpq*TmpArray1(1,6)
     tmpArray2(4,5) = Zpa*tmpArray1(1,5) + alphaZpq*TmpArray1(1,6)
     tmpArray2(2,6) = Xpa*tmpArray1(1,6) + alphaXpq*TmpArray1(1,7)
     tmpArray2(3,6) = Ypa*tmpArray1(1,6) + alphaYpq*TmpArray1(1,7)
     tmpArray2(4,6) = Zpa*tmpArray1(1,6) + alphaZpq*TmpArray1(1,7)
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
     TwoTerms(1) = inv2expP*(TmpArray1(1,4) + alphaP*TmpArray1(1,5))
     tmpArray3(5,4) = Xpa*tmpArray2(2,4) + alphaXpq*TmpArray2(2,5) + TwoTerms(1)
     tmpArray3(6,4) = Xpa*tmpArray2(3,4) + alphaXpq*TmpArray2(3,5)
     tmpArray3(7,4) = Xpa*tmpArray2(4,4) + alphaXpq*TmpArray2(4,5)
     tmpArray3(8,4) = Ypa*tmpArray2(3,4) + alphaYpq*TmpArray2(3,5) + TwoTerms(1)
     tmpArray3(9,4) = Ypa*tmpArray2(4,4) + alphaYpq*TmpArray2(4,5)
     tmpArray3(10,4) = Zpa*tmpArray2(4,4) + alphaZpq*TmpArray2(4,5) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,5) + alphaP*TmpArray1(1,6))
     tmpArray3(5,5) = Xpa*tmpArray2(2,5) + alphaXpq*TmpArray2(2,6) + TwoTerms(1)
     tmpArray3(6,5) = Xpa*tmpArray2(3,5) + alphaXpq*TmpArray2(3,6)
     tmpArray3(7,5) = Xpa*tmpArray2(4,5) + alphaXpq*TmpArray2(4,6)
     tmpArray3(8,5) = Ypa*tmpArray2(3,5) + alphaYpq*TmpArray2(3,6) + TwoTerms(1)
     tmpArray3(9,5) = Ypa*tmpArray2(4,5) + alphaYpq*TmpArray2(4,6)
     tmpArray3(10,5) = Zpa*tmpArray2(4,5) + alphaZpq*TmpArray2(4,6) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(TmpArray2(2,3) + alphaP*TmpArray2(2,4))
     TwoTerms(2) = inv2expP*(TmpArray2(3,3) + alphaP*TmpArray2(3,4))
     TwoTerms(3) = inv2expP*(TmpArray2(4,3) + alphaP*TmpArray2(4,4))
     tmpArray4(11,3) = Xpa*tmpArray3(5,3) + alphaXpq*TmpArray3(5,4) + 2*TwoTerms(1)
     tmpArray4(12,3) = Ypa*tmpArray3(5,3) + alphaYpq*TmpArray3(5,4)
     tmpArray4(13,3) = Zpa*tmpArray3(5,3) + alphaZpq*TmpArray3(5,4)
     tmpArray4(14,3) = Xpa*tmpArray3(8,3) + alphaXpq*TmpArray3(8,4)
     tmpArray4(15,3) = Xpa*tmpArray3(9,3) + alphaXpq*TmpArray3(9,4)
     tmpArray4(16,3) = Xpa*tmpArray3(10,3) + alphaXpq*TmpArray3(10,4)
     tmpArray4(17,3) = Ypa*tmpArray3(8,3) + alphaYpq*TmpArray3(8,4) + 2*TwoTerms(2)
     tmpArray4(18,3) = Zpa*tmpArray3(8,3) + alphaZpq*TmpArray3(8,4)
     tmpArray4(19,3) = Ypa*tmpArray3(10,3) + alphaYpq*TmpArray3(10,4)
     tmpArray4(20,3) = Zpa*tmpArray3(10,3) + alphaZpq*TmpArray3(10,4) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,4) + alphaP*TmpArray2(2,5))
     TwoTerms(2) = inv2expP*(TmpArray2(3,4) + alphaP*TmpArray2(3,5))
     TwoTerms(3) = inv2expP*(TmpArray2(4,4) + alphaP*TmpArray2(4,5))
     tmpArray4(11,4) = Xpa*tmpArray3(5,4) + alphaXpq*TmpArray3(5,5) + 2*TwoTerms(1)
     tmpArray4(12,4) = Ypa*tmpArray3(5,4) + alphaYpq*TmpArray3(5,5)
     tmpArray4(13,4) = Zpa*tmpArray3(5,4) + alphaZpq*TmpArray3(5,5)
     tmpArray4(14,4) = Xpa*tmpArray3(8,4) + alphaXpq*TmpArray3(8,5)
     tmpArray4(15,4) = Xpa*tmpArray3(9,4) + alphaXpq*TmpArray3(9,5)
     tmpArray4(16,4) = Xpa*tmpArray3(10,4) + alphaXpq*TmpArray3(10,5)
     tmpArray4(17,4) = Ypa*tmpArray3(8,4) + alphaYpq*TmpArray3(8,5) + 2*TwoTerms(2)
     tmpArray4(18,4) = Zpa*tmpArray3(8,4) + alphaZpq*TmpArray3(8,5)
     tmpArray4(19,4) = Ypa*tmpArray3(10,4) + alphaYpq*TmpArray3(10,5)
     tmpArray4(20,4) = Zpa*tmpArray3(10,4) + alphaZpq*TmpArray3(10,5) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(10) + alphaP*TmpArray3(10,2))
     TMPAuxArray(21) = Xpa*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Ypa*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zpa*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xpa*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Ypa*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xpa*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xpa*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xpa*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xpa*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xpa*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Ypa*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zpa*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Ypa*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Ypa*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zpa*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,2) + alphaP*TmpArray3(5,3))
     TwoTerms(2) = inv2expP*(TmpArray3(8,2) + alphaP*TmpArray3(8,3))
     TwoTerms(3) = inv2expP*(TmpArray3(10,2) + alphaP*TmpArray3(10,3))
     tmpArray5(21,2) = Xpa*tmpArray4(11,2) + alphaXpq*TmpArray4(11,3) + 3*TwoTerms(1)
     tmpArray5(22,2) = Ypa*tmpArray4(11,2) + alphaYpq*TmpArray4(11,3)
     tmpArray5(23,2) = Zpa*tmpArray4(11,2) + alphaZpq*TmpArray4(11,3)
     tmpArray5(24,2) = Xpa*tmpArray4(14,2) + alphaXpq*TmpArray4(14,3) + TwoTerms(2)
     tmpArray5(25,2) = Ypa*tmpArray4(13,2) + alphaYpq*TmpArray4(13,3)
     tmpArray5(26,2) = Xpa*tmpArray4(16,2) + alphaXpq*TmpArray4(16,3) + TwoTerms(3)
     tmpArray5(27,2) = Xpa*tmpArray4(17,2) + alphaXpq*TmpArray4(17,3)
     tmpArray5(28,2) = Xpa*tmpArray4(18,2) + alphaXpq*TmpArray4(18,3)
     tmpArray5(29,2) = Xpa*tmpArray4(19,2) + alphaXpq*TmpArray4(19,3)
     tmpArray5(30,2) = Xpa*tmpArray4(20,2) + alphaXpq*TmpArray4(20,3)
     tmpArray5(31,2) = Ypa*tmpArray4(17,2) + alphaYpq*TmpArray4(17,3) + 3*TwoTerms(2)
     tmpArray5(32,2) = Zpa*tmpArray4(17,2) + alphaZpq*TmpArray4(17,3)
     tmpArray5(33,2) = Ypa*tmpArray4(19,2) + alphaYpq*TmpArray4(19,3) + TwoTerms(3)
     tmpArray5(34,2) = Ypa*tmpArray4(20,2) + alphaYpq*TmpArray4(20,3)
     tmpArray5(35,2) = Zpa*tmpArray4(20,2) + alphaZpq*TmpArray4(20,3) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,3) + alphaP*TmpArray3(5,4))
     TwoTerms(2) = inv2expP*(TmpArray3(8,3) + alphaP*TmpArray3(8,4))
     TwoTerms(3) = inv2expP*(TmpArray3(10,3) + alphaP*TmpArray3(10,4))
     tmpArray5(21,3) = Xpa*tmpArray4(11,3) + alphaXpq*TmpArray4(11,4) + 3*TwoTerms(1)
     tmpArray5(22,3) = Ypa*tmpArray4(11,3) + alphaYpq*TmpArray4(11,4)
     tmpArray5(23,3) = Zpa*tmpArray4(11,3) + alphaZpq*TmpArray4(11,4)
     tmpArray5(24,3) = Xpa*tmpArray4(14,3) + alphaXpq*TmpArray4(14,4) + TwoTerms(2)
     tmpArray5(25,3) = Ypa*tmpArray4(13,3) + alphaYpq*TmpArray4(13,4)
     tmpArray5(26,3) = Xpa*tmpArray4(16,3) + alphaXpq*TmpArray4(16,4) + TwoTerms(3)
     tmpArray5(27,3) = Xpa*tmpArray4(17,3) + alphaXpq*TmpArray4(17,4)
     tmpArray5(28,3) = Xpa*tmpArray4(18,3) + alphaXpq*TmpArray4(18,4)
     tmpArray5(29,3) = Xpa*tmpArray4(19,3) + alphaXpq*TmpArray4(19,4)
     tmpArray5(30,3) = Xpa*tmpArray4(20,3) + alphaXpq*TmpArray4(20,4)
     tmpArray5(31,3) = Ypa*tmpArray4(17,3) + alphaYpq*TmpArray4(17,4) + 3*TwoTerms(2)
     tmpArray5(32,3) = Zpa*tmpArray4(17,3) + alphaZpq*TmpArray4(17,4)
     tmpArray5(33,3) = Ypa*tmpArray4(19,3) + alphaYpq*TmpArray4(19,4) + TwoTerms(3)
     tmpArray5(34,3) = Ypa*tmpArray4(20,3) + alphaYpq*TmpArray4(20,4)
     tmpArray5(35,3) = Zpa*tmpArray4(20,3) + alphaZpq*TmpArray4(20,4) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(11) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(14) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(16) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(17) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(19) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(20) + alphaP*TmpArray4(20,2))
     TMPAuxArray(36) = Xpa*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     TMPAuxArray(37) = Ypa*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     TMPAuxArray(38) = Zpa*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     TMPAuxArray(39) = Xpa*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     TMPAuxArray(40) = Ypa*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     TMPAuxArray(41) = Xpa*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     TMPAuxArray(42) = Xpa*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     TMPAuxArray(43) = Zpa*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     TMPAuxArray(44) = Ypa*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     TMPAuxArray(45) = Xpa*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     TMPAuxArray(46) = Xpa*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     TMPAuxArray(47) = Xpa*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     TMPAuxArray(48) = Xpa*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     TMPAuxArray(49) = Xpa*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     TMPAuxArray(50) = Xpa*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     TMPAuxArray(51) = Ypa*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     TMPAuxArray(52) = Zpa*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     TMPAuxArray(53) = Ypa*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     TMPAuxArray(54) = Ypa*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     TMPAuxArray(55) = Ypa*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     TMPAuxArray(56) = Zpa*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TmpArray4(11,2) + alphaP*TmpArray4(11,3))
     TwoTerms(2) = inv2expP*(TmpArray4(14,2) + alphaP*TmpArray4(14,3))
     TwoTerms(3) = inv2expP*(TmpArray4(16,2) + alphaP*TmpArray4(16,3))
     TwoTerms(4) = inv2expP*(TmpArray4(17,2) + alphaP*TmpArray4(17,3))
     TwoTerms(5) = inv2expP*(TmpArray4(19,2) + alphaP*TmpArray4(19,3))
     TwoTerms(6) = inv2expP*(TmpArray4(20,2) + alphaP*TmpArray4(20,3))
     tmpArray6(36,2) = Xpa*tmpArray5(21,2) + alphaXpq*TmpArray5(21,3) + 4*TwoTerms(1)
     tmpArray6(37,2) = Ypa*tmpArray5(21,2) + alphaYpq*TmpArray5(21,3)
     tmpArray6(38,2) = Zpa*tmpArray5(21,2) + alphaZpq*TmpArray5(21,3)
     tmpArray6(39,2) = Xpa*tmpArray5(24,2) + alphaXpq*TmpArray5(24,3) + 2*TwoTerms(2)
     tmpArray6(40,2) = Ypa*tmpArray5(23,2) + alphaYpq*TmpArray5(23,3)
     tmpArray6(41,2) = Xpa*tmpArray5(26,2) + alphaXpq*TmpArray5(26,3) + 2*TwoTerms(3)
     tmpArray6(42,2) = Xpa*tmpArray5(27,2) + alphaXpq*TmpArray5(27,3) + TwoTerms(4)
     tmpArray6(43,2) = Zpa*tmpArray5(24,2) + alphaZpq*TmpArray5(24,3)
     tmpArray6(44,2) = Ypa*tmpArray5(26,2) + alphaYpq*TmpArray5(26,3)
     tmpArray6(45,2) = Xpa*tmpArray5(30,2) + alphaXpq*TmpArray5(30,3) + TwoTerms(6)
     tmpArray6(46,2) = Xpa*tmpArray5(31,2) + alphaXpq*TmpArray5(31,3)
     tmpArray6(47,2) = Xpa*tmpArray5(32,2) + alphaXpq*TmpArray5(32,3)
     tmpArray6(48,2) = Xpa*tmpArray5(33,2) + alphaXpq*TmpArray5(33,3)
     tmpArray6(49,2) = Xpa*tmpArray5(34,2) + alphaXpq*TmpArray5(34,3)
     tmpArray6(50,2) = Xpa*tmpArray5(35,2) + alphaXpq*TmpArray5(35,3)
     tmpArray6(51,2) = Ypa*tmpArray5(31,2) + alphaYpq*TmpArray5(31,3) + 4*TwoTerms(4)
     tmpArray6(52,2) = Zpa*tmpArray5(31,2) + alphaZpq*TmpArray5(31,3)
     tmpArray6(53,2) = Ypa*tmpArray5(33,2) + alphaYpq*TmpArray5(33,3) + 2*TwoTerms(5)
     tmpArray6(54,2) = Ypa*tmpArray5(34,2) + alphaYpq*TmpArray5(34,3) + TwoTerms(6)
     tmpArray6(55,2) = Ypa*tmpArray5(35,2) + alphaYpq*TmpArray5(35,3)
     tmpArray6(56,2) = Zpa*tmpArray5(35,2) + alphaZpq*TmpArray5(35,3) + 4*TwoTerms(6)
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
      AuxArray(iTUV,IP) = AuxArray(iTUV,IP) + TMPAuxarray(iTUV)
     enddo
     AuxArray(57,IP) = AuxArray(57,IP) + Xpa*TMPAuxArray(36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,IP) = AuxArray(58,IP) + Ypa*TMPAuxArray(36) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,IP) = AuxArray(59,IP) + Zpa*TMPAuxArray(36) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,IP) = AuxArray(60,IP) + Xpa*TMPAuxArray(39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,IP) = AuxArray(61,IP) + Ypa*TMPAuxArray(38) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,IP) = AuxArray(62,IP) + Xpa*TMPAuxArray(41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,IP) = AuxArray(63,IP) + Xpa*TMPAuxArray(42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,IP) = AuxArray(64,IP) + Zpa*TMPAuxArray(39) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,IP) = AuxArray(65,IP) + Ypa*TMPAuxArray(41) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,IP) = AuxArray(66,IP) + Xpa*TMPAuxArray(45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,IP) = AuxArray(67,IP) + Xpa*TMPAuxArray(46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,IP) = AuxArray(68,IP) + Zpa*TMPAuxArray(42) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,IP) = AuxArray(69,IP) + Xpa*TMPAuxArray(48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,IP) = AuxArray(70,IP) + Ypa*TMPAuxArray(45) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,IP) = AuxArray(71,IP) + Xpa*TMPAuxArray(50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,IP) = AuxArray(72,IP) + Xpa*TMPAuxArray(51) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,IP) = AuxArray(73,IP) + Xpa*TMPAuxArray(52) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,IP) = AuxArray(74,IP) + Xpa*TMPAuxArray(53) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,IP) = AuxArray(75,IP) + Xpa*TMPAuxArray(54) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,IP) = AuxArray(76,IP) + Xpa*TMPAuxArray(55) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,IP) = AuxArray(77,IP) + Xpa*TMPAuxArray(56) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,IP) = AuxArray(78,IP) + Ypa*TMPAuxArray(51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,IP) = AuxArray(79,IP) + Zpa*TMPAuxArray(51) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,IP) = AuxArray(80,IP) + Ypa*TMPAuxArray(53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,IP) = AuxArray(81,IP) + Ypa*TMPAuxArray(54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,IP) = AuxArray(82,IP) + Ypa*TMPAuxArray(55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,IP) = AuxArray(83,IP) + Ypa*TMPAuxArray(56) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,IP) = AuxArray(84,IP) + Zpa*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSegQ7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:10,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Acenter(3)
  real(realk),intent(inout) :: AUXarray(  120,nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(   84)
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
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
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Ax = -Acenter(1)
  Ay = -Acenter(2)
  Az = -Acenter(3)
  DO iPassQ = 1,nPasses
   iP = (iPassQ-1)*nPrimP
   DO iPrimP=1, nPrimP
    iP = iP + 1
    DO iTUV=1,  120
     AUXarray(iTUV,iP)=0.0E0_realk
    ENDDO
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP) + Ax
    Ypa = Pcent(2,iPrimP) + Ay
    Zpa = Pcent(3,iPrimP) + Az
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
     tmpArray2(2,5) = Xpa*tmpArray1(1,5) + alphaXpq*TmpArray1(1,6)
     tmpArray2(3,5) = Ypa*tmpArray1(1,5) + alphaYpq*TmpArray1(1,6)
     tmpArray2(4,5) = Zpa*tmpArray1(1,5) + alphaZpq*TmpArray1(1,6)
     tmpArray2(2,6) = Xpa*tmpArray1(1,6) + alphaXpq*TmpArray1(1,7)
     tmpArray2(3,6) = Ypa*tmpArray1(1,6) + alphaYpq*TmpArray1(1,7)
     tmpArray2(4,6) = Zpa*tmpArray1(1,6) + alphaZpq*TmpArray1(1,7)
     tmpArray2(2,7) = Xpa*tmpArray1(1,7) + alphaXpq*TmpArray1(1,8)
     tmpArray2(3,7) = Ypa*tmpArray1(1,7) + alphaYpq*TmpArray1(1,8)
     tmpArray2(4,7) = Zpa*tmpArray1(1,7) + alphaZpq*TmpArray1(1,8)
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
     TwoTerms(1) = inv2expP*(TmpArray1(1,4) + alphaP*TmpArray1(1,5))
     tmpArray3(5,4) = Xpa*tmpArray2(2,4) + alphaXpq*TmpArray2(2,5) + TwoTerms(1)
     tmpArray3(6,4) = Xpa*tmpArray2(3,4) + alphaXpq*TmpArray2(3,5)
     tmpArray3(7,4) = Xpa*tmpArray2(4,4) + alphaXpq*TmpArray2(4,5)
     tmpArray3(8,4) = Ypa*tmpArray2(3,4) + alphaYpq*TmpArray2(3,5) + TwoTerms(1)
     tmpArray3(9,4) = Ypa*tmpArray2(4,4) + alphaYpq*TmpArray2(4,5)
     tmpArray3(10,4) = Zpa*tmpArray2(4,4) + alphaZpq*TmpArray2(4,5) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,5) + alphaP*TmpArray1(1,6))
     tmpArray3(5,5) = Xpa*tmpArray2(2,5) + alphaXpq*TmpArray2(2,6) + TwoTerms(1)
     tmpArray3(6,5) = Xpa*tmpArray2(3,5) + alphaXpq*TmpArray2(3,6)
     tmpArray3(7,5) = Xpa*tmpArray2(4,5) + alphaXpq*TmpArray2(4,6)
     tmpArray3(8,5) = Ypa*tmpArray2(3,5) + alphaYpq*TmpArray2(3,6) + TwoTerms(1)
     tmpArray3(9,5) = Ypa*tmpArray2(4,5) + alphaYpq*TmpArray2(4,6)
     tmpArray3(10,5) = Zpa*tmpArray2(4,5) + alphaZpq*TmpArray2(4,6) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,6) + alphaP*TmpArray1(1,7))
     tmpArray3(5,6) = Xpa*tmpArray2(2,6) + alphaXpq*TmpArray2(2,7) + TwoTerms(1)
     tmpArray3(6,6) = Xpa*tmpArray2(3,6) + alphaXpq*TmpArray2(3,7)
     tmpArray3(7,6) = Xpa*tmpArray2(4,6) + alphaXpq*TmpArray2(4,7)
     tmpArray3(8,6) = Ypa*tmpArray2(3,6) + alphaYpq*TmpArray2(3,7) + TwoTerms(1)
     tmpArray3(9,6) = Ypa*tmpArray2(4,6) + alphaYpq*TmpArray2(4,7)
     tmpArray3(10,6) = Zpa*tmpArray2(4,6) + alphaZpq*TmpArray2(4,7) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(TmpArray2(2,3) + alphaP*TmpArray2(2,4))
     TwoTerms(2) = inv2expP*(TmpArray2(3,3) + alphaP*TmpArray2(3,4))
     TwoTerms(3) = inv2expP*(TmpArray2(4,3) + alphaP*TmpArray2(4,4))
     tmpArray4(11,3) = Xpa*tmpArray3(5,3) + alphaXpq*TmpArray3(5,4) + 2*TwoTerms(1)
     tmpArray4(12,3) = Ypa*tmpArray3(5,3) + alphaYpq*TmpArray3(5,4)
     tmpArray4(13,3) = Zpa*tmpArray3(5,3) + alphaZpq*TmpArray3(5,4)
     tmpArray4(14,3) = Xpa*tmpArray3(8,3) + alphaXpq*TmpArray3(8,4)
     tmpArray4(15,3) = Xpa*tmpArray3(9,3) + alphaXpq*TmpArray3(9,4)
     tmpArray4(16,3) = Xpa*tmpArray3(10,3) + alphaXpq*TmpArray3(10,4)
     tmpArray4(17,3) = Ypa*tmpArray3(8,3) + alphaYpq*TmpArray3(8,4) + 2*TwoTerms(2)
     tmpArray4(18,3) = Zpa*tmpArray3(8,3) + alphaZpq*TmpArray3(8,4)
     tmpArray4(19,3) = Ypa*tmpArray3(10,3) + alphaYpq*TmpArray3(10,4)
     tmpArray4(20,3) = Zpa*tmpArray3(10,3) + alphaZpq*TmpArray3(10,4) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,4) + alphaP*TmpArray2(2,5))
     TwoTerms(2) = inv2expP*(TmpArray2(3,4) + alphaP*TmpArray2(3,5))
     TwoTerms(3) = inv2expP*(TmpArray2(4,4) + alphaP*TmpArray2(4,5))
     tmpArray4(11,4) = Xpa*tmpArray3(5,4) + alphaXpq*TmpArray3(5,5) + 2*TwoTerms(1)
     tmpArray4(12,4) = Ypa*tmpArray3(5,4) + alphaYpq*TmpArray3(5,5)
     tmpArray4(13,4) = Zpa*tmpArray3(5,4) + alphaZpq*TmpArray3(5,5)
     tmpArray4(14,4) = Xpa*tmpArray3(8,4) + alphaXpq*TmpArray3(8,5)
     tmpArray4(15,4) = Xpa*tmpArray3(9,4) + alphaXpq*TmpArray3(9,5)
     tmpArray4(16,4) = Xpa*tmpArray3(10,4) + alphaXpq*TmpArray3(10,5)
     tmpArray4(17,4) = Ypa*tmpArray3(8,4) + alphaYpq*TmpArray3(8,5) + 2*TwoTerms(2)
     tmpArray4(18,4) = Zpa*tmpArray3(8,4) + alphaZpq*TmpArray3(8,5)
     tmpArray4(19,4) = Ypa*tmpArray3(10,4) + alphaYpq*TmpArray3(10,5)
     tmpArray4(20,4) = Zpa*tmpArray3(10,4) + alphaZpq*TmpArray3(10,5) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,5) + alphaP*TmpArray2(2,6))
     TwoTerms(2) = inv2expP*(TmpArray2(3,5) + alphaP*TmpArray2(3,6))
     TwoTerms(3) = inv2expP*(TmpArray2(4,5) + alphaP*TmpArray2(4,6))
     tmpArray4(11,5) = Xpa*tmpArray3(5,5) + alphaXpq*TmpArray3(5,6) + 2*TwoTerms(1)
     tmpArray4(12,5) = Ypa*tmpArray3(5,5) + alphaYpq*TmpArray3(5,6)
     tmpArray4(13,5) = Zpa*tmpArray3(5,5) + alphaZpq*TmpArray3(5,6)
     tmpArray4(14,5) = Xpa*tmpArray3(8,5) + alphaXpq*TmpArray3(8,6)
     tmpArray4(15,5) = Xpa*tmpArray3(9,5) + alphaXpq*TmpArray3(9,6)
     tmpArray4(16,5) = Xpa*tmpArray3(10,5) + alphaXpq*TmpArray3(10,6)
     tmpArray4(17,5) = Ypa*tmpArray3(8,5) + alphaYpq*TmpArray3(8,6) + 2*TwoTerms(2)
     tmpArray4(18,5) = Zpa*tmpArray3(8,5) + alphaZpq*TmpArray3(8,6)
     tmpArray4(19,5) = Ypa*tmpArray3(10,5) + alphaYpq*TmpArray3(10,6)
     tmpArray4(20,5) = Zpa*tmpArray3(10,5) + alphaZpq*TmpArray3(10,6) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(10) + alphaP*TmpArray3(10,2))
     TMPAuxArray(21) = Xpa*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Ypa*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zpa*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xpa*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Ypa*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xpa*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xpa*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xpa*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xpa*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xpa*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Ypa*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zpa*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Ypa*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Ypa*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zpa*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,2) + alphaP*TmpArray3(5,3))
     TwoTerms(2) = inv2expP*(TmpArray3(8,2) + alphaP*TmpArray3(8,3))
     TwoTerms(3) = inv2expP*(TmpArray3(10,2) + alphaP*TmpArray3(10,3))
     tmpArray5(21,2) = Xpa*tmpArray4(11,2) + alphaXpq*TmpArray4(11,3) + 3*TwoTerms(1)
     tmpArray5(22,2) = Ypa*tmpArray4(11,2) + alphaYpq*TmpArray4(11,3)
     tmpArray5(23,2) = Zpa*tmpArray4(11,2) + alphaZpq*TmpArray4(11,3)
     tmpArray5(24,2) = Xpa*tmpArray4(14,2) + alphaXpq*TmpArray4(14,3) + TwoTerms(2)
     tmpArray5(25,2) = Ypa*tmpArray4(13,2) + alphaYpq*TmpArray4(13,3)
     tmpArray5(26,2) = Xpa*tmpArray4(16,2) + alphaXpq*TmpArray4(16,3) + TwoTerms(3)
     tmpArray5(27,2) = Xpa*tmpArray4(17,2) + alphaXpq*TmpArray4(17,3)
     tmpArray5(28,2) = Xpa*tmpArray4(18,2) + alphaXpq*TmpArray4(18,3)
     tmpArray5(29,2) = Xpa*tmpArray4(19,2) + alphaXpq*TmpArray4(19,3)
     tmpArray5(30,2) = Xpa*tmpArray4(20,2) + alphaXpq*TmpArray4(20,3)
     tmpArray5(31,2) = Ypa*tmpArray4(17,2) + alphaYpq*TmpArray4(17,3) + 3*TwoTerms(2)
     tmpArray5(32,2) = Zpa*tmpArray4(17,2) + alphaZpq*TmpArray4(17,3)
     tmpArray5(33,2) = Ypa*tmpArray4(19,2) + alphaYpq*TmpArray4(19,3) + TwoTerms(3)
     tmpArray5(34,2) = Ypa*tmpArray4(20,2) + alphaYpq*TmpArray4(20,3)
     tmpArray5(35,2) = Zpa*tmpArray4(20,2) + alphaZpq*TmpArray4(20,3) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,3) + alphaP*TmpArray3(5,4))
     TwoTerms(2) = inv2expP*(TmpArray3(8,3) + alphaP*TmpArray3(8,4))
     TwoTerms(3) = inv2expP*(TmpArray3(10,3) + alphaP*TmpArray3(10,4))
     tmpArray5(21,3) = Xpa*tmpArray4(11,3) + alphaXpq*TmpArray4(11,4) + 3*TwoTerms(1)
     tmpArray5(22,3) = Ypa*tmpArray4(11,3) + alphaYpq*TmpArray4(11,4)
     tmpArray5(23,3) = Zpa*tmpArray4(11,3) + alphaZpq*TmpArray4(11,4)
     tmpArray5(24,3) = Xpa*tmpArray4(14,3) + alphaXpq*TmpArray4(14,4) + TwoTerms(2)
     tmpArray5(25,3) = Ypa*tmpArray4(13,3) + alphaYpq*TmpArray4(13,4)
     tmpArray5(26,3) = Xpa*tmpArray4(16,3) + alphaXpq*TmpArray4(16,4) + TwoTerms(3)
     tmpArray5(27,3) = Xpa*tmpArray4(17,3) + alphaXpq*TmpArray4(17,4)
     tmpArray5(28,3) = Xpa*tmpArray4(18,3) + alphaXpq*TmpArray4(18,4)
     tmpArray5(29,3) = Xpa*tmpArray4(19,3) + alphaXpq*TmpArray4(19,4)
     tmpArray5(30,3) = Xpa*tmpArray4(20,3) + alphaXpq*TmpArray4(20,4)
     tmpArray5(31,3) = Ypa*tmpArray4(17,3) + alphaYpq*TmpArray4(17,4) + 3*TwoTerms(2)
     tmpArray5(32,3) = Zpa*tmpArray4(17,3) + alphaZpq*TmpArray4(17,4)
     tmpArray5(33,3) = Ypa*tmpArray4(19,3) + alphaYpq*TmpArray4(19,4) + TwoTerms(3)
     tmpArray5(34,3) = Ypa*tmpArray4(20,3) + alphaYpq*TmpArray4(20,4)
     tmpArray5(35,3) = Zpa*tmpArray4(20,3) + alphaZpq*TmpArray4(20,4) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,4) + alphaP*TmpArray3(5,5))
     TwoTerms(2) = inv2expP*(TmpArray3(8,4) + alphaP*TmpArray3(8,5))
     TwoTerms(3) = inv2expP*(TmpArray3(10,4) + alphaP*TmpArray3(10,5))
     tmpArray5(21,4) = Xpa*tmpArray4(11,4) + alphaXpq*TmpArray4(11,5) + 3*TwoTerms(1)
     tmpArray5(22,4) = Ypa*tmpArray4(11,4) + alphaYpq*TmpArray4(11,5)
     tmpArray5(23,4) = Zpa*tmpArray4(11,4) + alphaZpq*TmpArray4(11,5)
     tmpArray5(24,4) = Xpa*tmpArray4(14,4) + alphaXpq*TmpArray4(14,5) + TwoTerms(2)
     tmpArray5(25,4) = Ypa*tmpArray4(13,4) + alphaYpq*TmpArray4(13,5)
     tmpArray5(26,4) = Xpa*tmpArray4(16,4) + alphaXpq*TmpArray4(16,5) + TwoTerms(3)
     tmpArray5(27,4) = Xpa*tmpArray4(17,4) + alphaXpq*TmpArray4(17,5)
     tmpArray5(28,4) = Xpa*tmpArray4(18,4) + alphaXpq*TmpArray4(18,5)
     tmpArray5(29,4) = Xpa*tmpArray4(19,4) + alphaXpq*TmpArray4(19,5)
     tmpArray5(30,4) = Xpa*tmpArray4(20,4) + alphaXpq*TmpArray4(20,5)
     tmpArray5(31,4) = Ypa*tmpArray4(17,4) + alphaYpq*TmpArray4(17,5) + 3*TwoTerms(2)
     tmpArray5(32,4) = Zpa*tmpArray4(17,4) + alphaZpq*TmpArray4(17,5)
     tmpArray5(33,4) = Ypa*tmpArray4(19,4) + alphaYpq*TmpArray4(19,5) + TwoTerms(3)
     tmpArray5(34,4) = Ypa*tmpArray4(20,4) + alphaYpq*TmpArray4(20,5)
     tmpArray5(35,4) = Zpa*tmpArray4(20,4) + alphaZpq*TmpArray4(20,5) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(11) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(14) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(16) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(17) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(19) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(20) + alphaP*TmpArray4(20,2))
     TMPAuxArray(36) = Xpa*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     TMPAuxArray(37) = Ypa*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     TMPAuxArray(38) = Zpa*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     TMPAuxArray(39) = Xpa*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     TMPAuxArray(40) = Ypa*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     TMPAuxArray(41) = Xpa*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     TMPAuxArray(42) = Xpa*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     TMPAuxArray(43) = Zpa*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     TMPAuxArray(44) = Ypa*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     TMPAuxArray(45) = Xpa*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     TMPAuxArray(46) = Xpa*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     TMPAuxArray(47) = Xpa*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     TMPAuxArray(48) = Xpa*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     TMPAuxArray(49) = Xpa*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     TMPAuxArray(50) = Xpa*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     TMPAuxArray(51) = Ypa*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     TMPAuxArray(52) = Zpa*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     TMPAuxArray(53) = Ypa*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     TMPAuxArray(54) = Ypa*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     TMPAuxArray(55) = Ypa*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     TMPAuxArray(56) = Zpa*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TmpArray4(11,2) + alphaP*TmpArray4(11,3))
     TwoTerms(2) = inv2expP*(TmpArray4(14,2) + alphaP*TmpArray4(14,3))
     TwoTerms(3) = inv2expP*(TmpArray4(16,2) + alphaP*TmpArray4(16,3))
     TwoTerms(4) = inv2expP*(TmpArray4(17,2) + alphaP*TmpArray4(17,3))
     TwoTerms(5) = inv2expP*(TmpArray4(19,2) + alphaP*TmpArray4(19,3))
     TwoTerms(6) = inv2expP*(TmpArray4(20,2) + alphaP*TmpArray4(20,3))
     tmpArray6(36,2) = Xpa*tmpArray5(21,2) + alphaXpq*TmpArray5(21,3) + 4*TwoTerms(1)
     tmpArray6(37,2) = Ypa*tmpArray5(21,2) + alphaYpq*TmpArray5(21,3)
     tmpArray6(38,2) = Zpa*tmpArray5(21,2) + alphaZpq*TmpArray5(21,3)
     tmpArray6(39,2) = Xpa*tmpArray5(24,2) + alphaXpq*TmpArray5(24,3) + 2*TwoTerms(2)
     tmpArray6(40,2) = Ypa*tmpArray5(23,2) + alphaYpq*TmpArray5(23,3)
     tmpArray6(41,2) = Xpa*tmpArray5(26,2) + alphaXpq*TmpArray5(26,3) + 2*TwoTerms(3)
     tmpArray6(42,2) = Xpa*tmpArray5(27,2) + alphaXpq*TmpArray5(27,3) + TwoTerms(4)
     tmpArray6(43,2) = Zpa*tmpArray5(24,2) + alphaZpq*TmpArray5(24,3)
     tmpArray6(44,2) = Ypa*tmpArray5(26,2) + alphaYpq*TmpArray5(26,3)
     tmpArray6(45,2) = Xpa*tmpArray5(30,2) + alphaXpq*TmpArray5(30,3) + TwoTerms(6)
     tmpArray6(46,2) = Xpa*tmpArray5(31,2) + alphaXpq*TmpArray5(31,3)
     tmpArray6(47,2) = Xpa*tmpArray5(32,2) + alphaXpq*TmpArray5(32,3)
     tmpArray6(48,2) = Xpa*tmpArray5(33,2) + alphaXpq*TmpArray5(33,3)
     tmpArray6(49,2) = Xpa*tmpArray5(34,2) + alphaXpq*TmpArray5(34,3)
     tmpArray6(50,2) = Xpa*tmpArray5(35,2) + alphaXpq*TmpArray5(35,3)
     tmpArray6(51,2) = Ypa*tmpArray5(31,2) + alphaYpq*TmpArray5(31,3) + 4*TwoTerms(4)
     tmpArray6(52,2) = Zpa*tmpArray5(31,2) + alphaZpq*TmpArray5(31,3)
     tmpArray6(53,2) = Ypa*tmpArray5(33,2) + alphaYpq*TmpArray5(33,3) + 2*TwoTerms(5)
     tmpArray6(54,2) = Ypa*tmpArray5(34,2) + alphaYpq*TmpArray5(34,3) + TwoTerms(6)
     tmpArray6(55,2) = Ypa*tmpArray5(35,2) + alphaYpq*TmpArray5(35,3)
     tmpArray6(56,2) = Zpa*tmpArray5(35,2) + alphaZpq*TmpArray5(35,3) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TmpArray4(11,3) + alphaP*TmpArray4(11,4))
     TwoTerms(2) = inv2expP*(TmpArray4(14,3) + alphaP*TmpArray4(14,4))
     TwoTerms(3) = inv2expP*(TmpArray4(16,3) + alphaP*TmpArray4(16,4))
     TwoTerms(4) = inv2expP*(TmpArray4(17,3) + alphaP*TmpArray4(17,4))
     TwoTerms(5) = inv2expP*(TmpArray4(19,3) + alphaP*TmpArray4(19,4))
     TwoTerms(6) = inv2expP*(TmpArray4(20,3) + alphaP*TmpArray4(20,4))
     tmpArray6(36,3) = Xpa*tmpArray5(21,3) + alphaXpq*TmpArray5(21,4) + 4*TwoTerms(1)
     tmpArray6(37,3) = Ypa*tmpArray5(21,3) + alphaYpq*TmpArray5(21,4)
     tmpArray6(38,3) = Zpa*tmpArray5(21,3) + alphaZpq*TmpArray5(21,4)
     tmpArray6(39,3) = Xpa*tmpArray5(24,3) + alphaXpq*TmpArray5(24,4) + 2*TwoTerms(2)
     tmpArray6(40,3) = Ypa*tmpArray5(23,3) + alphaYpq*TmpArray5(23,4)
     tmpArray6(41,3) = Xpa*tmpArray5(26,3) + alphaXpq*TmpArray5(26,4) + 2*TwoTerms(3)
     tmpArray6(42,3) = Xpa*tmpArray5(27,3) + alphaXpq*TmpArray5(27,4) + TwoTerms(4)
     tmpArray6(43,3) = Zpa*tmpArray5(24,3) + alphaZpq*TmpArray5(24,4)
     tmpArray6(44,3) = Ypa*tmpArray5(26,3) + alphaYpq*TmpArray5(26,4)
     tmpArray6(45,3) = Xpa*tmpArray5(30,3) + alphaXpq*TmpArray5(30,4) + TwoTerms(6)
     tmpArray6(46,3) = Xpa*tmpArray5(31,3) + alphaXpq*TmpArray5(31,4)
     tmpArray6(47,3) = Xpa*tmpArray5(32,3) + alphaXpq*TmpArray5(32,4)
     tmpArray6(48,3) = Xpa*tmpArray5(33,3) + alphaXpq*TmpArray5(33,4)
     tmpArray6(49,3) = Xpa*tmpArray5(34,3) + alphaXpq*TmpArray5(34,4)
     tmpArray6(50,3) = Xpa*tmpArray5(35,3) + alphaXpq*TmpArray5(35,4)
     tmpArray6(51,3) = Ypa*tmpArray5(31,3) + alphaYpq*TmpArray5(31,4) + 4*TwoTerms(4)
     tmpArray6(52,3) = Zpa*tmpArray5(31,3) + alphaZpq*TmpArray5(31,4)
     tmpArray6(53,3) = Ypa*tmpArray5(33,3) + alphaYpq*TmpArray5(33,4) + 2*TwoTerms(5)
     tmpArray6(54,3) = Ypa*tmpArray5(34,3) + alphaYpq*TmpArray5(34,4) + TwoTerms(6)
     tmpArray6(55,3) = Ypa*tmpArray5(35,3) + alphaYpq*TmpArray5(35,4)
     tmpArray6(56,3) = Zpa*tmpArray5(35,3) + alphaZpq*TmpArray5(35,4) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TMPAuxArray(21) + alphaP*TmpArray5(21,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(24) + alphaP*TmpArray5(24,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(26) + alphaP*TmpArray5(26,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(27) + alphaP*TmpArray5(27,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(30) + alphaP*TmpArray5(30,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(31) + alphaP*TmpArray5(31,2))
     TwoTerms(7) = inv2expP*(TMPAuxArray(33) + alphaP*TmpArray5(33,2))
     TwoTerms(8) = inv2expP*(TMPAuxArray(34) + alphaP*TmpArray5(34,2))
     TwoTerms(9) = inv2expP*(TMPAuxArray(35) + alphaP*TmpArray5(35,2))
     TMPAuxArray(57) = Xpa*TMPAuxArray(36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     TMPAuxArray(58) = Ypa*TMPAuxArray(36) + alphaYpq*TmpArray6(36,2)
     TMPAuxArray(59) = Zpa*TMPAuxArray(36) + alphaZpq*TmpArray6(36,2)
     TMPAuxArray(60) = Xpa*TMPAuxArray(39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     TMPAuxArray(61) = Ypa*TMPAuxArray(38) + alphaYpq*TmpArray6(38,2)
     TMPAuxArray(62) = Xpa*TMPAuxArray(41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     TMPAuxArray(63) = Xpa*TMPAuxArray(42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     TMPAuxArray(64) = Zpa*TMPAuxArray(39) + alphaZpq*TmpArray6(39,2)
     TMPAuxArray(65) = Ypa*TMPAuxArray(41) + alphaYpq*TmpArray6(41,2)
     TMPAuxArray(66) = Xpa*TMPAuxArray(45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     TMPAuxArray(67) = Xpa*TMPAuxArray(46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     TMPAuxArray(68) = Zpa*TMPAuxArray(42) + alphaZpq*TmpArray6(42,2)
     TMPAuxArray(69) = Xpa*TMPAuxArray(48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     TMPAuxArray(70) = Ypa*TMPAuxArray(45) + alphaYpq*TmpArray6(45,2)
     TMPAuxArray(71) = Xpa*TMPAuxArray(50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     TMPAuxArray(72) = Xpa*TMPAuxArray(51) + alphaXpq*TmpArray6(51,2)
     TMPAuxArray(73) = Xpa*TMPAuxArray(52) + alphaXpq*TmpArray6(52,2)
     TMPAuxArray(74) = Xpa*TMPAuxArray(53) + alphaXpq*TmpArray6(53,2)
     TMPAuxArray(75) = Xpa*TMPAuxArray(54) + alphaXpq*TmpArray6(54,2)
     TMPAuxArray(76) = Xpa*TMPAuxArray(55) + alphaXpq*TmpArray6(55,2)
     TMPAuxArray(77) = Xpa*TMPAuxArray(56) + alphaXpq*TmpArray6(56,2)
     TMPAuxArray(78) = Ypa*TMPAuxArray(51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     TMPAuxArray(79) = Zpa*TMPAuxArray(51) + alphaZpq*TmpArray6(51,2)
     TMPAuxArray(80) = Ypa*TMPAuxArray(53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     TMPAuxArray(81) = Ypa*TMPAuxArray(54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     TMPAuxArray(82) = Ypa*TMPAuxArray(55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     TMPAuxArray(83) = Ypa*TMPAuxArray(56) + alphaYpq*TmpArray6(56,2)
     TMPAuxArray(84) = Zpa*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expP*(TmpArray5(21,2) + alphaP*TmpArray5(21,3))
     TwoTerms(2) = inv2expP*(TmpArray5(24,2) + alphaP*TmpArray5(24,3))
     TwoTerms(3) = inv2expP*(TmpArray5(26,2) + alphaP*TmpArray5(26,3))
     TwoTerms(4) = inv2expP*(TmpArray5(27,2) + alphaP*TmpArray5(27,3))
     TwoTerms(5) = inv2expP*(TmpArray5(30,2) + alphaP*TmpArray5(30,3))
     TwoTerms(6) = inv2expP*(TmpArray5(31,2) + alphaP*TmpArray5(31,3))
     TwoTerms(7) = inv2expP*(TmpArray5(33,2) + alphaP*TmpArray5(33,3))
     TwoTerms(8) = inv2expP*(TmpArray5(34,2) + alphaP*TmpArray5(34,3))
     TwoTerms(9) = inv2expP*(TmpArray5(35,2) + alphaP*TmpArray5(35,3))
     tmpArray7(57,2) = Xpa*tmpArray6(36,2) + alphaXpq*TmpArray6(36,3) + 5*TwoTerms(1)
     tmpArray7(58,2) = Ypa*tmpArray6(36,2) + alphaYpq*TmpArray6(36,3)
     tmpArray7(59,2) = Zpa*tmpArray6(36,2) + alphaZpq*TmpArray6(36,3)
     tmpArray7(60,2) = Xpa*tmpArray6(39,2) + alphaXpq*TmpArray6(39,3) + 3*TwoTerms(2)
     tmpArray7(61,2) = Ypa*tmpArray6(38,2) + alphaYpq*TmpArray6(38,3)
     tmpArray7(62,2) = Xpa*tmpArray6(41,2) + alphaXpq*TmpArray6(41,3) + 3*TwoTerms(3)
     tmpArray7(63,2) = Xpa*tmpArray6(42,2) + alphaXpq*TmpArray6(42,3) + 2*TwoTerms(4)
     tmpArray7(64,2) = Zpa*tmpArray6(39,2) + alphaZpq*TmpArray6(39,3)
     tmpArray7(65,2) = Ypa*tmpArray6(41,2) + alphaYpq*TmpArray6(41,3)
     tmpArray7(66,2) = Xpa*tmpArray6(45,2) + alphaXpq*TmpArray6(45,3) + 2*TwoTerms(5)
     tmpArray7(67,2) = Xpa*tmpArray6(46,2) + alphaXpq*TmpArray6(46,3) + TwoTerms(6)
     tmpArray7(68,2) = Zpa*tmpArray6(42,2) + alphaZpq*TmpArray6(42,3)
     tmpArray7(69,2) = Xpa*tmpArray6(48,2) + alphaXpq*TmpArray6(48,3) + TwoTerms(7)
     tmpArray7(70,2) = Ypa*tmpArray6(45,2) + alphaYpq*TmpArray6(45,3)
     tmpArray7(71,2) = Xpa*tmpArray6(50,2) + alphaXpq*TmpArray6(50,3) + TwoTerms(9)
     tmpArray7(72,2) = Xpa*tmpArray6(51,2) + alphaXpq*TmpArray6(51,3)
     tmpArray7(73,2) = Xpa*tmpArray6(52,2) + alphaXpq*TmpArray6(52,3)
     tmpArray7(74,2) = Xpa*tmpArray6(53,2) + alphaXpq*TmpArray6(53,3)
     tmpArray7(75,2) = Xpa*tmpArray6(54,2) + alphaXpq*TmpArray6(54,3)
     tmpArray7(76,2) = Xpa*tmpArray6(55,2) + alphaXpq*TmpArray6(55,3)
     tmpArray7(77,2) = Xpa*tmpArray6(56,2) + alphaXpq*TmpArray6(56,3)
     tmpArray7(78,2) = Ypa*tmpArray6(51,2) + alphaYpq*TmpArray6(51,3) + 5*TwoTerms(6)
     tmpArray7(79,2) = Zpa*tmpArray6(51,2) + alphaZpq*TmpArray6(51,3)
     tmpArray7(80,2) = Ypa*tmpArray6(53,2) + alphaYpq*TmpArray6(53,3) + 3*TwoTerms(7)
     tmpArray7(81,2) = Ypa*tmpArray6(54,2) + alphaYpq*TmpArray6(54,3) + 2*TwoTerms(8)
     tmpArray7(82,2) = Ypa*tmpArray6(55,2) + alphaYpq*TmpArray6(55,3) + TwoTerms(9)
     tmpArray7(83,2) = Ypa*tmpArray6(56,2) + alphaYpq*TmpArray6(56,3)
     tmpArray7(84,2) = Zpa*tmpArray6(56,2) + alphaZpq*TmpArray6(56,3) + 5*TwoTerms(9)
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
      AuxArray(iTUV,IP) = AuxArray(iTUV,IP) + TMPAuxarray(iTUV)
     enddo
     AuxArray(85,IP) = AuxArray(85,IP) + Xpa*TMPAuxArray(57) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     AuxArray(86,IP) = AuxArray(86,IP) + Ypa*TMPAuxArray(57) + alphaYpq*TmpArray7(57,2)
     AuxArray(87,IP) = AuxArray(87,IP) + Zpa*TMPAuxArray(57) + alphaZpq*TmpArray7(57,2)
     AuxArray(88,IP) = AuxArray(88,IP) + Xpa*TMPAuxArray(60) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     AuxArray(89,IP) = AuxArray(89,IP) + Ypa*TMPAuxArray(59) + alphaYpq*TmpArray7(59,2)
     AuxArray(90,IP) = AuxArray(90,IP) + Xpa*TMPAuxArray(62) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     AuxArray(91,IP) = AuxArray(91,IP) + Xpa*TMPAuxArray(63) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     AuxArray(92,IP) = AuxArray(92,IP) + Zpa*TMPAuxArray(60) + alphaZpq*TmpArray7(60,2)
     AuxArray(93,IP) = AuxArray(93,IP) + Ypa*TMPAuxArray(62) + alphaYpq*TmpArray7(62,2)
     AuxArray(94,IP) = AuxArray(94,IP) + Xpa*TMPAuxArray(66) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     AuxArray(95,IP) = AuxArray(95,IP) + Xpa*TMPAuxArray(67) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     AuxArray(96,IP) = AuxArray(96,IP) + Zpa*TMPAuxArray(63) + alphaZpq*TmpArray7(63,2)
     AuxArray(97,IP) = AuxArray(97,IP) + Xpa*TMPAuxArray(69) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     AuxArray(98,IP) = AuxArray(98,IP) + Ypa*TMPAuxArray(66) + alphaYpq*TmpArray7(66,2)
     AuxArray(99,IP) = AuxArray(99,IP) + Xpa*TMPAuxArray(71) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     AuxArray(100,IP) = AuxArray(100,IP) + Xpa*TMPAuxArray(72) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     AuxArray(101,IP) = AuxArray(101,IP) + Zpa*TMPAuxArray(67) + alphaZpq*TmpArray7(67,2)
     AuxArray(102,IP) = AuxArray(102,IP) + Xpa*TMPAuxArray(74) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     AuxArray(103,IP) = AuxArray(103,IP) + Xpa*TMPAuxArray(75) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     AuxArray(104,IP) = AuxArray(104,IP) + Ypa*TMPAuxArray(71) + alphaYpq*TmpArray7(71,2)
     AuxArray(105,IP) = AuxArray(105,IP) + Xpa*TMPAuxArray(77) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     AuxArray(106,IP) = AuxArray(106,IP) + Xpa*TMPAuxArray(78) + alphaXpq*TmpArray7(78,2)
     AuxArray(107,IP) = AuxArray(107,IP) + Xpa*TMPAuxArray(79) + alphaXpq*TmpArray7(79,2)
     AuxArray(108,IP) = AuxArray(108,IP) + Xpa*TMPAuxArray(80) + alphaXpq*TmpArray7(80,2)
     AuxArray(109,IP) = AuxArray(109,IP) + Xpa*TMPAuxArray(81) + alphaXpq*TmpArray7(81,2)
     AuxArray(110,IP) = AuxArray(110,IP) + Xpa*TMPAuxArray(82) + alphaXpq*TmpArray7(82,2)
     AuxArray(111,IP) = AuxArray(111,IP) + Xpa*TMPAuxArray(83) + alphaXpq*TmpArray7(83,2)
     AuxArray(112,IP) = AuxArray(112,IP) + Xpa*TMPAuxArray(84) + alphaXpq*TmpArray7(84,2)
     AuxArray(113,IP) = AuxArray(113,IP) + Ypa*TMPAuxArray(78) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     AuxArray(114,IP) = AuxArray(114,IP) + Zpa*TMPAuxArray(78) + alphaZpq*TmpArray7(78,2)
     AuxArray(115,IP) = AuxArray(115,IP) + Ypa*TMPAuxArray(80) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     AuxArray(116,IP) = AuxArray(116,IP) + Ypa*TMPAuxArray(81) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     AuxArray(117,IP) = AuxArray(117,IP) + Ypa*TMPAuxArray(82) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     AuxArray(118,IP) = AuxArray(118,IP) + Ypa*TMPAuxArray(83) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     AuxArray(119,IP) = AuxArray(119,IP) + Ypa*TMPAuxArray(84) + alphaYpq*TmpArray7(84,2)
     AuxArray(120,IP) = AuxArray(120,IP) + Zpa*TMPAuxArray(84) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
    ENDDO
   ENDDO
  ENDDO
 end subroutine

subroutine VerticalRecurrenceCPUSegQ8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  REAL(REALK),intent(in) :: TABFJW(0:11,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
  real(realk),intent(in) :: Acenter(3)
  real(realk),intent(inout) :: AUXarray(  165,nPrimP*nPasses)
  !local variables
  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV
  real(realk) :: TMPAUXarray(  120)
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
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
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  Ax = -Acenter(1)
  Ay = -Acenter(2)
  Az = -Acenter(3)
  DO iPassQ = 1,nPasses
   iP = (iPassQ-1)*nPrimP
   DO iPrimP=1, nPrimP
    iP = iP + 1
    DO iTUV=1,  165
     AUXarray(iTUV,iP)=0.0E0_realk
    ENDDO
    Pexpfac = PpreExpFac(iPrimP)
    mPX = -Pcent(1,iPrimP)
    mPY = -Pcent(2,iPrimP)
    mPZ = -Pcent(3,iPrimP)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP) + Ax
    Ypa = Pcent(2,iPrimP) + Ay
    Zpa = Pcent(3,iPrimP) + Az
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
     tmpArray2(2,5) = Xpa*tmpArray1(1,5) + alphaXpq*TmpArray1(1,6)
     tmpArray2(3,5) = Ypa*tmpArray1(1,5) + alphaYpq*TmpArray1(1,6)
     tmpArray2(4,5) = Zpa*tmpArray1(1,5) + alphaZpq*TmpArray1(1,6)
     tmpArray2(2,6) = Xpa*tmpArray1(1,6) + alphaXpq*TmpArray1(1,7)
     tmpArray2(3,6) = Ypa*tmpArray1(1,6) + alphaYpq*TmpArray1(1,7)
     tmpArray2(4,6) = Zpa*tmpArray1(1,6) + alphaZpq*TmpArray1(1,7)
     tmpArray2(2,7) = Xpa*tmpArray1(1,7) + alphaXpq*TmpArray1(1,8)
     tmpArray2(3,7) = Ypa*tmpArray1(1,7) + alphaYpq*TmpArray1(1,8)
     tmpArray2(4,7) = Zpa*tmpArray1(1,7) + alphaZpq*TmpArray1(1,8)
     tmpArray2(2,8) = Xpa*tmpArray1(1,8) + alphaXpq*TmpArray1(1,9)
     tmpArray2(3,8) = Ypa*tmpArray1(1,8) + alphaYpq*TmpArray1(1,9)
     tmpArray2(4,8) = Zpa*tmpArray1(1,8) + alphaZpq*TmpArray1(1,9)
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
     TwoTerms(1) = inv2expP*(TmpArray1(1,4) + alphaP*TmpArray1(1,5))
     tmpArray3(5,4) = Xpa*tmpArray2(2,4) + alphaXpq*TmpArray2(2,5) + TwoTerms(1)
     tmpArray3(6,4) = Xpa*tmpArray2(3,4) + alphaXpq*TmpArray2(3,5)
     tmpArray3(7,4) = Xpa*tmpArray2(4,4) + alphaXpq*TmpArray2(4,5)
     tmpArray3(8,4) = Ypa*tmpArray2(3,4) + alphaYpq*TmpArray2(3,5) + TwoTerms(1)
     tmpArray3(9,4) = Ypa*tmpArray2(4,4) + alphaYpq*TmpArray2(4,5)
     tmpArray3(10,4) = Zpa*tmpArray2(4,4) + alphaZpq*TmpArray2(4,5) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,5) + alphaP*TmpArray1(1,6))
     tmpArray3(5,5) = Xpa*tmpArray2(2,5) + alphaXpq*TmpArray2(2,6) + TwoTerms(1)
     tmpArray3(6,5) = Xpa*tmpArray2(3,5) + alphaXpq*TmpArray2(3,6)
     tmpArray3(7,5) = Xpa*tmpArray2(4,5) + alphaXpq*TmpArray2(4,6)
     tmpArray3(8,5) = Ypa*tmpArray2(3,5) + alphaYpq*TmpArray2(3,6) + TwoTerms(1)
     tmpArray3(9,5) = Ypa*tmpArray2(4,5) + alphaYpq*TmpArray2(4,6)
     tmpArray3(10,5) = Zpa*tmpArray2(4,5) + alphaZpq*TmpArray2(4,6) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,6) + alphaP*TmpArray1(1,7))
     tmpArray3(5,6) = Xpa*tmpArray2(2,6) + alphaXpq*TmpArray2(2,7) + TwoTerms(1)
     tmpArray3(6,6) = Xpa*tmpArray2(3,6) + alphaXpq*TmpArray2(3,7)
     tmpArray3(7,6) = Xpa*tmpArray2(4,6) + alphaXpq*TmpArray2(4,7)
     tmpArray3(8,6) = Ypa*tmpArray2(3,6) + alphaYpq*TmpArray2(3,7) + TwoTerms(1)
     tmpArray3(9,6) = Ypa*tmpArray2(4,6) + alphaYpq*TmpArray2(4,7)
     tmpArray3(10,6) = Zpa*tmpArray2(4,6) + alphaZpq*TmpArray2(4,7) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,7) + alphaP*TmpArray1(1,8))
     tmpArray3(5,7) = Xpa*tmpArray2(2,7) + alphaXpq*TmpArray2(2,8) + TwoTerms(1)
     tmpArray3(6,7) = Xpa*tmpArray2(3,7) + alphaXpq*TmpArray2(3,8)
     tmpArray3(7,7) = Xpa*tmpArray2(4,7) + alphaXpq*TmpArray2(4,8)
     tmpArray3(8,7) = Ypa*tmpArray2(3,7) + alphaYpq*TmpArray2(3,8) + TwoTerms(1)
     tmpArray3(9,7) = Ypa*tmpArray2(4,7) + alphaYpq*TmpArray2(4,8)
     tmpArray3(10,7) = Zpa*tmpArray2(4,7) + alphaZpq*TmpArray2(4,8) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(TmpArray2(2,3) + alphaP*TmpArray2(2,4))
     TwoTerms(2) = inv2expP*(TmpArray2(3,3) + alphaP*TmpArray2(3,4))
     TwoTerms(3) = inv2expP*(TmpArray2(4,3) + alphaP*TmpArray2(4,4))
     tmpArray4(11,3) = Xpa*tmpArray3(5,3) + alphaXpq*TmpArray3(5,4) + 2*TwoTerms(1)
     tmpArray4(12,3) = Ypa*tmpArray3(5,3) + alphaYpq*TmpArray3(5,4)
     tmpArray4(13,3) = Zpa*tmpArray3(5,3) + alphaZpq*TmpArray3(5,4)
     tmpArray4(14,3) = Xpa*tmpArray3(8,3) + alphaXpq*TmpArray3(8,4)
     tmpArray4(15,3) = Xpa*tmpArray3(9,3) + alphaXpq*TmpArray3(9,4)
     tmpArray4(16,3) = Xpa*tmpArray3(10,3) + alphaXpq*TmpArray3(10,4)
     tmpArray4(17,3) = Ypa*tmpArray3(8,3) + alphaYpq*TmpArray3(8,4) + 2*TwoTerms(2)
     tmpArray4(18,3) = Zpa*tmpArray3(8,3) + alphaZpq*TmpArray3(8,4)
     tmpArray4(19,3) = Ypa*tmpArray3(10,3) + alphaYpq*TmpArray3(10,4)
     tmpArray4(20,3) = Zpa*tmpArray3(10,3) + alphaZpq*TmpArray3(10,4) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,4) + alphaP*TmpArray2(2,5))
     TwoTerms(2) = inv2expP*(TmpArray2(3,4) + alphaP*TmpArray2(3,5))
     TwoTerms(3) = inv2expP*(TmpArray2(4,4) + alphaP*TmpArray2(4,5))
     tmpArray4(11,4) = Xpa*tmpArray3(5,4) + alphaXpq*TmpArray3(5,5) + 2*TwoTerms(1)
     tmpArray4(12,4) = Ypa*tmpArray3(5,4) + alphaYpq*TmpArray3(5,5)
     tmpArray4(13,4) = Zpa*tmpArray3(5,4) + alphaZpq*TmpArray3(5,5)
     tmpArray4(14,4) = Xpa*tmpArray3(8,4) + alphaXpq*TmpArray3(8,5)
     tmpArray4(15,4) = Xpa*tmpArray3(9,4) + alphaXpq*TmpArray3(9,5)
     tmpArray4(16,4) = Xpa*tmpArray3(10,4) + alphaXpq*TmpArray3(10,5)
     tmpArray4(17,4) = Ypa*tmpArray3(8,4) + alphaYpq*TmpArray3(8,5) + 2*TwoTerms(2)
     tmpArray4(18,4) = Zpa*tmpArray3(8,4) + alphaZpq*TmpArray3(8,5)
     tmpArray4(19,4) = Ypa*tmpArray3(10,4) + alphaYpq*TmpArray3(10,5)
     tmpArray4(20,4) = Zpa*tmpArray3(10,4) + alphaZpq*TmpArray3(10,5) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,5) + alphaP*TmpArray2(2,6))
     TwoTerms(2) = inv2expP*(TmpArray2(3,5) + alphaP*TmpArray2(3,6))
     TwoTerms(3) = inv2expP*(TmpArray2(4,5) + alphaP*TmpArray2(4,6))
     tmpArray4(11,5) = Xpa*tmpArray3(5,5) + alphaXpq*TmpArray3(5,6) + 2*TwoTerms(1)
     tmpArray4(12,5) = Ypa*tmpArray3(5,5) + alphaYpq*TmpArray3(5,6)
     tmpArray4(13,5) = Zpa*tmpArray3(5,5) + alphaZpq*TmpArray3(5,6)
     tmpArray4(14,5) = Xpa*tmpArray3(8,5) + alphaXpq*TmpArray3(8,6)
     tmpArray4(15,5) = Xpa*tmpArray3(9,5) + alphaXpq*TmpArray3(9,6)
     tmpArray4(16,5) = Xpa*tmpArray3(10,5) + alphaXpq*TmpArray3(10,6)
     tmpArray4(17,5) = Ypa*tmpArray3(8,5) + alphaYpq*TmpArray3(8,6) + 2*TwoTerms(2)
     tmpArray4(18,5) = Zpa*tmpArray3(8,5) + alphaZpq*TmpArray3(8,6)
     tmpArray4(19,5) = Ypa*tmpArray3(10,5) + alphaYpq*TmpArray3(10,6)
     tmpArray4(20,5) = Zpa*tmpArray3(10,5) + alphaZpq*TmpArray3(10,6) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray2(2,6) + alphaP*TmpArray2(2,7))
     TwoTerms(2) = inv2expP*(TmpArray2(3,6) + alphaP*TmpArray2(3,7))
     TwoTerms(3) = inv2expP*(TmpArray2(4,6) + alphaP*TmpArray2(4,7))
     tmpArray4(11,6) = Xpa*tmpArray3(5,6) + alphaXpq*TmpArray3(5,7) + 2*TwoTerms(1)
     tmpArray4(12,6) = Ypa*tmpArray3(5,6) + alphaYpq*TmpArray3(5,7)
     tmpArray4(13,6) = Zpa*tmpArray3(5,6) + alphaZpq*TmpArray3(5,7)
     tmpArray4(14,6) = Xpa*tmpArray3(8,6) + alphaXpq*TmpArray3(8,7)
     tmpArray4(15,6) = Xpa*tmpArray3(9,6) + alphaXpq*TmpArray3(9,7)
     tmpArray4(16,6) = Xpa*tmpArray3(10,6) + alphaXpq*TmpArray3(10,7)
     tmpArray4(17,6) = Ypa*tmpArray3(8,6) + alphaYpq*TmpArray3(8,7) + 2*TwoTerms(2)
     tmpArray4(18,6) = Zpa*tmpArray3(8,6) + alphaZpq*TmpArray3(8,7)
     tmpArray4(19,6) = Ypa*tmpArray3(10,6) + alphaYpq*TmpArray3(10,7)
     tmpArray4(20,6) = Zpa*tmpArray3(10,6) + alphaZpq*TmpArray3(10,7) + 2*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(5) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(8) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(10) + alphaP*TmpArray3(10,2))
     TMPAuxArray(21) = Xpa*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Ypa*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zpa*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xpa*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Ypa*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xpa*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xpa*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xpa*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xpa*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xpa*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Ypa*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zpa*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Ypa*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Ypa*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zpa*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,2) + alphaP*TmpArray3(5,3))
     TwoTerms(2) = inv2expP*(TmpArray3(8,2) + alphaP*TmpArray3(8,3))
     TwoTerms(3) = inv2expP*(TmpArray3(10,2) + alphaP*TmpArray3(10,3))
     tmpArray5(21,2) = Xpa*tmpArray4(11,2) + alphaXpq*TmpArray4(11,3) + 3*TwoTerms(1)
     tmpArray5(22,2) = Ypa*tmpArray4(11,2) + alphaYpq*TmpArray4(11,3)
     tmpArray5(23,2) = Zpa*tmpArray4(11,2) + alphaZpq*TmpArray4(11,3)
     tmpArray5(24,2) = Xpa*tmpArray4(14,2) + alphaXpq*TmpArray4(14,3) + TwoTerms(2)
     tmpArray5(25,2) = Ypa*tmpArray4(13,2) + alphaYpq*TmpArray4(13,3)
     tmpArray5(26,2) = Xpa*tmpArray4(16,2) + alphaXpq*TmpArray4(16,3) + TwoTerms(3)
     tmpArray5(27,2) = Xpa*tmpArray4(17,2) + alphaXpq*TmpArray4(17,3)
     tmpArray5(28,2) = Xpa*tmpArray4(18,2) + alphaXpq*TmpArray4(18,3)
     tmpArray5(29,2) = Xpa*tmpArray4(19,2) + alphaXpq*TmpArray4(19,3)
     tmpArray5(30,2) = Xpa*tmpArray4(20,2) + alphaXpq*TmpArray4(20,3)
     tmpArray5(31,2) = Ypa*tmpArray4(17,2) + alphaYpq*TmpArray4(17,3) + 3*TwoTerms(2)
     tmpArray5(32,2) = Zpa*tmpArray4(17,2) + alphaZpq*TmpArray4(17,3)
     tmpArray5(33,2) = Ypa*tmpArray4(19,2) + alphaYpq*TmpArray4(19,3) + TwoTerms(3)
     tmpArray5(34,2) = Ypa*tmpArray4(20,2) + alphaYpq*TmpArray4(20,3)
     tmpArray5(35,2) = Zpa*tmpArray4(20,2) + alphaZpq*TmpArray4(20,3) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,3) + alphaP*TmpArray3(5,4))
     TwoTerms(2) = inv2expP*(TmpArray3(8,3) + alphaP*TmpArray3(8,4))
     TwoTerms(3) = inv2expP*(TmpArray3(10,3) + alphaP*TmpArray3(10,4))
     tmpArray5(21,3) = Xpa*tmpArray4(11,3) + alphaXpq*TmpArray4(11,4) + 3*TwoTerms(1)
     tmpArray5(22,3) = Ypa*tmpArray4(11,3) + alphaYpq*TmpArray4(11,4)
     tmpArray5(23,3) = Zpa*tmpArray4(11,3) + alphaZpq*TmpArray4(11,4)
     tmpArray5(24,3) = Xpa*tmpArray4(14,3) + alphaXpq*TmpArray4(14,4) + TwoTerms(2)
     tmpArray5(25,3) = Ypa*tmpArray4(13,3) + alphaYpq*TmpArray4(13,4)
     tmpArray5(26,3) = Xpa*tmpArray4(16,3) + alphaXpq*TmpArray4(16,4) + TwoTerms(3)
     tmpArray5(27,3) = Xpa*tmpArray4(17,3) + alphaXpq*TmpArray4(17,4)
     tmpArray5(28,3) = Xpa*tmpArray4(18,3) + alphaXpq*TmpArray4(18,4)
     tmpArray5(29,3) = Xpa*tmpArray4(19,3) + alphaXpq*TmpArray4(19,4)
     tmpArray5(30,3) = Xpa*tmpArray4(20,3) + alphaXpq*TmpArray4(20,4)
     tmpArray5(31,3) = Ypa*tmpArray4(17,3) + alphaYpq*TmpArray4(17,4) + 3*TwoTerms(2)
     tmpArray5(32,3) = Zpa*tmpArray4(17,3) + alphaZpq*TmpArray4(17,4)
     tmpArray5(33,3) = Ypa*tmpArray4(19,3) + alphaYpq*TmpArray4(19,4) + TwoTerms(3)
     tmpArray5(34,3) = Ypa*tmpArray4(20,3) + alphaYpq*TmpArray4(20,4)
     tmpArray5(35,3) = Zpa*tmpArray4(20,3) + alphaZpq*TmpArray4(20,4) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,4) + alphaP*TmpArray3(5,5))
     TwoTerms(2) = inv2expP*(TmpArray3(8,4) + alphaP*TmpArray3(8,5))
     TwoTerms(3) = inv2expP*(TmpArray3(10,4) + alphaP*TmpArray3(10,5))
     tmpArray5(21,4) = Xpa*tmpArray4(11,4) + alphaXpq*TmpArray4(11,5) + 3*TwoTerms(1)
     tmpArray5(22,4) = Ypa*tmpArray4(11,4) + alphaYpq*TmpArray4(11,5)
     tmpArray5(23,4) = Zpa*tmpArray4(11,4) + alphaZpq*TmpArray4(11,5)
     tmpArray5(24,4) = Xpa*tmpArray4(14,4) + alphaXpq*TmpArray4(14,5) + TwoTerms(2)
     tmpArray5(25,4) = Ypa*tmpArray4(13,4) + alphaYpq*TmpArray4(13,5)
     tmpArray5(26,4) = Xpa*tmpArray4(16,4) + alphaXpq*TmpArray4(16,5) + TwoTerms(3)
     tmpArray5(27,4) = Xpa*tmpArray4(17,4) + alphaXpq*TmpArray4(17,5)
     tmpArray5(28,4) = Xpa*tmpArray4(18,4) + alphaXpq*TmpArray4(18,5)
     tmpArray5(29,4) = Xpa*tmpArray4(19,4) + alphaXpq*TmpArray4(19,5)
     tmpArray5(30,4) = Xpa*tmpArray4(20,4) + alphaXpq*TmpArray4(20,5)
     tmpArray5(31,4) = Ypa*tmpArray4(17,4) + alphaYpq*TmpArray4(17,5) + 3*TwoTerms(2)
     tmpArray5(32,4) = Zpa*tmpArray4(17,4) + alphaZpq*TmpArray4(17,5)
     tmpArray5(33,4) = Ypa*tmpArray4(19,4) + alphaYpq*TmpArray4(19,5) + TwoTerms(3)
     tmpArray5(34,4) = Ypa*tmpArray4(20,4) + alphaYpq*TmpArray4(20,5)
     tmpArray5(35,4) = Zpa*tmpArray4(20,4) + alphaZpq*TmpArray4(20,5) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TmpArray3(5,5) + alphaP*TmpArray3(5,6))
     TwoTerms(2) = inv2expP*(TmpArray3(8,5) + alphaP*TmpArray3(8,6))
     TwoTerms(3) = inv2expP*(TmpArray3(10,5) + alphaP*TmpArray3(10,6))
     tmpArray5(21,5) = Xpa*tmpArray4(11,5) + alphaXpq*TmpArray4(11,6) + 3*TwoTerms(1)
     tmpArray5(22,5) = Ypa*tmpArray4(11,5) + alphaYpq*TmpArray4(11,6)
     tmpArray5(23,5) = Zpa*tmpArray4(11,5) + alphaZpq*TmpArray4(11,6)
     tmpArray5(24,5) = Xpa*tmpArray4(14,5) + alphaXpq*TmpArray4(14,6) + TwoTerms(2)
     tmpArray5(25,5) = Ypa*tmpArray4(13,5) + alphaYpq*TmpArray4(13,6)
     tmpArray5(26,5) = Xpa*tmpArray4(16,5) + alphaXpq*TmpArray4(16,6) + TwoTerms(3)
     tmpArray5(27,5) = Xpa*tmpArray4(17,5) + alphaXpq*TmpArray4(17,6)
     tmpArray5(28,5) = Xpa*tmpArray4(18,5) + alphaXpq*TmpArray4(18,6)
     tmpArray5(29,5) = Xpa*tmpArray4(19,5) + alphaXpq*TmpArray4(19,6)
     tmpArray5(30,5) = Xpa*tmpArray4(20,5) + alphaXpq*TmpArray4(20,6)
     tmpArray5(31,5) = Ypa*tmpArray4(17,5) + alphaYpq*TmpArray4(17,6) + 3*TwoTerms(2)
     tmpArray5(32,5) = Zpa*tmpArray4(17,5) + alphaZpq*TmpArray4(17,6)
     tmpArray5(33,5) = Ypa*tmpArray4(19,5) + alphaYpq*TmpArray4(19,6) + TwoTerms(3)
     tmpArray5(34,5) = Ypa*tmpArray4(20,5) + alphaYpq*TmpArray4(20,6)
     tmpArray5(35,5) = Zpa*tmpArray4(20,5) + alphaZpq*TmpArray4(20,6) + 3*TwoTerms(3)
     TwoTerms(1) = inv2expP*(TMPAuxArray(11) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(14) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(16) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(17) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(19) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(20) + alphaP*TmpArray4(20,2))
     TMPAuxArray(36) = Xpa*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     TMPAuxArray(37) = Ypa*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     TMPAuxArray(38) = Zpa*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     TMPAuxArray(39) = Xpa*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     TMPAuxArray(40) = Ypa*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     TMPAuxArray(41) = Xpa*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     TMPAuxArray(42) = Xpa*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     TMPAuxArray(43) = Zpa*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     TMPAuxArray(44) = Ypa*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     TMPAuxArray(45) = Xpa*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     TMPAuxArray(46) = Xpa*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     TMPAuxArray(47) = Xpa*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     TMPAuxArray(48) = Xpa*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     TMPAuxArray(49) = Xpa*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     TMPAuxArray(50) = Xpa*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     TMPAuxArray(51) = Ypa*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     TMPAuxArray(52) = Zpa*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     TMPAuxArray(53) = Ypa*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     TMPAuxArray(54) = Ypa*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     TMPAuxArray(55) = Ypa*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     TMPAuxArray(56) = Zpa*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TmpArray4(11,2) + alphaP*TmpArray4(11,3))
     TwoTerms(2) = inv2expP*(TmpArray4(14,2) + alphaP*TmpArray4(14,3))
     TwoTerms(3) = inv2expP*(TmpArray4(16,2) + alphaP*TmpArray4(16,3))
     TwoTerms(4) = inv2expP*(TmpArray4(17,2) + alphaP*TmpArray4(17,3))
     TwoTerms(5) = inv2expP*(TmpArray4(19,2) + alphaP*TmpArray4(19,3))
     TwoTerms(6) = inv2expP*(TmpArray4(20,2) + alphaP*TmpArray4(20,3))
     tmpArray6(36,2) = Xpa*tmpArray5(21,2) + alphaXpq*TmpArray5(21,3) + 4*TwoTerms(1)
     tmpArray6(37,2) = Ypa*tmpArray5(21,2) + alphaYpq*TmpArray5(21,3)
     tmpArray6(38,2) = Zpa*tmpArray5(21,2) + alphaZpq*TmpArray5(21,3)
     tmpArray6(39,2) = Xpa*tmpArray5(24,2) + alphaXpq*TmpArray5(24,3) + 2*TwoTerms(2)
     tmpArray6(40,2) = Ypa*tmpArray5(23,2) + alphaYpq*TmpArray5(23,3)
     tmpArray6(41,2) = Xpa*tmpArray5(26,2) + alphaXpq*TmpArray5(26,3) + 2*TwoTerms(3)
     tmpArray6(42,2) = Xpa*tmpArray5(27,2) + alphaXpq*TmpArray5(27,3) + TwoTerms(4)
     tmpArray6(43,2) = Zpa*tmpArray5(24,2) + alphaZpq*TmpArray5(24,3)
     tmpArray6(44,2) = Ypa*tmpArray5(26,2) + alphaYpq*TmpArray5(26,3)
     tmpArray6(45,2) = Xpa*tmpArray5(30,2) + alphaXpq*TmpArray5(30,3) + TwoTerms(6)
     tmpArray6(46,2) = Xpa*tmpArray5(31,2) + alphaXpq*TmpArray5(31,3)
     tmpArray6(47,2) = Xpa*tmpArray5(32,2) + alphaXpq*TmpArray5(32,3)
     tmpArray6(48,2) = Xpa*tmpArray5(33,2) + alphaXpq*TmpArray5(33,3)
     tmpArray6(49,2) = Xpa*tmpArray5(34,2) + alphaXpq*TmpArray5(34,3)
     tmpArray6(50,2) = Xpa*tmpArray5(35,2) + alphaXpq*TmpArray5(35,3)
     tmpArray6(51,2) = Ypa*tmpArray5(31,2) + alphaYpq*TmpArray5(31,3) + 4*TwoTerms(4)
     tmpArray6(52,2) = Zpa*tmpArray5(31,2) + alphaZpq*TmpArray5(31,3)
     tmpArray6(53,2) = Ypa*tmpArray5(33,2) + alphaYpq*TmpArray5(33,3) + 2*TwoTerms(5)
     tmpArray6(54,2) = Ypa*tmpArray5(34,2) + alphaYpq*TmpArray5(34,3) + TwoTerms(6)
     tmpArray6(55,2) = Ypa*tmpArray5(35,2) + alphaYpq*TmpArray5(35,3)
     tmpArray6(56,2) = Zpa*tmpArray5(35,2) + alphaZpq*TmpArray5(35,3) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TmpArray4(11,3) + alphaP*TmpArray4(11,4))
     TwoTerms(2) = inv2expP*(TmpArray4(14,3) + alphaP*TmpArray4(14,4))
     TwoTerms(3) = inv2expP*(TmpArray4(16,3) + alphaP*TmpArray4(16,4))
     TwoTerms(4) = inv2expP*(TmpArray4(17,3) + alphaP*TmpArray4(17,4))
     TwoTerms(5) = inv2expP*(TmpArray4(19,3) + alphaP*TmpArray4(19,4))
     TwoTerms(6) = inv2expP*(TmpArray4(20,3) + alphaP*TmpArray4(20,4))
     tmpArray6(36,3) = Xpa*tmpArray5(21,3) + alphaXpq*TmpArray5(21,4) + 4*TwoTerms(1)
     tmpArray6(37,3) = Ypa*tmpArray5(21,3) + alphaYpq*TmpArray5(21,4)
     tmpArray6(38,3) = Zpa*tmpArray5(21,3) + alphaZpq*TmpArray5(21,4)
     tmpArray6(39,3) = Xpa*tmpArray5(24,3) + alphaXpq*TmpArray5(24,4) + 2*TwoTerms(2)
     tmpArray6(40,3) = Ypa*tmpArray5(23,3) + alphaYpq*TmpArray5(23,4)
     tmpArray6(41,3) = Xpa*tmpArray5(26,3) + alphaXpq*TmpArray5(26,4) + 2*TwoTerms(3)
     tmpArray6(42,3) = Xpa*tmpArray5(27,3) + alphaXpq*TmpArray5(27,4) + TwoTerms(4)
     tmpArray6(43,3) = Zpa*tmpArray5(24,3) + alphaZpq*TmpArray5(24,4)
     tmpArray6(44,3) = Ypa*tmpArray5(26,3) + alphaYpq*TmpArray5(26,4)
     tmpArray6(45,3) = Xpa*tmpArray5(30,3) + alphaXpq*TmpArray5(30,4) + TwoTerms(6)
     tmpArray6(46,3) = Xpa*tmpArray5(31,3) + alphaXpq*TmpArray5(31,4)
     tmpArray6(47,3) = Xpa*tmpArray5(32,3) + alphaXpq*TmpArray5(32,4)
     tmpArray6(48,3) = Xpa*tmpArray5(33,3) + alphaXpq*TmpArray5(33,4)
     tmpArray6(49,3) = Xpa*tmpArray5(34,3) + alphaXpq*TmpArray5(34,4)
     tmpArray6(50,3) = Xpa*tmpArray5(35,3) + alphaXpq*TmpArray5(35,4)
     tmpArray6(51,3) = Ypa*tmpArray5(31,3) + alphaYpq*TmpArray5(31,4) + 4*TwoTerms(4)
     tmpArray6(52,3) = Zpa*tmpArray5(31,3) + alphaZpq*TmpArray5(31,4)
     tmpArray6(53,3) = Ypa*tmpArray5(33,3) + alphaYpq*TmpArray5(33,4) + 2*TwoTerms(5)
     tmpArray6(54,3) = Ypa*tmpArray5(34,3) + alphaYpq*TmpArray5(34,4) + TwoTerms(6)
     tmpArray6(55,3) = Ypa*tmpArray5(35,3) + alphaYpq*TmpArray5(35,4)
     tmpArray6(56,3) = Zpa*tmpArray5(35,3) + alphaZpq*TmpArray5(35,4) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TmpArray4(11,4) + alphaP*TmpArray4(11,5))
     TwoTerms(2) = inv2expP*(TmpArray4(14,4) + alphaP*TmpArray4(14,5))
     TwoTerms(3) = inv2expP*(TmpArray4(16,4) + alphaP*TmpArray4(16,5))
     TwoTerms(4) = inv2expP*(TmpArray4(17,4) + alphaP*TmpArray4(17,5))
     TwoTerms(5) = inv2expP*(TmpArray4(19,4) + alphaP*TmpArray4(19,5))
     TwoTerms(6) = inv2expP*(TmpArray4(20,4) + alphaP*TmpArray4(20,5))
     tmpArray6(36,4) = Xpa*tmpArray5(21,4) + alphaXpq*TmpArray5(21,5) + 4*TwoTerms(1)
     tmpArray6(37,4) = Ypa*tmpArray5(21,4) + alphaYpq*TmpArray5(21,5)
     tmpArray6(38,4) = Zpa*tmpArray5(21,4) + alphaZpq*TmpArray5(21,5)
     tmpArray6(39,4) = Xpa*tmpArray5(24,4) + alphaXpq*TmpArray5(24,5) + 2*TwoTerms(2)
     tmpArray6(40,4) = Ypa*tmpArray5(23,4) + alphaYpq*TmpArray5(23,5)
     tmpArray6(41,4) = Xpa*tmpArray5(26,4) + alphaXpq*TmpArray5(26,5) + 2*TwoTerms(3)
     tmpArray6(42,4) = Xpa*tmpArray5(27,4) + alphaXpq*TmpArray5(27,5) + TwoTerms(4)
     tmpArray6(43,4) = Zpa*tmpArray5(24,4) + alphaZpq*TmpArray5(24,5)
     tmpArray6(44,4) = Ypa*tmpArray5(26,4) + alphaYpq*TmpArray5(26,5)
     tmpArray6(45,4) = Xpa*tmpArray5(30,4) + alphaXpq*TmpArray5(30,5) + TwoTerms(6)
     tmpArray6(46,4) = Xpa*tmpArray5(31,4) + alphaXpq*TmpArray5(31,5)
     tmpArray6(47,4) = Xpa*tmpArray5(32,4) + alphaXpq*TmpArray5(32,5)
     tmpArray6(48,4) = Xpa*tmpArray5(33,4) + alphaXpq*TmpArray5(33,5)
     tmpArray6(49,4) = Xpa*tmpArray5(34,4) + alphaXpq*TmpArray5(34,5)
     tmpArray6(50,4) = Xpa*tmpArray5(35,4) + alphaXpq*TmpArray5(35,5)
     tmpArray6(51,4) = Ypa*tmpArray5(31,4) + alphaYpq*TmpArray5(31,5) + 4*TwoTerms(4)
     tmpArray6(52,4) = Zpa*tmpArray5(31,4) + alphaZpq*TmpArray5(31,5)
     tmpArray6(53,4) = Ypa*tmpArray5(33,4) + alphaYpq*TmpArray5(33,5) + 2*TwoTerms(5)
     tmpArray6(54,4) = Ypa*tmpArray5(34,4) + alphaYpq*TmpArray5(34,5) + TwoTerms(6)
     tmpArray6(55,4) = Ypa*tmpArray5(35,4) + alphaYpq*TmpArray5(35,5)
     tmpArray6(56,4) = Zpa*tmpArray5(35,4) + alphaZpq*TmpArray5(35,5) + 4*TwoTerms(6)
     TwoTerms(1) = inv2expP*(TMPAuxArray(21) + alphaP*TmpArray5(21,2))
     TwoTerms(2) = inv2expP*(TMPAuxArray(24) + alphaP*TmpArray5(24,2))
     TwoTerms(3) = inv2expP*(TMPAuxArray(26) + alphaP*TmpArray5(26,2))
     TwoTerms(4) = inv2expP*(TMPAuxArray(27) + alphaP*TmpArray5(27,2))
     TwoTerms(5) = inv2expP*(TMPAuxArray(30) + alphaP*TmpArray5(30,2))
     TwoTerms(6) = inv2expP*(TMPAuxArray(31) + alphaP*TmpArray5(31,2))
     TwoTerms(7) = inv2expP*(TMPAuxArray(33) + alphaP*TmpArray5(33,2))
     TwoTerms(8) = inv2expP*(TMPAuxArray(34) + alphaP*TmpArray5(34,2))
     TwoTerms(9) = inv2expP*(TMPAuxArray(35) + alphaP*TmpArray5(35,2))
     TMPAuxArray(57) = Xpa*TMPAuxArray(36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     TMPAuxArray(58) = Ypa*TMPAuxArray(36) + alphaYpq*TmpArray6(36,2)
     TMPAuxArray(59) = Zpa*TMPAuxArray(36) + alphaZpq*TmpArray6(36,2)
     TMPAuxArray(60) = Xpa*TMPAuxArray(39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     TMPAuxArray(61) = Ypa*TMPAuxArray(38) + alphaYpq*TmpArray6(38,2)
     TMPAuxArray(62) = Xpa*TMPAuxArray(41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     TMPAuxArray(63) = Xpa*TMPAuxArray(42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     TMPAuxArray(64) = Zpa*TMPAuxArray(39) + alphaZpq*TmpArray6(39,2)
     TMPAuxArray(65) = Ypa*TMPAuxArray(41) + alphaYpq*TmpArray6(41,2)
     TMPAuxArray(66) = Xpa*TMPAuxArray(45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     TMPAuxArray(67) = Xpa*TMPAuxArray(46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     TMPAuxArray(68) = Zpa*TMPAuxArray(42) + alphaZpq*TmpArray6(42,2)
     TMPAuxArray(69) = Xpa*TMPAuxArray(48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     TMPAuxArray(70) = Ypa*TMPAuxArray(45) + alphaYpq*TmpArray6(45,2)
     TMPAuxArray(71) = Xpa*TMPAuxArray(50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     TMPAuxArray(72) = Xpa*TMPAuxArray(51) + alphaXpq*TmpArray6(51,2)
     TMPAuxArray(73) = Xpa*TMPAuxArray(52) + alphaXpq*TmpArray6(52,2)
     TMPAuxArray(74) = Xpa*TMPAuxArray(53) + alphaXpq*TmpArray6(53,2)
     TMPAuxArray(75) = Xpa*TMPAuxArray(54) + alphaXpq*TmpArray6(54,2)
     TMPAuxArray(76) = Xpa*TMPAuxArray(55) + alphaXpq*TmpArray6(55,2)
     TMPAuxArray(77) = Xpa*TMPAuxArray(56) + alphaXpq*TmpArray6(56,2)
     TMPAuxArray(78) = Ypa*TMPAuxArray(51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     TMPAuxArray(79) = Zpa*TMPAuxArray(51) + alphaZpq*TmpArray6(51,2)
     TMPAuxArray(80) = Ypa*TMPAuxArray(53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     TMPAuxArray(81) = Ypa*TMPAuxArray(54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     TMPAuxArray(82) = Ypa*TMPAuxArray(55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     TMPAuxArray(83) = Ypa*TMPAuxArray(56) + alphaYpq*TmpArray6(56,2)
     TMPAuxArray(84) = Zpa*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expP*(TmpArray5(21,2) + alphaP*TmpArray5(21,3))
     TwoTerms(2) = inv2expP*(TmpArray5(24,2) + alphaP*TmpArray5(24,3))
     TwoTerms(3) = inv2expP*(TmpArray5(26,2) + alphaP*TmpArray5(26,3))
     TwoTerms(4) = inv2expP*(TmpArray5(27,2) + alphaP*TmpArray5(27,3))
     TwoTerms(5) = inv2expP*(TmpArray5(30,2) + alphaP*TmpArray5(30,3))
     TwoTerms(6) = inv2expP*(TmpArray5(31,2) + alphaP*TmpArray5(31,3))
     TwoTerms(7) = inv2expP*(TmpArray5(33,2) + alphaP*TmpArray5(33,3))
     TwoTerms(8) = inv2expP*(TmpArray5(34,2) + alphaP*TmpArray5(34,3))
     TwoTerms(9) = inv2expP*(TmpArray5(35,2) + alphaP*TmpArray5(35,3))
     tmpArray7(57,2) = Xpa*tmpArray6(36,2) + alphaXpq*TmpArray6(36,3) + 5*TwoTerms(1)
     tmpArray7(58,2) = Ypa*tmpArray6(36,2) + alphaYpq*TmpArray6(36,3)
     tmpArray7(59,2) = Zpa*tmpArray6(36,2) + alphaZpq*TmpArray6(36,3)
     tmpArray7(60,2) = Xpa*tmpArray6(39,2) + alphaXpq*TmpArray6(39,3) + 3*TwoTerms(2)
     tmpArray7(61,2) = Ypa*tmpArray6(38,2) + alphaYpq*TmpArray6(38,3)
     tmpArray7(62,2) = Xpa*tmpArray6(41,2) + alphaXpq*TmpArray6(41,3) + 3*TwoTerms(3)
     tmpArray7(63,2) = Xpa*tmpArray6(42,2) + alphaXpq*TmpArray6(42,3) + 2*TwoTerms(4)
     tmpArray7(64,2) = Zpa*tmpArray6(39,2) + alphaZpq*TmpArray6(39,3)
     tmpArray7(65,2) = Ypa*tmpArray6(41,2) + alphaYpq*TmpArray6(41,3)
     tmpArray7(66,2) = Xpa*tmpArray6(45,2) + alphaXpq*TmpArray6(45,3) + 2*TwoTerms(5)
     tmpArray7(67,2) = Xpa*tmpArray6(46,2) + alphaXpq*TmpArray6(46,3) + TwoTerms(6)
     tmpArray7(68,2) = Zpa*tmpArray6(42,2) + alphaZpq*TmpArray6(42,3)
     tmpArray7(69,2) = Xpa*tmpArray6(48,2) + alphaXpq*TmpArray6(48,3) + TwoTerms(7)
     tmpArray7(70,2) = Ypa*tmpArray6(45,2) + alphaYpq*TmpArray6(45,3)
     tmpArray7(71,2) = Xpa*tmpArray6(50,2) + alphaXpq*TmpArray6(50,3) + TwoTerms(9)
     tmpArray7(72,2) = Xpa*tmpArray6(51,2) + alphaXpq*TmpArray6(51,3)
     tmpArray7(73,2) = Xpa*tmpArray6(52,2) + alphaXpq*TmpArray6(52,3)
     tmpArray7(74,2) = Xpa*tmpArray6(53,2) + alphaXpq*TmpArray6(53,3)
     tmpArray7(75,2) = Xpa*tmpArray6(54,2) + alphaXpq*TmpArray6(54,3)
     tmpArray7(76,2) = Xpa*tmpArray6(55,2) + alphaXpq*TmpArray6(55,3)
     tmpArray7(77,2) = Xpa*tmpArray6(56,2) + alphaXpq*TmpArray6(56,3)
     tmpArray7(78,2) = Ypa*tmpArray6(51,2) + alphaYpq*TmpArray6(51,3) + 5*TwoTerms(6)
     tmpArray7(79,2) = Zpa*tmpArray6(51,2) + alphaZpq*TmpArray6(51,3)
     tmpArray7(80,2) = Ypa*tmpArray6(53,2) + alphaYpq*TmpArray6(53,3) + 3*TwoTerms(7)
     tmpArray7(81,2) = Ypa*tmpArray6(54,2) + alphaYpq*TmpArray6(54,3) + 2*TwoTerms(8)
     tmpArray7(82,2) = Ypa*tmpArray6(55,2) + alphaYpq*TmpArray6(55,3) + TwoTerms(9)
     tmpArray7(83,2) = Ypa*tmpArray6(56,2) + alphaYpq*TmpArray6(56,3)
     tmpArray7(84,2) = Zpa*tmpArray6(56,2) + alphaZpq*TmpArray6(56,3) + 5*TwoTerms(9)
     TwoTerms(1) = inv2expP*(TmpArray5(21,3) + alphaP*TmpArray5(21,4))
     TwoTerms(2) = inv2expP*(TmpArray5(24,3) + alphaP*TmpArray5(24,4))
     TwoTerms(3) = inv2expP*(TmpArray5(26,3) + alphaP*TmpArray5(26,4))
     TwoTerms(4) = inv2expP*(TmpArray5(27,3) + alphaP*TmpArray5(27,4))
     TwoTerms(5) = inv2expP*(TmpArray5(30,3) + alphaP*TmpArray5(30,4))
     TwoTerms(6) = inv2expP*(TmpArray5(31,3) + alphaP*TmpArray5(31,4))
     TwoTerms(7) = inv2expP*(TmpArray5(33,3) + alphaP*TmpArray5(33,4))
     TwoTerms(8) = inv2expP*(TmpArray5(34,3) + alphaP*TmpArray5(34,4))
     TwoTerms(9) = inv2expP*(TmpArray5(35,3) + alphaP*TmpArray5(35,4))
     tmpArray7(57,3) = Xpa*tmpArray6(36,3) + alphaXpq*TmpArray6(36,4) + 5*TwoTerms(1)
     tmpArray7(58,3) = Ypa*tmpArray6(36,3) + alphaYpq*TmpArray6(36,4)
     tmpArray7(59,3) = Zpa*tmpArray6(36,3) + alphaZpq*TmpArray6(36,4)
     tmpArray7(60,3) = Xpa*tmpArray6(39,3) + alphaXpq*TmpArray6(39,4) + 3*TwoTerms(2)
     tmpArray7(61,3) = Ypa*tmpArray6(38,3) + alphaYpq*TmpArray6(38,4)
     tmpArray7(62,3) = Xpa*tmpArray6(41,3) + alphaXpq*TmpArray6(41,4) + 3*TwoTerms(3)
     tmpArray7(63,3) = Xpa*tmpArray6(42,3) + alphaXpq*TmpArray6(42,4) + 2*TwoTerms(4)
     tmpArray7(64,3) = Zpa*tmpArray6(39,3) + alphaZpq*TmpArray6(39,4)
     tmpArray7(65,3) = Ypa*tmpArray6(41,3) + alphaYpq*TmpArray6(41,4)
     tmpArray7(66,3) = Xpa*tmpArray6(45,3) + alphaXpq*TmpArray6(45,4) + 2*TwoTerms(5)
     tmpArray7(67,3) = Xpa*tmpArray6(46,3) + alphaXpq*TmpArray6(46,4) + TwoTerms(6)
     tmpArray7(68,3) = Zpa*tmpArray6(42,3) + alphaZpq*TmpArray6(42,4)
     tmpArray7(69,3) = Xpa*tmpArray6(48,3) + alphaXpq*TmpArray6(48,4) + TwoTerms(7)
     tmpArray7(70,3) = Ypa*tmpArray6(45,3) + alphaYpq*TmpArray6(45,4)
     tmpArray7(71,3) = Xpa*tmpArray6(50,3) + alphaXpq*TmpArray6(50,4) + TwoTerms(9)
     tmpArray7(72,3) = Xpa*tmpArray6(51,3) + alphaXpq*TmpArray6(51,4)
     tmpArray7(73,3) = Xpa*tmpArray6(52,3) + alphaXpq*TmpArray6(52,4)
     tmpArray7(74,3) = Xpa*tmpArray6(53,3) + alphaXpq*TmpArray6(53,4)
     tmpArray7(75,3) = Xpa*tmpArray6(54,3) + alphaXpq*TmpArray6(54,4)
     tmpArray7(76,3) = Xpa*tmpArray6(55,3) + alphaXpq*TmpArray6(55,4)
     tmpArray7(77,3) = Xpa*tmpArray6(56,3) + alphaXpq*TmpArray6(56,4)
     tmpArray7(78,3) = Ypa*tmpArray6(51,3) + alphaYpq*TmpArray6(51,4) + 5*TwoTerms(6)
     tmpArray7(79,3) = Zpa*tmpArray6(51,3) + alphaZpq*TmpArray6(51,4)
     tmpArray7(80,3) = Ypa*tmpArray6(53,3) + alphaYpq*TmpArray6(53,4) + 3*TwoTerms(7)
     tmpArray7(81,3) = Ypa*tmpArray6(54,3) + alphaYpq*TmpArray6(54,4) + 2*TwoTerms(8)
     tmpArray7(82,3) = Ypa*tmpArray6(55,3) + alphaYpq*TmpArray6(55,4) + TwoTerms(9)
     tmpArray7(83,3) = Ypa*tmpArray6(56,3) + alphaYpq*TmpArray6(56,4)
     tmpArray7(84,3) = Zpa*tmpArray6(56,3) + alphaZpq*TmpArray6(56,4) + 5*TwoTerms(9)
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
     TMPAuxArray(85) = Xpa*TMPAuxArray(57) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     TMPAuxArray(86) = Ypa*TMPAuxArray(57) + alphaYpq*TmpArray7(57,2)
     TMPAuxArray(87) = Zpa*TMPAuxArray(57) + alphaZpq*TmpArray7(57,2)
     TMPAuxArray(88) = Xpa*TMPAuxArray(60) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     TMPAuxArray(89) = Ypa*TMPAuxArray(59) + alphaYpq*TmpArray7(59,2)
     TMPAuxArray(90) = Xpa*TMPAuxArray(62) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     TMPAuxArray(91) = Xpa*TMPAuxArray(63) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     TMPAuxArray(92) = Zpa*TMPAuxArray(60) + alphaZpq*TmpArray7(60,2)
     TMPAuxArray(93) = Ypa*TMPAuxArray(62) + alphaYpq*TmpArray7(62,2)
     TMPAuxArray(94) = Xpa*TMPAuxArray(66) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     TMPAuxArray(95) = Xpa*TMPAuxArray(67) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     TMPAuxArray(96) = Zpa*TMPAuxArray(63) + alphaZpq*TmpArray7(63,2)
     TMPAuxArray(97) = Xpa*TMPAuxArray(69) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     TMPAuxArray(98) = Ypa*TMPAuxArray(66) + alphaYpq*TmpArray7(66,2)
     TMPAuxArray(99) = Xpa*TMPAuxArray(71) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     TMPAuxArray(100) = Xpa*TMPAuxArray(72) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     TMPAuxArray(101) = Zpa*TMPAuxArray(67) + alphaZpq*TmpArray7(67,2)
     TMPAuxArray(102) = Xpa*TMPAuxArray(74) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     TMPAuxArray(103) = Xpa*TMPAuxArray(75) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     TMPAuxArray(104) = Ypa*TMPAuxArray(71) + alphaYpq*TmpArray7(71,2)
     TMPAuxArray(105) = Xpa*TMPAuxArray(77) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     TMPAuxArray(106) = Xpa*TMPAuxArray(78) + alphaXpq*TmpArray7(78,2)
     TMPAuxArray(107) = Xpa*TMPAuxArray(79) + alphaXpq*TmpArray7(79,2)
     TMPAuxArray(108) = Xpa*TMPAuxArray(80) + alphaXpq*TmpArray7(80,2)
     TMPAuxArray(109) = Xpa*TMPAuxArray(81) + alphaXpq*TmpArray7(81,2)
     TMPAuxArray(110) = Xpa*TMPAuxArray(82) + alphaXpq*TmpArray7(82,2)
     TMPAuxArray(111) = Xpa*TMPAuxArray(83) + alphaXpq*TmpArray7(83,2)
     TMPAuxArray(112) = Xpa*TMPAuxArray(84) + alphaXpq*TmpArray7(84,2)
     TMPAuxArray(113) = Ypa*TMPAuxArray(78) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     TMPAuxArray(114) = Zpa*TMPAuxArray(78) + alphaZpq*TmpArray7(78,2)
     TMPAuxArray(115) = Ypa*TMPAuxArray(80) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     TMPAuxArray(116) = Ypa*TMPAuxArray(81) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     TMPAuxArray(117) = Ypa*TMPAuxArray(82) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     TMPAuxArray(118) = Ypa*TMPAuxArray(83) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     TMPAuxArray(119) = Ypa*TMPAuxArray(84) + alphaYpq*TmpArray7(84,2)
     TMPAuxArray(120) = Zpa*TMPAuxArray(84) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
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
     tmpArray8(85,2) = Xpa*tmpArray7(57,2) + alphaXpq*TmpArray7(57,3) + 6*TwoTerms(1)
     tmpArray8(86,2) = Ypa*tmpArray7(57,2) + alphaYpq*TmpArray7(57,3)
     tmpArray8(87,2) = Zpa*tmpArray7(57,2) + alphaZpq*TmpArray7(57,3)
     tmpArray8(88,2) = Xpa*tmpArray7(60,2) + alphaXpq*TmpArray7(60,3) + 4*TwoTerms(2)
     tmpArray8(89,2) = Ypa*tmpArray7(59,2) + alphaYpq*TmpArray7(59,3)
     tmpArray8(90,2) = Xpa*tmpArray7(62,2) + alphaXpq*TmpArray7(62,3) + 4*TwoTerms(3)
     tmpArray8(91,2) = Xpa*tmpArray7(63,2) + alphaXpq*TmpArray7(63,3) + 3*TwoTerms(4)
     tmpArray8(92,2) = Zpa*tmpArray7(60,2) + alphaZpq*TmpArray7(60,3)
     tmpArray8(93,2) = Ypa*tmpArray7(62,2) + alphaYpq*TmpArray7(62,3)
     tmpArray8(94,2) = Xpa*tmpArray7(66,2) + alphaXpq*TmpArray7(66,3) + 3*TwoTerms(5)
     tmpArray8(95,2) = Xpa*tmpArray7(67,2) + alphaXpq*TmpArray7(67,3) + 2*TwoTerms(6)
     tmpArray8(96,2) = Zpa*tmpArray7(63,2) + alphaZpq*TmpArray7(63,3)
     tmpArray8(97,2) = Xpa*tmpArray7(69,2) + alphaXpq*TmpArray7(69,3) + 2*TwoTerms(7)
     tmpArray8(98,2) = Ypa*tmpArray7(66,2) + alphaYpq*TmpArray7(66,3)
     tmpArray8(99,2) = Xpa*tmpArray7(71,2) + alphaXpq*TmpArray7(71,3) + 2*TwoTerms(8)
     tmpArray8(100,2) = Xpa*tmpArray7(72,2) + alphaXpq*TmpArray7(72,3) + TwoTerms(9)
     tmpArray8(101,2) = Zpa*tmpArray7(67,2) + alphaZpq*TmpArray7(67,3)
     tmpArray8(102,2) = Xpa*tmpArray7(74,2) + alphaXpq*TmpArray7(74,3) + TwoTerms(10)
     tmpArray8(103,2) = Xpa*tmpArray7(75,2) + alphaXpq*TmpArray7(75,3) + TwoTerms(11)
     tmpArray8(104,2) = Ypa*tmpArray7(71,2) + alphaYpq*TmpArray7(71,3)
     tmpArray8(105,2) = Xpa*tmpArray7(77,2) + alphaXpq*TmpArray7(77,3) + TwoTerms(13)
     tmpArray8(106,2) = Xpa*tmpArray7(78,2) + alphaXpq*TmpArray7(78,3)
     tmpArray8(107,2) = Xpa*tmpArray7(79,2) + alphaXpq*TmpArray7(79,3)
     tmpArray8(108,2) = Xpa*tmpArray7(80,2) + alphaXpq*TmpArray7(80,3)
     tmpArray8(109,2) = Xpa*tmpArray7(81,2) + alphaXpq*TmpArray7(81,3)
     tmpArray8(110,2) = Xpa*tmpArray7(82,2) + alphaXpq*TmpArray7(82,3)
     tmpArray8(111,2) = Xpa*tmpArray7(83,2) + alphaXpq*TmpArray7(83,3)
     tmpArray8(112,2) = Xpa*tmpArray7(84,2) + alphaXpq*TmpArray7(84,3)
     tmpArray8(113,2) = Ypa*tmpArray7(78,2) + alphaYpq*TmpArray7(78,3) + 6*TwoTerms(9)
     tmpArray8(114,2) = Zpa*tmpArray7(78,2) + alphaZpq*TmpArray7(78,3)
     tmpArray8(115,2) = Ypa*tmpArray7(80,2) + alphaYpq*TmpArray7(80,3) + 4*TwoTerms(10)
     tmpArray8(116,2) = Ypa*tmpArray7(81,2) + alphaYpq*TmpArray7(81,3) + 3*TwoTerms(11)
     tmpArray8(117,2) = Ypa*tmpArray7(82,2) + alphaYpq*TmpArray7(82,3) + 2*TwoTerms(12)
     tmpArray8(118,2) = Ypa*tmpArray7(83,2) + alphaYpq*TmpArray7(83,3) + TwoTerms(13)
     tmpArray8(119,2) = Ypa*tmpArray7(84,2) + alphaYpq*TmpArray7(84,3)
     tmpArray8(120,2) = Zpa*tmpArray7(84,2) + alphaZpq*TmpArray7(84,3) + 6*TwoTerms(13)
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
      AuxArray(iTUV,IP) = AuxArray(iTUV,IP) + TMPAuxarray(iTUV)
     enddo
     AuxArray(121,IP) = AuxArray(121,IP) + Xpa*TMPAuxArray(85) + alphaXpq*TmpArray8(85,2) + 7*TwoTerms(1)
     AuxArray(122,IP) = AuxArray(122,IP) + Ypa*TMPAuxArray(85) + alphaYpq*TmpArray8(85,2)
     AuxArray(123,IP) = AuxArray(123,IP) + Zpa*TMPAuxArray(85) + alphaZpq*TmpArray8(85,2)
     AuxArray(124,IP) = AuxArray(124,IP) + Xpa*TMPAuxArray(88) + alphaXpq*TmpArray8(88,2) + 5*TwoTerms(2)
     AuxArray(125,IP) = AuxArray(125,IP) + Ypa*TMPAuxArray(87) + alphaYpq*TmpArray8(87,2)
     AuxArray(126,IP) = AuxArray(126,IP) + Xpa*TMPAuxArray(90) + alphaXpq*TmpArray8(90,2) + 5*TwoTerms(3)
     AuxArray(127,IP) = AuxArray(127,IP) + Xpa*TMPAuxArray(91) + alphaXpq*TmpArray8(91,2) + 4*TwoTerms(4)
     AuxArray(128,IP) = AuxArray(128,IP) + Zpa*TMPAuxArray(88) + alphaZpq*TmpArray8(88,2)
     AuxArray(129,IP) = AuxArray(129,IP) + Ypa*TMPAuxArray(90) + alphaYpq*TmpArray8(90,2)
     AuxArray(130,IP) = AuxArray(130,IP) + Xpa*TMPAuxArray(94) + alphaXpq*TmpArray8(94,2) + 4*TwoTerms(5)
     AuxArray(131,IP) = AuxArray(131,IP) + Xpa*TMPAuxArray(95) + alphaXpq*TmpArray8(95,2) + 3*TwoTerms(6)
     AuxArray(132,IP) = AuxArray(132,IP) + Zpa*TMPAuxArray(91) + alphaZpq*TmpArray8(91,2)
     AuxArray(133,IP) = AuxArray(133,IP) + Xpa*TMPAuxArray(97) + alphaXpq*TmpArray8(97,2) + 3*TwoTerms(7)
     AuxArray(134,IP) = AuxArray(134,IP) + Ypa*TMPAuxArray(94) + alphaYpq*TmpArray8(94,2)
     AuxArray(135,IP) = AuxArray(135,IP) + Xpa*TMPAuxArray(99) + alphaXpq*TmpArray8(99,2) + 3*TwoTerms(8)
     AuxArray(136,IP) = AuxArray(136,IP) + Xpa*TMPAuxArray(100) + alphaXpq*TmpArray8(100,2) + 2*TwoTerms(9)
     AuxArray(137,IP) = AuxArray(137,IP) + Zpa*TMPAuxArray(95) + alphaZpq*TmpArray8(95,2)
     AuxArray(138,IP) = AuxArray(138,IP) + Xpa*TMPAuxArray(102) + alphaXpq*TmpArray8(102,2) + 2*TwoTerms(10)
     AuxArray(139,IP) = AuxArray(139,IP) + Xpa*TMPAuxArray(103) + alphaXpq*TmpArray8(103,2) + 2*TwoTerms(11)
     AuxArray(140,IP) = AuxArray(140,IP) + Ypa*TMPAuxArray(99) + alphaYpq*TmpArray8(99,2)
     AuxArray(141,IP) = AuxArray(141,IP) + Xpa*TMPAuxArray(105) + alphaXpq*TmpArray8(105,2) + 2*TwoTerms(12)
     AuxArray(142,IP) = AuxArray(142,IP) + Xpa*TMPAuxArray(106) + alphaXpq*TmpArray8(106,2) + TwoTerms(13)
     AuxArray(143,IP) = AuxArray(143,IP) + Zpa*TMPAuxArray(100) + alphaZpq*TmpArray8(100,2)
     AuxArray(144,IP) = AuxArray(144,IP) + Xpa*TMPAuxArray(108) + alphaXpq*TmpArray8(108,2) + TwoTerms(14)
     AuxArray(145,IP) = AuxArray(145,IP) + Xpa*TMPAuxArray(109) + alphaXpq*TmpArray8(109,2) + TwoTerms(15)
     AuxArray(146,IP) = AuxArray(146,IP) + Xpa*TMPAuxArray(110) + alphaXpq*TmpArray8(110,2) + TwoTerms(16)
     AuxArray(147,IP) = AuxArray(147,IP) + Ypa*TMPAuxArray(105) + alphaYpq*TmpArray8(105,2)
     AuxArray(148,IP) = AuxArray(148,IP) + Xpa*TMPAuxArray(112) + alphaXpq*TmpArray8(112,2) + TwoTerms(18)
     AuxArray(149,IP) = AuxArray(149,IP) + Xpa*TMPAuxArray(113) + alphaXpq*TmpArray8(113,2)
     AuxArray(150,IP) = AuxArray(150,IP) + Xpa*TMPAuxArray(114) + alphaXpq*TmpArray8(114,2)
     AuxArray(151,IP) = AuxArray(151,IP) + Xpa*TMPAuxArray(115) + alphaXpq*TmpArray8(115,2)
     AuxArray(152,IP) = AuxArray(152,IP) + Xpa*TMPAuxArray(116) + alphaXpq*TmpArray8(116,2)
     AuxArray(153,IP) = AuxArray(153,IP) + Xpa*TMPAuxArray(117) + alphaXpq*TmpArray8(117,2)
     AuxArray(154,IP) = AuxArray(154,IP) + Xpa*TMPAuxArray(118) + alphaXpq*TmpArray8(118,2)
     AuxArray(155,IP) = AuxArray(155,IP) + Xpa*TMPAuxArray(119) + alphaXpq*TmpArray8(119,2)
     AuxArray(156,IP) = AuxArray(156,IP) + Xpa*TMPAuxArray(120) + alphaXpq*TmpArray8(120,2)
     AuxArray(157,IP) = AuxArray(157,IP) + Ypa*TMPAuxArray(113) + alphaYpq*TmpArray8(113,2) + 7*TwoTerms(13)
     AuxArray(158,IP) = AuxArray(158,IP) + Zpa*TMPAuxArray(113) + alphaZpq*TmpArray8(113,2)
     AuxArray(159,IP) = AuxArray(159,IP) + Ypa*TMPAuxArray(115) + alphaYpq*TmpArray8(115,2) + 5*TwoTerms(14)
     AuxArray(160,IP) = AuxArray(160,IP) + Ypa*TMPAuxArray(116) + alphaYpq*TmpArray8(116,2) + 4*TwoTerms(15)
     AuxArray(161,IP) = AuxArray(161,IP) + Ypa*TMPAuxArray(117) + alphaYpq*TmpArray8(117,2) + 3*TwoTerms(16)
     AuxArray(162,IP) = AuxArray(162,IP) + Ypa*TMPAuxArray(118) + alphaYpq*TmpArray8(118,2) + 2*TwoTerms(17)
     AuxArray(163,IP) = AuxArray(163,IP) + Ypa*TMPAuxArray(119) + alphaYpq*TmpArray8(119,2) + TwoTerms(18)
     AuxArray(164,IP) = AuxArray(164,IP) + Ypa*TMPAuxArray(120) + alphaYpq*TmpArray8(120,2)
     AuxArray(165,IP) = AuxArray(165,IP) + Zpa*TMPAuxArray(120) + alphaZpq*TmpArray8(120,2) + 7*TwoTerms(18)
    ENDDO
   ENDDO
  ENDDO
 end subroutine
end module
