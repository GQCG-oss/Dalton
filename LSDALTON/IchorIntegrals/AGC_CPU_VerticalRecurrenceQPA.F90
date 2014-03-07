MODULE AGC_CPU_OBS_VERTICALRECURRENCEMODA
 use IchorPrecisionModule
  
 CONTAINS

subroutine VerticalRecurrenceCPU0(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,&
         & Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  REAL(REALK),intent(in) :: integralPrefactor(nprimQ,nPrimP)
  REAL(REALK),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  REAL(REALK),intent(in) :: TABFJW(0:3,0:1200)
  REAL(REALK),intent(in) :: QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: AUXarray(nPrimQ,nPrimP,nPasses)
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
  Integer :: IPNT,iPassP,iPrimP,iPrimQ,iPQ,iAtomA,iAtomB
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,px,py,pz,pqx,pqy,pqz,&
!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$OMP         Pexpfac,iPrimP,iPrimQ,&
!$OMP         RWVAL,GVAL,AUXarray) 
!!$OMP SHARED(nPasses,iAtomApass,iAtomBpass,PpreExpFac,&
!!$OMP        nPrimP,nPrimQ,&
!!$OMP        QpreExpFac,Pcent,Qcent,reducedExponents,TABFJW,&
!!$OMP        integralPrefactor)
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    px = Pcent(1,iPrimP,iAtomA,iAtomB)
    py = Pcent(2,iPrimP,iAtomA,iAtomB)
    pz = Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     pqx = px - Qcent(1,iPrimQ)
     pqy = py - Qcent(2,iPrimQ)
     pqz = pz - Qcent(3,iPrimQ)
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
     AUXarray(iPrimQ,iPrimP,iPassP)=integralPrefactor(iPrimQ,iPrimP)*&
          & QpreExpFac(iPrimQ)*Pexpfac*RJ000
    enddo
   enddo
  enddo
!$OMP END DO
end subroutine VerticalRecurrenceCPU0

subroutine VerticalRecurrenceCPU1A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: TABFJW(0:4,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA)
  real(realk),intent(inout) :: AUXarray(4,nPrimQ,nPrimP,nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,iAtomA,iAtomB
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
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,invexpP,&
!$OMP         Pexpfac,mPX,mPY,mPZ,RWVAL,GVAL,&
!$OMP         alphaXpq,alphaYpq,alphaZpq) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
    DO iPrimQ=1, nPrimQ
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
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     TMP1 = PREF*RJ000(0)
     TMP2 = PREF*RJ000(1)
     AUXarray(1,iPrimQ,iPrimP,iPassP) = TMP1
     AUXarray(2,iPrimQ,iPrimP,iPassP) = Xpa*TMP1 + alphaXpq*TMP2
     AUXarray(3,iPrimQ,iPrimP,iPassP) = Ypa*TMP1 + alphaYpq*TMP2
     AUXarray(4,iPrimQ,iPrimP,iPassP) = Zpa*TMP1 + alphaZpq*TMP2
    enddo
   enddo
  enddo
!$OMP END DO
end subroutine

subroutine BuildRJ000CPU2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: TABFJW(0: 5,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(inout) :: RJ000array(0: 5,nPrimQ,nPrimP,nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,iAtomA,iAtomB
  real(realk) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(realk) :: RJ000(0: 5)
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
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         mPX,mPY,mPZ,RWVAL,GVAL)
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   DO iPrimP=1, nPrimP
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
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
      RJ000Array( 0,iPrimQ,iPrimP,iPassP) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000( 1) = RWVAL*(( 1 - D05)*RJ000( 0)-REXPW)
      RJ000( 2) = RWVAL*(( 2 - D05)*RJ000( 1)-REXPW)
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
     ENDIF
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine BuildRJ000CPU3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: TABFJW(0: 6,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(inout) :: RJ000array(0: 6,nPrimQ,nPrimP,nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,iAtomA,iAtomB
  real(realk) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(realk) :: RJ000(0: 6)
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
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         mPX,mPY,mPZ,RWVAL,GVAL)
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   DO iPrimP=1, nPrimP
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
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
      RJ000Array( 0,iPrimQ,iPrimP,iPassP) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000( 1) = RWVAL*(( 1 - D05)*RJ000( 0)-REXPW)
      RJ000( 2) = RWVAL*(( 2 - D05)*RJ000( 1)-REXPW)
      RJ000( 3) = RWVAL*(( 3 - D05)*RJ000( 2)-REXPW)
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000( 3) = RWVAL*( 3 - D05)*RJ000( 2)
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
     ENDIF
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine BuildRJ000CPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: TABFJW(0: 7,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(inout) :: RJ000array(0: 7,nPrimQ,nPrimP,nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,iAtomA,iAtomB
  real(realk) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(realk) :: RJ000(0: 7)
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
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         mPX,mPY,mPZ,RWVAL,GVAL)
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   DO iPrimP=1, nPrimP
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
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
      RJ000Array( 0,iPrimQ,iPrimP,iPassP) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
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
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = RJ000( 4)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000( 3) = RWVAL*( 3 - D05)*RJ000( 2)
      RJ000( 4) = RWVAL*( 4 - D05)*RJ000( 3)
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = RJ000( 4)
     ENDIF
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine BuildRJ000CPU5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: TABFJW(0: 8,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(inout) :: RJ000array(0: 8,nPrimQ,nPrimP,nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,iAtomA,iAtomB
  real(realk) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(realk) :: RJ000(0: 8)
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
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         mPX,mPY,mPZ,RWVAL,GVAL)
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   DO iPrimP=1, nPrimP
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
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
      RJ000Array( 0,iPrimQ,iPrimP,iPassP) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
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
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = RJ000( 4)
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = RJ000( 5)
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
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = RJ000( 4)
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = RJ000( 5)
     ENDIF
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine BuildRJ000CPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: TABFJW(0: 9,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(inout) :: RJ000array(0: 9,nPrimQ,nPrimP,nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,iAtomA,iAtomB
  real(realk) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(realk) :: RJ000(0: 9)
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
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         mPX,mPY,mPZ,RWVAL,GVAL)
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   DO iPrimP=1, nPrimP
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
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
      RJ000Array( 0,iPrimQ,iPrimP,iPassP) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
      RJ000Array( 6,iPrimQ,iPrimP,iPassP) = TABFJW( 6,IPNT)-TABFJW( 7,IPNT)*WDIFF+TABFJW( 8,IPNT)*W2+TABFJW( 9,IPNT)*W3
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
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = RJ000( 4)
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = RJ000( 5)
      RJ000Array( 6,iPrimQ,iPrimP,iPassP) = RJ000( 6)
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
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = RJ000( 4)
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = RJ000( 5)
      RJ000Array( 6,iPrimQ,iPrimP,iPassP) = RJ000( 6)
     ENDIF
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine BuildRJ000CPU7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: TABFJW(0:10,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(inout) :: RJ000array(0:10,nPrimQ,nPrimP,nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,iAtomA,iAtomB
  real(realk) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(realk) :: RJ000(0:10)
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
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         mPX,mPY,mPZ,RWVAL,GVAL)
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   DO iPrimP=1, nPrimP
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
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
      RJ000Array( 0,iPrimQ,iPrimP,iPassP) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
      RJ000Array( 6,iPrimQ,iPrimP,iPassP) = TABFJW( 6,IPNT)-TABFJW( 7,IPNT)*WDIFF+TABFJW( 8,IPNT)*W2+TABFJW( 9,IPNT)*W3
      RJ000Array( 7,iPrimQ,iPrimP,iPassP) = TABFJW( 7,IPNT)-TABFJW( 8,IPNT)*WDIFF+TABFJW( 9,IPNT)*W2+TABFJW(10,IPNT)*W3
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
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = RJ000( 4)
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = RJ000( 5)
      RJ000Array( 6,iPrimQ,iPrimP,iPassP) = RJ000( 6)
      RJ000Array( 7,iPrimQ,iPrimP,iPassP) = RJ000( 7)
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
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = RJ000( 4)
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = RJ000( 5)
      RJ000Array( 6,iPrimQ,iPrimP,iPassP) = RJ000( 6)
      RJ000Array( 7,iPrimQ,iPrimP,iPassP) = RJ000( 7)
     ENDIF
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine BuildRJ000CPU8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: TABFJW(0:11,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(inout) :: RJ000array(0:11,nPrimQ,nPrimP,nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,iAtomA,iAtomB
  real(realk) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(realk) :: RJ000(0:11)
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
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         mPX,mPY,mPZ,RWVAL,GVAL)
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   DO iPrimP=1, nPrimP
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
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
      RJ000Array( 0,iPrimQ,iPrimP,iPassP) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
      RJ000Array( 6,iPrimQ,iPrimP,iPassP) = TABFJW( 6,IPNT)-TABFJW( 7,IPNT)*WDIFF+TABFJW( 8,IPNT)*W2+TABFJW( 9,IPNT)*W3
      RJ000Array( 7,iPrimQ,iPrimP,iPassP) = TABFJW( 7,IPNT)-TABFJW( 8,IPNT)*WDIFF+TABFJW( 9,IPNT)*W2+TABFJW(10,IPNT)*W3
      RJ000Array( 8,iPrimQ,iPrimP,iPassP) = TABFJW( 8,IPNT)-TABFJW( 9,IPNT)*WDIFF+TABFJW(10,IPNT)*W2+TABFJW(11,IPNT)*W3
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
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = RJ000( 4)
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = RJ000( 5)
      RJ000Array( 6,iPrimQ,iPrimP,iPassP) = RJ000( 6)
      RJ000Array( 7,iPrimQ,iPrimP,iPassP) = RJ000( 7)
      RJ000Array( 8,iPrimQ,iPrimP,iPassP) = RJ000( 8)
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
      RJ000Array(0,iPrimQ,iPrimP,iPassP) = RJ000(0)
      RJ000Array( 1,iPrimQ,iPrimP,iPassP) = RJ000( 1)
      RJ000Array( 2,iPrimQ,iPrimP,iPassP) = RJ000( 2)
      RJ000Array( 3,iPrimQ,iPrimP,iPassP) = RJ000( 3)
      RJ000Array( 4,iPrimQ,iPrimP,iPassP) = RJ000( 4)
      RJ000Array( 5,iPrimQ,iPrimP,iPassP) = RJ000( 5)
      RJ000Array( 6,iPrimQ,iPrimP,iPassP) = RJ000( 6)
      RJ000Array( 7,iPrimQ,iPrimP,iPassP) = RJ000( 7)
      RJ000Array( 8,iPrimQ,iPrimP,iPassP) = RJ000( 8)
     ENDIF
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU2A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & RJ000Array,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 5,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA)
  real(realk),intent(inout) :: AUXarray(   10,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(   1)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:3)
  real(realk) :: TMParray2(  2:  4,2:2)
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,invexpP,inv2expP,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
   iP = (iPassP-1)*nPrimQ*nPrimP
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xpa*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Ypa*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zpa*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpa*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypa*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpa*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,IP) = Xpa*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xpa*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xpa*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Ypa*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Ypa*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zpa*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU3A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & RJ000Array,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 6,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA)
  real(realk),intent(inout) :: AUXarray(   20,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(   3)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:4)
  real(realk) :: TMParray2(  2:  4,2:3)
  real(realk) :: TMParray3(  5: 10,2:2)
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,invexpP,inv2expP,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
   iP = (iPassP-1)*nPrimQ*nPrimP
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xpa*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Ypa*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zpa*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpa*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypa*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpa*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpa*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypa*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpa*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,IP) = Xpa*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xpa*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xpa*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Ypa*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Ypa*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zpa*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(TmpArray1(1,2) + alphaP*TmpArray1(1,3))
     tmpArray3(5,2) = Xpa*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xpa*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xpa*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Ypa*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Ypa*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zpa*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expP*(AuxArray(2,IP) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(3,IP) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(4,IP) + alphaP*TmpArray2(4,2))
     AuxArray(11,IP) = Xpa*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Ypa*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zpa*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xpa*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xpa*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xpa*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Ypa*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zpa*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Ypa*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zpa*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU4A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & RJ000Array,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 7,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA)
  real(realk),intent(inout) :: AUXarray(   35,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(   6)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:5)
  real(realk) :: TMParray2(  2:  4,2:4)
  real(realk) :: TMParray3(  5: 10,2:3)
  real(realk) :: TMParray4( 11: 20,2:2)
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,invexpP,inv2expP,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
   iP = (iPassP-1)*nPrimQ*nPrimP
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xpa*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Ypa*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zpa*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xpa*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Ypa*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zpa*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xpa*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Ypa*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zpa*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xpa*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Ypa*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zpa*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,IP) = Xpa*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xpa*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xpa*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Ypa*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Ypa*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zpa*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(AuxArray(2,IP) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(3,IP) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(4,IP) + alphaP*TmpArray2(4,2))
     AuxArray(11,IP) = Xpa*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Ypa*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zpa*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xpa*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xpa*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xpa*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Ypa*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zpa*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Ypa*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zpa*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(5,IP) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(AuxArray(8,IP) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(AuxArray(10,IP) + alphaP*TmpArray3(10,2))
     AuxArray(21,IP) = Xpa*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Ypa*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zpa*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xpa*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Ypa*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xpa*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xpa*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xpa*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xpa*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xpa*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Ypa*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zpa*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Ypa*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Ypa*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zpa*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU5A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & RJ000Array,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 8,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA)
  real(realk),intent(inout) :: AUXarray(   56,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(  10)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:6)
  real(realk) :: TMParray2(  2:  4,2:5)
  real(realk) :: TMParray3(  5: 10,2:4)
  real(realk) :: TMParray4( 11: 20,2:3)
  real(realk) :: TMParray5( 21: 35,2:2)
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,invexpP,inv2expP,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
   iP = (iPassP-1)*nPrimQ*nPrimP
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 6) = PREF*RJ000Array( 5,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xpa*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Ypa*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zpa*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,IP) = Xpa*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xpa*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xpa*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Ypa*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Ypa*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zpa*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(AuxArray(2,IP) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(3,IP) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(4,IP) + alphaP*TmpArray2(4,2))
     AuxArray(11,IP) = Xpa*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Ypa*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zpa*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xpa*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xpa*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xpa*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Ypa*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zpa*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Ypa*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zpa*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(5,IP) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(AuxArray(8,IP) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(AuxArray(10,IP) + alphaP*TmpArray3(10,2))
     AuxArray(21,IP) = Xpa*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Ypa*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zpa*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xpa*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Ypa*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xpa*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xpa*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xpa*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xpa*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xpa*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Ypa*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zpa*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Ypa*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Ypa*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zpa*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(11,IP) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(AuxArray(14,IP) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(AuxArray(16,IP) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(AuxArray(17,IP) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(AuxArray(19,IP) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(AuxArray(20,IP) + alphaP*TmpArray4(20,2))
     AuxArray(36,IP) = Xpa*AuxArray(21,IP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IP) = Ypa*AuxArray(21,IP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IP) = Zpa*AuxArray(21,IP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IP) = Xpa*AuxArray(24,IP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IP) = Ypa*AuxArray(23,IP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IP) = Xpa*AuxArray(26,IP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IP) = Xpa*AuxArray(27,IP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IP) = Zpa*AuxArray(24,IP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IP) = Ypa*AuxArray(26,IP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IP) = Xpa*AuxArray(30,IP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IP) = Xpa*AuxArray(31,IP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IP) = Xpa*AuxArray(32,IP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IP) = Xpa*AuxArray(33,IP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IP) = Xpa*AuxArray(34,IP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IP) = Xpa*AuxArray(35,IP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IP) = Ypa*AuxArray(31,IP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IP) = Zpa*AuxArray(31,IP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IP) = Ypa*AuxArray(33,IP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IP) = Ypa*AuxArray(34,IP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IP) = Ypa*AuxArray(35,IP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IP) = Zpa*AuxArray(35,IP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU6A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & RJ000Array,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 9,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA)
  real(realk),intent(inout) :: AUXarray(   84,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(  15)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:7)
  real(realk) :: TMParray2(  2:  4,2:6)
  real(realk) :: TMParray3(  5: 10,2:5)
  real(realk) :: TMParray4( 11: 20,2:4)
  real(realk) :: TMParray5( 21: 35,2:3)
  real(realk) :: TMParray6( 36: 56,2:2)
  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,invexpP,inv2expP,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
   iP = (iPassP-1)*nPrimQ*nPrimP
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 6) = PREF*RJ000Array( 5,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 7) = PREF*RJ000Array( 6,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xpa*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Ypa*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zpa*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,IP) = Xpa*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xpa*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xpa*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Ypa*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Ypa*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zpa*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(AuxArray(2,IP) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(3,IP) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(4,IP) + alphaP*TmpArray2(4,2))
     AuxArray(11,IP) = Xpa*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Ypa*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zpa*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xpa*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xpa*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xpa*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Ypa*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zpa*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Ypa*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zpa*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(5,IP) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(AuxArray(8,IP) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(AuxArray(10,IP) + alphaP*TmpArray3(10,2))
     AuxArray(21,IP) = Xpa*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Ypa*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zpa*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xpa*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Ypa*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xpa*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xpa*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xpa*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xpa*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xpa*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Ypa*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zpa*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Ypa*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Ypa*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zpa*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(11,IP) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(AuxArray(14,IP) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(AuxArray(16,IP) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(AuxArray(17,IP) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(AuxArray(19,IP) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(AuxArray(20,IP) + alphaP*TmpArray4(20,2))
     AuxArray(36,IP) = Xpa*AuxArray(21,IP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IP) = Ypa*AuxArray(21,IP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IP) = Zpa*AuxArray(21,IP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IP) = Xpa*AuxArray(24,IP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IP) = Ypa*AuxArray(23,IP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IP) = Xpa*AuxArray(26,IP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IP) = Xpa*AuxArray(27,IP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IP) = Zpa*AuxArray(24,IP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IP) = Ypa*AuxArray(26,IP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IP) = Xpa*AuxArray(30,IP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IP) = Xpa*AuxArray(31,IP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IP) = Xpa*AuxArray(32,IP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IP) = Xpa*AuxArray(33,IP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IP) = Xpa*AuxArray(34,IP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IP) = Xpa*AuxArray(35,IP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IP) = Ypa*AuxArray(31,IP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IP) = Zpa*AuxArray(31,IP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IP) = Ypa*AuxArray(33,IP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IP) = Ypa*AuxArray(34,IP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IP) = Ypa*AuxArray(35,IP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IP) = Zpa*AuxArray(35,IP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     TwoTerms(1) = inv2expP*(AuxArray(21,IP) + alphaP*TmpArray5(21,2))
     TwoTerms(2) = inv2expP*(AuxArray(24,IP) + alphaP*TmpArray5(24,2))
     TwoTerms(3) = inv2expP*(AuxArray(26,IP) + alphaP*TmpArray5(26,2))
     TwoTerms(4) = inv2expP*(AuxArray(27,IP) + alphaP*TmpArray5(27,2))
     TwoTerms(5) = inv2expP*(AuxArray(30,IP) + alphaP*TmpArray5(30,2))
     TwoTerms(6) = inv2expP*(AuxArray(31,IP) + alphaP*TmpArray5(31,2))
     TwoTerms(7) = inv2expP*(AuxArray(33,IP) + alphaP*TmpArray5(33,2))
     TwoTerms(8) = inv2expP*(AuxArray(34,IP) + alphaP*TmpArray5(34,2))
     TwoTerms(9) = inv2expP*(AuxArray(35,IP) + alphaP*TmpArray5(35,2))
     AuxArray(57,IP) = Xpa*AuxArray(36,IP) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,IP) = Ypa*AuxArray(36,IP) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,IP) = Zpa*AuxArray(36,IP) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,IP) = Xpa*AuxArray(39,IP) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,IP) = Ypa*AuxArray(38,IP) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,IP) = Xpa*AuxArray(41,IP) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,IP) = Xpa*AuxArray(42,IP) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,IP) = Zpa*AuxArray(39,IP) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,IP) = Ypa*AuxArray(41,IP) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,IP) = Xpa*AuxArray(45,IP) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,IP) = Xpa*AuxArray(46,IP) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,IP) = Zpa*AuxArray(42,IP) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,IP) = Xpa*AuxArray(48,IP) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,IP) = Ypa*AuxArray(45,IP) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,IP) = Xpa*AuxArray(50,IP) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,IP) = Xpa*AuxArray(51,IP) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,IP) = Xpa*AuxArray(52,IP) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,IP) = Xpa*AuxArray(53,IP) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,IP) = Xpa*AuxArray(54,IP) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,IP) = Xpa*AuxArray(55,IP) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,IP) = Xpa*AuxArray(56,IP) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,IP) = Ypa*AuxArray(51,IP) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,IP) = Zpa*AuxArray(51,IP) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,IP) = Ypa*AuxArray(53,IP) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,IP) = Ypa*AuxArray(54,IP) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,IP) = Ypa*AuxArray(55,IP) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,IP) = Ypa*AuxArray(56,IP) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,IP) = Zpa*AuxArray(56,IP) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU7A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & RJ000Array,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0:10,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA)
  real(realk),intent(inout) :: AUXarray(  120,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(  21)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
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
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,invexpP,inv2expP,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
   iP = (iPassP-1)*nPrimQ*nPrimP
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 6) = PREF*RJ000Array( 5,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 7) = PREF*RJ000Array( 6,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 8) = PREF*RJ000Array( 7,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xpa*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Ypa*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zpa*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,IP) = Xpa*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xpa*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xpa*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Ypa*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Ypa*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zpa*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(AuxArray(2,IP) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(3,IP) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(4,IP) + alphaP*TmpArray2(4,2))
     AuxArray(11,IP) = Xpa*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Ypa*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zpa*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xpa*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xpa*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xpa*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Ypa*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zpa*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Ypa*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zpa*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(5,IP) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(AuxArray(8,IP) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(AuxArray(10,IP) + alphaP*TmpArray3(10,2))
     AuxArray(21,IP) = Xpa*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Ypa*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zpa*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xpa*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Ypa*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xpa*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xpa*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xpa*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xpa*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xpa*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Ypa*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zpa*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Ypa*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Ypa*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zpa*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(11,IP) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(AuxArray(14,IP) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(AuxArray(16,IP) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(AuxArray(17,IP) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(AuxArray(19,IP) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(AuxArray(20,IP) + alphaP*TmpArray4(20,2))
     AuxArray(36,IP) = Xpa*AuxArray(21,IP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IP) = Ypa*AuxArray(21,IP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IP) = Zpa*AuxArray(21,IP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IP) = Xpa*AuxArray(24,IP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IP) = Ypa*AuxArray(23,IP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IP) = Xpa*AuxArray(26,IP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IP) = Xpa*AuxArray(27,IP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IP) = Zpa*AuxArray(24,IP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IP) = Ypa*AuxArray(26,IP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IP) = Xpa*AuxArray(30,IP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IP) = Xpa*AuxArray(31,IP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IP) = Xpa*AuxArray(32,IP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IP) = Xpa*AuxArray(33,IP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IP) = Xpa*AuxArray(34,IP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IP) = Xpa*AuxArray(35,IP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IP) = Ypa*AuxArray(31,IP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IP) = Zpa*AuxArray(31,IP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IP) = Ypa*AuxArray(33,IP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IP) = Ypa*AuxArray(34,IP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IP) = Ypa*AuxArray(35,IP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IP) = Zpa*AuxArray(35,IP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     TwoTerms(1) = inv2expP*(AuxArray(21,IP) + alphaP*TmpArray5(21,2))
     TwoTerms(2) = inv2expP*(AuxArray(24,IP) + alphaP*TmpArray5(24,2))
     TwoTerms(3) = inv2expP*(AuxArray(26,IP) + alphaP*TmpArray5(26,2))
     TwoTerms(4) = inv2expP*(AuxArray(27,IP) + alphaP*TmpArray5(27,2))
     TwoTerms(5) = inv2expP*(AuxArray(30,IP) + alphaP*TmpArray5(30,2))
     TwoTerms(6) = inv2expP*(AuxArray(31,IP) + alphaP*TmpArray5(31,2))
     TwoTerms(7) = inv2expP*(AuxArray(33,IP) + alphaP*TmpArray5(33,2))
     TwoTerms(8) = inv2expP*(AuxArray(34,IP) + alphaP*TmpArray5(34,2))
     TwoTerms(9) = inv2expP*(AuxArray(35,IP) + alphaP*TmpArray5(35,2))
     AuxArray(57,IP) = Xpa*AuxArray(36,IP) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,IP) = Ypa*AuxArray(36,IP) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,IP) = Zpa*AuxArray(36,IP) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,IP) = Xpa*AuxArray(39,IP) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,IP) = Ypa*AuxArray(38,IP) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,IP) = Xpa*AuxArray(41,IP) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,IP) = Xpa*AuxArray(42,IP) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,IP) = Zpa*AuxArray(39,IP) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,IP) = Ypa*AuxArray(41,IP) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,IP) = Xpa*AuxArray(45,IP) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,IP) = Xpa*AuxArray(46,IP) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,IP) = Zpa*AuxArray(42,IP) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,IP) = Xpa*AuxArray(48,IP) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,IP) = Ypa*AuxArray(45,IP) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,IP) = Xpa*AuxArray(50,IP) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,IP) = Xpa*AuxArray(51,IP) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,IP) = Xpa*AuxArray(52,IP) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,IP) = Xpa*AuxArray(53,IP) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,IP) = Xpa*AuxArray(54,IP) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,IP) = Xpa*AuxArray(55,IP) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,IP) = Xpa*AuxArray(56,IP) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,IP) = Ypa*AuxArray(51,IP) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,IP) = Zpa*AuxArray(51,IP) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,IP) = Ypa*AuxArray(53,IP) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,IP) = Ypa*AuxArray(54,IP) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,IP) = Ypa*AuxArray(55,IP) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,IP) = Ypa*AuxArray(56,IP) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,IP) = Zpa*AuxArray(56,IP) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
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
     TwoTerms(1) = inv2expP*(AuxArray(36,IP) + alphaP*TmpArray6(36,2))
     TwoTerms(2) = inv2expP*(AuxArray(39,IP) + alphaP*TmpArray6(39,2))
     TwoTerms(3) = inv2expP*(AuxArray(41,IP) + alphaP*TmpArray6(41,2))
     TwoTerms(4) = inv2expP*(AuxArray(42,IP) + alphaP*TmpArray6(42,2))
     TwoTerms(5) = inv2expP*(AuxArray(45,IP) + alphaP*TmpArray6(45,2))
     TwoTerms(6) = inv2expP*(AuxArray(46,IP) + alphaP*TmpArray6(46,2))
     TwoTerms(7) = inv2expP*(AuxArray(48,IP) + alphaP*TmpArray6(48,2))
     TwoTerms(8) = inv2expP*(AuxArray(50,IP) + alphaP*TmpArray6(50,2))
     TwoTerms(9) = inv2expP*(AuxArray(51,IP) + alphaP*TmpArray6(51,2))
     TwoTerms(10) = inv2expP*(AuxArray(53,IP) + alphaP*TmpArray6(53,2))
     TwoTerms(11) = inv2expP*(AuxArray(54,IP) + alphaP*TmpArray6(54,2))
     TwoTerms(12) = inv2expP*(AuxArray(55,IP) + alphaP*TmpArray6(55,2))
     TwoTerms(13) = inv2expP*(AuxArray(56,IP) + alphaP*TmpArray6(56,2))
     AuxArray(85,IP) = Xpa*AuxArray(57,IP) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     AuxArray(86,IP) = Ypa*AuxArray(57,IP) + alphaYpq*TmpArray7(57,2)
     AuxArray(87,IP) = Zpa*AuxArray(57,IP) + alphaZpq*TmpArray7(57,2)
     AuxArray(88,IP) = Xpa*AuxArray(60,IP) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     AuxArray(89,IP) = Ypa*AuxArray(59,IP) + alphaYpq*TmpArray7(59,2)
     AuxArray(90,IP) = Xpa*AuxArray(62,IP) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     AuxArray(91,IP) = Xpa*AuxArray(63,IP) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     AuxArray(92,IP) = Zpa*AuxArray(60,IP) + alphaZpq*TmpArray7(60,2)
     AuxArray(93,IP) = Ypa*AuxArray(62,IP) + alphaYpq*TmpArray7(62,2)
     AuxArray(94,IP) = Xpa*AuxArray(66,IP) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     AuxArray(95,IP) = Xpa*AuxArray(67,IP) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     AuxArray(96,IP) = Zpa*AuxArray(63,IP) + alphaZpq*TmpArray7(63,2)
     AuxArray(97,IP) = Xpa*AuxArray(69,IP) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     AuxArray(98,IP) = Ypa*AuxArray(66,IP) + alphaYpq*TmpArray7(66,2)
     AuxArray(99,IP) = Xpa*AuxArray(71,IP) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     AuxArray(100,IP) = Xpa*AuxArray(72,IP) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     AuxArray(101,IP) = Zpa*AuxArray(67,IP) + alphaZpq*TmpArray7(67,2)
     AuxArray(102,IP) = Xpa*AuxArray(74,IP) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     AuxArray(103,IP) = Xpa*AuxArray(75,IP) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     AuxArray(104,IP) = Ypa*AuxArray(71,IP) + alphaYpq*TmpArray7(71,2)
     AuxArray(105,IP) = Xpa*AuxArray(77,IP) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     AuxArray(106,IP) = Xpa*AuxArray(78,IP) + alphaXpq*TmpArray7(78,2)
     AuxArray(107,IP) = Xpa*AuxArray(79,IP) + alphaXpq*TmpArray7(79,2)
     AuxArray(108,IP) = Xpa*AuxArray(80,IP) + alphaXpq*TmpArray7(80,2)
     AuxArray(109,IP) = Xpa*AuxArray(81,IP) + alphaXpq*TmpArray7(81,2)
     AuxArray(110,IP) = Xpa*AuxArray(82,IP) + alphaXpq*TmpArray7(82,2)
     AuxArray(111,IP) = Xpa*AuxArray(83,IP) + alphaXpq*TmpArray7(83,2)
     AuxArray(112,IP) = Xpa*AuxArray(84,IP) + alphaXpq*TmpArray7(84,2)
     AuxArray(113,IP) = Ypa*AuxArray(78,IP) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     AuxArray(114,IP) = Zpa*AuxArray(78,IP) + alphaZpq*TmpArray7(78,2)
     AuxArray(115,IP) = Ypa*AuxArray(80,IP) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     AuxArray(116,IP) = Ypa*AuxArray(81,IP) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     AuxArray(117,IP) = Ypa*AuxArray(82,IP) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     AuxArray(118,IP) = Ypa*AuxArray(83,IP) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     AuxArray(119,IP) = Ypa*AuxArray(84,IP) + alphaYpq*TmpArray7(84,2)
     AuxArray(120,IP) = Zpa*AuxArray(84,IP) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU8A(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & RJ000Array,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0:11,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA)
  real(realk),intent(inout) :: AUXarray(  165,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa
  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP
  real(realk) :: TwoTerms(  28)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
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
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,invexpP,inv2expP,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Ax = -Acenter(1,iAtomA)
   Ay = -Acenter(2,iAtomA)
   Az = -Acenter(3,iAtomA)
   iP = (iPassP-1)*nPrimQ*nPrimP
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    invexpP = D1/Pexp(iPrimP)
    inv2expP = D05*invexpP
    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax
    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay
    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
     alphaXpq = -alphaP*Xpq
     alphaYpq = -alphaP*Ypq
     alphaZpq = -alphaP*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 6) = PREF*RJ000Array( 5,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 7) = PREF*RJ000Array( 6,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 8) = PREF*RJ000Array( 7,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 9) = PREF*RJ000Array( 8,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xpa*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Ypa*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zpa*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expP*(AuxArray(1,IP) + alphaP*TmpArray1(1,2))
     AuxArray(5,IP) = Xpa*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xpa*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xpa*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Ypa*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Ypa*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zpa*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expP*(AuxArray(2,IP) + alphaP*TmpArray2(2,2))
     TwoTerms(2) = inv2expP*(AuxArray(3,IP) + alphaP*TmpArray2(3,2))
     TwoTerms(3) = inv2expP*(AuxArray(4,IP) + alphaP*TmpArray2(4,2))
     AuxArray(11,IP) = Xpa*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Ypa*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zpa*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xpa*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xpa*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xpa*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Ypa*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zpa*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Ypa*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zpa*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(5,IP) + alphaP*TmpArray3(5,2))
     TwoTerms(2) = inv2expP*(AuxArray(8,IP) + alphaP*TmpArray3(8,2))
     TwoTerms(3) = inv2expP*(AuxArray(10,IP) + alphaP*TmpArray3(10,2))
     AuxArray(21,IP) = Xpa*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Ypa*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zpa*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xpa*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Ypa*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xpa*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xpa*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xpa*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xpa*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xpa*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Ypa*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zpa*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Ypa*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Ypa*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zpa*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expP*(AuxArray(11,IP) + alphaP*TmpArray4(11,2))
     TwoTerms(2) = inv2expP*(AuxArray(14,IP) + alphaP*TmpArray4(14,2))
     TwoTerms(3) = inv2expP*(AuxArray(16,IP) + alphaP*TmpArray4(16,2))
     TwoTerms(4) = inv2expP*(AuxArray(17,IP) + alphaP*TmpArray4(17,2))
     TwoTerms(5) = inv2expP*(AuxArray(19,IP) + alphaP*TmpArray4(19,2))
     TwoTerms(6) = inv2expP*(AuxArray(20,IP) + alphaP*TmpArray4(20,2))
     AuxArray(36,IP) = Xpa*AuxArray(21,IP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IP) = Ypa*AuxArray(21,IP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IP) = Zpa*AuxArray(21,IP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IP) = Xpa*AuxArray(24,IP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IP) = Ypa*AuxArray(23,IP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IP) = Xpa*AuxArray(26,IP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IP) = Xpa*AuxArray(27,IP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IP) = Zpa*AuxArray(24,IP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IP) = Ypa*AuxArray(26,IP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IP) = Xpa*AuxArray(30,IP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IP) = Xpa*AuxArray(31,IP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IP) = Xpa*AuxArray(32,IP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IP) = Xpa*AuxArray(33,IP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IP) = Xpa*AuxArray(34,IP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IP) = Xpa*AuxArray(35,IP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IP) = Ypa*AuxArray(31,IP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IP) = Zpa*AuxArray(31,IP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IP) = Ypa*AuxArray(33,IP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IP) = Ypa*AuxArray(34,IP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IP) = Ypa*AuxArray(35,IP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IP) = Zpa*AuxArray(35,IP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     TwoTerms(1) = inv2expP*(AuxArray(21,IP) + alphaP*TmpArray5(21,2))
     TwoTerms(2) = inv2expP*(AuxArray(24,IP) + alphaP*TmpArray5(24,2))
     TwoTerms(3) = inv2expP*(AuxArray(26,IP) + alphaP*TmpArray5(26,2))
     TwoTerms(4) = inv2expP*(AuxArray(27,IP) + alphaP*TmpArray5(27,2))
     TwoTerms(5) = inv2expP*(AuxArray(30,IP) + alphaP*TmpArray5(30,2))
     TwoTerms(6) = inv2expP*(AuxArray(31,IP) + alphaP*TmpArray5(31,2))
     TwoTerms(7) = inv2expP*(AuxArray(33,IP) + alphaP*TmpArray5(33,2))
     TwoTerms(8) = inv2expP*(AuxArray(34,IP) + alphaP*TmpArray5(34,2))
     TwoTerms(9) = inv2expP*(AuxArray(35,IP) + alphaP*TmpArray5(35,2))
     AuxArray(57,IP) = Xpa*AuxArray(36,IP) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,IP) = Ypa*AuxArray(36,IP) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,IP) = Zpa*AuxArray(36,IP) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,IP) = Xpa*AuxArray(39,IP) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,IP) = Ypa*AuxArray(38,IP) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,IP) = Xpa*AuxArray(41,IP) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,IP) = Xpa*AuxArray(42,IP) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,IP) = Zpa*AuxArray(39,IP) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,IP) = Ypa*AuxArray(41,IP) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,IP) = Xpa*AuxArray(45,IP) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,IP) = Xpa*AuxArray(46,IP) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,IP) = Zpa*AuxArray(42,IP) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,IP) = Xpa*AuxArray(48,IP) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,IP) = Ypa*AuxArray(45,IP) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,IP) = Xpa*AuxArray(50,IP) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,IP) = Xpa*AuxArray(51,IP) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,IP) = Xpa*AuxArray(52,IP) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,IP) = Xpa*AuxArray(53,IP) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,IP) = Xpa*AuxArray(54,IP) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,IP) = Xpa*AuxArray(55,IP) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,IP) = Xpa*AuxArray(56,IP) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,IP) = Ypa*AuxArray(51,IP) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,IP) = Zpa*AuxArray(51,IP) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,IP) = Ypa*AuxArray(53,IP) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,IP) = Ypa*AuxArray(54,IP) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,IP) = Ypa*AuxArray(55,IP) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,IP) = Ypa*AuxArray(56,IP) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,IP) = Zpa*AuxArray(56,IP) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
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
     TwoTerms(1) = inv2expP*(AuxArray(36,IP) + alphaP*TmpArray6(36,2))
     TwoTerms(2) = inv2expP*(AuxArray(39,IP) + alphaP*TmpArray6(39,2))
     TwoTerms(3) = inv2expP*(AuxArray(41,IP) + alphaP*TmpArray6(41,2))
     TwoTerms(4) = inv2expP*(AuxArray(42,IP) + alphaP*TmpArray6(42,2))
     TwoTerms(5) = inv2expP*(AuxArray(45,IP) + alphaP*TmpArray6(45,2))
     TwoTerms(6) = inv2expP*(AuxArray(46,IP) + alphaP*TmpArray6(46,2))
     TwoTerms(7) = inv2expP*(AuxArray(48,IP) + alphaP*TmpArray6(48,2))
     TwoTerms(8) = inv2expP*(AuxArray(50,IP) + alphaP*TmpArray6(50,2))
     TwoTerms(9) = inv2expP*(AuxArray(51,IP) + alphaP*TmpArray6(51,2))
     TwoTerms(10) = inv2expP*(AuxArray(53,IP) + alphaP*TmpArray6(53,2))
     TwoTerms(11) = inv2expP*(AuxArray(54,IP) + alphaP*TmpArray6(54,2))
     TwoTerms(12) = inv2expP*(AuxArray(55,IP) + alphaP*TmpArray6(55,2))
     TwoTerms(13) = inv2expP*(AuxArray(56,IP) + alphaP*TmpArray6(56,2))
     AuxArray(85,IP) = Xpa*AuxArray(57,IP) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     AuxArray(86,IP) = Ypa*AuxArray(57,IP) + alphaYpq*TmpArray7(57,2)
     AuxArray(87,IP) = Zpa*AuxArray(57,IP) + alphaZpq*TmpArray7(57,2)
     AuxArray(88,IP) = Xpa*AuxArray(60,IP) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     AuxArray(89,IP) = Ypa*AuxArray(59,IP) + alphaYpq*TmpArray7(59,2)
     AuxArray(90,IP) = Xpa*AuxArray(62,IP) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     AuxArray(91,IP) = Xpa*AuxArray(63,IP) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     AuxArray(92,IP) = Zpa*AuxArray(60,IP) + alphaZpq*TmpArray7(60,2)
     AuxArray(93,IP) = Ypa*AuxArray(62,IP) + alphaYpq*TmpArray7(62,2)
     AuxArray(94,IP) = Xpa*AuxArray(66,IP) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     AuxArray(95,IP) = Xpa*AuxArray(67,IP) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     AuxArray(96,IP) = Zpa*AuxArray(63,IP) + alphaZpq*TmpArray7(63,2)
     AuxArray(97,IP) = Xpa*AuxArray(69,IP) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     AuxArray(98,IP) = Ypa*AuxArray(66,IP) + alphaYpq*TmpArray7(66,2)
     AuxArray(99,IP) = Xpa*AuxArray(71,IP) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     AuxArray(100,IP) = Xpa*AuxArray(72,IP) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     AuxArray(101,IP) = Zpa*AuxArray(67,IP) + alphaZpq*TmpArray7(67,2)
     AuxArray(102,IP) = Xpa*AuxArray(74,IP) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     AuxArray(103,IP) = Xpa*AuxArray(75,IP) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     AuxArray(104,IP) = Ypa*AuxArray(71,IP) + alphaYpq*TmpArray7(71,2)
     AuxArray(105,IP) = Xpa*AuxArray(77,IP) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     AuxArray(106,IP) = Xpa*AuxArray(78,IP) + alphaXpq*TmpArray7(78,2)
     AuxArray(107,IP) = Xpa*AuxArray(79,IP) + alphaXpq*TmpArray7(79,2)
     AuxArray(108,IP) = Xpa*AuxArray(80,IP) + alphaXpq*TmpArray7(80,2)
     AuxArray(109,IP) = Xpa*AuxArray(81,IP) + alphaXpq*TmpArray7(81,2)
     AuxArray(110,IP) = Xpa*AuxArray(82,IP) + alphaXpq*TmpArray7(82,2)
     AuxArray(111,IP) = Xpa*AuxArray(83,IP) + alphaXpq*TmpArray7(83,2)
     AuxArray(112,IP) = Xpa*AuxArray(84,IP) + alphaXpq*TmpArray7(84,2)
     AuxArray(113,IP) = Ypa*AuxArray(78,IP) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     AuxArray(114,IP) = Zpa*AuxArray(78,IP) + alphaZpq*TmpArray7(78,2)
     AuxArray(115,IP) = Ypa*AuxArray(80,IP) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     AuxArray(116,IP) = Ypa*AuxArray(81,IP) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     AuxArray(117,IP) = Ypa*AuxArray(82,IP) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     AuxArray(118,IP) = Ypa*AuxArray(83,IP) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     AuxArray(119,IP) = Ypa*AuxArray(84,IP) + alphaYpq*TmpArray7(84,2)
     AuxArray(120,IP) = Zpa*AuxArray(84,IP) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
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
     TwoTerms(1) = inv2expP*(AuxArray(57,IP) + alphaP*TmpArray7(57,2))
     TwoTerms(2) = inv2expP*(AuxArray(60,IP) + alphaP*TmpArray7(60,2))
     TwoTerms(3) = inv2expP*(AuxArray(62,IP) + alphaP*TmpArray7(62,2))
     TwoTerms(4) = inv2expP*(AuxArray(63,IP) + alphaP*TmpArray7(63,2))
     TwoTerms(5) = inv2expP*(AuxArray(66,IP) + alphaP*TmpArray7(66,2))
     TwoTerms(6) = inv2expP*(AuxArray(67,IP) + alphaP*TmpArray7(67,2))
     TwoTerms(7) = inv2expP*(AuxArray(69,IP) + alphaP*TmpArray7(69,2))
     TwoTerms(8) = inv2expP*(AuxArray(71,IP) + alphaP*TmpArray7(71,2))
     TwoTerms(9) = inv2expP*(AuxArray(72,IP) + alphaP*TmpArray7(72,2))
     TwoTerms(10) = inv2expP*(AuxArray(74,IP) + alphaP*TmpArray7(74,2))
     TwoTerms(11) = inv2expP*(AuxArray(75,IP) + alphaP*TmpArray7(75,2))
     TwoTerms(12) = inv2expP*(AuxArray(77,IP) + alphaP*TmpArray7(77,2))
     TwoTerms(13) = inv2expP*(AuxArray(78,IP) + alphaP*TmpArray7(78,2))
     TwoTerms(14) = inv2expP*(AuxArray(80,IP) + alphaP*TmpArray7(80,2))
     TwoTerms(15) = inv2expP*(AuxArray(81,IP) + alphaP*TmpArray7(81,2))
     TwoTerms(16) = inv2expP*(AuxArray(82,IP) + alphaP*TmpArray7(82,2))
     TwoTerms(17) = inv2expP*(AuxArray(83,IP) + alphaP*TmpArray7(83,2))
     TwoTerms(18) = inv2expP*(AuxArray(84,IP) + alphaP*TmpArray7(84,2))
     AuxArray(121,IP) = Xpa*AuxArray(85,IP) + alphaXpq*TmpArray8(85,2) + 7*TwoTerms(1)
     AuxArray(122,IP) = Ypa*AuxArray(85,IP) + alphaYpq*TmpArray8(85,2)
     AuxArray(123,IP) = Zpa*AuxArray(85,IP) + alphaZpq*TmpArray8(85,2)
     AuxArray(124,IP) = Xpa*AuxArray(88,IP) + alphaXpq*TmpArray8(88,2) + 5*TwoTerms(2)
     AuxArray(125,IP) = Ypa*AuxArray(87,IP) + alphaYpq*TmpArray8(87,2)
     AuxArray(126,IP) = Xpa*AuxArray(90,IP) + alphaXpq*TmpArray8(90,2) + 5*TwoTerms(3)
     AuxArray(127,IP) = Xpa*AuxArray(91,IP) + alphaXpq*TmpArray8(91,2) + 4*TwoTerms(4)
     AuxArray(128,IP) = Zpa*AuxArray(88,IP) + alphaZpq*TmpArray8(88,2)
     AuxArray(129,IP) = Ypa*AuxArray(90,IP) + alphaYpq*TmpArray8(90,2)
     AuxArray(130,IP) = Xpa*AuxArray(94,IP) + alphaXpq*TmpArray8(94,2) + 4*TwoTerms(5)
     AuxArray(131,IP) = Xpa*AuxArray(95,IP) + alphaXpq*TmpArray8(95,2) + 3*TwoTerms(6)
     AuxArray(132,IP) = Zpa*AuxArray(91,IP) + alphaZpq*TmpArray8(91,2)
     AuxArray(133,IP) = Xpa*AuxArray(97,IP) + alphaXpq*TmpArray8(97,2) + 3*TwoTerms(7)
     AuxArray(134,IP) = Ypa*AuxArray(94,IP) + alphaYpq*TmpArray8(94,2)
     AuxArray(135,IP) = Xpa*AuxArray(99,IP) + alphaXpq*TmpArray8(99,2) + 3*TwoTerms(8)
     AuxArray(136,IP) = Xpa*AuxArray(100,IP) + alphaXpq*TmpArray8(100,2) + 2*TwoTerms(9)
     AuxArray(137,IP) = Zpa*AuxArray(95,IP) + alphaZpq*TmpArray8(95,2)
     AuxArray(138,IP) = Xpa*AuxArray(102,IP) + alphaXpq*TmpArray8(102,2) + 2*TwoTerms(10)
     AuxArray(139,IP) = Xpa*AuxArray(103,IP) + alphaXpq*TmpArray8(103,2) + 2*TwoTerms(11)
     AuxArray(140,IP) = Ypa*AuxArray(99,IP) + alphaYpq*TmpArray8(99,2)
     AuxArray(141,IP) = Xpa*AuxArray(105,IP) + alphaXpq*TmpArray8(105,2) + 2*TwoTerms(12)
     AuxArray(142,IP) = Xpa*AuxArray(106,IP) + alphaXpq*TmpArray8(106,2) + TwoTerms(13)
     AuxArray(143,IP) = Zpa*AuxArray(100,IP) + alphaZpq*TmpArray8(100,2)
     AuxArray(144,IP) = Xpa*AuxArray(108,IP) + alphaXpq*TmpArray8(108,2) + TwoTerms(14)
     AuxArray(145,IP) = Xpa*AuxArray(109,IP) + alphaXpq*TmpArray8(109,2) + TwoTerms(15)
     AuxArray(146,IP) = Xpa*AuxArray(110,IP) + alphaXpq*TmpArray8(110,2) + TwoTerms(16)
     AuxArray(147,IP) = Ypa*AuxArray(105,IP) + alphaYpq*TmpArray8(105,2)
     AuxArray(148,IP) = Xpa*AuxArray(112,IP) + alphaXpq*TmpArray8(112,2) + TwoTerms(18)
     AuxArray(149,IP) = Xpa*AuxArray(113,IP) + alphaXpq*TmpArray8(113,2)
     AuxArray(150,IP) = Xpa*AuxArray(114,IP) + alphaXpq*TmpArray8(114,2)
     AuxArray(151,IP) = Xpa*AuxArray(115,IP) + alphaXpq*TmpArray8(115,2)
     AuxArray(152,IP) = Xpa*AuxArray(116,IP) + alphaXpq*TmpArray8(116,2)
     AuxArray(153,IP) = Xpa*AuxArray(117,IP) + alphaXpq*TmpArray8(117,2)
     AuxArray(154,IP) = Xpa*AuxArray(118,IP) + alphaXpq*TmpArray8(118,2)
     AuxArray(155,IP) = Xpa*AuxArray(119,IP) + alphaXpq*TmpArray8(119,2)
     AuxArray(156,IP) = Xpa*AuxArray(120,IP) + alphaXpq*TmpArray8(120,2)
     AuxArray(157,IP) = Ypa*AuxArray(113,IP) + alphaYpq*TmpArray8(113,2) + 7*TwoTerms(13)
     AuxArray(158,IP) = Zpa*AuxArray(113,IP) + alphaZpq*TmpArray8(113,2)
     AuxArray(159,IP) = Ypa*AuxArray(115,IP) + alphaYpq*TmpArray8(115,2) + 5*TwoTerms(14)
     AuxArray(160,IP) = Ypa*AuxArray(116,IP) + alphaYpq*TmpArray8(116,2) + 4*TwoTerms(15)
     AuxArray(161,IP) = Ypa*AuxArray(117,IP) + alphaYpq*TmpArray8(117,2) + 3*TwoTerms(16)
     AuxArray(162,IP) = Ypa*AuxArray(118,IP) + alphaYpq*TmpArray8(118,2) + 2*TwoTerms(17)
     AuxArray(163,IP) = Ypa*AuxArray(119,IP) + alphaYpq*TmpArray8(119,2) + TwoTerms(18)
     AuxArray(164,IP) = Ypa*AuxArray(120,IP) + alphaYpq*TmpArray8(120,2)
     AuxArray(165,IP) = Zpa*AuxArray(120,IP) + alphaZpq*TmpArray8(120,2) + 7*TwoTerms(18)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine
end module
