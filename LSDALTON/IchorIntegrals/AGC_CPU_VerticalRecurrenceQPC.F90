MODULE AGC_CPU_OBS_VERTICALRECURRENCEMODC
 use IchorPrecisionModule
  
 CONTAINS

subroutine VerticalRecurrenceCPU1C(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Ccenter,Pcent,Qcent,&
         & integralPrefactor,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: TABFJW(0:4,0:1200)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(4,nPrimQ,nPrimP,nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,iAtomA,iAtomB
  real(realk) :: Cx,Cy,Cz,Xqc,Yqc,Zqc
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
  !ThetaAux(n,1,0,0) = Xqc*ThetaAux(n,0,0,0) + (-alpha/q*Xpq)*ThetaAux(n+1,0,0,0)
  !i = 0 last 2 term vanish
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Cx,Cy,Cz,Xqc,Yqc,Zqc,alphaQ,&
!$OMP         Pexpfac,mPX,mPY,mPZ,RWVAL,GVAL,&
!$OMP         alphaXpq,alphaYpq,alphaZpq) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   Cx = -Ccenter(1)
   Cy = -Ccenter(2)
   Cz = -Ccenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     Xqc = Qcent(1,iPrimQ) + Cx
     Yqc = Qcent(2,iPrimQ) + Cy
     Zqc = Qcent(3,iPrimQ) + Cz
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
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     TMP1 = PREF*RJ000(0)
     TMP2 = PREF*RJ000(1)
     AUXarray(1,iPrimQ,iPrimP,iPassP) = TMP1
     AUXarray(2,iPrimQ,iPrimP,iPassP) = Xqc*TMP1 + alphaXpq*TMP2
     AUXarray(3,iPrimQ,iPrimP,iPassP) = Yqc*TMP1 + alphaYpq*TMP2
     AUXarray(4,iPrimQ,iPrimP,iPassP) = Zqc*TMP1 + alphaZpq*TMP2
    enddo
   enddo
  enddo
!$OMP END DO
end subroutine

subroutine VerticalRecurrenceCPU2C(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 5,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(   10,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Cx,Cy,Cz,Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ
  real(realk) :: TwoTerms(   1)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:3)
  real(realk) :: TMParray2(  2:  4,2:2)
  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Cx,Cy,Cz,Xqc,Yqc,Zqc,alphaQ,inv2expQ,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   iP = (iPassP-1)*nPrimQ*nPrimP
   Cx = -Ccenter(1)
   Cy = -Ccenter(2)
   Cz = -Ccenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     Xqc = Qcent(1,iPrimQ) + Cx
     Yqc = Qcent(2,iPrimQ) + Cy
     Zqc = Qcent(3,iPrimQ) + Cz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xqc*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqc*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqc*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqc*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqc*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqc*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,IP) = Xqc*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqc*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqc*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqc*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqc*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqc*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU3C(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 6,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(   20,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Cx,Cy,Cz,Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ
  real(realk) :: TwoTerms(   3)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:4)
  real(realk) :: TMParray2(  2:  4,2:3)
  real(realk) :: TMParray3(  5: 10,2:2)
  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Cx,Cy,Cz,Xqc,Yqc,Zqc,alphaQ,inv2expQ,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   iP = (iPassP-1)*nPrimQ*nPrimP
   Cx = -Ccenter(1)
   Cy = -Ccenter(2)
   Cz = -Ccenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     Xqc = Qcent(1,iPrimQ) + Cx
     Yqc = Qcent(2,iPrimQ) + Cy
     Zqc = Qcent(3,iPrimQ) + Cz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xqc*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqc*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqc*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqc*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqc*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqc*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqc*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqc*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqc*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expQ*(AuxArray(1,IP) + alphaQ*TmpArray1(1,2))
     AuxArray(5,IP) = Xqc*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqc*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqc*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqc*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqc*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqc*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     AuxArray(11,IP) = Xqc*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Yqc*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zqc*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xqc*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xqc*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xqc*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Yqc*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zqc*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Yqc*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zqc*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU4C(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 7,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(   35,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Cx,Cy,Cz,Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ
  real(realk) :: TwoTerms(   6)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D1=1.0E0_realk
  real(realk) :: TMParray1(  1:  1,2:5)
  real(realk) :: TMParray2(  2:  4,2:4)
  real(realk) :: TMParray3(  5: 10,2:3)
  real(realk) :: TMParray4( 11: 20,2:2)
  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Cx,Cy,Cz,Xqc,Yqc,Zqc,alphaQ,inv2expQ,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   iP = (iPassP-1)*nPrimQ*nPrimP
   Cx = -Ccenter(1)
   Cy = -Ccenter(2)
   Cz = -Ccenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     Xqc = Qcent(1,iPrimQ) + Cx
     Yqc = Qcent(2,iPrimQ) + Cy
     Zqc = Qcent(3,iPrimQ) + Cz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xqc*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqc*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqc*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     AuxArray(5,IP) = Xqc*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqc*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqc*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqc*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqc*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqc*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     AuxArray(11,IP) = Xqc*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Yqc*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zqc*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xqc*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xqc*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xqc*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Yqc*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zqc*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Yqc*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zqc*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     AuxArray(21,IP) = Xqc*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Yqc*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zqc*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xqc*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Yqc*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xqc*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xqc*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xqc*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xqc*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xqc*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Yqc*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zqc*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Yqc*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Yqc*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zqc*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU5C(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 8,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(   56,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Cx,Cy,Cz,Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ
  real(realk) :: TwoTerms(  10)
  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
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
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Cx,Cy,Cz,Xqc,Yqc,Zqc,alphaQ,inv2expQ,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   iP = (iPassP-1)*nPrimQ*nPrimP
   Cx = -Ccenter(1)
   Cy = -Ccenter(2)
   Cz = -Ccenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     Xqc = Qcent(1,iPrimQ) + Cx
     Yqc = Qcent(2,iPrimQ) + Cy
     Zqc = Qcent(3,iPrimQ) + Cz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 6) = PREF*RJ000Array( 5,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xqc*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqc*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqc*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     AuxArray(5,IP) = Xqc*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqc*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqc*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqc*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqc*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqc*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     AuxArray(11,IP) = Xqc*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Yqc*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zqc*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xqc*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xqc*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xqc*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Yqc*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zqc*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Yqc*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zqc*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     AuxArray(21,IP) = Xqc*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Yqc*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zqc*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xqc*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Yqc*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xqc*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xqc*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xqc*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xqc*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xqc*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Yqc*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zqc*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Yqc*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Yqc*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zqc*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     AuxArray(36,IP) = Xqc*AuxArray(21,IP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IP) = Yqc*AuxArray(21,IP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IP) = Zqc*AuxArray(21,IP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IP) = Xqc*AuxArray(24,IP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IP) = Yqc*AuxArray(23,IP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IP) = Xqc*AuxArray(26,IP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IP) = Xqc*AuxArray(27,IP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IP) = Zqc*AuxArray(24,IP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IP) = Yqc*AuxArray(26,IP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IP) = Xqc*AuxArray(30,IP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IP) = Xqc*AuxArray(31,IP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IP) = Xqc*AuxArray(32,IP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IP) = Xqc*AuxArray(33,IP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IP) = Xqc*AuxArray(34,IP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IP) = Xqc*AuxArray(35,IP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IP) = Yqc*AuxArray(31,IP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IP) = Zqc*AuxArray(31,IP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IP) = Yqc*AuxArray(33,IP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IP) = Yqc*AuxArray(34,IP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IP) = Yqc*AuxArray(35,IP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IP) = Zqc*AuxArray(35,IP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU6C(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0: 9,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(   84,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Cx,Cy,Cz,Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ
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
  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Cx,Cy,Cz,Xqc,Yqc,Zqc,alphaQ,inv2expQ,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   iP = (iPassP-1)*nPrimQ*nPrimP
   Cx = -Ccenter(1)
   Cy = -Ccenter(2)
   Cz = -Ccenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     Xqc = Qcent(1,iPrimQ) + Cx
     Yqc = Qcent(2,iPrimQ) + Cy
     Zqc = Qcent(3,iPrimQ) + Cz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 6) = PREF*RJ000Array( 5,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 7) = PREF*RJ000Array( 6,iPrimQ,iPrimP,iPassP)
     AuxArray(2,IP) = Xqc*AuxArray(1,IP) + alphaXpq*TmpArray1(1,2)
     AuxArray(3,IP) = Yqc*AuxArray(1,IP) + alphaYpq*TmpArray1(1,2)
     AuxArray(4,IP) = Zqc*AuxArray(1,IP) + alphaZpq*TmpArray1(1,2)
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
     AuxArray(5,IP) = Xqc*AuxArray(2,IP) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(6,IP) = Xqc*AuxArray(3,IP) + alphaXpq*TmpArray2(3,2)
     AuxArray(7,IP) = Xqc*AuxArray(4,IP) + alphaXpq*TmpArray2(4,2)
     AuxArray(8,IP) = Yqc*AuxArray(3,IP) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(9,IP) = Yqc*AuxArray(4,IP) + alphaYpq*TmpArray2(4,2)
     AuxArray(10,IP) = Zqc*AuxArray(4,IP) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     AuxArray(11,IP) = Xqc*AuxArray(5,IP) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(12,IP) = Yqc*AuxArray(5,IP) + alphaYpq*TmpArray3(5,2)
     AuxArray(13,IP) = Zqc*AuxArray(5,IP) + alphaZpq*TmpArray3(5,2)
     AuxArray(14,IP) = Xqc*AuxArray(8,IP) + alphaXpq*TmpArray3(8,2)
     AuxArray(15,IP) = Xqc*AuxArray(9,IP) + alphaXpq*TmpArray3(9,2)
     AuxArray(16,IP) = Xqc*AuxArray(10,IP) + alphaXpq*TmpArray3(10,2)
     AuxArray(17,IP) = Yqc*AuxArray(8,IP) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(18,IP) = Zqc*AuxArray(8,IP) + alphaZpq*TmpArray3(8,2)
     AuxArray(19,IP) = Yqc*AuxArray(10,IP) + alphaYpq*TmpArray3(10,2)
     AuxArray(20,IP) = Zqc*AuxArray(10,IP) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     AuxArray(21,IP) = Xqc*AuxArray(11,IP) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(22,IP) = Yqc*AuxArray(11,IP) + alphaYpq*TmpArray4(11,2)
     AuxArray(23,IP) = Zqc*AuxArray(11,IP) + alphaZpq*TmpArray4(11,2)
     AuxArray(24,IP) = Xqc*AuxArray(14,IP) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(25,IP) = Yqc*AuxArray(13,IP) + alphaYpq*TmpArray4(13,2)
     AuxArray(26,IP) = Xqc*AuxArray(16,IP) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(27,IP) = Xqc*AuxArray(17,IP) + alphaXpq*TmpArray4(17,2)
     AuxArray(28,IP) = Xqc*AuxArray(18,IP) + alphaXpq*TmpArray4(18,2)
     AuxArray(29,IP) = Xqc*AuxArray(19,IP) + alphaXpq*TmpArray4(19,2)
     AuxArray(30,IP) = Xqc*AuxArray(20,IP) + alphaXpq*TmpArray4(20,2)
     AuxArray(31,IP) = Yqc*AuxArray(17,IP) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(32,IP) = Zqc*AuxArray(17,IP) + alphaZpq*TmpArray4(17,2)
     AuxArray(33,IP) = Yqc*AuxArray(19,IP) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(34,IP) = Yqc*AuxArray(20,IP) + alphaYpq*TmpArray4(20,2)
     AuxArray(35,IP) = Zqc*AuxArray(20,IP) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     AuxArray(36,IP) = Xqc*AuxArray(21,IP) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(37,IP) = Yqc*AuxArray(21,IP) + alphaYpq*TmpArray5(21,2)
     AuxArray(38,IP) = Zqc*AuxArray(21,IP) + alphaZpq*TmpArray5(21,2)
     AuxArray(39,IP) = Xqc*AuxArray(24,IP) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(40,IP) = Yqc*AuxArray(23,IP) + alphaYpq*TmpArray5(23,2)
     AuxArray(41,IP) = Xqc*AuxArray(26,IP) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(42,IP) = Xqc*AuxArray(27,IP) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(43,IP) = Zqc*AuxArray(24,IP) + alphaZpq*TmpArray5(24,2)
     AuxArray(44,IP) = Yqc*AuxArray(26,IP) + alphaYpq*TmpArray5(26,2)
     AuxArray(45,IP) = Xqc*AuxArray(30,IP) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(46,IP) = Xqc*AuxArray(31,IP) + alphaXpq*TmpArray5(31,2)
     AuxArray(47,IP) = Xqc*AuxArray(32,IP) + alphaXpq*TmpArray5(32,2)
     AuxArray(48,IP) = Xqc*AuxArray(33,IP) + alphaXpq*TmpArray5(33,2)
     AuxArray(49,IP) = Xqc*AuxArray(34,IP) + alphaXpq*TmpArray5(34,2)
     AuxArray(50,IP) = Xqc*AuxArray(35,IP) + alphaXpq*TmpArray5(35,2)
     AuxArray(51,IP) = Yqc*AuxArray(31,IP) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(52,IP) = Zqc*AuxArray(31,IP) + alphaZpq*TmpArray5(31,2)
     AuxArray(53,IP) = Yqc*AuxArray(33,IP) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(54,IP) = Yqc*AuxArray(34,IP) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(55,IP) = Yqc*AuxArray(35,IP) + alphaYpq*TmpArray5(35,2)
     AuxArray(56,IP) = Zqc*AuxArray(35,IP) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     AuxArray(57,IP) = Xqc*AuxArray(36,IP) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(58,IP) = Yqc*AuxArray(36,IP) + alphaYpq*TmpArray6(36,2)
     AuxArray(59,IP) = Zqc*AuxArray(36,IP) + alphaZpq*TmpArray6(36,2)
     AuxArray(60,IP) = Xqc*AuxArray(39,IP) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(61,IP) = Yqc*AuxArray(38,IP) + alphaYpq*TmpArray6(38,2)
     AuxArray(62,IP) = Xqc*AuxArray(41,IP) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(63,IP) = Xqc*AuxArray(42,IP) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(64,IP) = Zqc*AuxArray(39,IP) + alphaZpq*TmpArray6(39,2)
     AuxArray(65,IP) = Yqc*AuxArray(41,IP) + alphaYpq*TmpArray6(41,2)
     AuxArray(66,IP) = Xqc*AuxArray(45,IP) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(67,IP) = Xqc*AuxArray(46,IP) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(68,IP) = Zqc*AuxArray(42,IP) + alphaZpq*TmpArray6(42,2)
     AuxArray(69,IP) = Xqc*AuxArray(48,IP) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(70,IP) = Yqc*AuxArray(45,IP) + alphaYpq*TmpArray6(45,2)
     AuxArray(71,IP) = Xqc*AuxArray(50,IP) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(72,IP) = Xqc*AuxArray(51,IP) + alphaXpq*TmpArray6(51,2)
     AuxArray(73,IP) = Xqc*AuxArray(52,IP) + alphaXpq*TmpArray6(52,2)
     AuxArray(74,IP) = Xqc*AuxArray(53,IP) + alphaXpq*TmpArray6(53,2)
     AuxArray(75,IP) = Xqc*AuxArray(54,IP) + alphaXpq*TmpArray6(54,2)
     AuxArray(76,IP) = Xqc*AuxArray(55,IP) + alphaXpq*TmpArray6(55,2)
     AuxArray(77,IP) = Xqc*AuxArray(56,IP) + alphaXpq*TmpArray6(56,2)
     AuxArray(78,IP) = Yqc*AuxArray(51,IP) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(79,IP) = Zqc*AuxArray(51,IP) + alphaZpq*TmpArray6(51,2)
     AuxArray(80,IP) = Yqc*AuxArray(53,IP) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(81,IP) = Yqc*AuxArray(54,IP) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(82,IP) = Yqc*AuxArray(55,IP) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(83,IP) = Yqc*AuxArray(56,IP) + alphaYpq*TmpArray6(56,2)
     AuxArray(84,IP) = Zqc*AuxArray(56,IP) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
 end subroutine

subroutine VerticalRecurrenceCPU7C(nPasses,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPasses,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  REAL(REALK),intent(in) :: RJ000Array(0:10,nPrimQ,nPrimP,nPasses)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)
  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)
  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(  120,nPrimQ*nPrimP*nPasses)
  !local variables
  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: Cx,Cy,Cz,Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ
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
  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)
  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))
  !We include scaling of RJ000 
  DO iPrimQ=1, nPrimQ
     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)
  ENDDO
!$OMP DO &
!$OMP PRIVATE(iPassP,iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$OMP         iPrimP,iPrimQ,&
!$OMP         Cx,Cy,Cz,Xqc,Yqc,Zqc,alphaQ,inv2expQ,&
!$OMP         Pexpfac,mPX,mPY,mPZ,&
!$OMP         alphaXpq,alphaYpq,alphaZpq,TwoTerms,iP) 
  DO iPassP = 1,nPasses
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
   iP = (iPassP-1)*nPrimQ*nPrimP
   Cx = -Ccenter(1)
   Cy = -Ccenter(2)
   Cz = -Ccenter(3)
   DO iPrimP=1, nPrimP
    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)
    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)
    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)
    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)
    DO iPrimQ=1, nPrimQ
     iP = iP + 1
     Xpq = mPX + Qcent(1,iPrimQ)
     Ypq = mPY + Qcent(2,iPrimQ)
     Zpq = mPZ + Qcent(3,iPrimQ)
     Xqc = Qcent(1,iPrimQ) + Cx
     Yqc = Qcent(2,iPrimQ) + Cy
     Zqc = Qcent(3,iPrimQ) + Cz
     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)
     inv2expQ = D05*invexpQ(iPrimQ)
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*Pexpfac
     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 2) = PREF*RJ000Array( 1,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 3) = PREF*RJ000Array( 2,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 4) = PREF*RJ000Array( 3,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 5) = PREF*RJ000Array( 4,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 6) = PREF*RJ000Array( 5,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 7) = PREF*RJ000Array( 6,iPrimQ,iPrimP,iPassP)
     TMParray1(1, 8) = PREF*RJ000Array( 7,iPrimQ,iPrimP,iPassP)
