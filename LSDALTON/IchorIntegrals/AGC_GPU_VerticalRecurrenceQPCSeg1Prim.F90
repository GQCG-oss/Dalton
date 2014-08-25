MODULE AGC_GPU_OBS_VERTICALRECURRENCEMODCSeg1Prim
 use IchorPrecisionModule
  
 CONTAINS


subroutine VerticalRecurrenceGPUSeg1Prim1C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: TABFJW(0:4,0:1200)
  real(realk),intent(in) :: reducedExponents(1),Qexp(1)
  real(realk),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3),integralPrefactor(1),QpreExpFac(1),PpreExpFac(nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(nPassP,4)
  !local variables
  Integer :: iP,iPassP,ipnt,iAtomA,iAtomB
  real(realk) :: Xqc,Yqc,Zqc
  real(realk) :: mPX,mPY,mPZ,invexpQ,alphaQ,RJ000(0:1)
  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq
  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(realk),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk
  real(realk),parameter :: D2=2.0E0_realk
  real(realk),PARAMETER :: D2JP36=  3.8000000000000000E+01_realk
  real(realk),parameter :: D1=1.0E0_realk,D03333=1.0E0_realk/3.0E0_realk
  real(realk),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk
  real(realk),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
  real(realk),PARAMETER :: SMALL = 1E-15_realk,D12 = 12.0E0_realk
  real(realk), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
  real(realk), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
  real(realk), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
  real(realk), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
  real(realk), parameter :: PI=3.14159265358979323846E0_realk
  real(realk), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
  real(realk), PARAMETER :: SQRPIH = SQRTPI/D2
  real(realk), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
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
!$ACC         iP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,TABFJW,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Ccenter, &
!$ACC        nPrimP,nPrimQ,nPassP)
  DO iP = 1,nPassP
   iPassP = iP
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
    mPX = -Pcent(1,iAtomA,iAtomB)
    mPY = -Pcent(2,iAtomA,iAtomB)
    mPZ = -Pcent(3,iAtomA,iAtomB)
    invexpQ = D1/Qexp(1)
   Xpq = mPX + Qcent(1)
   Ypq = mPY + Qcent(2)
   Zpq = mPZ + Qcent(3)
     Xqc = Qcent(1) - Ccenter(1)
     Yqc = Qcent(2) - Ccenter(2)
     Zqc = Qcent(3) - Ccenter(3)
     alphaQ = -reducedExponents(1)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(1)*squaredDistance
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
     PREF = integralPrefactor(1)*QpreExpFac(1)*PpreExpFac(iAtomA,iAtomB)
     TMP1 = PREF*RJ000(0)
     TMP2 = PREF*RJ000(1)
     AUXarray(iP,1) = TMP1
     AUXarray(iP,2) = Xqc*TMP1 + alphaXpq*TMP2
     AUXarray(iP,3) = Yqc*TMP1 + alphaYpq*TMP2
     AUXarray(iP,4) = Zqc*TMP1 + alphaZpq*TMP2
  ENDDO !iP = 1,nPassP
end subroutine VerticalRecurrenceGPUSeg1Prim1C

subroutine VerticalRecurrenceGPUSeg1Prim2C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPassP,0: 2)
  real(realk),intent(in) :: reducedExponents(1),Qexp(1)
  real(realk),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3),integralPrefactor(1),QpreExpFac(1),PpreExpFac(nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(nPassP,   10)
  !local variables
  integer :: iPassP,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: TMPAUXarray(    4)
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
!$ACC         TMPAUXarray,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TwoTerms,&
!$ACC         iP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Ccenter, &
!$ACC        nPrimP,nPrimQ,nPassP)
  DO iP = 1,nPassP
   iPassP = iP
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     Xqc = Qcent(1) - Ccenter(1)
     Yqc = Qcent(2) - Ccenter(2)
     Zqc = Qcent(3) - Ccenter(3)
     alphaQ = -reducedExponents(1)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(1)*QpreExpFac(1)*PpreExpFac(iAtomA,iAtomB)
     TMPAuxarray(1) = PREF*RJ000Array(iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPassP, 2)
     TMPAuxArray(2) = Xqc*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqc*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqc*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqc*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqc*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqc*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(1) + alphaQ*TmpArray1(1,2))
     do iTUV = 1,    4
      AuxArray(iP,iTUV) = TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,5) = Xqc*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     AuxArray(iP,6) = Xqc*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     AuxArray(iP,7) = Xqc*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     AuxArray(iP,8) = Yqc*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     AuxArray(iP,9) = Yqc*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     AuxArray(iP,10) = Zqc*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
  ENDDO !iP = 1,nPassP
 end subroutine

subroutine VerticalRecurrenceGPUSeg1Prim3C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPassP,0: 3)
  real(realk),intent(in) :: reducedExponents(1),Qexp(1)
  real(realk),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3),integralPrefactor(1),QpreExpFac(1),PpreExpFac(nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(nPassP,   20)
  !local variables
  integer :: iPassP,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: TMPAUXarray(   10)
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
!$ACC         TMPAUXarray,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TwoTerms,&
!$ACC         iP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Ccenter, &
!$ACC        nPrimP,nPrimQ,nPassP)
  DO iP = 1,nPassP
   iPassP = iP
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     Xqc = Qcent(1) - Ccenter(1)
     Yqc = Qcent(2) - Ccenter(2)
     Zqc = Qcent(3) - Ccenter(3)
     alphaQ = -reducedExponents(1)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(1)*QpreExpFac(1)*PpreExpFac(iAtomA,iAtomB)
     TMPAuxarray(1) = PREF*RJ000Array(iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPassP, 3)
     TMPAuxArray(2) = Xqc*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqc*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqc*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqc*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqc*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqc*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqc*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqc*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqc*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(1) + alphaQ*TmpArray1(1,2))
     TMPAuxArray(5) = Xqc*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     TMPAuxArray(6) = Xqc*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     TMPAuxArray(7) = Xqc*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     TMPAuxArray(8) = Yqc*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     TMPAuxArray(9) = Yqc*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     TMPAuxArray(10) = Zqc*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TmpArray1(1,2) + alphaQ*TmpArray1(1,3))
     tmpArray3(5,2) = Xqc*tmpArray2(2,2) + alphaXpq*TmpArray2(2,3) + TwoTerms(1)
     tmpArray3(6,2) = Xqc*tmpArray2(3,2) + alphaXpq*TmpArray2(3,3)
     tmpArray3(7,2) = Xqc*tmpArray2(4,2) + alphaXpq*TmpArray2(4,3)
     tmpArray3(8,2) = Yqc*tmpArray2(3,2) + alphaYpq*TmpArray2(3,3) + TwoTerms(1)
     tmpArray3(9,2) = Yqc*tmpArray2(4,2) + alphaYpq*TmpArray2(4,3)
     tmpArray3(10,2) = Zqc*tmpArray2(4,2) + alphaZpq*TmpArray2(4,3) + TwoTerms(1)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(4) + alphaQ*TmpArray2(4,2))
     do iTUV = 1,   10
      AuxArray(iP,iTUV) = TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,11) = Xqc*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     AuxArray(iP,12) = Yqc*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     AuxArray(iP,13) = Zqc*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     AuxArray(iP,14) = Xqc*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     AuxArray(iP,15) = Xqc*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     AuxArray(iP,16) = Xqc*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     AuxArray(iP,17) = Yqc*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     AuxArray(iP,18) = Zqc*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     AuxArray(iP,19) = Yqc*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     AuxArray(iP,20) = Zqc*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
  ENDDO !iP = 1,nPassP
 end subroutine

subroutine VerticalRecurrenceGPUSeg1Prim4C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPassP,0: 4)
  real(realk),intent(in) :: reducedExponents(1),Qexp(1)
  real(realk),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3),integralPrefactor(1),QpreExpFac(1),PpreExpFac(nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(nPassP,   35)
  !local variables
  integer :: iPassP,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: TMPAUXarray(   20)
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
!$ACC         TMPAUXarray,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TMParray4,&
!$ACC         TwoTerms,&
!$ACC         iP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Ccenter, &
!$ACC        nPrimP,nPrimQ,nPassP)
  DO iP = 1,nPassP
   iPassP = iP
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     Xqc = Qcent(1) - Ccenter(1)
     Yqc = Qcent(2) - Ccenter(2)
     Zqc = Qcent(3) - Ccenter(3)
     alphaQ = -reducedExponents(1)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(1)*QpreExpFac(1)*PpreExpFac(iAtomA,iAtomB)
     TMPAuxarray(1) = PREF*RJ000Array(iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPassP, 3)
     TMParray1(1, 5) = PREF*RJ000Array(iPassP, 4)
     TMPAuxArray(2) = Xqc*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqc*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqc*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
     tmpArray2(2,2) = Xqc*tmpArray1(1,2) + alphaXpq*TmpArray1(1,3)
     tmpArray2(3,2) = Yqc*tmpArray1(1,2) + alphaYpq*TmpArray1(1,3)
     tmpArray2(4,2) = Zqc*tmpArray1(1,2) + alphaZpq*TmpArray1(1,3)
     tmpArray2(2,3) = Xqc*tmpArray1(1,3) + alphaXpq*TmpArray1(1,4)
     tmpArray2(3,3) = Yqc*tmpArray1(1,3) + alphaYpq*TmpArray1(1,4)
     tmpArray2(4,3) = Zqc*tmpArray1(1,3) + alphaZpq*TmpArray1(1,4)
     tmpArray2(2,4) = Xqc*tmpArray1(1,4) + alphaXpq*TmpArray1(1,5)
     tmpArray2(3,4) = Yqc*tmpArray1(1,4) + alphaYpq*TmpArray1(1,5)
     tmpArray2(4,4) = Zqc*tmpArray1(1,4) + alphaZpq*TmpArray1(1,5)
     TwoTerms(1) = inv2expQ*(TMPAuxArray(1) + alphaQ*TmpArray1(1,2))
     TMPAuxArray(5) = Xqc*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     TMPAuxArray(6) = Xqc*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     TMPAuxArray(7) = Xqc*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     TMPAuxArray(8) = Yqc*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     TMPAuxArray(9) = Yqc*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     TMPAuxArray(10) = Zqc*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(4) + alphaQ*TmpArray2(4,2))
     TMPAuxArray(11) = Xqc*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Yqc*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zqc*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xqc*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xqc*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xqc*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Yqc*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zqc*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Yqc*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zqc*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(5) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(8) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(10) + alphaQ*TmpArray3(10,2))
     do iTUV = 1,   20
      AuxArray(iP,iTUV) = TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,21) = Xqc*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     AuxArray(iP,22) = Yqc*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     AuxArray(iP,23) = Zqc*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     AuxArray(iP,24) = Xqc*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     AuxArray(iP,25) = Yqc*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     AuxArray(iP,26) = Xqc*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     AuxArray(iP,27) = Xqc*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     AuxArray(iP,28) = Xqc*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     AuxArray(iP,29) = Xqc*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     AuxArray(iP,30) = Xqc*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     AuxArray(iP,31) = Yqc*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     AuxArray(iP,32) = Zqc*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     AuxArray(iP,33) = Yqc*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     AuxArray(iP,34) = Yqc*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     AuxArray(iP,35) = Zqc*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
  ENDDO !iP = 1,nPassP
 end subroutine

subroutine VerticalRecurrenceGPUSeg1Prim5C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPassP,0: 5)
  real(realk),intent(in) :: reducedExponents(1),Qexp(1)
  real(realk),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3),integralPrefactor(1),QpreExpFac(1),PpreExpFac(nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(nPassP,   56)
  !local variables
  integer :: iPassP,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: TMPAUXarray(   35)
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
!$ACC         TMPAUXarray,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TMParray4,&
!$ACC         TMParray5,&
!$ACC         TwoTerms,&
!$ACC         iP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Ccenter, &
!$ACC        nPrimP,nPrimQ,nPassP)
  DO iP = 1,nPassP
   iPassP = iP
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     Xqc = Qcent(1) - Ccenter(1)
     Yqc = Qcent(2) - Ccenter(2)
     Zqc = Qcent(3) - Ccenter(3)
     alphaQ = -reducedExponents(1)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(1)*QpreExpFac(1)*PpreExpFac(iAtomA,iAtomB)
     TMPAuxarray(1) = PREF*RJ000Array(iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPassP, 3)
     TMParray1(1, 5) = PREF*RJ000Array(iPassP, 4)
     TMParray1(1, 6) = PREF*RJ000Array(iPassP, 5)
     TMPAuxArray(2) = Xqc*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqc*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqc*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(1) + alphaQ*TmpArray1(1,2))
     TMPAuxArray(5) = Xqc*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     TMPAuxArray(6) = Xqc*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     TMPAuxArray(7) = Xqc*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     TMPAuxArray(8) = Yqc*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     TMPAuxArray(9) = Yqc*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     TMPAuxArray(10) = Zqc*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(4) + alphaQ*TmpArray2(4,2))
     TMPAuxArray(11) = Xqc*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Yqc*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zqc*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xqc*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xqc*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xqc*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Yqc*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zqc*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Yqc*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zqc*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(5) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(8) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(10) + alphaQ*TmpArray3(10,2))
     TMPAuxArray(21) = Xqc*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Yqc*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zqc*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xqc*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Yqc*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xqc*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xqc*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xqc*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xqc*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xqc*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Yqc*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zqc*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Yqc*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Yqc*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zqc*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(11) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(14) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(16) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(17) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(19) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(20) + alphaQ*TmpArray4(20,2))
     do iTUV = 1,   35
      AuxArray(iP,iTUV) = TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,36) = Xqc*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     AuxArray(iP,37) = Yqc*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     AuxArray(iP,38) = Zqc*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     AuxArray(iP,39) = Xqc*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     AuxArray(iP,40) = Yqc*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     AuxArray(iP,41) = Xqc*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     AuxArray(iP,42) = Xqc*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     AuxArray(iP,43) = Zqc*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     AuxArray(iP,44) = Yqc*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     AuxArray(iP,45) = Xqc*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     AuxArray(iP,46) = Xqc*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     AuxArray(iP,47) = Xqc*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     AuxArray(iP,48) = Xqc*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     AuxArray(iP,49) = Xqc*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     AuxArray(iP,50) = Xqc*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     AuxArray(iP,51) = Yqc*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     AuxArray(iP,52) = Zqc*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     AuxArray(iP,53) = Yqc*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     AuxArray(iP,54) = Yqc*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     AuxArray(iP,55) = Yqc*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     AuxArray(iP,56) = Zqc*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
  ENDDO !iP = 1,nPassP
 end subroutine

subroutine VerticalRecurrenceGPUSeg1Prim6C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPassP,0: 6)
  real(realk),intent(in) :: reducedExponents(1),Qexp(1)
  real(realk),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3),integralPrefactor(1),QpreExpFac(1),PpreExpFac(nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(nPassP,   84)
  !local variables
  integer :: iPassP,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: TMPAUXarray(   56)
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
!$ACC         TMPAUXarray,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TMParray4,&
!$ACC         TMParray5,&
!$ACC         TMParray6,&
!$ACC         TwoTerms,&
!$ACC         iP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Ccenter, &
!$ACC        nPrimP,nPrimQ,nPassP)
  DO iP = 1,nPassP
   iPassP = iP
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     Xqc = Qcent(1) - Ccenter(1)
     Yqc = Qcent(2) - Ccenter(2)
     Zqc = Qcent(3) - Ccenter(3)
     alphaQ = -reducedExponents(1)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(1)*QpreExpFac(1)*PpreExpFac(iAtomA,iAtomB)
     TMPAuxarray(1) = PREF*RJ000Array(iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPassP, 3)
     TMParray1(1, 5) = PREF*RJ000Array(iPassP, 4)
     TMParray1(1, 6) = PREF*RJ000Array(iPassP, 5)
     TMParray1(1, 7) = PREF*RJ000Array(iPassP, 6)
     TMPAuxArray(2) = Xqc*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqc*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqc*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(1) + alphaQ*TmpArray1(1,2))
     TMPAuxArray(5) = Xqc*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     TMPAuxArray(6) = Xqc*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     TMPAuxArray(7) = Xqc*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     TMPAuxArray(8) = Yqc*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     TMPAuxArray(9) = Yqc*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     TMPAuxArray(10) = Zqc*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(4) + alphaQ*TmpArray2(4,2))
     TMPAuxArray(11) = Xqc*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Yqc*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zqc*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xqc*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xqc*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xqc*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Yqc*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zqc*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Yqc*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zqc*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(5) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(8) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(10) + alphaQ*TmpArray3(10,2))
     TMPAuxArray(21) = Xqc*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Yqc*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zqc*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xqc*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Yqc*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xqc*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xqc*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xqc*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xqc*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xqc*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Yqc*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zqc*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Yqc*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Yqc*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zqc*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(11) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(14) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(16) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(17) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(19) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(20) + alphaQ*TmpArray4(20,2))
     TMPAuxArray(36) = Xqc*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     TMPAuxArray(37) = Yqc*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     TMPAuxArray(38) = Zqc*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     TMPAuxArray(39) = Xqc*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     TMPAuxArray(40) = Yqc*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     TMPAuxArray(41) = Xqc*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     TMPAuxArray(42) = Xqc*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     TMPAuxArray(43) = Zqc*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     TMPAuxArray(44) = Yqc*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     TMPAuxArray(45) = Xqc*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     TMPAuxArray(46) = Xqc*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     TMPAuxArray(47) = Xqc*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     TMPAuxArray(48) = Xqc*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     TMPAuxArray(49) = Xqc*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     TMPAuxArray(50) = Xqc*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     TMPAuxArray(51) = Yqc*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     TMPAuxArray(52) = Zqc*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     TMPAuxArray(53) = Yqc*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     TMPAuxArray(54) = Yqc*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     TMPAuxArray(55) = Yqc*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     TMPAuxArray(56) = Zqc*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
      AuxArray(iP,iTUV) = TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,57) = Xqc*TMPAuxArray(36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     AuxArray(iP,58) = Yqc*TMPAuxArray(36) + alphaYpq*TmpArray6(36,2)
     AuxArray(iP,59) = Zqc*TMPAuxArray(36) + alphaZpq*TmpArray6(36,2)
     AuxArray(iP,60) = Xqc*TMPAuxArray(39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     AuxArray(iP,61) = Yqc*TMPAuxArray(38) + alphaYpq*TmpArray6(38,2)
     AuxArray(iP,62) = Xqc*TMPAuxArray(41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     AuxArray(iP,63) = Xqc*TMPAuxArray(42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     AuxArray(iP,64) = Zqc*TMPAuxArray(39) + alphaZpq*TmpArray6(39,2)
     AuxArray(iP,65) = Yqc*TMPAuxArray(41) + alphaYpq*TmpArray6(41,2)
     AuxArray(iP,66) = Xqc*TMPAuxArray(45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     AuxArray(iP,67) = Xqc*TMPAuxArray(46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     AuxArray(iP,68) = Zqc*TMPAuxArray(42) + alphaZpq*TmpArray6(42,2)
     AuxArray(iP,69) = Xqc*TMPAuxArray(48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     AuxArray(iP,70) = Yqc*TMPAuxArray(45) + alphaYpq*TmpArray6(45,2)
     AuxArray(iP,71) = Xqc*TMPAuxArray(50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     AuxArray(iP,72) = Xqc*TMPAuxArray(51) + alphaXpq*TmpArray6(51,2)
     AuxArray(iP,73) = Xqc*TMPAuxArray(52) + alphaXpq*TmpArray6(52,2)
     AuxArray(iP,74) = Xqc*TMPAuxArray(53) + alphaXpq*TmpArray6(53,2)
     AuxArray(iP,75) = Xqc*TMPAuxArray(54) + alphaXpq*TmpArray6(54,2)
     AuxArray(iP,76) = Xqc*TMPAuxArray(55) + alphaXpq*TmpArray6(55,2)
     AuxArray(iP,77) = Xqc*TMPAuxArray(56) + alphaXpq*TmpArray6(56,2)
     AuxArray(iP,78) = Yqc*TMPAuxArray(51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     AuxArray(iP,79) = Zqc*TMPAuxArray(51) + alphaZpq*TmpArray6(51,2)
     AuxArray(iP,80) = Yqc*TMPAuxArray(53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     AuxArray(iP,81) = Yqc*TMPAuxArray(54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     AuxArray(iP,82) = Yqc*TMPAuxArray(55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     AuxArray(iP,83) = Yqc*TMPAuxArray(56) + alphaYpq*TmpArray6(56,2)
     AuxArray(iP,84) = Zqc*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
  ENDDO !iP = 1,nPassP
 end subroutine

subroutine VerticalRecurrenceGPUSeg1Prim7C(nPassP,nPrimP,nPrimQ,&
         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&
         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: RJ000Array(nPassP,0: 7)
  real(realk),intent(in) :: reducedExponents(1),Qexp(1)
  real(realk),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3),integralPrefactor(1),QpreExpFac(1),PpreExpFac(nAtomsA,nAtomsB)
  real(realk),intent(in) :: Ccenter(3)
  real(realk),intent(inout) :: AUXarray(nPassP,  120)
  !local variables
  integer :: iPassP,ipnt,IP,iTUV,iAtomA,iAtomB
  real(realk) :: TMPAUXarray(   84)
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
!$ACC         TMPAUXarray,&
!$ACC         TMParray1,&
!$ACC         TMParray2,&
!$ACC         TMParray3,&
!$ACC         TMParray4,&
!$ACC         TMParray5,&
!$ACC         TMParray6,&
!$ACC         TMParray7,&
!$ACC         TwoTerms,&
!$ACC         iP,iPassP) &
!$ACC PRESENT(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&
!$ACC        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&
!$ACC        Qexp,Ccenter, &
!$ACC        nPrimP,nPrimQ,nPassP)
  DO iP = 1,nPassP
   iPassP = iP
     iAtomA = iAtomApass(iPassP)
     iAtomB = iAtomBpass(iPassP)
     invexpQ = D1/Qexp(1)
     inv2expQ = D05*invexpQ
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     Xqc = Qcent(1) - Ccenter(1)
     Yqc = Qcent(2) - Ccenter(2)
     Zqc = Qcent(3) - Ccenter(3)
     alphaQ = -reducedExponents(1)*invexpQ
     alphaXpq = alphaQ*Xpq
     alphaYpq = alphaQ*Ypq
     alphaZpq = alphaQ*Zpq
     PREF = integralPrefactor(1)*QpreExpFac(1)*PpreExpFac(iAtomA,iAtomB)
     TMPAuxarray(1) = PREF*RJ000Array(iPassP,0)
     TMParray1(1, 2) = PREF*RJ000Array(iPassP, 1)
     TMParray1(1, 3) = PREF*RJ000Array(iPassP, 2)
     TMParray1(1, 4) = PREF*RJ000Array(iPassP, 3)
     TMParray1(1, 5) = PREF*RJ000Array(iPassP, 4)
     TMParray1(1, 6) = PREF*RJ000Array(iPassP, 5)
     TMParray1(1, 7) = PREF*RJ000Array(iPassP, 6)
     TMParray1(1, 8) = PREF*RJ000Array(iPassP, 7)
     TMPAuxArray(2) = Xqc*TMPAuxArray(1) + alphaXpq*TmpArray1(1,2)
     TMPAuxArray(3) = Yqc*TMPAuxArray(1) + alphaYpq*TmpArray1(1,2)
     TMPAuxArray(4) = Zqc*TMPAuxArray(1) + alphaZpq*TmpArray1(1,2)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(1) + alphaQ*TmpArray1(1,2))
     TMPAuxArray(5) = Xqc*TMPAuxArray(2) + alphaXpq*TmpArray2(2,2) + TwoTerms(1)
     TMPAuxArray(6) = Xqc*TMPAuxArray(3) + alphaXpq*TmpArray2(3,2)
     TMPAuxArray(7) = Xqc*TMPAuxArray(4) + alphaXpq*TmpArray2(4,2)
     TMPAuxArray(8) = Yqc*TMPAuxArray(3) + alphaYpq*TmpArray2(3,2) + TwoTerms(1)
     TMPAuxArray(9) = Yqc*TMPAuxArray(4) + alphaYpq*TmpArray2(4,2)
     TMPAuxArray(10) = Zqc*TMPAuxArray(4) + alphaZpq*TmpArray2(4,2) + TwoTerms(1)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(2) + alphaQ*TmpArray2(2,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(3) + alphaQ*TmpArray2(3,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(4) + alphaQ*TmpArray2(4,2))
     TMPAuxArray(11) = Xqc*TMPAuxArray(5) + alphaXpq*TmpArray3(5,2) + 2*TwoTerms(1)
     TMPAuxArray(12) = Yqc*TMPAuxArray(5) + alphaYpq*TmpArray3(5,2)
     TMPAuxArray(13) = Zqc*TMPAuxArray(5) + alphaZpq*TmpArray3(5,2)
     TMPAuxArray(14) = Xqc*TMPAuxArray(8) + alphaXpq*TmpArray3(8,2)
     TMPAuxArray(15) = Xqc*TMPAuxArray(9) + alphaXpq*TmpArray3(9,2)
     TMPAuxArray(16) = Xqc*TMPAuxArray(10) + alphaXpq*TmpArray3(10,2)
     TMPAuxArray(17) = Yqc*TMPAuxArray(8) + alphaYpq*TmpArray3(8,2) + 2*TwoTerms(2)
     TMPAuxArray(18) = Zqc*TMPAuxArray(8) + alphaZpq*TmpArray3(8,2)
     TMPAuxArray(19) = Yqc*TMPAuxArray(10) + alphaYpq*TmpArray3(10,2)
     TMPAuxArray(20) = Zqc*TMPAuxArray(10) + alphaZpq*TmpArray3(10,2) + 2*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(5) + alphaQ*TmpArray3(5,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(8) + alphaQ*TmpArray3(8,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(10) + alphaQ*TmpArray3(10,2))
     TMPAuxArray(21) = Xqc*TMPAuxArray(11) + alphaXpq*TmpArray4(11,2) + 3*TwoTerms(1)
     TMPAuxArray(22) = Yqc*TMPAuxArray(11) + alphaYpq*TmpArray4(11,2)
     TMPAuxArray(23) = Zqc*TMPAuxArray(11) + alphaZpq*TmpArray4(11,2)
     TMPAuxArray(24) = Xqc*TMPAuxArray(14) + alphaXpq*TmpArray4(14,2) + TwoTerms(2)
     TMPAuxArray(25) = Yqc*TMPAuxArray(13) + alphaYpq*TmpArray4(13,2)
     TMPAuxArray(26) = Xqc*TMPAuxArray(16) + alphaXpq*TmpArray4(16,2) + TwoTerms(3)
     TMPAuxArray(27) = Xqc*TMPAuxArray(17) + alphaXpq*TmpArray4(17,2)
     TMPAuxArray(28) = Xqc*TMPAuxArray(18) + alphaXpq*TmpArray4(18,2)
     TMPAuxArray(29) = Xqc*TMPAuxArray(19) + alphaXpq*TmpArray4(19,2)
     TMPAuxArray(30) = Xqc*TMPAuxArray(20) + alphaXpq*TmpArray4(20,2)
     TMPAuxArray(31) = Yqc*TMPAuxArray(17) + alphaYpq*TmpArray4(17,2) + 3*TwoTerms(2)
     TMPAuxArray(32) = Zqc*TMPAuxArray(17) + alphaZpq*TmpArray4(17,2)
     TMPAuxArray(33) = Yqc*TMPAuxArray(19) + alphaYpq*TmpArray4(19,2) + TwoTerms(3)
     TMPAuxArray(34) = Yqc*TMPAuxArray(20) + alphaYpq*TmpArray4(20,2)
     TMPAuxArray(35) = Zqc*TMPAuxArray(20) + alphaZpq*TmpArray4(20,2) + 3*TwoTerms(3)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(11) + alphaQ*TmpArray4(11,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(14) + alphaQ*TmpArray4(14,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(16) + alphaQ*TmpArray4(16,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(17) + alphaQ*TmpArray4(17,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(19) + alphaQ*TmpArray4(19,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(20) + alphaQ*TmpArray4(20,2))
     TMPAuxArray(36) = Xqc*TMPAuxArray(21) + alphaXpq*TmpArray5(21,2) + 4*TwoTerms(1)
     TMPAuxArray(37) = Yqc*TMPAuxArray(21) + alphaYpq*TmpArray5(21,2)
     TMPAuxArray(38) = Zqc*TMPAuxArray(21) + alphaZpq*TmpArray5(21,2)
     TMPAuxArray(39) = Xqc*TMPAuxArray(24) + alphaXpq*TmpArray5(24,2) + 2*TwoTerms(2)
     TMPAuxArray(40) = Yqc*TMPAuxArray(23) + alphaYpq*TmpArray5(23,2)
     TMPAuxArray(41) = Xqc*TMPAuxArray(26) + alphaXpq*TmpArray5(26,2) + 2*TwoTerms(3)
     TMPAuxArray(42) = Xqc*TMPAuxArray(27) + alphaXpq*TmpArray5(27,2) + TwoTerms(4)
     TMPAuxArray(43) = Zqc*TMPAuxArray(24) + alphaZpq*TmpArray5(24,2)
     TMPAuxArray(44) = Yqc*TMPAuxArray(26) + alphaYpq*TmpArray5(26,2)
     TMPAuxArray(45) = Xqc*TMPAuxArray(30) + alphaXpq*TmpArray5(30,2) + TwoTerms(6)
     TMPAuxArray(46) = Xqc*TMPAuxArray(31) + alphaXpq*TmpArray5(31,2)
     TMPAuxArray(47) = Xqc*TMPAuxArray(32) + alphaXpq*TmpArray5(32,2)
     TMPAuxArray(48) = Xqc*TMPAuxArray(33) + alphaXpq*TmpArray5(33,2)
     TMPAuxArray(49) = Xqc*TMPAuxArray(34) + alphaXpq*TmpArray5(34,2)
     TMPAuxArray(50) = Xqc*TMPAuxArray(35) + alphaXpq*TmpArray5(35,2)
     TMPAuxArray(51) = Yqc*TMPAuxArray(31) + alphaYpq*TmpArray5(31,2) + 4*TwoTerms(4)
     TMPAuxArray(52) = Zqc*TMPAuxArray(31) + alphaZpq*TmpArray5(31,2)
     TMPAuxArray(53) = Yqc*TMPAuxArray(33) + alphaYpq*TmpArray5(33,2) + 2*TwoTerms(5)
     TMPAuxArray(54) = Yqc*TMPAuxArray(34) + alphaYpq*TmpArray5(34,2) + TwoTerms(6)
     TMPAuxArray(55) = Yqc*TMPAuxArray(35) + alphaYpq*TmpArray5(35,2)
     TMPAuxArray(56) = Zqc*TMPAuxArray(35) + alphaZpq*TmpArray5(35,2) + 4*TwoTerms(6)
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
     TwoTerms(1) = inv2expQ*(TMPAuxArray(21) + alphaQ*TmpArray5(21,2))
     TwoTerms(2) = inv2expQ*(TMPAuxArray(24) + alphaQ*TmpArray5(24,2))
     TwoTerms(3) = inv2expQ*(TMPAuxArray(26) + alphaQ*TmpArray5(26,2))
     TwoTerms(4) = inv2expQ*(TMPAuxArray(27) + alphaQ*TmpArray5(27,2))
     TwoTerms(5) = inv2expQ*(TMPAuxArray(30) + alphaQ*TmpArray5(30,2))
     TwoTerms(6) = inv2expQ*(TMPAuxArray(31) + alphaQ*TmpArray5(31,2))
     TwoTerms(7) = inv2expQ*(TMPAuxArray(33) + alphaQ*TmpArray5(33,2))
     TwoTerms(8) = inv2expQ*(TMPAuxArray(34) + alphaQ*TmpArray5(34,2))
     TwoTerms(9) = inv2expQ*(TMPAuxArray(35) + alphaQ*TmpArray5(35,2))
     TMPAuxArray(57) = Xqc*TMPAuxArray(36) + alphaXpq*TmpArray6(36,2) + 5*TwoTerms(1)
     TMPAuxArray(58) = Yqc*TMPAuxArray(36) + alphaYpq*TmpArray6(36,2)
     TMPAuxArray(59) = Zqc*TMPAuxArray(36) + alphaZpq*TmpArray6(36,2)
     TMPAuxArray(60) = Xqc*TMPAuxArray(39) + alphaXpq*TmpArray6(39,2) + 3*TwoTerms(2)
     TMPAuxArray(61) = Yqc*TMPAuxArray(38) + alphaYpq*TmpArray6(38,2)
     TMPAuxArray(62) = Xqc*TMPAuxArray(41) + alphaXpq*TmpArray6(41,2) + 3*TwoTerms(3)
     TMPAuxArray(63) = Xqc*TMPAuxArray(42) + alphaXpq*TmpArray6(42,2) + 2*TwoTerms(4)
     TMPAuxArray(64) = Zqc*TMPAuxArray(39) + alphaZpq*TmpArray6(39,2)
     TMPAuxArray(65) = Yqc*TMPAuxArray(41) + alphaYpq*TmpArray6(41,2)
     TMPAuxArray(66) = Xqc*TMPAuxArray(45) + alphaXpq*TmpArray6(45,2) + 2*TwoTerms(5)
     TMPAuxArray(67) = Xqc*TMPAuxArray(46) + alphaXpq*TmpArray6(46,2) + TwoTerms(6)
     TMPAuxArray(68) = Zqc*TMPAuxArray(42) + alphaZpq*TmpArray6(42,2)
     TMPAuxArray(69) = Xqc*TMPAuxArray(48) + alphaXpq*TmpArray6(48,2) + TwoTerms(7)
     TMPAuxArray(70) = Yqc*TMPAuxArray(45) + alphaYpq*TmpArray6(45,2)
     TMPAuxArray(71) = Xqc*TMPAuxArray(50) + alphaXpq*TmpArray6(50,2) + TwoTerms(9)
     TMPAuxArray(72) = Xqc*TMPAuxArray(51) + alphaXpq*TmpArray6(51,2)
     TMPAuxArray(73) = Xqc*TMPAuxArray(52) + alphaXpq*TmpArray6(52,2)
     TMPAuxArray(74) = Xqc*TMPAuxArray(53) + alphaXpq*TmpArray6(53,2)
     TMPAuxArray(75) = Xqc*TMPAuxArray(54) + alphaXpq*TmpArray6(54,2)
     TMPAuxArray(76) = Xqc*TMPAuxArray(55) + alphaXpq*TmpArray6(55,2)
     TMPAuxArray(77) = Xqc*TMPAuxArray(56) + alphaXpq*TmpArray6(56,2)
     TMPAuxArray(78) = Yqc*TMPAuxArray(51) + alphaYpq*TmpArray6(51,2) + 5*TwoTerms(6)
     TMPAuxArray(79) = Zqc*TMPAuxArray(51) + alphaZpq*TmpArray6(51,2)
     TMPAuxArray(80) = Yqc*TMPAuxArray(53) + alphaYpq*TmpArray6(53,2) + 3*TwoTerms(7)
     TMPAuxArray(81) = Yqc*TMPAuxArray(54) + alphaYpq*TmpArray6(54,2) + 2*TwoTerms(8)
     TMPAuxArray(82) = Yqc*TMPAuxArray(55) + alphaYpq*TmpArray6(55,2) + TwoTerms(9)
     TMPAuxArray(83) = Yqc*TMPAuxArray(56) + alphaYpq*TmpArray6(56,2)
     TMPAuxArray(84) = Zqc*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
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
      AuxArray(iP,iTUV) = TMPAuxarray(iTUV)
     enddo
     AuxArray(iP,85) = Xqc*TMPAuxArray(57) + alphaXpq*TmpArray7(57,2) + 6*TwoTerms(1)
     AuxArray(iP,86) = Yqc*TMPAuxArray(57) + alphaYpq*TmpArray7(57,2)
     AuxArray(iP,87) = Zqc*TMPAuxArray(57) + alphaZpq*TmpArray7(57,2)
     AuxArray(iP,88) = Xqc*TMPAuxArray(60) + alphaXpq*TmpArray7(60,2) + 4*TwoTerms(2)
     AuxArray(iP,89) = Yqc*TMPAuxArray(59) + alphaYpq*TmpArray7(59,2)
     AuxArray(iP,90) = Xqc*TMPAuxArray(62) + alphaXpq*TmpArray7(62,2) + 4*TwoTerms(3)
     AuxArray(iP,91) = Xqc*TMPAuxArray(63) + alphaXpq*TmpArray7(63,2) + 3*TwoTerms(4)
     AuxArray(iP,92) = Zqc*TMPAuxArray(60) + alphaZpq*TmpArray7(60,2)
     AuxArray(iP,93) = Yqc*TMPAuxArray(62) + alphaYpq*TmpArray7(62,2)
     AuxArray(iP,94) = Xqc*TMPAuxArray(66) + alphaXpq*TmpArray7(66,2) + 3*TwoTerms(5)
     AuxArray(iP,95) = Xqc*TMPAuxArray(67) + alphaXpq*TmpArray7(67,2) + 2*TwoTerms(6)
     AuxArray(iP,96) = Zqc*TMPAuxArray(63) + alphaZpq*TmpArray7(63,2)
     AuxArray(iP,97) = Xqc*TMPAuxArray(69) + alphaXpq*TmpArray7(69,2) + 2*TwoTerms(7)
     AuxArray(iP,98) = Yqc*TMPAuxArray(66) + alphaYpq*TmpArray7(66,2)
     AuxArray(iP,99) = Xqc*TMPAuxArray(71) + alphaXpq*TmpArray7(71,2) + 2*TwoTerms(8)
     AuxArray(iP,100) = Xqc*TMPAuxArray(72) + alphaXpq*TmpArray7(72,2) + TwoTerms(9)
     AuxArray(iP,101) = Zqc*TMPAuxArray(67) + alphaZpq*TmpArray7(67,2)
     AuxArray(iP,102) = Xqc*TMPAuxArray(74) + alphaXpq*TmpArray7(74,2) + TwoTerms(10)
     AuxArray(iP,103) = Xqc*TMPAuxArray(75) + alphaXpq*TmpArray7(75,2) + TwoTerms(11)
     AuxArray(iP,104) = Yqc*TMPAuxArray(71) + alphaYpq*TmpArray7(71,2)
     AuxArray(iP,105) = Xqc*TMPAuxArray(77) + alphaXpq*TmpArray7(77,2) + TwoTerms(13)
     AuxArray(iP,106) = Xqc*TMPAuxArray(78) + alphaXpq*TmpArray7(78,2)
     AuxArray(iP,107) = Xqc*TMPAuxArray(79) + alphaXpq*TmpArray7(79,2)
     AuxArray(iP,108) = Xqc*TMPAuxArray(80) + alphaXpq*TmpArray7(80,2)
     AuxArray(iP,109) = Xqc*TMPAuxArray(81) + alphaXpq*TmpArray7(81,2)
     AuxArray(iP,110) = Xqc*TMPAuxArray(82) + alphaXpq*TmpArray7(82,2)
     AuxArray(iP,111) = Xqc*TMPAuxArray(83) + alphaXpq*TmpArray7(83,2)
     AuxArray(iP,112) = Xqc*TMPAuxArray(84) + alphaXpq*TmpArray7(84,2)
     AuxArray(iP,113) = Yqc*TMPAuxArray(78) + alphaYpq*TmpArray7(78,2) + 6*TwoTerms(9)
     AuxArray(iP,114) = Zqc*TMPAuxArray(78) + alphaZpq*TmpArray7(78,2)
     AuxArray(iP,115) = Yqc*TMPAuxArray(80) + alphaYpq*TmpArray7(80,2) + 4*TwoTerms(10)
     AuxArray(iP,116) = Yqc*TMPAuxArray(81) + alphaYpq*TmpArray7(81,2) + 3*TwoTerms(11)
     AuxArray(iP,117) = Yqc*TMPAuxArray(82) + alphaYpq*TmpArray7(82,2) + 2*TwoTerms(12)
     AuxArray(iP,118) = Yqc*TMPAuxArray(83) + alphaYpq*TmpArray7(83,2) + TwoTerms(13)
     AuxArray(iP,119) = Yqc*TMPAuxArray(84) + alphaYpq*TmpArray7(84,2)
     AuxArray(iP,120) = Zqc*TMPAuxArray(84) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
  ENDDO !iP = 1,nPassP
 end subroutine
end module
