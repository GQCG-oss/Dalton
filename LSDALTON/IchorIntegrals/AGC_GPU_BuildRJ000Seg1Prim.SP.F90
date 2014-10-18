module SPAGC_GPU_OBS_BUILDRJ000MODSeg1Prim
 use IchorPrecisionMod
  
 CONTAINS

subroutine SPBuildRJ000GPUSeg1Prim2(nPassP,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: TABFJW(0: 5,0:1200)
  real(reals),intent(in) :: reducedExponents(1)
  real(reals),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3)
  real(reals),intent(inout) :: RJ000array(nPassP,0: 2)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iP,iPassP,ipnt,iAtomA,iAtomB
  real(reals) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(reals) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(reals) :: RJ000(0: 2)
  real(reals),PARAMETER :: TENTH = 0.01E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D2=2.0E0_reals
  real(reals),PARAMETER :: D2JP36=  4.0000000000000000E+01_reals
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
!$ACC parallel loop &
!$ACC present(nPassP,nPrimP,nPrimQ,IatomApass,IatomBpass,&
!$ACC         TABFJW,reducedExponents,Pcent,Qcent,RJ000array)       &
!$ACC private(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         iP,iPassP,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         mPX,mPY,mPZ,RWVAL,GVAL) ASYNC(iASync)
  DO iP = 1,nPassP
   iPassP = iP
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(1)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000Array(iP, 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array(iP, 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array(iP, 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000( 1) = RWVAL*(( 1 - D05)*RJ000( 0)-REXPW)
      RJ000( 2) = RWVAL*(( 2 - D05)*RJ000( 1)-REXPW)
      RJ000Array(iP,0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000Array(iP, 0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
     ENDIF
  ENDDO
 end subroutine

subroutine SPBuildRJ000GPUSeg1Prim3(nPassP,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: TABFJW(0: 6,0:1200)
  real(reals),intent(in) :: reducedExponents(1)
  real(reals),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3)
  real(reals),intent(inout) :: RJ000array(nPassP,0: 3)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iP,iPassP,ipnt,iAtomA,iAtomB
  real(reals) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(reals) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(reals) :: RJ000(0: 3)
  real(reals),PARAMETER :: TENTH = 0.01E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D2=2.0E0_reals
  real(reals),PARAMETER :: D2JP36=  4.2000000000000000E+01_reals
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
!$ACC parallel loop &
!$ACC present(nPassP,nPrimP,nPrimQ,IatomApass,IatomBpass,&
!$ACC         TABFJW,reducedExponents,Pcent,Qcent,RJ000array)       &
!$ACC private(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         iP,iPassP,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         mPX,mPY,mPZ,RWVAL,GVAL) ASYNC(iASync)
  DO iP = 1,nPassP
   iPassP = iP
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(1)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000Array(iP, 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array(iP, 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array(iP, 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array(iP, 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
     !  12 < WVAL <= (2J+36) 
     ELSE IF (WVAL.LE.D2JP36) THEN
      REXPW = D05*EXP(-WVAL)
      RWVAL = D1/WVAL
      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
      RJ000( 1) = RWVAL*(( 1 - D05)*RJ000( 0)-REXPW)
      RJ000( 2) = RWVAL*(( 2 - D05)*RJ000( 1)-REXPW)
      RJ000( 3) = RWVAL*(( 3 - D05)*RJ000( 2)-REXPW)
      RJ000Array(iP,0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000( 3) = RWVAL*( 3 - D05)*RJ000( 2)
      RJ000Array(iP, 0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
     ENDIF
  ENDDO
 end subroutine

subroutine SPBuildRJ000GPUSeg1Prim4(nPassP,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: TABFJW(0: 7,0:1200)
  real(reals),intent(in) :: reducedExponents(1)
  real(reals),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3)
  real(reals),intent(inout) :: RJ000array(nPassP,0: 4)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iP,iPassP,ipnt,iAtomA,iAtomB
  real(reals) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(reals) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(reals) :: RJ000(0: 4)
  real(reals),PARAMETER :: TENTH = 0.01E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D2=2.0E0_reals
  real(reals),PARAMETER :: D2JP36=  4.4000000000000000E+01_reals
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
!$ACC parallel loop &
!$ACC present(nPassP,nPrimP,nPrimQ,IatomApass,IatomBpass,&
!$ACC         TABFJW,reducedExponents,Pcent,Qcent,RJ000array)       &
!$ACC private(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         iP,iPassP,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         mPX,mPY,mPZ,RWVAL,GVAL) ASYNC(iASync)
  DO iP = 1,nPassP
   iPassP = iP
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(1)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000Array(iP, 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array(iP, 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array(iP, 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array(iP, 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000Array(iP, 4) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
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
      RJ000Array(iP,0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
      RJ000Array(iP, 4) = RJ000( 4)
     !  (2J+36) < WVAL 
     ELSE
      RWVAL = PID4/WVAL
      RJ000(0) = SQRT(RWVAL)
      RWVAL = RWVAL*PID4I
      RJ000( 1) = RWVAL*( 1 - D05)*RJ000( 0)
      RJ000( 2) = RWVAL*( 2 - D05)*RJ000( 1)
      RJ000( 3) = RWVAL*( 3 - D05)*RJ000( 2)
      RJ000( 4) = RWVAL*( 4 - D05)*RJ000( 3)
      RJ000Array(iP, 0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
      RJ000Array(iP, 4) = RJ000( 4)
     ENDIF
  ENDDO
 end subroutine

subroutine SPBuildRJ000GPUSeg1Prim5(nPassP,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: TABFJW(0: 8,0:1200)
  real(reals),intent(in) :: reducedExponents(1)
  real(reals),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3)
  real(reals),intent(inout) :: RJ000array(nPassP,0: 5)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iP,iPassP,ipnt,iAtomA,iAtomB
  real(reals) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(reals) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(reals) :: RJ000(0: 5)
  real(reals),PARAMETER :: TENTH = 0.01E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D2=2.0E0_reals
  real(reals),PARAMETER :: D2JP36=  4.6000000000000000E+01_reals
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
!$ACC parallel loop &
!$ACC present(nPassP,nPrimP,nPrimQ,IatomApass,IatomBpass,&
!$ACC         TABFJW,reducedExponents,Pcent,Qcent,RJ000array)       &
!$ACC private(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         iP,iPassP,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         mPX,mPY,mPZ,RWVAL,GVAL) ASYNC(iASync)
  DO iP = 1,nPassP
   iPassP = iP
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(1)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000Array(iP, 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array(iP, 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array(iP, 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array(iP, 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000Array(iP, 4) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000Array(iP, 5) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
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
      RJ000Array(iP,0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
      RJ000Array(iP, 4) = RJ000( 4)
      RJ000Array(iP, 5) = RJ000( 5)
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
      RJ000Array(iP, 0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
      RJ000Array(iP, 4) = RJ000( 4)
      RJ000Array(iP, 5) = RJ000( 5)
     ENDIF
  ENDDO
 end subroutine

subroutine SPBuildRJ000GPUSeg1Prim6(nPassP,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: TABFJW(0: 9,0:1200)
  real(reals),intent(in) :: reducedExponents(1)
  real(reals),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3)
  real(reals),intent(inout) :: RJ000array(nPassP,0: 6)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iP,iPassP,ipnt,iAtomA,iAtomB
  real(reals) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(reals) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(reals) :: RJ000(0: 6)
  real(reals),PARAMETER :: TENTH = 0.01E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D2=2.0E0_reals
  real(reals),PARAMETER :: D2JP36=  4.8000000000000000E+01_reals
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
!$ACC parallel loop &
!$ACC present(nPassP,nPrimP,nPrimQ,IatomApass,IatomBpass,&
!$ACC         TABFJW,reducedExponents,Pcent,Qcent,RJ000array)       &
!$ACC private(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         iP,iPassP,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         mPX,mPY,mPZ,RWVAL,GVAL) ASYNC(iASync)
  DO iP = 1,nPassP
   iPassP = iP
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(1)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000Array(iP, 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array(iP, 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array(iP, 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array(iP, 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000Array(iP, 4) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000Array(iP, 5) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
      RJ000Array(iP, 6) = TABFJW( 6,IPNT)-TABFJW( 7,IPNT)*WDIFF+TABFJW( 8,IPNT)*W2+TABFJW( 9,IPNT)*W3
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
      RJ000Array(iP,0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
      RJ000Array(iP, 4) = RJ000( 4)
      RJ000Array(iP, 5) = RJ000( 5)
      RJ000Array(iP, 6) = RJ000( 6)
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
      RJ000Array(iP, 0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
      RJ000Array(iP, 4) = RJ000( 4)
      RJ000Array(iP, 5) = RJ000( 5)
      RJ000Array(iP, 6) = RJ000( 6)
     ENDIF
  ENDDO
 end subroutine

subroutine SPBuildRJ000GPUSeg1Prim7(nPassP,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: TABFJW(0:10,0:1200)
  real(reals),intent(in) :: reducedExponents(1)
  real(reals),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3)
  real(reals),intent(inout) :: RJ000array(nPassP,0: 7)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iP,iPassP,ipnt,iAtomA,iAtomB
  real(reals) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(reals) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(reals) :: RJ000(0: 7)
  real(reals),PARAMETER :: TENTH = 0.01E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D2=2.0E0_reals
  real(reals),PARAMETER :: D2JP36=  5.0000000000000000E+01_reals
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
!$ACC parallel loop &
!$ACC present(nPassP,nPrimP,nPrimQ,IatomApass,IatomBpass,&
!$ACC         TABFJW,reducedExponents,Pcent,Qcent,RJ000array)       &
!$ACC private(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         iP,iPassP,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         mPX,mPY,mPZ,RWVAL,GVAL) ASYNC(iASync)
  DO iP = 1,nPassP
   iPassP = iP
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(1)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000Array(iP, 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array(iP, 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array(iP, 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array(iP, 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000Array(iP, 4) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000Array(iP, 5) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
      RJ000Array(iP, 6) = TABFJW( 6,IPNT)-TABFJW( 7,IPNT)*WDIFF+TABFJW( 8,IPNT)*W2+TABFJW( 9,IPNT)*W3
      RJ000Array(iP, 7) = TABFJW( 7,IPNT)-TABFJW( 8,IPNT)*WDIFF+TABFJW( 9,IPNT)*W2+TABFJW(10,IPNT)*W3
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
      RJ000Array(iP,0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
      RJ000Array(iP, 4) = RJ000( 4)
      RJ000Array(iP, 5) = RJ000( 5)
      RJ000Array(iP, 6) = RJ000( 6)
      RJ000Array(iP, 7) = RJ000( 7)
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
      RJ000Array(iP, 0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
      RJ000Array(iP, 4) = RJ000( 4)
      RJ000Array(iP, 5) = RJ000( 5)
      RJ000Array(iP, 6) = RJ000( 6)
      RJ000Array(iP, 7) = RJ000( 7)
     ENDIF
  ENDDO
 end subroutine

subroutine SPBuildRJ000GPUSeg1Prim8(nPassP,nPrimP,nPrimQ,reducedExponents,&
         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&
         & RJ000array,iASync)
  implicit none
  integer,intent(in) :: nPassP,nPrimP,nPrimQ
  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(reals),intent(in) :: TABFJW(0:11,0:1200)
  real(reals),intent(in) :: reducedExponents(1)
  real(reals),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3)
  real(reals),intent(inout) :: RJ000array(nPassP,0: 8)
  integer(kind=acckind),intent(in) :: iASync
  !local variables
  integer :: iP,iPassP,ipnt,iAtomA,iAtomB
  real(reals) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq
  real(reals) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL
  real(reals) :: RJ000(0: 8)
  real(reals),PARAMETER :: TENTH = 0.01E0_reals,D05 =0.5E0_reals
  real(reals),parameter :: D2=2.0E0_reals
  real(reals),PARAMETER :: D2JP36=  5.2000000000000000E+01_reals
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
!$ACC parallel loop &
!$ACC present(nPassP,nPrimP,nPrimQ,IatomApass,IatomBpass,&
!$ACC         TABFJW,reducedExponents,Pcent,Qcent,RJ000array)       &
!$ACC private(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!$ACC         iP,iPassP,&
!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!$ACC         mPX,mPY,mPZ,RWVAL,GVAL) ASYNC(iASync)
  DO iP = 1,nPassP
   iPassP = iP
   iAtomA = iAtomApass(iPassP)
   iAtomB = iAtomBpass(iPassP)
     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)
     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)
     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)
     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq
     WVAL = reducedExponents(1)*squaredDistance
     !  0 < WVAL < 12 
     IF (WVAL .LT. D12) THEN
      IPNT = NINT(D100*WVAL)
      WDIFF = WVAL - TENTH*IPNT
      W2    = WDIFF*WDIFF
      W3    = W2*WDIFF
      W2    = W2*D05
      W3    = W3*COEF3
      RJ000Array(iP, 0) = TABFJW( 0,IPNT)-TABFJW( 1,IPNT)*WDIFF+TABFJW( 2,IPNT)*W2+TABFJW( 3,IPNT)*W3
      RJ000Array(iP, 1) = TABFJW( 1,IPNT)-TABFJW( 2,IPNT)*WDIFF+TABFJW( 3,IPNT)*W2+TABFJW( 4,IPNT)*W3
      RJ000Array(iP, 2) = TABFJW( 2,IPNT)-TABFJW( 3,IPNT)*WDIFF+TABFJW( 4,IPNT)*W2+TABFJW( 5,IPNT)*W3
      RJ000Array(iP, 3) = TABFJW( 3,IPNT)-TABFJW( 4,IPNT)*WDIFF+TABFJW( 5,IPNT)*W2+TABFJW( 6,IPNT)*W3
      RJ000Array(iP, 4) = TABFJW( 4,IPNT)-TABFJW( 5,IPNT)*WDIFF+TABFJW( 6,IPNT)*W2+TABFJW( 7,IPNT)*W3
      RJ000Array(iP, 5) = TABFJW( 5,IPNT)-TABFJW( 6,IPNT)*WDIFF+TABFJW( 7,IPNT)*W2+TABFJW( 8,IPNT)*W3
      RJ000Array(iP, 6) = TABFJW( 6,IPNT)-TABFJW( 7,IPNT)*WDIFF+TABFJW( 8,IPNT)*W2+TABFJW( 9,IPNT)*W3
      RJ000Array(iP, 7) = TABFJW( 7,IPNT)-TABFJW( 8,IPNT)*WDIFF+TABFJW( 9,IPNT)*W2+TABFJW(10,IPNT)*W3
      RJ000Array(iP, 8) = TABFJW( 8,IPNT)-TABFJW( 9,IPNT)*WDIFF+TABFJW(10,IPNT)*W2+TABFJW(11,IPNT)*W3
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
      RJ000Array(iP,0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
      RJ000Array(iP, 4) = RJ000( 4)
      RJ000Array(iP, 5) = RJ000( 5)
      RJ000Array(iP, 6) = RJ000( 6)
      RJ000Array(iP, 7) = RJ000( 7)
      RJ000Array(iP, 8) = RJ000( 8)
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
      RJ000Array(iP, 0) = RJ000(0)
      RJ000Array(iP, 1) = RJ000( 1)
      RJ000Array(iP, 2) = RJ000( 2)
      RJ000Array(iP, 3) = RJ000( 3)
      RJ000Array(iP, 4) = RJ000( 4)
      RJ000Array(iP, 5) = RJ000( 5)
      RJ000Array(iP, 6) = RJ000( 6)
      RJ000Array(iP, 7) = RJ000( 7)
      RJ000Array(iP, 8) = RJ000( 8)
     ENDIF
  ENDDO
 end subroutine
end module
