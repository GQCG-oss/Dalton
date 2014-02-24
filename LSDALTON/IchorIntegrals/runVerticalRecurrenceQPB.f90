MODULE TESTMODULE
CONTAINS
  subroutine PASSsub
    IMPLICIT NONE
    INTEGER :: JMAX,nTUV,nTUVprev,ituvP,J,Tp,Up,Vp,N,N2,ntuvP,ituv,C,nTUVTMP,nTUVTMPprev
    Integer :: MaxAngmomQP,nTUVplus,JTMP,ntuvprev2,ntuvprev3
!    logical :: CREATED(-2:8,-2:8,-2:8)
    logical,pointer :: CREATED(:,:,:)
    logical :: TREC,UREC,VREC,TREC2,UREC2,VREC2
    integer,pointer :: TUVINDEX(:,:,:)
    integer,pointer :: TINDEX(:)
    integer,pointer :: UINDEX(:)
    integer,pointer :: VINDEX(:)
    integer,pointer :: JINDEX(:)
    integer :: nTUVLIST,nTUVLISTactual
    integer,pointer :: TwoTermTUVLIST(:)
    
    WRITE(*,'(A)')'MODULE AGC_OBS_VERTICALRECURRENCEMODB'
    WRITE(*,'(A)')' use IchorPrecisionModule'
    WRITE(*,'(A)')'  '
    WRITE(*,'(A)')' CONTAINS'
    MaxAngmomQP = 8

    WRITE(*,'(A)')''
    WRITE(*,'(A)')'subroutine VerticalRecurrence1B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
    WRITE(*,'(A)')'         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,&'
    WRITE(*,'(A)')'         & QpreExpFac,AUXarray)'
    WRITE(*,'(A)')'  implicit none'
    WRITE(*,'(A)')'  integer,intent(in) :: nPasses,nPrimP,nPrimQ'
    WRITE(*,'(A)')'  REAL(REALK),intent(in) :: TABFJW(0:4,0:1200)'
    WRITE(*,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)'
!    WRITE(*,'(A)')'!  REAL(REALK),intent(in) :: RJ000(0:1,nPrimQ,nPrimP,nPasses)'
    WRITE(*,'(A)')'  REAL(REALK),intent(in) :: integralPrefactor(nprimQ,nPrimP)'
    WRITE(*,'(A)')'  REAL(REALK),intent(in) :: QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)'
    WRITE(*,'(A)')'  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses),Bcenter(3)'
    WRITE(*,'(A)')'  real(realk),intent(inout) :: AUXarray(4,nPrimQ,nPrimP,nPasses)'
    WRITE(*,'(A)')'  !local variables'
    WRITE(*,'(A)')'  integer :: iPassQ,iPrimP,iPrimQ,ipnt'
    WRITE(*,'(A)')'  real(realk) :: Bx,By,Bz,Pexpfac,invexpP,mPX,mPY,mPZ,Xpb,Ypb,Zpb,RJ000(0:1)'
    WRITE(*,'(A)')'  real(realk) :: PREF,TMP1,TMP2,alphaP,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq'    
    WRITE(*,'(A)')'  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL' 
    WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk'
    WRITE(*,'(A)')'  real(realk),parameter :: D2=2.0E0_realk'
    WRITE(*,'(A,ES24.16,A)')'  REAL(REALK),PARAMETER :: D2JP36=',2.0d0 + 36.0d0,'_realk'
    WRITE(*,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D03333=1.0E0_realk/3.0E0_realk'
    WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk'
    WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk'
    WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: SMALL = 1E-15_realk,D12 = 12.0E0_realk'
    WRITE(*,'(A)')'  REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk'
    WRITE(*,'(A)')'  REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk'
    WRITE(*,'(A)')'  REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk'
    WRITE(*,'(A)')'  REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk'
    WRITE(*,'(A)')'  Real(realk), parameter :: PI=3.14159265358979323846E0_realk'
    WRITE(*,'(A)')'  REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk'
    WRITE(*,'(A)')'  REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2'
    WRITE(*,'(A)')'  REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI'
    WRITE(*,'(A)')'  !ThetaAux(n,1,0,0) = Xpb*ThetaAux(n,0,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,0,0,0)'
    WRITE(*,'(A)')'  !i = 0 last 2 term vanish'
    WRITE(*,'(A)')'  !We include scaling of RJ000 '
    WRITE(*,'(A)')'  Bx = -Bcenter(1)'
    WRITE(*,'(A)')'  By = -Bcenter(2)'
    WRITE(*,'(A)')'  Bz = -Bcenter(3)'
    WRITE(*,'(A)')'  DO iPassQ = 1,nPasses'
    WRITE(*,'(A)')'   DO iPrimP=1, nPrimP'
    WRITE(*,'(A)')'    Pexpfac = PpreExpFac(iPrimP)'
    WRITE(*,'(A)')'    invexpP = D1/Pexp(iPrimP)'
    WRITE(*,'(A)')'    mPX = -Pcent(1,iPrimP)'
    WRITE(*,'(A)')'    mPY = -Pcent(2,iPrimP)'
    WRITE(*,'(A)')'    mPZ = -Pcent(3,iPrimP)'
    WRITE(*,'(A)')'    Xpb = Pcent(1,iPrimP) + Bx'
    WRITE(*,'(A)')'    Ypb = Pcent(2,iPrimP) + By'
    WRITE(*,'(A)')'    Zpb = Pcent(3,iPrimP) + Bz'
    WRITE(*,'(A)')'    DO iPrimQ=1, nPrimQ'
    WRITE(*,'(A)')'     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)'
    WRITE(*,'(A)')'     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)'
    WRITE(*,'(A)')'     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)'
    WRITE(*,'(A)')'     alphaP = reducedExponents(iPrimQ,iPrimP)*invexpP'    
    WRITE(*,'(A)')'     alphaXpq = alphaP*Xpq'
    WRITE(*,'(A)')'     alphaYpq = alphaP*Ypq'
    WRITE(*,'(A)')'     alphaZpq = alphaP*Zpq'
    WRITE(*,'(A)')'     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq'
    WRITE(*,'(A)')'     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance'
    WRITE(*,'(A)')'     !  0 < WVAL < 0.000001'
    WRITE(*,'(A)')'!     IF (ABS(WVAL) .LT. SMALL) THEN'     
    WRITE(*,'(A)')'!      RJ000(0) = D1'
    WRITE(*,'(A)')'!      RJ000(1)= D03333 !THE BOYS FUNCTION FOR ZERO ARGUMENT'
    WRITE(*,'(A)')'!     !  0 < WVAL < 12 '
    WRITE(*,'(A)')'     IF (WVAL .LT. D12) THEN'
    WRITE(*,'(A)')'      IPNT = NINT(D100*WVAL)'
    WRITE(*,'(A)')'      WDIFF = WVAL - TENTH*IPNT'
    WRITE(*,'(A)')'      W2    = WDIFF*WDIFF'
    WRITE(*,'(A)')'      W3    = W2*WDIFF'
    WRITE(*,'(A)')'      W2    = W2*D05'
    WRITE(*,'(A)')'      W3    = W3*COEF3'
    WRITE(*,'(A)')'      RJ000(0)=TABFJW(0,IPNT)-TABFJW(1,IPNT)*WDIFF+TABFJW(2,IPNT)*W2+TABFJW(3,IPNT)*W3'
    WRITE(*,'(A)')'      RJ000(1)=TABFJW(1,IPNT)-TABFJW(2,IPNT)*WDIFF+TABFJW(3,IPNT)*W2+TABFJW(4,IPNT)*W3'
    WRITE(*,'(A)')'     !  12 < WVAL <= (2J+36) '
    WRITE(*,'(A)')'     ELSE IF (WVAL.LE.D2JP36) THEN'
    WRITE(*,'(A)')'      REXPW = D05*EXP(-WVAL)'
    WRITE(*,'(A)')'      RWVAL = D1/WVAL'
    WRITE(*,'(A)')'      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))'
    WRITE(*,'(A)')'      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
    WRITE(*,'(A)')'      RJ000(1) = RWVAL*(D05*RJ000(0)-REXPW)'
    WRITE(*,'(A)')'     !  (2J+36) < WVAL '
    WRITE(*,'(A)')'     ELSE'
    WRITE(*,'(A)')'      RWVAL = PID4/WVAL'
    WRITE(*,'(A)')'      RJ000(0) = SQRT(RWVAL)'
    WRITE(*,'(A)')'      RWVAL = RWVAL*PID4I'
    WRITE(*,'(A)')'      RJ000(1) = RWVAL*D05*RJ000(0)'
    WRITE(*,'(A)')'     ENDIF'
    WRITE(*,'(A)')'     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac'
    WRITE(*,'(A)')'     TMP1 = PREF*RJ000(0)'
    WRITE(*,'(A)')'     TMP2 = PREF*RJ000(1)'
    WRITE(*,'(A)')'     AUXarray(1,iPrimQ,iPrimP,iPassQ) = TMP1'
    WRITE(*,'(A)')'     AUXarray(2,iPrimQ,iPrimP,iPassQ) = Xpb*TMP1 + alphaXpq*TMP2'
    WRITE(*,'(A)')'     AUXarray(3,iPrimQ,iPrimP,iPassQ) = Ypb*TMP1 + alphaYpq*TMP2'
    WRITE(*,'(A)')'     AUXarray(4,iPrimQ,iPrimP,iPassQ) = Zpb*TMP1 + alphaZpq*TMP2'
    WRITE(*,'(A)')'    enddo'
    WRITE(*,'(A)')'   enddo'
    WRITE(*,'(A)')'  enddo'
    WRITE(*,'(A)')'end subroutine VerticalRecurrence1B'

    DO JMAX=2,MaxAngmomQP
       nTUV = (JMAX+1)*(JMAX+2)*(JMAX+3)/6   
       nTUVPLUS = (JMAX+2)*(JMAX+3)*(JMAX+4)/6   
       allocate(TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
       TUVINDEX = 0
       allocate(TINDEX(nTUVPLUS))
       allocate(UINDEX(nTUVPLUS))
       allocate(VINDEX(nTUVPLUS))
       allocate(JINDEX(nTUVPLUS))

       nTUVprev3 = (JMAX-2)*(JMAX-1)*(JMAX)/6
       nTUVprev2 = (JMAX-1)*(JMAX)*(JMAX+1)/6
       nTUVprev = (JMAX)*(JMAX+1)*(JMAX+2)/6
       ituv = 0 
       TUVINDEX = 0
       DO J = 0, JMAX
          DO Tp=J,0,-1       
             DO Up=J-Tp,0,-1
                Vp=J-Tp-Up
                ituv = ituv + 1 
                TUVINDEX(Tp,Up,Vp) = ituv
                TINDEX(iTUV) = Tp
                UINDEX(iTUV) = Up
                VINDEX(iTUV) = Vp
                JINDEX(iTUV) = J
             ENDDO
          ENDDO
       ENDDO

       WRITE(*,'(A)')''
       WRITE(*,'(A,I1,A)')'subroutine VerticalRecurrence',JMAX,'B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       WRITE(*,'(A)')'         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&'
       WRITE(*,'(A)')'         & AUXarray)'
       WRITE(*,'(A)')'  implicit none'
       WRITE(*,'(A)')'  integer,intent(in) :: nPasses,nPrimP,nPrimQ'
       WRITE(*,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)'
       WRITE(*,'(A,I1,A)')'!  REAL(REALK),intent(in) :: RJ000(0:',JMAX,',nPrimQ*nPrimP*nPasses)'
       WRITE(*,'(A)')'  REAL(REALK),intent(in) :: integralPrefactor(nprimQ,nPrimP)'
       WRITE(*,'(A)')'  REAL(REALK),intent(in) :: QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)'
       WRITE(*,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPrimQ*nPrimP*nPasses)'
       WRITE(*,'(A)')'  real(realk),intent(in) :: Pcent(3,nPrimP),Bcenter(3),Qcent(3,nPrimQ,nPasses)'
       WRITE(*,'(A,I2,A)')'  REAL(REALK),intent(in) :: TABFJW(0:',JMAX+3,',0:1200)'
       WRITE(*,'(A)')'  !Local variables'
       WRITE(*,'(A)')'  integer :: iPassQ,iPrimP,iPrimQ,IPNT,IP'
       WRITE(*,'(A)')'  real(realk) :: Bx,By,Bz,Pexpfac,invexpP,inv2expP,mPX,mPY,mPZ,Xpb,Ypb,Zpb,WVAL'
       WRITE(*,'(A)')'  real(realk) :: WDIFF,W2,W3,REXPW,RWVAL,GVAL'
       WRITE(*,'(A)')'  real(realk) :: PREF,alphaP,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq,squaredDistance'
       WRITE(*,'(A,i4,A)')'  real(realk) :: TwoTerms(',MAX(1,nTUVprev2-nTUVprev3),')'
       WRITE(*,'(A,i2,A)')'  real(realk) :: RJ000(0:',JMAX,')'
       C = JMAX+2
       DO JTMP=1,JMAX
          C = C-1
          nTUVTMPprev=(JTMP-1)*(JTMP)*(JTMP+1)/6
          nTUVTMP=(JTMP)*(JTMP+1)*(JTMP+2)/6
          if(JTMP.LT.10)THEN
             WRITE(*,'(A,I1,A,I3,A,I3,A,I1,A)')'  real(realk) :: TMParray',JTMP,'(',nTUVTMPprev+1,':',nTUVTMP,',2:',C,')'
          else
             WRITE(*,'(A,I2,A,I3,A,I3,A,I1,A)')'  real(realk) :: TMParray',JTMP,'(',nTUVTMPprev+1,':',nTUVTMP,',2:',C,')'
             
          endif
       ENDDO
       WRITE(*,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D03333=1.0E0_realk/3.0E0_realk'
       WRITE(*,'(A)')'  real(realk),parameter :: D2=2.0E0_realk,D4 = 4E0_realk, D100=100E0_realk'
       WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk'
       WRITE(*,'(A,ES24.16,A)')'  REAL(REALK),PARAMETER :: D2JP36=',2.d0*JMAX + 36.d0,'_realk'
       WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk'
       WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: SMALL = 1E-15_realk,D12 = 12.0E0_realk'
       WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk'
       WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk'
       WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: GFAC2 =  D2*0.321180909E0_realk'
       WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk'
       WRITE(*,'(A)')'  Real(realk),parameter :: PI=3.14159265358979323846E0_realk'
       WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk'
       WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: SQRPIH = SQRTPI/D2'
       WRITE(*,'(A)')'  REAL(REALK),PARAMETER :: PID4 = PI/D4, PID4I = D4/PI'
       WRITE(*,'(A)')'  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)'
       WRITE(*,'(A)')'  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))'
       WRITE(*,'(A)')'  !We include scaling of RJ000 '
       WRITE(*,'(A)')'  Bx = -Bcenter(1)'
       WRITE(*,'(A)')'  By = -Bcenter(2)'
       WRITE(*,'(A)')'  Bz = -Bcenter(3)'
       WRITE(*,'(A)')'  DO iPassQ = 1,nPasses'
       WRITE(*,'(A)')'   IP = (iPassQ-1)*nPrimQ*nPrimP'
       WRITE(*,'(A)')'   DO iPrimP=1, nPrimP'
       WRITE(*,'(A)')'    Pexpfac = PpreExpFac(iPrimP)'
       WRITE(*,'(A)')'    invexpP = D1/Pexp(iPrimP)'
       WRITE(*,'(A)')'    inv2expP = D05*invexpP'
       WRITE(*,'(A)')'    mPX = -Pcent(1,iPrimP)'
       WRITE(*,'(A)')'    mPY = -Pcent(2,iPrimP)'
       WRITE(*,'(A)')'    mPZ = -Pcent(3,iPrimP)'
       WRITE(*,'(A)')'    Xpb = Pcent(1,iPrimP) + Bx'
       WRITE(*,'(A)')'    Ypb = Pcent(2,iPrimP) + By'
       WRITE(*,'(A)')'    Zpb = Pcent(3,iPrimP) + Bz'
       WRITE(*,'(A)')'    DO iPrimQ=1, nPrimQ'
       WRITE(*,'(A)')'     IP = IP + 1'
       WRITE(*,'(A)')'     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP'
       WRITE(*,'(A)')'     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)'
       WRITE(*,'(A)')'     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)'
       WRITE(*,'(A)')'     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)'
       WRITE(*,'(A)')'     alphaXpq = -alphaP*Xpq'
       WRITE(*,'(A)')'     alphaYpq = -alphaP*Ypq'
       WRITE(*,'(A)')'     alphaZpq = -alphaP*Zpq'
       WRITE(*,'(A)')'     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq'
       WRITE(*,'(A)')'     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance'
       WRITE(*,'(A)')'     !  0 < WVAL < 0.000001'
       WRITE(*,'(A)')'!     IF (ABS(WVAL) .LT. SMALL) THEN'     
       WRITE(*,'(A)')'!      RJ000(0) = D1 !THE BOYS FUNCTION FOR ZERO ARGUMENT'
       DO J=1,JMAX
       WRITE(*,'(A,I1,A,ES23.16,A)')'!      RJ000(',J,')= ',1.0d0/(2*J + 1),'_realk'
       ENDDO
       WRITE(*,'(A)')'     !  0 < WVAL < 12 '
       WRITE(*,'(A)')'     IF (WVAL .LT. D12) THEN'
       WRITE(*,'(A)')'      IPNT = NINT(D100*WVAL)'
       WRITE(*,'(A)')'      WDIFF = WVAL - TENTH*IPNT'
       WRITE(*,'(A)')'      W2    = WDIFF*WDIFF'
       WRITE(*,'(A)')'      W3    = W2*WDIFF'
       WRITE(*,'(A)')'      W2    = W2*D05'
       WRITE(*,'(A)')'      W3    = W3*COEF3'
       DO J=0,JMAX
       WRITE(*,'(A,I2,A,I2,A,I2,A,I2,A,I2,A)')'      RJ000(',J,') = TABFJW(',J,',IPNT)-TABFJW(',J+1,',IPNT)*WDIFF+TABFJW(',J+2,',IPNT)*W2+TABFJW(',J+3,',IPNT)*W3'
       ENDDO
   !    WRITE(*,'(A)')'        DO J=0,JMAX'
  !     WRITE(*,'(A)')'           R = TABFJW(J,IPNT)'
  !     WRITE(*,'(A)')'           R = R -TABFJW(J+1,IPNT)*WDIFF'
  !     WRITE(*,'(A)')'           R = R + TABFJW(J+2,IPNT)*W2'
  !     WRITE(*,'(A)')'           R = R + TABFJW(J+3,IPNT)*W3'
  !     WRITE(*,'(A)')'           RJ000(J,ipq,ipassq) = R'
   !    WRITE(*,'(A)')'        ENDDO'
       WRITE(*,'(A)')'     !  12 < WVAL <= (2J+36) '
       WRITE(*,'(A)')'     ELSE IF (WVAL.LE.D2JP36) THEN'
       WRITE(*,'(A)')'      REXPW = D05*EXP(-WVAL)'
       WRITE(*,'(A)')'      RWVAL = D1/WVAL'
       WRITE(*,'(A)')'      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))'
          WRITE(*,'(A)')'      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
!       WRITE(*,'(A)')'        RJ000(0,ipq,ipassq) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
       DO J=1,JMAX
          WRITE(*,'(A,I2,A,I2,A,I2,A)')'      RJ000(',J,') = RWVAL*((',J,' - D05)*RJ000(',J-1,')-REXPW)'          
       ENDDO
!       WRITE(*,'(A)')'        DO J=1,JMAX'
!       WRITE(*,'(A)')'           RJ000(J,ipq,ipassq) = RWVAL*((J - D05)*RJ000(J-1,ipq,ipassq)-REXPW)'
!       WRITE(*,'(A)')'        ENDDO'
       WRITE(*,'(A)')'     !  (2J+36) < WVAL '
       WRITE(*,'(A)')'     ELSE'
       WRITE(*,'(A)')'      RWVAL = PID4/WVAL'
       WRITE(*,'(A)')'      RJ000(0) = SQRT(RWVAL)'
!       WRITE(*,'(A)')'        RJ000(0,ipq,ipassq) = SQRT(RWVAL)'
       WRITE(*,'(A)')'      RWVAL = RWVAL*PID4I'
       DO J=1,JMAX
          WRITE(*,'(A,I2,A,I2,A,I2,A)')'      RJ000(',J,') = RWVAL*(',J,' - D05)*RJ000(',J-1,')'
       ENDDO
!       WRITE(*,'(A)')'        DO J = 1, JMAX'
!       WRITE(*,'(A)')'           RJ000(J,ipq,ipassq) = RWVAL*(J - D05)*RJ000(J-1,ipq,ipassq)'
!       WRITE(*,'(A)')'        ENDDO'
       WRITE(*,'(A)')'     ENDIF'
       WRITE(*,'(A)')'     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac'
       WRITE(*,'(A,I1,A,I1)')'     Auxarray(1,IP) = PREF*RJ000(0)'
       DO J=2,JMAX+1
          WRITE(*,'(A,I2,A,I2,A)')'     TMParray1(1,',J,') = PREF*RJ000(',J-1,')'
       ENDDO

       allocate(CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
       CREATED  = .FALSE.
       CREATED(0,0,0) = .TRUE.
       DO JTMP=1,JMAX

          DO J = 1,JMAX-JTMP+1

          !determine twoterms
          nTUVLIST = (JTMP+1)*(JTMP+2)*(JTMP+3)/6   
          allocate(TwoTermTUVLIST(nTUVLIST))
          call TwoTerms(J,JTMP,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)

           DO Tp=JTMP,0,-1       
             DO Up=JTMP-Tp,0,-1
              Vp=JTMP-Tp-Up                
              iTUV = TUVINDEX(Tp,Up,Vp)

!                print*,'!iTUV=',iTUV
                !step 1.
                !how can the (Tp,Up,Vp) be built from lower
                TREC = CREATED(Tp-1,Up,Vp)
                UREC = CREATED(Tp,Up-1,Vp)
                VREC = CREATED(Tp,Up,Vp-1)
                N=0
                IF(TREC)N=N+1
                IF(UREC)N=N+1
                IF(VREC)N=N+1
                IF(N.EQ.1)THEN
                   !only one possible way to construct it
                   IF(TREC)THEN
!                      print*,'!A TRECURRENCE iTUV',iTUV
                      call TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
                   ELSEIF(UREC)THEN
!                      print*,'!A URECURRENCE iTUV',iTUV
                      call URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
                   ELSEIF(VREC)THEN
!                      print*,'!A VRECURRENCE iTUV',iTUV
                      call VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
                   ENDIF
                ELSE
                   !several ways to construct it
                   TREC2 = CREATED(Tp-2,Up,Vp)
                   UREC2 = CREATED(Tp,Up-2,Vp)
                   VREC2 = CREATED(Tp,Up,Vp-2)
                   N2=0
                   IF(TREC2)N2=N2+1
                   IF(UREC2)N2=N2+1
                   IF(VREC2)N2=N2+1
                   IF(N2.LT.N)THEN
                      !two term recurrence possible for one or more the possibilities 
                      !we chose the one term possibility
                      IF(.NOT.(TREC.AND.TREC2).AND.TREC)THEN
!                         print*,'!B TRECURRENCE iTUV',iTUV
                         CALL TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
                      ELSEIF(.NOT.(UREC.AND.UREC2).AND.UREC)THEN
!                         print*,'!B URECURRENCE iTUV',iTUV
                         CALL URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
                      ELSEIF(.NOT.(VREC.AND.VREC2).AND.VREC)THEN
!                         print*,'!B VRECURRENCE iTUV',iTUV
                         CALL VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
                      ENDIF
                   ELSE
                      !chose one of the possibilities
                      IF(TREC)THEN
!                         print*,'!C TRECURRENCE iTUV',iTUV
                         CALL TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
                      ELSEIF(UREC)THEN
!                         print*,'!C URECURRENCE iTUV',iTUV
                         call URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
                      ELSEIF(VREC)THEN
!                         print*,'!D VRECURRENCE iTUV',iTUV
                         call VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
                      ENDIF
                   ENDIF
                ENDIF
                CREATED(Tp,Up,Vp) = .TRUE.
              ENDDO
           ENDDO
          ENDDO
           deallocate(TwoTermTUVLIST)
          !WRITE(*,'(A,I4)')'           ENDDO'
       ENDDO
       WRITE(*,'(A)')'    ENDDO'
       WRITE(*,'(A)')'   ENDDO'
       WRITE(*,'(A)')'  ENDDO'
       WRITE(*,'(A,I1,A)')'end subroutine VerticalRecurrence',JMAX,'B'
       deallocate(TUVINDEX)
       deallocate(TINDEX)
       deallocate(UINDEX)
       deallocate(VINDEX)
       deallocate(JINDEX)
    ENDDO
    WRITE(*,'(A)')'end module'
  END subroutine PASSsub
  
subroutine TwoTerms(J1,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,JTMP,I,J1
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: nTUVLIST,nTUVLISTactual
integer :: TwoTermTUVLIST(nTUVLIST)
logical :: Unique,TREC,UREC,VREC
!TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
!               +T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
!TwoTerms are all  TUV(T-2) that is non zero ! to start with maybe not all needed
TwoTermTUVLIST = 0 
nTUVLISTactual = 0 

DO Tp=J,0,-1       
   DO Up=J-Tp,0,-1
      Vp=J-Tp-Up       
      IF(Tp-2.GE.0)THEN
         ituvP2 = TUVINDEX(Tp-2,Up,Vp)
         Unique = .TRUE.
         DO I = 1, nTUVLISTactual
            IF(ituvP2.EQ.TwoTermTUVLIST(I))Unique=.FALSE.
         ENDDO
         IF(Unique)THEN
            nTUVLISTactual = nTUVLISTactual + 1
            TwoTermTUVLIST(nTUVLISTactual) = iTUVP2
         ENDIF
      ENDIF
      IF(Up-2.GE.0)THEN
         ituvP2 = TUVINDEX(Tp,Up-2,Vp)
         Unique = .TRUE.
         DO I = 1, nTUVLISTactual
            IF(ituvP2.EQ.TwoTermTUVLIST(I))Unique=.FALSE.
         ENDDO
         IF(Unique)THEN
            nTUVLISTactual = nTUVLISTactual + 1
            TwoTermTUVLIST(nTUVLISTactual) = iTUVP2
         ENDIF
      ENDIF
      IF(Vp-2.GE.0)THEN
         ituvP2 = TUVINDEX(Tp,Up,Vp-2)
         Unique = .TRUE.
         DO I = 1, nTUVLISTactual
            IF(ituvP2.EQ.TwoTermTUVLIST(I))Unique=.FALSE.
         ENDDO
         IF(Unique)THEN
            nTUVLISTactual = nTUVLISTactual + 1
            TwoTermTUVLIST(nTUVLISTactual) = iTUVP2
         ENDIF
      ENDIF
   ENDDO
ENDDO

!!$DO ituvP2 = 1,nTUVLISTactual
!!$   IF(J1.EQ.1)THEN
!!$      IF(JTMP-1.LT.10)THEN
!!$         WRITE(*,'(A,I3,A,I3,A,I1,A,I3,A)')&
!!$              &'     TwoTerms(',ituvP2,') = invexpP*(AuxArray(',TwoTermTUVLIST(ituvP2),',IP) + alphaP*TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',2))'
!!$      ELSE
!!$         WRITE(*,'(A,I3,A,I3,A,I2,A,I3,A)')&
!!$              &'     TwoTerms(',ituvP2,') = invexpP*(AuxArray(',TwoTermTUVLIST(ituvP2),',IP) + alphaP*TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',2))'
!!$      ENDIF
!!$   ELSE
!!$      IF(JTMP-1.LT.10)THEN
!!$         WRITE(*,'(A,I3,A,I1,A,I3,A,I3,A,I1,A,I3,A,I3,A)')&
!!$              &'     TwoTerms(',ituvP2,') = invexpP*(TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',',J1,') + alphaP*TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',',J1+1,'))'
!!$      ELSE
!!$         WRITE(*,'(A,I3,A,I1,A,I3,A,I3,A,I2,A,I3,A,I3,A)')&
!!$              &'     TwoTerms(',ituvP2,') = invexpP*(TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',',J1,') + alphaP*TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',',J1+1,'))'
!!$      ENDIF
!!$   ENDIF
!!$ENDDO
DO ituvP2 = 1,nTUVLISTactual
   IF(J1.EQ.1)THEN
      IF(JTMP-1.LT.10)THEN
         WRITE(*,'(A,I3,A,I3,A,I1,A,I3,A)')&
              &'     TwoTerms(',ituvP2,') = inv2expP*(AuxArray(',TwoTermTUVLIST(ituvP2),',IP) + alphaP*TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',2))'
      ELSE
         WRITE(*,'(A,I3,A,I3,A,I2,A,I3,A)')&
              &'     TwoTerms(',ituvP2,') = inv2expP*(AuxArray(',TwoTermTUVLIST(ituvP2),',IP) + alphaP*TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',2))'
      ENDIF
   ELSE
      IF(JTMP-1.LT.10)THEN
         WRITE(*,'(A,I3,A,I1,A,I3,A,I3,A,I1,A,I3,A,I3,A)')&
              &'     TwoTerms(',ituvP2,') = inv2expP*(TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',',J1,') + alphaP*TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',',J1+1,'))'
      ELSE
         WRITE(*,'(A,I3,A,I1,A,I3,A,I3,A,I2,A,I3,A,I3,A)')&
              &'     TwoTerms(',ituvP2,') = inv2expP*(TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',',J1,') + alphaP*TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',',J1+1,'))'
      ENDIF
   ENDIF
ENDDO

end subroutine TwoTerms

subroutine TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: nTUVLIST,nTUVLISTactual,JTMP,TM2
integer :: TwoTermTUVLIST(nTUVLIST)

!TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
!               +T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
!TwoTerms are all  TUV(T-2) that is non zero ! to start with maybe not all needed
ituvP0 = TUVINDEX(Tp,Up,Vp)
ituvP1 = TUVINDEX(Tp-1,Up,Vp)

ituvP2 = 1
TM1 = Tp-1
TM2 = Tp-2

iTwoTerms = 1
IF(TM2.GE.0)THEN
   ituvP2 = TUVINDEX(Tp-2,Up,Vp)
   DO I=1,nTUVLISTactual
      IF(ituvP2.EQ.TwoTermTUVLIST(I))THEN
         iTwoTerms = I
      ENDIF
   ENDDO
ENDIF

IF(J.EQ.JMAX)THEN
   IF(J.EQ.1)THEN
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            STOP 'TREC D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A)')&
              &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2)'
      ENDIF
   ELSE
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            stop 'Trec B D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A)')&
              &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,')'
      ENDIF      
   ENDIF
ELSE
   IF(J.EQ.1)THEN
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2) + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            STOP 'Trec C D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A)')&
              &'     AuxArray(',ituvP0,',IP) = Xpb*AuxArray(',ituvP1,',IP) + alphaXpq*TmpArray',JTMP,'(',ituvP1,',2)'
      ENDIF
   ELSE
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            STOP 'TREC D D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A)')&
              &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Xpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaXpq*TmpArray',JTMP,'(',ituvP1,',',J+1,')'
      ENDIF
   ENDIF
ENDIF
end subroutine TRECURRENCE

subroutine URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: nTUVLIST,nTUVLISTactual,JTMP,TM2
integer :: TwoTermTUVLIST(nTUVLIST)

ituvP0 = TUVINDEX(Tp,Up,Vp)
ituvP1 = TUVINDEX(Tp,Up-1,Vp)

ituvP2 = 1
TM1 = Up-1
TM2 = Up-2

iTwoTerms = 1
IF(TM2.GE.0)THEN
   ituvP2 = TUVINDEX(Tp,Up-2,Vp)
   DO I=1,nTUVLISTactual
      IF(ituvP2.EQ.TwoTermTUVLIST(I))THEN
         iTwoTerms = I
      ENDIF
   ENDDO
ENDIF

IF(J.EQ.JMAX)THEN
   IF(J.EQ.1)THEN
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            STOP 'TREC D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A)')&
              &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2)'
      ENDIF
   ELSE
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            stop 'Trec B D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A)')&
              &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,')'
      ENDIF      
   ENDIF
ELSE
   IF(J.EQ.1)THEN
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2) + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            STOP 'Trec C D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A)')&
              &'     AuxArray(',ituvP0,',IP) = Ypb*AuxArray(',ituvP1,',IP) + alphaYpq*TmpArray',JTMP,'(',ituvP1,',2)'
      ENDIF
   ELSE
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            STOP 'TREC D D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A)')&
              &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Ypb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaYpq*TmpArray',JTMP,'(',ituvP1,',',J+1,')'
      ENDIF
   ENDIF
ENDIF
end subroutine URECURRENCE

subroutine VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: nTUVLIST,nTUVLISTactual,JTMP,TM2
integer :: TwoTermTUVLIST(nTUVLIST)

ituvP0 = TUVINDEX(Tp,Up,Vp)
ituvP1 = TUVINDEX(Tp,Up,Vp-1)

ituvP2 = 1
TM1 = Vp-1
TM2 = Vp-2

iTwoTerms = 1
IF(TM2.GE.0)THEN
   ituvP2 = TUVINDEX(Tp,Up,Vp-2)
   DO I=1,nTUVLISTactual
      IF(ituvP2.EQ.TwoTermTUVLIST(I))THEN
         iTwoTerms = I
      ENDIF
   ENDDO
ENDIF

IF(J.EQ.JMAX)THEN
   IF(J.EQ.1)THEN
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            STOP 'TREC D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A)')&
              &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2)'
      ENDIF
   ELSE
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            stop 'Trec B D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A)')&
              &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,')'
      ENDIF      
   ENDIF
ELSE
   IF(J.EQ.1)THEN
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A,I3,A)')&
                 &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2) + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            STOP 'Trec C D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I4,A,I4,A,I1,A,I4,A)')&
              &'     AuxArray(',ituvP0,',IP) = Zpb*AuxArray(',ituvP1,',IP) + alphaZpq*TmpArray',JTMP,'(',ituvP1,',2)'
      ENDIF
   ELSE
      IF(TM2.GE.0)THEN
         !four term relation
         IF(TM1.EQ.1)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.2)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 2.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.3)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 3.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.4)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 4.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.5)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 5.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.6)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 6.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.7)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 7.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.8)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 8.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.9)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 9.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSEIF(TM1.EQ.10)THEN
            WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I3,A)')&
                 &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,') + 10.0E0_realk*TwoTerms(',iTwoTerms,')'
         ELSE
            STOP 'TREC D D4'
         ENDIF
      ELSE
         !two term relation
         WRITE(*,'(A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A,I1,A,I4,A,I2,A)')&
              &'     tmpArray',JTMP+1,'(',ituvP0,',',J,') = Zpb*tmpArray',JTMP,'(',ituvP1,',',J,') + alphaZpq*TmpArray',JTMP,'(',ituvP1,',',J+1,')'
      ENDIF
   ENDIF
ENDIF
end subroutine VRECURRENCE



end MODULE TESTMODULE

PROGRAM PASSprogram
use TESTMODULE
INTEGER :: JMAX,JMIN

call PASSSUB

end PROGRAM PASSprogram
