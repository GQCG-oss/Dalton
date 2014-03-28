MODULE TESTMODULE
  use stringsMODULE
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
    integer :: nTUVLIST,nTUVLISTactual,I
    integer,pointer :: TwoTermTUVLIST(:)
    logical,pointer :: TwoTermsUsed(:)
    integer :: iseglabel,lufile,iseg,iPrim,nPrim
    integer :: non1Prim(16),pure1Prim(4),ia,ib,ic,id,GPUrun,K,center
    logical :: Gen,SegQ,SegP,Seg,Seg1Prim,CPU,nPrimnTUV,DoOpenMP
    logical :: Collapse,segwtuv,DoOpenACC
    Character(len=48) :: FileName    
    character(len=3) :: ARCSTRING,Xdir,Ydir,Zdir
    character(len=1) :: centerString
    Character(len=8)  :: SegLabel
    Character(len=20) :: PrimLabel,nPrimLabel

!=====================================================================================================0
! Vertical
!=====================================================================================================
    DO GPUrun = 1,2
       CPU = .TRUE.
       IF(GPUrun.EQ.2)CPU = .FALSE.
       DoOpenMP = .FALSE.
       DoOpenACC = .FALSE.
       IF(CPU)DoOpenMP = .TRUE.
       IF(.NOT.CPU)DoOpenACC = .TRUE.
       Collapse = .TRUE.
!       IF(.NOT.CPU)Collapse = .TRUE.
!       Collapse = .FALSE.
       IF(.NOT.CPU)nPrimnTUV = .FALSE.
       IF(CPU)THEN
          ARCSTRING = 'CPU'
       ELSE
          ARCSTRING = 'GPU'
       ENDIF
       do center=1,4
          IF(center.EQ.1)centerString='A'
          IF(center.EQ.2)centerString='B'
          IF(center.EQ.3)centerString='C'
          IF(center.EQ.4)centerString='D'
          DO iseg = 1,5
             Gen = .FALSE.; SegQ=.FALSE.; SegP=.FALSE.;Seg=.FALSE.;Seg1Prim=.FALSE.
             IF(iseg.EQ.1)THEN
                Gen = .TRUE.      ; SegLabel = 'Gen     '; iSegLabel = 3
             ELSEIF(iseg.EQ.2)THEN
                SegQ = .TRUE.     ; SegLabel = 'SegQ    '; iSegLabel = 4
             ELSEIF(iseg.EQ.3)THEN
                SegP = .TRUE.     ; SegLabel = 'SegP    '; iSegLabel = 4
             ELSEIF(iseg.EQ.4)THEN
                Seg = .TRUE.      ; SegLabel = 'Seg     '; iSegLabel = 3
             ELSEIF(iseg.EQ.5)THEN
                Seg1Prim = .TRUE. ; SegLabel = 'Seg1Prim'; iSegLabel = 8
             ENDIF
             DO I =1,48
                FileName(I:I) = ' '
             ENDDO
             WRITE(FileName,'(4A)')'runVerticalRecurrence'//ARCSTRING//'QP',centerString,SegLabel(1:iSegLabel),'.F90'
             print*,'FileName:',FileName
             LUFILE = 12
             open(unit = LUFILE, file=TRIM(FileName),status="unknown")

             IF(COLLAPSE)THEN
                IF(Gen)THEN
                   PrimLabel = 'iP'; iPrim = 2
                   nPrimLabel = 'nPrimQ*nPrimP*nPassP'; nPrim = 20 
                ELSEIF(SegQ)THEN
!                   PrimLabel = 'iPrimP,iPassP'; iPrim = 13 
!                   nPrimLabel = 'nPrimP,nPassP'; nPrim = 13 
                   PrimLabel = 'iP'; iPrim = 2
                   nPrimLabel = 'nPrimP*nPassP'; nPrim = 13 
                ELSEIF(SegP)THEN
!                   PrimLabel = 'iPrimQ,iPassP'; iPrim = 13 
!                   nPrimLabel = 'nPrimQ,nPassP'; nPrim = 13 
                   PrimLabel = 'iP'; iPrim = 2
                   nPrimLabel = 'nPrimQ*nPassP'; nPrim = 13 
                ELSEIF(Seg)THEN
                   PrimLabel = 'iPassP'; iPrim = 2
                   nPrimLabel = 'nPassP'; nPrim = 6 
                ELSE
                   PrimLabel = 'iPassP'; iPrim = 2
                   nPrimLabel = 'nPassP'; nPrim = 6 
                ENDIF
             ELSE
                IF(Gen)THEN
                   PrimLabel = 'iPrimQ,iPrimP,iPassP'; iPrim = 20 
                   nPrimLabel = 'nPrimQ,nPrimP,nPassP'; nPrim = 20 
                ELSEIF(SegQ)THEN
                   PrimLabel = 'iPrimP,iPassP'; iPrim = 13 
                   nPrimLabel = 'nPrimP,nPassP'; nPrim = 13 
                ELSEIF(SegP)THEN
                   PrimLabel = 'iPrimQ,iPassP'; iPrim = 13 
                   nPrimLabel = 'nPrimQ,nPassP'; nPrim = 13 
                ELSEIF(Seg)THEN
                   PrimLabel = 'iPassP'; iPrim = 6 
                   nPrimLabel = 'nPassP'; nPrim = 6 
                ELSE
                   PrimLabel = 'iPassP'; iPrim = 6 
                   nPrimLabel = 'nPassP'; nPrim = 6 
                ENDIF
             ENDIF

             WRITE(LUFILE,'(5A)')'MODULE AGC_',ARCSTRING,'_OBS_VERTICALRECURRENCEMOD',centerString,SegLabel(1:iSegLabel)
             MaxAngmomQP = 8 !currently only D functions
             IF((SegQ.OR.SegP).OR.Seg)THEN
                MaxAngmomQP = 4 
                !Highest possible is (DD|SS) otherwise
                !a General Vertical Recurrence is required followed by
                !a ElectronTransfer
             ENDIF
             WRITE(LUFILE,'(A)')' use IchorPrecisionModule'
             WRITE(LUFILE,'(A)')'  '
             WRITE(LUFILE,'(A)')' CONTAINS'

             !========================================================================================================
             !    VerticalRecurrence 0 only in Acenter
             !========================================================================================================
             WRITE(LUFILE,'(A)')''
             JMAX = 0
             IF(center.EQ.1)THEN
                WRITE(LUFILE,'(3A)')'subroutine VerticalRecurrence'//ARCSTRING,SegLabel(1:iSegLabel),'0(nPassP,nPrimP,nPrimQ,&'
                WRITE(LUFILE,'(A)')'         & reducedExponents,TABFJW,&'
                WRITE(LUFILE,'(A)')'         & Pcent,Qcent,integralPrefactor,&'
                WRITE(LUFILE,'(A)')'         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&'
                WRITE(LUFILE,'(A)')'         & AUXarray)'
                WRITE(LUFILE,'(A)')'  implicit none'
                WRITE(LUFILE,'(A)')'  integer,intent(in) :: nPassP,nPrimP,nPrimQ'
                WRITE(LUFILE,'(A)')'  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB'
                WRITE(LUFILE,'(A)')'  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)'
                WRITE(LUFILE,'(A)')'  REAL(REALK),intent(in) :: TABFJW(0:3,0:1200)'
                IF(Seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'  REAL(REALK),intent(in) :: reducedExponents(1)'
                   WRITE(LUFILE,'(A)')'  REAL(REALK),intent(in) :: integralPrefactor(1)'
                   WRITE(LUFILE,'(A)')'  REAL(REALK),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3)'
                   WRITE(LUFILE,'(A)')'  REAL(REALK),intent(in) :: QpreExpFac(1),PpreExpFac(nAtomsA,nAtomsB)'
                ELSE
                   WRITE(LUFILE,'(A)')'  REAL(REALK),intent(in) :: reducedExponents(nPrimQ,nPrimP)'
                   WRITE(LUFILE,'(A)')'  REAL(REALK),intent(in) :: integralPrefactor(nprimQ,nPrimP)'
                   WRITE(LUFILE,'(A)')'  REAL(REALK),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)'
                   WRITE(LUFILE,'(A)')'  REAL(REALK),intent(in) :: QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)'
                ENDIF
                WRITE(LUFILE,'(A,A,A)')'  real(realk),intent(inout) :: AUXarray(',nPrimLabel(1:nPrim),')'
                WRITE(LUFILE,'(A)')'  !local variables'
                WRITE(LUFILE,'(A,ES24.16,A)')'  REAL(REALK),PARAMETER :: D2JP36=',36.0d0,'_realk'
                WRITE(LUFILE,'(A)')'  real(realk),parameter :: D2=2.0E0_realk'
                WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: D05 =0.5E0_realk,D1=1E0_realk'
                WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk'
                WRITE(LUFILE,'(A)')'  Real(realk),parameter :: D12 = 12E0_realk, TENTH = 0.01E0_realk'
                WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk'
                WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk'
                WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk'
                WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk'
                WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: GFAC2 =  D2*0.321180909E0_realk'
                WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk'
                WRITE(LUFILE,'(A)')'  Real(realk),parameter :: PI=3.14159265358979323846E0_realk'
                WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk'
                WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: SQRPIH = SQRTPI/D2'
                WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: PID4 = PI/D4, PID4I = D4/PI'
                WRITE(LUFILE,'(A)')'!  REAL(REALK),PARAMETER :: SMALL = 1E-15_realk'
                WRITE(LUFILE,'(A)')'  Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA,WVAL'
                WRITE(LUFILE,'(A)')'  Real(realk) :: W2,W3,PX,PY,PZ,XPQ,YPQ,ZPQ,squaredDistance,RJ000'
                WRITE(LUFILE,'(A)')'  Integer :: IPNT,iAtomA,iAtomB'
                IF(COLLAPSE.AND.Seg)THEN
                   WRITE(LUFILE,'(A)')'  Integer :: iPrimQP'
                ENDIF
                WRITE(LUFILE,'(A)')'  Integer :: iP,iPrimQ,iPrimP,iPassP'
                IF(COLLAPSE)THEN
                   call PrintCollapseInitLoop(Gen,SegQ,SegP,Seg,seg1prim,LUFILE,0,nTUV,DoOpenMP)
                ENDIF
                IF(DoOpenMP.OR.DoOpenACC)THEN
                   call PrintOpenMP(Gen,SegQ,SegP,Seg,seg1prim,LUFILE,0,Collapse,Center,centerstring,DoOpenMP,DoOpenACC)
                ENDIF
!                         WRITE(LUFILE,'(A)')'  DO iP = 1,nPrimQ*nPrimP*nPassP'
!                         WRITE(LUFILE,'(A)')'   iPrimQ = mod(IP-1,nPrimQ)+1'
!                         WRITE(LUFILE,'(A)')'   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1'
!                         WRITE(LUFILE,'(A)')'   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1'
!                i1 = i12 - ((i12-1)/n1)*n1
!                i2 = (i12-1)/n1+1

!DO i123 = 1,n1*n2*n3'
!i1 = mod(I123-1,n1)+1'
!i2 = mod((I123-(mod(I123-1,n1)+1))/n1,n2)+1'
!i3 = (I123-1)/(n1*n2) + 1'


                IF(Collapse)THEN
                   call PrintCollapseLoopStart(Gen,SegQ,SegP,Seg,seg1prim,LUFILE)
                   WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
                   WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    px = Pcent(1,iPrimP,iAtomA,iAtomB)'
                      WRITE(LUFILE,'(A)')'    py = Pcent(2,iPrimP,iAtomA,iAtomB)'
                      WRITE(LUFILE,'(A)')'    pz = Pcent(3,iPrimP,iAtomA,iAtomB)'
                      WRITE(LUFILE,'(A)')'     Xpq = px - Qcent(1,iPrimQ)'
                      WRITE(LUFILE,'(A)')'     Ypq = py - Qcent(2,iPrimQ)'
                      WRITE(LUFILE,'(A)')'     Zpq = pz - Qcent(3,iPrimQ)'
                   ELSE
                      WRITE(LUFILE,'(A)')'     Xpq = Pcent(1,iAtomA,iAtomB) - Qcent(1)'
                      WRITE(LUFILE,'(A)')'     Ypq = Pcent(2,iAtomA,iAtomB) - Qcent(2)'
                      WRITE(LUFILE,'(A)')'     Zpq = Pcent(3,iAtomA,iAtomB) - Qcent(3)'
                   ENDIF
                ELSE
                   WRITE(LUFILE,'(A)')'  DO iPassP = 1,nPassP'
                   WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
                   WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
                   IF(Seg)WRITE(LUFILE,'(A)')'   AUXarray(iPassP)=0.0E0_realk'
                   IF(SegP)THEN
                      WRITE(LUFILE,'(A)')'   DO iPrimQ=1, nPrimQ'
                      WRITE(LUFILE,'(A)')'    AUXarray(iPrimQ,iPassP)=0.0E0_realk'
                      WRITE(LUFILE,'(A)')'   ENDDO'
                   ENDIF
                   IF(.NOT.Seg1Prim)WRITE(LUFILE,'(A)')'   DO iPrimP=1, nPrimP'
                   IF(SegQ)WRITE(LUFILE,'(A)')'    AUXarray(iPrimP,iPassP)=0.0E0_realk'
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    px = Pcent(1,iPrimP,iAtomA,iAtomB)'
                      WRITE(LUFILE,'(A)')'    py = Pcent(2,iPrimP,iAtomA,iAtomB)'
                      WRITE(LUFILE,'(A)')'    pz = Pcent(3,iPrimP,iAtomA,iAtomB)'
                      WRITE(LUFILE,'(A)')'    DO iPrimQ=1, nPrimQ'
                      WRITE(LUFILE,'(A)')'     Xpq = px - Qcent(1,iPrimQ)'
                      WRITE(LUFILE,'(A)')'     Ypq = py - Qcent(2,iPrimQ)'
                      WRITE(LUFILE,'(A)')'     Zpq = pz - Qcent(3,iPrimQ)'
                   ELSE
                      WRITE(LUFILE,'(A)')'     Xpq = Pcent(1,iAtomA,iAtomB) - Qcent(1)'
                      WRITE(LUFILE,'(A)')'     Ypq = Pcent(2,iAtomA,iAtomB) - Qcent(2)'
                      WRITE(LUFILE,'(A)')'     Zpq = Pcent(3,iAtomA,iAtomB) - Qcent(3)'
                   ENDIF
                ENDIF
                WRITE(LUFILE,'(A)')'     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq'
                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance'
                ELSE
                   WRITE(LUFILE,'(A)')'     WVAL = reducedExponents(1)*squaredDistance'
                ENDIF
                WRITE(LUFILE,'(A)')'     !  0 < WVAL < 0.000001'
                WRITE(LUFILE,'(A)')'!     IF (ABS(WVAL) .LT. SMALL) THEN'     
                WRITE(LUFILE,'(A)')'!      RJ000 = D1'
                WRITE(LUFILE,'(A)')'!     !  0 < WVAL < 12 '
                WRITE(LUFILE,'(A)')'     IF (WVAL .LT. D12) THEN'
                WRITE(LUFILE,'(A)')'      IPNT = NINT(D100*WVAL)'
                WRITE(LUFILE,'(A)')'      WDIFF = WVAL - TENTH*IPNT'
                WRITE(LUFILE,'(A)')'      W2    = WDIFF*WDIFF'
                WRITE(LUFILE,'(A)')'      W3    = W2*WDIFF'
                WRITE(LUFILE,'(A)')'      W2    = W2*D05'
                WRITE(LUFILE,'(A)')'      W3    = W3*COEF3'
                WRITE(LUFILE,'(A)')'      RJ000 = TABFJW(0,IPNT)-TABFJW(1,IPNT)*WDIFF+TABFJW(2,IPNT)*W2+TABFJW(3,IPNT)*W3'
                WRITE(LUFILE,'(A)')'     !  12 < WVAL <= (2J+36) '
                WRITE(LUFILE,'(A)')'     ELSE IF (WVAL.LE.D2JP36) THEN'
                WRITE(LUFILE,'(A)')'      REXPW = D05*EXP(-WVAL)'
                WRITE(LUFILE,'(A)')'      RWVAL = D1/WVAL'
                WRITE(LUFILE,'(A)')'      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))'
                WRITE(LUFILE,'(A)')'      RJ000 = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
                WRITE(LUFILE,'(A)')'     !  (2J+36) < WVAL '
                WRITE(LUFILE,'(A)')'     ELSE'
                WRITE(LUFILE,'(A)')'      RJ000 = SQRT(PID4/WVAL)'
                WRITE(LUFILE,'(A)')'     ENDIF'
                IF(Gen)THEN
                   WRITE(LUFILE,'(A,A,A)')'     AUXarray(',PrimLabel(1:iPrim),')=integralPrefactor(iPrimQ,iPrimP)*&'
                   WRITE(LUFILE,'(A)')'          & QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)*RJ000'
                ELSEIF(Seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'     AUXarray(iP)=integralPrefactor(1)*QpreExpFac(1)*PpreExpFac(iAtomA,iAtomB)*RJ000'
                ELSE
                   WRITE(LUFILE,'(A,A,A,A,A)')'     AUXarray(',PrimLabel(1:iPrim),')=AUXarray(',PrimLabel(1:iPrim),') + integralPrefactor(iPrimQ,iPrimP)*&'
                   WRITE(LUFILE,'(A)')'          & QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)*RJ000'
                ENDIF
                IF(Collapse)THEN
                   call PrintCollapseLoopEnd(Gen,SegQ,SegP,Seg,seg1prim,LUFILE)
                ELSE
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    enddo'
                      WRITE(LUFILE,'(A)')'   enddo'
                   ENDIF
                   WRITE(LUFILE,'(A)')'  enddo'
                ENDIF
!                IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END PARALLEL DO'
                IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END DO'
                WRITE(LUFILE,'(3A)')'end subroutine VerticalRecurrence'//ARCSTRING,SegLabel(1:iSegLabel),'0'
             endif
             !========================================================================================================
             !    VerticalRecurrence 1 
             !========================================================================================================
             JMAX = 1
             WRITE(LUFILE,'(A)')''
             WRITE(LUFILE,'(5A)')'subroutine VerticalRecurrence'//ARCSTRING,SegLabel(1:iSegLabel),'1',centerString,'(nPassP,nPrimP,nPrimQ,&'
             IF(center.EQ.1)THEN
                WRITE(LUFILE,'(A)')'         & reducedExponents,TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&'
             ELSEIF(center.EQ.2)THEN
                WRITE(LUFILE,'(A)')'         & reducedExponents,TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&'
             ELSEIF(center.EQ.3)THEN
                WRITE(LUFILE,'(A)')'         & reducedExponents,TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&'
             ELSEIF(center.EQ.4)THEN
                WRITE(LUFILE,'(A)')'         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&'
             ENDIF
             WRITE(LUFILE,'(A)')'         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&'
             WRITE(LUFILE,'(A)')'         & PpreExpFac,QpreExpFac,AUXarray)'
             WRITE(LUFILE,'(A)')'  implicit none'
             WRITE(LUFILE,'(A)')'  integer,intent(in) :: nPassP,nPrimP,nPrimQ'
             WRITE(LUFILE,'(A)')'  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB'
             WRITE(LUFILE,'(A)')'  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)'
             WRITE(LUFILE,'(A)')'  REAL(REALK),intent(in) :: TABFJW(0:4,0:1200)'
             IF(.NOT.Seg1Prim)THEN
                IF(center.LE.2)THEN
                   WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)'
                ELSE
                   WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)'
                ENDIF
             ELSE
                IF(center.LE.2)THEN
                   WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Pexp(1)'
                ELSE
                   WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Qexp(1)'
                ENDIF
             ENDIF
             IF(.NOT.Seg1Prim)THEN
                WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)'
                WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)'
             ELSE
                WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3),integralPrefactor(1),QpreExpFac(1),PpreExpFac(nAtomsA,nAtomsB)'
             ENDIF
             IF(center.LE.2)THEN
                WRITE(LUFILE,'(A,A1,A,A1,A1)')'  real(realk),intent(in) :: ',centerString,'center(3,nAtoms',centerString,')'
             ELSE
                WRITE(LUFILE,'(A,A1,A)')'  real(realk),intent(in) :: ',centerString,'center(3)'
             ENDIF
             WRITE(LUFILE,'(A,A,A)')'  real(realk),intent(inout) :: AUXarray(4,',nPrimLabel(1:nPrim),')'
             WRITE(LUFILE,'(A)')'  !local variables'
             IF(.NOT.Seg1Prim)THEN
                WRITE(LUFILE,'(A)')'  Integer :: iP,iPrimQ,iPrimP,iPassP,ipnt,iAtomA,iAtomB'
             ELSE
                WRITE(LUFILE,'(A)')'  Integer :: iP,iPassP,ipnt,iAtomA,iAtomB'
             ENDIF
             IF(COLLAPSE.AND.Seg)THEN
                WRITE(LUFILE,'(A)')'  Integer :: iPrimQP'
             ENDIF
             IF(center.EQ.1)WRITE(LUFILE,'(A)')'  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa'
             IF(center.EQ.2)WRITE(LUFILE,'(A)')'  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb'
             IF(center.EQ.3)WRITE(LUFILE,'(A)')'  real(realk) :: Xqc,Yqc,Zqc'
             IF(center.EQ.4)WRITE(LUFILE,'(A)')'  real(realk) :: Xqd,Yqd,Zqd'
             IF(center.LE.2)THEN
                WRITE(LUFILE,'(A)')'  real(realk) :: mPX,mPY,mPZ,invexpP,alphaP,RJ000(0:1)'
             else
                WRITE(LUFILE,'(A)')'  real(realk) :: mPX,mPY,mPZ,invexpQ,alphaQ,RJ000(0:1)'
             endif
             WRITE(LUFILE,'(A)')'  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq'
             WRITE(LUFILE,'(A)')'  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL' 
             WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk'
             WRITE(LUFILE,'(A)')'  real(realk),parameter :: D2=2.0E0_realk'
             WRITE(LUFILE,'(A,ES24.16,A)')'  REAL(REALK),PARAMETER :: D2JP36=',2.0d0 + 36.0d0,'_realk'
             WRITE(LUFILE,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D03333=1.0E0_realk/3.0E0_realk'
             WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk'
             WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk'
             WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: SMALL = 1E-15_realk,D12 = 12.0E0_realk'
             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk'
             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk'
             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk'
             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk'
             WRITE(LUFILE,'(A)')'  Real(realk), parameter :: PI=3.14159265358979323846E0_realk'
             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk'
             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2'
             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI'
             IF(center.EQ.1)WRITE(LUFILE,'(A)')'  !ThetaAux(n,1,0,0) = Xpa*ThetaAux(n,0,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,0,0,0)'
             IF(center.EQ.2)WRITE(LUFILE,'(A)')'  !ThetaAux(n,1,0,0) = Xpb*ThetaAux(n,0,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,0,0,0)'
             IF(center.EQ.3)WRITE(LUFILE,'(A)')'  !ThetaAux(n,1,0,0) = Xqc*ThetaAux(n,0,0,0) + (-alpha/q*Xpq)*ThetaAux(n+1,0,0,0)'
             IF(center.EQ.4)WRITE(LUFILE,'(A)')'  !ThetaAux(n,1,0,0) = Xqd*ThetaAux(n,0,0,0) + (-alpha/q*Xpq)*ThetaAux(n+1,0,0,0)'
             WRITE(LUFILE,'(A)')'  !i = 0 last 2 term vanish'
             WRITE(LUFILE,'(A)')'  !We include scaling of RJ000 '
             IF(Collapse)THEN
                call PrintCollapseInitLoop(Gen,SegQ,SegP,Seg,seg1prim,LUFILE,JMAX,nTUV,DoOpenMP)
             ENDIF
             IF(DoOpenMP.OR.DoOpenACC)THEN
                call PrintOpenMP(Gen,SegQ,SegP,Seg,seg1prim,LUFILE,JMAX,Collapse,Center,centerstring,DoOpenMP,DoOpenACC)
             ENDIF
             IF(Collapse)THEN
                call PrintCollapseLoopStart(Gen,SegQ,SegP,Seg,seg1prim,LUFILE)
                WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
                WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
                IF(center.EQ.1)THEN
                   WRITE(LUFILE,'(A)')'   Ax = -Acenter(1,iAtomA)'
                   WRITE(LUFILE,'(A)')'   Ay = -Acenter(2,iAtomA)'
                   WRITE(LUFILE,'(A)')'   Az = -Acenter(3,iAtomA)'
                ELSEIF(center.eq.2)THEN
                   WRITE(LUFILE,'(A)')'   Bx = -Bcenter(1,iAtomB)'
                   WRITE(LUFILE,'(A)')'   By = -Bcenter(2,iAtomB)'
                   WRITE(LUFILE,'(A)')'   Bz = -Bcenter(3,iAtomB)'
                ENDIF
                IF(seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'    mPX = -Pcent(1,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'    mPY = -Pcent(2,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'    mPZ = -Pcent(3,iAtomA,iAtomB)'
                ELSE
                   WRITE(LUFILE,'(A)')'    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)'
                ENDIF
                IF(center.LE.2)THEN
                   IF(seg1Prim)THEN 
                      WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(1)'
                   ELSE
                      WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(iPrimP)'
                   ENDIF
                ENDIF
                IF(center.GT.2)THEN
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    invexpQ = D1/Qexp(iPrimQ)'
                   ELSE
                      WRITE(LUFILE,'(A)')'    invexpQ = D1/Qexp(1)'
                   ENDIF
                ENDIF
                IF(center.EQ.1)THEN
                   IF(seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    Xpa = Pcent(1,iAtomA,iAtomB) + Ax'
                      WRITE(LUFILE,'(A)')'    Ypa = Pcent(2,iAtomA,iAtomB) + Ay'
                      WRITE(LUFILE,'(A)')'    Zpa = Pcent(3,iAtomA,iAtomB) + Az'
                   ELSE
                      WRITE(LUFILE,'(A)')'    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax'
                      WRITE(LUFILE,'(A)')'    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay'
                      WRITE(LUFILE,'(A)')'    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az'
                   ENDIF
                ELSEIF(center.EQ.2)THEN
                   IF(seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    Xpb = Pcent(1,iAtomA,iAtomB) + Bx'
                      WRITE(LUFILE,'(A)')'    Ypb = Pcent(2,iAtomA,iAtomB) + By'
                      WRITE(LUFILE,'(A)')'    Zpb = Pcent(3,iAtomA,iAtomB) + Bz'
                   ELSE
                      WRITE(LUFILE,'(A)')'    Xpb = Pcent(1,iPrimP,iAtomA,iAtomB) + Bx'
                      WRITE(LUFILE,'(A)')'    Ypb = Pcent(2,iPrimP,iAtomA,iAtomB) + By'
                      WRITE(LUFILE,'(A)')'    Zpb = Pcent(3,iPrimP,iAtomA,iAtomB) + Bz'
                   ENDIF
                ENDIF
                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'   Xpq = mPX + Qcent(1,iPrimQ)'
                   WRITE(LUFILE,'(A)')'   Ypq = mPY + Qcent(2,iPrimQ)'
                   WRITE(LUFILE,'(A)')'   Zpq = mPZ + Qcent(3,iPrimQ)'
                ELSE
                   WRITE(LUFILE,'(A)')'   Xpq = mPX + Qcent(1)'
                   WRITE(LUFILE,'(A)')'   Ypq = mPY + Qcent(2)'
                   WRITE(LUFILE,'(A)')'   Zpq = mPZ + Qcent(3)'
                ENDIF
             ELSE
                WRITE(LUFILE,'(A)')'  DO iPassP = 1,nPassP'
                WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
                WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
                IF(center.EQ.1)THEN
                   WRITE(LUFILE,'(A)')'   Ax = -Acenter(1,iAtomA)'
                   WRITE(LUFILE,'(A)')'   Ay = -Acenter(2,iAtomA)'
                   WRITE(LUFILE,'(A)')'   Az = -Acenter(3,iAtomA)'
                ELSEIF(center.eq.2)THEN
                   WRITE(LUFILE,'(A)')'   Bx = -Bcenter(1,iAtomB)'
                   WRITE(LUFILE,'(A)')'   By = -Bcenter(2,iAtomB)'
                   WRITE(LUFILE,'(A)')'   Bz = -Bcenter(3,iAtomB)'
                ENDIF
                IF(SegP)WRITE(LUFILE,'(A)')'   DO iPrimQ=1, nPrimQ'
                IF(Seg.OR.SegP)THEN
                   WRITE(LUFILE,'(A,A,A)')'   AUXarray(1,',PrimLabel(1:iPrim),')=0.0E0_realk'
                   WRITE(LUFILE,'(A,A,A)')'   AUXarray(2,',PrimLabel(1:iPrim),')=0.0E0_realk'
                   WRITE(LUFILE,'(A,A,A)')'   AUXarray(3,',PrimLabel(1:iPrim),')=0.0E0_realk'
                   WRITE(LUFILE,'(A,A,A)')'   AUXarray(4,',PrimLabel(1:iPrim),')=0.0E0_realk'
                ENDIF
                IF(SegP)WRITE(LUFILE,'(A)')'   ENDDO'
                IF(.NOT.seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'   DO iPrimP=1, nPrimP'
                ENDIF
                IF(SegQ)THEN
                   WRITE(LUFILE,'(A)')'    AUXarray(1,iPrimP,iPassP)=0.0E0_realk'
                   WRITE(LUFILE,'(A)')'    AUXarray(2,iPrimP,iPassP)=0.0E0_realk'
                   WRITE(LUFILE,'(A)')'    AUXarray(3,iPrimP,iPassP)=0.0E0_realk'
                   WRITE(LUFILE,'(A)')'    AUXarray(4,iPrimP,iPassP)=0.0E0_realk'
                ENDIF
                IF(seg1Prim)THEN
                   !          WRITE(LUFILE,'(A)')'    Pexpfac = PpreExpFac(iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'    mPX = -Pcent(1,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'    mPY = -Pcent(2,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'    mPZ = -Pcent(3,iAtomA,iAtomB)'
                ELSE
!                   WRITE(LUFILE,'(A)')'    Pexpfac = PpreExpFac(iPrimP,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)'
                ENDIF
                IF(center.LE.2)THEN
                   IF(seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(1)'
                   ELSE
                      WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(iPrimP)'
                   ENDIF
                ENDIF
                IF(center.GT.2)THEN
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    invexpQ = D1/Qexp(iPrimQ)'
                   ELSE
                      WRITE(LUFILE,'(A)')'    invexpQ = D1/Qexp(1)'
                   ENDIF
                ENDIF
                IF(center.EQ.1)THEN
                   IF(seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    Xpa = Pcent(1,iAtomA,iAtomB) + Ax'
                      WRITE(LUFILE,'(A)')'    Ypa = Pcent(2,iAtomA,iAtomB) + Ay'
                      WRITE(LUFILE,'(A)')'    Zpa = Pcent(3,iAtomA,iAtomB) + Az'
                   ELSE
                      WRITE(LUFILE,'(A)')'    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax'
                      WRITE(LUFILE,'(A)')'    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay'
                      WRITE(LUFILE,'(A)')'    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az'
                   ENDIF
                ELSEIF(center.EQ.2)THEN
                   IF(seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    Xpb = Pcent(1,iAtomA,iAtomB) + Bx'
                      WRITE(LUFILE,'(A)')'    Ypb = Pcent(2,iAtomA,iAtomB) + By'
                      WRITE(LUFILE,'(A)')'    Zpb = Pcent(3,iAtomA,iAtomB) + Bz'
                   ELSE
                      WRITE(LUFILE,'(A)')'    Xpb = Pcent(1,iPrimP,iAtomA,iAtomB) + Bx'
                      WRITE(LUFILE,'(A)')'    Ypb = Pcent(2,iPrimP,iAtomA,iAtomB) + By'
                      WRITE(LUFILE,'(A)')'    Zpb = Pcent(3,iPrimP,iAtomA,iAtomB) + Bz'
                   ENDIF
                ENDIF
                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'    DO iPrimQ=1, nPrimQ'
                   WRITE(LUFILE,'(A)')'     Xpq = mPX + Qcent(1,iPrimQ)'
                   WRITE(LUFILE,'(A)')'     Ypq = mPY + Qcent(2,iPrimQ)'
                   WRITE(LUFILE,'(A)')'     Zpq = mPZ + Qcent(3,iPrimQ)'
                ELSE
                   WRITE(LUFILE,'(A)')'     Xpq = mPX + Qcent(1)'
                   WRITE(LUFILE,'(A)')'     Ypq = mPY + Qcent(2)'
                   WRITE(LUFILE,'(A)')'     Zpq = mPZ + Qcent(3)'
                ENDIF
             ENDIF !collapse
             IF(center.EQ.3)THEN
                IF(seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'     Xqc = Qcent(1) - Ccenter(1)'
                   WRITE(LUFILE,'(A)')'     Yqc = Qcent(2) - Ccenter(2)'
                   WRITE(LUFILE,'(A)')'     Zqc = Qcent(3) - Ccenter(3)'
                ELSE
                   WRITE(LUFILE,'(A)')'     Xqc = Qcent(1,iPrimQ) - Ccenter(1)'
                   WRITE(LUFILE,'(A)')'     Yqc = Qcent(2,iPrimQ) - Ccenter(2)'
                   WRITE(LUFILE,'(A)')'     Zqc = Qcent(3,iPrimQ) - Ccenter(3)'
                ENDIF
             ELSEIF(center.EQ.4)THEN
                IF(seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'     Xqd = Qcent(1) - Dcenter(1)'
                   WRITE(LUFILE,'(A)')'     Yqd = Qcent(2) - Dcenter(2)'
                   WRITE(LUFILE,'(A)')'     Zqd = Qcent(3) - Dcenter(3)'
                ELSE
                   WRITE(LUFILE,'(A)')'     Xqd = Qcent(1,iPrimQ) - Dcenter(1)'
                   WRITE(LUFILE,'(A)')'     Yqd = Qcent(2,iPrimQ) - Dcenter(2)'
                   WRITE(LUFILE,'(A)')'     Zqd = Qcent(3,iPrimQ) - Dcenter(3)'
                ENDIF
             ENDIF
             IF(.NOT.Seg1Prim)THEN
                IF(center.LE.2)THEN
                   WRITE(LUFILE,'(A)')'     alphaP = reducedExponents(iPrimQ,iPrimP)*invexpP'    
                ELSE
                   WRITE(LUFILE,'(A)')'     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ'
                ENDIF
             ELSE
                IF(center.LE.2)THEN
                   WRITE(LUFILE,'(A)')'     alphaP = reducedExponents(1)*invexpP'    
                ELSE
                   WRITE(LUFILE,'(A)')'     alphaQ = -reducedExponents(1)*invexpQ'
                ENDIF
             ENDIF
             IF(center.LE.2)THEN
                WRITE(LUFILE,'(A)')'     alphaXpq = alphaP*Xpq'
                WRITE(LUFILE,'(A)')'     alphaYpq = alphaP*Ypq'
                WRITE(LUFILE,'(A)')'     alphaZpq = alphaP*Zpq'
             ELSE
                WRITE(LUFILE,'(A)')'     alphaXpq = alphaQ*Xpq'
                WRITE(LUFILE,'(A)')'     alphaYpq = alphaQ*Ypq'
                WRITE(LUFILE,'(A)')'     alphaZpq = alphaQ*Zpq'
             ENDIF
             WRITE(LUFILE,'(A)')'     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq'
             IF(.NOT.Seg1Prim)THEN
                WRITE(LUFILE,'(A)')'     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance'
             ELSE
                WRITE(LUFILE,'(A)')'     WVAL = reducedExponents(1)*squaredDistance'
             ENDIF
             !       WRITE(LUFILE,'(A)')'     !  0 < WVAL < 0.000001'
             !       WRITE(LUFILE,'(A)')'!     IF (ABS(WVAL) .LT. SMALL) THEN'     
             !       WRITE(LUFILE,'(A)')'!      RJ000(0) = D1'
             !       WRITE(LUFILE,'(A)')'!      RJ000(1)= D03333 !THE BOYS FUNCTION FOR ZERO ARGUMENT'
             !       WRITE(LUFILE,'(A)')'!     !  0 < WVAL < 12 '
             WRITE(LUFILE,'(A)')'     IF (WVAL .LT. D12) THEN'
             WRITE(LUFILE,'(A)')'      IPNT = NINT(D100*WVAL)'
             WRITE(LUFILE,'(A)')'      WDIFF = WVAL - TENTH*IPNT'
             WRITE(LUFILE,'(A)')'      W2    = WDIFF*WDIFF'
             WRITE(LUFILE,'(A)')'      W3    = W2*WDIFF'
             WRITE(LUFILE,'(A)')'      W2    = W2*D05'
             WRITE(LUFILE,'(A)')'      W3    = W3*COEF3'
             WRITE(LUFILE,'(A)')'      RJ000(0)=TABFJW(0,IPNT)-TABFJW(1,IPNT)*WDIFF+TABFJW(2,IPNT)*W2+TABFJW(3,IPNT)*W3'
             WRITE(LUFILE,'(A)')'      RJ000(1)=TABFJW(1,IPNT)-TABFJW(2,IPNT)*WDIFF+TABFJW(3,IPNT)*W2+TABFJW(4,IPNT)*W3'
             WRITE(LUFILE,'(A)')'     !  12 < WVAL <= (2J+36) '
             WRITE(LUFILE,'(A)')'     ELSE IF (WVAL.LE.D2JP36) THEN'
             WRITE(LUFILE,'(A)')'      REXPW = D05*EXP(-WVAL)'
             WRITE(LUFILE,'(A)')'      RWVAL = D1/WVAL'
             WRITE(LUFILE,'(A)')'      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))'
             WRITE(LUFILE,'(A)')'      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
             WRITE(LUFILE,'(A)')'      RJ000(1) = RWVAL*(D05*RJ000(0)-REXPW)'
             WRITE(LUFILE,'(A)')'     !  (2J+36) < WVAL '
             WRITE(LUFILE,'(A)')'     ELSE'
             WRITE(LUFILE,'(A)')'      RWVAL = PID4/WVAL'
             WRITE(LUFILE,'(A)')'      RJ000(0) = SQRT(RWVAL)'
             WRITE(LUFILE,'(A)')'      RJ000(1) = RWVAL*PID4I*D05*RJ000(0)'
             WRITE(LUFILE,'(A)')'     ENDIF'
             IF(.NOT.Seg1Prim)THEN
                WRITE(LUFILE,'(A)')'     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)'
             ELSE
                WRITE(LUFILE,'(A)')'     PREF = integralPrefactor(1)*QpreExpFac(1)*PpreExpFac(iAtomA,iAtomB)'
             ENDIF
             WRITE(LUFILE,'(A)')'     TMP1 = PREF*RJ000(0)'
             WRITE(LUFILE,'(A)')'     TMP2 = PREF*RJ000(1)'
             IF(Gen.OR.Seg1Prim) THEN
                WRITE(LUFILE,'(A,A,A)')'     AUXarray(1,',PrimLabel(1:iPrim),') = TMP1'
             ELSE
                WRITE(LUFILE,'(A,A,A,A,A)')'     AUXarray(1,',PrimLabel(1:iPrim),') = AUXarray(1,',PrimLabel(1:iPrim),') + TMP1'
             ENDIF
             IF(center.EQ.1)THEN
                Xdir='Xpa';Ydir='Ypa';Zdir='Zpa' 
             ENDIF
             IF(center.EQ.2)THEN
                Xdir='Xpb';Ydir='Ypb';Zdir='Zpb' 
             ENDIF
             IF(center.EQ.3)THEN
                Xdir='Xqc';Ydir='Yqc';Zdir='Zqc' 
             ENDIF
             IF(center.EQ.4)THEN
                Xdir='Xqd';Ydir='Yqd';Zdir='Zqd' 
             ENDIF
             IF(Gen.OR.Seg1prim)THEN
                WRITE(LUFILE,'(A,A,A,A,A)')'     AUXarray(2,',PrimLabel(1:iPrim),') = ',Xdir,'*TMP1 + alphaXpq*TMP2'
                WRITE(LUFILE,'(A,A,A,A,A)')'     AUXarray(3,',PrimLabel(1:iPrim),') = ',Ydir,'*TMP1 + alphaYpq*TMP2'
                WRITE(LUFILE,'(A,A,A,A,A)')'     AUXarray(4,',PrimLabel(1:iPrim),') = ',Zdir,'*TMP1 + alphaZpq*TMP2'
             ELSE
                WRITE(LUFILE,'(A,A,A,A,A,A,A)')'     AUXarray(2,',PrimLabel(1:iPrim),') = AUXarray(2,',PrimLabel(1:iPrim),') + ',Xdir,'*TMP1 + alphaXpq*TMP2'
                WRITE(LUFILE,'(A,A,A,A,A,A,A)')'     AUXarray(3,',PrimLabel(1:iPrim),') = AUXarray(3,',PrimLabel(1:iPrim),') + ',Ydir,'*TMP1 + alphaYpq*TMP2'
                WRITE(LUFILE,'(A,A,A,A,A,A,A)')'     AUXarray(4,',PrimLabel(1:iPrim),') = AUXarray(4,',PrimLabel(1:iPrim),') + ',Zdir,'*TMP1 + alphaZpq*TMP2'
             ENDIF
             IF(Collapse)THEN
                call PrintCollapseLoopEnd(Gen,SegQ,SegP,Seg,seg1prim,LUFILE)
             ELSE
                IF(.NOT.seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'    enddo'
                   WRITE(LUFILE,'(A)')'   enddo'
                ENDIF
                WRITE(LUFILE,'(A)')'  enddo'
             ENDIF
!             IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END PARALLEL DO'
             IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END DO'
             WRITE(LUFILE,'(4A)')'end subroutine VerticalRecurrence'//ARCSTRING,SegLabel(1:iSegLabel),'1',centerString

             !============================================================================================================
             !         VerticalRecurrence GENERAL
             !============================================================================================================

             DO JMAX=2,MaxAngmomQP                                  !Gen   !Seg
                IF(center.EQ.1.AND.JMAX.GT.MaxAngmomQP)CYCLE        !DDDD  DDSS
                IF(Gen.OR.Seg1Prim)THEN
                   IF(center.EQ.3.AND.JMAX.GT.MaxAngmomQP-1)CYCLE  !PDDD
                   IF(center.EQ.2.AND.JMAX.GT.MaxAngmomQP-2)CYCLE   !PDPD
                   IF(center.EQ.4.AND.JMAX.GT.MaxAngmomQP-3)CYCLE   !PPPD
                ELSE
                   IF(center.EQ.3.AND.JMAX.GT.MaxAngmomQP)CYCLE    !SSDD
                   IF(center.EQ.2.AND.JMAX.GT.MaxAngmomQP-1)CYCLE   !PDSS 
                   IF(center.EQ.4.AND.JMAX.GT.MaxAngmomQP-1)CYCLE   !SSPD                   
                ENDIF
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
                WRITE(LUFILE,'(A)')''
                WRITE(LUFILE,'(A,A,I1,A,A)')'subroutine VerticalRecurrence'//ARCSTRING,SegLabel(1:iSegLabel),JMAX,centerString,'(nPassP,nPrimP,nPrimQ,&'
                IF(center.EQ.1)THEN
                   WRITE(LUFILE,'(A)')'         & reducedExponents,RJ000Array,Pexp,Acenter,Pcent,Qcent,integralPrefactor,&'
                ELSEIF(center.EQ.2)THEN
                   WRITE(LUFILE,'(A)')'         & reducedExponents,RJ000Array,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,&'
                ELSEIF(center.EQ.3)THEN
                   WRITE(LUFILE,'(A)')'         & reducedExponents,RJ000Array,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&'
                ELSEIF(center.EQ.4)THEN
                   WRITE(LUFILE,'(A)')'         & reducedExponents,RJ000Array,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&'
                ENDIF
                WRITE(LUFILE,'(A)')'         & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,AUXarray)'
                WRITE(LUFILE,'(A)')'  implicit none'
                WRITE(LUFILE,'(A)')'  integer,intent(in) :: nPassP,nPrimP,nPrimQ'
                WRITE(LUFILE,'(A)')'  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB'
                WRITE(LUFILE,'(A)')'  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)'
                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(A,I2,A)')'  REAL(REALK),intent(in) :: RJ000Array(0:',JMAX,',nPrimQ,nPrimP,nPassP)'
                   IF(center.LE.2)THEN
                      WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)'
                   ELSE
                      WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)'
                   ENDIF
                   WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)'
                   WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ),PpreExpFac(nPrimP,nAtomsA,nAtomsB)'
                ELSE
                   WRITE(LUFILE,'(A,I2,A)')'  REAL(REALK),intent(in) :: RJ000Array(0:',JMAX,',nPassP)'
                   IF(center.LE.2)THEN
                      WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Pexp(1)'
                   ELSE
                      WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Qexp(1)'
                   ENDIF
                   WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3),integralPrefactor(1),QpreExpFac(1),PpreExpFac(nAtomsA,nAtomsB)'
                ENDIF
                IF(center.EQ.1)WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Acenter(3,nAtomsA)'
                IF(center.EQ.2)WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Bcenter(3,nAtomsB)'
                IF(center.EQ.3)WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Ccenter(3)'
                IF(center.EQ.4)WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Dcenter(3)'
                IF(COLLAPSE)THEN
                   IF(Gen) WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPrimQ*nPrimP*nPassP)'
                   IF(SegQ)WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPrimP*nPassP)'
                   IF(SegP)WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPrimQ*nPassP)'
                   IF(Seg) WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPassP)'
                   IF(Seg1Prim)WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPassP)'
                ELSE
                   IF(Gen) WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPrimQ,nPrimP,nPassP)'
                   IF(SegQ)WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPrimP,nPassP)'
                   IF(SegP)WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPrimQ,nPassP)'
                   IF(Seg) WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPassP)'
                   IF(Seg1Prim)WRITE(LUFILE,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPassP)'
                ENDIF
                WRITE(LUFILE,'(A)')'  !local variables'
                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'  integer :: iPassP,iPrimP,iPrimQ,ipnt,IP,iTUV,iAtomA,iAtomB'
                ELSE
                   WRITE(LUFILE,'(A)')'  integer :: iPassP,ipnt,IP,iTUV,iAtomA,iAtomB'
                ENDIF
                IF(COLLAPSE.AND.Seg)THEN
                   WRITE(LUFILE,'(A)')'  Integer :: iPrimQP'
                ENDIF
                IF(.NOT.Gen)THEN
                   WRITE(LUFILE,'(A,I5,A)')'  real(realk) :: TMPAUXarray(',nTUVprev,')'
                ENDIF
                IF(center.EQ.1)WRITE(LUFILE,'(A)')'  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa'
                IF(center.EQ.2)WRITE(LUFILE,'(A)')'  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb'
                IF(center.EQ.3)WRITE(LUFILE,'(A)')'  real(realk) :: Xqc,Yqc,Zqc'
                IF(center.EQ.4)WRITE(LUFILE,'(A)')'  real(realk) :: Xqd,Yqd,Zqd'
                IF(center.LE.2)THEN
                   WRITE(LUFILE,'(A,i2,A)')'  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP'
                ELSE
                   WRITE(LUFILE,'(A,i2,A)')'  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ'
                ENDIF
                WRITE(LUFILE,'(A,i4,A)')'  real(realk) :: TwoTerms(',MAX(1,nTUVprev2-nTUVprev3),')'
                WRITE(LUFILE,'(A)')'  real(realk) :: PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq'
                !             WRITE(LUFILE,'(A)')'  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL' 
                !             WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk'
                WRITE(LUFILE,'(A)')'  real(realk),parameter :: D2=2.0E0_realk,D05 =0.5E0_realk'
                !             WRITE(LUFILE,'(A,ES24.16,A)')'  REAL(REALK),PARAMETER :: D2JP36=',2.d0*JMAX + 36.d0,'_realk'
                WRITE(LUFILE,'(A)')'  real(realk),parameter :: D1=1.0E0_realk'
                !             WRITE(LUFILE,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D03333=1.0E0_realk/3.0E0_realk'
                !             WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk'
                !             WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk'
                !             WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: SMALL = 1E-15_realk,D12 = 12.0E0_realk'
                !             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk'
                !             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk'
                !             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk'
                !             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk'
                !             WRITE(LUFILE,'(A)')'  Real(realk), parameter :: PI=3.14159265358979323846E0_realk'
                !             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk'
                !             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2'
                !             WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI'
                C = JMAX+2
                DO JTMP=1,JMAX
                   C = C-1
                   nTUVTMPprev=(JTMP-1)*(JTMP)*(JTMP+1)/6
                   nTUVTMP=(JTMP)*(JTMP+1)*(JTMP+2)/6
                   if(JTMP.LT.10)THEN
                      WRITE(LUFILE,'(A,I1,A,I3,A,I3,A,I1,A)')'  real(realk) :: TMParray',JTMP,'(',nTUVTMPprev+1,':',nTUVTMP,',2:',C,')'
                   else
                      WRITE(LUFILE,'(A,I2,A,I3,A,I3,A,I1,A)')'  real(realk) :: TMParray',JTMP,'(',nTUVTMPprev+1,':',nTUVTMP,',2:',C,')'

                   endif
                ENDDO
                IF(center.Eq.1)then
                   WRITE(LUFILE,'(A)')'  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)'
                   WRITE(LUFILE,'(A)')'  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))'
                elseif(center.eq.2)then
                   WRITE(LUFILE,'(A)')'  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)'
                   WRITE(LUFILE,'(A)')'  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))'
                elseif(center.eq.3)then
                   WRITE(LUFILE,'(A)')'  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)'
                   WRITE(LUFILE,'(A)')'  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))'
                elseif(center.eq.4)then
                   WRITE(LUFILE,'(A)')'  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)'
                   WRITE(LUFILE,'(A)')'  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))'
                endif
                WRITE(LUFILE,'(A)')'  !We include scaling of RJ000 '
                IF(Collapse)THEN
                   call PrintCollapseInitLoop(Gen,SegQ,SegP,Seg,seg1prim,LUFILE,JMAX,nTUV,DoOpenMP)
                ENDIF
                IF(DoOpenMP.OR.DoOpenACC)THEN
                   call PrintOpenMP(Gen,SegQ,SegP,Seg,seg1prim,LUFILE,JMAX,Collapse,Center,centerstring,DoOpenMP,DoOpenACC)
                ENDIF
                ! ======================================================================
                !    iPassP Loop 
                ! ======================================================================
                IF(Collapse)THEN
                   call PrintCollapseLoopStart(Gen,SegQ,SegP,Seg,seg1prim,LUFILE)
                   WRITE(LUFILE,'(A)')'     iAtomA = iAtomApass(iPassP)'
                   WRITE(LUFILE,'(A)')'     iAtomB = iAtomBpass(iPassP)'
                   IF(center.EQ.1)THEN
                      WRITE(LUFILE,'(A)')'   Ax = -Acenter(1,iAtomA)'
                      WRITE(LUFILE,'(A)')'   Ay = -Acenter(2,iAtomA)'
                      WRITE(LUFILE,'(A)')'   Az = -Acenter(3,iAtomA)'
                   ELSEIF(center.eq.2)THEN
                      WRITE(LUFILE,'(A)')'   Bx = -Bcenter(1,iAtomB)'
                      WRITE(LUFILE,'(A)')'   By = -Bcenter(2,iAtomB)'
                      WRITE(LUFILE,'(A)')'   Bz = -Bcenter(3,iAtomB)'
                   ENDIF
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)'
                      WRITE(LUFILE,'(A)')'    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)'
                      WRITE(LUFILE,'(A)')'    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)'
                   ENDIF
                   IF(center.LE.2)THEN
                      IF(.NOT.Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(iPrimP)'
                         WRITE(LUFILE,'(A)')'    inv2expP = D05*invexpP'
                      ELSE
                         WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(1)'
                         WRITE(LUFILE,'(A)')'    inv2expP = D05*invexpP'
                      ENDIF                      
                   ENDIF
                   IF(center.EQ.1)THEN
                      IF(.NOT.Seg1Prim)THEN             
                         WRITE(LUFILE,'(A)')'    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax'
                         WRITE(LUFILE,'(A)')'    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay'
                         WRITE(LUFILE,'(A)')'    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az'
                      ELSE
                         WRITE(LUFILE,'(A)')'    Xpa = Pcent(1,iAtomA,iAtomB) + Ax'
                         WRITE(LUFILE,'(A)')'    Ypa = Pcent(2,iAtomA,iAtomB) + Ay'
                         WRITE(LUFILE,'(A)')'    Zpa = Pcent(3,iAtomA,iAtomB) + Az'
                      ENDIF
                   ENDIF
                   IF(center.EQ.2)THEN
                      IF(.NOT.Seg1Prim)THEN             
                         WRITE(LUFILE,'(A)')'    Xpb = Pcent(1,iPrimP,iAtomA,iAtomB) + Bx'
                         WRITE(LUFILE,'(A)')'    Ypb = Pcent(2,iPrimP,iAtomA,iAtomB) + By'
                         WRITE(LUFILE,'(A)')'    Zpb = Pcent(3,iPrimP,iAtomA,iAtomB) + Bz'
                      ELSE
                         WRITE(LUFILE,'(A)')'    Xpb = Pcent(1,iAtomA,iAtomB) + Bx'
                         WRITE(LUFILE,'(A)')'    Ypb = Pcent(2,iAtomA,iAtomB) + By'
                         WRITE(LUFILE,'(A)')'    Zpb = Pcent(3,iAtomA,iAtomB) + Bz'
                      ENDIF
                   ENDIF
                ELSE
                   WRITE(LUFILE,'(A)')'  DO iPassP = 1,nPassP'
                   WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
                   WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
                   IF(center.EQ.1)THEN
                      WRITE(LUFILE,'(A)')'   Ax = -Acenter(1,iAtomA)'
                      WRITE(LUFILE,'(A)')'   Ay = -Acenter(2,iAtomA)'
                      WRITE(LUFILE,'(A)')'   Az = -Acenter(3,iAtomA)'
                   ELSEIF(center.eq.2)THEN
                      WRITE(LUFILE,'(A)')'   Bx = -Bcenter(1,iAtomB)'
                      WRITE(LUFILE,'(A)')'   By = -Bcenter(2,iAtomB)'
                      WRITE(LUFILE,'(A)')'   Bz = -Bcenter(3,iAtomB)'
                   ENDIF
                   IF(Gen) WRITE(LUFILE,'(A)')'   iP = (iPassP-1)*nPrimQ*nPrimP'
                   IF(SegQ)WRITE(LUFILE,'(A)')'   iP = (iPassP-1)*nPrimP'
                   IF(SegP)THEN
                      WRITE(LUFILE,'(A)')'   DO iPrimQ=1, nPrimQ'
                      WRITE(LUFILE,'(A,i5)')'    DO iTUV=1,',nTUV
                      WRITE(LUFILE,'(A)')'     AUXarray(iTUV,iPrimQ,iPassP)=0.0E0_realk'
                      WRITE(LUFILE,'(A)')'    ENDDO'
                      WRITE(LUFILE,'(A)')'   ENDDO'
                      WRITE(LUFILE,'(A)')'   iP = (iPassP-1)*nPrimQ'
                   ENDIF
                   !seg
                   IF(Seg)THEN
                      WRITE(LUFILE,'(A)')'   iP = iPassP'
                      WRITE(LUFILE,'(A,i5)')'   DO iTUV=1,',nTUV
                      WRITE(LUFILE,'(A)')'    AUXarray(iTUV,iPassP)=0.0E0_realk'
                      WRITE(LUFILE,'(A)')'   ENDDO'
                   endif
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'   DO iPrimP=1, nPrimP'
                   ENDIF
                   IF(SegP)WRITE(LUFILE,'(A)')'    iP = (iPassP-1)*nPrimQ'
                   IF(SegQ)THEN
                      WRITE(LUFILE,'(A,i5)')'    DO iTUV=1,',nTUV
                      WRITE(LUFILE,'(A)')'     AUXarray(iTUV,iPrimP,iPassP)=0.0E0_realk'
                      WRITE(LUFILE,'(A)')'    ENDDO'
                   ENDIF
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)'
                      WRITE(LUFILE,'(A)')'    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)'
                      WRITE(LUFILE,'(A)')'    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)'
                   ENDIF
                   IF(center.LE.2)THEN
                      IF(.NOT.Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(iPrimP)'
                         WRITE(LUFILE,'(A)')'    inv2expP = D05*invexpP'
                      ELSE
                         WRITE(LUFILE,'(A)')'    invexpP = D1/Pexp(1)'
                         WRITE(LUFILE,'(A)')'    inv2expP = D05*invexpP'
                      ENDIF
                   ENDIF
                   IF(center.EQ.1)THEN
                      IF(.NOT.Seg1Prim)THEN             
                         WRITE(LUFILE,'(A)')'    Xpa = Pcent(1,iPrimP,iAtomA,iAtomB) + Ax'
                         WRITE(LUFILE,'(A)')'    Ypa = Pcent(2,iPrimP,iAtomA,iAtomB) + Ay'
                         WRITE(LUFILE,'(A)')'    Zpa = Pcent(3,iPrimP,iAtomA,iAtomB) + Az'
                      ELSE
                         WRITE(LUFILE,'(A)')'    Xpa = Pcent(1,iAtomA,iAtomB) + Ax'
                         WRITE(LUFILE,'(A)')'    Ypa = Pcent(2,iAtomA,iAtomB) + Ay'
                         WRITE(LUFILE,'(A)')'    Zpa = Pcent(3,iAtomA,iAtomB) + Az'
                      ENDIF
                   ENDIF
                   IF(center.EQ.2)THEN
                      IF(.NOT.Seg1Prim)THEN             
                         WRITE(LUFILE,'(A)')'    Xpb = Pcent(1,iPrimP,iAtomA,iAtomB) + Bx'
                         WRITE(LUFILE,'(A)')'    Ypb = Pcent(2,iPrimP,iAtomA,iAtomB) + By'
                         WRITE(LUFILE,'(A)')'    Zpb = Pcent(3,iPrimP,iAtomA,iAtomB) + Bz'
                      ELSE
                         WRITE(LUFILE,'(A)')'    Xpb = Pcent(1,iAtomA,iAtomB) + Bx'
                         WRITE(LUFILE,'(A)')'    Ypb = Pcent(2,iAtomA,iAtomB) + By'
                         WRITE(LUFILE,'(A)')'    Zpb = Pcent(3,iAtomA,iAtomB) + Bz'
                      ENDIF
                   ENDIF
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    DO iPrimQ=1, nPrimQ'
                   ENDIF
                   IF(Gen) WRITE(LUFILE,'(A)')'     iP = iP + 1'
                   IF(SegP)WRITE(LUFILE,'(A)')'   iP = iP + 1'
                ENDIF !collapse
                IF(center.gt.2)THEN
                   IF(Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'     invexpQ = D1/Qexp(1)'
                      WRITE(LUFILE,'(A)')'     inv2expQ = D05*invexpQ'
                   ELSE
                      WRITE(LUFILE,'(A)')'     invexpQ = D1/Qexp(iPrimQ)'
                      WRITE(LUFILE,'(A)')'     inv2expQ = D05*invexpQ'
                   ENDIF
                ENDIF
                !if Seg or seg1prim iP not used
                !if SegQ then IP should be increased after iPrimP only 
                !if Gen then IP should be increased after iPrimQ 
                !if SegP then IP should be increased after iPrimQ 
                IF(.NOT.Seg1Prim)THEN
                   WRITE(LUFILE,'(A)')'     Xpq = mPX + Qcent(1,iPrimQ)'
                   WRITE(LUFILE,'(A)')'     Ypq = mPY + Qcent(2,iPrimQ)'
                   WRITE(LUFILE,'(A)')'     Zpq = mPZ + Qcent(3,iPrimQ)'
                ELSE
                   WRITE(LUFILE,'(A)')'     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)'
                   WRITE(LUFILE,'(A)')'     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)'
                ENDIF
                IF(center.EQ.3)THEN
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'     Xqc = Qcent(1,iPrimQ) - Ccenter(1)'
                      WRITE(LUFILE,'(A)')'     Yqc = Qcent(2,iPrimQ) - Ccenter(2)'
                      WRITE(LUFILE,'(A)')'     Zqc = Qcent(3,iPrimQ) - Ccenter(3)'
                   ELSE
                      WRITE(LUFILE,'(A)')'     Xqc = Qcent(1) - Ccenter(1)'
                      WRITE(LUFILE,'(A)')'     Yqc = Qcent(2) - Ccenter(2)'
                      WRITE(LUFILE,'(A)')'     Zqc = Qcent(3) - Ccenter(3)'
                   ENDIF
                ENDIF
                IF(center.EQ.4)THEN
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'     Xqd = Qcent(1,iPrimQ) - Dcenter(1)'
                      WRITE(LUFILE,'(A)')'     Yqd = Qcent(2,iPrimQ) - Dcenter(2)'
                      WRITE(LUFILE,'(A)')'     Zqd = Qcent(3,iPrimQ) - Dcenter(3)'
                   ELSE
                      WRITE(LUFILE,'(A)')'     Xqd = Qcent(1) - Dcenter(1)'
                      WRITE(LUFILE,'(A)')'     Yqd = Qcent(2) - Dcenter(2)'
                      WRITE(LUFILE,'(A)')'     Zqd = Qcent(3) - Dcenter(3)'
                   ENDIF
                ENDIF
                IF(.NOT.Seg1Prim)THEN
                   IF(center.LE.2)THEN
                      WRITE(LUFILE,'(A)')'     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP'    
                   ELSE
                      WRITE(LUFILE,'(A)')'     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ'
                   ENDIF
                ELSE
                   IF(center.LE.2)THEN
                      WRITE(LUFILE,'(A)')'     alphaP = -reducedExponents(1)*invexpP'    
                   ELSE
                      WRITE(LUFILE,'(A)')'     alphaQ = -reducedExponents(1)*invexpQ'
                   ENDIF
                ENDIF
                IF(center.LE.2)THEN
                   WRITE(LUFILE,'(A)')'     alphaXpq = -alphaP*Xpq'
                   WRITE(LUFILE,'(A)')'     alphaYpq = -alphaP*Ypq'
                   WRITE(LUFILE,'(A)')'     alphaZpq = -alphaP*Zpq'
                ENDIF
!                IF(.NOT.Seg1Prim)THEN
!                   IF(center.GT.2)THEN
!                      WRITE(LUFILE,'(A)')'     inv2expQ = D05*invexpQ'
!                   ENDIF
!                ENDIF
                IF(center.GT.2)THEN
                   WRITE(LUFILE,'(A)')'     alphaXpq = alphaQ*Xpq'
                   WRITE(LUFILE,'(A)')'     alphaYpq = alphaQ*Ypq'
                   WRITE(LUFILE,'(A)')'     alphaZpq = alphaQ*Zpq'
                ENDIF
                IF(.NOT.Seg1Prim)THEN             
                   WRITE(LUFILE,'(A)')'     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ)*PpreExpFac(iPrimP,iAtomA,iAtomB)'
                ElSE
                   WRITE(LUFILE,'(A)')'     PREF = integralPrefactor(1)*QpreExpFac(1)*PpreExpFac(iAtomA,iAtomB)'
                ENDIF
                IF(COLLAPSE)THEN
                   IF(Gen) WRITE(LUFILE,'(A)')'     Auxarray(1,IP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)'
                ELSE
                   IF(Gen) WRITE(LUFILE,'(A)')'     Auxarray(1,iPrimQ,iPrimP,iPassP) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)'
                ENDIF
                IF(SegQ)WRITE(LUFILE,'(A)')'     TMPAuxarray(1) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)'
                IF(SegP)WRITE(LUFILE,'(A)')'     TMPAuxarray(1) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)'
                IF(Seg) WRITE(LUFILE,'(A)')'     TMPAuxarray(1) = PREF*RJ000Array(0,iPrimQ,iPrimP,iPassP)'
                IF(Seg1Prim) WRITE(LUFILE,'(A)')'     TMPAuxarray(1) = PREF*RJ000Array(0,iPassP)'
                DO J=2,JMAX+1
                   IF(.NOT.seg1Prim)THEN
                      WRITE(LUFILE,'(A,I2,A,I2,A)')'     TMParray1(1,',J,') = PREF*RJ000Array(',J-1,',iPrimQ,iPrimP,iPassP)'
                   ELSE
                      WRITE(LUFILE,'(A,I2,A,I2,A)')'     TMParray1(1,',J,') = PREF*RJ000Array(',J-1,',iPassP)'
                   ENDIF
                ENDDO

                SegWTUV =.FALSE.
                IF(Seg) SegWTUV =.TRUE.
                IF(Seg1Prim) SegWTUV =.TRUE.

                !Gen,SegQ,SegP,Seg
                allocate(CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1))
                CREATED  = .FALSE.
                CREATED(0,0,0) = .TRUE.
                DO JTMP=1,JMAX
                 DO J = 1,JMAX-JTMP+1
                  !determine twoterms
                  nTUVLIST = (JTMP+1)*(JTMP+2)*(JTMP+3)/6   
                  allocate(TwoTermTUVLIST(nTUVLIST))
                  allocate(TwoTermsUsed(nTUVLIST))
                  call TwoTerms1(J,JTMP,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
                  DO Tp=JTMP,0,-1       
                   DO Up=JTMP-Tp,0,-1
                    Vp=JTMP-Tp-Up                
                    iTUV = TUVINDEX(Tp,Up,Vp)
                    TREC = CREATED(Tp-1,Up,Vp); UREC = CREATED(Tp,Up-1,Vp);VREC = CREATED(Tp,Up,Vp-1)
                    N=0
                    IF(TREC)N=N+1
                    IF(UREC)N=N+1
                    IF(VREC)N=N+1
                    IF(N.EQ.1)THEN
                     !only one possible way to construct it
                     IF(TREC)THEN
                      call TRECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
                     ELSEIF(UREC)THEN
                      call URECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
                     ELSEIF(VREC)THEN
                      call VRECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
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
                           CALL TRECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
                        ELSEIF(.NOT.(UREC.AND.UREC2).AND.UREC)THEN
                                     !                         print*,'!B URECURRENCE iTUV',iTUV
                                     CALL URECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
                                  ELSEIF(.NOT.(VREC.AND.VREC2).AND.VREC)THEN
                                     !                         print*,'!B VRECURRENCE iTUV',iTUV
                                     CALL VRECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
                                  ENDIF
                               ELSE
                                  !chose one of the possibilities
                                  IF(TREC)THEN
                                     !                         print*,'!C TRECURRENCE iTUV',iTUV
                                     CALL TRECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
                                  ELSEIF(UREC)THEN
                                     !                         print*,'!C URECURRENCE iTUV',iTUV
                                     call URECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
                                  ELSEIF(VREC)THEN
                                     !                         print*,'!D VRECURRENCE iTUV',iTUV
                                     call VRECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
                                  ENDIF
                               ENDIF
                            ENDIF
                            CREATED(Tp,Up,Vp) = .TRUE.
                         ENDDO
                      ENDDO

                      call WriteTwoTerms(J,JTMP,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,&
                           & TwoTermTUVLIST,JTMP,TwoTermsUsed,LUFILE,Gen,SegQ,SegP,Segwtuv,&
                           & seg1prim,center,COLLAPSE,PrimLabel,iPrim)


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
                                  call TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,LUFILE,Gen,SegQ,SegP,Segwtuv,seg1prim,nTUVprev,center,PrimLabel,iPrim)
                               ELSEIF(UREC)THEN
                                  !                      print*,'!A URECURRENCE iTUV',iTUV
                                  call URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,LUFILE,Gen,SegQ,SegP,Segwtuv,seg1prim,nTUVprev,center,PrimLabel,iPrim)
                               ELSEIF(VREC)THEN
                                  !                      print*,'!A VRECURRENCE iTUV',iTUV
                                  call VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,LUFILE,Gen,SegQ,SegP,Segwtuv,seg1prim,nTUVprev,center,PrimLabel,iPrim)
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
                                     CALL TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,LUFILE,Gen,SegQ,SegP,Segwtuv,seg1prim,nTUVprev,center,PrimLabel,iPrim)
                                  ELSEIF(.NOT.(UREC.AND.UREC2).AND.UREC)THEN
                                     !                         print*,'!B URECURRENCE iTUV',iTUV
                                     CALL URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,LUFILE,Gen,SegQ,SegP,Segwtuv,seg1prim,nTUVprev,center,PrimLabel,iPrim)
                                  ELSEIF(.NOT.(VREC.AND.VREC2).AND.VREC)THEN
                                     !                         print*,'!B VRECURRENCE iTUV',iTUV
                                     CALL VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,LUFILE,Gen,SegQ,SegP,Segwtuv,seg1prim,nTUVprev,center,PrimLabel,iPrim)
                                  ENDIF
                               ELSE
                                  !chose one of the possibilities
                                  IF(TREC)THEN
                                     !                         print*,'!C TRECURRENCE iTUV',iTUV
                                     CALL TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,LUFILE,Gen,SegQ,SegP,Segwtuv,seg1prim,nTUVprev,center,PrimLabel,iPrim)
                                  ELSEIF(UREC)THEN
                                     !                         print*,'!C URECURRENCE iTUV',iTUV
                                     call URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,LUFILE,Gen,SegQ,SegP,Segwtuv,seg1prim,nTUVprev,center,PrimLabel,iPrim)
                                  ELSEIF(VREC)THEN
                                     !                         print*,'!D VRECURRENCE iTUV',iTUV
                                     call VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,LUFILE,Gen,SegQ,SegP,Segwtuv,seg1prim,nTUVprev,center,PrimLabel,iPrim)
                                  ENDIF
                               ENDIF
                            ENDIF
                            CREATED(Tp,Up,Vp) = .TRUE.
                         ENDDO
                      ENDDO
                   ENDDO
                   deallocate(TwoTermTUVLIST)
                   deallocate(TwoTermsUsed)
                ENDDO
                IF(Collapse)THEN
                   call PrintCollapseLoopEnd(Gen,SegQ,SegP,Seg,seg1prim,LUFILE)
                ELSE
                   IF(.NOT.seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'    ENDDO'
                      WRITE(LUFILE,'(A)')'   ENDDO'
                   ENDIF
                   WRITE(LUFILE,'(A)')'  ENDDO'
                ENDIF
!                IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END PARALLEL DO'
                IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END DO'
                WRITE(LUFILE,'(A)')' end subroutine'


                deallocate(TUVINDEX)
                deallocate(TINDEX)
                deallocate(UINDEX)
                deallocate(VINDEX)
                deallocate(JINDEX)
             ENDDO
             WRITE(LUFILE,'(A)')'end module'
             close(unit = LUFILE)
          END DO
       ENDDO
    ENDDO

!BUILDRJ000
    DO GPUrun = 1,2
       CPU = .TRUE.
       IF(GPUrun.EQ.2)CPU = .FALSE.
       DoOpenMP = .FALSE.
       DoOpenACC = .FALSE.
       IF(CPU)DoOpenMP = .TRUE.
       IF(.NOT.CPU)DoOpenACC = .TRUE.
       Collapse = .TRUE.
!       IF(.NOT.CPU)Collapse = .TRUE.
       IF(.NOT.CPU)nPrimnTUV = .FALSE.
       center=1
       centerString='A'
       IF(CPU)THEN
          ARCSTRING = 'CPU'
       ELSE
          ARCSTRING = 'GPU'
       ENDIF
       DO iseg = 1,5
          IF(iseg.GT.1.AND.iseg.LT.5)CYCLE
          Gen = .FALSE.; SegQ=.FALSE.; SegP=.FALSE.;Seg=.FALSE.;Seg1Prim=.FALSE.
          IF(iseg.EQ.1)THEN
             Gen = .TRUE.      ; SegLabel = 'Gen     '; iSegLabel = 3
          ELSE!IF(iseg.EQ.5)THEN
             Seg1Prim = .TRUE. ; SegLabel = 'Seg1Prim'; iSegLabel = 8
          ENDIF
          DO I =1,48
             FileName(I:I) = ' '
          ENDDO
          WRITE(FileName,'(4A)')'BUILDRJ000'//ARCSTRING//'QP',SegLabel(1:iSegLabel),'.F90'
          print*,'FileName:',FileName
          LUFILE = 12
          open(unit = LUFILE, file=TRIM(FileName),status="unknown")
          IF(COLLAPSE)THEN
             IF(Gen)THEN
                PrimLabel = 'iP'; iPrim = 2
                nPrimLabel = 'nPrimQ*nPrimP*nPassP'; nPrim = 20 
             ELSE
                PrimLabel = 'iPassP'; iPrim = 2
                nPrimLabel = 'nPassP'; nPrim = 6 
             ENDIF
          ELSE
             IF(Gen)THEN
                PrimLabel = 'iPrimQ,iPrimP,iPassP'; iPrim = 20 
                nPrimLabel = 'nPrimQ,nPrimP,nPassP'; nPrim = 20 
             ELSE
                PrimLabel = 'iPassP'; iPrim = 6 
                nPrimLabel = 'nPassP'; nPrim = 6 
             ENDIF
          ENDIF
          
          WRITE(LUFILE,'(5A)')'MODULE AGC_',ARCSTRING,'_OBS_BUILDRJ000MOD',SegLabel(1:iSegLabel)
          MaxAngmomQP = 8
          WRITE(LUFILE,'(A)')' use IchorPrecisionModule'
          WRITE(LUFILE,'(A)')'  '
          WRITE(LUFILE,'(A)')' CONTAINS'

          !============================================================================================================
          !         BuildRj000 angmom > 2 only for center A and seg1prim and Gen 
          !============================================================================================================
          IF(center.EQ.1)THEN
             IF(Gen.OR.Seg1Prim)THEN
                DO JMAX=2,MaxAngmomQP
                   IF(center.EQ.1.AND.JMAX.GT.MaxAngmomQP)CYCLE    !DDDD
                   IF(center.EQ.3.AND.JMAX.GT.MaxAngmomQP-1)CYCLE  !PDDD
                   IF(center.EQ.2.AND.JMAX.GT.MaxAngmomQP-2)CYCLE  !PDPD
                   IF(center.EQ.4.AND.JMAX.GT.MaxAngmomQP-3)CYCLE  !PPPD
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
                   
                   WRITE(LUFILE,'(A)')''
                   IF(JMAX.LT.10)THEN
                      WRITE(LUFILE,'(A,A,I1,A)')'subroutine BuildRJ000'//ARCSTRING,SegLabel(1:iSegLabel),JMAX,'(nPassP,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
                      WRITE(LUFILE,'(A,A,I2,A)')'subroutine BuildRJ000'//ARCSTRING,SegLabel(1:iSegLabel),JMAX,'(nPassP,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                   WRITE(LUFILE,'(A)')'         & TABFJW,Pcent,Qcent,IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&'
                   WRITE(LUFILE,'(A)')'         & RJ000array)'
                   WRITE(LUFILE,'(A)')'  implicit none'
                   WRITE(LUFILE,'(A)')'  integer,intent(in) :: nPassP,nPrimP,nPrimQ'
                   WRITE(LUFILE,'(A)')'  integer,intent(in) :: MaxPasses,nAtomsA,nAtomsB'
                   WRITE(LUFILE,'(A)')'  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)'
                   WRITE(LUFILE,'(A,I2,A)')'  REAL(REALK),intent(in) :: TABFJW(0:',JMAX+3,',0:1200)'
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)'
                      WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Pcent(3,nPrimP,nAtomsA,nAtomsB),Qcent(3,nPrimQ)'
                   ELSE
                      WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: reducedExponents(1)'
                      WRITE(LUFILE,'(A)')'  real(realk),intent(in) :: Pcent(3,nAtomsA,nAtomsB),Qcent(3)'
                   ENDIF
                   IF(.NOT.Seg1Prim)THEN
                      WRITE(LUFILE,'(A,I2,A)')'  real(realk),intent(inout) :: RJ000array(0:',JMAX,',nPrimQ*nPrimP*nPassP)'
                   ELSE
                      WRITE(LUFILE,'(A,I2,A)')'  real(realk),intent(inout) :: RJ000array(0:',JMAX,',nPassP)'
                   ENDIF
                   WRITE(LUFILE,'(A)')'  !local variables'
                   IF(COLLAPSE)THEN
                      IF(Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')'  integer :: iP,iPassP,ipnt,iAtomA,iAtomB'
                      ELSE
                         WRITE(LUFILE,'(A)')'  integer :: iP,iPrimQ,iPrimP,iPassP,ipnt,iAtomA,iAtomB'
                      ENDIF
                   ELSE
                      WRITE(LUFILE,'(A)')'  integer :: iPrimQ,iPrimP,iPassP,ipnt,iAtomA,iAtomB' 
                   ENDIF
                   WRITE(LUFILE,'(A)')'  real(realk) :: mPX,mPY,mPZ,Xpq,Ypq,Zpq'
                   WRITE(LUFILE,'(A)')'  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL' 
                   WRITE(LUFILE,'(A,I2,A)')'  real(realk) :: RJ000(0:',JMAX,')'
                   WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk'
                   WRITE(LUFILE,'(A)')'  real(realk),parameter :: D2=2.0E0_realk'
                   WRITE(LUFILE,'(A,ES24.16,A)')'  REAL(REALK),PARAMETER :: D2JP36=',2.d0*JMAX + 36.d0,'_realk'
                   WRITE(LUFILE,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D03333=1.0E0_realk/3.0E0_realk'
                   WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk'
                   WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk'
                   WRITE(LUFILE,'(A)')'  REAL(REALK),PARAMETER :: SMALL = 1E-15_realk,D12 = 12.0E0_realk'
                   WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk'
                   WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk'
                   WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk'
                   WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk'
                   WRITE(LUFILE,'(A)')'  Real(realk), parameter :: PI=3.14159265358979323846E0_realk'
                   WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk'
                   WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2'
                   WRITE(LUFILE,'(A)')'  REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI'
                   IF(DoOpenMP)THEN
                      !OPENMP
                      WRITE(LUFILE,'(A)')  '!$OMP DO &'
!                      WRITE(LUFILE,'(A)')  '!$OMP PARALLEL DO DEFAULT(none) &'
                      WRITE(LUFILE,'(A)')  '!$OMP PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&'
                      IF(COLLAPSE)THEN
                       IF(Seg1Prim)THEN
                        WRITE(LUFILE,'(A)')'!$OMP         iP,iPassP,&'
                       ELSE
                        WRITE(LUFILE,'(A)')'!$OMP         iP,iPrimQ,iPrimP,iPassP,&'
                       ENDIF
                      ELSE
                       WRITE(LUFILE,'(A)') '!$OMP         iPrimQ,iPrimP,iPassP,&'
                      ENDIF
                      WRITE(LUFILE,'(A)')  '!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&'
                      WRITE(LUFILE,'(A)')  '!$OMP         mPX,mPY,mPZ,RWVAL,GVAL) '
!                      WRITE(LUFILE,'(A)')  '!$OMP         mPX,mPY,mPZ,RWVAL,GVAL) &'
!                      WRITE(LUFILE,'(A)')  '!$OMP SHARED(nPassP,nPrimP,nPrimQ,IatomApass,IatomBpass,&'
!                      WRITE(LUFILE,'(A)')  '!$OMP        TABFJW,reducedExponents,Pcent,Qcent,RJ000array)'
                   ENDIF
                   IF(DoOpenACC)THEN
                      WRITE(LUFILE,'(A)')  '!$ACC parallel loop &'
                      WRITE(LUFILE,'(A)')  '!$ACC present(nPassP,nPrimP,nPrimQ,IatomApass,IatomBpass,&'
                      WRITE(LUFILE,'(A)')  '!$ACC         TABFJW,reducedExponents,Pcent,Qcent,RJ000array)       &'
                      WRITE(LUFILE,'(A)')  '!$ACC private(iAtomA,iAtomB,Xpq,Ypq,Zpq,&'
                      IF(COLLAPSE)THEN
                       IF(Seg1Prim)THEN
                        WRITE(LUFILE,'(A)')'!$ACC         iP,iPassP,&'
                       ELSE
                        WRITE(LUFILE,'(A)')'!$ACC         iP,iPrimQ,iPrimP,iPassP,&'
                       ENDIF
                      ELSE
                       WRITE(LUFILE,'(A)') '!$ACC         iPrimQ,iPrimP,iPassP,&'
                      ENDIF
                      WRITE(LUFILE,'(A)')  '!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&'
                      WRITE(LUFILE,'(A)')  '!$ACC         mPX,mPY,mPZ,RWVAL,GVAL)'
                   ENDIF
                   IF(Collapse)THEN
                      IF(Seg1Prim)THEN
                         WRITE(LUFILE,'(A)')'  DO iP = 1,nPassP'
                         WRITE(LUFILE,'(A)')'   iPassP = iP'
                      ELSE !gen
                         WRITE(LUFILE,'(A)')'  DO iP = 1,nPrimQ*nPrimP*nPassP'
                         WRITE(LUFILE,'(A)')'   iPrimQ = mod(IP-1,nPrimQ)+1'
                         WRITE(LUFILE,'(A)')'   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1'
                         WRITE(LUFILE,'(A)')'   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1'
                      ENDIF
                      WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
                      WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
                      IF(.NOT.seg1Prim)THEN
                         WRITE(LUFILE,'(A)')'    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)'
                         WRITE(LUFILE,'(A)')'    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)'
                         WRITE(LUFILE,'(A)')'    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)'
                         WRITE(LUFILE,'(A)')'     Xpq = mPX + Qcent(1,iPrimQ)'
                         WRITE(LUFILE,'(A)')'     Ypq = mPY + Qcent(2,iPrimQ)'
                         WRITE(LUFILE,'(A)')'     Zpq = mPZ + Qcent(3,iPrimQ)'
                      ELSE
                         WRITE(LUFILE,'(A)')'     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)'
                         WRITE(LUFILE,'(A)')'     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)'
                         WRITE(LUFILE,'(A)')'     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)'
                      ENDIF
                   ELSE
                      WRITE(LUFILE,'(A)')'  DO iPassP = 1,nPassP'
                      WRITE(LUFILE,'(A)')'   iAtomA = iAtomApass(iPassP)'
                      WRITE(LUFILE,'(A)')'   iAtomB = iAtomBpass(iPassP)'
                      IF(.NOT.seg1Prim)THEN
                         WRITE(LUFILE,'(A)')'   DO iPrimP=1, nPrimP'
                         WRITE(LUFILE,'(A)')'    mPX = -Pcent(1,iPrimP,iAtomA,iAtomB)'
                         WRITE(LUFILE,'(A)')'    mPY = -Pcent(2,iPrimP,iAtomA,iAtomB)'
                         WRITE(LUFILE,'(A)')'    mPZ = -Pcent(3,iPrimP,iAtomA,iAtomB)'
                         WRITE(LUFILE,'(A)')'    DO iPrimQ=1, nPrimQ'
                         WRITE(LUFILE,'(A)')'     Xpq = mPX + Qcent(1,iPrimQ)'
                         WRITE(LUFILE,'(A)')'     Ypq = mPY + Qcent(2,iPrimQ)'
                         WRITE(LUFILE,'(A)')'     Zpq = mPZ + Qcent(3,iPrimQ)'
                      ELSE
                         WRITE(LUFILE,'(A)')'     Xpq = Qcent(1)-Pcent(1,iAtomA,iAtomB)'
                         WRITE(LUFILE,'(A)')'     Ypq = Qcent(2)-Pcent(2,iAtomA,iAtomB)'
                         WRITE(LUFILE,'(A)')'     Zpq = Qcent(3)-Pcent(3,iAtomA,iAtomB)'
                      ENDIF
                   ENDIF
                   WRITE(LUFILE,'(A)')'     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq'
                   IF(.NOT.seg1Prim)THEN
                      WRITE(LUFILE,'(A)')'     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance'
                   ELSE
                      WRITE(LUFILE,'(A)')'     WVAL = reducedExponents(1)*squaredDistance'
                   ENDIF
                   WRITE(LUFILE,'(A)')'     !  0 < WVAL < 12 '
                   WRITE(LUFILE,'(A)')'     IF (WVAL .LT. D12) THEN'
                   WRITE(LUFILE,'(A)')'      IPNT = NINT(D100*WVAL)'
                   WRITE(LUFILE,'(A)')'      WDIFF = WVAL - TENTH*IPNT'
                   WRITE(LUFILE,'(A)')'      W2    = WDIFF*WDIFF'
                   WRITE(LUFILE,'(A)')'      W3    = W2*WDIFF'
                   WRITE(LUFILE,'(A)')'      W2    = W2*D05'
                   WRITE(LUFILE,'(A)')'      W3    = W3*COEF3'
                   DO J=0,JMAX
                      IF(COLLAPSE)THEN
                         WRITE(LUFILE,'(A,I2,A,I2,A,I2,A,I2,A,I2,A)')'      RJ000Array(',J,',iP) = TABFJW(',J,',IPNT)-TABFJW(',J+1,',IPNT)*WDIFF+TABFJW(',J+2,',IPNT)*W2+TABFJW(',J+3,',IPNT)*W3'
                      ELSE
                         WRITE(LUFILE,'(A,I2,A,I2,A,I2,A,I2,A,I2,A)')'      RJ000Array(',J,',iPrimQ,iPrimP,iPassP) = TABFJW(',J,',IPNT)-TABFJW(',J+1,',IPNT)*WDIFF+TABFJW(',J+2,',IPNT)*W2+TABFJW(',J+3,',IPNT)*W3'
                      ENDIF
                   ENDDO
                   WRITE(LUFILE,'(A)')'     !  12 < WVAL <= (2J+36) '
                   WRITE(LUFILE,'(A)')'     ELSE IF (WVAL.LE.D2JP36) THEN'
                   WRITE(LUFILE,'(A)')'      REXPW = D05*EXP(-WVAL)'
                   WRITE(LUFILE,'(A)')'      RWVAL = D1/WVAL'
                   WRITE(LUFILE,'(A)')'      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))'
                   WRITE(LUFILE,'(A)')'      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
                   DO J=1,JMAX
                      WRITE(LUFILE,'(A,I2,A,I2,A,I2,A)')'      RJ000(',J,') = RWVAL*((',J,' - D05)*RJ000(',J-1,')-REXPW)'          
                   ENDDO
                   IF(COLLAPSE)THEN
                      WRITE(LUFILE,'(A)')'      RJ000Array( 0,iP) = RJ000(0)'
                   ELSE
                      WRITE(LUFILE,'(A)')'      RJ000Array( 0,iPrimQ,iPrimP,iPassP) = RJ000(0)'
                   ENDIF
                   DO J=1,JMAX
                      IF(COLLAPSE)THEN
                         WRITE(LUFILE,'(A,I2,A,I2,A)')'      RJ000Array(',J,',iP) = RJ000(',J,')'
                      ELSE
                         WRITE(LUFILE,'(A,I2,A,I2,A)')'      RJ000Array(',J,',iPrimQ,iPrimP,iPassP) = RJ000(',J,')'
                      ENDIF
                   ENDDO
                   WRITE(LUFILE,'(A)')'     !  (2J+36) < WVAL '
                   WRITE(LUFILE,'(A)')'     ELSE'
                   WRITE(LUFILE,'(A)')'      RWVAL = PID4/WVAL'
                   WRITE(LUFILE,'(A)')'      RJ000(0) = SQRT(RWVAL)'
                   !        WRITE(LUFILE,'(A)')'        RJ000(0,ipq,ipassq) = SQRT(RWVAL)'
                   WRITE(LUFILE,'(A)')'      RWVAL = RWVAL*PID4I'
                   DO J=1,JMAX
                      WRITE(LUFILE,'(A,I2,A,I2,A,I2,A)')'      RJ000(',J,') = RWVAL*(',J,' - D05)*RJ000(',J-1,')'
                   ENDDO
                   IF(COLLAPSE)THEN
                      WRITE(LUFILE,'(A)')'      RJ000Array( 0,iP) = RJ000(0)'
                   ELSE
                      WRITE(LUFILE,'(A)')'      RJ000Array( 0,iPrimQ,iPrimP,iPassP) = RJ000(0)'
                   ENDIF
                   DO J=1,JMAX
                      IF(COLLAPSE)THEN
                         WRITE(LUFILE,'(A,I2,A,I2,A)')'      RJ000Array(',J,',iP) = RJ000(',J,')'
                      ELSE
                         WRITE(LUFILE,'(A,I2,A,I2,A)')'      RJ000Array(',J,',iPrimQ,iPrimP,iPassP) = RJ000(',J,')'
                      ENDIF
                   ENDDO
                   WRITE(LUFILE,'(A)')'     ENDIF'
                   IF(.NOT.seg1Prim)THEN
                      IF(.NOT.COLLAPSE)THEN
                         WRITE(LUFILE,'(A)')'    ENDDO'
                         WRITE(LUFILE,'(A)')'   ENDDO'
                      ENDIF
                   ENDIF
                   WRITE(LUFILE,'(A)')'  ENDDO'
                   IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END DO'
!                   IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END PARALLEL DO'
                   WRITE(LUFILE,'(A)')' end subroutine'
                   deallocate(TUVINDEX)
                   deallocate(TINDEX)
                   deallocate(UINDEX)
                   deallocate(VINDEX)
                   deallocate(JINDEX)
                ENDDO
             ENDIF
          ENDIF
          WRITE(LUFILE,'(A)')'end module'
          close(unit = LUFILE)
       ENDDO
    END DO
  END subroutine PASSsub

  Subroutine PrintCollapseLoopEnd(Gen,SegQ,SegP,Seg,seg1prim,LUFILE)
    implicit none
    logical,intent(in) :: Gen,SegQ,SegP,Seg,seg1prim
    integer,intent(in) :: LUFILE   
    IF(Gen)THEN
       WRITE(LUFILE,'(A)') '  ENDDO !iP = 1,nPrimQ*nPrimP*nPassP'
    ELSEIF(SegQ)THEN
       WRITE(LUFILE,'(A)') '   ENDDO !iPrimP=1, nPrimP'       
       WRITE(LUFILE,'(A)') '  ENDDO !iP = 1,nPrimQ*nPassP'
    ELSEIF(SegP)THEN
       WRITE(LUFILE,'(A)') '   ENDDO !iPrimQ=1, nPrimQ'
       WRITE(LUFILE,'(A)') '  ENDDO !iP = 1,nPrimP*nPassP'
    ELSEIF(Seg)THEN
       WRITE(LUFILE,'(A)') '   ENDDO !iPrimQP = 1,nPrimQ*nPrimP'
       WRITE(LUFILE,'(A)') '  ENDDO !iP = 1,nPassP'
    ELSEIF(Seg1Prim)THEN
       WRITE(LUFILE,'(A)') '  ENDDO !iP = 1,nPassP'
    ENDIF
  End Subroutine PrintCollapseLoopEnd

  Subroutine PrintCollapseLoopStart(Gen,SegQ,SegP,Seg,seg1prim,LUFILE)
    implicit none
    logical,intent(in) :: Gen,SegQ,SegP,Seg,seg1prim
    integer,intent(in) :: LUFILE   
    IF(Gen)THEN
       WRITE(LUFILE,'(A)') '  DO iP = 1,nPrimQ*nPrimP*nPassP'
       WRITE(LUFILE,'(A)') '   iPrimQ = mod(IP-1,nPrimQ)+1'
       WRITE(LUFILE,'(A)') '   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1'
       WRITE(LUFILE,'(A)') '   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1'
    ELSEIF(SegP)THEN
       WRITE(LUFILE,'(A)') '  DO iP = 1,nPrimQ*nPassP'
       WRITE(LUFILE,'(A)') '   DO iPrimP=1, nPrimP'
       WRITE(LUFILE,'(A)') '    iPrimQ = iP - ((iP-1)/nPrimQ)*nPrimQ'
       WRITE(LUFILE,'(A)') '    iPassP = (iP-1)/nPrimQ + 1'
    ELSEIF(SegQ)THEN
       WRITE(LUFILE,'(A)') '  DO iP = 1,nPrimP*nPassP'
       WRITE(LUFILE,'(A)') '   DO iPrimQ=1, nPrimQ'
       WRITE(LUFILE,'(A)') '    iPrimP = iP - ((iP-1)/nPrimP)*nPrimP'
       WRITE(LUFILE,'(A)') '    iPassP = (iP-1)/nPrimP + 1'
    ELSEIF(Seg)THEN
!       WRITE(LUFILE,'(A)') '  DO iP = 1,nPrimP*nPrimQ*nPassP'
!       WRITE(LUFILE,'(A)') '   iPrimQ = mod(IP-1,nPrimQ)+1'
!       WRITE(LUFILE,'(A)') '   iPrimP = mod((IP-(mod(IP-1,nPrimQ)+1))/nPrimQ,nPrimP)+1'
!       WRITE(LUFILE,'(A)') '   iPassP = (IP-1)/(nPrimQ*nPrimP) + 1'
       WRITE(LUFILE,'(A)') '  DO iP = 1,nPassP'
       WRITE(LUFILE,'(A)') '   DO iPrimQP=1,nPrimQ*nPrimP'
       WRITE(LUFILE,'(A)') '    iPrimQ = iPrimQP - ((iPrimQP-1)/nPrimQ)*nPrimQ'
       WRITE(LUFILE,'(A)') '    iPrimP = (iPrimQP-1)/nPrimQ + 1'       
       WRITE(LUFILE,'(A)') '    iPassP = iP'
    ELSEIF(Seg1Prim)THEN
       WRITE(LUFILE,'(A)') '  DO iP = 1,nPassP'
       WRITE(LUFILE,'(A)') '   iPassP = iP'
    ENDIF
  End Subroutine PrintCollapseLoopStart

  Subroutine PrintOpenMP(Gen,SegQ,SegP,Seg,seg1prim,LUFILE,JMAX,Collapse,Center,centerstring,OpenMP,OpenACC)
    implicit none
    logical,intent(in) :: Gen,SegQ,SegP,Seg,seg1prim,Collapse
    integer,intent(in) :: LUFILE,JMAX,Center                
    logical,intent(in) :: OpenMP,OpenACC
    integer :: JTMP,iSHARED
    logical :: INCLUDESHARED
    character(len=1) :: centerString
    character(len=5) :: DIR
    character(len=7) :: SHARED
    DIR = '!    '     
    SHARED = '       '
    iSHARED = 7
    IF(OpenMP)THEN
       DIR = '!$OMP'
       SHARED = 'SHARED'
       iSHARED = 6
!       WRITE(LUFILE,'(2A)')DIR,' PARALLEL DO DEFAULT(none) &'
       WRITE(LUFILE,'(2A)')DIR,' DO &'
    ENDIF
    IF(OpenACC)THEN
       DIR = '!$ACC' 
       SHARED = 'PRESENT'
       iSHARED = 7
       WRITE(LUFILE,'(2A)')DIR,' PARALLEL LOOP &'
    ENDIF
    INCLUDESHARED = OpenACC

    WRITE(LUFILE,'(2A)')DIR,' PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&'
    IF(JMAX.EQ.0.OR.JMAX.EQ.1)THEN !RJ000 calc
       WRITE(LUFILE,'(2A)')DIR,'         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&'
       WRITE(LUFILE,'(2A)')DIR,'         RWVAL,GVAL,&'
    ENDIF
    IF(JMAX.EQ.0)THEN
       WRITE(LUFILE,'(2A)')DIR,'         Px,Py,Pz,&'             
    ELSE
       WRITE(LUFILE,'(2A)')DIR,'         mPx,mPy,mPz,&'             
    ENDIF
    IF(JMAX.GE.1)THEN
       IF(center.EQ.1)WRITE(LUFILE,'(2A)')DIR,'         Ax,Ay,Az,Xpa,Ypa,Zpa,alphaP,&'             
       IF(center.EQ.2)WRITE(LUFILE,'(2A)')DIR,'         Bx,By,Bz,Xpb,Ypb,Zpb,alphaP,&'             
       IF(center.EQ.3)WRITE(LUFILE,'(2A)')DIR,'         Xqc,Yqc,Zqc,alphaQ,&'             
       IF(center.EQ.4)WRITE(LUFILE,'(2A)')DIR,'         Xqd,Yqd,Zqd,alphaQ,&'             
       WRITE(LUFILE,'(2A)')DIR,'         alphaXpq,alphaYpq,alphaZpq,&'
    ENDIF
    IF(center.GT.2)THEN
       IF(JMAX.GT.1)THEN
          WRITE(LUFILE,'(2A)')DIR,'         invexpQ,inv2expQ,&'   
       ELSE
          WRITE(LUFILE,'(2A)')DIR,'         invexpQ,&'   
       ENDIF
    ELSE
       IF(JMAX.GT.1)THEN
          WRITE(LUFILE,'(2A)')DIR,'         invexpP,inv2expP,&'  
       ELSEIF(JMAX.GT.0)THEN
          WRITE(LUFILE,'(2A)')DIR,'         invexpP,&'   
       ENDIF
    ENDIF
    IF(JMAX.GT.0)THEN    
       WRITE(LUFILE,'(2A)')DIR,'         PREF,&'   
    ENDIF
    IF(JMAX.EQ.1)THEN    
       WRITE(LUFILE,'(2A)')DIR,'         TMP1,TMP2,&'   
    ENDIF
    IF(JMAX.GT.1)THEN
       IF(.NOT.Gen)THEN
          WRITE(LUFILE,'(2A)')DIR,'         TMPAUXarray,&'
       ENDIF
       DO JTMP=1,JMAX
          WRITE(LUFILE,'(A,A,I1,A)')DIR,'         TMParray',JTMP,',&'
       ENDDO
       WRITE(LUFILE,'(2A)')DIR,'         TwoTerms,&'
    ENDIF
    IF(INCLUDESHARED)THEN
       IF(COLLAPSE)THEN
          IF(Seg1Prim)THEN
             WRITE(LUFILE,'(2A)')DIR,'         iP,iPassP) &'
          ELSEIF(Seg)THEN
             WRITE(LUFILE,'(2A)')DIR,'         iP,iPrimQ,iPrimP,iPrimQP,iPassP) &'
          ELSE !Gen
             WRITE(LUFILE,'(2A)')DIR,'         iP,iPrimQ,iPrimP,iPassP) &'
          ENDIF
       ELSE
          IF(.Not.Seg1Prim)THEN
             WRITE(LUFILE,'(2A)')DIR,'         iPrimQ,iPrimP,iPassP) &'
          ELSE
             WRITE(LUFILE,'(2A)')DIR,'         iPassP) &'
          ENDIF
       ENDIF
    ELSE
       IF(COLLAPSE)THEN
          IF(Seg1Prim)THEN
             WRITE(LUFILE,'(2A)')DIR,'         iP,iPassP)'
          ELSEIF(Seg)THEN
             WRITE(LUFILE,'(2A)')DIR,'         iP,iPrimQ,iPrimP,iPrimQP,iPassP)'
          ELSE !Gen
             WRITE(LUFILE,'(2A)')DIR,'         iP,iPrimQ,iPrimP,iPassP)'
          ENDIF
       ELSE
          IF(.Not.Seg1Prim)THEN
             WRITE(LUFILE,'(2A)')DIR,'         iPrimQ,iPrimP,iPassP)'
          ELSE
             WRITE(LUFILE,'(2A)')DIR,'         iPassP)'
          ENDIF
       ENDIF
    ENDIF
    IF(INCLUDESHARED)THEN
       IF(JMAX.LT.2)THEN
          WRITE(LUFILE,'(4A)')DIR,' ',SHARED(1:iSHARED),'(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,TABFJW,&'
       ELSE
          WRITE(LUFILE,'(4A)')DIR,' ',SHARED(1:iSHARED),'(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,RJ000Array,&'
       ENDIF
       WRITE(LUFILE,'(2A)')DIR,'        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray,&'
       IF(center.GT.2)THEN
          IF(JMAX.GT.0)THEN
             WRITE(LUFILE,'(4A)')DIR,'        Qexp,',centerstring,'center, &'
          ENDIF
       ELSE
          IF(JMAX.GT.0)THEN
             WRITE(LUFILE,'(4A)')DIR,'        Pexp,',centerstring,'center, &'
          ENDIF
       ENDIF
       WRITE(LUFILE,'(2A)')      DIR,'        nPrimP,nPrimQ,nPassP)'
    ENDIF
!!This is a CPU code so if OpenMP it uses OpenMP. However, I would like 
!!to test OpenACC so in case of no OpenMP and OpenACC it uses OpenACC)
!#ifdef VAR_OMP
!!$OMP DO &
!!$OMP PRIVATE(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!!$OMP         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!!$OMP         RWVAL,GVAL,&
!!$OMP         Px,Py,Pz,&
!!$OMP         iP,iPrimQ,iPrimP,iPassP) 
!!!$OMP SHARED(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,TABFJW,&
!!!$OMP        integralPrefactor,PpreExpFac,QpreExpFac,AUXarray)
!#else
!#ifdef VAR_OPENACC
!!$ACC parallel loop gang worker vector &
!!$ACC present(iAtomApass,iAtomBpass,Pcent,Qcent,reducedExponents,TABFJW,&
!!$ACC         integralPrefactor,PpreExpFac,QpreExpFac,AUXarray)&
!!$ACC private(iAtomA,iAtomB,Xpq,Ypq,Zpq,&
!!$ACC         squaredDistance,WVAL,IPNT,WDIFF,W2,W3,RJ000,REXPW,&
!!$ACC         RWVAL,GVAL,&
!!$ACC         Px,Py,Pz,&
!!$ACC         iP,iPrimQ,iPrimP,iPassP)  
!#endif
!#endif
  end Subroutine PrintOpenMP
  
  Subroutine PrintCollapseInitLoop(Gen,SegQ,SegP,Seg,seg1prim,LUFILE,JMAX,nTUV,DoOpenMP)
    implicit none
    logical,intent(in) :: Gen,SegQ,SegP,Seg,seg1prim,DoOpenMP
    integer,intent(in) :: LUFILE,JMAX,nTUV
    IF(SegQ.OR.SegP)THEN
       IF(JMAX.LT.2)THEN
          IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP DO PRIVATE(iP)'
       ELSE
          IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP DO PRIVATE(iP,iTUV)'          
       ENDIF
       IF(SegP)WRITE(LUFILE,'(A)')'  DO iP = 1,nPrimQ*nPassP'
       IF(SegQ)WRITE(LUFILE,'(A)')'  DO iP = 1,nPrimP*nPassP'
       IF(JMAX.EQ.0)THEN
          WRITE(LUFILE,'(A)')'    AUXarray(iP)=0.0E0_realk'
       ELSEIF(JMAX.EQ.1)THEN
          WRITE(LUFILE,'(A)')'    AUXarray(1,iP)=0.0E0_realk'
          WRITE(LUFILE,'(A)')'    AUXarray(2,iP)=0.0E0_realk'
          WRITE(LUFILE,'(A)')'    AUXarray(3,iP)=0.0E0_realk'
          WRITE(LUFILE,'(A)')'    AUXarray(4,iP)=0.0E0_realk'
       ELSE
          WRITE(LUFILE,'(A,I5)')'    DO iTUV=1,',nTUV
          WRITE(LUFILE,'(A)')'     AUXarray(iTUV,iP)=0.0E0_realk'
          WRITE(LUFILE,'(A)')'    ENDDO'
       ENDIF
       WRITE(LUFILE,'(A)')'  ENDDO'
       IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END DO'
    ELSEIF(Seg)THEN
       IF(JMAX.LT.2)THEN
          IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP DO PRIVATE(iPassP)'
       ELSE
          IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP DO PRIVATE(iPassP,iTUV)'
       ENDIF
       WRITE(LUFILE,'(A)')'  DO iPassP = 1,nPassP'
       IF(JMAX.EQ.0)THEN
          WRITE(LUFILE,'(A)')'   AUXarray(iPassP)=0.0E0_realk'
       ELSEIF(JMAX.EQ.1)THEN
          WRITE(LUFILE,'(A)')'   AUXarray(1,iPassP)=0.0E0_realk'
          WRITE(LUFILE,'(A)')'   AUXarray(2,iPassP)=0.0E0_realk'
          WRITE(LUFILE,'(A)')'   AUXarray(3,iPassP)=0.0E0_realk'
          WRITE(LUFILE,'(A)')'   AUXarray(4,iPassP)=0.0E0_realk'
       ELSE
          WRITE(LUFILE,'(A,I5)')'   DO iTUV=1,',nTUV
          WRITE(LUFILE,'(A)')'    AUXarray(iTUV,iPassP)=0.0E0_realk'
          WRITE(LUFILE,'(A)')'   ENDDO'
       ENDIF
       WRITE(LUFILE,'(A)')'  ENDDO'
       IF(DoOpenMP)WRITE(LUFILE,'(A)')'!$OMP END DO'
    ENDIF
  End Subroutine PrintCollapseInitLoop

subroutine TwoTerms1(J1,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,JTMP,I,J1
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: nTUVLIST,nTUVLISTactual
integer :: TwoTermTUVLIST(nTUVLIST)
logical :: TwoTermsUsed(nTUVLIST)
logical :: Unique,TREC,UREC,VREC
!TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
!               +T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
!TwoTerms are all  TUV(T-2) that is non zero ! to start with maybe not all needed
TwoTermTUVLIST = 0 
TwoTermsUsed = .FALSE.
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

end subroutine TwoTerms1

subroutine WriteTwoTerms(J1,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,&
     & TwoTermTUVLIST,JTMP,TwoTermsUsed,LUPRI,Gen,SegQ,SegP,Seg,seg1prim,center,&
     & COLLAPSE,PrimLabel,iPrim)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,JTMP,I,J1
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),center
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),Gen,SegQ,SegP,Seg,seg1prim
integer :: nTUVLIST,nTUVLISTactual,LUPRI
integer :: TwoTermTUVLIST(nTUVLIST),iPrim
logical :: Unique,TREC,UREC,VREC,TwoTermsUsed(nTUVLIST),COLLAPSE
Character(len=20) :: PrimLabel

i=0 
DO ituvP2 = 1,nTUVLISTactual
   IF(TwoTermsUsed(ituvP2))THEN
      i=i+1 
      call initString(5)
      IF(J1.EQ.1)THEN
         !         WRITE(*,'(A,I3,A,I3,A,I1,A,I3,A)')&
         !              &'     TwoTerms(',ituvP2,') = inv2expP*(AuxArray(',TwoTermTUVLIST(ituvP2),',IP) + alphaP*TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',2))'
         call AddToString('TwoTerms(')
!         call AddToString(ituvP2)
         call AddToString(i)
         IF(center.LT.3)THEN
            call AddToString(') = inv2expP*(')
         ELSEIF(center.GT.2)THEN
            call AddToString(') = inv2expQ*(')
         ENDIF
         IF(Gen)THEN
            call AddToString('AuxArray(')
         ELSE
            call AddToString('TMPAuxArray(')
         ENDIF
         call AddToString(TwoTermTUVLIST(ituvP2))
         IF(Gen)THEN
            IF(COLLAPSE)THEN
               call AddToString(',IP)')
            ELSE
               call AddToString(',')
               call AddToString(PrimLabel(1:iPrim))
               call AddToString(')')
            ENDIF
         ELSE
            call AddToString(')')
         ENDIF
         call AddToString(' + ')
         IF(center.LT.3)THEN
            call AddToString('alphaP*TmpArray')
         ELSEIF(center.GT.2)THEN
            call AddToString('alphaQ*TmpArray')
         ENDIF
         call AddToString(JTMP-1)
         call AddToString('(')
         call AddToString(TwoTermTUVLIST(ituvP2))
         call AddToString(',2))')
      ELSE
         call AddToString('TwoTerms(')
!         call AddToString(ituvP2)
         call AddToString(i)
         IF(center.LT.3)THEN
            call AddToString(') = inv2expP*(')
         ELSEIF(center.GT.2)THEN
            call AddToString(') = inv2expQ*(')
         ENDIF
         call AddToString('TmpArray')
         call AddToString(JTMP-1)
         call AddToString('(')
         call AddToString(TwoTermTUVLIST(ituvP2))
         call AddToString(',')
         call AddToString(J1)
         IF(center.LT.3)THEN
            call AddToString(') + alphaP*TmpArray')
         ELSEIF(center.GT.2)THEN
            call AddToString(') + alphaQ*TmpArray')
         ENDIF
         call AddToString(JTMP-1)
         call AddToString('(')
         call AddToString(TwoTermTUVLIST(ituvP2))
         call AddToString(',')
         call AddToString(J1+1)
         call AddToString('))')
!         WRITE(*,'(A,I3,A,I1,A,I3,A,I3,A,I1,A,I3,A,I3,A)')&
         !              &'     TwoTerms(',ituvP2,') = inv2expP*(TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',',J1,') + alphaP*TmpArray',JTMP-1,'(',TwoTermTUVLIST(ituvP2),',',J1+1,'))'
      ENDIF
      call writeString(LUPRI)
   ENDIF
ENDDO

i=0 
DO ituvP2 = 1,nTUVLISTactual
   IF(TwoTermsUsed(ituvP2))THEN
      i=i+1 
      TwoTermTUVLIST(i) = TwoTermTUVLIST(ituvP2)
   ENDIF
ENDDO
nTUVLISTactual = i
DO ituvP2 = i+1,nTUVLIST
   TwoTermTUVLIST(ituvP2) = 0
ENDDO
end subroutine WriteTwoTerms

subroutine TRECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,&
     & nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: nTUVLIST,nTUVLISTactual,JTMP,TM2
integer :: TwoTermTUVLIST(nTUVLIST)
logical :: TwoTermsUsed(nTUVLIST)
character(len=3) :: DIRECTIONSTRING
!TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
!               +T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
!TwoTerms are all  TUV(T-2) that is non zero ! to start with maybe not all needed
ituvP0 = TUVINDEX(Tp,Up,Vp)
ituvP1 = TUVINDEX(Tp-1,Up,Vp)

ituvP2 = 1
TM1 = Tp-1
TM2 = Tp-2

IF(TM2.GE.0)THEN
   ituvP2 = TUVINDEX(Tp-2,Up,Vp)
   DO I=1,nTUVLISTactual
      IF(ituvP2.EQ.TwoTermTUVLIST(I))THEN
         iTwoTerms = I
         TwoTermsUsed(iTwoTerms) = .TRUE.
      ENDIF
   ENDDO
ENDIF

end subroutine TRECURRENCETWOTERMS

subroutine URECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,&
     & nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: nTUVLIST,nTUVLISTactual,JTMP,TM2
logical :: TwoTermsUsed(nTUVLIST)
integer :: TwoTermTUVLIST(nTUVLIST)
character(len=3) :: DIRECTIONSTRING
!TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
!               +T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
!TwoTerms are all  TUV(T-2) that is non zero ! to start with maybe not all needed
ituvP0 = TUVINDEX(Tp,Up,Vp)
ituvP1 = TUVINDEX(Tp,Up-1,Vp)

ituvP2 = 1
TM1 = Up-1
TM2 = Up-2

IF(TM2.GE.0)THEN
   ituvP2 = TUVINDEX(Tp,Up-2,Vp)
   DO I=1,nTUVLISTactual
      IF(ituvP2.EQ.TwoTermTUVLIST(I))THEN
         iTwoTerms = I
         TwoTermsUsed(iTwoTerms) = .TRUE.
      ENDIF
   ENDDO
ENDIF

end subroutine URECURRENCETWOTERMS

subroutine VRECURRENCETWOTERMS(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,&
     & nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
integer :: nTUVLIST,nTUVLISTactual,JTMP,TM2
logical :: TwoTermsUsed(nTUVLIST)
integer :: TwoTermTUVLIST(nTUVLIST)
character(len=3) :: DIRECTIONSTRING
!TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
!               +T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))
!TwoTerms are all  TUV(T-2) that is non zero ! to start with maybe not all needed
ituvP0 = TUVINDEX(Tp,Up,Vp)
ituvP1 = TUVINDEX(Tp,Up,Vp-1)

ituvP2 = 1
TM1 = Vp-1
TM2 = Vp-2

IF(TM2.GE.0)THEN
   ituvP2 = TUVINDEX(Tp,Up,Vp-2)
   DO I=1,nTUVLISTactual
      IF(ituvP2.EQ.TwoTermTUVLIST(I))THEN
         iTwoTerms = I
         TwoTermsUsed(iTwoTerms) = .TRUE.
      ENDIF
   ENDDO
ENDIF

end subroutine VRECURRENCETWOTERMS

subroutine TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev,center,PrimLabel,iPrim)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms,lupri,nTUVprev
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),center
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),Gen,SegQ,SegP,Seg,seg1prim
integer :: nTUVLIST,nTUVLISTactual,JTMP,TM2
integer :: TwoTermTUVLIST(nTUVLIST)
character(len=3) :: DIRECTIONSTRING
integer :: iPrim
Character(len=20) :: PrimLabel

!TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)
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
IF(center.EQ.1)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Xpa','alphaXpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev,PrimLabel,iPrim)
ELSEIF(center.EQ.2)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Xpb','alphaXpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev,PrimLabel,iPrim)
ELSEIF(center.EQ.3)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Xqc','alphaXpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev,PrimLabel,iPrim)
ELSEIF(center.EQ.4)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Xqd','alphaXpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev,PrimLabel,iPrim)
ENDIF
end subroutine TRECURRENCE


subroutine URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev,center,PrimLabel,iPrim)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms,lupri,nTUVprev
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),center
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),Gen,SegQ,SegP,Seg,seg1prim
integer :: nTUVLIST,nTUVLISTactual,JTMP,TM2
integer :: TwoTermTUVLIST(nTUVLIST)
integer :: iPrim
Character(len=20) :: PrimLabel

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
IF(center.EQ.1)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Ypa','alphaYpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev,PrimLabel,iPrim)
ELSEIF(center.EQ.2)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Ypb','alphaYpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev,PrimLabel,iPrim)
ELSEIF(center.EQ.3)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Yqc','alphaYpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev,PrimLabel,iPrim)
ELSEIF(center.EQ.4)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Yqd','alphaYpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev,PrimLabel,iPrim)
ENDIF
end subroutine URECURRENCE

subroutine VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,lupri,Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev,center,PrimLabel,iPrim)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms,lupri,nTUVprev,center
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),Gen,SegQ,SegP,Seg,Seg1Prim
integer :: nTUVLIST,nTUVLISTactual,JTMP,TM2
integer :: TwoTermTUVLIST(nTUVLIST)
integer :: iPrim
Character(len=20) :: PrimLabel

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
IF(center.EQ.1)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Zpa','alphaZpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev,PrimLabel,iPrim)
ELSEIF(center.EQ.2)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Zpb','alphaZpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev,PrimLabel,iPrim)
ELSEIF(center.EQ.3)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Zqc','alphaZpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev,PrimLabel,iPrim)
ELSEIF(center.EQ.4)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Zqd','alphaZpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev,PrimLabel,iPrim)
ENDIF

end subroutine VRECURRENCE

SUBROUTINE XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
     & DIRECTIONSTRING1,DIRECTIONSTRING2,TM1,TM2,lupri,&
     & Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev,PrimLabel,iPrim)
implicit none
integer :: J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,TM1,TM2,lupri,nTUVprev
logical :: Gen,SegQ,SegP,Seg,Seg1Prim
character(len=3) :: DIRECTIONSTRING1 !Xpa
character(len=8) :: DIRECTIONSTRING2 !alphaXpq
integer :: iPrim
Character(len=20) :: PrimLabel
   
call initString(5)
IF(J.EQ.1)THEN
   IF(Gen)THEN
      !place in AuxArray and use AuxArray in X*Aux
      call AddToString('AuxArray(')
      call AddToString(ituvP0)
!      call AddToString(',IP) = ')
call AddToString(',')
call AddToString(PrimLabel(1:iPrim))
call AddToString(') = ')
      !      call AddToString('Xpa')
      call AddToString(DIRECTIONSTRING1) !Xpa
      call AddToString('*AuxArray(')
      call AddToString(ituvP1)
!      call AddToString(',IP)')
call AddToString(',')
call AddToString(PrimLabel(1:iPrim))
call AddToString(')')
   ELSE !segmentet 
      IF(ituvP0.LE.nTUVprev)THEN
         !place in TMPAuxArray and use TMPAuxArray in X*Aux
         call AddToString('TMPAuxArray(')
         call AddToString(ituvP0)
         call AddToString(') = ')
         !      call AddToString('Xpa')
         call AddToString(DIRECTIONSTRING1) !Xpa
         call AddToString('*TMPAuxArray(')
         call AddToString(ituvP1)
         call AddToString(')')
      ELSE         
         !add loop
         IF(ituvP0.EQ.nTUVprev+1)THEN
            WRITE(lupri,'(A,I5)')'     do iTUV = 1,',nTUVprev
            IF(Seg1Prim)THEN
               WRITE(lupri,'(A,A,A)')'      AuxArray(iTUV,',PrimLabel(1:iPrim),') = TMPAuxarray(iTUV)'
            ELSEIF(Seg)THEN
               WRITE(lupri,'(A,A,A,A,A)')'      AuxArray(iTUV,',PrimLabel(1:iPrim),') = AuxArray(iTUV,',PrimLabel(1:iPrim),') + TMPAuxarray(iTUV)'
            ELSE !segQ or SegP
               WRITE(lupri,'(A,A,A,A,A)')'      AuxArray(iTUV,',PrimLabel(1:iPrim),') = AuxArray(iTUV,',PrimLabel(1:iPrim),') + TMPAuxarray(iTUV)'
            ENDIF
            WRITE(lupri,'(A)')'     enddo'
         ENDIF
         !place in AuxArray and use AuxArray in X*Aux
         call AddToString('AuxArray(')
         call AddToString(ituvP0)
!         IF(Seg)THEN
!            call AddToString(',IPassP)')
!         ELSE !segQ or SegP or Gen
!            call AddToString(',IP)')
!         ENDIF
call AddToString(',')
call AddToString(PrimLabel(1:iPrim))
call AddToString(')')

         call AddToString(' = ')
         IF(.NOT.Seg1Prim)THEN
            call AddToString('AuxArray(')
            call AddToString(ituvP0)
!            IF(Seg)THEN
!               call AddToString(',IPassP)')
!            ELSE !segQ or SegP or Gen
!               call AddToString(',IP)')
!            ENDIF
call AddToString(',')
call AddToString(PrimLabel(1:iPrim))
call AddToString(')')
            call AddToString(' + ')
         ENDIF
         !      call AddToString('Xpa')
         call AddToString(DIRECTIONSTRING1) !Xpa
         call AddToString('*TMPAuxArray(')
         call AddToString(ituvP1)
         call AddToString(')')
         
!         it should be
!         ....
!         TMPAuxArray(84) = Zqd*TMPAuxArray(56) + alphaZpq*TmpArray6(56,2) + 5*TwoTerms(9)
!         TwoTerms(9) = inv2expQ*(TmpArray5(35,2) + alphaQ*TmpArray5(35,3))
!         tmpArray7(84,2) = Zqd*tmpArray6(56,2) + alphaZpq*TmpArray6(56,3) + 5*TwoTerms(9)
!         TwoTerms(9) = inv2expQ*(TmpArray5(35,3) + alphaQ*TmpArray5(35,4))
!         tmpArray7(84,3) = Zqd*tmpArray6(56,3) + alphaZpq*TmpArray6(56,4) + 5*TwoTerms(9)
!         TwoTerms(13) = inv2expQ*(TMPAuxArray(56) + alphaQ*TmpArray6(56,2))
!         TMPAuxArray(120) = Zqd*TMPAuxArray(84) + alphaZpq*TmpArray7(84,2) + 6*TwoTerms(13)
!         TwoTerms(13) = inv2expQ*(TmpArray6(56,2) + alphaQ*TmpArray6(56,3))
!         tmpArray8(120,2) = Zqd*tmpArray7(84,2) + alphaZpq*TmpArray7(84,3) + 6*TwoTerms(13)
!         do iTUV = 1,120
!            AuxArray(iTUV,Ipass) = AuxArray(iTUV,Ipass) + TMPAuxarray(iTUV)
!         enddo
!         AuxArray(121,Ipass) = AuxArray(121,Ipass) + Xqd*TMPAuxArray(85) + alphaXpq*TmpArray8(85,2) + 7*TwoTerms(1)
!         ....
!         AuxArray(165,Ipass) =  AuxArray(165,Ipass) + Zqd*TMPAuxArray(120) + alphaZpq*TmpArray8(120,2) + 7*TwoTerms(18)
!         
!         so have TMPAuxArray(1:120)
      ENDIF
   ENDIF
   call AddToString(' + ')
   !      call AddToString('alphaXpq')
   call AddToString(DIRECTIONSTRING2) !alphaXpq
   call AddToString('*TmpArray')
   call AddToString(JTMP)
   call AddToString('(')
   call AddToString(ituvP1)
   call AddToString(',2)')
   IF(TM2.GE.0)THEN
      !four term relation
      call AddToString(' + ')
      IF(TM1.EQ.1)THEN
         call AddToString('TwoTerms(')
         call AddToString(iTwoTerms)
         call AddToString(')')
      ELSE
         call AddToString(TM1)
         call AddToString('*')
         call AddToString('TwoTerms(')
         call AddToString(iTwoTerms)
         call AddToString(')')
      ENDIF
   ELSE
      !two term relation alrady done
   ENDIF
ELSE
   !place in tmpArray and use tmpArray in X*tmp
   call AddToString('tmpArray')
   call AddToString(JTMP+1)
   call AddToString('(')
   call AddToString(ituvP0)
   call AddToString(',')
   call AddToString(J)
   call AddToString(') = ')
   !      call AddToString('Xpa')
   call AddToString(DIRECTIONSTRING1) !Xpa
   call AddToString('*tmpArray')
   call AddToString(JTMP)
   call AddToString('(')
   call AddToString(ituvP1)
   call AddToString(',')
   call AddToString(J)
   call AddToString(') + ')
   !      call AddToString('alphaXpq')
   call AddToString(DIRECTIONSTRING2) !alphaXpq
   call AddToString('*TmpArray')
   call AddToString(JTMP)
   call AddToString('(')
   call AddToString(ituvP1)
   call AddToString(',')
   call AddToString(J+1)
   call AddToString(')')
   IF(TM2.GE.0)THEN
      !four term relation
      call AddToString(' + ')
      IF(TM1.EQ.1)THEN
         call AddToString('TwoTerms(')
         call AddToString(iTwoTerms)
         call AddToString(')')
      ELSE
         call AddToString(TM1)
         call AddToString('*')
         call AddToString('TwoTerms(')
         call AddToString(iTwoTerms)
         call AddToString(')')
      ENDIF
   ELSE
      !two term relation alrady done
   ENDIF
ENDIF
call writeString(lupri)
END SUBROUTINE XYZVERTICALRECURRENCE


end MODULE TESTMODULE

PROGRAM PASSprogram
use TESTMODULE
INTEGER :: JMAX,JMIN

call PASSSUB

end PROGRAM PASSprogram
