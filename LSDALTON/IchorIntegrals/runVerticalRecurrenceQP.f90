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
    integer :: LUMODA1,LUMODA2,LUMODA3,LUMODA4,LUMODA5
    integer :: LUMODB1,LUMODB2,LUMODB3,LUMODB4,LUMODB5
    integer :: LUMODC1,LUMODC2,LUMODC3,LUMODC4,LUMODC5
    integer :: LUMODD1,LUMODD2,LUMODD3,LUMODD4,LUMODD5
    integer :: non1Prim(16),pure1Prim(4),ia,ib,ic,id,GPUrun
    logical :: Gen,SegQ,SegP,Seg,Seg1Prim,CPU
    character(len=3) :: ARCSTRING
    DO GPUrun = 1,1!2
       CPU = .TRUE.
       IF(GPUrun.EQ.2)CPU = .FALSE.

       LUMODA1=1; LUMODA2=2; LUMODA3=3; LUMODA4=4; LUMODA5=5
       LUMODB1=6; LUMODB2=7; LUMODB3=8; LUMODB4=9; LUMODB5=10
       LUMODC1=11; LUMODC2=12; LUMODC3=13; LUMODC4=14; LUMODC5=15
       LUMODD1=16; LUMODD2=17; LUMODD3=18; LUMODD4=19; LUMODD5=20

       non1Prim(1) = 1; non1Prim(2) = 2; non1Prim(3) = 3; non1Prim(4) = 4
       non1Prim(5) = 6; non1Prim(6) = 7; non1Prim(7) = 8; non1Prim(8) = 9
       non1Prim(9) = 11; non1Prim(10) = 12; non1Prim(11) = 13; non1Prim(12) = 14
       non1Prim(13) = 16; non1Prim(14) = 17; non1Prim(15) = 18; non1Prim(16) = 19
       pure1Prim(1) = 5; pure1Prim(2) = 10; pure1Prim(3) = 15; pure1Prim(4) = 20
       IF(CPU)THEN
          ARCSTRING = 'CPU'
       ELSE
          ARCSTRING = 'GPU'
       ENDIF

       open(unit = LUMODA1, file="runVerticalRecurrence"//ARCSTRING//"QPA.F90",status="unknown")
       open(unit = LUMODA2, file="runVerticalRecurrence"//ARCSTRING//"QPASegQ.F90",status="unknown")
       open(unit = LUMODA3, file="runVerticalRecurrence"//ARCSTRING//"QPASegP.F90",status="unknown")
       open(unit = LUMODA4, file="runVerticalRecurrence"//ARCSTRING//"QPASeg.F90",status="unknown")
       open(unit = LUMODA5, file="runVerticalRecurrence"//ARCSTRING//"QPASeg1Prim.F90",status="unknown")

       open(unit = LUMODB1, file="runVerticalRecurrence"//ARCSTRING//"QPB.F90",status="unknown")
       open(unit = LUMODB2, file="runVerticalRecurrence"//ARCSTRING//"QPBSegQ.F90",status="unknown")
       open(unit = LUMODB3, file="runVerticalRecurrence"//ARCSTRING//"QPBSegP.F90",status="unknown")
       open(unit = LUMODB4, file="runVerticalRecurrence"//ARCSTRING//"QPBSeg.F90",status="unknown")
       open(unit = LUMODB5, file="runVerticalRecurrence"//ARCSTRING//"QPBSeg1Prim.F90",status="unknown")

       open(unit = LUMODC1, file="runVerticalRecurrence"//ARCSTRING//"QPC.F90",status="unknown")
       open(unit = LUMODC2, file="runVerticalRecurrence"//ARCSTRING//"QPCSegQ.F90",status="unknown")
       open(unit = LUMODC3, file="runVerticalRecurrence"//ARCSTRING//"QPCSegP.F90",status="unknown")
       open(unit = LUMODC4, file="runVerticalRecurrence"//ARCSTRING//"QPCSeg.F90",status="unknown")
       open(unit = LUMODC5, file="runVerticalRecurrence"//ARCSTRING//"QPCSeg1Prim.F90",status="unknown")

       open(unit = LUMODD1, file="runVerticalRecurrence"//ARCSTRING//"QPD.F90",status="unknown")
       open(unit = LUMODD2, file="runVerticalRecurrence"//ARCSTRING//"QPDSegQ.F90",status="unknown")
       open(unit = LUMODD3, file="runVerticalRecurrence"//ARCSTRING//"QPDSegP.F90",status="unknown")
       open(unit = LUMODD4, file="runVerticalRecurrence"//ARCSTRING//"QPDSeg.F90",status="unknown")
       open(unit = LUMODD5, file="runVerticalRecurrence"//ARCSTRING//"QPDSeg1Prim.F90",status="unknown")

       WRITE(LUMODA1,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODA'
       WRITE(LUMODA2,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASegQ'
       WRITE(LUMODA3,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASegP'
       WRITE(LUMODA4,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASeg'
       WRITE(LUMODA5,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASeg1Prim'

       WRITE(LUMODB1,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODB'
       WRITE(LUMODB2,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSegQ'
       WRITE(LUMODB3,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSegP'
       WRITE(LUMODB4,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSeg'
       WRITE(LUMODB5,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSeg1Prim'

       WRITE(LUMODC1,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODC'
       WRITE(LUMODC2,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSegQ'
       WRITE(LUMODC3,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSegP'
       WRITE(LUMODC4,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSeg'
       WRITE(LUMODC5,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSeg1Prim'

       WRITE(LUMODD1,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODD'
       WRITE(LUMODD2,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSegQ'
       WRITE(LUMODD3,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSegP'
       WRITE(LUMODD4,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSeg'
       WRITE(LUMODD5,'(A)')'MODULE AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSeg1Prim'

       MaxAngmomQP = 8
       do I=1,20
          WRITE(I,'(A)')' use IchorPrecisionModule'
          WRITE(I,'(A)')'  '
          WRITE(I,'(A)')' CONTAINS'
       enddo

       !========================================================================================================
       !    VerticalRecurrence 0
       !========================================================================================================

       do IA=1,5
          WRITE(IA,'(A)')''
       enddo
       WRITE(LUMODA1,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'0(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODA2,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegQ0(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODA3,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegP0(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODA4,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg0(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODA5,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1Prim0(nPasses,nPrimP,nPrimQ,&'
       do IA=1,5
          WRITE(IA,'(A)')'         & reducedExponents,TABFJW,&'
          WRITE(IA,'(A)')'         & Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&'
          WRITE(IA,'(A)')'         & AUXarray)'
          WRITE(IA,'(A)')'  implicit none'
          WRITE(IA,'(A)')'  integer,intent(in) :: nPasses,nPrimP,nPrimQ'
       enddo
       do IA=1,4
          WRITE(IA,'(A)')'  REAL(REALK),intent(in) :: reducedExponents(nPrimQ,nPrimP)'
          WRITE(IA,'(A)')'  REAL(REALK),intent(in) :: integralPrefactor(nprimQ,nPrimP)'
          WRITE(IA,'(A)')'  REAL(REALK),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)'
          WRITE(IA,'(A)')'  REAL(REALK),intent(in) :: TABFJW(0:3,0:1200)'
          WRITE(IA,'(A)')'  REAL(REALK),intent(in) :: QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)'
       enddo
       WRITE(LUMODA1,'(A)')'  real(realk),intent(inout) :: AUXarray(nPrimQ,nPrimP,nPasses)'
       WRITE(LUMODA2,'(A)')'  real(realk),intent(inout) :: AUXarray(nPrimP,nPasses)'
       WRITE(LUMODA3,'(A)')'  real(realk),intent(inout) :: AUXarray(nPrimQ,nPasses)'
       WRITE(LUMODA4,'(A)')'  real(realk),intent(inout) :: AUXarray(nPasses)'

       WRITE(LUMODA5,'(A)')'  REAL(REALK),intent(in) :: reducedExponents(1)'
       WRITE(LUMODA5,'(A)')'  REAL(REALK),intent(in) :: integralPrefactor(1)'
       WRITE(LUMODA5,'(A)')'  REAL(REALK),intent(in) :: Pcent(3),Qcent(3,nPasses)'
       WRITE(LUMODA5,'(A)')'  REAL(REALK),intent(in) :: TABFJW(0:3,0:1200)'
       WRITE(LUMODA5,'(A)')'  REAL(REALK),intent(in) :: QpreExpFac(nPasses),PpreExpFac(1)'
       WRITE(LUMODA5,'(A)')'  real(realk),intent(inout) :: AUXarray(nPasses)'

       do IA=1,5
          WRITE(IA,'(A)')'  !local variables'
          WRITE(IA,'(A,ES24.16,A)')'  REAL(REALK),PARAMETER :: D2JP36=',36.0d0,'_realk'
          WRITE(IA,'(A)')'  real(realk),parameter :: D2=2.0E0_realk'
          WRITE(IA,'(A)')'  REAL(REALK),PARAMETER :: D05 =0.5E0_realk,D1=1E0_realk'
          WRITE(IA,'(A)')'  REAL(REALK),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk'
          WRITE(IA,'(A)')'  Real(realk),parameter :: D12 = 12E0_realk, TENTH = 0.01E0_realk'
          WRITE(IA,'(A)')'  REAL(REALK),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk'
          WRITE(IA,'(A)')'  REAL(REALK),PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk'
          WRITE(IA,'(A)')'  REAL(REALK),PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk'
          WRITE(IA,'(A)')'  REAL(REALK),PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk'
          WRITE(IA,'(A)')'  REAL(REALK),PARAMETER :: GFAC2 =  D2*0.321180909E0_realk'
          WRITE(IA,'(A)')'  REAL(REALK),PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk'
          WRITE(IA,'(A)')'  Real(realk),parameter :: PI=3.14159265358979323846E0_realk'
          WRITE(IA,'(A)')'  REAL(REALK),PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk'
          WRITE(IA,'(A)')'  REAL(REALK),PARAMETER :: SQRPIH = SQRTPI/D2'
          WRITE(IA,'(A)')'  REAL(REALK),PARAMETER :: PID4 = PI/D4, PID4I = D4/PI'
          WRITE(IA,'(A)')'!  REAL(REALK),PARAMETER :: SMALL = 1E-15_realk'
          WRITE(IA,'(A)')'  Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA,WVAL,Pexpfac'
          WRITE(IA,'(A)')'  Real(realk) :: W2,W3,PX,PY,PZ,PQX,PQY,PQZ,squaredDistance,RJ000'
          WRITE(IA,'(A)')'  Integer :: IPNT,iPassQ,iPrimP,iPrimQ,iPQ'
       enddo

       do IA=1,4
          WRITE(IA,'(A)')'  DO iPassQ = 1,nPasses'
       enddo
       !seg
       WRITE(LUMODA4,'(A)')'   AUXarray(iPassQ)=0.0E0_realk'
       !segP
       WRITE(LUMODA3,'(A)')'   DO iPrimQ=1, nPrimQ'
       WRITE(LUMODA3,'(A)')'    AUXarray(iPrimQ,iPassQ)=0.0E0_realk'
       WRITE(LUMODA3,'(A)')'   ENDDO'
       do IA=1,4
          WRITE(IA,'(A)')'   DO iPrimP=1, nPrimP'
       enddo
       !segQ
       WRITE(LUMODA2,'(A)')'    AUXarray(iPrimP,iPassQ)=0.0E0_realk'
       do IA=1,4
          WRITE(IA,'(A)')'    Pexpfac = PpreExpFac(iPrimP)'
          WRITE(IA,'(A)')'    px = Pcent(1,iPrimP)'
          WRITE(IA,'(A)')'    py = Pcent(2,iPrimP)'
          WRITE(IA,'(A)')'    pz = Pcent(3,iPrimP)'
          WRITE(IA,'(A)')'    DO iPrimQ=1, nPrimQ'
          WRITE(IA,'(A)')'     pqx = px - Qcent(1,iPrimQ,iPassQ)'
          WRITE(IA,'(A)')'     pqy = py - Qcent(2,iPrimQ,iPassQ)'
          WRITE(IA,'(A)')'     pqz = pz - Qcent(3,iPrimQ,iPassQ)'
          WRITE(IA,'(A)')'     squaredDistance = pqx*pqx+pqy*pqy+pqz*pqz'
          WRITE(IA,'(A)')'     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance'
       enddo

       WRITE(LUMODA5,'(A)')'  DO iPassQ = 1,nPasses'
       WRITE(LUMODA5,'(A)')'   px = Pcent(1)'
       WRITE(LUMODA5,'(A)')'   py = Pcent(2)'
       WRITE(LUMODA5,'(A)')'   pz = Pcent(3)'
       WRITE(LUMODA5,'(A)')'   pqx = px - Qcent(1,iPassQ)'
       WRITE(LUMODA5,'(A)')'   pqy = py - Qcent(2,iPassQ)'
       WRITE(LUMODA5,'(A)')'   pqz = pz - Qcent(3,iPassQ)'
       WRITE(LUMODA5,'(A)')'   squaredDistance = pqx*pqx+pqy*pqy+pqz*pqz'
       WRITE(LUMODA5,'(A)')'   WVAL = reducedExponents(1)*squaredDistance'

       do IA=1,5
          WRITE(IA,'(A)')'     !  0 < WVAL < 0.000001'
          WRITE(IA,'(A)')'!     IF (ABS(WVAL) .LT. SMALL) THEN'     
          WRITE(IA,'(A)')'!      RJ000 = D1'
          WRITE(IA,'(A)')'!     !  0 < WVAL < 12 '
          WRITE(IA,'(A)')'     IF (WVAL .LT. D12) THEN'
          WRITE(IA,'(A)')'      IPNT = NINT(D100*WVAL)'
          WRITE(IA,'(A)')'      WDIFF = WVAL - TENTH*IPNT'
          WRITE(IA,'(A)')'      W2    = WDIFF*WDIFF'
          WRITE(IA,'(A)')'      W3    = W2*WDIFF'
          WRITE(IA,'(A)')'      W2    = W2*D05'
          WRITE(IA,'(A)')'      W3    = W3*COEF3'
          WRITE(IA,'(A)')'      RJ000 = TABFJW(0,IPNT)-TABFJW(1,IPNT)*WDIFF+TABFJW(2,IPNT)*W2+TABFJW(3,IPNT)*W3'
          WRITE(IA,'(A)')'     !  12 < WVAL <= (2J+36) '
          WRITE(IA,'(A)')'     ELSE IF (WVAL.LE.D2JP36) THEN'
          WRITE(IA,'(A)')'      REXPW = D05*EXP(-WVAL)'
          WRITE(IA,'(A)')'      RWVAL = D1/WVAL'
          WRITE(IA,'(A)')'      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))'
          WRITE(IA,'(A)')'      RJ000 = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
          WRITE(IA,'(A)')'     !  (2J+36) < WVAL '
          WRITE(IA,'(A)')'     ELSE'
          WRITE(IA,'(A)')'      RJ000 = SQRT(PID4/WVAL)'
          WRITE(IA,'(A)')'     ENDIF'
       enddo
       WRITE(LUMODA1,'(A)')'     AUXarray(iPrimQ,iPrimP,iPassQ)=integralPrefactor(iPrimQ,iPrimP)*&'
       WRITE(LUMODA1,'(A)')'          & QpreExpFac(iPrimQ,iPassQ)*Pexpfac*RJ000'

       WRITE(LUMODA2,'(A)')'     AUXarray(iPrimP,iPassQ)=AUXarray(iPrimP,iPassQ) + integralPrefactor(iPrimQ,iPrimP)*&'
       WRITE(LUMODA2,'(A)')'          & QpreExpFac(iPrimQ,iPassQ)*Pexpfac*RJ000'

       WRITE(LUMODA3,'(A)')'     AUXarray(iPrimQ,iPassQ)=AUXarray(iPrimQ,iPassQ)+integralPrefactor(iPrimQ,iPrimP)*&'
       WRITE(LUMODA3,'(A)')'          & QpreExpFac(iPrimQ,iPassQ)*Pexpfac*RJ000'

       WRITE(LUMODA4,'(A)')'     AUXarray(iPassQ)=AUXarray(iPassQ) + integralPrefactor(iPrimQ,iPrimP)*&'
       WRITE(LUMODA4,'(A)')'          & QpreExpFac(iPrimQ,iPassQ)*Pexpfac*RJ000'
       do IA=1,4
          WRITE(IA,'(A)')'    enddo'
          WRITE(IA,'(A)')'   enddo'
          WRITE(IA,'(A)')'  enddo'
       enddo
       WRITE(LUMODA5,'(A)')'     AUXarray(iPassQ)=integralPrefactor(1)*QpreExpFac(iPassQ)*PpreExpFac(1)*RJ000'
       WRITE(LUMODA5,'(A)')'  enddo'

       WRITE(LUMODA1,'(A)')'end subroutine VerticalRecurrence'//ARCSTRING//'0'
       WRITE(LUMODA2,'(A)')'end subroutine VerticalRecurrence'//ARCSTRING//'SegQ0'
       WRITE(LUMODA3,'(A)')'end subroutine VerticalRecurrence'//ARCSTRING//'SegP0'
       WRITE(LUMODA4,'(A)')'end subroutine VerticalRecurrence'//ARCSTRING//'Seg0'
       WRITE(LUMODA5,'(A)')'end subroutine VerticalRecurrence'//ARCSTRING//'Seg1Prim0'

       !========================================================================================================
       !    VerticalRecurrence 1 
       !========================================================================================================

       do IA=1,20
          WRITE(IA,'(A)')''
       enddo
       WRITE(LUMODA1,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'1A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       WRITE(LUMODA2,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegQ1A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       WRITE(LUMODA3,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegP1A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       WRITE(LUMODA4,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       WRITE(LUMODA5,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1Prim1A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       do IA=1,5
          WRITE(IA,'(A)')'         & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,&'
          WRITE(IA,'(A)')'         & QpreExpFac,AUXarray)'
       enddo
       WRITE(LUMODB1,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'1B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       WRITE(LUMODB2,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegQ1B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       WRITE(LUMODB3,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegP1B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       WRITE(LUMODB4,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       WRITE(LUMODB5,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1Prim1B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       do IB=6,10
          WRITE(IB,'(A)')'         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,&'
          WRITE(IB,'(A)')'         & QpreExpFac,AUXarray)'
       enddo

       WRITE(LUMODC1,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'1C(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODC2,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegQ1C(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODC3,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegP1C(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODC4,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1C(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODC5,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1Prim1C(nPasses,nPrimP,nPrimQ,&'

       do IC=11,15
          WRITE(IC,'(A)')'         & reducedExponents,TABFJW,Qexp,Ccenter,Pcent,Qcent,&'
          WRITE(IC,'(A)')'         & integralPrefactor,PpreExpFac,QpreExpFac,AUXarray)'
       enddo

       WRITE(LUMODD1,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'1D(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODD2,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegQ1D(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODD3,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegP1D(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODD4,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1D(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMODD5,'(A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1Prim1D(nPasses,nPrimP,nPrimQ,&'
       do ID=16,20
          WRITE(ID,'(A)')'         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,&'
          WRITE(ID,'(A)')'         & integralPrefactor,PpreExpFac,QpreExpFac,AUXarray)'
       enddo

       do I=1,20
          WRITE(I,'(A)')'  implicit none'
          WRITE(I,'(A)')'  integer,intent(in) :: nPasses,nPrimP,nPrimQ'
       enddo
       !C and D
       !    do I=11,20
       !       WRITE(I,'(A)')'  integer,intent(in) :: nAtomsC,nAtomsD'
       !    enddo
       do I=1,20
          WRITE(I,'(A)')'  REAL(REALK),intent(in) :: TABFJW(0:4,0:1200)'
       enddo
       do J=1,8
          I = non1Prim(J)
          WRITE(I,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)'
       enddo
       WRITE(LUMODA5,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Pexp(1)'
       WRITE(LUMODB5,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Pexp(1)'
       do J=9,16
          I = non1Prim(J)
          WRITE(I,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)'
       enddo
       WRITE(LUMODC5,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Qexp(1)'
       WRITE(LUMODD5,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Qexp(1)'
       do J=1,16
          I = non1Prim(J)
          WRITE(I,'(A)')'  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)'
          WRITE(I,'(A)')'  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)'
       enddo
       WRITE(LUMODA5,'(A)')'  real(realk),intent(in) :: Pcent(3),Qcent(3,nPasses),integralPrefactor(1),QpreExpFac(nPasses),PpreExpFac(1)'
       WRITE(LUMODB5,'(A)')'  real(realk),intent(in) :: Pcent(3),Qcent(3,nPasses),integralPrefactor(1),QpreExpFac(nPasses),PpreExpFac(1)'
       WRITE(LUMODC5,'(A)')'  real(realk),intent(in) :: Pcent(3),Qcent(3,nPasses),integralPrefactor(1),QpreExpFac(nPasses),PpreExpFac(1)'
       WRITE(LUMODD5,'(A)')'  real(realk),intent(in) :: Pcent(3),Qcent(3,nPasses),integralPrefactor(1),QpreExpFac(nPasses),PpreExpFac(1)'
       do IA=1,5
          WRITE(IA,'(A)')'  real(realk),intent(in) :: Acenter(3)'
       enddo
       do IB=6,10
          WRITE(IB,'(A)')'  real(realk),intent(in) :: Bcenter(3)'
       enddo
       !    do IC=11,15
       !       WRITE(IC,'(A)')'  real(realk),intent(in) :: Ccenter(3,nAtomsC)'
       !    enddo
       !    do ID=16,20
       !       WRITE(ID,'(A)')'  real(realk),intent(in) :: Dcenter(3,nAtomsD)'
       !    enddo
       do IC=11,15
          WRITE(IC,'(A)')'  real(realk),intent(in) :: Ccenter(3)'
       enddo
       do ID=16,20
          WRITE(ID,'(A)')'  real(realk),intent(in) :: Dcenter(3)'
       enddo
       !standard
       do I=1,16,5
          WRITE(I,'(A)')'  real(realk),intent(inout) :: AUXarray(4,nPrimQ,nPrimP,nPasses)'
       enddo
       !segQ
       do I=2,17,5
          WRITE(I,'(A)')'  real(realk),intent(inout) :: AUXarray(4,nPrimP,nPasses)'
       enddo
       !segP
       do I=3,18,5
          WRITE(I,'(A)')'  real(realk),intent(inout) :: AUXarray(4,nPrimQ,nPasses)'
       enddo
       !seg
       do I=4,19,5
          WRITE(I,'(A)')'  real(realk),intent(inout) :: AUXarray(4,nPasses)'
       enddo
       !seg1Prim
       do I=5,20,5
          WRITE(I,'(A)')'  real(realk),intent(inout) :: AUXarray(4,nPasses)'
       enddo
       !==================================================================================================
       do I=1,20
          WRITE(I,'(A)')'  !local variables'
          WRITE(I,'(A)')'  integer :: iPassQ,iPrimP,iPrimQ,ipnt'
       enddo
       !    do I=11,15
       !       WRITE(I,'(A)')'  integer :: iAtomC'
       !    enddo
       !    do I=16,20
       !       WRITE(I,'(A)')'  integer :: iAtomD'
       !    enddo

       do I=1,5
          WRITE(I,'(A)')'  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa'
       enddo
       do I=6,10
          WRITE(I,'(A)')'  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb'
       enddo
       do I=11,15
          WRITE(I,'(A)')'  real(realk) :: Cx,Cy,Cz,Xqc,Yqc,Zqc'
       enddo
       do I=16,20
          WRITE(I,'(A)')'  real(realk) :: Dx,Dy,Dz,Xqd,Yqd,Zqd'
       enddo
       do I=1,10
          WRITE(I,'(A)')'  real(realk) :: mPX,mPY,mPZ,invexpP,alphaP,RJ000(0:1)'
       enddo
       do J=9,16
          I=non1Prim(J)
          WRITE(I,'(A)')'  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),alphaQ,RJ000(0:1)'
       enddo
       WRITE(LUMODC5,'(A)')'  real(realk) :: mPX,mPY,mPZ,invexpQ,alphaQ,RJ000(0:1)'
       WRITE(LUMODD5,'(A)')'  real(realk) :: mPX,mPY,mPZ,invexpQ,alphaQ,RJ000(0:1)'
       do I=1,20
          WRITE(I,'(A)')'  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq'
          WRITE(I,'(A)')'  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL' 
          WRITE(I,'(A)')'  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk'
          WRITE(I,'(A)')'  real(realk),parameter :: D2=2.0E0_realk'
          WRITE(I,'(A,ES24.16,A)')'  REAL(REALK),PARAMETER :: D2JP36=',2.0d0 + 36.0d0,'_realk'
          WRITE(I,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D03333=1.0E0_realk/3.0E0_realk'
          WRITE(I,'(A)')'  REAL(REALK),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk'
          WRITE(I,'(A)')'  REAL(REALK),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk'
          WRITE(I,'(A)')'  REAL(REALK),PARAMETER :: SMALL = 1E-15_realk,D12 = 12.0E0_realk'
          WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk'
          WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk'
          WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk'
          WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk'
          WRITE(I,'(A)')'  Real(realk), parameter :: PI=3.14159265358979323846E0_realk'
          WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk'
          WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2'
          WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI'
       enddo
       do I=1,5
          WRITE(I,'(A)')'  !ThetaAux(n,1,0,0) = Xpa*ThetaAux(n,0,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,0,0,0)'
       enddo
       do I=6,10
          WRITE(I,'(A)')'  !ThetaAux(n,1,0,0) = Xpb*ThetaAux(n,0,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,0,0,0)'
       enddo
       do I=11,15
          WRITE(I,'(A)')'  !ThetaAux(n,1,0,0) = Xqc*ThetaAux(n,0,0,0) + (-alpha/q*Xpq)*ThetaAux(n+1,0,0,0)'
       enddo
       do I=16,20
          WRITE(I,'(A)')'  !ThetaAux(n,1,0,0) = Xqd*ThetaAux(n,0,0,0) + (-alpha/q*Xpq)*ThetaAux(n+1,0,0,0)'
       enddo
       do I=1,20
          WRITE(I,'(A)')'  !i = 0 last 2 term vanish'
          WRITE(I,'(A)')'  !We include scaling of RJ000 '
       enddo
       do I=1,5
          WRITE(I,'(A)')'  Ax = -Acenter(1)'
          WRITE(I,'(A)')'  Ay = -Acenter(2)'
          WRITE(I,'(A)')'  Az = -Acenter(3)'
       enddo
       do I=6,10
          WRITE(I,'(A)')'  Bx = -Bcenter(1)'
          WRITE(I,'(A)')'  By = -Bcenter(2)'
          WRITE(I,'(A)')'  Bz = -Bcenter(3)'
       enddo
       do J=9,16
          I=non1Prim(J)
          WRITE(I,'(A)')'  DO iPrimQ=1, nPrimQ'
          WRITE(I,'(A)')'     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)'
          WRITE(I,'(A)')'  ENDDO'
       enddo
       WRITE(LUMODC5,'(A)')'  invexpQ = D1/Qexp(1)'
       WRITE(LUMODD5,'(A)')'  invexpQ = D1/Qexp(1)'
       do I=1,20
          WRITE(I,'(A)')'  DO iPassQ = 1,nPasses'
       enddo
       !seg
       do I=4,19,5
          WRITE(I,'(A)')'   AUXarray(1,iPassQ)=0.0E0_realk'
          WRITE(I,'(A)')'   AUXarray(2,iPassQ)=0.0E0_realk'
          WRITE(I,'(A)')'   AUXarray(3,iPassQ)=0.0E0_realk'
          WRITE(I,'(A)')'   AUXarray(4,iPassQ)=0.0E0_realk'
       enddo
       !segP
       do I=3,18,5
          WRITE(I,'(A)')'   DO iPrimQ=1, nPrimQ'
          WRITE(I,'(A)')'    AUXarray(1,iPrimQ,iPassQ)=0.0E0_realk'
          WRITE(I,'(A)')'    AUXarray(2,iPrimQ,iPassQ)=0.0E0_realk'
          WRITE(I,'(A)')'    AUXarray(3,iPrimQ,iPassQ)=0.0E0_realk'
          WRITE(I,'(A)')'    AUXarray(4,iPrimQ,iPassQ)=0.0E0_realk'
          WRITE(I,'(A)')'   ENDDO'
       enddo
       do I=11,15
          !       WRITE(I,'(A)')'   iAtomC = iPassQ - ((iPassQ-1)/nAtomsC)*nAtomsC'
          !       WRITE(I,'(A)')'   Cx = -Ccenter(1,iAtomC)'
          !       WRITE(I,'(A)')'   Cy = -Ccenter(2,iAtomC)'
          !       WRITE(I,'(A)')'   Cz = -Ccenter(3,iAtomC)'
          WRITE(I,'(A)')'   Cx = -Ccenter(1)'
          WRITE(I,'(A)')'   Cy = -Ccenter(2)'
          WRITE(I,'(A)')'   Cz = -Ccenter(3)'
       enddo
       do I=16,20
          !       WRITE(I,'(A)')'   iAtomD = (iPassQ-1)/nAtomsC+1'
          !       WRITE(I,'(A)')'   Dx = -Dcenter(1,iAtomD)'
          !       WRITE(I,'(A)')'   Dy = -Dcenter(2,iAtomD)'
          !       WRITE(I,'(A)')'   Dz = -Dcenter(3,iAtomD)'
          WRITE(I,'(A)')'   Dx = -Dcenter(1)'
          WRITE(I,'(A)')'   Dy = -Dcenter(2)'
          WRITE(I,'(A)')'   Dz = -Dcenter(3)'
       enddo
       do J=1,16
          I = non1Prim(J)
          WRITE(I,'(A)')'   DO iPrimP=1, nPrimP'
       enddo
       !segQ
       do I=2,17,5
          WRITE(I,'(A)')'    AUXarray(1,iPrimP,iPassQ)=0.0E0_realk'
          WRITE(I,'(A)')'    AUXarray(2,iPrimP,iPassQ)=0.0E0_realk'
          WRITE(I,'(A)')'    AUXarray(3,iPrimP,iPassQ)=0.0E0_realk'
          WRITE(I,'(A)')'    AUXarray(4,iPrimP,iPassQ)=0.0E0_realk'
       enddo
       do J=1,16
          I = non1Prim(J)
          WRITE(I,'(A)')'    Pexpfac = PpreExpFac(iPrimP)'
          WRITE(I,'(A)')'    mPX = -Pcent(1,iPrimP)'
          WRITE(I,'(A)')'    mPY = -Pcent(2,iPrimP)'
          WRITE(I,'(A)')'    mPZ = -Pcent(3,iPrimP)'
       enddo
       do J=1,8
          I = non1Prim(J)
          WRITE(I,'(A)')'    invexpP = D1/Pexp(iPrimP)'
       enddo
       do J=1,4
          I = pure1Prim(J)
          WRITE(I,'(A)')'    Pexpfac = PpreExpFac(1)'
          WRITE(I,'(A)')'    mPX = -Pcent(1)'
          WRITE(I,'(A)')'    mPY = -Pcent(2)'
          WRITE(I,'(A)')'    mPZ = -Pcent(3)'
       enddo
       do J=1,2
          I = pure1Prim(J)
          WRITE(I,'(A)')'    invexpP = D1/Pexp(1)'
       enddo
       do I=1,4
          WRITE(I,'(A)')'    Xpa = Pcent(1,iPrimP) + Ax'
          WRITE(I,'(A)')'    Ypa = Pcent(2,iPrimP) + Ay'
          WRITE(I,'(A)')'    Zpa = Pcent(3,iPrimP) + Az'
       enddo
       WRITE(LUMODA5,'(A)')'    Xpa = Pcent(1) + Ax'
       WRITE(LUMODA5,'(A)')'    Ypa = Pcent(2) + Ay'
       WRITE(LUMODA5,'(A)')'    Zpa = Pcent(3) + Az'
       do I=6,9
          WRITE(I,'(A)')'    Xpb = Pcent(1,iPrimP) + Bx'
          WRITE(I,'(A)')'    Ypb = Pcent(2,iPrimP) + By'
          WRITE(I,'(A)')'    Zpb = Pcent(3,iPrimP) + Bz'
       enddo
       WRITE(LUMODB5,'(A)')'    Xpb = Pcent(1) + Bx'
       WRITE(LUMODB5,'(A)')'    Ypb = Pcent(2) + By'
       WRITE(LUMODB5,'(A)')'    Zpb = Pcent(3) + Bz'

       do J=1,16
          I=non1Prim(J)
          WRITE(I,'(A)')'    DO iPrimQ=1, nPrimQ'
          WRITE(I,'(A)')'     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)'
          WRITE(I,'(A)')'     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)'
          WRITE(I,'(A)')'     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)'
       enddo
       do J=1,4
          I=pure1Prim(J)
          WRITE(I,'(A)')'     Xpq = mPX + Qcent(1,iPassQ)'
          WRITE(I,'(A)')'     Ypq = mPY + Qcent(2,iPassQ)'
          WRITE(I,'(A)')'     Zpq = mPZ + Qcent(3,iPassQ)'
       enddo
       do J=9,12
          I=non1Prim(J)
          WRITE(I,'(A)')'     Xqc = Qcent(1,iPrimQ,iPassQ) + Cx'
          WRITE(I,'(A)')'     Yqc = Qcent(2,iPrimQ,iPassQ) + Cy'
          WRITE(I,'(A)')'     Zqc = Qcent(3,iPrimQ,iPassQ) + Cz'
       enddo
       I=pure1Prim(3)
       WRITE(I,'(A)')'     Xqc = Qcent(1,iPassQ) + Cx'
       WRITE(I,'(A)')'     Yqc = Qcent(2,iPassQ) + Cy'
       WRITE(I,'(A)')'     Zqc = Qcent(3,iPassQ) + Cz'
       do J=13,16
          I=non1Prim(J)
          WRITE(I,'(A)')'     Xqd = Qcent(1,iPrimQ,iPassQ) + Dx'
          WRITE(I,'(A)')'     Yqd = Qcent(2,iPrimQ,iPassQ) + Dy'
          WRITE(I,'(A)')'     Zqd = Qcent(3,iPrimQ,iPassQ) + Dz'
       enddo
       I=pure1Prim(4)
       WRITE(I,'(A)')'     Xqd = Qcent(1,iPassQ) + Dx'
       WRITE(I,'(A)')'     Yqd = Qcent(2,iPassQ) + Dy'
       WRITE(I,'(A)')'     Zqd = Qcent(3,iPassQ) + Dz'

       do J=1,8
          I=non1Prim(J)
          WRITE(I,'(A)')'     alphaP = reducedExponents(iPrimQ,iPrimP)*invexpP'    
       enddo
       WRITE(LUMODA5,'(A)')'     alphaP = reducedExponents(1)*invexpP'    
       WRITE(LUMODB5,'(A)')'     alphaP = reducedExponents(1)*invexpP'    
       do J=9,16
          I=non1Prim(J)
          WRITE(I,'(A)')'     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)'
       enddo
       WRITE(LUMODC5,'(A)')'     alphaQ = -reducedExponents(1)*invexpQ'
       WRITE(LUMODD5,'(A)')'     alphaQ = -reducedExponents(1)*invexpQ'
       do I=1,10
          WRITE(I,'(A)')'     alphaXpq = alphaP*Xpq'
          WRITE(I,'(A)')'     alphaYpq = alphaP*Ypq'
          WRITE(I,'(A)')'     alphaZpq = alphaP*Zpq'
       enddo
       do I=11,20
          WRITE(I,'(A)')'     alphaXpq = alphaQ*Xpq'
          WRITE(I,'(A)')'     alphaYpq = alphaQ*Ypq'
          WRITE(I,'(A)')'     alphaZpq = alphaQ*Zpq'
       enddo
       do I=1,200
          WRITE(I,'(A)')'     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq'
       enddo
       do J=1,16
          I=non1Prim(J)
          WRITE(I,'(A)')'     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance'
       enddo
       do J=1,4
          I=pure1Prim(J)
          WRITE(I,'(A)')'     WVAL = reducedExponents(1)*squaredDistance'
       enddo
       do I=1,20
          !       WRITE(I,'(A)')'     !  0 < WVAL < 0.000001'
          !       WRITE(I,'(A)')'!     IF (ABS(WVAL) .LT. SMALL) THEN'     
          !       WRITE(I,'(A)')'!      RJ000(0) = D1'
          !       WRITE(I,'(A)')'!      RJ000(1)= D03333 !THE BOYS FUNCTION FOR ZERO ARGUMENT'
          !       WRITE(I,'(A)')'!     !  0 < WVAL < 12 '
          WRITE(I,'(A)')'     IF (WVAL .LT. D12) THEN'
          WRITE(I,'(A)')'      IPNT = NINT(D100*WVAL)'
          WRITE(I,'(A)')'      WDIFF = WVAL - TENTH*IPNT'
          WRITE(I,'(A)')'      W2    = WDIFF*WDIFF'
          WRITE(I,'(A)')'      W3    = W2*WDIFF'
          WRITE(I,'(A)')'      W2    = W2*D05'
          WRITE(I,'(A)')'      W3    = W3*COEF3'
          WRITE(I,'(A)')'      RJ000(0)=TABFJW(0,IPNT)-TABFJW(1,IPNT)*WDIFF+TABFJW(2,IPNT)*W2+TABFJW(3,IPNT)*W3'
          WRITE(I,'(A)')'      RJ000(1)=TABFJW(1,IPNT)-TABFJW(2,IPNT)*WDIFF+TABFJW(3,IPNT)*W2+TABFJW(4,IPNT)*W3'
          WRITE(I,'(A)')'     !  12 < WVAL <= (2J+36) '
          WRITE(I,'(A)')'     ELSE IF (WVAL.LE.D2JP36) THEN'
          WRITE(I,'(A)')'      REXPW = D05*EXP(-WVAL)'
          WRITE(I,'(A)')'      RWVAL = D1/WVAL'
          WRITE(I,'(A)')'      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))'
          WRITE(I,'(A)')'      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
          WRITE(I,'(A)')'      RJ000(1) = RWVAL*(D05*RJ000(0)-REXPW)'
          WRITE(I,'(A)')'     !  (2J+36) < WVAL '
          WRITE(I,'(A)')'     ELSE'
          WRITE(I,'(A)')'      RWVAL = PID4/WVAL'
          WRITE(I,'(A)')'      RJ000(0) = SQRT(RWVAL)'
          WRITE(I,'(A)')'      RJ000(1) = RWVAL*PID4I*D05*RJ000(0)'
          WRITE(I,'(A)')'     ENDIF'
       enddo
       do J=1,16
          I=non1Prim(J)
          WRITE(I,'(A)')'     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac'
       enddo
       do J=1,4
          I=pure1Prim(J)
          WRITE(I,'(A)')'     PREF = integralPrefactor(1)*QpreExpFac(iPassQ)*PpreExpFac(1)'
       enddo
       do I=1,20
          WRITE(I,'(A)')'     TMP1 = PREF*RJ000(0)'
          WRITE(I,'(A)')'     TMP2 = PREF*RJ000(1)'
       enddo
       !standard
       do I=1,16,5
          WRITE(I,'(A)')'     AUXarray(1,iPrimQ,iPrimP,iPassQ) = TMP1'
       enddo
       !segQ
       do I=2,17,5
          WRITE(I,'(A)')'     AUXarray(1,iPrimP,iPassQ) = AUXarray(1,iPrimP,iPassQ) + TMP1'
       enddo
       !segP
       do I=3,18,5
          WRITE(I,'(A)')'     AUXarray(1,iPrimQ,iPassQ) = AUXarray(1,iPrimQ,iPassQ) + TMP1'
       enddo
       !seg
       do I=4,19,5
          WRITE(I,'(A)')'     AUXarray(1,iPassQ) = AUXarray(1,iPassQ) + TMP1'
       enddo
       !seg1Prim
       do I=5,20,5
          WRITE(I,'(A)')'     AUXarray(1,iPassQ) = TMP1'
       enddo

       WRITE(LUMODA1,'(A)')'     AUXarray(2,iPrimQ,iPrimP,iPassQ) = Xpa*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODA1,'(A)')'     AUXarray(3,iPrimQ,iPrimP,iPassQ) = Ypa*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODA1,'(A)')'     AUXarray(4,iPrimQ,iPrimP,iPassQ) = Zpa*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODA2,'(A)')'     AUXarray(2,iPrimP,iPassQ) = AUXarray(2,iPrimP,iPassQ) + Xpa*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODA2,'(A)')'     AUXarray(3,iPrimP,iPassQ) = AUXarray(3,iPrimP,iPassQ) + Ypa*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODA2,'(A)')'     AUXarray(4,iPrimP,iPassQ) = AUXarray(4,iPrimP,iPassQ) + Zpa*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODA3,'(A)')'     AUXarray(2,iPrimQ,iPassQ) = AUXarray(2,iPrimQ,iPassQ) + Xpa*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODA3,'(A)')'     AUXarray(3,iPrimQ,iPassQ) = AUXarray(3,iPrimQ,iPassQ) + Ypa*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODA3,'(A)')'     AUXarray(4,iPrimQ,iPassQ) = AUXarray(4,iPrimQ,iPassQ) + Zpa*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODA4,'(A)')'     AUXarray(2,iPassQ) = AUXarray(2,iPassQ) + Xpa*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODA4,'(A)')'     AUXarray(3,iPassQ) = AUXarray(3,iPassQ) + Ypa*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODA4,'(A)')'     AUXarray(4,iPassQ) = AUXarray(4,iPassQ) + Zpa*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODA5,'(A)')'     AUXarray(2,iPassQ) = Xpa*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODA5,'(A)')'     AUXarray(3,iPassQ) = Ypa*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODA5,'(A)')'     AUXarray(4,iPassQ) = Zpa*TMP1 + alphaZpq*TMP2'

       WRITE(LUMODB1,'(A)')'     AUXarray(2,iPrimQ,iPrimP,iPassQ) = Xpb*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODB1,'(A)')'     AUXarray(3,iPrimQ,iPrimP,iPassQ) = Ypb*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODB1,'(A)')'     AUXarray(4,iPrimQ,iPrimP,iPassQ) = Zpb*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODB2,'(A)')'     AUXarray(2,iPrimP,iPassQ) = AUXarray(2,iPrimP,iPassQ) + Xpb*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODB2,'(A)')'     AUXarray(3,iPrimP,iPassQ) = AUXarray(3,iPrimP,iPassQ) + Ypb*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODB2,'(A)')'     AUXarray(4,iPrimP,iPassQ) = AUXarray(4,iPrimP,iPassQ) + Zpb*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODB3,'(A)')'     AUXarray(2,iPrimQ,iPassQ) = AUXarray(2,iPrimQ,iPassQ) + Xpb*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODB3,'(A)')'     AUXarray(3,iPrimQ,iPassQ) = AUXarray(3,iPrimQ,iPassQ) + Ypb*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODB3,'(A)')'     AUXarray(4,iPrimQ,iPassQ) = AUXarray(4,iPrimQ,iPassQ) + Zpb*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODB4,'(A)')'     AUXarray(2,iPassQ) = AUXarray(2,iPassQ) + Xpb*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODB4,'(A)')'     AUXarray(3,iPassQ) = AUXarray(3,iPassQ) + Ypb*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODB4,'(A)')'     AUXarray(4,iPassQ) = AUXarray(4,iPassQ) + Zpb*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODB5,'(A)')'     AUXarray(2,iPassQ) = Xpb*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODB5,'(A)')'     AUXarray(3,iPassQ) = Ypb*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODB5,'(A)')'     AUXarray(4,iPassQ) = Zpb*TMP1 + alphaZpq*TMP2'

       WRITE(LUMODC1,'(A)')'     AUXarray(2,iPrimQ,iPrimP,iPassQ) = Xqc*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODC1,'(A)')'     AUXarray(3,iPrimQ,iPrimP,iPassQ) = Yqc*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODC1,'(A)')'     AUXarray(4,iPrimQ,iPrimP,iPassQ) = Zqc*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODC2,'(A)')'     AUXarray(2,iPrimP,iPassQ) = AUXarray(2,iPrimP,iPassQ) + Xqc*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODC2,'(A)')'     AUXarray(3,iPrimP,iPassQ) = AUXarray(3,iPrimP,iPassQ) + Yqc*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODC2,'(A)')'     AUXarray(4,iPrimP,iPassQ) = AUXarray(4,iPrimP,iPassQ) + Zqc*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODC3,'(A)')'     AUXarray(2,iPrimQ,iPassQ) = AUXarray(2,iPrimQ,iPassQ) + Xqc*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODC3,'(A)')'     AUXarray(3,iPrimQ,iPassQ) = AUXarray(3,iPrimQ,iPassQ) + Yqc*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODC3,'(A)')'     AUXarray(4,iPrimQ,iPassQ) = AUXarray(4,iPrimQ,iPassQ) + Zqc*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODC4,'(A)')'     AUXarray(2,iPassQ) = AUXarray(2,iPassQ) + Xqc*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODC4,'(A)')'     AUXarray(3,iPassQ) = AUXarray(3,iPassQ) + Yqc*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODC4,'(A)')'     AUXarray(4,iPassQ) = AUXarray(4,iPassQ) + Zqc*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODC5,'(A)')'     AUXarray(2,iPassQ) = Xqc*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODC5,'(A)')'     AUXarray(3,iPassQ) = Yqc*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODC5,'(A)')'     AUXarray(4,iPassQ) = Zqc*TMP1 + alphaZpq*TMP2'

       WRITE(LUMODD1,'(A)')'     AUXarray(2,iPrimQ,iPrimP,iPassQ) = Xqd*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODD1,'(A)')'     AUXarray(3,iPrimQ,iPrimP,iPassQ) = Yqd*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODD1,'(A)')'     AUXarray(4,iPrimQ,iPrimP,iPassQ) = Zqd*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODD2,'(A)')'     AUXarray(2,iPrimP,iPassQ) = AUXarray(2,iPrimP,iPassQ) + Xqd*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODD2,'(A)')'     AUXarray(3,iPrimP,iPassQ) = AUXarray(3,iPrimP,iPassQ) + Yqd*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODD2,'(A)')'     AUXarray(4,iPrimP,iPassQ) = AUXarray(4,iPrimP,iPassQ) + Zqd*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODD3,'(A)')'     AUXarray(2,iPrimQ,iPassQ) = AUXarray(2,iPrimQ,iPassQ) + Xqd*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODD3,'(A)')'     AUXarray(3,iPrimQ,iPassQ) = AUXarray(3,iPrimQ,iPassQ) + Yqd*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODD3,'(A)')'     AUXarray(4,iPrimQ,iPassQ) = AUXarray(4,iPrimQ,iPassQ) + Zqd*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODD4,'(A)')'     AUXarray(2,iPassQ) = AUXarray(2,iPassQ) + Xqd*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODD4,'(A)')'     AUXarray(3,iPassQ) = AUXarray(3,iPassQ) + Yqd*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODD4,'(A)')'     AUXarray(4,iPassQ) = AUXarray(4,iPassQ) + Zqd*TMP1 + alphaZpq*TMP2'
       WRITE(LUMODD5,'(A)')'     AUXarray(2,iPassQ) = Xqd*TMP1 + alphaXpq*TMP2'
       WRITE(LUMODD5,'(A)')'     AUXarray(3,iPassQ) = Yqd*TMP1 + alphaYpq*TMP2'
       WRITE(LUMODD5,'(A)')'     AUXarray(4,iPassQ) = Zqd*TMP1 + alphaZpq*TMP2'


       do J=1,16
          I=non1Prim(J)       
          WRITE(I,'(A)')'    enddo'
          WRITE(I,'(A)')'   enddo'
       enddo
       do I=1,20
          WRITE(I,'(A)')'  enddo'
          WRITE(I,'(A)')'end subroutine'
       enddo


       !============================================================================================================
       !         VerticalRecurrence GENERAL
       !============================================================================================================

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

          do IA=1,20
             WRITE(IA,'(A)')''
          enddo
          WRITE(LUMODA1,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING,JMAX,'A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
          WRITE(LUMODA2,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegQ',JMAX,'A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
          WRITE(LUMODA3,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegP',JMAX,'A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
          WRITE(LUMODA4,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg',JMAX,'A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
          WRITE(LUMODA5,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1Prim',JMAX,'A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
          do I=1,5
             WRITE(I,'(A)')'         & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&'
             WRITE(I,'(A)')'         & AUXarray)'
          enddo

          WRITE(LUMODB1,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING,JMAX,'B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
          WRITE(LUMODB2,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegQ',JMAX,'B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
          WRITE(LUMODB3,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegP',JMAX,'B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
          WRITE(LUMODB4,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg',JMAX,'B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
          WRITE(LUMODB5,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1Prim',JMAX,'B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
          do I=6,10
             WRITE(I,'(A)')'         & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,&'
             WRITE(I,'(A)')'         & AUXarray)'
          enddo

          WRITE(LUMODC1,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING,JMAX,'C(nPasses,nPrimP,nPrimQ,&'
          WRITE(LUMODC2,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegQ',JMAX,'C(nPasses,nPrimP,nPrimQ,&'
          WRITE(LUMODC3,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegP',JMAX,'C(nPasses,nPrimP,nPrimQ,&'
          WRITE(LUMODC4,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg',JMAX,'C(nPasses,nPrimP,nPrimQ,&'
          WRITE(LUMODC5,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1Prim',JMAX,'C(nPasses,nPrimP,nPrimQ,&'
          do I=11,15
             WRITE(I,'(A)')'         & reducedExponents,TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,&'
             WRITE(I,'(A)')'         & PpreExpFac,QpreExpFac,AUXarray)'
          enddo

          WRITE(LUMODD1,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING,JMAX,'D(nPasses,nPrimP,nPrimQ,&'
          WRITE(LUMODD2,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegQ',JMAX,'D(nPasses,nPrimP,nPrimQ,&'
          WRITE(LUMODD3,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'SegP',JMAX,'D(nPasses,nPrimP,nPrimQ,&'
          WRITE(LUMODD4,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg',JMAX,'D(nPasses,nPrimP,nPrimQ,&'
          WRITE(LUMODD5,'(A,I1,A)')'subroutine VerticalRecurrence'//ARCSTRING//'Seg1Prim',JMAX,'D(nPasses,nPrimP,nPrimQ,&'
          do I=16,20
             WRITE(I,'(A)')'         & reducedExponents,TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,&'
             WRITE(I,'(A)')'         & PpreExpFac,QpreExpFac,AUXarray)'
          enddo
          do I=1,20
             WRITE(I,'(A)')'  implicit none'
             WRITE(I,'(A)')'  integer,intent(in) :: nPasses,nPrimP,nPrimQ'
          enddo
          !C and D
          !       do I=11,20
          !          WRITE(I,'(A)')'  integer,intent(in) :: nAtomsC,nAtomsD'
          !       enddo
          do I=1,20
             WRITE(I,'(A,I2,A)')'  REAL(REALK),intent(in) :: TABFJW(0:',JMAX+3,',0:1200)'
          enddo
          do J=1,8
             I = non1Prim(J)
             WRITE(I,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP)'
          enddo
          WRITE(LUMODA5,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Pexp(1)'
          WRITE(LUMODB5,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Pexp(1)'
          do J=9,16
             I = non1Prim(J)
             WRITE(I,'(A)')'  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Qexp(nPrimQ)'
          enddo
          WRITE(LUMODC5,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Qexp(1)'
          WRITE(LUMODD5,'(A)')'  real(realk),intent(in) :: reducedExponents(1),Qexp(1)'
          do J=1,16
             I = non1Prim(J)
             WRITE(I,'(A)')'  real(realk),intent(in) :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)'
             WRITE(I,'(A)')'  real(realk),intent(in) :: integralPrefactor(nprimQ,nPrimP),QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)'
          enddo
          WRITE(LUMODA5,'(A)')'  real(realk),intent(in) :: Pcent(3),Qcent(3,nPasses),integralPrefactor(1),QpreExpFac(nPasses),PpreExpFac(1)'
          WRITE(LUMODB5,'(A)')'  real(realk),intent(in) :: Pcent(3),Qcent(3,nPasses),integralPrefactor(1),QpreExpFac(nPasses),PpreExpFac(1)'
          WRITE(LUMODC5,'(A)')'  real(realk),intent(in) :: Pcent(3),Qcent(3,nPasses),integralPrefactor(1),QpreExpFac(nPasses),PpreExpFac(1)'
          WRITE(LUMODD5,'(A)')'  real(realk),intent(in) :: Pcent(3),Qcent(3,nPasses),integralPrefactor(1),QpreExpFac(nPasses),PpreExpFac(1)'
          do IA=1,5
             WRITE(IA,'(A)')'  real(realk),intent(in) :: Acenter(3)'
          enddo
          do IB=6,10    
             WRITE(IB,'(A)')'  real(realk),intent(in) :: Bcenter(3)'
          enddo
          !       do IC=11,15
          !          WRITE(IC,'(A)')'  real(realk),intent(in) :: Ccenter(3,nAtomsC)'
          !       enddo
          !       do ID=16,20
          !          WRITE(ID,'(A)')'  real(realk),intent(in) :: Dcenter(3,nAtomsD)'
          !       enddo
          do IC=11,15
             WRITE(IC,'(A)')'  real(realk),intent(in) :: Ccenter(3)'
          enddo
          do ID=16,20
             WRITE(ID,'(A)')'  real(realk),intent(in) :: Dcenter(3)'
          enddo
          do I=1,16,5
             WRITE(I,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPrimQ*nPrimP*nPasses)'
          enddo
          !segQ
          do I=2,17,5
             WRITE(I,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPrimP*nPasses)'
          enddo
          !segP
          do I=3,18,5
             WRITE(I,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPrimQ*nPasses)'
          enddo
          !seg
          do I=4,19,5
             WRITE(I,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPasses)'
          enddo
          !seg1Prim
          do I=5,20,5
             WRITE(I,'(A,I5,A)')'  real(realk),intent(inout) :: AUXarray(',nTUV,',nPasses)'
          enddo
          !==================================================================================================
          do I=1,20
             WRITE(I,'(A)')'  !local variables'
             WRITE(I,'(A)')'  integer :: iPassQ,iPrimP,iPrimQ,ipnt,IP,iTUV'
          enddo
          !       do I=11,15
          !          WRITE(I,'(A)')'  integer :: iAtomC'
          !       enddo
          !       do I=16,20
          !          WRITE(I,'(A)')'  integer :: iAtomD'
          !       enddo
          do I=2,17,5
             WRITE(I,'(A,I5,A)')'  real(realk) :: TMPAUXarray(',nTUVprev,')'
          enddo
          do I=3,18,5
             WRITE(I,'(A,I5,A)')'  real(realk) :: TMPAUXarray(',nTUVprev,')'
          enddo
          do I=4,19,5
             WRITE(I,'(A,I5,A)')'  real(realk) :: TMPAUXarray(',nTUVprev,')'
          enddo
          do I=5,20,5
             WRITE(I,'(A,I5,A)')'  real(realk) :: TMPAUXarray(',nTUVprev,')'
          enddo
          do I=1,5
             WRITE(I,'(A)')'  real(realk) :: Ax,Ay,Az,Xpa,Ypa,Zpa'
          enddo
          do I=6,10
             WRITE(I,'(A)')'  real(realk) :: Bx,By,Bz,Xpb,Ypb,Zpb'
          enddo
          do I=11,15
             WRITE(I,'(A)')'  real(realk) :: Cx,Cy,Cz,Xqc,Yqc,Zqc'
          enddo
          do I=16,20
             WRITE(I,'(A)')'  real(realk) :: Dx,Dy,Dz,Xqd,Yqd,Zqd'
          enddo
          do I=1,10
             WRITE(I,'(A,i2,A)')'  real(realk) :: mPX,mPY,mPZ,invexpP,inv2expP,alphaP,RJ000(0:',JMAX,')'
          enddo
          do J=9,16
             I=non1Prim(J)
             WRITE(I,'(A,i2,A)')'  real(realk) :: mPX,mPY,mPZ,invexpQ(nPrimQ),inv2expQ,alphaQ,RJ000(0:',JMAX,')'
          enddo
          WRITE(LUMODC5,'(A,i2,A)')'  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ,RJ000(0:',JMAX,')'
          WRITE(LUMODD5,'(A,i2,A)')'  real(realk) :: mPX,mPY,mPZ,invexpQ,inv2expQ,alphaQ,RJ000(0:',JMAX,')'
          do I=1,20
             WRITE(I,'(A,i4,A)')'  real(realk) :: TwoTerms(',MAX(1,nTUVprev2-nTUVprev3),')'
             WRITE(I,'(A)')'  real(realk) :: Pexpfac,PREF,TMP1,TMP2,Xpq,Ypq,Zpq,alphaXpq,alphaYpq,alphaZpq'
             WRITE(I,'(A)')'  real(realk) :: squaredDistance,WVAL,WDIFF,W2,W3,REXPW,RWVAL,GVAL' 
             WRITE(I,'(A)')'  REAL(REALK),PARAMETER :: TENTH = 0.01E0_realk,D05 =0.5E0_realk'
             WRITE(I,'(A)')'  real(realk),parameter :: D2=2.0E0_realk'
             WRITE(I,'(A,ES24.16,A)')'  REAL(REALK),PARAMETER :: D2JP36=',2.d0*JMAX + 36.d0,'_realk'
             WRITE(I,'(A)')'  real(realk),parameter :: D1=1.0E0_realk,D03333=1.0E0_realk/3.0E0_realk'
             WRITE(I,'(A)')'  REAL(REALK),PARAMETER :: D4 = 4E0_realk, D100=100E0_realk'
             WRITE(I,'(A)')'  REAL(REALK),PARAMETER :: COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk'
             WRITE(I,'(A)')'  REAL(REALK),PARAMETER :: SMALL = 1E-15_realk,D12 = 12.0E0_realk'
             WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk'
             WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk'
             WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk'
             WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk'
             WRITE(I,'(A)')'  Real(realk), parameter :: PI=3.14159265358979323846E0_realk'
             WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk'
             WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2'
             WRITE(I,'(A)')'  REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI'
          enddo
          do I=1,20
             C = JMAX+2
             DO JTMP=1,JMAX
                C = C-1
                nTUVTMPprev=(JTMP-1)*(JTMP)*(JTMP+1)/6
                nTUVTMP=(JTMP)*(JTMP+1)*(JTMP+2)/6
                if(JTMP.LT.10)THEN
                   WRITE(I,'(A,I1,A,I3,A,I3,A,I1,A)')'  real(realk) :: TMParray',JTMP,'(',nTUVTMPprev+1,':',nTUVTMP,',2:',C,')'
                else
                   WRITE(I,'(A,I2,A,I3,A,I3,A,I1,A)')'  real(realk) :: TMParray',JTMP,'(',nTUVTMPprev+1,':',nTUVTMP,',2:',C,')'

                endif
             ENDDO
          enddo
          do I=1,5
             WRITE(I,'(A)')'  !TUV(T,0,0,N) = Xpa*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)'
             WRITE(I,'(A)')'  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))'
          enddo
          do I=6,10
             WRITE(I,'(A)')'  !TUV(T,0,0,N) = Xpb*TUV(T-1,0,0,N)-(alpha/p)*Xpq*TUV(T-1,0,0,N+1)'
             WRITE(I,'(A)')'  !             + T/(2p)*(TUV(T-2,0,0,N)-(alpha/p)*TUV(T-2,0,0,N+1))'
          enddo
          do I=11,15
             WRITE(I,'(A)')'  !TUV(T,0,0,N) = Xqc*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)'
             WRITE(I,'(A)')'  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))'
          enddo
          do I=16,20
             WRITE(I,'(A)')'  !TUV(T,0,0,N) = Xqd*TUV(T-1,0,0,N)-(alpha/q)*Xpq*TUV(T-1,0,0,N+1)'
             WRITE(I,'(A)')'  !             + T/(2q)*(TUV(T-2,0,0,N)-(alpha/q)*TUV(T-2,0,0,N+1))'
          enddo
          do I=1,20
             WRITE(I,'(A)')'  !We include scaling of RJ000 '
          enddo
          do I=1,5
             WRITE(I,'(A)')'  Ax = -Acenter(1)'
             WRITE(I,'(A)')'  Ay = -Acenter(2)'
             WRITE(I,'(A)')'  Az = -Acenter(3)'
          enddo
          do I=6,10
             WRITE(I,'(A)')'  Bx = -Bcenter(1)'
             WRITE(I,'(A)')'  By = -Bcenter(2)'
             WRITE(I,'(A)')'  Bz = -Bcenter(3)'
          enddo
          do J=9,16
             I=non1Prim(J)
             WRITE(I,'(A)')'  DO iPrimQ=1, nPrimQ'
             WRITE(I,'(A)')'     invexpQ(iPrimQ) = D1/Qexp(iPrimQ)'
             WRITE(I,'(A)')'  ENDDO'
          enddo
          WRITE(LUMODC5,'(A)')'  invexpQ = D1/Qexp(1)'
          WRITE(LUMODD5,'(A)')'  invexpQ = D1/Qexp(1)'
          ! ======================================================================
          !    iPassQ Loop 
          ! ======================================================================
          do I=1,20
             WRITE(I,'(A)')'  DO iPassQ = 1,nPasses'
          enddo
          !general
          do I=1,16,5
             WRITE(I,'(A)')'   iP = (iPassQ-1)*nPrimQ*nPrimP'
          enddo
          !segQ
          do I=2,17,5
             WRITE(I,'(A)')'   iP = (iPassQ-1)*nPrimP'
          enddo
          !segP
          do I=3,18,5
             WRITE(I,'(A)')'   iP = (iPassQ-1)*nPrimQ'
             WRITE(I,'(A)')'   DO iPrimQ=1, nPrimQ'
             WRITE(I,'(A)')'    iP = iP + 1'
             WRITE(I,'(A,i5)')'    DO iTUV=1,',nTUV
             WRITE(I,'(A)')'     AUXarray(iTUV,iP)=0.0E0_realk'
             WRITE(I,'(A)')'    ENDDO'
             WRITE(I,'(A)')'   ENDDO'
             WRITE(I,'(A)')'   iP = (iPassQ-1)*nPrimQ'
          enddo
          !seg
          do I=4,19,5
             WRITE(I,'(A)')'   iP = iPassQ'
             WRITE(I,'(A,i5)')'   DO iTUV=1,',nTUV
             WRITE(I,'(A)')'    AUXarray(iTUV,iPassQ)=0.0E0_realk'
             WRITE(I,'(A)')'   ENDDO'
          enddo
          do I=5,20,5
             WRITE(I,'(A)')'   iP = iPassQ'
          enddo
          do I=11,15
             !          WRITE(I,'(A)')'   iAtomC = iPassQ - ((iPassQ-1)/nAtomsC)*nAtomsC'
             !          WRITE(I,'(A)')'   Cx = -Ccenter(1,iAtomC)'
             !          WRITE(I,'(A)')'   Cy = -Ccenter(2,iAtomC)'
             !          WRITE(I,'(A)')'   Cz = -Ccenter(3,iAtomC)'
             WRITE(I,'(A)')'   Cx = -Ccenter(1)'
             WRITE(I,'(A)')'   Cy = -Ccenter(2)'
             WRITE(I,'(A)')'   Cz = -Ccenter(3)'
          enddo
          do I=16,20
             !          WRITE(I,'(A)')'   iAtomD = (iPassQ-1)/nAtomsC+1'
             !          WRITE(I,'(A)')'   Dx = -Dcenter(1,iAtomD)'
             !          WRITE(I,'(A)')'   Dy = -Dcenter(2,iAtomD)'
             !          WRITE(I,'(A)')'   Dz = -Dcenter(3,iAtomD)'
             WRITE(I,'(A)')'   Dx = -Dcenter(1)'
             WRITE(I,'(A)')'   Dy = -Dcenter(2)'
             WRITE(I,'(A)')'   Dz = -Dcenter(3)'
          enddo
          ! ======================================================================
          !    iPrimP Loop 
          ! ======================================================================
          do J=1,16
             I = non1Prim(J)
             WRITE(I,'(A)')'   DO iPrimP=1, nPrimP'
          enddo
          !segP
          do I=3,18,5
             WRITE(I,'(A)')'    iP = (iPassQ-1)*nPrimQ'
          enddo
          !segQ
          do I=2,17,5
             WRITE(I,'(A)')'    iP = iP + 1'
             WRITE(I,'(A,i5)')'    DO iTUV=1,',nTUV
             WRITE(I,'(A)')'     AUXarray(iTUV,iP)=0.0E0_realk'
             WRITE(I,'(A)')'    ENDDO'
          enddo
          do J=1,16
             I = non1Prim(J)
             WRITE(I,'(A)')'    Pexpfac = PpreExpFac(iPrimP)'
             WRITE(I,'(A)')'    mPX = -Pcent(1,iPrimP)'
             WRITE(I,'(A)')'    mPY = -Pcent(2,iPrimP)'
             WRITE(I,'(A)')'    mPZ = -Pcent(3,iPrimP)'
          enddo
          do J=1,8
             I = non1Prim(J)
             WRITE(I,'(A)')'    invexpP = D1/Pexp(iPrimP)'
             WRITE(I,'(A)')'    inv2expP = D05*invexpP'
          enddo
          do J=1,4
             I = pure1Prim(J)
             WRITE(I,'(A)')'    Pexpfac = PpreExpFac(1)'
             WRITE(I,'(A)')'    mPX = -Pcent(1)'
             WRITE(I,'(A)')'    mPY = -Pcent(2)'
             WRITE(I,'(A)')'    mPZ = -Pcent(3)'
          enddo
          do J=1,2
             I = pure1Prim(J)
             WRITE(I,'(A)')'    invexpP = D1/Pexp(1)'
             WRITE(I,'(A)')'    inv2expP = D05*invexpP'
          enddo

          do I=1,4
             WRITE(I,'(A)')'    Xpa = Pcent(1,iPrimP) + Ax'
             WRITE(I,'(A)')'    Ypa = Pcent(2,iPrimP) + Ay'
             WRITE(I,'(A)')'    Zpa = Pcent(3,iPrimP) + Az'
          enddo
          WRITE(LUMODA5,'(A)')'    Xpa = Pcent(1) + Ax'
          WRITE(LUMODA5,'(A)')'    Ypa = Pcent(2) + Ay'
          WRITE(LUMODA5,'(A)')'    Zpa = Pcent(3) + Az'
          do I=6,9
             WRITE(I,'(A)')'    Xpb = Pcent(1,iPrimP) + Bx'
             WRITE(I,'(A)')'    Ypb = Pcent(2,iPrimP) + By'
             WRITE(I,'(A)')'    Zpb = Pcent(3,iPrimP) + Bz'
          enddo
          WRITE(LUMODB5,'(A)')'    Xpb = Pcent(1) + Bx'
          WRITE(LUMODB5,'(A)')'    Ypb = Pcent(2) + By'
          WRITE(LUMODB5,'(A)')'    Zpb = Pcent(3) + Bz'

          ! ======================================================================
          !    iPrimQ Loop 
          ! ======================================================================
          do J=1,16
             I=non1Prim(J)
             WRITE(I,'(A)')'    DO iPrimQ=1, nPrimQ'
          enddo
          !if Seg or seg1prim iP not used
          !if SegQ then IP should be increased after iPrimP only 
          !if Gen then IP should be increased after iPrimQ 
          !if SegP then IP should be increased after iPrimQ 
          !general
          do I=1,16,5
             WRITE(I,'(A)')'     iP = iP + 1'
          enddo
          !segP
          do I=3,18,5
             WRITE(I,'(A)')'   iP = iP + 1'
          enddo

          do J=1,16
             I=non1Prim(J)
             WRITE(I,'(A)')'     Xpq = mPX + Qcent(1,iPrimQ,iPassQ)'
             WRITE(I,'(A)')'     Ypq = mPY + Qcent(2,iPrimQ,iPassQ)'
             WRITE(I,'(A)')'     Zpq = mPZ + Qcent(3,iPrimQ,iPassQ)'
          enddo
          do J=1,4
             I=pure1Prim(J)
             WRITE(I,'(A)')'     Xpq = mPX + Qcent(1,iPassQ)'
             WRITE(I,'(A)')'     Ypq = mPY + Qcent(2,iPassQ)'
             WRITE(I,'(A)')'     Zpq = mPZ + Qcent(3,iPassQ)'
          enddo
          do J=9,12
             I=non1Prim(J)
             WRITE(I,'(A)')'     Xqc = Qcent(1,iPrimQ,iPassQ) + Cx'
             WRITE(I,'(A)')'     Yqc = Qcent(2,iPrimQ,iPassQ) + Cy'
             WRITE(I,'(A)')'     Zqc = Qcent(3,iPrimQ,iPassQ) + Cz'
          enddo
          I=pure1Prim(3)
          WRITE(I,'(A)')'     Xqc = Qcent(1,iPassQ) + Cx'
          WRITE(I,'(A)')'     Yqc = Qcent(2,iPassQ) + Cy'
          WRITE(I,'(A)')'     Zqc = Qcent(3,iPassQ) + Cz'       
          do J=13,16
             I=non1Prim(J)
             WRITE(I,'(A)')'     Xqd = Qcent(1,iPrimQ,iPassQ) + Dx'
             WRITE(I,'(A)')'     Yqd = Qcent(2,iPrimQ,iPassQ) + Dy'
             WRITE(I,'(A)')'     Zqd = Qcent(3,iPrimQ,iPassQ) + Dz'
          enddo
          I=pure1Prim(4)
          WRITE(I,'(A)')'     Xqd = Qcent(1,iPassQ) + Dx'
          WRITE(I,'(A)')'     Yqd = Qcent(2,iPassQ) + Dy'
          WRITE(I,'(A)')'     Zqd = Qcent(3,iPassQ) + Dz'
          do J=1,8
             I=non1Prim(J)
             WRITE(I,'(A)')'     alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP'    
          enddo
          WRITE(LUMODA5,'(A)')'     alphaP = -reducedExponents(1)*invexpP'    
          WRITE(LUMODB5,'(A)')'     alphaP = -reducedExponents(1)*invexpP'    
          do J=9,16
             I=non1Prim(J)
             WRITE(I,'(A)')'     alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ(iPrimQ)'
          enddo
          WRITE(LUMODC5,'(A)')'     alphaQ = -reducedExponents(1)*invexpQ'
          WRITE(LUMODD5,'(A)')'     alphaQ = -reducedExponents(1)*invexpQ'
          do I=1,10
             WRITE(I,'(A)')'     alphaXpq = -alphaP*Xpq'
             WRITE(I,'(A)')'     alphaYpq = -alphaP*Ypq'
             WRITE(I,'(A)')'     alphaZpq = -alphaP*Zpq'
          enddo
          do J=9,16
             I=non1Prim(J)
             WRITE(I,'(A)')'     inv2expQ = D05*invexpQ(iPrimQ)'
          enddo
          do J=3,4
             I=pure1Prim(J)
             WRITE(I,'(A)')'     inv2expQ = D05*invexpQ'
          enddo
          do I=11,20
             WRITE(I,'(A)')'     alphaXpq = alphaQ*Xpq'
             WRITE(I,'(A)')'     alphaYpq = alphaQ*Ypq'
             WRITE(I,'(A)')'     alphaZpq = alphaQ*Zpq'
          enddo
          do I=1,200
             WRITE(I,'(A)')'     squaredDistance = Xpq*Xpq+Ypq*Ypq+Zpq*Zpq'
          enddo
          do J=1,16
             I=non1Prim(J)
             WRITE(I,'(A)')'     WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance'
          enddo
          do J=1,4
             I=pure1Prim(J)
             WRITE(I,'(A)')'     WVAL = reducedExponents(1)*squaredDistance'
          enddo
          do I=1,20
             !          WRITE(I,'(A)')'     !  0 < WVAL < 0.000001'
             !          WRITE(I,'(A)')'!     IF (ABS(WVAL) .LT. SMALL) THEN'     
             !          WRITE(I,'(A)')'!      RJ000(0) = D1 !THE BOYS FUNCTION FOR ZERO ARGUMENT'
             !          DO J=1,JMAX
             !             WRITE(I,'(A,I1,A,ES23.16,A)')'!      RJ000(',J,')= ',1.0d0/(2*J + 1),'_realk'
             !          ENDDO
             WRITE(I,'(A)')'     !  0 < WVAL < 12 '
             WRITE(I,'(A)')'     IF (WVAL .LT. D12) THEN'
             WRITE(I,'(A)')'      IPNT = NINT(D100*WVAL)'
             WRITE(I,'(A)')'      WDIFF = WVAL - TENTH*IPNT'
             WRITE(I,'(A)')'      W2    = WDIFF*WDIFF'
             WRITE(I,'(A)')'      W3    = W2*WDIFF'
             WRITE(I,'(A)')'      W2    = W2*D05'
             WRITE(I,'(A)')'      W3    = W3*COEF3'
             DO J=0,JMAX
                WRITE(I,'(A,I2,A,I2,A,I2,A,I2,A,I2,A)')'      RJ000(',J,') = TABFJW(',J,',IPNT)-TABFJW(',J+1,',IPNT)*WDIFF+TABFJW(',J+2,',IPNT)*W2+TABFJW(',J+3,',IPNT)*W3'
             ENDDO
             !     WRITE(I,'(A)')'        DO J=0,JMAX'
             !      WRITE(I,'(A)')'           R = TABFJW(J,IPNT)'
             !      WRITE(I,'(A)')'           R = R -TABFJW(J+1,IPNT)*WDIFF'
             !      WRITE(I,'(A)')'           R = R + TABFJW(J+2,IPNT)*W2'
             !      WRITE(I,'(A)')'           R = R + TABFJW(J+3,IPNT)*W3'
             !      WRITE(I,'(A)')'           RJ000(J,ipq,ipassq) = R'
             !     WRITE(I,'(A)')'        ENDDO'
             WRITE(I,'(A)')'     !  12 < WVAL <= (2J+36) '
             WRITE(I,'(A)')'     ELSE IF (WVAL.LE.D2JP36) THEN'
             WRITE(I,'(A)')'      REXPW = D05*EXP(-WVAL)'
             WRITE(I,'(A)')'      RWVAL = D1/WVAL'
             WRITE(I,'(A)')'      GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))'
             WRITE(I,'(A)')'      RJ000(0) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
             !        WRITE(I,'(A)')'        RJ000(0,ipq,ipassq) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
             DO J=1,JMAX
                WRITE(I,'(A,I2,A,I2,A,I2,A)')'      RJ000(',J,') = RWVAL*((',J,' - D05)*RJ000(',J-1,')-REXPW)'          
             ENDDO
             !        WRITE(I,'(A)')'        DO J=1,JMAX'
             !        WRITE(I,'(A)')'           RJ000(J,ipq,ipassq) = RWVAL*((J - D05)*RJ000(J-1,ipq,ipassq)-REXPW)'
             !        WRITE(I,'(A)')'        ENDDO'
             WRITE(I,'(A)')'     !  (2J+36) < WVAL '
             WRITE(I,'(A)')'     ELSE'
             WRITE(I,'(A)')'      RWVAL = PID4/WVAL'
             WRITE(I,'(A)')'      RJ000(0) = SQRT(RWVAL)'
             !        WRITE(I,'(A)')'        RJ000(0,ipq,ipassq) = SQRT(RWVAL)'
             WRITE(I,'(A)')'      RWVAL = RWVAL*PID4I'
             DO J=1,JMAX
                WRITE(I,'(A,I2,A,I2,A,I2,A)')'      RJ000(',J,') = RWVAL*(',J,' - D05)*RJ000(',J-1,')'
             ENDDO
             !        WRITE(I,'(A)')'        DO J = 1, JMAX'
             !        WRITE(I,'(A)')'           RJ000(J,ipq,ipassq) = RWVAL*(J - D05)*RJ000(J-1,ipq,ipassq)'
             !        WRITE(I,'(A)')'        ENDDO'
             WRITE(I,'(A)')'     ENDIF'
          enddo
          do J=1,16
             I=non1Prim(J)
             WRITE(I,'(A)')'     PREF = integralPrefactor(iPrimQ,iPrimP)*QpreExpFac(iPrimQ,iPassQ)*Pexpfac'
          enddo
          do J=1,4
             I=pure1Prim(J)
             WRITE(I,'(A)')'     PREF = integralPrefactor(1)*QpreExpFac(iPassQ)*PpreExpFac(1)'
          enddo
          !gen og seg
          DO I=1,16,5
             WRITE(I,'(A)')'     Auxarray(1,IP) = PREF*RJ000(0)'
             WRITE(I+1,'(A)')'     TMPAuxarray(1) = PREF*RJ000(0)'
             WRITE(I+2,'(A)')'     TMPAuxarray(1) = PREF*RJ000(0)'
             WRITE(I+3,'(A)')'     TMPAuxarray(1) = PREF*RJ000(0)'
             WRITE(I+4,'(A)')'     TMPAuxarray(1) = PREF*RJ000(0)'
          ENDDO
          !seg TMPAuxarray
          do I=1,20
             DO J=2,JMAX+1
                WRITE(I,'(A,I2,A,I2,A)')'     TMParray1(1,',J,') = PREF*RJ000(',J-1,')'
             ENDDO
          enddo
          do I=1,20
             Gen = .FALSE.; SegQ = .FALSE.; SegP = .FALSE.; Seg = .FALSE.; Seg1Prim = .FALSE.
             IF(MOD(I,5).EQ.1)Gen = .TRUE.
             IF(MOD(I,5).EQ.2)SegQ = .TRUE.
             IF(MOD(I,5).EQ.3)SegP = .TRUE.
             IF(MOD(I,5).EQ.4)Seg = .TRUE.
             IF(MOD(I,5).EQ.0)THEN
                Seg = .TRUE.
                Seg1Prim=.TRUE.
             ENDIF
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

                   call WriteTwoTerms(J,JTMP,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,TwoTermsUsed,I,Gen,SegQ,SegP,Seg,seg1prim)


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
                               call TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,I,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
                            ELSEIF(UREC)THEN
                               !                      print*,'!A URECURRENCE iTUV',iTUV
                               call URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,I,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
                            ELSEIF(VREC)THEN
                               !                      print*,'!A VRECURRENCE iTUV',iTUV
                               call VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,I,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
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
                                  CALL TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,I,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
                               ELSEIF(.NOT.(UREC.AND.UREC2).AND.UREC)THEN
                                  !                         print*,'!B URECURRENCE iTUV',iTUV
                                  CALL URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,I,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
                               ELSEIF(.NOT.(VREC.AND.VREC2).AND.VREC)THEN
                                  !                         print*,'!B VRECURRENCE iTUV',iTUV
                                  CALL VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,I,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
                               ENDIF
                            ELSE
                               !chose one of the possibilities
                               IF(TREC)THEN
                                  !                         print*,'!C TRECURRENCE iTUV',iTUV
                                  CALL TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,I,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
                               ELSEIF(UREC)THEN
                                  !                         print*,'!C URECURRENCE iTUV',iTUV
                                  call URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,I,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
                               ELSEIF(VREC)THEN
                                  !                         print*,'!D VRECURRENCE iTUV',iTUV
                                  call VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,I,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
                               ENDIF
                            ENDIF
                         ENDIF
                         CREATED(Tp,Up,Vp) = .TRUE.
                      ENDDO
                   ENDDO
                ENDDO
                deallocate(TwoTermTUVLIST)
                deallocate(TwoTermsUsed)
                ! WRITE(I,'(A,I4)')'           ENDDO'
             ENDDO
          enddo
          do J=1,16
             I=non1Prim(J)       
             WRITE(I,'(A)')'    ENDDO'
             WRITE(I,'(A)')'   ENDDO'
          enddo
          do I=1,20
             WRITE(I,'(A)')'  ENDDO'
             WRITE(I,'(A)')' end subroutine'
          enddo
          deallocate(TUVINDEX)
          deallocate(TINDEX)
          deallocate(UINDEX)
          deallocate(VINDEX)
          deallocate(JINDEX)
       ENDDO
       do I=1,20
          WRITE(I,'(A)')'end module'
       enddo


       close(unit = LUMODA1)
       close(unit = LUMODA2)
       close(unit = LUMODA3)
       close(unit = LUMODA4)
       close(unit = LUMODA5)

       close(unit = LUMODB1)
       close(unit = LUMODB2)
       close(unit = LUMODB3)
       close(unit = LUMODB4)
       close(unit = LUMODB5)

       close(unit = LUMODC1)
       close(unit = LUMODC2)
       close(unit = LUMODC3)
       close(unit = LUMODC4)
       close(unit = LUMODC5)

       close(unit = LUMODD1)
       close(unit = LUMODD2)
       close(unit = LUMODD3)
       close(unit = LUMODD4)
       close(unit = LUMODD5)
    END DO

  END subroutine PASSsub
  
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
     & TwoTermTUVLIST,JTMP,TwoTermsUsed,LUPRI,Gen,SegQ,SegP,Seg,seg1prim)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,JTMP,I,J1
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),Gen,SegQ,SegP,Seg,seg1prim
integer :: nTUVLIST,nTUVLISTactual,LUPRI
integer :: TwoTermTUVLIST(nTUVLIST)
logical :: Unique,TREC,UREC,VREC,TwoTermsUsed(nTUVLIST)

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
         IF(lupri.LT.11)THEN
            call AddToString(') = inv2expP*(')
         ELSEIF(lupri.GT.10)THEN
            call AddToString(') = inv2expQ*(')
         ENDIF
         IF(Gen)THEN
            call AddToString('AuxArray(')
         ELSE
            call AddToString('TMPAuxArray(')
         ENDIF
         call AddToString(TwoTermTUVLIST(ituvP2))
         IF(Gen)THEN
            call AddToString(',IP)')
         ELSE
            call AddToString(')')
         ENDIF
         call AddToString(' + ')
         IF(lupri.LT.11)THEN
            call AddToString('alphaP*TmpArray')
         ELSEIF(lupri.GT.10)THEN
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
         IF(lupri.LT.11)THEN
            call AddToString(') = inv2expP*(')
         ELSEIF(lupri.GT.10)THEN
            call AddToString(') = inv2expQ*(')
         ENDIF
         call AddToString('TmpArray')
         call AddToString(JTMP-1)
         call AddToString('(')
         call AddToString(TwoTermTUVLIST(ituvP2))
         call AddToString(',')
         call AddToString(J1)
         IF(lupri.LT.11)THEN
            call AddToString(') + alphaP*TmpArray')
         ELSEIF(lupri.GT.10)THEN
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

subroutine TRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms,lupri,nTUVprev
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),Gen,SegQ,SegP,Seg,seg1prim
integer :: nTUVLIST,nTUVLISTactual,JTMP,TM2
integer :: TwoTermTUVLIST(nTUVLIST)
character(len=3) :: DIRECTIONSTRING
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
IF(lupri.LT.6)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Xpa','alphaXpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
ELSEIF(lupri.LT.11)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Xpb','alphaXpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
ELSEIF(lupri.LT.16)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Xqc','alphaXpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
ELSEIF(lupri.LT.21)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Xqd','alphaXpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
ENDIF
end subroutine TRECURRENCE


subroutine URECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms,lupri,nTUVprev
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),Gen,SegQ,SegP,Seg,seg1prim
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
IF(lupri.LT.6)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Ypa','alphaYpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
ELSEIF(lupri.LT.11)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Ypb','alphaYpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
ELSEIF(lupri.LT.16)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Yqc','alphaYpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
ELSEIF(lupri.LT.21)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Yqd','alphaYpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,seg1prim,nTUVprev)
ENDIF
end subroutine URECURRENCE

subroutine VRECURRENCE(Tp,Up,Vp,J,TUVINDEX,CREATED,JMAX,nTUVLIST,nTUVLISTactual,TwoTermTUVLIST,JTMP,lupri,Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev)
implicit none
integer :: Tp,Up,Vp,J,ituvP,TM1,ituvP2,ituvP3,JMAX,ituvp0,ituvp1,I,iTwoTerms,lupri,nTUVprev
integer :: TUVINDEX(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1)
logical :: CREATED(-2:JMAX+1,-2:JMAX+1,-2:JMAX+1),Gen,SegQ,SegP,Seg,Seg1Prim
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
IF(lupri.LT.6)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Zpa','alphaZpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev)
ELSEIF(lupri.LT.11)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Zpb','alphaZpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev)
ELSEIF(lupri.LT.16)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Zqc','alphaZpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev)
ELSEIF(lupri.LT.21)THEN
   CALL XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
        & 'Zqd','alphaZpq',TM1,TM2,lupri,Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev)
ENDIF

end subroutine VRECURRENCE

SUBROUTINE XYZVERTICALRECURRENCE(J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,&
     & DIRECTIONSTRING1,DIRECTIONSTRING2,TM1,TM2,lupri,&
     & Gen,SegQ,SegP,Seg,Seg1prim,nTUVprev)
implicit none
integer :: J,JMAX,ituvP0,ituvP1,JTMP,iTwoTerms,TM1,TM2,lupri,nTUVprev
logical :: Gen,SegQ,SegP,Seg,Seg1Prim
character(len=3) :: DIRECTIONSTRING1 !Xpa
character(len=8) :: DIRECTIONSTRING2 !alphaXpq
   
call initString(5)
IF(J.EQ.1)THEN
   IF(Gen)THEN
      !place in AuxArray and use AuxArray in X*Aux
      call AddToString('AuxArray(')
      call AddToString(ituvP0)
      call AddToString(',IP) = ')
      !      call AddToString('Xpa')
      call AddToString(DIRECTIONSTRING1) !Xpa
      call AddToString('*AuxArray(')
      call AddToString(ituvP1)
      call AddToString(',IP)')
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
               WRITE(lupri,'(A)')'      AuxArray(iTUV,IPassQ) = TMPAuxarray(iTUV)'
            ELSEIF(Seg)THEN
               WRITE(lupri,'(A)')'      AuxArray(iTUV,IPassQ) = AuxArray(iTUV,IPassQ) + TMPAuxarray(iTUV)'
            ELSE !segQ or SegP
               WRITE(lupri,'(A)')'      AuxArray(iTUV,IP) = AuxArray(iTUV,IP) + TMPAuxarray(iTUV)'
            ENDIF
            WRITE(lupri,'(A)')'     enddo'
         ENDIF
         !place in AuxArray and use AuxArray in X*Aux
         call AddToString('AuxArray(')
         call AddToString(ituvP0)
         IF(Seg)THEN
            call AddToString(',IPassQ)')
         ELSE !segQ or SegP or Gen
            call AddToString(',IP)')
         ENDIF
         call AddToString(' = ')
         IF(.NOT.Seg1Prim)THEN
            call AddToString('AuxArray(')
            call AddToString(ituvP0)
            IF(Seg)THEN
               call AddToString(',IPassQ)')
            ELSE !segQ or SegP or Gen
               call AddToString(',IP)')
            ENDIF
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
