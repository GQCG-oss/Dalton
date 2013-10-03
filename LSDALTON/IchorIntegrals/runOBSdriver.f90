PROGRAM TUV
  use math
  implicit none
  integer,pointer :: TUVINDEX(:,:,:),TUVINDEXP(:,:,:)
  integer :: JMAX,J,JMAX1,JMAXP
  logical,pointer :: Enoscreen(:,:),EnoscreenS(:,:),zero(:)
  integer :: ijk1,ijk2,ijkcart,ijk,ijkcart1,ijkcart2,nTUV,ijkP
  integer :: iTUV,ilmP
  real(realk),pointer :: SCMAT1(:,:),SCMAT2(:,:)
  logical :: sph1,sph2,sphericalTrans,sphericalGTO,newparam,ELSESTATEMENT,Spherical
  integer :: LUMAIN,LUMOD1,LUMOD2,LUMOD3,LUMOD4,iparam,iparam2,nparam
  integer :: JMAX2,JMAX3,JMAX4,ituvP,jp,tp,up,vp,ntuvP,l1,l2,l12
  integer :: ip1,jp1,kp1,p1,ip2,jp2,kp2,p2,ijkcartP,AngmomA,AngmomB
  integer :: AngmomC,AngmomD,AngmomP,AngmomQ,nTUVQ,AngmomPQ
  integer :: nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec
  integer :: nTUVA,nTUVB,nTUVC,nTUVD
  integer :: nlmA,nlmB,nlmC,nlmD,angmomID
  real(realk),pointer :: uniqeparam(:)
  character(len=15),pointer :: uniqeparamNAME(:)
  character(len=9) :: STRINGIN,STRINGOUT,TMPSTRING
  character(len=4) :: SPEC
!TODO
!remove mem_alloc
!remove CDAB = TMParray2
!add PrimitiveContractionSeg to Transfer or Vertical 
!

  LUMOD3=3
  open(unit = LUMOD3, file="MAIN_OBS_DRIVER.f90",status="unknown")
  WRITE(LUMOD3,'(A)')'MODULE IchorEriCoulombintegralOBSGeneralMod'
  WRITE(LUMOD3,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD3,'(A)')'use IchorprecisionModule'
  WRITE(LUMOD3,'(A)')'use IchorCommonModule'
  WRITE(LUMOD3,'(A)')'use IchorMemory'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODA'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODD'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODC'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODAtoC'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODDtoA'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODCtoA'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_HorizontalRecurrenceLHSModAtoB'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_HorizontalRecurrenceRHSModCtoD'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_HorizontalRecurrenceRHSModDtoC'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_Sphcontract1Mod'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_Sphcontract2Mod'
  WRITE(LUMOD3,'(A)')'  '
  WRITE(LUMOD3,'(A)')'private   '
  WRITE(LUMOD3,'(A)')'public :: IchorCoulombIntegral_OBS_general,IchorCoulombIntegral_OBS_general_size  '
  WRITE(LUMOD3,'(A)')'  '
  WRITE(LUMOD3,'(A)')'CONTAINS'
  WRITE(LUMOD3,'(A)')'  '
  WRITE(LUMOD3,'(A)')'  '
  WRITE(LUMOD3,'(A)')'  subroutine IchorCoulombIntegral_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD3,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD3,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD3,'(A)')'       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
  WRITE(LUMOD3,'(A)')'       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&'
  WRITE(LUMOD3,'(A)')'       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&'
  WRITE(LUMOD3,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&'
  WRITE(LUMOD3,'(A)')'       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsC,nAtomsD,spherical,&'
  WRITE(LUMOD3,'(A)')'       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)'
  WRITE(LUMOD3,'(A)')'    implicit none'
  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD'
  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri'
  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2'
  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nAtomsC,nAtomsD'
  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)'
  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)'
  WRITE(LUMOD3,'(A)')'    logical,intent(in)     :: Qsegmented,Psegmented'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: qcent(3*nPrimQ*MaxPasses) !qcent(3,nPrimQ,MaxPasses)'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: QpreExpFac(nPrimQ*MaxPasses),PpreExpFac(nPrimP)'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)'
  WRITE(LUMOD3,'(A)')'    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
  WRITE(LUMOD3,'(A)')'    !    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)'
  WRITE(LUMOD3,'(A)')'    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
  WRITE(LUMOD3,'(A)')'    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: CDAB(:)'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: integralPrefactor(nPrimQP)'
  WRITE(LUMOD3,'(A)')'    logical,intent(in) :: PQorder'
  WRITE(LUMOD3,'(A)')'    !integralPrefactor(nPrimP,nPrimQ)'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: reducedExponents(nPrimQP)'
  WRITE(LUMOD3,'(A)')'    !reducedExponents(nPrimP,nPrimQ)'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: Qdistance12(3*MaxPasses) !Ccenter-Dcenter'
  WRITE(LUMOD3,'(A)')'    !Qdistance12(3,MaxPasses)'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter '
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: Acenter(3),Bcenter(3),Ccenter(3,nAtomsC),Dcenter(3,nAtomsD)'
  WRITE(LUMOD3,'(A)')'    logical,intent(in) :: spherical'
  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize'
  WRITE(LUMOD3,'(A)')'!   TMP variables - allocated outside'  
  WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)'

  WRITE(LUMOD3,'(A)')'!   Local variables '  
  WRITE(LUMOD3,'(A)')'    real(realk),pointer :: squaredDistance(:)'!,Rpq(:)'!,Rqc(:),Rpa(:)
  WRITE(LUMOD3,'(A)')'    integer :: AngmomPQ,AngmomP,AngmomQ,I,J,nContQP,la,lb,lc,ld,nsize,angmomid'
  WRITE(LUMOD3,'(A)')'    real(realk),pointer :: RJ000(:),OUTPUTinterest(:)'


  WRITE(LUMOD3,'(A)')'  '
  WRITE(LUMOD3,'(A)')' IF(nAtomsC*nAtomsD.NE.nPasses)Call ichorquit(''nPass error'',-1)'
  WRITE(LUMOD3,'(A)')'  '
  WRITE(LUMOD3,'(A)')'!IF(.TRUE.)THEN'
  WRITE(LUMOD3,'(A)')'!    call interest_initialize()'
  WRITE(LUMOD3,'(A)')'!'
  WRITE(LUMOD3,'(A)')'!    nsize = size(CDAB)*10'
  WRITE(LUMOD3,'(A)')'!    call mem_ichor_alloc(OUTPUTinterest,nsize)'
  WRITE(LUMOD3,'(A)')'!    nsize = nPrimQP'
  WRITE(LUMOD3,'(A)')'!    '
  WRITE(LUMOD3,'(A)')'!    '
  WRITE(LUMOD3,'(A)')'!    la = AngmomA+1'
  WRITE(LUMOD3,'(A)')'!    lb = AngmomB+1'
  WRITE(LUMOD3,'(A)')'!    lc = AngmomC+1'
  WRITE(LUMOD3,'(A)')'!    ld = AngmomD+1'
  WRITE(LUMOD3,'(A)')'!    call interest_eri(OUTPUTinterest,nsize,&'
  WRITE(LUMOD3,'(A)')'!       & la,Aexp,Acenter(1),Acenter(2),Acenter(3),ACC,&'
  WRITE(LUMOD3,'(A)')'!         & lb,Bexp,Bcenter(1),Bcenter(2),Bcenter(3),BCC,&'
  WRITE(LUMOD3,'(A)')'!         & lc,Cexp,Ccenter(1,1),Ccenter(2,1),Ccenter(3,1),CCC,&'
  WRITE(LUMOD3,'(A)')'!         & ld,Dexp,Dcenter(1,1),Dcenter(2,1),Dcenter(3,1),DCC,&'
  WRITE(LUMOD3,'(A)')'!         & lupri)!,&'
  WRITE(LUMOD3,'(A)')'!         !         & .false.)'
  WRITE(LUMOD3,'(A)')'!    write(lupri,*)''OUTPUTinterest'',OUTPUTinterest'
  WRITE(LUMOD3,'(A)')'!ENDIF'
  WRITE(LUMOD3,'(A)')' '
!  WRITE(LUMOD3,'(A)')'    !build the distance between center P and A in Rpa '
!  WRITE(LUMOD3,'(A)')'    !used in Vertical and Electron Transfer Recurrence Relations '
!  WRITE(LUMOD3,'(A)')'    call mem_ichor_alloc(Rpa,3*nPrimP)'
!  WRITE(LUMOD3,'(A)')'    call build_Rpa(nPrimP,Pcent,Acenter,Rpa)'
!  WRITE(LUMOD3,'(A)')'    !build the distance between center Q and C in Rqc used'
!  WRITE(LUMOD3,'(A)')'    !used in Electron Transfer Recurrence Relations '
!  WRITE(LUMOD3,'(A)')'    call mem_ichor_alloc(Rqc,3*nPrimQ)'
!  WRITE(LUMOD3,'(A)')'    call build_Rpa(nPrimQ,Qcent,Ccenter,Rqc)'
  WRITE(LUMOD3,'(A)')'    '
  WRITE(LUMOD3,'(A)')'    IF(PQorder)THEN'
  WRITE(LUMOD3,'(A)')'       call IchorQuit(''PQorder OBS general expect to get QP ordering'',-1)'
  WRITE(LUMOD3,'(A)')'    ENDIF'
  WRITE(LUMOD3,'(A)')'    IF(.NOT.spherical)THEN'
  WRITE(LUMOD3,'(A)')'       call IchorQuit(''cartesian not testet'',-1)'
  WRITE(LUMOD3,'(A)')'    ENDIF'
  WRITE(LUMOD3,'(A)')'    '
!  WRITE(LUMOD3,'(A)')'    call mem_ichor_alloc(squaredDistance,nPasses*nPrimQP)'
!  WRITE(LUMOD3,'(A)')'    call mem_ichor_alloc(Rpq,3*nPrimQP*nPasses)'
!  WRITE(LUMOD3,'(A)')'    !builds squaredDistance between center P and Q. Order(nPrimQ,nPrimP,nPassQ)'
!  WRITE(LUMOD3,'(A)')'    !builds Distance between center P and Q in Rpq. Order(3,nPrimQ,nPrimP,nPassQ)'
!  WRITE(LUMOD3,'(A)')'    call build_QP_squaredDistance_and_Rpq(nPrimQ,nPasses,nPrimP,Qcent,Pcent,&'
!  WRITE(LUMOD3,'(A)')'         & squaredDistance,Rpq)'
  WRITE(LUMOD3,'(A)')'    '
  WRITE(LUMOD3,'(A)')'    !Setup combined Angmom info'
  WRITE(LUMOD3,'(A)')'    AngmomP = AngmomA+AngmomB'
  WRITE(LUMOD3,'(A)')'    AngmomQ = AngmomC+AngmomD'
  WRITE(LUMOD3,'(A)')'    AngmomPQ  = AngmomP + AngmomQ'
!  WRITE(LUMOD3,'(A)')'    !Build the Boys Functions for argument squaredDistance*reducedExponents'
!  WRITE(LUMOD3,'(A)')'    !save in RJ000 ordering (AngmomPQ+1),nPrimQ,nPrimP,nPasses'
!  WRITE(LUMOD3,'(A)')'    call mem_ichor_alloc(RJ000,(AngmomPQ+1)*nPasses*nPrimQP)'
!  WRITE(LUMOD3,'(A)')'    call buildRJ000_general(nPasses,nPrimQ,nPrimP,nTABFJW1,nTABFJW2,reducedExponents,&'
!  WRITE(LUMOD3,'(A)')'         & TABFJW,RJ000,AngmomPQ,Pcent,Qcent)'
!  WRITE(LUMOD3,'(A)')'    IF (INTPRINT .GE. 10) THEN'
!  WRITE(LUMOD3,'(A)')'     WRITE(lupri,*)''Output from W000'''
!  WRITE(LUMOD3,'(A)')'     DO I=1,nPrimQ*nPrimP*nPasses'
!  WRITE(LUMOD3,'(A)')'       DO J=0,AngmomPQ'
!  WRITE(LUMOD3,'(A)')'          WRITE(LUPRI,''(2X,A6,I4,A1,I4,A2,ES16.8)'')''RJ000('',J,'','',I,'')='',RJ000(1+J+(I-1)*(AngmomPQ+1))'
!  WRITE(LUMOD3,'(A)')'       ENDDO'
!  WRITE(LUMOD3,'(A)')'     ENDDO'
!  WRITE(LUMOD3,'(A)')'    END IF'
  WRITE(LUMOD3,'(A)')'!    nTUV = (AngmomPQ+1)*(AngmomPQ+2)*(AngmomPQ+3)/6'
  WRITE(LUMOD3,'(A)')'!    nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6'
  WRITE(LUMOD3,'(A)')'!    nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6'
  WRITE(LUMOD3,'(A)')'!    nTUVC = (AngmomC+1)*(AngmomC+2)*(AngmomC+3)/6'
  WRITE(LUMOD3,'(A)')'!    nTUVD = (AngmomD+1)*(AngmomD+2)*(AngmomD+3)/6'
  WRITE(LUMOD3,'(A)')'!    nlmA = 2*AngmomA+1'
  WRITE(LUMOD3,'(A)')'!    nlmB = 2*AngmomB+1'
  WRITE(LUMOD3,'(A)')'!    nlmC = 2*AngmomC+1'
  WRITE(LUMOD3,'(A)')'!    nlmD = 2*AngmomD+1'
  WRITE(LUMOD3,'(A)')'    AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD'
       WRITE(LUMOD3,'(A)')'    SELECT CASE(AngmomID)'
  
  DO AngmomA = 0,2
!   DO AngmomB = 0,2!AngmomA
   DO AngmomB = 0,AngmomA
!    DO AngmomC = 0,2!AngmomA
    DO AngmomC = 0,AngmomA
!     DO AngmomD = 0,2!AngmomC
     DO AngmomD = 0,AngmomC
        AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD
        WRITE(LUMOD3,'(A,I4,A,I2,A,I2,A,I2,A,I2,A)')'    CASE(',AngmomID,')  !Angmom(A=',AngmomA,',B=',AngmomB,',C=',AngmomC,',D=',AngmomD,') combi'
      AngmomP = AngmomA + AngmomB
      AngmomQ = AngmomC + AngmomD
!      IF(AngmomQ.GT.AngmomP)CYCLE
      AngmomPQ = AngmomA + AngmomB + AngmomC + AngmomD
      nTUV = (AngmomPQ+1)*(AngmomPQ+2)*(AngmomPQ+3)/6
      nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
      nTUVQ = (AngmomQ+1)*(AngmomQ+2)*(AngmomQ+3)/6
      nTUVAspec = (AngmomA+1)*(AngmomA+2)/2
      nTUVBspec = (AngmomB+1)*(AngmomB+2)/2
      nTUVCspec = (AngmomC+1)*(AngmomC+2)/2
      nTUVDspec = (AngmomD+1)*(AngmomD+2)/2
      nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
      nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6
      nTUVC = (AngmomC+1)*(AngmomC+2)*(AngmomC+3)/6
      nTUVD = (AngmomD+1)*(AngmomD+2)*(AngmomD+3)/6
      IF((AngmomA.GT.1.OR.AngmomB.GT.1).OR.(AngmomC.GT.1.OR.AngmomD.GT.1))THEN
         spherical = .TRUE.
         nlmA = 2*AngmomA+1
         nlmB = 2*AngmomB+1
         nlmC = 2*AngmomC+1
         nlmD = 2*AngmomD+1
         STRINGIN(1:9)  = 'TMParray1'
         STRINGOUT(1:9) = 'TMParray2'
         TMPSTRING(1:9) = '         '
!         WRITE(LUMOD3,'(A)')'      IF(spherical)THEN'
         call subroutineMAIN(LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
              & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,spherical)
!         WRITE(LUMOD3,'(A)')'      ELSE'
!         spherical = .FALSE.
!         nlmA = nTUVAspec
!         nlmB = nTUVBspec
!         nlmC = nTUVCspec
!         nlmD = nTUVDspec
!         STRINGIN(1:9)  = 'TMParray1'
!         STRINGOUT(1:9) = 'TMParray2'
!         TMPSTRING(1:9) = '         '
!         call subroutineMAIN(LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
!              & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,spherical)
!         WRITE(LUMOD3,'(A)')'      ENDIF'
      ELSE
         spherical = .TRUE.
         nlmA = 2*AngmomA+1
         nlmB = 2*AngmomB+1
         nlmC = 2*AngmomC+1
         nlmD = 2*AngmomD+1
         STRINGIN(1:9)  = 'TMParray1'
         STRINGOUT(1:9) = 'TMParray2'
         TMPSTRING(1:9) = '         '
         call subroutineMAIN(LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
              & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,spherical)
      ENDIF
     ENDDO
!     WRITE(LUMOD3,'(A)')'       ENDIF ! D if statement'
    ENDDO
!    WRITE(LUMOD3,'(A)')'      ENDIF ! C if statement'
   ENDDO
!   WRITE(LUMOD3,'(A)')'     ENDIF ! B if statement'
  ENDDO
!  WRITE(LUMOD3,'(A)')'    ENDIF ! A if statement'
  WRITE(LUMOD3,'(A)')'    CASE DEFAULT'
  WRITE(LUMOD3,'(A)')'        CALL ICHORQUIT(''Unknown Case in IchorCoulombIntegral_OBS_general'',-1)'
  WRITE(LUMOD3,'(A)')'    END SELECT'
  
  WRITE(LUMOD3,'(A)')'  end subroutine IchorCoulombIntegral_OBS_general'

  WRITE(LUMOD3,'(A)')'  '
  WRITE(LUMOD3,'(A)')'  subroutine IchorCoulombIntegral_OBS_general_size(TMParray1maxsize,&'
  WRITE(LUMOD3,'(A)')'         &TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD3,'(A)')'         &nPrimQP,nContQP)'
  WRITE(LUMOD3,'(A)')'    implicit none'
  WRITE(LUMOD3,'(A)')'    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize'
  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD'
  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimQP,nContQP'
  WRITE(LUMOD3,'(A)')'    ! local variables'
  WRITE(LUMOD3,'(A)')'    integer :: AngmomID'
  WRITE(LUMOD3,'(A)')'    '
  WRITE(LUMOD3,'(A)')'    AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD'
  WRITE(LUMOD3,'(A)')'    TMParray2maxSize = 0'
  WRITE(LUMOD3,'(A)')'    TMParray1maxSize = 0'
  WRITE(LUMOD3,'(A)')'    SELECT CASE(AngmomID)'  
  DO AngmomA = 0,2
   DO AngmomB = 0,AngmomA
    DO AngmomC = 0,AngmomA
     DO AngmomD = 0,AngmomC
      AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD
      WRITE(LUMOD3,'(A,I4,A,I2,A,I2,A,I2,A,I2,A)')'    CASE(',AngmomID,')  !Angmom(A=',AngmomA,',B=',AngmomB,',C=',AngmomC,',D=',AngmomD,') combi'
      AngmomP = AngmomA + AngmomB
      AngmomQ = AngmomC + AngmomD
      AngmomPQ = AngmomA + AngmomB + AngmomC + AngmomD
      nTUV = (AngmomPQ+1)*(AngmomPQ+2)*(AngmomPQ+3)/6
      nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
      nTUVQ = (AngmomQ+1)*(AngmomQ+2)*(AngmomQ+3)/6
      nTUVAspec = (AngmomA+1)*(AngmomA+2)/2
      nTUVBspec = (AngmomB+1)*(AngmomB+2)/2
      nTUVCspec = (AngmomC+1)*(AngmomC+2)/2
      nTUVDspec = (AngmomD+1)*(AngmomD+2)/2
      nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
      nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6
      nTUVC = (AngmomC+1)*(AngmomC+2)*(AngmomC+3)/6
      nTUVD = (AngmomD+1)*(AngmomD+2)*(AngmomD+3)/6

      spherical = .TRUE.
      nlmA = 2*AngmomA+1
      nlmB = 2*AngmomB+1
      nlmC = 2*AngmomC+1
      nlmD = 2*AngmomD+1
      STRINGIN(1:9)  = 'TMParray1'
      STRINGOUT(1:9) = 'TMParray2'
      TMPSTRING(1:9) = '         '
      call determineSizes(LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
           & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,spherical)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
  WRITE(LUMOD3,'(A)')'    CASE DEFAULT'
  WRITE(LUMOD3,'(A)')'        CALL ICHORQUIT(''Unknown Case in IchorCoulombIntegral_OBS_general_size'',-1)'
  WRITE(LUMOD3,'(A)')'    END SELECT'
  WRITE(LUMOD3,'(A)')'  end subroutine IchorCoulombIntegral_OBS_general_size'


WRITE(LUMOD3,'(A)')'  subroutine PrimitiveContraction(AUXarray2,AUXarrayCont,nTUVfull,nPrimP,nPrimQ,nPasses,&'
WRITE(LUMOD3,'(A)')'       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD)'
WRITE(LUMOD3,'(A)')'    implicit none'
WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nTUVfull,nPrimP,nPrimQ,nPasses,nContP,nContQ'
WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)'
WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray2(nTUVfull,nPrimQ,nPrimP,nPasses)'
WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nTUVfull,nContQ,nContP,nPasses)'
WRITE(LUMOD3,'(A)')'    !'
WRITE(LUMOD3,'(A)')'    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD'
WRITE(LUMOD3,'(A)')'    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP'
WRITE(LUMOD3,'(A)')'    real(realk) :: B,ABCDTMP,ABDTMP,ABTMP'
WRITE(LUMOD3,'(A)')'    !all passes have same ACCs,BCCs,...'
WRITE(LUMOD3,'(A)')'    !maybe construct big CC(nPrimQP,nContQP) matrix and call dgemm nPass times'
WRITE(LUMOD3,'(A)')'    !the construction of CC should scale as c**4*p**4 and the '
WRITE(LUMOD3,'(A)')'    !dgemm should scale as c**4*p**4*L**6 but hopefully with efficient FLOP count, although not quadratic matrices....'
WRITE(LUMOD3,'(A)')'    !special for nContPQ = 1 '
WRITE(LUMOD3,'(A)')'    !special for nContP = 1'
WRITE(LUMOD3,'(A)')'    !special for nContQ = 1'
WRITE(LUMOD3,'(A)')'    !special for nContA = 1 ...'
WRITE(LUMOD3,'(A)')'    !memory should be c**4*p**4 + p**4*L**6 which is fine'
WRITE(LUMOD3,'(A)')'    !this would be a simple sum for segmentet! or maybe the sum can be moved into the previous electron transfer reccurence'
WRITE(LUMOD3,'(A)')'    do iPassQ = 1,nPasses'
WRITE(LUMOD3,'(A)')'       do iContB=1,nContB'
WRITE(LUMOD3,'(A)')'          do iContA=1,nContA'
WRITE(LUMOD3,'(A)')'             iContP = iContA+(iContB-1)*nContA'
WRITE(LUMOD3,'(A)')'             do iContD=1,nContD'
WRITE(LUMOD3,'(A)')'                do iContC=1,nContC'
WRITE(LUMOD3,'(A)')'                   iContQ = iContC+(iContD-1)*nContC'
WRITE(LUMOD3,'(A)')'                   do iTUV=1,nTUVfull'
WRITE(LUMOD3,'(A)')'                      AUXarrayCont(iTUV,iContQ,iContP,iPassQ) = 0.0E0_realk'
WRITE(LUMOD3,'(A)')'                   enddo'
WRITE(LUMOD3,'(A)')'                enddo'
WRITE(LUMOD3,'(A)')'             enddo'
WRITE(LUMOD3,'(A)')'          enddo'
WRITE(LUMOD3,'(A)')'       enddo'
WRITE(LUMOD3,'(A)')'    enddo'
WRITE(LUMOD3,'(A)')'    do iPassQ = 1,nPasses'
WRITE(LUMOD3,'(A)')'       do iContB=1,nContB'
WRITE(LUMOD3,'(A)')'          do iPrimB=1,nPrimB'
WRITE(LUMOD3,'(A)')'             B = BCC(iPrimB,iContB)'
WRITE(LUMOD3,'(A)')'             do iContA=1,nContA'
WRITE(LUMOD3,'(A)')'                iContP = iContA+(iContB-1)*nContA'
WRITE(LUMOD3,'(A)')'                do iPrimA=1,nPrimA'
WRITE(LUMOD3,'(A)')'                   iPrimP = iPrimA + (iPrimB-1)*nPrimA'
WRITE(LUMOD3,'(A)')'                   ABTMP = ACC(iPrimA,iContA)*B'
WRITE(LUMOD3,'(A)')'                   do iContD=1,nContD'
WRITE(LUMOD3,'(A)')'                      do iPrimD=1,nPrimD'
WRITE(LUMOD3,'(A)')'                         ABDTMP = DCC(iPrimD,iContD)*ABTMP'
WRITE(LUMOD3,'(A)')'                         do iContC=1,nContC'
WRITE(LUMOD3,'(A)')'                            iContQ = iContC+(iContD-1)*nContC'
WRITE(LUMOD3,'(A)')'                            iPrimQ = (iPrimD-1)*nPrimC'
WRITE(LUMOD3,'(A)')'                            do iPrimC=1,nPrimC'
WRITE(LUMOD3,'(A)')'                               ABCDTMP = CCC(iPrimC,iContC)*ABDTMP'
WRITE(LUMOD3,'(A)')'                               iPrimQ = iPrimQ + 1'
WRITE(LUMOD3,'(A)')'                               do iTUV=1,nTUVfull'
WRITE(LUMOD3,'(A)')'                                  AUXarrayCont(iTUV,iContQ,iContP,iPassQ) = AUXarrayCont(iTUV,iContQ,iContP,iPassQ) + &'
WRITE(LUMOD3,'(A)')'                                       & ABCDTMP*AUXarray2(iTUV,iPrimQ,iPrimP,iPassQ)'
WRITE(LUMOD3,'(A)')'                               enddo'
WRITE(LUMOD3,'(A)')'                            enddo'
WRITE(LUMOD3,'(A)')'                         enddo'
WRITE(LUMOD3,'(A)')'                      enddo'
WRITE(LUMOD3,'(A)')'                   enddo'
WRITE(LUMOD3,'(A)')'                enddo'
WRITE(LUMOD3,'(A)')'             enddo'
WRITE(LUMOD3,'(A)')'          enddo'
WRITE(LUMOD3,'(A)')'       enddo'
WRITE(LUMOD3,'(A)')'    enddo'
WRITE(LUMOD3,'(A)')'  end subroutine PrimitiveContraction'
WRITE(LUMOD3,'(A)')''
WRITE(LUMOD3,'(A)')'  subroutine PrimitiveContractionSeg(AUXarray2,AUXarrayCont,nTUVfull,nPrimP,nPrimQ,nPasses)'
WRITE(LUMOD3,'(A)')'    implicit none'
WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nTUVfull,nPrimP,nPrimQ,nPasses'
WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray2(nTUVfull,nPrimQ,nPrimP,nPasses)'
WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nTUVfull,nPasses)'
WRITE(LUMOD3,'(A)')'    !'
WRITE(LUMOD3,'(A)')'    integer :: iPassQ,iTUV,iPrimQ,iPrimP'
WRITE(LUMOD3,'(A)')'    !Maybe the sum can be moved into the previous electron transfer reccurence or Vertical '
WRITE(LUMOD3,'(A)')'    do iPassQ = 1,nPasses'
WRITE(LUMOD3,'(A)')'     do iTUV=1,nTUVfull'
WRITE(LUMOD3,'(A)')'      AUXarrayCont(iTUV,iPassQ) = 0.0E0_realk'
WRITE(LUMOD3,'(A)')'     enddo'
WRITE(LUMOD3,'(A)')'    enddo'
WRITE(LUMOD3,'(A)')'    do iPassQ = 1,nPasses'
WRITE(LUMOD3,'(A)')'     do iPrimP = 1,nPrimP'
WRITE(LUMOD3,'(A)')'      do iPrimQ = 1,nPrimQ'
WRITE(LUMOD3,'(A)')'       do iTUV=1,nTUVfull'
WRITE(LUMOD3,'(A)')'        AUXarrayCont(iTUV,iPassQ) = AUXarrayCont(iTUV,iPassQ) &'
WRITE(LUMOD3,'(A)')'                                  & + AUXarray2(iTUV,iPrimQ,iPrimP,iPassQ)'
WRITE(LUMOD3,'(A)')'       enddo'
WRITE(LUMOD3,'(A)')'      enddo'
WRITE(LUMOD3,'(A)')'     enddo'
WRITE(LUMOD3,'(A)')'    enddo'
WRITE(LUMOD3,'(A)')'  end subroutine PrimitiveContractionSeg'
WRITE(LUMOD3,'(A)')''
!!$WRITE(LUMOD3,'(A)')'  subroutine build_Rpa(nPrimP,Pcent,Acent,Rpa)'
!!$WRITE(LUMOD3,'(A)')'    implicit none'
!!$WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP'
!!$WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: Pcent(3,nPrimP),Acent(3)'
!!$WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: Rpa(3,nPrimP)'
!!$WRITE(LUMOD3,'(A)')'    !'
!!$WRITE(LUMOD3,'(A)')'    integer :: iPrimP'
!!$WRITE(LUMOD3,'(A)')'    real(realk) :: Ax,Ay,Az'
!!$WRITE(LUMOD3,'(A)')'    Ax = -Acent(1)'
!!$WRITE(LUMOD3,'(A)')'    Ay = -Acent(2)'
!!$WRITE(LUMOD3,'(A)')'    Az = -Acent(3)'
!!$WRITE(LUMOD3,'(A)')'    DO iPrimP=1, nPrimP'
!!$WRITE(LUMOD3,'(A)')'       Rpa(1,iPrimP) = Pcent(1,iPrimP) + Ax'
!!$WRITE(LUMOD3,'(A)')'       Rpa(2,iPrimP) = Pcent(2,iPrimP) + Ay'
!!$WRITE(LUMOD3,'(A)')'       Rpa(3,iPrimP) = Pcent(3,iPrimP) + Az'
!!$WRITE(LUMOD3,'(A)')'    ENDDO'
!!$WRITE(LUMOD3,'(A)')'  end subroutine build_Rpa'
!!$WRITE(LUMOD3,'(A)')''
!!$WRITE(LUMOD3,'(A)')'  subroutine build_QP_squaredDistance_and_Rpq(nPrimQ,nPasses,nPrimP,Qcent,Pcent,squaredDistance,Rpq)'
!!$WRITE(LUMOD3,'(A)')'    implicit none'
!!$WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimQ,nPasses,nPrimP'
!!$WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: Qcent(3,nPrimQ,nPasses),Pcent(3,nPrimP)'
!!$WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: squaredDistance(nPrimQ,nPrimP,nPasses)'
!!$WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: Rpq(3,nPrimQ,nPrimP,nPasses)'
!!$WRITE(LUMOD3,'(A)')'    !'
!!$WRITE(LUMOD3,'(A)')'    integer :: iPrimQ,iPassQ,iPrimP'
!!$WRITE(LUMOD3,'(A)')'    real(realk) :: px,py,pz,pqx,pqy,pqz'
!!$WRITE(LUMOD3,'(A)')'    DO iPassQ=1, nPasses'
!!$WRITE(LUMOD3,'(A)')'       DO iPrimP=1, nPrimP'
!!$WRITE(LUMOD3,'(A)')'          px = Pcent(1,iPrimP)'
!!$WRITE(LUMOD3,'(A)')'          py = Pcent(2,iPrimP)'
!!$WRITE(LUMOD3,'(A)')'          pz = Pcent(3,iPrimP)'
!!$WRITE(LUMOD3,'(A)')'          DO iPrimQ=1, nPrimQ'
!!$WRITE(LUMOD3,'(A)')'             pqx = px - Qcent(1,iPrimQ,iPassQ)'
!!$WRITE(LUMOD3,'(A)')'             pqy = py - Qcent(2,iPrimQ,iPassQ)'
!!$WRITE(LUMOD3,'(A)')'             pqz = pz - Qcent(3,iPrimQ,iPassQ)'
!!$WRITE(LUMOD3,'(A)')'             rPQ(1,iPrimQ,iPrimP,iPassQ) = pqx'
!!$WRITE(LUMOD3,'(A)')'             rPQ(2,iPrimQ,iPrimP,iPassQ) = pqy'
!!$WRITE(LUMOD3,'(A)')'             rPQ(3,iPrimQ,iPrimP,iPassQ) = pqz'
!!$WRITE(LUMOD3,'(A)')'             squaredDistance(iPrimQ,iPrimP,iPassQ) = pqx*pqx+pqy*pqy+pqz*pqz'
!!$WRITE(LUMOD3,'(A)')'          ENDDO'
!!$WRITE(LUMOD3,'(A)')'       ENDDO'
!!$WRITE(LUMOD3,'(A)')'    ENDDO'
!!$WRITE(LUMOD3,'(A)')'  end subroutine build_QP_squaredDistance_and_Rpq'
!WRITE(LUMOD3,'(A)')''
WRITE(LUMOD3,'(A)')'  SUBROUTINE buildRJ000_general(nPasses,nPrimQ,nPrimP,nTABFJW1,nTABFJW2,reducedExponents,&'
WRITE(LUMOD3,'(A)')'       & TABFJW,RJ000,JMAX,Pcent,Qcent)'
WRITE(LUMOD3,'(A)')'    IMPLICIT NONE'
WRITE(LUMOD3,'(A)')'    INTEGER,intent(in)         :: nPrimP,nPrimQ,Jmax,nTABFJW1,nTABFJW2,nPasses'
WRITE(LUMOD3,'(A)')'    REAL(REALK),intent(in)     :: reducedExponents(nPrimQ,nPrimP)'
WRITE(LUMOD3,'(A)')'    REAL(REALK),intent(in)     :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)'
WRITE(LUMOD3,'(A)')'    REAL(REALK),intent(in)     :: TABFJW(0:nTABFJW1,0:nTABFJW2)'
WRITE(LUMOD3,'(A)')'    REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrimQ*nPrimP,nPasses)'
WRITE(LUMOD3,'(A)')'    !'
WRITE(LUMOD3,'(A)')'    REAL(REALK)     :: D2JP36,WVAL'
WRITE(LUMOD3,'(A)')'    REAL(REALK),PARAMETER :: HALF =0.5E0_realk,D1=1E0_realk,D2 = 2E0_realk, D4 = 4E0_realk, D100=100E0_realk'
WRITE(LUMOD3,'(A)')'    Real(realk),parameter :: D12 = 12E0_realk, TENTH = 0.01E0_realk'
WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk'
WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk'
WRITE(LUMOD3,'(A)')'    Integer :: IPNT,J'
WRITE(LUMOD3,'(A)')'    Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA'
WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk'
WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk'
WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk'
WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk'
WRITE(LUMOD3,'(A)')'    Real(realk), parameter :: PI=3.14159265358979323846E0_realk'
WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk'
WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2'
WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI'
WRITE(LUMOD3,'(A)')'    Real(realk) :: W2,W3,R'
WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: SMALL = 1E-15_realk'
WRITE(LUMOD3,'(A)')'    REAL(REALK) :: PX,PY,PZ,PQX,PQY,PQZ,squaredDistance'
WRITE(LUMOD3,'(A)')'    integer :: iPassQ,iPQ,iPrimP,iPrimQ'
WRITE(LUMOD3,'(A)')'    !make for different values of JMAX => loop unroll  '
WRITE(LUMOD3,'(A)')'    !sorting? '
WRITE(LUMOD3,'(A)')'    D2JP36 = 2*JMAX + 36'
WRITE(LUMOD3,'(A)')'    DO iPassQ=1, nPasses'
WRITE(LUMOD3,'(A)')'      ipq = 0'
WRITE(LUMOD3,'(A)')'      DO iPrimP=1, nPrimP'
WRITE(LUMOD3,'(A)')'        px = Pcent(1,iPrimP)'
WRITE(LUMOD3,'(A)')'        py = Pcent(2,iPrimP)'
WRITE(LUMOD3,'(A)')'        pz = Pcent(3,iPrimP)'
WRITE(LUMOD3,'(A)')'        DO iPrimQ=1, nPrimQ'
WRITE(LUMOD3,'(A)')'          ipq = ipq + 1'
WRITE(LUMOD3,'(A)')'          pqx = px - Qcent(1,iPrimQ,iPassQ)'
WRITE(LUMOD3,'(A)')'          pqy = py - Qcent(2,iPrimQ,iPassQ)'
WRITE(LUMOD3,'(A)')'          pqz = pz - Qcent(3,iPrimQ,iPassQ)'
WRITE(LUMOD3,'(A)')'          squaredDistance = pqx*pqx+pqy*pqy+pqz*pqz'
WRITE(LUMOD3,'(A)')'          WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance'
WRITE(LUMOD3,'(A)')'          !  0 < WVAL < 0.000001'
WRITE(LUMOD3,'(A)')'          IF (ABS(WVAL) .LT. SMALL) THEN'         
WRITE(LUMOD3,'(A)')'             RJ000(0,ipq,ipassq) = D1'
WRITE(LUMOD3,'(A)')'             DO J=1,JMAX'
WRITE(LUMOD3,'(A)')'                RJ000(J,ipq,ipassq)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT'
WRITE(LUMOD3,'(A)')'             ENDDO'
WRITE(LUMOD3,'(A)')'             !  0 < WVAL < 12 '
WRITE(LUMOD3,'(A)')'          ELSE IF (WVAL .LT. D12) THEN'
WRITE(LUMOD3,'(A)')'             IPNT = NINT(D100*WVAL)'
WRITE(LUMOD3,'(A)')'             WDIFF = WVAL - TENTH*IPNT'
WRITE(LUMOD3,'(A)')'             W2    = WDIFF*WDIFF'
WRITE(LUMOD3,'(A)')'             W3    = W2*WDIFF'
WRITE(LUMOD3,'(A)')'             W2    = W2*COEF2'
WRITE(LUMOD3,'(A)')'             W3    = W3*COEF3'
WRITE(LUMOD3,'(A)')'             DO J=0,JMAX'
WRITE(LUMOD3,'(A)')'                R = TABFJW(J,IPNT)'
WRITE(LUMOD3,'(A)')'                R = R -TABFJW(J+1,IPNT)*WDIFF'
WRITE(LUMOD3,'(A)')'                R = R + TABFJW(J+2,IPNT)*W2'
WRITE(LUMOD3,'(A)')'                R = R + TABFJW(J+3,IPNT)*W3'
WRITE(LUMOD3,'(A)')'                RJ000(J,ipq,ipassq) = R'
WRITE(LUMOD3,'(A)')'             ENDDO'
WRITE(LUMOD3,'(A)')'             !  12 < WVAL <= (2J+36) '
WRITE(LUMOD3,'(A)')'          ELSE IF (WVAL.LE.D2JP36) THEN'
WRITE(LUMOD3,'(A)')'             REXPW = HALF*EXP(-WVAL)'
WRITE(LUMOD3,'(A)')'             RWVAL = D1/WVAL'
WRITE(LUMOD3,'(A)')'             GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))'
WRITE(LUMOD3,'(A)')'             RJ000(0,ipq,ipassq) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
WRITE(LUMOD3,'(A)')'             DO J=1,JMAX'
WRITE(LUMOD3,'(A)')'                RJ000(J,ipq,ipassq) = RWVAL*((J - HALF)*RJ000(J-1,ipq,ipassq)-REXPW)'
WRITE(LUMOD3,'(A)')'             ENDDO'
WRITE(LUMOD3,'(A)')'             !  (2J+36) < WVAL '
WRITE(LUMOD3,'(A)')'          ELSE'
WRITE(LUMOD3,'(A)')'             RWVAL = PID4/WVAL'
WRITE(LUMOD3,'(A)')'             RJ000(0,ipq,ipassq) = SQRT(RWVAL)'
WRITE(LUMOD3,'(A)')'             RWVAL = RWVAL*PID4I'
WRITE(LUMOD3,'(A)')'             DO J = 1, JMAX'
WRITE(LUMOD3,'(A)')'                RJ000(J,ipq,ipassq) = RWVAL*(J - HALF)*RJ000(J-1,ipq,ipassq)'
WRITE(LUMOD3,'(A)')'             ENDDO'
WRITE(LUMOD3,'(A)')'          ENDIF'
WRITE(LUMOD3,'(A)')'        ENDDO'
WRITE(LUMOD3,'(A)')'      ENDDO'
WRITE(LUMOD3,'(A)')'    ENDDO'
WRITE(LUMOD3,'(A)')'  END SUBROUTINE buildRJ000_general'
WRITE(LUMOD3,'(A)')'  '
  WRITE(LUMOD3,'(A)')'END MODULE IchorEriCoulombintegralOBSGeneralMod'
  
  close(unit = LUMOD3)
  
contains
  subroutine subroutineMain(LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
       & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,spherical)
    integer,intent(in) :: LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,nTUV,AngmomP,AngmomQ
    integer,intent(in) :: AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec
    character(len=9) :: STRINGIN,STRINGOUT,TMPSTRING
    logical :: spherical

    IF(nTUV.LT.10)THEN
!       WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ELSEIF(nTUV.LT.100)THEN
!       WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ELSEIF(nTUV.LT.1000)THEN
!       WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ELSEIF(nTUV.LT.10000)THEN
!       WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ENDIF
    IF(AngmomPQ.EQ.0)THEN
       WRITE(LUMOD3,'(A)')'        call VerticalRecurrence0(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMOD3,'(A)')'               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&'
       WRITE(LUMOD3,'(A,A,A)')'               & PpreExpFac,QpreExpFac,',STRINGOUT,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
       WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
    ELSE
       IF(AngmomP.GE.AngmomQ)THEN
          IF(AngmomA.GE.AngmomB)THEN
             !A Vertical recurrence
             IF(AngmomPQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ELSE
                WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ENDIF
             !swap 
             TMPSTRING = STRINGIN
             STRINGIN  = STRINGOUT
             STRINGOUT  = TMPSTRING
             !determine TransferRecurrence
             IF(AngmomQ.EQ.0)THEN
                WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
                !No reason for the Electron Transfer Recurrence Relation 
             ELSE
!                WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation '
!                call sub_alloc1(nTUVP,nTUVQ,STRINGOUT,LUMOD3) !call mem_alloc(STRINGOUT,nTUVP*nTUVQ,nPrimQP*nPasses)
                IF(AngmomC.GE.AngmomD)THEN
                   !A to C TransferRecurrence
                   SPEC = 'AtoC'
                ELSE
                   !A to D TransferRecurrence
                   SPEC = 'AtoD'
                ENDIF                
                IF(AngmomP.LT.10)THEN
                   IF(AngmomQ.LT.10)THEN
                      WRITE(LUMOD3,'(A,I1,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
                      WRITE(LUMOD3,'(A,I1,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ELSE
                   IF(AngmomQ.LT.10)THEN
                      WRITE(LUMOD3,'(A,I2,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
                      WRITE(LUMOD3,'(A,I2,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ENDIF
                WRITE(LUMOD3,'(A)')'               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&'                
                WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
!                WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
                !swap 
                TMPSTRING = STRINGIN
                STRINGIN  = STRINGOUT
                STRINGOUT  = TMPSTRING
             ENDIF
          ELSE
             !B Vertical recurrence
             IF(AngmomPQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ELSE
                WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ENDIF
             !swap 
             TMPSTRING = STRINGIN
             STRINGIN  = STRINGOUT
             STRINGOUT  = TMPSTRING
             !determine TransferRecurrence
             IF(AngmomQ.EQ.0)THEN
                WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
                !No reason for the Electron Transfer Recurrence Relation 
             ELSE
!                WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation '
!                call sub_alloc1(nTUVP,nTUVQ,STRINGOUT,LUMOD3) !call mem_alloc(STRINGOUT,nTUVP*nTUVQ,nPrimQP*nPasses)
                IF(AngmomC.GE.AngmomD)THEN
                   !B to C TransferRecurrence
                   SPEC = 'BtoC'
                ELSE
                   !B to D TransferRecurrence
                   SPEC = 'BtoD'
                ENDIF
                IF(AngmomP.LT.10)THEN
                   IF(AngmomQ.LT.10)THEN
                      WRITE(LUMOD3,'(A,I1,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
                      WRITE(LUMOD3,'(A,I1,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ELSE
                   IF(AngmomQ.LT.10)THEN
                      WRITE(LUMOD3,'(A,I2,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
                      WRITE(LUMOD3,'(A,I2,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ENDIF
                WRITE(LUMOD3,'(A)')'               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&'                
                WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
!                WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
                !swap 
                TMPSTRING = STRINGIN
                STRINGIN  = STRINGOUT
                STRINGOUT  = TMPSTRING
             ENDIF
          ENDIF
       ELSE
          IF(AngmomC.GE.AngmomD)THEN
             !C Vertical recurrence
             IF(AngmomPQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'C(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ELSE
                WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'C(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ENDIF
             !swap 
             TMPSTRING = STRINGIN
             STRINGIN  = STRINGOUT
             STRINGOUT  = TMPSTRING 
             !determine TransferRecurrence
             IF(AngmomP.EQ.0)THEN
                WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
                !No reason for the Electron Transfer Recurrence Relation 
             ELSE
!                WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation '
!                call sub_alloc1(nTUVP,nTUVQ,STRINGOUT,LUMOD3) !call mem_alloc(STRINGOUT,nTUVP*nTUVQ,nPrimQP*nPasses)
                IF(AngmomA.GE.AngmomB)THEN
                   !C to A TransferRecurrence
                   SPEC = 'CtoA'
                ELSE
                   !C to B TransferRecurrence
                   SPEC = 'CtoB'
                ENDIF
                IF(AngmomP.LT.10)THEN
                   IF(AngmomQ.LT.10)THEN
                      WRITE(LUMOD3,'(A,I1,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
                      WRITE(LUMOD3,'(A,I1,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ELSE
                   IF(AngmomQ.LT.10)THEN
                      WRITE(LUMOD3,'(A,I2,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
                      WRITE(LUMOD3,'(A,I2,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ENDIF
                WRITE(LUMOD3,'(A)')'               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&'                
                WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
!                WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
                !swap 
                TMPSTRING = STRINGIN
                STRINGIN  = STRINGOUT
                STRINGOUT  = TMPSTRING
             ENDIF
          ELSE
             !D Vertical recurrence
             IF(AngmomPQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'D(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ELSE
                WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'D(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ENDIF
             !swap 
             TMPSTRING = STRINGIN
             STRINGIN  = STRINGOUT
             STRINGOUT  = TMPSTRING             
             IF(AngmomP.EQ.0)THEN
                WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
                !No reason for the Electron Transfer Recurrence Relation 
             ELSE
!                WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation '
!                call sub_alloc1(nTUVP,nTUVQ,STRINGOUT,LUMOD3) !call mem_alloc(STRINGOUT,nTUVP*nTUVQ,nPrimQP*nPasses)
                IF(AngmomA.GE.AngmomB)THEN
                   !D to A TransferRecurrence
                   SPEC = 'DtoA'
                ELSE
                   !D to B TransferRecurrence
                   SPEC = 'DtoB'
                ENDIF
                IF(AngmomP.LT.10)THEN
                   IF(AngmomQ.LT.10)THEN
                      WRITE(LUMOD3,'(A,I1,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
                      WRITE(LUMOD3,'(A,I1,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ELSE
                   IF(AngmomQ.LT.10)THEN
                      WRITE(LUMOD3,'(A,I2,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
                      WRITE(LUMOD3,'(A,I2,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ENDIF
                WRITE(LUMOD3,'(A)')'               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&'                
                WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
!                WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
                !swap 
                TMPSTRING = STRINGIN
                STRINGIN  = STRINGOUT
                STRINGOUT  = TMPSTRING
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    WRITE(LUMOD3,'(A)')'        nContQP = nContQ*nContP'
    IF(nTUVP*nTUVQ.LT.10)THEN
!       WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.100)THEN
!       WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.1000)THEN
!       WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.10000)THEN
!       WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ENDIF

    WRITE(LUMOD3,'(A)')'        IF(Qsegmented.AND.Psegmented)THEN'
    WRITE(LUMOD3,'(A,A,A,A,A,I4,A)')'         call PrimitiveContractionSeg(',STRINGIN,',',STRINGOUT,',',nTUVQ*nTUVP,',nPrimP,nPrimQ,nPasses)'
    WRITE(LUMOD3,'(A)')'        ELSE'
    WRITE(LUMOD3,'(A,A,A,A,A,I4,A)')'         call PrimitiveContraction(',STRINGIN,',',STRINGOUT,',',nTUVQ*nTUVP,',nPrimP,nPrimQ,nPasses,&'
    WRITE(LUMOD3,'(A)')'              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&'
    WRITE(LUMOD3,'(A)')'              & nContC,nPrimD,nContD)'
    WRITE(LUMOD3,'(A)')'        ENDIF'
!    WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'

    !swap 
    TMPSTRING = STRINGIN
    STRINGIN  = STRINGOUT
    STRINGOUT  = TMPSTRING
    !     The horizontal recurrence also extract nTUVspec from nTUV       
    IF(AngmomP.EQ.0)THEN
       WRITE(LUMOD3,'(A)')'        !no need for LHS Horizontal recurrence relations a simply copy'
    ELSE
!     WRITE(LUMOD3,'(A)')'        !LHS Horizontal recurrence relations '
     IF(nTUVAspec*nTUVBspec*nTUVQ.LT.10)THEN
!       WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
     ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.100)THEN
!       WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
     ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.1000)THEN
!       WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
     ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.10000)THEN
!       WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
     ENDIF
     IF(AngmomA.GE.AngmomB)THEN
        SPEC = 'AtoB'
     ELSE
        SPEC = 'BtoA'
     ENDIF
     IF(AngmomP.LT.10)THEN
       IF(AngmomA.LT.10)THEN
          IF(AngmomB.LT.10)THEN
             WRITE(LUMOD3,'(A,I1,A,I1,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ELSEIF(AngmomB.LT.100)THEN
             WRITE(LUMOD3,'(A,I1,A,I1,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ENDIF
       ELSEIF(AngmomA.LT.100)THEN
          IF(AngmomB.LT.10)THEN
             WRITE(LUMOD3,'(A,I1,A,I2,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ELSEIF(AngmomB.LT.100)THEN
             WRITE(LUMOD3,'(A,I1,A,I2,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ENDIF
       ENDIF
     ELSEIF(AngmomP.LT.100)THEN
       IF(AngmomA.LT.10)THEN
          IF(AngmomB.LT.10)THEN
             WRITE(LUMOD3,'(A,I2,A,I1,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ELSEIF(AngmomB.LT.100)THEN
             WRITE(LUMOD3,'(A,I2,A,I1,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ENDIF
       ELSEIF(AngmomA.LT.100)THEN
          IF(AngmomB.LT.10)THEN
             WRITE(LUMOD3,'(A,I2,A,I2,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ELSEIF(AngmomB.LT.100)THEN
             WRITE(LUMOD3,'(A,I2,A,I2,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ENDIF
       ENDIF
     ENDIF
!     WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
     !swap 
     TMPSTRING = STRINGIN
     STRINGIN  = STRINGOUT
     STRINGOUT  = TMPSTRING
    ENDIF


    IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
!       WRITE(LUMOD3,'(A)')'        !Spherical Transformation LHS'         
       IF(nlmA*nlmB*nTUVQ.LT.10)THEN
!          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.100)THEN
!          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.1000)THEN
!          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.10000)THEN
!          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ENDIF
       IF(AngmomP.LT.10)THEN
          IF(AngmomA.LT.10)THEN
             WRITE(LUMOD3,'(A,I1,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomA.LT.100)THEN
             WRITE(LUMOD3,'(A,I1,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ELSEIF(AngmomP.LT.100)THEN
          IF(AngmomA.LT.10)THEN
             WRITE(LUMOD3,'(A,I2,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomA.LT.100)THEN
             WRITE(LUMOD3,'(A,I2,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ENDIF
!       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ELSE
       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation LHS needed'
    ENDIF


    IF(AngmomQ.EQ.0)THEN
       WRITE(LUMOD3,'(A)')'        !no need for RHS Horizontal recurrence relations '
    ELSE
!       WRITE(LUMOD3,'(A)')'        !RHS Horizontal recurrence relations '
       IF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.10)THEN
!          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.100)THEN
!          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.1000)THEN
!          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.10000)THEN
!          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ENDIF
       IF(AngmomC.GE.AngmomD)THEN
          SPEC = 'CtoD'
       ELSE
          SPEC = 'DtoC'
       ENDIF
       IF(AngmomQ.LT.10)THEN
          IF(AngmomC.LT.10)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A,I1,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I1,A,I1,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ELSEIF(AngmomC.LT.100)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A,I2,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I1,A,I2,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ENDIF
       ELSEIF(AngmomQ.LT.100)THEN
          IF(AngmomC.LT.10)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I2,A,I1,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I2,A,I1,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ELSEIF(AngmomC.LT.100)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I2,A,I2,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I2,A,I2,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ENDIF
       ENDIF
!       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF

    IF(Spherical.AND.(AngmomC.GT.1.OR.AngmomD.GT.1))THEN
!       WRITE(LUMOD3,'(A)')'        !Spherical Transformation RHS'
       STRINGOUT  = 'CDAB     '
       IF(nlmA*nlmB*nlmC*nlmD.LT.10)THEN
!          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.100)THEN
!          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.1000)THEN
!          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.10000)THEN
!          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ENDIF
       IF(AngmomQ.LT.10)THEN
          IF(AngmomC.LT.10)THEN
             WRITE(LUMOD3,'(A,I1,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomC.LT.100)THEN
             WRITE(LUMOD3,'(A,I1,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ELSEIF(AngmomQ.LT.100)THEN
          IF(AngmomC.LT.10)THEN
             WRITE(LUMOD3,'(A,I2,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomC.LT.100)THEN
             WRITE(LUMOD3,'(A,I2,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ENDIF
!       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
!       TMPSTRING = STRINGIN
!       STRINGIN  = STRINGOUT
!       STRINGOUT  = TMPSTRING
!       WRITE(LUMOD3,'(A,A)')'        CDAB = ',STRINGIN
    ELSE
       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation RHS needed'
       WRITE(LUMOD3,'(A,A)')'        CDAB = ',STRINGIN
    ENDIF

!    WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
  end subroutine subroutineMain

  subroutine sub_alloc1(nTUVP,nTUVQ,STRINGOUT,LUMOD3)
    implicit none
    integer :: nTUVP,nTUVQ,LUMOD3
    character(len=9) :: STRINGOUT                  
    IF(nTUVP*nTUVQ.LT.10)THEN
       WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.100)THEN
       WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.1000)THEN
       WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.10000)THEN
       WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
    ENDIF
  end subroutine sub_alloc1

  subroutine sub_maxalloc1(nTUVP,nTUVQ,STRINGOUT,LUMOD3)
    implicit none
    integer :: nTUVP,nTUVQ,LUMOD3
    character(len=9) :: STRINGOUT                  
    IF(nTUVP*nTUVQ.LT.10)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I1,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVP*nTUVQ,'*nPrimQP)'
    ELSEIF(nTUVP*nTUVQ.LT.100)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I2,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVP*nTUVQ,'*nPrimQP)'
    ELSEIF(nTUVP*nTUVQ.LT.1000)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I3,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVP*nTUVQ,'*nPrimQP)'
    ELSEIF(nTUVP*nTUVQ.LT.10000)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I4,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVP*nTUVQ,'*nPrimQP)'
    ENDIF
  end subroutine sub_maxalloc1
  
  subroutine determineSizes(LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
       & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,spherical)
    integer,intent(in) :: LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,nTUV,AngmomP,AngmomQ
    integer,intent(in) :: AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec
    character(len=9) :: STRINGIN,STRINGOUT,TMPSTRING
    logical :: spherical

    IF(nTUV.LT.10)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I1,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUV,'*nPrimQP)'
    ELSEIF(nTUV.LT.100)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I2,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUV,'*nPrimQP)'
    ELSEIF(nTUV.LT.1000)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I3,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUV,'*nPrimQP)'
    ELSEIF(nTUV.LT.10000)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I4,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUV,'*nPrimQP)'
    ENDIF
    IF(AngmomPQ.EQ.0)THEN
!       WRITE(LUMOD3,'(A)')'        call VerticalRecurrence0(nPasses,nPrimP,nPrimQ,&'
!       WRITE(LUMOD3,'(A)')'               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&'
!       WRITE(LUMOD3,'(A,A,A)')'               & PpreExpFac,QpreExpFac,',STRINGOUT,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
!       WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
    ELSE
       IF(AngmomP.GE.AngmomQ)THEN
          IF(AngmomA.GE.AngmomB)THEN
             !A Vertical recurrence
             IF(AngmomPQ.LT.10)THEN
!                WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
!                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ELSE
!                WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'A(nPasses,nPrimP,nPrimQ,reducedExponents,&'
!                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Pexp,Acenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ENDIF
             !swap 
             TMPSTRING = STRINGIN
             STRINGIN  = STRINGOUT
             STRINGOUT  = TMPSTRING
             !determine TransferRecurrence
             IF(AngmomQ.EQ.0)THEN
!                WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
                !No reason for the Electron Transfer Recurrence Relation 
             ELSE
!                WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation '
                call sub_maxalloc1(nTUVP,nTUVQ,STRINGOUT,LUMOD3) !call mem_alloc(STRINGOUT,nTUVP*nTUVQ,nPrimQP*nPasses)
                IF(AngmomC.GE.AngmomD)THEN
                   !A to C TransferRecurrence
                   SPEC = 'AtoC'
                ELSE
                   !A to D TransferRecurrence
                   SPEC = 'AtoD'
                ENDIF                
                IF(AngmomP.LT.10)THEN
                   IF(AngmomQ.LT.10)THEN
!                      WRITE(LUMOD3,'(A,I1,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
!                      WRITE(LUMOD3,'(A,I1,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ELSE
                   IF(AngmomQ.LT.10)THEN
!                      WRITE(LUMOD3,'(A,I2,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
!                      WRITE(LUMOD3,'(A,I2,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ENDIF
!                WRITE(LUMOD3,'(A)')'               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&'                
!                WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
!                WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
                !swap 
                TMPSTRING = STRINGIN
                STRINGIN  = STRINGOUT
                STRINGOUT  = TMPSTRING
             ENDIF
          ELSE
             !B Vertical recurrence
             IF(AngmomPQ.LT.10)THEN
!                WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
!                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ELSE
!                WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'B(nPasses,nPrimP,nPrimQ,reducedExponents,&'
!                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Pexp,Bcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ENDIF
             !swap 
             TMPSTRING = STRINGIN
             STRINGIN  = STRINGOUT
             STRINGOUT  = TMPSTRING
             !determine TransferRecurrence
             IF(AngmomQ.EQ.0)THEN
!                WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
                !No reason for the Electron Transfer Recurrence Relation 
             ELSE
!                WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation '
                call sub_maxalloc1(nTUVP,nTUVQ,STRINGOUT,LUMOD3) !call mem_alloc(STRINGOUT,nTUVP*nTUVQ,nPrimQP*nPasses)
                IF(AngmomC.GE.AngmomD)THEN
                   !B to C TransferRecurrence
                   SPEC = 'BtoC'
                ELSE
                   !B to D TransferRecurrence
                   SPEC = 'BtoD'
                ENDIF
                IF(AngmomP.LT.10)THEN
                   IF(AngmomQ.LT.10)THEN
 !                     WRITE(LUMOD3,'(A,I1,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
 !                     WRITE(LUMOD3,'(A,I1,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ELSE
                   IF(AngmomQ.LT.10)THEN
 !                     WRITE(LUMOD3,'(A,I2,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
 !                     WRITE(LUMOD3,'(A,I2,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ENDIF
 !               WRITE(LUMOD3,'(A)')'               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&'                
 !               WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
 !               WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
                !swap 
                TMPSTRING = STRINGIN
                STRINGIN  = STRINGOUT
                STRINGOUT  = TMPSTRING
             ENDIF
          ENDIF
       ELSE
          IF(AngmomC.GE.AngmomD)THEN
             !C Vertical recurrence
             IF(AngmomPQ.LT.10)THEN
  !              WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'C(nPasses,nPrimP,nPrimQ,reducedExponents,&'
  !              WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ELSE
  !              WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'C(nPasses,nPrimP,nPrimQ,reducedExponents,&'
  !              WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ENDIF
             !swap 
             TMPSTRING = STRINGIN
             STRINGIN  = STRINGOUT
             STRINGOUT  = TMPSTRING 
             !determine TransferRecurrence
             IF(AngmomP.EQ.0)THEN
!                WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
                !No reason for the Electron Transfer Recurrence Relation 
             ELSE
!                WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation '
                call sub_maxalloc1(nTUVP,nTUVQ,STRINGOUT,LUMOD3) !call mem_alloc(STRINGOUT,nTUVP*nTUVQ,nPrimQP*nPasses)
                IF(AngmomA.GE.AngmomB)THEN
                   !C to A TransferRecurrence
                   SPEC = 'CtoA'
                ELSE
                   !C to B TransferRecurrence
                   SPEC = 'CtoB'
                ENDIF
                IF(AngmomP.LT.10)THEN
                   IF(AngmomQ.LT.10)THEN
!                      WRITE(LUMOD3,'(A,I1,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
!                      WRITE(LUMOD3,'(A,I1,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ELSE
                   IF(AngmomQ.LT.10)THEN
!                      WRITE(LUMOD3,'(A,I2,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
!                      WRITE(LUMOD3,'(A,I2,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ENDIF
!                WRITE(LUMOD3,'(A)')'               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&'                
!                WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
!                WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
                !swap 
                TMPSTRING = STRINGIN
                STRINGIN  = STRINGOUT
                STRINGOUT  = TMPSTRING
             ENDIF
          ELSE
             !D Vertical recurrence
             IF(AngmomPQ.LT.10)THEN
!                WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'D(nPasses,nPrimP,nPrimQ,reducedExponents,&'
!                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ELSE
!                WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'D(nPasses,nPrimP,nPrimQ,reducedExponents,&'
!                WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
             ENDIF
             !swap 
             TMPSTRING = STRINGIN
             STRINGIN  = STRINGOUT
             STRINGOUT  = TMPSTRING             
             IF(AngmomP.EQ.0)THEN
!                WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
                !No reason for the Electron Transfer Recurrence Relation 
             ELSE
!                WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation '
                call sub_maxalloc1(nTUVP,nTUVQ,STRINGOUT,LUMOD3) !call mem_alloc(STRINGOUT,nTUVP*nTUVQ,nPrimQP*nPasses)
                IF(AngmomA.GE.AngmomB)THEN
                   !D to A TransferRecurrence
                   SPEC = 'DtoA'
                ELSE
                   !D to B TransferRecurrence
                   SPEC = 'DtoB'
                ENDIF
                IF(AngmomP.LT.10)THEN
                   IF(AngmomQ.LT.10)THEN
!                      WRITE(LUMOD3,'(A,I1,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
!                      WRITE(LUMOD3,'(A,I1,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ELSE
                   IF(AngmomQ.LT.10)THEN
!                      WRITE(LUMOD3,'(A,I2,A,I1,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ELSE
!                      WRITE(LUMOD3,'(A,I2,A,I2,A4,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,SPEC,'(nPasses,nPrimP,nPrimQ,reducedExponents,&'
                   ENDIF
                ENDIF
!                WRITE(LUMOD3,'(A)')'               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&'                
!                WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
!                WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
                !swap 
                TMPSTRING = STRINGIN
                STRINGIN  = STRINGOUT
                STRINGOUT  = TMPSTRING
             ENDIF
          ENDIF
       ENDIF
    ENDIF
!    WRITE(LUMOD3,'(A)')'        nContQP = nContQ*nContP'

    IF(nTUVP*nTUVQ.LT.10)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I1,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVP*nTUVQ,'*nPrimQP)'
    ELSEIF(nTUVP*nTUVQ.LT.100)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I2,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVP*nTUVQ,'*nContQP)'
    ELSEIF(nTUVP*nTUVQ.LT.1000)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I3,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVP*nTUVQ,'*nContQP)'
    ELSEIF(nTUVP*nTUVQ.LT.10000)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I4,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVP*nTUVQ,'*nContQP)'
    ENDIF

!    WRITE(LUMOD3,'(A)')'        IF(Qsegmented.AND.Psegmented)THEN'
!    WRITE(LUMOD3,'(A,A,A,A,A,I4,A)')'         call PrimitiveContractionSeg(',STRINGIN,',',STRINGOUT,',',nTUVQ*nTUVP,',nPrimP,nPrimQ,nPasses)'
!    WRITE(LUMOD3,'(A)')'        ELSE'
!    WRITE(LUMOD3,'(A,A,A,A,A,I4,A)')'         call PrimitiveContraction(',STRINGIN,',',STRINGOUT,',',nTUVQ*nTUVP,',nPrimP,nPrimQ,nPasses,&'
!    WRITE(LUMOD3,'(A)')'              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&'
!    WRITE(LUMOD3,'(A)')'              & nContC,nPrimD,nContD)'
!    WRITE(LUMOD3,'(A)')'        ENDIF'
!    WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'

    !swap 
    TMPSTRING = STRINGIN
    STRINGIN  = STRINGOUT
    STRINGOUT  = TMPSTRING
    !     The horizontal recurrence also extract nTUVspec from nTUV       
    IF(AngmomP.EQ.0)THEN
!       WRITE(LUMOD3,'(A)')'        !no need for LHS Horizontal recurrence relations a simply copy'
    ELSE
!     WRITE(LUMOD3,'(A)')'        !LHS Horizontal recurrence relations '
     IF(nTUVAspec*nTUVBspec*nTUVQ.LT.10)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I1,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP)'
     ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.100)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I2,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP)'
     ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.1000)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I3,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP)'
     ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.10000)THEN
       WRITE(LUMOD3,'(A7,A9,A,A9,A,I4,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP)'
     ENDIF
     IF(AngmomA.GE.AngmomB)THEN
        SPEC = 'AtoB'
     ELSE
        SPEC = 'BtoA'
     ENDIF
     IF(AngmomP.LT.10)THEN
       IF(AngmomA.LT.10)THEN
          IF(AngmomB.LT.10)THEN
!             WRITE(LUMOD3,'(A,I1,A,I1,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
!                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ELSEIF(AngmomB.LT.100)THEN
!             WRITE(LUMOD3,'(A,I1,A,I1,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
!                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ENDIF
       ELSEIF(AngmomA.LT.100)THEN
          IF(AngmomB.LT.10)THEN
!             WRITE(LUMOD3,'(A,I1,A,I2,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
!                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ELSEIF(AngmomB.LT.100)THEN
!             WRITE(LUMOD3,'(A,I1,A,I2,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
!                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ENDIF
       ENDIF
     ELSEIF(AngmomP.LT.100)THEN
       IF(AngmomA.LT.10)THEN
          IF(AngmomB.LT.10)THEN
!             WRITE(LUMOD3,'(A,I2,A,I1,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
!                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ELSEIF(AngmomB.LT.100)THEN
!             WRITE(LUMOD3,'(A,I2,A,I1,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
!                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ENDIF
       ELSEIF(AngmomA.LT.100)THEN
          IF(AngmomB.LT.10)THEN
!             WRITE(LUMOD3,'(A,I2,A,I2,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
!                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ELSEIF(AngmomB.LT.100)THEN
!             WRITE(LUMOD3,'(A,I2,A,I2,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
!                  &AngmomA,'B',AngmomB,SPEC,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
          ENDIF
       ENDIF
     ENDIF
!     WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
     !swap 
     TMPSTRING = STRINGIN
     STRINGIN  = STRINGOUT
     STRINGOUT  = TMPSTRING
    ENDIF


    IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
!       WRITE(LUMOD3,'(A)')'        !Spherical Transformation LHS'         
       IF(nlmA*nlmB*nTUVQ.LT.10)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I1,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nTUVQ,'*nContQP)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.100)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I2,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nTUVQ,'*nContQP)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.1000)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I3,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nTUVQ,'*nContQP)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.10000)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I4,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nTUVQ,'*nContQP)'
       ENDIF
       IF(AngmomP.LT.10)THEN
          IF(AngmomA.LT.10)THEN
!             WRITE(LUMOD3,'(A,I1,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomA.LT.100)THEN
!             WRITE(LUMOD3,'(A,I1,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ELSEIF(AngmomP.LT.100)THEN
          IF(AngmomA.LT.10)THEN
!             WRITE(LUMOD3,'(A,I2,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomA.LT.100)THEN
!             WRITE(LUMOD3,'(A,I2,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ENDIF
!       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ELSE
!       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation LHS needed'
    ENDIF


    IF(AngmomQ.EQ.0)THEN
!       WRITE(LUMOD3,'(A)')'        !no need for RHS Horizontal recurrence relations '
    ELSE
!       WRITE(LUMOD3,'(A)')'        !RHS Horizontal recurrence relations '
       IF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.10)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I1,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP)'
!          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.100)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I2,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.1000)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I3,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.10000)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I4,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP)'
       ENDIF
       IF(AngmomC.GE.AngmomD)THEN
          SPEC = 'CtoD'
       ELSE
          SPEC = 'DtoC'
       ENDIF
       IF(AngmomQ.LT.10)THEN
          IF(AngmomC.LT.10)THEN
             IF(AngmomD.LT.10)THEN
!                WRITE(LUMOD3,'(A,I1,A,I1,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
!                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
!                WRITE(LUMOD3,'(A,I1,A,I1,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
!                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ELSEIF(AngmomC.LT.100)THEN
             IF(AngmomD.LT.10)THEN
!                WRITE(LUMOD3,'(A,I1,A,I2,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
!                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
!                WRITE(LUMOD3,'(A,I1,A,I2,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
!                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ENDIF
       ELSEIF(AngmomQ.LT.100)THEN
          IF(AngmomC.LT.10)THEN
             IF(AngmomD.LT.10)THEN
!                WRITE(LUMOD3,'(A,I2,A,I1,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
!                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
!                WRITE(LUMOD3,'(A,I2,A,I1,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
!                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ELSEIF(AngmomC.LT.100)THEN
             IF(AngmomD.LT.10)THEN
!                WRITE(LUMOD3,'(A,I2,A,I2,A,I1,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
!                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
!                WRITE(LUMOD3,'(A,I2,A,I2,A,I2,A4,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
!                     &AngmomC,'D',AngmomD,SPEC,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ENDIF
       ENDIF
!       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF

    IF(Spherical.AND.(AngmomC.GT.1.OR.AngmomD.GT.1))THEN
!       WRITE(LUMOD3,'(A)')'        !Spherical Transformation RHS'

       IF(nlmA*nlmB*nlmC*nlmD.LT.10)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I1,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nlmC*nlmD,'*nContQP)'
!          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.100)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I2,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nlmC*nlmD,'*nContQP)'
!          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.1000)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I3,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nlmC*nlmD,'*nContQP)'
!          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.10000)THEN
          WRITE(LUMOD3,'(A7,A9,A,A9,A,I4,A)')'       ',STRINGOUT,'maxSize = MAX(',STRINGOUT,'maxSize,',nlmA*nlmB*nlmC*nlmD,'*nContQP)'
!          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ENDIF
       IF(AngmomQ.LT.10)THEN
          IF(AngmomC.LT.10)THEN
!             WRITE(LUMOD3,'(A,I1,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomC.LT.100)THEN
!             WRITE(LUMOD3,'(A,I1,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ELSEIF(AngmomQ.LT.100)THEN
          IF(AngmomC.LT.10)THEN
!             WRITE(LUMOD3,'(A,I2,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomC.LT.100)THEN
!             WRITE(LUMOD3,'(A,I2,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ENDIF
!       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ELSE
!       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation RHS needed'
    ENDIF

!    WRITE(LUMOD3,'(A,A)')'        CDAB = ',STRINGIN
!    WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
  end subroutine determineSizes

  !IF(AngmomD.GT.AngmomC.AND.(AngmomQ.GT.AngmomP))THEN
  subroutine subroutineMainD(LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
       & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,spherical)
    integer,intent(in) :: LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,nTUV,AngmomP,AngmomQ
    integer,intent(in) :: AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec
    character(len=9) :: STRINGIN,STRINGOUT,TMPSTRING
    logical :: spherical

    WRITE(LUMOD3,'(A,I2,A,I2,A,I2,A,I2,A)')'        !This is the Angmom(A=',AngmomA,',B=',AngmomB,',C=',AngmomC,',D=',AngmomD,') combi'
    IF(nTUV.LT.10)THEN
       WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ELSEIF(nTUV.LT.100)THEN
       WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ELSEIF(nTUV.LT.1000)THEN
       WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ELSEIF(nTUV.LT.10000)THEN
       WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ENDIF
    IF(AngmomPQ.EQ.0)THEN
       stop 'this should never happenA'
    ELSE
       IF(AngmomPQ.LT.10)THEN
          WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'D(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       ELSE
          WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'D(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       ENDIF
       WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Qexp,Dcenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
    ENDIF
!    WRITE(LUMOD3,'(A,I4,A)')'        call mem_ichor_dealloc(RJ000)'

    !swap 
    TMPSTRING = STRINGIN
    STRINGIN  = STRINGOUT
    STRINGOUT  = TMPSTRING

    IF(AngmomP.EQ.0)THEN
       WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
       !No reason for the Electron Transfer Recurrence Relation 
    ELSE
       IF(AngmomA.GE.AngmomB)THEN
          WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation from electron 4 to electron 1'
       ELSE
          WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation from electron 4 to electron 2'
       ENDIF
       IF(nTUVP*nTUVQ.LT.10)THEN
          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
       ELSEIF(nTUVP*nTUVQ.LT.100)THEN
          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
       ELSEIF(nTUVP*nTUVQ.LT.1000)THEN
          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
       ELSEIF(nTUVP*nTUVQ.LT.10000)THEN
          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
       ENDIF
       IF(AngmomA.GE.AngmomB)THEN
          IF(AngmomP.LT.10)THEN
             IF(AngmomQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A,I1,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(LUMOD3,'(A,I1,A,I2,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ELSE
             IF(AngmomQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I2,A,I1,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(LUMOD3,'(A,I2,A,I2,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'DtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ENDIF
          WRITE(LUMOD3,'(A)')'               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&'
          WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
          WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       ELSE
          WRITE(LUMOD3,'(A)')'        NOT IMPLEMENTED '
          WRITE(LUMOD3,'(A)')'        CALL ICHORQUIT(''NOT IMPLEMENTED'',-1)'          
          IF(AngmomP.LT.10)THEN
             IF(AngmomQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A,I1,A)')'   !     call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(LUMOD3,'(A,I1,A,I2,A)')'   !     call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ELSE
             IF(AngmomQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I2,A,I1,A)')'   !     call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(LUMOD3,'(A,I2,A,I2,A)')'   !     call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'DtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ENDIF
          WRITE(LUMOD3,'(A)')'         !      & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Dexp,nPrimA,nPrimB,nPrimC,nPrimD,&'
          WRITE(LUMOD3,'(A,A,A,A,A)')' !              & ',STRINGIN,',',STRINGOUT,')'
          WRITE(LUMOD3,'(A,A,A)')'     !   call mem_ichor_dealloc(',STRINGIN,')'
       ENDIF
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF

    WRITE(LUMOD3,'(A)')'        nContQP = nContQ*nContP'
    IF(nTUVP*nTUVQ.LT.10)THEN
       WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP)'
    ELSEIF(nTUVP*nTUVQ.LT.100)THEN
       WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP)'
    ELSEIF(nTUVP*nTUVQ.LT.1000)THEN
       WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP)'
    ELSEIF(nTUVP*nTUVQ.LT.10000)THEN
       WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP)'
    ENDIF

    WRITE(LUMOD3,'(A)')'        IF(Qsegmented.AND.Psegmented)THEN'
    WRITE(LUMOD3,'(A,A,A,A,A,I4,A)')'         call PrimitiveContractionSeg(',STRINGIN,',',STRINGOUT,',',nTUVQ*nTUVP,',nPrimP,nPrimQ,nPasses)'
    WRITE(LUMOD3,'(A)')'        ELSE'
    WRITE(LUMOD3,'(A,A,A,A,A,I4,A)')'         call PrimitiveContraction(',STRINGIN,',',STRINGOUT,',',nTUVQ*nTUVP,',nPrimP,nPrimQ,nPasses,&'
    WRITE(LUMOD3,'(A)')'              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&'
    WRITE(LUMOD3,'(A)')'              & nContC,nPrimD,nContD)'
    WRITE(LUMOD3,'(A)')'        ENDIF'
    WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'

    !swap 
    TMPSTRING = STRINGIN
    STRINGIN  = STRINGOUT
    STRINGOUT  = TMPSTRING
    !     The horizontal recurrence also extract nTUVspec from nTUV       
    IF(AngmomP.EQ.0)THEN
       WRITE(LUMOD3,'(A)')'        !no need for LHS Horizontal recurrence relations a simply copy'
    ELSE
      IF(AngmomA.GE.AngmomB)THEN
         WRITE(LUMOD3,'(A)')'        !LHS Horizontal recurrence relations A to B'
      ELSE !B gt A
         WRITE(LUMOD3,'(A)')'        !LHS Horizontal recurrence relations B to A'
      ENDIF
      IF(nTUVAspec*nTUVBspec*nTUVQ.LT.10)THEN
        WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
      ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.100)THEN
        WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
      ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.1000)THEN
        WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
      ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.10000)THEN
        WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
      ENDIF
      IF(AngmomP.LT.10)THEN
        IF(AngmomA.LT.10)THEN
           IF(AngmomB.LT.10)THEN
              WRITE(LUMOD3,'(A,I1,A,I1,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ELSEIF(AngmomB.LT.100)THEN
              WRITE(LUMOD3,'(A,I1,A,I1,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ENDIF
        ELSEIF(AngmomA.LT.100)THEN
           IF(AngmomB.LT.10)THEN
              WRITE(LUMOD3,'(A,I1,A,I2,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ELSEIF(AngmomB.LT.100)THEN
              WRITE(LUMOD3,'(A,I1,A,I2,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ENDIF
        ENDIF
      ELSEIF(AngmomP.LT.100)THEN
        IF(AngmomA.LT.10)THEN
           IF(AngmomB.LT.10)THEN
              WRITE(LUMOD3,'(A,I2,A,I1,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ELSEIF(AngmomB.LT.100)THEN
              WRITE(LUMOD3,'(A,I2,A,I1,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ENDIF
        ELSEIF(AngmomA.LT.100)THEN
           IF(AngmomB.LT.10)THEN
              WRITE(LUMOD3,'(A,I2,A,I2,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ELSEIF(AngmomB.LT.100)THEN
              WRITE(LUMOD3,'(A,I2,A,I2,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ENDIF
        ENDIF
      ENDIF
      WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
      !swap 
      TMPSTRING = STRINGIN
      STRINGIN  = STRINGOUT
      STRINGOUT  = TMPSTRING
    ENDIF


    IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
       WRITE(LUMOD3,'(A)')'        !Spherical Transformation LHS'         
       IF(nlmA*nlmB*nTUVQ.LT.10)THEN
          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.100)THEN
          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.1000)THEN
          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.10000)THEN
          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ENDIF
       IF(AngmomP.LT.10)THEN
          IF(AngmomA.LT.10)THEN
             WRITE(LUMOD3,'(A,I1,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomA.LT.100)THEN
             WRITE(LUMOD3,'(A,I1,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ELSEIF(AngmomP.LT.100)THEN
          IF(AngmomA.LT.10)THEN
             WRITE(LUMOD3,'(A,I2,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomA.LT.100)THEN
             WRITE(LUMOD3,'(A,I2,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ENDIF
       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ELSE
       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation LHS needed'
    ENDIF


    IF(AngmomQ.EQ.0)THEN
       stop 'this should never happenD'
    ELSE
       WRITE(LUMOD3,'(A)')'        !RHS Horizontal recurrence relations D to C'
       IF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.10)THEN
          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.100)THEN
          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.1000)THEN
          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.10000)THEN
          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ENDIF
       IF(AngmomQ.LT.10)THEN
          IF(AngmomC.LT.10)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A,I1,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'DtoC(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I1,A,I1,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'DtoC(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ELSEIF(AngmomC.LT.100)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A,I2,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'DtoC(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I1,A,I2,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'DtoC(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ENDIF
       ELSEIF(AngmomQ.LT.100)THEN
          IF(AngmomC.LT.10)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I2,A,I1,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'DtoC(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I2,A,I1,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'DtoC(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ELSEIF(AngmomC.LT.100)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I2,A,I2,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'DtoC(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I2,A,I2,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'DtoC(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ENDIF
       ENDIF
       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF

    IF(Spherical.AND.(AngmomC.GT.1.OR.AngmomD.GT.1))THEN
       WRITE(LUMOD3,'(A)')'        !Spherical Transformation RHS'

       IF(nlmA*nlmB*nlmC*nlmD.LT.10)THEN
          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.100)THEN
          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.1000)THEN
          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.10000)THEN
          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ENDIF
       IF(AngmomQ.LT.10)THEN
          IF(AngmomC.LT.10)THEN
             WRITE(LUMOD3,'(A,I1,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomC.LT.100)THEN
             WRITE(LUMOD3,'(A,I1,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ELSEIF(AngmomQ.LT.100)THEN
          IF(AngmomC.LT.10)THEN
             WRITE(LUMOD3,'(A,I2,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomC.LT.100)THEN
             WRITE(LUMOD3,'(A,I2,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ENDIF
       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ELSE
       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation RHS needed'
    ENDIF

    WRITE(LUMOD3,'(A,A)')'        CDAB = ',STRINGIN
    WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
  end subroutine subroutineMainD

  !IF(AngmomD.GT.AngmomC.AND.(AngmomQ.GT.AngmomP))THEN
  subroutine subroutineMainC(LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
       & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,spherical)
    integer,intent(in) :: LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,nTUV,AngmomP,AngmomQ
    integer,intent(in) :: AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec
    character(len=9) :: STRINGIN,STRINGOUT,TMPSTRING
    logical :: spherical

    WRITE(LUMOD3,'(A,I2,A,I2,A,I2,A,I2,A)')'        !This is the Angmom(A=',AngmomA,',B=',AngmomB,',C=',AngmomC,',D=',AngmomD,') combi'
    IF(nTUV.LT.10)THEN
       WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ELSEIF(nTUV.LT.100)THEN
       WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ELSEIF(nTUV.LT.1000)THEN
       WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ELSEIF(nTUV.LT.10000)THEN
       WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUV,'*nPrimQP*nPasses)'
    ENDIF
    IF(AngmomPQ.EQ.0)THEN
       stop 'this should never happenA'
    ELSE
       IF(AngmomPQ.LT.10)THEN
          WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'C(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       ELSE
          WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'C(nPasses,nPrimP,nPrimQ,reducedExponents,&'
       ENDIF
       WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Qexp,Ccenter,Pcent,Qcent,integralPrefactor,PpreExpFac,QpreExpFac,',STRINGOUT,')'
    ENDIF
!    WRITE(LUMOD3,'(A,I4,A)')'        call mem_ichor_dealloc(RJ000)'

    !swap 
    TMPSTRING = STRINGIN
    STRINGIN  = STRINGOUT
    STRINGOUT  = TMPSTRING

    IF(AngmomP.EQ.0)THEN
       WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
       !No reason for the Electron Transfer Recurrence Relation 
    ELSE
       IF(AngmomA.GE.AngmomB)THEN
          WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation from electron 3 to electron 1'
       ELSE
          WRITE(LUMOD3,'(A)')'        !Electron Transfer Recurrence Relation from electron 3 to electron 2'
       ENDIF
       IF(nTUVP*nTUVQ.LT.10)THEN
          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
       ELSEIF(nTUVP*nTUVQ.LT.100)THEN
          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
       ELSEIF(nTUVP*nTUVQ.LT.1000)THEN
          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
       ELSEIF(nTUVP*nTUVQ.LT.10000)THEN
          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nPrimQP*nPasses)'
       ENDIF
       IF(AngmomA.GE.AngmomB)THEN
          IF(AngmomP.LT.10)THEN
             IF(AngmomQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A,I1,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(LUMOD3,'(A,I1,A,I2,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ELSE
             IF(AngmomQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I2,A,I1,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(LUMOD3,'(A,I2,A,I2,A)')'        call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'CtoA(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ENDIF
          WRITE(LUMOD3,'(A)')'               & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&'
          WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
          WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       ELSE
          WRITE(LUMOD3,'(A)')'        NOT IMPLEMENTED '
          WRITE(LUMOD3,'(A)')'        CALL ICHORQUIT(''NOT IMPLEMENTED'',-1)'          
          IF(AngmomP.LT.10)THEN
             IF(AngmomQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A,I1,A)')'   !     call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(LUMOD3,'(A,I1,A,I2,A)')'   !     call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ELSE
             IF(AngmomQ.LT.10)THEN
                WRITE(LUMOD3,'(A,I2,A,I1,A)')'   !     call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ELSE
                WRITE(LUMOD3,'(A,I2,A,I2,A)')'   !     call TransferRecurrenceP',AngmomP,'Q',AngmomQ,'CtoB(nPasses,nPrimP,nPrimQ,reducedExponents,&'
             ENDIF
          ENDIF
          WRITE(LUMOD3,'(A)')'         !      & Pexp,Qexp,Pdistance12,Qdistance12,Bexp,Cexp,nPrimA,nPrimB,nPrimC,nPrimD,&'
          WRITE(LUMOD3,'(A,A,A,A,A)')' !              & ',STRINGIN,',',STRINGOUT,')'
          WRITE(LUMOD3,'(A,A,A)')'     !   call mem_ichor_dealloc(',STRINGIN,')'
       ENDIF
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF

    WRITE(LUMOD3,'(A)')'        nContQP = nContQ*nContP'
    IF(nTUVP*nTUVQ.LT.10)THEN
       WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.100)THEN
       WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.1000)THEN
       WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.10000)THEN
       WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ENDIF

    WRITE(LUMOD3,'(A)')'        IF(Qsegmented.AND.Psegmented)THEN'
    WRITE(LUMOD3,'(A,A,A,A,A,I4,A)')'         call PrimitiveContractionSeg(',STRINGIN,',',STRINGOUT,',',nTUVQ*nTUVP,',nPrimP,nPrimQ,nPasses)'
    WRITE(LUMOD3,'(A)')'        ELSE'
    WRITE(LUMOD3,'(A,A,A,A,A,I4,A)')'         call PrimitiveContraction(',STRINGIN,',',STRINGOUT,',',nTUVQ*nTUVP,',nPrimP,nPrimQ,nPasses,&'
    WRITE(LUMOD3,'(A)')'              & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,&'
    WRITE(LUMOD3,'(A)')'              & nContC,nPrimD,nContD)'
    WRITE(LUMOD3,'(A)')'        ENDIF'
    WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'

    !swap 
    TMPSTRING = STRINGIN
    STRINGIN  = STRINGOUT
    STRINGOUT  = TMPSTRING
    !     The horizontal recurrence also extract nTUVspec from nTUV       
    IF(AngmomP.EQ.0)THEN
       WRITE(LUMOD3,'(A)')'        !no need for LHS Horizontal recurrence relations a simply copy'
    ELSE
      IF(AngmomA.GE.AngmomB)THEN
         WRITE(LUMOD3,'(A)')'        !LHS Horizontal recurrence relations A to B'
      ELSE !B gt A
         WRITE(LUMOD3,'(A)')'        !LHS Horizontal recurrence relations B to A'
      ENDIF
      IF(nTUVAspec*nTUVBspec*nTUVQ.LT.10)THEN
        WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
      ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.100)THEN
        WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
      ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.1000)THEN
        WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
      ELSEIF(nTUVAspec*nTUVBspec*nTUVQ.LT.10000)THEN
        WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVAspec*nTUVBspec*nTUVQ,'*nContQP*nPasses)'
      ENDIF
      IF(AngmomP.LT.10)THEN
        IF(AngmomA.LT.10)THEN
           IF(AngmomB.LT.10)THEN
              WRITE(LUMOD3,'(A,I1,A,I1,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ELSEIF(AngmomB.LT.100)THEN
              WRITE(LUMOD3,'(A,I1,A,I1,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ENDIF
        ELSEIF(AngmomA.LT.100)THEN
           IF(AngmomB.LT.10)THEN
              WRITE(LUMOD3,'(A,I1,A,I2,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ELSEIF(AngmomB.LT.100)THEN
              WRITE(LUMOD3,'(A,I1,A,I2,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ENDIF
        ENDIF
      ELSEIF(AngmomP.LT.100)THEN
        IF(AngmomA.LT.10)THEN
           IF(AngmomB.LT.10)THEN
              WRITE(LUMOD3,'(A,I2,A,I1,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ELSEIF(AngmomB.LT.100)THEN
              WRITE(LUMOD3,'(A,I2,A,I1,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ENDIF
        ELSEIF(AngmomA.LT.100)THEN
           IF(AngmomB.LT.10)THEN
              WRITE(LUMOD3,'(A,I2,A,I2,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ELSEIF(AngmomB.LT.100)THEN
              WRITE(LUMOD3,'(A,I2,A,I2,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_LHS_P',AngmomP,'A',&
                   &AngmomA,'B',AngmomB,'(nContQP*nPasses,',nTUVQ,',Pdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
           ENDIF
        ENDIF
      ENDIF
      WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
      !swap 
      TMPSTRING = STRINGIN
      STRINGIN  = STRINGOUT
      STRINGOUT  = TMPSTRING
    ENDIF


    IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
       WRITE(LUMOD3,'(A)')'        !Spherical Transformation LHS'         
       IF(nlmA*nlmB*nTUVQ.LT.10)THEN
          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.100)THEN
          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.1000)THEN
          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVQ.LT.10000)THEN
          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVQ,'*nContQP*nPasses)'
       ENDIF
       IF(AngmomP.LT.10)THEN
          IF(AngmomA.LT.10)THEN
             WRITE(LUMOD3,'(A,I1,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomA.LT.100)THEN
             WRITE(LUMOD3,'(A,I1,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ELSEIF(AngmomP.LT.100)THEN
          IF(AngmomA.LT.10)THEN
             WRITE(LUMOD3,'(A,I2,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomA.LT.100)THEN
             WRITE(LUMOD3,'(A,I2,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS1_maxAngP',AngmomP,'_maxAngA',AngmomA,'(',nTUVQ,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ENDIF
       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ELSE
       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation LHS needed'
    ENDIF


    IF(AngmomQ.EQ.0)THEN
       stop 'this should never happenD'
    ELSE
       WRITE(LUMOD3,'(A)')'        !RHS Horizontal recurrence relations D to C'
       IF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.10)THEN
          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.100)THEN
          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.1000)THEN
          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nTUVCspec*nTUVDspec.LT.10000)THEN
          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nTUVCspec*nTUVDspec,'*nContQP*nPasses)'
       ENDIF
       IF(AngmomQ.LT.10)THEN
          IF(AngmomC.LT.10)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A,I1,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I1,A,I1,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ELSEIF(AngmomC.LT.100)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I1,A,I2,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I1,A,I2,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ENDIF
       ELSEIF(AngmomQ.LT.100)THEN
          IF(AngmomC.LT.10)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I2,A,I1,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I2,A,I1,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ELSEIF(AngmomC.LT.100)THEN
             IF(AngmomD.LT.10)THEN
                WRITE(LUMOD3,'(A,I2,A,I2,A,I1,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ELSEIF(AngmomD.LT.100)THEN
                WRITE(LUMOD3,'(A,I2,A,I2,A,I2,A,I4,A,A,A,A,A)')'        call HorizontalRR_RHS_Q',AngmomQ,'C',&
                     &AngmomC,'D',AngmomD,'(nContQP,nPasses,',nlmA*nlmB,',Qdistance12,',STRINGIN,',',STRINGOUT,',lupri)'
             ENDIF
          ENDIF
       ENDIF
       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF

    IF(Spherical.AND.(AngmomC.GT.1.OR.AngmomD.GT.1))THEN
       WRITE(LUMOD3,'(A)')'        !Spherical Transformation RHS'

       IF(nlmA*nlmB*nlmC*nlmD.LT.10)THEN
          WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.100)THEN
          WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.1000)THEN
          WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ELSEIF(nlmA*nlmB*nlmC*nlmD.LT.10000)THEN
          WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nlmA*nlmB*nlmC*nlmD,'*nContQP*nPasses)'
       ENDIF
       IF(AngmomQ.LT.10)THEN
          IF(AngmomC.LT.10)THEN
             WRITE(LUMOD3,'(A,I1,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomC.LT.100)THEN
             WRITE(LUMOD3,'(A,I1,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ELSEIF(AngmomQ.LT.100)THEN
          IF(AngmomC.LT.10)THEN
             WRITE(LUMOD3,'(A,I2,A,I1,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ELSEIF(AngmomC.LT.100)THEN
             WRITE(LUMOD3,'(A,I2,A,I2,A,I4,A,A,A,A,A)')'        call SphericalContractOBS2_maxAngQ',AngmomQ,'_maxAngC',AngmomC,'(',nlmA*nlmB,',nContQP*nPasses,',STRINGIN,',',STRINGOUT,')'
          ENDIF
       ENDIF
       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ELSE
       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation RHS needed'
    ENDIF

    WRITE(LUMOD3,'(A,A)')'        CDAB = ',STRINGIN
    WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
  end subroutine subroutineMainC

END PROGRAM TUV

!contractecoeff_gen

