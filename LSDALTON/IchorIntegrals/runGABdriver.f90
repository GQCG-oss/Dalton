PROGRAM TUV
  use math
  use stringsMODULE
  implicit none
  integer,pointer :: TUVINDEX(:,:,:),TUVINDEXP(:,:,:)
  integer :: JMAX,J,JMAX1,JMAXP
  logical,pointer :: Enoscreen(:,:),EnoscreenS(:,:),zero(:)
  integer :: ijk1,ijk2,ijkcart,ijk,ijkcart1,ijkcart2,nTUV,ijkP
  integer :: iTUV,ilmP
  real(realk),pointer :: SCMAT1(:,:),SCMAT2(:,:)
  logical :: sph1,sph2,sphericalTrans,sphericalGTO,newparam,ELSESTATEMENT,Spherical
  integer :: LUMAIN,LUMOD1,LUMOD2,LUMOD3,LUMOD4,iparam,iparam2,nparam
  integer :: LUMOD5,LUMOD6,LUMOD7
  integer :: JMAX2,JMAX3,JMAX4,ituvP,jp,tp,up,vp,ntuvP,l1,l2,l12
  integer :: ip1,jp1,kp1,p1,ip2,jp2,kp2,p2,ijkcartP,AngmomA,AngmomB
  integer :: AngmomC,AngmomD,AngmomP,AngmomQ,nTUVQ,AngmomPQ
  integer :: nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec
  integer :: nTUVA,nTUVB,nTUVC,nTUVD
  integer :: nlmA,nlmB,nlmC,nlmD,angmomID,iseg,ILUMOD,I,nUniquenTUVs
  real(realk),pointer :: uniqeparam(:)
  character(len=15),pointer :: uniqeparamNAME(:)
  character(len=9) :: STRINGIN,STRINGOUT,TMPSTRING
  character(len=4) :: SPEC
  logical :: BUILD(0:2,0:2),Gen,Seg,SegP,segQ,Seg1Prim,UNIQUE
  integer,pointer :: UniquenTUVs(:)
  !TODO
  !remove mem_alloc
  !remove CDAB = TMParray2
  !add PrimitiveContractionSeg to Transfer or Vertical 
  !
  LUMOD2=2
  open(unit = LUMOD2, file="GAB_OBS_DRIVER.f90",status="unknown")
  WRITE(LUMOD2,'(A)')'MODULE IchorEriGabintegralOBSGeneralMod'
  WRITE(LUMOD2,'(A)')'!Automatic Generated Code (AGC) by runGABdriver.f90 in tools directory'

  LUMOD3=3
  open(unit = LUMOD3, file="GAB_OBS_DRIVERGen.f90",status="unknown")
  WRITE(LUMOD3,'(A)')'MODULE IchorEriGabintegralOBSGeneralModGen'
  WRITE(LUMOD3,'(A)')'!Automatic Generated Code (AGC) by runGABdriver.f90 in tools directory'
  WRITE(LUMOD3,'(A)')'!Contains routines for General Contracted Basisset '

  LUMOD6=4
  open(unit = LUMOD6, file="GAB_OBS_DRIVERSeg.f90",status="unknown")
  WRITE(LUMOD6,'(A)')'MODULE IchorEriGabintegralOBSGeneralModSeg'
  WRITE(LUMOD6,'(A)')'!Automatic Generated Code (AGC) by runGABdriver.f90 in tools directory'
  WRITE(LUMOD6,'(A)')'!Contains routines for Segmented contracted Basisset '

  WRITE(LUMOD2,'(A)')'use IchorEriGabintegralOBSGeneralModGen'
  WRITE(LUMOD2,'(A)')'use IchorEriGabintegralOBSGeneralModSeg'
  DO ILUMOD=2,4
     WRITE(ILUMOD,'(A)')'use IchorprecisionModule'
     WRITE(ILUMOD,'(A)')'use IchorCommonModule'
     WRITE(ILUMOD,'(A)')'use IchorMemory'
  ENDDO
  !Since we need primitivecontractiongenXXX
  WRITE(LUMOD3,'(A)')'use IchorEriCoulombintegralOBSGeneralModGen'

  WRITE(LUMOD3,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODA'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODB'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODD'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODC'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODA'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODB'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODD'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODC'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODASeg'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODBSeg'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODDSeg'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_VERTICALRECURRENCEMODCSeg'

  WRITE(LUMOD3,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODAtoCGen'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODAtoDGen'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODBtoCGen'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODBtoDGen'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODCtoAGen'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODDtoAGen'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODCtoBGen'
  WRITE(LUMOD3,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODDtoBGen'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODAtoCSeg'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODAtoDSeg'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODBtoCSeg'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODBtoDSeg'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODCtoASeg'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODDtoASeg'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODCtoBSeg'
  WRITE(LUMOD6,'(A)')'use AGC_OBS_TRANSFERRECURRENCEMODDtoBSeg'
  DO ILUMOD=3,4
     WRITE(ILUMOD,'(A)')'use AGC_OBS_HorizontalRecurrenceLHSModAtoB'
     WRITE(ILUMOD,'(A)')'use AGC_OBS_HorizontalRecurrenceLHSModBtoA'
     WRITE(ILUMOD,'(A)')'use AGC_OBS_HorizontalRecurrenceRHSModCtoD'
     WRITE(ILUMOD,'(A)')'use AGC_OBS_HorizontalRecurrenceRHSModDtoC'
     WRITE(ILUMOD,'(A)')'use AGC_OBS_Sphcontract1Mod'
     WRITE(ILUMOD,'(A)')'use AGC_OBS_Sphcontract2Mod'
     WRITE(ILUMOD,'(A)')'  '
     WRITE(ILUMOD,'(A)')'private   '
  ENDDO
  WRITE(LUMOD2,'(A)')'public :: IchorGabIntegral_OBS_general,IchorGabIntegral_OBS_general_size  '
  WRITE(LUMOD3,'(A)')'public :: IchorGabIntegral_OBS_Gen,IchorGabIntegral_OBS_general_sizeGen  '
  WRITE(LUMOD6,'(A)')'public :: IchorGabIntegral_OBS_Seg,IchorGabIntegral_OBS_general_sizeSeg  '

  DO ILUMOD=2,4
     WRITE(ILUMOD,'(A)')'  '
     WRITE(ILUMOD,'(A)')'CONTAINS'
     WRITE(ILUMOD,'(A)')'  '
     WRITE(ILUMOD,'(A)')'  '
  ENDDO
  WRITE(LUMOD2,'(A)')'  subroutine IchorGabIntegral_OBS_general(nPrimA,nPrimB,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContP,pexp,ACC,BCC,&'
  WRITE(LUMOD2,'(A)')'       & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
  WRITE(LUMOD2,'(A)')'       & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&'
  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,Pdistance12,PQorder,CDAB,Acenter,Bcenter,&'
  WRITE(LUMOD2,'(A)')'       & spherical,TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)'
  WRITE(LUMOD2,'(A)')'    implicit none'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nPrimP,nPasses,nPrimA,nPrimB'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: MaxPasses,IntPrint,lupri'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nContA,nContB,nContP,nTABFJW1,nTABFJW2'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: AngmomA,AngmomB'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB)'
  WRITE(LUMOD2,'(A)')'    logical,intent(in)     :: Psegmented'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: pexp(nPrimP)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: PpreExpFac(nPrimP)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)'
  WRITE(LUMOD2,'(A)')'    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
  WRITE(LUMOD2,'(A)')'    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(inout) :: CDAB(nPasses)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: integralPrefactor(nPrimP*nPrimP)'
  WRITE(LUMOD2,'(A)')'    logical,intent(in) :: PQorder'
  WRITE(LUMOD2,'(A)')'    !integralPrefactor(nPrimP,nPrimP)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: reducedExponents(nPrimP*nPrimP)'
  WRITE(LUMOD2,'(A)')'    !reducedExponents(nPrimP,nPrimP)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter '
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: Acenter(3),Bcenter(3)'
  WRITE(LUMOD2,'(A)')'    logical,intent(in) :: spherical'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize'
  WRITE(LUMOD2,'(A)')'!   TMP variables - allocated outside'  
  WRITE(LUMOD2,'(A)')'    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)'
  WRITE(LUMOD2,'(A)')'    IF(PQorder)THEN'
  WRITE(LUMOD2,'(A)')'       call IchorQuit(''PQorder OBS general expect to get QP ordering'',-1)'
  WRITE(LUMOD2,'(A)')'    ENDIF'
  WRITE(LUMOD2,'(A)')'    IF(.NOT.spherical)THEN'
  WRITE(LUMOD2,'(A)')'       call IchorQuit(''cartesian not testet'',-1)'
  WRITE(LUMOD2,'(A)')'    ENDIF'
  WRITE(LUMOD2,'(A)')'    '
  WRITE(LUMOD2,'(A)')'   IF(Psegmented)THEN'
  WRITE(LUMOD2,'(A)')'    call IchorGabIntegral_OBS_Seg(nPrimA,nPrimB,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContP,pexp,ACC,BCC,&'
  WRITE(LUMOD2,'(A)')'       & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
  WRITE(LUMOD2,'(A)')'       & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&'
  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,Pdistance12,PQorder,CDAB,Acenter,Bcenter,&'
  WRITE(LUMOD2,'(A)')'       & spherical,TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)'
  WRITE(LUMOD2,'(A)')'   ELSE'
  WRITE(LUMOD2,'(A)')'    call IchorGabIntegral_OBS_Gen(nPrimA,nPrimB,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContP,pexp,ACC,BCC,&'
  WRITE(LUMOD2,'(A)')'       & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
  WRITE(LUMOD2,'(A)')'       & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&'
  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,Pdistance12,PQorder,CDAB,Acenter,Bcenter,&'
  WRITE(LUMOD2,'(A)')'       & spherical,TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)'
  WRITE(LUMOD2,'(A)')'   ENDIF'
  WRITE(LUMOD2,'(A)')'  end subroutine IchorGabIntegral_OBS_general'
  WRITE(LUMOD2,'(A)')'  '

  DO ISEG = 1,2
     Gen=.FALSE.; Seg=.FALSE.;
     IF(ISEG.EQ.1)THEN
        Gen=.TRUE.       
     ELSEIF(ISEG.EQ.2)THEN
        Seg=.TRUE.
     ENDIF

     IF(Gen)THEN
        WRITE(LUMOD3,'(A)')'  subroutine IchorGabIntegral_OBS_Gen(nPrimA,nPrimB,&'
        ILUMOD = 3
     ELSEIF(Seg)THEN
        WRITE(LUMOD6,'(A)')'  subroutine IchorGabIntegral_OBS_Seg(nPrimA,nPrimB,&'
        ILUMOD = LUMOD6
     ENDIF
     WRITE(ILUMOD,'(A)')'       & nPrimP,nPasses,MaxPasses,IntPrint,lupri,&'
     WRITE(ILUMOD,'(A)')'       & nContA,nContB,nContP,pexp,ACC,BCC,&'
     WRITE(ILUMOD,'(A)')'       & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
     WRITE(ILUMOD,'(A)')'       & Aexp,Bexp,Psegmented,reducedExponents,integralPrefactor,&'
     WRITE(ILUMOD,'(A)')'       & AngmomA,AngmomB,Pdistance12,PQorder,CDAB,Acenter,Bcenter,&'
     WRITE(ILUMOD,'(A)')'       & spherical,TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize)'
     WRITE(ILUMOD,'(A)')'    implicit none'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nPrimP,nPasses,nPrimA,nPrimB'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: MaxPasses,IntPrint,lupri'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nContA,nContB,nContP,nTABFJW1,nTABFJW2'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: AngmomA,AngmomB'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB)'
     WRITE(ILUMOD,'(A)')'    logical,intent(in)     :: Psegmented'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: pexp(nPrimP)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: PpreExpFac(nPrimP)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)'
     WRITE(ILUMOD,'(A)')'    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
     WRITE(ILUMOD,'(A)')'    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(inout) :: CDAB(nPasses)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: integralPrefactor(nPrimP*nPrimP)'
     WRITE(ILUMOD,'(A)')'    logical,intent(in) :: PQorder'
     WRITE(ILUMOD,'(A)')'    !integralPrefactor(nPrimP,nPrimP)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: reducedExponents(nPrimP*nPrimP)'
     WRITE(ILUMOD,'(A)')'    !reducedExponents(nPrimP,nPrimP)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter '
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: Acenter(3),Bcenter(3)'
     WRITE(ILUMOD,'(A)')'    logical,intent(in) :: spherical'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize'
     WRITE(ILUMOD,'(A)')'!   TMP variables - allocated outside'  
     WRITE(ILUMOD,'(A)')'    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)'

     WRITE(ILUMOD,'(A)')'!   Local variables '  
     !  WRITE(ILUMOD,'(A)')'    real(realk),pointer :: squaredDistance(:)'!,Rpq(:)'!,Rqc(:),Rpa(:)
     WRITE(ILUMOD,'(A)')'    integer :: AngmomP,I,J,nContQP,la,lb,lc,ld,nsize,angmomid'


     WRITE(ILUMOD,'(A)')'    '
     WRITE(ILUMOD,'(A)')'    !Setup combined Angmom info'
     WRITE(ILUMOD,'(A)')'    AngmomP = AngmomA+AngmomB'
     !  WRITE(ILUMOD,'(A)')'    !Build the Boys Functions for argument squaredDistance*reducedExponents'
     !  WRITE(ILUMOD,'(A)')'    !save in RJ000 ordering (AngmomPQ+1),nPrimQ,nPrimP,nPasses'
     !  WRITE(ILUMOD,'(A)')'    call mem_ichor_alloc(RJ000,(AngmomPQ+1)*nPasses*nPrimQP)'
     !  WRITE(ILUMOD,'(A)')'    call buildRJ000_general(nPasses,nPrimQ,nPrimP,nTABFJW1,nTABFJW2,reducedExponents,&'
     !  WRITE(ILUMOD,'(A)')'         & TABFJW,RJ000,AngmomPQ,Pcent,Qcent)'
     !  WRITE(ILUMOD,'(A)')'    IF (INTPRINT .GE. 10) THEN'
     !  WRITE(ILUMOD,'(A)')'     WRITE(lupri,*)''Output from W000'''
     !  WRITE(ILUMOD,'(A)')'     DO I=1,nPrimQ*nPrimP*nPasses'
     !  WRITE(ILUMOD,'(A)')'       DO J=0,AngmomPQ'
     !  WRITE(ILUMOD,'(A)')'          WRITE(LUPRI,''(2X,A6,I4,A1,I4,A2,ES16.8)'')''RJ000('',J,'','',I,'')='',RJ000(1+J+(I-1)*(AngmomPQ+1))'
     !  WRITE(ILUMOD,'(A)')'       ENDDO'
     !  WRITE(ILUMOD,'(A)')'     ENDDO'
     !  WRITE(ILUMOD,'(A)')'    END IF'
     WRITE(ILUMOD,'(A)')'!    nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6'
     WRITE(ILUMOD,'(A)')'!    nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6'
     WRITE(ILUMOD,'(A)')'!    nlmA = 2*AngmomA+1'
     WRITE(ILUMOD,'(A)')'!    nlmB = 2*AngmomB+1'
     WRITE(ILUMOD,'(A)')'    AngmomID = 10*AngmomA+AngmomB'
     WRITE(ILUMOD,'(A)')'    SELECT CASE(AngmomID)'

     DO AngmomA = 0,2
        DO AngmomB = 0,2
           BUILD(AngmomA,AngmomB) = .TRUE.
        ENDDO
     ENDDO
     !the order may be important and should be 
     DO AngmomA = 0,2
        DO AngmomB = 0,AngmomA
           BUILD(AngmomA,AngmomB) = .FALSE.
           !==========================00=====================================

           AngmomID = 10*AngmomA+AngmomB
           WRITE(ILUMOD,'(A,I4,A,I2,A,I2,A,I2,A,I2,A)')'    CASE(',AngmomID,')  !Angmom(A=',AngmomA,',B=',AngmomB,',C=',AngmomA,',D=',AngmomB,') combi'
           AngmomP = AngmomA + AngmomB
           nTUV = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
           nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
           nTUVQ = nTUVP
           nTUVAspec = (AngmomA+1)*(AngmomA+2)/2
           nTUVBspec = (AngmomB+1)*(AngmomB+2)/2
           nTUVCspec = nTUVAspec
           nTUVDspec = nTUVBspec
           nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
           nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6
           nTUVC = nTUVA
           nTUVD = nTUVB
           IF((AngmomA.GT.1.OR.AngmomB.GT.1))THEN
              spherical = .TRUE.
              nlmA = 2*AngmomA+1
              nlmB = 2*AngmomB+1
              nlmC = nlmA
              nlmD = nlmB
              STRINGIN(1:9)  = 'TMParray1'
              STRINGOUT(1:9) = 'TMParray2'
              TMPSTRING(1:9) = '         '
              !         WRITE(ILUMOD,'(A)')'      IF(spherical)THEN'
              call subroutineMAIN(ILUMOD,AngmomA,AngmomB,STRINGIN,STRINGOUT,TMPSTRING,AngmomP,&
                   & nTUVP,nTUVAspec,nTUVBspec,spherical,Gen,Seg)
           ELSE
              spherical = .TRUE.
              nlmA = 2*AngmomA+1
              nlmB = 2*AngmomB+1
              nlmC = nlmA
              nlmD = nlmB
              STRINGIN(1:9)  = 'TMParray1'
              STRINGOUT(1:9) = 'TMParray2'
              TMPSTRING(1:9) = '         '
              call subroutineMAIN(ILUMOD,AngmomA,AngmomB,STRINGIN,STRINGOUT,TMPSTRING,AngmomP,&
                   & nTUVP,nTUVAspec,nTUVBspec,spherical,Gen,Seg)
           ENDIF
        ENDDO
     ENDDO
     DO AngmomA = 0,2
        DO AngmomB = 0,2
           IF(BUILD(AngmomA,AngmomB))THEN
              !==========================00=====================================              
              AngmomID = 10*AngmomA+AngmomB
              WRITE(ILUMOD,'(A,I4,A,I2,A,I2,A,I2,A,I2,A)')'    CASE(',AngmomID,')  !Angmom(A=',AngmomA,',B=',AngmomB,',C=',AngmomA,',D=',AngmomB,') combi'
              AngmomP = AngmomA + AngmomB
              nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
              nTUVAspec = (AngmomA+1)*(AngmomA+2)/2
              nTUVBspec = (AngmomB+1)*(AngmomB+2)/2
              nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
              nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6
              IF((AngmomA.GT.1.OR.AngmomB.GT.1))THEN
                 spherical = .TRUE.
                 nlmA = 2*AngmomA+1
                 nlmB = 2*AngmomB+1
                 STRINGIN(1:9)  = 'TMParray1'
                 STRINGOUT(1:9) = 'TMParray2'
                 TMPSTRING(1:9) = '         '
                 !         WRITE(ILUMOD,'(A)')'      IF(spherical)THEN'
                 call subroutineMAIN(ILUMOD,AngmomA,AngmomB,STRINGIN,STRINGOUT,TMPSTRING,AngmomP,&
                      & nTUVP,nTUVAspec,nTUVBspec,spherical,Gen,Seg)
              ELSE
                 spherical = .TRUE.
                 nlmA = 2*AngmomA+1
                 nlmB = 2*AngmomB+1
                 STRINGIN(1:9)  = 'TMParray1'
                 STRINGOUT(1:9) = 'TMParray2'
                 TMPSTRING(1:9) = '         '
                 call subroutineMAIN(ILUMOD,AngmomA,AngmomB,STRINGIN,STRINGOUT,TMPSTRING,AngmomP,&
                      & nTUVP,nTUVAspec,nTUVBspec,spherical,Gen,Seg)
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     WRITE(ILUMOD,'(A)')'    CASE DEFAULT'

     IF(Gen)THEN
        WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in IchorGabIntegral_OBS_Gen'',-1)'
     ELSEIF(Seg)THEN
        WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in IchorGabIntegral_OBS_Seg'',-1)'
     ENDIF
     WRITE(ILUMOD,'(A)')'    END SELECT'
     IF(Gen)THEN
        WRITE(ILUMOD,'(A)')'  end subroutine IchorGabIntegral_OBS_Gen'
     ELSEIF(Seg)THEN
        WRITE(ILUMOD,'(A)')'  end subroutine IchorGabIntegral_OBS_Seg'
     ENDIF
     WRITE(ILUMOD,'(A)')'  '
  ENDDO


  WRITE(LUMOD2,'(A)')'  '
  WRITE(LUMOD2,'(A)')'  subroutine IchorGabIntegral_OBS_general_size(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         &TMParray2maxsize,AngmomA,AngmomB,nPrimP,nContP,Psegmented)'
  WRITE(LUMOD2,'(A)')'    implicit none'
  WRITE(LUMOD2,'(A)')'    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: AngmomA,AngmomB'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nPrimP,nContP'
  WRITE(LUMOD2,'(A)')'    logical,intent(in) :: Psegmented'
  WRITE(LUMOD2,'(A)')'    IF(Psegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call IchorGabIntegral_OBS_general_sizeSeg(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         &TMParray2maxsize,AngmomA,AngmomB,nPrimP,nContP)'
  WRITE(LUMOD2,'(A)')'    ELSE'
  WRITE(LUMOD2,'(A)')'     call IchorGabIntegral_OBS_general_sizeGen(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         &TMParray2maxsize,AngmomA,AngmomB,nPrimP,nContP)'
  WRITE(LUMOD2,'(A)')'    ENDIF'
  WRITE(LUMOD2,'(A)')'  end subroutine IchorGabIntegral_OBS_general_size'
  WRITE(LUMOD2,'(A)')'  '


  DO ISEG = 1,2
     Gen=.FALSE.; Seg=.FALSE.
     IF(ISEG.EQ.1)THEN
        Seg=.TRUE.
     ELSE
        Gen=.TRUE.       
     ENDIF

     WRITE(ISEG+2,'(A)')'  '
     IF(Gen)THEN
        WRITE(LUMOD3,'(A)')'  subroutine IchorGabIntegral_OBS_general_sizeGen(TMParray1maxsize,&'
        ILUMOD = LUMOD3
     ELSEIF(Seg)THEN
        WRITE(LUMOD6,'(A)')'  subroutine IchorGabIntegral_OBS_general_sizeSeg(TMParray1maxsize,&'
        ILUMOD = LUMOD6
     ENDIF
     WRITE(ILUMOD,'(A)')'         &TMParray2maxsize,AngmomA,AngmomB,nPrimP,nContP)'
     WRITE(ILUMOD,'(A)')'    implicit none'
     WRITE(ILUMOD,'(A)')'    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: AngmomA,AngmomB'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nPrimP,nContP'
     WRITE(ILUMOD,'(A)')'    ! local variables'
     WRITE(ILUMOD,'(A)')'    integer :: AngmomID'
     WRITE(ILUMOD,'(A)')'    '
     WRITE(ILUMOD,'(A)')'    AngmomID = 10*AngmomA+AngmomB'
     WRITE(ILUMOD,'(A)')'    TMParray2maxSize = 1'
     WRITE(ILUMOD,'(A)')'    TMParray1maxSize = 1'
     WRITE(ILUMOD,'(A)')'    SELECT CASE(AngmomID)'  
     DO AngmomA = 0,2
        DO AngmomB = 0,2
           AngmomID = 10*AngmomA+AngmomB
           WRITE(ILUMOD,'(A,I4,A,I2,A,I2,A,I2,A,I2,A)')'    CASE(',AngmomID,')  !Angmom(A=',AngmomA,',B=',AngmomB,',C=',AngmomA,',D=',AngmomB,') combi'

           AngmomP = AngmomA + AngmomB
           nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
           nTUVAspec = (AngmomA+1)*(AngmomA+2)/2
           nTUVBspec = (AngmomB+1)*(AngmomB+2)/2
           nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
           nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6           
           spherical = .TRUE.
           nlmA = 2*AngmomA+1
           nlmB = 2*AngmomB+1
           STRINGIN(1:9)  = 'TMParray1'
           STRINGOUT(1:9) = 'TMParray2'
           TMPSTRING(1:9) = '         '
           call determineSizes(ILUMOD,AngmomA,AngmomB,STRINGIN,STRINGOUT,TMPSTRING,AngmomP,&
                      & nTUVP,nTUVAspec,nTUVBspec,spherical,Gen,Seg)
        ENDDO
     ENDDO
     WRITE(ILUMOD,'(A)')'    CASE DEFAULT'
     WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in IchorGabIntegral_OBS_general_size'',-1)'
     WRITE(ILUMOD,'(A)')'    END SELECT'

     IF(Gen)THEN
        WRITE(LUMOD3,'(A)')'  end subroutine IchorGabIntegral_OBS_general_sizeGen'
        ILUMOD = LUMOD3
     ELSEIF(Seg)THEN
        WRITE(LUMOD6,'(A)')'  end subroutine IchorGabIntegral_OBS_general_sizeSeg'
        ILUMOD = LUMOD6
     ENDIF
  ENDDO

!!$  WRITE(LUMOD3,'(A)')''
!!$  WRITE(LUMOD3,'(A)')'  subroutine PrimitiveGabContractionGen1(AUXarrayCont,AUXarrayCont2,nContP,nPasses)'
!!$  WRITE(LUMOD3,'(A)')'    implicit none'
!!$  WRITE(LUMOD3,'(A)')'    !Only intrested in diagonal terms iContQ = iContP'
!!$  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPasses,nContP'
!!$  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarrayCont(nContP,nContP,nPasses)'
!!$  WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont2(nContP,nPasses)'
!!$  WRITE(LUMOD3,'(A)')'    !'
!!$  WRITE(LUMOD3,'(A)')'    integer :: iPassQ,iContP'
!!$  WRITE(LUMOD3,'(A)')'    do iPassQ = 1,nPasses'
!!$  WRITE(LUMOD3,'(A)')'     do iContP=1,nContP'
!!$  WRITE(LUMOD3,'(A)')'      AUXarrayCont2(iContP,iPassQ) = AUXarrayCont(iContP,iContP,iPassQ)'
!!$  WRITE(LUMOD3,'(A)')'     enddo'
!!$  WRITE(LUMOD3,'(A)')'    enddo'
!!$  WRITE(LUMOD3,'(A)')'  end subroutine PrimitiveGabContractionGen1'
!!$  WRITE(LUMOD3,'(A)')''
!!$  WRITE(LUMOD3,'(A)')''
  WRITE(LUMOD3,'(A)')'  subroutine GabPrimitiveContractionGen1(AUXarray2,AUXarrayCont,nPrimP,nPasses,&'
  WRITE(LUMOD3,'(A)')'       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)'
  WRITE(LUMOD3,'(A)')'    implicit none'
  WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPasses,nContP'
  WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray2(nPrimA,nPrimB,nPrimA,nPrimB,nPasses)'
  WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContA,nContB,nPasses)'
  WRITE(LUMOD3,'(A)')'    !'
  WRITE(LUMOD3,'(A)')'    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD'
  WRITE(LUMOD3,'(A)')'    real(realk) :: TMP,TMP1'
  WRITE(LUMOD3,'(A)')'    real(realk) :: tmpArray1(nPrimB,nPrimA,nPrimB)'
  WRITE(LUMOD3,'(A)')'    real(realk) :: tmpArray2(nPrimA,nPrimB)'
  WRITE(LUMOD3,'(A)')'    real(realk) :: tmpArray3(nPrimB)'
  WRITE(LUMOD3,'(A)')'    !Scaling p**4*c*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nPassQ'
  WRITE(LUMOD3,'(A)')'    do iPassQ = 1,nPasses'
  WRITE(LUMOD3,'(A)')'     do iContC=1,nContA'
  WRITE(LUMOD3,'(A)')'      do iPrimB=1,nPrimB'
  WRITE(LUMOD3,'(A)')'       do iPrimA=1,nPrimA'
  WRITE(LUMOD3,'(A)')'        do iPrimD=1,nPrimB'
  WRITE(LUMOD3,'(A)')'         TMP = 0.0E0_realk'
  WRITE(LUMOD3,'(A)')'         do iPrimC=1,nPrimA'
  WRITE(LUMOD3,'(A)')'          TMP = TMP + ACC(iPrimC,iContC)*AUXarray2(iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)'
  WRITE(LUMOD3,'(A)')'         enddo'
  WRITE(LUMOD3,'(A)')'         tmpArray1(iPrimD,iPrimA,iPrimB) = TMP'
  WRITE(LUMOD3,'(A)')'        enddo'
  WRITE(LUMOD3,'(A)')'       enddo'
  WRITE(LUMOD3,'(A)')'      enddo'
  WRITE(LUMOD3,'(A)')'      do iContD=1,nContB'
  WRITE(LUMOD3,'(A)')'       do iPrimB=1,nPrimB'
  WRITE(LUMOD3,'(A)')'        do iPrimA=1,nPrimA'
  WRITE(LUMOD3,'(A)')'         TMP = 0.0E0_realk'
  WRITE(LUMOD3,'(A)')'         do iPrimD=1,nPrimB'
  WRITE(LUMOD3,'(A)')'          TMP = TMP + BCC(iPrimD,iContD)*tmpArray1(iPrimD,iPrimA,iPrimB)'
  WRITE(LUMOD3,'(A)')'         enddo'
  WRITE(LUMOD3,'(A)')'         tmpArray2(iPrimA,iPrimB) = TMP'
  WRITE(LUMOD3,'(A)')'        enddo'
  WRITE(LUMOD3,'(A)')'       enddo'
  WRITE(LUMOD3,'(A)')'       TMP = 0.0E0_realk'
  WRITE(LUMOD3,'(A)')'       do iPrimB=1,nPrimB'
  WRITE(LUMOD3,'(A)')'        TMP1 = 0.0E0_realk'
  WRITE(LUMOD3,'(A)')'        do iPrimA=1,nPrimA'
  WRITE(LUMOD3,'(A)')'         TMP1 = TMP1 + ACC(iPrimA,iContC)*tmpArray2(iPrimA,iPrimB)'
  WRITE(LUMOD3,'(A)')'        enddo'
  WRITE(LUMOD3,'(A)')'        TMP = TMP + BCC(iPrimB,iContD)*TMP1'
  WRITE(LUMOD3,'(A)')'       enddo'
  WRITE(LUMOD3,'(A)')'       AUXarrayCont(iContC,iContD,iPassQ) = TMP'
  WRITE(LUMOD3,'(A)')'      enddo'
  WRITE(LUMOD3,'(A)')'     enddo'
  WRITE(LUMOD3,'(A)')'    enddo'
  WRITE(LUMOD3,'(A)')'  end subroutine GabPrimitiveContractionGen1'


  allocate(UniquenTUVs(3*3*3*3))
  UniquenTUVs(1) = 1
  nUniquenTUVs = 1
  DO AngmomA = 0,2
   DO AngmomB = 0,2
    AngmomP = AngmomA + AngmomB
    nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
    UNIQUE = .TRUE.
    DO I=1,nUniquenTUVs
       IF(nTUVP*nTUVP.EQ.UniquenTUVs(I))THEN
          UNIQUE = .FALSE.
       ENDIF
    ENDDO
    IF(UNIQUE)THEN
       nUniquenTUVs = nUniquenTUVs + 1
       UniquenTUVs(nUniquenTUVs) = nTUVP*nTUVP
!GENERAL 
         WRITE(LUMOD3,'(A)')''
         IF(nTUVP*nTUVP.LT.10)THEN
            WRITE(LUMOD3,'(A,I1,A)')'  subroutine GabPrimitiveContractionGen',nTUVP*nTUVP,'(AUXarray2,AUXarrayCont,nPrimP,nPasses,&'
         ELSEIF(nTUVP*nTUVP.LT.100)THEN
            WRITE(LUMOD3,'(A,I2,A)')'  subroutine GabPrimitiveContractionGen',nTUVP*nTUVP,'(AUXarray2,AUXarrayCont,nPrimP,nPasses,&'
         ELSEIF(nTUVP*nTUVP.LT.1000)THEN
            WRITE(LUMOD3,'(A,I3,A)')'  subroutine GabPrimitiveContractionGen',nTUVP*nTUVP,'(AUXarray2,AUXarrayCont,nPrimP,nPasses,&'
         ELSEIF(nTUVP*nTUVP.LT.10000)THEN
            WRITE(LUMOD3,'(A,I4,A)')'  subroutine GabPrimitiveContractionGen',nTUVP*nTUVP,'(AUXarray2,AUXarrayCont,nPrimP,nPasses,&'
         ELSEIF(nTUVP*nTUVP.LT.100000)THEN
            WRITE(LUMOD3,'(A,I5,A)')'  subroutine GabPrimitiveContractionGen',nTUVP*nTUVP,'(AUXarray2,AUXarrayCont,nPrimP,nPasses,&'
         ELSE
            STOP 'Primitive contraction'
         ENDIF
         WRITE(LUMOD3,'(A)')'       & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)'
         WRITE(LUMOD3,'(A)')'    implicit none'
         WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
         WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPasses,nContP'
         WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB'
         WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
         WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(in) :: AUXarray2(',nTUVP*nTUVP,',nPrimA,nPrimB,nPrimA,nPrimB,nPasses)'
         WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(inout) :: AUXarrayCont(',nTUVP*nTUVP,',nContA,nContB,nPasses)'
         WRITE(LUMOD3,'(A)')'    !'
         WRITE(LUMOD3,'(A)')'    integer :: iPassQ,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD'
         WRITE(LUMOD3,'(A)')'    integer :: iTUV'
         WRITE(LUMOD3,'(A,I5,A)')'    real(realk) :: TMPArray(',nTUVP*nTUVP,')'
         WRITE(LUMOD3,'(A,I5,A)')'    real(realk) :: tmpArray1(',nTUVP*nTUVP,',nPrimB,nPrimA,nPrimB)'
         WRITE(LUMOD3,'(A,I5,A)')'    real(realk) :: tmpArray2(',nTUVP*nTUVP,',nPrimA,nPrimB)'
         WRITE(LUMOD3,'(A,I5,A)')'    real(realk) :: tmpArray3(',nTUVP*nTUVP,')'
         WRITE(LUMOD3,'(A)')'    real(realk) :: ACCTMP,BCCTMP,CCCTMP,DCCTMP'
         WRITE(LUMOD3,'(A)')'    do iPassQ = 1,nPasses'
         WRITE(LUMOD3,'(A)')'     do iContC=1,nContA'
         WRITE(LUMOD3,'(A)')'      do iPrimB=1,nPrimB'
         WRITE(LUMOD3,'(A)')'       do iPrimA=1,nPrimA'
         WRITE(LUMOD3,'(A)')'        do iPrimD=1,nPrimB'
      WRITE(LUMOD3,'(A,I5)')'         do iTUV=1,',nTUVP*nTUVP
         WRITE(LUMOD3,'(A)')'          TMPArray(iTUV) = 0.0E0_realk'
         WRITE(LUMOD3,'(A)')'         enddo'
         WRITE(LUMOD3,'(A)')'         do iPrimC=1,nPrimA'
         WRITE(LUMOD3,'(A)')'          CCCTMP = ACC(iPrimC,iContC)'
         WRITE(LUMOD3,'(A)')'          !Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ'
      WRITE(LUMOD3,'(A,I5)')'          do iTUV=1,',nTUVP*nTUVP
         WRITE(LUMOD3,'(A)')'           TMPArray(iTUV) = TMPArray(iTUV) + CCCTMP*AUXarray2(iTUV,iPrimC,iPrimD,iPrimA,iPrimB,iPassQ)'
         WRITE(LUMOD3,'(A)')'          enddo'
         WRITE(LUMOD3,'(A)')'         enddo'
      WRITE(LUMOD3,'(A,I5)')'         do iTUV=1,',nTUVP*nTUVP
         WRITE(LUMOD3,'(A)')'          tmpArray1(iTUV,iPrimD,iPrimA,iPrimB) = TMPArray(iTUV)'
         WRITE(LUMOD3,'(A)')'         enddo'
         WRITE(LUMOD3,'(A)')'        enddo'
         WRITE(LUMOD3,'(A)')'       enddo'
         WRITE(LUMOD3,'(A)')'      enddo'
         WRITE(LUMOD3,'(A)')'      do iContD=1,nContB'
         WRITE(LUMOD3,'(A)')'       do iPrimB=1,nPrimB'
         WRITE(LUMOD3,'(A)')'        do iPrimA=1,nPrimA'
      WRITE(LUMOD3,'(A,I5)')'         do iTUV=1,',nTUVP*nTUVP
         WRITE(LUMOD3,'(A)')'          TMPArray(iTUV) = 0.0E0_realk'
         WRITE(LUMOD3,'(A)')'         enddo'
         WRITE(LUMOD3,'(A)')'         do iPrimD=1,nPrimB'
         WRITE(LUMOD3,'(A)')'          DCCTMP = BCC(iPrimD,iContD)'
         WRITE(LUMOD3,'(A)')'          !Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ'
      WRITE(LUMOD3,'(A,I5)')'          do iTUV=1,',nTUVP*nTUVP
         WRITE(LUMOD3,'(A)')'           TMPArray(iTUV) = TMPArray(iTUV) + DCCTMP*tmpArray1(iTUV,iPrimD,iPrimA,iPrimB)'
         WRITE(LUMOD3,'(A)')'          enddo'
         WRITE(LUMOD3,'(A)')'         enddo'
      WRITE(LUMOD3,'(A,I5)')'         do iTUV=1,',nTUVP*nTUVP
         WRITE(LUMOD3,'(A)')'          tmpArray2(iTUV,iPrimA,iPrimB) = TMPArray(iTUV)'
         WRITE(LUMOD3,'(A)')'         enddo'
         WRITE(LUMOD3,'(A)')'        enddo'
         WRITE(LUMOD3,'(A)')'       enddo'
      WRITE(LUMOD3,'(A,I5)')'       do iTUV=1,',nTUVP*nTUVP
         WRITE(LUMOD3,'(A)')'        TMPArray3(iTUV) = 0.0E0_realk'
         WRITE(LUMOD3,'(A)')'       enddo'
         WRITE(LUMOD3,'(A)')'       do iPrimB=1,nPrimB'
      WRITE(LUMOD3,'(A,I5)')'        do iTUV=1,',nTUVP*nTUVP
         WRITE(LUMOD3,'(A)')'         TMPArray(iTUV) = 0.0E0_realk'
         WRITE(LUMOD3,'(A)')'        enddo'
         WRITE(LUMOD3,'(A)')'        do iPrimA=1,nPrimA'
         WRITE(LUMOD3,'(A)')'         ACCTMP = ACC(iPrimA,iContC)'
      WRITE(LUMOD3,'(A,I5)')'         do iTUV=1,',nTUVP*nTUVP
         WRITE(LUMOD3,'(A)')'          TMPArray(iTUV) = TMPArray(iTUV) + ACCTMP*tmpArray2(iTUV,iPrimA,iPrimB)'
         WRITE(LUMOD3,'(A)')'         enddo'
         WRITE(LUMOD3,'(A)')'        enddo'
         WRITE(LUMOD3,'(A)')'        BCCTMP = BCC(iPrimB,iContD)'
      WRITE(LUMOD3,'(A,I5)')'        do iTUV=1,',nTUVP*nTUVP
         WRITE(LUMOD3,'(A)')'         TMPArray3(iTUV) = TMPArray3(iTUV) + BCCTMP*tmpArray(iTUV)'
         WRITE(LUMOD3,'(A)')'        enddo'
         WRITE(LUMOD3,'(A)')'       enddo'
      WRITE(LUMOD3,'(A,I5)')'       do iTUV=1,',nTUVP*nTUVP
         WRITE(LUMOD3,'(A)')'        AUXarrayCont(iTUV,iContC,iContD,iPassQ) = TMPArray3(iTUV)'
         WRITE(LUMOD3,'(A)')'       enddo'
         WRITE(LUMOD3,'(A)')'      enddo'
         WRITE(LUMOD3,'(A)')'     enddo'
         WRITE(LUMOD3,'(A)')'    enddo'
         IF(nTUVP*nTUVP.LT.10)THEN
            WRITE(LUMOD3,'(A,I1)')'  end subroutine GabPrimitiveContractionGen',nTUVP*nTUVP
         ELSEIF(nTUVP*nTUVP.LT.100)THEN
            WRITE(LUMOD3,'(A,I2)')'  end subroutine GabPrimitiveContractionGen',nTUVP*nTUVP
         ELSEIF(nTUVP*nTUVP.LT.1000)THEN
            WRITE(LUMOD3,'(A,I3)')'  end subroutine GabPrimitiveContractionGen',nTUVP*nTUVP
         ELSEIF(nTUVP*nTUVP.LT.10000)THEN
            WRITE(LUMOD3,'(A,I4)')'  end subroutine GabPrimitiveContractionGen',nTUVP*nTUVP
         ELSEIF(nTUVP*nTUVP.LT.100000)THEN
            WRITE(LUMOD3,'(A,I5)')'  end subroutine GabPrimitiveContractionGen',nTUVP*nTUVP
         ELSE
            STOP 'Primitive contraction'
         ENDIF

!!$
!!$       WRITE(LUMOD3,'(A)')''
!!$       IF(nTUVP*nTUVP.LT.10)THEN
!!$          WRITE(LUMOD3,'(A,I1,A)')'  subroutine PrimitiveGabContractionGen',nTUVP*nTUVP,'(AUXarrayCont,AUXarrayCont2,nContP,nPasses)'
!!$       ELSEIF(nTUVP*nTUVP.LT.100)THEN
!!$          WRITE(LUMOD3,'(A,I2,A)')'  subroutine PrimitiveGabContractionGen',nTUVP*nTUVP,'(AUXarrayCont,AUXarrayCont2,nContP,nPasses)'
!!$       ELSEIF(nTUVP*nTUVP.LT.1000)THEN
!!$          WRITE(LUMOD3,'(A,I3,A)')'  subroutine PrimitiveGabContractionGen',nTUVP*nTUVP,'(AUXarrayCont,AUXarrayCont2,nContP,nPasses)'
!!$       ELSEIF(nTUVP*nTUVP.LT.10000)THEN
!!$          WRITE(LUMOD3,'(A,I4,A)')'  subroutine PrimitiveGabContractionGen',nTUVP*nTUVP,'(AUXarrayCont,AUXarrayCont2,nContP,nPasses)'
!!$       ELSEIF(nTUVP*nTUVP.LT.100000)THEN
!!$          WRITE(LUMOD3,'(A,I5,A)')'  subroutine PrimitiveGabContractionGen',nTUVP*nTUVP,'(AUXarrayCont,AUXarrayCont2,nContP,nPasses)'
!!$       ELSE
!!$          STOP 'Primitive contraction'
!!$       ENDIF
!!$       WRITE(LUMOD3,'(A)')'    implicit none'
!!$       WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
!!$       WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPasses,nContP'
!!$       WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(in) :: AUXarrayCont(',nTUVP*nTUVP,',nContP,nContP,nPasses)'
!!$       WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(inout) :: AUXarrayCont2(',nTUVP*nTUVP,',nContP,nPasses)'
!!$       WRITE(LUMOD3,'(A)')'    !'
!!$       WRITE(LUMOD3,'(A)')'    integer :: iPassQ,iContP,iTUV'
!!$       WRITE(LUMOD3,'(A)')'    do iPassQ = 1,nPasses'
!!$       WRITE(LUMOD3,'(A)')'     do iContP=1,nContP'
!!$       WRITE(LUMOD3,'(A,I5)')'      do iTUV=1,',nTUVP*nTUVP
!!$       WRITE(LUMOD3,'(A)')'       AUXarrayCont2(iTUV,iContP,iPassQ) = AUXarrayCont(iTUV,iContP,iContP,iPassQ)'
!!$       WRITE(LUMOD3,'(A)')'      enddo'
!!$       WRITE(LUMOD3,'(A)')'     enddo'
!!$       WRITE(LUMOD3,'(A)')'    enddo'
!!$       IF(nTUVP*nTUVP.LT.10)THEN
!!$          WRITE(LUMOD3,'(A,I1)')'  end subroutine PrimitiveGabContractionGen',nTUVP*nTUVP
!!$       ELSEIF(nTUVP*nTUVP.LT.100)THEN
!!$          WRITE(LUMOD3,'(A,I2)')'  end subroutine PrimitiveGabContractionGen',nTUVP*nTUVP
!!$       ELSEIF(nTUVP*nTUVP.LT.1000)THEN
!!$          WRITE(LUMOD3,'(A,I3)')'  end subroutine PrimitiveGabContractionGen',nTUVP*nTUVP
!!$       ELSEIF(nTUVP*nTUVP.LT.10000)THEN
!!$          WRITE(LUMOD3,'(A,I4)')'  end subroutine PrimitiveGabContractionGen',nTUVP*nTUVP
!!$       ELSEIF(nTUVP*nTUVP.LT.100000)THEN
!!$          WRITE(LUMOD3,'(A,I5)')'  end subroutine PrimitiveGabContractionGen',nTUVP*nTUVP
!!$       ELSE
!!$          STOP 'Primitive contraction'
!!$       ENDIF
    ENDIF
 ENDDO
ENDDO
deallocate(UniquenTUVs)


WRITE(LUMOD3,'(A)')''
WRITE(LUMOD3,'(A)')'  subroutine ExtractGabElmP1Gen(AUXarray,Output,nContP,nPasses)'
WRITE(LUMOD3,'(A)')'    implicit none'
WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nContP,nPasses'
WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray(nContP,nPasses)'
WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: Output(nPasses)'
WRITE(LUMOD3,'(A)')'    !'
WRITE(LUMOD3,'(A)')'    integer :: iPassQ,iContP'
WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'
WRITE(LUMOD3,'(A)')'    do iPassQ = 1,nPasses'
WRITE(LUMOD3,'(A)')'     TMP = ABS(AUXarray(1,iPassQ))'
WRITE(LUMOD3,'(A)')'     do iContP=2,nContP'
WRITE(LUMOD3,'(A)')'      IF(ABS(AUXarray(iContP,iPassQ)).GT.TMP)THEN'
WRITE(LUMOD3,'(A)')'       TMP = ABS(AUXarray(iContP,iPassQ))'
WRITE(LUMOD3,'(A)')'      ENDIF'
WRITE(LUMOD3,'(A)')'     enddo'
WRITE(LUMOD3,'(A)')'     Output(iPassQ) = SQRT(TMP)'
WRITE(LUMOD3,'(A)')'    enddo'
WRITE(LUMOD3,'(A)')'  end subroutine ExtractGabElmP1Gen'
WRITE(LUMOD3,'(A)')''
WRITE(LUMOD6,'(A)')'  subroutine ExtractGabElmP1Seg(AUXarray,Output,nPasses)'
WRITE(LUMOD6,'(A)')'    implicit none'
WRITE(LUMOD6,'(A)')'    integer,intent(in) :: nPasses'
WRITE(LUMOD6,'(A)')'    real(realk),intent(in) :: AUXarray(nPasses)'
WRITE(LUMOD6,'(A)')'    real(realk),intent(inout) :: Output(nPasses)'
WRITE(LUMOD6,'(A)')'    !'
WRITE(LUMOD6,'(A)')'    integer :: iPassQ'
WRITE(LUMOD6,'(A)')'    do iPassQ = 1,nPasses'
WRITE(LUMOD6,'(A)')'     Output(iPassQ) = SQRT(ABS(AUXarray(iPassQ)))'
WRITE(LUMOD6,'(A)')'    enddo'
WRITE(LUMOD6,'(A)')'  end subroutine ExtractGabElmP1Seg'

allocate(UniquenTUVs(3*3*3*3))
UniquenTUVs(1) = 1
nUniquenTUVs = 1
DO AngmomA = 0,2
   DO AngmomB = 0,2
      AngmomP = AngmomA + AngmomB
      nlmA = 2*AngmomA+1
      nlmB = 2*AngmomB+1
      UNIQUE = .TRUE.
      DO I=1,nUniquenTUVs
         IF(nlmA*nlmB.EQ.UniquenTUVs(I))THEN
            UNIQUE = .FALSE.
         ENDIF
      ENDDO
      IF(UNIQUE)THEN
         nUniquenTUVs = nUniquenTUVs + 1
         UniquenTUVs(nUniquenTUVs) = nlmA*nlmB
         !GENERAL 
         WRITE(LUMOD3,'(A)')''
         IF(nlmA*nlmB.LT.10)THEN
            WRITE(LUMOD3,'(A,I1,A)')'  subroutine ExtractGabElmP',nlmA*nlmB,'Gen(AUXarray,Output,nContP,nPasses)'
         ELSEIF(nlmA*nlmB.LT.100)THEN
            WRITE(LUMOD3,'(A,I2,A)')'  subroutine ExtractGabElmP',nlmA*nlmB,'Gen(AUXarray,Output,nContP,nPasses)'
         ELSEIF(nlmA*nlmB.LT.1000)THEN
            WRITE(LUMOD3,'(A,I3,A)')'  subroutine ExtractGabElmP',nlmA*nlmB,'Gen(AUXarray,Output,nContP,nPasses)'
         ELSEIF(nlmA*nlmB.LT.10000)THEN
            WRITE(LUMOD3,'(A,I4,A)')'  subroutine ExtractGabElmP',nlmA*nlmB,'Gen(AUXarray,Output,nContP,nPasses)'
         ELSEIF(nlmA*nlmB.LT.100000)THEN
            WRITE(LUMOD3,'(A,I5,A)')'  subroutine ExtractGabElmP',nlmA*nlmB,'Gen(AUXarray,Output,nContP,nPasses)'
         ELSE
            STOP 'Primitive contraction'
         ENDIF
         WRITE(LUMOD3,'(A)')'    implicit none'
         WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nContP,nPasses'
         WRITE(LUMOD3,'(A,I5,A,I5,A)')'    real(realk),intent(in) :: AUXarray(',nlmA*nlmB,',',nlmA*nlmB,',nContP,nPasses)'
         WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: Output(nPasses)'
         WRITE(LUMOD3,'(A)')'    !'
         WRITE(LUMOD3,'(A)')'    integer :: iPassQ,iContP,i'
         WRITE(LUMOD3,'(A,I5,A)')'    real(realk) :: TMP(',nlmA*nlmB,')'
         WRITE(LUMOD3,'(A)')'    real(realk) :: MaxValue,TotalMaxValue'
         WRITE(LUMOD3,'(A)')'    do iPassQ = 1,nPasses'
         WRITE(LUMOD3,'(A,I5)')'     do i=1,',nlmA*nlmB
         WRITE(LUMOD3,'(A)')'      TMP(i) = ABS(AUXarray(i,i,1,iPassQ))'
         WRITE(LUMOD3,'(A)')'     enddo'
         WRITE(LUMOD3,'(A)')'     TotalMaxValue = MAXVAL(TMP)'
         WRITE(LUMOD3,'(A)')'     do iContP=2,nContP'
         WRITE(LUMOD3,'(A,I5)')'      do i=1,',nlmA*nlmB
         WRITE(LUMOD3,'(A)')'       TMP(i) = ABS(AUXarray(i,i,iContP,iPassQ))'
         WRITE(LUMOD3,'(A)')'      enddo'
         WRITE(LUMOD3,'(A)')'      maxvalue = MAXVAL(TMP)'
         WRITE(LUMOD3,'(A)')'      IF(MaxValue.GT.TotalMaxValue)THEN'
         WRITE(LUMOD3,'(A)')'       TotalMaxValue = MaxValue'       
         WRITE(LUMOD3,'(A)')'      ENDIF'
         WRITE(LUMOD3,'(A)')'     enddo'
         WRITE(LUMOD3,'(A)')'     Output(iPassQ) = SQRT(TotalMaxValue)'
         WRITE(LUMOD3,'(A)')'    enddo'
         IF(nlmA*nlmB.LT.10)THEN
            WRITE(LUMOD3,'(A,I1,A)')'  end subroutine ExtractGabElmP',nlmA*nlmB,'Gen'
         ELSEIF(nlmA*nlmB.LT.100)THEN
            WRITE(LUMOD3,'(A,I2,A)')'  end subroutine ExtractGabElmP',nlmA*nlmB,'Gen'
         ELSEIF(nlmA*nlmB.LT.1000)THEN
            WRITE(LUMOD3,'(A,I3,A)')'  end subroutine ExtractGabElmP',nlmA*nlmB,'Gen'
         ELSEIF(nlmA*nlmB.LT.10000)THEN
            WRITE(LUMOD3,'(A,I4,A)')'  end subroutine ExtractGabElmP',nlmA*nlmB,'Gen'
         ELSEIF(nlmA*nlmB.LT.100000)THEN
            WRITE(LUMOD3,'(A,I5,A)')'  end subroutine ExtractGabElmP',nlmA*nlmB,'Gen'
         ELSE
            STOP 'Primitive contraction'
         ENDIF

         !Segmented
         WRITE(LUMOD6,'(A)')''
         IF(nlmA*nlmB.LT.10)THEN
            WRITE(LUMOD6,'(A,I1,A)')'  subroutine ExtractGabElmP',nlmA*nlmB,'Seg(AUXarray,Output,nPasses)'
         ELSEIF(nlmA*nlmB.LT.100)THEN
            WRITE(LUMOD6,'(A,I2,A)')'  subroutine ExtractGabElmP',nlmA*nlmB,'Seg(AUXarray,Output,nPasses)'
         ELSEIF(nlmA*nlmB.LT.1000)THEN
            WRITE(LUMOD6,'(A,I3,A)')'  subroutine ExtractGabElmP',nlmA*nlmB,'Seg(AUXarray,Output,nPasses)'
         ELSEIF(nlmA*nlmB.LT.10000)THEN
            WRITE(LUMOD6,'(A,I4,A)')'  subroutine ExtractGabElmP',nlmA*nlmB,'Seg(AUXarray,Output,nPasses)'
         ELSEIF(nlmA*nlmB.LT.100000)THEN
            WRITE(LUMOD6,'(A,I5,A)')'  subroutine ExtractGabElmP',nlmA*nlmB,'Seg(AUXarray,Output,nPasses)'
         ELSE
            STOP 'Primitive contraction'
         ENDIF
         WRITE(LUMOD6,'(A)')'    implicit none'
         WRITE(LUMOD6,'(A)')'    integer,intent(in) :: nPasses'
         WRITE(LUMOD6,'(A,I5,A,I5,A)')'    real(realk),intent(in) :: AUXarray(',nlmA*nlmB,',',nlmA*nlmB,',nPasses)'
         WRITE(LUMOD6,'(A)')'    real(realk),intent(inout) :: Output(nPasses)'
         WRITE(LUMOD6,'(A)')'    !'
         WRITE(LUMOD6,'(A)')'    integer :: iPassQ,i'
         WRITE(LUMOD6,'(A,I5,A)')'    real(realk) :: TMP(',nlmA*nlmB,')'
         WRITE(LUMOD6,'(A)')'    real(realk) :: MaxValue,TotalMaxValue'
         WRITE(LUMOD6,'(A)')'    do iPassQ = 1,nPasses'
      WRITE(LUMOD6,'(A,I5)')'     do i=1,',nlmA*nlmB
         WRITE(LUMOD6,'(A)')'      TMP(i) = ABS(AUXarray(i,i,iPassQ))'
         WRITE(LUMOD6,'(A)')'     enddo'
         WRITE(LUMOD6,'(A)')'     Output(iPassQ) = SQRT(MAXVAL(TMP))'
         WRITE(LUMOD6,'(A)')'    enddo'
         IF(nlmA*nlmB.LT.10)THEN
            WRITE(LUMOD6,'(A,I1,A)')'  end subroutine ExtractGabElmP',nlmA*nlmB,'Seg'
         ELSEIF(nlmA*nlmB.LT.100)THEN
            WRITE(LUMOD6,'(A,I2,A)')'  end subroutine ExtractGabElmP',nlmA*nlmB,'Seg'
         ELSEIF(nlmA*nlmB.LT.1000)THEN
            WRITE(LUMOD6,'(A,I3,A)')'  end subroutine ExtractGabElmP',nlmA*nlmB,'Seg'
         ELSEIF(nlmA*nlmB.LT.10000)THEN
            WRITE(LUMOD6,'(A,I4,A)')'  end subroutine ExtractGabElmP',nlmA*nlmB,'Seg'
         ELSEIF(nlmA*nlmB.LT.100000)THEN
            WRITE(LUMOD6,'(A,I5,A)')'  end subroutine ExtractGabElmP',nlmA*nlmB,'Seg'
         ELSE
            STOP 'Primitive contraction'
         ENDIF
      ENDIF
   ENDDO
ENDDO
deallocate(UniquenTUVs)

WRITE(LUMOD2,'(A)')'END MODULE IchorEriGabintegralOBSGeneralMod'
WRITE(LUMOD3,'(A)')'END MODULE IchorEriGabintegralOBSGeneralModGen'
WRITE(LUMOD6,'(A)')'END MODULE IchorEriGabintegralOBSGeneralModSeg'

close(unit = LUMOD2)
close(unit = LUMOD3)
close(unit = LUMOD6)

contains
  subroutine subroutineMain(LUMOD3,AngmomA,AngmomB,STRINGIN,STRINGOUT,TMPSTRING,AngmomP,&
       & nTUVP,nTUVAspec,nTUVBspec,spherical,Gen,Seg)
    implicit none
    integer,intent(in) :: LUMOD3,AngmomA,AngmomB,AngmomP
    integer,intent(in) :: nTUVP,nTUVAspec,nTUVBspec
    character(len=9) :: STRINGIN,STRINGOUT,TMPSTRING
    logical :: spherical,OutputSet,Gen,Seg,Contracted
    character(len=8) :: BASISSPEC
    integer :: iBasisSpec
    AngmomPQ = AngmomP + AngmomP
    IF(Gen)THEN
       iBasisSpec = 3
       BASISSPEC = 'Gen     '
    ELSEIF(Seg)THEN
       iBasisSpec = 3
       BASISSPEC = 'Seg     '
    ENDIF
    Contracted = .FALSE.
    IF(AngmomP.EQ.0)THEN
       IF(Gen)THEN
          WRITE(LUMOD3,'(A)')'        call VerticalRecurrence0(nPasses,nPrimP,nPrimP,&'
       ELSEIF(Seg)THEN
          WRITE(LUMOD3,'(A)')'        call VerticalRecurrenceSeg0(nPasses,nPrimP,nPrimP,&'
          Contracted = .TRUE.
       ENDIF
       WRITE(LUMOD3,'(A)')'               & reducedExponents,TABFJW,Pcent,Pcent,integralPrefactor,&'
       WRITE(LUMOD3,'(A,A,A)')'               & PpreExpFac,PpreExpFac,',STRINGOUT,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
       WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
    ELSE
       IF(AngmomA.GE.AngmomB)THEN
          !A Vertical recurrence
          IF(AngmomPQ.LT.10)THEN
             WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'A(nPasses,nPrimP,nPrimP,reducedExponents,&'
          ELSE
             WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'A(nPasses,nPrimP,nPrimP,reducedExponents,&'
          ENDIF
          WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Pexp,Acenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,',STRINGOUT,')'
          !swap 
          TMPSTRING = STRINGIN
          STRINGIN  = STRINGOUT
          STRINGOUT  = TMPSTRING
          !determine TransferRecurrence
          !A to C TransferRecurrence
          SPEC = 'AtoC'
          call initString(8)
          call AddToString('call TransferRecurrenceP')
          call AddToString(AngmomP)
          call AddToString('Q')
          call AddToString(AngmomP)
          call AddToString(SPEC)                
          call AddToString(BASISSPEC(1:iBasisSpec))
          call AddToString('(nPasses,nPrimP,nPrimP,reducedExponents,&')
          call writeString(LUMOD3)
          IF(.NOT.Gen)Contracted = .TRUE.
          WRITE(LUMOD3,'(A)')'               & Pexp,Pexp,Pdistance12,Pdistance12,Bexp,Bexp,nPrimA,nPrimB,nPrimA,nPrimB,&'                
          WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
          !swap 
          TMPSTRING = STRINGIN
          STRINGIN  = STRINGOUT
          STRINGOUT  = TMPSTRING
       ELSE
          !B Vertical recurrence
          IF(AngmomPQ.LT.10)THEN
             WRITE(LUMOD3,'(A,I1,A)')'        call VerticalRecurrence',AngmomPQ,'B(nPasses,nPrimP,nPrimP,reducedExponents,&'
          ELSE
             WRITE(LUMOD3,'(A,I2,A)')'        call VerticalRecurrence',AngmomPQ,'B(nPasses,nPrimP,nPrimP,reducedExponents,&'
          ENDIF
          WRITE(LUMOD3,'(A,A,A)')'               & TABFJW,Pexp,Bcenter,Pcent,Pcent,integralPrefactor,PpreExpFac,PpreExpFac,',STRINGOUT,')'
          !swap 
          TMPSTRING = STRINGIN
          STRINGIN  = STRINGOUT
          STRINGOUT  = TMPSTRING
          !determine TransferRecurrence
          !B to D TransferRecurrence
          SPEC = 'BtoD'
          call initString(8)
          call AddToString('call TransferRecurrenceP')
          call AddToString(AngmomP)
          call AddToString('Q')
          call AddToString(AngmomP)
          call AddToString(SPEC)                
          call AddToString(BASISSPEC(1:iBasisSpec))
          call AddToString('(nPasses,nPrimP,nPrimP,reducedExponents,&')
          call writeString(LUMOD3)
          IF(.NOT.Gen)Contracted = .TRUE.
          !B to D TransferRecurrence
          WRITE(LUMOD3,'(A)')'               & Pexp,Pexp,Pdistance12,Pdistance12,Aexp,Aexp,nPrimA,nPrimB,nPrimA,nPrimB,&'                
          WRITE(LUMOD3,'(A,A,A,A,A)')'               & ',STRINGIN,',',STRINGOUT,')'
          !swap 
          TMPSTRING = STRINGIN
          STRINGIN  = STRINGOUT
          STRINGOUT  = TMPSTRING
       ENDIF
    ENDIF
    !================= DONE WITH VERTICAL AND TRANSFER ================================================================
!    IF(Gen)THEN
!       WRITE(LUMOD3,'(A)')'        nContQP = nContP*nContP'
!    ENDIF
    IF(Contracted)THEN
       WRITE(LUMOD3,'(A)')'        !Primitive Contraction have already been done'
    ELSE
       IF(nTUVP*nTUVP.LT.10)THEN       
          WRITE(LUMOD3,'(A,I1,A,A,A,A,A)')'         call GabPrimitiveContractionGen',nTUVP*nTUVP,'(',STRINGIN,',',STRINGOUT,',nPrimP,nPasses,&'
       ELSEIF(nTUVP*nTUVP.LT.100)THEN
          WRITE(LUMOD3,'(A,I2,A,A,A,A,A)')'         call GabPrimitiveContractionGen',nTUVP*nTUVP,'(',STRINGIN,',',STRINGOUT,',nPrimP,nPasses,&'
       ELSEIF(nTUVP*nTUVP.LT.1000)THEN
          WRITE(LUMOD3,'(A,I3,A,A,A,A,A)')'         call GabPrimitiveContractionGen',nTUVP*nTUVP,'(',STRINGIN,',',STRINGOUT,',nPrimP,nPasses,&'
       ELSEIF(nTUVP*nTUVP.LT.10000)THEN
          WRITE(LUMOD3,'(A,I4,A,A,A,A,A)')'         call GabPrimitiveContractionGen',nTUVP*nTUVP,'(',STRINGIN,',',STRINGOUT,',nPrimP,nPasses,&'
       ELSE
          STOP 'GenPrimCont'
       ENDIF
       WRITE(LUMOD3,'(A)')'              & nContP,ACC,BCC,nPrimA,nContA,nPrimB,nContB)'
       !swap 
!!$       TMPSTRING = STRINGIN
!!$       STRINGIN  = STRINGOUT
!!$       STRINGOUT  = TMPSTRING
!!$       IF(nTUVP*nTUVP.LT.10)THEN       
!!$          WRITE(LUMOD3,'(A,I1,A,A,A,A,A)')'         call PrimitiveGabContractionGen',nTUVP*nTUVP,'(',STRINGIN,',',STRINGOUT,',nContP,nPasses)'
!!$       ELSEIF(nTUVP*nTUVP.LT.100)THEN
!!$          WRITE(LUMOD3,'(A,I2,A,A,A,A,A)')'         call PrimitiveGabContractionGen',nTUVP*nTUVP,'(',STRINGIN,',',STRINGOUT,',nContP,nPasses)'
!!$       ELSEIF(nTUVP*nTUVP.LT.1000)THEN
!!$          WRITE(LUMOD3,'(A,I3,A,A,A,A,A)')'         call PrimitiveGabContractionGen',nTUVP*nTUVP,'(',STRINGIN,',',STRINGOUT,',nContP,nPasses)'
!!$       ELSEIF(nTUVP*nTUVP.LT.10000)THEN
!!$          WRITE(LUMOD3,'(A,I4,A,A,A,A,A)')'         call PrimitiveGabContractionGen',nTUVP*nTUVP,'(',STRINGIN,',',STRINGOUT,',nContP,nPasses)'
!!$       ELSE
!!$          STOP 'GenPrimCont'
!!$       ENDIF
!!$       WRITE(LUMOD3,'(A)')'        nContQP = nContP ! The diagonal elements have been extracted'       
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF
    !     The horizontal recurrence also extract nTUVspec from nTUV       
    IF(AngmomP.EQ.0)THEN
       WRITE(LUMOD3,'(A)')'        !no need for LHS Horizontal recurrence relations, it would be a simply copy'
       !there will not be need for spherical transformation afterwards
       !there will not be need for RHS Horizontal recurrence relations nor Spherical Transformation'
    ELSE
       IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
          !need for Spherical Transformation so we cannot place output in CDAB yet
       ELSE
          !no need for LHS Spherical Transformation
          !need for RHS Horizontal recurrence
       ENDIF
       IF(AngmomA.GE.AngmomB)THEN
          SPEC = 'AtoB'
       ELSE
          SPEC = 'BtoA'
       ENDIF
       call initString(8)
       call AddToString('call HorizontalRR_LHS_P')
       call AddToString(AngmomP)
       call AddToString('A')
       call AddToString(AngmomA)
       call AddToString('B')
       call AddToString(AngmomB)
       call AddToString(SPEC)
       IF(Gen)THEN
          call AddToString('(nContP*nPasses,')
       ELSE
          call AddToString('(nPasses,')
       ENDIF
       call AddToString(nTUVP)
       call AddToString(',Pdistance12,')
       call AddToString(STRINGIN)
       call AddToString(',')
       call AddToString(STRINGOUT)
       call AddToString(',lupri)')
       call writeString(LUMOD3)
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF


    IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
       !need for RHS horizontal recurrence
       call initString(8)
       call AddToString('call SphericalContractOBS1_maxAngP')
       call AddToString(AngmomP)
       call AddToString('_maxAngA')
       call AddToString(AngmomA)
       call AddToString('(')
       call AddToString(nTUVP)
       IF(Gen)THEN
          call AddToString(',nContP*nPasses,')
       ELSE
          call AddToString(',nPasses,')
       ENDIF
       call AddToString(STRINGIN)
       call AddToString(',')
       call AddToString(STRINGOUT)
       call AddToString(')')
       call writeString(LUMOD3)
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ELSE
       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation LHS needed'
    ENDIF

    IF(AngmomP.EQ.0)THEN
       WRITE(LUMOD3,'(A)')'        !no need for RHS Horizontal recurrence relations '
       !there will not be need for Spherical either 
    ELSE
       IF(AngmomA.GE.AngmomB)THEN
          SPEC = 'CtoD'
       ELSE
          SPEC = 'DtoC'
       ENDIF
       call initString(8)
       call AddToString('call HorizontalRR_RHS_Q')
       call AddToString(AngmomP)
       call AddToString('C')
       call AddToString(AngmomA)
       call AddToString('D')
       call AddToString(AngmomB)
       call AddToString(SPEC)
       IF(Gen)THEN
          call AddToString('(nContP,nPasses,')
       ELSE
          call AddToString('(1,nPasses,')
       ENDIF
       call AddToString(nlmA*nlmB)
       call AddToString(',Pdistance12,')
       call AddToString(STRINGIN)
       call AddToString(',')
       call AddToString(STRINGOUT)
       call AddToString(',lupri)')
       call writeString(LUMOD3)
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF

    IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
       !       WRITE(LUMOD3,'(A)')'        !Spherical Transformation RHS'
       call initString(8)
       call AddToString('call SphericalContractOBS2_maxAngQ')
       call AddToString(AngmomP)
       call AddToString('_maxAngC')
       call AddToString(AngmomA)
       call AddToString('(')
       call AddToString(nlmA*nlmB)
       IF(Gen)THEN
          call AddToString(',nContP*nPasses,')
       ELSE
          call AddToString(',nPasses,')
       ENDIF
       call AddToString(STRINGIN)
       call AddToString(',')
       call AddToString(STRINGOUT)
       call AddToString(')')
       call writeString(LUMOD3)
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ELSE
       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation RHS needed'
    ENDIF

    call initString(8)
    call AddToString('call ExtractGabElmP')
    call AddToString(nlmA*nlmB)
    call AddToString(BASISSPEC(1:iBasisSpec))
    call AddToString('(')
    call AddToString(STRINGIN)
    call AddToString(',')
    call AddToString('CDAB,')
    IF(Gen)THEN
       call AddToString('nContP,nPasses)')
    ELSE
       call AddToString('nPasses)')
    ENDIF
    call writeString(LUMOD3)

  end subroutine subroutineMain

  subroutine determineSizes(LUMOD3,AngmomA,AngmomB,STRINGIN,STRINGOUT,TMPSTRING,AngmomP,&
       & nTUVP,nTUVAspec,nTUVBspec,spherical,Gen,Seg)
    implicit none
    integer,intent(in) :: LUMOD3,AngmomA,AngmomB,AngmomP
    integer,intent(in) :: nTUVP,nTUVAspec,nTUVBspec
    character(len=9) :: STRINGIN,STRINGOUT,TMPSTRING
    logical :: spherical,OutputSet,Gen,Seg,Contracted
    character(len=8) :: BASISSPEC
    integer :: iBasisSpec,AngmomQ,AngmomC,AngmomD,nTUVQ,AngmomPQ,nTUVCspec,nTUVDspec
    logical :: PerformSphericaQAndPlaceInTmp
    logical :: PerformHorizontQAndPlaceInTmp
    logical :: PerformSphericaPAndPlaceInTmp
    logical :: PerformHorizontPAndPlaceInTmp
    logical :: PerformBasisContAndPlaceInTmp
    logical :: PerformGabBasisContAndPlaceInTmp
    logical :: PerformTranserAndPlaceInTmp
    logical :: PerformVerticalAndPlaceInTmp
    AngmomQ = AngmomP
    AngmomC=AngmomA
    AngmomD=AngmomB
    nTUVQ=nTUVP
    AngmomPQ=AngmomP+AngmomP
    nTUVCspec=nTUVAspec
    nTUVDspec=nTUVBspec
    nTUV = nTUVP*nTUVQ
    IF(Gen)THEN
       iBasisSpec = 3
       BASISSPEC = 'Gen     '
    ELSEIF(Seg)THEN
       iBasisSpec = 3
       BASISSPEC = 'Seg     '
    ENDIF
    OutputSet = .FALSE.
    Contracted = .FALSE.

    !CDAB always output from ExtractGabElmP 
    IF(Spherical.AND.(AngmomC.GT.1.OR.AngmomD.GT.1))THEN
       PerformSphericaQAndPlaceInTmp = .TRUE. 
    ELSE
       PerformSphericaQAndPlaceInTmp = .FALSE. 
    ENDIF
    IF(AngmomC.EQ.0.AND.AngmomD.EQ.0)THEN
       !no RHS Horizontal Recurrence Relation
       PerformHorizontQAndPlaceInTmp = .FALSE. 
    ELSE
       !RHS Horizontal Recurrence Relation
       PerformHorizontQAndPlaceInTmp = .TRUE. 
    ENDIF
    IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
       PerformSphericaPAndPlaceInTmp = .TRUE. 
    ELSE
       PerformSphericaPAndPlaceInTmp = .FALSE. 
    ENDIF
    IF(AngmomA.EQ.0.AND.AngmomB.EQ.0)THEN
       !no LHS Horizontal Recurrence Relation
       PerformHorizontPAndPlaceInTmp = .FALSE. 
    ELSE
       !LHS Horizontal Recurrence Relation
       PerformHorizontPAndPlaceInTmp = .TRUE. 
    ENDIF
    IF(Seg)THEN
       PerformBasisContAndPlaceInTmp = .FALSE. 
       PerformGabBasisContAndPlaceInTmp = .FALSE. 
    ELSE
       PerformBasisContAndPlaceInTmp = .TRUE. 
       PerformGabBasisContAndPlaceInTmp = .TRUE. 
    ENDIF
    IF(AngmomP.GT.0.AND.AngmomQ.GT.0)THEN
       PerformTranserAndPlaceInTmp = .TRUE. 
    ELSE
       PerformTranserAndPlaceInTmp = .FALSE. 
    ENDIF
    PerformVerticalAndPlaceInTmp = .TRUE. 

    IF(PerformVerticalAndPlaceInTmp)THEN
!       WRITE(LUMOD3,'(A)')'       ! Vertical Recurrence'
       call initString(7)
       call AddToString(STRINGOUT)
       call AddToString('maxSize = MAX(')
       call AddToString(STRINGOUT)
       call AddToString('maxSize,')
       call AddToString(nTUV)       
       IF(PerformTranserAndPlaceInTmp)THEN
          call AddToString('*nPrimP*nPrimP)')
       ELSE 
          IF(Gen)THEN
             call AddToString('*nPrimP*nPrimP)')
          ELSEIF(Seg)THEN
             call AddToString(')')
          ENDIF
       ENDIF
       call writeString(LUMOD3)                
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF
    IF(PerformTranserAndPlaceInTmp)THEN
       call initString(7)
       call AddToString(STRINGOUT)
       call AddToString('maxSize = MAX(')
       call AddToString(STRINGOUT)
       call AddToString('maxSize,')
       call AddToString(nTUVQ*nTUVP)       
       IF(Gen)THEN
          call AddToString('*nPrimP*nPrimP)')
       ELSEIF(Seg)THEN
          call AddToString(')')
       ENDIF
       call writeString(LUMOD3)                
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF
    IF(PerformBasisContAndPlaceInTmp)THEN
!       WRITE(LUMOD3,'(A)')'       ! BasisCont Recurrence'
       call initString(7)
       call AddToString(STRINGOUT)
       call AddToString('maxSize = MAX(')
       call AddToString(STRINGOUT)
       call AddToString('maxSize,')
       call AddToString(nTUVQ*nTUVP)       
       IF(Gen)THEN
!          call AddToString('*nContP*nContP)')
          call AddToString('*nContP)')
       ELSE
          call AddToString(')')
       ENDIF
       call writeString(LUMOD3)                
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF
!!$    IF(PerformGabBasisContAndPlaceInTmp)THEN
!!$       !additional exctraction
!!$       IF(Gen)THEN
!!$          call initString(7)
!!$          call AddToString(STRINGOUT)
!!$          call AddToString('maxSize = MAX(')
!!$          call AddToString(STRINGOUT)
!!$          call AddToString('maxSize,')
!!$          call AddToString(nTUVQ*nTUVP)       
!!$          call AddToString('*nContP)')
!!$          call writeString(LUMOD3)                
!!$          !swap 
!!$          TMPSTRING = STRINGIN
!!$          STRINGIN  = STRINGOUT
!!$          STRINGOUT  = TMPSTRING
!!$       ENDIF
!!$    ENDIF
    IF(PerformHorizontPAndPlaceInTmp)THEN
!       WRITE(LUMOD3,'(A)')'       ! Horizontal LHS Recurrence'
       call initString(7)
       call AddToString(STRINGOUT)
       call AddToString('maxSize = MAX(')
       call AddToString(STRINGOUT)
       call AddToString('maxSize,')
       call AddToString(nTUVQ*nTUVAspec*nTUVBspec)       
       IF(Gen)THEN
          call AddToString('*nContP)')
       ELSE
          call AddToString(')')
       ENDIF
       call writeString(LUMOD3)                
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING       
    ENDIF
    IF(PerformSphericaPAndPlaceInTmp)THEN
!       WRITE(LUMOD3,'(A)')'       ! Spherical LHS'
       call initString(7)
       call AddToString(STRINGOUT)
       call AddToString('maxSize = MAX(')
       call AddToString(STRINGOUT)
       call AddToString('maxSize,')
       call AddToString(nTUVQ*nlmA*nlmB)       
       IF(Gen)THEN
          call AddToString('*nContP)')
       ELSE
          call AddToString(')')
       ENDIF
       call writeString(LUMOD3)                
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF
    IF(PerformHorizontQAndPlaceInTmp)THEN
!       WRITE(LUMOD3,'(A)')'       ! Horizontal RHS Recurrence'
       call initString(7)
       call AddToString(STRINGOUT)
       call AddToString('maxSize = MAX(')
       call AddToString(STRINGOUT)
       call AddToString('maxSize,')
       call AddToString(nTUVCspec*nTUVDspec*nlmA*nlmB)       
       IF(Gen)THEN
          call AddToString('*nContP)')
       ELSE
          call AddToString(')')
       ENDIF
       call writeString(LUMOD3)                
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF
    IF(PerformSphericaQAndPlaceInTmp)THEN
!       WRITE(LUMOD3,'(A)')'       ! Spherical LHS'
       call initString(7)
       call AddToString(STRINGOUT)
       call AddToString('maxSize = MAX(')
       call AddToString(STRINGOUT)
       call AddToString('maxSize,')
       call AddToString(nlmA*nlmB*nlmA*nlmB)       
       IF(Gen)THEN
          call AddToString('*nContP)')
       ELSE
          call AddToString(')')
       ENDIF
       call writeString(LUMOD3)                
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF
  end subroutine determineSizes

END PROGRAM TUV

!contractecoeff_gen

