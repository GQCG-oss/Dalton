PROGRAM TUV
  use math
  use stringsMODULE
  implicit none
  logical :: nPrimLast
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
  integer :: nTUVA,nTUVB,nTUVC,nTUVD,GPUrun
  integer :: nlmA,nlmB,nlmC,nlmD,angmomID,iseg,ILUMOD,I,nUniquenTUVs
  real(realk),pointer :: uniqeparam(:)
  character(len=15),pointer :: uniqeparamNAME(:)
  character(len=9) :: STRINGIN,STRINGOUT,TMPSTRING
  character(len=4) :: SPEC
  character(len=3) :: ARCSTRING
  logical :: BUILD(0:2,0:2,0:2,0:2),Gen,Seg,SegP,segQ,Seg1Prim,UNIQUE,CPU
  integer,pointer :: UniquenTUVs(:)
  !TODO
  !remove mem_alloc
  !remove LOCALINTS = TMParray2
  !add PrimitiveContractionSeg to Transfer or Vertical 
  !
DO GPUrun=1,2
   CPU = .TRUE.
   IF(GPUrun.EQ.2)CPU = .FALSE.
   nPrimLAST = .FALSE.
   IF(CPU)nPrimLAST = .TRUE.
   IF(CPU)THEN
      ARCSTRING = 'CPU'
   ELSE
      ARCSTRING = 'GPU'
   ENDIF
  LUMOD2=2
  open(unit = LUMOD2, file='MAIN_'//ARCSTRING//'_OBS_DRIVER.f90',status="unknown")
  WRITE(LUMOD2,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralMod'
  WRITE(LUMOD2,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'

  LUMOD3=3
  open(unit = LUMOD3, file='MAIN_'//ARCSTRING//'_OBS_DRIVERGen.f90',status="unknown")
  WRITE(LUMOD3,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGen'
  WRITE(LUMOD3,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD3,'(A)')'!Contains routines for General Contracted Basisset '

  LUMOD4=4
  open(unit = LUMOD4, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSegQ.f90',status="unknown")
  WRITE(LUMOD4,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQ'
  WRITE(LUMOD4,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD4,'(A)')'!Contains routines for a General Contracted LHS Segmented contracted RHS and Basisset '

  LUMOD5=5
  open(unit = LUMOD5, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSegP.f90',status="unknown")
  WRITE(LUMOD5,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegP'
  WRITE(LUMOD5,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD5,'(A)')'!Contains routines for Segmented contracted LHS and a General Contracted RHS Basisset '

  LUMOD6=6
  open(unit = LUMOD6, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSeg.f90',status="unknown")
  WRITE(LUMOD6,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg'
  WRITE(LUMOD6,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD6,'(A)')'!Contains routines for Segmented contracted Basisset '

  LUMOD7=7
  open(unit = LUMOD7, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSeg1Prim.f90',status="unknown")
  WRITE(LUMOD7,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1Prim'
  WRITE(LUMOD7,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD7,'(A)')'!Contains routines for Segmented contracted Basisset containing a single primitive'

  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGen'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQ'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegP'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1Prim'

  IF(CPU)THEN
     DO ILUMOD=2,7
        WRITE(ILUMOD,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'McMGeneralMod'
     ENDDO
  ENDIF

  DO ILUMOD=2,7
     WRITE(ILUMOD,'(A)')'use IchorprecisionModule'
     WRITE(ILUMOD,'(A)')'use IchorCommonModule'
     WRITE(ILUMOD,'(A)')'use IchorMemory'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_BUILDRJ000ModGen'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_BUILDRJ000ModSeg1Prim'
  ENDDO
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODAGen'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASegQ'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASegP'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASeg'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASeg1Prim'

  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBGen'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSegQ'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSegP'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSeg'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSeg1Prim'

  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDGen'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSegQ'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSegP'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSeg'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSeg1Prim'

  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCGen'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSegQ'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSegP'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSeg'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSeg1Prim'
  !also needed by the Segs
  DO ILUMOD=4,7
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODAGen'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBGen'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCGen'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDGen'
  ENDDO
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCGen'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDGen'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCGen'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDGen'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoAGen'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoAGen'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBGen'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBGen'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegQ'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSegQ'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSegQ'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSegQ'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASegQ'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASegQ'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSegQ'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSegQ'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegP'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSegP'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSegP'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSegP'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASegP'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASegP'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSegP'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSegP'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSeg'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSeg'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSeg'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASeg'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASeg'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSeg'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSeg'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg1Prim'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSeg1Prim'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSeg1Prim'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSeg1Prim'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASeg1Prim'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASeg1Prim'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSeg1Prim'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSeg1Prim'
  DO ILUMOD=3,7
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_HorizontalRecurrenceLHSModAtoB'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_HorizontalRecurrenceLHSModBtoA'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_HorizontalRecurrenceRHSModCtoD'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_HorizontalRecurrenceRHSModDtoC'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_Sphcontract1Mod'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_Sphcontract2Mod'
     WRITE(ILUMOD,'(A)')'  '
     WRITE(ILUMOD,'(A)')'private   '
  ENDDO
  WRITE(LUMOD2,'(A)')'public :: IchorCoulombIntegral_'//ARCSTRING//'_OBS_general,IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_size  '
  WRITE(LUMOD3,'(A)')'public :: IchorCoulombIntegral_'//ARCSTRING//'_OBS_Gen,IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeGen  '
  WRITE(LUMOD4,'(A)')'public :: IchorCoulombIntegral_'//ARCSTRING//'_OBS_SegQ,IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSegQ  '
  WRITE(LUMOD5,'(A)')'public :: IchorCoulombIntegral_'//ARCSTRING//'_OBS_SegP,IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSegP  '
  WRITE(LUMOD6,'(A)')'public :: IchorCoulombIntegral_'//ARCSTRING//'_OBS_Seg,IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSeg  '
  WRITE(LUMOD7,'(A)')'public :: IchorCoulombIntegral_'//ARCSTRING//'_OBS_Seg1Prim,IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSeg1Prim  '

  DO ILUMOD=2,7
     WRITE(ILUMOD,'(A)')'  '
     WRITE(ILUMOD,'(A)')'CONTAINS'
     WRITE(ILUMOD,'(A)')'  '
     WRITE(ILUMOD,'(A)')'  '
  ENDDO
  WRITE(LUMOD2,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&'
  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&'
  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&'
  WRITE(LUMOD2,'(A)')'       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&'
  WRITE(LUMOD2,'(A)')'       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&'
  IF(CPU)THEN
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)'
  ELSE
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)'
  ENDIF
  WRITE(LUMOD2,'(A)')'    implicit none'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nAtomsA,nAtomsB'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ)'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)'
  WRITE(LUMOD2,'(A)')'    logical,intent(in)     :: Qsegmented,Psegmented'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: pcent(3*nPrimP*nAtomsA*nAtomsB)   !pcent(3,nPrimP)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: qcent(3*nPrimQ)           !qcent(3,nPrimQ)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: QpreExpFac(nPrimQ),PpreExpFac(nPrimP*nAtomsA*nAtomsB)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)'
  WRITE(LUMOD2,'(A)')'    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
  WRITE(LUMOD2,'(A)')'    !    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)'
  WRITE(LUMOD2,'(A)')'    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
  WRITE(LUMOD2,'(A)')'    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: localintsmaxsize'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(inout) :: LOCALINTS(localintsmaxsize)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: integralPrefactor(nPrimQP)'
  WRITE(LUMOD2,'(A)')'    logical,intent(in) :: PQorder'
  WRITE(LUMOD2,'(A)')'    !integralPrefactor(nPrimP,nPrimQ)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: reducedExponents(nPrimQP)'
  WRITE(LUMOD2,'(A)')'    !reducedExponents(nPrimP,nPrimQ)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: Qdistance12(3) !Ccenter-Dcenter'
  WRITE(LUMOD2,'(A)')'    !Qdistance12(3)'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: Pdistance12(3*nAtomsA*nAtomsB)  !Acenter-Bcenter '
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: Acenter(3,nAtomsA),Bcenter(3,nAtomsB),Ccenter(3),Dcenter(3)'
  WRITE(LUMOD2,'(A)')'    logical,intent(in) :: spherical'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize'
  WRITE(LUMOD2,'(A)')'!   TMP variables - allocated outside'  
  WRITE(LUMOD2,'(A)')'    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)'
  IF(.NOT.CPU)WRITE(LUMOD2,'(A)')'    integer(kind=acckind),intent(in) :: iASync'


  !  WRITE(LUMOD2,'(A)')'!   Local variables '  
  !  WRITE(LUMOD2,'(A)')'    real(realk),pointer :: squaredDistance(:)'!,Rpq(:)'!,Rqc(:),Rpa(:)
  !  WRITE(LUMOD2,'(A)')'    integer :: AngmomPQ,AngmomP,AngmomQ,I,J,nContQP,la,lb,lc,ld,nsize,angmomid'
  !  WRITE(LUMOD2,'(A)')'    real(realk),pointer :: RJ000(:),OUTPUTinterest(:)'


!  WRITE(LUMOD2,'(A)')'  '
!  WRITE(LUMOD2,'(A)')'    IF(nAtomsC*nAtomsD.NE.nPasses)Call ichorquit(''nPass error'',-1)'
  WRITE(LUMOD2,'(A)')'  '
!!$  WRITE(LUMOD2,'(A)')'!IF(.TRUE.)THEN'
!!$  WRITE(LUMOD2,'(A)')'!    call interest_initialize()'
!!$  WRITE(LUMOD2,'(A)')'!'
!!$  WRITE(LUMOD2,'(A)')'!    nsize = size(LOCALINTS)*10'
!!$  WRITE(LUMOD2,'(A)')'!    call mem_ichor_alloc(OUTPUTinterest,nsize)'
!!$  WRITE(LUMOD2,'(A)')'!    nsize = nPrimQP'
!!$  WRITE(LUMOD2,'(A)')'!    '
!!$  WRITE(LUMOD2,'(A)')'!    '
!!$  WRITE(LUMOD2,'(A)')'!    la = AngmomA+1'
!!$  WRITE(LUMOD2,'(A)')'!    lb = AngmomB+1'
!!$  WRITE(LUMOD2,'(A)')'!    lc = AngmomC+1'
!!$  WRITE(LUMOD2,'(A)')'!    ld = AngmomD+1'
!!$  WRITE(LUMOD2,'(A)')'!    call interest_eri(OUTPUTinterest,nsize,&'
!!$  WRITE(LUMOD2,'(A)')'!       & la,Aexp,Acenter(1),Acenter(2),Acenter(3),ACC,&'
!!$  WRITE(LUMOD2,'(A)')'!         & lb,Bexp,Bcenter(1),Bcenter(2),Bcenter(3),BCC,&'
!!$  WRITE(LUMOD2,'(A)')'!         & lc,Cexp,Ccenter(1,1),Ccenter(2,1),Ccenter(3,1),CCC,&'
!!$  WRITE(LUMOD2,'(A)')'!         & ld,Dexp,Dcenter(1,1),Dcenter(2,1),Dcenter(3,1),DCC,&'
!!$  WRITE(LUMOD2,'(A)')'!         & lupri)!,&'
!!$  WRITE(LUMOD2,'(A)')'!         !         & .false.)'
!!$  WRITE(LUMOD2,'(A)')'!    write(lupri,*)''OUTPUTinterest'',OUTPUTinterest'
!!$  WRITE(LUMOD2,'(A)')'!ENDIF'
  !  WRITE(LUMOD2,'(A)')' '
  !  WRITE(LUMOD2,'(A)')'    !build the distance between center P and A in Rpa '
  !  WRITE(LUMOD2,'(A)')'    !used in Vertical and Electron Transfer Recurrence Relations '
  !  WRITE(LUMOD2,'(A)')'    call mem_ichor_alloc(Rpa,3*nPrimP)'
  !  WRITE(LUMOD2,'(A)')'    call build_Rpa(nPrimP,Pcent,Acenter,Rpa)'
  !  WRITE(LUMOD2,'(A)')'    !build the distance between center Q and C in Rqc used'
  !  WRITE(LUMOD2,'(A)')'    !used in Electron Transfer Recurrence Relations '
  !  WRITE(LUMOD2,'(A)')'    call mem_ichor_alloc(Rqc,3*nPrimQ)'
  !  WRITE(LUMOD2,'(A)')'    call build_Rpa(nPrimQ,Qcent,Ccenter,Rqc)'
  !  WRITE(LUMOD2,'(A)')'    '
  WRITE(LUMOD2,'(A)')'    IF(PQorder)THEN'
  WRITE(LUMOD2,'(A)')'       call IchorQuit(''PQorder OBS general expect to get QP ordering'',-1)'
  WRITE(LUMOD2,'(A)')'    ENDIF'
  WRITE(LUMOD2,'(A)')'    IF(.NOT.spherical)THEN'
  WRITE(LUMOD2,'(A)')'       call IchorQuit(''cartesian not testet'',-1)'
  WRITE(LUMOD2,'(A)')'    ENDIF'
  WRITE(LUMOD2,'(A)')'    '
  !  WRITE(LUMOD2,'(A)')'    call mem_ichor_alloc(squaredDistance,nPasses*nPrimQP)'
  !  WRITE(LUMOD2,'(A)')'    call mem_ichor_alloc(Rpq,3*nPrimQP*nPasses)'
  !  WRITE(LUMOD2,'(A)')'    !builds squaredDistance between center P and Q. Order(nPrimQ,nPrimP,nPassQ)'
  !  WRITE(LUMOD2,'(A)')'    !builds Distance between center P and Q in Rpq. Order(3,nPrimQ,nPrimP,nPassQ)'
  !  WRITE(LUMOD2,'(A)')'    call build_QP_squaredDistance_and_Rpq(nPrimQ,nPasses,nPrimP,Qcent,Pcent,&'
  !  WRITE(LUMOD2,'(A)')'         & squaredDistance,Rpq)'


  WRITE(LUMOD2,'(A)')'   IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN'
  WRITE(LUMOD2,'(A)')'    call IchorCoulombIntegral_'//ARCSTRING//'_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&'
  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&'
  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&'
  WRITE(LUMOD2,'(A)')'       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&'
  WRITE(LUMOD2,'(A)')'       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&'
  IF(CPU)THEN
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
  ELSE
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
  ENDIF
  WRITE(LUMOD2,'(A)')'   ELSEIF(Psegmented.AND.Qsegmented)THEN'
  WRITE(LUMOD2,'(A)')'    call IchorCoulombIntegral_'//ARCSTRING//'_OBS_Seg(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&'
  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&'
  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&'
  WRITE(LUMOD2,'(A)')'       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&'
  WRITE(LUMOD2,'(A)')'       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&'
  IF(CPU)THEN
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
  ELSE
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
  ENDIF
  WRITE(LUMOD2,'(A)')'   ELSEIF(Psegmented)THEN'
  WRITE(LUMOD2,'(A)')'    call IchorCoulombIntegral_'//ARCSTRING//'_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&'
  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&'
  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&'
  WRITE(LUMOD2,'(A)')'       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&'
  WRITE(LUMOD2,'(A)')'       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&'
  IF(CPU)THEN
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
  ELSE
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
  ENDIF
  WRITE(LUMOD2,'(A)')'   ELSEIF(Qsegmented)THEN'
  WRITE(LUMOD2,'(A)')'    call IchorCoulombIntegral_'//ARCSTRING//'_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&'
  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&'
  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&'
  WRITE(LUMOD2,'(A)')'       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&'
  WRITE(LUMOD2,'(A)')'       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&'
  IF(CPU)THEN
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
  ELSE
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
  ENDIF

  WRITE(LUMOD2,'(A)')'   ELSE'
  WRITE(LUMOD2,'(A)')'    call IchorCoulombIntegral_'//ARCSTRING//'_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&'
  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&'
  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&'
  WRITE(LUMOD2,'(A)')'       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&'
  WRITE(LUMOD2,'(A)')'       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&'
  IF(CPU)THEN
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
  ELSE
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
  ENDIF

  WRITE(LUMOD2,'(A)')'   ENDIF'
  WRITE(LUMOD2,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general'
  WRITE(LUMOD2,'(A)')'  '

  DO ISEG = 1,5
     Gen=.FALSE.; SegQ=.FALSE.; SegP=.FALSE.; Seg=.FALSE.; Seg1Prim=.FALSE.
     IF(ISEG.EQ.1)THEN
        Seg1Prim=.TRUE.
     ELSEIF(ISEG.EQ.2)THEN
        Seg=.TRUE.
     ELSEIF(ISEG.EQ.3)THEN
        SegP=.TRUE.
     ELSEIF(ISEG.EQ.4)THEN
        SegQ=.TRUE.
     ELSEIF(ISEG.EQ.5)THEN
        Gen=.TRUE.       
     ENDIF

     IF(Gen)THEN
        WRITE(LUMOD3,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&'
        ILUMOD = 3
     ELSEIF(SegQ)THEN
        WRITE(LUMOD4,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&'
        ILUMOD = 4
     ELSEIF(SegP)THEN
        WRITE(LUMOD5,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&'
        ILUMOD = 5
     ELSEIF(Seg)THEN
        WRITE(LUMOD6,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_Seg(nPrimA,nPrimB,nPrimC,nPrimD,&'
        ILUMOD = 6
     ELSEIF(Seg1Prim)THEN
        WRITE(LUMOD7,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&'
        ILUMOD = 7
     ENDIF
     WRITE(ILUMOD,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
     WRITE(ILUMOD,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
     WRITE(ILUMOD,'(A)')'       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
     WRITE(ILUMOD,'(A)')'       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&'
     WRITE(ILUMOD,'(A)')'       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&'
     WRITE(ILUMOD,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&'
     WRITE(ILUMOD,'(A)')'       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&'
     WRITE(ILUMOD,'(A)')'       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&'
     IF(CPU)THEN
        WRITE(ILUMOD,'(A)')'       & IatomAPass,iatomBPass)'
     ELSE
        WRITE(ILUMOD,'(A)')'       & IatomAPass,iatomBPass,iASync)'
     ENDIF
     WRITE(ILUMOD,'(A)')'    implicit none'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nAtomsA,nAtomsB'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ)'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)'
     WRITE(ILUMOD,'(A)')'    logical,intent(in)     :: Qsegmented,Psegmented'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: pcent(3*nPrimP*nAtomsA*nAtomsB)           !qcent(3,nPrimP)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: qcent(3*nPrimQ)           !qcent(3,nPrimQ)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: QpreExpFac(nPrimQ),PpreExpFac(nPrimP*nAtomsA*nAtomsB)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)'
     WRITE(ILUMOD,'(A)')'    !    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
     WRITE(ILUMOD,'(A)')'    !    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)'
     WRITE(ILUMOD,'(A)')'    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)'
     WRITE(ILUMOD,'(A)')'    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: localintsmaxsize'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(inout) :: LOCALINTS(localintsmaxsize)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: integralPrefactor(nPrimQP)'
     WRITE(ILUMOD,'(A)')'    logical,intent(in) :: PQorder'
     WRITE(ILUMOD,'(A)')'    !integralPrefactor(nPrimP,nPrimQ)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: reducedExponents(nPrimQP)'
     WRITE(ILUMOD,'(A)')'    !reducedExponents(nPrimP,nPrimQ)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: Qdistance12(3) !Ccenter-Dcenter'
     WRITE(ILUMOD,'(A)')'    !Qdistance12(3)'
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: Pdistance12(3*nAtomsA*nAtomsB)           !Acenter-Bcenter '
     WRITE(ILUMOD,'(A)')'    real(realk),intent(in) :: Acenter(3,nAtomsA),Bcenter(3,nAtomsB),Ccenter(3),Dcenter(3)'
     WRITE(ILUMOD,'(A)')'    logical,intent(in) :: spherical'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: TMParray1maxsize,TMParray2maxsize'
     WRITE(ILUMOD,'(A)')'!   TMP variables - allocated outside'  
     WRITE(ILUMOD,'(A)')'    real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)'
     IF(.NOT.CPU)WRITE(ILUMOD,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(ILUMOD,'(A)')'!   Local variables '  
     !  WRITE(ILUMOD,'(A)')'    real(realk),pointer :: squaredDistance(:)'!,Rpq(:)'!,Rqc(:),Rpa(:)
     WRITE(ILUMOD,'(A)')'    integer :: AngmomPQ,AngmomP,AngmomQ,I,J,nContQP,la,lb,lc,ld,nsize,angmomid'


     WRITE(ILUMOD,'(A)')'    '
     WRITE(ILUMOD,'(A)')'    !Setup combined Angmom info'
     WRITE(ILUMOD,'(A)')'    AngmomP = AngmomA+AngmomB'
     WRITE(ILUMOD,'(A)')'    AngmomQ = AngmomC+AngmomD'
     WRITE(ILUMOD,'(A)')'    AngmomPQ  = AngmomP + AngmomQ'
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
     WRITE(ILUMOD,'(A)')'!    nTUV = (AngmomPQ+1)*(AngmomPQ+2)*(AngmomPQ+3)/6'
     WRITE(ILUMOD,'(A)')'!    nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6'
     WRITE(ILUMOD,'(A)')'!    nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6'
     WRITE(ILUMOD,'(A)')'!    nTUVC = (AngmomC+1)*(AngmomC+2)*(AngmomC+3)/6'
     WRITE(ILUMOD,'(A)')'!    nTUVD = (AngmomD+1)*(AngmomD+2)*(AngmomD+3)/6'
     WRITE(ILUMOD,'(A)')'!    nlmA = 2*AngmomA+1'
     WRITE(ILUMOD,'(A)')'!    nlmB = 2*AngmomB+1'
     WRITE(ILUMOD,'(A)')'!    nlmC = 2*AngmomC+1'
     WRITE(ILUMOD,'(A)')'!    nlmD = 2*AngmomD+1'
     WRITE(ILUMOD,'(A)')'    AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD'
     WRITE(ILUMOD,'(A)')'    SELECT CASE(AngmomID)'

     DO AngmomA = 0,2
        DO AngmomB = 0,2
           DO AngmomC = 0,2
              DO AngmomD = 0,2
                 BUILD(AngmomA,AngmomB,AngmomC,AngmomD) = .TRUE.
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !the order may be important and should be 
     DO AngmomA = 0,2
        DO AngmomB = 0,AngmomA
           DO AngmomC = 0,AngmomA
              DO AngmomD = 0,AngmomC

                 BUILD(AngmomA,AngmomB,AngmomC,AngmomD) = .FALSE.

                 !==========================00=====================================

                 AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD
                 WRITE(ILUMOD,'(A,I4,A,I2,A,I2,A,I2,A,I2,A)')'    CASE(',AngmomID,')  !Angmom(A=',AngmomA,',B=',AngmomB,',C=',AngmomC,',D=',AngmomD,') combi'
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
                    !         WRITE(ILUMOD,'(A)')'      IF(spherical)THEN'
                    call subroutineMAIN(ILUMOD,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
                         & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,nlmA,nlmB,nlmC,nlmD,spherical,Gen,SegQ,SegP,Seg,Seg1Prim)
                    !         WRITE(ILUMOD,'(A)')'      ELSE'
                    !         spherical = .FALSE.
                    !         nlmA = nTUVAspec
                    !         nlmB = nTUVBspec
                    !         nlmC = nTUVCspec
                    !         nlmD = nTUVDspec
                    !         STRINGIN(1:9)  = 'TMParray1'
                    !         STRINGOUT(1:9) = 'TMParray2'
                    !         TMPSTRING(1:9) = '         '
                    !         call subroutineMAIN(ILUMOD,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
                    !              & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,spherical)
                    !         WRITE(ILUMOD,'(A)')'      ENDIF'
                 ELSE
                    spherical = .TRUE.
                    nlmA = 2*AngmomA+1
                    nlmB = 2*AngmomB+1
                    nlmC = 2*AngmomC+1
                    nlmD = 2*AngmomD+1
                    STRINGIN(1:9)  = 'TMParray1'
                    STRINGOUT(1:9) = 'TMParray2'
                    TMPSTRING(1:9) = '         '
                    call subroutineMAIN(ILUMOD,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
                         & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,nlmA,nlmB,nlmC,nlmD,spherical,Gen,SegQ,SegP,Seg,Seg1Prim)
                 ENDIF


                 !==========================00=====================================

              ENDDO
           ENDDO
        ENDDO
     ENDDO


     DO AngmomA = 0,2
        DO AngmomB = 0,2
           DO AngmomC = 0,2
              DO AngmomD = 0,2
                 IF(BUILD(AngmomA,AngmomB,AngmomC,AngmomD))THEN

                    !==========================00=====================================

                    AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD
                    WRITE(ILUMOD,'(A,I4,A,I2,A,I2,A,I2,A,I2,A)')'    CASE(',AngmomID,')  !Angmom(A=',AngmomA,',B=',AngmomB,',C=',AngmomC,',D=',AngmomD,') combi'
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
                       !         WRITE(ILUMOD,'(A)')'      IF(spherical)THEN'
                       call subroutineMAIN(ILUMOD,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
                            & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,nlmA,nlmB,nlmC,nlmD,spherical,Gen,SegQ,SegP,Seg,Seg1Prim)
                       !         WRITE(ILUMOD,'(A)')'      ELSE'
                       !         spherical = .FALSE.
                       !         nlmA = nTUVAspec
                       !         nlmB = nTUVBspec
                       !         nlmC = nTUVCspec
                       !         nlmD = nTUVDspec
                       !         STRINGIN(1:9)  = 'TMParray1'
                       !         STRINGOUT(1:9) = 'TMParray2'
                       !         TMPSTRING(1:9) = '         '
                       !         call subroutineMAIN(ILUMOD,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
                       !              & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,spherical)
                       !         WRITE(ILUMOD,'(A)')'      ENDIF'
                    ELSE
                       spherical = .TRUE.
                       nlmA = 2*AngmomA+1
                       nlmB = 2*AngmomB+1
                       nlmC = 2*AngmomC+1
                       nlmD = 2*AngmomD+1
                       STRINGIN(1:9)  = 'TMParray1'
                       STRINGOUT(1:9) = 'TMParray2'
                       TMPSTRING(1:9) = '         '
                       call subroutineMAIN(ILUMOD,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
                            & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,nlmA,nlmB,nlmC,nlmD,spherical,Gen,SegQ,SegP,Seg,Seg1Prim)
                    ENDIF

                    !==========================00=====================================

                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
        WRITE(ILUMOD,'(A)')'    CASE DEFAULT'
     IF(CPU)THEN
        WRITE(ILUMOD,'(A)')'        call IchorCoulombIntegral_'//ARCSTRING//'_McM_general(nPrimA,nPrimB,nPrimC,nPrimD,&'
        WRITE(ILUMOD,'(A)')'           & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
        WRITE(ILUMOD,'(A)')'           & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
        WRITE(ILUMOD,'(A)')'           & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
        WRITE(ILUMOD,'(A)')'           & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&'
        WRITE(ILUMOD,'(A)')'           & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&'
        WRITE(ILUMOD,'(A)')'           & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&'
        WRITE(ILUMOD,'(A)')'           & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&'
        WRITE(ILUMOD,'(A)')'           & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&'
        IF(CPU)THEN
           WRITE(ILUMOD,'(A)')'           & IatomAPass,iatomBPass)' 
        ELSE
           WRITE(ILUMOD,'(A)')'           & IatomAPass,iatomBPass,iASync)' 
        ENDIF
     ELSE
        IF(Gen)THEN
           WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in IchorCoulombIntegral_'//ARCSTRING//'_OBS_Gen'',-1)'
        ELSEIF(SegQ)THEN
           WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in IchorCoulombIntegral_'//ARCSTRING//'_OBS_SegQ'',-1)'
        ELSEIF(SegP)THEN
           WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in IchorCoulombIntegral_'//ARCSTRING//'_OBS_SegP'',-1)'
        ELSEIF(Seg)THEN
           WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in IchorCoulombIntegral_'//ARCSTRING//'_OBS_Seg'',-1)'
        ELSEIF(Seg1Prim)THEN
           WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in IchorCoulombIntegral_'//ARCSTRING//'_OBS_Seg1Prim'',-1)'
        ENDIF
     ENDIF
     WRITE(ILUMOD,'(A)')'    END SELECT'
     IF(Gen)THEN
        WRITE(ILUMOD,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_Gen'
     ELSEIF(SegQ)THEN
        WRITE(ILUMOD,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_SegQ'
     ELSEIF(SegP)THEN
        WRITE(ILUMOD,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_SegP'
     ELSEIF(Seg)THEN
        WRITE(ILUMOD,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_Seg'
     ELSEIF(Seg1Prim)THEN
        WRITE(ILUMOD,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_Seg1Prim'
     ENDIF
     WRITE(ILUMOD,'(A)')'  '
  ENDDO

  WRITE(LUMOD2,'(A)')'  '
  WRITE(LUMOD2,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_size(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD2,'(A)')'         & nContA,nContB,nContC,nContD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,&'
  WRITE(LUMOD2,'(A)')'         & nContP,nContQ,nPrimQP,nContQP,Psegmented,Qsegmented)'
  WRITE(LUMOD2,'(A)')'    implicit none'
  WRITE(LUMOD2,'(A)')'    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nContA,nContB,nContC,nContD'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nPrimA,nPrimB,nPrimC,nPrimD'
  WRITE(LUMOD2,'(A)')'    logical,intent(in) :: Psegmented,Qsegmented'
  WRITE(LUMOD2,'(A)')'    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN'
  WRITE(LUMOD2,'(A)')'     call IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSeg1Prim(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD2,'(A)')'         & nContA,nContB,nContC,nContD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)'
  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented.AND.Qsegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSeg(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD2,'(A)')'         & nContA,nContB,nContC,nContD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)'
  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSegP(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD2,'(A)')'         & nContA,nContB,nContC,nContD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)'
  WRITE(LUMOD2,'(A)')'    ELSEIF(Qsegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSegQ(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD2,'(A)')'         & nContA,nContB,nContC,nContD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)'
  WRITE(LUMOD2,'(A)')'    ELSE'
  WRITE(LUMOD2,'(A)')'     call IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeGen(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD2,'(A)')'         & nContA,nContB,nContC,nContD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)'
  WRITE(LUMOD2,'(A)')'    ENDIF'
  WRITE(LUMOD2,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_size'
  WRITE(LUMOD2,'(A)')'  '


  DO ISEG = 1,5
     Gen=.FALSE.; SegQ=.FALSE.; SegP=.FALSE.; Seg=.FALSE.; Seg1Prim=.FALSE.
     IF(ISEG.EQ.1)THEN
        Seg1Prim=.TRUE.
     ELSEIF(ISEG.EQ.2)THEN
        Seg=.TRUE.
     ELSEIF(ISEG.EQ.3)THEN
        SegP=.TRUE.
     ELSEIF(ISEG.EQ.4)THEN
        SegQ=.TRUE.
     ELSEIF(ISEG.EQ.5)THEN
        Gen=.TRUE.       
     ENDIF

     WRITE(ISEG+2,'(A)')'  '
     IF(Gen)THEN
        WRITE(LUMOD3,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeGen(TMParray1maxsize,&'
        ILUMOD = LUMOD3
     ELSEIF(SegQ)THEN
        WRITE(LUMOD4,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSegQ(TMParray1maxsize,&'
        ILUMOD = LUMOD4
     ELSEIF(SegP)THEN
        WRITE(LUMOD5,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSegP(TMParray1maxsize,&'
        ILUMOD = LUMOD5
     ELSEIF(Seg)THEN
        WRITE(LUMOD6,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSeg(TMParray1maxsize,&'
        ILUMOD = LUMOD6
     ELSEIF(Seg1Prim)THEN
        WRITE(LUMOD7,'(A)')'  subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSeg1Prim(TMParray1maxsize,&'
        ILUMOD = LUMOD7
     ENDIF
     WRITE(ILUMOD,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,nContA,nContB,nContC,nContD,&'
     WRITE(ILUMOD,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)'
     WRITE(ILUMOD,'(A)')'    implicit none'
     WRITE(ILUMOD,'(A)')'    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nContA,nContB,nContC,nContD'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nPrimA,nPrimB,nPrimC,nPrimD'
     WRITE(ILUMOD,'(A)')'    ! local variables'
     WRITE(ILUMOD,'(A)')'    integer :: AngmomID'
     WRITE(ILUMOD,'(A)')'    '
     WRITE(ILUMOD,'(A)')'    AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD'
     WRITE(ILUMOD,'(A)')'    TMParray2maxSize = 1'
     WRITE(ILUMOD,'(A)')'    TMParray1maxSize = 1'

     WRITE(ILUMOD,'(A)')'    SELECT CASE(AngmomID)'  
     DO AngmomA = 0,2
        DO AngmomB = 0,2
           DO AngmomC = 0,2
              DO AngmomD = 0,2
                 !   DO AngmomB = 0,AngmomA
                 !    DO AngmomC = 0,AngmomA
                 !     DO AngmomD = 0,AngmomC
                 AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD
                 WRITE(ILUMOD,'(A,I4,A,I2,A,I2,A,I2,A,I2,A)')'    CASE(',AngmomID,')  !Angmom(A=',AngmomA,',B=',AngmomB,',C=',AngmomC,',D=',AngmomD,') combi'
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
                 call determineSizes(ILUMOD,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
                      & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,nlmA,nlmB,nlmC,nlmD,spherical,Gen,SegQ,SegP,Seg,Seg1Prim)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     WRITE(ILUMOD,'(A)')'    CASE DEFAULT'
     IF(CPU)THEN
        WRITE(ILUMOD,'(A)')'     call IchorCoulombIntegral_'//ARCSTRING//'_McM_general_size(TMParray1maxsize,&'
        WRITE(ILUMOD,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
        WRITE(ILUMOD,'(A)')'         & nContA,nContB,nContC,nContD,&'
        WRITE(ILUMOD,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,&'
        WRITE(ILUMOD,'(A)')'         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP,&'
        IF(Gen)THEN
           WRITE(ILUMOD,'(A)')'         & .FALSE.,.FALSE.)'
        ELSEIF(SegQ)THEN
           WRITE(ILUMOD,'(A)')'         & .FALSE.,.TRUE.)'
        ELSEIF(SegP)THEN
           WRITE(ILUMOD,'(A)')'         & .TRUE.,.FALSE.)'
        ELSE !Seg or Seg1Prim
           WRITE(ILUMOD,'(A)')'         & .TRUE.,.TRUE.)'
        ENDIF
     ELSE
        WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in IchorCoulombIntegral_OBS_general_size'',-1)'
     ENDIF
     WRITE(ILUMOD,'(A)')'    END SELECT'

     IF(Gen)THEN
        WRITE(LUMOD3,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeGen'
        ILUMOD = LUMOD3
     ELSEIF(SegQ)THEN
        WRITE(LUMOD4,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSegQ'
        ILUMOD = LUMOD4
     ELSEIF(SegP)THEN
        WRITE(LUMOD5,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSegP'
        ILUMOD = LUMOD5
     ELSEIF(Seg)THEN
        WRITE(LUMOD6,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSeg'
        ILUMOD = LUMOD6
     ELSEIF(Seg1Prim)THEN
        WRITE(LUMOD7,'(A)')'  end subroutine IchorCoulombIntegral_'//ARCSTRING//'_OBS_general_sizeSeg1Prim'
        ILUMOD = LUMOD7
     ENDIF
  ENDDO
!====================================QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ

  IF(.NOT.nPrimLast)THEN
     !Only General
     WRITE(LUMOD3,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionC'//ARCSTRING//'Gen')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD3)
     WRITE(LUMOD3,'(A)')'       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,nTUVP,nTUVQ,iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'    implicit none'
     WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ,nTUVP,nTUVQ'
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: CCC(nPrimC,nContC)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD*nPrimA*nPrimB*nPasses*nTUVP*nTUVQ)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPrimA*nPrimB*nPasses*nTUVP*nTUVQ)'
     IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD3,'(A)')'    !'
     WRITE(LUMOD3,'(A)')'    integer :: iP,iContC,iPrimC'

     WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'
     WRITE(LUMOD3,'(A)')'!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP DO &'
        WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iP,iContC,iPrimC,TMP)'
!        WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
!        WRITE(LUMOD3,'(A)')'!$OMP        CCC,AUXarrayCont,AUXarray2,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iP,iContC,iPrimC,TMP) &'
        WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
        WRITE(LUMOD3,'(A)')'!$ACC        CCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'   do iP=1,nTUVP*nTUVQ*nPrimA*nPrimB*nPrimD*nPasses'
     WRITE(LUMOD3,'(A)')'    do iContC=1,nContC'               
     WRITE(LUMOD3,'(A)')'     TMP = 0.0E0_realk'
     WRITE(LUMOD3,'(A)')'     do iPrimC=1,nPrimC'
     WRITE(LUMOD3,'(A)')'      TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)'
     WRITE(LUMOD3,'(A)')'     enddo'     
     WRITE(LUMOD3,'(A)')'     AUXarrayCont(iContC,iP)=TMP'
     WRITE(LUMOD3,'(A)')'    enddo'
     WRITE(LUMOD3,'(A)')'   enddo'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionC'//ARCSTRING//'Gen')
     call writeString(LUMOD3)


     WRITE(LUMOD3,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionD'//ARCSTRING//'Gen')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD3)
     WRITE(LUMOD3,'(A)')'       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,nTUVP,nTUVQ,iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'    implicit none'
     WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ,nTUVP,nTUVQ'
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: DCC(nPrimD,nContD)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray2(nContC,nPrimD,nPrimA*nPrimB*nPasses*nTUVP*nTUVQ)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPrimA*nPrimB*nPasses*nTUVP*nTUVQ)'
     IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD3,'(A)')'    !'
     WRITE(LUMOD3,'(A)')'    integer :: iP,iContC,iContD,iPrimD'
     WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'
     WRITE(LUMOD3,'(A)')'!Scaling p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP DO &'
        WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iP,iContC,iContD,iPrimD,TMP) '
!        WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
!        WRITE(LUMOD3,'(A)')'!$OMP        DCC,AUXarrayCont,AUXarray2,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iP,iContC,iContD,iPrimD,TMP) &'
        WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
        WRITE(LUMOD3,'(A)')'!$ACC         DCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     call initString(3)
     WRITE(LUMOD3,'(A)')'   do iP=1,nTUVP*nTUVQ*nPrimA*nPrimB*nPasses'
     WRITE(LUMOD3,'(A)')'    do iContD=1,nContD'               
     WRITE(LUMOD3,'(A)')'     do iContC=1,nContC'               
     WRITE(LUMOD3,'(A)')'      AUXarrayCont(iContC,iContD,iP)=0.0E0_realk'
     WRITE(LUMOD3,'(A)')'     enddo'
     WRITE(LUMOD3,'(A)')'     do iPrimD=1,nPrimD'
     WRITE(LUMOD3,'(A)')'      TMP = DCC(iPrimD,iContD)'
     WRITE(LUMOD3,'(A)')'      do iContC=1,nContC'               
     WRITE(LUMOD3,'(A)')'       AUXarrayCont(iContC,iContD,iP)=AUXarrayCont(iContC,iContD,iP)+TMP*AUXarray2(iContC,iPrimD,iP)'
     WRITE(LUMOD3,'(A)')'      enddo'     
     WRITE(LUMOD3,'(A)')'     enddo'
     WRITE(LUMOD3,'(A)')'    enddo'
     WRITE(LUMOD3,'(A)')'   enddo'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionD'//ARCSTRING//'Gen')
     call writeString(LUMOD3)

     WRITE(LUMOD3,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionA'//ARCSTRING//'Gen')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD3)
     WRITE(LUMOD3,'(A)')'       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,nTUVP,nTUVQ,iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'    implicit none'
     WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ,nTUVP,nTUVQ'
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: ACC(nPrimA,nContA)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray2(nContC*nContD,nPrimA,nPrimB*nPasses*nTUVP*nTUVQ)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC*nContD,nContA,nPrimB*nPasses*nTUVP*nTUVQ)'
     IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'

     WRITE(LUMOD3,'(A)')'    !'
     WRITE(LUMOD3,'(A)')'    integer :: iP,iC,iContA,iPrimA'
     WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'
     WRITE(LUMOD3,'(A)')'!Scaling p**2*c**3*nTUV*nPassQ: nPrimA*nPrimB*nContC*nContD*nContA*nTUV*nPassQ'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP DO &'
        WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iP,iContA,iPrimA,TMP,iC)'
!        WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContC,nContD,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
!        WRITE(LUMOD3,'(A)')'!$OMP        ACC,AUXarrayCont,AUXarray2,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iP,iContA,iPrimA,TMP,iC) &'
        WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
        WRITE(LUMOD3,'(A)')'!$ACC         ACC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'   do iP=1,nTUVP*nTUVQ*nPrimB*nPasses'
     WRITE(LUMOD3,'(A)')'    do iContA=1,nContA'
     WRITE(LUMOD3,'(A)')'     do iC=1,nContC*nContD'
     WRITE(LUMOD3,'(A)')'      AUXarrayCont(iC,iContA,iP)=0.0E0_realk'
     WRITE(LUMOD3,'(A)')'     enddo'
     WRITE(LUMOD3,'(A)')'     do iPrimA=1,nPrimA'
     WRITE(LUMOD3,'(A)')'      TMP = ACC(iPrimA,iContA)'
     WRITE(LUMOD3,'(A)')'      do iC=1,nContC*nContD'
     WRITE(LUMOD3,'(A)')'       AUXarrayCont(iC,iContA,iP)=AUXarrayCont(iC,iContA,iP)+TMP*AUXarray2(iC,iPrimA,iP)'
     WRITE(LUMOD3,'(A)')'      enddo'     
     WRITE(LUMOD3,'(A)')'     enddo'
     WRITE(LUMOD3,'(A)')'    enddo'
     WRITE(LUMOD3,'(A)')'   enddo'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionA'//ARCSTRING//'Gen')
     call writeString(LUMOD3)

     WRITE(LUMOD3,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionB'//ARCSTRING//'Gen')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD3)
     WRITE(LUMOD3,'(A)')'       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,nTUVP,nTUVQ,iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'    implicit none'
     WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ,nTUVP,nTUVQ'
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: BCC(nPrimB,nContB)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray2(nContC*nContD*nContA,nPrimB,nPasses*nTUVP*nTUVQ)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC*nContD*nContA,nContB,nPasses*nTUVP*nTUVQ)'
     IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD3,'(A)')'    !'
     WRITE(LUMOD3,'(A)')'    integer :: iP,iC,iContB,iPrimB'
     WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'
     WRITE(LUMOD3,'(A)')'    !Scaling p*c**4*nTUV*nPassQ: nPrimB*nContC*nContD*nContA*nContB*nTUV*nPassQ'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP DO &'
        WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iP,iContB,iPrimB,TMP,iC) '
!        WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContC,nContD,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
!        WRITE(LUMOD3,'(A)')'!$OMP        BCC,AUXarrayCont,AUXarray2,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iP,iContB,iPrimB,TMP,iC) &'
        WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
        WRITE(LUMOD3,'(A)')'!$ACC        BCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'   do iP=1,nTUVP*nTUVQ*nPasses'
     WRITE(LUMOD3,'(A)')'    do iContB=1,nContB'
     WRITE(LUMOD3,'(A)')'     do iC=1,nContC*nContD*nContA'
     WRITE(LUMOD3,'(A)')'      AUXarrayCont(iC,iContB,iP)=0.0E0_realk'
     WRITE(LUMOD3,'(A)')'     enddo'
     WRITE(LUMOD3,'(A)')'     do iPrimB=1,nPrimB'
     WRITE(LUMOD3,'(A)')'      TMP = BCC(iPrimB,iContB)'
     WRITE(LUMOD3,'(A)')'      do iC=1,nContC*nContD*nContA'
     WRITE(LUMOD3,'(A)')'       AUXarrayCont(iC,iContB,iP)=AUXarrayCont(iC,iContB,iP)+TMP*AUXarray2(iC,iPrimB,iP)'
     WRITE(LUMOD3,'(A)')'      enddo'     
     WRITE(LUMOD3,'(A)')'     enddo'
     WRITE(LUMOD3,'(A)')'    enddo'
     WRITE(LUMOD3,'(A)')'   enddo'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionB'//ARCSTRING//'Gen')
     call writeString(LUMOD3)


     WRITE(LUMOD5,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionC'//ARCSTRING//'SegP')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD5)
     IF(CPU)THEN
        WRITE(LUMOD5,'(A)')'       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD5,'(A)')'       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,nTUVP,nTUVQ,iASync)'
     ENDIF

     WRITE(LUMOD5,'(A)')'    implicit none'
     WRITE(LUMOD5,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ,nTUVP,nTUVQ'
     WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(in) :: CCC(nPrimC,nContC)'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD*nPasses*nTUVP*nTUVQ)'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPasses*nTUVP*nTUVQ)'
     IF(.NOT.CPU)WRITE(LUMOD5,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD5,'(A)')'    !'
     WRITE(LUMOD5,'(A)')'    integer :: iP,iContC,iPrimC'
     WRITE(LUMOD5,'(A)')'    real(realk) :: TMP'
     IF(CPU)THEN
        WRITE(LUMOD5,'(A)')'!$OMP DO &'
        WRITE(LUMOD5,'(A)')'!$OMP PRIVATE(iP,iContC,iPrimC,TMP) '
!        WRITE(LUMOD5,'(A)')'!$OMP SHARED(nContC,nPrimC,nPasses,nPrimD,nTUVP,nTUVQ,&'
!        WRITE(LUMOD5,'(A)')'!$OMP        CCC,AUXarrayCont,AUXarray2,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD5,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD5,'(A)')'!$ACC PRIVATE(iP,iContC,iPrimC,TMP) &'
        WRITE(LUMOD5,'(A)')'!$ACC PRESENT(nContC,nPrimC,nPasses,nPrimD,&'
        WRITE(LUMOD5,'(A)')'!$ACC        CCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD5,'(A)')'    do iP = 1,nPasses*nTUVP*nTUVQ*nPrimD'
     WRITE(LUMOD5,'(A)')'     do iContC=1,nContC'
     WRITE(LUMOD5,'(A)')'      tmp = 0.0E0_realk'
     WRITE(LUMOD5,'(A)')'      do iPrimC=1,nPrimC'
     WRITE(LUMOD5,'(A)')'       tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)'
     WRITE(LUMOD5,'(A)')'      enddo'
     WRITE(LUMOD5,'(A)')'      AUXarrayCont(iContC,iP) = tmp'  
     WRITE(LUMOD5,'(A)')'     enddo'
     WRITE(LUMOD5,'(A)')'    enddo'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionC'//ARCSTRING//'SegP')
     call writeString(LUMOD5)
     WRITE(LUMOD5,'(A)')''

     call initString(1)
     call AddToString('  subroutine PrimitiveContractionD'//ARCSTRING//'SegP')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD5)
     IF(CPU)THEN
        WRITE(LUMOD5,'(A)')'       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD5,'(A)')'       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,nTUVP,nTUVQ,iASync)'
     ENDIF
     WRITE(LUMOD5,'(A)')'    implicit none'
     WRITE(LUMOD5,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ,nTUVP,nTUVQ'
     WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(in) :: DCC(nPrimD,nContD)'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(in) :: AUXarray2(nContC,nPrimD,nPasses*nTUVP*nTUVQ)'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPasses*nTUVP*nTUVQ)'
     IF(.NOT.CPU)WRITE(LUMOD5,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD5,'(A)')'    !'
     WRITE(LUMOD5,'(A)')'    integer :: iP,iContC,iContD,iPrimD'
     WRITE(LUMOD5,'(A)')'    real(realk) :: TMP'
     IF(CPU)THEN
        WRITE(LUMOD5,'(A)')'!$OMP DO &'
        WRITE(LUMOD5,'(A)')'!$OMP PRIVATE(iP,iContD,iPrimD,iContC,TMP) '
!        WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContC,nPrimD,nPasses,nContD,nTUVP,nTUVQ,&'
!        WRITE(LUMOD3,'(A)')'!$OMP        DCC,AUXarrayCont,AUXarray2)'
     ELSE
        WRITE(LUMOD5,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD5,'(A)')'!$ACC PRIVATE(iP,iContD,iPrimD,iContC,TMP) &'
        WRITE(LUMOD5,'(A)')'!$ACC PRESENT(nContC,nPrimD,nPasses,nContD,&'
        WRITE(LUMOD5,'(A)')'!$ACC         DCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD5,'(A)')'    do iP = 1,nPasses*nTUVP*nTUVQ'
     WRITE(LUMOD5,'(A)')'     do iContD=1,nContD'
     WRITE(LUMOD5,'(A)')'      do iContC=1,nContC'
     WRITE(LUMOD5,'(A)')'       AUXarrayCont(iContC,iContD,iP) = 0.0E0_realk'
     WRITE(LUMOD5,'(A)')'      enddo'
     WRITE(LUMOD5,'(A)')'      do iPrimD=1,nPrimD'
     WRITE(LUMOD5,'(A)')'       tmp = DCC(nPrimD,nContD)'
     WRITE(LUMOD5,'(A)')'       do iContC=1,nContC'
     WRITE(LUMOD5,'(A)')'        AUXarrayCont(iContC,iContD,iP) = AUXarrayCont(iContC,iContD,iP) + TMP*AUXarray2(iContC,iPrimD,iP)'
     WRITE(LUMOD5,'(A)')'       enddo'
     WRITE(LUMOD5,'(A)')'      enddo'
     WRITE(LUMOD5,'(A)')'     enddo'
     WRITE(LUMOD5,'(A)')'    enddo'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionD'//ARCSTRING//'SegP')
     call writeString(LUMOD5)


     WRITE(LUMOD4,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionA'//ARCSTRING//'SegQ')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD4)
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'       & nContP,ACC,nPrimA,nContA,nPrimB,nContB,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD4,'(A)')'       & nContP,ACC,nPrimA,nContA,nPrimB,nContB,nTUVP,nTUVQ,iASync)'
     ENDIF

     WRITE(LUMOD4,'(A)')'    implicit none'
     WRITE(LUMOD4,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nTUVP,nTUVQ'
     WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(in) :: ACC(nPrimA,nContA)'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(in) :: AUXarray2(nPrimA,nPrimB*nPasses*nTUVP*nTUVQ)'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContA,nPrimB*nPasses*nTUVP*nTUVQ)'
     IF(.NOT.CPU)WRITE(LUMOD4,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD4,'(A)')'    !'
     WRITE(LUMOD4,'(A)')'    integer :: iP,iContA,iPrimA'
     WRITE(LUMOD4,'(A)')'    real(realk) :: TMP'
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'!$OMP DO &'
        WRITE(LUMOD4,'(A)')'!$OMP PRIVATE(iP,iContA,iPrimA,TMP) '
!        WRITE(LUMOD3,'(A)')'!$OMP SHARED(nPasses,nTUVP,nTUVQ,nPrimA,nContA,nPrimB,&'
!        WRITE(LUMOD3,'(A)')'!$OMP        ACC,AUXarrayCont,AUXarray2,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD4,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD4,'(A)')'!$ACC PRIVATE(iP,iContA,iPrimA,TMP) &'
        WRITE(LUMOD4,'(A)')'!$ACC PRESENT(nPasses,nPrimA,nContA,nPrimB,&'
        WRITE(LUMOD4,'(A)')'!$ACC         ACC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD4,'(A)')'    do iP=1,nPasses*nTUVP*nTUVQ*nPrimB'
     WRITE(LUMOD4,'(A)')'      do iContA=1,nContA'
     WRITE(LUMOD4,'(A)')'        TMP = 0.0E0_realk'
     WRITE(LUMOD4,'(A)')'        do iPrimA=1,nPrimA'
     WRITE(LUMOD4,'(A)')'         TMP = TMP + ACC(iPrimA,iContA)*AUXarray2(iPrimA,iP)'
     WRITE(LUMOD4,'(A)')'        enddo'
     WRITE(LUMOD4,'(A)')'        AUXarrayCont(iContA,iP) = TMP'
     WRITE(LUMOD4,'(A)')'      enddo'
     WRITE(LUMOD4,'(A)')'    enddo'
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionA'//ARCSTRING//'SegQ')
     call writeString(LUMOD4)
     WRITE(LUMOD4,'(A)')''

     WRITE(LUMOD4,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionB'//ARCSTRING//'SegQ')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD4)
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'       & nContP,BCC,nPrimA,nContA,nPrimB,nContB,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD4,'(A)')'       & nContP,BCC,nPrimA,nContA,nPrimB,nContB,nTUVP,nTUVQ,iASync)'
     ENDIF

     WRITE(LUMOD4,'(A)')'    implicit none'
     WRITE(LUMOD4,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nTUVP,nTUVQ'
     WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(in) :: BCC(nPrimB,nContB)'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(in) :: AUXarray2(nContA,nPrimB,nPasses*nTUVP*nTUVQ)'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContA,nContB,nPasses*nTUVP*nTUVQ)'
     IF(.NOT.CPU)WRITE(LUMOD4,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD4,'(A)')'    !'
     WRITE(LUMOD4,'(A)')'    integer :: iP,iContB,iContA,iPrimB'
     WRITE(LUMOD4,'(A)')'    real(realk) :: TMP'
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'!$OMP DO &'
        WRITE(LUMOD4,'(A)')'!$OMP PRIVATE(iP,iContA,iContB,iPrimB,TMP) '
!        WRITE(LUMOD3,'(A)')'!$OMP SHARED(nPasses,nTUVP,nTUVQ,nPrimB,nContB,nContA,&'
!        WRITE(LUMOD3,'(A)')'!$OMP        BCC,AUXarrayCont,AUXarray2,nTUVP,nTUVQ)'
     ELSE
        WRITE(LUMOD4,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD4,'(A)')'!$ACC PRIVATE(iP,iContA,iContB,iPrimB,TMP) &'
        WRITE(LUMOD4,'(A)')'!$ACC PRESENT(nPasses,nPrimB,nContB,nContA,&'
        WRITE(LUMOD4,'(A)')'!$ACC         BCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD4,'(A)')'    do iP=1,nPasses*nTUVP*nTUVQ'
     WRITE(LUMOD4,'(A)')'     do iContB=1,nContB'
     WRITE(LUMOD4,'(A)')'      do iContA=1,nContA'
     WRITE(LUMOD4,'(A)')'       AUXarrayCont(iContA,iContB,iP) = 0.0E0_realk'
     WRITE(LUMOD4,'(A)')'      enddo'
     WRITE(LUMOD4,'(A)')'      do iPrimB=1,nPrimB'
     WRITE(LUMOD4,'(A)')'       TMP = BCC(iPrimB,iContB)'
     WRITE(LUMOD4,'(A)')'       do iContA=1,nContA'
     WRITE(LUMOD4,'(A)')'        AUXarrayCont(iContA,iContB,iP) = AUXarrayCont(iContA,iContB,iP) + TMP*AUXarray2(iContA,iPrimB,iP)'
     WRITE(LUMOD4,'(A)')'       enddo'
     WRITE(LUMOD4,'(A)')'      enddo'
     WRITE(LUMOD4,'(A)')'     enddo'
     WRITE(LUMOD4,'(A)')'    enddo'
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionB'//ARCSTRING//'SegQ')
     call writeString(LUMOD4)
     WRITE(LUMOD4,'(A)')''
  ELSE

!Special for nTUVP*nTUVQ=1

     WRITE(LUMOD3,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionC'//ARCSTRING//'Gen1')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD3)
     WRITE(LUMOD3,'(A)')'       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD)'
     ELSE
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,iASync)'
     ENDIF

     WRITE(LUMOD3,'(A)')'    implicit none'
     WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ'
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: CCC(nPrimC,nContC)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD*nPrimA*nPrimB*nPasses)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPrimA*nPrimB*nPasses)'
     IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD3,'(A)')'    !'
     WRITE(LUMOD3,'(A)')'    integer :: iP,iContC,iPrimC'
     WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP DO &'
        WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iP,iPrimC,iContC,TMP) '
!        WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
!        WRITE(LUMOD3,'(A)')'!$OMP        CCC,AUXarrayCont,AUXarray2)'
     ELSE
        WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iP,iPrimC,iContC,TMP) &'
        WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
        WRITE(LUMOD3,'(A)')'!$ACC         CCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses'
     WRITE(LUMOD3,'(A)')'     do iContC=1,nContC'
     WRITE(LUMOD3,'(A)')'      TMP = 0.0E0_realk'
     WRITE(LUMOD3,'(A)')'      do iPrimC=1,nPrimC'
     WRITE(LUMOD3,'(A)')'       TMP = TMP + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)'
     WRITE(LUMOD3,'(A)')'      enddo'
     WRITE(LUMOD3,'(A)')'      AUXarrayCont(iContC,iP)=TMP'
     WRITE(LUMOD3,'(A)')'     enddo'
     WRITE(LUMOD3,'(A)')'    enddo'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionC'//ARCSTRING//'Gen1')
     call writeString(LUMOD3)
          
     WRITE(LUMOD3,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionD'//ARCSTRING//'Gen1')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD3)
     WRITE(LUMOD3,'(A)')'       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD)'
     ELSE
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'    implicit none'
     WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ'
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: DCC(nPrimD,nContD)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray2(nContC,nPrimD,nPrimA*nPrimB*nPasses)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPrimA*nPrimB*nPasses)'
     IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD3,'(A)')'    !'
     WRITE(LUMOD3,'(A)')'    integer :: iP,iContD,iContC,iPrimD'
     WRITE(LUMOD3,'(A)')'    integer :: iTUV'
     WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'            
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP DO &'
        WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iP,iContD,iContC,iPrimD,TMP) '
!        WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContC,nContD,nPasses,nPrimD,nPrimA,nPrimB,&'
!        WRITE(LUMOD3,'(A)')'!$OMP        DCC,AUXarrayCont,AUXarray2)'
     ELSE
        WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iP,iContD,iContC,iPrimD,TMP) &'
        WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nContC,nContD,nPasses,nPrimD,nPrimA,nPrimB,&'
        WRITE(LUMOD3,'(A)')'!$ACC        DCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'    do iP = 1,nPrimA*nPrimB*nPasses'
     WRITE(LUMOD3,'(A)')'     do iContD=1,nContD'
     WRITE(LUMOD3,'(A)')'      do iContC=1,nContC'
     WRITE(LUMOD3,'(A)')'       AUXarrayCont(iContC,iContD,iP) = 0.0E0_realk'
     WRITE(LUMOD3,'(A)')'      enddo'
     WRITE(LUMOD3,'(A)')'      do iPrimD=1,nPrimD'
     WRITE(LUMOD3,'(A)')'       TMP = DCC(iPrimD,iContD)'
     WRITE(LUMOD3,'(A)')'       do iContC=1,nContC'
     WRITE(LUMOD3,'(A)')'        AUXarrayCont(iContC,iContD,iP)=AUXarrayCont(iContC,iContD,iP)+TMP*AUXarray2(iContC,iPrimD,iP)'
     WRITE(LUMOD3,'(A)')'       enddo'
     WRITE(LUMOD3,'(A)')'      enddo'
     WRITE(LUMOD3,'(A)')'     enddo'
     WRITE(LUMOD3,'(A)')'    enddo'               
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionD'//ARCSTRING//'Gen1')
     call writeString(LUMOD3)
          
     WRITE(LUMOD3,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionA'//ARCSTRING//'Gen1')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD3)
     WRITE(LUMOD3,'(A)')'       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD)'
     ELSE
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,iASync)'
     ENDIF

     WRITE(LUMOD3,'(A)')'    implicit none'
     WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ'
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: ACC(nPrimA,nContA)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray2(nContC*nContD,nPrimA,nPrimB*nPasses)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC*nContD,nContA,nPrimB*nPasses)'
     IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD3,'(A)')'    !'
     WRITE(LUMOD3,'(A)')'    integer :: iC,iP,iContA,iPrimA'
     WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'            
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP DO &'
        WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iC,iP,iContA,iPrimA,TMP) '
!        WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContA,nContC,nContD,nPasses,nPrimA,nPrimB,&'
!        WRITE(LUMOD3,'(A)')'!$OMP        ACC,AUXarrayCont,AUXarray2)'
     ELSE
        WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iC,iP,iContA,iPrimD,TMP) &'
        WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nContA,nContC,nContD,nPasses,nPrimA,nPrimB,&'
        WRITE(LUMOD3,'(A)')'!$ACC        ACC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'    do iP = 1,nPrimB*nPasses'
     WRITE(LUMOD3,'(A)')'     do iContA=1,nContA'
     WRITE(LUMOD3,'(A)')'      do iC=1,nContC*nContD'
     WRITE(LUMOD3,'(A)')'       AUXarrayCont(iC,iContA,iP) = 0.0E0_realk'
     WRITE(LUMOD3,'(A)')'      enddo'
     WRITE(LUMOD3,'(A)')'      do iPrimA=1,nPrimA'
     WRITE(LUMOD3,'(A)')'       TMP = ACC(iPrimA,iContA)'
     WRITE(LUMOD3,'(A)')'       do iC=1,nContC*nContD'
     WRITE(LUMOD3,'(A)')'        AUXarrayCont(iC,iContA,iP)=AUXarrayCont(iC,iContA,iP)+TMP*AUXarray2(iC,iPrimA,iP)'
     WRITE(LUMOD3,'(A)')'       enddo'
     WRITE(LUMOD3,'(A)')'      enddo'
     WRITE(LUMOD3,'(A)')'     enddo'
     WRITE(LUMOD3,'(A)')'    enddo'               
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionA'//ARCSTRING//'Gen1')
     call writeString(LUMOD3)
          
     WRITE(LUMOD3,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionB'//ARCSTRING//'Gen1')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD3)
     WRITE(LUMOD3,'(A)')'       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD)'
     ELSE
        WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,iASync)'
     ENDIF

     WRITE(LUMOD3,'(A)')'    implicit none'
     WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ'
     WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: BCC(nPrimB,nContB)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: AUXarray2(nContC*nContD*nContA,nPrimB,nPasses)'
     WRITE(LUMOD3,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC*nContD*nContA,nContB,nPasses)'
     IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD3,'(A)')'    !'
     WRITE(LUMOD3,'(A)')'    integer :: iC,iP,iContB,iPrimB'
     WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'            
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP DO &'
        WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iC,iP,iContB,iPrimB,TMP) '
!        WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContB,nContA,nContC,nContD,nPasses,nPrimB,&'
!        WRITE(LUMOD3,'(A)')'!$OMP        BCC,AUXarrayCont,AUXarray2)'
     ELSE
        WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iC,iP,iContA,iPrimD,TMP) &'
        WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nCOntB,nContA,nContC,nContD,nPasses,nPrimB,&'
        WRITE(LUMOD3,'(A)')'!$ACC         BCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD3,'(A)')'    do iP = 1,nPasses'
     WRITE(LUMOD3,'(A)')'     do iContB=1,nContB'
     WRITE(LUMOD3,'(A)')'      do iC=1,nContC*nContD*nContA'
     WRITE(LUMOD3,'(A)')'       AUXarrayCont(iC,iContB,iP) = 0.0E0_realk'
     WRITE(LUMOD3,'(A)')'      enddo'
     WRITE(LUMOD3,'(A)')'      do iPrimB=1,nPrimB'
     WRITE(LUMOD3,'(A)')'       TMP = BCC(iPrimB,iContB)'
     WRITE(LUMOD3,'(A)')'       do iC=1,nContC*nContD*nContA'
     WRITE(LUMOD3,'(A)')'        AUXarrayCont(iC,iContB,iP)=AUXarrayCont(iC,iContB,iP)+TMP*AUXarray2(iC,iPrimB,iP)'
     WRITE(LUMOD3,'(A)')'       enddo'
     WRITE(LUMOD3,'(A)')'      enddo'
     WRITE(LUMOD3,'(A)')'     enddo'
     WRITE(LUMOD3,'(A)')'    enddo'               
     IF(CPU)THEN
        WRITE(LUMOD3,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionB'//ARCSTRING//'Gen1')
     call writeString(LUMOD3)

! PSEGMENTED
     WRITE(LUMOD5,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionC'//ARCSTRING//'SegP1')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD5)
     IF(CPU)THEN
        WRITE(LUMOD5,'(A)')'       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)'
     ELSE
        WRITE(LUMOD5,'(A)')'       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,iASync)'
     ENDIF

     WRITE(LUMOD5,'(A)')'    implicit none'
     WRITE(LUMOD5,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ'
     WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(in) :: CCC(nPrimC,nContC)'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(in) :: AUXarray2(nPrimC,nPrimD*nPasses)'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC,nPrimD*nPasses)'
     IF(.NOT.CPU)WRITE(LUMOD5,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD5,'(A)')'    !'
     WRITE(LUMOD5,'(A)')'    integer :: iP,iContC,iPrimC'
     WRITE(LUMOD5,'(A)')'    real(realk) :: TMP'
     IF(CPU)THEN
        WRITE(LUMOD5,'(A)')'!$OMP DO &'
        WRITE(LUMOD5,'(A)')'!$OMP PRIVATE(iP,iContC,iPrimC,TMP) '
!        WRITE(LUMOD5,'(A)')'!$OMP SHARED(nContC,nContD,nPasses,nPrimC,nPrimD,&'
!        WRITE(LUMOD5,'(A)')'!$OMP        CCC,AUXarrayCont,AUXarray2)'
     ELSE
        WRITE(LUMOD5,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD5,'(A)')'!$ACC PRIVATE(iP,iContC,iPrimC,TMP) &'
        WRITE(LUMOD5,'(A)')'!$ACC PRESENT(nContC,nContD,nPasses,nPrimC,nPrimD,&'
        WRITE(LUMOD5,'(A)')'!$ACC        CCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD5,'(A)')'    do iP = 1,nPasses*nPrimD'
     WRITE(LUMOD5,'(A)')'     do iContC=1,nContC'
     WRITE(LUMOD5,'(A)')'      tmp = 0.0E0_realk'
     WRITE(LUMOD5,'(A)')'      do iPrimC=1,nPrimC'
     WRITE(LUMOD5,'(A)')'       tmp = tmp + CCC(iPrimC,iContC)*AUXarray2(iPrimC,iP)'
     WRITE(LUMOD5,'(A)')'      enddo'
     WRITE(LUMOD5,'(A)')'      AUXarrayCont(iContC,iP) = tmp'
     WRITE(LUMOD5,'(A)')'     enddo'
     WRITE(LUMOD5,'(A)')'    enddo'
     IF(CPU)THEN
        WRITE(LUMOD5,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionC'//ARCSTRING//'SegP1')
     call writeString(LUMOD5)

     WRITE(LUMOD5,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionD'//ARCSTRING//'SegP1')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD5)
     IF(CPU)THEN
        WRITE(LUMOD5,'(A)')'       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)'
     ELSE
        WRITE(LUMOD5,'(A)')'       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,iASync)'
     ENDIF

     WRITE(LUMOD5,'(A)')'    implicit none'
     WRITE(LUMOD5,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ'
     WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(in) :: DCC(nPrimD,nContD)'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(in) :: AUXarray2(nContC,nPrimD,nPasses)'
     WRITE(LUMOD5,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContC,nContD,nPasses)'
     IF(.NOT.CPU)WRITE(LUMOD5,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD5,'(A)')'    !'
     WRITE(LUMOD5,'(A)')'    integer :: iPassP,iContD,iPrimD,iContC'
     WRITE(LUMOD5,'(A)')'    real(realk) :: TMP'
     IF(CPU)THEN
        WRITE(LUMOD5,'(A)')'!$OMP DO &'
        WRITE(LUMOD5,'(A)')'!$OMP PRIVATE(iPassP,iContD,iPrimD,iContC,TMP) '
!        WRITE(LUMOD5,'(A)')'!$OMP SHARED(nContC,nContD,nPasses,nPrimC,nPrimD,&'
!        WRITE(LUMOD5,'(A)')'!$OMP        DCC,AUXarrayCont,AUXarray2)'
     ELSE
        WRITE(LUMOD5,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD5,'(A)')'!$ACC PRIVATE(iPassP,iContD,iPrimD,iContC,TMP) &'
        WRITE(LUMOD5,'(A)')'!$ACC PRESENT(nContC,nContD,nPasses,nPrimD,&'
        WRITE(LUMOD5,'(A)')'!$ACC        DCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD5,'(A)')'    do iPassP = 1,nPasses'
     WRITE(LUMOD5,'(A)')'     do iContD=1,nContD'
     WRITE(LUMOD5,'(A)')'      do iContC=1,nContC'
     WRITE(LUMOD5,'(A)')'        AUXarrayCont(iContC,iContD,iPassP) = 0.0E0_realk'
     WRITE(LUMOD5,'(A)')'       enddo'
     WRITE(LUMOD5,'(A)')'      do iPrimD=1,nPrimD'
     WRITE(LUMOD5,'(A)')'       tmp = DCC(iPrimD,iContD)'
     WRITE(LUMOD5,'(A)')'       do iContC=1,nContC'
     WRITE(LUMOD5,'(A)')'        AUXarrayCont(iContC,iContD,iPassP) = AUXarrayCont(iContC,iContD,iPassP) + tmp*AUXarray2(iContC,iPrimD,iPassP)'
     WRITE(LUMOD5,'(A)')'       enddo'
     WRITE(LUMOD5,'(A)')'      enddo'
     WRITE(LUMOD5,'(A)')'     enddo'
     WRITE(LUMOD5,'(A)')'    enddo'
     IF(CPU)THEN
        WRITE(LUMOD5,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionD'//ARCSTRING//'SegP1')
     call writeString(LUMOD5)
     
     ! QSEGMENTED
     WRITE(LUMOD4,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionA'//ARCSTRING//'SegQ1')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD4)
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)'
     ELSE
        WRITE(LUMOD4,'(A)')'       & nContP,ACC,nPrimA,nContA,nPrimB,nContB,iASync)'
     ENDIF

     WRITE(LUMOD4,'(A)')'    implicit none'
     WRITE(LUMOD4,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP'
     WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(in) :: ACC(nPrimA,nContA)'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(in) :: AUXarray2(nPrimA,nPrimB*nPasses)'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContA,nPrimB*nPasses)'
     IF(.NOT.CPU)WRITE(LUMOD4,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD4,'(A)')'    !'
     WRITE(LUMOD4,'(A)')'    integer :: iP,iContA,iPrimA'
     WRITE(LUMOD4,'(A)')'    real(realk) :: TMP'
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'!$OMP DO &'
        WRITE(LUMOD4,'(A)')'!$OMP PRIVATE(iP,iContA,iPrimA,TMP) '
!        WRITE(LUMOD4,'(A)')'!$OMP SHARED(nContA,nPasses,nPrimA,nPrimB,&'
!        WRITE(LUMOD4,'(A)')'!$OMP        ACC,AUXarrayCont,AUXarray2)'
     ELSE
        WRITE(LUMOD4,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD4,'(A)')'!$ACC PRIVATE(iP,iContA,iPrimA,TMP) &'
        WRITE(LUMOD4,'(A)')'!$ACC PRESENT(nContA,nPasses,nPrimA,nPrimB,&'
        WRITE(LUMOD4,'(A)')'!$ACC         ACC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD4,'(A)')'   do iP = 1,nPasses*nPrimB'
     WRITE(LUMOD4,'(A)')'    do iContA=1,nContA'
     WRITE(LUMOD4,'(A)')'     TMP = 0.0E0_realk'
     WRITE(LUMOD4,'(A)')'     do iPrimA=1,nPrimA'
     WRITE(LUMOD4,'(A)')'      TMP = TMP + ACC(iPrimA,iContA)*AUXarray2(iPrimA,iP)'
     WRITE(LUMOD4,'(A)')'     enddo'
     WRITE(LUMOD4,'(A)')'     AUXarrayCont(iContA,iP) = tmp'
     WRITE(LUMOD4,'(A)')'    enddo'
     WRITE(LUMOD4,'(A)')'   enddo'
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionA'//ARCSTRING//'SegQ1')
     call writeString(LUMOD4)
     WRITE(LUMOD4,'(A)')''

     WRITE(LUMOD4,'(A)')''
     call initString(1)
     call AddToString('  subroutine PrimitiveContractionB'//ARCSTRING//'SegQ1')
     call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
     call writeString(LUMOD4)
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)'
     ELSE
        WRITE(LUMOD4,'(A)')'       & nContP,BCC,nPrimA,nContA,nPrimB,nContB,iASync)'
     ENDIF

     WRITE(LUMOD4,'(A)')'    implicit none'
     WRITE(LUMOD4,'(A)')'    !Warning Primitive screening modifies this!!! '
     WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP'
     WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(in) :: BCC(nPrimB,nContB)'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(in) :: AUXarray2(nContA,nPrimB,nPasses)'
     WRITE(LUMOD4,'(A)')'    real(realk),intent(inout) :: AUXarrayCont(nContA,nContB,nPasses)'
     IF(.NOT.CPU)WRITE(LUMOD4,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
     WRITE(LUMOD4,'(A)')'    !'
     WRITE(LUMOD4,'(A)')'    integer :: iPassP,iContB,iPrimB,iContA'
     WRITE(LUMOD4,'(A)')'    real(realk) :: TMP'
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'!$OMP DO &'
        WRITE(LUMOD4,'(A)')'!$OMP PRIVATE(iContA,iPassP,iContB,iPrimB,TMP) '
!        WRITE(LUMOD4,'(A)')'!$OMP SHARED(nContA,nContB,nPasses,nPrimB,&'
!        WRITE(LUMOD4,'(A)')'!$OMP        BCC,AUXarrayCont,AUXarray2)'
     ELSE
        WRITE(LUMOD4,'(A)')'!$ACC PARALLEL LOOP &'
        WRITE(LUMOD4,'(A)')'!$ACC PRIVATE(iContA,iPassP,iContB,iPrimB,TMP) &'
        WRITE(LUMOD4,'(A)')'!$ACC PRESENT(nContA,nContB,nPasses,nPrimB,&'
        WRITE(LUMOD4,'(A)')'!$ACC         BCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
     ENDIF
     WRITE(LUMOD4,'(A)')'   do iPassP = 1,nPasses'
     WRITE(LUMOD4,'(A)')'    do iContB=1,nContB'
     WRITE(LUMOD4,'(A)')'     do iContA=1,nContA'
     WRITE(LUMOD4,'(A)')'       AUXarrayCont(iContA,iContB,iPassP)=0.0E0_realk'
     WRITE(LUMOD4,'(A)')'     enddo'
     WRITE(LUMOD4,'(A)')'     do iPrimB=1,nPrimB'
     WRITE(LUMOD4,'(A)')'      TMP = BCC(iPrimB,iContB)'
     WRITE(LUMOD4,'(A)')'      do iContA=1,nContA'
     WRITE(LUMOD4,'(A)')'       AUXarrayCont(iContA,iContB,iPassP)=AUXarrayCont(iContA,iContB,iPassP)+tmp*AUXarray2(iContA,iPrimB,iPassP)'
     WRITE(LUMOD4,'(A)')'      enddo'
     WRITE(LUMOD4,'(A)')'     enddo'
     WRITE(LUMOD4,'(A)')'    enddo'
     WRITE(LUMOD4,'(A)')'   enddo'
     IF(CPU)THEN
        WRITE(LUMOD4,'(A)')'!$OMP END DO'
     ENDIF
     call initString(1)
     call AddToString('  end subroutine PrimitiveContractionB'//ARCSTRING//'SegQ1')
     call writeString(LUMOD4)
     WRITE(LUMOD4,'(A)')''
     
     allocate(UniquenTUVs(3*3*3*3))
     UniquenTUVs(1) = 1
     nUniquenTUVs = 1
     DO AngmomA = 0,2
      DO AngmomB = 0,2
       DO AngmomC = 0,2
        DO AngmomD = 0,2
         AngmomP = AngmomA + AngmomB
         AngmomQ = AngmomC + AngmomD
         nTUV = (AngmomPQ+1)*(AngmomPQ+2)*(AngmomPQ+3)/6
         nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
         nTUVQ = (AngmomQ+1)*(AngmomQ+2)*(AngmomQ+3)/6
         UNIQUE = .TRUE.
         DO I=1,nUniquenTUVs
            IF(nTUVP*nTUVQ.EQ.UniquenTUVs(I))THEN
               UNIQUE = .FALSE.
            ENDIF
         ENDDO
         IF(UNIQUE)THEN
          nUniquenTUVs = nUniquenTUVs + 1
          UniquenTUVs(nUniquenTUVs) = nTUVP*nTUVQ
!GENERAL 

          WRITE(LUMOD3,'(A)')''
          call initString(1)
          call AddToString('  subroutine PrimitiveContractionC'//ARCSTRING//'Gen')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
          call writeString(LUMOD3)
          WRITE(LUMOD3,'(A)')'       & nContP,nContQ,CCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'       & nPrimD,nContD)'
          ELSE
             WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,iASync)'
          ENDIF
          WRITE(LUMOD3,'(A)')'    implicit none'
          WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
          WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ'
          WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
          WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: CCC(nPrimC,nContC)'
          WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(in) :: AUXarray2(',nTUVP*nTUVQ,',nPrimC,nPrimD*nPrimA*nPrimB*nPasses)'
          WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(inout) :: AUXarrayCont(',nTUVP*nTUVQ,',nContC,nPrimD*nPrimA*nPrimB*nPasses)'
          IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
          WRITE(LUMOD3,'(A)')'    !'
          WRITE(LUMOD3,'(A)')'    integer :: iP,iPrimC,iContC,iTUV'
          WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'
          WRITE(LUMOD3,'(A)')'!Scaling p**4*c*nTUV*nPassQ: nPrimA*nPrimB*nPrimC*nPrimD*nContC*nTUV*nPassQ'
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'!$OMP DO &'
             WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iTUV,iP,iPrimC,iContC,TMP) '
!             WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
!             WRITE(LUMOD3,'(A)')'!$OMP        CCC,AUXarrayCont,AUXarray2)'
          ELSE
             WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
             WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iTUV,iP,iPrimC,iContC,TMP) &'
             WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nContC,nPasses,nPrimC,nPrimD,nPrimA,nPrimB,&'
             WRITE(LUMOD3,'(A)')'!$ACC        CCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
          ENDIF
          WRITE(LUMOD3,'(A)')'    do iP = 1,nPrimD*nPrimA*nPrimB*nPasses'
          WRITE(LUMOD3,'(A)')'     do iContC=1,nContC'
          WRITE(LUMOD3,'(A,I5)')'      do iTUV=1,',nTUVP*nTUVQ
          WRITE(LUMOD3,'(A)')'       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk'
          WRITE(LUMOD3,'(A)')'      enddo'
          WRITE(LUMOD3,'(A)')'      do iPrimC=1,nPrimC'
          WRITE(LUMOD3,'(A)')'       TMP = CCC(iPrimC,iContC)'
          WRITE(LUMOD3,'(A,I5)')'       do iTUV=1,',nTUVP*nTUVQ
          WRITE(LUMOD3,'(A)')'        AUXarrayCont(iTUV,iContC,iP)=AUXarrayCont(iTUV,iContC,iP)+TMP*AUXarray2(iTUV,iPrimC,iP)'
          WRITE(LUMOD3,'(A)')'       enddo'
          WRITE(LUMOD3,'(A)')'      enddo'
          WRITE(LUMOD3,'(A)')'     enddo'
          WRITE(LUMOD3,'(A)')'    enddo'
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'!$OMP END DO'
          ENDIF
          call initString(1)
          call AddToString('  end subroutine PrimitiveContractionC'//ARCSTRING//'Gen')
          call AddToString(nTUVP*nTUVQ)
          call writeString(LUMOD3)
          
          WRITE(LUMOD3,'(A)')''
          call initString(1)
          call AddToString('  subroutine PrimitiveContractionD'//ARCSTRING//'Gen')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
          call writeString(LUMOD3)
          WRITE(LUMOD3,'(A)')'       & nContP,nContQ,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'       & nPrimD,nContD)'
          ELSE
             WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,iASync)'
          ENDIF
          WRITE(LUMOD3,'(A)')'    implicit none'
          WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
          WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ'
          WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
          WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: DCC(nPrimD,nContD)'
          WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(in) :: AUXarray2(',nTUVP*nTUVQ,'*nContC,nPrimD,nPrimA*nPrimB*nPasses)'
          WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(inout) :: AUXarrayCont(',nTUVP*nTUVQ,'*nContC,nContD,nPrimA*nPrimB*nPasses)'
          IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
          WRITE(LUMOD3,'(A)')'    !'
          WRITE(LUMOD3,'(A)')'    integer :: iTUV,iP,iContD,iPrimD'
          WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'            
          WRITE(LUMOD3,'(A)')'!Scaling  p**3*c**2*nTUV*nPassQ: nPrimA*nPrimB*nPrimD*nContC*nContD*nTUV*nPassQ'
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'!$OMP DO &'
             WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iTUV,iP,iContD,iPrimD,TMP) '
!             WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContC,nContD,nPasses,nPrimD,nPrimA,nPrimB,&'
!             WRITE(LUMOD3,'(A)')'!$OMP        DCC,AUXarrayCont,AUXarray2)'
          ELSE
             WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
             WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iTUV,iP,iContD,iPrimD,TMP) &'
             WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nContC,nContD,nPasses,nPrimD,nPrimA,nPrimB,&'
             WRITE(LUMOD3,'(A)')'!$ACC        DCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
          ENDIF
          WRITE(LUMOD3,'(A)')'    do iP = 1,nPrimA*nPrimB*nPasses'
          WRITE(LUMOD3,'(A)')'     do iContD=1,nContD'
          call initString(7)
          call AddToString('do iTUV=1,')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('*nContC')
          call writeString(LUMOD3)               
          WRITE(LUMOD3,'(A)')'       AUXarrayCont(iTUV,iContD,iP) = 0.0E0_realk'
          WRITE(LUMOD3,'(A)')'      enddo'
          WRITE(LUMOD3,'(A)')'      do iPrimD=1,nPrimD'
          WRITE(LUMOD3,'(A)')'       TMP = DCC(iPrimD,iContD)'
          call initString(8)
          call AddToString('do iTUV=1,')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('*nContC')
          call writeString(LUMOD3)               
          WRITE(LUMOD3,'(A)')'        AUXarrayCont(iTUV,iContD,iP)=AUXarrayCont(iTUV,iContD,iP)+TMP*AUXarray2(iTUV,iPrimD,iP)'
          WRITE(LUMOD3,'(A)')'       enddo'
          WRITE(LUMOD3,'(A)')'      enddo'
          WRITE(LUMOD3,'(A)')'     enddo'
          WRITE(LUMOD3,'(A)')'    enddo'               
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'!$OMP END DO'
          ENDIF
          call initString(1)
          call AddToString('  end subroutine PrimitiveContractionD'//ARCSTRING//'Gen')
          call AddToString(nTUVP*nTUVQ)
          call writeString(LUMOD3)
          
          WRITE(LUMOD3,'(A)')''
          call initString(1)
          call AddToString('  subroutine PrimitiveContractionA'//ARCSTRING//'Gen')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
          call writeString(LUMOD3)
          WRITE(LUMOD3,'(A)')'       & nContP,nContQ,ACC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'       & nPrimD,nContD)'
          ELSE
             WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,iASync)'
          ENDIF
          WRITE(LUMOD3,'(A)')'    implicit none'
          WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
          WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ'
          WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
          WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: ACC(nPrimA,nContA)'
          WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(in) :: AUXarray2(',nTUVP*nTUVQ,'*nContC*nContD,nPrimA,nPrimB*nPasses)'
          WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(inout) :: AUXarrayCont(',nTUVP*nTUVQ,'*nContC*nContD,nContA,nPrimB*nPasses)'
          IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
          WRITE(LUMOD3,'(A)')'    !'
          WRITE(LUMOD3,'(A)')'    integer :: iTUV,iP,iContA,iPrimA'
          WRITE(LUMOD3,'(A)')'    real(realk) :: TMP'            
          WRITE(LUMOD3,'(A)')'!Scaling  p**2*c**3*nTUV*nPassQ: nPrimB*nPrimA*nContC*nContD*nContA*nTUV*nPassQ'
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'!$OMP DO &'
             WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) '
!             WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContA,nContC,nContD,nPasses,nPrimA,nPrimB,&'
!             WRITE(LUMOD3,'(A)')'!$OMP        ACC,AUXarrayCont,AUXarray2)'
          ELSE
             WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
             WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iTUV,iP,iContA,iPrimA,TMP) &'
             WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nContA,nContC,nContD,nPasses,nPrimA,nPrimB,&'
             WRITE(LUMOD3,'(A)')'!$ACC        ACC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
          ENDIF
          WRITE(LUMOD3,'(A)')'    do iP = 1,nPrimB*nPasses'
          WRITE(LUMOD3,'(A)')'     do iContA=1,nContA'
          call initString(7)
          call AddToString('do iTUV=1,')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('*nContC*nContD')
          call writeString(LUMOD3)               
          WRITE(LUMOD3,'(A)')'       AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk'
          WRITE(LUMOD3,'(A)')'      enddo'
          WRITE(LUMOD3,'(A)')'      do iPrimA=1,nPrimA'
          WRITE(LUMOD3,'(A)')'       TMP = ACC(iPrimA,iContA)'
          call initString(8)
          call AddToString('do iTUV=1,')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('*nContC*nContD')
          call writeString(LUMOD3)               
          WRITE(LUMOD3,'(A)')'        AUXarrayCont(iTUV,iContA,iP)=AUXarrayCont(iTUV,iContA,iP)+TMP*AUXarray2(iTUV,iPrimA,iP)'
          WRITE(LUMOD3,'(A)')'       enddo'
          WRITE(LUMOD3,'(A)')'      enddo'
          WRITE(LUMOD3,'(A)')'     enddo'
          WRITE(LUMOD3,'(A)')'    enddo'               
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'!$OMP END DO'
          ENDIF
          call initString(1)
          call AddToString('  end subroutine PrimitiveContractionA'//ARCSTRING//'Gen')
          call AddToString(nTUVP*nTUVQ)
          call writeString(LUMOD3)
          
          WRITE(LUMOD3,'(A)')''
          call initString(1)
          call AddToString('  subroutine PrimitiveContractionB'//ARCSTRING//'Gen')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
          call writeString(LUMOD3)
          WRITE(LUMOD3,'(A)')'       & nContP,nContQ,BCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,&'
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'       & nPrimD,nContD)'
          ELSE
             WRITE(LUMOD3,'(A)')'       & nPrimD,nContD,iASync)'
          ENDIF
          WRITE(LUMOD3,'(A)')'    implicit none'
          WRITE(LUMOD3,'(A)')'    !Warning Primitive screening modifies this!!! '
          WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP,nContQ'
          WRITE(LUMOD3,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD'
          WRITE(LUMOD3,'(A)')'    real(realk),intent(in) :: BCC(nPrimB,nContB)'
          WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(in) :: AUXarray2(',nTUVP*nTUVQ,'*nContC*nContD*nContA,nPrimB,nPasses)'
          WRITE(LUMOD3,'(A,I5,A)')'    real(realk),intent(inout) :: AUXarrayCont(',nTUVP*nTUVQ,'*nContC*nContD*nContA,nContB,nPasses)'
          IF(.NOT.CPU)WRITE(LUMOD3,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
          WRITE(LUMOD3,'(A)')'    !'
          WRITE(LUMOD3,'(A)')'    integer :: iTUV,iP,iContB,iPrimB'
          WRITE(LUMOD3,'(A,I5,A)')'    real(realk) :: TMP'            
          WRITE(LUMOD3,'(A)')'!Scaling  p*c**4*nTUV*nPassQ: nPrimB*nContA*nContC*nContD*nContB*nTUV*nPassQ'
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'!$OMP DO &'
             WRITE(LUMOD3,'(A)')'!$OMP PRIVATE(iTUV,iP,iContB,iPrimB,TMP) '
!             WRITE(LUMOD3,'(A)')'!$OMP SHARED(nContB,nContA,nContC,nContD,nPasses,nPrimB,&'
!             WRITE(LUMOD3,'(A)')'!$OMP        BCC,AUXarrayCont,AUXarray2)'
          ELSE
             WRITE(LUMOD3,'(A)')'!$ACC PARALLEL LOOP &'
             WRITE(LUMOD3,'(A)')'!$ACC PRIVATE(iTUV,iP,iContB,iPrimB,TMP) &'
             WRITE(LUMOD3,'(A)')'!$ACC PRESENT(nCOntB,nContA,nContC,nContD,nPasses,nPrimB,&'
             WRITE(LUMOD3,'(A)')'!$ACC        BCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
          ENDIF
          WRITE(LUMOD3,'(A)')'    do iP = 1,nPasses'
          WRITE(LUMOD3,'(A)')'     do iContB=1,nContB'
          call initString(7)
          call AddToString('do iTUV=1,')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('*nContC*nContD*nContA')
          call writeString(LUMOD3)               
          WRITE(LUMOD3,'(A)')'       AUXarrayCont(iTUV,iContB,iP) = 0.0E0_realk'
          WRITE(LUMOD3,'(A)')'      enddo'
          WRITE(LUMOD3,'(A)')'      do iPrimB=1,nPrimB'
          WRITE(LUMOD3,'(A)')'       TMP = BCC(iPrimB,iContB)'
          call initString(8)
          call AddToString('do iTUV=1,')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('*nContC*nContD*nContA')
          call writeString(LUMOD3)               
          WRITE(LUMOD3,'(A)')'        AUXarrayCont(iTUV,iContB,iP)=AUXarrayCont(iTUV,iContB,iP)+TMP*AUXarray2(iTUV,iPrimB,iP)'
          WRITE(LUMOD3,'(A)')'       enddo'
          WRITE(LUMOD3,'(A)')'      enddo'
          WRITE(LUMOD3,'(A)')'     enddo'
          WRITE(LUMOD3,'(A)')'    enddo'               
          IF(CPU)THEN
             WRITE(LUMOD3,'(A)')'!$OMP END DO'
          ENDIF
          call initString(1)
          call AddToString('  end subroutine PrimitiveContractionB'//ARCSTRING//'Gen')
          call AddToString(nTUVP*nTUVQ)
          call writeString(LUMOD3)

! PSEGMENTED
          WRITE(LUMOD5,'(A)')''
          call initString(1)
          call AddToString('  subroutine PrimitiveContractionC'//ARCSTRING//'SegP')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
          call writeString(LUMOD5)
          IF(CPU)THEN
             WRITE(LUMOD5,'(A)')'       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD)'
          ELSE
             WRITE(LUMOD5,'(A)')'       & nContQ,CCC,nPrimC,nContC,nPrimD,nContD,iASync)'
          ENDIF
          WRITE(LUMOD5,'(A)')'    implicit none'
          WRITE(LUMOD5,'(A)')'    !Warning Primitive screening modifies this!!! '
          WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ'
          WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD'
          WRITE(LUMOD5,'(A)')'    real(realk),intent(in) :: CCC(nPrimC,nContC)'
          WRITE(LUMOD5,'(A,I5,A)')'    real(realk),intent(in) :: AUXarray2(',nTUVP*nTUVQ,',nPrimC,nPrimD*nPasses)'
          WRITE(LUMOD5,'(A,I5,A)')'    real(realk),intent(inout) :: AUXarrayCont(',nTUVP*nTUVQ,',nContC,nPrimD*nPasses)'
          IF(.NOT.CPU)WRITE(LUMOD5,'(A)')'    integer(kind=acckind),intent(in) :: iASync'         
          WRITE(LUMOD5,'(A)')'    !'
          WRITE(LUMOD5,'(A)')'    integer :: iP,iContC,iPrimC,iTUV'
          WRITE(LUMOD5,'(A,I5,A)')'    real(realk) :: TMP'
          IF(CPU)THEN
             WRITE(LUMOD5,'(A)')'!$OMP DO &'
             WRITE(LUMOD5,'(A)')'!$OMP PRIVATE(iP,iContC,iPrimC,iTUV,TMP) '
!             WRITE(LUMOD5,'(A)')'!$OMP SHARED(nContC,nContD,nPasses,nPrimC,nPrimD,&'
!             WRITE(LUMOD5,'(A)')'!$OMP        CCC,AUXarrayCont,AUXarray2)'
          ELSE
             WRITE(LUMOD5,'(A)')'!$ACC PARALLEL LOOP &'
             WRITE(LUMOD5,'(A)')'!$ACC PRIVATE(iP,iContC,iPrimC,iTUV,TMP) &'
             WRITE(LUMOD5,'(A)')'!$ACC PRESENT(nContC,nContD,nPasses,nPrimC,nPrimD,&'
             WRITE(LUMOD5,'(A)')'!$ACC        CCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
          ENDIF
          WRITE(LUMOD5,'(A)')'    do iP = 1,nPasses*nPrimD'
          WRITE(LUMOD5,'(A)')'     do iContC=1,nContC'
       WRITE(LUMOD5,'(A,I5)')'      do iTUV=1,',nTUVP*nTUVQ
          WRITE(LUMOD5,'(A)')'       AUXarrayCont(iTUV,iContC,iP) = 0.0E0_realk'
          WRITE(LUMOD5,'(A)')'      enddo'
          WRITE(LUMOD5,'(A)')'      do iPrimC=1,nPrimC'
          WRITE(LUMOD5,'(A)')'       tmp = CCC(iPrimC,iContC)'
       WRITE(LUMOD5,'(A,I5)')'       do iTUV=1,',nTUVP*nTUVQ
          WRITE(LUMOD5,'(A)')'        AUXarrayCont(iTUV,iContC,iP) = AUXarrayCont(iTUV,iContC,iP) + TMP*AUXarray2(iTUV,iPrimC,iP)'
          WRITE(LUMOD5,'(A)')'       enddo'
          WRITE(LUMOD5,'(A)')'      enddo'
          WRITE(LUMOD5,'(A)')'     enddo'
          WRITE(LUMOD5,'(A)')'    enddo'
          IF(CPU)THEN
             WRITE(LUMOD5,'(A)')'!$OMP END DO'
          ENDIF
          call initString(1)
          call AddToString('  end subroutine PrimitiveContractionC'//ARCSTRING//'SegP')
          call AddToString(nTUVP*nTUVQ)
          call writeString(LUMOD5)

          WRITE(LUMOD5,'(A)')''
          call initString(1)
          call AddToString('  subroutine PrimitiveContractionD'//ARCSTRING//'SegP')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
          call writeString(LUMOD5)
          IF(CPU)THEN
             WRITE(LUMOD5,'(A)')'       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD)'
          ELSE
             WRITE(LUMOD5,'(A)')'       & nContQ,DCC,nPrimC,nContC,nPrimD,nContD,iASync)'
          ENDIF
          WRITE(LUMOD5,'(A)')'    implicit none'
          WRITE(LUMOD5,'(A)')'    !Warning Primitive screening modifies this!!! '
          WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContQ'
          WRITE(LUMOD5,'(A)')'    integer,intent(in) :: nPrimC,nContC,nPrimD,nContD'
          WRITE(LUMOD5,'(A)')'    real(realk),intent(in) :: DCC(nPrimD,nContD)'
          WRITE(LUMOD5,'(A,I5,A)')'    real(realk),intent(in) :: AUXarray2(',nTUVP*nTUVQ,'*nContC,nPrimD,nPasses)'
          WRITE(LUMOD5,'(A,I5,A)')'    real(realk),intent(inout) :: AUXarrayCont(',nTUVP*nTUVQ,'*nContC,nContD,nPasses)'
          IF(.NOT.CPU)WRITE(LUMOD5,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
          WRITE(LUMOD5,'(A)')'    !'
          WRITE(LUMOD5,'(A)')'    integer :: iPassP,iContD,iPrimD,iTUV'
          WRITE(LUMOD5,'(A)')'    real(realk) :: TMP'
          IF(CPU)THEN
             WRITE(LUMOD5,'(A)')'!$OMP DO &'
             WRITE(LUMOD5,'(A)')'!$OMP PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) '
!             WRITE(LUMOD5,'(A)')'!$OMP SHARED(nContC,nContD,nPasses,nPrimC,nPrimD,&'
!             WRITE(LUMOD5,'(A)')'!$OMP        DCC,AUXarrayCont,AUXarray2)'
          ELSE
             WRITE(LUMOD5,'(A)')'!$ACC PARALLEL LOOP &'
             WRITE(LUMOD5,'(A)')'!$ACC PRIVATE(iPassP,iContD,iPrimD,iTUV,TMP) &'
             WRITE(LUMOD5,'(A)')'!$ACC PRESENT(nContC,nContD,nPasses,nPrimD,&'
             WRITE(LUMOD5,'(A)')'!$ACC        DCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
          ENDIF
          WRITE(LUMOD5,'(A)')'    do iPassP = 1,nPasses'
          WRITE(LUMOD5,'(A)')'     do iContD=1,nContD'
          WRITE(LUMOD5,'(A,I5,A)')'       do iTUV=1,',nTUVP*nTUVQ,'*nContC'
          WRITE(LUMOD5,'(A)')'        AUXarrayCont(iTUV,iContD,iPassP) = 0.0E0_realk'
          WRITE(LUMOD5,'(A)')'       enddo'
          WRITE(LUMOD5,'(A)')'      do iPrimD=1,nPrimD'
          WRITE(LUMOD5,'(A)')'       tmp = DCC(iPrimD,iContD)'
          WRITE(LUMOD5,'(A,I5,A)')'       do iTUV=1,',nTUVP*nTUVQ,'*nContC'
          WRITE(LUMOD5,'(A)')'        AUXarrayCont(iTUV,iContD,iPassP) = AUXarrayCont(iTUV,iContD,iPassP) + tmp*AUXarray2(iTUV,iPrimD,iPassP)'
          WRITE(LUMOD5,'(A)')'       enddo'
          WRITE(LUMOD5,'(A)')'      enddo'
          WRITE(LUMOD5,'(A)')'     enddo'
          WRITE(LUMOD5,'(A)')'    enddo'
          IF(CPU)THEN
             WRITE(LUMOD5,'(A)')'!$OMP END DO'
          ENDIF
          call initString(1)
          call AddToString('  end subroutine PrimitiveContractionD'//ARCSTRING//'SegP')
          call AddToString(nTUVP*nTUVQ)
          call writeString(LUMOD5)

         ! QSEGMENTED
          WRITE(LUMOD4,'(A)')''
          call initString(1)
          call AddToString('  subroutine PrimitiveContractionA'//ARCSTRING//'SegQ')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
          call writeString(LUMOD4)
          IF(CPU)THEN
             WRITE(LUMOD4,'(A)')'       & nContP,ACC,nPrimA,nContA,nPrimB,nContB)'
          ELSE
             WRITE(LUMOD4,'(A)')'       & nContP,ACC,nPrimA,nContA,nPrimB,nContB,iASync)'
          ENDIF
          WRITE(LUMOD4,'(A)')'    implicit none'
          WRITE(LUMOD4,'(A)')'    !Warning Primitive screening modifies this!!! '
          WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP'
          WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB'
          WRITE(LUMOD4,'(A)')'    real(realk),intent(in) :: ACC(nPrimA,nContA)'
          WRITE(LUMOD4,'(A,I5,A)')'    real(realk),intent(in) :: AUXarray2(',nTUVP*nTUVQ,',nPrimA,nPrimB*nPasses)'
          WRITE(LUMOD4,'(A,I5,A)')'    real(realk),intent(inout) :: AUXarrayCont(',nTUVP*nTUVQ,',nContA,nPrimB*nPasses)'
          IF(.NOT.CPU)WRITE(LUMOD4,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
          WRITE(LUMOD4,'(A)')'    !'
          WRITE(LUMOD4,'(A)')'    integer :: iP,iContA,iPrimA,iTUV'
          WRITE(LUMOD4,'(A)')'    real(realk) :: TMP'
          IF(CPU)THEN
             WRITE(LUMOD4,'(A)')'!$OMP DO &'
             WRITE(LUMOD4,'(A)')'!$OMP PRIVATE(iTUV,iP,iContA,iPrimA,TMP) '
!             WRITE(LUMOD4,'(A)')'!$OMP SHARED(nContA,nPasses,nPrimA,nPrimB,&'
!             WRITE(LUMOD4,'(A)')'!$OMP        ACC,AUXarrayCont,AUXarray2)'
          ELSE
             WRITE(LUMOD4,'(A)')'!$ACC PARALLEL LOOP &'
             WRITE(LUMOD4,'(A)')'!$ACC PRIVATE(iTUV,iP,iContA,iPrimA,TMP) &'
             WRITE(LUMOD4,'(A)')'!$ACC PRESENT(nContA,nPasses,nPrimA,nPrimB,&'
             WRITE(LUMOD4,'(A)')'!$ACC         ACC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
          ENDIF
          WRITE(LUMOD4,'(A)')'   do iP = 1,nPasses*nPrimB'
          WRITE(LUMOD4,'(A)')'    do iContA=1,nContA'
       WRITE(LUMOD4,'(A,I5)')'     do iTUV=1,',nTUVP*nTUVQ
          WRITE(LUMOD4,'(A)')'      AUXarrayCont(iTUV,iContA,iP) = 0.0E0_realk'
          WRITE(LUMOD4,'(A)')'     enddo'
          WRITE(LUMOD4,'(A)')'     do iPrimA=1,nPrimA'
          WRITE(LUMOD4,'(A)')'      TMP = ACC(iPrimA,iContA)'
       WRITE(LUMOD4,'(A,I5)')'      do iTUV=1,',nTUVP*nTUVQ
          WRITE(LUMOD4,'(A)')'       AUXarrayCont(iTUV,iContA,iP) = AUXarrayCont(iTUV,iContA,iP) + tmp*AUXarray2(iTUV,iPrimA,iP)'
          WRITE(LUMOD4,'(A)')'      enddo'
          WRITE(LUMOD4,'(A)')'     enddo'
          WRITE(LUMOD4,'(A)')'    enddo'
          WRITE(LUMOD4,'(A)')'   enddo'
          IF(CPU)THEN
             WRITE(LUMOD4,'(A)')'!$OMP END DO'
          ENDIF
          call initString(1)
          call AddToString('  end subroutine PrimitiveContractionA'//ARCSTRING//'SegQ')
          call AddToString(nTUVP*nTUVQ)
          call writeString(LUMOD4)
          WRITE(LUMOD4,'(A)')''

          WRITE(LUMOD4,'(A)')''
          call initString(1)
          call AddToString('  subroutine PrimitiveContractionB'//ARCSTRING//'SegQ')
          call AddToString(nTUVP*nTUVQ)
          call AddToString('(AUXarray2,AUXarrayCont,nPrimP,nPrimQ,nPasses,&')
          call writeString(LUMOD4)
          IF(CPU)THEN
             WRITE(LUMOD4,'(A)')'       & nContP,BCC,nPrimA,nContA,nPrimB,nContB)'
          ELSE
             WRITE(LUMOD4,'(A)')'       & nContP,BCC,nPrimA,nContA,nPrimB,nContB,iASync)'
          ENDIF
          WRITE(LUMOD4,'(A)')'    implicit none'
          WRITE(LUMOD4,'(A)')'    !Warning Primitive screening modifies this!!! '
          WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimP,nPrimQ,nPasses,nContP'
          WRITE(LUMOD4,'(A)')'    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB'
          WRITE(LUMOD4,'(A)')'    real(realk),intent(in) :: BCC(nPrimB,nContB)'
          WRITE(LUMOD4,'(A,I5,A)')'    real(realk),intent(in) :: AUXarray2(',nTUVP*nTUVQ,'*nContA,nPrimB,nPasses)'
          WRITE(LUMOD4,'(A,I5,A)')'    real(realk),intent(inout) :: AUXarrayCont(',nTUVP*nTUVQ,'*nContA,nContB,nPasses)'
          IF(.NOT.CPU)WRITE(LUMOD4,'(A)')'    integer(kind=acckind),intent(in) :: iASync'
          WRITE(LUMOD4,'(A)')'    !'
          WRITE(LUMOD4,'(A)')'    integer :: iPassP,iContB,iPrimB,iTUV'
          WRITE(LUMOD4,'(A)')'    real(realk) :: TMP'
          IF(CPU)THEN
             WRITE(LUMOD4,'(A)')'!$OMP DO &'
             WRITE(LUMOD4,'(A)')'!$OMP PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) '
!             WRITE(LUMOD4,'(A)')'!$OMP SHARED(nContA,nContB,nPasses,nPrimB,&'
!             WRITE(LUMOD4,'(A)')'!$OMP        BCC,AUXarrayCont,AUXarray2)'
          ELSE
             WRITE(LUMOD4,'(A)')'!$ACC PARALLEL LOOP &'
             WRITE(LUMOD4,'(A)')'!$ACC PRIVATE(iTUV,iPassP,iContB,iPrimB,TMP) &'
             WRITE(LUMOD4,'(A)')'!$ACC PRESENT(nContA,nContB,nPasses,nPrimB,&'
             WRITE(LUMOD4,'(A)')'!$ACC         BCC,AUXarrayCont,AUXarray2) ASYNC(iASync)'
          ENDIF
          WRITE(LUMOD4,'(A)')'   do iPassP = 1,nPasses'
          WRITE(LUMOD4,'(A)')'    do iContB=1,nContB'
     WRITE(LUMOD4,'(A,I5,A)')'     do iTUV=1,',nTUVP*nTUVQ,'*nContA'
          WRITE(LUMOD4,'(A)')'      AUXarrayCont(iTUV,iContB,iPassP)=0.0E0_realk'
          WRITE(LUMOD4,'(A)')'     enddo'
          WRITE(LUMOD4,'(A)')'     do iPrimB=1,nPrimB'
          WRITE(LUMOD4,'(A)')'      TMP = BCC(iPrimB,iContB)'
     WRITE(LUMOD4,'(A,I5,A)')'      do iTUV=1,',nTUVP*nTUVQ,'*nContA'
          WRITE(LUMOD4,'(A)')'       AUXarrayCont(iTUV,iContB,iPassP)=AUXarrayCont(iTUV,iContB,iPassP)+tmp*AUXarray2(iTUV,iPrimB,iPassP)'
          WRITE(LUMOD4,'(A)')'      enddo'
          WRITE(LUMOD4,'(A)')'     enddo'
          WRITE(LUMOD4,'(A)')'    enddo'
          WRITE(LUMOD4,'(A)')'   enddo'
          IF(CPU)THEN
             WRITE(LUMOD4,'(A)')'!$OMP END DO'
          ENDIF
          call initString(1)
          call AddToString('  end subroutine PrimitiveContractionB'//ARCSTRING//'SegQ')
          call AddToString(nTUVP*nTUVQ)
          call writeString(LUMOD4)
          WRITE(LUMOD4,'(A)')''
       ENDIF
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   deallocate(UniquenTUVs)
ENDIF

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
!!$WRITE(LUMOD3,'(A)')'    integer :: iPrimQ,iPassP,iPrimP'
!!$WRITE(LUMOD3,'(A)')'    real(realk) :: px,py,pz,pqx,pqy,pqz'
!!$WRITE(LUMOD3,'(A)')'    DO iPassP=1, nPasses'
!!$WRITE(LUMOD3,'(A)')'       DO iPrimP=1, nPrimP'
!!$WRITE(LUMOD3,'(A)')'          px = Pcent(1,iPrimP)'
!!$WRITE(LUMOD3,'(A)')'          py = Pcent(2,iPrimP)'
!!$WRITE(LUMOD3,'(A)')'          pz = Pcent(3,iPrimP)'
!!$WRITE(LUMOD3,'(A)')'          DO iPrimQ=1, nPrimQ'
!!$WRITE(LUMOD3,'(A)')'             pqx = px - Qcent(1,iPrimQ,iPassP)'
!!$WRITE(LUMOD3,'(A)')'             pqy = py - Qcent(2,iPrimQ,iPassP)'
!!$WRITE(LUMOD3,'(A)')'             pqz = pz - Qcent(3,iPrimQ,iPassP)'
!!$WRITE(LUMOD3,'(A)')'             rPQ(1,iPrimQ,iPrimP,iPassP) = pqx'
!!$WRITE(LUMOD3,'(A)')'             rPQ(2,iPrimQ,iPrimP,iPassP) = pqy'
!!$WRITE(LUMOD3,'(A)')'             rPQ(3,iPrimQ,iPrimP,iPassP) = pqz'
!!$WRITE(LUMOD3,'(A)')'             squaredDistance(iPrimQ,iPrimP,iPassP) = pqx*pqx+pqy*pqy+pqz*pqz'
!!$WRITE(LUMOD3,'(A)')'          ENDDO'
!!$WRITE(LUMOD3,'(A)')'       ENDDO'
!!$WRITE(LUMOD3,'(A)')'    ENDDO'
!!$WRITE(LUMOD3,'(A)')'  end subroutine build_QP_squaredDistance_and_Rpq'
  !WRITE(LUMOD3,'(A)')''
!!$WRITE(LUMOD3,'(A)')'  SUBROUTINE buildRJ000_general(nPasses,nPrimQ,nPrimP,nTABFJW1,nTABFJW2,reducedExponents,&'
!!$WRITE(LUMOD3,'(A)')'       & TABFJW,RJ000,JMAX,Pcent,Qcent)'
!!$WRITE(LUMOD3,'(A)')'    IMPLICIT NONE'
!!$WRITE(LUMOD3,'(A)')'    INTEGER,intent(in)         :: nPrimP,nPrimQ,Jmax,nTABFJW1,nTABFJW2,nPasses'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK),intent(in)     :: reducedExponents(nPrimQ,nPrimP)'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK),intent(in)     :: Pcent(3,nPrimP),Qcent(3,nPrimQ,nPasses)'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK),intent(in)     :: TABFJW(0:nTABFJW1,0:nTABFJW2)'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrimQ*nPrimP,nPasses)'
!!$WRITE(LUMOD3,'(A)')'    !'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK)     :: D2JP36,WVAL'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK),PARAMETER :: HALF =0.5E0_realk,D1=1E0_realk,D2 = 2E0_realk, D4 = 4E0_realk, D100=100E0_realk'
!!$WRITE(LUMOD3,'(A)')'    Real(realk),parameter :: D12 = 12E0_realk, TENTH = 0.01E0_realk'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk'
!!$WRITE(LUMOD3,'(A)')'    Integer :: IPNT,J'
!!$WRITE(LUMOD3,'(A)')'    Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk'
!!$WRITE(LUMOD3,'(A)')'    Real(realk), parameter :: PI=3.14159265358979323846E0_realk'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI'
!!$WRITE(LUMOD3,'(A)')'    Real(realk) :: W2,W3,R'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK), PARAMETER :: SMALL = 1E-15_realk'
!!$WRITE(LUMOD3,'(A)')'    REAL(REALK) :: PX,PY,PZ,PQX,PQY,PQZ,squaredDistance'
!!$WRITE(LUMOD3,'(A)')'    integer :: iPassP,iPQ,iPrimP,iPrimQ'
!!$WRITE(LUMOD3,'(A)')'    !make for different values of JMAX => loop unroll  '
!!$WRITE(LUMOD3,'(A)')'    !sorting? '
!!$WRITE(LUMOD3,'(A)')'    D2JP36 = 2*JMAX + 36'
!!$WRITE(LUMOD3,'(A)')'    DO iPassP=1, nPasses'
!!$WRITE(LUMOD3,'(A)')'      ipq = 0'
!!$WRITE(LUMOD3,'(A)')'      DO iPrimP=1, nPrimP'
!!$WRITE(LUMOD3,'(A)')'        px = Pcent(1,iPrimP)'
!!$WRITE(LUMOD3,'(A)')'        py = Pcent(2,iPrimP)'
!!$WRITE(LUMOD3,'(A)')'        pz = Pcent(3,iPrimP)'
!!$WRITE(LUMOD3,'(A)')'        DO iPrimQ=1, nPrimQ'
!!$WRITE(LUMOD3,'(A)')'          ipq = ipq + 1'
!!$WRITE(LUMOD3,'(A)')'          pqx = px - Qcent(1,iPrimQ,iPassP)'
!!$WRITE(LUMOD3,'(A)')'          pqy = py - Qcent(2,iPrimQ,iPassP)'
!!$WRITE(LUMOD3,'(A)')'          pqz = pz - Qcent(3,iPrimQ,iPassP)'
!!$WRITE(LUMOD3,'(A)')'          squaredDistance = pqx*pqx+pqy*pqy+pqz*pqz'
!!$WRITE(LUMOD3,'(A)')'          WVAL = reducedExponents(iPrimQ,iPrimP)*squaredDistance'
!!$WRITE(LUMOD3,'(A)')'          !  0 < WVAL < 0.000001'
!!$WRITE(LUMOD3,'(A)')'          IF (ABS(WVAL) .LT. SMALL) THEN'         
!!$WRITE(LUMOD3,'(A)')'             RJ000(0,ipq,ipassq) = D1'
!!$WRITE(LUMOD3,'(A)')'             DO J=1,JMAX'
!!$WRITE(LUMOD3,'(A)')'                RJ000(J,ipq,ipassq)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT'
!!$WRITE(LUMOD3,'(A)')'             ENDDO'
!!$WRITE(LUMOD3,'(A)')'             !  0 < WVAL < 12 '
!!$WRITE(LUMOD3,'(A)')'          ELSE IF (WVAL .LT. D12) THEN'
!!$WRITE(LUMOD3,'(A)')'             IPNT = NINT(D100*WVAL)'
!!$WRITE(LUMOD3,'(A)')'             WDIFF = WVAL - TENTH*IPNT'
!!$WRITE(LUMOD3,'(A)')'             W2    = WDIFF*WDIFF'
!!$WRITE(LUMOD3,'(A)')'             W3    = W2*WDIFF'
!!$WRITE(LUMOD3,'(A)')'             W2    = W2*COEF2'
!!$WRITE(LUMOD3,'(A)')'             W3    = W3*COEF3'
!!$WRITE(LUMOD3,'(A)')'             DO J=0,JMAX'
!!$WRITE(LUMOD3,'(A)')'                R = TABFJW(J,IPNT)'
!!$WRITE(LUMOD3,'(A)')'                R = R -TABFJW(J+1,IPNT)*WDIFF'
!!$WRITE(LUMOD3,'(A)')'                R = R + TABFJW(J+2,IPNT)*W2'
!!$WRITE(LUMOD3,'(A)')'                R = R + TABFJW(J+3,IPNT)*W3'
!!$WRITE(LUMOD3,'(A)')'                RJ000(J,ipq,ipassq) = R'
!!$WRITE(LUMOD3,'(A)')'             ENDDO'
!!$WRITE(LUMOD3,'(A)')'             !  12 < WVAL <= (2J+36) '
!!$WRITE(LUMOD3,'(A)')'          ELSE IF (WVAL.LE.D2JP36) THEN'
!!$WRITE(LUMOD3,'(A)')'             REXPW = HALF*EXP(-WVAL)'
!!$WRITE(LUMOD3,'(A)')'             RWVAL = D1/WVAL'
!!$WRITE(LUMOD3,'(A)')'             GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))'
!!$WRITE(LUMOD3,'(A)')'             RJ000(0,ipq,ipassq) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL'
!!$WRITE(LUMOD3,'(A)')'             DO J=1,JMAX'
!!$WRITE(LUMOD3,'(A)')'                RJ000(J,ipq,ipassq) = RWVAL*((J - HALF)*RJ000(J-1,ipq,ipassq)-REXPW)'
!!$WRITE(LUMOD3,'(A)')'             ENDDO'
!!$WRITE(LUMOD3,'(A)')'             !  (2J+36) < WVAL '
!!$WRITE(LUMOD3,'(A)')'          ELSE'
!!$WRITE(LUMOD3,'(A)')'             RWVAL = PID4/WVAL'
!!$WRITE(LUMOD3,'(A)')'             RJ000(0,ipq,ipassq) = SQRT(RWVAL)'
!!$WRITE(LUMOD3,'(A)')'             RWVAL = RWVAL*PID4I'
!!$WRITE(LUMOD3,'(A)')'             DO J = 1, JMAX'
!!$WRITE(LUMOD3,'(A)')'                RJ000(J,ipq,ipassq) = RWVAL*(J - HALF)*RJ000(J-1,ipq,ipassq)'
!!$WRITE(LUMOD3,'(A)')'             ENDDO'
!!$WRITE(LUMOD3,'(A)')'          ENDIF'
!!$WRITE(LUMOD3,'(A)')'        ENDDO'
!!$WRITE(LUMOD3,'(A)')'      ENDDO'
!!$WRITE(LUMOD3,'(A)')'    ENDDO'
!!$WRITE(LUMOD3,'(A)')'  END SUBROUTINE buildRJ000_general'
!  WRITE(LUMOD3,'(A)')'  '
!  WRITE(LUMOD3,'(A)')'END MODULE IchorEriCoulombintegralOBSGeneralMod'

  WRITE(LUMOD2,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralMod'
  WRITE(LUMOD3,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGen'
  WRITE(LUMOD4,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQ'
  WRITE(LUMOD5,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegP'
  WRITE(LUMOD6,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg'
  WRITE(LUMOD7,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1Prim'

  close(unit = LUMOD2)
  close(unit = LUMOD3)
  close(unit = LUMOD4)
  close(unit = LUMOD5)
  close(unit = LUMOD6)
  close(unit = LUMOD7)

enddo

contains
  subroutine CallPrimitiveContractionString1(LUMOD3,ARCSTRING2,BASISSPEC,iBasisSpec,&
       & nTUVQ,nTUVP,STRINGIN,STRINGOUT,CenterString)
    implicit none
    character(len=3) :: ARCSTRING2
    character(len=1) :: CenterString
    character(len=8) :: BASISSPEC
    character(len=9) :: STRINGIN,STRINGOUT
    integer :: iBasisSpec,nTUVP,nTUVQ,LUMOD3
    call initString(9)
    call AddToString('call PrimitiveContraction')
    call AddToString(CenterString)
    call AddToString(ARCSTRING2)
    call AddToString(BASISSPEC(1:iBasisSpec))
    IF(nPrimLast)THEN
       call AddToString(nTUVQ*nTUVP)
    ENDIF
    call AddToString('(')
    call AddToString(STRINGIN)
    call AddToString(',')
    call AddToString(STRINGOUT)
    call AddToString(',nPrimP,nPrimQ,nPasses,&')
    call writeString(LUMOD3)             
    call initString(15)
    IF(CenterString.EQ.'C'.OR.CenterString.EQ.'D')THEN
       call AddToString('& nContQ,')
       call AddToString(CenterString)
       IF(nPrimLast)THEN
          call AddToString('CC,nPrimC,nContC,nPrimD,nContD')
       ELSE
          call AddToString('CC,nPrimC,nContC,nPrimD,nContD,')
          call AddToString(nTUVP)
          call AddToString(',')
          call AddToString(nTUVQ)
       ENDIF
    ELSE !A or B
       call AddToString('& nContP,')
       call AddToString(CenterString)
       IF(nPrimLast)THEN
          call AddToString('CC,nPrimA,nContA,nPrimB,nContB')
       ELSE
          call AddToString('CC,nPrimA,nContA,nPrimB,nContB,')
          call AddToString(nTUVP)
          call AddToString(',')
          call AddToString(nTUVQ)
       ENDIF
    ENDIF
    IF(CPU)THEN
       call AddToString(')')
    ELSE
       call AddToString(',iASync)')
    ENDIF
    
    call writeString(LUMOD3)             
  end subroutine CallPrimitiveContractionString1

  subroutine CallPrimitiveContractionString2(LUMOD3,ARCSTRING2,BASISSPEC,iBasisSpec,&
       & nTUVQ,nTUVP,STRINGIN,STRINGOUT,CenterString)
    implicit none
    character(len=3) :: ARCSTRING2
    character(len=1) :: CenterString
    character(len=8) :: BASISSPEC
    character(len=9) :: STRINGIN,STRINGOUT
    integer :: iBasisSpec,nTUVP,nTUVQ,LUMOD3
    call initString(9)
    call AddToString('call PrimitiveContraction')
    call AddToString(CenterString)
    call AddToString(ARCSTRING2)
    call AddToString(BASISSPEC(1:iBasisSpec))
    IF(nPrimLast)THEN
       call AddToString(nTUVQ*nTUVP)
    ENDIF
    call AddToString('(')
    call AddToString(STRINGIN)
    call AddToString(',')
    call AddToString(STRINGOUT)
    call AddToString(',nPrimP,nPrimQ,nPasses,&')
    call writeString(LUMOD3)             
    call initString(15)
    call AddToString('& nContP,nContQ,')
    call AddToString(CenterString)
    call AddToString('CC,nPrimA,nContA,nPrimB,nContB,nPrimC,&')
    call writeString(LUMOD3)             
    IF(nPrimLast)THEN
       IF(CPU)THEN
          WRITE(LUMOD3,'(A)')'              & nContC,nPrimD,nContD)'
       ELSE
          WRITE(LUMOD3,'(A)')'              & nContC,nPrimD,nContD,iASync)'
       ENDIF
    ELSE
       call initString(15)
       call AddToString('& nContC,nPrimD,nContD,')
       call AddToString(nTUVP)
       call AddToString(',')
       call AddToString(nTUVQ)
       IF(CPU)THEN
          call AddToString(')')
       ELSE
          call AddToString(',iASync)')
       ENDIF
       call writeString(LUMOD3)             
    ENDIF
  end subroutine CallPrimitiveContractionString2

  subroutine subroutineMain(LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
       & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,nlmA,nlmB,nlmC,nlmD,&
       & spherical,Gen,SegQ,SegP,Seg,Seg1Prim)
    implicit none
    integer,intent(in) :: LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,nTUV,AngmomP,AngmomQ
    integer,intent(in) :: nlmA,nlmB,nlmC,nlmD
    integer,intent(in) :: AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec
    character(len=9) :: STRINGIN,STRINGOUT,TMPSTRING
    logical :: spherical,OutputSet,Gen,SegQ,SegP,Seg,Seg1Prim,Contracted,noTransfer
    character(len=8) :: BASISSPEC
    integer :: iBasisSpec
    Character(len=1)  :: FromLabel,ToLabel,FromExpLabel,ToExpLabel
    Character(len=1)  :: CenterString,QPCenterString
    IF(Gen)THEN
       iBasisSpec = 3
       BASISSPEC = 'Gen     '
    ELSEIF(SegQ)THEN
       iBasisSpec = 4
       BASISSPEC = 'SegQ    '
    ELSEIF(SegP)THEN
       iBasisSpec = 4
       BASISSPEC = 'SegP    '
    ELSEIF(Seg)THEN
       iBasisSpec = 3
       BASISSPEC = 'Seg     '
    ELSEIF(Seg1Prim)THEN
       iBasisSpec = 8
       BASISSPEC = 'Seg1Prim'
    ENDIF
    OutputSet = .FALSE.
    Contracted = .FALSE.
    IF(AngmomPQ.EQ.0)THEN
       IF(Gen)THEN
          call DebugMemoryTest(STRINGOUT,'nPrimP*nPrimQ*nPasses',nTUV,LUMOD3)
       ELSEIF(Seg)THEN
          call DebugMemoryTest(STRINGOUT,'nPasses',nTUV,LUMOD3)
       ELSEIF(SegQ)THEN
          call DebugMemoryTest(STRINGOUT,'nPrimP*nPasses',nTUV,LUMOD3)
       ELSEIF(SegP)THEN
          call DebugMemoryTest(STRINGOUT,'nPrimQ*nPasses',nTUV,LUMOD3)
       ELSEIF(Seg1Prim)THEN
          call DebugMemoryTest(STRINGOUT,'nPasses',nTUV,LUMOD3)
       ENDIF
       WRITE(LUMOD3,'(A)')'        call VerticalRecurrence'//ARCSTRING//BASISSPEC(1:iBasisSpec)//'0(nPasses,nPrimP,nPrimQ,&'
       WRITE(LUMOD3,'(A)')'               & reducedExponents,TABFJW,Pcent,Qcent,integralPrefactor,&'
       WRITE(LUMOD3,'(A)')'               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,&'
       call initString(15)
       call AddToString('& PpreExpFac,QpreExpFac,')
       IF(Gen)THEN
          call AddToString(STRINGOUT)
          !          call AddToString('(1:nPrimP*nPrimQ*nPasses*')
          !       call AddToString(nTUV)
          !       call AddToString(')')                
       ELSEIF(Seg)THEN
          call AddToString('LOCALINTS')
          OutputSet = .TRUE.
          !          call AddToString('(1:nPasses*')
          !       call AddToString(nTUV)
          !       call AddToString(')')                
          Contracted = .TRUE.
       ELSEIF(SegQ)THEN
          call AddToString(STRINGOUT)
          !          call AddToString('(1:nPrimP*nPasses*')
          !       call AddToString(nTUV)
          !       call AddToString(')')                
!          Contracted = .TRUE.
       ELSEIF(SegP)THEN
          call AddToString(STRINGOUT)
!          Contracted = .TRUE.
       ELSEIF(Seg1Prim)THEN
          call AddToString('LOCALINTS')
          OutputSet = .TRUE.
          Contracted = .TRUE.
       ENDIF
       IF(CPU)THEN
          call AddToString(')')                
       ELSE
          call AddToString(',iASync)')
       ENDIF
       call writeString(LUMOD3)
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
       WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
    ELSE
       !build RJ000 
       IF(AngmomPQ.GT.1)THEN
          call DebugMemoryTest(STRINGOUT,'nPrimP*nPrimQ*nPasses',AngmomPQ+1,LUMOD3)
          call initString(8)
          call AddToString('call BuildRJ000')
          call AddToString(ARCSTRING)
          IF(Seg1Prim)THEN
             call AddToString(BASISSPEC(1:iBasisSpec))
          ELSE
             call AddToString('Gen')
          ENDIF
          call AddToString(AngmomPQ)
!          call AddToString(centerstring)                
          call AddToString('(nPasses,nPrimP,nPrimQ,reducedExponents,&')
          call writeString(LUMOD3)

          call initString(15)
          call AddToString('& TABFJW,Pcent,Qcent,IatomApass,IatomBpass,&')
          call writeString(LUMOD3)

          call initString(15)
          call AddToString('& MaxPasses,nAtomsA,nAtomsB,')
          call AddToString(STRINGOUT)
          IF(CPU)THEN
             call AddToString(')')
          ELSE
             call AddToString(',iASync)')
          ENDIF
          call writeString(LUMOD3)
          
          TMPSTRING = STRINGIN
          STRINGIN  = STRINGOUT
          STRINGOUT  = TMPSTRING
       ENDIF

       IF(AngmomP.GE.AngmomQ)THEN
          noTransfer = AngmomQ.EQ.0
          IF(AngmomA.GE.AngmomB)THEN
             !A Vertical recurrence
             centerString='A'
             QPcenterString='P'
             IF(AngmomC.GE.AngmomD)THEN
                !A to C TransferRecurrence
                FromLabel = 'A'; ToLabel = 'C'; FromExpLabel = 'B'; ToExpLabel = 'D'
                SPEC = 'AtoC'
             ELSE
                !A to D TransferRecurrence
                FromLabel = 'A'; ToLabel = 'D'; FromExpLabel = 'B'; ToExpLabel = 'C' 
                SPEC = 'AtoD'
             ENDIF
          ELSE
             !B Vertical recurrence
             centerString='B'
             QPcenterString='P'
             IF(AngmomC.GE.AngmomD)THEN
                !A to C TransferRecurrence
FromLabel = 'B'; ToLabel = 'C'; FromExpLabel = 'A'; ToExpLabel = 'D' 
                SPEC = 'BtoC'
             ELSE
                !A to D TransferRecurrence
                FromLabel = 'B'; ToLabel = 'D'; FromExpLabel = 'A'; ToExpLabel = 'C' 
                SPEC = 'BtoD'
             ENDIF
          ENDIF
       ELSE
          noTransfer = AngmomP.EQ.0
          IF(AngmomC.GE.AngmomD)THEN
             !C Vertical recurrence
             centerString='C'
             QPcenterString='Q'
             IF(AngmomA.GE.AngmomB)THEN
                !C to A TransferRecurrence
                FromLabel = 'C'; ToLabel = 'A'; FromExpLabel = 'D'; ToExpLabel = 'B'
                SPEC = 'CtoA'
             ELSE
                !C to B TransferRecurrence
                FromLabel = 'C'; ToLabel = 'B'; FromExpLabel = 'D'; ToExpLabel = 'A' 
                SPEC = 'CtoB'
             ENDIF
          ELSE
             !D Vertical recurrence
             centerString='D'
             QPcenterString='Q'
             IF(AngmomA.GE.AngmomB)THEN
                !C to A TransferRecurrence
                FromLabel = 'D'; ToLabel = 'A'; FromExpLabel = 'C'; ToExpLabel = 'B' 
                SPEC = 'DtoA'
             ELSE
                !C to B TransferRecurrence
                FromLabel = 'D'; ToLabel = 'B'; FromExpLabel = 'C'; ToExpLabel = 'A' 
                SPEC = 'DtoB'
             ENDIF
          ENDIF
       ENDIF

       !VERTICAL
       IF(noTransfer)THEN
          !no Electron Transfer Recurrence Relation so in case of segmented we can add together now
          IF(Gen)THEN
             call DebugMemoryTest(STRINGOUT,'nPrimP*nPrimQ*nPasses',nTUV,LUMOD3)
          ELSEIF(Seg)THEN
             call DebugMemoryTest(STRINGOUT,'nPasses',nTUV,LUMOD3)
             Contracted = .TRUE.
          ELSEIF(SegQ)THEN
             call DebugMemoryTest(STRINGOUT,'nPrimP*nPasses',nTUV,LUMOD3)
             Contracted = .TRUE.
          ELSEIF(SegP)THEN
             call DebugMemoryTest(STRINGOUT,'nPrimQ*nPasses',nTUV,LUMOD3)
             Contracted = .TRUE.
          ELSEIF(Seg1Prim)THEN
             call DebugMemoryTest(STRINGOUT,'nPasses',nTUV,LUMOD3)
             Contracted = .TRUE.
          ENDIF
          call initString(8)
          call AddToString('call VerticalRecurrence')
          call AddToString(ARCSTRING)
          call AddToString(BASISSPEC(1:iBasisSpec))
       ELSE
          call DebugMemoryTest(STRINGOUT,'nPrimP*nPrimQ*nPasses',nTUV,LUMOD3)
          call initString(8)
          call AddToString('call VerticalRecurrence')
          call AddToString(ARCSTRING)
          IF(Seg1Prim)THEN
             call AddToString(BASISSPEC(1:iBasisSpec))
          ELSE
             call AddToString('Gen')
          ENDIF
       ENDIF
       call AddToString(AngmomPQ)
       call AddToString(centerstring)                
       call AddToString('(nPasses,nPrimP,nPrimQ,reducedExponents,&')
       call writeString(LUMOD3)

       call initString(15)
       call AddToString('& ')
       IF(AngmomPQ.EQ.1)THEN
          call AddToString('TABFJW')
       ELSE
          call AddToString(STRINGIN)
       ENDIF
       call AddToString(',')
       call AddToString(QPcenterString)       
       call AddToString('exp,')
       call AddToString(centerString)
       call AddToString('center,Pcent,Qcent,integralPrefactor,&')
       call writeString(LUMOD3)
       WRITE(LUMOD3,'(A)')'               & IatomApass,IatomBpass,MaxPasses,nAtomsA,nAtomsB,PpreExpFac,QpreExpFac,&'
       call initString(15)
       call AddToString('& ')
       call AddToString(STRINGOUT)
!             IF(Gen)THEN
!                call AddToString('(1:nPrimP*nPrimQ*nPasses*')
!             ELSEIF(Seg)THEN
!                call AddToString('(1:nPasses*')
!             ELSEIF(SegQ)THEN
!                call AddToString('(1:nPrimP*nPasses*')
!             ELSEIF(SegP)THEN
!                call AddToString('(1:nPrimQ*nPasses*')
!             ELSEIF(Seg1Prim)THEN
!                call AddToString('(1:nPasses*')
!             ENDIF
!             call AddToString(nTUV)
!             call AddToString(')')                
       IF(CPU)THEN
          call AddToString(')')
       ELSE
          call AddToString(',iASync)')
       ENDIF         
       call writeString(LUMOD3)             

       !WRITE(LUMOD3,'(A,A,A)')'               & ',STRINGOUT,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
       
       !determine TransferRecurrence
       IF(noTransfer)THEN
          WRITE(LUMOD3,'(A)')'        !No reason for the Electron Transfer Recurrence Relation '
          !No reason for the Electron Transfer Recurrence Relation 
       ELSE
          IF(Gen)THEN
             call DebugMemoryTest(STRINGOUT,'nPrimP*nPrimQ*nPasses',nTUVP*nTUVQ,LUMOD3)
          ELSEIF(Seg)THEN
             call DebugMemoryTest(STRINGOUT,'nPasses',nTUVP*nTUVQ,LUMOD3)
          ELSEIF(SegP)THEN
             call DebugMemoryTest(STRINGOUT,'nPrimQ*nPasses',nTUVP*nTUVQ,LUMOD3)
          ELSEIF(SegQ)THEN
             call DebugMemoryTest(STRINGOUT,'nPrimP*nPasses',nTUVP*nTUVQ,LUMOD3)
          ELSEIF(Seg1Prim)THEN
             call DebugMemoryTest(STRINGOUT,'nPasses',nTUVP*nTUVQ,LUMOD3)
          ENDIF
          call initString(8)
          call AddToString('call TransferRecurrence'//ARCSTRING//'P')
          call AddToString(AngmomP)
          call AddToString('Q')
          call AddToString(AngmomQ)
          call AddToString(SPEC)                
          call AddToString(BASISSPEC(1:iBasisSpec))
          call AddToString('(nPasses,nPrimP,nPrimQ,reducedExponents,&')
          call writeString(LUMOD3)
          IF(.NOT.Gen)Contracted = .TRUE.
          call initString(15)
          call AddToString('& Pexp,Qexp,Pdistance12,Qdistance12,')          
          call AddToString(FromExpLabel)
          call AddToString('exp,')
          call AddToString(ToExpLabel)
          call AddToString('exp,nPrimA,nPrimB,nPrimC,nPrimD,&')
          call writeString(LUMOD3)             
          call initString(15)
          call AddToString('& MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,&')
          call writeString(LUMOD3)             
          call initString(15)
          call AddToString('& ')
          call AddToString(STRINGIN)
!          IF(Gen)THEN
!             call AddToString('(1:nPrimP*nPrimQ*nPasses*')
!          ELSEIF(Seg)THEN
!             call AddToString('(1:nPasses*')
!          ELSEIF(SegQ)THEN
!             call AddToString('(1:nPrimP*nPasses*')
!          ELSEIF(SegP)THEN
!             call AddToString('(1:nPrimQ*nPasses*')
!          ELSEIF(Seg1Prim)THEN
!             call AddToString('(1:nPasses*')
!          ENDIF
!          call AddToString(nTUV)
!                call AddToString(')')
          call AddToString(',')
          call AddToString(STRINGOUT)
!          IF(Gen)THEN
!             call AddToString('(1:nPrimP*nPrimQ*nPasses*')
!          ELSEIF(Seg)THEN
!             call AddToString('(1:nPasses*')
!          ELSEIF(SegQ)THEN
!             call AddToString('(1:nPrimP*nPasses*')
!          ELSEIF(SegP)THEN
!             call AddToString('(1:nPrimQ*nPasses*')
!          ELSEIF(Seg1Prim)THEN
!             call AddToString('(1:nPasses*')
!          ENDIF
!          call AddToString(nTUVP*nTUVQ)
!          call AddToString(')')
          IF(CPU)THEN
             call AddToString(')')
          ELSE
             call AddToString(',iASync)')
          ENDIF
          call writeString(LUMOD3)             
          !swap 
          TMPSTRING = STRINGIN
          STRINGIN  = STRINGOUT
          STRINGOUT  = TMPSTRING
       ENDIF
    ENDIF
    !================= DONE WITH VERTICAL AND TRANSFER ================================================================
    IF(Gen)THEN
       WRITE(LUMOD3,'(A)')'        nContQP = nContQ*nContP'
    ENDIF
    IF(nTUVP*nTUVQ.LT.10)THEN
       !       WRITE(LUMOD3,'(A,A,A,I1,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.100)THEN
       !       WRITE(LUMOD3,'(A,A,A,I2,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.1000)THEN
       !       WRITE(LUMOD3,'(A,A,A,I3,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ELSEIF(nTUVP*nTUVQ.LT.10000)THEN
       !       WRITE(LUMOD3,'(A,A,A,I4,A)')'        call mem_ichor_alloc(',STRINGOUT,',',nTUVP*nTUVQ,'*nContQP*nPasses)'
    ENDIF

    IF(Contracted)THEN
       IF(SegP.OR.SegQ)THEN          
          IF(SegP)THEN
             call DebugMemoryTest(STRINGOUT,'nContC*nPrimD*nPasses',nTUVP*nTUVQ,LUMOD3)
             call CallPrimitiveContractionString1(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'C')
             TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING          
             IF(AngmomP.EQ.0.AND.AngmomQ.EQ.0)THEN
                !no subsequent Horizontal or spherical transformations 
                IF(.NOT.OutputSet)THEN
                   STRINGOUT  = 'LOCALINTS     ';OutputSet = .TRUE.
                ELSE
                   STOP 'LOCALINTS already set MAJOER PROBLEM A2'
                ENDIF
             ELSE
                call DebugMemoryTest(STRINGOUT,'nContC*nContD*nPasses',nTUVP*nTUVQ,LUMOD3)
             ENDIF
             call CallPrimitiveContractionString1(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'D')
          ELSE !segQ
             call DebugMemoryTest(STRINGOUT,'nContA*nPrimB*nPasses',nTUVP*nTUVQ,LUMOD3)
             call CallPrimitiveContractionString1(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'A')
             TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING          
             IF(AngmomP.EQ.0.AND.AngmomQ.EQ.0)THEN
                !no subsequent Horizontal or spherical transformations 
                IF(.NOT.OutputSet)THEN
                   STRINGOUT  = 'LOCALINTS     ';OutputSet = .TRUE.
                ELSE
                   STOP 'LOCALINTS already set MAJOER PROBLEM A2'
                ENDIF
             ELSE
                call DebugMemoryTest(STRINGOUT,'nContA*nContB*nPasses',nTUVP*nTUVQ,LUMOD3)
             ENDIF
             call CallPrimitiveContractionString1(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'B')
          ENDIF
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING          
       ELSE
          WRITE(LUMOD3,'(A)')'        !Primitive Contraction have already been done'
       ENDIF
    ELSE
       IF(SegP)THEN
          call DebugMemoryTest(STRINGOUT,'nContC*nPrimD*nPasses',nTUVP*nTUVQ,LUMOD3)
          call CallPrimitiveContractionString1(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'C')
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING          
          IF(AngmomP.EQ.0.AND.AngmomQ.EQ.0)THEN
             !no subsequent Horizontal or spherical transformations 
             IF(.NOT.OutputSet)THEN
                STRINGOUT  = 'LOCALINTS     ';OutputSet = .TRUE.
             ELSE
                STOP 'LOCALINTS already set MAJOER PROBLEM A2'
             ENDIF
          ELSE
             call DebugMemoryTest(STRINGOUT,'nContC*nContD*nPasses',nTUVP*nTUVQ,LUMOD3)
          ENDIF
          call CallPrimitiveContractionString1(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'D')
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING          
       ELSEIF(SegQ)THEN
          call DebugMemoryTest(STRINGOUT,'nContA*nPrimB*nPasses',nTUVP*nTUVQ,LUMOD3)
          call CallPrimitiveContractionString1(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'A')
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING          
          IF(AngmomP.EQ.0.AND.AngmomQ.EQ.0)THEN
             !no subsequent Horizontal or spherical transformations 
             IF(.NOT.OutputSet)THEN
                STRINGOUT  = 'LOCALINTS     ';OutputSet = .TRUE.
             ELSE
                STOP 'LOCALINTS already set MAJOER PROBLEM A2'
             ENDIF
          ELSE
             call DebugMemoryTest(STRINGOUT,'nContA*nContB*nPasses',nTUVP*nTUVQ,LUMOD3)
          ENDIF
          call CallPrimitiveContractionString1(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'B')
       ELSE
          call DebugMemoryTest(STRINGOUT,'nPrimA*nPrimB*nContC*nPrimD*nPasses',nTUVP*nTUVQ,LUMOD3)
          call CallPrimitiveContractionString2(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'C')
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING          
          call DebugMemoryTest(STRINGOUT,'nPrimA*nPrimB*nContC*nContD*nPasses',nTUVP*nTUVQ,LUMOD3)
          call CallPrimitiveContractionString2(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'D')
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING          
          call DebugMemoryTest(STRINGOUT,'nContA*nPrimB*nContC*nContD*nPasses',nTUVP*nTUVQ,LUMOD3)
          call CallPrimitiveContractionString2(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'A')
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING          
          IF(AngmomP.EQ.0.AND.AngmomQ.EQ.0)THEN
             !no subsequent Horizontal or spherical transformations 
             IF(.NOT.OutputSet)THEN
                STRINGOUT  = 'LOCALINTS     ';OutputSet = .TRUE.
             ELSE
                STOP 'LOCALINTS already set MAJOER PROBLEM A2'
             ENDIF
          ELSE
             call DebugMemoryTest(STRINGOUT,'nContA*nContB*nContC*nContD*nPasses',nTUVP*nTUVQ,LUMOD3)
          ENDIF
          call CallPrimitiveContractionString2(LUMOD3,ARCSTRING,BASISSPEC,iBasisSpec,nTUVQ,nTUVP,STRINGIN,STRINGOUT,'B')
       ENDIF
       TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING          
    ENDIF

    !================= DONE WITH PRIMITIVE CONTRACTION ================================================================

    !     The horizontal recurrence also extract nTUVspec from nTUV       
    IF(AngmomP.EQ.0)THEN
       WRITE(LUMOD3,'(A)')'        !no need for LHS Horizontal recurrence relations, it would be a simply copy'
       !there will not be need for spherical transformation afterwards
       IF(AngmomQ.EQ.0)THEN
          !there will not be need for RHS Horizontal recurrence relations nor Spherical Transformation'
          IF(.NOT.OutputSet)THEN
             WRITE(LUMOD3,'(A,A)')'        LOCALINTS = ',STRINGIN
             OutputSet = .TRUE.
          ENDIF
       ELSE
          !need for RHS Horizontal
       ENDIF
    ELSE
       IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
          !need for Spherical Transformation so we cannot place output in LOCALINTS yet
       ELSE
          !no need for LHS Spherical Transformation
          IF(AngmomQ.EQ.0)THEN
             !there will not be need for RHS Horizontal recurrence relations nor Spherical Transformation'
             IF(.NOT.OutputSet)THEN
                STRINGOUT  = 'LOCALINTS     '
                OutputSet = .TRUE.
             ELSE
                STOP 'LOCALINTS already set MAJOER PROBLEM A2'
             ENDIF
          ELSE
             !need for RHS Horizontal recurrence
          ENDIF
       ENDIF
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
       IF(Gen)THEN
          call DebugMemoryTest(STRINGOUT,'nContQP*nPasses',nTUVAspec*nTUVBspec*nTUVQ,LUMOD3)
       ELSEIF(SegQ)THEN
          call DebugMemoryTest(STRINGOUT,'nContP*nPasses',nTUVAspec*nTUVBspec*nTUVQ,LUMOD3)
       ELSEIF(SegP)THEN
          call DebugMemoryTest(STRINGOUT,'nContQ*nPasses',nTUVAspec*nTUVBspec*nTUVQ,LUMOD3)
       ELSE
          call DebugMemoryTest(STRINGOUT,'nPasses',nTUVAspec*nTUVBspec*nTUVQ,LUMOD3)
       ENDIF
       call initString(8)
       call AddToString('call HorizontalRR_'//ARCSTRING//'_LHS_P')
       call AddToString(AngmomP)
       call AddToString('A')
       call AddToString(AngmomA)
       call AddToString('B')
       call AddToString(AngmomB)
       call AddToString(SPEC)
       IF(Gen)THEN
          call AddToString('(nContQP,nPasses,')
       ELSEIF(SegQ)THEN
          call AddToString('(nContP,nPasses,')
       ELSEIF(SegP)THEN
          call AddToString('(nContQ,nPasses,')
       ELSE
          call AddToString('(1,nPasses,')
       ENDIF
       call AddToString(nTUVQ)
       call AddToString(',&')
       call writeString(LUMOD3)
       call initString(12)
       call AddToString('& Pdistance12,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass,')
       call AddToString(STRINGIN)
!!$       IF(Gen)THEN
!!$          call AddToString('(1:nContQP*nPasses*')
!!$          call AddToString(nTUVP*nTUVQ)
!!$          call AddToString(')')
!!$       ELSEIF(SegQ)THEN
!!$          call AddToString('(1:nContP*nPasses*')
!!$          call AddToString(nTUVP*nTUVQ)
!!$          call AddToString(')')
!!$       ELSEIF(SegP)THEN
!!$          call AddToString('(1:nContQ*nPasses*')
!!$          call AddToString(nTUVP*nTUVQ)
!!$          call AddToString(')')
!!$       ELSE
!!$          call AddToString('(1:nPasses*')
!!$          call AddToString(nTUVP*nTUVQ)
!!$          call AddToString(')')
!!$       ENDIF
       call AddToString(',&')
       call writeString(LUMOD3)
       call initString(12)
       call AddToString('& ')
       call AddToString(STRINGOUT)
!!$       IF(Gen)THEN
!!$          call AddToString('(1:nContQP*nPasses*')
!!$          call AddToString(nTUVAspec*nTUVBspec*nTUVQ)
!!$          call AddToString(')')
!!$       ELSEIF(SegQ)THEN
!!$          call AddToString('(1:nContP*nPasses*')
!!$          call AddToString(nTUVAspec*nTUVBspec*nTUVQ)
!!$          call AddToString(')')
!!$       ELSEIF(SegP)THEN
!!$          call AddToString('(1:nContQ*nPasses*')
!!$          call AddToString(nTUVAspec*nTUVBspec*nTUVQ)
!!$          call AddToString(')')
!!$       ELSE
!!$          call AddToString('(1:nPasses*')
!!$          call AddToString(nTUVAspec*nTUVBspec*nTUVQ)
!!$          call AddToString(')')
!!$       ENDIF
       call AddToString(',lupri')
       IF(CPU)THEN
          call AddToString(')')
       ELSE
          call AddToString(',iASync)')
       ENDIF
       call writeString(LUMOD3)
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF

    IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
       IF(AngmomQ.EQ.0)THEN
          !no need for RHS horizontal transfer nor spherical transformation
          IF(.NOT.OutputSet)THEN
             !there will not be need for RHS Horizontal recurrence relations nor Spherical Transformation'
             STRINGOUT  = 'LOCALINTS     '
             OutputSet = .TRUE.
          ELSE
             STOP 'LOCALINTS already set MAJOER PROBLEM B1'
          ENDIF
       ELSE
          !need for RHS horizontal recurrence
       ENDIF
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
       IF(Gen)THEN
          call DebugMemoryTest(STRINGOUT,'nContQP*nPasses',nlmA*nlmB*nTUVQ,LUMOD3)
       ELSEIF(SegQ)THEN
          call DebugMemoryTest(STRINGOUT,'nContP*nPasses',nlmA*nlmB*nTUVQ,LUMOD3)
       ELSEIF(SegP)THEN
          call DebugMemoryTest(STRINGOUT,'nContQ*nPasses',nlmA*nlmB*nTUVQ,LUMOD3)
       ELSE
          call DebugMemoryTest(STRINGOUT,'nPasses',nlmA*nlmB*nTUVQ,LUMOD3)
       ENDIF
       call initString(8)
       call AddToString('call SphericalContractOBS1_'//ARCSTRING//'_maxAngP')
       call AddToString(AngmomP)
       call AddToString('_maxAngA')
       call AddToString(AngmomA)
       call AddToString('(')
       call AddToString(nTUVQ)
       IF(Gen)THEN
          call AddToString(',nContQP*nPasses,')
       ELSEIF(SegP)THEN
          call AddToString(',nContQ*nPasses,')
       ELSEIF(SegQ)THEN
          call AddToString(',nContP*nPasses,')
       ELSE
          call AddToString(',nPasses,')
       ENDIF
       call AddToString(STRINGIN)
!       IF(Gen)THEN
!          call AddToString('(1:nContQP*nPasses*')
!          call AddToString(nTUVAspec*nTUVBspec*nTUVQ)
!          call AddToString(')')
!       ELSEIF(SegQ)THEN
!          call AddToString('(1:nContP*nPasses*')
!          call AddToString(nTUVAspec*nTUVBspec*nTUVQ)
!          call AddToString(')')
!       ELSEIF(SegP)THEN
!          call AddToString('(1:nContQ*nPasses*')
!          call AddToString(nTUVAspec*nTUVBspec*nTUVQ)
!          call AddToString(')')
!       ELSE
!          call AddToString('(1:nPasses*')
!          call AddToString(nTUVAspec*nTUVBspec*nTUVQ)
!          call AddToString(')')
!       ENDIF
       call AddToString(',&')
       call writeString(LUMOD3)
       call initString(12)
       call AddToString('& ')
       call AddToString(STRINGOUT)
!       IF(Gen)THEN
!          call AddToString('(1:nContQP*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVQ)
!          call AddToString(')')
!       ELSEIF(SegQ)THEN
!          call AddToString('(1:nContP*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVQ)
!          call AddToString(')')
!       ELSEIF(SegP)THEN
!          call AddToString('(1:nContQ*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVQ)
!          call AddToString(')')
!       ELSE
!          call AddToString('(1:nPasses*')
!          call AddToString(nlmA*nlmB*nTUVQ)
!          call AddToString(')')
!       ENDIF
       IF(CPU)THEN
          call AddToString(')')
       ELSE
          call AddToString(',iASync)')
       ENDIF

       call writeString(LUMOD3)
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ELSE
       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation LHS needed'
       IF(AngmomQ.EQ.0)THEN
          !there will not be need for RHS Horizontal recurrence relations nor Spherical Transformation'
          !afterwards which means we can 
          IF(.NOT.OutputSet)THEN
             WRITE(LUMOD3,'(A,A)')'        LOCALINTS = ',STRINGIN
             OutputSet = .TRUE.
          ENDIF
       ELSE
          !need for RHS Horizontal so no copy
       ENDIF
    ENDIF

    IF(AngmomQ.EQ.0)THEN
       WRITE(LUMOD3,'(A)')'        !no need for RHS Horizontal recurrence relations '
       !there will not be need for Spherical either 
       IF(.NOT.OutputSet)THEN
          WRITE(LUMOD3,'(A,A)')'        LOCALINTS = ',STRINGIN
          OutputSet = .TRUE.
          print*,'ouptutsat1'
       ENDIF
    ELSE
       !       WRITE(LUMOD3,'(A)')'        !RHS Horizontal recurrence relations '
       IF(.NOT.OutputSet)THEN
          IF(.NOT.(Spherical.AND.(AngmomC.GT.1.OR.AngmomD.GT.1)))THEN
             !no Spherical afterwards which means we can 
             STRINGOUT  = 'LOCALINTS     '
             OutputSet = .TRUE.
          ENDIF
       ELSE
          STOP 'MAJOR ERROR OUTPUT SET BUT RHS HORIZONTAL NEEDED C1'
       ENDIF
       IF(AngmomC.GE.AngmomD)THEN
          SPEC = 'CtoD'
       ELSE
          SPEC = 'DtoC'
       ENDIF
       IF(Gen)THEN
          call DebugMemoryTest(STRINGOUT,'nContQP*nPasses',nlmA*nlmB*nTUVCspec*nTUVDspec,LUMOD3)
       ELSEIF(SegQ)THEN
          call DebugMemoryTest(STRINGOUT,'nContP*nPasses',nlmA*nlmB*nTUVCspec*nTUVDspec,LUMOD3)
       ELSEIF(SegP)THEN
          call DebugMemoryTest(STRINGOUT,'nContQ*nPasses',nlmA*nlmB*nTUVCspec*nTUVDspec,LUMOD3)
       ELSE
          call DebugMemoryTest(STRINGOUT,'nPasses',nlmA*nlmB*nTUVCspec*nTUVDspec,LUMOD3)
       ENDIF
       call initString(8)
       call AddToString('call HorizontalRR_'//ARCSTRING//'_RHS_Q')
       call AddToString(AngmomQ)
       call AddToString('C')
       call AddToString(AngmomC)
       call AddToString('D')
       call AddToString(AngmomD)
       call AddToString(SPEC)
       IF(Gen)THEN
          call AddToString('(nContQP,nPasses,')
       ELSEIF(SegP)THEN
          call AddToString('(nContQ,nPasses,')
       ELSEIF(SegQ)THEN
          call AddToString('(nContP,nPasses,')
       ELSE
          call AddToString('(1,nPasses,')
       ENDIF
       call AddToString(nlmA*nlmB)
       call AddToString(',Qdistance12,')
       call AddToString(STRINGIN)
!       IF(Gen)THEN
!          call AddToString('(1:nContQP*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVQ)
!          call AddToString(')')
!       ELSEIF(SegQ)THEN
!          call AddToString('(1:nContP*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVQ)
!          call AddToString(')')
!       ELSEIF(SegP)THEN
!          call AddToString('(1:nContQ*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVQ)
!          call AddToString(')')
!       ELSE
!          call AddToString('(1:nPasses*')
!          call AddToString(nlmA*nlmB*nTUVQ)
!          call AddToString(')')
!       ENDIF
       call AddToString(',&')
       call writeString(LUMOD3)
       call initString(12)
       call AddToString('& ')
       call AddToString(STRINGOUT)
!       IF(Gen)THEN
!          call AddToString('(1:nContQP*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVCspec*nTUVDspec)
!          call AddToString(')')
!       ELSEIF(SegQ)THEN
!          call AddToString('(1:nContP*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVCspec*nTUVDspec)
!          call AddToString(')')
!       ELSEIF(SegP)THEN
!          call AddToString('(1:nContQ*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVCspec*nTUVDspec)
!          call AddToString(')')
!       ELSE
!          call AddToString('(1:nPasses*')
!          call AddToString(nlmA*nlmB*nTUVCspec*nTUVDspec)
!          call AddToString(')')
!       ENDIF
       call AddToString(',lupri')
       IF(CPU)THEN
          call AddToString(')')
       ELSE
          call AddToString(',iASync)')
       ENDIF
       call writeString(LUMOD3)
       !       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF

    IF(Spherical.AND.(AngmomC.GT.1.OR.AngmomD.GT.1))THEN
       !       WRITE(LUMOD3,'(A)')'        !Spherical Transformation RHS'
       IF(.NOT.OutputSet)THEN
          STRINGOUT  = 'LOCALINTS     '
       ELSE
          STOP 'MAJOR ERROR OUTPUT SET BUT RHS HORIZONTAL NEEDED D1'
       ENDIF
       IF(Gen)THEN
          call DebugMemoryTest(STRINGOUT,'nContQP*nPasses',nlmA*nlmB*nlmC*nlmD,LUMOD3)
       ELSEIF(SegQ)THEN
          call DebugMemoryTest(STRINGOUT,'nContP*nPasses',nlmA*nlmB*nlmC*nlmD,LUMOD3)
       ELSEIF(SegP)THEN
          call DebugMemoryTest(STRINGOUT,'nContQ*nPasses',nlmA*nlmB*nlmC*nlmD,LUMOD3)
       ELSE
          call DebugMemoryTest(STRINGOUT,'nPasses',nlmA*nlmB*nlmC*nlmD,LUMOD3)
       ENDIF
       call initString(8)
       call AddToString('call SphericalContractOBS2_'//ARCSTRING//'_maxAngQ')
       call AddToString(AngmomQ)
       call AddToString('_maxAngC')
       call AddToString(AngmomC)
       call AddToString('(')
       call AddToString(nlmA*nlmB)
       IF(Gen)THEN
          call AddToString(',nContQP*nPasses,')
       ELSEIF(SegQ)THEN
          call AddToString(',nContP*nPasses,')
       ELSEIF(SegP)THEN
          call AddToString(',nContQ*nPasses,')
       ELSE
          call AddToString(',nPasses,')
       ENDIF
       call AddToString(STRINGIN)
!       IF(Gen)THEN
!          call AddToString('(1:nContQP*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVCspec*nTUVDspec)
!          call AddToString(')')
!       ELSEIF(SegQ)THEN
!          call AddToString('(1:nContP*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVCspec*nTUVDspec)
!          call AddToString(')')
!       ELSEIF(SegP)THEN
!          call AddToString('(1:nContQ*nPasses*')
!          call AddToString(nlmA*nlmB*nTUVCspec*nTUVDspec)
!          call AddToString(')')
!       ELSE
!          call AddToString('(1:nPasses*')
!          call AddToString(nlmA*nlmB*nTUVCspec*nTUVDspec)
!          call AddToString(')')
!       ENDIF
       call AddToString(',&')
       call writeString(LUMOD3)
       call initString(12)
       call AddToString('& ')
       call AddToString(STRINGOUT)
 !      IF(Gen)THEN
 !         call AddToString('(1:nContQP*nPasses*')
 !         call AddToString(nlmA*nlmB*nlmC*nlmD)
 !         call AddToString(')')
 !      ELSEIF(SegQ)THEN
 !         call AddToString('(1:nContP*nPasses*')
 !         call AddToString(nlmA*nlmB*nlmC*nlmD)
 !         call AddToString(')')
 !      ELSEIF(SegP)THEN
 !         call AddToString('(1:nContQ*nPasses*')
 !         call AddToString(nlmA*nlmB*nlmC*nlmD)
 !         call AddToString(')')
 !      ELSE
 !         call AddToString('(1:nPasses*')
 !         call AddToString(nlmA*nlmB*nlmC*nlmD)
 !         call AddToString(')')
 !      ENDIF
       IF(CPU)THEN
          call AddToString(')')
       ELSE
          call AddToString(',iASync)')
       ENDIF
       call writeString(LUMOD3)
       !       WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
       !swap 
       !       TMPSTRING = STRINGIN
       !       STRINGIN  = STRINGOUT
       !       STRINGOUT  = TMPSTRING
       !       WRITE(LUMOD3,'(A,A)')'        LOCALINTS = ',STRINGIN
    ELSE
       WRITE(LUMOD3,'(A)')'        !no Spherical Transformation RHS needed'
    ENDIF

    !    WRITE(LUMOD3,'(A,A,A)')'        call mem_ichor_dealloc(',STRINGIN,')'
  end subroutine subroutineMain

  subroutine DebugMemoryTest(STRING,DIMSTRING,DIMINT,LUPRI2)
    implicit none
    integer :: LUPRI2
    character(len=9),intent(in) :: STRING
    character*(*) :: DIMSTRING
    integer :: DIMINT
    IF(STRING(1:4).NE.'LOCALINTS')THEN
       WRITE(LUPRI2,'(A)')'#ifdef VAR_DEBUGICHOR'
       call initString(8)
       call AddToString('IF(')
       call AddToString(DIMSTRING)
       call AddToString('*')
       call AddToString(DIMINT)
       call AddToString('.GT.')                
       call AddToString(STRING)
       call AddToString('maxsize)THEN')
       call writeString(LUPRI2)
       
       call initString(10)
       call AddToString('call ichorquit(''')
       call AddToString(DIMSTRING)
       call AddToString(' too small'',-1)')
       call writeString(LUPRI2)
       WRITE(LUPRI2,'(A)')'        ENDIF'
       WRITE(LUPRI2,'(A)')'#endif'
    ENDIF
  end subroutine DebugMemoryTest

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
       & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,nlmA,nlmB,nlmC,nlmD,spherical,Gen,SegQ,SegP,Seg,Seg1Prim)
    implicit none
    integer,intent(in) :: LUMOD3,AngmomA,AngmomB,AngmomC,AngmomD,nTUV,AngmomP,AngmomQ
    integer,intent(in) :: AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec
    integer,intent(in) :: nlmA,nlmB,nlmC,nlmD
    character(len=9) :: STRINGIN,STRINGOUT,TMPSTRING
    logical :: spherical,OutputSet,Gen,SegQ,SegP,Seg,Seg1Prim,Contracted
    character(len=8) :: BASISSPEC
    integer :: iBasisSpec
    logical :: PerformSphericaQAndPlaceInTmp
    logical :: PerformHorizontQAndPlaceInTmp
    logical :: PerformSphericaPAndPlaceInTmp
    logical :: PerformHorizontPAndPlaceInTmp
    logical :: PerformBasisContAndPlaceInTmp
    logical :: PerformTranserAndPlaceInTmp
    logical :: PerformVerticalAndPlaceInTmp
    logical :: PerformBasisContLAST

    IF(Gen)THEN
       iBasisSpec = 3
       BASISSPEC = 'Gen     '
    ELSEIF(Seg)THEN
       iBasisSpec = 4
       BASISSPEC = 'SegQ    '
    ELSEIF(SegQ)THEN
       iBasisSpec = 4
       BASISSPEC = 'SegP    '
    ELSEIF(SegP)THEN
       iBasisSpec = 3
       BASISSPEC = 'Seg     '
    ELSEIF(Seg1Prim)THEN
       iBasisSpec = 8
       BASISSPEC = 'Seg1Prim'
    ENDIF
    OutputSet = .FALSE.
    Contracted = .FALSE.
    PerformBasisContLAST = .FALSE. 

    PerformSphericaQAndPlaceInTmp = .FALSE. !always false as placed in LOCALINTS 
    IF(Spherical.AND.(AngmomC.GT.1.OR.AngmomD.GT.1))THEN
       !RHS Spherical Transformation placed in LOCALINTS
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
       IF(Seg.OR.Seg1Prim)THEN
          PerformBasisContAndPlaceInTmp = .FALSE. 
       ELSE
          PerformBasisContAndPlaceInTmp = .TRUE. 
       ENDIF
       IF(AngmomP.GT.0.AND.AngmomQ.GT.0)THEN
          PerformTranserAndPlaceInTmp = .TRUE. 
       ELSE
          PerformTranserAndPlaceInTmp = .FALSE. 
       ENDIF
       PerformVerticalAndPlaceInTmp = .TRUE. 
    ELSE
       PerformHorizontQAndPlaceInTmp = .FALSE.        
       IF(.NOT.(AngmomC.EQ.0.AND.AngmomD.EQ.0))THEN
          !RHS Horizontal Recurrence Relation place in LOCALINTS
          IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
             PerformSphericaPAndPlaceInTmp = .TRUE. 
          ELSE
             PerformSphericaPAndPlaceInTmp = .FALSE. 
          ENDIF
          IF(AngmomA.EQ.0.AND.AngmomB.EQ.0)THEN
             PerformHorizontPAndPlaceInTmp = .FALSE. 
          ELSE
             !LHS Horizontal Recurrence Relation
             PerformHorizontPAndPlaceInTmp = .TRUE. 
          ENDIF
          IF(Seg.OR.Seg1Prim)THEN
             PerformBasisContAndPlaceInTmp = .FALSE. 
          ELSE
             PerformBasisContAndPlaceInTmp = .TRUE. 
          ENDIF
          IF(AngmomP.GT.0.AND.AngmomQ.GT.0)THEN
             PerformTranserAndPlaceInTmp = .TRUE. 
          ELSE
             PerformTranserAndPlaceInTmp = .FALSE. 
          ENDIF
          PerformVerticalAndPlaceInTmp = .TRUE. 
       ELSE
          PerformSphericaPAndPlaceInTmp = .FALSE. 
          IF(Spherical.AND.(AngmomA.GT.1.OR.AngmomB.GT.1))THEN
             !LHS Spherical Transformation placed in LOCALINTS
             IF(AngmomA.EQ.0.AND.AngmomB.EQ.0)THEN
                PerformHorizontPAndPlaceInTmp = .FALSE. 
             ELSE
                !LHS Horizontal Recurrence Relation
                PerformHorizontPAndPlaceInTmp = .TRUE. 
             ENDIF
             IF(Seg.OR.Seg1Prim)THEN
                PerformBasisContAndPlaceInTmp = .FALSE. 
             ELSE
                PerformBasisContAndPlaceInTmp = .TRUE. 
             ENDIF
             IF(AngmomP.GT.0.AND.AngmomQ.GT.0)THEN
                PerformTranserAndPlaceInTmp = .TRUE. 
             ELSE
                PerformTranserAndPlaceInTmp = .FALSE. 
             ENDIF
             PerformVerticalAndPlaceInTmp = .TRUE.             
          ELSE
             PerformHorizontPAndPlaceInTmp = .FALSE. 
             IF(.NOT.(AngmomA.EQ.0.AND.AngmomB.EQ.0))THEN
                !LHS Horizontal Recurrence Relation placed in LOCALINTS
                IF(Seg.OR.Seg1Prim)THEN
                   PerformBasisContAndPlaceInTmp = .FALSE. 
                ELSE
                   PerformBasisContAndPlaceInTmp = .TRUE. 
                ENDIF
                IF(AngmomP.GT.0.AND.AngmomQ.GT.0)THEN
                   PerformTranserAndPlaceInTmp = .TRUE. 
                ELSE
                   PerformTranserAndPlaceInTmp = .FALSE. 
                ENDIF
                PerformVerticalAndPlaceInTmp = .TRUE.                 
             ELSE
                PerformBasisContAndPlaceInTmp = .FALSE. 
                PerformBasisContLAST = .TRUE. 
                !Basis set contraction last step 
                IF(.NOT.(Seg.OR.Seg1Prim))THEN
                   !Basis set Contraction placed in LOCALINTS
                   IF(AngmomP.GT.0.AND.AngmomQ.GT.0)THEN
                      PerformTranserAndPlaceInTmp = .TRUE. 
                   ELSE
                      PerformTranserAndPlaceInTmp = .FALSE. 
                   ENDIF
                   PerformVerticalAndPlaceInTmp = .TRUE.                    
                ELSE
                   PerformTranserAndPlaceInTmp = .FALSE. 
                   IF(AngmomP.GT.0.AND.AngmomQ.GT.0)THEN
                      !Electron Transfer placed in LOCALINTS
                      PerformVerticalAndPlaceInTmp = .TRUE.                    
                   ELSE
                      !Vertical Transfer placed in LOCALINTS
                      PerformVerticalAndPlaceInTmp = .FALSE.
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    IF(PerformVerticalAndPlaceInTmp)THEN
       !       WRITE(LUMOD3,'(A)')'       ! Vertical Recurrence'
       IF(AngmomPQ.EQ.0)THEN
          call initString(7)
          call AddToString(STRINGOUT)
          call AddToString('maxSize = MAX(')
          call AddToString(STRINGOUT)
          call AddToString('maxSize,')
          call AddToString(nTUV)       
          !          IF(PerformTranserAndPlaceInTmp)THEN
          !             call AddToString('*nPrimQP)')
          !          ELSE 
          IF(Gen)THEN
             call AddToString('*nPrimQP)')
          ELSEIF(SegQ)THEN
             call AddToString('*nPrimP)')
          ELSEIF(SegP)THEN
             call AddToString('*nPrimQ)')
          ELSEIF(Seg)THEN
             call AddToString(')')
          ELSEIF(Seg1Prim)THEN
             call AddToString(')')
          ENDIF
          !          ENDIF
          call writeString(LUMOD3)                
          !swap 
          TMPSTRING = STRINGIN
          STRINGIN  = STRINGOUT
          STRINGOUT  = TMPSTRING
       ELSE
          IF(AngmomPQ.GT.1)THEN
             call initString(7)
             call AddToString(STRINGOUT)
             call AddToString('maxSize = MAX(')
             call AddToString(STRINGOUT)
             call AddToString('maxSize,')
             call AddToString(AngmomPQ+1)       
             IF(PerformTranserAndPlaceInTmp)THEN
                call AddToString('*nPrimQP)')
             ELSE 
                IF(Seg1Prim)THEN
                   call AddToString(')')
                ELSE
                   call AddToString('*nPrimQP)')
                ENDIF
             ENDIF
             call writeString(LUMOD3)                
             !swap 
             TMPSTRING = STRINGIN
             STRINGIN  = STRINGOUT
             STRINGOUT  = TMPSTRING
          ENDIF
          call initString(7)
          call AddToString(STRINGOUT)
          call AddToString('maxSize = MAX(')
          call AddToString(STRINGOUT)
          call AddToString('maxSize,')
          call AddToString(nTUV)
          IF(PerformTranserAndPlaceInTmp)THEN
             call AddToString('*nPrimQP)')
          ELSE 
             IF(Gen)THEN
                call AddToString('*nPrimQP)')
             ELSEIF(SegQ)THEN
                call AddToString('*nPrimP)')
             ELSEIF(SegP)THEN
                call AddToString('*nPrimQ)')
             ELSEIF(Seg)THEN
                call AddToString(')')
             ELSEIF(Seg1Prim)THEN
                call AddToString(')')
             ENDIF
          ENDIF
          call writeString(LUMOD3)                
          !swap 
          TMPSTRING = STRINGIN
          STRINGIN  = STRINGOUT
          STRINGOUT  = TMPSTRING
       ENDIF
    ENDIF
    IF(PerformTranserAndPlaceInTmp)THEN
!       WRITE(LUMOD3,'(A)')'       ! Transfer Recurrence'
       call initString(7)
       call AddToString(STRINGOUT)
       call AddToString('maxSize = MAX(')
       call AddToString(STRINGOUT)
       call AddToString('maxSize,')
       call AddToString(nTUVQ*nTUVP)       
       IF(Gen)THEN
          call AddToString('*nPrimQP)')
       ELSEIF(SegQ)THEN
          call AddToString('*nPrimP)')
       ELSEIF(SegP)THEN
          call AddToString('*nPrimQ)')
       ELSEIF(Seg)THEN
          call AddToString(')')
       ELSEIF(Seg1Prim)THEN
          call AddToString(')')
       ENDIF
       call writeString(LUMOD3)                
       !swap 
       TMPSTRING = STRINGIN
       STRINGIN  = STRINGOUT
       STRINGOUT  = TMPSTRING
    ENDIF
    IF(PerformBasisContAndPlaceInTmp.OR.PerformBasisContLAST)THEN
       !       WRITE(LUMOD3,'(A)')'       ! BasisCont Recurrence'
       IF(Gen)THEN
          call initString(4)
          call AddToString(STRINGOUT);call AddToString('maxSize = MAX(')
          call AddToString(STRINGOUT);call AddToString('maxSize,');call AddToString(nTUVQ*nTUVP)       
          call AddToString('*nContC*nPrimD*nPrimA*nPrimB)')          
          call writeString(LUMOD3)                
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING
          call initString(4)
          call AddToString(STRINGOUT);call AddToString('maxSize = MAX(')
          call AddToString(STRINGOUT);call AddToString('maxSize,');call AddToString(nTUVQ*nTUVP)       
          call AddToString('*nContC*nContD*nPrimA*nPrimB)')          
          call writeString(LUMOD3)                
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING
          call initString(4)
          call AddToString(STRINGOUT);call AddToString('maxSize = MAX(')
          call AddToString(STRINGOUT);call AddToString('maxSize,');call AddToString(nTUVQ*nTUVP)       
          call AddToString('*nContC*nContD*nContA*nPrimB)')          
          call writeString(LUMOD3)                
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING
          IF(.NOT.PerformBasisContLAST)THEN
             !place in Output not Tmparray
             call initString(4)
             call AddToString(STRINGOUT);call AddToString('maxSize = MAX(')
             call AddToString(STRINGOUT);call AddToString('maxSize,');call AddToString(nTUVQ*nTUVP)       
             call AddToString('*nContC*nContD*nContA*nContB)')          
             call writeString(LUMOD3)                
             TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING
          ENDIF
       ELSEIF(SegQ)THEN
          call initString(4)
          call AddToString(STRINGOUT);call AddToString('maxSize = MAX(')
          call AddToString(STRINGOUT);call AddToString('maxSize,');call AddToString(nTUVQ*nTUVP)       
          call AddToString('*nContA*nPrimB)')
          call writeString(LUMOD3)                
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING
          IF(.NOT.PerformBasisContLAST)THEN
             !place in Output not Tmparray
             call initString(4)
             call AddToString(STRINGOUT);call AddToString('maxSize = MAX(')
             call AddToString(STRINGOUT);call AddToString('maxSize,');call AddToString(nTUVQ*nTUVP)       
             call AddToString('*nContA*nContB)')
             call writeString(LUMOD3)                
             TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING
          ENDIF
       ELSEIF(SegP)THEN
          call initString(4)
          call AddToString(STRINGOUT);call AddToString('maxSize = MAX(')
          call AddToString(STRINGOUT);call AddToString('maxSize,');call AddToString(nTUVQ*nTUVP)       
          call AddToString('*nContC*nPrimD)')
          call writeString(LUMOD3)                
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING
          IF(.NOT.PerformBasisContLAST)THEN
             !place in Output not Tmparray
             call initString(4)
             call AddToString(STRINGOUT);call AddToString('maxSize = MAX(')
             call AddToString(STRINGOUT);call AddToString('maxSize,');call AddToString(nTUVQ*nTUVP)       
             call AddToString('*nContC*nContD)')
             call writeString(LUMOD3)                
             TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING
          ENDIF
       ELSE
!          STOP 'NOT CORRECT basis should nbot be called PerformBasisContAndPlaceInTmp'
          call initString(4)
          call AddToString(STRINGOUT);call AddToString('maxSize = MAX(')
          call AddToString(STRINGOUT);call AddToString('maxSize,');call AddToString(nTUVQ*nTUVP)       
          call AddToString(')')
          call writeString(LUMOD3)                
          TMPSTRING = STRINGIN; STRINGIN  = STRINGOUT; STRINGOUT  = TMPSTRING
       ENDIF
    ENDIF
    IF(PerformHorizontPAndPlaceInTmp)THEN
!       WRITE(LUMOD3,'(A)')'       ! Horizontal LHS Recurrence'
       call initString(7)
       call AddToString(STRINGOUT)
       call AddToString('maxSize = MAX(')
       call AddToString(STRINGOUT)
       call AddToString('maxSize,')
       call AddToString(nTUVQ*nTUVAspec*nTUVBspec)       
       IF(Gen)THEN
          call AddToString('*nContQP)')
       ELSEIF(SegQ)THEN
          call AddToString('*nContP)')
       ELSEIF(SegP)THEN
          call AddToString('*nContQ)')
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
          call AddToString('*nContQP)')
       ELSEIF(SegQ)THEN
          call AddToString('*nContP)')
       ELSEIF(SegP)THEN
          call AddToString('*nContQ)')
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
          call AddToString('*nContQP)')
       ELSEIF(SegQ)THEN
          call AddToString('*nContP)')
       ELSEIF(SegP)THEN
          call AddToString('*nContQ)')
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

