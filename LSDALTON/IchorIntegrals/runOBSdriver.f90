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
  integer :: LUMOD5,LUMOD6,LUMOD7,LUMOD8,LUMOD9,LUMOD10,LUMOD11,LUMOD12
  integer :: JMAX2,JMAX3,JMAX4,ituvP,jp,tp,up,vp,ntuvP,l1,l2,l12
  integer :: ip1,jp1,kp1,p1,ip2,jp2,kp2,p2,ijkcartP,AngmomA,AngmomB
  integer :: AngmomC,AngmomD,AngmomP,AngmomQ,nTUVQ,AngmomPQ,IFILES
  integer :: nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec
  integer :: nTUVA,nTUVB,nTUVC,nTUVD,GPUrun,MaxAngmomSpecial
  integer :: nlmA,nlmB,nlmC,nlmD,angmomID,iseg,ILUMOD,I,nUniquenTUVs
  real(realk),pointer :: uniqeparam(:)
  character(len=15),pointer :: uniqeparamNAME(:)
  character(len=12) :: STRINGIN,STRINGOUT,TMPSTRING
  character(len=9) :: STRINGIN2,STRINGOUT2,TMPSTRING2
  character(len=4) :: SPEC
  character(len=3) :: ARCSTRING
  logical :: BUILD(0:2,0:2,0:2,0:2),Gen,Seg,SegP,segQ,Seg1Prim,UNIQUE,CPU
  integer,pointer :: UniquenTUVs(:)

  MaxAngmomSpecial = 2 !D functions
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

  LUMOD8=8
  open(unit = LUMOD8, file='MAIN_'//ARCSTRING//'_OBS_DRIVERGen2.f90',status="unknown")
  WRITE(LUMOD8,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGen2'
  WRITE(LUMOD8,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD8,'(A)')'!Contains routines for General Contracted Basisset '

  LUMOD9=9
  open(unit = LUMOD9, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSegQ2.f90',status="unknown")
  WRITE(LUMOD9,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQ2'
  WRITE(LUMOD9,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD9,'(A)')'!Contains routines for a General Contracted LHS Segmented contracted RHS and Basisset '

  LUMOD10=10
  open(unit = LUMOD10, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSegP2.f90',status="unknown")
  WRITE(LUMOD10,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegP2'
  WRITE(LUMOD10,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD10,'(A)')'!Contains routines for Segmented contracted LHS and a General Contracted RHS Basisset '

  LUMOD11=11
  open(unit = LUMOD11, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSeg2.f90',status="unknown")
  WRITE(LUMOD11,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg2'
  WRITE(LUMOD11,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD11,'(A)')'!Contains routines for Segmented contracted Basisset '

  LUMOD12=12
  open(unit = LUMOD12, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSeg1Prim2.f90',status="unknown")
  WRITE(LUMOD12,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1Prim2'
  WRITE(LUMOD12,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD12,'(A)')'!Contains routines for Segmented contracted Basisset containing a single primitive'

  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGen'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQ'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegP'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1Prim'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGen2'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQ2'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegP2'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg2'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1Prim2'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGenSize'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQSize'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegPSize'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegSize'
  WRITE(LUMOD2,'(A)')'use IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1PrimSize'

!  WRITE(LUMOD2,'(A)')'use SPIchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGen'
!  WRITE(LUMOD2,'(A)')'use SPIchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQ'
!  WRITE(LUMOD2,'(A)')'use SPIchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegP'
!  WRITE(LUMOD2,'(A)')'use SPIchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg'
!  WRITE(LUMOD2,'(A)')'use SPIchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1Prim'
!  WRITE(LUMOD2,'(A)')'use SPIchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGen2'
!  WRITE(LUMOD2,'(A)')'use SPIchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQ2'
!  WRITE(LUMOD2,'(A)')'use SPIchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegP2'
!  WRITE(LUMOD2,'(A)')'use SPIchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg2'
!  WRITE(LUMOD2,'(A)')'use SPIchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1Prim2'

!  IF(CPU)THEN
     DO ILUMOD=8,12
        WRITE(ILUMOD,'(A)')'use IchorEriCoulombintegralCPUMcMGeneralMod'
     ENDDO
!  ENDIF

  DO ILUMOD=2,12
     WRITE(ILUMOD,'(A)')'use IchorprecisionMod'
     WRITE(ILUMOD,'(A)')'use IchorCommonMod'
     WRITE(ILUMOD,'(A)')'use IchorMemory'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_BUILDRJ000ModGen'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_BUILDRJ000ModSeg1Prim'
  ENDDO
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODAGen'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASegQ'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASegP'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASeg'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASeg1Prim'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODAGen'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASegQ'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASegP'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASeg'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODASeg1Prim'

  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_PrimContractGenMod'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_PrimContractSegQMod'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_PrimContractSegPMod'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_PrimContractGenMod'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_PrimContractSegQMod'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_PrimContractSegPMod'

  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBGen'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSegQ'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSegP'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSeg'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSeg1Prim'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBGen'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSegQ'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSegP'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSeg'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBSeg1Prim'


  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDGen'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSegQ'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSegP'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSeg'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSeg1Prim'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDGen'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSegQ'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSegP'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSeg'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDSeg1Prim'

  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCGen'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSegQ'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSegP'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSeg'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSeg1Prim'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCGen'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSegQ'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSegP'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSeg'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCSeg1Prim'
  !also needed by the Segs
  DO ILUMOD=4,12
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODAGen'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODBGen'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODCGen'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_VERTICALRECURRENCEMODDGen'
  ENDDO
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCGen1'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCGen2'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDGen1'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDGen2'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCGen1'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDGen1'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoAGen'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoAGen'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBGen'
  WRITE(LUMOD3,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBGen'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegQ1'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegQ2'
!  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegQ3'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSegQ1'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSegQ2'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSegQ1'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSegQ1'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASegQ'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASegQ'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSegQ'
  WRITE(LUMOD4,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSegQ'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegP1'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegP2'
!  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegP3'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSegP1'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSegP2'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSegP1'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSegP1'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASegP'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASegP'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSegP'
  WRITE(LUMOD5,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSegP'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg1'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg2'
!  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg3'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSeg1'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSeg2'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSeg1'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSeg1'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASeg'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASeg'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSeg'
  WRITE(LUMOD6,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSeg'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg1Prim1'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg1Prim2'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSeg1Prim1'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSeg1Prim2'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSeg1Prim1'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSeg1Prim1'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASeg1Prim'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASeg1Prim'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSeg1Prim'
  WRITE(LUMOD7,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSeg1Prim'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCGen1'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCGen2'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDGen1'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDGen2'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCGen1'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDGen1'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoAGen'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoAGen'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBGen'
  WRITE(LUMOD8,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBGen'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegQ1'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegQ2'
!  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegQ3'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSegQ1'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSegQ2'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSegQ1'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSegQ1'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASegQ'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASegQ'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSegQ'
  WRITE(LUMOD9,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSegQ'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegP1'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegP2'
!  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSegP3'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSegP1'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSegP2'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSegP1'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSegP1'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASegP'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASegP'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSegP'
  WRITE(LUMOD10,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSegP'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg1'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg2'
!  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg3'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSeg1'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSeg2'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSeg1'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSeg1'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASeg'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASeg'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSeg'
  WRITE(LUMOD11,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSeg'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg1Prim1'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoCSeg1Prim2'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSeg1Prim1'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODAtoDSeg1Prim2'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoCSeg1Prim1'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODBtoDSeg1Prim1'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoASeg1Prim'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoASeg1Prim'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODCtoBSeg1Prim'
  WRITE(LUMOD12,'(A)')'use AGC_'//ARCSTRING//'_OBS_TRMODDtoBSeg1Prim'
  DO ILUMOD=3,12
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_HorizontalRecurrenceLHSModAtoB'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_HorizontalRecurrenceLHSModBtoA'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_HorizontalRecurrenceRHSModCtoD'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_HorizontalRecurrenceRHSModDtoC'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_Sphcontract1Mod'
     WRITE(ILUMOD,'(A)')'use AGC_'//ARCSTRING//'_OBS_Sphcontract2Mod'
     WRITE(ILUMOD,'(A)')'  '
     WRITE(ILUMOD,'(A)')'private   '
  ENDDO
  WRITE(LUMOD2,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_general'
  WRITE(LUMOD3,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_Gen'
  WRITE(LUMOD4,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_SegQ'
  WRITE(LUMOD5,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_SegP'
  WRITE(LUMOD6,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_Seg'
  WRITE(LUMOD7,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_Seg1Prim'
  WRITE(LUMOD8,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_Gen2'
  WRITE(LUMOD9,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_SegQ2'
  WRITE(LUMOD10,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_SegP2'
  WRITE(LUMOD11,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_Seg2'
  WRITE(LUMOD12,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_Seg1Prim2'

  DO ILUMOD=2,12
     WRITE(ILUMOD,'(A)')'  '
     WRITE(ILUMOD,'(A)')'CONTAINS'
     WRITE(ILUMOD,'(A)')'  '
     WRITE(ILUMOD,'(A)')'  '
  ENDDO
  WRITE(LUMOD2,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
  WRITE(LUMOD2,'(A)')'       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&'
  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&'
  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&'
  WRITE(LUMOD2,'(A)')'       & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&'
  WRITE(LUMOD2,'(A)')'       & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&'
  IF(CPU)THEN
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,UseSP)'
  ELSE
     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync,UseSP)'
  ENDIF
  WRITE(LUMOD2,'(A)')'    implicit none'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nPrimQP,MaxPasses,IntPrint,lupri'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nAtomsA,nAtomsB'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ)'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD'
  WRITE(LUMOD2,'(A)')'    integer,intent(in) :: nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV'
  WRITE(LUMOD2,'(A)')'    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)'
  WRITE(LUMOD2,'(A)')'    logical,intent(in)     :: Qsegmented,Psegmented,UseSP'
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
!!$  WRITE(LUMOD2,'(A)')'    real(reals),allocatable :: Aexp_SP(:),Bexp_SP(:),Cexp_SP(:),Dexp_SP(:)'
!!$  WRITE(LUMOD2,'(A)')'    real(reals),allocatable :: pexp_SP(:),qexp_SP(:),pcent_SP(:),qcent_SP(:)'
!!$  WRITE(LUMOD2,'(A)')'    real(reals),allocatable :: QpreExpFac_SP(:),PpreExpFac_SP(:),TABFJW_SP(:,:)'
!!$  WRITE(LUMOD2,'(A)')'    real(reals),allocatable :: ACC_SP(:,:),BCC_SP(:,:),CCC_SP(:,:),DCC_SP(:,:)'
!!$  WRITE(LUMOD2,'(A)')'    real(reals),allocatable :: LOCALINTS_SP(:),integralPrefactor_SP(:)'
!!$  WRITE(LUMOD2,'(A)')'    real(reals),allocatable :: reducedExponents_SP(:)'
!!$  WRITE(LUMOD2,'(A)')'    real(reals) :: Qdistance12_SP(3),Ccenter_SP(3),Dcenter_SP(3)'
!!$  WRITE(LUMOD2,'(A)')'    real(reals),allocatable :: Pdistance12_SP(:)'
!!$  WRITE(LUMOD2,'(A)')'    real(reals),allocatable :: Acenter_SP(:,:),Bcenter_SP(:,:)'
!!$  WRITE(LUMOD2,'(A)')'    real(reals),allocatable :: TmpArray1_SP(:),TmpArray2_SP(:)'


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
!  WRITE(LUMOD2,'(A)')'  IF(.NOT.UseSP)THEN'
  WRITE(LUMOD2,'(A,I1,A)')'   IF(.NOT.UseGeneralCode.AND.(((AngmomA.LE.',MaxAngmomSpecial,').AND.(AngmomA.GE.AngmomB)).AND.&'
  WRITE(LUMOD2,'(A)')'       ((AngmomA.GE.AngmomC).AND.(AngmomC.GE.AngmomD))))THEN'
  WRITE(LUMOD2,'(A)')'    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
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
  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented.AND.Qsegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_Seg(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
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
  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
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
  WRITE(LUMOD2,'(A)')'    ELSEIF(Qsegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
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

  WRITE(LUMOD2,'(A)')'    ELSE'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
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

  WRITE(LUMOD2,'(A)')'    ENDIF'
  WRITE(LUMOD2,'(A)')'   ELSE' !OTHER OPTIONS
  WRITE(LUMOD2,'(A)')'    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_Seg1Prim2(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
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
  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented.AND.Qsegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_Seg2(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
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
  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_SegP2(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
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
  WRITE(LUMOD2,'(A)')'    ELSEIF(Qsegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_SegQ2(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
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

  WRITE(LUMOD2,'(A)')'    ELSE'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_Gen2(nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
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

  WRITE(LUMOD2,'(A)')'    ENDIF'
  WRITE(LUMOD2,'(A)')'   ENDIF'
!!$  WRITE(LUMOD2,'(A)')'   ELSE !Single Precision Routines!'
!!$  WRITE(LUMOD2,'(A)')'  allocate(Aexp_SP(nPrimA))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(Bexp_SP(nPrimB))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(Cexp_SP(nPrimC))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(Dexp_SP(nPrimD))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(pexp_SP(nPrimP))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(qexp_SP(nPrimQ))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(pcent_SP(3*nPrimP*nAtomsA*nAtomsB))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(qcent_SP(3*nPrimQ))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(QpreExpFac_SP(nPrimQ))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(PpreExpFac_SP(nPrimP*nAtomsA*nAtomsB))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(TABFJW_SP(0:nTABFJW1,0:nTABFJW2))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(ACC_SP(nPrimA,nContA))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(BCC_SP(nPrimB,nContB))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(CCC_SP(nPrimC,nContC))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(DCC_SP(nPrimD,nContD))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(LOCALINTS_SP(localintsmaxsize))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(integralPrefactor_SP(nPrimQP))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(reducedExponents_SP(nPrimQP))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(Pdistance12_SP(3*nAtomsA*nAtomsB))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(Acenter_SP(3,nAtomsA))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(Bcenter_SP(3,nAtomsB))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(TmpArray1_SP(TMParray1maxsize))'
!!$  WRITE(LUMOD2,'(A)')'  allocate(TmpArray2_SP(TMParray2maxsize))'
!!$
!!$  WRITE(LUMOD2,'(A)')'  Aexp_SP = Real(Aexp,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  Bexp_SP = Real(Bexp,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  Cexp_SP = Real(Cexp,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  Dexp_SP = Real(Dexp,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  pexp_SP = Real(pexp,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  qexp_SP = Real(qexp,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  pcent_SP = Real(pcent,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  qcent_SP = Real(qcent,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  QpreExpFac_SP = Real(QpreExpFac,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  PpreExpFac_SP = Real(PpreExpFac,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  TABFJW_SP = Real(TABFJW,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  ACC_SP = Real(ACC,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  BCC_SP = Real(BCC,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  CCC_SP = Real(CCC,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  DCC_SP = Real(DCC,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  integralPrefactor_SP = Real(integralPrefactor,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  reducedExponents_SP = Real(reducedExponents,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  Pdistance12_SP = Real(Pdistance12,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  Acenter_SP = Real(Acenter,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  Bcenter_SP = Real(Bcenter,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  Qdistance12_SP = Real(Qdistance12,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  Ccenter_SP = Real(Ccenter,kind=reals)'
!!$  WRITE(LUMOD2,'(A)')'  Dcenter_SP = Real(Dcenter,kind=reals)'
!!$  WRITE(LUMOD2,'(A,I1,A)')'   IF(.NOT.UseGeneralCode.AND.(((AngmomA.LE.',MaxAngmomSpecial,').AND.(AngmomA.GE.AngmomB)).AND.&'
!!$  WRITE(LUMOD2,'(A)')'       ((AngmomA.GE.AngmomC).AND.(AngmomC.GE.AngmomD))))THEN'
!!$  WRITE(LUMOD2,'(A)')'    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN'
!!$  WRITE(LUMOD2,'(A)')'     call SPICI_'//ARCSTRING//'_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
!!$  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
!!$  WRITE(LUMOD2,'(A)')'       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & PQorder,LOCALINTS_SP,localintsmaxsize,&'
!!$  WRITE(LUMOD2,'(A)')'       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&'
!!$  WRITE(LUMOD2,'(A)')'       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&'
!!$  IF(CPU)THEN
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
!!$  ELSE
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
!!$  ENDIF
!!$  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented.AND.Qsegmented)THEN'
!!$  WRITE(LUMOD2,'(A)')'     call SPICI_'//ARCSTRING//'_OBS_Seg(nPrimA,nPrimB,nPrimC,nPrimD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
!!$  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
!!$  WRITE(LUMOD2,'(A)')'       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & PQorder,LOCALINTS_SP,localintsmaxsize,&'
!!$  WRITE(LUMOD2,'(A)')'       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&'
!!$  WRITE(LUMOD2,'(A)')'       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&'
!!$  IF(CPU)THEN
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
!!$  ELSE
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
!!$  ENDIF
!!$  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented)THEN'
!!$  WRITE(LUMOD2,'(A)')'     call SPICI_'//ARCSTRING//'_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
!!$  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
!!$  WRITE(LUMOD2,'(A)')'       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & PQorder,LOCALINTS_SP,localintsmaxsize,&'
!!$  WRITE(LUMOD2,'(A)')'       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&'
!!$  WRITE(LUMOD2,'(A)')'       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&'
!!$  IF(CPU)THEN
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
!!$  ELSE
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
!!$  ENDIF
!!$  WRITE(LUMOD2,'(A)')'    ELSEIF(Qsegmented)THEN'
!!$  WRITE(LUMOD2,'(A)')'     call SPICI_'//ARCSTRING//'_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
!!$  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
!!$  WRITE(LUMOD2,'(A)')'       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & PQorder,LOCALINTS_SP,localintsmaxsize,&'
!!$  WRITE(LUMOD2,'(A)')'       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&'
!!$  WRITE(LUMOD2,'(A)')'       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&'
!!$  IF(CPU)THEN
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
!!$  ELSE
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
!!$  ENDIF
!!$
!!$  WRITE(LUMOD2,'(A)')'    ELSE'
!!$  WRITE(LUMOD2,'(A)')'     call SPICI_'//ARCSTRING//'_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
!!$  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
!!$  WRITE(LUMOD2,'(A)')'       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & PQorder,LOCALINTS_SP,localintsmaxsize,&'
!!$  WRITE(LUMOD2,'(A)')'       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&'
!!$  WRITE(LUMOD2,'(A)')'       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&'
!!$  IF(CPU)THEN
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
!!$  ELSE
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
!!$  ENDIF
!!$
!!$  WRITE(LUMOD2,'(A)')'    ENDIF'
!!$  WRITE(LUMOD2,'(A)')'   ELSE' !OTHER OPTIONS
!!$  WRITE(LUMOD2,'(A)')'    IF((Psegmented.AND.Qsegmented).AND.(nPrimQP.EQ.1))THEN'
!!$  WRITE(LUMOD2,'(A)')'     call SPICI_'//ARCSTRING//'_OBS_Seg1Prim2(nPrimA,nPrimB,nPrimC,nPrimD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
!!$  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
!!$  WRITE(LUMOD2,'(A)')'       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & PQorder,LOCALINTS_SP,localintsmaxsize,&'
!!$  WRITE(LUMOD2,'(A)')'       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&'
!!$  WRITE(LUMOD2,'(A)')'       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&'
!!$  IF(CPU)THEN
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
!!$  ELSE
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
!!$  ENDIF
!!$  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented.AND.Qsegmented)THEN'
!!$  WRITE(LUMOD2,'(A)')'     call SPICI_'//ARCSTRING//'_OBS_Seg2(nPrimA,nPrimB,nPrimC,nPrimD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
!!$  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
!!$  WRITE(LUMOD2,'(A)')'       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & PQorder,LOCALINTS_SP,localintsmaxsize,&'
!!$  WRITE(LUMOD2,'(A)')'       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&'
!!$  WRITE(LUMOD2,'(A)')'       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&'
!!$  IF(CPU)THEN
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
!!$  ELSE
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
!!$  ENDIF
!!$  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented)THEN'
!!$  WRITE(LUMOD2,'(A)')'     call SPICI_'//ARCSTRING//'_OBS_SegP2(nPrimA,nPrimB,nPrimC,nPrimD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
!!$  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
!!$  WRITE(LUMOD2,'(A)')'       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & PQorder,LOCALINTS_SP,localintsmaxsize,&'
!!$  WRITE(LUMOD2,'(A)')'       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&'
!!$  WRITE(LUMOD2,'(A)')'       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&'
!!$  IF(CPU)THEN
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
!!$  ELSE
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
!!$  ENDIF
!!$  WRITE(LUMOD2,'(A)')'    ELSEIF(Qsegmented)THEN'
!!$  WRITE(LUMOD2,'(A)')'     call SPICI_'//ARCSTRING//'_OBS_SegQ2(nPrimA,nPrimB,nPrimC,nPrimD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
!!$  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
!!$  WRITE(LUMOD2,'(A)')'       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & PQorder,LOCALINTS_SP,localintsmaxsize,&'
!!$  WRITE(LUMOD2,'(A)')'       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&'
!!$  WRITE(LUMOD2,'(A)')'       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&'
!!$  IF(CPU)THEN
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
!!$  ELSE
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
!!$  ENDIF
!!$
!!$  WRITE(LUMOD2,'(A)')'    ELSE'
!!$  WRITE(LUMOD2,'(A)')'     call SPICI_'//ARCSTRING//'_OBS_Gen2(nPrimA,nPrimB,nPrimC,nPrimD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
!!$  WRITE(LUMOD2,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp_SP,qexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & ACC_SP,BCC_SP,CCC_SP,DCC_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
!!$  WRITE(LUMOD2,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
!!$  WRITE(LUMOD2,'(A)')'       & pcent_SP,qcent_SP,Ppreexpfac_SP,Qpreexpfac_SP,nTABFJW1,nTABFJW2,TABFJW_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qiprim1,Qiprim2,Aexp_SP,Bexp_SP,Cexp_SP,Dexp_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & Qsegmented,Psegmented,reducedExponents_SP,integralPrefactor_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12_SP,Qdistance12_SP,&'
!!$  WRITE(LUMOD2,'(A)')'       & PQorder,LOCALINTS_SP,localintsmaxsize,&'
!!$  WRITE(LUMOD2,'(A)')'       & Acenter_SP,Bcenter_SP,Ccenter_SP,Dcenter_SP,nAtomsA,nAtomsB,spherical,&'
!!$  WRITE(LUMOD2,'(A)')'       & TmpArray1_SP,TMParray1maxsize,TmpArray2_SP,TMParray2maxsize,&'
!!$  IF(CPU)THEN
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass)' 
!!$  ELSE
!!$     WRITE(LUMOD2,'(A)')'       & IatomAPass,iatomBPass,iASync)' 
!!$  ENDIF
!!$
!!$  WRITE(LUMOD2,'(A)')'    ENDIF'
!!$  WRITE(LUMOD2,'(A)')'   ENDIF'
!!$  WRITE(LUMOD2,'(A)')'   !COPY TO DP'
!!$  WRITE(LUMOD2,'(A)')'   LOCALINTS  = REAL(LOCALINTS_SP,KIND=realk)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(Aexp_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(Bexp_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(Cexp_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(Dexp_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(pexp_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(qexp_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(pcent_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(qcent_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(QpreExpFac_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(PpreExpFac_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(TABFJW_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(ACC_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(BCC_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(CCC_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(DCC_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(LOCALINTS_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(integralPrefactor_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(reducedExponents_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(Pdistance12_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(Acenter_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(Bcenter_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(TmpArray1_SP)'
!!$  WRITE(LUMOD2,'(A)')'  deallocate(TmpArray2_SP)'
!!$  WRITE(LUMOD2,'(A)')'  ENDIF'
  WRITE(LUMOD2,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_general'
  WRITE(LUMOD2,'(A)')'  '

  DO IFILES = 1,2
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

     IF(IFILES.EQ.1)THEN
        IF(Gen)THEN
           WRITE(LUMOD3,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_Gen(nPrimA,nPrimB,nPrimC,nPrimD,&'
           ILUMOD = 3
        ELSEIF(SegQ)THEN
           WRITE(LUMOD4,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_SegQ(nPrimA,nPrimB,nPrimC,nPrimD,&'
           ILUMOD = 4
        ELSEIF(SegP)THEN
           WRITE(LUMOD5,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_SegP(nPrimA,nPrimB,nPrimC,nPrimD,&'
           ILUMOD = 5
        ELSEIF(Seg)THEN
           WRITE(LUMOD6,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_Seg(nPrimA,nPrimB,nPrimC,nPrimD,&'
           ILUMOD = 6
        ELSEIF(Seg1Prim)THEN
           WRITE(LUMOD7,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_Seg1Prim(nPrimA,nPrimB,nPrimC,nPrimD,&'
           ILUMOD = 7
        ENDIF
     ELSE
        IF(Gen)THEN
           WRITE(LUMOD8,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_Gen2(nPrimA,nPrimB,nPrimC,nPrimD,&'
           ILUMOD = 8
        ELSEIF(SegQ)THEN
           WRITE(LUMOD9,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_SegQ2(nPrimA,nPrimB,nPrimC,nPrimD,&'
           ILUMOD = 9
        ELSEIF(SegP)THEN
           WRITE(LUMOD10,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_SegP2(nPrimA,nPrimB,nPrimC,nPrimD,&'
           ILUMOD = 10
        ELSEIF(Seg)THEN
           WRITE(LUMOD11,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_Seg2(nPrimA,nPrimB,nPrimC,nPrimD,&'
           ILUMOD = 11
        ELSEIF(Seg1Prim)THEN
           WRITE(LUMOD12,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_Seg1Prim2(nPrimA,nPrimB,nPrimC,nPrimD,&'
           ILUMOD = 12
        ENDIF
     ENDIF
     WRITE(ILUMOD,'(A)')'       & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
     WRITE(ILUMOD,'(A)')'       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
     WRITE(ILUMOD,'(A)')'       & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
     WRITE(ILUMOD,'(A)')'       & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
     WRITE(ILUMOD,'(A)')'       & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
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
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD'
     WRITE(ILUMOD,'(A)')'    integer,intent(in) :: nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV'
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
     WRITE(ILUMOD,'(A)')'    IF(UseGeneralCode) AngmomID = AngmomID + 10000 !force to use general code'
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
     DO AngmomA = 0,MaxAngmomSpecial
        DO AngmomB = 0,AngmomA
           DO AngmomC = 0,AngmomA
              DO AngmomD = 0,AngmomC

               BUILD(AngmomA,AngmomB,AngmomC,AngmomD) = .FALSE.
               IF(IFILES.EQ.1)THEN

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
                    STRINGIN(1:12)  = 'TMParray1(1)'
                    STRINGOUT(1:12) = 'TMParray2(1)'                                      
                    TMPSTRING(1:12) = '            '
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
                    STRINGIN(1:12)  = 'TMParray1(1)'
                    STRINGOUT(1:12) = 'TMParray2(1)'
                    TMPSTRING(1:12) = '            '
                    call subroutineMAIN(ILUMOD,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
                         & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,nlmA,nlmB,nlmC,nlmD,spherical,Gen,SegQ,SegP,Seg,Seg1Prim)
                 ENDIF

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

                  IF(IFILES.EQ.2)THEN
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
                       STRINGIN(1:12)  = 'TMParray1(1)'
                       STRINGOUT(1:12) = 'TMParray2(1)'
                       TMPSTRING(1:12) = '            '
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
                       STRINGIN(1:12)  = 'TMParray1(1)'
                       STRINGOUT(1:12) = 'TMParray2(1)'
                       TMPSTRING(1:12) = '            '
                       call subroutineMAIN(ILUMOD,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN,STRINGOUT,TMPSTRING,nTUV,AngmomP,AngmomQ,&
                            & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,nlmA,nlmB,nlmC,nlmD,spherical,Gen,SegQ,SegP,Seg,Seg1Prim)
                    ENDIF

                    !==========================00=====================================
                  ENDIF !IFILES
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     WRITE(ILUMOD,'(A)')'    CASE DEFAULT'
     IF(IFILES.EQ.2)THEN
        WRITE(ILUMOD,'(A)')'#ifdef VAR_OPENACC'
        WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''ICI_CPU_McM_general called with OpenACC'',-1)'
        WRITE(ILUMOD,'(A)')'#endif'
        WRITE(ILUMOD,'(A)')'        call ICI_CPU_McM_general(nPrimA,nPrimB,nPrimC,nPrimD,&'
        WRITE(ILUMOD,'(A)')'           & nPrimP,nPrimQ,nPrimQP,nPasses,MaxPasses,IntPrint,lupri,&'
        WRITE(ILUMOD,'(A)')'           & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&'
        WRITE(ILUMOD,'(A)')'           & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&'
        WRITE(ILUMOD,'(A)')'           & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&'
        WRITE(ILUMOD,'(A)')'           & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&'
        WRITE(ILUMOD,'(A)')'           & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&'
        WRITE(ILUMOD,'(A)')'           & Qiprim1,Qiprim2,Aexp,Bexp,Cexp,Dexp,&'
        WRITE(ILUMOD,'(A)')'           & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&'
        WRITE(ILUMOD,'(A)')'           & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,LOCALINTS,localintsmaxsize,&'
        WRITE(ILUMOD,'(A)')'           & Acenter,Bcenter,Ccenter,Dcenter,nAtomsA,nAtomsB,spherical,&'
        WRITE(ILUMOD,'(A)')'           & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&'
        WRITE(ILUMOD,'(A)')'           & IatomAPass,iatomBPass)' 
     ELSE
        IF(Gen)THEN
           WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in ICI_'//ARCSTRING//'_OBS_Gen'',-1)'
        ELSEIF(SegQ)THEN
           WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in ICI_'//ARCSTRING//'_OBS_SegQ'',-1)'
        ELSEIF(SegP)THEN
           WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in ICI_'//ARCSTRING//'_OBS_SegP'',-1)'
        ELSEIF(Seg)THEN
           WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in ICI_'//ARCSTRING//'_OBS_Seg'',-1)'
        ELSEIF(Seg1Prim)THEN
           WRITE(ILUMOD,'(A)')'        CALL ICHORQUIT(''Unknown Case in ICI_'//ARCSTRING//'_OBS_Seg1Prim'',-1)'
        ENDIF
     ENDIF
     WRITE(ILUMOD,'(A)')'    END SELECT'
     IF(IFILES.EQ.1)THEN
        IF(Gen)THEN
           WRITE(ILUMOD,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_Gen'
        ELSEIF(SegQ)THEN
           WRITE(ILUMOD,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_SegQ'
        ELSEIF(SegP)THEN
           WRITE(ILUMOD,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_SegP'
        ELSEIF(Seg)THEN
           WRITE(ILUMOD,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_Seg'
        ELSEIF(Seg1Prim)THEN
           WRITE(ILUMOD,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_Seg1Prim'
        ENDIF
     ELSE
        IF(Gen)THEN
           WRITE(ILUMOD,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_Gen2'
        ELSEIF(SegQ)THEN
           WRITE(ILUMOD,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_SegQ2'
        ELSEIF(SegP)THEN
           WRITE(ILUMOD,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_SegP2'
        ELSEIF(Seg)THEN
           WRITE(ILUMOD,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_Seg2'
        ELSEIF(Seg1Prim)THEN
           WRITE(ILUMOD,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_Seg1Prim2'
        ENDIF
     ENDIF
     WRITE(ILUMOD,'(A)')'  '
  ENDDO
  ENDDO !FILES

  WRITE(LUMOD2,'(A)')'  '
  WRITE(LUMOD2,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_general_size(TMParray1maxsize,&'
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
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_general_sizeSeg1Prim(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD2,'(A)')'         & nContA,nContB,nContC,nContD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)'
  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented.AND.Qsegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_general_sizeSeg(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD2,'(A)')'         & nContA,nContB,nContC,nContD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)'
  WRITE(LUMOD2,'(A)')'    ELSEIF(Psegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_general_sizeSegP(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD2,'(A)')'         & nContA,nContB,nContC,nContD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)'
  WRITE(LUMOD2,'(A)')'    ELSEIF(Qsegmented)THEN'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_general_sizeSegQ(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD2,'(A)')'         & nContA,nContB,nContC,nContD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)'
  WRITE(LUMOD2,'(A)')'    ELSE'
  WRITE(LUMOD2,'(A)')'     call ICI_'//ARCSTRING//'_OBS_general_sizeGen(TMParray1maxsize,&'
  WRITE(LUMOD2,'(A)')'         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&'
  WRITE(LUMOD2,'(A)')'         & nContA,nContB,nContC,nContD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimA,nPrimB,nPrimC,nPrimD,&'
  WRITE(LUMOD2,'(A)')'         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)'
  WRITE(LUMOD2,'(A)')'    ENDIF'
  WRITE(LUMOD2,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_general_size'
  WRITE(LUMOD2,'(A)')'  '

  WRITE(LUMOD2,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralMod'
  WRITE(LUMOD3,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGen'
  WRITE(LUMOD4,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQ'
  WRITE(LUMOD5,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegP'
  WRITE(LUMOD6,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg'
  WRITE(LUMOD7,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1Prim'
  WRITE(LUMOD8,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGen2'
  WRITE(LUMOD9,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQ2'
  WRITE(LUMOD10,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegP2'
  WRITE(LUMOD11,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg2'
  WRITE(LUMOD12,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1Prim2'

  close(unit = LUMOD2)
  close(unit = LUMOD3)
  close(unit = LUMOD4)
  close(unit = LUMOD5)
  close(unit = LUMOD6)
  close(unit = LUMOD7)
  close(unit = LUMOD8)
  close(unit = LUMOD9)
  close(unit = LUMOD10)
  close(unit = LUMOD11)
  close(unit = LUMOD12)

  LUMOD3=3
  open(unit = LUMOD3, file='MAIN_'//ARCSTRING//'_OBS_DRIVERGenSize.f90',status="unknown")
  WRITE(LUMOD3,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGenSize'
  WRITE(LUMOD3,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD3,'(A)')'!Contains routines for General Contracted Basisset '

  LUMOD4=4
  open(unit = LUMOD4, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSegQSize.f90',status="unknown")
  WRITE(LUMOD4,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQSize'
  WRITE(LUMOD4,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD4,'(A)')'!Contains routines for a General Contracted LHS Segmented contracted RHS and Basisset '

  LUMOD5=5
  open(unit = LUMOD5, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSegPSize.f90',status="unknown")
  WRITE(LUMOD5,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegPSize'
  WRITE(LUMOD5,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD5,'(A)')'!Contains routines for Segmented contracted LHS and a General Contracted RHS Basisset '

  LUMOD6=6
  open(unit = LUMOD6, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSegSize.f90',status="unknown")
  WRITE(LUMOD6,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegSize'
  WRITE(LUMOD6,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD6,'(A)')'!Contains routines for Segmented contracted Basisset '

  LUMOD7=7
  open(unit = LUMOD7, file='MAIN_'//ARCSTRING//'_OBS_DRIVERSeg1PrimSize.f90',status="unknown")
  WRITE(LUMOD7,'(A)')'MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1PrimSize'
  WRITE(LUMOD7,'(A)')'!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory'
  WRITE(LUMOD7,'(A)')'!Contains routines for Segmented contracted Basisset containing a single primitive'
  
  DO ILUMOD=3,7
!     WRITE(ILUMOD,'(A)')'use IchorprecisionMod'
     WRITE(ILUMOD,'(A)')'use IchorCommonMod'
     WRITE(ILUMOD,'(A)')'use IchorEriCoulombintegralCPUMcMGeneralMod'
     WRITE(ILUMOD,'(A)')'  '
     WRITE(ILUMOD,'(A)')'private   '
  ENDDO
  WRITE(LUMOD2,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_general_size  '
  WRITE(LUMOD3,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_general_sizeGen  '
  WRITE(LUMOD4,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_general_sizeSegQ  '
  WRITE(LUMOD5,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_general_sizeSegP  '
  WRITE(LUMOD6,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_general_sizeSeg  '
  WRITE(LUMOD7,'(A)')'public :: ICI_'//ARCSTRING//'_OBS_general_sizeSeg1Prim  '
  DO ILUMOD=2,7
     WRITE(ILUMOD,'(A)')'  '
     WRITE(ILUMOD,'(A)')'CONTAINS'
  ENDDO

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
        WRITE(LUMOD3,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_general_sizeGen(TMParray1maxsize,&'
        ILUMOD = LUMOD3
     ELSEIF(SegQ)THEN
        WRITE(LUMOD4,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_general_sizeSegQ(TMParray1maxsize,&'
        ILUMOD = LUMOD4
     ELSEIF(SegP)THEN
        WRITE(LUMOD5,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_general_sizeSegP(TMParray1maxsize,&'
        ILUMOD = LUMOD5
     ELSEIF(Seg)THEN
        WRITE(LUMOD6,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_general_sizeSeg(TMParray1maxsize,&'
        ILUMOD = LUMOD6
     ELSEIF(Seg1Prim)THEN
        WRITE(LUMOD7,'(A)')'  subroutine ICI_'//ARCSTRING//'_OBS_general_sizeSeg1Prim(TMParray1maxsize,&'
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
     WRITE(ILUMOD,'(A)')'    IF(UseGeneralCode) AngmomID = AngmomID + 10000 !force to use general code'
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
                 STRINGIN2(1:9)  = 'TMParray1'
                 STRINGOUT2(1:9) = 'TMParray2'
                 TMPSTRING2(1:9) = '         '
                 call determineSizes(ILUMOD,AngmomA,AngmomB,AngmomC,AngmomD,STRINGIN2,STRINGOUT2,TMPSTRING2,nTUV,AngmomP,AngmomQ,&
                      & AngmomPQ,nTUVP,nTUVQ,nTUVAspec,nTUVBspec,nTUVCspec,nTUVDspec,nlmA,nlmB,nlmC,nlmD,spherical,Gen,SegQ,SegP,Seg,Seg1Prim)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     WRITE(ILUMOD,'(A)')'    CASE DEFAULT'
     WRITE(ILUMOD,'(A)')'     call ICI_CPU_McM_general_size(TMParray1maxsize,&'
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
     WRITE(ILUMOD,'(A)')'    END SELECT'

     IF(Gen)THEN
        WRITE(LUMOD3,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_general_sizeGen'
        ILUMOD = LUMOD3
     ELSEIF(SegQ)THEN
        WRITE(LUMOD4,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_general_sizeSegQ'
        ILUMOD = LUMOD4
     ELSEIF(SegP)THEN
        WRITE(LUMOD5,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_general_sizeSegP'
        ILUMOD = LUMOD5
     ELSEIF(Seg)THEN
        WRITE(LUMOD6,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_general_sizeSeg'
        ILUMOD = LUMOD6
     ELSEIF(Seg1Prim)THEN
        WRITE(LUMOD7,'(A)')'  end subroutine ICI_'//ARCSTRING//'_OBS_general_sizeSeg1Prim'
        ILUMOD = LUMOD7
     ENDIF
  ENDDO

  WRITE(LUMOD3,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModGenSize'
  WRITE(LUMOD4,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegQSize'
  WRITE(LUMOD5,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegPSize'
  WRITE(LUMOD6,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSegSize'
  WRITE(LUMOD7,'(A)')'END MODULE IchorEriCoulombintegral'//ARCSTRING//'OBSGeneralModSeg1PrimSize'

!  close(unit = LUMOD2)
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
    character(len=12) :: STRINGIN,STRINGOUT,TMPSTRING
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
!!$    IF(STRING(1:4).NE.'LOCALINTS')THEN
!!$       WRITE(LUPRI2,'(A)')'#ifdef VAR_DEBUGICHOR'
!!$       call initString(8)
!!$       call AddToString('IF(')
!!$       call AddToString(DIMSTRING)
!!$       call AddToString('*')
!!$       call AddToString(DIMINT)
!!$       call AddToString('.GT.')                
!!$       call AddToString(STRING)
!!$       call AddToString('maxsize)THEN')
!!$       call writeString(LUPRI2)
!!$       
!!$       call initString(10)
!!$       call AddToString('call ichorquit(''')
!!$       call AddToString(DIMSTRING)
!!$       call AddToString(' too small'',-1)')
!!$       call writeString(LUPRI2)
!!$       WRITE(LUPRI2,'(A)')'        ENDIF'
!!$       WRITE(LUPRI2,'(A)')'#endif'
!!$    ENDIF
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

