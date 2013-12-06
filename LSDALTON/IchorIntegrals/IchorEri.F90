!> @file
!> Contains the main Ichor integral drivers for calculation electron repulsion integrals
!> Ichor is the "Integral Code Hand Optimized for Rapid evaluation" 
!> Ichor is the ethereal golden fluid that is the blood of the greek gods

!> \brief Main Ichor drivers for the calculation of integrals 
!> based on the Obara Saika(OS)/Head-Gordon-Pople(HGP)
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorErimodule
  use IchorprecisionModule
  use IchorCommonModule
  use IchorBatchToolsModule
!  use IchorEriCoulombintegralMod, only: IchorCoulombIntegral_SSSS
!  use IchorEriCoulombintegralGeneralMod, only: IchorCoulombIntegral_McM_general
  use IchorEriCoulombintegralOBSGeneralMod, only: IchorCoulombIntegral_OBS_general, &
       & IchorCoulombIntegral_OBS_general_size
  use IchorCoulombIntegral_seg_seg_SSSS_mod, only: IchorCoulombIntegral_seg_seg_SSSS
  use IchorMemory
  use IchorGammaTabulationModule
  use IchorParametersModule
  use IchorSaveGabModule

public:: IchorEri
private
CONTAINS
subroutine IchorEri(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
     & AngmomOfTypeA,nAtomsOfTypeA,nPrimOfTypeA,nContOfTypeA,&
     & startOrbitalOfTypeA,Acenters,exponentsOfTypeA,ContractCoeffOfTypeA,&
     & startBatchA,endBatchA,&
     & nTypesB,MaxNatomsB,MaxnPrimB,MaxnContB,&
     & AngmomOfTypeB,nAtomsOfTypeB,nPrimOfTypeB,nContOfTypeB,&
     & startOrbitalOfTypeB,Bcenters,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & startBatchB,endBatchB,&
     & nTypesC,MaxNatomsC,MaxnPrimC,MaxnContC,&
     & AngmomOfTypeC,nAtomsOfTypeC,nPrimOfTypeC,nContOfTypeC,&
     & startOrbitalOfTypeC,Ccenters,exponentsOfTypeC,ContractCoeffOfTypeC,&
     & startBatchC,endBatchC,&
     & nTypesD,MaxNatomsD,MaxnPrimD,MaxnContD,&
     & AngmomOfTypeD,nAtomsOfTypeD,nPrimOfTypeD,nContOfTypeD,&
     & startOrbitalOfTypeD,Dcenters,exponentsOfTypeD,ContractCoeffOfTypeD,&
     & startBatchD,endBatchD,&
     & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorInputDim1,&
     & IchorInputDim2,IchorInputDim3,&
     & InputStorage,IchorParSpec,IchorScreenSpec,THRESHOLD_OD,&
     & THRESHOLD_CS,THRESHOLD_QQR,&
     & IchorGabID1,IchorGabID2,IchorDebugSpec,&
     & IchorAlgoSpec,IchorPermuteSpec,filestorageIdentifier,MaxMem,&
     & MaxFileStorage,MaxMemAllocated,MemAllocated,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
     & OutputStorage,lupri)
implicit none
!> nTypesA is the number of different types of shells, each type is defined by 
!> an angular momentum, a number of primitives(nPrim), a number of contracted functions
!> (nCont) a set of exponents and a set of contraction coefficients. [For Center A]
integer,intent(in) :: nTypesA
!> MaxnAtomsA is the maximum number of Atoms that have the same type. [For Center A]
integer,intent(in) :: MaxnAtomsA
!> MaxnPrim is the maximum number of Primitives among all the types. [For Center A]
integer,intent(in) :: MaxnPrimA
!> MaxnPrim is the maximum number of Contracted functions among all the types. [For Center A]
integer,intent(in) :: MaxnContA
!> AngmomOfTypeA is the angular momentum for each type. [For Center A]
Integer,intent(in) :: AngmomOfTypeA(ntypesA)
!> nAtomsOfTypeA is the number of atoms that have the given type. [For Center A]
Integer,intent(in) :: nAtomsOfTypeA(ntypesA)
!> nPrimOfTypeA is the number of primitive for the given type. [For Center A]
Integer,intent(in) :: nPrimOfTypeA(ntypesA)
!> nContOfTypeA is the number of contracted function for the given type. [For Center A]
Integer,intent(in) :: nContOfTypeA(ntypesA)
!> startorbital for this atom of this type
Integer,intent(in) :: startOrbitalOfTypeA(MaxNatomsA,ntypesA)
!> Acenters is the centers of all the atoms that have the given type. [For Center A]
Real(realk),intent(in) :: Acenters(3,MaxNatomsA,ntypesA)
!> exponentsOfTypeA is the nPrim exponents for the given type. [For Center A]
Real(realk),intent(in) :: exponentsOfTypeA(MaxnprimA,ntypesA)
!> ContractCoeffOfTypeA is the contrction coefficient matrix (nPrim,nCont) for the given type. [For Center A]
Real(realk),intent(in) :: ContractCoeffOfTypeA(MaxnprimA,MaxnContA,ntypesA)
!> startindex of first Batch
Integer,intent(in) :: startBatchA
!> endindex of last Batch
Integer,intent(in) :: endBatchA
!
! Same for Center B
!
integer,intent(in) :: nTypesB,MaxnAtomsB,MaxnPrimB,MaxnContB,startBatchB,endBatchB
Integer,intent(in) :: AngmomOfTypeB(ntypesB),nAtomsOfTypeB(ntypesB)
Integer,intent(in) :: nContOfTypeB(ntypesB),nPrimOfTypeB(ntypesB),startOrbitalOfTypeB(MaxNatomsB,ntypesB)
Real(realk),intent(in) :: Bcenters(3,MaxNatomsB,ntypesB),exponentsOfTypeB(MaxnprimB,ntypesB)
Real(realk),intent(in) :: ContractCoeffOfTypeB(MaxnprimB,MaxnContB,ntypesB)
!
! Same for Center C
!
integer,intent(in) :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC,startBatchC,endBatchC
Integer,intent(in) :: AngmomOfTypeC(ntypesC),nAtomsOfTypeC(ntypesC)
Integer,intent(in) :: nContOfTypeC(ntypesC),nPrimOfTypeC(ntypesC),startOrbitalOfTypeC(MaxNatomsC,ntypesC)
Real(realk),intent(in) :: Ccenters(3,MaxNatomsC,ntypesC),exponentsOfTypeC(MaxnprimC,ntypesC)
Real(realk),intent(in) :: ContractCoeffOfTypeC(MaxnprimC,MaxnContC,ntypesC)
!
! Same for Center D
!
integer,intent(in) :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD,startBatchD,endBatchD
Integer,intent(in) :: AngmomOfTypeD(ntypesD),nAtomsOfTypeD(ntypesD)
Integer,intent(in) :: nContOfTypeD(ntypesD),nPrimOfTypeD(ntypesD),startOrbitalOfTypeD(MaxNatomsD,ntypesD)
Real(realk),intent(in) :: Dcenters(3,MaxNatomsD,ntypesD),exponentsOfTypeD(MaxnprimD,ntypesD)
Real(realk),intent(in) :: ContractCoeffOfTypeD(MaxnprimD,MaxnContD,ntypesD)
!
!> Spherical Specification (SphericalSpec = SphericalParam = 1) means to use Spherical Harmonic basis functions
Integer,intent(in) :: SphericalSpec
!> Job Specification (IcorJob = IcorJobEri = 1) means that the 4 center 2 electron repulsion integrals
!> should be calculated. 
Integer,intent(in) :: IchorJobSpec
!> Input Specification (IchorInputSpec = IcorInputNoInput = 1) means no Input have been provided
Integer,intent(in) :: IchorInputSpec
!> Input dimensions assuming InputStorage(IchorInputDim1,IchorInputDim2,IchorInputDim3)
Integer,intent(in) :: IchorInputDim1,IchorInputDim2,IchorInputDim3
!> InputStorage
real(realk),intent(in) :: InputStorage(IchorInputDim1,IchorInputDim2,IchorInputDim3)
!> Parallelization specification (communicator and other stuff should be set up with other call) 
!> IchorParSpec = IchorParNone = 1 means no parallelization - no OpenMP, no MPI, no GPU
Integer,intent(in) :: IchorParSpec
!> Screening specification 
!> IchorScreenSpec = IchorScreen = 1 means default screening including Cauchy-Schwarz screening and QQR and OD
!> IchorScreenSpec = IchorScreenNone = 2 means no screening
Integer,intent(in) :: IchorScreenSpec
!> Overlap Density Screening 
real(realk),intent(in) :: THRESHOLD_OD
!> Cauchy-Schwarz screening Threshold only used if Cauchy-Schwarz screening is activated
real(realk),intent(in) :: THRESHOLD_CS
!> Cauchy-Schwarz screening with distance dependence (QQR) Threshold
real(realk),intent(in) :: THRESHOLD_QQR
!> Screening Matrix LHS Identification, used in connection with IchorSaveGabModule
Integer,intent(in) :: IchorGabID1
!> Screening Matrix RHS Identification, used in connection with IchorSaveGabModule
Integer,intent(in) :: IchorGabID2
!> Debug info specification - The print Level or IchorDebugNone=0
!> IchorDebugSpec = IchorDebugNone = 0 means no debug info
Integer,intent(in) :: IchorDebugSpec
!> Integral Algorithm specification 
!> IchorAlgoSpec = IchorAlgoOS = 1 means Obara-Saika (Head-Gordon Pople)
Integer,intent(in) :: IchorAlgoSpec
!> Permutation specification (SameLHSaos, SameRHSaos, SameODs) 
!> IchorPermuteSpec = IchorPermuteTTT = 1 means (SameLHSaos=.TRUE., SameRHSaos=.TRUE., SameODs=.TRUE.) 
Integer,intent(in) :: IchorPermuteSpec
!> Identifier to determine which file should be used to save integrals on disk, This 
!> should be logical unit number, if IchorNofilestorage=0 no file is open. 
Integer,intent(in) :: filestorageIdentifier
!> Maximum Memory that the integral program is allowed to use. Zero means no restrictions 
Integer(kind=long),intent(in) :: MaxMem
!> Maximum File size, if zero - no file will be written or read. 
Integer,intent(in) :: MaxFileStorage
!> Maximum Memory used in the program. The Ichor program adds to this value.  
Integer(kind=long),intent(inout) :: MaxMemAllocated
!> Memory allocated in the program. If the input value is not equal to the output value the
!> there is a memory leak inside the program. 
Integer(kind=long),intent(inout) :: MemAllocated
!> Output dimensions assuming 
!> OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5)
Integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5
!> OutputStorage
real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5)
!> Logical unit number of output file.
Integer,intent(in) :: lupri

! MISSING BATCHINFORMATION IF NOT THE FULL SET SHOULD BE CALCULATED

! Local variables
integer :: nPrimP,nContP,nPrimQ,nContQ
integer :: nTABFJW1,nTABFJW2,i1,i2,i3,i4,AngmomQ,TotalAngmom
integer :: i12,offset,K,I,iPrimQ,iPrimP,icont
integer :: oldmaxangmomABCD
integer :: ItypeA,ItypeB,itypeC,itypeD,AngmomA,AngmomB,AngmomC,AngmomD
integer :: ItypeAnon,ItypeBnon,itypeCnon,itypeDnon,nLocalInt
integer :: nPrimA,nPrimB,nContA,nAtomsA,nAtomsB,nAtomsC,nAtomsD,nOrbA
integer :: nDimA,nOrbCompA,nContB,nOrbB,nDimB
integer :: nPrimC,nContC,nPrimD,nContD,nOrbC,nOrbD,nOrbCompB,nOrbCompC,nDimC
integer :: nDimD,nOrbCompD,INTPRINT,maxangmomABCD
integer,pointer :: Piprim1(:),Piprim2(:),Qiprim1(:),Qiprim2(:)
integer,pointer :: OrderdListA(:),OrderdListB(:),OrderdListC(:),OrderdListD(:)
logical :: Psegmented,Qsegmented,PQorder,Spherical,SameLHSaos,SameRHSaos,SameODs
logical :: TriangularLHSAtomLoop,TriangularRHSAtomLoop,PermuteRHS,CSScreen
logical :: TriangularODAtomLoop
logical :: NOTDoSSSS,Segmented,PermuteLHSTypes,PermuteRHSTypes,PermuteODTypes
real(realk),pointer :: expP(:),Pcent(:),PpreExpFac(:),Pdistance12Pass(:,:),LocalInt(:)
real(realk),pointer :: expQ(:),PcentPass(:,:),PpreExpFacPass(:,:),inversexpP(:),LocalIntPass(:,:)
real(realk),pointer :: QpreExpFac(:),Qcent(:)
REAL(realk),pointer :: TABFJW(:,:),reducedExponents(:),integralPrefactor(:)
real(realk) :: AcenterSpec(3),BcenterSpec(3),CcenterSpec(3),DcenterSpec(3)
real(realk) :: Pdistance12(3),Qdistance12(3)
real(realk),pointer :: expA(:),ContractCoeffA(:,:),Acenter(:,:)
real(realk),pointer :: expB(:),ContractCoeffB(:,:),Bcenter(:,:)
real(realk),pointer :: expC(:),ContractCoeffC(:,:),Ccenter(:,:)
real(realk),pointer :: expD(:),ContractCoeffD(:,:),Dcenter(:,:)
integer,pointer :: startOrbitalA(:),startOrbitalB(:)
integer,pointer :: startOrbitalC(:),startOrbitalD(:)
!Tmporary array used in the innermost routines 
integer :: TMParray1maxsize,TMParray2maxsize
real(realk),pointer :: TmpArray1(:)
real(realk),pointer :: TmpArray2(:)
!Batch info
integer :: iBatchC,iBatchD,nBatchA,nBatchB,nBatchC,nBatchD
integer,pointer :: BatchIndexOfTypeA(:)
integer,pointer :: BatchIndexOfTypeB(:)
integer,pointer :: BatchIndexOfTypeC(:)
integer,pointer :: BatchIndexOfTypeD(:)
integer :: iBatchIndexOfTypeA,iBatchIndexOfTypeB,iBatchIndexOfTypeC,iBatchIndexOfTypeD
!Screening variables
real(realk),pointer :: MaxGabForTypeAB(:,:),MaxGabForTypeCD(:,:)
real(realk),pointer :: BATCHGAB(:,:),BATCHGCD(:,:),BATCHGAB2(:,:),BATCHGCD2(:,:)
logical,pointer :: noScreenAB(:,:),noScreenCD(:,:),noScreenCD2(:,:)
integer :: nBatchAGAB,nBatchBGAB,nBatchCGCD,nBatchDGCD
real(realk) :: GABELM
!Passes
integer :: nPasses,iPass,MaxPasses!,nPassLoop,nPassesInPassLoop,iPassTMP !,ndimPass
!integer :: iPassLoop
integer,pointer :: IatomAPass(:),IatomBPass(:)
integer :: MaxTotalAngmom,IAngmomTypes
!ODscreening
real(realk) :: ExtentOfTypeA(ntypesA),ExtentOfTypeB(ntypesB),extentA,extentB,extentC,extentD
real(realk) :: ExtentOfTypeC(ntypesC),ExtentOfTypeD(ntypesD)
real(realk) :: sumExtent2AB,sumExtent2CD
logical :: ODscreen,QQRscreen
! POSSIBLE IMPROVEMENTS / TO DO LIST:
! FFFF
! FUSE LOOPS IN TRANSFER
! Replace HARDCODED LOOPUNROLLED WITH LOOP?
! MKL primitives VML functions - google "List of VML Functions" ?
! Get CPI rate down! below 0.5!
! FIX noscreen not pretty to have 2 noscreen!
! nPasses collect - do not collect? PROFILE


call determineScreening(IchorScreenSpec,CSscreen,ODscreen,QQRscreen)
call build_ichor_AOextent(MaxnAtomsA,MaxnprimA,MaxnContA,ntypesA,exponentsOfTypeA,ODscreen,&
     & nContOfTypeA,nPrimOfTypeA,ContractCoeffOfTypeA,AngmomOfTypeA,Threshold_OD,ExtentOfTypeA)
call build_ichor_AOextent(MaxnAtomsB,MaxnprimB,MaxnContB,ntypesB,exponentsOfTypeB,ODscreen,&
     & nContOfTypeB,nPrimOfTypeB,ContractCoeffOfTypeB,AngmomOfTypeB,Threshold_OD,ExtentOfTypeB)
call build_ichor_AOextent(MaxnAtomsC,MaxnprimC,MaxnContC,ntypesC,exponentsOfTypeC,ODscreen,&
     & nContOfTypeC,nPrimOfTypeC,ContractCoeffOfTypeC,AngmomOfTypeC,Threshold_OD,ExtentOfTypeC)
call build_ichor_AOextent(MaxnAtomsD,MaxnprimD,MaxnContD,ntypesD,exponentsOfTypeD,ODscreen,&
     & nContOfTypeD,nPrimOfTypeD,ContractCoeffOfTypeD,AngmomOfTypeD,Threshold_OD,ExtentOfTypeD)
call set_ichor_memvar(MaxMemAllocated,MemAllocated)
INTPRINT=IchorDebugSpec
IF(INTPRINT.GT.50)WRITE(LUPRI,'(A)') ' IchorEri'
call mem_ichor_alloc(OrderdListA,nTypesA)
call GenerateOrderdListOfTypes(lupri,nTypesA,AngmomOfTypeA,OrderdListA)
call mem_ichor_alloc(OrderdListB,nTypesB)
call GenerateOrderdListOfTypes(lupri,nTypesB,AngmomOfTypeB,OrderdListB)
call mem_ichor_alloc(OrderdListC,nTypesC)
call GenerateOrderdListOfTypes(lupri,nTypesC,AngmomOfTypeC,OrderdListC)
call mem_ichor_alloc(OrderdListD,nTypesD)
call GenerateOrderdListOfTypes(lupri,nTypesD,AngmomOfTypeD,OrderdListD)
PQorder=.FALSE.
IF(CSScreen)THEN
   call mem_ichor_alloc(BatchIndexOfTypeA,nTypesA)
   call ConstructBatchIndexOfType(BatchIndexOfTypeA,nTypesA,nAtomsOfTypeA,nBatchA)
   call mem_ichor_alloc(BatchIndexOfTypeB,nTypesB)
   call ConstructBatchIndexOfType(BatchIndexOfTypeB,nTypesB,nAtomsOfTypeB,nBatchB)
   call RetrieveGabDimFromIchorSaveGabModule(nBatchAGAB,nBatchBGAB,IchorGabID1)
   call mem_ichor_alloc(BATCHGAB,nBatchA,nBatchB)
   IF(nBatchA.NE.endBatchA-startBatchA+1)call ichorQuit('Screening Dim Mismatch1',-1)
   IF(nBatchB.NE.endBatchB-startBatchB+1)call ichorQuit('Screening Dim Mismatch2',-1)
   IF(nBatchAGAB.NE.nBatchA.OR.nBatchBGAB.NE.nBatchB)THEN
      !WARNING BATCHCALC NOT FULL INTEGRAL IS CALCULATED
      call mem_ichor_alloc(BATCHGAB2,nBatchAGAB,nBatchBGAB)
      call RetrieveGabFromIchorSaveGabModule(nBatchAGAB,nBatchBGAB,IchorGabID1,BATCHGAB2)
      call ExtractBatchGabFromFullGab(nBatchA,nBatchB,BATCHGAB,nBatchAGAB,nBatchBGAB,BATCHGAB2,&
           & startBatchA,endBatchA,startBatchB,endBatchB)
      call mem_ichor_dealloc(BATCHGAB2)
   ELSE
      call RetrieveGabFromIchorSaveGabModule(nBatchA,nBatchB,IchorGabID1,BATCHGAB)
   ENDIF

   call mem_ichor_alloc(MaxGabForTypeAB,nTypesB,nTypesA)
   call ObtainMaxGabForType(MaxGabForTypeAB,nTypesA,nTypesB,nAtomsOfTypeA,&
        & nAtomsOfTypeB,BatchIndexOfTypeA,BatchIndexOfTypeB,BATCHGAB,&
        & nBatchA,nBatchB)
!   print*,'IchorEri: MaxGabForTypeBA:',MaxGabForTypeBA
   IF(IchorGabID1.EQ.IchorGabID2)THEN
      BATCHGCD => BATCHGAB
      MaxGabForTypeCD => MaxGabForTypeAB
      BatchIndexOfTypeC => BatchIndexOfTypeA
      BatchIndexOfTypeD => BatchIndexOfTypeB
      nBatchC = nBatchA
      nBatchD = nBatchB
   ELSE
      call mem_ichor_alloc(BatchIndexOfTypeC,nTypesC)
      call ConstructBatchIndexOfType(BatchIndexOfTypeC,nTypesC,nAtomsOfTypeC,nBatchC)
      call mem_ichor_alloc(BatchIndexOfTypeD,nTypesD)
      call ConstructBatchIndexOfType(BatchIndexOfTypeD,nTypesD,nAtomsOfTypeD,nBatchD)
      call mem_ichor_alloc(BATCHGCD,nBatchC,nBatchD)
      IF(nBatchC.NE.endBatchC-startBatchC+1)call ichorQuit('Screening Dim Mismatch3',-1)
      IF(nBatchD.NE.endBatchD-startBatchD+1)call ichorQuit('Screening Dim Mismatch4',-1)
      call RetrieveGabDimFromIchorSaveGabModule(nBatchCGCD,nBatchDGCD,IchorGabID2)      
      IF(nBatchCGCD.NE.nBatchC.OR.nBatchDGCD.NE.nBatchD)THEN
         !WARNING BATCHCALC NOT FULL INTEGRAL IS CALCULATED
         call mem_ichor_alloc(BATCHGCD2,nBatchCGCD,nBatchDGCD)
         call RetrieveGabFromIchorSaveGabModule(nBatchCGCD,nBatchDGCD,IchorGabID1,BATCHGCD2)
         call ExtractBatchGabFromFullGab(nBatchC,nBatchD,BATCHGCD,nBatchCGCD,nBatchDGCD,BATCHGCD2,&
              & startBatchC,endBatchC,startBatchD,endBatchD)
         call mem_ichor_dealloc(BATCHGCD2)
      ELSE
         call RetrieveGabFromIchorSaveGabModule(nBatchC,nBatchD,IchorGabID2,BATCHGCD)
      ENDIF
      call mem_ichor_alloc(MaxGabForTypeCD,nTypesC,nTypesD)
      call ObtainMaxGabForType(MaxGabForTypeCD,nTypesC,nTypesD,nAtomsOfTypeC,&
           & nAtomsOfTypeD,BatchIndexOfTypeC,BatchIndexOfTypeD,BATCHGCD,&
           & nBatchC,nBatchD)
   ENDIF
ELSE
   nBatchA = 1
   nBatchB = 1
   nBatchC = 1
   nBatchD = 1
   iBatchIndexOfTypeA = 0
   iBatchIndexOfTypeB = 0
   iBatchIndexOfTypeD = 0
   iBatchIndexOfTypeC = 0
   call mem_ichor_alloc(BATCHGAB,nBatchA,nBatchB)
   BATCHGCD => BATCHGAB
ENDIF
IF(CSScreen.OR.ODscreen)THEN
   !it is possible that we skip integrals
   call ichorzero2(OutputStorage,OutputDim1*OutputDim2,OutputDim3*OutputDim4*OutputDim5)
ENDIF

call determinePermuteSym(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
IF(INTPRINT.GT.50)write(lupri,*)'SameLHSaos,SameRHSaos,SameODs',SameLHSaos,SameRHSaos,SameODs
Spherical = SphericalSpec.EQ.SphericalParam
oldmaxangmomABCD = -25

MaxTotalAngmom = MAXVAL(AngmomOfTypeA) + MAXVAL(AngmomOfTypeB) &
     & + MAXVAL(AngmomOfTypeC) + MAXVAL(AngmomOfTypeD)

!we loop over Total angmom in order to ensure that we first do all
!SSSS integrals then PSSS,SPSS,SSPS,SSSP, ...
!this mean we call GAMMATABULATION a limited number of times and
!we reduce branch misprediction inside the code. 
DO IAngmomTypes = 0,MaxTotalAngmom
 DO ItypeDnon=1,nTypesD
  ! TYPE D CALC ===========================
  ItypeD = OrderdListD(ItypeDnon)  
  nAtomsD = nAtomsOfTypeD(ItypeD)
  IF(nAtomsD.EQ.0)CYCLE
  AngmomD = AngmomOfTypeD(ItypeD)
  nPrimD = nPrimOfTypeD(ItypeD)
  nContD = nContOfTypeD(ItypeD)
  extentD = ExtentOfTypeD(ItypeD)
  IF (spherical) THEN
   nOrbCompD = 2*(AngmomD+1)-1
  ELSE
   nOrbCompD = (AngmomD+1)*(AngmomD+2)/2
  ENDIF
  nOrbD = nContD*nOrbCompD
  nDimD = nContD*nOrbCompD*nAtomsD
  call mem_ichor_alloc(expD,nPrimD)
  call mem_ichor_alloc(ContractCoeffD,nPrimD,nContD)
  call mem_ichor_alloc(Dcenter,3,nAtomsD)
  call mem_ichor_alloc(StartOrbitalD,nAtomsD)
  call build_exp_ContractCoeff_center(nPrimD,nContD,nAtomsD,ntypesD,iTypeD,&
       & exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters,startOrbitalOfTypeD,&
       & expD,ContractCoeffD,Dcenter,StartOrbitalD,MaxnAtomsD,MaxnprimD,MaxnContD,lupri)
  ! DONE TYPE D CALC ===========================
  DO ItypeCnon=1,nTypesC
   ! TYPE C CALC ================================
   ItypeC = OrderdListC(ItypeCnon)  
   IF(SameRHSaos .AND. ItypeD.GT.ItypeC)CYCLE
   TriangularRHSAtomLoop = SameRHSaos .AND. ItypeD.EQ.ItypeC
   PermuteRHSTypes = SameRHSaos .AND. ItypeD.LT.ItypeC
   nAtomsC = nAtomsOfTypeC(ItypeC)
   IF(nAtomsC.EQ.0)CYCLE
   AngmomC = AngmomOfTypeC(ItypeC)
   nPrimC = nPrimOfTypeC(ItypeC)
   nContC = nContOfTypeC(ItypeC)
   extentC = ExtentOfTypeC(ItypeC)
   IF (spherical) THEN
      nOrbCompC = 2*(AngmomC+1)-1
   ELSE
      nOrbCompC = (AngmomC+1)*(AngmomC+2)/2
   ENDIF
   nOrbC = nContC*nOrbCompC
   nDimC = nOrbC*nAtomsC
   call mem_ichor_alloc(expC,nPrimC)
   call mem_ichor_alloc(ContractCoeffC,nPrimC,nContC)
   call mem_ichor_alloc(Ccenter,3,nAtomsC)
   call mem_ichor_alloc(StartOrbitalC,nAtomsC)
   call build_exp_ContractCoeff_center(nPrimC,nContC,nAtomsC,ntypesC,iTypeC,&
        & exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters,startOrbitalOfTypeC,&
        & expC,ContractCoeffC,Ccenter,StartOrbitalC,MaxnAtomsC,MaxnprimC,MaxnContC,lupri)
   ! DONE TYPE C CALC ===========================
   ! TYPE Q CALC ================================
   sumExtent2CD = (extentC + extentD)*(extentC + extentD)
   
   AngmomQ = AngmomC + AngmomD
   nPrimQ = nPrimC*nPrimD
   nContQ = nContC*nContD
   call mem_ichor_alloc(Qiprim1,nPrimQ) !not used
   call mem_ichor_alloc(Qiprim2,nPrimQ) !not used
   IF (nContQ.EQ. 1)THEN
      Qsegmented = .TRUE.
   ELSE
      Qsegmented = .FALSE.
   ENDIF
   call mem_ichor_alloc(expQ,nPrimQ)
   call build_expP(nPrimC,nPrimD,expC,expD,expQ)

   call mem_ichor_alloc(noScreenCD,nAtomsC,nAtomsD)
   call mem_ichor_alloc(noScreenCD2,nAtomsC,nAtomsD)
   IF(ODScreen)THEN
      call ODscreen_noScreen(nAtomsC,nAtomsD,Ccenter,Dcenter,sumextent2CD,noScreenCD)
!      call loutput(noscreenCD,nAtomsC,nAtomsD,6)
   ELSE
      call build_EmptynoScreen1(nAtomsC,nAtomsD,noScreenCD)
   ENDIF
   
   ! DONE TYPE Q CALC ===========================
   DO ItypeBnon=1,nTypesB
    ! TYPE B CALC ================================
    ItypeB = OrderdListB(ItypeBnon)
    nAtomsB = nAtomsOfTypeB(ItypeB)
    IF(nAtomsB.EQ.0)CYCLE
    AngmomB = AngmomOfTypeB(ItypeB)
    nPrimB = nPrimOfTypeB(ItypeB)
    nContB = nContOfTypeB(ItypeB)
    extentB = ExtentOfTypeB(ItypeB)
    IF (spherical) THEN
       nOrbCompB = 2*(AngmomB+1)-1
    ELSE
       nOrbCompB = (AngmomB+1)*(AngmomB+2)/2
    ENDIF
    nOrbB = nContB*nOrbCompB
    nDimB = nOrbB*nAtomsB
    call mem_ichor_alloc(expB,nPrimB)
    call mem_ichor_alloc(ContractCoeffB,nPrimB,nContB)
    call mem_ichor_alloc(Bcenter,3,nAtomsB)
    call mem_ichor_alloc(StartOrbitalB,nAtomsB)
    call build_exp_ContractCoeff_center(nPrimB,nContB,nAtomsB,ntypesB,iTypeB,&
         & exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters,startOrbitalOfTypeB,&
         & expB,ContractCoeffB,Bcenter,StartOrbitalB,MaxnAtomsB,MaxnprimB,MaxnContB,lupri)
    ! DONE TYPE B CALC ===========================
    DO ItypeAnon=1,nTypesA             !Ordered according to Angmom (highest angular momentum first)
     ! TYPE A CALC ================================
     ItypeA = OrderdListA(ItypeAnon)   !Due to Permutational Symmetry! (FD|PS) most eff. not (SP|DF)
     IF(SameLHSaos .AND. ItypeB.GT.ItypeA)CYCLE
     TriangularLHSAtomLoop = SameLHSaos .AND. ItypeB.EQ.ItypeA
     PermuteLHSTypes = SameLHSaos .AND. ItypeB.LT.ItypeA
     IF(SameODs .AND. ((ItypeC.GT.ItypeA).OR.((ItypeC.EQ.ItypeA).AND.(ItypeD.GT.ItypeB))))CYCLE
     PermuteODTypes = SameODs
!     IF(ItypeC.EQ.ItypeA.AND.ItypeD.EQ.ItypeB)PermuteODTypes=.FALSE.
     TriangularODAtomLoop = SameODs .AND. ((ItypeC.EQ.ItypeA).AND.(ItypeD.EQ.ItypeB))

     nAtomsA = nAtomsOfTypeA(ItypeA)
     IF(nAtomsA.EQ.0)CYCLE
     AngmomA = AngmomOfTypeA(ItypeA)
     TotalAngmom = AngmomA + AngmomB + AngmomQ
     IF(TotalAngmom.EQ.IAngmomTypes)THEN
      !This if statement ensures that we call GAMMATABULATION a very limited number of times
      !and it reduceses branch mispredictions inside the code,...
      nPrimA = nPrimOfTypeA(ItypeA)
      nContA = nContOfTypeA(ItypeA)
      extentA = ExtentOfTypeA(ItypeA)
      IF (spherical) THEN
         nOrbCompA = 2*(AngmomA+1)-1
      ELSE
         nOrbCompA = (AngmomA+1)*(AngmomA+2)/2
      ENDIF
      nOrbA = nContA*nOrbCompA
      nDimA = nContA*nOrbCompA*nAtomsA
      call mem_ichor_alloc(expA,nPrimA)
      call mem_ichor_alloc(ContractCoeffA,nPrimA,nContA)
      call mem_ichor_alloc(Acenter,3,nAtomsA)
      call mem_ichor_alloc(StartOrbitalA,nAtomsA)
      call build_exp_ContractCoeff_center(nPrimA,nContA,nAtomsA,ntypesA,iTypeA,&
           & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,startOrbitalOfTypeA,&
           & expA,ContractCoeffA,Acenter,StartOrbitalA,MaxnAtomsA,MaxnprimA,MaxnContA,lupri)
      ! DONE TYPE A CALC ===========================
      ! TYPE P CALC ================================
      sumExtent2AB = (extentA + extentB)*(extentA + extentB)
      nPrimP = nPrimA*nPrimB
      nContP = nContA*nContB
      call mem_ichor_alloc(Piprim1,nPrimP) !not used
      call mem_ichor_alloc(Piprim2,nPrimP) !not used
      IF (nContP.EQ. 1)THEN
         Psegmented = .TRUE.
      ELSE
         Psegmented = .FALSE.
      ENDIF
      call mem_ichor_alloc(expP,nPrimP)
      call mem_ichor_alloc(inversexpP,nPrimP)
      call build_expQ_inverseexpQ(nPrimA,nPrimB,expA,expB,expP,inversexpP)
      ! DONE TYPE P CALC ===========================
      
      ! TYPE PQ CALC ================================
      TotalAngmom = AngmomA + AngmomB + AngmomQ
      maxangmomABCD = AngmomA + AngmomB + AngmomQ
      IF(maxangmomABCD.NE.oldmaxangmomABCD)THEN
       IF(oldmaxangmomABCD.NE.-25)THEN
          call mem_ichor_dealloc(TABFJW)
       ENDIF
       nTABFJW1 = AngmomA + AngmomB + AngmomQ + 3 
       !only need + 3 after Branos change in BUILD_RJ000 
       nTABFJW2 = 1200
       !TABFJW(0:nTABFJW1,0:nTABFJW2)
       call mem_ichor_alloc(TABFJW,nTABFJW1,nTABFJW2,.TRUE.,.TRUE.)
       CALL GAMMATABULATION(lupri,maxangmomABCD,nTABFJW1,nTABFJW2,TABFJW)  
       oldmaxangmomABCD = maxangmomABCD
      ENDIF
      call mem_ichor_alloc(reducedExponents,nPrimP*nPrimQ)
      call mem_ichor_alloc(integralPrefactor,nPrimP*nPrimQ)
      call build_reducedExponents_integralPrefactorQP(nPrimP,nPrimQ,expQ,expP,&
           & reducedExponents,integralPrefactor)
      ! DONE TYPE PQ CALC ===========================
      
      NOTDoSSSS = .NOT.(TotalAngmom.EQ.0.AND.(Psegmented.AND.Qsegmented))
      IF(NOTDoSSSS)THEN
         !Determine Sizes of TmpArrays and MaxPasses
         call IchorCoulombIntegral_OBS_general_size(TMParray1maxsize,&
              &TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
              &nPrimP,nPrimQ,nContP,nContQ,nPrimQ*nPrimP,nContQ*nContP,&
              &Psegmented,Qsegmented)
         TMParray1maxsize = TMParray1maxsize
         TMParray2maxsize = TMParray2maxsize
         call mem_ichor_alloc(TmpArray1,TMParray1maxsize)
         call mem_ichor_alloc(TmpArray2,TMParray2maxsize)     
         CALL MEM_ICHOR_ALLOC(LocalInt,nOrbA*nOrbB*nOrbC*nOrbD)
         nLocalInt = nOrbA*nOrbB*nOrbC*nOrbD
         call mem_ichor_alloc(PpreExpFac,nPrimP)
         call mem_ichor_alloc(Pcent,3*nPrimP)
      ENDIF
      call mem_ichor_alloc(QpreExpFac,nPrimQ)
      call mem_ichor_alloc(Qcent,3*nPrimQ)      
      CALL MEM_ICHOR_ALLOC(LocalIntPass,nOrbA*nAtomsA*nOrbB*nAtomsB,nOrbC*nOrbD)
      !calc       
      IF (INTPRINT .GE. 10) THEN
         call PrintTypeExpInfo(nPrimP,nPrimQ,reducedExponents,integralPrefactor,lupri)
         call PrintTypeInfo(AngmomA,AngmomB,AngmomC,AngmomD,nPrimA,nPrimB,nPrimC,nPrimD,&
              & nContA,nContB,nContC,nContD,expA,ContractCoeffA,expB,ContractCoeffB,&
              & expC,ContractCoeffC,expD,ContractCoeffD,&
              & nAtomsA,nAtomsB,nAtomsC,nAtomsD,Acenter,Bcenter,Ccenter,Dcenter,lupri)
      ENDIF
      call mem_ichor_alloc(noScreenAB,nAtomsA,nAtomsB)
      IF(ODScreen)THEN
         call ODscreen_noScreen(nAtomsA,nAtomsB,Acenter,Bcenter,sumextent2AB,noScreenAB)
      ELSE
         call build_EmptynoScreen1(nAtomsA,nAtomsB,noScreenAB)
         GABELM = 0.0E0_realk
      ENDIF
      IF(CSScreen)THEN
         call build_noScreen1(ItypeC,ItypeD,ntypesC,ntypesD,nAtomsC,nAtomsD,nBatchC,&
              & nBatchD,BATCHGCD,MaxGabForTypeAB(ItypeA,ItypeB),THRESHOLD_CS,&
              & BatchIndexOfTypeC,BatchIndexOfTypeD,noScreenCD,noScreenCD2)
         iBatchIndexOfTypeA = BatchIndexOfTypeA(ItypeA)
         iBatchIndexOfTypeB = BatchIndexOfTypeB(ItypeB)
!         print*,'CS screenCD:'
!         call loutput(noscreenCD2,nAtomsC,nAtomsD,6)
      ELSE
         noscreenCD2 = noscreenCD
      ENDIF 
      
      !At this point we know all info about the integral except for the coordinates. 
      !Can be used for screening? Different screening for different types (SS|SS) 
      !APE type screening? 
      !For SSSS it should be exact to do
      !exp(-mu_AB R_AB^2)exp(-mu_CD R_CD^2) R000 and use
      !R000(alpha*R_pq^2) = F0(alpha*R_pq^2) = sqrt(pi/(4*x))erf(x) ;x=alpha*R_pq^2
      !R000(alpha*R_pq^2) > sqrt(pi/(4*x)) ;x=alpha*R_pq^2
      !look at possible things for SSSP type things.... 
      
      !LinK: Form active Dmat 1,2 or 4 (using screening) matrices (D_BD,D_AD,D_BC,D_AC)
      call mem_ichor_alloc(IatomAPass,nAtomsA*nAtomsB)
      call mem_ichor_alloc(IatomBPass,nAtomsA*nAtomsB)
      call mem_ichor_alloc(PpreExpFacPass,nPrimP,nAtomsA*nAtomsB)
      call mem_ichor_alloc(PcentPass,3*nPrimP,nAtomsA*nAtomsB)
      call mem_ichor_alloc(Pdistance12Pass,3,nAtomsA*nAtomsB)
      CALL Build_pcent_Pdistance12_PpreExpFac(nPrimA,nPrimB,natomsA,natomsB,nContA,nContB,&
           & inversexpP,expA,expB,Acenter,Bcenter,ContractCoeffA,ContractCoeffB,PSegmented,&
           & PcentPass,Pdistance12Pass,PpreExpFacPass,INTPRINT)
      IF(CSScreen)THEN
         iBatchIndexOfTypeD = BatchIndexOfTypeD(ItypeD)
         iBatchIndexOfTypeC = BatchIndexOfTypeC(ItypeC)
      ENDIF
      IF(NOTDoSSSS)THEN
         !     call IchorTimer('START',TSTART,TEND,LUPRI)
         
         !      WRITE(lupri,*)'CALC TYPE:',ItypeA,ItypeB,ItypeC,ItypeD
         !      WRITE(lupri,*)'CALC  Ang:',AngmomA,AngmomB,AngmomC,AngmomD
         call IchorTypeIntegralLoop(nAtomsA,nPrimA,nContA,nOrbCompA,startOrbitalA,&
              & iBatchIndexOfTypeA,expA,ContractCoeffA,AngmomA,Acenter,nBatchA,nOrbA,&
              & nAtomsB,nPrimB,nContB,nOrbCompB,startOrbitalB,&
              & iBatchIndexOfTypeB,expB,ContractCoeffB,AngmomB,Bcenter,nBatchB,nOrbB,&
              & nAtomsC,nPrimC,nContC,nOrbCompC,startOrbitalC,&
              & iBatchIndexOfTypeC,expC,ContractCoeffC,AngmomC,Ccenter,nBatchC,nOrbC,&
              & nAtomsD,nPrimD,nContD,nOrbCompD,startOrbitalD,&
              & iBatchIndexOfTypeD,expD,ContractCoeffD,AngmomD,Dcenter,nBatchD,nOrbD,&
              & nPrimP,nPrimQ,nContP,nContQ,expP,expQ,&
              & pcent,qcent,Ppreexpfac,Qpreexpfac,&
              & Qiprim1,Qiprim2,Piprim1,Piprim2,&
              & nTABFJW1,nTABFJW2,TABFJW,TotalAngmom,&
              & CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
              & Qsegmented,Psegmented,TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
              & reducedExponents,integralPrefactor,&
              & PcentPass,Pdistance12Pass,PpreExpFacPass,&
              & Pdistance12,Qdistance12,PQorder,&
              & BATCHGCD,IatomAPass,IatomBPass,BATCHGAB,&
              & LocalInt,nLocalInt,LocalIntPass,Spherical,TmpArray1,TMParray1maxsize,TmpArray2,&
              & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,OutputStorage,&
              & TMParray2maxsize,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri)
         call mem_ichor_dealloc(TmpArray1)
         call mem_ichor_dealloc(TmpArray2)
         CALL Mem_ichor_dealloc(LocalInt)
         call mem_ichor_dealloc(pcent)
         call mem_ichor_dealloc(PpreExpFac)
      ELSE
         call IchorsegsegSSSSIntegralLoop(nAtomsA,nPrimA,startOrbitalA,&
              & iBatchIndexOfTypeA,nBatchA,&
              & nAtomsB,nPrimB,startOrbitalB,iBatchIndexOfTypeB,nBatchB,&
              & nAtomsC,nPrimC,startOrbitalC,iBatchIndexOfTypeC,&
              & expC,ContractCoeffC,Ccenter,nBatchC,&
              & nAtomsD,nPrimD,startOrbitalD,iBatchIndexOfTypeD,&
              & expD,ContractCoeffD,Dcenter,nBatchD,&
              & nPrimP,nPrimQ,qcent,Qpreexpfac,&
              & TABFJW,CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
              & TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
              & reducedExponents,integralPrefactor,&
              & PcentPass,PpreExpFacPass,BATCHGCD,IatomAPass,IatomBPass,BATCHGAB,&
              & LocalIntPass,OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
              & OutputStorage,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri)
      ENDIF
      call mem_ichor_dealloc(Qcent)
      call mem_ichor_dealloc(QpreExpFac)
      call mem_ichor_dealloc(PpreExpFacPass)
      call mem_ichor_dealloc(PcentPass)
      call mem_ichor_dealloc(Pdistance12Pass)
      call mem_ichor_dealloc(IatomBPass) 
      call mem_ichor_dealloc(IatomAPass) 
      IF(PermuteODTypes)THEN
         !place (ABCD) in (CDAB) and possible (DCAB),(CDBA),(DCBA)
         call PermuteODtypesSub(OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputStorage,&
              & nAtomsA,nAtomsB,startOrbitalA,startOrbitalB,nOrbA,nOrbB,&
              & nAtomsC,nAtomsD,startOrbitalC,startOrbitalD,nOrbC,nOrbD,&
              & TriangularODAtomLoop,TriangularRHSAtomLoop,TriangularLHSAtomLoop,&
              & SameLHSaos,SameRHSaos,lupri)
      ENDIF

      call mem_ichor_dealloc(reducedExponents)
      call mem_ichor_dealloc(integralPrefactor)
      call mem_ichor_dealloc(noScreenAB)
      CALL MEM_ICHOR_deALLOC(LocalIntPass)
      !       call IchorTimer('TYPE',TSTART,TEND,LUPRI)
      call mem_ichor_dealloc(expP)
      call mem_ichor_dealloc(Piprim1)
      call mem_ichor_dealloc(Piprim2)
      call mem_ichor_dealloc(inversexpP)    
      call mem_ichor_dealloc(expA)
      call mem_ichor_dealloc(ContractCoeffA)
      call mem_ichor_dealloc(Acenter)
      call mem_ichor_dealloc(StartOrbitalA)
     ENDIF
    ENDDO !typeA
    call mem_ichor_dealloc(expB)
    call mem_ichor_dealloc(ContractCoeffB)
    call mem_ichor_dealloc(Bcenter)
    call mem_ichor_dealloc(StartOrbitalB)
   ENDDO !typeB
   IF(PermuteRHSTypes)THEN
      !place (ABCD) in (ABDC)
      Call PermuteRHStypesSub(OutputDim1*OutputDim2,OutputDim3,OutputDim4,OutputStorage,&
           & nAtomsC,nAtomsD,startOrbitalC,startOrbitalD,nOrbC,nOrbD,&
           & TriangularRHSAtomLoop,lupri)
   ENDIF
   call mem_ichor_dealloc(noScreenCD)
   call mem_ichor_dealloc(noScreenCD2)
   call mem_ichor_dealloc(expQ)
   call mem_ichor_dealloc(Qiprim1) 
   call mem_ichor_dealloc(Qiprim2) 
   
   call mem_ichor_dealloc(expC)
   call mem_ichor_dealloc(ContractCoeffC)
   call mem_ichor_dealloc(Ccenter)
   call mem_ichor_dealloc(StartOrbitalC)
  ENDDO !typeC
  call mem_ichor_dealloc(expD)
  call mem_ichor_dealloc(ContractCoeffD)
  call mem_ichor_dealloc(Dcenter)
  call mem_ichor_dealloc(StartOrbitalD)
 ENDDO !typeD
ENDDO 
call mem_ichor_dealloc(TABFJW)

call mem_ichor_dealloc(BATCHGAB)
IF(CSScreen)THEN
   call mem_ichor_dealloc(BatchIndexOfTypeA)
   call mem_ichor_dealloc(BatchIndexOfTypeB)
   call mem_ichor_dealloc(MaxGabForTypeAB)
   IF(IchorGabID1.EQ.IchorGabID2)THEN
      !nothing
   ELSE
      call mem_ichor_dealloc(BatchIndexOfTypeC)
      call mem_ichor_dealloc(BatchIndexOfTypeD)
      call mem_ichor_dealloc(BATCHGCD)
      call mem_ichor_dealloc(MaxGabForTypeCD)
   ENDIF
ENDIF

call mem_ichor_dealloc(OrderdListA)
call mem_ichor_dealloc(OrderdListB)
call mem_ichor_dealloc(OrderdListC)
call mem_ichor_dealloc(OrderdListD)
call retrieve_ichor_memvar(MaxMemAllocated,MemAllocated)
IF(INTPRINT.GT.3)THEN
   call stats_ichor_mem(lupri)
ENDIF

end subroutine IchorEri

subroutine IchorTypeIntegralLoop(nAtomsA,nPrimA,nContA,nOrbCompA,startOrbitalA,&
     & iBatchIndexOfTypeA,expA,ContractCoeffA,AngmomA,Acenter,nBatchA,nOrbA,&
     & nAtomsB,nPrimB,nContB,nOrbCompB,startOrbitalB,&
     & iBatchIndexOfTypeB,expB,ContractCoeffB,AngmomB,Bcenter,nBatchB,nOrbB,&
     & nAtomsC,nPrimC,nContC,nOrbCompC,startOrbitalC,&
     & iBatchIndexOfTypeC,expC,ContractCoeffC,AngmomC,Ccenter,nBatchC,nOrbC,&
     & nAtomsD,nPrimD,nContD,nOrbCompD,startOrbitalD,&
     & iBatchIndexOfTypeD,expD,ContractCoeffD,AngmomD,Dcenter,nBatchD,nOrbD,&
     & nPrimP,nPrimQ,nContP,nContQ,expP,expQ,&
     & pcent,qcent,Ppreexpfac,Qpreexpfac,&
     & Qiprim1,Qiprim2,Piprim1,Piprim2,&
     & nTABFJW1,nTABFJW2,TABFJW,TotalAngmom,&
     & CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
     & Qsegmented,Psegmented,TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
     & reducedExponents,integralPrefactor,&
     & PcentPass,Pdistance12Pass,PpreExpFacPass,&
     & Pdistance12,Qdistance12,PQorder,&
     & BATCHGCD,IatomAPass,IatomBPass,BATCHGAB,&
     & LocalInt,nLocalInt,LocalIntPass,Spherical,TmpArray1,TMParray1maxsize,TmpArray2,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,OutputStorage,&
     & TMParray2maxsize,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri)
  implicit none
  logical,intent(in) :: TriangularRHSAtomLoop,TriangularLHSAtomLoop,PermuteLHSTypes
  logical,intent(in) :: Qsegmented,Psegmented,PQorder,Spherical,CSscreen
  logical,intent(in) :: TriangularODAtomLoop
  integer,intent(in) :: nTABFJW1,nTABFJW2,intprint,lupri
  real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2),THRESHOLD_CS
  integer,intent(in) :: Piprim1(nPrimP),Piprim2(nPrimP),Qiprim1(nPrimQ),Qiprim2(nPrimQ)
  !D
  integer,intent(in) :: nAtomsD,nPrimD,nContD,nOrbCompD,AngmomD,nBatchD,nOrbD
  integer,intent(in) :: iBatchIndexOfTypeD
  integer,intent(in) :: startOrbitalD(nAtomsD)
  real(realk),intent(in) :: Dcenter(3,nAtomsD),expD(nPrimD),ContractCoeffD(nPrimD,nContD)
  !C
  integer,intent(in) :: nAtomsC,nPrimC,nContC,nOrbCompC,AngmomC,nBatchC,nOrbC
  integer,intent(in) :: iBatchIndexOfTypeC
  integer,intent(in) :: startOrbitalC(nAtomsC)
  real(realk),intent(in) :: Ccenter(3,nAtomsC),expC(nPrimC),ContractCoeffC(nPrimC,nContC)
  logical,intent(in) :: noScreenCD2(nAtomsC,nAtomsD)
  logical,intent(in) :: noScreenAB(nAtomsA,nAtomsB)
  !Q
  integer,intent(in) :: nContQ,nPrimQ
  real(realk),intent(inout) :: Qcent(3,nPrimQ),Qdistance12(3),QpreExpFac(nPrimQ),expQ(nPrimP)
  !P
  integer,intent(in) :: nContP,nPrimP
  real(realk),intent(inout) :: Pcent(3,nPrimP),Pdistance12(3),PpreExpFac(nPrimP)
  real(realk),intent(in) :: PcentPass(3,nPrimP,natomsA*natomsB),Pdistance12Pass(3,natomsA*natomsB)
  real(realk),intent(in) :: PpreExpFacPass(nPrimP,natomsA*natomsB),expP(nPrimP)
  !
  real(realk),intent(in) :: BATCHGCD(nBatchC,nBatchD)
  !A & B
  integer,intent(in) :: iBatchIndexOfTypeA,iBatchIndexOfTypeB,AngmomB,AngmomA,nOrbA,nOrbB,nPrimA,nPrimB
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbCompA,nOrbCompB,nBatchA,nBatchB,nLocalInt,nContA,nContB
  integer,intent(inout) :: IatomAPass(natomsA*natomsB)
  integer,intent(inout) :: IatomBPass(natomsA*natomsB)
  integer,intent(in) :: startOrbitalA(nAtomsA)
  integer,intent(in) :: startOrbitalB(nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA),expA(nPrimA),ContractCoeffA(nPrimA,nContA)
  real(realk),intent(in) :: Bcenter(3,nAtomsB),expB(nPrimB),ContractCoeffB(nPrimB,nContB)
  real(realk),intent(in) ::  BATCHGAB(nBatchA*nBatchB)
  !collected
  real(realk),intent(inout) :: LocalInt(nLocalInt)
  real(realk),intent(inout) :: LocalIntPass(nOrbCompA*nContA*nAtomsA,nOrbCompB*nContB*nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,TotalAngmom
  real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5)
  integer,intent(in) :: TMParray1maxsize,TMParray2maxsize
  real(realk),intent(inout) :: TmpArray1(TMParray1maxsize),TmpArray2(TMParray2maxsize)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: integralPrefactor(nPrimQ,nPrimP)
  !local variables
  integer :: iBatchD,IatomD,IatomC,startC,IPass,IatomB,IatomA,nPasses,startD
  real(realk) :: DcenterSpec(3),CcenterSpec(3),AcenterSpec(3),BcenterSpec(3),GABELM
  logical :: PermuteRHS

  DO IatomD = 1,nAtomsD
   GABELM = 0.0E0_realk 
   startD = startOrbitalD(iAtomD)
   iBatchD = iBatchIndexOfTypeD + IatomD
   DcenterSpec(1) = Dcenter(1,IAtomD)
   DcenterSpec(2) = Dcenter(2,IAtomD)
   DcenterSpec(3) = Dcenter(3,IAtomD)
   DO IatomC = 1,nAtomsC
    IF(TriangularRHSAtomLoop.AND.IatomD.GT.IatomC)CYCLE
    PermuteRHS = TriangularRHSAtomLoop.AND.IatomD.LT.IatomC
    startC = startOrbitalC(iAtomC)
    IF(noScreenCD2(IatomC,IatomD))THEN
     CcenterSpec(1) = Ccenter(1,IAtomC)
     CcenterSpec(2) = Ccenter(2,IAtomC)
     CcenterSpec(3) = Ccenter(3,IAtomC)
     nPasses = nAtomsA*nAtomsB
     IF(CSscreen)GABELM = BATCHGCD(iBatchIndexOfTypeC+IatomC,iBatchD)
     CALL BUILD_noScreen2(CSscreen,nAtomsA,nAtomsB,&
          & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,&
          & BATCHGAB,THRESHOLD_CS,GABELM,nPasses,IatomAPass,IatomBPass,&
          & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,IatomD,noScreenAB) !Determines nPasses = COUNT(noScreenDC)
!     IF(CSscreen)THEN
!        print*,'CS screenAB:'
!        call loutput(noscreenAB,nAtomsA,nAtomsB,6)
!     ENDIF
     IF(nPasses.NE.nAtomsA*nAtomsB)THEN
        call ichorzero(LocalIntPass,SIZE(LocalIntPass))
     ENDIF
     IF(nPasses.EQ.0)CYCLE
     CALL Build_qcent_Qdistance12_QpreExpFac(nPrimC,nPrimD,nContC,nContD,&
          & expC,expD,CcenterSpec,DcenterSpec,ContractCoeffC,ContractCoeffD,Qsegmented,&
          & Qcent,Qdistance12,QpreExpFac,INTPRINT)
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(Ipass,IatomB,IatomA,AcenterSpec,BcenterSpec,&
!$OMP TmpArray1,TmpArray2,LocalInt,Pcent,Pdistance12,PpreExpFac) SHARED(Qcent,&
!$OMP Qdistance12,QpreExpFac,intprint,lupri,&
!$OMP nPasses,IatomAPass,IatomBPass,LocalIntPass,PcentPass,Pdistance12Pass,&
!$OMP PpreExpFacPass,Acenter,Bcenter,expA,expB,expP,expQ,CONTRACTCOEFFA,&
!$OMP CONTRACTCOEFFB,CONTRACTCOEFFC,CONTRACTCOEFFD,TABFJW,QIPRIM1,QIPRIM2,&
!$OMP PIPRIM1,PIPRIM2,expC,expD,QSEGMENTED,PSEGMENTED,REDUCEDEXPONENTS,&
!$OMP INTEGRALPREFACTOR,PQORDER,CCENTERSPEC,DCENTERSPEC,&
!$OMP SPHERICAL) FIRSTPRIVATE(nPrimP,nAtomsA,nAtomsB,nPrimA,nPrimB,&
!$OMP nPrimC,nPrimD,nPrimQ,nCOntA,nContB,nContC,nContD,nContQ,NTABFJW1,&
!$OMP NTABFJW2,nContP,AngmomA,AngmomB,AngmomC,AngmomD,nLocalInt,TMParray1maxsize,&
!$OMP TMParray2maxsize,TOTALANGMOM,NORBCOMPA,NORBCOMPB,NORBCOMPC,NORBCOMPD) SCHEDULE(DYNAMIC,1)
     DO IPass = 1,nPasses
        IatomB = IatomBPass(IPass)
        IatomA = IatomAPass(IPass)
        AcenterSpec(1) = Acenter(1,IAtomA)
        AcenterSpec(2) = Acenter(2,IAtomA)
        AcenterSpec(3) = Acenter(3,IAtomA)
        BcenterSpec(1) = Bcenter(1,IAtomB)
        BcenterSpec(2) = Bcenter(2,IAtomB)
        BcenterSpec(3) = Bcenter(3,IAtomB)
        CALL Build_pcent_Pdistance12_PpreExpFac2(nPrimP,nPasses,PcentPass,Pdistance12Pass,&
             & PpreExpFacPass,Pcent,Pdistance12,PpreExpFac,nAtomsA,nAtomsB,IatomA,iatomB)
!        write(lupri,'(A,4I)')'CALC atoms:',IatomA,IatomB,IatomC,IatomD
!        write(lupri,'(A,I2,A,I2,A,I2,A,I2,A,I2,A,I2,A,I2,A,I2,A)')&
!             &'(',startOrbitalA(iAtomA)+1,':',startOrbitalA(iAtomA)+nOrbA,&
!             & ',',startOrbitalB(iAtomB)+1,':',startOrbitalB(iAtomB)+nOrbB,&
!             & ',',startC+1,':',startC+nOrbC,&
!             & ',',startD+1,':',startD+nOrbD,')'
        call IchorCoulombIntegral_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,&
             & nPrimQ,nPrimP*nPrimQ,1,1,intprint,lupri,&
             & nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,&
             & ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
             & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
             & Qiprim1,Qiprim2,Piprim1,Piprim2,expA,expB,expC,expD,&
             & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
             & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,&
             & LocalInt,nLocalInt,AcenterSpec,BcenterSpec,CcenterSpec,DcenterSpec,&
             & nAtomsA,nAtomsB,Spherical,TmpArray1,TMParray1maxsize,TmpArray2,&
             & TMParray2maxsize)                   
        !output LocalInt(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP)
        IF(Qsegmented.AND.Psegmented)THEN
         IF(TotalAngmom.NE.0)THEN
          call DistributeToLocalIntPassSeg(LocalInt,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
               & nAtomsA,nAtomsB,iAtomA,iAtomB,LocalIntPass)
         ELSE !TotalAngmom=0
          call DistributeToLocalIntPassSeg0000(LocalInt,nAtomsA,nAtomsB,iAtomA,iAtomB,LocalIntPass)
         ENDIF
        ELSE
         IF(TotalAngmom.NE.0)THEN
          call DistributeToLocalIntPass(LocalInt,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
               & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,iAtomA,iAtomB,LocalIntPass)
         ELSE !TotalAngmom=0
          call DistributeToLocalIntPass0000(LocalInt,&
               & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,iAtomA,iAtomB,LocalIntPass)
         ENDIF
        ENDIF
        !output LocalIntPass(nOrbA*nAtomsA,nOrbB*nAtomsB,nOrbC,nOrbD)
     ENDDO
!$OMP END PARALLEL DO
     !symmetrize the LHS
     IF(TriangularLHSAtomLoop)THEN
        IF(Qsegmented.AND.Psegmented)THEN
         IF(TotalAngmom.NE.0)THEN
          call IchorPermuteLHS(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD,LocalIntPass)
         ELSE
          call IchorPermuteLHSSeg0000(nAtomsA,nAtomsB,LocalIntPass)
         ENDIF
        ELSE
         call IchorPermuteLHS(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD,LocalIntPass)
        ENDIF
     ENDIF
     call IchorDistribute(nAtomsA,nAtomsB,startOrbitalA,startOrbitalB,&
          & startC,startD,OutputStorage,Outputdim1,Outputdim2,Outputdim3,&
          & Outputdim4,&
          & LocalIntPass,nOrbA,nOrbB,nOrbC,nOrbD,PermuteRHS,PermuteLHSTypes,&
          & TriangularODAtomLoop,lupri)         
    ENDIF !ScreenCD
   ENDDO !IatomC
  ENDDO !iAtomD
end subroutine IchorTypeIntegralLoop

subroutine IchorsegsegSSSSIntegralLoop(nAtomsA,nPrimA,startOrbitalA,&
     & iBatchIndexOfTypeA,nBatchA,&
     & nAtomsB,nPrimB,startOrbitalB,iBatchIndexOfTypeB,nBatchB,&
     & nAtomsC,nPrimC,startOrbitalC,iBatchIndexOfTypeC,&
     & expC,ContractCoeffC,Ccenter,nBatchC,&
     & nAtomsD,nPrimD,startOrbitalD,iBatchIndexOfTypeD,&
     & expD,ContractCoeffD,Dcenter,nBatchD,&
     & nPrimP,nPrimQ,qcent,Qpreexpfac,&
     & TABFJW,CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
     & TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
     & reducedExponents,integralPrefactor,&
     & PcentPass,PpreExpFacPass,BATCHGCD,IatomAPass,IatomBPass,BATCHGAB,&
     & LocalIntPass,OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
     & OutputStorage,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri)
  implicit none
  logical,intent(in) :: TriangularRHSAtomLoop,TriangularLHSAtomLoop,PermuteLHSTypes
  logical,intent(in) :: CSscreen,TriangularODAtomLoop
  integer,intent(in) :: nPrimQ,nPrimP,intprint,lupri
  integer,intent(in) :: nAtomsD,nPrimD,nBatchD,iBatchIndexOfTypeD
  integer,intent(in) :: nAtomsC,nPrimC,nBatchC,iBatchIndexOfTypeC
  integer,intent(in) :: startOrbitalD(nAtomsD),startOrbitalC(nAtomsC)
  integer,intent(in) :: startOrbitalA(nAtomsA),startOrbitalB(nAtomsB)
  logical,intent(in) :: noScreenCD2(nAtomsC,nAtomsD),noScreenAB(nAtomsA,nAtomsB)
  integer,intent(in) :: iBatchIndexOfTypeA,iBatchIndexOfTypeB,nPrimA,nPrimB
  integer,intent(in) :: nAtomsA,nAtomsB,nBatchA,nBatchB
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5
  real(realk),intent(in) :: TABFJW(0:3,0:1200),THRESHOLD_CS
  real(realk),intent(in) :: Dcenter(3,nAtomsD),expD(nPrimD),ContractCoeffD(nPrimD)
  real(realk),intent(in) :: Ccenter(3,nAtomsC),expC(nPrimC),ContractCoeffC(nPrimC)
  real(realk),intent(in) :: PcentPass(3,nPrimP,natomsA*natomsB)
  real(realk),intent(in) :: PpreExpFacPass(nPrimP,natomsA*natomsB)
  real(realk),intent(in) :: BATCHGCD(nBatchC,nBatchD)
  real(realk),intent(in) ::  BATCHGAB(nBatchA*nBatchB)
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: integralPrefactor(nPrimQ,nPrimP)
  !used as TMP
  real(realk),intent(inout) :: Qcent(3,nPrimQ),QpreExpFac(nPrimQ)
  integer,intent(inout) :: IatomAPass(natomsA*natomsB)
  integer,intent(inout) :: IatomBPass(natomsA*natomsB)
  real(realk),intent(inout) :: LocalIntPass(nAtomsA,nAtomsB)
  !Output
  real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5)
  !local variables
  integer :: iBatchD,IatomD,IatomC,startC,IPass,IatomB,IatomA,nPasses,startD
  real(realk) :: DcenterSpec(3),CcenterSpec(3),GABELM
  logical :: PermuteRHS
  DO IatomD = 1,nAtomsD
   GABELM = 0.0E0_realk 
   startD = startOrbitalD(iAtomD)
   iBatchD = iBatchIndexOfTypeD + IatomD
   DcenterSpec(1) = Dcenter(1,IAtomD)
   DcenterSpec(2) = Dcenter(2,IAtomD)
   DcenterSpec(3) = Dcenter(3,IAtomD)
   DO IatomC = 1,nAtomsC
    IF(TriangularRHSAtomLoop.AND.IatomD.GT.IatomC)CYCLE
    PermuteRHS = TriangularRHSAtomLoop.AND.IatomD.LT.IatomC
    startC = startOrbitalC(iAtomC)
    IF(noScreenCD2(IatomC,IatomD))THEN
     CcenterSpec(1) = Ccenter(1,IAtomC)
     CcenterSpec(2) = Ccenter(2,IAtomC)
     CcenterSpec(3) = Ccenter(3,IAtomC)
     nPasses = nAtomsA*nAtomsB
     IF(CSscreen)GABELM = BATCHGCD(iBatchIndexOfTypeC+IatomC,iBatchD)
     CALL BUILD_noScreen2(CSscreen,nAtomsA,nAtomsB,&
          & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,&
          & BATCHGAB,THRESHOLD_CS,GABELM,nPasses,IatomAPass,IatomBPass,&
          & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,IatomD,noScreenAB) 
     !Determines nPasses = COUNT(noScreenDC)
     IF(nPasses.NE.nAtomsA*nAtomsB)THEN
        call ichorzero(LocalIntPass,SIZE(LocalIntPass))
     ENDIF
     IF(nPasses.EQ.0)CYCLE
     CALL Build_seg_qcent_QpreExpFac(nPrimC,nPrimD,expC,expD,CcenterSpec,DcenterSpec,&
          & ContractCoeffC,ContractCoeffD,Qcent,QpreExpFac,INTPRINT)
     call IchorCoulombIntegral_seg_seg_SSSS(nPrimP,nPrimQ,nPasses,nAtomsA,nAtomsB,&
          & IatomAPass,iatomBPass,pcentPass,qcent,PpreexpfacPass,Qpreexpfac,TABFJW,&
          & reducedExponents,integralPrefactor,LocalIntPass)
     !symmetrize the LHS
     IF(TriangularLHSAtomLoop)THEN
        call IchorPermuteLHSSeg0000(nAtomsA,nAtomsB,LocalIntPass)
     ENDIF
     call Distribute_seg_seg_SSSS(nAtomsA,nAtomsB,startOrbitalA,startOrbitalB,&
          & startC,startD,OutputStorage,Outputdim1,Outputdim2,Outputdim3,Outputdim4,&
          & LocalIntPass,PermuteRHS,PermuteLHSTypes,TriangularODAtomLoop,lupri)
    ENDIF !ScreenCD
   ENDDO !IatomC
  ENDDO !iAtomD
end subroutine IchorsegsegSSSSIntegralLoop

subroutine determinePermuteSym(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
implicit none
integer,intent(in) :: IchorPermuteSpec
logical,intent(inout) ::  SameLHSaos,SameRHSaos,SameODs
SELECT CASE(IchorPermuteSpec)
CASE(IchorPermuteTTT)
   SameLHSaos=.TRUE.;  SameRHSaos=.TRUE. ; SameODs = .TRUE.
CASE(IchorPermuteFFT)
   SameLHSaos=.FALSE.; SameRHSaos=.FALSE.; SameODs = .TRUE.
CASE(IchorPermuteTTF)
   SameLHSaos=.TRUE.;  SameRHSaos=.TRUE. ; SameODs = .FALSE.
CASE(IchorPermuteTFF)
   SameLHSaos=.TRUE.;  SameRHSaos=.FALSE.; SameODs = .FALSE.
CASE(IchorPermuteFTF)
   SameLHSaos=.FALSE.; SameRHSaos=.TRUE. ; SameODs = .FALSE.
CASE(IchorPermuteFFF)
   SameLHSaos=.FALSE.; SameRHSaos=.FALSE.; SameODs = .FALSE.
CASE DEFAULT
   call ichorquit('unknown case in determinePermuteSym',-1)
END SELECT
end subroutine determinePermuteSym

subroutine PermuteRHStypesSub(dim1dim2,dim3,dim4,OutputStorage,nAtomsC,nAtomsD,&
     & startOrbitalC,startOrbitalD,nOrbC,nOrbD,TriangularRHSAtomLoop,lupri)
  implicit none
  integer,intent(in) :: dim1dim2,dim3,dim4,nAtomsC,nAtomsD,nOrbC,nOrbD,lupri
  integer,intent(in) :: startOrbitalC(nAtomsC),startOrbitalD(nAtomsD)
  real(realk),intent(inout) :: OutputStorage(dim1dim2,dim3,dim4)
  logical,intent(in) :: TriangularRHSAtomLoop
  !
  integer :: I,C,D,startC,startD,IatomC,IatomD,A,B
  IF(TriangularRHSAtomLoop)THEN
!$OMP PARALLEL DEFAULT(none) PRIVATE(IatomD,startD,IatomC,startC,C,D,&
!$OMP I) SHARED(startOrbitalD,startOrbitalC,OutputStorage) FIRSTPRIVATE(nAtomsD,&
!$OMP nAtomsC,nOrbC,nOrbD,dim1dim2) 
   DO IatomC = 1,nAtomsC
    startC = startOrbitalC(iAtomC)
    DO IatomD = 1,IatomC
     startD = startOrbitalD(iAtomD)
!$OMP DO SCHEDULE(DYNAMIC,1)
     DO D=startD+1,startD+nOrbD
      DO C=startC+1,startC+nOrbC
       DO I=1,dim1dim2
        OutputStorage(I,D,C) = OutputStorage(I,C,D)
       ENDDO
      ENDDO
     ENDDO
!$OMP END DO NOWAIT
    ENDDO
   ENDDO
!$OMP END PARALLEL
  ELSE
!$OMP PARALLEL DEFAULT(none) PRIVATE(IatomD,startD,IatomC,startC,C,D,&
!$OMP I) SHARED(startOrbitalD,startOrbitalC,OutputStorage) FIRSTPRIVATE(nAtomsD,&
!$OMP nAtomsC,nOrbC,nOrbD,dim1dim2)
   DO IatomD = 1,nAtomsD
    startD = startOrbitalD(iAtomD)
!$OMP DO SCHEDULE(DYNAMIC,1)
    DO IatomC = 1,nAtomsC
     startC = startOrbitalC(iAtomC)
     DO D=startD+1,startD+nOrbD
      DO C=startC+1,startC+nOrbC
       DO I=1,dim1dim2
        OutputStorage(I,D,C) = OutputStorage(I,C,D)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
!$OMP END DO NOWAIT
   ENDDO
!$OMP END PARALLEL
  ENDIF
end subroutine PermuteRHStypesSub

subroutine PermuteODtypesSub(dim1,dim2,dim3,dim4,OutputStorage,&
     & nAtomsA,nAtomsB,startOrbitalA,startOrbitalB,nOrbA,nOrbB,&
     & nAtomsC,nAtomsD,startOrbitalC,startOrbitalD,nOrbC,nOrbD,&
     & TriangularODAtomLoop,TriangularRHSAtomLoop,&
     & TriangularLHSAtomLoop,SameLHSaos,SameRHSaos,lupri)
  implicit none
  integer,intent(in) :: dim1,dim2,dim3,dim4,lupri
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbA,nOrbB
  integer,intent(in) :: nAtomsC,nAtomsD,nOrbC,nOrbD
  integer,intent(in) :: startOrbitalA(nAtomsA),startOrbitalB(nAtomsB)
  integer,intent(in) :: startOrbitalC(nAtomsC),startOrbitalD(nAtomsD)
  real(realk),intent(inout) :: OutputStorage(dim1,dim2,dim3,dim4)
  logical,intent(in) :: TriangularODAtomLoop,TriangularRHSAtomLoop
  logical,intent(in) :: TriangularLHSAtomLoop,SameLHSaos,SameRHSaos
  !
  integer :: A,B,C,D,startC,startD,IatomC,IatomD,IatomB,IatomA,startA,startB
  integer :: IatomBend
!$OMP PARALLEL DEFAULT(none) PRIVATE(IatomB,IatomBend,iatomD,IatomC,startC,startD,&
!$OMP IatomA,startA,startB,A,B,C,D) SHARED(startOrbitalD,startOrbitalC,startOrbitalA,&
!$OMP startOrbitalB,TriangularODAtomLoop,OutputStorage) FIRSTPRIVATE(nOrbA,nOrbB,nOrbC,&
!$OMP nOrbD,nAtomsC,nAtomsD,nAtomsB,nAtomsA,TriangularLHSAtomLoop,&
!$OMP TriangularRHSAtomLoop,SameRHSaos,SameLHSaos)
  IF(SameLHSaos.AND.SameRHSaos)THEN
   IatomBend = nAtomsB
   DO IatomD = 1,nAtomsD
    startD = startOrbitalD(iAtomD)
    DO IatomC = 1,nAtomsC
     startC = startOrbitalC(iAtomC)
     IF(TriangularRHSAtomLoop.AND.IatomD.GT.IatomC)CYCLE
!$OMP DO SCHEDULE(DYNAMIC,1)
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      IF(TriangularLHSAtomLoop)IatomBend = IatomA
      DO IatomB = 1,IatomBend       
       startB = startOrbitalB(iAtomB)
       IF(TriangularODAtomLoop)THEN
        IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
         DO D=startD+1,startD+nOrbD
          DO C=startC+1,startC+nOrbC
           DO B=startB+1,startB+nOrbB
            DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,B,A) = OutputStorage(A,B,C,D)
              OutputStorage(C,D,B,A) = OutputStorage(A,B,C,D) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ELSE
        DO D=startD+1,startD+nOrbD
         DO C=startC+1,startC+nOrbC
          DO B=startB+1,startB+nOrbB
           DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,B,A) = OutputStorage(A,B,C,D)
              OutputStorage(C,D,B,A) = OutputStorage(A,B,C,D)
           ENDDO !A
          ENDDO !B
         ENDDO !C
        ENDDO !D
       ENDIF
      ENDDO
     ENDDO
!$OMP END DO NOWAIT
    ENDDO
   ENDDO
 ELSEIF(SameLHSaos)THEN
   IatomBend = nAtomsB
   DO IatomD = 1,nAtomsD
    startD = startOrbitalD(iAtomD)
    DO IatomC = 1,nAtomsC
     startC = startOrbitalC(iAtomC)
!$OMP DO SCHEDULE(DYNAMIC,1)
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      IF(TriangularLHSAtomLoop)IatomBend = IatomA
      DO IatomB = 1,IatomBend       
       startB = startOrbitalB(iAtomB)
       IF(TriangularODAtomLoop)THEN
        IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
         DO D=startD+1,startD+nOrbD
          DO C=startC+1,startC+nOrbC
           DO B=startB+1,startB+nOrbB
            DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(C,D,B,A) = OutputStorage(A,B,C,D)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ELSE
        DO D=startD+1,startD+nOrbD
         DO C=startC+1,startC+nOrbC
          DO B=startB+1,startB+nOrbB
           DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(C,D,B,A) = OutputStorage(A,B,C,D)
           ENDDO !A
          ENDDO !B
         ENDDO !C
        ENDDO !D
       ENDIF
      ENDDO
     ENDDO
!$OMP END DO NOWAIT
    ENDDO
   ENDDO
 ELSEIF(SameRHSaos)THEN
   DO IatomD = 1,nAtomsD
    startD = startOrbitalD(iAtomD)
    DO IatomC = 1,nAtomsC
     startC = startOrbitalC(iAtomC)
     IF(TriangularRHSAtomLoop.AND.IatomD.GT.IatomC)CYCLE
!$OMP DO SCHEDULE(DYNAMIC,1)
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      DO IatomB = 1,nAtomsB
       startB = startOrbitalB(iAtomB)
       IF(TriangularODAtomLoop)THEN
        IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
         DO D=startD+1,startD+nOrbD
          DO C=startC+1,startC+nOrbC
           DO B=startB+1,startB+nOrbB
            DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,A,B) = OutputStorage(A,B,C,D)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ELSE
        DO D=startD+1,startD+nOrbD
         DO C=startC+1,startC+nOrbC
          DO B=startB+1,startB+nOrbB
           DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
              OutputStorage(D,C,A,B) = OutputStorage(A,B,C,D)
           ENDDO !A
          ENDDO !B
         ENDDO !C
        ENDDO !D
       ENDIF
      ENDDO
     ENDDO
!$OMP END DO NOWAIT
    ENDDO
   ENDDO
 ELSE
   DO IatomD = 1,nAtomsD
    startD = startOrbitalD(iAtomD)
    DO IatomC = 1,nAtomsC
     startC = startOrbitalC(iAtomC)
!$OMP DO SCHEDULE(DYNAMIC,1)
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      DO IatomB = 1,nAtomsB
       startB = startOrbitalB(iAtomB)
       IF(TriangularODAtomLoop)THEN
        IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
         DO D=startD+1,startD+nOrbD
          DO C=startC+1,startC+nOrbC
           DO B=startB+1,startB+nOrbB
            DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ELSE
        DO D=startD+1,startD+nOrbD
         DO C=startC+1,startC+nOrbC
          DO B=startB+1,startB+nOrbB
           DO A=startA+1,startA+nOrbA
!              WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,F16.8)')&
!                   & 'CALC: (',A,',',B,',',C,',',D,')=',OutputStorage(A,B,C,D)
              OutputStorage(C,D,A,B) = OutputStorage(A,B,C,D)
           ENDDO !A
          ENDDO !B
         ENDDO !C
        ENDDO !D
       ENDIF
      ENDDO
     ENDDO
!$OMP END DO NOWAIT
    ENDDO
   ENDDO
 ENDIF
!$OMP END PARALLEL
end subroutine PermuteODtypesSub

subroutine ExtractBatchGabFromFullGab(nBatchASmall,nBatchBSmall,&
     & BATCHGABSmall,nBatchABIG,nBatchBBIG,BATCHGABBig,&
     & startBatchA,endBatchA,startBatchB,endBatchB)
implicit none
integer,intent(in) :: nBatchASmall,nBatchBSmall
integer,intent(in) :: nBatchABIG,nBatchBBIG
integer,intent(in) :: startBatchA,endBatchA,startBatchB,endBatchB
real(realk),intent(in) :: BATCHGABBig(nBatchABIG,nBatchBBIG) 
real(realk),intent(inout) :: BATCHGABSmall(nBatchASmall,nBatchBSmall) 
!local
integer :: sA,sB,IBB,IA,IB
sA = startBatchA-1
sB = startBatchB-1
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(IB,IBB,&
!$OMP IA) SHARED(BATCHGABSmall,BATCHGABBig) FIRSTPRIVATE(sA,&
!$OMP sB,nBatchASmall,nBatchBSmall) SCHEDULE(DYNAMIC,3)
DO IB=1,nBatchBSmall
   IBB = sB+IB
   DO IA=1,nBatchASmall
      BATCHGABSmall(IA,IB) = BATCHGABBig(sA+IA,IBB)
   ENDDO
ENDDO 
!$OMP END PARALLEL DO 
end subroutine ExtractBatchGabFromFullGab

subroutine Distribute_seg_seg_SSSS(nAtomsA,nAtomsB,startOrbitalA,startOrbitalB,&
     & startC,startD,OutputStorage,dim1,dim2,dim3,dim4,LocalIntPassSSSS,PermuteRHS,&
     & PermuteLHSTypes,TriangularODAtomLoop,lupri)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,startC,startD,dim1,dim2,dim3,dim4
  integer,intent(in) :: lupri
  integer,intent(in) :: startOrbitalA(nAtomsA),startOrbitalB(nAtomsB)
  logical,intent(in) :: PermuteRHS,PermuteLHSTypes,TriangularODAtomLoop
  real(realk),intent(inout) :: OutputStorage(dim1,dim2,dim3,dim4)
  real(realk),intent(in) :: LocalIntPassSSSS(nAtomsA,nAtomsB)
  !local variables
  integer :: i4,i3,iatomb,startb,i2,iatomA,startA
  I4 = startD + 1
  I3 = startC + 1
  IF(PermuteLHSTypes)THEN
   IF(PermuteRHS)THEN
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(iatomB,I2,&
!$OMP IatomA,startA) FIRSTPRIVATE(I4,I3,nAtomsA,&
!$OMP nAtomsB,PermuteLHSTypes,PermuteRHS) SHARED(startOrbitalB,&
!$OMP startOrbitalA,LocalIntPassSSSS,OutputStorage) SCHEDULE(DYNAMIC,1)
    DO IatomB = 1,nAtomsB
     I2 = startOrbitalB(iAtomB) + 1
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      OutputStorage(1 + startA,I2,I3,I4)=LocalIntPassSSSS(iAtomA,iAtomB)
      OutputStorage(1 + startA,I2,I4,I3)=LocalIntPassSSSS(iAtomA,iAtomB)
      OutputStorage(I2,1 + startA,I3,I4)=LocalIntPassSSSS(iAtomA,iAtomB)
      OutputStorage(I2,1 + startA,I4,I3)=LocalIntPassSSSS(iAtomA,iAtomB)
     ENDDO
    ENDDO
!$OMP END PARALLEL DO
   ELSE  !PermuteLHSTypes NOT PermuteRHS
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(iatomB,I2,&
!$OMP IatomA,startA) FIRSTPRIVATE(I4,I3,nAtomsA,&
!$OMP nAtomsB,PermuteLHSTypes,PermuteRHS) SHARED(startOrbitalB,&
!$OMP startOrbitalA,LocalIntPassSSSS,OutputStorage) SCHEDULE(DYNAMIC,1)
    DO IatomB = 1,nAtomsB
     I2 = startOrbitalB(iAtomB) + 1
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      OutputStorage(1 + startA,I2,I3,I4)=LocalIntPassSSSS(iAtomA,iAtomB)
      OutputStorage(I2,1 + startA,I3,I4)=LocalIntPassSSSS(iAtomA,iAtomB)
     ENDDO
    ENDDO
!$OMP END PARALLEL DO
   ENDIF
  ELSE !NOT PermuteLHSTypes
   IF(PermuteRHS)THEN
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(iatomB,I2,&
!$OMP IatomA,startA) FIRSTPRIVATE(I4,I3,nAtomsA,&
!$OMP nAtomsB,PermuteLHSTypes,PermuteRHS) SHARED(startOrbitalB,&
!$OMP startOrbitalA,LocalIntPassSSSS,OutputStorage) SCHEDULE(DYNAMIC,1)
    DO IatomB = 1,nAtomsB
     I2 = startOrbitalB(iAtomB) + 1
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      OutputStorage(1 + startA,I2,I3,I4)=LocalIntPassSSSS(iAtomA,iAtomB)
      OutputStorage(1 + startA,I2,I4,I3)=LocalIntPassSSSS(iAtomA,iAtomB)
     ENDDO
    ENDDO
!$OMP END PARALLEL DO
   ELSE !NOT PermuteLHSTypes NOT PermuteRHS 
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(iatomB,I2,&
!$OMP IatomA,startA) FIRSTPRIVATE(I4,I3,nAtomsA,&
!$OMP nAtomsB,PermuteLHSTypes,PermuteRHS) SHARED(startOrbitalB,&
!$OMP startOrbitalA,LocalIntPassSSSS,OutputStorage) SCHEDULE(DYNAMIC,1)
    DO IatomB = 1,nAtomsB
     I2 = startOrbitalB(iAtomB) + 1
     DO IatomA = 1,nAtomsA
      startA = startOrbitalA(iAtomA)
      OutputStorage(1 + startA,I2,I3,I4)=LocalIntPassSSSS(iAtomA,iAtomB)
     ENDDO
    ENDDO
!$OMP END PARALLEL DO
   ENDIF
  ENDIF
end subroutine Distribute_seg_seg_SSSS

subroutine DistributeToLocalIntPass(LocalInt,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,iAtomA,iAtomB,LocalIntPass1)
  implicit none 
  integer,intent(in) :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
  integer,intent(in) :: nContA,nContB,nContC,nContD,iAtomA,iAtomB
  real(realk),intent(in) :: LocalInt(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContC*nContD,nContA*nContB)
  real(realk),intent(inout) :: LocalIntPass1(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
  !local variables
  integer :: iContQ,iContA,iContB,iContC,iContD,iContP,offset
  integer :: iAngA,iAngB,iAngC,iAngD,I2,I3,I4
  iContQ = 0
  DO iContD = 1,nContD
   DO iContC = 1,nContC
    iContQ = iContQ + 1
    DO iAngD = 1,nOrbCompD
     I4 = iAngD + (iContD-1)*nOrbCompD
     DO iAngC = 1,nOrbCompC
      I3 = iAngC + (iContC-1)*nOrbCompC
      DO iContB = 1,nContB
       DO iAngB = 1,nOrbCompB
        I2 = iAngB + (iContB-1)*nOrbCompB
        DO iContA = 1,nContA
         iContP = iContA+(iContB-1)*nContA
         offset = (iContA-1)*nOrbCompA
         DO iAngA = 1,nOrbCompA
          LocalIntPass1(iAngA + offset,iatomA,I2,iatomB,I3,I4) = LocalInt(iAngA,iAngB,iAngC,iAngD,iContQ,iContP)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
end subroutine DistributeToLocalIntPass

subroutine DistributeToLocalIntPass0000(LocalInt,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,iAtomA,iAtomB,LocalIntPass2)
  implicit none 
  integer,intent(in) :: nAtomsA,nAtomsB
  integer,intent(in) :: nContA,nContB,nContC,nContD,iAtomA,iAtomB
  real(realk),intent(in) :: LocalInt(nContC*nContD,nContA,nContB)
  real(realk),intent(inout) :: LocalIntPass2(nContA,nAtomsA,nContB,nAtomsB,nContC*nContD)
  !local variables
  integer :: iContQ,iContA,iContB
  DO iContQ = 1,nContD*nContC
   DO iContB = 1,nContB
    DO iContA = 1,nContA
     LocalIntPass2(iContA,iatomA,iContB,iatomB,iContQ) = LocalInt(iContQ,iContA,iContB)
    ENDDO
   ENDDO
  ENDDO
end subroutine DistributeToLocalIntPass0000

subroutine DistributeToLocalIntPassSeg(LocalInt,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nAtomsA,nAtomsB,iAtomA,iAtomB,LocalIntPass3)
  implicit none 
  integer,intent(in) :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
  integer,intent(in) :: iAtomA,iAtomB
  real(realk),intent(in) :: LocalInt(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
  real(realk),intent(inout) :: LocalIntPass3(nOrbCompA,nAtomsA,nOrbCompB,nAtomsB,nOrbCompC,nOrbCompD)
  !local variables
  integer :: iAngA,iAngB,iAngC,iAngD
  DO iAngD = 1,nOrbCompD
   DO iAngC = 1,nOrbCompC
    DO iAngB = 1,nOrbCompB
     DO iAngA = 1,nOrbCompA
      LocalIntPass3(iAngA,iatomA,iAngB,iatomB,iAngC,iAngD) = LocalInt(iAngA,iAngB,iAngC,iAngD)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
end subroutine DistributeToLocalIntPassSeg

subroutine DistributeToLocalIntPassSeg0000(LocalInt,nAtomsA,nAtomsB,iAtomA,iAtomB,LocalIntPass4)
  implicit none 
  integer,intent(in) :: nAtomsA,nAtomsB
  integer,intent(in) :: iAtomA,iAtomB
  real(realk),intent(in) :: LocalInt(1)
  real(realk),intent(inout) :: LocalIntPass4(nAtomsA,nAtomsB)
  LocalIntPass4(iatomA,iatomB) = LocalInt(1)
end subroutine DistributeToLocalIntPassSeg0000

subroutine IchorPermuteLHS(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD,LocalIntPass5)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB
  integer,intent(in) :: nOrbA,nOrbB,nOrbC,nOrbD
  real(realk),intent(inout) :: LocalIntPass5(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC*nOrbD)  
  !local variables
  integer :: iOrbQ,iAtomB,iOrbB,iatomA,iOrbA
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(iOrbQ,iatomB,iOrbB,IatomA,&
!$OMP iOrbA) FIRSTPRIVATE(nAtomsB,nAtomsA,nOrbA,nOrbB,nOrbC,&
!$OMP nOrbD) SHARED(LocalIntPass5) SCHEDULE(DYNAMIC,1)
  DO iOrbQ = 1,nOrbC*nOrbD
   DO IatomB = 1,nAtomsB
    DO iOrbB = 1,nOrbB
     DO IatomA = 1,IatomB-1
      DO iOrbA = 1,nOrbA
       LocalIntPass5(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)=LocalIntPass5(iOrbB,iAtomB,iOrbA,iAtomA,iOrbQ)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
end subroutine IchorPermuteLHS

subroutine IchorPermuteLHSSeg0000(nAtomsA,nAtomsB,LocalIntPass5)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB
  real(realk),intent(inout) :: LocalIntPass5(nAtomsA,nAtomsB)
  !local variables
  integer :: iAtomB,iatomA
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(iatomB,&
!$OMP IatomA) FIRSTPRIVATE(nAtomsB) SHARED(LocalIntPass5) SCHEDULE(DYNAMIC,1)
  DO IatomB = 1,nAtomsB
     DO IatomA = 1,IatomB-1
        LocalIntPass5(iAtomA,iAtomB)=LocalIntPass5(iAtomB,iAtomA)
     ENDDO
  ENDDO
!$OMP END PARALLEL DO
end subroutine IchorPermuteLHSSeg0000

subroutine IchorDistribute(nAtomsA,nAtomsB,startOrbitalA,startOrbitalB,&
     & startC,startD,OutputStorage,dim1,dim2,dim3,dim4,&
     & LocalIntPass6,nOrbA,nOrbB,nOrbC,nOrbD,PermuteRHS,PermuteLHSTypes,&
     & TriangularODAtomLoop,lupri)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,startC,startD,dim1,dim2,dim3,dim4
  integer,intent(in) :: nOrbA,nOrbB,nOrbC,nOrbD,lupri
  integer,intent(in) :: startOrbitalA(nAtomsA),startOrbitalB(nAtomsB)
  logical,intent(in) :: PermuteRHS,PermuteLHSTypes,TriangularODAtomLoop
  real(realk),intent(inout) :: OutputStorage(dim1,dim2,dim3,dim4)
  real(realk),intent(in) :: LocalIntPass6(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
  !local variables
  integer :: iOrbD,i4,iOrbC,i3,iatomb,startb,iOrbB,i2,iatomA,startA,iOrbA
!$OMP PARALLEL DEFAULT(none) PRIVATE(iOrbD,I4,iOrbC,I3,iatomB,startB,iOrbB,I2,IatomA,&
!$OMP iOrbA,startA) FIRSTPRIVATE(nOrbD,nOrbC,nOrbB,nOrbA,nAtomsA,nAtomsB,PermuteLHSTypes,&
!$OMP PermuteRHS) SHARED(startOrbitalB,startOrbitalA,startC,startD,LocalIntPass6,OutputStorage)
  IF(PermuteLHSTypes)THEN
   IF(PermuteRHS)THEN
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!$OMP DO SCHEDULE(DYNAMIC,1)
       DO IatomB = 1,nAtomsB
        startB = startOrbitalB(iAtomB)
        DO iOrbB = 1,nOrbB
         I2 = startB + iOrbB
         DO IatomA = 1,nAtomsA
          startA = startOrbitalA(iAtomA)
          DO iOrbA = 1,nOrbA
           OutputStorage(iOrbA + startA,I2,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
           OutputStorage(iOrbA + startA,I2,I4,I3)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
           OutputStorage(I2,iOrbA + startA,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
           OutputStorage(I2,iOrbA + startA,I4,I3)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ELSE  !PermuteLHSTypes NOT PermuteRHS
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!$OMP DO SCHEDULE(DYNAMIC,1)
       DO IatomB = 1,nAtomsB
        startB = startOrbitalB(iAtomB)
        DO iOrbB = 1,nOrbB
         I2 = startB + iOrbB
         DO IatomA = 1,nAtomsA
          startA = startOrbitalA(iAtomA)
          DO iOrbA = 1,nOrbA
           OutputStorage(iOrbA + startA,I2,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
           OutputStorage(I2,iOrbA + startA,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ENDIF
  ELSE !NOT PermuteLHSTypes
   IF(PermuteRHS)THEN
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!$OMP DO SCHEDULE(DYNAMIC,1)
       DO IatomB = 1,nAtomsB
        startB = startOrbitalB(iAtomB)
        DO iOrbB = 1,nOrbB
         I2 = startB + iOrbB
         DO IatomA = 1,nAtomsA
          startA = startOrbitalA(iAtomA)
          DO iOrbA = 1,nOrbA
           OutputStorage(iOrbA + startA,I2,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
           OutputStorage(iOrbA + startA,I2,I4,I3)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ELSE !NOT PermuteLHSTypes NOT PermuteRHS 
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!$OMP DO SCHEDULE(DYNAMIC,1)
       DO IatomB = 1,nAtomsB
        startB = startOrbitalB(iAtomB)
        DO iOrbB = 1,nOrbB
         I2 = startB + iOrbB
         DO IatomA = 1,nAtomsA
          startA = startOrbitalA(iAtomA)
          DO iOrbA = 1,nOrbA
           OutputStorage(iOrbA + startA,I2,I3,I4)=LocalIntPass6(iOrbA,iAtomA,iOrbB,iAtomB,iOrbC,iOrbD)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ENDIF
  ENDIF
!$OMP END PARALLEL
 end subroutine IchorDistribute

END MODULE IchorErimodule
