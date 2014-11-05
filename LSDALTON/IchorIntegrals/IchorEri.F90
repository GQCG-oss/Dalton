!> @file
!> Contains the main Ichor integral drivers for calculation electron repulsion integrals
!> Ichor is the "Integral Code Hand Optimized for Rapid evaluation" 
!> Ichor is the ethereal golden fluid that is the blood of the greek gods
!> The code was written by Thomas Kjaergaard with inspiration from the 
!> McMurchie-Davidson algorithm based Thermite code written by Simen Reine and Thomas Kjaergaard and 
!> as well as the Obara-Saika algorithm based Interest code written by Michal Repisky.
!> This code is based on Obara-Saika (and Head-Gordon-Pople) 
!> The Interest code was also used for debugging purposes
!> \brief Main Ichor drivers for the calculation of integrals 
!> based on the Obara Saika(OS)/Head-Gordon-Pople(HGP)
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorErimod
  use IchorprecisionMod
  use IchorCommonMod
  use IchorBatchToolsMod
  use IchorEriCoulombintegralCPUMcMGeneralMod, only: TmpArray3,TmpArray4,DetermineSizeTmpArray34,&
       & precalcichorsphmat, freeichorsphmat,nTmpArray3,nTmpArray4
  use IchorEriCoulombintegralCPUOBSGeneralMod, only: ICI_CPU_OBS_general, &
       & ICI_CPU_OBS_general_size
  use IchorEriCoulombintegralGPUOBSGeneralMod, only: ICI_GPU_OBS_general, &
       & ICI_GPU_OBS_general_size
  use ICI_seg_seg_SSSS_mod, only: ICI_seg_seg_SSSS
  use IchorMemory
  use IchorGammaTabulationMod
  use IchorParametersMod
  use IchorSaveGabMod
  use IchorInputInfoMod
  use IchorEriToolsmod
  use IchorEriLinkmod
  use IchorEriDistmod
#ifdef VAR_OPENACC
  !OpenACC libary routines  
  use openacc!, only: acc_async_test
#endif
  use AGC_GPU_OBS_TRParamMod
  use IchorGaussianGeminalMod, only: set_GGem, free_GGem, GGemOperatorCalc

  public:: IchorEri,IchorEriMem
  private

logical,parameter :: UseCPU = .TRUE.

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
     & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorOperatorSpec,&
     & IchorInputDim1,IchorInputDim2,IchorInputDim3,&
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
Integer,intent(in) :: nContOfTypeB(ntypesB),nPrimOfTypeB(ntypesB)
Integer,intent(in) :: startOrbitalOfTypeB(MaxNatomsB,ntypesB)
Real(realk),intent(in) :: Bcenters(3,MaxNatomsB,ntypesB),exponentsOfTypeB(MaxnprimB,ntypesB)
Real(realk),intent(in) :: ContractCoeffOfTypeB(MaxnprimB,MaxnContB,ntypesB)
!
! Same for Center C
!
integer,intent(in) :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC,startBatchC,endBatchC
Integer,intent(in) :: AngmomOfTypeC(ntypesC),nAtomsOfTypeC(ntypesC)
Integer,intent(in) :: nContOfTypeC(ntypesC),nPrimOfTypeC(ntypesC)
Integer,intent(in) :: startOrbitalOfTypeC(MaxNatomsC,ntypesC)
Real(realk),intent(in) :: Ccenters(3,MaxNatomsC,ntypesC),exponentsOfTypeC(MaxnprimC,ntypesC)
Real(realk),intent(in) :: ContractCoeffOfTypeC(MaxnprimC,MaxnContC,ntypesC)
!
! Same for Center D
!
integer,intent(in) :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD,startBatchD,endBatchD
Integer,intent(in) :: AngmomOfTypeD(ntypesD),nAtomsOfTypeD(ntypesD)
Integer,intent(in) :: nContOfTypeD(ntypesD),nPrimOfTypeD(ntypesD)
Integer,intent(in) :: startOrbitalOfTypeD(MaxNatomsD,ntypesD)
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
!> Operator Specification (IchorOperatorSpec = CoulombOperator = 1) means Coulomb operator
Integer,intent(in) :: IchorOperatorSpec
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
real(realk),intent(inout)::OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5)
!> Logical unit number of output file.
Integer,intent(in) :: lupri
! Local variables
integer :: nPrimP,nContP,nPrimQ,nContQ
integer :: nTABFJW1,nTABFJW2,i1,i2,i3,i4,AngmomQ,TotalAngmom
integer :: i12,offset,K,I,iPrimQ,iPrimP,icont,AngmomP
integer :: oldmaxangmomABCD
integer :: ItypeA,ItypeB,itypeC,itypeD,AngmomA,AngmomB,AngmomC,AngmomD
integer :: ItypeAnon,ItypeBnon,itypeCnon,itypeDnon,nLocalInt
integer :: nPrimA,nPrimB,nContA,nAtomsA,nAtomsB,nAtomsC,nAtomsD,nOrbA
integer :: nDimA,nOrbCompA,nContB,nOrbB,nDimB
integer :: nPrimC,nContC,nPrimD,nContD,nOrbC,nOrbD,nOrbCompB,nOrbCompC,nDimC
integer :: nDimD,nOrbCompD,INTPRINT,maxangmomABCD
integer :: nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD
integer :: nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV
integer,allocatable :: Qiprim1(:),Qiprim2(:)
integer,allocatable :: OrderdListA(:),OrderdListB(:),OrderdListC(:),OrderdListD(:)
logical :: Psegmented,Qsegmented,PQorder,Spherical,SameLHSaos,SameRHSaos,SameODs
logical :: TriangularLHSAtomLoop,TriangularRHSAtomLoop,PermuteRHS,CSScreen
logical :: TriangularODAtomLoop
logical :: NOTDoSSSS,Segmented,PermuteLHSTypes,PermuteRHSTypes,PermuteODTypes
real(realk),allocatable :: expP(:),Pcent(:),PpreExpFac(:),Pdistance12Pass(:,:)!,LocalInt(:)
real(realk),allocatable :: expQ(:),PcentPass(:,:),PpreExpFacPass(:,:)
real(realk),allocatable :: inversexpP(:),LocalIntPass(:,:)
real(realk),allocatable :: QpreExpFac(:),Qcent(:)
REAL(realk),allocatable :: TABFJW(:,:),reducedExponents(:),integralPrefactor(:)
real(realk) :: AcenterSpec(3),BcenterSpec(3),CcenterSpec(3),DcenterSpec(3)
real(realk) :: Qdistance12(3)
real(realk),allocatable :: expA(:),ContractCoeffA(:,:),Acenter(:,:)
real(realk),allocatable :: expB(:),ContractCoeffB(:,:),Bcenter(:,:)
real(realk),allocatable :: expC(:),ContractCoeffC(:,:),Ccenter(:,:)
real(realk),allocatable :: expD(:),ContractCoeffD(:,:),Dcenter(:,:)
integer,allocatable :: startOrbitalA(:),startOrbitalB(:)
integer,allocatable :: startOrbitalC(:),startOrbitalD(:)
!Tmporary array used in the innermost routines 
integer :: TMParray1maxsize,TMParray2maxsize
!real(realk),allocatable :: TmpArray1(:)
!real(realk),allocatable :: TmpArray2(:)
!Batch info
integer :: iBatchC,iBatchD,nBatchA,nBatchB,nBatchC,nBatchD
integer,allocatable :: BatchIndexOfTypeA(:)
integer,allocatable :: BatchIndexOfTypeB(:)
integer,allocatable :: BatchIndexOfTypeC(:)
integer,allocatable :: BatchIndexOfTypeD(:)
integer :: iBatchIndexOfTypeA,iBatchIndexOfTypeB,iBatchIndexOfTypeC,iBatchIndexOfTypeD
!Screening variables
real(realk),allocatable :: MaxGabForTypeAB(:,:),MaxGabForTypeCD(:,:)
real(realk),allocatable,target :: BATCHGAB(:,:)
real(realk),allocatable :: BATCHGAB2(:,:),BATCHGCD2(:,:)
real(realk),pointer :: BATCHGCD(:,:)
logical,allocatable :: noScreenAB(:,:),noScreenCD(:,:),noScreenCD2(:,:)
integer :: nBatchAGAB,nBatchBGAB,nBatchCGCD,nBatchDGCD
real(realk) :: GABELM
!Passes
integer :: nPasses,iPass,MaxPasses!,nPassLoop,nPassesInPassLoop,iPassTMP !,ndimPass
!integer :: iPassLoop
integer :: MaxTotalAngmom,IAngmomTypes
!ODscreening
real(realk) :: ExtentOfTypeA(ntypesA),ExtentOfTypeB(ntypesB),extentA,extentB,extentC,extentD
real(realk) :: ExtentOfTypeC(ntypesC),ExtentOfTypeD(ntypesD)
real(realk) :: sumExtent2AB,sumExtent2CD,Tstart,Tend
logical :: ODscreen,QQRscreen,dolink
character(len=16)  :: TYPESTRING
!Link specific stuff
real(realk) :: MaxGabLHS,MaxGabRHS !MaximumAllover
real(realk) :: MaxAtomGabelmLHS,MaxAtomGabelmRHS !MaximumAllover
real(realk),allocatable :: AtomGCD(:,:),MaxGCDvec(:),AtomGAB(:,:),MaxGABvec(:)
real(realk),allocatable :: ReducedDmatBD(:,:),ReducedDmatBC(:,:)
real(realk),allocatable :: DmatBD(:,:,:),DmatAD(:,:,:),DmatBC(:,:,:),DmatAC(:,:,:)
real(realk),allocatable :: KmatBD(:,:,:),KmatAD(:,:,:),KmatBC(:,:,:),KmatAC(:,:,:)
integer,allocatable :: nKetList(:),KetList(:,:),nBraList(:),BraList(:,:)
integer,allocatable :: nBraketList(:),BraketList(:,:)
logical :: doMoTrans
!MOtrans specific stuff
integer :: nCMO1,nCMO2,nCMO3,nCMO4
integer :: nStaticParamIfac
real(realk),allocatable :: CMO1A(:,:),CMO2B(:,:)
real(realk),allocatable :: CMO1B(:,:),CMO2A(:,:)
real(realk),allocatable :: CMO3C(:,:),CMO4D(:,:)
real(realk),allocatable :: CMO3D(:,:),CMO4C(:,:)

IF(.NOT.UseCPU)THEN
   Write(lupri,'(A,F10.3,A)')'Ichor: GPU Maximum Memory : ', IchorGPUMAXMEM, ' GB'
ENDIF
IF(.NOT.UseCPU) call Init_AGC_TransferParam()
#ifdef VAR_OPENACC
!$ACC DATA COPYIN(TUVindexX1_35,TUVindexX2_35,TUVindexX3_35,&
!$ACC             TUVindexX1_56,TUVindexX2_56,TUVindexX3_56,&
!$ACC             TUVindexX1_84,TUVindexX2_84,TUVindexX3_84,&
!$ACC             TUVindexX1_120,TUVindexX2_120,TUVindexX3_120,&
!$ACC             IfacX1_20,IfacX2_20,IfacX3_20,&
!$ACC             IfacX1_35,IfacX2_35,IfacX3_35,&
!$ACC             IfacX1_56,IfacX2_56,IfacX3_56,&
!$ACC             IfacX1_84,IfacX2_84,IfacX3_84)
nStaticParamIfac = 3*(35 + 20 + 56 + 35 + 84 + 56 + 120 + 84)
#endif

doMOtrans = IchorJobSpec.EQ.IchorJobMOtrans
IF(doMOtrans)Then
   nCMO1 = InputM(InputMspec(1))%n2
   nCMO2 = InputM(InputMspec(2))%n2
   nCMO3 = InputM(InputMspec(3))%n2
   nCMO4 = InputM(InputMspec(4))%n2
   IF(nCMO1.NE.OutputDim1)call ichorQuit('MoTrans Dim Mismatch1',-1)
   IF(nCMO2.NE.OutputDim2)call ichorQuit('MoTrans Dim Mismatch2',-1)
   IF(nCMO3.NE.OutputDim3)call ichorQuit('MoTrans Dim Mismatch3',-1)
   IF(nCMO4.NE.OutputDim4)call ichorQuit('MoTrans Dim Mismatch4',-1)
   IF(OutputDim5.NE.1)call ichorQuit('MoTrans Dim Mismatch5',-1)
ENDIF

! FIXME 
! LinK requires GAB - if not CSscreen fix
doLink = IchorJobSpec.EQ.IchorJobLink
IF(doLink)Then
   IF(OutputDim5.NE.IchorInputDim3)call ichorQuit('Link Dim Mismatch5',-1)
ENDIF

! POSSIBLE IMPROVEMENTS / TO DO / TODO LIST:
! FFFF
! FUSE LOOPS IN TRANSFER
! Replace HARDCODED LOOPUNROLLED WITH LOOP?
! MKL primitives VML functions - google "List of VML Functions" ?
! Get CPI rate down! below 0.5!
! FIX noscreen not pretty to have 2 noscreen!
! nPasses collect - do not collect? PROFILE
! Try TABFJW(0:3,1200) => TABFJW(0:6,120) Meaning size(4800) => size(840) 
! At the expense of RJ000(2,I) = TABFJW(0,I) + .. +  TABFJW(6,I) (6 TERMS)
! TWO IchorTypeIntegralLoop subs with and without OMP (or ifdefs)
! in primitiveContraction4Gen loop unroll iTUV=1,4!
! S only primitive screening?
call determineScreening(IchorScreenSpec,CSscreen,ODscreen,QQRscreen)
call build_ichor_AOextent(MaxnAtomsA,MaxnprimA,MaxnContA,ntypesA,exponentsOfTypeA,ODscreen,&
     & nContOfTypeA,nPrimOfTypeA,ContractCoeffOfTypeA,AngmomOfTypeA,Threshold_OD,ExtentOfTypeA)
call build_ichor_AOextent(MaxnAtomsB,MaxnprimB,MaxnContB,ntypesB,exponentsOfTypeB,ODscreen,&
     & nContOfTypeB,nPrimOfTypeB,ContractCoeffOfTypeB,AngmomOfTypeB,Threshold_OD,ExtentOfTypeB)
call build_ichor_AOextent(MaxnAtomsC,MaxnprimC,MaxnContC,ntypesC,exponentsOfTypeC,ODscreen,&
     & nContOfTypeC,nPrimOfTypeC,ContractCoeffOfTypeC,AngmomOfTypeC,Threshold_OD,ExtentOfTypeC)
call build_ichor_AOextent(MaxnAtomsD,MaxnprimD,MaxnContD,ntypesD,exponentsOfTypeD,ODscreen,&
     & nContOfTypeD,nPrimOfTypeD,ContractCoeffOfTypeD,AngmomOfTypeD,Threshold_OD,ExtentOfTypeD)
call set_ichor_memvar(MaxMemAllocated,MemAllocated,MaxMem)
INTPRINT=IchorDebugSpec
IF(INTPRINT.GT.50)WRITE(LUPRI,'(A)') ' IchorEri'
allocate(OrderdListA(nTypesA))
call mem_ichor_alloc(OrderdListA)
call GenerateOrderdListOfTypes(lupri,nTypesA,AngmomOfTypeA,OrderdListA)
allocate(OrderdListB(nTypesB))
call mem_ichor_alloc(OrderdListB)
call GenerateOrderdListOfTypes(lupri,nTypesB,AngmomOfTypeB,OrderdListB)
allocate(OrderdListC(nTypesC))
call mem_ichor_alloc(OrderdListC)
call GenerateOrderdListOfTypes(lupri,nTypesC,AngmomOfTypeC,OrderdListC)
allocate(OrderdListD(nTypesD))
call mem_ichor_alloc(OrderdListD)
call GenerateOrderdListOfTypes(lupri,nTypesD,AngmomOfTypeD,OrderdListD)
PQorder=.FALSE.
IF(CSScreen)THEN
   allocate(BatchIndexOfTypeA(nTypesA))
   call mem_ichor_alloc(BatchIndexOfTypeA)
   call ConstructBatchIndexOfType(BatchIndexOfTypeA,nTypesA,nAtomsOfTypeA,nBatchA)
   allocate(BatchIndexOfTypeB(nTypesB))
   call mem_ichor_alloc(BatchIndexOfTypeB)
   call ConstructBatchIndexOfType(BatchIndexOfTypeB,nTypesB,nAtomsOfTypeB,nBatchB)
   call RetrieveGabDimFromIchorSaveGabModule(nBatchAGAB,nBatchBGAB,IchorGabID1)
   allocate(BATCHGAB(nBatchA,nBatchB))
   call mem_ichor_alloc(BATCHGAB)
   IF(nBatchA.NE.endBatchA-startBatchA+1)call ichorQuit('Screening Dim Mismatch1',-1)
   IF(nBatchB.NE.endBatchB-startBatchB+1)call ichorQuit('Screening Dim Mismatch2',-1)
   IF(nBatchAGAB.NE.nBatchA.OR.nBatchBGAB.NE.nBatchB)THEN
      !WARNING BATCHCALC NOT FULL INTEGRAL IS CALCULATED
      allocate(BATCHGAB2(nBatchAGAB,nBatchBGAB))
      call mem_ichor_alloc(BATCHGAB2)
      call RetrieveGabFromIchorSaveGabModule(nBatchAGAB,nBatchBGAB,IchorGabID1,BATCHGAB2)
      call ExtractBatchGabFromFullGab(nBatchA,nBatchB,BATCHGAB,nBatchAGAB,nBatchBGAB,BATCHGAB2,&
           & startBatchA,endBatchA,startBatchB,endBatchB)
      call mem_ichor_dealloc(BATCHGAB2)
      deallocate(BATCHGAB2)
   ELSE
      call RetrieveGabFromIchorSaveGabModule(nBatchA,nBatchB,IchorGabID1,BATCHGAB)
   ENDIF

   allocate(MaxGabForTypeAB(nTypesA,nTypesB))
   call mem_ichor_alloc(MaxGabForTypeAB)
   call ObtainMaxGabForType(MaxGabForTypeAB,nTypesA,nTypesB,nAtomsOfTypeA,&
        & nAtomsOfTypeB,BatchIndexOfTypeA,BatchIndexOfTypeB,BATCHGAB,&
        & nBatchA,nBatchB)
   MaxGabLHS = MAXVAL(MaxGabForTypeAB)
   allocate(BatchIndexOfTypeC(nTypesC))
   call mem_ichor_alloc(BatchIndexOfTypeC)
   allocate(BatchIndexOfTypeD(nTypesD))
   call mem_ichor_alloc(BatchIndexOfTypeD)
   allocate(MaxGabForTypeCD(nTypesC,nTypesD))
   call mem_ichor_alloc(MaxGabForTypeCD)
   IF(IchorGabID1.EQ.IchorGabID2)THEN
      BATCHGCD => BATCHGAB
      MaxGabForTypeCD = MaxGabForTypeAB
      MaxGabRHS = MaxGabLHS
      BatchIndexOfTypeC = BatchIndexOfTypeA
      BatchIndexOfTypeD = BatchIndexOfTypeB
      nBatchC = nBatchA
      nBatchD = nBatchB
   ELSE
      call ConstructBatchIndexOfType(BatchIndexOfTypeC,nTypesC,nAtomsOfTypeC,nBatchC)
      call ConstructBatchIndexOfType(BatchIndexOfTypeD,nTypesD,nAtomsOfTypeD,nBatchD)
      allocate(BATCHGCD(nBatchC,nBatchD))
      call mem_ichor_alloc(BATCHGCD)
      IF(nBatchC.NE.endBatchC-startBatchC+1)call ichorQuit('Screening Dim Mismatch3',-1)
      IF(nBatchD.NE.endBatchD-startBatchD+1)call ichorQuit('Screening Dim Mismatch4',-1)
      call RetrieveGabDimFromIchorSaveGabModule(nBatchCGCD,nBatchDGCD,IchorGabID2)      
      IF(nBatchCGCD.NE.nBatchC.OR.nBatchDGCD.NE.nBatchD)THEN
         !WARNING BATCHCALC NOT FULL INTEGRAL IS CALCULATED
         allocate(BATCHGCD2(nBatchCGCD,nBatchDGCD))
         call mem_ichor_alloc(BATCHGCD2)
         call RetrieveGabFromIchorSaveGabModule(nBatchCGCD,nBatchDGCD,IchorGabID2,BATCHGCD2)
         call ExtractBatchGabFromFullGab(nBatchC,nBatchD,BATCHGCD,nBatchCGCD,nBatchDGCD,&
              & BATCHGCD2,startBatchC,endBatchC,startBatchD,endBatchD)
         call mem_ichor_dealloc(BATCHGCD2)
         deallocate(BATCHGCD2)
      ELSE
         call RetrieveGabFromIchorSaveGabModule(nBatchC,nBatchD,IchorGabID2,BATCHGCD)
      ENDIF
      call ObtainMaxGabForType(MaxGabForTypeCD,nTypesC,nTypesD,nAtomsOfTypeC,&
           & nAtomsOfTypeD,BatchIndexOfTypeC,BatchIndexOfTypeD,BATCHGCD,&
           & nBatchC,nBatchD)
      MaxGabRHS = MAXVAL(MaxGabForTypeCD)
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
   allocate(BATCHGAB(nBatchA,nBatchB))
   call mem_ichor_alloc(BATCHGAB)
   BATCHGCD => BATCHGAB
ENDIF
IF(CSScreen.OR.ODscreen)THEN
   !it is possible that we skip integrals
   call ichorzero2(OutputStorage,OutputDim1*OutputDim2,OutputDim3*OutputDim4*OutputDim5)
ENDIF

call determinePermuteSym(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
IF(INTPRINT.GT.50)write(lupri,*)'SameLHSaos,SameRHSaos,SameODs',SameLHSaos,SameRHSaos,SameODs
IF(DoLink)THEN
   IF((SameLHSaos.OR.SameRHSaos).OR.SameODs)call ichorQuit('Ichor Link Symmetry not enabled yet',-1)
ENDIF
IF(doMOtrans)THEN
   !have to deactivate SameOD at least for now. 
   SameODs = .FALSE.
ENDIF
Spherical = SphericalSpec.EQ.SphericalParam
oldmaxangmomABCD = -25

MaxTotalAngmom = MAXVAL(AngmomOfTypeA) + MAXVAL(AngmomOfTypeB) &
     & + MAXVAL(AngmomOfTypeC) + MAXVAL(AngmomOfTypeD)
call set_GGem(IchorOperatorSpec,MaxTotalAngmom)
!we loop over Total angmom in order to ensure that we first do all
!SSSS integrals then PSSS,SPSS,SSPS,SSSP, ...
!this mean we call GAMMATABULATION a limited number of times and
!we reduce branch misprediction inside the code and reuse instruction cache. 
DO IAngmomTypes = 0,MaxTotalAngmom
 DO ItypeDnon=1,nTypesD
  ! TYPE D CALC ===========================
  call ObtainTypeInfo(nTypesD,ItypeDnon,OrderdListD,nAtomsOfTypeD,AngmomOfTypeD,nPrimOfTypeD,&
       & nContOfTypeD,ExtentOfTypeD,ItypeD,nAtomsD,AngmomD,nPrimD,nContD,extentD,nOrbCompD,&
       & nOrbD,nDimD,nCartOrbCompD,spherical)
  IF(nAtomsD.EQ.0)CYCLE
  allocate(expD(nPrimD))
  call mem_ichor_alloc(expD)
  allocate(ContractCoeffD(nPrimD,nContD))
  call mem_ichor_alloc(ContractCoeffD)
  allocate(Dcenter(3,nAtomsD))
  call mem_ichor_alloc(Dcenter)
  allocate(StartOrbitalD(nAtomsD))
  call mem_ichor_alloc(StartOrbitalD)
  call build_exp_ContractCoeff_center(nPrimD,nContD,nAtomsD,ntypesD,iTypeD,&
       & exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters,startOrbitalOfTypeD,&
       & expD,ContractCoeffD,Dcenter,StartOrbitalD,MaxnAtomsD,MaxnprimD,MaxnContD,lupri)

  IF(doMOtrans)Then
     allocate(CMO4D(nDimD,nCMO4))
     call mem_ichor_alloc(CMO4D)
     call FormActiveCmat(nDimD,nCMO4,CMO4D,nAtomsD,startOrbitalD,nOrbD,&
          & InputM(InputMspec(4))%M,InputM(InputMspec(4))%n1)
     allocate(CMO3D(nDimD,nCMO3))
     call mem_ichor_alloc(CMO3D)
     call FormActiveCmat(nDimD,nCMO3,CMO3D,nAtomsD,startOrbitalD,nOrbD,&
          & InputM(InputMspec(3))%M,InputM(InputMspec(3))%n1)
  ENDIF
  ! DONE TYPE D CALC ===========================
  DO ItypeCnon=1,nTypesC
   ! TYPE C CALC ================================
   ItypeC = OrderdListC(ItypeCnon)  
   IF(SameRHSaos .AND. ItypeD.GT.ItypeC)CYCLE
   TriangularRHSAtomLoop = SameRHSaos .AND. ItypeD.EQ.ItypeC
   PermuteRHSTypes = SameRHSaos .AND. ItypeD.LT.ItypeC
   call ObtainTypeInfo(nTypesC,ItypeCnon,OrderdListC,nAtomsOfTypeC,AngmomOfTypeC,nPrimOfTypeC,&
        & nContOfTypeC,ExtentOfTypeC,ItypeC,nAtomsC,AngmomC,nPrimC,nContC,extentC,nOrbCompC,&
        & nOrbC,nDimC,nCartOrbCompC,spherical)   
   IF(nAtomsC.EQ.0)CYCLE
   allocate(expC(nPrimC))
   call mem_ichor_alloc(expC)
   allocate(ContractCoeffC(nPrimC,nContC))
   call mem_ichor_alloc(ContractCoeffC)
   allocate(Ccenter(3,nAtomsC))
   call mem_ichor_alloc(Ccenter)
   allocate(StartOrbitalC(nAtomsC))
   call mem_ichor_alloc(StartOrbitalC)
   call build_exp_ContractCoeff_center(nPrimC,nContC,nAtomsC,ntypesC,iTypeC,&
        & exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters,startOrbitalOfTypeC,&
        & expC,ContractCoeffC,Ccenter,StartOrbitalC,MaxnAtomsC,MaxnprimC,MaxnContC,lupri)
   IF(doMOtrans)Then
      allocate(CMO3C(nDimC,nCMO3))
      call mem_ichor_alloc(CMO3C)
      call FormActiveCmat(nDimC,nCMO3,CMO3C,nAtomsC,startOrbitalC,nOrbC,&
           & InputM(InputMspec(3))%M,InputM(InputMspec(3))%n1)
      allocate(CMO4C(nDimC,nCMO4))
      call mem_ichor_alloc(CMO4C)
      call FormActiveCmat(nDimC,nCMO4,CMO4C,nAtomsC,startOrbitalC,nOrbC,&
           & InputM(InputMspec(4))%M,InputM(InputMspec(4))%n1)
   ENDIF
   ! DONE TYPE C CALC ===========================
   ! TYPE Q CALC ================================
   sumExtent2CD = (extentC + extentD)*(extentC + extentD)
   
   AngmomQ = AngmomC + AngmomD
   nPrimQ = nPrimC*nPrimD
   nContQ = nContC*nContD
   nOrbCompQ = nOrbCompC*nOrbCompD
   nCartOrbCompQ = nCartOrbCompC*nCartOrbCompD
   nTUVQ = (AngmomQ+1)*(AngmomQ+2)*(AngmomQ+3)/6
   allocate(Qiprim1(nPrimQ))
   call mem_ichor_alloc(Qiprim1) !not used
   allocate(Qiprim2(nPrimQ))
   call mem_ichor_alloc(Qiprim2) !not used
   IF (nContQ.EQ. 1)THEN
      Qsegmented = .TRUE.
   ELSE
      Qsegmented = .FALSE.
   ENDIF
   allocate(expQ(nPrimQ))
   call mem_ichor_alloc(expQ)
   call build_expP(nPrimC,nPrimD,expC,expD,expQ)

   allocate(noScreenCD(nAtomsC,nAtomsD))
   call mem_ichor_alloc(noScreenCD)
   allocate(noScreenCD2(nAtomsC,nAtomsD))
   call mem_ichor_alloc(noScreenCD2)
   IF(ODScreen)THEN
      call ODscreen_noScreen(nAtomsC,nAtomsD,Ccenter,Dcenter,sumextent2CD,noScreenCD)
!      call loutput(noscreenCD,nAtomsC,nAtomsD,6)
   ELSE
      call build_EmptynoScreen1(nAtomsC,nAtomsD,noScreenCD)
   ENDIF
   !LinK: 
   IF(DoLink)THEN !Make Atomic Screening matric AtomGCD
      !maybe this should be used in general - smaller than BATCHGCD ad more straightforward 
      allocate(AtomGCD(nAtomsC,nAtomsD))
      call mem_ichor_alloc(AtomGCD)
      allocate(MaxGCDVec(nAtomsD))
      call mem_ichor_alloc(MaxGCDVec)
      call Build_AtomGAB(nAtomsC,nAtomsD,nBatchC,nBatchD,ntypesC,ntypesD,&
           & ItypeC,ItypeD,BatchIndexOfTypeC,BatchIndexOfTypeD,BATCHGCD,&
           & AtomGCD,MaxGCDVec,MaxAtomGabelmRHS)
      !Make from Screening CD (GCD) a list of contributing C atoms for given D atom)
      !no restriction on atom C for given atom D (unlike A,B)
      allocate(nKetList(nAtomsD))
      allocate(KetList(nAtomsC,nAtomsD))
      call mem_ichor_alloc(nKetList)
      call mem_ichor_alloc(KetList)
      call ConstructBraList(nKetList,KetList,nAtomsC,nAtomsD,THRESHOLD_CS,MaxGabLHS,.FALSE.,&
           & atomGCD)
   ENDIF

   ! DONE TYPE Q CALC ===========================

   DO ItypeBnon=1,nTypesB
    ! TYPE B CALC ================================
    call ObtainTypeInfo(nTypesB,ItypeBnon,OrderdListB,nAtomsOfTypeB,AngmomOfTypeB,nPrimOfTypeB,&
         & nContOfTypeB,ExtentOfTypeB,ItypeB,nAtomsB,AngmomB,nPrimB,nContB,extentB,nOrbCompB,&
         & nOrbB,nDimB,nCartOrbCompB,spherical)
    IF(nAtomsB.EQ.0)CYCLE
    allocate(expB(nPrimB))
    call mem_ichor_alloc(expB)
    allocate(ContractCoeffB(nPrimB,nContB))
    call mem_ichor_alloc(ContractCoeffB)
    allocate(Bcenter(3,nAtomsB))
    call mem_ichor_alloc(Bcenter)
    allocate(StartOrbitalB(nAtomsB))
    call mem_ichor_alloc(StartOrbitalB)
    call build_exp_ContractCoeff_center(nPrimB,nContB,nAtomsB,ntypesB,iTypeB,&
         & exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters,startOrbitalOfTypeB,&
         & expB,ContractCoeffB,Bcenter,StartOrbitalB,MaxnAtomsB,MaxnprimB,MaxnContB,lupri)
    IF(doMOtrans)Then
       allocate(CMO2B(nDimB,nCMO2))
       call mem_ichor_alloc(CMO2B)
       call FormActiveCmat(nDimB,nCMO2,CMO2B,nAtomsB,startOrbitalB,nOrbB,&
            & InputM(InputMspec(2))%M,InputM(InputMspec(2))%n1)
       allocate(CMO1B(nDimB,nCMO1))
       call mem_ichor_alloc(CMO1B)
       call FormActiveCmat(nDimB,nCMO1,CMO1B,nAtomsB,startOrbitalB,nOrbB,&
            & InputM(InputMspec(1))%M,InputM(InputMspec(1))%n1)
    ENDIF
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
     call ObtainTypeInfo(nTypesA,ItypeAnon,OrderdListA,nAtomsOfTypeA,AngmomOfTypeA,nPrimOfTypeA,&
          & nContOfTypeA,ExtentOfTypeA,ItypeA,nAtomsA,AngmomA,nPrimA,nContA,extentA,nOrbCompA,&
          & nOrbA,nDimA,nCartOrbCompA,spherical)
     IF(nAtomsA.EQ.0)CYCLE
     TotalAngmom = AngmomA + AngmomB + AngmomQ
     IF(TotalAngmom.EQ.IAngmomTypes)THEN
      !This if statement ensures that we call GAMMATABULATION a very limited number of times
      !and it reduceses branch mispredictions inside the code,...
      allocate(expA(nPrimA))
      call mem_ichor_alloc(expA)
      allocate(ContractCoeffA(nPrimA,nContA))
      call mem_ichor_alloc(ContractCoeffA)
      allocate(Acenter(3,nAtomsA))
      call mem_ichor_alloc(Acenter)
      allocate(StartOrbitalA(nAtomsA))
      call mem_ichor_alloc(StartOrbitalA)
      call build_exp_ContractCoeff_center(nPrimA,nContA,nAtomsA,ntypesA,iTypeA,&
           & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,startOrbitalOfTypeA,&
           & expA,ContractCoeffA,Acenter,StartOrbitalA,MaxnAtomsA,MaxnprimA,MaxnContA,lupri)

      IF(doMOtrans)Then
         allocate(CMO1A(nDimA,nCMO1))
         call mem_ichor_alloc(CMO1A)
         call FormActiveCmat(nDimA,nCMO1,CMO1A,nAtomsA,startOrbitalA,nOrbA,&
              & InputM(InputMspec(1))%M,InputM(InputMspec(1))%n1)
         allocate(CMO2A(nDimA,nCMO2))
         call mem_ichor_alloc(CMO2A)
         call FormActiveCmat(nDimA,nCMO2,CMO2A,nAtomsA,startOrbitalA,nOrbA,&
              & InputM(InputMspec(2))%M,InputM(InputMspec(2))%n1)
      ENDIF

      ! DONE TYPE A CALC ===========================
      ! TYPE P CALC ================================
      AngmomP = AngmomA + AngmomB
      sumExtent2AB = (extentA + extentB)*(extentA + extentB)
      nPrimP = nPrimA*nPrimB
      nContP = nContA*nContB
      nOrbCompP = nOrbCompA*nOrbCompB
      nCartOrbCompP = nCartOrbCompA*nCartOrbCompB
      nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
      nTUV = (TotalAngmom+1)*(TotalAngmom+2)*(TotalAngmom+3)/6
      IF (nContP.EQ. 1)THEN
         Psegmented = .TRUE.
      ELSE
         Psegmented = .FALSE.
      ENDIF
      allocate(expP(nPrimP))
      call mem_ichor_alloc(expP)
      allocate(inversexpP(nPrimP))
      call mem_ichor_alloc(inversexpP)
      call build_expQ_inverseexpQ(nPrimA,nPrimB,expA,expB,expP,inversexpP)

      !LinK: 
      IF(DoLink)THEN !Make Atomic Screening matric AtomGAB
         !maybe this should be used in general - smaller than BATCHGCD ad more straightforward 
         allocate(AtomGAB(nAtomsA,nAtomsB))
         call mem_ichor_alloc(AtomGAB)
         allocate(MaxGABVec(nAtomsB))
         call mem_ichor_alloc(MaxGABVec)
         call Build_AtomGAB(nAtomsA,nAtomsB,nBatchA,nBatchB,ntypesA,ntypesB,&
              & ItypeA,ItypeB,BatchIndexOfTypeA,BatchIndexOfTypeB,BATCHGAB,&
              & AtomGAB,MaxGABVec,MaxAtomGabelmLHS)
      ENDIF
      ! DONE TYPE P CALC ===========================
      
      ! TYPE PQ CALC ================================
      TotalAngmom = AngmomA + AngmomB + AngmomQ
      maxangmomABCD = AngmomA + AngmomB + AngmomQ
      UseGeneralCode = .FALSE. !Use Specialized code when appropriate. 
      IF(GGemOperatorCalc)UseGeneralCode = .TRUE.
      IF(MAX(AngmomA,AngmomB,AngmomC,AngmomD).GT.MaxSpecialAngmom)UseGeneralCode = .TRUE.
      IF(maxangmomABCD.NE.oldmaxangmomABCD)THEN
       IF(oldmaxangmomABCD.NE.-25)THEN
          call mem_ichor_dealloc(TABFJW)
          deallocate(TABFJW)
       ENDIF
       nTABFJW1 = AngmomA + AngmomB + AngmomQ + 3 
       !only need + 3 after Branos change in BUILD_RJ000 
       nTABFJW2 = 1200
       !TABFJW(0:nTABFJW1,0:nTABFJW2)
       allocate(TABFJW(0:nTABFJW1,0:nTABFJW2))
       call mem_ichor_alloc(TABFJW)
       CALL GAMMATABULATION(lupri,maxangmomABCD,nTABFJW1,nTABFJW2,TABFJW)  
       oldmaxangmomABCD = maxangmomABCD
      ENDIF
      allocate(reducedExponents(nPrimP*nPrimQ))
      call mem_ichor_alloc(reducedExponents)
      allocate(integralPrefactor(nPrimP*nPrimQ))
      call mem_ichor_alloc(integralPrefactor)
      call build_reducedExponents_integralPrefactorQP(nPrimP,nPrimQ,expQ,expP,&
           & reducedExponents,integralPrefactor)
      ! DONE TYPE PQ CALC ===========================
      
      NOTDoSSSS = .NOT.(TotalAngmom.EQ.0.AND.(Psegmented.AND.Qsegmented))
      IF(DoLink) NOTDoSSSS=.TRUE.
      IF(DoMoTrans) NOTDoSSSS=.TRUE.
      IF(UseGeneralCode)NOTDoSSSS = .TRUE.
      IF(.NOT.UseCPU)NOTDoSSSS = .TRUE.
      IF(NOTDoSSSS)THEN
         !Determine Sizes of TmpArrays and MaxPasses
         IF(UseCPU)THEN
            call ICI_CPU_OBS_general_size(TMParray1maxsize,&
                 & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
                 & nContA,nContB,nContC,nContD,&
                 & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,nContP,&
                 & nContQ,nPrimQ*nPrimP,nContQ*nContP,Psegmented,Qsegmented)
         ELSE !use GPU code
            call ICI_GPU_OBS_general_size(TMParray1maxsize,&
                 & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
                 & nContA,nContB,nContC,nContD,&
                 & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,nContP,&
                 & nContQ,nPrimQ*nPrimP,nContQ*nContP,Psegmented,Qsegmented)
         ENDIF
         nLocalInt = nOrbA*nOrbB*nOrbC*nOrbD
      ENDIF
      IF(NOTDoSSSS)THEN !FIXME A MESS
         IF(DoLink)THEN         
            allocate(QpreExpFac(nPrimQ))
            call mem_ichor_alloc(QpreExpFac)
            allocate(Qcent(3*nPrimQ))
            call mem_ichor_alloc(Qcent)
         ELSEIF(doMOtrans)THEN
            allocate(QpreExpFac(nPrimQ))
            call mem_ichor_alloc(QpreExpFac)
            allocate(Qcent(3*nPrimQ))
            call mem_ichor_alloc(Qcent)
         ELSEIF(UseCPU)THEN
            allocate(QpreExpFac(nPrimQ))
            call mem_ichor_alloc(QpreExpFac)
            allocate(Qcent(3*nPrimQ))
            call mem_ichor_alloc(Qcent)
         ELSE
            !not needed 
         ENDIF
      ELSE
         allocate(QpreExpFac(nPrimQ))
         call mem_ichor_alloc(QpreExpFac)
         allocate(Qcent(3*nPrimQ))
         call mem_ichor_alloc(Qcent)
      ENDIF
      !calc       
      IF (INTPRINT .GE. 10) THEN
         call PrintTypeExpInfo(nPrimP,nPrimQ,reducedExponents,integralPrefactor,lupri)
         call PrintTypeInfo(AngmomA,AngmomB,AngmomC,AngmomD,nPrimA,nPrimB,nPrimC,nPrimD,&
              & nContA,nContB,nContC,nContD,expA,ContractCoeffA,expB,ContractCoeffB,&
              & expC,ContractCoeffC,expD,ContractCoeffD,&
              & nAtomsA,nAtomsB,nAtomsC,nAtomsD,Acenter,Bcenter,Ccenter,Dcenter,lupri)
      ENDIF
      allocate(noScreenAB(nAtomsA,nAtomsB))
      call mem_ichor_alloc(noScreenAB)
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
      !Can be used for screening? 
      !It should be possible to make sort of a 3D volume - Determining 
      !For max(Rcd) = 0,max(Rpq) = 0 which max(Rab) will make the Integral vanish?  
      !For max(Rpq) = 0,max(Rab) = 0 which max(Rcd) will make the Integral vanish?  
      !For max(Rab) = 0,max(Rcd) = 0 which max(Rpq) will make the Integral vanish?  
      !Different screening for different types    (SS|SS) special 
      !APE type screening? 
      !For SSSS it should be exact to do
      !exp(-mu_AB R_AB^2)exp(-mu_CD R_CD^2) R000 and use
      !R000(alpha*R_pq^2) = F0(alpha*R_pq^2) = sqrt(pi/(4*x))erf(x) ;x=alpha*R_pq^2
      !R000(alpha*R_pq^2) > sqrt(pi/(4*x)) ;x=alpha*R_pq^2
      !look at possible things for SSSP type things.... 
      
      !LinK: 
      IF(DoLink)THEN !Form active Dmat: DmatBD,DmatAD,DmatBC,DmatAC
         allocate(DmatBD(nDimB,nDimD,IchorInputDim3))
         allocate(DmatAD(nDimA,nDimD,IchorInputDim3))
         allocate(DmatBC(nDimB,nDimC,IchorInputDim3))
         allocate(DmatAC(nDimA,nDimC,IchorInputDim3))
         call mem_ichor_alloc(DmatBD)
         call mem_ichor_alloc(DmatAD)
         call mem_ichor_alloc(DmatBC)
         call mem_ichor_alloc(DmatAC)
         Call FormActiveDmat(nDimB,nDimD,DmatBD,nAtomsB,nAtomsD,startOrbitalB,startOrbitalD,&
              & nOrbB,nOrbD,InputStorage,IchorInputDim1,IchorInputDim2,IchorInputDim3)
!         Call FormActiveDmat(nDimA,nDimD,DmatAD,nAtomsA,nAtomsD,startOrbitalA,startOrbitalD,&
!              & nOrbA,nOrbD,InputStorage,IchorInputDim1,IchorInputDim2,IchorInputDim3)
!         Call FormActiveDmat(nDimB,nDimC,DmatBC,nAtomsB,nAtomsC,startOrbitalB,startOrbitalC,&
!              & nOrbB,nOrbC,InputStorage,IchorInputDim1,IchorInputDim2,IchorInputDim3)
!         Call FormActiveDmat(nDimA,nDimC,DmatAC,nAtomsA,nAtomsC,startOrbitalA,startOrbitalC,&
!              & nOrbA,nOrbC,InputStorage,IchorInputDim1,IchorInputDim2,IchorInputDim3)
         allocate(KmatAC(nDimA,nDimC,IchorInputDim3))
         allocate(KmatBC(nDimB,nDimC,IchorInputDim3))
         allocate(KmatAD(nDimA,nDimD,IchorInputDim3))
         allocate(KmatBD(nDimB,nDimD,IchorInputDim3))
         call mem_ichor_alloc(KmatAC)
         call mem_ichor_alloc(KmatBC)
         call mem_ichor_alloc(KmatAD)
         call mem_ichor_alloc(KmatBD)
         call ichorzero(KmatAC,SIZE(KmatAC))
         call ichorzero(KmatBC,SIZE(KmatBC))
         call ichorzero(KmatAD,SIZE(KmatAD))
         call ichorzero(KmatBD,SIZE(KmatBD))
         allocate(ReducedDmatBD(nAtomsB,nAtomsD))
         allocate(ReducedDmatBC(nAtomsB,nAtomsC))
         call mem_ichor_alloc(ReducedDmatBD)
         call mem_ichor_alloc(ReducedDmatBC)
         call FormReducedDmat(DmatBD,ReducedDmatBD,nOrbB,nOrbD,nAtomsB,nAtomsD,IchorInputDim3)
!         call FormReducedDmat(DmatBC,ReducedDmatBC,nOrbB,nOrbC,nAtomsB,nAtomsC,IchorInputDim3)

         !Make from Screening AB (GAB) a list of contributing A atoms for given B atom)
         !Including permutation/symmetry restrictions on atom A,B
         allocate(nBraList(nAtomsB))
         allocate(BraList(nAtomsA,nAtomsB))
         call mem_ichor_alloc(nBraList)
         call mem_ichor_alloc(BraList)
         !note we here use MaxAtomGabelmRHS instead of MaxGabRHS because we are constructing
         !this BraList for each type of CD - and can use the more accurate MaxAtomGabelmRHS
         call ConstructBraList(nBraList,BraList,nAtomsA,nAtomsB,THRESHOLD_CS,MaxAtomGabelmRHS,&
              & TriangularLHSAtomLoop,atomGAB)
            
         !Make from Screening info a list of contributing B atoms for given D atom)
         allocate(nBraketList(nAtomsD))
         allocate(BraketList(nAtomsB,nAtomsD))
         call mem_ichor_alloc(nBraketList)
         call mem_ichor_alloc(BraketList)
         call ConstructBraketList(nAtomsB,nAtomsD,nBraketList,BraketList,TriangularODAtomLoop,&
              & ReducedDmatBD,MaxGABVec,MaxGCDVec,THRESHOLD_CS)
      ENDIF

      allocate(PpreExpFacPass(nPrimP,nAtomsA*nAtomsB))
      call mem_ichor_alloc(PpreExpFacPass)
      allocate(PcentPass(3*nPrimP,nAtomsA*nAtomsB))
      call mem_ichor_alloc(PcentPass)
      IF(CSScreen)THEN
         iBatchIndexOfTypeD = BatchIndexOfTypeD(ItypeD)
         iBatchIndexOfTypeC = BatchIndexOfTypeC(ItypeC)
      ENDIF
!      call IchorTimer('START',TSTART,TEND,LUPRI)
      IF(NOTDoSSSS)THEN
         allocate(Pdistance12Pass(3,nAtomsA*nAtomsB))
         call mem_ichor_alloc(Pdistance12Pass)
         CALL Build_pcent_Pdistance12_PpreExpFac(nPrimA,nPrimB,natomsA,natomsB,nContA,nContB,&
              & inversexpP,expA,expB,Acenter,Bcenter,ContractCoeffA,ContractCoeffB,PSegmented,&
              & PcentPass,Pdistance12Pass,PpreExpFacPass,INTPRINT)
         allocate(PpreExpFac(nPrimP))
         call mem_ichor_alloc(PpreExpFac)

         IF(DoLink)THEN         
            !Link type Loop over the atoms of this type 
            !The integrals are contracted with densities to yield exchange matrices
            call IchorTypeLinKLoop(nAtomsA,nPrimA,nContA,nOrbCompA,&
                 & expA,ContractCoeffA,AngmomA,Acenter,nOrbA,&
                 & nAtomsB,nPrimB,nContB,nOrbCompB,&
                 & expB,ContractCoeffB,AngmomB,Bcenter,nOrbB,&
                 & nAtomsC,nPrimC,nContC,nOrbCompC,&
                 & expC,ContractCoeffC,AngmomC,Ccenter,nOrbC,&
                 & nAtomsD,nPrimD,nContD,nOrbCompD,&
                 & expD,ContractCoeffD,AngmomD,Dcenter,nOrbD,&
                 & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
                 & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
                 & nPrimP,nPrimQ,nContP,nContQ,expP,expQ,&
                 & qcent,Ppreexpfac,Qpreexpfac,&
                 & Qiprim1,Qiprim2,&
                 & nTABFJW1,nTABFJW2,TABFJW,TotalAngmom,&
                 & CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
                 & Qsegmented,Psegmented,TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
                 & reducedExponents,integralPrefactor,&
                 & PcentPass,Pdistance12Pass,PpreExpFacPass,&
                 & Qdistance12,PQorder,Spherical,TMParray1maxsize,nLocalInt,&
                 & TMParray2maxsize,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri,&
                 & DmatBD,DmatAD,DmatBC,&
                 & DmatAC,KmatBD,KmatAD,KmatBC,KmatAC,nDimA,nDimB,nDimC,nDImD,IchorInputDim3,&
                 & ReducedDmatBD,ReducedDmatBC,AtomGAB,AtomGCD,&
                 & nKetList,KetList,nBraList,BraList,nBraketList,BraketList,&
                 & SameRHSaos,SameLHSaos,SameODs)

            Call addActiveKmatToOutput(nDimA,nDimB,nDimC,nDimD,nAtomsA,nAtomsB,nAtomsC,nAtomsD,&
                 & nOrbA,nOrbB,nOrbC,nOrbD,startOrbitalA,startOrbitalB,startOrbitalC,startOrbitalD,&
                 & OutputDim1,OutputDim2,OutputDim5,OutputStorage,KmatAC)
           
         ELSEIF(doMOtrans)THEN

            call IchorTypeMOtransLoop(nAtomsA,nPrimA,nContA,nOrbCompA,startOrbitalA,&
                 & iBatchIndexOfTypeA,expA,ContractCoeffA,AngmomA,Acenter,nBatchA,nOrbA,&
                 & nAtomsB,nPrimB,nContB,nOrbCompB,startOrbitalB,&
                 & iBatchIndexOfTypeB,expB,ContractCoeffB,AngmomB,Bcenter,nBatchB,nOrbB,&
                 & nAtomsC,nPrimC,nContC,nOrbCompC,startOrbitalC,&
                 & iBatchIndexOfTypeC,expC,ContractCoeffC,AngmomC,Ccenter,nBatchC,nOrbC,&
                 & nAtomsD,nPrimD,nContD,nOrbCompD,startOrbitalD,&
                 & iBatchIndexOfTypeD,expD,ContractCoeffD,AngmomD,Dcenter,nBatchD,nOrbD,&
                 & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
                 & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
                 & nPrimP,nPrimQ,nContP,nContQ,expP,expQ,&
                 & qcent,Ppreexpfac,Qpreexpfac,&
                 & Qiprim1,Qiprim2,&
                 & nTABFJW1,nTABFJW2,TABFJW,TotalAngmom,&
                 & CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
                 & Qsegmented,Psegmented,TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
                 & reducedExponents,integralPrefactor,&
                 & PcentPass,Pdistance12Pass,PpreExpFacPass,&
                 & Qdistance12,PQorder,&
                 & BATCHGCD,BATCHGAB,Spherical,TMParray1maxsize,nLocalInt,&
                 & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,OutputStorage,&
                 & TMParray2maxsize,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri,&
                 & nCMO1,nCMO2,nCMO3,nCMO4,CMO1A,CMO2B,CMO3C,CMO4D,&
                 & CMO1B,CMO2A,CMO3D,CMO4C,PermuteRHSTypes)
         ELSEIF(UseCPU)THEN
            call IchorTypeIntegralLoopCPU(nAtomsA,nPrimA,nContA,nOrbCompA,startOrbitalA,&
                 & iBatchIndexOfTypeA,expA,ContractCoeffA,AngmomA,Acenter,nBatchA,nOrbA,&
                 & nAtomsB,nPrimB,nContB,nOrbCompB,startOrbitalB,&
                 & iBatchIndexOfTypeB,expB,ContractCoeffB,AngmomB,Bcenter,nBatchB,nOrbB,&
                 & nAtomsC,nPrimC,nContC,nOrbCompC,startOrbitalC,&
                 & iBatchIndexOfTypeC,expC,ContractCoeffC,AngmomC,Ccenter,nBatchC,nOrbC,&
                 & nAtomsD,nPrimD,nContD,nOrbCompD,startOrbitalD,&
                 & iBatchIndexOfTypeD,expD,ContractCoeffD,AngmomD,Dcenter,nBatchD,nOrbD,&
                 & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
                 & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
                 & nPrimP,nPrimQ,nContP,nContQ,expP,expQ,&
                 & qcent,Ppreexpfac,Qpreexpfac,&
                 & Qiprim1,Qiprim2,&
                 & nTABFJW1,nTABFJW2,TABFJW,TotalAngmom,&
                 & CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
                 & Qsegmented,Psegmented,TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
                 & reducedExponents,integralPrefactor,&
                 & PcentPass,Pdistance12Pass,PpreExpFacPass,&
                 & Qdistance12,PQorder,&
                 & BATCHGCD,BATCHGAB,Spherical,TMParray1maxsize,nLocalInt,&
                 & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,OutputStorage,&
                 & TMParray2maxsize,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri)
         ELSE !GPU
            call IchorTypeIntegralLoopGPU(nAtomsA,nPrimA,nContA,nOrbCompA,startOrbitalA,&
                 & iBatchIndexOfTypeA,expA,ContractCoeffA,AngmomA,Acenter,nBatchA,nOrbA,&
                 & nAtomsB,nPrimB,nContB,nOrbCompB,startOrbitalB,&
                 & iBatchIndexOfTypeB,expB,ContractCoeffB,AngmomB,Bcenter,nBatchB,nOrbB,&
                 & nAtomsC,nPrimC,nContC,nOrbCompC,startOrbitalC,&
                 & iBatchIndexOfTypeC,expC,ContractCoeffC,AngmomC,Ccenter,nBatchC,nOrbC,&
                 & nAtomsD,nPrimD,nContD,nOrbCompD,startOrbitalD,&
                 & iBatchIndexOfTypeD,expD,ContractCoeffD,AngmomD,Dcenter,nBatchD,nOrbD,&
                 & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
                 & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
                 & nPrimP,nPrimQ,nContP,nContQ,expP,expQ,&
                 & Ppreexpfac,Qiprim1,Qiprim2,&
                 & nTABFJW1,nTABFJW2,TABFJW,TotalAngmom,&
                 & CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
                 & Qsegmented,Psegmented,TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
                 & reducedExponents,integralPrefactor,PcentPass,Pdistance12Pass,PpreExpFacPass,&
                 & PQorder,BATCHGCD,BATCHGAB,Spherical,TMParray1maxsize,nLocalInt,&
                 & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,OutputStorage,&
                 & TMParray2maxsize,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri,nStaticParamIfac)
         ENDIF
!         call mem_ichor_dealloc(pcent)
!         deallocate(pcent)
         call mem_ichor_dealloc(PpreExpFac)
         deallocate(PpreExpFac)
         call mem_ichor_dealloc(Pdistance12Pass)
         deallocate(Pdistance12Pass)
      ELSE
         CALL Build_pcent_PpreExpFac(nPrimA,nPrimB,natomsA,natomsB,nContA,nContB,&
              & inversexpP,expA,expB,Acenter,Bcenter,ContractCoeffA,ContractCoeffB,PSegmented,&
              & PcentPass,PpreExpFacPass,INTPRINT)
         allocate(LocalIntPass(nOrbA*nAtomsA*nOrbB*nAtomsB,nOrbC*nOrbD))
         CALL MEM_ICHOR_ALLOC(LocalIntPass)
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
              & PcentPass,PpreExpFacPass,BATCHGCD,BATCHGAB,&
              & LocalIntPass,OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
              & OutputStorage,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri)
         CALL MEM_ICHOR_deALLOC(LocalIntPass)
         deallocate(LocalIntPass)
      ENDIF
!      call build_TYPESTRING(TYPESTRING,AngmomA,AngmomB,AngmomC,AngmomD,&
!           & nTypesA,nTypesB,nTypesC,nTypesD,ItypeA,ItypeB,ItypeC,ItypeD)
!      call IchorTimer(TYPESTRING,TSTART,TEND,LUPRI)

      IF(NOTDoSSSS)THEN !FIXME A MESS
         IF(DoLink)THEN         
            call mem_ichor_dealloc(Qcent)
            deallocate(Qcent)
            call mem_ichor_dealloc(QpreExpFac)
            deallocate(QpreExpFac)
         ELSEIF(doMOtrans)THEN
            call mem_ichor_dealloc(Qcent)
            deallocate(Qcent)
            call mem_ichor_dealloc(QpreExpFac)
            deallocate(QpreExpFac)
         ELSEIF(UseCPU)THEN
            call mem_ichor_dealloc(Qcent)
            deallocate(Qcent)
            call mem_ichor_dealloc(QpreExpFac)
            deallocate(QpreExpFac)
         ELSE
            !not needed 
         ENDIF
      ELSE
         call mem_ichor_dealloc(Qcent)
         deallocate(Qcent)
         call mem_ichor_dealloc(QpreExpFac)
         deallocate(QpreExpFac)
      ENDIF

      call mem_ichor_dealloc(PpreExpFacPass)
      deallocate(PpreExpFacPass)
      call mem_ichor_dealloc(PcentPass)
      deallocate(PcentPass)
      IF(PermuteODTypes)THEN
         !place (ABCD) in (CDAB) and possible (DCAB),(CDBA),(DCBA)
         call PermuteODtypesSub(OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputStorage,&
              & nAtomsA,nAtomsB,startOrbitalA,startOrbitalB,nOrbA,nOrbB,&
              & nAtomsC,nAtomsD,startOrbitalC,startOrbitalD,nOrbC,nOrbD,&
              & TriangularODAtomLoop,TriangularRHSAtomLoop,TriangularLHSAtomLoop,&
              & SameLHSaos,SameRHSaos,lupri)
      ENDIF
      IF(DoLink)THEN         
         call mem_ichor_dealloc(DmatBD)
         call mem_ichor_dealloc(DmatAD)
         call mem_ichor_dealloc(DmatBC)
         call mem_ichor_dealloc(DmatAC)
         call mem_ichor_dealloc(KmatBD)
         call mem_ichor_dealloc(KmatAD)
         call mem_ichor_dealloc(KmatBC)
         call mem_ichor_dealloc(KmatAC)
         call mem_ichor_dealloc(ReducedDmatBD)
         call mem_ichor_dealloc(ReducedDmatBC)
         call mem_ichor_dealloc(nBraList)
         call mem_ichor_dealloc(nBraketList)
         call mem_ichor_dealloc(BraketList)
         call mem_ichor_dealloc(BraList)
         deallocate(DmatBD)
         deallocate(DmatAD)
         deallocate(DmatBC)
         deallocate(DmatAC)
         deallocate(KmatBD)
         deallocate(KmatAD)
         deallocate(KmatBC)
         deallocate(KmatAC)
         deallocate(ReducedDmatBD)
         deallocate(ReducedDmatBC)
         deallocate(nBraList)
         deallocate(BraList)
         deallocate(nBraketList)
         deallocate(BraketList)
      ENDIF
      call mem_ichor_dealloc(reducedExponents)
      deallocate(reducedExponents)
      call mem_ichor_dealloc(integralPrefactor)
      deallocate(integralPrefactor)
      call mem_ichor_dealloc(noScreenAB)
      deallocate(noScreenAB)
      call mem_ichor_dealloc(expP)
      deallocate(expP)
      call mem_ichor_dealloc(inversexpP)    
      deallocate(inversexpP)    
      call mem_ichor_dealloc(expA)
      deallocate(expA)
      call mem_ichor_dealloc(ContractCoeffA)
      deallocate(ContractCoeffA)
      call mem_ichor_dealloc(Acenter)
      deallocate(Acenter)
      call mem_ichor_dealloc(StartOrbitalA)
      deallocate(StartOrbitalA)

      IF(doMOtrans)Then
         call mem_ichor_dealloc(CMO1A)
         deallocate(CMO1A)
         call mem_ichor_dealloc(CMO2A)
         deallocate(CMO2A)
      ENDIF

      IF(DoLink)THEN
         call mem_ichor_dealloc(AtomGAB)
         deallocate(AtomGAB)
         call mem_ichor_dealloc(MaxGABvec)
         deallocate(MaxGABvec)
      ENDIF
     ENDIF !correct angmom type
    ENDDO !typeA
    call mem_ichor_dealloc(expB)
    deallocate(expB)
    call mem_ichor_dealloc(ContractCoeffB)
    deallocate(ContractCoeffB)
    call mem_ichor_dealloc(Bcenter)
    deallocate(Bcenter)
    call mem_ichor_dealloc(StartOrbitalB)
    deallocate(StartOrbitalB)
    IF(doMOtrans)Then
       call mem_ichor_dealloc(CMO2B)
       deallocate(CMO2B)
       call mem_ichor_dealloc(CMO1B)
       deallocate(CMO1B)
    ENDIF
   ENDDO !typeB
   IF(PermuteRHSTypes.AND.(.NOT.doMOtrans))THEN
      !place (ABCD) in (ABDC)
      Call PermuteRHStypesSub(OutputDim1*OutputDim2,OutputDim3,OutputDim4,OutputStorage,&
           & nAtomsC,nAtomsD,startOrbitalC,startOrbitalD,nOrbC,nOrbD,&
           & TriangularRHSAtomLoop,lupri)
      !In case of doMOtrans the PermuteRHSTypes have already been performed
   ENDIF
   IF(DoLink)THEN
      call mem_ichor_dealloc(AtomGCD)
      deallocate(AtomGCD)
      call mem_ichor_dealloc(MaxGCDvec)
      deallocate(MaxGCDvec)
      call mem_ichor_dealloc(nKetList)
      deallocate(nKetList)
      call mem_ichor_dealloc(KetList)
      deallocate(KetList)
   ENDIF
   call mem_ichor_dealloc(noScreenCD)
   deallocate(noScreenCD)
   call mem_ichor_dealloc(noScreenCD2)
   deallocate(noScreenCD2)
   call mem_ichor_dealloc(expQ)
   deallocate(expQ)
   call mem_ichor_dealloc(Qiprim1) 
   deallocate(Qiprim1) 
   call mem_ichor_dealloc(Qiprim2) 
   deallocate(Qiprim2) 
   
   call mem_ichor_dealloc(expC)
   deallocate(expC)
   call mem_ichor_dealloc(ContractCoeffC)
   deallocate(ContractCoeffC)
   call mem_ichor_dealloc(Ccenter)
   deallocate(Ccenter)
   call mem_ichor_dealloc(StartOrbitalC)
   deallocate(StartOrbitalC)   

   IF(doMOtrans)Then
      call mem_ichor_dealloc(CMO3C)
      deallocate(CMO3C)
      call mem_ichor_dealloc(CMO4C)
      deallocate(CMO4C)
   ENDIF
  ENDDO !typeC
  call mem_ichor_dealloc(expD)
  deallocate(expD)
  call mem_ichor_dealloc(ContractCoeffD)
  deallocate(ContractCoeffD)
  call mem_ichor_dealloc(Dcenter)
  deallocate(Dcenter)
  call mem_ichor_dealloc(StartOrbitalD)
  deallocate(StartOrbitalD)
  IF(doMOtrans)Then
     call mem_ichor_dealloc(CMO4D)
     deallocate(CMO4D)
     call mem_ichor_dealloc(CMO3D)
     deallocate(CMO3D)
  ENDIF
 ENDDO !typeD
ENDDO 
call mem_ichor_dealloc(TABFJW)
deallocate(TABFJW)

call mem_ichor_dealloc(BATCHGAB)
deallocate(BATCHGAB)
IF(CSScreen)THEN
   call mem_ichor_dealloc(BatchIndexOfTypeA)
   deallocate(BatchIndexOfTypeA)
   call mem_ichor_dealloc(BatchIndexOfTypeB)
   deallocate(BatchIndexOfTypeB)
   call mem_ichor_dealloc(MaxGabForTypeAB)
   deallocate(MaxGabForTypeAB)
   call mem_ichor_dealloc(BatchIndexOfTypeC)
   deallocate(BatchIndexOfTypeC)
   call mem_ichor_dealloc(BatchIndexOfTypeD)
   deallocate(BatchIndexOfTypeD)
   call mem_ichor_dealloc(MaxGabForTypeCD)
   deallocate(MaxGabForTypeCD)
   IF(IchorGabID1.EQ.IchorGabID2)THEN
      !nothing
   ELSE
      call mem_ichor_dealloc(BATCHGCD)
      deallocate(BATCHGCD)
   ENDIF
ENDIF
call mem_ichor_dealloc(OrderdListA)
deallocate(OrderdListA)
call mem_ichor_dealloc(OrderdListB)
deallocate(OrderdListB)
call mem_ichor_dealloc(OrderdListC)
deallocate(OrderdListC)
call mem_ichor_dealloc(OrderdListD)
deallocate(OrderdListD)
call free_GGem()

call retrieve_ichor_memvar(MaxMemAllocated,MemAllocated)
IF(INTPRINT.GT.3)THEN
   call stats_ichor_mem(lupri)
ENDIF
IF(doMOtrans)THEN
   call stats_ichor_mem(lupri)
ENDIF
IF(MemAllocated.NE.0)THEN
   call stats_ichor_mem(lupri)
   call ichorquit('MemoryLeak in IchorEri',lupri)
ENDIF
#ifdef VAR_OPENACC
!$ACC END DATA
#endif

end subroutine IchorEri

subroutine IchorEriMem(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
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
     & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorOperatorSpec,&
     & IchorInputDim1,IchorInputDim2,IchorInputDim3,&
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
Integer,intent(in) :: nContOfTypeB(ntypesB),nPrimOfTypeB(ntypesB)
Integer,intent(in) :: startOrbitalOfTypeB(MaxNatomsB,ntypesB)
Real(realk),intent(in) :: Bcenters(3,MaxNatomsB,ntypesB),exponentsOfTypeB(MaxnprimB,ntypesB)
Real(realk),intent(in) :: ContractCoeffOfTypeB(MaxnprimB,MaxnContB,ntypesB)
!
! Same for Center C
!
integer,intent(in) :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC,startBatchC,endBatchC
Integer,intent(in) :: AngmomOfTypeC(ntypesC),nAtomsOfTypeC(ntypesC)
Integer,intent(in) :: nContOfTypeC(ntypesC),nPrimOfTypeC(ntypesC)
Integer,intent(in) :: startOrbitalOfTypeC(MaxNatomsC,ntypesC)
Real(realk),intent(in) :: Ccenters(3,MaxNatomsC,ntypesC),exponentsOfTypeC(MaxnprimC,ntypesC)
Real(realk),intent(in) :: ContractCoeffOfTypeC(MaxnprimC,MaxnContC,ntypesC)
!
! Same for Center D
!
integer,intent(in) :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD,startBatchD,endBatchD
Integer,intent(in) :: AngmomOfTypeD(ntypesD),nAtomsOfTypeD(ntypesD)
Integer,intent(in) :: nContOfTypeD(ntypesD),nPrimOfTypeD(ntypesD)
Integer,intent(in) :: startOrbitalOfTypeD(MaxNatomsD,ntypesD)
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
!> Operator Specification (IchorOperatorSpec = CoulombOperator = 1) means Coulomb operator
Integer,intent(in) :: IchorOperatorSpec
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
real(realk),intent(inout)::OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5)
!> Logical unit number of output file.
Integer,intent(in) :: lupri
! Local variables
integer :: nPrimP,nContP,nPrimQ,nContQ
integer :: nTABFJW1,nTABFJW2,i1,i2,i3,i4,AngmomQ,TotalAngmom
integer :: i12,offset,K,I,iPrimQ,iPrimP,icont,AngmomP
integer :: oldmaxangmomABCD
integer :: ItypeA,ItypeB,itypeC,itypeD,AngmomA,AngmomB,AngmomC,AngmomD
integer :: ItypeAnon,ItypeBnon,itypeCnon,itypeDnon,nLocalInt
integer :: nPrimA,nPrimB,nContA,nAtomsA,nAtomsB,nAtomsC,nAtomsD,nOrbA
integer :: nDimA,nOrbCompA,nContB,nOrbB,nDimB
integer :: nPrimC,nContC,nPrimD,nContD,nOrbC,nOrbD,nOrbCompB,nOrbCompC,nDimC
integer :: nDimD,nOrbCompD,INTPRINT,maxangmomABCD
integer :: nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD
integer :: nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV
logical :: Psegmented,Qsegmented,PQorder,Spherical,SameLHSaos,SameRHSaos,SameODs
logical :: TriangularLHSAtomLoop,TriangularRHSAtomLoop,PermuteRHS,CSScreen
logical :: TriangularODAtomLoop
logical :: NOTDoSSSS,Segmented,PermuteLHSTypes,PermuteRHSTypes,PermuteODTypes
!Tmporary array used in the innermost routines 
integer :: TMParray1maxsize,TMParray2maxsize
!Batch info
integer :: iBatchC,iBatchD,nBatchA,nBatchB,nBatchC,nBatchD
integer :: iBatchIndexOfTypeA,iBatchIndexOfTypeB,iBatchIndexOfTypeC,iBatchIndexOfTypeD
!Screening variables
integer :: nBatchAGAB,nBatchBGAB,nBatchCGCD,nBatchDGCD
real(realk) :: GABELM
!Passes
integer :: nPasses,iPass,MaxPasses!,nPassLoop,nPassesInPassLoop,iPassTMP !,ndimPass
!integer :: iPassLoop
integer :: MaxTotalAngmom,IAngmomTypes
!ODscreening
real(realk) :: sumExtent2AB,sumExtent2CD,Tstart,Tend
logical :: ODscreen,QQRscreen,dolink
character(len=16)  :: TYPESTRING
!Link specific stuff
real(realk) :: MaxGabLHS,MaxGabRHS !MaximumAllover
real(realk) :: MaxAtomGabelmLHS,MaxAtomGabelmRHS !MaximumAllover
logical :: doMoTrans
!MOtrans specific stuff
integer :: nCMO1,nCMO2,nCMO3,nCMO4

doMOtrans = IchorJobSpec.EQ.IchorJobMOtrans
IF(doMOtrans)Then
   nCMO1 = InputM(InputMspec(1))%n2
   nCMO2 = InputM(InputMspec(2))%n2
   nCMO3 = InputM(InputMspec(3))%n2
   nCMO4 = InputM(InputMspec(4))%n2
   IF(nCMO1.NE.OutputDim1)call ichorQuit('MoTrans Dim Mismatch1',-1)
   IF(nCMO2.NE.OutputDim2)call ichorQuit('MoTrans Dim Mismatch2',-1)
   IF(nCMO3.NE.OutputDim3)call ichorQuit('MoTrans Dim Mismatch3',-1)
   IF(nCMO4.NE.OutputDim4)call ichorQuit('MoTrans Dim Mismatch4',-1)
   IF(OutputDim5.NE.1)call ichorQuit('MoTrans Dim Mismatch5',-1)
ENDIF

doLink = IchorJobSpec.EQ.IchorJobLink
IF(doLink)Then
   IF(OutputDim5.NE.IchorInputDim3)call ichorQuit('Link Dim Mismatch5',-1)
ENDIF

call set_ichor_memvar(MaxMemAllocated,MemAllocated,MaxMem)
INTPRINT=IchorDebugSpec
call mem_ichor_alloc_dryrun(nTypesA) !OrderdListA
call mem_ichor_alloc_dryrun(nTypesB) !OrderdListB
call mem_ichor_alloc_dryrun(nTypesC) !OrderdListC
call mem_ichor_alloc_dryrun(nTypesD) !OrderdListD
PQorder=.FALSE.
call determineScreening(IchorScreenSpec,CSscreen,ODscreen,QQRscreen)
IF(CSScreen)THEN
   call mem_ichor_alloc_dryrun(nTypesA) !BatchIndexOfTypeA
   call mem_ichor_alloc_dryrun(nTypesB) !BatchIndexOfTypeB
   call mem_ichor_alloc_dryrun(nBatchA*nBatchB) !BATCHGAB
   IF(nBatchA.NE.endBatchA-startBatchA+1)call ichorQuit('Screening Dim Mismatch1',-1)
   IF(nBatchB.NE.endBatchB-startBatchB+1)call ichorQuit('Screening Dim Mismatch2',-1)
   call RetrieveGabDimFromIchorSaveGabModule(nBatchAGAB,nBatchBGAB,IchorGabID1)
   IF(nBatchAGAB.NE.nBatchA.OR.nBatchBGAB.NE.nBatchB)THEN
      call mem_ichor_alloc_dryrun(nBatchAGAB*nBatchBGAB) !BATCHGAB2
      call mem_ichor_dealloc_dryrun(nBatchAGAB*nBatchBGAB) !BATCHGAB2
   ENDIF

   call mem_ichor_alloc_dryrun(nTypesA*nTypesB) !MaxGabForTypeAB
   call mem_ichor_alloc_dryrun(nTypesC) !BatchIndexOfTypeC    
   call mem_ichor_alloc_dryrun(nTypesD) !BatchIndexOfTypeD
   call mem_ichor_alloc_dryrun(nTypesC,nTypesD) !MaxGabForTypeCD
   IF(IchorGabID1.EQ.IchorGabID2)THEN
      nBatchC = nBatchA
      nBatchD = nBatchB
   ELSE
      call mem_ichor_alloc_dryrun(nBatchC,nBatchD) !BATCHGCD
      IF(nBatchC.NE.endBatchC-startBatchC+1)call ichorQuit('Screening Dim Mismatch3',-1)
      IF(nBatchD.NE.endBatchD-startBatchD+1)call ichorQuit('Screening Dim Mismatch4',-1)
      call RetrieveGabDimFromIchorSaveGabModule(nBatchCGCD,nBatchDGCD,IchorGabID2)
      IF(nBatchCGCD.NE.nBatchC.OR.nBatchDGCD.NE.nBatchD)THEN
         !WARNING BATCHCALC NOT FULL INTEGRAL IS CALCULATED
         call mem_ichor_alloc_dryrun(nBatchCGCD*nBatchDGCD)   !BATCHGCD2
         call mem_ichor_dealloc_dryrun(nBatchCGCD*nBatchDGCD) !BATCHGCD2
      ENDIF
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
   call mem_ichor_alloc_dryrun(nBatchA,nBatchB) !BATCHGAB
ENDIF

!CSScreen = .FALSE.
!ODscreen = .FALSE.
call determinePermuteSym(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
IF(DoLink)THEN
   IF((SameLHSaos.OR.SameRHSaos).OR.SameODs)call ichorQuit('Ichor Link Symmetry not enabled yet',-1)
ENDIF
IF(doMOtrans)THEN
   !have to deactivate SameOD at least for now. 
   SameODs = .FALSE.
ENDIF
Spherical = SphericalSpec.EQ.SphericalParam
oldmaxangmomABCD = -25
MaxTotalAngmom = MAXVAL(AngmomOfTypeA) + MAXVAL(AngmomOfTypeB) &
     & + MAXVAL(AngmomOfTypeC) + MAXVAL(AngmomOfTypeD)
!we loop over Total angmom in order to ensure that we first do all
!SSSS integrals then PSSS,SPSS,SSPS,SSSP, ...
!this mean we call GAMMATABULATION a limited number of times and
!we reduce branch misprediction inside the code and reuse instruction cache. 
DO IAngmomTypes = 0,MaxTotalAngmom
 DO ItypeD=1,nTypesD
  ! TYPE D CALC ===========================
  call ObtainTypeInfoNoExtent(nTypesD,nAtomsOfTypeD,AngmomOfTypeD,nPrimOfTypeD,&
       & nContOfTypeD,ItypeD,nAtomsD,AngmomD,nPrimD,nContD,nOrbCompD,&
       & nOrbD,nDimD,nCartOrbCompD,spherical)
  IF(nAtomsD.EQ.0)CYCLE
  call mem_ichor_alloc_dryrun(nPrimD) !expD
  call mem_ichor_alloc_dryrun(nPrimD,nContD) !ContractCoeffD
  call mem_ichor_alloc_dryrun(3*nAtomsD) !Dcenter
  call mem_ichor_alloc_dryrun(nAtomsD) !StartOrbitalD
  IF(doMOtrans)Then
     call mem_ichor_alloc_dryrun(nDimD*nCMO4) !CMO4D
     call mem_ichor_alloc_dryrun(nDimD*nCMO3) !CMO3D
  ENDIF
  ! DONE TYPE D CALC ===========================
  DO ItypeC=1,nTypesC
   ! TYPE C CALC ================================
   IF(SameRHSaos .AND. ItypeD.GT.ItypeC)CYCLE
   TriangularRHSAtomLoop = SameRHSaos .AND. ItypeD.EQ.ItypeC
   PermuteRHSTypes = SameRHSaos .AND. ItypeD.LT.ItypeC
   call ObtainTypeInfoNoExtent(nTypesC,nAtomsOfTypeC,AngmomOfTypeC,nPrimOfTypeC,&
        & nContOfTypeC,ItypeC,nAtomsC,AngmomC,nPrimC,nContC,nOrbCompC,&
        & nOrbC,nDimC,nCartOrbCompC,spherical)   
   IF(nAtomsC.EQ.0)CYCLE
   call mem_ichor_alloc_dryrun(nPrimC) !expC
   call mem_ichor_alloc_dryrun(nPrimC,nContC) !ContractCoeffC
   call mem_ichor_alloc_dryrun(3*nAtomsC) !Ccenter
   call mem_ichor_alloc_dryrun(nAtomsC) !StartOrbitalC
   IF(doMOtrans)Then
      call mem_ichor_alloc_dryrun(nDimC*nCMO3) !CMO3C
      call mem_ichor_alloc_dryrun(nDimC*nCMO4) !CMO4C
   ENDIF
   ! DONE TYPE C CALC ===========================
   ! TYPE Q CALC ================================
   AngmomQ = AngmomC + AngmomD
   nPrimQ = nPrimC*nPrimD
   nContQ = nContC*nContD
   nOrbCompQ = nOrbCompC*nOrbCompD
   nCartOrbCompQ = nCartOrbCompC*nCartOrbCompD
   nTUVQ = (AngmomQ+1)*(AngmomQ+2)*(AngmomQ+3)/6
   call mem_ichor_alloc_dryrun(nPrimQ) !Qiprim1
   call mem_ichor_alloc_dryrun(nPrimQ) !Qiprim2
   IF (nContQ.EQ. 1)THEN
      Qsegmented = .TRUE.
   ELSE
      Qsegmented = .FALSE.
   ENDIF
   call mem_ichor_alloc_dryrun(nPrimQ) !expQ

   call mem_ichor_alloc_dryrun(nAtomsC*nAtomsD) !noScreenCD
   call mem_ichor_alloc_dryrun(nAtomsC*nAtomsD) !noScreenCD2
   !LinK: 
   IF(DoLink)THEN !Make Atomic Screening matric AtomGCD
      !maybe this should be used in general - smaller than BATCHGCD ad more straightforward 
      call mem_ichor_alloc_dryrun(nAtomsC*nAtomsD) !AtomGCD
      call mem_ichor_alloc_dryrun(nAtomsD)         !MaxGCDVec
      call mem_ichor_alloc_dryrun(nAtomsD)         !nKetList
      call mem_ichor_alloc_dryrun(nAtomsC*nAtomsD) !KetList
   ENDIF   
   ! DONE TYPE Q CALC ===========================
   DO ItypeB=1,nTypesB
    ! TYPE B CALC ================================
    call ObtainTypeInfoNoExtent(nTypesB,nAtomsOfTypeB,AngmomOfTypeB,nPrimOfTypeB,&
         & nContOfTypeB,ItypeB,nAtomsB,AngmomB,nPrimB,nContB,nOrbCompB,&
         & nOrbB,nDimB,nCartOrbCompB,spherical)
    IF(nAtomsB.EQ.0)CYCLE
    call mem_ichor_alloc_dryrun(nPrimB) !expB
    call mem_ichor_alloc_dryrun(nPrimB*nContB) !ContractCoeffB
    call mem_ichor_alloc_dryrun(3*nAtomsB) !Bcenter
    call mem_ichor_alloc_dryrun(nAtomsB) !StartOrbitalB
    IF(doMOtrans)Then
       call mem_ichor_alloc_dryrun(nDimB*nCMO2) !CMO2B
       call mem_ichor_alloc_dryrun(nDimB*nCMO1) !CMO1B
    ENDIF
    ! DONE TYPE B CALC ===========================
    DO ItypeA=1,nTypesA 
     ! TYPE A CALC ================================
     IF(SameLHSaos .AND. ItypeB.GT.ItypeA)CYCLE
     TriangularLHSAtomLoop = SameLHSaos .AND. ItypeB.EQ.ItypeA
     PermuteLHSTypes = SameLHSaos .AND. ItypeB.LT.ItypeA
     IF(SameODs .AND. ((ItypeC.GT.ItypeA).OR.((ItypeC.EQ.ItypeA).AND.(ItypeD.GT.ItypeB))))CYCLE
     PermuteODTypes = SameODs
!     IF(ItypeC.EQ.ItypeA.AND.ItypeD.EQ.ItypeB)PermuteODTypes=.FALSE.
     TriangularODAtomLoop = SameODs .AND. ((ItypeC.EQ.ItypeA).AND.(ItypeD.EQ.ItypeB))
     call ObtainTypeInfoNoExtent(nTypesA,nAtomsOfTypeA,AngmomOfTypeA,nPrimOfTypeA,&
          & nContOfTypeA,ItypeA,nAtomsA,AngmomA,nPrimA,nContA,nOrbCompA,&
          & nOrbA,nDimA,nCartOrbCompA,spherical)
     IF(nAtomsA.EQ.0)CYCLE
     TotalAngmom = AngmomA + AngmomB + AngmomQ
     IF(TotalAngmom.EQ.IAngmomTypes)THEN
      !This if statement ensures that we call GAMMATABULATION a very limited number of times
      !and it reduceses branch mispredictions inside the code,...
      call mem_ichor_alloc_dryrun(nPrimA) !expA
      call mem_ichor_alloc_dryrun(nPrimA*nContA) !ContractCoeffA
      call mem_ichor_alloc_dryrun(3*nAtomsA) !Acenter
      call mem_ichor_alloc_dryrun(nAtomsA) !StartOrbitalA
      IF(doMOtrans)Then
         call mem_ichor_alloc_dryrun(nDimA*nCMO1) !CMO1A
         call mem_ichor_alloc_dryrun(nDimA*nCMO2) !CMO2A
      ENDIF
      ! DONE TYPE A CALC ===========================
      ! TYPE P CALC ================================
      AngmomP = AngmomA + AngmomB
      nPrimP = nPrimA*nPrimB
      nContP = nContA*nContB
      nOrbCompP = nOrbCompA*nOrbCompB
      nCartOrbCompP = nCartOrbCompA*nCartOrbCompB
      nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
      nTUV = (TotalAngmom+1)*(TotalAngmom+2)*(TotalAngmom+3)/6
      IF (nContP.EQ. 1)THEN
         Psegmented = .TRUE.
      ELSE
         Psegmented = .FALSE.
      ENDIF
      call mem_ichor_alloc_dryrun(nPrimP) !expP
      call mem_ichor_alloc_dryrun(nPrimP) !inversexpP
      !LinK: 
      IF(DoLink)THEN !Make Atomic Screening matric AtomGAB
         !maybe this should be used in general - smaller than BATCHGCD ad more straightforward 
         call mem_ichor_alloc_dryrun(nAtomsA*nAtomsB) !AtomGAB
         call mem_ichor_alloc_dryrun(nAtomsB) !MaxGABVec
      ENDIF
      ! DONE TYPE P CALC ===========================
      
      ! TYPE PQ CALC ================================
      TotalAngmom = AngmomA + AngmomB + AngmomQ
      maxangmomABCD = AngmomA + AngmomB + AngmomQ

      UseGeneralCode = .FALSE. !Use Specialized code when appropriate. 
      IF(GGemOperatorCalc)UseGeneralCode = .TRUE.
      IF(MAX(AngmomA,AngmomB,AngmomC,AngmomD).GT.MaxSpecialAngmom)UseGeneralCode = .TRUE.
      IF(maxangmomABCD.NE.oldmaxangmomABCD)THEN
       IF(oldmaxangmomABCD.NE.-25)THEN
          call mem_ichor_dealloc_dryrun((nTABFJW1+1)*(nTABFJW2+1)) !TABFJW
       ENDIF
       nTABFJW1 = AngmomA + AngmomB + AngmomQ + 3 
       !only need + 3 after Branos change in BUILD_RJ000 
       nTABFJW2 = 1200
       !TABFJW(0:nTABFJW1,0:nTABFJW2)
       call mem_ichor_alloc_dryrun((nTABFJW1+1)*(nTABFJW2+1)) !TABFJW
       oldmaxangmomABCD = maxangmomABCD
      ENDIF
      call mem_ichor_alloc_dryrun(nPrimP*nPrimQ) !reducedExponents
      call mem_ichor_alloc_dryrun(nPrimP*nPrimQ) !integralPrefactor
      ! DONE TYPE PQ CALC ===========================      
      NOTDoSSSS = .NOT.(TotalAngmom.EQ.0.AND.(Psegmented.AND.Qsegmented))
      IF(DoLink) NOTDoSSSS=.TRUE.
      IF(DoMoTrans) NOTDoSSSS=.TRUE.
      IF(NOTDoSSSS)THEN
         !Determine Sizes of TmpArrays and MaxPasses
         IF(UseCPU)THEN
            call ICI_CPU_OBS_general_size(TMParray1maxsize,&
                 & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
                 & nContA,nContB,nContC,nContD,&
                 & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,nContP,&
                 & nContQ,nPrimQ*nPrimP,nContQ*nContP,Psegmented,Qsegmented)
         ELSE !use GPU code
            call ICI_GPU_OBS_general_size(TMParray1maxsize,&
                 & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
                 & nContA,nContB,nContC,nContD,&
                 & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,nContP,&
                 & nContQ,nPrimQ*nPrimP,nContQ*nContP,Psegmented,Qsegmented)
         ENDIF
         nLocalInt = nOrbA*nOrbB*nOrbC*nOrbD
      ENDIF
      IF(NOTDoSSSS)THEN !FIXME A MESS
         IF(DoLink)THEN         
            call mem_ichor_alloc_dryrun(nPrimQ) !QpreExpFac
            call mem_ichor_alloc_dryrun(3*nPrimQ) !Qcent
         ELSEIF(doMOtrans)THEN
            call mem_ichor_alloc_dryrun(nPrimQ) !QpreExpFac
            call mem_ichor_alloc_dryrun(3*nPrimQ) !Qcent
         ELSEIF(UseCPU)THEN
            call mem_ichor_alloc_dryrun(nPrimQ) !QpreExpFac
            call mem_ichor_alloc_dryrun(3*nPrimQ) !Qcent
         ELSE
            !not needed 
         ENDIF
      ELSE
         call mem_ichor_alloc_dryrun(nPrimQ) !QpreExpFac
         call mem_ichor_alloc_dryrun(3*nPrimQ) !Qcent
      ENDIF
      !calc       
      call mem_ichor_alloc_dryrun(nAtomsA*nAtomsB) !noScreenAB
      
      !LinK: 
      IF(DoLink)THEN !Form active Dmat: DmatBD,DmatAD,DmatBC,DmatAC
         call mem_ichor_alloc_dryrun(nDimB*nDimD*IchorInputDim3) !DmatBD
         call mem_ichor_alloc_dryrun(nDimA*nDimD*IchorInputDim3) !DmatAD
         call mem_ichor_alloc_dryrun(nDimB*nDimC*IchorInputDim3) !DmatBC
         call mem_ichor_alloc_dryrun(nDimA*nDimC*IchorInputDim3) !DmatAC
         call mem_ichor_alloc_dryrun(nDimA,nDimC,IchorInputDim3) !KmatAC
         call mem_ichor_alloc_dryrun(nDimB,nDimC,IchorInputDim3) !KmatBC
         call mem_ichor_alloc_dryrun(nDimA,nDimD,IchorInputDim3) !KmatAD
         call mem_ichor_alloc_dryrun(nDimB,nDimD,IchorInputDim3) !KmatBD
         call mem_ichor_alloc_dryrun(nAtomsB,nAtomsD) !ReducedDmatBD
         call mem_ichor_alloc_dryrun(nAtomsB,nAtomsC) !ReducedDmatBC
         call mem_ichor_alloc_dryrun(nAtomsB) !nBraList
         call mem_ichor_alloc_dryrun(nAtomsA,nAtomsB) !BraList
         call mem_ichor_alloc_dryrun(nAtomsD) !nBraketList
         call mem_ichor_alloc_dryrun(nAtomsB,nAtomsD) !BraketList
      ENDIF
      call mem_ichor_alloc_dryrun(nPrimP,nAtomsA*nAtomsB) !PpreExpFacPass
      call mem_ichor_alloc_dryrun(3*nPrimP,nAtomsA*nAtomsB) !PcentPass
      IF(NOTDoSSSS)THEN
         call mem_ichor_alloc_dryrun(3,nAtomsA*nAtomsB) !Pdistance12Pass
         call mem_ichor_alloc_dryrun(nPrimP) !PpreExpFac
         IF(DoLink)THEN !call IchorTypeLinKLoop
            !WARNING assumes MaxPasses = nAtomsA*nAtomsB
            MaxPasses = nAtomsA*nAtomsB
            call mem_ichor_alloc_dryrun(TMParray1maxsize*MaxPasses) !TmpArray1
            call mem_ichor_alloc_dryrun(TMParray2maxsize*MaxPasses) !TmpArray2
            call mem_ichor_alloc_dryrun(MaxPasses) !IatomAPass
            call mem_ichor_alloc_dryrun(MaxPasses) !IatomBPass
            call mem_ichor_alloc_dryrun(nLocalint*MaxPasses) !LocalIntPass1
            call mem_ichor_alloc_dryrun(nLocalInt*nAtomsA*nAtomsB) !LocalIntPass2
            call mem_ichor_alloc_dryrun(nAtomsA*nAtomsB) !DoINT
            IF(UseGeneralCode)THEN
               call DetermineSizeTmpArray34(nTUVQ,nCartOrbCompQ,nPrimQ,nTUVP,nCartOrbCompP,&
                    & nPrimP,MaxPasses,AngmomA,AngmomB,AngmomC,AngmomD,&
                    & AngmomA+AngmomB,AngmomC+AngmomD,TotalAngmom)
               call mem_ichor_alloc_dryrun(nTmpArray3) !TmpArray3
               call mem_ichor_alloc_dryrun(nTmpArray4) !TmpArray4
               !     CALL PreCalciChorSPHMAT(MAX(AngmomA,AngmomB,AngmomC,AngmomD))
            ENDIF
            !call IchorTypeLinKLoop
            IF(UseGeneralCode)THEN
               call mem_ichor_dealloc_dryrun(nTmpArray3) !TmpArray3
               call mem_ichor_dealloc_dryrun(nTmpArray4) !TmpArray4
            ENDIF
            call mem_ichor_dealloc_dryrun(nAtomsA*nAtomsB) !DoINT
            call mem_ichor_dealloc_dryrun(TMParray1maxsize*MaxPasses) !TmpArray1
            call mem_ichor_dealloc_dryrun(TMParray2maxsize*MaxPasses) !TmpArray2
            call mem_ichor_dealloc_dryrun(MaxPasses) !IatomAPass
            call mem_ichor_dealloc_dryrun(MaxPasses) !IatomBPass
            call mem_ichor_dealloc_dryrun(nLocalInt*MaxPasses) !LocalIntPass1
            call mem_ichor_dealloc_dryrun(nLocalInt*nAtomsA*nAtomsB) !LocalIntPass2
         ELSEIF(doMOtrans)THEN!call IchorTypeMOtransLoop
            !WARNING assumes MaxPasses = nAtomsA*nAtomsB
            MaxPasses = nAtomsA*nAtomsB
            call mem_ichor_alloc_dryrun(TMParray1maxsize*MaxPasses) !TmpArray1
            call mem_ichor_alloc_dryrun(TMParray2maxsize*MaxPasses) !TmpArray2
            call mem_ichor_alloc_dryrun(MaxPasses) !IatomAPass
            call mem_ichor_alloc_dryrun(MaxPasses) !IatomBPass
            call mem_ichor_alloc_dryrun(nLocalint*MaxPasses) !LocalIntPass1
            call mem_ichor_alloc_dryrun(nLocalInt*nAtomsA*nAtomsB) !LocalIntPass2
            IF(UseGeneralCode)THEN
               call DetermineSizeTmpArray34(nTUVQ,nCartOrbCompQ,nPrimQ,nTUVP,nCartOrbCompP,&
                    & nPrimP,MaxPasses,AngmomA,AngmomB,AngmomC,AngmomD,&
                    & AngmomA+AngmomB,AngmomC+AngmomD,TotalAngmom)
               call mem_ichor_alloc_dryrun(nTmpArray3) !TmpArray3
               call mem_ichor_alloc_dryrun(nTmpArray4) !TmpArray4
               !     CALL PreCalciChorSPHMAT(MAX(AngmomA,AngmomB,AngmomC,AngmomD))
            ENDIF
            call mem_ichor_alloc_dryrun(nCMO1,MAX(ndimA,ndimB),nOrbC*nOrbD) !OutputA
            call mem_ichor_alloc_dryrun(nCMO1,nCMO2,ndimC*nDimD) !OutputCD
            !call IchorTypeMOtransLoop
            IF(UseGeneralCode)THEN
               call mem_ichor_dealloc_dryrun(nTmpArray3) !TmpArray3
               call mem_ichor_dealloc_dryrun(nTmpArray4) !TmpArray4
            ENDIF
            call mem_ichor_dealloc_dryrun(nCMO1,MAX(ndimA,ndimB),nOrbC*nOrbD) !OutputA
            call mem_ichor_dealloc_dryrun(TMParray1maxsize*MaxPasses) !TmpArray1
            call mem_ichor_dealloc_dryrun(TMParray2maxsize*MaxPasses) !TmpArray2
            call mem_ichor_dealloc_dryrun(MaxPasses) !IatomAPass
            call mem_ichor_dealloc_dryrun(MaxPasses) !IatomBPass
            call mem_ichor_dealloc_dryrun(nLocalInt*MaxPasses) !LocalIntPass1
            call mem_ichor_dealloc_dryrun(nLocalInt*nAtomsA*nAtomsB) !LocalIntPass2            
            IF(PermuteRHSTypes)THEN
               call mem_ichor_alloc_dryrun(nCMO1,nCMO2,nCMO3*MAX(nDimD,nDimC)) !OutputC
               call mem_ichor_dealloc_dryrun(nCMO1,nCMO2,ndimC*nDimD) !OutputCD
               call mem_ichor_dealloc_dryrun(nCMO1,nCMO2,nCMO3*MAX(nDimD,nDimC)) !OutputC
            ELSE
               call mem_ichor_alloc_dryrun(nCMO1,nCMO2,nCMO3*nDimD) !OutputC
               call mem_ichor_dealloc_dryrun(nCMO1,nCMO2,ndimC*nDimD) !OutputCD
               call mem_ichor_dealloc_dryrun(nCMO1,nCMO2,nCMO3*nDimD) !OutputC
            ENDIF
         ELSEIF(UseCPU)THEN !IchorTypeIntegralLoopCPU
            !WARNING assumes MaxPasses = nAtomsA*nAtomsB
            MaxPasses = nAtomsA*nAtomsB
            call mem_ichor_alloc_dryrun(TMParray1maxsize*MaxPasses) !TmpArray1
            call mem_ichor_alloc_dryrun(TMParray2maxsize*MaxPasses) !TmpArray2
            call mem_ichor_alloc_dryrun(MaxPasses) !IatomAPass
            call mem_ichor_alloc_dryrun(MaxPasses) !IatomBPass
            call mem_ichor_alloc_dryrun(nLocalint*MaxPasses) !LocalIntPass1
            call mem_ichor_alloc_dryrun(nLocalInt*nAtomsA*nAtomsB) !LocalIntPass2
            IF(UseGeneralCode)THEN
               call DetermineSizeTmpArray34(nTUVQ,nCartOrbCompQ,nPrimQ,nTUVP,nCartOrbCompP,&
                    & nPrimP,MaxPasses,AngmomA,AngmomB,AngmomC,AngmomD,&
                    & AngmomA+AngmomB,AngmomC+AngmomD,TotalAngmom)
               call mem_ichor_alloc_dryrun(nTmpArray3) !TmpArray3
               call mem_ichor_alloc_dryrun(nTmpArray4) !TmpArray4
               !     CALL PreCalciChorSPHMAT(MAX(AngmomA,AngmomB,AngmomC,AngmomD))
            ENDIF
            !IchorTypeIntegralLoopCPU
            IF(UseGeneralCode)THEN
               call mem_ichor_dealloc_dryrun(nTmpArray3) !TmpArray3
               call mem_ichor_dealloc_dryrun(nTmpArray4) !TmpArray4
            ENDIF
            call mem_ichor_dealloc_dryrun(TMParray1maxsize*MaxPasses) !TmpArray1
            call mem_ichor_dealloc_dryrun(TMParray2maxsize*MaxPasses) !TmpArray2
            call mem_ichor_dealloc_dryrun(MaxPasses) !IatomAPass
            call mem_ichor_dealloc_dryrun(MaxPasses) !IatomBPass
            call mem_ichor_dealloc_dryrun(nLocalInt*MaxPasses) !LocalIntPass1
            call mem_ichor_dealloc_dryrun(nLocalInt*nAtomsA*nAtomsB) !LocalIntPass2
         ELSE ! IchorTypeIntegralLoopGPU
            !WARNING assumes MaxPasses = nAtomsA*nAtomsB
            !WARNING assumes nAsyncHandles = maxnAsyncHandles
            call mem_ichor_alloc_dryrun(nLocalInt*nAtomsA*nAtomsB) !LocalIntPass2
            call mem_ichor_alloc_dryrun(3*nPrimQ*maxnAsyncHandles) !Qcent
            call mem_ichor_alloc_dryrun(3*maxnAsyncHandles) !Qdistance12
            call mem_ichor_alloc_dryrun(nPrimQ,maxnAsyncHandles) !QpreExpFac
            call mem_ichor_alloc_dryrun(TMParray1maxsize*MaxPasses,maxnAsyncHandles) !TmpArray1
            call mem_ichor_alloc_dryrun(TMParray2maxsize*MaxPasses,maxnAsyncHandles) !TmpArray2
            call mem_ichor_alloc_dryrun(MaxPasses*maxnAsyncHandles) !IatomAPass
            call mem_ichor_alloc_dryrun(MaxPasses*maxnAsyncHandles) !IatomBPass
            call mem_ichor_alloc_dryrun(nLocalInt*nAtomsA*nAtomsB*maxnAsyncHandles) !LocalIntPass1
            !IchorTypeIntegralLoopGPU
            call mem_ichor_dealloc_dryrun(nLocalInt*nAtomsA*nAtomsB) !LocalIntPass2
            call mem_ichor_dealloc_dryrun(3*nPrimQ*maxnAsyncHandles) !Qcent
            call mem_ichor_dealloc_dryrun(3*maxnAsyncHandles) !Qdistance12
            call mem_ichor_dealloc_dryrun(nPrimQ,maxnAsyncHandles) !QpreExpFac
            call mem_ichor_dealloc_dryrun(TMParray1maxsize*MaxPasses,maxnAsyncHandles) !TmpArray1
            call mem_ichor_dealloc_dryrun(TMParray2maxsize*MaxPasses,maxnAsyncHandles) !TmpArray2
            call mem_ichor_dealloc_dryrun(MaxPasses*maxnAsyncHandles) !IatomAPass
            call mem_ichor_dealloc_dryrun(MaxPasses*maxnAsyncHandles) !IatomBPass
            call mem_ichor_dealloc_dryrun(nLocalInt*nAtomsA*nAtomsB*maxnAsyncHandles) !LocalIntPass1
         ENDIF
         call mem_ichor_dealloc_dryrun(nPrimP) !PpreExpFac
         call mem_ichor_dealloc_dryrun(3,nAtomsA*nAtomsB) !Pdistance12Pass
      ELSE
         call mem_ichor_alloc_dryrun(nOrbA*nAtomsA*nOrbB*nAtomsB,nOrbC*nOrbD) !LocalIntPass
         call mem_ichor_alloc_dryrun(nAtomsA*nAtomsB) !IatomAPass
         call mem_ichor_alloc_dryrun(nAtomsA*nAtomsB) !IatomBPass
!         call IchorsegsegSSSSIntegralLoop
         call mem_ichor_dealloc_dryrun(nAtomsA*nAtomsB) !IatomAPass
         call mem_ichor_dealloc_dryrun(nAtomsA*nAtomsB) !IatomBPass
         call mem_ichor_dealloc_dryrun(nOrbA*nAtomsA*nOrbB*nAtomsB,nOrbC*nOrbD) !LocalIntPass
      ENDIF
      IF(NOTDoSSSS)THEN !FIXME A MESS
         IF(DoLink)THEN         
            call mem_ichor_dealloc_dryrun(3*nPrimQ) !Qcent
            call mem_ichor_dealloc_dryrun(nPrimQ) !QpreExpFac
         ELSEIF(doMOtrans)THEN
            call mem_ichor_dealloc_dryrun(3*nPrimQ) !Qcent
            call mem_ichor_dealloc_dryrun(nPrimQ) !QpreExpFac
         ELSEIF(UseCPU)THEN
            call mem_ichor_dealloc_dryrun(3*nPrimQ) !Qcent
            call mem_ichor_dealloc_dryrun(nPrimQ) !QpreExpFac
         ELSE
            !not needed 
         ENDIF
      ELSE
         call mem_ichor_dealloc_dryrun(3*nPrimQ) !Qcent
         call mem_ichor_dealloc_dryrun(nPrimQ) !QpreExpFac
      ENDIF
      call mem_ichor_dealloc_dryrun(nPrimP,nAtomsA*nAtomsB) !PpreExpFacPass
      call mem_ichor_dealloc_dryrun(3*nPrimP,nAtomsA*nAtomsB) !PcentPass
      IF(DoLink)THEN         
         call mem_ichor_dealloc_dryrun(nDimB*nDimD*IchorInputDim3) !DmatBD
         call mem_ichor_dealloc_dryrun(nDimA*nDimD*IchorInputDim3) !DmatAD
         call mem_ichor_dealloc_dryrun(nDimB*nDimC*IchorInputDim3) !DmatBC
         call mem_ichor_dealloc_dryrun(nDimA*nDimC*IchorInputDim3) !DmatAC
         call mem_ichor_dealloc_dryrun(nDimA,nDimC,IchorInputDim3) !KmatAC
         call mem_ichor_dealloc_dryrun(nDimB,nDimC,IchorInputDim3) !KmatBC
         call mem_ichor_dealloc_dryrun(nDimA,nDimD,IchorInputDim3) !KmatAD
         call mem_ichor_dealloc_dryrun(nDimB,nDimD,IchorInputDim3) !KmatBD
         call mem_ichor_dealloc_dryrun(nAtomsB,nAtomsD) !ReducedDmatBD
         call mem_ichor_dealloc_dryrun(nAtomsB,nAtomsC) !ReducedDmatBC
         call mem_ichor_dealloc_dryrun(nAtomsB) !nBraList
         call mem_ichor_dealloc_dryrun(nAtomsA,nAtomsB) !BraList
         call mem_ichor_dealloc_dryrun(nAtomsD) !nBraketList
         call mem_ichor_dealloc_dryrun(nAtomsB,nAtomsD) !BraketList
      ENDIF
      call mem_ichor_dealloc_dryrun(nPrimP*nPrimQ) !reducedExponents
      call mem_ichor_dealloc_dryrun(nPrimP*nPrimQ) !integralPrefactor
      call mem_ichor_dealloc_dryrun(nAtomsA*nAtomsB) !noScreenAB

      call mem_ichor_dealloc_dryrun(nPrimP) !expP
      call mem_ichor_dealloc_dryrun(nPrimP) !inversexpP

      call mem_ichor_dealloc_dryrun(nPrimA) !expA
      call mem_ichor_dealloc_dryrun(nPrimA*nContA) !ContractCoeffA
      call mem_ichor_dealloc_dryrun(3*nAtomsA) !Acenter
      call mem_ichor_dealloc_dryrun(nAtomsA) !StartOrbitalA
      IF(doMOtrans)Then
         call mem_ichor_dealloc_dryrun(nDimA*nCMO1) !CMO1A
         call mem_ichor_dealloc_dryrun(nDimA*nCMO2) !CMO2A
      ENDIF
      IF(DoLink)THEN
         call mem_ichor_dealloc_dryrun(nAtomsA*nAtomsB) !AtomGAB
         call mem_ichor_dealloc_dryrun(nAtomsB) !MaxGABVec
      ENDIF
     ENDIF !correct angmom type
    ENDDO !typeA
    call mem_ichor_dealloc_dryrun(nPrimB) !expB
    call mem_ichor_dealloc_dryrun(nPrimB*nContB) !ContractCoeffB
    call mem_ichor_dealloc_dryrun(3*nAtomsB) !Bcenter
    call mem_ichor_dealloc_dryrun(nAtomsB) !StartOrbitalB
    IF(doMOtrans)Then
       call mem_ichor_dealloc_dryrun(nDimB*nCMO1) !CMO1B
       call mem_ichor_dealloc_dryrun(nDimB*nCMO2) !CMO2B
    ENDIF
   ENDDO !typeB
   IF(DoLink)THEN
      call mem_ichor_dealloc_dryrun(nAtomsC*nAtomsD) !AtomGCD
      call mem_ichor_dealloc_dryrun(nAtomsD)         !MaxGCDVec
      call mem_ichor_dealloc_dryrun(nAtomsD)         !nKetList
      call mem_ichor_dealloc_dryrun(nAtomsC*nAtomsD) !KetList
   ENDIF   
   call mem_ichor_dealloc_dryrun(nAtomsC*nAtomsD) !noScreenCD
   call mem_ichor_dealloc_dryrun(nAtomsC*nAtomsD) !noScreenCD2
   call mem_ichor_dealloc_dryrun(nPrimQ) !expQ
   call mem_ichor_dealloc_dryrun(nPrimQ) !Qiprim1
   call mem_ichor_dealloc_dryrun(nPrimQ) !Qiprim2

   call mem_ichor_dealloc_dryrun(nPrimC) !expC
   call mem_ichor_dealloc_dryrun(nPrimC,nContC) !ContractCoeffC
   call mem_ichor_dealloc_dryrun(3*nAtomsC) !Ccenter
   call mem_ichor_dealloc_dryrun(nAtomsC) !StartOrbitalC
   IF(doMOtrans)Then
      call mem_ichor_dealloc_dryrun(nDimC*nCMO3) !CMO3C
      call mem_ichor_dealloc_dryrun(nDimC*nCMO4) !CMO4C
   ENDIF
  ENDDO !typeC
  call mem_ichor_dealloc_dryrun(nPrimD) !expD
  call mem_ichor_dealloc_dryrun(nPrimD,nContD) !ContractCoeffD
  call mem_ichor_dealloc_dryrun(3*nAtomsD) !Dcenter
  call mem_ichor_dealloc_dryrun(nAtomsD) !StartOrbitalD
  IF(doMOtrans)Then
     call mem_ichor_dealloc_dryrun(nDimD*nCMO4) !CMO4D
     call mem_ichor_dealloc_dryrun(nDimD*nCMO3) !CMO3D
  ENDIF
 ENDDO !typeD
ENDDO 
call mem_ichor_dealloc_dryrun((nTABFJW1+1)*(nTABFJW2+1)) !TABFJW
call mem_ichor_dealloc_dryrun(nBatchA,nBatchB) !BATCHGAB

IF(CSScreen)THEN
   call mem_ichor_alloc_dryrun(nTypesA) !BatchIndexOfTypeA
   call mem_ichor_alloc_dryrun(nTypesB) !BatchIndexOfTypeB
   call mem_ichor_alloc_dryrun(nBatchA*nBatchB) !BATCHGAB
   call mem_ichor_alloc_dryrun(nTypesA*nTypesB) !MaxGabForTypeAB
   call mem_ichor_alloc_dryrun(nTypesC) !BatchIndexOfTypeC    
   call mem_ichor_alloc_dryrun(nTypesD) !BatchIndexOfTypeD
   call mem_ichor_alloc_dryrun(nTypesC,nTypesD) !MaxGabForTypeCD
   IF(IchorGabID1.EQ.IchorGabID2)THEN
   ELSE
      call mem_ichor_alloc_dryrun(nBatchC,nBatchD) !BATCHGCD
   ENDIF
ENDIF
call mem_ichor_dealloc_dryrun(nTypesA) !OrderdListA
call mem_ichor_dealloc_dryrun(nTypesB) !OrderdListB
call mem_ichor_dealloc_dryrun(nTypesC) !OrderdListC
call mem_ichor_dealloc_dryrun(nTypesD) !OrderdListD

call retrieve_ichor_memvar(MaxMemAllocated,MemAllocated)
IF(MemAllocated.NE.0)THEN
   call ichorquit('MemoryLeak in IchorEri',lupri)
ENDIF
end subroutine IchorEriMem

subroutine IchorTypeIntegralLoopCPU(nAtomsA,nPrimA,nContA,nOrbCompA,startOrbitalA,&
     & iBatchIndexOfTypeA,expA,ContractCoeffA,AngmomA,Acenter,nBatchA,nOrbA,&
     & nAtomsB,nPrimB,nContB,nOrbCompB,startOrbitalB,&
     & iBatchIndexOfTypeB,expB,ContractCoeffB,AngmomB,Bcenter,nBatchB,nOrbB,&
     & nAtomsC,nPrimC,nContC,nOrbCompC,startOrbitalC,&
     & iBatchIndexOfTypeC,expC,ContractCoeffC,AngmomC,Ccenter,nBatchC,nOrbC,&
     & nAtomsD,nPrimD,nContD,nOrbCompD,startOrbitalD,&
     & iBatchIndexOfTypeD,expD,ContractCoeffD,AngmomD,Dcenter,nBatchD,nOrbD,&
     & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
     & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
     & nPrimP,nPrimQ,nContP,nContQ,expP,expQ,&
     & qcent,Ppreexpfac,Qpreexpfac,&
     & Qiprim1,Qiprim2,&
     & nTABFJW1,nTABFJW2,TABFJW,TotalAngmom,&
     & CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
     & Qsegmented,Psegmented,TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
     & reducedExponents,integralPrefactor,&
     & PcentPass,Pdistance12Pass,PpreExpFacPass,&
     & Qdistance12,PQorder,&
     & BATCHGCD,BATCHGAB,Spherical,TMParray1maxsize,nLocalInt,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,OutputStorage,&
     & TMParray2maxsize,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri)
  implicit none
  logical,intent(in) :: TriangularRHSAtomLoop,TriangularLHSAtomLoop,PermuteLHSTypes
  logical,intent(in) :: Qsegmented,Psegmented,PQorder,Spherical,CSscreen
  logical,intent(in) :: TriangularODAtomLoop
  integer,intent(in) :: nTABFJW1,nTABFJW2,lupri
  real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2),THRESHOLD_CS
  integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ)
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
  real(realk),intent(inout) :: Qcent(3,nPrimQ),Qdistance12(3),QpreExpFac(nPrimQ),expQ(nPrimQ)
  !P
  integer,intent(in) :: nContP,nPrimP
  real(realk),intent(inout) :: PpreExpFac(nPrimP)
  real(realk),intent(in) :: PcentPass(3,nPrimP,natomsA*natomsB)
  real(realk),intent(in) :: Pdistance12Pass(3,natomsA*natomsB)
  real(realk),intent(in) :: PpreExpFacPass(nPrimP,natomsA*natomsB),expP(nPrimP)
  !
  real(realk),intent(in) :: BATCHGCD(nBatchC,nBatchD)
  !A & B
  integer,intent(in) :: nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV
  integer,intent(in) :: nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD
  integer,intent(in) :: iBatchIndexOfTypeA,iBatchIndexOfTypeB,AngmomB,AngmomA,nOrbA,nOrbB
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbCompA,nOrbCompB,nBatchA,nBatchB,nContA,nContB
  integer,intent(in) :: nPrimA,nPrimB
  integer,intent(in) :: startOrbitalA(nAtomsA)
  integer,intent(in) :: startOrbitalB(nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA),expA(nPrimA),ContractCoeffA(nPrimA,nContA)
  real(realk),intent(in) :: Bcenter(3,nAtomsB),expB(nPrimB),ContractCoeffB(nPrimB,nContB)
  real(realk),intent(in) ::  BATCHGAB(nBatchA*nBatchB)
  !collected
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,TotalAngmom
  real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  integer,intent(in) :: TMParray1maxsize,TMParray2maxsize,nLocalInt
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: integralPrefactor(nPrimQ,nPrimP)
  !local variables
  integer :: iBatchD,IatomD,IatomC,startC,IatomB,nPasses,startD,intprint
  real(realk) :: DcenterSpec(3),CcenterSpec(3),AcenterSpec(3),BcenterSpec(3),GABELM
  logical :: PermuteRHS
  real(realk),allocatable :: TmpArray1(:),TmpArray2(:)
  real(realk),allocatable :: LocalIntPass1(:),LocalIntPass2(:)
  integer,allocatable :: IatomAPass(:),IatomBPass(:)
  integer :: iOrbQ,iOrbB,iOrbA,iOrbD,iOrbC,I4,I3,I2
  integer :: startA,startB,ndim,nOrbQ,MaxPasses
  integer :: TMParray1maxsizePass,TMParray2maxsizePass,nLocalIntPass
#ifdef VAR_OMP
  integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
  nLocalIntPass = nOrbA*nAtomsA*nOrbB*nAtomsB*nOrbC*nOrbD
  call DetermineMaxPasses(nAtomsD,iBatchIndexOfTypeD,nAtomsC,nAtomsA,nAtomsB,&
       & iBatchIndexOfTypeC,iBatchIndexOfTypeA,nBatchB,nBatchA,iBatchIndexOfTypeB,&
       & TriangularRHSAtomLoop,CSscreen,TriangularLHSAtomLoop,TriangularODAtomLoop,&
       & noScreenCD2,BATCHGAB,THRESHOLD_CS,noScreenAB,BATCHGCD,nBatchC,nBatchD,MaxPasses)
  TMParray1maxsizePass = TMParray1maxsize*MaxPasses
  TMParray2maxsizePass = TMParray2maxsize*MaxPasses

  ndim = nOrbA*nOrbB*nAtomsA*nAtomsB
  nOrbQ = nOrbC*nOrbD
  allocate(TmpArray1(TMParray1maxsize*MaxPasses))
  call mem_ichor_alloc(TmpArray1)
  allocate(TmpArray2(TMParray2maxsize*MaxPasses))     
  call mem_ichor_alloc(TmpArray2)
  allocate(IatomAPass(MaxPasses))
  call mem_ichor_alloc(IatomAPass)  
  allocate(IatomBPass(MaxPasses))
  call mem_ichor_alloc(IatomBPass)

  nLocalIntPass = nLocalint*MaxPasses
  allocate(LocalIntPass1(nLocalIntPass))
  CALL Mem_ichor_alloc(LocalIntPass1)
  allocate(LocalIntPass2(nLocalint*nAtomsA*nAtomsB))
  CALL Mem_ichor_alloc(LocalIntPass2)
  IF(UseGeneralCode)THEN
     call DetermineSizeTmpArray34(nTUVQ,nCartOrbCompQ,nPrimQ,nTUVP,nCartOrbCompP,nPrimP,MaxPasses,&
          & AngmomA,AngmomB,AngmomC,AngmomD,AngmomA+AngmomB,AngmomC+AngmomD,TotalAngmom)
     allocate(TmpArray3(nTmpArray3))
     call mem_ichor_alloc(TmpArray3)
     allocate(TmpArray4(nTmpArray4))
     call mem_ichor_alloc(TmpArray4)
     CALL PreCalciChorSPHMAT(MAX(AngmomA,AngmomB,AngmomC,AngmomD))
  ENDIF
!$OMP PARALLEL DEFAULT(none) &
!$OMP PRIVATE(iAtomD,iAtomC,GABELM,startD,iBatchD,DcenterSpec,PermuteRHS,startC,&
!$OMP         CcenterSpec,iOrbQ,I3,startB,I4,iOrbD,iOrbC,iAtomB) &
!$OMP SHARED(nAtomsD,startOrbitalD,iBatchIndexOfTypeD,Dcenter,nAtomsC,&
!$OMP        TriangularRHSAtomLoop,startOrbitalC,noScreenCD2,Ccenter,&
!$OMP        nAtomsA,nAtomsB,BATCHGCD,iBatchIndexOfTypeC,CSscreen,&
!$OMP        nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,&
!$OMP        BATCHGAB,THRESHOLD_CS,IatomAPass,IatomBPass,startOrbitalB,&
!$OMP        MaxPasses,TriangularLHSAtomLoop,TriangularODAtomLoop,nOrbB,&
!$OMP        Qsegmented,nPasses,noScreenAB,nLocalInt,TotalAngmom,nOrbA,&
!$OMP        nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,intprint,lupri,&
!$OMP        nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,TABFJW,&
!$OMP        ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
!$OMP        pcentPass,qcent,PpreexpfacPass,Qpreexpfac,nTABFJW1,nTABFJW2,&
!$OMP        Qiprim1,Qiprim2,expA,expB,expC,expD,Psegmented,reducedExponents,&
!$OMP        integralPrefactor,AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12Pass,&
!$OMP        Qdistance12,PQorder,LocalIntPass1,LocalIntPass2,nLocalIntPass,&
!$OMP        Spherical,TmpArray1,TMParray1maxsizePass,TmpArray2,Bcenter,nOrbQ,&
!$OMP        TMParray2maxsizePass,Acenter,nTmpArray3,nTmpArray4,&
!$OMP        nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,PermuteLHSTypes,nOrbD,nOrbC,&
!$OMP        startOrbitalA,OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputStorage,&
!$OMP        nTUVQ,nCartOrbCompQ,nTUVP,nCartOrbCompP,TmpArray3,TmpArray4,nTUV,&
!$OMP        nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,nOrbCompP,nOrbCompQ)
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
     IF(CSscreen)GABELM = BATCHGCD(iBatchIndexOfTypeC+IatomC,iBatchD)
     !output: IatomAPass,IatomBPass,nPasses
!$OMP SINGLE
     nPasses = nAtomsA*nAtomsB
     CALL BUILD_noScreen2(CSscreen,nAtomsA,nAtomsB,&
          & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,&
          & BATCHGAB,THRESHOLD_CS,GABELM,nPasses,IatomAPass,IatomBPass,MaxPasses,&
          & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,IatomD,noScreenAB) 
!$OMP END SINGLE
!$OMP BARRIER
     IF(nPasses.EQ.0)CYCLE
     IF(nPasses.NE.nAtomsA*nAtomsB)THEN
!$OMP DO PRIVATE(I4)
        do I4 = 1,nLocalIntPass
           LocalIntPass1(I4) = 0.0E0_realk
        enddo
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(I4)
        do I4 = 1,nLocalint*nAtomsA*nAtomsB
           LocalIntPass2(I4) = 0.0E0_realk
        enddo
!$OMP END DO
     ENDIF
!$OMP SINGLE
     !output: Qcent,Qdistance12,QpreExpFac
     CALL Build_qcent_Qdistance12_QpreExpFac(nPrimC,nPrimD,nContC,nContD,&
          & expC,expD,CcenterSpec,DcenterSpec,ContractCoeffC,ContractCoeffD,&
          & Qsegmented,Qcent,Qdistance12,QpreExpFac,INTPRINT)
!$OMP END SINGLE
!$OMP BARRIER
     !Unique for each iPassQ (iAtomC,iAtomD) iteration: qcent,qdistance12,qpreexpfac, Qiprim1(nPrimQ), output:
     !LocalIntPass(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,MaxPasses)
     !IatomAPass,iatomBPass changes and 
!     IF(iAtomC.EQ.1.AND.iAtomD.EQ.1)INTPRINT=1000
     call ICI_CPU_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,&
          & nPrimQ,nPrimP*nPrimQ,nPasses,MaxPasses,intprint,lupri,&
          & nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,&
          & ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
          & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
          & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
          & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
          & pcentPass,qcent,PpreexpfacPass,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
          & Qiprim1,Qiprim2,expA,expB,expC,expD,&
          & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
          & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12Pass,Qdistance12,PQorder,&
          & LocalIntPass1,nLocalIntPass,Acenter,Bcenter,CcenterSpec,DcenterSpec,&
          & nAtomsA,nAtomsB,Spherical,TmpArray1,TMParray1maxsizePass,TmpArray2,&
          & TMParray2maxsizePass,IatomAPass,iatomBPass,UseSP)
     !output private LocalIntPass(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)
     !reorder (including LHS permute) to LocalIntPass(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
     !this can be done on the accelerator
     call MainTriDistributetoLocalIntPass2CPU(TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,&
          & nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD,nContA,nContB,nContC,nContD,&
          & MaxPasses,nPasses,TriangularLHSAtomLoop,Qsegmented,Psegmented,LocalIntPass1,&
          & LocalIntPass2,IatomAPass,iatomBPass,nContQ,nContP)
     
     call TypeDistribution(PermuteLHSTypes,PermuteRHS,nOrbD,nOrbC,startC,startD,&
          & nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startOrbitalB,OutputDim1,OutputDim2,&
          & OutputDim3,OutputDim4,OutputStorage,LocalIntPass2,nOrbQ)

    ENDIF !noscreenCD2
   ENDDO !IatomC
  ENDDO !iAtomD
!$OMP END PARALLEL
  IF(UseGeneralCode)THEN
    call mem_ichor_dealloc(TmpArray3)
    deallocate(TmpArray3)
    call mem_ichor_dealloc(TmpArray4)
    deallocate(TmpArray4)
    call FreeIchorSPHMAT()
  ENDIF
  call mem_ichor_dealloc(TmpArray1)
  deallocate(TmpArray1)
  call mem_ichor_dealloc(TmpArray2)
  deallocate(TmpArray2)    
  call mem_ichor_dealloc(IatomBPass) 
  deallocate(IatomBPass) 
  call mem_ichor_dealloc(IatomAPass) 
  deallocate(IatomAPass) 
  CALL Mem_ichor_dealloc(LocalIntPass1)
  deallocate(LocalIntPass1)
  CALL Mem_ichor_dealloc(LocalIntPass2)
  deallocate(LocalIntPass2)
end subroutine IchorTypeIntegralLoopCPU

subroutine IchorTypeIntegralLoopGPU(nAtomsA,nPrimA,nContA,nOrbCompA,startOrbitalA,&
     & iBatchIndexOfTypeA,expA,ContractCoeffA,AngmomA,Acenter,nBatchA,nOrbA,&
     & nAtomsB,nPrimB,nContB,nOrbCompB,startOrbitalB,&
     & iBatchIndexOfTypeB,expB,ContractCoeffB,AngmomB,Bcenter,nBatchB,nOrbB,&
     & nAtomsC,nPrimC,nContC,nOrbCompC,startOrbitalC,&
     & iBatchIndexOfTypeC,expC,ContractCoeffC,AngmomC,Ccenter,nBatchC,nOrbC,&
     & nAtomsD,nPrimD,nContD,nOrbCompD,startOrbitalD,&
     & iBatchIndexOfTypeD,expD,ContractCoeffD,AngmomD,Dcenter,nBatchD,nOrbD,&
     & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
     & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
     & nPrimP,nPrimQ,nContP,nContQ,expP,expQ,&
     & Ppreexpfac,Qiprim1,Qiprim2,&
     & nTABFJW1,nTABFJW2,TABFJW,TotalAngmom,&
     & CSscreen,noScreenCD2,noScreenAB,THRESHOLD_CS,&
     & Qsegmented,Psegmented,TriangularRHSAtomLoop,TRIANGULARLHSATOMLOOP,&
     & reducedExponents,integralPrefactor,&
     & PcentPass,Pdistance12Pass,PpreExpFacPass,&
     & PQorder,&
     & BATCHGCD,BATCHGAB,Spherical,TMParray1maxsize,nLocalInt,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,OutputStorage,&
     & TMParray2maxsize,PermuteLHSTypes,TriangularODAtomLoop,intprint,lupri,nStaticParamIfac)
  implicit none
  logical,intent(in) :: TriangularRHSAtomLoop,TriangularLHSAtomLoop,PermuteLHSTypes
  logical,intent(in) :: Qsegmented,Psegmented,PQorder,Spherical,CSscreen
  logical,intent(in) :: TriangularODAtomLoop
  integer,intent(in) :: nTABFJW1,nTABFJW2,lupri
  real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2),THRESHOLD_CS
  integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),nStaticParamIfac
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
  real(realk),intent(inout) :: expQ(nPrimQ)
  !P
  integer,intent(in) :: nContP,nPrimP
  real(realk),intent(inout) :: PpreExpFac(nPrimP)
  real(realk),intent(in) :: PcentPass(3,nPrimP,natomsA*natomsB)
  real(realk),intent(in) :: Pdistance12Pass(3,natomsA*natomsB)
  real(realk),intent(in) :: PpreExpFacPass(nPrimP,natomsA*natomsB),expP(nPrimP)
  !
  real(realk),intent(in) :: BATCHGCD(nBatchC,nBatchD)
  !A & B
  integer,intent(in) :: nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV
  integer,intent(in) :: iBatchIndexOfTypeA,iBatchIndexOfTypeB,AngmomB,AngmomA,nOrbA,nOrbB
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbCompA,nOrbCompB,nBatchA,nBatchB,nContA,nContB
  integer,intent(in) :: nPrimA,nPrimB,nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD
  integer,intent(in) :: startOrbitalA(nAtomsA)
  integer,intent(in) :: startOrbitalB(nAtomsB)
  real(realk),intent(in) :: Acenter(3,nAtomsA),expA(nPrimA),ContractCoeffA(nPrimA,nContA)
  real(realk),intent(in) :: Bcenter(3,nAtomsB),expB(nPrimB),ContractCoeffB(nPrimB,nContB)
  real(realk),intent(in) ::  BATCHGAB(nBatchA*nBatchB)
  !collected
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,TotalAngmom
  real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  integer,intent(in) :: TMParray1maxsize,TMParray2maxsize,nLocalInt
  real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP)
  real(realk),intent(in) :: integralPrefactor(nPrimQ,nPrimP)
  !local variables
  integer :: iBatchD,IatomD,IatomC,startC,IatomB,startD,intprint
  real(realk) :: AcenterSpec(3),BcenterSpec(3),GABELM
  logical :: PermuteRHS
  integer :: iOrbQ,iOrbB,iOrbA,iOrbD,iOrbC,I4,I3,I2
  integer :: startA,startB,ndim,nOrbQ,MaxPasses,iAtomCcurr,iAtomDcurr,iCAH
  integer :: TMParray1maxsizePass,TMParray2maxsizePass,nLocalIntPass
  real(realk),allocatable :: LocalIntPass2(:)
#ifdef VAR_OMP
  integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif  
  integer(kind=acckind) :: iSync(maxnAsyncHandles) !The Async handle! (0=free)
  integer(kind=long) :: MaxGPUmemory
  integer(kind=long) :: nSizeStatic,nSizeAsync
  integer :: iAsyncHandles,nAsyncHandles
  !Variables unique for each Async Handle
  integer :: nPasses(maxnAsyncHandles)
  integer,allocatable :: IatomAPass(:,:),IatomBPass(:,:)
  real(realk) :: DcenterSpec(3,maxnAsyncHandles),CcenterSpec(3,maxnAsyncHandles)
  real(realk),allocatable :: TmpArray1(:,:),TmpArray2(:,:)
  real(realk),allocatable :: LocalIntPass1(:,:)
  real(realk),allocatable :: Qcent(:,:,:),Qdistance12(:,:),QpreExpFac(:,:)
! real(realk),intent(inout):: Qcent(3,nPrimQ),Qdistance12(3),QpreExpFac(nPrimQ)
  nLocalIntPass = nOrbA*nAtomsA*nOrbB*nAtomsB*nOrbC*nOrbD
  call DetermineMaxPasses(nAtomsD,iBatchIndexOfTypeD,nAtomsC,nAtomsA,nAtomsB,&
       & iBatchIndexOfTypeC,iBatchIndexOfTypeA,nBatchB,nBatchA,iBatchIndexOfTypeB,&
       & TriangularRHSAtomLoop,CSscreen,TriangularLHSAtomLoop,TriangularODAtomLoop,&
       & noScreenCD2,BATCHGAB,THRESHOLD_CS,noScreenAB,BATCHGCD,nBatchC,nBatchD,MaxPasses)
  TMParray1maxsizePass = TMParray1maxsize*MaxPasses
  TMParray2maxsizePass = TMParray2maxsize*MaxPasses

  !Determine number of Async handles (related to size of memory required)
#ifdef VAR_OPENACC
  MaxGPUmemory = FLOOR(IchorGPUMAXMEM,kind=long)*1000_long*1000_long*1000_long !given by input
  call DeterminenAsyncHandles(nAsyncHandles,MaxGPUmemory,maxnAsyncHandles,nPrimP,&
       & nPrimQ,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD,nTABFJW1,&
       & nTABFJW2,natomsA,natomsB,TMParray1maxsize,TMParray2maxsize,MaxPasses,nLocalIntPass,nStaticParamIfac)
!  WRITE(lupri,*)'nAsyncHandles',nAsyncHandles
  IF(nAsyncHandles.EQ.0)call ichorquit('GPU Memory Error. Calc require too much memory transported to device',-1)
#else
  nAsyncHandles = 1
#endif
  DO iAsyncHandles=1,maxnAsyncHandles
     iSync(iAsyncHandles) = 0          !All Async handles are not in use (set to zero) 
  ENDDO

  ndim = nOrbA*nOrbB*nAtomsA*nAtomsB
  nOrbQ = nOrbC*nOrbD

  allocate(LocalIntPass2(nLocalint*nAtomsA*nAtomsB))
  CALL Mem_ichor_alloc(LocalIntPass2)

! VARIABLES THAT CHANGES FOR EACH ASYNC HANDLE
! SO THAT THEY NEED TO BE UPDATED ON THE DEVICE OR HOST

  allocate(Qcent(3,nPrimQ,nAsyncHandles))
  call mem_ichor_alloc(Qcent)
  allocate(Qdistance12(3,nAsyncHandles))
  call mem_ichor_alloc(Qdistance12)
  allocate(QpreExpFac(nPrimQ,nAsyncHandles))
  call mem_ichor_alloc(QpreExpFac)

  allocate(TmpArray1(TMParray1maxsize*MaxPasses,nAsyncHandles))
  call mem_ichor_alloc(TmpArray1)
  allocate(TmpArray2(TMParray2maxsize*MaxPasses,nAsyncHandles))     
  call mem_ichor_alloc(TmpArray2)
  allocate(IatomAPass(MaxPasses,nAsyncHandles))
  call mem_ichor_alloc(IatomAPass)  
  allocate(IatomBPass(MaxPasses,nAsyncHandles))
  call mem_ichor_alloc(IatomBPass)

  nLocalIntPass = nLocalint*MaxPasses
  allocate(LocalIntPass1(nLocalIntPass,nAsyncHandles))
  CALL Mem_ichor_alloc(LocalIntPass1)
!$ACC DATA COPYIN(nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,&
!$ACC             nPrimQ,nPasses,MaxPasses,intprint,lupri,&
!$ACC             nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,&
!$ACC             ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
!$ACC             nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
!$ACC             nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
!$ACC             nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
!$ACC             pcentPass,Qcent,PpreexpfacPass,Qpreexpfac,&
!$ACC             nTABFJW1,nTABFJW2,TABFJW,Qiprim1,Qiprim2,expA,expB,expC,expD,&
!$ACC             Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
!$ACC             AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12Pass,Qdistance12,PQorder,&
!$ACC             nLocalIntPass,Acenter,Bcenter,&
!$ACC             CcenterSpec,DcenterSpec,nAtomsA,nAtomsB,Spherical,&
!$ACC             TMParray1maxsizePass,TMParray2maxsizePass,&
!$ACC             IatomAPass,iatomBPass) &
!$ACC CREATE(LocalIntPass1) &
!$ACC CREATE(TmpArray1,TmpArray2)
  DO IatomD = 1,nAtomsD
   GABELM = 0.0E0_realk 
   iBatchD = iBatchIndexOfTypeD + IatomD
   DO IatomC = 1,nAtomsC
    IF(TriangularRHSAtomLoop.AND.IatomD.GT.IatomC)CYCLE
    PermuteRHS = TriangularRHSAtomLoop.AND.IatomD.LT.IatomC
    IF(noScreenCD2(IatomC,IatomD))THEN
     IF(CSscreen)GABELM = BATCHGCD(iBatchIndexOfTypeC+IatomC,iBatchD)

     !Find a stream that is available 
     iCAH = -1
     DO iAsyncHandles=1,nAsyncHandles
        IF(iSync(iAsyncHandles).EQ.0_acckind)THEN
           !iAsyncHandles is free 
           iSync(iAsyncHandles) = IatomC*1_acckind + (IatomD-1_acckind)*nAtomsC !unique handle 
           iCAH = iAsyncHandles
           !CurrentAsyncHandles 
           EXIT
        ENDIF
     ENDDO
     IF(iCAH .EQ. -1)call ichorquit('iCAH .EQ. -1',-1)

     nPasses(iCAH) = nAtomsA*nAtomsB
     DcenterSpec(1,iCAH) = Dcenter(1,IAtomD)
     DcenterSpec(2,iCAH) = Dcenter(2,IAtomD)
     DcenterSpec(3,iCAH) = Dcenter(3,IAtomD)
     CcenterSpec(1,iCAH) = Ccenter(1,IAtomC)
     CcenterSpec(2,iCAH) = Ccenter(2,IAtomC)
     CcenterSpec(3,iCAH) = Ccenter(3,IAtomC)
     
     !output: IatomAPass,IatomBPass,nPasses
     CALL BUILD_noScreen2(CSscreen,nAtomsA,nAtomsB,&
          & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,&
          & BATCHGAB,THRESHOLD_CS,GABELM,nPasses(iCAH),&
          & IatomAPass(:,iCAH),IatomBPass(:,iCAH),MaxPasses,&
          & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,IatomD,noScreenAB) 
     IF(nPasses(iCAH).EQ.0)THEN
        iSync(iCAH) = 0 ! setting Async handle free  
        CYCLE
     ENDIF
     IF(nPasses(iCAH).NE.nAtomsA*nAtomsB)THEN
        !$ACC PARALLEL LOOP PRESENT(LocalIntPass1) ASYNC(iSync(iCAH))
        do I4 = 1,nLocalIntPass
           LocalIntPass1(I4,iCAH) = 0.0E0_realk
        enddo
     ENDIF

     !output: Qcent,Qdistance12,QpreExpFac
     CALL Build_qcent_Qdistance12_QpreExpFac(nPrimC,nPrimD,nContC,nContD,&
          & expC,expD,CcenterSpec(:,iCAH),DcenterSpec(:,iCAH),&
          & ContractCoeffC,ContractCoeffD,Qsegmented,&
          & Qcent(:,:,iCAH),Qdistance12(:,iCAH),QpreExpFac(:,iCAH),INTPRINT)


!$ACC UPDATE DEVICE(CcenterSpec(:,iCAH),DcenterSpec(:,iCAH),IatomAPass(:,iCAH),&
!$ACC               iatomBPass(:,iCAH),nPasses(iCAH),Qcent(:,:,iCAH),&
!$ACC               Qpreexpfac(:,iCAH),Qdistance12(:,iCAH)) ASYNC(iSync(iCAH))

     call ICI_GPU_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,&
          & nPrimQ,nPrimP*nPrimQ,nPasses(iCAH),MaxPasses,intprint,lupri,&
          & nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,&
          & ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
          & nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
          & nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,&
          & nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV,&
          & pcentPass,Qcent(:,:,iCAH),PpreexpfacPass,Qpreexpfac(:,iCAH),&
          & nTABFJW1,nTABFJW2,TABFJW,Qiprim1,Qiprim2,expA,expB,expC,expD,&
          & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
          & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12Pass,Qdistance12(:,iCAH),PQorder,&
          & LocalIntPass1(:,iCAH),nLocalIntPass,Acenter,Bcenter,&
          & CcenterSpec(:,iCAH),DcenterSpec(:,iCAH),nAtomsA,nAtomsB,Spherical,&
          & TmpArray1(:,iCAH),TMParray1maxsizePass,TmpArray2(:,iCAH),&
          & TMParray2maxsizePass,IatomAPass(:,iCAH),iatomBPass(:,iCAH),iSync(iCAH),UseSP)
     
!$ACC UPDATE HOST(LocalIntPass1(:,iCAH)) ASYNC(iSync(iCAH))

     !If all iSync(:) are not zero then all streams have been engaged
     !We need to extract result in order to assign that stream new jobs. 
     DO WHILE(MINVAL(iSync(1:nAsyncHandles)).NE.0)
        
        DO iAsyncHandles=1,nAsyncHandles
#ifdef VAR_OPENACC
           IF(acc_async_test(iSync(iAsyncHandles)))THEN
#endif
              !The iAsyncHandles stream is done with the last async task (updating LocalIntPass1)
              iCAH = iAsyncHandles !CurrentAsyncHandles 
              IatomCcurr = mod(iSync(iCAH)-1,nAtomsC)+1
              IatomDcurr = (iSync(iCAH)-1)/nAtomsC + 1
              IF(nPasses(iCAH).NE.nAtomsA*nAtomsB)THEN
                 do I4 = 1,nLocalint*nAtomsA*nAtomsB
                    LocalIntPass2(I4) = 0.0E0_realk
                 enddo
              ENDIF

              !output private LocalIntPass(nContQ,nContP,nPasses,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
              !reorder (including LHS permute) to LocalIntPass(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
              !this can be done on the accelerator
              call MainTriDistributetoLocalIntPass2GPU(TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,&
                   & nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD,nContA,nContB,nContC,nContD,&
                   & MaxPasses,nPasses(iCAH),TriangularLHSAtomLoop,Qsegmented,Psegmented,LocalIntPass1(:,iCAH),&
                   & LocalIntPass2,IatomAPass(:,iCAH),iatomBPass(:,iCAH),nContQ,nContP)
              
              startD = startOrbitalD(iAtomDcurr)
              startC = startOrbitalC(iAtomCcurr)

              call TypeDistribution(PermuteLHSTypes,PermuteRHS,nOrbD,nOrbC,startC,startD,&
                   & nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startOrbitalB,OutputDim1,OutputDim2,&
                   & OutputDim3,OutputDim4,OutputStorage,LocalIntPass2,nOrbQ)
              iSync(iCAH) = 0 !Free for new assignments 
#ifdef VAR_OPENACC
           ENDIF
#endif
        ENDDO
     ENDDO
     
    ENDIF !noscreenCD2
   ENDDO !IatomC
  ENDDO !iAtomD

  !All the IatomC,IatomD have been assigned to streams and started to compute
  !we must wait for all streams to finish - must wait for all iSync to be zero
  DO WHILE(MAXVAL(iSync(1:nAsyncHandles)).NE.0)
        
     DO iAsyncHandles=1,nAsyncHandles
        !Wait untill the iAsyncHandles stream is done with the last async task (updating LocalIntPass1)
        iCAH = iAsyncHandles !CurrentAsyncHandles 
        IF(iSync(iCAH).EQ.0)CYCLE !already done
        !$ACC WAIT(iSync(iCAH))

        IatomCcurr = mod(iSync(iCAH)-1,nAtomsC)+1
        IatomDcurr = (iSync(iCAH)-1)/nAtomsC + 1

        IF(nPasses(iCAH).NE.nAtomsA*nAtomsB)THEN
           do I4 = 1,nLocalint*nAtomsA*nAtomsB
              LocalIntPass2(I4) = 0.0E0_realk
           enddo
        ENDIF
           
        !output private LocalIntPass(nContQ,nContP,nPasses,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
        !reorder (including LHS permute) to LocalIntPass(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
        !this can be done on the accelerator
        call MainTriDistributetoLocalIntPass2GPU(TotalAngmom,nOrbCompA,nOrbCompB,nOrbCompC,&
             & nOrbCompD,nAtomsA,nAtomsB,nOrbA,nOrbB,nOrbC,nOrbD,nContA,nContB,nContC,nContD,&
             & MaxPasses,nPasses(iCAH),TriangularLHSAtomLoop,Qsegmented,Psegmented,LocalIntPass1(:,iCAH),&
             & LocalIntPass2,IatomAPass(:,iCAH),iatomBPass(:,iCAH),nContQ,nContP)
           
        startD = startOrbitalD(iAtomDcurr)
        startC = startOrbitalC(iAtomCcurr)
           
        call TypeDistribution(PermuteLHSTypes,PermuteRHS,nOrbD,nOrbC,startC,startD,&
             & nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startOrbitalB,OutputDim1,OutputDim2,&
             & OutputDim3,OutputDim4,OutputStorage,LocalIntPass2,nOrbQ)
        
        iSync(iCAH) = 0 !Free for new assignments 
     ENDDO
  ENDDO
!$ACC END DATA

  call mem_ichor_dealloc(Qcent)
  deallocate(Qcent)
  call mem_ichor_dealloc(Qdistance12)
  deallocate(Qdistance12)
  call mem_ichor_dealloc(QpreExpFac)
  deallocate(QpreExpFac)

  call mem_ichor_dealloc(TmpArray1)
  deallocate(TmpArray1)
  call mem_ichor_dealloc(TmpArray2)
  deallocate(TmpArray2)    
  call mem_ichor_dealloc(IatomBPass) 
  deallocate(IatomBPass) 
  call mem_ichor_dealloc(IatomAPass) 
  deallocate(IatomAPass) 
  CALL Mem_ichor_dealloc(LocalIntPass1)
  deallocate(LocalIntPass1)
  CALL Mem_ichor_dealloc(LocalIntPass2)
  deallocate(LocalIntPass2)
end subroutine IchorTypeIntegralLoopGPU

!  GPU OpenACC Memory calculation 
subroutine DeterminenAsyncHandles(nAsyncHandles,MaxGPUmemory,maxnAsyncHandles,nPrimP,&
     & nPrimQ,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD,nTABFJW1,&
     & nTABFJW2,natomsA,natomsB,TMParray1maxsize,TMParray2maxsize,MaxPasses,&
     & nLocalIntPass,nStaticParamIfac)
  implicit none
  integer(kind=long),intent(in) :: MaxGPUmemory
  integer,intent(in) :: maxnAsyncHandles,nPrimP,nPrimQ,nPrimA,nContA,nPrimB,nContB
  integer,intent(in) :: nPrimC,nContC,nPrimD,nContD,nTABFJW1,nTABFJW2,natomsA
  integer,intent(in) :: natomsB,TMParray1maxsize,MaxPasses,nLocalIntPass
  integer,intent(in) :: TMParray2maxsize,nStaticParamIfac
  integer,intent(inout) :: nAsyncHandles
  !local variables
  integer(kind=long) :: nSizeStatic,nSizeAsync
  integer :: iAsyncHandles

  !Calculate the size of the memory independent on the number of streams
  !                                                                     
  !Copy of parameters in AGC_TransferRecurrenceParam.F90
  nSizeStatic = nStaticParamIfac*mem_intsize
  ! ACC DATA COPYIN OF                                                
  ! nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,                               
  nSizeStatic = nSizeStatic + mem_intsize*5 
  ! nPrimQ,nPasses,MaxPasses,intprint,lupri,                          
  nSizeStatic = nSizeStatic + mem_intsize*(4 + maxnAsyncHandles)  !nPasses is static allocated
  ! nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,              
  nSizeStatic = nSizeStatic + mem_intsize*6 + mem_realsize*(nPrimP+nPrimQ)
  ! ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,      
  nSizeStatic = nSizeStatic + mem_realsize*(nPrimA*nContA+nPrimB*nContB+nPrimC*nContC+nPrimD*nContD)
  ! nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,                          
  nSizeStatic = nSizeStatic + mem_intsize*4
  ! nCartOrbCompA,nCartOrbCompB,nCartOrbCompC,nCartOrbCompD,          
  nSizeStatic = nSizeStatic + mem_intsize*4
  ! nCartOrbCompP,nCartOrbCompQ,nOrbCompP,nOrbCompQ,nTUVP,nTUVQ,nTUV, 
  nSizeStatic = nSizeStatic + mem_intsize*7
  ! pcentPass,Qcent,PpreexpfacPass,Qpreexpfac,
  nSizeAsync = mem_realsize*4*nPrimQ   !Qcent,Qpreexpfac is nAsyncDependent
  ! nTABFJW1,nTABFJW2,TABFJW,Qiprim1,Qiprim2
  nSizeStatic = nSizeStatic + mem_intsize*(2+2*nPrimQ) + mem_realsize*(nTABFJW1+1)*(nTABFJW2+1) 
  ! expA,expB,expC,expD,     
  nSizeStatic = nSizeStatic + mem_realsize*(nPrimA+nPrimB+nPrimC+nPrimD) 
  ! Qsegmented,Psegmented,reducedExponents,integralPrefactor,         
  nSizeStatic = nSizeStatic + mem_intsize*2 + mem_realsize*nPrimQ*nPrimP
  ! AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12Pass,Qdistance12,PQorder,
  nSizeStatic = nSizeStatic + mem_intsize*5 + mem_realsize*3*natomsA*natomsB
  nSizeAsync = nSizeAsync + mem_realsize*3   !Qdistance12 is nAsyncDependent
  ! nLocalIntPass,Acenter,Bcenter,                                      
  nSizeStatic = nSizeStatic + mem_intsize + mem_realsize*3*(nAtomsA+nAtomsB)
  ! CcenterSpec,DcenterSpec,nAtomsA,nAtomsB,Spherical,                  
  nSizeStatic = nSizeStatic + mem_intsize*3 + mem_realsize*2*3*maxnAsyncHandles
  ! TMParray1maxsizePass,TMParray2maxsizePass,                          
  nSizeStatic = nSizeStatic + mem_intsize*2 + mem_realsize*2*3*maxnAsyncHandles
  ! IatomAPass,iatomBPass                                               
  nSizeAsync = nSizeAsync + mem_intsize*2*MaxPasses   !IatomAPass,iatomBPass is nAsyncDependent
  !ACC CREATE(LocalIntPass1)                                            
  nSizeAsync = nSizeAsync + mem_realsize*nLocalIntPass !LocalIntPass1 is nAsyncDependent
    !ACC CREATE(TmpArray1,TmpArray2)                                      
  nSizeAsync = nSizeAsync + mem_realsize*(TMParray1maxsize*MaxPasses+TMParray2maxsize*MaxPasses)

!  WRITE(lupri,*)'nSizeStatic',nSizeStatic
!  WRITE(lupri,*)'nSizeAsync',nSizeAsync
!  WRITE(lupri,*)'MaxGPUmemory',MaxGPUmemory
  DO iAsyncHandles=1,maxnAsyncHandles
     IF(nSizeStatic+iAsyncHandles*nSizeAsync.LT.MaxGPUmemory)THEN
        nAsyncHandles = iAsyncHandles
     ENDIF
  ENDDO
END subroutine DeterminenAsyncHandles

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
     & PcentPass,PpreExpFacPass,BATCHGCD,BATCHGAB,&
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
  real(realk),intent(inout) :: LocalIntPass(nAtomsA,nAtomsB)
  !Output
  real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5)
  !local variables
  integer :: iBatchD,IatomD,IatomC,startC,IPass,IatomB,IatomA,nPasses,startD
  real(realk) :: DcenterSpec(3),CcenterSpec(3),GABELM
  logical :: PermuteRHS
  integer :: MaxPasses 
  integer :: IatomAPass(natomsA*natomsB)
  integer :: IatomBPass(natomsA*natomsB)
  MaxPasses = natomsA*natomsB
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
          & BATCHGAB,THRESHOLD_CS,GABELM,nPasses,IatomAPass,IatomBPass,MaxPasses,&
          & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,IatomD,noScreenAB) 
     !Determines nPasses = COUNT(noScreenDC)
     IF(nPasses.NE.nAtomsA*nAtomsB)THEN
        call ichorzero(LocalIntPass,SIZE(LocalIntPass))
     ENDIF
     IF(nPasses.EQ.0)CYCLE
     CALL Build_seg_qcent_QpreExpFac(nPrimC,nPrimD,expC,expD,CcenterSpec,DcenterSpec,&
          & ContractCoeffC,ContractCoeffD,Qcent,QpreExpFac,INTPRINT)
     call ICI_seg_seg_SSSS(nPrimP,nPrimQ,nPasses,nAtomsA,nAtomsB,&
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

END MODULE IchorErimod
