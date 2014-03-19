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
MODULE IchorErimodule
  use IchorprecisionModule
  use IchorCommonModule
  use IchorBatchToolsModule
  use IchorEriCoulombintegralCPUOBSGeneralMod, only: IchorCoulombIntegral_CPU_OBS_general, &
       & IchorCoulombIntegral_CPU_OBS_general_size
  use IchorEriCoulombintegralGPUOBSGeneralMod, only: IchorCoulombIntegral_GPU_OBS_general, &
       & IchorCoulombIntegral_GPU_OBS_general_size
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
integer :: i12,offset,K,I,iPrimQ,iPrimP,icont
integer :: oldmaxangmomABCD
integer :: ItypeA,ItypeB,itypeC,itypeD,AngmomA,AngmomB,AngmomC,AngmomD
integer :: ItypeAnon,ItypeBnon,itypeCnon,itypeDnon,nLocalInt
integer :: nPrimA,nPrimB,nContA,nAtomsA,nAtomsB,nAtomsC,nAtomsD,nOrbA
integer :: nDimA,nOrbCompA,nContB,nOrbB,nDimB
integer :: nPrimC,nContC,nPrimD,nContD,nOrbC,nOrbD,nOrbCompB,nOrbCompC,nDimC
integer :: nDimD,nOrbCompD,INTPRINT,maxangmomABCD
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
integer :: BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize
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
logical :: ODscreen,QQRscreen
character(len=16)  :: TYPESTRING
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
call set_ichor_memvar(MaxMemAllocated,MemAllocated)
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

   allocate(MaxGabForTypeAB(nTypesB,nTypesA))
   call mem_ichor_alloc(MaxGabForTypeAB)
   call ObtainMaxGabForType(MaxGabForTypeAB,nTypesA,nTypesB,nAtomsOfTypeA,&
        & nAtomsOfTypeB,BatchIndexOfTypeA,BatchIndexOfTypeB,BATCHGAB,&
        & nBatchA,nBatchB)
   IF(IchorGabID1.EQ.IchorGabID2)THEN
      BATCHGCD => BATCHGAB
      allocate(MaxGabForTypeCD(nTypesC,nTypesD))
      call mem_ichor_alloc(MaxGabForTypeCD)
      MaxGabForTypeCD = MaxGabForTypeAB
      allocate(BatchIndexOfTypeC(nTypesC))
      call mem_ichor_alloc(BatchIndexOfTypeC)
      BatchIndexOfTypeC = BatchIndexOfTypeA
      allocate(BatchIndexOfTypeD(nTypesD))
      call mem_ichor_alloc(BatchIndexOfTypeD)
      BatchIndexOfTypeD = BatchIndexOfTypeB
      nBatchC = nBatchA
      nBatchD = nBatchB
   ELSE
      allocate(BatchIndexOfTypeC(nTypesC))
      call mem_ichor_alloc(BatchIndexOfTypeC)
      call ConstructBatchIndexOfType(BatchIndexOfTypeC,nTypesC,nAtomsOfTypeC,nBatchC)
      allocate(BatchIndexOfTypeD(nTypesD))
      call mem_ichor_alloc(BatchIndexOfTypeD)
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
         call RetrieveGabFromIchorSaveGabModule(nBatchCGCD,nBatchDGCD,IchorGabID1,BATCHGCD2)
         call ExtractBatchGabFromFullGab(nBatchC,nBatchD,BATCHGCD,nBatchCGCD,nBatchDGCD,&
              & BATCHGCD2,startBatchC,endBatchC,startBatchD,endBatchD)
         call mem_ichor_dealloc(BATCHGCD2)
         deallocate(BATCHGCD2)
      ELSE
         call RetrieveGabFromIchorSaveGabModule(nBatchC,nBatchD,IchorGabID2,BATCHGCD)
      ENDIF
      allocate(MaxGabForTypeCD(nTypesC,nTypesD))
      call mem_ichor_alloc(MaxGabForTypeCD)
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
Spherical = SphericalSpec.EQ.SphericalParam
oldmaxangmomABCD = -25

MaxTotalAngmom = MAXVAL(AngmomOfTypeA) + MAXVAL(AngmomOfTypeB) &
     & + MAXVAL(AngmomOfTypeC) + MAXVAL(AngmomOfTypeD)

!we loop over Total angmom in order to ensure that we first do all
!SSSS integrals then PSSS,SPSS,SSPS,SSSP, ...
!this mean we call GAMMATABULATION a limited number of times and
!we reduce branch misprediction inside the code and reuse instruction cache. 
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
   ! DONE TYPE C CALC ===========================
   ! TYPE Q CALC ================================
   sumExtent2CD = (extentC + extentD)*(extentC + extentD)
   
   AngmomQ = AngmomC + AngmomD
   nPrimQ = nPrimC*nPrimD
   nContQ = nContC*nContD
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
      ! DONE TYPE A CALC ===========================
      ! TYPE P CALC ================================
      sumExtent2AB = (extentA + extentB)*(extentA + extentB)
      nPrimP = nPrimA*nPrimB
      nContP = nContA*nContB
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
      ! DONE TYPE P CALC ===========================
      
      ! TYPE PQ CALC ================================
      TotalAngmom = AngmomA + AngmomB + AngmomQ
      maxangmomABCD = AngmomA + AngmomB + AngmomQ
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
      IF(NOTDoSSSS)THEN
         !Determine Sizes of TmpArrays and MaxPasses
         call IchorCoulombIntegral_CPU_OBS_general_size(TMParray1maxsize,&
              & TMParray2maxsize,BasisCont1maxsize,BasisCont2maxsize,&
              & BasisCont3maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
              & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,nContP,&
              & nContQ,nPrimQ*nPrimP,nContQ*nContP,Psegmented,Qsegmented)
         nLocalInt = nOrbA*nOrbB*nOrbC*nOrbD
      ENDIF
      allocate(QpreExpFac(nPrimQ))
      call mem_ichor_alloc(QpreExpFac)
      allocate(Qcent(3*nPrimQ))
      call mem_ichor_alloc(Qcent)
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
      
      !LinK: Form active Dmat 1,2 or 4 matrices (D_BD(nDimB,nDimD),D_AD(nDimA,nDimD),D_BC(nDimB,nDimC),D_AC(nDimA,nDimC))
      allocate(PpreExpFacPass(nPrimP,nAtomsA*nAtomsB))
      call mem_ichor_alloc(PpreExpFacPass)
      allocate(PcentPass(3*nPrimP,nAtomsA*nAtomsB))
      call mem_ichor_alloc(PcentPass)
      IF(CSScreen)THEN
         iBatchIndexOfTypeD = BatchIndexOfTypeD(ItypeD)
         iBatchIndexOfTypeC = BatchIndexOfTypeC(ItypeC)
      ENDIF
      call IchorTimer('START',TSTART,TEND,LUPRI)
      IF(NOTDoSSSS)THEN
         allocate(Pdistance12Pass(3,nAtomsA*nAtomsB))
         call mem_ichor_alloc(Pdistance12Pass)
         CALL Build_pcent_Pdistance12_PpreExpFac(nPrimA,nPrimB,natomsA,natomsB,nContA,nContB,&
              & inversexpP,expA,expB,Acenter,Bcenter,ContractCoeffA,ContractCoeffB,PSegmented,&
              & PcentPass,Pdistance12Pass,PpreExpFacPass,INTPRINT)
         allocate(PpreExpFac(nPrimP))
         call mem_ichor_alloc(PpreExpFac)
!         allocate(Pcent(3*nPrimP))
!         call mem_ichor_alloc(Pcent)
!         WRITE(lupri,*)'CALC TYPE:',ItypeA,ItypeB,ItypeC,ItypeD
!         WRITE(lupri,*)'CALC  Ang:',AngmomA,AngmomB,AngmomC,AngmomD
         call IchorTypeIntegralLoop(nAtomsA,nPrimA,nContA,nOrbCompA,startOrbitalA,&
              & iBatchIndexOfTypeA,expA,ContractCoeffA,AngmomA,Acenter,nBatchA,nOrbA,&
              & nAtomsB,nPrimB,nContB,nOrbCompB,startOrbitalB,&
              & iBatchIndexOfTypeB,expB,ContractCoeffB,AngmomB,Bcenter,nBatchB,nOrbB,&
              & nAtomsC,nPrimC,nContC,nOrbCompC,startOrbitalC,&
              & iBatchIndexOfTypeC,expC,ContractCoeffC,AngmomC,Ccenter,nBatchC,nOrbC,&
              & nAtomsD,nPrimD,nContD,nOrbCompD,startOrbitalD,&
              & iBatchIndexOfTypeD,expD,ContractCoeffD,AngmomD,Dcenter,nBatchD,nOrbD,&
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
              & BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize)

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
      call build_TYPESTRING(TYPESTRING,AngmomA,AngmomB,AngmomC,AngmomD,&
           & nTypesA,nTypesB,nTypesC,nTypesD,ItypeA,ItypeB,ItypeC,ItypeD)
      call IchorTimer(TYPESTRING,TSTART,TEND,LUPRI)
      call mem_ichor_dealloc(Qcent)
      deallocate(Qcent)
      call mem_ichor_dealloc(QpreExpFac)
      deallocate(QpreExpFac)
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
     ENDIF
    ENDDO !typeA
    call mem_ichor_dealloc(expB)
    deallocate(expB)
    call mem_ichor_dealloc(ContractCoeffB)
    deallocate(ContractCoeffB)
    call mem_ichor_dealloc(Bcenter)
    deallocate(Bcenter)
    call mem_ichor_dealloc(StartOrbitalB)
    deallocate(StartOrbitalB)
   ENDDO !typeB
   IF(PermuteRHSTypes)THEN
      !place (ABCD) in (ABDC)
      Call PermuteRHStypesSub(OutputDim1*OutputDim2,OutputDim3,OutputDim4,OutputStorage,&
           & nAtomsC,nAtomsD,startOrbitalC,startOrbitalD,nOrbC,nOrbD,&
           & TriangularRHSAtomLoop,lupri)
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
  ENDDO !typeC
  call mem_ichor_dealloc(expD)
  deallocate(expD)
  call mem_ichor_dealloc(ContractCoeffD)
  deallocate(ContractCoeffD)
  call mem_ichor_dealloc(Dcenter)
  deallocate(Dcenter)
  call mem_ichor_dealloc(StartOrbitalD)
  deallocate(StartOrbitalD)
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
   IF(IchorGabID1.EQ.IchorGabID2)THEN
      !nothing
   ELSE
      call mem_ichor_dealloc(BatchIndexOfTypeC)
      deallocate(BatchIndexOfTypeC)
      call mem_ichor_dealloc(BatchIndexOfTypeD)
      deallocate(BatchIndexOfTypeD)
      call mem_ichor_dealloc(BATCHGCD)
      deallocate(BATCHGCD)
      call mem_ichor_dealloc(MaxGabForTypeCD)
      deallocate(MaxGabForTypeCD)
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
call retrieve_ichor_memvar(MaxMemAllocated,MemAllocated)
IF(INTPRINT.GT.3)THEN
   call stats_ichor_mem(lupri)
ENDIF

end subroutine IchorEri

subroutine build_TYPESTRING(TYPESTRING,AngmomA,AngmomB,AngmomC,AngmomD,&
     & nTypesA,nTypesB,nTypesC,nTypesD,ItypeA,ItypeB,ItypeC,ItypeD)
  implicit none
  character(len=16),intent(inout)  :: TYPESTRING
  integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
  integer,intent(in) :: nTypesA,nTypesB,nTypesC,nTypesD
  integer,intent(in) :: ItypeA,ItypeB,ItypeC,ItypeD
  integer :: n  
  TYPESTRING(1:3) = 'ANG'
  WRITE(TYPESTRING(4:4),'(A1)') AngmomA
  WRITE(TYPESTRING(5:5),'(A1)') AngmomB
  WRITE(TYPESTRING(6:6),'(A1)') AngmomC
  WRITE(TYPESTRING(7:7),'(A1)') AngmomD
  TYPESTRING(8:8) = 'T'
  n=9
  call typeAspec(ItypeA,nTypesA,n,TYPESTRING)
  call typeAspec(ItypeB,nTypesB,n,TYPESTRING)
  call typeAspec(ItypeC,nTypesC,n,TYPESTRING)
  call typeAspec(ItypeD,nTypesD,n,TYPESTRING)
end subroutine build_TYPESTRING

subroutine typeAspec(ItypeA,nTypesA,n,TYPESTRING)
  integer,intent(in) :: ItypeA,nTypesA
  integer,intent(inout) :: n
  character(len=16),intent(inout)  :: TYPESTRING
  IF(nTypesA.LT.10)THEN
     WRITE(TYPESTRING(n:n),'(A1)') ItypeA
     n=n+1
  ELSEIF(nTypesA.LT.100)THEN
     IF(ItypeA.LT.10)THEN
        TYPESTRING(n:n) = '0'
        WRITE(TYPESTRING(n+1:n+1),'(A1)') ItypeA
     ELSE
        WRITE(TYPESTRING(n:n+1),'(A1)') ItypeA
     ENDIF
     n=n+2
  ELSE
     TYPESTRING(n:n+1) = '  '
     n=n+2
  ENDIF
end subroutine typeAspec

subroutine IchorTypeIntegralLoop(nAtomsA,nPrimA,nContA,nOrbCompA,startOrbitalA,&
     & iBatchIndexOfTypeA,expA,ContractCoeffA,AngmomA,Acenter,nBatchA,nOrbA,&
     & nAtomsB,nPrimB,nContB,nOrbCompB,startOrbitalB,&
     & iBatchIndexOfTypeB,expB,ContractCoeffB,AngmomB,Bcenter,nBatchB,nOrbB,&
     & nAtomsC,nPrimC,nContC,nOrbCompC,startOrbitalC,&
     & iBatchIndexOfTypeC,expC,ContractCoeffC,AngmomC,Ccenter,nBatchC,nOrbC,&
     & nAtomsD,nPrimD,nContD,nOrbCompD,startOrbitalD,&
     & iBatchIndexOfTypeD,expD,ContractCoeffD,AngmomD,Dcenter,nBatchD,nOrbD,&
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
     & BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize)
  implicit none
  logical,intent(in) :: TriangularRHSAtomLoop,TriangularLHSAtomLoop,PermuteLHSTypes
  logical,intent(in) :: Qsegmented,Psegmented,PQorder,Spherical,CSscreen
  logical,intent(in) :: TriangularODAtomLoop
  integer,intent(in) :: nTABFJW1,nTABFJW2,lupri
  integer,intent(in) :: BasisCont1maxsize,BasisCont2maxsize,BasisCont3maxsize
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
  real(realk),intent(inout) :: Qcent(3,nPrimQ),Qdistance12(3),QpreExpFac(nPrimQ),expQ(nPrimP)
  !P
  integer,intent(in) :: nContP,nPrimP
  real(realk),intent(inout) :: PpreExpFac(nPrimP)
  real(realk),intent(in) :: PcentPass(3,nPrimP,natomsA*natomsB)
  real(realk),intent(in) :: Pdistance12Pass(3,natomsA*natomsB)
  real(realk),intent(in) :: PpreExpFacPass(nPrimP,natomsA*natomsB),expP(nPrimP)
  !
  real(realk),intent(in) :: BATCHGCD(nBatchC,nBatchD)
  !A & B
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
  real(realk),allocatable :: BasisCont1(:),BasisCont2(:),BasisCont3(:)
  real(realk),allocatable :: LocalIntPass1(:),LocalIntPass2(:)
  integer,allocatable :: IatomAPass(:),IatomBPass(:)
  integer :: iOrbQ,iOrbB,iOrbA,iOrbD,iOrbC,I4,I3,I2
  integer :: startA,startB,ndim,nOrbQ,MaxPasses
  integer :: TMParray1maxsizePass,TMParray2maxsizePass,nLocalIntPass
  integer :: BasisCont1maxsizePass,BasisCont2maxsizePass,BasisCont3maxsizePass
#ifdef VAR_OMP
  integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
  nLocalIntPass = nOrbA*nAtomsA*nOrbB*nAtomsB*nOrbC*nOrbD
  MaxPasses = 0
  DO IatomD = 1,nAtomsD
   GABELM = 0.0E0_realk 
   startD = startOrbitalD(iAtomD)
   iBatchD = iBatchIndexOfTypeD + IatomD
   DO IatomC = 1,nAtomsC
    IF(TriangularRHSAtomLoop.AND.IatomD.GT.IatomC)CYCLE
    PermuteRHS = TriangularRHSAtomLoop.AND.IatomD.LT.IatomC
    startC = startOrbitalC(iAtomC)
    IF(noScreenCD2(IatomC,IatomD))THEN
     nPasses = nAtomsA*nAtomsB
     IF(CSscreen)GABELM = BATCHGCD(iBatchIndexOfTypeC+IatomC,iBatchD)
     !output: IatomAPass,IatomBPass,nPasses
     CALL BUILD_noScreenRed(CSscreen,nAtomsA,nAtomsB,&
          & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,&
          & BATCHGAB,THRESHOLD_CS,GABELM,nPasses,&
          & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,IatomD,noScreenAB)
     MaxPasses = MAX(MaxPasses,nPasses)
    ENDIF
   ENDDO
  ENDDO
  TMParray1maxsizePass = TMParray1maxsize*MaxPasses
  TMParray2maxsizePass = TMParray2maxsize*MaxPasses
  BasisCont1maxsizePass =  BasisCont1maxsize*MaxPasses
  BasisCont2maxsizePass =  BasisCont2maxsize*MaxPasses
  BasisCont3maxsizePass =  BasisCont3maxsize*MaxPasses

  ndim = nOrbA*nOrbB*nAtomsA*nAtomsB
  nOrbQ = nOrbC*nOrbD
  allocate(TmpArray1(TMParray1maxsize*MaxPasses))
  call mem_ichor_alloc(TmpArray1)
  allocate(TmpArray2(TMParray2maxsize*MaxPasses))     
  call mem_ichor_alloc(TmpArray2)
  allocate(BasisCont1(BasisCont1maxsize*MaxPasses))
  call mem_ichor_alloc(BasisCont1)
  allocate(BasisCont2(BasisCont2maxsize*MaxPasses))
  call mem_ichor_alloc(BasisCont2)
  allocate(BasisCont3(BasisCont3maxsize*MaxPasses))
  call mem_ichor_alloc(BasisCont3)
  allocate(IatomAPass(MaxPasses))
  call mem_ichor_alloc(IatomAPass)  
  allocate(IatomBPass(MaxPasses))
  call mem_ichor_alloc(IatomBPass)

  nLocalIntPass = nLocalint*MaxPasses
  allocate(LocalIntPass1(nLocalIntPass))
  CALL Mem_ichor_alloc(LocalIntPass1)
  allocate(LocalIntPass2(nLocalint*nAtomsA*nAtomsB))
  CALL Mem_ichor_alloc(LocalIntPass2)

!!$OMP PARALLEL DEFAULT(none) &
!!$OMP PRIVATE(iAtomD,iAtomC,GABELM,startD,iBatchD,DcenterSpec,PermuteRHS,startC,&
!!$OMP         CcenterSpec,iOrbQ,I3,startB,I4,iOrbD,iOrbC,iAtomB) &
!!$OMP SHARED(nAtomsD,startOrbitalD,iBatchIndexOfTypeD,Dcenter,nAtomsC,&
!!$OMP        TriangularRHSAtomLoop,startOrbitalC,noScreenCD2,Ccenter,&
!!$OMP        nAtomsA,nAtomsB,BATCHGCD,iBatchIndexOfTypeC,CSscreen,&
!!$OMP        nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,&
!!$OMP        BATCHGAB,THRESHOLD_CS,IatomAPass,IatomBPass,startOrbitalB,&
!!$OMP        MaxPasses,TriangularLHSAtomLoop,TriangularODAtomLoop,nOrbB,&
!!$OMP        Qsegmented,nPasses,noScreenAB,nLocalInt,TotalAngmom,nOrbA,&
!!$OMP        nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,intprint,lupri,&
!!$OMP        nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,TABFJW,&
!!$OMP        ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
!!$OMP        pcentPass,qcent,PpreexpfacPass,Qpreexpfac,nTABFJW1,nTABFJW2,&
!!$OMP        Qiprim1,Qiprim2,expA,expB,expC,expD,Psegmented,reducedExponents,&
!!$OMP        integralPrefactor,AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12Pass,&
!!$OMP        Qdistance12,PQorder,LocalIntPass1,LocalIntPass2,nLocalIntPass,&
!!$OMP        Spherical,TmpArray1,TMParray1maxsizePass,TmpArray2,Bcenter,nOrbQ,&
!!$OMP        TMParray2maxsizePass,BasisCont1maxsizePass,BasisCont2maxsizePass,&
!!$OMP        BasisCont3maxsizePass,BasisCont1,BasisCont2,BasisCont3,Acenter,&
!!$OMP        nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,PermuteLHSTypes,nOrbD,nOrbC,&
!!$OMP        startOrbitalA,OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputStorage)
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
!!$OMP SINGLE
     nPasses = nAtomsA*nAtomsB
     CALL BUILD_noScreen2(CSscreen,nAtomsA,nAtomsB,&
          & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,&
          & BATCHGAB,THRESHOLD_CS,GABELM,nPasses,IatomAPass,IatomBPass,MaxPasses,&
          & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,IatomD,noScreenAB) 
!!$OMP END SINGLE
!!$OMP BARRIER
     IF(nPasses.EQ.0)CYCLE
     IF(nPasses.NE.nAtomsA*nAtomsB)THEN
!!$OMP DO 
        do I4 = 1,nLocalIntPass
           LocalIntPass1(I4) = 0.0E0_realk
        enddo
!!$OMP END DO NOWAIT
!!$OMP DO 
        do I4 = 1,nLocalint*nAtomsA*nAtomsB
           LocalIntPass2(I4) = 0.0E0_realk
        enddo
!!$OMP END DO
     ENDIF
!!$OMP SINGLE
     !output: Qcent,Qdistance12,QpreExpFac
     CALL Build_qcent_Qdistance12_QpreExpFac(nPrimC,nPrimD,nContC,nContD,&
          & expC,expD,CcenterSpec,DcenterSpec,ContractCoeffC,ContractCoeffD,&
          & Qsegmented,Qcent,Qdistance12,QpreExpFac,INTPRINT)
!!$OMP END SINGLE
!!$OMP BARRIER
     !Unique for each iPassQ (iAtomC,iAtomD) iteration: qcent,qdistance12,qpreexpfac, Qiprim1(nPrimQ), output:
     !LocalIntPass(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,MaxPasses)
     !IatomAPass,iatomBPass changes and 
!     IF(iAtomC.EQ.1.AND.iAtomD.EQ.1)INTPRINT=1000
     call IchorCoulombIntegral_CPU_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,&
             & nPrimQ,nPrimP*nPrimQ,nPasses,MaxPasses,intprint,lupri,&
             & nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,&
             & ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
             & pcentPass,qcent,PpreexpfacPass,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
             & Qiprim1,Qiprim2,expA,expB,expC,expD,&
             & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
             & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12Pass,Qdistance12,PQorder,&
             & LocalIntPass1,nLocalIntPass,Acenter,Bcenter,CcenterSpec,DcenterSpec,&
             & nAtomsA,nAtomsB,Spherical,TmpArray1,TMParray1maxsizePass,TmpArray2,&
             & TMParray2maxsizePass,BasisCont1maxsizePass,BasisCont2maxsizePass,&
             & BasisCont3maxsizePass,BasisCont1,BasisCont2,BasisCont3,IatomAPass,iatomBPass)
!     IF(iAtomC.EQ.1.AND.iAtomD.EQ.1)INTPRINT=0
     !output private LocalIntPass(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nPasses)
     !reorder (including LHS permute) to LocalIntPass(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD)
     !this can be done on the accelerator
     IF(TriangularLHSAtomLoop)THEN
        IF(Qsegmented.AND.Psegmented)THEN
           IF(TotalAngmom.NE.0)THEN
            call TriDistributeToLocalIntPassSeg(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
              & nOrbCompD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,IatomAPass,iatomBPass,nPasses)
           ELSE !TotalAngmom=0
              !use LocalIntPass (pure copy) but this should not be called
              call ichorquit('DistributeToLocalIntPassSeg0000 not needed - use diff path',-1)
           ENDIF
        ELSE
           IF(TotalAngmom.NE.0)THEN
              call TriDistributeToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
                   & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
                   & MaxPasses,IatomAPass,iatomBPass,nPasses)
           ELSE !TotalAngmom=0
              call TriDistributeToLocalIntPass0000(LocalIntPass1,&
                   & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,&
                   & IatomAPass,iatomBPass,nPasses)
           ENDIF
        ENDIF
     ELSE
        IF(Qsegmented.AND.Psegmented)THEN
           IF(TotalAngmom.NE.0)THEN
              call DistributeToLocalIntPassSeg(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
                   & nOrbCompD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,IatomAPass,iatomBPass,&
                   & nPasses)
           ELSE !TotalAngmom=0
              !use LocalIntPass (pure copy) but this should not be called
              call ichorquit('DistributeToLocalIntPassSeg0000 not needed - use diff path',-1)
           ENDIF
        ELSE
           IF(TotalAngmom.NE.0)THEN
              call DistributeToLocalIntPass(LocalIntPass1,nOrbCompA,nOrbCompB,nOrbCompC,&
                   & nOrbCompD,nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,&
                   & MaxPasses,IatomAPass,iatomBPass,nPasses)
           ELSE !TotalAngmom=0
              call DistributeToLocalIntPass0000(LocalIntPass1,&
                   & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LocalIntPass2,MaxPasses,&
                   & IatomAPass,iatomBPass,nPasses)
           ENDIF
        ENDIF
     ENDIF
!================================================================================
     IF(PermuteLHSTypes)THEN
      IF(PermuteRHS)THEN
!!$OMP DO 
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!$OMP PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB) &
!$OMP SHARED(startOrbitalB,nOrbD,nOrbC,nAtomsB,startD,&
!$OMP        startC,nAtomsA,nOrbA,nOrbB,startOrbitalA,&
!$OMP        OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
!$OMP        OutputStorage,LocalIntPass2,nOrbQ)
       DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
         DO IatomB = 1,nAtomsB
          iOrbQ = iOrbC+(iOrbD-1)*nOrbC
          I4 = startD + iOrbD
          I3 = startC + iOrbC
          startB = startOrbitalB(iAtomB)
          CALL DIST1(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
               & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
               & OutputStorage,LocalIntPass2,nOrbQ,I3,I4,iOrbQ,iatomB)
         ENDDO
        ENDDO
       ENDDO
!$OMP END PARALLEL DO
!!$OMP END DO 
      ELSE  !PermuteLHSTypes NOT PermuteRHS
!!$OMP DO 
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!$OMP PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB) &
!$OMP SHARED(startOrbitalB,nOrbD,nOrbC,nAtomsB,startD,&
!$OMP        startC,nAtomsA,nOrbA,nOrbB,startOrbitalA,&
!$OMP        OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
!$OMP        OutputStorage,LocalIntPass2,nOrbQ)
       DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
         DO IatomB = 1,nAtomsB
          iOrbQ = iOrbC+(iOrbD-1)*nOrbC
          I4 = startD + iOrbD
          I3 = startC + iOrbC
          startB = startOrbitalB(iAtomB)
          CALL DIST2(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
               & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
               & OutputStorage,LocalIntPass2,nOrbQ,I3,I4,iOrbQ,iatomB)
         ENDDO
        ENDDO
       ENDDO
!$OMP END PARALLEL DO
!!$OMP END DO 
      ENDIF
     ELSE !NOT PermuteLHSTypes
      IF(PermuteRHS)THEN
!!$OMP DO 
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!$OMP PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB) &
!$OMP SHARED(startOrbitalB,nOrbD,nOrbC,nAtomsB,startD,&
!$OMP        startC,nAtomsA,nOrbA,nOrbB,startOrbitalA,&
!$OMP        OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
!$OMP        OutputStorage,LocalIntPass2,nOrbQ)
       DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
         DO IatomB = 1,nAtomsB
          iOrbQ = iOrbC+(iOrbD-1)*nOrbC
          I4 = startD + iOrbD
          I3 = startC + iOrbC
          startB = startOrbitalB(iAtomB)
          CALL DIST3(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
               & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
               & OutputStorage,LocalIntPass2,nOrbQ,I3,I4,iOrbQ,iatomB)
         ENDDO
        ENDDO
       ENDDO
!$OMP END PARALLEL DO
!!$OMP END DO 
      ELSE !NOT PermuteLHSTypes NOT PermuteRHS 
!!$OMP DO 
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!$OMP PRIVATE(iOrbD,iOrbC,IatomB,IorbQ,I4,I3,startB) &
!$OMP SHARED(startOrbitalB,nOrbD,nOrbC,nAtomsB,startD,&
!$OMP        startC,nAtomsA,nOrbA,nOrbB,startOrbitalA,&
!$OMP        OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
!$OMP        OutputStorage,LocalIntPass2,nOrbQ)
       DO iOrbD = 1,nOrbD
        DO iOrbC = 1,nOrbC
         DO IatomB = 1,nAtomsB
          iOrbQ = iOrbC+(iOrbD-1)*nOrbC
          I4 = startD + iOrbD
          I3 = startC + iOrbC
          startB = startOrbitalB(iAtomB)
          CALL DIST4(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
               & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
               & OutputStorage,LocalIntPass2,nOrbQ,I3,I4,iOrbQ,iatomB)
         ENDDO
        ENDDO
       ENDDO
!$OMP END PARALLEL DO
!!$OMP END DO 
      ENDIF
     ENDIF
!=========================================================
    ENDIF !noscreenCD2
   ENDDO !IatomC
  ENDDO !iAtomD
!!$OMP END PARALLEL

  call mem_ichor_dealloc(TmpArray1)
  deallocate(TmpArray1)
  call mem_ichor_dealloc(TmpArray2)
  deallocate(TmpArray2)    
  call mem_ichor_dealloc(BasisCont1)
  deallocate(BasisCont1)
  call mem_ichor_dealloc(BasisCont2)
  deallocate(BasisCont2)
  call mem_ichor_dealloc(BasisCont3)
  deallocate(BasisCont3)
  call mem_ichor_dealloc(IatomBPass) 
  deallocate(IatomBPass) 
  call mem_ichor_dealloc(IatomAPass) 
  deallocate(IatomAPass) 
  CALL Mem_ichor_dealloc(LocalIntPass1)
  deallocate(LocalIntPass1)
  CALL Mem_ichor_dealloc(LocalIntPass2)
  deallocate(LocalIntPass2)
end subroutine IchorTypeIntegralLoop

SUBROUTINE BUILD_noScreen2(CSscreen,nAtomsA,nAtomsB,&
     & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,BATCHGAB,&
     & THRESHOLD_CS,GABELM,nPasses,IatomAPass,IatomBPass,MaxPasses,&
     & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,iAtomD,noScreenABin) 
  implicit none
  logical,intent(in) :: CSScreen,TriangularLHSAtomLoop,TriangularODAtomLoop
  integer,intent(in) :: nAtomsA,nAtomsB,iAtomC,iAtomD,MaxPasses
  integer,intent(in) :: iBatchIndexOfTypeA,nBatchB,nBatchA
  integer,intent(in) :: iBatchIndexOfTypeB
  real(realk),intent(in) :: BATCHGAB(nBatchA,nBatchB),THRESHOLD_CS,GABELM
  logical,intent(in) :: noScreenABin(natomsA,natomsB)
  integer,intent(inout) :: nPasses
  integer,intent(inout) :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  !local variables
  integer :: iBatchA,IatomA,iBatchB,IatomB,iPass,IatomAstart,IatomBend
  iPass=0
!  IF(TriangularODAtomLoop)THEN
!     IatomAstart = iAtomC !Restrict AtomC =< AtomA
!  ELSE
     IatomAstart = 1     
!  ENDIF
  IatomBend = nAtomsB
  IF(CSScreen)THEN
   DO IatomA = IatomAstart,nAtomsA
    iBatchA = iBatchIndexOfTypeA + IatomA
    IF(TriangularLHSAtomLoop)IatomBend = IatomA !Restrict AtomB =< AtomA
    DO IatomB = 1,IatomBend
     IF(TriangularODAtomLoop)THEN !If AtomC=AtomA restrict AtomD =< AtomB
      IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
       IF(noScreenABin(IatomA,IatomB))THEN
        IF(GABELM*BATCHGAB(iBatchA,iBatchIndexOfTypeB + IatomB).GT.THRESHOLD_CS)THEN
           iPass = iPass + 1
           IatomAPass(iPass) = IatomA
           IatomBPass(iPass) = IatomB
        ENDIF
       ENDIF
      ENDIF
     ELSE
      IF(noScreenABin(IatomA,IatomB))THEN
       IF(GABELM*BATCHGAB(iBatchA,iBatchIndexOfTypeB + IatomB).GT.THRESHOLD_CS)THEN
          iPass = iPass + 1
          IatomAPass(iPass) = IatomA
          IatomBPass(iPass) = IatomB
       ENDIF
      ENDIF
     ENDIF
    ENDDO
   ENDDO
  ELSE
   DO IatomA = IatomAstart,nAtomsA
    IF(TriangularLHSAtomLoop)IatomBend = IatomA !Restrict AtomB =< AtomA
    DO IatomB = 1,IatomBend
     IF(TriangularODAtomLoop)THEN !If AtomC=AtomA restrict AtomD =< AtomB
      IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
       IF(noScreenABin(IatomA,IatomB))THEN
          iPass = iPass + 1
          IatomAPass(iPass) = IatomA
          IatomBPass(iPass) = IatomB
       ENDIF
      ENDIF
     ELSE
      IF(noScreenABin(IatomA,IatomB))THEN
         iPass = iPass + 1
         IatomAPass(iPass) = IatomA
         IatomBPass(iPass) = IatomB
      ENDIF
     ENDIF
    ENDDO
   ENDDO
  ENDIF
  nPasses = iPass
END SUBROUTINE BUILD_NOSCREEN2

SUBROUTINE BUILD_noScreenRed(CSscreen,nAtomsA,nAtomsB,&
     & nBatchB,nBatchA,iBatchIndexOfTypeA,iBatchIndexOfTypeB,BATCHGAB,&
     & THRESHOLD_CS,GABELM,nPasses,&
     & TriangularLHSAtomLoop,TriangularODAtomLoop,iAtomC,iAtomD,noScreenABin) 
  implicit none
  logical,intent(in) :: CSScreen,TriangularLHSAtomLoop,TriangularODAtomLoop
  integer,intent(in) :: nAtomsA,nAtomsB,iAtomC,iAtomD
  integer,intent(in) :: iBatchIndexOfTypeA,nBatchB,nBatchA
  integer,intent(in) :: iBatchIndexOfTypeB
  real(realk),intent(in) :: BATCHGAB(nBatchA,nBatchB),THRESHOLD_CS,GABELM
  logical,intent(in) :: noScreenABin(natomsA,natomsB)
  integer,intent(inout) :: nPasses
  !local variables
  integer :: iBatchA,IatomA,iBatchB,IatomB,iPass,IatomAstart,IatomBend
  iPass=0
!  IF(TriangularODAtomLoop)THEN
!     IatomAstart = iAtomC !Restrict AtomC =< AtomA
!  ELSE
     IatomAstart = 1     
!  ENDIF
  IatomBend = nAtomsB
  IF(CSScreen)THEN
     DO IatomA = IatomAstart,nAtomsA
        iBatchA = iBatchIndexOfTypeA + IatomA
        IF(TriangularLHSAtomLoop)IatomBend = IatomA !Restrict AtomB =< AtomA
        DO IatomB = 1,IatomBend
         IF(TriangularODAtomLoop)THEN !If AtomC=AtomA restrict AtomD =< AtomB
          IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
           IF(noScreenABin(IatomA,IatomB))THEN
            IF(GABELM*BATCHGAB(iBatchA,iBatchIndexOfTypeB + IatomB).GT.THRESHOLD_CS)THEN
               iPass = iPass + 1
            ENDIF
           ENDIF
          ENDIF
         ELSE
          IF(noScreenABin(IatomA,IatomB))THEN
           IF(GABELM*BATCHGAB(iBatchA,iBatchIndexOfTypeB + IatomB).GT.THRESHOLD_CS)THEN
              iPass = iPass + 1
           ENDIF
          ENDIF
         ENDIF
        ENDDO
     ENDDO
  ELSE
     DO IatomA = IatomAstart,nAtomsA
        IF(TriangularLHSAtomLoop)IatomBend = IatomA !Restrict AtomB =< AtomA
        DO IatomB = 1,IatomBend
         IF(TriangularODAtomLoop)THEN !If AtomC=AtomA restrict AtomD =< AtomB
          IF(IatomA.GT.iAtomC.OR.((IatomA.EQ.iAtomC).AND.(IatomB.GE.IatomD)))THEN
           IF(noScreenABin(IatomA,IatomB))THEN
              iPass = iPass + 1
           ENDIF
          ENDIF
         ELSE
          IF(noScreenABin(IatomA,IatomB))THEN
             iPass = iPass + 1
          ENDIF
         ENDIF
        ENDDO
     ENDDO
  ENDIF
  nPasses = iPass
END SUBROUTINE BUILD_NOSCREENRed

subroutine DIST1(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
     & OutputStorage1,LocalIntPass1,nOrbQ,I3,I4,iOrbQ,iAtomB)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbA,nOrbB,I3,I4,iOrbQ
  integer,intent(in) :: startOrbitalA(nAtomsA),startB,iAtomB,nOrbQ
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4
  real(realk),intent(inout) :: OutputStorage1(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  real(realk),intent(in) :: LocalIntPass1(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbQ)
  !local
  integer :: IatomA,startA,iOrbB,I2,iOrbA
  DO IatomA = 1,nAtomsA
   startA = startOrbitalA(iAtomA)
   DO iOrbB = 1,nOrbB
    I2 = startB + iOrbB
    DO iOrbA = 1,nOrbA
     OutputStorage1(iOrbA + startA,I2,I3,I4)=LocalIntPass1(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     OutputStorage1(iOrbA + startA,I2,I4,I3)=LocalIntPass1(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     OutputStorage1(I2,iOrbA + startA,I3,I4)=LocalIntPass1(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     OutputStorage1(I2,iOrbA + startA,I4,I3)=LocalIntPass1(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
    ENDDO
   ENDDO
  ENDDO
end subroutine DIST1

subroutine DIST2(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
     & OutputStorage2,LocalIntPass2,nOrbQ,I3,I4,iOrbQ,iAtomB)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbA,nOrbB,I3,I4,iOrbQ
  integer,intent(in) :: startOrbitalA(nAtomsA),startB,iAtomB,nOrbQ
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4
  real(realk),intent(inout) :: OutputStorage2(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  real(realk),intent(in) :: LocalIntPass2(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbQ)
  !local
  integer :: IatomA,startA,iOrbB,I2,iOrbA
  DO IatomA = 1,nAtomsA
   startA = startOrbitalA(iAtomA)
   DO iOrbB = 1,nOrbB
    DO iOrbA = 1,nOrbA
     OutputStorage2(iOrbA + startA,iOrbB + startB,I3,I4)=LocalIntPass2(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     OutputStorage2(iOrbB + startB,iOrbA + startA,I3,I4)=LocalIntPass2(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
    ENDDO
   ENDDO
  ENDDO
end subroutine DIST2

subroutine DIST3(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
     & OutputStorage3,LocalIntPass3,nOrbQ,I3,I4,iOrbQ,iAtomB)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbA,nOrbB,I3,I4,iOrbQ
  integer,intent(in) :: startOrbitalA(nAtomsA),startB,iAtomB,nOrbQ
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4
  real(realk),intent(inout) :: OutputStorage3(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  real(realk),intent(in) :: LocalIntPass3(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbQ)
  !local
  integer :: IatomA,startA,iOrbB,I2,iOrbA
  DO IatomA = 1,nAtomsA
   startA = startOrbitalA(iAtomA)
   DO iOrbB = 1,nOrbB
    DO iOrbA = 1,nOrbA
     OutputStorage3(iOrbA + startA,iOrbB + startB,I3,I4)=LocalIntPass3(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     OutputStorage3(iOrbA + startA,iOrbB + startB,I4,I3)=LocalIntPass3(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
    ENDDO
   ENDDO
  ENDDO
end subroutine DIST3

subroutine DIST4(nAtomsA,nAtomsB,nOrbA,nOrbB,startOrbitalA,startB,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,&
     & OutputStorage4,LocalIntPass4,nOrbQ,I3,I4,iOrbQ,iAtomB)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,nOrbA,nOrbB,I3,I4,iOrbQ
  integer,intent(in) :: startOrbitalA(nAtomsA),startB,iAtomB,nOrbQ
  integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4
  real(realk),intent(inout) :: OutputStorage4(OutputDim1,OutputDim2,OutputDim3,OutputDim4)
  real(realk),intent(in) :: LocalIntPass4(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbQ)
  !local
  integer :: IatomA,startA,iOrbB,I2,iOrbA
  DO IatomA = 1,nAtomsA
   startA = startOrbitalA(iAtomA)
   DO iOrbB = 1,nOrbB
    DO iOrbA = 1,nOrbA
     OutputStorage4(iOrbA + startA,iOrbB + startB,I3,I4)=LocalIntPass4(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
    ENDDO
   ENDDO
  ENDDO
end subroutine DIST4

subroutine DistributeIchor1(startA,startB,startC,startD,nOrbD,nOrbC,nOrbB,nOrbA,&
     & OutputStorage,LocalInt,Outputdim1,Outputdim2,Outputdim3,Outputdim4)
  implicit none
  integer,intent(in) :: startA,startB,startC,startD,nOrbA,nOrbB,nOrbC,nOrbD
  integer,intent(in) :: Outputdim1,Outputdim2,Outputdim3,Outputdim4
  real(realk),intent(in) :: LocalInt(nOrbA,nOrbB,nOrbC,nOrbD)
  real(realk),intent(inout) :: OutputStorage(Outputdim1,Outputdim2,Outputdim3,Outputdim4)
  !
  integer :: I2,I3,I4,iOrbA,iOrbB,iOrbC,iOrbD
  DO iOrbD = 1,nOrbD
   I4 = startD + iOrbD
   DO iOrbC = 1,nOrbC
    I3 = startC + iOrbC
    DO iOrbB = 1,nOrbB
     I2 = startB + iOrbB
     DO iOrbA = 1,nOrbA
      OutputStorage(iOrbA + startA,I2,I3,I4)=LocalInt(iOrbA,iOrbB,iOrbC,iOrbD)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
end subroutine DistributeIchor1

subroutine IchorPermuteLHS1(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD,LocalIntPass5,iAtomA,iAtomB,iOrbQ)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB,iAtomA,iAtomB,iOrbQ
  integer,intent(in) :: nOrbA,nOrbB,nOrbC,nOrbD
  real(realk),intent(inout) :: LocalIntPass5(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC*nOrbD)  
  !local variables
  integer :: iOrbB,iOrbA
  DO iOrbB = 1,nOrbB
     DO iOrbA = 1,nOrbA
        LocalIntPass5(iOrbB,iAtomB,iOrbA,iAtomA,iOrbQ)=LocalIntPass5(iOrbA,iAtomA,iOrbB,iAtomB,iOrbQ)
     ENDDO
  ENDDO
end subroutine IchorPermuteLHS1

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

subroutine DistributeToLocalIntPass(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
implicit none 
integer,intent(in)    :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
integer,intent(in)    :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
integer,intent(in)    :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
real(realk),intent(in)::LP1(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContC*nContD,nContA*nContB,MaxPasses)
real(realk),intent(inout)::LP2(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
!local variables
integer :: iContQ,iContA,iContB,iContC,iContD,iContP,offsetA,iAngA,iAngB
integer :: iAngC,iAngD,I3,I4,offsetB,iPass,iAtomA,iAtomB
!!$OMP DO PRIVATE(iContQ,iContA,iContB,iContC,iContD,iContP,offsetA,iAngA,iAngB,iAngC,iAngD,I3,I4,offsetB,iPass,iAtomA,iAtomB)

!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!$OMP PRIVATE(iContQ,iContA,iContB,iContC,iContD,iContP,offsetA,iAngA,iAngB,&
!$OMP         iAngC,iAngD,I3,I4,offsetB,iPass,iAtomA,iAtomB) &
!$OMP SHARED(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB,&
!$OMP        nContA,nContB,nContC,nContD,MaxPasses,nPasses,&
!$OMP        IatomAPass,IatomBPass,LP1,LP2)
  DO IPass = 1,nPasses
   DO iContD = 1,nContD
    DO iContC = 1,nContC
     DO iAngD = 1,nOrbCompD
      DO iAngC = 1,nOrbCompC
       IatomB = IatomBPass(IPass)
       IatomA = IatomAPass(IPass)
       iContQ = iContC + (iContD-1)*nContC
       I4 = iAngD + (iContD-1)*nOrbCompD
       I3 = iAngC + (iContC-1)*nOrbCompC
       DO iContB = 1,nContB
        offsetB = (iContB-1)*nOrbCompB
        DO iContA = 1,nContA
         iContP = iContA+(iContB-1)*nContA
         offsetA = (iContA-1)*nOrbCompA
         DO iAngB = 1,nOrbCompB
          DO iAngA = 1,nOrbCompA
           LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) =&
                & LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
!!$OMP END DO
end subroutine DistributeToLocalIntPass

subroutine TriDistributeToLocalIntPass(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
implicit none 
integer,intent(in)        :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB
integer,intent(in)        :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
real(realk),intent(in)   :: LP1(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContC*nContD,nContA*nContB,MaxPasses)
real(realk),intent(inout):: LP2(nOrbCompA*nContA,nAtomsA,nOrbCompB*nContB,nAtomsB,nOrbCompC*nContC,nOrbCompD*nContD)
!local variables
integer :: iPass,iContP,iContQ,iAngD,iAngC,iAtomA,iAtomB,iContA,iContB
integer :: iContC,iContD,i4,i3,offsetA,offsetB,iAngB,iAngA
 !!!  DO IatomB = 1,nAtomsB
 !!!   DO IatomA = 1,IatomB-1
!!$OMP DO COLLAPSE(3) PRIVATE(iPass,iContP,iContQ,iAngD,iAngC,iAtomA,iAtomB,&
!!$OMP & iContA,iContB,iContC,iContD,i4,i3,offsetA,offsetB,iAngB,iAngA)

!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!$OMP PRIVATE(iPass,iContP,iContQ,iAngD,iAngC,iAtomA,iAtomB,iContA,iContB,iContC,&
!$OMP         iContD,i4,i3,offsetA,offsetB,iAngB,iAngA) &
!$OMP SHARED(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB,&
!$OMP        nContA,nContB,nContC,nContD,MaxPasses,nPasses,IatomAPass,IatomBPass,LP1,LP2)
  DO IPass = 1,nPasses
   DO iContP = 1,nContA*nContA
    DO iContQ = 1,nContC*nContD
     DO iAngD = 1,nOrbCompD
      DO iAngC = 1,nOrbCompC
       IatomB = IatomBPass(IPass)
       IatomA = IatomAPass(IPass)

       !Ordering of Ipass is 
       !iPass = 0 
       !DO IatomA = 1,natomsA
       ! DO IatomB = 1,IatomBend
       !   iPass = iPass + 1
       ! ENDDO
       !ENDDO
       !Where IatomBend=IatomA for triangularLHSatomLoop
       !   or IatomBend=natomsB 

       iContA = iContP - ((iContP-1)/nContA)*nContA
       iContB = (iContP-1)/nContA+1
       iContC = iContQ - ((iContQ-1)/nContC)*nContC
       iContD = (iContQ-1)/nContC+1
       I4 = iAngD + (iContD-1)*nOrbCompD
       I3 = iAngC + (iContC-1)*nOrbCompC
       offsetB = (iContB-1)*nOrbCompA
       offsetA = (iContA-1)*nOrbCompA
       DO iAngB = 1,nOrbCompA
        DO iAngA = 1,nOrbCompA
         LP2(iAngA + offsetA,iatomA,iAngB + offsetB,iatomB,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
         LP2(iAngB + offsetB,iatomB,iAngA + offsetA,iatomA,I3,I4) = LP1(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IPass)
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
!!$OMP END DO
end subroutine TriDistributeToLocalIntPass

subroutine TriDistributeToLocalIntPassSeg(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB,MaxPasses,nPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  real(realk),intent(in)    :: LP1(nOrbCompA,nOrbCompB,nOrbCompC*nOrbCompD,MaxPasses)
  real(realk),intent(inout) :: LP2(nOrbCompA,nAtomsA,nOrbCompB,nAtomsB,nOrbCompC*nOrbCompD)
  !local variables
  integer :: iPass,iAngQ,iAtomA,iAtomB,i4,i3,offsetA,offsetB,iAngB,iAngA
 !!!  DO IatomB = 1,nAtomsB
 !!!   DO IatomA = 1,IatomB-1
!!$OMP DO COLLAPSE(2) PRIVATE(iPass,iAngQ,iAtomA,iAtomB,i4,i3,offsetA,offsetB,iAngB,iAngA)

!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(2) &
!$OMP PRIVATE(iPass,iAngQ,iAtomA,iAtomB,i4,i3,offsetA,offsetB,iAngB,iAngA) &
!$OMP SHARED(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB,MaxPasses,nPasses,&
!$OMP        IatomAPass,IatomBPass,LP1,LP2)
  DO IPass = 1,nPasses
   DO iAngQ = 1,nOrbCompD*nOrbCompC
    IatomB = IatomBPass(IPass)
    IatomA = IatomAPass(IPass)    
    DO iAngB = 1,nOrbCompA
     DO iAngA = 1,nOrbCompA
      LP2(iAngA,iatomA,iAngB,iatomB,iAngQ) = LP1(iAngA,iAngB,iAngQ,IPass)
      LP2(iAngB,iatomB,iAngA,iatomA,iAngQ) = LP1(iAngA,iAngB,iAngQ,IPass)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
!!$OMP END DO
end subroutine TriDistributeToLocalIntPassSeg

subroutine TriDistributeToLocalIntPass0000(LP1,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nAtomsA,nAtomsB
  integer,intent(in)        :: nContA,nContB,nContC,nContD,MaxPasses,nPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  real(realk),intent(in)    :: LP1(nContC*nContD,nContA,nContB,MaxPasses)
  real(realk),intent(inout) :: LP2(nContA,nAtomsA,nContB,nAtomsB,nContC*nContD)
  !local variables
  integer :: iPass,iContP,iContQ,iAtomA,iAtomB,iContA,iContB,iContC,iContD,i4,i3,offsetA,offsetB
 !!!  DO IatomB = 1,nAtomsB
 !!!   DO IatomA = 1,IatomB-1
!!$OMP DO PRIVATE(iPass,iContP,iContQ,iAtomA,iAtomB,iContA,iContB,iContC,iContD,i4,i3,offsetA,offsetB)
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!$OMP PRIVATE(iPass,iContP,iContQ,iAtomA,iAtomB,iContA,iContB, &
!$OMP         iContC,iContD,i4,i3,offsetA,offsetB) &
!$OMP SHARED(nAtomsA,nAtomsB,nContA,nContB,nContC,nContD,MaxPasses,nPasses,&
!$OMP        IatomAPass,IatomBPass,LP1,LP2)
  DO IPass = 1,nPasses
   DO iContB = 1,nContA
    DO iContA = 1,nContA
     IatomB = IatomBPass(IPass)
     IatomA = IatomAPass(IPass)    
     DO iContQ = 1,nContC*nContD
      LP2(iContA,iatomA,iContB,iatomB,iContQ) = LP1(iContQ,iContA,iContB,IPass)
      LP2(iContB,iatomB,iContA,iatomA,iContQ) = LP1(iContQ,iContA,iContB,IPass)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
!!$OMP END DO
end subroutine TriDistributeToLocalIntPass0000

subroutine DistributeToLocalIntPass0000(LP1,&
     & nContA,nContB,nContC,nContD,nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nAtomsA,nAtomsB,nContA,nContB,nContC,nContD,MaxPasses,nPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses)
  real(realk),intent(in)    :: LP1(nContC*nContD,nContA,nContB,MaxPasses)
  real(realk),intent(inout) :: LP2(nContA,nAtomsA,nContB,nAtomsB,nContC*nContD)
  !local variables
  integer :: iContA,iContB,iContQ,ipass,iAtomA,iAtomB
!!$OMP DO PRIVATE(iContA,iContB,iContQ,ipass,iAtomA,iAtomB)
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!$OMP PRIVATE(iContA,iContB,iContQ,ipass,iAtomA,iAtomB) &
!$OMP SHARED(nAtomsA,nAtomsB,nContA,nContB,nContC,nContD,&
!$OMP        MaxPasses,nPasses,IatomAPass,IatomBPass,LP1,LP2)
  DO IPass = 1,nPasses
   DO iContB = 1,nContB
    DO iContA = 1,nContA
     IatomB = IatomBPass(IPass)
     IatomA = IatomAPass(IPass)
     DO iContQ = 1,nContD*nContC
      LP2(iContA,iatomA,iContB,iatomB,iContQ) = LP1(iContQ,iContA,iContB,IPass)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
!!$OMP END DO
end subroutine DistributeToLocalIntPass0000

subroutine DistributeToLocalIntPassSeg(LP1,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,&
     & nAtomsA,nAtomsB,LP2,MaxPasses,IatomAPass,iatomBPass,nPasses)
  implicit none 
  integer,intent(in)        :: nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB,MaxPasses
  integer,intent(in)        :: IatomAPass(MaxPasses),IatomBPass(MaxPasses),nPasses
  real(realk),intent(in)    :: LP1(nOrbCompA,nOrbCompB,nOrbCompC*nOrbCompD,MaxPasses)
  real(realk),intent(inout) :: LP2(nOrbCompA,nAtomsA,nOrbCompB,nAtomsB,nOrbCompC*nOrbCompD)
  !local variables
  integer :: iAngA,iAngB,iAngQ,ipass,iAtomA,iAtomB
!!$OMP DO COLLAPSE(2) PRIVATE(iAngA,iAngB,iAngQ,ipass,iAtomA,iAtomB)
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) &
!$OMP PRIVATE(iAngA,iAngB,iAngQ,ipass,iAtomA,iAtomB) &
!$OMP SHARED(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nAtomsA,nAtomsB,MaxPasses,&
!$OMP        nPasses,IatomAPass,IatomBPass,LP1,LP2)
  DO IPass = 1,nPasses
   DO iAngQ = 1,nOrbCompD*nOrbCompC
    DO iAngB = 1,nOrbCompB
     IatomB = IatomBPass(IPass)
     IatomA = IatomAPass(IPass)
     DO iAngA = 1,nOrbCompA
      LP2(iAngA,iAtomA,iAngB,iAtomB,iAngQ) = LP1(iAngA,iAngB,iAngQ,IPass)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
!!$OMP END DO
end subroutine DistributeToLocalIntPassSeg

subroutine IchorPermuteLHS(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC,nOrbD,LocalIntPass5)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB
  integer,intent(in) :: nOrbA,nOrbB,nOrbC,nOrbD
  real(realk),intent(inout) :: LocalIntPass5(nOrbA,nAtomsA,nOrbB,nAtomsB,nOrbC*nOrbD)  
  !local variables
  integer :: iOrbQ,iAtomB,iOrbB,iatomA,iOrbA
!!$OMP PARALLEL DO DEFAULT(none) PRIVATE(iOrbQ,iatomB,iOrbB,IatomA,&
!!$OMP iOrbA) FIRSTPRIVATE(nAtomsB,nAtomsA,nOrbA,nOrbB,nOrbC,&
!!$OMP nOrbD) SHARED(LocalIntPass5) SCHEDULE(DYNAMIC,1)
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
!!$OMP END PARALLEL DO
end subroutine IchorPermuteLHS

subroutine IchorPermuteLHSSeg0000(nAtomsA,nAtomsB,LocalIntPass5)
  implicit none
  integer,intent(in) :: nAtomsA,nAtomsB
  real(realk),intent(inout) :: LocalIntPass5(nAtomsA,nAtomsB)
  !local variables
  integer :: iAtomB,iatomA
!!$OMP PARALLEL DO DEFAULT(none) PRIVATE(iatomB,&
!!$OMP IatomA) FIRSTPRIVATE(nAtomsB) SHARED(LocalIntPass5) SCHEDULE(DYNAMIC,1)
  DO IatomB = 1,nAtomsB
     DO IatomA = 1,IatomB-1
        LocalIntPass5(iAtomA,iAtomB)=LocalIntPass5(iAtomB,iAtomA)
     ENDDO
  ENDDO
!!$OMP END PARALLEL DO
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
!!$OMP PARALLEL DEFAULT(none) PRIVATE(iOrbD,I4,iOrbC,I3,iatomB,startB,iOrbB,I2,IatomA,&
!!$OMP iOrbA,startA) FIRSTPRIVATE(nOrbD,nOrbC,nOrbB,nOrbA,nAtomsA,nAtomsB,PermuteLHSTypes,&
!!$OMP PermuteRHS) SHARED(startOrbitalB,startOrbitalA,startC,startD,LocalIntPass6,OutputStorage)
  IF(PermuteLHSTypes)THEN
   IF(PermuteRHS)THEN
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!!$OMP DO SCHEDULE(DYNAMIC,1)
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
!!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ELSE  !PermuteLHSTypes NOT PermuteRHS
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!!$OMP DO SCHEDULE(DYNAMIC,1)
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
!!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ENDIF
  ELSE !NOT PermuteLHSTypes
   IF(PermuteRHS)THEN
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!!$OMP DO SCHEDULE(DYNAMIC,1)
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
!!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ELSE !NOT PermuteLHSTypes NOT PermuteRHS 
     DO iOrbD = 1,nOrbD
      I4 = startD + iOrbD
      DO iOrbC = 1,nOrbC
       I3 = startC + iOrbC
!!$OMP DO SCHEDULE(DYNAMIC,1)
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
!!$OMP END DO NOWAIT
      ENDDO
     ENDDO
   ENDIF
  ENDIF
!!$OMP END PARALLEL
 end subroutine IchorDistribute

END MODULE IchorErimodule
