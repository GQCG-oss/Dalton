!> @file
!> Contains the main Ichor screening integral drivers for calculation electron repulsion screening integrals
!> Ichor is the "Integral Code Hand Optimized for Rapid evaluation" 
!> Ichor is the ethereal golden fluid that is the blood of the greek gods

!> \brief Main Ichor drivers for the calculation of screening integrals 
!> based on the Obara Saika(OS)/Head-Gordon-Pople(HGP)
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorGabmod
  use IchorprecisionMod
  use IchorBatchToolsMod
  use IchorCommonMod
  use IchorEriGabintegralOBSGeneralMod, only: IGI_OBS_general, &
       & IGI_OBS_general_size
  use IchorEriCoulombintegralCPUMcMGeneralMod, only: TmpArray3,TmpArray4,&
       & DetermineSizeTmpArray34,precalcichorsphmat,freeichorsphmat,&
       & nTmpArray3,nTmpArray4
  use IchorMemory
  use IchorGammaTabulationMod
  use IchorParametersMod
!debugging
  use IchorEriCoulombintegralCPUOBSGeneralMod, only: ICI_CPU_OBS_general
  use IchorGaussianGeminalMod, only: set_GGem, free_GGem, GGemOperatorCalc
public :: IchorGab
private
CONTAINS

subroutine IchorGab(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
     & AngmomOfTypeA,nAtomsOfTypeA,nPrimOfTypeA,nContOfTypeA,&
     & startOrbitalOfTypeA,Acenters,exponentsOfTypeA,ContractCoeffOfTypeA,&
     & nTypesB,MaxNatomsB,MaxnPrimB,MaxnContB,&
     & AngmomOfTypeB,nAtomsOfTypeB,nPrimOfTypeB,nContOfTypeB,&
     & startOrbitalOfTypeB,Bcenters,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorOperatorSpec,&
     & IchorInputDim1,IchorInputDim2,IchorInputDim3,&
     & InputStorage,IchorParSpec,IchorScreenSpec,IchorDebugSpec,&
     & IchorAlgoSpec,SameLHSaos,filestorageIdentifier,MaxMem,&
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
!
! Same for Center B
!
integer,intent(in) :: nTypesB,MaxnAtomsB,MaxnPrimB,MaxnContB
Integer,intent(in) :: AngmomOfTypeB(ntypesB),nAtomsOfTypeB(ntypesB)
Integer,intent(in) :: nContOfTypeB(ntypesB),nPrimOfTypeB(ntypesB),startOrbitalOfTypeB(MaxNatomsB,ntypesB)
Real(realk),intent(in) :: Bcenters(3,MaxNatomsB,ntypesB),exponentsOfTypeB(MaxnprimB,ntypesB)
Real(realk),intent(in) :: ContractCoeffOfTypeB(MaxnprimB,MaxnContB,ntypesB)
!
!> Spherical Specification (SphericalSpec = SphericalParam = 1) means to use Spherical Harmonic basis functions
Integer,intent(in) :: SphericalSpec
!> Job Specification (IcorJob = IcorJobEri = 1) means that the 4 center 2 electron repulsion integrals
!> should be calculated. 
Integer,intent(in) :: IchorJobSpec
!> Input Specification (IchorInputSpec = IcorInputNoInput = 1) means no Input have been provided
Integer,intent(in) :: IchorInputSpec
!> Operator Specification (IchorOperatorSpec = 1) means Coulomb Operator
Integer,intent(in) :: IchorOperatorSpec
!> Input dimensions assuming InputStorage(IchorInputDim1,IchorInputDim2,IchorInputDim3)
Integer,intent(in) :: IchorInputDim1,IchorInputDim2,IchorInputDim3
!> InputStorage
real(realk),intent(in) :: InputStorage(IchorInputDim1,IchorInputDim2,IchorInputDim3)
!> Parallelization specification (communicator and other stuff should be set up with other call) 
!> IchorParSpec = IchorParNone = 1 means no parallelization - no OpenMP, no MPI, no GPU
Integer,intent(in) :: IchorParSpec
!> Screening specification 
!> IchorScreenSpec = IchorScreen = 1 means default screening including Cauchy-Schwarz screening and QQR
!> IchorScreenSpec = IchorScreenNone = 2 means no screening
Integer,intent(in) :: IchorScreenSpec
!> Debug info specification 
!> IchorDebugSpec = IchorDebugNone = 1 means no debug info
Integer,intent(in) :: IchorDebugSpec
!> Integral Algorithm specification 
!> IchorAlgoSpec = IchorAlgoOS = 1 means Obara-Saika (Head-Gordon Pople)
Integer,intent(in) :: IchorAlgoSpec
!> Permutation specification (SameLHSaos) 
logical,intent(in) :: SameLHSaos
!> Identifier to determine which file should be used to save integrals on disk, This 
!> should be logical unit number, if zero no file is open. 
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
! Local variables
integer :: nPrimP,nContP
integer :: nTABFJW1,nTABFJW2,TotalAngmom
integer :: i12,nPasses,offset,K,I,iPass,MaxPasses,iPrimQ,iPrimP,icont,i2,i1
integer :: ndimPass,oldmaxangmomABCD
integer :: ItypeA,ItypeB,itypeC,itypeD,AngmomA,AngmomB
integer :: ItypeAnon,ItypeBnon,nDimB
integer :: nPrimA,nPrimB,nContA,nAtomsA,nAtomsB,nOrbA,nContB,nOrbCompB
integer :: nDimA,nOrbCompA,iPrimA,iContA,iAtomA,iPrimB,iContB,iAtomB
integer :: nCartOrbCompA,nCartOrbCompB
integer :: INTPRINT,nBatchA,nBatchB,nOrbB
!
integer :: iBatchB,iBatchA
real(realk) :: BcenterSpec(3),AcenterSpec(3),Pdistance12(3),CDAB(1)
!
integer :: iBatchIndexOfTypeB,iBatchIndexOfTypeA,AngmomP
integer :: nOrbCompP,nCartOrbCompP,nTUVP,nTUV
integer,allocatable :: OrderdListA(:),OrderdListB(:)
integer,allocatable :: BatchIndexOfTypeA(:),BatchIndexOfTypeB(:)
logical :: Psegmented,Spherical
logical :: TriangularLHSAtomLoop, PermuteLHS,PQorder
real(realk),allocatable :: expP(:),Pcent(:),PpreExpFac(:)!,CDAB2(:)!,CDAB(:)
REAL(realk),allocatable :: TABFJW(:,:),reducedExponents(:),integralPrefactor(:)
real(realk),allocatable :: expA(:),ContractCoeffA(:,:),Acenter(:,:)
real(realk),allocatable :: expB(:),ContractCoeffB(:,:),Bcenter(:,:)
integer,allocatable :: StartOrbitalA(:),StartOrbitalB(:)
integer :: TMParray1maxsize,TMParray2maxsize,IAngmomTypes,MaxTotalAngmomAB
integer :: BasisContmaxsize
!real(realk) :: SYMFAC,SYMFAC2
real(realk),allocatable :: TmpArray1(:)
real(realk),allocatable :: TmpArray2(:)
INTPRINT=IchorDebugSpec
call set_ichor_memvar(MaxMemAllocated,MemAllocated,MaxMem)
allocate(OrderdListA(nTypesA))
call mem_ichor_alloc(OrderdListA)
call GenerateOrderdListOfTypes(lupri,nTypesA,AngmomOfTypeA,OrderdListA)
allocate(OrderdListB(nTypesB))
call mem_ichor_alloc(OrderdListB)
call GenerateOrderdListOfTypes(lupri,nTypesB,AngmomOfTypeB,OrderdListB)

!TODO 
! Can you set Rpq = 0 always? does it affect Gab matrix?

! GAMMATABULATION 
!     is this needed for (SSSS) ? is it better to build it several times for 
!     different Angmom combis ? 
allocate(BatchIndexOfTypeA(nTypesA))
call mem_ichor_alloc(BatchIndexOfTypeA)
call ConstructBatchIndexOfType(BatchIndexOfTypeA,nTypesA,nAtomsOfTypeA,nBatchA)
IF(nBatchA.NE.OutputDim1)THEN
   CALL ICHORQUIT('Error BATCHGAB dim1 not consistent',-1)
ENDIF

allocate(BatchIndexOfTypeB(nTypesB))
call mem_ichor_alloc(BatchIndexOfTypeB)
call ConstructBatchIndexOfType(BatchIndexOfTypeB,nTypesB,nAtomsOfTypeB,nBatchB)
IF(nBatchB.NE.OutputDim2)THEN
   CALL ICHORQUIT('Error BATCHGAB dim2 not consistent',-1)
ENDIF

IF(SameLHSaos)THEN
   IF(.NOT.(nBatchA.EQ.nBatchB))THEN
      CALL ICHORQUIT('Error IchorGab Screening symmetry error ',-1)
   ENDIF
ENDIF
!call ichorzero2(OutputStorage,OutputDim1*OutputDim2,OutputDim3*OutputDim4*OutputDim5)
!INTRODUCE iAngmomType Loop Like IchorEri.F90 !!!
MaxTotalAngmomAB = MAXVAL(AngmomOfTypeA) + MAXVAL(AngmomOfTypeB)
UseGeneralCode = .FALSE. !Use Specialized code when appropriate. 
call set_GGem(IchorOperatorSpec,2*MaxTotalAngmomAB)
IF(GGemOperatorCalc)UseGeneralCode = .TRUE.

Spherical = SphericalSpec.EQ.SphericalParam
oldmaxangmomABCD = -25
DO IAngmomTypes = 0,MaxTotalAngmomAB
 DO ItypeBnon=1,nTypesB
    call ObtainTypeInfoGab(nTypesB,ItypeBnon,OrderdListB,nAtomsOfTypeB,AngmomOfTypeB,nPrimOfTypeB,&
         & nContOfTypeB,ItypeB,nAtomsB,AngmomB,nPrimB,nContB,nOrbCompB,&
         & nOrbB,nDimB,nCartOrbCompB,spherical)
  IF(nAtomsB.EQ.0)CALL ICHORQUIT('Gab cyle atomB',-1)!CYCLE
  iBatchIndexOfTypeB = BatchIndexOfTypeB(ItypeB)

  allocate(expB(nPrimB))
  call mem_ichor_alloc(expB)
  allocate(ContractCoeffB(nPrimB,nContB))
  call mem_ichor_alloc(ContractCoeffB)
  allocate(Bcenter(3,nAtomsB))
  call mem_ichor_alloc(Bcenter)
  allocate(StartOrbitalB(nAtomsB))
  call mem_ichor_alloc(StartOrbitalB)
  call build_exp_ContractCoeff_center(nPrimB,nContB,nAtomsB,ntypesB,iTypeB,&
       & exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters,StartOrbitalOfTypeB,&
       & expB,ContractCoeffB,Bcenter,StartOrbitalB,MaxnAtomsB,MaxnprimB,MaxnContB,lupri)
  call mem_ichor_dealloc(startOrbitalB) !do not need it
  deallocate(startOrbitalB)
  DO ItypeAnon=1,nTypesA !non ordered loop
   call ObtainTypeInfoGab(nTypesA,ItypeAnon,OrderdListA,nAtomsOfTypeA,AngmomOfTypeA,nPrimOfTypeA,&
        & nContOfTypeA,ItypeA,nAtomsA,AngmomA,nPrimA,nContA,nOrbCompA,&
        & nOrbA,nDimA,nCartOrbCompA,spherical)
   IF(AngmomA+AngmomB.EQ.IAngmomTypes)THEN
    IF(SameLHSaos .AND. ItypeB.GT.ItypeA)CYCLE
    TriangularLHSAtomLoop = SameLHSaos .AND. ItypeB.EQ.ItypeA
    IF(nAtomsA.EQ.0)CALL ICHORQUIT('Gab cyle atomA',-1)!CYCLE
    iBatchIndexOfTypeA = BatchIndexOfTypeA(ItypeA)
    allocate(expA(nPrimA))
    call mem_ichor_alloc(expA)
    allocate(ContractCoeffA(nPrimA,nContA))
    call mem_ichor_alloc(ContractCoeffA)
    allocate(Acenter(3,nAtomsA))
    call mem_ichor_alloc(Acenter)
    allocate(StartOrbitalA(nAtomsA))
    call mem_ichor_alloc(StartOrbitalA)
    call build_exp_ContractCoeff_center(nPrimA,nContA,nAtomsA,ntypesA,iTypeA,&
         & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,StartOrbitalOfTypeA,&
         & expA,ContractCoeffA,Acenter,StartOrbitalA,MaxnAtomsA,MaxnprimA,MaxnContA,lupri) 
    call mem_ichor_dealloc(startOrbitalA) !do not need it
    deallocate(startOrbitalA)
    TotalAngmom = AngmomA + AngmomB + AngmomA + AngmomB 
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
    allocate(expP(nPrimP))
    call mem_ichor_alloc(expP)
    call build_expP(nPrimA,nPrimB,expA,expB,expP)
    IF(TotalAngmom.NE.oldmaxangmomABCD)THEN
       IF(oldmaxangmomABCD.NE.-25)THEN
          call mem_ichor_dealloc(TABFJW)
          deallocate(TABFJW)
       ENDIF
       nTABFJW1 = AngmomA + AngmomB + AngmomA + AngmomB + 3 
       !only need + 3 after Branos change in BUILD_RJ000 
       nTABFJW2 = 1200
       !TABFJW(0:nTABFJW1,0:nTABFJW2)
       allocate(TABFJW(0:nTABFJW1,0:nTABFJW2))
       call mem_ichor_alloc(TABFJW)
       CALL GAMMATABULATION(lupri,TotalAngmom,nTABFJW1,nTABFJW2,TABFJW)  
       oldmaxangmomABCD = TotalAngmom
    ENDIF
    !it may be easy to include primitive screening on the LHS here
    allocate(reducedExponents(nPrimP*nPrimP))
    call mem_ichor_alloc(reducedExponents)
    allocate(integralPrefactor(nPrimP*nPrimP))
    call mem_ichor_alloc(integralPrefactor)
    PQorder = .FALSE.
    call build_reducedExponents_integralPrefactorQP(nPrimP,nPrimP,expP,expP,&
         & reducedExponents,integralPrefactor)
    allocate(pcent(3*nPrimP))
    call mem_ichor_alloc(pcent)
    allocate(PpreExpFac(nPrimP))
    call mem_ichor_alloc(PpreExpFac)
    !       call IchorTimer('START',TSTART,TEND,LUPRI)
    call IGI_OBS_general_size(TMParray1maxsize,&
        & TMParray2maxsize,BasisContmaxsize,AngmomA,AngmomB,&
        & nPrimP,nContP,nPrimB,Psegmented)
    MaxPasses = 1
    !possibly change MaxPasses according to Sizes!
    call GabIntLoop(nPrimA,nPrimB,nPrimP,intprint,lupri,nContA,&
     & nAtomsA,nAtomsB,nContB,nContP,nTABFJW1,nTABFJW2,AngmomA,AngmomB,&
     & TMParray1maxsize,TMParray2maxsize,iBatchIndexOfTypeA,&
     & iBatchIndexOfTypeB,OutputDim1,OutputDim2,expB,ContractCoeffB,&
     & nCartOrbCompA,nCartOrbCompB,nOrbCompA,nOrbCompB,nCartOrbCompP,&
     & nOrbCompP,nTUVP,nTUV,&
     & Bcenter,expA,ContractCoeffA,Acenter,expP,TABFJW,&
     & reducedExponents,integralPrefactor,pcent,PpreExpFac,&
     & OutputStorage,Psegmented,PQorder,Spherical,TriangularLHSAtomLoop,&
     & BasisContmaxsize)
    call mem_ichor_dealloc(pcent)
    deallocate(pcent)
    call mem_ichor_dealloc(PpreExpFac)
    deallocate(PpreExpFac)
    call mem_ichor_dealloc(reducedExponents)
    deallocate(reducedExponents)
    call mem_ichor_dealloc(integralPrefactor)
    deallocate(integralPrefactor)
    call mem_ichor_dealloc(expP)
    deallocate(expP)

    call mem_ichor_dealloc(expA)
    deallocate(expA)
    call mem_ichor_dealloc(ContractCoeffA)
    deallocate(ContractCoeffA)
    call mem_ichor_dealloc(Acenter)
    deallocate(Acenter)
   ENDIF
  ENDDO !typeA
  call mem_ichor_dealloc(expB)
  deallocate(expB)
  call mem_ichor_dealloc(ContractCoeffB)
  deallocate(ContractCoeffB)
  call mem_ichor_dealloc(Bcenter)
  deallocate(Bcenter)
 ENDDO !typeB
ENDDO 
!symmetrize BATCHGAB!
!output requires to zero first otherwise a floating invalid
!WRITE(lupri,*)'The BatchGab NON SYMM'
!call ls_output(OutputStorage,1,nBatchA,1,nBatchB,nBatchA,nBatchB,1,lupri)

IF(SameLHSaos)THEN
   call AddUpperTriAngular(OutputStorage,OutputDim1,lupri)
!   print*,'SameLHSaos GAB OutputStorage'
!   call ls_output(OutputStorage,1,OutputDim1,1,OutputDim1,OutputDim1,OutputDim1,1,6)
ENDIF

call free_GGem()
call mem_ichor_dealloc(TABFJW)
deallocate(TABFJW)
call mem_ichor_dealloc(BatchIndexOfTypeA)
deallocate(BatchIndexOfTypeA)
call mem_ichor_dealloc(BatchIndexOfTypeB)
deallocate(BatchIndexOfTypeB)
call mem_ichor_dealloc(OrderdListA)
deallocate(OrderdListA)
call mem_ichor_dealloc(OrderdListB)
deallocate(OrderdListB)

!WRITE(lupri,*)'The BatchGab'
!call ls_output(OutputStorage,1,nBatchA,1,nBatchB,nBatchA,nBatchB,1,lupri)
!call ls_output(OutputStorage,1,nBatchA,1,nBatchB,nBatchA,nBatchB,1,6)
call retrieve_ichor_memvar(MaxMemAllocated,MemAllocated)
IF(INTPRINT.GT.3)THEN
   call stats_ichor_mem(lupri)
ENDIF

IF(MemAllocated.NE.0)THEN
   call stats_ichor_mem(lupri)
   call ichorquit('MemoryLeak in IchorGab',lupri)
ENDIF
end subroutine IchorGab

subroutine ObtainTypeInfoGab(nTypesD,ItypeDnon,OrderdListD,nAtomsOfTypeD,AngmomOfTypeD,nPrimOfTypeD,&
     & nContOfTypeD,ItypeD,nAtomsD,AngmomD,nPrimD,nContD,nOrbCompD,nOrbD,&
     & nDimD,nCartOrbCompD,spherical)
  implicit none
  logical,intent(in) :: spherical
  integer,intent(in) :: nTypesD,ItypeDnon
  integer,intent(inout) :: nAtomsD,AngmomD,nPrimD,nContD,nCartOrbCompD
  integer,intent(inout) :: nOrbCompD,nOrbD,nDimD,ItypeD
  integer,intent(in) :: OrderdListD(nTypesD),nAtomsOfTypeD(nTypesD),AngmomOfTypeD(nTypesD)
  integer,intent(in) :: nPrimOfTypeD(nTypesD),nContOfTypeD(nTypesD)
  
  ItypeD = OrderdListD(ItypeDnon)  
  nAtomsD = nAtomsOfTypeD(ItypeD)
  AngmomD = AngmomOfTypeD(ItypeD)
  nPrimD = nPrimOfTypeD(ItypeD)
  nContD = nContOfTypeD(ItypeD)
  nCartOrbCompD = (AngmomD+1)*(AngmomD+2)/2
  IF (spherical) THEN
     nOrbCompD = 2*(AngmomD+1)-1
  ELSE
     nOrbCompD = nCartOrbCompD
  ENDIF
  nOrbD = nContD*nOrbCompD
  nDimD = nContD*nOrbCompD*nAtomsD
end subroutine ObtainTypeInfoGab

subroutine GabIntLoop(nPrimA,nPrimB,nPrimP,intprint,lupri,nContA,&
     & nAtomsA,nAtomsB,nContB,nContP,nTABFJW1,nTABFJW2,AngmomA,AngmomB,&
     & TMParray1maxsize,TMParray2maxsize,iBatchIndexOfTypeA,&
     & iBatchIndexOfTypeB,OutputDim1,OutputDim2,expB,ContractCoeffB,&
     & nCartOrbCompA,nCartOrbCompB,nOrbCompA,nOrbCompB,nCartOrbCompP,&
     & nOrbCompP,nTUVP,nTUV,&
     & Bcenter,expA,ContractCoeffA,Acenter,expP,TABFJW,&
     & reducedExponents,integralPrefactor,pcent,PpreExpFac,&
     & OutputStorage,Psegmented,PQorder,Spherical,TriangularLHSAtomLoop,&
     & BasisContmaxsize)
  implicit none
  integer,intent(in) :: nPrimA,nPrimB,nPrimP,intprint,lupri,nContA,nAtomsA
  integer,intent(in) :: nAtomsB,nContB,nContP,nTABFJW1,nTABFJW2,AngmomA,AngmomB,TMParray1maxsize,TMParray2maxsize
  integer,intent(in) :: iBatchIndexOfTypeA,iBatchIndexOfTypeB
  integer,intent(in) :: OutputDim1,OutputDim2,BasisContmaxsize
  integer,intent(in) :: nCartOrbCompA,nCartOrbCompB,nOrbCompA,nOrbCompB,nCartOrbCompP
  integer,intent(in) :: nOrbCompP,nTUVP,nTUV
  real(realk),intent(in) :: expB(nPrimB),ContractCoeffB(nPrimB,nContB),Bcenter(3,nAtomsB)
  real(realk),intent(in) :: expA(nPrimA),ContractCoeffA(nPrimA,nContA),Acenter(3,nAtomsA)
  real(realk),intent(in) :: expP(nPrimP),TABFJW(0:nTABFJW1,0:nTABFJW2)
  real(realk),intent(in) :: reducedExponents(nPrimP*nPrimP)
  real(realk),intent(in) :: integralPrefactor(nPrimP*nPrimP)
  real(realk),intent(inout) :: pcent(3*nPrimP),PpreExpFac(nPrimP)
!  real(realk),intent(inout) :: TmpArray1(TMParray1maxsize)
!  real(realk),intent(inout) :: TmpArray2(TMParray2maxsize)
  real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2)
  logical,intent(in) :: Psegmented,PQorder,Spherical,TriangularLHSAtomLoop
  !local variables
  real(realk) :: CDAB(1),BcenterSpec(3),AcenterSpec(3),Pdistance12(3)
  integer :: IatomB,iBatchA,IatomA,iBatchB,nPasses,MaxPasses,iPass,nPass
  real(realk),allocatable :: TmpArray1(:),TmpArray2(:),BasisCont(:)
  integer,allocatable :: iPassA(:),iPassB(:)
  IF(TriangularLHSAtomLoop)THEN
     allocate(iPassA(nAtomsA*nAtomsB+1/2))
     call mem_ichor_alloc(iPassA)
     allocate(iPassB(nAtomsA*nAtomsB+1/2))
     call mem_ichor_alloc(iPassB)
     iPass = 0
     DO IatomB = 1,nAtomsB
        DO IatomA = IatomB,nAtomsA
           iPass = iPass+1
           iPassA(iPass) = IatomA
           iPassB(iPass) = IatomB
        ENDDO
     ENDDO
     nPass = iPass
  ELSE
     nPass = nAtomsA*nAtomsB
  ENDIF
  allocate(TmpArray1(TMParray1maxsize))
  allocate(TmpArray2(TMParray2maxsize))
  allocate(BasisCont(BasisContmaxsize))
  call mem_ichor_alloc(TmpArray1) 
  call mem_ichor_alloc(TmpArray2) 
  call mem_ichor_alloc(BasisCont) 
  MaxPasses = 1
  nPasses = 1
  IF(MAX(AngmomA,AngmomB).GT.MaxSpecialAngmom.OR.GGemOperatorCalc)THEN
     call DetermineSizeTmpArray34(nTUVP,nCartOrbCompP,nPrimP,nTUVP,&
          & nCartOrbCompP,nPrimP,1,&
          & AngmomA,AngmomB,AngmomA,AngmomB,AngmomA+AngmomB,AngmomA+AngmomB,&
          & AngmomA+AngmomB+AngmomA+AngmomB)
     allocate(TmpArray3(nTmpArray3))
     call mem_ichor_alloc(TmpArray3)
     allocate(TmpArray4(nTmpArray4))
     call mem_ichor_alloc(TmpArray4)
     CALL PreCalciChorSPHMAT(MAX(AngmomA,AngmomB))
  ENDIF
!$OMP PARALLEL DEFAULT(none) &
!$OMP PRIVATE(iPass,iAtomA,iAtomB,iBatchB,BcenterSpec,iBatchA,AcenterSpec) &
!$OMP SHARED(nPrimA,nPrimB,nPrimP,intprint,lupri,nContA,nAtomsA,&
!$OMP        nAtomsB,nContB,nContP,nTABFJW1,nTABFJW2,AngmomA,AngmomB,&
!$OMP        TMParray1maxsize,TMParray2maxsize,iBatchIndexOfTypeA,&
!$OMP        iBatchIndexOfTypeB,OutputDim1,OutputDim2,BasisContmaxsize,&
!$OMP        expB,ContractCoeffB,Bcenter,expA,ContractCoeffA,Acenter,&
!$OMP        expP,TABFJW,reducedExponents,integralPrefactor,nPass,&
!$OMP        pcent,PpreExpFac,OutputStorage,iPassA,iPassB,CDAB,&
!$OMP        Psegmented,PQorder,Spherical,TriangularLHSAtomLoop,&
!$OMP        nCartOrbCompA,nCartOrbCompB,nOrbCompA,nOrbCompB,&
!$OMP        nCartOrbCompP,nOrbCompP,nTUVP,nTUV,&
!$OMP        Pdistance12,nPasses,MaxPasses,TmpArray1,TmpArray2,BasisCont)
  DO IPass = 1,nPass
     IF(TriangularLHSAtomLoop)THEN
        IatomA = iPassA(iPass)
        IatomB = iPassB(iPass)
     ELSE
        IatomA = iPass - ((iPass-1)/nAtomsA)*nAtomsA
        IatomB = (iPass-1)/nAtomsA+1
     ENDIF
     iBatchB = iBatchIndexOfTypeB + IatomB
     BcenterSpec(1) = Bcenter(1,IatomB)
     BcenterSpec(2) = Bcenter(2,IatomB)
     BcenterSpec(3) = Bcenter(3,IatomB)
     iBatchA = iBatchIndexOfTypeA+IatomA
     AcenterSpec(1) = Acenter(1,IatomA)
     AcenterSpec(2) = Acenter(2,IatomA)
     AcenterSpec(3) = Acenter(3,IatomA)
!$OMP SINGLE
     CALL Build_qcent_Qdistance12_QpreExpFac(nPrimA,nPrimB,&
          & nContA,nContB,expA,expB,AcenterSpec,BcenterSpec,ContractCoeffA,&
          & ContractCoeffB,PSegmented,&
          & pcent,Pdistance12,PpreExpFac,INTPRINT)
!$OMP END SINGLE
!$OMP BARRIER
     call IGI_OBS_general(nPrimA,nPrimB,nPrimP,&
          & intprint,lupri,nContA,nContB,nContP,expP,&
          & ContractCoeffA,ContractCoeffB,&
          & nOrbCompA,nOrbCompB,nCartOrbCompA,nCartOrbCompB,&
          & nCartOrbCompP,nOrbCompP,nTUVP,nTUV,&
          & pcent,Ppreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
          & expA,expB,Psegmented,reducedExponents,integralPrefactor,&
          & AngmomA,AngmomB,Pdistance12,PQorder,CDAB,&
          & AcenterSpec,BcenterSpec,Spherical,&
          & TmpArray1,TMParray1maxsize,TmpArray2,TMParray2maxsize,&
          & BasisContmaxsize,BasisCont)
!$OMP SINGLE
     OutputStorage(iBatchA,iBatchB) = CDAB(1)
!$OMP END SINGLE
  ENDDO
!$OMP END PARALLEL
  IF(MAX(AngmomA,AngmomB).GT.MaxSpecialAngmom.OR.GGemOperatorCalc)THEN
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
  call mem_ichor_dealloc(BasisCont)
  deallocate(BasisCont)
  IF(TriangularLHSAtomLoop)THEN
     call mem_ichor_dealloc(iPassA)
     deallocate(iPassA)
     call mem_ichor_dealloc(iPassB)
     deallocate(iPassB)
  ENDIF
end subroutine GabIntLoop

subroutine AddUpperTriAngular(MAT,nBatchA,lupri)
implicit none
integer,intent(in) :: nBatchA,lupri
real(realk),intent(inout) :: MAT(nBatchA,nBatchA)
!
integer :: IA,IB
DO IA = 1,nBatchA
   DO IB = 1,IA-1
      MAT(IB,IA) = MAT(IA,IB)
   ENDDO
ENDDO
end subroutine AddUpperTriAngular
END MODULE IchorGabmod
