
Subroutine IchorEriInterface(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
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
  use IchorErimodule
  use IchorPrecisionModule
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

call IchorEri(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
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

end Subroutine IchorEriInterface

subroutine IchorGabInterface(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
     & AngmomOfTypeA,nAtomsOfTypeA,nPrimOfTypeA,nContOfTypeA,&
     & startOrbitalOfTypeA,Acenters,exponentsOfTypeA,ContractCoeffOfTypeA,&
     & nTypesB,MaxNatomsB,MaxnPrimB,MaxnContB,&
     & AngmomOfTypeB,nAtomsOfTypeB,nPrimOfTypeB,nContOfTypeB,&
     & startOrbitalOfTypeB,Bcenters,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorInputDim1,&
     & IchorInputDim2,IchorInputDim3,&
     & InputStorage,IchorParSpec,IchorScreenSpec,IchorDebugSpec,&
     & IchorAlgoSpec,SameLHSaos,filestorageIdentifier,MaxMem,&
     & MaxFileStorage,MaxMemAllocated,MemAllocated,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
     & OutputStorage,lupri)
  use IchorGabmodule
  use IchorPrecisionModule
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

call IchorGab(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
     & AngmomOfTypeA,nAtomsOfTypeA,nPrimOfTypeA,nContOfTypeA,&
     & startOrbitalOfTypeA,Acenters,exponentsOfTypeA,ContractCoeffOfTypeA,&
     & nTypesB,MaxNatomsB,MaxnPrimB,MaxnContB,&
     & AngmomOfTypeB,nAtomsOfTypeB,nPrimOfTypeB,nContOfTypeB,&
     & startOrbitalOfTypeB,Bcenters,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorInputDim1,&
     & IchorInputDim2,IchorInputDim3,&
     & InputStorage,IchorParSpec,IchorScreenSpec,IchorDebugSpec,&
     & IchorAlgoSpec,SameLHSaos,filestorageIdentifier,MaxMem,&
     & MaxFileStorage,MaxMemAllocated,MemAllocated,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
     & OutputStorage,lupri)

end subroutine IchorGabInterface

subroutine GET_IchorGabIDInterface(IchorGabID1,IchorGabID2)
  use IchorSaveGabModule
implicit none
integer,intent(inout) :: IchorGabID1,IchorGabID2
call GET_IchorGabID(IchorGabID1,IchorGabID2)
end subroutine GET_IchorGabIDInterface

subroutine SET_IchorGabIDInterface(IchorGabID1,IchorGabID2)
  use IchorSaveGabModule
implicit none
integer,intent(in) :: IchorGabID1,IchorGabID2
call SET_IchorGabID(IchorGabID1,IchorGabID2)
end subroutine SET_IchorGabIDInterface

subroutine GetIchorScreeningParameter(IchorScreenSpec,CSscreen,&
     & ODscreen,QQRscreen)
  use IchorParametersModule
  use IchorCommonModule
  implicit none
  integer,intent(inout) :: IchorScreenSpec
  logical,intent(in)    :: CSscreen,ODscreen,QQRscreen

  IF(CSscreen.AND.ODscreen.AND.QQRscreen)THEN
     IchorScreenSpec = IchorScreen
  ELSEIF(CSscreen.AND.ODscreen)THEN
     IchorScreenSpec = IchorScreenCSOD
  ELSEIF(CSscreen)THEN
     IchorScreenSpec = IchorScreenCS
  ELSEIF(ODscreen)THEN
     IchorScreenSpec = IchorScreenOD
  ELSEIF(QQRscreen)THEN
     call IchorQuit('QQR screening not implemented in ichor',-1)
     IchorScreenSpec = IchorScreenQQR
  ELSE
     IchorScreenSpec = IchorScreenNone
  ENDIF

END subroutine GetIchorScreeningParameter

subroutine GetIchorPermuteParameter(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
  use IchorParametersModule
  implicit none
  integer,intent(inout) :: IchorPermuteSpec
  logical,intent(in)    :: SameLHSaos,SameRHSaos,SameODs
  IF(SameODs.AND.SameRHSaos.AND.SameLHSaos)THEN
     IchorPermuteSpec= IchorPermuteTTT
  ELSEIF(SameRHSaos.AND.SameLHSaos)THEN
     IchorPermuteSpec=IchorPermuteTTF 
  ELSEIF(SameODs)THEN
     IchorPermuteSpec=IchorPermuteFFT 
  ELSEIF(SameLHSaos)THEN
     IchorPermuteSpec=IchorPermuteTFF
  ELSEIF(SameRHSaos)THEN
     IchorPermuteSpec=IchorPermuteFTF
  ELSE
     IchorPermuteSpec=IchorPermuteFFF
  ENDIF
end subroutine GetIchorPermuteParameter

subroutine GetIchorFilestorageIdentifier(filestorageIdentifier)
  use IchorParametersModule
  implicit none
  integer,intent(inout) :: filestorageIdentifier
  filestorageIdentifier = IchorNofilestorage
end subroutine GetIchorFilestorageIdentifier

subroutine GetIchorSphericalParamIdentifier(Identifier)
  use IchorParametersModule
  implicit none
  integer,intent(inout) :: Identifier
  Identifier = SphericalParam
end subroutine GetIchorSphericalParamIdentifier

subroutine GetIchorJobEriIdentifier(Identifier,doLink)
  use IchorParametersModule
  implicit none
  integer,intent(inout) :: Identifier
  logical,intent(in)    :: doLink
  IF(doLink)THEN
     Identifier = IchorJobLink
  ELSE
     Identifier = IchorJobEri
  ENDIF
end subroutine GetIchorJobEriIdentifier

subroutine GetIchorJobMOtransIdentifier(Identifier)
  use IchorParametersModule
  implicit none
  integer,intent(inout) :: Identifier
  Identifier = IchorJobMoTrans
end subroutine GetIchorJobMOtransIdentifier

subroutine GetIchorInputIdentifier(Identifier,rhsDmat)
  use IchorParametersModule
  implicit none
  integer,intent(inout) :: Identifier
  logical,intent(in)    :: rhsDmat
  IF(rhsDmat)THEN
     Identifier = IchorInputDmat
  ELSE
     Identifier = IchorInputNoInput
  ENDIF
end subroutine GetIchorInputIdentifier

subroutine GetIchorParallelSpecIdentifier(Identifier)
  use IchorParametersModule
  implicit none
  integer,intent(inout) :: Identifier
  Identifier = IchorParNone
end subroutine GetIchorParallelSpecIdentifier

subroutine GetIchorDebugIdentifier(Identifier,iprint)
  use IchorParametersModule
  implicit none
  integer,intent(inout) :: Identifier
  integer,intent(in)    :: IPRINT
  Identifier = Iprint
end subroutine GetIchorDebugIdentifier

subroutine GetIchorAlgorithmSpecIdentifier(Identifier)
  use IchorParametersModule
  implicit none
  integer,intent(inout) :: Identifier
  Identifier = IchorAlgoOS
end subroutine GetIchorAlgorithmSpecIdentifier

subroutine FreeIchorSaveGabModuleInterface()
use IchorSaveGabModule
implicit none
call FreeIchorSaveGabModule
end subroutine FreeIchorSaveGabModuleInterface

subroutine InitIchorSaveGabModuleInterface()
use IchorSaveGabModule
implicit none
call InitIchorSaveGabModule
end subroutine InitIchorSaveGabModuleInterface

subroutine AddGabToIchorSaveGabModuleInterface(nBatch1,nBatch2,GabIdentifier,BATCHGAB)
use IchorSaveGabModule
use IchorPrecisionModule
implicit none
integer :: nBatch1,nBatch2,GabIdentifier
real(realk) :: BATCHGAB(nBatch1*nBatch2) 
call AddGabToIchorSaveGabModule(nBatch1,nBatch2,GabIdentifier,BATCHGAB)
end subroutine AddGabToIchorSaveGabModuleInterface

subroutine RetrieveGabFromIchorSaveGabModuleInterface(nBatch1,nBatch2,&
     &GabIdentifier,BATCHGAB)
use IchorSaveGabModule
use IchorPrecisionModule
implicit none
integer :: nBatch1,nBatch2,GabIdentifier
real(realk) :: BATCHGAB(nBatch1*nBatch2) 

call RetrieveGabFromIchorSaveGabModule(nBatch1,nBatch2,&
     &GabIdentifier,BATCHGAB)

end subroutine RetrieveGabFromIchorSaveGabModuleInterface

subroutine RetrieveGabDimFromIchorSaveGabModuleInterface(&
     & nBatch1,nBatch2,GabIdentifier)
use IchorSaveGabModule
implicit none
integer :: nBatch1,nBatch2,GabIdentifier
call RetrieveGabDimFromIchorSaveGabModule(nBatch1,nBatch2,GabIdentifier)
end subroutine RetrieveGabDimFromIchorSaveGabModuleInterface

subroutine InitIchorInputInfo()
  use IchorInputInfoModule
  implicit none
  call InitIchorInputInfoModule
end subroutine InitIchorInputInfo

subroutine FreeIchorInputInfo()
  use IchorInputInfoModule
  implicit none
  call FreeIchorInputInfoModule()
end subroutine FreeIchorInputInfo

subroutine IchorInputM1(M,n1,n2)
  use IchorInputInfoModule
  use IchorPrecisionModule
  implicit none
  integer,intent(in) :: n1,n2
  real(realk),intent(in) :: M(n1,n2)
  call IchorInputInfoM1(M,n1,n2)
end subroutine IchorInputM1

subroutine IchorInputM2(M1,n11,n12,M2,n21,n22)
  use IchorInputInfoModule
  use IchorPrecisionModule
  implicit none
  integer,intent(in) :: n11,n12,n21,n22
  real(realk),intent(in) :: M1(n11,n12),M2(n21,n22)
  call IchorInputInfoM2(M1,n11,n12,M2,n21,n22)
end subroutine IchorInputM2

subroutine IchorInputM3(M1,n11,n12,M2,n21,n22,M3,n31,n32)
  use IchorInputInfoModule
  use IchorPrecisionModule
  implicit none
  integer,intent(in) :: n11,n12,n21,n22,n31,n32
  real(realk),intent(in) :: M1(n11,n12),M2(n21,n22),M3(n31,n32)
  call IchorInputInfoM3(M1,n11,n12,M2,n21,n22,M3,n31,n32)
end subroutine IchorInputM3

subroutine IchorInputM4(M1,n11,n12,M2,n21,n22,M3,n31,n32,M4,n41,n42)
  use IchorInputInfoModule
  use IchorPrecisionModule
  implicit none
  integer,intent(in) :: n11,n12,n21,n22,n31,n32,n41,n42
  real(realk),intent(in) :: M1(n11,n12),M2(n21,n22),M3(n31,n32),M4(n41,n42)
  call IchorInputInfoM4(M1,n11,n12,M2,n21,n22,M3,n31,n32,M4,n41,n42)
end subroutine IchorInputM4

subroutine IchorInputSpec(CenterA,CenterB,CenterC,CenterD)
  use IchorInputInfoModule
  implicit none
  integer :: CenterA,CenterB,CenterC,CenterD
  call IchorInputInfoSpec(CenterA,CenterB,CenterC,CenterD)
end subroutine IchorInputSpec

subroutine WriteCenterInfo1(luoutput,nTypesA,nBatchesA,&
     & MaxnAtomsA,MaxnPrimA,MaxnContA,spherical)
implicit none
integer,intent(in) :: LUOUTPUT !logical unit number of file to write
logical,intent(in) :: spherical
integer,intent(in) :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchesA

WRITE(LUOUTPUT,'(I9)') nTypesA
WRITE(LUOUTPUT,'(I9)') MaxnAtomsA
WRITE(LUOUTPUT,'(I9)') MaxnPrimA
WRITE(LUOUTPUT,'(I9)') MaxnContA
WRITE(LUOUTPUT,'(I9)') nbatchesA
WRITE(LUOUTPUT,'(L1)') spherical

end subroutine WriteCenterInfo1

subroutine WriteCenterInfo2(luoutput,nTypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,MaxnAtomsA,MaxnPrimA,MaxnContA,&
     & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,spherical)
  use IchorPrecisionModule
implicit none
integer,intent(in) :: LUOUTPUT !logical unit number of file to write
logical,intent(in) :: spherical
integer,intent(in) :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchesA
integer,intent(in) :: nAtomsOfTypeA(ntypesA),nPrimOfTypeA(ntypesA) 
integer,intent(in) :: nContOfTypeA(ntypesA),AngmomOfTypeA(ntypesA) 
integer,intent(in) :: startOrbitalOfTypeA(MaxNatomsA,ntypesA)
real(realk),intent(in) :: exponentsOfTypeA(MaxnPrimA,ntypesA),Acenters(3,MaxnAtomsA,nTypesA)
real(realk),intent(in) :: ContractCoeffOfTypeA(MaxnPrimA,MaxnContA,ntypesA)
!
integer :: I,J,K
DO I = 1, ntypesA
   WRITE(LUOUTPUT,'(I9)') nAtomsOfTypeA(I)
ENDDO
DO I = 1, ntypesA
   WRITE(LUOUTPUT,'(I9)') nPrimOfTypeA(I)
ENDDO
DO I = 1, ntypesA
   WRITE(LUOUTPUT,'(I9)') nContOfTypeA(I)
ENDDO
DO I = 1, ntypesA
   WRITE(LUOUTPUT,'(I9)') AngmomOfTypeA(I)
ENDDO
DO J = 1, ntypesA
   DO I = 1, MaxNatomsA
      WRITE(LUOUTPUT,'(I9)') startOrbitalOfTypeA(I,J)
   ENDDO
ENDDO
DO J = 1, ntypesA
   DO I = 1, MaxnPrimA
      WRITE(LUOUTPUT,'(F25.15)') exponentsOfTypeA(I,J)
   ENDDO
ENDDO
DO K = 1, ntypesA
   DO J = 1, MaxnAtomsA
      DO I = 1, 3
         WRITE(LUOUTPUT,'(F25.15)') Acenters(I,J,K)
      ENDDO
   ENDDO
ENDDO
DO K = 1, ntypesA
   DO J = 1, MaxnContA
      DO I = 1, MaxnPrimA
         WRITE(LUOUTPUT,'(F25.15)') ContractCoeffOfTypeA(I,J,K)
      ENDDO
   ENDDO
ENDDO
end subroutine WriteCenterInfo2

subroutine ReadCenterInfo1(luoutput,nTypesA,nBatchesA,&
     & MaxnAtomsA,MaxnPrimA,MaxnContA,spherical)
implicit none
integer,intent(in) :: LUOUTPUT !logical unit number of file to read from
logical,intent(inout) :: spherical
integer,intent(inout) :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchesA

READ(LUOUTPUT,'(I9)') nTypesA
READ(LUOUTPUT,'(I9)') MaxnAtomsA
READ(LUOUTPUT,'(I9)') MaxnPrimA
READ(LUOUTPUT,'(I9)') MaxnContA
READ(LUOUTPUT,'(I9)') nbatchesA
READ(LUOUTPUT,'(L1)') spherical

end subroutine ReadCenterInfo1

subroutine ReadCenterInfo2(luoutput,nTypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,MaxnAtomsA,MaxnPrimA,MaxnContA,&
     & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,spherical)
  use IchorPrecisionModule
implicit none
integer,intent(in) :: LUOUTPUT !logical unit number of file to read from
logical,intent(in) :: spherical
integer,intent(in) :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchesA
integer,intent(inout) :: nAtomsOfTypeA(ntypesA),nPrimOfTypeA(ntypesA) 
integer,intent(inout) :: nContOfTypeA(ntypesA),AngmomOfTypeA(ntypesA) 
integer,intent(inout) :: startOrbitalOfTypeA(MaxNatomsA,ntypesA)
real(realk),intent(inout) :: exponentsOfTypeA(MaxnPrimA,ntypesA),Acenters(3,MaxnAtomsA,nTypesA)
real(realk),intent(inout) :: ContractCoeffOfTypeA(MaxnPrimA,MaxnContA,ntypesA)

integer :: I,J,K
DO I = 1, ntypesA
   READ(LUOUTPUT,'(I9)') nAtomsOfTypeA(I)
ENDDO
DO I = 1, ntypesA
   READ(LUOUTPUT,'(I9)') nPrimOfTypeA(I)
ENDDO
DO I = 1, ntypesA
   READ(LUOUTPUT,'(I9)') nContOfTypeA(I)
ENDDO
DO I = 1, ntypesA
   READ(LUOUTPUT,'(I9)') AngmomOfTypeA(I)
ENDDO
DO J = 1, ntypesA
   DO I = 1, MaxNatomsA
      READ(LUOUTPUT,'(I9)') startOrbitalOfTypeA(I,J)
   ENDDO
ENDDO
DO J = 1, ntypesA
   DO I = 1, MaxnPrimA
      READ(LUOUTPUT,'(F25.15)') exponentsOfTypeA(I,J)
   ENDDO
ENDDO
DO K = 1, ntypesA
   DO J = 1, MaxnAtomsA
      DO I = 1, 3
         READ(LUOUTPUT,'(F25.15)') Acenters(I,J,K)
      ENDDO
   ENDDO
ENDDO
DO K = 1, ntypesA
   DO J = 1, MaxnContA
      DO I = 1, MaxnPrimA
         READ(LUOUTPUT,'(F25.15)') ContractCoeffOfTypeA(I,J,K)
      ENDDO
   ENDDO
ENDDO

end subroutine ReadCenterInfo2
