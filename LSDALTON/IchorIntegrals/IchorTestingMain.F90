PROGRAM IchorErimoduleTEST
implicit none
integer        :: LUPRI,IPRINT
integer        :: dim1,dim2,dim3,dim4 !AO dim
integer        :: OutDim1,OutDim2,OutDim3,OutDim4 !output dim
!
integer                :: nTypes
!A
integer :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchA
integer,pointer     :: nAtomsOfTypeA(:),nPrimOfTypeA(:)
integer,pointer     :: nContOfTypeA(:),AngmomOfTypeA(:),startOrbitalOfTypeA(:,:)
real(8),pointer :: exponentsOfTypeA(:,:),Acenters(:,:,:),ContractCoeffOfTypeA(:,:,:)
!B
integer :: nTypesB,MaxnAtomsB,MaxnPrimB,MaxnContB,nbatchB
integer,pointer     :: nAtomsOfTypeB(:),nPrimOfTypeB(:)
integer,pointer     :: nContOfTypeB(:),AngmomOfTypeB(:),startOrbitalOfTypeB(:,:)
real(8),pointer :: exponentsOfTypeB(:,:),Bcenters(:,:,:),ContractCoeffOfTypeB(:,:,:)
!C
integer :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC,nbatchC
integer,pointer     :: nAtomsOfTypeC(:),nPrimOfTypeC(:)
integer,pointer     :: nContOfTypeC(:),AngmomOfTypeC(:),startOrbitalOfTypeC(:,:)
real(8),pointer :: exponentsOfTypeC(:,:),Ccenters(:,:,:),ContractCoeffOfTypeC(:,:,:)
!D
integer :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD,nbatchD
integer,pointer     :: nAtomsOfTypeD(:),nPrimOfTypeD(:)
integer,pointer     :: nContOfTypeD(:),AngmomOfTypeD(:),startOrbitalOfTypeD(:,:)
real(8),pointer :: exponentsOfTypeD(:,:),Dcenters(:,:,:),ContractCoeffOfTypeD(:,:,:)
!job specification
integer :: SphericalSpec,IchorJobSpec,IchorInputSpec,IchorParSpec,IchorScreenSpec
Integer :: IchorInputDim1,IchorInputDim2,IchorInputDim3,IchorDebugSpec,IchorAlgoSpec
integer :: filestorageIdentifier,MaxFileStorage,IchorPermuteSpec
Integer :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5
Integer :: GabIdentifier, IchorGabID1, IchorGabID2
logical :: SameLHSaos
real(8) :: THRESHOLD_CS,THRESHOLD_QQR,THRESHOLD_OD
real(8),pointer :: InputStorage(:)
real(8),pointer :: BATCHGAB(:),BATCHGCD(:)
Integer(kind=8) :: MaxMem,MaxMemAllocated,MemAllocated
Integer :: nBatchesA,nBatchesB,nBatchesC,nBatchesD,a,b,c,d
logical :: spherical
integer :: nbatchAstart2,nbatchAend2,nbatchBstart2,nbatchBend2
integer :: nbatchCstart2,nbatchCend2,nbatchDstart2,nbatchDend2
logical :: SameRHSaos,SameODs,CRIT1,CRIT2,CRIT3,CRIT4,doLink,rhsDmat,CRIT5
real(8),pointer :: integrals(:,:,:,:)
!FULLABATCH,FULLBBATCH,FULLCBATCH,FULLDBATCH

LUPRI = 6
IPRINT = 0
spherical = .TRUE.
!A
Call BuildCenterAndTypeInfo(1,ntypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,MaxnAtomsA,MaxnPrimA,MaxnContA,&
     & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,&
     & Acenters,spherical,Outdim1)
!B
Call BuildCenterAndTypeInfo(2,ntypesB,nBatchesB,nAtomsOfTypeB,AngmomOfTypeB,&
     & nPrimOfTypeB,nContOfTypeB,MaxnAtomsB,MaxnPrimB,MaxnContB,&
     & startOrbitalOfTypeB,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & Bcenters,spherical,Outdim2)
!C
Call BuildCenterAndTypeInfo(3,ntypesC,nBatchesC,nAtomsOfTypeC,AngmomOfTypeC,&
     & nPrimOfTypeC,nContOfTypeC,MaxnAtomsC,MaxnPrimC,MaxnContC,&
     & startOrbitalOfTypeC,exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters,&
     & spherical,Outdim3)
!D
Call BuildCenterAndTypeInfo(4,ntypesD,nBatchesD,nAtomsOfTypeD,AngmomOfTypeD,&
     & nPrimOfTypeD,nContOfTypeD,MaxnAtomsD,MaxnPrimD,MaxnContD,&
     & startOrbitalOfTypeD,exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters,spherical,Outdim4)

call GetIchorSphericalParamIdentifier(SphericalSpec)
doLink = .FALSE.
call GetIchorJobEriIdentifier(IchorJobSpec,doLink)
rhsDmat = .FALSE. !no rhs density matrix supplied as input 
call GetIchorInputIdentifier(IchorInputSpec,rhsDmat)
IchorInputDim1=1                 !not used since   IcorInputNoInput
IchorInputDim2=1                 !not used since   IcorInputNoInput
IchorInputDim3=1                 !not used since   IcorInputNoInput
allocate(InputStorage(1))
call GetIchorParallelSpecIdentifier(IchorParSpec)   !no parallelization
call GetIchorDebugIdentifier(IchorDebugSpec,iprint) !Debug PrintLevel
call GetIchorAlgorithmSpecIdentifier(IchorAlgoSpec)

SameLHSaos = .FALSE. 
SameRHSaos = .FALSE. 
SameODs = .FALSE. 

call GetIchorPermuteParameter(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
call GetIchorFileStorageIdentifier(filestorageIdentifier)

MaxMem=0         !Maximum Memory Ichor is allowed to use. Zero = no restrictions
MaxFileStorage=0 !Maximum File size, if zero - no file will be written or read. 
MaxMemAllocated=0!Maximum Memory used in the program. Ichor adds to this value
MemAllocated = 0 !Memory allocated in the Ichor program

call GetIchorScreeningParameter(IchorScreenSpec,.FALSE.,&
     & .FALSE.,.FALSE.)
IchorGabID1=0 !screening Matrix Identifier, not used if IchorScreenNone
IchorGabID2=0 !screening Matrix Identifier, not used if IchorScreenNone
allocate(integrals(Outdim1,Outdim2,Outdim3,Outdim4))
OutputDim1=Outdim1
OutputDim2=Outdim2
OutputDim3=Outdim3
OutputDim4=Outdim4
OutputDim5=1
THRESHOLD_OD = 1.0d-10
THRESHOLD_CS = 1.0d-10
THRESHOLD_QQR = 1.0d-10
!print*,'THRESHOLD_CS',THRESHOLD_CS,'THRESHOLD_OD',THRESHOLD_OD
!=====================================================================
!  Main Call
!=====================================================================
call IchorEriInterface(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
     & AngmomOfTypeA,nAtomsOfTypeA,nPrimOfTypeA,nContOfTypeA,&
     & startOrbitalOfTypeA,Acenters,exponentsOfTypeA,ContractCoeffOfTypeA,&
     & nbatchAstart2,nbatchAend2,&
     & nTypesB,MaxNatomsB,MaxnPrimB,MaxnContB,&
     & AngmomOfTypeB,nAtomsOfTypeB,nPrimOfTypeB,nContOfTypeB,&
     & startOrbitalOfTypeB,Bcenters,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & nbatchBstart2,nbatchBend2,&
     & nTypesC,MaxNatomsC,MaxnPrimC,MaxnContC,&
     & AngmomOfTypeC,nAtomsOfTypeC,nPrimOfTypeC,nContOfTypeC,&
     & startOrbitalOfTypeC,Ccenters,exponentsOfTypeC,ContractCoeffOfTypeC,&
     & nbatchCstart2,nbatchCend2,&
     & nTypesD,MaxNatomsD,MaxnPrimD,MaxnContD,&
     & AngmomOfTypeD,nAtomsOfTypeD,nPrimOfTypeD,nContOfTypeD,&
     & startOrbitalOfTypeD,Dcenters,exponentsOfTypeD,ContractCoeffOfTypeD,&
     & nbatchDstart2,nbatchDend2,&
     & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorInputDim1,IchorInputDim2,&
     & IchorInputDim3,&
     & InputStorage,IchorParSpec,IchorScreenSpec,THRESHOLD_OD,THRESHOLD_CS,&
     & THRESHOLD_QQR,IchorGabID1,IchorGabID2,IchorDebugSpec,&
     & IchorAlgoSpec,IchorPermuteSpec,filestorageIdentifier,MaxMem,&
     & MaxFileStorage,MaxMemAllocated,MemAllocated,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
     & integrals,lupri)

do d=1,OutputDim4
   do c=1,OutputDim3
      do b=1,OutputDim2
         do a=1,OutputDim1            
            write(lupri,*)'int',integrals(a,b,c,d)
         enddo
      enddo
   enddo
enddo
!call Mem_Add_external_memory(MaxMemAllocated)
!call mem_dealloc(InputStorage)
!=====================================================================


!=====================================================================
!free space
!call FreeCenterAndTypeInfo(nAtomsOfTypeA,AngmomOfTypeA,&
!           & nPrimOfTypeA,nContOfTypeA,startOrbitalOfTypeA,&
!           & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters)
!call FreeCenterAndTypeInfo(nAtomsOfTypeB,AngmomOfTypeB,&
!           & nPrimOfTypeB,nContOfTypeB,startOrbitalOfTypeB,&
!           & exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters)
!call FreeCenterAndTypeInfo(nAtomsOfTypeC,AngmomOfTypeC,&
!           & nPrimOfTypeC,nContOfTypeC,startOrbitalOfTypeC,&
!           & exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters)
!call FreeCenterAndTypeInfo(nAtomsOfTypeD,AngmomOfTypeD,&
!           & nPrimOfTypeD,nContOfTypeD,startOrbitalOfTypeD,&
!           & exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters)
CONTAINS
Subroutine BuildCenterAndTypeInfo(Center,ntypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,MaxnAtomsA,MaxnPrimA,MaxnContA,&
     & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,spherical,dim)
implicit none
integer,intent(IN)      :: Center
logical,intent(in) :: spherical
integer,intent(inout) :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchesA,dim
integer,pointer     :: nAtomsOfTypeA(:),nPrimOfTypeA(:)      !intent(inout)
integer,pointer     :: nContOfTypeA(:),AngmomOfTypeA(:)      !intent(inout)
integer,pointer     :: startOrbitalOfTypeA(:,:)              !intent(inout)
real(8),pointer :: exponentsOfTypeA(:,:),Acenters(:,:,:)
real(8),pointer :: ContractCoeffOfTypeA(:,:,:)
!
ntypesA = 1           !before(2 of S type 2 of P and 2 of D)
nBatchesA = 2*ntypesA !2 atoms of each
allocate(nAtomsOfTypeA(ntypesA))
allocate(AngmomOfTypeA(ntypesA))
allocate(nPrimOfTypeA(ntypesA))
allocate(nContOfTypeA(ntypesA))

call build_TypeInfo1(nTypesA,&
     & nAtomsOfTypeA,AngmomOfTypeA,nPrimOfTypeA,nContOfTypeA,&
    & MaxnAtomsA,MaxnPrimA,MaxnContA,startOrbitalOfTypeA,&
     & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,dim)

END Subroutine BUILDCENTERANDTYPEINFO

Subroutine build_TypeInfo1(nTypes,nAtomsOfType,&
     & AngmomOfType,nPrimOfType,nContOfType,MaxnAtoms,&
     & MaxnPrim,MaxnCont,startOrbitalOfTypeA,exponentsOfTypeA,&
     & ContractCoeffOfTypeA,CentersOfTypeA,dim)
implicit none
!> the number of different types of batches
INTEGER,intent(in)           :: ntypes
!> the number of atoms for each type of batches
integer,intent(inout)    :: nAtomsOfType(nTypes)
!> the angmom for each type of batches (0 for S, 1 for P,...)
integer,intent(inout)    :: AngmomOfType(nTypes)
!> the number of primitive functions for each type of batches 
integer,intent(inout)    :: nPrimOfType(nTypes)
!> the number of contracted functions for each type of batches 
integer,intent(inout)    :: nContOfType(nTypes)
!> the maximum number of Atoms with a given shell type
integer,intent(inout)    :: MaxnAtoms
!> the maximum number of Primitive in a given shell type
integer,intent(inout)    :: MaxnPrim
!> the maximum number of Contracted functions in a given shell type
integer,intent(inout)    :: MaxnCont,dim
!> the start orbital of the atoms with this type
integer,pointer :: startOrbitalOfTypeA(:,:)
!> the exponents for each type of batches 
real(8),pointer :: exponentsOfTypeA(:,:)
!> the contraction coefficients for each type of batches 
real(8),pointer :: ContractCoeffOfTypeA(:,:,:)
!> the centers of the atoms with this type
real(8),pointer :: CentersOfTypeA(:,:,:)
!
integer :: nBatchType,iBatchType,i,j,k,orbitalindex
integer :: ncol,nOrbComp
nBatchType = ntypes
do iBatchType=1,ntypes 
   nAtomsOfType(iBatchType) = 2 !2
   nPrimOfType(iBatchType) = 2  
   nContOfType(iBatchType) = 1
   AngmomOfType(iBatchType) = 0 !Stype!(iBatchType-1)/2 !0,0,1,1,2,2
enddo
MaxNatoms = 0 
do iBatchType=1,nBatchType
   MaxNatoms = MAX(MaxNatoms,nAtomsOfType(iBatchType))
enddo
MaxnPrim = 0 
do iBatchType=1,nBatchType
   MaxnPrim = MAX(MaxnPrim,nPrimOfType(iBatchType))
enddo
MaxnCont = 0 
do iBatchType=1,nBatchType
   MaxnCont = MAX(MaxnCont,nContOfType(iBatchType))
enddo

allocate(exponentsOfTypeA(MaxnPrim,ntypes))
allocate(startOrbitalOfTypeA(MaxNatoms,ntypes))
allocate(ContractCoeffOfTypeA(MaxnPrim,MaxnCont,ntypes))
allocate(centersOfTypeA(3,MaxnAtoms,nTypes))
do J=1,ntypes
   do I=1,MaxnPrim
      exponentsOfTypeA(I,J) = 1.25d0*I+0.03333d0*J
   enddo
   do K=1,MaxnCont
      do I=1,MaxnPrim
         ContractCoeffOfTypeA(I,K,J) = 0.1d0*I+0.01d0*K+0.001d0*J
      enddo
   enddo
enddo
orbitalindex = 0
do J=1,ntypes
   K = AngmomOfType(J)+1
   ncol = nContOfType(J)
   nOrbComp = 2*K-1
   do I=1,MaxNatoms
      startOrbitalOfTypeA(I,J) = orbitalindex
      orbitalindex = orbitalindex + nOrbComp*ncol
      centersOfTypeA(1,I,J) = 2.25d0*I+0.03333d0*J
      centersOfTypeA(2,I,J) = 0.00001d0*I+0.12345d0*J
      centersOfTypeA(3,I,J) = 0.999d0+I-0.03333d0*J
   enddo
enddo
dim = orbitalindex
end Subroutine Build_TypeInfo1


END PROGRAM IchorErimoduleTEST

