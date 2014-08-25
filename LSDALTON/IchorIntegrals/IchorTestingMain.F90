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
integer :: nbatchCstart2,nbatchCend2,nbatchDstart2,nbatchDend2,LUOUTPUT
logical :: SameRHSaos,SameODs,CRIT1,CRIT2,CRIT3,CRIT4,doLink,rhsDmat,CRIT5,FAIL
real(8),pointer :: integrals(:,:,:,:),ComparisonInt(:,:,:,:)
character(len=100) :: filename
do A=1,100
   filename(A:A) = ' '
enddo

!FULLABATCH,FULLBBATCH,FULLCBATCH,FULLDBATCH

LUPRI = 6
IPRINT = 0
spherical = .TRUE.
FileName = 'IchorUnitTestUnitTest_segD1pUnitTest_segDUnitTest_genDUnitTest_segD1p'

LUOUTPUT = 12
print*,'FileName',FileName
open(unit = LUOUTPUT, file=TRIM(FileName),status='OLD',FORM='FORMATTED')

!A
call ReadCenterInfo1(luoutput,nTypesA,nBatchesA,MaxnAtomsA,MaxnPrimA,MaxnContA,spherical)
allocate(nAtomsOfTypeA(ntypesA))
allocate(AngmomOfTypeA(ntypesA))
allocate(nPrimOfTypeA(ntypesA))
allocate(nContOfTypeA(ntypesA))
allocate(startOrbitalOfTypeA(MaxNatomsA,ntypesA))
allocate(exponentsOfTypeA(MaxnPrimA,ntypesA))
allocate(ContractCoeffOfTypeA(MaxnPrimA,MaxnContA,ntypesA))
allocate(Acenters(3,MaxnAtomsA,nTypesA))
Call ReadCenterInfo2(luoutput,ntypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,MaxnAtomsA,MaxnPrimA,MaxnContA,&
     & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,&
     & Acenters,spherical)
!B
call ReadCenterInfo1(luoutput,nTypesB,nBatchesB,MaxnAtomsB,MaxnPrimB,MaxnContB,spherical)
allocate(nAtomsOfTypeB(ntypesB))
allocate(AngmomOfTypeB(ntypesB))
allocate(nPrimOfTypeB(ntypesB))
allocate(nContOfTypeB(ntypesB))
allocate(startOrbitalOfTypeB(MaxNatomsB,ntypesB))
allocate(exponentsOfTypeB(MaxnPrimB,ntypesB))
allocate(ContractCoeffOfTypeB(MaxnPrimB,MaxnContB,ntypesB))
allocate(Bcenters(3,MaxnAtomsB,nTypesB))
Call ReadCenterInfo2(luoutput,ntypesB,nBatchesB,nAtomsOfTypeB,AngmomOfTypeB,&
     & nPrimOfTypeB,nContOfTypeB,MaxnAtomsB,MaxnPrimB,MaxnContB,&
     & startOrbitalOfTypeB,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & Bcenters,spherical)
!C
call ReadCenterInfo1(luoutput,nTypesC,nBatchesC,MaxnAtomsC,MaxnPrimC,MaxnContC,spherical)
allocate(nAtomsOfTypeC(ntypesC))
allocate(AngmomOfTypeC(ntypesC))
allocate(nPrimOfTypeC(ntypesC))
allocate(nContOfTypeC(ntypesC))
allocate(startOrbitalOfTypeC(MaxNatomsC,ntypesC))
allocate(exponentsOfTypeC(MaxnPrimC,ntypesC))
allocate(ContractCoeffOfTypeC(MaxnPrimC,MaxnContC,ntypesC))
allocate(Ccenters(3,MaxnAtomsC,nTypesC))
Call ReadCenterInfo2(luoutput,ntypesC,nBatchesC,nAtomsOfTypeC,AngmomOfTypeC,&
     & nPrimOfTypeC,nContOfTypeC,MaxnAtomsC,MaxnPrimC,MaxnContC,&
     & startOrbitalOfTypeC,exponentsOfTypeC,ContractCoeffOfTypeC,&
     & Ccenters,spherical)
!D
call ReadCenterInfo1(luoutput,nTypesD,nBatchesD,MaxnAtomsD,MaxnPrimD,MaxnContD,spherical)
allocate(nAtomsOfTypeD(ntypesD))
allocate(AngmomOfTypeD(ntypesD))
allocate(nPrimOfTypeD(ntypesD))
allocate(nContOfTypeD(ntypesD))
allocate(startOrbitalOfTypeD(MaxNatomsD,ntypesD))
allocate(exponentsOfTypeD(MaxnPrimD,ntypesD))
allocate(ContractCoeffOfTypeD(MaxnPrimD,MaxnContD,ntypesD))
allocate(Dcenters(3,MaxnAtomsD,nTypesD))
Call ReadCenterInfo2(luoutput,ntypesD,nBatchesD,nAtomsOfTypeD,AngmomOfTypeD,&
     & nPrimOfTypeD,nContOfTypeD,MaxnAtomsD,MaxnPrimD,MaxnContD,&
     & startOrbitalOfTypeD,exponentsOfTypeD,ContractCoeffOfTypeD,&
     & Dcenters,spherical)

READ(LUOUTPUT,*) Outdim1,Outdim2,Outdim3,Outdim4

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

!     print*,'OutputDim1:',OutputDim1,OutputDim2,OutputDim3,OutputDim4
!     do d=1,OutputDim4
!        do c=1,OutputDim3
!           do b=1,OutputDim2
!              do a=1,OutputDim1            
!                print*,'int',integrals(a,b,c,d)
!             enddo
!          enddo
!       enddo
!    enddo
     
     allocate(ComparisonInt(Outdim1,Outdim2,Outdim3,Outdim4))
     READ(LUOUTPUT,*) ComparisonInt
     
     FAIL = .FALSE.
     do d=1,OutputDim4
        do c=1,OutputDim3
           do b=1,OutputDim2
              do a=1,OutputDim1   
                 IF(ABS(integrals(A,B,C,D)).GT.1.0E-10_8)THEN     
                    IF(ABS(integrals(a,b,c,d)-ComparisonInt(a,b,c,d)).GT.1.0E-10_8/ABS(ComparisonInt(a,b,c,d)))THEN
                       print*,'ERROR a,b,c,d = ',a,b,c,d
                       print*,'integrals(a,b,c,d)    ',integrals(a,b,c,d)
                       print*,'ComparisonInt(a,b,c,d)',ComparisonInt(a,b,c,d)
                       print*,'DIFF ',integrals(a,b,c,d)-ComparisonInt(a,b,c,d)
                       FAIL = .TRUE.
                    ENDIF
                 ELSE                    
                    IF(ABS(integrals(a,b,c,d)-ComparisonInt(a,b,c,d)).GT.1.0E-10_8)THEN
                       print*,'ERROR a,b,c,d = ',a,b,c,d
                       print*,'integrals(a,b,c,d)    ',integrals(a,b,c,d)
                       print*,'ComparisonInt(a,b,c,d)',ComparisonInt(a,b,c,d)
                       print*,'DIFF ',integrals(a,b,c,d)-ComparisonInt(a,b,c,d)
                       FAIL = .TRUE.
                    ENDIF
                 ENDIF
              enddo
           enddo
        enddo
     enddo
     IF(FAIL)THEN
        print*,'TEST FAILED '
     ELSE
        print*,'TEST SUCCEEDED '
     ENDIF
     deallocate(integrals)
     deallocate(ComparisonInt)

close(unit = LUOUTPUT)

deallocate(nAtomsOfTypeA)
deallocate(AngmomOfTypeA)
deallocate(nPrimOfTypeA)
deallocate(nContOfTypeA)
deallocate(startOrbitalOfTypeA)
deallocate(exponentsOfTypeA)
deallocate(ContractCoeffOfTypeA)
deallocate(Acenters)

deallocate(nAtomsOfTypeB)
deallocate(AngmomOfTypeB)
deallocate(nPrimOfTypeB)
deallocate(nContOfTypeB)
deallocate(startOrbitalOfTypeB)
deallocate(exponentsOfTypeB)
deallocate(ContractCoeffOfTypeB)
deallocate(Bcenters)

deallocate(nAtomsOfTypeC)
deallocate(AngmomOfTypeC)
deallocate(nPrimOfTypeC)
deallocate(nContOfTypeC)
deallocate(startOrbitalOfTypeC)
deallocate(exponentsOfTypeC)
deallocate(ContractCoeffOfTypeC)
deallocate(Ccenters)

deallocate(nAtomsOfTypeD)
deallocate(AngmomOfTypeD)
deallocate(nPrimOfTypeD)
deallocate(nContOfTypeD)
deallocate(startOrbitalOfTypeD)
deallocate(exponentsOfTypeD)
deallocate(ContractCoeffOfTypeD)
deallocate(Dcenters)

END PROGRAM IchorErimoduleTEST

