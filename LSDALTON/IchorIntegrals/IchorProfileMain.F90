PROGRAM IchorErimoduleTEST
use IchorPrecisionMod
#ifdef VAR_OPENACC
#ifdef VAR_PGF90
use openacc, only: acc_get_device_type,ACC_DEVICE_NONE,&
     & ACC_DEVICE_DEFAULT,ACC_DEVICE_HOST,ACC_DEVICE_NOT_HOST,&
     & acc_get_num_devices,acc_set_device_num, acc_init, acc_shutdown,&
     & ACC_DEVICE_NVIDIA
!ACC_DEVICE_NVIDIA only defined with PGI
#else
use openacc, only: acc_get_device_type,ACC_DEVICE_NONE,&
     & ACC_DEVICE_DEFAULT,ACC_DEVICE_HOST,ACC_DEVICE_NOT_HOST,&
     & acc_get_num_devices,acc_set_device_num, acc_init, acc_shutdown
#endif
#endif
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
Integer :: GabIdentifier, IchorGabID1, IchorGabID2,IchorOperatorSpec
logical :: SameLHSaos
real(8) :: THRESHOLD_CS,THRESHOLD_QQR,THRESHOLD_OD
real(8),pointer :: InputStorage(:)
real(8),pointer :: BATCHGAB(:),BATCHGCD(:)
Integer(kind=8) :: MaxMem,MaxMemAllocated,MemAllocated
Integer :: nBatchesA,nBatchesB,nBatchesC,nBatchesD,a,b,c,d
logical :: spherical
integer :: nbatchAstart2,nbatchAend2,nbatchBstart2,nbatchBend2
integer :: nbatchCstart2,nbatchCend2,nbatchDstart2,nbatchDend2,LUOUTPUT
logical :: SameRHSaos,SameODs,CRIT1,CRIT2,CRIT3,CRIT4,doLink,rhsDmat,CRIT5
integer :: nBasisA,nBasisB,nBasisC,nBasisD,iPassStart,iPassEnd
integer :: Ipass,iBasis1,iBasis2,iBasis3,iBasis4,ibasiselm(4)
real(8),pointer :: integrals(:,:,:,:),ComparisonInt(:,:,:,:)
integer :: iBasis1Q,iBasis2Q,iBasis3Q,iBasis4Q
integer,pointer :: iBasisA(:),iBasisB(:),iBasisC(:),iBasisD(:)
CHARACTER(len=20)    :: BASISTYPE(13)
integer :: iBASISTYPE(13),DebugIchorOption,ifilename,nKK,KK
integer :: DebugIchorOption2,nRepetitions,L,K,I,J,selected_device_number
character(len=100) :: filename
logical      :: SpecialPass,FAIL(13,13,13,13),ALLPASS
real(8) :: WALLTIMEFULL,WALLTIMECASE,WALLTIMEseg,WALLTIMEsegP,WALLTIMEsegQ
real(8) :: WALLTIMEseg1Prim,WALLTIMEGen,TIME1,TIME2,DELTAWALL,CPUTIME,WALLTIME,GPUMAXMEM
integer(kind=accdevkind) :: acc_device_type

#ifdef VAR_OPENACC
acc_device_type = acc_get_device_type()
print*,'acc_get_device_type = ',acc_device_type
print*,'================================'
!4 device types are always supported 
print*,'ACC_DEVICE_NONE     = ',ACC_DEVICE_NONE
print*,'ACC_DEVICE_DEFAULT  = ',ACC_DEVICE_DEFAULT
print*,'ACC_DEVICE_HOST     = ',ACC_DEVICE_HOST
print*,'ACC_DEVICE_NOT_HOST = ',ACC_DEVICE_NOT_HOST
#ifdef VAR_PGF90
!PGF90 supports 
print*,'ACC_DEVICE_NVIDIA   = ',ACC_DEVICE_NVIDIA
#endif
print*,'================================'
#ifdef VAR_PGF90
IF(acc_device_type.EQ.ACC_DEVICE_NVIDIA)print*,'ACC_DEVICE_NVIDIA have been selected'
#endif
IF(acc_device_type.EQ.ACC_DEVICE_NONE)THEN
   print*,'ACC_DEVICE_NONE have been selected'
   print*,'please use the command'
   print*,'export ACC_DEVICE=NVIDIA'
   print*,'to chose the NVIDIA graphics card or '
   print*,'export ACC_DEVICE=HOST'
   print*,'to chose to run the GPU kernel on the CPU'
ENDIF
#ifdef VAR_PGF90
print*,'There are ',acc_get_num_devices(acc_device_type),'devices'
print*,'The Program only support 1 device at present'
selected_device_number = 2
IF(acc_get_num_devices(acc_device_type).GT.1)THEN
   IF(selected_device_number .EQ. 0)THEN
      print*,'Using default behavior of which device to use'
      print*,'Change this behavior by choosing device number in input or'
      print*,'export ACC_DEVICE_NUM=1'
   ELSE
      print*,'The selected device number=',selected_device_number
      call acc_set_device_num(2,0)
   ENDIF
ENDIF
#endif
call acc_init(acc_device_type)
#endif


SpecialPass = .FALSE.
DebugIchorOption = 9
LUPRI = 6
IPRINT = 0
spherical = .TRUE.
nKK = 5
IF(DebugIchorOption.EQ.9)nKK = 1

WALLTIMEFULL=0.0E0_8

DO KK=1,nKK      
WALLTIMECASE=0.0E0_8
DebugIchorOption2 = kk
IF(DebugIchorOption.EQ.9)DebugIchorOption2 = DebugIchorOption
SELECT CASE(DebugIchorOption2)
CASE(1)
   !Special types
   !   A          B          C          D          OVERALL
   !Seg1Prim     Seg        Gen      Seg1Prim      SegP
   nbasisA = 3; nbasisB = 3; nbasisC = 3; nbasisD = 3
   SpecialPass = .TRUE.
   allocate(iBasisA(nBasisA))
   allocate(iBasisB(nBasisB))
   allocate(iBasisC(nBasisC))
   allocate(iBasisD(nBasisD))
   iBasisA(1) = 1; iBasisA(2) = 2; iBasisA(3) = 3  !Seg1Prim
   iBasisB(1) = 4; iBasisB(2) = 5; iBasisB(3) = 6  !Seg
   iBasisC(1) = 7; iBasisC(2) = 8; iBasisC(3) = 9  !Gen
   iBasisD(1) = 1; iBasisD(2) = 2; iBasisD(3) = 3  !Seg1Prim
CASE(2)
   !Special types
   !   A          B          C           D          OVERALL
   !  Seg        Seg1Prim   Seg1Prim  Seg1Prim      Seg
   nbasisA = 3; nbasisB = 3; nbasisC = 3; nbasisD = 3
   allocate(iBasisA(nBasisA))
   allocate(iBasisB(nBasisB))
   allocate(iBasisC(nBasisC))
   allocate(iBasisD(nBasisD))
   !Pure Seg (S,P,D ; S,P,D | S,P,D ; S,P,D)
   iBasisA(1) = 4; iBasisA(2) = 5; iBasisA(3) = 6
   iBasisB(1) = 1; iBasisB(2) = 2; iBasisB(3) = 3
   iBasisC(1) = 1; iBasisC(2) = 2; iBasisC(3) = 3
   iBasisD(1) = 1; iBasisD(2) = 2; iBasisD(3) = 3
CASE(3)
   !Special types
   !   A          B          C          D          OVERALL
   !Seg1Prim   Seg1Prim   Seg1Prim   Seg1Prim      Seg1Prim 
   nbasisA = 3; nbasisB = 3; nbasisC = 3; nbasisD = 3
   allocate(iBasisA(nBasisA))
   allocate(iBasisB(nBasisB))
   allocate(iBasisC(nBasisC))
   allocate(iBasisD(nBasisD))
   !Pure Seg1Prim (S,P,D ; S,P,D | S,P,D ; S,P,D)
   iBasisA(1) = 1; iBasisA(2) = 2; iBasisA(3) = 3
   iBasisB(1) = 1; iBasisB(2) = 2; iBasisB(3) = 3
   iBasisC(1) = 1; iBasisC(2) = 2; iBasisC(3) = 3
   iBasisD(1) = 1; iBasisD(2) = 2; iBasisD(3) = 3
CASE(4)
   !Special types
   !   A          B          C          D          OVERALL
   !Seg1Prim     Gen        Seg      Seg1Prim      SegQ
   nbasisA = 3; nbasisB = 3; nbasisC = 3; nbasisD = 3
   SpecialPass = .TRUE.
   allocate(iBasisA(nBasisA))
   allocate(iBasisB(nBasisB))
   allocate(iBasisC(nBasisC))
   allocate(iBasisD(nBasisD))
   iBasisA(1) = 1; iBasisA(2) = 2; iBasisA(3) = 3  !Seg1Prim
   iBasisB(1) = 7; iBasisB(2) = 8; iBasisB(3) = 9  !Gen
   iBasisC(1) = 4; iBasisC(2) = 5; iBasisC(3) = 6  !Seg
   iBasisD(1) = 1; iBasisD(2) = 2; iBasisD(3) = 3  !Seg1Prim
CASE(5)
   !Gen1
   !Special types
   !   A          B          C          D          OVERALL
   !  Gen        Gen        Gen        Gen         Gen    
   !Pure Gen (S,P,D ; S,P,D | S,P,D ; S,P,D)
   nbasisA = 3; nbasisB = 3; nbasisC = 3; nbasisD = 3
   SpecialPass = .TRUE. !otherwise too expensive
   allocate(iBasisA(nBasisA))
   allocate(iBasisB(nBasisB))
   allocate(iBasisC(nBasisC))
   allocate(iBasisD(nBasisD))
   iBasisA(1) = 7; iBasisA(2) = 8; iBasisA(3) = 9
   iBasisB(1) = 7; iBasisB(2) = 8; iBasisB(3) = 9
   iBasisC(1) = 7; iBasisC(2) = 8; iBasisC(3) = 9
   iBasisD(1) = 7; iBasisD(2) = 8; iBasisD(3) = 9
CASE DEFAULT
   FileName = 'IchorProfTestUnitTest_genDUnitTest_genDUnitTest_genDUnitTest_genD1'
   nBasisA = 1
   nBasisB = 1
   nBasisC = 1
   nBasisD = 1
   nbasisA = 1; nbasisB = 1; nbasisC = 1; nbasisD = 1
   allocate(iBasisA(nBasisA))
   allocate(iBasisB(nBasisB))
   allocate(iBasisC(nBasisC))
   allocate(iBasisD(nBasisD))
   do iBasis1Q=1,nbasisA
      iBasisA(iBasis1Q) = iBasis1Q
   enddo
   do iBasis2Q=1,nbasisB
      iBasisB(iBasis2Q) = iBasis2Q
   enddo
   do iBasis3Q=1,nbasisC
      iBasisC(iBasis3Q) = iBasis3Q
   enddo
   do iBasis4Q=1,nbasisD
      iBasisD(iBasis4Q) = iBasis4Q
   enddo
END SELECT
IpassStart = 1; IpassEnd = 1
BASISTYPE(1) = 'UnitTest_segS1p     '; iBASISTYPE(1) = 15
BASISTYPE(2) = 'UnitTest_segP1p     '; iBASISTYPE(2) = 15
BASISTYPE(3) = 'UnitTest_segD1p     '; iBASISTYPE(3) = 15
BASISTYPE(4) = 'UnitTest_segS       '; iBASISTYPE(4) = 13
BASISTYPE(5) = 'UnitTest_segP       '; iBASISTYPE(5) = 13
BASISTYPE(6) = 'UnitTest_segD       '; iBASISTYPE(6) = 13
BASISTYPE(7) = 'UnitTest_genS       '; iBASISTYPE(7) = 13
BASISTYPE(8) = 'UnitTest_genP       '; iBASISTYPE(8) = 13
BASISTYPE(9) = 'UnitTest_genD       '; iBASISTYPE(9) = 13
BASISTYPE(10) = 'UnitTest_segSP      '; iBASISTYPE(10) = 14
FAIL = .FALSE.
do Ipass = IpassStart,IpassEnd
 do iBasis1Q = 1,nBasisA
  iBasis1 = iBasisA(iBasis1Q)
  do iBasis2Q = 1,nBasisB
   iBasis2 = iBasisB(iBasis2Q)
   do iBasis3Q = 1,nBasisC
    iBasis3 = iBasisC(iBasis3Q)
    do iBasis4Q = 1,nBasisD
     iBasis4 = iBasisD(iBasis4Q)

     ibasiselm(1) = iBasis1
     ibasiselm(2) = iBasis2
     ibasiselm(3) = iBasis3
     ibasiselm(4) = iBasis4

     IF(DebugIchorOption.NE.9)THEN
        print*,'DebugIchorOption'
        do A=1,100
           filename(A:A) = ' '
        enddo
        filename(1:13) = 'IchorProfTest'
        ifilename = 14
        do A = 1,4       
           filename(ifilename:ifilename+iBASISTYPE(iBasiselm(A))-1) =  BASISTYPE(iBasiselm(A))(1:iBASISTYPE(iBasiselm(A)))
           ifilename = ifilename+iBASISTYPE(iBasiselm(A)) 
        enddo
        WRITE(filename(ifilename:ifilename),'(I1)') Ipass
        ifilename = ifilename + 1
        LUOUTPUT = 12
        !        print*,'FileName',FileName
     ELSE
        print*,'DebugIchorOption = 9 ',Filename
     ENDIF
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
     GPUMAXMEM = 2.0E0_8
     call SetIchorGPUMaxMem(GPUMAXMEM)
     
     SameLHSaos = .FALSE. 
     SameRHSaos = .FALSE. 
     SameODs = .FALSE. 
     
     call GetIchorPermuteParameter(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
     call GetIchorFileStorageIdentifier(filestorageIdentifier)
     call GetIchorOpereratorIntSpec('C',IchorOperatorSpec) !Coulomb Operator
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

!     IF(DebugIchorOption.EQ.9)THEN
!        nRepetitions = 5
!     ELSEIF(DebugIchorOption2.EQ.3)THEN
!        nRepetitions = 3
!     ELSE
        nRepetitions = 1
!     ENDIF

     allocate(ComparisonInt(Outdim1,Outdim2,Outdim3,Outdim4))

     IF(MIN(Outputdim1,Outputdim2,Outputdim3,Outputdim4).LT.7)THEN
        DO L=1,7
           DO K=1,7
              DO J=1,7
                 DO I=1,7
                    READ(LUOUTPUT,*) ComparisonInt(I,J,K,L)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDIF
     DO L=1,Outputdim4,10
        DO K=1,Outputdim3,10
           DO J=1,Outputdim2,10
              DO I=1,Outputdim1,10
                 READ(LUOUTPUT,*) ComparisonInt(I,J,K,L)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!     READ(LUOUTPUT,*) ComparisonInt
     DO I=1,nRepetitions
        !print*,'THRESHOLD_CS',THRESHOLD_CS,'THRESHOLD_OD',THRESHOLD_OD
        !=====================================================================
        !  Main Call
        !=====================================================================
        
        CALL ICHOR_GETTIM(CPUTIME,WALLTIME)

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
             & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorOperatorSpec,&
             & IchorInputDim1,IchorInputDim2,IchorInputDim3,&
             & InputStorage,IchorParSpec,IchorScreenSpec,THRESHOLD_OD,THRESHOLD_CS,&
             & THRESHOLD_QQR,IchorGabID1,IchorGabID2,IchorDebugSpec,&
             & IchorAlgoSpec,IchorPermuteSpec,filestorageIdentifier,MaxMem,&
             & MaxFileStorage,MaxMemAllocated,MemAllocated,&
             & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
             & integrals,lupri)


        CALL ICHOR_GETTIM(TIME1,TIME2)
        DELTAWALL=TIME2-WALLTIME
        write(lupri,'(A,A,A,A,A,A,A,A,A,F16.8)')'BASIS(',TRIM(BASISTYPE(iBasis1)),',',TRIM(BASISTYPE(iBasis2)),&
             & ',',TRIM(BASISTYPE(iBasis3)),',',TRIM(BASISTYPE(iBasis4)),') Wall Time=',DELTAWALL
        WALLTIMECASE = WALLTIMECASE + DELTAWALL
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
        IF(I.EQ.1.OR.MOD(I,5).EQ.0)THEN
           FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .FALSE.

           IF(MIN(Outputdim1,Outputdim2,Outputdim3,Outputdim4).LT.7)THEN
              do d=1,7
                 do c=1,7
                    do b=1,7
                       do a=1,7
                         IF(ABS(integrals(A,B,C,D)).GT.1.0E-10_8)THEN     
                          IF(ABS(integrals(a,b,c,d)-ComparisonInt(a,b,c,d)).GT.1.0E-10_8/ABS(ComparisonInt(a,b,c,d)))THEN
                             FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .TRUE.
                             print*,'ERROR a,b,c,d = ',a,b,c,d
                             print*,'integrals(a,b,c,d)    ',integrals(a,b,c,d)
                             print*,'ComparisonInt(a,b,c,d)',ComparisonInt(a,b,c,d)
                             print*,'DIFF ',integrals(a,b,c,d)-ComparisonInt(a,b,c,d)
                             FAIL = .TRUE.
                          ENDIF
                         ELSE                    
                          IF(ABS(integrals(a,b,c,d)-ComparisonInt(a,b,c,d)).GT.1.0E-10_8)THEN
                             FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .TRUE.
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
           ENDIF
           do d=1,OutputDim4,10
              do c=1,OutputDim3,10
                 do b=1,OutputDim2,10
                    do a=1,OutputDim1,10   
                       IF(ABS(integrals(A,B,C,D)).GT.1.0E-10_8)THEN     
                          IF(ABS(integrals(a,b,c,d)-ComparisonInt(a,b,c,d)).GT.1.0E-10_8/ABS(ComparisonInt(a,b,c,d)))THEN
                             FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .TRUE.
                             print*,'ERROR a,b,c,d = ',a,b,c,d
                             print*,'integrals(a,b,c,d)    ',integrals(a,b,c,d)
                             print*,'ComparisonInt(a,b,c,d)',ComparisonInt(a,b,c,d)
                             print*,'DIFF ',integrals(a,b,c,d)-ComparisonInt(a,b,c,d)
                             FAIL = .TRUE.
                          ENDIF
                       ELSE                    
                          IF(ABS(integrals(a,b,c,d)-ComparisonInt(a,b,c,d)).GT.1.0E-10_8)THEN
                             FAIL(iBasis1,ibasis2,ibasis3,ibasis4) = .TRUE.
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
           IF(FAIL(iBasis1,ibasis2,ibasis3,ibasis4))THEN
              write(*,'(A)')'TEST FAILED '
!           ELSE
!              write(*,'(A)')'TEST SUCCEEDED '
           ENDIF
        ENDIF
     enddo
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
    enddo
   enddo
  enddo
 enddo
 ALLPASS = .TRUE.
 write(*,'(A)')'Summary'
 do iBasis1Q = 1,nBasisA
  iBasis1 = iBasisA(iBasis1Q)
  do iBasis2Q = 1,nBasisB
   iBasis2 = iBasisB(iBasis2Q)
   do iBasis3Q = 1,nBasisC
    iBasis3 = iBasisC(iBasis3Q)
    do iBasis4Q = 1,nBasisD
     iBasis4 = iBasisD(iBasis4Q)

     IF(FAIL(iBasis1,ibasis2,ibasis3,ibasis4)) THEN
        write(*,'(A,A,A,A,A,A,A,A,A,I1,A)')'BASIS(',BASISTYPE(iBasis1)(10:15),',',&
             & BASISTYPE(iBasis2)(10:15),',',BASISTYPE(iBasis3)(10:15),',',&
             & BASISTYPE(iBasis4)(10:15),',',ipass,') FAILED'
        ALLPASS = .FALSE.
!     ELSE
!        write(*,'(A,A,A,A,A,A,A,A,A,I1,A)')'BASIS(',BASISTYPE(iBasis1)(10:15),',',&
!             & BASISTYPE(iBasis2)(10:15),',',BASISTYPE(iBasis3)(10:15),',',&
!             & BASISTYPE(iBasis4)(10:15),',',ipass,') SUCCESSFUL'
     ENDIF     
    enddo
   enddo
  enddo
 enddo
enddo

IF(ALLPASS)THEN
   WRITE(*,'(A)')'Ichor Integral UnitTest: SUCCESSFUL'
ELSE
   WRITE(*,'(A)')'Ichor Integrals UnitTest: FAILED'
ENDIF

WALLTIMEFULL=WALLTIMEFULL+WALLTIMECASE

SELECT CASE(DebugIchorOption2)
CASE(1)
   WALLTIMEsegP = WALLTIMECASE
   WRITE(*,'(A,F16.8)')'SegP Wall Time     =',WALLTIMEsegP
CASE(2)
   WALLTIMEseg = WALLTIMECASE
   WRITE(*,'(A,F16.8)')'Seg Wall Time      =',WALLTIMEseg
CASE(3)
   WALLTIMEseg1Prim = WALLTIMECASE
   WRITE(*,'(A,F16.8)')'Seg1Prim Wall Time =',WALLTIMEseg1Prim
CASE(4)
   WALLTIMEsegQ = WALLTIMECASE
   WRITE(*,'(A,F16.8)')'SegQ Wall Time     =',WALLTIMEsegQ
CASE(5)
   WALLTIMEGen = WALLTIMECASE
   WRITE(*,'(A,F16.8)')'Gen Wall Time      =',WALLTIMEGen
END SELECT

ENDDO

WRITE(*,'(A,F16.8)')'SegP Wall Time     =',WALLTIMEsegP
WRITE(*,'(A,F16.8)')'Seg Wall Time      =',WALLTIMEseg
WRITE(*,'(A,F16.8)')'Seg1Prim Wall Time =',WALLTIMEseg1Prim
WRITE(*,'(A,F16.8)')'SegQ Wall Time     =',WALLTIMEsegQ
WRITE(*,'(A,F16.8)')'Gen Wall Time      =',WALLTIMEGen
WRITE(*,'(A,F16.8)')'Total Wall Time    =',WALLTIMEFULL

call acc_shutdown(acc_device_type)

CONTAINS
subroutine GetIchorOpereratorIntSpec(intSpec,IchorOperatorSpec)
  implicit none
  character :: intspec 
  integer,intent(inout) :: IchorOperatorSpec
  IF (intSpec.EQ.'C') THEN
     ! Regular Coulomb operator 1/r12                                                                   
     call GetIchorOpererator('Coulomb',IchorOperatorSpec)
  ELSE IF (intSpec.EQ.'G') THEN
     ! The Gaussian geminal operator g                                                                  
     call GetIchorOpererator('GGem   ',IchorOperatorSpec)
  ELSE IF (intSpec.EQ.'F') THEN
     ! The Gaussian geminal divided by the Coulomb operator g/r12                                       
     call GetIchorOpererator('GGemCou',IchorOperatorSpec)
  ELSE IF (intSpec.EQ.'D') THEN
     ! The double commutator [[T,g],g]                                                                  
     call GetIchorOpererator('GGemGrd',IchorOperatorSpec)
  ELSE IF (intSpec.EQ.'2') THEN
     ! The Gaussian geminal operator squared g^2                                                        
     call GetIchorOpererator('GGemSq ',IchorOperatorSpec)
  ELSE
     STOP 'Error in specification of operator in GetIchorOpereratorIntSpec'
  ENDIF
end subroutine GetIchorOpereratorIntSpec

END PROGRAM IchorErimoduleTEST

