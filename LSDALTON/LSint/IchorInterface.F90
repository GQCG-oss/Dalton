!> @file
!> Contains the main Ichor integral drivers
!> Ichor is the "Integral Code Hand Optimized for Rapid evaluation" 
!> Ichor is the ethereal golden fluid that is the blood of the greek gods

!> \brief Main Ichor drivers for the calculation of integrals 
!> based on the McMurchie-Davidson(McM)/Obara Saika(OS)/Head-Gordon-Pople(HGP)/Rys 
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorErimoduleHost
  use precision
  use TYPEDEFTYPE,only: lssetting,BASISSETINFO,MOLECULEINFO
  use memory_handling, only: mem_alloc,mem_dealloc
  use IchorErimodule !IchorEri plus the parameters
public:: MAIN_ICHORERI_DRIVER
private
CONTAINS
SUBROUTINE MAIN_ICHORERI_DRIVER(LUPRI,IPRINT,setting,dim1,dim2,dim3,dim4,integrals,Spherical)
implicit none
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: LUPRI,IPRINT,dim1,dim2,dim3,dim4
real(realk),intent(inout) :: integrals(dim1,dim2,dim3,dim4)
logical,intent(in)        :: spherical
!
integer                :: nTypes
!A
integer :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA
integer,pointer     :: nAtomsOfTypeA(:),nPrimOfTypeA(:)
integer,pointer     :: nContOfTypeA(:),AngmomOfTypeA(:),startOrbitalOfTypeA(:,:)
real(realk),pointer :: exponentsOfTypeA(:,:),Acenters(:,:,:),ContractCoeffOfTypeA(:,:,:)
!B
integer :: nTypesB,MaxnAtomsB,MaxnPrimB,MaxnContB
integer,pointer     :: nAtomsOfTypeB(:),nPrimOfTypeB(:)
integer,pointer     :: nContOfTypeB(:),AngmomOfTypeB(:),startOrbitalOfTypeB(:,:)
real(realk),pointer :: exponentsOfTypeB(:,:),Bcenters(:,:,:),ContractCoeffOfTypeB(:,:,:)
!C
integer :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC
integer,pointer     :: nAtomsOfTypeC(:),nPrimOfTypeC(:)
integer,pointer     :: nContOfTypeC(:),AngmomOfTypeC(:),startOrbitalOfTypeC(:,:)
real(realk),pointer :: exponentsOfTypeC(:,:),Ccenters(:,:,:),ContractCoeffOfTypeC(:,:,:)
!D
integer :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD
integer,pointer     :: nAtomsOfTypeD(:),nPrimOfTypeD(:)
integer,pointer     :: nContOfTypeD(:),AngmomOfTypeD(:),startOrbitalOfTypeD(:,:)
real(realk),pointer :: exponentsOfTypeD(:,:),Dcenters(:,:,:),ContractCoeffOfTypeD(:,:,:)
!job specification
integer :: SphericalSpec,IchorJobSpec,IchorInputSpec,IchorParSpec,IchorScreenSpec
Integer :: IchorInputDim1,IchorInputDim2,IchorInputDim3,IchorDebugSpec,IchorAlgoSpec
integer :: filestorageIdentifier,MaxFileStorage,IchorPermuteSpec
Integer :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5
real(realk),pointer :: InputStorage(:)
Integer(kind=long) :: MaxMem,MaxMemAllocated,MemAllocated

!A
call Determine_nTypesForBatch(setting%BASIS(1)%p%REGULAR,ntypesA)
call mem_alloc(nAtomsOfTypeA,ntypesA)
call mem_alloc(AngmomOfTypeA,ntypesA)
call mem_alloc(nPrimOfTypeA,ntypesA)
call mem_alloc(nContOfTypeA,ntypesA)

call build_TypeInfo1(setting%MOLECULE(1)%p,setting%BASIS(1)%p%REGULAR,nTypesA,&
     & nAtomsOfTypeA,AngmomOfTypeA,nPrimOfTypeA,nContOfTypeA,MaxnAtomsA,MaxnPrimA,MaxnContA)

call mem_alloc(startOrbitalOfTypeA,MaxNatomsA,ntypesA)
call mem_alloc(exponentsOfTypeA,MaxnPrimA,ntypesA)
call mem_alloc(ContractCoeffOfTypeA,MaxnPrimA,MaxnContA,ntypesA)
call mem_alloc(Acenters,3,MaxnAtomsA,nTypesA)

call build_TypeInfo2(setting%MOLECULE(1)%p,setting%BASIS(1)%p%REGULAR,nTypesA,&
     & spherical,MaxnAtomsA,MaxnPrimA,MaxnContA,startOrbitalOfTypeA,&
     & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters)

!B
call Determine_nTypesForBatch(setting%BASIS(2)%p%REGULAR,ntypesB)
call mem_alloc(nAtomsOfTypeB,ntypesB)
call mem_alloc(AngmomOfTypeB,ntypesB)
call mem_alloc(nPrimOfTypeB,ntypesB)
call mem_alloc(nContOfTypeB,ntypesB)

call build_TypeInfo1(setting%MOLECULE(2)%p,setting%BASIS(2)%p%REGULAR,nTypesB,&
     & nAtomsOfTypeB,AngmomOfTypeB,nPrimOfTypeB,nContOfTypeB,MaxnAtomsB,MaxnPrimB,MaxnContB)

call mem_alloc(startOrbitalOfTypeB,MaxNatomsB,ntypesB)
call mem_alloc(exponentsOfTypeB,MaxnPrimB,ntypesB)
call mem_alloc(ContractCoeffOfTypeB,MaxnPrimB,MaxnContB,ntypesB)
call mem_alloc(Bcenters,3,MaxnAtomsB,nTypesB)

call build_TypeInfo2(setting%MOLECULE(2)%p,setting%BASIS(2)%p%REGULAR,nTypesB,&
     & spherical,MaxnAtomsB,MaxnPrimB,MaxnContB,startOrbitalOfTypeB,&
     & exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters)

!C
call Determine_nTypesForBatch(setting%BASIS(3)%p%REGULAR,ntypesC)
call mem_alloc(nAtomsOfTypeC,ntypesC)
call mem_alloc(AngmomOfTypeC,ntypesC)
call mem_alloc(nPrimOfTypeC,ntypesC)
call mem_alloc(nContOfTypeC,ntypesC)

call build_TypeInfo1(setting%MOLECULE(3)%p,setting%BASIS(3)%p%REGULAR,nTypesC,&
     & nAtomsOfTypeC,AngmomOfTypeC,nPrimOfTypeC,nContOfTypeC,MaxnAtomsC,MaxnPrimC,MaxnContC)

call mem_alloc(startOrbitalOfTypeC,MaxNatomsC,ntypesC)
call mem_alloc(exponentsOfTypeC,MaxnPrimC,ntypesC)
call mem_alloc(ContractCoeffOfTypeC,MaxnPrimC,MaxnContC,ntypesC)
call mem_alloc(Ccenters,3,MaxnAtomsC,nTypesC)

call build_TypeInfo2(setting%MOLECULE(3)%p,setting%BASIS(3)%p%REGULAR,nTypesC,&
     & spherical,MaxnAtomsC,MaxnPrimC,MaxnContC,startOrbitalOfTypeC,&
     & exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters)

!D
call Determine_nTypesForBatch(setting%BASIS(4)%p%REGULAR,ntypesD)
call mem_alloc(nAtomsOfTypeD,ntypesD)
call mem_alloc(AngmomOfTypeD,ntypesD)
call mem_alloc(nPrimOfTypeD,ntypesD)
call mem_alloc(nContOfTypeD,ntypesD)
call build_TypeInfo1(setting%MOLECULE(4)%p,setting%BASIS(4)%p%REGULAR,nTypesD,&
     & nAtomsOfTypeD,AngmomOfTypeD,nPrimOfTypeD,nContOfTypeD,MaxnAtomsD,MaxnPrimD,MaxnContD)

call mem_alloc(startOrbitalOfTypeD,MaxNatomsD,ntypesD)
call mem_alloc(exponentsOfTypeD,MaxnPrimD,ntypesD)
call mem_alloc(ContractCoeffOfTypeD,MaxnPrimD,MaxnContD,ntypesD)
call mem_alloc(Dcenters,3,MaxnAtomsD,nTypesD)

call build_TypeInfo2(setting%MOLECULE(4)%p,setting%BASIS(4)%p%REGULAR,nTypesD,&
     & spherical,MaxnAtomsD,MaxnPrimD,MaxnContD,startOrbitalOfTypeD,&
     & exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters)

SphericalSpec=SphericalParam
IchorJobSpec=IcorJobEri
IchorInputSpec=1
IchorInputDim1=1
IchorInputDim2=1
IchorInputDim3=1
call mem_alloc(InputStorage,1)
IchorParSpec=IchorParNone
IchorScreenSpec=IchorScreenNone
IchorDebugSpec=iprint!IchorDebugNone
IchorAlgoSpec=IchorAlgoOS
IchorPermuteSpec=IchorPermuteFFF
filestorageIdentifier=1
MaxMem = 0 
MaxFileStorage = 0
MaxMemAllocated = 0 
MemAllocated = 0
OutputDim1=dim1
OutputDim2=dim2
OutputDim3=dim3
OutputDim4=dim4
OutputDim5=1
!=====================================================================
!  Main Call
!=====================================================================
 
call IchorEri(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
     & AngmomOfTypeA,nAtomsOfTypeA,nPrimOfTypeA,nContOfTypeA,&
     & startOrbitalOfTypeA,Acenters,exponentsOfTypeA,ContractCoeffOfTypeA,&
     & nTypesB,MaxNatomsB,MaxnPrimB,MaxnContB,&
     & AngmomOfTypeB,nAtomsOfTypeB,nPrimOfTypeB,nContOfTypeB,&
     & startOrbitalOfTypeB,Bcenters,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & nTypesC,MaxNatomsC,MaxnPrimC,MaxnContC,&
     & AngmomOfTypeC,nAtomsOfTypeC,nPrimOfTypeC,nContOfTypeC,&
     & startOrbitalOfTypeC,Ccenters,exponentsOfTypeC,ContractCoeffOfTypeC,&
     & nTypesD,MaxNatomsD,MaxnPrimD,MaxnContD,&
     & AngmomOfTypeD,nAtomsOfTypeD,nPrimOfTypeD,nContOfTypeD,&
     & startOrbitalOfTypeD,Dcenters,exponentsOfTypeD,ContractCoeffOfTypeD,&
     & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorInputDim1,IchorInputDim2,&
     & IchorInputDim3,&
     & InputStorage,IchorParSpec,IchorScreenSpec,IchorDebugSpec,&
     & IchorAlgoSpec,IchorPermuteSpec,filestorageIdentifier,MaxMem,&
     & MaxFileStorage,MaxMemAllocated,MemAllocated,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
     & integrals,lupri)
!=====================================================================

call mem_dealloc(InputStorage)

!=====================================================================

!free space
!A
call mem_dealloc(nAtomsOfTypeA)
call mem_dealloc(AngmomOfTypeA)
call mem_dealloc(nPrimOfTypeA)
call mem_dealloc(nContOfTypeA)
call mem_dealloc(startOrbitalOfTypeA)
call mem_dealloc(exponentsOfTypeA)
call mem_dealloc(ContractCoeffOfTypeA)
call mem_dealloc(Acenters)

!B
call mem_dealloc(nAtomsOfTypeB)
call mem_dealloc(AngmomOfTypeB)
call mem_dealloc(nPrimOfTypeB)
call mem_dealloc(nContOfTypeB)
call mem_dealloc(startOrbitalOfTypeB)
call mem_dealloc(exponentsOfTypeB)
call mem_dealloc(ContractCoeffOfTypeB)
call mem_dealloc(Bcenters)
!C
call mem_dealloc(nAtomsOfTypeC)
call mem_dealloc(AngmomOfTypeC)
call mem_dealloc(nPrimOfTypeC)
call mem_dealloc(nContOfTypeC)
call mem_dealloc(startOrbitalOfTypeC)
call mem_dealloc(exponentsOfTypeC)
call mem_dealloc(ContractCoeffOfTypeC)
call mem_dealloc(Ccenters)

!D
call mem_dealloc(nAtomsOfTypeD)
call mem_dealloc(AngmomOfTypeD)
call mem_dealloc(nPrimOfTypeD)
call mem_dealloc(nContOfTypeD)
call mem_dealloc(startOrbitalOfTypeD)
call mem_dealloc(exponentsOfTypeD)
call mem_dealloc(ContractCoeffOfTypeD)
call mem_dealloc(Dcenters)

END SUBROUTINE MAIN_ICHORERI_DRIVER

Subroutine Determine_nTypesForBatch(BASISINFO,ntypes)
implicit none
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO),intent(in):: BASISINFO
!> the number of different types of batches
INTEGER,intent(inout)        :: ntypes
!
INTEGER                   :: SUM1,I,K
IF(BASISINFO%natomtypes.EQ. 0)&
     & CALL LSQUIT('Error Determine_nTypesForBatch called with empty basis',-1)
ntypes=0
DO I=1,BASISINFO%natomtypes
   SUM1=0
   DO K=1,BASISINFO%ATOMTYPE(I)%nAngmom
      SUM1=SUM1+BASISINFO%ATOMTYPE(I)%SHELL(K)%nsegments
   ENDDO
   ntypes=ntypes + SUM1
ENDDO
END Subroutine DETERMINE_NTYPESFORBATCH

Subroutine build_TypeInfo1(MOLECULE,BASISINFO,nTypes,nAtomsOfType,&
     & AngmomOfType,nPrimOfType,nContOfType,MaxnAtoms,MaxnPrim,MaxnCont)
implicit none
!> the number of different types of batches
INTEGER,intent(in)           :: ntypes
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO),intent(in) :: MOLECULE
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO),intent(in):: BASISINFO
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
integer,intent(inout)    :: MaxnCont
!
INTEGER,pointer          :: MODELTYPES(:),MODELBATCHTYPES(:,:,:)
INTEGER                  :: maxseg,maxang,I,K,L,iBatchType,R,icharge,type,iseg
INTEGER                  :: nBatchType,nrow,ncol,irow,icol,orbitalindex,nOrbComp
!assume nofamily
call mem_alloc(MODELTYPES,BASISINFO%natomtypes)

maxSeg = 0
maxAng = 0
DO I=1,BASISINFO%natomtypes
   maxAng = MAX(maxAng,BASISINFO%ATOMTYPE(I)%nAngmom)
   DO K=1,BASISINFO%ATOMTYPE(I)%nAngmom
      maxSeg = MAX(maxseg,BASISINFO%ATOMTYPE(I)%SHELL(K)%nsegments)
   ENDDO
ENDDO

call mem_alloc(MODELBATCHTYPES,maxSeg,maxAng,BASISINFO%natomtypes)
L=0
iBatchType =  0 
DO I=1,BASISINFO%natomtypes
   L=L+1
   MODELTYPES(I)=L
   DO K=1,BASISINFO%ATOMTYPE(I)%nAngmom
      DO iseg = 1,BASISINFO%ATOMTYPE(I)%SHELL(K)%nsegments
         iBatchType = iBatchType + 1 
         MODELBATCHTYPES(iseg,K,L) = iBatchType
      ENDDO
   ENDDO
ENDDO
nBatchType = iBatchType
IF(nBatchType.NE.ntypes)call lsquit('dim mismatch in Determine_nAtomsOfTypes',-1)
!init
do iBatchType=1,nBatchType
   nAtomsOfType(iBatchType) = 0
enddo
!build nAtomsOfType
R = BASISINFO%Labelindex
DO I=1,MOLECULE%natoms  
 IF(MOLECULE%ATOM(I)%pointcharge)CYCLE !no basis functions on this point charge
 IF(R.EQ. 0)THEN
    ICHARGE = INT(MOLECULE%ATOM(I)%CHARGE)
    type = BASISINFO%CHARGEINDEX(ICHARGE)
 ELSE
    type = MOLECULE%ATOM(I)%IDtype(R)
 ENDIF
 L=MODELTYPES(type)
 IF(BASISINFO%ATOMTYPE(type)%nAngmom.EQ.0) CYCLE
 DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
  DO iseg = 1,BASISINFO%ATOMTYPE(type)%SHELL(K)%nsegments
   iBatchType = MODELBATCHTYPES(iseg,K,L)
   nAtomsOfType(iBatchType) = nAtomsOfType(iBatchType) + 1
   AngmomOfType(iBatchType) = K-1
   nrow = BASISINFO%ATOMTYPE(type)%SHELL(K)%segment(iseg)%nrow
   ncol = BASISINFO%ATOMTYPE(type)%SHELL(K)%segment(iseg)%ncol
   nPrimOfType(iBatchType) = nrow
   nContOfType(iBatchType) = ncol
  ENDDO
 ENDDO
ENDDO
call mem_dealloc(MODELTYPES)
call mem_dealloc(MODELBATCHTYPES)

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

end Subroutine Build_TypeInfo1

Subroutine Build_TypeInfo2(MOLECULE,BASISINFO,nTypes,spherical,MaxnAtoms,MaxnPrim,&
     & MaxnCont,startOrbitalOfType,exponentsOfType,ContractCoeffOfType,CentersOfType)
implicit none
!> if it is spherical or cartesian GTOs
LOGICAL,intent(in)           :: spherical
!> the number of different types of batches
INTEGER,intent(in)           :: ntypes
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO),intent(in) :: MOLECULE
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO),intent(in):: BASISINFO
!> the maximum number of Atoms with a given shell type
integer,intent(in)    :: MaxnAtoms
!> the maximum number of Primitive in a given shell type
integer,intent(in)    :: MaxnPrim
!> the maximum number of Contracted functions in a given shell type
integer,intent(in)    :: MaxnCont
!> the start orbital of the atoms with this type
integer,intent(inout) :: startOrbitalOfType(MaxnAtoms,nTypes)
!> the exponents for each type of batches 
real(realk),intent(inout) :: exponentsOfType(MaxnPrim,nTypes)
!> the contraction coefficients for each type of batches 
real(realk),intent(inout) :: ContractCoeffOfType(MaxnPrim,MaxnCont,nTypes)
!> the centers of the atoms with this type
real(realk),intent(inout) :: CentersOfType(3,MaxnAtoms,nTypes)
!
INTEGER,pointer          :: MODELTYPES(:),MODELBATCHTYPES(:,:,:)
INTEGER                  :: maxseg,maxang,I,K,L,iBatchType,R,icharge,type,iseg
INTEGER                  :: nBatchType,nrow,ncol,irow,icol,orbitalindex,nOrbComp
!tmp
integer   :: nAtomsOfType(nTypes)
!assume nofamily
call mem_alloc(MODELTYPES,BASISINFO%natomtypes)

maxSeg = 0
maxAng = 0
DO I=1,BASISINFO%natomtypes
   maxAng = MAX(maxAng,BASISINFO%ATOMTYPE(I)%nAngmom)
   DO K=1,BASISINFO%ATOMTYPE(I)%nAngmom
      maxSeg = MAX(maxseg,BASISINFO%ATOMTYPE(I)%SHELL(K)%nsegments)
   ENDDO
ENDDO

call mem_alloc(MODELBATCHTYPES,maxSeg,maxAng,BASISINFO%natomtypes)
L=0
iBatchType =  0 
DO I=1,BASISINFO%natomtypes
   L=L+1
   MODELTYPES(I)=L
   DO K=1,BASISINFO%ATOMTYPE(I)%nAngmom
      DO iseg = 1,BASISINFO%ATOMTYPE(I)%SHELL(K)%nsegments
         iBatchType = iBatchType + 1 
         MODELBATCHTYPES(iseg,K,L) = iBatchType
      ENDDO
   ENDDO
ENDDO
nBatchType = iBatchType
IF(nBatchType.NE.ntypes)call lsquit('dim mismatch in Determine_nAtomsOfTypes',-1)
!init
do iBatchType=1,nBatchType
   nAtomsOfType(iBatchType) = 0
enddo
orbitalindex = 0
R = BASISINFO%Labelindex
DO I=1,MOLECULE%natoms  
 IF(MOLECULE%ATOM(I)%pointcharge)CYCLE !no basis functions on this point charge
 IF(R.EQ. 0)THEN
    ICHARGE = INT(MOLECULE%ATOM(I)%CHARGE)
    type = BASISINFO%CHARGEINDEX(ICHARGE)
 ELSE
    type = MOLECULE%ATOM(I)%IDtype(R)
 ENDIF
 L=MODELTYPES(type)
 IF(BASISINFO%ATOMTYPE(type)%nAngmom.EQ.0) CYCLE
 DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
  DO iseg = 1,BASISINFO%ATOMTYPE(type)%SHELL(K)%nsegments
   iBatchType = MODELBATCHTYPES(iseg,K,L)
   nrow = BASISINFO%ATOMTYPE(type)%SHELL(K)%segment(iseg)%nrow
   ncol = BASISINFO%ATOMTYPE(type)%SHELL(K)%segment(iseg)%ncol
   do irow = 1,nrow
      exponentsOfType(irow,iBatchType) = &
           & BASISINFO%ATOMTYPE(type)%SHELL(K)%SEGMENT(iseg)%Exponents(irow)
   enddo
   do icol = 1,ncol
    do irow = 1,nrow
       ContractCoeffOfType(irow,icol,iBatchType) = &
          & BASISINFO%ATOMTYPE(type)%SHELL(K)%SEGMENT(iseg)%elms(irow+(icol-1)*nrow)
    enddo
   enddo
   nAtomsOfType(iBatchType) = nAtomsOfType(iBatchType) + 1
   CentersOfType(1,nAtomsOfType(iBatchType),iBatchType) = MOLECULE%ATOM(I)%Center(1)
   CentersOfType(2,nAtomsOfType(iBatchType),iBatchType) = MOLECULE%ATOM(I)%Center(2)
   CentersOfType(3,nAtomsOfType(iBatchType),iBatchType) = MOLECULE%ATOM(I)%Center(3)
   startorbitalOfType(nAtomsOfType(iBatchType),iBatchType) = orbitalindex
   IF (spherical) THEN
      nOrbComp = 2*K-1
   ELSE
      nOrbComp = K*(K+1)/2
   ENDIF
   orbitalindex = orbitalindex + nOrbComp*ncol
  ENDDO
 ENDDO
ENDDO
call mem_dealloc(MODELTYPES)
call mem_dealloc(MODELBATCHTYPES)

end Subroutine Build_TypeInfo2

END MODULE IchorErimoduleHost
