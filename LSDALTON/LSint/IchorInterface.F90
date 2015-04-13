!> @file
!> Contains the main LSDALTON to Ichor integral Interfaces 
!> Ichor is the "Integral Code Hand Optimized for Rapid evaluation" 
!> Ichor is the ethereal golden fluid that is the blood of the greek gods
!> The Ichor code is in the IchorIntegrals directory. 

!> \brief Main Ichor Interface for the calculation of integrals 
!> based on the Obara Saika(OS)/Head-Gordon-Pople(HGP)
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorErimoduleHost
  use precision
  use TYPEDEFTYPE,only: lssetting,BASISSETINFO,MOLECULEINFO
  use memory_handling, only: mem_alloc,mem_dealloc, mem_add_external_memory,&
       & mem_allocated_global
  use basis_typetype,only: RegBasParam,CABBasParam
  use dec_typedef_module, only: DecAObatchinfo
public:: MAIN_ICHORERI_DRIVER, SCREEN_ICHORERI_DRIVER, &
     & determine_MinimumAllowedAObatchSize, &
     & determine_Ichor_nbatchesofAOS, determine_Ichor_batchesofAOS,&
     & determine_Ichor_nAObatches, FREE_SCREEN_ICHORERI,&
     & screen_ichoreri_retrieve_gabdim,screen_ichoreri_retrieve_gab,&
     & MAIN_LINK_ICHORERI_DRIVER, MAIN_ICHORERI_MOTRANS_DRIVER, &
     & write_ichoreri_info, main_ichoreri_readdriver,&
     & determine_Ichor_ActualDim
private
CONTAINS
SUBROUTINE determine_MinimumAllowedAObatchSize(setting,iAO,AOSPEC,MinimumAllowedAObatchSize)
implicit none
character(len=1),intent(in) :: AOspec
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: iAO
integer,intent(inout)     :: MinimumAllowedAObatchSize
!
#ifdef VAR_ICHOR
TYPE(BASISSETINFO),pointer :: AObasis
integer :: nAOBatches
integer,pointer :: OrbSizeOfAOBatches(:)
logical   :: spherical
spherical = .TRUE.
IF(AOspec.EQ.'R')THEN
   !   The regular AO-basis
   AObasis => setting%basis(iAO)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis => setting%basis(iAO)%p%BINFO(CABBasParam)   
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF
call Determine_nBatches(setting%MOLECULE(iAO)%p,AObasis,nAOBatches)
call mem_alloc(OrbSizeOfAOBatches,nAOBatches)
call Determine_OrbSizeOfBatches(setting%MOLECULE(iAO)%p,AObasis,&
     & nAOBatches,OrbSizeOfAOBatches,Spherical)
MinimumAllowedAObatchSize = MAXVAL(OrbSizeOfAOBatches)
call mem_dealloc(OrbSizeOfAOBatches)
#endif
end SUBROUTINE determine_MinimumAllowedAObatchSize

SUBROUTINE determine_Ichor_nbatchesofAOS(setting,iAO,AOSPEC,&
     & RequestedOrbitalDimOfAObatch,nbatchesofAOS,lupri)
implicit none
character(len=1),intent(in) :: AOspec
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: iAO,lupri
integer,intent(inout)     :: nbatchesofAOS
integer,intent(inout)     :: RequestedOrbitalDimOfAObatch
!
#ifdef VAR_ICHOR
TYPE(BASISSETINFO),pointer :: AObasis
integer     :: MaxOrbitalDimOfAObatch,nBatches,MinimumAllowedAObatchSize
integer     :: MinOrbitalDimOfAObatch,J,n,k
real(realk) :: MaxOrbitalDimOfAObatchR,MinOrbitalDimOfAObatchR
real(realk) :: Ratio1
integer,pointer :: OrbSizeOfBatches(:),MaxOrbitalDimOfAObatch2(:)
real(realk),pointer :: Ratio2(:)
logical   :: spherical
spherical = .TRUE.
IF(AOspec.EQ.'R')THEN
   !   The regular AO-basis
   AObasis => setting%basis(iAO)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis => setting%basis(iAO)%p%BINFO(CABBasParam)   
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF
call Determine_nBatches(setting%MOLECULE(iAO)%p,AObasis,nBatches)
call mem_alloc(OrbSizeOfBatches,nBatches)
call Determine_OrbSizeOfBatches(setting%MOLECULE(iAO)%p,AObasis,&
     & nBatches,OrbSizeOfBatches,Spherical)
MinimumAllowedAObatchSize = MAXVAL(OrbSizeOfBatches)
IF(RequestedOrbitalDimOfAObatch.LT.MAXVAL(OrbSizeOfBatches))THEN
   WRITE(lupri,'(A)')'Error In Determine_Ichor_nbatchesofAOS RequestedMaximumOrbitalDimOfAObatch too small'
   WRITE(lupri,'(A,I6)')'RequestedMaximumOrbitalDimOfAObatch have to be at least',MAXVAL(OrbSizeOfBatches)
   WRITE(lupri,'(A)')'You can call determine_MinimumAllowedAObatchSize(setting,spherical,iAO)'
   WRITE(lupri,'(A)')'to determine this size before this call'
   call LSQUIT('Error In Determine_Ichor_nbatchesofAOS RequestedMaximumOrbitalDimOfAObatch too small',-1)
ENDIF

n = RequestedOrbitalDimOfAObatch-MinimumAllowedAObatchSize+1
call mem_alloc(ratio2,n)
call mem_alloc(MaxOrbitalDimOfAObatch2,n)
ratio2 = 0
n=0
DO J=RequestedOrbitalDimOfAObatch,MinimumAllowedAObatchSize,-1
   call loop1(nbatchesofAOS,nBatches,OrbSizeOfBatches,&
        & MaxOrbitalDimOfAObatch,MinOrbitalDimOfAObatch,J)
   n=n+1
   MaxOrbitalDimOfAObatchR = MaxOrbitalDimOfAObatch
   MinOrbitalDimOfAObatchR = MinOrbitalDimOfAObatch
   Ratio2(n) = MaxOrbitalDimOfAObatchR/MinOrbitalDimOfAObatchR
   MaxOrbitalDimOfAObatch2(n) = MaxOrbitalDimOfAObatch
ENDDO
Ratio1 = Ratio2(1)
k=1
DO J=2,n
   IF(ABS(Ratio2(J)-1.0E0_realk).LT.ABS(Ratio1-1.0E0_realk))THEN
      Ratio1 = Ratio2(J)
      k=J
   ELSEIF(ABS(Ratio2(J)-Ratio1).LT.1.0E-8_realk)THEN
      !same Ratio
      IF(RequestedOrbitalDimOfAObatch+1-J.EQ.MaxOrbitalDimOfAObatch2(J))THEN
         !A desireable quanity that the RequestedOrbitalDimOfAObatch=MaxOrbitalDimOfAObatch2
         IF(RequestedOrbitalDimOfAObatch+1-k.EQ.MaxOrbitalDimOfAObatch2(k))THEN
            !This desireable quanity already fulfilled for larger RequestedOrbitalDimOfAObatch
            !do nothing
         ELSE
            Ratio1 = Ratio2(J)
            k=J
         ENDIF
      ENDIF
   ENDIF
ENDDO
RequestedOrbitalDimOfAObatch = RequestedOrbitalDimOfAObatch+1-k
!determine nbatchesofAOS for the right RequestedOrbitalDimOfAObatch
call loop1(nbatchesofAOS,nBatches,OrbSizeOfBatches,&
     & MaxOrbitalDimOfAObatch,MinOrbitalDimOfAObatch,&
     & RequestedOrbitalDimOfAObatch)

call mem_dealloc(ratio2)
call mem_dealloc(MaxOrbitalDimOfAObatch2)
call mem_dealloc(OrbSizeOfBatches)
#endif
end SUBROUTINE determine_Ichor_nbatchesofAOS

SUBROUTINE determine_Ichor_nAObatches(setting,iAO,AOSPEC,nAObatches,lupri)
implicit none
character(len=1),intent(in) :: AOspec
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: iAO,lupri
integer,intent(inout)     :: nAObatches
!
#ifdef VAR_ICHOR
TYPE(BASISSETINFO),pointer :: AObasis
IF(AOspec.EQ.'R')THEN
   !   The regular AO-basis
   AObasis => setting%basis(iAO)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis => setting%basis(iAO)%p%BINFO(CABBasParam)   
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF
call Determine_nBatches(setting%MOLECULE(iAO)%p,AObasis,nAOBatches)
#endif
end SUBROUTINE determine_Ichor_nAObatches

subroutine loop1(nbatchesofAOS,nBatches,OrbSizeOfBatches,&
     & MaxOrbitalDimOfAObatch,MinOrbitalDimOfAObatch,&
     & RequestedOrbitalDimOfAObatch)
implicit none
integer,intent(in) :: RequestedOrbitalDimOfAObatch,nBatches
integer,intent(inout) :: nbatchesofAOS,MaxOrbitalDimOfAObatch
integer,intent(inout) :: MinOrbitalDimOfAObatch 
integer,intent(in) :: OrbSizeOfBatches(nBatches)
!local
#ifdef VAR_ICHOR
integer :: I,DIM
nbatchesofAOS=1
DIM = 0
MaxOrbitalDimOfAObatch = 0
MinOrbitalDimOfAObatch = HUGE(DIM)
DO I=1,nBatches
   IF(DIM+OrbSizeOfBatches(I).LE.RequestedOrbitalDimOfAObatch)THEN
      !BatchOrbitaldimension smaller than allowed
      DIM = DIM + OrbSizeOfBatches(I)
   ELSE
      MaxOrbitalDimOfAObatch = MAX(MaxOrbitalDimOfAObatch,DIM)
      MinOrbitalDimOfAObatch = MIN(MinOrbitalDimOfAObatch,DIM)
      nbatchesofAOS=nbatchesofAOS+1
      DIM = OrbSizeOfBatches(I)      
   ENDIF
ENDDO
MaxOrbitalDimOfAObatch = MAX(MaxOrbitalDimOfAObatch,DIM)
MinOrbitalDimOfAObatch = MIN(MinOrbitalDimOfAObatch,DIM)
#endif
end subroutine loop1

SUBROUTINE determine_Ichor_batchesofAOS(setting,iAO,AOSPEC,&
     & RequestedOrbitalDimOfAObatch,nbatchesofAOS,AObatchinfo,&
     & MaxOrbitalDimOfAObatch,lupri)
implicit none
character(len=1),intent(in) :: AOspec
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: iAO,nbatchesofAOS,lupri
type(DecAObatchinfo)      :: AObatchinfo(nbatchesofAOS)
integer,intent(in)        :: RequestedOrbitalDimOfAObatch
integer,intent(inout)     :: MaxOrbitalDimOfAObatch
!
#ifdef VAR_ICHOR
TYPE(BASISSETINFO),pointer :: AObasis
integer,pointer :: OrbSizeOfBatches(:)
integer :: nBatches,MinimumAllowedAObatchSize,ibatchesofAOS
integer :: I,DIM,ORBINDEX,MinOrbitalDimOfAObatch
logical   :: spherical
spherical = .TRUE.
IF(AOspec.EQ.'R')THEN
   !   The regular AO-basis
   AObasis => setting%basis(iAO)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis => setting%basis(iAO)%p%BINFO(CABBasParam)   
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF
call Determine_nBatches(setting%MOLECULE(iAO)%p,AObasis,nBatches)
call mem_alloc(OrbSizeOfBatches,nBatches)
call Determine_OrbSizeOfBatches(setting%MOLECULE(iAO)%p,AObasis,&
     & nBatches,OrbSizeOfBatches,Spherical)
MinimumAllowedAObatchSize = MAXVAL(OrbSizeOfBatches)
IF(RequestedOrbitalDimOfAObatch.LT.MAXVAL(OrbSizeOfBatches))THEN
   WRITE(lupri,'(A)')'Error In Determine_Ichor_batchesofAOS RequestedMaximumOrbitalDimOfAObatch too small'
   WRITE(lupri,'(A,I6)')'RequestedMaximumOrbitalDimOfAObatch have to be at least',MAXVAL(OrbSizeOfBatches)
   WRITE(lupri,'(A)')'You can call determine_MinimumAllowedAObatchSize(setting,spherical,iAO)'
   WRITE(lupri,'(A)')'to determine this size before this call'
   call LSQUIT('Error In Determine_Ichor_batchesofAOS RequestedMaximumOrbitalDimOfAObatch too small',-1)
ENDIF
ibatchesofAOS=1
AObatchinfo(ibatchesofAOS)%OrbStart = 1
AObatchinfo(ibatchesofAOS)%AOStart = 1
AObatchinfo(ibatchesofAOS)%OrbEnd = 0
AObatchinfo(ibatchesofAOS)%AOEnd = 0
ORBINDEX = 0
DIM = 0 
MaxOrbitalDimOfAObatch = 0
MinOrbitalDimOfAObatch = HUGE(ORBINDEX)
DO I=1,nBatches-1
   IF(DIM+OrbSizeOfBatches(I).LE.RequestedOrbitalDimOfAObatch)THEN
      !BatchOrbitaldimension smaller than allowed
      !we add to Batch of batches
      DIM = DIM + OrbSizeOfBatches(I)
      AObatchinfo(ibatchesofAOS)%OrbEnd = AObatchinfo(ibatchesofAOS)%OrbEnd+OrbSizeOfBatches(I)
      AObatchinfo(ibatchesofAOS)%AOEnd = AObatchinfo(ibatchesofAOS)%AOEnd+1
      AObatchinfo(ibatchesofAOS)%DIM = DIM
   ELSE
      !we start a new batch of batches
      MaxOrbitalDimOfAObatch = MAX(MaxOrbitalDimOfAObatch,DIM)
      MinOrbitalDimOfAObatch = MIN(MinOrbitalDimOfAObatch,DIM)
      ibatchesofAOS=ibatchesofAOS+1
      DIM = OrbSizeOfBatches(I)
      AObatchinfo(ibatchesofAOS)%OrbStart = ORBINDEX+1      
      AObatchinfo(ibatchesofAOS)%OrbEnd = ORBINDEX+OrbSizeOfBatches(I)
      AObatchinfo(ibatchesofAOS)%AOStart = I
      AObatchinfo(ibatchesofAOS)%AOEnd = I
      AObatchinfo(ibatchesofAOS)%DIM = DIM
   ENDIF
   ORBINDEX = ORBINDEX+OrbSizeOfBatches(I)
ENDDO
I=nBatches
IF(DIM+OrbSizeOfBatches(I).LE.RequestedOrbitalDimOfAObatch)THEN
   !BatchOrbitaldimension smaller than allowed
   !we add to final batch to current ibatchesofAOS
   DIM = DIM + OrbSizeOfBatches(I)
   AObatchinfo(ibatchesofAOS)%OrbEnd = AObatchinfo(ibatchesofAOS)%OrbEnd+OrbSizeOfBatches(I)
   AObatchinfo(ibatchesofAOS)%AOEnd = AObatchinfo(ibatchesofAOS)%AOEnd+1
   AObatchinfo(ibatchesofAOS)%DIM = DIM
ELSE
   !we start a new batch of just this final batch
   ibatchesofAOS=ibatchesofAOS+1
   DIM = OrbSizeOfBatches(I)
   AObatchinfo(ibatchesofAOS)%DIM = OrbSizeOfBatches(I)
   AObatchinfo(ibatchesofAOS)%OrbStart = ORBINDEX+1      
   AObatchinfo(ibatchesofAOS)%AOStart = I
   AObatchinfo(ibatchesofAOS)%OrbEnd = ORBINDEX+OrbSizeOfBatches(I)
   AObatchinfo(ibatchesofAOS)%AOEnd = I
ENDIF
MaxOrbitalDimOfAObatch = MAX(MaxOrbitalDimOfAObatch,DIM)
MinOrbitalDimOfAObatch = MIN(MinOrbitalDimOfAObatch,DIM)
call mem_dealloc(OrbSizeOfBatches)
#endif
end SUBROUTINE determine_Ichor_batchesofAOS

SUBROUTINE determine_Ichor_ActualDim(setting,iAO,AOSPEC,&
     & RequestedOrbitalDimOfAObatch,MaxOrbitalDimOfAObatch,lupri)
implicit none
character(len=1),intent(in) :: AOspec
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: iAO,lupri
integer,intent(in)        :: RequestedOrbitalDimOfAObatch
integer,intent(inout)     :: MaxOrbitalDimOfAObatch
!
#ifdef VAR_ICHOR
TYPE(BASISSETINFO),pointer :: AObasis
integer,pointer :: OrbSizeOfBatches(:)
integer :: nBatches,MinimumAllowedAObatchSize,ibatchesofAOS
integer :: I,DIM,ORBINDEX,MinOrbitalDimOfAObatch
logical   :: spherical
spherical = .TRUE.
IF(AOspec.EQ.'R')THEN
   !   The regular AO-basis
   AObasis => setting%basis(iAO)%p%BINFO(RegBasParam)
ELSEIF(AOspec.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis => setting%basis(iAO)%p%BINFO(CABBasParam)   
ELSE
   call lsquit('Unknown specification in build_batchesOfAOs',-1)
ENDIF
call Determine_nBatches(setting%MOLECULE(iAO)%p,AObasis,nBatches)
call mem_alloc(OrbSizeOfBatches,nBatches)
call Determine_OrbSizeOfBatches(setting%MOLECULE(iAO)%p,AObasis,&
     & nBatches,OrbSizeOfBatches,Spherical)
MinimumAllowedAObatchSize = MAXVAL(OrbSizeOfBatches)
IF(RequestedOrbitalDimOfAObatch.LT.MAXVAL(OrbSizeOfBatches))THEN
   WRITE(lupri,'(A)')'Error In Determine_Ichor_batchesofAOS RequestedMaximumOrbitalDimOfAObatch too small'
   WRITE(lupri,'(A,I6)')'RequestedMaximumOrbitalDimOfAObatch have to be at least',MAXVAL(OrbSizeOfBatches)
   WRITE(lupri,'(A)')'You can call determine_MinimumAllowedAObatchSize(setting,spherical,iAO)'
   WRITE(lupri,'(A)')'to determine this size before this call'
   call LSQUIT('Error In Determine_Ichor_batchesofAOS RequestedMaximumOrbitalDimOfAObatch too small',-1)
ENDIF
ibatchesofAOS=1
DIM = 0 
MaxOrbitalDimOfAObatch = 0
DO I=1,nBatches-1
   IF(DIM+OrbSizeOfBatches(I).LE.RequestedOrbitalDimOfAObatch)THEN
      DIM = DIM + OrbSizeOfBatches(I)
   ELSE
      MaxOrbitalDimOfAObatch = MAX(MaxOrbitalDimOfAObatch,DIM)
      DIM = OrbSizeOfBatches(I)
   ENDIF
ENDDO
I=nBatches
IF(DIM+OrbSizeOfBatches(I).LE.RequestedOrbitalDimOfAObatch)THEN
   !BatchOrbitaldimension smaller than allowed
   !we add to final batch to current ibatchesofAOS
   DIM = DIM + OrbSizeOfBatches(I)
ELSE
   !we start a new batch of just this final batch
   ibatchesofAOS=ibatchesofAOS+1
   DIM = OrbSizeOfBatches(I)
ENDIF
MaxOrbitalDimOfAObatch = MAX(MaxOrbitalDimOfAObatch,DIM)
call mem_dealloc(OrbSizeOfBatches)

#endif
end SUBROUTINE determine_Ichor_ActualDim

!dim1,dim2,dim3,dim4 are the AO dimensions 
SUBROUTINE MAIN_ICHORERI_MOTRANS_DRIVER(LUPRI,IPRINT,setting,dim1,dim2,dim3,dim4,integrals,intspec,FullBatch,&
     & nbatchAstart,nbatchAend,nbatchBstart,nbatchBend,nbatchCstart,nbatchCend,nbatchDstart,nbatchDend,&
     & Cocc,Cvirt,nOcc,nVirt,intThreshold)
implicit none
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: LUPRI,IPRINT,dim1,dim2,dim3,dim4,nOcc,nVirt
logical,intent(in)        :: FullBatch  !full batches are assumed and the nbatchXXX arguments are not used 
integer,intent(in)        :: nbatchAstart,nbatchAend,nbatchBstart,nbatchBend
integer,intent(in)        :: nbatchCstart,nbatchCend,nbatchDstart,nbatchDend
real(realk),intent(in)    :: Cvirt(dim1,nVirt)
real(realk),intent(in)    :: Cocc(dim2,nOcc),intThreshold
real(realk),intent(inout) :: integrals(nVirt,nOcc,nVirt,nOcc)
Character,intent(IN)      :: intSpec(5)
!
#ifdef VAR_ICHOR
logical :: MoTrans,NoSymmetry
integer :: A,I,B,J
!DO J=1,nOcc
!   DO B=1,nVirt
!      DO I=1,nOcc
!         DO A=1,nVirt
!            integrals(A,I,B,J) = 0.0E0_realk
!         ENDDO
!      ENDDO
!   ENDDO
!ENDDO
MoTrans = .TRUE.
NoSymmetry = .FALSE.
!all these subroutines are in the file IchorIntegrals/MainIchorInterface.F90
call InitIchorInputInfo
call IchorInputM2(Cvirt,dim1,nVirt,Cocc,dim2,nOcc)
call IchorInputSpec(1,2,1,2)
CALL MAIN_ICHORERI_DRIVER(LUPRI,IPRINT,setting,dim1,dim2,dim3,dim4,integrals,intspec,FullBatch,&
     & nbatchAstart,nbatchAend,nbatchBstart,nbatchBend,nbatchCstart,nbatchCend,nbatchDstart,&
     & nbatchDend,MoTrans,nVirt,nOcc,nVirt,nOcc,NoSymmetry,intThreshold)
call FreeIchorInputInfo()
#endif
END SUBROUTINE MAIN_ICHORERI_MOTRANS_DRIVER

SUBROUTINE MAIN_ICHORERI_DRIVER(LUPRI,IPRINT,setting,dim1,dim2,dim3,dim4,integrals,intspec,FullBatch,&
     & nbatchAstart,nbatchAend,nbatchBstart,nbatchBend,nbatchCstart,nbatchCend,nbatchDstart,nbatchDend,&
     & MoTrans,OutDim1,OutDim2,OutDim3,OutDim4,NoSymmetry,intThreshold)
implicit none
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: LUPRI,IPRINT
integer,intent(in)        :: dim1,dim2,dim3,dim4 !AO dim
integer,intent(in)        :: OutDim1,OutDim2,OutDim3,OutDim4 !output dim
logical,intent(in)        :: FullBatch!full batches are assumed and the nbatchXXX arguments are not used 
logical,intent(in)        :: MoTrans,NoSymmetry
integer,intent(in)        :: nbatchAstart,nbatchAend,nbatchBstart,nbatchBend
integer,intent(in)        :: nbatchCstart,nbatchCend,nbatchDstart,nbatchDend
real(realk),intent(inout) :: integrals(OutDim1,OutDim2,OutDim3,OutDim4)
real(realk),intent(in)    :: intThreshold
Character,intent(IN)      :: intSpec(5)
!
#ifdef VAR_ICHOR
integer                :: nTypes
!A
integer :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchA
integer,pointer     :: nAtomsOfTypeA(:),nPrimOfTypeA(:)
integer,pointer     :: nContOfTypeA(:),AngmomOfTypeA(:),startOrbitalOfTypeA(:,:)
real(realk),pointer :: exponentsOfTypeA(:,:),Acenters(:,:,:),ContractCoeffOfTypeA(:,:,:)
!B
integer :: nTypesB,MaxnAtomsB,MaxnPrimB,MaxnContB,nbatchB
integer,pointer     :: nAtomsOfTypeB(:),nPrimOfTypeB(:)
integer,pointer     :: nContOfTypeB(:),AngmomOfTypeB(:),startOrbitalOfTypeB(:,:)
real(realk),pointer :: exponentsOfTypeB(:,:),Bcenters(:,:,:),ContractCoeffOfTypeB(:,:,:)
!C
integer :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC,nbatchC
integer,pointer     :: nAtomsOfTypeC(:),nPrimOfTypeC(:)
integer,pointer     :: nContOfTypeC(:),AngmomOfTypeC(:),startOrbitalOfTypeC(:,:)
real(realk),pointer :: exponentsOfTypeC(:,:),Ccenters(:,:,:),ContractCoeffOfTypeC(:,:,:)
!D
integer :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD,nbatchD
integer,pointer     :: nAtomsOfTypeD(:),nPrimOfTypeD(:)
integer,pointer     :: nContOfTypeD(:),AngmomOfTypeD(:),startOrbitalOfTypeD(:,:)
real(realk),pointer :: exponentsOfTypeD(:,:),Dcenters(:,:,:),ContractCoeffOfTypeD(:,:,:)
!job specification
integer :: SphericalSpec,IchorJobSpec,IchorInputSpec,IchorParSpec,IchorScreenSpec
Integer :: IchorInputDim1,IchorInputDim2,IchorInputDim3,IchorDebugSpec,IchorAlgoSpec
integer :: filestorageIdentifier,MaxFileStorage,IchorPermuteSpec
Integer :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5
Integer :: GabIdentifier, IchorGabID1, IchorGabID2
logical :: SameLHSaos
real(realk) :: THRESHOLD_CS,THRESHOLD_QQR,THRESHOLD_OD
real(realk),pointer :: InputStorage(:)
real(realk),pointer :: BATCHGAB(:),BATCHGCD(:)
Integer(kind=long) :: MaxMem,MaxMemAllocated,MemAllocated
Integer :: nBatchesA,nBatchesB,nBatchesC,nBatchesD
logical :: spherical
TYPE(BASISSETINFO),pointer :: AObasis
integer :: nbatchAstart2,nbatchAend2,nbatchBstart2,nbatchBend2
integer :: nbatchCstart2,nbatchCend2,nbatchDstart2,nbatchDend2
integer :: IchorOperatorSpec
logical :: SameRHSaos,SameODs,CRIT1,CRIT2,CRIT3,CRIT4,doLink,rhsDmat,CRIT5
logical :: ForceCPU,ForceGPU
!FULLABATCH,FULLBBATCH,FULLCBATCH,FULLDBATCH
spherical = .TRUE.

!A
Call BuildCenterAndTypeInfo(1,intSpec(1),setting,ntypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,nbatchAstart2,nbatchAend2,MaxnAtomsA,MaxnPrimA,MaxnContA,&
     & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,FullBatch,&
     & nbatchAstart,nbatchAend,spherical)
!B
Call BuildCenterAndTypeInfo(2,intSpec(2),setting,ntypesB,nBatchesB,nAtomsOfTypeB,AngmomOfTypeB,&
     & nPrimOfTypeB,nContOfTypeB,nbatchBstart2,nbatchBend2,MaxnAtomsB,MaxnPrimB,MaxnContB,&
     & startOrbitalOfTypeB,exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters,FullBatch,&
     & nbatchBstart,nbatchBend,spherical)
!C
Call BuildCenterAndTypeInfo(3,intSpec(3),setting,ntypesC,nBatchesC,nAtomsOfTypeC,AngmomOfTypeC,&
     & nPrimOfTypeC,nContOfTypeC,nbatchCstart2,nbatchCend2,MaxnAtomsC,MaxnPrimC,MaxnContC,&
     & startOrbitalOfTypeC,exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters,FullBatch,&
     & nbatchCstart,nbatchCend,spherical)
!D
Call BuildCenterAndTypeInfo(4,intSpec(4),setting,ntypesD,nBatchesD,nAtomsOfTypeD,AngmomOfTypeD,&
     & nPrimOfTypeD,nContOfTypeD,nbatchDstart2,nbatchDend2,MaxnAtomsD,MaxnPrimD,MaxnContD,&
     & startOrbitalOfTypeD,exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters,FullBatch,&
     & nbatchDstart,nbatchDend,spherical)

call GetIchorSphericalParamIdentifier(SphericalSpec)
IF(MoTrans)THEN
   call GetIchorJobMOtransIdentifier(IchorJobSpec)   
ELSE
   doLink = .FALSE.
   call GetIchorJobEriIdentifier(IchorJobSpec,doLink)
ENDIF
rhsDmat = .FALSE. !no rhs density matrix supplied as input 
call GetIchorInputIdentifier(IchorInputSpec,rhsDmat)
IchorInputDim1=1                 !not used since   IcorInputNoInput
IchorInputDim2=1                 !not used since   IcorInputNoInput
IchorInputDim3=1                 !not used since   IcorInputNoInput
call mem_alloc(InputStorage,1)
call GetIchorParallelSpecIdentifier(IchorParSpec)   !no parallelization
call GetIchorDebugIdentifier(IchorDebugSpec,iprint) !Debug PrintLevel
call GetIchorAlgorithmSpecIdentifier(IchorAlgoSpec)
call SetIchorGPUMaxMem(setting%GPUMAXMEM)

IF(FullBatch)THEN
   SameLHSaos = (intSpec(1).EQ.intSpec(2).AND.Setting%sameMol(1,2)).AND.Setting%sameBas(1,2)
   SameRHSaos = (intSpec(3).EQ.intSpec(4).AND.Setting%sameMol(3,4)).AND.Setting%sameBas(3,4)
   SameODs = ((intSpec(1).EQ.intSpec(3)).AND.(intSpec(2).EQ.intSpec(4))&
        & .AND.(Setting%sameMol(1,3).AND.Setting%sameMol(3,4)))&
        &.AND.(Setting%sameBAS(1,3).AND.Setting%sameBAS(3,4))
!   SameLHSaos = intSpec(1).EQ.intSpec(2)
!   SameRHSaos = intSpec(3).EQ.intSpec(4)
!   SameODs = (intSpec(1).EQ.intSpec(3)).AND.(intSpec(2).EQ.intSpec(4))
ELSE
!   FULLABATCH = nbatchAstart2.EQ.1.AND.nbatchAend2.EQ.nBatchesA
!   FULLBBATCH = nbatchBstart2.EQ.1.AND.nbatchBend2.EQ.nBatchesB
!   FULLCBATCH = nbatchCstart2.EQ.1.AND.nbatchCend2.EQ.nBatchesC
!   FULLDBATCH = nbatchDstart2.EQ.1.AND.nbatchDend2.EQ.nBatchesD
   SameLHSaos = (intSpec(1).EQ.intSpec(2).AND.Setting%sameMol(1,2)).AND.&
        & ((nbatchAstart2.EQ.nbatchBstart2).AND.(nbatchAend2.EQ.nbatchBend2))
   SameLHSaos = SameLHSaos.AND.Setting%sameBas(1,2)
   SameRHSaos = (intSpec(3).EQ.intSpec(4).AND.Setting%sameMol(3,4)).AND.&
        & ((nbatchCstart2.EQ.nbatchDstart2).AND.(nbatchCend2.EQ.nbatchDend2))
   SameRHSaos = SameRHSaos.AND.Setting%sameBas(3,4)
   CRIT1 = (intSpec(1).EQ.intSpec(3)).AND.(intSpec(2).EQ.intSpec(4))
   CRIT2 = Setting%sameMol(1,3).AND.Setting%sameMol(2,4)
   CRIT3 = (nbatchAstart2.EQ.nbatchCstart2).AND.(nbatchAend2.EQ.nbatchCend2)
   CRIT4 = (nbatchBstart2.EQ.nbatchDstart2).AND.(nbatchBend2.EQ.nbatchDend2)
   CRIT5 = Setting%sameBas(1,3).AND.Setting%sameBas(2,4)
   SameODs = ((CRIT1.AND.CRIT2).AND.(CRIT3.AND.CRIT4)).AND.CRIT5

ENDIF
IF(NoSymmetry)THEN
   SameLHSaos = .FALSE. 
   SameRHSaos = .FALSE. 
   SameODs = .FALSE. 
ENDIF
call GetIchorPermuteParameter(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
call GetIchorFileStorageIdentifier(filestorageIdentifier)
call GetIchorOpereratorIntSpec(intSpec(5),IchorOperatorSpec)

MaxMem=0         !Maximum Memory Ichor is allowed to use. Zero = no restrictions
MaxFileStorage=0 !Maximum File size, if zero - no file will be written or read. 
MaxMemAllocated=0!Maximum Memory used in the program. Ichor adds to this value
MemAllocated = 0 !Memory allocated in the Ichor program

call GetIchorScreeningParameter(IchorScreenSpec,SETTING%SCHEME%CS_SCREEN,&
     & SETTING%SCHEME%OD_SCREEN,.FALSE.)
IF(SETTING%SCHEME%CS_SCREEN.AND.SETTING%SCHEME%OD_SCREEN)THEN
   CALL GET_IchorGabIDInterface(IchorGabID1,IchorGabID2)
ELSEIF(SETTING%SCHEME%CS_SCREEN)THEN
   CALL GET_IchorGabIDInterface(IchorGabID1,IchorGabID2)
ELSE   
   IchorGabID1=0 !screening Matrix Identifier, not used if IchorScreenNone
   IchorGabID2=0 !screening Matrix Identifier, not used if IchorScreenNone
ENDIF

OutputDim1=Outdim1
OutputDim2=Outdim2
OutputDim3=Outdim3
OutputDim4=Outdim4
OutputDim5=1
THRESHOLD_OD = intThreshold*1.0E-1_realk
THRESHOLD_CS = intThreshold
THRESHOLD_QQR = intThreshold
ForceCPU = SETTING%SCHEME%IchorForceCPU
ForceGPU = SETTING%SCHEME%IchorForceGPU
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
     & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorOperatorSpec,&
     & IchorInputDim1,IchorInputDim2,IchorInputDim3,&
     & InputStorage,IchorParSpec,IchorScreenSpec,THRESHOLD_OD,THRESHOLD_CS,&
     & THRESHOLD_QQR,IchorGabID1,IchorGabID2,IchorDebugSpec,&
     & IchorAlgoSpec,IchorPermuteSpec,filestorageIdentifier,MaxMem,&
     & MaxFileStorage,MaxMemAllocated,MemAllocated,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
     & integrals,ForceCPU,ForceGPU,lupri)

call Mem_Add_external_memory(MaxMemAllocated)
call mem_dealloc(InputStorage)
!=====================================================================


!=====================================================================
!free space
call FreeCenterAndTypeInfo(nAtomsOfTypeA,AngmomOfTypeA,&
           & nPrimOfTypeA,nContOfTypeA,startOrbitalOfTypeA,&
           & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeB,AngmomOfTypeB,&
           & nPrimOfTypeB,nContOfTypeB,startOrbitalOfTypeB,&
           & exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeC,AngmomOfTypeC,&
           & nPrimOfTypeC,nContOfTypeC,startOrbitalOfTypeC,&
           & exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeD,AngmomOfTypeD,&
           & nPrimOfTypeD,nContOfTypeD,startOrbitalOfTypeD,&
           & exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters)
#endif
END SUBROUTINE MAIN_ICHORERI_DRIVER

subroutine GetIchorOpereratorIntSpec(intSpec,IchorOperatorSpec)
  implicit none
  character :: intspec 
  integer,intent(inout) :: IchorOperatorSpec
#ifdef VAR_ICHOR
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
     call lsquit('Error in specification of operator in GetIchorOpereratorIntSpec',-1)
  ENDIF
#endif
end subroutine GetIchorOpereratorIntSpec

SUBROUTINE MAIN_ICHORERIMEM_DRIVER(LUPRI,IPRINT,setting,dim1,dim2,dim3,dim4,integrals,intspec,FullBatch,&
     & nbatchAstart,nbatchAend,nbatchBstart,nbatchBend,nbatchCstart,nbatchCend,nbatchDstart,nbatchDend,&
     & MoTrans,OutDim1,OutDim2,OutDim3,OutDim4,NoSymmetry,MaxMemoryUsage,intThreshold)
implicit none
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: LUPRI,IPRINT
integer,intent(in)        :: dim1,dim2,dim3,dim4 !AO dim
integer,intent(in)        :: OutDim1,OutDim2,OutDim3,OutDim4 !output dim
logical,intent(in)        :: FullBatch!full batches are assumed and the nbatchXXX arguments are not used 
logical,intent(in)        :: MoTrans,NoSymmetry
integer,intent(in)        :: nbatchAstart,nbatchAend,nbatchBstart,nbatchBend
integer,intent(in)        :: nbatchCstart,nbatchCend,nbatchDstart,nbatchDend
real(realk),intent(inout) :: integrals(OutDim1,OutDim2,OutDim3,OutDim4)
Character,intent(IN)      :: intSpec(5)
integer(kind=long),intent(inout) :: MaxMemoryUsage !only actual output  
real(realk),intent(in) :: intThreshold
!
#ifdef VAR_ICHOR
integer                :: nTypes
!A
integer :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchA
integer,pointer     :: nAtomsOfTypeA(:),nPrimOfTypeA(:)
integer,pointer     :: nContOfTypeA(:),AngmomOfTypeA(:),startOrbitalOfTypeA(:,:)
real(realk),pointer :: exponentsOfTypeA(:,:),Acenters(:,:,:),ContractCoeffOfTypeA(:,:,:)
!B
integer :: nTypesB,MaxnAtomsB,MaxnPrimB,MaxnContB,nbatchB
integer,pointer     :: nAtomsOfTypeB(:),nPrimOfTypeB(:)
integer,pointer     :: nContOfTypeB(:),AngmomOfTypeB(:),startOrbitalOfTypeB(:,:)
real(realk),pointer :: exponentsOfTypeB(:,:),Bcenters(:,:,:),ContractCoeffOfTypeB(:,:,:)
!C
integer :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC,nbatchC
integer,pointer     :: nAtomsOfTypeC(:),nPrimOfTypeC(:)
integer,pointer     :: nContOfTypeC(:),AngmomOfTypeC(:),startOrbitalOfTypeC(:,:)
real(realk),pointer :: exponentsOfTypeC(:,:),Ccenters(:,:,:),ContractCoeffOfTypeC(:,:,:)
!D
integer :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD,nbatchD
integer,pointer     :: nAtomsOfTypeD(:),nPrimOfTypeD(:)
integer,pointer     :: nContOfTypeD(:),AngmomOfTypeD(:),startOrbitalOfTypeD(:,:)
real(realk),pointer :: exponentsOfTypeD(:,:),Dcenters(:,:,:),ContractCoeffOfTypeD(:,:,:)
!job specification
integer :: SphericalSpec,IchorJobSpec,IchorInputSpec,IchorParSpec,IchorScreenSpec
Integer :: IchorInputDim1,IchorInputDim2,IchorInputDim3,IchorDebugSpec,IchorAlgoSpec
integer :: filestorageIdentifier,MaxFileStorage,IchorPermuteSpec
Integer :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5
Integer :: GabIdentifier, IchorGabID1, IchorGabID2
logical :: SameLHSaos
real(realk) :: THRESHOLD_CS,THRESHOLD_QQR,THRESHOLD_OD
real(realk),pointer :: InputStorage(:)
real(realk),pointer :: BATCHGAB(:),BATCHGCD(:)
Integer(kind=long) :: MaxMem,MaxMemAllocated,MemAllocated
Integer :: nBatchesA,nBatchesB,nBatchesC,nBatchesD,IchorOperatorSpec
logical :: spherical
TYPE(BASISSETINFO),pointer :: AObasis
integer :: nbatchAstart2,nbatchAend2,nbatchBstart2,nbatchBend2
integer :: nbatchCstart2,nbatchCend2,nbatchDstart2,nbatchDend2
logical :: SameRHSaos,SameODs,CRIT1,CRIT2,CRIT3,CRIT4,doLink,rhsDmat,CRIT5
integer(kind=long) :: mem_allocated_global_current
logical :: ForceCPU,ForceGPU
!FULLABATCH,FULLBBATCH,FULLCBATCH,FULLDBATCH
spherical = .TRUE.
!IF (intSpec(5).NE.'C') CALL LSQUIT('MAIN_ICHORERI_DRIVER limited to Coulomb Integrals for now',-1)
mem_allocated_global_current = mem_allocated_global
!A
Call BuildCenterAndTypeInfo(1,intSpec(1),setting,ntypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,nbatchAstart2,nbatchAend2,MaxnAtomsA,MaxnPrimA,MaxnContA,&
     & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,FullBatch,&
     & nbatchAstart,nbatchAend,spherical)
!B
Call BuildCenterAndTypeInfo(2,intSpec(2),setting,ntypesB,nBatchesB,nAtomsOfTypeB,AngmomOfTypeB,&
     & nPrimOfTypeB,nContOfTypeB,nbatchBstart2,nbatchBend2,MaxnAtomsB,MaxnPrimB,MaxnContB,&
     & startOrbitalOfTypeB,exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters,FullBatch,&
     & nbatchBstart,nbatchBend,spherical)
!C
Call BuildCenterAndTypeInfo(3,intSpec(3),setting,ntypesC,nBatchesC,nAtomsOfTypeC,AngmomOfTypeC,&
     & nPrimOfTypeC,nContOfTypeC,nbatchCstart2,nbatchCend2,MaxnAtomsC,MaxnPrimC,MaxnContC,&
     & startOrbitalOfTypeC,exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters,FullBatch,&
     & nbatchCstart,nbatchCend,spherical)
!D
Call BuildCenterAndTypeInfo(4,intSpec(4),setting,ntypesD,nBatchesD,nAtomsOfTypeD,AngmomOfTypeD,&
     & nPrimOfTypeD,nContOfTypeD,nbatchDstart2,nbatchDend2,MaxnAtomsD,MaxnPrimD,MaxnContD,&
     & startOrbitalOfTypeD,exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters,FullBatch,&
     & nbatchDstart,nbatchDend,spherical)

call GetIchorSphericalParamIdentifier(SphericalSpec)
IF(MoTrans)THEN
   call GetIchorJobMOtransIdentifier(IchorJobSpec)   
ELSE
   doLink = .FALSE.
   call GetIchorJobEriIdentifier(IchorJobSpec,doLink)
ENDIF
rhsDmat = .FALSE. !no rhs density matrix supplied as input 
call GetIchorInputIdentifier(IchorInputSpec,rhsDmat)
IchorInputDim1=1                 !not used since   IcorInputNoInput
IchorInputDim2=1                 !not used since   IcorInputNoInput
IchorInputDim3=1                 !not used since   IcorInputNoInput
call mem_alloc(InputStorage,1)
call GetIchorParallelSpecIdentifier(IchorParSpec)   !no parallelization
call GetIchorDebugIdentifier(IchorDebugSpec,iprint) !Debug PrintLevel
call GetIchorAlgorithmSpecIdentifier(IchorAlgoSpec)
call SetIchorGPUMaxMem(setting%GPUMAXMEM)
IF(FullBatch)THEN
   SameLHSaos = (intSpec(1).EQ.intSpec(2).AND.Setting%sameMol(1,2)).AND.Setting%sameBas(1,2)
   SameRHSaos = (intSpec(3).EQ.intSpec(4).AND.Setting%sameMol(3,4)).AND.Setting%sameBas(3,4)
   SameODs = ((intSpec(1).EQ.intSpec(3)).AND.(intSpec(2).EQ.intSpec(4))&
        & .AND.(Setting%sameMol(1,3).AND.Setting%sameMol(3,4)))&
        &.AND.(Setting%sameBAS(1,3).AND.Setting%sameBAS(3,4))
!   SameLHSaos = intSpec(1).EQ.intSpec(2)
!   SameRHSaos = intSpec(3).EQ.intSpec(4)
!   SameODs = (intSpec(1).EQ.intSpec(3)).AND.(intSpec(2).EQ.intSpec(4))
ELSE
!   FULLABATCH = nbatchAstart2.EQ.1.AND.nbatchAend2.EQ.nBatchesA
!   FULLBBATCH = nbatchBstart2.EQ.1.AND.nbatchBend2.EQ.nBatchesB
!   FULLCBATCH = nbatchCstart2.EQ.1.AND.nbatchCend2.EQ.nBatchesC
!   FULLDBATCH = nbatchDstart2.EQ.1.AND.nbatchDend2.EQ.nBatchesD
   SameLHSaos = (intSpec(1).EQ.intSpec(2).AND.Setting%sameMol(1,2)).AND.&
        & ((nbatchAstart2.EQ.nbatchBstart2).AND.(nbatchAend2.EQ.nbatchBend2))
   SameLHSaos = SameLHSaos.AND.Setting%sameBas(1,2)
   SameRHSaos = (intSpec(3).EQ.intSpec(4).AND.Setting%sameMol(3,4)).AND.&
        & ((nbatchCstart2.EQ.nbatchDstart2).AND.(nbatchCend2.EQ.nbatchDend2))
   SameRHSaos = SameRHSaos.AND.Setting%sameBas(3,4)
   CRIT1 = (intSpec(1).EQ.intSpec(3)).AND.(intSpec(2).EQ.intSpec(4))
   CRIT2 = Setting%sameMol(1,3).AND.Setting%sameMol(2,4)
   CRIT3 = (nbatchAstart2.EQ.nbatchCstart2).AND.(nbatchAend2.EQ.nbatchCend2)
   CRIT4 = (nbatchBstart2.EQ.nbatchDstart2).AND.(nbatchBend2.EQ.nbatchDend2)
   CRIT5 = Setting%sameBas(1,3).AND.Setting%sameBas(2,4)
   SameODs = ((CRIT1.AND.CRIT2).AND.(CRIT3.AND.CRIT4)).AND.CRIT5
ENDIF
IF(NoSymmetry)THEN
   SameLHSaos = .FALSE. 
   SameRHSaos = .FALSE. 
   SameODs = .FALSE. 
ENDIF
call GetIchorPermuteParameter(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
call GetIchorFileStorageIdentifier(filestorageIdentifier)
call GetIchorOpereratorIntSpec(intSpec(5),IchorOperatorSpec)

MaxMem=0         !Maximum Memory Ichor is allowed to use. Zero = no restrictions
MaxFileStorage=0 !Maximum File size, if zero - no file will be written or read. 
MaxMemAllocated=0!Maximum Memory used in the program. Ichor adds to this value
MemAllocated = 0 !Memory allocated in the Ichor program

call GetIchorScreeningParameter(IchorScreenSpec,SETTING%SCHEME%CS_SCREEN,&
     & SETTING%SCHEME%OD_SCREEN,.FALSE.)
IF(SETTING%SCHEME%CS_SCREEN.AND.SETTING%SCHEME%OD_SCREEN)THEN
   CALL GET_IchorGabIDInterface(IchorGabID1,IchorGabID2)
ELSEIF(SETTING%SCHEME%CS_SCREEN)THEN
   CALL GET_IchorGabIDInterface(IchorGabID1,IchorGabID2)
ELSE   
   IchorGabID1=0 !screening Matrix Identifier, not used if IchorScreenNone
   IchorGabID2=0 !screening Matrix Identifier, not used if IchorScreenNone
ENDIF

OutputDim1=Outdim1
OutputDim2=Outdim2
OutputDim3=Outdim3
OutputDim4=Outdim4
OutputDim5=1
THRESHOLD_OD = intThreshold*1.0E-1_realk
THRESHOLD_CS = intThreshold
THRESHOLD_QQR = intThreshold
ForceCPU = SETTING%SCHEME%IchorForceCPU
ForceGPU = SETTING%SCHEME%IchorForceGPU
!print*,'THRESHOLD_CS',THRESHOLD_CS,'THRESHOLD_OD',THRESHOLD_OD
!=====================================================================
!  Main Call
!=====================================================================
call IchorEriMemInterface(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
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
     & integrals,ForceCPU,ForceGPU,lupri)

MaxMemoryUsage = MaxMemAllocated + (mem_allocated_global-mem_allocated_global_current)
print*,'MaxMemoryUsage',MaxMemoryUsage

if (MaxMemoryUsage < 100_long) then !Divide by 1 to typecast real
   WRITE(LUPRI,'(A, f10.3,A)') 'MaxMemoryUsage ',MaxMemoryUsage/1.0E0_realk,' Byte'
else if (MaxMemoryUsage < 1000000_long) then
   WRITE(LUPRI,'(A, f10.3,A)') 'MaxMemoryUsage ',MaxMemoryUsage*1.0E-3_realk,' kB'
else if (MaxMemoryUsage < 1000000000_long) then
   WRITE(LUPRI,'(A, f10.3,A)') 'MaxMemoryUsage ',MaxMemoryUsage*1.E-6_realk,' MB'
else
   WRITE(LUPRI,'(A, f10.3,A)') 'MaxMemoryUsage ',MaxMemoryUsage*1.E-9_realk,' GB'
endif

!call Mem_Add_external_memory(MaxMemAllocated)
call mem_dealloc(InputStorage)
!=====================================================================


!=====================================================================
!free space
call FreeCenterAndTypeInfo(nAtomsOfTypeA,AngmomOfTypeA,&
           & nPrimOfTypeA,nContOfTypeA,startOrbitalOfTypeA,&
           & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeB,AngmomOfTypeB,&
           & nPrimOfTypeB,nContOfTypeB,startOrbitalOfTypeB,&
           & exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeC,AngmomOfTypeC,&
           & nPrimOfTypeC,nContOfTypeC,startOrbitalOfTypeC,&
           & exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeD,AngmomOfTypeD,&
           & nPrimOfTypeD,nContOfTypeD,startOrbitalOfTypeD,&
           & exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters)
#endif
END SUBROUTINE MAIN_ICHORERIMEM_DRIVER

Subroutine BuildCenterAndTypeInfo(Center,intSpecA,setting,ntypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,nbatchAstart2,nbatchAend2,MaxnAtomsA,MaxnPrimA,MaxnContA,&
     & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,FullBatch,&
     & nbatchAstart,nbatchAend,spherical)
implicit none
integer,intent(IN)      :: Center,nbatchAstart,nbatchAend
Character,intent(IN)    :: intSpecA
TYPE(lssetting),intent(in):: setting
logical,intent(in) :: FullBatch,spherical
integer,intent(inout) :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchesA
integer,intent(inout) :: nbatchAstart2,nbatchAend2
integer,pointer     :: nAtomsOfTypeA(:),nPrimOfTypeA(:)      !intent(inout)
integer,pointer     :: nContOfTypeA(:),AngmomOfTypeA(:)      !intent(inout)
integer,pointer     :: startOrbitalOfTypeA(:,:)              !intent(inout)
real(realk),pointer :: exponentsOfTypeA(:,:),Acenters(:,:,:)
real(realk),pointer :: ContractCoeffOfTypeA(:,:,:)
!
#ifdef VAR_ICHOR
TYPE(BASISSETINFO),pointer :: AObasis

IF (intSpecA.EQ.'R') THEN
   !   The regular AO-basis
   AObasis => setting%basis(Center)%p%BINFO(RegBasParam)
ELSEIF(intspecA.EQ.'C')THEN
   !   The CABS AO-type basis
   AObasis => setting%basis(Center)%p%BINFO(CABBasParam)   
ELSE
   call lsquit('Unknown specification in MAIN_ICHORERI_DRIVER',-1)
ENDIF
call Determine_nTypesForBatch(AObasis,ntypesA)
call Determine_nBatches(setting%MOLECULE(Center)%p,AObasis,nBatchesA)
call mem_alloc(nAtomsOfTypeA,ntypesA)
call mem_alloc(AngmomOfTypeA,ntypesA)
call mem_alloc(nPrimOfTypeA,ntypesA)
call mem_alloc(nContOfTypeA,ntypesA)

IF(FullBatch)THEN
   nbatchAstart2 = 1
   nbatchAend2 = nBatchesA
ELSE
   nbatchAstart2 = nbatchAstart
   nbatchAend2 = nbatchAend
ENDIF
call build_TypeInfo1(setting%MOLECULE(Center)%p,AObasis,nTypesA,&
     & nAtomsOfTypeA,AngmomOfTypeA,nPrimOfTypeA,nContOfTypeA,&
     & MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchAstart2,nbatchAend2)

call mem_alloc(startOrbitalOfTypeA,MaxNatomsA,ntypesA)
call mem_alloc(exponentsOfTypeA,MaxnPrimA,ntypesA)
call mem_alloc(ContractCoeffOfTypeA,MaxnPrimA,MaxnContA,ntypesA)
call mem_alloc(Acenters,3,MaxnAtomsA,nTypesA)

call build_TypeInfo2(setting%MOLECULE(Center)%p,AObasis,nTypesA,&
     & spherical,MaxnAtomsA,MaxnPrimA,MaxnContA,startOrbitalOfTypeA,&
     & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,nbatchAstart2,nbatchAend2)

#endif
END Subroutine BUILDCENTERANDTYPEINFO

subroutine FreeCenterAndTypeInfo(nAtomsOfTypeA,AngmomOfTypeA,&
           & nPrimOfTypeA,nContOfTypeA,startOrbitalOfTypeA,&
           & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters)
implicit none
integer,pointer     :: nAtomsOfTypeA(:),nPrimOfTypeA(:)      !intent(inout)
integer,pointer     :: nContOfTypeA(:),AngmomOfTypeA(:)      !intent(inout)
integer,pointer     :: startOrbitalOfTypeA(:,:)              !intent(inout)
real(realk),pointer :: exponentsOfTypeA(:,:),Acenters(:,:,:)
real(realk),pointer :: ContractCoeffOfTypeA(:,:,:)
call mem_dealloc(nAtomsOfTypeA)
call mem_dealloc(AngmomOfTypeA)
call mem_dealloc(nPrimOfTypeA)
call mem_dealloc(nContOfTypeA)
call mem_dealloc(startOrbitalOfTypeA)
call mem_dealloc(exponentsOfTypeA)
call mem_dealloc(ContractCoeffOfTypeA)
call mem_dealloc(Acenters)
end subroutine FreeCenterAndTypeInfo


SUBROUTINE SCREEN_ICHORERI_DRIVER(LUPRI,IPRINT,setting,intspec,SameMOL)
implicit none
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: LUPRI,IPRINT
Character,intent(IN)      :: intSpec(5)
logical,intent(IN) :: SameMOL
!
#ifdef VAR_ICHOR
integer                :: nTypes
!A
integer :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchesA
integer,pointer     :: nAtomsOfTypeA(:),nPrimOfTypeA(:)
integer,pointer     :: nContOfTypeA(:),AngmomOfTypeA(:),startOrbitalOfTypeA(:,:)
real(realk),pointer :: exponentsOfTypeA(:,:),Acenters(:,:,:),ContractCoeffOfTypeA(:,:,:)
!B
integer :: nTypesB,MaxnAtomsB,MaxnPrimB,MaxnContB,nbatchesB
integer,pointer     :: nAtomsOfTypeB(:),nPrimOfTypeB(:)
integer,pointer     :: nContOfTypeB(:),AngmomOfTypeB(:),startOrbitalOfTypeB(:,:)
real(realk),pointer :: exponentsOfTypeB(:,:),Bcenters(:,:,:),ContractCoeffOfTypeB(:,:,:)
!C
integer :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC,nbatchesC
integer,pointer     :: nAtomsOfTypeC(:),nPrimOfTypeC(:)
integer,pointer     :: nContOfTypeC(:),AngmomOfTypeC(:),startOrbitalOfTypeC(:,:)
real(realk),pointer :: exponentsOfTypeC(:,:),Ccenters(:,:,:),ContractCoeffOfTypeC(:,:,:)
!D
integer :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD,nbatchesD
integer,pointer     :: nAtomsOfTypeD(:),nPrimOfTypeD(:)
integer,pointer     :: nContOfTypeD(:),AngmomOfTypeD(:),startOrbitalOfTypeD(:,:)
real(realk),pointer :: exponentsOfTypeD(:,:),Dcenters(:,:,:),ContractCoeffOfTypeD(:,:,:)
!job specification
integer :: SphericalSpec,IchorJobSpec,IchorInputSpec,IchorParSpec,IchorScreenSpec
Integer :: IchorInputDim1,IchorInputDim2,IchorInputDim3,IchorDebugSpec,IchorAlgoSpec
integer :: filestorageIdentifier,MaxFileStorage,IchorPermuteSpec
Integer :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5
Integer :: GabIdentifier, IchorGabID1, IchorGabID2,IchorOperatorSpec
logical :: SameLHSaos
real(realk) :: THRESHOLD_CS,THRESHOLD_QQR
real(realk),pointer :: InputStorage(:)
real(realk),pointer :: BATCHGAB(:),BATCHGCD(:)
Integer(kind=long) :: MaxMem,MaxMemAllocated,MemAllocated
!Integer :: nBatchesA,nBatchesB,nBatchesC,nBatchesD
logical   :: spherical,doLink,rhsDmat,FullBatch
TYPE(BASISSETINFO),pointer :: AObasis
integer :: nbatchAstart2,nbatchAend2,nbatchBstart2,nbatchBend2
integer :: nbatchCstart2,nbatchCend2,nbatchDstart2,nbatchDend2
call FreeIchorSaveGabModuleInterface
call InitIchorSaveGabModuleInterface

FullBatch = .TRUE.

spherical = .TRUE.
!IF (intSpec(5).NE.'C') CALL LSQUIT('MAIN_ICHORERI_DRIVER limited to Coulomb Integrals for now',-1)

IF(SETTING%SCHEME%CS_SCREEN)THEN
   !A
   Call BuildCenterAndTypeInfo(1,intSpec(1),setting,ntypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
        & nPrimOfTypeA,nContOfTypeA,nbatchAstart2,nbatchAend2,MaxnAtomsA,MaxnPrimA,MaxnContA,&
        & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,FullBatch,&
        & 1,nBatchesA,spherical)
   !B
   Call BuildCenterAndTypeInfo(2,intSpec(2),setting,ntypesB,nBatchesB,nAtomsOfTypeB,AngmomOfTypeB,&
        & nPrimOfTypeB,nContOfTypeB,nbatchBstart2,nbatchBend2,MaxnAtomsB,MaxnPrimB,MaxnContB,&
        & startOrbitalOfTypeB,exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters,FullBatch,&
        & 1,nBatchesB,spherical)
   
   call GetIchorSphericalParamIdentifier(SphericalSpec)
   doLink = .FALSE.
   call GetIchorJobEriIdentifier(IchorJobSpec,doLink)
   rhsDmat = .FALSE. !no rhs density matrix supplied as input 
   call GetIchorInputIdentifier(IchorInputSpec,rhsDmat)
   call GetIchorOpereratorIntSpec(intSpec(5),IchorOperatorSpec)
   IchorInputDim1=1                 !not used since   IcorInputNoInput
   IchorInputDim2=1                 !not used since   IcorInputNoInput
   IchorInputDim3=1                 !not used since   IcorInputNoInput
   call mem_alloc(InputStorage,1)
   call GetIchorParallelSpecIdentifier(IchorParSpec)   !no parallelization
   call GetIchorDebugIdentifier(IchorDebugSpec,iprint) !Debug PrintLevel
   call GetIchorAlgorithmSpecIdentifier(IchorAlgoSpec)
   call SetIchorGPUMaxMem(setting%GPUMAXMEM)
   !no Permutational Sym: SameLHSaos=SameRHSaos=SameODs=.FALSE. 
   call GetIchorPermuteParameter(IchorPermuteSpec,.FALSE.,.FALSE.,.FALSE.)
   call GetIchorFileStorageIdentifier(filestorageIdentifier)
   MaxMem = 0                       !Maximum Memory Ichor is allowed to use. Zero means no restrictions
   MaxFileStorage = 0               !Maximum File size, if zero - no file will be written or read. 
   MaxMemAllocated = 0              !Maximum Memory used in the program. Ichor adds to this value
   MemAllocated = 0                 !Memory allocated in the Ichor program
   OutputDim1 = nBatchesA
   OutputDim2 = nBatchesB
   OutputDim3 = 1
   OutputDim4 = 1
   OutputDim5 = 1
   call GetIchorScreeningParameter(IchorScreenSpec,.TRUE.,.TRUE.,.FALSE.)
   !LHS
   IF (intSpec(1).EQ.intSpec(2).AND.SAMEMOL)THEN
      SameLHSaos = .TRUE.
   ELSE
      SameLHSaos = .FALSE.
   ENDIF
   call mem_alloc(BATCHGAB,nBatchesA*nBatchesB)
   call IchorGabInterface(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
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
        & BATCHGAB,lupri)
   call Mem_Add_external_memory(MaxMemAllocated)

   CALL GenerateIdentifier(INTSPEC,GabIdentifier)
   call AddGabToIchorSaveGabModuleInterface(nBatchesA,nBatchesB,&
        & GabIdentifier,BATCHGAB)

   IF(IPRINT.GT.2)THEN
      WRITE(lupri,*)'The LHS Ichor GAB Matrix with Identifier:',GabIdentifier
      call LS_Output(BATCHGAB,1,nBatchesA,1,nBatchesB,nBatchesA,nBatchesB,1,lupri)
   ENDIF
   call mem_dealloc(BATCHGAB)   
   IchorGabID1=GabIdentifier !screening Matrix Identifier
   
   !A
   call FreeCenterAndTypeInfo(nAtomsOfTypeA,AngmomOfTypeA,&
        & nPrimOfTypeA,nContOfTypeA,startOrbitalOfTypeA,&
        & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters)
   call FreeCenterAndTypeInfo(nAtomsOfTypeB,AngmomOfTypeB,&
        & nPrimOfTypeB,nContOfTypeB,startOrbitalOfTypeB,&
        & exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters)

   IF(SameMOL.AND.(intSpec(1).EQ.intspec(3).AND.intSpec(2).EQ.intspec(4)))THEN
      IchorGabID2=IchorGabID1
   ELSE
      !calc RHS
      
      !C
      Call BuildCenterAndTypeInfo(3,intSpec(3),setting,ntypesC,nBatchesC,nAtomsOfTypeC,AngmomOfTypeC,&
           & nPrimOfTypeC,nContOfTypeC,nbatchCstart2,nbatchCend2,MaxnAtomsC,MaxnPrimC,MaxnContC,&
           & startOrbitalOfTypeC,exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters,FullBatch,&
           & 1,nBatchesC,spherical)
      !D
      Call BuildCenterAndTypeInfo(4,intSpec(4),setting,ntypesD,nBatchesD,nAtomsOfTypeD,AngmomOfTypeD,&
           & nPrimOfTypeD,nContOfTypeD,nbatchDstart2,nbatchDend2,MaxnAtomsD,MaxnPrimD,MaxnContD,&
           & startOrbitalOfTypeD,exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters,FullBatch,&
           & 1,nBatchesD,spherical)

      !RHS
      OutputDim1 = nBatchesC
      OutputDim2 = nBatchesD
      OutputDim3 = 1
      OutputDim4 = 1
      OutputDim5 = 1
      IF (intSpec(3).EQ.intSpec(4).AND.SAMEMOL)THEN
         SameLHSaos = .TRUE.
      ELSE
         SameLHSaos = .FALSE.
      ENDIF
      call mem_alloc(BATCHGCD,nBatchesC*nBatchesD)
      call IchorGabInterface(nTypesC,MaxNatomsC,MaxnPrimC,MaxnContC,&
           & AngmomOfTypeC,nAtomsOfTypeC,nPrimOfTypeC,nContOfTypeC,&
           & startOrbitalOfTypeC,Ccenters,exponentsOfTypeC,ContractCoeffOfTypeC,&
           & nTypesD,MaxNatomsD,MaxnPrimD,MaxnContD,&
           & AngmomOfTypeD,nAtomsOfTypeD,nPrimOfTypeD,nContOfTypeD,&
           & startOrbitalOfTypeD,Dcenters,exponentsOfTypeD,ContractCoeffOfTypeD,&
           & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorOperatorSpec,&
           & IchorInputDim1,IchorInputDim2,IchorInputDim3,&
           & InputStorage,IchorParSpec,IchorScreenSpec,IchorDebugSpec,&
           & IchorAlgoSpec,SameLHSaos,filestorageIdentifier,MaxMem,&
           & MaxFileStorage,MaxMemAllocated,MemAllocated,&
           & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
           & BATCHGCD,lupri)
      call mem_Add_external_memory(MaxMemAllocated)

      CALL GenerateIdentifier(INTSPEC,GabIdentifier)
      IF(GabIdentifier.EQ.IchorGabID1)THEN
         GabIdentifier = GabIdentifier + 53210
      ENDIF
      call AddGabToIchorSaveGabModuleInterface(nBatchesC,nBatchesD,GabIdentifier,BATCHGCD)

      IF(IPRINT.GT.2)THEN
         WRITE(lupri,*)'The RHS Ichor GAB Matrix with Identifier:',GabIdentifier
         call LS_Output(BATCHGCD,1,nBatchesC,1,nBatchesD,nBatchesC,nBatchesD,1,lupri)
      ENDIF
      call mem_dealloc(BATCHGCD)
      IchorGabID2=GabIdentifier !screening Matrix Identifier   

      !free space
      call FreeCenterAndTypeInfo(nAtomsOfTypeC,AngmomOfTypeC,&
           & nPrimOfTypeC,nContOfTypeC,startOrbitalOfTypeC,&
           & exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters)
      call FreeCenterAndTypeInfo(nAtomsOfTypeD,AngmomOfTypeD,&
           & nPrimOfTypeD,nContOfTypeD,startOrbitalOfTypeD,&
           & exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters)
   ENDIF
   call GetIchorScreeningParameter(IchorScreenSpec,.TRUE.,.FALSE.,.FALSE.)
   call mem_dealloc(InputStorage)
ELSE   
   call GetIchorScreeningParameter(IchorScreenSpec,.FALSE.,.FALSE.,.FALSE.)
   IchorGabID1=0 !screening Matrix Identifier, not used if IchorScreenNone
   IchorGabID2=0 !screening Matrix Identifier, not used if IchorScreenNone
ENDIF
CALL SET_IchorGabIDinterface(IchorGabID1,IchorGabID2)

#endif
END SUBROUTINE SCREEN_ICHORERI_DRIVER

!K_(A,C) = (AB|CD) D_(BD)
SUBROUTINE MAIN_LINK_ICHORERI_DRIVER(LUPRI,IPRINT,setting,dim1,dim2,dim3,dim4,nDmat,Kmat,Dmat,intspec,FullBatch,&
     & nbatchAstart,nbatchAend,nbatchBstart,nbatchBend,nbatchCstart,nbatchCend,nbatchDstart,nbatchDend,intThreshold)
implicit none
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: LUPRI,IPRINT,dim1,dim2,dim3,dim4,nDmat
logical,intent(in)        :: FullBatch!full batches are assumed and the nbatchXXX arguments are not used 
integer,intent(in)        :: nbatchAstart,nbatchAend,nbatchBstart,nbatchBend
integer,intent(in)        :: nbatchCstart,nbatchCend,nbatchDstart,nbatchDend
real(realk),intent(in)    :: Dmat(dim2,dim4,nDmat),intThreshold
real(realk),intent(inout) :: Kmat(dim1,dim3,nDmat)
Character,intent(IN)      :: intSpec(5)
!
#ifdef VAR_ICHOR
integer                :: nTypes
!A
integer :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchA
integer,pointer     :: nAtomsOfTypeA(:),nPrimOfTypeA(:)
integer,pointer     :: nContOfTypeA(:),AngmomOfTypeA(:),startOrbitalOfTypeA(:,:)
real(realk),pointer :: exponentsOfTypeA(:,:),Acenters(:,:,:),ContractCoeffOfTypeA(:,:,:)
!B
integer :: nTypesB,MaxnAtomsB,MaxnPrimB,MaxnContB,nbatchB
integer,pointer     :: nAtomsOfTypeB(:),nPrimOfTypeB(:)
integer,pointer     :: nContOfTypeB(:),AngmomOfTypeB(:),startOrbitalOfTypeB(:,:)
real(realk),pointer :: exponentsOfTypeB(:,:),Bcenters(:,:,:),ContractCoeffOfTypeB(:,:,:)
!C
integer :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC,nbatchC
integer,pointer     :: nAtomsOfTypeC(:),nPrimOfTypeC(:)
integer,pointer     :: nContOfTypeC(:),AngmomOfTypeC(:),startOrbitalOfTypeC(:,:)
real(realk),pointer :: exponentsOfTypeC(:,:),Ccenters(:,:,:),ContractCoeffOfTypeC(:,:,:)
!D
integer :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD,nbatchD
integer,pointer     :: nAtomsOfTypeD(:),nPrimOfTypeD(:)
integer,pointer     :: nContOfTypeD(:),AngmomOfTypeD(:),startOrbitalOfTypeD(:,:)
real(realk),pointer :: exponentsOfTypeD(:,:),Dcenters(:,:,:),ContractCoeffOfTypeD(:,:,:)
!job specification
integer :: SphericalSpec,IchorJobSpec,IchorInputSpec,IchorParSpec,IchorScreenSpec
Integer :: IchorInputDim1,IchorInputDim2,IchorInputDim3,IchorDebugSpec,IchorAlgoSpec
integer :: filestorageIdentifier,MaxFileStorage,IchorPermuteSpec
Integer :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5
Integer :: GabIdentifier, IchorGabID1, IchorGabID2
logical :: SameLHSaos
real(realk) :: THRESHOLD_CS,THRESHOLD_QQR,THRESHOLD_OD
!real(realk),pointer :: InputStorage(:)
real(realk),pointer :: BATCHGAB(:),BATCHGCD(:)
Integer(kind=long) :: MaxMem,MaxMemAllocated,MemAllocated
Integer :: nBatchesA,nBatchesB,nBatchesC,nBatchesD,IchorOperatorSpec
logical :: spherical
TYPE(BASISSETINFO),pointer :: AObasis
integer :: nbatchAstart2,nbatchAend2,nbatchBstart2,nbatchBend2
integer :: nbatchCstart2,nbatchCend2,nbatchDstart2,nbatchDend2
logical :: SameRHSaos,SameODs,CRIT1,CRIT2,CRIT3,CRIT4,doLink,rhsDmat
logical :: ForceCPU,ForceGPU

spherical = .TRUE.
!IF (intSpec(5).NE.'C') CALL LSQUIT('MAIN_LINK_ICHORERI_DRIVER limited to Coulomb Integrals for now',-1)
!A
Call BuildCenterAndTypeInfo(1,intSpec(1),setting,ntypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,nbatchAstart2,nbatchAend2,MaxnAtomsA,MaxnPrimA,MaxnContA,&
     & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,FullBatch,&
     & nbatchAstart,nbatchAend,spherical)
!B
Call BuildCenterAndTypeInfo(2,intSpec(2),setting,ntypesB,nBatchesB,nAtomsOfTypeB,AngmomOfTypeB,&
     & nPrimOfTypeB,nContOfTypeB,nbatchBstart2,nbatchBend2,MaxnAtomsB,MaxnPrimB,MaxnContB,&
     & startOrbitalOfTypeB,exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters,FullBatch,&
     & nbatchBstart,nbatchBend,spherical)
!C
Call BuildCenterAndTypeInfo(3,intSpec(3),setting,ntypesC,nBatchesC,nAtomsOfTypeC,AngmomOfTypeC,&
     & nPrimOfTypeC,nContOfTypeC,nbatchCstart2,nbatchCend2,MaxnAtomsC,MaxnPrimC,MaxnContC,&
     & startOrbitalOfTypeC,exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters,FullBatch,&
     & nbatchCstart,nbatchCend,spherical)
!D
Call BuildCenterAndTypeInfo(4,intSpec(4),setting,ntypesD,nBatchesD,nAtomsOfTypeD,AngmomOfTypeD,&
     & nPrimOfTypeD,nContOfTypeD,nbatchDstart2,nbatchDend2,MaxnAtomsD,MaxnPrimD,MaxnContD,&
     & startOrbitalOfTypeD,exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters,FullBatch,&
     & nbatchDstart,nbatchDend,spherical)

call GetIchorSphericalParamIdentifier(SphericalSpec)
doLink = .TRUE.
call GetIchorJobEriIdentifier(IchorJobSpec,DoLink)
rhsDmat = .TRUE.
call GetIchorInputIdentifier(IchorInputSpec,rhsDmat)
IchorInputDim1=dim2
IchorInputDim2=dim4        
IchorInputDim3=nDmat       
call GetIchorParallelSpecIdentifier(IchorParSpec)   !no parallelization
call GetIchorDebugIdentifier(IchorDebugSpec,iprint) !Debug PrintLevel
call GetIchorAlgorithmSpecIdentifier(IchorAlgoSpec)
call SetIchorGPUMaxMem(setting%GPUMAXMEM)

SameLHSaos = .FALSE.!intSpec(1).EQ.intSpec(2).AND.Setting%sameMol(1,2)
SameRHSaos = .FALSE.!intSpec(3).EQ.intSpec(4).AND.Setting%sameMol(3,4)
SameODs = .FALSE.!(intSpec(1).EQ.intSpec(3)).AND.(intSpec(2).EQ.intSpec(4)).AND.(Setting%sameMol(1,3).AND.Setting%sameMol(3,4))

!!$IF(FullBatch)THEN
!!$   SameLHSaos = intSpec(1).EQ.intSpec(2).AND.Setting%sameMol(1,2)
!!$   SameRHSaos = intSpec(3).EQ.intSpec(4).AND.Setting%sameMol(3,4)
!!$   SameODs = (intSpec(1).EQ.intSpec(3)).AND.(intSpec(2).EQ.intSpec(4)).AND.(Setting%sameMol(1,3).AND.Setting%sameMol(3,4))
!!$ELSE
!!$   SameLHSaos = (intSpec(1).EQ.intSpec(2).AND.Setting%sameMol(1,2)).AND.&
!!$        & ((nbatchAstart2.EQ.nbatchBstart2).AND.(nbatchAend2.EQ.nbatchBend2))
!!$   SameRHSaos = (intSpec(3).EQ.intSpec(4).AND.Setting%sameMol(3,4)).AND.&
!!$        & ((nbatchCstart2.EQ.nbatchDstart2).AND.(nbatchCend2.EQ.nbatchDend2))
!!$   CRIT1 = (intSpec(1).EQ.intSpec(3)).AND.(intSpec(2).EQ.intSpec(4))
!!$   CRIT2 = Setting%sameMol(1,3).AND.Setting%sameMol(2,4)
!!$   CRIT3 = (nbatchAstart2.EQ.nbatchCstart2).AND.(nbatchAend2.EQ.nbatchCend2)
!!$   CRIT4 = (nbatchBstart2.EQ.nbatchDstart2).AND.(nbatchBend2.EQ.nbatchDend2)
!!$   SameODs = (CRIT1.AND.CRIT2).AND.(CRIT3.AND.CRIT4)
!!$ENDIF

call GetIchorPermuteParameter(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
call GetIchorFileStorageIdentifier(filestorageIdentifier)
call GetIchorOpereratorIntSpec(intSpec(5),IchorOperatorSpec)

MaxMem=0         !Maximum Memory Ichor is allowed to use. Zero = no restrictions
MaxFileStorage=0 !Maximum File size, if zero - no file will be written or read. 
MaxMemAllocated=0!Maximum Memory used in the program. Ichor adds to this value
MemAllocated = 0 !Memory allocated in the Ichor program

call GetIchorScreeningParameter(IchorScreenSpec,SETTING%SCHEME%CS_SCREEN,&
     & SETTING%SCHEME%OD_SCREEN,.FALSE.)
IF(SETTING%SCHEME%CS_SCREEN.AND.SETTING%SCHEME%OD_SCREEN)THEN
   CALL GET_IchorGabIDInterface(IchorGabID1,IchorGabID2)
ELSEIF(SETTING%SCHEME%CS_SCREEN)THEN
   CALL GET_IchorGabIDInterface(IchorGabID1,IchorGabID2)
ELSE   
   IchorGabID1=0 !screening Matrix Identifier, not used if IchorScreenNone
   IchorGabID2=0 !screening Matrix Identifier, not used if IchorScreenNone
ENDIF

OutputDim1=dim1
OutputDim2=dim3
OutputDim3=1
OutputDim4=1
OutputDim5=nDmat
THRESHOLD_OD = intThreshold*1.0E-1_realk
THRESHOLD_CS = intThreshold
THRESHOLD_QQR = intThreshold
ForceCPU = SETTING%SCHEME%IchorForceCPU
ForceGPU = SETTING%SCHEME%IchorForceGPU
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
     & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorOperatorSpec,&
     & IchorInputDim1,IchorInputDim2,IchorInputDim3,&
     & Dmat,IchorParSpec,IchorScreenSpec,THRESHOLD_OD,THRESHOLD_CS,&
     & THRESHOLD_QQR,IchorGabID1,IchorGabID2,IchorDebugSpec,&
     & IchorAlgoSpec,IchorPermuteSpec,filestorageIdentifier,MaxMem,&
     & MaxFileStorage,MaxMemAllocated,MemAllocated,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
     & Kmat,ForceCPU,ForceGPU,lupri)

call Mem_Add_external_memory(MaxMemAllocated)
!=====================================================================


!=====================================================================
!free space
call FreeCenterAndTypeInfo(nAtomsOfTypeA,AngmomOfTypeA,&
           & nPrimOfTypeA,nContOfTypeA,startOrbitalOfTypeA,&
           & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeB,AngmomOfTypeB,&
           & nPrimOfTypeB,nContOfTypeB,startOrbitalOfTypeB,&
           & exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeC,AngmomOfTypeC,&
           & nPrimOfTypeC,nContOfTypeC,startOrbitalOfTypeC,&
           & exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeD,AngmomOfTypeD,&
           & nPrimOfTypeD,nContOfTypeD,startOrbitalOfTypeD,&
           & exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters)

#endif
END SUBROUTINE MAIN_LINK_ICHORERI_DRIVER

SUBROUTINE SCREEN_ICHORERI_RETRIEVE_GABDIM(LUPRI,IPRINT,setting,nBatchA,nBatchB,LHS)
implicit none
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: LUPRI,IPRINT
logical,intent(IN) :: LHS
integer,intent(inout) :: nBatchA,nBatchB
! local variables
#ifdef VAR_ICHOR
integer :: IchorGabID1,IchorGabID2
CALL GET_IchorGabIDinterface(IchorGabID1,IchorGabID2)
IF(LHS)THEN
   call RetrieveGabDIMFromIchorSaveGabModuleInterface(&
        & nBatchA,nBatchB,IchorGabID1)
ELSE
   call RetrieveGabDIMFromIchorSaveGabModuleInterface(&
        & nBatchA,nBatchB,IchorGabID2)   
ENDIF
#endif
END SUBROUTINE SCREEN_ICHORERI_RETRIEVE_GABDIM

SUBROUTINE SCREEN_ICHORERI_RETRIEVE_GAB(LUPRI,IPRINT,setting,nBatchA,nBatchB,LHS,BATCHGAB)
implicit none
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: LUPRI,IPRINT
logical,intent(IN) :: LHS
integer,intent(IN) :: nBatchA,nBatchB
real(realk),intent(inout) :: BATCHGAB(nBatchA*nBatchB)
! local variables
#ifdef VAR_ICHOR
integer :: IchorGabID1,IchorGabID2
CALL GET_IchorGabIDInterface(IchorGabID1,IchorGabID2)
IF(LHS)THEN
   call RetrieveGabFromIchorSaveGabModuleInterface(nBatchA,nBatchB,IchorGabID1,BATCHGAB)
ELSE
   call RetrieveGabFromIchorSaveGabModuleInterface(nBatchA,nBatchB,IchorGabID2,BATCHGAB)
ENDIF
#endif
END SUBROUTINE SCREEN_ICHORERI_RETRIEVE_GAB

SUBROUTINE FREE_SCREEN_ICHORERI()
implicit none
#ifdef VAR_ICHOR
CALL SET_IchorGabIDInterface(0,0)
call FreeIchorSaveGabModuleInterface
#endif
END SUBROUTINE FREE_SCREEN_ICHORERI

SUBROUTINE GenerateIdentifier(INTSPEC,GabIdentifier)
  implicit none
  Character,intent(IN)      :: intSpec(5)
  integer,intent(inout) :: GabIdentifier
  !
#ifdef VAR_ICHOR
  Integer :: I
  GabIdentifier = 0 
  DO I = 1,4
     IF(intSpec(I).EQ.'R')THEN
        GabIdentifier = GabIdentifier+10**I         
     ELSEIF(intSpec(I).EQ.'C')THEN
        GabIdentifier = GabIdentifier+2*10**I         
     ELSE
        call lsquit('unknown spec in GENERATEIDENTIFIER',-1)
     ENDIF
  ENDDO
  GabIdentifier = 0 
  IF(intSpec(5).EQ.'C')THEN
     GabIdentifier = GabIdentifier+100000         
  ELSEIF(intSpec(5).EQ.'G')THEN
     GabIdentifier = GabIdentifier+200000         
  ELSEIF(intSpec(5).EQ.'F')THEN
     GabIdentifier = GabIdentifier+300000         
  ELSEIF(intSpec(5).EQ.'D')THEN
     GabIdentifier = GabIdentifier+400000         
  ELSEIF(intSpec(5).EQ.'2')THEN
     GabIdentifier = GabIdentifier+500000         
  ELSE
     call lsquit('unknown spec in GENERATEIDENTIFIER',-1)
  ENDIF
#endif
END SUBROUTINE GENERATEIDENTIFIER

Subroutine Determine_nTypesForBatch(BASISINFO,ntypes)
implicit none
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO),intent(in):: BASISINFO
!> the number of different types of batches
INTEGER,intent(inout)        :: ntypes
!
#ifdef VAR_ICHOR
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
#endif
END Subroutine DETERMINE_NTYPESFORBATCH

Subroutine Determine_nBatches(MOLECULE,BASISINFO,nBatches)
implicit none
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO),intent(in) :: MOLECULE
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO),intent(in):: BASISINFO
!> the number of Batches
integer,intent(inout)    :: nBatches
#ifdef VAR_ICHOR
!local variables
integer :: R,I,ICHARGE,TYPE,K,iseg
!build nAtomsOfType
nBatches = 0
R = BASISINFO%Labelindex
DO I=1,MOLECULE%natoms  
 IF(MOLECULE%ATOM(I)%pointcharge)CYCLE !no basis functions on this point charge
 IF(R.EQ. 0)THEN
    ICHARGE = INT(MOLECULE%ATOM(I)%CHARGE)
    type = BASISINFO%CHARGEINDEX(ICHARGE)
 ELSE
    type = MOLECULE%ATOM(I)%IDtype(R)
 ENDIF
 IF(BASISINFO%ATOMTYPE(type)%nAngmom.EQ.0) CYCLE
 DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
  DO iseg = 1,BASISINFO%ATOMTYPE(type)%SHELL(K)%nsegments
   nBatches = nBatches + 1 
  ENDDO
 ENDDO
ENDDO
#endif
END Subroutine DETERMINE_NBATCHES

Subroutine Determine_OrbSizeOfBatches(MOLECULE,BASISINFO,nBatches,OrbSizeOfBatches,Spherical)
implicit none
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO),intent(in) :: MOLECULE
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO),intent(in):: BASISINFO
!> if it is spherical or cartesian GTOs
LOGICAL,intent(in)           :: spherical
!> the number of Batches
integer,intent(in)    :: nBatches
!> the orbital size of Batches
integer,intent(inout) :: OrbSizeOfBatches(nBatches)
#ifdef VAR_ICHOR
!local variables
integer :: R,I,ICHARGE,TYPE,K,iseg,ncol,nOrbComp,iBatches
!build nAtomsOfType
iBatches = 0
R = BASISINFO%Labelindex
DO I=1,MOLECULE%natoms  
 IF(MOLECULE%ATOM(I)%pointcharge)CYCLE !no basis functions on this point charge
 IF(R.EQ. 0)THEN
    ICHARGE = INT(MOLECULE%ATOM(I)%CHARGE)
    type = BASISINFO%CHARGEINDEX(ICHARGE)
 ELSE
    type = MOLECULE%ATOM(I)%IDtype(R)
 ENDIF
 IF(BASISINFO%ATOMTYPE(type)%nAngmom.EQ.0) CYCLE
 DO K=1,BASISINFO%ATOMTYPE(type)%nAngmom
  DO iseg = 1,BASISINFO%ATOMTYPE(type)%SHELL(K)%nsegments
   iBatches = iBatches + 1 

   ncol = BASISINFO%ATOMTYPE(type)%SHELL(K)%segment(iseg)%ncol
   IF (spherical) THEN
      nOrbComp = 2*K-1
   ELSE
      nOrbComp = K*(K+1)/2
   ENDIF
   OrbSizeOfBatches(iBatches) = nOrbComp*ncol
  ENDDO
 ENDDO
ENDDO
#endif
END Subroutine DETERMINE_ORBSIZEOFBATCHES

Subroutine build_TypeInfo1(MOLECULE,BASISINFO,nTypes,nAtomsOfType,&
     & AngmomOfType,nPrimOfType,nContOfType,MaxnAtoms,MaxnPrim,MaxnCont,iBatchStart,iBatchEnd)
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
!> the starting Batch index (normally 1 else not the full set is used)
integer,intent(in)      :: iBatchStart
!> the end Batch index (normally nBatches else not the full set is used)
integer,intent(in)      :: iBatchEnd
#ifdef VAR_ICHOR
!
INTEGER,pointer          :: MODELTYPES(:),MODELBATCHTYPES(:,:,:)
INTEGER                  :: maxseg,maxang,I,K,L,iBatchType,R,icharge,type,iseg
INTEGER                  :: nBatchType,nrow,ncol,irow,icol,orbitalindex,nOrbComp
INTEGER                  :: iBatches
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
   nPrimOfType(iBatchType) = 0
   nContOfType(iBatchType) = 0
   AngmomOfType(iBatchType) = 0
enddo
!build nAtomsOfType
iBatches = 0
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
   iBatches = iBatches + 1   
   iBatchType = MODELBATCHTYPES(iseg,K,L)
   IF((iBatches.GE.iBatchStart).AND.(iBatches.LE.iBatchEnd))THEN
      AngmomOfType(iBatchType) = K-1
      nAtomsOfType(iBatchType) = nAtomsOfType(iBatchType) + 1
      nrow = BASISINFO%ATOMTYPE(type)%SHELL(K)%segment(iseg)%nrow
      ncol = BASISINFO%ATOMTYPE(type)%SHELL(K)%segment(iseg)%ncol
      nPrimOfType(iBatchType) = nrow
      nContOfType(iBatchType) = ncol
   ENDIF
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

#endif
end Subroutine Build_TypeInfo1

Subroutine Build_TypeInfo2(MOLECULE,BASISINFO,nTypes,spherical,MaxnAtoms,MaxnPrim,&
     & MaxnCont,startOrbitalOfType,exponentsOfType,ContractCoeffOfType,CentersOfType,iBatchStart,iBatchEnd)
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
!> the starting Batch index (normally 1 else not the full set is used)
integer,intent(in)      :: iBatchStart
!> the end Batch index (normally nBatches else not the full set is used)
integer,intent(in)      :: iBatchEnd
#ifdef VAR_ICHOR
!
INTEGER,pointer          :: MODELTYPES(:),MODELBATCHTYPES(:,:,:)
INTEGER                  :: maxseg,maxang,I,K,L,iBatchType,R,icharge,type,iseg
INTEGER                  :: nBatchType,nrow,ncol,irow,icol,orbitalindex,nOrbComp
!tmp
integer   :: nAtomsOfType(nTypes),iBatches
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
iBatches = 0
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
   iBatches = iBatches + 1   
   iBatchType = MODELBATCHTYPES(iseg,K,L)
   IF((iBatches.GE.iBatchStart).AND.(iBatches.LE.iBatchEnd))THEN
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
   ENDIF
  ENDDO
 ENDDO
ENDDO
call mem_dealloc(MODELTYPES)
call mem_dealloc(MODELBATCHTYPES)

#endif
end Subroutine Build_TypeInfo2

SUBROUTINE WRITE_ICHORERI_INFO(LUPRI,IPRINT,setting,dim1,dim2,dim3,dim4,integrals,profile,LUOUTPUT)
implicit none
TYPE(lssetting),intent(in):: setting
integer,intent(in)        :: LUPRI,IPRINT,LUOUTPUT
integer,intent(in)        :: dim1,dim2,dim3,dim4 
real(realk),intent(inout) :: integrals(dim1,dim2,dim3,dim4)
logical,intent(in)        :: profile
#ifdef VAR_ICHOR
!
integer                :: nTypes
!A
integer :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchesA
integer,pointer     :: nAtomsOfTypeA(:),nPrimOfTypeA(:)
integer,pointer     :: nContOfTypeA(:),AngmomOfTypeA(:),startOrbitalOfTypeA(:,:)
real(realk),pointer :: exponentsOfTypeA(:,:),Acenters(:,:,:),ContractCoeffOfTypeA(:,:,:)
!B
integer :: nTypesB,MaxnAtomsB,MaxnPrimB,MaxnContB,nbatchesB
integer,pointer     :: nAtomsOfTypeB(:),nPrimOfTypeB(:)
integer,pointer     :: nContOfTypeB(:),AngmomOfTypeB(:),startOrbitalOfTypeB(:,:)
real(realk),pointer :: exponentsOfTypeB(:,:),Bcenters(:,:,:),ContractCoeffOfTypeB(:,:,:)
!C
integer :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC,nbatchesC
integer,pointer     :: nAtomsOfTypeC(:),nPrimOfTypeC(:)
integer,pointer     :: nContOfTypeC(:),AngmomOfTypeC(:),startOrbitalOfTypeC(:,:)
real(realk),pointer :: exponentsOfTypeC(:,:),Ccenters(:,:,:),ContractCoeffOfTypeC(:,:,:)
!D
integer :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD,nbatchesD
integer,pointer     :: nAtomsOfTypeD(:),nPrimOfTypeD(:)
integer,pointer     :: nContOfTypeD(:),AngmomOfTypeD(:),startOrbitalOfTypeD(:,:)
real(realk),pointer :: exponentsOfTypeD(:,:),Dcenters(:,:,:),ContractCoeffOfTypeD(:,:,:)
logical :: spherical,FullBatch
integer :: nbatchAstart2,nbatchAend2,nbatchAstart,nbatchAend,I,J,K,L
integer :: nbatchBstart2,nbatchBend2,nbatchBstart,nbatchBend
integer :: nbatchCstart2,nbatchCend2,nbatchCstart,nbatchCend
integer :: nbatchDstart2,nbatchDend2,nbatchDstart,nbatchDend
FullBatch = .TRUE.
nbatchAstart2=1; nbatchAend2=1; nbatchBstart2=1; nbatchBend2=1
nbatchCstart2=1; nbatchCend2=1; nbatchDstart2=1; nbatchDend2=1
spherical = .TRUE.
!A
Call BuildCenterAndTypeInfo(1,'R',setting,ntypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,nbatchAstart2,nbatchAend2,MaxnAtomsA,MaxnPrimA,MaxnContA,&
     & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,Acenters,FullBatch,&
     & nbatchAstart,nbatchAend,spherical)
Call WriteCenterInfo1(luoutput,ntypesA,nBatchesA,MaxnAtomsA,MaxnPrimA,MaxnContA,spherical)
Call WriteCenterInfo2(luoutput,ntypesA,nBatchesA,nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,MaxnAtomsA,MaxnPrimA,MaxnContA,&
     & startOrbitalOfTypeA,exponentsOfTypeA,ContractCoeffOfTypeA,&
     & Acenters,spherical)
call FreeCenterAndTypeInfo(nAtomsOfTypeA,AngmomOfTypeA,&
     & nPrimOfTypeA,nContOfTypeA,startOrbitalOfTypeA,&
     & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters)

!B
Call BuildCenterAndTypeInfo(2,'R',setting,ntypesB,nBatchesB,nAtomsOfTypeB,AngmomOfTypeB,&
     & nPrimOfTypeB,nContOfTypeB,nbatchBstart2,nbatchBend2,MaxnAtomsB,MaxnPrimB,MaxnContB,&
     & startOrbitalOfTypeB,exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters,FullBatch,&
     & nbatchBstart,nbatchBend,spherical)
Call WriteCenterInfo1(luoutput,ntypesB,nBatchesB,MaxnAtomsB,MaxnPrimB,MaxnContB,spherical)
Call WriteCenterInfo2(luoutput,ntypesB,nBatchesB,nAtomsOfTypeB,AngmomOfTypeB,&
     & nPrimOfTypeB,nContOfTypeB,MaxnAtomsB,MaxnPrimB,MaxnContB,&
     & startOrbitalOfTypeB,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & Bcenters,spherical)
call FreeCenterAndTypeInfo(nAtomsOfTypeB,AngmomOfTypeB,&
     & nPrimOfTypeB,nContOfTypeB,startOrbitalOfTypeB,&
     & exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters)

!C
Call BuildCenterAndTypeInfo(3,'R',setting,ntypesC,nBatchesC,nAtomsOfTypeC,AngmomOfTypeC,&
     & nPrimOfTypeC,nContOfTypeC,nbatchCstart2,nbatchCend2,MaxnAtomsC,MaxnPrimC,MaxnContC,&
     & startOrbitalOfTypeC,exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters,FullBatch,&
     & nbatchCstart,nbatchCend,spherical)
Call WriteCenterInfo1(luoutput,ntypesC,nBatchesC,MaxnAtomsC,MaxnPrimC,MaxnContC,spherical)
Call WriteCenterInfo2(luoutput,ntypesC,nBatchesC,nAtomsOfTypeC,AngmomOfTypeC,&
     & nPrimOfTypeC,nContOfTypeC,MaxnAtomsC,MaxnPrimC,MaxnContC,&
     & startOrbitalOfTypeC,exponentsOfTypeC,ContractCoeffOfTypeC,&
     & Ccenters,spherical)
call FreeCenterAndTypeInfo(nAtomsOfTypeC,AngmomOfTypeC,&
     & nPrimOfTypeC,nContOfTypeC,startOrbitalOfTypeC,&
     & exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters)
!D
Call BuildCenterAndTypeInfo(4,'R',setting,ntypesD,nBatchesD,nAtomsOfTypeD,AngmomOfTypeD,&
     & nPrimOfTypeD,nContOfTypeD,nbatchDstart2,nbatchDend2,MaxnAtomsD,MaxnPrimD,MaxnContD,&
     & startOrbitalOfTypeD,exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters,FullBatch,&
     & nbatchDstart,nbatchDend,spherical)
Call WriteCenterInfo1(luoutput,ntypesD,nBatchesD,MaxnAtomsD,MaxnPrimD,MaxnContD,spherical)
Call WriteCenterInfo2(luoutput,ntypesD,nBatchesD,nAtomsOfTypeD,AngmomOfTypeD,&
     & nPrimOfTypeD,nContOfTypeD,MaxnAtomsD,MaxnPrimD,MaxnContD,&
     & startOrbitalOfTypeD,exponentsOfTypeD,ContractCoeffOfTypeD,&
     & Dcenters,spherical)
call FreeCenterAndTypeInfo(nAtomsOfTypeD,AngmomOfTypeD,&
     & nPrimOfTypeD,nContOfTypeD,startOrbitalOfTypeD,&
     & exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters)

IF(.NOT.profile)THEN
   WRITE(LUOUTPUT,*) dim1,dim2,dim3,dim4 
   WRITE(LUOUTPUT,*) integrals
ELSE
   !dimensions huge so we only write a few elements.
   WRITE(LUOUTPUT,*) dim1,dim2,dim3,dim4 
   IF(MIN(dim1,dim2,dim3,dim4).LT.7)THEN
    DO L=1,7
     DO K=1,7
      DO J=1,7
       DO I=1,7
        WRITE(LUOUTPUT,*) integrals(I,J,K,L)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDIF
   DO L=1,dim4,10
    DO K=1,dim3,10
     DO J=1,dim2,10
      DO I=1,dim1,10
       WRITE(LUOUTPUT,*) integrals(I,J,K,L)
      ENDDO
     ENDDO
    ENDDO
   ENDDO
ENDIF
#endif
END SUBROUTINE WRITE_ICHORERI_INFO

SUBROUTINE MAIN_ICHORERI_READDRIVER(LUPRI,IPRINT,LUOUTPUT,CS_SCREEN,OD_SCREEN,integrals,&
           & dim1,dim2,dim3,dim4)
implicit none
integer,intent(in)        :: LUPRI,IPRINT,LUOUTPUT
logical,intent(in)        :: CS_SCREEN,OD_SCREEN
real(realk),pointer       :: integrals(:,:,:,:)
integer,intent(inout)     :: dim1,dim2,dim3,dim4
#ifdef VAR_ICHOR
!
integer                :: nTypes
!A
integer :: nTypesA,MaxnAtomsA,MaxnPrimA,MaxnContA,nbatchA
integer,pointer     :: nAtomsOfTypeA(:),nPrimOfTypeA(:)
integer,pointer     :: nContOfTypeA(:),AngmomOfTypeA(:),startOrbitalOfTypeA(:,:)
real(realk),pointer :: exponentsOfTypeA(:,:),Acenters(:,:,:),ContractCoeffOfTypeA(:,:,:)
!B
integer :: nTypesB,MaxnAtomsB,MaxnPrimB,MaxnContB,nbatchB
integer,pointer     :: nAtomsOfTypeB(:),nPrimOfTypeB(:)
integer,pointer     :: nContOfTypeB(:),AngmomOfTypeB(:),startOrbitalOfTypeB(:,:)
real(realk),pointer :: exponentsOfTypeB(:,:),Bcenters(:,:,:),ContractCoeffOfTypeB(:,:,:)
!C
integer :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC,nbatchC
integer,pointer     :: nAtomsOfTypeC(:),nPrimOfTypeC(:)
integer,pointer     :: nContOfTypeC(:),AngmomOfTypeC(:),startOrbitalOfTypeC(:,:)
real(realk),pointer :: exponentsOfTypeC(:,:),Ccenters(:,:,:),ContractCoeffOfTypeC(:,:,:)
!D
integer :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD,nbatchD
integer,pointer     :: nAtomsOfTypeD(:),nPrimOfTypeD(:)
integer,pointer     :: nContOfTypeD(:),AngmomOfTypeD(:),startOrbitalOfTypeD(:,:)
real(realk),pointer :: exponentsOfTypeD(:,:),Dcenters(:,:,:),ContractCoeffOfTypeD(:,:,:)
!job specification
integer :: SphericalSpec,IchorJobSpec,IchorInputSpec,IchorParSpec,IchorScreenSpec
Integer :: IchorInputDim1,IchorInputDim2,IchorInputDim3,IchorDebugSpec,IchorAlgoSpec
integer :: filestorageIdentifier,MaxFileStorage,IchorPermuteSpec
Integer :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5
Integer :: GabIdentifier, IchorGabID1, IchorGabID2
logical :: SameLHSaos
real(realk) :: THRESHOLD_CS,THRESHOLD_QQR,THRESHOLD_OD
real(realk),pointer :: InputStorage(:)
real(realk),pointer :: BATCHGAB(:),BATCHGCD(:)
Integer(kind=long) :: MaxMem,MaxMemAllocated,MemAllocated
Integer :: nBatchesA,nBatchesB,nBatchesC,nBatchesD
logical :: spherical
TYPE(BASISSETINFO),pointer :: AObasis
logical :: SameRHSaos,SameODs,CRIT1,CRIT2,CRIT3,CRIT4,doLink,rhsDmat,CRIT5
integer :: nbatchAstart2,nbatchAend2
integer :: nbatchBstart2,nbatchBend2
integer :: nbatchCstart2,nbatchCend2
integer :: nbatchDstart2,nbatchDend2
logical :: ForceCPU,ForceGPU
nbatchAstart2=1; nbatchAend2=1; nbatchBstart2=1; nbatchBend2=1
nbatchCstart2=1; nbatchCend2=1; nbatchDstart2=1; nbatchDend2=1
Call ReadCenterInfo1(luoutput,ntypesB,nBatchesB,MaxnAtomsB,MaxnPrimB,MaxnContB,spherical)
call mem_alloc(nAtomsOfTypeA,ntypesA)
call mem_alloc(AngmomOfTypeA,ntypesA)
call mem_alloc(nPrimOfTypeA,ntypesA)
call mem_alloc(nContOfTypeA,ntypesA)
call mem_alloc(startOrbitalOfTypeA,MaxNatomsA,ntypesA)
call mem_alloc(exponentsOfTypeA,MaxnPrimA,ntypesA)
call mem_alloc(ContractCoeffOfTypeA,MaxnPrimA,MaxnContA,ntypesA)
call mem_alloc(Acenters,3,MaxnAtomsA,nTypesA)
Call ReadCenterInfo2(luoutput,ntypesB,nBatchesB,nAtomsOfTypeB,AngmomOfTypeB,&
     & nPrimOfTypeB,nContOfTypeB,MaxnAtomsB,MaxnPrimB,MaxnContB,&
     & startOrbitalOfTypeB,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & Bcenters,spherical)

Call ReadCenterInfo1(luoutput,ntypesB,nBatchesB,MaxnAtomsB,MaxnPrimB,MaxnContB,spherical)
call mem_alloc(nAtomsOfTypeB,ntypesB)
call mem_alloc(AngmomOfTypeB,ntypesB)
call mem_alloc(nPrimOfTypeB,ntypesB)
call mem_alloc(nContOfTypeB,ntypesB)
call mem_alloc(startOrbitalOfTypeB,MaxNatomsB,ntypesB)
call mem_alloc(exponentsOfTypeB,MaxnPrimB,ntypesB)
call mem_alloc(ContractCoeffOfTypeB,MaxnPrimB,MaxnContB,ntypesB)
call mem_alloc(Bcenters,3,MaxnAtomsB,nTypesB)
Call ReadCenterInfo2(luoutput,ntypesB,nBatchesB,nAtomsOfTypeB,AngmomOfTypeB,&
     & nPrimOfTypeB,nContOfTypeB,MaxnAtomsB,MaxnPrimB,MaxnContB,&
     & startOrbitalOfTypeB,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & Bcenters,spherical)

!C
Call ReadCenterInfo1(luoutput,ntypesC,nBatchesC,MaxnAtomsC,MaxnPrimC,MaxnContC,spherical)
call mem_alloc(nAtomsOfTypeC,ntypesC)
call mem_alloc(AngmomOfTypeC,ntypesC)
call mem_alloc(nPrimOfTypeC,ntypesC)
call mem_alloc(nContOfTypeC,ntypesC)
call mem_alloc(startOrbitalOfTypeC,MaxNatomsC,ntypesC)
call mem_alloc(exponentsOfTypeC,MaxnPrimC,ntypesC)
call mem_alloc(ContractCoeffOfTypeC,MaxnPrimC,MaxnContC,ntypesC)
call mem_alloc(Ccenters,3,MaxnAtomsC,nTypesC)
Call ReadCenterInfo2(luoutput,ntypesC,nBatchesC,nAtomsOfTypeC,AngmomOfTypeC,&
     & nPrimOfTypeC,nContOfTypeC,MaxnAtomsC,MaxnPrimC,MaxnContC,&
     & startOrbitalOfTypeC,exponentsOfTypeC,ContractCoeffOfTypeC,&
     & Ccenters,spherical)

!D
Call ReadCenterInfo1(luoutput,ntypesD,nBatchesD,MaxnAtomsD,MaxnPrimD,MaxnContD,spherical)
call mem_alloc(nAtomsOfTypeD,ntypesD)
call mem_alloc(AngmomOfTypeD,ntypesD)
call mem_alloc(nPrimOfTypeD,ntypesD)
call mem_alloc(nContOfTypeD,ntypesD)
call mem_alloc(startOrbitalOfTypeD,MaxNatomsD,ntypesD)
call mem_alloc(exponentsOfTypeD,MaxnPrimD,ntypesD)
call mem_alloc(ContractCoeffOfTypeD,MaxnPrimD,MaxnContD,ntypesD)
call mem_alloc(Dcenters,3,MaxnAtomsD,nTypesD)
Call ReadCenterInfo2(luoutput,ntypesD,nBatchesD,nAtomsOfTypeD,AngmomOfTypeD,&
     & nPrimOfTypeD,nContOfTypeD,MaxnAtomsD,MaxnPrimD,MaxnContD,&
     & startOrbitalOfTypeD,exponentsOfTypeD,ContractCoeffOfTypeD,&
     & Dcenters,spherical)


READ(LUOUTPUT,*) dim1,dim2,dim3,dim4 
call mem_alloc(integrals,dim1,dim2,dim3,dim4)


!===================================================================0

call GetIchorSphericalParamIdentifier(SphericalSpec)
doLink = .FALSE.
call GetIchorJobEriIdentifier(IchorJobSpec,doLink)
rhsDmat = .FALSE. !no rhs density matrix supplied as input 
call GetIchorInputIdentifier(IchorInputSpec,rhsDmat)
IchorInputDim1=1                 !not used since   IcorInputNoInput
IchorInputDim2=1                 !not used since   IcorInputNoInput
IchorInputDim3=1                 !not used since   IcorInputNoInput
call mem_alloc(InputStorage,1)
call GetIchorParallelSpecIdentifier(IchorParSpec)   !no parallelization
call GetIchorDebugIdentifier(IchorDebugSpec,iprint) !Debug PrintLevel
call GetIchorAlgorithmSpecIdentifier(IchorAlgoSpec)
call SetIchorGPUMaxMem(0.0E0_realk)
SameLHSaos = .FALSE. 
SameRHSaos = .FALSE. 
SameODs = .FALSE. 
call GetIchorPermuteParameter(IchorPermuteSpec,SameLHSaos,SameRHSaos,SameODs)
call GetIchorFileStorageIdentifier(filestorageIdentifier)
MaxMem=0         !Maximum Memory Ichor is allowed to use. Zero = no restrictions
MaxFileStorage=0 !Maximum File size, if zero - no file will be written or read. 
MaxMemAllocated=0!Maximum Memory used in the program. Ichor adds to this value
MemAllocated = 0 !Memory allocated in the Ichor program
call GetIchorScreeningParameter(IchorScreenSpec,CS_SCREEN,OD_SCREEN,.FALSE.)
IF(CS_SCREEN.AND.OD_SCREEN)THEN
   CALL GET_IchorGabIDInterface(IchorGabID1,IchorGabID2)
ELSEIF(CS_SCREEN)THEN
   CALL GET_IchorGabIDInterface(IchorGabID1,IchorGabID2)
ELSE   
   IchorGabID1=0 !screening Matrix Identifier, not used if IchorScreenNone
   IchorGabID2=0 !screening Matrix Identifier, not used if IchorScreenNone
ENDIF
OutputDim1=dim1
OutputDim2=dim2
OutputDim3=dim3
OutputDim4=dim4
OutputDim5=1
THRESHOLD_OD = 1.0E-8_realk*1.0E-1_realk
THRESHOLD_CS = 1.0E-8_realk*1.0E-2_realk
THRESHOLD_QQR = 1.0E-8_realk*1.0E+0_realk*0.50E0_realk
ForceCPU = .FALSE.!SETTING%SCHEME%IchorForceCPU
ForceGPU = .FALSE.!SETTING%SCHEME%IchorForceGPU
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
     & integrals,ForceCPU,ForceGPU,lupri)

call Mem_Add_external_memory(MaxMemAllocated)
call mem_dealloc(InputStorage)
!=====================================================================


!=====================================================================
!free space
call FreeCenterAndTypeInfo(nAtomsOfTypeA,AngmomOfTypeA,&
           & nPrimOfTypeA,nContOfTypeA,startOrbitalOfTypeA,&
           & exponentsOfTypeA,ContractCoeffOfTypeA,Acenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeB,AngmomOfTypeB,&
           & nPrimOfTypeB,nContOfTypeB,startOrbitalOfTypeB,&
           & exponentsOfTypeB,ContractCoeffOfTypeB,Bcenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeC,AngmomOfTypeC,&
           & nPrimOfTypeC,nContOfTypeC,startOrbitalOfTypeC,&
           & exponentsOfTypeC,ContractCoeffOfTypeC,Ccenters)
call FreeCenterAndTypeInfo(nAtomsOfTypeD,AngmomOfTypeD,&
           & nPrimOfTypeD,nContOfTypeD,startOrbitalOfTypeD,&
           & exponentsOfTypeD,ContractCoeffOfTypeD,Dcenters)
#endif
END SUBROUTINE MAIN_ICHORERI_READDRIVER

END MODULE IchorErimoduleHost
