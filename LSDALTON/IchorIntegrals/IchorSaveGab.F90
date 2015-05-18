MODULE IchorSaveGabMod
  use IchorprecisionMod
  use IchorCommonMod
!  use IchorMemory
  public :: InitIchorSaveGabModule, AddGabToIchorSaveGabModule,&
       & RetrieveGabFromIchorSaveGabModule,FreeIchorSaveGabModule,&
       & SET_IchorGabID,GET_IchorGabID,RetrieveGabDimFromIchorSaveGabModule
  private
  
TYPE GABM
   real(realk),pointer :: BATCHGAB(:) 
   integer :: GabIdentifier
   integer :: nBatch1
   integer :: nBatch2
END type GABM

TYPE GABsaveItem
   integer :: nGABM
   type(GABM),pointer :: GABM(:) 
END type GABsaveItem

integer(KIND=long),save :: mem_allocated_screenichor
integer(KIND=long),parameter :: memrealsz=8
logical,save :: allocated
type(GABsaveItem),save :: GabItem
Integer,save :: IchorGabID1current
Integer,save :: IchorGabID2current
CONTAINS
subroutine InitIchorSaveGabModule()
implicit none
mem_allocated_screenichor = 0 
allocated = .FALSE.
nullify(GABItem%GABM)
end subroutine InitIchorSaveGabModule

subroutine AddGabToIchorSaveGabModule(nBatch1,nBatch2,GabIdentifier,BATCHGAB)
implicit none
integer,intent(in) :: nBatch1,nBatch2,GabIdentifier
real(realk),intent(in) :: BATCHGAB(nBatch1*nBatch2) 
! local variables
integer :: I,J,K,n
real(realk),pointer :: elms(:)
type(GABM),pointer :: GABMTMP(:) 

IF(allocated)THEN
   !copy to tmp
   allocate(GABMTMP(GABItem%nGABM))
   do K = 1,GABItem%nGABM
    allocate(GABMTMP(K)%BATCHGAB(GABItem%GABM(K)%nBatch1*GABItem%GABM(K)%nBatch2))
    mem_allocated_screenichor = mem_allocated_screenichor + &
         & size(GABMTMP(K)%BATCHGAB,KIND=long)*memrealsz
    n = GABItem%GABM(K)%nBatch1*GABItem%GABM(K)%nBatch2
    call gabdcopy(n,GABItem%GABM(K)%BATCHGAB,GABMTMP(K)%BATCHGAB)
    GABMTMP(K)%GabIdentifier = GABItem%GABM(K)%GabIdentifier
    GABMTMP(K)%nBatch1 = GABItem%GABM(K)%nBatch1
    GABMTMP(K)%nBatch2 = GABItem%GABM(K)%nBatch2
    mem_allocated_screenichor = mem_allocated_screenichor - &
         & size(GABItem%GABM(K)%BATCHGAB,KIND=long)*memrealsz
    deallocate(GABItem%GABM(K)%BATCHGAB)
   enddo
   deallocate(GABItem%GABM)
   allocate(GABItem%GABM(GABItem%nGABM+1))
   do K = 1,GABItem%nGABM
    allocate(GABItem%GABM(K)%BATCHGAB(GABMTMP(K)%nBatch1*GABMTMP(K)%nBatch2))
    mem_allocated_screenichor = mem_allocated_screenichor + &
         & size(GABItem%GABM(K)%BATCHGAB,KIND=long)*memrealsz
    n = GABMTMP(K)%nBatch1*GABMTMP(K)%nBatch2
    call gabdcopy(n,GABMTMP(K)%BATCHGAB,GABItem%GABM(K)%BATCHGAB)
    GABItem%GABM(K)%GabIdentifier = GABMTMP(K)%GabIdentifier
    GABItem%GABM(K)%nBatch1 = GABMTMP(K)%nBatch1
    GABItem%GABM(K)%nBatch2 = GABMTMP(K)%nBatch2
    mem_allocated_screenichor = mem_allocated_screenichor - &
         & size(GABMTMP(K)%BATCHGAB,KIND=long)*memrealsz
    deallocate(GABMTMP(K)%BATCHGAB)
   enddo

   K = GABItem%nGABM + 1
   GABItem%nGABM = K
   allocate(GABItem%GABM(K)%BATCHGAB(nBatch1*nBatch2))
   mem_allocated_screenichor = mem_allocated_screenichor + &
        & size(GABItem%GABM(K)%BATCHGAB,KIND=long)*memrealsz
   n = nBatch2*nBatch1
   call gabdcopy(n,BATCHGAB,GABItem%GABM(K)%BATCHGAB)
   GABItem%GABM(K)%GabIdentifier = GabIdentifier
   GABItem%GABM(K)%nBatch1 = nBatch1
   GABItem%GABM(K)%nBatch2 = nBatch2
ELSE
   allocate(GABItem%GABM(1))
   GABItem%nGABM = 1
   allocate(GABItem%GABM(1)%BATCHGAB(nBatch1*nBatch2))
   mem_allocated_screenichor = mem_allocated_screenichor + &
        & size(GABItem%GABM(1)%BATCHGAB,KIND=long)*memrealsz
   n = nBatch2*nBatch1
   call gabdcopy(n,BATCHGAB,GABItem%GABM(1)%BATCHGAB)
   GABItem%GABM(1)%GabIdentifier = GabIdentifier
   GABItem%GABM(1)%nBatch1 = nBatch1
   GABItem%GABM(1)%nBatch2 = nBatch2
   allocated = .TRUE.
ENDIF
end subroutine AddGabToIchorSaveGabModule

subroutine RetrieveGabFromIchorSaveGabModule(nBatch1,nBatch2,GabIdentifier,BATCHGAB)
implicit none
integer,intent(in) :: nBatch1,nBatch2,GabIdentifier
real(realk),intent(inout) :: BATCHGAB(nBatch1*nBatch2) 
! local variables
integer :: I,J,K
real(realk),pointer :: elms(:)
IF(allocated)THEN
   gabllop: DO K=1,GABItem%nGABM
      !Sanity check
      IF(GABItem%GABM(K)%GabIdentifier .EQ. GabIdentifier)THEN
         IF(GABItem%GABM(K)%nBatch1 .NE. nBatch1)THEN
            CALL ICHORQUIT('RetrieveGab: Dim1 wrong ',-1)
         ENDIF
         IF(GABItem%GABM(K)%nBatch2 .NE. nBatch2)THEN
            CALL ICHORQUIT('RetrieveGab: Dim2 wrong ',-1)
         ENDIF
         call gabdcopy(nBatch2*nBatch1,GABItem%GABM(K)%BATCHGAB,BATCHGAB)
         EXIT gabllop
      ENDIF
   ENDDO gabllop
ELSE
   CALL ICHORQUIT('RetrieveGab: GABitem Not allocated ',-1)
ENDIF
end subroutine RetrieveGabFromIchorSaveGabModule

subroutine RetrieveGabDimFromIchorSaveGabModule(nBatch1,nBatch2,GabIdentifier)
implicit none
integer,intent(in) :: GabIdentifier
integer,intent(inout) :: nBatch1,nBatch2
! local variables
integer :: I,J,K
real(realk),pointer :: elms(:)
IF(allocated)THEN
   gabllop: DO K=1,GABItem%nGABM
      !Sanity check
      IF(GABItem%GABM(K)%GabIdentifier .EQ. GabIdentifier)THEN
         nBatch1 = GABItem%GABM(K)%nBatch1
         nBatch2 = GABItem%GABM(K)%nBatch2
         EXIT gabllop
      ENDIF
   ENDDO gabllop
ELSE
   CALL ICHORQUIT('RetrieveGabDim: GABitem Not allocated ',-1)
ENDIF
end subroutine RetrieveGabDimFromIchorSaveGabModule

subroutine FreeIchorSaveGabModule()
implicit none
integer :: I
IF(allocated)THEN
   DO I=1,GABItem%nGABM
      mem_allocated_screenichor = mem_allocated_screenichor - &
           & size(GABItem%GABM(I)%BATCHGAB,KIND=long)*memrealsz
      deallocate(GABItem%GABM(I)%BATCHGAB)
   ENDDO
   deallocate(GABItem%GABM)
ENDIF
nullify(GABItem%GABM)
allocated = .FALSE.
IF(mem_allocated_screenichor.NE.0)THEN
   CALL ICHORQUIT('FreeIchorSaveGabModule: still mem allocated ',-1)
ENDIF
end subroutine FreeIchorSaveGabModule

subroutine SET_IchorGabID(IchorGabID1,IchorGabID2)
implicit none
integer,intent(in) :: IchorGabID1,IchorGabID2
IchorGabID1current = IchorGabID1
IchorGabID2current = IchorGabID2
end subroutine SET_IchorGabID

subroutine GET_IchorGabID(IchorGabID1,IchorGabID2)
implicit none
integer,intent(inout) :: IchorGabID1,IchorGabID2
IchorGabID1 = IchorGabID1current
IchorGabID2 = IchorGabID2current
end subroutine GET_IchorGabID

subroutine gabdcopy(n,dx,dy)
  implicit none
  integer :: n
  real(realk), intent(in) :: dx(n)
  real(realk), intent(inout):: dy(n)
  integer :: i         
  do I = 1,n
     dy(I) = dx(I)
  enddo
end subroutine gabdcopy

END MODULE ICHORSAVEGABMOD
