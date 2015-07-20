!> @file
!> Contains the Atomic Orbital batch structure and associated subroutines

MODULE AO_Type
use AO_typeType
use precision
use LSmatrix_type
use lsmatrix_operations_dense
use memory_handling
private 
public :: SET_EMPTY_AO, FREE_EMPTY_AO, copy_aobatch, copy_aobatch2,&
     & PRINT_AO, PRINT_AOBATCH, FREE_AOITEM, initBatchOrbitalInfo,&
     & setBatchOrbitalInfo, freeBatchOrbitalInfo, SET_MAXJ, &
     & nullifyAOBATCH, nullifyAOITEM, COPY_AOITEM
Contains

subroutine nullifyAOITEM(AO)
implicit none
type(aoitem) :: AO

AO%EMPTY = .TRUE.
nullify(AO%BATCH)
AO%nbatches = 0
nullify(AO%CC)
AO%nCC=0
nullify(AO%Exponents)
nullify(AO%Angmom)
AO%nExp=0
AO%natoms=0
AO%ntype=0
AO%nredtype=0
AO%nbast=0
AO%nprimbast=0
AO%maxJ=0
nullify(AO%ATOMICnORB)
nullify(AO%ATOMICnBATCH)
end subroutine nullifyAOITEM

subroutine copy_aoitem(AO,AONEW)
implicit none
type(aoitem) :: AO,AONEW
integer :: I,nrow,ncol,iangmom

AONEW%EMPTY = AO%EMPTY
nullify(AONEW%BATCH)
call mem_alloc(AONEW%BATCH,AO%nbatches)
AONEW%nbatches = AO%nbatches
DO I=1,AONEW%nbatches
   call copy_aobatch2(AO%BATCH(I),AONEW%BATCH(I))
ENDDO
nullify(AONEW%CC)
call mem_alloc(AONEW%CC,AO%nCC)
AONEW%nCC=AO%nCC
DO I=1,AONEW%nCC
   nrow = AO%CC(I)%nrow
   ncol = AO%CC(I)%ncol
   CALL lsmat_dense_init(AONEW%CC(I),AO%CC(I)%nrow,AO%CC(I)%ncol)
   call dcopy(nrow*ncol,AO%CC(I)%elms,1,AONEW%CC(I)%elms,1)   
ENDDO
nullify(AONEW%Angmom)
call mem_alloc(AONEW%Angmom,SIZE(AO%Angmom))
DO I=1,SIZE(AO%Angmom)
   AONEW%angmom(I) = AO%angmom(I)
ENDDO
AONEW%nExp=AO%nExp
nullify(AONEW%Exponents)
call mem_alloc(AONEW%Exponents,AO%nExp)
DO I=1,AONEW%nEXP
   nrow = AO%Exponents(I)%nrow
   CALL lsmat_dense_init(AONEW%Exponents(I),AO%Exponents(I)%nrow,1)
   call dcopy(nrow,AO%Exponents(I)%elms,1,AONEW%Exponents(I)%elms,1)
ENDDO

DO I=1,AONEW%nbatches
   nullify(AONEW%BATCH(I)%pExponents)
   AONEW%BATCH(I)%pExponents => AONEW%Exponents(AONEW%BATCH(I)%ExpIndex)
   do iangmom = 1, AONEW%BATCH(I)%nAngmom
      nullify(AONEW%BATCH(I)%pCC(iangmom)%p)
      AONEW%BATCH(I)%pCC(iangmom)%p => AONEW%CC(AONEW%BATCH(I)%CCindex(iangmom))
   enddo
ENDDO

AONEW%natoms=AO%natoms
AONEW%ntype=AO%ntype
AONEW%nredtype=AO%nredtype
AONEW%nbast=AO%nbast
AONEW%nprimbast=AO%nprimbast
AONEW%maxJ=AO%maxJ
nullify(AONEW%ATOMICnORB)
call mem_alloc(AONEW%ATOMICnORB,SIZE(AO%ATOMICnORB))
DO I=1,SIZE(AO%ATOMICnORB)
   AONEW%ATOMICnORB(I) = AO%ATOMICnORB(I)
ENDDO
nullify(AONEW%ATOMICnBATCH)
call mem_alloc(AONEW%ATOMICnBATCH,SIZE(AO%ATOMICnBATCH))
DO I=1,SIZE(AO%ATOMICnBATCH)
   AONEW%ATOMICnBATCH(I) = AO%ATOMICnBATCH(I)
ENDDO
end subroutine copy_aoitem

subroutine nullifyAOBATCH(AO)
implicit none
type(AOITEM) :: AO
!
integer :: i,n,J
n = size(AO%BATCH)
do I=1,n
   AO%BATCH(I)%TYPE_Empty = .TRUE.
   AO%BATCH(I)%TYPE_Nucleus = .FALSE.
   AO%BATCH(I)%TYPE_pCharge = .FALSE.
   AO%BATCH(I)%TYPE_elField = .FALSE.
   AO%BATCH(I)%spherical = .FALSE.
   AO%BATCH(I)%CENTER(1) = 0.0E0_realk
   AO%BATCH(I)%CENTER(2) = 0.0E0_realk
   AO%BATCH(I)%CENTER(3) = 0.0E0_realk
   AO%BATCH(I)%nPrimitives = 0
   AO%BATCH(I)%atom = 0
   AO%BATCH(I)%molecularIndex = 0
   AO%BATCH(I)%batch = 0
   AO%BATCH(I)%maxContracted = 0
   AO%BATCH(I)%maxAngmom = 0
   nullify(AO%BATCH(I)%pExponents)
   AO%BATCH(I)%nAngmom = 0
   AO%BATCH(I)%extent = 0.0E0_realk
   do J = 1,maxAOangmom
      AO%BATCH(I)%ANGMOM(J) = 0
      AO%BATCH(I)%nContracted(J) = 0
      AO%BATCH(I)%startOrbital(J) = 0
      AO%BATCH(I)%startprimOrbital(J) = 0
      AO%BATCH(I)%nOrbComp(J) = 0
      AO%BATCH(I)%nPrimOrbComp(J) = 0
      AO%BATCH(I)%nOrbitals(J) = 0
      nullify(AO%BATCH(I)%pCC(J)%p) 
      AO%BATCH(I)%CCindex(J)  = 0
   ENDDO
   AO%BATCH(I)%Expindex = 0
   AO%BATCH(I)%itype = 0
   AO%BATCH(I)%redtype = 0
enddo
end subroutine nullifyAOBATCH

!> \brief make an empty AOITEM
!> \author T. Kjaergaard
!> \date 2010
!> \param AO the AOITEM to be built
SUBROUTINE SET_EMPTY_AO(AO)
IMPLICIT NONE
TYPE(AOITEM)  :: AO

AO%natoms = 1
CALL MEM_ALLOC(AO%ATOMICnORB,AO%natoms)  
CALL MEM_ALLOC(AO%ATOMICnBATCH,AO%natoms)
AO%ATOMICnORB(1)=1
AO%ATOMICnBATCH(1)=1
AO%nbast = 1                        
AO%nprimbast = 1
AO%empty=.TRUE.
AO%nbatches=1
AO%nCC=1
AO%nExp=1
CALL MEM_ALLOC(AO%BATCH,1)   
AO%BATCH(1)%batch=1
AO%BATCH(1)%nAngmom=1
AO%BATCH(1)%nContracted(1)=1
AO%BATCH(1)%startOrbital(1)=1
AO%BATCH(1)%nOrbComp(1)=1
AO%BATCH(1)%nPrimitives=1
AO%BATCH(1)%type_empty=.TRUE.

END SUBROUTINE SET_EMPTY_AO

!> \brief free an empty AOITEM
!> \author T. Kjaergaard
!> \date 2010
!> \param AO free the AOITEM
SUBROUTINE FREE_EMPTY_AO(AO)
IMPLICIT NONE
TYPE(AOITEM)  :: AO

CALL MEM_DEALLOC(AO%ATOMICnORB)  
CALL MEM_DEALLOC(AO%ATOMICnBATCH)
CALL MEM_DEALLOC(AO%BATCH)   

END SUBROUTINE FREE_EMPTY_AO

subroutine copy_aobatch(AO1,AO2)
implicit none
type(AOBATCH) :: AO1,AO2
!
integer :: iangmom,nrow,ncont

call copy_aobatch2(AO1,AO2)

nullify(AO2%pExponents)
allocate(AO2%pExponents)
nrow = AO1%pExponents%nrow
CALL lsmat_dense_init(AO2%pExponents,nrow,1)
AO2%pExponents%elms(1:nrow) = AO1%pExponents%elms(1:nrow)

do iangmom = 1, AO1%nAngmom
   nrow = AO1%pCC(iangmom)%p%nrow
   ncont = AO1%pCC(iangmom)%p%ncol
   nullify(AO2%pCC(iangmom)%p)
   allocate(AO2%pCC(iangmom)%p)
   CALL lsmat_dense_init(AO2%pCC(iangmom)%p,nrow,ncont)
   call dcopy(nrow*nCont,AO1%pCC(iangmom)%p%elms,1,AO2%pCC(iangmom)%p%elms,1)
enddo

end subroutine copy_aobatch

subroutine copy_aobatch2(AO1,AO2)
implicit none
type(AOBATCH),intent(in) :: AO1
type(AOBATCH),intent(inout) :: AO2
AO2%TYPE_Empty = AO1%TYPE_Empty
AO2%TYPE_pCharge = AO1%TYPE_pCharge
AO2%TYPE_Nucleus = AO1%TYPE_Nucleus
AO2%TYPE_elField = AO1%TYPE_elField
AO2%spherical = AO1%spherical
AO2%CENTER = AO1%CENTER
AO2%nPrimitives = AO1%nPrimitives
AO2%atom= AO1%atom
AO2%molecularIndex = AO1% molecularIndex
AO2%batch = AO1%batch
AO2%maxContracted = AO1%maxContracted
AO2%maxAngmom = AO1%maxAngmom
AO2%nAngmom = AO1%nAngmom
AO2%extent = AO1%extent
AO2%ANGMOM = AO1%ANGMOM
AO2%nContracted = AO1%nContracted
AO2%startOrbital = AO1%startOrbital
AO2%startprimOrbital = AO1%startprimOrbital
AO2%nOrbComp = AO1%nOrbComp
AO2%nPrimOrbComp = AO1%nPrimOrbComp
AO2%nOrbitals = AO1%nOrbitals
AO2%CCindex = AO1%CCindex
AO2%Expindex = AO1%Expindex
AO2%itype = AO1%itype
AO2%redtype = AO1%redtype
end subroutine copy_aobatch2
!> \brief Print the AOITEM structure
!> \author T. Kjaergaard
!> \date 2010
!> \param LUPRI the logical unit number for the output file
!> \param AO the Atomic Orbital item to be printet
SUBROUTINE PRINT_AO(LUPRI,AO)
implicit none
TYPE(AOITEM)        :: AO
INTEGER             :: I
INTEGER             :: LUPRI

   WRITE(LUPRI,*) '                     '
   WRITE(LUPRI,'(A)')'PRINTING AO BATCH '
   WRITE(LUPRI,'(3X,A,I7)')'# ATOMS             ',AO%natoms
   WRITE(LUPRI,'(3X,A,I7)')'# BASISFUNCTIONS    ',AO%nbast
   WRITE(LUPRI,'(3X,A,I7)')'# PRIMITIVEFUNCTIONS',AO%nprimbast
   WRITE(LUPRI,'(3X,A,I7)')'# BATCHES           ',AO%nbatches
   WRITE(LUPRI,'(3X,A6,2X,A10,2X,A11)')'ATOM  ','# ORBITALS','# AOBATCHES'
   DO I=1,AO%natoms
      WRITE(LUPRI,'(3X,I6,6X,I6,7X,I6)')I,AO%ATOMICnORB(I),AO%ATOMICnBATCH(I)
   ENDDO
   DO I=1,AO%nbatches
      WRITE(LUPRI,'(5X,A,I4)')'AOBATCH NUMBER',I
      CALL PRINT_AOBATCH(AO%BATCH(I),LUPRI)
   ENDDO
   WRITE(LUPRI,'(3X,A18,L7)')'Empty             ',AO%Empty
   WRITE(LUPRI,'(3X,A18,I7)')'maxJ              ',AO%maxJ
   WRITE(LUPRI,'(3X,A18,I7)')'nredtype          ',AO%nredtype
   WRITE(LUPRI,'(3X,A18,I7)')'ntype             ',AO%ntype
   WRITE(LUPRI,*)' '

END SUBROUTINE PRINT_AO

!> \brief Print the AOBATCH structure
!> \author T. Kjaergaard
!> \date 2010
!> \param IUNIT the logical unit number for the output file
!> \param AOB the Atomic Orbital batch to be printet
SUBROUTINE PRINT_AOBATCH(AOB,IUNIT)
 IMPLICIT NONE
 TYPE(AOBATCH)  :: AOB
 Integer        :: IUNIT
!
 Integer        :: i,k,nangmom,nprimitives,ncol
!
 WRITE(IUNIT,'(5X,A,L1)') 'Type of basis-function: Empty    :',AOB%type_empty
 WRITE(IUNIT,'(5X,A,L1)') 'Type of basis-function: Nucleus  :',AOB%type_nucleus
 WRITE(IUNIT,'(5X,A,L1)') 'Type of basis-function: pCharge  :',AOB%type_pCharge
 IF (AOB%spherical) THEN
   WRITE(IUNIT,'(5X,A)') 'The orbitals are spherical'
 ELSE
   WRITE(IUNIT,'(5X,A)') 'The orbitals are cartesian'
 ENDIF
 WRITE(IUNIT,'(5X,A,2X,3F12.8)')'Cartesian center (A):',&
       &AOB%CENTER(1),AOB%CENTER(2),AOB%CENTER(3)
 nPrimitives=AOB%nPrimitives
 WRITE(IUNIT,'(5X,A,I3)')    'Number of primitives,       nPrimitives   = ', AOB%nPrimitives
 WRITE(IUNIT,'(5X,A,I3)')    'Max. # of contracted,       maxContracted = ', AOB%maxContracted
 WRITE(IUNIT,'(5X,A,F12.8)') 'The max. extent of the AO-batch,   extent = ', AOB%extent
 WRITE(IUNIT,'(5X,A,I3)')    'Max. angular momentum,      maxAngmom     = ', AOB%maxAngmom
 WRITE(IUNIT,'(5X,A,I3)')    '# of ang. mom. blocks,      nAngmom       = ', AOB%nAngmom
 WRITE(IUNIT,'(5X,A,I3)')    'Atom index                                = ', AOB%atom
 WRITE(IUNIT,'(5X,A,I3)')    'batch index                               = ', AOB%batch
 WRITE(IUNIT,'(5X,A,I3)')    'type index                                = ', AOB%itype
 WRITE(IUNIT,'(5X,A,I3)')    'reduced type index                        = ', AOB%redtype
 WRITE(IUNIT,'(3X,A)')       '------------------ Orbital block information -------------------------------'
 WRITE(IUNIT,'(3X,A)')       '-    Block#  Ang.mom.  #cont.   1.orb.   #orb.   #orb.c.  #pri.o.c  1.porb -'
 DO I=1,AOB%nAngmom
   WRITE(IUNIT,'(3X,8I9)') I,AOB%angmom(I),AOB%nContracted(I),AOB%startOrbital(I),&
      & AOB%nOrbitals(I),AOB%nOrbComp(I),AOB%nPrimOrbComp(I),AOB%startprimOrbital(I)
 ENDDO
 WRITE(IUNIT,'(3X,A)')       '----------------------------------------------------------------------------'
 WRITE(IUNIT,'(3X,A)')       '----------- Exponents -----------'
 call LSMAT_dense_PRINT(AOB%pExponents, 1, nPrimitives, 1, 1,IUNIT)
 WRITE(IUNIT,'(3X,A)')       '---------------------------------'
 DO K=1,AOB%nAngmom
    nAngmom=AOB%Angmom(K)
    WRITE(IUNIT,'(3X,A,I3,A,I3,A)') '--- Contraction coefficient for ang.mom: ',nAngmom,' - CCindex:',AOB%CCindex(K),' ---'
    ncol=AOB%pCC(K)%p%ncol
    CALL LSMAT_DENSE_PRINT(AOB%pCC(K)%p,1,nPrimitives,1,ncol,IUNIT)
 WRITE(IUNIT,'(3X,A)')         '-----------------------------------------------------------------------'
 ENDDO

END SUBROUTINE PRINT_AOBATCH

!> \brief free the AOitem structure
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE FREE_AOITEM(LUPRI,AO)
IMPLICIT NONE
TYPE(AOITEM)              :: AO
INTEGER                   :: I,J,LUPRI

CALL MEM_DEALLOC(AO%ATOMICnORB)
CALL MEM_DEALLOC(AO%ATOMICnBATCH)

DO I=1,AO%nCC
   CALL LSMAT_DENSE_FREE(AO%CC(I))
ENDDO
DO I=1,AO%nExp
   CALL LSMAT_DENSE_FREE(AO%Exponents(I))
ENDDO
DO I=1,AO%nbatches
   DO J=1,AO%BATCH(I)%nAngmom
      NULLIFY(AO%BATCH(I)%pCC(J)%p)
   ENDDO
   NULLIFY(AO%BATCH(I)%pExponents)
ENDDO
CALL MEM_DEALLOC(AO%Angmom)
CALL MEM_DEALLOC(AO%BATCH)
CALL MEM_DEALLOC(AO%CC)
CALL MEM_DEALLOC(AO%Exponents)

END SUBROUTINE FREE_AOITEM

!> \brief Initialize a batchorbitalinfo type
!> \author S. Reine
!> \date 18-03-2010
!> \param BO The batchorbitalinfo
!> \param nbast The number of orbitals/basis functions
SUBROUTINE initBatchOrbitalInfo(BO,nBast)
use memory_handling
implicit none
TYPE(BATCHORBITALINFO),intent(INOUT) :: BO
!this is really intent out but this will diassociate the pointer
Integer,intent(IN)                 :: nBast

BO%nBatches = 0
BO%nBast = nBast
call mem_alloc(BO%orbToBatch,nBast)

END SUBROUTINE initBatchOrbitalInfo

!> \brief Set up the batchorbitalinfo
!> \author S. Reine
!> \date 18-03-2010
!> \param BO The batchorbitalinfo
!> \param AO The AO-batch
!> \param lupri Deafult output unit
SUBROUTINE setBatchOrbitalInfo(BO,AO,lupri)
implicit none
TYPE(BATCHORBITALINFO),intent(INOUT) :: BO 
!this is really intent out but this will diassociate the pointer
TYPE(AOITEM),intent(IN)            :: AO
Integer,intent(IN)                 :: lupri
!
Integer :: iBatch,iAngmom,startOrb,numOrb,iOrb

IF (BO%nBast.NE.AO%nBast) THEN
  WRITE(LUPRI,'(1X,A,I8,A,I8)') 'Error in setBatchOrbitalInfo. AO%nBast = ',AO%nBast,&
    & ' and BO%nBast = ',BO%nBast
  CALL LSQUIT('Error in setBatchOrbitalInfo. nBast mismatch!',-1)
ELSE
  BO%nBatches = AO%nBatches
  DO iBatch=1,AO%nBatches
     !iatom = AO%Batch(iBatch)%atom
    DO iAngmom=1,AO%Batch(iBatch)%nAngmom
      startOrb = AO%Batch(iBatch)%startOrbital(iAngmom)
      numOrb   = AO%Batch(iBatch)%nOrbitals(iAngmom)
      DO iOrb=startOrb,startOrb+numOrb-1
        BO%orbToBatch(iOrb) = iBatch
      ENDDO
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE setBatchOrbitalInfo

!> \brief Free a batchorbitalinfo type
!> \author S. Reine
!> \date 18-03-2010
!> \param BO The batchorbitalinfo
SUBROUTINE freeBatchOrbitalInfo(BO)
use memory_handling
implicit none
TYPE(BATCHORBITALINFO),intent(INOUT) :: BO

IF (.NOT.ASSOCIATED(BO%orbToBatch)) THEN
  CALL LSQUIT('Error in freeBatchOrbitalInfo. orbToBatch not associated!',-1)
ELSE
  BO%nBatches = 0
  BO%nBast = 0
  call mem_dealloc(BO%orbToBatch)
ENDIF
END SUBROUTINE freeBatchOrbitalInfo

!> \brief set the maximum angularmoment of this AOtype
!> \author T. Kjaergaard
!> \date 2011
!> \param LUPRI the logical unit number for the output file
!> \param AO the Atomic Orbital item
SUBROUTINE SET_MAXJ(AO,LUPRI)
implicit none
INTEGER             :: LUPRI
TYPE(AOITEM)        :: AO
!
INTEGER             :: I,J,maxJ
maxJ = 0
DO I=1,AO%nbatches
   DO J=1,AO%BATCH(I)%nAngmom
      maxJ = MAX(maxJ,AO%BATCH(I)%angmom(J))
   ENDDO
ENDDO
AO%maxJ = maxJ
END SUBROUTINE SET_MAXJ

END MODULE AO_Type

