!> @file 
!> contains the overlap distribution batch structures
MODULE OD_Type
 use precision
 use OD_TypeType
 use AO_typetype
 use AO_type
 use memory_handling
 private
 public :: getODscreening, FREE_ODitem, PRINT_OD, PRINT_ODBATCH
! use ls_util
 !************************************************
 !*
 !* OBJECT CONTAINING INFORMATION ABOUT OD-BATCHES
 !*
 !************************************************
 ! An OD-batch is a set of:
 !  i)   product overlaps between two batches of basis-functions -
 !       where each AO-batch belong to a given center, and
 !       share a common set of primitive basis-functions
 !  ii)  AO basis-functions (beloning to one
 !       center and with a set of primitive functions).
 !  iii) an empty batch (used for debugging purposes only -
 !       useful when testing contruction of E-coefficients).
 !
 ! Type containing collection of OD-batches
Contains
!> \brief determine if the integral should be screening based on overlap distance screening
!> \author T. Kjaergaard
!> \date 2010
!> \param BatchA the 1. AO batch
!> \param BatchB the 2. AO batch
!> \param screen true if the integral should be screened away
SUBROUTINE getODscreening(BATCHA,BATCHB,screen) 
implicit none
type(AOBATCH) :: batchA,batchB
logical :: screen
!
real(realk) :: sumExtent2,distance,distance2
integer :: I

IF ((.NOT.BATCHA%type_Empty).AND.(.NOT.BATCHB%type_Empty)) THEN
   sumExtent2 = BATCHA%extent + BATCHB%extent
   sumExtent2 = sumExtent2*sumExtent2
   distance2 = 0E0_realk
   DO I=1,3
      distance  = BATCHA%CENTER(I) - BATCHB%CENTER(I)
      distance2 = distance2 + distance*distance
   ENDDO
   IF (distance2.GT.sumExtent2) screen = .TRUE.
ELSE
   screen = .FALSE.
ENDIF
END SUBROUTINE getODscreening

!> \brief free OD-item
!> \author T. Kjaergaard
!> \date 2010
!> \param OD the OD item
SUBROUTINE FREE_ODitem(OD)
  use memory_handling
IMPLICIT NONE
 TYPE(ODITEM)   :: OD
 Integer        :: SIZE
 integer (kind=long) :: nsize

 SIZE = 5*mem_realsize+5*mem_intsize+1*mem_logicalsize
 nsize = SIZE*OD%nbatches 
 call mem_deallocated_mem_ODitem(nsize)
 call mem_dealloc(OD%BATCH)
END SUBROUTINE FREE_ODitem

!> \brief Print OD-item
!> \author T. Kjaergaard
!> \date 2010
!> \param OD the OD item
!> \param IUNIT the logical unit number for the output 
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE PRINT_OD(OD,IUNIT,IPRINT)
IMPLICIT NONE
  TYPE(ODITEM)      :: OD
  INTEGER           :: IUNIT,IPRINT
!
  INTEGER  :: IOD
!
  IF (IPRINT.GT. 2) THEN
!    WRITE(TEXT,'(2A)') 'ODBATCHES ',OD%IDENTIFIER
!    CALL HEADER(TEXT,-1)
    WRITE(IUNIT,'(1X,A,I8)') 'Number of batches:',OD%nbatches
    IF (OD%sameAOs) THEN
      WRITE(IUNIT,'(1X,A)') 'Batches built from identical AO-batches'
    ELSE
      WRITE(IUNIT,'(1X,A)') 'Batches built from different AO-batches'
    ENDIF
    IF (IPRINT.GT. 3) THEN
      DO IOD=1,OD%nbatches
        WRITE(IUNIT,'(1X,A,I8,A)') '*** BATCH number ',IOD,' ***'
        CALL PRINT_ODBATCH(OD%BATCH(IOD),IUNIT,IPRINT)
        WRITE(IUNIT,'(1X,A)') '************************'
      ENDDO
    ENDIF
  ENDIF
END SUBROUTINE PRINT_OD
!
!> \brief Print induvidual OD-batch
!> \author T. Kjaergaard
!> \date 2010
!> \param ODB the OD batch
!> \param IUNIT the logical unit number for the output 
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE PRINT_ODBATCH(ODB,IUNIT,IPRINT)
IMPLICIT NONE
  TYPE(ODBATCH)     :: ODB
  INTEGER           :: IUNIT,IPRINT
!
  WRITE(IUNIT,'(5X,A,I3)') 'Number of primitives,      nPrimitives   = ', ODB%nPrimitives
  WRITE(IUNIT,'(5X,A,I3)') 'Number of angular blocks,  nAngmom       = ', ODB%nAngmom
  WRITE(IUNIT,'(5X,A,I3)') 'Max. number of contracted, maxContracted = ', ODB%maxContracted
  WRITE(IUNIT,'(5X,A,I4)')'Max. Schwarz(GAB) matrix element         = ', ODB%maxGAB
  WRITE(IUNIT,'(5X,A,L1)')'Same center and primitives = ', ODB%sameAO

  IF (IPRINT.GT. 5) THEN
    WRITE(IUNIT,'(3X,A)') '* First AO-batch  *'
    CALL PRINT_AOBATCH(ODB%AO(1)%p,IUNIT)
    WRITE(IUNIT,'(3X,A)') '* Second AO-batch *'
    CALL PRINT_AOBATCH(ODB%AO(2)%p,IUNIT)
  ENDIF

END SUBROUTINE PRINT_ODBATCH

END MODULE OD_Type

