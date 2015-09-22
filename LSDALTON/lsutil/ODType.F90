!> @file 
!> contains the overlap distribution batch structures
MODULE OD_TypeType
 use precision
 use AO_typetype
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
 TYPE ODITEM
  Integer               :: nbatches
  Integer               :: maxJ
  Logical               :: sameAOs
  TYPE(ODBATCH),pointer :: BATCH(:)
 END TYPE ODITEM
 !
 ! Type containing each induvidual OD-batch
 TYPE ODBATCH  !size  5*8+5*4+1*4+2*8 = 80
  TYPE(AOBATCHPOINTER)  :: AO(2)
  INTEGER               :: nAngmom
  INTEGER               :: nPrimitives
  INTEGER               :: maxContracted
  INTEGER(kind=short)   :: maxGAB
  LOGICAL               :: sameAO
  INTEGER               :: IA !THE AOBATCH INDEX
  INTEGER               :: IB !THE AOBATCH INDEX
  INTEGER               :: ITYPE
  INTEGER               :: redTYPE
  REAL(REALK)           :: ODcenter(3)
  REAL(REALK)           :: ODextent
 END TYPE ODBATCH
public :: ODBATCH,ODITEM
private 
contains

!Added to avoid "has no symbols" linking warning
subroutine OD_TypeType_void()
end subroutine OD_TypeType_void

END MODULE OD_TypeType

