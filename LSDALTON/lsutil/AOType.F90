!> @file
!> Contains the Atomic Orbital batch structure and associated subroutines

MODULE AO_TypeType
use precision
use LSmatrix_type
!*****************************************
!*
!* 
!*
!*****************************************
!> maximum AO angular moments
INTEGER, PARAMETER      :: maxAOangmom=10
real(realk),parameter   :: ExpThr = 10E-9_realk
TYPE AOBATCHPOINTER 
TYPE(AOBATCH),pointer  :: p
END TYPE AOBATCHPOINTER

!> OBJECT CONTAINING INFORMATION ABOUT THE AOBATCHES
TYPE AOBATCH
!The type specifies the kind of primitive basis-functions
Logical                 :: TYPE_Empty
Logical                 :: TYPE_Nucleus
Logical                 :: TYPE_pCharge
Logical                 :: TYPE_elField
Logical                 :: spherical
real(realk)             :: CENTER(3)
INTEGER                 :: nPrimitives
INTEGER                 :: atom
INTEGER                 :: molecularIndex
INTEGER                 :: batch
INTEGER                 :: maxContracted
integer                 :: maxAngmom
TYPE(LSMatrix),pointer  :: pExponents
!nAngmom greater than one is for family basis-sets
integer                 :: nAngmom     
real(realk)             :: extent
INTEGER                 :: ANGMOM(maxAOangmom)
INTEGER                 :: nContracted(maxAOangmom)
INTEGER                 :: startOrbital(maxAOangmom)
INTEGER                 :: startprimOrbital(maxAOangmom)
INTEGER                 :: nOrbComp(maxAOangmom)
INTEGER                 :: nPrimOrbComp(maxAOangmom)
INTEGER                 :: nOrbitals(maxAOangmom)
TYPE(LSMatrixpointer)   :: pCC(maxAOangmom)    
INTEGER                 :: CCindex(maxAOangmom) 
INTEGER                 :: Expindex
INTEGER                 :: itype
INTEGER                 :: redtype
END TYPE AOBATCH

TYPE AOITEMPOINTER 
TYPE(AOITEM),pointer  :: p
END TYPE AOITEMPOINTER

!> OBJECT CONTAINING collection of AObatches 
TYPE AOITEM
LOGICAL                      :: EMPTY
TYPE(AOBATCH),pointer        :: BATCH(:)
INTEGER                      :: nbatches
TYPE(LSMatrix),pointer         :: CC(:) !contractioncoefficients 
INTEGER                      :: nCC
TYPE(LSMatrix),pointer         :: Exponents(:) !contractioncoefficients 
INTEGER,pointer              :: Angmom(:)
INTEGER                      :: nExp
INTEGER                      :: natoms
INTEGER                      :: ntype
INTEGER                      :: nredtype
INTEGER                      :: nbast
INTEGER                      :: nprimbast
INTEGER                      :: maxJ
INTEGER,pointer              :: ATOMICnORB(:)!size = natoms
INTEGER,pointer              :: ATOMICnBATCH(:)!size = natoms
END TYPE AOITEM

TYPE BATCHORBITALINFO
Integer :: nBatches
Integer :: nBast
Integer,pointer :: orbToBatch(:)
END TYPE BATCHORBITALINFO

public :: BATCHORBITALINFO,AOITEM,AOBATCH,AOITEMPOINTER,&
     & AOBATCHPOINTER,maxAOangmom,ExpThr

private 
contains

!Added to avoid "has no symbols" linking warning
subroutine AO_TypeType_void()
end subroutine AO_TypeType_void

END MODULE AO_TypeType

