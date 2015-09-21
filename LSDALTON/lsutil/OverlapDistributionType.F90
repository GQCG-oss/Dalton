!> @file
!> Contains the orbital and overlap types, and subroutines to set these up
MODULE OverlapType
use precision
use LSmatrix_type
use AO_TypeType

TYPE SPHMAT
REAL(REALK),pointer  :: elms(:)
END TYPE SPHMAT

TYPE SPHMATPOINTER
TYPE(SPHMAT),pointer :: p
END TYPE SPHMATPOINTER

TYPE Orbital
integer,pointer               :: angmom(:) !1:nangmom
integer,pointer               :: nContracted(:)
integer,pointer               :: startOrbital(:)
integer,pointer               :: startLocOrb(:)
integer,pointer               :: nOrbComp(:)
integer,pointer               :: nOrbitals(:)
!Dimensions according to nPrimitives
real(realk),pointer           :: exponents(:)
!Dimensions according to nPrimitives,maxContracted,nAngmom
!type(lsmatrixpointer),pointer :: CC(:)
type(lsmatrixpointer)         :: CC(maxAOangmom)
!Only one of these types can be true
Integer                       :: nAngmom
Integer                       :: nPrimitives
Integer                       :: totOrbitals
Logical                       :: TYPE_Empty
Logical                       :: TYPE_Hermite
Logical                       :: TYPE_Cartesian
Logical                       :: TYPE_Nucleus
Logical                       :: TYPE_pCharge
Logical                       :: spherical
Logical                       :: FTUVorb
Integer                       :: maxAngmom
Integer                       :: nPasses
Integer                       :: maxContracted
Integer                       :: CCidentifier
Real(realk)                   :: center(3)       !Cartesian center
END TYPE Orbital

!****** OVERLAPPOINTER
TYPE Overlappointer
TYPE(Overlap),pointer :: p
END TYPE Overlappointer

!****** OVERLAP
TYPE Overlap
!Dimensions according to nPrimitives*nPasses
Real(realk),pointer :: center(:)
Real(realk),pointer :: distance12(:)
Real(realk),pointer :: preExpFac(:)
Integer,pointer     :: iprim1(:)    !poniter to primitive of orbital 1
Integer,pointer     :: iprim2(:)    !poniter to primitive of orbital 2
integer,pointer     :: Orb1atom(:)
integer,pointer     :: Orb2atom(:)
integer,pointer     :: Orb1batch(:)
integer,pointer     :: Orb2batch(:)
integer,pointer     :: Orb1mol(:)
integer,pointer     :: Orb2mol(:)
Integer             :: nAngmom
Integer             :: nPrimitives
Logical             :: ContractBasis 
Logical             :: single
!only one of these types can be true
Logical             :: TYPE_Empty
Logical             :: TYPE_Hermite
Logical             :: TYPE_Hermite_single
Logical             :: TYPE_Cartesian
Logical             :: TYPE_Cartesian_single
Logical             :: TYPE_Nucleus
Logical             :: TYPE_FTUV
!
Logical             :: segmented
real(realk)         :: ODcenter(3)
real(realk)         :: ODextent
Logical             :: sphericalEcoeff
TYPE(Orbital)       :: orbital1,orbital2
logical             :: sameAO
Integer,pointer     :: totOrbitals(:) !(startGeoOrder:endGeoOrder)
Integer             :: nTUV
Integer             :: maxContracted
integer(kind=short) :: maxGab
Real(realk)         :: mbie(2)
!minAngmom and maxAngmom represents the actual angular momenta of the orbitals
Integer             :: minAngmom
Integer             :: maxAngmom
!The order of differentiation (0 for undifferentiated, 1 for first-order differentiation etc.)
Integer             :: startGeoOrder
Integer             :: endGeoOrder
!The number of differentiated components (for gradients: 3 for a single AO gradient 
!and 6 for a product of two AOs)
Integer             :: nGeoDerivComp 
Integer             :: nCartesianMomentComp
!startAngmom and endAngmom represents the angular momenta of the inner 
!Hermite integrals (will typically differ from min and max when dealing
!with differentiated integrals, etc.)
Integer             :: startAngmom
Integer             :: endAngmom
!ToDo: for ETUV's not attached, remember to increase dimension here, and
!      update the building of Ecoefficients
Real(realk)         :: squaredDistance
Real(realk),pointer :: exponents(:)
Real(realk),pointer :: reducedExponents(:)
!Dimensions according to nAngmom
Integer,pointer     :: angmom(:)
Integer,pointer     :: indexAng1(:)
Integer,pointer     :: indexAng2(:)
Integer,pointer     :: nOrbitals(:)
Integer,pointer     :: nOrbComp(:)
Integer,pointer     :: nContracted(:)
Integer             :: batchindex1
Integer             :: batchindex2
!McMurcie-Davidson E-coefficients (and/or F-coefficients in case of J-engine) 
Logical             :: FTUVisAlloc
Integer,pointer     :: lenETUV(:)
Real(realk),pointer :: FTUV(:)
integer             :: nPrimAlloc
integer             :: nTUVAlloc
integer             :: nAngAlloc
!Pass specific variables (passes are used to increase perforance by 
!collecting OD-batches of similar structure - this decrease the time
!spent in the different contraction/matrix multiplies steps, as 
!well as reducing the overhead related to memory/cache loadings and
!subroutine calls).
Integer             :: nPasses
Integer             :: passType
!
Integer             :: cmorder ! order of cartesianmoments
Integer             :: mmorder ! order of multipolmoments
Logical             :: do_mulmom
Integer             :: magderiv
END TYPE Overlap

TYPE TUVitem
INTEGER,pointer        :: TUVindex(:,:,:)
INTEGER                :: nTABFJW1
INTEGER                :: nTABFJW2
Real(realk),pointer    :: TABFJW(:,:)
INTEGER,pointer        :: Tindex(:)
INTEGER,pointer        :: Uindex(:)
INTEGER,pointer        :: Vindex(:)
TYPE(SPHMAT),pointer   :: SPH_MAT(:)
INTEGER                :: nSPHMAT
INTEGER,pointer        :: SYMindex(:)
INTEGER                :: iPQxyz
!integer                :: maxJ
!integer                :: ntuvmax
END TYPE TUVitem

TYPE IntegerPointer 
INTEGER                :: dim 
integer,pointer        :: elms(:) 
END type IntegerPointer

public :: SPHMAT,SPHMATPOINTER,Orbital,Overlappointer,&
     & Overlap,TUVitem,IntegerPointer 

private

contains

!Added to avoid "has no symbols" linking warning
subroutine OverlapType_void()
end subroutine OverlapType_void

END MODULE OverlapType
