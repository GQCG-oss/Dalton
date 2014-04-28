!> @file
!> basisinfo type module, contains also standard operation subroutines for this type
!> \brief basisinfo type module and associated subroutines for this type 
!> \author T. Kjaergaard
!> \date 2010
MODULE basis_typetype
 use precision
 use AO_typeType
 INTEGER, PARAMETER      :: maxBASISsegment=25
 INTEGER, PARAMETER      :: maxBasisSetInLIB=10
 INTEGER, PARAMETER      :: maxNumberOfChargesinLIB=10

TYPE segment
INTEGER                 :: nrow
INTEGER                 :: ncol
real(realk),pointer     :: elms(:)     
real(realk),pointer     :: Exponents(:)  !length=nrow
real(realk),pointer     :: UCCelms(:) !Unmodified ContractionCoefficients
END TYPE segment

TYPE SHELL
INTEGER                 :: nprim  !=#primitives
INTEGER                 :: norb   !=#orbitals
!TYPE(segment),pointer :: segment(:) 
TYPE(segment)           :: segment(maxBASISsegment) 
INTEGER                 :: nsegments
END TYPE SHELL

TYPE ATOMTYPEITEM
LOGICAL                  :: FAMILY !TYPE BASISSET
INTEGER                  :: nAngmom
INTEGER                  :: ToTnorb !total number of orbitals
INTEGER                  :: ToTnprim!total number of primitives
INTEGER                  :: Charge
!TYPE(ANGMOMITEM),pointer :: ANGMOM(:)
TYPE(SHELL)              :: SHELL(maxAOangmom)
Character(len=80)        :: NAME   
END TYPE ATOMTYPEITEM

TYPE BASISSETINFO
LOGICAL                    :: DunningsBasis
LOGICAL                    :: SPHERICAL
LOGICAL                    :: GCbasis
LOGICAL                    :: GCONT
INTEGER                    :: natomtypes
TYPE(ATOMTYPEITEM),pointer :: ATOMTYPE(:)
INTEGER                    :: labelindex 
!  Labelindex indicate the ordering of the different types in ATOMTYPEITEM
!
!  labelindex=n for n>0 indicate MoleculeSpecific ordering which means that 
!  the molecule%ATOM(iatom)%IDtype(n) determines which ATOMTYPE the given atom has
!  labelindex=1 is for regularMoleculeSpecific ordering 
!  labelindex=2 is for auxiliaryMoleculeSpecific ordering
!  labelindex=3 is for Complementary auxiliaryMoleculeSpecific ordering
!  labelindex=4 is for JKauxiliaryMoleculeSpecific ordering
!  labelindex=5 is for ValenceMoleculeSpecific ordering
!
!  This is the case when ATOMBASIS is used in MOLECULE.INP
!   
!  labelindex=0 indicate charge based indexing which means that all atoms 
!  share the same basisset like 6-31G or cc-pVTZ - this is the case when 
!  BASIS is used in MOLECULE.INP, This means that the chargeindex is used
!  to acces the given ATOMTYPE
INTEGER,pointer            :: Chargeindex(:)!length nChargeindex=maxcharge
!only allocated when labelindex=0
INTEGER                    :: nChargeindex
INTEGER                    :: nbast,nprimbast
Character(len=9)           :: label !REGULAR for a regular basis
END TYPE BASISSETINFO

TYPE BASISSET_PT
TYPE(BASISSETINFO),pointer :: p
END TYPE BASISSET_PT

TYPE BASISINFO
LOGICAL                   :: GCtransAlloc 
TYPE(BASISSETINFO)        :: REGULAR
TYPE(BASISSETINFO)        :: GCtrans
TYPE(BASISSETINFO)        :: AUXILIARY
TYPE(BASISSETINFO)        :: CABS
TYPE(BASISSETINFO)        :: JK
TYPE(BASISSETINFO)        :: VALENCE 
END TYPE BASISINFO

TYPE BASIS_PT
TYPE(BASISINFO),pointer :: p
END TYPE BASIS_PT

TYPE BASINF
LOGICAL             :: spherical
INTEGER             :: maxnshell !number of shells for regular basis
! so an atom with 10s8p5d1f1g would have 25 shells (75 basisfunctions)
INTEGER,pointer     :: shellangmom(:) !size kmax, angmom +1 for this shell
INTEGER,pointer     :: shellnprim(:) !size kmax, Nr. primitives
INTEGER,pointer     :: shell2atom(:) !size kmax, which atom the shell belongs to
INTEGER,pointer     :: NSTART(:)!size kmax, the shells start contracted Orbital 
INTEGER,pointer     :: priexpstart(:) !size kmax, the accumulated number of primitives
INTEGER,pointer     :: CCSTART(:)!size ushells, where the CCmatrix starts 
INTEGER,pointer     :: CCINDEX(:)!size kmax, where the CCmatrix starts 
real(realk),pointer :: X(:) !size natoms
real(realk),pointer :: Y(:) !size natoms
real(realk),pointer :: Z(:) !size natoms
real(realk),pointer :: PRIEXP(:) !size mxprim, primitive exponents 
TYPE(LSMatrix),pointer  :: CC(:) !size ushells, ContractionCoefficients for shell
TYPE(LSMatrix),pointer  :: priexpM(:) !size ushells, Exponents for shell
real(realk),pointer :: RSHEL(:) !size kmax, radii of shells
real(realk),pointer :: CENT(:,:) !size 3,kmax x,y,z coord of shell
INTEGER             :: maxAngmom,nAtoms,GRDONE,MXPRIM,Ushells 
INTEGER,pointer     :: CHARGE(:)!size natoms
END TYPE BASINF

contains
subroutine nullifyBasis(BAS)
  implicit none
  TYPE(BASISSETINFO) :: BAS
  BAS%DunningsBasis = .FALSE.
  BAS%SPHERICAL = .FALSE.
  BAS%GCbasis = .FALSE.
  BAS%GCONT = .FALSE.
  BAS%nAtomtypes=0
  nullify(BAS%ATOMTYPE)
  BAS%labelindex = 0  
  nullify(BAS%Chargeindex)
  BAS%nChargeindex = 0  
  BAS%nbast = -1
  BAS%nprimbast = -1
  BAS%label = 'Empty Bas'
end subroutine nullifyBasis

END MODULE basis_typetype

