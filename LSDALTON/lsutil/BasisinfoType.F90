!> @file
!> basisinfo type module, contains also standard operation subroutines for this type
!> \brief basisinfo type module and associated subroutines for this type 
!> \author T. Kjaergaard
!> \date 2010
MODULE basis_typetype
 use precision
 use AO_typeType
 use LSmatrix_type
 INTEGER, PARAMETER      :: maxBASISsegment=25
 INTEGER, PARAMETER      :: maxBasisSetInLIB=10
 INTEGER, PARAMETER      :: maxNumberOfChargesinLIB=10
!Remember to modify nBasisBasParam if you add a basis set to the list. 
 Integer,parameter :: nBasisBasParam=7
 Integer,parameter :: RegBasParam=1
 Integer,parameter :: AUXBasParam=2
 Integer,parameter :: CABBasParam=3
 Integer,parameter :: JKBasParam=4
 Integer,parameter :: VALBasParam=5
 Integer,parameter :: GCTBasParam=6
 Integer,parameter :: ADMBasParam=7

 character(len=9),parameter :: BasParamLABEL(nBasisBasParam) = &
      & (/'REGULAR  ','AUXILIARY','CABS     ','JKAUX    ',&
      & 'VALENCE  ','GCTRANS  ','ADMM     '/)

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
REAL(REALK)                :: GeminalScalingFactor
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

!todo change to BASIS(nBasis) using Regularbasis=1,GCtransbasis=2,..
TYPE BASISINFO
LOGICAL                 :: WBASIS(nBasisBasParam)!Which basis is used/allocated
TYPE(BASISSETINFO)      :: BINFO(nBasisBasParam)
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

TYPE BRAKEBASINFO
INTEGER             :: OriginalType
INTEGER             :: OriginalnAngmom
INTEGER             :: nAngmom
INTEGER,pointer     :: OriginalnOrb(:)
INTEGER,pointer     :: nOrb(:)
INTEGER,pointer     :: Orb(:,:) 
INTEGER,pointer     :: angmom(:)
END TYPE BRAKEBASINFO

public ::  maxBASISsegment,maxBasisSetInLIB,maxNumberOfChargesinLIB,&
     & nBasisBasParam,RegBasParam,AUXBasParam,&
     & CABBasParam,JKBasParam,VALBasParam,&
     & GCTBasParam,ADMBasParam,BasParamLABEL,&
     & BRAKEBASINFO,BASINF,BASISINFO,BASIS_PT,&
     & BASISSETINFO,BASISSET_PT,SHELL,ATOMTYPEITEM,&
     & segment,nullifyBasisset,nullifyMainBasis,&
     & nullifyAtomType,nullifyShell,nullifySegment

private
contains
subroutine nullifyBasisset(BAS)
  implicit none
  TYPE(BASISSETINFO) :: BAS
  BAS%DunningsBasis = .FALSE.
  BAS%GeminalScalingFactor = 1.0E0_realk
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
end subroutine nullifyBasisset

subroutine nullifyMainBasis(BAS)
  implicit none
  TYPE(BASISINFO) :: BAS
  integer :: I
  BAS%WBASIS = .FALSE.
  DO I=1,nBasisBasParam
     call nullifyBasisset(BAS%BINFO(I))
  ENDDO
end subroutine nullifyMainBasis

subroutine nullifyAtomType(AT)
implicit none
type(atomtypeitem) :: AT
!
integer :: j
AT%FAMILY = .FALSE.
AT%nAngmom=0
AT%ToTnorb=0
AT%ToTnprim=0
AT%Charge=0
do J=1,80
   AT%NAME(J:J) = ' '
enddo
do J=1,maxAOangmom
   call nullifySHELL(AT%SHELL(J))
end do
end subroutine nullifyAtomType

subroutine nullifyShell(sh)
implicit none
type(SHELL) :: SH
!
integer :: I

SH%nprim = 0
SH%norb = 0 
SH%nsegments = 0 
do I=1,maxBASISsegment
   call nullifySegment(SH%segment(I))
enddo
end subroutine nullifyShell

subroutine nullifySegment(se)
implicit none
type(Segment) :: Se
Se%nrow = 0
Se%ncol = 0 
nullify(Se%elms)
nullify(Se%Exponents)
nullify(Se%UCCelms)
end subroutine nullifySegment

END MODULE basis_typetype

