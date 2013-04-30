!> @file
!> molecule type module, contains also standard operation subroutines for this type
!> \author S. Reine and T.Kjaergaard
!> \date 2010-02-21
MODULE molecule_typetype
use precision
!*****************************************
!*
!* OBJECT CONTAINING INFORMATION ABOUT THE MOLECULE
!*
!*****************************************
TYPE ATOMITEM
Integer           :: Isotope !Isotope number acc.to abundancy 
Character(len=4)  :: Name    !Name of atom
real(realk)       :: Mass    !Atomic mass 
real(realk)       :: CovRad  ! Covalent radius
real(realk)       :: Frag    ! Assigned fragment
real(realk)       :: CENTER(3)
Integer           :: Atomic_number !Atomic number
real(realk)       :: Charge        !Atomic Charge
integer           :: nbasis
Character(len=9)  :: basislabel(5) !REGULAR,AUXILIARY,CABS,JK,VALENCE
INTEGER           :: Basisindex(5) !which set in the BASISSETLIBRARY
INTEGER           :: IDtype(5) !A unique identifier - identifying the type
LOGICAL           :: Phantom !true if basisfunction but no actual atom 
LOGICAL           :: Pointcharge !true if no basis functions
! THE FOLLOWING ARE ADDED AFTER BUILD BASIS TO DO PARALLELIZATION
!INTEGER           :: SphOrbREG !Spherical orbitals for Regular basis
!INTEGER           :: CAROrbREG !Cartesian orbitals
INTEGER           :: nContOrbREG !# contracted orbitals
INTEGER           :: nPrimOrbREG !# primitives orbitals
!INTEGER           :: SphOrbAUX !Spherical orbitals for Aux bas
!INTEGER           :: CAROrbAUX !Cartesian orbitals
INTEGER           :: nContOrbAUX !# contracted orbitals
INTEGER           :: nPrimOrbAUX !# primitives orbitals
INTEGER           :: nContOrbCABS !# contracted orbitals
INTEGER           :: nPrimOrbCABS !# primitives orbitals
INTEGER           :: nContOrbJK !# contracted orbitals
INTEGER           :: nPrimOrbJK !# primitives orbitals
INTEGER           :: nContOrbVAL !# contracted orbitals for valence basis
INTEGER           :: nPrimOrbVAL !# primitives orbitals for valence basis
END TYPE ATOMITEM

TYPE MOLECULEINFO
Character(len=22)    :: label
TYPE(ATOMITEM), pointer  :: ATOM(:) !length = nAtomtypes
INTEGER              :: nAtoms
INTEGER              :: nAtomsNPC !nAtoms NOT including pointcharges
INTEGER              :: nelectrons
real(realk)          :: charge !molecular charge
! THE FOLLOWING ARE ADDED AFTER BUILD BASIS
INTEGER              :: nbastREG
INTEGER              :: nbastAUX
INTEGER              :: nbastCABS
INTEGER              :: nbastJK
INTEGER              :: nbastVAL
INTEGER              :: nprimbastREG
INTEGER              :: nprimbastAUX
INTEGER              :: nprimbastCABS
INTEGER              :: nprimbastJK
INTEGER              :: nprimbastVAL
logical              :: pointMolecule
END TYPE MOLECULEINFO

TYPE MOLECULE_PT
TYPE(MOLECULEINFO),pointer :: p
END TYPE MOLECULE_PT

TYPE MOLECULARORBITALINFO
INTEGER         :: nAtoms
INTEGER         :: nBastReg
INTEGER         :: nBastAux
INTEGER,POINTER :: numAtomicOrbitalsReg(:)
INTEGER,POINTER :: startAtomicOrbitalsReg(:)
INTEGER,POINTER :: endAtomicOrbitalsReg(:)
INTEGER,POINTER :: numAtomicOrbitalsAux(:)
INTEGER,POINTER :: startAtomicOrbitalsAux(:)
INTEGER,POINTER :: endAtomicOrbitalsAux(:)
END TYPE MOLECULARORBITALINFO

contains

!Added to avoid "has no symbols" linking warning
subroutine molecule_typetype_void()
end subroutine molecule_typetype_void

END MODULE molecule_typetype
