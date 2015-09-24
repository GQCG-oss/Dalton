!> @file
!> molecule type module, contains also standard operation subroutines for this type
!> \author S. Reine and T.Kjaergaard
!> \date 2010-02-21
MODULE molecule_typetype
use precision
use basis_typetype
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
!integer           :: nbasis
!REGULAR,AUXILIARY,CABS,JK,ADMM,VALENCE,...
Character(len=9)  :: basislabel(nBasisBasParam) 
!which set in the BASISSETLIBRARY
INTEGER           :: Basisindex(nBasisBasParam) 
!A unique identifier - identifying the type
INTEGER           :: IDtype(nBasisBasParam) 
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
INTEGER           :: nContOrbADMM !# contracted orbitals
INTEGER           :: nPrimOrbADMM !# primitives orbitals
INTEGER           :: nContOrbVAL !# contracted orbitals for valence basis
INTEGER           :: nPrimOrbVAL !# primitives orbitals for valence basis
INTEGER           :: molecularIndex !# consecutive atomic index (for full the molecule)
INTEGER           :: SubSystemIndex !(index in Moleculeinfo%SubsystemLabel
END TYPE ATOMITEM

TYPE MOLECULEINFO
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
INTEGER              :: nbastADMM
INTEGER              :: nbastVAL
INTEGER              :: nprimbastREG
INTEGER              :: nprimbastAUX
INTEGER              :: nprimbastCABS
INTEGER              :: nprimbastJK
INTEGER              :: nprimbastADMM
INTEGER              :: nprimbastVAL
logical              :: pointMolecule
INTEGER              :: nSubSystems
Character(len=22)    :: label
Character(len=80),pointer :: SubsystemLabel(:)
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

public :: MOLECULARORBITALINFO,MOLECULE_PT,MOLECULEINFO,&
     & ATOMITEM

private

contains

!Added to avoid "has no symbols" linking warning
subroutine molecule_typetype_void()
end subroutine molecule_typetype_void

END MODULE molecule_typetype
