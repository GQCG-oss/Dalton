!> @file
!> molecule type module, contains also standard operation subroutines for this type
!> \author S. Reine and T.Kjaergaard
!> \date 2010-02-21
MODULE molecule_type
use molecule_typetype
use precision
use memory_handling
use basis_typetype
private
public ::  DETERMINE_MAXCOOR, PRINT_MOL, &
     & write_atom_to_disk, &
     & read_atom_from_disk, &
     & write_moleculeinfo_to_disk, &
     & read_moleculeinfo_from_disk, &
     & init_Moleculeinfo, &
     & free_Moleculeinfo, &
     & build_atomicmolecule, &
     & build_empty_molecule, &
     & build_pointmolecule, &
     & build_molecule_from_coordlist, &
     & copy_molecule,copy_atom, &
     & molecule_get_num_atoms, &
     & molecule_get_atom_coord, &
     & molecule_get_atom_charge

CONTAINS
!#################################################################
!#
!# SUBROUTINES THAT ONLY AFFECT THE TYPES DEFINED IN THIS FILE
!# init_moleculeinfo
!# free_Moleculeinfo
!# build_atomicmolecule
!# build_pointmolecule
!# COPY_ATOM
!# COPY_MOLECULE
!# MPICOPY_ATOM
!# MPICOPY_MOLECULE
!# write_moleculeinfo_to_disk
!# read_moleculeinfo_from_disk
!# 
!################################################################

!> \brief print bare molecule
!> \author T. Kjaergaard
!> \date  2010-02-24
SUBROUTINE determine_maxCoor(molecule,maxcoor)
TYPE(MOLECULEINFO),intent(in) :: MOLECULE
real(realk),intent(out) :: maxcoor
!
INTEGER :: I 
real(realk) :: t

maxcoor = TINY(t)
DO I=1,MOLECULE%nAtoms
 maxcoor = MAX(maxcoor,MOLECULE%ATOM(I)%CENTER(1),MOLECULE%ATOM(I)%CENTER(2),MOLECULE%ATOM(I)%CENTER(3))
ENDDO

END SUBROUTINE DETERMINE_MAXCOOR

!> \brief print bare molecule
!> \author T. Kjaergaard
!> \date  2010-02-24
SUBROUTINE print_mol(molecule,lupri)
TYPE(MOLECULEINFO) :: MOLECULE
INTEGER :: LUPRI
!
INTEGER :: I 
WRITE(LUPRI,*) '                     '
WRITE(LUPRI,'(A)')'THE BARE MOLECULE:',MOLECULE%label
WRITE(LUPRI,*) '--------------------------------------------------------------------'
WRITE(LUPRI,'(A38,2X,I7)')'Molecular nAtoms                 :',MOLECULE%nAtoms
WRITE(LUPRI,'(A38,2X,I7)')'Molecular nAtomsNPC              :',MOLECULE%nAtomsNPC
WRITE(LUPRI,'(A38,2X,I7)')'Molecular nelec                  :',MOLECULE%nelectrons
WRITE(LUPRI,'(A38,2X,F8.4)')'Molecular Charge                 :',MOLECULE%charge
WRITE(LUPRI,'(2X,A38,2X,I7)')'Regular basisfunctions             :',MOLECULE%nbastREG
WRITE(LUPRI,'(2X,A38,2X,I7)')'Auxiliary basisfunctions           :',MOLECULE%nbastAUX
WRITE(LUPRI,'(2X,A38,2X,I7)')'CABS basisfunctions                :',MOLECULE%nbastCABS
WRITE(LUPRI,'(2X,A38,2X,I7)')'JK-fit basisfunctions              :',MOLECULE%nbastJK
WRITE(LUPRI,'(2X,A38,2X,I7)')'ADMM basisfunctions                :',MOLECULE%nbastADMM
WRITE(LUPRI,'(2X,A38,2X,I7)')'Valence basisfunctions             :',MOLECULE%nbastVAL
WRITE(LUPRI,'(2X,A38,2X,I7)')'Primitive Regular basisfunctions   :',MOLECULE%nprimbastREG
WRITE(LUPRI,'(2X,A38,2X,I7)')'Primitive Auxiliary basisfunctions :',MOLECULE%nprimbastAUX
WRITE(LUPRI,'(2X,A38,2X,I7)')'Primitive CABS basisfunctions      :',MOLECULE%nprimbastCABS
WRITE(LUPRI,'(2X,A38,2X,I7)')'Primitive JK-fit basisfunctions    :',MOLECULE%nprimbastJK
WRITE(LUPRI,'(2X,A38,2X,I7)')'Primitive ADMM basisfunctions      :',MOLECULE%nprimbastADMM
WRITE(LUPRI,'(2X,A38,2X,I7)')'Primitive Valence basisfunctions   :',MOLECULE%nprimbastVAL
IF(MOLECULE%nSubSystems.NE.0)THEN
   WRITE(LUPRI,'(2X,A38,2X,I7)')'number of Subsystem Labels         :',MOLECULE%nSubSystems
   DO I=1,MOLECULE%nSubSystems
      WRITE(LUPRI,'(2X,A,I3,A,A80)')'Label(',I,'):',MOLECULE%SubsystemLabel(I)
   ENDDO
ENDIF
WRITE(LUPRI,*) '--------------------------------------------------------------------'
WRITE(LUPRI,*) '                     '

DO I=1,MOLECULE%nAtoms
   IF(MOLECULE%ATOM(I)%pointcharge)THEN
      WRITE(LUPRI,'(2X,I4,2X,A4,2X,I7,2X,F16.8,2X,F16.8,2X,F16.8)') I,&
           & MOLECULE%ATOM(I)%Name,&
           & MOLECULE%ATOM(I)%Isotope,&
           & MOLECULE%ATOM(I)%CENTER(1),&
           & MOLECULE%ATOM(I)%CENTER(2),&
           & MOLECULE%ATOM(I)%CENTER(3)
   ELSE
      WRITE(LUPRI,'(2X,I4,2X,A4,2X,I7,2X,F16.8,2X,F16.8,2X,F16.8,A12)') I,&
           & MOLECULE%ATOM(I)%Name,&
           & MOLECULE%ATOM(I)%Isotope,&
           & MOLECULE%ATOM(I)%CENTER(1),&
           & MOLECULE%ATOM(I)%CENTER(2),&
           & MOLECULE%ATOM(I)%CENTER(3),' pointcharge'
   ENDIF
ENDDO
WRITE(LUPRI,*) '                     '

END SUBROUTINE PRINT_MOL

!> \brief Initializes the molecular info
!> \author S. Host
!> \date 2010-02-21
!> \param molecule Contains the information about the molecule
!> \param nAtoms The number of atoms
!> \param label  Label to uniqely specify the molecule (use for instance for screening matrix storage)
SUBROUTINE init_Moleculeinfo(Molecule,nAtoms,label)
implicit none
TYPE(MOLECULEINFO),intent(inout) :: MOLECULE
INTEGER,intent(IN)               :: nAtoms
Character(len=22),intent(IN)     :: label

NULLIFY(MOLECULE%ATOM)
MOLECULE%nAtoms = nAtoms
MOLECULE%nAtomsNPC = nAtoms
MOLECULE%nelectrons = 0
IF (nAtoms.GT. 0) call mem_alloc(MOLECULE%ATOM,nAtoms)
call nullifyAtoms(MOLECULE%ATOM)
MOLECULE%label = label

MOLECULE%nbastREG = 0
MOLECULE%nbastAUX = 0
MOLECULE%nbastCABS = 0
MOLECULE%nbastJK = 0
MOLECULE%nbastADMM = 0
MOLECULE%nbastVAL = 0
MOLECULE%nprimbastREG = 0
MOLECULE%nprimbastAUX = 0
MOLECULE%nprimbastCABS = 0
MOLECULE%nprimbastJK = 0
MOLECULE%nprimbastADMM = 0
MOLECULE%nprimbastVAL = 0
MOLECULE%pointmolecule = .false.
MOLECULE%nSubSystems = 0
NULLIFY(MOLECULE%SubsystemLabel)
END SUBROUTINE init_Moleculeinfo

subroutine nullifyAtoms(ATOMS)
implicit none
type(atomitem) :: ATOMS(:)
!
integer :: n,i
n=size(ATOMS)
do i=1,n
   ATOMS(I)%Isotope = 0 
   ATOMS(I)%Name = '    '  
   ATOMS(I)%Mass = 0.0E0_realk 
   ATOMS(I)%CovRad = 0.0E0_realk 
   ATOMS(I)%Frag  = 0.0E0_realk  
   ATOMS(I)%CENTER(1) = 0.0E0_realk 
   ATOMS(I)%CENTER(2) = 0.0E0_realk 
   ATOMS(I)%CENTER(3) = 0.0E0_realk 
   ATOMS(I)%Atomic_number = 0 
   ATOMS(I)%Charge = 0.0E0_realk        
   ATOMS(I)%basislabel = '         '
   ATOMS(I)%Basisindex = 0  
   ATOMS(I)%IDtype = 0
   ATOMS(I)%Phantom = .FALSE.
   ATOMS(I)%Pointcharge = .FALSE.
   ATOMS(I)%nContOrbREG = 0
   ATOMS(I)%nPrimOrbREG = 0
   ATOMS(I)%nContOrbAUX = 0
   ATOMS(I)%nPrimOrbAUX = 0
   ATOMS(I)%nContOrbCABS = 0
   ATOMS(I)%nPrimOrbCABS = 0
   ATOMS(I)%nContOrbJK = 0
   ATOMS(I)%nPrimOrbJK = 0
   ATOMS(I)%nContOrbADMM = 0
   ATOMS(I)%nPrimOrbADMM = 0
   ATOMS(I)%nContOrbVAL = 0
   ATOMS(I)%nPrimOrbVAL = 0
   ATOMS(I)%molecularIndex = 0
   ATOMS(I)%SubSystemIndex = 0
enddo

end subroutine nullifyAtoms

!> \brief Frees the molecular info
!> \author S. Host
!> \date 2010-02-21
!> \param molecule Contains the information about the molecule
SUBROUTINE free_Moleculeinfo(Molecule)
implicit none
TYPE(MOLECULEINFO) :: MOLECULE

if (MOLECULE%nAtoms.GT. 0) then
  if (.not.ASSOCIATED(MOLECULE%ATOM)) then
    print*,'memory previously released!!'
    print*,'MOLECULE%nAtoms',MOLECULE%nAtoms
    call lsquit('Error in Molecule_free - memory previously released',-1)
  endif
  call mem_dealloc(MOLECULE%ATOM)
endif
IF(Molecule%nSubSystems.NE.0)THEN
   if (.not.ASSOCIATED(Molecule%SubSystemLabel)) then
      print*,'memory previously released!!'
      print*,'Molecule%nSubSystems',Molecule%nSubSystems
      call lsquit('Error in Molecule_free -SubSystemLabel memory previously released',-1)
   endif   
   call mem_dealloc(Molecule%SubSystemLabel)
ENDIF

END SUBROUTINE free_Moleculeinfo

!> \brief build moleculeinfo containing only 1 atom from full molecule
!> \author T. Kjaergaard
!> \date 2010-03-23
!> \param Molecule Contains the information about the original molecule
!> \param atomicMolecule the moleculeinfo to be built
!> \param iatom atom in full molecule
SUBROUTINE build_atomicmolecule(molecule,atomicmolecule,iatom,lupri)
IMPLICIT NONE
TYPE(MOLECULEINFO),intent(INOUT) :: atomicmolecule
TYPE(MOLECULEINFO),intent(IN)  :: MOLECULE
integer,intent(in) :: iatom,lupri
character(len=22) :: label

call mem_alloc(atomicmolecule%ATOM,1)
call nullifyAtoms(atomicmolecule%ATOM)
atomicmolecule%nAtoms=1
atomicmolecule%nAtomsNPC=1
write(label,'(A6,I16)') 'GCATOM',iatom
atomicmolecule%label=label

atomicmolecule%nbastREG     = 0
atomicmolecule%nPrimbastREG = 0
atomicmolecule%nbastAUX     = 0
atomicmolecule%nPrimbastAUX = 0
atomicmolecule%nbastCABS     = 0
atomicmolecule%nPrimbastCABS = 0
atomicmolecule%nbastJK     = 0
atomicmolecule%nPrimbastJK = 0
atomicmolecule%nbastADMM     = 0
atomicmolecule%nPrimbastADMM = 0
atomicmolecule%nbastVAL     = 0
atomicmolecule%nPrimbastVAL = 0

call copy_atom(MOLECULE,Iatom,atomicmolecule,1,lupri)
atomicmolecule%nelectrons=INT(atomicmolecule%ATOM(1)%Charge)
atomicmolecule%charge=0
atomicmolecule%pointmolecule = .false.
atomicmolecule%nSubSystems = 0
NULLIFY(atomicmolecule%SubsystemLabel)

END SUBROUTINE build_atomicmolecule

!> \brief build an 'empty'  moleculeinfo containing only 1 atom
!> \author S. Reine
!> \date 2012-02-07
!> \param Molecule Contains the information about the original molecule
SUBROUTINE build_empty_molecule(molecule)
IMPLICIT NONE
TYPE(MOLECULEINFO),intent(INOUT)  :: MOLECULE
!
character(len=22) :: label
integer :: i

call mem_alloc(molecule%ATOM,1)
call nullifyAtoms(molecule%ATOM)
molecule%nAtoms=1
molecule%nAtomsNPC=1
write(label,'(A5,17X)') 'Empty'
molecule%label=label

molecule%nbastREG     = 0
molecule%nPrimbastREG = 0
molecule%nbastAUX     = 0
molecule%nPrimbastAUX = 0
molecule%nbastCABS     = 0
molecule%nPrimbastCABS = 0
molecule%nbastJK     = 0
molecule%nPrimbastJK = 0
molecule%nbastADMM     = 0
molecule%nPrimbastADMM = 0
molecule%nbastVAL     = 0
molecule%nPrimbastVAL = 0
molecule%nelectrons=0
molecule%charge=0
molecule%pointmolecule = .false.
molecule%nSubSystems = 0
NULLIFY(molecule%SubsystemLabel)

molecule%ATOM(1)%Isotope = 0
molecule%ATOM(1)%Name='Empt'
molecule%ATOM(1)%Mass=0_realk
molecule%ATOM(1)%CovRad=0_realk
molecule%ATOM(1)%Frag=0_realk
molecule%ATOM(1)%CENTER(1)=0_realk
molecule%ATOM(1)%CENTER(2)=0_realk
molecule%ATOM(1)%CENTER(3)=0_realk
molecule%ATOM(1)%Atomic_number=0
molecule%ATOM(1)%Charge=0_realk
do i=1,nBasisBasParam
   molecule%ATOM(1)%basislabel(i)='None'
   molecule%ATOM(1)%Basisindex(i)=0
   molecule%ATOM(1)%IDtype(i)=0
enddo
molecule%ATOM(1)%phantom=.FALSE.
molecule%ATOM(1)%Pointcharge=.FALSE.
molecule%ATOM(1)%molecularIndex=0
molecule%ATOM(1)%SubSystemIndex=-1
molecule%ATOM(1)%nContOrbREG=0
molecule%ATOM(1)%nPrimOrbREG=0
molecule%ATOM(1)%nContOrbAUX=0
molecule%ATOM(1)%nPrimOrbAUX=0
molecule%ATOM(1)%nContOrbCABS=0
molecule%ATOM(1)%nPrimOrbCABS=0
molecule%ATOM(1)%nContOrbJK=0
molecule%ATOM(1)%nPrimOrbJK=0
molecule%ATOM(1)%nContOrbADMM=0
molecule%ATOM(1)%nPrimOrbADMM=0
molecule%ATOM(1)%nContOrbVAL=0
molecule%ATOM(1)%nPrimOrbVAL=0

END SUBROUTINE build_empty_molecule

!> \brief build moleculeinfo containing N points
!> \author T. Kjaergaard
!> \date 2010-03-23
!> \param pointMolecule Contains the information about constructed molecule
!> \param R the array of point coordinates
!> \param N the number of points
!> \param lupri logical unit number of output file
SUBROUTINE build_pointmolecule(pointmolecule,R,N,lupri,charge)
IMPLICIT NONE
TYPE(MOLECULEINFO),intent(INOUT) :: pointmolecule
integer,intent(in) :: lupri,N
real(realk) :: R(3,N)
real(realk),optional :: charge(N)
!
Integer :: I,J
Logical :: unitCharge

unitCharge = .NOT.present(charge)

call mem_alloc(pointmolecule%ATOM,N)
call nullifyAtoms(pointmolecule%ATOM)
pointmolecule%nAtoms=N
pointmolecule%nAtomsNPC=0
pointmolecule%label='Point                 '
pointmolecule%nbastREG     = 0
pointmolecule%nPrimbastREG = 0
pointmolecule%nbastAUX     = 0
pointmolecule%nPrimbastAUX = 0
pointmolecule%nbastCABS     = 0
pointmolecule%nPrimbastCABS = 0
pointmolecule%nbastJK     = 0
pointmolecule%nPrimbastJK = 0
pointmolecule%nbastADMM     = 0
pointmolecule%nPrimbastADMM = 0
pointmolecule%nbastVAL     = 0
pointmolecule%nPrimbastVAL = 0
pointmolecule%nelectrons=0
pointmolecule%charge=0
pointmolecule%pointmolecule = .true.
pointmolecule%nSubSystems = 0
NULLIFY(pointmolecule%SubsystemLabel)

!the Points
DO I=1,N
   pointmolecule%ATOM(I)%Isotope = 0
   pointmolecule%ATOM(I)%Name='XXXX'
   pointmolecule%ATOM(I)%Mass=0.d0   
   pointmolecule%ATOM(I)%CovRad=0.d0   
   pointmolecule%ATOM(I)%Frag=0.d0   
   pointmolecule%ATOM(I)%CENTER(1)=R(1,I)
   pointmolecule%ATOM(I)%CENTER(2)=R(2,I)
   pointmolecule%ATOM(I)%CENTER(3)=R(3,I)
   pointmolecule%ATOM(I)%Atomic_number=0
   IF (unitCharge) THEN
     pointmolecule%ATOM(I)%Charge=-1.d0
   ELSE
     pointmolecule%ATOM(I)%Charge=-Charge(I)
   ENDIF
   do j=1,nBasisBasParam
      pointmolecule%ATOM(I)%basislabel(j)='XXXXXXXXX'
      pointmolecule%ATOM(I)%Basisindex(j)=0
      pointmolecule%ATOM(I)%IDtype(j)=0
   enddo
   pointmolecule%ATOM(I)%phantom=.FALSE.
   pointmolecule%ATOM(I)%Pointcharge=.FALSE.
   pointmolecule%ATOM(I)%molecularIndex=1
   pointmolecule%ATOM(I)%SubSystemIndex=-1
   pointmolecule%ATOM(I)%nContOrbREG=0
   pointmolecule%ATOM(I)%nPrimOrbREG=0
   pointmolecule%ATOM(I)%nContOrbAUX=0
   pointmolecule%ATOM(I)%nPrimOrbAUX=0
   pointmolecule%ATOM(I)%nContOrbCABS=0
   pointmolecule%ATOM(I)%nPrimOrbCABS=0
   pointmolecule%ATOM(I)%nContOrbJK=0
   pointmolecule%ATOM(I)%nPrimOrbJK=0
   pointmolecule%ATOM(I)%nContOrbADMM=0
   pointmolecule%ATOM(I)%nPrimOrbADMM=0
   pointmolecule%ATOM(I)%nContOrbVAL=0
   pointmolecule%ATOM(I)%nPrimOrbVAL=0
enddo

END SUBROUTINE build_pointmolecule

!> \brief build moleculeinfo containing N points
!> \author T. Kjaergaard
!> \date 2010-03-23
!> \param pointMolecule Contains the information about constructed molecule
!> \param R the array of point coordinates
!> \param N the number of points
!> \param lupri logical unit number of output file
SUBROUTINE build_molecule_from_coordlist(molecule,coord,N,Charge,lupri)
IMPLICIT NONE
TYPE(MOLECULEINFO),intent(INOUT) :: molecule
integer,intent(in) :: lupri,N
real(realk) :: coord(3,N),charge(N)
!
Integer :: I,j
character(len=22) :: string22
character(len=9) :: stringA9,stringB9
character(len=4) :: string4
string22='ListMol               '
stringA9 = 'REGULAR  '
stringB9 = 'NOT SET  '
string4 = 'XXXX'
call init_Moleculeinfo(Molecule,N,string22)
molecule%nelectrons=0
molecule%charge=0
!the Points
DO I=1,N
   molecule%ATOM(I)%Isotope = 0
   molecule%ATOM(I)%Name=string4
   molecule%ATOM(I)%Mass=0.d0   
   molecule%ATOM(I)%CovRad=0.d0   
   molecule%ATOM(I)%Frag=0.d0   
   molecule%ATOM(I)%CENTER(1)=coord(1,I)
   molecule%ATOM(I)%CENTER(2)=coord(2,I)
   molecule%ATOM(I)%CENTER(3)=coord(3,I)
   molecule%ATOM(I)%Atomic_number=0
   molecule%ATOM(I)%Charge=charge(I)
   do j=1,nBasisBasParam
      molecule%ATOM(I)%basislabel(j)=stringB9
      molecule%ATOM(I)%Basisindex(j)=1
      molecule%ATOM(I)%IDtype(j)=0
   enddo
   molecule%ATOM(I)%basislabel(RegBasParam)=stringA9
   molecule%ATOM(I)%Basisindex(RegBasParam)=1
   molecule%ATOM(I)%IDtype(RegBasParam)=0
   molecule%ATOM(I)%phantom=.FALSE.
   molecule%ATOM(I)%Pointcharge=.FALSE.
   molecule%ATOM(I)%molecularIndex=I
   molecule%ATOM(I)%SubSystemIndex=-1
   molecule%ATOM(I)%nContOrbREG=0
   molecule%ATOM(I)%nPrimOrbREG=0
   molecule%ATOM(I)%nContOrbAUX=0
   molecule%ATOM(I)%nPrimOrbAUX=0
   molecule%ATOM(I)%nContOrbCABS=0
   molecule%ATOM(I)%nPrimOrbCABS=0
   molecule%ATOM(I)%nContOrbJK=0
   molecule%ATOM(I)%nPrimOrbJK=0
   molecule%ATOM(I)%nContOrbADMM=0
   molecule%ATOM(I)%nPrimOrbADMM=0
   molecule%ATOM(I)%nContOrbVAL=0
   molecule%ATOM(I)%nPrimOrbVAL=0
enddo

END SUBROUTINE build_molecule_from_coordlist

!> \brief Copies a oldMOLECULE to newMOLECULE
!> \author T. Kjaergaard
!> \date 2010-05-31
!> \param oldMolecule Contains the information about the original molecule
!> \param newMolecule new molecule - same as oldMolecule afterwards
!> \param LUPRI LOGICAL UNIT NUMBER FOR OUTPUT
subroutine copy_molecule(oldMolecule,newMolecule,lupri)
implicit none
TYPE(MOLECULEINFO),intent(INOUT) :: newMOLECULE
TYPE(MOLECULEINFO),intent(IN)    :: oldMOLECULE
integer :: lupri
!
integer :: I

newMOLECULE%nAtoms = oldMOLECULE%nAtoms
newMOLECULE%nAtomsNPC = oldMOLECULE%nAtomsNPC
newMOLECULE%label = oldMOLECULE%label
newMOLECULE%nbastREG = 0
newMOLECULE%nbastAUX = 0
newMOLECULE%nbastCABS = 0
newMOLECULE%nbastJK = 0
newMOLECULE%nbastADMM = 0
newMOLECULE%nbastVAL = 0
newMOLECULE%nprimbastREG = 0
newMOLECULE%nprimbastAUX = 0
newMOLECULE%nprimbastCABS = 0
newMOLECULE%nprimbastJK = 0
newMOLECULE%nprimbastADMM = 0
newMOLECULE%nprimbastVAL = 0

call mem_alloc(newMOLECULE%ATOM,newMOLECULE%nAtoms)
call nullifyAtoms(newmolecule%ATOM)
do I = 1, newMOLECULE%nAtoms
   call copy_atom(oldMOLECULE,I,newMOLECULE,I,LUPRI)
enddo
newMOLECULE%nelectrons = oldMOLECULE%nelectrons
newMOLECULE%charge = oldMOLECULE%charge
newMolecule%pointmolecule = oldMolecule%pointmolecule
newMolecule%nSubSystems = oldMolecule%nSubSystems
IF(newMolecule%nSubSystems.NE.0)THEN
   call mem_alloc(newMolecule%SubSystemLabel,newMolecule%nSubSystems)
   do I = 1, newMolecule%nSubSystems
      newMolecule%SubSystemLabel(I) = oldMolecule%SubSystemLabel(I)
   enddo
ELSE
   NULLIFY(newMolecule%SubSystemLabel)
ENDIF
end subroutine copy_molecule

!> \brief Copies an atom from MOLECULE to FRAGMENT
!> \author S. Reine
!> \date 2010-02-21
!> \param Molecule Contains the information about the original molecule
!> \param I Atomic number of atom to be copied
!> \param Fragment Contains the information about the molecule copy
!> \param J Atomic number to be given for the copied atom
!> \param LUPRI LOGICAL UNIT NUMBER FOR OUTPUT
SUBROUTINE copy_atom(MOLECULE,I,FRAGMENT,J,LUPRI)
IMPLICIT NONE
INTEGER,intent(IN)               :: I,J,LUPRI
TYPE(MOLECULEINFO),intent(INOUT) :: FRAGMENT
TYPE(MOLECULEINFO),intent(IN)    :: MOLECULE
integer :: K

IF ((J.GT.FRAGMENT%nAtoms).OR.(J.LT. 0)) THEN
  WRITE(LUPRI,'(1X,A,I4,A,I4)') 'Error in COPY_ATOM. Fragment atomic number ',J,&
     &                          ' is wrong - natoms = ',FRAGMENT%nAtoms
  CALL LSQUIT('Input in-consistency in COPY_ATOM',lupri)
ELSEIF ((I.GT.MOLECULE%nAtoms).OR.(I.LT. 0)) THEN
  WRITE(LUPRI,'(1X,A,I4,A,I4)') 'Error in COPY_ATOM. Molecular atomic number ',I,&
     &                          ' is wrong - natoms = ',MOLECULE%nAtoms
  CALL LSQUIT('Input in-consistency in COPY_ATOM',lupri)
ENDIF

FRAGMENT%ATOM(J)%Isotope=MOLECULE%ATOM(I)%Isotope
FRAGMENT%ATOM(J)%Name=MOLECULE%ATOM(I)%Name
FRAGMENT%ATOM(J)%Mass=MOLECULE%ATOM(I)%Mass
FRAGMENT%ATOM(J)%CovRad=MOLECULE%ATOM(I)%CovRad
FRAGMENT%ATOM(J)%Frag=MOLECULE%ATOM(I)%Frag
FRAGMENT%ATOM(J)%CENTER(1)=MOLECULE%ATOM(I)%CENTER(1)
FRAGMENT%ATOM(J)%CENTER(2)=MOLECULE%ATOM(I)%CENTER(2)
FRAGMENT%ATOM(J)%CENTER(3)=MOLECULE%ATOM(I)%CENTER(3)
FRAGMENT%ATOM(J)%Atomic_number=MOLECULE%ATOM(I)%Atomic_number
FRAGMENT%ATOM(J)%Charge=MOLECULE%ATOM(I)%Charge
do K = 1,nBasisBasParam
   FRAGMENT%ATOM(J)%basislabel(K)=MOLECULE%ATOM(I)%basislabel(K)
   FRAGMENT%ATOM(J)%Basisindex(K)=MOLECULE%ATOM(I)%Basisindex(K)
   FRAGMENT%ATOM(J)%IDtype(K)=MOLECULE%ATOM(I)%IDtype(K)
enddo
FRAGMENT%ATOM(J)%phantom=MOLECULE%ATOM(I)%phantom
FRAGMENT%ATOM(J)%PointCharge=MOLECULE%ATOM(I)%PointCharge
FRAGMENT%ATOM(J)%molecularIndex=MOLECULE%ATOM(I)%molecularIndex
FRAGMENT%ATOM(J)%SubSystemIndex=MOLECULE%ATOM(I)%SubSystemIndex
FRAGMENT%ATOM(J)%nContOrbREG=MOLECULE%ATOM(I)%nContOrbREG
FRAGMENT%ATOM(J)%nPrimOrbREG=MOLECULE%ATOM(I)%nPrimOrbREG
FRAGMENT%ATOM(J)%nContOrbAUX=MOLECULE%ATOM(I)%nContOrbAUX
FRAGMENT%ATOM(J)%nPrimOrbAUX=MOLECULE%ATOM(I)%nPrimOrbAUX
FRAGMENT%ATOM(J)%nContOrbCABS=MOLECULE%ATOM(I)%nContOrbCABS
FRAGMENT%ATOM(J)%nPrimOrbCABS=MOLECULE%ATOM(I)%nPrimOrbCABS
FRAGMENT%ATOM(J)%nContOrbJK=MOLECULE%ATOM(I)%nContOrbJK
FRAGMENT%ATOM(J)%nPrimOrbJK=MOLECULE%ATOM(I)%nPrimOrbJK
FRAGMENT%ATOM(J)%nContOrbADMM=MOLECULE%ATOM(I)%nContOrbADMM
FRAGMENT%ATOM(J)%nPrimOrbADMM=MOLECULE%ATOM(I)%nPrimOrbADMM
FRAGMENT%ATOM(J)%nContOrbVAL=MOLECULE%ATOM(I)%nContOrbVAL
FRAGMENT%ATOM(J)%nPrimOrbVAL=MOLECULE%ATOM(I)%nPrimOrbVAL

FRAGMENT%nbastREG     = FRAGMENT%nbastREG     + MOLECULE%ATOM(I)%nContOrbREG
FRAGMENT%nPrimbastREG = FRAGMENT%nPrimbastREG + MOLECULE%ATOM(I)%nPrimOrbREG
FRAGMENT%nbastAUX     = FRAGMENT%nbastAUX     + MOLECULE%ATOM(I)%nContOrbAUX
FRAGMENT%nPrimbastAUX = FRAGMENT%nPrimbastAUX + MOLECULE%ATOM(I)%nPrimOrbAUX
FRAGMENT%nbastCABS     = FRAGMENT%nbastCABS     + MOLECULE%ATOM(I)%nContOrbCABS
FRAGMENT%nPrimbastCABS = FRAGMENT%nPrimbastCABS + MOLECULE%ATOM(I)%nPrimOrbCABS
FRAGMENT%nbastJK     = FRAGMENT%nbastJK     + MOLECULE%ATOM(I)%nContOrbJK
FRAGMENT%nPrimbastJK = FRAGMENT%nPrimbastJK + MOLECULE%ATOM(I)%nPrimOrbJK
FRAGMENT%nbastADMM     = FRAGMENT%nbastADMM     + MOLECULE%ATOM(I)%nContOrbADMM
FRAGMENT%nPrimbastADMM = FRAGMENT%nPrimbastADMM + MOLECULE%ATOM(I)%nPrimOrbADMM
FRAGMENT%nbastVAL     = FRAGMENT%nbastVAL     + MOLECULE%ATOM(I)%nContOrbVAL
FRAGMENT%nPrimbastVAL = FRAGMENT%nPrimbastVAL + MOLECULE%ATOM(I)%nPrimOrbVAL

END SUBROUTINE copy_atom

!> \brief gets the number of atoms
!> \author Bin Gao
!> \date 2012-05-24
!> \param mol_info contains the information of molecule
!> \return num_atoms is the number of atoms
subroutine molecule_get_num_atoms(mol_info, num_atoms)
  implicit none
  type(MOLECULEINFO), intent(in) :: mol_info
  integer, intent(out) :: num_atoms
  num_atoms = mol_info%nAtoms
  return
end subroutine molecule_get_num_atoms

!> \brief gets the coordinates of atoms
!> \author Bin Gao
!> \date 2012-05-24
!> \param mol_info contains the information of molecule
!> \param num_atoms is the number of atoms
!> \param lupri is the logical unit number of the standard output
!> \return coord_atoms contains the coordinates of atoms
subroutine molecule_get_atom_coord(mol_info, num_atoms, coord_atoms, lupri)
  implicit none
  type(MOLECULEINFO), intent(in) :: mol_info
  integer, intent(in) :: num_atoms
  real(realk), intent(out) :: coord_atoms(3,num_atoms)
  integer, intent(in) :: lupri
  integer iatom  !incremental recorder over atoms
  if (num_atoms/=mol_info%nAtoms) then
    call lsquit("molecule_get_atom_coord>> incorrect number of atoms", lupri)
  else
    do iatom = 1, mol_info%nAtoms
      coord_atoms(:,iatom) = mol_info%ATOM(iatom)%CENTER
    end do
  end if
  return
end subroutine molecule_get_atom_coord

!> \brief gets the charges of atoms
!> \author Bin Gao
!> \date 2012-05-24
!> \param mol_info contains the information of molecule
!> \param num_atoms is the number of atoms
!> \param lupri is the logical unit number of the standard output
!> \return charge_atoms contains the charges of atoms
subroutine molecule_get_atom_charge(mol_info, num_atoms, charge_atoms, lupri)
  implicit none
  type(MOLECULEINFO), intent(in) :: mol_info
  integer, intent(in) :: num_atoms
  real(realk), intent(out) :: charge_atoms(num_atoms)
  integer, intent(in) :: lupri
  integer iatom  !incremental recorder over atoms
  if (num_atoms/=mol_info%nAtoms) then
    call lsquit("molecule_get_atom_charge>> incorrect number of atoms", lupri)
  else
    do iatom = 1, mol_info%nAtoms
      charge_atoms(iatom) = -mol_info%ATOM(iatom)%Charge
    end do
  end if
  return
end subroutine molecule_get_atom_charge

END MODULE molecule_type
