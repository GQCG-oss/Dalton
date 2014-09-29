!> @file 
!> Contains AtomSparse module

!> \brief Contains structure and routines for a pair-wise sparse matrix structure
!> \author S. Reine and P. Merlot
!> \date 2010-07-20
!> (Developed for the PARI scheme - but can be extended for other purposes)
MODULE AtomSparse
use molecule_type
use molecule_typetype
use memory_handling
use precision

!> \brief Contains the pair-atomic matrix subblock (i.e. the matrix elements between 
!>        the basis-functions centered on the two atoms A and B)
!> \author S. Reine and P. Merlot
!> \date 2010-07-20
!> \param atomB      Atomic index of atom B
!> \param nB         Number of basis functions on atom B
!> \param elms       The matrix elements corresponding to the AB-block of the matrix
!> \param calculated Specifying whether the elements have been calculated (or only allocated)
!> \param next       Pointer to the next ABblock in the ABlist
TYPE ABblock
integer                 :: atomB
integer                 :: nB
integer                 :: nMat
real(realk), pointer    :: elms(:,:,:) !nA,nB,nMat
logical                 :: calculated
type(ABblock), pointer  :: next
END TYPE ABblock

!> \brief For given atom A this list contains the ABblocks (for list of atoms B)
!> \author S. Reine and P. Merlot
!> \date 2010-07-20
!> Note that this is a nested list (where one element points to the next)
!> \param nAtomsB Number of atoms B in the list
!> \param nA      Number of basis functions on atom A
!> \param first   Points to the first ABblock in the list
TYPE ABlist
integer                 :: nAtomsB
integer                 :: nA
type(ABblock), pointer  :: first
END TYPE ABlist

!> \brief The atomic-sparse matrix (storage by pair-atomic matrix subblocks - corresponding 
!>         to the interaction between the basis functions centered on two atoms)
!> \author S. Reine and P. Merlot
!> \date 2010-07-20
!> \param nAtomsA     The number of atoms A
!> \param nAtomsBfull The number of atoms B
!> \param Alist       The list of pairs AB for each atom A
!> \param AO1         Specified the basis employed for the first AO
!> \param AO2         Specified the basis employed for the second AO
!> \param Mol1        The molecule of the first AO
!> \param Mol2        The molecule of the second AO
!> \param primitive   Specifies whether a primitive or contracted basis is used
!> \param allocated   Specifies if the Alist have been allocated
TYPE AtomSparseMat
integer                    :: nAtomsA
integer                    :: nAtomsBfull
type(ABlist), pointer      :: Alist(:)
Character*(80)             :: AO1
Character*(80)             :: AO2
TYPE(moleculeInfo),pointer :: Mol1
TYPE(moleculeInfo),pointer :: Mol2
logical                    :: primitive
logical                    :: sameDims
logical                    :: allocated
END TYPE AtomSparseMat

CONTAINS

!> \brief This subroutine allocates memory of a list of atom pairs containing (alpha|beta)
!> \author P. Merlot, S. Reine
!> \date 2010-07-20
!> \param ABsparseMat The atomic-sparse matrix
!> \param Mol1        The molecule of the first AO
!> \param Mol2        The molecule of the second AO
!> \param AO1         Specification of the first AO ('Regular','DF-Aux', etc.)
!> \param AO2         Specification of the second AO 
!> \param intType     Specifies whether we use contracted or primitive basis
!> \param lupri       Default print unit
SUBROUTINE init_AtomSparseMat(ABsparseMat,Mol1,Mol2,AO1,AO2,intType,lupri)
implicit none
type(AtomSparseMat),intent(inout)    :: ABsparseMat
integer,intent(in)                   :: lupri
type(moleculeinfo),intent(in),target :: Mol1
type(moleculeinfo),intent(in),target :: Mol2
character*(*),intent(in)             :: AO1,AO2,intType
!
Integer :: nAtoms1,iAtom,nBastA
!
nAtoms1 = Mol1%nAtoms
!
ABsparseMat%nAtomsA     = nAtoms1
ABsparseMat%nAtomsBfull = Mol2%nAtoms
nullify(ABsparseMat%Alist)
allocate(ABsparseMat%Alist(nAtoms1))

! Select primitive or contracted basis
SELECT CASE (intType)
CASE ('Primitive')
  ABsparseMat%primitive = .TRUE.
CASE ('Contracted')
  ABsparseMat%primitive = .FALSE.
CASE DEFAULT
  write(lupri,'(1X,2A)') 'Error in init_AtomSparseMat. intType not idenified:',intType
  CALL LSQUIT('intType not idenfified in init_AtomSparseMat',lupri)
END SELECT

DO iAtom=1,nAtoms1
!   Select AO1 type
  SELECT CASE(AO1)
  CASE('Regular')
   nBastA = Mol1%Atom(iAtom)%nContOrbREG
   IF (ABsparseMat%primitive) nBastA = Mol1%Atom(iAtom)%nPrimOrbREG
  CASE('DF-Aux')
   nBastA = Mol1%Atom(iAtom)%nContOrbAUX
   IF (ABsparseMat%primitive) nBastA = Mol1%Atom(iAtom)%nPrimOrbAUX
  CASE DEFAULT
    write(lupri,'(1X,2A)') 'Error in init_AtomSparseMat. AO1 not identified:',AO1
    CALL LSQUIT('AO1 not idenfified in init_AtomSparseMat',lupri)
  END SELECT
  call init_ablist(ABsparseMat%Alist(iAtom),nBastA,lupri)
ENDDO
ABsparseMat%allocated = .TRUE.

ABsparseMat%AO1     = AO1
ABsparseMat%AO2     = AO2
ABsparseMat%MOL1    => Mol1
ABsparseMat%MOL2    => Mol2

END SUBROUTINE init_AtomSparseMat

!> \brief This subroutine deallocates memory of a list of atom pairs containing (alpha|beta)
!> \author P. Merlot, S. Reine
!> \date 2010-07-20
!> \param ABsparseMat The atomic-sparse matrix
 SUBROUTINE free_AtomSparseMat(ABsparseMat)
  implicit none
  type(AtomSparseMat),intent(inout) :: ABsparseMat
  !
  Integer :: iAtom
  !
  IF (.NOT.ABsparseMat%allocated) &
     & CALL LSQUIT('Error in free_AtomSparseMat: ABsparseMat not allocated',-1)
  !
  DO iAtom=1,ABsparseMat%nAtomsA
!    call free_ablist(ABsparseMat%Alist(iAtom))
  ENDDO
  ABsparseMat%nAtomsA     = 0
  ABsparseMat%nAtomsBfull = 0
  deallocate(ABsparseMat%Alist)
  nullify(ABsparseMat%Alist)
  ABsparseMat%allocated = .FALSE.
END SUBROUTINE free_AtomSparseMat


!> \brief Initializes an ABlist item
!> \author P. Merlot, S. Reine
!> \date 2010-07-20
!> \param list The list of atoms B for a given atom A, for which the AB-pair exists
!> \param nBasis The number of basis functions for atom A
!> \param lupri Default print unit
SUBROUTINE init_ablist(list,nBasis,lupri)
implicit none
TYPE(ABlist),intent(inout) :: list
Integer,intent(in)         :: nBasis,lupri
!
IF (nBasis.LE. 0) THEN
  WRITE(LUPRI,'(1X,A,I10)') 'Error in init_ablist. nBasis =',nBasis
  CALL LSQUIT('Error in init_ablist. Zero or negative nBasis',lupri)
ENDIF
!
list%nAtomsB = 0
list%nA      = nBasis
nullify(list%first)
!
END SUBROUTINE init_ablist

!> \brief Frees an ABlist item
!> \author P. Merlot, S. Reine
!> \date 2010-07-20
!> \param list The list of atoms B for a given atom A, for which the AB-pair exists
!> \param lupri Default print unit
SUBROUTINE free_ablist(list,lupri)
implicit none
TYPE(ABlist),intent(inout) :: list
Integer,intent(in)         :: lupri
!
Integer :: iAtom
TYPE(ABblock),pointer :: current,next
!
! Free each allocated ABblock
IF (list%nAtomsB.GT. 0) THEN
  current => list%first
  next    => current%next
  call free_abblock(current,lupri)
  DO iAtom=2,list%nAtomsB
    current => next
    next    => current%next
    call free_abblock(current,lupri)
  ENDDO
  IF (associated(next)) THEN
    WRITE(LUPRI,'(1X,A)') 'Error in free_ablist. Last item have a non-empty next'
    CALL LSQUIT('Error in free_ablist. Non-empty next in last item',lupri)
  ENDIF
ENDIF
!
list%nAtomsB = 0
list%nA      = 0
nullify(list%first)
!
END SUBROUTINE free_ablist

!SUBROUTINE add_abblock(block,lupri)
!implicit none
!TYPE(ABblock),intent(inout) :: block
!Integer,intent(in)          :: lupri
!!
!block%iAtomB     = 0
!block%nB         = 0
!block%calculated = .FALSE.
!nullify(block%next)
!!
!IF (.NOT.associated(block%elms)) THEN
!  WRITE(LUPRI,'(1X,A)') 'Error in free_abblock. Elements not associated'
!  CALL LSQUIT('Error in free_abblock. Elements not associated',lupri)
!ENDIF
!call mem_dealloc(block%elms)
!!
!END SUBROUTINE add_abblock

!> \brief Initializes an abblock item (for a pair of atoms AB)
!> \author S. Reine and P. merlot
!> \date 2010-07-20
!> \param pairBlock Abblock
!> \param Atom2 The atom B
!> \param AO2   The orbital type of B ('Regular', 'DF-aux' etc.)
!> \param intType Specifies whether we use primitive or contracted (or other) basis
!> \param iAtom2 The atomic number of atom B
!> \param nextBlock The next block in the (ordered) list (can be empty)
!> \param nA The number of basis functions for atom A
!> \param nMat The number of AB matrices
!> \param lupri Default print unit number
SUBROUTINE init_abblock(pairBlock,Atom2,AO2,primitive,iAtom2,nextBlock,nA,nMat,lupri)
use molecule_type
implicit none
TYPE(ABblock),intent(inout) :: pairBlock
type(atomitem),intent(in)       :: Atom2
character*(*),intent(in)    :: AO2
logical,intent(in)          :: primitive
Integer,intent(in)          :: iAtom2,lupri,nA,nMat
TYPE(ABblock),pointer       :: nextBlock
!
Integer :: nBastB
!
! Select AO2 type
SELECT CASE(AO2)
CASE('Regular')
 nBastB = Atom2%nContOrbREG
 IF (primitive) nBastB = Atom2%nPrimOrbREG
CASE('DF-Aux')
 nBastB = Atom2%nContOrbAUX
 IF (primitive) nBastB = Atom2%nPrimOrbAUX
CASE DEFAULT
  write(lupri,'(1X,2A)') 'Error in init_abblock. AO2 not identified:',AO2
  CALL LSQUIT('AO1 not idenfified in init_abblock',lupri)
END SELECT

pairBlock%atomB      = iAtom2
pairBlock%nB         = nBastB
pairBlock%nMat       = nMat
call mem_alloc(pairBlock%elms,nA,pairBlock%nB,nMat)
pairBlock%calculated = .FALSE.
pairBlock%next       = nextBlock
!
END SUBROUTINE init_abblock

!> \brief Frees an abblock item (for a pair of atoms AB)
!> \author S. Reine and P. merlot
!> \date 2010-07-20
!> \param pairBlock Abblock
!> \param lupri Default print unit number
SUBROUTINE free_abblock(pairBlock,lupri)
implicit none
TYPE(ABblock),intent(inout) :: pairBlock
Integer,intent(in)          :: lupri
!
pairBlock%atomB     = 0
pairBlock%nB         = 0
pairBlock%calculated = .FALSE.
nullify(pairBlock%next)
!
IF (.NOT.associated(pairBlock%elms)) THEN
  WRITE(LUPRI,'(1X,A)') 'Error in free_abblock. Elements not associated'
  CALL LSQUIT('Error in free_abblock. Elements not associated',lupri)
ENDIF
call mem_dealloc(pairBlock%elms)
!
END SUBROUTINE free_abblock

!> \brief Extracts the ABblock (the AB pair-atomic subblock) from the metrix
!> \author S. Reine and P. Merlot
!> \date 2010-07-20
!> \param ABsparseMat Atomic sparse matrix (to extract from)
!> \param atomA       Atomic index of atom A
!> \param atomB       Atomic index of atom B
!> \param nA          Number of basis functions on atom A
!> \param nB          Number of basis functions on atom B
!> \param nMat        The number of AB matrices
!> \param lupri       Default print unit
FUNCTION extract_ABblock(ABsparseMat,atomA,atomB,nA,nB,nMat,lupri)
implicit none
TYPE(AtomSparseMat),intent(in) :: ABsparseMat
Integer,intent(in)             :: atomA,atomB,nA,nB,nMat,lupri
Real(realk)                    :: extract_ABblock(nA,nB,nMat)
!
Type(ABblock),pointer :: pairBlock

pairBlock => get_block(ABsparseMat,atomA,AtomB)
IF (.NOT.associated(pairBlock)) THEN
  WRITE(LUPRI,'(1X,A,2I7)') 'Error in extract_ABblock. Block does not exist for atomA,atomB =',atomA,atomB
  CALL LSQUIT('Error in extract_ABblock. Block does not exist',lupri)
ENDIF

extract_ABblock = pairBlock%elms 

END FUNCTION extract_ABblock

!> \brief Puts an AB block into the matrix
!> \author S. Reine and P. Merlot
!> \date 2010-07-20
!> \param ABsparseMat Atomic sparse matrix
!> \param ABelms      The elements of the subblock AB to put into the ABsparseMat
!> \param atomA       Atomic index of atom A
!> \param atomB       Atomic index of atom B
!> \param nA          Number of basis functions on atom A
!> \param nB          Number of basis functions on atom B
!> \param nMat        The number of AB matrices
!> \param lupri       Default print unit
SUBROUTINE put_ABblock(ABsparseMat,ABelms,atomA,atomB,nA,nB,nMat,lupri)
implicit none
TYPE(AtomSparseMat)            :: ABsparseMat
Integer,intent(in)             :: atomA,atomB,nA,nB,nMat,lupri
Real(realk),intent(in)         :: ABelms(nA,nB,nMat)
!
Type(ABblock),pointer :: beforeBlock,afterBlock,insertBlock
Integer               :: iAtomB

! Consistency testing of nA
IF (nA.NE.ABsparseMat%Alist(atomA)%nA) THEN
  WRITE(LUPRI,'(1X,A,2I7)') &
     &'Error in put_ABblock. Mismatching dimensions nA =',&
     &nA,ABsparseMat%Alist(atomA)%nA
  CALL LSQUIT('Error in put_ABblock. Mismatching dimensions nA',lupri)
ENDIF

! We will not overwrite an existing block
IF (block_exists(ABsparseMat,atomA,AtomB)) THEN
  WRITE(LUPRI,'(1X,A,2I7)') &
     &'Error in put_ABblock. Trying to add an exsisting block for atomA,atomB =',&
     &atomA,atomB
  CALL LSQUIT('Error in put_ABblock. Trying to add an exsisting block',lupri)
ENDIF

! Find the blocks to insert the new block in between
IF (ABsparseMat%Alist(atomA)%nAtomsB.GT. 0) THEN
  beforeBlock => ABsparseMat%Alist(atomA)%first
  afterBlock  => beforeBlock%next
  DO iAtomB=2,ABsparseMat%Alist(atomA)%nAtomsB
    IF (afterBlock%atomB.GT.atomB) EXIT
    beforeBlock => afterBlock
    afterBlock  => beforeBlock%next
  ENDDO
! Special case where the new block is the first element
  IF (beforeBlock%atomB.GT.atomB) THEN
    NULLIFY(ABsparseMat%Alist(atomA)%first)
    ALLOCATE(ABsparseMat%Alist(atomA)%first)
    insertBlock => ABsparseMat%Alist(atomA)%first
! General case
  ELSE
   NULLIFY(beforeBlock%next)
   ALLOCATE(beforeBlock%next)
   insertBlock => beforeBlock%next
  ENDIF
! Special case where we create the first block
ELSE
  NULLIFY(ABsparseMat%Alist(atomA)%first)
  ALLOCATE(ABsparseMat%Alist(atomA)%first)
  insertBlock => ABsparseMat%Alist(atomA)%first
  NULLIFY(afterBlock)
ENDIF

!
CALL init_ABblock(insertBlock,ABsparseMat%MOL2%ATOM(atomB),ABsparseMat%AO2,&
     &            ABsparseMat%primitive,atomB,afterBlock,nA,nMat,lupri)
!
! Consistency testing of nB
IF (insertBlock%nB.NE.nB) THEN
  WRITE(LUPRI,'(1X,A,2I6)') 'Error in put_ABblock. Mismatching dimensions nB =',&
     &nB,insertBlock%nB
  CALL LSQUIT('Error in put_ABblock. Mismatching dimensions nB',lupri)
ENDIF

insertBlock%elms       = ABelms
insertBlock%calculated = .true.

END SUBROUTINE put_ABblock

!> \brief Finds the ABblock corresponding to atomA and atomB if it exists - otherwise it returns a null pointer
!> \author S. Reine and P. Merlot
!> \date 2010-7-20
!> \param ABsparseMat The atomic sparse matrix
!> \param atomA The first atom
!> \param atomB The second atom
FUNCTION get_block(ABsparseMat,atomA,atomB)
implicit none
TYPE(AtomSparseMat),intent(in) :: ABsparseMat
Integer,intent(in)             :: atomA,atomB
TYPE(ABblock),pointer          :: get_block
!
Integer :: iAtom
!
IF (atomA.GT.ABsparseMat%nAtomsA) CALL LSQUIT('Error in get_block. atomA too large',-1)
IF (ABsparseMat%Alist(atomA)%nAtomsB.EQ. 0) CALL LSQUIT('Error in get_block. No atoms B',-1)

get_block => ABsparseMat%Alist(atomA)%first
DO iAtom=2,ABsparseMat%Alist(atomA)%nAtomsB
  IF (get_block%atomB.GE.atomB) THEN
    EXIT
  ENDIF
  get_block => get_block%next
ENDDO
! Check to see if block has been found - otherwise nullify it
IF (.not.(get_block%atomB.EQ.atomB)) THEN
  nullify(get_block)
ENDIF
END FUNCTION get_block

!> \brief Checks if the ABblock corresponding to atomA and atomB exists
!> \author S. Reine and P. Merlot
!> \date 2010-7-20
!> \param ABsparseMat The atomic sparse matrix
!> \param atomA The first atom
!> \param atomB The second atom
LOGICAL FUNCTION block_exists(ABsparseMat,atomA,atomB)
implicit none
TYPE(AtomSparseMat),intent(in) :: ABsparseMat
Integer,intent(in)             :: atomA,atomB
!
TYPE(ABblock),pointer :: pairBlock
!
pairBlock => get_block(ABsparseMat,atomA,atomB)
block_exists = associated(pairBlock)
!
END FUNCTION block_exists

!> \brief Prints information about the atomic-sparse matrix based on the value of iprint
!> \author S. Reine
!> \date 2010-07-21
!> \param ABsparseMat The atomic-sparse matrix
!> \param iunit       The output unit number
!> \param iprint      The print level
!> \param label       Idenfifier to be printed
SUBROUTINE print_atomicSparseMat(ABsparseMat,iunit,iprint,label)
implicit none
TYPE(AtomSparseMat),intent(in) :: ABsparseMat
Integer,intent(in)             :: iunit,iprint
Character*(*),intent(in)       :: label

! Do not print anything if iprint < 2
IF (iprint.GT. 1) THEN
  write(iunit,'(1X,2A)')    'Specification for atomic-sparse matrix ', label
  write(iunit,'(3X,A,I5)')  'Number of atoms A ', ABsparseMat%nAtomsA
  write(iunit,'(3X,A,I5)')  'Number of atoms B ', ABsparseMat%nAtomsBfull
  IF (ABsparseMat%primitive) THEN
    write(iunit,'(3X,A)')  'Primitive basis set used'
  ELSE
    write(iunit,'(3X,A)')  'Contracted basis set used'
  ENDIF   ! primitive
  write(iunit,'(3X,2A)') 'AO1 specification', ABsparseMat%AO1
  write(iunit,'(3X,2A)') 'AO2 specification', ABsparseMat%AO2
  IF (iprint.GT. 10) call mat_print_atomicSparseMat(ABsparseMat,iunit)
ENDIF   ! iprint
END SUBROUTINE print_atomicSparseMat

!> \brief Prints the atomic-sparse matrix elements block-by-block
!> \author S. Reine
!> \date 2010-07-21
!> \param ABsparseMat The atomic-sparse matrix
!> \param iunit       The output unit number
SUBROUTINE mat_print_atomicSparseMat(ABsparseMat,iunit)
implicit none
TYPE(AtomSparseMat),intent(in) :: ABsparseMat
Integer,intent(in)             :: iunit
!
Integer                :: iMat,iAtomA,iAtomB,nA,nB
TYPE(ABblock), pointer :: currentblock
Logical                :: printed

printed = .false.
DO iAtomA=1,ABsparseMat%nAtomsA
  currentBlock => ABsparseMat%Alist(iAtomA)%first
  nA = ABsparseMat%Alist(iAtomA)%nA
  DO iAtomB=1,ABsparseMat%Alist(iAtomA)%nAtomsB
    nB = currentBlock%nB
    WRITE(IUNIT,'(5X,A,I6,A,I6)') 'ABblock for atom A =',iAtomA, ' and atom B =',iAtomB
    DO iMat=1,currentBlock%nMat
      WRITE(IUNIT,'(5X,A,I6,A,I6)') 'Matrix number ',iMat, ' out of a total of',currentBlock%nMat
          CALL ls_output(currentBlock%elms(1,1,iMat),1,nA,1,nB,nA,nB,1,iunit)
    ENDDO
    currentBlock => currentBlock%next
    printed = .true.
  ENDDO
ENDDO

IF (.not. printed) WRITE(IUNIT,'(3X,A)') 'Zero matrix'

END SUBROUTINE mat_print_atomicSparseMat

END MODULE AtomSparse
