!> @file
!> Module containing subroutines related to the molecule

!> Standard molecule module. Contains also the orbital information.
!> \author S. Reine
!> \date 2010-02-21
MODULE molecule_module
use precision
use typedeftype
use typedef
use memory_handling
use molecule_type
use molecule_typetype
use LSparameters
CONTAINS


!*****************************************
!*
!*  BASISSETINFO INITIATION ROUTINES
!*
!*****************************************

!> \brief Returns number of atoms and regular and auxiliary orbitals of given molecule
!> \author S. Reine
!> \date 2010-02-21
!> \param molecule Contains the information about the molecule
!> \param nAtoms The number of atoms
!> \param nBastReg The number of regular orbitals
!> \param nBastAux The number of auxiliary orbitals
SUBROUTINE getMolecularDimensions(MOLECULE,nAtoms,nBast,nBastAux)
implicit none
TYPE(MOLECULEINFO),intent(in) :: MOLECULE
Integer,intent(out)           :: nAtoms,nBast,nBastAux
!AORdefault   !AODFdefault
nAtoms   = MOLECULE%nAtoms
IF(AORdefault.EQ.AOregular)THEN
   nbast = MOLECULE%nbastREG
ELSEIF(AORdefault.EQ.AOVAL)THEN
   nbast = MOLECULE%nbastVAL
ELSEIF(AORdefault.EQ.AOdfAux)THEN
   nbast = MOLECULE%nBastAUX
ELSEIF(AORdefault.EQ.AOdfCABS)THEN
   nbast = MOLECULE%nBastCABS + MOLECULE%nbastREG
ELSEIF(AORdefault.EQ.AOdfCABO)THEN !CABS only
   nbast = MOLECULE%nBastCABS
ELSEIF(AORdefault.EQ.AOdfJK)THEN
   nbast = MOLECULE%nBastJK
ELSEIF(AORdefault.EQ.AOadmm)THEN
   nbast = MOLECULE%nBastADMM
ELSE   
   CALL LSQUIT('ERROR IN NBASIS DETERMINATION in getMolecularDimensions',-1)
ENDIF

IF(AODFdefault.EQ.AOregular)THEN
   nbastAux = MOLECULE%nbastREG
ELSEIF(AODFdefault.EQ.AOVAL)THEN
   nbastAux = MOLECULE%nbastVAL
ELSEIF(AODFdefault.EQ.AOdfAux)THEN
   nBastAux = MOLECULE%nBastAUX
ELSEIF(AODFdefault.EQ.AOdfCABS)THEN
   nBastAux = MOLECULE%nBastCABS + MOLECULE%nbastREG
ELSEIF(AODFdefault.EQ.AOdfCABO)THEN !CABS only
   nBastAux = MOLECULE%nBastCABS
ELSEIF(AODFdefault.EQ.AOdfJK)THEN
   nBastAux = MOLECULE%nBastJK
ELSEIF(AODFdefault.EQ.AOadmm)THEN
   nBastAux = MOLECULE%nBastADMM
ELSE
   CALL LSQUIT('ERROR IN NAUX DETERMINATION in getMolecularDimensions',-1)
ENDIF
!
END SUBROUTINE getMolecularDimensions

!> \brief Sets up the orbital information for a given molecule
!> \author S. Reine
!> \date 2010-02-21
!> \param Molecule Contains the information about the molecule
!> \param orbitalInfo Contains orbital indeces for the different atoms
SUBROUTINE setMolecularOrbitalInfo(MOLECULE,orbitalInfo)
implicit none
TYPE(MOLECULEINFO),intent(in)            :: MOLECULE
TYPE(MOLECULARORBITALINFO),intent(inout) :: orbitalInfo
!
integer :: iAtom,iReg,iAux,nOrbReg,nOrbAux


CALL getMolecularDimensions(MOLECULE,orbitalInfo%nAtoms,orbitalInfo%nBastReg,&
   &                        orbitalInfo%nBastAux)

CALL initMolecularOrbitalInfo(orbitalInfo,orbitalInfo%nAtoms)

iReg = 0
iAux = 0
DO iAtom=1,orbitalInfo%nAtoms
   IF(AORdefault.EQ.AOregular)THEN
      nOrbReg = MOLECULE%ATOM(iAtom)%nContOrbREG
   ELSEIF(AORdefault.EQ.AOVAL)THEN
      nOrbReg = MOLECULE%ATOM(iAtom)%nContOrbVAL
   ELSEIF(AORdefault.EQ.AOdfAux)THEN
      nOrbReg = MOLECULE%ATOM(iAtom)%nContOrbAUX
   ELSEIF(AORdefault.EQ.AOdfCABS)THEN
      nOrbReg = MOLECULE%ATOM(iAtom)%nContOrbCABS+MOLECULE%ATOM(iAtom)%nContOrbREG
   ELSEIF(AORdefault.EQ.AOdfCABO)THEN !CABS only
      nOrbReg = MOLECULE%ATOM(iAtom)%nContOrbCABS
   ELSEIF(AORdefault.EQ.AOADMM)THEN
      nOrbReg = MOLECULE%ATOM(iAtom)%nContOrbADMM
   ELSEIF(AORdefault.EQ.AOdfJK)THEN
      nOrbReg = MOLECULE%ATOM(iAtom)%nContOrbJK
   ENDIF
   IF(AODFdefault.EQ.AOregular)THEN
      nOrbAux = MOLECULE%ATOM(iAtom)%nContOrbREG
   ELSEIF(AODFdefault.EQ.AOVAL)THEN
      nOrbAux = MOLECULE%ATOM(iAtom)%nContOrbVAL
   ELSEIF(AODFdefault.EQ.AOdfAux)THEN
      nOrbAux = MOLECULE%ATOM(iAtom)%nContOrbAUX
   ELSEIF(AODFdefault.EQ.AOdfCABS)THEN
      nOrbAux = MOLECULE%ATOM(iAtom)%nContOrbCABS+MOLECULE%ATOM(iAtom)%nContOrbREG
   ELSEIF(AODFdefault.EQ.AOdfCABO)THEN !CABS only
      nOrbAux = MOLECULE%ATOM(iAtom)%nContOrbCABS
   ELSEIF(AODFdefault.EQ.AOADMM)THEN
      nOrbAux = MOLECULE%ATOM(iAtom)%nContOrbADMM
   ELSEIF(AODFdefault.EQ.AOdfJK)THEN
      nOrbAux = MOLECULE%ATOM(iAtom)%nContOrbJK
   ENDIF
   orbitalInfo%numAtomicOrbitalsReg(iAtom)=nOrbReg
   orbitalInfo%numAtomicOrbitalsAux(iAtom)=nOrbAux
   orbitalInfo%startAtomicOrbitalsReg(iAtom) = iReg+1
   orbitalInfo%startAtomicOrbitalsAux(iAtom) = iAux+1
   iReg = iReg + nOrbReg
   iAux = iAux + nOrbAux
   orbitalInfo%endAtomicOrbitalsReg(iAtom) = iReg
   orbitalInfo%endAtomicOrbitalsAux(iAtom) = iAux
ENDDO
END SUBROUTINE setMolecularOrbitalInfo

!> \brief Returns the orbital information of a given atom
!> \author S. Reine
!> \date 2010-02-19
!> \param orbitalInfo Contains the orbital-information of a given molecule
!> \param iAtom Specifies the atomic number in question
!> \param nReg The number of regular basis functions on given atom
!> \param startReg The starting orbital index of given atom
!> \param endReg The last orbital index of given atom
!> \param nAux The number of auxiliary basis functions on given atom
!> \param startAux The starting auxiliary orbital index of given atom
!> \param endAux The last auxiliary orbital index of given atom
SUBROUTINE getAtomicOrbitalInfo(orbitalInfo,iAtom,nReg,startReg,endReg,nAux,startAux,endAux)
use typedef
implicit none
TYPE(MOLECULARORBITALINFO),intent(IN) :: orbitalInfo
Integer,intent(IN)  :: iAtom
Integer,intent(OUT) :: nReg,startReg,endReg,nAux,startAux,endAux
!
nReg     = orbitalInfo%numAtomicOrbitalsReg(iAtom)
startReg = orbitalInfo%startAtomicOrbitalsReg(iAtom)
endReg   = orbitalInfo%endAtomicOrbitalsReg(iAtom)
nAux     = orbitalInfo%numAtomicOrbitalsAux(iAtom)
startAux = orbitalInfo%startAtomicOrbitalsAux(iAtom)
endAux   = orbitalInfo%endAtomicOrbitalsAux(iAtom)
!
END SUBROUTINE getAtomicOrbitalInfo

!> \brief Initialize orbitalInfo type
!> \author S. Reine
!> \date 2010-02-21
!> \param orbitalInfo Contains orbital indeces for the different atoms
!> \param nAtoms The number of atoms
SUBROUTINE initMolecularOrbitalInfo(orbitalInfo,nAtoms)
implicit none
TYPE(MOLECULARORBITALINFO),intent(inout) :: orbitalInfo
integer,intent(in)                       :: nAtoms

CALL mem_alloc(orbitalInfo%numAtomicOrbitalsReg,nAtoms)
CALL mem_alloc(orbitalInfo%startAtomicOrbitalsReg,nAtoms)
CALL mem_alloc(orbitalInfo%endAtomicOrbitalsReg,nAtoms)
CALL mem_alloc(orbitalInfo%numAtomicOrbitalsAux,nAtoms)
CALL mem_alloc(orbitalInfo%startAtomicOrbitalsAux,nAtoms)
CALL mem_alloc(orbitalInfo%endAtomicOrbitalsAux,nAtoms)

END SUBROUTINE initMolecularOrbitalInfo

 
!> \brief Frees orbitalInfo type
!> \author S. Reine
!> \date 2010-02-21
!> \param orbitalInfo Contains orbital indeces for the different atoms
SUBROUTINE freeMolecularOrbitalInfo(orbitalInfo)
implicit none
TYPE(MOLECULARORBITALINFO),INTENT(INOUT) :: orbitalInfo
!
CALL mem_dealloc(orbitalInfo%numAtomicOrbitalsReg)
CALL mem_dealloc(orbitalInfo%startAtomicOrbitalsReg)
CALL mem_dealloc(orbitalInfo%endAtomicOrbitalsReg)
CALL mem_dealloc(orbitalInfo%numAtomicOrbitalsAux)
CALL mem_dealloc(orbitalInfo%startAtomicOrbitalsAux)
CALL mem_dealloc(orbitalInfo%endAtomicOrbitalsAux)
END SUBROUTINE freeMolecularOrbitalInfo

!> \brief Determined the number of electrons for a given molecule
!> \author T. Kjaergaard
!> \date 2010-02-21
!> \param Molecule Contains the information about the molecule
!> \param Moleculecharge The charge of the molecule
!> \param nelectrons The number of electrons
SUBROUTINE DETERMINE_NELECTRONS(Molecule,Moleculecharge,nelectrons)
implicit none
TYPE(MOLECULEINFO),intent(IN) :: MOLECULE
real(realk),intent(IN)        :: Moleculecharge
integer,intent(OUT)           :: Nelectrons
!
integer :: I,NCHARGE

NCHARGE = 0
DO I = 1,MOLECULE%NATOMS
   IF(MOLECULE%ATOM(I)%phantom)CYCLE !no electrons on this  
   NCHARGE = INT(MOLECULE%ATOM(I)%CHARGE)+NCHARGE
ENDDO

NELECTRONS = NCHARGE - INT(Moleculecharge)

END SUBROUTINE DETERMINE_NELECTRONS

!!$!> \brief Divide a molecule into molecular fragments (by setting up indices)
!!$!> \author S. Reine
!!$!> \date 2010-02-05
!!$!> \param Molecule The molecule to be fragmented
!!$!> \param fragmentIndex Indices specifying for each atom in Molecule which fragment it belongs to
!!$!> \param numFragments The number of fragments the molecule should be divied into
!!$!> \param lupri Default output unit
!!$SUBROUTINE fragmentMolecule(Molecule,fragmentIndex,numFragments,lupri)
!!$implicit none
!!$TYPE(MOLECULEINFO),intent(in) :: MOLECULE
!!$Integer,intent(in)            :: numFragments,lupri
!!$Integer,intent(inout)         :: fragmentIndex(MOLECULE%nAtoms)
!!$!
!!$Integer :: numOrbitals,numFragOrbitals,iFragment,I,totOrb
!!$logical :: Increased
!!$
!!$
!!$IF (numFragments.GT.MOLECULE%nAtoms) THEN
!!$  CALL LSQUIT('ERROR: fragmentMolecule entered with numFragments > nAtoms',lupri)
!!$ELSEIF (numFragments.EQ.MOLECULE%nAtoms) THEN
!!$  TOTorb=0
!!$  iFragment = 0
!!$  DO I=1,MOLECULE%nAtoms
!!$    numOrbitals = MOLECULE%ATOM(I)%nContOrbREG
!!$    TOTorb = TOTorb + MOLECULE%ATOM(I)%nContOrbREG
!!$    iFragment = iFragment + 1
!!$    fragmentIndex(I) = iFragment
!!$  ENDDO
!!$ELSE
!!$! Divide the molecule into fragments with approximately similar number of
!!$! orbitals
!!$
!!$! First get the average number of orbitals per fragment
!!$  numFragOrbitals = MOLECULE%nbastREG/numFragments
!!$
!!$! Then partition the molecule into fragments of about this number of orbitails
!!$  TOTorb=0
!!$  numOrbitals = 0
!!$  iFragment   = 1
!!$  Increased = .FALSE.
!!$  DO I=1,MOLECULE%nAtoms
!!$    numOrbitals = numOrbitals + MOLECULE%ATOM(I)%nContOrbREG
!!$    TOTorb = TOTorb + MOLECULE%ATOM(I)%nContOrbREG
!!$    fragmentIndex(I) = iFragment
!!$    Increased = .TRUE.
!!$    IF((TOTorb .GE. iFragment*numFragOrbitals .AND. .NOT. (ifragment.EQ.numFragments) ))THEN
!!$      iFragment = iFragment + 1
!!$      numOrbitals = 0
!!$      Increased = .FALSE.
!!$    ENDIF
!!$  ENDDO
!!$  IF(.NOT.Increased) ifragment=ifragment-1
!!$  IF(iFragment .NE. numFragments) THEN
!!$     TOTorb=0
!!$     numOrbitals = 0
!!$     iFragment   = numFragments
!!$     Increased = .FALSE.
!!$     DO I=MOLECULE%nAtoms,1,-1
!!$        numOrbitals = numOrbitals + MOLECULE%ATOM(I)%nContOrbREG
!!$        TOTorb = TOTorb + MOLECULE%ATOM(I)%nContOrbREG
!!$        fragmentIndex(I) = iFragment
!!$        Increased = .TRUE.
!!$        IF (TOTorb .GE. (numFragments-iFragment+1)*numFragOrbitals .AND. .NOT. ((numFragments-ifragment+1).EQ.numFragments)) THEN
!!$           iFragment = iFragment - 1
!!$           numOrbitals = 0
!!$           Increased = .FALSE.
!!$        ENDIF
!!$     ENDDO
!!$     IF(.NOT.Increased)ifragment=ifragment+1
!!$     ifragment = numFragments-iFragment+1
!!$  ENDIF
!!$  IF(iFragment .NE. numFragments) THEN
!!$     WRITE(LUPRI,*)'FRAGEMENT ',iFragment
!!$     WRITE(LUPRI,*)'NODES     ',numFragments
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     WRITE(LUPRI,*)'WARNING WARNING WANING '
!!$     CALL LSQUIT('ifrag not equal to number of nodes',lupri)
!!$  ENDIF
!!$ENDIF
!!$
!!$END SUBROUTINE fragmentMolecule

!> \brief Builds a molecular fragment from a subset of the atoms in the original molecule
!> \author S. Reine and T. Kjaergaard
!> \date 2010-02-21
!> \param DALMOL The original molecule
!> \param FRAGMOL The fragment molecule
!> \param FRAGBASIS the basisinfo 
!> \param ATOMS List of atoms to be included in the fragment
!> \param nATOMS The number of atoms to be included
!> \param lupri Default output unit
SUBROUTINE BUILD_FRAGMENT(DALMOL,FRAGMOL,FRAGBASIS,ATOMS,nATOMS,lupri)
implicit none
INTEGER,intent(IN)                    :: NATOMS,lupri
INTEGER,intent(IN)                    :: ATOMS(NATOMS)
TYPE(MOLECULEINFO),intent(IN)         :: DALMOL
TYPE(MOLECULEINFO),intent(INOUT)      :: FRAGMOL
TYPE(BASISINFO),intent(INOUT)         :: FRAGBASIS
call BUILD_FRAGMENT2(DALMOL,FRAGMOL,FRAGBASIS,ATOMS,nATOMS,lupri)
call DETERMINE_FRAGMENTNBAST(DALMOL,FRAGMOL,FRAGBASIS,lupri)
END SUBROUTINE BUILD_FRAGMENT

!> \brief Sets info about nbast
!> \author S. Reine and T. Kjaergaard
!> \date 2010-02-21
!> \param DALMOL The original molecule
!> \param FRAGMOL The fragment molecule
!> \param FRAGBASIS the basisinfo 
!> \param lupri Default output unit
subroutine DETERMINE_FRAGMENTNBAST(DALMOL,FRAGMOL,FRAGBASIS,lupri)
implicit none
INTEGER,intent(IN)                    :: lupri
TYPE(MOLECULEINFO),intent(IN)         :: DALMOL
TYPE(MOLECULEINFO),intent(INOUT)      :: FRAGMOL
TYPE(BASISINFO),intent(INOUT)         :: FRAGBASIS
!
INTEGER            :: I

do I=1,nBasisBasParam
   IF(FRAGBASIS%WBASIS(I))THEN
      CALL DETERMINE_NBAST(FRAGMOL,FRAGBASIS%BINFO(I))
   ENDIF
enddo
CALL DETERMINE_NBAST(FRAGMOL,FRAGBASIS%BINFO(RegBasParam))
end subroutine DETERMINE_FRAGMENTNBAST

!> \brief Builds a molecular fragment from a subset of the atoms in the original molecule
!> \author S. Reine and T. Kjaergaard
!> \date 2010-02-21
!> \param DALMOL The original molecule
!> \param FRAGMOL The fragment molecule
!> \param FRAGBASIS the basisinfo 
!> \param ATOMS List of atoms to be included in the fragment
!> \param nATOMS The number of atoms to be included
!> \param lupri Default output unit
SUBROUTINE BUILD_FRAGMENT2(DALMOL,FRAGMOL,FRAGBASIS,ATOMS,nATOMS,lupri)
implicit none
INTEGER,intent(IN)                    :: NATOMS,lupri
INTEGER,intent(IN)                    :: ATOMS(NATOMS)
TYPE(MOLECULEINFO),intent(IN)         :: DALMOL
TYPE(MOLECULEINFO),intent(INOUT)      :: FRAGMOL
TYPE(BASISINFO),intent(INOUT)         :: FRAGBASIS
!
INTEGER            :: I,nelectrons
Character(len=22)  :: FRAGMENTNAME

IF ((NATOMS.GT. 999999).OR.(ATOMS(1).GT. 999999).OR.(ATOMS(NATOMS).GT. 999999)) &
     & CALL LSQUIT('Error in BUILD_FRAGMENT -> FRAGMENTNAME',-1)
write(FRAGMENTNAME,'(A4,3I6)') 'FRAG',NATOMS,ATOMS(1),ATOMS(NATOMS)
CALL init_MoleculeInfo(FRAGMOL,natoms,FRAGMENTNAME)
FRAGMOL%charge = DALMOL%charge
nelectrons = 0
DO I = 1,nAtoms
   CALL COPY_ATOM(DALMOL,ATOMS(I),FRAGMOL,I,lupri)
   nelectrons = nelectrons + INT(FRAGMOL%ATOM(I)%Charge)
ENDDO
FRAGMOL%nelectrons = nelectrons
FRAGMOL%nbastREG=0
FRAGMOL%nprimbastREG=0
FRAGMOL%nbastAUX=0
FRAGMOL%nprimbastAUX=0
FRAGMOL%nbastCABS=0
FRAGMOL%nprimbastCABS=0
FRAGMOL%nbastADMM=0
FRAGMOL%nprimbastADMM=0
FRAGMOL%nbastJK=0
FRAGMOL%nprimbastJK=0
FRAGMOL%nbastVAL=0
FRAGMOL%nprimbastVAL=0

FRAGMOL%nSubSystems = DALMOL%nSubSystems
IF(FRAGMOL%nSubSystems.NE.0)THEN
   call mem_alloc(FRAGMOL%SubSystemLabel,FRAGMOL%nSubSystems)
   do I = 1, FRAGMOL%nSubSystems
      FRAGMOL%SubSystemLabel(I) = DALMOL%SubSystemLabel(I)
   enddo
ELSE
   NULLIFY(FRAGMOL%SubSystemLabel)
ENDIF

END SUBROUTINE BUILD_FRAGMENT2

!!$!> \brief 
!!$!> \author
!!$!> \date
!!$!> \param 
!!$SUBROUTINE buildFragmentFromFragmentIndex(FRAGMENT,MOLECULE,FragmentIndex,iFrag,lupri)
!!$implicit none
!!$TYPE(MOLECULEINFO) :: FRAGMENT
!!$TYPE(MOLECULEINFO),intent(IN)  :: MOLECULE
!!$Integer,intent(IN)             :: FragmentIndex(MOLECULE%nAtoms)
!!$Integer,intent(IN)             :: iFrag
!!$Integer,intent(IN)             :: lupri
!!$!
!!$Integer :: iAtom
!!$Integer :: nAtoms
!!$!
!!$nAtoms=0
!!$Do iAtom=1,MOLECULE%nATOMS
!!$  IF(FragmentIndex(iAtom) .EQ. iFrag)THEN
!!$    nAtoms=nAtoms+1
!!$    CALL COPY_ATOM(MOLECULE,iAtom,FRAGMENT,nAtoms,lupri)
!!$  ENDIF
!!$ENDDO
!!$FRAGMENT%nAtoms = nAtoms
!!$!
!!$END SUBROUTINE buildFragmentFromFragmentIndex


!> \brief Free the dalton-fragments
!> \autor S. Reine
!> \param Setting Contains the integral settings
SUBROUTINE freeDaltonFragments(SETTING)
implicit none
Type(LSSETTING),intent(inout) :: SETTING
!
Integer :: indAO,I

CALL freeFragments(SETTING%FRAGMENT,SETTING%fragBuild,setting%nAO)

! Restore to defaul settings
DO indAO=1,setting%nAO
 setting%fragment(indAO)%p => setting%molecule(indAO)%p
 IF (associated(setting%fragment(indAO)%p)) THEN
   DO I=1,nBasisBasParam
      call DETERMINE_NBAST(setting%MOLECULE(indAO)%p,setting%BASIS(indAO)%p%BINFO(I),&
           & setting%scheme%DoSpherical,setting%scheme%uncont)
   ENDDO
 ENDIF
ENDDO
!Set up fragments
END SUBROUTINE FreeDaltonFragments

!> \brief free the molecule fragments 
!> \author S. Reine
!> \date 2010
SUBROUTINE freeFragments(FRAGMENT,fragBuild,nAO)
implicit none
INTEGER, intent(in)             :: nAO
TYPE(MOLECULE_PT),intent(inout) :: FRAGMENT(nAO)
LOGICAL,intent(inout)           :: fragBuild(nAO)
!
integer :: indAO

DO indAO = 1,nAO
  IF (fragBuild(indAO)) THEN
    CALL free_MoleculeInfo(FRAGMENT(indAO)%p)
    DEALLOCATE(FRAGMENT(indAO)%p)
    NULLIFY(FRAGMENT(indAO)%p)
   fragBuild(indAO) = .FALSE.
  ENDIF
ENDDO
END SUBROUTINE freeFragments


!> \brief 
!> \author
!> \date
!> \param 
SUBROUTINE DETERMINE_NBAST(MOLECULE,BASINFO,spherical,UNCONTRACTED)
implicit none
TYPE(BASISSETINFO)  :: BASINFO
TYPE(MOLECULEINFO)  :: MOLECULE
INTEGER             :: I,TOTcont,TOTprim,R,K,type,lupri,TEMP1,TEMP2,icharge
LOGICAL,OPTIONAL    :: spherical,UNCONTRACTED
!
Logical :: spher, uncont,REG,AUX,VAL,JKAUX,CABS,ADMM

IF(BASINFO%natomtypes.NE.0)THEN
 IF(BASINFO%label(1:9) .NE. 'GCTRANS  ')THEN 
   ! Defaults
   spher  = .true.
   uncont = .false.
   ! Optional settings
   IF (present(spherical)) spher = spherical
   IF (present(UNCONTRACTED)) uncont = UNCONTRACTED
   REG = .FALSE.
   AUX = .FALSE.
   CABS = .FALSE.
   JKAUX = .FALSE.
   VAL = .FALSE.
   ADMM = .FALSE.
   IF(BASINFO%label(1:9) .EQ. 'REGULAR  ') REG = .TRUE.
   IF(BASINFO%label(1:9) .EQ. 'AUXILIARY') AUX = .TRUE.
   IF(BASINFO%label(1:9) .EQ. 'CABS     ') CABS = .TRUE.
   IF(BASINFO%label(1:9) .EQ. 'JKAUX    ') JKAUX = .TRUE.
   IF(BASINFO%label(1:9) .EQ. 'ADMM     ') ADMM = .TRUE.
   IF(BASINFO%label(1:9) .EQ. 'VALENCE  ') VAL = .TRUE.
   IF(.NOT.MOLECULE%pointMolecule)THEN
      TOTcont=0
      TOTprim=0
      R = BASINFO%Labelindex
      DO I=1,MOLECULE%nAtoms
         IF(R.EQ. 0)THEN
            icharge = INT(MOLECULE%ATOM(I)%charge)
            type = BASINFO%chargeindex(icharge) 
         ELSE
            type=MOLECULE%ATOM(I)%IDtype(R)
         ENDIF
         IF(.NOT.MOLECULE%ATOM(I)%Pointcharge)THEN
            IF(uncont)THEN
               TOTcont=TOTcont+BASINFO%ATOMTYPE(type)%Totnprim
               IF(REG) MOLECULE%ATOM(I)%nprimOrbREG=BASINFO%ATOMTYPE(type)%Totnprim
               IF(AUX) MOLECULE%ATOM(I)%nprimOrbAUX=BASINFO%ATOMTYPE(type)%Totnprim
               IF(CABS) MOLECULE%ATOM(I)%nprimOrbCABS=BASINFO%ATOMTYPE(type)%Totnprim
               IF(JKAUX) MOLECULE%ATOM(I)%nprimOrbJK=BASINFO%ATOMTYPE(type)%Totnprim
               IF(ADMM) MOLECULE%ATOM(I)%nprimOrbADMM=BASINFO%ATOMTYPE(type)%Totnprim
               IF(VAL) MOLECULE%ATOM(I)%nprimOrbVAL=BASINFO%ATOMTYPE(type)%Totnprim
               IF(REG) MOLECULE%ATOM(I)%ncontOrbREG=BASINFO%ATOMTYPE(type)%Totnprim
               IF(AUX) MOLECULE%ATOM(I)%ncontOrbAUX=BASINFO%ATOMTYPE(type)%Totnprim
               IF(CABS) MOLECULE%ATOM(I)%ncontOrbCABS=BASINFO%ATOMTYPE(type)%Totnprim
               IF(JKAUX) MOLECULE%ATOM(I)%ncontOrbJK=BASINFO%ATOMTYPE(type)%Totnprim
               IF(ADMM) MOLECULE%ATOM(I)%ncontOrbADMM=BASINFO%ATOMTYPE(type)%Totnprim
               IF(VAL) MOLECULE%ATOM(I)%ncontOrbVAL=BASINFO%ATOMTYPE(type)%Totnprim
               TOTprim=TOTcont
            ELSE !DEFAULT
               TOTcont=TOTcont+BASINFO%ATOMTYPE(type)%Totnorb      
               TOTprim=TOTprim+BASINFO%ATOMTYPE(type)%Totnprim      
               IF(REG) MOLECULE%ATOM(I)%nprimOrbREG=BASINFO%ATOMTYPE(type)%Totnprim      
               IF(AUX) MOLECULE%ATOM(I)%nprimOrbAUX=BASINFO%ATOMTYPE(type)%Totnprim      
               IF(CABS) MOLECULE%ATOM(I)%nprimOrbCABS=BASINFO%ATOMTYPE(type)%Totnprim
               IF(JKAUX) MOLECULE%ATOM(I)%nprimOrbJK=BASINFO%ATOMTYPE(type)%Totnprim
               IF(ADMM) MOLECULE%ATOM(I)%nprimOrbADMM=BASINFO%ATOMTYPE(type)%Totnprim
               IF(VAL) MOLECULE%ATOM(I)%nprimOrbVAL=BASINFO%ATOMTYPE(type)%Totnprim      
               IF(REG) MOLECULE%ATOM(I)%ncontOrbREG=BASINFO%ATOMTYPE(type)%Totnorb      
               IF(AUX) MOLECULE%ATOM(I)%ncontOrbAUX=BASINFO%ATOMTYPE(type)%Totnorb      
               IF(CABS) MOLECULE%ATOM(I)%ncontOrbCABS=BASINFO%ATOMTYPE(type)%Totnorb      
               IF(JKAUX) MOLECULE%ATOM(I)%ncontOrbJK=BASINFO%ATOMTYPE(type)%Totnorb      
               IF(ADMM) MOLECULE%ATOM(I)%ncontOrbADMM=BASINFO%ATOMTYPE(type)%Totnorb      
               IF(VAL) MOLECULE%ATOM(I)%ncontOrbVAL=BASINFO%ATOMTYPE(type)%Totnorb      
            ENDIF
         ELSE
            IF(REG) MOLECULE%ATOM(I)%nprimOrbREG=0
            IF(AUX) MOLECULE%ATOM(I)%nprimOrbAUX=0
            IF(CABS) MOLECULE%ATOM(I)%nprimOrbCABS=0
            IF(JKAUX) MOLECULE%ATOM(I)%nprimOrbJK=0
            IF(ADMM) MOLECULE%ATOM(I)%nprimOrbADMM=0
            IF(VAL) MOLECULE%ATOM(I)%nprimOrbVAL=0
            IF(REG) MOLECULE%ATOM(I)%ncontOrbREG=0
            IF(AUX) MOLECULE%ATOM(I)%ncontOrbAUX=0
            IF(CABS) MOLECULE%ATOM(I)%ncontOrbCABS=0
            IF(JKAUX) MOLECULE%ATOM(I)%ncontOrbJK=0
            IF(ADMM) MOLECULE%ATOM(I)%ncontOrbADMM=0
            IF(VAL) MOLECULE%ATOM(I)%ncontOrbVAL=0
         ENDIF
      ENDDO
      BASINFO%nbast=TOTcont
      BASINFO%nprimbast=TOTprim
      IF(REG)MOLECULE%nbastREG=TOTcont
      IF(REG)MOLECULE%nprimbastREG=TOTprim
      
      IF(AUX)MOLECULE%nbastAUX=TOTcont
      IF(AUX)MOLECULE%nprimbastAUX=TOTprim
      
      IF(CABS)MOLECULE%nbastCABS=TOTcont
      IF(CABS)MOLECULE%nprimbastCABS=TOTprim
      
      IF(JKAUX)MOLECULE%nbastJK=TOTcont
      IF(JKAUX)MOLECULE%nprimbastJK=TOTprim
      
      IF(ADMM)MOLECULE%nbastADMM=TOTcont
      IF(ADMM)MOLECULE%nprimbastADMM=TOTprim

      IF(VAL)MOLECULE%nbastVAL=TOTcont
      IF(VAL)MOLECULE%nprimbastVAL=TOTprim
   ENDIF
 ENDIF
ENDIF

END SUBROUTINE DETERMINE_NBAST

SUBROUTINE DETERMINE_NBAST2(MOLECULE,BASINFO,spherical,UNCONTRACTED,nbast)
implicit none
INTEGER,intent(inout) :: nbast
TYPE(BASISSETINFO),intent(in)  :: BASINFO
TYPE(MOLECULEINFO),intent(in)  :: MOLECULE
LOGICAL,OPTIONAL    :: spherical,UNCONTRACTED
!
INTEGER             :: I,TOTcont,R,K,type,icharge
Logical             :: spher, uncont
!
! Defaults
spher  = .true.
uncont = .false.
! Optional settings
IF (present(spherical)) spher = spherical
IF (present(UNCONTRACTED)) uncont = UNCONTRACTED

TOTcont=0
R = BASINFO%Labelindex
DO I=1,MOLECULE%nAtoms
   IF(R.EQ. 0)THEN
      icharge = INT(MOLECULE%ATOM(I)%charge)
      type = BASINFO%chargeindex(icharge) 
   ELSE
      type=MOLECULE%ATOM(I)%IDtype(R)
   ENDIF
   IF(.NOT.MOLECULE%ATOM(I)%Pointcharge)THEN
      IF(uncont)THEN
         TOTcont=TOTcont+BASINFO%ATOMTYPE(type)%Totnprim
      ELSE !DEFAULT
         TOTcont=TOTcont+BASINFO%ATOMTYPE(type)%Totnorb      
      ENDIF
   ENDIF
ENDDO
nbast=TOTcont

END SUBROUTINE DETERMINE_NBAST2

SUBROUTINE GET_GEOMETRY(LUPRI,IPRINT,MOLECULE,natoms,X,Y,Z)
IMPLICIT NONE
INTEGER            :: LUPRI,IPRINT,natoms
TYPE(MOLECULEINFO) :: MOLECULE
REAL(REALK)        :: X(natoms),Y(natoms),Z(natoms)
!
integer :: I

DO I=1,nAtoms
     X(I) = MOLECULE%ATOM(I)%CENTER(1)
     Y(I) = MOLECULE%ATOM(I)%CENTER(2)
     Z(I) = MOLECULE%ATOM(I)%CENTER(3)
ENDDO

END SUBROUTINE GET_GEOMETRY

SUBROUTINE PRINT_GEOMETRY(MOLECULE,LUPRI)
IMPLICIT NONE
INTEGER :: LUPRI
INTEGER :: I
CHARACTER(len=1)   :: CHRXYZ(3)=(/'x','y','z'/)
TYPE(MOLECULEINFO) :: MOLECULE
   WRITE (LUPRI,'(2X,A,I3)')' Total number of coordinates:',3*MOLECULE%natoms
   WRITE (LUPRI,'(2X,A)')' Written in atomic units    '
   DO I=1,MOLECULE%nAtoms
      WRITE (LUPRI,'(A)')' '
      IF(MOLECULE%ATOM(I)%phantom)THEN
         WRITE (LUPRI,'(2X,A)')' Phantom Atom (Only Basis Functions)'
      ENDIF
      WRITE (LUPRI,'(I4,3X,A,5X,A,3X,F15.10)')&
           &  (3*I-2), MOLECULE%ATOM(I)%Name,CHRXYZ(1),&
           & MOLECULE%ATOM(I)%CENTER(1)
      WRITE (LUPRI,'(I4,12X,A,3X,F15.10)')&
           &  3*I-1, CHRXYZ(2), MOLECULE%ATOM(I)%CENTER(2),&
           &  3*I, CHRXYZ(3), MOLECULE%ATOM(I)%CENTER(3)
   ENDDO
END SUBROUTINE PRINT_GEOMETRY 

END MODULE molecule_module
