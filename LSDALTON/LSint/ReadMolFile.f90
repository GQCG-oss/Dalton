!> @file
!> Read molecule input
MODULE READMOLEFILE
  use fundamental
  use precision
  use TYPEDEF
  use molecule_module
  use fundamental
  use lattice_type
  use lattice_vectors
contains
!> \brief read the molecule file and build the molecule structure 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param molecule the molecule structure to be built
!> \param BASISSETLIBRARY the info about basisset read from mol file
!> \param doprint if the information should be printet
!> \param iprint the printlevel integer, determining how much output should be generated
!> \param DoSpherical if the basis should be a Spherical or cartesian basis
!> \param auxbasis if an auxilliary basis is specified in mol file
!> \param cabsbasis if an Complementary auxilliary basis is specified in molfile
!> \param JKbasis if an JK auxilliary basis is specified in mol file
!>
!> THIS IS THE MAIN DRIVER ROUTINE WHICH READS THE MOLECULE INPUT FILE
!> AND BUILDS THE MOLECULE OBJECT 
!>
!> THE RECIPE:
!> 1. OPEN MOLECULE.INP FILE -> LUINFO
!> 2. READ THE FIRST LINE CONTAINING KEYWORD BASIS,ATOMBASIS,...
!> 3. IN CASE OF BASIS READ THE SECOND LINE CONTAINING BASISSET
!> 4. READ THE 2 LINES OF COMMENTS
!> 5. READ THE FOURTH LINE CONTAINING NUMBER OF ATOMS, MOLECULARCHARGE,..    
!> 6. CALL 'READ_GEOMETRY' WHICH READS THE FIFTH LINE A NUMBER OF TIME 
!>    AND READS THE X,Y,Z COORDINATES OF THE INDIVIDUAL ATOMS
!>
SUBROUTINE READ_MOLFILE_AND_BUILD_MOLECULE(LUPRI,MOLECULE,&
     &BASISSETLIBRARY,doprint,iprint,DoSpherical,auxbasis,&
	 &cabsbasis,jkbasis,latt_config)
implicit none
INTEGER,intent(in)               :: LUPRI,iprint
TYPE(MOLECULEINFO),intent(inout) :: MOLECULE
LOGICAL,intent(inout) :: AUXBASIS,DoSpherical
LOGICAL,intent(in) :: doprint
TYPE(BASISSETLIBRARYITEM),intent(inout) :: BASISSETLIBRARY
TYPE(lvec_list_t),INTENT(INOUT) :: latt_config
!
integer            :: LUINFO
logical            :: PRINTATOMCOORD,file_exist,Angstrom,Symmetry,dopbc
logical            :: ATOMBASIS,BASIS,CABSBASIS,JKBASIS
CHARACTER(len=80)  :: BASISSET,AUXBASISSET,CABSBASISSET,JKBASISSET
integer            :: MolecularCharge,Atomtypes,Totalnatoms,I,IPOS

IF (doprint) WRITE(LUPRI,'(2X,A18,2X,I6)')'MOLPRINT IS SET TO',IPRINT
BASIS=.FALSE.
ATOMBASIS=.FALSE.
AUXBASIS=.FALSE.
CABSBASIS=.FALSE.
JKBASIS=.FALSE.
do I=1,80
   AUXBASISSET(I:I)=' '
enddo
do I=1,80
   CABSBASISSET(I:I)=' '
enddo
do I=1,80
   JKBASISSET(I:I)=' '
enddo

INQUIRE(file='MOLECULE.INP',EXIST=file_exist) 
IF(file_exist)THEN
  LUINFO=-1
  CALL LSOPEN(LUINFO,'MOLECULE.INP','OLD','FORMATTED')
ELSE
  CALL LSQUIT('MOLECULE.INP does not exist',lupri)
ENDIF

CALL Obtain_Totalnatoms(LUPRI,LUINFO,Totalnatoms)

CALL init_MoleculeInfo(Molecule,Totalnatoms,'INPUT-Molecule________')

PRINTATOMCOORD = .FALSE.
IF(Totalnatoms .LE. 11)PRINTATOMCOORD = .TRUE.
CALL READ_LINE1(LUPRI,LUINFO,BASIS,ATOMBASIS)

IF(BASIS)THEN
  CALL READ_LINE2(LUPRI,LUINFO,IPRINT,BASISSET,AUXBASISSET,AUXBASIS,&
       &CABSBASISSET,CABSBASIS,JKBASISSET,JKBASIS)
ENDIF

BASISSETLIBRARY%DunningsBasis = .FALSE.
IF(BASIS)THEN
   BASISSETLIBRARY%BASISSETNAME(1)=BASISSET
   IPOS = INDEX(BASISSET,'cc-pV')
   IF (IPOS .NE. 0) THEN 
      BASISSETLIBRARY%DunningsBasis = .TRUE.
   ENDIF
   BASISSETLIBRARY%BASISSETNAME(2)=AUXBASISSET
   BASISSETLIBRARY%BASISSETNAME(3)=CABSBASISSET
   BASISSETLIBRARY%BASISSETNAME(4)=JKBASISSET
   BASISSETLIBRARY%nbasissets=1
   IF(AUXBASIS) THEN
      BASISSETLIBRARY%nbasissets=2
   ENDIF
   IF(CABSBASIS) THEN
      IF(.NOT.AUXBASIS)THEN
         CALL LSQUIT('Aux is required When supplying a CABS basis set, even if not used',-1)
      ENDIF
      BASISSETLIBRARY%nbasissets=3
   ENDIF
   IF(JKBASIS) THEN
      IF(.NOT.CABSBASIS)THEN
         CALL LSQUIT('CABS is required When supplying a JK basis set, even if not used',-1)
      ENDIF
      BASISSETLIBRARY%nbasissets=4
   ENDIF
ELSE !ATOMBASIS
   !We dont know yet 
ENDIF

CALL READ_COMMENTS(LUPRI,LUINFO,.FALSE.)

CALL READ_LINE4(LUPRI,LUINFO,Atomtypes,DoSpherical,MolecularCharge&
                           &,Angstrom,Symmetry,doprint,&
						   &latt_config%setup_pbclatt)

Molecule%charge = MolecularCharge

IF (Angstrom.AND.doprint)THEN
  WRITE (LUPRI,'(2X,A)')'Coordinates are entered in Angstroms&
                         & and converted to atomic units.'
  WRITE (LUPRI,'(2X,A,F11.8,A2)')'Conversion factor : 1 bohr ='&
                               &,bohr_to_angstrom,' A'
ENDIF

!************************************************
!*****       Read geometry input data       *****
!************************************************
DOPBC=.FALSE.

CALL READ_GEOMETRY(LUPRI,LUINFO,IPRINT,BASISSETLIBRARY,Atomtypes,dopbc,&
     &ATOMBASIS,BASIS,AUXBASIS,CABSBASIS,JKBASIS,Angstrom,MOLECULE,&
     &PRINTATOMCOORD,doprint,latt_config)

CALL DETERMINE_NELECTRONS(Molecule,Molecule%charge,Molecule%nelectrons)

!CALL PRINT_MOLECULEINFO(LUPRI,MOLECULE,BASISSETLIBRARY,1)

!CALL PRINT_BASISSETLIBRARY(LUPRI,BASISSETLIBRARY)

CALL LSCLOSE(LUINFO,'KEEP')

END SUBROUTINE READ_MOLFILE_AND_BUILD_MOLECULE

!> \brief determine the total number of atoms
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param totalnatoms the number of total atoms to be determined
SUBROUTINE Obtain_Totalnatoms(LUPRI,LUINFO,Totalnatoms)
implicit none
integer          :: Totalnatoms,ipos,ipos2,natoms,I,J
integer          :: LUINFO,Atomtypes,numbertypes
CHARACTER(len=80):: LINE
LOGICAL          :: FOUND,BASIS
INTEGER          :: LUPRI
INTEGER          :: MolecularCharge,ios
CHARACTER(len=2) :: SYMTXT,ID3
CHARACTER(len=1) :: KASYM(3,3),CRT
real(realk)      :: AtomicCharge

rewind(LUINFO)
BASIS=.FALSE.

DO
  READ (LUINFO,'(a80)') LINE
  IF (LINE(1:1) == '!' .or. LINE(1:1) == '#') CYCLE
  EXIT
ENDDO

IPOS = INDEX(LINE,'ATOM')
IF (IPOS .NE. 0) THEN !ATOMBASIS
  IPOS2 = INDEX(LINE(IPOS:80),'BASIS')
  IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. (IPOS+10) )) THEN
     WRITE (LUPRI,'(2X,a10,a70)') ' Keyword ',LINE
     WRITE (LUPRI,'(2X,a80)') ' not recognized in READ_MOLECULE_FILE.'
     WRITE (LUPRI,'(2X,a80)') ' maybe it is not implemented yet.'
     WRITE (LUPRI,'(2X,a80)') ' format is ATOMBASIS'
     CALL LSQUIT('Only ATOMBASIS and BASIS is implemented',lupri)
  ENDIF
ELSE
  IPOS2 = INDEX(LINE,'BASIS')
  IF (IPOS2 .EQ. 0) THEN
     WRITE (LUPRI,'(2X,a10,a70)') ' Keyword ',LINE
     WRITE (LUPRI,'(2X,a80)') ' not recognized in READ_MOLECULE_FILE.'
     WRITE (LUPRI,'(2X,a80)') ' maybe it is not implemented yet.'
     WRITE (LUPRI,'(2X,a80)') ' format is BASIS'
     CALL LSQUIT('Only ATOMBASIS and BASIS is implemented',lupri)
  ELSE
    BASIS=.TRUE.
  ENDIF
ENDIF

IF(BASIS) THEN
   DO
      READ (LUINFO,'(a80)') LINE  !BASISSETLINE
      IF (LINE(1:1) == '!' .or. LINE(1:1) == '#') CYCLE
      EXIT
   ENDDO
ENDIF

CALL READ_COMMENTS(LUPRI,LUINFO,.FALSE.)
Atomtypes = 0
Totalnatoms=0
FOUND=.FALSE.
DO   
  IF (FOUND) EXIT
  READ (LUINFO,'(a80)') LINE
  IPOS = INDEX(LINE,'Atomtypes')
  IF (IPOS .EQ. 0) IPOS = INDEX(LINE,'AtomTypes')
  IF (IPOS .NE. 0) THEN
    IPOS2 = INDEX(LINE(IPOS:80),'=')
    IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 10)) THEN
      WRITE (LUPRI,'(2X,A40)') 'Incorrect input for # atomtypes'
      WRITE (LUPRI,'(2X,A40)') 'Format is "Atomtypes=?"'
      CALL LSQUIT('Incorrect input for # atomtypes',lupri)
    ELSE
      READ (LINE((IPOS+IPOS2):80),*) Atomtypes
      FOUND=.TRUE.
    ENDIF
  ELSE !ASSUME OLD INPUT
     READ (LINE,'(BN,A1,I4,I3,A2,10A1,D10.2,6I5)',IOSTAT=ios) CRT,&
          & Atomtypes,MolecularCharge,SYMTXT,((KASYM(I,J),I=1,3),J=1,3),&
          & ID3
     IF(IOS .NE. 0) THEN
        WRITE (LUPRI,'(2X,A40)') ' Error in the determination of the number &
             & of atomic types'
        WRITE (LUPRI,'(2X,A60)') "Correct input structure is: Atomtypes=???"
        CALL LSQUIT('Error in determining the number of different &
             & atom types',lupri)
     ELSE
     ENDIF
     EXIT
  ENDIF
ENDDO

numbertypes=1
DO 
  IF(numbertypes>Atomtypes) EXIT
  READ (LUINFO,'(a80)') LINE
  IPOS = INDEX(LINE,'Atoms')
  IF (IPOS .NE. 0) THEN
    IPOS2 = INDEX(LINE(IPOS:),'=')
    IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
      WRITE (LUPRI,'(2X,A40)') 'Obtain_Totalnatoms: Incorrect input for # of atoms'
      WRITE (LUPRI,'(2X,A40)') 'Obtain_Totalnatoms: Format is "Atoms=?"'
      CALL LSQUIT('Obtain_Totalnatoms: Incorrect input for # of atoms',lupri)
    ELSE
      READ (LINE((IPOS+IPOS2):),*) nAtoms
      numbertypes=numbertypes+1
    ENDIF
    Totalnatoms=Totalnatoms+nAtoms
    DO J=1,nAtoms
       READ (LUINFO, '(a80)') LINE
    ENDDO
  ELSE !OLD INPUT
     READ (LINE,'(BN,F10.0,I5)',IOSTAT=ios) AtomicCharge, nAtoms
     IF(IOS .NE. 0) THEN
        WRITE (LUPRI,'(2X,A40)') 'Obtain_Totalnatoms: Error in the determination of the number of atoms'
        WRITE (LUPRI,'(2X,A40)') "Obtain_Totalnatoms: Correct input structure is: Atoms=???"
        CALL LSQUIT('Obtain_Totalnatoms: Error in determining the number of atoms',lupri)
     ENDIF
     numbertypes=numbertypes+1
     Totalnatoms=Totalnatoms+nAtoms
     DO J=1,nAtoms
        READ (LUINFO, '(a80)') LINE
     ENDDO
  ENDIF
ENDDO
!CALL LSCLOSE(LUINFO,'KEEP')

END SUBROUTINE Obtain_Totalnatoms

!> \brief read the first line of molecule file (BASIS or ATOMBASIS)
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param basis if the BASIS is used
!> \param atombasis if the ATOMBASIS is used
SUBROUTINE READ_LINE1(LUPRI,LUINFO,BASIS,ATOMBASIS)
implicit none
CHARACTER(len=80)  :: WORD
LOGICAL            :: BASIS,ATOMBASIS
INTEGER            :: LUINFO,IPOS,IPOS2,LUPRI

rewind(LUINFO)
DO
  READ (LUINFO,'(a80)') WORD
  IF(WORD(1:1) == '!' .OR. WORD(1:1) == '#') CYCLE
  EXIT
ENDDO

BASIS=.FALSE.
ATOMBASIS=.FALSE.

IPOS = INDEX(WORD,'ATOM')
IF (IPOS .NE. 0) THEN !ATOMBASIS
  IPOS2 = INDEX(WORD(IPOS:80),'BASIS')
  IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. (IPOS+10) )) THEN
     WRITE (LUPRI,'(2X,a10,a70)') ' Keyword ',WORD
     WRITE (LUPRI,'(2X,a80)') ' not recognized in READ_MOLECULE_FILE.'
     WRITE (LUPRI,'(2X,a80)') ' maybe it is not implemented yet.'
     WRITE (LUPRI,'(2X,a80)') ' format is ATOMBASIS'
     CALL LSQUIT('Only ATOMBASIS and BASIS is implemented',lupri)
  ELSE
     ATOMBASIS=.TRUE.
  ENDIF
ELSE
  IPOS2 = INDEX(WORD,'BASIS')
  IF (IPOS2 .NE. 0) THEN
     BASIS=.TRUE.
  ELSE
     WRITE (LUPRI,'(2X,a10,a70)') ' Keyword ',WORD
     WRITE (LUPRI,'(2X,a80)') ' not recognized in READ_MOLECULE_FILE.'
     WRITE (LUPRI,'(2X,a80)') ' maybe it is not implemented yet.'
     WRITE (LUPRI,'(2X,a80)') ' format is BASIS'
     CALL LSQUIT('Only ATOMBASIS and BASIS is implemented',lupri)
  ENDIF
ENDIF

END SUBROUTINE READ_LINE1

!> \brief read the second line of molecule file (the name of basisset) if BASIS is used 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param iprint the printlevel integer, determining how much output should be generated
!> \param basisset the name of the basisset 
!> \param auxbasisset the name of the auxilliary basisset 
!> \param auxbasis if an auxilliary basisset is given
SUBROUTINE READ_LINE2(LUPRI,LUINFO,IPRINT,BASISSET,AUXBASISSET,AUXBASIS,&
     & CABSBASISSET,CABSBASIS,JKBASISSET,JKBASIS)
implicit none
CHARACTER(len=80)  :: WORD,BASISSET,AUXBASISSET,CABSBASISSET,JKBASISSET
INTEGER            :: N,I,K,Itmp,M,L
LOGICAL            :: NEXT,AUXBASIS,CABSBASIS,JKBASIS
INTEGER            :: LUINFO
INTEGER            :: LUPRI,IPRINT,BAS
AUXBASIS=.FALSE.
CABSBASIS=.FALSE.
JKBASIS=.FALSE.
N=0
DO I=1,80
  BASISSET(I:I) = ' '
  AUXBASISSET(I:I) = ' '
ENDDO
READ (LUINFO, '(a80)') WORD
BAS=0
NEXT = .FALSE.
DO I=1,80
  IF (WORD(I:I) .NE. ' ' .AND. .NOT. NEXT) THEN
    NEXT=.TRUE.  !INSIDE WORD
    k=i          !Start of word
  ELSEIF(WORD(I:I) .EQ. ' ' .AND. NEXT) THEN
    L=I-1        !end of WORD
    NEXT=.FALSE. !OUT OF WORD
    N=N+1
    IF(N==1)THEN
      BASISSET(1:I-K)=WORD(K:L) !FOUND BASISSET
    ELSE !FOUND ANOTHER BASISSET
       IF(WORD(K:K+2) .EQ. 'Aux')THEN
          AUXBASIS=.TRUE.
          Itmp = K + 3
          BAS=2
       ELSEIF(WORD(K:K+3) .EQ. 'CABS')THEN
          CABSBASIS=.TRUE.
          Itmp = K + 4
          BAS=3
       ELSEIF(WORD(K:K+1) .EQ. 'JK')THEN
          JKBASIS=.TRUE.
          Itmp = K + 2
          BAS=4
       ELSEIF(WORD(K:K) .EQ. '!')THEN !comment
          EXIT
       ELSEIF(WORD(K:K) .EQ. '#')THEN !comment
          EXIT
       ENDIF
       !     Read a portion of the form "Aux = AUXBASNAM", with arbitrarily many blanks
       !     or a portion of the form "CABS = AUXBASNAM", with arbitrarily many blanks
       !     or a portion of the form "JK = AUXBASNAM", with arbitrarily many blanks
       !     1. find '='
       DO WHILE (WORD(Itmp:Itmp) .ne. '=' .and. Itmp .lt. 80)
          Itmp = itmp + 1
       ENDDO
       Itmp = Itmp + 1
       IF (Itmp .ge. 80)THEN
          IF(BAS.EQ.2)CALL LSquit('"Aux" given but no basis provided!',lupri)
          IF(BAS.EQ.3)CALL LSquit('"CABS" given but no basis provided!',lupri)
          IF(BAS.EQ.4)CALL LSquit('"JK" given but no basis provided!',lupri)
       ENDIF
       !     2. find non-blank
       DO WHILE (WORD(Itmp:Itmp) .eq. ' ' .and. itmp .lt. 80)
          Itmp = Itmp + 1
       ENDDO
       If (itmp .ge. 80) CALL LSquit('"Aux" given but no basis provided!',lupri)
       !     3. copy to next blank
       M = 1
       DO WHILE(WORD(itmp:itmp) .ne. ' ' .and. itmp .lt. 80)
          IF(BAS.EQ.2)AUXBASISSET(M:M) = WORD(itmp:itmp)
          IF(BAS.EQ.3)CABSBASISSET(M:M) = WORD(itmp:itmp)
          IF(BAS.EQ.4)JKBASISSET(M:M) = WORD(itmp:itmp)
          M=M+1
          Itmp = Itmp + 1
       ENDDO
    ENDIF
    IF(N>4)THEN !only Basis,Aux,CABS,JK is implemented
        WRITE (LUPRI,'(/A)') ' Found string"'//WORD(K:L)&
        &//'" but not correct format in MOLECULE.INP.'
    ENDIF
  ENDIF
ENDDO
IF (IPRINT .GT. 1) THEN
  WRITE (LUPRI,'(/A)') ' Basis set "'//BASISSET(1:I-K)&
        &//'" from the basis set library will be used.'
  IF(AUXBASIS)THEN
     WRITE (LUPRI,'(/A)') ' AuxBasis set "'//AUXBASISSET(1:M)&
          &//'" from the basis set library will be used.'
  ENDIF
  IF(CABSBASIS)THEN
     WRITE (LUPRI,'(/A)') ' CABSBasis set "'//CABSBASISSET(1:M)&
          &//'" from the basis set library will be used.'
  ENDIF
  IF(JKBASIS)THEN
     WRITE (LUPRI,'(/A)') ' JKBasis set "'//JKBASISSET(1:M)&
          &//'" from the basis set library will be used.'
  ENDIF
END IF
END SUBROUTINE READ_LINE2

!> \brief find N strings (seperated by empty ' ' characters) in the word
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param N the number of strings in the word
!> \param word the character string
SUBROUTINE FIND_N_STRING(LUPRI,N,WORD)
CHARACTER(len=80)  :: RWORD,WORD
INTEGER            :: I,K,L,N,M
LOGICAL            :: NEXT
INTEGER            :: LUPRI

DO I=1,80
  RWORD(I:I) = WORD(I:I)
  WORD(I:I) = ' '
ENDDO
M=0
NEXT=.FALSE.
DO I=1,80
  IF (RWORD(I:I) .NE. ' ' .AND. .NOT. NEXT) THEN
    NEXT=.TRUE.  !INSIDE WORD
    K=I          !Start of word
  ELSEIF(RWORD(I:I) .EQ. ' ' .AND. NEXT) THEN
    L=I-1        !end of WORD
    NEXT=.FALSE. !OUT OF WORD
    M=M+1   
    IF(N==M)THEN
      WORD(1:I-K)=RWORD(K:L)
    ENDIF
    RETURN
  ENDIF
END DO
IF (N<M) THEN
 WRITE(LUPRI,*)'TRYING TO READ MORE LABELS THEN IN STRING'
 CALL LSQUIT('TO FEW LABELS IN FIND WORD.',lupri)
ENDIF

END SUBROUTINE FIND_N_STRING

!> \brief read the 2 comment lines in the molecule input file
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param printing should the comment lines be printed to output file
SUBROUTINE READ_COMMENTS(LUPRI,LUINFO,PRINTING)
implicit none
INTEGER  :: LUINFO,LUPRI
LOGICAL  :: PRINTING
CHARACTER(len=60):: LINE 

READ (LUINFO, '(2X,a80)') LINE
IF (PRINTING) WRITE(LUPRI,*)LINE
READ (LUINFO, '(2X,a80)') LINE
IF (PRINTING) WRITE(LUPRI,*)LINE

END SUBROUTINE READ_COMMENTS

!> \brief read lin 4 in the molecule input file (Atomtypes= Charge= and so on)
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param Atomtypes the number of different atomtypes
!> \param DoSpherical should the basis be spherical or cartesian
!> \param molecularcharge the molecular charge
!> \param angstrom is the coordinates given in angstrom or bohr (au)
!> \param symmetry if we should use symmetry - obsolete 
!> \param doprint if we should print this to output files
SUBROUTINE READ_LINE4(LUPRI,LUINFO,Atomtypes,DoSpherical,MolecularCharge&
                                    &,Angstrom,Symmetry,doprint,setup_pbclatt)
! Read in the fourth molecule.inp line using the new (and old) input scheme
implicit none
INTEGER          :: IPOS,IPOS2,Atomtypes,MolecularCharge,IPSO
Real(realk)      :: Charge
LOGICAL          :: DoSpherical,Angstrom,Symmetry,doprint
CHARACTER(len=2) :: SYMTXT,ID3
CHARACTER(len=80):: LINE 
INTEGER          :: LUINFO,IOS,i,j
CHARACTER(len=1) :: KASYM(3,3),CRT
INTEGER          :: LUPRI
CHARACTER(len=3) :: pbc_check
LOGICAL,INTENT(INOUT) :: setup_pbclatt

READ (LUINFO, '(a80)') LINE
  ! Number of different atom types
  IPOS = INDEX(LINE,'Ato')
  IF (IPOS .NE. 0) THEN
    IPOS2 = INDEX(LINE(IPOS:80),'=')
    IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 10)) THEN
      WRITE (LUPRI,'(2X,A40)') 'Incorrect input for # atomtypes'
      WRITE (LUPRI,'(2X,A40)') 'Format is "Atomtypes=?"'
      CALL LSQUIT('Incorrect input for # atomtypes',lupri)
    ELSE
      READ (LINE((IPOS+IPOS2):80),*) Atomtypes
    ENDIF

    IF (Atomtypes.EQ. 0) THEN
      WRITE (LUPRI,'(A)')&
      &' You have specified a molecule with zero atoms,&
      & thus all answers to all your input are zero! &
      & (or you made an input error in the .mol file)'
      CALL LSQUIT('No atoms according to .mol input!',lupri)
    ELSE IF (Atomtypes.LT. 0) THEN
      WRITE (LUPRI,'(/A,I6)')&
      &' >>> READI1 error, no. of atomic types negative:',Atomtypes
      CALL LSQUIT('Negative number of atoms according to .mol input',lupri)
    ENDIF

!
!     Kind of basis set
!     

    DoSpherical = .TRUE.
    IPOS = INDEX(LINE,'Car')
    IF (IPOS .NE. 0) DoSpherical = .FALSE.
!    DOSPHERICAL = .TRUE.
!    IPOS = INDEX(LINE,'NoSph')
!    IF (IPOS .EQ. 0) THEN
!      DOSPHERICAL = .FALSE.
!    ENDIF


!
!     Charge of molecule
!     
    IPOS = INDEX(LINE,'Cha')
    IF (IPOS .NE. 0) THEN
      IPOS2 = INDEX(LINE(IPOS:80),'=')
      IF ((IPOS2 .EQ. 0) .OR. (IPOS2 .GT. 7)) THEN
        WRITE (LUPRI,'(2X,A40)') 'Incorrect input for molecular charge'
        WRITE (LUPRI,'(2X,A40)') 'Format is "Charge=?"'
        CALL LSQUIT('Incorrect input for molecular charge',lupri)
      ELSE
        READ (LINE(IPOS+IPOS2:80),*) Charge
        MolecularCharge = nint(charge + 1E-8_realk)
        IF (abs(charge - MolecularCharge).GT. 1E-7_realk) THEN
          WRITE(LUPRI,'(2X,A60)') 'Only integer charges allowed in current lsdalton version'
          CALL LSQUIT('Non-integer molecular charge',lupri)
        ENDIF
      ENDIF
    ELSE
        MolecularCharge = 0
    ENDIF

!     
!     Angstrom?
!
    IPOS = INDEX(LINE,'Ang') 
    IF (IPOS .NE. 0) THEN
      Angstrom = .TRUE.
    ELSE
      Angstrom = .FALSE.
    ENDIF

!     
!     Symmetry generators is disregarded at the moment but the old code to read
!     is in abacus/herrdn.F in subroutine LINE4
!   Symmetry turned off always
    SYMMETRY = .FALSE.
!
!     Change of integral threshold?
!
!    IPOS = INDEX(LINE,'Int')
!    IF (IPOS .NE. 0) THEN
!      IPOS2 = INDEX(LINE(IPOS:80),'=')
!      IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 11)) THEN
!        WRITE (LUPRI,'(2X,A40)') 'Incorrect input for integral threshold'
!        WRITE (LUPRI,'(2X,A40)') 'Format is "Integrals=?"'
!        CALL LSQUIT('Incorrect input for integral threshold',lupri)
!      ELSE
!        !SUBROUTINE BY KENNETH THAT SIMULATES A FREEFORMAT READ
!        CALL FREFRM(LINE,IPOS+IPOS2,IDUMMY,IntegralThreshold,'REA',IERR)
!        READ (LINE((IPOS+IPOS2):),*) IntegralThreshold
!        IF(IntegralThreshold>1.00E+0_realk)THEN
!         WRITE (LUPRI,'(2X,A40)') 'Incorrect input for integral threshold'
!         WRITE (LUPRI,'(2X,A40)') 'Format is "Integrals=1.00E-15_realk"'
!         WRITE (LUPRI,'(2X,A40)') 'Your choice was: ',IntegralThreshold
!         CALL LSQUIT('Incorrect input for integral threshold',lupri)
!        ENDIF       
!      ENDIF
!    ELSE
!      IF (doprint) WRITE(LUPRI,'(2X,A49)')'Integral threshold is set to the default 1.00E-15_realk'
!      IntegralThreshold=1.00E-15_realk
!    ENDIF

  !johannesfor reading lattice vectors in pbc
    IPOS = INDEX(LINE,'PBC')
    IF (IPOS .NE. 0) THEN
      setup_pbclatt=.TRUE.
    ENDIF
    IPOS = INDEX(LINE,'pbc')
    IF (IPOS .NE. 0) THEN
      setup_pbclatt=.TRUE.
    ENDIF

  ELSE
!*******************************************************************
!*
!*  OLD INDPUT FORMAT
!*
!******************************************************************

    READ (LINE,'(BN,A1,I4,I3,A2,10A1,D10.2,6I5,A3)',IOSTAT=ios) CRT,&
    & Atomtypes,MolecularCharge,SYMTXT,((KASYM(I,J),I=1,3),J=1,3),&
    & ID3!,pbc_check
    IF(IOS .NE. 0) THEN
      WRITE (LUPRI,'(2X,A)') ' Error in the determination of the number &
     & of atomic types'
      WRITE (LUPRI,*) "Correct input structure is: Atomtypes=???"
      CALL LSQUIT('Error in determining the number of different &
     & atom types',lupri)
    ENDIF
!    IF (LINE(21:30) .EQ. '          ') THEN
!      IF (doprint) WRITE(LUPRI,'(2X,A55)') &
!     &'Integral threshold is set to the default 1.00E-15_realk'
!      IntegralThreshold=1.00E-15_realk
!    ENDIF
    IF (SYMTXT(1:1) .EQ. '0')THEN
       SYMMETRY = .FALSE.
    ELSE IF(SYMTXT(2:2) .EQ. '0') THEN
       SYMMETRY = .FALSE.
    ELSE
       SYMMETRY = .FALSE.       
    ENDIF
    IF(CRT .EQ. 'C' .OR. CRT .EQ. 'c') THEN
      DoSpherical = .FALSE.
    ELSE
      DoSpherical = .true.
    ENDIF
    IF((ID3 .EQ. 'A' .OR. ID3 .EQ. 'X').OR. ID3 .EQ. '*' ) THEN
      Angstrom = .TRUE.
    ELSE
      Angstrom = .FALSE.
    ENDIF

    !johannesfor reading lattice vectors in pbc
   ! IPOS = INDEX(LINE,'PBC')
   ! IF (pbc_check .eq. 'PBC' .or.pbc_check .eq. 'pbc') THEN
   !   setup_pbclatt=.TRUE.
   ! else
   !   setup_pbclatt=.false.
   ! ENDIF
  ENDIF
!  IF(.NOT.DoSpherical)call lsquit('at the moment only spherical harmonic basis functionas have been implemented.',lupri)

END SUBROUTINE READ_LINE4

!> \brief read line 5 in the molecule input file CONTAINING ATOMICCHARGE,NUMBER OF ATOMS,ATOMICBASISSET and then x,y,z coordinates.
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param iprint the printlevel, determining how much output should be generated
!> \param BASISSETLIBRARY the info about the basisset required
!> \param Atomtypes the number of different atomtypes
!> \param if the periodic boundry conditions is to be used
!> \param atombasis is the atombasis keyword used in 1. line
!> \param basis is the basis keyword used in 1. line
!> \param auxbasis is the auxilliary basis given in 1. line
!> \param angstrom is the coordinates given in angstrom or bohr (au)
!> \param molecule the molecule to be built
!> \param PRINTATOMCOORD should the coordinates be printed
!> \param doprint if we should print this to output files
!>
!> THIS IS THE MAIN ROUTINE WHICH READS THE ATOMIC INFORMATION 
!>
!> THE RECIPE:
!> 1. READ LINE5 CONTAINING ATOMICCHARGE,NUMBER OF ATOMS,ATOMICBASISSET,..
!> 2. READS THE X,Y,Z COORDINATES OF THE INDIVIDUAL ATOMS
!> 3. IF THE FILE IS WRITTEN IN ANGSTROM WE CONVERT TO ATOMIC UNIT 
!>
SUBROUTINE READ_GEOMETRY(LUPRI,LUINFO,IPRINT,BASISSETLIBRARY,Atomtypes,dopbc,&
     &ATOMBASIS,BASIS,AUXBASIS,CABSBASIS,JKBASIS,Angstrom,MOLECULE,PRINTATOMCOORD,doprint,&
     &latt_config)
  use ls_util
implicit none
TYPE(MOLECULEINFO) :: MOLECULE
TYPE(BASISSETLIBRARYITEM) :: BASISSETLIBRARY
TYPE(lvec_list_t)  :: latt_config
INTEGER            :: Atomtypes,LUINFO,IPOS,atomnumber
INTEGER            :: I,J,natoms
CHARACTER(len=4)   :: StringFormat
LOGICAL            :: ATOMBASIS,BASIS,Angstrom,dopbc,AUXBASIS,doprint
LOGICAL            :: CABSBASIS,JKBASIS
real(realk)        :: AtomicCharge
CHARACTER(len=80)  :: LINE,Atomicbasisset,Auxbasisset,Cabsbasisset,JKbasisset  
CHARACTER(len=1)   :: CHRXYZ(3)=(/'x','y','z'/)
INTEGER            :: LUPRI,basissetnumber,basissetnumber1,basissetnumber2
INTEGER            :: unique1,unique2,uniqueCharge,IPRINT
LOGICAL            :: PRINTATOMCOORD,pointcharge,phantom
INTEGER :: auxunique2,auxbasissetnumber2,cabsunique2,cabsbasissetnumber2
INTEGER :: jkunique2,jkbasissetnumber2
atomnumber=0
basissetnumber=0
DO I=1,Atomtypes
 CALL READ_LINE5(LUPRI,LUINFO,AtomicCharge,nAtoms,AtomicBasisset,ATOMBASIS,&
      & BASIS,AUXBASIS,AUXBASISSET,CABSBASIS,CABSbasisset,JKBASIS,JKbasisset,pointcharge,phantom)

 IF(ATOMBASIS)THEN
    CALL DETERMINE_UNIQUE_BASIS(BASISSETLIBRARY,ATOMICBASISSET,&
                                             &unique1,basissetnumber)
    IF(unique1 == 0)THEN !found new basisset
       basissetnumber=basissetnumber+1
       IF(basissetnumber .GT. maxBasisSetInLIB) THEN
          WRITE(LUPRI,*)'You use many different basisset (which is okay),&
               & but you need to increase maxBasisSetInLIB in TYPE-DEF.f90&
               & to a number equal to or greater than the number of different&
               & basissets (including auxiliary basissets). At the moment it&
               & is set to ', maxBasisSetInLIB,'Which is clearly not enough for you.&
               & Thomas Kjaergaard'
          CALL LSQUIT('Increase maxBasisSetInLIB in TYPE-DEF.f90',lupri)
       ENDIF
       BASISSETLIBRARY%BASISSETNAME(basissetnumber)=ATOMICBASISSET       
       BASISSETLIBRARY%Charges(basissetnumber,1)=AtomicCharge
       BASISSETLIBRARY%pointcharges(basissetnumber,1)=pointcharge
       BASISSETLIBRARY%phantom(basissetnumber,1)=phantom
       BASISSETLIBRARY%nCharges(basissetnumber)=1 
       basissetnumber1=basissetnumber
       BASISSETLIBRARY%nbasissets=basissetnumber
    ELSE !old basisset
       CALL DETERMINE_UNIQUE_CHARGE(BASISSETLIBRARY,unique1,AtomicCharge,uniqueCharge,pointcharge)
       IF(uniqueCharge /= 0)THEN !found new charge
          IF(uniqueCharge .GT. maxNumberOfChargesinLIB) THEN
             WRITE(LUPRI,*)'You use many different charges (which is okay),&
                  & but you need to increase maxNumberOfChargesinLIB in TYPE-DEF.f90&
                  & to a number equal to or greater than the number of different&
                  & charges in the same basisset. So if you have 18 different atoms&
                  & with 6-31G basis and 15 different atoms with STO-3G you need&
                  & to increase maxNumberOfChargesinLIB to 18 or higher. &
                  & maxNumberOfChargesinLIB is set to ', maxNumberOfChargesinLIB,&
                  &'Which is clearly not enough for you.&
                  & Thomas Kjaergaard'
             CALL LSQUIT('Increase maxBasisSetInLIB in TYPE-DEF.f90',lupri)
          ENDIF

          BASISSETLIBRARY%Charges(unique1,uniqueCharge)=AtomicCharge
          BASISSETLIBRARY%pointcharges(unique1,uniqueCharge)=pointcharge
          BASISSETLIBRARY%phantom(unique1,uniqueCharge)=phantom
          BASISSETLIBRARY%nCharges(unique1)=uniqueCharge
!       ELSE
       ENDIF !old basisset old charge
    ENDIF
    IF(AUXBASIS)THEN
       CALL DETERMINE_UNIQUE_AND_SET_BASIS(LUPRI,BASISSETLIBRARY,AUXBASISSET,&
                &basissetnumber,AtomicCharge,uniqueCharge,auxunique2,pointcharge,phantom)
       auxbasissetnumber2=basissetnumber
    ENDIF
    IF(CABSBASIS)THEN
       IF(.NOT.AUXBASIS)THEN
          CALL LSQUIT('Aux is required When supplying a CABS basis set, even if not used',-1)
       ENDIF
       CALL DETERMINE_UNIQUE_AND_SET_BASIS(LUPRI,BASISSETLIBRARY,CABSBASISSET,&
                &basissetnumber,AtomicCharge,uniqueCharge,cabsunique2,pointcharge,phantom)
       cabsbasissetnumber2=basissetnumber
    ENDIF
    IF(JKBASIS)THEN
       IF(.NOT.CABSBASIS)THEN
          CALL LSQUIT('CABS is required When supplying a JK basis set, even if not used',-1)
       ENDIF
       CALL DETERMINE_UNIQUE_AND_SET_BASIS(LUPRI,BASISSETLIBRARY,JKBASISSET,&
                &basissetnumber,AtomicCharge,uniqueCharge,jkunique2,pointcharge,phantom)
       jkbasissetnumber2=basissetnumber
    ENDIF
 ENDIF


 IF(BASIS)THEN
    unique1=1
    basissetnumber1=1
    unique2=2
    basissetnumber2=2
    IF(AUXBASIS)THEN
       auxunique2 = unique2
       auxbasissetnumber2=basissetnumber2
       unique2=unique2+1
       basissetnumber2=basissetnumber2+1
    ENDIF
    IF(CABSBASIS)THEN
       CABSunique2 = unique2
       CABSbasissetnumber2=basissetnumber2
       unique2=unique2+1
       basissetnumber2=basissetnumber2+1
    ENDIF
    IF(JKBASIS)THEN
       JKunique2 = unique2
       JKbasissetnumber2=basissetnumber2
    ENDIF
    IF(I==1)THEN
       BASISSETLIBRARY%Charges(unique1,1)=AtomicCharge
       BASISSETLIBRARY%pointcharges(unique1,1)=pointcharge
       BASISSETLIBRARY%phantom(unique1,1)=phantom
       BASISSETLIBRARY%nCharges(unique1)=1
       IF(AUXBASIS)THEN
          BASISSETLIBRARY%Charges(auxunique2,1)=AtomicCharge
          BASISSETLIBRARY%pointcharges(auxunique2,1)=pointcharge
          BASISSETLIBRARY%phantom(auxunique2,1)=phantom
          BASISSETLIBRARY%nCharges(auxunique2)=1
       ENDIF
       IF(CABSBASIS)THEN
          BASISSETLIBRARY%Charges(cabsunique2,1)=AtomicCharge
          BASISSETLIBRARY%pointcharges(cabsunique2,1)=pointcharge
          BASISSETLIBRARY%phantom(cabsunique2,1)=phantom
          BASISSETLIBRARY%nCharges(cabsunique2)=1
       ENDIF
       IF(JKBASIS)THEN
          BASISSETLIBRARY%Charges(jkunique2,1)=AtomicCharge
          BASISSETLIBRARY%pointcharges(jkunique2,1)=pointcharge
          BASISSETLIBRARY%phantom(jkunique2,1)=phantom
          BASISSETLIBRARY%nCharges(jkunique2)=1
       ENDIF
    ELSE
       CALL DETERMINE_UNIQUE_CHARGE(BASISSETLIBRARY,unique1,AtomicCharge,uniqueCharge,pointcharge)
       IF(uniqueCharge /= 0)THEN !found new charge
          IF(uniqueCharge .GT. maxNumberOfChargesinLIB) THEN
             WRITE(LUPRI,*)'You use many different atoms/charges (which is okay),&
                  & but you need to increase maxNumberOfChargesinLIB in TYPE-DEF.f90&
                  & to a number equal to or greater than the number of different&
                  & charges. So if you have 18 different atoms&
                  & with 6-31G basis, you need&
                  & to increase maxNumberOfChargesinLIB to 18 or higher. &
                  & maxNumberOfChargesinLIB is set to ', maxNumberOfChargesinLIB,&
                  &'Which is clearly not enough for you.&
                  & Thomas Kjaergaard'
             CALL LSQUIT('Increase maxBasisSetInLIB in TYPE-DEF.f90',lupri)
          ENDIF
          BASISSETLIBRARY%Charges(unique1,uniqueCharge)=AtomicCharge
          BASISSETLIBRARY%pointcharges(unique1,uniqueCharge)=pointcharge
          BASISSETLIBRARY%phantom(unique1,uniqueCharge)=phantom
          BASISSETLIBRARY%nCharges(unique1)=uniqueCharge
          IF(AUXBASIS)THEN
             BASISSETLIBRARY%Charges(auxunique2,uniqueCharge)=AtomicCharge
             BASISSETLIBRARY%pointcharges(auxunique2,uniqueCharge)=pointcharge
             BASISSETLIBRARY%phantom(auxunique2,uniqueCharge)=phantom
             BASISSETLIBRARY%nCharges(auxunique2)=uniqueCharge
          ENDIF
          IF(CABSBASIS)THEN
             BASISSETLIBRARY%Charges(cabsunique2,uniqueCharge)=AtomicCharge
             BASISSETLIBRARY%pointcharges(cabsunique2,uniqueCharge)=pointcharge
             BASISSETLIBRARY%phantom(cabsunique2,uniqueCharge)=phantom
             BASISSETLIBRARY%nCharges(cabsunique2)=uniqueCharge
          ENDIF
          IF(JKBASIS)THEN
             BASISSETLIBRARY%Charges(jkunique2,uniqueCharge)=AtomicCharge
             BASISSETLIBRARY%pointcharges(jkunique2,uniqueCharge)=pointcharge
             BASISSETLIBRARY%phantom(jkunique2,uniqueCharge)=phantom
             BASISSETLIBRARY%nCharges(jkunique2)=uniqueCharge
          ENDIF
       ENDIF
    ENDIF

 ENDIF
 DO J=1,nAtoms
    atomnumber=atomnumber+1
    MOLECULE%ATOM(atomnumber)%phantom = phantom
    MOLECULE%ATOM(atomnumber)%pointCharge = pointCharge
    MOLECULE%ATOM(atomnumber)%Charge=AtomicCharge
    CALL ATTACH_BASISINDEX_AND_BASISLABEL(LUPRI,MOLECULE,&
         &atomnumber,basissetnumber1,auxbasissetnumber2,CABSbasissetnumber2,&
         &JKbasissetnumber2,unique1,auxunique2,CABSunique2,JKunique2,&
         &BASIS,ATOMBASIS,AUXBASIS,CABSBASIS,JKBASIS)
    ! Set atomic numbers, isotopes and masses
    MOLECULE%ATOM(atomnumber)%Isotope = 1
    MOLECULE%ATOM(atomnumber)%Atomic_number = NINT(AtomicCharge)
    MOLECULE%ATOM(atomnumber)%Mass = &
    Isotopes(MOLECULE%ATOM(Atomnumber)%Atomic_number, &
    MOLECULE%ATOM(AtomNumber)%Isotope,'MASS',LUPRI)
    !READ_ATOMCOORD 

    READ (LUINFO, '(a80)') LINE

    IPOS = INDEX(LINE,' ')
    IF (IPOS .LT. 2) THEN
      WRITE (LUPRI,*) 'Atom name must start in first column'
      CALL LSQUIT('Error in placement of atom name. See output',lupri)
    ELSEIF (IPOS .EQ. 2) THEN
      StringFormat = '(A1)'
    ELSEIF (IPOS .EQ. 3) THEN
      StringFormat = '(A2)'
    ELSEIF (IPOS .EQ. 4) THEN
      StringFormat = '(A3)'
    ELSEIF (IPOS .EQ. 5) THEN
      StringFormat = '(A4)'
    ELSE
      WRITE (LUPRI,*) 'Atom name must be less then 4 letters'
      CALL LSQUIT('Error in atom name. See output',lupri)
    END IF

    READ (LINE,StringFormat) MOLECULE%ATOM(atomnumber)%Name

!Read x,y and z coordinate 
    READ (LINE(IPOS:80),*) MOLECULE%ATOM(atomnumber)%CENTER(1), &
    & MOLECULE%ATOM(atomnumber)%CENTER(2), MOLECULE%ATOM(atomnumber)%CENTER(3)

!#if defined(SYS_CRAY) && defined (VAR_NOFREE) && defined (SYS_T3D) && defined (SYS_T90)
!         READ (LINE(IPOS:80),*,ERR=101) &
!                & MOLECULE%ATOM(atomnumber)%CENTER(1)&  !X
!                & MOLECULE%ATOM(atomnumber)%CENTER(2)&  !Y
!                & MOLECULE%ATOM(atomnumber)%CENTER(3)   !Z
!         GOTO 104
!  101    READ (LINE(IPOS:80),'(BN,3F20.0)',ERR=102) & 
!                & MOLECULE%ATOM(atomnumber)%CENTER(1)&  !X
!                & MOLECULE%ATOM(atomnumber)%CENTER(2)&  !Y
!                & MOLECULE%ATOM(atomnumber)%CENTER(3)   !Z
!         GO TO 104
!  102    READ (LINE(IPOS:80),'(BN,3F10.0)',ERR=103)& 
!                & MOLECULE%ATOM(atomnumber)%CENTER(1)&  !X
!                & MOLECULE%ATOM(atomnumber)%CENTER(2)&  !Y
!                & MOLECULE%ATOM(atomnumber)%CENTER(3)   !Z
!         GO TO 104
!  103    CONTINUE
!            WRITE(LUPRI,'(/A,I5/A,I5,A)')&
!     &      ' ERROR: Unable to read Cartesian coordinates of atom no.',&
!     &      atomnumber,' from the following line in the MOLECULE input file:'
!            WRITE(LUPRI,'(A)') LINE
!         CALL QUIT('ERROR reading atomic coordinates in MOLECULE input')
!  104    CONTINUE
!#else
!    CALL FREFRM(LINE,IPOS,IDUMMY,MOLECULE%ATOM(atomnumber)%CENTER(1),'REA',IERR)
!    CALL FREFRM(LINE,IPOS,IDUMMY,MOLECULE%ATOM(atomnumber)%CENTER(2),'REA',IERR)
!    CALL FREFRM(LINE,IPOS,IDUMMY,MOLECULE%ATOM(atomnumber)%CENTER(3),'REA',IERR)
!#endif

    IF (dopbc) THEN
      write(LUPRI,*) 'debug: read atom position'
      write(LUPRI,*) 'PBC not implemented'
!     Here we need to do two things: First, check that
!     the given atom position is within the cell. Second,
!     if the position was given in lattice coordinates we
!     need to transform to standard coordinates, since
!     subsequent calculations assume standard coordinates.
!     ErikT

!     call pbc_transf_atom_coord(NUCIND,CORD(1,NUCIND),atom_coord_lat)
    ENDIF

    IPOS = INDEX(LINE,'Isotope=')
    IF (IPOS .NE. 0) THEN
      IPOS = IPOS + 8
      READ (LINE(IPOS:80),'(I3)') MOLECULE%ATOM(atomnumber)%Isotope
    ELSE
      MOLECULE%ATOM(atomnumber)%Isotope = 1
    END IF

    IF(Angstrom) THEN
      MOLECULE%ATOM(atomnumber)%CENTER(1) = MOLECULE%ATOM(atomnumber)%CENTER(1)/bohr_to_angstrom 
      MOLECULE%ATOM(atomnumber)%CENTER(2) = MOLECULE%ATOM(atomnumber)%CENTER(2)/bohr_to_angstrom 
      MOLECULE%ATOM(atomnumber)%CENTER(3) = MOLECULE%ATOM(atomnumber)%CENTER(3)/bohr_to_angstrom 
    ENDIF

! IF(ABS(CORDinates).GT.CORMAX) THEN
!        WRITE (LUPRI,'(A,E12.6,A/A/A,E12.6)')
!     &    ' Atomic coordinate ',CORD(J,NUCIND),
!     &    ' too large in CNTINP.',
!     &    ' Note: Program is unstable for large coordinates.',
!     &    ' Maximum coordinate value:',CORMAX
!        CALL QUIT('*** ERROR: Atomic coordinate too large in CNTINP')
!ENDIF

!IF(IPRINT) THEN
!  WRITE(LUPRI,'(6X,A4,3F20.15)') NAMN(NUCIND),(CORD(J,NUCIND), J = 1,3)
!ENDIF

!!$C
!!$C     To avoid problems in later stages, we now rewrite the coordinates
!!$C     of the nuclei in traditional Dalton format, with 4 characters for the
!!$C     name of the atom
!!$C
!!$         IF (IPOS .EQ. 0) THEN
!!$            WRITE (MLINE(NMLINE),'(A4,3F20.10,2X,A14)') NAMN(NUCIND),
!!$     &           (CORD(J,NUCIND),J = 1, 3), '                '
!!$         ELSE
!!$            IF (MASSNM .LT. 10) THEN
!!$               WRITE (MLINE(NMLINE),'(A4,3F20.10,A10,I1,A5)') 
!!$     &              NAMN(NUCIND),(CORD(J,NUCIND),J = 1, 3), 
!!$     &              '  Isotope=', MASSNM, '     '
!!$            ELSE IF (MASSNM .LT. 100) THEN
!!$               WRITE (MLINE(NMLINE),'(A4,3F20.10,A10,I2,A4)') 
!!$     &              NAMN(NUCIND),(CORD(J,NUCIND),J = 1, 3), 
!!$     &              '  Isotope=', MASSNM, '    '
!!$            ELSE
!!$               WRITE (MLINE(NMLINE),'(A4,3F20.10,A10,I3,A3)') 
!!$     &              NAMN(NUCIND),(CORD(J,NUCIND),J = 1, 3), 
!!$     &              '  Isotope=', MASSNM, '   '
!!$            END IF
!!$         END IF


!C*****************************************************************************
!C       MULBSI  - basis-set identifier (WK/UniKA/04-11-2002).
!C       CHARGE  - charge of center
!C       NOORBT  - TRUE: no orbitals on this center
!C       GNUEXP  - exponent of Gaussian nuclear charge distribution
!C*****************************************************************************

!         MULBSI(NUCIND) = MAX(MBSI,1)
!         IF (MBSI .GT. 1) THEN
!            CHARGE(NUCIND) = D0
!            IZATOM(NUCIND) = -1
!C           ... code -1 to tell MBSI (value of 0 is floating orbitals)
!            GNUEXP(NUCIND) = D0
!         ELSE
!            CHARGE(NUCIND) = Q
!            IZATOM(NUCIND) = NQ
!C           ... if point charge this is reset outside to -2
!C               to tell this is a point charge /hjaaj Nov 2003
!            GNUEXP(NUCIND) = GEXP
!         END IF
  ENDDO
ENDDO

J=0
DO I = 1,MOLECULE%natoms
   IF(MOLECULE%ATOM(I)%pointCharge)CYCLE
   J=J+1
ENDDO
MOLECULE%natomsNPC = J

IF ((IPRINT .GT. 0).AND.doprint) THEN
   WRITE (LUPRI,'(2X,A,I3)')' Total number of atoms:',MOLECULE%natoms
   IF(MOLECULE%natomsNPC.NE.MOLECULE%natoms)&
        &WRITE (LUPRI,'(2X,A,I3)')' Total number of atoms NOT including PointCharges:',MOLECULE%natomsNPC
ENDIF

IF(IPRINT .GT. -1 .AND. PRINTATOMCOORD.AND.DOPRINT) THEN
   CALL LSHEADER(LUPRI,'Cartesian Coordinates Linsca (au)')
   CALL PRINT_GEOMETRY(MOLECULE,LUPRI)
ENDIF
IF(latt_config%setup_pbclatt) THEN
  !READ lattice vectors
  CALL READ_LATT_VECTORS(LUPRI,LUINFO,latt_config)
ENDIF
END SUBROUTINE READ_GEOMETRY

!> \brief determines if the basisset is unique and set the basissetlibrary accordingly
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param BASISSETLIBRARY the info about the basisset required
!> \param auxbasisset the auxilliary basisset for this atom
!> \param basissetnumber the counter that counts the number of basisset
!> \param AtomicCharge the charge of the atom
!> \param UniqueCharge the index if the charge is unique
!> \param unique2 zero if unique, otherwise the index to oldbasis
SUBROUTINE DETERMINE_UNIQUE_AND_SET_BASIS(LUPRI,BASISSETLIBRARY,AUXBASISSET,&
     &basissetnumber,AtomicCharge,uniqueCharge,unique2,pointcharge,phantom)
  implicit none
  TYPE(BASISSETLIBRARYITEM) :: BASISSETLIBRARY
  real(realk)        :: AtomicCharge
  INTEGER            :: LUPRI,basissetnumber,unique2
  INTEGER            :: uniqueCharge
  CHARACTER(len=80)  :: Auxbasisset
  LOGICAL            :: pointcharge,phantom
  CALL DETERMINE_UNIQUE_BASIS(BASISSETLIBRARY,AUXBASISSET,&
       &unique2,basissetnumber)
  IF(unique2 == 0)THEN !found new basisset
     basissetnumber=basissetnumber+1
     IF(basissetnumber .GT. maxBasisSetInLIB) THEN
        WRITE(LUPRI,*)'You use many different basisset (which is okay),&
             & but you need to increase maxBasisSetInLIB in TYPE-DEF.f90&
             & to a number equal to or greater than the number of different&
             & basissets (including auxiliary basissets). At the moment it&
             & is set to ', maxBasisSetInLIB,'Which is clearly not enough for you.&
             & Thomas Kjaergaard'
        CALL LSQUIT('Increase maxBasisSetInLIB in TYPE-DEF.f90',lupri)
     ENDIF
     BASISSETLIBRARY%BASISSETNAME(basissetnumber)=AUXBASISSET
     BASISSETLIBRARY%Charges(basissetnumber,1)=AtomicCharge
     BASISSETLIBRARY%pointcharges(basissetnumber,1)=pointcharge
     BASISSETLIBRARY%phantom(basissetnumber,1)=phantom
     BASISSETLIBRARY%nCharges(basissetnumber)=1 
     BASISSETLIBRARY%nbasissets=basissetnumber
  ELSE !old basisset
     CALL DETERMINE_UNIQUE_CHARGE(BASISSETLIBRARY,unique2,AtomicCharge,uniqueCharge,pointcharge)
     
     IF(uniqueCharge /= 0)THEN !found new charge
        IF(uniqueCharge .GT. maxNumberOfChargesinLIB) THEN
           WRITE(LUPRI,*)'You use many different charges (which is okay),&
                & but you need to increase maxNumberOfChargesinLIB in TYPE-DEF.f90&
                & to a number equal to or greater than the number of different&
                & charges in the same basisset. So if you have 18 different atoms&
                & with 6-31G basis and 15 different atoms with STO-3G you need&
                & to increase maxNumberOfChargesinLIB to 18 or higher. &
                & maxNumberOfChargesinLIB is set to ', maxNumberOfChargesinLIB,&
                &'Which is clearly not enough for you.&
                & Thomas Kjaergaard'
           CALL LSQUIT('Increase maxBasisSetInLIB in TYPE-DEF.f90',lupri)
        ENDIF
        BASISSETLIBRARY%Charges(unique2,uniqueCharge)=AtomicCharge
        BASISSETLIBRARY%pointcharges(unique2,uniqueCharge)=pointcharge
        BASISSETLIBRARY%phantom(unique2,uniqueCharge)=phantom
        BASISSETLIBRARY%nCharges(unique2)=uniqueCharge
     ELSE
     ENDIF !old basisset old charge
  ENDIF
  
END SUBROUTINE DETERMINE_UNIQUE_AND_SET_BASIS

!> \brief read line 5 in the molecule input file: Read the atom-specific information in new input style
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param luinfo the logical unit number for the molecule input file
!> \param AtomicCharge the charge of the atom
!> \param natoms the number of atoms of this type
!> \param Atomicbasisset the basisset for this atom
!> \param atombasis if ATOMBASIS keyword is given in line 1
!> \param basis if BASIS keyword is given in line 1
!> \param auxbasis is the auxilliary basis given in 1. line
!> \param auxbasisset the auxilliary basisset for this atom
SUBROUTINE READ_LINE5(LUPRI,LUINFO,AtomicCharge,nAtoms,AtomicBasisset,&
     &ATOMBASIS,BASIS,AuxBASIS,Auxbasisset,CABSBASIS,CABSbasisset,JKBASIS,JKbasisset,pointcharge,phantom)
implicit none
real(realk)        :: AtomicCharge
INTEGER            :: IPOS,IPOS2,IPOS3,nAtoms,LUINFO
CHARACTER(len=80)  :: TEMPLINE,Atomicbasisset,Auxbasisset,CABSbasisset,JKbasisset
LOGICAL            :: ATOMBASIS,BASIS,AUXBASIS,pointcharge,phantom,CABSBASIS,JKBASIS
CHARACTER(len=5)   :: StringFormat
INTEGER            :: LUPRI,ios,I

READ (LUINFO, '(a80)') TEMPLINE
IPOS = INDEX(TEMPLINE,'Cha')
IF (IPOS .NE. 0 ) THEN
   IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
   IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 7)) THEN
      WRITE (LUPRI,'(2X,A40)') 'Incorrect input for atomic charge'
      WRITE (LUPRI,'(2X,A40)') 'Format is "Charge=?"'
      CALL LSQUIT('Incorrect input for atomic charge',lupri)
   ELSE
      READ (TEMPLINE((IPOS+IPOS2):),*) AtomicCharge
      IPOS = INDEX(TEMPLINE,'Ato')
      IF (IPOS .NE. 0) THEN
         IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
         IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
            WRITE (LUPRI,'(2X,A40)') 'READ_LINE5: Incorrect input for # of atoms'
            WRITE (LUPRI,'(2X,A40)') 'READ_LINE5: Format is "Atoms=?"'
            CALL LSQUIT('READ_LINE5: Incorrect input for # of atoms',lupri)
         ELSE
            READ (TEMPLINE((IPOS+IPOS2):),*) nAtoms
         ENDIF
      ENDIF
   ENDIF
ELSE !OLD INPUT STYLE
   READ (TEMPLINE,'(BN,F10.0,I5)',IOSTAT=ios) AtomicCharge, nAtoms
   IF(IOS .NE. 0) THEN
      WRITE (LUPRI,'(2X,A50)') 'READ_LINE5: OLD: Error in the determination of the number of atoms'
      WRITE (LUPRI,'(2X,A50)') "READ_LINE5: OLD: Correct input structure is: Atoms=???"
      CALL LSQUIT('READ_LINE5: OLD: Error in determining the number of atoms',lupri)
   ENDIF
ENDIF

IPOS = INDEX(TEMPLINE,'pointcharge')
pointcharge = .FALSE.
IF (IPOS .NE. 0) THEN
   pointcharge = .TRUE.
ENDIF

IPOS = INDEX(TEMPLINE,'phantom')
phantom = .FALSE.
IF (IPOS .NE. 0) THEN
   phantom = .TRUE.
ENDIF

!     Multiple basis sets used?
!         
!  IF (LMULBS) THEN
!    IPOS = INDEX(TEMPLINE,'Set')
!    IF (IPOS .NE. 0) THEN
!      IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
!      IF (IPOS2 .EQ. 0 .OR. (IPOS2. GT. 5)) THEN
!        WRITE (2,*) 'Incorrect input for # of basis sets'
!        WRITE (2,*) 'Format is "Sets=?"'
!        CALL QUIT('Incorrect input for # of basis sets')
!      ELSE
!        READ (TEMPLINE((IPOS+IPOS2):),*) MBSI
!      END IF
!    END IF
!  END IF
! Read in basis set information
! Integral ignored for now
IF (.NOT. (BASIS .OR. ATOMBASIS)) THEN
   IPOS = INDEX(TEMPLINE,'Blo')
   IF (IPOS .NE. 0) THEN
      IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
      IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 7)) THEN
         WRITE (LUPRI,'(2X,A40)') 'Incorrect input for # of integral blocks'
         WRITE (LUPRI,'(2X,A40)') 'Format is "Blocks=? ? ?"'
         CALL LSQUIT('Incorrect input for # of integral blocks',lupri)
      ELSE
         !        ISTART = IPOS + IPOS2
         !        CALL FREFRM(TEMPLINE,ISTART,IQM,DUMMY,'INT',IERR)
         !        DO IQMLOP = 1, IQM
         !          CALL FREFRM(TEMPLINE,ISTART,JCO(IQMLOP),DUMMY,'INT',IERR)
         !        END DO
         WRITE(LUPRI,*) 'integral blocks not implemented'
         CALL LSQUIT('# of integral blocks not implemented',lupri)       
      END IF
   END IF
ENDIF

IF (ATOMBASIS) THEN
   AUXBASIS=.FALSE.
   CABSBASIS=.FALSE.
   JKBASIS=.FALSE.
   IPOS = INDEX(TEMPLINE,'Bas')
   IF (IPOS .NE. 0) THEN
      IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
      IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
         WRITE (LUPRI,'(2X,A40)') 'Incorrect input for choice of atomic basis set'
         WRITE (LUPRI,'(2X,A40)') 'Format is "Basis=? ? ?"'
         CALL LSQUIT('Incorrect input for choice of atomic basis set',lupri)
      ELSE
         IPOS3 = INDEX(TEMPLINE((IPOS+IPOS2):),' ')
         IF (IPOS3 .LT. 10) THEN
            WRITE (StringFormat,'(A2,I1,A1,1X)') '(A',IPOS3 - 1,')'
         ELSE
            WRITE (StringFormat,'(A2,I2,A1)') '(A',(IPOS3 - 1),')'
         ENDIF
         READ (TEMPLINE((IPOS + IPOS2):),StringFormat) AtomicBasisset
      ENDIF
   ELSE
      IF(pointcharge)THEN
         DO I=1,80
            AtomicBasisset(I:I)=' '
         ENDDO
         AtomicBasisset(1:11)='pointcharge'
      ELSE
         WRITE (LUPRI,*) 'ATOMBASIS selected, but no atomic basis&
              & set specified for one atom type'
         CALL LSQUIT( 'ATOMBASIS selected, but no atomic basis set &
              &specified for one atom type',lupri)
      ENDIF
   ENDIF

   !Auxiliary basisset
   IPOS = INDEX(TEMPLINE,'Aux')
   IF (IPOS .NE. 0) THEN
      AUXBASIS=.TRUE.
      IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
      IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 4)) THEN
         WRITE (LUPRI,*) 'Incorrect input for choice of auxiliary atomic basis set'
         WRITE (LUPRI,*) ' Format is "Aux=? ? ?"'
         CALL LSQUIT('Incorrect input for choice of auxiliary atomic basis set',lupri)
      ELSE
         IPOS3 = INDEX(TEMPLINE((IPOS+IPOS2):),' ')
         IF (IPOS3 .LT. 10) THEN
            WRITE (StringFormat,'(A2,I1,A1,1X)') '(A',IPOS3 - 1,')'
         ELSE
            WRITE (StringFormat,'(A2,I2,A1)') '(A',(IPOS3 - 1),')'
         ENDIF
         READ (TEMPLINE((IPOS + IPOS2):),StringFormat) AuxBasisset
      ENDIF
   ENDIF

   !CABS basisset - Complementary Auxiliary basis set
   IPOS = INDEX(TEMPLINE,'CABS')
   IF (IPOS .NE. 0) THEN
      CABSBASIS=.TRUE.
      IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
      IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 5)) THEN
         WRITE (LUPRI,*) 'Incorrect input for choice of Complementary auxiliary atomic basis set'
         WRITE (LUPRI,*) ' Format is "CABS=? ? ?"'
         CALL LSQUIT('Incorrect input for choice of Complementary auxiliary atomic basis set',lupri)
      ELSE
         IPOS3 = INDEX(TEMPLINE((IPOS+IPOS2):),' ')
         IF (IPOS3 .LT. 10) THEN
            WRITE (StringFormat,'(A2,I1,A1,1X)') '(A',IPOS3 - 1,')'
         ELSE
            WRITE (StringFormat,'(A2,I2,A1)') '(A',(IPOS3 - 1),')'
         ENDIF
         READ (TEMPLINE((IPOS + IPOS2):),StringFormat) CABSBasisset
      ENDIF
   ENDIF
   
   !JK basisset
   IPOS = INDEX(TEMPLINE,'JK')
   IF (IPOS .NE. 0) THEN
      JKBASIS=.TRUE.
      IPOS2 = INDEX(TEMPLINE(IPOS:),'=')
      IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 3)) THEN
         WRITE (LUPRI,*) 'Incorrect input for choice of JK atomic basis set'
         WRITE (LUPRI,*) ' Format is "JK=? ? ?"'
         CALL LSQUIT('Incorrect input for choice of JK atomic basis set',lupri)
      ELSE
         IPOS3 = INDEX(TEMPLINE((IPOS+IPOS2):),' ')
         IF (IPOS3 .LT. 10) THEN
            WRITE (StringFormat,'(A2,I1,A1,1X)') '(A',IPOS3 - 1,')'
         ELSE
            WRITE (StringFormat,'(A2,I2,A1)') '(A',(IPOS3 - 1),')'
         ENDIF
         READ (TEMPLINE((IPOS + IPOS2):),StringFormat) JKBasisset
      ENDIF
   ENDIF
ENDIF
END SUBROUTINE READ_LINE5

!> \brief determine if the basis is unique
!> \author T. Kjaergaard
!> \date 2010
!> \param BASISSETLIBRARY the info about the basisset required
!> \param BASISSET the basis that may or may not be unique
!> \param unique output determining if the basis is unique
!> \param nsets the number of different basissets in library
!>
!> if the basiset is not contained in basissetlibrary unique is set to
!> zero otherwise it is set to the index of the basisset in the basisset 
!> library 
!>
SUBROUTINE DETERMINE_UNIQUE_BASIS(BASISSETLIBRARY,BASISSET,unique,nsets)
implicit none
TYPE(BASISSETLIBRARYITEM) :: BASISSETLIBRARY
CHARACTER(len=80)         :: basisset
INTEGER                   :: unique,nsets,I

unique=0

IF(nsets/=0)THEN
   DO I=1,nsets 
      IF(BASISSETLIBRARY%BASISSETNAME(I) .EQ. basisset)THEN
         unique=I
      ENDIF
   ENDDO
ENDIF

END SUBROUTINE DETERMINE_UNIQUE_BASIS

!> \brief determine if the charge is unique
!> \author T. Kjaergaard
!> \date 2010
!> \param BASISSETLIBRARY the info about the basisset required
!> \param basisnumber the number in the basislibrary
!> \param atomicCharge that may or may not be unique
!> \param uniqueCharge 0 if unique
!>
!> if the charge is unique in the basissetlibrary
!> uniqueCharge is set to zero otherwise it is set to the index
!>
SUBROUTINE DETERMINE_UNIQUE_CHARGE(BASISSETLIBRARY,basisnumber,&
     &atomicCharge,uniqueCharge,pointcharge)
implicit none
TYPE(BASISSETLIBRARYITEM) :: BASISSETLIBRARY
INTEGER               :: basisnumber,uniqueCharge,I
real(realk)           :: atomicCharge
LOGICAL               :: FOUND,pointcharge
uniqueCharge=0
IF(.NOT.pointcharge)THEN
   FOUND=.FALSE.
   DO I=1,BASISSETLIBRARY%nCharges(basisnumber)
      IF(BASISSETLIBRARY%pointcharges(basisnumber,I))CYCLE 
      IF(ABS(BASISSETLIBRARY%Charges(basisnumber,I)-atomicCharge).LT.1.0E-12_realk)THEN
         FOUND=.TRUE.
      ENDIF
   ENDDO
   IF(.NOT. FOUND) uniqueCharge=BASISSETLIBRARY%nCharges(basisnumber)+1
ENDIF
END SUBROUTINE DETERMINE_UNIQUE_CHARGE

!> \brief ATTACH_BASISINDEX AND BASISLABEL to the atom
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param moleclue the molecule structure
!> \param I the atom index in the molecule structure
!> \param basissetnumber1 the basisindex
!> \param auxbasissetnumber2 the auxilliary basisindex
!> \param CABSbasissetnumber2 the complementary auxilliary basisindex
!> \param JKbasissetnumber2 the JK auxilliary basisindex
!> \param unique1 basis index 
!> \param auxunique2 auxilliary basis index
!> \param CABSunique2 complementary auxilliary basis index
!> \param JKunique2 JK auxilliary basis index
!> \param basis if BASIS keyword is given in line 1
!> \param atombasis if ATOMBASIS keyword is given in line 1
!> \param auxbasis if the auxilliary basis given in 1. line
!> \param CABSbasis if the complementary auxilliary basis given in 1. line
!> \param JKbasis if the JK auxilliary basis given in 1. line
!>
!> IF YOU ADD MORE BASISSET TO MOLECULE FILE YOU NEED TO 
!> INCREASE THE SIZE IN THE ATOM DERIVED TYPE IN TYPE-DEF.f90. 
!>
SUBROUTINE ATTACH_BASISINDEX_AND_BASISLABEL(LUPRI,MOLECULE,&
     &I,basissetnumber1,auxbasissetnumber2,CABSbasissetnumber2,&
     &JKbasissetnumber2,unique1,auxunique2,CABSunique2,JKunique2,&
     &BASIS,ATOMBASIS,AUXBASIS,CABSBASIS,JKBASIS)
implicit none
TYPE(MOLECULEINFO) :: MOLECULE 
LOGICAL            :: BASIS,ATOMBASIS,AUXBASIS,CABSBASIS,JKBASIS
INTEGER            :: unique1,I,LUPRI,auxbasissetnumber2 !I=atomnumber
INTEGER            :: basissetnumber1,CABSbasissetnumber2,JKbasissetnumber2
INTEGER            :: auxunique2,CABSunique2,JKunique2

IF(BASIS)THEN
   MOLECULE%ATOM(I)%basisindex(1)=1
   MOLECULE%ATOM(I)%basislabel(1)='REGULAR  '
   MOLECULE%ATOM(I)%nbasis=1
   IF(AUXBASIS)THEN
      MOLECULE%ATOM(I)%basisindex(2)=2
      MOLECULE%ATOM(I)%basislabel(2)='AUXILIARY'
      MOLECULE%ATOM(I)%nbasis=2
   ENDIF
   IF(CABSBASIS)THEN
      MOLECULE%ATOM(I)%basisindex(3)=3
      MOLECULE%ATOM(I)%basislabel(3)='CABS     '
      MOLECULE%ATOM(I)%nbasis=3
   ENDIF
   IF(JKBASIS)THEN
      MOLECULE%ATOM(I)%basisindex(4)=4
      MOLECULE%ATOM(I)%basislabel(4)='JKAUX    '
      MOLECULE%ATOM(I)%nbasis=4
   ENDIF   
   
ELSEIF(ATOMBASIS)THEN
   IF(unique1 /= 0)THEN
      MOLECULE%ATOM(I)%basisindex(1)=unique1
      MOLECULE%ATOM(I)%basislabel(1)='REGULAR  '
   ELSE
      MOLECULE%ATOM(I)%basisindex(1)=basissetnumber1
      MOLECULE%ATOM(I)%basislabel(1)='REGULAR  '
   ENDIF
   MOLECULE%ATOM(I)%nbasis=1
   IF(AUXBASIS)THEN
      IF(AUXunique2 /= 0)THEN
         MOLECULE%ATOM(I)%basisindex(2)=AUXunique2
         MOLECULE%ATOM(I)%basislabel(2)='AUXILIARY'
      ELSE
         MOLECULE%ATOM(I)%basisindex(2)=AUXbasissetnumber2
         MOLECULE%ATOM(I)%basislabel(2)='AUXILIARY'
      ENDIF
      MOLECULE%ATOM(I)%nbasis=2
   ENDIF
   IF(CABSBASIS)THEN
      IF(CABSunique2 /= 0)THEN
         MOLECULE%ATOM(I)%basisindex(3)=CABSunique2
         MOLECULE%ATOM(I)%basislabel(3)='CABS     '
      ELSE
         MOLECULE%ATOM(I)%basisindex(3)=CABSbasissetnumber2
         MOLECULE%ATOM(I)%basislabel(3)='CABS     '
      ENDIF
      MOLECULE%ATOM(I)%nbasis=3
   ENDIF
   IF(JKBASIS)THEN
      IF(JKunique2 /= 0)THEN
         MOLECULE%ATOM(I)%basisindex(4)=JKunique2
         MOLECULE%ATOM(I)%basislabel(4)='JKAUX    '
      ELSE
         MOLECULE%ATOM(I)%basisindex(4)=JKbasissetnumber2
         MOLECULE%ATOM(I)%basislabel(4)='JKAUX    '
      ENDIF
      MOLECULE%ATOM(I)%nbasis=4
   ENDIF
ELSE
   CALL LSQUIT('Something wrong in determining basisetsets in READ_GEOMETRY',-1)
ENDIF

END SUBROUTINE ATTACH_BASISINDEX_AND_BASISLABEL

!> \brief build distinct atoms 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri the logical unit number for the output file
!> \param ndatoms the number of distinct atoms
!> \param natoms the total number of atoms
!> \param natomlist for a given distinct atom a atomindex i full molecule
!> \param uatomtype for a given distinct atom the atomtype in full molecule
!> \param ls contains molecule and basis structure
!> \param iprint the printlevel, determining how much output should be generate
SUBROUTINE BUILD_DISTINCT_ATOMS(LUPRI,NDATOMS,NATOMS,NATOMLIST,UATOMTYPE,LS,IPRINT)
use BUILDBASISSET
use molecule_module
IMPLICIT NONE

INTEGER  :: LUPRI,NDATOMS,NATOMS,IPRINT,I
TYPE(LSITEM)        :: LS
INTEGER             :: CHARGE(NDATOMS),NATOMLIST(NDATOMS),UATOMTYPE(NDATOMS)
INTEGER             :: NDATOMS2
Character(len=80)   :: NAME(NDATOMS)
INTEGER             :: itype,TYPEN(NDATOMS),icharge
INTEGER             :: R
INTEGER,pointer :: NEWATOM(:)

R = LS%INPUT%BASIS%REGULAR%Labelindex
call mem_alloc(NEWATOM,LS%INPUT%BASIS%REGULAR%natomtypes)
NEWATOM=0
NDATOMS2 = 0
IF(R.EQ. 0)THEN
   DO I=1,LS%INPUT%MOLECULE%nATOMs
      IF(LS%INPUT%MOLECULE%ATOM(I)%pointcharge)CYCLE
      IF(LS%INPUT%MOLECULE%ATOM(I)%phantom)CYCLE
      icharge = INT(LS%INPUT%MOLECULE%ATOM(I)%Charge)
      itype = LS%INPUT%BASIS%REGULAR%Chargeindex(icharge)
      IF(NEWATOM(itype).EQ. 0)THEN
         NEWATOM(itype)=1

         NDATOMS2 = NDATOMS2+1
         TYPEN(NDATOMS2) = itype
         CHARGE(NDATOMS2) = icharge
         NAME(NDATOMS2)= LS%INPUT%BASIS%REGULAR%ATOMTYPE(itype)%NAME
         NATOMLIST(NDATOMS2) = I
         UATOMTYPE(NDATOMS2) = itype
!         UATOMLIST(I)=NDATOMS2
      ENDIF
   ENDDO
ELSE
   DO I=1,LS%INPUT%MOLECULE%nATOMs
      IF(LS%INPUT%MOLECULE%ATOM(I)%pointcharge)CYCLE
      IF(LS%INPUT%MOLECULE%ATOM(I)%phantom)CYCLE
      icharge = INT(LS%INPUT%MOLECULE%ATOM(I)%Charge)
      itype = LS%INPUT%MOLECULE%ATOM(I)%IDtype(R)
      IF(NEWATOM(itype).EQ. 0)THEN
         NEWATOM(itype)=1

         NDATOMS2 = NDATOMS2+1
         TYPEN(NDATOMS2) = itype
         CHARGE(NDATOMS2) = icharge
         NAME(NDATOMS2)= LS%INPUT%BASIS%REGULAR%ATOMTYPE(itype)%NAME
         NATOMLIST(NDATOMS2) = I
         UATOMTYPE(NDATOMS2) = itype
!         UATOMLIST(I)=NDATOMS2        
      ENDIF
   ENDDO
ENDIF

IF(NDATOMS .LT. NDATOMS2)CALL LSQUIT('DISTINCT_ATOMS and BUILD_DISTINCT_ATOMS not consistent',lupri)
NDATOMS=NDATOMS2
call mem_dealloc(NEWATOM)

END SUBROUTINE BUILD_DISTINCT_ATOMS
!!$FUNCTION ISOMAS(CHARGE,MASSNumber)
!!$!  Function to switch from mass number to isotope number sorted
!!$!  according to abundance
!!$implicit none
!!$INTEGER    ::  IORD,I,QMASS,MASSNumber
!!$IORD = 0
!!$DO I = 1, 5
!!$  QMASS = DISOTP(ICHARG,I,'MASS')
!!$  IF (ANINT(QMASS) .EQ. MASSNM) IORD = I
!!$END DO
!!$IF (IORD .EQ. 0) THEN
!!$  WRITE (LUPRI,'(/A,I4,A,I4)') 'ERROR: unknown mass',MASSNM,
!!$  &' for atom with charge ',ICHARG
!!$  CALL QUIT('Unknown mass for chosen atomic charge')
!!$ELSE
!!$  ISOMAS = IORD
!!$END IF
!!$END FUNCTION ISOMAS
!!$
!!$FUNCTION DISOTP(IATOM,ISOTOP,TYPE)
!!$! NOTE: Isotopes are sorted according to abundance,
!!$! i.e. DISOTP(IATOM,1,TYPE) will return the most abundant
!!$! isotope etc.
!!$! Proton mass and electron charge: 1986 CODATA Recommended Values
!!$!
!!$! Nuclear masses:
!!$!   A. H. Wapstra and K. Bos,
!!$!   Atomic Data and Nuclear Tables 19 (1977) 177
!!$!
!!$! Abundancies:
!!$!   Handbook of Chemistry and Physics, 73rd Edition
!!$!        
!!$! Nuclear moments and spins:
!!$!   P. Raghavan,
!!$!   Atomic Data and Nuclear Data Tables 42 (1989) 189
!!$!
!!$! Quadrupole moments:
!!$!   P.Pykkoe and J.Li
!!$!   Report HUKI 1-92
!!$!   Updated  to P.Pyykkoe, Mol.Phys. (2001) by K.Ruud, Aug.2001
!!$!
!!$!   Nuclear masses, Abundancies, nuclear moments, spins 
!!$!   and quadrupole moments for Z= 55 to Z = 86:
!!$!
!!$!   I. Mills, T. Cvitas, K. Homann, N. Kallay, and K. Kuchitsu
!!$!   Quantities, Units and Symbols in Physical Chemistry
!!$!   (IUPAC, Blackwell Scientific, Oxford, 1988)
!!$!
!!$CHARACTER*(*) TYPE
!!$INTEGER           :: MAXISO,MAXCHR
!!$real(realk)       :: DATNUC(5,6,86),PMASS,EMASS,XFAMU,THRESH,DMP
!!$real(realk)       :: FOUR,ZERO
!!$PMASS = 1.007276470E0_realk
!!$EMASS = 9.10938188E-31_realk
!!$XFAMU = 1822.88848E0_realk
!!$THRESH = 1.0E-10_realk
!!$DMP = PMASS*XFAMU*EMASS
!!$MAXISO = 6
!!$MAXCHR = 86
!!$FOUR = 4.0E0_realk
!!$ZERO = 0.0E0_realk

!!$! H - Ne
!!$! ======
!!$     DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=0,10) /
!!$! Dummy:
!!$     &  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,  0.000000E0_realk,
!!$! H:
!!$     &  1.007825E0_realk, 99.985000E0_realk,   .500000E0_realk,  2.792847E0_realk,   .000000E0_realk,
!!$     &  2.014102E0_realk,   .015000E0_realk,  1.000000E0_realk,   .857438E0_realk,   .002860E0_realk,
!!$     &  3.016049E0_realk,   .000000E0_realk,   .500000E0_realk,  2.978962E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! He:
!!$     &  4.002603E0_realk, 99.999870E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &  3.016029E0_realk,   .000130E0_realk,   .500000E0_realk, -2.127625E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! Li:
!!$     &  7.016005E0_realk, 92.500000E0_realk,  1.500000E0_realk,  3.256427E0_realk,  -.040100E0_realk,
!!$     &  6.015123E0_realk,  7.500000E0_realk,  1.000000E0_realk,   .822047E0_realk,  -.000808E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! Be:
!!$     &  9.012183E0_realk,100.000000E0_realk,  1.500000E0_realk, -1.177800E0_realk,   .052880E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! B:
!!$     & 11.009305E0_realk, 80.100000E0_realk,  1.500000E0_realk,  2.688649E0_realk,   .040590E0_realk,
!!$     & 10.012938E0_realk, 19.900000E0_realk,  3.000000E0_realk,  1.800645E0_realk,   .084590E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$!  C:
!!$     & 12.000000E0_realk, 98.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 13.003355E0_realk,  1.100000E0_realk,   .500000E0_realk,   .702412E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! N:
!!$     & 14.003074E0_realk, 99.630000E0_realk,  1.000000E0_realk,   .403761E0_realk,   .020440E0_realk,
!!$     & 15.000109E0_realk,   .370000E0_realk,   .500000E0_realk,  -.283189E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 12.000000E0_realk, 98.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! O:
!!$     & 15.994915E0_realk, 99.760000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 17.999159E0_realk,   .200000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 16.999131E0_realk,   .040000E0_realk,  2.500000E0_realk, -1.893790E0_realk,  -.025580E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 14.003074E0_realk, 99.630000E0_realk,  1.000000E0_realk,   .403761E0_realk,   .020200E0_realk,
!!$! F:
!!$     & 18.998403E0_realk,100.000000E0_realk,   .500000E0_realk,  2.628868E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 15.994915E0_realk, 99.760000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! Ne:
!!$     & 19.992439E0_realk, 90.480000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 21.991384E0_realk,  9.250000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 20.993845E0_realk,   .270000E0_realk,  1.500000E0_realk,  -.661797E0_realk,   .101550E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!!$! Na - Ar
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=11,18) /
!!$! Na:
!!$     & 22.989770E0_realk,100.000000E0_realk,  1.500000E0_realk,  2.217656E0_realk,   .104000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! Mg:
!!$     & 23.985045E0_realk, 78.990000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 25.982595E0_realk, 11.010000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 24.985839E0_realk, 10.000000E0_realk,  2.500000E0_realk,  -.855450E0_realk,   .199400E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! Al:
!!$     & 26.981541E0_realk,100.000000E0_realk,  2.500000E0_realk,  3.641507E0_realk,   .146600E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! Si:
!!$     & 27.976928E0_realk, 92.230000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 28.976496E0_realk,  4.670000E0_realk,   .500000E0_realk,  -.555290E0_realk,   .000000E0_realk,
!!$     & 29.973772E0_realk,  3.100000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! P:
!!$     & 30.973763E0_realk,100.000000E0_realk,   .500000E0_realk,  1.131600E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! S:
!!$     & 31.972072E0_realk, 95.020000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 33.967868E0_realk,  4.210000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 32.971459E0_realk,   .750000E0_realk,  1.500000E0_realk,   .643821E0_realk,  -.067800E0_realk,
!!$     & 35.967079E0_realk,   .020000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! Cl:
!!$     & 34.968853E0_realk, 75.770000E0_realk,  1.500000E0_realk,   .821874E0_realk,  -.081650E0_realk,
!!$     & 36.965903E0_realk, 24.230000E0_realk,  1.500000E0_realk,   .684124E0_realk,  -.064350E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! Ar:
!!$     & 39.962383E0_realk, 99.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 35.967546E0_realk,   .337000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 37.962732E0_realk,   .063000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!!$! K - Ca
!!$! ======
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=19,20) /
!!$! K:
!!$     & 38.963708E0_realk, 93.258100E0_realk,  1.500000E0_realk,   .391507E0_realk,   .058500E0_realk,
!!$     & 40.961825E0_realk,  6.730200E0_realk,  1.500000E0_realk,   .214893E0_realk,   .071100E0_realk,
!!$     & 39.963999E0_realk,   .011700E0_realk,  4.000000E0_realk, -1.298100E0_realk,  -.073000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$! Ca:
!!$     & 39.962591E0_realk, 96.941000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 43.955485E0_realk,  2.086000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 41.958622E0_realk,   .647000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 47.952532E0_realk,   .187000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 42.958770E0_realk,   .135000E0_realk,  3.500000E0_realk, -1.317643E0_realk,  -.040800E0_realk,
!!$     & 45.953689E0_realk,   .004000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!!$! Sc - Zn
!!$! =======
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=21,30) /
!!$! Sc:
!!$     & 44.955914E0_realk,100.000000E0_realk,  3.500000E0_realk,  4.756487E0_realk,  -.220000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Ti:
!!$C
!!$     & 47.947947E0_realk, 73.800000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 45.952633E0_realk,  8.000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 46.951765E0_realk,  7.300000E0_realk,  2.500000E0_realk,  -.788480E0_realk,   .302000E0_realk,
!!$     & 48.947871E0_realk,  5.500000E0_realk,  3.500000E0_realk, -1.104170E0_realk,   .247000E0_realk,
!!$     & 49.944786E0_realk,  5.400000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     V:
!!$C
!!$     & 50.943963E0_realk, 99.750000E0_realk,  3.500000E0_realk,  5.148706E0_realk,  -.052000E0_realk,
!!$     & 49.947161E0_realk,   .250000E0_realk,  6.000000E0_realk,  3.345689E0_realk,   .210000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Cr:
!!$C
!!$     & 51.940510E0_realk, 83.790000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 52.940651E0_realk,  9.500000E0_realk,  1.500000E0_realk,  -.474540E0_realk,  -.150000E0_realk,
!!$     & 49.946046E0_realk,  4.345000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 53.938882E0_realk,  2.365000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Mn:
!!$C
!!$     & 54.938046E0_realk,100.000000E0_realk,  2.500000E0_realk,  3.468719E0_realk,   .330000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Fe:
!!$C
!!$     & 55.934939E0_realk, 91.720000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 53.939612E0_realk,  5.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 56.935396E0_realk,  2.100000E0_realk,   .500000E0_realk,   .090623E0_realk,   .000000E0_realk,
!!$     & 57.933278E0_realk,   .280000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Co:
!!$C
!!$     & 58.933198E0_realk,100.000000E0_realk,  3.500000E0_realk,  4.627000E0_realk,   .420000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Ni:
!!$C
!!$     & 57.935347E0_realk, 68.077000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 59.930789E0_realk, 26.223000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 61.928346E0_realk,  3.634000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 60.931059E0_realk,  1.140000E0_realk,  1.500000E0_realk,  -.750020E0_realk,   .162000E0_realk,
!!$     & 63.927968E0_realk,  0.926000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Cu:
!!$C
!!$     & 62.929599E0_realk, 69.170000E0_realk,  1.500000E0_realk,  2.227206E0_realk,  -.220000E0_realk,
!!$     & 64.927792E0_realk, 30.830000E0_realk,  1.500000E0_realk,  2.381610E0_realk,  -.204000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Zn:
!!$C
!!$     & 63.929145E0_realk, 48.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 65.926035E0_realk, 27.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 67.924846E0_realk, 18.800000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 66.927129E0_realk,  4.100000E0_realk,  2.500000E0_realk,   .875479E0_realk,   .150000E0_realk,
!!$     & 69.925325E0_realk,   .600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!!$C
!!$C     Ga - Kr
!!$C     =======
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=31,36) /
!!$C
!!$C     Ga:
!!$C
!!$     & 68.925581E0_realk, 60.108000E0_realk,  1.500000E0_realk,  2.016589E0_realk,   .171000E0_realk,
!!$     & 70.924701E0_realk, 39.892000E0_realk,  1.500000E0_realk,  2.562266E0_realk,   .107000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Ge:
!!$C
!!$     & 73.921179E0_realk, 35.940000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 71.922080E0_realk, 27.660000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 69.924250E0_realk, 21.240000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 72.923464E0_realk,  7.720000E0_realk,  4.500000E0_realk,  -.879468E0_realk,  -.196000E0_realk,
!!$     & 75.921403E0_realk,  7.440000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     As:
!!$C
!!$     & 74.921596E0_realk,100.000000E0_realk,  1.500000E0_realk,  1.439475E0_realk,   .314000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Se:
!!$C
!!$     & 79.916521E0_realk, 49.610000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 77.917304E0_realk, 23.770000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 75.919207E0_realk,  9.360000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 81.916709E0_realk,  8.740000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 76.919908E0_realk,  7.630000E0_realk,   .500000E0_realk,   .535042E0_realk,   .000000E0_realk,
!!$     & 73.922477E0_realk,  0.890000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Br:
!!$C
!!$     & 78.918336E0_realk, 50.690000E0_realk,  1.500000E0_realk,  2.106400E0_realk,   .313000E0_realk,
!!$     & 80.916290E0_realk, 49.310000E0_realk,  1.500000E0_realk,  2.270562E0_realk,   .261500E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Kr:
!!$C
!!$     & 83.911506E0_realk, 57.000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 85.910614E0_realk, 17.300000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 81.913483E0_realk, 11.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 82.914134E0_realk, 11.500000E0_realk,  4.500000E0_realk,  -.970669E0_realk,   .259000E0_realk,
!!$     & 79.916375E0_realk,  2.250000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 77.920397E0_realk,  0.350000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!!$C
!!$C     Rb:
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=37,45) /
!!$C
!!$     & 84.911800E0_realk, 72.170000E0_realk,  2.500000E0_realk,  1.353352E0_realk,   .276000E0_realk,
!!$     & 86.909184E0_realk, 27.830000E0_realk,  1.500000E0_realk,  2.751818E0_realk,   .133500E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Sr:
!!$C
!!$     & 87.905625E0_realk, 82.580000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 85.909273E0_realk,  9.860000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 86.908890E0_realk,  7.000000E0_realk,  4.500000E0_realk, -1.093603E0_realk,   .335000E0_realk,
!!$     & 83.913428E0_realk,   .560000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Y:
!!$C
!!$     & 88.905856E0_realk,100.000000E0_realk,   .500000E0_realk,  -.137415E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Zr:
!!$C
!!$     & 89.904708E0_realk, 51.450000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 93.906319E0_realk, 17.380000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 91.905039E0_realk, 17.150000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 90.905644E0_realk, 11.220000E0_realk,  2.500000E0_realk, -1.303620E0_realk,  -.176000E0_realk,
!!$     & 95.908272E0_realk,  2.800000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Nb:
!!$C
!!$     & 92.906378E0_realk,100.000000E0_realk,  4.500000E0_realk,  6.170500E0_realk,  -.320000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Mo:
!!$C
!!$     & 97.905405E0_realk, 24.130000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 95.904676E0_realk, 16.680000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 94.905838E0_realk, 15.920000E0_realk,  2.500000E0_realk,  -.914200E0_realk,  -.022000E0_realk,
!!$     & 93.905086E0_realk, 14.840000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 99.907473E0_realk,  9.630000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 96.906018E0_realk,  9.550000E0_realk,  2.500000E0_realk,  -.933500E0_realk,  0.255000E0_realk,
!!$C
!!$C     Tc:
!!$C
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Ru:
!!$C
!!$     &101.904348E0_realk, 31.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &103.905422E0_realk, 18.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &100.905581E0_realk, 17.100000E0_realk,  2.500000E0_realk,  -.718800E0_realk,   .457000E0_realk,
!!$     & 98.905937E0_realk, 12.700000E0_realk,  2.500000E0_realk,  -.641300E0_realk,   .079000E0_realk,
!!$     & 99.904218E0_realk, 12.600000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     & 95.907596E0_realk,  5.540000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Rh:
!!$C
!!$     &102.905503E0_realk,100.000000E0_realk,   .500000E0_realk,  -.088400E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=46,54) /
!!$C
!!$C     Pd:
!!$C
!!$     &105.903475E0_realk, 27.330000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &107.903894E0_realk, 26.460000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &104.905075E0_realk, 22.330000E0_realk,  2.500000E0_realk,  -.642000E0_realk,   .660000E0_realk,
!!$     &109.905169E0_realk, 11.720000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &103.904026E0_realk, 11.140000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &101.905609E0_realk,  1.020000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Ag:
!!$C
!!$     &106.905095E0_realk, 51.839000E0_realk,   .500000E0_realk,  -.113570E0_realk,   .000000E0_realk,
!!$     &108.904754E0_realk, 48.161000E0_realk,   .500000E0_realk,   .130563E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Cd:
!!$C
!!$     &113.903361E0_realk, 28.730000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &111.902761E0_realk, 24.130000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &110.904182E0_realk, 12.800000E0_realk,   .500000E0_realk,  -.594886E0_realk,   .000000E0_realk,
!!$     &109.903007E0_realk, 12.490000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &112.904401E0_realk, 12.220000E0_realk,   .500000E0_realk,  -.622301E0_realk,   .000000E0_realk,
!!$     &115.904758E0_realk,  7.490000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     In:
!!$C
!!$     &114.903875E0_realk, 95.700000E0_realk,  4.500000E0_realk,  5.540800E0_realk,   .810000E0_realk,
!!$     &112.904056E0_realk,  4.300000E0_realk,  4.500000E0_realk,  5.528900E0_realk,   .799000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Sn:
!!$C
!!$     &119.902199E0_realk, 32.590000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &117.901607E0_realk, 24.220000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &115.901744E0_realk, 14.530000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &118.903310E0_realk,  8.580000E0_realk,   .500000E0_realk, -1.047280E0_realk,   .000000E0_realk,
!!$     &116.902954E0_realk,  7.680000E0_realk,   .500000E0_realk, -1.001040E0_realk,   .000000E0_realk,
!!$     &123.905271E0_realk,  5.790000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Sb:
!!$C
!!$     &120.903824E0_realk, 57.360000E0_realk,  2.500000E0_realk,  3.363400E0_realk,  -.360000E0_realk,
!!$     &122.904222E0_realk, 42.640000E0_realk,  3.500000E0_realk,  2.549800E0_realk,  -.490000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Te:
!!$C
!!$     &129.906229E0_realk, 33.870000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &127.904464E0_realk, 31.700000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &125.903310E0_realk, 18.930000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &124.904435E0_realk,  7.120000E0_realk,   .500000E0_realk,  -.888505E0_realk,   .000000E0_realk,
!!$     &123.902825E0_realk,  4.790000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &121.903055E0_realk,  2.590000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     I:
!!$C
!!$     &126.904477E0_realk,100.000000E0_realk,  2.500000E0_realk,  2.813273E0_realk,  -.710000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$C
!!$C     Xe:
!!$C
!!$     &131.904148E0_realk, 26.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &128.904780E0_realk, 26.400000E0_realk,   .500000E0_realk,  -.777976E0_realk,   .000000E0_realk,
!!$     &130.905076E0_realk, 21.200000E0_realk,  1.500000E0_realk,   .691862E0_realk,  -.114000E0_realk,
!!$     &133.905395E0_realk, 10.400000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &135.907219E0_realk,  8.900000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk,
!!$     &129.903510E0_realk,  4.100000E0_realk,   .000000E0_realk,   .000000E0_realk,   .000000E0_realk/
!!$C
!!$C
!!$C
!!$C
!!$C    Cs - Rn*
!!$C    =======
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=55,64) /
!!$C
!!$C     Cs:
!!$C
!!$     &132.905429E0_realk,100.000000E0_realk,  3.50000E0_realk,   2.582025E0_realk,  -0.00343E0_realk,
!!$     &  0.00000E0_realk,   0.000000E0_realk,  0.00000E0_realk,   0.000000E0_realk,   0.00000E0_realk,
!!$     &  0.00000E0_realk,   0.000000E0_realk,  0.00000E0_realk,   0.000000E0_realk,   0.00000E0_realk,
!!$     &  0.00000E0_realk,   0.000000E0_realk,  0.00000E0_realk,   0.000000E0_realk,   0.00000E0_realk,
!!$     &  0.00000E0_realk,   0.000000E0_realk,  0.00000E0_realk,   0.000000E0_realk,   0.00000E0_realk,
!!$     &  0.00000E0_realk,   0.000000E0_realk,  0.00000E0_realk,   0.000000E0_realk,   0.00000E0_realk,
!!$C
!!$C     Ba:
!!$C   
!!$     &137.905232E0_realk, 71.70000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &136.905812E0_realk, 11.23000E0_realk,   1.50000E0_realk,  0.937365E0_realk,   0.245000E0_realk,
!!$     &135.904553E0_realk,  7.85400E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &134.905665E0_realk,  6.59200E0_realk,   1.50000E0_realk,  0.837943E0_realk,   0.160000E0_realk,
!!$     &133.904486E0_realk,  2.41700E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &131.905042E0_realk,  0.10100E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     La:
!!$C
!!$     &138.906347E0_realk, 99.90980E0_realk,   3.50000E0_realk,  2.7830455E0_realk,  0.200000E0_realk,
!!$     &137.907105E0_realk,  0.09020E0_realk,   5.00000E0_realk,  3.7136460E0_realk,  0.450000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Ce
!!$C
!!$     &139.905433E0_realk, 88.48000E0_realk,   0.00000E0_realk,  0.000000E0_realk,    0.00000E0_realk,
!!$     &141.909241E0_realk, 11.08000E0_realk,   0.00000E0_realk,  0.000000E0_realk,    0.00000E0_realk,
!!$     &137.905985E0_realk,  0.25000E0_realk,   0.00000E0_realk,  0.000000E0_realk,    0.00000E0_realk,
!!$     &135.907140E0_realk,  0.19000E0_realk,   0.00000E0_realk,  0.000000E0_realk,    0.00000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Pr:
!!$C
!!$     &140.907647E0_realk,100.0000E0_realk,    2.50000E0_realk,  4.275400E0_realk,  -0.058900E0_realk,
!!$     & 0.000000E0_realk,  0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     & 0.000000E0_realk,  0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     & 0.000000E0_realk,  0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     & 0.000000E0_realk,  0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     & 0.000000E0_realk,  0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Nd:
!!$C
!!$     &141.907719E0_realk, 27.130000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &143.910083E0_realk, 23.800000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &145.913113E0_realk, 17.190000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &142.909810E0_realk, 12.180000E0_realk,   3.50000E0_realk, -0.065000E0_realk,  -0.630000E0_realk,
!!$     &144.912570E0_realk,  8.300000E0_realk,   3.50000E0_realk, -1.065000E0_realk,  -0.330000E0_realk,
!!$     &147.916889E0_realk,  5.760000E0_realk,   0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Pm:
!!$C
!!$     &144.912743E0_realk,100.000000E0_realk,   2.50000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Sm:
!!$C
!!$     &151.919728E0_realk, 26.700000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &153.922205E0_realk, 22.700000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &146.914894E0_realk, 15.000000E0_realk,   3.500000E0_realk,-0.814800E0_realk,  -0.259000E0_realk,
!!$     &148.917180E0_realk, 13.800000E0_realk,   3.500000E0_realk,-0.671700E0_realk,   0.075000E0_realk,
!!$     &147.914819E0_realk, 11.300000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &149.917273E0_realk,  7.400000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Eu:
!!$C
!!$     &152.921225E0_realk, 52.200000E0_realk,   2.500000E0_realk, 1.533000E0_realk,   2.412000E0_realk,
!!$     &150.919702E0_realk, 47.800000E0_realk,   2.500000E0_realk, 3.471700E0_realk,   0.903000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.00000E0_realk,   0.00000E0_realk,    0.00000E0_realk,  0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Gd:
!!$C
!!$     &157.924019E0_realk, 24.840000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &159.927049E0_realk, 21.860000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &155.922118E0_realk, 20.470000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &156.923956E0_realk, 15.650000E0_realk,   1.500000E0_realk,-0.337260E0_realk,   1.350000E0_realk,
!!$     &154.922618E0_realk, 14.800000E0_realk,   1.500000E0_realk,-0.257230E0_realk,   1.270000E0_realk,
!!$     &153.920861E0_realk,  2.180000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk/
!!$C
!!$C     Tb:
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=65,74) /
!!$C
!!$     &158.925342E0_realk,100.000000E0_realk,   1.500000E0_realk, 2.014000E0_realk,   1.432000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Dy:
!!$C
!!$     &163.929171E0_realk, 28.200000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &161.926795E0_realk, 25.500000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &162.928728E0_realk, 24.900000E0_realk,   2.500000E0_realk, 0.672600E0_realk,   2.648000E0_realk,
!!$     &160.926930E0_realk, 18.900000E0_realk,   2.500000E0_realk,-0.480300E0_realk,   2.507000E0_realk,
!!$     &159.925193E0_realk,  2.340000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &157.924277E0_realk,  0.100000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Ho:
!!$C
!!$     &164.930319E0_realk,100.000000E0_realk,   3.500000E0_realk, 4.173000E0_realk,   3.580000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Er:
!!$C
!!$     &165.930290E0_realk, 33.600000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &167.932368E0_realk, 26.800000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &166.932368E0_realk, 22.950000E0_realk,   3.500000E0_realk,-0.563850E0_realk,   3.565000E0_realk,
!!$     &169.935461E0_realk, 14.900000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &163.929198E0_realk,  1.610000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &161.928775E0_realk,  0.140000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Tm:
!!$C
!!$     &168.934212E0_realk,100.000000E0_realk,   0.500000E0_realk,-0.231600E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Yb:
!!$C
!!$     &173.938859E0_realk, 31.800000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &171.936378E0_realk, 21.900000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &172.938208E0_realk, 16.120000E0_realk,   2.500000E0_realk,-0.679890E0_realk,   2.800000E0_realk,
!!$     &170.936323E0_realk, 14.300000E0_realk,   0.500000E0_realk, 0.493670E0_realk,   0.000000E0_realk,
!!$     &175.942564E0_realk, 12.700000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &169.934759E0_realk,  3.050000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Lu:
!!$C
!!$     &174.940770E0_realk, 97.410000E0_realk,   3.500000E0_realk, 2.232700E0_realk,   3.490000E0_realk,
!!$     &175.942679E0_realk,  2.590000E0_realk,   7.000000E0_realk, 3.169200E0_realk,   4.970000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Hf:
!!$C
!!$     &179.9465457E0_realk,35.100000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &177.943696E0_realk, 27.297000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &176.943217E0_realk, 18.606000E0_realk,   3.500000E0_realk, 0.793500E0_realk,   3.365000E0_realk,
!!$     &178.9458122E0_realk,13.629000E0_realk,   4.500000E0_realk,-0.640900E0_realk,   3.793000E0_realk,
!!$     &175.941406E0_realk,  5.206000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &173.940044E0_realk,  0.162000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Ta:
!!$C
!!$     &180.947992E0_realk, 99.988000E0_realk,   3.500000E0_realk, 2.370500E0_realk,   3.170000E0_realk,
!!$     &179.947462E0_realk,  0.012000E0_realk,   8.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     W:
!!$C
!!$     &183.950928E0_realk, 30.670000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &185.954357E0_realk, 28.600000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &181.948202E0_realk, 26.300000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &182.950928E0_realk, 14.300000E0_realk,   0.500000E0_realk, 0.11778476,   0.000000E0_realk,
!!$     &179.947462E0_realk,  0.162000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk/
!!$C
!!$C     Re:
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=75,80) /
!!$C
!!$     &186.955744E0_realk, 62.600000E0_realk,   2.500000E0_realk, 3.219700E0_realk,   2.070000E0_realk,
!!$     &184.952951E0_realk, 37.400000E0_realk,   2.500000E0_realk, 3.187100E0_realk,   2.180000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Os:
!!$C
!!$     &191.961467E0_realk, 41.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &189.958436E0_realk, 26.400000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &188.958436E0_realk, 16.100000E0_realk,   1.500000E0_realk, 0.659933E0_realk,   0.856000E0_realk,
!!$     &187.955830E0_realk, 13.300000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &186.955741E0_realk,  1.600000E0_realk,   0.500000E0_realk, 0.06465189,   0.000000E0_realk,
!!$     &185.953830E0_realk,  1.580000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Ir:
!!$C
!!$     &192.962917E0_realk, 62.600000E0_realk,   1.500000E0_realk, 0.163700E0_realk,   0.751000E0_realk,
!!$     &190.960584E0_realk, 37.400000E0_realk,   1.500000E0_realk, 0.150700E0_realk,   0.816000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Pt:
!!$C
!!$     &194.964766E0_realk, 33.800000E0_realk,   0.500000E0_realk, 0.609520E0_realk,   0.000000E0_realk,
!!$     &193.962655E0_realk, 32.900000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &195.964926E0_realk, 25.300000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &197.967869E0_realk,  7.200000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &191.961019E0_realk,  0.790000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &189.959917E0_realk,  0.010000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Au:
!!$C   
!!$     &196.966543E0_realk,100.000000E0_realk,   1.500000E0_realk, 0.148158E0_realk,   0.547000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Hg:
!!$C
!!$     &201.970617E0_realk, 29.860000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &199.968300E0_realk, 23.100000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &198.968254E0_realk, 16.870000E0_realk,   0.500000E0_realk, 0.50588549,   0.000000E0_realk,
!!$     &200.970277E0_realk, 13.180000E0_realk,   1.500000E0_realk,-0.5602257 ,   0.386000E0_realk,
!!$     &197.966743E0_realk,  9.970000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &203.973467E0_realk,  6.870000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk/
!!$C
!!$C     Tl:
!!$C
!!$      DATA (((DATNUC(I,J,K),I=1,5),J=1,MAXISO),K=81,86) /
!!$C
!!$     &204.974401E0_realk, 70.476000E0_realk,   0.500000E0_realk, 1.63831461E0_realk, 0.000000E0_realk,
!!$     &202.972320E0_realk, 29.524000E0_realk,   0.500000E0_realk, 1.62225787E0_realk, 0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Pb:
!!$C   
!!$     &207.976627E0_realk, 52.400000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &205.975872E0_realk, 24.100000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &206.975872E0_realk, 22.100000E0_realk,   0.500000E0_realk, 0.582583E0_realk,   0.000000E0_realk,
!!$     &203.973020E0_realk,  1.400000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Bi:
!!$C
!!$     &208.980374E0_realk,100.000000E0_realk,   4.500000E0_realk, 4.110600E0_realk,  -0.516000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Po:
!!$C
!!$     &208.982404E0_realk,  0.000000E0_realk,   0.500000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     At:
!!$C
!!$     &209.987126E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$C
!!$C     Rn:
!!$C
!!$     &222.017571E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk,
!!$     &  0.000000E0_realk,  0.000000E0_realk,   0.000000E0_realk, 0.000000E0_realk,   0.000000E0_realk/

!!$IF (ISOTOP .GT. 6) THEN
!!$  WRITE (LUPRI,'(//,A,2(/,A,I5),A)')'ISOTOP too large in DISOTP. ',&
!!$  &'Input value: ',ISOTOP,' Maximum value:',6,' Program cannot continue.'
!!$  CALL QUIT('MAXISO exceeded in DISOTP')
!!$END IF
!!$IF (IATOM .GT. 86) THEN
!!$  WRITE (LUPRI,'(//,A,2(/,A,I5),A)')' IATOM too large in DISOTP. ',&
!!$  &'Input value: ',IATOM,' Maximum value:',86,' Program cannot continue.'
!!$  CALL QUIT('MAXCHR exceeded in DISOTP')
!!$END IF
!!$
!!$IF (IATOM .LE. 0) THEN
!!$!  This is a floating orbital, a point charge,
!!$!   or an auxiliary basis set /Mar 2004 hjaaj
!!$  DISOTP = D0
!!$ELSE IF (TYPE .EQ. 'MASS') THEN
!!$  DISOTP = DATNUC(1,ISOTOP,IATOM)
!!$ELSE IF (TYPE .EQ. 'A') THEN
!!$  DISOTP = NINT(DATNUC(1,ISOTOP,IATOM))
!!$ELSE IF (TYPE .EQ. 'ABUNDANCE') THEN
!!$  DISOTP = DATNUC(LUPRI,ISOTOP,IATOM)
!!$ELSE IF (TYPE .EQ. 'SPIN') THEN
!!$  DISOTP = DATNUC(3,ISOTOP,IATOM)
!!$ELSE IF (TYPE .EQ. 'MMOM') THEN
!!$  DISOTP = DATNUC(4,ISOTOP,IATOM)
!!$ELSE IF (TYPE .EQ. 'GVAL') THEN
!!$  SPIN = DATNUC(3,ISOTOP,IATOM)
!!$  IF (SPIN .GT. THRESH) THEN
!!$    DISOTP = DATNUC(4,ISOTOP,IATOM)/SPIN
!!$  ELSE
!!$    DISOTP = 0.0E0_realk
!!$  END IF
!!$ELSE IF (TYPE .EQ. 'LARMOR') THEN
!!$  SPIN = DATNUC(3,ISOTOP,IATOM)
!!$  IF (SPIN .GT. THRESH) THEN
!!$    DISOTP = ABS(ECHARGE*DATNUC(4,ISOTOP,IATOM)/(D4*PI*SPIN*DMP))
!!$  ELSE
!!$    DISOTP = 0.0E0_realk
!!$  END IF
!!$ELSE IF (TYPE .EQ. 'QMOM') THEN
!!$  DISOTP = DATNUC(5,ISOTOP,IATOM)
!!$ELSE IF (TYPE .EQ. 'NEUTRONS') THEN
!!$  DISOTP = FLOAT(NINT(DATNUC(1,ISOTOP,IATOM)-IATOM))
!!$ELSE
!!$  WRITE (LUPRI,'(//,3A,/,A)')'Keyword ',TYPE,' unknown in DISOTP. ',&
!!$  &'Program cannot continue.'
!!$  CALL QUIT('Illegal keyword in DISOTP')
!!$END IF
!!$
!!$END FUNCTION DISOTP

END MODULE