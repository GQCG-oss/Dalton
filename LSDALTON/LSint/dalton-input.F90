!> @file
!> Contains module that reads the dalton input and initializes lsitem
MODULE DALTONINFO
use precision
use Integralparameters
use files, only: lsopen,lsclose
use typedeftype, only: LSITEM, integralconfig, DALTONINPUT,&
     & BASISSETLIBRARYITEM
use typedef, only: typedef_init_setting, getNbasis, PRINT_MOLECULEINFO,&
     & PRINT_MOLECULE_AND_BASIS, typedef_set_default_setting, &
     & typedef_free_setting, print_basissetlibrary
use lattice_type, only: lvec_list_t
#ifdef MOD_UNRELEASED
 use lattice_vectors, only: pbc_setup_default
#endif
use basis_type, only: copy_basissetinfo, free_basissetinfo
use basis_typetype,only: nullifyBasis
use io, only: io_init, io_free
use molecule_type, only: free_moleculeinfo
use READMOLEFILE, only: read_molfile_and_build_molecule
use BuildBasisSet, only: Build_BASIS
use lstiming, only: lstimer
use molecule_module, only: build_fragment
use screen_mod, only: screen_init, screen_free
private
public ::  ls_init, ls_free, build_ccfragmentlsitem, dalton_finalize
CONTAINS 
!> \brief initiate the lsitem structure which contain all info about integralevaluation schemes, molecule and basisset.
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ls_init(ls,lupri,luerr,nbast,integral,dodft,doDEC,doprint)
implicit none
!> contains everything needed for integralevaluation (incl. molecule,basis)
TYPE(LSITEM) :: ls
!> the logical unit number for the output file
Integer,intent(in)  :: lupri
!> the logical unit number for the error file
Integer,intent(in)  :: luerr
!> the number of basisfunctions
integer,intent(out) :: nbast
!> a collection of logicals read from the LSDALTON.INP file
type(integralconfig) :: integral 
!> Should we print stuff to output file
Logical      :: doprint
!> is this a dft calculation or a HF calculation
Logical      :: dodft
!> is this a DEC calculation or a HF calculation
Logical      :: doDEC

ls%lupri = lupri
ls%luerr = luerr
ls%optlevel = 3 
call init_AO_parameters
CALL dalton_init(ls%input,lupri,luerr,nbast,integral,dodft,doDEC,doprint)
CALL typedef_init_setting(ls%setting)
call screen_init()
CALL typedef_set_default_setting(ls%setting,ls%input)

END SUBROUTINE ls_init

!> \brief frees the lsitem structure.
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE ls_free(ls)
implicit none
!> contains everything needed for integralevaluation (incl. molecule,basis)
TYPE(LSITEM) :: ls
CALL dalton_finalize(ls%input,ls%lupri,ls%luerr)
CALL typedef_free_setting(ls%setting)
call screen_free()
END SUBROUTINE ls_free

!> \brief build the daltoninput structure associated with the lsitem 
!> \author T. Kjaergaard
!> \date 2010
!>
!> build the daltoninput structure associated with the lsitem, which means 
!> initialize the IOitem (filenames) read MOLECULE.INP file and build molecule
!> and finally build the basis. These informations should not be changed and
!> the name ls%input indicate that it is input information, while
!> ls%setting is the current info about molecule,basis, etc.
!>
SUBROUTINE dalton_init(intinp,LUPRI,LUERR,nbast,integral,dodft,doDEC,doprint)
#ifdef VAR_MPI
use infpar_module
#endif
implicit none
!> contains info about input from LSDALTON.INP(only integral info) and MOLECULE.INP
TYPE(DALTONINPUT)    :: intinp
!> the logical unit number for the output file
integer,intent(in)   :: LUPRI
!> the logical unit number for the error file
integer,intent(in)   :: LUERR
!> number of basisfunctions
integer, intent(out) :: nbast
!> information about the integral evaluation info read from LSDALTON.INP
type(integralconfig)      :: integral
!> is it a DFT run
logical              :: dodft
!> is it a DEC run
logical              :: doDEC
!> should we print to output file
LOGICAL              :: doprint
TYPE(lvec_list_t)    :: latt_config
!
TYPE(BASISSETLIBRARYITEM) :: LIBRARY
integer              :: LUCME,IDUMMY,I
CHARACTER(len=80)    :: BASISSETNAME
CHARACTER(len=9)     :: BASISLABEL
real(realk)          :: Tim1,Tim2,TIME_BUILD_BASIS,TIME_READ_DALTON
real(realk)          :: TIME_AUXBUILD_BASIS,TIME_READ_MOL,TSTART
Integer              :: numAtoms, numNodes
logical              :: lopen
LUCME=-1
!CAREFUL

intinp%nfock = 0 !number of fock matrices calculated

NULLIFY(intinp%MOLECULE)
NULLIFY(intinp%BASIS)
ALLOCATE(intinp%MOLECULE)
ALLOCATE(intinp%BASIS)

IF (doprint) THEN
   !CALL PRINT_INTRO(LUPRI)
   CALL LS_TSTAMP('Start simulation',LUPRI)
   CALL PRINT_MOL_FILE(LUPRI)
   CALL PRINT_DALTON_FILE(LUPRI)
ENDIF

!Initialize IO
CALL io_init(intinp%IO)
!*************************************************
!*
!*  READ THE LINSCA INTEGRAL SECTION OF THE DALTON
!*  INPUT-FILE 
!*
!*************************************************
CALL LSTIMER('START',TIM1,TIM2,LUPRI)

intinp%DALTON = integral

IF(intinp%DALTON%TIMINGS) CALL LSTIMER('READ DALTONFILE',TIM1,TIM2,LUPRI)
!*************************************************
!*
!*  READ THE MOLECULE.INP FILE AND BUILD THE MOLECULE 
!*  STRUCTURE
!*
!*************************************************
#ifdef MOD_UNRELEASED
  call pbc_setup_default(latt_config)
#endif
CALL READ_MOLFILE_AND_BUILD_MOLECULE(LUPRI,intinp%MOLECULE,LIBRARY,doprint,&
     & intinp%dalton%molprint,intinp%dalton%DoSpherical,intinp%dalton%auxbasis,&
     & intinp%dalton%CABSbasis,intinp%dalton%JKbasis,latt_config)

integral%nelectrons = intinp%MOLECULE%nelectrons 
integral%molcharge = INT(intinp%MOLECULE%charge)
numAtoms = intinp%MOLECULE%nAtoms
#ifdef VAR_MPI
  intinp%node         = infpar%mynum
  intinp%numNodes     = infpar%nodtot
  intinp%numFragments = min(numAtoms,intinp%numNodes)
#else
  intinp%node         = 0
  intinp%numNodes     = 1
  intinp%numFragments = 1
#endif
IF (doprint) THEN
  IF (intinp%numFragments.GT. 1) WRITE(LUPRI,*) 'Integrals calculated using ',intinp%numFragments, ' fragments'
ENDIF

IF(intinp%DALTON%TIMINGS) CALL LSTIMER('READ MOLFILE',TIM1,TIM2,LUPRI)

!*************************************************
!*
!*  CALL LINSCF WITH DEFAULT BASIS AS DEFINED IN INPUT
!*
!*************************************************

IF (doprint) THEN
  WRITE(LUPRI,'(2X,A)')' '
  WRITE(LUPRI,'(2X,A)')'CALLING BUILD BASIS WITH DEFAULT REGULAR BASIS'
  WRITE(LUPRI,'(2X,A)')' '
ENDIF
BASISLABEL='REGULAR  '
IF(doDEC)intinp%DALTON%NOFAMILY=.TRUE.
!this is due to the ordering of basisfunctions in the batches and the 
!minAObatch sizes calculation. 
!IF(doDEC)intinp%DALTON%NOSEGMENT=.TRUE. 
IF(doDEC)intinp%DALTON%NOGCINTEGRALTRANSFORM=.TRUE. 
!as we perform calculations of explicit 4 center integrals we cannot
!at this time perform the calculation in the input basis and transform to GC
!but this could be done I think. 
!added complications if the 2int screening matrix should be calculated. 
!it needs to be possible for a single integral to change basis in the integral call

CALL Build_basis(LUPRI,intinp%DALTON%BASPRINT,&
     &intinp%MOLECULE,intinp%BASIS%REGULAR,LIBRARY,&
     &BASISLABEL,intinp%DALTON%UNCONT,intinp%DALTON%NOSEGMENT,&
     &doprint,intinp%DALTON%DOSPHERICAL)
nbast = getNbasis(AORdefault,Contractedinttype,intinp%MOLECULE,LUPRI)

IF(intinp%DALTON%TIMINGS) CALL LSTIMER('BUILD BASIS',TIM1,TIM2,LUPRI)

IF (doprint) THEN
  CALL PRINT_BASISSETLIBRARY(LUPRI,LIBRARY)
ENDIF

IF(intinp%DALTON%AUXBASIS)THEN

   IF (doprint) THEN
      WRITE(LUPRI,'(2X,A)')' '
      WRITE(LUPRI,'(2X,A)')'CALLING BUILD BASIS WITH DEFAULT AUXILIARY BASIS'
      WRITE(LUPRI,'(2X,A)')' '
   ENDIF
   BASISLABEL='AUXILIARY'
   CALL Build_basis(LUPRI,intinp%DALTON%BASPRINT,&
     &intinp%MOLECULE,intinp%BASIS%AUXILIARY,LIBRARY,&
     &BASISLABEL,intinp%DALTON%UNCONT,intinp%DALTON%NOSEGMENT,&
     &doprint,intinp%DALTON%DOSPHERICAL)
   IF(intinp%DALTON%TIMINGS) CALL LSTIMER('BUILD AUXBASIS',TIM1,TIM2,LUPRI)
   IF (doprint) THEN
      CALL PRINT_MOLECULE_AND_BASIS(LUPRI,intinp%MOLECULE,intinp%BASIS%AUXILIARY)
   ENDIF
   
ELSE
intinp%BASIS%AUXILIARY%natomtypes = 0
intinp%BASIS%AUXILIARY%nChargeindex = 0  
intinp%BASIS%AUXILIARY%labelindex = 0  
nullify(intinp%BASIS%AUXILIARY%Chargeindex)
ENDIF

IF(intinp%DALTON%CABSBASIS)THEN
   IF (doprint) THEN
      WRITE(LUPRI,'(2X,A)')' '
      WRITE(LUPRI,'(2X,A)')'CALLING BUILD BASIS WITH COMPLEMENTARY AUXILIARY BASIS'
      WRITE(LUPRI,'(2X,A)')' '
   ENDIF
   BASISLABEL='CABS     '
   CALL Build_basis(LUPRI,intinp%DALTON%BASPRINT,&
     &intinp%MOLECULE,intinp%BASIS%CABS,LIBRARY,&
     &BASISLABEL,intinp%DALTON%UNCONT,intinp%DALTON%NOSEGMENT,&
     &doprint,intinp%DALTON%DOSPHERICAL)
   IF(intinp%DALTON%TIMINGS) CALL LSTIMER('BUILD CABSBASIS',TIM1,TIM2,LUPRI)
   IF (doprint) THEN
      CALL PRINT_MOLECULE_AND_BASIS(LUPRI,intinp%MOLECULE,intinp%BASIS%CABS)
   ENDIF
ELSE
   intinp%BASIS%CABS%natomtypes = 0
   intinp%BASIS%CABS%nChargeindex = 0  
   intinp%BASIS%CABS%labelindex = 0  
   nullify(intinp%BASIS%CABS%Chargeindex)
ENDIF

IF(intinp%DALTON%JKBASIS)THEN
   IF (doprint) THEN
      WRITE(LUPRI,'(2X,A)')' '
      WRITE(LUPRI,'(2X,A)')'CALLING BUILD BASIS WITH COMPLEMENTARY AUXILIARY BASIS'
      WRITE(LUPRI,'(2X,A)')' '
   ENDIF
   BASISLABEL='JKAUX    '
   CALL Build_basis(LUPRI,intinp%DALTON%BASPRINT,&
     &intinp%MOLECULE,intinp%BASIS%JK,LIBRARY,&
     &BASISLABEL,intinp%DALTON%UNCONT,intinp%DALTON%NOSEGMENT,&
     &doprint,intinp%DALTON%DOSPHERICAL)
   IF(intinp%DALTON%TIMINGS) CALL LSTIMER('BUILD JKAUXBASIS',TIM1,TIM2,LUPRI)
   IF (doprint) THEN
      CALL PRINT_MOLECULE_AND_BASIS(LUPRI,intinp%MOLECULE,intinp%BASIS%JK)
   ENDIF
ELSE
   intinp%BASIS%JK%natomtypes = 0
   intinp%BASIS%JK%nChargeindex = 0  
   intinp%BASIS%JK%labelindex = 0  
   nullify(intinp%BASIS%JK%Chargeindex)
ENDIF

intinp%BASIS%VALENCE%natomtypes = 0
intinp%BASIS%VALENCE%nChargeindex = 0  
intinp%BASIS%VALENCE%labelindex = 0  
nullify(intinp%BASIS%VALENCE%Chargeindex)

nullify(intinp%BASIS%GCtrans%ATOMTYPE)
intinp%BASIS%GCtrans%natomtypes = 0
intinp%BASIS%GCtrans%labelindex = 0
intinp%BASIS%GCtrans%nChargeindex = 0
intinp%BASIS%GCtransAlloc = .FALSE.
IF (doprint) THEN
   CALL PRINT_MOLECULEINFO(LUPRI,intinp%MOLECULE,intinp%BASIS,intinp%DALTON%MOLPRINT)
   CALL PRINT_MOLECULE_AND_BASIS(LUPRI,intinp%MOLECULE,intinp%BASIS%REGULAR)
ENDIF

intinp%DO_DFT = dodft
!CALL TSTAMP(' ',LUPRI)
!CALL LSCLOSE(LUCME,'KEEP')
!CALL LSCLOSE(LUMOLDEN,'KEEP')

END SUBROUTINE dalton_init

!> \brief finalize the daltoninput structure
!> \author T. Kjaergaard
!> \date 2010
!>
!> free molecule info
!> free basis
!> free io
!>
SUBROUTINE dalton_finalize(intinp,LUPRI,LUERR)
implicit none
!> contains info about input from LSDALTON.INP(only integral info) and MOLECULE.INP
TYPE(DALTONINPUT)    :: intinp
!> the logical unit number for the output file
integer   :: LUPRI
!> the logical unit number for the error file
integer   :: LUERR

call free_Moleculeinfo(intinp%MOLECULE)
IF(intinp%BASIS%REGULAR%natomtypes .NE. 0)THEN
   call free_basissetinfo(intinp%BASIS%REGULAR)
ENDIF
IF(intinp%BASIS%AUXILIARY%natomtypes .NE. 0)THEN
   call free_basissetinfo(intinp%BASIS%AUXILIARY)
ENDIF
IF(intinp%BASIS%CABS%natomtypes .NE. 0)THEN
   call free_basissetinfo(intinp%BASIS%CABS)
ENDIF
IF(intinp%BASIS%JK%natomtypes .NE. 0)THEN
   call free_basissetinfo(intinp%BASIS%JK)
ENDIF
IF(intinp%BASIS%VALENCE%natomtypes .NE. 0)THEN
   call free_basissetinfo(intinp%BASIS%VALENCE)
ENDIF
IF(intinp%BASIS%GCtransAlloc)THEN
   IF(intinp%BASIS%GCtrans%natomtypes .NE. 0)THEN
      call free_basissetinfo(intinp%BASIS%GCtrans)
   ENDIF
ENDIF
DEALLOCATE(intinp%MOLECULE)
DEALLOCATE(intinp%BASIS)
call io_free(intinp%io)

END SUBROUTINE dalton_finalize

!!$!> \brief free the daltoninput structure
!!$!> \author T. Kjaergaard
!!$!> \date 2010
!!$SUBROUTINE dalton_free(intinp)
!!$implicit none
!!$!> contains info about input from LSDALTON.INP(only integral info) and MOLECULE.INP
!!$TYPE(DALTONINPUT) :: intinp
!!$
!!$call free_Moleculeinfo(intinp%MOLECULE)
!!$call free_basissetinfo(intinp%BASIS%REGULAR)
!!$IF(intinp%DALTON%AUXBASIS)THEN
!!$   call free_basissetinfo(intinp%BASIS%AUXILIARY)
!!$ENDIF
!!$DEALLOCATE(intinp%MOLECULE)
!!$DEALLOCATE(intinp%BASIS)
!!$NULLIFY(intinp%MOLECULE)
!!$NULLIFY(intinp%BASIS)
!!$call io_free(intinp%io)
!!$
!!$END SUBROUTINE dalton_free
 
!> \brief print dalton file
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE PRINT_DALTON_FILE(LUPRI)
!> the logical unit number for the output file
implicit none
INTEGER            :: LUCMD !Logical unit number for the daltoninput
INTEGER            :: IDUMMY,LUPRI
character(len=100)  :: WORD
character(len=2)   :: PROMPT
LOGICAL            :: DONE_LINSCA,file_exist

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,*) '-----------------------------------------'
WRITE(LUPRI,*) '         PRINTING THE LSDALTON.INP FILE '
WRITE(LUPRI,*) '-----------------------------------------'
WRITE(LUPRI,*) '                     '

INQUIRE(file='LSDALTON.INP',EXIST=file_exist) 
IF(file_exist)THEN
  LUCMD=-1
  CALL LSOPEN(LUCMD,'LSDALTON.INP','OLD','FORMATTED')
ELSE
  CALL LSQUIT('LSDALTON.INP does not exist',lupri)
ENDIF
rewind(LUCMD)
DO
   READ (LUCMD, '(A100)') WORD
   PROMPT = WORD(1:2)
   IF(PROMPT(1:1) .EQ. '!' .OR. PROMPT .EQ. '#') CYCLE
   SELECT CASE(WORD) 
   CASE ('*END OF INPUT')
      WRITE(LUPRI,'(2X,A100)') WORD   
      EXIT
   CASE DEFAULT
      WRITE(LUPRI,'(2X,A100)') WORD
   END SELECT
ENDDO
WRITE(LUPRI,*)' '

CALL LSCLOSE(LUCMD,'KEEP')

END SUBROUTINE PRINT_DALTON_FILE

!> \brief print molecule file
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE PRINT_MOL_FILE(LUPRI)
implicit none
!> the logical unit number for the output file
integer            :: lupri
!
integer            :: LUINFO,IDUMMY,I,J,IPOS,IPOS2
logical            :: file_exist,Angstrom,Symmetry,dopbc
logical            :: ATOMBASIS,BASIS,ATOMDF,AUXBASIS
CHARACTER(len=80)  :: WORD
integer            :: Atomtypes,natoms,Molecularcharge,ios
integer,ALLOCATABLE:: BasisSetCharge(:) !Charge of each type 
real(realk)        :: Q
logical            :: DoOwn,DoCartesian
CHARACTER(len=1) :: KASYM(3,3),CRT
CHARACTER(len=80),ALLOCATABLE  :: ATOMBASISSET(:)
CHARACTER(len=2) :: SYMTXT,ID3

WRITE(LUPRI,*) '                     '
WRITE(LUPRI,*) '-----------------------------------------'
WRITE(LUPRI,*) '         PRINTING THE MOLECULE.INP FILE '
WRITE(LUPRI,*) '-----------------------------------------'
WRITE(LUPRI,*) '                     '
nAtoms = 0
INQUIRE(file='MOLECULE.INP',EXIST=file_exist) 
IF(file_exist)THEN
  LUINFO=-1
  CALL LSOPEN(LUINFO,'MOLECULE.INP','OLD','FORMATTED')
ELSE
  CALL LSQUIT('MOLECULE.INP does not exist',lupri)
ENDIF

rewind(LUINFO)
DO
   READ (LUINFO, '(A80)') WORD
   IF (WORD(1:1) == '!' .or. WORD(1:1) == '#') CYCLE
   SELECT CASE(WORD)
   CASE('BASIS') 
      WRITE(LUPRI,'(2X,A40)') WORD 
      DO
         READ (LUINFO, '(A80)') WORD
         IF (WORD(1:1) == '!' .or. WORD(1:1) == '#') CYCLE
         WRITE(LUPRI,'(2X,A40)') WORD 
         EXIT
      ENDDO
   CASE DEFAULT
      WRITE(LUPRI,'(2X,A40)') WORD     
   END SELECT
   READ (LUINFO, '(A80)') WORD
   WRITE(LUPRI,'(2X,A40)') WORD     
   READ (LUINFO, '(A80)') WORD
   WRITE(LUPRI,'(2X,A40)') WORD     
   EXIT
ENDDO

READ (LUINFO, '(a80)') WORD
! Number of different atom types
IPOS = INDEX(WORD,'Ato')
IF (IPOS .NE. 0) THEN
   IPOS2 = INDEX(WORD(IPOS:80),'=')
   IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 10)) THEN
      WRITE (LUPRI,'(2X,A50)') 'Incorrect input for # atomtypes'
      WRITE (LUPRI,'(2X,A50)') 'Format is "Atomtypes=?"'
      CALL LSQUIT('Incorrect input for # atomtypes',lupri)
   ELSE
      READ (WORD((IPOS+IPOS2):80),*) Atomtypes
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
ELSE
   !*******************************************************************
   !*  OLD INDPUT FORMAT
   !******************************************************************
   READ (WORD,'(BN,A1,I4,I3,A2,10A1,D10.2,6I5)',IOSTAT=ios) CRT,&
        & Atomtypes,MolecularCharge,SYMTXT,((KASYM(I,J),I=1,3),J=1,3),&
        & ID3
   IF(IOS .NE. 0) THEN
      WRITE (LUPRI,'(A)') ' Error in the determination of the number &
           & of atomic types'
      WRITE (LUPRI,*) "Correct input structure is: Atomtypes=???"
      CALL LSQUIT('Error in determining the number of different &
           & atom types',lupri)
   ENDIF
ENDIF

   WRITE(LUPRI,'(2X,A60)') WORD
   
   DO I=1,Atomtypes
      
      READ (LUINFO, '(a80)') WORD
      IPOS = INDEX(WORD,'Ato')
      IF (IPOS .NE. 0) THEN
         IPOS2 = INDEX(WORD(IPOS:),'=')
         IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
            WRITE (LUPRI,'(2X,A50)') 'Incorrect input for # of atoms'
            WRITE (LUPRI,'(2X,A50)') 'Format is "Atoms=?"'
            CALL LSQUIT('Incorrect input for # of atoms',lupri)
         ELSE
            READ (WORD((IPOS+IPOS2):),*) nAtoms
         ENDIF
      ELSE ! OLD INPUT STYLE
         READ (WORD,'(BN,F10.0,I5)',IOSTAT=ios) Q,nAtoms
         IF(IOS .NE. 0) THEN
            WRITE (LUPRI,'(2X,A50,I4)') ' Error in the determination of the number &
                 & of atoms  for type ',I
            WRITE (LUPRI,'(2X,A50)') "Correct input structure is: Atomtypes=???"
            CALL LSQUIT('Error in determining the number of atoms',lupri)
         ENDIF
      ENDIF
      WRITE(LUPRI,'(2X,A60)')  WORD
      
!       IF(nAtoms.GT. 10)THEN
!          DO J=1,10
!             READ (LUINFO, '(a80)') WORD
!             WRITE(LUPRI,'(2X,A60)') WORD
!          ENDDO
!          WRITE(LUPRI,'(2X,A60)') 'WARNING: SINCE YOU HAVE MORE THAN 10 ATOMS&
!               & OF THIS'
!          WRITE(LUPRI,'(2X,A60)') 'TYPE I WILL NOT PRINT THEM ALL TO LIMIT OUTPUT'
!          DO J=11,nAtoms
!             READ (LUINFO, '(a80)') WORD
!          ENDDO
!       ELSE
         DO J=1,nAtoms
            READ (LUINFO, '(a80)') WORD
            WRITE(LUPRI,'(2X,A60)') WORD
         ENDDO
!      ENDIF
   ENDDO
   
   CALL LSCLOSE(LUINFO,'KEEP')

END SUBROUTINE PRINT_MOL_FILE

!> \brief build a lsitem as a fragment of a full lsitem, used in DEC
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE build_ccfragmentlsitem(LSFULL,FRAGMENT,ATOMS,NATOMS,LUPRI,IPRINT)
#ifdef VAR_MPI
use infpar_module
#endif
implicit none
!> full lsitem structure 
TYPE(LSITEM)       :: LSFULL
!> fragment lsitem structure to be build 
TYPE(LSITEM)       :: FRAGMENT
!> number of atoms in the fragment
integer            :: NATOMS
!> list of atoms in the fragment
integer            :: ATOMS(NATOMS)
!> the logical unit number for the output file
integer            :: LUPRI
!> the printlevel integer, determining how much output should be generated
integer            :: IPRINT
!
integer            :: IAO
CHARACTER(len=9)   :: BASISLABEL

FRAGMENT%lupri = lsfull%lupri
FRAGMENT%luerr = lsfull%luerr
FRAGMENT%optlevel = lsfull%optlevel
FRAGMENT%fopen = lsfull%fopen

CALL io_init(FRAGMENT%INPUT%IO)
!no pointers so a simple copy is fine
FRAGMENT%INPUT%DALTON = lsfull%INPUT%DALTON
FRAGMENT%INPUT%DALTON%DFT = lsfull%INPUT%DALTON%DFT

NULLIFY(FRAGMENT%INPUT%BASIS)
ALLOCATE(FRAGMENT%INPUT%BASIS)
CALL Copy_basissetinfo(lsfull%INPUT%BASIS%REGULAR,fragment%INPUT%BASIS%REGULAR)

!The values of daltonitem which are set in readmolfile
fragment%input%DALTON%AUXBASIS = lsfull%input%DALTON%AUXBASIS
fragment%input%DALTON%CABSBASIS = lsfull%input%DALTON%CABSBASIS
fragment%input%DALTON%JKBASIS = lsfull%input%DALTON%JKBASIS
FRAGMENT%input%dalton%DoSpherical  = lsfull%input%DALTON%DoSpherical 
FRAGMENT%input%numFragments  = lsfull%input%numFragments 
FRAGMENT%input%DALTON%NOGCINTEGRALTRANSFORM=.TRUE. 

IF(FRAGMENT%INPUT%DALTON%AUXBASIS)THEN
   CALL Copy_basissetinfo(lsfull%input%BASIS%AUXILIARY,fragment%input%BASIS%AUXILIARY)
ELSE
   call nullifyBasis(fragment%input%BASIS%AUXILIARY)
ENDIF

IF(FRAGMENT%INPUT%DALTON%CABSBASIS)THEN
   CALL Copy_basissetinfo(lsfull%input%BASIS%CABS,fragment%input%BASIS%CABS)
ELSE
   call nullifyBasis(fragment%input%BASIS%CABS)
ENDIF

IF(FRAGMENT%INPUT%DALTON%JKBASIS)THEN
   CALL Copy_basissetinfo(lsfull%input%BASIS%JK,fragment%input%BASIS%JK)
ELSE
   call nullifyBasis(fragment%input%BASIS%JK)
ENDIF

call nullifyBasis(fragment%input%BASIS%VALENCE)
call nullifyBasis(fragment%input%BASIS%GCtrans)
fragment%INPUT%BASIS%GCtransAlloc = .FALSE.

NULLIFY(FRAGMENT%INPUT%MOLECULE)
ALLOCATE(FRAGMENT%INPUT%MOLECULE)
CALL BUILD_FRAGMENT(lsfull%input%MOLECULE,FRAGMENT%input%MOLECULE,fragment%input%BASIS,&
     &FRAGMENT%INPUT%DALTON%AUXBASIS,FRAGMENT%INPUT%DALTON%CABSBASIS,&
     &FRAGMENT%INPUT%DALTON%JKBASIS,ATOMS,nATOMS,LUPRI)

IF(IPRINT .GT. 0)THEN
   CALL PRINT_MOLECULEINFO(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS,IPRINT)
   CALL PRINT_MOLECULE_AND_BASIS(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS%REGULAR)
   IF(FRAGMENT%INPUT%DALTON%AUXBASIS)THEN
      CALL PRINT_MOLECULE_AND_BASIS(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS%AUXILIARY)
   ENDIF
ENDIF

CALL typedef_init_setting(FRAGMENT%setting)
CALL typedef_set_default_setting(FRAGMENT%SETTING,FRAGMENT%INPUT)
FRAGMENT%SETTING%IntegralTransformGC =.FALSE.

! Set setting communicator to be local group communicator
#ifdef VAR_MPI
fragment%setting%comm = infpar%lg_comm
fragment%setting%node = infpar%lg_mynum
fragment%setting%numnodes = infpar%lg_nodtot
#else
fragment%setting%comm = 0
fragment%setting%node = 0
fragment%setting%numnodes = 1
#endif
END SUBROUTINE build_ccfragmentlsitem

END MODULE DALTONINFO
