!> @file
!> Contains module that reads the dalton input and initializes lsitem
MODULE DALTONINFO
use memory_handling
use precision
use LSparameters
use files, only: lsopen,lsclose
use typedeftype, only: LSITEM, integralconfig, DALTONINPUT,&
     & BASISSETLIBRARYITEM
use typedef, only: typedef_init_setting, getNbasis, PRINT_MOLECULEINFO,&
     & PRINT_MOLECULE_AND_BASIS, typedef_set_default_setting, &
     & typedef_free_setting, print_basissetlibrary
use lattice_type, only: lvec_list_t
use lattice_vectors, only: pbc_setup_default
use basis_type, only: copy_basissetinfo, free_basissetinfo,&
     & freeBrakebasinfo,initBrakebasinfo,&
     & copybrakebasinfo,buildbasisfrombrakebasinf,print_basissetinfo,&
     & print_brakebas, add_basissetinfo
use basis_typetype,only: nullifyBasisset,nullifyMainBasis,&
     & BasParamLABEL,nBasisBasParam,GCTBasParam,BRAKEBASINFO,&
     & RegBasParam,CABBasParam
use io, only: io_init, io_free
use molecule_type, only: free_moleculeinfo
use READMOLEFILE, only: read_molfile_and_build_molecule,Geometry_analysis  
use BuildBasisSet, only: Build_BASIS,copy_Molecule_IDTYPE
use lstiming, only: lstimer
use molecule_module, only: build_fragment,build_fragment2,&
     & DETERMINE_FRAGMENTNBAST,determine_nbast
use screen_mod, only: screen_init, screen_free
private
public ::  ls_init, ls_free, build_ccfragmentlsitem, dalton_finalize,&
     & build_atomspecfragmentlsitem
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
!!!!!!!!!!!!!!!!!!!!!!!!
TYPE(LSITEM) :: FRAGMENT
integer :: nATOMS,nBASIS,nbasisfull,I
integer,pointer :: ATOMS(:),BASIS(:)
ls%lupri = lupri
ls%luerr = luerr
ls%optlevel = 3 
call init_AO_parameters
CALL dalton_init(ls%input,lupri,luerr,nbast,integral,dodft,doDEC,doprint)
CALL typedef_init_setting(ls%setting)
call screen_init()
CALL typedef_set_default_setting(ls%setting,ls%input)


!!$WRITE(lupri,*)'LS%input%BASIS%BINFO(1)'
!!$call print_basissetinfo(LUPRI,LS%input%BASIS%BINFO(1))
!!$
!!$nATOMS = 7
!!$nbasis = 19
!!$nbasisfull = nbast
!!$call mem_alloc(ATOMS,nATOMS)
!!$do I=1,nATOMS
!!$   ATOMS(I) = I
!!$enddo
!!$call mem_alloc(BASIS,nBASIS)
!!$do I=1,nBASIS
!!$   BASIS(I) = I
!!$enddo
!!$call build_ccfragmentlsitem(LS,FRAGMENT,ATOMS,NATOMS,BASIS,NBASIS,NBASISFULL,LUPRI,1000)!IPRINT)
!!$call mem_dealloc(ATOMS)
!!$call mem_dealloc(BASIS)
!!$
!!$WRITE(lupri,*)'THE MODIFIED CC FRAGMENT BASIS'
!!$IF (doprint) THEN
!!$   WRITE(lupri,*)'FRAGMENT%input%BASIS%BINFO(1)'
!!$   call print_basissetinfo(LUPRI,FRAGMENT%input%BASIS%BINFO(1))
!!$
!!$   CALL PRINT_MOLECULEINFO(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS,1000)
!!$   do i=1,nBasisBasParam
!!$      IF(ls%input%BASIS%WBASIS(I))THEN
!!$         CALL PRINT_MOLECULE_AND_BASIS(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS%BINFO(I))
!!$      ENDIF
!!$   enddo
!!$ENDIF
!!$
!!$call lsquit('test done',-1)
!!$

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
TYPE(BASISSETLIBRARYITEM) :: LIBRARY(nBasisBasParam)
integer              :: LUCME,IDUMMY,I
CHARACTER(len=80)    :: BASISSETNAME
CHARACTER(len=9)     :: BASISLABEL
real(realk)          :: Tim1,Tim2,TIME_BUILD_BASIS,TIME_READ_DALTON
real(realk)          :: TIME_AUXBUILD_BASIS,TIME_READ_MOL,TSTART
Integer              :: numAtoms, numNodes
integer,pointer :: atoms(:)
logical              :: lopen
LUCME=-1
!CAREFUL

intinp%nfock = 0 !number of fock matrices calculated

NULLIFY(intinp%MOLECULE)
NULLIFY(intinp%AUXMOLECULE)
NULLIFY(intinp%BASIS)
ALLOCATE(intinp%MOLECULE)
ALLOCATE(intinp%AUXMOLECULE)
intinp%AUXMOLECULE%nAtoms=0
intinp%AUXMOLECULE%nSubSystems=0
ALLOCATE(intinp%BASIS)
call nullifyMainBasis(intinp%BASIS)

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
  call pbc_setup_default(latt_config)
  CALL READ_MOLFILE_AND_BUILD_MOLECULE(LUPRI,intinp%MOLECULE,LIBRARY,doprint,&
       & intinp%dalton%molprint,intinp%dalton%DoSpherical,intinp%dalton%basis,&
       & latt_config,intinp%dalton%atombasis)

CALL Geometry_analysis(intinp%MOLECULE,LUPRI)  

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

IF (doprint) THEN
   do i=1,nBasisBasParam
      IF(intinp%DALTON%BASIS(i))THEN
         WRITE(lupri,*)'BASISSETLIBRARY  : ',BasParamLABEL(i)
         CALL PRINT_BASISSETLIBRARY(LUPRI,LIBRARY(i))
      ENDIF
   enddo
ENDIF

IF(intinp%DALTON%TIMINGS) CALL LSTIMER('BUILD BASIS',TIM1,TIM2,LUPRI)

call nullifyMainBasis(intinp%BASIS)
intinp%BASIS%WBASIS = intinp%DALTON%BASIS
do i=1,nBasisBasParam
   IF(intinp%BASIS%WBASIS(I))THEN
      IF (doprint) THEN
         WRITE(LUPRI,'(2X,A)')' '
         WRITE(LUPRI,'(2X,3A)')'CALLING BUILD BASIS WITH DEFAULT ',BasParamLABEL(i),' BASIS'
         WRITE(LUPRI,'(2X,A)')' '
      ENDIF
      CALL Build_basis(LUPRI,intinp%DALTON%BASPRINT,&
           & intinp%MOLECULE,intinp%BASIS%BINFO(I),LIBRARY,&
           & BasParamLABEL(I),intinp%DALTON%UNCONT,intinp%DALTON%NOSEGMENT,&
           & doprint,intinp%DALTON%DOSPHERICAL,I)
      IF(i.EQ.RegBasParam)THEN !regular basis
         nbast = getNbasis(AORdefault,Contractedinttype,intinp%MOLECULE,LUPRI)
      ENDIF
      IF(intinp%DALTON%TIMINGS) CALL LSTIMER('BUILD '//BasParamLABEL(I),TIM1,TIM2,LUPRI)
   ENDIF
ENDDO

IF (doprint) THEN
   WRITE(lupri,*)'PRINT MOLECULE AND BASIS IN FULL INPUT'
   CALL PRINT_MOLECULEINFO(LUPRI,intinp%MOLECULE,intinp%BASIS,intinp%DALTON%MOLPRINT)
   do i=1,nBasisBasParam
      IF(intinp%BASIS%WBASIS(I))THEN
         CALL PRINT_MOLECULE_AND_BASIS(LUPRI,intinp%MOLECULE,intinp%BASIS%BINFO(I))
      ENDIF
   enddo
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
integer :: ibas

call free_Moleculeinfo(intinp%MOLECULE)
call free_Moleculeinfo(intinp%AUXMOLECULE)
DO ibas=1,nBasisBasParam
   IF(intinp%BASIS%WBASIS(ibas))THEN
      IF(intinp%BASIS%BINFO(ibas)%natomtypes .NE. 0)THEN
         call free_basissetinfo(intinp%BASIS%BINFO(ibas))
      ELSE
         call lsquit('Error in dalton_finalize',-1)
      ENDIF
   ENDIF
ENDDO
DEALLOCATE(intinp%MOLECULE)
DEALLOCATE(intinp%AUXMOLECULE)
DEALLOCATE(intinp%BASIS)
call io_free(intinp%io)

END SUBROUTINE dalton_finalize

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
logical            :: BASIS,ATOMDF,AUXBASIS
CHARACTER(len=80)  :: WORD
integer            :: Atomtypes,natoms,Molecularcharge,ios
integer,ALLOCATABLE:: BasisSetCharge(:) !Charge of each type 
real(realk)        :: Q
logical            :: DoOwn,DoCartesian
CHARACTER(len=1) :: KASYM(3,3),CRT
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
SUBROUTINE build_AtomSpecfragmentlsitem(LSFULL,FRAGMENT,ATOMS,NATOMS,LUPRI,IPRINT)
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
integer            :: IAO,I
CHARACTER(len=9)   :: BASISLABEL
integer,pointer :: FullAtomList(:)

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
call nullifyMainBasis(FRAGMENT%INPUT%BASIS)
FRAGMENT%INPUT%BASIS%WBASIS = lsfull%INPUT%BASIS%WBASIS
do i=1,nBasisBasParam
 IF(lsfull%INPUT%BASIS%WBASIS(I))THEN
  CALL Copy_basissetinfo(lsfull%INPUT%BASIS%BINFO(i),fragment%INPUT%BASIS%BINFO(I))
 ENDIF
enddo

IF(.NOT.lsfull%input%DALTON%NOGCINTEGRALTRANSFORM)THEN
   CALL LSQUIT('NOGCINTEGRALTRANSFORM false in build_AtomSpecfragmentlsitem',-1)
ENDIF

!The values of daltonitem which are set in readmolfile
FRAGMENT%input%dalton%DoSpherical  = lsfull%input%DALTON%DoSpherical 
FRAGMENT%input%numFragments  = lsfull%input%numFragments 

NULLIFY(FRAGMENT%INPUT%MOLECULE)
ALLOCATE(FRAGMENT%INPUT%MOLECULE)
CALL BUILD_FRAGMENT(lsfull%input%MOLECULE,FRAGMENT%input%MOLECULE,&
     & fragment%input%BASIS,ATOMS,nATOMS,LUPRI)
NULLIFY(FRAGMENT%INPUT%AUXMOLECULE)
ALLOCATE(FRAGMENT%INPUT%AUXMOLECULE)
call mem_alloc(FullAtomList,lsfull%input%MOLECULE%natoms)
do I=1,lsfull%input%MOLECULE%natoms
   FullAtomList(I) = I
enddo
CALL BUILD_FRAGMENT(lsfull%input%MOLECULE,FRAGMENT%input%AUXMOLECULE,&
     & fragment%input%BASIS,FullAtomList,lsfull%input%MOLECULE%natoms,LUPRI)
call mem_dealloc(FullAtomList)

IF(IPRINT .GT. 0)THEN
   CALL PRINT_MOLECULEINFO(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS,IPRINT)
   do i=1,nBasisBasParam
    IF(lsfull%input%BASIS%WBASIS(I))THEN
     CALL PRINT_MOLECULE_AND_BASIS(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS%BINFO(I))
    ENDIF
   enddo
ENDIF

CALL typedef_init_setting(FRAGMENT%setting)
CALL typedef_set_default_setting(FRAGMENT%SETTING,FRAGMENT%INPUT)
IF(FRAGMENT%SETTING%IntegralTransformGC)THEN
   CALL LSQUIT('IntegralTransformGC true in build_AtomSpecfragmentlsitem',-1)
ENDIF
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
END SUBROUTINE build_AtomSpecfragmentlsitem

!> \brief build a lsitem as a fragment of a full lsitem, used in DEC
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE build_ccfragmentlsitem(LSFULL,FRAGMENT,ATOMS,NATOMS,BASIS,NBASIS,NBASISFULL,LUPRI,IPRINT)
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
!> number of aos (basis functions) in the fragment
integer            :: NBASIS
!> number of aos (basis functions) in the molecule
integer            :: NBASISFULL
!> list of basis functions  in the fragment
integer            :: BASIS(NBASIS)
!> the logical unit number for the output file
integer            :: LUPRI
!> the printlevel integer, determining how much output should be generated
integer            :: IPRINT
!
integer            :: IAO,I,NATOMSFULL
CHARACTER(len=9)   :: BASISLABEL
logical            :: WhichAos(NBASISFULL)

NATOMSFULL = lsfull%input%MOLECULE%nAtoms
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
call nullifyMainBasis(FRAGMENT%INPUT%BASIS)
FRAGMENT%INPUT%BASIS%WBASIS = lsfull%INPUT%BASIS%WBASIS
Call BuildNewBasInfoStructure(LSFULL,FRAGMENT,ATOMS,NATOMS,NATOMSFULL,BASIS,NBASIS,NBASISFULL,LUPRI,IPRINT)

IF(.NOT.lsfull%input%DALTON%NOGCINTEGRALTRANSFORM)THEN
   CALL LSQUIT('NOGCINTEGRALTRANSFORM false in build_ccfragmentlsitem',-1)
ENDIF

!The values of daltonitem which are set in readmolfile
FRAGMENT%input%dalton%DoSpherical  = lsfull%input%DALTON%DoSpherical 
FRAGMENT%input%numFragments  = lsfull%input%numFragments 

IF(IPRINT .GT. 0)THEN
   WRITE(lupri,*)'FRAGMENT%input%BASIS%BINFO(1)'
   call print_basissetinfo(LUPRI,FRAGMENT%input%BASIS%BINFO(1))
   WRITE(lupri,*)'PRINT MOLECULE AND BASIS IN CCFRAGMENT'
   CALL PRINT_MOLECULEINFO(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS,IPRINT)
   do i=1,nBasisBasParam
    IF(lsfull%input%BASIS%WBASIS(I))THEN
     CALL PRINT_MOLECULE_AND_BASIS(LUPRI,FRAGMENT%input%MOLECULE,FRAGMENT%input%BASIS%BINFO(I))
    ENDIF
   enddo
ENDIF

CALL typedef_init_setting(FRAGMENT%setting)
CALL typedef_set_default_setting(FRAGMENT%SETTING,FRAGMENT%INPUT)
IF(FRAGMENT%SETTING%IntegralTransformGC)THEN
   CALL LSQUIT('IntegralTransformGC true in build_ccfragmentlsitem',-1)
ENDIF
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

SUBROUTINE BuildNewBasInfoStructure(LSFULL,FRAGMENT,ATOMS,NATOMS,NATOMSFULL,BASIS,NBASIS,NBASISFULL,LUPRI,IPRINT)
implicit none
!> full lsitem structure 
TYPE(LSITEM)       :: LSFULL
!> fragment lsitem structure to be build 
TYPE(LSITEM)       :: FRAGMENT
!> number of atoms in the fragment
integer            :: NATOMS
!> list of atoms in the fragment
integer            :: ATOMS(NATOMS)
!> number of atoms in the full
integer            :: NATOMSFULL
!> number of aos (basis functions) in the fragment
integer            :: NBASIS
!> list of basis functions  in the fragment
integer            :: BASIS(NBASIS)
!> number of aos (basis functions) in the full molecule
integer            :: NBASISFULL
!> the logical unit number for the output file
integer            :: LUPRI
!> the printlevel integer, determining how much output should be generated
integer            :: IPRINT
!
TYPE(BRAKEBASINFO) :: BASINFO
TYPE(BRAKEBASINFO) :: BASINFOARRAY(NATOMS)
logical :: WhichAos(NBASISFULL),ContainOrb,WhichAtoms(nAtomsFull)
integer :: iB,r,iOrbitalIndex,i,icharge,itype,MaxnOrb,nangmom,kmult,iAngNew,ang
integer :: nOrb,iOrbNew,nBASINFOARRAY,itypeOld,iOrb,iK,unique,iAtom,ik1,ik2
integer :: iOrbitalIndexSS
integer,pointer :: newType(:),FullAtomList(:)
call mem_alloc(newType,natoms)
WhichAos = .FALSE.
DO I=1,NBASIS
   WhichAos(BASIS(I)) = .TRUE.
ENDDO
WhichAtoms = .FALSE.
DO I=1,NATOMS
   WhichAtoms(ATOMS(I)) = .TRUE.
ENDDO

nBASINFOARRAY = 0
iB=RegBasParam
r = lsfull%input%basis%binfo(IB)%labelindex
iOrbitalIndex = 0
itypeOld = -24
iAtom = 0 
DO i=1,natomsfull
 if(r == 0) then            
    icharge = int(lsfull%input%molecule%atom(i)%charge)
    itype = lsfull%input%basis%binfo(iB)%chargeindex(icharge)
 else
    itype = lsfull%input%molecule%atom(i)%idtype(r)
 end if

 IF(WhichAtoms(i))THEN
   iAtom = iAtom + 1 
   nangmom = lsfull%input%basis%binfo(iB)%ATOMTYPE(itype)%nAngmom
   MaxnOrb = 0 
   do ang = 1,nAngmom
      MaxnOrb = MAX(MaxnOrb,lsfull%input%basis%binfo(iB)%ATOMTYPE(itype)%SHELL(ang)%norb)
   enddo
!   IF(itype.NE.itypeOld)THEN
      IF(itypeOld.NE.-24)call freeBrakebasinfo(BASINFO)
      call initBrakebasinfo(BASINFO,itype,nangmom,MaxnOrb)
      itypeOld = itype
!   ENDIF
   kmult = 1
   iAngNew = 1
   do ang = 0,nAngmom-1
      nOrb = lsfull%input%basis%binfo(iB)%ATOMTYPE(itype)%SHELL(ang+1)%norb
      BASINFO%OriginalnOrb(ang+1) = nOrb      
      iOrbNew = 1 
      do iOrb = 1,nOrb
         ContainOrb = .TRUE.
         iOrbitalIndexSS = iOrbitalIndex
         do iK = 1,kmult
            iOrbitalIndex = iOrbitalIndex + 1
            IF(.NOT.WhichAos(iOrbitalIndex))THEN !at least one orb comp is false
               !ensure WhichAos are false for all Orbital components
               do iK1 = 1,kmult
                  IF(WhichAos(iOrbitalIndexSS+iK1))THEN !At least one orb comp is true                      
                     CALL LSQUIT('Only part of Orbital componets include ',-1)
                  ENDIF
               enddo
               ContainOrb = .FALSE.
            ENDIF
         enddo
         IF(ContainOrb)THEN
            BASINFO%Orb(iOrbNew,iAngNew) = iOrb
            iOrbNew = iOrbNew +  1
         ENDIF
      enddo
      BASINFO%nOrb(iAngNew) = iOrbNew-1
      IF(BASINFO%nOrb(iAngNew).NE.0)THEN
         BASINFO%angmom(iAngNew) = ang 
         iAngNew = iAngNew + 1 
      ENDIF
      kmult = kmult + 2
   enddo
   BASINFO%nangmom = iAngNew-1 
   IF(iAngNew.EQ.1)THEN
      !atom not included 
      call lsquit('atom not included in BASIS list, but included in ATOMS list',-1)
   ELSE
      !Determine if unique BASINFO structure
      call DetermineIfUniqueBASINFO(BASINFO,BASINFOARRAY,NATOMS,nBASINFOARRAY,unique)
      IF(unique.EQ.0)THEN
!         call print_Brakebas(BASINFO,6)
         call copyBrakebasinfo(BASINFO,BASINFOARRAY(nBASINFOARRAY+1))
         nBASINFOARRAY = nBASINFOARRAY + 1 
         newType(iAtom) = nBASINFOARRAY
      ELSE
         newType(iAtom) = unique
      ENDIF         
   ENDIF
 ELSE
    iOrbitalIndex = iOrbitalIndex + lsfull%input%basis%binfo(iB)%ATOMTYPE(itype)%TotnOrb
 ENDIF
ENDDO
IF(itypeOld.NE.-24)call freeBrakebasinfo(BASINFO)

Call BuildBasisFromBRAKEBASINF(BASINFOARRAY,NATOMS,nBASINFOARRAY,&
     & lsfull%INPUT%BASIS%BINFO(iB),fragment%INPUT%BASIS%BINFO(IB),RegBasParam)

DO I=1,nBASINFOARRAY
   call freeBrakebasinfo(BASINFOARRAY(I))
ENDDO

!copy the rest (for now the full Aux/CABS AO set is used for each atom)
do iB=2,nBasisBasParam
   IF(lsfull%INPUT%BASIS%WBASIS(IB))THEN
      CALL Copy_basissetinfo(lsfull%INPUT%BASIS%BINFO(iB),fragment%INPUT%BASIS%BINFO(IB))
   ENDIF
enddo
NULLIFY(FRAGMENT%INPUT%MOLECULE)
ALLOCATE(FRAGMENT%INPUT%MOLECULE)
CALL BUILD_FRAGMENT2(lsfull%input%MOLECULE,FRAGMENT%input%MOLECULE,&
     & fragment%input%BASIS,ATOMS,nATOMS,LUPRI)

call mem_alloc(FullAtomList,natomsfull)
do I=1,natomsfull
   FullAtomList(I) = I
enddo
NULLIFY(FRAGMENT%INPUT%AUXMOLECULE)
ALLOCATE(FRAGMENT%INPUT%AUXMOLECULE)
CALL BUILD_FRAGMENT2(lsfull%input%MOLECULE,FRAGMENT%input%AUXMOLECULE,&
     & lsfull%input%BASIS,FullAtomList,natomsfull,LUPRI)
call mem_dealloc(FullAtomList)

DO i=1,natoms
   fragment%input%molecule%atom(i)%idtype(RegBasParam) = newType(i)
ENDDO
call mem_dealloc(newType)

call DETERMINE_FRAGMENTNBAST(lsfull%input%MOLECULE,FRAGMENT%input%AUXMOLECULE,&
     & lsfull%input%BASIS,LUPRI)
call DETERMINE_FRAGMENTNBAST(lsfull%input%MOLECULE,FRAGMENT%input%MOLECULE,&
     & fragment%input%BASIS,LUPRI)

do iB=nBasisBasParam,1,-1
   IF(lsfull%INPUT%BASIS%WBASIS(IB))THEN
      CALL DETERMINE_NBAST(FRAGMENT%INPUT%AUXMOLECULE,lsfull%INPUT%BASIS%BINFO(iB),&
           & FRAGMENT%input%dalton%DoSpherical,.FALSE.)
   ENDIF
enddo
do iB=nBasisBasParam,1,-1
   IF(lsfull%INPUT%BASIS%WBASIS(IB))THEN
      CALL DETERMINE_NBAST(FRAGMENT%INPUT%MOLECULE,FRAGMENT%INPUT%BASIS%BINFO(iB),&
           & FRAGMENT%input%dalton%DoSpherical,.FALSE.)
   ENDIF
enddo

end SUBROUTINE BuildNewBasInfoStructure

subroutine DetermineIfUniqueBASINFO(BASINFO,BASINFOARRAY,NATOMS,nBASINFOARRAY,unique)
implicit none
integer,intent(inout) :: unique
integer,intent(in) :: NATOMS,nBASINFOARRAY
TYPE(BRAKEBASINFO),intent(in) :: BASINFO
TYPE(BRAKEBASINFO) :: BASINFOARRAY(NATOMS)
!
integer :: i
!logical,external :: IdenticalBrakebasinfo

unique = 0
IF(nBASINFOARRAY.EQ.0)THEN
   unique = 0
ELSE
   basinfoloop: DO i = 1,nBASINFOARRAY
      IF(IdenticalBrakebasinfo(BASINFO,BASINFOARRAY(i)))THEN
         unique = i
         exit basinfoloop         
      ENDIF
   ENDDO basinfoloop
ENDIF
CONTAINS
logical function IdenticalBrakebasinfo(BBASINFO,BBASINFO2)
implicit none
TYPE(BRAKEBASINFO) :: BBASINFO,BBASINFO2
!
integer :: I,J
IdenticalBrakebasinfo = .TRUE.
IF(BBASINFO%OriginalType.NE.BBASINFO2%OriginalType)THEN
   IdenticalBrakebasinfo = .FALSE. !not same original type
ELSE
   IF(BBASINFO%nAngmom.NE.BBASINFO2%nAngmom)THEN
      IdenticalBrakebasinfo = .FALSE. !not same number of angular moments
   ELSE
      DO I=1,BBASINFO%nAngmom
         IF(BBASINFO%nOrb(I).NE.BBASINFO2%nOrb(I).OR.&
              & BBASINFO%Angmom(I).NE.BBASINFO2%Angmom(I))THEN
            IdenticalBrakebasinfo = .FALSE. !not same number of orbitals or angmom
         ELSE
            DO J=1,BBASINFO%nOrb(I)
               IF(BBASINFO%Orb(J,I).NE.BBASINFO2%Orb(J,I))THEN
                  IdenticalBrakebasinfo = .FALSE. !not same orbitals
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDIF
ENDIF
end function IdenticalBrakebasinfo

end subroutine DetermineIfUniqueBASINFO



END MODULE DALTONINFO
