!> @file
!> Module containing subroutines for building basisSetInfo and basisSetLibrary
MODULE BUILDBASISSET
  use precision
  use TYPEDEF
  use LSMatrix_Operations_dense
  use lstiming
  use memory_handling
  use AO_TypeType, only: ExpThr
  use molecule_type
  use files
  use molecule_module
  TYPE CAR200
    CHARACTER(len=200) :: p
  END TYPE CAR200
contains
!> \brief builds the basis structure
!> \author T. Kjaergaard
!> \date 2008
!>
!> THIS IS THE MAIN DRIVER WHICH BUILDS THE BASISSETINFO OBJECT 
!> (See TYPE-DEF.f90) IF The BASISSETNAME IS GIVEN IT WILL BUILD 
!> THE BASISINFO IS BUILD FROM THIS. OTHERWISE THE BASISINFO IS 
!> BUILD FROM WHAT IS SPECIFIED FROM IN THE MOLECULEINFO OBJECT 
!> WHICH MEANS THAT IT IS BUILD FROM INFORMATION GIVEN IN 
!> THE MOLECULE INPUT FILE
!> 
!> THE DRIVER IS DEVIDED INTO 3 SPECIAL CASES.
!> 1. BASISSETNAME IS GIVEN AND WE BUILD BASISSETINFO FROM THIS  
!> 2. 'BASIS' KEYWORD IS GIVEN IN MOLECULE.INP  
!> 3. 'ATOMBASIS' KEYWORD IS GIVEN IN MOLECULE.INP  
!> THE RECIPE IS
!> 1 DETERMINE UNIQUE BASISSETS IN CASE OF MULTIPLE BASISSETS
!> 2 DETERMINE UNIQUE CHARGES FOR EACH BASISSET
!> 3 DETERMINE AUGMENTATION OF THE BASISSET
!> 4 DETERMINE POLARIZATION FUNCTIONS IF ANY SPECIFIED
!> 5 OPEN BASISSET FILE
!> 6 CALLS DRIVER WHICH READS THE BASISSET FILE AND BUILDS
!>   THE BASISSETINFO 
!> 7 CLOSING THE BASISSETFILE  
!> 8 SYNC THE MOLECULE OBJECT WITH THE BASISSET OBJECT
!>
!> THE BASISSETINFO CAN BE PRINTED WITH 'PRINT_BASISSETINFO'
!>
SUBROUTINE Build_BASIS(LUPRI,IPRINT,MOLECULE,BASINFO,BASISSETLIBRARY,&
     & BASISLABEL,UNCONTRACTED,SINGLESEGMENT,doprint,DOSPHERICAL,iBAS,&
     & BASISSETNAME)
implicit none
!> the logical unit number for the output file
INTEGER            :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER            :: IPRINT
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO) :: MOLECULE
!> the basisinfo to be build
TYPE(BASISSETINFO) :: BASINFO
!> the basissetlibrary contains info on basisset from molecule file.
TYPE(BASISSETLIBRARYITEM) :: BASISSETLIBRARY(nBasisBasParam)
!> label 'REGULAR  ' or 'AUXILIARY' or 'CABS     ' or 'JKAUX    '
CHARACTER(len=9)   :: BASISLABEL
!> if the AOitem should be uncontracted
LOGICAL            :: UNCONTRACTED
!> no segmentation of the basis is used if this is true
LOGICAL            :: SINGLESEGMENT
!> same as iprint?
LOGICAL            :: doprint
!> if you want to use spherical basis function (default) or cartesian
LOGICAL            :: dospherical
! which basis in BASISSETLIBRARY(nBasisBasParam)
INTEGER :: iBAS
!> if you want to disregard the input basis and use own basis use this option
CHARACTER(len=80),OPTIONAL  :: BASISSETNAME
!
!INTEGER,pointer         :: CHARGES(:)
REAL(REALK),pointer         :: CHARGES(:)
LOGICAL,pointer    :: uCHARGES(:)
TYPE(CAR200),pointer :: BASISDIR(:)
INTEGER,pointer    :: LEN_BASISDIR(:)
INTEGER            :: LUBAS,NEND,i,j,k,IPOLST,IAUG,NBASDIR,IPOS,iBAS2
LOGICAL            :: FILE_EXIST,POLFUN,CONTRACTED,GCONT,pointcharge,BASIS_EXIST
TYPE(CAR200),pointer :: STRING(:)
!INTEGER,pointer  :: BINDEXES(:)
real(realk)        :: tstart,tend
integer            :: IT,II,natomtypes,atomtype,MaxCharge,iatom,icharge
logical,pointer    :: POINTCHARGES(:)

CONTRACTED=.NOT.UNCONTRACTED

CALL LSTIMER('START ',Tstart,Tend,LUPRI)

CALL GET_BASISSET_LIBS(LUPRI,BASISDIR,LEN_BASISDIR,NBASDIR,IPRINT.GT.0)
ALLOCATE(STRING(NBASDIR))

BASINFO%GCbasis = .FALSE.
BASINFO%SPHERICAL = DOSPHERICAL
IF(present(BASISSETNAME))THEN
 BASINFO%Labelindex = 0
! call mem_alloc(BINDEXES,1)
! BINDEXES(1)=1
 J=1
 BASINFO%DunningsBasis = .FALSE.
 IF(IPRINT .GT. 5)WRITE(LUPRI,*)'UNIQUE BASISSETS',J
!ELSEIF(BASISSETLIBRARY(iBas)%nbasissets == 1)THEN
! !This is the case of BASIS OR ALL BASISSETS USING ATOMBASIS IS THE SAME
! BASINFO%Labelindex = 0
!! call mem_alloc(BINDEXES,1)
!! BINDEXES(1)=1
! J=1
! BASINFO%DunningsBasis = BASISSETLIBRARY(iBas)%DunningsBasis
! IF(IPRINT .GT. 5)WRITE(LUPRI,*)'UNIQUE BASISSETS',J
ELSE
! call mem_alloc(BINDEXES,BASISSETLIBRARY(iBas)%nbasissets)
 BASINFO%DunningsBasis = BASISSETLIBRARY(iBas)%DunningsBasis
 SELECT CASE(BASISLABEL)
 CASE('REGULAR  ')
    BASINFO%Labelindex = RegBasParam
 CASE('AUXILIARY')
    BASINFO%Labelindex = AUXBasParam
 CASE('CABS     ')
    BASINFO%Labelindex = CABBasParam
 CASE('JKAUX    ')
    BASINFO%Labelindex = JKBasParam
 CASE('VALENCE  ')
    BASINFO%Labelindex = VALBasParam
    call LSQUIT('VALENCE no legal keyword in Build_Basis.',lupri)  
 CASE DEFAULT
    WRITE (LUPRI,'(A5,2X,A9)') 'LABEL',BASISLABEL
    WRITE (LUPRI,'(a80)') ' not recognized in Build_basis.'
    WRITE (LUPRI,'(a80)') ' maybe it is not implemented yet.'
    CALL LSQUIT('Illegal keyword in Build_Basis.',lupri)  
 END SELECT
 J = BASISSETLIBRARY(iBas)%nbasissets
 IF(IPRINT .GT. 5)WRITE(LUPRI,*)'UNIQUE BASISSETS',J
 IF(IPRINT .GT. 5)THEN
  DO I=1,J  
   WRITE(LUPRI,'(A18,2X,A20)')'BUILD BASIS FOR ',BASISSETLIBRARY(iBas)%BASISSETNAME(I)
  ENDDO
 ENDIF
ENDIF
BASINFO%LABEL=BASISLABEL 

natomtypes = 0
DO I=1,J  
 IF(present(BASISSETNAME))THEN
!  CALL UNIQUE_CHARGES2(LUPRI,MOLECULE,CHARGES,k)
    MaxCharge = 0
    DO IATOM=1,MOLECULE%natoms  
       ICHARGE = INT(MOLECULE%ATOM(IATOM)%CHARGE)
       MaxCharge = MAX(Icharge,MaxCharge)
    ENDDO
    call mem_alloc(uCHARGES,MaxCharge)
    uCHARGES = .FALSE.
    DO IATOM=1,MOLECULE%natoms  
       ICHARGE = INT(MOLECULE%ATOM(IATOM)%CHARGE)
       uCHARGES(ICHARGE) = .TRUE.
    ENDDO
    k = 0
    DO ICharge=1,MaxCharge
       IF(uCHARGES(ICharge)) k = k +1
    ENDDO
    call mem_alloc(CHARGES,k)
    k = 0
    DO ICharge=1,MaxCharge
       IF(uCHARGES(ICharge))THEN
          k = k + 1
          CHARGES(k) = ICharge
       ENDIF
    ENDDO
    call mem_dealloc(uCHARGES)    
 ELSE
   k=BASISSETLIBRARY(iBas)%nCharges(I)
 ENDIF
 natomtypes = natomtypes + k
ENDDO
IF(natomtypes.EQ. 0)CALL LSQUIT('Error trying to build basis but found no atoms',lupri)
CALL INIT_BASISSETINFO_TYPES(BASINFO,natomtypes)

atomtype = 0
DO I=1,J  
 IF(present(BASISSETNAME))THEN
!  call mem_alloc(CHARGES,MOLECULE%nAtoms)
!  CALL UNIQUE_CHARGES(LUPRI,BASISSETLIBRARY(iBas),CHARGES,k)
!  CALL UNIQUE_CHARGES(LUPRI,BASISSETLIBRARY(iBas),CHARGES,k)
  DO II=1,k
   BASINFO%ATOMTYPE(atomtype+II)%NAME = BASISSETNAME
  ENDDO
  IF(IPRINT .GT. 5)THEN 
   WRITE(LUPRI,*)'UNIQUE CHARGES FOR BASISSET',BASISSETNAME
   DO IT=1,K
    WRITE(LUPRI,*)'CHARGE(IT):',CHARGES(IT)
   ENDDO
  ENDIF
 ELSE
  k=BASISSETLIBRARY(iBas)%nCharges(I)
  DO II=1,k
   BASINFO%ATOMTYPE(atomtype+II)%NAME = BASISSETLIBRARY(iBas)%BASISSETNAME(I)
  ENDDO
  IF(IPRINT .GT. 5)THEN 
   WRITE(LUPRI,*)'UNIQUE CHARGES FOR BASISSET',BASINFO%ATOMTYPE(atomtype+k)%NAME
   DO IT=1,K
    WRITE(LUPRI,*)'CHARGE(IT):',BASISSETLIBRARY(iBas)%CHARGES(I,IT)
   ENDDO
  ENDIF
 ENDIF
 IF(present(BASISSETNAME))THEN
  CALL DETERMINE_AUGMENTATION(LUPRI,BASISSETNAME,IAUG)
  CALL DETERMINE_POLARIZATION(LUPRI,BASISSETNAME,POLFUN,IPOLST)
  NEND = INDEX(BASISSETNAME,' ') - 1
  DO IBAS2=1,NBASDIR
    STRING(IBAS2)%p = BASISDIR(IBAS2)%p(1:LEN_BASISDIR(IBAS2))//BASISSETNAME(1:NEND)//' '
  ENDDO
  IF(BASISSETNAME(1:11).EQ.'pointcharge')THEN
     pointcharge=.TRUE.
  ELSE
     pointcharge=.FALSE.
  ENDIF
 ELSE
  CALL DETERMINE_AUGMENTATION(LUPRI,BASISSETLIBRARY(iBas)%BASISSETNAME(I),IAUG)
  CALL DETERMINE_POLARIZATION(LUPRI,BASISSETLIBRARY(iBas)%BASISSETNAME(I),POLFUN,IPOLST)
  NEND = INDEX(BASISSETLIBRARY(iBas)%BASISSETNAME(I),' ') - 1
  DO IBAS2=1,NBASDIR
    STRING(IBAS2)%p = BASISDIR(IBAS2)%p(1:LEN_BASISDIR(IBAS2))//BASISSETLIBRARY(iBas)%BASISSETNAME(I)(1:NEND)//' '
  ENDDO
  IF(BASISSETLIBRARY(iBas)%BASISSETNAME(I)(1:11).EQ.'pointcharge')THEN
     pointcharge=.TRUE.
  ELSE
     pointcharge=.FALSE.
  ENDIF
 ENDIF
 IF(pointcharge)THEN
    IF(present(BASISSETNAME))THEN
       DO II=1,k
          BASINFO%ATOMTYPE(atomtype+II)%Charge = NINT(CHARGES(II))
          BASINFO%ATOMTYPE(atomtype+II)%nAngmom=0
       ENDDO
    ELSE
       DO II=1,k
          BASINFO%ATOMTYPE(atomtype+II)%Charge = NINT(BASISSETLIBRARY(iBas)%CHARGES(I,II))
          BASINFO%ATOMTYPE(atomtype+II)%nAngmom=0
       ENDDO
    ENDIF
 ELSE
  NEND= INDEX(STRING(1)%p(1:),' ') - 1
  LUBAS=-1
  IF(STRING(1)%p(NEND-10:NEND).EQ.'USERDEFINED')THEN
     IBAS2 = 1
     CALL LSOPEN(LUBAS,'MOLECULE.INP','OLD','FORMATTED')
     call FIND_BASIS_IN_MOLFILE(LUPRI,LUBAS) 
  ELSE
     BASIS_EXIST = .FALSE.
     DO II=1,NBASDIR
       NEND= INDEX(STRING(II)%p(1:),' ') - 1
       INQUIRE (FILE = STRING(II)%p(1:NEND), EXIST = FILE_EXIST)
       IF (FILE_EXIST) THEN
         BASIS_EXIST = .TRUE.
         IBAS2 = II
         EXIT
       ENDIF
     ENDDO
     IF(BASIS_EXIST) THEN
        IF(doprint) WRITE(LUPRI,*)'OPENING FILE',STRING(IBAS2)%p(1:NEND)
        CALL LSOPEN(LUBAS,STRING(IBAS2)%p(1:NEND),'OLD','FORMATTED')
     ELSE
        IPOS = INDEX(STRING(1)%p,'/',.TRUE.)+1
        WRITE (LUPRI,'(/,3A)')'Basis "',STRING(1)%p(IPOS:NEND),&
     &    '" doesn''t exist in any of the provided basis-set libraries:'
        DO IBAS2=1,NBASDIR
          WRITE(LUPRI,'(3X,A)') BASISDIR(IBAS2)%p
        ENDDO
        CALL LSQUIT('Non-existing basis set in Linsca Integral',lupri)
     ENDIF
  ENDIF
  call mem_alloc(POINTCHARGES,k)
  IF(present(BASISSETNAME))THEN
     DO IT=1,k
        POINTCHARGES(IT) = .FALSE.
     ENDDO
  ELSE
     DO IT=1,k
        POINTCHARGES(IT) = BASISSETLIBRARY(iBas)%POINTCHARGES(I,IT)
     ENDDO     
  ENDIF
  IF(present(BASISSETNAME))THEN
   CALL READ_BASISSET_FILE_AND_BUILD_BASISINFO(LUPRI,IPRINT,LUBAS,&
        &CONTRACTED,STRING(IBAS2)%p,IAUG,POLFUN,IPOLST,BASINFO,CHARGES,&
        &k,atomtype,SINGLESEGMENT,DOSPHERICAL,POINTCHARGES)
   call mem_dealloc(CHARGES)
  ELSE
   call mem_alloc(CHARGES,k)
   DO IT=1,k
    CHARGES(IT) = BASISSETLIBRARY(iBas)%CHARGES(I,IT)
!    BASINFO%ATOMTYPE(atomtype+IT)%Charge = CHARGES(IT)
   ENDDO
   CALL READ_BASISSET_FILE_AND_BUILD_BASISINFO(LUPRI,IPRINT,LUBAS,&
        &CONTRACTED,STRING(IBAS2)%p,IAUG,POLFUN,IPOLST,BASINFO,&
        &CHARGES,k,atomtype,SINGLESEGMENT,DOSPHERICAL,POINTCHARGES)
   call mem_dealloc(CHARGES)
  ENDIF
  call mem_dealloc(POINTCHARGES)
  CALL LSCLOSE(LUBAS,'KEEP')
  atomtype = atomtype+k
 ENDIF
ENDDO
BASINFO%nChargeindex = 0  
CALL ATTACH_Chargeindex_IDTYPE(MOLECULE,BASINFO,BASISSETLIBRARY(iBas))
CALL DETERMINE_NBAST(MOLECULE,BASINFO,DOSPHERICAL,.FALSE.)

IF(IPRINT .GT. 5) CALL PRINT_BASISSETINFO(LUPRI,BASINFO)

!call mem_dealloc(BINDEXES)

CALL DETERMINE_FAMILYTYPEBASISSET(lupri,iprint,BASINFO)

CALL DETERMINE_GENERALCONTRACTED(lupri,iprint,BASINFO,GCONT)
BASINFO%GCONT = GCONT
DEALLOCATE(BASISDIR)
DEALLOCATE(STRING)
call mem_dealloc(LEN_BASISDIR)
END SUBROUTINE Build_BASIS

!> \brief determine if the different atoms have a family type basisset
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE DETERMINE_FAMILYTYPEBASISSET(LUPRI,IPRINT,BASINFO)
implicit none
!> the logical unit number for the output file
INTEGER :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER :: IPRINT
!> contains all info about the basisset 
TYPE(BASISSETINFO)        :: BASINFO
!
INTEGER :: type,B1,J1,B2,J2,nprim1,nprim2,i,iprim
real(realk) :: exp1,exp2
LOGICAL :: FAMILY

do type = 1,BASINFO%nATOMTYPES
 FAMILY = .FALSE.
 DO B1=1,BASINFO%ATOMTYPE(type)%nAngmom
  DO J1=1,BASINFO%ATOMTYPE(type)%SHELL(B1)%nsegments
   nprim1=BASINFO%ATOMTYPE(type)%SHELL(B1)%segment(J1)%nrow
   DO B2=B1+1,BASINFO%ATOMTYPE(type)%nAngmom
    DO J2=1,BASINFO%ATOMTYPE(type)%SHELL(B2)%nsegments
     nprim2=BASINFO%ATOMTYPE(type)%SHELL(B2)%segment(J2)%nrow
     IF(nprim1.EQ.nprim2)THEN
      i = 0
      do iprim = 1,nprim1
       exp1 =BASINFO%ATOMTYPE(type)%SHELL(B1)%segment(J1)%Exponents(iprim)
       exp2 =BASINFO%ATOMTYPE(type)%SHELL(B2)%segment(J2)%Exponents(iprim)
       if(ABS(exp1-exp2) < ExpThr) i = i+1
      enddo
      if(i.EQ.nprim1)FAMILY=.TRUE.
     endif
     IF(FAMILY)EXIT
    enddo
    IF(FAMILY)EXIT
   ENDDO
   IF(FAMILY)EXIT
  enddo
  IF(FAMILY)EXIT
 ENDDO
 BASINFO%ATOMTYPE(type)%FAMILY = FAMILY
enddo

END SUBROUTINE DETERMINE_FAMILYTYPEBASISSET

!> \brief determine if the basis set is general contracted
!> \author T. Kjaergaard
!> \date 2011
SUBROUTINE DETERMINE_GENERALCONTRACTED(LUPRI,IPRINT,BASINFO,GCONT)
implicit none
!> the logical unit number for the output file
INTEGER,intent(in) :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER,intent(in) :: IPRINT
!> contains all info about the basisset 
TYPE(BASISSETINFO),intent(in) :: BASINFO
!> general contracted basis ?
LOGICAL,intent(inout) :: GCONT
!
INTEGER :: type,B1

GCONT = .TRUE.
do type = 1,BASINFO%nATOMTYPES
 DO B1=1,BASINFO%ATOMTYPE(type)%nAngmom
  IF(BASINFO%ATOMTYPE(type)%SHELL(B1)%nprim.NE.0)THEN
   IF(BASINFO%ATOMTYPE(type)%SHELL(B1)%segment(1)%ncol.NE.&
        & BASINFO%ATOMTYPE(type)%SHELL(B1)%norb)THEN
      GCONT = .FALSE.
      RETURN
   ENDIF
  ENDIF
 ENDDO
enddo

END SUBROUTINE DETERMINE_GENERALCONTRACTED

!> \brief attach charge index to basisinfo
!> \author T. Kjaergaard
!> \date 2010
!>
!> BASINFO%Labelindex = 0 :
!> indicate charge based indexing which means that all atoms 
!> share the same basisset like 6-31G or cc-pVTZ - this is the case when 
!> BASIS is used in MOLECULE.INP, This means that the chargeindex is used
!> to acces the given ATOMTYPE
!> BASINFO%Labelindex = n:
!> for n>0 indicate MoleculeSpecific ordering which means that 
!> the molecule%ATOM(iatom)%IDtype(n) determines which ATOMTYPE the given atom has
!> labelindex=1 is for regularMoleculeSpecific ordering 
!> labelindex=2 is for auxiliaryMoleculeSpecific ordering
!> labelindex=3 is for cabs MoleculeSpecific ordering
!> labelindex=4 is for JKaux MoleculeSpecific ordering
!> This is the case when ATOMBASIS is used in MOLECULE.INP
!>
SUBROUTINE ATTACH_Chargeindex_IDTYPE(MOLECULE,BASINFO,BASISSETLIBRARY)
implicit none
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO)        :: BASINFO
!> contains all info about the input basis from the molecule input file
TYPE(BASISSETLIBRARYITEM) :: BASISSETLIBRARY
INTEGER                   :: R,itype,J,I,MAXCHARGE,icharge,IR
IF(BASINFO%Labelindex .EQ. 0)THEN
!  labelindex=0 indicate charge based indexing which means that all atoms 
!  share the same basisset like 6-31G or cc-pVTZ - this is the case when 
!  BASIS is used in MOLECULE.INP, This means that the chargeindex is used
!  to acces the given ATOMTYPE
   MAXCHARGE = 0
   DO ITYPE = 1,BASINFO%natomtypes
      ICHARGE = BASINFO%ATOMTYPE(ITYPE)%CHARGE
      MAXCHARGE = MAX(MAXCHARGE,ICHARGE)
   ENDDO
!   call mem_alloc(BASINFO%Chargeindex,MaxCharge)   
   call mem_alloc(BASINFO%Chargeindex,MaxCharge,.TRUE.)
   BASINFO%nChargeindex=MaxCharge   
   BASINFO%Chargeindex = 0
   DO ITYPE = 1,BASINFO%natomtypes
      ICHARGE = BASINFO%ATOMTYPE(ITYPE)%CHARGE
      BASINFO%ChargeIndex(ICHARGE) = ITYPE
   ENDDO
ELSE
!  labelindex=n for n>0 indicate MoleculeSpecific ordering which means that 
!  the molecule%ATOM(iatom)%IDtype(n) determines which ATOMTYPE the given atom has
!  labelindex=1 is for regularMoleculeSpecific ordering 
!  labelindex=2 is for auxiliaryMoleculeSpecific ordering
!  labelindex=3 is for cabsMoleculeSpecific ordering
!  labelindex=4 is for JKauxMoleculeSpecific ordering
!  This is the case when ATOMBASIS is used in MOLECULE.INP
   R = BASINFO%Labelindex
   DO I=1,MOLECULE%natoms
      MOLECULE%ATOM(I)%IDtype(R) = 0
      IR = MOLECULE%ATOM(I)%basisindex(R)
      IF(MOLECULE%ATOM(I)%pointcharge)THEN
         MOLECULE%ATOM(I)%IDtype(R) = -56
      ELSE       
         DO J=1,BASINFO%natomtypes
!            print*,'iatomtype=',J,'basisindex=IR=',IR,'IATOM=',I
!            print*,'BASINFO%ATOMTYPE(J)%NAME',BASINFO%ATOMTYPE(J)%NAME
!            print*,'BASISSETLIBRARY%BASISSETNAME(IR)',BASISSETLIBRARY%BASISSETNAME(IR)
!            print*,'BASINFO%ATOMTYPE(J)%NAME .EQ.BASISSETLIBRARY%BASISSETNAME(IR)',&
!                 &BASINFO%ATOMTYPE(J)%NAME .EQ.BASISSETLIBRARY%BASISSETNAME(IR)
            IF(BASINFO%ATOMTYPE(J)%NAME .EQ. &
                 &BASISSETLIBRARY%BASISSETNAME(IR)) THEN
!               print*,'BASINFO%ATOMTYPE(J)%Charge',BASINFO%ATOMTYPE(J)%Charge
!               print*,'INT(Molecule%Atom(I)%Charge)',INT(Molecule%Atom(I)%Charge)
               IF(BASINFO%ATOMTYPE(J)%Charge .EQ. INT(Molecule%Atom(I)%Charge))THEN
                  MOLECULE%ATOM(I)%IDtype(R) = J
!                  print*,'MOLECULE%ATOM(I)%IDtype(R) = J = ',J
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      IF(MOLECULE%ATOM(I)%IDtype(R).EQ. 0)&
           & CALL LSQUIT('ERROR IN IDtype assignment',-1)
   ENDDO
ENDIF

END SUBROUTINE ATTACH_CHARGEINDEX_IDTYPE

!> \brief get the path to the basisset library
!> \author T. Kjaergaard
!> \date 2010
!> Modified by S. Reine to handle multiple library directories NBASDIR, Feb 3, 2014
SUBROUTINE GET_BASISSET_LIBS(LUPRI,BASDIR,LENBAS,NBASDIR,doprint)
!****************************************************************
!* OBTAIN A STRING CONTAINING THE BASISSET LIBRARY AND LENGTH=LENBAS 
!****************************************************************
implicit none
TYPE(CAR200),pointer :: BASDIR(:)
CHARACTER(len=1000) :: BASISDIR
INTEGER,pointer  :: LENBAS(:)
INTEGER          :: IDUMMY,IERR,I,LUPRI,NBASDIR,IPOS,IBAS,IPOSNEW
logical          :: doprint,lastbasis

DO I = 1,1000
  BASISDIR(I:I)=' '
ENDDO

#if defined(SYS_AIX)||defined(SYS_CRAY)||defined (SYS_LINUX)
      CALL GETENV ('BASDIR',BASISDIR)
!CALL GET_ENVIRONMENT_VARIABLE(NAME[, VALUE, LENGTH, STATUS, TRIM_NAME) 
!NAME 	Shall be a scalar of type CHARACTER and of default kind.
!VALUE 	(Optional) Shall be a scalar of type CHARACTER and of default kind.
!LENGTH 	(Optional) Shall be a scalar of type INTEGER and of default kind.
!STATUS 	(Optional) Shall be a scalar of type INTEGER and of default kind.
!TRIM_NAME 	(Optional) Shall be a scalar of type LOGICAL and of default kind. 
#endif

IF (BASISDIR(1:1) .NE. '/') THEN

#if defined (VAR_ABSOFT)
         BASISDIR="Insert basis-directory path here"
!        Example:
!        BASDIR="/Users/tkjaer/Dalton/basis/"
#else
!Warning using gfortran it may put in the value of INSTALL_BASDIR and then if the 
!        resulting fortran line is longer is longer then 132 characters it will 
!        truncate the line.
!        the following line therefor restrict the INSTALL_BASDIR size
!                                                         BASISDIR=INSTALL_BASDIR
!        the following line is better
!        BASISDIR=INSTALL_BASDIR
!        while the following line make sure that INSTALL_BASDIR can be 131 long
BASISDIR=&
&INSTALL_BASDIR
#endif
     IF (doprint) WRITE(LUPRI,'(A,A)') ' Default basis set library used:',TRIM(BASISDIR)
ELSE
     IF (doprint) WRITE(LUPRI,'(A,A)') ' User supplied basis set directory :',TRIM(BASISDIR)
END IF

NBASDIR = 1
IPOS    = 1
lastbasis = .FALSE.
DO WHILE (.NOT.lastbasis)
  IPOS = INDEX(BASISDIR(IPOS:1000),':')
  IF (IPOS.EQ.0) THEN
    lastbasis = .true.
  ELSE
    IPOS = IPOS + 1
    NBASDIR = NBASDIR + 1
  ENDIF
ENDDO

ALLOCATE(BASDIR(NBASDIR))
call mem_alloc(LENBAS,NBASDIR)

IF (doprint) write(lupri,'(A,I3)') 'Number of basis set directories in BASDIR is',NBASDIR

IPOS    = 1
DO IBAS=1,NBASDIR
  DO I=1,200
    BASDIR(IBAS)%p(I:I) = ' '
  ENDDO
  IPOSNEW = INDEX(BASISDIR(IPOS:1000),':')
  IF (IPOSNEW.EQ.0) THEN
    IPOSNEW=MIN(IPOS+200,1001)
    LENBAS(IBAS) = 0
    LENBAS(IBAS) = INDEX(BASISDIR(IPOS:IPOSNEW),' ') - 1
  ELSE
    LENBAS(IBAS) = IPOSNEW - IPOS
  ENDIF
  BASDIR(IBAS)%p(1:LENBAS(IBAS)) = BASISDIR(IPOS:IPOSNEW-1)
  IF (BASDIR(IBAS)%p(LENBAS(IBAS):LENBAS(IBAS)) .NE. '/') THEN
    LENBAS(IBAS) = LENBAS(IBAS) + 1
    BASDIR(IBAS)%p(LENBAS(IBAS):LENBAS(IBAS)) = '/'
  END IF
  IPOS = IPOSNEW+1
  IF (doprint) WRITE(LUPRI,'(3X,A)') BASDIR(IBAS)%p
ENDDO
END SUBROUTINE GET_BASISSET_LIBS

!> \brief DETERMINE NUMBER OF UNIQUE CHARGES
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE UNIQUE_CHARGES(LUPRI,BASISSETLIBRARY,CHARGES,k)
implicit none
INTEGER  :: LUPRI,k,I,IT,J
real(realk) :: CHARGES(:)
TYPE(BASISSETLIBRARYITEM) :: BASISSETLIBRARY
LOGICAL  :: NEWCHARGE
DO I=1,BASISSETLIBRARY%nbasissets
   DO J=1,BASISSETLIBRARY%nCharges(I)
      IF(I==1 .AND. J==1)THEN
         CHARGES(1)=BASISSETLIBRARY%Charges(I,J)
         k=1
      ELSE
         NEWCHARGE=.TRUE.
         DO IT=1,k
            IF(ABS(BASISSETLIBRARY%Charges(I,J)-CHARGES(IT)).LT.1.0E-12_realk)THEN
               NEWCHARGE=.FALSE.
            ENDIF
         ENDDO
         IF(NEWCHARGE) THEN
            k=k+1
            CHARGES(k)=BASISSETLIBRARY%Charges(I,J)
         ENDIF
      ENDIF
   ENDDO
ENDDO

END SUBROUTINE UNIQUE_CHARGES

!> \brief PRODUCES THE AUGMENTATIONLEVEL AND NEW BASISNAME WITHOUT AUGMENTATION
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE DETERMINE_AUGMENTATION(LUPRI,BASNAM,IAUG)
implicit none
CHARACTER(len=80)  :: BASNAM,BASTMP
INTEGER            :: IAUG,LUPRI

IF (BASNAM(2:9) .EQ. 'aug-cc-p') THEN
  BASTMP = BASNAM
  BASNAM(1:79) = BASTMP(2:80)
  IF (BASTMP(1:1) .EQ. 'd') THEN
    IAUG = 1
  ELSEIF (BASTMP(1:1) .EQ. 't') THEN
    IAUG = 2
  ELSE IF (BASTMP(1:1) .EQ. 'q') THEN
    IAUG = 3
  ELSE
    WRITE (LUPRI,'(/A1,A)') BASTMP(1:1),' is an unknown augmentation level'
    CALL LSQUIT('Illegal augmentation level in LINSCA',lupri)
  END IF
ELSEIF (BASNAM(1:7) .EQ. 'aug-ecp') THEN
  BASTMP = BASNAM
  BASNAM(1:75) = BASTMP(5:80)
  IAUG = 1
ELSEIF (BASNAM(2:8) .EQ. 'aug-ecp') THEN
  BASTMP = BASNAM
  BASNAM(1:74) = BASTMP(6:80)
  IF (BASTMP(1:1) .EQ. 'd') THEN
    IAUG = 1
  ELSE IF (BASTMP(1:1) .EQ. 't') THEN
    IAUG = 2
  ELSE IF (BASTMP(1:1) .EQ. 'q') THEN
    IAUG = 3
  ELSE
    WRITE (LUPRI,'(/A1,A)') BASTMP(1:1),' is an unknown augmentation level'
    CALL LSQUIT('Illegal augmentation level in LINSCA',lupri)
  END IF
ELSE
  IAUG=0
ENDIF

END SUBROUTINE DETERMINE_AUGMENTATION

!> \brief PRODUCES THE POLARIZATIONLEVEL = IPOLST AND LOGICAL POLFUN 
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE DETERMINE_POLARIZATION(LUPRI,BASNAM,POLFUN,IPOLST)
implicit none
CHARACTER(len=80)  :: BASNAM
INTEGER            :: IPOLST,LUPRI
LOGICAL            :: POLFUN

IPOLST = INDEX(BASNAM,'Pol')
IF (IPOLST .GT. 0) THEN
   POLFUN = .TRUE.
   IPOLST = IPOLST + 3
   CALL LSQUIT('Polarization functions are not yet implemented',lupri)
ELSE
   POLFUN = .FALSE.
END IF
END SUBROUTINE DETERMINE_POLARIZATION

!> \brief read the basisset file and build the basisinfo
!> \author T. Kjaergaard
!> \date 2010
!>
!> THIS IS THE ROUTINE WHICH READS THE BASISSET FILE
!> WITH LOGICAL UNIT NUMBER LUBAS AND THEN BUILDS THE BASISINFO
!> Written by Thomas Kjaergaard  
!>
!> THE RECIPE IS
!> 1 SORT THE CHARGES - by simple bubblesort
!> 2 RUN OVER ALL DIFFERENT CHARGES   
!> 3 FIND POSITION OF THIS CHARGE ON FILE 
!> 4 READ # Exponents = nprim
!>   AND # Primitive orbitals = nOrbitals FOR
!>   EACH ANGULAR MOMENTUM. WHEN IT REACH THE END OF THE ELEMENT
!>   NewElement is set to .TRUE. AND ELEMENT IS THE NEXT ELEMENT  
!> 5 THE EXPONENTS AND CONTRACTIONCOEFFICIENTS ARE READ IN THE
!>   SUBROUTINE 'READ_COEFFICIENT_AND_EXPONENTS'
!> 6 UNLESS UNCONTRACTED IS SPECIFIED THE SUBROUTINE 
!>   'ANALYSE_CONTRACTIONMATRIX' IS CALLED TO FIND THE NONZERO 
!>   BLOCKS OF THE CONTRACTIONCOEFFICIENTMATRIX   
!>
!> S.Reine 2010-02-23 segmentedFormat added
!>
!>
SUBROUTINE READ_BASISSET_FILE_AND_BUILD_BASISINFO(LUPRI,IPRINT,LUBAS,&
     &CONTRACTED,BASISSETNAME,IAUG,POLFUN,IPOLST,&
     &BASINFO,CHARGE,L,atomtype,SINGLESEGMENT,dospherical,POINTCHARGES)
implicit none
!> the logical unit number for the output file
INTEGER                   :: LUPRI
!> the printlevel integer, determining how much output should be generated
INTEGER                   :: IPRINT
!> the logical unit number for the basisset file
INTEGER                   :: LUBAS
!> if the basis should be contracted
LOGICAL                   :: CONTRACTED
!> basisset name
CHARACTER(len=200)   :: BASISSETNAME
!> IAUG      : AUGMENTATION LEVEL  
INTEGER  :: IAUG
!> POLFUN  : TRUE IF POLARIZATION FUNCTIONS ARE SPECIFIED  
LOGICAL  :: POLFUN
!> IPOLST 
INTEGER  :: IPOLST
!> contains all info about the basisset (exponents,number of primitives,...)
TYPE(BASISSETINFO)  :: BASINFO
!> CHARGE    : A VECTOR WITH LENGTH L CONTAINING THE CHARGES FOR 
!>           : WHICH THE BASISSET SHOULD BE USED.
!INTEGER         :: CHARGE(L)
REAL(REALK) :: CHARGE(L)
LOGICAL :: POINTCHARGES(L)
!> size of CHARGE 
INTEGER             :: L
!>atomtype
INTEGER :: atomtype
!> IF NO SEGMENTATION SHOULD BE DONE OF THE basis
LOGICAL :: SINGLESEGMENT
!> if you want to use spherical basis function (default) or cartesian
LOGICAL            :: dospherical
!
LOGICAL             :: SEARCHING,NewElements
INTEGER             :: nprim, nOrbitals,K,nAngmom
INTEGER             :: ELEMENT,Totalnprim,TotalnOrbitals,AA,I,TESTICHARGE
Integer,pointer :: icharge(:)
TYPE(Lsmatrix)        :: ContractionMatrix
TYPE(Lsmatrix)        :: ContractionMatrixNORM
TYPE(Lsmatrix)        :: Exponents
logical             :: segmentedFormat,HFORMAT
CHARACTER(len=200)  :: STRING
IF(IPOLST.GT. 0)call lsquit('polarization functions not implemented',lupri)
SEARCHING=.TRUE.
ELEMENT=0

call mem_alloc(ICHARGE,L)
DO I=1,L !number of unique charges in basisset
   IF(POINTCHARGES(I))CYCLE
   ICHARGE(I) = NINT(CHARGE(I))
   TESTICHARGE = NINT(CHARGE(I) + 1E-8_realk)
   IF(abs(charge(I) - TESTICHARGE).GT. 1E-7_realk)THEN
      CALL LSQUIT('Non-integer charge in READ_BASISSET_FILE_AND_BUILD_BASISINFO',lupri)
   ENDIF
ENDDO

! Sort the Charges so that only read file one time
Call bubbleSort(ICHARGE,L,CHARGE,POINTCHARGES)
HFORMAT = .FALSE.
DO I=1,L !number of unique charges in basisset
  IF(POINTCHARGES(I))THEN
     BASINFO%ATOMTYPE(atomtype+I)%nAngmom=0
     BASINFO%ATOMTYPE(atomtype+I)%Charge=0
     CYCLE
  ENDIF
  IF(IPRINT .GT. 5) WRITE(LUPRI,*)'FIND CHARGE',ICHARGE(I)
  !SEARCH TO FIND ELEMENT WITH CHARGE(I) 
  CALL FIND_POSITION_ON_FILE(LUPRI,LUBAS,ICHARGE(I),BASISSETNAME,ELEMENT,HFORMAT,STRING) 
  BASINFO%ATOMTYPE(atomtype+I)%Charge=NINT(CHARGE(I))
  BASINFO%ATOMTYPE(atomtype+I)%nAngmom=0
  NewElements=.FALSE.
  nAngmom=0
  TotalnOrbitals=0   !Total number of contracted orbitals and
  Totalnprim=0 !Total number of primitives for this charge
  
  DO
  !READ number of primitives = number of exponents
  !And number of orbitals = number of contracted basisfunctions
  CALL TEST_IF_NEXT_ANGMOM_OR_NEW_ELM(LUPRI,LUBAS, nprim, nOrbitals,&
                                      & NewElements,ELEMENT,HFORMAT,STRING)
  IF(NewElements) EXIT

  segmentedFormat = nOrbitals.EQ. 0
  IF (segmentedFormat) nOrbitals = nprim
  nAngmom=nAngmom+1
  IF(nAngmom .GT. maxAOangmom)THEN
       WRITE(LUPRI,*)'You have at least one atom with many different orbitals&
            & - more than just the normal s, p, d, f, g, h, i, j, k and l&
            & orbitals so you need to increase the maxAOangmom parameter&
            & in TYPE-DEF.f90 file and recompile. It should be even to or&
            & greater than the number of different orbitals you need.&
            & Thomas Kjaergaard'
       CALL LSQUIT('Increase maxAOangmom in TYPE-DEF.f90',lupri)
  ENDIF

  IF(nprim .NE. 0)THEN
    IF(dospherical)THEN !default
       ! 2*(nAngmom)-1 is the degenracy for given angularmomentum 
       Totalnprim=(nprim+IAUG)*(2*(nAngmom)-1)+Totalnprim
       TotalnOrbitals=(nOrbitals+IAUG)*(2*(nAngmom)-1)+TotalnOrbitals
    ELSE
       !for cartesian the degeneracy is nAngmom*(nAngmom+1)/2 for given angularmomentum 
       Totalnprim=(nprim+IAUG)*(nAngmom*(nAngmom+1)/2)+Totalnprim
       TotalnOrbitals=(nOrbitals+IAUG)*(nAngmom*(nAngmom+1)/2)+TotalnOrbitals
    ENDIF
    IF ((.NOT.CONTRACTED).OR.segmentedFormat) THEN 
      BASINFO%ATOMTYPE(atomtype+I)%SHELL(nAngmom)%norb=nprim+IAUG
      BASINFO%ATOMTYPE(atomtype+I)%SHELL(nAngmom)%nprim=nprim+IAUG
      !If the basis is uncontracted or forced uncontracted, we know that each
      !block is a 1 by 1 matrix and we can imidiate initialise the exponents 
      !and contraction matrices
      CALL INIT_BASISSETINFO_ContractionM(BasInfo,atomtype+I,&
                                               &nAngmom,nprim+IAUG)
      DO K=1,nprim+IAUG
        CALL INIT_BASISSETINFO_elms(BasInfo,atomtype+I,nAngmom,K,1,1)
      ENDDO
      !WHEN THE .UNCONT OPTION IS USED THE ORBITALS ARE NORMALIZED
      !IN READ_COEFFICIENT_AND_EXPONENTS
      CALL READ_COEFFICIENT_AND_EXPONENTS(LUPRI,IPRINT,LUBAS,BASINFO,&
           &atomtype+I,nAngmom,nprim,nOrbitals,IAUG,POLFUN,CONTRACTED,&
           &ContractionMatrix,Exponents,segmentedFormat)
      !maybe ContractionMatrix optional
    ELSE
      BASINFO%ATOMTYPE(atomtype+I)%SHELL(nAngmom)%norb=nOrbitals+IAUG
      BASINFO%ATOMTYPE(atomtype+I)%SHELL(nAngmom)%nprim=nprim+IAUG
      !As we do not yet know how many blocks we need we initialise a 
      !full Matrix and a full vector and read the exponents and coefficients
      !afterwards we then call the ANALYSE_CONTRACTIONMATRIX
      !To determine the number of blocks and then initiate and construct the 
      !proper blocks in BASINFO
      CALL LSMAT_DENSE_INIT(ContractionMatrix,nprim+IAUG,nOrbitals+IAUG)
      CALL LSMAT_DENSE_INIT(Exponents,nprim+IAUG,1)
      CALL LSMAT_DENSE_ZERO(ContractionMatrix)
      CALL READ_COEFFICIENT_AND_EXPONENTS(LUPRI,IPRINT,LUBAS,BASINFO,&
           &atomtype+I,nAngmom,nprim,nOrbitals,IAUG,POLFUN,CONTRACTED,&
           &ContractionMatrix,Exponents,segmentedFormat)
      !FIXED PROBLEMS WITH AUGMENTATION AND DIMENSIONS
      !AFTER READ_COEFFICIENT_AND_EXPONENTS 
      ! nprim+iaug -> nprim  AND  nOrbitals+iaug -> nOrbitals
      CALL LSMAT_DENSE_INIT(ContractionMatrixNorm,nprim+IAUG,nOrbitals+IAUG)
      CALL LSMAT_DENSE_ZERO(ContractionMatrixNorm)
      CALL NORMALIZE_ORBITALS(LUPRI,IPRINT,nangmom,nOrbitals,nprim,&
           &Exponents,ContractionMatrix,ContractionMatrixNorm)

      IF(SINGLESEGMENT)THEN
         WRITE(LUPRI,*)'NO SEGMENTATION IS USED'
         CALL INIT_BASISSETINFO_ContractionM(BasInfo,atomtype+I,&
                                                                &nAngmom,1)
         CALL INIT_BASISSETINFO_elms(BasInfo,atomtype+I,nAngmom,&
                                                   &1,nprim,nOrbitals)
         DO AA=1,nprim*nOrbitals
            BASINFO%ATOMTYPE(atomtype+I)%SHELL(nAngmom)%segment(1)&
                 &%elms(AA)=ContractionMatrixNORM%elms(AA)
         ENDDO
         DO AA=1,nprim*nOrbitals
            BASINFO%ATOMTYPE(atomtype+I)%SHELL(nAngmom)%segment(1)&
                 &%UCCelms(AA)=ContractionMatrix%elms(AA)
         ENDDO
         DO AA=1,nprim
            BASINFO%ATOMTYPE(atomtype+I)%SHELL(nAngmom)%&
                 &segment(1)%Exponents(AA)=Exponents%elms(AA)
         ENDDO
      ELSE !DEFAULT
         CALL ANALYSE_CONTRACTIONMATRIX(LUPRI,IPRINT,BASINFO,atomtype+I,&
              &nAngmom,nprim,nOrbitals,ContractionMatrix,ContractionMatrixNORM,Exponents)   
      ENDIF
      
      CALL LSMAT_DENSE_FREE(ContractionMatrix)
      CALL LSMAT_DENSE_FREE(ContractionMatrixNORM)
      CALL LSMAT_DENSE_FREE(Exponents)
    ENDIF

  ELSE
     BASINFO%ATOMTYPE(atomtype+I)%SHELL(nAngmom)%norb=0
     BASINFO%ATOMTYPE(atomtype+I)%SHELL(nAngmom)%nprim=0
     BASINFO%ATOMTYPE(atomtype+I)%SHELL(nAngmom)%nsegments=0
  ENDIF
  STRING(1:12) = '            '
  ENDDO

  BASINFO%ATOMTYPE(atomtype+I)%Totnprim=Totalnprim
  IF(CONTRACTED)THEN !DEFAULT
     BASINFO%ATOMTYPE(atomtype+I)%Totnorb=TotalnOrbitals
  ELSE
     BASINFO%ATOMTYPE(atomtype+I)%Totnorb=Totalnprim
  ENDIF
ENDDO
call mem_dealloc(ICHARGE)
END SUBROUTINE READ_BASISSET_FILE_AND_BUILD_BASISINFO

!> \brief basic bubblesort
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE bubbleSort(B,N,C,P)
!**********************************************************************
!* SORTING ROUTINE 
!*********************************************************************
implicit none
integer :: N, newB, B(N), k, j
real(realk) :: C(N),newC
logical :: P(N),newP
do k=1, N
   do j=1, N-1
    if (B(j+1)<B(j)) then
     newB=B(j)   
     B(j)=B(j+1)
     B(j+1)=newB
     newC=C(j)   
     C(j)=C(j+1)
     C(j+1)=newC
     newP=P(j)   
     P(j)=P(j+1)
     P(j+1)=newP
    endif
   enddo
 enddo
end subroutine bubbleSort

!> \brief test if the entire line is blank
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE TEST_BLANK_LINE(BLANK, STRING)
!**********************************************************************
!* TEST IF THE STRING IS BLANK
!*********************************************************************
implicit none
CHARACTER*(*) :: STRING
LOGICAL       :: BLANK
INTEGER       :: J

BLANK = .TRUE.
DO J = 1, LEN(STRING)
   BLANK = BLANK .AND. (STRING(J:J) .EQ. ' ')
ENDDO

END SUBROUTINE TEST_BLANK_LINE

subroutine FIND_BASIS_IN_MOLFILE(LUPRI,LUBAS) 
implicit none
integer,intent(in) :: lupri,lubas
!
logical :: SEARCHING
INTEGER :: IOS
CHARACTER(len=200)  :: STRING
SEARCHING =.TRUE.
DO WHILE(SEARCHING) 
   READ(LUBAS,'(A200)', IOSTAT=ios) STRING
   IF(STRING(1:17).EQ.'USERDEFINED BASIS')THEN
      SEARCHING=.FALSE.
   ENDIF
ENDDO
end subroutine FIND_BASIS_IN_MOLFILE

!> \brief find the correct position on file for the given charge
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE FIND_POSITION_ON_FILE(LUPRI,LUBAS,CHARGE,BASISSETNAME,ELEMENT,HFORMAT,STRING) 
implicit none 
CHARACTER(len=200)  :: STRING
CHARACTER(len=1)    :: SIGN
INTEGER             :: CHARGE,NUCLEARCHARGE,ios,ELEMENT,LUBAS,LUPRI,IPOS
LOGICAL             :: SEARCHING,SEARCHING_H,HFORMAT
CHARACTER(len=200)   :: BASISSETNAME
CHARACTER(len=10)   :: atname(74)
atname = (/'HYDROGEN  ','HELIUM    ','LITHIUM   ',&
     & 'BERYLLIUM ','BORON     ','CARBON    ','NITROGEN  ',&
     & 'OXYGEN    ','FLUORINE  ','NEON      ','SODIUM    ',&
     & 'MAGNESIUM ','ALUMINUM  ','SILICON   ','PHOSPHORUS',&
     & 'SULFUR    ','CHLORINE  ','ARGON     ','POTASSIUM ',&
     & 'CALCIUM   ','SCANDIUM  ','TITANIUM  ','VANADIUM  ',&
     & 'CHROMIUM  ','MANGANESE ','IRON      ','COBALT    ',&
     & 'NICKEL    ','COPPER    ','ZINC      ','GALLIUM   ',&
     & 'GERMANIUM ','ARSENIC   ','SELENIUM  ','BROMINE   ',&
     & 'KRYPTON   ','RUBIDIUM  ','STRONTIUM ','YTTRIUM   ',&
     & 'ZIRCONIUM ','NIOBIUM   ','MOLYBDENUM','TECHNETIUM',&
     & 'RUTHENIUM ','RHODIUM   ','PALLADIUM ','SILVER    ',&
     & 'CADMIUM   ','INDIUM    ','TIN       ','ANTIMONY  ',&
     & 'TELLERIUM ','IODINE    ','XENON     ','CESIUM    ',&
     & 'BARIUM    ','LANTHANUM ','HAFNIUM   ','TANTALUM  ',&
     & 'WOLFFRAM  ','RHENIUM   ','OSMIUM    ','IRIDIUM   ',&
     & 'PLATINUM  ','GOLD      ','MERCURY   ','THALLIUM  ',&
     & 'LEAD      ','BISMUTH   ','POLONIUM  ','ASTATINE  ',&
     & 'RADON     ','FRANCIUM  ','RADIUM    '/)
! Searching the file for the element.
SEARCHING=.TRUE.
IF((ELEMENT .EQ. 0) .OR. (ELEMENT.NE.CHARGE))THEN
  DO WHILE(SEARCHING) 
    READ(LUBAS,'(A88)', IOSTAT=ios) STRING
    IF(ios /= 0)THEN
      WRITE (LUPRI,'(/I3,2A)') CHARGE&
      & ,' is an unsupported element for basis ',BASISSETNAME
      WRITE (LUPRI,'(A)') 'You need to choose a basis set that support this element'
      WRITE (LUPRI,'(A)') 'You can look at EMSL to find a suitable basis set that support this element'
      WRITE (LUPRI,'(A)') 'If you download a new basis set from EMSL please read the EMSL section in the manual!'
      WRITE (LUPRI,'(A)') 'Note that LSDALTON do NOT support Effective Core Potentials (ECP)'
      CALL LSQUIT('Unsupported element in basis set, choose a proper basis set',lupri)
    ELSE
      READ (STRING, '(A1)') SIGN
      IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
        READ (STRING, '(A1, I4)') SIGN, NUCLEARCHARGE
        IF (CHARGE .EQ. NUCLEARCHARGE) THEN
          SEARCHING=.FALSE.
        ENDIF     
      ELSE
         IPOS = INDEX(STRING,atname(charge))
         IF (IPOS .NE. 0) THEN   
            HFORMAT = .TRUE.
            SEARCHING=.FALSE.
            SEARCHING_H=.TRUE.
            DO WHILE(SEARCHING_H) 
               READ(LUBAS,'(A88)', IOSTAT=ios) STRING
               READ (STRING, '(A1)') SIGN
               IF (SIGN .EQ. 'H')THEN
                  SEARCHING_H = .FALSE.
               ELSEIF ((SIGN .EQ. '!') .OR. (SIGN .EQ. '$'))THEN
                  !comment line
               ELSE
                  print*,'error in basis set input expected a comment'
                  print*,'or a line starting with H, but instead got the string'
                  print*,STRING
                  call lsquit('error in basis set input ',-1)
               ENDIF
            ENDDO
         ENDIF
      ENDIF
    ENDIF
  ENDDO
ELSE ! ELEMENT .EQ. CHARGE
!ALREADY AT THE CORRECT SPOT IN THE FILE
ENDIF

END SUBROUTINE FIND_POSITION_ON_FILE

!> \brief TEST IF THE NEXT LINE IS THE NEXT ANGULAR MOMENTUM OR WE HAVE REACH THE NEXT ELEMENT
!> \author T. Kjaergaard
!> \date 2010
!> 
!> TEST IF THE NEXT LINE IS THE NEXT ANGULAR MOMENTUM OR WE HAVE REACH THE NEXT ELEMENT
!> IN CASE OF FINDING THE NEXT ANGULAR MOMENTUM WE READ 
!> NUMBER OF EXPONENTS AND NUMBER OF PRIMITIVE ORBITALS = nExponential, nOrbitals
!> IF WE INSTEAD FIND A NEW ELEMENT WE SET NEWEL=.TRUE.
!> AND ELEMENT EQUAL TO THE NEW ELEMENT 
!> 
SUBROUTINE TEST_IF_NEXT_ANGMOM_OR_NEW_ELM(LUPRI,LUBAS, nExponential, nOrbitals, NEWEL,ELEMENT,HFORMAT,STRING)
implicit none
INTEGER             :: LUBAS,nExponential, nOrbitals,ELEMENT,ios,LUPRI,IPOS,IATOM
LOGICAL             :: BLANK,NEWEL,HFORMAT,START
CHARACTER(len=200)  :: STRING
CHARACTER(len=1)    :: SIGN
CHARACTER(len=10)   :: atname(74)
atname = (/'HYDROGEN  ','HELIUM    ','LITHIUM   ',&
     & 'BERYLLIUM ','BORON     ','CARBON    ','NITROGEN  ',&
     & 'OXYGEN    ','FLUORINE  ','NEON      ','SODIUM    ',&
     & 'MAGNESIUM ','ALUMINUM  ','SILICON   ','PHOSPHORUS',&
     & 'SULFUR    ','CHLORINE  ','ARGON     ','POTASSIUM ',&
     & 'CALCIUM   ','SCANDIUM  ','TITANIUM  ','VANADIUM  ',&
     & 'CHROMIUM  ','MANGANESE ','IRON      ','COBALT    ',&
     & 'NICKEL    ','COPPER    ','ZINC      ','GALLIUM   ',&
     & 'GERMANIUM ','ARSENIC   ','SELENIUM  ','BROMINE   ',&
     & 'KRYPTON   ','RUBIDIUM  ','STRONTIUM ','YTTRIUM   ',&
     & 'ZIRCONIUM ','NIOBIUM   ','MOLYBDENUM','TECHNETIUM',&
     & 'RUTHENIUM ','RHODIUM   ','PALLADIUM ','SILVER    ',&
     & 'CADMIUM   ','INDIUM    ','TIN       ','ANTIMONY  ',&
     & 'TELLERIUM ','IODINE    ','XENON     ','CESIUM    ',&
     & 'BARIUM    ','LANTHANUM ','HAFNIUM   ','TANTALUM  ',&
     & 'WOLFFRAM  ','RHENIUM   ','OSMIUM    ','IRIDIUM   ',&
     & 'PLATINUM  ','GOLD      ','MERCURY   ','THALLIUM  ',&
     & 'LEAD      ','BISMUTH   ','POLONIUM  ','ASTATINE  ',&
     & 'RADON     ','FRANCIUM  ','RADIUM    '/)
NEWEL = .FALSE.
ELEMENT=0
START=.TRUE.
ios=0
FINPOSLOOP: DO
  IF(START)THEN
     !the first time you enter this subroutine you should only read a newline 
     !if this is not the EMSL H format.
     IF(.NOT.HFORMAT)THEN
        READ(LUBAS,'(A88)', IOSTAT=ios) STRING
     ENDIF
  ELSE
     READ(LUBAS,'(A88)', IOSTAT=ios) STRING
  ENDIF
  START=.FALSE.
  IF(ios /= 0)THEN
     NEWEL=.TRUE.
     nExponential=0
     nOrbitals=0
     EXIT FINPOSLOOP
  ELSE
    READ (STRING, '(A1)') SIGN
    IF (SIGN .EQ. ' ')THEN
      CALL TEST_BLANK_LINE(BLANK, STRING)
      IF (BLANK) THEN
        cycle FINPOSLOOP
      ELSE
        READ (STRING,*) nExponential, nOrbitals
        EXIT FINPOSLOOP
      END IF
    ELSEIF(SIGN .EQ. 'H') THEN
       READ (STRING,*) SIGN,nExponential, nOrbitals
       EXIT FINPOSLOOP
    ELSEIF((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
      READ (STRING, '(A1, I4)') SIGN, ELEMENT
      NEWEL = .TRUE.
      nExponential=0
      nOrbitals=0
      EXIT FINPOSLOOP
    ELSEIF(HFORMAT) THEN
       DO IATOM=1,size(atname)
          IPOS = INDEX(STRING,atname(iatom))
          IF (IPOS .NE. 0) THEN   
             NEWEL = .TRUE.
             nExponential=0
             nOrbitals=0
             EXIT FINPOSLOOP
          ENDIF
       ENDDO
    ELSEIF(SIGN .EQ. '!') THEN !Comment line
       IF(HFORMAT) THEN
          DO IATOM=1,size(atname)
             IPOS = INDEX(STRING,atname(iatom))
             IF (IPOS .NE. 0) THEN   
                NEWEL = .TRUE.
                nExponential=0
                nOrbitals=0
                EXIT FINPOSLOOP
             ENDIF
          ENDDO
       ELSE
          cycle FINPOSLOOP
       ENDIF
    ELSEIF(SIGN .EQ. '$') THEN !Comment line
       IF(HFORMAT) THEN
          DO IATOM=1,size(atname)
             IPOS = INDEX(STRING,atname(iatom))
             IF (IPOS .NE. 0) THEN   
                NEWEL = .TRUE.
                nExponential=0
                nOrbitals=0
                EXIT FINPOSLOOP
             ENDIF
          ENDDO
       ELSE
          cycle FINPOSLOOP
       ENDIF
    ENDIF
  ENDIF
ENDDO FINPOSLOOP

END SUBROUTINE TEST_IF_NEXT_ANGMOM_OR_NEW_ELM

!> \brief READ_COEFFICIENT_AND_EXPONENTS
!> \author T. Kjaergaard
!> \date 2010
!> 
!> WE READ THE NUMBER OF EXPONENTS  =  nExponential, 
!> AND NUMBER OF PRIMITIVE ORBITALS =  nOrbitals
!> THE CONTRACTION COEFFICIENTS CAN BE WRITTEN ON A NUMBER OF LINES 
!> 6 ON THE FIRST LINE AND 7 ON THE REMAINING LINES
!> 
!> AFTERWARDS WE AUGMENT AND ADD POLARIZED FUNCTIONS IF NECCESARY
!> S.Reine 2010-02-23 segmentedFormat added
!>
SUBROUTINE READ_COEFFICIENT_AND_EXPONENTS(LUPRI,IPRINT,LUBAS,BASINFO,&
     &atype,nAng,nprim,nOrbital,IAUG,POLFUN,CONTRACTED,&
     &Contractionmatrix,Exponents,segmentedFormat)
  implicit none
  TYPE(BASISSETINFO)    :: BASINFO
  INTEGER               :: LUBAS,ios,IEX,I,M,LUPRI,IPRINT
  LOGICAL               :: POLFUN,CONTRACTED,segmentedFormat,BLANK
  INTEGER               :: atype,nang,nprim,nOrbital,IAUG,NUMNUMOLD
  INTEGER               :: J,NUMBER_OF_LINES,KNTORB,NUMNUM,KAUG,nNumbers
  CHARACTER(len=200)    :: STRING
  CHARACTER(len=1)      :: SIGN
  real(realk)           :: exmin2,exmin1,PI,Exp,PIPPI
  CHARACTER(len=1)      :: SPDFGH(10)=(/'S','P','D','F','G','H','I','J','K','L'/) 
  TYPE(Lsmatrix)          :: ContractionMatrix  !Not initiated if UNCONTRACTED
  TYPE(Lsmatrix)          :: Exponents          !Not initiated if UNCONTRACTED
  PI=3.14159265358979323846E0_realk
  PIPPI = (0.5E0_realk/PI)**(0.75E0_realk)
  J = 0
  DO WHILE( J .LT. nprim) 
     ! Reading the primitive and contracted coeffecients
     READ(LUBAS, '(A200)', IOSTAT = IOS) STRING
     IF(ios /= 0)THEN
        WRITE (LUPRI,'(2A)') ' Error in basisset file'
        CALL LSQUIT('Error in basisset file',lupri)
     ELSE
        READ (STRING, '(A1)') SIGN
        IF (SIGN .EQ. ' ') THEN
           CALL TEST_BLANK_LINE(BLANK, STRING)
           IF (.NOT. BLANK) THEN
              !We have found a line with a primitive and some coeffecients
              J=J+1
              !LINES_OF_CONTRACTION returns the number of lines 
              !the contraction coeffecients written on.
              CALL determine_nNumbers_in_string(STRING,nNUMBERS)
              IF(nNUMBERS.GT.7)THEN
                 WRITE(lupri,*)'WARNING This basis format violte the old dalton basis format'
                 WRITE(lupri,*)'Which consist of 1 Exponent F16.9 and up to 6 contraction'
                 WRITE(lupri,*)'coefficients on the first line followed by a up to 7 contraction coefficients'
                 WRITE(lupri,*)'on the following lines until the full number of contraction coefficients are given'
                 WRITE(lupri,*)'We will try to this basis set, but this code is not very well testet. TK'
              ENDIF
              IF(nNUMBERS-1.LT.nOrbital)THEN
                 CALL LINES_OF_CONTRACTION(nOrbital,nNUMBERS-1,NUMBER_OF_LINES,segmentedFormat)
              ELSEIF(nNUMBERS-1.EQ.nOrbital)THEN
                 NUMBER_OF_LINES=1
              ELSE
                 CALL LSQUIT('Error in determining number of continuation lines',-1)
              ENDIF
              IF (.NOT.CONTRACTED.OR.segmentedFormat) THEN
                 IF(J .GT. maxBASISsegment)THEN
                    WRITE(LUPRI,*)'J=',J
                    WRITE(LUPRI,*)'nang=',nang
                    WRITE(LUPRI,*)'nang=',nang
                    WRITE(LUPRI,*)'You have alot of segments&
                         & This is okay, but you need to increase the parameter&
                         & maxBASISsegment in the TYPE-DEF.f90 file and recompile.&
                         & Currently the parameter is set to ',maxBASISsegment,&
                         &'Thomas Kjaergaard'
                    CALL LSQUIT('Increase maxBASISsegment in TYPE-DEF.f90 file',lupri)
                 ENDIF
                 !WRITE(LUPRI,*) 'INSIDE UNCONTRACTED'
                 !uncontracted basis set are forced uncontracted with .UNCONT

                 READ (STRING, '(F16.9)') Exp
                 !skip contiuation lines for contraction coeff.
                 DO I = 2,NUMBER_OF_LINES
                    READ (LUBAS,*)
                 END DO
                 BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(J)%Exponents(1) = Exp
                 BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(J)%elms(1)= &
                      &(4*Exp)**(0.5E0_realk*nang+0.25E0_realk)*PIPPI 
              ELSE
                 !         Getting the format for the read statement right.
                 IF (NUMBER_OF_LINES .EQ. 1) THEN
                    KNTORB = nOrbital
                 ELSE
                    KNTORB = nNUMBERS-1
                 END IF
                 READ (STRING,*) Exponents%elms(J),&
                      &(ContractionMatrix%elms(J+(I-1)*(nprim+IAUG)),I = 1, KNTORB)
                 !         Reading the first line with exponents and contractioncoeffecients
                 !     READ (STRING, '(F16.9, 6F12.9)') Exponents%elms(J),&
                 !          &(ContractionMatrix%elms(J+(I-1)*(nprim+IAUG)),I = 1, KNTORB)
                 !         If there are more lines with contraction-coeffecients
                 !         they will be read here.
                 NUMNUM = nNUMBERS-1
                 NUMNUMOLD = nNUMBERS-1
                 DO I=2, NUMBER_OF_LINES
                    CALL determine_nNumbers_in_string(STRING,nNUMBERS)
                    IF(nNUMBERS.GT.7)THEN
                       WRITE(lupri,*)'WARNING This is a continuation line, which violte the old dalton basis format'
                       WRITE(lupri,*)'The continuation line should contain at most 7 contraction coefficients'
                       WRITE(lupri,*)'We will try to this basis set, but this code is not very well testet. TK'           
                    ENDIF
                    NUMNUM = NUMNUM + nNUMBERS
                    KNTORB = MIN(NUMNUM, nOrbital)
                    !Getting the format for the read-stat right.
                    !Making the usual safety-precautions before we read the 
                    !contraction-coeffecients.
                    READ(LUBAS, '(A200)', IOSTAT = IOS) STRING
                    IF(ios /= 0)THEN
                       WRITE (LUPRI,'(2A)') ' Error in basisset file'
                       CALL LSQUIT('Error in basisset file',lupri)
                    ELSE
                       READ (STRING, '(A1)') SIGN
                       IF (SIGN .EQ. ' ') THEN
                          CALL TEST_BLANK_LINE(BLANK, STRING)
                          IF (.NOT. BLANK) THEN
                             !We now have a line with contraction-coeffecients.
                             !         READ (STRING,'(F16.9,6F12.9)')&
                             !              &(Contractionmatrix%elms(J+(M-1)*(nprim+IAUG)),&
                             !              & M = 6 + (I-2)*7 +1, KNTORB)
                             READ (STRING,*)&
                                  &(Contractionmatrix%elms(J+(M-1)*(nprim+IAUG)),&
                                  & M = NUMNUMOLD+1, KNTORB)
                          END IF
                       ELSE
                          cycle
                       ENDIF
                    ENDIF
                    NUMNUMOLD = NUMNUM
                    !just in case NUMBER_OF_LINES was wrong
                    IF(KNTORB.EQ.nOrbital)EXIT
                 ENDDO
              ENDIF
           ENDIF
        ENDIF
     ENDIF
  ENDDO

  ! If this is some kind of augmented basis set, we augment it here, kr-96
  !
  ! Fix jan -01 VB: For certain elements some of the cc basis sets do 
  ! not contain exponents that are monotonically decreasing.
  ! Before augmenting we have to determine the two lowest exponents.
  IF(IAUG .GT. 0) THEN
     IF(.NOT.CONTRACTED)THEN
        IF(BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(1)%exponents(1)&
             &.GT.BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(2)%exponents(1))THEN
           EXMIN1=BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(2)%exponents(1)
           EXMIN2=BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(1)%exponents(1)
        ELSE
           EXMIN1=BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(1)%exponents(1)
           EXMIN2=BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(2)%exponents(1)
        ENDIF
        DO IEX = 3, nprim
           IF(BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(IEX)%exponents(1).LT.EXMIN1)THEN
              EXMIN2 = EXMIN1
              EXMIN1 = BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(IEX)%exponents(1)
           ELSEIF(BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(IEX)%exponents(1).LT.EXMIN2)THEN
              EXMIN2=BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(IEX)%exponents(1)
           ENDIF
        ENDDO
        DO KAUG = 1, IAUG
           Exp = EXMIN1*EXMIN1/EXMIN2
           BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(nprim+kaug)%exponents(1)=Exp
           BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(nprim+kaug)%elms(1)=&
                &(4*Exp)**(0.5E0_realk*nang+0.25E0_realk)*PIPPI 
           BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(nprim+kaug)%UCCelms(1)=&
                &(4*Exp)**(0.5E0_realk*nang+0.25E0_realk)*PIPPI 
           EXMIN2=EXMIN1
           EXMIN1=Exp
        ENDDO
     ELSE !THE NORMALE CASE OF A SEGMENTED OR GENERAL CONTRACTED
        IF(Exponents%elms(1) .GT. Exponents%elms(2)) THEN
           EXMIN1 = Exponents%elms(2)
           EXMIN2 = Exponents%elms(1)
        ELSE
           EXMIN1 = Exponents%elms(1)
           EXMIN2 = Exponents%elms(2)
        END IF
        DO IEX = 3, nprim
           IF(Exponents%elms(IEX) .LT. EXMIN1) THEN
              EXMIN2 = EXMIN1
              EXMIN1 = Exponents%elms(IEX)
           ELSEIF(Exponents%elms(IEX) .LT. EXMIN2) THEN
              EXMIN2 = Exponents%elms(IEX)
           ENDIF
        ENDDO
        DO KAUG = 1, IAUG
           Exp = EXMIN1*EXMIN1/EXMIN2
           Exponents%elms(nprim+kaug) = Exp
           Contractionmatrix%elms(nprim+kaug+(nOrbital+kaug-1)*(nprim+iaug))=1.0E0_realk 
           EXMIN2=EXMIN1
           EXMIN1=Exp
        ENDDO
     ENDIF
     nprim = nprim + iaug
     nOrbital = nOrbital + iaug
  ENDIF

  !****************************************************************
  !*
  !* PRINT SECTION
  !*
  !***************************************************************
  IF(IPRINT .GT. 5 .AND. CONTRACTED)THEN
     WRITE(LUPRI,*)SPDFGH(nang),'-TYPE FUNCTIONS'
     WRITE(LUPRI,*)'Exponents    nprim',nprim
     call LSMAT_DENSE_PRINT(Exponents, 1, nprim, 1, 1, LUPRI)
     WRITE(LUPRI,*)'Contractionmatrix        IAUG=',IAUG,'  nOrbital',nOrbital
     CALL LSMAT_DENSE_PRINT(Contractionmatrix,1,nprim,1,nOrbital,LUPRI)
  ENDIF

  ! NOT TESTET
  IF (POLFUN) THEN
     CALL LSQUIT('polarized functions not yet implemented',lupri)
  ENDIF

END SUBROUTINE READ_COEFFICIENT_AND_EXPONENTS

subroutine determine_nNumbers_in_string(STRING,nNUMBERS)
  implicit none
  CHARACTER(len=200)    :: STRING
  integer :: nNUMBERS
  !
  logical :: INSIDENUMBER,SCIENTIFIC
  integer :: I
  nNUMBERS=0
  INSIDENUMBER=.FALSE.      
  SCIENTIFIC = .FALSE.
  DO I=1,LEN(STRING)
!    The '-' allows fixed format type numbers with no space separation
     IF((STRING(I:I).EQ.' ').OR.(STRING(I:I).EQ.'-').AND.INSIDENUMBER)THEN
!       In case of scientific number representation 1.2345678D-09 we do not 
!       accept '-' to separate two nnumbers
        INSIDENUMBER=SCIENTIFIC
     ELSEIF(STRING(I:I).EQ.' '.AND..NOT.INSIDENUMBER)THEN
        !still outside number but not yet inside new number
!     ELSEIF(STRING(I:I).EQ.'H'.AND..NOT.INSIDENUMBER)THEN
!        !still outside number but not yet inside new number
     ELSE
        IF(.NOT.INSIDENUMBER)THEN
           nNUMBERS=nNUMBERS+1
           INSIDENUMBER=.TRUE.
        ENDIF
     ENDIF
     SCIENTIFIC = STRING(I:I).EQ.'D'
  ENDDO
END subroutine DETERMINE_NNUMBERS_IN_STRING

!> \brief determine how many lines the contraction matrix is distributed over
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE LINES_OF_CONTRACTION(nOrbital,nCont,NUMBER_OF_LINES,segmentedFormat)
!*********************************************************************
!* CALCULATE ON HOW MANY LINES THE CONTRACTION COEFFICIENTS ARE 
!* WRITTEN ON
!*********************************************************************
implicit none
INTEGER     :: NUMBER_OF_LINES,nOrbital,nCont
REAL(realk) :: B,C
LOGICAL     :: segmentedFormat
!The intrisic functions DBLE makes a souble precision reak number of an integer.
IF (segmentedFormat) THEN
  NUMBER_OF_LINES = 1
  RETURN
ENDIF
IF (MOD(nOrbital,nCont).EQ.0) THEN
  NUMBER_OF_LINES = nOrbital/nCont
ELSE
  NUMBER_OF_LINES = nOrbital/nCont + 1
ENDIF
END SUBROUTINE LINES_OF_CONTRACTION

!> \brief analyze the contraction matrix
!> \author T. Kjaergaard
!> \date 2010
!>
!> ANALYSE THE CONTRACTIONCOEFFICIENTMATRIX IN ORDER TO CONSTRUCT 
!> SEGMENTED BLOCKS SO NOT TO STORE ZERO ELEMENTS AND PREPARE TO CONSTRUCT
!> AOBATCHES
!>
SUBROUTINE ANALYSE_CONTRACTIONMATRIX(LUPRI,IPRINT,BASINFO,atomtype,&
          &nAngmom,nprim,nOrbital,Contractionmatrix,ContractionmatrixNORM,Exponents)
implicit none
TYPE(BASISSETINFO)  :: BASINFO
TYPE(lsmatrix) :: Contractionmatrix,ContractionmatrixNORM
TYPE(lsmatrix) :: Exponents
INTEGER      :: atomtype,nAngmom,nprim,nOrbital,ELEMENTS
INTEGER      :: ncol,nrow,K,L,LUPRI,J,istart,iend,IPRINT,EXTRAROWS
INTEGER      :: start(maxBASISsegment),end(maxBASISsegment),I,segments,Nsegments
INTEGER      :: Nend(maxBASISsegment),Nstart(maxBASISsegment)
LOGICAL      :: NEWBLOCK,INSIDEBLOCK
INTEGER,pointer :: SEGMENTrow(:),SEGMENTcol(:)
REAL(REALK),pointer :: TEMPEXP1(:),TEMPEXP2(:)
IF (IPRINT .GT. 200)THEN
WRITE(LUPRI,*)'ANALYSE CONTRACTION MATRIX'
CALL LSMAT_DENSE_PRINT(Contractionmatrix,1,Contractionmatrix%nrow,1,&
                                &Contractionmatrix%ncol,lupri)
ENDIF
!**********************************************************
!*
!* STEP 1. DETERMINE SECTIONS OF NONZERO ELEMENTS
!*
!**********************************************************
start=0
end=0
start(1)=1
segments=0
NEWBLOCK=.FALSE.
DO I=1,nprim*nOrbital
  IF(NEWBLOCK)THEN
    IF(ABS(Contractionmatrix%elms(I)) .LT. 1.0E-30_realk) THEN 
      segments=segments+1      
      end(segments)=I-1        
      NEWBLOCK=.FALSE.
      start(segments+1)=I+1    
    ENDIF
  ELSE
    IF(ABS(Contractionmatrix%elms(I)) .LT. 1.0E-30_realk) THEN 
      start(segments+1)=start(segments+1)+1 
    ELSE
      NEWBLOCK=.TRUE.
    ENDIF
  ENDIF
  IF(I==nprim*nOrbital) then
    segments=segments+1
    end(segments)=nprim*nOrbital 
  ENDIF
ENDDO

DO I=1,segments
   IF(start(I).GT.end(I))THEN
      write(lupri,*)'Error in basis set, the last primitive most contribute to the last contracted function'
      print*,'Error in basis set, the last primitive most contribute to the last contracted function'
      call lsquit('Error in basis set, the last primitive most contribute to the last contracted function',-1)
   ENDIF
ENDDO

IF (IPRINT .GT. 200)THEN
   WRITE(lupri,*)'STEP 1'
   WRITE(lupri,*)'SEGMENT    START     END'
   DO I=1,segments
      WRITE(lupri,'(I5,2X,I5,2X,I5)') I,start(I),end(I)
   ENDDO
ENDIF
IF (segments.GT.maxBASISsegment) THEN
  WRITE(LUPRI,'(1X,A,I2,A,I2)') 'Error in ANALYSE_CONTRACTIONMATRIX. segments =', &
     & segments,' > maxBASISsegment = ',maxBASISsegment
  CALL LSQUIT('Error in ANALYSE_CONTRACTIONMATRIX. segments > maxBASISsegment',lupri)
ENDIF
!**********************************************************
!*
!* STEP 2. IF NUMBER OF NONZERO ELEMENTS IN A SECTION IS
!* LARGER THAN # EXPONENT THE NUMBER OF COLLUMS IS 
!* DETERMINED AND THE NUMBER OF BLOCKS IS DETERMINED
!*
!**********************************************************

istart=1
iend=1

DO I=1,segments
   IF (segments==1)THEN
     Nstart(1)=1
     Nend(1)=nprim*nOrbital
   ELSE
      IF(I==1)THEN
         Nstart(1)=1
         istart=istart+1
         IF(START(I+1)-END(I) > nprim)THEN  
            !THEY BELONG TO 2 DIFFERENT BLOCKS MEANING THAT I CAN 
            !END FIRST BLOCK
            Nend(1)=end(I)
            iend=2
            INSIDEBLOCK=.FALSE.
         ELSE !5-9=4
            !THEY BELONG TO SAME BLOCK
            INSIDEBLOCK=.TRUE.
            Nend(1)=end(1)+nprim  !10
            end(2)=end(1)+nprim    !10
         ENDIF
      ELSEIF(I==segments) THEN
         IF(.NOT.INSIDEBLOCK) THEN
            Nstart(istart)=START(I)  
         ENDIF
         Nend(iend)=nprim*nOrbital
      ELSE
         IF(INSIDEBLOCK)THEN
            IF(START(I+1)-END(I) > nprim)THEN  !15-10=5
               !THEY BELONG TO 2 DIFFERENT BLOCKS MEANING THAT I CAN 
               !END ONE BLOCK
               Nend(iend)=end(I)  
               iend=iend+1 
               INSIDEBLOCK=.FALSE.
            ELSE
               end(I+1)=end(I)+nprim  
               !THEY BELONG TO SAME BLOCK
            ENDIF
         ELSE
            !I START A NEW BLOCK
            IF(START(I+1)-END(I) > nprim)THEN
               !THEY BELONG TO 2 DIFFERENT BLOCKS MEANING THAT I CAN 
               !END FIRST BLOCK
               Nstart(istart)=start(I)
               Nend(iend)=end(I)
               iend=iend+1 
               istart=istart+1
               INSIDEBLOCK=.FALSE.
            ELSE
               !THEY BELONG TO SAME BLOCK
               Nstart(istart)=start(I)
               istart=istart+1
               INSIDEBLOCK=.TRUE.
            ENDIF
         ENDIF
      ENDIF
   ENDIF
ENDDO

Nsegments=iend

IF (IPRINT .GT. 200)THEN
   WRITE(lupri,*)'STEP 2'
   WRITE(lupri,*)'SEGMENT    START     END'
   DO I=1,Nsegments
      WRITE(lupri,'(I5,2X,I5,2X,I5)') I,Nstart(I),Nend(I)
   ENDDO
ENDIF

CALL INIT_BASISSETINFO_ContractionM(BasInfo,atomtype,&
     &nAngmom,Nsegments)

!**********************************************************
!*
!* STEP 3. DETERMINE nROW AND nCOL FOR EACH SEGMENT
!*
!**********************************************************
call mem_alloc(SEGMENTROW,Nsegments)
call mem_alloc(SEGMENTCOL,Nsegments)

DO I=1,Nsegments
  ELEMENTS=Nend(I)-Nstart(I)+1
  IF(ELEMENTS < nprim)THEN
    SEGMENTROW(I)=ELEMENTS
    SEGMENTcol(I)=1
  ELSE
    SEGMENTcol(I)=ELEMENTS/nprim
    IF(MOD(ELEMENTS,nprim).NE. 0)THEN
       SEGMENTcol(I)=SEGMENTcol(I)+1
       SEGMENTrow(I)=ELEMENTS-(SEGMENTcol(I)-1)*nprim
       DO K=1,SEGMENTrow(I)
          IF(Nstart(I)+K+(SEGMENTcol(I)+1-1)*nprim-1.LE.nprim*norbital)THEN
             IF(ABS(Contractionmatrix%elms(Nstart(I)+K+(SEGMENTcol(I)+1-1)*nprim-1)) .GT. 1.0E-30_realk)THEN
                WRITE(LUPRI,*)'CC(',Nstart(I)+K+(SEGMENTcol(I)+1-1)*nprim-1,')=',&
                     &Contractionmatrix%elms(Nstart(I)+K+(SEGMENTcol(I)+1-1)*nprim-1)
                CALL LSQUIT('something is wrong in ANALYSE_CONTRACTIONMATRIX',lupri)
             ENDIF
          ENDIF
       ENDDO
       EXTRAROWS=0
       DO L=1,SEGMENTcol(I)
          IF(Nstart(I)+SEGMENTrow(I)+(L-1)*nprim.LE.nprim*norbital)THEN
!          WRITE(LUPRI,*)'CC(',Nstart(I)+SEGMENTrow(I)+(L-1)*nprim,')=',Contractionmatrix%elms(Nstart(I)+SEGMENTrow(I)+(L-1)*nprim)
             IF(ABS(Contractionmatrix%elms(Nstart(I)+SEGMENTrow(I)+(L-1)*nprim)) .GT. 1.0E-30_realk)EXTRAROWS=EXTRAROWS+1          
          ENDIF
       ENDDO
       IF(EXTRAROWS .GT. 0) SEGMENTrow(I)=SEGMENTrow(I)+1
    ELSE
       SEGMENTrow(I)=ELEMENTS-(SEGMENTcol(I)-1)*nprim
    ENDIF
  ENDIF
ENDDO

IF (IPRINT .GT. 200)THEN
   WRITE(lupri,*)'STEP 3'
   WRITE(lupri,*)'SEGMENT    ROW     COL'
   DO I=1,Nsegments
      WRITE(lupri,'(I5,2X,I5,2X,I5)') I,SEGMENTrow(I),SEGMENTcol(I)
   ENDDO
ENDIF

!**********************************************************
!*
!* STEP 4. REORDER PRIMITIVES
!*
!**********************************************************
!DO I=1,Nsegments
!   nrow=SEGMENTrow(I)
!   ncol=SEGMENTcol(I)
!   WRITE(lupri,*)'SECTION =',I,'nrow',nrow,'ncol',ncol
!ENDDO

K=1
I=1
DO 
   IF(I.EQ.Nsegments)EXIT
   ELEMENTS=Nend(I)-Nstart(I)+1
   nrow=SEGMENTrow(I)
   ncol=SEGMENTcol(I)
   
   IF(Exponents%elms(K) .LE. Exponents%elms(nrow+K))THEN
!      WRITE(lupri,*)'EXP',Exponents%elms(K),'is less than',Exponents%elms(nrow+K)
!      WRITE(lupri,*)'SO WE SWAP'
!      WRITE(lupri,*)'THE ORIGINAL VERSION'
!      DO J=K,K+SEGMENTrow(I)+SEGMENTrow(I+1)-1
!         WRITE(lupri,*)'EXP(',J,')=',Exponents%elms(J)
!      ENDDO
      !SWAP EXPONENTS AND SWAP THE SEGMENTS
      call mem_alloc(TEMPEXP1,SEGMENTrow(I))
      DO J=1,SEGMENTrow(I)
         TEMPEXP1(J)=Exponents%elms(K+J-1)
      ENDDO
!      DO J=1,SEGMENTrow(I)
!         WRITE(lupri,*)'TEMPEXP1(',J,')=',TEMPEXP1(J)
!      ENDDO
      call mem_alloc(TEMPEXP2,(SEGMENTrow(I+1)))
      DO J=1,SEGMENTrow(I+1)
         TEMPEXP2(J)=Exponents%elms(K+SEGMENTrow(I)+J-1)
      ENDDO
!      DO J=1,SEGMENTrow(I+1)
!         WRITE(lupri,*)'TEMPEXP2(',J,')=',TEMPEXP2(J)
!      ENDDO

      DO J=1,SEGMENTrow(I+1)
         Exponents%elms(K+J-1)=TEMPEXP2(J)
      ENDDO
      DO J=1,SEGMENTrow(I)
         Exponents%elms(K+SEGMENTrow(I+1)+J-1)=TEMPEXP1(J)
      ENDDO

      !SWAP SEGMENTS
      nrow=SEGMENTrow(I+1)
      SEGMENTrow(I+1)=SEGMENTrow(I)
      SEGMENTrow(I)=nrow
      ncol=SEGMENTcol(I+1)
      SEGMENTcol(I+1)=SEGMENTcol(I)
      SEGMENTcol(I)=ncol     
      nrow=Nend(I+1)
      Nend(I+1)=Nend(I)
      Nend(I)=nrow
      nrow=Nstart(I+1)
      Nstart(I+1)=Nstart(I)
      Nstart(I)=nrow
      call mem_dealloc(TEMPEXP1)
      call mem_dealloc(TEMPEXP2)
!      WRITE(lupri,*)'THE SWAPED VERSION'
!      DO J=K,K+SEGMENTrow(I+1)+SEGMENTrow(I)-1
!         WRITE(lupri,*)'EXP(',J,')=',Exponents%elms(J)
!      ENDDO
      K=K+SEGMENTrow(I+1)
      I=I+1
   ELSE     
      K=K+nrow
      I=I+1
   ENDIF
ENDDO

!**********************************************************
!*
!* STEP 5. BUILD CONTRACTIONMATRIX
!*
!**********************************************************

DO I=1,Nsegments
   nrow=SEGMENTrow(I)
   ncol=SEGMENTcol(I)
   CALL INIT_BASISSETINFO_elms(BasInfo,atomtype,nAngmom,I,nrow,ncol)
   DO K=1,nrow
      DO L=1,ncol
         BASINFO%ATOMTYPE(atomtype)%SHELL(nAngmom)%&
              &segment(I)%elms(K+(L-1)*nrow)=&
              &ContractionmatrixNORM%elms(Nstart(I)+K+(L-1)*nprim-1)
      ENDDO
   ENDDO
   DO K=1,nrow
      DO L=1,ncol
         BASINFO%ATOMTYPE(atomtype)%SHELL(nAngmom)%&
              &segment(I)%UCCelms(K+(L-1)*nrow)=&
              &Contractionmatrix%elms(Nstart(I)+K+(L-1)*nprim-1)
      ENDDO
   ENDDO
   IF(IPRINT .GT. 200)THEN
      WRITE(LUPRI,*)'Exponents nr.',I,'nrow',nrow
      call OUTPUT(BASINFO%ATOMTYPE(atomtype)%SHELL(nAngmom)%&
           &segment(I)%elms, 1, nrow, 1, ncol, nrow, ncol, 1, LUPRI)
   ENDIF
   IF(IPRINT .GT. 200)THEN
      WRITE(LUPRI,*)'UCC Exponents nr.',I,'nrow',nrow
      call OUTPUT(BASINFO%ATOMTYPE(atomtype)%SHELL(nAngmom)%&
           &segment(I)%UCCelms, 1, nrow, 1, ncol, nrow, ncol, 1, LUPRI)
   ENDIF
ENDDO

!**********************************************************
!*
!* STEP 6. BUILD EXPONENTS
!*
!**********************************************************

k=1
DO I=1,Nsegments
  nrow=BASINFO%ATOMTYPE(atomtype)%SHELL(nAngmom)%&
        &segment(I)%nrow
  DO J=1,nrow
     BASINFO%ATOMTYPE(atomtype)%SHELL(nAngmom)%&
     &segment(I)%Exponents(J)=Exponents%elms(k)
     k=k+1
  ENDDO 
  IF(IPRINT .GT. 200)THEN
     WRITE(LUPRI,*)'Exponents nr.',I,'nrow',nrow
     call OUTPUT(BASINFO%ATOMTYPE(atomtype)%SHELL(nAngmom)%&
          &segment(I)%Exponents, 1, nrow, 1, 1, nrow, 1, 1, LUPRI)
  ENDIF
ENDDO

call mem_dealloc(SEGMENTROW)
call mem_dealloc(SEGMENTCOL)

END SUBROUTINE ANALYSE_CONTRACTIONMATRIX

!> \brief Normalize the orbitals
!> \author T. Kjaergaard
!> \date 2010
SUBROUTINE NORMALIZE_ORBITALS(LUPRI,IPRINT,angmom,nOrbitals,nExp,Exp,CmatrixIN,CmatrixOUT)
implicit none
real(realk) :: PI,PIPPI,D0,D1,D2,D4,DP25,DP5,DP75,THRMIN,EXPMIN
TYPE(lsmatrix):: Exp       !Vector of exponents
TYPE(lsmatrix):: CmatrixIN   !matrix of non-Normalized ContractionCoefficients
TYPE(lsmatrix):: CmatrixOUT  !OUTPUT: Normalised matrix of C Coefficients
INTEGER     :: nEXP      !number of exponents
INTEGER     :: nOrbitals !number of contracted orbitals  
INTEGER     :: angmom    !angular momentum=angmom-1
INTEGER     :: LUPRI,M,L,N,IPRINT
real(realk) :: SUM,T
IF(IPRINT .GT. 10) THEN
   WRITE(LUPRI,*)'THE UNNORMALISED Cmatrix'
   CALL LSMAT_DENSE_PRINT(CmatrixIN,1,CmatrixIN%nrow,1,CmatrixIN%ncol,lupri)
ENDIF
PI=3.14159265358979323846E0_realk
D0 = 0.0E0_realk
D1 = 1.0E0_realk
D2 = 2.0E0_realk
D4 = 4.0E0_realk
DP25 = 0.25E0_realk
DP5 = 0.5E0_realk
DP75 = 0.75E0_realk
THRMIN = 1E-17_realk
EXPMIN = 1E-10_realk
PIPPI = (0.5E0_realk/PI)**(0.75E0_realk) !=\sqrt{\frac{1}{(2\pi)^{3/2}}} !3 dimensions

IF (nExp.EQ. 1 .AND. ABS(Exp%elms(1)).LE.EXPMIN) THEN
  CmatrixOUT%elms(1) = D1 
ELSE 
  DO N = 1, nOrbitals
    SUM = D0
    DO L = 1, nExp
      DO M = 1, nExp
        T = D2*SQRT(Exp%elms(L)*Exp%elms(M))/(Exp%elms(L)+Exp%elms(M))
        SUM = SUM + CmatrixIN%elms(L+(N-1)*nExp)*CmatrixIN%elms(M+(N-1)*nExp)&
                                              &*(T**(angmom + DP5))
      END DO
    END DO
    IF (SQRT(SUM) .LT. THRMIN) THEN
      WRITE (LUPRI,*)&
      & 'INPUT ERROR: CGTO no.',N,'of angular momentum',angmom,&
      & 'has zero norm.'
      CALL LSQUIT('CGTO with zero norm.',lupri)
    ENDIF
    SUM=D1/SQRT(SUM)
    DO L=1,nExp
      CmatrixOUT%elms(L+(N-1)*nExp)=CmatrixIN%elms(L+(N-1)*nExp)*SUM&
                                 &*(D4*Exp%elms(L))**(DP5*angmom+DP25)*PIPPI
    END DO
  END DO
END IF
IF(IPRINT .GT. 10) THEN
   WRITE(LUPRI,*)'THE NORMALISED Cmatrix'
   CALL LSMAT_DENSE_PRINT(CmatrixOUT,1,CmatrixOUT%nrow,1,CmatrixOUT%ncol,lupri)
ENDIF
END SUBROUTINE NORMALIZE_ORBITALS

END MODULE BUILDBASISSET



