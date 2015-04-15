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
!> label 'REGULAR  ' or 'AUXILIARY' or 'CABS     ' or 'JKAUX    ' or 'ADMM     '
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
BASINFO%GeminalScalingFactor = 1.0E0_realk
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

 !Setting the F12 geminal scaling factor 
 !see J. Chem. Phys 128, 084102
 IPOS = INDEX(BASISSETNAME,'aug-cc-pVDZ')
 IF (IPOS .NE. 0) BASINFO%GeminalScalingFactor = 1.1E0_realk
 IPOS = INDEX(BASISSETNAME,'aug-cc-pVTZ')
 IF (IPOS .NE. 0) BASINFO%GeminalScalingFactor = 1.2E0_realk
 IPOS = INDEX(BASISSETNAME,'aug-cc-pVQZ')
 IF (IPOS .NE. 0) BASINFO%GeminalScalingFactor = 1.4E0_realk
 IPOS = INDEX(BASISSETNAME,'aug-cc-pV5Z')
 IF (IPOS .NE. 0) BASINFO%GeminalScalingFactor = 1.4E0_realk
 IPOS = INDEX(BASISSETNAME,'cc-pVDZ-F12')
 IF (IPOS .NE. 0) BASINFO%GeminalScalingFactor = 0.9E0_realk
 IPOS = INDEX(BASISSETNAME,'cc-pVTZ-F12')
 IF (IPOS .NE. 0) BASINFO%GeminalScalingFactor = 1.0E0_realk
 IPOS = INDEX(BASISSETNAME,'cc-pVQZ-F12')
 IF (IPOS .NE. 0) BASINFO%GeminalScalingFactor = 1.1E0_realk
ELSE
! call mem_alloc(BINDEXES,BASISSETLIBRARY(iBas)%nbasissets)
 BASINFO%DunningsBasis = BASISSETLIBRARY(iBas)%DunningsBasis
 BASINFO%GeminalScalingFactor = BASISSETLIBRARY(iBas)%GeminalScalingFactor
 SELECT CASE(BASISLABEL)
 CASE('REGULAR  ') !identical to BasParamLABEL(RegBasParam)
    BASINFO%Labelindex = RegBasParam
 CASE('AUXILIARY') !identical to BasParamLABEL(AUXBasParam)
    BASINFO%Labelindex = AUXBasParam
 CASE('CABS     ') !identical to BasParamLABEL(CABBasParam)
    BASINFO%Labelindex = CABBasParam
 CASE('JKAUX    ') !identical to BasParamLABEL(JKBasParam)
    BASINFO%Labelindex = JKBasParam
 CASE('ADMM     ') !identical to BasParamLABEL(ADMBasParam)
    BASINFO%Labelindex = ADMBasParam
 CASE('VALENCE  ') !identical to BasParamLABEL(VALBasParam)
    !build as subset of regular in trilevel algorithm
    call LSQUIT('VALENCE no legal keyword in Build_Basis.',lupri)  
 CASE('GCTRANS  ')!identical to BasParamLABEL(GCTBasParam)
    !build as a transformation basis in trilevel algorithm
    call LSQUIT('GCTRANS no legal keyword in Build_Basis.',lupri)  
 CASE('CABSP    ') !identical to BasParamLABEL(CAPBasParam)
    BASINFO%Labelindex = CAPBasParam
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
  IF(IAUG.EQ.0)THEN
     DO IBAS2=1,NBASDIR
        STRING(IBAS2)%p = BASISDIR(IBAS2)%p(1:LEN_BASISDIR(IBAS2))//BASISSETLIBRARY(iBas)%BASISSETNAME(I)(1:NEND)//' '
     ENDDO
  ELSE
     DO IBAS2=1,NBASDIR
        STRING(IBAS2)%p = BASISDIR(IBAS2)%p(1:LEN_BASISDIR(IBAS2))//BASISSETLIBRARY(iBas)%BASISSETNAME(I)(2:NEND)//' '
     ENDDO
  ENDIF
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
!> labelindex=RegBasParam is for regular MoleculeSpecific ordering 
!> labelindex=AuxBasParam is for auxiliary MoleculeSpecific ordering
!> labelindex=CABBasParam is for cabs MoleculeSpecific ordering
!> labelindex=JKBasParam is for JK MoleculeSpecific ordering
!> labelindex=ADMBasParam is for ADMM MoleculeSpecific ordering
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
!  labelindex=RegBasParam is for regular MoleculeSpecific ordering 
!  labelindex=AuxBasParam is for auxiliary MoleculeSpecific ordering
!  labelindex=CABBasParam is for cabs MoleculeSpecific ordering
!  labelindex=JKBasParam is for JKaux MoleculeSpecific ordering
!  labelindex=ADMBasParam is for ADMM aux MoleculeSpecific ordering
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

SUBROUTINE COPY_MOLECULE_IDTYPE(MOLECULE,iR,iC)
implicit none
!> contains all info about the molecule (atoms, charge,...)
TYPE(MOLECULEINFO)        :: MOLECULE
!> labelindexes
INTEGER,intent(in)        :: iR,iC
!
INTEGER                   :: I

DO I=1,MOLECULE%natoms
   MOLECULE%ATOM(I)%IDtype(iC) = MOLECULE%ATOM(I)%IDtype(iR)
ENDDO

END SUBROUTINE COPY_MOLECULE_IDTYPE

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
!  BASNAM(1:79) = BASTMP(2:80)
!  BASNAM(80:80) = ' '
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
   CALL LSQUIT('ecp not supported',-1)
   BASTMP = BASNAM
   !  BASNAM(1:75) = BASTMP(5:80)
   !  BASNAM(76:80) = '     '
   IAUG = 1
ELSEIF (BASNAM(2:8) .EQ. 'aug-ecp') THEN
   CALL LSQUIT('ecp not supported',-1)
   BASTMP = BASNAM
   !  BASNAM(1:74) = BASTMP(6:80)
   !  BASNAM(75:80) = '      '
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
INTEGER             :: nprim, nOrbitals,K,nAngmom,JJ
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
  DO JJ=1,200
     STRING(JJ:JJ) = ' '
  ENDDO
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
        READ (STRING, '(A1, I8)') SIGN, NUCLEARCHARGE
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
  LOGICAL               :: POLFUN,CONTRACTED,segmentedFormat,BLANK,FOUNDnewContractionCoeffLine
  INTEGER               :: atype,nang,nprim,nOrbital,IAUG,NUMNUMOLD
  INTEGER               :: J,NUMBER_OF_LINES,KNTORB,NUMNUM,KAUG,nNumbers,LIST(2,20)
  CHARACTER(len=280)    :: STRING
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
     READ(LUBAS, '(A280)', IOSTAT = IOS) STRING
     IF(ios /= 0)THEN
        WRITE (LUPRI,'(A)') ' Error in basisset file'
        WRITE(lupri,'(A)')'This could mean that the line containing exponents'
        WRITE(lupri,'(A)')'and contraction coefficients fill more than 280 characters'
        WRITE(lupri,'(A)')'Which means you need to manually split the line'
        CALL LSQUIT('Error in basisset file',lupri)
     ELSE
        READ (STRING, '(A1)') SIGN
        IF (SIGN .EQ. ' ') THEN
           CALL TEST_BLANK_LINE(BLANK, STRING)
           IF (.NOT. BLANK) THEN
              !We have found a line with a primitive and some coeffecients
              J=J+1 !The primitive index 
              NUMBER_OF_LINES=1           
              KNTORB = 0 !number of Contraction Coefficients read 
              DO WHILE(KNTORB.NE.nOrbital)
                 !LINES_OF_CONTRACTION returns the number of lines 
                 !the contraction coeffecients written on.
                 CALL determine_nNumbers_in_string(STRING,nNUMBERS,LIST)
                 IF(nNUMBERS.GT.7)THEN
                    IF(NUMBER_OF_LINES.EQ.1)THEN
                       WRITE(lupri,*)'WARNING This basis format violate the old dalton basis format'
                       WRITE(lupri,*)'Which consist of 1 Exponent F16.9 and up to 6 contraction'
                       WRITE(lupri,*)'coefficients on the first line followed by a up to 7 contraction coefficients'
                       WRITE(lupri,*)'on the following lines until the full number of contraction coefficients are given'
                       WRITE(lupri,*)'We will try to this basis set, but this code is not very well testet. TK'
                    ELSE
                       WRITE(lupri,*)'WARNING This is a continuation line, which violte the old dalton basis format'
                       WRITE(lupri,*)'The continuation line should contain at most 7 contraction coefficients'
                       WRITE(lupri,*)'We will try to this basis set, but this code is not very well testet. TK'           
                    ENDIF
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
                    IF(NUMBER_OF_LINES.EQ.1)THEN
                       READ (STRING, '(F16.9)') Exp
                       BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(J)%Exponents(1) = Exp
                       BASINFO%ATOMTYPE(atype)%SHELL(nang)%segment(J)%elms(1)= &
                            &(4*Exp)**(0.5E0_realk*nang+0.25E0_realk)*PIPPI 
                       KNTORB = nNUMBERS-1
                    ELSE
                       !continuation line with no exponents
                       KNTORB = KNTORB + nNUMBERS
                    ENDIF
                 ELSE
                    IF(NUMBER_OF_LINES.EQ.1)THEN                    
                       READ (STRING(LIST(1,1):LIST(2,1)),*) Exponents%elms(J)
                       DO I = 1, nNUMBERS-1
                          READ (STRING(LIST(1,I+1):LIST(2,I+1)),*) ContractionMatrix%elms(J+(I-1)*(nprim+IAUG))
                       ENDDO
                       !         Reading the first line with exponents and contractioncoeffecients
                       !     READ (STRING, '(F16.9, 6F12.9)') Exponents%elms(J),&
                       !          &(ContractionMatrix%elms(J+(I-1)*(nprim+IAUG)),I = 1, KNTORB)
                       !         If there are more lines with contraction-coeffecients
                       !         they will be read here.
                       KNTORB = nNUMBERS-1
                    ELSE
                       !We now have a line with contraction-coeffecients.
                       !         READ (STRING,'(F16.9,6F12.9)')&
                       !              &(Contractionmatrix%elms(J+(M-1)*(nprim+IAUG)),&
                       !              & M = 6 + (I-2)*7 +1, KNTORB)
                       DO I = 1, nNUMBERS
                          READ (STRING(LIST(1,I):LIST(2,I)),*) ContractionMatrix%elms(J+(I+KNTORB-1)*(nprim+IAUG))
                       ENDDO
                       KNTORB = KNTORB + nNUMBERS
                    ENDIF
                 ENDIF
                 IF(KNTORB.NE.nOrbital)THEN
                    IF (segmentedFormat) exit
                    FOUNDnewContractionCoeffLine = .FALSE.
                    DO WHILE(.NOT.FOUNDnewContractionCoeffLine)
                       READ(LUBAS, '(A280)', IOSTAT = IOS) STRING
                       IF(ios /= 0)THEN
                          WRITE (LUPRI,'(2A)') ' Error in basisset file'
                          CALL LSQUIT('Error in basisset file',lupri)
                       ELSE
                          READ (STRING, '(A1)') SIGN
                          IF (SIGN .EQ. ' ') THEN
                             CALL TEST_BLANK_LINE(BLANK, STRING)
                             IF (.NOT. BLANK) THEN
                                NUMBER_OF_LINES = NUMBER_OF_LINES + 1
                                FOUNDnewContractionCoeffLine = .TRUE.
                             END IF
                          ELSE
                             cycle
                          ENDIF
                       ENDIF
                    ENDDO
                 ENDIF
                 IF(KNTORB.GT.nOrbital)CALL LSQUIT('Error in BuildBasis: too many coefficients',-1)
              ENDDO
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

subroutine determine_nNumbers_in_string(STRING,nNUMBERS,LIST)
  implicit none
  CHARACTER(len=280)    :: STRING
  integer :: nNUMBERS,LIST(2,20)
  !
  logical :: INSIDENUMBER,SCIENTIFIC
  integer :: I
  nNUMBERS=0
  INSIDENUMBER=.FALSE.      
  SCIENTIFIC = .FALSE.
  DO I=1,LEN(STRING)
!    The '-' allows fixed format type numbers with no space separation
     IF(INSIDENUMBER.AND.STRING(I:I).EQ.' ')THEN
        INSIDENUMBER=.FALSE.
        LIST(2,nNUMBERS) = I-1 !Last charater in String
     ELSEIF((STRING(I:I).EQ.'-'.OR.STRING(I:I).EQ.'+').AND.INSIDENUMBER)THEN
!       In case of scientific number representation 1.2345678D-09 we do not 
!       accept '-' to separate two nnumbers
        IF(SCIENTIFIC)THEN
           !still inside a number 
           INSIDENUMBER=.TRUE.
        ELSE
           LIST(2,nNUMBERS) = I-1 !Last charater in String
           nNUMBERS=nNUMBERS+1    !new number starting with - or plus
           LIST(1,nNUMBERS) = I   !First charater in String
           INSIDENUMBER=.TRUE.
        ENDIF
     ELSEIF(STRING(I:I).EQ.' '.AND..NOT.INSIDENUMBER)THEN
        !still outside number but not yet inside new number
!     ELSEIF(STRING(I:I).EQ.'H'.AND..NOT.INSIDENUMBER)THEN
!        !still outside number but not yet inside new number
     ELSE
        IF(.NOT.INSIDENUMBER)THEN
           nNUMBERS=nNUMBERS+1
           LIST(1,nNUMBERS) = I !First charater in String
           INSIDENUMBER=.TRUE.
        ENDIF
     ENDIF
     SCIENTIFIC = STRING(I:I).EQ.'D' .OR. STRING(I:I).EQ.'E'
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
SUBROUTINE ANALYSE_CONTRACTIONMATRIX(LUPRI,IPRINT,BASINFO,at,&
          & nAngmom,nprim,nOrbital,Contractionmatrix,&
          & ContractionmatrixNORM,Exponents)
implicit none
TYPE(BASISSETINFO),intent(inout)  :: BASINFO
TYPE(lsmatrix),intent(in) :: Contractionmatrix,ContractionmatrixNORM
TYPE(lsmatrix),intent(inout) :: Exponents
INTEGER,intent(in)      :: at,nAngmom,nprim,nOrbital
!local variables
INTEGER      :: ncol,nrow,K,L,LUPRI,J,IPRINT
INTEGER      :: I,Nsegments,kk,icol
INTEGER      :: SEGMENTrow(nOrbital),SEGMENTcol(nOrbital)
INTEGER      :: nOrb(nprim),mPrim(nOrbital)
LOGICAL      :: Segmented,Prim(nPrim,nOrbital),PerfomMerge
LOGICAL      :: merged(nOrbital),MergedPrim(nPrim,nOrbital)
LOGICAL      :: SEGMENTCOLID(nOrbital,nOrbital)
real(realk),pointer  :: CCN(:,:),CC(:,:) 
real(realk)  :: newExponents,newCC,newCCN

IF (IPRINT .GT. 200)THEN
   WRITE(LUPRI,*)'ANALYSE CONTRACTION MATRIX'
   CALL LSMAT_DENSE_PRINT(ContractionmatrixNorm,1,ContractionmatrixNorm%nrow,1,&
        &ContractionmatrixNorm%ncol,lupri)
ENDIF

nrow = Contractionmatrix%nrow
ncol = Contractionmatrix%ncol
call mem_alloc(CCN,nrow,ncol)
call dcopy (nrow*ncol,ContractionmatrixNorm%elms,1,CCN,1)
call mem_alloc(CC,nrow,ncol)
call dcopy (nrow*ncol,Contractionmatrix%elms,1,CC,1)
IF (IPRINT .GT. 200)THEN
   WRITE(LUPRI,*)'nrow,ncol,nOrbital,nprim',nrow,ncol,nOrbital,nprim
   WRITE(LUPRI,*)'CCN'
   call ls_output(CCN,1,nrow,1,ncol,nrow,ncol,1,LUPRI)
   WRITE(LUPRI,*)'CC'
   call ls_output(CC,1,nrow,1,ncol,nrow,ncol,1,LUPRI)
ENDIF
!reorder primitives at this level using basic bubblesort
DO K=1,nrow
   DO J=1,nrow-1
      IF(Exponents%elms(J+1).GT. Exponents%elms(J))THEN
         newExponents=Exponents%elms(J)
         Exponents%elms(J)=Exponents%elms(J+1)
         Exponents%elms(J+1)=newExponents
         do I=1,nOrbital
            newCC=CC(J,I)
            CC(J,I)=CC(J+1,I)
            CC(J+1,I)=newCC
            newCCN=CCN(J,I)
            CCN(J,I)=CCN(J+1,I)
            CCN(J+1,I)=newCCN
         enddo
      ENDIF
   ENDDO
ENDDO

!Determine the number of orbitals each primitive function contributes to
do i=1,nPrim
   nOrb(i) = 0 
   do j=1,nOrbital
      IF(ABS(CC(i,j)).GT.1.0E-12_realk) nOrb(i) = nOrb(i) + 1
   enddo
enddo

Segmented = .TRUE.
do i=1,nPrim
   IF(nOrb(i).NE.1)THEN
      !at least one primitive contribute to more than 1 contracted function
      !hence the basis set is general contracted
      Segmented = .FALSE.
   ENDIF
enddo

IF(Segmented)THEN
 IF (IPRINT .GT. 2)WRITE(LUPRI,*)'This BasisBlock is a Segmented Basis  '
 !all primitive orbitals only contribute to one contracted function
 nSegments = ncol
 IF (IPRINT .GT. 5)WRITE(LUPRI,*)'nSegments',nSegments
 CALL INIT_BASISSETINFO_ContractionM(BasInfo,at,nAngmom,Nsegments)

 !DETERMINE number of Primitives FOR EACH SEGMENT
 do j=1,nOrbital
    mPrim(j) = 0 
    do i=1,nPrim
       IF(ABS(CC(i,j)).GT.1.0E-12_realk) mPrim(j) = mPrim(j) + 1
    enddo
 enddo

 DO J=1,Nsegments
   !For each segment copy info into structure
   nrow=mPrim(J)
   ncol=1
   CALL INIT_BASISSETINFO_elms(BasInfo,at,nAngmom,J,nrow,ncol)
   K=0
   do i=1,nPrim
      IF(ABS(CC(i,j)).GT.1.0E-12_realk)THEN
         K = K + 1
         BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%elms(K)=CCN(I,J)
         BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%UCCelms(K)=CC(I,J)
         BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%Exponents(K)=Exponents%elms(I)
      ENDIF
   enddo

   !Print Stuff
   IF(IPRINT .GT. 200)THEN
    WRITE(LUPRI,*)'Exponents nr.',J,'nrow',nrow
    call LS_OUTPUT(BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%Exponents,1,nrow,1,1,nrow,1,1,LUPRI)
   ENDIF
   IF(IPRINT .GT. 200)THEN
    WRITE(LUPRI,*)'Coefficients nr.',J,'nrow',nrow
    call LS_OUTPUT(BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%elms,1,nrow,1,ncol,nrow,ncol,1,LUPRI)
   ENDIF
   IF(IPRINT .GT. 200)THEN
    WRITE(LUPRI,*)'Unnormalizes Coefficients nr.',J,'nrow',nrow
    call LS_OUTPUT(BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%UCCelms,1,nrow,1,ncol,nrow,ncol,1,LUPRI)
   ENDIF
 ENDDO
ELSE
   IF(IPRINT .GT. 2)WRITE(LUPRI,*)'This BasisBlock is a General Contracted Basis'
   !primitive orbitals can contribute to more than one contracted function
   
   !for each contracted function determine which Primitives contribute. 
   do j=1,nOrbital
      do i=1,nPrim
         Prim(i,j) = ABS(CC(i,j)).GT.1.0E-12_realk
      enddo
   enddo
   !loop over contracted functions if j can be merged with j+1 do so
   !Determine the number of segments and sizes of the segments
   merged = .FALSE.
   Nsegments = 0
   SEGMENTCOLID = .FALSE.
   do j=1,nOrbital
      IF(.NOT.merged(j))THEN
         Nsegments = Nsegments + 1  
         do i=1,nPrim
            MergedPrim(i,Nsegments) = Prim(i,j) 
         enddo
         SEGMENTROW(Nsegments) = COUNT(MergedPrim(:,Nsegments))
         SEGMENTCOL(Nsegments) = 1
         SEGMENTCOLID(j,Nsegments) = .TRUE.
         do k=j+1,nOrbital
            IF(.NOT.merged(k))THEN
               PerfomMerge = TrueSubset(MergedPrim(:,Nsegments),Prim(:,k),nrow)
               IF(PerfomMerge)THEN
                  !merge k into j 
                  merged(k) = .TRUE.
                  do i=1,nPrim
                     MergedPrim(i,Nsegments) = Prim(i,k) .OR. MergedPrim(i,Nsegments)
                  enddo
                  SEGMENTROW(Nsegments) = COUNT(MergedPrim(:,Nsegments))
                  SEGMENTCOL(Nsegments) = SEGMENTCOL(Nsegments) + 1  
                  SEGMENTCOLID(k,Nsegments) = .TRUE.
               ENDIF
            ENDIF
         enddo
      ENDIF
   enddo
   IF(IPRINT.GT.10)THEN
      WRITE(lupri,*)'Nsegments',Nsegments
      WRITE(lupri,*)'SEGMENTROW:',SEGMENTROW(1:Nsegments)
      WRITE(lupri,*)'SEGMENTCOL:',SEGMENTCOL(1:Nsegments)
      do k=1,nOrbital
         WRITE(lupri,*)'SEGMENTCOLID(',k,'):',SEGMENTCOLID(k,1:Nsegments)
      enddo
   ENDIF
   CALL INIT_BASISSETINFO_ContractionM(BasInfo,at,nAngmom,Nsegments)
   DO J=1,Nsegments
      nrow=SEGMENTROW(J)
      ncol=SEGMENTCOL(J)
      !For each segment copy info into structure
      CALL INIT_BASISSETINFO_elms(BasInfo,at,nAngmom,J,nrow,ncol)
      kk=0
      do i=1,nPrim
         IF(MergedPrim(i,J))THEN
            kk=kk+1           
            BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%Exponents(kk)=Exponents%elms(I)
         ENDIF
      enddo
      icol = 0
      do k=1,nOrbital
         IF(SEGMENTCOLID(k,J))THEN
            icol = icol + 1
            kk=0
            do i=1,nPrim
               IF(MergedPrim(i,J))THEN
                  kk=kk+1
                  BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%elms(kk+(icol-1)*nrow)=CCN(I,k)
                  BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%UCCelms(kk+(icol-1)*nrow)=CC(I,k)
               ENDIF
            enddo
         ENDIF
      enddo
      !Print Stuff
      IF(IPRINT .GT. 200)THEN
         WRITE(LUPRI,*)'Exponents nr.',J,'nrow',nrow
         call LS_OUTPUT(BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%Exponents,1,nrow,1,1,nrow,1,1,LUPRI)
      ENDIF
      IF(IPRINT .GT. 200)THEN
         WRITE(LUPRI,*)'Coefficients nr.',J,'nrow',nrow,'ncol',ncol
         call LS_OUTPUT(BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%elms,1,nrow,1,ncol,nrow,ncol,1,LUPRI)
      ENDIF
      IF(IPRINT .GT. 200)THEN
         WRITE(LUPRI,*)'Unnormalizes Coefficients nr.',J,'nrow',nrow,'ncol',ncol
         call LS_OUTPUT(BASINFO%ATOMTYPE(at)%SHELL(nAngmom)%segment(J)%UCCelms,1,nrow,1,ncol,nrow,ncol,1,LUPRI)
      ENDIF
   enddo
ENDIF

call mem_dealloc(CCN)
call mem_dealloc(CC)

END SUBROUTINE ANALYSE_CONTRACTIONMATRIX

logical function TrueSubset(PrimJ,PrimK,ndim)
implicit none
integer :: ndim
logical :: PrimJ(ndim),PrimK(ndim)
!
integer :: I 

TrueSubset = .FALSE.
do I = 1,ndim
   IF(PrimJ(I).AND.PrimK(I))THEN
      !share primitive orbital
      TrueSubset = .TRUE.
      EXIT
   ENDIF
enddo
end function TrueSubset

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



