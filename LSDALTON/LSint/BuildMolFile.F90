!> @file
!> Write a molecule input/output file
MODULE WRITEMOLEFILE
use files
  use precision
  use TYPEDEF
  use molecule_module
  use fundamental
  use lattice_type
contains
!> \brief builds a molecule input or output file
!> \author T. Kjaergaard
!> \date 2010
!> \param filename name of the molecule file to be written
!> \param molecule the molecule structure to be built
!> \param basis the basisinfo structure - info on basis sets 
!> \param lupri the logical unit number for the output file
SUBROUTINE WRITE_MOLECULE_OUTPUT(filename,MOLECULE,BASIS,lupri)
implicit none
Character*(*)        :: filename
INTEGER,intent(in)            :: LUPRI
TYPE(MOLECULEINFO),intent(in) :: MOLECULE
TYPE(BASISINFO),intent(in)    :: BASIS
!
integer :: LUMOL,R,A,IATOMTYPE,natomsOfAtomtype,I,type,natomtypes,icharge
logical :: moleculefile_exsist,AUXBASIS,ADMMBASIS
real(realk) :: CHARGE
character(len=5) :: natomsstring
character(len=6) :: chargestring
character(len=3) :: chargestring2

INQUIRE(file=filename,EXIST=moleculefile_exsist) 
IF(moleculefile_exsist)THEN
   WRITE(lupri,*)'WARNING WRITE_MOLECULE_OUTPUT subroutine is overwriting file'
   WRITE(lupri,*)'with filename:',filename
ENDIF
LUMOL=-1
CALL LSOPEN(LUMOL,filename,'UNKNOWN','FORMATTED')

!we use ATOMBASIS format as default as all BASIS formats can be written 
!in ATOMBASIS format, but not all ATOMBASIS formats can be written in BASIS
WRITE(lumol,'(A)')'ATOMBASIS' 
WRITE(lumol,'(A)')'This file have been generated using the WRITE_MOLECULE_OUTPUT subroutine'
WRITE(lumol,'(A)')'Output is in Bohr'
natomtypes = BASIS%BINFO(RegBasParam)%natomtypes
IF(natomtypes.LT. 10)THEN
   WRITE(natomsstring,'(I1,A)')natomtypes,'    '
ELSEIF(natomtypes.LT. 100)THEN
   WRITE(natomsstring,'(I2,A)')natomtypes,'   '
ELSEIF(natomtypes.LT. 1000)THEN
   WRITE(natomsstring,'(I3,A)')natomtypes,'  '
ELSEIF(natomtypes.LT. 10000)THEN
   WRITE(natomsstring,'(I4,A)')natomtypes,' '
ELSEIF(natomtypes.LT. 100000)THEN
   WRITE(natomsstring,'(I5)')natomtypes,' '
ELSE
   WRITE(lupri,*)'Error in WRITE_MOLECULE_OUTPUT 100.000 atoms or more of a single type.'
ENDIF
iCharge=INT(MOLECULE%charge)
IF(icharge.LT.-9)THEN
   WRITE(chargestring2,'(A,I2)') '-',ABS(icharge)
ELSEIF(icharge.LT. 0)THEN
   WRITE(chargestring2,'(A,I1,A)') '-',ABS(icharge),' '
ELSEIF(icharge.LT. 10)THEN
   WRITE(chargestring2,'(I1,A)') icharge,'  '
ELSEIF(icharge.LT. 100)THEN
   WRITE(chargestring2,'(I2,A)') icharge,' '
ELSEIF(icharge.LT. 1000)THEN
   WRITE(chargestring2,'(I3)') icharge
ELSE
   WRITE(lupri,*)'Error in WRITE_MOLECULE_OUTPUT charge .gt. 1000?'
ENDIF

WRITE(lumol,'(A,A5,A,A3,A)')'Atomtypes=',natomsstring,'Charge=',chargestring2,' Nosymmetry'

R = BASIS%BINFO(RegBasParam)%Labelindex
A = BASIS%BINFO(AuxBasParam)%Labelindex
IF(BASIS%BINFO(AuxBasParam)%natomtypes.EQ. 0)THEN
   AUXBASIS = .FALSE.
ELSE
   AUXBASIS = .TRUE.
ENDIF
IF(BASIS%BINFO(ADMBasParam)%natomtypes.EQ. 0)THEN
   ADMMBASIS = .FALSE.
ELSE
   ADMMBASIS = .TRUE.
ENDIF

DO IATOMTYPE=1,BASIS%BINFO(REGBASPARAM)%natomtypes
   natomsOfAtomtype=0
   DO I=1,MOLECULE%natoms   
      IF(R.EQ. 0)THEN
         ICHARGE = INT(MOLECULE%ATOM(I)%CHARGE)
         type = BASIS%BINFO(REGBASPARAM)%CHARGEINDEX(ICHARGE)
      ELSE
         type = MOLECULE%ATOM(I)%IDtype(R)
      ENDIF
      IF(type.EQ.IATOMTYPE)THEN
         CHARGE = MOLECULE%ATOM(I)%CHARGE
         natomsOfAtomtype=natomsOfAtomtype+1
      ENDIF
   ENDDO
   IF(natomsOfAtomtype.LT. 10)THEN
      WRITE(natomsstring,'(I1,A)')natomsOfAtomtype,'    '
   ELSEIF(natomsOfAtomtype.LT. 100)THEN
      WRITE(natomsstring,'(I2,A)')natomsOfAtomtype,'   '
   ELSEIF(natomsOfAtomtype.LT. 1000)THEN
      WRITE(natomsstring,'(I3,A)')natomsOfAtomtype,'  '
   ELSEIF(natomsOfAtomtype.LT. 10000)THEN
      WRITE(natomsstring,'(I4,A)')natomsOfAtomtype,' '
   ELSEIF(natomsOfAtomtype.LT. 100000)THEN
      WRITE(natomsstring,'(I5)')natomsOfAtomtype,' '
   ELSE
      WRITE(lupri,*)'Error in WRITE_MOLECULE_OUTPUT 100.000 atoms or more of a single type.'
   ENDIF
   IF(charge.LT. 100)THEN
      WRITE(chargestring,'(F5.2,A)') charge,' '
   ELSEIF(charge.LT. 1000)THEN
      WRITE(chargestring,'(F6.2)') charge
   ELSE
      WRITE(lupri,*)'Error in WRITE_MOLECULE_OUTPUT charge .gt. 1000?'
   ENDIF

   IF(AUXBASIS)THEN
      IF (ADMMBASIS) THEN
          WRITE(lumol,'(A,A6,A,A5,6A)')'Charge=',chargestring,' Atoms=',natomsString,&
           & ' Basis=',TRIM(BASIS%BINFO(REGBASPARAM)%ATOMTYPE(type)%NAME),&
           & ' Aux=',TRIM(BASIS%BINFO(AUXBASPARAM)%ATOMTYPE(type)%NAME),&
           & ' ADMM=',TRIM(BASIS%BINFO(ADMBASPARAM)%ATOMTYPE(type)%NAME)
      ELSE
        WRITE(lumol,'(A,A6,A,A5,4A)')'Charge=',chargestring,' Atoms=',natomsString,&
             & ' Basis=',TRIM(BASIS%BINFO(REGBASPARAM)%ATOMTYPE(type)%NAME),&
             & ' Aux=',TRIM(BASIS%BINFO(AUXBASPARAM)%ATOMTYPE(type)%NAME)
      ENDIF
   ELSE
      IF (ADMMBASIS) THEN
        WRITE(lumol,'(A,A6,A,A5,4A)')'Charge=',Chargestring,' Atoms=',natomsString,&
           & ' Basis=',TRIM(BASIS%BINFO(REGBASPARAM)%ATOMTYPE(type)%NAME),&
           & ' ADMM=',TRIM(BASIS%BINFO(ADMBASPARAM)%ATOMTYPE(type)%NAME)
      ELSE
        WRITE(lumol,'(A,A6,A,A5,2A)')'Charge=',Chargestring,' Atoms=',natomsString,&
             & ' Basis=',TRIM(BASIS%BINFO(REGBASPARAM)%ATOMTYPE(type)%NAME)
      ENDIF
   ENDIF
   DO I=1,MOLECULE%natoms   
      CHARGE = MOLECULE%ATOM(I)%CHARGE
      IF(R.EQ. 0)THEN
         type = BASIS%BINFO(REGBASPARAM)%CHARGEINDEX(INT(CHARGE))
      ELSE
         type = MOLECULE%ATOM(I)%IDtype(R)
      ENDIF
      IF(type.EQ.IATOMTYPE)THEN
         WRITE(lumol,'(A4,1X,3F16.8)') MOLECULE%ATOM(I)%NAME,&
              & MOLECULE%ATOM(I)%CENTER(1),&
              & MOLECULE%ATOM(I)%CENTER(2),MOLECULE%ATOM(I)%CENTER(3)
      ENDIF
   ENDDO
ENDDO

CALL LSCLOSE(LUMOL,'KEEP')

END SUBROUTINE WRITE_MOLECULE_OUTPUT

END MODULE
