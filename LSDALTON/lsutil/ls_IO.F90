!> @file 
!> Contain the IOITEM which keeps track of files written to disk
MODULE io
use io_type
use files
use precision
use Integralparameters
use matrix_module
use matrix_operations
use memory_handling
use mat3d_mod, only: read_mat3d_from_disk, mat3d
use molecule_type
use molecule_typetype
use LSTENSOR_OPERATIONSMOD, only: LSTENSOR
use mat3d_mod, only: write_mat3d_to_disk, mat3d
CONTAINS
!> \brief initialise the IOitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
SUBROUTINE io_init(IO)
implicit none
TYPE(IOITEM)  :: IO
IO%numFiles = 0
call io_alloc(IO,increment)
IO%isopen = .FALSE.
END SUBROUTINE io_init

!> \brief free the IOitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
SUBROUTINE io_free(IO)
implicit none
TYPE(IOITEM)  :: IO
IO%numFiles = 0
call io_dealloc(IO)
END SUBROUTINE io_free

!> \brief allocate the IOitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
!> \param nsize the number of files that should be allocated
SUBROUTINE io_alloc(IO,nsize)
implicit none
integer       :: nsize
TYPE(IOITEM)  :: IO
IO%nallocFiles = nsize
call mem_alloc(IO%filename,nsize)
call mem_alloc(IO%IUNIT,nsize)
call mem_alloc(IO%isopen,nsize)
END SUBROUTINE io_alloc

!> \brief deallocate the IOitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
SUBROUTINE io_dealloc(IO)
implicit none
TYPE(IOITEM)  :: IO
integer :: i
IO%nallocFiles = 0
do i=1,size(IO%isopen)
   IF(IO%isopen(i))THEN
      CALL LSCLOSE(IO%IUNIT(I),'KEEP')
      WRITE(*,'(1X,3A)') 'Error in io_dealloc. File: ',TRIM(IO%Filename(I)),' still open.'
      CALL lsQUIT('Error in io_dealloc. file still open!',-1)
   ENDIF
enddo
call mem_dealloc(IO%filename)
call mem_dealloc(IO%IUNIT)
call mem_dealloc(IO%isopen)
END SUBROUTINE io_dealloc

!> \brief copy the IOitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM, to be copied
!> \param IONEW the new copied IOITEM
SUBROUTINE COPY_IOITEM(IO,IONEW)
IMPLICIT NONE
TYPE(IOITEM),intent(in) :: IO
TYPE(IOITEM),intent(inout) :: IONEW
INTEGER      :: I
IOnew%numfiles=IO%numfiles
IF(IOnew%numfiles.GT.IOnew%nallocFiles)&
     &call lsquit('Programming Error in copy_ioitem, not initialized IONEW',-1)
DO I=1,IO%numfiles
 IOnew%filename(I)=IO%filename(I)
 IOnew%IUNIT(I)=IO%IUNIT(I)
 IOnew%isopen(I)=IO%isopen(I)
ENDDO

END SUBROUTINE COPY_IOITEM

!> \brief copy the IOitem
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM, to be copied
!> \param IONEW the new copied IOITEM
SUBROUTINE COPY_AND_ALLOC_IOITEM(newIO,IO)
IMPLICIT NONE
TYPE(IOITEM),intent(inout) :: newIO
TYPE(IOITEM),intent(in) :: IO
INTEGER      :: I

call io_alloc(newIO,IO%nallocFiles)
newIO%isopen = .FALSE.
call COPY_IOITEM(IO,newIO)

END SUBROUTINE COPY_AND_ALLOC_IOITEM

!> \brief add a filename to the IOITEM structue
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
!> \param Filename the filename to be added
!> \param LUPRI the logical unit number for the output file
SUBROUTINE io_add_filename(IO,Filename,LUPRI)
implicit none
TYPE(IOITEM)  :: IO
Character(80) :: Filename
Integer       :: LUPRI
!
Logical :: fileFound
Integer :: iFile,oldnallocFiles
TYPE(IOITEM)  :: TEMPIO

IF (io_file_exist(Filename,IO)) THEN
  WRITE(LUPRI,'(1X,3A)') 'Error in io_add_filename. File ',TRIM(Filename),&
   &                      ' allready exist in io-list!'
  CALL lsQUIT('Error in io_add_filename. Trying to add an existing filename.',lupri)
ENDIF

IF(IO%numFiles+1.GT.IO%nallocFiles)THEN
   oldnallocFiles = IO%nallocFiles
   call IO_alloc(TEMPIO,oldnallocFiles)
   CALL COPY_IOITEM(IO,TEMPIO)
   call IO_dealloc(IO)
   call IO_alloc(IO,oldnallocFiles+increment)
   CALL COPY_IOITEM(TEMPIO,IO)
   call IO_dealloc(TEMPIO)
   IF(IO%numFiles+1.GT.IO%nallocFiles)THEN
      WRITE(LUPRI,'(1X,2A)') 'Error in io_add_filename. something &
           & strange happend, last add file =',TRIM(Filename)
      CALL lsQUIT('Error in io_add_filename',lupri)
   ENDIF
endif

IO%numFiles = IO%numFiles + 1
IO%filename(IO%numFiles) = Filename
IO%IUNIT(IO%numFiles) = -1
IO%isOpen(IO%numFiles) = .FALSE.

END SUBROUTINE io_add_filename

!> \brief determines if the file exist in the IOITEM
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Filename the filename to be added
!> \param IO the IOITEM
LOGICAL FUNCTION io_file_exist(Filename,IO)
implicit none
Character(80) :: Filename
TYPE(IOITEM)  :: IO
!
Integer :: iFile
io_file_exist = .FALSE.
DO iFile=1,IO%numFiles
  IF (Filename.EQ.IO%Filename(iFile)) THEN
    io_file_exist = .TRUE.
    RETURN
  ENDIF
ENDDO
END FUNCTION io_file_exist

!> \brief determines logical unit number of filename
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param IO the IOITEM
!> \param Filename the filename
INTEGER FUNCTION io_iunit(IO,Filename)
implicit none
Character(80) :: Filename
TYPE(IOITEM)  :: IO
!
Integer :: iFile

IF (.not.io_file_exist(Filename,IO)) THEN
  CALL lsQUIT('Error in io_iunit, file does not exist!',-1)
ENDIF

DO iFile=1,IO%numFiles
  IF (Filename.EQ.IO%Filename(iFile)) EXIT
ENDDO
io_iunit = IO%iUnit(iFile)
END FUNCTION io_iunit

!> \brief open a filename
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_open(Filename,IO,LUPRI,LUERR,fileopen)
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR
TYPE(IOITEM)  :: IO
logical,optional :: fileopen
!
Integer :: IDUM,IUNIT,iFile
Logical :: fileFound,isOpen,LDUM,fileopen2,NewFile
fileopen2 = .TRUE.
IF(present(fileopen))THEN
   fileopen2=fileopen
ENDIF
fileFound = .FALSE.
DO iFile=1,IO%numFiles
  IF (Filename.EQ.IO%Filename(iFile)) THEN
    fileFound = .TRUE.
    IUNIT     = IO%IUNIT(iFile)
    isOpen    = IO%isOpen(iFile)
    EXIT
  ENDIF
ENDDO

IF (.NOT. fileFound) THEN
  WRITE(LUPRI,'(1X,2A)') 'Error in io_open. Could not find file: ',TRIM(Filename)
  CALL lsQUIT('Error in io_open. File not found!',lupri)
ENDIF

IF (isOpen) THEN
  WRITE(LUPRI,'(1X,3A)') 'Error in io_open. File: ',TRIM(Filename),' allready opened.'
  CALL lsQUIT('Error in io_open. Trying to open an allready opened file!',lupri)
ENDIF

IF(fileopen2)THEN
   CALL LSOPEN(IUNIT,FILENAME,'UNKNOWN','UNFORMATTED')
ELSE
   call matrixmembuf_Open(iunit,filename)
ENDIF
IO%isOpen(iFile) = .TRUE.
IO%IUNIT(iFile)  = IUNIT

END SUBROUTINE io_open

!> \brief close filename
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_close(Filename,IO,LUPRI,LUERR,fileclose)
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR
TYPE(IOITEM)  :: IO
logical,optional :: fileclose
!
Integer :: IDUM,IUNIT,iFile
Logical :: fileFound,isOpen,LDUM,fileExist,fileclose2
fileclose2 = .TRUE.
IF(present(fileclose))THEN
   fileclose2=fileclose
ENDIF

fileFound = .FALSE.
DO iFile=1,IO%numFiles
  IF (Filename.EQ.IO%Filename(iFile)) THEN
    fileFound = .TRUE.
    IUNIT     = IO%IUNIT(iFile)
    isOpen    = IO%isOpen(iFile)
    EXIT
  ENDIF
ENDDO

IF (.NOT. fileFound) THEN
  WRITE(LUPRI,'(1X,2A)') 'Error in io_close. Could not find file: ',TRIM(Filename)
  CALL lsQUIT('Error in io_close. File not found!',lupri)
ENDIF

IF (.NOT.isOpen) THEN
  WRITE(LUPRI,'(1X,3A)') 'Error in io_close. File: ',TRIM(Filename),' not opened.'
  CALL lsQUIT('Error in io_close. Trying to close a file that is not open!',lupri)
ENDIF

IF(fileclose2)THEN
   INQUIRE(file=FILENAME,exist=fileExist)
   IF (.NOT.fileExist) THEN
      WRITE(LUPRI,'(1X,3A)') 'Error in io_close. File: ',TRIM(Filename),' does not exist on disk!'
      CALL lsQUIT('Error in io_close. File missing!',lupri)
   ENDIF
   CALL LSCLOSE(IUNIT,'KEEP')
ELSE
   call matrixmembuf_Close(iunit,.FALSE.)
ENDIF
IO%isOpen(iFile) = .FALSE.
IO%IUNIT(iFile)  = -1

END SUBROUTINE io_close

SUBROUTINE io_get_CSidentifier(identifier,THR,mol1,mol2,CS,PS)
implicit none
Character(53)              :: identifier
Integer,intent(IN)         :: THR
TYPE(MOLECULEINFO),pointer :: mol1,mol2
LOGICAL,intent(IN)         :: CS,PS
!
Character(22) :: label1,label2

IF (ASSOCIATED(mol1)) THEN
   label1 = mol1%label
ELSE
   label1 = 'Empty_________________'
ENDIF
IF (ASSOCIATED(mol2)) THEN
   label2 = mol2%label
ELSE
   label2 = 'Empty_________________'
ENDIF
!$OMP CRITICAL (ifortwrite)
IF(THR.GT.9)THEN
  write(identifier,'(A2,I2,A1,A22,A1,A22,A1,L1,L1)')'CS',THR,'_',label1,'_',label2,'_',CS,PS
ELSE
  write(identifier,'(A3,I1,A1,A22,A1,A22,A1,L1,L1)')'CS0',THR,'_',label1,'_',label2,'_',CS,PS
ENDIF
!$OMP END CRITICAL (ifortwrite)

END SUBROUTINE io_get_CSidentifier


!> \brief obtain a filename based on input values used for screening
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Filename the filename
!> \param Identifier Character string to identify the filename
!> \param AO1 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 1
!> \param AO2 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 2
!> \param AO3 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 3
!> \param AO4 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 4
!> \param start1 start index for AO1
!> \param start2 start index for AO2
!> \param Oper operator label
!> \param intType 'Contracted or Primitive
!> \param FRAGMENT is this a fragment calculation
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_get_filename(Filename,Identifier,AO1,AO2,AO3,AO4,start1,start2,Oper,intType,FRAGMENT,LUPRI,LUERR)
implicit none
Character*(*)        :: Identifier
integer        :: AO1,AO2,AO3,AO4,Oper,intType
Integer        :: LUPRI,LUERR,start1,start2
Character(80)  :: Filename
logical        :: FRAGMENT
!
integer  :: AOstring(4)
Integer       :: iLen,iFilename,IAO,i,tmp(2)
Character(len=1)  :: STRING1
Character(len=2)  :: STRING2
Character(len=3)  :: STRING3
Character(len=4)  :: STRING4
Character(len=5)  :: STRING5
Character(len=6)  :: STRING6
Character(len=7)  :: OperatorString 
Character(len=10) :: InttypeString
call param_oper_Stringfromparam(OperatorString,Oper)
call param_inttype_Stringfromparam(InttypeString,IntType)
AOstring(1) = AO1
AOstring(2) = AO2
AOstring(3) = AO3
AOstring(4) = AO4
!$OMP CRITICAL (ifortwrite)
iLen = LEN(TRIM(Identifier))
Filename(1:iLen) = Identifier(1:iLen)
iFilename = iLen + 1
iLen = LEN(TRIM(InttypeString))
iLen = min(iLen,3)
Filename(iFilename:iFilename+iLen-1) = InttypeString(1:iLen)
iFilename = iFilename + iLen

iLen = LEN(TRIM(OperatorString))
iLen = min(iLen,3)
Filename(iFilename:iFilename+iLen-1) = OperatorString(1:iLen)
iFilename = iFilename + iLen
DO IAO=1,4
  IF (AOstring(IAO).EQ.AORegular) THEN
    Filename(iFilename:iFilename) = 'r'
    iFilename = iFilename + 1
  ELSEIF (AOstring(IAO).EQ.AOdfAux) THEN
    Filename(iFilename:iFilename) = 'd'
    iFilename = iFilename + 1
  ELSEIF (AOstring(IAO).EQ.AOdfCABS) THEN
    Filename(iFilename:iFilename) = 'c'
    iFilename = iFilename + 1
  ELSEIF (AOstring(IAO).EQ.AOdfJK) THEN
    Filename(iFilename:iFilename+1) = 'jk'
    iFilename = iFilename + 2
  ELSEIF (AOstring(IAO).EQ.AOVAL) THEN
    Filename(iFilename:iFilename) = 'v'
    iFilename = iFilename + 1
  ELSEIF (AOstring(IAO).EQ.AOEmpty) THEN
    Filename(iFilename:iFilename) = 'e'
    iFilename = iFilename + 1
  ELSEIF (AOstring(IAO).EQ.AONuclear) THEN
    Filename(iFilename:iFilename) = 'n'
    iFilename = iFilename + 1
  ELSEIF (AOstring(IAO).EQ.AOpCharge) THEN
    Filename(iFilename:iFilename) = 'p'
    iFilename = iFilename + 1
  ELSEIF (AOstring(IAO).EQ.AOelField) THEN
    Filename(iFilename:iFilename) = 'f'
    iFilename = iFilename + 1
  ELSE
    WRITE(LUPRI,'(1X,A,I1,2A)') 'Error! Wrong AOstring(',IAO,')=',AOstring(IAO),' in io_get_filename.'
    CALL lsQUIT('Error in io_get_filename Wrong AOstring',lupri)
  ENDIF
ENDDO

IF(FRAGMENT)THEN
   Filename(iFilename:iFilename) = 'F'
   iFilename = iFilename + 1
ENDIF

tmp(1)=start1
tmp(2)=start2
do I=1,2
   Filename(iFilename:iFilename) = 'S'
   iFilename = iFilename + 1
   IF(tmp(I) .LT. 10)THEN
      WRITE(STRING1,'(I1)') tmp(I)
      Filename(iFilename:iFilename) = STRING1
      iFilename = iFilename + 1
   ELSEIF(tmp(I) .LT. 100)THEN
      WRITE(STRING2,'(I2)') tmp(I)
      Filename(iFilename:iFilename+1) = STRING2
      iFilename = iFilename + 2
   ELSEIF(tmp(I) .LT. 1000)THEN
      WRITE(STRING3,'(I3)') tmp(I)
      Filename(iFilename:iFilename+2) = STRING3
      iFilename = iFilename + 3
   ELSEIF(tmp(I) .LT. 10000)THEN
      WRITE(STRING4,'(I4)') tmp(I)
      Filename(iFilename:iFilename+3) = STRING4
      iFilename = iFilename + 4
   ELSEIF(tmp(I) .LT. 100000)THEN
      WRITE(STRING5,'(I5)') tmp(I)
      Filename(iFilename:iFilename+4) = STRING5
      iFilename = iFilename + 5
   ELSEIF(tmp(I) .LT. 1000000)THEN
      WRITE(STRING6,'(I6)') tmp(I)
      Filename(iFilename:iFilename+5) = STRING6
      iFilename = iFilename + 6
   ELSE
      WRITE(LUPRI,'(1X,A,I1,2A)') 'startbasisfunction in io_get_filename is larger than 100000 not implemented'
      CALL lsQUIT('Error in io_get_filename startbasisfunction too large',lupri)
   ENDIF
enddo

iFilename=iFilename-1
IF (iFilename.GT. 80) THEN
  CALL lsQUIT('Error: iFilename > 80 in io_get_filename.',lupri)
ELSE
  DO i=iFilename+1,80
    Filename(i:i) = ' '
  ENDDO
ENDIF
!$OMP END CRITICAL (ifortwrite)
END SUBROUTINE io_get_filename

!> \brief write mat to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat to be written to disk
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_write_mat(Mat,Filename,IO,OnMaster,LUPRI,LUERR)
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR,n1,n2,n3,n4,n5
TYPE(matrix)  :: Mat
TYPE(IOITEM)  :: IO
logical       :: OnMaster
CALL io_open(Filename,IO,LUPRI,LUERR,OnMaster)
IF(.NOT.OnMaster)call matrixmembuf_Overwrite(io_iunit(IO,Filename),filename)
CALL mat_write_to_disk(io_iunit(IO,Filename),Mat,OnMaster)
CALL io_close(Filename,IO,LUPRI,LUERR,OnMaster)
END SUBROUTINE io_write_mat

!> \brief write tensor to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat 5 dim array to be written to disk
!> \param n1 the size of dimension 1
!> \param n2 the size of dimension 2
!> \param n3 the size of dimension 3
!> \param n4 the size of dimension 4
!> \param n5 the size of dimension 5
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_write(Mat,n1,n2,n3,n4,n5,Filename,IO,LUPRI,LUERR)
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR,n1,n2,n3,n4,n5
Real(realk)   :: Mat(n1,n2,n3,n4,n5)
TYPE(IOITEM)  :: IO
CALL io_open(Filename,IO,LUPRI,LUERR)
CALL io_write_tensor(Mat,n1,n2,1,1,1,io_iunit(IO,Filename))
CALL io_close(Filename,IO,LUPRI,LUERR)
END SUBROUTINE io_write

!> \brief write tensor to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat 5 dim array to be written to disk
!> \param n1 the size of dimension 1
!> \param n2 the size of dimension 2
!> \param n3 the size of dimension 3
!> \param n4 the size of dimension 4
!> \param n5 the size of dimension 5
!> \param IUNIT the logical unit number for the file to write to
SUBROUTINE io_write_tensor(Mat,n1,n2,n3,n4,n5,IUNIT)
implicit none
Integer       :: n1,n2,n3,n4,n5,IUNIT
Real(realk)   :: Mat(n1,n2,n3,n4,n5)
!
Integer :: i1,i2,i3,i4,i5
Integer :: IDUM
Logical :: LDUM
WRITE(IUNIT) n1,n2,n3,n4,n5
DO i5=1,n5
  DO i4=1,n4
    DO i3=1,n3
      DO i2=1,n2
        WRITE(IUNIT) (Mat(i1,i2,i3,i4,i5),i1=1,n1)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE io_write_tensor

!!$!> \brief write lstensor to disk
!!$!> \author S. Reine and T. Kjaergaard
!!$!> \date 2010
!!$!> \param tensor the lstensor to be written to disk  
!!$!> \param Filename the filename
!!$!> \param IO the IOITEM
!!$!> \param LUPRI the logical unit number for the output file
!!$!> \param LUERR the logical unit number for the error file
!!$SUBROUTINE io_write_lstensor(tensor,Filename,IO,LUPRI,LUERR)
!!$implicit none
!!$Character(80) :: Filename
!!$Integer       :: LUPRI,LUERR
!!$type(lstensor) :: tensor
!!$TYPE(IOITEM)  :: IO
!!$CALL io_open(Filename,IO,LUPRI,LUERR)
!!$CALL write_lstensor_to_disk(tensor,io_iunit(IO,Filename),lupri)
!!$CALL io_close(Filename,IO,LUPRI,LUERR)
!!$END SUBROUTINE io_write_lstensor

!> \brief read mat to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat to be read from file 
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_read_mat(Mat,Filename,IO,OnMaster,LUPRI,LUERR)
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR,n1,n2,n3,n4,n5
TYPE(matrix)  :: Mat
TYPE(IOITEM)  :: IO
logical       :: OnMaster
CALL io_open(Filename,IO,LUPRI,LUERR,OnMaster)
CALL mat_read_from_disk(io_iunit(IO,Filename),Mat,OnMaster)
CALL io_close(Filename,IO,LUPRI,LUERR,OnMaster)
END SUBROUTINE io_read_mat

!> \brief read 5 dim array to disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat to be read from file 
!> \param n1 the size of dimension 1
!> \param n2 the size of dimension 2
!> \param n3 the size of dimension 3
!> \param n4 the size of dimension 4
!> \param n5 the size of dimension 5
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_read(Mat,n1,n2,n3,n4,n5,Filename,IO,LUPRI,LUERR)
implicit none
Character(80) :: Filename
Integer       :: LUPRI,LUERR,n1,n2,n3,n4,n5
Real(realk)   :: Mat(n1,n2,n3,n4,n5)
TYPE(IOITEM)  :: IO
CALL io_open(Filename,IO,LUPRI,LUERR)
CALL io_read_tensor(Mat,n1,n2,1,1,1,io_iunit(IO,Filename))
CALL io_close(Filename,IO,LUPRI,LUERR)
END SUBROUTINE io_read

!> \brief read 5 dim array from disk
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param mat to be read from file 
!> \param n1 the size of dimension 1
!> \param n2 the size of dimension 2
!> \param n3 the size of dimension 3
!> \param n4 the size of dimension 4
!> \param n5 the size of dimension 5
!> \param IUNIT the logical unit number for file from which to read
SUBROUTINE io_read_tensor(Mat,n1,n2,n3,n4,n5,IUNIT)
implicit none
Integer       :: n1,n2,n3,n4,n5,IUNIT
Real(realk)   :: Mat(n1,n2,n3,n4,n5)
!
Integer :: i1,i2,i3,i4,i5
Integer :: IDUM
Logical :: LDUM
READ(IUNIT) i1,i2,i3,i4,i5
IF ((i1.NE.n1).OR.(i2.NE.n2).OR.(i3.NE.n3).OR.(i4.NE.n4).OR.(i5.NE.n5)) THEN
  CALL lsQUIT('Dimension error in io_read_tensor!',-1)
ENDIF
DO i5=1,n5
  DO i4=1,n4
    DO i3=1,n3
      DO i2=1,n2
        READ(IUNIT) (Mat(i1,i2,i3,i4,i5),i1=1,n1)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE io_read_tensor

!!$!> \brief read lstensor from disk
!!$!> \author S. Reine and T. Kjaergaard
!!$!> \date 2010
!!$!> \param tensor the lstensor to be written to disk  
!!$!> \param Filename the filename
!!$!> \param IO the IOITEM
!!$!> \param LUPRI the logical unit number for the output file
!!$!> \param LUERR the logical unit number for the error file
!!$SUBROUTINE io_read_lstensor(tensor,Filename,IO,LUPRI,LUERR)
!!$implicit none
!!$Character(80) :: Filename
!!$Integer       :: LUPRI,LUERR
!!$type(lstensor) :: tensor
!!$TYPE(IOITEM)  :: IO
!!$!$OMP CRITICAL
!!$CALL io_open(Filename,IO,LUPRI,LUERR)
!!$CALL read_lstensor_from_disk(tensor,io_iunit(IO,Filename),lupri)
!!$CALL io_close(Filename,IO,LUPRI,LUERR)
!!$!$OMP END CRITICAL
!!$END SUBROUTINE io_read_lstensor

!Type(IOitem)         :: IO
!IUNIT = -1
!CALL LSOPEN(IUNIT,FILENAME,'UNKNOWN','SEQUENTIAL','UNFORMATTED',IDUM,LDUM)

!> \brief read mat3d from disk
!> \author S. Reine
!> \date 2011-01-13
!> \param mat  The MAT3D's to be written to disk  
!> \param nmat The number of MAT3D's to be written to disk
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_read_mat3d(mat,nmat,Filename,IO,LUPRI,LUERR)
implicit none
Character(80)  :: Filename
Integer        :: LUPRI,LUERR,nmat
type(MAT3D)    :: mat(nmat)
TYPE(IOITEM)   :: IO
!
Integer :: imat,nmat_file
CALL io_open(Filename,IO,LUPRI,LUERR)
read(io_iunit(IO,Filename)) nmat_file
IF (nmat_file.NE.nmat) CALL LSQUIT('io_read_mat3d',LUPRI)
DO imat=1,nmat
  CALL read_mat3d_from_disk(mat(imat),io_iunit(IO,Filename),lupri)
ENDDO
CALL io_close(Filename,IO,LUPRI,LUERR)
END SUBROUTINE io_read_mat3d

!> \brief write mat3d from disk
!> \author S. Reine
!> \date 2011-01-13
!> \param mat the MAT3D to be written to disk  
!> \param Filename the filename
!> \param IO the IOITEM
!> \param LUPRI the logical unit number for the output file
!> \param LUERR the logical unit number for the error file
SUBROUTINE io_write_mat3d(mat,nmat,Filename,IO,LUPRI,LUERR)
implicit none
Character(80)  :: Filename
Integer        :: LUPRI,LUERR,nmat
type(MAT3D)    :: mat(nmat)
TYPE(IOITEM)   :: IO
!
Integer :: imat
CALL io_open(Filename,IO,LUPRI,LUERR)
write(io_iunit(IO,Filename)) nmat
DO imat=1,nmat
  CALL write_mat3d_to_disk(mat(imat),io_iunit(IO,Filename))
ENDDO
CALL io_close(Filename,IO,LUPRI,LUERR)
END SUBROUTINE io_write_mat3d

  
END MODULE io


!> \brief write vector 
!> \author S. Reine
!> \date 2010
!> \param iunit the logical unit number for the file, from which to read
!> \param vector the vectro th write
!> \param N the size of the vector
SUBROUTINE ls_write(IUNIT,vector,N)
use io, only: maxRecord
use precision
implicit none
Integer,intent(IN)     :: IUNIT,N
Real(realk),intent(IN) :: vector(N)
!
Integer :: nBUF,startBUF,endBUF,I
nBUF     = n/maxRecord + 1
startBUF = 1
DO I=1,nBUF
  endBUF = MIN(n,startBUF+maxRecord-1)
  WRITE(IUNIT) vector(startBUF:endBUF)
  startBUF = startBUF + maxRecord
ENDDO
END SUBROUTINE ls_write

!> \brief read vector
!> \author S. Reine
!> \date 2010
!> \param iunit the logical unit number for the file, from which to read
!> \param vector the vectro th write
!> \param N the size of the vector
SUBROUTINE ls_read(IUNIT,vector,N)
use io, only: maxRecord
use precision
implicit none
Integer,intent(IN)      :: IUNIT,N
Real(realk),intent(OUT) :: vector(N)
!
Integer :: nBUF,startBUF,endBUF,I
nBUF     = n/maxRecord + 1
startBUF = 1
DO I=1,nBUF
  endBUF = MIN(n,startBUF+maxRecord-1)
  READ(IUNIT) vector(startBUF:endBUF)
  startBUF = startBUF + maxRecord
ENDDO
END SUBROUTINE ls_read

