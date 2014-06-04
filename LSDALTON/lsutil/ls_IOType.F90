!> @file 
!> Contain the IOITEM which keeps track of files written to disk
MODULE io_type
use precision
use Matrix_module
!Integer,parameter :: maxFiles = 30000
Integer,parameter :: maxRecord = 2**28
Integer,parameter :: increment = 50

TYPE memMat
  TYPE(matrix)              :: mat
  Character(len=80)         :: filename
END TYPE memMat

TYPE memMatP
  TYPE(memMat),pointer  :: p
  TYPE(memMatP),pointer :: next
  TYPE(memMatP),pointer :: previous
END TYPE memMatP

TYPE IOITEM
Integer                   :: nallocFiles
Integer                   :: numFiles
Character(len=80),pointer :: filename(:)
Integer,pointer           :: IUNIT(:)
Logical,pointer           :: isOpen(:)
Integer                   :: nMemMat
TYPE(memMatP),pointer     :: first
TYPE(memMatP),pointer     :: current
Logical                   :: saveInMem
END TYPE IOITEM

contains

!Added to avoid "has no symbols" linking warning
subroutine io_TYPE_void()
end subroutine io_TYPE_void

END MODULE io_TYPE
