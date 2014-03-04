!> @file 
!> Contain the IOITEM which keeps track of files written to disk
MODULE io_type
use precision
!Integer,parameter :: maxFiles = 30000
Integer,parameter :: maxRecord = 2**28
Integer,parameter :: increment = 50

TYPE IOITEM
Integer       :: nallocFiles
Integer       :: numFiles
Character(len=80),pointer :: filename(:)
Integer,pointer       :: IUNIT(:)
Logical,pointer       :: isOpen(:)
END TYPE IOITEM

contains

!Added to avoid "has no symbols" linking warning
subroutine io_TYPE_void()
end subroutine io_TYPE_void

END MODULE io_TYPE
