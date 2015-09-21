!> @file 
!> Contains the LSmatrix structure 
MODULE LSmatrix_type
use precision
TYPE LSMATRIX
INTEGER :: nrow
INTEGER :: ncol
real(realk), pointer   :: elms(:)
END TYPE LSMATRIX

type LSMatrixpointer
TYPE(LSMatrix), pointer :: p
end type LSMatrixpointer

private
public :: LSMatrixpointer,LSMATRIX
contains

!Added to avoid "has no symbols" linking warning
subroutine LSmatrix_type_void()
end subroutine LSmatrix_type_void

END MODULE LSmatrix_type
