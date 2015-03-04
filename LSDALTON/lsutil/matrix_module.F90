!> @file
!> Contains matrix type definition module.

!> \brief Matrix type definitions.
!>
!> General rules:
!> NEVER put a type(Matrix) as intent(out), this will on some platforms
!>       make the pointer 
!>       disassociated entering the routine, and memory already 
!>       allocated for the matrix will be lost. \n
!> NEVER think that e.g. A%elms = matrix will copy the matrix elements from 
!>       matrix to A%elms, it will only associate the pointer with the array
!>       matrix. \n
!> Use mat_assign for A = B operation
!> ALWAYS and ONLY call mat_free on a matrix you have initialized with mat_init
!>
MODULE Matrix_module  
   use precision
   private
   public :: mat_lu,mat_info,mat_mem_monitor, Matrix, Matrixp, MAT_NULLIFY, &
        & mat_init_magic_value
   !> Logical unit number of file LSDALTON.OUT used by matrix operation modules
   integer, save :: mat_lu
   !> True if various info from matrix module should be printed
   logical, save :: mat_info
   !> True if number of allocated matrices should be monitored
   logical, save :: mat_mem_monitor
   !> The general matrix type that is used throughout LSDALTON
   TYPE Matrix
      !> number of rows
      INTEGER :: nrow
      !> number of columns
      INTEGER :: ncol
!      !> room for various flags/information for the matrix
!      integer, dimension(:),pointer :: idata
!      !> bsm permutation pointer
!      integer, dimension(:),pointer :: permutation

      !> pointer to double precision matrix element storage
      real(realk), pointer   :: elms(:)
      !> pointer to double precision matrix element storage for beta part
      real(realk), pointer   :: elmsb(:)
!      !> pointer to complex matrix element storage
!      complex(realk), pointer :: celms(:)
!      !> pointer to complex matrix element storage for beta part
!      complex(realk), pointer :: celmsb(:)
      !> room for any integer auxiliary information
      integer, pointer     :: iaux(:)
      !> room for any real auxiliary information
      real(realk), pointer :: raux(:)
!CSR INFO 
      real(realk),pointer :: val(:)
      integer,pointer :: col(:)
      integer,pointer :: row(:)
      integer :: nnz

!      !> flag for complex elements
      logical :: complex

      !> tag used to spot accidental use of uninitialized and memory-corrupted
      !> matrices, and fail with something other than 'segmentation fault'.
      !> Tag is set by mat_init and cleared by mat_free
      integer :: init_magic_tag
      !> pointer to self, set by init to mark the matrix' correct location.
      !> Used to distinguish init'ed matrices from copy-in'ed matrices like
      !> "call subr((/A,B,C/))". Set by mat_init and cleared by mat_free
      type(Matrix), pointer :: init_self_ptr
!SCALAPACK INFO
#ifdef VAR_SCALAPACK
      !> number of rows
      INTEGER :: localnrow
      !> number of columns
      INTEGER :: localncol
      integer, pointer :: addr_on_grid(:,:)
      real(REALK), pointer :: p(:,:)
#endif
      !> Parallel Distributed Memory Matrix Identification 
      INTEGER :: PDMID
   END TYPE Matrix

   !> Pointer to a type(matrix). Necessary if we want arrays of derived types!
   type Matrixp
      TYPE(Matrix), pointer :: p
   end type Matrixp

   !> Large random number >-2^31, <+2^31 (but also works on 64bit),
   !> not likely to be found in memory by accident. Mat_init sets
   !> type(matrix)%init_magic_tag to this value, all matrix operations
   !> will expect it to be there, and mat_free clears it
   integer, parameter :: mat_init_magic_value = -1981879812
   CONTAINS

     SUBROUTINE MAT_NULLIFY(MAT)
       type(matrix) :: MAT
!       nullify(MAT%idata)
!       nullify(MAT%permutation)
       nullify(MAT%elms)
       nullify(MAT%elmsb)
!       nullify(MAT%celms)
!       nullify(MAT%celmsb)
       nullify(MAT%iaux)
       nullify(MAT%raux)
       nullify(MAT%val)
       nullify(MAT%col)
       nullify(MAT%row)
       nullify(MAT%init_self_ptr)
#ifdef VAR_SCALAPACK
       nullify(MAT%addr_on_grid)
       nullify(MAT%p)
#endif
       MAT%init_magic_tag = 0 !clear magic tag
     END SUBROUTINE MAT_NULLIFY

END MODULE Matrix_module

