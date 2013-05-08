!> @file
!> Contains queue type definition module.

!> \brief Queue type definition.
!> \author Stinne Host
!> \date 2007
module queue_module
use matrix_module

   !> Pointer to an integer. Necessary if we want arrays of pointers!
   type integerp
      integer, pointer :: p
   end type integerp

   !> \brief Queue type definition.
   !> \author Stinne Host
   !> \date 2007
   !>
   !> modFIFO stands for modified First-In-First-Out. FIFO is a standard way to
   !> store a number of old matrices (or vectors), where the oldest matrix is
   !> thrown away when there is no more space. This is a modified type of FIFO,
   !> because it contains two queues of matrices instead of one, and treats one
   !> member (usually the newest) in a special way (the expansion point). \n
   !> The matrices are labelled F and D for Fock and density, but the queue is
   !> completely general and could be used for anything. \n
   !> We have to use the type(matrixp) which only contains one item, a type(matrix) (see matrix_module): \n
   !>   type Matrixp \n
   !>      TYPE(Matrix), pointer :: p  \n
   !>   end type Matrixp \n
   !> This is a workaround because we are not allowed by fortran90 to make a
   !> pointer to an array of type(matrix), because it contains multiple items.
   !> The type(matrixp), on the other hand, contains only one member - a type(matrix).\n
   !> The use of Dmatrices, Fmatrices and qcounter is not pretty! This is the only
   !> way I have been able to make the pointer queue work. I suspect it is because
   !> pointers is really just an add-on to fortran, and they don't work optimally.
   !> But it is of course also possible that I'm just an idiot. Suggestions are
   !> welcome!
   !>
   type modFIFO 
      !> Indicates how many matrices we have in the queue. Eq. to queuesize-1 if queue is full
      integer :: offset
      !> Indicates how many matrices can be saved, including expansion point
      integer :: queuesize
      !> The first queue (densities, or "array1")
      TYPE(Matrixp), pointer :: Darray(:)
      !> The second queue (Fock matrices, or "array2")
      TYPE(Matrixp), pointer :: Farray(:)
      !> Expansion point for the first queue
      TYPE(Matrix), pointer  :: D_exp
      !> Expansion point for the second queue
      TYPE(Matrix), pointer  :: F_exp
      !> A real to hold info for expansion point (used for energy, but could be anything).
      real(realk)            :: energy
      !> f90 workaround to make the pointer queue work!!
      type(Matrix) :: Dmatrices(200)
      !> f90 workaround to make the pointer queue work!!
      type(Matrix) :: Fmatrices(200)
      !> f90 workaround to make the pointer queue work!!
      integer      :: qcounter
   !Disk variables:
   !===============
      !> True if queue should be stored on disk
      logical          :: disk
      !> LUNs for the first queue (densities, or "array1") - if storing on disk
      type(integerp), pointer :: iDarray(:)
      !> LUNs for the second queue (Fock matrices, or "array2") - if storing on disk
      type(integerp), pointer :: iFarray(:)
      !> LUN for expansion point for the first queue - if storing on disk
      integer, pointer :: iD_exp
      !> LUN for expansion point for the second queue - if storing on disk
      integer, pointer :: iF_exp
      !> Temporary storage
      TYPE(Matrix) :: tempD
      !> Temporary storage
      TYPE(Matrix) :: tempF
      !> f90 workaround to make the pointer queue work!!
      integer :: iDmatrices(200)
      !> f90 workaround to make the pointer queue work!!
      integer :: iFmatrices(200)
   end type modFIFO

contains

!Added to avoid "has no symbols" linking warning
subroutine queue_module_void()
end subroutine queue_module_void

end module queue_module
