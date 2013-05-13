module queue4_module
use Matrix_module
!   type Matrixp
!      TYPE(Matrix), pointer :: p moved to matrix_module
!   end type Matrixp

   type modFIFO4 !FIXME: could be done more elegant - make just one array where each element contains F, D, E, etx.
      integer :: offset
      integer :: queuesize
      TYPE(Matrixp), pointer :: xarray(:)
      TYPE(Matrixp), pointer :: parray(:)
      TYPE(Matrixp), pointer :: Aparray(:)
      TYPE(Matrixp), pointer :: Apparray(:)
      TYPE(Matrix), pointer  :: x_exp
      TYPE(Matrix), pointer  :: p_exp
      TYPE(Matrix), pointer  :: Ap_exp
      TYPE(Matrix), pointer  :: App_exp
      real(realk)            :: energy
      !f90 Workaround to make the pointer queue work!!
      type(Matrix) :: xmatrices(500), pmatrices(500), Apmatrices(500),Appmatrices(500)
      integer                 :: qcounter
   end type modFIFO4

end module queue4_module
