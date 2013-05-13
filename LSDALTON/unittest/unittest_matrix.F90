module unittest_mat
use unittest_mat_dense

contains

   subroutine MatrixTest

   implicit none

      !Test the different matrix types:
      !1. Dense
         call MatrixTestDense

      !2. Unrestricted dense
      !  call MatrixTestUnres 

      !3. CSR
      !  call MatrixTestCSR         

   end subroutine MatrixTest

end module unittest_mat



