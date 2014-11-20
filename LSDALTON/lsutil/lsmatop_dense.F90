!> @file
!> Contains a dense matrix module used in LSint directory 
!> \brief Contains matrix operation routines for type(lsmatrix).
!>
module lsmatrix_operations_dense
  use precision
  use LSmatrix_type

CONTAINS
!> \brief Initialize a type(lsmatrix)
!> \author T. Kjaergaard
!> \date 2010
!> \param a type(lsmatrix) that should be initialized
!> \param nrow Number of rows for a
!> \param ncol Number of columns for a
  subroutine lsmat_dense_init(A,nrow,ncol)
    use memory_handling
    implicit none
    TYPE(Lsmatrix) :: A
    integer, intent(in) :: nrow, ncol
    integer(kind=long) :: nsize
    !A%elms should not be associated with anything.
    !if so, the memory will be lost!!
    A%nrow = nrow
    A%ncol = ncol
    allocate(A%elms(A%nrow * A%ncol))
    nsize = size(A%elms,KIND=long)*mem_realsize
    call mem_allocated_mem_type_matrix(nsize)
  end subroutine lsmat_dense_init

!> \brief Free a type(lsmatrix) that has been initialized
!> \author T. Kjaergaard
!> \date 2010
!> \param a type(lsmatrix) that should be freed
  subroutine lsmat_dense_free(a)
    use memory_handling
    implicit none
    TYPE(Lsmatrix) :: a
    integer(kind=long) :: nsize
    if (.not.ASSOCIATED(a%elms)) then
       print*,'memory previously released!!'
       STOP 'Error in lsmat_dense_free - memory previously released'
    endif
    nsize = SIZE(a%elms,KIND=long)*mem_realsize
    call mem_deallocated_mem_type_matrix(nsize)
    DEALLOCATE(a%elms)
    NULLIFY(a%elms)
  end subroutine lsmat_dense_free

!> \brief Print a type(lsmatrix) to file in pretty format
!> \author T. Kjaergaard
!> \date 2010
!> \param a The type(lsmatrix) that should be printed
!> \param i_row1 Print starting from this row
!> \param i_rown Print ending at this row
!> \param j_col1 Print starting from this column
!> \param j_coln Print ending at this column
!> \param lu Print to file with this logical unit number
  subroutine lsmat_dense_print(a, i_row1, i_rown, j_col1, j_coln, lu)
    implicit none
    TYPE(Lsmatrix),intent(in) :: a
    integer, intent(in)     :: i_row1, i_rown, j_col1, j_coln, lu
    !to be found in lsutil/common_utilities.F90
    call LS_OUTPUT(a%elms, i_row1, i_rown, j_col1, j_coln, a%nrow, a%ncol, 1, lu)
  end subroutine lsmat_dense_print

!> \brief Set a type(lsmatrix) A to zero
!> \author T. Kjaergaard
!> \date 2010
!> \param A Input/output matrix which should be set to zero
  subroutine lsmat_dense_zero(A)
    implicit none
    type(Lsmatrix), intent(inout) :: A
    integer :: i

    do i = 1,A%nrow*A%ncol
       A%elms(i) = 0.0E0_realk
    enddo

  end subroutine lsmat_dense_zero

end module lsmatrix_operations_dense
