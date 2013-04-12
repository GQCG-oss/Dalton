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
    !to be found in pdpack/printpkg.F
    call OUTPUT(a%elms, i_row1, i_rown, j_col1, j_coln, a%nrow, a%ncol, 1, lu)
  end subroutine lsmat_dense_print

!> \brief Transpose a type(lsmatrix).
!> \author T. Kjaergaard
!> \date 2010
!> \param a The type(lsmatrix) that should be transposed
!> \param b The transposed output type(lsmatrix).
!>
!> Usage discouraged! If what you want is to multiply your transposed
!> matrix with something else, you should instead use mat_mul with the
!> transpose flag 'T'. This is much more efficient than transposing first 
!> and then multiplying.
  subroutine lsmat_dense_trans(a,b)
    implicit none
    type(lsmatrix),intent(in)   :: a
    type(lsmatrix)              :: b   !output  
    integer                   :: i, j

    do j = 1,a%ncol
       do i = 1,a%nrow
          b%elms(a%nrow*(i-1)+j) = a%elms(a%nrow*(j-1)+i)
       enddo
    enddo
  end subroutine lsmat_dense_trans

!> \brief Copy a type(lsmatrix).
!> \author T. Kjaergaard
!> \date 2010
!> \param a The type(lsmatrix) that should be copied
!> \param b The copy output type(lsmatrix)
  subroutine lsmat_dense_assign(a,b)
    implicit none
    TYPE(Lsmatrix), INTENT(INOUT) :: a
    TYPE(Lsmatrix), INTENT(IN)    :: b
    integer                     :: i

    !    do i = 1,a%nrow*a%ncol
    !      a%elms(i) = b%elms(i)
    !    enddo

    i = a%nrow*a%ncol
    call dcopy (i,b%elms,1,a%elms,1)

  end subroutine lsmat_dense_assign

!> \brief Copy and scale a type(lsmatrix).
!> \author T. Kjaergaard
!> \date 2010
!> \param alpha The scaling parameter
!> \param a The type(lsmatrix) that should be copied
!> \param b The scaled output type(lsmatrix).
  subroutine lsmat_dense_copy(alpha,a,b)
    implicit none
    REAL(REALK),  INTENT(IN)    :: alpha
    TYPE(Lsmatrix), INTENT(IN)    :: a
    TYPE(Lsmatrix), INTENT(INOUT) :: b

    call lsmat_dense_assign(b,a)
    if (ABS(alpha-1.0E0_realk).GT.1.0E-15_realk) call lsmat_dense_scal(alpha,b)

  end subroutine lsmat_dense_copy

!> \brief Makes the trace of a type(lsmatrix).
!> \param a The type(lsmatrix) we want the trace of
!> \return The trace of a
!> \author T. Kjaergaard
!> \date 2010
  function lsmat_dense_Tr(a)
    implicit none
    TYPE(Lsmatrix), intent(IN) :: a
    REAL(realk) :: lsmat_dense_tr 
    integer :: i

    lsmat_dense_tr = 0E0_realk
    do i = 1,a%nrow
       lsmat_dense_tr = lsmat_dense_tr + a%elms((a%Nrow+1)*i-a%Nrow)
    enddo

  end function lsmat_dense_Tr

!> \brief Make the trace of the product of type(lsmatrix) A and B.
!> \author T. Kjaergaard
!> \date 2010
!> \param a The first type(lsmatrix) factor
!> \param b The second type(lsmatrix) factor
!> \return Tr(a*b)
  function lsmat_dense_TrAB(a,b)
    implicit none
    TYPE(Lsmatrix), intent(IN) :: a,b
    REAL(realk) :: lsmat_dense_trAB 
    real(realk), external :: ddot
    integer :: j

    lsmat_dense_TrAB = 0.0E0_realk
    do j = 1,a%ncol
       !      do i = 1,a%nrow
       !        lsmat_dense_TrAB = lsmat_dense_TrAB + a%elms(a%nrow*(j-1)+i)*&
       !                                         &b%elms(b%nrow*(i-1)+j)
       !      enddo
       lsmat_dense_TrAB = lsmat_dense_TrAB + &
            & ddot(a%nrow,a%elms(a%nrow*(j-1)+1),1,b%elms(j),a%ncol)
    enddo

  end function lsmat_dense_TrAB

!> \brief Make c = alpha*ab + beta*c, where a,b,c are type(lsmatrix) and alpha,beta are parameters
!> \author T. Kjaergaard
!> \date 2010
!> \param a The first type(lsmatrix) factor
!> \param b The second type(lsmatrix) factor
!> \param transa 'T'/'t' if a should be transposed, 'N'/'n' otherwise
!> \param transb 'T'/'t' if b should be transposed, 'N'/'n' otherwise
!> \param alpha The alpha parameter
!> \param beta The beta parameter
!> \param c The output type(lsmatrix)
  subroutine lsmat_dense_mul(a,b,transa,transb,alpha,beta,c) 
    implicit none
    TYPE(Lsmatrix), intent(IN) :: a, b
    character, intent(in)    :: transa, transb
    REAL(realk), INTENT(IN)  :: alpha, beta
    TYPE(Lsmatrix)             :: c  
    INTEGER                  :: m,n,k

    if (transa == 'n' .or. transa == 'N') then
       m = a%nrow
       k = a%ncol
    elseif (transa == 't' .or. transa == 'T') then
       m = a%ncol
       k = a%nrow
    else
       print*,'unknown format in lsmat_dense_mul'
       STOP 'unknown format in lsmat_dense_mul'
    endif
    if (transb == 'n' .or. transb == 'N') then
       n = b%ncol
       if (b%nrow /= k) STOP 'weird'
    elseif (transb == 't' .or. transb == 'T') then
       n = b%nrow
       if (b%ncol /= k) STOP 'weird'
    endif
    !write (lsmat_lu,*) 'lsmatrix a'
    !call lsmat_dense_print(a, 1, a%nrow, 1, a%ncol, lsmat_lu)
    !write (lsmat_lu,*) 'lsmatrix b'
    !call lsmat_dense_print(b, 1, b%nrow, 1, b%ncol, lsmat_lu)
    !write (lsmat_lu,*) 'lsmatrix c'
    !call lsmat_dense_print(c, 1, c%nrow, 1, c%ncol, lsmat_lu)
    call DGEMM(transa,transb,m,n,k,alpha,&
         &a%elms,a%nrow,b%elms,b%nrow,beta,c%elms,c%nrow)

  end subroutine lsmat_dense_mul

!> \brief Make c = alpha*a + beta*b, where a,b are type(lsmatrix) and alpha,beta are parameters
!> \author T. Kjaergaard
!> \date 2010
!> \param a The first type(lsmatrix) 
!> \param alpha The alpha parameter
!> \param b The second type(lsmatrix) 
!> \param beta The beta parameter
!> \param c The output type(lsmatrix)
  subroutine lsmat_dense_add(alpha,a,beta,b,c)
    implicit none
    TYPE(Lsmatrix), intent(IN) :: a, b
    REAL(realk), INTENT(IN)  :: alpha, beta
    TYPE(Lsmatrix)             :: c
    integer                  :: i

    !    do i = 1,a%nrow*a%ncol
    !      c%elms(i) = alpha*a%elms(i) + beta*b%elms(i)
    !    enddo

    call lsmat_dense_copy(alpha,a,c)
    call lsmat_dense_daxpy(beta,b,c)

  end subroutine lsmat_dense_add

!> \brief Make Y = alpha*X + Y where X,Y are type(lsmatrix) and a is a parameter
!> \author T. Kjaergaard
!> \date 2010
!> \param alpha The alpha parameter
!> \param X The input type(lsmatrix) 
!> \param Y The input/output type(lsmatrix) 
  subroutine lsmat_dense_daxpy(alpha,x,y)
    implicit none
    real(realk),intent(in)       :: alpha
    TYPE(Lsmatrix), intent(IN)     :: X
    TYPE(Lsmatrix), intent(INOUT)  :: Y
    integer                      :: i

    i = x%nrow*x%ncol
    call daxpy(i,alpha,x%elms,1,y%elms,1)

  end subroutine lsmat_dense_daxpy

!> \brief Make the dot product of type(lsmatrix) a and b.
!> \author T. Kjaergaard
!> \date 2010
!> \param a The first type(lsmatrix) factor
!> \param b The second type(lsmatrix) factor
!> \return The dot product of a and b
  function lsmat_dense_dotproduct(a,b)
    implicit none
    TYPE(Lsmatrix), intent(IN) :: a,b
    REAL(realk) :: lsmat_dense_dotproduct
    real(realk), external :: ddot
    integer     :: i

    !    lsmat_dense_dotproduct = 0.0E0_realk
    !    do i = 1,a%nrow*a%ncol
    !      lsmat_dense_dotproduct = lsmat_dense_dotproduct + a%elms(i)*b%elms(i)
    !    enddo

    i = a%nrow*a%ncol
    lsmat_dense_dotproduct = ddot(i,a%elms,1,b%elms,1)

  end function lsmat_dense_dotproduct

!> \brief Make the dot product of type(lsmatrix) a with itself.
!> \author T. Kjaergaard
!> \date 2010
!> \param a The type(lsmatrix) input
!> \return The dot product of a with itself
  function lsmat_dense_sqnorm2(a)
    implicit none
    TYPE(Lsmatrix), intent(IN) :: a
    REAL(realk) :: lsmat_dense_sqnorm2
    !    integer     :: i

    !    lsmat_dense_sqnorm2 = 0.0E0_realk
    !    do i = 1,a%nrow*a%ncol
    !      lsmat_dense_sqnorm2 = lsmat_dense_sqnorm2 + a%elms(i)*a%elms(i)
    !    enddo
    lsmat_dense_sqnorm2 = lsmat_dense_dotproduct(a,a)

  end function lsmat_dense_sqnorm2

!> \brief Scale a type(lsmatrix) A by a scalar alpha
!> \author T. Kjaergaard
!> \date 2010
!> \param alpha Scaling parameter
!> \param A Input/output matrix which we want to scale
  subroutine lsmat_dense_scal(alpha,A)
    implicit none
    real(realk), intent(in) :: alpha
    type(Lsmatrix), intent(inout) :: A
    integer :: i

    !   do i = 1,A%nrow*A%ncol
    !     A%elms(i) = A%elms(i) * alpha
    !   enddo
    i = A%nrow*A%ncol
    call dscal(i,alpha,A%elms,1)

  end subroutine lsmat_dense_scal

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

!> \brief Write a type(lsmatrix) to disk.
!> \author T. Kjaergaard
!> \date 2010
!> \param iunit Logical unit number of file which matrix should be written to
!> \param A Matrix which should be written on disk
  subroutine lsmat_dense_WRITE_TO_DISK(iunit,A)
    implicit none
    integer, intent(in) :: iunit
    type(Lsmatrix), intent(in) :: A
    integer :: i

    if (.not.ASSOCIATED(A%elms)) STOP 'A in lsmat_dense_WRITE_TO_DISK non-existant'

    WRITE(iunit) A%Nrow, A%Ncol
    WRITE(iunit)(A%elms(I),I=1,A%nrow*A%ncol)

  end subroutine lsmat_dense_WRITE_TO_DISK

!> \brief Read a type(lsmatrix) from disk.
!> \author T. Kjaergaard
!> \date 2010
!> \param iunit Logical unit number of file from which matrix should be read
!> \param A Output matrix which is read from disk
  subroutine lsmat_dense_READ_FROM_DISK(iunit,A)
    implicit none
    integer, intent(in) :: iunit
    type(Lsmatrix), intent(inout) :: A
    integer :: i
    READ(iunit) A%Nrow, A%Ncol
    READ(iunit)(A%elms(I),I=1,A%nrow*A%ncol)
  end subroutine lsmat_dense_READ_FROM_DISK

!> \brief Returns sum of all elements of matrix.
!> \param A The input lsmatrix 
!> \return The sum of all elements of matrix
  function lsmat_dense_sum(A)
    implicit none
    type(Lsmatrix), intent(in) :: A
    real(realk) :: lsmat_dense_sum
    integer :: i

    lsmat_dense_sum=0E0_realk
    do i=1,A%nrow*A%ncol
       lsmat_dense_sum=lsmat_dense_sum+a%elms(i)
    enddo

    return
  end function lsmat_dense_sum

end module lsmatrix_operations_dense
