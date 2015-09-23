!> @file
!> Contains csr (compressed sparse row) matrix module.
!> 
!> All functions are derived from their dense counterparts (mat_dense.f90) 
!> and adapted to the Compressed Sparse Row matrix format. The module
!> requires Intel MKL.
!> Rasmus Andersen <raand@chem.au.dk>, June 2010

module matrix_operations_csr
  use matrix_module
  use memory_handling
  use precision
! HACK
  use matrix_operations_dense
! HACK
  real(realk), save :: zeroCSR = 1e-20_realk
!  integer, save :: nnzdiff
! this parameter is used in the integral program (lstensor_operations.f90) and 
! should always be equal to the tol parameter in mat_csr_cleanup (matop_csr_aux.c)
  public :: zeroCSR
  public :: mat_csr_init
  public :: mat_csr_allocate
  public :: mat_csr_copy
  public :: mat_csr_TrAB
  public :: mat_csr_identity
  public :: mat_csr_add_identity
  public :: mat_csr_assign
  public :: mat_csr_dotproduct
  public :: mat_csr_sqnorm2
  public :: mat_csr_outdia_sqnorm2
  public :: mat_csr_abs_max_elm
  public :: mat_csr_max_elm
  public :: mat_csr_trans
  public :: mat_csr_Tr
  public :: mat_csr_scal
  public :: mat_csr_add
  public :: mat_csr_daxpy
  public :: mat_csr_mul
  public :: mat_csr_ao_precond
  public :: mat_csr_set_from_full
  public :: mat_csr_to_full
  public :: mat_csr_to_full3d
  public :: mat_csr_retrieve_block_full
  public :: mat_csr_report_sparsity
  public :: mat_csr_write_to_disk
  public :: mat_csr_read_from_disk
  public :: mat_csr_get_elem
  public :: mat_csr_print
  public :: mat_csr_free
  public :: mat_csr_zero
  public :: mat_csr_inquire_cutoff
  public :: mat_csr_zero_cutoff

  
contains
  !> \brief See mat_init in mat-operations.f90
  
! When initialized, a CSR matrix gets its row array filled with 1's.
! :
  subroutine mat_csr_init(a,nrow,ncol)
    implicit none
    TYPE(Matrix) :: a 
    integer, intent(in) :: nrow, ncol
#ifdef VAR_CSR
    integer(kind=long) :: nsize
    NULLIFY(a%val)
    NULLIFY(a%row)
    NULLIFY(a%col)

    !call mem_alloc(a%row,nrow+1)
    allocate(a%row(nrow+1)) !It's more practical to count this as type(matrix) instead of real and integer /Stinne
    nsize = mem_intsize*size(a%row,kind=long)
    call mem_allocated_mem_type_matrix(nsize)

    a%row = 1 
    a%nrow = nrow
    a%ncol = ncol
    a%nnz = 0 
#endif
  end subroutine mat_csr_init

  subroutine mat_csr_allocate(a,nnz)
    implicit none
    TYPE(Matrix) :: a 
    integer, intent(in) :: nnz
#ifdef VAR_CSR
    integer(kind=long) :: nsize

    if (a%nnz .ne. 0) then
       !call mem_dealloc(a%val)
       nsize = mem_realsize*size(a%val,kind=long)
       deallocate(a%val) !It's more practical to count this as type(matrix) instead of real and integer /Stinne
       call mem_deallocated_mem_type_matrix(nsize)

       !call mem_dealloc(a%col)
       nsize = mem_intsize*size(a%col,kind=long)
       deallocate(a%col) !It's more practical to count this as type(matrix) instead of real and integer /Stinne
       call mem_deallocated_mem_type_matrix(nsize)
    endif
    !call mem_alloc(a%val,nnz)
    allocate(a%val(nnz)) !It's more practical to count this as type(matrix) instead of real and integer /Stinne
    nsize = mem_realsize*size(a%val,kind=long)
    call mem_allocated_mem_type_matrix(nsize)

    !call mem_alloc(a%col,nnz)
    allocate(a%col(nnz)) !It's more practical to count this as type(matrix) instead of real and integer /Stinne
    nsize = mem_intsize*size(a%col,kind=long)
    call mem_allocated_mem_type_matrix(nsize)

    a%nnz = nnz 
#endif
  end subroutine mat_csr_allocate

  subroutine mat_csr_copy(alpha,a,b)
     implicit none
     REAL(REALK),  INTENT(IN)    :: alpha
     TYPE(Matrix), INTENT(IN)    :: a
     TYPE(Matrix), INTENT(INOUT) :: b
#ifdef VAR_CSR
     call mat_csr_assign(b,a)
     if (ABS(alpha-1.0E0_realk).GT.1.0E-14_realk) then
        call mat_csr_scal(alpha,b)
     endif
#endif
  end subroutine mat_csr_copy

  function mat_csr_TrAB(a,b)
     implicit none
     TYPE(Matrix), intent(IN) :: a,b
     REAL(realk) :: mat_csr_trAB 
#ifdef VAR_CSR
     TYPE(Matrix) :: b_temp 
!     real(realk), external :: ddot
!     print *,"in mat_csr nnz's are: ", a%nnz, b%nnz
     if ((a%nnz .eq. 0) .or. (b%nnz .eq. 0)) then
        mat_csr_TrAB = 0
     else
        call mat_csr_init(b_temp, b%nrow, b%ncol)
        call mat_csr_trans(b,b_temp)
        mat_csr_TrAB = mat_csr_dotproduct(a, b_temp)
        call mat_csr_free(b_temp)
     endif
#else
     mat_csr_TrAB = 0.0E0_realk
#endif
  end function mat_csr_TrAB

  subroutine mat_csr_identity(a)
    implicit none
    type(matrix), intent(inout) :: a
#ifdef VAR_CSR
    integer :: i

    ! Zero the matrix
    call mat_csr_zero(a)
    ! Allocate for the diagonal, i.e. nrow entries
    call mat_csr_allocate(a, a%nrow)
    do i = 1,a%nrow
       a%val(i) = 1.0E0_realk
       a%col(i) = i
       a%row(i) = i
    enddo
    ! Set last row element
    a%row(a%nrow+1) = a%nrow+1
#endif
  end subroutine mat_csr_identity
   
  !C = alpha*I + beta*B
  subroutine mat_csr_add_identity(alpha, beta, B, C)
    TYPE(Matrix), intent(IN) :: B
    REAL(realk), INTENT(IN)  :: alpha, beta
    TYPE(Matrix)             :: C
#ifdef VAR_CSR
    TYPE(Matrix)             :: I

    call mat_csr_init(I, B%nrow, B%ncol)
    call mat_csr_identity(I)
    call mat_csr_add(alpha,I,beta,b,c)
    call mat_csr_free(I)
#endif
  end subroutine mat_csr_add_identity

  subroutine mat_csr_assign(a,b)
    implicit none
    TYPE(Matrix), INTENT(INOUT) :: a
    TYPE(Matrix), INTENT(IN)    :: b
#ifdef VAR_CSR
    
    if (a%nnz .ne. 0) then
       call mat_csr_zero(a)
    endif
    if (b%nnz .eq. 0) then
       return
    endif
    call mat_csr_allocate(a,b%nnz)
    a%val = b%val
    a%col = b%col
    a%row = b%row  
#endif
  end subroutine mat_csr_assign

  function mat_csr_dotproduct(a,b)
     implicit none
     TYPE(Matrix), intent(IN) :: a,b
     REAL(realk) :: mat_csr_dotproduct
     real(realk), external :: ddot
     integer     :: i,j,k
 
     mat_csr_dotproduct = 0E0_realk
#ifdef VAR_CSR
     if ((a%nnz .eq. 0) .or. (b%nnz .eq. 0)) then
        return
     endif
     do i=1,a%nrow
        do j=a%row(i),a%row(i+1)-1
           do k=b%row(i),b%row(i+1)-1
              if (a%col(j).eq.b%col(k)) then
                 mat_csr_dotproduct = mat_csr_dotproduct + (a%val(j) * b%val(k))
                 exit
              endif
           enddo
        enddo
     enddo
#endif
  end function mat_csr_dotproduct
  
  function mat_csr_sqnorm2(a)
    implicit none
    TYPE(Matrix), intent(IN) :: a
    REAL(realk) :: mat_csr_sqnorm2
    
    mat_csr_sqnorm2 = mat_csr_dotproduct(a,a)
  end function mat_csr_sqnorm2
  
! Returns sum of the squares of off-diagonal elements in a matrix
  function mat_csr_outdia_sqnorm2(a)
    implicit none
    TYPE(Matrix), intent(IN) :: a
    REAL(realk) :: mat_csr_outdia_sqnorm2
    integer     :: i,j
    
    mat_csr_outdia_sqnorm2 = 0.0E0_realk
#ifdef VAR_CSR
    do i = 1, a%nrow
       do j = a%row(i), a%row(i+1)-1
          if (a%col(j) .ne. i) then
             mat_csr_outdia_sqnorm2 = mat_csr_outdia_sqnorm2 + a%val(j)*a%val(j)
             exit
          endif
       enddo
    enddo
#endif
  end function mat_csr_outdia_sqnorm2

  subroutine mat_csr_abs_max_elm(a,val)  
    implicit none
    type(matrix),intent(in)  :: a
    real(realk), intent(inout) :: val
#ifdef VAR_CSR
    integer                  :: i
    
    if (a%nnz .eq. 0) then 
       val = 0E0_realk
    else
       val = abs(a%val(1))
       do i = 2, a%nnz 
          if (abs(a%val(i)) > val) then
             val = abs(a%val(i))
          endif
       enddo
    endif
#endif
  end subroutine mat_csr_abs_max_elm

  subroutine mat_csr_max_elm(a,val)  
    implicit none
    type(matrix),intent(in)  :: a
    real(realk), intent(inout) :: val
#ifdef VAR_CSR
    integer                  :: i
    
    if (a%nnz .eq. 0) then 
       val = 0E0_realk
    else
       val = a%val(1)
       do i = 2, a%nnz 
          if (a%val(i) > val) then
             val = a%val(i)
          endif
       enddo
    endif
#endif
  end subroutine mat_csr_max_elm

!RA: This one just converts to CSC, consider more optimal solution 
  subroutine mat_csr_trans(a, b)
    implicit none
    TYPE(Matrix), intent(IN) :: a
    TYPE(Matrix), intent(INOUT) :: b 
#ifdef VAR_CSR
    integer :: info
    integer :: job(8)
    job(1) = 0
    job(2) = 1
    job(3) = 1
!RA: Error in mkl doc; job(6) should not be 0 
    job(6) = 1
    
    if (a%nnz .eq. 0) then
       return
    endif
    call mat_csr_allocate(b, a%nnz)
    call mkl_dcsrcsc(job, a%nrow, a%val, a%col, a%row, b%val, b%col,b%row, info)
!RA: Currently (mkl v. 10.2.1.017), info is not used in csrcsc as result status!
!RA: When it comes, check it for the exit status
#endif
  end subroutine mat_csr_trans

  function mat_csr_Tr(a)
    implicit none
    TYPE(Matrix), intent(IN) :: a
    REAL(realk) :: mat_csr_tr 
    integer :: i,j
    
    mat_csr_tr = 0E0_realk
#ifdef VAR_CSR
    do i = 1, a%nrow
       do j = a%row(i), a%row(i+1)-1
          if (a%col(j) == i) then
             mat_csr_tr = mat_csr_tr + a%val(j)
             exit
          endif
       enddo
    enddo
#endif
  end function mat_csr_Tr
  
  subroutine mat_csr_scal(alpha, a)
    !A = alpha * A
    implicit none
    real(realk), intent(in) :: alpha
    type(Matrix), intent(inout) :: a
#ifdef VAR_CSR
    integer :: i
    
    do i=1, a%nnz
       a%val(i) = a%val(i) * alpha
    enddo
#endif
  end subroutine mat_csr_scal

  subroutine mat_csr_add(alpha, a, beta, b, c)
    !C = A + beta*op(B)
    implicit none
    TYPE(Matrix),intent(in)    :: a,b 
    REAL(realk), INTENT(in)    :: alpha,beta
    TYPE(Matrix), intent(inout):: c
#ifdef VAR_CSR
!
    TYPE(Matrix) :: a_temp 
    integer :: i
    !request calculation of nzmax
    integer :: request 
    !sort specifies reordering of the input matrices, none necessary here
    integer :: sort
    !dimension (m=rows of A, n=cols of A)
    integer :: m,n
    !info is the status of the calls to csradd 
    integer :: info
    integer:: idummy
    real(realk) :: ddummy

    info = 0
    sort = 0
    request = 1
    m = a%nrow
    n = a%ncol

    if (a%nnz .eq. 0) then
       call mat_csr_assign(c, b)
       call mat_csr_scal(beta, c)
       return
    elseif (b%nnz .eq. 0) then
       call mat_csr_assign(c, a)
       call mat_csr_scal(alpha, c)
       return
    endif
    ! RA, ugly hack: make a copy of A to avoid messing 
    ! up all "intent(in)" declarations
    call mat_csr_init(a_temp, a%nrow, a%ncol)
    call mat_csr_assign(a_temp,a)

    !scale matrix a_temp
    call mat_csr_scal(alpha, a_temp)

    !zero matrix c if already in use 
    if (c%nnz .ne. 0) then
       call mat_csr_zero(c)
    endif
    !calculate nnz
    call mkl_dcsradd('n', request, sort, m, n, a_temp%val, a_temp%col, a_temp%row, beta, &
         & b%val, b%col, b%row, ddummy, idummy, c%row, idummy, info)
    
    !assert(info!=0)
    !nnz can now be found as the last element in c's row - 1.
    !So use this to allocate for c's val and col arrays.
    call mat_csr_allocate(c, c%row(n+1)-1)
    !alter request to initiate addition
    request = 2
    call mkl_dcsradd('n', request, sort, m, n, a_temp%val, a_temp%col, a_temp%row, beta, &
         & b%val, b%col, b%row, c%val, c%col, c%row, c%nnz, info)
!assert(info!=0)
    call mat_csr_free(a_temp)
#endif
  end subroutine mat_csr_add

!RA: csr_add disallows reuse of a matrix for input and output, 
!so we have to make a temporary copy of y. Consider manual addition
  subroutine mat_csr_daxpy(a,x,y)
    !Y = aX+Y
    implicit none
    real(realk),intent(in)       :: a
    TYPE(Matrix), intent(IN)     :: x
    TYPE(Matrix), intent(INOUT)  :: y
#ifdef VAR_CSR
    TYPE(Matrix) :: y_temp
    real(realk) :: beta
    
    beta = 1.0E0_realk
    call mat_csr_init(y_temp, y%nrow, y%ncol)
    call mat_csr_assign(y_temp, y)
    call mat_csr_add(a, x, beta, y_temp, y) 
    call mat_csr_free(y_temp)
#endif
  end subroutine mat_csr_daxpy

  subroutine mat_csr_mul(a,b, transa, transb, alpha, beta, c)
    !C = op(A) * B
    implicit none
    TYPE(Matrix),intent(in)    :: a,b 
    character, intent(in)      :: transa, transb
    REAL(realk), INTENT(in)    :: alpha, beta
    TYPE(Matrix), intent(inout):: c
    TYPE(Matrix) :: b_temp, c_temp
#ifdef VAR_CSR
!
    integer :: i,j,nnz
    !request calculation of nnz
    integer :: request
    !sort specifies reordering of the input matrices, none necessary here
    integer :: sort 
    !dimension (m=rows of A, n=cols of A, k=cols of B)
    integer :: m,n,k 
    !info is the status of the calls to csrmultcsr 
    integer :: info
    !dummy variables for first call to csrmultcsr
    REAL(realk) :: ddummy
    integer :: idummy 
    
    idummy = 0
    ddummy = 0E0_realk
    request = 1
    sort = 0
    m = a%nrow
    n = a%ncol
    k = b%ncol

    !check for zero-matrices and just return a zero-matrix
    if ((a%nnz .eq. 0) .or. (b%nnz .eq. 0)) then
       !Thomas K: Hack to allow for non zero beta
       IF(ABS(beta).LT.1.0E-25_realk)THEN
          if (c%nnz .ne. 0) then
             call mat_csr_zero(c)
          endif
       ELSE
          call mat_csr_scal(beta,c)
       ENDIF
       return
    endif

    !RA: hack to allow transb, consider better solution
    call mat_csr_init(b_temp, b%nrow, b%ncol)
    if (transb == 't' .or. transb == 'T') then
       call mat_csr_allocate(b_temp, b%nnz)
       call mat_csr_trans(b,b_temp)
    else
       call mat_csr_assign(b_temp, b)
    endif

    IF(ABS(beta).GT.1.0E-25_realk)THEN
       !Thomas K: Hack to allow for non zero beta
       call mat_csr_init(c_temp, c%nrow, c%ncol)
       call mat_csr_assign(c_temp, c)
    ENDIF

    !zero matrix c if already in use 
    if (c%nnz .ne. 0) then
       call mat_csr_zero(c)
    endif
    !calculate nnz
    call mkl_dcsrmultcsr(transa, request, sort, m, n, k, a%val, a%col, a%row, &
         & b_temp%val, b_temp%col, b_temp%row, ddummy, idummy, c%row, idummy, info)
!RA: assert(info!=0)
    !nnz can now be found as the last element in c's row - 1
    nnz = c%row(n+1)-1
    !allocate for c's column and values arrays
    call mat_csr_allocate(c, nnz)
    !alter request to initiate multiplication
    request = 2

    call mkl_dcsrmultcsr(transa, request, sort, m, n, k, a%val, a%col, a%row, &
         & b_temp%val, b_temp%col, b_temp%row, c%val, c%col, c%row, c%nnz, info)
    !Stinne fix: multiply by alpha if different from one!
    if (abs(1.0E0_realk-alpha) > 1.0E-7_realk) then
       call mat_csr_scal(alpha,c)
    endif
    !Thomas K: Hack to allow for non zero beta
    IF(ABS(beta).GT.1.0E-25_realk)THEN
       call mat_csr_daxpy(beta,c_temp,c)
       call mat_csr_free(c_temp)
    ENDIF

!RA: assert(info!=0)
    !print *, "Will clean result matrix, nnz is ", c%nnz
    call mat_csr_cleanup(c%val, c%col, c%row, n, c%nnz, zeroCSR)
    !nnzdiff = nnzdiff + c%nnz-nnz
    !RA ugly hack: for the internal memory management, 
    !register how many elements were removed from the 
    !col and val arrays during the cleanup phase
!    call mem_deallocated_mem_integer((nnz-c%nnz)*mem_intsize)
!    call mem_deallocated_mem_real((nnz-c%nnz)*mem_realsize)
    !print *, "Done cleaning result matrix in matrix mult, nnz is now", c%nnz, ", result matrix:"
    ! RA: free b_temp hack
    call mat_csr_free(b_temp)
#endif
  end subroutine mat_csr_mul
  
  subroutine mat_csr_ao_precond(symmetry,omega,FUP,FUQ,DU,X_AO)
    implicit none
    integer, intent(in) :: symmetry
    real(realk), intent(in) :: omega
    type(Matrix), intent(in) :: FUP, FUQ, DU
    type(Matrix), intent(inout) :: X_AO
#ifdef VAR_CSR
    real(realk) :: denom, err, fup_i, fup_j, fuq_i, fuq_j, du_i, du_j
    integer :: ndim,i,j,k

    do i = 1, X_AO%nrow !row i
       do k = X_AO%row(i), X_AO%row(i+1)-1 
          j = X_AO%col(k) !column j
          denom = 1
          fuq_j = mat_csr_get_elem(FUQ, j,j)
          fuq_i = mat_csr_get_elem(FUQ, i,i)
          fup_j = mat_csr_get_elem(FUP, j,j)
          fup_i = mat_csr_get_elem(FUP, i,i)
          if (symmetry == 1 .or. symmetry == 2) then 
             !Symmetric or antisymmetric X_AO
             denom = fuq_j - fup_j + fuq_i - fup_i - omega
             if (ABS(denom) > 1.0E-10_realk) then
                X_AO%val(k) = X_AO%val(k)/denom 
             endif
          else 
             !X_AO not symmetric in any way
             du_j = mat_csr_get_elem(DU, j,j)
             du_i = mat_csr_get_elem(DU, i,i)
             denom = fuq_j - fup_j + fuq_i - fup_i - omega*(du_j - du_i)
             if (ABS(denom) > 1.0E-10_realk) then
                X_AO%val(k) = X_AO%val(k)/denom 
             endif
          endif
       enddo
    enddo
#endif
  end subroutine mat_csr_ao_precond

  subroutine mat_csr_set_from_full(afull,alpha,a)
    implicit none
    real(realk), INTENT(IN) :: afull(*)
    real(realk), intent(in) :: alpha
    TYPE(Matrix)            :: a 
#ifdef VAR_CSR
    INTEGER       job(8)
    INTEGER       m, n, lda, info, nnz
    !dummy variables for first call to ddnscsr
    REAL(realk) :: ddummy
    integer :: idummy     
    !    call mat_csr_init(a,afull%nrow, afull%ncol)
    job(1) = 0
    job(2) = 1
    job(3) = 1
    job(4) = 2
    job(5) = 0
    job(6) = 0
    lda = a%nrow
    m = a%nrow
    n = a%ncol
    call mkl_ddnscsr(job, m, n, afull, lda, ddummy, idummy, a%row, info)
    !nnz can now be found as the last element -1 in the row array
    nnz = a%row(n+1)-1
    job(5) = nnz
    !the column and values arrays can now be allocated 
    if (a%nnz .ne. 0) then
       call mat_csr_zero(a)
    endif
    call mat_csr_allocate(a, nnz)
    !alter request to initiate actual conversion
    job(6) = 1
    call mkl_ddnscsr(job, m, n, afull, lda, a%val, a%col, a%row, info)
     ! assert info     
     ! RA: maybe we should clean the matrix here?
#endif
  end subroutine mat_csr_set_from_full

!> \brief See mat_to_full in mat-operations.f90
  subroutine mat_csr_to_full(a, alpha, afull)
     implicit none
     TYPE(Matrix), intent(in) :: a
     real(realk), intent(in) :: alpha
     real(realk), intent(inout):: afull(*)  
#ifdef VAR_CSR
     integer                 :: i

     INTEGER       job(8)
     INTEGER       m, n, lda, info
     afull(1) = 0.0E0_realk !Stinne fix - some compilers complain about intent(out)
                      !without MKL present
     job(1) = 1
     job(2) = 1
     job(3) = 1
     job(4) = 2
     job(5) = a%nnz
     job(6) = 2
     lda = a%nrow
     m = a%nrow
     n = a%ncol
     !RA: hack to handle zero matrix
     if (a%nnz .eq. 0) then
        return
     endif
     call mkl_ddnscsr(job, m, n, afull, lda, a%val, a%col, a%row, info)
    if (abs(1.0E0_realk-alpha) > 1.0E-7_realk) then
       do i=1, m*n
          afull(i) = afull(i) * alpha
       enddo
    endif
#endif
  end subroutine mat_csr_to_full

!> \brief See mat_to_full in mat-operations.f90
  subroutine mat_csr_to_full3D(a, alpha, afull,n1,n2,n3,i3)
     implicit none
     integer, INTENT(IN)           :: n1,n2,n3,i3
     TYPE(Matrix), intent(in) :: a
     real(realk), intent(in) :: alpha
     real(realk), intent(inout):: afull(n1,n2,n3)  
#ifdef VAR_CSR
     INTEGER      ::  job(8)
     INTEGER      ::  m, n, lda, info,j,i,offset,mp1
     real(realk),pointer  :: Afulltmp(:)
     call mem_alloc(Afulltmp,n1*n2)
     job(1) = 1
     job(2) = 1
     job(3) = 1
     job(4) = 2
     job(5) = a%nnz
     job(6) = 2
     lda = a%nrow
     m = a%nrow
     n = a%ncol
     !RA: hack to handle zero matrix
     if (a%nnz .eq. 0) then
        return
     endif
     call mkl_ddnscsr(job, m, n, afulltmp, lda, a%val, a%col, a%row, info)
     N = a%nrow    !change diff
     M = MOD(N,7)  !change diff
     IF (M.NE.0) THEN
        do j = 1,a%ncol
           offset = (j-1)*N
           DO I = 1,M
              afull(i,j,i3) = alpha*afulltmp(i+offset)              
           ENDDO
        enddo
        MP1 = M + 1
        IF (N.GE.7)THEN
           do j = 1,a%ncol
              offset = (j-1)*N
              DO I = MP1,N,7
                 afull(i,j,i3) = alpha*afulltmp(i+offset)
                 afull(i+1,j,i3) = alpha*afulltmp(i+1+offset)
                 afull(i+2,j,i3) = alpha*afulltmp(i+2+offset)
                 afull(i+3,j,i3) = alpha*afulltmp(i+3+offset)
                 afull(i+4,j,i3) = alpha*afulltmp(i+4+offset)
                 afull(i+5,j,i3) = alpha*afulltmp(i+5+offset)
                 afull(i+6,j,i3) = alpha*afulltmp(i+6+offset)
              END DO
           enddo
        ENDIF
     ELSE
        do j = 1,a%ncol
           offset = (j-1)*N
           DO I = 1,N,7
              afull(i,j,i3) = alpha*afulltmp(i+offset)
              afull(i+1,j,i3) = alpha*afulltmp(i+1+offset)
              afull(i+2,j,i3) = alpha*afulltmp(i+2+offset)
              afull(i+3,j,i3) = alpha*afulltmp(i+3+offset)
              afull(i+4,j,i3) = alpha*afulltmp(i+4+offset)
              afull(i+5,j,i3) = alpha*afulltmp(i+5+offset)
              afull(i+6,j,i3) = alpha*afulltmp(i+6+offset)
           END DO
        ENDDO
     ENDIF
     call mem_dealloc(Afulltmp)
#endif
  end subroutine mat_csr_to_full3D

  subroutine mat_csr_retrieve_block_full(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in)         :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(inout)     :: fullmat(fullrow,fullcol)
    type(Matrix), intent(in) :: A
#ifdef VAR_CSR
    integer                     :: i, j, k
    
    fullmat = 0E0_realk
    do i = insertrow, insertrow+fullrow-1   !rows
       do k = a%row(i), a%row(i+1)-1 
          j = a%col(k) !column indexes for this row
          if ((j .ge. insertcol) .and. (j .lt. insertcol+fullcol)) then
             fullmat(j-insertcol+1, i-insertrow+1) = a%val(k)
          endif
       enddo
    enddo
#endif
  end subroutine mat_csr_retrieve_block_full

! RA: todo: implement this one...
  subroutine mat_csr_retrieve_block_csr(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in)         :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(inout)     :: fullmat(fullrow,fullcol)
    type(Matrix), intent(in) :: A
#ifdef VAR_CSR
    integer                     :: i, j, k
    
    fullmat = 0E0_realk
    do i = fullrow, fullrow+insertrow-1   !rows
       do k = a%row(i), a%row(i+1)-1 
          j = a%col(k) !column indexes for this row
!          print *,"searching col ", j, " in row ",i 
          if ((j .ge. fullcol) .and. (j.lt.fullcol+insertcol)) then
!             fullmat(i-fullrow+1, j-fullcol+1) = a%val(k)
!             print *, "using val: ", i-fullrow+1, j-fullcol+1, a%val(k)
          endif
       enddo
    enddo
!    do i=a%row(fullrow), a%row(fullrow+insertrow)
!       do j=a%col()
    !do j = insertcol, insertcol+fullcol-1
       !do i = insertrow, insertrow+fullrow-1
!          fullmat(i-insertrow+1,j-insertcol+1) = A%elms((j-1)*A%nrow+i)
#endif
  end subroutine mat_csr_retrieve_block_csr
  
  subroutine mat_csr_report_sparsity(A,sparsity)
    implicit none
    type(Matrix) :: A
    real(realk)  :: sparsity
#ifdef VAR_CSR
    integer      :: nnz
    
    if (A%nnz == 0) then
       ! matrix not associated
       sparsity = 0.0E0_realk
    else
       nnz = a%nnz
       call mat_csr_cleanup(A%val, A%col, A%row, A%nrow, A%nnz, zeroCSR)
       !if (nnz /= a%nnz) then
       !   call mem_deallocated_mem_integer((nnz-a%nnz)*mem_intsize)
       !   call mem_deallocated_mem_real((nnz-a%nnz)*mem_realsize)
       !endif

       sparsity = A%nnz/(A%ncol*A%nrow)
    endif
#endif
  end subroutine mat_csr_report_sparsity

  subroutine mat_csr_write_to_disk(iunit,A)
    implicit none
    integer, intent(in) :: iunit
    type(Matrix), intent(in) :: A
#ifdef VAR_CSR
    integer :: i
    integer(kind=long) :: ncol,nrow,nnz 

    nnz = A%nnz
    nrow = A%Nrow
    ncol = A%Ncol
    WRITE(iunit) nnz, nrow, ncol
    WRITE(iunit) (A%val(I),I=1,A%nnz)
    WRITE(iunit) (A%col(I),I=1,A%nnz)
    WRITE(iunit) (A%row(I),I=1,A%nrow+1)
#endif
  end subroutine mat_csr_write_to_disk

  subroutine mat_csr_read_from_disk(iunit,A)
    implicit none
    integer, intent(in) :: iunit
    type(Matrix), intent(inout) :: A
#ifdef VAR_CSR
    integer :: i
    integer(kind=long) :: ncol,nrow,nnz 

    READ(iunit) nnz, nrow, ncol
    IF(Nrow.EQ.A%Nrow.AND.Ncol.EQ.A%Ncol)THEN
       if (nnz .ne. A%nnz) then
          call mat_csr_zero(A)
          call mat_csr_allocate(A, nnz)
       endif
       READ(iunit) (A%val(I),I=1,A%nnz)
       READ(iunit) (A%col(I),I=1,A%nnz)
       READ(iunit) (A%row(I),I=1,A%nrow+1)
    ELSE
       print*,'Error in reading CSR matrix from disk. Dimension mismatch'
       print*,'Dimensions of the CSR matrix on Disk  :',Nrow,Ncol
       print*,'Allocated Dimensions of the CSR matrix:',A%Nrow,A%Ncol
       CALL LSQUIT('Error in reading CSR matrix from disk. Dimension mismatch',-1)
    ENDIF
#endif
  end subroutine mat_csr_read_from_disk

  function mat_csr_get_elem(a, i, j)
     TYPE(Matrix), intent(IN) :: a
     INTEGER, intent(IN)  :: i,j 
#ifdef VAR_CSR
     REAL(realk) :: mat_csr_get_elem
     integer column

     mat_csr_get_elem = 0.0
     ! loop through all column indexes belonging to row i
     do column=a%row(i), a%row(i+1)-1
        if (a%col(column) .eq. j) then
           mat_csr_get_elem = a%val(column)
        endif
     enddo
#else
     mat_csr_get_elem = 0.0
#endif
  end function mat_csr_get_elem

  !> \brief See mat_print in mat-operations.f90
  subroutine mat_csr_print(a, lu)
    implicit none
    TYPE(Matrix), intent(in) :: a
    integer, intent(in) :: lu 
    real(realk), allocatable :: afull(:,:)
#ifdef VAR_CSR
    real(realk) :: alpha
!    integer :: i,j,k,cur_row,cur_col, cur_val
!    WRITE(lu,'(2X,A4,5E13.3,/(6X,5E13.3))')'VAL:',(a%val(j),j=1,a%nnz)
!    WRITE(lu,'(2X,A4,15I4,/(6X,15I4))')'COL:',(a%col(j),j=1,a%nnz)
!    WRITE(lu,'(2X,A4,15I4,/(6X,15I4))')'ROW:',(a%row(j),j=1,a%nrow+1)
! Hack to handle zero-matrices
!    if (a%nnz .eq. 0) then
!       call mem_alloc(a%val, 1)
!       call mem_alloc(a%col, 1)
!       a%val = 0
!       a%col = 0
!    endif
!    call mat_csr_pretty_print(a%val, a%col, a%row, a%nrow)
!    fd = fnum(lu)
!    call mat_csr_column_print(fd, a%val, a%col, a%row, a%nrow)
    allocate(afull(a%nrow,a%ncol))
    afull = 0.0E0_realk
    alpha = 1.0E0_realk
    call mat_csr_to_full(a, alpha, afull)
    call LS_OUTPUT(afull, 1, a%nrow, 1, a%ncol, a%nrow, a%ncol, 1, lu)
    !call LS_OUTPUT(afull, 1, a%nrow, 1, a%ncol, a%nrow, a%ncol, 1, 6)
    deallocate(afull)
#endif
  end subroutine mat_csr_print
  
  !> \brief See mat_free in mat-operations.f90
  subroutine mat_csr_free(a)
    implicit none
    TYPE(Matrix) :: a 
#ifdef VAR_CSR
    integer(kind=long) :: nsize

    if (a%nnz .ne. 0) then
       !call mem_dealloc(a%val)
       nsize = mem_realsize*size(a%val,kind=long)
       deallocate(a%val) !It's more practical to count this as type(matrix) instead of real and integer /Stinne
       call mem_deallocated_mem_type_matrix(nsize)

       !call mem_dealloc(a%col)
       nsize = mem_intsize*size(a%col,kind=long)
       deallocate(a%col) !It's more practical to count this as type(matrix) instead of real and integer /Stinne
       call mem_deallocated_mem_type_matrix(nsize)

       a%nnz = 0
    endif
    if (associated(a%row)) then
       !call mem_dealloc(a%row)
       nsize = mem_intsize*size(a%row,kind=long)
       deallocate(a%row) !It's more practical to count this as type(matrix) instead of real and integer /Stinne
       call mem_deallocated_mem_type_matrix(nsize)
    endif
    a%nrow = 0
    a%ncol = 0
#endif
  end subroutine mat_csr_free

  
! A zero matrix in CSR is expressed as:
! val = []
! col = []
! row = [1,1,1,..., 1] ; nrow+1 entries
! nnz = 0
! nrow = n
! ncol = n
  subroutine mat_csr_zero(a)
    implicit none
    TYPE(Matrix) :: a 
#ifdef VAR_CSR
    integer :: i
    integer(kind=long) :: nsize

    if (a%nnz .ne. 0) then
       !call mem_dealloc(a%val)
       nsize = mem_realsize*size(a%val,kind=long)
       deallocate(a%val) !It's more practical to count this as type(matrix) instead of real and integer /Stinne
       call mem_deallocated_mem_type_matrix(nsize)

       !call mem_dealloc(a%col)
       nsize = mem_intsize*size(a%col,kind=long)
       deallocate(a%col) !It's more practical to count this as type(matrix) instead of real and integer /Stinne
       call mem_deallocated_mem_type_matrix(nsize)

       NULLIFY(a%val)
       NULLIFY(a%col)
    endif
    do i=1, a%nrow+1
       a%row(i) = 1
    enddo
    a%nnz = 0
    !call mat_csr_cleanup(a%val, a%col, a%row, a%nrow, a%nnz, zeroCSR)
#endif
  end subroutine mat_csr_zero

!> See mat_inquire_cutoff in mat-operations.f90
  subroutine mat_csr_inquire_cutoff(cutoff)
    implicit none
    real(realk), intent(inout) :: cutoff

    cutoff = zeroCSR

  end subroutine mat_csr_inquire_cutoff

!> See mat_zero_cutoff in mat-operations.f90
  subroutine mat_csr_zero_cutoff(cutoff)
    implicit none
    real(realk), intent(in) :: cutoff
    zeroCSR = cutoff
  end subroutine mat_csr_zero_cutoff

end module matrix_operations_csr
