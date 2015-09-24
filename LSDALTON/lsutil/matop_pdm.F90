!> @file
!> Contains parallel distributed memory (PDM) matrix module.
!> 
!> All functions are derived from their dense counterparts (mat_dense.f90) 
!> and adapted to the PDM matrix format. 
!>  
!> The thinking is 
!> 1. Memory is only distributed on the slaves. The master does not contain 
!>    any part of the matrices. This is to reduce memory on master, which most likely
!>    have other things.  
!> 2. When doing LAPACK calls the main matrices are collected on node pdmm_collector. 
!>    node pdmm_collector, is special in the sense that it is used to collect the full 
!>    matrix. Perform LAPACK and redestribute the resulting matrix. Normally this is 
!>    done on rank 1 (and on the master rank in other codes) to ensure reduced memory usage 
!>    on master
!> 3. The matrices are not necesarrily distributed across all nodes, this depend on memory
!>    matrices are not necessarily distributed in the same manner 
!>      Mat A could be distributed like (  1   2  )
!>                                      (  3   4  )
!>      Mat B could be distributed like (  5   6  )
!>                                      (  7   8  )
!>      Mat C could be distributed like (  1   2  )
!>                                      (  3   4  )
!>      so A + C is easy A + B 
!>      For now use same. for all. 
!>    In general the focus is to distribute the mem not on the parallelization.
!>    so if all mats are 5 GB on option is to have full matrices on different nodes.
!>    the CPU LAPACK calls would be fast - at least make sure that N is big enough to 
!>    hide memory latency - so around N=2000 but preferably bigger.                                     
!> 4. Use LSGEMM wrapper to use GPU for large matrices 
!> 5.  
!>  
module matrix_operations_pdmm
#ifdef VAR_ENABLE_TENSORS
  use lsparameters
  use memory_handling
  use matrix_module
  use LSmatrix_type
  use precision
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
  use lsmpi_param
#endif
  use tensor_interface_module

  !> Pointer to a type(tensor). Necessary if we want arrays of derived types!
  type MatrixPDMp
     type(tensor), pointer :: p
  end type MatrixPDMp

  integer(kind=ls_mpik) :: pdmm_master
  integer(kind=ls_mpik) :: pdmm_nodtot
  integer(kind=ls_mpik) :: pdmm_mynum
  integer(kind=ls_mpik) :: pdmm_comm
  integer(kind=ls_mpik) :: pdmm_collector
  logical :: pdmm_member
  save
  Type(MatrixPDMp) :: pdmm_Marray(500)
  integer :: pdmm_MaxAllocatedMat
  integer :: BLOCK_SIZE_PDM
!  integer(kind=8) :: iwrkTensor
!  real(realk),pointer :: wrkTensor(:)
  logical :: pdmm_mpi_set
contains
  subroutine InitMatrixModulePDM()
    implicit none
    pdmm_MaxAllocatedMat = 0
    call MARRAY_NULLIFY
!    iwrkTensor = 4*BLOCK_SIZE_PDM*BLOCK_SIZE_PDM
!    call mem_alloc(wrkTensor,iwrkTensor)
  end subroutine InitMatrixModulePDM

  subroutine FreeMatrixModulePDM()
    implicit none
    integer :: I
    MarrayList: DO I=1,size(pdmm_Marray)
       IF(ASSOCIATED(pdmm_Marray(I)%p))THEN
          call lsquit('Memory Leak in PDM',-1)
       ENDIF
    ENDDO MarrayList
!    IF(iwrkTensor.EQ.-1)CALL LSQUIT('iwrkTensor error in PDMM',-1)
!    call mem_dealloc(wrkTensor)
!    iwrkTensor = -1 
  end subroutine FreeMatrixModulePDM

  SUBROUTINE MARRAY_NULLIFY
    IMPLICIT NONE
    INTEGER(KIND=LS_MPIK) :: I
    DO I=1, size(pdmm_Marray)
       NULLIFY(pdmm_Marray(I)%p)
    ENDDO
  END SUBROUTINE MARRAY_NULLIFY
  
  !> \brief See mat_init in mat-operations.f90
  subroutine mat_pdmm_init(A,nrow,ncol)
    implicit none
    TYPE(Matrix) :: A
    integer, intent(in) :: nrow, ncol
    integer :: i,j,k,ID
    integer(kind=long) :: nsizeFULL,nsizeLOCAL
    integer :: mode,offsetpdmm
    A%nrow = nrow
    A%ncol = ncol
    mode = 2
    !Step 1: Put into list - find PDMID  
    MarrayList: DO I=1,pdmm_MaxAllocatedMat+1
       IF(.NOT.ASSOCIATED(pdmm_Marray(I)%p))THEN
          Allocate(pdmm_Marray(I)%p)
          ID = I          
          pdmm_MaxAllocatedMat = MAX(pdmm_MaxAllocatedMat,I)
          Exit MarrayList
       ENDIF
    ENDDO MarrayList
    A%PDMID = ID
    !minit  master calls only while slaves are asleep
    offsetpdmm = 0
#ifdef VAR_MPI
    IF(infpar%nodtot.GT.1)offsetpdmm=1
#else
    call lsquit('PDMM requires MPI',-1)
#endif
    call tensor_minit(pdmm_Marray(ID)%p,[nrow,ncol],mode,tdims=[BLOCK_SIZE_PDM,BLOCK_SIZE_PDM],atype="TDAR",fo=offsetpdmm)

  end subroutine mat_pdmm_init

  !> \brief See mat_free in mat-operations.f90
  subroutine mat_pdmm_free(A)
    implicit none
    TYPE(Matrix) :: A
    call tensor_free(pdmm_Marray(A%PDMID)%p)
    deallocate(pdmm_Marray(A%PDMID)%p)
    nullify(pdmm_Marray(A%PDMID)%p)
   end subroutine mat_pdmm_free

   !> \brief See mat_set_from_full in mat-operations.f90
   subroutine mat_pdmm_set_from_full(afull,alpha,a)
     implicit none
     TYPE(Matrix)               :: a 
     real(realk), INTENT(IN)    :: afull(A%nrow,A%ncol)
     real(realk), intent(in)    :: alpha
     
     !  convert from arg1 to arg2 
     call tensor_convert(Afull,pdmm_Marray(A%PDMID)%p)!,&
!          & wrk=wrkTensor,iwrk=iwrkTensor)
     if (ABS(ALPHA-1.0E0_realk).GT.1.0E-15_realk)THEN
        call mat_pdmm_scal(alpha,A)
     ENDIF
   end subroutine mat_pdmm_set_from_full

!> \brief See mat_to_full in mat-operations.f90
   subroutine mat_pdmm_to_full(a,alpha,afull)
     implicit none
     TYPE(Matrix),intent(in) :: a 
     real(realk), intent(in) :: alpha
     real(realk)  :: afull(A%nrow,A%ncol)     
     call tensor_convert(pdmm_Marray(A%PDMID)%p,Afull)!,wrk=wrkTensor,iwrk=iwrkTensor)
     if (ABS(ALPHA-1.0E0_realk).GT.1.0E-15_realk)THEN
        call dscal(A%nrow*A%ncol,alpha,afull,1)
     endif
   end subroutine mat_pdmm_to_full

   !> \brief See mat_print in mat-operations.f90
   subroutine mat_pdmm_print(a, i_row1, i_rown, j_col1, j_coln, lu)
     implicit none
     TYPE(Matrix),intent(in) :: a
     integer, intent(in)     :: i_row1, i_rown, j_col1, j_coln, lu
     real(realk) :: alpha
     real(realk),pointer :: Afull(:,:)
     alpha=1.0E0_realk
     call mem_alloc(Afull,A%nrow,A%ncol)
     call mat_pdmm_to_full(A,alpha,Afull)
     call LS_OUTPUT(Afull, i_row1, i_rown, j_col1, j_coln, a%nrow, a%ncol, 1, lu)
     call mem_dealloc(Afull)
   end subroutine mat_pdmm_print

  !> \brief See mat_trans in mat-operations.f90
  subroutine mat_pdmm_trans(A,B) 
    implicit none
    TYPE(Matrix) :: A,B !B=A^T
    integer :: order(2)
    order(1) = 2
    order(2) = 1
    call tensor_cp_data(pdmm_Marray(A%PDMID)%p,pdmm_Marray(B%PDMID)%p,order)
    
  end subroutine mat_pdmm_trans

  !> \brief See mat_dpotrf in mat-operations.f90
  SUBROUTINE mat_pdmm_dpotrf(a) 
    implicit none
    TYPE(Matrix),intent(inout)     :: a !overwrites matrix A
    integer :: info
    real(realk) :: alpha
    real(realk),pointer :: Afull(:,:)
    INFO = 0
    alpha=1.0E0_realk
    call mem_alloc(Afull,A%nrow,A%ncol)
    call mat_pdmm_to_full(A,alpha,Afull)
    CALL DPOTRF('U',A%nrow,Afull,A%nrow,INFO)
    IF(INFO.NE.0)CALL LSQUIT('DPOTRF ERROR',-1)
    call mat_pdmm_set_from_full(afull,alpha,A)
    call mem_dealloc(Afull)
  end SUBROUTINE mat_pdmm_dpotrf

  !> \brief See mat_dpotrs in mat-operations.f90
  SUBROUTINE mat_pdmm_dpotrs(A,B) 
    implicit none
    TYPE(Matrix),intent(in)    :: a
    TYPE(Matrix),intent(inout) :: b
    integer :: info
    real(realk) :: alpha
    real(realk),pointer :: Afull(:,:),Bfull(:,:)
    INFO = 0
    alpha=1.0E0_realk
    call mem_alloc(Afull,A%nrow,A%ncol)
    call mat_pdmm_to_full(A,alpha,Afull)
    call mem_alloc(Bfull,B%nrow,B%ncol)
    call mat_pdmm_to_full(B,alpha,Bfull)
    CALL DPOTRS('U',A%nrow,B%ncol,Afull,A%nrow,Bfull,B%nrow,INFO)
    call mem_dealloc(Afull)
    call mat_pdmm_set_from_full(Bfull,alpha,B)
    call mem_dealloc(Bfull)
  end SUBROUTINE mat_pdmm_dpotrs

  !> \brief See mat_dpotri in mat-operations.f90
  SUBROUTINE mat_pdmm_dpotri(A) 
    implicit none
    TYPE(Matrix),intent(inout)     :: a
    integer :: info
    real(realk) :: alpha
    real(realk),pointer :: Afull(:,:)
    INFO = 0
    alpha=1.0E0_realk
    call mem_alloc(Afull,A%nrow,A%ncol)
    call mat_pdmm_to_full(A,alpha,Afull)
    CALL DPOTRI('U',A%nrow,Afull,A%nrow,INFO)           
    IF(INFO.NE.0)CALL LSQUIT('DPOTRI ERROR',-1)
    call mat_pdmm_set_from_full(Afull,alpha,A)
    call mem_dealloc(Afull)
  end SUBROUTINE mat_pdmm_dpotri

  !> \brief  See mat_inv in mat-operations.f90
  SUBROUTINE mat_pdmm_inv(A, A_inv) 
    implicit none
    TYPE(Matrix),intent(in)     :: A
    TYPE(Matrix)                :: A_inv !output
    real(realk), pointer   :: work1(:)
    integer,pointer    :: IPVT(:)
    real(realk)            :: RCOND, dummy(2), tmstart, tmend
    integer                :: IERR, i, j, fulldim, ndim
    real(realk),pointer :: A_inv_full(:,:)
    fulldim = a%nrow
    call mem_alloc(A_inv_full,fulldim,fulldim) 
    call mem_alloc(work1,fulldim)
    call mem_alloc(IPVT,fulldim)
    !Invert U and Ut:
    IPVT = 0 ; RCOND = 0.0E0_realk  
    call mat_pdmm_to_full(A,1.0E0_realk,A_inv_full)
    call DGECO(A_inv_full,fulldim,fulldim,IPVT,RCOND,work1)
    call DGEDI(A_inv_full,fulldim,fulldim,IPVT,dummy,work1,01)
    !Convert framework:
    call mat_pdmm_set_from_full(A_inv_full,1.0E0_realk,A_inv)
    call mem_dealloc(A_inv_full) 
    call mem_dealloc(work1)
    call mem_dealloc(IPVT)
  END SUBROUTINE mat_pdmm_inv

  !> \brief See mat_assign in mat-operations.f90
  !B = A
  subroutine mat_pdmm_assign(B,A) 
    implicit none
    TYPE(Matrix) :: A,B     
    call tensor_cp_data(pdmm_Marray(A%PDMID)%p,pdmm_Marray(B%PDMID)%p)
  end subroutine mat_pdmm_assign
  
  subroutine mat_pdmm_copy(alpha,A,B) 
    implicit none
    TYPE(Matrix) :: A,B
    real(realk) :: Alpha
    call mat_pdmm_assign(B,A)
    if (ABS(alpha-1.0E0_realk).GT.1.0E-15_realk)call mat_pdmm_scal(alpha,b)
  end subroutine mat_pdmm_copy

  !> \brief See mat_tr in mat-operations.f90
  function mat_pdmm_tr(A) 
    implicit none
    TYPE(Matrix) :: A !B=A^T
    real(realk) :: tracefunc,mat_pdmm_tr,alpha
    real(realk),pointer :: Afull(:,:)
    integer :: i
    alpha=1.0E0_realk
    call mem_alloc(Afull,A%nrow,A%ncol)
    call mat_pdmm_to_full(A,alpha,Afull)
    mat_pdmm_tr = 0E0_realk
    do i = 1,a%nrow
       mat_pdmm_tr = mat_pdmm_tr + Afull(i,i)
    enddo
    call mem_dealloc(Afull)
  end function mat_pdmm_tr
  
  function mat_pdmm_trAB(a,b)
    implicit none
    TYPE(Matrix), intent(IN) :: a,b
    REAL(realk) :: mat_pdmm_trAB
    !
    integer :: i,j
    real(realk),pointer :: Afull(:,:),Bfull(:,:)
    call mem_alloc(Afull,A%nrow,A%ncol)
    call mat_pdmm_to_full(A,1.0E0_realk,Afull)
    call mem_alloc(Bfull,B%nrow,B%ncol)
    call mat_pdmm_to_full(B,1.0E0_realk,Bfull)
    mat_pdmm_TrAB = 0.0E0_realk
    do j = 1,a%ncol
       do i = 1,a%nrow
          mat_pdmm_TrAB = mat_pdmm_TrAB + Afull(i,j)*Bfull(j,i)
       enddo
    enddo    
    call mem_dealloc(Afull)
    call mem_dealloc(Bfull)
  end function mat_pdmm_trAB

  !> \brief See mat_zero in mat-operations.f90
  subroutine mat_pdmm_mul(a, b, ta, tb, alpha, beta, c)
    implicit none
    TYPE(Matrix) :: A,B,C
    CHARACTER    :: TA,TB
    REAL(REALK)  :: Alpha,Beta    
    !local
    !PDMM fallback
    type(matrix) :: TMP
    integer :: M,N,K,m2cA(1),m2cB(1),nmodes2c,order(2)
    real(realk),pointer :: Afull(:,:),Bfull(:,:),Cfull(:,:)    
!!$    call mem_alloc(Afull,A%nrow,A%ncol)
!!$    call mat_pdmm_to_full(A,1.0E0_realk,Afull)
!!$    call mem_alloc(Bfull,B%nrow,B%ncol)
!!$    call mat_pdmm_to_full(B,1.0E0_realk,Bfull)
!!$    call mem_alloc(Cfull,C%nrow,C%ncol)
!!$    IF(ABS(beta).GT.1.0E-14_realk)THEN
!!$       call mat_pdmm_to_full(C,1.0E0_realk,Cfull)
!!$    ENDIF
!!$    m = C%nrow
!!$    n = C%ncol
!!$    if (ta == 'n' .or. ta == 'N') then
!!$       k = a%ncol
!!$    elseif(ta == 't' .or. ta == 'T') then
!!$       k = a%nrow
!!$    else
!!$       print*,'unknown format in mat_pdmm_mul'
!!$       CALL LSQUIT( 'unknown format in mat_pdmm_mul',-1)
!!$    endif
!!$    call DGEMM(ta,tb,m,n,k,alpha,&
!!$         &afull,a%nrow,bfull,b%nrow,beta,cfull,c%nrow)
!!$    call mat_pdmm_set_from_full(Cfull,1.0E0_realk,C)
!!$    call mem_dealloc(Cfull)
!!$    call mem_dealloc(Bfull)
!!$    call mem_dealloc(Afull)

    if (ta == 'n' .or. ta == 'N') then
       m2cA(1) = 2 !mode on A to contract 
    else
       m2cA(1) = 1 !mode on A to contract 
    endif
    if (tb == 'n' .or. tb == 'N') then
       m2cB(1) = 1 !mode on B to contract 
    else
       m2cB(1) = 2 !mode on B to contract 
    endif
    order(1) = 1 !C(n,m) = A(n,k)*B(k,m)
    order(2) = 2
    nmodes2c = 1   !number of modes to contract 
    IF(A%PDMID.EQ.B%PDMID)THEN
       !Workaround for identical data for A and B
       call mat_pdmm_init(TMP,B%nrow,B%ncol)
       call mat_pdmm_assign(TMP,B) !TMP=B
       call tensor_contract(alpha,pdmm_Marray(A%PDMID)%p,pdmm_Marray(TMP%PDMID)%p,&
            & m2cA,m2cB,nmodes2c,beta,pdmm_Marray(C%PDMID)%p,order)!,&
!            & wrk=wrkTensor,iwrk=iwrkTensor)
       call mat_pdmm_free(TMP)
    ELSE
       call tensor_contract(alpha,pdmm_Marray(A%PDMID)%p,pdmm_Marray(B%PDMID)%p,&
            & m2cA,m2cB,nmodes2c,beta,pdmm_Marray(C%PDMID)%p,order)!,&
!            & wrk=wrkTensor,iwrk=iwrkTensor)
    ENDIF
  end subroutine mat_pdmm_mul

  !> \brief See mat_max_elm in mat-operations.f90
  subroutine mat_pdmm_max_elm(a, val, pos)
    implicit none
    type(matrix),intent(in)  :: a
    real(realk), intent(inout) :: val
    integer    , intent(out)   :: pos(2)
    integer                  :: i,j
    real(realk) ,pointer :: Afull(:,:)
    call mem_alloc(Afull,A%nrow,A%ncol)
    call mat_pdmm_to_full(A,1.0E0_realk,Afull)
    val = afull(1,1); pos(1) = 1; pos(2) = 1 
    do j = 1, a%ncol 
       do i = 1, a%nrow
          if (afull(i,j) > val) then
             val = afull(i,j)
             pos(1) = i
             pos(2) = j
          endif
       enddo
    enddo
    call mem_dealloc(Afull)
  end subroutine mat_pdmm_max_elm

  !> \brief See mat_min_elm in mat-operations.f90
  subroutine mat_pdmm_min_elm(a, val, pos)
    implicit none
    type(matrix),intent(in)  :: a
    real(realk), intent(inout) :: val
    integer    , intent(out)   :: pos(2)
    integer                  :: i,j
    real(realk) ,pointer :: Afull(:,:)
    call mem_alloc(Afull,A%nrow,A%ncol)
    call mat_pdmm_to_full(A,1.0E0_realk,Afull)
    val = afull(1,1); pos(1) = 1; pos(2) = 1 
    do j = 1, a%ncol 
       do i = 1, a%nrow
          if (afull(i,j) < val) then
             val = afull(i,j)
             pos(1) = i
             pos(2) = j
          endif
       enddo
    enddo
    call mem_dealloc(Afull)
  end subroutine mat_pdmm_min_elm

!> \brief See mat_dmul in mat-operations.f90
!> Make c = alpha*diag(a)b + beta*c, 
!> where a is realk(:) b,c are type(matrix) and alpha,beta are parameters
  subroutine mat_pdmm_dmul(a,b,transb,alpha,beta,c) 
    implicit none
    real(realk), intent(in)  :: a(:)
    TYPE(Matrix), intent(IN) :: b
    character, intent(in)    :: transb
    REAL(realk), INTENT(IN)  :: alpha, beta
    TYPE(Matrix)             :: c  
    INTEGER                  :: m,n,i
    real(realk),pointer :: Bfull(:),Cfull(:)   
    call mem_alloc(Bfull,B%nrow*B%ncol)
    call mat_pdmm_to_full(B,1.0E0_realk,Bfull)
    call mem_alloc(Cfull,C%nrow*C%ncol)
    call mat_pdmm_to_full(C,1.0E0_realk,Cfull)

    n = b%nrow
    m = b%ncol
    
    if (ABS(beta).LT.1.0E-15_realk)then 
       call ls_dzero(n*m,cfull)
    elseif(ABS(beta-1.0E0_realk).GT.1.0E-15_realk) then
       call dscal(n*m,beta,cfull,1)
    endif
        
    if (transb == 'n' .or. transb == 'N') then       
       do i=1,n
          call daxpy(m,alpha*a(i),bfull(i),n,cfull(i),n)
       enddo       
    elseif (transb == 't' .or. transb == 'T') then       
       do i=1,m
          call daxpy(n,alpha*a(i),bfull(n*(i-1)+1),1,cfull(i),m)
       enddo       
    else
       print*,'unknown format in mat_pdmm_dmul'
       CALL LSQUIT( 'unknown format in mat_pdmm_dmul',-1)
    endif
    call mem_dealloc(Bfull)
    call mat_pdmm_set_from_full(Cfull,1.0E0_realk,C)    
    call mem_dealloc(Cfull)
    
  end subroutine mat_pdmm_dmul

!> \brief See mat_hmul in mat-operations.f90
!> Make Cij = alpha*Aij*Bij+beta*Cij (Hadamard product)
  subroutine mat_pdmm_hmul(alpha,A,B,beta,C)
    implicit none
    real(realk) :: alpha,beta
    type(matrix),intent(in) :: A,B
    type(matrix),intent(inout) :: C

    call tensor_hmul(alpha,pdmm_Marray(A%PDMID)%p,pdmm_Marray(B%PDMID)%p,&
         & beta,pdmm_Marray(C%PDMID)%p)
    
!    INTEGER                  :: i
!    real(realk),pointer :: Afull(:),Bfull(:),Cfull(:)    
!    call mem_alloc(Afull,A%nrow*A%ncol)
!    call mat_pdmm_to_full(A,1.0E0_realk,Afull)
!    call mem_alloc(Bfull,B%nrow*B%ncol)
!    call mat_pdmm_to_full(B,1.0E0_realk,Bfull)
!    call mem_alloc(Cfull,C%nrow*C%ncol)
!    call mat_pdmm_to_full(C,1.0E0_realk,Cfull)
!    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(Afull,Bfull,Cfull,alpha,beta,A)
!    do i=1,A%nrow*A%ncol
!       Cfull(i)= alpha*Afull(i)*Bfull(i)+beta*Cfull(i)
!    enddo
!    !$OMP END PARALLEL DO
!    call mem_dealloc(Afull)
!    call mem_dealloc(Bfull)
!    call mat_pdmm_set_from_full(Cfull,1.0E0_realk,C)    
!    call mem_dealloc(Cfull)
    
  end subroutine mat_pdmm_hmul
  
  !> \brief See mat_add in mat-operations.f90
  !B = A
  subroutine mat_pdmm_add(alpha,A,beta,B,C) 
    implicit none
    TYPE(Matrix) :: A,B,C
    real(realk) :: Alpha,Beta    
    call mat_pdmm_copy(alpha,a,c)
    call mat_pdmm_daxpy(beta,b,c)
  end subroutine mat_pdmm_add

  !> \brief See mat_add_identity in mat-operations.f90
  !> A = alpha*I + beta*B
  subroutine mat_pdmm_add_identity(alpha, beta, B, A)
    implicit none
    TYPE(Matrix), intent(IN) :: B
    REAL(realk), INTENT(IN)  :: alpha, beta
    TYPE(Matrix)             :: A
    !
    call mat_pdmm_identity(A)
    call mat_pdmm_scal(alpha,A) 
    IF(ABS(beta).GT. 1.0E-17_realk)then
       call mat_pdmm_daxpy(beta,B,A)
    ENDIF
  end subroutine mat_pdmm_add_identity

  subroutine mat_pdmm_identity(I)
    implicit none
    type(Matrix), intent(inout) :: I
    real(realk), pointer    :: ifull(:,:)    
    integer :: j
    call mem_alloc(Ifull,I%nrow,I%ncol)
    call ls_dzero(Ifull,I%nrow*I%ncol)
    do j=1,I%nrow
       ifull(j,j) = 1.0E0_realk
    enddo
    call mat_pdmm_set_from_full(Ifull,1.0E0_realk,I)
    call mem_dealloc(Ifull)

  end subroutine mat_pdmm_identity

  !B = alpha*A + B
  !> \brief See mat_daxpy in mat-operations.f90
  subroutine mat_pdmm_daxpy(alpha,A,B)
    implicit none
    TYPE(Matrix) :: A,B
    REAL(REALK)  :: Alpha
    call tensor_add(pdmm_Marray(B%PDMID)%p,alpha,pdmm_Marray(A%PDMID)%p)
  end subroutine mat_pdmm_daxpy

  !> \brief See mat_dotproduct in mat-operations.f90
  function mat_pdmm_dotproduct(a,b)
     implicit none
     TYPE(Matrix), intent(IN) :: a,b
     REAL(realk) :: mat_pdmm_dotproduct
     mat_pdmm_dotproduct = tensor_ddot(pdmm_Marray(A%PDMID)%p,pdmm_Marray(B%PDMID)%p)
  end function mat_pdmm_dotproduct

  !> \brief See mat_sqnorm2 in mat-operations.f90
  function mat_pdmm_sqnorm2(a)
    implicit none
    TYPE(Matrix), intent(IN) :: a
    REAL(realk) :: mat_pdmm_sqnorm2
    call print_norm(pdmm_Marray(A%PDMID)%p,mat_pdmm_sqnorm2,.TRUE.)
  end function mat_pdmm_sqnorm2

  !> \brief See mat_outdia_sqnorm2 in mat-operations.f90
  function mat_pdmm_outdia_sqnorm2(a)
     implicit none
     TYPE(Matrix), intent(IN) :: a
     REAL(realk) :: mat_pdmm_outdia_sqnorm2
     integer     :: i,n
     mat_pdmm_outdia_sqnorm2 = 0.0E0_realk
     call lsquit('mat_pdmm_outdia_sqnorm2 not implemented',-1)
  end function mat_pdmm_outdia_sqnorm2

  !> \brief See mat_diag_f in mat-operations.f90
  subroutine mat_pdmm_diag_f(F,S,eival,Cmo)
    !solves FC = SCe 
    implicit none
    TYPE(Matrix), intent(IN) :: F,S
    type(matrix)             :: Cmo  !output
    real(realk), intent(inout) :: eival(:)
    real(realk), pointer :: tmp(:,:),CMOfull(:,:)
    integer :: ndim,i
    call mem_alloc(CMOfull,F%nrow,F%ncol)
    call mat_pdmm_to_full(F,1.0E0_realk,CMOfull)
    call mem_alloc(tmp,S%nrow,S%ncol)
    call mat_pdmm_to_full(S,1.0E0_realk,tmp)
    ndim = s%nrow
    call my_DSYGV(ndim,Cmofull,tmp,eival,"PDMM                ")
    call mem_dealloc(tmp)
    call mat_pdmm_set_from_full(CMOfull,1.0E0_realk,CMO)
    call mem_dealloc(CMOfull)
  end subroutine mat_pdmm_diag_f

  !> \brief See mat_dsyev in mat-operations.f90
  subroutine mat_pdmm_dsyev(S,eival,ndim)
    implicit none
    TYPE(Matrix), intent(INOUT) :: S
    real(realk), intent(INOUT) :: eival(ndim)
    real(realk), pointer :: Sfull(:,:)
    integer :: ndim,i
    real(realk),pointer :: work(:)
    integer :: infdiag,lwork
#ifdef VAR_LSESSL
    integer :: liwork
    integer, pointer:: iwork(:)
    liwork=-1
#endif
    infdiag=0
    lwork = -1

    call mem_alloc(Sfull,S%nrow,S%ncol)
    call mat_pdmm_to_full(S,1.0E0_realk,Sfull)
 
    ! we inquire the size of lwork
#ifdef VAR_LSESSL
    call mem_alloc(work,1)
    call mem_alloc(iwork,1)
    iwork = 0.0E0_realk
    !print *,"here1 and trialrun",lwork,liwork,work(1),iwork(1)
    call DSYEVD('V','U',ndim,Sfull,ndim,eival,work,lwork,iwork,liwork,infdiag)
    !print *,work,lwork,iwork,liwork
    liwork = iwork(1)
    call mem_dealloc(iwork)
    call mem_alloc(iwork,liwork)
#else
    call mem_alloc(work,5)
    call DSYEV('V','U',ndim,Sfull,ndim,eival,work,lwork,infdiag)
#endif
    if(infdiag.ne. 0) then
       print*,'mat_dsyev: dsyev query failed, info=',infdiag
       call lsquit('mat_dsyev: query failed.',-1)
    end if

    lwork = NINT(work(1))
    call mem_dealloc(work)
    call mem_alloc(work,lwork)
#ifdef VAR_LSESSL
    !print *,"here1 and run",lwork,liwork,work(1),iwork(1)
    call DSYEVD('V','U',ndim,Sfull,ndim,eival,work,lwork,iwork,liwork,infdiag)
    call mem_dealloc(iwork)
#else
    call DSYEV('V','U',ndim,Sfull,ndim,eival,work,lwork,infdiag)
#endif
    call mem_dealloc(work)
    if(infdiag.ne. 0) then
       print*,'mat_dsyev: dsyev failed, info=',infdiag
       call lsquit('mat_dsyev: diagonalization failed.',-1)
    end if
    call mat_pdmm_set_from_full(Sfull,1.0E0_realk,S)
    call mem_dealloc(Sfull)
  end subroutine mat_pdmm_dsyev

  !> \brief See mat_dsyevx in mat-operations.f90
  subroutine mat_pdmm_dsyevx(S,eival,ieig)
    implicit none
    TYPE(Matrix), intent(INOUT) :: S
    real(realk), intent(INOUT) :: eival
    integer,intent(in) :: ieig
    !
    real(realk),pointer :: work(:),eivec(:,:),elms(:,:)
    integer,pointer     :: icholtemp(:),IFAIL(:)
    integer             :: neig,info,lwork,i,j,ndim
    real(realk),pointer :: eivalTmp(:),B(:,:),C(:,:)
    real(realk)         :: abstol,VL,VU
    real(realk),external :: DLAMCH
    call mem_alloc(elms,S%nrow,S%ncol)
    call mat_pdmm_to_full(S,1.0E0_realk,elms)

    ndim = S%nrow
    abstol =  2.0E0_realk*DLAMCH('S')
    neig = 1
    VL = 0.0E0_realk
    VU = 0.0E0_realk
    INFO = 0
    call mem_alloc(work,8)
    call mem_alloc(icholtemp,5*ndim)
    call mem_alloc(ifail,ndim)
    call mem_alloc(eivalTmp,ndim)
    call mem_alloc(eivec,1,1)
    lwork = -1
    call DSYEVX('N', 'I', 'U', ndim, elms, ndim, VL, VU, ieig, ieig, &
       &  abstol, neig, eivalTmp, eivec, 1, work, lwork, icholtemp, &
       &  IFAIL, INFO )
    lwork = work(1)
    call mem_dealloc(work)
    call mem_alloc(work,lwork)
    call DSYEVX('N', 'I', 'U', ndim, elms, ndim, VL, VU, ieig, ieig, &
       &  abstol, neig, eivalTmp, eivec, 1, work, lwork, icholtemp, &
       &  IFAIL, INFO )
    eival = eivalTmp(1)
    if(info.ne. 0) then
       print*,'mat_pdmm_dsyevx_aux: dsyevx failed, info=',info
       print*,'abstol',abstol
       print*,'lwork',lwork
       call mem_alloc(B,ndim,ndim)
       call mem_alloc(C,ndim,ndim)
       do j = 1,ndim
          do i = 1,ndim
             C(J,I) = elms(I,J)
          enddo
       enddo
       call dcopy(ndim*ndim,elms,1,B,1)
       call dscal(ndim*ndim,0.5E0_realk,B,1)
       call daxpy(ndim*ndim,0.5E0_realk,C,1,B,1)
       call mem_dealloc(C)

       call DSYEVX('N', 'I', 'U', ndim, B, ndim, VL, VU, ieig, ieig, &
            &  abstol, neig, eivalTmp, eivec, 1, work, lwork, icholtemp, &
            &  IFAIL, INFO )
       call mem_dealloc(B)
       if(info.eq. 0) then
          print*,'looks like the matrix was not fully symmetric'
          print*,'dsyevx was succesfull for a symmetriced version'
       else
          print*,'dsyevx also fails for a symmetriced version'
          call lsquit('mat_pdmm_dsyevx_aux: diagonalization failed.',-1)
       endif
    end if
    call mem_dealloc(elms)
    call mem_dealloc(eivec)
    call mem_dealloc(work)
    call mem_dealloc(icholtemp)
    call mem_dealloc(ifail)
    call mem_dealloc(eivalTmp)

  end subroutine mat_pdmm_dsyevx

  !> \brief See mat_assign in mat-operations.f90
  !A = alpha*A
  subroutine mat_pdmm_scal(alpha,A) 
    implicit none
    TYPE(Matrix) :: A
    REAL(REALK),intent(in)  :: alpha
    real(realk),pointer :: AF(:,:)
    call tensor_scale(pdmm_Marray(A%PDMID)%p,alpha)
  end subroutine mat_pdmm_scal

   !> \brief See mat_zero in mat-operations.f90
   subroutine mat_pdmm_zero(A)
     implicit none
     TYPE(Matrix) :: A
     call tensor_zero(pdmm_Marray(A%PDMID)%p)
   end subroutine mat_pdmm_zero

   SUBROUTINE mat_pdmm_dposv(A,b,lupri)
     implicit none
     TYPE(Matrix), intent(INOUT)  :: A
     TYPE(Matrix), intent(INOUT)  :: b
     INTEGER,INTENT(IN)           :: lupri
     Real(Realk),pointer          :: Af(:,:)
     Real(Realk),pointer          :: bf(:,:)
     INTEGER                      :: dim1,dim2,dim3,dim4
     INTEGER                      :: info
     
     dim1=A%nrow
     dim2=A%ncol
     dim3=b%nrow
     dim4=b%ncol
     if(dim3 .ne. A%nrow) then
        call lsquit('mat_pdmm_dposv, Reason: Wrong dim3',lupri)
     endif
     call mem_alloc(Af,dim1,dim2)
     call mem_alloc(bf,dim3,dim4)
     call mat_pdmm_to_full(A,1D0,Af)
     call mat_pdmm_to_full(b,1D0,bf)
     
     call dposv('U',dim1,dim4,Af,dim1,bf,dim3,info)
     
     if(info .ne. 0) then
        call lsquit('mat_pdmm_dposv, Reason: info not 0',lupri)
     endif
     
     call mat_pdmm_set_from_full(Af,1.0E0_realk,A)
     call mat_pdmm_set_from_full(bf,1.0E0_realk,b)
     call mem_dealloc(Af)
     call mem_dealloc(bf)
     
   END SUBROUTINE mat_pdmm_dposv

   !> \brief Compute the Cholesky decomposition factors of a positive definite matrix
   !> \author T. Kjærgaard
   !> \date 2012
   !> \param a The type(matrix) that should be decomposed
   !> \param b The output type(matrix) that contains the cholesky factors.
   SUBROUTINE mat_pdmm_chol(a, b) 
     implicit none
     TYPE(Matrix),intent(in)     :: a
     TYPE(Matrix),intent(inout)  :: b !output
     real(realk), pointer   :: work1(:)
     real(realk), pointer   :: U_full(:,:)
     integer,pointer        :: JPVT(:)
     real(realk)            :: RCOND, dummy(2), tmstart, tmend
     integer                :: IERR, i, j, fulldim, ndim
     fulldim = a%nrow
     call mem_alloc(U_full,fulldim,fulldim) 
     call mem_alloc(work1,fulldim)
     call mem_alloc(JPVT,fulldim)
     call mat_pdmm_to_full(A,1.0E0_realk,U_full)
     !Set lower half of U = 0:
     do i = 1, fulldim
        do j = 1, i-1
           U_full(i,j) = 0.0E0_realk
        enddo
     enddo
     JPVT(1) = 0
     call dchdc(U_full,fulldim,fulldim,work1,JPVT,0,IERR)           
     call mat_pdmm_set_from_full(U_full,1.0E0_realk,B)
     call mem_dealloc(U_full) 
     call mem_dealloc(work1)
     call mem_dealloc(JPVT)
   end SUBROUTINE mat_pdmm_chol

   !> \brief See mat_retrieve_block in mat-operations.f90
  subroutine mat_pdmm_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in)         :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(inout)     :: fullmat(fullrow,fullcol)
    type(Matrix), intent(in) :: A
    integer                     :: i, j
    real(realk),pointer :: Afull(:,:)
    call mem_alloc(Afull,A%nrow,A%ncol)
    call mat_pdmm_to_full(A,1.0E0_realk,Afull)
    do j = insertcol, insertcol+fullcol-1
       do i = insertrow, insertrow+fullrow-1
          fullmat(i-insertrow+1,j-insertcol+1) = Afull(i,j)
       enddo
    enddo
    call mem_dealloc(Afull)
  end subroutine mat_pdmm_retrieve_block

  !> \brief Extract diagonal
  subroutine mat_pdmm_extract_diagonal(diag,A)
    implicit none
    real(realk), intent(inout) :: diag(:)
    type(Matrix), intent(in) :: A
    integer :: n
    real(realk),pointer :: Afull(:)
    call mem_alloc(Afull,A%nrow*A%ncol)
    call mat_pdmm_to_full(A,1.0E0_realk,Afull)
    n=A%nrow
    call dcopy(n,Afull,n+1,diag,1)
    call mem_dealloc(Afull)
  end subroutine mat_pdmm_extract_diagonal

  !> \brief See mat_dger in mat-operations.f90
  SUBROUTINE mat_pdmm_dger(alpha,x,y,A)
    implicit none
    real(realk) :: alpha
    type(matrix) :: A
    real(realk) :: x(:),y(:)
    real(realk),pointer :: Afull(:,:) 
    call mem_alloc(Afull,A%nrow,A%ncol)
    call mat_pdmm_to_full(A,1.0E0_realk,Afull)
    call DGER(A%nrow,A%ncol,alpha,x,1,y,1,Afull,A%nrow)
    call mat_pdmm_set_from_full(Afull,1.0E0_realk,A)
    call mem_dealloc(Afull)
  end SUBROUTINE mat_pdmm_dger

  !> \brief See mat_create_block in mat-operations.f90
  subroutine mat_pdmm_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(in) :: fullmat(fullrow,fullcol)
    type(Matrix), intent(inout) :: A
    real(realk),pointer :: Afull(:,:) 
    integer :: i,j
    call mem_alloc(Afull,A%nrow,A%ncol)
    call mat_pdmm_to_full(A,1.0E0_realk,Afull)
    do j = insertcol, insertcol+fullcol-1
       do i = insertrow, insertrow+fullrow-1
          Afull(i,j) = fullmat(i-insertrow+1,j-insertcol+1)
       enddo
    enddo
    call mat_pdmm_set_from_full(Afull,1.0E0_realk,A)
    call mem_dealloc(Afull)
  end subroutine mat_pdmm_create_block
#else
  contains
    subroutine dummy_no_use
    end subroutine dummy_no_use
#endif

end module matrix_operations_pdmm

#ifdef VAR_ENABLE_TENSORS
SUBROUTINE PDMM_GRIDINIT(NBAST)
  use matrix_operations_pdmm
  use matrix_module
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif
  IMPLICIT NONE
  INTEGER(KIND=LS_MPIK) :: PR, PC, TMP(2), IERR, I, NBAST,nSlaves
#ifdef VAR_MPI
  pdmm_collector=0 !collect on master, rank = 0
  call ls_mpibcast(PDMMGRIDINIT,infpar%master,MPI_COMM_LSDALTON)   
  IF(infpar%inputblocksize.EQ.0)THEN
     IF(NBAST.GT.4096)THEN
        infpar%inputblocksize = 2048
     ELSE !debug 
        infpar%inputblocksize = NBAST/(infpar%nodtot)
     ENDIF
  ENDIF
  ENTRY PDMM_GRIDINIT_SLAVE
  call ls_mpibcast(infpar%inputblocksize,infpar%master,pdmm_comm)
  BLOCK_SIZE_PDM = infpar%inputblocksize
#endif      
  call InitMatrixModulePDM()
END SUBROUTINE PDMM_GRIDINIT

SUBROUTINE PDMM_GRIDEXIT()
   use matrix_operations_pdmm
#ifdef VAR_MPI
   use infpar_module
   use lsmpi_type
#endif
   IMPLICIT NONE
#ifdef VAR_MPI
   call ls_mpibcast(PDMMGRIDEXIT,infpar%master,MPI_COMM_LSDALTON)
#endif
   ENTRY PDMM_GRIDEXIT_SLAVE
   call FreeMatrixModulePDM
 END SUBROUTINE PDMM_GRIDEXIT
#endif

