!> @file
!> Contains dense matrix module.

!> \brief Contains matrix operation routines for type(matrix) = dense.
!>
!> Modified to use blas and lapack routines by B. Jansik, Arhus, Jun 2006. \n
!> ajt mar09 Generalize to complex matrices. \n
!> ajt FIXME Currently both re and im in elms(:). Should use celms/ielms
!>
module matrix_operations_dense
  use memory_handling
  use matrix_module
#ifdef VAR_ENABLE_TENSORS
  use reorder_frontend_module
#endif
  use precision
  contains
!> \brief See mat_init in mat-operations.f90
  subroutine mat_dense_init(A,nrow,ncol)
     implicit none
     TYPE(Matrix) :: A
     integer, intent(in) :: nrow, ncol
     integer(kind=long) :: nsize
     !A%elms should not be associated with anything.
     !if so, the memory will be lost!!
     A%nrow = nrow
     A%ncol = ncol
     allocate(A%elms(A%nrow * A%ncol &
                   & * merge(2,1,A%complex)))
     nsize = size(A%elms,KIND=long)*mem_realsize
     call mem_allocated_mem_type_matrix(nsize)
  end subroutine mat_dense_init

!> \brief See mat_free in mat-operations.f90
  subroutine mat_dense_free(a)
     implicit none
     TYPE(Matrix) :: a
     integer(kind=long) :: nsize
     if (.not.ASSOCIATED(a%elms)) then
       print*,'memory previously released!!'
       CALL LSQUIT( 'Error in mat_dense_free - memory previously released',-1)
     endif
     nsize = SIZE(a%elms,KIND=long)*mem_realsize
     call mem_deallocated_mem_type_matrix(nsize)
     DEALLOCATE(a%elms)
     NULLIFY(a%elms)
  end subroutine mat_dense_free

!> \brief See mat_set_from_full in mat-operations.f90
  subroutine mat_dense_set_from_full(afull,alpha,a)
     implicit none
     real(realk), INTENT(IN) :: afull(*)
     real(realk), intent(in) :: alpha
     TYPE(Matrix)            :: a 
     integer                 :: i

!    do i = 1,a%nrow*a%ncol
!      a%elms(i) = alpha*afull(i)
!    enddo

     i = a%nrow*a%ncol
     call dcopy (i,afull,1,a%elms,1)
     if (ABS(alpha-1.0E0_realk).GT.1.0E-15_realk) call mat_dense_scal(alpha,a)

  end subroutine mat_dense_set_from_full

!> \brief See mat_to_full in mat-operations.f90
  subroutine mat_dense_to_full(a, alpha, afull)
     implicit none
     TYPE(Matrix), intent(in):: a
     real(realk), intent(in) :: alpha
     real(realk), intent(inout):: afull(*)  
     integer                 :: i

!    do i = 1,a%nrow*a%ncol
!      afull(i) = alpha*a%elms(i)
!    enddo

     i = a%nrow*a%ncol
     call dcopy (i,a%elms,1,afull,1)
     if (ABS(alpha-1.0E0_realk).GT.1.0E-15_realk) call dscal(i,alpha,afull,1)

  end subroutine mat_dense_to_full

!> \brief See mat_to_full3D in mat-operations.f90
  subroutine mat_dense_to_full3D(a, alpha, afull,n1,n2,n3,i3)
     implicit none
     integer, INTENT(IN)           :: n1,n2,n3,i3
     TYPE(Matrix), intent(in):: a
     real(realk), intent(in) :: alpha
     real(realk), intent(inout):: afull(n1,n2,n3)  
     integer                 :: i,N,M,MP1,j,offset
     N = a%nrow
     M = MOD(N,7)
     IF (M.NE.0) THEN
        do j = 1,a%ncol
           offset = (j-1)*N
           DO I = 1,M
              afull(i,j,i3) = alpha*a%elms(i+offset)              
           ENDDO
        enddo
        MP1 = M + 1
        IF (N.GE.7)THEN
           do j = 1,a%ncol
              offset = (j-1)*N
              DO I = MP1,N,7
                 afull(i,j,i3) = alpha*a%elms(i+offset)
                 afull(i+1,j,i3) = alpha*a%elms(i+1+offset)
                 afull(i+2,j,i3) = alpha*a%elms(i+2+offset)
                 afull(i+3,j,i3) = alpha*a%elms(i+3+offset)
                 afull(i+4,j,i3) = alpha*a%elms(i+4+offset)
                 afull(i+5,j,i3) = alpha*a%elms(i+5+offset)
                 afull(i+6,j,i3) = alpha*a%elms(i+6+offset)
              END DO
           enddo
        ENDIF
     ELSE
        do j = 1,a%ncol
           offset = (j-1)*N
           DO I = 1,N,7
              afull(i,j,i3) = alpha*a%elms(i+offset)
              afull(i+1,j,i3) = alpha*a%elms(i+1+offset)
              afull(i+2,j,i3) = alpha*a%elms(i+2+offset)
              afull(i+3,j,i3) = alpha*a%elms(i+3+offset)
              afull(i+4,j,i3) = alpha*a%elms(i+4+offset)
              afull(i+5,j,i3) = alpha*a%elms(i+5+offset)
              afull(i+6,j,i3) = alpha*a%elms(i+6+offset)
           END DO
        ENDDO
     ENDIF
  end subroutine mat_dense_to_full3D

!> \brief See mat_print in mat-operations.f90
  subroutine mat_dense_print(a, i_row1, i_rown, j_col1, j_coln, lu)
     implicit none
     TYPE(Matrix),intent(in) :: a
     integer, intent(in)     :: i_row1, i_rown, j_col1, j_coln, lu
     !to be found in pdpack/printpkg.F
     call LS_OUTPUT(a%elms, i_row1, i_rown, j_col1, j_coln, a%nrow, a%ncol, 1, lu)
  end subroutine mat_dense_print

!> \brief See mat_trans in mat-operations.f90
  subroutine mat_dense_trans(a,b)
     implicit none
     type(matrix),intent(in)   :: a
     type(matrix)              :: b   !output  
     integer                   :: i, j
!!#ifdef VAR_MKL
!!     call mkl_domatcopy('R', 'T', A%nrow, A%ncol, 1.0E0_realk, A%elms, 1, B%elms,1)
!!#else
!     call ls_transpose(A%elms,B%elms,A%nrow)
#ifdef VAR_ENABLE_TENSORS
     call mat_transpose(a%nrow,a%ncol,1.0E0_realk,a%elms,0.0E0_realk,b%elms)
#else
     do j = 1,a%ncol
       do i = 1,a%nrow
         b%elms(b%nrow*(i-1)+j) = a%elms(a%nrow*(j-1)+i)
       enddo
     enddo
#endif
  end subroutine mat_dense_trans

!> \brief See mat_assign in mat-operations.f90
  subroutine mat_dense_assign(a,b)
     implicit none
     TYPE(Matrix), INTENT(INOUT) :: a
     TYPE(Matrix), INTENT(IN)    :: b
     integer                     :: i

!    do i = 1,a%nrow*a%ncol
!      a%elms(i) = b%elms(i)
!    enddo

     i = a%nrow*a%ncol
     call dcopy (i,b%elms,1,a%elms,1)

  end subroutine mat_dense_assign


#ifndef UNITTEST
!> \brief See mat_mpicopy in mat-operations.f90
  subroutine mat_dense_mpicopy(a,slave,master)
#ifdef VAR_MPI
    use lsmpi_type
#endif
     implicit none
     TYPE(Matrix), TARGET, INTENT(INOUT) :: a
     logical,intent(in)               :: slave
     integer(kind=ls_mpik),intent(in) :: master
     integer                          :: i
#ifdef VAR_MPI
     call LS_MPI_BUFFER(a%nrow,Master)
     call LS_MPI_BUFFER(a%ncol,Master)
     call LS_MPI_BUFFER(a%complex,Master)
     IF(SLAVE)THEN
        NULLIFY(a%iaux, a%raux)
        nullify(A%elms)
        nullify(A%elmsb)                
        call mat_dense_init(a,a%nrow,a%ncol)
        a%init_self_ptr => a
        a%init_magic_tag = mat_init_magic_value
     ENDIF
     i = size(a%elms) !a%nrow*a%ncol
     call LS_MPI_BUFFER(a%elms,i,Master)
#else
     call lsquit('call mat_dense_mpicopy, without LSMPI',-1)
#endif
   end subroutine mat_dense_mpicopy
#endif

!> \brief See mat_copy in mat-operations.f90
  subroutine mat_dense_copy(alpha,a,b)
     implicit none
     REAL(REALK),  INTENT(IN)    :: alpha
     TYPE(Matrix), INTENT(IN)    :: a
     TYPE(Matrix), INTENT(INOUT) :: b

!    do i = 1,a%nrow*a%ncol
!      b%elms(i) = alpha*a%elms(i)
!    enddo

     call mat_dense_assign(b,a)
     if (ABS(alpha-1.0E0_realk).GT.1.0E-15_realk) call mat_dense_scal(alpha,b)

  end subroutine mat_dense_copy

!> \brief See mat_tr in mat-operations.f90
  function mat_dense_Tr(a)
     implicit none
     TYPE(Matrix), intent(IN) :: a
     REAL(realk) :: mat_dense_tr 
     integer :: i
  
     mat_dense_tr = 0E0_realk
     do i = 1,a%nrow
       mat_dense_tr = mat_dense_tr + a%elms((a%Nrow+1)*i-a%Nrow)
     enddo

  end function mat_dense_Tr
 
!> \brief See mat_TrAB in mat-operations.f90
  function mat_dense_TrAB(a,b)
     implicit none
     TYPE(Matrix), intent(IN) :: a,b
     REAL(realk) :: mat_dense_trAB 
     real(realk), external :: ddot
     integer :: j

     mat_dense_TrAB = 0.0E0_realk
     do j = 1,a%ncol
!      do i = 1,a%nrow
!        mat_dense_TrAB = mat_dense_TrAB + a%elms(a%nrow*(j-1)+i)*&
!                                         &b%elms(b%nrow*(i-1)+j)
!      enddo
       mat_dense_TrAB = mat_dense_TrAB + &
           & ddot(a%nrow,a%elms(a%nrow*(j-1)+1),1,b%elms(j),a%ncol)
     enddo

  end function mat_dense_TrAB 

!> \brief See mat_mul in mat-operations.f90
  subroutine mat_dense_mul(a,b,transa,transb,alpha,beta,c) 
     implicit none
     TYPE(Matrix), intent(IN) :: a, b
     character, intent(in)    :: transa, transb
     REAL(realk), INTENT(IN)  :: alpha, beta
     TYPE(Matrix)             :: c  
     INTEGER                  :: m,n,k

     if (transa == 'n' .or. transa == 'N') then
       m = a%nrow
       k = a%ncol
     elseif (transa == 't' .or. transa == 'T') then
       m = a%ncol
       k = a%nrow
     else
       print*,'unknown format in mat_dense_mul'
       CALL LSQUIT( 'unknown format in mat_dense_mul',-1)
     endif
     if (transb == 'n' .or. transb == 'N') then
       n = b%ncol
       if (b%nrow /= k) CALL LSQUIT( 'weird',-1)
     elseif (transb == 't' .or. transb == 'T') then
       n = b%nrow
       if (b%ncol /= k) CALL LSQUIT( 'weird',-1)
     endif
!write (mat_lu,*) 'matrix a'
!call mat_dense_print(a, 1, a%nrow, 1, a%ncol, mat_lu)
!write (mat_lu,*) 'matrix b'
!call mat_dense_print(b, 1, b%nrow, 1, b%ncol, mat_lu)
!write (mat_lu,*) 'matrix c'
!call mat_dense_print(c, 1, c%nrow, 1, c%ncol, mat_lu)
     call DGEMM(transa,transb,m,n,k,alpha,&
               &a%elms,a%nrow,b%elms,b%nrow,beta,c%elms,c%nrow)

  end subroutine mat_dense_mul

!> \brief See mat_dmul in mat-operations.f90
  subroutine mat_dense_dmul(a,b,transb,alpha,beta,c) 
     implicit none
     real(realk), intent(in)  :: a(:)
     TYPE(Matrix), intent(IN) :: b
     character, intent(in)    :: transb
     REAL(realk), INTENT(IN)  :: alpha, beta
     TYPE(Matrix)             :: c  
     INTEGER                  :: m,n,i

     n = b%nrow
     m = b%ncol
       
     if (ABS(beta).LT.1.0E-15_realk)then 
        c%elms=0E0_realk
     elseif (ABS(beta-1.0E0_realk).GT.1.0E-15_realk) then
        call dscal(n*m,beta,c%elms,1)
     endif


     if (transb == 'n' .or. transb == 'N') then

       do i=1,n
        call daxpy(m,alpha*a(i),b%elms(i),n,c%elms(i),n)
       enddo

     elseif (transb == 't' .or. transb == 'T') then

       do i=1,m
        call daxpy(n,alpha*a(i),b%elms(n*(i-1)+1),1,c%elms(i),m)
       enddo

     else
       print*,'unknown format in mat_dense_dmul'
       CALL LSQUIT( 'unknown format in mat_dense_dmul',-1)
     endif
  end subroutine mat_dense_dmul

!> \brief Extract diagonal
  subroutine mat_dense_extract_diagonal(diag,A)
      implicit none
      real(realk), intent(inout) :: diag(:)
      type(Matrix), intent(in) :: A
      integer :: n

      n=A%nrow
      call dcopy(n,A%elms,n+1,diag,1)

  end subroutine mat_dense_extract_diagonal

!> \brief See mat_add in mat-operations.f90
  subroutine mat_dense_add(alpha,a,beta,b,c)
     implicit none
     TYPE(Matrix), intent(IN) :: a, b
     REAL(realk), INTENT(IN)  :: alpha, beta
     TYPE(Matrix)             :: c

!    do i = 1,a%nrow*a%ncol
!      c%elms(i) = alpha*a%elms(i) + beta*b%elms(i)
!    enddo

     call mat_dense_copy(alpha,a,c)
     call mat_dense_daxpy(beta,b,c)
     
  end subroutine mat_dense_add

!> \brief See mat_daxpy in mat-operations.f90
  subroutine mat_dense_daxpy(alpha,x,y)
     implicit none
     real(realk),intent(in)       :: alpha
     TYPE(Matrix), intent(IN)     :: X
     TYPE(Matrix), intent(INOUT)  :: Y
     integer                      :: i

!    do i = 1,x%nrow*x%ncol
!      y%elms(i) = y%elms(i) + alpha*x%elms(i)
!    enddo

     i = x%nrow*x%ncol
     call daxpy(i,alpha,x%elms,1,y%elms,1)

  end subroutine mat_dense_daxpy

!> \brief See mat_dotproduct in mat-operations.f90
  function mat_dense_dotproduct(a,b)
     implicit none
     TYPE(Matrix), intent(IN) :: a,b
     REAL(realk) :: mat_dense_dotproduct
     real(realk), external :: ddot
     integer     :: i

!     mat_dense_dotproduct = 0.0E0_realk
!    do i = 1,a%nrow*a%ncol
!      mat_dense_dotproduct = mat_dense_dotproduct + a%elms(i)*b%elms(i)
!    enddo

     i = a%nrow*a%ncol
     mat_dense_dotproduct = ddot(i,a%elms,1,b%elms,1)

  end function mat_dense_dotproduct

!> \brief See mat_sqnorm2 in mat-operations.f90
  function mat_dense_sqnorm2(a)
     implicit none
     TYPE(Matrix), intent(IN) :: a
     REAL(realk) :: mat_dense_sqnorm2
!    integer     :: i

!     mat_dense_sqnorm2 = 0.0E0_realk
!    do i = 1,a%nrow*a%ncol
!      mat_dense_sqnorm2 = mat_dense_sqnorm2 + a%elms(i)*a%elms(i)
!    enddo
     mat_dense_sqnorm2 = mat_dense_dotproduct(a,a)
      
  end function mat_dense_sqnorm2

!> \brief See mat_outdia_sqnorm2 in mat-operations.f90
  function mat_dense_outdia_sqnorm2(a)
     implicit none
     TYPE(Matrix), intent(IN) :: a
     REAL(realk) :: mat_dense_outdia_sqnorm2
     integer     :: i,n

     mat_dense_outdia_sqnorm2 = 0.0E0_realk
     n = 1
     do i = 1,a%nrow*a%ncol
       if (i /= (n-1)*a%nrow+n) then
         mat_dense_outdia_sqnorm2 = mat_dense_outdia_sqnorm2 + a%elms(i)*a%elms(i)
       else
         n = n+1
       endif
     enddo
  end function mat_dense_outdia_sqnorm2

!> \brief See mat_dsyev in mat-operations.f90
  SUBROUTINE mat_dense_dsyev(S,eival,ndim)
    implicit none
    TYPE(Matrix), intent(INOUT) :: S
    integer,intent(in) :: ndim
    real(realk), intent(INOUT) :: eival(ndim)
!
    real(realk),pointer :: work(:)
    integer :: infdiag,lwork
#ifdef VAR_LSESSL
    integer :: liwork
    integer, pointer:: iwork(:)
    liwork=-1
#endif
    infdiag=0
    lwork = -1
 
    ! we inquire the size of lwork
#ifdef VAR_LSESSL
    call mem_alloc(work,1)
    call mem_alloc(iwork,1)
    iwork = 0.0E0_realk
    !print *,"here1 and trialrun",lwork,liwork,work(1),iwork(1)
    call DSYEVD('V','U',ndim,S%elms,ndim,eival,work,lwork,iwork,liwork,infdiag)
    !print *,work,lwork,iwork,liwork
    liwork = iwork(1)
    call mem_dealloc(iwork)
    call mem_alloc(iwork,liwork)
#else
    call mem_alloc(work,5)
    call DSYEV('V','U',ndim,S%elms,ndim,eival,work,lwork,infdiag)
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
    call DSYEVD('V','U',ndim,S%elms,ndim,eival,work,lwork,iwork,liwork,infdiag)
    call mem_dealloc(iwork)
#else
    call DSYEV('V','U',ndim,S%elms,ndim,eival,work,lwork,infdiag)
#endif
    call mem_dealloc(work)
    if(infdiag.ne. 0) then
       print*,'mat_dsyev: dsyev failed, info=',infdiag
       call lsquit('mat_dsyev: diagonalization failed.',-1)
    end if
  END SUBROUTINE mat_dense_dsyev

!> \brief See mat_dsyevx in mat-operations.f90
  SUBROUTINE mat_dense_dsyevx(S,eival,ieig)
    implicit none
    TYPE(Matrix), intent(INOUT) :: S
    real(realk), intent(INOUT) :: eival
    integer,intent(in) :: ieig
!
    integer :: ndim

    ndim = S%nrow
    call mat_dense_dsyevx_aux(S%elms,eival,ndim,ieig)

  END SUBROUTINE mat_dense_dsyevx

!> \brief See mat_dsyevx in mat-operations.f90
  SUBROUTINE mat_dense_dsyevx_aux(elms,eival,ndim,ieig)
    implicit none
    integer,intent(in)         :: ndim,ieig
    Real(realk), intent(INOUT) :: elms(ndim,ndim)
    real(realk), intent(INOUT) :: eival
!
    real(realk),pointer :: work(:),eivec(:,:)
    integer,pointer     :: icholtemp(:),IFAIL(:)
    integer             :: neig,info,lwork,i,j
    real(realk),pointer :: eivalTmp(:),B(:,:),C(:,:)
    real(realk)         :: abstol,VL,VU
    real(realk),external :: DLAMCH

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
       print*,'mat_dense_dsyevx_aux: dsyevx failed, info=',info
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
          call lsquit('mat_dense_dsyevx_aux: diagonalization failed.',-1)
       endif
    end if
    call mem_dealloc(eivec)
    call mem_dealloc(work)
    call mem_dealloc(icholtemp)
    call mem_dealloc(ifail)
    call mem_dealloc(eivalTmp)
  END SUBROUTINE mat_dense_dsyevx_aux


!> \brief See mat_diag_f in mat-operations.f90
  subroutine mat_dense_diag_f(F,S,eival,Cmo)
    !solves FC = SCe 
    implicit none
    TYPE(Matrix), intent(IN) :: F,S
    type(matrix)             :: Cmo  !output
    real(realk), intent(inout) :: eival(:)
    real(realk), allocatable :: tmp(:)
    integer :: ndim,i
    ndim = s%nrow
    ALLOCATE(tmp(Ndim*Ndim))
    do i = 1,ndim*ndim
       tmp(i) = S%elms(i)
    enddo
    call mat_dense_assign(Cmo, F)
    call my_DSYGV(ndim,Cmo%elms,tmp,eival,"DENSE_GET_DENS_SOL  ")
    DEALLOCATE(tmp)
  end subroutine mat_dense_diag_f

!> \brief See mat_abs_max_elm in mat-operations.f90
  subroutine mat_dense_abs_max_elm(a, val)
  implicit none
     type(matrix),intent(in)  :: a
     real(realk), intent(inout) :: val
     integer                  :: i

     val = abs(a%elms(1))
     do i = 1, a%nrow*a%ncol 
        if (abs(a%elms(i)) > val) then
           val = abs(a%elms(i))
        endif
     enddo
  end subroutine mat_dense_abs_max_elm

!> \brief See mat_max_elm in mat-operations.f90
  subroutine mat_dense_max_elm(a, val, pos)
  implicit none
     type(matrix),intent(in)  :: a
     real(realk), intent(inout) :: val
     integer    , intent(out)   :: pos(2)
     integer                  :: i, imax

     val = a%elms(1); imax=1;
     do i = 2, a%nrow*a%ncol 
        if (a%elms(i) > val) then
           val = a%elms(i)
           imax=i
        endif
     enddo

     !in base 0
     imax   = imax -1
     pos(1) = mod(imax,a%nrow)
     pos(2) = (imax - pos(1))/a%nrow

     !in base 1
     pos = pos + 1
     
  end subroutine mat_dense_max_elm

  subroutine mat_dense_min_elm(a, val, pos)
  implicit none
     type(matrix),intent(in)  :: a
     real(realk), intent(inout) :: val
     integer    , intent(out)   :: pos(2)
     integer                  :: i, imax
    
    val = a%elms(1); imax=1;
     do i = 2, a%nrow*a%ncol 
        if (a%elms(i) < val) then
           val = a%elms(i)
           imax=i
        endif
     enddo

     !in base 0
     imax   = imax -1
     pos(1) = mod(imax,a%nrow)
     pos(2) = (imax - pos(1))/a%nrow

     !in base 1
     pos = pos + 1
     
  end subroutine mat_dense_min_elm


!> \brief See mat_max_diag_elm in mat-operations.f90
  subroutine mat_dense_max_diag_elm(a,pos,val)
  implicit none
     type(matrix),intent(in)  :: a
     real(realk), intent(inout) :: val
     integer                  :: i
     integer, intent(inout)     :: pos

     val = abs(a%elms(1))
     pos = 1
     do i = 2, a%nrow
        if (abs(a%elms(a%nrow*(i-1)+i)) > abs(val)) then
           val = a%elms(a%nrow*(i-1)+i)
           pos = i
        endif
     enddo

  end subroutine mat_dense_max_diag_elm

!> \brief See mat_dE_dmu in mat-operations.f90
  FUNCTION mat_dense_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
     !Computes dE_SCF/dmu where mu is the damping in damped roothan
     !Find dE/dmu = sum_nu,I dE/dC_nu,I dC_nu,I/dmu
     implicit none
     type(matrix),intent(in) :: Fnew,Fdamp,SDS,Cmo
     integer, intent(in)     :: nocc
     real(realk)             :: mat_dense_dE_dmu
     type(matrix)            :: CTFC, tSDS, Cmo_occ, X
     real(realk)             :: tFdamp_dia(Fnew%nrow)
     integer                 :: ndim,i,k,m

     ndim = Fnew%nrow
   ![dE/dC]ij = 4[FC]ij  
     call mat_dense_init(CTFC,ndim,nocc)
     call mat_dense_init(X,ndim,nocc)
     call mat_dense_init(cmo_occ,ndim,nocc)
     call mat_dense_section(Cmo,1,ndim,1,nocc,Cmo_occ)
     call mat_dense_mul(Fnew,Cmo_occ,'n','n',1E0_realk,0E0_realk,X)
     call mat_dense_mul(Cmo,X,'t','n',4.0E0_realk,0E0_realk,CTFC)
   ![dC/dmu]ij = -[dC/dmu]ji => [dC/dmu]ij = 0 for i=j
   ![dC/dmu]ij = [C^TSDSC]ij/(C^TFdampC_ii - C^TFdampC_jj))
     call mat_dense_init(tSDS,ndim,nocc)
     call mat_dense_mul(SDS,Cmo_occ,'n','n',1E0_realk,0E0_realk,X)
     call mat_dense_mul(cmo,X,'t','n',1E0_realk,0E0_realk,tSDS)
     call mat_dense_free(X)
     call mat_dense_free(Cmo_occ)
     !construct the C^TFdampC diagonal - maybe some day a mat_get_dia routine 
     !will be made and some of this would be replaced with a call for that
     !tFdamp_ii = sum_km (Fdamp_mk*C_mi*C_ki)
     do i = 1,ndim
       tFdamp_dia(i) = 0.0E0_realk
       do k = 1,ndim  !col
         do m = 1,ndim  !row
           tFdamp_dia(i) = tFdamp_dia(i) + Fdamp%elms((k-1)*Fdamp%nrow+m)*&
                       &Cmo%elms((i-1)*Cmo%nrow+m)*Cmo%elms((i-1)*Cmo%nrow+k)
         enddo
       enddo
     enddo
!FIXME: make only sum over virt,occ
   !dE/dmu = sum_nu,I dE/dC_nu,I dC_nu,I/dmu
     mat_dense_dE_dmu = 0.0E0_realk
     do i = 1,nocc  !I
       do m = 1,ndim  !J
         if (i /= m) then                            !CTFC(m,i)
           mat_dense_dE_dmu = mat_dense_dE_dmu + CTFC%elms((i-1)*CTFC%nrow+m)*&
                 &tSDS%elms((i-1)*tSDS%nrow+m)/(tFdamp_dia(m)-tFdamp_dia(i))
                      !tSDS(m,i)
         endif
       enddo
     enddo     
     call mat_dense_free(CTFC)
     call mat_dense_free(tSDS)
  END FUNCTION mat_dense_dE_dmu

!> \brief See mat_column_norm in mat-operations.f90
  FUNCTION mat_dense_column_norm(Mat,ncol,from_row,to_row)
  !Returns the sum of the elements A(from_row:to_row,col_num) squared
     implicit none
     type(Matrix), intent(in) :: Mat
     integer, intent(in) :: ncol, from_row, to_row
     real(realk) :: mat_dense_column_norm
     integer :: i
  
     mat_dense_column_norm = 0E0_realk
     do i = from_row,to_row
       mat_dense_column_norm = mat_dense_column_norm + &
                       & (Mat%elms((ncol-1)*mat%nrow+i))**2
     enddo  
  end function mat_dense_column_norm

!> \brief See mat_section in mat-operations.f90
  subroutine mat_dense_section(A,from_row,to_row,from_col,to_col,Asec)
     implicit none
     type(Matrix), intent(in) :: A
     integer, intent(in) :: from_row, to_row, from_col, to_col
     type(Matrix), intent(inout) :: Asec  !output
     integer :: i, index, irow, icol

     Asec%nrow = to_row - from_row + 1
     Asec%ncol = to_col - from_col + 1
     if (SIZE(Asec%elms) < Asec%nrow*Asec%ncol) CALL LSQUIT( 'too little memory allocated in mat_dense_section',-1)
     index = (from_col-1)*A%nrow + from_row
     irow = 0
     icol = from_col
     do i = 1,Asec%nrow*Asec%ncol
       Asec%elms(i) = A%elms(index)
       irow = irow + 1
       if (irow == Asec%nrow) then
         icol = icol + 1
         index = (icol - 1)*A%nrow + from_row
         irow = 0 !Added 17/11-05. Stinne thinks she has found a bug :-)
       else
         index = index + 1
       endif
     enddo
  end subroutine mat_dense_section

!> \brief See mat_insert_section in mat-operations-essentials.f90
subroutine mat_dense_insert_section (Asec, from_row, to_row, from_col, to_col, A)

implicit none

type(matrix), intent(in)    :: Asec
integer, intent(in)         :: from_row, to_row, from_col, to_col
type(matrix), intent(inout) :: A
integer                     :: i, index, irow, icol

index = (from_col-1)*A%nrow + from_row
irow = 0
icol = from_col
do i = 1,Asec%nrow*Asec%ncol
  A%elms(index) = Asec%elms(i)
  irow = irow + 1
  if (irow == Asec%nrow) then
    icol = icol + 1
    index = (icol - 1)*A%nrow + from_row
    irow = 0
  else
    index = index + 1
  endif
enddo

end subroutine mat_dense_insert_section
  
!> \brief See mat_section2 in mat-operations.f90
  subroutine mat_dense_section2(A,from_row,to_row,from_col,to_col,B)
     implicit none
     type(Matrix), intent(in) :: A
     integer, intent(in) :: from_row, to_row, from_col, to_col
     type(Matrix), intent(inout) :: B  !output
     integer :: i, j, k, ndim

     ndim = B%nrow

     k = 1
     do i = from_col, to_col
        do j = from_row, to_row
           B%elms((i-1)*ndim+j) = A%elms(k)
           k = k + 1
        enddo
     enddo 

  end subroutine mat_dense_section2

!> \brief See mat_precond in mat-operations.f90
  subroutine mat_dense_precond(M,x,xprec)
    implicit none
    type(Matrix), intent(inout) :: xprec
    type(Matrix), intent(in)    :: M, x
    integer :: i

    do i = 1, x%nrow
       if (abs(M%elms((i-1)*M%nrow+i)) > 1.0E-9_realk) then
          xprec%elms(i) = x%elms(i)/M%elms((i-1)*M%nrow+i)
       endif
    enddo

  end subroutine mat_dense_precond

!> \brief See mat_mo_precond in mat-operations.f90
  subroutine mat_dense_mo_precond(nocc,omega,Eorb_final,X_MO)
    implicit none
    integer, intent(in) :: nocc
    real(realk), intent(in) :: omega,  Eorb_final(:)
    type(Matrix), intent(inout) :: X_MO
    real(realk) :: dia
    integer :: i,j

    !modify lower left block (ai block)
    do j = 1,nocc  !columns
      do i = nocc+1,X_MO%nrow  !rows
        !dia = E[2]_dia + omega*S[2]_dia
        !NOTE TO SELF: WHY the 2.d0s and not 1E0_realk????
        dia = 2.0E0_realk*(Eorb_final(i) - Eorb_final(j)) + omega * 2E0_realk
        X_MO%elms((j-1)*X_MO%nrow+i) = X_MO%elms((j-1)*X_MO%nrow+i) /dia
      enddo
    enddo
    !modify upper right block (ia block) 
    do j = nocc+1,X_MO%ncol  !columns
      do i = 1,nocc          !rows
        !dia = E[2]_dia + omega*S[2]_dia
        dia = 2.0E0_realk*(Eorb_final(j) - Eorb_final(i)) - omega * 2E0_realk
        X_MO%elms((j-1)*X_MO%nrow+i) = X_MO%elms((j-1)*X_MO%nrow+i)/dia
      enddo
    enddo

  end subroutine mat_dense_mo_precond

  subroutine mat_dense_new_mo_precond(nocc,omega,Eorb_final,Xp_MO,Xm_mo)
    implicit none
    integer, intent(in) :: nocc
    real(realk), intent(in) :: omega,  Eorb_final(:)
    type(Matrix), intent(inout) :: Xp_MO,xm_mo
    real(realk) :: dia,aa
    integer :: i,j
    
    !modify lower left block (ai block)
    do j = 1,nocc  !columns
      do i = nocc+1,Xp_MO%nrow  !rows
        !dia = E[2]_dia + omega*S[2]_dia
        !NOTE TO SELF: WHY the 2.d0s and not 1E0_realk????
        aa=2E0_realk*(Eorb_final(i) - Eorb_final(j))*(Eorb_final(i)- Eorb_final(j))-2E0_realk*omega*omega
        !dia = 2.0E0_realk*(Eorb_final(i) - Eorb_final(j)) + omega * 2E0_realk
        dia=Xp_mo%elms((j-1)*Xp_mo%nrow+i)
        xp_mo%elms((j-1)*Xp_MO%nrow+i) = ((Eorb_final(i) - Eorb_final(j))*(Xp_mo%elms((j-1)*Xp_mo%nrow+i))&
        &-omega*Xm_mo%elms((j-1)*Xm_mo%nrow+i)) /aa
        xm_mo%elms((j-1)*Xm_MO%nrow+i) = ((Eorb_final(i) - Eorb_final(j))*(Xm_mo%elms((j-1)*Xm_mo%nrow+i))&
        &-omega*dia) /aa
      enddo
    enddo
    !modify upper right block (ia block) 
    do j = nocc+1,Xp_MO%ncol  !columns
      do i = 1,nocc          !rows
        !dia = E[2]_dia + omega*S[2]_dia
        aa=2E0_realk*(Eorb_final(j) - Eorb_final(i))*(Eorb_final(j)- Eorb_final(i))-2E0_realk*omega*omega
        !dia = 2.0E0_realk*(Eorb_final(j) - Eorb_final(i)) - omega * 2E0_realk
        dia=Xp_mo%elms((j-1)*Xp_mo%nrow+i)
        xp_mo%elms((j-1)*Xp_MO%nrow+i) = ((Eorb_final(j) - Eorb_final(i))*(Xp_mo%elms((j-1)*Xp_mo%nrow+i))&
        &+omega*Xm_mo%elms((j-1)*Xm_mo%nrow+i)) /aa
        xm_mo%elms((j-1)*Xm_MO%nrow+i) = ((Eorb_final(j)- Eorb_final(i))*(Xm_mo%elms((j-1)*Xm_mo%nrow+i))+omega*dia)/aa
      enddo
    enddo
  end subroutine mat_dense_new_mo_precond


  subroutine mat_dense_new_complex_precond(nocc,omega,gammma,Eorb_final,Xp_MO,Xm_mo,Xpi_mo,xmi_mo)
    implicit none
    integer, intent(in) :: nocc
    real(realk), intent(in) :: omega,  Eorb_final(:),gammma
    type(Matrix), intent(inout) :: Xp_MO,xm_mo,xpi_mo,xmi_mo
    type(matrix)           :: xp,xm,xpi,xmi
    real(realk) :: aa,a,b,c,d,ab
    integer :: i,j,ndim

    ndim=Xp_MO%nrow 
    call mat_dense_init(xp,ndim,ndim)
    call mat_dense_init(xm,ndim,ndim)
    call mat_dense_init(xpi,ndim,ndim)
    call mat_dense_init(xmi,ndim,ndim)
    call mat_dense_assign(xp,xp_mo)
    call mat_dense_assign(xm,xm_mo)
    call mat_dense_assign(xpi,xpi_mo)
    call mat_dense_assign(xmi,xmi_mo)
    
    
    
    !modify lower left block (ai block)
    do j = 1,nocc  !columns
      do i = nocc+1,Xp_MO%nrow  !rows
        !dia = E[2]_dia + omega*S[2]_dia
        !NOTE TO SELF: WHY the 2.d0s and not 1E0_realk????
       ab=(Eorb_final(i) - Eorb_final(j))*(Eorb_final(i)-Eorb_final(j))-(omega*omega-gammma*gammma) 
        aa=2E0_realk*ab*ab+8E0_realk*(omega*omega*gammma*gammma)
        !dia = 2.0E0_realk*(Eorb_final(i) - Eorb_final(j)) + omega * 2E0_realk
        !dia=Xp_mo%elms((j-1)*Xp_mo%nrow+i)
        a=(Eorb_final(i) - Eorb_final(j))*ab
        b=-omega*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
        &-(omega*omega+gammma*gammma))
        c=gammma*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
        &+(omega*omega+gammma*gammma))
        d=-2E0_realk*omega*gammma*(Eorb_final(i) - Eorb_final(j))
        
        xp%elms((j-1)*Xp_MO%nrow+i) = (a*(Xp_mo%elms((j-1)*Xp_mo%nrow+i))+b*Xm_mo%elms((j-1)*Xm_mo%nrow+i)&
        &-d*(Xpi_mo%elms((j-1)*Xpi_mo%nrow+i))-c*Xmi_mo%elms((j-1)*Xmi_mo%nrow+i)) /aa
        xm%elms((j-1)*Xm_MO%nrow+i) = (b*(Xp_mo%elms((j-1)*Xp_mo%nrow+i))+a*Xm_mo%elms((j-1)*Xm_mo%nrow+i)&
        &-c*(Xpi_mo%elms((j-1)*Xpi_mo%nrow+i))-d*Xmi_mo%elms((j-1)*Xmi_mo%nrow+i)) /aa
        xpi%elms((j-1)*Xpi_MO%nrow+i) = (d*(Xp_mo%elms((j-1)*Xp_mo%nrow+i))+c*Xm_mo%elms((j-1)*Xm_mo%nrow+i)&
        &+a*(Xpi_mo%elms((j-1)*Xpi_mo%nrow+i))+b*Xmi_mo%elms((j-1)*Xmi_mo%nrow+i)) /aa
        xmi%elms((j-1)*Xpi_MO%nrow+i) = (c*(Xp_mo%elms((j-1)*Xp_mo%nrow+i))+d*Xm_mo%elms((j-1)*Xm_mo%nrow+i)&
        &+b*(Xpi_mo%elms((j-1)*Xpi_mo%nrow+i))+a*Xmi_mo%elms((j-1)*Xmi_mo%nrow+i)) /aa
      enddo
    enddo
    !modify upper right block (ia block) 
    do j = nocc+1,Xp_MO%ncol  !columns
      do i = 1,nocc          !rows
       ab=(Eorb_final(i) - Eorb_final(j))*(Eorb_final(i)-Eorb_final(j))-(omega*omega-gammma*gammma) 
        aa=2E0_realk*ab*ab+8E0_realk*(omega*omega*gammma*gammma)
        !dia = 2.0E0_realk*(Eorb_final(i) - Eorb_final(j)) + omega * 2E0_realk
        !dia=Xp_mo%elms((j-1)*Xp_mo%nrow+i)
        a=(Eorb_final(j) - Eorb_final(i))*ab
        b=omega*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
        &-(omega*omega+gammma*gammma))
        c=-gammma*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
        &+(omega*omega+gammma*gammma))
        d=-2E0_realk*omega*gammma*(Eorb_final(j) - Eorb_final(i))
        
        xp%elms((j-1)*Xp_MO%nrow+i) = (a*(Xp_mo%elms((j-1)*Xp_mo%nrow+i))+b*Xm_mo%elms((j-1)*Xm_mo%nrow+i)&
        &-d*(Xpi_mo%elms((j-1)*Xpi_mo%nrow+i))-c*Xmi_mo%elms((j-1)*Xmi_mo%nrow+i)) /aa
        xm%elms((j-1)*Xm_MO%nrow+i) = (b*(Xp_mo%elms((j-1)*Xp_mo%nrow+i))+a*Xm_mo%elms((j-1)*Xm_mo%nrow+i)&
        &-c*(Xpi_mo%elms((j-1)*Xpi_mo%nrow+i))-d*Xmi_mo%elms((j-1)*Xmi_mo%nrow+i)) /aa
        xpi%elms((j-1)*Xpi_MO%nrow+i) = (d*(Xp_mo%elms((j-1)*Xp_mo%nrow+i))+c*Xm_mo%elms((j-1)*Xm_mo%nrow+i)&
        &+a*(Xpi_mo%elms((j-1)*Xpi_mo%nrow+i))+b*Xmi_mo%elms((j-1)*Xmi_mo%nrow+i)) /aa
        xmi%elms((j-1)*Xpi_MO%nrow+i) = (c*(Xp_mo%elms((j-1)*Xp_mo%nrow+i))+d*Xm_mo%elms((j-1)*Xm_mo%nrow+i)&
        &+b*(Xpi_mo%elms((j-1)*Xpi_mo%nrow+i))+a*Xmi_mo%elms((j-1)*Xmi_mo%nrow+i)) /aa
      enddo
    enddo
    call mat_dense_assign(xp_mo,xp)
    call mat_dense_assign(xm_mo,xm)
    call mat_dense_assign(xpi_mo,xpi)
    call mat_dense_assign(xmi_mo,xmi)
    call mat_dense_free(xp)
    call mat_dense_free(xm)
    call mat_dense_free(xpi)
    call mat_dense_free(xmi)
  end subroutine mat_dense_new_complex_precond

!> \brief See mat_mo_precond_complex in mat-operations.f90
  subroutine mat_dense_mo_precond_complex(nocc,Eorb_final,X_MO)
    implicit none
    integer, intent(in) :: nocc
    real(realk), intent(in) :: Eorb_final(:)
    type(Matrix), intent(inout) :: X_MO
    real(realk) :: dia
    integer :: i,j

    !modify lower left block (ai block)
    do j = 1,nocc  !columns
      do i = nocc+1,X_MO%nrow  !rows
        !dia = E[2]_dia + omega*S[2]_dia
        !NOTE TO SELF: WHY the 2.d0s and not 1E0_realk????
        dia = 2.0E0_realk*(Eorb_final(i) - Eorb_final(j))
        X_MO%elms((j-1)*X_MO%nrow+i) = X_MO%elms((j-1)*X_MO%nrow+i) /dia
      enddo
    enddo
    !modify upper right block (ia block) 
    do j = nocc+1,X_MO%ncol  !columns
      do i = 1,nocc          !rows
        !dia = E[2]_dia
        dia = 2.0E0_realk*(Eorb_final(j) - Eorb_final(i))
        X_MO%elms((j-1)*X_MO%nrow+i) = X_MO%elms((j-1)*X_MO%nrow+i)/dia
      enddo
    enddo

  end subroutine mat_dense_mo_precond_complex

!> \brief See mat_ao_precond in mat-operations.f90
  subroutine mat_dense_ao_precond(symmetry,omega,FUP,FUQ,DU,X_AO)
    implicit none
    integer, intent(in) :: symmetry
    real(realk), intent(in) :: omega
    type(Matrix), intent(in) :: FUP, FUQ, DU
    type(Matrix), intent(inout) :: X_AO
    real(realk) :: denom
    integer :: ndim,i,j

    ndim = FUP%nrow
    if (symmetry == 1 .or. symmetry == 2) then !Symmetric or antisymmetric X_AO
       do j = 1,ndim   !columns
         do i = 1,ndim  !rows
           denom = FUQ%elms((j-1)*ndim+j) + FUQ%elms((i-1)*ndim+i)  &
                &- FUP%elms((i-1)*ndim+i) - FUP%elms((j-1)*ndim+j)  &  
                &- omega
           if (dabs(denom) > 1.0E-8_realk) then
              X_AO%elms((j-1)*ndim+i) = X_AO%elms((j-1)*ndim+i)/(denom)
           endif
         enddo
       enddo
       !DEBUG:
       !if (symmetry == 1) then
       !   do j = 1, ndim
       !      do i = 1, j
       !         err = X_AO%elms((j-1)*ndim+i) - X_AO%elms((i-1)*ndim+j)
       !         if (ABS(err) > 1.0E-4_realk) then
       !            print *, 'Warning: Matrix is not symmetric, err =', err 
       !            !CALL LSQUIT(
       !         endif
       !      enddo
       !   enddo
       !else
       !   do j = 1, ndim
       !      do i = 1, j
       !         err = X_AO%elms((j-1)*ndim+i) + X_AO%elms((i-1)*ndim+j)
       !         if (ABS(err) > 1.0E-4_realk) then
       !            print *, 'Warning: Matrix is not antisymmetric, err =', err 
       !            !CALL LSQUIT(
       !         endif
       !      enddo
       !   enddo
       !endif
       !end DEBUG
    else  !X_AO not symmetric in any way
       do j = 1,ndim   !columns
         do i = 1,ndim  !rows
           denom = FUQ%elms((j-1)*ndim+j) + FUQ%elms((i-1)*ndim+i)  &
                &- FUP%elms((i-1)*ndim+i) - FUP%elms((j-1)*ndim+j)  &
                &- omega*(DU%elms((j-1)*ndim+j) - DU%elms((i-1)*ndim+i))
           if (abs(denom) < 1.0E-1_realk) then !Do not divide by too small elements
              denom = denom*1.0E-1_realk/(abs(denom)) !Keep the sign on denominator
           endif
           X_AO%elms((j-1)*ndim+i) = X_AO%elms((j-1)*ndim+i)/(denom)
         enddo
       enddo
    endif
  end subroutine mat_dense_ao_precond

!  !Find lowest hessian diagonal element
!  subroutine mat_dense_hessian_diag(FUP,FUQ,alfa,beta)
!    implicit none
!    type(Matrix), intent(in) :: FUP, FUQ
!    integer, intent(out)     :: alfa, beta
!    real(realk)              :: diag, diag_old
!    integer :: ndim,i,j
!
!    ndim = FUP%nrow
!    diag_old =  FUQ%elms((2-1)*ndim+2) + FUQ%elms((1-1)*ndim+1)  & !1st elm is initial guess
!             &- FUP%elms((1-1)*ndim+1) - FUP%elms((2-1)*ndim+2) 
!    alfa = 1 ; beta = 1
!    do j = 1,ndim   
!      do i = 1, j-1  
!        diag =  FUQ%elms((j-1)*ndim+j) + FUQ%elms((i-1)*ndim+i)  &
!             &- FUP%elms((i-1)*ndim+i) - FUP%elms((j-1)*ndim+j) 
!             !print *, 'diag', diag 
!        if (diag < diag_old) then
!           alfa = j ; beta = i
!           !print *, 'new diag chosen', diag, 'alfa, beta', alfa, beta
!           diag_old = diag
!        endif
!      enddo
!    enddo
!  end subroutine mat_dense_hessian_diag

  !Find all hessian diagonal elementss and alpha/beta indices
  !and sort them in increasing order. 
  !SONIA
!  subroutine mat_dense_hessian_diags(FUP,FUQ,indici)
!    implicit none
!    type(Matrix), intent(in) :: FUP, FUQ
!    integer, intent(out)     :: indici(FUP%nrow*(FUP%nrow-1)/2,2)
!    real(realk)              :: diag_old
!    integer                  :: ndim,i,j,xv(2),n
!    real(realk)              :: x
!    real(realk), allocatable :: diag(:)
!
!    allocate(diag(FUP%nrow*(FUP%nrow-1)/2))
!    ndim = FUP%nrow
!    n=0
!    do j = 1,ndim   
!      do i = 1, j-1  
!        n = n + 1
!        diag(n) =  FUQ%elms((j-1)*ndim+j) + FUQ%elms((i-1)*ndim+i)  &
!                &- FUP%elms((i-1)*ndim+i) - FUP%elms((j-1)*ndim+j) 
!        indici(n,1) = j
!        indici(n,2) = i
!      enddo
!    enddo
!
!    !print *, 'INDICI in mat_dense_hessian_diags before sorting'
!    !do i=1,n
!    !  print *, 'Indici alpha,beta ', indici(i,1), indici(i,2)
!    !end do
!    !do i=1,n
!    !  print *, 'Diag elements ', i, diag(i)
!    !end do
!
!    do i = 1,n
!       do j = 1,n-1
!         if (diag(j) > diag(j+1)) then
!           x = diag(j+1)
!           diag(j+1) = diag(j)
!           diag(j) = x
!           xv = indici(j+1,:)
!           indici(j+1,:) = indici(j,:)
!           indici(j,:) = xv
!         endif
!       enddo
!    enddo
!    deallocate(diag)
!    !print *, 'INDICI in mat_dense_hessian_diags AFTER sorting'
!    !do i=1,n
!    !  print *, 'Indici alpha,beta ', indici(i,1), indici(i,2)
!    !end do
!    !do i=1,n
!    !  print *, 'Diag elements ', i, diag(i)
!    !end do
!  end subroutine mat_dense_hessian_diags
!
!  subroutine mat_dense_hes_ao_precond(FPPU,FQQU,SPPU,SQQU,X_AO)
!    implicit none
!    !integer, intent(in) :: symmetry
!    !real(realk), intent(in) :: omega
!    type(Matrix), intent(in) :: FPPU, FQQU, SPPU, SQQU
!    type(Matrix)             :: hes
!    type(Matrix), intent(inout) :: X_AO !On input, X_AO has the value of the residual
!    real(realk) :: value
!    real(realk), allocatable :: hes_full(:,:), res(:,:), x_full(:)
!    integer, allocatable :: IPIV(:)
!!    real(realk), parameter :: alpha = 4.20E0_realk
!    integer :: ndim, hesdim, mu, nu, alfa, beta, INFO
!    integer :: index1, index2, col_index
!
!    ndim = FPPU%nrow
!    hesdim = ndim*ndim
!    call mat_dense_init(hes,hesdim,hesdim)
!    ALLOCATE(hes_full(hesdim,hesdim),res(hesdim,1),IPIV(hesdim), x_full(hesdim))
!
!    do mu = 1, ndim
!       do nu = 1, ndim
!          index1 = ndim*(nu-1) + mu
!          do alfa = 1, ndim
!             do beta = 1, ndim
!             index2 = ndim*(beta-1) + alfa
!             value = SPPU%elms(ndim*(alfa-1)+mu)*FQQU%elms(ndim*(nu-1)+beta) &
!                 & + FQQU%elms(ndim*(alfa-1)+mu)*SPPU%elms(ndim*(nu-1)+beta) &
!                 & - SQQU%elms(ndim*(alfa-1)+mu)*FPPU%elms(ndim*(nu-1)+beta) &
!                 & - FPPU%elms(ndim*(alfa-1)+mu)*SQQU%elms(ndim*(nu-1)+beta)
!             col_index = hesdim*(index2-1)+index1
!             hes%elms(col_index) = value
!             enddo
!          enddo
!       enddo
!    enddo
!    write (mat_lu,*) 'Hessian from mat_dense_hes_ao_precond:' 
!    call mat_dense_print(hes, 1, hes%nrow, 1, hes%ncol, mat_lu)
!
!    !Solve linear equation hes*X_AO = R:
!    
!    !Convert framework:
!    call mat_dense_to_full(hes, 1.0E0_realk, hes_full)
!!call LS_OUTPUT (hes_full,1,hesdim,1,hesdim,hesdim,hesdim,1,mat_lu)
!!    write (mat_lu,*) 'Incoming residual:' 
!!    call mat_dense_print(X_AO, 1, X_AO%nrow, 1, X_AO%ncol, mat_lu)
!
!    call mat_dense_to_full(X_AO, 1.0E0_realk, res)
!!    write (mat_lu,*) 'Residual converted to full:'
!!call LS_OUTPUT (res,1,hesdim,1,1,hesdim,1,1,mat_lu)
!
!    call LIN_EQ(hes_full,res,hesdim,x_full,INFO)
!
!    print *, "Did solution of linear equations go okay?", INFO
!
!!CALL LSQUIT(
!!        DGETRS(TRANS,  N,  NRHS, A,       LDA,   IPIV, B,    LDB,   INFO)
!!   call DGETRS('N', hesdim, 1, hes_full, hesdim, IPIV, res, hesdim, INFO)
!    !call DGESV(ndim*ndim, 1, hes_full, ndim*ndim, IPIV, res, ndim*ndim, INFO)
!!CALL LSQUIT(
!!    print *, "Did solution of linear equations go okay?", INFO
!
!!CALL LSQUIT(
!    !Convert framework:
!    call mat_dense_set_from_full(x_full,1.0E0_realk,X_AO)
!CALL LSQUIT(
!  call mat_dense_free(hes)
!  DEALLOCATE(res)
!  DEALLOCATE(hes_full)
!  DEALLOCATE(IPIV)
!  DEALLOCATE(x_full)
!  end subroutine mat_dense_hes_ao_precond

!########################## Stinne lineq solver!!! ########################################
!1. Solving linear equations A*x = b
!   by QR decomposition/Gram-Schmidt and
!   back-substitution.
!
!   Input: A(n,n), real double precision matrix
!          b(n), real double precision vector
!   Output: x(n), solution vector, double precision
!
!   subroutine LIN_EQ(A,b,n,x,INFO)
!   implicit none
!   
!        integer, intent(in)     :: n
!        integer, intent(out)    :: INFO
!        real(realk)             :: val
!        real(realk), intent(in) :: A(n,n), b(n)
!        real(realk), allocatable, dimension(:,:) :: Q, R, QR, A_change
!        real(realk), allocatable, dimension(:)   :: test1, test2 
!        real(realk), intent(out) :: x(n)
!        integer :: i
!
!   INFO = 0
!
!   ALLOCATE(QR(n,n), Q(n,n), R(n,n), test1(n), test2(n), A_change(n,n))
!   
!   A_change = A
!!print *, 'test1'   
!   CALL QR_GS(A_change, n, Q, R)
!!print *, 'test2'
!!Test:
!   QR = matmul(Q,R)
!!print *, 'test3'
!write (mat_lu,*) 'Q*R:'
!call LS_OUTPUT (QR,1,n,1,n,n,n,1,mat_lu)
!write (mat_lu,*) 'Q:'
!call LS_OUTPUT (Q,1,n,1,n,n,n,1,mat_lu)
!write (mat_lu,*) 'R:'
!call LS_OUTPUT (R,1,n,1,n,n,n,1,mat_lu)
!
!!print *, 'test4'    
!   CALL QR_BACK(Q, R, n, b, x)
!write (mat_lu,*) 'Incoming b in lin_eq:'
!call LS_OUTPUT (b,1,n,1,1,n,1,1,mat_lu)
!
!   !Test if solution of linear equations went okay:
!      test1 = matmul(A,x)
!write (mat_lu,*) 'A*x:'
!call LS_OUTPUT (test1,1,n,1,1,n,1,1,mat_lu)
!
!      test2 = b - test1
!      val = 0.0E0_realk
!      do i = 1, n
!         val = val + test2(i)**2
!      enddo                                                                 
!      val=SQRT(val)
!      print *, "Norm square of b - x =", val
!      if (val > 1.0E-8_realk) then
!         INFO = -1
!      endif
!   DEALLOCATE(Q,R,QR,test1,test2,A_change)
!   end subroutine LIN_EQ
!   
!   !Subroutine for QR decomposition by a modified Gram-Schmidt method
!   !THIS SUBROUTINE CHANGES A!!!!!!!!!!!!!!!!!!!!!!
!   subroutine QR_GS(A, n, Q, R)
!   Implicit none
!   
!        integer, intent(in) :: n
!        integer             :: i, j, k, l
!        real(realk), intent(inout) :: A(n,n)
!        real(realk), intent(out):: Q(n,n), R(n,n)
!   
!   do i = 1, n
!      R(i,i)=SQRT(DOT_PRODUCT(A(1:n,i),A(1:n,i)))
!      do j = 1, n
!         Q(j,i) = A(j,i)/R(i,i)
!      enddo
!      do j = i+1, n
!         R(i,j)=DOT_PRODUCT(Q(1:n,i), A(1:n,j))
!         do k = 1, n
!            A(k,j) = A(k,j)-Q(k,i)*R(i,j)
!         enddo
!      enddo
!   enddo   
!   end subroutine QR_GS
!   
!   !Subroutine for solving linear equations Ax=b => QRx=b by QR back-substitution.
!   subroutine QR_BACK(Q, R, n, b, x)
!   implicit none
!   
!        integer, intent(in)     :: n
!        real(realk), intent(in) :: Q(n,n), R(n,n), b(n)
!        real(realk),allocatable :: Q_T(:,:), b_ny(:)
!        real(realk), intent(out):: x(n)
!        integer :: i, k
!   
!   ALLOCATE(Q_T(n,n), b_ny(n))
!   Q_T = TRANSPOSE(Q)
!   b_ny = MATMUL(Q_T, b)
!
!write (mat_lu,*) 'b_ny = Q_T*b:'
!call LS_OUTPUT (b_ny,1,n,1,1,n,1,1,mat_lu)
!   
!   x = 0.0
!   
!   do i = n, 1, -1
!      do k = i+1, n
!         x(i) = x(i) - R(i,k)*x(k)
!      enddo
!      x(i) = (x(i) + b_ny(i))/R(i,i)
!   enddo
!   DEALLOCATE(Q_T)
!   DEALLOCATE(b_ny)
!   end subroutine QR_BACK
!########################## end of Stinne lineq solver!!! ########################################

!> \brief See mat_identity in mat-operations.f90
  subroutine mat_dense_identity(Id)
    implicit none
    type(matrix), intent(inout) :: Id
    integer :: i,j

    do j = 1,Id%ncol
      do i = 1,Id%nrow
        if (i == j) then
          Id%elms((j-1)*Id%nrow+i) = 1.0E0_realk
        else
          Id%elms((j-1)*Id%nrow+i) = 0.0E0_realk
        endif
      enddo
    enddo
  end subroutine mat_dense_identity

!> \brief See mat_create_elm in mat-operations.f90
  subroutine mat_dense_create_elm(i,j,val,A)
    implicit none
    integer, intent(in) :: i,j
    real(realk), intent(in) :: val
    type(Matrix), intent(inout) :: A

    A%elms((j-1)*A%nrow+i) = val

  end subroutine mat_dense_create_elm

!> \brief See mat_get_elm in mat-operations.f90
!=======================================================================
  subroutine mat_dense_get_elm (A, r, c, elm)

  implicit none

  type(matrix), intent(in) :: A
  integer, intent(in)      :: r, c
  real(realk), intent(out) :: elm

  elm = A%elms(A%nrow*(c-1)+r)

  end subroutine mat_dense_get_elm
!=======================================================================

!> \brief See mat_create_block in mat-operations.f90
  subroutine mat_dense_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in)         :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(in)     :: fullmat(fullrow,fullcol)
    type(Matrix), intent(inout) :: A
    integer                     :: i, j

    do j = insertcol, insertcol+fullcol-1
       do i = insertrow, insertrow+fullrow-1
          A%elms((j-1)*A%nrow+i) = fullmat(i-insertrow+1,j-insertcol+1)
       enddo
    enddo

  end subroutine mat_dense_create_block

!> \brief See mat_add_block in mat-operations.f90
  subroutine mat_dense_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in)         :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(in)     :: fullmat(fullrow,fullcol)
    type(Matrix), intent(inout) :: A
    integer                     :: i, j

    do j = insertcol, insertcol+fullcol-1
       do i = insertrow, insertrow+fullrow-1
          A%elms((j-1)*A%nrow+i) = A%elms((j-1)*A%nrow+i) + fullmat(i-insertrow+1,j-insertcol+1)
       enddo
    enddo

  end subroutine mat_dense_add_block

!> \brief See mat_retrieve_block in mat-operations.f90
  subroutine mat_dense_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in)         :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(inout)     :: fullmat(fullrow,fullcol)
    type(Matrix), intent(in) :: A
    integer                     :: i, j

    do j = insertcol, insertcol+fullcol-1
       do i = insertrow, insertrow+fullrow-1
          fullmat(i-insertrow+1,j-insertcol+1) = A%elms((j-1)*A%nrow+i)
       enddo
    enddo

  end subroutine mat_dense_retrieve_block

!> \brief See mat_scal in mat-operations.f90
  subroutine mat_dense_scal(alpha,A)
    implicit none
    real(realk), intent(in) :: alpha
    type(Matrix), intent(inout) :: A
    integer :: i

!   do i = 1,A%nrow*A%ncol
!     A%elms(i) = A%elms(i) * alpha
!   enddo
    i = A%nrow*A%ncol
    call dscal(i,alpha,A%elms,1)

  end subroutine mat_dense_scal

!> \brief See mat_scal_dia in mat-operations.f90
  subroutine mat_dense_scal_dia(alpha,A)
    implicit none
    real(realk), intent(in) :: alpha
    type(Matrix), intent(inout) :: A

!   do i = 1,A%nrow
!     A%elms((i-1)*A%nrow+i) = A%elms((i-1)*A%nrow+i) * alpha
!   enddo

    call dscal(A%nrow,alpha,A%elms,A%nrow+1)

  end subroutine mat_dense_scal_dia

!> \brief See mat_scal_dia in mat-operations.f90
  subroutine mat_dense_scal_dia_vec(alpha,A,ndim)
    implicit none
    integer, intent(in)     :: ndim
    real(realk), intent(in) :: alpha(ndim)
    type(Matrix), intent(inout) :: A
    integer :: i
    do i = 1,A%nrow
       A%elms((i-1)*A%nrow+i) = A%elms((i-1)*A%nrow+i) * alpha(i)
    enddo
  end subroutine mat_dense_scal_dia_vec

!> \brief See mat_zero in mat-operations.f90
  subroutine mat_dense_zero(A)
    implicit none
    type(Matrix), intent(inout) :: A
    integer :: i

    do i = 1,A%nrow*A%ncol
      A%elms(i) = 0.0E0_realk
    enddo
 
  end subroutine mat_dense_zero

!> \brief See mat_zerohalf in mat-operations.f90
  subroutine mat_dense_zerohalf(part,A)
    implicit none
    character(len=2), intent(in) :: part
    type(Matrix), intent(inout) :: A
    integer :: i,j

    if (part == 'ut' .or. part == 'UT') then
      !set the upper triangle to zero - diagonal is kept
      do j = 2,A%ncol
        do i = 1,j-1
          A%elms((j-1)*A%nrow+i) = 0.0E0_realk
        enddo
      enddo
    else
      !set the lower triangle to zero - diagonal is also zeroed
      do j = 1,A%ncol
        do i = j,A%nrow
          A%elms((j-1)*A%nrow+i) = 0.0E0_realk
        enddo
      enddo
    endif

  end subroutine mat_dense_zerohalf

!> \brief See mat_write_to_disk in mat-operations.f90
  subroutine mat_dense_write_to_disk(iunit,A)
  !===========================================================================
  ! Writes matrix A to file
  !           INPUT: iunit (integer) - The unit number
  !                  A (matrix) - The matrix that should be written to disk
    implicit none
    integer, intent(in) :: iunit
    type(Matrix), intent(in) :: A
!    character(len=1), allocatable :: cmpd(:)
    integer :: i
    integer(kind=long) :: ncol,nrow 

    if (.not.ASSOCIATED(A%elms)) CALL LSQUIT( 'A in mat_dense_WRITE_TO_DISK non-existant',-1)

 !   INQUIRE(iunit,NAME=filename)     !Find the filename of the file
 !   print*,filename
 !   CLOSE(iunit)
!TODO: Check this
!    OPEN(iunit, FILE=filename, STATUS="old", FORM="unformatted", POSITION="append" )
!    OPEN(iunit, FILE=filename, STATUS="old", FORM="unformatted")
 !   REWIND(iunit)
!   if (.not.cfg_compress) then

    nrow = A%Nrow
    ncol = A%Ncol
    write(iunit) Nrow, Ncol
    WRITE(iunit)(A%elms(I),I=1,A%nrow*A%ncol)
    !WRITE(iunit) A%Nrow, A%Ncol
    !WRITE(iunit) A%elms

!    else
!       slen = A%Nrow*A%Ncol*8
!       dlen = (slen*1.001) + 24
!       allocate(cmpd(dlen))        

!       call compress(A%elms,slen,cmpd,dlen)
!       WRITE(iunit) A%Nrow, A%Ncol, dlen
!       WRITE(iunit) (cmpd(i),i=1,dlen)
!       
!       deallocate(cmpd)
!   endif


 !   ENDFILE iunit
  end subroutine mat_dense_WRITE_TO_DISK

!> \brief See mat_read_from_disk in mat-operations.f90
  subroutine mat_dense_read_from_disk(iunit,A)
  !===============================================================================
  ! Reads a matrix from disk
  !            INPUT: iunit (integer) - The unit number
  !           OUTPUT: A (matrix) - The matrix wanted.
    implicit none
    integer, intent(in) :: iunit
    type(Matrix), intent(inout) :: A
    integer :: i
    integer(kind=long) :: Nrow,Ncol
!   logical :: cfg_compress = .false.

!    REWIND iunit

!   if (.not.cfg_compress) then
    READ(iunit) Nrow, Ncol
    IF(Nrow.EQ.A%Nrow.AND.Ncol.EQ.A%Ncol)THEN
       READ(iunit)(A%elms(I),I=1,A%nrow*A%ncol)
    ELSE
       print*,'Error in reading matrix from disk. Dimension mismatch'
       print*,'Dimensions of the matrix on Disk  :',Nrow,Ncol
       print*,'Allocated Dimensions of the matrix:',A%Nrow,A%Ncol
       CALL LSQUIT('Error in reading matrix from disk. Dimension mismatch',-1)
    ENDIF
!    READ(iunit) A%Nrow, A%Ncol         
!    READ(iunit) A%elms
!   else
!       READ(iunit) A%Nrow, A%Ncol, slen
!       allocate(cmpd(slen))
!       READ(iunit) (cmpd(i),i=1,slen)
!       dlen=A%Nrow*A%Ncol*8
!       call uncompress(cmpd,slen,A%elms,dlen)
!       deallocate(cmpd)
!   endif

  end subroutine mat_dense_READ_FROM_DISK

!> \brief See mat_write_to_disk2 in mat-operations.f90
  subroutine mat_dense_write_to_disk2(iunit,A)
  !Hack routine - see debug_convert_density in debug.f90
    implicit none
    integer, intent(in) :: iunit
    type(Matrix), intent(in) :: A
!    character(len=20) :: filename
!    character(len=1), allocatable :: cmpd(:)
    integer :: i

    if (.not.ASSOCIATED(A%elms)) CALL LSQUIT( 'A in mat_dense_WRITE_TO_DISK non-existant',-1)

 !   INQUIRE(iunit,NAME=filename)     !Find the filename of the file
 !   print*,filename
 !   CLOSE(iunit)
!TODO: Check this
!    OPEN(iunit, FILE=filename, STATUS="old", FORM="unformatted", POSITION="append" )
!    OPEN(iunit, FILE=filename, STATUS="old", FORM="unformatted")
 !   REWIND(iunit)
!   if (.not.cfg_compress) then

    WRITE(iunit,*) A%Nrow, A%Ncol
    !WRITE(iunit,*) A%elms
    WRITE(iunit,*) (A%elms(I),I=1,A%nrow*A%ncol)

!    else
!       slen = A%Nrow*A%Ncol*8
!       dlen = (slen*1.001) + 24
!       allocate(cmpd(dlen))        

!       call compress(A%elms,slen,cmpd,dlen)
!       WRITE(iunit) A%Nrow, A%Ncol, dlen
!       WRITE(iunit) (cmpd(i),i=1,dlen)
!       
!       deallocate(cmpd)
!   endif


 !   ENDFILE iunit
  end subroutine mat_dense_WRITE_TO_DISK2

!> \brief See mat_from_disk2 in mat-operations.f90
  subroutine mat_dense_read_from_disk2(iunit,A)
  !===============================================================================
  ! Reads a matrix from disk
  !            INPUT: iunit (integer) - The unit number
  !           OUTPUT: A (matrix) - The matrix wanted.
  !Hack routine - see debug_convert_density in debug.f90
    implicit none
    integer, intent(in) :: iunit
    type(Matrix), intent(inout) :: A
!   logical :: cfg_compress = .false.

!    REWIND iunit

!   if (.not.cfg_compress) then
    READ(iunit,*) A%Nrow, A%Ncol         
    READ(iunit,*) A%elms
!   else
!       READ(iunit) A%Nrow, A%Ncol, slen
!       allocate(cmpd(slen))
!       READ(iunit) (cmpd(i),i=1,slen)
!       dlen=A%Nrow*A%Ncol*8
!       call uncompress(cmpd,slen,A%elms,dlen)
!       deallocate(cmpd)
!   endif

  end subroutine mat_dense_READ_FROM_DISK2

!> \brief See mat_vec_to_mat in mat-operations.f90
    subroutine mat_dense_vec_to_mat(symmetry, VEC, MAT)
    !Change a vector to a matrix. The matrix is symmetric or antisymmetric.
    implicit none
                                           
         character, intent(in)     :: symmetry                                   
         integer                   :: n, m, i, k
         TYPE(matrix), intent(in)  :: VEC
         TYPE(matrix)              :: MAT   !Intent(out)
                                       
 
    if (symmetry == 's' .OR. symmetry  == 'S') then
       do n = 1, MAT%nrow
          do m = 1, n
             i = n*(n-1)/2 + m
             MAT%elms((m-1)*MAT%nrow + n) = VEC%elms(i)
             MAT%elms((n-1)*MAT%nrow + m) = VEC%elms(i)
          enddo
       enddo
    endif

 !      do n = 1, dm_mat
 !         do m = 1, n
 !            i = n*(n-1)/2 + m
 !            MAT(n,m) = VEC(i)
 !            MAT(m,n) = -VEC(i)
 !         enddo
 !      enddo
    k = 0
    if (symmetry == 'a' .OR. symmetry == 'A') then
       do n = 1, MAT%nrow
          do m = 1, n
             if (n == m) then
                MAT%elms((n-1)*MAT%nrow + m) = 0.0E0_realk !Diagonal element are 0 (antisym. matrix)
             else
                k = k + 1
                MAT%elms((m-1)*MAT%nrow + n) = VEC%elms(k)
                MAT%elms((n-1)*MAT%nrow + m) = -VEC%elms(k)
             endif
          enddo
       enddo
    endif
    end subroutine mat_dense_vec_to_mat
 
!> \brief See mat_to_vec in mat-operations.f90
    subroutine mat_dense_mat_to_vec(symmetry, MAT, VEC)
    implicit none
 
         character, intent(in)    :: symmetry
         integer                  :: n, m, i, k
         TYPE(matrix)             :: VEC        !Intent(out)
         TYPE(matrix), intent(in) :: MAT

    if (symmetry == 's' .OR. symmetry == 'S') then
       do n = 1, MAT%nrow
          do m = 1, n
             i = n*(n-1)/2 + m
             VEC%elms(i) = MAT%elms((m-1)*MAT%nrow + n)
          enddo
       enddo
    endif

    k = 0
    if (symmetry == 'a' .OR. symmetry == 'A') then
       do n = 1, MAT%nrow
          do m = 1, n
             if (n /= m) then
                k = k + 1
                VEC%elms(k) = MAT%elms((m-1)*MAT%nrow + n)
             endif
          enddo
       enddo
    endif

    end subroutine mat_dense_mat_to_vec

    SUBROUTINE mat_dense_dposv(A,b,lupri)
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
         call lsquit('mat_dense_dposv, Reason: Wrong dim3',lupri)
      endif
      call mem_alloc(Af,dim1,dim2)
      call mem_alloc(bf,dim3,dim4)
      call mat_dense_to_full(A,1D0,Af)
      call mat_dense_to_full(b,1D0,bf)

      call dposv('U',dim1,dim4,Af,dim1,bf,dim3,info)
      
      if(info .ne. 0) then
         call lsquit('mat_dense_dposv, Reason: info not 0',lupri)
      endif
      
      call mat_dense_set_from_full(Af,1D0,A)
      call mat_dense_set_from_full(bf,1D0,b)
      call mem_dealloc(Af)
      call mem_dealloc(bf)
      
    END SUBROUTINE mat_dense_dposv

    !> \brief creates the inverse matrix of type(matrix). mat_inv
    !> \author T. Kjærgaard
    !> \date 2012
    !> \param a The type(matrix) that should be inversed
    !> \param chol The type(matrix) that contains cholesky factors (from mat_chol)
    !> \param c The inverse output type(matrix).
    SUBROUTINE mat_dense_inv(A, A_inv) 
      implicit none
      TYPE(Matrix),intent(in)     :: A
      TYPE(Matrix)                :: A_inv !output
      real(realk), pointer   :: work1(:)
      real(realk), pointer   :: A_inv_full(:,:) 
      integer,pointer    :: IPVT(:)
      real(realk)            :: RCOND, dummy(2), tmstart, tmend
      integer                :: IERR, i, j, fulldim, ndim
      fulldim = a%nrow
      call mem_alloc(A_inv_full,fulldim,fulldim) 
      call mem_alloc(work1,fulldim)
      call mem_alloc(IPVT,fulldim)
      !Invert U and Ut:
      IPVT = 0 ; RCOND = 0.0E0_realk  
      call mat_dense_to_full(A,1.0E0_realk,A_inv_full)
      call DGECO(A_inv_full,fulldim,fulldim,IPVT,RCOND,work1)
      call DGEDI(A_inv_full,fulldim,fulldim,IPVT,dummy,work1,01)
      !Convert framework:
      call mat_dense_set_from_full(A_inv_full,1.0E0_realk,A_inv)
      call mem_dealloc(A_inv_full) 
      call mem_dealloc(work1)
      call mem_dealloc(IPVT)
    END SUBROUTINE mat_dense_inv

!> \brief Compute the Cholesky decomposition factors of a positive definite matrix
!> \author T. Kjærgaard
!> \date 2012
!> \param a The type(matrix) that should be decomposed
!> \param b The output type(matrix) that contains the cholesky factors.
      SUBROUTINE mat_dense_chol(a, b) 
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
        call mat_dense_to_full(A,1.0E0_realk,U_full)
           !Set lower half of U = 0:
        do i = 1, fulldim
           do j = 1, i-1
              U_full(i,j) = 0.0E0_realk
           enddo
        enddo
        JPVT(1) = 0
        call dchdc(U_full,fulldim,fulldim,work1,JPVT,0,IERR)           
        call mat_dense_set_from_full(U_full,1.0E0_realk,B)
        call mem_dealloc(U_full) 
        call mem_dealloc(work1)
        call mem_dealloc(JPVT)
      end SUBROUTINE mat_dense_chol

!Routines needed for purification
! - commented out because purification is not documented and no one really knows
! if purification works!

!!> \brief See mat_cholesky in mat-operations.f90
!  subroutine mat_dense_cholesky(A,B)
!  !returns cholesky decomposition of matrix A in B
!    implicit none
!    type(Matrix), intent(in) :: A
!    type(Matrix), intent(inout) :: B
!    integer :: info,i,j
!
!!   call mat_dense_copy(1E0_realk,A,B)
!    call mat_dense_assign(B,A)
!    call dpotrf('L',B%nrow,B%elms,B%nrow,info)
!    if(info.ne. 0) call lsquit( 'exit in mat_dense_cholesky!'
!    call mat_dense_zerohalf('ut',B)
!
!    return
!  end subroutine mat_dense_cholesky
!
!!> \brief See mat_inverse_triang in mat-operations.f90
!  subroutine mat_dense_inverse_triang(A,B)
!  !returns inverse of lower triangular matrix A in B
!    implicit none
!    type(Matrix), intent(in) :: A
!    type(Matrix), intent(inout) :: B
!    integer :: info
!
!!   call mat_dense_copy(1E0_realk,A,B)
!    call mat_dense_assign(B,A)
!    call dtrtri('L','N',B%nrow,B%elms,B%nrow,info)
!    if(info.ne. 0) call lsquit( 'exit in mat_dense_inverse_triang!'
!
!    return
!  end subroutine mat_dense_inverse_triang
!
!!> \brief See mat_simtran in mat-operations.f90
!  subroutine mat_dense_simtran(A,B,transb,C)
!  !returns similarity transformation B.A.B^T in C
!    implicit none
!    character, intent(in) :: transb
!    type(Matrix), intent(in) :: A,B
!    type(Matrix), intent(inout) :: C
!    type(Matrix) :: T
!
!    call mat_dense_init(T,A%nrow,A%ncol)
!
!    if(transb.eq.'t'.or.transb.eq.'T') then
!       call mat_dense_mul(B,A,'T','N',1E0_realk,0E0_realk,T)
!       call mat_dense_mul(T,B,'N','N',1E0_realk,0E0_realk,C)
!    else
!       call mat_dense_mul(B,A,'N','N',1E0_realk,0E0_realk,T)
!       call mat_dense_mul(T,B,'N','T',1E0_realk,0E0_realk,C)
!    endif
!
!    call mat_dense_free(T)
!
!    return
!  end subroutine mat_dense_simtran
!
!!> \brief See mat_gershgorin_minmax in mat-operations.f90
!  subroutine mat_dense_gershgorin_minmax(A,min,max)
!    implicit none
!    type(Matrix), intent(in) :: A
!    real(realk), intent(out) :: min,max
!    integer :: i,j
!    real(realk) :: offsum,tmin,tmax
!
!    do i=1,A%nrow
!       offsum=0E0_realk
!       do j=1,A%ncol
!          if(i.eq.j) cycle
!          offsum=offsum+abs(A%elms(a%nrow*(i-1)+j))
!       enddo
!       tmin=A%elms(a%nrow*(i-1)+i)-offsum
!       tmax=A%elms(a%nrow*(i-1)+i)+offsum
!       if(i.eq. 1) then
!          min=tmin
!          max=tmax
!       else
!          if(min.gt.tmin) min=tmin
!          if(max.lt.tmax) max=tmax
!       endif
!    enddo
!
!    return
!  end subroutine mat_dense_gershgorin_minmax
!
!!> \brief See mat_clean in mat-operations.f90
!  subroutine mat_dense_clean(A,lowcut)
!    !symmetrise and filter small values of matrix
!    implicit none
!    type(Matrix), intent(inout) :: A
!    real(realk), intent(in) :: lowcut
!    integer :: i,j
!    real(realk) :: aij
!
!    do i=1,A%nrow
!       do j=1,i
!          aij=A%elms(a%nrow*(i-1)+j)
!          if(abs(aij).lt.lowcut) aij=0E0_realk
!          A%elms(a%nrow*(i-1)+j)=aij
!          A%elms(a%nrow*(j-1)+i)=aij
!       enddo
!    enddo
!
!    return
!  end subroutine mat_dense_clean
!
!!> \brief See mat_to_minus_one_half in mat-operations.f90
!  subroutine mat_dense_to_minus_one_half(A,B)
!    implicit none
!    type(Matrix), intent(in) :: A
!    type(Matrix), intent(inout) :: B
!    type(Matrix) :: e
!    INTEGER :: lwork, info
!    real(realk), dimension(:), allocatable :: work
!
!    lwork=6*A%nrow
!    allocate(work(lwork))
!    call mat_dense_init(e,A%nrow,1)
!
!!   call mat_dense_copy(1E0_realk,A,B)
!    call mat_dense_assign(B,A)
!    call dsyev('V','U',A%nrow,B%elms,A%nrow,e%elms,work,lwork,info)
!    if(info.ne. 0) call lsquit( 'exit in mat_dense_to_minus_one_half!'
!    call s_to_minus_one_half(B%elms,A%nrow,e%elms,A%nrow)
!
!    deallocate(work)
!    call mat_dense_free(e)
!
!    return
!  end subroutine mat_dense_to_minus_one_half
!
!!----------------------------------------------------------------------
!subroutine s_to_minus_one_half(u,nu,e,n)
!!----------------------------------------------------------------------
!! forms u=s**(-1/2)
!! dimension u(nu,n),e(n)
!! eigenvectors and eigenvalues of s must be provided in u and e,
!! respectively. u and e are overwritten.
!   implicit none
!   integer, intent(in) :: nu,n
!   real(8), intent(inout) :: u(nu,n),e(n)
!   integer :: i,j,k
!   real(8) :: fak
!   real(8), dimension(:), allocatable :: v
!
!   if(n.gt.nu) call LSquit('DIMENSION ERROR in SMH', mat_lu)
!   allocate(v(n))
!
!   do i=1,n
!      if(e(i).gt. 1E-10_realk) then
!        e(i)=1.0E0_realk/sqrt(e(i))
!      else
!        e(i)=0
!      end if
!   end do
!   do 50 j=1,n
!   do 20 i=j,n
!20 v(i)=0
!   do 30 k=1,n
!   fak=e(k)*u(j,k)
!   do 30 i=j,n
!30 v(i)=v(i)+u(i,k)*fak
!   do 50 i=j,n
!50 u(j,i)=v(i)
!   do 60 i=2,n
!   do 60 j=1,i-1
!60 u(i,j)=u(j,i)
!   
!   return
!end subroutine s_to_minus_one_half
!
!> \brief See mat_sum in mat-operations.f90
  function mat_dense_sum(A)
    !returns sum of all elements of matrix
    implicit none
    type(Matrix), intent(in) :: A
    real(realk) :: mat_dense_sum
    integer :: i

    mat_dense_sum=0E0_realk
    do i=1,A%nrow*A%ncol
          mat_dense_sum=mat_dense_sum+a%elms(i)
    enddo

    return
  end function mat_dense_sum

!> \brief See mat_report_sparsity in mat-operations.f90
  subroutine mat_dense_report_sparsity(A,sparsity)
  implicit none
  type(Matrix) :: A
  real(realk)  :: nnz,sparsity
  integer      :: i

       nnz=0E0_realk
       do i=1, A%ncol*A%nrow
         if (abs(A%elms(i)).gt. 1E-9_realk) nnz = nnz+1E0_realk
       enddo

       sparsity = nnz/(A%ncol*A%nrow)
  end subroutine mat_dense_report_sparsity

end module matrix_operations_dense
