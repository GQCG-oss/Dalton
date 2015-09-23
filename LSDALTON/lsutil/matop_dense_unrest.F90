!> @file
!> Contains unrestricted dense matrix module.

!> \brief Contains matrix operation routines for type(matrix) = unrestricted dense.
module matrix_op_unres_dense
  use matrix_module
  use memory_handling
  use precision
public :: mat_unres_dense_init
public :: mat_unres_dense_free
public :: mat_unres_dense_set_from_full
public :: mat_unres_dense_to_full
public :: mat_unres_dense_to_full3d
public :: mat_unres_dense_print
public :: mat_unres_dense_trans
public :: mat_unres_dense_assign
public :: mat_unres_dense_mpicopy
public :: mat_unres_dense_copy
public :: mat_unres_dense_Tr
public :: mat_unres_dense_TrAB 
public :: mat_unres_dense_mul
public :: mat_unres_dense_dmul
public :: mat_unres_dense_extract_diagonal
public :: mat_unres_dense_add
public :: mat_unres_dense_daxpy
public :: mat_unres_dense_dotproduct
public :: mat_unres_dense_sqnorm2
public :: mat_unres_dense_outdia_sqnorm2
public :: mat_unres_dense_dsyev
public :: mat_unres_dense_dsyevx
public :: mat_unres_dense_dsyevx_aux
public :: mat_unres_dense_diag_f
public :: mat_unres_dense_abs_max_elm
public :: mat_unres_dense_max_elm
public :: mat_unres_dense_min_elm
public :: mat_unres_dense_max_diag_elm
public :: mat_unres_dense_dE_dmu
public :: mat_unres_dense_column_norm
public :: mat_unres_dense_section
public :: mat_unres_dense_insert_sectio
public :: mat_unres_dense_section2
public :: mat_unres_dense_precond
public :: mat_unres_dense_mo_precond
public :: mat_unres_dense_new_mo_precond
public :: mat_unres_dense_new_complex_precond
public :: mat_unres_dense_mo_precond_complex
public :: mat_unres_dense_ao_precond
public :: mat_unres_dense_identity
public :: mat_unres_dense_create_elm
public :: mat_unres_dense_get_el
public :: mat_unres_dense_create_block
public :: mat_unres_dense_add_block
public :: mat_unres_dense_retrieve_block
public :: mat_unres_dense_scal
public :: mat_unres_dense_scal_dia
public :: mat_unres_dense_scal_dia_vec
public :: mat_unres_dense_zero
public :: mat_unres_dense_zerohalf
public :: mat_unres_dense_write_to_disk
public :: mat_unres_dense_read_from_disk
public :: mat_unres_dense_write_to_disk2
public :: mat_unres_dense_read_from_disk2
public :: mat_unres_dense_vec_to_mat
public :: mat_unres_dense_mat_to_vec
public :: mat_unres_dense_dposv
public :: mat_unres_dense_inv
public :: mat_unres_dense_chol
public :: mat_unres_dense_sum
public :: mat_unres_dense_report_sparsity
public :: mat_unres_dense_insert_section
public :: mat_unres_dense_get_elm
public :: mat_unres_dense_create_elm_alph
public :: mat_unres_dense_create_elm_beta
public :: mat_unres_dense_part_from_full
public :: mat_unres_dense_create_ab_elms
public :: mat_unres_dense_mix_homolumo
public :: mat_unres_dense_ab_daxpy
public :: mat_unres_dense_trab_ab
public :: mat_unres_dense_to_full2
public :: mat_unres_dense_set_from_full2
public :: mat_unres_dense_add_to_fullunres
public :: mat_unres_dense_set_from_full3
public :: unres_dens_stat_deallocated_mem
public :: unres_dens_stat_allocated_mem
private

  contains
!> \brief See mat_init in mat-operations.f90
  subroutine mat_unres_dense_init(a,nrow,ncol)
     implicit none
     TYPE(Matrix) :: a
     integer, intent(in) :: nrow, ncol
     integer :: nsize
!A%elms should not be associated with anything.
!if so, the memory will be lost!!
     NULLIFY(a%elms)
     NULLIFY(a%elmsb)
     a%nrow = nrow
     a%ncol = ncol
     nsize = A%nrow*A%ncol
     ALLOCATE(a%elms(nsize),a%elmsb(nsize))
     nsize = 2*nsize*mem_realsize
     call unres_dens_stat_allocated_mem(nsize)
  end subroutine mat_unres_dense_init
  
!> \brief See mat_sum in mat-operations.f90
  function mat_unres_dense_sum(A)
    !returns sum of all elements of matrix
    implicit none
    type(Matrix), intent(in) :: A
    real(realk) :: mat_unres_dense_sum
    integer :: i,j

    mat_unres_dense_sum=0E0_realk
    do i=1,A%nrow*A%ncol
          mat_unres_dense_sum=mat_unres_dense_sum+a%elms(i)
    enddo
    do i=1,A%nrow*A%ncol
          mat_unres_dense_sum=mat_unres_dense_sum+a%elmsb(i)
    enddo

    return
  end function mat_unres_dense_sum

!> \brief See mat_free in mat-operations.f90
  subroutine mat_unres_dense_free(a)
     implicit none
     TYPE(Matrix) :: a
     integer :: nsize
     if (.not.ASSOCIATED(a%elms)) then
       print*,'memory previously released!!'
       STOP 'Error in mat_dense_free - memory previously released'
     endif
     nsize = (SIZE(a%elms)+SIZE(a%elmsb))*mem_realsize
     call unres_dens_stat_deallocated_mem(nsize)
     DEALLOCATE(a%elms,a%elmsb)
     NULLIFY(a%elms)
     NULLIFY(a%elmsb)
  end subroutine mat_unres_dense_free
  
!> \brief See mat_set_from_full in mat-operations.f90
!>
!> This routine works only for square matrices - 
!> for a(nbas,nbas), afull must be (2*nbas,2*nbas)
!>
  subroutine mat_unres_dense_set_from_full(afull,alpha,a)
     implicit none
     real(realk), INTENT(IN) :: afull(*)
     real(realk), intent(in) :: alpha
     TYPE(Matrix)            :: a 
     integer                 :: i,j,k

     if (a%nrow /= a%ncol) STOP 'non-square matrix in mat_unres_dense_set_from_full!'
     j = 1
     k = 2*a%nrow*a%nrow + a%nrow + 1
     do i = 1,a%nrow*a%nrow
        a%elms(i) = alpha*afull(j)
        a%elmsb(i) = alpha*afull(k)
        if (mod(j,a%nrow)==0) then
           j = j + a%nrow + 1
           k = k + a%nrow + 1
        else
           j = j + 1 
           k = k + 1
        endif
     enddo
  end subroutine mat_unres_dense_set_from_full

!> \brief Sets either alfa or beta part of of the type(matrix) equal to afull
!> \param afull Standard fortran matrix that should be converted (n x n)
!> \param a The output type(matrix), (n x n)
!> \param part If 'a'/'A', set alpha part equal to afull, if 'b'/'B', set beta part
  subroutine mat_unres_dense_part_from_full(afull,part,a)
     implicit none
     real(realk), INTENT(IN)    :: afull(*)
     character, intent(in)      :: part
     TYPE(Matrix),intent(inout) :: a 
     integer                    :: i

     if (part=='a' .or. part=='A') then
        do i = 1,a%nrow*a%nrow
           a%elms(i) = afull(i)
        enddo
     else if (part=='b' .or. part=='B') then
        do i = 1,a%nrow*a%nrow
           a%elmsb(i) = afull(i)
        enddo
     endif
  end subroutine mat_unres_dense_part_from_full

!> \brief See mat_to_full in mat-operations.f90
!>
!> This routine works only for square matrices -
!> for a(nbas,nbas), afull must be (2*nbas,2*nbas)
!>
  subroutine mat_unres_dense_to_full(a, alpha, afull)
     implicit none
     TYPE(Matrix), intent(in):: a
     real(realk), intent(in) :: alpha
     real(realk), intent(out):: afull(*)  
     integer                 :: i,j,k

     do i = 1,(2*a%nrow)*(2*a%nrow)
        afull(i) = 0E0_realk
     enddo
     if (a%nrow /= a%ncol) STOP 'non-square matrix in mat_unres_dense_to_full!'
     j = 1
     k = 2*a%nrow*a%nrow + a%nrow + 1
     do i = 1,a%nrow*a%nrow
        afull(j) = alpha*a%elms(i)
        afull(k) = alpha*a%elmsb(i)
        if (mod(j,a%nrow)==0) then
           j = j + a%nrow + 1
           k = k + a%nrow + 1
        else
           j = j + 1
           k = k + 1
        endif
     enddo

     end subroutine mat_unres_dense_to_full

!> \brief See mat_set_from_full2 in mat-operations.f90
  subroutine mat_unres_dense_set_from_full2(afull,alpha,a)
     implicit none
     real(realk), INTENT(IN) :: afull(*)
     real(realk), intent(in) :: alpha
     TYPE(Matrix)            :: a 
     integer                 :: i, ndim2

     ndim2=a%nrow*a%ncol
     do i = 1,ndim2
        a%elms(i)  = alpha*afull(i)
        a%elmsb(i) = alpha*afull(ndim2+i)
     enddo

  end subroutine mat_unres_dense_set_from_full2

!> \brief See mat_set_from_full2 in mat-operations.f90
  subroutine mat_unres_dense_add_to_fullunres(a,alpha,afull,spin)
     implicit none
     real(realk), INTENT(IN) :: afull(*)
     real(realk), intent(in) :: alpha
     TYPE(Matrix)            :: a 
     integer                 :: i, ndim2,spin

     ndim2=a%nrow*a%ncol
     IF(spin.EQ.1)THEN
        CALL DAXPY(ndim2,alpha,a%elms,1,Afull,1)
     ELSEIF(spin.EQ.2)THEN
        CALL DAXPY(ndim2,alpha,a%elmsb,1,Afull,1)
     ELSE
        call lsquit('unknown spin in mat_unres_dense_to_full2',-1)
     ENDIF
   end subroutine mat_unres_dense_add_to_fullunres

!> \brief See mat_to_full in mat-operations.f90
!>
!> This routine works only for square matrices -
!> for a(nbas,nbas), afull must be (2*nbas,2*nbas)
!>
  subroutine mat_unres_dense_to_full3D(a, alpha, afull,n1,n2,n3,i3,i4)
     implicit none
     integer, INTENT(IN)     :: n1,n2,n3,i3,i4
     TYPE(Matrix), intent(in):: a
     real(realk), intent(in) :: alpha
     real(realk), intent(out):: afull(n1,n2,n3)  
     integer                 :: i,j,k,n,mp1,m,offset
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

     IF(i3.NE.i4)THEN
        N = a%nrow
        M = MOD(N,7)
        IF (M.NE.0) THEN
           do j = 1,a%ncol
              offset = (j-1)*N
              DO I = 1,M
                 afull(i,j,i4) = alpha*a%elmsb(i+offset)              
              ENDDO
           enddo
           MP1 = M + 1
           IF (N.GE.7)THEN
              do j = 1,a%ncol
                 offset = (j-1)*N
                 DO I = MP1,N,7
                    afull(i,j,i4) = alpha*a%elmsb(i+offset)
                    afull(i+1,j,i4) = alpha*a%elmsb(i+1+offset)
                    afull(i+2,j,i4) = alpha*a%elmsb(i+2+offset)
                    afull(i+3,j,i4) = alpha*a%elmsb(i+3+offset)
                    afull(i+4,j,i4) = alpha*a%elmsb(i+4+offset)
                    afull(i+5,j,i4) = alpha*a%elmsb(i+5+offset)
                    afull(i+6,j,i4) = alpha*a%elmsb(i+6+offset)
                 END DO
              enddo
           ENDIF
        ELSE
           do j = 1,a%ncol
              offset = (j-1)*N
              DO I = 1,N,7
                 afull(i,j,i4) = alpha*a%elmsb(i+offset)
                 afull(i+1,j,i4) = alpha*a%elmsb(i+1+offset)
                 afull(i+2,j,i4) = alpha*a%elmsb(i+2+offset)
                 afull(i+3,j,i4) = alpha*a%elmsb(i+3+offset)
                 afull(i+4,j,i4) = alpha*a%elmsb(i+4+offset)
                 afull(i+5,j,i4) = alpha*a%elmsb(i+5+offset)
                 afull(i+6,j,i4) = alpha*a%elmsb(i+6+offset)
              END DO
           ENDDO
        ENDIF
     ELSE !i3=i4  => ADD
        N = a%nrow
        M = MOD(N,7)
        IF (M.NE.0) THEN
           do j = 1,a%ncol
              offset = (j-1)*N
              DO I = 1,M
                 afull(i,j,i4) = afull(i,j,i4) + alpha*a%elmsb(i+offset)
              ENDDO
           enddo
           MP1 = M + 1
           IF (N.GE.7)THEN
              do j = 1,a%ncol
                 offset = (j-1)*N
                 DO I = MP1,N,7
                    afull(i,j,i4) = afull(i,j,i4) + alpha*a%elmsb(i+offset)
                    afull(i+1,j,i4) = afull(i+1,j,i4) + alpha*a%elmsb(i+1+offset)
                    afull(i+2,j,i4) = afull(i+2,j,i4) + alpha*a%elmsb(i+2+offset)
                    afull(i+3,j,i4) = afull(i+3,j,i4) + alpha*a%elmsb(i+3+offset)
                    afull(i+4,j,i4) = afull(i+4,j,i4) + alpha*a%elmsb(i+4+offset)
                    afull(i+5,j,i4) = afull(i+5,j,i4) + alpha*a%elmsb(i+5+offset)
                    afull(i+6,j,i4) = afull(i+6,j,i4) + alpha*a%elmsb(i+6+offset)
                 END DO
              enddo
           ENDIF
        ELSE
           do j = 1,a%ncol
              offset = (j-1)*N
              DO I = 1,N,7
                 afull(i,j,i4) = afull(i,j,i4) + alpha*a%elmsb(i+offset)
                 afull(i+1,j,i4) = afull(i+1,j,i4) + alpha*a%elmsb(i+1+offset)
                 afull(i+2,j,i4) = afull(i+2,j,i4) + alpha*a%elmsb(i+2+offset)
                 afull(i+3,j,i4) = afull(i+3,j,i4) + alpha*a%elmsb(i+3+offset)
                 afull(i+4,j,i4) = afull(i+4,j,i4) + alpha*a%elmsb(i+4+offset)
                 afull(i+5,j,i4) = afull(i+5,j,i4) + alpha*a%elmsb(i+5+offset)
                 afull(i+6,j,i4) = afull(i+6,j,i4) + alpha*a%elmsb(i+6+offset)
              END DO
           ENDDO
        ENDIF
     ENDIF

   end subroutine mat_unres_dense_to_full3D

!> \brief See mat_set_from_full in mat-operations.f90
!>
!> For a(nbas,nbas), afull must be (nbas,nbas)
!>
  subroutine mat_unres_dense_set_from_full3(afull,alpha,a)
    implicit none
    real(realk), INTENT(IN) :: afull(*)
    real(realk), intent(in) :: alpha
    TYPE(Matrix)            :: a 
    integer                 :: i,j,k
    
    i = a%nrow*a%ncol
    call dcopy (i,afull,1,a%elms,1)
    call dcopy (i,afull,1,a%elmsb,1)
    if (alpha.ne. 1.0E0_realk) call mat_unres_dense_scal(alpha,a)
    
  end subroutine mat_unres_dense_set_from_full3
  
!> \brief See mat_to_full2 in mat-operations.f90
  subroutine mat_unres_dense_to_full2(a, alpha, afull)
     implicit none
     TYPE(Matrix), intent(in):: a
     real(realk), intent(in) :: alpha
     real(realk), intent(out):: afull(*)  
     integer                 :: i, ndim2

     ndim2=a%nrow*a%ncol
     do i = 1,ndim2
        afull(i)       = alpha*a%elms(i)
        afull(ndim2+i) = alpha*a%elmsb(i)
     enddo

  end subroutine mat_unres_dense_to_full2

!> \brief See mat_print in mat-operations.f90
  subroutine mat_unres_dense_print(a, i_row1, i_rown, j_col1, j_coln, lu)
     implicit none
     TYPE(Matrix),intent(in) :: a
     integer, intent(in)     :: i_row1, i_rown, j_col1, j_coln, lu

     !to be found in pdpack/printpkg.F
     write(lu,'(1x,a)') 'ALPHA part:'
     call LS_OUTPUT(a%elms, i_row1, i_rown, j_col1, j_coln, a%nrow, a%ncol, 1, lu)
     write(lu,'(1x,a)') 'BETA part:'
     call LS_OUTPUT(a%elmsb, i_row1, i_rown, j_col1, j_coln, a%nrow, a%ncol, &
               & 1, lu)
  end subroutine mat_unres_dense_print

!> \brief See mat_trans in mat-operations.f90
  subroutine mat_unres_dense_trans(a,b)
     implicit none
     type(matrix),intent(in)   :: a
     type(matrix)              :: b   !output  
     integer                   :: i, j

     do j = 1,a%ncol
       do i = 1,a%nrow
          b%elms(a%nrow*(i-1)+j) = a%elms(a%nrow*(j-1)+i)
          b%elmsb(a%nrow*(i-1)+j) = a%elmsb(a%nrow*(j-1)+i)
       enddo
     enddo

  end subroutine mat_unres_dense_trans

!> \brief See mat_assign in mat-operations.f90
  subroutine mat_unres_dense_assign(a,b)
     implicit none
     TYPE(Matrix), INTENT(INOUT) :: a
     TYPE(Matrix), INTENT(IN)    :: b
     integer                     :: i

     !do i = 1,a%nrow*a%ncol
     !   a%elms(i) = b%elms(i)
     !   a%elmsb(i) = b%elmsb(i)
     !enddo
     i = a%nrow*a%ncol
     call dcopy (i,b%elms,1,a%elms,1)
     call dcopy (i,b%elmsb,1,a%elmsb,1)
  end subroutine mat_unres_dense_assign

#ifndef UNITTEST
!> \brief See mat_mpicopy in mat-operations.f90
  subroutine mat_unres_dense_mpicopy(a,slave,master)
    use lsmpi_type
     implicit none
     TYPE(Matrix), TARGET, INTENT(INOUT) :: a
     logical                     :: slave
     integer(kind=ls_mpik)       :: master
     integer                     :: i

     call LS_MPI_BUFFER(a%nrow,Master)
     call LS_MPI_BUFFER(a%ncol,Master)
     IF(Slave)THEN
        NULLIFY(a%iaux, a%raux)
        nullify(A%elms)
        nullify(A%elmsb)
        call mat_unres_dense_init(a,a%nrow,a%ncol)
        a%init_self_ptr => a
        a%init_magic_tag = mat_init_magic_value
     ENDIF
     i = size(a%elms)!a%nrow*a%ncol
     call LS_MPI_BUFFER(a%elms,i,Master)
     i = size(a%elmsb)!a%nrow*a%ncol
     call LS_MPI_BUFFER(a%elmsb,i,Master)

   end subroutine mat_unres_dense_mpicopy
#endif

!> \brief See mat_copy in mat-operations.f90
  subroutine mat_unres_dense_copy(alpha,a,b)
     implicit none
     REAL(REALK),  INTENT(IN)    :: alpha
     TYPE(Matrix), INTENT(IN)    :: a
     TYPE(Matrix), INTENT(INOUT) :: b
     integer                     :: i

     !do i = 1,a%nrow*a%ncol
     !   b%elms(i) = alpha*a%elms(i)
     !   b%elmsb(i) = alpha*a%elmsb(i)
     !enddo
     call mat_unres_dense_assign(b,a)
     if (alpha.ne. 1.0E0_realk) call mat_unres_dense_scal(alpha,b)
  end subroutine mat_unres_dense_copy

!> \brief See mat_tr in mat-operations.f90
  function mat_unres_dense_Tr(a)
     implicit none
     TYPE(Matrix), intent(IN) :: a
     REAL(realk) :: mat_unres_dense_tr 
     integer :: i

     mat_unres_dense_tr = 0E0_realk
     do i = 1,a%nrow
        mat_unres_dense_tr = mat_unres_dense_tr + a%elms((a%Nrow+1)*i-a%Nrow) + &
                                              &   a%elmsb((a%Nrow+1)*i-a%Nrow)
    enddo

  end function mat_unres_dense_Tr
 
!> \brief See mat_TrAB in mat-operations.f90
  function mat_unres_dense_TrAB(a,b)
     implicit none
     TYPE(Matrix), intent(IN) :: a,b
     REAL(realk) :: mat_unres_dense_trAB 
     real(realk), external :: ddot
     integer :: i, j

     mat_unres_dense_TrAB = 0.0E0_realk
     !do i = 1,a%nrow
       do j = 1,a%ncol
         !mat_unres_dense_TrAB = mat_unres_dense_TrAB + &
         !                               &  a%elms(a%nrow*(j-1)+i)*  &
         !                               &  b%elms(b%nrow*(i-1)+j) + &
         !                               &  a%elmsb(a%nrow*(j-1)+i)*  &
         !                               &  b%elmsb(b%nrow*(i-1)+j) 
          mat_unres_dense_TrAB = mat_unres_dense_TrAB + &
           & ddot(a%nrow,a%elms(a%nrow*(j-1)+1),1,b%elms(j),a%ncol) + &
           & ddot(a%nrow,a%elmsb(a%nrow*(j-1)+1),1,b%elmsb(j),a%ncol)
       enddo
     !enddo

  end function mat_unres_dense_TrAB 

!> \brief See mat_TrAB_ab in mat-operations.f90
subroutine mat_unres_dense_TrAB_ab (A, B, trace)

implicit none

type(matrix), intent(in) :: A, B
real(realk), intent(out) :: trace(:)

integer                  :: r, c, i
real(realk), external    :: ddot

r = A%nrow ; c = A%ncol
trace(1) = 0.0E0_realk ; trace(2) = 0.0E0_realk
do i=1,c
  trace(1) = trace(1) + ddot(r, A%elms(r*(i-1)+1), 1, B%elms(i), c)
  trace(2) = trace(2) + ddot(r, A%elmsb(r*(i-1)+1), 1, B%elmsb(i), c)
enddo

end subroutine mat_unres_dense_TrAB_ab

!> \brief See mat_mul in mat-operations.f90
  subroutine mat_unres_dense_mul(a,b,transa,transb,alpha,beta,c) 
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
       STOP 'unknown format in mat_dense_mul'
     endif
     if (transb == 'n' .or. transb == 'N') then
       n = b%ncol
       if (b%nrow /= k) STOP 'weird'
     elseif (transb == 't' .or. transb == 'T') then
       n = b%nrow
       if (b%ncol /= k) STOP 'weird'
     endif
     call DGEMM(transa,transb,m,n,k,alpha,&
               & a%elms,a%nrow,b%elms,b%nrow,beta,c%elms,c%nrow)
     call DGEMM(transa,transb,m,n,k,alpha,&
               & a%elmsb,a%nrow,b%elmsb,b%nrow,beta,c%elmsb,c%nrow)

  end subroutine mat_unres_dense_mul

!> \brief See mat_add in mat-operations.f90
  subroutine mat_unres_dense_add(alpha,a,beta,b,c)
     implicit none
     TYPE(Matrix), intent(IN) :: a, b
     REAL(realk), INTENT(IN)  :: alpha, beta
     TYPE(Matrix)             :: c
     integer                  :: i

     !do i = 1,a%nrow*a%ncol
     !   c%elms(i) = alpha*a%elms(i) + beta*b%elms(i)
     !   c%elmsb(i) = alpha*a%elmsb(i) + beta*b%elmsb(i)
     !enddo

     call mat_unres_dense_copy(alpha,a,c)
     call mat_unres_dense_daxpy(beta,b,c)
     
  end subroutine mat_unres_dense_add

!> \brief See mat_daxpy in mat-operations.f90
  subroutine mat_unres_dense_daxpy(alpha,x,y)
     implicit none
     real(realk),intent(in)       :: alpha
     TYPE(Matrix), intent(IN)     :: X
     TYPE(Matrix), intent(INOUT)  :: Y
     integer                      :: i

     !do i = 1,x%nrow*x%ncol
     !   y%elms(i)  = y%elms(i) + alpha*x%elms(i)
     !   y%elmsb(i) = y%elmsb(i) + alpha*x%elmsb(i)
     !enddo
     i = x%nrow*x%ncol
     call daxpy(i,alpha,x%elms,1,y%elms,1)
     call daxpy(i,alpha,x%elmsb,1,y%elmsb,1)

  end subroutine mat_unres_dense_daxpy

  !> \brief creates the inverse matrix of type(matrix). mat_inv
  !> \author T. Kjærgaard
  !> \date 2012
  !> \param a The type(matrix) that should be inversed
  !> \param chol The type(matrix) that contains cholesky factors (from mat_chol)
  !> \param c The inverse output type(matrix).
  SUBROUTINE mat_unres_dense_inv(A, A_inv) 
    implicit none
    TYPE(Matrix),intent(in)     :: A
    TYPE(Matrix)                :: A_inv !output
    real(realk), pointer   :: work1(:)
    real(realk), pointer   :: A_inv_full(:,:) 
    integer,pointer    :: IPVT(:)
    real(realk)            :: RCOND, dummy(2), tmstart, tmend
    integer                :: IERR, i, j, fulldim, ndim
    fulldim = 2*a%nrow
    call mem_alloc(A_inv_full,fulldim,fulldim) 
    call mem_alloc(work1,fulldim)
    call mem_alloc(IPVT,fulldim)
    !Invert U and Ut:
    IPVT = 0 ; RCOND = 0.0E0_realk  
    call mat_unres_dense_to_full(A,1.0E0_realk,A_inv_full)
    call DGECO(A_inv_full,fulldim,fulldim,IPVT,RCOND,work1)
    call DGEDI(A_inv_full,fulldim,fulldim,IPVT,dummy,work1,01)
    !Convert framework:
    call mat_unres_dense_set_from_full(A_inv_full,1.0E0_realk,A_inv)
    call mem_dealloc(A_inv_full) 
    call mem_dealloc(work1)
    call mem_dealloc(IPVT)
  END SUBROUTINE mat_unres_dense_inv
  

!> \brief See mat_ab_daxpy in mat-operations.f90
!=======================================================================
subroutine mat_unres_dense_ab_daxpy (alpha, X, Y)

implicit none

real(realk), intent(in)     :: alpha(2)
type(matrix), intent(in)    :: X
type(matrix), intent(inout) :: Y
integer                     :: i

i = X%nrow * X%ncol

call daxpy (i, alpha(1), X%elms, 1, Y%elms, 1)
call daxpy (i, alpha(2), X%elmsb, 1, Y%elmsb, 1)

end subroutine mat_unres_dense_ab_daxpy
!=======================================================================

!> \brief See mat_dotproduct in mat-operations.f90
  function mat_unres_dense_dotproduct(a,b)
     implicit none
     TYPE(Matrix), intent(IN) :: a,b
     REAL(realk) :: mat_unres_dense_dotproduct
     real(realk), external :: ddot
     integer     :: i

     mat_unres_dense_dotproduct = 0.0E0_realk
     !do i = 1,a%nrow*a%ncol
     !   mat_unres_dense_dotproduct = mat_unres_dense_dotproduct + &
     !        & a%elms(i)*b%elms(i) + a%elmsb(i)*b%elmsb(i)
     !enddo
     i = a%nrow*a%ncol
     mat_unres_dense_dotproduct = ddot(i,a%elms,1,b%elms,1) + ddot(i,a%elmsb,1,b%elmsb,1)

  end function mat_unres_dense_dotproduct

!> \brief See mat_sqnorm2 in mat-operations.f90
  function mat_unres_dense_sqnorm2(a)
     implicit none
     TYPE(Matrix), intent(IN) :: a
     REAL(realk) :: mat_unres_dense_sqnorm2
     integer     :: i

     mat_unres_dense_sqnorm2 = 0.0E0_realk
     !do i = 1,a%nrow*a%ncol
     !   mat_unres_dense_sqnorm2 = mat_unres_dense_sqnorm2 + &
     !        & a%elms(i)*a%elms(i) +  a%elmsb(i)*a%elmsb(i)
     !enddo
     mat_unres_dense_sqnorm2 = mat_unres_dense_dotproduct(a,a)
  end function mat_unres_dense_sqnorm2

!> \brief See mat_outdia_sqnorm2 in mat-operations.f90
  function mat_unres_dense_outdia_sqnorm2(a)
     implicit none
     TYPE(Matrix), intent(IN) :: a
     REAL(realk) :: mat_unres_dense_outdia_sqnorm2
     integer     :: i,n

     mat_unres_dense_outdia_sqnorm2 = 0.0E0_realk
     n = 1
     do i = 1,a%nrow*a%ncol
       if (i /= (n-1)*a%nrow+n) then
         mat_unres_dense_outdia_sqnorm2 = mat_unres_dense_outdia_sqnorm2 + &
                      & a%elms(i)*a%elms(i) + a%elmsb(i)*a%elmsb(i)
       else
         n = n+1
       endif
     enddo

  end function mat_unres_dense_outdia_sqnorm2

!> \brief See mat_dsyev in mat-operations.f90
  SUBROUTINE mat_unres_dense_dsyev(S,eival,ndim)
    implicit none
    TYPE(Matrix), intent(INOUT) :: S
    real(realk), intent(INOUT) :: eival(2*ndim)
    integer,intent(in) :: ndim
!
    real(realk),pointer :: work(:)
    real(realk), allocatable :: eivala(:), eivalb(:)
    integer :: infdiag,lwork
    infdiag=0

    allocate (eivala(ndim), eivalb(ndim))

    lwork = -1
    call mem_alloc(work,5)
    ! we inquire the size of lwork
    call DSYEV('V','U',ndim,S%elms,ndim,eivala,work,lwork,infdiag)
    lwork = work(1)
    call mem_dealloc(work)
    call mem_alloc(work,lwork)
    call DSYEV('V','U',ndim,S%elms,ndim,eivala,work,lwork,infdiag)
    call mem_dealloc(work)
    if(infdiag.ne. 0) then
       print*,'mat_dsyev: dsyev failed, info=',infdiag
       call lsquit('mat_dsyev: diagonalization failed.',-1)
    end if

    lwork = -1
    call mem_alloc(work,5)
    ! we inquire the size of lwork
    call DSYEV('V','U',ndim,S%elmsb,ndim,eivalb,work,lwork,infdiag)
    lwork = work(1)
    call mem_dealloc(work)
    call mem_alloc(work,lwork)
    call DSYEV('V','U',ndim,S%elmsb,ndim,eivalb,work,lwork,infdiag)
    call mem_dealloc(work)
    if(infdiag.ne. 0) then
       print*,'mat_dsyev: dsyev failed, info=',infdiag
       call lsquit('mat_dsyev: diagonalization failed.',-1)
    end if

    eival(1:ndim) = eivala
    eival(ndim+1:2*ndim) = eivalb

    deallocate (eivala, eivalb)

  END SUBROUTINE mat_unres_dense_dsyev

!> \brief See mat_diag_f in mat-operations.f90
  subroutine mat_unres_dense_diag_f(F,S,eival,Cmo)
    !solves FC = SCe 
    implicit none
    TYPE(Matrix), intent(IN) :: F,S
    real(realk), intent(OUT) :: eival(:) ! eigenvalues
    type(matrix)             :: Cmo  !output
    real(realk), allocatable :: wrk(:),tmp(:)
    integer :: ndim,i
   
    ndim = S%nrow
    ALLOCATE(tmp(Ndim*Ndim))
    do i = 1,ndim*ndim
       tmp(i) = S%elms(i)
    enddo
    call  mat_unres_dense_assign(Cmo, F)
    !diagonalise alpha part
    call my_DSYGV(ndim,Cmo%elms,tmp,eival(1:ndim),&
                & "DENSE_GET_DENS_SOL,a")

    do i = 1,ndim*ndim
       tmp(i) = S%elmsb(i)
    enddo
    !diagonalise beta part
    call my_DSYGV(ndim,Cmo%elmsb,tmp,eival(ndim+1:2*ndim),&
                & "DENSE_GET_DENS_SOL,b")


    DEALLOCATE(tmp)
  end subroutine mat_unres_dense_diag_f

!> \brief See mat_max_elm in mat-operations.f90
  subroutine mat_unres_dense_max_elm(a, val)
  implicit none
     type(matrix),intent(in)  :: a
     real(realk), intent(out) :: val
     integer                  :: i

     val = a%elms(1)
     do i = 1, a%nrow*a%ncol 
        if (a%elms(i) > val) then
           val = a%elms(i)
        endif
        if (a%elmsb(i) > val) then
           val = a%elmsb(i)
        endif
     enddo
  end subroutine mat_unres_dense_max_elm

!> \brief See mat_abs_max_elm in mat-operations.f90
  subroutine mat_unres_dense_abs_max_elm(a, val)
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
     val = MAX(val,abs(a%elmsb(1)))
     do i = 1, a%nrow*a%ncol 
        if (abs(a%elmsb(i)) > val) then
           val = abs(a%elmsb(i))
        endif
     enddo
   end subroutine mat_unres_dense_abs_max_elm

!> \brief See mat_max_diag_elm in mat-operations.f90
  subroutine mat_unres_dense_max_diag_elm(a,pos,val)
  implicit none
     type(matrix),intent(in)  :: a
     real(realk), intent(out) :: val
     integer                  :: i
     integer, intent(out)     :: pos

     val = abs(a%elms(1))
     pos = 1
     do i = 2, a%nrow
        if (abs(a%elms(a%nrow*(i-1)+i)) > abs(val)) then
           val = a%elms(a%nrow*(i-1)+i)
           pos = i
        endif
        if (abs(a%elmsb(a%nrow*(i-1)+i)) > abs(val)) then
           val = a%elmsb(a%nrow*(i-1)+i)
           pos = i
        endif
     enddo

  end subroutine mat_unres_dense_max_diag_elm

!> \brief See mat_column_norm in mat-operations.f90
  FUNCTION mat_unres_dense_column_norm(Mat,ncol,from_row,to_row)
    !Returns the sum of the elements A(from_row:to_row,col_num) squared
     implicit none
     type(Matrix), intent(in) :: Mat
     integer, intent(in) :: ncol, from_row, to_row
     real(realk) :: mat_unres_dense_column_norm
     integer :: i

     mat_unres_dense_column_norm = 0E0_realk
     do i = from_row,to_row
        mat_unres_dense_column_norm = mat_unres_dense_column_norm + &
             & (Mat%elms((ncol-1)*mat%nrow+i))**2     + &
             & (Mat%elmsb((ncol-1)*mat%nrow+i))**2    
     enddo  
  end function mat_unres_dense_column_norm

!> \brief See mat_section in mat-operations.f90
  subroutine mat_unres_dense_section(A,from_row,to_row,from_col,to_col,Asec)
     implicit none
     type(Matrix), intent(in) :: A
     integer, intent(in) :: from_row, to_row, from_col, to_col
     type(Matrix), intent(inout) :: Asec  !output
     integer :: i, index, irow, icol

!     Asec%nrow = to_row - from_row + 1
!     Asec%ncol = to_col - from_col + 1
!     if (SIZE(Asec%elms) < Asec%nrow*Asec%ncol*2) &
!          & STOP 'too little memory allocated in mat_dense_section'
     index = (from_col-1)*A%nrow + from_row
     irow = 0
     icol = from_col
     do i = 1,Asec%nrow*Asec%ncol
       Asec%elms(i) = A%elms(index)
       Asec%elmsb(i) = A%elmsb(index)
       irow = irow + 1
       if (irow == Asec%nrow) then
         icol = icol + 1
         index = (icol - 1)*A%nrow + from_row
         irow = 0
       else
         index = index + 1
       endif
     enddo

  end subroutine mat_unres_dense_section

!> \brief See mat_insert_section in mat-operations-essentials.f90
subroutine mat_unres_dense_insert_section (Asec, from_row, to_row, from_col, to_col, A)

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
  A%elmsb(index) = Asec%elmsb(i)
  irow = irow + 1
  if (irow == Asec%nrow) then
    icol = icol + 1
    index = (icol - 1)*A%nrow + from_row
    irow = 0
  else
    index = index + 1
  endif
enddo

end subroutine mat_unres_dense_insert_section
 
!> \brief See mat_mix_homolumo in mat_operations.f90
subroutine mat_unres_dense_mix_homolumo (plus, minus, homo, lumo)

implicit none

type(matrix), intent(in)    :: plus, minus
type(matrix), intent(inout) :: homo, lumo

homo%elms = plus%elms
homo%elmsb = minus%elmsb

lumo%elms = minus%elms
lumo%elmsb = plus%elmsb

end subroutine mat_unres_dense_mix_homolumo
 
!> \brief See mat_precond in mat-operations.f90
  subroutine mat_unres_dense_precond(M,x,xprec)
    implicit none
    type(Matrix), intent(inout) :: xprec
    type(Matrix), intent(in)    :: M, x
    integer :: i

    do i = 1, x%nrow
       if (abs(M%elms((i-1)*M%nrow+i)) > 1.0E-9_realk) then
          xprec%elms(i) = x%elms(i)/M%elms((i-1)*M%nrow+i)
       endif
       if (abs(M%elmsb((i-1)*M%nrow+i)) > 1.0E-9_realk) then
          xprec%elmsb(i) = x%elmsb(i)/M%elmsb((i-1)*M%nrow+i)
       endif
    enddo

  end subroutine mat_unres_dense_precond

!> \brief See mat_mo_precond in mat-operations.f90
  subroutine mat_unres_dense_mo_precond(nocc,omega,Eorb_final,X_MO)
   !which orbital energies should be used for precond?? and how?
    implicit none
    integer, intent(in) :: nocc
    real(realk), intent(in) :: omega, Eorb_final(:)
    type(Matrix), intent(inout) :: X_MO
    real(realk) :: dia
    integer :: i,j,nvirt

    !modify lower left block (ai block)
    do j = 1,nocc  !columns
      do i = nocc+1,X_MO%nrow  !rows
        !dia = E[2]_dia + omega*S[2]_dia
        !NOTE TO SELF: WHY the 2.d0s and not 1E0_realk????
        dia = 2.0E0_realk*(Eorb_final(i) - Eorb_final(j)) + omega * 2E0_realk
        X_MO%elms((j-1)*X_MO%nrow+i) = X_MO%elms((j-1)*X_MO%nrow+i) /dia
        X_MO%elmsb((j-1)*X_MO%nrow+i) = X_MO%elmsb((j-1)*X_MO%nrow+i) /dia
      enddo
    enddo
    !modify upper right block (ia block) 
    do j = nocc+1,X_MO%ncol  !columns
      do i = 1,nocc          !rows
        !dia = E[2]_dia + omega*S[2]_dia
        dia = 2.0E0_realk*(Eorb_final(j) - Eorb_final(i)) - omega * 2E0_realk
        X_MO%elms((j-1)*X_MO%nrow+i) = X_MO%elms((j-1)*X_MO%nrow+i)/dia
        X_MO%elmsb((j-1)*X_MO%nrow+i) = X_MO%elmsb((j-1)*X_MO%nrow+i)/dia
      enddo
    enddo

  end subroutine mat_unres_dense_mo_precond

!> \brief See mat_ao_precond in mat-operations.f90
  subroutine mat_unres_dense_ao_precond(symmetry,omega,FUP,FUQ,DU,X_AO)
    implicit none
    integer, intent(in) :: symmetry
    real(realk), intent(in) :: omega
    type(Matrix), intent(in) :: FUP, FUQ, DU
    type(Matrix), intent(inout) :: X_AO
    real(realk) :: denom,err
    integer :: ndim,i,j

    ndim = FUP%nrow
    if (symmetry == 1 .or. symmetry == 2) then !Symmetric or antisymmetric X_AO
       do j = 1,ndim   !columns
         do i = 1,ndim  !rows
           !alpha part
           denom = FUQ%elms((j-1)*ndim+j) + FUQ%elms((i-1)*ndim+i)  &
                &- FUP%elms((i-1)*ndim+i) - FUP%elms((j-1)*ndim+j)  &  
                &- omega
           if (ABS(denom) > 1.0E-10_realk) then
              X_AO%elms((j-1)*ndim+i) = X_AO%elms((j-1)*ndim+i)/(denom)
           endif
           !beta part
           denom = FUQ%elmsb((j-1)*ndim+j) + FUQ%elmsb((i-1)*ndim+i)  &
                &- FUP%elmsb((i-1)*ndim+i) - FUP%elmsb((j-1)*ndim+j)  &  
                &- omega
           if (ABS(denom) > 1.0E-10_realk) then
              X_AO%elmsb((j-1)*ndim+i) = X_AO%elmsb((j-1)*ndim+i)/(denom)
           endif
         enddo
       enddo
    else  !X_AO not symmetric in any way
       do j = 1,ndim   !columns
         do i = 1,ndim  !rows
           !alpha part
           denom = FUQ%elms((j-1)*ndim+j) + FUQ%elms((i-1)*ndim+i)  &
                &- FUP%elms((i-1)*ndim+i) - FUP%elms((j-1)*ndim+j)  &
                &- omega*(DU%elms((j-1)*ndim+j) - DU%elms((i-1)*ndim+i))
           if (abs(denom) < 1.0E-1_realk) then !Do not divide by too small elements
              denom = denom*1.0E-1_realk/(abs(denom)) !Keep the sign on denominator
           endif
           X_AO%elms((j-1)*ndim+i) = X_AO%elms((j-1)*ndim+i)/(denom)
           !beta part
           denom = FUQ%elmsb((j-1)*ndim+j) + FUQ%elmsb((i-1)*ndim+i)  &
                &- FUP%elmsb((i-1)*ndim+i) - FUP%elmsb((j-1)*ndim+j)  &
                &- omega*(DU%elmsb((j-1)*ndim+j) - DU%elmsb((i-1)*ndim+i))
           if (abs(denom) < 1.0E-1_realk) then !Do not divide by too small elements
              denom = denom*1.0E-1_realk/(abs(denom)) !Keep the sign on denominator
           endif
           X_AO%elmsb((j-1)*ndim+i) = X_AO%elmsb((j-1)*ndim+i)/(denom)
         enddo
       enddo
    endif
  end subroutine mat_unres_dense_ao_precond

!> \brief See mat_identity in mat-operations.f90
  subroutine mat_unres_dense_identity(Id)
    implicit none
    type(matrix), intent(inout) :: Id
    integer :: i,j

    do j = 1,Id%ncol
      do i = 1,Id%nrow
        if (i == j) then
          Id%elms((j-1)*Id%nrow+i) = 1.0E0_realk
          Id%elmsb((j-1)*Id%nrow+i) = 1.0E0_realk
        else
          Id%elms((j-1)*Id%nrow+i) = 0.0E0_realk
          Id%elmsb((j-1)*Id%nrow+i) = 0.0E0_realk 
        endif
      enddo
    enddo
  end subroutine mat_unres_dense_identity

!> \brief See mat_create_elm in mat-operations.f90
  subroutine mat_unres_dense_create_elm(i,j,val,A)
    implicit none
    integer, intent(in) :: i,j
    real(realk), intent(in) :: val
    type(Matrix), intent(inout) :: A

    A%elms((j-1)*A%nrow+i) = val
    A%elmsb((j-1)*A%nrow+i) = val

  end subroutine mat_unres_dense_create_elm

!> \brief Create only element in alpha part.
  subroutine mat_unres_dense_create_elm_alph(i,j,val,A)
    implicit none
    integer, intent(in) :: i,j
    real(realk), intent(in) :: val
    type(Matrix), intent(inout) :: A

    A%elms((j-1)*A%nrow+i) = val

  end subroutine mat_unres_dense_create_elm_alph

!> \brief Create only element in beta part.
  subroutine mat_unres_dense_create_elm_beta(i,j,val,A)
    implicit none
    integer, intent(in) :: i,j
    real(realk), intent(in) :: val
    type(Matrix), intent(inout) :: A

    A%elmsb((j-1)*A%nrow+i) = val

  end subroutine mat_unres_dense_create_elm_beta


!> \brief See mat_create_ab_elm in mat-operations.f90
!=======================================================================
subroutine mat_unres_dense_create_ab_elms (r, c, elm, A)

implicit none

integer, intent(in)         :: r, c
real(realk), intent(in)     :: elm(:)
type(matrix), intent(inout) :: A

A%elms((c-1)*A%nrow+r)  = elm(1)
A%elmsb((c-1)*A%nrow+r) = elm(2)

end subroutine mat_unres_dense_create_ab_elms
!=======================================================================

!> \brief See mat_get_ab_elms in mat-operations.f90
!=======================================================================
  subroutine mat_unres_dense_get_elm (A, r, c, elm)

  implicit none

  type(matrix), intent(in) :: A
  integer, intent(in)      :: r, c
  real(realk), intent(out) :: elm(:)

  if (size(elm) /= 2)THEN
     stop 'wrong dimension of elm in mat_unres_dense_get_elm'
  ENDIF
  elm(1) = A%elms(A%nrow*(c-1)+r)
  elm(2) = A%elmsb(A%nrow*(c-1)+r)

  end subroutine mat_unres_dense_get_elm
!=======================================================================

!> \brief See mat_create_block in mat-operations.f90
  subroutine mat_unres_dense_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
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
    do j = insertcol, insertcol+fullcol-1
       do i = insertrow, insertrow+fullrow-1
          A%elmsb((j-1)*A%nrow+i) = fullmat(i-insertrow+1,j-insertcol+1)
       enddo
    enddo
  end subroutine mat_unres_dense_create_block

!> \brief See mat_add_block in mat-operations.f90
  subroutine mat_unres_dense_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in)         :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(in)     :: fullmat(fullrow,fullcol)
    type(Matrix), intent(inout) :: A
    integer                     :: i, j
    
    do j = insertcol, insertcol+fullcol-1
       do i = insertrow, insertrow+fullrow-1
          A%elms((j-1)*A%nrow+i) = A%elms((j-1)*A%nrow+i)+fullmat(i-insertrow+1,j-insertcol+1)
       enddo
    enddo
    do j = insertcol, insertcol+fullcol-1
       do i = insertrow, insertrow+fullrow-1
          A%elmsb((j-1)*A%nrow+i) = A%elmsb((j-1)*A%nrow+i)+fullmat(i-insertrow+1,j-insertcol+1)
       enddo
    enddo
  end subroutine mat_unres_dense_add_block
  
!> \brief See mat_retrieve_block in mat-operations.f90
  subroutine mat_unres_dense_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in)         :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(out)    :: fullmat(fullrow,fullcol)
    type(Matrix), intent(in) :: A
    integer                     :: i, j
    do j = insertcol, insertcol+fullcol-1
       do i = insertrow, insertrow+fullrow-1
           fullmat(i-insertrow+1,j-insertcol+1) = 0.5E0_realk*(A%elms((j-1)*A%nrow+i)+A%elmsb((j-1)*A%nrow+i))
       enddo
    enddo
  end subroutine mat_unres_dense_retrieve_block

!> \brief See mat_scal in mat-operations.f90
  subroutine mat_unres_dense_scal(alpha,A)
    implicit none
    real(realk), intent(in) :: alpha
    type(Matrix), intent(inout) :: A
    integer :: i

    !do i = 1,A%nrow*A%ncol
    !   A%elms(i) = A%elms(i) * alpha
    !   A%elmsb(i) = A%elmsb(i) * alpha
    !enddo
    i = A%nrow*A%ncol
    call dscal(i,alpha,A%elms,1)
    call dscal(i,alpha,A%elmsb,1)
  end subroutine mat_unres_dense_scal

!> \brief See mat_scal_dia in mat-operations.f90
  subroutine mat_unres_dense_scal_dia(alpha,A)
    implicit none
    real(realk), intent(in) :: alpha
    type(Matrix), intent(inout) :: A
    integer :: i

    !do i = 1,A%nrow
    !   A%elms((i-1)*A%nrow+i) = A%elms((i-1)*A%nrow+i) * alpha
    !   A%elmsb((i-1)*A%nrow+i) = A%elmsb((i-1)*A%nrow+i) * alpha
    !enddo    
    call dscal(A%nrow,alpha,A%elms,A%nrow+1)
    call dscal(A%nrow,alpha,A%elmsb,A%nrow+1)
  end subroutine mat_unres_dense_scal_dia

!> \brief See mat_zero in mat-operations.f90
  subroutine mat_unres_dense_zero(A)
    implicit none
    type(Matrix), intent(inout) :: A
    integer :: i

    do i = 1,A%nrow*A%ncol
       A%elms(i) = 0.0E0_realk
       A%elmsb(i) = 0.0E0_realk
    enddo
 
  end subroutine mat_unres_dense_zero

!> \brief See mat_zerohalf in mat-operations.f90
  subroutine mat_unres_dense_zerohalf(part,A)
    implicit none
    character(len=2), intent(in) :: part
    type(Matrix), intent(inout) :: A
    integer :: i,j

    if (part == 'ut' .or. part == 'UT') then
      !set the upper triangle to zero - diagonal is kept
      do j = 2,A%ncol
        do i = 1,j-1
           A%elms((j-1)*A%nrow+i) = 0.0E0_realk
           A%elmsb((j-1)*A%nrow+i) = 0.0E0_realk
        enddo
      enddo
    else
      !set the lower triangle to zero - diagonal is also zeroed
      do j = 1,A%ncol
        do i = j,A%nrow
           A%elms((j-1)*A%nrow+i) = 0.0E0_realk
           A%elmsb((j-1)*A%nrow+i) = 0.0E0_realk
        enddo
      enddo
    endif

  end subroutine mat_unres_dense_zerohalf

!> \brief See mat_write_to_disk in mat-operations.f90
  subroutine mat_unres_dense_write_to_disk(iunit,A)
  !===========================================================================
  ! Writes matrix A to file
  !           INPUT: iunit (integer) - The unit number
  !                  A (matrix) - The matrix that should be written to disk
    implicit none
    integer, intent(in) :: iunit
    integer :: i
    type(Matrix), intent(in) :: A
    character(len=20) :: filename
    integer(kind=long) :: ncol,nrow 
    if (.not.ASSOCIATED(A%elms)) &
         & STOP 'A in mat_unres_dense_WRITE_TO_DISK non-existant'

 !   INQUIRE(iunit,NAME=filename)     !Find the filename of the file
 !   print*,filename
 !   CLOSE(iunit)
!TODO: Check this
!    OPEN(iunit, FILE=filename, STATUS="old", FORM="unformatted", POSITION="append" )
!    OPEN(iunit, FILE=filename, STATUS="old", FORM="unformatted")
 !   REWIND(iunit)
    nrow = A%Nrow
    ncol = A%Ncol
    write(iunit) Nrow, Ncol
    WRITE(iunit)(A%elms(I),I=1,A%nrow*A%ncol)
    WRITE(iunit)(A%elmsb(I),I=1,A%nrow*A%ncol)
 !   WRITE(iunit) A%elms
 !   WRITE(iunit) A%elmsb
 !   ENDFILE iunit
  end subroutine mat_unres_dense_WRITE_TO_DISK

!> \brief See mat_read_from_disk in mat-operations.f90
  subroutine mat_unres_dense_read_from_disk(iunit,A)
  !===============================================================================
  ! Reads a matrix from disk
  !            INPUT: iunit (integer) - The unit number
  !           OUTPUT: A (matrix) - The matrix wanted.
    implicit none
    integer, intent(in) :: iunit
    type(Matrix), intent(inout) :: A
    integer :: i
    integer(kind=long) :: ncol,nrow 

!    REWIND iunit
    READ(iunit) Nrow, Ncol
    IF(nrow .NE. A%Nrow) call lsquit( 'mat_unres_dense_read_from_disk: Nrow /= A%nrow',-1)
    IF(ncol .NE. A%Ncol) call lsquit( 'mat_unres_dense_read_from_disk: Ncol /= A%ncol',-1)
    READ(iunit)(A%elms(I),I=1,A%nrow*A%ncol)
    READ(iunit)(A%elmsb(I),I=1,A%nrow*A%ncol)
    !READ(iunit) A%elms
    !READ(iunit) A%elmsb

  end subroutine mat_unres_dense_READ_FROM_DISK

!> \brief See mat_vec_to_mat in mat-operations.f90
    subroutine mat_unres_dense_vec_to_mat(symmetry, VEC, MAT)
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
             MAT%elmsb((m-1)*MAT%nrow + n) = VEC%elmsb(i)
             MAT%elmsb((n-1)*MAT%nrow + m) = VEC%elmsb(i)
          enddo
       enddo
    else if (symmetry == 'a' .OR. symmetry == 'A') then
       k = 0
       do n = 1, MAT%nrow
          do m = 1, n
             if (n == m) then
                MAT%elms((n-1)*MAT%nrow + m) = 0.0E0_realk !Diagonal element are 0 (antisym. matrix)
                MAT%elmsb((n-1)*MAT%nrow + m) = 0.0E0_realk !Diagonal element are 0 (antisym. matrix)
             else
                k = k + 1
                MAT%elms((m-1)*MAT%nrow + n) = VEC%elms(k)
                MAT%elms((n-1)*MAT%nrow + m) = -VEC%elms(k)
                MAT%elmsb((m-1)*MAT%nrow + n) = VEC%elmsb(k)
                MAT%elmsb((n-1)*MAT%nrow + m) = -VEC%elmsb(k)
             endif
          enddo
       enddo
    endif
    end subroutine mat_unres_dense_vec_to_mat
 
!> \brief See mat_to_vec in mat-operations.f90
    subroutine mat_unres_dense_mat_to_vec(symmetry, MAT, VEC)
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
             VEC%elmsb(i) = MAT%elmsb((m-1)*MAT%nrow + n)
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
                VEC%elmsb(k) = MAT%elmsb((m-1)*MAT%nrow + n)
             endif
          enddo
       enddo
    endif
    
    end subroutine mat_unres_dense_mat_to_vec

!> \brief See stat_allocated_memory in mat-operations.f90
  subroutine unres_dens_stat_allocated_mem(size)
     implicit none
     integer, intent(in) :: size
     integer(kind=long) :: nsize
     nsize = size
    
     call mem_allocated_mem_type_matrix(nsize)
  
  end subroutine unres_dens_stat_allocated_mem

!> \brief See stat_deallocated_memory in mat-operations.f90
  subroutine unres_dens_stat_deallocated_mem(size)
     implicit none
     integer, intent(in) :: size
     integer(kind=long) :: nsize
     nsize = size
     call mem_deallocated_mem_type_matrix(nsize)

  end subroutine unres_dens_stat_deallocated_mem

end module matrix_op_unres_dense

