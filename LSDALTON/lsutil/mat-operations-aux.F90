!> @file 
!> Contains additional non essentail matrix operations module

!> Contains wrapper routines that branch out to matrix routine for chosen matrix type.
!> \author L. Thogersen
!> \date 2003
!>
!> General rules:
!> NEVER put a type(Matrix) as intent(out), this will on some platforms
!>       make the pointer
!>       disassociated entering the routine, and memory already
!>       allocated for the matrix will be lost. \n
!> NEVER think that e.g. A%elms = matrix will copy the matrix elements from
!>       matrix to A%elms, it will only associate the pointer with the array
!>       matrix. \n
!> BUT type(Matrix) :: A,B; A = B SHOULD copy the matrix elements from matrix B to A 
!>     (see mat_assign). \n
!> ALWAYS and ONLY call mat_free on a matrix you have initialized with mat_init.
!>
MODULE matrix_operations_aux
!FIXME: order routines alphabetically
!   use lstiming
  use precision
  use lstiming
  use matrix_operations
   use matrix_module
!   Use matrix_operations_symm_dense
   use matrix_operations_dense
   use matrix_operations_scalapack
#ifdef VAR_ENABLE_TENSORS
   use matrix_operations_pdmm
#endif
   use matrix_operations_csr
!   Use matrix_operations_unres_symm_dense
   use matrix_op_unres_dense
   
   private
   public :: mat_density_from_orbs, mat_eigenvalues_to_aux, &
        & mat_set_from_full2, mat_to_full2, mat_TrAB_ab, mat_ab_daxpy,&
        & mat_section2, mat_mix_homolumo, mat_precond, mat_mo_precond,&
        & mat_new_mo_precond, mat_new_complex_precond, mat_mo_precond_complex,&
        & mat_ao_precond_fallback, mat_ao_precond, mat_create_elm, &
        & mat_create_ab_elms, mat_get_elm, mat_get_ab_elms, mat_zerohalf,&
        & mat_write_to_disk2, mat_read_from_disk2, mat_VEC_TO_MAT,&
        & MAT_TO_VEC, mat_report_sparsity, mat_inquire_cutoff, mat_dmul,&
        & mat_hmul, mat_hdiv, mat_dger, mat_dhmul, mat_zero_cutoff,&
        & mat_dE_dmu, mat_column_norm, mat_sum, mat_Condition_Number
   contains

SUBROUTINE Mat_Condition_Number(A,ConditionNumber)
  implicit none
  TYPE(Matrix), INTENT(in) :: A
  Real(realk), INTENT(inout) :: ConditionNumber
  !
  integer :: nrow,ncol
  TYPE(Matrix) :: A2
  Real(realk) :: minEigv,maxEigv
  Real(realk),pointer :: eival(:) 
  nrow = A%nrow
  ncol = A%ncol
  IF(nrow.NE.ncol)then
     call lsquit('MAT_CONDITION_NUMBER Require a square matrix',-1)
  ENDIF
  call mem_alloc(eival,nrow)
  call mat_init(A2,nrow,ncol)
  call mat_assign(A2,A)
  call mat_dsyev(A2,eival,nrow)
  call mat_free(A2)
  minEigV = MINVAL(eival)
  maxEigV = MAXVAL(eival)
  call mem_dealloc(eival)
  conditionNumber = abs(maxEigV)/abs(minEigV)
END SUBROUTINE MAT_CONDITION_NUMBER

!> \brief Make c = a(1:ndim,1:nocc)*a^T(1:nocc,1:ndim) where a and c are type(matrix) 
!> \author T. Kjærgaard
!> \date 2012
!> \param a The first type(matrix) factor
!> \param c The output type(matrix)
SUBROUTINE mat_density_from_orbs(a, c,nocc,nocca,noccb,orbfree)
  !c = a(n,m)*a^T(m,n)
  implicit none
  TYPE(Matrix), intent(IN) :: a
  TYPE(Matrix), intent(inout):: c
  integer,intent(in) :: nocc
  integer,optional :: nocca,noccb
  logical,optional :: orbfree
  logical :: doOrbFree

  doOrbFree = .false.
  IF (present(orbfree)) doOrbFree = orbFree

  IF (doOrbFree) THEN
!   Special case for orbital free DFT. Only the lowest occopied orbital is included, with occupation number N
    call mat_density_from_orbs_nocc(a, c,1,1,1)
    call mat_scal(1.e0_realk*nocc,c)
  ELSE
!   Regular case
    call mat_density_from_orbs_nocc(a, c,nocc,nocca,noccb)
  ENDIF
END SUBROUTINE mat_density_from_orbs

!> \brief Make c = a(1:ndim,1:nocc)*a^T(1:nocc,1:ndim) where a and c are type(matrix) 
!> \author T. Kjærgaard
!> \date 2012
!> \param a The first type(matrix) factor
!> \param c The output type(matrix)
SUBROUTINE mat_density_from_orbs_nocc(a, c,nocc,nocca,noccb)
  !c = a(n,m)*a^T(m,n)
  implicit none
  TYPE(Matrix), intent(IN) :: a
  TYPE(Matrix), intent(inout):: c
  integer,intent(in) :: nocc
  integer,optional :: nocca,noccb
  integer :: ndim
  real(realk), allocatable :: orb(:),den(:),orb_ab(:,:),den_ab(:,:)
  ndim=a%nrow
  call time_mat_operations1
  no_of_matmuls = no_of_matmuls + 1
  select case(matrix_type)
  case(mtype_unres_dense)
     if(.not.present(nocca))call lsquit('no nocca',-1)
     if(.not.present(noccb))call lsquit('no noccb',-1)
     allocate(orb_ab(ndim*ndim,2),den_ab(ndim*ndim,2))
     call mat_to_full2(A,1E0_realk,orb_ab)
     CALL DGEMM('N','T',ndim,ndim,nocca,1E0_realk,orb_ab(:,1),ndim,orb_ab(:,1),ndim, &
          & 0E0_realk,den_ab(:,1),ndim)
     CALL DGEMM('N','T',ndim,ndim,noccb,1E0_realk,orb_ab(:,2),ndim,orb_ab(:,2),ndim, &
          & 0E0_realk,den_ab(:,2),ndim)
     call mat_set_from_full2(den_ab,1E0_realk,C)
     deallocate(orb_ab,den_ab)
  case(mtype_dense)
     call DGEMM('N','T',ndim,ndim,nocc,1E0_realk,&
          &A%elms,A%nrow,A%elms,A%nrow,0E0_realk,C%elms,C%nrow)
  case(mtype_scalapack)
     call mat_scalapack_density_from_orbs(a,c,ndim,nocc)
  case default
     allocate(den(ndim*ndim),orb(ndim*ndim))
     call mat_to_full(A,1E0_realk,orb)
     CALL DGEMM('N','T',ndim,ndim,nocc,1E0_realk,orb,ndim,orb,ndim, &
          & 0E0_realk,den,ndim)
     call mat_set_from_full(den,1E0_realk,C)
     deallocate(orb,den)
  end select
  call time_mat_operations2(JOB_mat_density_from_orbs)
  
END SUBROUTINE mat_density_from_orbs_nocc

!> \brief Compute maximum and minimum eigenvalue of matrix A
!> \author B. Jansik
!> \date March 2006
!> 
!> Store the results to A%raux(1) for emax and A%raux(2) for emin
!>
SUBROUTINE mat_eigenvalues_to_aux(unres,A)
  implicit none
   !> True if unrestricted calculation
   logical, intent(in)      :: unres
   !> Calculate max and min eigenvalue of this matrix
   Type(Matrix),intent(inout) :: A !Output because eigenvalues are stored in A%raux
   Type(Matrix)             :: I,Ix,IT
   real(realk), allocatable :: T(:,:), eval(:), wrk(:),TMPE(:)
   integer,     allocatable :: iwrk(:) 
   integer                  :: info, n, n2, lwrk, ilaenv, m
   real(realk)              :: NULLR,Z(1,1),eival
   integer,     allocatable :: NULLI(:)
   info=0
   NULLR = 0.0E0_realk
   if (.not.associated(A%raux)) allocate(A%raux(2))
   select case(matrix_type)
   case(mtype_unres_dense)
      n = 2*A%ncol; n2=n*n          
      !get optimal size of wrk array
      lwrk = (ilaenv(1,'DSYTRD','U',n,n,n,n)+3)*n
      lwrk = max(lwrk,8*n)
      allocate(T(n,n),wrk(lwrk),iwrk(5*n),TMPE(n))
      allocate(NULLI(n))
      !compute extremal eigenvalues
      !compute max eigenvalue
      call mat_to_full(A,1.0E0_realk,T)
      call dsyevx('N','I','U',n,T,n,NULLR,NULLR,n,n,0.0E0_realk,m,TMPE,&
           & Z,1,wrk,lwrk,iwrk,NULLI,info)
      A%raux(1)=TMPE(1)
      if (info .ne. 0) then
         write(*,*) &
              & 'DSYEVX routine failed to compute max eigenvalue. INFO=',info, &
              & 'See man dsyevx for more information.'
         call lsquit('DSYEVX failed in mat_eigenvalues_to_aux',-1)
      endif
      
      !compute min eigenvalue
      call mat_to_full(A,1.0E0_realk,T)
      
      call dsyevx('N','I','U',n,T,n,NULLR,NULLR,1,1,0.0E0_realk,m,TMPE,&
           & Z,1,wrk,lwrk,iwrk,NULLI,info)
      A%raux(2) = TMPE(1)
      deallocate(T,iwrk,wrk,TMPE) 
      deallocate(NULLI)
      if (info .ne. 0) then
         write(*,*) &
              & 'DSYEVX routine failed to compute min eigenvalue. INFO=',info, &
              & 'See man dsyevx for more information.'
         call lsquit('DSYEVX failed in mat_eigenvalues_to_aux',-1)
      endif
      !      print *, 'emax: ', A%raux(1), ' emin: ', A%raux(2)          
   case(mtype_dense)
      call mat_init(I,A%nrow,A%ncol)
      call mat_assign(I,A)
      !symmetrize 
      call mat_init(Ix,I%nrow,I%ncol)
      call mat_init(IT,I%nrow,I%ncol)
      call mat_assign(Ix,I)
      call mat_trans(I,IT)
      call mat_add(0.5E0_realk,Ix,0.5E0_realk,IT,I)
      call mat_free(Ix)
      call mat_free(IT)

      !mat_dsyevx changes TMP
      !compute extremal eigenvalues
      !compute max eigenvalue
      n = A%ncol
      call mat_dsyevx(I,eival,n)
      A%raux(1) = eival
      !compute min eigenvalue
      n = 1
      call mat_dsyevx(I,eival,n)
      A%raux(2) = eival
      call mat_free(I)
   case default
      call lsquit('mat_eigenvalues_to_aux not implemted for this type',-1)
   end select
  
END SUBROUTINE mat_eigenvalues_to_aux

!> \brief Unrestricted only! Convert a standard fortran matrix to a type(matrix) - USAGE DISCOURAGED!
!> \author L. Thogersen
!> \date 2003
!> \param afull (2n x 2n) standard fortran matrix that should be converted
!> \param alpha The output type(matrix) is multiplied by alpha
!> \param a The output (n x n) type(matrix)
SUBROUTINE mat_set_from_full2(afull,alpha, a)
  implicit none
  real(realk), INTENT(IN) :: afull(*)
  real(realk), intent(in) :: alpha
  TYPE(Matrix)            :: a  !output
  select case(matrix_type)
  case(mtype_unres_dense)
     call mat_unres_dense_set_from_full2(afull,alpha,a)
  case default
     call lsquit("mat_set_from_full2 not implemented for this type of matrix",-1)
  end select
END SUBROUTINE mat_set_from_full2

!> \brief Unrestricted only! Convert a type(matrix) to a standard fortran matrix - USAGE DISCOURAGED!
!> \author L. Thogersen
!> \date 2003
!> \param a The (n x n) type(matrix) that should be converted
!> \param afull The (2n x 2n) output standard fortran matrix 
!> \param alpha The output standard matrix is multiplied by alpha
SUBROUTINE mat_to_full2(a, alpha, afull)
  implicit none
  TYPE(Matrix), intent(in):: a
  real(realk), intent(in) :: alpha
  real(realk), intent(out):: afull(*)  !output  
  select case(matrix_type)
  case(mtype_unres_dense)
     call mat_unres_dense_to_full2(a, alpha, afull)
  case default
     call lsquit("mat_to_full2 not implemented for this type of matrix",-1)
  end select
END SUBROUTINE mat_to_full2
   
!> \brief Traces of alpha- and beta-part of matrix-product AB
!> \author C. Nygaard
!> \date 2010-07-02
!> \param A The first matrix
!> \param B The second matrix
!> \param trace The trace of the alpha- and beta-part of the matrix-product
subroutine mat_TrAB_ab (A, B, trace)
implicit none
type(matrix), intent(in) :: A, B
real(realk), intent(out) :: trace(:)

if (size(trace) == 1) then
  trace(1) = mat_TrAB (A, B)
elseif (size(trace) == 2) then
  select case (matrix_type)
  case (mtype_dense)
    trace(1) = mat_dense_TrAB(A,B)
    trace(2) = 0.0E0_realk
    print *, 'Warning: mat_TrAB_ab is used with mtype_dense, trace(2) = 0'
  case (mtype_unres_dense)
    call mat_unres_dense_TrAB_ab (A, B, trace)
  case default
    call lsquit ('mat_TrAB_ab not implemented for this type of matrix',-1)
  end select
else
  call lsquit ('Wrong dimension of trace in mat_TrAB_ab',-1)
endif
end subroutine mat_TrAB_ab

!> \brief mat_daxpy if alpha is different for the two matrix parts (unres)
!> \author C. Nygaard
!> \date 2010
!> \param alpha The alpha parameters
!> \param X The input type(matrix) 
!> \param Y The input/output type(matrix)
subroutine mat_ab_daxpy (alpha, X, Y)
implicit none
real(realk), intent(in)     :: alpha(:)
type(matrix), intent(in)    :: X
type(matrix), intent(inout) :: Y

if (size(alpha) == 1) then
  call mat_daxpy (alpha(1), X, Y)
elseif (size(alpha) == 2) then
  select case (matrix_type)
  case (mtype_unres_dense)
    call mat_unres_dense_ab_daxpy (alpha, X, Y)
  case default
    stop 'mat_ab_daxpy only works for mtype_unres_dense'
  end select
else
  stop 'Wrong dimension of alpha in mat_ab_daxpy'
endif

end subroutine mat_ab_daxpy

!> \brief Computes dE_SCF/dmu where mu is the damping in damped roothan
!> \author L. Thogersen
!> \date 2003
!> \param Fnew New Fock/Kohn-Sham matrix
!> \param Fdamp Damped Fock/Kohn-Sham matrix
!> \param SDS S*D*S, S = overlap matrix, D = density matrix
!> \param Cmo C coefficients
!> \param nocc Number of occupied orbitals
!> \return dE_SCF/dmu
FUNCTION mat_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
  !Find dE/dmu = sum_nu,I dE/dC_nu,I dC_nu,I/dmu
  implicit none
  type(matrix),intent(in) :: Fnew,Fdamp,SDS,Cmo
  integer, intent(in)     :: nocc
  real(realk)             :: mat_dE_dmu
  
  select case(matrix_type)
     !         case(mtype_symm_dense)
     !             mat_dE_dmu = mat_symm_dense_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
  case(mtype_dense)
     mat_dE_dmu = mat_dense_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
!         case(mtype_unres_symm_dense)
!             mat_dE_dmu = mat_unres_symm_dense_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
!         case(mtype_unres_dense)
!             mat_dE_dmu = mat_unres_dense_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
         case default
              call lsquit("mat_dE_dmu not implemented for this type of matrix",-1)
         end select
      END FUNCTION mat_dE_dmu
 
!> \brief Returns the sum of the elements Mat(from_row:to_row,ncol) squared
!> \author L. Thogersen
!> \date 2003
!> \param Mat Input type(matrix)
!> \param from_row Begin at this row
!> \param to_row End at this row
!> \param ncol Number of column to use
      FUNCTION mat_column_norm(Mat,ncol,from_row,to_row)
         implicit none
         type(Matrix), intent(in) :: Mat
         integer, intent(in) :: ncol, from_row, to_row
         real(realk) :: mat_column_norm

         if (to_row > Mat%nrow .or. from_row < 1 .or. ncol < 1 .or. Mat%ncol < ncol) then
           STOP 'wrong dimensions in mat_column_norm'
         endif
         if (from_row > to_row) then
           STOP 'from_row > to_row in mat_column_norm'
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             mat_column_norm = mat_symm_dense_column_norm(Mat,ncol,from_row,to_row)
         case(mtype_dense)
             mat_column_norm = mat_dense_column_norm(Mat,ncol,from_row,to_row)
!         case(mtype_unres_symm_dense)
!             mat_column_norm = mat_unres_symm_dense_column_norm(Mat,ncol,from_row,to_row)
         case(mtype_unres_dense)
             mat_column_norm = mat_unres_dense_column_norm(Mat,ncol,from_row,to_row)
         case default
              call lsquit("mat_column_norm not implemented for this type of matrix",-1)
         end select
      END FUNCTION mat_column_norm

!> \brief Puts smaller matrix A in selected position in larger matrix B
!> \author L. Thogersen
!> \date 2003
!> \param A The input type(matrix)
!> \param from_row Begin at this row in B
!> \param to_row End at this row in B
!> \param from_col Begin at this column in B
!> \param to_col End at this column in B
!> \param B The output type(matrix)
      subroutine mat_section2(A,from_row,to_row,from_col,to_col,B)
         implicit none
         type(Matrix), intent(in) :: A
         integer, intent(in) :: from_row, to_row, from_col, to_col
         type(Matrix), intent(inout) :: B  !output

         if (to_row > B%nrow .or. from_row < 1 .or. from_col < 1 .or. B%ncol < to_col) then
           STOP 'wrong dimensions in mat_section2'
         endif
         if (from_row > to_row .or. from_col > to_col) then
           STOP 'from_row or from_col > to_row or to_col in mat_section2'
         endif
         if ((to_row - from_row + 1) /= A%nrow .or. (to_col - from_col + 1) /= A%ncol) then
           STOP 'wrong dimensions in mat_section2'
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_section(A,from_row,to_row,from_col,to_col,Asec)
         case(mtype_dense)
             call mat_dense_section2(A,from_row,to_row,from_col,to_col,B)
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_section(A,from_row,to_row,from_col,to_col,Asec)
!         case(mtype_unres_dense)
!             call mat_unres_dense_section(A,from_row,to_row,from_col,to_col,Asec)
         case default
              call lsquit("mat_section2 not implemented for this type of matrix",-1)
         end select
      END SUBROUTINE mat_section2

!> \brief Mix two matrices into a plus and a minus combination
!> \author C. Nygaard
!> \date 2010
!> \param homo One of the matrices (+ in both combinations)
!> \param lumo The other matrix (different signs in the combinations)
subroutine mat_mix_homolumo (homo, lumo)

!Made for the unrestricted part,
! removes the symmetry of the starting guess
!Mix homo=(homo^alpha , homo^beta) and lumo=(lumo^alpha , lumo^beta)
! into newhomo=(plus^alpha , minus^beta)
!  and newlumo=(minus^alpha , plus^beta)
! where plus=homo+lumo and minus=homo-lumo

implicit none

type(matrix), intent(inout) :: homo, lumo
type(matrix)                :: plus, minus
integer                     :: r, c
real(realk)                 :: Na, Nb, ya, yb

r = homo%nrow ; c = homo%ncol

if (lumo%nrow /= r .or. lumo%ncol /= c) then
  write (mat_lu, *) 'Different dimensions of homo and lumo in mat_mix_homolumo!'
  call lsquit ('Wrong dimensions in mat_mix_homolumo',mat_lu)
endif

call mat_init (plus, r, c)
call mat_init (minus, r, c)

ya = 0.1E0_realk ; yb = 0.05E0_realk
Na = 1.0E0_realk/SQRT((1.0E0_realk-ya)**2 + ya**2)
Nb = 1.0E0_realk/SQRT((1.0E0_realk-yb)**2 + yb**2)

call mat_add (Na*(1.0E0_realk-ya), homo, Na*ya, lumo, plus)
call mat_add (Nb*(1.0E0_realk-yb), homo, -Nb*yb, lumo, minus)

select case (matrix_type)
  case (mtype_unres_dense)
    call mat_unres_dense_mix_homolumo (plus, minus, homo, lumo)
  case default
    stop 'mat_mix_homolumo not implemented for this type of matrix'
end select

call mat_free (plus)
call mat_free (minus)

end subroutine mat_mix_homolumo

!> \brief Preconditioning of vector x by matrix M: xprec(i) = x(i)/M(i,i)
!> \author S. Host
!> \date 2005
!> \param M The preconditioner
!> \param x Input vector to be preconditioned
!> \param xprec Preconditioned output vector
      subroutine mat_precond(M,x,xprec)
         implicit none
         type(Matrix), intent(inout) :: xprec
         type(Matrix), intent(in)    :: M, x

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_precond(M,x,xprec)
         case(mtype_dense)
             call mat_dense_precond(M,x,xprec)
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_precond(M,x,xprec)
         case(mtype_unres_dense)
             call mat_unres_dense_precond(M,x,xprec)
         case default
            call lsquit("mat_precond not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MATPRE',mat_TSTR,mat_TEN,mat_lu)
      END SUBROUTINE mat_precond

!> \brief Divide the occ-virt part of the matrix X_MO with the orbital energy difference E_a - E_i
!> \author L. Thogersen
!> \date 2005
!> \param nocc Number of occupied orbitals
!> \param omega Levelshift
!> \param Eorb_final Final orbital energies
!> \param X_MO Input/output - matrix to be preconditioned
      subroutine mat_mo_precond(nocc,omega,Eorb_final,X_MO)
         implicit none
         integer, intent(in) :: nocc
         real(realk), intent(in) :: omega
         real(realk), intent(in) :: Eorb_final(:)
         type(Matrix), intent(inout) :: X_MO
         integer :: i, j
         real(realk), ALLOCATABLE :: XMO_full(:,:)
         real(realk) :: dia

         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_mo_precond(nocc,omega,Eorb_final,X_MO)
         case(mtype_dense)
             call mat_dense_mo_precond(nocc,omega,Eorb_final,X_MO)
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_mo_precond(nocc,omega,Eorb_final,X_MO)
         case(mtype_unres_dense)
             call mat_unres_dense_mo_precond(nocc,omega,Eorb_final,X_MO)
         case default
            print *, "FALLBACK mo_precond"
            ! see mat_dense_mo_precond
            allocate(xmo_full(X_MO%nrow, X_MO%ncol))
            call mat_to_full(X_MO,1E0_realk,xmo_full)
            do j = 1,nocc  !columns
               do i = nocc+1,X_MO%nrow  !rows
                  dia = 2.0E0_realk*(Eorb_final(i) - Eorb_final(j)) + omega * 2E0_realk
                  xmo_full(i,j) = xmo_full(i,j) /dia
               enddo
            enddo
            !modify upper right block (ia block) 
            do j = nocc+1,X_MO%ncol  !columns
               do i = 1,nocc          !rows
                  !dia = E[2]_dia + omega*S[2]_dia
                  dia = 2.0E0_realk*(Eorb_final(j) - Eorb_final(i)) - omega * 2E0_realk
                  xmo_full(i,j) = xmo_full(i,j)/dia
               enddo
            enddo
            call mat_set_from_full(xmo_full, 1E0_realk, X_MO)
            DEALLOCATE(XMO_full)
         end select
      END SUBROUTINE mat_mo_precond

!vide the occ-virt part of the matrix X_MO with the orbital energy difference E_a - E_i
 !> \author J. Kauczor
 !> \date 2010
 !> \param nocc Number of occupied orbitals
 !> \param omega Levelshift
 !> \param Eorb_final Final orbital energies
 !> \param X_MO Input/output - matrix to be preconditioned
       subroutine mat_new_mo_precond(nocc,omega,Eorb_final,Xp_MO,xm_mo)
          implicit none
          integer, intent(in) :: nocc
          real(realk), intent(in) :: omega
          real(realk), intent(in) :: Eorb_final(:)
          type(Matrix), intent(inout) :: Xp_MO,xm_mo
          integer :: i, j,ndim
          real(realk), ALLOCATABLE :: XpMO_full(:,:),XmMO_full(:,:)
          real(realk) :: dia,aa
 
          ndim=xp_mo%nrow
          select case(matrix_type)
          case(mtype_dense)
              call mat_dense_new_mo_precond(nocc,omega,Eorb_final,Xp_MO,xm_mo)
          case default
             print *, "FALLBACK mo_precond"
             ! see mat_dense_mo_precond
             allocate(xpmo_full(Xp_MO%nrow, Xp_MO%ncol))
             allocate(xmmo_full(Xm_MO%nrow, Xm_MO%ncol))
             call mat_to_full(Xp_MO,1E0_realk,xpmo_full)
             call mat_to_full(Xm_MO,1E0_realk,xmmo_full)
             do j = 1,nocc  !columns
                do i = nocc+1,Xp_MO%nrow  !rows
                   aa=2E0_realk*(Eorb_final(i) - Eorb_final(j))*(Eorb_final(i)- Eorb_final(j))-2E0_realk*omega*omega
                 !  dia = 2.0E0_realk*(Eorb_final(i) - Eorb_final(j)) + omega * 2E0_realk
                   dia=xpmo_full(i,j)
                   xpmo_full(i,j) = ((Eorb_final(i) - Eorb_final(j))*xpmo_full(i,j)-omega*xmmo_full(i,j)) /aa
                   xmmo_full(i,j) = ((Eorb_final(i) - Eorb_final(j))*xmmo_full(i,j)-omega*dia) /aa
                enddo
             enddo
             !modify upper right block (ia block) 
             do j = nocc+1,Xp_MO%ncol  !columns
                do i = 1,nocc          !rows
                   !dia = E[2]_dia + omega*S[2]_dia
                   aa=2E0_realk*(Eorb_final(i) - Eorb_final(j))*(Eorb_final(i)-Eorb_final(j))-2E0_realk*omega*omega
                   !dia = 2.0E0_realk*(Eorb_final(j) - Eorb_final(i)) - omega * 2E0_realk
                   dia=xpmo_full(i,j)
                   xpmo_full(i,j) = ((Eorb_final(j) - Eorb_final(i))*xpmo_full(i,j)+omega*xmmo_full(i,j)) /aa
                   xmmo_full(i,j) = ((Eorb_final(j) - Eorb_final(i))*xmmo_full(i,j)+omega*dia) /aa
                  ! xmo_full(i,j) = xmo_full(i,j)/dia
                enddo
             enddo
             call mat_set_from_full(xpmo_full, 1E0_realk, Xp_MO)
             call mat_set_from_full(xmmo_full, 1E0_realk, Xm_MO)
             DEALLOCATE(XpMO_full)
             DEALLOCATE(XmMO_full)
          end select
       END SUBROUTINE mat_new_mo_precond
 
 
 !> \brief Divide the occ-virt part of the matrix X_MO with the orbital energy difference E_a - E_i (complex orbitals)
 !> \author J. Kauczor
 !> \date 2010
 !> \param nocc Number of occupied orbitals
 !> \param Eorb_final Final orbital energies
 !> \param X_MO Input/output - matrix to be preconditioned
       subroutine mat_new_complex_precond(nocc,omega,gammma,Eorb_final,Xp_MO,xm_mo,xpi_mo,xmi_mo)
          implicit none
          integer, intent(in) :: nocc
          real(realk), intent(in) :: omega,gammma
          real(realk), intent(in) :: Eorb_final(:)
          type(Matrix), intent(inout) :: Xp_MO,xm_mo,xpi_mo,xmi_mo
          integer :: i, j,ndim
          real(realk), ALLOCATABLE :: XpMO_full(:,:),XmMO_full(:,:)
          real(realk), ALLOCATABLE :: XpiMO_full(:,:),XmiMO_full(:,:)
          real(realk), ALLOCATABLE :: Xp(:,:),Xm(:,:)
          real(realk), ALLOCATABLE :: Xpi(:,:),Xmi(:,:)
          real(realk) :: dia,aa,ab,a,b,c,d
 
          ndim=xp_mo%nrow
          select case(matrix_type)
          case(mtype_dense)
              call mat_dense_new_complex_precond(nocc,omega,gammma,Eorb_final,Xp_MO,xm_mo,xpi_mo,xmi_mo)
          case default
             print *, "FALLBACK mo_precond"
             ! see mat_dense_mo_precond
             allocate(xpmo_full(Xp_MO%nrow, Xp_MO%ncol))
             allocate(xmmo_full(Xm_MO%nrow, Xm_MO%ncol))
             allocate(xpimo_full(Xpi_MO%nrow, Xpi_MO%ncol))
             allocate(xmimo_full(Xmi_MO%nrow, Xmi_MO%ncol))
             allocate(xp(Xp_MO%nrow, Xp_MO%ncol))
             allocate(xm(Xm_MO%nrow, Xm_MO%ncol))
             allocate(xpi(Xpi_MO%nrow, Xpi_MO%ncol))
             allocate(xmi(Xmi_MO%nrow, Xmi_MO%ncol))
             xp=0E0_realk; xm=0E0_realk; xpi=0E0_realk; xmi=0E0_realk
             call mat_to_full(Xp_MO,1E0_realk,xpmo_full)
             call mat_to_full(Xm_MO,1E0_realk,xmmo_full)
             call mat_to_full(Xpi_MO,1E0_realk,xpimo_full)
             call mat_to_full(Xmi_MO,1E0_realk,xmimo_full)
             do i=1,Xp_MO%nrow
                xp(i,i)=xpmo_full(i,i)
                xpi(i,i)=xpimo_full(i,i)
               ! xm(i,i)=xpmo_full(i,i)
               ! xmi(i,i)=xpimo_full(i,i)
             enddo
             do j = 1,nocc  !columns
                do i = nocc+1,Xp_MO%nrow  !rows
                   ab=(Eorb_final(i) - Eorb_final(j))*(Eorb_final(i)-Eorb_final(j))-(omega*omega-gammma*gammma) 
                   aa=2E0_realk*ab*ab+8E0_realk*(omega*omega*gammma*gammma)
                 !  dia=xpmo_full(i,j)
                   a=(Eorb_final(i) - Eorb_final(j))*ab
                   b=-omega*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
                      &-(omega*omega+gammma*gammma))
                   c=-gammma*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
                      &+(omega*omega+gammma*gammma))
                   d=2E0_realk*omega*gammma*(Eorb_final(i) - Eorb_final(j))
                   
                   xp(i,j)  = (a*xpmo_full(i,j)+b*xmmo_full(i,j)-d*xpimo_full(i,j)-c*xmimo_full(i,j)) /aa
                   xm(i,j)  = (b*xpmo_full(i,j)+a*xmmo_full(i,j)-c*xpimo_full(i,j)-d*xmimo_full(i,j))/aa
                   xpi(i,j) = (d*xpmo_full(i,j)+c*xmmo_full(i,j)+a*xpimo_full(i,j)+b*xmimo_full(i,j))/aa
                   xmi(i,j) = (c*xpmo_full(i,j)+d*xmmo_full(i,j)+b*xpimo_full(i,j)+a*xmimo_full(i,j))/aa
                enddo
             enddo
             !modify upper right block (ia block) 
             do j = nocc+1,Xp_MO%ncol  !columns
                do i = 1,nocc          !rows
        ab=(Eorb_final(i) - Eorb_final(j))*(Eorb_final(i)-Eorb_final(j))-(omega*omega-gammma*gammma) 
         aa=2E0_realk*ab*ab+8E0_realk*(omega*omega*gammma*gammma)
         
         a=(Eorb_final(j) - Eorb_final(i))*ab
         b=omega*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
         &-(omega*omega+gammma*gammma))
         c=gammma*((Eorb_final(i) - Eorb_final(j))*(Eorb_final(i) - Eorb_final(j))&
         &+(omega*omega+gammma*gammma))
         d=2E0_realk*omega*gammma*(Eorb_final(j) - Eorb_final(i))
                   
                   xp(i,j)  = (a*xpmo_full(i,j)+b*xmmo_full(i,j)-d*xpimo_full(i,j)-c*xmimo_full(i,j)) /aa
                   xm(i,j)  = (b*xpmo_full(i,j)+a*xmmo_full(i,j)-c*xpimo_full(i,j)-d*xmimo_full(i,j))/aa
                   xpi(i,j) = (d*xpmo_full(i,j)+c*xmmo_full(i,j)+a*xpimo_full(i,j)+b*xmimo_full(i,j))/aa
                   xmi(i,j) = (c*xpmo_full(i,j)+d*xmmo_full(i,j)+b*xpimo_full(i,j)+a*xmimo_full(i,j))/aa
                enddo
             enddo
     
             call mat_set_from_full(xp, 1E0_realk, Xp_MO)
             call mat_set_from_full(xm, 1E0_realk, Xm_MO)
             call mat_set_from_full(xpi, 1E0_realk, Xpi_MO)
             call mat_set_from_full(xmi, 1E0_realk, Xmi_MO)
             DEALLOCATE(XpMO_full)
             DEALLOCATE(XmMO_full)
             DEALLOCATE(XpiMO_full)
             DEALLOCATE(XmiMO_full)
             DEALLOCATE(Xp)
             DEALLOCATE(Xm)
             DEALLOCATE(Xpi)
             DEALLOCATE(Xmi)
          end select
       END SUBROUTINE mat_new_complex_precond

!> \author J. Kauczor
!> \date 2009
!> \param nocc Number of occupied orbitals
!> \param Eorb_final Final orbital energies
!> \param X_MO Input/output - matrix to be preconditioned
      subroutine mat_mo_precond_complex(nocc,Eorb_final,X_MO)
         implicit none
         integer, intent(in) :: nocc
         real(realk), intent(in) :: Eorb_final(:)
         type(Matrix), intent(inout) :: X_MO
         integer :: i, j,ndim
         real(realk), ALLOCATABLE :: XMO_full(:,:)
         real(realk) :: dia

         ndim=X_MO%nrow
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_mo_precond_complex(nocc,Eorb_final,X_MO)
         case default
            print *, "FALLBACK mo_precond"
            ! see mat_dense_mo_precond
            allocate(xmo_full(X_MO%nrow, X_MO%ncol))
            call mat_to_full(X_MO,1E0_realk,xmo_full)
            do j = 1,nocc  !columns
               do i = nocc+1,X_MO%nrow  !rows
                  dia = 2.0E0_realk*(Eorb_final(i) - Eorb_final(j))
                  xmo_full(i,j) = xmo_full(i,j) /dia
               enddo
            enddo
            !modify upper right block (ia block) 
            do j = nocc+1,X_MO%ncol  !columns
               do i = 1,nocc          !rows
                  !dia = E[2]_dia + omega*S[2]_dia
                  dia = 2.0E0_realk*(Eorb_final(j) - Eorb_final(i))
                  xmo_full(i,j) = xmo_full(i,j)/dia
               enddo
            enddo
            call mat_set_from_full(xmo_full, 1E0_realk, X_MO)
            DEALLOCATE(XMO_full)
         end select
      END SUBROUTINE mat_mo_precond_complex

!> \brief Diagonal preconditioning in orthonormal AO basis 
!> \author S. Host
!> \date 2005
!> \param symm Symmetry indicator: symmetric = 1, antisymmetric = 2, nonsymmetric = 0
!> \param omega Level shift
!> \param FUP Fock/KS matrix in OAO basis, occupied part (virtual part projected out)
!> \param FUQ Fock/KS matrix in OAO basis, virtual part (occupied part projected out)
!> \param DU Density matrix in OAO basis
!> \param X_AO Matrix to be preconditioned
!> 
!> Preconditioning with orbital energy difference. This is done by using the occupied and virtual 
!> parts of the Fock/KS matrix: 
!> X_prec(i,j) = X(i,j) / [FUQ(j,j)-FUP(j,j) + FUQ(i,i)-FUP(i,i) - omega*(DU(j,j)-DU(i,i))]
!> taking care not to divide by zero and exploiting symmetry if X is symm or antisymm
!> 
      subroutine mat_ao_precond(symm,omega,FUP,FUQ,DU,X_AO)
         implicit none
         integer, intent(in) :: symm
         real(realk), intent(in) :: omega
         type(Matrix), intent(in) :: FUP, FUQ,DU
         type(Matrix), intent(inout) :: X_AO

         if (info_memory) write(mat_lu,*) 'Before mat_ao_precond: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
         case(mtype_dense)
            call mat_dense_ao_precond(symm,omega,FUP,FUQ,DU,X_AO)
         case(mtype_csr)
            call mat_csr_ao_precond(symm,omega,FUP,FUQ,DU,X_AO)
         case(mtype_unres_dense)
            call mat_unres_dense_ao_precond(symm,omega,FUP,FUQ,DU,X_AO)
         case(mtype_scalapack)
            call mat_scalapack_ao_precond(symm,omega,FUP,FUQ,DU,X_AO)
         case default
            print *, "FALLBACK: mat_ao_precond"
            call mat_ao_precond_fallback
         end select
         if (info_memory) write(mat_lu,*) 'After mat_ao_precond: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('AOPREC',mat_TSTR,mat_TEN,mat_lu)
       contains
         subroutine mat_ao_precond_fallback
           implicit none
           real(realk), allocatable :: FPPd(:,:), FQQd(:,:)
           real(realk), allocatable :: Dd(:,:)
           real(realk), allocatable :: X(:,:)
           integer :: n, i, j
           real(realk) :: denom

           n = FUP%nrow
           allocate(FPPd(n,n),FQQd(n,n),Dd(n,n),X(n,n))
           call mat_to_full(FUP,1E0_realk,FPPd)
           call mat_to_full(FUQ,1E0_realk,FQQd)
           call mat_to_full(DU,1E0_realk,Dd)
           !call mat_to_full(SQQ,1E0_realk,SQQd)
           call mat_to_full(X_AO,1E0_realk,X)
              ! extract diagonals and use them instead of inefficient
              ! stride (n+1) access!
           if(symm .eq. 1 .or. symm .eq. 2) THEN
              DO i = 1, n
                 DO j = 1, n
                    denom = FQQd(j,j) + FQQd(i,i) &
                         &- FPPd(i,i) - FPPd(j,j) &
                         &- omega
                    IF(ABS(denom)>1E-9_realk) X(i,j) = X(i,j)/(denom)  !12/10-09
                                                                 !Stinne removed factor 2 to match matop_dense
                 END DO
              END DO
           ELSE
              DO i = 1, n
                 DO j = 1, n
                    denom = FQQd(j,j) + FQQd(i,i) &
                         &- FPPd(i,i) - FPPd(j,j) &
                         &- omega*(Dd(i,i)-Dd(j,j)) 
                    IF(ABS(denom)>1E-9_realk) X(i,j) = X(i,j)/(denom) !12/10-09
                                                                 !Stinne removed factor 2 to match matop_dense 
                 END DO
              END DO
           end if
           call mat_set_from_full(X,1E0_realk, X_AO)
           deallocate(FPPd,FQQd,Dd,X)
         end subroutine mat_ao_precond_fallback

      END SUBROUTINE mat_ao_precond


!> \brief Create or overwrite element in matrix
!> \author S. Host
!> \date 2005
!> \param i Row index for element to be created
!> \param j Column index for element to be created
!> \param val Value of element to be created
!> \param A Input/output matrix
      subroutine mat_create_elm(i,j,val,A)
         implicit none
         integer, intent(in) :: i,j
         real(Realk), intent(in) :: val
         type(Matrix), intent(inout) :: A

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (i > A%nrow .or. j > A%ncol .or. i < 0 .or. j < 0) then
           WRITE(mat_lu,*) 'cannot create element, the indexes',i,j,&
                        & 'are out of the bounds nrow,ncol ',A%nrow,A%ncol
           STOP 'cannot create element, the indexes are out of bounds'
         endif
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_create_elm(i,j,val,A)
         case(mtype_dense)
             call mat_dense_create_elm(i,j,val,A)
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_create_elm(i,j,val,A)
         case(mtype_unres_dense)
             call mat_unres_dense_create_elm(i,j,val,A)
         case default
              call lsquit("mat_create_elm not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('CREATE',mat_TSTR,mat_TEN,mat_lu)
      END SUBROUTINE mat_create_elm

!> \brief mat_create_elm for different alpha and beta parts
!> \author C. Nygaard
!> \date 2010
!> \param r Row index for element to be created
!> \param c Column index for element to be created
!> \param elm Values of elements to be created
!> \param A Input/output matrix
!=======================================================================
  subroutine mat_create_ab_elms (r, c, elm, A)

  implicit none

  integer, intent(in)         :: r, c
  real(realk), intent(in)     :: elm(:)
  type(matrix), intent(inout) :: A

  if (size(elm) == 1) then
    call mat_create_elm (r, c, elm(1), A)
  elseif (size(elm) == 2) then
    select case (matrix_type)
    case (mtype_unres_dense)
      call mat_unres_dense_create_ab_elms (r, c, elm, A)
    case (mtype_dense)
      call mat_dense_create_elm (r, c, elm(1), A)
    case default
      stop 'mat_create_ab_elm only works for unrestricted, you want&
           & mat_create_elm'
    end select
  else
    print *, 'expected dimension = 1 or 2, actual dimension =', size(elm)
    stop 'wrong dimension of elm in mat_create_ab_elms'
  endif

  end subroutine mat_create_ab_elms
!=======================================================================

!> \brief Get element from matrix
!> \author C. Nygaard
!> \date 2010
!> \param A Input matrix
!> \param r Row index for element
!> \param c Column index for element
!> \param elm Value of element (output)
!=======================================================================
      subroutine mat_get_elm (A, r, c, elm)

      implicit none

      type(matrix), intent(in) :: A
      integer, intent(in)      :: r, c
      real(realk), intent(out) :: elm
      real(realk)              :: tmp(2)

      select case (matrix_type)
      case (mtype_dense)
        call mat_dense_get_elm (A, r, c, elm)
      case (mtype_unres_dense)
        call mat_unres_dense_get_elm (A, r, c, tmp)
        elm = tmp(1)
      case default
        call lsquit("mat_get_elm not implemented for this matrix-type",-1)
      end select

      end subroutine mat_get_elm
!=======================================================================

!> \brief mat_get_elm for different alpha and beta parts
!> \author C. Nygaard
!> \date 2010
!> \param A Input matrix
!> \param r Row index for element
!> \param c Column index for element
!> \param elm Value of elements (output)
!=======================================================================
subroutine mat_get_ab_elms (A, r, c, elm)

implicit none

type(matrix), intent(in) :: A
integer, intent(in)      :: r, c
real(realk), intent(out) :: elm(:)

if (size(elm) == 1) then
  call mat_get_elm (A, r, c, elm(1))
elseif (size(elm) == 2) then
  select case (matrix_type)
  case (mtype_dense)
    call mat_get_elm (A, r, c, elm(1))
    elm(2) = 0.0E0_realk
  case (mtype_unres_dense)
    call mat_unres_dense_get_elm (A, r, c, elm)
  case default
    stop 'mat_get_ab_elm is only implemented for mtype_unres_dense'
  end select
else
  stop 'Wrong dimensions of elm in mat_get_ab_elms!'
endif

end subroutine mat_get_ab_elms

!> \brief Set upper or lower triangle of a type(matrix) to zero. Diagonal belongs to lower triangle.
!> \author L. Thogersen
!> \date 2003
!> \param part Indicates which part should be set to zero. 'UT'/'ut' = upper, 'LT'/'lt' = lower triangle
!> \param A Input/output matrix
subroutine mat_zerohalf(part,A)
  implicit none
  character(len=2), intent(in) :: part
  type(Matrix), intent(inout) :: A
  real(realk), allocatable :: Afull(:,:)
  integer  i, j
  
  !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
  if (part /= 'UT' .and. part /= 'ut' .and. part /= 'LT' .and. part /= 'lt') then
     STOP 'unknown part of matrix to zero'
  elseif (A%nrow /= A%ncol) then
     STOP 'cannot define triangles since nrow /= ncol'
  endif
  select case(matrix_type)
     !         case(mtype_symm_dense)
     !             call mat_symm_dense_zerohalf(part,A)
  case(mtype_dense)
     call mat_dense_zerohalf(part,A)
     !         case(mtype_unres_symm_dense)
     !             call mat_unres_symm_dense_zerohalf(part,A)
  case(mtype_unres_dense)
     call mat_unres_dense_zerohalf(part,A)
  case default
     PRINT *, "FALLBACK zerohalf ", part
     allocate(afull(a%nrow, a%ncol))
     call mat_to_full(a,1E0_realk,afull)
     if (part == 'ut' .or. part == 'UT') then
        !set the upper triangle to zero - diagonal is kept
        do j = 2,A%ncol
           afull(1:j-1,j) = 0E0_realk
        enddo
     else
        !set the lower triangle to zero - diagonal is also zeroed
        do j = 1,A%ncol
           afull(j:A%nrow,j) = 0E0_realk
        enddo
     endif
     call mat_set_from_full(afull, 1E0_realk, a)
     DEALLOCATE(afull)
  end select
  !if (INFO_TIME_MAT) CALL LSTIMER('ZEROHA',mat_TSTR,mat_TEN,mat_lu)
end subroutine mat_zerohalf

!> \brief Write a type(matrix) to disk in formatted form.
!> \author S. Host
!> \date 2007
!> \param iunit Logical unit number of file which matrix should be written to
!> \param A Matrix which should be written on disk
!>
!>  DEBUG ROUTINE (see debug_convert_density in debug.f90).
!>
subroutine mat_write_to_disk2(iunit,A)
  implicit none
  integer, intent(in) :: iunit
  type(Matrix), intent(in) :: A
  real(realk), allocatable :: afull(:,:)
  !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
  select case(matrix_type)
     !         case(mtype_symm_dense)
     !             call mat_symm_dense_write_to_disk(iunit,A)
  case(mtype_dense)
     call mat_dense_write_to_disk2(iunit,A)
     !         case(mtype_unres_symm_dense)
     !             call mat_unres_symm_dense_write_to_disk(iunit,A)
  case(mtype_unres_dense)
     call mat_unres_dense_write_to_disk(iunit,A)
  case default
     !print *, "FALLBACK: mat_write_to_disk"
     allocate(afull(a%nrow, a%ncol))
     call mat_to_full(a,1E0_realk,afull)
     write(iunit) A%Nrow, A%Ncol
     write(iunit) afull
     deallocate(afull)
     
  end select
  !if (INFO_TIME_MAT) CALL LSTIMER('WRITE ',mat_TSTR,mat_TEN,mat_lu)
end subroutine mat_write_to_disk2

!> \brief Read a type(matrix) from disk in formatted form.
!> \author S. Host
!> \date 2007
!> \param iunit Logical unit number of file from which matrix should be read
!> \param A Output matrix which should be read from disk
!>
!>  DEBUG ROUTINE (see debug_convert_density in debug.f90).
!>
subroutine mat_read_from_disk2(iunit,A)
  implicit none
  integer, intent(in) :: iunit
  type(Matrix), intent(inout) :: A  !output
  real(realk), allocatable :: afull(:,:)
  integer                  :: nrow, ncol
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
!         case(mtype_symm_dense)
!             call mat_symm_dense_read_from_disk(iunit,A)
         case(mtype_dense)
             call mat_dense_read_from_disk2(iunit,A)
!         case(mtype_unres_symm_dense)
!             call mat_unres_symm_dense_read_from_disk(iunit,A)
         case(mtype_unres_dense)
             call mat_unres_dense_read_from_disk(iunit,A)
         case default
            !print *, "FALLBACK: mat_read_from_disk"
            allocate(afull(a%nrow, a%ncol))
            READ(iunit) Nrow, Ncol
            if(Nrow /= A%nrow) stop 'mat_read_from_disk: Nrow /= A%nrow'
            if(Ncol /= A%ncol) stop 'mat_read_from_disk: Ncol /= A%ncol'
            read(iunit) afull
            call mat_set_from_full(afull,1E0_realk,a)
            deallocate(afull)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('READ  ',mat_TSTR,mat_TEN,mat_lu)
      end subroutine mat_read_from_disk2

!> \brief Change a vector to a matrix. The matrix must be symmetric or antisymmetric.
!> \author S. Host
!> \date 2005
!> \param symmetry Indicates symmetry of matrix, 's'/'S' = symmetric, 'a'/'A' = !antisymmetric
!> \param vec Input vector that should be transformed
!> \param mat Output matrix
!>
!> For debug purposes.
!> Symmetric: diagonal included, vecdim = matdim(matdim+1)/2
!> Antisym:   diagonal excluded, vecdim = matdim(matdim+1)/2 - matdim 
!>
    subroutine mat_vec_to_mat(symmetry, vec, mat)
    implicit none
                                                                                 
         integer :: n, m, i
         character, intent(in)      :: symmetry
         TYPE(matrix),intent(in)    :: vec
         TYPE(matrix),intent(inout) :: mat
     
    if  (symmetry /= 'a' .AND. symmetry /= 'A' .AND.       &
       & symmetry /= 's' .AND. symmetry /= 'S') then
          STOP 'Unknown symmetry possibility in mat_VEC_TO_MAT'
    endif

    if (symmetry == 's' .OR. symmetry == 'S') then     
       if (VEC%nrow /= MAT%nrow*(MAT%nrow+1)/2 .OR. &
         & VEC%nrow /= MAT%ncol*(MAT%ncol+1)/2 .OR. &
         & VEC%ncol /= 1) then
          STOP 'Wrong dimensions in mat_VEC_TO_MAT'
       endif
    endif

    if (symmetry == 'a' .OR. symmetry == 'A') then     
       if (VEC%nrow /= MAT%nrow*(MAT%nrow+1)/2 - MAT%nrow .OR. &
         & VEC%nrow /= MAT%ncol*(MAT%ncol+1)/2 - MAT%nrow .OR. &
         & VEC%ncol /= 1) then
          STOP 'Wrong dimensions in mat_VEC_TO_MAT'
       endif
    endif

    select case(matrix_type)
 !         case(mtype_symm_dense)
 !             call mat_symm_dense_vec_to_mat()
          case(mtype_dense)
              call mat_dense_vec_to_mat(symmetry, VEC, MAT)
#ifndef UNITTEST
#endif
 !         case(mtype_unres_symm_dense)
 !             call mat_unres_symm_dense_vec_to_mat()
          case(mtype_unres_dense)
              call mat_unres_dense_vec_to_mat(symmetry, VEC, MAT)
          case default
               CALL LSQUIT("mat_VEC_TO_MAT not implemented for this type of matrix",-1)
          end select
    end subroutine mat_VEC_TO_MAT
 
!> \brief Change a matrix to a vector. The matrix must be symmetric or antisymmetric.
!> \author S. Host
!> \date 2005
!> \param symmetry Indicates symmetry of matrix, 's'/'S' = symmetric, 'a'/'A' = !antisymmetric
!> \param mat Input matrix that should be transformed
!> \param vec Output vector
!>
!> For debug purposes.
!> Symmetric: diagonal included, vecdim = matdim(matdim+1)/2
!> Antisym:   diagonal excluded, vecdim = matdim(matdim+1)/2 - matdim 
!>
    subroutine mat_to_vec(symmetry, mat, vec)
    implicit none
                                    
         character, intent(in)      :: symmetry                                             
         TYPE(matrix),intent(inout) :: vec
         TYPE(matrix), intent(in)   :: mat

    if ( symmetry /= 'a' .AND. symmetry /= 'A' .AND.       &
       & symmetry /= 's' .AND. symmetry /= 'S') then
          STOP 'Unknown symmetry possibility in MAT_TO_VEC'
    endif

    if (symmetry == 's' .OR. symmetry == 'S') then     
       if (VEC%nrow /= MAT%nrow*(MAT%nrow+1)/2 .OR. &
         & VEC%nrow /= MAT%ncol*(MAT%ncol+1)/2 .OR. &
         & VEC%ncol /= 1) then
          STOP 'Wrong dimensions in MAT_TO_VEC'
       endif
    endif

    if (symmetry == 'a' .OR. symmetry == 'A') then     
       if (VEC%nrow /= MAT%nrow*(MAT%nrow+1)/2 - MAT%nrow .OR. &
         & VEC%nrow /= MAT%ncol*(MAT%ncol+1)/2 - MAT%nrow .OR. &
         & VEC%ncol /= 1) then
          STOP 'Wrong dimensions in MAT_TO_VEC'
       endif
    endif

    select case(matrix_type)
 !         case(mtype_symm_dense)
 !             call mat_symm_dense_mat_to_vec()
          case(mtype_dense)
              call mat_dense_mat_to_vec(symmetry, MAT, VEC)
#ifndef UNITTEST
#endif
 !         case(mtype_unres_symm_dense)
 !             call mat_unres_symm_dense_mat_to_vec()
          case(mtype_unres_dense)
              call mat_unres_dense_mat_to_vec(symmetry, MAT, VEC)
          case default
               CALL LSQUIT("MAT_TO_VEC not implemented for this type of matrix",-1)
          end select
    end subroutine MAT_TO_VEC

!> \brief Reports number of non-zero elements and sparsity of a type(matrix) A.
!> \author P. Salek
!> \date 2003
!> \param A Matrix for which number of non-zero elements i requested
!> \param mat_label The matrix is identified by this character string
!> \param iunit Logical unit number of file to which sparsity should be printed
!>
!> The data is printed to stream specified by iunit.
!> The matrix is identified by the mat_label string.
!>
    SUBROUTINE mat_report_sparsity(A,mat_label,nnz,iunit)
      implicit none
      TYPE(Matrix), INTENT(in) :: A
      CHARACTER(*), INTENT(IN) :: mat_label
      INTEGER, INTENT(in) :: iunit
      INTEGER, INTENT(out) :: nnz
      real(realk) :: sparsity
      
      if (info_memory) write(mat_lu,*) 'Before mat_report_sparsity: mem_allocated_global =', mem_allocated_global
      SELECT CASE(matrix_type)
      CASE(mtype_dense)
         CALL mat_dense_report_sparsity(A,sparsity)
      CASE(mtype_csr)
         CALL mat_csr_report_sparsity(A,sparsity)
      CASE default
         return
      END SELECT

#ifndef UNITTEST
      WRITE(iunit,'("Matrix ",A," has nnz=",I10," sparsity: ",F10.3," %")')&
           &mat_label, INT(sparsity*A%nrow*A%ncol), sparsity*100.0
#endif
      nnz = INT(sparsity*A%nrow*A%ncol)
      if (info_memory) write(mat_lu,*) 'After mat_report_sparsity: mem_allocated_global =', mem_allocated_global
   END SUBROUTINE mat_report_sparsity

!> \brief Returns sum of all elements of matrix.
!> \param A The input matrix 
!> \return The sum of all elements of matrix
    function mat_sum(A)
      implicit none
      type(Matrix), intent(in) :: A
      real(realk) :: mat_sum

      select case(matrix_type)
!       case(mtype_symm_dense)
!           mat_sum=mat_symm_dense_sum(A)
        case(mtype_dense)
            mat_sum=mat_dense_sum(A)
#ifndef UNITTEST
#endif
!       case(mtype_unres_symm_dense_cholseky)
!           mat_sum=mat_unres_symm_dense_sum(A)
       case(mtype_unres_dense)
           mat_sum=mat_unres_dense_sum(A)
        case default
            call lsquit("mat_sum not implemented for this type of matrix",-1)
      end select

      return
    end function mat_sum

#ifndef UNITTEST
!> \brief Inquire zero cutoff - for csr matrices only!! 
!> \author S. Host
!> \date 2009
!> \param cutoff The zero cutoff forcsr  matrices
    subroutine mat_inquire_cutoff(cutoff)
      implicit none
      real(realk), intent(out) :: cutoff

      select case(matrix_type)
!       case(mtype_symm_dense)
!           call mat_symm_dense_zero_cutoff(cutoff)
!       case(mtype_dense)
!           call mat_dense_zero_cutoff(cutoff)
!       case(mtype_unres_symm_dense_cholseky)
!           call mat_unres_symm_dense_zero_cutoff(cutoff)
!       case(mtype_unres_dense_zero_cutoff)
!           call mat_unres_dense_zero_cutoff(cutoff)
        case(mtype_csr)
             call mat_csr_inquire_cutoff(cutoff)
        case default
            call lsquit("mat_zero_cutoff not implemented for this type of matrix",-1)
      end select

      return
    end subroutine mat_inquire_cutoff
#endif

!> \brief Make c = alpha*diag(a)b + beta*c, where a is realk(:) b,c are type(matrix) and alpha,beta are parameters
!> \author B. Jansik
!> \date 2010
!> \param a The first realk(:) diagonal 
!> \param b The second type(matrix) factor
!> \param transb 'T'/'t' if b should be transposed, 'N'/'n' otherwise
!> \param alpha The alpha parameter
!> \param beta The beta parameter
!> \param c The output type(matrix)
    subroutine mat_dmul (a, b, transb, alpha, beta, c)
         implicit none
         real(realk), intent(in)  :: a(:)
         TYPE(Matrix), intent(IN) :: b
         character, intent(in)    :: transb
         REAL(realk), INTENT(IN)  :: alpha, beta
         TYPE(Matrix), intent(inout):: c
         integer :: ak, bk, ci, cj
         call time_mat_operations1
         select case(matrix_type)
         case(mtype_dense)
            call mat_dense_dmul(a,b,transb,alpha,beta,c)
         case(mtype_scalapack)
            call mat_scalapack_dmul(a,b,transb,alpha,beta,c)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            call mat_pdmm_dmul(a,b,transb,alpha,beta,c)
#endif
         case default
              call lsquit("mat_dmul not implemented for this type of matrix",-1)
         end select
         call time_mat_operations2(JOB_mat_dmul)
    end subroutine mat_dmul

    !> \brief Make Cij = alpha*Aij*Bij+beta*Cij (Hadamard product), where alpha is a realk and A,B,C are type(matrix)
    !> \author I.M Hoyvik
    !> \date 2012
    !> param alpha Scaling factor
    !> param beta Scaling factor
    !> param A type(matrix) input
    !> param B type(matrix) input
    !> param C type(matrix) output
    subroutine mat_hmul(alpha,A,B,beta,C)
    implicit none
    real(realk) :: alpha,beta
    type(matrix),intent(in) :: A,B
    type(matrix),intent(inout) :: C
    integer :: i
    call time_mat_operations1
 
    select case(matrix_type)
    case (mtype_dense)
       !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(A,B,C,alpha,beta)
       do i=1,A%nrow*A%ncol
          C%elms(i)= alpha*A%elms(i)*B%elms(i)+beta*C%elms(i)
       end do
       !$OMP END PARALLEL DO
    case (mtype_scalapack)
      call mat_scalapack_hmul(alpha,A,B,beta,C)
#ifdef VAR_ENABLE_TENSORS
    case (mtype_pdmm)
      call mat_pdmm_hmul(alpha,A,B,beta,C)
#endif
    case default
       call lsquit("mat_hmul not implemented for this type of matrix",-1)
    end select
    call time_mat_operations2(JOB_mat_hmul)
    end subroutine mat_hmul

    !> \brief Make Aij = Aij/(Bij -mu)
    !> \author B. Jansik
    !> \date 2012
    !> param alpha Scaling factor
    !> param A type(matrix) inout
    !> param B type(matrix) input
    !> param mu type(matrix) input
    subroutine mat_hdiv(A,B,mu)
    implicit none
    type(matrix),intent(inout) :: A
    type(matrix),intent(in) :: B
    real(realk), intent(in) :: mu

    real(realk) :: denom
    integer :: i
    call time_mat_operations1
   
    if ((A%nrow.ne.B%nrow).or.(A%ncol.ne.B%ncol)) &
    &  call lsquit('Incompatible matrix dimensions in mat_hdiv',-1)

    select case(matrix_type)
    case (mtype_dense)
       !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,denom) SHARED(A,B,mu)
       do i=1,A%nrow*A%ncol
          denom = B%elms(i) - mu
          if (denom.ne.0.0E0_realk) A%elms(i)= A%elms(i)/denom
       enddo
       !$OMP END PARALLEL DO
    case (mtype_scalapack)
      call mat_scalapack_hdiv(A,B,mu)
    case default
      call lsquit('mat_hdiv not implemented for this type of matrix',-1)
    end select
    call time_mat_operations2(JOB_mat_hdiv)
    end subroutine mat_hdiv


!> \brief   A = alpha*x*y^T+A
!> \author: Ida-Marie Hoeyvik
!> \date May 2012
!> \param alpha; realk scalar
!> \param x realk(:)
!> \param y realk(:)
!> \param A type(matrix)
  subroutine mat_dger(alpha,x,y,A)
    implicit none
    real(realk) :: alpha
    type(matrix) :: A
    real(realk) :: x(A%nrow),y(A%ncol)
    real(realk),pointer :: Afull(:,:)
    call time_mat_operations1

    select case(matrix_type)
    case(mtype_dense)
       call dger(A%nrow,A%ncol,alpha,x,1,y,1,A%elms,A%nrow)       
!    case(mtype_scalapack)
!       call mat_scalapack_dger(alpha,x,y,A)
#ifdef VAR_ENABLE_TENSORS
    case(mtype_pdmm)
       call mat_pdmm_dger(alpha,x,y,A)
#endif
    case default
       !FALLBACK 
       call mem_alloc(Afull,A%nrow,A%ncol)
       call mat_to_full(A,1.0E0_realk,Afull)
       call DGER(A%nrow,A%ncol,alpha,x,1,y,1,Afull,A%nrow)
       call mat_set_from_full(Afull,1.0E0_realk,A)
       call mem_dealloc(Afull)
!       call lsquit("mat_dger not implemented for this type of matrix",-1)
   end select
   call time_mat_operations2(JOB_mat_dger)
 end subroutine mat_dger


!> \brief Make C = alpha*diag(x)*A*B+beta*C (diagonal matrix multiplied with Hadamard product). Diagonal matrix is represented by vector a realk vector.
!> \author I.M Hoyvik
!> \date 2012
!> \param x realk(:)
!> \param A type(matrix) input
!> \param B type(matrix) input
!> \param transa 'T'/'t' if A should be transposed, 'N'/'n' otherwise
!> \param alpha realk scalar
!> \param beta realk scalar
!> \param C type(matrix) output
subroutine mat_dhmul(x,A,B,transa,transb,alpha,beta,C)
implicit none
type(matrix),intent(in) :: A,B
type(matrix),intent(inout) :: C
character, intent(in)    ::transa, transb
type(matrix) :: temp
real(realk)  :: x(A%nrow),alpha,beta

select case(matrix_type)
case default
  call mat_init(temp,A%nrow,A%ncol)
  call mat_dmul(x,A,transa,alpha,0d0,temp)
  call mat_hmul(1d0,temp,B,1d0,C)
  call mat_free(temp)
end select
end subroutine mat_dhmul





#ifndef UNITTEST
!> \brief Set zero cutoff - for CSR matrices only!! 
!> \author S. Host
!> \date 2009
!> \param zero The desired zero cutoff for CSR matrices
    subroutine mat_zero_cutoff(zero)
      implicit none
      real(realk), intent(in) :: zero

      select case(matrix_type)
!       case(mtype_symm_dense)
!           call mat_symm_dense_zero_cutoff(zero)
!       case(mtype_dense)
!           call mat_dense_zero_cutoff(zero)
!       case(mtype_unres_symm_dense_cholseky)
!           call mat_unres_symm_dense_zero_cutoff(zero)
!       case(mtype_unres_dense_zero_cutoff)
!           call mat_unres_dense_zero_cutoff(zero)
        case(mtype_csr)
             call mat_csr_zero_cutoff(zero)
        case default
            call lsquit("mat_zero_cutoff not implemented for this type of matrix",-1)
      end select

      return
    end subroutine mat_zero_cutoff
#endif

  END MODULE Matrix_operations_aux

