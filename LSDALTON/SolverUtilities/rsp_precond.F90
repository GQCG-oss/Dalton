!> @file
!> Contains module rspPrecond

!> \brief Contains routines for preconditioning of rsp equations in OAO basis.
!> \author S. Host
!> \date 2007
!>
!> See Coriani et al., JCP 126, 154108 (2007)
!>
module rspPrecond
use decompMod
use files
use precision
use matrix_module
use matrix_operations
use matrix_operations_aux,only:mat_zerohalf,mat_ao_precond

private
public :: rsp_Ared, rsp_Gred, rsp_Sred, rsp_boverlaps, lusigmarsp, lubrsp, lurhors,&
     & oao_rsp_solver, oao_rsp_setup_redspace, oao_rsp_solve_in_red_space, oao_red_space_EIGEN,&
     & oao_rsp_lintrans, oao_rsp_orthonormalize, oao_rsp_symm_orthonormalize, &
     & oao_rsp_test_for_symmetry, rsp_diag_precond
!For reduced space iterative solver. Should be put into structure
!instead of being global!
real(realk), allocatable, save :: rsp_Ared(:,:), rsp_Gred(:), rsp_Sred(:,:), rsp_boverlaps(:,:)
integer, save      :: lusigmarsp, lubrsp, lurhorsp

contains

   !> \brief Driver for solving rsp preconditioning equations in orthonormal AO basis.
   !> \author S. Host
   !> \date 2007
   subroutine oao_rsp_solver(decomp, res_in, omega, x)
      implicit none
      !> Contains matrices from OAO decomposition of overlap matrix
      TYPE(decompItem), intent(in)  :: decomp
      !> Residual to be preconditioned
      type(Matrix), intent(in)   :: res_in
      !> Solution vector, i.e. preconditioned residual (output)
      type(Matrix),intent(inout) :: x 
      !> Damping factor. Only output if debug_rsp_linsca is specified
      real(realk),intent(inout)  :: omega
      type(Matrix)             :: scrmat(1), b_current, sigma, rho, res
      integer                  :: i, j, k, matdim, max_it, iter, red_space_dim, IERR
      real(realk)              :: err, thresh, norm, test, fac
      real(realk)              :: t1, t2, t_lintra, t_ortho, t_setup, t_solve, t_write, & 
                                & t_precond, t_testforsym, t_init, t_end
      real(realk), allocatable :: testmat(:,:), fv1(:), fv2(:), eigenval(:), eigenvec(:,:)
      logical :: OnMaster
   !call CPU_TIME(t1)
   max_it = 20  !Stinne's first guess as to how many iterations are needed
                !FIXME: 30 is too many - for test calculations only. Reduce to, say, 15? 
   i = 1
   OnMaster =.TRUE.
   matdim = x%nrow

   allocate(rsp_Ared(2*max_it,2*max_it), rsp_Gred(2*max_it), rsp_Sred(2*max_it,2*max_it), rsp_boverlaps(2*max_it,2*max_it))
   call mat_init(b_current,matdim,matdim)
   call mat_init(scrmat(1),matdim,matdim)
   call mat_init(res,matdim,matdim)
   call mat_init(sigma,matdim,matdim)
   call mat_init(rho,matdim,matdim)
   rsp_Ared = 0.0E0_realk ; rsp_Gred = 0.0E0_realk ; rsp_Sred = 0.0E0_realk

   lusigmarsp = -1 ; lubrsp = -1 ; lurhorsp = -1
   CALL LSOPEN(lurhorsp,'rhovecs','unknown','UNFORMATTED')
   CALL LSOPEN(lusigmarsp,'sigmavecs','unknown','UNFORMATTED')
   CALL LSOPEN(lubrsp,'bvecs','unknown','UNFORMATTED')

   !if (decomp%debug_rsp_linsca) then
   !   if (decomp%info_rsp_precond) write (lupri,*) 'In oao_rsp_solver - HACK CODE TO INIGUESS'
   !   thresh = 1.0E-6_realk
   !   call oao_rsp_firstguess(1,scrmat(1),0) 
   !   !Generate starting guess???
   !        !write(lupri,*) 'firstguess:'
   !        !!call mat_print(scrmat(1),1,matdim,1,matdim,lupri)
   !else
      if (decomp%info_rsp_precond) write (decomp%lupri,*) 'In oao_rsp_solver - preconditioning of response equations'
      norm = sqrt(mat_sqnorm2(res_in))
      if (decomp%info_rsp_precond) write (decomp%lupri,*) 'oao_rsp_solver: norm of res in', norm
      thresh = norm*decomp%cfg_micro_thresh
      !call mat_copy(-1.0E0_realk,res_in,scrmat(1))
      call mat_assign(scrmat(1),res_in)
      call mat_scal(-1.0E0_realk, scrmat(1))
   !endif
   !call flshfo(lupri)

   !Test
   norm = sqrt(mat_sqnorm2(scrmat(1)))
   if (decomp%info_rsp_precond) write (decomp%lupri,*) 'oao_rsp_solver: norm of iniguess', norm
   !end test

   call oao_rsp_test_for_symmetry(decomp,scrmat(1))

   !call flshfo(decomp%lupri)
   !Test 
   norm = sqrt(mat_sqnorm2(scrmat(1)))
   if (decomp%info_rsp_precond) write (decomp%lupri,*) 'oao_rsp_solver: norm of iniguess after symcheck', norm
   !end test

   call project_oao_basis(decomp, scrmat(1), 0, b_current) !0 => neither symm nor antisymm
   !call flshfo(decomp%lupri)

   !Test 
   norm = sqrt(mat_sqnorm2(b_current))
   if (decomp%info_rsp_precond) write (decomp%lupri,*) 'oao_rsp_solver: norm of b_current', norm
   !end test
   !call flshfo(decomp%lupri)

   !Perform symmetric orthonormalization of (Z Y) and (Y Z) pair
   !and save first trial b vector (= incoming residual):
   call oao_rsp_symm_orthonormalize(decomp,b_current)
           !write(decomp%lupri,*) 'b_current after orthonormalize:'
           !call mat_print(b_current,1,matdim,1,matdim,decomp%lupri)
   
   call oao_rsp_lintrans(decomp,b_current,sigma,rho)
           !write(decomp%lupri,*) 'rho:'
           !call mat_print(rho,1,matdim,1,matdim,decomp%lupri)

   !Test 
   !norm = sqrt(mat_sqnorm2(rho))
   !if (decomp%info_rsp_precond) write (decomp%lupri,*) 'oao_rsp_solver: norm of rho', norm
   !norm = sqrt(mat_sqnorm2(sigma))
   !if (decomp%info_rsp_precond) write (decomp%lupri,*) 'oao_rsp_solver: norm of sigma', norm
   !norm = mat_dotproduct(b_current,rho)
   !if (decomp%info_rsp_precond) write (decomp%lupri,*) 'oao_rsp_solver: b_current dot rho', norm
   !end test
   !call flshfo(decomp%lupri)

   call mat_write_to_disk(lubrsp,b_current,OnMaster)
   call mat_write_to_disk(lusigmarsp,sigma,OnMaster)
   call mat_write_to_disk(lurhorsp,rho,OnMaster)

           !write(decomp%lupri,*) 'sigma_current:' 
           !call mat_print(sigma,1,matdim,1,matdim,decomp%lupri)
           !write(decomp%lupri,*) 'rho_current:'
           !call mat_print(rho,1,matdim,1,matdim,decomp%lupri)
           !write(decomp%lupri,*) 'b_current:'
           !call mat_print(b_current,1,matdim,1,matdim,decomp%lupri)

   !Setup reduced residual
   if (.not. decomp%debug_rsp_linsca) then
      rsp_Gred(1) = mat_dotproduct(b_current,res_in)
      call mat_trans(b_current, scrmat(1))
      rsp_Gred(2) = mat_dotproduct(scrmat(1),res_in)
      !write (decomp%lupri,*) 'gred 1, 2',  rsp_Gred(1), rsp_Gred(2)
   endif

   call oao_rsp_setup_redspace(decomp,b_current,sigma,rho,res_in,i)

   if (decomp%debug_rsp_linsca) then
      call oao_red_space_EIGEN(decomp,i,omega,x,res)
   else
      call oao_rsp_solve_in_red_space(decomp,i,res_in,omega,x,res)
   endif

        !write(decomp%LUPRI,*) 'res test'
        !call mat_print(res,1,matdim,1,matdim,decomp%lupri)
   !t_lintra = 0.0E0_realk ; t_ortho = 0.0E0_realk ; t_setup = 0.0E0_realk
   !t_solve = 0.0E0_realk  ; t_write = 0.0E0_realk ; t_precond = 0.0E0_realk ; t_testforsym = 0.0E0_realk
   !call CPU_TIME(t2)
   !t_init = t2 -t1

   do
      i = i + 1
      if (i == max_it + 1) then
         write (decomp%lupri,*) 'Warning: Max. no of iterations reached, exiting'
         exit
         !WRITE(decomp%LUPRI,'(/A)') &
         !&     'Linear equations not converged'
         !CALL QUIT(' Linear equations not converged')
      endif
      !Convergence?
      err = SQRT(mat_sqnorm2(res))
      !if (INFO_RSP) then
         write (decomp%LUPRI, '("          Preconditioning: Error of red space iteration", i3, ":  ", E16.10)') i-1, err
      !endif
      if (err < thresh) then
         !if (INFO_RSP) then
            write (decomp%lupri,'("     Preconditioning: Linear equations converged in", i3, " iterations!")') i-1
         !endif
         exit
      endif 
      !Preconditioning of residual = resP and b(iter+1) = resP:
      !call CPU_TIME(t1)
      call oao_rsp_test_for_symmetry(decomp,res)
      !call CPU_TIME(t2)
      !t_testforsym = t_testforsym + (t2 -t1)

      !call CPU_TIME(t1)
      call rsp_diag_precond(decomp,res,0,-omega,b_current)
      !call chol_AO_precond2(res,0,-omega,b_current) !Stinne: Change sign on omega 10/7-2009
                                                    !It should actually be done inside prec routine
                                                    !but that causes problems since it is also used for
                                                    !density optimization!
      !call CPU_TIME(t2)
      !t_precond = t_precond + (t2 -t1)
      
      !Orthonormalize b_current against set of previous b vectors and save:
      !call CPU_TIME(t1)
      call oao_rsp_orthonormalize(decomp,i-1,b_current)
      !call CPU_TIME(t2)
      !t_ortho = t_ortho + (t2 -t1)

      !if (linear_dep) exit
      !Save current on disc:
      !call CPU_TIME(t1)
      call mat_write_to_disk(lubrsp,b_current,OnMaster)
      !call CPU_TIME(t2)
      !t_write = t_write + (t2 -t1)

      !Calculate and save new sigma vector on disc:
      !call CPU_TIME(t1)
      call oao_rsp_lintrans(decomp, b_current,sigma,rho)
      !call CPU_TIME(t2)
      !t_lintra = t_lintra + (t2 -t1)

      !call CPU_TIME(t1)
      call mat_write_to_disk(lusigmarsp,sigma,OnMaster)
      call mat_write_to_disk(lurhorsp,rho,OnMaster)
      !call CPU_TIME(t2)
      !t_write = t_write + (t2 -t1)

      !call CPU_TIME(t1)
      call oao_rsp_setup_redspace(decomp,b_current, sigma, rho, res_in, i)
      !call CPU_TIME(t2)
      !t_setup = t_setup + (t2 -t1)

      if (decomp%debug_rsp_linsca) then
         call oao_red_space_EIGEN(decomp,i,omega,x,res)
      else
         !call CPU_TIME(t1)
         call oao_rsp_solve_in_red_space(decomp,i,res_in,omega,x,res)
         !call CPU_TIME(t2)
         !t_solve = t_solve + (t2 -t1)
      endif
   enddo

   !call CPU_TIME(t1)

   if (decomp%debug_rsp_linsca) then !Normalize X*S*X=1
      call oao_rsp_lintrans(decomp,x,scrmat(1),rho)
      fac = abs(mat_dotproduct(rho,x))
      write (decomp%lupri,*) 'X*S*X before normalization', fac
           !write(decomp%lupri,*) 'x before scaling:'
           !call mat_print(x,1,matdim,1,matdim,decomp%lupri)

      call mat_scal(1.0E0_realk/sqrt(fac),x)
           !write(decomp%lupri,*) 'x after scaling:'
           !call mat_print(x,1,matdim,1,matdim,decomp%lupri)

      !TEST
      call oao_rsp_lintrans(decomp,x,scrmat(1),rho)
      test = mat_dotproduct(rho,x) 
      write (decomp%lupri,*) 'X*S*X after normalization', test
      !END TEST
   endif

   CALL LSCLOSE(lusigmarsp,'DELETE')
   CALL LSCLOSE(lubrsp,'DELETE')
   CALL LSCLOSE(lurhorsp,'DELETE')

   deallocate(rsp_Ared, rsp_Gred, rsp_Sred, rsp_boverlaps)
   call mat_free(b_current)
   call mat_free(scrmat(1))
   call mat_free(res)
   call mat_free(sigma)
   call mat_free(rho)

   !call CPU_TIME(t2)
   !t_end = t2 -t1

   !write(decomp%lupri,*) 'Time used in init:', t_init
   !write(decomp%lupri,*) 'Time used in test for symmetry:', t_testforsym
   !write(decomp%lupri,*) 'Time used in diag precond + projection:', t_precond
   !write(decomp%lupri,*) 'Time used in orthonormalize:', t_ortho
   !write(decomp%lupri,*) 'Time used in writing to disk:', t_write
   !write(decomp%lupri,*) 'Time used in linear transformation:', t_lintra
   !write(decomp%lupri,*) 'Time used in setting up reduced space:', t_setup
   !write(decomp%lupri,*) 'Time used in solving in reduced space:', t_solve
   !write(decomp%lupri,*) 'Time used in shutdown:', t_end
   !t1 = t_init + t_testforsym + t_precond + t_ortho + t_write + t_lintra + t_setup + t_solve + t_end
   !write(decomp%lupri,*) 'Total time used in preconditioning:', t1

   end subroutine oao_rsp_solver

   !> \brief Set up reduced space.
   !> \author S. Host
   !> \date 2007
   subroutine oao_rsp_setup_redspace(decomp,b_current,sigma_current,rho_current,res_in,iter)
      implicit none

      !> Contains matrices from OAO decomposition of overlap matrix
      TYPE(decompItem), intent(in)  :: decomp
      !> Current trial vector
      type(Matrix), intent(in) :: b_current
      !> Sigma part of current linearly transformed trial vector
      type(Matrix), intent(in) :: sigma_current
      !> Rho part of current linearly transformed trial vector
      type(Matrix), intent(in) :: rho_current 
      !> The residual which is preconditioned in this module. Only used for testing element of reduced residual.
      type(Matrix), intent(in) :: res_in
      !> Iteration number => dim of red space = 2*iter (paired vectors)
      integer, intent(in)      :: iter
      type(Matrix)             :: vec_i, vec_i_trans, b_cur_trans, sigma_cur_trans, rho_cur_trans
      integer                  :: i, matdim
      real(realk)              :: error, test
      logical                  :: OnMaster
      OnMaster = .TRUE.
   matdim = b_current%nrow

   call mat_init(vec_i,matdim,matdim)
   call mat_init(vec_i_trans,matdim,matdim)
   call mat_init(b_cur_trans,matdim,matdim)
   call mat_init(sigma_cur_trans,matdim,matdim)
   call mat_init(rho_cur_trans,matdim,matdim)

   rewind(lusigmarsp)
   rewind(lurhorsp)
   rewind(lubrsp)

           !write(decomp%lupri,*) 'b_current:' 
           !call mat_print(b_current,1,matdim,1,matdim,decomp%lupri)
           !write(decomp%lupri,*) 'rho_current:'
           !call mat_print(rho_current,1,matdim,1,matdim,decomp%lupri)

   !Setup the 4 'diagonal' elements of reduced E2 and S2
   call mat_trans(b_current, b_cur_trans)
   call mat_trans(sigma_current, sigma_cur_trans)
   call mat_trans(rho_current, rho_cur_trans)

   rsp_Ared(2*iter-1,2*iter-1) =  mat_dotproduct(b_current,sigma_current)
   rsp_Ared(2*iter,2*iter-1)   =  mat_dotproduct(b_cur_trans,sigma_current)
   rsp_Ared(2*iter-1,2*iter)   =  mat_dotproduct(b_current,sigma_cur_trans)
   rsp_Ared(2*iter,2*iter)     =  mat_dotproduct(b_cur_trans,sigma_cur_trans)

   rsp_Sred(2*iter-1,2*iter-1) =  mat_dotproduct(b_current,rho_current)
   rsp_Sred(2*iter,2*iter-1)   =  mat_dotproduct(b_cur_trans,rho_current)
   rsp_Sred(2*iter-1,2*iter)   = -mat_dotproduct(b_current,rho_cur_trans)
   rsp_Sred(2*iter,2*iter)     = -mat_dotproduct(b_cur_trans,rho_cur_trans)

   rsp_boverlaps(2*iter-1,2*iter-1) =  mat_dotproduct(b_current,b_current)
   rsp_boverlaps(2*iter,2*iter-1)   =  mat_dotproduct(b_cur_trans,b_current)
   rsp_boverlaps(2*iter-1,2*iter)   =  mat_dotproduct(b_current,b_cur_trans)
   rsp_boverlaps(2*iter,2*iter)     =  mat_dotproduct(b_cur_trans,b_cur_trans)

   do i = 1, iter-1  !Setup lower half of reduced E2 and S2
      call mat_read_from_disk(lusigmarsp,vec_i,OnMaster)
      call mat_trans(vec_i, vec_i_trans)
      rsp_Ared(2*iter-1,2*i-1) = mat_dotproduct(b_current,vec_i)
      rsp_Ared(2*iter,2*i-1)   = mat_dotproduct(b_cur_trans,vec_i)
      rsp_Ared(2*iter-1,2*i)   = mat_dotproduct(b_current,vec_i_trans)
      rsp_Ared(2*iter,2*i)     = mat_dotproduct(b_cur_trans,vec_i_trans)

      call mat_read_from_disk(lurhorsp,vec_i,OnMaster)
      call mat_trans(vec_i, vec_i_trans)
      rsp_Sred(2*iter-1,2*i-1) =  mat_dotproduct(b_current,vec_i)
      rsp_Sred(2*iter,2*i-1)   =  mat_dotproduct(b_cur_trans,vec_i)
      rsp_Sred(2*iter-1,2*i)   = -mat_dotproduct(b_current,vec_i_trans)
      rsp_Sred(2*iter,2*i)     = -mat_dotproduct(b_cur_trans,vec_i_trans)
   enddo 

   do i = 1, iter-1  !Setup upper half of reduced E2 and S2
      call mat_read_from_disk(lubrsp,vec_i,OnMaster)
      call mat_trans(vec_i, vec_i_trans)
      rsp_Ared(2*i-1,2*iter-1) = mat_dotproduct(vec_i,sigma_current)
      rsp_Ared(2*i,2*iter-1)   = mat_dotproduct(vec_i_trans,sigma_current)
      rsp_Ared(2*i-1,2*iter)   = mat_dotproduct(vec_i,sigma_cur_trans)
      rsp_Ared(2*i,2*iter)     = mat_dotproduct(vec_i_trans,sigma_cur_trans)

      rsp_Sred(2*i-1,2*iter-1) =  mat_dotproduct(vec_i,rho_current)
      rsp_Sred(2*i,2*iter-1)   =  mat_dotproduct(vec_i_trans,rho_current)
      rsp_Sred(2*i-1,2*iter)   = -mat_dotproduct(vec_i,rho_cur_trans)
      rsp_Sred(2*i,2*iter)     = -mat_dotproduct(vec_i_trans,rho_cur_trans)

      rsp_boverlaps(2*i-1,2*iter-1) = mat_dotproduct(vec_i,b_current)
      rsp_boverlaps(2*iter-1,2*i-1) = mat_dotproduct(b_current, vec_i)

      rsp_boverlaps(2*i,2*iter-1)   = mat_dotproduct(vec_i_trans,b_current)
      rsp_boverlaps(2*iter-1,2*i)   = mat_dotproduct(b_current,vec_i_trans)

      rsp_boverlaps(2*i-1,2*iter)   = mat_dotproduct(vec_i,b_cur_trans)
      rsp_boverlaps(2*iter,2*i-1)   = mat_dotproduct(b_cur_trans,vec_i)

      rsp_boverlaps(2*i,2*iter)     = mat_dotproduct(vec_i_trans,b_cur_trans)
      rsp_boverlaps(2*iter,2*i)     = mat_dotproduct(vec_i_trans,b_cur_trans)
   enddo

   if (decomp%info_rsp_precond) then
      write (decomp%lupri,*) 'b vector overlaps, response solver:'
      call LS_OUTPUT(rsp_boverlaps, 1, 2*iter, 1, 2*iter, 200, 200, 1, decomp%lupri)

      write (decomp%lupri,*) 'Reduced E2, response solver:'
      call LS_OUTPUT(rsp_Ared, 1, 2*iter, 1, 2*iter, 200, 200, 1, decomp%lupri)

      write (decomp%lupri,*) 'Reduced S2, response solver:'
      call LS_OUTPUT(rsp_Sred, 1, 2*iter, 1, 2*iter, 200, 200, 1, decomp%lupri)
   endif

   error = abs(rsp_Ared(1,2*iter) - rsp_Ared(2*iter,1))
   if (abs(error) > 1.0E-10_realk) write (decomp%lupri,*) 'WARNING: Error on non-diag elms of reduced E2 ', error

   error = abs(rsp_Sred(1,2*iter) - rsp_Sred(2*iter,1))
   if (abs(error) > 1.0E-10_realk) write (decomp%lupri,*) 'WARNING: Error on non-diag elms of reduced S2 ', error

  !TEST: remaining elms of reduced residual should be zero:
   if (iter > 1 .and. .not. decomp%debug_rsp_linsca) then
      test = mat_dotproduct(b_current,res_in)
      !write (decomp%lupri,*) 'elm', 2*iter-1, 'of reduced residual:', test
      if (abs(test) > 1.0E-7_realk) then
         write (decomp%lupri,*) 'WARNING: Remaining elms of reduced residual should be zero'
         write (decomp%lupri,*) 'but element ', 2*iter-1, 'of reduced residual is ', test
         !STOP
      endif 
      call mat_trans(b_current, vec_i_trans)
      test = mat_dotproduct(vec_i_trans,res_in)
      !write (decomp%lupri,*) 'elm', 2*iter, 'of reduced residual:', test
      if (abs(test) > 1.0E-7_realk) then
         write (decomp%lupri,*) 'WARNING: Remaining elms of reduced residual should be zero'
         write (decomp%lupri,*) 'but element ', 2*iter, 'of reduced residual is ', test
         !STOP 'error in oao_rsp_setup_redspace'
      endif 
   endif
   !END TEST

   call mat_free(vec_i)
   call mat_free(vec_i_trans)
   call mat_free(b_cur_trans)
   call mat_free(sigma_cur_trans)
   call mat_free(rho_cur_trans)
   end subroutine oao_rsp_setup_redspace

   !> \brief Solve eigenvalue equation in already set up reduced space.
   !> \author S. Host
   !> \date 2007
   subroutine oao_rsp_solve_in_red_space(decomp,iter,res_in,omega,x,res)
      implicit none

      !> Contains matrices from OAO decomposition of overlap matrix
      TYPE(decompItem), intent(in)  :: decomp
      !> Iteration number => Dim of red space = 2*iter (paired vectors)
      integer, intent(in)      :: iter
      !> The residual which is preconditioned by this module
      type(Matrix), intent(in) :: res_in
      !> Frequency/excitation energy, used as shift in solution of equations
      real(realk), intent(in)  :: omega
      !> Full space solution vector constructed from red. space solution vector (output)
      type(Matrix),intent(inout) :: x
      !> Full space residual of linear equations 
      type(Matrix),intent(inout) :: res
      type(Matrix)             :: sigma_i, b_i, rho_i, scrmat, vec_i_trans, vec_i_trans2
      integer                  :: i, IERR, matdim
      integer, allocatable, dimension(:)      :: IPIV
      real(realk), allocatable, dimension(:)  :: RHS, xred
      real(realk),allocatable, dimension(:,:) :: A, S
      logical :: OnMaster
      IERR=0
      OnMaster = .TRUE.
   matdim = res_in%nrow

   allocate(A(2*iter,2*iter), S(2*iter,2*iter))
   allocate(RHS(2*iter), IPIV(2*iter))
   allocate(xred(2*iter))
   call mat_init(scrmat,matdim,matdim)
   call mat_init(b_i,matdim,matdim)
   call mat_init(sigma_i,matdim,matdim)
   call mat_init(rho_i,matdim,matdim)
   call mat_init(vec_i_trans,matdim,matdim)
   call mat_init(vec_i_trans2,matdim,matdim)

   !Setup reduced E2, S2, and right hand side with proper dimension
   A(1:2*iter,1:2*iter) = rsp_Ared(1:2*iter,1:2*iter)
   S(1:2*iter,1:2*iter) = rsp_Sred(1:2*iter,1:2*iter)
   RHS(1:2*iter) = rsp_Gred(1:2*iter)

   !write (decomp%lupri,*) "E2:"
   !call LS_OUTPUT(A, 1, 2*iter, 1, 2*iter, 2*iter, 2*iter, 1, decomp%lupri)
   !write (decomp%lupri,*) "S2:"
   !call LS_OUTPUT(S, 1, 2*iter, 1, 2*iter, 2*iter, 2*iter, 1, decomp%lupri)
   A = A + omega*S  !Stinne: SIGN CHANGED ON OMEGA!!! 10/7-2009 
                    !NB: comments not changed!
   !write (decomp%lupri,*) 'E2 - omega*S2:'
   !call LS_OUTPUT(A, 1, 2*iter, 1, 2*iter, 2*iter, 2*iter, 1, decomp%lupri)
  
   !write (decomp%lupri,*) 'RHS:'
   !call LS_OUTPUT(RHS, 1, 2*iter, 1, 1, 2*iter, 1, 1, decomp%lupri)
   !Solve set of linear equations Ax = b:
   call DGESV(2*iter, 1, A, 2*iter, IPIV, RHS, 2*iter, IERR) !Solution vector is found in RHS.
   if (IERR /= 0) then
      WRITE(decomp%LUPRI,'(/A, i4)') &
      &     'Problem in DGESV, IERR = ', IERR
      CALL lsQUIT(' Problem in DGESV',decomp%lupri)
   endif
   !write (decomp%lupri,*) 'Solution vector:'
   !call LS_OUTPUT(RHS, 1, 2*iter, 1, 1, 2*iter, 1, 1, decomp%lupri)
   xred = RHS
   !Now obtain x and residual in real space
   call mat_zero(res) ; call mat_zero(x)
   rewind(lusigmarsp) ; rewind(lubrsp) ; rewind(lurhorsp)
   do i = 1, iter
      call mat_read_from_disk(lusigmarsp,sigma_i,OnMaster)
      call mat_read_from_disk(lubrsp,b_i,OnMaster)
      call mat_read_from_disk(lurhorsp,rho_i,OnMaster)
   !x = xred(1)*b1 + xred(2)*b1_T  + xred(3)*b2 + xred*b2_T + ...
      call mat_daxpy(xred(2*i-1), b_i, x)
      call mat_trans(b_i, vec_i_trans)
      call mat_daxpy(xred(2*i), vec_i_trans, x)
   !res = E2*x - omega*S2*x - G 
   !    = xred(1)*(sigma1 - omega*rho1) + xred(2)*(sigma1_T + omega*rho1_T)
   !    + xred(3)*(sigma2 - omega*rho2) + xred(4)*(sigma2_T + omega*rho2_T)
   !    - res_in
      call MAT_ADD(1.0E0_realk, sigma_i, omega, rho_i, scrmat)  !Stinne: SIGN CHANGED ON OMEGA!!! 10/7-2009
      call mat_daxpy(xred(2*i-1), scrmat, res)            !NB: comments not changed!

      call mat_trans(sigma_i, vec_i_trans)
      call mat_trans(rho_i, vec_i_trans2)

      call MAT_ADD(1.0E0_realk, vec_i_trans, -omega, vec_i_trans2, scrmat)  !Stinne: SIGN CHANGED ON OMEGA!!! 10/7-2009
      call mat_daxpy(xred(2*i), scrmat, res)                          !NB: comments not changed!
   enddo
   call mat_daxpy(-1.0E0_realk, res_in, res) !Subtract right-hand-side from residual

   deallocate(A,RHS,IPIV,S,xred)
   call mat_free(scrmat)
   call mat_free(b_i)
   call mat_free(sigma_i)
   call mat_free(rho_i)
   call mat_free(vec_i_trans)
   call mat_free(vec_i_trans2)
   end subroutine oao_rsp_solve_in_red_space

   !> \brief Test routine - for debugging mini_solver for getting the initial guess
   !> \author S. Host
   !> \date 2007
   subroutine oao_red_space_EIGEN(decomp,iter,omega,x,res)
      implicit none

      !> Contains matrices from OAO decomposition of overlap matrix
      TYPE(decompItem), intent(in)  :: decomp
      !> Iteration number => Dim of red space = 2*iter (paired vectors)
      integer, intent(in)         :: iter
      !> Frequency/excitation energy, used as shift in solution of equations
      real(realk), intent(inout)  :: omega
      !> Full space solution vector constructed from red. space solution vector (output)
      type(Matrix),intent(inout) :: x
      !> Full space residual of linear equations 
      type(Matrix),intent(inout) :: res
      integer                  :: ISNDX(3), chosen_eival
      type(Matrix)             :: sigma_i, b_i, rho_i, scrmat, vec_i_trans, vec_i_trans2
      integer                  :: i, IERR, matdim
      real(realk), allocatable, dimension(:)  :: ALFR, ALFI, BETA, xred
      real(realk),allocatable, dimension(:,:) :: A, S, EIGENVEC
      logical :: eigenval_found,OnMaster
      OnMaster = .TRUE.
   eigenval_found = .false.

   matdim = x%nrow

   allocate(A(2*iter,2*iter), S(2*iter,2*iter), EIGENVEC(2*iter,2*iter))
   allocate(alfr(2*iter), alfi(2*iter), beta(2*iter))
   allocate(xred(2*iter))
   call mat_init(scrmat,matdim,matdim)
   call mat_init(b_i,matdim,matdim)
   call mat_init(sigma_i,matdim,matdim)
   call mat_init(rho_i,matdim,matdim)
   call mat_init(vec_i_trans,matdim,matdim)
   call mat_init(vec_i_trans2,matdim,matdim)

   !Setup reduced E2, S2, and right hand side with proper dimension
   A(1:2*iter,1:2*iter) = rsp_Ared(1:2*iter,1:2*iter)
   S(1:2*iter,1:2*iter) = rsp_Sred(1:2*iter,1:2*iter)

   !write (decomp%lupri,*) "E2:"
   !call LS_OUTPUT(A, 1, 2*iter, 1, 2*iter, 2*iter, 2*iter, 1, decomp%lupri)
   !write (decomp%lupri,*) "S2:"
   !call LS_OUTPUT(S, 1, 2*iter, 1, 2*iter, 2*iter, 2*iter, 1, decomp%lupri)

   !write (decomp%lupri,*) 'E2 - omega*S2:'
   !call LS_OUTPUT(A, 1, 2*iter, 1, 2*iter, 2*iter, 2*iter, 1, decomp%lupri)
  
   !Solve eigenvalue equation AX = omegaBX:
   call RGG(2*iter,2*iter,A,S,ALFR,ALFI,BETA,1,EIGENVEC,IERR)
   if (IERR /= 0) then
      WRITE(decomp%LUPRI,'(/A, i4)') &
      &     'Problem in RGG, IERR = ', IERR
      CALL lsQUIT(' Problem in RGG',decomp%lupri)
   endif
   !write (decomp%lupri,*) 'Numerator: Real part of eigenval:'
   !call LS_OUTPUT(ALFR, 1, 2*iter, 1, 1, 2*iter, 1, 1, decomp%lupri)
   !write (decomp%lupri,*) 'Numerator: Im part of eigenval:'
   !call LS_OUTPUT(ALFI, 1, 2*iter, 1, 1, 2*iter, 1, 1, decomp%lupri)
   !write (decomp%lupri,*) 'Denominator: Real part of eigenval:'
   !call LS_OUTPUT(BETA, 1, 2*iter, 1, 1, 2*iter, 1, 1, decomp%lupri)
   !write (decomp%lupri,*) 'Eigenvectors:'
   !call LS_OUTPUT(EIGENVEC, 1, 2*iter, 1, 2*iter, 2*iter, 2*iter, 1, decomp%lupri)
  
   do i = 1, 2*iter
      if (abs(ALFI(i)) > 1.0E-10_realk) then
         STOP 'Imaginary eigenvalue in oao_rsp_solve_in_red_space_EIGEN'
      endif
      ALFR(i) = ALFR(i)/BETA(i)
   enddo 

   !write (decomp%lupri,*) 'Eigenvalues:'
   !call LS_OUTPUT(ALFR, 1, 2*iter, 1, 1, 2*iter, 1, 1, decomp%lupri)

   !Find lowest positive eigenvalue:
   omega = 1.0E6_realk
   do i = 1, 2*iter
      if (ALFR(i) > 0.0E0_realk .and. ALFR(i) < omega) then
         omega = ALFR(i)
         chosen_eival = i
         eigenval_found = .true.
      endif
   enddo 
   !Find highest negative eigenvalue
   !omega = -1.0E6_realk
   !do i = 1, 2*iter
   !   if (ALFR(i) < 0.0E0_realk .and. ALFR(i) > omega) then
   !      omega = ALFR(i)
   !      chosen_eival = i
   !      eigenval_found = .true.
   !   endif
   !enddo 
   if (.not. eigenval_found) STOP 'Error in oao_red_space_EIGEN'

   write (decomp%lupri,*) 'Chosen eigenval: ', omega
   xred =  EIGENVEC(1:2*iter,chosen_eival)
   
   !write (decomp%lupri,*) 'Chosen eigenvec:'
   !call LS_OUTPUT(xred, 1, 2*iter, 1, 1, 2*iter, 1, 1, decomp%lupri)
   
   !Now obtain x and residual in real space
   call mat_zero(res) ; call mat_zero(x)
   rewind(lusigmarsp) ; rewind(lubrsp) ; rewind(lurhorsp)
   do i = 1, iter
      call mat_read_from_disk(lusigmarsp,sigma_i,OnMaster)
      call mat_read_from_disk(lubrsp,b_i,OnMaster)
      call mat_read_from_disk(lurhorsp,rho_i,OnMaster)
   !x = xred(1)*b1 + xred(2)*b1_T  + xred(3)*b2 + xred*b2_T + ...
      call mat_daxpy(xred(2*i-1), b_i, x)
      call mat_trans(b_i, vec_i_trans)
      call mat_daxpy(xred(2*i), vec_i_trans, x)
        !write(decomp%LUPRI,*) 'x in real space:'
        !call mat_print(x,1,matdim,1,matdim,decomp%lupri)

   !res = E2*x - omega*S2*x 
   !    = xred(1)*(sigma1 - omega*rho1) + xred(2)*(sigma1_T + omega*rho1_T)
   !    + xred(3)*(sigma2 - omega*rho2) + xred(4)*(sigma2_T + omega*rho2_T)
   !    - res_in
      call MAT_ADD(1.0E0_realk, sigma_i, -omega, rho_i, scrmat) 
      call mat_daxpy(xred(2*i-1), scrmat, res)

      call mat_trans(sigma_i, vec_i_trans)
      call mat_trans(rho_i, vec_i_trans2)

      call MAT_ADD(1.0E0_realk, vec_i_trans, omega, vec_i_trans2, scrmat) 
      call mat_daxpy(xred(2*i), scrmat, res)
   enddo

   deallocate(A,S,EIGENVEC,ALFR,ALFI,BETA,xred)
   call mat_free(scrmat)
   call mat_free(b_i)
   call mat_free(sigma_i)
   call mat_free(rho_i)
   call mat_free(vec_i_trans)
   call mat_free(vec_i_trans2)
   end subroutine oao_red_space_EIGEN

   !> \brief Linear transformation.
   !> \author S. Host
   !> \date 2007
   subroutine oao_rsp_lintrans(decomp,x_in,sigma,rho)
   implicit none  
         !> Contains matrices from OAO decomposition of overlap matrix
         TYPE(decompItem), intent(in) :: decomp
         !> Trial vector to be linearly transformed
         TYPE(matrix), intent(in)     :: x_in
         !> Sigma part of linear transformation, (FUQ-FUP)*X + X*(FUQ-FUP) (output)
         TYPE(matrix), intent(inout)  :: sigma
         !> Rho part of linear transformation, -omega*(DX - XD) (output)
         TYPE(matrix), intent(inout)  :: rho
         TYPE(matrix)              :: FX, XF, FQminusFP, DX, XD

     call MAT_INIT(FQminusFP,x_in%nrow,x_in%ncol)
     call MAT_INIT(FX,x_in%nrow,x_in%ncol)
     call MAT_INIT(XF,x_in%nrow,x_in%ncol)
     call MAT_INIT(DX,x_in%nrow,x_in%ncol)
     call MAT_INIT(XD,x_in%nrow,x_in%ncol)

     call MAT_ADD(1.0E0_realk,decomp%FUQ,-1.0E0_realk,decomp%FUP,FQminusFP)
     !Sigma part of linear transformation, (FUQ-FUP)*X + X*(FUQ-FUP):
     call MAT_MUL(FQminusFP,x_in,'n','n',1.0E0_realk,0.0E0_realk,FX)
     call MAT_MUL(x_in,FQminusFP,'n','n',1.0E0_realk,0.0E0_realk,XF)
     call MAT_ADD(1.0E0_realk,FX,1.0E0_realk,XF,sigma)
     !Rho part, -omega*(DX - XD)
     call MAT_MUL(decomp%DU,x_in,'n','n',1.0E0_realk,0.0E0_realk,DX)
     call MAT_MUL(x_in,decomp%DU,'n','n',1.0E0_realk,0.0E0_realk,XD)
  !NOTE TO STINNE AND POUL FROM STINNE: CHANGE OF SIGN ON THIS LINEAR
  !TRANSFORMATION DESTROYS THE OVERALL RESPONSE CONVERGENCE (WILL NOT
  !BE THE SAME AS THE CONVERGENCE WITH MO PRECONDITIONING)!!!!! 16/8 - 2006
     call MAT_ADD(1.0E0_realk,DX,-1.0E0_realk,XD,rho) !Sign changed 10/8-06 Stinne !Changed back 16/8-06!!!!
                                          
     call MAT_FREE(FX)
     call MAT_FREE(XF)
     call MAT_FREE(FQminusFP)
     call MAT_FREE(DX)
     call MAT_FREE(XD)
   end subroutine oao_rsp_lintrans

  !> \brief Orthonormalization routine.
  !> \author S. Host
  !> \date 2007
  !> 
  !>  Modified 'orthonormalize' from module RESPONSE_SOLVER. \
  !>  Purpose: \n
  !>   Orthogonalize new b-vector against all previous b-vectors
  !>   and among themselves, and renormalize.
  !>   The b-vectors have the form \n
  !>         ( Y_dia    Z   )   ,    Z_mu_nu  mu < nu ; Y_dia = Y_mu_mu \n
  !>         (  Y     Y_dia )   ,    Y_mu_nu  mu > nu                   \n
  !>   Each b-vector is in (Z, Y) form, the (Y, Z) vector is obtained
  !>   by transposing.
  !>   (Orthogonalization is performed twice if round-off is large,
  !>    if larger than THRRND). \n
  !> 
  subroutine oao_rsp_orthonormalize(decomp,iter,b_current)
    implicit none
    !> Contains matrices from OAO decomposition of overlap matrix
    TYPE(decompItem), intent(in)  :: decomp
    !> Number of previous trial vectors (stored on file lubrsp)
    integer, intent(in)     :: iter
    !> The vector to be orthogonalized against set of previous trial vectors
    type(Matrix),intent(inout) :: b_current
    integer :: i, matdim
    type(matrix) :: B_scr, b_i
    real(realk) :: TT,T1,T2
    logical :: OnMaster
    OnMaster = .TRUE.
    matdim = b_current%nrow

    call mat_init(B_scr,matdim,matdim)
    call mat_init(b_i,matdim,matdim)

    rewind(lubrsp)
    do i = 1,iter
          call mat_read_from_disk(lubrsp,b_i,OnMaster)
          TT = mat_dotproduct(b_i,b_current)
          call mat_daxpy(-TT,b_i,b_current)
          call mat_trans(b_i,b_scr)
          TT = mat_dotproduct(B_scr,b_current)
          call mat_daxpy(-TT,B_scr,b_current)
    enddo

    call oao_rsp_symm_orthonormalize(decomp,b_current)

    call mat_free(B_scr)
    call mat_free(b_i)
  end subroutine oao_rsp_orthonormalize

  !> \brief Symmmetric orthonormalization routine.
  !> \author S. Host
  !> \date 2007
  !>
  !>  Perform symmetric orthonormalization of (Z Y) and (Y Z) pair
  !>  for vector bvec \n
  !>    -1/2       ( C1   C2 )               (  1     OVLPI ) \n
  !>   S      =    (         )   where S  =  (              ) \n
  !>               ( C2   C1 )               ( OVLPI     1  ) \n
  !>  - i.e. bvec should be orthonormal to its own transpose
  !> 
  subroutine oao_rsp_symm_orthonormalize(decomp,bvec)
  use matrix_util
    implicit none
    !> Contains matrices from OAO decomposition of overlap matrix
    TYPE(decompItem), intent(in)  :: decomp
    !> The trial bvector to be symmetrically orthonormalized
    type(Matrix), intent(inout) :: bvec
    type(Matrix) :: b_scr,b_scr2
    real(realk) :: ovlpi,x1,x2,c1,c2
    integer :: i,ndim

    ndim = bvec%nrow
    call mat_init(b_scr,ndim,ndim)
    call mat_init(b_scr2,ndim,ndim)

    call normalize(bvec)

!do it twice for numerical stability
    do i = 1,2
      call mat_trans(bvec,b_scr)
      ovlpi = mat_dotproduct(b_scr,bvec)!???what should be put here??? see rsp_solver.f90 l. 697)
      x1 = mat_sqnorm2(bvec)
      if (decomp%info_rsp_precond) then
         WRITE (decomp%LUPRI,'(A,1P,D14.6,D19.10)') &
         &         ' oao_rsp_symm_orth, norm squared and <ZY/YZ> overlap before S-1/2',X1,OVLPI
      endif
      X1 = 1E0_realk+OVLPI
      X2 = 1E0_realk-OVLPI
      if (decomp%info_rsp_precond) write (decomp%lupri,*) 'X1, X2: ', X1, X2
      IF (ABS(X1 - X2) > 1E-20_realk) THEN
         X1 = 0.5E0_realk / SQRT(X1)
         X2 = 0.5E0_realk / SQRT(X2)
         C1 = X1 + X2
         C2 = X1 - X2
         call mat_assign(b_scr2,bvec)
         call mat_add(c1,b_scr2,c2,b_scr,bvec)
         call normalize(bvec)
      endif
    enddo   
    call mat_free(b_Scr)
    call mat_free(b_scr2)
   end subroutine oao_rsp_symm_orthonormalize

  !> \brief Test trial vectors for symmetry.
  !> \author S. Host
  !> \date 2007
  !>
  !> The trial vectors can't be symmetric or antisymmetric,
  !> this causes linear dependencies. Here, we test for this and 
  !> correct if necessary (change a (z z) or (z -z) vector to a (0 z) vector). 
  !>
  subroutine oao_rsp_test_for_symmetry(decomp,bvec)
    implicit none
    !> Contains matrices from OAO decomposition of overlap matrix
    type(decompItem),intent(in) :: decomp
    !> The trial bvector to be tested for symmetry
    type(Matrix), intent(inout) :: bvec
    type(Matrix) :: b_trans,test
    real(realk) :: norm1, norm2
    integer :: i,ndim

    ndim = bvec%nrow
    call mat_init(b_trans,ndim,ndim)
    call mat_init(test,ndim,ndim)

    call mat_trans(bvec,b_trans)
    !Test for symmetry:
    call mat_add(1.0E0_realk, bvec, -1.0E0_realk, b_trans, test)
    norm1 = sqrt(mat_sqnorm2(test))
    !Test for antisymmetry:
    call mat_add(1.0E0_realk, bvec, 1.0E0_realk, b_trans, test)
    norm2 = sqrt(mat_sqnorm2(test))
    if (norm1 < 1.0E-4_realk) then
       if (decomp%info_rsp_precond) write (decomp%lupri,*) 'Oops. Trial vector is symmetric. I take action!'
       call mat_zerohalf('ut',bvec)
       call mat_scal_dia(0.5E0_realk,bvec)
       call project_oao_basis(decomp, bvec, 0, test)
       call mat_assign(bvec,test)
    else if (norm2 < 1.0E-4_realk) then
       if(decomp%info_rsp_precond) write (decomp%lupri,*) 'Oops. Trial vector is antisymmetric. I take action!'
       call mat_zerohalf('lt',bvec)
       call project_oao_basis(decomp, bvec, 0, test)
       call mat_assign(bvec,test)
    else if (decomp%info_rsp_precond) then
       write (decomp%lupri,*) 'Symmetry test OK - trial vector neither symmetric nor antisymmetric' 
    endif
    call mat_free(b_trans)
    call mat_free(test)
   end subroutine oao_rsp_test_for_symmetry

  !> \brief Wrapper routine for preconditioning of linear equations A Gnt = Gn
  !> \author S. Host
  !> \date 2007
   subroutine rsp_diag_precond(decomp,Gn,symm,omega,Gnt)
      !Preconditioning of Gn by solving A Gnt = Gn
      implicit none
      !> Contains matrices from OAO decomposition of overlap matrix
      type(decompItem),intent(in) :: decomp
      !> Residual to be preconditioned
      type(Matrix), intent(in) :: Gn
      !> Symmetry of residual, 1 = symmetric, 2 = antisymmetric
      integer, intent(in) :: symm
      !> Level shift
      real(realk), intent(in) :: omega
      !> Preconditioned residual (output)
      type(Matrix), intent(inout) :: Gnt
      type(Matrix) :: wrk, id
      real(realk)  :: norm
      integer      :: ndim

   ndim = Gn%nrow

   if (decomp%cfg_NOPREC) then
      call mat_assign(Gnt,Gn)
   else
      call mat_init(wrk,ndim,ndim)

      call mat_assign(wrk,Gn)
      call mat_ao_precond(symm,omega,decomp%FUP,decomp%FUQ,decomp%DU,wrk)
      call project_oao_basis(decomp, wrk, symm, Gnt)

      call mat_free(wrk)
   endif

    end subroutine rsp_diag_precond

end module rspPrecond
