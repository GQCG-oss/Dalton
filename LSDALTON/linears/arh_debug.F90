!> @file 
!> Contains debug module for augmented Roothaan-Hall code.

!> Debug routines for ARH module. Mainly different ways to compute Hessian and eigenvalues.
!> \author S. Host
!> \date 2005
module arh_debugging
use precision
use matrix_module
use matrix_operations
use matrix_operations_aux
use matrix_util
use memory_handling
use matrix_op_unres_dense, only: mat_unres_dense_create_elm_alph,&
     & mat_unres_dense_create_elm_beta
!> True if the debug header should be printed. Set to false after 1st iteration.
logical, save :: print_header = .true.

contains

!> \brief Compare numerical and analytical ARH Hessians and gradients.
!> \author S. Host
!> \date 2005
!> \param fifoqueue Contains densities and Fock/KS matrices from previous iterations
!> \param gradmat Analytical gradient (e.g. obtained from linear transformation). 
!> \param hes1 Analytical hessian (e.g. obtained from linear transformation).
!>
!> Very expensive, so don't use this for more than about 20 basisfunctions or so, at most.
!>
!   subroutine debug_hes(fifoqueue, gradmat,hes1)
!   implicit none
!           real(realk)             :: test1, test2
!           TYPE(matrix), intent(in):: gradmat, hes1 
!           TYPE(matrix)            :: grad1, grad2, gradtest, hes2, hestest, gradvec
!           type(modFIFO)           :: fifoqueue
!           integer                 :: matdim, vecdim, graddim
!
!   call MAT_INIT(grad1,vecdim,1)
!   call MAT_INIT(grad2,vecdim,1)
!   call MAT_INIT(gradtest,vecdim,1)
!   call MAT_INIT(hes2,vecdim,vecdim)
!   call MAT_INIT(hestest,vecdim,vecdim)
!
!   call MAT_TO_VEC('a', gradmat, grad1)
!   !write (LUPRI,*) "Before finite diff" ; call flshfo(lupri)
!   call finite_diff(fifoqueue, grad2, hes2)
!   !write (LUPRI,*) "After finite diff" ; call flshfo(lupri)
!   call MAT_ADD(1.0E0_realk, grad1, 1.0E0_realk, grad2, gradtest)  !The gradients have opposite signs
!      !write (LUPRI,*) "Gradient, linear transformation:"
!      !call MAT_PRINT(grad1, 1, grad1%nrow, 1, grad1%ncol, LUPRI)
!      !write (LUPRI,*) "Gradient, finite difference:"
!      !call MAT_PRINT(grad2, 1, grad2%nrow, 1, grad2%ncol, LUPRI)
!
!   test1 = sqrt(MAT_SQNORM2(gradtest))
!   write (LUPRI,*) "Diff. between lintrans and finite difference gradients:", test1 ; call flshfo(lupri)
!
!   call MAT_ADD(1.0E0_realk, hes1, -1.0E0_realk, hes2, hestest)
!      !write (LUPRI,*) "Hessian, linear transformation:"
!      !call MAT_PRINT(hes1, 1, hestest%nrow, 1, hestest%ncol, LUPRI)
!      !write (LUPRI,*) "Hessian, finite difference:"
!      !call MAT_PRINT(hes2, 1, hestest%nrow, 1, hestest%ncol, LUPRI)
!
!   test2 = sqrt(MAT_SQNORM2(hestest))
!   write (LUPRI,*) "Diff. between lintrans and finite difference Hessians:", test2
!
!   !Diagonalize Hessian:
!   write(lupri,*) 'Diagonalizing ARH Hessian:'
!   call util_diag(hes2,.false.,0.25E0_realk,'ARH lowest hes eigenval: ')
!
!   call MAT_FREE(grad1)
!   call MAT_FREE(grad2)
!   call MAT_FREE(gradtest)
!   call MAT_FREE(hes2)
!   call MAT_FREE(hestest)
!   end subroutine debug_hes

!> \brief Set up and diagonalize full exact second order Hessian.
!> \author S. Host
!> \date 2005
!>
!> Very expensive, so don't use this for more than about 20 basisfunctions or so, at most.
!>
   subroutine debug_diag_full_hessian(solver,decomp)
   use arhDensity
   use decompMod
   implicit none
           !> Contains info and setting for linear equations solver
           type(solverItem),intent(inout) :: solver
           !> Contains decomposed overlap matrix (Löwdin, Cholesky or other) 
           type(decompItem),intent(in)    :: decomp
           integer                 :: ndim, m, l, hesdim
           TYPE(matrix)            :: hessian, x_trial, x_trial_mat, column, scr
           logical                 :: temp

      ndim = decomp%U%nrow

      if (decomp%cfg_unres) then
         call debug_diag_full_hessian_unres(solver,decomp)
      else
         hesdim = ndim*(ndim+1)/2 - ndim
         solver%set_do_2nd_order = .true.

         call mat_init(hessian,hesdim,hesdim)
         call MAT_INIT(scr,ndim,ndim)
         call MAT_INIT(x_trial,hesdim,1)
         call MAT_INIT(x_trial_mat,ndim,ndim)
         call MAT_INIT(column,hesdim,1)

         do m = 1, hesdim
            call MAT_ZERO(x_trial)
            call MAT_ZERO(x_trial_mat)
            call mat_create_elm(m, 1, 1.0E0_realk, x_trial)
            call mat_VEC_TO_MAT('a', x_trial, scr)
            !write (LUPRI,*) "xmat:"
            !call MAT_PRINT(scr, 1, scr%nrow, 1, scr%ncol, LUPRI)
            call project_oao_basis(decomp, scr, arh_antisymmetric, x_trial_mat)
            !write (LUPRI,*) "xmat projected:"
            !call MAT_PRINT(x_trial_mat, 1, x_trial_mat%nrow, 1, x_trial_mat%ncol, LUPRI)
            !column = m'th column of hessian:
            call arh_lintrans(solver,decomp,x_trial_mat,arh_antisymmetric,0.0E0_realk,scr)
            call MAT_TO_VEC('a', scr, column)
            do l = 1, hesdim !Put elements of column in m'th column of hessian
               call mat_create_elm(l,m,column%elms(l),hessian)
            enddo         
         enddo

         !call mat_scal(4.0E0_realk, hes)

         !write (decomp%LUPRI,*) "Hessian:"
         !call MAT_PRINT(hessian, 1, hessian%nrow, 1, hessian%ncol, decomp%LUPRI)

         !Diagonalize Hessian:
         write(decomp%lupri,*) 'Diagonalizing full, exact Hessian:'
         call util_diag(decomp%lupri,hessian,.false.,1.0E0_realk,'Lowest Hessian eigenvalue: ')

         solver%set_do_2nd_order = solver%cfg_do_2nd_order

         call MAT_FREE(hessian)
         call MAT_FREE(x_trial)
         call MAT_FREE(x_trial_mat)
         call MAT_FREE(column)
         call MAT_FREE(scr)
      endif
   end subroutine debug_diag_full_hessian

!> \brief Set up and diagonalize full exact second order Hessian (for unrestricted).
!> \author S. Host
!> \date 2005
!>
!> Very expensive, so don't use this for more than about 20 basisfunctions or so, at most.
!>
   subroutine debug_diag_full_hessian_unres(solver, decomp)
   use arhDensity
   use decompMod
   implicit none
           !> Contains info and setting for linear equations solver
           type(solverItem),intent(inout) :: solver
           !> Contains decomposed overlap matrix (Löwdin, Cholesky or other) 
           type(decompItem),intent(in)    :: decomp
           integer                 :: ndim, m, l, hesdim, k
           TYPE(matrix)            :: hessian, x_trial, x_trial_mat, column, scr
           logical                 :: temp, alpha
           real(realk),pointer :: fullhes(:,:) , parthes(:,:)
           interface
              subroutine dsyevx_interface(A,ndim,print_eivecs,lupri,string)
                use precision
                implicit none
                character(6), intent(in), optional :: string
                integer, intent(in)      :: ndim, lupri
                real(realk), intent(in)  :: A(ndim,ndim)
                logical, intent(in)      :: print_eivecs 
              end subroutine dsyevx_interface
           end interface

      ndim = decomp%U%nrow

      hesdim = ndim*(ndim+1)/2 - ndim
      solver%set_do_2nd_order = .true.

      alpha = .true.

      !call mat_init(hessian,hesdim,hesdim)
      call MAT_INIT(scr,ndim,ndim)
      call MAT_INIT(x_trial,hesdim,1)
      call MAT_INIT(x_trial_mat,ndim,ndim)
      call MAT_INIT(column,hesdim,1)
      call mem_alloc(fullhes,2*hesdim,2*hesdim)

      do m = 1, 2*hesdim
         call MAT_ZERO(x_trial)
         call MAT_ZERO(x_trial_mat)
         if (m <= hesdim) then
            call mat_unres_dense_create_elm_alph(m, 1, 1.0E0_realk, x_trial)
         else
            call mat_unres_dense_create_elm_beta(m-hesdim, 1, 1.0E0_realk, x_trial)
         endif
         !write (decomp%LUPRI,*) "xvec:"
         !call MAT_PRINT(x_trial, 1, x_trial%nrow, 1, x_trial%ncol, decomp%LUPRI)
         call mat_VEC_TO_MAT('a', x_trial, scr)
         !write (decomp%LUPRI,*) "xmat:"
         !call MAT_PRINT(scr, 1, scr%nrow, 1, scr%ncol, decomp%LUPRI)
         call project_oao_basis(decomp, scr, arh_antisymmetric, x_trial_mat)
         !write (decomp%LUPRI,*) "xmat projected:"
         !call MAT_PRINT(x_trial_mat, 1, x_trial_mat%nrow, 1, x_trial_mat%ncol, decomp%LUPRI)
         !column = m'th column of hessian:
         call arh_lintrans(solver,decomp,x_trial_mat,arh_antisymmetric,0.0E0_realk,scr)
         !write (decomp%LUPRI,*) "Linear transformation:"
         !call MAT_PRINT(scr, 1, scr%nrow, 1, scr%ncol, decomp%LUPRI)
         call MAT_TO_VEC('a', scr, column)
         !write (decomp%LUPRI,*) "Column:"
         !call MAT_PRINT(column, 1, column%nrow, 1, column%ncol, decomp%LUPRI)
         do l = 1, hesdim !Put elements of column in m'th column of hessian
            fullhes(m,l)        = column%elms(l)
            fullhes(m,l+hesdim) = column%elmsb(l)
         enddo         
      enddo

      !write (decomp%LUPRI,*) "Unrestricted Hessian:"
      !call LS_OUTPUT(fullhes, 1, 2*hesdim, 1, 2*hesdim, 2*hesdim, 2*hesdim, 1, decomp%lupri)

      call dsyevx_interface(fullhes,2*hesdim,.false.,decomp%lupri,'HesA+B')

      solver%set_do_2nd_order = solver%cfg_do_2nd_order

      call MAT_FREE(x_trial)
      call MAT_FREE(x_trial_mat)
      call MAT_FREE(column)
      call MAT_FREE(scr)
      call mem_dealloc(fullhes)
   end subroutine debug_diag_full_hessian_unres

!Commented out because it it not yet changed to fit new configuration structure. /Stinne
!!> \brief Debug ARH Hessian by using energy function and comparing to linear transformation.
!!> \author S. Host
!!> \date 2005
!!> \param lintrans Linear transformation of x
!!> \param x Trial vector for debugging linear transf. against energy function
!!> \param fifoqueue Contains densities and Fock/KS matrices from previous iterations
!!>
!!> Very expensive, so don't use this for more than about 20 basisfunctions or so, at most.
!!>
!   subroutine debug_hes_via_energy(lintrans, x, fifoqueue)
!   !Check linear transformation by 
!   !(E[2]*x)_i = {E(beta*x + alfa*delta_i) - E(beta*x - alfa*delta_i)+ E(-beta*x - alfa*delta_i)- E(-beta*x + alfa*delta_i)}
!   !              ___________________________________________________________________________________________________________
!   !                                               4*alfa*beta
!   !delta_i is an antisymmetric matrix corresponding to a unit vector with 1 one the i'th position
!   !DIRECT MEDDLING WITH MATRIX ELEMENTS - USE ONLY FOR DENSE MATRICES
!   implicit none
!
!           TYPE(matrix), intent(in):: lintrans !E[2]*x obtained from linear transformation
!           TYPE(matrix), intent(in):: x        !The x used in linear transformation
!           !integer, intent(in)     :: pos      !position in queue to be used
!           real(realk)             :: alfa, beta, norm, e1, e2, e3, e4, lintra_elm, diff
!           TYPE(matrix)            :: delta_i, testmat
!           !TYPE(util_HistoryStore), intent(in) :: queue
!           type(modFIFO), intent(in)           :: fifoqueue
!           integer                             :: ndim, i, j, nelms, counter, do_check
!
!   ndim = x%nrow
!   nelms = ndim*ndim
!   do_check = nelms/10
!   
!   call MAT_INIT(delta_i,ndim,ndim)
!   call MAT_INIT(testmat,ndim,ndim)
!
!   alfa = 1.0E-4_realk
!   norm = sqrt(mat_sqnorm2(x))
!   !write(lupri,*) 'norm in debug_hes_via_energy', norm
!   beta = 1.0E-4_realk/norm
!   !write(lupri,*) 'alfa, beta in debug_hes_via_energy', alfa, beta
!
!   counter = 0
!   do i = 1, ndim
!      do j = 1, i
!         counter = counter + 1
!         !Random check:
!         if (counter == do_check) then
!            counter = 0
!            !Create delta_i*alfa
!            call mat_zero(delta_i)
!            call mat_create_elm(i, j, 1.0E0_realk, delta_i)
!            call mat_create_elm(j, i, -1.0E0_realk, delta_i)
!            if (i == j) call mat_zero(delta_i)
!            !write (LUPRI,*) "delta_i in debug_hes_via_energy:"
!            !call MAT_PRINT(delta_i, 1, ndim, 1, ndim, LUPRI)
!            call mat_scal(alfa, delta_i)
!
!            !Calculate first energy contribution
!            call MAT_ADD(beta, x, 1.0E0_realk, delta_i, testmat)
!            !write (LUPRI,*) "testmat in debug_hes_via_energy:"
!            !call MAT_PRINT(testmat, 1, ndim, 1, ndim, LUPRI)
!            call e_appr(testmat, fifoqueue, e1)
!
!            !Calculate second energy contribution
!            call MAT_ADD(beta, x, -1.0E0_realk, delta_i, testmat)
!            !write (LUPRI,*) "testmat in debug_hes_via_energy:"
!            !call MAT_PRINT(testmat, 1, ndim, 1, ndim, LUPRI)
!            call e_appr(testmat, fifoqueue, e2)
!
!            !Calculate third energy contribution
!            call MAT_ADD(-beta, x, -1.0E0_realk, delta_i, testmat)
!            !write (LUPRI,*) "testmat in debug_hes_via_energy:"
!            !call MAT_PRINT(testmat, 1, ndim, 1, ndim, LUPRI)
!            call e_appr(testmat, fifoqueue, e3)
!
!            !Calculate fourth energy contribution
!            call MAT_ADD(-beta, x, 1.0E0_realk, delta_i, testmat)
!            !write (LUPRI,*) "testmat in debug_hes_via_energy:"
!            !call MAT_PRINT(testmat, 1, ndim, 1, ndim, LUPRI)
!            call e_appr(testmat, fifoqueue, e4)
!
!            lintra_elm = (e1-e2+e3-e4)/(4.0E0_realk*alfa*beta)
!            lintra_elm = lintra_elm/4E0_realk !linear transformation is divided by 4
!            diff = lintrans%elms((j-1)*ndim+i)-lintra_elm
!            if (abs(diff) > 1.0E-6_realk) then
!               write(LUPRI,*) 'WARNING: found inconsistency in linear transformation, error =', diff
!               write(LUPRI,*) 'in element (i,j) =', i, j, ' of linear transformation'
!            endif
!            !write(lupri,*) 'Diff between elm from lintrans and elm from energy calculation:', lintrans%elms((j-1)*ndim+i)-lintra_elm
!            !write(lupri,*) 'Energy elms 1, 2, 3, 4:', e1,e2,e3,e4
!            !write(lupri,*) 'Element from linear transformed vector:', lintrans%elms((j-1)*ndim+i)
!            !write(lupri,*) 'Element from energy calculation       :', lintra_elm
!         endif
!      enddo
!   enddo 
!
!   call MAT_FREE(delta_i)
!   call MAT_FREE(testmat)
!   end subroutine debug_hes_via_energy

!Commented out because it it not yet changed to fit new configuration structure. /Stinne
!   subroutine arh_debug_print(reject)
!   implicit none
!        logical, intent(in) :: reject !if info is from a rejection step, we do not want to print all info again
!                                                                                                               
!   if (reject) then
!      write (lupri, "(i5, F36.6, F12.6, '                ///')") &
!           & debug_arh_scfit, final_redspace_eival, debug_arh_arheival
!   else          
!      write (lupri, "(i5, F12.6, F12.6, F12.6, F12.6, F12.6, '    ///')") &
!           & debug_arh_scfit, debug_arh_diag_hlgap, debug_arh_iter_hlgap, final_redspace_eival, &
!           & debug_arh_arheival, debug_arh_heseival
!   endif
!   end subroutine arh_debug_print

!Commented out because it it not yet changed to fit new configuration structure. /Stinne
!   subroutine debug_arh_quasinewton(x,lintra, arhnormal, G1, G2)
!   implicit none
!        TYPE(matrix), intent(inout) :: x,lintra, G1, G2, arhnormal
!        type(matrix)             :: scr1, scr2
!        integer                  :: ndim
!        real(realk)              :: xnorm, sec_norm, arh_norm, arhnormal_norm, graddiff, lintradiff
!
!   if (print_header) then
!      write(lupri, "('    2nd order   ARH          xnorm       graddiff     :::')")
!      print_header = .false.
!   endif
!
!   ndim = G1%nrow
!      
!   call mat_init(scr1, ndim, ndim)
!   call mat_init(scr2, ndim, ndim)
!   
!   call mat_scal(0.25E0_realk,g1) ; call mat_scal(0.25E0_realk,g2)
!   call mat_add(1E0_realk,g1,-1.0E0_realk, g2, scr1)
!
!   graddiff = sqrt(mat_sqnorm2(scr1))
!   !write (lupri,*) 'Grad difference:'
!   !call MAT_PRINT(scr1, 1, scr1%nrow, 1, scr1%ncol, LUPRI)
!   !write (lupri,*) 'Lintra:'
!   !call MAT_PRINT(lintra, 1, scr1%nrow, 1, scr1%ncol, LUPRI)
!
!   call mat_add(1E0_realk,scr1,-1.0E0_realk, lintra, scr2)
!   sec_norm = sqrt(mat_sqnorm2(scr2))
!   call mat_add(1E0_realk,scr1,-1.0E0_realk, arhnormal, scr2)
!   arhnormal_norm = sqrt(mat_sqnorm2(scr2))
!
!   !call mat_add(1E0_realk,scr1,1.0E0_realk, arhlintra, scr2)
!   !arh_norm = sqrt(mat_sqnorm2(scr2))
!
!   !call mat_add(1E0_realk,lintra, 1.0E0_realk, arhlintra, scr2)
!   !lintradiff = sqrt(mat_sqnorm2(scr2))
!
!   xnorm = sqrt(mat_sqnorm2(x))
!
!   write(lupri, "(4E13_realk.5, '  :::')") sec_norm, arhnormal_norm, xnorm, graddiff
!
!   call mat_free(scr1)
!   call mat_free(scr2)
!                                                                                                           
!   end subroutine debug_arh_quasinewton

!Commented out because it it not yet changed to fit new configuration structure. /Stinne
!   subroutine debug_arh_superlin_conv(x,fifoqueue,symm,SCF_it)
!   !Test superlinear convergence
!   !lim k->inf || (HesExact(k) - HesARH(k))X(k) || / ||X(k)|| = 0
!   implicit none
!        TYPE(matrix), intent(in) :: x
!        integer, intent(in)      :: symm, SCF_it
!        TYPE(modFIFO),intent(in) :: fifoqueue
!        type(matrix)             :: lintraARH, lintraExact, lintraRH, scr1
!        integer                  :: ndim
!        real(realk)              :: xnorm, hesnormARH, hesnormRH, ratioARH, ratioRH
!        logical                  :: dummy
!
!      ndim = x%nrow
!      call mat_init(lintraRH, ndim,ndim)
!      call mat_init(lintraARH, ndim,ndim)
!      call mat_init(lintraExact, ndim,ndim)
!      call mat_init(scr1, ndim,ndim)
!
!      !ARH lintrans:
!      call arh_lintrans(x,symm,0.0E0_realk,lintraARH,fifoqueue)
!
!      !Turn off ARH terms and get RH lintrans:
!      dummy = cfg_arhterms 
!      cfg_arhterms = .false.
!      call arh_lintrans(x,symm,0.0E0_realk,lintraRH,fifoqueue)
!      cfg_arhterms = dummy
!
!      !Turn on sec. order terms and get sec. order lintrans:
!      dummy = cfg_do_2nd_order
!      cfg_do_2nd_order = .true.
!      call chol_hessian_times_vector(x,symm,0.0E0_realk,lintraExact) !lintra is E[2](x)
!      cfg_do_2nd_order = dummy
!
!      xnorm = SQRT(mat_sqnorm2(x))
!
!      CALL mat_add(1E0_realk,lintraExact,-1.0E0_realk,lintraARH,scr1)
!      hesnormARH = SQRT(mat_sqnorm2(scr1))
!      ratioARH = hesnormARH/xnorm
!
!      CALL mat_add(1E0_realk,lintraExact,-1.0E0_realk,lintraRH,scr1)
!      hesnormRH = SQRT(mat_sqnorm2(scr1))
!      ratioRH = hesnormRH/xnorm
!
!      write(lupri,"('||ExH - apprH)X||/||X|| = ', F12.8, ' (ARH) ', F12.8, ' (RH)   k = ', i3)") ratioARH, ratioRH, SCF_it
!
!      call mat_free(lintraRH)
!      call mat_free(lintraARH)
!      call mat_free(lintraExact)
!      call mat_free(scr1)
!   end subroutine debug_arh_superlin_conv

end module arh_debugging
