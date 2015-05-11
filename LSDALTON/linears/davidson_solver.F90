MODULE  davidson_solv_mod
use precision
use matrix_module
use matrix_operations
use davidson_settings
use typedef
use decompMod
use arhDensity
use kurtosis
use direct_dens_util
use memory_handling
use orbspread_hess_prec_mod,only: orbspread_hesslin, orbspread_precond
CONTAINS

!> \brief Contains reduced space solver which solves the (level-shifted) Newton equation (Hoeyvik et al., JCTC, 8)
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param CFG  contains information needed to call correct linear transforms etc.
!> \param grad Gradient (type(matrix)) 
!> \param X solution to (H-muI)X = -G as found from reduced space solver
  subroutine davidson_solver(CFG,grad,x) 
    implicit none
    type(matrix),intent(in)  :: grad
    real(realk), pointer :: xred(:)
    type(RedSpaceItem) :: CFG
    type(matrix)       :: x,b_current,res  
    type(matrix)       :: sigma_temp
    logical            :: converged
    real(realk)        :: stepsize,alpha
    real(realk)        :: CurrentResNorm,GdotG,test
    real(realk)        :: resnorm_2D 
    integer            :: coldim,rowdim,iter,i,j
    logical            :: Levelshift
    real(realk)        :: xHx,gx
    real(realk)        :: initial_stepsize
    integer            :: MoreRed
    real(realk) :: minel,arh_thresh=0_realk
    integer :: minel_pos(2)

    ! Initialize reduced space matrices
    call mem_alloc(CFG%Ared,CFG%max_it,CFG%max_it) 
    call mem_alloc(CFG%Gred,CFG%max_it)
    call mem_alloc(CFG%xred,CFG%max_it)
    stepsize = CFG%stepsize
    coldim = grad%ncol
    rowdim = grad%nrow

    initial_stepsize=CFG%stepsize 
    ! Find lowest hessian diagonal, use extra start vector if negative 
    CFG%start_it = 3 !default init
    if ((.not. CFG%arh_davidson)) then
       call mat_min_elm(CFG%P,minel,minel_pos)
       if (minel < 0_realk) CFG%start_it = 4
    else
       if (CFG%arh_extravec) CFG%start_it = 4
       if (CFG%arh_gradnorm .ge. 1E-3_realk) arh_thresh=0.0001_realk
       if (CFG%arh_gradnorm < 1E-3_realk) arh_thresh=0.00001_realk
    end if

    if (grad%ncol .lt. 6) CFG%start_it = 2 ! reduced for few orbs 

    call mem_alloc(CFG%Allb,CFG%max_it)
    call mem_alloc(CFG%AllSigma,CFG%max_it)
    do i=1,CFG%start_it-1
       call mat_init(CFG%Allb(i),grad%nrow,grad%ncol)
       call mat_init(CFG%AllSigma(i),grad%nrow,grad%ncol)
    end do
    call mat_init(b_current,rowdim,coldim)
    call mat_init(sigma_temp,rowdim,coldim)
    call FirstTrialVec(grad,CFG,b_current,sigma_temp)
    call InitializeCFG(grad,CFG,b_current,sigma_temp)
    if (CFG%start_it .gt. 2) then 
       call SecondTrialVec(grad,CFG,b_current,sigma_temp) 
       !Construct rest of Ared in 3D starting space
       call IncreaseDimRedSpace(grad,CFG,b_current,2)
       if (CFG%start_it ==4) then
          call ThirdTrialVec(grad,CFG,b_current,sigma_temp,minel_pos,CFG%start_it)
          ! Check if coupling to b1 and b2 is strong enough
          if (CFG%singularity) then
             CFG%start_it = 3
             call mat_free(CFG%Allb(3))
             call mat_free(CFG%AllSigma(3))
             CFG%singularity = .false.
          endif
       end if
    endif
    call mat_free(b_current)
    call mat_free(sigma_temp)
    !************************************************************
    !*             Start Reduced Space Loop                     *
    !************************************************************
    converged =.false.
    call mat_init(res,rowdim,coldim)
    do MoreRed=1,5
       resnorm_2D =0_realk
       ReducedSpaceLoop: do iter=CFG%start_it,CFG%max_it-1
          call mem_alloc(xred,iter-1)
          call mat_init(CFG%Allb(iter),grad%nrow,grad%ncol)
          call mat_init(CFG%AllSigma(iter),grad%nrow,grad%ncol)

          ! Determine level shift 
          call get_levelshift(grad,CFG,X,xred,iter,alpha,Levelshift)

          ! Solve linear equations in reduced space
          if (Levelshift) then
             call GetResidual(grad,CFG%mu,iter,Res,X,xred,CFG)
          else
             call SolveNewton_ReducedSpace(grad,CFG,X,xred,res,iter)
          end if
          CFG%xred(1:iter-1)=xred(:)

          CurrentResNorm = sqrt(mat_dotproduct(res,res))
          ! PRINT INFO FOR MICRO ITERATION
          if (CFG%arh_davidson .and. CFG%arh_davidson_debug) then
             write(CFG%lupri,'(a,i3,a,ES13.5,a,ES13.5,a,ES13.5)') "iter :",iter, "    mu :",&
                  &CFG%mu,  "    ResNorm :",CurrentResNorm,&
                  &" Gradient norm",CFG%arh_gradnorm
          elseif (.not. CFG%arh_davidson) then
             if (iter > 3) then
                call test_convergence(CFG,grad,resnorm_2d)
             else
                resnorm_2d = sqrt(mat_sqnorm2(grad))
             endif
             if (CFG%orb_debug) then
                write(CFG%lupri,'(a,i3,a,ES13.5,a,ES13.5,a,ES13.5)') "iter :",iter, "    mu :",&
                     &CFG%mu,  "    ResNorm :",CurrentResNorm,&
                     &"   2D ResNorm",resnorm_2d
             end if
          end if

          !TEST CONVERGENCE
          if (CFG%arh_davidson .and. (CurrentResNorm < arh_thresh*CFG%arh_gradnorm .or. CurrentResNorm<1.0E-8_realk)) then
             call mem_dealloc(xred)
             CFG%it =iter
             converged =.true.
             exit
          elseif ( (CurrentResNorm < CFG%conv_thresh*resnorm_2D).or. (CurrentResNorm < 1E-6_realk .and.&
               &(.not. CFG%arh_davidson))) then
             call mem_dealloc(xred)
             CFG%it =iter
             converged =.true.
             exit    
          end if


          call mat_init(b_current,rowdim,coldim)
          !b_current=res
          if (CFG%precond) then
             call Precondition(CFG,res,b_current,grad)
          else
             call mat_assign(b_current,res)
          end if
          call orthonormalize(b_current,iter-1,CFG)
          call mat_init(sigma_temp,rowdim,coldim)
          call LinearTransformations(CFG,b_current,sigma_temp,grad)
          !Save matrices 
          call mat_assign(CFG%AllSigma(iter),sigma_temp)
          call mat_free(sigma_temp)
          call mat_assign(CFG%Allb(iter),b_current)
          !Increase dimension
          call IncreaseDimRedSpace(grad,CFG,b_current,iter)
          call mat_free(b_current)

          call mem_dealloc(xred)
          CFG%it = iter


       end do ReducedSpaceLoop
       ! Check if micro iterations have converged
       if (converged) then
          if (CFG%orb_debug) write(CFG%lupri,'(a,i4)') &
               & 'iter : Micro iterations converged in RedSpaceLoop number :', MoreRed
          exit
       elseif ((.not. converged) .and. CFG%stepsize > 0.0009_realk .and. (MoreRed .ne. 5)) then 
          CFG%stepsize = 0.5_realk*CFG%stepsize
          if (CFG%orb_debug) write(CFG%lupri,'(a,i4,a,f9.4)') &
               &'iter : Micro iterations not converged for RedSpaceLoop:', MoreRed, &
               & " new stepsize : ",CFG%stepsize
          do i=CFG%start_it,CFG%it
             call mat_free(CFG%AllSigma(i))
             call mat_free(CFG%Allb(i))
          end do
       elseif ((CFG%stepsize <0.0009_realk).or. (MoreRed==5)) then
          if (CFg%arh_davidson) resnorm_2d=CFG%arh_gradnorm
          write(CFG%lupri,*) 
          write(CFG%lupri,*) 'WARNING!  Did not converge micro iterations. Residual norm reduction is ',&
               &CurrentresNorm/resnorm_2d
          write(CFG%lupri,*) 
          exit
       end if
    end do
    call mat_free(res)
    ! ********************************
    !    END OF REDUCED SPACE LOOP   *
    ! ********************************



    !Compute RedHessian eigenvalues and print
    if (CFG%arh_davidson .and. CFG%arh_davidson_debug) call RedHessian_eigvals(CFG) 

    call predicted_change(CFG,grad,x)

    call mem_dealloc(CFG%Ared)
    call mem_dealloc(CFG%Gred)
    call mem_dealloc(CFG%xred)

    do i=1,CFG%it
       call mat_free(CFG%Allb(i))
    end do
    do i=1,CFG%it
       call mat_free(CFG%AllSigma(i))
    end do
    call mem_dealloc(CFG%Allb)
    call mem_dealloc(CFG%AllSigma)

  end subroutine davidson_solver


!> \brief Subroutine that takes in type(Matrix) and normalize_mats it
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param M matrix (type(matrix)) to be normalized
  subroutine normalize_mat(M)
    implicit none
    type(matrix), intent(inout) :: M
    real(realk)		      :: norm

    norm = sqrt(mat_sqnorm2(M))
    call mat_scal(1.0_realk/norm,M)

  end subroutine normalize_mat

!> \brief Subroutine that performs linear transformation HX 
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param X trial vector (input)
!> \param HX resulting linear transformation (output)
!> \param G gradient (input)
!> \param CFG reduced space item (input) 
  subroutine LinearTransformations(CFG,X,HX,G)
    implicit none
    type(RedSpaceItem) :: CFG
    type(MATRIX),intent(in) :: X,G
    type(matrix)            :: HX
    integer :: n
    if (CFG%arh_davidson .and. CFG%arh_lintrans) then
       call arh_lintrans(CFG%arh,CFG%decomp,X,CFG%symm,0.0_realk,HX,CFG%fifoqueue)
       return
    end if

    if (CFG%orbspread) then
       call orbspread_hesslin(HX,X,0.0_realk,X%nrow,CFG%orbspread_inp)
       return
    end if

    if (CFG%PFM) then
       call compute_lintrans(CFG%PFM_input,X,G,HX)
       return
    end if


    if (CFG%PM_input%PipekMezeyLowdin) then
       call PMLowdin_LinTra(G,HX,X,CFG%PM_input)
    elseif (CFG%PM_input%PipekMezeyMull) then
       call PMMull_LinTra(G,HX,X,CFG%PM_input) 
    end if


  end subroutine LinearTransformations

!> \brief Compute linear equations in already set up reduced space (used for zero level shift)
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param grad gradient (input)
!> \param CFG Reduced space item (input)
!> \param X solution vector in full space (output)
!> \param xred solution in reduced space (output)
!> \param red residual (output)
!> \param iter dimension of augmented reduced Hessian (input)
  subroutine SolveNewton_ReducedSpace(grad,CFG,X,xred,res,iter)
    implicit none
    type(RedSpaceItem) :: CFG
    integer            :: iter
    type(matrix)       :: x,bn,sigma_n,res
    real(realk)        :: xred(iter-1)
    type(matrix),intent(in)  :: grad
    integer            :: i
    real(realk),pointer:: A(:,:)
    real(realk)        :: RHS(iter-1)
    integer            :: IPIV(iter-1),IERR
    IERR=0

    call mem_alloc(A,iter-1,iter-1)
    A(1:iter-1,1:iter-1)= CFG%Ared(2:iter,2:iter)
    RHS(1:iter-1) = -CFG%Ared(2:iter,1)



    !solve linear equations Hredxred=-Gred
    call DGESV(iter-1,1,A,iter-1,IPIV,RHS,iter-1,IERR)
    if (IERR .ne. 0) STOP 'DGESV in SolveNewton'

    xred=RHS

    !Obtain x and residual in real space
    !Res = Hx+G = sum_i xred(i)*sigma_i + G
    call mat_zero(res); call mat_zero(X)

    do i=1,iter-1
       !obtain x
       call mat_daxpy(xred(i),CFG%Allb(i),x)
       !obtain residual
       call mat_daxpy(xred(i),CFG%AllSigma(i),res)
    end do

    !Add G to residual
    call mat_daxpy(1.0_realk,grad,res)
    call mem_dealloc(A)


  end subroutine SolveNewton_ReducedSpace


!> \brief Increase dimension of reduced space augmented Hessian
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param grad gradient (input)
!> \param CFG Reduced space input
!> \param b_current trial vector for current iteration (input)
!> \param iter current iteration (input)
  subroutine IncreaseDimRedSpace(grad,CFG,b_current,iter)
    implicit none
    type(RedSpaceItem) :: CFG
    integer,intent(in) :: iter
    type(matrix),intent(in) :: b_current
    type(matrix)       :: bn
    type(matrix),intent(in)::grad
    type(matrix)       :: sigma_n
    integer            :: i,j

    CFG%Ared(1,iter+1) = 0.0_realk
    CFG%Ared(iter+1,1) = CFG%Ared(1,iter+1)
    !Set up reduced augmented hessian
    do i=1,iter
       CFG%Ared(iter+1,i+1) = mat_dotproduct(b_current,CFG%AllSigma(i)) 
       CFG%Ared(i+1,iter+1) = mat_dotproduct(CFG%Allb(i),CFG%AllSigma(iter))
    end do


    if (CFG%start_it == 4 .and. iter==3) then
       if ((dabs(CFG%Ared(2,4))+dabs(CFG%Ared(3,4))) < 1.0E-3_realk) then  
          !Must abandon b3  due to no coupling to b1 and b2
          CFG%singularity=.true.
       endif
    end if


  end subroutine IncreaseDimRedSpace

!> \brief Routine which calls correct preconditioning routines
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param res Residual to be preconditioned (input)
!> \param b_current preconditioned residual (output)
!> \param grad gradient (input)
  subroutine Precondition(CFG,res,b_current,grad)
    implicit none
    type(matrix)       :: H,res,b_current
    type(matrix),intent(in)::grad
    type(RedSpaceItem) :: CFG

    if (CFG%arh_davidson .and. CFG%arh_precond) then
       call arh_precond(CFG%arh,CFG%decomp,res,CFG%symm,CFG%mu,b_current)
       return
    end if

    if (CFG%PFM) then 
       call kurtosis_precond(b_current,res,CFG%mu,CFG%PFM_input)
       return
    end if

    if (CFG%orbspread) then
       call orbspread_precond(b_current,res,CFG%mu,CFG%orbspread_inp)
    elseif (CFG%PM) then
       call charge_precond(b_current,res,CFG%mu,CFG%PM_input)
    end if

  end subroutine Precondition


!> \brief Subroutine that takes V and projects out U component.
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param V input/output matrix 
!> \param U component to be projected out
  subroutine Orthogonalize(V,U)
    implicit none
    type(matrix),intent(in) :: U
    type(matrix),intent(inout) :: V
    real(realk) :: factor,norm2

    norm2 = mat_sqnorm2(U)
    factor = mat_dotproduct(V,U)/norm2
    !Make V = V -factor*U
    call mat_daxpy(-factor,U,V)

  end subroutine Orthogonalize

!> \brief Computes residual from reduced space solution
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param grad gradient (input)
!> \param mu level shift (input)
!> \param iter current iteration (dimension of augmented RS Hessian) (input)
!> \param Res  residual (output)
!> \param X solution in full space (output)
!> \param Xred solution in reduced space (input)
  subroutine GetResidual(grad,mu,iter,Res,X,xred,CFG)
    implicit none
    type(RedSpaceItem)  :: CFG
    integer, intent(in) :: iter
    real(realk),intent(in) :: xred(iter-1)
    real(realk),intent(in) :: mu
    type(matrix) :: Res,X
    type(matrix),intent(in)::grad
    integer :: i

    call mat_zero(Res);call mat_zero(X)

    do i=1,iter-1
       !Obtain x
       call mat_daxpy(xred(i),CFG%Allb(i),x)
       !obtain residual
       call mat_daxpy(xred(i),CFG%AllSigma(i),res)
       call mat_daxpy(-xred(i)*mu,CFG%Allb(i),res)
    end do

    !Add G to residual
    call mat_daxpy(1.0_realk,grad,Res)


  end subroutine GetResidual

!> \brief Subroutine that generates first trial vector for reduced space solver
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param grad gradient (input)
!> \param b_current 1st trial vector, normalized grad. (output)
!> \param sigma_temp 1st linear transformation (output)
  subroutine FirstTrialVec(grad,CFG,b_current,sigma_temp)
    implicit none
    type(RedSpaceItem) :: CFG
    type(matrix)       :: b_current, sigma_temp,temp
    type(matrix),intent(in)  :: grad

    call mat_assign(b_current ,grad)
    call normalize_mat(b_current)
    call LinearTransformations(CFG,b_current,sigma_temp,grad)
    call mat_assign(CFG%Allb(1),b_current)
    call mat_assign(CFG%AllSigma(1),sigma_temp)


  end subroutine FirstTrialVec

!> \brief  Subroutine that generates second trial vector
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param grad gradient (input)
!> \param b_current 2nd trial vector (output)
!> \param sigma_temp 2nd linear transformation (output)
  subroutine SecondTrialVec(grad,CFG,b_current,sigma_temp)
    implicit none
    type(RedSpaceItem) :: CFG
    type(matrix),intent(in)       :: grad
    type(matrix)       :: b_current,sigma_temp,temp
    type(matrix) :: scrmat
    real(realk)  :: tcpu,twall

    call mat_init(temp,grad%nrow,grad%ncol)
    call mat_copy(-1.0_realk,grad,temp)
    call LinearTransformations(CFG,temp,b_current,grad) 
    call mat_daxpy(-1.0_realk,grad,b_current)

    call orthonormalize(b_current,1,CFG)
    call LinearTransformations(CFG,b_current,sigma_temp,grad) 

    call mat_assign(CFG%Allb(2),b_current)
    call mat_assign(CFG%AllSigma(2),sigma_temp)

    call mat_free(temp)

  end subroutine SecondTrialVec

!> \brief Creat third trial vector for reduced space
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param G gradient (input)
!> \param b 3rd trial vector  (output)
!> \param sigma 3rd linear transformation
  subroutine ThirdTrialVec(G,CFG,b,sigma,pos,start_it)
    implicit none
    type(RedSpaceItem) :: CFG
    type(matrix) ::b,G,sigma
    integer :: indx(2)
    integer :: pos(2)
    integer :: i,j,norb,ndim,start_it
    real(realk) :: min_diag
    real(realk),pointer :: trial(:,:)

    norb= G%ncol

    call mat_zero(b)
    ndim=b%ncol

    if (CFG%arh_davidson) then

       call ComputeFock_eigvecs(CFG%decomp,b,CFG%Allb(1),ndim,CFG%lupri,CFG%arh_dodft,CFG%arh_extravec) 
       if (.not. CFG%arh_extravec) then
          call mat_free(CFG%Allb(3))
          call mat_free(CFG%AllSigma(3))
          start_it = 3
          return
       endif
       call orthonormalize(b,2,CFG)
       call LinearTransformations(CFG,b,sigma,G)
       call mat_assign(CFG%Allb(3),b)
       call mat_assign(CFG%AllSigma(3),sigma)
       call IncreaseDimRedSpace(G,CFG,b,3)

    else
       call mem_alloc(trial,1,1)
       trial(1,1)=1.0_realk
       call mat_create_block(b,trial,1,1,Pos(1),Pos(2))
       trial(1,1)=-1.0_realk
       call mat_create_block(b,trial,1,1,Pos(2),Pos(1))
       call orthonormalize(b,2,CFG)
       call LinearTransformations(CFG,b,sigma,G)
       call mat_assign(CFG%Allb(3),b)
       call mat_assign(CFG%AllSigma(3),sigma)
       call IncreaseDimRedSpace(G,CFG,b,3)
       call mem_dealloc(trial)
    endif


  end subroutine ThirdTrialVec


!> \brief Subroutine that initializes matrices Ared,Sred and Gred
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param grad (input)
!> \param b_current 1st trial vector
!> \param sigma_temp 1st linear transformation
  subroutine InitializeCFG(grad,CFG,b_current,sigma_temp)
    implicit none
    type(matrix)        :: sigma_temp,b_current
    type(RedSpaceItem)  :: CFG
    type(matrix),intent(in)        :: grad


    CFG%Ared = 0.0_realk; CFG%Gred = 0.0_realk;CFG%xred=0.0_realk
    !Construct 1st part of reduced augmented hessian
    CFG%Gred(1) = mat_dotproduct(b_current,grad)
    CFG%Ared(1,2) = CFG%Gred(1)
    CFG%Ared(2,1) = CFG%Ared(1,2)
    CFG%Ared(2,2) = mat_dotproduct(b_current,sigma_temp)


  end subroutine InitializeCFG


!> \brief Subroutine that determines levelshift by using a bisectional algorithm
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param grad gradient (input)
!> \param X solution in full space (output)
!> \param Xred solution in reduced space (output)
!> \param alpha parameter to determine levelshift 
!> \param CFG%mu level shift, main output
!> \param Levelshift logical, F is no levelshift, T if level shift is needed 
  subroutine get_levelshift(grad,CFG,X,xred,iter,alpha,Levelshift)
    implicit none
    integer, intent(in) ::iter
    type(RedSpaceItem)  :: CFG
    real(realk),intent(inout)  :: xred(iter-1)
    type(matrix) :: X
    type(matrix),intent(in) ::grad
    real(realk)  :: alpha1,alpha2,alpha
    logical      :: done,Levelshift
    !initialize levelshift
    Levelshift = .true.

    !Find interval which contains alpha giving the correct stepsize
    call BracketAlpha(CFG,X,xred,iter,alpha1,alpha2,done,Levelshift)
    if (.not. done) STOP 'Something went wrong in BracketAlpha'

    !USrede interval to determine alpha and thus mu
    if (Levelshift) then
       call LocateAlpha(CFG,X,xred,iter,alpha1,alpha2,alpha)
    else
       alpha = 0.0_realk
    end if


  end subroutine get_levelshift



!> \brief Bracket interval of alpha in which alpha corresponding to correct stepsize 
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param alpha1, alpha2 (output), brackets interval with correct alpha
!> \param Levelshift (output) returned as F is Hessian is positive definite
  subroutine BracketAlpha(CFG,X,xred,iter,alpha1,alpha2,done,Levelshift)
    implicit none
    type(RedSpaceItem)     :: CFG
    integer,intent(in)     :: iter
    real(realk),intent(inout) :: xred(iter-1)
    real(realk)            :: alpha1,alpha2,alpha
    logical                :: done
    type(matrix)           :: X
    integer,parameter :: N=500
    integer           :: j
    !f1,f2: Function values at alpha1 and alpha2
    real(realk)       :: f1,f2,factor,alpha_max
    logical :: Levelshift


    alpha_max=10000000.0_realk
    factor = 1.0_realk

    !Initialize values for alpha1 and alpha2
    alpha1=0.00001_realk

    alpha2= alpha1+factor 

    done = .false.
    call NormXminusStep(CFG,X,xred,iter,alpha1,f1)
    call NormXminusStep(CFG,X,xred,iter,alpha2,f2)

    do j=1,N
       if (f1*f2 < 0.0_realk) then
          done = .true.
          exit
       end if

       alpha2 = alpha2+factor
       call NormXminusStep(CFG,X,xred,iter,alpha2,f2)
       if (alpha2 > alpha_max) then
          if  (f1<0.0_realk .and. f2<0.0_realk)  then
             CFG%mu = 0.0_realk
             Levelshift =.false.
             done = .true.
             exit
          else
             Stop 'Something is wrong. f1 and f2 are both positive in (0.001  50000)'
          end if
       end if

       factor = factor*10.0_realk
    end do


  end subroutine BracketAlpha

!> \brief Function that determined correct alpha
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param X solution in full space (output)
!> \param xred solution in reduced space (output)
!> \param alpha1 & alpha2 (input) brackets interval with correct alpha
!> \param alpha correct alpha giving the correct mu (output)
  subroutine LocateAlpha(CFG,X,xred,iter,alpha1,alpha2,alpha) 
    implicit none
    type(RedSpaceItem)     :: CFG
    integer, intent(in)    :: iter
    real(realk)            :: alpha,alpha1,alpha2, func
    type(matrix)           :: X
    real(realk)            :: xred(iter)
    !> MaxBS is maximum allowed number of bisetions
    integer, parameter     :: MaxBS = 100 
    !> Thresh defines the accuracy to which alpha is found (alpha +/- thresh)
    real(realk),parameter  :: thresh=0.0001_realk
    real(realk)            :: dx,f,f_mid,alpha_mid,tcpu1,twall1
    integer :: j
    call NormXminusStep(CFG,X,xred,iter,alpha1,f)
    call NormXminusStep(CFG,X,xred,iter,alpha2,f_mid)

    if (f*f_mid > 0.0_realk) STOP 'Alpha is not bracketed by alpha1,alpha2'

    ! Orient search so that f>0 lies at x+dx
    if (f < 0.0_realk) then
       alpha = alpha1
       dx    = alpha2-alpha1
    else
       alpha = alpha2
       dx    = alpha1-alpha2
    end if

    !Bisectional loop
    do j = 1,MaxBS
       dx        = dx*0.5_realk
       alpha_mid = alpha + dx
       call NormXminusStep(CFG,X,xred,iter,alpha_mid,f_mid)

       if (f_mid < 0.0_realk)  alpha = alpha_mid
       if ((dabs(dx) < thresh) .or. (dabs(f_mid) < thresh)) exit

       if (j==MaxBS) then
          STOP  'Bisectional algorithm did not converge in maximum iterations'
       end if
    end do

  end subroutine LocateAlpha

!> \brief Solve the eigen value equation in the already set up reduced space 
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param X solution in full space (output)
!> \param xred solution in reduced space (output)
  subroutine Solve_EigVal_RedSpace(CFG,X,xred,iter,alpha)
    implicit none
    type(RedSpaceItem)  :: CFG
    integer, intent(in) :: iter
    real(realk)   :: alpha
    type(Matrix)  :: bn,X
    real(realk)   :: xred(iter-1)
    real(realk),pointer  :: A(:,:)
    real(realk)   :: EigValues(iter),mu_Vec(iter),dummy(iter)
    real(realk),pointer :: WORK(:),VL(:,:),VR(:,:)
    integer :: LWORK,i,IERR,j,indx(1)
    character(len=1) :: V,L
    IERR=0

    call mem_alloc(A,iter,iter)
    !Set up reduced Hessian and overlap  with correct dimension
    A(1:iter,1:iter) = CFG%Ared(1:iter,1:iter)
    !Scale first column&row of aug hessian by alpha
    A(1,2) = alpha*A(1,2)
    A(2,1) = alpha*A(2,1)
    LWORK = 5*(iter)
    call mem_alloc(WORK,LWORK)
    call mem_alloc(VL,iter,iter)
    call mem_alloc(VR,iter,iter)
    call DGEEV('N','V',iter,A,iter,EigValues,dummy,VL,iter,VR,iter,WORK,LWORK,IERR)
    if (IERR .ne. 0) call lsquit('Something went wrong in DGEEV in orbital localization',CFG%lupri)
    call mem_dealloc(WORK)
    call mem_dealloc(VL)

    CFG%mu = minval(EigValues)
    indx= minloc(EigValues)
    !Store corresponding eigenvector (eig.vecs is given i A in output)
    mu_Vec(:) = VR(:,indx(1))
    call mem_dealloc(VR)
    !scale vector such that first element is 1
    if (dabs(mu_Vec(1)).le. 1E-15_realk) then
       if (CFG%arh_davidson) then
          write(CFG%lupri,*)
          write(*,*)
          write(CFG%lupri,*) ' ***** WARNING ***** '
          write(*,*) ' ***** WARNING ***** '
          write(CFG%lupri,*) '  A trial vector contains zero gradient component' 
          write(*,*) '  A trial vector contains zero gradient component' 
          write(*,*) '  Start calculation using keyword .ARH instead of .ARH DAVID /.ARH(LS) DAVID' 
          write(CFG%lupri,*) '  Start calculation using keyword .ARH instead of .ARH DAVID /.ARH(LS) DAVID' 
          call lsquit('Singularities due to trial vector containg zero gradient component',CFG%lupri)
       end if
    end if
    mu_Vec = mu_Vec/mu_Vec(1)

    !Xred is eigenvector  without first element
    xred(1:iter-1) = mu_Vec(2:iter) 
    xred = xred/alpha

    !Obtain X and in real space by X = sum_n xred(n)*bn 
    call mat_zero(X)
    do i=1,iter-1
       call mat_daxpy(xred(i),CFG%Allb(i),X)
    end do
    call mem_dealloc(A)

  end subroutine Solve_EigVal_RedSpace

!> \brief Compute  Norm(X)-stepsize (want to find Norm(X) = stepsize)
!> \author Ida-Marie Hoeyvik 
!> \date 2012 
!> \param X solution in full space (input)
!> \param xred solution in reduced space (output)
!> \param StepDev How far from Norm(X) = stepsize (output)
  subroutine NormXminusStep(CFG,X,xred,iter,alpha,StepDev)
    implicit none
    type(RedSpaceItem) :: CFG
    integer,intent(in) :: iter
    real(realk)        :: alpha
    real(realk)        :: xred(iter-1),StepDev
    type(matrix)       :: X

    call Solve_EigVal_RedSpace(CFG,X,xred,iter,alpha)

    StepDev = dsqrt(mat_sqnorm2(X))-CFG%stepsize


  end subroutine NormXminusStep


!> \brief Gram-Schmidt scheme to orthonormalize vectors
!> \author Ida-Marie Hoeyvik
!> \date 2012
!> \param b_current Vector to be orthogonalized against the n vectors in CFG
  subroutine orthonormalize(b_current,n,CFG)
    implicit none
    type(RedSpaceItem) :: CFG
    type(matrix) :: b_current
    integer      :: i,j 
    real(realk)  :: dotprod
    integer, intent(in) :: n

    call normalize_mat(b_current)

    do j=1,2
       do i=1,n
          dotprod = mat_dotproduct(CFG%Allb(i),b_current)
          call mat_daxpy(-dotprod,CFG%Allb(i),b_current)
       end do
       call normalize_mat(b_current)
    end do

  end subroutine orthonormalize

!> \brief Solve (H-mu)X=-g ---> Ax=b in 2d space with given mu
  subroutine test_convergence(CFG,grad,resnorm_2D)
    implicit none
    integer,parameter    :: ndim=2
    type(RedSpaceItem)  :: CFG
    real(realk)         :: A(ndim,ndim),b(ndim),id(ndim,ndim)
    integer             :: IPIV(ndim),IERR,i,j
    type(matrix)        :: residual_2D
    type(matrix),intent(in) :: grad
    real(realk)         :: resnorm_2D
    IERR=0

    do j=1,ndim
       do i=1,ndim
          A(i,j) = CFG%Ared(i+1,j+1)
       end do
    end do

    id = 0.0_realk
    do i=1,ndim
       id(i,i) = 1.0_realk
    end do

    !Construct H-mu*I
    A=A-CFG%mu*id

    !construct b=-Gred
    b=0.0_realk
    b(1)=-CFG%Gred(1)

    call DGESV(ndim,1,A,ndim,IPIV,b,ndim,IERR)
    if (IERR .ne. 0) call lsquit('Something wrong, DGESV, test_convergence, orbital localization',CFG%lupri)


    !Compute residual
    call mat_init(residual_2D,grad%nrow,grad%ncol)
    call mat_zero(residual_2D)
    do i=1,2
       !obtain residual
       call mat_daxpy(b(i),CFG%AllSigma(i),residual_2D)
       call mat_daxpy(-b(i)*CFG%mu,CFG%Allb(i),residual_2D)
    end do
    call mat_daxpy(1.0_realk,grad,residual_2D)
    resnorm_2D = sqrt(mat_dotproduct(residual_2D,residual_2D))

    call mat_free(residual_2D)
  end subroutine test_convergence



!> \brief Compute xHx in reduced space
!> \author Ida-Marie Hoeyvik
!> \param CFG  contains Xred, Ared and dimension (input)
!> \param xHx  output
!> \param D  dimension of reduced space
!> \param  Hred  Hessian in basis of D trialvectors
  subroutine Get_xHx_RedSpace(CFG,x,xHx,grad)
    implicit none
    type(RedSpaceItem)       :: CFG
    real(realk)              :: xHx
    integer                  :: i,D
    real(realk),pointer      :: xH(:),xred(:),Hred(:,:)
    type(matrix)             :: grad,X,HX
    integer :: n

    D = CFG%it-1
    call mem_alloc(xH,D)
    call mem_alloc(xred,D)
    call mem_alloc(Hred,D,D)
    do i=1,D
       xred(i) = CFG%xred(i)
    end do
    Hred(1:D,1:D) = CFG%Ared(2:D+1,2:D+1)
    call mat_init(Hx,x%nrow,x%ncol)
    call LinearTransformations(CFG,x,Hx,grad)
    xHx=mat_dotproduct(X,Hx)
    call mat_free(Hx)
    call mem_dealloc(xH)
    call mem_dealloc(Xred)
    call mem_dealloc(Hred)
  end subroutine Get_xHx_RedSpace

!> \brief Compute gx in reduced space
!> \author Ida-Marie Hoeyvik
!> \param CFG  contains Xred, Ared and dimension (input)
!> \param gx  output
  subroutine Get_gx_RedSpace(CFG,X,grad,gx)
    implicit none
    type(RedSpaceItem)       :: CFG
    real(realk), intent(out) :: gx
    integer                  :: i,D
    real(realk),pointer      :: xred(:),gred(:)
    real(8)                  :: ddot 
    type(matrix),intent(in)  :: grad,x
    integer :: n


    D=CFG%it-1
    call mem_alloc(gred,D)
    call mem_alloc(xred,D)

    do i=1,D
       xred(i)=CFG%xred(i)
    end do
    gred(1:D)=CFG%Gred(1:D)

    gx= mat_dotproduct(grad,x)

    call mem_dealloc(xred)
    call mem_dealloc(gred)

  end subroutine Get_gx_RedSpace


  subroutine predicted_change(CFG,grad,x)
    implicit none
    type(RedSpaceItem)  :: CFG
    real(realk)         :: xHx,gx
    type(matrix),intent(in) ::x,grad

    call Get_xHx_RedSpace(CFG,x,xHx,grad)
    call Get_gx_RedSpace(CFG,x,grad,gx)

    CFG%r_denom = gx + 0.5_realk*xHx

  end subroutine predicted_change

!> \brief Compute and print eigenvalues for reduced Hessian 
!> \author Ida-Marie Hoeyvik
!> \date 2012
  subroutine RedHessian_eigvals(CFG)
    implicit none
    type(RedSpaceItem) :: CFG
    real(realk),pointer :: RedH(:,:),eigvals(:),wrk(:)
    integer :: lwrk,IERR,i,j,n
    IERR=0

    n=CFG%it-1
    lwrk=10*n
    call mem_alloc(RedH,n,n)
    call mem_alloc(wrk,lwrk)
    call mem_alloc(eigvals,n)

    do i=1,n
       do j=1,n
          RedH(j,i)=CFG%Ared(j+1,i+1)
       end do
    end do

    call dsyev('N','L',n,RedH,n,eigvals,wrk,lwrk,IERR)
    if (IERR .ne. 0)  call lsquit('Something went wrong in DSYEV in orbital localization',CFG%lupri)


    write(CFG%lupri,*)
    write(CFG%lupri,'(a)')        'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    write(CFG%lupri,'(a,ES13.4,a)') 'X LOWEST EIGENVALUE OF REDUCED HESSIAN: ', eigvals(1), ' X'
    if (CFG%debug_info) print*, 'Lowest eigval redspace :', eigvals(1)
    write(CFG%lupri,'(a)')        'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    write(CFG%lupri,*)

    call mem_dealloc(eigvals)
    call mem_dealloc(RedH)
    call mem_dealloc(wrk)

  end subroutine RedHessian_eigvals


  subroutine ComputeFock_eigvecs(decomp,b,G,ndim,lupri,DFT,success)
    implicit none
    type(DecompItem) :: decomp
    integer          :: ndim
    type(matrix)  :: b
    type(matrix),intent(in) :: G
    real(realk),pointer    :: eigvalsO(:),eigvalsV(:)
    real(realk),pointer    :: FmatO(:,:),FmatV(:,:)
    real(realk),pointer    :: wrk(:)
    real(realk) :: sortocc(100,100),sortvirt(100,100),gap(100,100)
    integer :: lwrk,IERR,loc(2),locO(1),locV(1)
    integer :: i,j,indxO,indxV,lupri
    logical :: success,DFT
    real(realk) ::thr,nrm
    integer :: maxdim,maxdimV,maxdimO

    call mem_alloc(FmatO,ndim,ndim)
    call mem_alloc(FmatV,ndim,ndim)
    call mat_to_full(decomp%FUP,1.0_realk,FmatO)
    call mat_to_full(decomp%FUQ,1.0_realk,FmatV)
    call mem_alloc(eigvalsO,ndim)
    call mem_alloc(eigvalsV,ndim)


    lwrk=-1
    call mem_alloc(wrk,5)
    call dsyev('V','L',ndim,FmatO,ndim,eigvalsO,wrk,lwrk,IERR)

    lwrk=wrk(1)
    call mem_dealloc(wrk)
    call mem_alloc(wrk,lwrk)

    ! DIAGONALIZE OCCUPIED PART
    call dsyev('V','L',ndim,FmatO,ndim,eigvalsO,wrk,lwrk,IERR)
    if (IERR .ne. 0) STOP 'diag FUP error'

    ! COMPUTE VIRT EIGENVECTORS
    call dsyev('V','L',ndim,FmatV,ndim,eigvalsV,wrk,lwrk,IERR)
    if (IERR .ne. 0) STOP 'diag FUQ error'
    call mem_dealloc(wrk)

    ! Sort ten lowest virt eigvals and ten highest occ eigvals
    ! Save their locations, not their values
    thr=1.0E-10_realk
    maxdimO=ndim
    maxdimV=ndim
    do i=1,100
       locO=maxloc(eigvalsO,MASK=abs(eigvalsO)>thr)
       if (locO(1) == 0 .and. i>1) then
          maxdimO=i-1 
          maxdimV=i-1
          exit
       endif
       sortocc(i,1)=dble(locO(1))
       sortocc(i,2)=eigvalsO(locO(1))
       eigvalsO(locO(1))=0.0_realk
       locV=minloc(eigvalsV,MASK=abs(eigvalsV)>thr)
       if (locV(1) == 0 .and. i>1) then
          maxdimV=i-1 
          maxdimO=i-1 
          exit
       endif
       sortvirt(i,1)=dble(locV(1))
       sortvirt(i,2)=eigvalsV(locV(1))
       eigvalsV(locV(1))=0.0_realk
    end do
    ! Find combinations e_virt-e_occ in matrix
    do i=1,min(maxdimV,100)
       do j=1,min(maxdimO,100)
          gap(i,j)=sortvirt(i,2)-sortocc(j,2)
       end do
    end do
    write(lupri,*) ' HOMO-LUMO gap calc. by diagonalization ', minval(gap)
    nrm = 0.0_realk
    maxdim=min(100*100,maxdimV*maxdimO)
    do i=1,maxdim
       loc=minloc(gap)
       indxV=nint(sortvirt(loc(1),1))
       indxO=nint(sortocc(loc(2),1))
       call mat_zero(b)

       ! In case of DFT, must perturb vectors
       if (DFT) then
          do j=1,ndim
             if (FmatO(j,indxO)>1.0E-12_realk) FmatO(j,indxO)=FmatO(j,indxO) +0.001
             if (FmatO(j,indxO)<1.0E-12_realk) FmatO(j,indxO)=FmatO(j,indxO) -0.001
          end do
       end if
       call mat_dger(1.0_realk,FmatO(:,indxO),FmatV(:,indxV),b)
       call mat_dger(-1.0_realk,FmatV(:,indxV),FmatO(:,indxO),b)


       nrm=abs(mat_dotproduct(b,G))
       if (nrm> 1.0E-4_realk) then
          write(lupri,*) 'Special trial vector: orbital energy difference', minval(gap)
          write(lupri,*) "Special trial vector: successfull! Non-orthogonal to gradient", i
          exit
       else
          gap(loc(1),loc(2))=1.0E10_realk
       endif
       if (i==maxdim) write(lupri,*)"Special trial vector: FAILED TO FIND NON-ORTHOGONAL VEC"
       success = .false.
    end do


    call mem_dealloc(FmatO)
    call mem_dealloc(FmatV)
    call mem_dealloc(eigvalsO)
    call mem_dealloc(eigvalsV)


  end subroutine ComputeFock_eigvecs


END MODULE davidson_solv_mod
