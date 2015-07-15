!> @file 
!> Solver for response equations (linear equations or eigenvalue problem).

module RSPsolver
  use matrix_util
  use rspPrecond
  use rsp_util, only: util_scriptPx, get_rsp_trials_from_MO, MO_precond
  use TYPEDEFTYPE, only: LSSETTING !due to rsp_molcfg type def.
  use response_wrapper_type_module, only: RSPSOLVERinputitem
  use lstiming
  use matrix_operations_aux, only: mat_zerohalf
  use decompMod
  use dal_interface
  use direct_dens_util
  use direct_dens_util_unres
  use memory_handling
  use II_XC_interfaceModule
  private
  public :: rsp_init, rsp_solver, make_lintran_vecs,&
       &get_1st_orth_trial_lineq, transform_vectors, expand_on_basis, &
       &expand_on_basis_minus, precond_residual, verify_conv_prop,&
       & remove_initial_lindep, rsp_normalize, symm_orthonormalize, get_rho,&
       & make_rhos, rsp_molcfg, get_1st_rsp_trials_unres, get_1st_rsp_trials,&
       & init_rsp_molcfg

  logical, save :: LINEQ
  logical, save :: RSPOnMaster
  integer, save :: rsp_number_of_rhs, rsp_number_of_sols, &
                 & rsp_number_of_omegas, lusigma_rsp, lurho_rsp, lub_rsp, &
                 & lub_mcdvec, nmcdvec
  real(realk), allocatable, save :: red_E(:,:), red_S(:,:), red_GD(:,:), red_GDI(:,:)
  real(realk), allocatable, save :: b_overlaps(:,:),bS_overlaps(:,:)  !Mostly for debug

   !> Molecule configuration type, abstracting data and settings
   !> to be passed to solver and integral routines.
   !> Should eventually be moved to separate program-specific
   !> interface modules
   type rsp_molcfg
      !> basic/prototype zero matrix, such as overlap or diplen
      !> Must have shape set, but may (should) not be allocated.
      !> Used to create/initialize other matrices, and thus hide
      !> program-specific information, like sparse/full,
      !> 1-comp/2-comp/4-comp, real/complex, etc. from response
      !> code.
      type(matrix)              :: zeromat
      !> number of atoms
      integer,          pointer :: natoms
      !> unit number for printing output
      integer,          pointer :: lupri
      !> unit number for printing errors
      integer,          pointer :: luerr
      !> integral program settings
      type(LSSETTING),  pointer :: setting
      !> decomposition settings and precomputed data
      type(decompItem), pointer :: decomp
      !> RSPsolver settings
      type(RSPSOLVERinputitem),pointer::solver
   end type

contains

  !> \brief Initializations of rsp_molcfg
  !> \author T. Kjaergaard
  !> \date 2011
  !> 
  subroutine init_rsp_molcfg(molcfg,inputmatrix,natoms,lupri,luerr,setting,decomp,solver)
    implicit none
    type(rsp_molcfg) :: molcfg
    !> number of atoms
    type(matrix)             :: inputmatrix
    integer,          target :: natoms
    integer,          target :: lupri
    integer,          target :: luerr
    type(LSSETTING),  target :: setting
    type(decompItem), target :: decomp
    type(RSPSOLVERinputitem),target::solver

    molcfg%zeromat%nrow = inputmatrix%nrow
    molcfg%zeromat%ncol = inputmatrix%ncol
    call mat_nullify(molcfg%zeromat)
    molcfg%zeromat%nnz  = inputmatrix%nnz
    molcfg%zeromat%complex = inputmatrix%complex
    molcfg%zeromat%init_magic_tag = inputmatrix%init_magic_tag

    molcfg%natoms => natoms
    molcfg%lupri => lupri
    molcfg%luerr => luerr
    molcfg%setting => setting
    molcfg%decomp => decomp
    molcfg%solver => solver

  end subroutine init_rsp_molcfg

  !> \brief Initializations for solver
  !> \author S. Host, S. Coriani
  !> \date 2006
  !>
  subroutine rsp_init(ntrial, nrhs, nsol, nomega, nstart)
  implicit none
      !> (Max.) number of trial vectors in a given iteration
      integer, intent(in) :: ntrial
      !> Number of right-hand sides. Only relevant for linear equations (always 1 for eigenvalue problem)
      integer, intent(in) :: nrhs
      !> Number of solution (output) vectors
      integer, intent(in) :: nsol
      !> If LINEQ, number of laser freq.s (input). Otherwise number of excitation energies (output) 
      integer, intent(in) :: nomega
      !> Number of start vectors. Only relevant for eigenvalue problem
      integer, intent(in) :: nstart

      rsp_number_of_rhs = nrhs
      rsp_number_of_sols = nsol
      rsp_number_of_omegas = nomega
      RSPOnMaster=.TRUE.
  end subroutine rsp_init

!> \brief The linear scaling response solver
!> \author S. Host, S. Coriani
!> \date 2006
!> 
!> LINEQ == FALSE  DIRECT THE SOLUTION THE GENERALIZED RSP
!>                 EIGENVALUE PROBLEM
!>                  ( E(2)-w(I)S(2))X(I) = 0
!> LINEQ == TRUE  DIRECT THE SOLUTION THE GENERALIZED RSP
!>                LINEAR EQUATIONS
!>                  ( E(2)-w(I)S(2))X(I) - GD = 0
!>
  subroutine rsp_solver(molcfg,D,S,F,LINEQ_x,n_gd_or_exci,GD,EIVAL,eivecs,Xproject)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Density matrix
    type(Matrix), intent(in)    :: D
    !> Overlap matrix
    type(Matrix), intent(in)    :: S
    !> Fock/KS matrix
    type(Matrix), intent(in)    :: F
    !> If true, solve linear equations. Otherwise, solve eigenvalue problem
    logical, intent(in)         :: LINEQ_x
    !> Number of gradients or excitations 
    integer, intent(in)         :: n_gd_or_exci 
    !> Right hand sides = property gradients
    type(Matrix), intent(inout) :: gd(rsp_number_of_rhs)
    !if LINEQ_x=true, laser frequencies (input). Otherwise, excitation energies (output)
    real(realk), intent(inout)  :: eival(rsp_number_of_omegas)  
    !> Output solution vectors
    type(Matrix), intent(inout) :: eivecs(rsp_number_of_sols)
!local
    type(Matrix),optional,intent(inout) :: Xproject(:) 
    type(Matrix),pointer :: Bvecs(:), rhos(:),sigmas(:)
!    type(Matrix) :: Bvecs(rsp_bvec_dim), rhos(rsp_bvec_dim), &
!                  & sigmas(rsp_bvec_dim)
    real(realk), allocatable :: red_X(:,:)
    integer :: itmic,ndim,i,ndim_red,Nb_new,l,j,max_it, max_red
    logical :: conv,exit_loop, make_rhos !make_rhos should control whether the S[2]b vectors are required
    type(Matrix) :: Sigma_scr, Rho_scr  !output solution vectors
    logical :: cov,OnMaster
    real(realk)  :: TSTR, TEN, TIMSTR,TIMEND, rsp_thresh
!for restart
    integer :: lurestart,nXproject
    logical :: UseExcitationVecs
    nullify(Bvecs)
    nullify(rhos)
    nullify(sigmas)
    max_it= molcfg%solver%rsp_maxit
    max_red= molcfg%solver%rsp_maxred
    rsp_thresh= molcfg%solver%rsp_thresh

    !the number of vectors Xproject to project out of the equation
    nmcdvec = 0
    if (present(Xproject)) nmcdvec = size(Xproject)
    if(nmcdvec.EQ.1)WRITE(molcfg%lupri,*)'The equation will be solved orthogonal to the ',nmcdvec,' singular component(s)'
    if(nmcdvec.EQ.2 .AND. (.NOT.molcfg%solver%degeneratestates))&
     &CALL LSQUIT('The solver is projecting against 2 vectors for a nondegenerate molecule?',molcfg%lupri)
    if(nmcdvec.GT.1)THEN
       UseExcitationVecs = molcfg%solver%UseExcitationVecs
       molcfg%solver%UseExcitationVecs = .FALSE.
    endif

    ndim = S%nrow
    call mat_init(Sigma_scr,ndim,ndim)
    call mat_init(Rho_scr,ndim,ndim)
    LINEQ = LINEQ_x

    !The reduced E2, S2 and gradient are global variables
    allocate(red_E(2*max_red,2*max_red), red_S(2*max_red,2*max_red))
    allocate(red_X(2*max_red,n_gd_or_exci))
    allocate(b_overlaps(2*max_red,2*max_red)) !Mostly for debug
    allocate(bS_overlaps(2*max_red,2*max_red)) !Mostly for debug
    IF (LINEQ) then
       !allocate space for reduced gradient
       ALLOCATE(RED_GD(2*max_red,molcfg%solver%rsp_maxgd))
       red_GD = 0.0E0_realk
    endif
    red_E = 0.0E0_realk ; red_S = 0.0E0_realk ; red_X = 0.0E0_realk
    if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'RSP_SOLVER: LINEQ, LINEQ_x',  LINEQ, LINEQ_x
    IF (LINEQ) THEN
       WRITE (molcfg%lupri,'(///2A//A/)') &
  &     ' <<<  SOLVING SETS OF LINEAR EQUATIONS ', &
  &     'FOR LINEAR RESPONSE PROPERTIES >>>'
    ELSE
       WRITE (molcfg%lupri,'(///2A//A/)') &
  &     ' <<< EXCITATION ENERGIES', &
  &     ' AND TRANSITION MOMENT CALCULATION (LSTDHF) >>>'
    END IF
    lusigma_rsp = -1 ; lurho_rsp = -1 ; lub_rsp = -1
    CALL LSOPEN(lusigma_rsp,'sigmarsp','unknown','UNFORMATTED')
    CALL LSOPEN(lurho_rsp,'rhorsp','unknown','UNFORMATTED')
    CALL LSOPEN(lub_rsp,'brsp','unknown','UNFORMATTED')


!    allocate(rhos(rsp_bvec_dim))
!    allocate(sigmas(rsp_bvec_dim))
!    allocate(bvecs(rsp_bvec_dim))
!    do i = 1, rsp_bvec_dim
!      call mat_init(rhos(i),ndim,ndim)
!      call mat_init(sigmas(i),ndim,ndim)
!      call mat_init(Bvecs(i),ndim,ndim)
!    enddo

    !Find initial trial vectors
    !Nb_new to the number of acceptable start trials 

    ! if MCD the exication vectors are the first trial vector 
    if (nmcdvec /= 0) then
       !AndreasJT: Test redundant, as nmcdvec is set from size(Xproject) above
       !if(.NOT.present(Xproject))then
       !   CALL LSQUIT('when projecting in rsp_solver the Xproject must be given as input',molcfg%lupri)
       !endif
       lub_mcdvec=-1
       CALL LSOPEN(lub_mcdvec,'mcd_rsp_vecs','NEW','UNFORMATTED')
       if (molcfg%solver%info_rsp) WRITE(molcfg%lupri,*)'RSP_SOLVER -> MCD_get_1st_trials '
       CALL MCD_get_1st_vec1(molcfg,F,D,S,n_gd_or_exci,gd,eival,bvecs,Nb_new,&
       &                     lub_mcdvec,Xproject,nmcdvec)
       Nb_new = nmcdvec
       Ndim_red = 0
       make_rhos = .true.
       if (ndim>2) then
         call transform_vectors(molcfg,D,S,F,Nb_new,bvecs,sigmas,rhos,make_rhos)
         call extend_red_matrices(molcfg,ndim_red,nb_new,n_gd_or_exci,gd,sigmas,rhos,bvecs)
         do i=1,size(sigmas)
            call mat_free(bvecs(i))
            call mat_free(sigmas(i))
            call mat_free(rhos(i))
         enddo
         call mem_dealloc(bvecs)
         call mem_dealloc(sigmas)
         call mem_dealloc(rhos)
         call MCD_get_1st_vec2(molcfg,F,D,S,n_gd_or_exci,gd,EIVAL,bvecs,Nb_new,conv)
         IF(conv)then
            !already done due to small dimensions of space.
            do i = 1,n_gd_or_exci
               call mat_zero(eivecs(i))               
               call expand_on_basis(molcfg,ndim_red,RED_X(1:ndim_red*2,i),lub_rsp,eivecs(i))
            enddo
            do i = 1,size(Bvecs)
!               call mat_free(Bvecs(i))
               call mat_free(Bvecs(i))
            enddo
            call mem_dealloc(Bvecs)
            IF(associated(rhos))THEN
               do i = 1,size(rhos)
                  call mat_free(rhos(i))
                  call mat_free(sigmas(i))
               enddo
               call mem_dealloc(rhos)
               call mem_dealloc(sigmas)
            ENDIF
            call mat_free(Sigma_scr)
            call mat_free(Rho_scr)
            DEALLOCATE(RED_X,RED_E,RED_S,b_overlaps,bS_overlaps)
            if (LINEQ) then 
               DEALLOCATE(RED_GD)
            endif
            CALL LSCLOSE(lusigma_rsp,'DELETE')
            CALL LSCLOSE(lurho_rsp,'DELETE')
            CALL LSCLOSE(lub_rsp,'DELETE')
            if(nmcdvec /= 0) then
               CALL LSCLOSE(lub_mcdvec,'DELETE')
            endif
            return
         endif
       endif
    else
       if (molcfg%solver%info_rsp) WRITE(molcfg%lupri,*)'RSP_SOLVER -> get_1st_trials '
       call get_1st_orth_trials(molcfg,D,S,n_gd_or_exci,gd,EIVAL,Nb_new,bvecs)
       Ndim_red = 0
    endif

    !loop until convergence
    do itmic = 1,max_it
       if (ndim_red > max_red) then
          WRITE(molcfg%lupri,'(/A)') &
          &     'rsp_maxred too small - recompile with larger rsp_maxred parameter'
          CALL LSQUIT('rsp_maxred parameter too small',molcfg%lupri)
       endif

       if (.not. molcfg%solver%rsp_quiet) CALL LSTIMER('START ',TIMSTR,TIMEND,molcfg%lupri)

       if(.NOT. molcfg%solver%rsp_quiet) write(molcfg%lupri,*) ' ** RSP_SOLVER MICROITERATION NUMBER',ITMIC

       !find sigma=E[2]bvec and rho=S[2]bvec of trial vectors
       if (molcfg%solver%info_rsp) WRITE(molcfg%lupri,*)'RSP_SOLVER -> transform_vectors, ITMIC: ', itmic
       if(.not. molcfg%solver%rsp_quiet) CALL LSTIMER('START ',TSTR,TEN,molcfg%lupri)
       make_rhos = .true.
       call transform_vectors(molcfg,D,S,F,nb_new,bvecs,sigmas,rhos,make_rhos)
       if(.not. molcfg%solver%rsp_quiet) CALL LSTIMER('LINTRA',TSTR,TEN,molcfg%lupri)

       !construct the problem in the reduced space and solve it
       if (molcfg%solver%info_rsp) WRITE(molcfg%LUPRI,*)'RSP_SOLVER -> build_red_space_solve, ITMIC: ', itmic
       call build_reduced_space_and_solve(molcfg,ndim_red,nb_new,n_gd_or_exci,&
                                        & gd,sigmas,rhos,bvecs,eival,RED_X)
       !compute residuals and check for converged vectors
       !make new trials and orthogonalize them versus the others
       if (molcfg%solver%info_rsp) WRITE(molcfg%LUPRI,*)'RSP_SOLVER -> get_residuals_and_new_b, ITMIC: ', itmic
       call get_residuals_and_new_b(molcfg,D,S,itmic,ndim_red,n_gd_or_exci, &
           & gd,RED_X,eival,Nb_new,bvecs,conv,exit_loop) 

       if (conv) then
          IF(associated(bvecs))then
             do i=1,size(bvecs)
                call mat_free(bvecs(i))
             enddo
             call mem_dealloc(bvecs)
          ENDIF
          do i = 1,n_gd_or_exci
             if (i > molcfg%solver%rsp_maxgd) then
                WRITE(molcfg%LUPRI,'(/A)') &
                &     'molcfg%solver%rsp_maxgd too small - recompile with larger molcfg%solver%rsp_maxgd parameter'
                CALL LSQUIT('molcfg%solver%rsp_maxgd parameter too small',molcfg%lupri)
             endif
             call mat_zero(eivecs(i))

             !Construct final solution vector by expanding red space solution onto trial vector space 
             if (molcfg%solver%info_rsp) WRITE(molcfg%LUPRI,*)'RSP_SOLVER -> expand_on_basis (final s), ITMIC: ', itmic
             call expand_on_basis(molcfg,ndim_red,RED_X(1:ndim_red*2,i),lub_rsp,eivecs(i))
             !WRITE(molcfg%LUPRI,*)'RSP_SOLVER -> Final solution i''th', i
             !call mat_print(eivecs(i),1,ndim,1,ndim,molcfg%lupri)
          enddo
          IF(.NOT.LINEQ_x)then
             lurestart = -1
             CALL LSOPEN(lurestart,'rsp_eigenvecs','unknown','UNFORMATTED')
             rewind(lurestart)             
             WRITE(molcfg%lupri,*)'Write ',n_gd_or_exci,' excitation vectors to disk'
             do i = 1,n_gd_or_exci
                write(lurestart) i, eival(i)
                OnMaster=.TRUE.
                call mat_write_to_disk(lurestart,eivecs(i),OnMaster)
             enddo
             molcfg%solver%rsp_eigenvecs = n_gd_or_exci
             CALL LSCLOSE(lurestart,'KEEP')
          endif
          exit
       else if (exit_loop) then
          do i = 1,n_gd_or_exci
             if (i > molcfg%solver%rsp_maxgd) then
                WRITE(molcfg%LUPRI,'(/A)') &
                &     'molcfg%solver%rsp_maxgd too small - recompile with larger molcfg%solver%rsp_maxgd parameter'
                CALL LSQUIT('molcfg%solver%rsp_maxgd parameter too small',molcfg%lupri)
             endif

             call mat_zero(eivecs(i))

             !Construct final solution vector by expanding red space solution onto trial vector space
             if (molcfg%solver%info_rsp) WRITE(molcfg%LUPRI,*)'RSP_SOLVER -> expand_on_basis (exit_loop), ITMIC: ', itmic
             call expand_on_basis(molcfg,ndim_red,RED_X(1:ndim_red*2,i),lub_rsp,eivecs(i))
          enddo
          exit
      endif
      if(.not. molcfg%solver%rsp_quiet) CALL LSTIMER('RSP_IT',TIMSTR,TIMEND,molcfg%lupri)
    enddo
!    do i = 1,rsp_bvec_dim
!       call mat_free(Bvecs(i))
!       call mat_free(rhos(i))
!       call mat_free(sigmas(i))
!    enddo
!    deallocate(rhos)
!    deallocate(sigmas)
!    deallocate(bvecs)
    call mat_free(Sigma_scr)
    call mat_free(Rho_scr)
    DEALLOCATE(RED_X,RED_E,RED_S,b_overlaps,bS_overlaps)
    if (LINEQ) then 
       DEALLOCATE(RED_GD)
    endif

    CALL LSCLOSE(lusigma_rsp,'DELETE')
    CALL LSCLOSE(lurho_rsp,'DELETE')
    CALL LSCLOSE(lub_rsp,'DELETE')
    if(nmcdvec .GT. 0) then
       CALL LSCLOSE(lub_mcdvec,'DELETE')
       if(nmcdvec .GT. 1)molcfg%solver%UseExcitationVecs = UseExcitationVecs
    endif

    IF(.NOT.LINEQ_x)THEN
       call orthogonalizeDegenerate(molcfg,D,S,F,eivecs,eival,n_gd_or_exci)
    ENDIF

  end subroutine rsp_solver

!> \brief Wrapper for setting up initial trial vectors
!> \author S. Host, S. Coriani
!> \date 2006
  subroutine get_1st_orth_trials(molcfg,D,S,n_gd_or_exci,gd,eival,Nb_new,bvecs)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Density matrix
    type(Matrix), intent(in)   :: D
    !> Overlap matrix
    type(Matrix), intent(in)   :: S
    !> Number of gradients or excitations 
    integer, intent(in)        :: n_gd_or_exci
    !> Right hand sides = property gradients
    type(Matrix), intent(in)   :: gd(rsp_number_of_rhs)
    !if LINEQ=true, laser frequencies (input). Otherwise, excitation energies (output)
    real(realk), intent(inout) :: eival(rsp_number_of_omegas)
    !> Number of generated trial vectors
    integer, intent(out)       :: Nb_new
    !> The generated trial vectors (output)
    type(Matrix),pointer :: bvecs(:)

    IF (LINEQ) THEN
       if (molcfg%solver%info_rsp) WRITE(molcfg%LUPRI,*) 'Now enter 1st guess for linear equations '
       CALL get_1st_orth_trial_lineq(molcfg,D,S,n_gd_or_exci,gd,eival,bvecs,Nb_new)
    ELSE
       if (molcfg%solver%info_rsp) WRITE(molcfg%LUPRI,*) 'Now enter 1st guess for excitation energies'
       call get_1st_orth_trial_eigen(molcfg,D,S,n_gd_or_exci,bvecs,Nb_new)
    END IF

  end subroutine get_1st_orth_trials


!> \brief Set up initial trial vectors (linear equations)
!> \author S. Coriani, P. Jorgensen
!> \date June 2003
!> 
!>  Generate start vector(s) for solution of a linear
!>  set of equations from the preconditioned property
!>  gradient vector(s) 
!>  The new, linearly independent orthogonalized vectors
!>  to be used as trials are returned in BVECS
!> 
!>  Modified version of LRST for Linear Scaling. 
!> 
  subroutine get_1st_orth_trial_lineq(molcfg,D,S,ngd,gd,eival,bvecs,Nb_new)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Density matrix
    type(Matrix),intent(in)    :: D
    !> Overlap matrix
    type(Matrix),intent(in)    :: S
    !> Number of gradients
    integer,intent(in)         :: ngd
    !> Right hand sides = property gradients
    type(Matrix),intent(in)    :: gd(rsp_number_of_rhs)
    !> Laser frequencies (input)
    real(realk), intent(inout) :: eival(rsp_number_of_omegas)
    !> The generated trial vectors (output)
    type(Matrix),pointer :: bvecs(:)
    !> Number of generated trial vectors
    integer, intent(out)       :: Nb_new
    type(Matrix) :: grad_x,Pgrad
    type(matrix),pointer :: Bvec_tmp(:)
    real(realk) :: norm,eivaldummy
    integer :: ndim,nb_prev,ibx,i,lurestart,k,nextra,ibx2
    logical :: fileexists,OnMaster,UseExcitationVecs
    UseExcitationVecs = molcfg%solver%UseExcitationVecs
    ndim = S%nrow
    IF(UseExcitationVecs)THEN
       INQUIRE(file='rsp_eigenvecs',EXIST=fileexists)
       if (fileexists) then
          nextra = molcfg%solver%rsp_eigenvecs
       else
          nextra = 0
       endif
    ELSE
       nextra = 0
    ENDIF    
    call mem_alloc(Bvec_tmp,ngd+nextra)
!
! Generate preconditioned gradients as start vectors
!
    if (ngd > rsp_number_of_rhs) then
         WRITE(molcfg%LUPRI,'(/A)') &
         &     'ngd > rsp_number_of_rhs'
         CALL LSQUIT('ngd > rsp_number_of_rhs',molcfg%lupri)
    endif
    if (ngd > rsp_number_of_omegas) then
         WRITE(molcfg%LUPRI,'(/A)') &
         &     'ngd > rsp_number_of_omegas'
         CALL LSQUIT('ngd > rsp_number_of_omegas',molcfg%lupri)
    endif
    ibx2 = 0
    Nb_new = 0
    call mat_init(grad_x,ndim,ndim)
    do ibx = 1, ngd
        call mat_assign(grad_x,gd(ibx))  !STINNE - we do now not divide in symm and antisymm parts

       !Sonia_mag
       !WRITE(molcfg%LUPRI,*) '-----------------------------------------------'
       !WRITE(molcfg%LUPRI,*) 'The first trial of LINEQ before preconditioning'
       !call mat_print(grad_x,1,grad_x%nrow,1,grad_x%nrow,molcfg%lupri)
       !WRITE(molcfg%LUPRI,*) '-----------------------------------------------'
       norm = mat_sqnorm2(grad_x)                  
       if (norm > molcfg%solver%rsp_tolerance) then          
          !Precondition if the matrix exists  by
          !solving A Pgrad = grad
          !print *, 'precond 1st trial'
          call mat_init(Pgrad,ndim,ndim)
          call rsp_AB_precond(molcfg,grad_x,S,EIVAL(ibx),Pgrad) !FIXME: needs to be changed to work for more than one frequency at a time
          ibx2 = ibx2+1
          call mat_init(Bvec_tmp(ibx2),ndim,ndim)
          call mat_assign(Bvec_tmp(ibx2),Pgrad)
          call mat_free(Pgrad)
          !Sonia_mag
          !WRITE(molcfg%LUPRI,*) '-----------------------------------------------'
          !WRITE(molcfg%LUPRI,*) 'The first trial of LINEQ after preconditioning'
          !call mat_print(Pgrad,1,grad_x%nrow,1,grad_x%nrow,molcfg%lupri)
          !WRITE(molcfg%LUPRI,*) '-----------------------------------------------'
          !Remove preconditioning!!!
          !Bvec_tmp(ibx) = gd(ibx)
          !end Remove preconditioning!!!
          Nb_new = Nb_new + 1
       ELSE
          ibx2 = ibx2+1
          call mat_init(Bvec_tmp(ibx2),ndim,ndim)
          call mat_assign(Bvec_tmp(ibx2),grad_x)
          Nb_new = Nb_new + 1
       endif
    enddo
    call mat_free(grad_x)

!    WRITE(molcfg%lupri,*)'Nb_new',Nb_new
    Nb_prev = 0

    IF(UseExcitationVecs)THEN
       if (fileexists) then
          WRITE(molcfg%lupri,'(A,i5,A)')'Use the ',molcfg%solver%rsp_eigenvecs,' excitation vectors from the rsp_eigenvecs file.'
          WRITE(molcfg%lupri,'(A,i5,A)')'ngd = ',ngd
          OnMaster = .TRUE.
          lurestart = -1
          CALL LSOPEN(lurestart,'rsp_eigenvecs','old','UNFORMATTED')
          rewind(lurestart)
          do i = 1, molcfg%solver%rsp_eigenvecs
             read(lurestart) k, eivaldummy
             call mat_init(Bvec_tmp(Nb_new+i),ndim,ndim)
             call mat_read_from_disk(lurestart,Bvec_tmp(Nb_new+i),OnMaster)
          enddo
          call LSCLOSE(lurestart,'KEEP') 
          Nb_new = Nb_new + molcfg%solver%rsp_eigenvecs
          ibx2 = Nb_new
       endif
    ENDIF

!
! orthogonalize and remove linear dependent vectors
! Nb_new is updated inside orthonormalize
!
    if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Nb_new to orthonormalize:',Nb_new
    if (Nb_new == 0) then
       WRITE (molcfg%LUPRI,'(//A,I5)') &
                    &  'get_1st_orth_trial_lineq: No start vectors present!'
       CALL LSQUIT('get_1st_orth_trial_lineq: No start vectors present!',molcfg%lupri)
    endif
    !
    call orthonormalize(molcfg,D,S,Bvec_tmp,Nb_new,Nb_prev,bvecs)
!
! Number of start trial vectors = Nb_new 
!
    if (Nb_new <= 0) then
       WRITE (molcfg%LUPRI,'(//A,I5)') &
                    &  'get_1st_orth_trial_lineq: START VECTOR IS NOT LINEAR INDEPENDENT.'
      CALL LSQUIT('START VECTOR IS NOT LINEAR INDEPENDENT',molcfg%lupri)
    endif

!FREE USED MATRICES!!
    do i = 1,ibx2
       call mat_free(Bvec_tmp(i))
    enddo
    call mem_dealloc(Bvec_tmp)
  end subroutine get_1st_orth_trial_lineq

!> \brief Set up initial trial vectors (eigenvalue problem)
!> \author L. Thogersen, S. Coriani, S. Host
!> \date November 2005. Modified October 2006
!> 
!>  Based on old PPST routine.
!>  Generate start vector(s) for solution of generalized eigenvalue equations via either:
!>  - backtranforming the MO unit vectors corresponding to the resorted
!>             roots E[2]_ii/|S[2]_ii| into AO
!>  or
!>  - create start vectors based on the solution vectors from iterative
!>    determination of the orbital energies
!>
  subroutine get_1st_orth_trial_eigen(molcfg,D,S,nexcit,bvecs,Nb_new)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Density matrix
    type(Matrix),intent(in)    :: D
    !> Overlap matrix
    type(Matrix),intent(in)    :: S
    !> Number of eigenvalues/roots (required in input)
    integer,intent(in)         :: nexcit
    !> The generated trial vectors (output)
    type(Matrix),pointer :: bvecs(:)
    !> Number of generated trial vectors
    integer, intent(out)       :: Nb_new
    type(Matrix),allocatable   :: Bvec_tmp(:) !it is possible to have more start vecs than nexcit
    real(realk)                :: norm
    integer                    :: ndim,nb_prev,ibx,i,howmany
!RESTART variables:
    logical                    :: fileexists,OnMaster
    integer                    :: j, k, lurestart
    real(realk)                :: eival
    OnMaster = .TRUE.
    ndim = S%nrow

    if (molcfg%solver%rsp_startvectors) then 
       howmany = molcfg%solver%rsp_no_of_startvectors
    else
       howmany = nexcit
    endif

    allocate(Bvec_tmp(howmany))
    do i = 1,howmany
       call mat_init(Bvec_tmp(i),ndim,ndim)
    enddo
!
! Generate start vectors
!
    if (molcfg%solver%cfg_unres) then
       call get_1st_rsp_trials_unres(molcfg,howmany,Bvec_tmp)
    else
       call get_1st_rsp_trials(molcfg,howmany,Bvec_tmp)
    endif
!
! Stinne May 2008: Possibility to restart exci calculation from file rsp_eigenvecs!
!
    if (molcfg%solver%rsp_restart_exci) then
       INQUIRE(file='rsp_eigenvecs',EXIST=fileexists)
       if (.not. fileexists) then
          WRITE(molcfg%LUPRI,'(/A)') &
               &     'File "rsp_eigenvecs" does not exist! Must be present with .RESTART'
          CALL LSQUIT('File "rsp_eigenvecs" does not exist!',molcfg%lupri)
       endif
       WRITE(molcfg%lupri,'(A,i5,A)')'Restart from ',molcfg%solver%rsp_restart_nexci,' vectors from the rsp_eigenvecs file.'
       lurestart = -1
       CALL LSOPEN(lurestart,'rsp_eigenvecs','old','UNFORMATTED')
       rewind(lurestart)
       do i = 1, molcfg%solver%rsp_restart_nexci
          read(lurestart) k, eival
          call mat_read_from_disk(lurestart,Bvec_tmp(i),OnMaster)
       enddo
       print*,'REMOVE solutions to disk for restart molcfg%solver%rsp_restart_exci'
       call LSCLOSE(lurestart,'DELETE') !Delete so it doesn't conflict with creation of
    endif                               !new "rsp_eigenvecs" after the calculation is done
!
! orthogonalize and remove linear dependent vectors
! Nb_new is updated inside orthonormalize
!
    Nb_prev = 0
    Nb_new = howmany
    if(molcfg%solver%info_rsp) write(molcfg%lupri,*) 'GET_1ST_TRIALS_EIGEN: Nb_new to orthonormalize:',Nb_new
    if (Nb_new == 0) then
       WRITE (molcfg%LUPRI,'(//A,I5)') &
                    &  'get_1st_orth_trial_lineq: No start vectors present!'
       CALL LSQUIT('get_1st_orth_trial_lineq: No start vectors present!',molcfg%lupri)
    endif

    if (molcfg%solver%info_rsp) then
       do i = 1, Nb_new
          norm = mat_sqnorm2(Bvec_tmp(i))
          if(molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Before orthonormalize |btm| = ', i , norm
       end do
    end if
    call orthonormalize(molcfg,D,S,Bvec_tmp,Nb_new,Nb_prev,Bvecs)
!
! Number of start trial vectors = Nb_new  
!
    if (Nb_new <= 0) then
       WRITE (molcfg%LUPRI,'(//A,I5)') &
&       'get_1st_orth_trial_eigen: START VECTOR IS NOT LINEAR INDEPENDENT.'
      CALL LSQUIT('START VECTOR IS NOT LINEAR INDEPENDENT',molcfg%lupri)
    endif
    if (Nb_new < howmany) then
       !if(molcfg%solver%info_rsp)write(molcfg%lupri,*) 'Inside get_1st_orth_trial_eigen'
       if(molcfg%solver%info_rsp)write(molcfg%lupri,*) '**WARNING: Less start vectors than required excitations!'
    end if
!
! Free memory
!
    do i = 1,howmany
       call mat_free(Bvec_tmp(i))
    enddo
    deallocate(Bvec_tmp)

  end subroutine get_1st_orth_trial_eigen

!> \brief Orthogonalize new b-vectors against all previous b-vectors and among themselves, and renormalize
!> \author L. Thogersen, S. Coriani
!> \date November 2005
!> 
!> Based on old subroutine RSPORT
!> Purpose:
!>  Orthogonalize new b-vectors against all previous b-vectors
!>  and among themselves, and renormalize.
!>  The b-vectors have the form
!>        ( Y_dia    Z   )       ! Z_mu_nu  mu < nu ; Y_dia = Y_mu_mu 
!>        (  Y     Y_dia )       ! Y_mu_nu  mu > nu    
!>  Each b-vector is in (Z, Y) form, the (Y, Z) vector is obtained
!>  by transposing.  Each of the new b-vectors is
!>  orthogonalized both against (Z, Y) and (Y, Z) for each previous
!>  b-vectors.
!>  (Orthogonalization is performed twice if round-off is large,
!>   see normalize)
!> 
!> In input:
!>  previous vectors are on disk in lub_rsp (global variable, see top of file) 
!> 
  subroutine orthonormalize(molcfg,D,S,Bvec_tmp,Nb_new,Nb_prev,bvecs)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Density matrix
    type(Matrix),intent(in) :: D
    !> Overlap matrix
    type(Matrix),intent(in) :: S
    !> New non-orthogonal vectors
    type(Matrix),intent(inout) :: Bvec_tmp(:)
    !> Input: Number of new non-orthogonal b-vectors in Bvec_tmp. Output: Number of acceptable new trial vectors after orthonormalization
    integer, intent(inout)  :: Nb_new
    !> Number of previous b-vectors (on disk in lub_rsp)
    integer, intent(in)     :: Nb_prev
    !> New orthogonalized vectors (output)
    type(Matrix),pointer :: bvecs(:)
!local 
    integer :: irx,i,iturn,ndim,k,ibvec,jbvec,jrx,no_of_new,lub_rsp_vec,dummy_i,max_red
    integer,allocatable :: lin_depend(:) !linear dependency index array
    type(matrix) :: B_scr, b_k, Xf, rho_k
    type(matrix) :: orthovec !Will properly be changed
    real(realk) :: TT,T1,T2,dummy_real
    logical :: run_ortho_again

    max_red= molcfg%solver%rsp_maxred
    if ((Nb_prev + Nb_new) > max_red) then
      WRITE (molcfg%LUPRI,*) 'ERROR in orthonormalize: NB_PREV + NB_NEW > rsp_maxred'
      WRITE (molcfg%LUPRI,*) 'NB_PREV ,NB_NEW  =',NB_PREV,NB_NEW
      WRITE (molcfg%LUPRI,*) 'rsp_maxred =',max_red
      WRITE (molcfg%LUPRI,*) 'Reduce problem size or recompile orthonormalize '// &
      &                  'with larger rsp_maxred parameter'
      WRITE (*,*) 'ERROR in orthonormalize: NB_PREV + NB_NEW > rsp_maxred'
      WRITE (*,*) 'NB_PREV ,NB_NEW  =',NB_PREV,NB_NEW
      WRITE (*,*) 'rsp_maxred =',max_red
      WRITE (*,*) 'Reduce problem size or recompile orthonormalize '// &
      &                  'with larger rsp_maxred parameter'
      CALL LSQUIT('std orthonormalize error: NB_PREV + NB_NEW > rsp_maxred parameter',molcfg%lupri)
    END IF
    
    ndim = S%nrow
    ALLOCATE(lin_depend(Nb_new))
    call mat_init(B_scr,ndim,ndim)
    call mat_init(b_k,ndim,ndim)
    if(nmcdvec /= 0) then
      call mat_init(Xf,ndim,ndim)
      call mat_init(rho_k,ndim,ndim)
    endif

    ! STEP 1: check initial linear dependencies among new vectors
    if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'orthonormalize -> remove_initial_lindep'
    lin_depend(1:Nb_new) = 1 !Initialize to non-dependent
    call remove_initial_lindep(molcfg,Nb_new,Bvec_tmp,B_scr) 

   iturn = 0
   run_ortho_again = .false.
   do   !It might be necessary to run through everything twice - decided in normalize
   iturn = iturn + 1
   !
   !Project out
   !
      if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'orthonormalize -> util_scriptPx'
      do i = 1,Nb_new
         if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Current vector is number = ', i
         call util_scriptPx('N',D,S,Bvec_tmp(i))
      enddo
      !
      ! Orthogonalize new b-vectors agains previous (old) b-vectors
      ! (both (Z, Y) and (Y, Z))
      !
      if (molcfg%solver%info_rsp) then
         write(molcfg%lupri,*) 'Orthogonalize new b-vectors agains previous b'
         write(molcfg%lupri,*) 'NB_PREV, NB_NEW ', NB_PREV, NB_NEW
      endif 
   
      if (nmcdvec /= 0) rewind(lub_mcdvec)
      rewind(lub_rsp)
      do k = 1,Nb_prev
         call mat_read_from_disk(lub_rsp,b_k,RSPOnMaster)
         if (nmcdvec /= 0 .and. k <= nmcdvec) then
            call mat_read_from_disk(lub_mcdvec,Xf,RSPOnMaster)
            call mat_read_from_disk(lub_mcdvec,rho_k,RSPOnMaster)
         endif
         do irx = 1,Nb_new
            !if lin_depend(irx) == 0, the vector is skipped
            !because of linear dependencies
            if (lin_depend(irx) /= 0) then
               if (nmcdvec /= 0 .and. k <= nmcdvec) then
                  TT = mat_dotproduct(rho_k,Bvec_tmp(irx))
                  call mat_daxpy(-TT,Xf,Bvec_tmp(irx))
                  call mat_trans(rho_k,b_scr)
                  TT = -mat_dotproduct(B_scr,Bvec_tmp(irx))
                  call mat_trans(Xf,b_scr)
                  call mat_daxpy(TT,b_scr,Bvec_tmp(irx))
               else
                  !Orthogonalize to (Z Y)_old
                  TT = mat_dotproduct(b_k,Bvec_tmp(irx))
                  call mat_daxpy(-TT,b_k,Bvec_tmp(irx))
                  !Orthogonalize to (Y Z)_old
                  !find the paired
                  call mat_trans(b_k,b_scr)
                  TT = mat_dotproduct(B_scr,Bvec_tmp(irx))
                  call mat_daxpy(-TT,B_scr,Bvec_tmp(irx))
               endif
            endif
         enddo
      enddo

      if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Orthogonalize new b-vectors against each other '
      !
      ! Orthogonalize new vectors against each other
      !
      do ibvec = 1,Nb_new !index for current bvector 
         if (molcfg%solver%info_rsp) write(molcfg%lupri,*)'Current ibvec index is ', ibvec, lin_depend(ibvec)
         if (lin_depend(ibvec) /= 0) then
            jbvec = 1    !index for another bvector
            do jrx = 1,(ibvec-1)
               if (lin_depend(jrx) == 0) then
                  jbvec = jbvec + 1 !skip
               else
                  T1 = mat_sqnorm2(Bvec_tmp(jbvec))
                  if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Norm of jbvec', jbvec, T1
                  T2 = mat_dotproduct(Bvec_tmp(jbvec),Bvec_tmp(ibvec))
                  if (molcfg%solver%info_rsp) write(molcfg%lupri,*) '<bjbvec|bibvec>', T2
                  TT = -T2/T1
                  if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'TT = -T2/T1', TT
   
                  call mat_daxpy(TT,Bvec_tmp(jbvec),Bvec_tmp(ibvec))
   
                  call mat_trans(Bvec_tmp(jbvec),B_scr)
                  TT = mat_dotproduct(B_scr,Bvec_tmp(ibvec))
                  TT = -TT/T1
   
                  call mat_daxpy(TT,B_scr,Bvec_tmp(ibvec))
   
                  jbvec = jbvec + 1
               endif
            enddo
         endif
         !
         ! NORMALIZE VECTOR NUMBER Ibvec
         !
         if (lin_depend(ibvec) /= 0) then !not linear dependent
            if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'RSP_SOLVER: Orthonormalize -> normalize'
            call rsp_normalize(molcfg,iturn,ibvec,run_ortho_again,lin_depend(ibvec),Bvec_tmp(ibvec))
         endif

         ! Perform symmetric orthonormalization of the pair ibvec, ibvec^T
         if (lin_depend(ibvec) /= 0) then !not linear dependent
            if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'RSP_SOLVER: Orthonormalize -> Symm_orthonormalize'
            call symm_orthonormalize(molcfg,Bvec_tmp(ibvec),ibvec,lin_depend(ibvec))
         endif 

      enddo !ibvec
   
      if (run_ortho_again) then
         run_ortho_again = .false.
         if (iturn == 2) then !This is redundant since run_ortho_again is only set .true. if iturn = 1
            WRITE(molcfg%LUPRI,'(/A)') &
            &     'Error: Already ran twice through orthonormalize!'
            CALL LSQUIT('Error: Already ran twice through orthonormalize!',molcfg%lupri)
         else
            cycle
         endif
      else
         exit
      endif
   enddo !Maybe do twice

   if (nmcdvec /= 0) then
      call mat_free(Xf)
      call mat_free(rho_k)
   endif

   ! Add new vectors to file 
   ! ib counts how many "acceptable" vectors there are in total
   !    
   no_of_new = 0
   do irx = 1,Nb_new
      if (lin_depend(irx) /= 0) then
         no_of_new = no_of_new + 1
      endif
   enddo
   call mem_alloc(Bvecs,no_of_new)
   no_of_new = 0
   do irx = 1,Nb_new
      if (lin_depend(irx) /= 0) then
         no_of_new = no_of_new + 1
         call mat_init(bvecs(no_of_new),Bvec_tmp(irx)%nrow,Bvec_tmp(irx)%ncol)
         call mat_assign(bvecs(no_of_new),Bvec_tmp(irx))
      endif
   enddo
   !
   ! Set NB_NEW to actual number of acceptable new trial vectors
   !
   Nb_new = no_of_new
     
   DEALLOCATE(lin_depend)
   call mat_free(B_scr)
   call mat_free(b_k)

  end subroutine orthonormalize

!> \brief Remove linear dependence caused by symmetric or antisymmetric trial vectors
!> \author L. Thogersen, S. Coriani
!> \date November 2005
!>
!> For zero frequency linear response, the result will
!> be Z =-Y for real perturbations (antisymmetric matrix) and 
!> Z = Y for imaginary perturbations (symmetric matrix) 
!> ; trial vectors will have the same structure.
!>
!> if Z = -Y  (Y_dia=0) then the paired (Z Y) and (Y Z) are linear dep.
!> we set -Y = 0,  Y_dia is already zero 
!> if Z = Y then (Z Y) and (Y Z) are linear dependent
!> we set  Z = 0  and Y_dia = Y_dia/2 
!>
  subroutine remove_initial_lindep(molcfg,Nb_new,Bvec_tmp,B_scr)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Number of trial vectors in Bvec_tmp
    integer, intent(in) :: Nb_new
    !> The trial vectors
    type(Matrix), intent(inout) :: Bvec_tmp(:)
    !> Scratch space
    type(Matrix), intent(inout) :: B_scr
    integer :: irx,ndim
    real(realk) :: norm,norm_symm,norm_antisymm, lea
    logical :: sym_or_antisym

    sym_or_antisym= .false.
    ndim = B_scr%nrow
    do irx = 1,Nb_new    
       !norm^2 of current matrix
       norm = mat_sqnorm2(bvec_tmp(irx))
       !Find symmetric part of matrix
       call mat_assign(B_scr,Bvec_tmp(irx))
       call util_get_symm_part(B_scr)
       !Find norm of symmetric part
       !if norm = 0, matrix is antisymmetric
       norm_symm = mat_outdia_sqnorm2(B_scr)
       if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'remove_initial_lindep: Norm of symm part = ', norm_symm
       !test if Z + Y is zero
       if (norm_symm <= molcfg%solver%rsp_ovlmin*norm) then
          !antisymmetric matrix
          if (molcfg%solver%info_rsp) then
             WRITE(molcfg%LUPRI,*) '                                       '
             WRITE(molcfg%LUPRI,*) ' Z = -Y in trial vector no.',IRX
             WRITE(molcfg%LUPRI,*) ' Y (or Z?) component removed and Ydia=Ydia/2.'
          endif
          sym_or_antisym = .true.
       else
          !Find antisymmetric part of matrix
          call util_get_antisymm_part(Bvec_tmp(irx),B_scr)
          !Find norm of antisymmmetric part
          !if norm = 0, matrix is symmetric
          norm_antisymm = mat_sqnorm2(B_scr)
          if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'remove_initial_lindep: ', &
                  & '         Norm of asymm part = ', norm_antisymm
          if (norm_antisymm <= molcfg%solver%rsp_ovlmin*norm) then
             !symmetric matrix
             if (molcfg%solver%info_rsp) then
               WRITE(molcfg%LUPRI,*) '                                       '
               WRITE(molcfg%LUPRI,*) ' Z = Y in trial vector no.',IRX
               WRITE(molcfg%LUPRI,*) ' Z component removed and Y_dia= Y_dia/2.'
             endif
             sym_or_antisym = .true.
          endif
       endif
       if (sym_or_antisym) then
          !set the upper triangle to zero and divide the diagonal w. 2
          call mat_zerohalf('UT',Bvec_tmp(irx))
          call mat_scal_dia(0.5E0_realk,Bvec_tmp(irx))
          sym_or_antisym = .false.
       endif
    enddo
        
  end subroutine remove_initial_lindep

!> \brief Perform symmetric orthonormalization of trial vectors
!> \author L. Thogersen, S. Coriani
!> \date November 2005
!>
!> Perform symmetric orthonormalization of (Z Y) and (Y Z) pair
!> for vector number i
!>
!>   -1/2       ( C1   C2 )               (  1     OVLPI )
!>  S      =    (         )   where S  =  (              )
!>              ( C2   C1 )               ( OVLPI     1  )
!>------------------------------------------------------------------
!> - i.e. bvec should be orthonormal to its own transpose
!>
  subroutine symm_orthonormalize(molcfg,bvec_tmp_i,i,lin_depend_i)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Input: The i'th bvector to be symmetrically orthonormalized. Output: The i'th bvector, now symmetrically orthonormalized
    type(Matrix), intent(inout) :: Bvec_tmp_i
    !> Index of Bvec_temp_i among the new bvectors (used only for printout if linearly dependent)
    integer, intent(in) :: i 
    !> Linear dependence indicator: 1 if non-dependent, 0 if linearly dependent 
    integer, intent(inout) :: lin_depend_i
    type(Matrix) :: B_scr,B_scr2
    real(realk) :: ovlpi,x1,x2,c1,c2
    
    integer :: j,ndim

    ndim = Bvec_tmp_i%nrow
    call mat_init(B_scr,ndim,ndim)
    call mat_init(B_scr2,ndim,ndim)

    !do it twice for numerical stability
    do j = 1,2
       call mat_trans(Bvec_tmp_i,B_scr)
       ovlpi = mat_dotproduct(B_scr,Bvec_tmp_i) 
       x1 = mat_sqnorm2(Bvec_tmp_i)

       if (molcfg%solver%info_rsp) WRITE (molcfg%lupri,'(A,1P,D14.6,D19.10)') &
       &          ' symm_orth, norm and <ZY/YZ> overlap before S-1/2',X1,OVLPI

       X2 = ABS(OVLPI)

       IF (ABS(1E0_realk-X2) < molcfg%solver%rsp_OVLMIN) THEN
          write(molcfg%lupri,"('bvector number ', i3, ' is removed because of linear dependence')") i
          write(molcfg%lupri,"('since <Z Y|Y Z> overlap is ', d12.5, ' which is ~1')") X2
          lin_depend_i = 0
          exit
       END IF

       X1 = 1E0_realk+OVLPI
       X2 = 1E0_realk-OVLPI
       IF (ABS(X1 - X2) > molcfg%solver%rsp_thr_lin_depend) THEN
          X1 = 0.5E0_realk / SQRT(X1)
          X2 = 0.5E0_realk / SQRT(X2)
          C1 = X1 + X2
          C2 = X1 - X2
          call mat_assign(B_scr2,Bvec_tmp_i)
          call mat_add(c1,B_scr2,c2,B_scr,Bvec_tmp_i)
       endif
    enddo   
    call mat_free(B_Scr)
    call mat_free(B_scr2)

  end subroutine symm_orthonormalize

!> \brief Normalize the i'th trial vector
!> \author L. Thogersen, S. Coriani
!> \date November 2005
  subroutine rsp_normalize(molcfg,iturn,i,run_ortho_again,lin_depend_i,Bvec_tmp_i)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Which time we are running through rsp_orthonormalize, which calls normalize
    integer, intent(in) :: iturn
    !> Index of Bvec_temp_i among the new trial vectors (used for printout)
    integer, intent(in) :: i
    !> If linear dependencies are suspected, this is set to true, and orthonormalize is run again
    logical, intent(out) :: run_ortho_again
    !> Linear dependence indicator: 1 if non-dependent, set to 0 if linearly dependent
    integer, intent(inout) :: lin_depend_i
    !> Input: The i'th trial vector to be normalized. Output: The normalized i'th trial vector
    type(Matrix), intent(inout) :: Bvec_tmp_i
    real(realk) :: TT

    TT = mat_sqnorm2(Bvec_tmp_i)
    if (TT < molcfg%solver%rsp_thr_lin_depend) then
       write(molcfg%lupri,"('bvector number ', i3, ' is removed because its norm**2 :', d12.5, &
              & ' is < molcfg%solver%rsp_thr_lin_depend ', d12.5 )") i,TT,molcfg%solver%rsp_thr_lin_depend
       lin_depend_i = 0  !linear dependent
    else if (TT < molcfg%solver%rsp_thr_round) then
       if (iturn == 1) then
          run_ortho_again = .true.
       else
          write(molcfg%lupri,"('bvector number ', i3, ' is removed because its norm**2 :', d12.5, &
              & ' is < molcfg%solver%rsp_thr_round ', d12.5 ,' after 2nd Gram-Schmidt orth.')") i,TT,molcfg%solver%rsp_thr_round
          lin_depend_i = 0 !linear dependent
       endif
    else
       !I am not sure this is correct but this lin_depend_i
       !was uninitialised if non of the 2 cases were true
       lin_depend_i = 1  !not linear dependent
    endif
    !Vector is not linearly dependent - normalize it:
    IF (lin_depend_i /= 0) THEN
       IF (molcfg%solver%info_rsp) then
          write(molcfg%lupri,"('Now normalize bvector number ', i3, ' with initial norm**2 :', d12.5)") i,TT
       endif
       !If norm is very small, increase it by scaling bvector
       IF (TT < molcfg%solver%rsp_T1MIN) THEN
          TT = 1E0_realk / SQRT(TT)
          call mat_SCAL(TT,BVEC_tmp_i)
          TT = mat_sqnorm2(Bvec_tmp_i)
       END IF
       !Normalization
       TT = 1E0_realk / SQRT(TT)
       call mat_SCAL(TT,BVEC_tmp_i) 
    END IF

  end subroutine rsp_normalize

!> \brief Interface routine to carry out linear transformations
!> \author S. Coriani, T. Kjaergaard 
!> \date June 2003
!>  
!>   SIGMAS(I) = E[2]*N(I) (sigma)
!>   RHOS(I) = S[2]*N(I) (rho) > check omega of RPA
!>
  subroutine transform_vectors(molcfg,D,S,F,nnew,bvecs,sigmas,rhos,make_rhos)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Density matrix
    type(Matrix), intent(in) :: D
    !> Overlap matrix
    type(Matrix), intent(in) :: S
    !> Fock/KS matrix
    type(Matrix), intent(in) :: F
    !> Number of trial vectors in bvecs to be linearly transformed
    integer, intent(in) :: nnew
    !> Trial vectors in bvecs to be linearly transformed
    type(Matrix), intent(in) :: bvecs(:)
    !> Sigma part of linear transformations (output)
    type(Matrix), pointer :: sigmas(:)
    !> Rho part of linear transformations, if make_rhos = true (output)
    type(Matrix), pointer :: rhos(:)
    !> True if the rho part of linear transformation should be calculated
    logical, intent(in) :: make_rhos
    type(Matrix) :: scr
    integer :: i,ndim,l 
    ndim = S%nrow
!    if(make_rhos)then
       call mem_alloc(rhos,nnew)
!    endif
    call mem_alloc(sigmas,nnew)
    do i=1,nnew
!    if(make_rhos)then
       call mat_init(rhos(i),Bvecs(1)%nrow,Bvecs(1)%ncol)
!    endif
       call mat_init(sigmas(i),Bvecs(1)%nrow,Bvecs(1)%ncol)
    enddo

    if((nnew .GT. 1).AND. .NOT.molcfg%solver%cfg_unres)then
      call make_lintran_vecsArray(molcfg,D,S,F,Bvecs,sigmas,rhos,make_rhos,nnew)
    else
     do i = 1, nnew
      call make_lintran_vecs(molcfg,D,S,F,Bvecs(i),sigmas(i),rhos(i),make_rhos)
      !if (molcfg%solver%info_rsp) then
      !  write(molcfg%lupri,*) 'After make_lintran_vecs: SIGMA Matrix ', i
      !  call mat_print(sigmas(i),1,S%nrow,1,S%nrow,molcfg%lupri)
      !  write(molcfg%lupri,*) 'After make_lintran_vecs: RHO Matrix ', i
      !  call mat_print(rhos(i),1,S%nrow,1,S%nrow,molcfg%lupri)
      !endif
      !call mat_write_to_disk(lusigma_rsp,sigmas(i))
      !call mat_write_to_disk(lurho_rsp,rhos(i))
     enddo
    endif

  end subroutine transform_vectors

!> \brief Interface routine to carry out linear transformations
!> \author S. Coriani, T. Kjaergaard 
!> \date June 2003
!>  
!>   SIGMAS(I) = E[2]*N(I) (sigma)
!>   RHOS(I) = S[2]*N(I) (rho) > check omega of RPA
!>
  subroutine transform_vectors_symAsym(molcfg,D,S,F,&
       & nnew_sym,bvecs_sym,sigmas_sym,rhos_sym,&
       & nnew_asym,bvecs_asym,sigmas_asym,rhos_asym,&
       & make_rhos)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Density matrix
    type(Matrix), intent(in) :: D
    !> Overlap matrix
    type(Matrix), intent(in) :: S
    !> Fock/KS matrix
    type(Matrix), intent(in) :: F
    !> Number of trial vectors in bvecs to be linearly transformed
    integer, intent(in) :: nnew_sym
    !> Trial vectors in bvecs to be linearly transformed
    type(Matrix), intent(in) :: bvecs_sym(:)
    !> Sigma part of linear transformations (output)
    type(Matrix), pointer :: sigmas_sym(:)
    !> Rho part of linear transformations, if make_rhos = true (output)
    type(Matrix), pointer :: rhos_sym(:)
    !> Number of trial vectors in bvecs to be linearly transformed
    integer, intent(in) :: nnew_asym
    !> Trial vectors in bvecs to be linearly transformed
    type(Matrix), intent(in) :: bvecs_asym(:)
    !> Sigma part of linear transformations (output)
    type(Matrix), pointer :: sigmas_asym(:)
    !> Rho part of linear transformations, if make_rhos = true (output)
    type(Matrix), pointer :: rhos_asym(:)
    !> True if the rho part of linear transformation should be calculated
    logical, intent(in) :: make_rhos
    type(Matrix) :: scr
    integer :: i,ndim,l,nnew 
    ndim = S%nrow
!    if(make_rhos)then
       call mem_alloc(rhos_sym,nnew_sym)
!    endif
    call mem_alloc(sigmas_sym,nnew_sym)
    do i=1,nnew_sym
!    if(make_rhos)then
       call mat_init(rhos_sym(i),Bvecs_sym(1)%nrow,Bvecs_sym(1)%ncol)
!    endif
       call mat_init(sigmas_sym(i),Bvecs_sym(1)%nrow,Bvecs_sym(1)%ncol)
    enddo
!    if(make_rhos)then
       call mem_alloc(rhos_asym,nnew_asym)
!    endif
    call mem_alloc(sigmas_asym,nnew_asym)
    do i=1,nnew_asym
!    if(make_rhos)then
       call mat_init(rhos_asym(i),Bvecs_asym(1)%nrow,Bvecs_asym(1)%ncol)
!    endif
       call mat_init(sigmas_asym(i),Bvecs_asym(1)%nrow,Bvecs_asym(1)%ncol)
    enddo

    if(molcfg%solver%cfg_unres)then
       call lsquit('transform_vectors_symAsym not implemented for unres',-1)
    endif

    call make_lintran_vecsArray_symAsym(molcfg,D,S,F,Bvecs_sym,sigmas_sym,rhos_sym,&
         & Bvecs_asym,sigmas_asym,rhos_asym,make_rhos,nnew_sym,nnew_asym)

!    alternative
!    do i = 1, nnew_sym
!       call make_lintran_vecs(molcfg,D,S,F,Bvecs_sym(i),sigmas_sym(i),rhos_sym(i),make_rhos)
!    enddo
!    do i = 1, nnew_sym
!       call make_lintran_vecs(molcfg,D,S,F,Bvecs_asym(i),sigmas_asym(i),rhos_asym(i),make_rhos)
!    enddo

  end subroutine transform_vectors_symAsym
  
!> \brief Linear transformation routine
!> \author S. Coriani
!> \date June 2003
!>  
!> Calculate the rho and sigma parts of linear transformation:
!> rho = S[2]b = - S[b,D]_s S = - S(bSD-DSb)S = -SbSDS+SDSbS
!> sigma = E[2]b = F[b,D]_s S-S[b,D]_s F+G([b,D])DS-SDG([b,D])
!> for one b at a time.
!>
  subroutine make_lintran_vecs(molcfg,D,S,F,Bvecs_i,sigma_i,rho_i,make_rhos)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Density matrix
    type(Matrix), intent(in) :: D
    !> Overlap matrix
    type(Matrix), intent(in) :: S
    !> Fock/KS matrix
    type(Matrix), intent(in) :: F
    !> Trial vector to be linearly transformed
    type(Matrix), intent(in) :: Bvecs_i
    !> Sigma part of linear transformation (output)
    type(Matrix), intent(inout) :: sigma_i
    !> Rho part of linear transformation, if make_rhos = true (output)
    type(Matrix), intent(inout) :: rho_i
    !> If true, construct rho part of linear transformation
    logical, intent(in) :: make_rhos
    type(matrix) :: prod(1),prod2(1),GbDs(1)
    integer :: ndim
    logical :: cov
!Test: declarations used in testing
    real(realk), allocatable :: CMOfull(:), DiffEn(:), WrkScr(:)
    real(realk), allocatable :: Sigmafull(:)
    real(realk) :: T1, T2
  
    molcfg%solver%rsp_nlintra = molcfg%solver%rsp_nlintra + 1 !Count no of linear transformations
    ndim = S%nrow
    call mat_init(prod(1),ndim,ndim)
    call mat_init(prod2(1),ndim,ndim)
    call mat_init(GbDs(1),ndim,ndim)

    call ABCcommutator(ndim,bvecs_i,D,S,prod2(1))
    if (make_rhos) then
       call mat_mul(prod2(1),S,'n','n',1E0_realk,0E0_realk,prod(1))
       call mat_mul(S,prod(1),'n','n',2E0_realk,0E0_realk,rho_i) !Sign changed 15/7-09!
    endif
    call ABCcommutator(ndim,F,S,prod2(1),sigma_i)
    
    call di_GET_GbDs_and_XC_linrsp(GbDs,prod,molcfg%lupri,&
         & molcfg%luerr,prod2,1,ndim,D,&
         & molcfg%setting%do_dft,molcfg%setting)

!    call di_GET_GbDs(molcfg%lupri,molcfg%luerr,& 
!         &prod2(1),GbDs(1),molcfg%setting)
!   	 if (molcfg%setting%do_dft) THEN 
!       !Add extra G contributions 
!       call II_get_xc_linrsp(molcfg%lupri,molcfg%luerr,& 
!            &molcfg%setting,ndim,prod2,D,prod,1) 
!       call mat_daxpy(1E0_realk,prod(1),GbDs(1)) 
!    endif
    
    call ABCcommutator(ndim,GbDs(1),S,D,prod(1))
    call mat_DAXPY(1E0_realk,prod(1),sigma_i)
    call mat_scal(2.0E0_realk,sigma_i)
    
    ! Project out redundancies
    call util_scriptPx('T',D,S,sigma_i)
    
    if (make_rhos) then
       call util_scriptPx('T',D,S,rho_i)
    end if

    !FREE matrices!
    call mat_free(prod(1))
    call mat_free(prod2(1))
    call mat_free(GbDs(1))
  end subroutine make_lintran_vecs

!> \brief See make_lintran_vecs - this one is just for many at a time
!> \author T. Kjaergaard
!> \date 2009
  subroutine make_lintran_vecsArray_symAsym(molcfg,D,S,F,&
       & Bvecs_sym,sigma_sym,rho_sym,&
       & Bvecs_asym,sigma_asym,rho_asym,&
       & make_rhos,nnew_sym,nnew_asym)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Density matrix
    type(Matrix), intent(in) :: D
    !> Overlap matrix
    type(Matrix), intent(in) :: S
    !> Fock/KS matrix
    type(Matrix), intent(in) :: F
    !> Number of trial vectors in Bvecs
    integer, intent(in)      :: nnew_sym
    !> Trial vector to be linearly transformed
    type(Matrix), intent(in) :: Bvecs_sym(nnew_sym)
    !> Sigma parts of linear transformations (output)
    type(Matrix), intent(inout) :: sigma_sym(nnew_sym)
    !> Rho parts of linear transformations, if make_rhos = true (output)
    type(Matrix), intent(inout) :: rho_sym(nnew_sym)  !output
    !> Number of trial vectors in Bvecs
    integer, intent(in)      :: nnew_asym
    !> Trial vector to be linearly transformed
    type(Matrix), intent(in) :: Bvecs_asym(nnew_asym)
    !> Sigma parts of linear transformations (output)
    type(Matrix), intent(inout) :: sigma_asym(nnew_asym)
    !> Rho parts of linear transformations, if make_rhos = true (output)
    type(Matrix), intent(inout) :: rho_asym(nnew_asym)  !output
    !> If true, construct rho parts of linear transformations
    logical, intent(in) :: make_rhos
    type(matrix) :: prod1
    type(matrix),pointer :: prod2(:),GbDs(:),Sigma(:)
    logical :: cov
    integer :: ndim,i,nnew

    nnew = nnew_sym+nnew_asym
    ndim = S%nrow

    nullify(prod2)
    allocate(prod2(nnew))
    nullify(GbDs)
    allocate(GbDs(nnew))
    nullify(Sigma)
    allocate(Sigma(nnew))
    do i=1,nnew
       call mat_init(prod2(i),ndim,ndim)
       call mat_init(GbDs(i),ndim,ndim)
       call mat_init(Sigma(i),ndim,ndim)
    enddo

    call mat_init(prod1,ndim,ndim)
    molcfg%solver%rsp_nlintra = molcfg%solver%rsp_nlintra + 1 !Count no of linear transformations

    do i=1,nnew_sym
       call ABCcommutator(ndim,bvecs_sym(i),D,S,prod2(i))
       if (make_rhos) then
          call mat_mul(prod2(i),S,'n','n',1E0_realk,0E0_realk,prod1)
          call mat_mul(S,prod1,'n','n',2E0_realk,0E0_realk,rho_sym(i))    !Sign changed 15/7-09!
       endif
    enddo
    do i=1,nnew_asym
       call ABCcommutator(ndim,bvecs_asym(i),D,S,prod2(nnew_sym+i))
       if (make_rhos) then
          call mat_mul(prod2(nnew_sym+i),S,'n','n',1E0_realk,0E0_realk,prod1)
          call mat_mul(S,prod1,'n','n',2E0_realk,0E0_realk,rho_asym(i))    !Sign changed 15/7-09!
       endif
    enddo
    
    call di_GET_GbDs_and_XC_linrsp(GbDs,sigma,molcfg%lupri,&
         & molcfg%luerr,prod2,nnew,ndim,D,&
         & molcfg%setting%do_dft,molcfg%setting)
    do i=1,nnew
       call mat_free(sigma(i))
    enddo
    deallocate(Sigma)
    nullify(Sigma)
    
    do i=1,nnew_sym
       call ABCcommutator(ndim,F,S,prod2(i),sigma_sym(i))    
       call ABCcommutator(ndim,GbDs(i),S,D,prod1)
       call mat_DAXPY(1E0_realk,prod1,sigma_sym(i))
       call mat_scal(2.0E0_realk,sigma_sym(i))
       ! Project out redundancies
       call util_scriptPx('T',D,S,sigma_sym(i))
    
       if (make_rhos) then
          call util_scriptPx('T',D,S,rho_sym(i))
       end if
    enddo
    do i=1,nnew_asym
       call ABCcommutator(ndim,F,S,prod2(nnew_sym+i),sigma_asym(i))    
       call ABCcommutator(ndim,GbDs(nnew_sym+i),S,D,prod1)
       call mat_DAXPY(1E0_realk,prod1,sigma_asym(i))
       call mat_scal(2.0E0_realk,sigma_asym(i))
       ! Project out redundancies
       call util_scriptPx('T',D,S,sigma_asym(i))
    
       if (make_rhos) then
          call util_scriptPx('T',D,S,rho_asym(i))
       end if
    enddo

    !FREE matrices!
    do i=1,nnew
       call mat_free(prod2(i))
       call mat_free(GbDs(i))
    enddo
    call mat_free(prod1)
    deallocate(prod2)
    deallocate(GbDs)
    nullify(GbDs)
    nullify(prod2)

  end subroutine make_lintran_vecsArray_symAsym

!> \brief See make_lintran_vecs - this one is just for many at a time
!> \author T. Kjaergaard
!> \date 2009
  subroutine make_lintran_vecsArray(molcfg,D,S,F,Bvecs,sigma,rho,make_rhos,nnew)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Density matrix
    type(Matrix), intent(in) :: D
    !> Overlap matrix
    type(Matrix), intent(in) :: S
    !> Fock/KS matrix
    type(Matrix), intent(in) :: F
    !> Number of trial vectors in Bvecs
    integer, intent(in)      :: nnew
    !> Trial vector to be linearly transformed
    type(Matrix), intent(in) :: Bvecs(nnew)
    !> Sigma parts of linear transformations (output)
    type(Matrix), intent(inout) :: sigma(nnew)
    !> Rho parts of linear transformations, if make_rhos = true (output)
    type(Matrix), intent(inout) :: rho(nnew)  !output
    !> If true, construct rho parts of linear transformations
    logical, intent(in) :: make_rhos
    type(matrix) :: prod1
    type(matrix),pointer :: prod2(:),GbDs(:)
    logical :: cov
    integer :: ndim,i
!Test: declarations used in testing
    real(realk), allocatable :: CMOfull(:), DiffEn(:), WrkScr(:)
    real(realk), allocatable :: Sigmafull(:)
    real(realk) :: T1, T2

    ndim = S%nrow
    nullify(prod2)
    allocate(prod2(nnew))
    nullify(GbDs)
    allocate(GbDs(nnew))
    do i=1,nnew
       call mat_init(prod2(i),ndim,ndim)
       call mat_init(GbDs(i),ndim,ndim)
    enddo
    call mat_init(prod1,ndim,ndim)
    molcfg%solver%rsp_nlintra = molcfg%solver%rsp_nlintra + 1 !Count no of linear transformations

    do i=1,nnew
       call ABCcommutator(ndim,bvecs(i),D,S,prod2(i))
       if (make_rhos) then
          call mat_mul(prod2(i),S,'n','n',1E0_realk,0E0_realk,prod1)
          call mat_mul(S,prod1,'n','n',2E0_realk,0E0_realk,rho(i)) !Sign changed 15/7-09!
       endif
    enddo
    
	call di_GET_GbDs_and_XC_linrsp(GbDs,sigma,molcfg%lupri,&
						& molcfg%luerr,prod2,nnew,ndim,D,&
						& molcfg%setting%do_dft,molcfg%setting)
!    call di_GET_GbDs(molcfg%lupri,molcfg%luerr,& 
!         &prod2,GbDs,nnew,molcfg%setting) 
!    if(molcfg%setting%do_dft)then
!       call II_get_xc_linrsp(molcfg%lupri,molcfg%luerr,& 
!            &molcfg%setting,ndim,prod2,D,sigma,nnew)   !sigma used as temp mat 
!       do i=1,nnew
!          call mat_daxpy(1E0_realk,sigma(i),GbDs(i))
!       enddo
!    endif

    do i=1,nnew
       call ABCcommutator(ndim,F,S,prod2(i),sigma(i))    
       call ABCcommutator(ndim,GbDs(i),S,D,prod1)
       call mat_DAXPY(1E0_realk,prod1,sigma(i))
       call mat_scal(2.0E0_realk,sigma(i))
       ! Project out redundancies
       call util_scriptPx('T',D,S,sigma(i))
    
       if (make_rhos) then
          call util_scriptPx('T',D,S,rho(i))
       end if
    enddo

    !FREE matrices!
    do i=1,nnew
       call mat_free(prod2(i))
       call mat_free(GbDs(i))
    enddo
    call mat_free(prod1)
    deallocate(prod2)
    deallocate(GbDs)
    nullify(GbDs)
    nullify(prod2)

  end subroutine make_lintran_vecsArray

!> \brief Solve rpa problem in reduced subspace.
!> \author S. Coriani
!> \date June 2003
!> 
!> GD(gradient),BVECS(new subspace basis vectors),SIGMAS (new lin.transf.trial vecs),
!> RHOS (new rho vectors) are sparse (the old ones are on disk!)
!> RED_GD,RED_E,RED_S,EIVEC are not
!> N_GD_OR_EXCI is the number of response vectors or of excitations/eigenvectors to be found
!> EIVAL is an array of size N_GD_OR_EXCI:
!>       lineq=true: eival contains input frequencies 
!>       lineq=false: array for excitation energies
!>
!> PURPOSE:
!> Solve rpa problem in reduced subspace.
!> the structure of E[2] and S[2] is used to obtain
!> a reduced problem with pairing of eigensolutions. 
!> Although only ndim_red vectors are known, a reduced
!> problem with 2*(ndim_red + nb_new) vectors is used, since linear transformations 
!> for E[2] BVEC implicitly gives the linear transformation for (BVEC)^T.
!>
  subroutine build_reduced_space_and_solve(molcfg,ndim_red,nb_new,n_gd_or_exci,gd,sigmas,rhos,bvecs,&
                                          & eival,red_eivec)
    implicit none 
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Input: Number of vectors on disk and current half size of reduced space. \n
    !> Output: Updated number of vectors on disk and new half size of reduced space (updated in extend_red_matrices2)
    integer, intent(inout) :: ndim_red
    !> Number of new trial vectors and corresponding linear transformations
    integer, intent(in) :: nb_new
    !> Number of response vectors or of excitations/eigenvectors to be found
    integer, intent(in) :: n_gd_or_exci
    !> Right hand sides = property gradients
    type(Matrix), intent(in) :: gd(:)
    !> E[2] transformed trial vectors (old ones on disk, lusigma_rsp)
    type(Matrix), pointer :: sigmas(:)
    !> S[2] transformed trial vectors (old ones on disk, lurho_rsp)
    type(Matrix), pointer :: rhos(:)
    !> New trial vectors (old ones on disk, lub_rsp)
    type(Matrix), pointer :: bvecs(:)
    !> Reduced space excitation energies, or input frequencies. Dimension(n_gd_or_exci). (Input if LINEQ=true, otherwise output)
    real(realk), intent(inout) :: eival(:)
    !> Eigenvectors/solutions of the new reduced rpa equations. Dimension(ndim_RED,ngrad)
    real(realk), intent(inout) :: red_eivec(:,:)
    integer :: ndim_red_mat,i

    ndim_red_mat = 2*ndim_red
    if (ndim_red_mat > 2*molcfg%solver%rsp_maxred) then
      WRITE(molcfg%LUPRI,'(//A/A,I5,/A,I5)') ' >>> ERROR IN build_reduced.. >>>', &
&     ' DIMENSION OF REDUCED SPACE IS  ',ndim_red_mat, &
&     ' WHICH EXCEEDS ALLOWED DIMENSION',2*molcfg%solver%rsp_MAXRED
      STOP 'build_reduced: TOO LARGE DIMENSION OF REDUCED SPACE'
    END IF
!***************************************************************
! Section 1: extend reduced E[2]-AND S[2]-matrices and gradients 
! with N new b-vectors
!***************************************************************
    if (nb_new > 0) then  
       call extend_red_matrices(molcfg,ndim_red,nb_new,n_gd_or_exci,gd,sigmas,rhos,bvecs)
       do i=1,size(bvecs)
          call mat_free(bvecs(i))
          call mat_free(sigmas(i))
          call mat_free(rhos(i))
       enddo
       call mem_dealloc(bvecs)
       call mem_dealloc(sigmas)
       call mem_dealloc(rhos)
    endif
! ************************************************************
! Section 2: find eigenvalues and -vectors of reduced L-matrix
! ************************************************************
    if (LINEQ) then 
       if (.NOT. molcfg%solver%rsp_quiet) write(molcfg%lupri,*)'omega in build_reduc.. ',eival(1) 
       if (nmcdvec /= 0) then
          call solve_red_lineqMCD(molcfg,ndim_red,n_gd_or_exci,eival,red_eivec)
       else
          call solve_red_lineq(molcfg,ndim_red,n_gd_or_exci,eival,red_eivec)
       endif
    else
       call solve_red_eigen(molcfg,ndim_red,n_gd_or_exci,eival,red_eivec)
    endif

  end subroutine build_reduced_space_and_solve

!> \brief Extend reduced space with new vectors.
!> \author S. Host, S. Coriani
!> \date October 2006
!> 
!> The reduced E2 and S2 matrices (red_E and red_S) are global variables
!> (see top of file)
!>
  subroutine extend_red_matrices(molcfg,ndim_red,nb_new,ngd,gd,sigmas,rhos,bvecs)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Input: Number of vectors on disk and current half size of reduced space \n
    !> Output: Updated number of vectors on disk and new half size of reduced space
    integer, intent(inout) :: ndim_red
    !> Number of new trial vectors and corresponding linear transformations
    integer, intent(in) :: nb_new
    !> The number of response vectors or of excitations/eigenvectors to be found
    integer, intent(in) :: ngd
    !> if(lineq): matrix array containing ngd right hand sides (gradients). Not referenced if lineq=false
    type(Matrix), intent(in) :: gd(:)
    !> Sigma parts of linear transformations. Dimension(Nb_new)
    type(Matrix), intent(in) :: sigmas(:)
    !> Rho parts of linear transformations. Dimension(Nb_new)
    type(Matrix), intent(in) :: rhos(:)
    !> Trial vectors. Dimension(Nb_new)
    type(Matrix), intent(in) :: bvecs(:)
    type(Matrix) :: b_j, bT, b_jT, sigma_j, sigmaT, rho_j, rhoT,blt
    integer :: i,j,k,l,ndim
    ndim = bvecs(1)%nrow
    call mat_init(bT,ndim,ndim)
    call mat_init(rhoT,ndim,ndim)
    call mat_init(sigmaT,ndim,ndim)
    call mat_init(b_j,ndim,ndim)
    call mat_init(rho_j,ndim,ndim)
    call mat_init(sigma_j,ndim,ndim)
    call mat_init(b_jT,ndim,ndim)
    call mat_init(blT,ndim,ndim)

    if (LINEQ) then
      !EXTEND REDUCED GRADIENTS 
      call extend_gradients(molcfg,ndim_red,nb_new,ngd,gd,bvecs)
    endif

    !Calculate all new elements from new vectors 
    k = 0
    do i = 2*ndim_red+1, 2*ndim_red+2*nb_new, 2  !Run over 2 elements at a time because of pairing
       !k = i - 2*ndim_red
       k = k + 1
       l = 0
       do j = 2*ndim_red+1, 2*ndim_red+2*nb_new, 2
          !l = j - 2*ndim_red
          l = l + 1
          !write(molcfg%lupri,*) '1: k index', k 
          call mat_trans(bvecs(k), bT) 
          call mat_trans(sigmas(l), sigmaT)
          call mat_trans(rhos(l), rhoT)
          call mat_trans(bvecs(l), blT)
          red_E(i,j)     = mat_dotproduct(bvecs(k),sigmas(l))
          red_E(i+1,j)   = mat_dotproduct(bT,sigmas(l))
          red_E(i,j+1)   = mat_dotproduct(bvecs(k),sigmaT)
          red_E(i+1,j+1) = mat_dotproduct(bT,sigmaT)

          red_S(i,j)     = mat_dotproduct(bvecs(k),rhos(l))
          red_S(i+1,j)   = mat_dotproduct(bT,rhos(l))
          red_S(i,j+1)   = -mat_dotproduct(bvecs(k),rhoT)
          red_S(i+1,j+1) = -mat_dotproduct(bT,rhoT)

          b_overlaps(i,j)     = mat_dotproduct(bvecs(k),bvecs(k))
          b_overlaps(i+1,j)   = mat_dotproduct(bT,bvecs(l))
          b_overlaps(i,j+1)   = mat_dotproduct(bvecs(k),blT)
          b_overlaps(i+1,j+1) = mat_dotproduct(blT,blT)

          call mat_trans(rhos(k), rhoT)
          bS_overlaps(i,j)     = mat_dotproduct(rhos(k),bvecs(k))
          bS_overlaps(i+1,j)   = -mat_dotproduct(rhoT,bvecs(k))
          bS_overlaps(i,j+1)   = mat_dotproduct(rhos(k),bT)
          bS_overlaps(i+1,j+1) = -mat_dotproduct(rhoT,bT)
       enddo      
    enddo

    !Setup lower half of E2 and S2:
    rewind(lusigma_rsp) ; rewind(lurho_rsp)
    do j = 1, 2*ndim_red, 2
       call mat_read_from_disk(lusigma_rsp,sigma_j,RSPOnMaster)
       call mat_read_from_disk(lurho_rsp,rho_j,RSPOnMaster)
       call mat_trans(sigma_j, sigmaT) 
       call mat_trans(rho_j, rhoT)
       k = 0
       do i = 2*ndim_red+1, 2*ndim_red+2*nb_new, 2  !Run over 2 elements at a time because of pairing
          k = k + 1
          !writemolcfg%(lupri,*) '2: k index', k
          call mat_trans(bvecs(k), bT)
          red_E(i,j)     = mat_dotproduct(bvecs(k),sigma_j) 
          red_E(i+1,j)   = mat_dotproduct(bT,sigma_j)
          red_E(i,j+1)   = mat_dotproduct(bvecs(k),sigmaT)
          red_E(i+1,j+1) = mat_dotproduct(bT,sigmaT)

          red_S(i,j)     = mat_dotproduct(bvecs(k),rho_j)
          red_S(i+1,j)   = mat_dotproduct(bT,rho_j)
          red_S(i,j+1)   = -mat_dotproduct(bvecs(k),rhoT)
          red_S(i+1,j+1) = -mat_dotproduct(bT,rhoT)
       enddo
    enddo

    rewind(lub_rsp)
    !Explicitly calculate upper half of E2 and S2:
    do i = 1, 2*ndim_red, 2
       call mat_read_from_disk(lub_rsp,b_j,RSPOnMaster)
       call mat_trans(b_j, b_jT)
       k = 0
       do j = 2*ndim_red+1, 2*ndim_red+2*nb_new, 2  !Run over 2 elements at a time because of pairing
          k = k + 1
          call mat_trans(bvecs(k), bT)
          call mat_trans(sigmas(k), sigmaT) 
          call mat_trans(rhos(k), rhoT)

          red_E(i,j)     = mat_dotproduct(b_j,sigmas(k)) 
          red_E(i+1,j)   = mat_dotproduct(b_jT,sigmas(k))
          red_E(i,j+1)   = mat_dotproduct(b_j,sigmaT)
          red_E(i+1,j+1) = mat_dotproduct(b_jT,sigmaT)

          red_S(i,j)     = mat_dotproduct(b_j,rhos(k)) 
          red_S(i+1,j)   = mat_dotproduct(b_jT,rhos(k))
          red_S(i,j+1)   = -mat_dotproduct(b_j,rhoT)
          red_S(i+1,j+1) = -mat_dotproduct(b_jT,rhoT)

          b_overlaps(i,j)     = mat_dotproduct(bvecs(k),b_j) 
          b_overlaps(i+1,j)   = mat_dotproduct(bT,b_j)
          b_overlaps(i,j+1)   = mat_dotproduct(bvecs(k),b_jT)
          b_overlaps(i+1,j+1) = mat_dotproduct(bT,b_jT)
          b_overlaps(j,i)     = b_overlaps(i,j)     
          b_overlaps(j+1,i)   = b_overlaps(i+1,j)   
          b_overlaps(j,i+1)   = b_overlaps(i,j+1)   
          b_overlaps(j+1,i+1) = b_overlaps(i+1,j+1) 

          bS_overlaps(i,j)     = mat_dotproduct(rhos(k),b_j) 
          bS_overlaps(i+1,j)   = mat_dotproduct(rhoT,b_j)
          bS_overlaps(i,j+1)   = mat_dotproduct(rhos(k),b_jT)
          bS_overlaps(i+1,j+1) = mat_dotproduct(rhoT,b_jT)
          bS_overlaps(j,i)     = bS_overlaps(i,j)     
          bS_overlaps(j+1,i)   = bS_overlaps(i+1,j)   
          bS_overlaps(j,i+1)   = bS_overlaps(i,j+1)   
          bS_overlaps(j+1,i+1) = bS_overlaps(i+1,j+1) 
       enddo
    enddo

    ndim_red = ndim_red + nb_new !update reduced space dimension
    !Remember that the row-dimension of the reduced matrices 
    !is 2*ndim_red, where ndim_red is the number of trials

    if (MOLCFG%SOLVER%INFO_RSP_REDSPACE) then
       write (molcfg%lupri,*) 'b vector overlaps, extend_red_matrices2:'
       call LS_OUTPUT(b_overlaps, 1, 2*ndim_red, 1, 2*ndim_red, 2*molcfg%solver%rsp_maxred, &
            &2*molcfg%solver%rsp_maxred, 1, molcfg%lupri)
       write (molcfg%lupri,*) 'b vector overlaps over the METRIC, extend_red_matrices2:'
       call LS_OUTPUT(bS_overlaps, 1, 2*ndim_red, 1, 2*ndim_red, 2*molcfg%solver%rsp_maxred,& 
       &2*molcfg%solver%rsp_maxred, 1, molcfg%lupri)
       write (molcfg%lupri,*) 'Reduced E2, extend_red_matrices2:'
       call LS_OUTPUT(red_E, 1, 2*ndim_red, 1, 2*ndim_red, 2*molcfg%solver%rsp_maxred, &
            &2*molcfg%solver%rsp_maxred, 1, molcfg%lupri)
       write (molcfg%lupri,*) 'Reduced S2, extend_red_matrices2:'
       call LS_OUTPUT(red_S, 1, 2*ndim_red, 1, 2*ndim_red, 2*molcfg%solver%rsp_maxred, &
            &2*molcfg%solver%rsp_maxred, 1, molcfg%lupri)
    endif

    !Add the new sigma, rho, and b vectors to file:
    do i = 1,Nb_new       
       call mat_write_to_disk(lub_rsp,bvecs(i),RSPOnMaster)
       call mat_write_to_disk(lusigma_rsp,sigmas(i),RSPOnMaster)
       call mat_write_to_disk(lurho_rsp,rhos(i),RSPOnMaster)
       if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Extend_red_matrices: set of sigma, rho, and bvector is written to disk'
    enddo

    call mat_free(bT)
    call mat_free(rhoT)
    call mat_free(sigmaT)
    call mat_free(b_j)
    call mat_free(rho_j)
    call mat_free(sigma_j)
    call mat_free(b_jT)
    call mat_free(blT)

  end subroutine extend_red_matrices

!> \brief Extend reduced gradient with new vectors.
!> \author S. Coriani
!> \date June 2003
!> 
!> The reduced gradient (red_GD) is a global variable
!> (see top of file)
!> No output in argument list, but global variable red_GD is updated
!> 
  subroutine extend_gradients(molcfg,ndim_red,nb_new,ngd,gd,bvecs)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Number of vectors on disk and current half size of reduced space
    integer,intent(in) :: ndim_red
    !> Number of new trial vectors
    integer,intent(in) :: nb_new
    !> Number of response vectors to be found
    integer,intent(in) :: ngd
    !> Matrix array containing ngd right hand sides (gradients)
    type(Matrix), intent(in) :: gd(:)
    !> Matrix arrays of the Nb_new new trial vectors
    type(Matrix), intent(in) :: bvecs(:)
    type(Matrix) :: bT

    integer :: ndim, i, j, k

    ndim = gd(1)%nrow 
    call mat_init(bT,ndim,ndim)

    do j = 1, ngd  !Loop over number of gradients
       k=0
       do i = 2*ndim_red+1, 2*ndim_red+2*nb_new, 2 !Loop in steps of 2 because of pairing
          k = k+1
          call mat_trans(bvecs(k),bT)
          red_gd(i,j) = mat_dotproduct(bvecs(k),gd(j))
          red_gd(i+1,j) = mat_dotproduct(bT,gd(j))
       enddo
    enddo

    call mat_free(bT)

    if (molcfg%solver%info_rsp) then
      write(molcfg%lupri,*) 'extend_gradients: after gradient extension'
      call LS_OUTPUT(red_gd, 1, 2*(ndim_red+nb_new), 1, ngd, 2*molcfg%solver%rsp_maxred, &
           &2*molcfg%solver%rsp_maxgd, 1, molcfg%lupri)
    endif

  end subroutine extend_gradients

!> \brief Solve reduced linear response problem in subspace.
!> \author S. Coriani
!> \date June 2003
!> 
!> We consider one gradient/solution at a time.
!>
  subroutine solve_red_lineq(molcfg,ndim_red,ngd,freq,red_X)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Half size of reduced space
    integer, intent(in) :: ndim_red
    !> Number of right hand sides
    integer, intent(in) :: ngd
    !> The ngd frequencies
    real(realk),intent(in) :: freq(:)
    !> The ngd reduced space solution vectors (output)
    real(realk),intent(inout) :: red_X(:,:)
    real(realk),allocatable :: E2(:,:), S2(:,:), RHS(:) 
    integer,allocatable :: IPIV(:)
    integer :: igd, ierr
    ierr=0

    allocate(E2(2*ndim_red,2*ndim_red), S2(2*ndim_red,2*ndim_red))
    allocate(RHS(2*ndim_red), IPIV(2*ndim_red))


    do igd = 1, ngd
       !Setup reduced E2, S2, and right hand side with proper dimension
       !This must be done inside the igd loop because DGESV destroys it. TK
       E2 = red_E(1:2*ndim_red,1:2*ndim_red)
       S2 = red_S(1:2*ndim_red,1:2*ndim_red)

       RHS = red_GD(1:2*ndim_red,igd)

       if (MOLCFG%SOLVER%INFO_RSP_REDSPACE) then
          write (molcfg%lupri,*) "E2, solve_red_lineq2:"
          call LS_OUTPUT(E2, 1, 2*ndim_red, 1, 2*ndim_red, 2*ndim_red, 2*ndim_red, 1, molcfg%lupri)
          write (molcfg%lupri,*) "S2, solve_red_lineq2:"
          call LS_OUTPUT(S2, 1, 2*ndim_red, 1, 2*ndim_red, 2*ndim_red, 2*ndim_red, 1, molcfg%lupri)
       endif
       E2 = E2 - freq(igd)*S2   

       if (MOLCFG%SOLVER%INFO_RSP_REDSPACE) then
          write (molcfg%lupri,*) 'E2 - omega*S2:'
          call LS_OUTPUT(E2, 1, 2*ndim_red, 1, 2*ndim_red, 2*ndim_red, 2*ndim_red, 1, molcfg%lupri)
  
          write (molcfg%lupri,*) 'RHS:'
          call LS_OUTPUT(RHS, 1, 2*ndim_red, 1, 1, 2*ndim_red, 1, 1, molcfg%lupri)
       endif

       !Solve set of linear equations Ax = b:
       call DGESV(2*ndim_red, 1, E2, 2*ndim_red, IPIV, RHS, 2*ndim_red, IERR) !Solution vector is found in RHS.
       if (IERR /= 0) then
          WRITE(molcfg%LUPRI,'(/A, i4)') &
          &     'Problem in DGESV, IERR = ', IERR
          CALL LSQUIT(' Problem in DGESV',molcfg%lupri)
       endif

       if (MOLCFG%SOLVER%INFO_RSP_REDSPACE) then
          write (molcfg%lupri,*) 'Solution vector, solve_red_lineq:'
          call LS_OUTPUT(RHS, 1, 2*ndim_red, 1, 1, 2*ndim_red, 1, 1, molcfg%lupri)
       endif
       red_X(1:2*ndim_red,igd) = RHS
    enddo

    deallocate(E2, S2)
    deallocate(RHS, IPIV)

  end subroutine solve_red_lineq

!> \brief Construct residual, test convergence, and generate new trial vectors.
!> \author S. Coriani
!> \date June 2003
!> 
!> Purpose 
!> 1) Construct residual R=(E(2)-W(I)*S(2))*X(I) (- GD, if lineq)
!>    for NGD eigenvectors/rsp vectors X(I) of reduced rsp problem
!> 2) Test for convergence of ngd eigenvectors,
!>    Convergence criterion:
!>    ||(E(2)-W(I)*S(2))*X(I)|| < molcfg%solver%rsp_conv_thr * ||X(I)|| !!!IS THIS RIGHT????
!> 3) Generate new trial vectors
!>    from PCG residuals plus orthogonalization
!>
!> Based on old subroutine RSPNEX.
!>
  subroutine get_residuals_and_new_b(molcfg,D,S,itmic,ndim_red,ngd,gd,&
          &red_eivec,eival,Nb_new,bvecs,conv,exit_loop)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg),intent(inout) :: molcfg
    !> Density matrix
    type(Matrix),intent(in) :: D
    !> Overlap matrix
    type(Matrix),intent(in) :: S
    !> Number of the microiteration (for printout)
    integer,intent(in) :: itmic
    !> Half size of reduced space
    integer,intent(in) :: ndim_red
    !> Number of right hand sides
    integer,intent(in) :: ngd
    !> Matrix array containing ngd right hand sides (gradients)
    type(Matrix),intent(in) :: gd(:)
    !> The ngd reduced space solution vectors
    real(realk), intent(in) :: red_eivec(:,:)
    !> The ngd frequencies if lineq=true, otherwise the ngd excitation energies 
    real(realk), intent(inout) :: eival(:)
    !> Number of new trial vectors
    integer, intent(inout)  :: Nb_new
    !> Matrix arrays of the Nb_new new trial vectors
    type(Matrix),pointer :: bvecs(:) !output
    !> Set true if converged, and microiterations should be stopped
    logical, intent(out)    :: conv
    !> Set true if we have no new vectors due to linear dependencies, and microiterations should be stopped  
    logical, intent(out)    :: exit_loop
!local
    type(Matrix) :: residuals(ngd)
    type(Matrix),pointer :: residuals2(:)
    real(realk) :: wibx,res_norm,x_red_norm,conv_test,ddot
    integer :: i,n_not_conv,ndim_red_mat,ndim,index_conv(ngd),ibx,max_it,max_red

    type(Matrix) :: Sigma_scr, Rho_scr, X_scr
    real(realk)  :: av_norm,rsp_thresh
    logical :: conv_root, all_conv,DoSVD
!
    type(Matrix) :: tmp,tmp2,tmp3,sigma
    real(realk),pointer :: SV(:),optwrk(:)
    integer,pointer :: IWORK(:)
    integer :: ibx2,kk,info,lwork,nSVD,MAXnSVD
    info=0
    DoSVD = molcfg%solver%doSVD
    IF(DoSVD)WRITE(molcfg%lupri,*)'Use SVD for the response solver. TK'
!for restart
!    integer :: lurestart, j
 
    max_it= molcfg%solver%rsp_maxit
    max_red= molcfg%solver%rsp_maxred
    rsp_thresh= molcfg%solver%rsp_thresh

    ndim = S%nrow
    ndim_red_mat = 2*ndim_red
    n_not_conv = 0
    
    call mat_init(Sigma_scr,ndim,ndim)
    call mat_init(Rho_scr,ndim,ndim)
    call mat_init(X_scr,ndim,ndim)
    call mat_zero(X_scr)
    call mat_zero(Rho_scr)
    call mat_zero(Sigma_scr)
    !
    ! find E[2]*X(i)-w(i) S[2]*X(i) as linear combination of 
    ! previous Sigmas and RHOs. 
    ! X(i) = sum_j Bvecold_j REDX_j(i)
    ! Sigmanew(i) = E[2]X(i) = sum_j REDX_j(i) E[2]Bvecold_j =
    !             = sum_j REDX_j(i) Sigmaold_j
    ! and add to previous -w_i S[2]*X(I) vector
    ! 
    !    lurestart = -1 
    !    CALL LSOPEN(lurestart,'rsp_eigenvecs','unknown','UNFORMATTED')
    !    rewind(lurestart)
    do ibx = 1,ngd
       !get the frequency
       wibx = eival(ibx)
       !write(molcfg%lupri,*) 'eival(ibx), 1', eival(ibx)
       call mat_init(residuals(ibx),ndim,ndim)
       call mat_zero(residuals(ibx))
       if (ABS(wibx) > 1E-20_realk) then
          !for nozero frequencies fill in with -wi*S[2]X(ibx) contrib
          call expand_on_basis_minus(molcfg,ndim_red,red_eivec(1:ndim_red*2,ibx),lurho_rsp,residuals(ibx))
          !res is zero for ibx=2 after this call
          call mat_scal(-wibx,residuals(ibx))
       endif
       !add  sum_j ^RX_j(ibx) sigma_j
       if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Get_Residuals calls Expand_on_basis'
       call expand_on_basis(molcfg,ndim_red,red_eivec(1:ndim_red*2,ibx),lusigma_rsp,residuals(ibx))
       !Dump solutions to disk to make it possible to restart:
       !write(lurestart) ibx, eival(ibx)
       !call mat_write_to_disk(lurestart,residuals(ibx))
       !End Dump
       !TEST TEST TEST
       !Tjek at E[2]X giver Sigmanew
       if (LINEQ) then
          !If we are solving linear equations subtract gradient
          call mat_daxpy(-1E0_realk,gd(ibx),residuals(ibx))
       endif
    enddo

    IF(doSVD)THEN
       IF(LINEQ.AND.ngd.EQ.1)THEN
          MAXnSVD = 6
          nSVD = 0
          call mem_alloc(residuals2,MAXnSVD*ngd)
       ELSE
          MAXnSVD = 3
          nSVD = 3
          call mem_alloc(residuals2,MAXnSVD*ngd)
       ENDIF
       call mat_init(tmp,ndim,ndim)
       call mat_init(tmp2,ndim,ndim)
       call mat_init(tmp3,ndim,ndim)
       call mat_init(SIGMA,ndim,ndim)
       ibx2=0
       do ibx = 1,ngd
          !       WRITE(molcfg%lupri,*)'residual  '
          !       call mat_print(residuals(ibx),1,ndim,1,ndim,molcfg%lupri)
          call mat_assign(tmp3,residuals(ibx))
          
          call mem_alloc(SV,ndim)
          SV=0.0E0_realk
          call mem_alloc(optwrk,5) 
          call dgesvd('A','A',ndim,ndim,tmp3%elms,ndim, &
               & SV,tmp%elms,ndim,tmp2%elms,ndim,optwrk,-1,INFO)
          lwork = INT(optwrk(1))
          call mat_assign(tmp3,residuals(ibx))
          call mem_dealloc(optwrk) 
          call mem_alloc(optwrk,lwork) 
          call dgesvd('A','A',ndim,ndim,tmp3%elms,ndim, &
               & SV,tmp%elms,ndim,tmp2%elms,ndim,optwrk,lwork,INFO)
          IF( INFO.GT.0 ) THEN
             WRITE(*,*)'The algorithm computing SVD failed to converge.'
             call lsquit('The algorithm computing SVD failed to converge.',-1)
          ENDIF
          WRITE(molcfg%lupri,*)'SVD',SV
          IF(LINEQ.AND.ngd.EQ.1)THEN
             nSVD = 0
             do i=1,ndim
                IF(ABS(SV(i)).GT.1.0E-3)THEN
                   nSVD = nSVD + 1
                ENDIF
             enddo
             nSVD = MIN(MAXnSVD,nSVD)
             IF(nSVD.EQ.0)nSVD = 1
          ENDIF
          do kk=1,nSVD-1
             ibx2 = ibx2 + 1
             call mat_init(residuals2(ibx2),ndim,ndim)
             call mat_zero(SIGMA)
             SIGMA%elms(kk+(kk-1)*ndim) = SV(kk)
             call mat_mul(tmp,SIGMA,'N','N',1.0E0_realk,0.0E0_realk,tmp3)
             call mat_mul(tmp3,tmp2,'N','N',1.0E0_realk,0.0E0_realk,residuals2(ibx2))
          enddo

          ibx2 = ibx2 + 1
          call mat_init(residuals2(ibx2),ndim,ndim)
          call mat_zero(SIGMA)
          do kk=nSVD,ndim
             SIGMA%elms(kk+(kk-1)*ndim) = SV(kk)
          enddo
          call mat_mul(tmp,SIGMA,'N','N',1.0E0_realk,0.0E0_realk,tmp3)
          call mat_mul(tmp3,tmp2,'N','N',1.0E0_realk,0.0E0_realk,residuals2(ibx2))
          !       WRITE(molcfg%lupri,*)'sigma 3  U S V'
          !       call mat_print(residuals2(ibx2),1,ndim,1,ndim,molcfg%lupri)
          
!          call mat_add(1E0_realk,residuals2(ibx2-2),1E0_realk,residuals2(ibx2-1),tmp3)
!          call mat_add(1E0_realk,tmp3,1E0_realk,residuals2(ibx2),tmp2)
          
          !       WRITE(molcfg%lupri,*)'sigma 1+2+3'
          !       call mat_print(tmp2,1,ndim,1,ndim,molcfg%lupri)
!          call VerifyMatrices(tmp2,residuals(ibx),'residue',1.0E-14_realk,molcfg%lupri)
          
          call mem_dealloc(optwrk) 
          call mem_dealloc(SV)        
       enddo
       call mat_free(tmp)
       call mat_free(tmp2)
       call mat_free(tmp3)
       call mat_free(SIGMA)
    endif
!    call LSCLOSE(lurestart,'KEEP')

! the residual(s) is now done: test for convergence. If not converged
! form new trial(s)  in next subroutine
! New trial(s) is equal to preconditioned residual(s)

    !Stinne & Poul, 28/9-06
    !If dynamic convergence threshold is requested: Determine average value of initial residuals
    !and set molcfg%solver%rsp_conv_thr to this value times the factor requested in input (1.0E-1_realk, 1.0E-2_realk, or 1.0E-3_realk).
    if (molcfg%solver%rsp_convdyn) then
       if (ITMIC == 1) then
          molcfg%solver%rsp_dyn_thresh=0E0_realk
          av_norm = 0.0E0_realk
          do i = 1, ngd
             res_norm = sqrt(mat_sqnorm2(residuals(i)))
             av_norm = av_norm + res_norm
          enddo
          av_norm = av_norm/ngd
          molcfg%solver%rsp_dyn_thresh = av_norm*molcfg%solver%rsp_conv_factor
          rsp_thresh= molcfg%solver%rsp_dyn_thresh
          WRITE(molcfg%lupri,*) 'Dynamic response convergence threshold set to', molcfg%solver%rsp_dyn_thresh
       else
          rsp_thresh= molcfg%solver%rsp_dyn_thresh
       endif
    else
       if (ITMIC == 1) then
          WRITE(molcfg%lupri,*) 'Response convergence threshold is', molcfg%solver%rsp_thresh
       endif
       rsp_thresh= molcfg%solver%rsp_thresh
    endif
    
    do ibx = 1,ngd       
       conv_root = .false.
       res_norm = sqrt(mat_sqnorm2(residuals(ibx)))
       if (res_norm < rsp_thresh) conv_root = .true.

       if (lineq) then
          write (molcfg%lupri, '("Residual norm for vector ", i3, " is: ", E14.6, "&
               &    and frequency = " , F12.6, "   It = ", i3, " CONV =", l2)')&
               & ibx, res_norm, EIVAL(ibx), itmic, conv_root
       else
          if(.NOT. molcfg%solver%rsp_quiet)write (molcfg%lupri, '("Residual norm for vector ", i3, " is: ", E14.6, "&
               &    and eigenvalue = ", F12.8, "   It = ", i3,  " CONV =", l2)') &
               &ibx, res_norm, EIVAL(ibx), itmic, conv_root
       endif
       x_red_norm = sqrt(ddot(ndim_red_mat,red_eivec(1,ibx),1,red_eivec(1,ibx),1))
       !
       ! Insert back print statements
       !
       if (lineq) then !Stinne & Poul, 27/9-06: For eigenvalues, it doesn't make sense to divide by x_red_norm
                       !since it is normalized over S2 and therefore is arbitrary
          !conv_test = res_norm / x_red_norm
          conv_test = res_norm !Stinne 1/6-07 otherwise convergence is too sloppy!!!
       else
          conv_test = res_norm
       endif
       !if (molcfg%solver%info_rsp) then
       !   write(molcfg%lupri,*) 'conv_test',conv_test
       !   write(molcfg%lupri,*) 'res_norm',res_norm
       !   write(molcfg%lupri,*) 'x_red_norm',x_red_norm
       !   write(molcfg%lupri,*) 'molcfg%solver%rsp_conv_thr',molcfg%solver%rsp_conv_thr
       !endif
 
       if (conv_test <= rsp_thresh) then
          !this root is converged
          index_conv(ibx) = -1
       else
          !this root is NOT converged
          n_not_conv = n_not_conv + 1
          index_conv(ibx) = 0
          !form new trial from preconditioned residual
          !write(molcfg%lupri,*) 'eival(ibx)', eival(ibx)
          !write(molcfg%lupri,*) 'res in, get_res_and...:'
          !call mat_print(residuals(ibx),1,ndim,1,ndim,molcfg%lupri)
          IF(doSVD)THEN
             do kk=1,nSVD
                call precond_residual(molcfg,S,EIVAL(ibx),residuals2(kk+(ibx-1)*3))
             enddo
          ELSE
             call precond_residual(molcfg,S,EIVAL(ibx),residuals(ibx))
          ENDIF
       endif
    enddo
    !The residuals are now preconditioned!!!
 
    !Remove converged vectors
    nb_new = 0
    all_conv = .true.
    do ibx = 1,ngd
       if (index_conv(ibx) == 0) then !vec not converged
          IF(doSVD)THEN
             do kk=1,nSVD
                nb_new = nb_new + 1
                if (nb_new /= kk+(ibx-1)*3) call mat_assign(residuals2(nb_new),residuals2(kk+(ibx-1)*3))
             enddo
          ELSE
             nb_new = nb_new + 1
             if (nb_new /= ibx) call mat_assign(residuals(nb_new),residuals(ibx))
          endif
          all_conv = .false.
       endif
    enddo
!
!        ORTHOGONALIZE TRIAL VECTORS AND EXAMINE FOR LINEAR DEPENDENCE
!
    !Nb_new is changed 
    if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Nb_new to orthonormalize:',Nb_new
    if (molcfg%solver%info_rsp .and. Nb_new == 0 .and. .not. all_conv) then
       WRITE(molcfg%lupri,*) 'WARNING: no new vectors to orthonormalize!!?'
    endif
    if (Nb_new /= 0) then
       IF(doSVD)THEN
          call orthonormalize(molcfg,D,S,residuals2,Nb_new,Ndim_red,Bvecs)
       ELSE
          call orthonormalize(molcfg,D,S,residuals,Nb_new,Ndim_red,Bvecs)
       ENDIF
    endif
    call verify_conv_prop(molcfg,ngd,n_not_conv,Nb_new,itmic,ndim_red,conv,exit_loop)
!
! Finished, deallocate local arrays
    IF(doSVD)THEN
       do ibx = 1,ibx2
          call mat_free(residuals2(ibx))
       enddo
       call mem_dealloc(residuals2)
    ENDIF
    do ibx = 1,ngd
       call mat_free(residuals(ibx))
    enddo
    call mat_free(Sigma_scr)
    call mat_free(Rho_scr)
    call mat_free(X_scr)

  end subroutine get_residuals_and_new_b

!> \brief Construct new trial vector from preconditioned residual
!> \author S. Coriani
!> \date June 2003
  subroutine precond_residual(molcfg,S,omega,Pres_jr)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Overlap matrix
    type(matrix), intent(in) :: S
    !> Frequency/excitation energy
    real(realk), intent(inout) :: omega
    !> Input: Residual to be preconditioned. Output: New trial vector
    type(Matrix),intent(inout) :: Pres_jr
    type(Matrix) :: scr
    real(realk) :: norm
    integer :: ndim

    ndim = Pres_jr%nrow
    call mat_init(scr,ndim,ndim)
    call mat_zero(scr)
    call mat_assign(scr,Pres_jr) 
    norm = sqrt(mat_sqnorm2(scr))                      
    if (norm > molcfg%solver%rsp_tolerance) then                 
       call rsp_AB_precond(molcfg,scr,S,omega,Pres_jr)
    endif
    call mat_free(scr)
  end subroutine precond_residual

!> \brief Evaluates and prints convergence status.
!> \author L. Thogersen, S. Coriani
!> \date June 2003
!>  
!> Preeetty ugly routine!
!> return conv = .true. if all vectors are converged and exit_loop = .true. if
!> we instead hit some boundry/threshold
!> 
  subroutine verify_conv_prop(molcfg,ngd,n_not_conv,Nb_new,itmic,ndim_red,conv,exit_loop)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg)    :: molcfg
    !> Number of requested solution vectors
    integer, intent(in) :: ngd
    !> Number of not converged vectors
    integer, intent(in) :: n_not_conv
    !> Number of new trial vectors generated in this iteration
    integer, intent(in) :: Nb_new
    !> Iteration number
    integer, intent(in) :: itmic
    !> Half size of reduced space
    integer, intent(in) :: ndim_red
    !> Set true if convergence, and microiterations should be stopped
    logical, intent(out) :: conv
    !> Set true if we have no new vectors due to linear dependencies, and microiterations should be stopped
    logical, intent(out) :: exit_loop
    integer :: kconv, lupri

lupri = molcfg%lupri
!
! ---
! Output:
!
    IF (NGD > molcfg%solver%rsp_MAXVEC ) then
      WRITE (LUPRI,95) NGD,molcfg%solver%rsp_MAXVEC
 95   FORMAT(/,'  ** RSPNEX ** WARNING, NUMBER OF SIMULTANEOUS TRIAL', &
     & ' VECTORS',I5,/,' EXCEEDS THE PRINT MAXIMUM ',I5)
    endif

! transfer outside convergence informations

    if (n_not_conv == 0) then
      !EIGENVECTORS CONVERGED
      !IF (molcfg%solver%info_rsp) then
        WRITE (LUPRI,5030) NGD
 5030 FORMAT(/' *** THE REQUESTED',I15,' SOLUTION VECTORS CONVERGED')
      !endif
      kconv = 1
    else if (NB_NEW == 0) then
      !LINEAR DEPENDENCIES BETWEEN TRIAL VECTORS
      IF (ITMIC >= molcfg%solver%rsp_maxit) THEN
         WRITE(LUPRI,5020)
         KCONV=0
      ELSE
         WRITE(LUPRI,5010)
         KCONV=-1
      END IF
 5010    FORMAT(/' *** MICROITERATIONS STOPPED DUE TO LINEAR', &
     &           ' DEPENDENCE BETWEEN NEW TRIAL VECTORS')
 5020    FORMAT(/' *** MICROITERATIONS STOPPED DUE TO MAX', &
     &           ' ITERATIONS REACHED.')
         WRITE (LUPRI,5031) NGD 
 5031    FORMAT(/' *** REQUESTED',I5,' SOLUTION VECTORS NOT CONVERGED')

    else if (2*ndim_red + 2*nb_new > 2*molcfg%solver%rsp_maxred) then
      !MAXIMUM DIMENSION OF REDUCED SPACE EXCEEDED
      WRITE (LUPRI,'(/A/A,I5,A/)')&
     &   ' *** MICROITERATIONS STOPPED BEFORE CONVERGENCE BECAUSE', &
     &   '     MAXIMUM DIMENSION OF REDUCED SPACE',2*molcfg%solver%rsp_MAXRED,' EXCEEDED.'
         WRITE (LUPRI,5031) NGD
      kconv = -2

    else 
      !MICROITERATIONS NOT CONVERGED
      kconv = 0
    end if

    conv = .false.
    exit_loop = .false.
    IF (KCONV == -2) THEN
       !(MAXIMUM DIMENSION OF REDUCED SPACE EXCEEDED)
       WRITE (LUPRI,'(/A/A)') &
     &          ' *** RSPCTL WARNING-MICROITERATIONS STOPPED BECAUSE', &
     &          '     MAXIMUM DIMENSION OF REDUCED SPACE EXCEEDED.'
       exit_loop = .true.
    ELSE IF (KCONV < 0) THEN
       !(LINEAR DEPENDENCE BETWEEN NEW TRIAL VECTOR )
       WRITE (LUPRI,'(/A/A)') &
&               ' *** RSPCTL WARNING-MICROITERATIONS STOPPED BECAUSE', &
&               '     OF LINEAR DEPENDENCE BETWEEN NEW TRIAL VECTORS'
       exit_loop = .true.
    ELSE IF (KCONV > 0)THEN
      !(CONVERGED)
      IF (molcfg%solver%info_rsp) &
&           WRITE(LUPRI,'(/A)')' *** RSPCTL MICROITERATIONS CONVERGED'
      conv = .true.
    ELSEIF (ITMIC >= molcfg%solver%rsp_maxit) THEN
       !(MAX NO OF MICROITERATIONS REACHED)
       WRITE(LUPRI,'(/A,I4,A)') &
&      ' *** RSPCTL WARNING-MAXIMUM NUMBER OF MICROITERATIONS,', &
&      ITMIC,', REACHED'
       exit_loop = .true.
    ENDIF

end subroutine verify_conv_prop

!> \brief Expand reduced space solution in full basis.
!> \author L. Thogersen, S. Coriani
!> \date June 2003
!> 
!> Expanded_vec must be zeroed outside first time it is used!!!
!>
  subroutine expand_on_basis(molcfg,ndim_red,red_eivec_i,lu_basis,expanded_vec)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Number of vectors on disk
    integer, intent(in) :: ndim_red
    !> Contains the expansion coefficients 
    real(realk), intent(in) :: red_eivec_i(:)
    !> Logical unit number for the file where the basis vectors are found
    integer, intent(in) :: lu_basis
    !> Vector expanded in the set of basis vectors (output)
    type(Matrix), intent(inout) :: expanded_vec  
    integer :: ndim
    type(matrix) :: B_scr, basisvec_k
    real(realk) :: fac1,fac2,bnorm,dprodu
    integer :: k
    logical :: OnMaster

    ndim = expanded_vec%nrow

    call mat_init(B_scr,ndim,ndim)
    call mat_init(basisvec_k,ndim,ndim)

    rewind(lu_basis)
    do k = 1,ndim_red !length of red. solution for each grad
                      !and # of bf's (no pair)
       OnMaster = .TRUE.
      call mat_read_from_disk(lu_basis,basisvec_k,OnMaster) 
      fac1 = red_eivec_i(2*k-1)
      fac2 = red_eivec_i(2*k)
      bnorm = mat_sqnorm2(Basisvec_k)
      if (molcfg%solver%info_rsp) then
        write(molcfg%lupri,*) 'expand_on_basis: bnorm^2 = ', bnorm
        !write(molcfg%lupri,*) 'Vettore di base nro: ', k
        !call mat_print(Basis(k),1,ndim,1,ndim,molcfg%lupri)
      endif
      call mat_daxpy(fac1,Basisvec_k,expanded_vec)
      call mat_trans(Basisvec_k,B_scr)      
      if (molcfg%solver%info_rsp) then
        dprodu = mat_dotproduct(Basisvec_k,B_scr)
        write(molcfg%lupri,*) 'expand_on_basis: bxb^T = ',dprodu
      endif
      call mat_daxpy(fac2,b_scr,expanded_vec)
      !if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'After summation, expanded_vec, k = ',k
      !call mat_print(expanded_vec,1,ndim,1,ndim,molcfg%lupri)
    enddo
    call mat_free(B_scr)
    call mat_free(basisvec_k)

  end subroutine expand_on_basis

!> \brief Expand reduced space solution rho in full basis 
!> \author L. Thogersen, S. Coriani
!> \date June 2003
!> 
!> Expanded_vec must be zeroed outside first time it is used!!!
!>
  subroutine expand_on_basis_minus(molcfg,ndim_red,red_eivec_i,lu_basis,expanded_vec)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Number of vectors on disk
    integer, intent(in) :: ndim_red
    !> Contains the expansion coefficients 
    real(realk), intent(in) :: red_eivec_i(:)
    !> Logical unit number for the file where the basis vectors are found
    integer, intent(in) :: lu_basis
    !> Vector expanded in the set of basis vectors (output)
    type(Matrix), intent(inout) :: expanded_vec  

    type(matrix) :: B_scr, basisvec_k
    real(realk) :: fac1,fac2,bnorm,dprodu
    integer :: k, ndim
    logical :: OnMaster
    OnMaster = .TRUE.
    ndim = expanded_vec%nrow
    call mat_init(B_scr,ndim,ndim)
    call mat_init(basisvec_k,ndim,ndim)

    rewind(lu_basis)
    do k = 1,ndim_red !length of red. solution for each grad
                      !and # of bf's (no pair)
      call mat_read_from_disk(lu_basis,basisvec_k,OnMaster)

      fac1 = red_eivec_i(2*k-1)
      fac2 = red_eivec_i(2*k)

      bnorm = mat_sqnorm2(basisvec_k)
      call mat_daxpy(fac1,basisvec_k,expanded_vec)

      call mat_trans(basisvec_k,B_scr)      
      if (molcfg%solver%info_rsp) then
         dprodu = mat_dotproduct(basisvec_k,B_scr)
         write(molcfg%lupri,*) 'expand_on_basis_minus: bxb^T = ',dprodu
      endif
      call mat_daxpy(-fac2,b_scr,expanded_vec)
    enddo
    call mat_free(B_scr)
    call mat_free(basisvec_k)

  end subroutine expand_on_basis_minus

!> \brief Solve reduced eigenvalue equation in subspace.
!> \author L. Thogersen, S. Coriani
!> \date June 2003
  subroutine solve_red_eigen(molcfg,ndim_red,nexci,eival,red_eivec)
    implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Half size of reduced space
    integer, intent(in) :: ndim_red
    !> Requested number of excitation energies
    integer, intent(in) :: nexci
    !> Array containing the nexci excitation energies
    real(realk),intent(inout) :: eival(:)
    !> The nexci reduced space eigenvectors
    real(realk),intent(inout) :: red_eivec(:,:) !org dim:(2*molcfg%solver%rsp_maxred,2*molcfg%solver%rsp_maxred)
    integer :: ndim_red_mat,i,igd,ij,j,isndx(3),ineg,nsim,kk1
    real(realk) :: ddot, freq
    integer :: ije,jie,kzyrdq, ipos, ierr,matz
    real(realk), allocatable :: alfi(:),beta(:),alfr(:)
    real(realk), allocatable :: E2FULL(:,:),S2FULL(:,:)
    real(realk), allocatable :: red_eivec_scr(:)
    real(realk), allocatable :: scr1(:,:),scr2(:,:)
    real(realk), parameter :: drtest = 1.0E-5_realk
    ndim_red_mat = 2*ndim_red
    allocate(E2FULL(2*ndim_red,2*ndim_red), S2FULL(2*ndim_red,2*ndim_red))
    ALLOCATE(alfi(2*ndim_red),alfr(2*ndim_red),beta(2*ndim_red))
    ALLOCATE(red_eivec_scr(2*ndim_red*2*ndim_red)) !FIXME: check this dimension...
    allocate(scr1(2*ndim_red,2*ndim_red),scr2(2*ndim_red,2*ndim_red))

    E2FULL = red_E(1:2*ndim_red,1:2*ndim_red)
    S2FULL = red_S(1:2*ndim_red,1:2*ndim_red)

!************************************************************
!        Use EISPACK routine for real general matrices in
!        generalized eigenvalue problem
!
!        CALL RGG(NM,N,A,B,ALFR,ALFI,BETA,MATZ,Z,IERR)
!***********************************************************
    matz=1
    call RGG(2*ndim_red,2*ndim_red,E2FULL,S2FULL,ALFR,ALFI,BETA,matz,red_eivec_scr,IERR)

    IF ( IERR /= 0 ) THEN
       WRITE(molcfg%LUPRI,'(/A,I5)') 'solve_red_eigen: Problem in RGG, IERR = ', IERR
       CALL LSQUIT('Problem in RGG',molcfg%lupri)
    END IF
!
! Reorder the reduced eigenvectors
!
    !SONIA: S2 is regenerated because it is destroyed in RGG
    S2FULL = red_S(1:2*ndim_red,1:2*ndim_red)

    CALL ReorderEigenvalues(molcfg,2*ndim_red,S2FULL,RED_EIVEC_scr, & 
    &                       alfr,ALFI,BETA,scr1,scr2,ISNDX)
    if (molcfg%solver%info_rsp_redspace) then
        WRITE (molcfg%LUPRI,'(/A)') ' REDUCED EIGENVECTORS AFTER ORDERING:'
        CALL LS_OUTPUT(RED_EIVEC_scr,1,ndim_red_mat,1,ndim_red_mat, &
    &                     ndim_red_mat,ndim_red_mat,1,molcfg%LUPRI)
    endif

!following printout is also testing stuff
!    if (molcfg%solver%info_rsp_redspace) then
!         WRITE(LUPRI,'(//2A,3I5//A/A)') &
!     &   ' NUMBER OF EIGENVALUES WITH POSITIVE,', &
!     &   ' ZERO, AND NEGATIVE METRIC:',ISNDX, &
!     &   ' THE EIGENVALUES WITH POSITIVE METRIC:', &
!     &   ' NUMBER         EIGENVALUE'
!         IPOS = ISNDX(1)
!         DO I=1,IPOS
!            IF (ABS(alfr(I)) .LT. DRTEST) THEN
!               WRITE(LUPRI,'(I10,1P,D15.8,A/)') I,alfr(I), &
!     &         ' *** RSPRED  WARNING **** ZERO EIGENVALUE.'
!            ELSE
!               WRITE(LUPRI,'(I10,1P,D15.8)') I,alfr(I)
!            END IF
!         END DO
!    endif
    if(2*ndim_red .GT. nexci)then
       eival(1:nexci) = alfr(1:nexci) !Stinne
    else if(2*ndim_red .LE. nexci)then
       eival(1:2*ndim_red) = alfr !Stinne
    endif
! Put residual results in our RED_EIVEC matrix
!
    if (ndim_red_mat > 2*molcfg%solver%rsp_maxred) then !Stinne fix 30/9-06
       WRITE (molcfg%LUPRI,*) 'molcfg%solver%rsp_maxred, ndim_red_mat:', molcfg%solver%rsp_maxred, ndim_red_mat
       WRITE (molcfg%LUPRI,*) 'molcfg%solver%rsp_maxred too small (must be >= 1/2*ndim_red_mat)'
       CALL LSQUIT('molcfg%solver%rsp_maxred too small',molcfg%lupri)
    endif
    if(nexci .LE. ndim_red_mat)then
       do i=1,ndim_red_mat
          do j=1,nexci
             red_eivec(i,j) = red_eivec_scr((j-1)*ndim_red_mat+i)
          enddo
       enddo
    else
       do i=1,ndim_red_mat
          do j=1,ndim_red_mat
             red_eivec(i,j) = red_eivec_scr((j-1)*ndim_red_mat+i)
          enddo
       enddo
    endif
    DEALLOCATE(red_eivec_scr)
    DEALLOCATE(scr1,scr2)
    DEALLOCATE(E2FULL,S2FULL)
    DEALLOCATE(alfi,beta,alfr)
    
  end subroutine solve_red_eigen

!> \brief Analyze and reorder eigenvectors from obtained from subroutine solve_red_eigen.
!> \author L. Thogersen, S. Coriani
!> \date June 2003
!>  
!> Based on old subroutine RSPORD
!> Eigenvalue k is (ALFAR(k)+i*ALFAI(k))/BETA(k)
!> 
  subroutine ReorderEigenvalues(molcfg,NDIM,SRED,EIVEC,ALFAR,ALFAI,BETA,WRK1,WRK2,ISNDX)
  implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg) :: molcfg
    !> Number of eigenvalues/size of reduced space
    integer, intent(in) :: NDIM
    !> Reduced S(2) matrix
    real(realk),intent(in) :: SRED(NDIM,NDIM) 
    !> Eigenvectors in reduced basis
    !real(realk),intent(in) :: EIVEC(NDIM,NDIM)  !SONIA: FIXME
    real(realk),intent(inout) :: EIVEC(NDIM,NDIM)
    !> Real part of numerator of eigenvalues
    real(realk),intent(inout) :: ALFAR(:)
    !> Imaginary part of numerator of eigenvalues
    real(realk),intent(inout) :: ALFAI(:)
    !> Denominator of eigenvalues
    real(realk),intent(inout) :: BETA(:)
    !> Scratch space
    real(realk),intent(inout) :: WRK1(NDIM,NDIM)
    !> Scratch space
    real(realk),intent(inout) :: WRK2(NDIM,NDIM)
    !> Indicates: isndx(1) = Eigenvectors with positive metric
    !>            isndx(2) = Eigenvectors with zero metric
    !>            isndx(3) = Eigenvectors with negative metric
    integer, intent(inout) :: ISNDX(3)

  integer :: ipos, izer, ineg, i, j, jmin, nneg
  real(realk) :: xscale,xsave, amin
  real(realk), PARAMETER :: ZEROT=1.0E-9_realk, COMPLX = 1.0E7_realk
  real(realk), PARAMETER :: D0=0.0E0_realk, D1=1.0E0_realk
!
! IF (IPRRSP .GE. 81) THEN
  if (molcfg%solver%info_rsp) then
         WRITE (MOLCFG%LUPRI,'(//A/A)') &
     &      '  (ALFAR + i ALFAI) / BETA are eigenvalues;', &
     &      '     ALFAR           ALFAI          BETA' 
         WRITE (MOLCFG%LUPRI,'(1P,3D15.6)')(ALFAR(I),ALFAI(I),BETA(I),I=1,NDIM)
  end if
! END IF
  DO I=1,NDIM
     IF(ABS(BETA(I)).GE.ZEROT) THEN
       ALFAR(I)=ALFAR(I)/BETA(I)
       ALFAI(I)=ALFAI(I)/BETA(I)
     ELSE
!      singularities
       WRITE(MOLCFG%LUPRI,1010)I, ALFAR(I),ALFAI(I),BETA(I)
 1010  FORMAT(/' *** SINGULARITY IN REDUCED MCRSP :', &
     &             /'     I,ALFAR(I),ALFAI(I),BETA(I):',I6,1P,3D13.6)
     END IF
!        IF(ABS(ALFAR(I)).GT.ZEROT) THEN
!    &   .OR.(ABS(ALFAR(I)).LE.ZEROT.AND.ABS(ALFAI(I)).GT.ZCRIT)) THEN
!           complex eigenvalues
!        RATIO= ABS(ALFAI(I)/ALFAR(I))
!        IF(RATIO.GT.ZCRIT) WRITE(MOLCFG%LUPRI,1020) ALFAR(I),ALFAI(I)
     IF(ABS(ALFAI(I)).GT.D0) THEN
       WRITE(MOLCFG%LUPRI,1020) I,ALFAR(I),ALFAI(I)
 1020  FORMAT(/' *** COMPLEX EIGENVALUE IN REDUCED MCRSP :', &
     &          /'     REAL AND IMAGINARY PART :   ',I6,1P,2D13.6)
!
! SET EIGENVALUE EQUAL TO COMPLX IN ORDER TO BE ABLE TO SKIP
! CONTRIBUTIONS FROM THIS ROOT WHEN SUMMING UP TERMS
! IN THE CALCULATION OF THE EFFECTIVE SPECTRUM IN C6 CALCULATIONS
!
        ALFAR(I) = COMPLX
     END IF
  END DO
  if (molcfg%solver%info_rsp) then
         WRITE (MOLCFG%LUPRI,'(/A)') ' Eigenvalues of E(2)  :'
         WRITE (MOLCFG%LUPRI,'(I10,1P,D12.2)') (I,ALFAR(I),I=1,NDIM)
  endif
!
!     reduced S(2) in eigenvector basis
!
  call ls_dzero(WRK1,NDIM*NDIM)
  call ls_dzero(WRK2,NDIM*NDIM)
  if (molcfg%solver%info_rsp) then
    WRITE (MOLCFG%LUPRI,*) ' Reduced S(2) before diagonal basis '
    CALL LS_OUTPUT(SRED,1,NDIM,1,NDIM,NDIM,NDIM,1,MOLCFG%LUPRI)
    WRITE (MOLCFG%LUPRI,*) ' Reduced Xvec before diagonal basis '
    CALL LS_OUTPUT(EIVEC,1,NDIM,1,NDIM,NDIM,NDIM,1,MOLCFG%LUPRI)
  endif

  CALL DGEMM('N','N',NDIM,NDIM,NDIM,1E0_realk,SRED,NDIM, &
     &           EIVEC,NDIM,0E0_realk,WRK1,NDIM)
  if (molcfg%solver%info_rsp) then
    WRITE (MOLCFG%LUPRI,*) ' REDS[2]*X'
    CALL LS_OUTPUT(WRK1,1,NDIM,1,NDIM,NDIM,NDIM,1,MOLCFG%LUPRI)
  endif
  CALL DGEMM('T','N',NDIM,NDIM,NDIM,1E0_realk,EIVEC,NDIM, &
     &           WRK1,NDIM,0E0_realk,WRK2,NDIM)
  if (molcfg%solver%info_rsp) then
         WRITE (MOLCFG%LUPRI,'(/A)') ' Reduced S(2) in diagonal basis :'
         WRITE (MOLCFG%LUPRI,'(I10,1P,D12.2)') (I,WRK2(I,I),I=1,NDIM)
!        CALL LS_OUTPUT(WRK2,1,KZYRED,1,KZYRED,KZYRED,KZYRED,1,MOLCFG%LUPRI)
  endif
!
!     select eigenvectors with positive normalization and store them
!     as the first eigenvectors.
!
   IPOS = 0
   IZER = 0
   INEG = 0
   DO I=1,NDIM
      IF (BETA(I) .LE. ZEROT) THEN
         IZER = IZER + 1
      ELSE IF (WRK2(I,I) .GT. D0) THEN
         IPOS = IPOS + 1
         XSCALE= D1/SQRT(WRK2(I,I))
         CALL DSCAL(NDIM,XSCALE,EIVEC(1,I),1)
         IF (IPOS.NE.I) THEN
           CALL DSWAP(NDIM,EIVEC(1,I),1,EIVEC(1,IPOS),1)
           XSAVE = ALFAR(IPOS)
           ALFAR(IPOS)=ALFAR(I)
           ALFAR(I)   =XSAVE
           XSAVE = WRK2(IPOS,IPOS)
           WRK2(IPOS,IPOS) = WRK2(I,I)
           WRK2(I,I)       = XSAVE
           XSAVE = BETA(IPOS)
           BETA(IPOS) = BETA(I)
           BETA(I) = XSAVE
         END IF
      ELSE IF (WRK2(I,I) .LT. -D0) THEN
         INEG = INEG + 1
         XSCALE= D1/SQRT(ABS(WRK2(I,I)))
         CALL DSCAL(NDIM,XSCALE,EIVEC(1,I),1)
         IF (ALFAR(I) .EQ. COMPLX) ALFAR(I) = -COMPLX
      END IF
   END DO
   NNEG = IPOS
   DO I=IPOS+1,NDIM
      IF (ABS(BETA(I)) .GT. ZEROT) THEN
         NNEG = NNEG + 1
         IF (NNEG.NE.I) THEN
            CALL DSWAP(NDIM,EIVEC(1,I),1,EIVEC(1,NNEG),1)
            XSAVE = ALFAR(NNEG)
            ALFAR(NNEG)=ALFAR(I)
            ALFAR(I)   =XSAVE
            XSAVE = WRK2(NNEG,NNEG)
            WRK2(NNEG,NNEG) = WRK2(I,I)
            WRK2(I,I)       = XSAVE
         END IF
      END IF
   END DO
   ISNDX(1) = IPOS
   ISNDX(2) = IZER
   ISNDX(3) = INEG
   IF (IPOS.NE.(NDIM/2)) WRITE (MOLCFG%LUPRI,2020) IPOS,IZER,INEG
 2020 FORMAT(/' *** EIGENVECTORS WITH POSITIVE METRIC:',I6, &
     &       /'     EIGENVECTORS WITH ZERO METRIC:    ',I6, &
     &       /'     EIGENVECTORS WITH NEGATIVE METRIC:',I6)
!
! order eigensolutions in ascending order of eigenvalues.
!
   DO I=1,IPOS
      JMIN = I
      AMIN = ALFAR(I)
      DO J=I+1,IPOS
         IF(ALFAR(J).LT.AMIN) THEN
            AMIN = ALFAR(J)
            JMIN = J
         ENDIF
      END DO
      IF (JMIN.NE.I) THEN
         ALFAR(JMIN)=ALFAR(I)
         ALFAR(I)=AMIN
         CALL DSWAP(NDIM,EIVEC(1,I),1,EIVEC(1,JMIN),1)
      ENDIF
   END DO
   DO I=IPOS+1,IPOS+INEG
      JMIN = I
      AMIN = ALFAR(I)
      DO J=I+1,IPOS+INEG
         IF(ALFAR(J).GT.AMIN) THEN
            AMIN = ALFAR(J)
            JMIN = J
         ENDIF
      END DO
      IF (JMIN.NE.I) THEN
         ALFAR(JMIN)=ALFAR(I)
         ALFAR(I)=AMIN
         CALL DSWAP(NDIM,EIVEC(1,I),1,EIVEC(1,JMIN),1)
      ENDIF
   END DO
   IF (INEG.NE.IPOS) THEN
      WRITE(MOLCFG%LUPRI,'(3(/A))')' ***** WARNING *********' &
     &   ,' number of eigenvalues with negative metric differ from' &
     &   ,' number with positive metric'
      !IF (IPRRSP.GT. 20) THEN
      if (molcfg%solver%info_rsp) then
      WRITE(MOLCFG%LUPRI,'(/A)')'   NUMBER    EIGENVALUE '
      DO I=1,NDIM
         WRITE(MOLCFG%LUPRI,'(I10,5X,1P,D15.8)') I,ALFAR(I)
      END DO
      end if
      !END IF
   ELSE
     if (molcfg%solver%info_rsp) WRITE(MOLCFG%LUPRI,'(/A)') &
     &      '   NUMBER    EIGENVALUE    PAIRED EIGENVALUE'
     DO I=1,IPOS
        IF (ABS(ABS(ALFAR(I))-ABS(ALFAR(IPOS+I))).GT.ZEROT) THEN
            WRITE(MOLCFG%LUPRI,'(/A)')' **WARNING** EIGENVALUES NOT PAIRED'
        END IF
        if (molcfg%solver%info_rsp) &
        !IF (IPRRSP.GT. 20)
     &      WRITE(MOLCFG%LUPRI,'(I10,5X,1P,D15.8,5X,D15.8)') &
     &                      I,ALFAR(I),ALFAR(IPOS+I)
     END DO
     IF (IZER .GT. 0 ) THEN
        WRITE(MOLCFG%LUPRI,'(/A)') &
     &         ' ZERO METRIC EIGENVALUE AND ZERO EIGENVALUE'
        WRITE(MOLCFG%LUPRI,'(/A)') 'NUMBER    EIGENVALUE'
        DO I=1,IZER
           WRITE(MOLCFG%LUPRI,'(I10,5X,1P,D15.8)') I,ALFAR(IPOS+INEG+I)
        END DO
     END IF
   END IF
!
   end subroutine ReorderEigenvalues
!---------------
!***********************************************************
!
!  Subroutines to be used when singularities occur in the 
!  response equations. e.i. MCD and HerzbergTeller 
!  written by Thomas Kjaergaard 
!
!************************************************************
  subroutine MCD_get_1st_vec1(molcfg,F,D,S,ngd,GD,freq,bvecs,Nb_new,lub_mcdvec,Bvec_tmp,nBvec_tmp)
!*******************************************************
! Purpose: 
! generate start vector(s) for solution of a linear
! set of equations from the preconditioned property
! gradient vector(s) 
! The first being the excitation vector
! the second obtaind from the RightHand Side of the eq
! The new,linearly independent orthogonalized vectors
! to be used as trials are returned in BVECS
!*******************************************************
! Modified version of LRST for Linear Scaling. 
! Sonia & Poul, June 2003
! MCD modification Thomas
!*******************************************************
    implicit none
!    type(decompItem),intent(in) :: decomp
    type(rsp_molcfg), intent(inout) :: molcfg
    integer,intent(in)         :: ngd,nBvec_tmp   !number of gradients
    type(Matrix),intent(in)    :: F,D,S
    real(realk), intent(in)    :: freq(rsp_number_of_omegas)
    type(Matrix),intent(inout) :: gd(rsp_number_of_rhs) 
    type(Matrix),pointer :: bvecs(:)
    integer, intent(out)       :: Nb_new
    type(Matrix) :: Bvec_tmp(nBvec_tmp)
    type(Matrix) :: sigma_i,rho_i
    real(realk) :: norm,temps1,temps2,temps3,temps4,q,w,p,g
    integer :: ndim,nb_prev,ibx,i
! thomas
    real(realk) :: dummy_real,S12,S13,S14
    integer     :: dummy_i,lub_rsp_vec,k,lub_mcdvec    

    !print*,'inside DOUBLE TROUBLE'
    !write(molcfg%lupri,*)'inside DOUBLE TROUBLE'
    ndim = S%nrow
    Call Projection(GD,D,S,ngd,freq,Bvec_tmp,nBvec_tmp)

    if (molcfg%solver%info_rsp) then
       write(molcfg%lupri,*)'THE RHS after projection'
       call mat_print(GD(1),1,ndim,1,ndim,molcfg%lupri)
    endif

    call mat_init(sigma_i,ndim,ndim)
    call mat_init(rho_i,ndim,ndim)

    !transform projection vectors and write them to disk
    do i=1,nBvec_tmp
       call mat_write_to_disk(lub_mcdvec,Bvec_tmp(i),RSPOnMaster)
       call make_lintran_vecs(molcfg,D,S,F,Bvec_tmp(i),sigma_i,rho_i,.true.)
       call mat_write_to_disk(lub_mcdvec,rho_i,RSPOnMaster)
    enddo

    call mat_free(sigma_i)
    call mat_free(rho_i)

    k = 0
    Nb_new = nBvec_tmp
    call orthonormalizeMCD(molcfg,D,S,Bvec_tmp,Nb_new,k,bvecs)

  end subroutine MCD_get_1st_vec1

  subroutine MCD_get_1st_vec2(molcfg,F,D,S,ngd,GD,freq,bvecs,Nb_new,cov)
!*******************************************************
! Purpose: 
! generate start vector(s) for solution of a linear
! set of equations from the preconditioned property
! gradient vector(s) 
! The new,linearly independent orthogonalized vectors
! to be used as trials are returned in BVECS
!*******************************************************
! Modified version of LRST for Linear Scaling. 
! Sonia & Poul, June 2003
! MCD modification Thomas
!*******************************************************
    implicit none
    type(rsp_molcfg), intent(inout) :: molcfg
    integer,intent(in)         :: ngd   !number of gradients
    type(Matrix),intent(in)    :: F,D,S
    real(realk), intent(inout) :: freq(rsp_number_of_omegas)
    type(Matrix),pointer :: bvecs(:)
    type(Matrix),intent(inout) :: gd(rsp_number_of_rhs)  
    integer, intent(out)       :: Nb_new
    type(Matrix) :: grad_x,Pgrad, Bvec_tmp(ngd)
    real(realk) :: norm
    integer :: ndim,nb_prev,ibx,i
    logical :: cov

    cov = .FALSE.
    ndim = S%nrow
    call mat_init(grad_x,ndim,ndim)
    call mat_init(Pgrad,ndim,ndim)
    do i = 1,ngd
      call mat_init(Bvec_tmp(i),ndim,ndim)
    enddo
!
! Generate preconditioned gradients as start vectors
!
    if (ngd > rsp_number_of_rhs) then
         WRITE(MOLCFG%LUPRI,'(/A)') &
         &     'ngd > rsp2_number_of_rhs'
         CALL LSQUIT('ngd > rsp2_number_of_rhs', molcfg%lupri)
    endif
    if (ngd > rsp_number_of_omegas) then
         WRITE(MOLCFG%LUPRI,'(/A)') &
         &     'ngd > rsp2_number_of_omegas'
         CALL LSQUIT('ngd > rsp2_number_of_omegas',molcfg%lupri)
    endif
    do ibx = 1, ngd
       call mat_assign(grad_x,gd(ibx))  !STINNE - we do now not divide in symm and antisymm parts
       norm = mat_sqnorm2(grad_x)                  
       if (norm > molcfg%solver%rsp_tolerance) then          
          !Precondition if the matrix exists  by
          !solving A Pgrad = grad
          call rsp_AB_precond(molcfg,grad_x,S,freq(ibx),Pgrad)
          call mat_assign(Bvec_tmp(ibx),Pgrad)
       else
          call mat_assign(Bvec_tmp(ibx),grad_x)
       endif
    enddo
!
! orthogonalize and remove linear dependent vectors
! Nb_new is updated inside orthonormalize
!
    Nb_prev = nmcdvec
    Nb_new = ngd
    if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Nb_new to orthonormalize:',Nb_new

    if (Nb_new == 0) then
       WRITE (MOLCFG%LUPRI,'(//A,I5)') &
                    &  'get_1st_orth_trial_lineq: No start vectors present!'
       CALL LSQUIT('get_1st_orth_trial_lineq: No start vectors present!',molcfg%lupri)
    endif
    !
    call orthonormalize(molcfg,D,S,Bvec_tmp,Nb_new,Nb_prev,bvecs)
  
!
! Number of start trial vectors = Nb_new 
!
    if (Nb_new <= 0) then
       IF(ndim .EQ. Nb_prev+1.OR.ndim .EQ. Nb_prev)then
          cov = .TRUE.
       else
          WRITE (MOLCFG%LUPRI,'(//A,I5)') &
               &  'get_1st_orth_trial_lineq2: START VECTOR IS NOT LINEAR INDEPENDENT.'
          CALL LSQUIT('START VECTOR IS NOT LINEAR INDEPENDENT',molcfg%lupri)
       endif
    endif

!FREE USED MATRICES!!
    call mat_free(grad_x)
    call mat_free(Pgrad)
    do i = 1,ngd
      call mat_free(Bvec_tmp(i))
    enddo
  
  end subroutine MCD_get_1st_vec2

  subroutine orthonormalizeMCD(molcfg,D,S,Bvec_tmp,Nb_new,Nb_prev,bvecs)
!*******************************************************************
! related to old SUBROUTINE RSPORT
! Purpose:
!  Orthogonalize new b-vectors against all previous b-vectors
!  and among themselves, and renormalize.
!  The b-vectors have the form
!        ( Y_dia    Z   )       ! Z_mu_nu  mu < nu ; Y_dia = Y_mu_mu 
!        (  Y     Y_dia )       ! Y_mu_nu  mu > nu    

!  Each b-vector is in (Z, Y) form, the (Y, Z) vector is obtained
!  by transposing.  Each of the new b-vectors is
!  orthogonalized both against (Z, Y) and (Y, Z) for each previous
!  b-vectors.
!  (Orthogonalization is performed twice if round-off is large,
!   see normalize)
!
! In input:
!  previous vectors are on disk in lub_rsp (global variable, see top of file) 
!  BVEC_tmp,  new non-orthogonal vectors
!  NB_NEW,  number of new non-orthogonal b-vectors in BVEC_tmp
!  NB_PREV, number of previous b-vectors (on disk in lub_rsp)
!
! IN output:
!  NB_NEW,  number of acceptable new trial vectors after orthonormalization
!  BVECS, new orthogonalized vectors 
!
! MCD modification Thomas
!**************************************************************************
    implicit none
    type(rsp_molcfg),intent(inout) :: molcfg
    type(Matrix),intent(in) :: D,S
    integer, intent(in)     :: Nb_prev
    integer, intent(inout)  :: Nb_new
    type(Matrix),intent(inout) :: Bvec_tmp(:)
    type(Matrix),pointer :: bvecs(:)
!local 
    integer :: irx,i,iturn,ndim,k,ibvec,jbvec,jrx,no_of_new,lub_rsp_vec,dummy_i
    integer,allocatable :: lin_depend(:) !linear dependency index array
    type(matrix) :: B_scr, b_k
    type(matrix) :: orthovec !Will properly be changed
    real(realk) :: TT,T1,T2,dummy_real
    logical :: run_ortho_again

    if ((Nb_prev + Nb_new) > molcfg%solver%rsp_maxred) then
      WRITE (MOLCFG%LUPRI,*) 'ERROR in orthonormalize: NB_PREV + NB_NEW > molcfg%solver%rsp_maxred'
      WRITE (MOLCFG%LUPRI,*) 'NB_PREV ,NB_NEW  =',NB_PREV,NB_NEW
      WRITE (MOLCFG%LUPRI,*) 'molcfg%solver%rsp_maxred =',molcfg%solver%rsp_maxred
      WRITE (MOLCFG%LUPRI,*) 'Reduce problem size or recompile orthonormalize '// &
      &                  'with larger molcfg%solver%rsp_maxred parameter'
      CALL LSQUIT('orthonormalize error: NB_PREV + NB_NEW > molcfg%solver%rsp_maxred parameter',molcfg%lupri)
    END IF
    
    ndim = S%nrow
    ALLOCATE(lin_depend(Nb_new))
    call mat_init(B_scr,ndim,ndim)
    call mat_init(b_k,ndim,ndim)   

! STEP 1: check initial linear dependencies among new vectors
    if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'orthonormalize -> remove_initial_lindep'
    lin_depend = 1 !Initialize to non-dependent
    call remove_initial_lindep(molcfg,Nb_new,Bvec_tmp,B_scr) 

   iturn = 0
   run_ortho_again = .false.

   iturn = iturn + 1
   !
   !Project out
   !
       if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'orthonormalize -> util_scriptPx'
       do i = 1,Nb_new
          if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Current vector is number = ', i
          call util_scriptPx('N',D,S,Bvec_tmp(i))
       enddo
  
   !
   ! NORMALIZE VECTOR NUMBER Ibvec
   !
       do ibvec = 1,Nb_new !index for current bvector 
         if (lin_depend(ibvec) /= 0) then
           if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'RSP_SOLVER: Orthonormalize -> normalize'
           call rsp_normalize(molcfg,iturn,ibvec,run_ortho_again,lin_depend(ibvec),Bvec_tmp(ibvec))
         endif
       enddo !ibvec
!Sonia and Stinne new

    if(.NOT. molcfg%solver%rsp_quiet)write(molcfg%lupri,*)'done normaizing'

! Add new vectors to file 
! ib counts how many "acceptable" vectors there are in total
!    
    no_of_new = 0
    do irx = 1,Nb_new
       if (lin_depend(irx) /= 0) then
         no_of_new = no_of_new + 1
       endif
    enddo
    call mem_alloc(bvecs,no_of_new)
    no_of_new = 0
    do irx = 1,Nb_new
       if (lin_depend(irx) /= 0) then
         no_of_new = no_of_new + 1
         call mat_init(bvecs(no_of_new),Bvec_tmp(irx)%nrow,Bvec_tmp(irx)%ncol)
         call mat_assign(bvecs(no_of_new),Bvec_tmp(irx))
       endif
    enddo
    if(.NOT. molcfg%solver%rsp_quiet)write(molcfg%lupri,*)'done inside orthonormalize'
!
!     Set NB_NEW to actual number of acceptable new trial vectors
!
     
    Nb_new = no_of_new
     
    DEALLOCATE(lin_depend)
    call mat_free(B_scr)
    call mat_free(b_k)

  end subroutine orthonormalizeMCD

  subroutine Projection(GD,D,S,ngd,freq,Bvecs,NrP)
!  Project the Gradient
!  Should be called before a symmetric orthonomalization
!  ajt Generalization to treat negative frequencies too
    implicit none
    integer,intent(in)         :: ngd   !number of gradients
    real(realk), intent(in)    :: freq(rsp_number_of_omegas)
    type(Matrix),intent(in)    :: D,S
    type(Matrix),intent(in)    :: bvecs(:)
    type(Matrix),intent(inout) :: gd(rsp_number_of_rhs) !output   
    integer, intent(in)        :: NrP
    real(realk) :: temps1
    type(Matrix) :: tempm1,tempm2
    integer :: ndim,i,k,igd
    ndim = Bvecs(1)%nrow
    call mat_init(tempm1,ndim,ndim)
    call mat_init(tempm2,ndim,ndim)

    !loop over bvecs to project out
    do i=1,NrP

       !if(.NOT. molcfg%solver%rsp_quiet)write(molcfg%lupri,*)'norm',mat_dotproduct(Bvecs(I),Bvecs(I))
       call ABCcommutator(ndim,Bvecs(i),D,S,tempm1)
       call mat_mul(tempm1,S,'n','n',1E0_realk,0E0_realk,tempm2)
       call mat_mul(S,tempm2,'n','n',2E0_realk,0E0_realk,tempm1)

       !if an equation has negative frequency, it is resonant with
       !de-excitation, and the de-excitation vector (trans) should be projected
       if (any(freq(1:ngd) < 0.0E0_realk)) &
          call mat_trans(tempm1,tempm2)

       !loop over gradients to project
       do igd=1,ngd
          !if excitation resonant, project out Bvec,
          !if de-excitation resonant, project out trans(Bvec)
          if (freq(igd) >= 0.0E0_realk) then
             temps1 = mat_dotproduct(Bvecs(i),GD(igd))
             call mat_daxpy(-temps1,tempm1,GD(igd))
          else
             temps1 = mat_trab(Bvecs(i),GD(igd))
             call mat_daxpy(-temps1,tempm2,GD(igd))
          endif
          !if done with prjection, ensure OO and UU blocks are zero
          if (i==NrP) call util_scriptPx('T',D,S,GD(igd))
       enddo
    enddo

    call mat_free(tempm1)
    call mat_free(tempm2)

  end subroutine

  subroutine solve_red_lineqMCD(molcfg,ndim_red,ngd,freq,red_X)
!
!  Solve reduced linear response problem in subspace
!  Now finally we consider one gradient/solution at a time
!
!  INPUT:
!  ndim_red: half size of reduced space
!  ngd:      number of right hand sides
!  freq:     the ngd frequencies
!
!  OUTPUT:
!  red_X:    the ngd reduced space solution vectors
! MCD modification Thomas
! Generalization to equations with negative frequency AndreasJT
    implicit none
    type(rsp_molcfg),intent(inout) :: molcfg
    integer, intent(in) :: ndim_red, ngd
    real(realk),intent(in) :: freq(:)
    real(realk),intent(inout) :: red_X(:,:)
    real(realk),allocatable :: E2(:,:), S2(:,:), RHS(:,:)
    integer,allocatable :: IPIV(:)
    integer :: igd, ierr, rdim, nump, i, j, k, l
    ierr=0
    !number of vectors to project out (degeneracy of (de-)excitation(s))
    nump = nmcdvec
    !dimension of redsp. matrices after out-projection of (de-)excitation(s)
    rdim = 2*ndim_red - nump

    allocate(E2(rdim,rdim))
    allocate(S2(rdim,rdim))
    allocate(RHS(rdim,ngd))
    allocate(IPIV(rdim))

    !loop over equations (gradients and frequencies)
    do igd=1, ngd

       !only print E2 and S2 once (for the first equation)
       if (igd==1.and.molcfg%solver%info_rsp) then
          write (molcfg%lupri,*) "UNMOD E2 :"
          call LS_OUTPUT(red_E(1:2*ndim_red,1:2*ndim_red), 1, 2*ndim_red, &
          &           1, 2*ndim_red, 2*ndim_red, 2*ndim_red, 1, molcfg%lupri)
          write (molcfg%lupri,*) "UNMOD S2 :"
          call LS_OUTPUT(red_S(1:2*ndim_red,1:2*ndim_red), 1, 2*ndim_red, &
          &           1, 2*ndim_red, 2*ndim_red, 2*ndim_red, 1, molcfg%lupri)
       endif

       if (molcfg%solver%info_rsp) then
          write (molcfg%lupri,*) "FREQUENCY IN SOLVE RED SPACE:",freq(igd)
          write (molcfg%lupri,*) "UNMOD E2 - omega S2:"
          call LS_OUTPUT(red_E(1:2*ndim_red,1:2*ndim_red) - freq(igd)* &
          &           red_S(1:2*ndim_red,1:2*ndim_red), 1, 2*ndim_red, &
          &           1, 2*ndim_red, 2*ndim_red, 2*ndim_red, 1, molcfg%lupri)
          write (molcfg%lupri,*) 'UNMOD RHS:'
          call LS_OUTPUT(red_GD(1:2*ndim_red,igd), 1, 2*ndim_red, 1, 1, 2*ndim_red, 1, 1, molcfg%lupri)
       endif

       !copy from red_E/S/GD into E2/S2/RHS, skipping resonant rows/cols
       do j = 1, rdim
          l = min(j+nump,merge(2*j-1,2*j,freq(igd)<0.0E0_realk))
          do i = 1, rdim
             k = min(i+nump,merge(2*i-1,2*i,freq(igd)<0.0E0_realk))
             E2(i,j) = red_E(k,l)
             S2(i,j) = red_S(k,l)
          enddo
          RHS(j,igd) = red_GD(l,igd)
       enddo

       if (molcfg%solver%info_rsp) then
          write (molcfg%lupri,*) "REDUCED E2, solve_red_lineq:"
          call LS_OUTPUT(E2, 1, rdim, 1, rdim, rdim, rdim, 1, molcfg%lupri)
          write (molcfg%lupri,*) "REDUCED S2, solve_red_lineq:"
          call LS_OUTPUT(S2, 1, rdim, 1, rdim, rdim, rdim, 1, molcfg%lupri)
          write (molcfg%lupri,*) 'REDUCED E2 - omega*S2:'
          call LS_OUTPUT(E2 - freq(igd)*S2, 1, rdim, 1, rdim, rdim, rdim, 1, molcfg%lupri)
          write (molcfg%lupri,*) 'REDUCED RHS:'
          call LS_OUTPUT(RHS(:,igd), 1, rdim, 1, 1, rdim, 1, 1, molcfg%lupri)
       endif

       !Solve set of linear equations Ax = b: (Solution vector is found in RHS)
       call DGESV(rdim, 1, E2-freq(igd)*S2, rdim, IPIV, RHS(:,igd), rdim, IERR)
       if (IERR /= 0) then
          WRITE(MOLCFG%LUPRI,'(/A, i4)') &
          &     'Problem in DGESV, IERR = ', IERR
          CALL LSQUIT(' Problem in DGESV',molcfg%lupri)
       endif

       if (molcfg%solver%info_rsp) then
          write (molcfg%lupri,*) 'REDUCED Solution vector - solve_red_lineq:'
          call LS_OUTPUT(RHS(:,igd), 1, rdim, 1, 1, rdim, 1, 1, molcfg%lupri)
       endif

       !copy from RHS over to red_X, filling in zeros for the resonant indices
       i = merge(2, 1, freq(igd) < 0.0E0_realk)
       red_X(i:i+2*(nump-1):2,igd) = 0.0E0_realk
       i = merge(1, 2, freq(igd) < 0.0E0_realk)
       red_X(i:i+2*(nump-1):2,igd) = RHS(1:nump,igd)
       red_X(2*nump+1:2*ndim_red,igd) = RHS(nump+1:,igd)

       if (molcfg%solver%info_rsp) then
          write (molcfg%lupri,*) 'Solution vector - MODIFIED - solve_red_lineq:'
          call LS_OUTPUT(red_X(1:2*ndim_red,igd), 1, 2*ndim_red, 1, 1, 2*ndim_red, 1, 1, molcfg%lupri)
       endif

    enddo

    deallocate(E2, S2)
    deallocate(RHS, IPIV)

  end subroutine solve_red_lineqMCD
 
subroutine get_rho(molcfg,D,S,nnew,bvecs,rhos)
!******************************************************
! Interface routine to carry out linear transformations
!   RHOS(I) = S[2]*N(I) (rho) > check omega
!******************************************************
    implicit none
    type(rsp_molcfg) :: molcfg 
    
!    type(decompItem),intent(in) :: decomp
    type(Matrix), intent(in) :: D,S,Bvecs(:)
    integer, intent(in) :: nnew
    type(Matrix), intent(inout) :: rhos(:)
    integer :: i,ndim 
    ndim = S%nrow

    do i = 1, nnew
      call make_rhos(molcfg,D,S,Bvecs(i),rhos(i))
    enddo
  end subroutine get_rho
!====================================================  
 subroutine make_rhos(molcfg,D,S,Bvecs_i,rho_i)
!*******************************************************
!  Carry out linear transformation
!  for one trial vector b(i). Loop on nr b's outside 
!*******************************************************
! Calculate the linear transformed vector rho of RSP
! rho = S[2]b = - S[b,D]_s S = - S(bSD-DSb)S = -SbSDS+SDSbS
! For one b at a time.
!************************************************************

    implicit none
    type(rsp_molcfg) :: molcfg 
    type(Matrix), intent(in) :: D,S,Bvecs_i
    type(Matrix), intent(inout) :: rho_i  !output
    type(matrix) :: prod,prod1,prod2
    integer :: ndim
  
    molcfg%solver%rsp_nlintra = molcfg%solver%rsp_nlintra + 1 !Count no of linear transformations
    ndim = S%nrow
    call mat_init(prod,ndim,ndim)
    call mat_init(prod1,ndim,ndim)
    call mat_init(prod2,ndim,ndim)

    call ABCcommutator(ndim,bvecs_i,D,S,prod2)
    call mat_mul(prod2,S,'n','n',1E0_realk,0E0_realk,prod)
    call mat_mul(S,prod,'n','n',-2E0_realk,0E0_realk,rho_i)
    call util_scriptPx('T',D,S,rho_i)

    !FREE matrices!
    call mat_free(prod)
    call mat_free(prod1)
    call mat_free(prod2)
  end subroutine make_rhos
!=================================================

subroutine real_solver_check(molcfg,F,D,S,gd,XR,omega)
!***************************************************************
!checking if real equation solved to the right solution
!***************************************************************
implicit none
     type(rsp_molcfg), intent(inout) :: molcfg
!    type(decompItem),intent(in) :: decomp
    type(Matrix), intent(in) :: F,D,S,gd
    type(Matrix)   :: XR
    type(Matrix)   :: E2XR,S2XR,res_real
    real(realk) :: gammma, a
    real(realk),intent(in) :: omega
    integer :: i,ndim
     logical  :: make_rhos 
    ndim = S%nrow

call mat_init(E2XR,ndim,ndim)
call mat_init(S2XR,ndim,ndim)
    
call make_lintran_vecs(molcfg,D,S,F,XR,E2XR,S2XR,.true.)
call mat_daxpy(-omega,S2XR,E2XR)
call mat_free(S2XR)

call mat_init(res_real,ndim,ndim)
call mat_add(1E0_realk,gd,-1E0_realk,E2XR,res_real)

a=mat_sqnorm2(res_real)

write(molcfg%lupri,*) 'sqnorm res_real:', a

call mat_free(res_real)
call mat_free(E2XR)

end subroutine real_solver_check

!> \brief Wrapper for calling the proper preconditioning routine
!> \author S. Host
!> \date October 2006
subroutine rsp_AB_precond(molcfg,Gn,S,omega,Gnt)
  implicit none
    !> Contains settings for response solver 
    type(rsp_molcfg), intent(inout) :: molcfg
    !> Matrix to be preconditioned
    type(Matrix), intent(in) :: Gn
    !> Overlap matrix
    type(Matrix), intent(in) :: S
    !> Frequency/eigenvalue
    real(realk), intent(inout)  :: omega
    !> Preconditioned matrix (output)
    type(Matrix), intent(inout) :: Gnt
    integer                     :: ndim
    type(Matrix)                :: GnU, GntU

     ndim = Gn%nrow

     if (molcfg%solver%rsp_MO_precond) then 
        if(.not. molcfg%solver%rsp_quiet)write (molcfg%lupri,*) 'MO preconditioning'
        call MO_precond(Gn,S,omega,Gnt,molcfg%decomp%nocc) 
     else if (molcfg%solver%rsp_no_precond) then
        !do nothing
        call mat_assign(Gnt,Gn)
     else 
        call mat_init(GnU,ndim,ndim)
        call mat_init(GntU,ndim,ndim)
        !Preconditioning of Gn by solving A Gnt = Gn
        !Convert residual to orthonormal basis:
        call res_to_oao_basis(molcfg%decomp,Gn, GnU)
        call oao_rsp_solver(molcfg%decomp, GnU, omega, GntU)
        !Convert back to non-orthonormal basis:
        call x_from_oao_basis(molcfg%decomp,GntU, Gnt)
        call mat_free(GnU)
        call mat_free(GntU)
     endif

end subroutine rsp_AB_precond

!> \brief Wrapper for calling the chosen starting guess for eigenvalue problem.
!> \author S. Host
!> \date October 2006
!>
!> If (molcfg%solver%rsp_mostart):
!> In the response (E[2]-omega*S[2])b = 0 case (Excitation energies), 
!> the initial trial vector b is found using information about orbital
!> energies if they exists 
!> 
!> otherwise: Use eigenvectors corresponding to lowest eigenvalues of FUP and FUQ
!> They are calculated in dd_utilities.f90
!> 
    subroutine get_1st_rsp_trials(molcfg,nroots,bvec_ao)
      implicit none
      !> Contains settings for response solver 
      type(rsp_molcfg), intent(inout) :: molcfg
      !> Requested number of excitation energies
      integer, intent(in) :: nroots
      !> The nroots starting guesses (output)
      type(matrix), intent(inout) :: bvec_ao(nroots)
      type(debugItem)          :: debug
      type(DDitem)             :: DD

      if (molcfg%solver%rsp_mostart) then
         write(molcfg%lupri,*) 'Use MO orbital energy differences to find &
                      & first guess for eigenvectors/values'
         call get_rsp_trials_from_MO(nroots,bvec_ao,molcfg%decomp%nocc)
      else 
         call dd_init(molcfg%decomp,DD)
         call DD_homolumo_and_heseigen(DD,molcfg%decomp,debug,.true.,nroots,rsp_iniguess=bvec_ao)
         call dd_shutdown(molcfg%decomp,DD)
         !write(molcfg%lupri,*) 'Exci starting guess:'
         !call mat_print(bvec_ao(1),1,bvec_ao(1)%nrow,1,bvec_ao(1)%nrow,molcfg%lupri)
      endif

    end subroutine get_1st_rsp_trials

!> \brief As get_1st_rsp_trials, but for unrestricted calculations.
!> \author S. Host
!> \date October 2006
    subroutine get_1st_rsp_trials_unres(molcfg,nroots,bvec_ao)
      implicit none
      !> Contains settings for response solver 
      type(rsp_molcfg), intent(inout) :: molcfg
      !> Requested number of excitation energies
      integer, intent(in) :: nroots
      !> The nroots starting guesses (output)
      type(matrix), intent(inout) :: bvec_ao(nroots)
      type(debugItem)          :: debug
      type(DDitem)             :: DD

      if (molcfg%solver%rsp_mostart) then
         STOP 'MO starting guess for excitation energies currently not implemented for open shell'
         write(molcfg%lupri,*) 'Use MO orbital energy differences to find &
                      & first guess for eigenvectors/values'
         call get_rsp_trials_from_MO(nroots,bvec_ao,molcfg%decomp%nocc)
      else 
         call dd_init(molcfg%decomp,DD)
         call DD_homolumo_and_heseigen_unres(DD,molcfg%decomp,debug,.true.,nroots,bvec_ao)
         call dd_shutdown(molcfg%decomp,DD)
      endif

    end subroutine get_1st_rsp_trials_unres

subroutine orthogonalizeDegenerate(molcfg,D,S,F,eivecs,eival,nexci)
  implicit none
  type(rsp_molcfg), intent(inout) :: molcfg
!  type(decompItem),intent(inout) :: decomp
  type(Matrix), intent(in)    :: D,S,F
  integer, intent(in)      :: nexci
  type(matrix),intent(inout) :: eivecs(nexci)
  real(realk), intent(in)  :: eival(nexci)
  type(matrix) :: rho1,rho2,sigma,newxsol(2)
  integer :: iindex,i,j,degeneratestates,state1(nexci),state2(nexci),kk
  integer :: jindex,ndim,k
  real(realk) :: freqi,freqj,TT,EMAT(2,2),SMAT(2,2),THRESHOLD
  real(realk) :: EMAT2(nexci,nexci),SMAT2(nexci,nexci)
  
  ndim = F%nrow
  CALL MAT_INIT(rho1,ndim,ndim)
  CALL MAT_INIT(rho2,ndim,ndim)
  CALL MAT_INIT(sigma,ndim,ndim)     
  CALL MAT_INIT(NEWXSOL(1),ndim,ndim)
  CALL MAT_INIT(NEWXSOL(2),ndim,ndim)

  iindex=0
  i=0
  degeneratestates=0
  freqi=10000
  j=1
  THRESHOLD = molcfg%solver%rsp_conv_thr*1E-2_realk
!  print*,'orthogonalizeDegenerate',THRESHOLD
  do while (j .LE. nexci)
     freqj=eival(j)
!     print*,'freqi',freqi
!     print*,'freqj',freqj
!     print*,'abs(freqi-freqj)',abs(freqi-freqj),'thr:',THRESHOLD
     if(abs(freqi-freqj)<THRESHOLD) then
        if(i==0) then !new excited state
           i=2
           degeneratestates=degeneratestates+1
           state1(degeneratestates)=iindex
           state2(degeneratestates)=j
!           print*,'stat1',state1(degeneratestates)
!           print*,'stat2',state2(degeneratestates)
!        else !degeneracy higher than 2
!           print*,'warning orthogonalizeDegenerate degeneracy in higher than 2'
        endif
        freqi=freqj
        iindex=j
     else
        !nondegenerate
        i=0
        freqi=freqj
        iindex=j
     endif
     j=j+1
  enddo
  
!  print*,'number of degenerate states',degeneratestates
  do KK=1,degeneratestates
     !test if a orthogonalization is required
     
     
     call mat_assign(NEWXSOL(1),eivecs(state1(KK)))
     call mat_assign(NEWXSOL(2),eivecs(state2(KK)))
     !ONLY NEED RHO
     call make_lintran_vecs(molcfg,D,S,F,eivecs(state1(KK)),sigma,rho1,.true.)
     TT = mat_dotproduct(eivecs(state2(KK)),rho1)
     call mat_daxpy(-TT,eivecs(state2(KK)),NEWXSOL(1))      
     
     call make_lintran_vecs(molcfg,D,S,F,newxsol(1),sigma,rho1,.true.)
     TT = 1E0_realk/sqrt(ABS(mat_dotproduct(newXSOL(1),rho1)))
     call mat_scal(TT,NEWXSOL(1))
     
     call mat_assign(eivecs(state1(KK)),NEWXSOL(1))
     call mat_assign(eivecs(state2(KK)),NEWXSOL(2))
     
     !--------TEST--------------
     !EMAT=0E0_realk
     !SMAT=0E0_realk
     !do K=1,2
     !   call make_lintran_vecs(decomp,D,S,F,newXsol(k),sigma,rho1,.true.)
     !   EMAT(k,k)=mat_dotproduct(newxsol(k),sigma)
     !   SMAT(k,k)=mat_dotproduct(newxsol(k),rho1)
     !   do i=1,2
     !      EMAT(i,k)=mat_dotproduct(newxsol(i),sigma)
     !      SMAT(i,k)=mat_dotproduct(newxsol(i),rho1)
     !   enddo
     !enddo
     !WRITE(molcfg%lupri,*)'QQQ THE TRANSFORMED bE[2]b matrix deg=',KK
     !WRITE (MOLCFG%LUPRI,'(10X,F12.8,2X,F12.8)')  EMAT(1,1),EMAT(1,2)
     !WRITE (MOLCFG%LUPRI,'(10X,F12.8,2X,F12.8)')  EMAT(2,1),EMAT(2,2)
     !WRITE(molcfg%lupri,*)' '
     !WRITE(molcfg%lupri,*)'    THE TRANSFORMED bS[2]b matrix '
     !WRITE (MOLCFG%LUPRI,'(10X,F12.8,2X,F12.8)')  SMAT(1,1),SMAT(1,2)
     !WRITE (MOLCFG%LUPRI,'(10X,F12.8,2X,F12.8)')  SMAT(2,1),SMAT(2,2)
  enddo
  
  !--------------------------------------------------------------------
     !EMAT2=0E0_realk
     !SMAT2=0E0_realk
     !do K=1,nexci
     !   call make_lintran_vecs(decomp,D,S,F,eivecs(k),sigma,rho1,.true.)
     !   do i=1,nexci
     !      EMAT2(i,k)=mat_dotproduct(eivecs(i),sigma)
     !      SMAT2(i,k)=mat_dotproduct(eivecs(i),rho1)
     !   enddo
     !enddo
     !WRITE(molcfg%lupri,*)'QQQ THE TRANSFORMED bE[2]b matrix'
     !call ls_output(EMAT2,1,nexci,1,nexci,nexci,nexci,1,lupri)
     !WRITE(molcfg%lupri,*)' '
     !WRITE(molcfg%lupri,*)'    THE TRANSFORMED bS[2]b matrix'
     !call ls_output(EMAT2,1,nexci,1,nexci,nexci,nexci,1,lupri)
!----------------------------------------------------------------------

   CALL MAT_free(rho1)
   CALL MAT_free(rho2)
   CALL MAT_free(sigma)     
   CALL MAT_free(NEWXSOL(1))
   CALL MAT_free(NEWXSOL(2))

 end subroutine orthogonalizeDegenerate
end module RSPsolver
