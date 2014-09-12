module COMPLEXSOLVER
  use decompMod
  use files
  use memory_handling
  use matrix_module
  use matrix_operations
  use RSPsolver
  use rsp_util
  use precision
  private
  public ::  rsp_complex_solver,complex_solver_check,complex_solver, &
          &  rsp_complex_init
  integer, save :: rsp_number_of_current_trial, rsp_number_of_rhs, rsp_number_of_sols, &
                 & rsp_number_of_omegas, rsp_number_of_startvecs, lusigma_rsp, lurho_rsp, lub_rsp, &
                 & rsp_bvec_dim, lu_x_rsp, lu_sigma_rsp, lu_rho_rsp
                 !to complex solver
  real(realk), allocatable, save :: red_E(:,:), red_S(:,:), red_GD(:,:), red_GDI(:,:)
  real(realk), allocatable, save :: red_E_glob(:,:), red_S_glob(:,:), red_GD_glob(:,:), red_GDI_glob(:,:)
  logical,save :: RSPOnMaster
contains
  subroutine rsp_complex_init(ntrial, nrhs, nsol, nomega, nstart)
  implicit none
      integer, intent(in) :: ntrial, nrhs, nsol, nomega, nstart
      !nstart only relevant for eigenvalue problem
      !nrhs only relevant for linear equations (always 1 for eigenvalue problem)

      rsp_number_of_current_trial = ntrial 
      rsp_number_of_rhs = nrhs
      rsp_number_of_sols = nsol
      rsp_number_of_omegas = nomega
      rsp_number_of_startvecs = nstart
      rsp_bvec_dim = max(rsp_number_of_current_trial,rsp_number_of_startvecs)+3
      RSPonMaster = .TRUE.
  end subroutine rsp_complex_init

 !======================================================================== 
  subroutine rsp_complex_solver(molcfg,F,D,S,ngd,GD,EIVAL,ifreq,XSOL,XSOLimg,gd_complex,gdi)
 !=======================================================================
 ! The linear scaling complex response solver
 ! solving a complex equation problem
 ! (E[2]-(w+igamma)S[2])(X[R]+iX[I])=GD[R]+iGD[I]
 ! (a) if gamma=0 and GD[I]=0, then solving the standard response eqs
 ! (b) solvng damped response eq.
 !
 ! INPUT FILE
 ! .COMPLEX (damping parameter - gamma)
 ! gamma given in input as: cfg_rsp_gamma
 ! .DAMP_FREQ (defining the frequency spectrum we want to scan)
 ! number_of_frequencies, minimal_frequency, maximal_requency
 !
 ! input:
 ! ngd    number of RHS
 ! GD     RHS
 ! EIVAL  frequency
 ! ifreq  no. of frequency
 !
 ! output:
 ! XSOL    real solution vector
 ! XSOLimg imaginary solution vector
 !
 ! Joanna, January 2009
 !=========================================================================
    use RSPsolver
    implicit none
!    type(decompItem),intent(inout) :: decomp
    type(rsp_molcfg), intent(inout) :: molcfg
    integer, intent(in)         :: ngd,ifreq 
    type(Matrix), intent(in)    :: D,S,F
    logical,intent(in)          :: gd_complex
    type(Matrix)                :: gd(rsp_number_of_rhs)
    type(Matrix), optional      :: gdi(rsp_number_of_rhs)
    real(realk), intent(inout)  :: eival(rsp_number_of_omegas)  !if LINEQ, eival = laser freqs  => input
    type(Matrix), intent(inout) :: XSOL(rsp_number_of_sols) !output solution vectors
    type(Matrix), intent(inout) :: XSOLimg(rsp_number_of_sols) !output img solution vectors
    logical                     :: conv
    integer                     :: i

write(molcfg%lupri,*) 'Entering the complex rsp solver, gamma=', molcfg%solver%rsp_gamma
if ((abs(molcfg%solver%rsp_gamma) .LT. 1E-8_realk) .and. (.not. gd_complex)) then
   IF(molcfg%solver%UseExcitationVecs)then
      call rsp_init(rsp_number_of_current_trial,rsp_number_of_rhs,&
           & rsp_number_of_sols,rsp_number_of_omegas,&
           & rsp_number_of_startvecs+molcfg%solver%rsp_eigenvecs)
   ELSE
      call rsp_init(rsp_number_of_current_trial,rsp_number_of_rhs,&
           & rsp_number_of_sols,rsp_number_of_omegas,rsp_number_of_startvecs)
   ENDIF
   call rsp_solver(molcfg,D,S,F,.true.,1,GD,EIVAL,XSOL)
   do i=1,rsp_number_of_sols
      call mat_zero(xsolimg(i))
   enddo
else
   call complex_solver(molcfg,F,D,S,1,GD,GDI,EIVAL,XSOL,XSOLimg,gd_complex)
   !     call complex_solver_check(F,D,S,gd(1),XSOLimg(1),XSOL(1),EIVAL(ifreq),&
   !        & gd_complex,conv,gdi(1))
endif
    
end subroutine rsp_complex_solver
!==========================================================
!==============================================================
  subroutine complex_solver(molcfg,F,D,S,ngd,GDB,GDI,EIVAL,XSOL,XSOLimg,gd_complex)
!==============================================================
 ! The linear scaling complex response solver
 ! solving a complex equation problem
 ! (E[2]-(w+igamma)S[2])(X[R]+iX[I])=GD[R]+iGD[I]
 !
 ! input:
 ! ngd    number of RHS
 ! GDB    RHS
 ! EIVAL  frequency
 !
 ! output:
 ! XSOL    real solution vector
 ! XSOLimg imaginary solution vector
 !   
 ! Joanna, December 2008
 !=========================================================================
!
implicit none
!type(decompItem),intent(inout)      :: decomp
type(rsp_molcfg), intent(inout)     :: molcfg
type(Matrix),intent(in)             :: F,D,S
type(Matrix),intent(inout)          :: GDB(rsp_number_of_rhs)
type(Matrix),intent(inout),optional :: GDI(rsp_number_of_rhs)
real(realk),intent(inout)           :: eival(rsp_number_of_omegas)
type(Matrix),intent(inout)          :: XSOL(rsp_number_of_sols),XSOLimg(rsp_number_of_sols)
logical, intent(in)                 :: gd_complex
integer, intent(in)                 :: ngd
type(Matrix)                        :: RHS_real(rsp_number_of_rhs),RHS_img(rsp_number_of_rhs)
type(Matrix)                        :: gdb_prec(rsp_number_of_rhs),gdbi_prec(rsp_number_of_rhs)
type(Matrix),pointer                :: sig_prec(:),rho_prec(:)
type(Matrix)                        :: Xsolvec(2*rsp_number_of_rhs),sigmas(2*rsp_number_of_rhs)
type(Matrix)                        :: rhos(2*rsp_number_of_rhs),Xsolvector(2)
real(realk)                         :: gammma,res_norm_tot,r1,r2,r3
integer                             :: i,j,ndim,k,n_red,nx_new,m,l,n_v,max_it,max_red
logical                             :: conv,make_rhos,conv_real,conv_img
real(realk),allocatable             :: red_X_glob(:,:),red_Ximg_glob(:,:)
real(realk)                         :: absorp, disper, rsp_thresh
 nullify(sig_prec)
 nullify(rho_prec)
 ndim = S%nrow
 gammma=molcfg%solver%rsp_gamma
 max_it=molcfg%solver%rsp_maxit
 max_red=molcfg%solver%rsp_maxred
 rsp_thresh=molcfg%solver%rsp_thresh
   lu_x_rsp= -1 ; lu_sigma_rsp = -1 ; lu_rho_rsp = -1
   CALL LSOPEN(lu_x_rsp,'xrsp','unknown','UNFORMATTED')
   CALL LSOPEN(lu_sigma_rsp,'sigma_rsp','unknown','UNFORMATTED')
   CALL LSOPEN(lu_rho_rsp,'rho_rsp','unknown','UNFORMATTED')
   allocate(red_E_glob(2*max_red,2*max_red), red_S_glob(2*max_red,2*max_red))
   allocate(red_X_glob(2*max_red,ngd),red_Ximg_glob(2*max_red,ngd))
   ALLOCATE(RED_GD_glob(2*max_red,ngd))
   if (gd_complex) then
       allocate(red_gdi_glob(2*max_red,ngd))
       red_gdi_glob=0E0_realk
   endif
   red_GD_glob = 0.0E0_realk
   red_E_glob = 0.0E0_realk ; red_S_glob = 0.0E0_realk ; red_X_glob = 0.0E0_realk ; red_Ximg_glob = 0.0E0_realk
  
   do i=1,ngd
      call mat_init(RHS_real(i),ndim,ndim)
      call mat_init(RHS_img(i),ndim,ndim)
      call mat_zero(RHS_img(i))
      call mat_zero(RHS_real(i))
   enddo
   do i=1,ngd
      call mat_init(gdb_prec(i),ndim,ndim)
      call mat_zero(gdb_prec(i))
      if (gd_complex) then
         call mat_init(gdbi_prec(i),ndim,ndim)
         call mat_zero(gdbi_prec(i))
      endif
       
      call mat_assign(RHS_real(i),GDB(i))
      if (gd_complex) then
          call mat_assign(RHS_img(i),GDI(i))
      endif
      IF(molcfg%solver%rsp_MO_precond) then
         call MO_precond_complex(gdb(i),D,S,gdb_prec(i),molcfg%decomp%nocc)
         if (gd_complex) then
            call MO_precond_complex(gdi(i),D,S,gdbi_prec(i),molcfg%decomp%nocc)
         endif
      ELSEIF(molcfg%solver%rsp_no_precond) then
         call mat_assign(gdb_prec(i),gdb(i))
         if (gd_complex) then
            call mat_assign(gdbi_prec(i),gdi(i))
         endif
      ELSE
         call lsquit('complex AO precond not implemented ',-1)
      ENDIF
   enddo
   n_red=0
   conv_real=.false.
   conv_img=.false.
   conv=.false.
   nx_new=0

   do j=1,2
      call mat_init(xsolvector(j),ndim,ndim)
      call mat_zero(xsolvector(j))
   enddo

   ! get the start vectors
   ! if (molcfg%solver%rsp_damp_2start) (.TWOSTART) start vectors obtained by solving 2 sets of std response eqs
   ! default -> start vectors as preconditioned RHSs.
   call get_start_vectors(molcfg,F,D,S,RHS_real,RHS_img,ngd,eival,nx_new,xsolvector)
  
    do i = 1,max_it
     write(molcfg%lupri,*) '------------------'
     write(molcfg%lupri,*) 'Start macroiteration:',i
     write(molcfg%lupri,*) '------------------'
     n_v=nx_new 
     do j=1,nx_new
        call mat_init(xsolvec(j),ndim,ndim)
        call mat_init(sigmas(j),ndim,ndim)
        call mat_init(rhos(j),ndim,ndim)
!        call mat_init(sig_prec(j),ndim,ndim)
!        call mat_init(rho_prec(j),ndim,ndim)
      enddo
      
      call orthonormalize22(molcfg,D,S,xsolvector,nx_new,n_red,xsolvec)
      call transform_vectors(molcfg,D,S,F,nx_new,xsolvec,sig_prec,rho_prec,.true.)
      IF(molcfg%solver%rsp_MO_precond) then
         do j=1,nx_new
            call MO_precond_complex(sig_prec(j),D,S,sigmas(j),molcfg%decomp%nocc)
            call MO_precond_complex(rho_prec(j),D,S,rhos(j),molcfg%decomp%nocc)
         enddo
      ELSEIF(molcfg%solver%rsp_no_precond) then
         do j=1,nx_new
            call mat_assign(sigmas(j),sig_prec(j))
            call mat_assign(rhos(j),rho_prec(j))
         enddo
      ELSE
         call lsquit('complex AO precond not implemented ',-1)
      ENDIF
      print*,'sig_prec(j) FREE'
      do j=1,nx_new
         call mat_free(sig_prec(j))
         call mat_free(rho_prec(j))
      enddo
      call mem_dealloc(sig_prec)
      call mem_dealloc(rho_prec)

      do j=1,2
         call mat_free(xsolvector(j))
      enddo
        
      if (gd_complex) then
        !build a reduce space
         call extend_complex_reduced_matrices(molcfg,n_red,nx_new,ngd,gdb_prec,xsolvec,sigmas,rhos,.true.,gdbi_prec)
        
          !solve complex reduced equation
         call solve_complex(molcfg,ndim,n_red,ngd,eival,red_X_glob,red_Ximg_glob,.true.)
    
        ! check for convergece and get new trialvectors
         do j=1,2
             call mat_init(xsolvector(j),ndim,ndim)
             call mat_zero(xsolvector(j))
         enddo
         call get_complex_res(molcfg,F,D,S,i,n_red,ngd,gdb_prec,red_X_glob,red_Ximg_glob,eival, &
        &    xsolvector,nx_new,conv,conv_real,conv_img,res_norm_tot,.true.,gdbi_prec)

      else
        !build a reduce space
           call extend_complex_reduced_matrices(molcfg,n_red,nx_new,ngd,gdb_prec,xsolvec,sigmas,rhos,.false.)

        !solve complex reduced equation
            call solve_complex(molcfg,ndim,n_red,ngd,eival,red_X_glob,red_Ximg_glob,.false.)

        ! check for convergece and get new trialvectors
        do j=1,2
           call mat_init(xsolvector(j),ndim,ndim)
           call mat_zero(xsolvector(j))
        enddo
        call get_complex_res(molcfg,F,D,S,i,n_red,ngd,gdb_prec,red_X_glob,red_Ximg_glob,eival, &
        &    xsolvector,nx_new,conv,conv_real,conv_img,res_norm_tot,.false.)
    endif

    if (conv)  then
        do l=1,ngd
           rewind(lu_x_rsp)
           call mat_zero(XSOL(l))
           call expand_on_basis(molcfg,n_red,RED_X_glob(1:n_red*2,l),lu_x_rsp,XSOL(l))
           rewind(lu_x_rsp)
           call mat_zero(XSOLimg(l))
           call expand_on_basis(molcfg,n_red,RED_Ximg_glob(1:n_red*2,l),lu_x_rsp,XSOLimg(l))
        enddo
        do j=1,n_v
           call mat_free(xsolvec(j))
           call mat_free(sigmas(j))
           call mat_free(rhos(j))
        enddo

           exit
    else         
        do l=1,ngd
           rewind(lu_x_rsp)
           call mat_zero(XSOL(l))
           call expand_on_basis(molcfg,n_red,RED_X_glob(1:n_red*2,l),lu_x_rsp,XSOL(l))
           absorp=mat_dotproduct(gdb(l),xsol(l))
           write(molcfg%lupri,*) "Dispersion value macroiteration:",i, "for &
               & frequency",eival(l),"=", absorp
      
           rewind(lu_x_rsp)
           call mat_zero(XSOLimg(l))
           call expand_on_basis(molcfg,n_red,RED_Ximg_glob(1:n_red*2,l),lu_x_rsp,XSOLimg(l))
           disper=mat_dotproduct(gdb(l),xsolimg(l))
           write(molcfg%lupri,*) "Absorption value in macroiteration:", i, "for &
               & frequency ", eival(l),"=", disper
         enddo

         if (i==1) then
             r1=res_norm_tot
          elseif (i==2) then
             r2=res_norm_tot
          elseif (i==3) then
             r3=res_norm_tot
             if ((abs(r1-r2)<rsp_thresh) .and. (abs(r2-r3)<rsp_thresh)) then
                do j=1,n_v
                   call mat_free(xsolvec(j))
                   call mat_free(sigmas(j))
                   call mat_free(rhos(j))
                 enddo

                    exit
              endif    
           else
              r1=r2
              r2=r3
              r3=res_norm_tot
              if ((abs(r1-r2)<rsp_thresh) .and. (abs(r2-r3)<rsp_thresh)) then
                 do j=1,n_v
                    call mat_free(xsolvec(j))
                    call mat_free(sigmas(j))
                    call mat_free(rhos(j))
                 enddo

                    exit
                endif    
            endif
      endif
      do j=1,n_v
          call mat_free(xsolvec(j))
          call mat_free(sigmas(j))
          call mat_free(rhos(j))
      enddo
   enddo
   do j=1,2
      call mat_free(xsolvector(j))
   enddo

   do i=1,ngd
      call mat_free(RHS_real(i))
      call mat_free(RHS_img(i))
   enddo
   do i=1,ngd
      call mat_free(gdb_prec(i))
      if (gd_complex) then
         call mat_free(gdbi_prec(i))
      endif
    enddo
    
    deallocate(red_E_glob,red_S_glob,red_gd_glob,red_X_glob,red_Ximg_glob)
    if (gd_complex) deallocate(red_gdi_glob)
    
    CALL LSCLOSE(lu_sigma_rsp,'DELETE')
    CALL LSCLOSE(lu_rho_rsp,'DELETE')
    CALL LSCLOSE(lu_x_rsp,'DELETE')
      
!================================================================  
end subroutine complex_solver
!================================================================  
!================================================================  
subroutine extend_complex_reduced_matrices(molcfg,n_red,nx_new,ngd,gd,xsolvec,sigmas,rhos,gd_complex,gdi)
!================================================================  
! Purpose: extend global reduced space with new vectors
!          The reduced E2, S2,G matrices (red_E_glob, red_S_glob, red_G_glob) are global variables
!          (see top of file) 
!
! INPUT: 
! n_red, number of vectors on disk and current half size of reduced space
! nx_NEW, number of new vectors
! ngd, the number of response vectors to be found
! gd(ngd), matrix array containing ngd right hand sides (gradients)
! sigmas,rhos,xsolvec:  matrix arrays of the nb_new new vectors
!
! OUTPUT:
! n_red, updated number of vecs on disk and new half size of reduced space
!
!================================================================  
implicit none
type(rsp_molcfg), intent(inout)  :: molcfg
integer,intent(inout)            :: n_red
integer, intent(in)              :: nx_new,ngd
type(Matrix),intent(in)          :: xsolvec(:),sigmas(:),rhos(:)
type(Matrix),intent(in)          :: gd(:)
type(Matrix),intent(in),optional :: gdi(:)
logical,intent(in)               :: gd_complex
type(Matrix)                     :: xT,sigmaT,rhoT,x_j,rho_j,sigma_j,x_jT
integer                          :: ndim,k,j,i,l,max_red

max_red=molcfg%solver%rsp_maxred

if ((N_red + Nx_new) > max_red) then
   WRITE (molcfg%lupri,*) 'ERROR in extend_complex_reduced_matrices: N_RED + N_NEW > rsp_maxred'
   WRITE (molcfg%lupri,*) 'N_RED ,Nx_NEW  =',N_red,Nx_NEW
   WRITE (molcfg%lupri,*) 'rsp_maxred =',max_red
   WRITE (molcfg%lupri,*) 'Reduce problem size or recompile orthonormalize '// &
        &                  'with larger rsp_maxred parameter'
   WRITE (*,*) 'ERROR in extend_complex_reduced_matrices: N_RED + N_NEW > rsp_maxred'
   WRITE (*,*) 'N_RED ,Nx_NEW  =',N_red,Nx_NEW
   WRITE (*,*) 'rsp_maxred =',max_red
   WRITE (*,*) 'Reduce problem size or recompile orthonormalize '// &
        &                  'with larger rsp_maxred parameter'
   CALL lsQUIT('complex orthonormalize error: N_RED + Nx_NEW > rsp_maxred parameter',molcfg%lupri)
END IF

ndim = gd(1)%nrow
call mat_init(xT,ndim,ndim)
call mat_init(sigmaT,ndim,ndim)
call mat_init(rhoT,ndim,ndim)
call mat_init(sigma_j,ndim,ndim)
call mat_init(rho_j,ndim,ndim)
call mat_init(x_j,ndim,ndim)
call mat_init(x_jT,ndim,ndim)

    k = 0
    do i = 2*n_red+1, 2*n_red+2*nx_new, 2  !Run over 2 elements at a time because of pairing
       k = k + 1
       l = 0
       do j = 2*n_red+1, 2*n_red+2*nx_new, 2
          l = l + 1
          call mat_trans(xsolvec(k), xT) 
          call mat_trans(sigmas(l), sigmaT)
          call mat_trans(rhos(l), rhoT)
          red_E_glob(i,j)     = mat_dotproduct(xsolvec(k),sigmas(l))
          red_E_glob(i+1,j)   = mat_dotproduct(xT,sigmas(l))
          red_E_glob(i,j+1)   = mat_dotproduct(xsolvec(k),sigmaT)
          red_E_glob(i+1,j+1) = mat_dotproduct(xT,sigmaT)

          red_S_glob(i,j)     = mat_dotproduct(xsolvec(k),rhos(l))
          red_S_glob(i+1,j)   = mat_dotproduct(xT,rhos(l))
          red_S_glob(i,j+1)   = -mat_dotproduct(xsolvec(k),rhoT)
          red_S_glob(i+1,j+1) = -mat_dotproduct(xT,rhoT)
       enddo      
    enddo

    !Setup lower half of E2 and S2:
    rewind(lu_sigma_rsp) ; rewind(lu_rho_rsp)
    do j = 1, 2*n_red, 2
       call mat_read_from_disk(lu_sigma_rsp,sigma_j,RSPonMaster)
       call mat_read_from_disk(lu_rho_rsp,rho_j,RSPonMaster)
       call mat_trans(sigma_j, sigmaT) 
       call mat_trans(rho_j, rhoT)
       k = 0
       do i = 2*n_red+1, 2*n_red+2*nx_new, 2  !Run over 2 elements at a time because of pairing
          k = k + 1
          call mat_trans(xsolvec(k), xT)
          red_E_glob(i,j)     = mat_dotproduct(xsolvec(k),sigma_j) 
          red_E_glob(i+1,j)   = mat_dotproduct(xT,sigma_j)
          red_E_glob(i,j+1)   = mat_dotproduct(xsolvec(k),sigmaT)
          red_E_glob(i+1,j+1) = mat_dotproduct(xT,sigmaT)

          red_S_glob(i,j)     = mat_dotproduct(xsolvec(k),rho_j)
          red_S_glob(i+1,j)   = mat_dotproduct(xT,rho_j)
          red_S_glob(i,j+1)   = -mat_dotproduct(xsolvec(k),rhoT)
          red_S_glob(i+1,j+1) = -mat_dotproduct(xT,rhoT)
       enddo
    enddo

    rewind(lu_x_rsp)
    !Explicitly calculate upper half of E2 and S2:
    do i = 1, 2*n_red, 2
       call mat_read_from_disk(lu_x_rsp,x_j,RSPonMASTER)
       call mat_trans(x_j, x_jT)
       k = 0
       do j = 2*n_red+1, 2*n_red+2*nx_new, 2  !Run over 2 elements at a time because of pairing
          k = k + 1
          call mat_trans(xsolvec(k), xT)
          call mat_trans(sigmas(k), sigmaT) 
          call mat_trans(rhos(k), rhoT)

          red_E_glob(i,j)     = mat_dotproduct(x_j,sigmas(k)) 
          red_E_glob(i+1,j)   = mat_dotproduct(x_jT,sigmas(k))
          red_E_glob(i,j+1)   = mat_dotproduct(x_j,sigmaT)
          red_E_glob(i+1,j+1) = mat_dotproduct(x_jT,sigmaT)

          red_S_glob(i,j)     = mat_dotproduct(x_j,rhos(k)) 
          red_S_glob(i+1,j)   = mat_dotproduct(x_jT,rhos(k))
          red_S_glob(i,j+1)   = -mat_dotproduct(x_j,rhoT)
          red_S_glob(i+1,j+1) = -mat_dotproduct(x_jT,rhoT)
       enddo
    enddo
   
    do j = 1, ngd  !Loop over number of gradients
       k=0
       do i = 2*n_red+1, 2*n_red+2*nx_new, 2 !Loop in steps of 2 because of pairing
          k = k+1
          call mat_trans(xsolvec(k),xT)
           red_gd_glob(i,j) = mat_dotproduct(xsolvec(k),gd(j))
           red_gd_glob(i+1,j) = mat_dotproduct(xT,gd(j))
           if (gd_complex) then
              red_gdi_glob(i,j) = mat_dotproduct(xsolvec(k),gdi(j))
              red_gdi_glob(i+1,j) = mat_dotproduct(xT,gdi(j))
           endif 
       enddo
    enddo

        !save on disk
    do i = 1,nx_new
       call mat_write_to_disk(lu_x_rsp,xsolvec(i),RSPonMaster)
       call mat_write_to_disk(lu_sigma_rsp,sigmas(i),RSPonMaster)
       call mat_write_to_disk(lu_rho_rsp,rhos(i),RSPonMaster)
       if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Extend_red_matrices2: set of sigma, rho, and bvector is written to disk'
    enddo
    n_red = n_red + nx_new !update reduced space dimension
    !Remember that the row-dimension of the reduced matrices 
    !is 2*ndim_red, where ndim_red is the number of trials

    if (MOLCFG%SOLVER%INFO_RSP_REDSPACE) then

       write (molcfg%lupri,*) 'Reduced E2, extend_red_matrices2:'
       call LS_OUTPUT(red_E_glob, 1, 2*n_red, 1, 2*n_red, 2*max_red, 2*max_red, 1, molcfg%lupri)

       write (molcfg%lupri,*) 'Reduced S2, extend_red_matrices2:'
       call LS_OUTPUT(red_S_glob, 1, 2*n_red, 1, 2*n_red, 2*max_red, 2*max_red, 1, molcfg%lupri)
    endif

    if (molcfg%solver%info_rsp) then
            write(molcfg%lupri,*) 'extend_gradients: after gradient extension'; 
            call LS_OUTPUT(red_gd_glob, 1, 2*n_red, 1, ngd, 2*max_red, 2*molcfg%solver%rsp_maxgd, 1, molcfg%lupri)
    endif

call mat_free(xT)
call mat_free(x_jT)
call mat_free(x_j)
call mat_free(rho_j)
call mat_free(sigma_j)
call mat_free(rhoT)
call mat_free(sigmaT)
    end subroutine extend_complex_reduced_matrices

!================================================================  

  subroutine orthonormalize22(molcfg,D,S,Bvec_tmp,Nb_new,Nb_prev,bvecs)
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
!**************************************************************************
    implicit none
    type(rsp_molcfg), intent(inout) :: molcfg
    type(Matrix),intent(in)    :: D,S
    integer, intent(in)        :: Nb_prev
    integer, intent(inout)     :: Nb_new
    type(Matrix),intent(inout) :: Bvec_tmp(:), bvecs(:)
!local 
    integer                    :: irx,i,iturn,ndim,k,ibvec,jbvec,jrx,no_of_new,lub_rsp_vec,dummy_i,max_red
    integer,allocatable        :: lin_depend(:) !linear dependency index array
    type(matrix)               :: B_scr, b_k, Xf, rho_k
    type(matrix)               :: orthovec !Will properly be changed
    real(realk)                :: TT,T1,T2,dummy_real
    logical                    :: run_ortho_again
     

max_red=molcfg%solver%rsp_maxred
    if ((Nb_prev + Nb_new) > max_red) then
      WRITE (molcfg%lupri,*) 'ERROR in orthonormalize: NB_PREV + NB_NEW > rsp_maxred'
      WRITE (molcfg%lupri,*) 'NB_PREV ,NB_NEW  =',NB_PREV,NB_NEW
      WRITE (molcfg%lupri,*) 'rsp_maxred =',max_red
      WRITE (molcfg%lupri,*) 'Reduce problem size or recompile orthonormalize '// &
      &                  'with larger rsp_maxred parameter'
      CALL lsQUIT('orthonormalize error: NB_PREV + NB_NEW > rsp_maxred parameter',molcfg%lupri)
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
   
       rewind(lu_x_rsp)
       do k = 1,Nb_prev
                  call mat_read_from_disk(lu_x_rsp,b_k,RSPOnMaster)
          do irx = 1,Nb_new
             !if lin_depend(irx) == 0, the vector is skipped
             !because of linear dependencies
             if (lin_depend(irx) /= 0) then
                 !Orthogonalize to (Z Y)_old
                 TT = mat_dotproduct(b_k,Bvec_tmp(irx))
                 call mat_daxpy(-TT,b_k,Bvec_tmp(irx))
                 !Orthogonalize to (Y Z)_old
                 !find the paired
                 call mat_trans(b_k,b_scr)
                 TT = mat_dotproduct(B_scr,Bvec_tmp(irx))
                 call mat_daxpy(-TT,B_scr,Bvec_tmp(irx))
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
   !
               call mat_daxpy(TT,Bvec_tmp(jbvec),Bvec_tmp(ibvec))
   !
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
             WRITE(molcfg%lupri,'(/A)') &
             &     'Error: Already ran twice through orthonormalize!'
             CALL lsQUIT('Error: Already ran twice through orthonormalize!',molcfg%lupri)
          else
             cycle
          endif
       else
          exit
       endif
   enddo !Maybe do twice


!Sonia and Stinne new
!
! Add new vectors to file 
! ib counts how many "acceptable" vectors there are in total
!    
    no_of_new = 0
    do irx = 1,Nb_new
       if (lin_depend(irx) /= 0) then
         no_of_new = no_of_new + 1
         call mat_assign(bvecs(no_of_new),Bvec_tmp(irx))
       endif
    enddo
!
!     Set NB_NEW to actual number of acceptable new trial vectors
!
     
     Nb_new = no_of_new
     
    DEALLOCATE(lin_depend)
    call mat_free(B_scr)
    call mat_free(b_k)

  end subroutine orthonormalize22
 
 !================================================================  
subroutine complex_solver_check(molcfg,F,D,S,gd,XI,XR,omega,gd_complex,conv,gdi)
!**************************************************************
!Purpose:
!Checking if complex equation converged to the right solution
!*************************************************************
implicit none
!    type(decompItem),intent(in) :: decomp
    type(rsp_molcfg), intent(inout) :: molcfg
    type(Matrix), intent(in)        :: F,D,S,gd
    type(Matrix)                    :: XR,XI
    type(Matrix),intent(in),optional:: gdi
    type(Matrix)                    :: E2XR,S2XR,E2XI, S2XI,res_real,res_img
    real(realk)                     :: gammma, a, b, valthresh
    real(realk),intent(in)          :: omega
    integer                         :: i,ndim
    logical,intent(in)              :: gd_complex
    logical,intent(out)             :: conv
    ndim = S%nrow

    valthresh=1e-5
    gammma=molcfg%solver%rsp_gamma
    
call mat_init(E2XI,ndim,ndim)
call mat_init(S2XR,ndim,ndim)
call mat_init(E2XR,ndim,ndim)
call mat_init(S2XI,ndim,ndim)
call mat_init(res_real,ndim,ndim)
call mat_init(res_img,ndim,ndim)
 
call make_lintran_vecs(molcfg,D,S,F,XR,E2XR,S2XR,.true.)
call make_lintran_vecs(molcfg,D,S,F,XI,E2XI,S2XI,.true.)

call mat_daxpy(-omega,S2XR,E2XR)
call mat_daxpy(-omega,S2XI,E2XI)
call mat_scal(-gammma,S2XI)
call mat_scal(gammma,S2XR)
call mat_daxpy(1E0_realk,gd,S2XI)
if (gd_complex) then
call mat_daxpy(1E0_realk,gdi,S2XR)
endif

call mat_add(1E0_realk,S2XI,-1E0_realk,E2XR,res_real)
call mat_add(1E0_realk,S2XR,-1E0_realk,E2XI,res_img)

a=mat_sqnorm2(res_real)
b=mat_sqnorm2(res_img)

write(molcfg%lupri,*) 'Sqnorm of res_real:', a, 'and res_img:', b

if ((a>valthresh) .or. (b>valthresh)) then
   WRITE(molcfg%lupri,"('WARNING: it converged to a wrong solution!')")
   conv=.false.
else
   conv=.true.
endif

call mat_free(E2XR)
call mat_free(S2XR)
call mat_free(E2XI)
call mat_free(S2XI)
call mat_free(res_real)
call mat_free(res_img)

end subroutine complex_solver_check
  

subroutine solve_complex(molcfg,ndim,ndim_red,ngd,freq,red_X,red_Xi,gd_complex)
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
!  red_X:    the ngd reduced space real solution vectors
!  red_X:    the ngd reduced space real solution vectors
    implicit none
    type(rsp_molcfg), intent(inout) :: molcfg
    integer, intent(in)        :: ndim_red, ngd,ndim
    real(realk),intent(in)     :: freq(:)
    logical                    :: gd_complex
    real(realk),intent(inout)  :: red_X(:,:)
    real(realk),intent(inout)  :: red_Xi(:,:)
    real(realk),allocatable    :: E2(:,:), S2(:,:)
    complex(realk),allocatable :: A(:,:),KHS(:)
    integer,allocatable        :: IPIV(:)
    integer                    :: igd, ierr, i, j, k
    integer                    :: lda, Nb,nmax
    real(realk)                :: gammma

    gammma=molcfg%solver%rsp_gamma
    Nb=64
    nmax=4*ndim_red
    lda=nmax

    allocate(E2(2*ndim_red,2*ndim_red), S2(2*ndim_red,2*ndim_red))
    allocate(KHS(nmax), IPIV(nmax))
    allocate(A(lda,nmax))
    E2=0E0_realk
    S2=0E0_realk
    KHS=0E0_realk
    IPIV=0E0_realk
    A=0E0_realk

    !Setup reduced E2, S2, and right hand side with proper dimension
    E2 = red_E_glob(1:2*ndim_red,1:2*ndim_red)
    S2 = red_S_glob(1:2*ndim_red,1:2*ndim_red)

    do igd = 1,ngd
       if (gd_complex) then
           KHS(1:2*ndim_red)= CMPLX(red_GD_glob(1:2*ndim_red,igd),red_GDI_glob(1:2*ndim_red,igd))
       else        
           KHS(1:2*ndim_red) = CMPLX(red_GD_glob(1:2*ndim_red,igd),0E0_realk)
       endif
       if (molcfg%solver%info_rsp) then
          write (molcfg%lupri,*) "E2, solve_red_lineq2:"
          call LS_OUTPUT(E2, 1, 2*ndim_red, 1, 2*ndim_red, 2*ndim_red, 2*ndim_red, 1, molcfg%lupri)
          write (molcfg%lupri,*) "S2, solve_red_lineq2:"
          call LS_OUTPUT(S2, 1, 2*ndim_red, 1, 2*ndim_red, 2*ndim_red, 2*ndim_red, 1, molcfg%lupri)
       endif
       E2 = E2 - freq(igd)*S2   
       S2 = -gammma*S2

       do j=1,2*ndim_red
          do k=1,2*ndim_red
              A(j,k)=CMPLX(E2(j,k),S2(j,k))
          enddo
       enddo
       
       if (molcfg%solver%info_rsp) then
          write(molcfg%lupri,*) 'A:'
          write(molcfg%lupri,*) ((A(i,j),j=1,2*ndim_red),i=1,2*ndim_red)  
          write(molcfg%lupri,*) 'KHS:'
          write(molcfg%lupri,*) KHS(1:2*ndim_red)
       endif

       !Solve set of linear equations Ax = b:
       call  zgesv(2*ndim_red,1,A,lda,IPIV,KHS,2*ndim_red,IERR)
       !Solution vector is found in RHS.
       if (IERR /= 0) then
          WRITE(molcfg%lupri,'(/A, i4)') &
          &     'Problem in ZSYSV, IERR = ', IERR
          CALL lsQUIT(' Problem in ZSYSV',molcfg%lupri)
       endif
       
       red_X(1:2*ndim_red,igd) = real(KHS(1:2*ndim_red))
       red_Xi(1:2*ndim_red,igd) = aimag(KHS(1:2*ndim_red))
    enddo

    deallocate(E2, S2)
    deallocate(KHS, A)
    deallocate(IPIV)
  end subroutine solve_complex
  
  
  subroutine get_complex_res(molcfg,F,D,S,itmic,ndim_red,ngd,gd,red_X_glob,red_Ximg_glob,eival, &
             &    xsolvector,nx_new,conv,conv_root,conv_root_img,res_norm_tot,gd_complex,gdi)
    implicit none
!******************************************************************
! Purpose: 
! 1)Construct residuals:
!   R^R=(E(2)-W(I)*S(2))*X^R(I) - GD^R + gamma*S[2]X^I
!   R^I=(E(2)-W(I)*S(2))*X^I(I) - GD^I - gamma*S[2]X^R
!   for NGD rsp vectors X(I) of reduced rsp problem
! 2)Test for convergence of ngd eigenvectors,
!   Convergence criterium:
!   ||R^R|| < cfg_rsp_conv_thr
!   ||R^I|| < cfg_rsp_conv_thr 
! 3) contruct new trial vectors as preconditioned residuals
!
!  INPUT:
!  itmic:    number of the microiteration (for printout)
!  ndim_red: half size of reduced space
!  ngd:      number of right hand sides
!  gd(ngd):  matrix array containing ngd right hand sides (gradients)
!  red_X_glob:    the ngd reduced space real solution vectors (same as RED_X)
!  red_Ximg_glob:    the ngd reduced space img solution vectors (same as RED_Xi)
!  eival:    the ngd frequencies
!  

!  OUTPUT:
!  nx_NEW:   number of new vectors
!  conv:     set true if converged, and microiterations should be stopped
!  xsolvector(nx_new) new trial vectors
!******************************************************************

    type(rsp_molcfg), intent(inout)  :: molcfg
    integer,intent(in)               :: itmic,ndim_red,ngd
    type(Matrix),intent(in)          :: F,D,S
    type(Matrix),intent(in)          :: gd(:)
    type(Matrix),intent(in),optional :: gdi(:)
    logical,intent(in)               :: gd_complex
    real(realk), intent(in)          :: red_X_glob(:,:),eival(:)
    real(realk), intent(in)          :: red_Ximg_glob(:,:)
    type(Matrix),intent(inout)       :: xsolvector(2)
    logical, intent(out)             :: conv
    real(realk),intent(out)          :: res_norm_tot
    logical,intent(out)              :: conv_root, conv_root_img
    integer, intent(out)             :: nx_new
!local
    type(Matrix)                     :: S2XR(ngd), S2XI(ngd)
    type(Matrix)                     :: residuals(ngd),residuals_img(ngd),rhs(ngd),rhsi(ngd),residuals2(2*ngd)
    real(realk)                      :: wibx,res_norm,x_red_norm,conv_test,ddot,diff,x_red_norm_img
    integer                          :: i,n_not_conv,ndim_red_mat,ndim,index_conv(ngd),ibx,j,k,max_red
    real(realk)                      :: a, b,c,rsp_thresh
    real(realk)                      :: av_norm,res_norm_img, conv_test_img,gammma

    gammma=molcfg%solver%rsp_gamma
    max_red=molcfg%solver%rsp_maxred

    ndim = S%nrow
    ndim_red_mat = 2*ndim_red
    n_not_conv = 0
    conv=.false.
 
    do i=1,ngd
       call mat_init(S2XR(i),ndim,ndim)
       call mat_init(S2XI(i),ndim,ndim)
       call mat_zero(S2XR(i))
       call mat_zero(S2XI(i))
    enddo
!
! find E[2]*X(i)-w(i) S[2]*X(i) as linear combination of 
! previous Sigmas and RHOs. 
! X(i) = sum_j Bvecold_j REDX_j(i)
! Sigmanew(i) = E[2]X(i) = sum_j REDX_j(i) E[2]Bvecold_j =
!             = sum_j REDX_j(i) Sigmaold_j
! and add to previous -w_i S[2]*X(I) vector
! 
    do ibx = 1,ngd
      !get the frequency
       wibx = eival(ibx)
       call mat_init(residuals(ibx),ndim,ndim)
       call mat_zero(residuals(ibx))
       call mat_init(rhs(ibx),ndim,ndim)
       
       rewind(lu_rho_rsp)
       rewind(lu_sigma_rsp)
       call expand_on_basis_minus(molcfg,ndim_red,red_X_glob(1:ndim_red*2,ibx),lu_rho_rsp,residuals(ibx))
       call mat_assign(S2XR(ibx),residuals(ibx))
       call mat_scal(-wibx,residuals(ibx))
       !add  sum_j ^RX_j(ibx) sigma_j
       call expand_on_basis(molcfg,ndim_red,red_X_glob(1:ndim_red*2,ibx),lu_sigma_rsp,residuals(ibx))
       
       call mat_init(residuals_img(ibx),ndim,ndim)
       call mat_zero(residuals_img(ibx))
       rewind(lu_rho_rsp)
       rewind(lu_sigma_rsp)
       call expand_on_basis_minus(molcfg,ndim_red,red_Ximg_glob(1:ndim_red*2,ibx),lu_rho_rsp,residuals_img(ibx))
       call mat_assign(S2XI(ibx),residuals_img(ibx))
       call mat_scal(-wibx,residuals_img(ibx))
       !add  sum_j ^RX_j(ibx) sigma_j
       call expand_on_basis(molcfg,ndim_red,red_Ximg_glob(1:ndim_red*2,ibx),lu_sigma_rsp,residuals_img(ibx))
     
      ! subtract gradient
      call mat_assign(rhs(ibx),gd(ibx))
      call mat_daxpy(-gammma,S2XI(ibx),rhs(ibx))
      call mat_daxpy(-1.0E0_realk,rhs(ibx),residuals(ibx))
      
      if (gd_complex) then
          call mat_init(rhsi(ibx),ndim,ndim)
          call mat_assign(rhsi(ibx),gdi(ibx))
          call mat_daxpy(-1.0E0_realk,rhsi(ibx),residuals_img(ibx))
      endif
      call mat_daxpy(-gammma,S2XR(ibx),residuals_img(ibx))
      
      call util_scriptPx('N',D,S,residuals(ibx))
      call util_scriptPx('N',D,S,residuals_img(ibx))
   enddo

! the residual(s) is now done: test for convergence. If not converged
! form new trial(s)  in next subroutine
! New trial(s) is equal to preconditioned residual(s)
   if (molcfg%solver%rsp_convdyn) then !default
      if (ITMIC == 1) then
         molcfg%solver%rsp_dyn_thresh=0E0_realk
         av_norm = 0.0E0_realk
         do i = 1, ngd
            res_norm = sqrt(mat_sqnorm2(residuals(i)))
            res_norm_img = sqrt(mat_sqnorm2(residuals_img(i)))
            av_norm = av_norm + res_norm + res_norm_img
         enddo
         av_norm = 0.5E0_realk*av_norm
         rsp_thresh = av_norm*molcfg%solver%rsp_conv_factor
         molcfg%solver%rsp_dyn_thresh=max(rsp_thresh,molcfg%solver%rsp_thresh)
         rsp_thresh=molcfg%solver%rsp_dyn_thresh
         WRITE(molcfg%lupri,*) 'Dynamic response convergence threshold set to', rsp_thresh
        else
         rsp_thresh=molcfg%solver%rsp_dyn_thresh
        endif
   else
      rsp_thresh=molcfg%solver%rsp_thresh
   endif

   do ibx = 1,ngd
      conv_root = .false.
      conv_root_img=.false.
      if (molcfg%solver%rsp_single_norm) then
          res_norm=mat_sqnorm2(residuals(ibx))
         res_norm_img=mat_sqnorm2(residuals_img(ibx))
      else ! double norm
         res_norm= sqrt(mat_sqnorm2(residuals(ibx)))
         res_norm_img=sqrt(mat_sqnorm2(residuals_img(ibx)))
      endif
      if (res_norm < rsp_thresh) conv_root = .true.
      if (res_norm_img < rsp_thresh) conv_root_img = .true.
      res_norm_tot=res_norm+res_norm_img

     ! Insert back print statements
      write (molcfg%lupri, '("Residual norm for real vector ", i3, " is: ", E14.6, "    and frequency = " , F12.6, "   It = ", i3, &
                        & " CONV =", l2)') ibx, res_norm, EIVAL(ibx), itmic, conv_root
    
      write (molcfg%lupri, '("Residual norm for img vector ", i3, " is: ", E14.6, "    and frequency = " , F12.6, "   It = ", i3, &
                        & " CONV =", l2)') ibx, res_norm_img, EIVAL(ibx), itmic, conv_root_img
         
         conv_test = res_norm
         conv_test_img = res_norm_img
      if (molcfg%solver%info_rsp) then
        write(molcfg%lupri,*) 'conv_test',conv_test
        write(molcfg%lupri,*) 'res_norm',res_norm
        write(molcfg%lupri,*) 'rsp_conv_thr',rsp_thresh
        write(molcfg%lupri,*) 'conv_test_img',conv_test_img
        write(molcfg%lupri,*) 'res_norm_img',res_norm_img
      endif
      if (conv_root .and. conv_root_img) then 
          conv=.true.
      endif
   enddo  
     
   do j=1,2
      call mat_zero(xsolvector(j))  
   enddo 
   nx_new=0
   if (.not. conv) then     
       do j=1,ngd
          if (.not. conv_root_img) then 
             nx_new=nx_new+1
             call mat_assign(xsolvector(nx_new),residuals_img(j))
          endif    
          if (.not. conv_root) then 
             nx_new=nx_new+1
             call mat_assign(xsolvector(nx_new),residuals(j))
          endif
       enddo
   endif

! Finished, deallocate local arrays
   do ibx=1,ngd
        call mat_free(residuals(ibx))
        call mat_free(residuals_img(ibx))
        call mat_free(S2XR(ibx))
        call mat_free(S2XI(ibx))
        call mat_free(rhs(ibx))
        if (gd_complex) then
            call mat_free(rhsi(ibx))
        endif
   enddo

  end subroutine get_complex_res
!================================================================  

  subroutine get_start_vectors(molcfg,F,D,S,RHS_real,RHS_img,ngd,eival,nx_new,xsolvector)
    implicit none
!================================================================ 
! A subroutine to get start vectors for damped response eqs.
! by default:
! X0[R]=prec(G[R])
! X0[I]=prec(G[I]+gammaS[2]X0[R]
!
!if a keyword .TWOSTART used in input, start vectors obtained by solving two sets of standard response eqs 
!(E[2]-wS[2])X0[R]=G[R] 
!(E[2]-wS[2])X0[I]=G[I]+gammaS[2]X0[R]
!================================================================ 
!    type(decompItem),intent(inout):: decomp
    type(rsp_molcfg), intent(inout) :: molcfg
    type(Matrix),intent(in)    :: F,D,S
    type(Matrix),intent(inout) :: RHS_real(:),RHS_img(:)
    type(Matrix),intent(inout) :: xsolvector(:)
    real(realk), intent(inout) :: eival(:) 
    integer,intent(inout)      :: nx_new
    integer,intent(in)         :: ngd
    type(Matrix)               :: S2XR(rsp_number_of_rhs),xsol(rsp_number_of_rhs)
    integer                    :: j,ndim
    real(realk)                :: gammma

 ndim=D%nrow
 gammma=molcfg%solver%rsp_gamma

do j=1,ngd
   call mat_init(S2XR(j),ndim,ndim)
   call mat_init(xsol(j),ndim,ndim)
   call mat_zero(S2XR(j))
   call mat_zero(xsol(j))
enddo

 if (molcfg%solver%rsp_damp_2start) then
    IF(molcfg%solver%UseExcitationVecs)then
       call rsp_init(rsp_number_of_current_trial,rsp_number_of_rhs, &
            &rsp_number_of_sols,rsp_number_of_omegas,&
            &rsp_number_of_startvecs+molcfg%solver%rsp_eigenvecs)
    ELSE
       call rsp_init(rsp_number_of_current_trial,rsp_number_of_rhs, &
            &rsp_number_of_sols,rsp_number_of_omegas,rsp_number_of_startvecs)
    ENDIF
    call rsp_solver(molcfg,D,S,F,.true.,1,RHS_real,EIVAL,xsol)
    nx_new=nx_new+1
    call mat_assign(xsolvector(nx_new),XSOL(1))
    call get_rho(molcfg,D,S,1,XSOL,S2XR)
    do j=1,ngd
       call mat_daxpy(gammma,S2XR(j),RHS_img(j))
       call mat_zero(xsol(j))
    enddo
    call rsp_solver(molcfg,D,S,F,.true.,1,RHS_img,EIVAL,XSOL)
    nx_new=nx_new+1
    call mat_assign(xsolvector(nx_new),XSOL(1))
 else
    do j=1,ngd
       nx_new=nx_new+1
       IF(molcfg%solver%rsp_MO_precond) then
          call MO_precond_complex(RHS_real(j),D,S,xsolvector(nx_new),molcfg%decomp%nocc)          
       ELSEIF(molcfg%solver%rsp_no_precond) then
          call mat_assign(xsolvector(nx_new),RHS_real(j))
       ELSE
          call lsquit('complex AO precond not implemented ',-1)
       ENDIF
    enddo
    call get_rho(molcfg,D,S,1,xsolvector,S2XR)
    do j=1,ngd
       call mat_daxpy(gammma,S2XR(j),RHS_img(j))
       nx_new=nx_new+1
       IF(molcfg%solver%rsp_MO_precond) then
          call MO_precond_complex(RHS_img(j),D,S,xsolvector(nx_new),molcfg%decomp%nocc)          
       ELSEIF(molcfg%solver%rsp_no_precond) then
          call mat_assign(xsolvector(nx_new),RHS_img(j))
       ELSE
          call lsquit('complex AO precond not implemented ',-1)
       ENDIF
    enddo
 endif

do j=1,ngd
   call mat_free(S2XR(j))
   call mat_free(xsol(j))
enddo
end subroutine get_start_vectors

 end module COMPLEXSOLVER
