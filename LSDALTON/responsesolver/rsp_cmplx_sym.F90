module COMPLEXSYMSOLVER
  use memory_handling  
  use decompMod
  use precision
  use matrix_module
  use matrix_operations
  use RSPSOLVER
  use RSPSYMSOLVER
  use COMPLEXSOLVER
  use COMPLEXNEWSYMSOLVER
  use rsp_util
  use files
  !use configuration
  private
  public ::  rsp_sym_complex_solver,rsp_sym_complex_init
  integer, save :: rsp_number_of_current_trial, rsp_number_of_rhs, rsp_number_of_sols, &
                 & rsp_number_of_omegas, rsp_number_of_startvecs,rsp_bvec_dim, &
                 & lu_x_rsp, lu_sigma_rsp, lu_rho_rsp,&
                 & lu_xi_rsp, lu_sigmai_rsp, lu_rhoi_rsp
                 !to complex solver
  real(realk),allocatable,save :: E1(:,:), E2(:,:), S1(:,:), S2(:,:), G1(:,:), G2(:,:)
  real(realk),allocatable,save :: E3(:,:), E4(:,:), S7(:,:), S8(:,:), G3(:,:), G4(:,:)
  logical :: RSPonMaster
contains
  subroutine rsp_sym_complex_init(ntrial, nrhs, nsol, nomega, nstart)
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
  end subroutine rsp_sym_complex_init

 !======================================================================== 
  subroutine rsp_sym_complex_solver(molcfg,F,D,S,ngd,GD,EIVAL,ifreq,XSOL,XSOLimg,gd_complex,gdi)
 !=======================================================================
 ! The complex response solver
 ! solving a complex equation problem using the symmetrized vectors
 ! (E[2]-(w+igamma)S[2])(X^R_g+X^R_u+iX^I_g+iX^I_u)=GD[R]_g+GD[R]_u+iGD[I]_g+iGD[I]_u
 ! (a) if gamma=0 and GD[I]=0, then solving the standard response eqs
 ! (b) solvng damped response eq. where the vectors are split due to its symmetry
 !
 ! INPUT FILE
 ! .COMPLEX (damping parameter - gamma)
 ! gamma given in input as: config%response%solver%rsp_gamma
 !
 ! input:
 ! ngd    number of RHS
 ! GD     RHS
 ! EIVAL  frequency
 ! ifreq  no. of frequency
 ! gd_complex .true. if imaginary RHS present
 !            .false. if not
 !gdi imaginary RHS
 !
 ! output:
 ! XSOL    real solution vector
 ! XSOLimg imaginary solution vector
 !
 ! Joanna, July 2010
 !=========================================================================
    implicit none
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
    integer                     :: i,nfreqs

write(molcfg%lupri,*) 'Entering the complex sym rsp solver, gamma=', molcfg%solver%rsp_gamma
if ((abs(molcfg%solver%rsp_gamma) .LT. 1E-8_realk) .and. (.not. gd_complex)) then
    if (molcfg%solver%rsp_stdnew) then  
       call rsp_sym_init(rsp_number_of_current_trial,rsp_number_of_rhs,&
                        &rsp_number_of_sols,rsp_number_of_omegas,rsp_number_of_startvecs)
       call rsp_sym_solver(molcfg,D,S,F,.true.,1,GD,EIVAL,XSOL)
    else
       if (molcfg%solver%UseExcitationVecs) then
          call rsp_init(rsp_number_of_current_trial,rsp_number_of_rhs,&
                       & rsp_number_of_sols,rsp_number_of_omegas,&
                       & rsp_number_of_startvecs+molcfg%solver%rsp_eigenvecs)
       else
          call rsp_init(rsp_number_of_current_trial,rsp_number_of_rhs,&
                       & rsp_number_of_sols,rsp_number_of_omegas,rsp_number_of_startvecs)
       endif
       call rsp_solver(molcfg,D,S,F,.true.,1,GD,EIVAL,XSOL)
    endif
    do i=1,rsp_number_of_sols
       call mat_zero(xsolimg(i))
    enddo
else
   if (molcfg%solver%rsp_cmplxnew) then
      call symcomplex_solver(molcfg,F,D,S,1,GD,GDI,EIVAL,XSOL,XSOLimg,gd_complex)
 !  call complex_solver_check(F,D,S,gd(1),XSOLimg(1),XSOL(1),EIVAL(ifreq),&
  !                         & gd_complex,conv,gdi(1))
   elseif (molcfg%solver%rsp_cpp) then
      nfreqs=rsp_number_of_omegas
      call new_symcomplex_solver(molcfg,F,D,S,1,nfreqs,GD,GDI,EIVAL,XSOL,XSOLimg,gd_complex)
   else
      write(molcfg%solver%rsp_cpp,*) 'The code should terminate'
   endif
endif
    
    end subroutine rsp_sym_complex_solver   
!==========================================================
!==============================================================
    subroutine symcomplex_solver(molcfg,F,D,S,ngd,GDB,GDI,EIVAL,XSOL,XSOLimg,gd_complex)
      !==============================================================
      ! The linear scaling complex response solver
      ! solving a complex equation problem using the symmetrized trial vectors algorithm
      ! (E[2]-(w+igamma)S[2])(X^R_g+X^R_u+iX^I_g+iX^I_u)=GD[R]_g+GD[R]_u+iGD[I]_g+iGD[I]_u
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
      ! Joanna, May 2010
      !=========================================================================
      !
      use RSPSYMSOLVER
      implicit none
      type(rsp_molcfg), intent(inout)     :: molcfg
      type(Matrix),intent(in)             :: F,D,S
      type(Matrix),intent(inout)          :: GDB(rsp_number_of_rhs)
      type(Matrix),intent(inout),optional :: GDI(rsp_number_of_rhs)
      real(realk),intent(inout)           :: eival(rsp_number_of_omegas)
      type(Matrix),intent(inout)          :: XSOL(rsp_number_of_sols),XSOLimg(rsp_number_of_sols)
      logical, intent(in)                 :: gd_complex
      integer, intent(in)                 :: ngd
      type(Matrix)                        :: xp(rsp_number_of_rhs),xm(rsp_number_of_rhs),x
      type(Matrix)                        :: xpi(rsp_number_of_rhs),xmi(rsp_number_of_rhs),scrT
      real(realk)                         :: gammma,res_norm_tot,r1,r2,r3
      integer                             :: i,j,ndim,k,n_red,nx_new,m,l,n_v,nm_new,nm_red,n_i
      integer                             ::nxi_new,nmi_new,ni_red,nmi_red,ni_v,ni_i,nt_red,nti_red
      logical                             :: conv,make_rhos,conv_real,conv_img
      real(realk),allocatable             :: red_Xp_glob(:,:),red_Xm_glob(:,:)
      real(realk),allocatable             :: red_Xpi_glob(:,:),red_Xmi_glob(:,:)
      real(realk)                         :: absorp, disper
      type(Matrix)                        :: RHS_real(rsp_number_of_rhs),RHS_img(rsp_number_of_rhs)
      type(Matrix)                        :: Xsolvecp(rsp_number_of_rhs),sigmasp(rsp_number_of_rhs)
      type(Matrix)                        :: rhosm(rsp_number_of_rhs),rhosp(rsp_number_of_rhs)
      type(Matrix)                        :: Xsolvecm(rsp_number_of_rhs),sigmasm(rsp_number_of_rhs)
      type(Matrix)                        :: Xsolvecpi(rsp_number_of_rhs),sigmaspi(rsp_number_of_rhs)
      type(Matrix)                        :: rhosmi(rsp_number_of_rhs),rhospi(rsp_number_of_rhs)
      type(Matrix)                        :: Xsolvecmi(rsp_number_of_rhs),sigmasmi(rsp_number_of_rhs)
      type(Matrix)                        :: Xsolvec(rsp_number_of_rhs)
      !type(Matrix)                :: sigmas(rsp_number_of_rhs),rhos(rsp_number_of_rhs)
      type(Matrix),pointer                :: sigmas(:),rhos(:)
      type(Matrix)                        :: Xsolveci(rsp_number_of_rhs)
      !type(Matrix)                        :: sigmasi(rsp_number_of_rhs),rhosi(rsp_number_of_rhs)
      type(Matrix),pointer                :: sigmasi(:),rhosi(:)
      integer                             :: max_red, max_it, rsp_thresh
      ! call mem_alloc(rhos,rsp_number_of_rhs)
      ! call mem_alloc(sigmas,rsp_number_of_rhs)
      ! call mem_alloc(rhosi,rsp_number_of_rhs)
      ! call mem_alloc(sigmasi,rsp_number_of_rhs)
      max_red= molcfg%solver%rsp_maxred
      max_it= molcfg%solver%rsp_maxit
      rsp_thresh= molcfg%solver%rsp_thresh
      gammma=molcfg%solver%rsp_gamma

      ndim = S%nrow
      lu_x_rsp= -1 ; lu_sigma_rsp = -1 ; lu_rho_rsp = -1
      lu_xi_rsp= -1 ; lu_sigmai_rsp = -1 ; lu_rhoi_rsp = -1
      CALL LSOPEN(lu_x_rsp,'xrsp','unknown','UNFORMATTED')
      CALL LSOPEN(lu_sigma_rsp,'sigma_rsp','unknown','UNFORMATTED')
      CALL LSOPEN(lu_rho_rsp,'rho_rsp','unknown','UNFORMATTED')
      CALL LSOPEN(lu_xi_rsp,'xirsp','unknown','UNFORMATTED')
      CALL LSOPEN(lu_sigmai_rsp,'sigmai_rsp','unknown','UNFORMATTED')
      CALL LSOPEN(lu_rhoi_rsp,'rhoi_rsp','unknown','UNFORMATTED')

      allocate(red_Xp_glob(max_red,ngd),red_Xm_glob(max_red,ngd))
      allocate(red_Xpi_glob(max_red,ngd),red_Xmi_glob(max_red,ngd))
      ALLOCATE(E1(max_red,max_red),E2(max_red,max_red))
      ALLOCATE(E3(max_red,max_red),E4(max_red,max_red))
      allocate(S1(max_red,max_red),S2(max_red,max_red))
      allocate(S7(max_red,max_red),S8(max_red,max_red))
      allocate(G1(max_red,ngd))  
      allocate(G2(max_red,ngd))  
      allocate(G3(max_red,ngd))  
      allocate(G4(max_red,ngd))  

      red_Xp_glob = 0.0E0_realk ; red_Xm_glob = 0.0E0_realk
      red_Xpi_glob = 0.0E0_realk ; red_Xmi_glob = 0.0E0_realk
      E1=0E0_realk; E2=0E0_realk; S1=0E0_realk; S2=0E0_realk; G1=0E0_realk; G2=0E0_realk
      E3=0E0_realk; E4=0E0_realk; S7=0E0_realk; S8=0E0_realk; G3=0E0_realk; G4=0E0_realk

      n_red=0
      nm_red=0
      ni_red=0
      nmi_red=0
      nt_red=0
      nti_red=0
      conv_real=.false.
      conv_img=.false.
      conv=.false.
      nx_new=0
      nm_new=0
      nxi_new=0
      nmi_new=0

      call mat_zero(xsolimg(1))   
      do j=1,ngd
         call mat_init(xp(j),ndim,ndim)
         call mat_init(xm(j),ndim,ndim)
         call mat_init(xpi(j),ndim,ndim)
         call mat_init(xmi(j),ndim,ndim)
         call mat_init(RHS_real(j),ndim,ndim)
         call mat_init(RHS_img(j),ndim,ndim)
         call mat_zero(xp(j))
         call mat_zero(xm(j))
         call mat_zero(xpi(j))
         call mat_zero(xmi(j))
         call mat_zero(RHS_real(j))
         call mat_zero(RHS_img(j))
         call mat_assign(RHS_real(j),GDb(j))
         if (gd_complex) then
            call mat_assign(RHS_img(j),GDi(j))
         endif
      enddo

      !get start vectors by preconditioning RHS with efficient preconditoner by Villaume et al.
      call get_start_vectors(molcfg,F,D,S,RHS_real,RHS_img,ngd,eival,nx_new,xp,xm,xpi,xmi)


      do i = 1,max_it
         write(molcfg%lupri,*) '------------------'
         write(molcfg%lupri,*) 'Start macroiteration:',i
         write(molcfg%lupri,*) '------------------'
         nm_new  = nx_new
         nxi_new = nx_new
         nmi_new = nx_new
         n_v  = nx_new 
         n_i  = nm_new
         ni_v = nxi_new 
         ni_i = nmi_new
         do j=1,nx_new
            call mat_init(xsolvecp(j),ndim,ndim)
            call mat_init(sigmasp(j),ndim,ndim)
            call mat_init(rhosm(j),ndim,ndim)
            call mat_init(xsolvecpi(j),ndim,ndim)
            call mat_init(sigmaspi(j),ndim,ndim)
            call mat_init(rhosmi(j),ndim,ndim)
            call mat_zero(xsolvecp(j))
            call mat_zero(xsolvecpi(j))
         enddo
         do j=1,nm_new
            call mat_init(xsolvecm(j),ndim,ndim)
            call mat_init(sigmasm(j),ndim,ndim)
            call mat_init(rhosp(j),ndim,ndim)
            call mat_init(xsolvecmi(j),ndim,ndim)
            call mat_init(sigmasmi(j),ndim,ndim)
            call mat_init(rhospi(j),ndim,ndim)
            call mat_zero(xsolvecm(j))
            call mat_zero(xsolvecmi(j))
         enddo
         do j=1,max(nx_new,nm_new)
            call mat_init(xsolvec(j),ndim,ndim)
            !        call mat_init(sigmas(j),ndim,ndim)
            !        call mat_init(rhos(j),ndim,ndim)
            call mat_init(xsolveci(j),ndim,ndim)
            !        call mat_init(sigmasi(j),ndim,ndim)
            !        call mat_init(rhosi(j),ndim,ndim)
            call mat_zero(xsolvec(j))
            call mat_zero(xsolveci(j))
         enddo
         do j=1,nx_new
            call mat_add(1E0_realk,xp(j),1E0_realk,xm(j),xsolvecp(j))
            call mat_add(1E0_realk,xpi(j),1E0_realk,xmi(j),xsolvecpi(j))
         enddo

         call orthonormalize22(molcfg,D,S,xsolvecp,nx_new,nt_red,xsolvec,lu_x_rsp)
         call orthonormalize22(molcfg,D,S,xsolvecpi,nxi_new,nt_red,xsolveci,lu_xi_rsp)



         !  call orthonormalize21(molcfg,D,S,xp,xm,nx_new,nm_new,nt_red,xsolvecp,xsolvecm,lu_x_rsp,.true.)
         !  call orthonormalize21(molcfg,D,S,xpi,xmi,nxi_new,nmi_new,nti_red,xsolvecpi,xsolvecmi,lu_xi_rsp,.true.)

         !the number of symmetric and antisymmetric components might be different
         !   do j=1,min(nx_new,nm_new)
         !   call mat_add(1E0_realk,xsolvecp(j),1E0_realk,xsolvecm(j),xsolvec(j))
         !   enddo
         !   if (nx_new>nm_new) then
         !      do j=nm_new+1,nx_new
         !         xsolvec(j)=xsolvecp(j)
         !      enddo
         !   elseif (nm_new>nx_new) then
         !      do j=nx_new+1,nm_new
         !         xsolvec(j)=xsolvecm(j)
         !      enddo
         !  endif
         !   do j=1,min(nxi_new,nmi_new)
         !   call mat_add(1E0_realk,xsolvecpi(j),1E0_realk,xsolvecmi(j),xsolveci(j))
         !   enddo
         !   if (nxi_new>nmi_new) then
         !      do j=nmi_new+1,nxi_new
         !         xsolveci(j)=xsolvecpi(j)
         !      enddo
         !   elseif (nmi_new>nxi_new) then
         !      do j=nxi_new+1,nmi_new
         !         xsolveci(j)=xsolvecmi(j)
         !      enddo
         !  endif

         call transform_vectors(molcfg,D,S,F,max(nx_new,nm_new),xsolvec,sigmas,rhos,.true.)
         call transform_vectors(molcfg,D,S,F,max(nxi_new,nmi_new),xsolveci,sigmasi,rhosi,.true.)

         do j=1,nx_new
            call mat_init(scrT,ndim,ndim)
            call mat_trans(sigmas(j),scrT)  

            call mat_add(0.5E0_realk,sigmas(j),0.5E0_realk,scrT,sigmasp(j))
            call mat_add(0.5E0_realk,sigmas(j),-0.5E0_realk,scrT,sigmasm(j))

            call mat_trans(rhos(j),scrT)  

            call mat_add(0.5E0_realk,rhos(j),0.5E0_realk,scrT,rhosp(j))
            call mat_add(0.5E0_realk,rhos(j),-0.5E0_realk,scrT,rhosm(j))

            call mat_trans(xsolvec(j),scrT)

            call mat_add(0.5E0_realk,xsolvec(j),0.5E0_realk,scrT,xsolvecp(j))
            call mat_add(0.5E0_realk,xsolvec(j),-0.5E0_realk,scrT,xsolvecm(j))
            call mat_free(scrT)
         enddo
         !    if (nx_new>nm_new) then
         !      do j=nm_new+1,nx_new
         !         sigmasp(j)=sigmas(j)
         !         rhosm(j)=rhos(j)
         !      enddo
         !    elseif (nm_new>nx_new) then
         !      do j=nx_new+1,nm_new
         !         sigmasm(j)=sigmas(j)
         !         rhosp(j)=rhos(j)
         !      enddo
         !   endif 
         do j=1,min(nxi_new,nmi_new)
            call mat_init(scrT,ndim,ndim)
            call mat_trans(sigmasi(j),scrT)  

            call mat_add(0.5E0_realk,sigmasi(j),0.5E0_realk,scrT,sigmaspi(j))
            call mat_add(0.5E0_realk,sigmasi(j),-0.5E0_realk,scrT,sigmasmi(j))

            call mat_trans(rhosi(j),scrT)  

            call mat_add(0.5E0_realk,rhosi(j),0.5E0_realk,scrT,rhospi(j))
            call mat_add(0.5E0_realk,rhosi(j),-0.5E0_realk,scrT,rhosmi(j))

            call mat_trans(xsolveci(j),scrT)  

            call mat_add(0.5E0_realk,xsolveci(j),0.5E0_realk,scrT,xsolvecpi(j))
            call mat_add(0.5E0_realk,xsolveci(j),-0.5E0_realk,scrT,xsolvecmi(j))


            call mat_free(scrT)
         enddo
         if (nxi_new>nmi_new) then
            do j=nmi_new+1,nxi_new
               call mat_assign(sigmaspi(j),sigmasi(j))
               call mat_assign(rhosmi(j),rhosi(j))
            enddo
         elseif (nmi_new>nxi_new) then
            do j=nxi_new+1,nmi_new
               call mat_assign(sigmasmi(j),sigmasi(j))
               call mat_assign(rhospi(j),rhosi(j))
            enddo
         endif

         if (gd_complex) then
            !        !build a reduce space
            call extend_new_complex_reduced_matrices(molcfg,n_red,nm_red,ni_red,nmi_red,nt_red,nti_red,nx_new,nm_new,nxi_new,&
                 &nmi_new,ngd,RHS_real,xsolvecp,sigmasp,rhosm,xsolvecpi,sigmaspi,rhosmi,xsolvecm,sigmasm,rhosp,&
                 &xsolvecmi,sigmasmi,rhospi,.true.,RHS_img)
            !write on disk. The sum of symm and antisymm component is stored together
            do j = 1,max(nx_new,nm_new)
               call mat_write_to_disk(lu_x_rsp,xsolvec(j),RSPonMaster)
               call mat_write_to_disk(lu_sigma_rsp,sigmas(j),RSPonMaster)
               call mat_write_to_disk(lu_rho_rsp,rhos(j),RSPonMaster)
            enddo

            n_red = n_red + nx_new !update reduced space dimension
            nm_red = nm_red + nm_new !update reduced space dimension
            nt_red=nt_red+max(nx_new,nm_new)

            do j = 1,max(nxi_new,nmi_new)
               call mat_write_to_disk(lu_xi_rsp,xsolveci(j),RSPonMaster)
               call mat_write_to_disk(lu_sigmai_rsp,sigmasi(j),RSPonMaster)
               call mat_write_to_disk(lu_rhoi_rsp,rhosi(j),RSPonMaster)
            enddo

            ni_red = ni_red + nxi_new !update reduced space dimension
            nmi_red = nmi_red + nmi_new !update reduced space dimension
            nti_red=nti_red+max(nxi_new,nmi_new)
            !        
            !          !solve complex reduced equation
            call solve_complex(molcfg,ndim,n_red,nm_red,ni_red,nmi_red,ngd,eival,&
                 &red_Xp_glob,red_Xpi_glob,red_Xm_glob,red_Xmi_glob,.true.)
            !    
            !        ! check for convergece and get new trialvectors
            do j=1,ngd
               call mat_zero(xp(j))
               call mat_zero(xm(j))
               call mat_zero(xpi(j))
               call mat_zero(xmi(j))
            enddo
            call get_complex_res(molcfg,F,D,S,i,n_red,nm_red,ni_red,nmi_red,nt_red,&
                 &nti_red,ngd,RHS_real,red_Xp_glob,red_Xpi_glob,red_Xm_glob,&
                 &red_Xmi_glob,eival,xp,xm,xpi,xmi,nx_new,conv,res_norm_tot,&
                 &.true.,RHS_img)
         else
            !build a reduce space
            call extend_new_complex_reduced_matrices(molcfg,n_red,nm_red,ni_red,&
                 &nmi_red,nt_red,nti_red,nx_new,nm_new,nxi_new,nmi_new,ngd,&
                 &RHS_real,xsolvecp,sigmasp,rhosm,xsolvecpi,sigmaspi,rhosmi,&
                 &xsolvecm,sigmasm,rhosp,xsolvecmi,sigmasmi,rhospi,.false.)

            do j = 1,max(nx_new,nm_new)
               call mat_write_to_disk(lu_x_rsp,xsolvec(j),RSPonMaster)
               call mat_write_to_disk(lu_sigma_rsp,sigmas(j),RSPonMaster)
               call mat_write_to_disk(lu_rho_rsp,rhos(j),RSPonMaster)
            enddo

            n_red = n_red + nx_new !update reduced space dimension
            nm_red = nm_red + nm_new !update reduced space dimension
            nt_red=nt_red+max(nx_new,nm_new)

            do j = 1,max(nxi_new,nmi_new)
               call mat_write_to_disk(lu_xi_rsp,xsolveci(j),RSPonMaster)
               call mat_write_to_disk(lu_sigmai_rsp,sigmasi(j),RSPonMaster)
               call mat_write_to_disk(lu_rhoi_rsp,rhosi(j),RSPonMaster)
            enddo

            ni_red = ni_red + nxi_new !update reduced space dimension
            nmi_red = nmi_red + nmi_new !update reduced space dimension
            nti_red=nti_red+max(nxi_new,nmi_new)

            !solve complex reduced equation
            call solve_complex(molcfg,ndim,n_red,nm_red,ni_red,nmi_red,ngd,&
                 &eival,red_Xp_glob,red_Xpi_glob,red_Xm_glob,red_Xmi_glob,.false.)

            ! check for convergece and get new trialvectors
            do j=1,ngd
               call mat_zero(xp(j))
               call mat_zero(xm(j))
               call mat_zero(xpi(j))
               call mat_zero(xmi(j))
            enddo
            call get_complex_res(molcfg,F,D,S,i,n_red,nm_red,ni_red,nmi_red,nt_red,&
                 &nti_red,ngd,RHS_real,red_Xp_glob,red_Xpi_glob,red_Xm_glob,&
                 &red_Xmi_glob,eival,xp,xm,xpi,xmi,nx_new,conv,res_norm_tot,.false.)
         endif

         if (conv)  then
            do l=1,ngd
               rewind(lu_x_rsp)
               call expand_on_basis4(molcfg,nt_red,ndim,RED_Xp_glob(1:n_red,l),RED_Xm_glob(1:nm_red,l),lu_x_rsp,Xp(l),Xm(l))
               rewind(lu_xi_rsp)
               call expand_on_basis4(molcfg,nti_red,ndim,RED_Xpi_glob(1:ni_red,l),RED_Xmi_glob(1:nmi_red,l),lu_xi_rsp,Xpi(l),Xmi(l))
               call mat_zero(xsol(l))
               call mat_zero(xsolimg(l))
               call mat_add(1E0_realk,xp(l),1E0_realk,xm(l),xsol(l))
               call mat_add(1E0_realk,xpi(l),1E0_realk,xmi(l),xsolimg(l))

            enddo
            do j=1,n_v
               call mat_free(xsolvecp(j))
               call mat_free(sigmasp(j))
               call mat_free(rhosm(j))
            enddo
            do j=1,ni_v
               call mat_free(xsolvecpi(j))
               call mat_free(sigmaspi(j))
               call mat_free(rhosmi(j))
            enddo
            do j=1,n_i
               call mat_free(xsolvecm(j))
               call mat_free(sigmasm(j))
               call mat_free(rhosp(j))
            enddo
            do j=1,ni_i
               call mat_free(xsolvecmi(j))
               call mat_free(sigmasmi(j))
               call mat_free(rhospi(j))
            enddo
            do j=1,ngd
               call mat_free(xp(j))
               call mat_free(xm(j))
               call mat_free(xpi(j))
               call mat_free(xmi(j))
            enddo
            do j=1,max(n_i,n_v)
               call mat_free(xsolvec(j))
               call mat_free(sigmas(j))
               call mat_free(rhos(j))
            enddo
            call mem_dealloc(rhos)
            call mem_dealloc(sigmas)
            do j=1,max(ni_i,ni_v)
               call mat_free(xsolveci(j))
               call mat_free(sigmasi(j))
               call mat_free(rhosi(j))
            enddo
            call mem_dealloc(rhosi)
            call mem_dealloc(sigmasi)
            exit
         else         
            do l=1,ngd
               rewind(lu_x_rsp)
               call expand_on_basis4(molcfg,nt_red,ndim,RED_Xp_glob(1:n_red,l),RED_Xm_glob(1:nm_red,l),lu_x_rsp,Xp(l),Xm(l))
               rewind(lu_xi_rsp)
               call expand_on_basis4(molcfg,nti_red,ndim,RED_Xpi_glob(1:ni_red,l),RED_Xmi_glob(1:nmi_red,l),lu_xi_rsp,Xpi(l),Xmi(l))
               call mat_zero(xsol(l))
               call mat_zero(xsolimg(l))
               call mat_add(1E0_realk,xp(l),1E0_realk,xm(l),xsol(l))
               call mat_add(1E0_realk,xpi(l),1E0_realk,xmi(l),xsolimg(l))

               absorp=mat_dotproduct(gdb(l),xsol(l))
               write(molcfg%lupri,*) "Dispersion value macroiteration:",i, "for &
                    & frequency",eival(l),"=", absorp

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
                     call mat_free(xsolvecp(j))
                     call mat_free(sigmasp(j))
                     call mat_free(rhosm(j))
                  enddo
                  do j=1,ni_v
                     call mat_free(xsolvecpi(j))
                     call mat_free(sigmaspi(j))
                     call mat_free(rhosmi(j))
                  enddo
                  do j=1,n_i
                     call mat_free(xsolvecm(j))
                     call mat_free(sigmasm(j))
                     call mat_free(rhosp(j))
                  enddo
                  do j=1,ni_i
                     call mat_free(xsolvecmi(j))
                     call mat_free(sigmasmi(j))
                     call mat_free(rhospi(j))
                  enddo
                  do j=1,max(n_i,n_v)
                     call mat_free(xsolvec(j))
                     call mat_free(sigmas(j))
                     call mat_free(rhos(j))
                  enddo
                  call mem_dealloc(rhos)
                  call mem_dealloc(sigmas)
                  do j=1,max(ni_i,ni_v)
                     call mat_free(xsolveci(j))
                     call mat_free(sigmasi(j))
                     call mat_free(rhosi(j))
                  enddo
                  call mem_dealloc(rhosi)
                  call mem_dealloc(sigmasi)
                  do j=1,ngd
                     call mat_free(xp(j))
                     call mat_free(xm(j))
                     call mat_free(xpi(j))
                     call mat_free(xmi(j))
                  enddo
                  exit
               endif
            else
               r1=r2
               r2=r3
               r3=res_norm_tot
               if ((abs(r1-r2)<rsp_thresh) .and. (abs(r2-r3)<rsp_thresh)) then
                  do j=1,n_v
                     call mat_free(xsolvecp(j))
                     call mat_free(sigmasp(j))
                     call mat_free(rhosm(j))
                  enddo
                  do j=1,ni_v
                     call mat_free(xsolvecpi(j))
                     call mat_free(sigmaspi(j))
                     call mat_free(rhosmi(j))
                  enddo
                  do j=1,n_i
                     call mat_free(xsolvecm(j))
                     call mat_free(sigmasm(j))
                     call mat_free(rhosp(j))
                  enddo
                  do j=1,ni_i
                     call mat_free(xsolvecmi(j))
                     call mat_free(sigmasmi(j))
                     call mat_free(rhospi(j))
                  enddo
                  do j=1,max(n_i,n_v)
                     call mat_free(xsolvec(j))
                     call mat_free(sigmas(j))
                     call mat_free(rhos(j))
                  enddo
                  call mem_dealloc(rhos)
                  call mem_dealloc(sigmas)
                  do j=1,max(ni_i,ni_v)
                     call mat_free(xsolveci(j))
                     call mat_free(sigmasi(j))
                     call mat_free(rhosi(j))
                  enddo
                  call mem_dealloc(rhosi)
                  call mem_dealloc(sigmasi)
                  do j=1,ngd
                     call mat_free(xp(j))
                     call mat_free(xm(j))
                     call mat_free(xpi(j))
                     call mat_free(xmi(j))
                  enddo
                  exit
               endif
            endif
         endif
         do j=1,n_v
            call mat_free(xsolvecp(j))
            call mat_free(sigmasp(j))
            call mat_free(rhosm(j))
         enddo
         do j=1,ni_v
            call mat_free(xsolvecpi(j))
            call mat_free(sigmaspi(j))
            call mat_free(rhosmi(j))
         enddo
         do j=1,n_i
            call mat_free(xsolvecm(j))
            call mat_free(sigmasm(j))
            call mat_free(rhosp(j))
         enddo
         do j=1,ni_i
            call mat_free(xsolvecmi(j))
            call mat_free(sigmasmi(j))
            call mat_free(rhospi(j))
         enddo
         do j=1,max(n_i,n_v)
            call mat_free(xsolvec(j))
            call mat_free(sigmas(j))
            call mat_free(rhos(j))
         enddo
         call mem_dealloc(rhos)
         call mem_dealloc(sigmas)
         do j=1,max(ni_i,ni_v)
            call mat_free(xsolveci(j))
            call mat_free(sigmasi(j))
            call mat_free(rhosi(j))
         enddo
         call mem_dealloc(rhosi)
         call mem_dealloc(sigmasi)
      enddo
      do i=1,ngd
         call mat_free(RHS_real(i))
         call mat_free(RHS_img(i))
      enddo
      deallocate(E1,E2,E3,E4,S1,S2,G1)
      deallocate(S7,S8,G2,G3,G4)
      deallocate(red_Xp_glob,red_Xm_glob)
      deallocate(red_Xpi_glob,red_Xmi_glob)

      CALL lsCLOSE(lu_sigma_rsp,'DELETE')
      CALL lsCLOSE(lu_rho_rsp,'DELETE')
      CALL lsCLOSE(lu_x_rsp,'DELETE')
      CALL lsCLOSE(lu_sigmai_rsp,'DELETE')
      CALL lsCLOSE(lu_rhoi_rsp,'DELETE')
      CALL lsCLOSE(lu_xi_rsp,'DELETE')

      !================================================================  
    end subroutine symcomplex_solver
!================================================================  
!================================================================  
subroutine extend_new_complex_reduced_matrices(molcfg,n_red,nm_red,ni_red,&
     &nmi_red,nt_red,nti_red,nx_new,nm_new,nxi_new,nmi_new,ngd,gd,&
     &xsolvecp,sigmasp,rhosm, xsolvecpi,sigmaspi,rhosmi,xsolvecm,&
     &sigmasm,rhosp,xsolvecmi,sigmasmi,rhospi,gd_complex,gdi)
!================================================================  
! Purpose: extend reduced space with new vectors
!          A reduced space for each component is build separtely
!
!           E1     -wS1    gammaS2  0          Xg^R      Gg^R
!          -wS3     E2        0   gammaS4      Xu^R      Gu^R
!       (                                  ) (      )= (    )
!           gammaS5 0       -E3     wS6        Xu^I      Gu^I
!               0  gammaS7   wS8    -E4        Xg^I      Gg^I
!       where
!          E1 = Xg^R E[2] Xg^R
!          E2 = Xu^R E[2] Xu^R
!          E3 = Xu^I E[2] Xu^I
!          E4 = Xg^I E[2] Xg^I
!
!          the reduced spaces E1, E2, .. are global variables  
!          (see top of file) 
! INPUT: 
! n_red, number of g^R vectors on disk and current half size of reduced space
! nm_red, number of u^R vectors on disk and current half size of reduced space
! n_red, number of g^I vectors on disk and current half size of reduced space
! nm_red, number of u^I vectors on disk and current half size of reduced space
! nt_red, number of r vectors stored on disk (vectors are stored as g+u, 
!         in case when the trial vector is either symmetric or antisymmetric
!          the total number of trial vectors may be different than number of g/u)
! nti_red, number of i vectors stored on disk
! nx_NEW, number of new g^R vectors
! nm_NEW, number of new u^R vectors
! nx_NEW, number of new g^I vectors
! nm_NEW, number of new u^I vectors
! ngd, the number of response vectors to be found
! gd(ngd), matrix array containing ngd real right hand sides (gradients)
! gdi(ngd, matrix array containing ngd img right hand sides (gradients)
! sigmasp,rhosp,xsolvecp:  matrix arrays of the new g^R vectors
! sigmasm,rhosm,xsolvecm:  matrix arrays of the new u^R vectors
! sigmaspi,rhospi,xsolvecpi:  matrix arrays of the new g^I vectors
! sigmasmi,rhosmi,xsolvecmi:  matrix arrays of the new u^I vectors
!
!Joanna, July 2010
!================================================================  
implicit none
type(rsp_molcfg), intent(inout) :: molcfg
integer,intent(in)            :: n_red,nm_red,nt_red
integer,intent(in)            :: ni_red,nmi_red,nti_red
integer, intent(in)              :: nx_new,nm_new,ngd
integer, intent(in)              :: nxi_new,nmi_new
type(Matrix),intent(in)          :: xsolvecp(:),sigmasp(:),rhosm(:)
type(Matrix),intent(in)          :: xsolvecpi(:),sigmaspi(:),rhosmi(:)
type(Matrix),intent(in)          :: xsolvecm(:),sigmasm(:),rhosp(:)
type(Matrix),intent(in)          :: xsolvecmi(:),sigmasmi(:),rhospi(:)
type(Matrix),intent(in)          :: gd(:)
type(Matrix),intent(in),optional :: gdi(:)
type(Matrix)                     :: gdp,gdpi,gdm,gdmi
logical,intent(in)               :: gd_complex
integer                          :: ndim,k,j,i,l,max_red
type(Matrix)                     ::xi_j,sigmai_j,rhoi_j,x_j,rho_j,sigma_j,gT
type(Matrix)                     ::xpi_j,sigmapi_j,xp_j,sigmap_j,scrT
type(Matrix)                     ::xmi_j,sigmami_j,rhopi_j,xm_j,rhop_j,sigmam_j

 max_red= molcfg%solver%rsp_maxred
if ((N_red + Nx_new) > max_red) then
   WRITE (molcfg%LUPRI,*) 'ERROR in extend_complex_reduced_matrices: N_RED + N_NEW > rsp_maxred'
   WRITE (molcfg%LUPRI,*) 'N_RED ,Nx_NEW  =',N_red,Nx_NEW
   WRITE (molcfg%LUPRI,*) 'rsp_maxred =',max_red
   WRITE (molcfg%LUPRI,*) 'Reduce problem size or recompile orthonormalize '// &
        &                  'with larger rsp_maxred parameter'
   CALL lsQUIT('orthonormalize error: N_RED + Nx_NEW > rsp_maxred parameter',molcfg%lupri)
END IF
if ((nm_red+nm_new) > max_red) then
   WRITE (molcfg%LUPRI,*) 'ERROR in extend_complex_reduced_matrices: Nm_red + Nm_NEW > rsp_maxred'
   WRITE (molcfg%LUPRI,*) 'Nm_RED ,Nm_NEW  =',Nm_red,Nm_NEW
   WRITE (molcfg%LUPRI,*) 'rsp_maxred =',max_red
   WRITE (molcfg%LUPRI,*) 'Reduce problem size or recompile orthonormalize '// &
        &                  'with larger rsp_maxred parameter'
   CALL lsQUIT('orthonormalize error: Nm_red + Nm_NEW > rsp_maxred parameter',molcfg%lupri)
END IF

ndim = gd(1)%nrow
call mat_init(sigma_j,ndim,ndim)
call mat_init(rho_j,ndim,ndim)
call mat_init(x_j,ndim,ndim)
call mat_init(sigmai_j,ndim,ndim)
call mat_init(rhoi_j,ndim,ndim)
call mat_init(xi_j,ndim,ndim)
call mat_init(sigmap_j,ndim,ndim)
call mat_init(rhop_j,ndim,ndim)
call mat_init(sigmam_j,ndim,ndim)
call mat_init(xp_j,ndim,ndim)
call mat_init(sigmapi_j,ndim,ndim)
call mat_init(xpi_j,ndim,ndim)
call mat_init(xm_j,ndim,ndim)
call mat_init(sigmami_j,ndim,ndim)
call mat_init(rhopi_j,ndim,ndim)
call mat_init(xmi_j,ndim,ndim)
call mat_init(scrT,ndim,ndim)

if (nx_new>0) then
    k = 0
    do i = n_red+1, n_red+nx_new 
       k = k + 1
       l = 0
       do j = n_red+1, n_red+nx_new
          l = l + 1
          E1(i,j)   =  mat_dotproduct(xsolvecp(k),sigmasp(l))
       enddo      
    enddo
endif
if (nm_new>0) then
    k = 0
    do i = nm_red+1, nm_red+nm_new 
       k = k + 1
       l = 0
       do j = nm_red+1, nm_red+nm_new
          l = l + 1
          E2(i,j)   =  mat_dotproduct(xsolvecm(k),sigmasm(l))
       enddo      
    enddo
endif
if (nmi_new>0) then
    k = 0
    do i = nmi_red+1, nmi_red+nmi_new 
       k = k + 1
       l = 0
       do j = nmi_red+1, nmi_red+nmi_new
          l = l + 1
          E3(i,j)   =  mat_dotproduct(xsolvecmi(k),sigmasmi(l))
       enddo      
    enddo
endif
if (nxi_new>0) then
    k = 0
    do i = ni_red+1, ni_red+nxi_new 
       k = k + 1
       l = 0
       do j = ni_red+1, ni_red+nxi_new
          l = l + 1
          E4(i,j)   =  mat_dotproduct(xsolvecpi(k),sigmaspi(l))
       enddo      
    enddo
endif


if (nx_new>0) then
    k = 0
    do i = n_red+1, n_red+nx_new 
       k = k + 1
    if (nm_new>0) then
       l = 0
       do j = nm_red+1, nm_red+nm_new
          l = l + 1
          S1(i,j)   =  mat_dotproduct(xsolvecp(k),rhosp(l))
       enddo      
    endif
   if (nmi_new>0) then
       l = 0
       do j = nmi_red+1, nmi_red+nmi_new
          l = l + 1
          S2(i,j)   =  mat_dotproduct(xsolvecp(k),rhospi(l))
       enddo      
    endif
    enddo
endif
if (nxi_new>0) then
    k = 0
    do i = ni_red+1, ni_red+nxi_new 
       k = k + 1
    if (nm_new>0) then
       l = 0
       do j = nm_red+1, nm_red+nm_new
          l = l + 1
          S7(i,j)   =  mat_dotproduct(xsolvecpi(k),rhosp(l))
       enddo      
   endif
   if (nmi_new>0) then
       l = 0
       do j = nmi_red+1, nmi_red+nmi_new
          l = l + 1
          S8(i,j)   =  mat_dotproduct(xsolvecpi(k),rhospi(l))
       enddo      
    endif
    enddo
endif

    !Setup lower half of E2 and S2:
    rewind(lu_sigma_rsp) 
    do j = 1, nt_red
       call mat_read_from_disk(lu_sigma_rsp,sigma_j,RSPonMaster)
       call mat_trans(sigma_j,scrT)
       call mat_add(0.5E0_realk,sigma_j,0.5E0_realk,scrT,sigmap_j)
       call mat_add(0.5E0_realk,sigma_j,-0.5E0_realk,scrT,sigmam_j)
       
       if ((nx_new>0) .and. (sqrt(mat_sqnorm2(sigmap_j))>1E-9_realk)) then
       k = 0
       do i = n_red+1, n_red+nx_new
          k = k + 1
          E1(i,j)   =  mat_dotproduct(xsolvecp(k),sigmap_j)
       enddo
      endif
      if ((nm_new>0) .and. (sqrt(mat_sqnorm2(sigmam_j))>1E-9_realk)) then
           k = 0
           do i = nm_red+1, nm_red+nm_new
              k = k + 1
             E2(i,j)   =  mat_dotproduct(xsolvecm(k),sigmam_j)
          enddo
      endif
   enddo
   
    rewind(lu_sigmai_rsp) 
    do j = 1, nti_red
       call mat_read_from_disk(lu_sigmai_rsp,sigmai_j,RSPonMaster)
       call mat_trans(sigmai_j,scrT)
       call mat_add(0.5E0_realk,sigmai_j,0.5E0_realk,scrT,sigmapi_j)
       call mat_add(0.5E0_realk,sigmai_j,-0.5E0_realk,scrT,sigmami_j)
       
       if ((nxi_new>0).and. (sqrt(mat_sqnorm2(sigmapi_j))>1E-9_realk)) then
       k = 0
       do i = ni_red+1, ni_red+nxi_new
          k = k + 1
          E4(i,j)   =  mat_dotproduct(xsolvecpi(k),sigmapi_j)
       enddo
      endif
       if ((nmi_new>0).and. (sqrt(mat_sqnorm2(sigmami_j))>1E-9_realk)) then
       k = 0
       do i = nmi_red+1, nmi_red+nmi_new
          k = k + 1
          E3(i,j)   =  mat_dotproduct(xsolvecmi(k),sigmami_j) 
       enddo
      endif
    enddo
    
    rewind(lu_rho_rsp)
    call mat_zero(rho_j)
    do j = 1, nt_red
       call mat_read_from_disk(lu_rho_rsp,rho_j,RSPonMaster)
       call mat_trans(rho_j,scrT)
       call mat_add(0.5E0_realk,rho_j,0.5E0_realk,scrT,rhop_j)
       
       if (sqrt(mat_sqnorm2(rhop_j))>1E-9_realk) then
      if (nx_new>0) then 
       k = 0
       do i = n_red+1, n_red+nx_new
          k = k + 1
          S1(i,j)   =  mat_dotproduct(xsolvecp(k),rhop_j)
       enddo
      endif
      if (nxi_new>0) then 
       k = 0
       do i = ni_red+1, ni_red+nxi_new
          k = k + 1
          S7(i,j)   =  mat_dotproduct(xsolvecpi(k),rhop_j)
       enddo
      endif
      endif
    enddo
    
    rewind(lu_rhoi_rsp)
    call mat_zero(rhoi_j)
    do j = 1, nti_red
       call mat_read_from_disk(lu_rhoi_rsp,rhoi_j,RSPonMaster)
       call mat_trans(rhoi_j,scrT)
       call mat_add(0.5E0_realk,rhoi_j,0.5E0_realk,scrT,rhopi_j)
       
       if (sqrt(mat_sqnorm2(rhopi_j))>1E-9_realk) then
     if (nx_new>0) then
          k = 0
         do i = n_red+1, n_red+nx_new
          k = k + 1
          S2(i,j)   =  mat_dotproduct(xsolvecp(k),rhopi_j)
       enddo
     endif
     if (nxi_new>0) then
          k = 0
         do i = ni_red+1, ni_red+nxi_new
          k = k + 1
          S8(i,j)   =  mat_dotproduct(xsolvecpi(k),rhopi_j) 
       enddo
     endif
     endif
    enddo
    
    !Explicitly calculate upper half of E2 and S2:
    rewind(lu_x_rsp); rewind(lu_xi_rsp)
    do i = 1, nt_red
       call mat_read_from_disk(lu_x_rsp,x_j,RSPonMaster)
       call mat_trans(x_j,scrT)
       call mat_add(0.5E0_realk,x_j,0.5E0_realk,scrT,xp_j)
       call mat_add(0.5E0_realk,x_j,-0.5E0_realk,scrT,xm_j)
    
  if (sqrt(mat_sqnorm2(xp_j))>1E-9_realk) then
    if (nx_new>0) then
       k = 0
       do j = n_red+1, n_red+nx_new  !Run over 2 elements at a time because of pairing
          k = k + 1
          E1(i,j)   =  mat_dotproduct(xp_j,sigmasp(k))
       enddo
     endif
     if (nm_new>0) then
       k = 0
       do j = nm_red+1, nm_red+nm_new
          k = k + 1
          S1(i,j)   =  mat_dotproduct(xp_j,rhosp(k))
       enddo
     endif
     if (nmi_new>0) then
       k = 0
       do j = nmi_red+1, nmi_red+nmi_new
          k = k + 1
          S2(i,j)   =  mat_dotproduct(xp_j,rhospi(k))
       enddo
     endif
     endif
    if ((nm_new>0) .and. (sqrt(mat_sqnorm2(xm_j))>1E-9_realk)) then
       k = 0
       do j = nm_red+1, nm_red+nm_new  !Run over 2 elements at a time because of pairing
          k = k + 1
          E2(i,j)   =  mat_dotproduct(xm_j,sigmasm(k))
       enddo
     endif
    enddo
    
    do i = 1, nti_red
       call mat_read_from_disk(lu_xi_rsp,xi_j,RSPonMaster)
       call mat_trans(xi_j,scrT)
       call mat_add(0.5E0_realk,xi_j,0.5E0_realk,scrT,xpi_j)
       call mat_add(0.5E0_realk,xi_j,-0.5E0_realk,scrT,xmi_j)
    
  if (sqrt(mat_sqnorm2(xpi_j))>1E-9_realk) then
    if (nxi_new>0) then
       k = 0
       do j = ni_red+1, ni_red+nxi_new  !Run over 2 elements at a time because of pairing
          k = k + 1
          E4(i,j)   =  mat_dotproduct(xpi_j,sigmaspi(k))
       enddo
     endif
      if (nm_new>0) then
       k = 0
       do j = nm_red+1, nm_red+nm_new
          k = k + 1
          S7(i,j)   =  mat_dotproduct(xpi_j,rhosp(k))
       enddo
       endif
      if (nmi_new>0) then
       k = 0
       do j = nmi_red+1, nmi_red+nmi_new
          k = k + 1
          S8(i,j)   =  mat_dotproduct(xpi_j,rhospi(k))
       enddo
       endif
       endif
    if ((nmi_new>0).and. (sqrt(mat_sqnorm2(xmi_j))>1E-9_realk)) then
       k = 0
       do j = nmi_red+1, nmi_red+nmi_new  !Run over 2 elements at a time because of pairing
          k = k + 1
          E3(i,j)   =  mat_dotproduct(xmi_j,sigmasmi(k))
       enddo
     endif
    enddo
 
    do j=1,ngd
    call mat_init(gdp,ndim,ndim)  
    call mat_init(gdm,ndim,ndim) 
    call mat_init(gT,ndim,ndim)
    call mat_trans(gd(j),gT)  
    
    call mat_add(0.5E0_realk,gd(j),0.5E0_realk,gT,gdp)
    call mat_add(0.5E0_realk,gd(j),-0.5E0_realk,gT,gdm)
 
    
    if (gd_complex) then 
    call mat_init(gdpi,ndim,ndim)  
    call mat_init(gdmi,ndim,ndim)
    call mat_zero(gT)
    call mat_trans(gdi(j),gT)
    call mat_add(0.5E0_realk,gdi(j),0.5E0_realk,gT,gdpi)
    call mat_add(0.5E0_realk,gdi(j),-0.5E0_realk,gT,gdmi)
   endif
   call mat_free(gT)
   
   enddo  
    do j = 1, ngd  !Loop over number of gradients
       k=0
       do i = n_red+1, n_red+nx_new !Loop in steps of 2 because of pairing
          k = k+1
           G1(i,j) = mat_dotproduct(xsolvecp(k),gdp)
       enddo
       k=0
       do i = nm_red+1, nm_red+nm_new !Loop in steps of 2 because of pairing
          k = k+1
           G2(i,j) = mat_dotproduct(xsolvecm(k),gdm)
       enddo
          
   if (gd_complex) then 
       k=0
       do i = nmi_red+1, nmi_red+nmi_new !Loop in steps of 2 because of pairing
          k = k+1
           G3(i,j) = mat_dotproduct(xsolvecmi(k),gdmi)
       enddo
       k=0
       do i = ni_red+1, ni_red+nxi_new !Loop in steps of 2 because of pairing
          k = k+1
           G4(i,j) = mat_dotproduct(xsolvecpi(k),gdpi)
       enddo
    endif 
       call mat_free(gdp)
       call mat_free(gdm)
       if (gd_complex) then
       call mat_free(gdpi)
       call mat_free(gdmi)

       endif
    enddo
        !save on disk
    if (molcfg%solver%info_rsp_redspace) then
          write(molcfg%lupri,*) 'Reduced E1:'
          call OUTPUT(E1, 1, n_red+nx_new, 1, n_red+nx_new, max_red, max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced E2:'
          call OUTPUT(E2, 1, nm_red+nm_new, 1, nm_red+nm_new, max_red, max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced E3:'
          call OUTPUT(E3, 1, nmi_red+nmi_new, 1, nmi_red+nmi_new, max_red, max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced E4:'
          call OUTPUT(E4, 1, ni_red+nxi_new, 1, ni_red+nxi_new, max_red, max_red, 1, molcfg%lupri)
          
          write(molcfg%lupri,*) 'Reduced S1:'
          call OUTPUT(S1, 1, n_red+nx_new, 1, nm_red+nm_new, max_red, max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced S2:'
          call OUTPUT(S2, 1, n_red+nx_new, 1, nmi_red+nmi_new, max_red, max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced S7:'
          call OUTPUT(S7, 1, ni_red+nxi_new, 1, nm_red+nm_new, max_red, max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced S8:'
          call OUTPUT(S8, 1, ni_red+nxi_new, 1, nmi_red+nmi_new, max_red, max_red, 1, molcfg%lupri)
          
          write(molcfg%lupri,*) 'Reduced G1:'
          call OUTPUT(G1, 1, n_red+nx_new, 1, ngd, max_red, max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced G2:'
          call OUTPUT(G2, 1, nm_red+nm_new, 1, ngd, max_red, max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced G3:'
          call OUTPUT(G3, 1, nmi_red+nmi_new, 1, ngd, max_red, max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced G4:'
          call OUTPUT(G4, 1, ni_red+nxi_new, 1, ngd, max_red, max_red, 1, molcfg%lupri)
      endif

call mat_free(scrT)      
call mat_free(x_j)
call mat_free(rho_j)
call mat_free(sigma_j)
call mat_free(xi_j)
call mat_free(rhoi_j)
call mat_free(sigmai_j)
call mat_free(xp_j)
call mat_free(sigmap_j)
call mat_free(xm_j)
call mat_free(rhop_j)
call mat_free(sigmam_j)
call mat_free(xpi_j)
call mat_free(sigmapi_j)
call mat_free(xmi_j)
call mat_free(rhopi_j)
call mat_free(sigmami_j)
    end subroutine extend_new_complex_reduced_matrices

!================================================================  
 !================================================================  
subroutine solve_complex(molcfg,ndim,ndim_red,nm_red,ni_red,nmi_red,ngd,freq,red_X,red_Xi,red_Xm,red_Xmi,gd_complex)
!
!  Solve reduced linear response problem in subspace
!  we consider one gradient/solution at a time
!
!  INPUT:
!  ndim_red: half size of g^R reduced space (E1)
!  nm_red: half size of  u^R reduced space (E2)
!  nmi_red: half size of  u^I reduced space (E3)
!  ni_red: half size of  g^I reduced space (E4)
!  ngd:      number of right hand sides
!  freq:     the ngd frequencies
!
!  OUTPUT:
!  red_X:    the ngd reduced space real g solution vectors
!  red_Xm:    the ngd reduced space real u solution vectors
!  red_Xmi:    the ngd reduced space img u solution vectors
!  red_Xi:    the ngd reduced space img g solution vectors
    implicit none
    type(rsp_molcfg), intent(inout) :: molcfg
    integer, intent(in)      :: ndim_red, ngd,ndim,nm_red
    integer,intent(in)       :: ni_red, nmi_red
    real(realk),intent(in)   :: freq(:)
    logical                  :: gd_complex
    real(realk),intent(inout):: red_X(:,:),red_Xm(:,:)
    real(realk),intent(inout):: red_Xi(:,:),red_Xmi(:,:)
    real(realk),allocatable  :: S1p(:,:), S2p(:,:)
    real(realk),allocatable  :: S3m(:,:),S4m(:,:)
    real(realk),allocatable  :: S5p(:,:), S7m(:,:)
    real(realk),allocatable  :: S6p(:,:), S8m(:,:)
    real(realk),allocatable   :: A(:,:),KHS(:)
    integer,allocatable      :: IPIV(:)
    integer                  :: igd, ierr, i, j, k,n_red,max_red
    real(realk)              :: gammma
    ierr=0
    gammma=molcfg%solver%rsp_gamma
    max_red=molcfg%solver%rsp_maxred
    n_red=ndim_red+nm_red+ni_red+nmi_red

    allocate(S2p(ndim_red,nmi_red), S5p(nmi_red,ndim_red))
    allocate(S1p(ndim_red,nm_red), S3m(nm_red,ndim_red))
    allocate(S4m(nm_red,ni_red), S7m(ni_red,nm_red))
    allocate(S6p(nmi_red,ni_red), S8m(ni_red,nmi_red))
    allocate(KHS(n_red), IPIV(n_red))
    allocate(A(n_red,n_red))
    S1p=0E0_realk; S2p=0E0_realk; S3m=0E0_realk; S4m=0E0_realk
    S5p=0E0_realk; S6p=0E0_realk; S7m=0E0_realk; S8m=0E0_realk
    
    KHS=0E0_realk
    IPIV=0E0_realk
    A=0E0_realk

    !Setup reduced E2, S2, and right hand side with proper dimension
    S2p(1:ndim_red,1:nmi_red) = S2(1:ndim_red,1:nmi_red)
    S7m(1:nm_red,1:ni_red) = S7(1:nm_red,1:ni_red)
    S1p(1:ndim_red,1:nm_red) = S1(1:ndim_red,1:nm_red)
    S8m(1:nmi_red,1:ni_red) = S8(1:nmi_red,1:ni_red)
     
       S2p = gammma*S2p
       S7m = gammma*S7m
   
       do i=1,nmi_red
       do j=1,ndim_red
           S5p(i,j) = S2p(j,i)
       enddo
    enddo
   do i=1,nm_red
       do j=1,ni_red
           S4m(i,j) = S7m(j,i)
       enddo
    enddo

       A(1:ndim_red,1:ndim_red)                                                              =  E1(1:ndim_red,1:ndim_red)
       A(ndim_red+1:ndim_red+nm_red,ndim_red+1:ndim_red+nm_red)                              =  E2(1:nm_red,1:nm_red)
       A(ndim_red+nm_red+1:ndim_red+nm_red+nmi_red,ndim_red+nm_red+1:ndim_red+nm_red+nmi_red)= -E3(1:nm_red,1:nm_red)
       A(ndim_red+nm_red+nmi_red+1:n_red,ndim_red+nm_red+nmi_red+1:n_red)                    = -E4(1:ndim_red,1:ndim_red)

       A(1:ndim_red,ndim_red+nm_red+1:ndim_red+nm_red+nmi_red)                            =  S2p(1:ndim_red,1:nmi_red)
       A(ndim_red+nm_red+1:ndim_red+nm_red+nmi_red,1:ndim_red)                            =  S5p(1:nmi_red,1:ndim_red)
       A(ndim_red+1:ndim_red+nm_red,ndim_red+nm_red+nmi_red+1:n_red)                      =  S4m(1:nm_red,1:ni_red)
       A(ndim_red+nm_red+nmi_red+1:n_red,ndim_red+1:ndim_red+nm_red)                      =  S7m(1:ni_red,1:nm_red)
       
       do igd=1,ngd
       S1p = - freq(igd)*S1p
       S8m =  freq(igd)*S8m
   do i=1,nm_red
       do j=1,ndim_red
           S3m(i,j) = S1p(j,i)
       enddo
    enddo
   do i=1,nmi_red
       do j=1,ni_red
           S6p(i,j) = S8m(j,i)
       enddo
    enddo
    

       A(1:ndim_red,ndim_red+1:ndim_red+nm_red)                                           =  S1p(1:ndim_red,1:nm_red)
       A(ndim_red+1:ndim_red+nm_red,1:ndim_red)                                           =  S3m(1:nm_red,1:ndim_red)
       A(ndim_red+nm_red+1:ndim_red+nm_red+nmi_red,ndim_red+nm_red+nmi_red+1:n_red)       =  S6p(1:nmi_red,1:ni_red)
       A(ndim_red+nm_red+nmi_red+1:n_red,ndim_red+nm_red+1:ndim_red+nm_red+nmi_red)       =  S8m(1:ni_red,1:nmi_red)

       KHS(1:ndim_red)                                                                    =  G1(1:ndim_red,igd)
       KHS(ndim_red+1:ndim_red+nm_red)                                                    =  G2(1:nm_red,igd)
       if (gd_complex) then
           KHS(ndim_red+nm_red+1:ndim_red+nm_red+nmi_red)                                 = -G3(1:nmi_red,igd)
           KHS(ndim_red+nm_red+nmi_red+1:n_red)                                           = -G4(1:ni_red,igd)
       endif       

       if (molcfg%solver%info_rsp_redspace) then
          write(molcfg%lupri,*) 'Reduced A:'
          call OUTPUT(A, 1, n_red, 1, n_red, 2*max_red, 2*max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced RHS:'
          call OUTPUT(KHS, 1, n_red, 1, 1, 2*max_red, 1, 1, molcfg%lupri)
       endif
         
       !Solve set of linear equations Ax = b:
       call DGESV(n_red, 1, A, n_red, IPIV, KHS, n_red, IERR)
       !Solution vector is found in RHS.
       if (IERR /= 0) then
          WRITE(molcfg%LUPRI,'(/A, i4)') &
          &     'Problem in DGESV, IERR = ', IERR
          CALL lsQUIT(' Problem in DGESV',molcfg%lupri)
       endif
       
       red_X(1:ndim_red,igd) = KHS(1:ndim_red)
       red_Xm(1:nm_red,igd) = KHS(ndim_red+1:ndim_red+nm_red)
       red_Xmi(1:nmi_red,igd) = KHS(ndim_red+nm_red+1:ndim_red+nm_red+nmi_red)
       red_Xi(1:ni_red,igd) = KHS(ndim_red+nm_red+nmi_red+1:n_red)
    enddo

    deallocate(S1p, S2p)
    deallocate(S4m, S7m)
    deallocate(S5p, S3m)
    deallocate(S6p, S8m)
    deallocate(KHS, A)
    deallocate(IPIV)
  end subroutine solve_complex
  
  
  subroutine get_complex_res(molcfg,F,D,S,itmic,ndim_red,nm_red,ni_red,nmi_red,nt_red,nti_red,ngd,gd,red_Xp_glob,&
             &red_Xpi_glob,red_Xm_glob, red_Xmi_glob,eival,xp,xm,xpi,xmi,nx_new,conv,res_norm_tot,gd_complex,gdi)
  use RSPSYMSOLVER
    implicit none
!******************************************************************
! Purpose: 
! 1)Construct residuals:
!   Rg^R=E(2)Xg^R(I)-W(I)*S(2)Xu^R(I) - GDg^R + gamma*S[2]Xu^I
!   Ru^R=E(2)Xu^R-W*S(2)Xg^R - GDu^R + gamma*S[2]Xg^I
!   Ru^I=-E(2)Xu^I+W*S(2)Xg^I + GDu^I + gamma*S[2]Xg^R
!   Rg^I=-E(2)Xg^I+W*S(2)Xu^I + GDg^I + gamma*S[2]Xu^R
!   for NGD rsp vectors X(I) of reduced rsp problem
! 2)Test for convergence of ngd eigenvectors,
!   Convergence criterium:
!   ||Rg^R||+||Rg^I|| < rsp_thresh
!   ||Ru^R||+||Ru^I|| < rsp_thresh 
! 3) contruct new trial vectors as preconditioned residuals
!
!  INPUT:
!  itmic:    number of the microiteration (for printout)
!  ndim_red: half size of g^R reduced space (E1)
!  nm_red:   half size of u^R reduced space (E2)
!  ni_red: half size of g^I reduced space (E4)
!  nmi_red:   half size of u^I reduced space (E3)
!  nt_red:   number of vectors real stored on disk (vectors are stored as g+u, 
!            in case when the trial vector is either symmetric or antisymmetric
!            the total number of trial vectors may be different than number of g/u)
!  nti_red:   number of vectors img stored on disk (vectors are stored as g+u)
!  ngd:      number of right hand sides
!  gd(ngd):  matrix array containing ngd real right hand sides (gradients)
!  gdi(ngd):  matrix array containing ngd img right hand sides (gradients)
!  red_Xp_glob:    the ngd reduced space g^R solution vectors (same as RED_X)
!  red_Xm_glob:    the ngd reduced space u^R solution vectors (same as RED_Xm)
!  red_Xpi_glob:    the ngd reduced space g^I solution vectors (same as RED_Xi)
!  red_Xmi_glob:    the ngd reduced space u^I solution vectors (same as RED_Xmi)
!  eival:    the ngd frequencies
!  
!  OUTPUT:
!  nx_NEW:   number of new g^R vectors
!  nm_NEW:   number of new u^R vectors
!  nxi_NEW:   number of new g^I vectors
!  nmi_NEW:   number of new u^I vectors
!  xp(nx_new) new g^R trial vector
!  xm(nm_new) new u^R trial vector
!  xpi(nxi_new) new g^I trial vector
!  xmi(nmi_new) new u^I trial vector
!  conv:     set true if converged
! res_norm_tot : norm of residual to check if convergence not "stuck"
!******************************************************************

    type(rsp_molcfg), intent(inout)  :: molcfg
    integer,intent(in)               :: itmic,ndim_red,ngd,nm_red,ni_red,nmi_red,nt_red,nti_red
    type(Matrix),intent(in)          :: F,D,S
    type(Matrix),intent(in)          :: gd(:)
    type(Matrix),intent(in),optional :: gdi(:)
    logical,intent(in)               :: gd_complex
    real(realk), intent(in)          :: red_Xp_glob(:,:)
    real(realk), intent(inout)          :: eival(:)
    real(realk), intent(in)          :: red_Xpi_glob(:,:)
    real(realk), intent(in)          :: red_Xm_glob(:,:)
    real(realk), intent(in)          :: red_Xmi_glob(:,:)
    type(Matrix),intent(inout)       :: xp(:),xm(:),xpi(:),xmi(:)
    logical, intent(out)             :: conv
    real(realk),intent(out)          :: res_norm_tot
    integer, intent(out)             :: nx_new
!local
    type(Matrix)                     :: gdp,gdpi,gdm,gdmi,gT
    logical                          :: conv_root, conv_root_img
    type(Matrix)                     :: S2XR(ngd), S2XI(ngd),S2XRi(ngd),S2XIi(ngd),xr,xi
    type(Matrix)                     :: residualsp(ngd),residualsm(ngd)
    type(Matrix)                     :: residualspi(ngd),residualsmi(ngd)
    real(realk)                      :: wibx,res_norm,x_red_norm,conv_test,ddot,diff,x_red_norm_img
    integer                          :: i,n_not_conv,ndim_red_mat,ndim,index_conv(ngd),ibx,j,k,max_red
    real(realk)                      :: a, b,c,rsp_conv_thr,rsp_thresh
    real(realk)                      :: av_norm,res_norm_img, conv_test_img,gammma

    gammma=molcfg%solver%rsp_gamma
    max_red=molcfg%solver%rsp_maxred
    rsp_thresh= molcfg%solver%rsp_thresh

    ndim = S%nrow
    n_not_conv = 0
    conv=.false.
 
    do i=1,ngd
       call mat_init(S2XR(i),ndim,ndim)
       call mat_init(S2XI(i),ndim,ndim)
       call mat_init(S2XRi(i),ndim,ndim)
       call mat_init(S2XIi(i),ndim,ndim)
       call mat_zero(S2XR(i))
       call mat_zero(S2XI(i))
       call mat_zero(S2XRi(i))
       call mat_zero(S2XIi(i))
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
       call mat_init(residualsp(ibx),ndim,ndim)
       call mat_init(residualsm(ibx),ndim,ndim)
       call mat_init(residualspi(ibx),ndim,ndim)
       call mat_init(residualsmi(ibx),ndim,ndim)
       call mat_zero(residualsp(ibx))
       call mat_zero(residualsm(ibx))
       call mat_zero(residualspi(ibx))
       call mat_zero(residualsmi(ibx))
       
       rewind(lu_rho_rsp)
       rewind(lu_sigma_rsp)
       rewind(lu_sigmai_rsp)
       rewind(lu_rhoi_rsp)
call expand_on_basis4(molcfg,nt_red,ndim,red_Xm_glob(1:nm_red,ibx),red_Xp_glob(1:ndim_red,ibx),lu_rho_rsp,&
                                  &residualsp(ibx),residualsm(ibx))
call expand_on_basis4(molcfg,nti_red,ndim,red_Xmi_glob(1:nmi_red,ibx),red_Xpi_glob(1:ni_red,ibx),lu_rhoi_rsp,&
                                  &residualspi(ibx),residualsmi(ibx))
       
       call mat_assign(S2XR(ibx),residualsm(ibx))
       call mat_assign(S2XRi(ibx),residualsmi(ibx))
       call mat_assign(S2XI(ibx),residualsp(ibx))
       call mat_assign(S2XIi(ibx),residualspi(ibx))
       
       call mat_scal(-wibx,residualsp(ibx))
       call mat_scal(-wibx,residualsm(ibx))
       call mat_scal(-wibx,residualspi(ibx))
       call mat_scal(-wibx,residualsmi(ibx))
       !add  sum_j ^RX_j(ibx) sigma_j
       
 call expand_on_basis4(molcfg,nt_red,ndim,red_Xp_glob(1:ndim_red,ibx),red_Xm_glob(1:nm_red,ibx),lu_sigma_rsp,&
                       &residualsp(ibx),residualsm(ibx))
 call expand_on_basis4(molcfg,nti_red,ndim,red_Xpi_glob(1:ni_red,ibx),red_Xmi_glob(1:nmi_red,ibx),lu_sigmai_rsp,&
                       &residualspi(ibx),residualsmi(ibx))

       call mat_scal(-1E0_realk,residualspi(ibx))
       call mat_scal(-1E0_realk,residualsmi(ibx))
      ! subtract gradient
    call mat_init(gdp,ndim,ndim)  
    call mat_init(gdm,ndim,ndim) 
    call mat_init(gT,ndim,ndim)
    call mat_trans(gd(ibx),gT)  
    
    call mat_add(0.5E0_realk,gd(ibx),0.5E0_realk,gT,gdp)
    call mat_add(0.5E0_realk,gd(ibx),-0.5E0_realk,gT,gdm)
   
    call mat_init(gdpi,ndim,ndim)  
    call mat_init(gdmi,ndim,ndim)
    call mat_zero(gdpi)
    call mat_zero(gdmi)
   
    if (gd_complex) then 
    call mat_zero(gT)
    call mat_trans(gdi(ibx),gT)
    call mat_add(0.5E0_realk,gdi(ibx),0.5E0_realk,gT,gdpi)
    call mat_add(0.5E0_realk,gdi(ibx),-0.5E0_realk,gT,gdmi)
   endif
   
   call mat_free(gT)
   
      call mat_daxpy(gammma,S2XIi(ibx),residualsp(ibx))
      call mat_daxpy(gammma,S2XRi(ibx),residualsm(ibx))
      call mat_daxpy(gammma,S2XI(ibx),residualspi(ibx))
      call mat_daxpy(gammma,S2XR(ibx),residualsmi(ibx))
  
      call mat_daxpy(-1.0E0_realk,gdp,residualsp(ibx))
      call mat_daxpy(-1E0_realk,gdm,residualsm(ibx))
      if (gd_complex) then
      call mat_daxpy(1E0_realk,gdpi,residualspi(ibx))
      call mat_daxpy(1E0_realk,gdmi,residualsmi(ibx))
      endif
      
call mat_free(gdp)
call mat_free(gdm)
call mat_free(gdpi)
call mat_free(gdmi)
enddo
! the residual(s) is now done: test for convergence. If not converged
! form new trial(s)  in next subroutine
! New trial(s) is equal to preconditioned residual(s)
   if (molcfg%solver%rsp_convdyn) then
       if (ITMIC == 1) then
          molcfg%solver%rsp_dyn_thresh=0E0_realk
          av_norm = 0.0E0_realk
          do i = 1, ngd
             res_norm     = sqrt(mat_sqnorm2(residualsp(i)))+sqrt(mat_sqnorm2(residualspi(i)))
             res_norm_img = sqrt(mat_sqnorm2(residualsm(i)))+sqrt(mat_sqnorm2(residualsmi(i)))
             av_norm = av_norm + res_norm + res_norm_img
          enddo
          av_norm = 0.25E0_realk*av_norm
          rsp_conv_thr = av_norm*molcfg%solver%rsp_conv_factor
          molcfg%solver%rsp_dyn_thresh=max(rsp_conv_thr,rsp_thresh)
          rsp_thresh= molcfg%solver%rsp_dyn_thresh
          WRITE(molcfg%LUPRI,*) 'Dynamic response convergence threshold set to', rsp_thresh 
       else
          rsp_thresh= molcfg%solver%rsp_dyn_thresh
       endif
   else
          rsp_thresh= molcfg%solver%rsp_thresh
   endif

   do ibx = 1,ngd
      conv_root = .false.
      conv_root_img=.false.
      if (molcfg%solver%rsp_single_norm) then
          res_norm=mat_sqnorm2(residualsp(ibx))+mat_sqnorm2(residualspi(ibx))
         res_norm_img=mat_sqnorm2(residualsm(ibx))+mat_sqnorm2(residualsmi(ibx))
      else ! double norm
         res_norm    =sqrt(mat_sqnorm2(residualsp(ibx)))+sqrt(mat_sqnorm2(residualspi(ibx)))
         res_norm_img=sqrt(mat_sqnorm2(residualsm(ibx)))+sqrt(mat_sqnorm2(residualsmi(ibx)))
      endif
      if (res_norm < rsp_thresh) conv_root = .true.
      if (res_norm_img < rsp_thresh) conv_root_img = .true.
      res_norm_tot=res_norm+res_norm_img

     ! Insert back print statements
      write (molcfg%LUPRI, '("Residual norm for herm vector ", i3, " is: ", E14.6, "    and frequency = " , F12.6, "   It = ", i3, &
                        & " CONV =", l2)') ibx, res_norm, EIVAL(ibx), itmic, conv_root
    
      write (molcfg%LUPRI,'(A,I3,A,E14.6,A,F12.6,A,I3,A,l2)') &
           & 'Residual norm for antiherm vector',ibx,' is: ',res_norm_img,&
           & '    and frequency = ',EIVAL(ibx),'   It = ',itmic,' CONV =',conv_root_img
         
         conv_test = res_norm
         conv_test_img = res_norm_img
      if (molcfg%solver%info_rsp) then
        write(molcfg%lupri,*) 'conv_test',conv_test
        write(molcfg%lupri,*) 'res_norm',res_norm
        write(molcfg%lupri,*) 'rsp_thresh',rsp_thresh
        write(molcfg%lupri,*) 'conv_test_img',conv_test_img
        write(molcfg%lupri,*) 'res_norm_img',res_norm_img
      endif
      if (conv_root .and. conv_root_img) then 
          conv=.true.
      endif
   enddo  
   nx_new=0
if (.not. conv) then     
      do j=1,ngd
    call mat_init(xr,ndim,ndim)  
    call mat_init(xi,ndim,ndim) 
    call mat_init(gT,ndim,ndim)

          nx_new=nx_new+1
          call new_complex_precond(residualsp(j),residualsm(j),residualspi(j),residualsmi(j),S,-eival(j),-gammma, &
               & xp(nx_new),xm(nx_new),xpi(nx_new),xmi(nx_new),molcfg%decomp%nocc) 

          call util_scriptPx('N',D,S,xp(nx_new))
          call util_scriptPx('N',D,S,xm(nx_new))
          call util_scriptPx('N',D,S,xpi(nx_new))
          call util_scriptPx('N',D,S,xmi(nx_new))
     call mat_scal(-1E0_realk,xpi(nx_new))
     call mat_scal(-1E0_realk,xmi(nx_new))
    
      call mat_free(gT)
      call mat_free(xr)
      call mat_free(xi)
      enddo
endif
          
! Finished, deallocate local arrays
   do ibx=1,ngd
        call mat_free(residualsp(ibx))
        call mat_free(residualsm(ibx))
        call mat_free(residualspi(ibx))
        call mat_free(residualsmi(ibx))
        call mat_free(S2XR(ibx))
        call mat_free(S2XI(ibx))
        call mat_free(S2XRi(ibx))
        call mat_free(S2XIi(ibx))
   enddo

  end subroutine get_complex_res
!================================================================  

  subroutine get_start_vectors(molcfg,F,D,S,RHS_real,RHS_img,ngd,eival,nx_new,xp,xm,xpi,xmi)
    implicit none
!================================================================ 
! A subroutine to get start vectors for damped response eqs.
! by preconditioning the RHS vectors.
! g and u components obtianed directly.
!================================================================ 
    type(rsp_molcfg), intent(inout) :: molcfg
    type(Matrix),intent(in)    :: F,D,S
    type(Matrix),intent(inout) :: RHS_real(:),RHS_img(:)
    type(Matrix),intent(inout) :: xp(:),xm(:),xpi(:),xmi(:)
    real(realk), intent(inout) :: eival(:) 
    integer,intent(inout)      :: nx_new
    integer,intent(in)         :: ngd
    type(matrix)               :: gdp,gdm,gdpi,gdmi,gT
    integer                    :: j,ndim
    real(realk)                :: gammma

 ndim=D%nrow
 gammma=molcfg%solver%rsp_gamma

  do j=1,ngd
     nx_new=nx_new+1
      
     call mat_init(gT,ndim,ndim)
     call mat_init(gdp,ndim,ndim)
     call mat_init(gdm,ndim,ndim)
     call mat_init(gdpi,ndim,ndim)
     call mat_init(gdmi,ndim,ndim)
    
     call mat_trans(RHS_real(j),gT)  
     call mat_add(0.5E0_realk,RHS_real(j),0.5E0_realk,gT,gdp)
     call mat_add(0.5E0_realk,RHS_real(j),-0.5E0_realk,gT,gdm)
     call mat_trans(RHS_img(j),gT)  
     call mat_add(0.5E0_realk,RHS_img(j),0.5E0_realk,gT,gdpi)
     call mat_add(0.5E0_realk,RHS_img(j),-0.5E0_realk,gT,gdmi)

     call new_complex_precond(gdp,gdm,gdpi,gdmi,S,-eival(j),-gammma, &
               & xp(nx_new),xm(nx_new),xpi(nx_new),xmi(nx_new),molcfg%decomp%nocc) 

     call util_scriptPx('N',D,S,xp(nx_new))
     call util_scriptPx('N',D,S,xm(nx_new))
     call util_scriptPx('N',D,S,xpi(nx_new))
     call util_scriptPx('N',D,S,xmi(nx_new))

     call mat_scal(-1E0_realk,xpi(nx_new))
     call mat_scal(-1E0_realk,xmi(nx_new))
          
     call mat_free(gdp)
     call mat_free(gdm)
     call mat_free(gdpi)
     call mat_free(gdmi)
     call mat_free(gT)
 enddo
end subroutine get_start_vectors
 end module COMPLEXSYMSOLVER
