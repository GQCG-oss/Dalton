module RSPSYMSOLVER
  use precision
  use matrix_module
  use matrix_util
  use decompMod
  use matrix_operations
  use RSPSOLVER
  use rsp_util
  use files
  use memory_handling
  !use configuration
  private
  public ::  rsp_sym_solver,rsp_sym_init,expand_on_basis3,orthonormalize22,&
           &orthonormalize21,expand_on_basis4
  integer, save :: rsp_number_of_current_trial, rsp_number_of_rhs, rsp_number_of_sols, &
                 & rsp_number_of_omegas, rsp_number_of_startvecs,rsp_bvec_dim, &
                 & lu_x_rsp, lu_sigma_rsp, lu_rho_rsp
  real(realk),allocatable,save :: E1(:,:), E2(:,:), S1(:,:), S2(:,:), G1(:,:), G2(:,:)
  logical,save :: RSPonMaster
contains
  subroutine rsp_sym_init(ntrial, nrhs, nsol, nomega, nstart)
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
      RSPOnMaster=.TRUE.
  end subroutine rsp_sym_init
  subroutine rsp_sym_solver(molcfg,D,S,F,lineq,ngd,GDB,EIVAL,XSOL)
!==============================================================
 ! The response solver using the subspace algorithm with symmetrized trial vectors
 ! for solving the standard response equation and the eigenvalue equation
 ! (E[2]-wS[2])(X[g]+X[u])=GD[g]+GD[u]
 !
 ! input:
 ! lineq if .true. std response eq is solved
 !       if .false. eigenvalue eq is solved
 !
 ! if (lineq) then
 ! ngd   number of RHS else number of excitations
 ! GDB   RHS 
 ! EIVAL frequency
 ! else
 ! ngd   number of excitations
 !
 ! output:
 ! XSOL    solution vector
 ! if (lineq=.false.)
 ! eival excitantion energies
 !   
 ! Joanna, July 2010
 !=========================================================================
!
implicit none
type(rsp_molcfg), intent(inout)     :: molcfg
logical, intent(in)                 :: LINEQ
type(Matrix),intent(in)             :: F,D,S
type(Matrix),intent(inout)          :: GDB(rsp_number_of_rhs)
real(realk),intent(inout)           :: eival(rsp_number_of_omegas)
type(Matrix),intent(inout)          :: XSOL(rsp_number_of_sols)
integer, intent(in)                 :: ngd
type(Matrix)                        :: xp(rsp_bvec_dim),xm(rsp_bvec_dim)
real(realk)                         :: res_norm_tot,r1,r2,r3
integer                             :: i,j,ndim,k,n_red,nx_new,m,l,n_v,nm_new,nm_red,n_i,nt_red
logical                             :: conv
real(realk),allocatable             :: red_Xp_glob(:,:),red_Xm_glob(:,:)
real(realk)                         :: absorp, disper
type(Matrix)                        :: RHS_real(rsp_number_of_rhs),RHS_img(rsp_number_of_rhs)
type(Matrix)                        :: Xsolvecp(rsp_bvec_dim),sigmasp(rsp_bvec_dim)
type(Matrix)                        :: Xsolvec(rsp_bvec_dim)
type(Matrix),pointer                :: sigmas(:),rhos(:)
type(Matrix)                        :: rhosm(rsp_bvec_dim),rhosp(rsp_bvec_dim),gT
type(Matrix)                        :: Xsolvecm(rsp_bvec_dim),sigmasm(rsp_bvec_dim)
!debug
real(realk)                         :: rsp_thresh   
integer                             :: max_it, max_red
! call mem_alloc(sigmas,rsp_bvec_dim)
! call mem_alloc(rhos,rsp_bvec_dim)
 ndim = S%nrow
   lu_x_rsp= -1 ; lu_sigma_rsp = -1 ; lu_rho_rsp = -1
  CALL LSOPEN(lu_x_rsp,'xrsp','unknown','UNFORMATTED')
  CALL LSOPEN(lu_sigma_rsp,'sigma_rsp','unknown','UNFORMATTED')
  CALL LSOPEN(lu_rho_rsp,'rho_rsp','unknown','UNFORMATTED')
  
 max_it= molcfg%solver%rsp_maxit
 max_red= molcfg%solver%rsp_maxred
 rsp_thresh= molcfg%solver%rsp_thresh

   allocate(red_Xp_glob(max_red,ngd),red_Xm_glob(max_red,ngd))
   ALLOCATE(E1(max_red,max_red),E2(max_red,max_red))
   allocate(S1(max_red,max_red),S2(max_red,max_red))
   if (LINEQ) then
       !allocate space for reduced gradient
   allocate(G1(max_red,ngd))  
   allocate(G2(max_red,ngd))  
       G1=0E0_realk; G2=0E0_realk
   endif

  red_Xp_glob = 0.0E0_realk ; red_Xm_glob = 0.0E0_realk
  E1=0E0_realk; E2=0E0_realk; S1=0E0_realk; S2=0E0_realk
  
    if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'RSP_SOLVER: LINEQ',  LINEQ
    IF (LINEQ) THEN
       WRITE (molcfg%LUPRI,'(///2A//A/)') &
  &     ' <<<  SOLVING SETS OF LINEAR EQUATIONS ', &
  &     'FOR LINEAR RESPONSE PROPERTIES >>>'
    ELSE
       WRITE (molcfg%LUPRI,'(///2A//A/)') &
  &     ' <<< EXCITATION ENERGIES', &
  &     ' AND TRANSITION MOMENT CALCULATION (LSTDHF) >>>'
    END IF
   
   n_red=0
   nm_red=0
   nt_red=0
   conv=.false.
   nx_new=0
   nm_new=0

   do j=1,ngd
      call mat_init(xp(j),ndim,ndim)
      call mat_init(xm(j),ndim,ndim)
      call mat_zero(xp(j))
      call mat_zero(xm(j))
      if (lineq) then
      call mat_init(RHS_real(j),ndim,ndim)
      call mat_init(RHS_img(j),ndim,ndim)
      call mat_zero(RHS_real(j))
      call mat_zero(RHS_img(j))
      endif
   enddo
   if (lineq) then
     do j=1,ngd
         call mat_init(gT,ndim,ndim)
         call mat_trans(gdb(j),gT)  
    
         call mat_add(0.5E0_realk,gdb(j),0.5E0_realk,gT,RHS_real(j))
         call mat_add(0.5E0_realk,gdb(j),-0.5E0_realk,gT,RHS_img(j))
         call mat_free(gT)
      enddo
    endif

    IF (LINEQ) THEN              
       if (molcfg%solver%info_rsp) WRITE(molcfg%LUPRI,*) 'Now enter 1st guess for linear equations '
   call get_start_vectors(molcfg,F,D,S,RHS_real,RHS_img,ngd,eival,nx_new,nm_new,xp,xm)
    ELSE
       if (molcfg%solver%info_rsp) WRITE(molcfg%LUPRI,*) 'Now enter 1st guess for excitation energies'
       call get_1st_orth_trial_eigen(molcfg,D,S,ngd,xp,xm,nx_new,nm_new)
    END IF


    do i = 1,max_it
      write(molcfg%lupri,*) '------------------------------------------'
      write(molcfg%lupri,*) 'Start iteration in the rsp_std_sym solver:',i
      write(molcfg%lupri,*) '------------------------------------------'
       if (n_red > max_red) then
          WRITE(molcfg%LUPRI,'(/A)') &
          &     'rsp_maxred too small - recompile with larger rsp_maxred parameter'
          CALL lsQUIT('rsp_maxred parameter too small',molcfg%lupri)
       endif

       n_v=nx_new 
       n_i=nm_new
       do j=1,nx_new
          call mat_init(xsolvecp(j),ndim,ndim)
          call mat_init(sigmasp(j),ndim,ndim)
          call mat_init(rhosm(j),ndim,ndim)
       enddo
       do j=1,nm_new
          call mat_init(xsolvecm(j),ndim,ndim)
          call mat_init(sigmasm(j),ndim,ndim)
          call mat_init(rhosp(j),ndim,ndim)
        enddo
        do j=1,max(nx_new,nm_new)
            call mat_init(xsolvec(j),ndim,ndim)
!            call mat_init(sigmas(j),ndim,ndim)
!            call mat_init(rhos(j),ndim,ndim)
!            call mat_zero(sigmas(j))
!            call mat_zero(rhos(j))
        enddo
         if (lineq) then
             call orthonormalize21(molcfg,D,S,xp,xm,nx_new,nm_new,nt_red,xsolvecp,xsolvecm,lu_x_rsp,lineq)
         else 
            do j=1,nx_new
               call mat_add(0.5E0_realk,xp(j),0.5E0_realk,xm(j),xsolvec(j))
            enddo
               call orthonormalize22(molcfg,D,S,xsolvec,nx_new,nt_red,xsolvecp,lu_x_rsp)
               if (nx_new>0) then
                 do j=1,nx_new
                   call mat_init(gT,ndim,ndim)
                   call mat_trans(xsolvecp(j),gT)
                   call mat_add(0.5E0_realk,xsolvecp(j),0.5E0_realk,gT,xp(j))
                   call mat_add(0.5E0_realk,xsolvecp(j),-0.5E0_realk,gT,xm(j))
                   call mat_free(gT)
                 enddo
                 call orthonormalize21(molcfg,D,S,xp,xm,nx_new,nm_new,nt_red,xsolvecp,xsolvecm,lu_x_rsp,lineq)
              !   nm_new=nx_new
                 !debug
      !           do j=1,nx_new
      !              call mat_zero(xsolvecp(j))
      !              call mat_zero(xsolvecm(j))
      !              xsolvecp(j)=xp(j)
      !              xsolvecm(j)=xm(j)
      !           enddo
                 !debug
               else
                 write(molcfg%lupri,*) 'No new trial vectors'
                  exit        
                endif
         endif


       !  if (lineq) then 
          if ((nx_new>0) .and. (nm_new)>0) then    
            do j=1,min(nx_new,nm_new)
               call mat_add(1E0_realk,xsolvecp(j),1E0_realk,xsolvecm(j),xsolvec(j))
            enddo
            if (nx_new>nm_new) then
               do j=nm_new+1,nx_new
                   call mat_assign(xsolvec(j),xsolvecp(j))
               enddo
            elseif (nm_new>nx_new) then
               do j=nx_new+1,nm_new
                  call mat_assign(xsolvec(j),xsolvecm(j))
               enddo
             endif
          else
            if (nm_new==0) then 
               do j=1,nx_new
                   call mat_assign(xsolvec(j),xsolvecp(j))
               enddo
            elseif (nx_new==0) then
               do j=1,nm_new
                   call mat_assign(xsolvec(j),xsolvecm(j))
               enddo
           endif
          endif
        !  else
         !    do j=1,nx_new
         !       call mat_add(1E0_realk,xsolvecp(j),1E0_realk,xsolvecm(j),xsolvec(j))
         !    enddo
         !  endif
           call transform_vectors(molcfg,D,S,F,max(nx_new,nm_new),xsolvec,sigmas,rhos,.true.)
           do j=1,min(nx_new,nm_new)
              call mat_init(gT,ndim,ndim)
              call mat_trans(sigmas(j),gT)  
              call mat_add(0.5E0_realk,sigmas(j),0.5E0_realk,gT,sigmasp(j))
              call mat_add(0.5E0_realk,sigmas(j),-0.5E0_realk,gT,sigmasm(j))
    
              call mat_trans(rhos(j),gT)  
              call mat_add(0.5E0_realk,rhos(j),0.5E0_realk,gT,rhosp(j))
              call mat_add(0.5E0_realk,rhos(j),-0.5E0_realk,gT,rhosm(j))
              call mat_free(gT)
           enddo
          ! if (lineq) then
          if ((nx_new>0) .and. (nm_new)>0) then    
              if (nx_new>nm_new) then
                 do j=nm_new+1,nx_new
                    call mat_assign(sigmasp(j),sigmas(j))
                    call mat_assign(rhosm(j),rhos(j))
                 enddo
               elseif (nm_new>nx_new) then
                 do j=nx_new+1,nm_new
                      call mat_assign(sigmasm(j),sigmas(j))
                      call mat_assign(rhosp(j),rhos(j))
                 enddo
               endif 
           else
            if (nm_new==0) then 
               do j=1,nx_new
                   call mat_assign(sigmasp(j),sigmas(j))
                   call mat_assign(rhosm(j),rhos(j))
               enddo
            elseif (nx_new==0) then
               do j=1,nm_new
                   call mat_assign(sigmasm(j),sigmas(j))
                   call mat_assign(rhosp(j),rhos(j))
               enddo
           endif
           endif
       if (lineq) then
           call extend_new_reduced_matrices(molcfg,n_red,nm_red,nt_red,nx_new,nm_new,ngd,&
                &xsolvecp,sigmasp,rhosm, xsolvecm,sigmasm,rhosp,.true.,RHS_real,RHS_img)
       else
           call extend_new_reduced_matrices(molcfg,n_red,nm_red,nt_red,nx_new,nm_new,ngd,&
                &xsolvecp,sigmasp,rhosm, xsolvecm,sigmasm,rhosp,.false.)
       endif

        !save on disk
           do j = 1,max(nx_new,nm_new)
              call mat_write_to_disk(lu_x_rsp,xsolvec(j),RSPonMaster)
              call mat_write_to_disk(lu_sigma_rsp,sigmas(j),RSPonMaster)
              call mat_write_to_disk(lu_rho_rsp,rhos(j),RSPonMaster)
           enddo
    
           n_red = n_red + nx_new !update reduced space dimension
           nm_red = nm_red + nm_new !update reduced space dimension
           nt_red=nt_red+max(nx_new,nm_new)
   
        if (lineq) then
        !solve complex reduced equation
            call solve_std(molcfg,ndim,n_red,nm_red,ngd,eival,red_Xp_glob,red_Xm_glob)
        else
           call solve_red_eigen(molcfg,ndim,n_red,nm_red,ngd,eival,red_Xp_glob,red_Xm_glob)
        endif

        ! check for convergece and get new trialvectors
        do j=1,ngd
           call mat_zero(xp(j))
           call mat_zero(xm(j))
        enddo
if (lineq) then
        call get_new_res(molcfg,F,D,S,i,n_red,nm_red,nt_red,ngd,red_Xp_glob,red_Xm_glob,eival, &
             &xp,xm,nx_new,nm_new,conv,res_norm_tot,.true.,RHS_real,RHS_img)
 else
        call get_new_res(molcfg,F,D,S,i,n_red,nm_red,nt_red,ngd,red_Xp_glob,red_Xm_glob,eival, &
             &xp,xm,nx_new,nm_new,conv,res_norm_tot,.false.)
endif

         if (conv)  then
            do l=1,ngd
               rewind(lu_x_rsp)
               call expand_on_basis4(molcfg,nt_red,ndim,RED_Xp_glob(1:n_red,l),RED_Xm_glob(1:nm_red,l),&
                    &lu_x_rsp,Xp(l),Xm(l))
               call mat_zero(xsol(l))
               call mat_add(1E0_realk,xp(l),1E0_realk,xm(l),xsol(l))
            enddo
           do j=1,n_v
              call mat_free(xsolvecp(j))
              call mat_free(sigmasp(j))
              call mat_free(rhosm(j))
           enddo
           do j=1,n_i
             call mat_free(xsolvecm(j))
             call mat_free(sigmasm(j))
             call mat_free(rhosp(j))
           enddo
           do j=1,max(n_v,n_i)
              call mat_free(xsolvec(j))
              call mat_free(sigmas(j))
              call mat_free(rhos(j))
          enddo
          call mem_dealloc(sigmas)
          call mem_dealloc(rhos)
         do j=1,ngd
            call mat_free(xp(j))
            call mat_free(xm(j))
         enddo

           exit
    else         
        do l=1,ngd
           rewind(lu_x_rsp)
           call expand_on_basis4(molcfg,nt_red,ndim,RED_Xp_glob(1:n_red,l),RED_Xm_glob(1:nm_red,l),lu_x_rsp,&
                &Xp(l),Xm(l))
           call mat_zero(xsol(l))
           call mat_add(1E0_realk,xp(l),1E0_realk,xm(l),xsol(l))
          
           if (lineq) then 
              absorp=mat_dotproduct(gdb(l),xsol(l))
              write(molcfg%lupri,*) "Dispersion value macroiteration:",i, "for &
               & frequency",eival(l),"=", absorp
           endif
       enddo
   
  if (lineq) then
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
               do j=1,n_i
                  call mat_free(xsolvecm(j))
                  call mat_free(sigmasm(j))
                  call mat_free(rhosp(j))
               enddo
               do j=1,max(n_v,n_i)
                  call mat_free(xsolvec(j))
                  call mat_free(sigmas(j))
                  call mat_free(rhos(j))
                enddo
                call mem_dealloc(sigmas)
                call mem_dealloc(rhos)                
                do j=1,ngd
                   call mat_free(xp(j))
                   call mat_free(xm(j))
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
                 do j=1,n_i
                    call mat_free(xsolvecm(j))
                    call mat_free(sigmasm(j))
                    call mat_free(rhosp(j))
                 enddo
                 do j=1,max(n_v,n_i)
                    call mat_free(xsolvec(j))
                    call mat_free(sigmas(j))
                    call mat_free(rhos(j))
                 enddo
                 call mem_dealloc(sigmas)
                 call mem_dealloc(rhos)
                 do j=1,ngd
                    call mat_free(xp(j))
                    call mat_free(xm(j))
                 enddo
                    exit
             endif    
          endif
      endif
      endif

      do j=1,n_v
           call mat_free(xsolvecp(j))
           call mat_free(sigmasp(j))
           call mat_free(rhosm(j))
      enddo
        do j=1,n_i
           call mat_free(xsolvecm(j))
           call mat_free(sigmasm(j))
           call mat_free(rhosp(j))
        enddo
         do j=1,max(n_v,n_i)
           call mat_free(xsolvec(j))
           call mat_free(sigmas(j))
           call mat_free(rhos(j))
         enddo
         call mem_dealloc(sigmas)
         call mem_dealloc(rhos)
    enddo
if (lineq) then
  do j=1,ngd 
      call mat_free(RHS_real(j))
      call mat_free(RHS_img(j))
    enddo
endif
    deallocate(E1,E2,S1,S2)
    deallocate(red_Xp_glob,red_Xm_glob)
    if (lineq) then
        deallocate(G1,G2)
    endif    
    
    CALL lsCLOSE(lu_sigma_rsp,'DELETE')
    CALL lsCLOSE(lu_rho_rsp,'DELETE')
    CALL lsCLOSE(lu_x_rsp,'DELETE')
end subroutine rsp_sym_solver
!================================================================  
!================================================================  
subroutine extend_new_reduced_matrices(molcfg,n_red,nm_red,nt_red,nx_new,nm_new,ngd,xsolvecp,&
                                       & sigmasp,rhosm,xsolvecm,sigmasm,rhosp,lineq,gd,gdi)
!================================================================  
! Purpose: extend reduced subspace with new vectors
!          a whole reduced space has a form
!           E1  wS1    Xg      Gg
!       (           )(    )= (    )
!           wS2  E2    Xu      Gu
!       where
!          E1 = Xg E[2] Xg
!          E2 = Xu E[2] Xu
!          S1 = Xg S[2] Xu
!          S2 = Xu S[2] Xg
!          The reduced E2, S2, G matrices are global variables
!          (see top of file) 
!
! INPUT: 
! n_red, number of g vectors on disk and current half size of reduced space
! nm_red, number of u vectors on disk and current half size of reduced space
! nt_red, number of vectors stored on disk (vectors are stored as g+u, 
!         in case when the trial vector is either symmetric or antisymmetric
!          the total number of trial vectors may be different than number of g/u)
! nx_NEW, number of new g vectors
! nm_NEW, number of new u vectors
! ngd, the number of response vectors to be found
! gd(ngd), matrix array containing ngd g right hand sides (gradients)
! gdi(ngd), matrix array containing ngd u right hand sides (gradients)
! sigmasp,rhosp,xsolvecp:  matrix arrays of the new g vectors
! sigmasm,rhosm,xsolvecm:  matrix arrays of the new u vectors
!
!================================================================  
implicit none
type(rsp_molcfg), intent(inout)  :: molcfg
integer,intent(in)               :: n_red,nm_red,nt_red
integer, intent(in)              :: nx_new,nm_new,ngd
logical                          :: lineq
type(Matrix),intent(in)          :: xsolvecp(:),sigmasp(:),rhosm(:)
type(Matrix),intent(in)          :: xsolvecm(:),sigmasm(:),rhosp(:)
type(Matrix),intent(in),optional :: gd(:),gdi(:)
integer                          :: ndim,k,j,i,l,max_red
type(Matrix)                     ::x_j,rho_j,sigma_j,gT
type(Matrix)                     ::xp_j,rhop_j,sigmap_j,scrT
type(Matrix)                     ::xm_j,sigmam_j
real(realk)                      :: a,b

 max_red= molcfg%solver%rsp_maxred
if ((N_red + Nx_new) > max_red) then
   WRITE (molcfg%LUPRI,*) 'ERROR in extend_complex_reduced_matrices: N_RED + N_NEW > rsp_maxred'
   WRITE (molcfg%LUPRI,*) 'N_RED ,Nx_NEW  =',N_red,Nx_NEW
   WRITE (molcfg%LUPRI,*) 'rsp_maxred =',max_red
   WRITE (molcfg%LUPRI,*) 'Reduce problem size or recompile orthonormalize '// &
        &                  'with larger rsp_maxred parameter'
   WRITE (*,*) 'ERROR in extend_complex_reduced_matrices: N_RED + N_NEW > rsp_maxred'
   WRITE (*,*) 'N_RED ,Nx_NEW  =',N_red,Nx_NEW
   WRITE (*,*) 'rsp_maxred =',max_red
   WRITE (*,*) 'Reduce problem size or recompile orthonormalize '// &
        &                  'with larger rsp_maxred parameter'
   CALL lsQUIT('sym1 orthonormalize error: N_RED + Nx_NEW > rsp_maxred parameter',molcfg%lupri)
END IF
if ((nm_red+nm_new) > max_red) then
   WRITE (molcfg%LUPRI,*) 'ERROR in extend_complex_reduced_matrices: Nm_red + Nm_NEW > max_red'
   WRITE (molcfg%LUPRI,*) 'Nm_RED ,Nm_NEW  =',Nm_red,Nm_NEW
   WRITE (molcfg%LUPRI,*) 'rsp_maxred =',max_red
   WRITE (molcfg%LUPRI,*) 'Reduce problem size or recompile orthonormalize '// &
        &                  'with larger rsp_maxred parameter'
   WRITE (*,*) 'ERROR in extend_complex_reduced_matrices: Nm_red + Nm_NEW > max_red'
   WRITE (*,*) 'Nm_RED ,Nm_NEW  =',Nm_red,Nm_NEW
   WRITE (*,*) 'rsp_maxred =',max_red
   WRITE (*,*) 'Reduce problem size or recompile orthonormalize '// &
        &                  'with larger rsp_maxred parameter'
   CALL lsQUIT('sym2 orthonormalize error: Nm_red + Nm_NEW > rsp_maxred parameter',molcfg%lupri)
END IF

ndim = xsolvecp(1)%nrow
call mat_init(sigma_j,ndim,ndim)
call mat_init(rho_j,ndim,ndim)
call mat_init(x_j,ndim,ndim)
call mat_init(sigmap_j,ndim,ndim)
call mat_init(xp_j,ndim,ndim)
call mat_init(sigmam_j,ndim,ndim)
call mat_init(rhop_j,ndim,ndim)
call mat_init(xm_j,ndim,ndim)
call mat_init(scrT,ndim,ndim)

if (nx_new>0) then
    k = 0
    do i = n_red+1, n_red+nx_new 
       k = k + 1
       l = 0
       do j = n_red+1, n_red+nx_new
          l = l + 1
          E1(i,j)     = mat_dotproduct(xsolvecp(k),sigmasp(l))
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
          E2(i,j)     = mat_dotproduct(xsolvecm(k),sigmasm(l))
       enddo      
    enddo
endif
if ((nx_new>0) .and. (nm_new>0)) then
    k = 0
    do i = n_red+1, n_red+nx_new 
       k = k + 1
       l = 0
       do j = nm_red+1, nm_red+nm_new
          l = l + 1
          S1(i,j)     = mat_dotproduct(xsolvecp(k),rhosp(l))
       enddo      
    enddo
endif

!Setup lower half of E2 and S2:
 rewind(lu_sigma_rsp)
 do j = 1, nt_red
    call mat_zero(sigma_j)
    call mat_read_from_disk(lu_sigma_rsp,sigma_j,RSPonMaster)
    call mat_trans(sigma_j,scrT)
    call mat_add(0.5E0_realk,sigma_j,0.5E0_realk,scrT,sigmap_j)
    call mat_add(0.5E0_realk,sigma_j,-0.5E0_realk,scrT,sigmam_j)
       
    if ((nx_new>0) .and. (sqrt(mat_sqnorm2(sigmap_j))>1E-9_realk)) then
         k = 0
         do i = n_red+1, n_red+nx_new
            k = k + 1
            E1(i,j)     = mat_dotproduct(xsolvecp(k),sigmap_j)
         enddo
     endif
     if ((nm_new>0) .and. (sqrt(mat_sqnorm2(sigmam_j))>1E-9_realk)) then
         k = 0
         do i = nm_red+1, nm_red+nm_new
            k = k + 1
            E2(i,j)    = mat_dotproduct(xsolvecm(k),sigmam_j)
         enddo
     endif
 enddo
 
 rewind(lu_rho_rsp)
 call mat_zero(rho_j)
 do j = 1, nt_red
    call mat_read_from_disk(lu_rho_rsp,rho_j,RSPonMaster)
    call mat_trans(rho_j,scrT)
    call mat_add(0.5E0_realk,rho_j,0.5E0_realk,scrT,rhop_j)
    if ((nx_new>0) .and. (sqrt(mat_sqnorm2(rhop_j))>1E-9_realk)) then
        k = 0
        do i = n_red+1, n_red+nx_new
           k = k + 1
           S1(i,j)=mat_dotproduct(xsolvecp(k),rhop_j)
        enddo
     endif
 enddo
 
 rewind(lu_x_rsp)
 !Explicitly calculate upper half of E2 and S2:
 do i = 1, nt_red
     call mat_read_from_disk(lu_x_rsp,x_j,RSPonMaster)
     call mat_trans(x_j,scrT)
     call mat_add(0.5E0_realk,x_j,0.5E0_realk,scrT,xp_j)
     call mat_add(0.5E0_realk,x_j,-0.5E0_realk,scrT,xm_j)
    
    if ((nx_new>0) .and. (sqrt(mat_sqnorm2(xp_j))>1E-9_realk))then
       k = 0
       do j = n_red+1, n_red+nx_new  !Run over 2 elements at a time because of pairing
          k = k + 1
          E1(i,j)   = mat_dotproduct(xp_j,sigmasp(k))
       enddo
     endif
    if ((nm_new>0) .and. (sqrt(mat_sqnorm2(xm_j))>1E-9_realk)) then
       k = 0
       do j = nm_red+1, nm_red+nm_new  !Run over 2 elements at a time because of pairing
          k = k + 1
          E2(i,j)   = mat_dotproduct(xm_j,sigmasm(k))
       enddo
    endif
     
    if ((nm_new>0) .and. (sqrt(mat_sqnorm2(xp_j))>1E-9_realk)) then
       k = 0
       do j = nm_red+1, nm_red+nm_new
          k = k + 1
          S1(i,j)   = mat_dotproduct(xp_j,rhosp(k))
       enddo
    endif
 enddo

 if (lineq) then
    do j = 1, ngd  !Loop over number of gradients
       k=0
       do i = n_red+1, n_red+nx_new !Loop in steps of 2 because of pairing
          k = k+1
            if (sqrt(mat_sqnorm2(gd(j)))>1E-9_realk) then    
           G1(i,j) = mat_dotproduct(xsolvecp(k),gd(j))
            else
           G1(i,j)=0E0_realk
           endif
       enddo
       k=0
       do i = nm_red+1, nm_red+nm_new !Loop in steps of 2 because of pairing
          k = k+1
            if (sqrt(mat_sqnorm2(gdi(j)))>1E-9_realk) then    
           G2(i,j) = mat_dotproduct(xsolvecm(k),gdi(j))
            else
           G2(i,j)=0E0_realk
           endif
       enddo
    enddo
 endif

 if (molcfg%solver%info_rsp_redspace) then 
          write(molcfg%lupri,*) 'Reduced E1:'
          call LS_OUTPUT(E1, 1, n_red+nx_new, 1, n_red+nx_new, max_red, max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced E2:'
          call LS_OUTPUT(E2, 1, nm_red+nm_new, 1, nm_red+nm_new, max_red, max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced S1:'
          call LS_OUTPUT(S1, 1, n_red+nx_new, 1, nm_red+nm_new, max_red, max_red, 1, molcfg%lupri)
  endif

call mat_free(x_j)
call mat_free(rho_j)
call mat_free(sigma_j)
call mat_free(xp_j)
call mat_free(sigmap_j)
call mat_free(xm_j)
call mat_free(rhop_j)
call mat_free(sigmam_j)
call mat_free(scrT)

    end subroutine extend_new_reduced_matrices

!================================================================  
 !================================================================  
subroutine solve_std(molcfg,ndim,ndim_red,nm_red,ngd,freq,red_X,red_Xm)
!
!  Solve reduced linear response problem in subspace
!  Now finally we consider one gradient/solution at a time
!
!  INPUT:
!  ndim:     size of the full matrix
!  ndim_red: size of the g reduced space (E1)
!  nm_red:   size of the u reduced space (E2)
!  ngd:      number of right hand sides
!  freq:     the ngd frequencies
!
!  OUTPUT:
!  red_X:    the ngd reduced space real g solution vectors
!  red_Xm:   the ngd reduced space real u solution vectors
!
!------------------------------------------------------------------
    implicit none
    type(rsp_molcfg), intent(inout) :: molcfg
    integer, intent(in)      :: ndim_red, ngd,ndim,nm_red
    real(realk),intent(in)   :: freq(:)
    real(realk),intent(inout):: red_X(:,:),red_Xm(:,:)
    real(realk),allocatable  :: S1p(:,:), S2p(:,:)
    real(realk),allocatable   :: A(:,:),KHS(:)
    integer,allocatable      :: IPIV(:)
    integer                  :: igd, ierr, i, j, k,max_red

    ierr=0

    max_red= molcfg%solver%rsp_maxred

    allocate(S2p(nm_red,ndim_red), S1p(ndim_red,nm_red))
    allocate(KHS(ndim_red+nm_red), IPIV(ndim_red+nm_red))
    allocate(A(ndim_red+nm_red,ndim_red+nm_red))
    S1p=0E0_realk; S2p=0E0_realk
    
    KHS=0E0_realk
    IPIV=0E0_realk
    A=0E0_realk

    !Setup reduced E2, S2, and right hand side with proper dimension
    S1p(1:ndim_red,1:nm_red)                                 = S1(1:ndim_red,1:nm_red)

    A(1:ndim_red,1:ndim_red)                                 = E1(1:ndim_red,1:ndim_red)
    A(ndim_red+1:ndim_red+nm_red,ndim_red+1:ndim_red+nm_red) = E2(1:nm_red,1:nm_red)
       
    do igd=1,ngd
       S1p = - freq(igd)*S1p
       do i=1,ndim_red
          do j=1,nm_red
             S2p(j,i)=S1p(i,j)
          enddo
        enddo 
       A(1:ndim_red,ndim_red+1:ndim_red+nm_red)   = S1p(1:ndim_red,1:nm_red)
       A(ndim_red+1:ndim_red+nm_red,1:ndim_red)   = S2p(1:nm_red,1:ndim_red)

       KHS(1:ndim_red)                            = G1(1:ndim_red,igd)
       KHS(ndim_red+1:ndim_red+nm_red)            = G2(1:nm_red,igd)
         
       if (molcfg%solver%info_rsp_redspace) then
          write(molcfg%lupri,*) 'Reduced A:'
          call LS_OUTPUT(A, 1, ndim_red+nm_red, 1, ndim_red+nm_red, 2*max_red, 2*max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced RHS:'
          call LS_OUTPUT(KHS, 1, ndim_red+nm_red, 1, 1, 2*max_red, 1, 1, molcfg%lupri)
       endif

       !Solve set of linear equations Ax = b:
       call DGESV(ndim_red+nm_red, 1, A, ndim_red+nm_red, IPIV, KHS, ndim_red+nm_red, IERR)
       !Solution vector is found in RHS.
       if (IERR /= 0) then
          WRITE(molcfg%LUPRI,'(/A, i4)') &
          &     'Problem in DGESV, IERR = ', IERR
          CALL lsQUIT(' Problem in DGESV',molcfg%lupri)
       endif
       
       red_X(1:ndim_red,igd) = KHS(1:ndim_red)
       red_Xm(1:nm_red,igd)  = KHS(ndim_red+1:ndim_red+nm_red)
    enddo

    deallocate(S1p, S2p)
    deallocate(KHS, A)
    deallocate(IPIV)
  end subroutine solve_std
  
  
  subroutine get_new_res(molcfg,F,D,S,itmic,ndim_red,nm_red,nt_red,ngd,red_Xp_glob,red_Xm_glob, &
             &    eival,xp,xm,nx_new,nm_new,conv,res_norm_tot,lineq,gd,gdi)
    implicit none
!******************************************************************
! Purpose: 
! 1)Construct residuals:
!   Rg=E(2)Xg-W(I)*S(2)Xu(I) - GDg
!   Ru=E(2)Xu-W(I)*S(2)Xg(I) - GDu
!   for NGD rsp vectors X(I) of reduced rsp problem
! 2)Test for convergence of ngd eigenvectors,
!   Convergence criterium:
!   ||Rg|| < rsp_thresh
!   ||Ru|| < rsp_thresh 
! 3) contruct new trial vectors as preconditioned residuals
!
!  INPUT:
!  itmic:    number of the microiteration (for printout)
!  ndim_red: half size of g reduced space (E1)
!  nm_red:   half size of u reduced space (E2)
!  nt_red:   number of vectors stored on disk (vectors are stored as g+u, 
!            in case when the trial vector is either symmetric or antisymmetric
!            the total number of trial vectors may be different than number of g/u)
!  ngd:      number of right hand sides
!  gd(ngd):  matrix array containing ngd g right hand sides (gradients)
!  gdi(ngd):  matrix array containing ngd u right hand sides (gradients)
!  red_Xp_glob:    the ngd reduced space g solution vectors (same as RED_X)
!  red_Xm_glob:    the ngd reduced space u solution vectors (same as RED_Xm)
!  eival:    the ngd frequencies
!  
!  OUTPUT:
!  nx_NEW:   number of new g vectors
!  nm_NEW:   number of new u vectors
!  xp(nx_new) new g trial vector
!  xm(nm_new) new u trial vector
!  conv:     set true if converged
! res_norm_tot : norm of residual to check if convergence not "stuck"
!******************************************************************

    type(rsp_molcfg), intent(inout)  :: molcfg
    logical,intent(in)               :: lineq
    integer,intent(in)               :: itmic,ngd,nt_red,ndim_red,nm_red
    type(Matrix),intent(in)          :: F,D,S
    type(Matrix),intent(in),optional :: gd(:),gdi(:)
    real(realk), intent(in)          :: red_Xp_glob(:,:)
    real(realk), intent(inout)          :: eival(:)
    real(realk), intent(in)          :: red_Xm_glob(:,:)
    type(Matrix),intent(inout)       :: xp(:),xm(:)
    logical, intent(out)             :: conv
    real(realk),intent(out)          :: res_norm_tot
    integer, intent(out)             :: nx_new,nm_new
!local
    type(Matrix)                     :: gdp,gdpi,gdm,gdmi,gT
    type(Matrix)                     :: residualsp(ngd),residualsm(ngd),rhsp(ngd),rhsm(ngd)
    type(Matrix)                     :: residualspi(ngd),residualsmi(ngd),rhspi(ngd),rhsmi(ngd)
    real(realk)                      :: wibx,res_norm,x_red_norm,conv_test,ddot,diff,x_red_norm_img
    integer                          :: i,n_not_conv,ndim_red_mat,ndim,index_conv(ngd),ibx,j,k
    real(realk)                      :: a, b,c
    real(realk)                      :: av_norm,res_norm_img, conv_test_img
    logical                         :: conv_root, conv_root_img,all_conv
    !debug
    type(matrix)                     :: r,xa(ngd),xb(ngd),Lxa(ngd),Lxb(ngd)
    type(matrix)                     :: S2Xp(ngd),S2Xm(ngd),scr
    real(realk)                      :: omega1(ngd),rsp_thresh,rsp_conv_thr


       rsp_thresh= molcfg%solver%rsp_thresh
    ndim = S%nrow
    n_not_conv = 0
    conv=.false.
    conv_root=.false.
    conv_root_img=.false.
    omega1=0E0_realk
    all_conv=.false.
   
! find E[2]*X(i)-w(i) S[2]*X(i) as linear combination of 
! previous Sigmas and RHOs. 
! X(i) = sum_j Bvecold_j REDX_j(i)
! Sigmanew(i) = E[2]X(i) = sum_j REDX_j(i) E[2]Bvecold_j =
!             = sum_j REDX_j(i) Sigmaold_j
! and add to previous -w_i S[2]*X(I) vector
! 
call mat_init(scr,ndim,ndim)
call mat_zero(scr)
    do ibx = 1,ngd
      !get the frequency
       wibx = eival(ibx)
       call mat_init(residualsp(ibx),ndim,ndim)
       call mat_init(residualsm(ibx),ndim,ndim)
       call mat_zero(residualsp(ibx))
       call mat_zero(residualsm(ibx))
      if (molcfg%solver%rsp_olsen) then
          call mat_init(S2Xp(ibx),ndim,ndim)
          call mat_init(S2Xm(ibx),ndim,ndim)
          call mat_init(xa(ibx),ndim,ndim)
          call mat_init(xb(ibx),ndim,ndim)
          call mat_init(Lxa(ibx),ndim,ndim)
          call mat_init(Lxb(ibx),ndim,ndim)
       endif       

       rewind(lu_rho_rsp)
       rewind(lu_sigma_rsp)
       call expand_on_basis4(molcfg,nt_red,ndim,red_Xm_glob(1:nm_red,ibx),red_Xp_glob(1:ndim_red,ibx),&
                            &lu_rho_rsp,residualsp(ibx),residualsm(ibx))
       
      if (molcfg%solver%rsp_olsen) then
       call mat_assign(S2Xp(ibx),residualsm(ibx))
       call mat_assign(S2Xm(ibx),residualsp(ibx))
       endif
       call mat_scal(-wibx,residualsp(ibx))
       call mat_scal(-wibx,residualsm(ibx))
       !add  sum_j ^RX_j(ibx) sigma_j
       call expand_on_basis4(molcfg,nt_red,ndim,red_Xp_glob(1:ndim_red,ibx),red_Xm_glob(1:nm_red,ibx),&
                            &lu_sigma_rsp,residualsp(ibx),residualsm(ibx))
      
       ! subtract gradient
      if (lineq) then
         call mat_daxpy(-1.0E0_realk,gd(ibx),residualsp(ibx))
         call mat_daxpy(-1E0_realk,gdi(ibx),residualsm(ibx))
      endif
      call util_scriptPx('T',D,S,residualsp(ibx))
      call util_scriptPx('T',D,S,residualsm(ibx))
    enddo

 call mat_add(1E0_realk,residualsp(1),1E0_realk,residualsm(1),scr)
  write(molcfg%lupri,*) 'residual tot norm in it',itmic, sqrt(mat_sqnorm2(scr))
! the residual(s) is now done: test for convergence. If not converged
! form new trial(s)  in next subroutine
! New trial(s) is equal to preconditioned residual(s)
   if (molcfg%solver%rsp_convdyn) then !default
      if (ITMIC == 1) then
          av_norm = 0.0E0_realk
          molcfg%solver%rsp_dyn_thresh=0E0_realk
          do i = 1, ngd
             res_norm     = sqrt(mat_sqnorm2(residualsp(i)))
             res_norm_img = sqrt(mat_sqnorm2(residualsm(i)))
             av_norm = av_norm + res_norm + res_norm_img
          enddo
          av_norm = 0.5E0_realk*av_norm
          rsp_conv_thr = av_norm*molcfg%solver%rsp_conv_factor
          molcfg%solver%rsp_dyn_thresh=rsp_conv_thr
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
         res_norm=mat_sqnorm2(residualsp(ibx))
         res_norm_img=mat_sqnorm2(residualsm(ibx))
      else ! double norm
         res_norm    =sqrt(mat_sqnorm2(residualsp(ibx)))
         res_norm_img=sqrt(mat_sqnorm2(residualsm(ibx)))
      endif
      if (res_norm < rsp_thresh) conv_root = .true.
      if (res_norm_img < rsp_thresh) conv_root_img = .true.
      res_norm_tot=res_norm+res_norm_img

     ! Insert back print statements
      write (molcfg%LUPRI,'(A,I3,A,E14.6,A,F12.6,A,I3,A,L2)') &
           &'Residual norm for herm vector ',ibx,' is: ',res_norm,&
           &'    and frequency = ',EIVAL(ibx),'   It = ', itmic,&
           &' CONV =',conv_root
    
      write (molcfg%LUPRI,'(A,I3,A,E14.6,A,F12.6,A,I3,A,L2)') &
           &'Residual norm for antiherm vector ',ibx,' is: ',res_norm_img,&
           &'    and frequency = ',EIVAL(ibx),'   It = ', itmic,&
           &' CONV =',conv_root_img
         
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
          !this root is converged
          index_conv(ibx) = -1
       else
          !this root is NOT converged
          n_not_conv = n_not_conv + 1
          index_conv(ibx) = 0
       endif
    enddo
 
   if ((ngd==1) .and. (conv_root .and. conv_root_img)) then 
       conv=.true.
   else
      if (molcfg%solver%rsp_olsen) then
       write(molcfg%lupri,*) 'olsen algorithm is used'
       do ibx=1,ngd
        call expand_on_basis4(molcfg,nt_red,ndim,red_Xp_glob(1:ndim_red,ibx),red_Xm_glob(1:ndim_red,ibx),lu_x_rsp,&
                             &xa(ibx),xb(ibx))
        call new_std_MO_precond(xa(ibx),xb(ibx),S,eival(ibx),Lxa(ibx),Lxb(ibx),molcfg%decomp%nocc)
        omega1(ibx)=(mat_dotproduct(Lxa(ibx),residualsp(ibx))+mat_dotproduct(Lxb(ibx),residualsm(ibx)))/&
                   &(mat_dotproduct(Lxa(ibx),S2Xm(ibx))+mat_dotproduct(Lxb(ibx),S2Xp(ibx)))
        call mat_daxpy(-omega1(ibx),S2Xm(ibx),residualsp(ibx))
        call mat_daxpy(-omega1(ibx),S2Xp(ibx),residualsm(ibx))
       enddo
      endif
    !Remove converged vectors
      nx_new=0
      nm_new=0
      all_conv = .true.
      do ibx = 1,ngd
         if (index_conv(ibx) == 0) then !vec not converged
             nm_new=nm_new+1
             nx_new=nx_new+1
             if (nx_new /= ibx) residualsp(nx_new) = residualsp(ibx)
             if (nm_new /= ibx) residualsm(nm_new) = residualsm(ibx)
             call new_std_MO_precond(residualsp(nx_new),residualsm(nm_new),S,eival(ibx),xp(nx_new),xm(nm_new),molcfg%decomp%nocc)
             call util_scriptPx('N',D,S,xp(nx_new))
             call util_scriptPx('N',D,S,xm(nm_new))
             all_conv = .false.
         endif
     enddo
  endif
         
if (all_conv) conv=.true.

! Finished, deallocate local arrays
   do ibx=1,ngd
        call mat_free(residualsp(ibx))
        call mat_free(residualsm(ibx))
      if (molcfg%solver%rsp_olsen) then
          call mat_free(S2Xp(ibx))
          call mat_free(S2Xm(ibx))
          call mat_free(xa(ibx))
          call mat_free(xb(ibx))
          call mat_free(Lxa(ibx))
          call mat_free(Lxb(ibx))
        endif
   enddo

call mat_free(scr)


  end subroutine get_new_res
!================================================================  

  subroutine get_start_vectors(molcfg,F,D,S,RHS_real,RHS_img,ngd,eival,nx_new,nm_new,xp,xm)
    implicit none
!================================================================ 
! A subroutine to get start vectors 
! by preconditioning the RHS vectors.
! g and u components obtianed directly.
!================================================================ 
    type(rsp_molcfg), intent(inout) :: molcfg
    type(Matrix),intent(in)    :: S,D,F
    type(Matrix),intent(inout) :: RHS_real(:),RHS_img(:)
    type(Matrix),intent(inout) :: xp(:),xm(:)
    real(realk), intent(inout) :: eival(:) 
    integer,intent(inout)      :: nx_new,nm_new
    integer,intent(in)         :: ngd
    integer                    :: i,ndim

ndim=S%nrow
do i=1,ngd
          nx_new=nx_new+1
          nm_new=nm_new+1
         call new_std_MO_precond(RHS_real(i),RHS_img(i),S,eival(i),xp(nx_new),xm(nm_new),molcfg%decomp%nocc)
         call util_scriptPx('N',D,S,xp(nx_new))
         call util_scriptPx('N',D,S,xm(nm_new))
enddo
end subroutine get_start_vectors

  subroutine orthonormalize22(molcfg,D,S,Bvec_tmp,Nb_new,Nb_prev,bvecs,lu_basis)
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
    integer                    :: irx,i,iturn,ndim,k,ibvec,jbvec,jrx,no_of_new,lub_rsp_vec,dummy_i
    integer,allocatable        :: lin_depend(:) !linear dependency index array
    type(matrix)               :: B_scr, b_k, Xf, rho_k
    type(matrix)               :: orthovec !Will properly be changed
    real(realk)                :: TT,T1,T2,dummy_real
    logical                    :: run_ortho_again,OnMaster
    integer,intent(in)                :: lu_basis
    integer                    ::max_red

    OnMaster = .TRUE.
 max_red= molcfg%solver%rsp_maxred

    if ((Nb_prev + Nb_new) > max_red) then
      WRITE (molcfg%LUPRI,*) 'ERROR in orthonormalize: NB_PREV + NB_NEW > rsp_maxred'
      WRITE (molcfg%LUPRI,*) 'NB_PREV ,NB_NEW  =',NB_PREV,NB_NEW
      WRITE (molcfg%LUPRI,*) 'rsp_maxred =',max_red
      WRITE (molcfg%LUPRI,*) 'Reduce problem size or recompile orthonormalize '// &
      &                  'with larger rsp_maxred parameter'
      CALL lsQUIT('orthonormalize error: NB_PREV + NB_NEW > rsp_maxred parameter',molcfg%lupri)
    END IF
    
    ndim = S%nrow
    ALLOCATE(lin_depend(Nb_new))
    call mat_init(B_scr,ndim,ndim)
    call mat_init(b_k,ndim,ndim)   

    iturn=0
    ! STEP 1: check initial linear dependencies among new vectors
    if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'orthonormalize -> remove_initial_lindep'
    lin_depend = 1 !Initialize to non-dependent
    do   !It might be necessary to run through everything twice - decided in normalize
       iturn = iturn + 1
   !
   !Project out
   ! Orthogonalize new b-vectors agains previous (old) b-vectors
   ! (both (Z, Y) and (Y, Z))
   !
       if (molcfg%solver%info_rsp) then
         write(molcfg%lupri,*) 'Orthogonalize new b-vectors agains previous b'
         write(molcfg%lupri,*) 'NB_PREV, NB_NEW ', NB_PREV, NB_NEW
       endif 
   
       rewind(lu_basis)
       do k = 1,Nb_prev
                  call mat_read_from_disk(lu_basis,b_k,OnMaster)
          do irx = 1,Nb_new
             !if lin_depend(irx) == 0, the vector is skipped
             !because of linear dependencies
             if (lin_depend(irx) /= 0) then
                 !Orthogonalize to (Z Y)_old
                 TT = mat_dotproduct(b_k,Bvec_tmp(irx))
                 call mat_daxpy(-TT,b_k,Bvec_tmp(irx))
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
       enddo !ibvec
   
       if (run_ortho_again) then
          run_ortho_again = .false.
          if (iturn == 2) then !This is redundant since run_ortho_again is only set .true. if iturn = 1
             WRITE(molcfg%LUPRI,'(/A)') &
             &     'Error: Already ran twice through orthonormalize22!'
             CALL lsQUIT('Error: Already ran twice through orthonormalize!22',molcfg%lupri)
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

  subroutine expand_on_basis3(molcfg,ndim_red,ndim,red_eivec_i,lu_basis,expanded_vec)
!********************************************************************
! PURPOSE:
!   Expand reduced space solution on full basis 
!
!   INPUT:
!   ndim:     size of matrices/number of basis functions
!   ndim_red: number of vectors on disk
!   lu_basis: the lun for the file where the basis vectors are found
!   red_EIVEC_i:    contains the expansion coefficients 
!
!   OUTPUT:
!   expanded_vec: vector expanded in the set of basis vectors
!   
!    
!******************************************************
    !expanded_vec must be zeroed outside first time it is used!!!
    implicit none
    type(rsp_molcfg), intent(inout) :: molcfg
    integer, intent(in) :: ndim_red, lu_basis, ndim
    real(realk), intent(in) :: red_eivec_i(:)
    type(Matrix), intent(inout) :: expanded_vec  

    type(matrix) :: basisvec_k
    real(realk) :: fac1,bnorm,dprodu
    integer :: k
    logical :: OnMaster
    OnMaster = .TRUE.
    call mat_init(basisvec_k,ndim,ndim)

    rewind(lu_basis)
    do k = 1,ndim_red !length of red. solution for each grad
                      !and # of bf's (no pair)
      call mat_read_from_disk(lu_basis,basisvec_k,OnMaster) 
     
      fac1 = red_eivec_i(k)
      bnorm = mat_sqnorm2(Basisvec_k)
      if (molcfg%solver%info_rsp) then
        write(molcfg%lupri,*) 'expand_on_basis2: bnorm^2 = ', bnorm
        !write(molcfg%lupri,*) 'Vettore di base nro: ', k
        !call mat_print(Basis(k),1,ndim,1,ndim,molcfg%lupri)
      endif
      call mat_daxpy(fac1,Basisvec_k,expanded_vec)
      !if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'After summation, expanded_vec, k = ',k
      !call mat_print(expanded_vec,1,ndim,1,ndim,molcfg%lupri)
    enddo
    call mat_free(basisvec_k)

  end subroutine expand_on_basis3

  subroutine get_1st_trials(molcfg,F,D,S,ngd,eival,nx_new,nm_new,xp,xm,lineq,gdp,gdm)
!--------------------------------------
! Purpose: set up initial trial vectors 
!--------------------------------------
    implicit none
    type(rsp_molcfg), intent(inout)      :: molcfg
    logical                              :: lineq
    type(Matrix), intent(in)             :: F,D,S
    type(Matrix), intent(inout),optional :: gdp(rsp_number_of_rhs),gdm(rsp_number_of_rhs)
    real(realk), intent(inout)           :: eival(rsp_number_of_omegas)
    integer, intent(in)                  :: ngd
    type(Matrix),intent(inout)           :: xp(rsp_bvec_dim),xm(rsp_bvec_dim)  !output   
    integer, intent(out)                 :: nx_new,nm_new

    IF (LINEQ) THEN              
       if (molcfg%solver%info_rsp) WRITE(molcfg%LUPRI,*) 'Now enter 1st guess for linear equations '
   call get_start_vectors(molcfg,F,D,S,gdp,gdm,ngd,eival,nx_new,nm_new,xp,xm)
    ELSE
       if (molcfg%solver%info_rsp) WRITE(molcfg%LUPRI,*) 'Now enter 1st guess for excitation energies'
       call get_1st_orth_trial_eigen(molcfg,D,S,ngd,xp,xm,nx_new,nm_new)
    END IF

  end subroutine get_1st_trials
  subroutine get_1st_orth_trial_eigen(molcfg,D,S,nexcit,bvecs,bvecsm,Nb_new,nm_new)
!*******************************************************
! Purpose (~OLD PPST routine):
! generate start vector(s) for solution of generalized eigenvalue equations via:
! MO ->backtranforming the MO unit vectors corresponding to the resorted
!            roots E[2]_ii/|S[2]_ii| into AO
! CHO -> backtranforming the Cholesky unit vectors corresponding to the
!        resorted A_ab,ab diagonal elements into AO
!        and solving an eigenvalue problem in the Cholesky space
!        to generate orbital energy differences  
! ÅLea and Sonia, Feb. 2005/Sonia, Nov. 2005
!*******************************************************
    implicit none
    type(rsp_molcfg), intent(inout) :: molcfg
    integer,intent(in)         :: nexcit     !number of roots (required in input)
    type(Matrix),intent(in)    :: D,S
    type(Matrix),intent(inout) :: bvecs(rsp_bvec_dim),bvecsm(rsp_bvec_dim)  !output   
    integer, intent(out)       :: Nb_new,nm_new
    type(Matrix),allocatable   :: Bvec_tmp(:) !it is possible to have more start vecs than nexcit
    type(matrix)               :: bT
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

    if (howmany > rsp_bvec_dim) then
       WRITE(molcfg%LUPRI,'(/A)') &
       &     'More start vectors than b vectors allowed in core)'
       CALL lsQUIT('More start vectors than b vectors allowed in core',molcfg%lupri)
    endif

    allocate(Bvec_tmp(howmany))
    do i = 1,howmany
       call mat_init(Bvec_tmp(i),ndim,ndim)
    enddo
     call mat_init(bT,ndim,ndim)
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
           CALL lsQUIT('File "rsp_eigenvecs" does not exist!',molcfg%lupri)
       endif
       lurestart = -1
       CALL LSOPEN(lurestart,'rsp_eigenvecs','old','UNFORMATTED')
       do i = 1, molcfg%solver%rsp_restart_nexci
          read(lurestart) k, eival
          call mat_read_from_disk(lurestart,Bvec_tmp(i),OnMaster)          
       enddo
       call lsCLOSE(lurestart,'DELETE') !Delete so it doesn't conflict with creation of
    endif                               !new "rsp_eigenvecs" after the calculation is done
!
! orthogonalize and remove linear tors
! Nb_new is updated inside orthonormalize
!
    Nb_prev = 0
    Nb_new = howmany
    if (Nb_new == 0) then
       WRITE (molcfg%LUPRI,'(//A,I5)') &
                    &  'get_1st_orth_trial_lineq: No start vectors present!'
       CALL lsQUIT('get_1st_orth_trial_lineq: No start vectors present!',molcfg%lupri)
    endif

! Number of start trial vectors = Nb_new  
!
    nm_new=nb_new
    do i=1,nb_new
       call mat_trans(Bvec_tmp(i),bT)
       call mat_add(0.5E0_realk,Bvec_tmp(i),0.5E0_realk,bT,bvecs(i))
       call mat_add(0.5E0_realk,Bvec_tmp(i),-0.5E0_realk,bT,bvecsm(i))
    enddo


    if (Nb_new <= 0) then
       WRITE (molcfg%LUPRI,'(//A,I5)') &
&       'get_1st_orth_trial_eigen: START VECTOR IS NOT LINEAR INDEPENDENT.'
      CALL lsQUIT('START VECTOR IS NOT LINEAR INDEPENDENT',molcfg%lupri)
    endif
    if (Nb_new < howmany) then
       if(molcfg%solver%info_rsp)write(molcfg%lupri,*) 'Inside get_1st_orth_trial_eigen'
       if(molcfg%solver%info_rsp)write(molcfg%lupri,*) '**WARNING: Less start vectors than required excitations!'
    end if
!
! Free memory
!
call mat_free(bT)
    do i = 1,howmany
       call mat_free(Bvec_tmp(i))
    enddo
    deallocate(Bvec_tmp)

  end subroutine get_1st_orth_trial_eigen
subroutine solve_red_eigen(molcfg,ndim,ndim_red,nm_red,nexci,eival,red_X,red_Xm)
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
    integer, intent(in)      :: ndim_red, nexci,ndim,nm_red
    real(realk),intent(inout)   :: eival(:)
    real(realk),intent(inout):: red_X(:,:),red_Xm(:,:)
    real(realk),allocatable  :: S1p(:,:), S2p(:,:)
    real(realk),allocatable   :: A(:,:),S2FULL(:,:)
    real(realk), allocatable :: alfi(:),beta(:),alfr(:)
    real(realk), allocatable :: red_eivec_scr(:)
    integer                  :: igd, ierr, i, j, k,ndim_red_mat,l,isndx(3),max_red
    real(realk), allocatable :: scr1(:,:),scr2(:,:)
    real(realk), parameter :: drtest = 1.0E-5_realk

    ndim_red_mat = ndim_red+nm_red

 max_red= molcfg%solver%rsp_maxred

    allocate(S2p(nm_red,ndim_red), S1p(ndim_red,nm_red))
    allocate(A(ndim_red_mat,ndim_red_mat))
    allocate(S2FULL(ndim_red_mat,ndim_red_mat))
    ALLOCATE(alfi(ndim_red_mat),alfr(ndim_red_mat),beta(ndim_red_mat))
    ALLOCATE(red_eivec_scr(ndim_red_mat*ndim_red_mat)) !FIXME: check this dimension...
    allocate(scr1(ndim_red_mat,ndim_red_mat),scr2(ndim_red_mat,ndim_red_mat))
    S1p=0E0_realk; S2p=0E0_realk
    S2FULL=0E0_realk 
    A=0E0_realk

    !Setup reduced E2, S2, and right hand side with proper dimension
    S1p(1:ndim_red,1:nm_red) = S1(1:ndim_red,1:nm_red)
    do i=1,ndim_red
       do j=1,nm_red
           S2p(j,i)=S1p(i,j)
       enddo
    enddo 

       A(1:ndim_red,1:ndim_red)                       = E1(1:ndim_red,1:ndim_red)
       A(ndim_red+1:2*ndim_red,ndim_red+1:2*ndim_red) = E2(1:ndim_red,1:ndim_red)
       
       S2FULL(1:ndim_red,ndim_red+1:ndim_red_mat)       = S1p(1:ndim_red,1:nm_red)
       S2FULL(ndim_red+1:ndim_red_mat,1:ndim_red)       = S2p(1:nm_red,1:ndim_red)

      if (molcfg%solver%info_rsp_redspace) then
          write(molcfg%lupri,*) 'Reduced E:'
          call LS_OUTPUT(A, 1, ndim_red_mat, 1, ndim_red_mat, 2*max_red, 2*max_red, 1, molcfg%lupri)
          write(molcfg%lupri,*) 'Reduced S:'
          call LS_OUTPUT(S2FULL, 1, ndim_red_mat, 1, ndim_red_mat, 2*max_red, 2*max_red, 1, molcfg%lupri)
      endif

!************************************************************
!        Use EISPACK routine for real general matrices in
!        generalized eigenvalue problem
!
!        CALL RGG(NM,N,A,B,ALFR,ALFI,BETA,MATZ,Z,IERR)
!***********************************************************
    call RGG(ndim_red_mat,ndim_red_mat,A,S2FULL,ALFR,ALFI,BETA,1,red_eivec_scr,IERR)

    IF ( IERR /= 0 ) THEN
       WRITE(molcfg%LUPRI,'(/A)') 'solve_red_eigen: Problem in RGG, IERR = ', IERR
       CALL lsQUIT('Problem in RGG',molcfg%lupri)
    END IF
!
! Reorder the reduced eigenvectors
!
    !S2 is regenerated because it is destroyed in RGG
       S2FULL = 0E0_realk
       S2FULL(1:ndim_red,ndim_red+1:ndim_red_mat)=S1p(1:ndim_red,1:nm_red)
       S2FULL(ndim_red+1:ndim_red_mat,1:ndim_red)=S2p(1:nm_red,1:ndim_red)

    if (molcfg%solver%info_rsp_redspace) then
        WRITE (molcfg%LUPRI,'(/A)') ' REDUCED EIGENVECTORS BEFORE ORDERING:'
        CALL LS_OUTPUT(RED_EIVEC_scr,1,ndim_red_mat,1,ndim_red_mat, &
    &                     ndim_red_mat,ndim_red_mat,1,molcfg%LUPRI)
    endif
    CALL ReorderEigenvalues2(molcfg,ndim_red_mat,S2FULL,RED_EIVEC_scr, & 
    &                       alfr,ALFI,BETA,scr1,scr2,ISNDX)
    if (molcfg%solver%info_rsp_redspace) then
        WRITE (molcfg%LUPRI,'(/A)') ' REDUCED EIGENVECTORS AFTER ORDERING:'
        CALL LS_OUTPUT(RED_EIVEC_scr,1,ndim_red_mat,1,ndim_red_mat, &
    &                     ndim_red_mat,ndim_red_mat,1,molcfg%LUPRI)
    endif

    if ((ndim_red_mat) .GT. nexci) then
       eival(1:nexci) = alfr(1:nexci) 

    elseif ((ndim_red_mat) .LE. nexci) then
       eival(1:ndim_red_mat) = alfr 
    endif

! Put residual results in our RED_EIVEC matrix
!
    if (ndim_red_mat > 2*max_red) then !Stinne fix 30/9-06
       WRITE (molcfg%LUPRI,*) 'rsp_maxred, ndim_red_mat:', max_red, ndim_red_mat
       WRITE (molcfg%LUPRI,*) 'rsp_maxred too small (must be >= 1/2*ndim_red_mat)'
       CALL lsQUIT('rsp_maxred too small',molcfg%lupri)
    endif
    if (nexci .LE. ndim_red_mat) then
       do i=1,ndim_red
          do j=1,nexci
             red_X(i,j) = red_eivec_scr((j-1)*ndim_red_mat+i)
          enddo
       enddo
       l=0
       do j=1,nexci
          l=l+1
       k=0
          do i=ndim_red+1,ndim_red_mat
             k=k+1
             red_Xm(k,l) = red_eivec_scr((j-1)*ndim_red_mat+i)
          enddo
       enddo
    else
       do i=1,ndim_red
          do j=1,ndim_red_mat
             red_X(i,j) = red_eivec_scr((j-1)*ndim_red_mat+i)
          enddo
       enddo
       l=0
       do j=1,ndim_red_mat
          l=l+1
       k=0
          do i=ndim_red+1,ndim_red_mat
              k=k+1
             red_Xm(k,l) = red_eivec_scr((j-1)*ndim_red_mat+i)
          enddo
       enddo
    endif
    DEALLOCATE(red_eivec_scr)
    DEALLOCATE(scr1,scr2)
    DEALLOCATE(A)
    deallocate(S2FULL)
    DEALLOCATE(alfi,beta,alfr)
  end subroutine solve_red_eigen
      
  subroutine expand_on_basis4(molcfg,nt_red,ndim,red_eivec_i,red_eivec_m,lu_basis,expanded_vec,expanded_vecm)
!********************************************************************
! PURPOSE:
!   Expand reduced space solution on full basis for the g-u solver
! vectors stored on disk in a form g+u
!
!   INPUT:
!   ndim:     size of matrices/number of basis functions
!   ndim_red: number of vectors on disk
!   lu_basis: the lun for the file where the basis vectors are found
!   red_EIVEC_i:    contains the expansion coefficients 
!
!   OUTPUT:
!   expanded_vec: vector expanded in the set of basis vectors
!   
!    
!******************************************************
    !expanded_vec must be zeroed outside first time it is used!!!
    implicit none
    type(rsp_molcfg), intent(inout) :: molcfg
    integer, intent(in) :: nt_red, lu_basis, ndim
    real(realk), intent(in) :: red_eivec_i(:),red_eivec_m(:)
    type(Matrix), intent(inout) :: expanded_vec,expanded_vecm  

    type(matrix) :: basisvec_k,bT,bp,bm
    real(realk) :: fac1,fac2,bnorm,dprodu
    integer :: k,i,j
    logical :: OnMaster
    OnMaster = .TRUE.
    call mat_init(basisvec_k,ndim,ndim)
    call mat_init(bT,ndim,ndim)
    call mat_init(bp,ndim,ndim)
    call mat_init(bm,ndim,ndim)

    i=0
    j=0
    rewind(lu_basis)
    do k = 1,nt_red !length of red. solution for each grad
                      !and # of bf's (no pair)
      call mat_read_from_disk(lu_basis,basisvec_k,OnMaster)
      call mat_trans(basisvec_k,bT)
     call mat_add(0.5E0_realk,basisvec_k,0.5E0_realk,bT,bp) 
     call mat_add(0.5E0_realk,basisvec_k,-0.5E0_realk,bT,bm) 
    
     if (sqrt(mat_sqnorm2(bp))>1E-9_realk) then
      i=i+1
      fac1 = red_eivec_i(i)
      call mat_daxpy(fac1,bp,expanded_vec)
      endif
     if (sqrt(mat_sqnorm2(bm))>1E-9_realk) then
      j=j+1
      fac2 = red_eivec_m(j)
      call mat_daxpy(fac2,bm,expanded_vecm)
      endif
    enddo
    call mat_free(basisvec_k)
    call mat_free(bp)
    call mat_free(bm)
    call mat_free(bT)

  end subroutine expand_on_basis4
  
  subroutine orthonormalize21(molcfg,D,S,Bvec_tmp,bm_tmp,Nb_new,nm_new,Nb_prev,bvecs,bvecsm,lu_basis,lineq)
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
    integer, intent(inout)     :: Nb_new,nm_new
    type(Matrix),intent(inout) :: Bvec_tmp(:), bvecs(:)
    type(Matrix),intent(inout) :: Bm_tmp(:), bvecsm(:)
    logical                    :: lineq
    !local 
    integer                    :: irx,i,iturn,ndim,k,ibvec,jbvec,jrx,no_of_new,lub_rsp_vec,dummy_i,max_red
    integer,allocatable        :: lin_depend(:) !linear dependency index array
    type(matrix)               :: B_scr, b_k, Xf, rho_k,bT,bp,bm
    type(matrix)               :: orthovec !Will properly be changed
    real(realk)                :: TT,T1,T2,dummy_real
    logical                    :: run_ortho_again,OnMaster
    integer,intent(in)                :: lu_basis
    OnMaster=.TRUE.
    max_red= molcfg%solver%rsp_maxred

    if ((Nb_prev + Nb_new) > max_red) then
       WRITE (molcfg%LUPRI,*) 'ERROR in orthonormalize: NB_PREV + NB_NEW > rsp_maxred'
       WRITE (molcfg%LUPRI,*) 'NB_PREV ,NB_NEW  =',NB_PREV,NB_NEW
       WRITE (molcfg%LUPRI,*) 'rsp_maxred =',max_red
       WRITE (molcfg%LUPRI,*) 'Reduce problem size or recompile orthonormalize '// &
            &                  'with larger rsp_maxred parameter'
       CALL lsQUIT('orthonormalize error: NB_PREV + NB_NEW > rsp_maxred parameter',molcfg%lupri)
    END IF

    ndim = S%nrow
    ALLOCATE(lin_depend(Nb_new))
    call mat_init(B_scr,ndim,ndim)
    call mat_init(b_k,ndim,ndim)   
    call mat_init(bT,ndim,ndim)   
    call mat_init(bp,ndim,ndim)   
    call mat_init(bm,ndim,ndim)   
    iturn=0

    ! STEP 1: check initial linear dependencies among new vectors
    if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'orthonormalize -> remove_initial_lindep'
    !  lin_depend = 1 !Initialize to non-dependent
    do   !It might be necessary to run through everything twice - decided in normalize
       iturn = iturn + 1
       !
       !Project out
       ! Orthogonalize new b-vectors agains previous (old) b-vectors
       ! (both (Z, Y) and (Y, Z))
       !
       if (molcfg%solver%info_rsp) then
          write(molcfg%lupri,*) 'Orthogonalize new b-vectors agains previous b'
          write(molcfg%lupri,*) 'NB_PREV, NB_NEW ', NB_PREV, NB_NEW
       endif

       rewind(lu_basis)
       do k = 1,Nb_prev
          call mat_read_from_disk(lu_basis,b_k,OnMaster)
          call mat_trans(b_k,bT)
          call mat_add(0.5E0_realk,b_k,0.5E0_realk,bT,bp)
          call mat_add(0.5E0_realk,b_k,-0.5E0_realk,bT,bm)
          if (sqrt(mat_sqnorm2(bp)) >1E-9_realk) then
             do irx = 1,Nb_new
                !if lin_depend(irx) == 0, the vector is skipped
                !because of linear dependencies
                !          if (lin_depend(irx) /= 0) then
                !Orthogonalize to (Z Y)_old
                TT = mat_dotproduct(bp,Bvec_tmp(irx))
                call mat_daxpy(-TT,bp,Bvec_tmp(irx))
                !         endif
             enddo
          endif
          if (sqrt(mat_sqnorm2(bm))>1E-9_realk) then
             do irx = 1,Nm_new
                !if lin_depend(irx) == 0, the vector is skipped
                !because of linear dependencies
                !        if (lin_depend(irx) /= 0) then
                !Orthogonalize to (Z Y)_old
                TT = mat_dotproduct(bm,Bm_tmp(irx))
                call mat_daxpy(-TT,bm,Bm_tmp(irx))
                !       endif
             enddo
          endif
       enddo


       if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Orthogonalize new b-vectors against each other '
       !
       ! Orthogonalize new vectors against each other
       !
       do ibvec = 1,Nb_new !index for current bvector 
          if (molcfg%solver%info_rsp) write(molcfg%lupri,*)'Current ibvec index is ', ibvec, lin_depend(ibvec)
          !   if (lin_depend(ibvec) /= 0) then
          jbvec = 1    !index for another bvector
          do jrx = 1,(ibvec-1)
             !      if (lin_depend(jrx) == 0) then
             !       jbvec = jbvec + 1 !skip
             !     else
             T1 = mat_sqnorm2(Bvec_tmp(jbvec))
             if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Norm of jbvec', jbvec, T1
             T2 = mat_dotproduct(Bvec_tmp(jbvec),Bvec_tmp(ibvec))
             if (molcfg%solver%info_rsp) write(molcfg%lupri,*) '<bjbvec|bibvec>', T2
             TT = -T2/T1
             if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'TT = -T2/T1', TT
             !
             call mat_daxpy(TT,Bvec_tmp(jbvec),Bvec_tmp(ibvec))
             jbvec = jbvec + 1
             !    endif
          enddo
          !    endif
       enddo
       do ibvec = 1,Nm_new !index for current bvector 
          if (molcfg%solver%info_rsp) write(molcfg%lupri,*)'Current ibvec index is ', ibvec, lin_depend(ibvec)
          !   if (lin_depend(ibvec) /= 0) then
          jbvec = 1    !index for another bvector
          do jrx = 1,(ibvec-1)
             !      if (lin_depend(jrx) == 0) then
             !       jbvec = jbvec + 1 !skip
             !     else
             T1 = mat_sqnorm2(Bm_tmp(jbvec))
             if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Norm of jbvec', jbvec, T1
             T2 = mat_dotproduct(Bm_tmp(jbvec),Bm_tmp(ibvec))
             if (molcfg%solver%info_rsp) write(molcfg%lupri,*) '<bjbvec|bibvec>', T2
             TT = -T2/T1
             if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'TT = -T2/T1', TT
             !
             call mat_daxpy(TT,Bm_tmp(jbvec),Bm_tmp(ibvec))
             jbvec = jbvec + 1
             !    endif
          enddo
          ! endif
       enddo
       !
       ! NORMALIZE VECTOR NUMBER Ibvec
       !
       do ibvec=1,nb_new
          !  if (lin_depend(ibvec) /= 0) then !not linear dependent
          if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'RSP_SOLVER: Orthonormalize -> normalize'
          if (lineq) then
             call rsp_normalize(molcfg,iturn,ibvec,run_ortho_again,lin_depend(ibvec),Bvec_tmp(ibvec))
          else
             call normalize3(molcfg,iturn,ibvec,run_ortho_again,lin_depend(ibvec),Bvec_tmp(ibvec),.true.)
          endif

          !  endif
       enddo !ibvec
       do ibvec=1,nm_new
          ! if (lin_depend(ibvec) /= 0) then !not linear dependent
          if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'RSP_SOLVER: Orthonormalize -> normalize'
          if (lineq) then
             call rsp_normalize(molcfg,iturn,ibvec,run_ortho_again,lin_depend(ibvec),Bm_tmp(ibvec))
          else
             call normalize3(molcfg,iturn,ibvec,run_ortho_again,lin_depend(ibvec),Bm_tmp(ibvec),.false.)
          endif
          ! endif
       enddo !ibvec

       if (run_ortho_again) then
          run_ortho_again = .false.
          if (iturn == 2) then !This is redundant since run_ortho_again is only set .true. if iturn = 1
             WRITE(molcfg%LUPRI,'(/A)') &
                  &     'Error: Already ran twice through orthonormalize!21'
             CALL lsQUIT('Error: Already ran twice through orthonormalize21!',molcfg%lupri)
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
       !  if (lin_depend(irx) /= 0) then
       no_of_new = no_of_new + 1
       call mat_assign(bvecs(no_of_new),Bvec_tmp(irx))
       !  endif
    enddo
    !
    !     Set NB_NEW to actual number of acceptable new trial vectors
    !

    Nb_new = no_of_new


    no_of_new = 0
    do irx = 1,Nm_new
       !   if (lin_depend(irx) /= 0) then
       no_of_new = no_of_new + 1
       call mat_assign(bvecsm(no_of_new),Bm_tmp(irx))
       !   endif
    enddo
    !
    !     Set NB_NEW to actual number of acceptable new trial vectors
    !

    Nm_new = no_of_new

    DEALLOCATE(lin_depend)
    call mat_free(B_scr)
    call mat_free(b_k)
    call mat_free(bT)
    call mat_free(bp)
    call mat_free(bm)

  end subroutine orthonormalize21
  subroutine normalize3(molcfg,iturn,i,run_ortho_again,lin_depend_i,Bvec_tmp_i,symm)
!Purpose:
!   Normalize the i'th bvector
!
!   iturn: Which time we are running through orthonormalize2, which calls normalize2 (IN)
!   i    :  index of Bvec_temp_i among the new bvectors (for printout)               (IN)
!   run_ortho_again: If linear dependencies are suspected, run_ortho_again is set    (OUT)
!                    to true, and orthonormalize2 is run again
!   lin_depend_i: linear dependence indicator: 1 if non-dep, set to 0 if lin dep     (OUT)
!   Bvec_temp_i,  the i'th bvector to be normalized				     (INOUT)
    implicit none
    type(rsp_molcfg), intent(inout) :: molcfg
    integer, intent(in) :: iturn,i
    integer, intent(inout) :: lin_depend_i
    type(Matrix), intent(inout) :: Bvec_tmp_i
    real(realk) :: TT
    logical    :: symm
    logical, intent(inout) :: run_ortho_again
    integer  :: j,k

    TT = mat_sqnorm2(Bvec_tmp_i)
    if (TT==0E0_realk) then
            if (symm) then
             do j = 1,Bvec_tmp_i%ncol
                do k = 1,Bvec_tmp_i%nrow
                   if (k == j) then
                       Bvec_tmp_i%elms((j-1)*Bvec_tmp_i%nrow+k) = 0E0_realk
                   else 
                       if (k>j) then
                        Bvec_tmp_i%elms((j-1)*Bvec_tmp_i%nrow+k) = 1.0E0_realk
                        Bvec_tmp_i%elms((k-1)*Bvec_tmp_i%nrow+j)=Bvec_tmp_i%elms((j-1)*Bvec_tmp_i%nrow+k)
                       endif
                   endif
                 enddo
               enddo
            else
             do j = 1,Bvec_tmp_i%ncol
                do k = 1,Bvec_tmp_i%nrow
                   if (k == j) then
                       Bvec_tmp_i%elms((j-1)*Bvec_tmp_i%nrow+k) = 0E0_realk
                   else 
                       if (i>j) then
                        Bvec_tmp_i%elms((j-1)*Bvec_tmp_i%nrow+k) = 1.0E0_realk
                        Bvec_tmp_i%elms((k-1)*Bvec_tmp_i%nrow+j)=-Bvec_tmp_i%elms((j-1)*Bvec_tmp_i%nrow+k)
                       endif
                   endif
                 enddo
               enddo
            endif
    call mat_identity(Bvec_tmp_i) 
 !   elseif (TT < molcfg%solver%rsp_thr_round) then
!      if (lineq) then
!   call mat_scal(1/TT,Bvec_tmp_i)
!         run_ortho_again=.true.
!       endif
endif
     !if (molcfg%solver%info_rsp) then
!        write(molcfg%lupri,*)'bvector number ', i, ' is removed because its norm**2 :', TT, &
!              & ' is < cfg_rsp_thr_lin_depend '!, molcfg%rsp_thr_lin_depend
      !endif
 !     lin_depend_i = 0  !linear dependent
 !   elseif(TT < molcfg%solver%rsp_thr_round) then
 !     if (iturn == 1) then
 !        run_ortho_again = .true.
 !    else
 !   call mat_scal(1E0_realk/sqrt(TT),Bvec_tmp_i)
 !        run_ortho_again = .true.
        !if (molcfg%solver%info_rsp) then
  !        write(molcfg%lupri,*) 'bvector number ', i, ' is removed because its norm**2 :', TT, &
  !            & ' is < cfg_rsp_thr_round after 2nd Gram-Schmidt orth.'
!, molcfg%rsp_thr_round ,' 
!        endif
!       lin_depend_i = 0
!      endif
    !END IF
    !Vector is not linearly denpendent - normalize it:
 !   IF (lin_depend_i /= 0) THEN
      IF (molcfg%solver%info_rsp) then
          write(molcfg%lupri,*)'Now normalize bvector number ', i, ' with initial norm**2 :', TT
      endif
     !If norm is very small, increase it by scaling bvector
      IF (TT < 1E-5_realk) THEN
         TT = 1E0_realk / SQRT(TT)
         call mat_SCAL(TT,BVEC_tmp_i)
         TT = mat_sqnorm2(Bvec_tmp_i)
      END IF
      !Normalization
      TT = 1E0_realk / SQRT(TT)
      call mat_SCAL(TT,BVEC_tmp_i) 
!    END IF

  end subroutine normalize3

  subroutine ReorderEigenvalues2(molcfg,NDIM,SRED,EIVEC,ALFAR,ALFAI,BETA,WRK1,WRK2,ISNDX)
!
! Analyze and reorder eigenvectors from solve_red_eigen
! Based on RSPORD
!
  implicit none
  type(rsp_molcfg), intent(inout) :: molcfg
  integer, intent(in) :: NDIM
  integer, intent(inout) :: ISNDX(3)
  real(realk),intent(in) :: SRED(NDIM,NDIM) 
  !real(realk),intent(in) :: EIVEC(NDIM,NDIM)  !SONIA: FIXME
  real(realk),intent(inout) :: EIVEC(NDIM,NDIM)
  real(realk),intent(inout) :: WRK1(NDIM,NDIM), WRK2(NDIM,NDIM)
  real(realk),intent(inout) :: ALFAR(:),ALFAI(:),BETA(:)
!
! eigenvalue k is (ALFAR(k)+i*ALFAI(k))/BETA(k)
!
  integer :: ipos, izer, ineg, i, j, jmin, nneg
  real(realk) :: xscale,xsave, amin
  real(realk), PARAMETER :: ZEROT=1.0E-9_realk, COMPLX = 1.0E7_realk
  real(realk), PARAMETER :: D0=0.0E0_realk, D1=1.0E0_realk
!
  if (molcfg%solver%info_rsp) then
         WRITE (molcfg%LUPRI,'(//A/A)') &
     &      '  (ALFAR + i ALFAI) / BETA are eigenvalues;', &
     &      '     ALFAR           ALFAI          BETA' 
         WRITE (molcfg%LUPRI,'(1P,3D15.6)')(ALFAR(I),ALFAI(I),BETA(I),I=1,NDIM)
  end if
  DO I=1,NDIM
     IF(ABS(BETA(I)).GE.ZEROT) THEN
       ALFAR(I)=ALFAR(I)/BETA(I)
       ALFAI(I)=ALFAI(I)/BETA(I)
     ELSE
!      singularities
       WRITE(molcfg%LUPRI,1010)I, ALFAR(I),ALFAI(I),BETA(I)
 1010  FORMAT(/' *** SINGULARITY IN REDUCED MCRSP :', &
     &             /'     I,ALFAR(I),ALFAI(I),BETA(I):',I6,1P,3D13.6)
     END IF
!        IF(ABS(ALFAR(I)).GT.ZEROT) THEN
!    &   .OR.(ABS(ALFAR(I)).LE.ZEROT.AND.ABS(ALFAI(I)).GT.ZCRIT)) THEN
!           complex eigenvalues
!        RATIO= ABS(ALFAI(I)/ALFAR(I))
!        IF(RATIO.GT.ZCRIT) WRITE(molcfg%LUPRI,1020) ALFAR(I),ALFAI(I)
     IF(ABS(ALFAI(I)).GT.D0) THEN
       WRITE(molcfg%LUPRI,1020) I,ALFAR(I),ALFAI(I)
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
         WRITE (molcfg%LUPRI,'(/A)') ' Eigenvalues of E(2)  :'
         WRITE (molcfg%LUPRI,'(I10,1P,D12.2)') (I,ALFAR(I),I=1,NDIM)
  endif
!
!     reduced S(2) in eigenvector basis
!
  call ls_dzero(WRK1,NDIM*NDIM)
  call ls_dzero(WRK2,NDIM*NDIM)
  if (molcfg%solver%info_rsp) then
    WRITE (molcfg%LUPRI,*) ' Reduced S(2) before diagonal basis '
    CALL LS_OUTPUT(SRED,1,NDIM,1,NDIM,NDIM,NDIM,1,molcfg%LUPRI)
    WRITE (molcfg%LUPRI,*) ' Reduced Xvec before diagonal basis '
    CALL LS_OUTPUT(EIVEC,1,NDIM,1,NDIM,NDIM,NDIM,1,molcfg%LUPRI)
  endif

  CALL DGEMM('N','N',NDIM,NDIM,NDIM,1E0_realk,SRED,NDIM, &
     &           EIVEC,NDIM,0E0_realk,WRK1,NDIM)
  if (molcfg%solver%info_rsp) then
    WRITE (molcfg%LUPRI,*) ' REDS[2]*X'
    CALL LS_OUTPUT(WRK1,1,NDIM,1,NDIM,NDIM,NDIM,1,molcfg%LUPRI)
  endif
  CALL DGEMM('T','N',NDIM,NDIM,NDIM,1E0_realk,EIVEC,NDIM, &
     &           WRK1,NDIM,0E0_realk,WRK2,NDIM)
  if (molcfg%solver%info_rsp) then
         WRITE (molcfg%LUPRI,'(/A)') ' Reduced S(2) in diagonal basis :'
         WRITE (molcfg%LUPRI,'(I10,1P,D12.2)') (I,WRK2(I,I),I=1,NDIM)
!        CALL LS_OUTPUT(WRK2,1,KZYRED,1,KZYRED,KZYRED,KZYRED,1,molcfg%LUPRI)
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
   IF (IPOS.NE.(NDIM/2)) WRITE (molcfg%LUPRI,2020) IPOS,IZER,INEG
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
      WRITE(molcfg%LUPRI,'(3(/A))')' ***** WARNING *********' &
     &   ,' number of eigenvalues with negative metric differ from' &
     &   ,' number with positive metric'
      !IF (IPRRSP.GT. 20) THEN
      if (molcfg%solver%info_rsp) then
      WRITE(molcfg%LUPRI,'(/A)')'   NUMBER    EIGENVALUE '
      DO I=1,NDIM
         WRITE(molcfg%LUPRI,'(I10,5X,1P,D15.8)') I,ALFAR(I)
      END DO
      end if
      !END IF
   ELSE
     if (molcfg%solver%info_rsp) WRITE(molcfg%LUPRI,'(/A)') &
     &      '   NUMBER    EIGENVALUE    PAIRED EIGENVALUE'
     DO I=1,IPOS
        IF (ABS(ABS(ALFAR(I))-ABS(ALFAR(IPOS+I))).GT.ZEROT) THEN
            WRITE(molcfg%LUPRI,'(/A)')' **WARNING** EIGENVALUES NOT PAIRED'
        END IF
        if (molcfg%solver%info_rsp) &
        !IF (IPRRSP.GT. 20)
     &      WRITE(molcfg%LUPRI,'(I10,5X,1P,D15.8,5X,D15.8)') &
     &                      I,ALFAR(I),ALFAR(IPOS+I)
     END DO
     IF (IZER .GT. 0 ) THEN
        WRITE(molcfg%LUPRI,'(/A)') &
     &         ' ZERO METRIC EIGENVALUE AND ZERO EIGENVALUE'
        WRITE(molcfg%LUPRI,'(/A)') 'NUMBER    EIGENVALUE'
        DO I=1,IZER
           WRITE(molcfg%LUPRI,'(I10,5X,1P,D15.8)') I,ALFAR(IPOS+INEG+I)
        END DO
     END IF
   END IF
!
   end subroutine ReorderEigenvalues2
 end module RSPSYMSOLVER
