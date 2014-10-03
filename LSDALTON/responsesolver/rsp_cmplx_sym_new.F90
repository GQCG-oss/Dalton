module COMPLEXNEWSYMSOLVER
  use memory_handling  
  use decompMod
  use precision
  use matrix_module
  use matrix_operations
  use RSPSOLVER
  use RSPSYMSOLVER
  use COMPLEXSOLVER
  use rsp_util
  use files
  !use configuration
  private
  public ::  new_symcomplex_solver
  integer, save :: rsp_number_of_current_trial, rsp_number_of_rhs, rsp_number_of_sols, &
                 & rsp_number_of_omegas, rsp_number_of_startvecs,rsp_bvec_dim, &
                 & lu_x_s_rsp, lu_sigma_s_rsp, lu_rho_a_rsp,&
                 & lu_x_a_rsp, lu_sigma_a_rsp, lu_rho_s_rsp
                 !to complex solver
  real(realk),allocatable,save :: E1(:,:), E2(:,:), S1(:,:),S2(:,:), G1(:,:), G2(:,:)
  real(realk),allocatable,save :: G3(:,:), G4(:,:)
  logical :: RSPonMaster
contains
!==========================================================
!==============================================================
    subroutine new_symcomplex_solver(molcfg,F,D,S,ngd,nfreq,GDB,GDI,EIVAL,XSOL,XSOLimg,gd_complex)
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
      ! Joanna, October 2013
      !=========================================================================
      !
      use RSPSYMSOLVER
      implicit none
      type(rsp_molcfg), intent(inout)    :: molcfg
      type(Matrix),intent(in)             :: F,D,S
      integer, intent(in)                 :: ngd,nfreq
      type(Matrix),intent(inout)          :: GDB(ngd)
      type(Matrix),intent(inout),optional :: GDI(ngd)
      real(realk),intent(inout)           :: eival(nfreq)
      type(Matrix),intent(inout)          :: XSOL(nfreq,ngd),XSOLimg(nfreq,ngd)
      logical, intent(in)                 :: gd_complex
!
      real(realk),allocatable             :: red_Xp_glob(:,:,:),red_Xm_glob(:,:,:)
      real(realk),allocatable             :: red_Xpi_glob(:,:,:),red_Xmi_glob(:,:,:)
      type(Matrix)                        :: xp(nfreq,ngd),xm(nfreq,ngd),xpi(nfreq,ngd),xmi(nfreq,ngd)
      type(Matrix)                        :: RHS_real(ngd),RHS_img(ngd)
      type(Matrix),pointer                :: sigmass(:), rhosa(:), sigmasa(:), rhoss(:)
      type(Matrix)                        :: xsvec(2*nfreq*ngd),xavec(2*nfreq*ngd)
      logical                             :: conv(nfreq,ngd)
      integer                             :: max_red, max_it, ndim, i, j, nx_new, ns_red, na_red, nxs_new, nxa_new
      integer                             :: itmic,conv_tot
      real(realk)                         :: gammma, rsp_thresh
!
      max_red= molcfg%solver%rsp_maxred
      max_it= molcfg%solver%rsp_maxit
      rsp_thresh= molcfg%solver%rsp_thresh
      gammma=molcfg%solver%rsp_gamma

      ndim = S%nrow
      lu_x_s_rsp= -1 ; lu_sigma_s_rsp = -1 ; lu_rho_a_rsp = -1
      lu_x_a_rsp= -1 ; lu_sigma_a_rsp = -1 ; lu_rho_s_rsp = -1
      CALL LSOPEN(lu_x_s_rsp,'xsrsp','unknown','UNFORMATTED')
      CALL LSOPEN(lu_sigma_s_rsp,'sigmas_rsp','unknown','UNFORMATTED')
      CALL LSOPEN(lu_rho_a_rsp,'rhoa_rsp','unknown','UNFORMATTED')
      CALL LSOPEN(lu_x_a_rsp,'xarsp','unknown','UNFORMATTED')
      CALL LSOPEN(lu_sigma_a_rsp,'sigmaa_rsp','unknown','UNFORMATTED')
      CALL LSOPEN(lu_rho_s_rsp,'rhos_rsp','unknown','UNFORMATTED')

      allocate(red_Xp_glob(max_red,nfreq,ngd),red_Xm_glob(max_red,nfreq,ngd))
      allocate(red_Xpi_glob(max_red,nfreq,ngd),red_Xmi_glob(max_red,nfreq,ngd))
      ALLOCATE(E1(max_red,max_red),E2(max_red,max_red))
      allocate(S1(max_red,max_red),S2(max_red,max_red))
      allocate(G1(max_red,ngd))  
      allocate(G2(max_red,ngd))  
      allocate(G3(max_red,ngd))  
      allocate(G4(max_red,ngd))  

      red_Xp_glob = 0.0E0_realk ; red_Xm_glob = 0.0E0_realk
      red_Xpi_glob = 0.0E0_realk ; red_Xmi_glob = 0.0E0_realk
      E1=0E0_realk; E2=0E0_realk; S1=0E0_realk; S2=0E0_realk; G1=0E0_realk; G2=0E0_realk
      G3=0E0_realk; G4=0E0_realk
      
      nx_new=0
      ns_red=0
      na_red=0
      
      do j=1,ngd
         call mat_init(RHS_real(j),ndim,ndim)
         call mat_init(RHS_img(j),ndim,ndim)
         call mat_zero(RHS_real(j))
         call mat_zero(RHS_img(j))
         call mat_assign(RHS_real(j),GDb(j))
         if (gd_complex) then
            call mat_assign(RHS_img(j),GDi(j))
         endif
         do i=1,nfreq
           conv(i,j)=.false.
           call mat_init(xp(i,j),ndim,ndim)
           call mat_init(xm(i,j),ndim,ndim)
           call mat_init(xpi(i,j),ndim,ndim)
           call mat_init(xmi(i,j),ndim,ndim)
           call mat_zero(xp(i,j))
           call mat_zero(xm(i,j))
           call mat_zero(xpi(i,j))
           call mat_zero(xmi(i,j))
           call mat_zero(xsol(i,j))   
           call mat_zero(xsolimg(i,j))   
         enddo
      enddo

      call get_start_vectors2(molcfg,F,D,S,RHS_real,RHS_img,ngd,nfreq,eival,nx_new,xp,xm,xpi,xmi)

      do itmic = 1,max_it
         write(molcfg%lupri,*) '------------------'
         write(molcfg%lupri,*) 'Start macroiteration:',itmic
         write(molcfg%lupri,*) '------------------'
         do j=1,2*nx_new
            call mat_init(xsvec(j),ndim,ndim)
            call mat_init(xavec(j),ndim,ndim)
            call mat_zero(xsvec(j))
            call mat_zero(xavec(j))
         enddo
         nxs_new=nx_new
         nxa_new=nx_new
         call orthonormalize_sym(molcfg,D,S,xp,xpi,nxs_new,ns_red,xsvec,lu_x_s_rsp)
         call orthonormalize_sym(molcfg,D,S,xm,xmi,nxa_new,na_red,xavec,lu_x_a_rsp)
          nullify(sigmass)
          nullify(sigmasa)
          nullify(rhosa)
          nullify(rhoss)

!FIXME Thomas & Joanna
         call transform_vectors(molcfg,D,S,F,nxs_new,xsvec,sigmass,rhosa,.true.)
         call transform_vectors(molcfg,D,S,F,nxa_new,xavec,sigmasa,rhoss,.true.)
!end FIXME Thomas & Joanna

        if (gd_complex) then
         !build a reduce space
          call extend_new_complex_reduced_matrices2(molcfg,nxs_new,nxa_new,ns_red,na_red,&
                 &ngd,RHS_real,xsvec,sigmass,rhosa,xavec,sigmasa,rhoss,.true.,RHS_img)
            !write on disk. The sum of symm and antisymm component is stored together
            do j = 1,nxs_new
               call mat_write_to_disk(lu_x_s_rsp,xsvec(j),RSPonMaster)
               call mat_write_to_disk(lu_sigma_s_rsp,sigmass(j),RSPonMaster)
               call mat_write_to_disk(lu_rho_a_rsp,rhosa(j),RSPonMaster)
            enddo
            ns_red=ns_red+nxs_new
            do j = 1,nxa_new
               call mat_write_to_disk(lu_x_a_rsp,xavec(j),RSPonMaster)
               call mat_write_to_disk(lu_sigma_a_rsp,sigmasa(j),RSPonMaster)
               call mat_write_to_disk(lu_rho_s_rsp,rhoss(j),RSPonMaster)
            enddo
            na_red=na_red+nxa_new
            do j=1,2*nx_new
               call mat_free(xsvec(j))
               call mat_free(xavec(j))
             enddo
            do j=1,nxs_new
               call mat_free(sigmass(j))
               call mat_free(rhosa(j))
            enddo
            do j=1,nxa_new
               call mat_free(sigmasa(j))
               call mat_free(rhoss(j))
            enddo
            call mem_dealloc(sigmass)
            call mem_dealloc(sigmasa)
            call mem_dealloc(rhoss)
            call mem_dealloc(rhosa)
            !        
            !          !solve complex reduced equation
            call solve_complex2(molcfg,ndim,ns_red,na_red,ngd,nfreq,eival,&
                 &red_Xp_glob,red_Xpi_glob,red_Xm_glob,red_Xmi_glob,.true.)
            !    
            !        ! check for convergece and get new trialvectors
            do j=1,ngd
              do i=1,nfreq
               call mat_zero(xp(i,j))
               call mat_zero(xm(i,j))
               call mat_zero(xpi(i,j))
               call mat_zero(xmi(i,j))
              enddo
            enddo
            call get_complex_res2(molcfg,F,D,S,itmic,ns_red,na_red,nfreq,&
                 &ngd,RHS_real,red_Xp_glob,red_Xpi_glob,red_Xm_glob,&
                 &red_Xmi_glob,eival,xp,xm,xpi,xmi,nx_new,conv,conv_tot,xsol,xsolimg,&
                 &.true.,RHS_img)
         else
            !build a reduce space
          call extend_new_complex_reduced_matrices2(molcfg,nxs_new,nxa_new,ns_red,na_red,&
                 &ngd,RHS_real,xsvec,sigmass,rhosa,xavec,sigmasa,rhoss,.false.)

            do j = 1,nxs_new
               call mat_write_to_disk(lu_x_s_rsp,xsvec(j),RSPonMaster)
               call mat_write_to_disk(lu_sigma_s_rsp,sigmass(j),RSPonMaster)
               call mat_write_to_disk(lu_rho_a_rsp,rhosa(j),RSPonMaster)
            enddo
            ns_red=ns_red+nxs_new
            do j = 1,nxa_new
               call mat_write_to_disk(lu_x_a_rsp,xavec(j),RSPonMaster)
               call mat_write_to_disk(lu_sigma_a_rsp,sigmasa(j),RSPonMaster)
               call mat_write_to_disk(lu_rho_s_rsp,rhoss(j),RSPonMaster)
            enddo
            na_red=na_red+nxa_new
            do j=1,2*nx_new
               call mat_free(xsvec(j))
               call mat_free(xavec(j))
             enddo
            do j=1,nxs_new
               call mat_free(sigmass(j))
               call mat_free(rhosa(j))
            enddo
            do j=1,nxa_new
               call mat_free(sigmasa(j))
               call mat_free(rhoss(j))
            enddo

            !solve complex reduced equation
            call solve_complex2(molcfg,ndim,ns_red,na_red,ngd,nfreq,eival,&
                 &red_Xp_glob,red_Xpi_glob,red_Xm_glob,red_Xmi_glob,.false.)

            ! check for convergece and get new trialvectors
            do j=1,ngd
              do i=1,nfreq
                call mat_zero(xp(i,j))
                call mat_zero(xm(i,j))
                call mat_zero(xpi(i,j))
                call mat_zero(xmi(i,j))
              enddo
            enddo
            call get_complex_res2(molcfg,F,D,S,itmic,ns_red,na_red,nfreq,&
                 &ngd,RHS_real,red_Xp_glob,red_Xpi_glob,red_Xm_glob,&
                 &red_Xmi_glob,eival,xp,xm,xpi,xmi,nx_new,conv,conv_tot,xsol,xsolimg,&
                 &.false.)
         endif

         if (conv_tot==1)  then
            do j=1,ngd
              do i=1,nfreq
               call mat_free(xp(i,j))
               call mat_free(xpi(i,j))
               call mat_free(xm(i,j))
               call mat_free(xmi(i,j))
              enddo
            enddo
            exit
         endif         
      enddo
      do i=1,ngd
         call mat_free(RHS_real(i))
         call mat_free(RHS_img(i))
      enddo
      deallocate(E1,E2,S1,S2,G1,G2,G3,G4)
      deallocate(red_Xp_glob,red_Xm_glob)
      deallocate(red_Xpi_glob,red_Xmi_glob)

      CALL lsCLOSE(lu_sigma_s_rsp,'DELETE')
      CALL lsCLOSE(lu_rho_a_rsp,'DELETE')
      CALL lsCLOSE(lu_x_s_rsp,'DELETE')
      CALL lsCLOSE(lu_sigma_a_rsp,'DELETE')
      CALL lsCLOSE(lu_rho_s_rsp,'DELETE')
      CALL lsCLOSE(lu_x_a_rsp,'DELETE')

      !================================================================  
    end subroutine new_symcomplex_solver
!================================================================  
!================================================================  
subroutine extend_new_complex_reduced_matrices2(molcfg,nxs_new,nxa_new,ns_red,na_red,&
          &ngd,gd,xsvec,sigmass,rhosa, xavec,sigmasa,rhoss,gd_complex,gdi)
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
type(rsp_molcfg), intent(inout)    :: molcfg
integer,intent(in)            :: ns_red,na_red
integer, intent(in)              :: nxs_new,nxa_new,ngd
type(Matrix),intent(in)          :: xsvec(:),sigmass(:),rhosa(:)
type(Matrix),intent(in)          :: xavec(:),sigmasa(:),rhoss(:)
type(Matrix),intent(in)          :: gd(:)
type(Matrix),intent(in),optional :: gdi(:)
!
type(Matrix)                     :: gdp,gdpi,gdm,gdmi
logical,intent(in)               :: gd_complex
integer                          :: ndim,k,j,i,l,max_red
type(Matrix)                     ::x_i,sigmai_j,rhoi_j,x_j,rho_j,sigma_j,gT
type(Matrix)                     ::xpi_j,sigmapi_j,xp_j,sigmap_j,scrT
type(Matrix)                     ::xmi_j,sigmami_j,rhopi_j,xm_j,rhop_j,sigmam_j

 max_red= molcfg%solver%rsp_maxred
if ((ns_red + nxs_new) > max_red) then
   WRITE (molcfg%LUPRI,*) 'ERROR in extend_complex_reduced_matrices: ns_red + nxs_new > rsp_maxred'
   WRITE (molcfg%LUPRI,*) 'ns_red ,nxs_new  =',ns_red,nxs_new
   WRITE (molcfg%LUPRI,*) 'rsp_maxred =',max_red
   WRITE (molcfg%LUPRI,*) 'Reduce problem size or recompile orthonormalize '// &
        &                  'with larger rsp_maxred parameter'
   CALL lsQUIT('orthonormalize error: ns_red + nxs_new > rsp_maxred parameter',molcfg%lupri)
END IF
if ((na_red+nxa_new) > max_red) then
   WRITE (molcfg%LUPRI,*) 'ERROR in extend_complex_reduced_matrices: na_red + nxa_new > rsp_maxred'
   WRITE (molcfg%LUPRI,*) 'na_red ,nxa_new  =',na_red,nxa_NEW
   WRITE (molcfg%LUPRI,*) 'rsp_maxred =',max_red
   WRITE (molcfg%LUPRI,*) 'Reduce problem size or recompile orthonormalize '// &
        &                  'with larger rsp_maxred parameter'
   CALL lsQUIT('orthonormalize error: na_red + nxa_NEW > rsp_maxred parameter',molcfg%lupri)
END IF

ndim = gd(1)%nrow
call mat_init(rho_j,ndim,ndim)
call mat_init(x_i,ndim,ndim)

if (nxs_new>0) then
    do i = 1, nxs_new 
       do j = i, nxs_new
          E1(ns_red+i,ns_red+j)   =  mat_dotproduct(xsvec(i),sigmass(j))
          if (i .ne. j) E1(ns_red+j,ns_red+i)=E1(ns_red+i,ns_red+j)
       enddo
       if (nxa_new>0) then      
         do j = 1, nxa_new
           S1(ns_red+i,na_red+j)   =  mat_dotproduct(xsvec(i),rhoss(j))
         enddo
       endif      
    enddo
endif
if (nxa_new>0) then
    do i = 1, nxa_new 
       do j = i, nxa_new
          E2(na_red+i,na_red+j)   =  mat_dotproduct(xavec(i),sigmasa(j))
          if (i .ne. j) E2(na_red+j,na_red+i)=E2(na_red+i,na_red+j)
       enddo      
       if (nxs_new>0) then      
         do j = 1, nxs_new
           S2(na_red+i,ns_red+j)   =  mat_dotproduct(xavec(i),rhosa(j))
         enddo
       endif      
    enddo
endif
if (ns_red>0) then
    !Setup lower half of E2 and S2:
    rewind(lu_x_s_rsp) 
    do i = 1, ns_red
       call mat_read_from_disk(lu_x_s_rsp,x_i,RSPonMaster)
       if (nxs_new>0) then
         do j = 1, nxs_new
           E1(i,ns_red+j)   =  mat_dotproduct(x_i,sigmass(j))
           E1(ns_red+j,i) = E1(i,ns_red+j)
         enddo
       endif
       if (nxa_new>0) then
         do j = 1, nxa_new
           S1(i,na_red+j)   =  mat_dotproduct(x_i,rhoss(j))
         enddo
       endif
   enddo
    rewind(lu_rho_s_rsp)
    call mat_zero(rho_j)
    do j = 1, na_red
       call mat_read_from_disk(lu_rho_s_rsp,rho_j,RSPonMaster)
      if (nxs_new>0) then 
       do i = 1,nxs_new
          S1(ns_red+i,j)   =  mat_dotproduct(xsvec(i),rho_j)
       enddo
      endif
    enddo
endif
if (na_red>0) then
    rewind(lu_x_a_rsp) 
    do i = 1, na_red
       call mat_read_from_disk(lu_x_a_rsp,x_i,RSPonMaster)
       
       if (nxa_new>0) then
         do j = 1, nxa_new
           E2(i,na_red+j)   =  mat_dotproduct(x_i,sigmasa(j))
           E2(na_red+j,i) = E2(i,na_red+j)
         enddo
       endif
       if (nxs_new>0) then
         do j = 1, nxs_new
           S2(i,ns_red+j)   =  mat_dotproduct(x_i,rhosa(j))
         enddo
       endif
   enddo
    rewind(lu_rho_a_rsp)
    call mat_zero(rho_j)
    do j = 1, ns_red
       call mat_read_from_disk(lu_rho_a_rsp,rho_j,RSPonMaster)
      if (nxa_new>0) then 
       do i = 1,nxa_new
          S2(na_red+i,j)   =  mat_dotproduct(xavec(i),rho_j)
       enddo
      endif
    enddo
  endif 
    
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
       do i = 1,nxs_new
         G1(ns_red+i,j) = mat_dotproduct(xsvec(i),gdp)
         if (gd_complex) then
           G4(ns_red+i,j) = mat_dotproduct(xsvec(i),gdpi)
         endif
       enddo
       do i = 1,nxa_new
         G2(na_red+i,j) = mat_dotproduct(xavec(i),gdm)
         if (gd_complex) then
           G3(na_red+i,j) = mat_dotproduct(xavec(i),gdmi)
         endif
       enddo
    enddo
       call mat_free(gdp)
       call mat_free(gdm)
       if (gd_complex) then
       call mat_free(gdpi)
       call mat_free(gdmi)

       endif
        !save on disk
    if (molcfg%solver%info_rsp_redspace) then
     write(molcfg%lupri,*) 'Reduced E1:'
     call LS_OUTPUT(E1, 1, ns_red+nxs_new, 1, ns_red+nxs_new, max_red, max_red, 1, molcfg%lupri)
     write(molcfg%lupri,*) 'Reduced E2:'
     call LS_OUTPUT(E2, 1, na_red+nxa_new, 1, na_red+nxa_new, max_red, max_red, 1, molcfg%lupri)
     
     write(molcfg%lupri,*) 'Reduced S1:'
     call LS_OUTPUT(S1, 1, ns_red+nxs_new, 1, na_red+nxa_new, max_red, max_red, 1, molcfg%lupri)
     write(molcfg%lupri,*) 'Reduced S2:'
     call LS_OUTPUT(S2, 1, na_red+nxa_new, 1, ns_red+nxs_new, max_red, max_red, 1, molcfg%lupri)
     
     write(molcfg%lupri,*) 'Reduced G1:'
     call LS_OUTPUT(G1, 1, ns_red+nxs_new, 1, ngd, max_red, max_red, 1, molcfg%lupri)
     write(molcfg%lupri,*) 'Reduced G2:'
     call LS_OUTPUT(G2, 1, na_red+nxa_new, 1, ngd, max_red, max_red, 1, molcfg%lupri)
     if (gd_complex) then
        write(molcfg%lupri,*) 'Reduced G3:'
        call LS_OUTPUT(G3, 1, na_red+nxa_new, 1, ngd, max_red, max_red, 1, molcfg%lupri)
        write(molcfg%lupri,*) 'Reduced G4:'
        call LS_OUTPUT(G4, 1, ns_red+nxs_new, 1, ngd, max_red, max_red, 1, molcfg%lupri)
     endif
   endif

call mat_free(x_i)
call mat_free(rho_j)
    end subroutine extend_new_complex_reduced_matrices2

!================================================================  
 !================================================================  
subroutine solve_complex2(molcfg,ndim,ns_red,na_red,ngd,nfreq,freq,red_X,red_Xi,red_Xm,red_Xmi,gd_complex)
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
type(rsp_molcfg), intent(inout)    :: molcfg
    integer, intent(in)      :: ns_red, ngd,ndim,na_red,nfreq
    real(realk),intent(in)   :: freq(:)
    logical                  :: gd_complex
    real(realk),intent(inout):: red_X(:,:,:),red_Xm(:,:,:)
    real(realk),intent(inout):: red_Xi(:,:,:),red_Xmi(:,:,:)
    real(realk),allocatable   :: A(:,:),KHS(:)
    integer,allocatable      :: IPIV(:)
    integer                  :: igd, ifreq,ierr,n_red,max_red
    real(realk)              :: gammma

    gammma=molcfg%solver%rsp_gamma
    max_red=molcfg%solver%rsp_maxred
    n_red=2*ns_red+2*na_red

    allocate(KHS(n_red), IPIV(n_red))
    allocate(A(n_red,n_red))
   
  do igd=1,ngd 
    do ifreq=1,nfreq   
      KHS=0E0_realk
      IPIV=0E0_realk
      A=0E0_realk

      A(1:ns_red,1:ns_red)                                                       =  E1(1:ns_red,1:ns_red)
      A(ns_red+1:ns_red+na_red,ns_red+1:ns_red+na_red)                           =  E2(1:na_red,1:na_red)
      A(ns_red+na_red+1:ns_red+2*na_red,ns_red+na_red+1:ns_red+2*na_red)         = -E2(1:na_red,1:na_red)
      A(ns_red+2*na_red+1:2*ns_red+2*na_red,ns_red+2*na_red+1:2*ns_red+2*na_red) = -E1(1:ns_red,1:ns_red)

      A(1:ns_red,ns_red+1:ns_red+na_red)                                         = -freq(ifreq)*S1(1:ns_red,1:na_red)
      A(1:ns_red,ns_red+na_red+1:ns_red+2*na_red)                                =  gammma*S1(1:ns_red,1:na_red)
      A(ns_red+1:ns_red+na_red,1:ns_red)                                         = -freq(ifreq)*S2(1:na_red,1:ns_red)
      A(ns_red+1:ns_red+na_red,ns_red+2*na_red+1:2*ns_red+2*na_red)              =  gammma*S2(1:na_red,1:ns_red)
      A(ns_red+na_red+1:ns_red+2*na_red,1:ns_red)                                =  gammma*S2(1:na_red,1:ns_red)
      A(ns_red+na_red+1:ns_red+2*na_red,ns_red+2*na_red+1:2*ns_red+2*na_red)   =  freq(ifreq)*S2(1:na_red,1:ns_red)
      A(ns_red+2*na_red+1:2*ns_red+2*na_red,ns_red+1:ns_red+na_red)              =  gammma*S1(1:ns_red,1:na_red)
      A(ns_red+2*na_red+1:2*ns_red+2*na_red,ns_red+na_red+1:ns_red+2*na_red)     =  freq(ifreq)*S1(1:ns_red,1:na_red)
       
      KHS(1:ns_red)                                                              =  G1(1:ns_red,igd)
      KHS(ns_red+1:ns_red+na_red)                                                =  G2(1:na_red,igd)
      if (gd_complex) then
        KHS(ns_red+na_red+1:ns_red+2*na_red)                                     = -G3(1:ns_red,igd)
        KHS(ns_red+2*na_red+1:2*ns_red+2*na_red)                                 = -G4(1:ns_red,igd)
      endif       

      if (molcfg%solver%info_rsp_redspace) then
        write(molcfg%lupri,*) 'Reduced A:'
        call LS_OUTPUT(A, 1, n_red, 1, n_red, n_red, n_red, 1, molcfg%lupri)
        write(molcfg%lupri,*) 'Reduced RHS:'
        call LS_OUTPUT(KHS, 1, n_red, 1, 1, n_red, 1, 1, molcfg%lupri)
      endif
         
       !Solve set of linear equations Ax = b:
       call DGESV(n_red, 1, A, n_red, IPIV, KHS, n_red, IERR)
       !Solution vector is found in RHS.
       if (IERR /= 0) then
          WRITE(molcfg%LUPRI,'(/A, i4)') &
          &     'Problem in DGESV, IERR = ', IERR
          CALL lsQUIT(' Problem in DGESV',molcfg%lupri)
       endif
      if (molcfg%solver%info_rsp_redspace) then
        write(molcfg%lupri,*) 'Reduced solution:'
        call LS_OUTPUT(KHS, 1, n_red, 1, 1, n_red, 1, 1, molcfg%lupri)
      endif
       
       red_X(1:ns_red,ifreq,igd) = KHS(1:ns_red)
       red_Xm(1:na_red,ifreq,igd) = KHS(ns_red+1:ns_red+na_red)
       red_Xmi(1:ns_red,ifreq,igd) = KHS(ns_red+na_red+1:ns_red+2*na_red)
       red_Xi(1:ns_red,ifreq,igd) = KHS(ns_red+2*na_red+1:n_red)
    enddo
   enddo

    deallocate(KHS, A)
    deallocate(IPIV)
  end subroutine solve_complex2
  
  subroutine get_complex_res2(molcfg,F,D,S,itmic,ns_red,na_red,nfreq,ngd,gd,red_Xp_glob,&
             &red_Xpi_glob,red_Xm_glob, red_Xmi_glob,eival,xp,xm,xpi,xmi,nx_new,conv_root,conv,xsol,xsolimg,gd_complex,gdi)
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

type(rsp_molcfg), intent(inout)    :: molcfg
    integer,intent(in)               :: itmic,ns_red,ngd,na_red,nfreq
    type(Matrix),intent(in)          :: F,D,S
    type(Matrix),intent(in)          :: gd(:)
    type(Matrix),intent(in),optional :: gdi(:)
    type(Matrix), intent(inout)      :: xsol(nfreq,ngd),xsolimg(nfreq,ngd)
    logical,intent(in)               :: gd_complex
    real(realk), intent(in)          :: red_Xp_glob(:,:,:)
    real(realk), intent(inout)       :: eival(:)
    real(realk), intent(in)          :: red_Xpi_glob(:,:,:)
    real(realk), intent(in)          :: red_Xm_glob(:,:,:)
    real(realk), intent(in)          :: red_Xmi_glob(:,:,:)
    type(Matrix),intent(inout)       :: xp(nfreq,ngd),xm(nfreq,ngd),xpi(nfreq,ngd),xmi(nfreq,ngd)
    logical, intent(inout)           :: conv_root(:,:)
    integer,intent(out)              :: conv
    integer, intent(out)             :: nx_new
!local
    real(realk)                      :: gd_norm, res_norm_real, res_norm_tot, res_norm_img, gammma
    integer                          :: n_not_converged
!local
    type(Matrix)                     :: gdp,gdpi,gdm,gdmi,gT, scrR, scrI
    type(Matrix)                     :: S2Xp(nfreq,ngd), S2Xm(nfreq,ngd),S2Xpi(nfreq,ngd),S2Xmi(nfreq,ngd),xr,xi
    type(Matrix)                     :: residualsp(nfreq,ngd),residualsm(nfreq,ngd)
    type(Matrix)                     :: residualspi(nfreq,ngd),residualsmi(nfreq,ngd)
    real(realk)                      :: res_norm,x_red_norm,conv_test,ddot,diff,x_red_norm_img
    integer                          :: i,n_not_conv,ndim_red_mat,ndim,index_conv(ngd),j,k,max_red
    real(realk)                      :: a, b,c,rsp_conv_thr,rsp_thresh
    real(realk)                      :: av_norm, conv_test_img

    gammma=molcfg%solver%rsp_gamma
    max_red=molcfg%solver%rsp_maxred
    rsp_thresh= molcfg%solver%rsp_thresh

    ndim = S%nrow
    n_not_converged = 0
    conv=1
     
    do j=1,ngd
      do i=1,nfreq
        call mat_init(residualsp(i,j),ndim,ndim)
        call mat_init(residualsm(i,j),ndim,ndim)
        call mat_init(residualspi(i,j),ndim,ndim)
        call mat_init(residualsmi(i,j),ndim,ndim)
        call mat_init(S2Xp(i,j),ndim,ndim)
        call mat_init(S2Xm(i,j),ndim,ndim)
        call mat_init(S2Xpi(i,j),ndim,ndim)
        call mat_init(S2Xmi(i,j),ndim,ndim)
        call mat_zero(residualsp(i,j))
        call mat_zero(residualsm(i,j))
        call mat_zero(residualspi(i,j))
        call mat_zero(residualsmi(i,j))
        call mat_zero(S2Xp(i,j))
        call mat_zero(S2Xm(i,j))
        call mat_zero(S2Xpi(i,j))
        call mat_zero(S2Xmi(i,j))
      enddo
    enddo
!
! find E[2]*X(i)-w(i) S[2]*X(i) as linear combination of 
! previous Sigmas and RHOs. 
! X(i) = sum_j Bvecold_j REDX_j(i)
! Sigmanew(i) = E[2]X(i) = sum_j REDX_j(i) E[2]Bvecold_j =
!             = sum_j REDX_j(i) Sigmaold_j
! and add to previous -w_i S[2]*X(I) vector
! 
       rewind(lu_rho_s_rsp)
       rewind(lu_sigma_s_rsp)
       rewind(lu_sigma_a_rsp)
       rewind(lu_rho_a_rsp)
       call expand_on_basis_sym(molcfg,na_red,ndim,ngd,nfreq,red_Xm_glob,red_Xmi_glob,lu_rho_s_rsp,&
                               &S2Xm,S2Xmi)
       call expand_on_basis_sym(molcfg,ns_red,ndim,ngd,nfreq,red_Xp_glob,red_Xpi_glob,lu_rho_a_rsp,&
                               &S2Xp,S2Xpi)
      do j=1,ngd
        do i=1,nfreq
          call mat_add(-eival(i),S2Xm(i,j),gammma,S2Xmi(i,j),residualsp(i,j)) 
          call mat_add(-eival(i),S2Xp(i,j),gammma,S2Xpi(i,j),residualsm(i,j)) 
          call mat_add(eival(i),S2Xpi(i,j),gammma,S2Xp(i,j),residualsmi(i,j)) 
          call mat_add(eival(i),S2Xmi(i,j),gammma,S2Xm(i,j),residualspi(i,j))
        enddo
      enddo
      do j=1,ngd
        do i=1,nfreq
          call mat_zero(S2Xp(i,j))
          call mat_zero(S2Xm(i,j))
          call mat_zero(S2Xpi(i,j))
          call mat_zero(S2Xmi(i,j))
        enddo
      enddo
      call expand_on_basis_sym(molcfg,ns_red,ndim,ngd,nfreq,red_Xp_glob,red_Xpi_glob,lu_sigma_s_rsp,&
                              &S2Xp,S2Xpi)
      call expand_on_basis_sym(molcfg,na_red,ndim,ngd,nfreq,red_Xm_glob,red_Xmi_glob,lu_sigma_a_rsp,&
                              &S2Xm,S2Xmi)
      call mat_init(gdp,ndim,ndim)  
      call mat_init(gdm,ndim,ndim) 
      call mat_init(gT,ndim,ndim)
      if (gd_complex) then
          call mat_init(gdpi,ndim,ndim)  
          call mat_init(gdmi,ndim,ndim)
      endif
      do j=1,ngd
        ! subtract gradient
        call mat_trans(gd(j),gT)  
    
        call mat_add(0.5E0_realk,gd(j),0.5E0_realk,gT,gdp)
        call mat_add(0.5E0_realk,gd(j),-0.5E0_realk,gT,gdm)
   
        if (gd_complex) then 
          call mat_zero(gdpi)
          call mat_zero(gdmi)
   
          call mat_zero(gT)
          call mat_trans(gdi(j),gT)
          call mat_add(0.5E0_realk,gdi(j),0.5E0_realk,gT,gdpi)
          call mat_add(0.5E0_realk,gdi(j),-0.5E0_realk,gT,gdmi)
        endif
   
        do i=1,nfreq
          call mat_daxpy(1.0E0_realk,S2Xp(i,j),residualsp(i,j))
          call mat_daxpy(1.0E0_realk,S2Xm(i,j),residualsm(i,j))
          call mat_daxpy(-1.0E0_realk,S2Xmi(i,j),residualsmi(i,j))
          call mat_daxpy(-1.0E0_realk,S2Xpi(i,j),residualspi(i,j))
          call mat_daxpy(-1.0E0_realk,gdp,residualsp(i,j))
          call mat_daxpy(-1.0E0_realk,gdm,residualsm(i,j))
          if (gd_complex) then
            call mat_daxpy(1.0E0_realk,gdpi,residualspi(i,j))
            call mat_daxpy(1.0E0_realk,gdmi,residualsmi(i,j))
          endif
        enddo
      enddo
      
call mat_free(gT)
call mat_free(gdp)
call mat_free(gdm)
if (gd_complex) then
  call mat_free(gdpi)
  call mat_free(gdmi)
endif
! the residual(s) is now done: test for convergence. If not converged
! form new trial(s)  in next subroutine
! New trial(s) is equal to preconditioned residual(s)
   if (molcfg%solver%rsp_convdyn) then
          rsp_thresh= 1.0E-2_realk*molcfg%solver%rsp_dyn_thresh
!          WRITE(molcfg%LUPRI,*) 'Dynamic response convergence threshold set to', rsp_thresh
   else
          rsp_thresh= molcfg%solver%rsp_thresh
   endif
   call mat_init(scrR,ndim,ndim)
   call mat_init(scrI,ndim,ndim)
   call mat_zero(scrR)
   call mat_zero(scrI)
   do j = 1,ngd
     if (molcfg%solver%rsp_convdyn) then
       gd_norm     = sqrt(mat_sqnorm2(gd(j))+mat_sqnorm2(gdi(j)))
     endif
      do i = 1,nfreq
        call mat_add(1.0E0_realk,residualsp(i,j),1.0E0_realk,residualsm(i,j),scrR)
        call mat_add(1.0E0_realk,residualspi(i,j),1.0E0_realk,residualsmi(i,j),scrI)
        res_norm_real = mat_sqnorm2(scrR)
        res_norm_img  = mat_sqnorm2(scrI)
        res_norm_tot  = sqrt(res_norm_real*res_norm_real+res_norm_img*res_norm_img) 
        if (molcfg%solver%rsp_convdyn) then
          res_norm_tot=res_norm_tot/gd_norm
        endif
        if (res_norm_tot < rsp_thresh) then
          conv_root(i,j) = .true.
          call mat_zero(xsol(i,j))
          call mat_zero(xsolimg(i,j))
          call  expand_on_basis_full(molcfg,ns_red,na_red,ndim,red_Xp_glob(1:ns_red,i,j),red_Xm_glob(1:na_red,i,j),&
                                    &red_Xpi_glob(1:ns_red,i,j),red_Xmi_glob(1:na_red,i,j),&
                                    &lu_x_s_rsp,lu_x_a_rsp,xsol(i,j),xsolimg(i,j))
        else
          n_not_converged=n_not_converged+1
          conv=0
        endif

     ! Insert back print statements
      write (molcfg%LUPRI, '("Residual norm: rhs ", i3, " n ", i3, " is: ", E14.6, "    and frequency = " , F12.6, "   It = ", i3, &
                        & " CONV =", l2)') j, i, res_norm_tot, EIVAL(i), itmic, conv_root(i,j)
        enddo
      enddo   
      call mat_free(scrR) 
      call mat_free(scrI) 
!         conv_test = res_norm
!         conv_test_img = res_norm_img
!      if (molcfg%solver%info_rsp) then
!        write(molcfg%lupri,*) 'conv_test',conv_test
!        write(molcfg%lupri,*) 'res_norm',res_norm
!        write(molcfg%lupri,*) 'rsp_thresh',rsp_thresh
!        write(molcfg%lupri,*) 'conv_test_img',conv_test_img
!        write(molcfg%lupri,*) 'res_norm_img',res_norm_img
!      endif
!      if (conv_root .and. conv_root_img) then 
!          conv=.true.
!      endif
!   enddo  
   nx_new=0
if (n_not_converged>0) then     
      do j=1,ngd
        do i=1,nfreq
          call mat_scal(-1E0_realk,residualspi(i,j))
          call mat_scal(-1E0_realk,residualsmi(i,j))
          nx_new=nx_new+1
          call new_complex_precond(residualsp(i,j),residualsm(i,j),residualspi(i,j),residualsmi(i,j),S,-eival(j),-gammma, &
               & xp(i,j),xm(i,j),xpi(i,j),xmi(i,j),molcfg%decomp%nocc) 

          call util_scriptPx('N',D,S,xp(i,j))
          call util_scriptPx('N',D,S,xm(i,j))
          call util_scriptPx('N',D,S,xpi(i,j))
          call util_scriptPx('N',D,S,xmi(i,j))
          call mat_scal(-1E0_realk,xpi(i,j))
          call mat_scal(-1E0_realk,xmi(i,j))
        enddo
      enddo
endif
          
! Finished, deallocate local arrays
   do j=1,ngd
     do i=1,nfreq
        call mat_free(residualsp(i,j))
        call mat_free(residualsm(i,j))
        call mat_free(residualspi(i,j))
        call mat_free(residualsmi(i,j))
        call mat_free(S2Xp(i,j))
        call mat_free(S2Xm(i,j))
        call mat_free(S2Xpi(i,j))
        call mat_free(S2Xmi(i,j))
   enddo
   enddo

  end subroutine get_complex_res2
!================================================================  

  subroutine get_start_vectors2(molcfg,F,D,S,RHS_real,RHS_img,ngd,nfreq,eival,nx_new,xp,xm,xpi,xmi)
    implicit none
!================================================================ 
! A subroutine to get start vectors for damped response eqs.
! by preconditioning the RHS vectors.
! g and u components obtianed directly.
!================================================================ 
type(rsp_molcfg), intent(inout)    :: molcfg
    type(Matrix),intent(in)    :: F,D,S
    type(Matrix),intent(inout) :: RHS_real(:),RHS_img(:)
    type(Matrix),intent(inout) :: xp(:,:),xm(:,:),xpi(:,:),xmi(:,:)
    real(realk), intent(inout) :: eival(:) 
    integer,intent(inout)      :: nx_new
    integer,intent(in)         :: ngd,nfreq
    type(matrix)               :: gdp,gdm,gdpi,gdmi,gT
    integer                    :: i,j,ndim
    real(realk)                :: gammma

 ndim=D%nrow
 gammma=molcfg%solver%rsp_gamma

 
  call mat_init(gT,ndim,ndim)
  call mat_init(gdp,ndim,ndim)
  call mat_init(gdm,ndim,ndim)
  call mat_init(gdpi,ndim,ndim)
  call mat_init(gdmi,ndim,ndim)
  do j=1,ngd
     call mat_trans(RHS_real(j),gT)  
     call mat_add(0.5E0_realk,RHS_real(j),0.5E0_realk,gT,gdp)
     call mat_add(0.5E0_realk,RHS_real(j),-0.5E0_realk,gT,gdm)
     call mat_trans(RHS_img(j),gT)  
     call mat_add(0.5E0_realk,RHS_img(j),0.5E0_realk,gT,gdpi)
     call mat_add(0.5E0_realk,RHS_img(j),-0.5E0_realk,gT,gdmi)

     do i=1,nfreq  
       nx_new=nx_new+1
       call new_complex_precond(gdp,gdm,gdpi,gdmi,S,-eival(i),-gammma, &
               & xp(i,j),xm(i,j),xpi(i,j),xmi(i,j),molcfg%decomp%nocc) 

       call util_scriptPx('N',D,S,xp(i,j))
       call util_scriptPx('N',D,S,xm(i,j))
       call util_scriptPx('N',D,S,xpi(i,j))
       call util_scriptPx('N',D,S,xmi(i,j))
       call mat_scal(-1E0_realk,xpi(i,j))
       call mat_scal(-1E0_realk,xmi(i,j))
     enddo 
 enddo
 call mat_free(gdp)
 call mat_free(gdm)
 call mat_free(gdpi)
 call mat_free(gdmi)
 call mat_free(gT)
end subroutine get_start_vectors2
!================================================================  

  subroutine orthonormalize_sym(molcfg,D,S,br,bi,Nb_new,Nb_prev,bvecs,lu)
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
    type(rsp_molcfg),intent(inout) :: molcfg
    type(Matrix),intent(in)    :: D,S
    integer, intent(in)        :: Nb_prev
    integer, intent(inout)     :: Nb_new
    type(Matrix),intent(inout) :: br(Nb_new), bi(Nb_new)
    type(Matrix),intent(inout) :: bvecs(:)
    integer, intent(in)        :: lu
!local 
    integer                    :: irx,i,j,iturn,ndim,k,ibvec,jbvec,jrx,no_of_new,lub_rsp_vec,dummy_i,max_red
    logical,allocatable        :: lin_depend(:) !linear dependency index array
    type(matrix)               :: B_scr, b_k, Xf, rho_k
    type(matrix)               :: orthovec !Will properly be changed
    real(realk)                :: TT,T1,T2,dummy_real,pr1,pr2,thr_ld
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
    thr_ld=1E-15_realk
    allocate(lin_depend(2*Nb_new))
    do i=1,2*Nb_new
      lin_depend(i)=.false.
    enddo    
    call mat_init(b_k,ndim,ndim)
    call mat_zero(b_k)
    run_ortho_again = .false.
    do iturn = 1, 2 
    if (Nb_prev > 0) then
    !
    ! Orthogonalize new b-vectors agains previous (old) b-vectors
    ! (both (Z, Y) and (Y, Z))
    !
       if (molcfg%solver%info_rsp) then
         write(molcfg%lupri,*) 'Orthogonalize new b-vectors agains previous b'
         write(molcfg%lupri,*) 'NB_PREV, NB_NEW ', NB_PREV, NB_NEW
       endif 
   
       rewind(lu)
       do k = 1,Nb_prev
          call mat_read_from_disk(lu,b_k,RSPOnMaster)
          do irx = 1,Nb_new
             !if lin_depend(irx) == 0, the vector is skipped
             !because of linear dependencies
             if (.not.lin_depend(irx)) then
                 TT = mat_dotproduct(b_k,br(irx))
                 call mat_daxpy(-TT,b_k,br(irx))
             endif
             if (.not.lin_depend(Nb_new+irx)) then
                 TT = mat_dotproduct(b_k,bi(irx))
                 call mat_daxpy(-TT,b_k,bi(irx))
             endif
          enddo
       enddo
    do irx=1,Nb_new
      if (.not.lin_depend(irx)) then
        pr1=sqrt(mat_dotproduct(br(irx),br(irx)))
        if (pr1 < thr_ld) then
           lin_depend(irx)=.true.
        else
           pr2=1E0_realk/pr1
           call mat_scal(pr2,br(irx)) 
        endif
      endif
      if (.not.lin_depend(Nb_new+irx)) then
        pr1=sqrt(mat_dotproduct(bi(irx),bi(irx)))
        if (pr1 < thr_ld) then
           lin_depend(Nb_new+irx)=.true.
        else
           pr2=1E0_realk/pr1
           call mat_scal(pr2,bi(irx)) 
        endif
      endif
    enddo
    endif
      
    if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Orthogonalize new b-vectors against each other '
   !
   ! Orthogonalize new vectors against each other
   !
    do ibvec = 1,Nb_new !index for current bvector 
!       if (molcfg%solver%info_rsp) write(molcfg%lupri,*)'Current ibvec index is ', ibvec, lin_depend(ibvec), lin_depend(Nb_new+ibvec)
       if (.not.lin_depend(ibvec)) then
         do jrx = 1,(ibvec-1)
           if (.not.lin_depend(jrx)) then
             T1 = mat_sqnorm2(br(jrx))
!               if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Norm of jbvec', jbvec, T1
             if (T1 > thr_ld) then 
               T2 = mat_dotproduct(br(jrx),br(ibvec))
!               if (molcfg%solver%info_rsp) write(molcfg%lupri,*) '<bjbvec|bibvec>', T2
               TT = -T2/T1
!               if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'TT = -T2/T1', TT
               call mat_daxpy(TT,br(jrx),br(ibvec))
             endif
           endif
         enddo
       endif
    enddo
    do ibvec = 1,Nb_new !index for current bvector 
!       if (molcfg%solver%info_rsp) write(molcfg%lupri,*)'Current ibvec index is ', ibvec, lin_depend(ibvec), lin_depend(Nb_new+ibvec)
       if (.not.lin_depend(Nb_new+ibvec)) then
         do jrx = 1,Nb_new
           if (.not.lin_depend(jrx)) then
             T1 = mat_sqnorm2(br(jrx))
!               if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Norm of jbvec', jbvec, T1
             if (T1 > thr_ld) then 
               T2 = mat_dotproduct(br(jrx),bi(ibvec))
!               if (molcfg%solver%info_rsp) write(molcfg%lupri,*) '<bjbvec|bibvec>', T2
               TT = -T2/T1
!               if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'TT = -T2/T1', TT
               call mat_daxpy(TT,br(jrx),bi(ibvec))
             endif
           endif
         enddo
         do jrx = 1,(ibvec-1)
           if (.not.lin_depend(Nb_new+jrx)) then
             T1 = mat_sqnorm2(bi(jrx))
!               if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'Norm of jbvec', jbvec, T1
             if (T1 > thr_ld) then 
               T2 = mat_dotproduct(bi(jrx),bi(ibvec))
!               if (molcfg%solver%info_rsp) write(molcfg%lupri,*) '<bjbvec|bibvec>', T2
               TT = -T2/T1
!               if (molcfg%solver%info_rsp) write(molcfg%lupri,*) 'TT = -T2/T1', TT
               call mat_daxpy(TT,bi(jrx),bi(ibvec))
             endif
           endif
         enddo
       endif
     enddo
    do irx=1,Nb_new
      if (.not.lin_depend(irx)) then
        pr1=sqrt(mat_dotproduct(br(irx),br(irx)))
        if (pr1 < thr_ld) then
          if (iturn == 1) then
            lin_depend(irx)=.true.
          else
            stop 'error in R ortho'
          endif  
        else
           pr2=1E0_realk/pr1
           call mat_scal(pr2,br(irx)) 
        endif
      endif
      if (.not.lin_depend(Nb_new+irx)) then
        pr1=sqrt(mat_dotproduct(bi(irx),bi(irx)))
        if (pr1 < thr_ld) then
          if (iturn == 1) then
            lin_depend(Nb_new+irx)=.true.
          else
            stop 'error in I ortho'
          endif  
        else
           pr2=1E0_realk/pr1
           call mat_scal(pr2,bi(irx)) 
        endif
      endif
    enddo
 enddo

!             WRITE(molcfg%lupri,'(/A)') &
!             &     'Error: Already ran twice through orthonormalize!'
!             CALL lsQUIT('Error: Already ran twice through orthonormalize!',molcfg%lupri)


! Add new vectors
!    
    no_of_new = 0
    do irx = 1,Nb_new
       if (.not.lin_depend(irx)) then
         no_of_new = no_of_new + 1
         call mat_assign(bvecs(no_of_new),br(irx))
       endif
    enddo
    do irx = 1,Nb_new
       if (.not.lin_depend(Nb_new+irx)) then
         no_of_new = no_of_new + 1
         call mat_assign(bvecs(no_of_new),bi(irx))
       endif
    enddo
!
!     Set NB_NEW to actual number of acceptable new trial vectors
!
     
     Nb_new = no_of_new
     
!     do i=1,Nb_new
!       do j=1,Nb_new
!        pr1=mat_dotproduct(bvecs(i),bvecs(j))
!        write(molcfg%lupri,*) 'ortho czesc',i,j,pr1
!       enddo
!     enddo
     
    DEALLOCATE(lin_depend)
    call mat_free(b_k)

  end subroutine orthonormalize_sym
  subroutine expand_on_basis_sym(molcfg,n_red,ndim,ngd,nfreq,red_X,red_Xi,lu_basis,x,xi)
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
    type(rsp_molcfg), intent(inout)    :: molcfg
    integer, intent(in) :: n_red, lu_basis, ndim,ngd,nfreq
    real(realk), intent(in) :: red_X(:,:,:),red_Xi(:,:,:)
    type(Matrix), intent(inout) :: x(:,:),xi(:,:)  

    type(matrix) :: b_k
    integer :: k,i,j
    logical :: OnMaster
    OnMaster = .TRUE.

    call mat_init(b_k,ndim,ndim)

    rewind(lu_basis)
    do k = 1,n_red
      call mat_read_from_disk(lu_basis,b_k,OnMaster)
       do j=1,ngd
          do i=1,nfreq
            call mat_daxpy(red_X(k,i,j),b_k,x(i,j))
            call mat_daxpy(red_Xi(k,i,j),b_k,xi(i,j))
          enddo
       enddo
    enddo
    call mat_free(b_k)

  end subroutine expand_on_basis_sym
  subroutine expand_on_basis_full(molcfg,ns_red,na_red,ndim,red_Xp,red_Xm,red_Xpi,red_Xmi,lu_s,lu_a,x,xi)
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
    type(rsp_molcfg), intent(inout)    :: molcfg
    integer, intent(in) :: ns_red, na_red, lu_s, lu_a, ndim
    real(realk), intent(in) :: red_Xp(:),red_Xpi(:),red_Xm(:),red_Xmi(:)
    type(Matrix), intent(inout) :: x,xi  

    type(matrix) :: b_k
    integer :: k,i,j
    logical :: OnMaster
    OnMaster = .TRUE.

    call mat_init(b_k,ndim,ndim)

    rewind(lu_s)
    do k = 1,ns_red
      call mat_read_from_disk(lu_s,b_k,OnMaster)
            call mat_daxpy(red_Xp(k),b_k,x)
            call mat_daxpy(red_Xpi(k),b_k,xi)
    enddo
    rewind(lu_a)
    do k = 1,na_red
      call mat_read_from_disk(lu_a,b_k,OnMaster)
            call mat_daxpy(red_Xm(k),b_k,x)
            call mat_daxpy(red_Xmi(k),b_k,xi)
    enddo
    call mat_free(b_k)

  end subroutine expand_on_basis_full

 end module COMPLEXNEWSYMSOLVER
