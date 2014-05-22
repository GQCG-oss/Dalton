!@file: this file contains the pno-ccsd residual routine and all the relevant
!PNO-specific routines
!author: Patrick Ettenhuber

module pno_ccsd_module

  use precision
  use typedef
  use typedeftype
  use dec_typedef_module
  use screen_mod
  use tensor_interface_module
  use IntegralInterfaceDEC
  
  
  use dec_fragment_utils
  
  public :: get_ccsd_residual_pno_style, &
     & get_pno_trafo_matrices,      & 
     & get_pno_overlap_matrices,    &
     & successive_4ao_mo_trafo, free_PNOSpaceInfo
  
  private
  
  contains
  
  !> \author Patrick Ettenhuber
  !> \date November 2013
  !> \brief this subroutine calculates the ccsd residual by transforming each
  !>        doubles amplitudes to their respective set of PNO's and then
  !>        transforming the result vector back to the reference basis. 
  subroutine get_ccsd_residual_pno_style(t1,t2,o1,o2,no,nv,nb,xo,xv,yo,yv,&
        &mylsitem,fj,pno_cv,pno_S,nspaces,oof,vvf,ifo,iter,f)
     implicit none
     !ARGUMENTS
     integer, intent(in) :: no, nv, nb,iter,nspaces
     logical, intent(in) :: fj
     real(realk), intent(inout) :: t1(nv,no), t2(nv,no,nv,no)
     real(realk), intent(inout) :: o1(nv,no), o2(nv,no,nv,no)
     real(realk), intent(in) :: xo(nb,no), xv(nb,nv), yo(nb,no), yv(nb,nv),ifo(nb,nb)
     type(lsitem), intent(inout) :: mylsitem
     real(realk), intent(inout) :: oof(no,no),vvf(nv,nv)
     type(decfrag),intent(in),optional :: f
     type(PNOSpaceInfo),intent(inout) :: pno_cv(nspaces)
     type(PNOSpaceInfo),intent(inout) :: pno_S(nspaces*(nspaces-1)/2)
     !INTERNAL VARIABLES
     type(array),pointer :: pno_o2(:),pno_t2(:),pno_gvvvv(:),pno_govov(:),pno_gvovo(:)
     integer :: ns,c,nc,nc2
     real(realk),pointer :: w1(:),w2(:),w3(:), w4(:),w5(:)
     real(realk),pointer :: p1(:,:,:,:), p2(:,:,:,:), p3(:,:,:,:), p4(:,:,:,:),h1(:), h2(:), r1(:,:),r2(:,:)
     real(realk),pointer :: gvvvv(:), gvovo(:), govov(:), goooo(:), goovv(:), gvvov(:), gooov(:)
     integer :: i, j, a, b, i_idx, la, lg, fa, fg
     integer(kind=8) :: o2v2
     character(ARR_MSG_LEN) :: msg
     real(realk),pointer :: d(:,:), d1(:,:), d2(:,:),t(:), t22(:), t21(:), o(:),vof(:),ovf(:)
     real(realk),pointer :: Lvoov(:)  
     real(realk) :: nnorm, norm 
     integer, pointer :: idx(:),idx1(:),idx2(:), p_idx(:,:), p_nidx(:), oidx1(:,:),oidx2(:,:)
     integer, pointer :: s_idx(:,:,:), s_nidx(:)
     integer :: pno,pno1,pno2,pnv,pnv1,pnv2, k, l, nidx1, nidx2, spacemax
     logical :: skiptrafo,skiptrafo2,save_gvvvv_is,cyc,use_triangular,PS
     real(realk), pointer :: iFock(:,:), Dens(:,:)
     integer(kind=8) :: maxsize, myload
     integer :: pair,paircontribs,paircontrib(2,2),rpd
     integer :: order1(4)
     integer :: suborder(2)
     logical :: master
     !Integral stuff
     integer :: alphaB,gammaB,dimAlpha,dimGamma
     integer :: dim1,dim2,dim3,MinAObatch
     integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd
     integer :: iorb,nthreads
     type(batchtoorb), pointer :: batch2orbAlpha(:)
     type(batchtoorb), pointer :: batch2orbGamma(:)
     Character(80)        :: FilenameCS,FilenamePS
     Character(80),pointer:: BatchfilenamesCS(:,:)
     Character(80),pointer:: BatchfilenamesPS(:,:)
     Character            :: INTSPEC(5)
     logical :: FoundInMem,fullRHS, doscreen
     integer :: MaxAllowedDimAlpha,MaxActualDimAlpha,nbatchesAlpha
     integer :: MaxAllowedDimGamma,MaxActualDimGamma,nbatchesGamma
     integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
     integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
     TYPE(DECscreenITEM)    :: DecScreen
     real(realk) :: MemFree
     real(realk) :: ref(no*nv*nv*no), ref1(no*nv), u(nv,no,nv,no)
     real(realk), parameter :: p20 = 2.0E0_realk
     real(realk), parameter :: p10 = 1.0E0_realk
     real(realk), parameter :: m10 = -1.0E0_realk
     real(realk), parameter :: m05 = -0.5E0_realk
     real(realk), parameter :: p05 = 0.5E0_realk
     real(realk), parameter :: nul = 0.0E0_realk
  
  
     o2v2 = (i8*no**2)*nv**2
     use_triangular = .true.
     master         = .true.

     paircontribs = 2
     paircontrib(1:2,1) = [1,2]
     paircontrib(1:2,2) = [2,1]

     if(fj.and..not.present(f))call lsquit("ERROR(get_ccsd_residual_pno_style):wrong input fj without f",-1)


     call mem_alloc( pno_t2, nspaces  )
     call mem_alloc( pno_o2, nspaces  )
     call mem_alloc( gvvvv,  nv**4    )
     call mem_alloc( gvovo,  o2v2     )
     call mem_alloc( govov,  o2v2     )
     call mem_alloc( goooo,  no**4    )
     call mem_alloc( goovv,  o2v2     )
     call mem_alloc( Lvoov,  o2v2     )
     call mem_alloc( gvvov,  nv**3*no )
     call mem_alloc( gooov,  no**3*nv )
     call mem_alloc( p_nidx, nspaces  )
     call mem_alloc( p_idx,  nspaces  , nspaces )
     call mem_alloc( s_nidx, no )
     call mem_alloc( s_idx,  2       , nspaces , no )
     call mem_alloc( ovf,    no*nv    )
     call mem_alloc( vof,    nv*no    )

     maxsize=nb**4
     call mem_alloc( w1, maxsize )
     call mem_alloc( w2, nb**3*max(nv,no) )
     call mem_alloc( w3, maxsize )

     ! if  we have a fragment job the basis of the atomic (pair) site has to be
     ! treated in a special way, either LO or FNO basis, this will be element 1
     ! in all the array arrays
     spacemax = 2
     if(present(f))spacemax = max(f%noccEOS,spacemax)
     spacemax = spacemax + 2 + 2

     !Get pair interaction space information
     call get_pair_space_info(pno_cv,p_idx,p_nidx,s_idx,s_nidx,nspaces,no)

     !Get all the pno amplitudes with index restrictions i<=j
     call get_pno_amplitudes(t2,pno_cv,pno_t2,nspaces,no,nv)

     !initialize the pno_residual according to the allocated pno_cv
     call init_pno_residual(pno_cv,pno_o2,nspaces)



     !call II_get_AbsoluteValueOcc_overlap(DECinfo%output,DECinfo%output,setting,nb,no,out)

     !INTEGRAL DIRECT STUFF HAPPENINING HERE
     ! Set integral info
     ! *****************
     INTSPEC(1)               = 'R' !R = Regular Basis set on the 1th center 
     INTSPEC(2)               = 'R' !R = Regular Basis set on the 2th center 
     INTSPEC(3)               = 'R' !R = Regular Basis set on the 3th center 
     INTSPEC(4)               = 'R' !R = Regular Basis set on the 4th center 
     INTSPEC(5)               = 'C' !C = Coulomb operator
     doscreen                 = MyLsItem%setting%scheme%cs_screen.OR.MyLsItem%setting%scheme%ps_screen

     !==================================================
     !                  Batch construction             !
     !==================================================


     ! Get free memory and determine maximum batch sizes
     ! -------------------------------------------------
     if(master)then
        call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch,'R')
        call get_currently_available_memory(MemFree)
        !call get_max_batch_sizes(scheme,nb,nv,no,MaxAllowedDimAlpha,MaxAllowedDimGamma,&
        !   &MinAObatch,DECinfo%manual_batchsizes,iter,MemFree,.true.,els2add,local)
        MaxAllowedDimAlpha = nb
        MaxAllowedDimGamma = nb
     endif
     ! ************************************************
     ! * Determine batch information for Gamma batch  *
     ! ************************************************

     ! Orbital to batch information
     ! ----------------------------
     call mem_alloc(orb2batchGamma,nb)
     call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma,&
        & nb,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,&
        &nbatchesGamma,orb2BatchGamma,'R')
     if(master.and.DECinfo%PL>1)&
        &write(DECinfo%output,*) 'BATCH: Number of Gamma batches   = ', nbatchesGamma,&
        & 'with maximum size',MaxActualDimGamma

     ! Translate batchindex to orbital index
     ! -------------------------------------
     call mem_alloc(batch2orbGamma,nbatchesGamma)
     do i=1,nbatchesGamma
        call mem_alloc(batch2orbGamma(i)%orbindex,batchdimGamma(i))
        batch2orbGamma(i)%orbindex = 0
        batch2orbGamma(i)%norbindex = 0
     end do
     do iorb=1,nb
        i = orb2batchGamma(iorb)
        batch2orbGamma(i)%norbindex = batch2orbGamma(i)%norbindex+1
        K = batch2orbGamma(i)%norbindex
        batch2orbGamma(i)%orbindex(K) = iorb
     end do


     ! ************************************************
     ! * Determine batch information for Alpha batch  *
     ! ************************************************

     ! Orbital to batch information
     ! ----------------------------
     call mem_alloc(orb2batchAlpha,nb)
     call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha,&
        & nb,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha,'R')
     if(master.and.DECinfo%PL>1)&
        &write(DECinfo%output,*) 'BATCH: Number of Alpha batches   = ', nbatchesAlpha&
        &, 'with maximum size',MaxActualDimAlpha

     !Translate batchindex to orbital index
     !-------------------------------------
     call mem_alloc(batch2orbAlpha,nbatchesAlpha)
     do i=1,nbatchesAlpha
        call mem_alloc(batch2orbAlpha(i)%orbindex,batchdimAlpha(i) )
        batch2orbAlpha(i)%orbindex = 0
        batch2orbAlpha(i)%norbindex = 0
     end do
     do iorb=1,nb
        i = orb2batchAlpha(iorb)
        batch2orbAlpha(i)%norbindex = batch2orbAlpha(i)%norbindex+1
        K = batch2orbAlpha(i)%norbindex
        batch2orbAlpha(i)%orbindex(K) = iorb
     end do

     ! ************************************************
     ! *  Allocate matrices used in the batched loop  *
     ! ************************************************

     ! ************************************************
     ! *  precalculate the full schreening matrix     *
     ! ************************************************

     ! This subroutine builds the full screening matrix.
     call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting,&
        & nbatchesAlpha,nbatchesGamma,INTSPEC)
     if(mylsitem%setting%scheme%cs_screen.or.mylsitem%setting%scheme%ps_screen)then
        call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
           & nb,nbatchesAlpha,nbatchesGamma,&
           & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
           & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
        call II_getBatchOrbitalScreenK(DecScreen,mylsitem%setting,&
           & nb,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
           & batchindexAlpha,batchindexGamma,&
           & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
     endif

     myload = 0

     fullRHS = nbatchesGamma.EQ.1.AND.nbatchesAlpha.EQ.1

     BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
        dimGamma   = batchdimGamma(gammaB)                         ! Dimension of gamma batch
        GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
        GammaEnd   = batch2orbGamma(gammaB)%orbindex(dimGamma)     ! Last index in gamma batch
        !short hand notation
        fg         = GammaStart
        lg         = dimGamma


        BatchAlpha: do alphaB = 1, nbatchesAlpha

           !check if the current job is to be done by current node
           !call check_job(scheme,first_round,dynamic_load,alphaB,gammaB,nbatchesAlpha,&
           !  &nbatchesGamma,tasks,tasksw,print_debug)
           !break the loop if alpha become too large, necessary to account for all
           !of the mpi and non mpi schemes, this is accounted for, because static,
           !and dynamic load balancing are enabled
           if(alphaB>nbatchesAlpha) exit

           dimAlpha   = batchdimAlpha(alphaB)                              ! Dimension of alpha batch
           AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
           AlphaEnd   = batch2orbAlpha(alphaB)%orbindex(dimAlpha)          ! Last index in alpha batch

           !short hand notation
           fa         = AlphaStart
           la         = dimAlpha
           myload     = myload + la * lg

           IF(doscreen)Mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGabKLHS(alphaB)%p
           IF(doscreen)Mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGabKRHS(gammaB)%p

           call II_GET_DECPACKED4CENTER_K_ERI(DECinfo%output,DECinfo%output, &
              & Mylsitem%setting,w1,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
              & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),dimAlpha,nb,dimGamma,nb,INTSPEC,fullRHS)

           w3 = w1
           !gvvvv
           call successive_4ao_mo_trafo(nb,w1,xv,nv,yv,nv,xv,nv,yv,nv,w2)
           call dcopy(nv**4,w1,1,gvvvv,1)
           !goooo
           w1 = w3
           call successive_4ao_mo_trafo(nb,w1,xo,no,yo,no,xo,no,yo,no,w2)
           call dcopy(no**4,w1,1,goooo,1)
           !govov
           w1 = w3
           call successive_4ao_mo_trafo(nb,w1,xo,no,yv,nv,xo,no,yv,nv,w2)
           call dcopy(nv**2*no**2,w1,1,govov,1)
           !goovv
           w1 = w3
           call successive_4ao_mo_trafo(nb,w1,xo,no,yo,no,xv,nv,yv,nv,w2)
           call dcopy(nv**2*no**2,w1,1,goovv,1)
           !Lvoov = 2gvoov - gvvoo
           w1 = w3
           call successive_4ao_mo_trafo(nb,w1,xv,nv,yo,no,xo,no,yv,nv,w2)
           call array_reorder_4d( p20, w1, nv, no ,no, nv, [1,2,3,4], nul, Lvoov)
           call array_reorder_4d( m10, goovv, no, no ,nv, nv, [3,2,1,4], p10, Lvoov)
           !gvvov
           w1 = w3
           call successive_4ao_mo_trafo(nb,w1,xv,nv,yv,nv,xo,no,yv,nv,w2)
           call dcopy(nv**3*no,w1,1,gvvov,1)
           !gooov
           w1 = w3
           call successive_4ao_mo_trafo(nb,w1,xo,no,yo,no,xo,no,yv,nv,w2)
           call dcopy(nv*no**3,w1,1,gooov,1)
           !gvovo
           w1 = w3
           call successive_4ao_mo_trafo(nb,w1,xv,nv,yo,no,xv,nv,yo,no,w2)
           call dcopy(nv**2*no**2,w1,1,gvovo,1)

           if(fa<=fg+lg-1)then
              do ns = 1, nspaces
                 !  no_pair
                 !  nv_pair
                 !  xo_pair
                 !  xv_pair
                 !  yo_pair
                 !  yv_pair
                 !  call get_a22_and_prepb22_terms_ex(w0,w1,w2,w3,tpl(ns),tmi(ns),no_pair,nv_pair,nb,fa,fg,la,lg,&
                 !    &xo_pair,yo_pair,xv_pair,yv_pair,omega2(ns),sio4(ns),4,[w0%n,w1%n,w2%n,w3%n],.false.,.false.,.false.)
              enddo
           endif
        enddo BatchAlpha
     enddo BatchGamma
     

     ! Free integral stuff
     ! *******************
     nullify(Mylsitem%setting%LST_GAB_LHS)
     nullify(Mylsitem%setting%LST_GAB_RHS)
     call free_decscreen(DECSCREEN)

     ! Free gamma stuff
     call mem_dealloc(orb2batchGamma)
     call mem_dealloc(batchdimGamma)
     call mem_dealloc(batchsizeGamma)
     call mem_dealloc(batchindexGamma)
     do i=1,nbatchesGamma
        call mem_dealloc(batch2orbGamma(i)%orbindex)
        batch2orbGamma(i)%orbindex => null()
     end do
     call mem_dealloc(batch2orbGamma)

     ! Free alpha stuff
     call mem_dealloc(orb2batchAlpha)
     call mem_dealloc(batchdimAlpha)
     call mem_dealloc(batchsizeAlpha)
     call mem_dealloc(batchindexAlpha)
     do i=1,nbatchesAlpha
        call mem_dealloc(batch2orbAlpha(i)%orbindex)
        batch2orbAlpha(i)%orbindex => null()
     end do
     call mem_dealloc(batch2orbAlpha)

     call mem_dealloc( w1 )
     call mem_dealloc( w2 )
     call mem_dealloc( w3 )


     maxsize=max(max(i8*no,i8*nv)**4,i8*nb*max(nv,nb))
     call mem_alloc( w1, maxsize )
     call mem_alloc( w2, maxsize )
     call mem_alloc( w3, maxsize )
     call mem_alloc( w4, maxsize )
     call mem_alloc( w5, maxsize )

     !!!!!!!!!!!!!!!!!!!
     !GET FOCK MATRICES!
     !!!!!!!!!!!!!!!!!!!

     !allocate the density matrix
     call mem_alloc(iFock,nb,nb)
     call mem_alloc(Dens,nb,nb)
     !calculate inactive fock matrix in ao basis
     call dgemm('n','t',nb,nb,no,1.0E0_realk,yo,nb,xo,nb,0.0E0_realk,Dens,nb)
     iFock = 0.0E0_realk
     call II_get_fock_mat_full(DECinfo%output,DECinfo%output,MyLsItem%setting,nb,&
        & Dens,.false.,iFock)
     !use dens as temporay array 
     call ii_get_h1_mixed_full(DECinfo%output,DECinfo%output,MyLsItem%setting,&
        & Dens,nb,nb,AORdefault,AORdefault)
     ! Add one- and two-electron contributions to Fock matrix
     call daxpy(nb**2,1.0E0_realk,Dens,1,iFock,1)
     call daxpy(nb**2,1.0E0_realk,ifo,1,iFock,1)
     !Free the density matrix
     call mem_dealloc(Dens)
     !Transform inactive Fock matrix into the different mo subspaces
     ! -> Foo
     call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,iFock,nb,0.0E0_realk,w1,no)
     call dgemm('n','n',no,no,nb,1.0E0_realk,w1,no,yo,nb,0.0E0_realk,oof,no)
     ! -> Fov
     call dgemm('n','n',no,nv,nb,1.0E0_realk,w1,no,yv,nb,0.0E0_realk,ovf,no)
     ! -> Fvo
     call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,iFock,nb,0.0E0_realk,w1,nv)
     call dgemm('n','n',nv,no,nb,1.0E0_realk,w1,nv,yo,nb,0.0E0_realk,vof,nv)
     ! -> Fvv
     call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1,nv,yv,nb,0.0E0_realk,vvf,nv)
     call mem_dealloc(iFock)


     !!DEBUG: A2 term
     !!**************
     !u = p20*t2
     !call array_reorder_4d( m10, t2,   nv, no, nv, no, [1,4,3,2], p10, u  )
     !ref = gvovo

     !!A2.2 contribution
     !call array_reorder_4d( p10, gvvvv, nv, nv, nv, nv, [1,3,2,4], nul, w1  )
     !call array_reorder_4d( p10, t2,    nv, no, nv, no, [1,3,2,4], nul, w2  )
     !call dgemm( 'n', 'n', nv**2, no**2, nv**2, p10, w1, nv**2, w2, nv**2, nul, w3, nv**2 )
     !call array_reorder_4d( p10, w3,    nv, nv, no, no, [1,3,2,4], p10, ref )

     !call print_norm(w3,o2v2,nnorm,.true.)
     !call print_norm(ref,o2v2,norm,.true.)
     !write(*,*)' DEBUG A2/TOT:',sqrt(nnorm),sqrt(norm)

     !!DEBUG: B2 term
     !!**************
     !call array_reorder_4d( p10, govov, no, nv, no, nv, [2,4,1,3], nul, w1  ) ! kcld -> cdkl
     !call array_reorder_4d( p10, t2,    nv, no, nv, no, [2,4,1,3], nul, w2  ) ! cidj -> ijcd
     !call array_reorder_4d( p10, goooo, no, no, no, no, [2,4,1,3], nul, w3  ) ! kilj -> ijkl
     !! g(ijkl) + t(ijcd) g(cdkl) = B(klij)
     !call dgemm( 'n', 'n', no**2, no**2, nv**2, p10, w2, no**2, w1, nv**2, p10, w3, no**2)
     !! B(ijkl) t(klab) = B2(ijab) 
     !call dgemm( 'n', 'n', no**2, nv**2, no**2, p10, w3, no**2, w2, no**2, nul, w1, no**2)
     !call array_reorder_4d( p10, w1, no, no, nv, nv, [3,1,4,2], p10, ref )

     !call print_norm(w1,o2v2,nnorm,.true.)
     !call print_norm(ref,o2v2,norm,.true.)
     !write (*,*)' DEBUG B2/TOT:',sqrt(nnorm),sqrt(norm)


     !!DEBUG: C2 term
     !!**************
     !ref = 0.0E0_realk
     !call array_reorder_4d( p10, goovv, no, no, nv, nv, [4,1,2,3], nul, w1 ) ! kjac -> ckja
     !call array_reorder_4d( p10, t2,    nv, no, nv, no, [2,3,4,1], nul, w2 ) ! aldj -> ldja
     !call array_reorder_4d( p10, govov, no, nv, no, nv, [4,1,3,2], nul, w3 ) ! kdlc -> ckld
     !!C intermediate w1(ckja) -  0.5 w3(ckld) w2(ldja)
     !call dgemm( 'n', 'n', no*nv, no*nv, no*nv, m05, w3, no*nv, w2, no*nv, p10, w1, no*nv)
     !call array_reorder_4d( p10, t2,    nv, no, nv, no, [1,4,3,2], nul, w2 ) ! bkci -> bick
     !!call print_tensor_unfolding_with_labels(w2,&
     !!  &[nv,no],'bi',2,[nv,no],'ck',2,'FULL T2')
     !!call print_tensor_unfolding_with_labels(w1,&
     !!  &[nv,no],'ck',2,[no,nv],'ja',2,'FULL Goovv')
     !call dgemm( 'n', 'n', no*nv, no*nv, no*nv, m10, w2, no*nv, w1, no*nv, nul, w3, no*nv)

     !!call print_tensor_unfolding_with_labels(w3,&
     !!  &[nv,no],'bi',2,[no,nv],'ja',2,'FULL ERG GEMM')
     !!USE THE SYMMETRIZED CONTRIBUTION, i.e. P_{ij}^{ab} (1+0.5P_{ij}) * w3
     !call array_reorder_4d( p10, w3, nv, no, no, nv, [4,2,1,3], nul, w2 ) ! bija -> aibj
     !call array_reorder_4d( p05, w3, nv, no, no, nv, [4,3,1,2], p10, w2 ) ! bjia -> aibj
     !call array_reorder_4d( p10, w3, nv, no, no, nv, [1,3,4,2], p10, w2 ) ! ajib -> aibj
     !call array_reorder_4d( p05, w3, nv, no, no, nv, [1,2,4,3], p10, w2 ) ! aijb -> aibj

     !ref = ref + w2(1:o2v2)

     !!call print_tensor_unfolding_with_labels(ref,&
     !!  &[nv,no],'ai',2,[nv,no],'bj',2,'FULL ERG C')
     !call print_norm(w2,o2v2,nnorm,.true.)
     !call print_norm(ref,o2v2,norm,.true.)
     !write (*,*)' DEBUG C2/TOT:',sqrt(nnorm),sqrt(norm)

     !!DEBUG: D2 term
     !!**************
     !call array_reorder_4d( p10, Lvoov, nv, no, no, nv, [4,3,1,2], nul, w1 ) ! aikc -> ckai
     !call array_reorder_4d( p10, u,     nv, no, nv, no, [4,3,1,2], nul, w2 ) ! aidl -> ldai
     !call array_reorder_4d( p20, govov, no, nv, no, nv, [4,3,1,2], nul, w3 ) ! ldkc -> ckld
     !call array_reorder_4d( m10, govov, no, nv, no, nv, [4,1,3,2], p10, w3 ) ! ldkc -> clkd
     !call dgemm('n','n',nv*no,nv*no,no*nv, p05, w3,nv*no,w2,no*nv, p10, w1, nv*no)
     !call dgemm('n','n',nv*no,nv*no,nv*no, p05, u ,nv*no,w1,nv*no, nul, w2, nv*no)
     !w3(1:o2v2) = w2(1:o2v2)
     !call array_reorder_4d( p10, w2, nv, no, nv, no, [3,4,1,2], p10, w3 )

     !!ref = ref + w3(1:o2v2)

     !call print_norm(w3,o2v2,nnorm,.true.)
     !call print_norm(ref,o2v2,norm,.true.)
     !!write (*,*)' DEBUG D2/TOT:',sqrt(nnorm),sqrt(norm)


     !!DEBUG: E2 term
     !!**************
     !!part 1 vv
     !call ass_D2to1(vvf,h1,[nv,nv])
     !w1(1:nv**2) = h1(1:nv**2)
     !h1 => null()
     !call array_reorder_4d( p10, govov, no, nv, no, nv, [3,2,1,4], nul, w3 )
     !call dgemm('n','n',nv,nv,no*nv*no, m10, u,nv,w3,no*nv*no, p10, w1, nv)
     !call array_reorder_4d( p10, t2, nv, no, nv, no, [3,4,1,2], nul, w2)
     !call dgemm( 'n','n',nv,no*nv*no,nv,p10,w1,nv,w2,nv,nul,w3,nv)
     !w2(1:o2v2) = w3(1:o2v2)
     !call array_reorder_4d( p10, w3, nv, no, nv, no, [3,4,1,2], p10, w2)

     !!ref = ref + w2(1:o2v2)

     !call print_norm(w2,o2v2,nnorm,.true.)
     !call print_norm(ref,o2v2,norm,.true.)
     !write (*,*)' DEBUG E21/TOT:',sqrt(nnorm),sqrt(norm)

     !!part 2
     !call ass_D2to1(oof,h1,[no,no])
     !w1(1:no**2) = h1(1:no**2)
     !h1 => null()
     !call array_reorder_4d( p10, govov, no, nv, no, nv, [1,4,3,2], nul, w2)
     !call dgemm('n','n',no,no,nv*no*nv, p10, w2,no,u,nv*no*nv, p10, w1, no)
     !call dgemm('n','n',nv*no*nv,no,no, m10, t2, nv*no*nv,w1,no, nul, w2, nv*no*nv)
     !w3(1:o2v2) = w2(1:o2v2)
     !call array_reorder_4d( p10, w3, nv, no, nv, no, [3,4,1,2], p10, w2)

     !!ref = ref + w2(1:o2v2)

     !call print_norm(w2,o2v2,nnorm,.true.)
     !call print_norm(ref,o2v2,norm,.true.)
     !!write (*,*)' DEBUG E22/TOT:',sqrt(nnorm),sqrt(norm)


     !ref1 = vof

     !!DEBUG SINGLES A1
     !!****************
     !call array_reorder_4d( p10, u, nv, no, nv, no, [3,2,1,4], nul, w1) ! ckdi -> dkci
     !call dgemm( 'n','n',nv, no, nv**2*no, p10,gvvov,nv,w1,nv**2*no,nul,w2,nv)

     !ref1 = ref1 + w2(1:nv*no)

     !call print_norm(w2,i8*nv*no,nnorm,.true.)
     !call print_norm(ref1,i8*nv*no,norm,.true.)
     !!write (*,*)' DEBUG A1/TOT:',sqrt(nnorm),sqrt(norm)


     !!DEBUG SINGLES B1
     !!****************
     !call array_reorder_4d( p10, gooov, no, no, no, nv, [1,4,3,2], nul, w1)
     !call dgemm('n','n',nv, no,nv*no**2,m10,u,nv,w1,nv*no**2,nul,w2,nv)

     !ref1 = ref1 + w2(1:nv*no)

     !call print_norm(w2,i8*nv*no,nnorm,.true.)
     !call print_norm(ref1,i8*nv*no,norm,.true.)
     !!write (*,*)' DEBUG B1/TOT:',sqrt(nnorm),sqrt(norm)

     !!DEBUG SINGLES C1
     !!****************
     !call array_reorder_2d( p10, ovf, no, nv, [2,1], nul, w2)
     !call dgemv('n',no*nv,no*nv,p10,u,no*nv,w2,1, nul,w1,1)

     !ref1 = ref1 + w1(1:nv*no)

     !!call print_norm(w1,i8*nv*no,nnorm,.true.)
     !!call print_norm(ref1,i8*nv*no,norm,.true.)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! ref is not written after this point!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call mem_dealloc( w1 )
     call mem_dealloc( w2 )
     call mem_dealloc( w3 )
     call mem_dealloc( w4 )
     call mem_dealloc( w5 )

     call ass_D2to1(o1,h1,[nv,no])
     h1 = vof
     h1 => null()

     call mem_TurnONThread_Memory()
     !$OMP PARALLEL DEFAULT(NONE) PRIVATE(d,t,idx,pnv,pno,a,i,b,j,ns,pnv1,pnv2,pno1,pno2,&
     !$OMP d1,d2,t21,t22,w1,w2,w3,w4,w5,o,idx1,idx2,p1,p2,p3,p4,h1,h2,&
     !$OMP skiptrafo, skiptrafo2,oidx1,nidx1,oidx2,nidx2,i_idx,r1,r2,cyc,& 
     !$OMP nc,nc2,rpd,PS,ref) SHARED(pno_cv,pno_s,pno_t2,gvovo,goovv,gvvvv,&
     !$OMP vvf,goooo,Lvoov,pno_o2,govov,&
     !$OMP oof, maxsize, nspaces, ovf, gvvov, s_idx,o1,&
     !$OMP s_nidx,gooov, no, nv, p_idx, p_nidx,spacemax) 
     call init_threadmemvar()

     call mem_alloc( w1, maxsize )
     call mem_alloc( w2, maxsize )
     call mem_alloc( w3, maxsize )
     call mem_alloc( w4, maxsize )
     call mem_alloc( w5, maxsize )
     call mem_alloc( oidx1, spacemax, 3)
     call mem_alloc( oidx2, spacemax, 3)

     !$OMP DO SCHEDULE(DYNAMIC)
     LoopContribs:do ns = 1, nspaces


        if(.not.pno_cv(ns)%allocd)then

           cycle LoopContribs

        endif

        !The original space quantities carry no numbering
        d   => pno_cv(ns)%d
        t   => pno_t2(ns)%elm1
        idx => pno_cv(ns)%iaos
        pnv =  pno_cv(ns)%ns2
        pno =  pno_cv(ns)%n
        rpd =  pno_cv(ns)%rpd
        PS  =  pno_cv(ns)%PS

        o   => pno_o2(ns)%elm1

        o = 0.0E0_realk

        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  A2 Term !!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!

        !A2.1
        !Get the integral contribution, sort it first like the integrals then transform it
        call extract_from_gvovo_transform_add(gvovo,w1,w2,d,o,idx,rpd,pno,no,pnv,nv,&
           &pno_cv(ns)%n,PS)

        !A2.2
        call add_A22_contribution_simple(gvvvv,w1,w2,w3,d,t,o,pno,no,pnv,nv,pno_cv(ns)%n,PS)


        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  E2 Term part1, B2 !! 
        !!!!!!!!!!!!!!!!!!!!!!!!!

        call get_free_summation_for_current_aibj(no,ns,pno_cv,pno_S,pno_t2,o,&
           &w1,w2,w3,w4,w5,goooo,govov,vvf,nspaces)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  E2 Term part2, C2, D2 !! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Loop only over the indices which have a common index with the current
        !pair index, this could in principle also be solved with if statements
        !in the previous full loop and only doing the following work in a subset
        call get_common_idx_summation_for_current_aibj(no,ns,pno_cv,pno_S,pno_t2,&
           &o,w1,w2,w3,w4,w5,goovv,govov,Lvoov,oof,p_idx,p_nidx,oidx1,oidx2,nspaces,ref,nv)



        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  B1 Term !!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!


        !get gooov(kilc) as klic and transform c to pno basis
        call ass_D1to4( w3,    p3, [pno,pno,no, nv] )
        call ass_D1to4( gooov, p2, [no, no, no, nv] )
        do b=1,nv
           do i=1,no
              do j=1,pno
                 do a=1,pno
                    p3(a,j,i,b) = p2(idx(a), i, idx(j),b)
                 enddo
              enddo
           enddo
        enddo

        ! transform c such that d(c\bar{c})^T w(kli,c)^T = w1(\bar{c}kli)
        call dgemm('t','t',pnv,pno**2*no,nv,p10,d,nv,w3,pno**2*no,nul,w1,pnv)

        !get u from current amplitudes contract and transform back to local basis, as u(\bar{a}\bar{c}kl) from akcl
        call array_reorder_4d( p20, t, pnv,pno,pnv,pno, [1,3,2,4], nul,w2)
        call array_reorder_4d( m10, t, pnv,pno,pnv,pno, [1,3,4,2], p10,w2)


        ! carry out w2(\bar{a}\bar{c} kl) w1(\bar{c} kl i) = omega1{\bar{a}i}
        call dgemm('n','n',pnv,no,pnv*pno**2,m10,w2,pnv,w1,pnv*pno**2, nul, w3,pnv)
        !transform d(a\bar{a}) omega1{\bar{a} i} -> o1(a,i)
        !$OMP CRITICAL
        call dgemm('n','n',nv, no,pnv, p10,d, nv, w3, pnv, p10,o1,nv)
        !$OMP END CRITICAL


        d   => null()
        t   => null()
        idx => null()
        o   => null()
        pnv =  0
        pno =  0 
     enddo LoopContribs
     !$OMP END DO NOWAIT


     ! Add the missing singles contributions
     !$OMP DO SCHEDULE(DYNAMIC)
     LoopSingles: do nc=1,no


        !loop over all spaces in which the corresponding occupied index occurs
        OverlapLoop: do nc2=1,s_nidx(nc)


           ns    = s_idx(1,nc2,nc)
           i_idx = s_idx(2,nc2,nc)


           if(.not.pno_cv(ns)%allocd)then

              cycle OverlapLoop

           endif


           d   => pno_cv(ns)%d
           t   => pno_t2(ns)%elm1
           idx => pno_cv(ns)%iaos
           pnv =  pno_cv(ns)%ns2
           pno =  pno_cv(ns)%n

           !!!!!!!!!!!!!!!!!!!!!!!!!
           !!!  A1 Term !!!!!!!!!!!!
           !!!!!!!!!!!!!!!!!!!!!!!!!

           !extract gvvov(adkc) as dakc and transform d and c, keep a since singles
           !are constructed fully -> akcd
           call ass_D1to4( w1,    p1, [nv, nv,pno,nv] )
           call ass_D1to4( gvvov, p2, [nv, nv, no, nv] )
           do b=1,nv
              do j=1,pno
                 do i=1,nv
                    do a=1,nv
                       p1(i,a,j,b) = p2(a,i, idx(j),b)
                    enddo
                 enddo
              enddo
           enddo
           p1 => null()
           p2 => null()

           call dgemm('n','n',nv*nv*pno,pnv,nv,p10,w1,nv*nv*pno,d,nv,nul,w2,nv*nv*pno)
           call dgemm('t','n',nv*pno*pnv,pnv,nv,p10,w2,nv,d,nv,nul,w1,nv*pno*pnv)

           !extract amplitudes as u kcdi with i = nc2
           call ass_D1to4( w3, p3, [pno,pnv,pnv,1] )
           call ass_D1to4( t,  p2, [pnv,pno,pnv,pno] )
           do b=1,pnv
              do j=1,pno
                 do a=1,pnv
                    p3(j,a,b,1) = p20 * p2(a,j,b,i_idx) - p2(a,i_idx,b,j)
                 enddo
              enddo
           enddo
           p2 => null()
           p3 => null()

           call dgemv('n',nv,pno*pnv*pnv,p10,w1,nv,w3, 1, nul,w2,1)

           !$OMP CRITICAL
           o1(:,nc) = o1(:,nc) + w2(1:nv)
           !$OMP END CRITICAL

           !!!!!!!!!!!!!!!!!!!!!!!!!
           !!!  C1 Term !!!!!!!!!!!!
           !!!!!!!!!!!!!!!!!!!!!!!!!

           !extract amplitudes as u aick with i = nc2
           call ass_D1to4( w3, p3, [pnv,1,pnv,pno] )
           call ass_D1to4( t,  p2, [pnv,pno,pnv,pno] )
           do b=1,pnv
              do j=1,pno
                 do a=1,pnv
                    p3(a,1,b,j) = p20 * p2(a,i_idx,b,j) - p2(a,j,b,i_idx)
                 enddo
              enddo
           enddo
           p2 => null()
           p3 => null()

           !get fock matrix in the current space and transform virtual idx
           call ass_D1to2( w1, r1, [pno,nv] )
           call ass_D1to2( ovf, r2,[ no,nv] )
           do b=1,nv
              do j=1,pno
                 r1(j,b) = r2(idx(j),b)
              enddo
           enddo
           r1 => null()
           r2 => null()
           call dgemm('t','t',pnv,pno,nv,p10,d,nv,w1,pno,nul,w2,pnv)

           !contract amplitudes with fock matrix and transform back
           call dgemv('n',pnv,pno*pnv,p10,w3,pnv,w2,1,nul,w1,1)

           !transform back
           !$OMP CRITICAL
           call dgemv('n',nv,pnv,p10,d,nv,w1,1,p10,o1(1,nc),1)
           !$OMP END CRITICAL


           d     => null()
           t     => null()
           idx   => null()
           pnv   =  0
           pno   =  0 
           ns    =  0
           i_idx =  0
        enddo OverlapLoop
     enddo LoopSingles
     !$OMP END DO NOWAIT

     call mem_dealloc( w1 )
     call mem_dealloc( w2 )
     call mem_dealloc( w3 )
     call mem_dealloc( w4 )
     call mem_dealloc( w5 )
     call mem_dealloc( oidx1 )
     call mem_dealloc( oidx2 )
     o => null()

     call collect_thread_memory()
     !$OMP END PARALLEL
     call mem_TurnOffThread_Memory()


     !this subroutine assumes that symmetrization has already occured and only a
     !backtransformation to the original space is carried out
     !call print_norm(gvovo,o2v2)
     !do ns =1 ,nspaces
     !   call print_norm(pno_o2(ns))
     !enddo
     call backtransform_omegas(pno_o2,pno_cv,o2,nspaces,no,nv)

     !call print_norm(o2,o2v2)
     !call print_norm(o1,i8*no*nv)

     do ns = 1, nspaces

        if(  pno_cv(ns)%allocd )then
           call array_free( pno_t2(ns) )
           call array_free( pno_o2(ns) )
        endif

     enddo

     call mem_dealloc( pno_t2 )
     call mem_dealloc( pno_o2 )
     call mem_dealloc( gvvvv )
     call mem_dealloc( gvovo )
     call mem_dealloc( govov )
     call mem_dealloc( goooo )
     call mem_dealloc( goovv )
     call mem_dealloc( Lvoov )
     call mem_dealloc( gvvov )
     call mem_dealloc( gooov )
     call mem_dealloc( p_idx )
     call mem_dealloc( p_nidx )
     call mem_dealloc( s_nidx )
     call mem_dealloc( s_idx )
     call mem_dealloc( vof )
     call mem_dealloc( ovf )

     !o2 = 0.0E0_realk
     !o1 = 0.0E0_realk

  end subroutine get_ccsd_residual_pno_style

  subroutine get_overlap_idx(n1,n2,cv,idx,nidx,ndidx1,ndidx2)
     implicit none
     integer,intent(in) :: n1,n2
     type(PNOSpaceInfo), intent(in) :: cv(:)
     integer, intent(out) :: idx(:,:),nidx 
     integer, intent(out), optional :: ndidx1,ndidx2
     integer :: n,nc1,nc2,pos

     if(present(ndidx1)) ndidx1 = 0
     if(present(ndidx2)) ndidx2 = 0

     if(n1/=n2)then
        nidx = 0
        pos  = 0
        do nc1=1,cv(n1)%n
           do nc2=1,cv(n2)%n
              if(cv(n1)%iaos(nc1) == cv(n2)%iaos(nc2))then
                 pos = pos + 1
                 ! pos in first that equals second 
                 idx(pos,1) = nc1
                 ! pos in second that equals first
                 idx(pos,2) = nc2
                 ! the aos idx they refer to
                 idx(pos,3) = cv(n2)%iaos(nc2)
              endif
           enddo
        enddo
        nidx = pos


        if(present(ndidx1))then
           ndidx1 = 0
           do nc1 = 1,cv(n1)%n
              do n=1, nidx
                 if(cv(n1)%iaos(nc1) /= idx(n,3))then
                    pos    = pos    + 1
                    ndidx1 = ndidx1 + 1
                    ! pos in first that equals second 
                    idx(pos,1) = nc1
                    ! pos in second that equals first - set invalid since there is none
                    idx(pos,2) = -1
                    ! the aos idx they refer to
                    idx(pos,3) = cv(n1)%iaos(nc1)
                 endif
              enddo
           enddo
        endif

        if(present(ndidx2))then
           ndidx2 = 0
           do nc2 = 1,cv(n2)%n
              do n=1, nidx
                 if(cv(n2)%iaos(nc2) /= idx(n,3))then
                    pos    = pos    + 1
                    ndidx2 = ndidx2 + 1
                    ! pos in first that equals second - set invalid since there is none
                    idx(pos,1) = -1
                    ! pos in second that equals first
                    idx(pos,2) = nc2
                    ! the aos idx they refer to
                    idx(pos,3) = cv(n2)%iaos(nc2)
                 endif
              enddo
           enddo
        endif

     else
        !just copy if they are the same
        nidx = cv(n1)%n
        do nc1=1,cv(n1)%n
           ! pos in first that equals second 
           idx(nc1,1) = nc1
           ! pos in second that equals first
           idx(nc1,2) = nc1
           ! the aos idx they refer to
           idx(nc1,3) = cv(n1)%iaos(nc1)
        enddo
        pos = nidx
     endif

     if(nidx==0)then
        print *, nidx,n1,n2
        print *, cv(n1)%iaos
        print *, cv(n2)%iaos
        print *,"ONLY CALL THIS IF THERE ARE COMMON INDICES:something wrong check subroutine and input"
        stop 0 
     endif


  end subroutine get_overlap_idx

  !>\brief extract the information from the SpaceInfo structure in a nice and useful shape outside of this routine
  !>\author Patrick Ettenhuber
  !>\date december 2013
  subroutine get_overlap_ptr(n1,n2,pS,tr1,tr2,st,S,ldS,U,ldU,VT,ldVT,red1,red2,ns1,ns2)
     implicit none
     integer,intent(in) :: n1,n2
     type(PNOSpaceInfo), intent(in) :: pS(:)
     character, intent(out) :: tr1,tr2
     logical, intent(out) :: st
     real(realk), pointer, intent(out) :: S(:,:)
     integer, intent(out) :: ldS
     real(realk), pointer, intent(out) :: U(:,:),VT(:,:)
     integer, intent(out) :: ldU, ldVT
     integer, intent(out),optional :: red1, red2, ns1,ns2
     integer :: Sidx

     !Get the overlap identificaton and transformation props

     if(n1>n2)then

        !trafo from n1 to n2 
        Sidx      =  (n2 - n1 + 1) + n1 * (n1 - 1 )/2
        tr1       =  't'
        tr2       =  'n'
        st        =  .false.
        S         => pS(Sidx)%d
        ldS       =  pS(Sidx)%red1
        U         => pS(Sidx)%s1
        VT        => pS(Sidx)%s2
        ldU       =  pS(Sidx)%ns1
        ldVT      =  pS(Sidx)%red2

        if(present(red1))red1 = pS(Sidx)%red1
        if(present(red2))red2 = pS(Sidx)%red2
        if(present(ns1)) ns1 = pS(Sidx)%ns1
        if(present(ns2)) ns2 = pS(Sidx)%ns2
        !check if correct matrix was chosen
        if(pS(Sidx)%iaos(1)/=n1.or.pS(Sidx)%iaos(2)/=n2)then
           print *,"S mat wrong",pS(Sidx)%iaos(1),n1,pS(Sidx)%iaos(2),n2
        endif

     else if(n2>n1)then

        !trafo from n2 to n1
        Sidx      =  (n1 - n2 + 1) + n2 * (n2 - 1 )/2
        tr1       =  'n'
        tr2       =  't'
        st        =  .false.
        S         => pS(Sidx)%d
        ldS       =  pS(Sidx)%red1
        U         => pS(Sidx)%s1
        VT        => pS(Sidx)%s2
        ldU       =  pS(Sidx)%ns1
        ldVT      =  pS(Sidx)%red2

        if(present(red1))red1 = pS(Sidx)%red1
        if(present(red2))red2 = pS(Sidx)%red2
        if(present(ns1)) ns1 = pS(Sidx)%ns1
        if(present(ns2)) ns2 = pS(Sidx)%ns2
        !check if correct matrix was chosen
        if(pS(Sidx)%iaos(1)/=n2.or.pS(Sidx)%iaos(2)/=n1)then
           print *,"S mat wrong",pS(Sidx)%iaos(1),n2,pS(Sidx)%iaos(2),n1
        endif

     else

        !skip the transformation if the amplitudes reference the same space
        st        =  .true.
        S         => null()
        ldS       =  0

     endif
  end subroutine get_overlap_ptr


  subroutine get_pair_space_info(cv,p_idx,p_nidx,s_idx,s_nidx,ns,no)
     implicit none
     integer, intent(in) :: ns,no
     type(PNOSpaceInfo),intent(in) :: cv(ns)
     integer,intent(inout) :: p_idx(ns,ns),p_nidx(ns),s_idx(2,ns,no),s_nidx(no)
     integer :: cntr,n1,n2,k,l
     p_idx  = -1
     p_nidx = -1
     s_idx  = -1
     s_nidx = -1
     SpaceLoop:do n1 = 1, ns
        !set the first index to be the space itself
        cntr        = 1
        p_nidx(n1)  = 1
        p_idx(1,n1) = n1
        if(cv(n1)%allocd)then
           idxloop: do k = 1,cv(n1)%n
              !search for indices in the pair space
              SpaceLoop2:do n2 = 1, ns
                 if(n1/=n2.and.cv(n2)%allocd)then
                    do l = 1, cv(n2)%n
                       if(cv(n1)%iaos(k)==cv(n2)%iaos(l))then
                          cntr           = cntr+1
                          p_nidx(n1)     = cntr
                          p_idx(cntr,n1) = n2
                          cycle SpaceLoop2
                       endif
                    enddo
                 endif
              enddo SpaceLoop2
           enddo idxloop
        endif
     enddo SpaceLoop

     !print *,"DEBUG: print pair info"
     !do n1=1,ns
     !  print *,(p_idx(k,n1),k=1,p_nidx(n1))
     !enddo

     !search for the occupied indices in the spaces
     occupiedloop: do n1=1,no
        cntr = 0
        SpaceLoop3:do n2 = 1, ns
           if(cv(n2)%allocd)then
              do k = 1,cv(n2)%n
                 if(cv(n2)%iaos(k) == n1 )then
                    cntr = cntr + 1
                    s_nidx(n1) = cntr
                    s_idx(1,cntr,n1) = n2
                    s_idx(2,cntr,n1) = k
                 endif
              enddo
           endif
        enddo SpaceLoop3
     enddo occupiedloop
  end subroutine get_pair_space_info


  subroutine backtransform_omegas(pno_o2,pno_cv,o2,n,no,nv)
     implicit none
     integer,intent(in) :: n,no,nv
     type(array), intent(in) :: pno_o2(n)
     type(PNOSpaceInfo),intent(in) :: pno_cv(n)
     real(realk), intent(inout) :: o2(nv,no,nv,no)
     integer :: ns,pno,pnv,i,j,a,b, rpd
     real(realk), pointer :: tmp1(:),tmp2(:),d(:,:),po2(:,:,:,:), w1(:,:,:,:)
     real(realk) :: one,nul
     integer, pointer :: idx(:)
     logical :: PS

     one = 1.0E0_realk
     nul = 0.0E0_realk

     call mem_alloc(tmp1,nv**2*no**2)
     call mem_alloc(tmp2,nv**2*no**2)

     do ns = 1, n

        if(pno_cv(ns)%allocd)then

           pno =  pno_cv(ns)%n
           rpd =  pno_cv(ns)%rpd
           pnv =  pno_cv(ns)%ns2
           PS  =  pno_cv(ns)%PS
           d   => pno_cv(ns)%d
           po2 => pno_o2(ns)%elm4
           idx => pno_cv(ns)%iaos


           !For both cases the order should be correct after the
           !backtransformation. only in the first case the order is important,
           !since the transposition in the second part leads to the symmetric
           !contribution
           if( PS )then

              call dgemm( 'n', 'n', nv, pnv, pnv, one, d, nv, po2, pnv, nul, tmp1, nv )
              call dgemm( 'n', 't', nv, nv, pnv, one, tmp1, nv, d, nv, nul, tmp2, nv )
              tmp1(1:nv**2) = tmp2(1:nv**2)

           else

              call dgemm( 'n', 'n', nv, pno**2*pnv, pnv, one, d, nv, po2, pnv, nul, tmp1, nv )
              call array_reorder_4d( one, tmp1, nv, pno, pnv, pno, [3,4,1,2], nul, tmp2)
              call dgemm( 'n', 'n', nv, pno**2*nv, pnv, one, d, nv, tmp2, pnv, nul, tmp1, nv )

           endif


           !sort the contribution back and add up, again, because we assume a
           !symmetrized contribution in pno_o2 we can add up without taking care
           call ass_D1to4(tmp1,w1,[nv,rpd,nv,rpd])

           if(pno/=2.or.pno_cv(ns)%is_FA_space)then
              do j = 1, pno
                 do b = 1, nv
                    do i = 1, pno
                       do a = 1, nv
                          o2(a,idx(i),b,idx(j)) = o2(a,idx(i),b,idx(j)) + w1(a,i,b,j)
                       enddo
                    enddo
                 enddo
              enddo
           else
              i = 1
              j = 2

              if( PS )then
                 do b = 1, nv
                    do a = 1, nv
                       o2(a,idx(i),b,idx(j)) = o2(a,idx(i),b,idx(j)) + w1(a,1,b,1)
                       o2(a,idx(j),b,idx(i)) = o2(a,idx(j),b,idx(i)) + w1(b,1,a,1)
                    enddo
                 enddo
              else
                 do b = 1, nv
                    do a = 1, nv
                       o2(a,idx(i),b,idx(j)) = o2(a,idx(i),b,idx(j)) + w1(a,i,b,j)
                       o2(a,idx(j),b,idx(i)) = o2(a,idx(j),b,idx(i)) + w1(b,i,a,j)
                    enddo
                 enddo
              endif
           endif

           pno =  0
           pnv =  0
           d   => null()
           po2 => null()
           idx => null()
           w1  => null()

        endif
     enddo

     call mem_dealloc(tmp1)
     call mem_dealloc(tmp2)

  end subroutine backtransform_omegas

  subroutine check_if_contributes(n1,n2,cv,S,cyc)
     implicit none
     integer,intent(in) :: n1,n2
     type(PNOSpaceInfo),intent(in) :: cv(:),S(:)
     logical, intent(out) :: cyc
     integer :: Sidx

     cyc = .false.    

     !if the trafo matrix has been screened away, cycle
     if(n1>n2)then

        !trafo from n1 to n2 
        Sidx      =  (n2 - n1 + 1) + n1 * (n1 - 1 )/2
        cyc       =  .not. S(Sidx)%allocd

     else if(n2>n1)then

        !trafo from n2 to n1
        Sidx      =  (n1 - n2 + 1) + n2 * (n2 - 1 )/2
        cyc       =  .not. S(Sidx)%allocd

     else

        !do not skip prematurely
        cyc       =  .false.

     endif

     !or if the contribution n2 does not exist ( only important if n1==n2, else
     !the overlap will not exist in first place )
     cyc = ( cyc .or. .not. cv(n2)%allocd )

  end subroutine check_if_contributes

  !> \brief This routine calculates the overlap transformation from one PNO
  !space to another, the spaces are specified by their identification numbers n1
  !and n2. The overlap obtained from the SVD is saved in s and the necessary
  !information is extracted from S by the call to get_overlap_ptr. Pos of
  !overlap specifies wheter the overlap is the first or the second argument in
  !an untransposed gemm.
  !right now the routine assumes that A is never transposed. optionally
  !the pointers ptr and ptr2 can be assoicated with the transformed data and the
  !matrix not containing relevant data. this simplifies the code in the loops a
  !bit
  !> \author Patrick Ettenhuber
  !> \date december 2013
  ! 
  ! TODO: make this routine independent of the internal mem_allocations
  subroutine do_overlap_trafo(ns1,ns2,pos_of_overlap,S,m,n,k,A,C,ptr,ptr2,pC)
     implicit none
     integer, intent(in) :: pos_of_overlap,m,n,k,ns1,ns2
     type(PNOSpaceInfo),intent(inout) :: S(:)
     real(realk),intent(in),target :: A(:)
     real(realk),intent(inout),target :: C(:)
     real(realk),pointer,optional :: ptr(:),ptr2(:)
     real(realk),optional :: pC
     real(realk),pointer :: S1(:,:)
     real(realk),pointer :: tmp1(:),tmp2(:), U(:,:), VT(:,:)
     integer :: ldS1,ldU,ldVT
     logical :: skiptrafo
     character :: tr1,tr2
     real(realk), parameter :: nul = 0.0E0_realk
     real(realk), parameter :: p10 = 1.0E0_realk
     real(realk) :: preC


     preC = 0.0E0_realk
     if(present(pC))preC = pC

     call get_overlap_ptr(ns1,ns2,S,tr1,tr2,skiptrafo,S1,ldS1,U=U,ldU=ldU,VT=VT,ldVT=ldVT)


     if(.not. skiptrafo)then

        select case(pos_of_overlap)
        case(1)

           call mem_alloc(tmp1,ldS1 * n)
           call mem_alloc(tmp2,ldVT * n)
           if(tr2=='t')then
              call dgemm('t','n', ldS1, n, ldU, p10, U, ldU, A, k, nul, tmp1, ldS1 )
              call dgemm('t','n', ldVT, n, ldS1,p10, S1,ldS1,tmp1,ldS1,nul,tmp2,ldVT)
              call dgemm('t','n',m,n,ldVT,p10, VT,ldVT,tmp2,ldVT,preC,C,m)
           else if(tr2=='n')then
              call dgemm('n','n',ldVT, n,k, p10, VT,ldVT,A,k,nul,tmp2,ldVT)
              call dgemm('n','n',ldS1,n,ldVT,p10, S1,ldS1,tmp2,ldVT,nul,tmp1,ldS1)
              call dgemm('n','n',m,n,ldS1, p10, U,ldU,tmp1,ldS1,preC,C,m)
           else
              call lsquit("ERROR(do_overlap_trafo):this should never happen, check get_overlap_ptr",-1)
           endif
           call mem_dealloc(tmp1)
           call mem_dealloc(tmp2)

        case(2)
           call mem_alloc(tmp1,ldS1 * m)
           call mem_alloc(tmp2,ldVT * m)
           if(tr1=='t')then
              call dgemm('n','t',m,ldVT,k, p10, A,m,VT,ldVT,nul,tmp2,m)
              call dgemm('n','t',m,ldS1,ldVT,p10, tmp2,m,S1,ldS1,nul,tmp1,m)
              call dgemm('n','t',m, n,ldS1,p10, tmp1,m,U,ldU,preC,C,m)
           else if(tr1=='n')then
              call dgemm('n','n',m, ldS1,k, p10, A,m,U,ldU,nul,tmp1,m)
              call dgemm('n','n',m,ldVT,ldS1,p10, tmp1,m,S1,ldS1,nul,tmp2,m)
              call dgemm('n','n',m,n,ldVT, p10, tmp2,m,VT,ldVT,preC,C,m)
           else
              call lsquit("ERROR(do_overlap_trafo):this should never happen, check get_overlap_ptr",-1)
           endif
           call mem_dealloc(tmp1)
           call mem_dealloc(tmp2)

        case default

           call lsquit("ERROR(do_overlap_trafo): wrong selection of pos_of_overlap",-1)

        end select

        !associate the pointer to the result
        if(present(ptr)) ptr  => C
        !associate the pointer to the input matrix
        if(present(ptr2))ptr2 => A

     else

        !Do the association the other way round, since the data are in the input
        !matrix
        if(present(ptr)) ptr  => A
        if(present(ptr2))ptr2 => C

     endif

  end subroutine do_overlap_trafo


  !> \brief Get the overlap matrices specifying the transformation from one PNO
  !space to another. Here the overlap screening with the singular value
  !decomposition is used to screen away some contributions.
  !> \author Patrick Ettenhuber
  !> \date december 2013
  subroutine get_pno_overlap_matrices(no,nv,pno_cv,pno_S,n,with_svd)
     implicit none
     integer :: no, nv, n
     type(PNOSpaceInfo),intent(in) :: pno_cv(n)
     type(PNOSpaceInfo),intent(inout) :: pno_S(n*(n-1)/2)
     logical, intent(in) :: with_svd
     integer :: i, j, c, t1,t2,dg, n1
     integer :: ns1,ns2,INFO,lwork,mindim,maxdim,red1,red2,kerdim,diag,remove
     real(realk),pointer:: s1(:,:), s2(:,:), sv(:),U(:,:), VT(:,:),work(:)
     real(realk) :: norm,thr
     real(realk),parameter :: p10 = 1.0E0_realk
     real(realk),parameter :: nul = 0.0E0_realk
     logical :: keep_pair
     integer :: allremoved, ofmindim, ofmaxdim, allocpcount

     if( DECinfo%noPNOoverlaptrunc ) then
        thr = -1.0E0_realk * huge(thr)
     else
        thr = DECinfo%PNOoverlapthr
     endif

     allocpcount= 0
     allremoved = 0
     ofmindim   = 0
     ofmaxdim   = 0

     call mem_TurnONThread_Memory()
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP REDUCTION(+:allremoved,ofmindim,ofmaxdim,allocpcount)&
     !$OMP SHARED(pno_cv,pno_S,n,no,nv,with_svd,thr)&
     !$OMP PRIVATE(ns1,ns2,i,j,c,s1,s2,norm,sv,U,VT,work,remove,&
     !$OMP lwork,info,diag,kerdim,red1,red2,maxdim,mindim,dg,&
     !$OMP keep_pair)
     call init_threadmemvar()

     !CURRENT HACK FOR PGI COMPILER, SOMETHING WITH THE ALLOCATIONS IN THE LOOP
     !(AND MAYBE STACK MEMORY), FIXME: MOVE ALLOCATION OUTSIDE OF PARALLEL REGION

     !$OMP SINGLE

     !OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
     do i=1,n
        do j=1,n

           if(j>=i) cycle

           ! COUNT UPPER TRIANGULAR ELEMENTS WITHOUT DIAGONAL ELEMENTS
           c = (j - i + 1) + i*(i-1)/2

           ns1 = pno_cv(i)%ns2
           ns2 = pno_cv(j)%ns2

           pno_S(c)%ns1 = ns1
           pno_S(c)%ns2 = ns2

           keep_pair = ( pno_cv(i)%allocd .and. pno_cv(j)%allocd )

           if( keep_pair )then


              pno_S(c)%s_associated = .false.
              call mem_alloc(pno_S(c)%d,ns1,ns2)

              s1 => pno_cv(i)%d
              s2 => pno_cv(j)%d

              call dgemm('t','n',ns1,ns2,nv,p10,s1,nv,s2,nv,nul,pno_S(c)%d,ns1)

              if(with_svd)then
                 allocpcount = allocpcount + 1
                 !get the minimum dimension, the maximum dimension and the dimension
                 !of the kernel of the transformation
                 mindim = min(ns1,ns2)
                 maxdim = max(ns1,ns2)
                 kerdim = maxdim - mindim

                 !Characterize the type of overlap just produced, does it need to be
                 !considered at all -> calculate the singular values for checking
                 lwork = max(1,3*mindim+maxdim,5*mindim)
                 call mem_alloc(sv,min(ns1,ns2))
                 call mem_alloc(work,lwork)
                 call mem_alloc(U,ns1,ns1)
                 call mem_alloc(VT,ns2,ns2)
                 sv = 0.0E0_realk; INFO=0
                 call dgesvd('A','A',ns1,ns2,pno_S(c)%d,ns1,sv,U,ns1,VT,ns2,work,lwork,INFO)
                 if(INFO/=0)call &
                    &lsquit("ERROR(get_pno_overlap_matrices): dgesvd failed",-1)

                 !screen singular values according to a predefined threshold
                 do diag=1,mindim
                    if(sv(diag)<thr)then
                       exit
                    endif
                 enddo
                 ! go one step back to where the singular value is still above thr
                 ! and calculate the number of elements to remove additionally to the
                 ! kernel dimension
                 if(diag/=mindim)diag = diag - 1
                 remove = mindim - diag

                 ! Find the number of elements in the trafo, note that kerdim == 0 if
                 ! ns1 == ns2, therefore an if else is enough
                 if(ns1>ns2)then
                    red1 = ns1 - remove - kerdim
                    red2 = ns2 - remove
                 else
                    red1 = ns1 - remove
                    red2 = ns2 - remove - kerdim
                 endif

                 ofmindim   = ofmindim   + mindim
                 ofmaxdim   = ofmaxdim   + maxdim
                 allremoved = allremoved + remove

                 keep_pair = ( red1 > 0 .and. red2 > 0 )

                 if(red1/=diag.or.red2/=diag)call &
                    &lsquit("ERROR(get_pno_overlap_matrices)calculated wrong dimensions",-1)


                 call mem_dealloc( pno_S(c)%d )

                 if( keep_pair )then
                    pno_S(c)%s_associated = .true.
                    call mem_alloc( pno_S(c)%s1, ns1,  red1 )
                    call mem_alloc( pno_S(c)%s2, red2, ns2  )
                    call mem_alloc( pno_S(c)%d,  red1, red2 )
                    pno_S(c)%red1 = red1
                    pno_S(c)%red2 = red2

                    pno_S(c)%s1 = U(:,1:red1)
                    pno_S(c)%s2 = VT(1:red2,:)
                    pno_S(c)%d  = nul
                    do dg = 1, diag
                       pno_S(c)%d(dg,dg) = sv(dg)
                    enddo

                 endif


                 call mem_dealloc( U )
                 call mem_dealloc( VT )
                 call mem_dealloc( work )
                 call mem_dealloc( sv )
              endif


              if( keep_pair ) then

                 pno_S(c)%n = 2
                 call mem_alloc(pno_S(c)%iaos,pno_S(c)%n)
                 pno_S(c)%iaos = [i,j]

              endif

              pno_S(c)%allocd = keep_pair

           else

              pno_S(c)%allocd = .false.

           endif


        enddo
     enddo
     !$OMP END SINGLE
     !OMP END DO NOWAIT
     call collect_thread_memory()
     !$OMP END PARALLEL
     call mem_TurnOffThread_Memory()

     print *,"overlapscreening removed ",allremoved," of ",ofmindim,ofmaxdim,"in",allocpcount,"of",n*(n-1)/2,"pairs"

  end subroutine get_pno_overlap_matrices

  subroutine init_pno_residual(cv,o2,n)
     implicit none
     integer, intent(in) :: n
     type(PNOSpaceInfo),intent(in) :: cv(n)
     type(array),intent(inout) :: o2(n)
     integer :: nn, pnv,rpd

     do nn=1,n

        pnv = cv(nn)%ns2
        rpd = cv(nn)%rpd

        if(cv(nn)%allocd)then

           o2(nn) = array_init([pnv,rpd,pnv,rpd],4)

        endif
     enddo
  end subroutine init_pno_residual

  subroutine get_pno_amplitudes(t2,cv,pno_t2,n,no,nv)
     implicit none
     integer, intent(in) :: n,no,nv
     real(realk),intent(in) :: t2(:,:,:,:)
     type(PNOSpaceInfo), intent(in) :: cv(n)
     type(array), intent(inout) :: pno_t2(n)
     real(realk), pointer :: tmp1(:),tmp2(:)
     real(realk), pointer :: w1(:,:,:,:)
     integer :: nn, pnv, pno, a, b, i, j, rpd
     logical :: PS, FAspace
     call mem_alloc(tmp1,no**2*nv**2)
     call mem_alloc(tmp2,no**2*nv**2)

     do nn=1,n

        pnv     = cv(nn)%ns2
        pno     = cv(nn)%n
        rpd     = cv(nn)%rpd
        PS      = cv(nn)%PS
        FAspace = cv(nn)%is_FA_space

        if(cv(nn)%allocd)then

           pno_t2(nn) = array_init([pnv,rpd,pnv,rpd],4)

           call ass_D1to4(tmp1,w1,[nv,rpd,nv,rpd])

           if( PS )then
              i           = cv(nn)%iaos(1)
              j           = cv(nn)%iaos(2)
              w1(:,1,:,1) = t2(:,i,:,j)
           else
              do j=1,pno
                 do i=1,pno
                    w1(:,i,:,j) = t2(:,cv(nn)%iaos(i),:,cv(nn)%iaos(j))
                 enddo
              enddo
           endif

           if ( PS )then
              call dgemm('t','n',pnv,nv,nv,1.0E0_realk,cv(nn)%d,nv,w1,nv,0.0E0_realk,tmp2,pnv)
              call dgemm('n','n',pnv,pnv,nv,1.0E0_realk,tmp2,pnv,cv(nn)%d,nv,0.0E0_realk,pno_t2(nn)%elm1,pnv)
           else
              call dgemm('t','n',pnv,pno**2*nv,nv,1.0E0_realk,cv(nn)%d,nv,w1,nv,0.0E0_realk,tmp2,pnv)
              !the amplitudes are symmetric, also after the transformation,
              !therefore it is not important in which order they are stored
              call array_reorder_4d(1.0E0_realk,tmp2,pnv,pno,nv,pno,[3,4,1,2],0.0E0_realk,tmp1)
              call dgemm('t','n',pnv,pno**2*pnv,nv,1.0E0_realk,cv(nn)%d,nv,tmp1,nv,0.0E0_realk,pno_t2(nn)%elm1,pnv)
           endif

           w1 => null()

           !To avoid double counting -> this can be removed and savings introduced in
           !the algorithm if restrictions are put on the pair indices, but a bit more
           !complicated than this
           if( rpd==2 .and. .not.  cv(nn)%is_FA_space)then
              do i = 1, pno
                 pno_t2(nn)%elm4(:,i,:,i) = 0.0E0_realk
              enddo
           endif
        endif

     enddo

     call mem_dealloc(tmp1)
     call mem_dealloc(tmp2)
  end subroutine get_pno_amplitudes

  subroutine get_pno_trafo_matrices(no,nv,nb,t_mp2,cv,n,fj,f)
     implicit none
     !ARGUMENTS
     integer, intent(in) :: no, nv, nb, n
     real(realk), intent(in) :: t_mp2(nv,no,nv,no)
     type(PNOSpaceInfo),pointer :: cv(:)
     logical,intent(in) :: fj
     type(decfrag),intent(in),optional :: f
     !INTERNAL
     real(realk) :: virteival(nv),U(nv,nv),PD(nv,nv)
     integer :: i,j,oi,oj,counter, calc_parameters,det_parameters
     integer :: pno_gvvvv_size,find_pos(no,no)
     logical :: doit

     find_pos = -1

     counter = 0
     if(fj) counter = 1
     doi1 :do i = 1, no
        doj1: do j = i, no

           doit=.true.
           !check if both indices occur in occ EOS, if yes -> skip
           if(fj)then
              oiloop2: do oi = 1, f%noccEOS
                 if(f%idxo(oi) == i)then
                    do oj = 1, f%noccEOS
                       if(f%idxo(oj) == j)then
                          doit = .false.
                          exit oiloop2
                       endif
                    enddo
                 endif
              enddo oiloop2
           endif

           if(doit)then
              counter = counter + 1
              find_pos(i,j) = counter
           endif
        enddo doj1
     enddo doi1

     if(counter /= n )then
        call lsquit("ERROR(get_pno_trafo_matrices):wrong counter",-1)
     endif

     calc_parameters = 0
     det_parameters  = 0
     pno_gvvvv_size  = 0

     call mem_TurnONThread_Memory()
     !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:calc_parameters,det_parameters,pno_gvvvv_size)&
     !$OMP SHARED(no,nv,nb,n,fj,f,DECinfo,cv,find_pos,t_mp2)&
     !$OMP PRIVATE(counter,virteival,U,PD,doit)
     call init_threadmemvar()

     if(fj)then

        if(.not.associated(f%VirtMat))then
           call lsquit("Error(get_pno_trafo_matrices)Fragment Correlation density matrix not allocated",-1)
        endif


        !$OMP SINGLE
        call solve_eigenvalue_problem_unitoverlap(nv,f%VirtMat,virteival,U)
        call truncate_trafo_mat_from_EV(U,virteival,nv,cv(1),ext_thr=DECinfo%EOSPNOthr)
        call mem_alloc(cv(1)%iaos,f%noccEOS)
        cv(1)%n            = f%noccEOS
        cv(1)%rpd          = f%noccEOS
        cv(1)%iaos         = f%idxo
        cv(1)%is_FA_space  = .true.
        cv(1)%PS           = .false.
        counter            = 1

        calc_parameters = calc_parameters + cv(1)%ns2**2*cv(1)%n**2
        det_parameters  = det_parameters  + cv(1)%ns2**2*cv(1)%n**2
        pno_gvvvv_size  = pno_gvvvv_size  + cv(1)%ns2**4

        if(.not.cv(1)%allocd)then
           call lsquit("ERROR(get_pno_trafo_matrices):EOS does not contribute&
              & according to the current threshold, skipping this fragment should be&
              & implemented",-1)
        endif
        !$OMP END SINGLE

        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
        doi :do i = 1, no
           doj: do j = 1, no

              if(j<i.or.(i==1.and.j==1)) cycle doj

              !check if both indices occur in occ EOS, if yes -> skip
              doit = ( find_pos(i,j)/= -1 )

              !calculate the pair density matrix, diagonalize it, truncate the
              !respective transformation matrix and save it in c
              if(doit)then

                 counter = find_pos(i,j)
                 call calculate_pair_density_matrix(PD,t_mp2(:,i,:,j),nv,(i==j))
                 call solve_eigenvalue_problem_unitoverlap(nv,PD,virteival,U)
                 call truncate_trafo_mat_from_EV(U,virteival,nv,cv(counter))
                 cv(counter)%is_FA_space  = .false.

                 if(cv(counter)%allocd)then
                    if(i==j)then

                       cv(counter)%n    = 1
                       cv(counter)%rpd  = 1
                       call mem_alloc(cv(counter)%iaos,cv(counter)%n)
                       cv(counter)%iaos = [i]
                       cv(counter)%PS   = .false.

                       det_parameters   = det_parameters + cv(counter)%ns2*cv(counter)%ns2*cv(counter)%n**2

                    else

                       cv(counter)%n    = 2
                       if(DECinfo%PNOtriangular)then
                          cv(counter)%rpd  = 1
                          cv(counter)%PS   = .true.
                       else
                          cv(counter)%rpd  = 2
                          cv(counter)%PS   = .false.
                       endif

                       call mem_alloc(cv(counter)%iaos,cv(counter)%n)
                       cv(counter)%iaos = [i,j]

                       det_parameters   = det_parameters + cv(counter)%ns2*cv(counter)%ns2*2

                    endif

                    calc_parameters = calc_parameters + cv(counter)%ns2*cv(counter)%ns2*cv(counter)%rpd**2
                    pno_gvvvv_size  = pno_gvvvv_size  + cv(counter)%ns2**4

                 endif
              endif
           enddo doj
        enddo doi
        !$OMP END DO NOWAIT
     else
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
        doiful :do i = 1, no
           dojful: do j = 1, no

              if(j<i) cycle dojful

              counter = find_pos(i,j)

              call calculate_pair_density_matrix(PD,t_mp2(:,i,:,j),nv,(i==j))
              call solve_eigenvalue_problem_unitoverlap(nv,PD,virteival,U)
              call truncate_trafo_mat_from_EV(U,virteival,nv,cv(counter))
              cv(counter)%is_FA_space  = .false.

              if(cv(counter)%allocd)then

                 if(i==j)then

                    cv(counter)%n    = 1
                    cv(counter)%rpd  = 1
                    call mem_alloc(cv(counter)%iaos,cv(counter)%n)
                    cv(counter)%iaos = [i]
                    cv(counter)%PS   = .false.

                    det_parameters = det_parameters + cv(counter)%ns2*cv(counter)%ns2

                 else

                    cv(counter)%n = 2
                    if(DECinfo%PNOtriangular)then
                       cv(counter)%rpd  = 1
                       cv(counter)%PS   = .true.
                    else
                       cv(counter)%rpd  = 2
                       cv(counter)%PS   = .false.
                    endif

                    call mem_alloc(cv(counter)%iaos,cv(counter)%n)
                    cv(counter)%iaos = [i,j]

                    det_parameters = det_parameters + cv(counter)%ns2*cv(counter)%ns2*2

                 endif

                 calc_parameters = calc_parameters + cv(counter)%ns2*cv(counter)%ns2*cv(counter)%rpd**2
                 pno_gvvvv_size  = pno_gvvvv_size  + cv(counter)%ns2**4

              endif

           enddo dojful
        enddo doiful
        !$OMP END DO NOWAIT
     endif

     call collect_thread_memory()
     !$OMP END PARALLEL
     call mem_TurnOffThread_Memory()

     print *,"I have to determine",det_parameters," of ",no**2*nv**2," using ",calc_parameters
     !print *,"full gvvvv",nv**4," vs ",pno_gvvvv_size


  end subroutine get_pno_trafo_matrices

  !\brief Calculation of the pair density matrix from a set of MP2 amplitudes
  !for a given pair. The input amplitudes are a virt-virt block for indices (ij)
  !and the routine has to know wheter i==j
  !\author Patrick Ettenhuber
  subroutine calculate_pair_density_matrix(PD,tvv,nv,ieqj)
     implicit none
     !ARGUMENTS
     integer, intent(in) :: nv
     real(realk),intent(inout) :: PD(nv,nv)
     real(realk),intent(in)    :: tvv(nv,nv)
     logical, intent(in) :: ieqj
     !INTERNAL
     real(realk) :: tildetvv(nv,nv), fact

     !build tilde t from mp2 amplitudes, set the prefactor correct, i.e. if i==j
     !a factor of one half is needed in the construction of the pair density matrix
     fact = 2.0E0_realk
     if(ieqj) fact = 0.5E0_realk * fact

     tildetvv = 2.0E0_realk * fact * tvv
     call array_reorder_2d(-1.0E0_realk*fact,tvv,nv,nv,[2,1],1.0E0_realk,tildetvv)


     ! do the contractions of tilde t with t and save them in the output matrix
     call dgemm('n','t',nv,nv,nv,1.0E0_realk,tildetvv,nv,tvv,nv,0.0E0_realk,PD,nv)
     call dgemm('t','n',nv,nv,nv,1.0E0_realk,tildetvv,nv,tvv,nv,1.0E0_realk,PD,nv)

  end subroutine calculate_pair_density_matrix

  !\brief 
  subroutine truncate_trafo_mat_from_EV(U,EV,n,NU,ext_thr)
     implicit none
     !ARGUMENTS
     integer,intent(in) :: n
     real(realk), intent(in) :: U(n,n),EV(n)
     type(PNOSpaceInfo),intent(inout) :: NU
     real(realk), intent(in),optional :: ext_thr
     !INTERNAL
     integer :: i,nn
     real(realk) :: thr

     thr = DECinfo%simplePNOthr
     if(present(ext_thr)) thr = ext_thr

     !on finishing the loop i contains the position of the first element that should be in the
     !transformation
     do i = 1, n
        if(EV(i)>thr)then
           exit
        endif
     enddo
     ! n elements in the transformation
     nn = n - i + 1

     !DEBUGGING KEYWORD to avoid truncation
     if(DECinfo%noPNOtrunc.and..not.present(ext_thr))then 
        nn = n
        i  = 1
     endif
     if(DECinfo%noFAtrunc .and.     present(ext_thr))then 
        nn = n
        i  = 1
     endif

     if(DECinfo%PL>2.and.present(ext_thr))write(DECinfo%output,'("The FO trafo  matrix has dims",2I4)')n,nn


     NU%ns1 = n
     NU%ns2 = nn

     if(nn<=0)then

        NU%allocd = .false.

     else

        NU%s_associated = .false.
        call mem_alloc(NU%d,n,nn)

        if((DECinfo%noPNOtrafo.and..not.present(ext_thr)).or.&
           (DECinfo%noFAtrafo .and.     present(ext_thr))       )then

           NU%d = 0.0E0_realk
           do i = 1, n
              NU%d(i,i) = 1.0E0_realk
           enddo

        else

           NU%d = U(1:n,i:n)

        endif

        NU%allocd = .true.

     endif


  end subroutine truncate_trafo_mat_from_EV

  subroutine extract_from_gvovo_transform_add(gvovo,w1,w2,d,o,idx,rpd,pno,no,pnv,nv,n,PS)
     implicit none
     integer,intent(in) :: rpd,pno,no,pnv,nv,n
     integer,intent(in) :: idx(n)
     real(realk),intent(in) ::gvovo(nv,no,nv,no), d(nv,pnv)
     real(realk),intent(inout) :: o(pnv*pnv*pno*pno)
     real(realk),intent(inout) :: w1(nv,rpd,nv,rpd)
     real(realk),pointer,intent(inout) :: w2(:)
     logical, intent(in) :: PS
     real(realk),pointer :: check(:,:,:,:)
     integer :: i,j,a,b
     real(realk), parameter :: nul = 0.0E0_realk
     real(realk), parameter :: p10 = 1.0E0_realk

     if( PS )then
        do b=1,nv
           do a=1,nv
              w1(a,1,b,1) = gvovo(a,idx(1),b,idx(2))
           enddo
        enddo
     else
        do j=1,pno
           do b=1,nv
              do i=1,pno
                 do a=1,nv
                    w1(a,i,b,j) = gvovo(a,idx(i),b,idx(j))
                 enddo
              enddo
           enddo
        enddo
     endif

     if( PS )then
        !transform integral contribution, no symmetry here
        call dgemm( 't', 'n', pnv, nv,  nv, p10, d,  nv, w1, nv, nul, w2, pnv )
        call dgemm( 'n', 'n', pnv, pnv, nv, p10, w2, pnv, d, nv, nul, o, pnv )
     else
        !transform integral contribution, use symmetry  gvovo(aibj) => gvovo(\bar{b} j \bar{a} i)
        call dgemm( 't', 'n', pnv, pno**2*nv, nv, p10, d, nv, w1, nv, nul, w2, pnv )
        call array_reorder_4d( p10, w2, pnv, pno, nv, pno, [3,4,1,2], nul, w1 )
        call dgemm( 't', 'n', pnv, pno**2*pnv, nv, p10, d, nv, w1, nv, nul, o, pnv )
     endif

  end subroutine extract_from_gvovo_transform_add


  subroutine add_A22_contribution_simple(gvvvv,w1,w2,w3,d,t,o,pno,no,pnv,nv,n,PS)
     implicit none
     integer,intent(in) :: pno,no,pnv,nv,n
     real(realk), intent(in)    :: gvvvv(nv**4), t(pnv**2*pno**2), d(nv,pnv)
     real(realk), intent(out)   :: w1(:),w2(:),w3(:)
     real(realk), intent(inout) :: o(pnv**2*pno**2)
     logical, intent(in)        :: PS
     real(realk), parameter     :: nul = 0.0E0_realk
     real(realk), parameter     :: p10 = 1.0E0_realk
     !A2.2
     !transform to basis of space gvvvv(acbd) => gvvvv(\bar{a}\bar{c}\bar{b}\bar{d})
     call dgemm( 't', 'n', nv**3,     pnv, nv, p10, gvvvv, nv, d, nv, nul, w1, nv**3     )
     call dgemm( 't', 'n', nv**2*pnv, pnv, nv, p10, w1   , nv, d, nv, nul, w2, nv**2*pnv )
     call dgemm( 't', 'n', nv*pnv**2, pnv, nv, p10, w2   , nv, d, nv, nul, w1, nv*pnv**2 )
     call dgemm( 't', 'n', pnv**3,    pnv, nv, p10, w1   , nv, d, nv, nul, w2, pnv**3    )
     !end transformation, integrals are now in the order gvvvv(\bar{a}\bar{c}\bar{b}\bar{d})
     !reorder corresponding integrals (w1) to gvvvv(\bar{a}\bar{b}\bar{c}\bar{d})and 
     call array_reorder_4d( p10, w2, pnv, pnv, pnv, pnv, [1,3,2,4], nul, w1 )

     if( PS )then
        call dgemv('n',pnv**2,pnv**2,p10,w1,pnv**2,t,1,p10,o,1)
     else
        !amplitudes to t((\bar{c}\bar{d} i j)
        call array_reorder_4d( p10, t,  pnv, pno, pnv, pno, [1,3,2,4], nul, w2 )
        !contract the amplitudes and integrals to get the A2.2 contribution
        call dgemm( 'n', 'n', pnv**2, pno**2, pnv**2, p10, w1, pnv**2, w2, pnv**2, nul, w3, pnv**2 )
        call array_reorder_4d( p10, w3, pnv, pnv, pno, pno, [1,3,2,4], p10, o )
     endif
  end subroutine add_A22_contribution_simple


  subroutine get_free_summation_for_current_aibj(no,ns,pno_cv,pno_S,pno_t2,o2_space,&
        &w1,w2,w3,w4,w5,goooo,govov,vvf,nspaces)
     implicit none
     integer, intent(in) :: no,ns,nspaces
     type(PNOSpaceInfo), intent(inout) :: pno_cv(nspaces),pno_S(nspaces*(nspaces-1)/2)
     type(array), intent(in) :: pno_t2(nspaces)
     real(realk),pointer,intent(inout) :: o2_space(:)
     real(realk),pointer,intent(inout) :: w1(:),w2(:),w3(:),w4(:),w5(:)
     real(realk),intent(in) :: goooo(:),govov(:),vvf(:,:)
     character :: tr11,tr12,tr21,tr22
     real(realk),pointer :: p1(:,:,:,:), p2(:,:,:,:), p3(:,:,:,:), p4(:,:,:,:),h1(:), h2(:), r1(:,:),r2(:,:),d(:,:),d1(:,:)
     real(realk),pointer :: o(:),t(:),S1(:,:), t21(:)
     logical :: skiptrafo, cyc, PS, PS1
     integer :: space, a,i,b,j, nv,pno,pnv, ns2
     integer :: pno1,pnv1,Sidx1,ldS1, rpd1,rpd
     integer :: pair,paircontrib(2,2)
     integer :: order1(4)
     integer :: suborder(2)
     integer,pointer :: idx(:),idx1(:)
     integer(kind=8) :: o2v2
     integer, parameter :: paircontribs = 2
     real(realk), parameter :: nul =  0.0E0_realk
     real(realk), parameter :: p10 =  1.0E0_realk
     real(realk), parameter :: p20 =  2.0E0_realk
     real(realk), parameter :: m10 = -1.0E0_realk

     nv  =  pno_cv(ns)%ns1
     d   => pno_cv(ns)%d
     idx => pno_cv(ns)%iaos
     pnv =  pno_cv(ns)%ns2
     pno =  pno_cv(ns)%n
     rpd =  pno_cv(ns)%rpd
     PS  =  pno_cv(ns)%PS
     t   => pno_t2(ns)%elm1
     o   => o2_space

     o2v2 = (i8*no**2)*nv**2

     paircontrib(1:2,1) = [1,2]
     paircontrib(1:2,2) = [2,1]

     !!!!!!!!!!!!!!!!!!!!!!!!!
     !!!  E2 Term part 1!!!!!! -- continued in the following loop and after the loop
     !!!!!!!!!!!!!!!!!!!!!!!!!

     !the first step is transforming the vv Fock matrix to the pno space of (ij), instead of constructing the
     !full fock matrix and doing the trafo here, it might already be
     !constructed in the pno basis, probably at the expense of memory,
     !depending on the sizes
     call dgemm('t','n',pnv,nv,nv,p10, d, nv,vvf,nv,nul,w1,pnv)
     call dgemm('n','n',pnv,pnv,nv,p10, w1,pnv,d,nv,nul,w5,pnv)

     FullSpaceLoop: do ns2 = 1, nspaces

        call check_if_contributes(ns,ns2,pno_cv,pno_S,cyc)

        if(cyc)then

           cycle FullSpaceLoop

        endif

        d1   => pno_cv(ns2)%d
        t21  => pno_t2(ns2)%elm1
        idx1 => pno_cv(ns2)%iaos
        pnv1 =  pno_cv(ns2)%ns2
        pno1 =  pno_cv(ns2)%n
        rpd1 =  pno_cv(ns2)%rpd
        PS1  =  pno_cv(ns2)%PS



        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  E2 Term part1!!!!!!! - quadratic contribution
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !Get the integral contribution, sort it first like the integrals then transform it, govov
        call ass_D1to4( w1,    p1, [rpd1,nv,rpd1,nv] )
        call ass_D1to4( govov, p2, [no,   nv,no, nv] )

        if( PS1 )then

           do pair = 1, paircontribs


              p1(1,:,1,:) = p2(idx1(paircontrib(1,pair)),:,idx1(paircontrib(2,pair)),:)


              !transform integral contribution govov(ldkc) => govov(\bar{d} k l  \bar{c}) to the space of (ij) -> w2
              call dgemm( 'n', 'n', nv, pnv, nv, p10, w1, nv, d, nv, nul, w2, nv)
              call dgemm( 't', 'n', pnv1, pnv, nv, p10, d1, nv, w2, nv, nul, w1, pnv1 )

              ! Quadratic part of the E2 term use u^{bd}_{kl} (bkdl) as b,dkl
              call array_reorder_2d( p20, t21, pnv1, pnv1, paircontrib(:,3-pair), nul, w3)
              call array_reorder_2d( m10, t21, pnv1, pnv1, paircontrib(:,pair), p10, w3)
              call do_overlap_trafo(ns,ns2,1,pno_S,pnv,pnv1,pnv1,w3,w2,ptr=h1)

              !contract amplitudes in h1 with integrals in w1 and add to w4 : -1 * h1(bdkl) w1(dlkc) += w4(bc)
              call dgemm('n','n',pnv,pnv,pnv1,m10, h1,pnv,w1,pnv1,p10,w5,pnv)

           enddo
        else


           do j=1,rpd1
              do i=1,rpd1
                 p1(i,:,j,:) = p2(idx1(i),:,idx1(j),:)
              enddo
           enddo
           p1 => null()
           p2 => null()

           !transform integral contribution, use symmetry  govov(ldkc) => govov(\bar{d} k l  \bar{c}) to the space of (ij) -> w2
           call dgemm( 'n', 'n', nv*rpd1**2, pnv, nv, p10, w1, nv*rpd1**2, d, nv, nul, w2, nv*rpd1**2 )
           call array_reorder_4d( p10, w2, rpd1, nv, rpd1, pnv, [2,3,1,4], nul, w3 )
           call dgemm( 't', 'n', pnv1, rpd1**2*pnv, nv, p10, d1, nv, w3, nv, nul, w1, pnv1 )

           ! Quadratic part of the E2 term use u^{bd}_{kl} (bkdl) as b,dkl
           call array_reorder_4d( p20, t21, pnv1, rpd1, pnv1, rpd1, [1,3,2,4], nul, w3)
           call array_reorder_4d( m10, t21, pnv1, rpd1, pnv1, rpd1, [1,3,4,2], p10, w3)
           call do_overlap_trafo(ns,ns2,1,pno_S,pnv,rpd1*pnv1*rpd1,pnv1,w3,w2,ptr=h1)


           !contract amplitudes in h1 with integrals in w1 and add to w4 : -1 * h1(bdkl) w1(dlkc) += w4(bc)
           call dgemm('n','n',pnv,pnv,pnv1*rpd1**2,m10, h1,pnv,w1,pnv1*rpd1**2,p10,w5,pnv)


        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  B2 Term !!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!

        !prepare 4 occupied integral goooo for B2 term
        call ass_D1to4( w3(1:rpd**2*rpd1**2), p1, [rpd,rpd,rpd1,rpd1] )
        call ass_D1to4( goooo, p2, [ no, no,  no,  no] )
        !Get the integral contribution, sort it first like the integrals then transform it, govov
        call ass_D1to4( w1,    p3, [rpd1,nv,rpd1,nv] )
        call ass_D1to4( govov, p4, [no,  nv,no,  nv] )


        !CODE FOR PAIR SPACES WITH RESTRICTIONS
        if( PS1 )then

           !loop over pair contributions kl and lk
           do pair = 1, paircontribs

              p3(1,:,1,:) = p4(idx1(paircontrib(1,pair)),:,idx1(paircontrib(2,pair)),:)

              !transform integral contribution
              !govov(kcld) => govov(\bar{c} \bar{d} k l) or lk to the space of (ij) -> w2
              call dgemm( 'n', 'n', nv, pnv, nv, p10, w1, nv, d, nv, nul, w3, nv )
              call dgemm( 't', 'n', pnv, pnv, nv, p10, d, nv, w3, nv, nul, w4, pnv )

              if( PS )then
                 p1(1,1,1,1) = p2(idx1(paircontrib(1,pair)),idx(1),idx1(paircontrib(2,pair)),idx(2))
              else
                 do j = 1,pno
                    do i = 1,pno
                       p1(i,j,1,1) = p2(idx1(paircontrib(1,pair)),idx(i),idx1(paircontrib(2,pair)),idx(j))
                    enddo
                 enddo
              endif

              !sort the amplitudes and contract cidj -> ijcd, ijcd cdkl + ijkl = ijkl
              call array_reorder_4d( p10, t, pnv, rpd, pnv, rpd, [2,4,1,3], nul, w1 )
              call dgemm( 'n', 'n', rpd**2, rpd1**2, pnv**2, p10, w1, rpd**2, w4, pnv**2, p10, w3, rpd**2 )

              !contract the B intermediate in w3 with the amplitudes (kl) from the
              !inner loop and use the overlap to transform to the omega space, ijkl klab
              call array_reorder_2d( p10, t21, pnv1, pnv1 , paircontrib(:,pair), nul, w1 )

              !print *,"pair",pair,"contirb", paircontrib(:,3-pair)
              call dgemm( 'n', 'n', rpd**2, pnv1**2, rpd1**2, p10, w3, rpd**2, w1, rpd1**2, dble(pair-1), w2, rpd**2 )

           enddo

           !CODE FOR RECTANGULAR SPACE IN \sum_kl
        else
           do j=1,pno1
              do i=1,pno1
                 p3(i,:,j,:) = p4(idx1(i),:,idx1(j),:)
              enddo
           enddo

           !transform integral contribution
           !govov(kcld) => govov(\bar{c} \bar{d} k l) to the space of (ij) -> w2
           call dgemm( 'n', 'n', nv*rpd1**2, pnv, nv, p10, w1, nv*rpd1**2, d, nv, nul, w2, nv*rpd1**2 )
           call array_reorder_4d( p10, w2, rpd1, nv, rpd1, pnv, [2,4,1,3], nul, w1 )
           call dgemm( 't', 'n', pnv, rpd1**2*pnv, nv, p10, d, nv, w1, nv, nul, w4, pnv )

           !prepare 4 occupied integral goooo for B2 term
           if( PS )then
              do j=1,pno1
                 do i=1,pno1
                    p1(1,1,i,j) = p2(idx1(i),idx(1),idx1(j),idx(2))
                 enddo
              enddo
           else
              do j=1,pno1
                 do i=1,pno1
                    do b=1,pno
                       do a=1,pno
                          p1(a,b,i,j) = p2(idx1(i),idx(a),idx1(j),idx(b))
                       enddo
                    enddo
                 enddo
              enddo
           endif
           !sort the amplitudes and contract cidj -> ijcd, ijcd cdkl + ijkl = ijkl
           call array_reorder_4d( p10, t, pnv, rpd, pnv, rpd, [2,4,1,3], nul, w1 )
           call dgemm( 'n', 'n', rpd**2, rpd1**2, pnv**2, p10, w1, rpd**2, w4, pnv**2, p10, w3, rpd**2 )

           !contract the B intermediate in w3 with the amplitudes (kl) from the
           !inner loop and use the overlap to transform to the omega space, ijkl klab
           call array_reorder_4d( p10, t21, pnv1, rpd1, pnv1, rpd1, [2,4,1,3], nul, w1 )

           call dgemm( 'n', 'n', rpd**2, pnv1**2, rpd1**2, p10, w3, rpd**2, w1, rpd1**2, nul, w2, rpd**2 )

        endif

        ! transform back, or in the case of ns==ns2 just order correctly and add
        ! the contributions directly to the residual
        if(ns==ns2)then
           call array_reorder_4d( p10, w2, rpd, rpd, pnv, pnv, [3,1,4,2], p10, o )
        else
           call do_overlap_trafo(ns,ns2,2,pno_S,rpd**2*pnv1, pnv, pnv1,w2,w1)
           call array_reorder_4d( p10, w1, rpd, rpd, pnv1, pnv, [3,1,4,2], nul, w3 )
           call do_overlap_trafo(ns,ns2,1,pno_S,pnv,rpd**2*pnv, pnv1,w3,o,pC=p10)
        endif



     enddo FullSpaceLoop


     !Add the E21 contribution
     !CODE FOR PAIR SPACES WITH RESTRICTIONS
     if( PS )then

        !contrib 1
        call dgemm('n','t',pnv,pnv,pnv,p10,t,pnv,w5,pnv,p10,o,pnv)
        !contrib 2
        call dgemm('n','n',pnv,pnv,pnv,p10,w5,pnv,t,pnv,p10,o,pnv)

     else

        call array_reorder_4d( p10, t, pnv, rpd, pnv, rpd, [3,4,1,2], nul, w1)
        call dgemm('n','n',pnv,rpd*pnv*rpd,pnv,p10,w5,pnv,w1,pnv,nul,w2,pnv)
        o = o + w2(1:pnv*rpd*pnv*rpd)
        call array_reorder_4d( p10, w2, pnv, rpd, pnv, rpd, [3,4,1,2], p10, o )

     endif



  end subroutine get_free_summation_for_current_aibj


  subroutine get_common_idx_summation_for_current_aibj(no,ns,pno_cv,pno_S,pno_t2,o2_space,&
        &w1,w2,w3,w4,w5,goovv,govov,Lvoov,oof,p_idx,p_nidx,oidx1,oidx2,nspaces,reference,yep)
     implicit none
     integer, intent(in) :: no,ns,nspaces
     type(PNOSpaceInfo), intent(inout) :: pno_cv(nspaces),pno_S(nspaces*(nspaces-1)/2)
     type(array), intent(in) :: pno_t2(nspaces)
     real(realk),pointer,intent(inout) :: o2_space(:)
     real(realk),pointer,intent(inout) :: w1(:),w2(:),w3(:),w4(:),w5(:)
     real(realk),intent(in) :: goovv(:),govov(:),Lvoov(:),oof(:,:)
     integer,intent(in)    :: p_idx(:,:),p_nidx(:)
     integer,intent(inout) :: oidx1(:,:),oidx2(:,:)
     character :: tr11,tr12,tr21,tr22,TRamp_pos1(2),TRamp_pos2(2),trh1,trh2
     real(realk),pointer :: p1(:,:,:,:), p2(:,:,:,:), p3(:,:,:,:), p4(:,:,:,:)
     real(realk),pointer :: p5(:,:,:,:), p6(:,:,:,:)
     real(realk),pointer :: pij(:,:,:,:), pji(:,:,:,:)
     real(realk),pointer :: h1(:), h2(:), h3(:), h4(:)
     real(realk),pointer :: r1(:,:),r2(:,:),d(:,:),d1(:,:),d2(:,:)
     real(realk),pointer :: o(:),t(:),S1(:,:), t21(:), t22(:), f1(:,:,:)
     logical :: skiptrafo,cyc
     integer :: space,nv,rpd,pno,pnv,nc,nc2,nidx1,nidx2,ns1,ns2,kp(2)
     integer :: a,i,ic,b,j,jc,k,kc,l,lc,c
     integer :: rpd1,pno1,pnv1,Sidx1,ldS1
     integer :: rpd2,pno2,pnv2,Sidx2,ldS2
     integer :: pair1,pair2,paircontrib(2,2)
     integer :: order1(4), ldh1,ldh2,ldh3
     integer :: suborder(2),i_idx,j_idx,pos_in_res,ipos_in_res
     integer :: nidx_h1, nidx_h2
     integer,pointer :: idx(:),idx1(:),idx2(:)
     integer(kind=8) :: o2v2
     integer,parameter :: paircontribs = 2
     real(realk), parameter :: nul =  0.0E0_realk
     real(realk), parameter :: p20 =  2.0E0_realk
     real(realk), parameter :: p10 =  1.0E0_realk
     real(realk), parameter :: p05 =  0.5E0_realk
     real(realk), parameter :: m05 = -0.5E0_realk
     real(realk), parameter :: m10 = -1.0E0_realk
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! DEBUG VARIABLES, REMOVE FOR WORKING CODE
     integer, intent(in) :: yep
     real(realk), intent(in) :: reference(yep,no,yep,no)
     real(realk), pointer :: check_ref(:,:,:,:), phelp(:,:,:,:)
     integer :: diff11,diff12,diff21,diff22, id(2),kcount
     integer :: bpc,epc
     !reduce the no
     logical :: add_contrib,add_contrib1,add_contrib2,FAspace,FAspace1,PS,PS1,PS2
     real(realk) :: pref
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     nv      =  pno_cv(ns)%ns1
     d       => pno_cv(ns)%d
     idx     => pno_cv(ns)%iaos
     pnv     =  pno_cv(ns)%ns2
     pno     =  pno_cv(ns)%n
     rpd     =  pno_cv(ns)%rpd
     FAspace =  pno_cv(ns)%is_FA_space
     PS      =  pno_cv(ns)%PS
     t       => pno_t2(ns)%elm1

     o       => o2_space

     paircontrib(1:2,1) = [1,2]
     paircontrib(1:2,2) = [2,1]

     OneIdxSpaceLoop1: do nc = 1, p_nidx(ns)
        ! extract indices:
        ns1 = p_idx(nc,ns)

        call check_if_contributes(ns,ns1,pno_cv,pno_S,cyc)

        !print *,""
        !print *,""
        !print *,""
        !print *,""
        !print *,"ENTERING LOOP",ns,ns1,cyc,pno,pno_cv(ns1)%n

        if( cyc )then

           cycle OneIdxSpaceLoop1

        endif

        call get_overlap_idx(ns,ns1,pno_cv,oidx1,nidx1,ndidx1=diff11,ndidx2=diff12)

        d1       => pno_cv(ns1)%d
        t21      => pno_t2(ns1)%elm1
        idx1     => pno_cv(ns1)%iaos
        pnv1     =  pno_cv(ns1)%ns2
        pno1     =  pno_cv(ns1)%n
        rpd1     =  pno_cv(ns1)%rpd
        FAspace1 =  pno_cv(ns1)%is_FA_space
        PS1      =  pno_cv(ns1)%PS


        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  C2 Term !!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!

        !Transform the integral contribution to the space of the current amps
        !get the order kjac -> ckja, please note, that this loop is only
        !inside OneIdxSpaceLoop2 because I only use 3 working matrices, p10
        !might easily move the following part outside the loop and add stuff
        !up during the loops
        call ass_D1to4( goovv, p2, [no, no,nv,  nv] )
        if( PS )then

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !CODE FOR TRIANGULAR RESIDUAL CONTRIBUTIONS!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !associate pointers such as to write k last in this case, does not
           !change anything for the pno1==2 case but helps a lot otherwise
           call ass_D1to4( w1,    p1, [nv,rpd,nv,rpd1] )
           !store also the part for Pij
           call ass_D1to4( w1(nv*rpd1*rpd*nv+1:2*nv*rpd1*rpd*nv), p3, [nv,rpd,nv,rpd1] )

           !FIND OUT WHICH CONTRIBUTION TO ADD AND SKIP THE OTHER, FROM THE
           !BEGINNING bpc and epc are set such that no looping occurs, both are
           !then replaced in the loops according to the memreqs, there should be
           !no case where no loop is required, so the choice of begin pair
           !contribuion (bpc) and end pair contribution (epc) should be valid in
           !the end
           bpc = 3
           epc = 2
           do pair1 = 1, paircontribs

              j = idx(paircontrib(2,pair1))

              if( PS1 )then

                 i = idx(paircontrib(1,pair1))

                 if(ns==ns1)then
                    k = j
                 else
                    k = oidx1(nidx1+diff11+1,3)
                 endif

#ifdef VAR_LSDEBUG
                 if(k==i)call lsquit("ERROR(get_common_idx_summation_for_current_aibj): &
                    &this should never happen i==j",-1)
#endif

              else

#ifdef VAR_LSDEBUG
                 if(nidx1/=1)call lsquit("ERROR(get_common_idx_summation_for_current_aibj): &
                    &this should never occur nidx1/=1",-1)
#endif
                 i = oidx1(1,3)
                 k = idx1(1)

              endif
              add_contrib = ((PS1 .and. ((i<j.and.oidx1(1,3)==idx(1)).or.(i>j.and.oidx1(1,3)==idx(2)))  ) .or.ns==ns1 ) &
                 &.or.(.not. PS1 .and. i/=k)

              if(add_contrib     .and.pair1==1)bpc = 1
              if(.not.add_contrib.and.pair1==1)bpc = 2
              if(add_contrib     .and.pair1==2)epc = 2
              if(.not.add_contrib.and.pair1==2)epc = 1
           enddo

           bpc = 1
           epc = 2

           !LOOP OVER THE CONTRIBUTIONS
           do pair1 = bpc, epc

              i = idx(paircontrib(1,pair1))
              j = idx(paircontrib(2,pair1))
              if(pair1==1)then
                 p4 => p1
                 h4 => w4(1:pnv1*rpd1*rpd*pnv)
              else if(pair1==2)then
                 p4 => p3
                 h4 => w4(pnv1*rpd1*rpd*pnv+1:2*pnv1*rpd1*rpd*pnv)
              endif

              if( PS1 ) then

                 !find index k, k=j if ns=ns1 else the non-overlapping index
                 if(ns==ns1)then
                    k = idx1(paircontrib(2,pair1))
#ifdef VAR_LSDEBUG
                    if(k /= j)call lsquit("ERROR():something wrong with &
                       &obtaining the indices",-1)
#endif
                 else
                    k = oidx1(nidx1+diff11+1,3)
                 endif

                 p4(:,1,:,1) = p2(k,j,:,:)

                 ! transform c to \bar(c} of (ki)  and a to \bar{a} of (ij)
                 call dgemm('t','n', pnv, nv, nv,  p10, d, nv, p4, nv, nul, w2, pnv)
                 call dgemm('n','n', pnv, pnv1,nv, p10, w2, pnv, d1, nv, nul, h4, pnv)

              else

                 do kc=1,pno1

                    p4(:,1,:,kc) = p2(idx1(kc),j,:,:)
                    call dgemm('t','n', pnv, nv, nv, p10, d, nv, p4(1,1,1,kc), nv, nul, w2, pnv)
                    call dgemm('n','n', pnv, pnv1,nv, p10, w2, pnv, d1, nv, nul, h4(1+(kc-1)*pnv*pnv1), pnv)

                 enddo

              endif

              !!THE INNER CONTRACTION LOOP - BUILDING THE C INTERMEDIATE
              OneIdxSpaceLoop21: do nc2=1, p_nidx(ns)
                 ! extract indices:
                 ns2 = p_idx(nc2,ns)

                 call check_if_contributes(ns,ns2,pno_cv,pno_S,cyc)

                 if( cyc )then

                    cycle OneIdxSpaceLoop21

                 endif


                 call get_overlap_idx(ns,ns2,pno_cv,oidx2,nidx2,ndidx1=diff21,ndidx2=diff22)

                 d2   => pno_cv(ns2)%d
                 t22  => pno_t2(ns2)%elm1
                 idx2 => pno_cv(ns2)%iaos
                 pnv2 =  pno_cv(ns2)%ns2
                 pno2 =  pno_cv(ns2)%n
                 rpd2 =  pno_cv(ns2)%rpd
                 PS2  =  pno_cv(ns2)%PS


                 if( PS2 )then

                    if(ns==ns2)then
                       l = idx2(paircontrib(1,pair1))
                    else
                       l = oidx2(nidx2+diff21+1,3)
                    endif

                    elseif( rpd2 == 1 )then
                    l = idx2(1)

                    if( l == i ) cycle
                 endif
                 !print *,""
                 !print *,"PAIR CONTRIB LOOP",ns,ns1,ns2,PS1,PS2,i,j,k,l

                 !Get the integrals kdlc -> ckld and transform c and d to their
                 !corresponding spaces, (iajb -> bija) 
                 call ass_D1to4( govov, p2, [no, nv, no, nv]  )
                 call ass_D1to4( w1,    p1, [rpd1,rpd2,nv,nv] )
                 if( PS1 ) then

                    if( PS2 )then

                       !find index k, k=j if ns=ns1 else the non-overlapping index
                       p1(1,1,:,:) = p2(k,:,l,:)

                    else

                       do lc=1,pno2
                          p1(1,lc,:,:) = p2(k,:,idx2(lc),:)
                       enddo

                    endif

                 else

                    if( PS2 ) then

                       do kc=1,pno1
                          p1(kc,1,:,:) = p2(idx1(kc),:,l,:)
                       enddo

                    else

                       do lc=1,pno2
                          do kc=1,pno1
                             p1(kc,lc,:,:) = p2(idx1(kc),:,idx2(lc),:)
                          enddo
                       enddo

                    endif

                 endif

                 ! transform c to \bar(c} in (ki) and d to \bar{d} in (lj)
                 call dgemm('t','t',pnv1,rpd1*rpd2*nv,nv,  p10, d1, nv, w1, rpd1*rpd2*nv, nul, w3, pnv1)
                 call dgemm('n','n',pnv1*rpd1*rpd2,pnv2,nv,p10, w3, pnv1*rpd1*rpd2, d2, nv, nul, w1, pnv1*rpd1*rpd2)

                 !get the amplitudes in the correct order eldj -> ldje transform to a and contract to
                 ! -0.5 w1(ckld) w2(ldja) += w4(ckja)
                 call ass_D1to4( t22, p2, [pnv2,rpd2,pnv2,rpd2] )
                 if( PS2 )then
                    if(l<j)call array_reorder_2d(p10,t22,pnv2,pnv2,[1,2],nul,w3)
                    if(j>l)call array_reorder_2d(p10,t22,pnv2,pnv2,[2,1],nul,w3)
                    nidx_h2 = 1
                 else
                    call ass_D1to4( w3,  p3, [rpd2,pnv2,nidx_h2,pnv2] )
                    !print *,"#CHECK",nidx2,pno2,oidx2(1,2)
                    do b=1,pnv2
                       do jc=1,nidx2
                          do a=1,pnv2
                             do ic=1,pno2
                                p3(ic,b,jc,a) = p2(a,ic,b,oidx2(jc,2))
                             enddo
                          enddo
                       enddo
                    enddo
                    nidx_h2 = nidx2
                 endif

                 call do_overlap_trafo(ns,ns2,2,pno_S,rpd2*pnv2*nidx_h2,pnv,pnv2,w3,w2,ptr=h1,ptr2=h2)

                 !call print_tensor_unfolding_with_labels(w1,&
                 !   &[pnv1,rpd1],'ck',2,[rpd2,pnv2],'ld',2,'inner pair mat 1')
                 !call print_tensor_unfolding_with_labels(h1,&
                 !   &[rpd2,pnv2],'ld',2,[nidx_h2,pnv],'ja',2,'inner pair mat 2')

                 call dgemm('n','n', pnv1*rpd1, nidx_h2*pnv, rpd2*pnv2, m05, w1, pnv1*rpd1, h1, rpd2*pnv2, nul, h2, pnv1*rpd1)

                 !call print_tensor_unfolding_with_labels(h2,&
                 !   &[pnv1,rpd1],'ck',2,[nidx_h2,pnv],'ja',2,'inner pair result')

                 call ass_D1to4( h2, p3, [pnv1,rpd1,nidx_h2,pnv] )
                 call ass_D1to4( h4, p4, [rpd,pnv,pnv1,rpd1] )

                 !call print_tensor_unfolding_with_labels(h4,&
                 !   &[rpd,pnv],'ja',2,[pnv1,rpd1],'ck',2,'pair mat 2 - before add')

                 if( PS1 )then

                    do c=1,pnv1
                       do a=1,pnv
                          p4(1,a,c,1) = p4(1,a,c,1) + p3(c,1,1,a)
                       enddo
                    enddo

                 else

                    do kc=1,rpd1
                       do c=1,pnv1
                          do a=1,pnv
                             p4(1,a,c,kc) = p4(1,a,c,kc) + p3(c,kc,1,a)
                          enddo
                       enddo
                    enddo

                    !do kc=1,rpd1
                    !   do c=1,pnv1
                    !      do a=1,pnv
                    !         do jc=1,nidx_h2
                    !            p4(oidx2(jc,1),a,c,kc) = p4(oidx2(jc,1),a,c,kc) + p3(c,kc,jc,a)
                    !         enddo
                    !      enddo
                    !   enddo
                    !enddo
                 endif

                 !call print_tensor_unfolding_with_labels(h4,&
                 !   &[rpd,pnv],'ja',2,[pnv1,rpd1],'ck',2,'pair mat 2 - after add')

              enddo OneIdxSpaceLoop21

              !print *,""
              !print *,"EXIT PAIR CONTRIB LOOP",pair1


              !get the amplitudes, extract the necessary indices, 
              !reorder dkci -> dick :D transform to current space (bick) and do the contraction,
              !bick ja,ck^T = bija, do the permutation and addition of the contribution

              call ass_D1to4( t21, p6, [pnv1,rpd1,pnv1, rpd1] )
              call ass_D1to4( w1,  p5, [pnv1,1   ,pnv1, rpd1] )

              if( PS1 )then

                 if(k<i)call array_reorder_2d(p10,p6,pnv1,pnv1,[1,2],nul,w1)
                 if(k>i)call array_reorder_2d(p10,p6,pnv1,pnv1,[2,1],nul,w1)

              else

                 do kc=1,pno1
                    p5(:,1,:,kc) = p6(:,kc,:,oidx1(1,2))
                 enddo

              endif


              do kc = 1,rpd1

                 call do_overlap_trafo(ns,ns1,1,pno_S, pnv,pnv1, pnv1,w1(1+(kc-1)*pnv1**2:),w2,ptr=h1,ptr2=h2)

                 !call print_tensor_unfolding_with_labels(h1,&
                 !   &[pnv,rpd],'bi',2,[pnv1,1],'ck',2,'pair mat 1')
                 !call print_tensor_unfolding_with_labels(h4(1+(kc-1)*pnv1*pnv:),&
                 !   &[pnv1,1],'ck',2,[rpd,pnv],'ja',2,'pair mat 2')


                 call dgemm('n','t', pnv, pnv, pnv1, m10, h1, pnv, h4(1+(kc-1)*pnv1*pnv), pnv, nul, h2, pnv)


                 !call print_tensor_unfolding_with_labels(h2,&
                 !   &[pnv],'b',1,[pnv],'a',1,'pair result after adding')


                 !call print_tensor_unfolding_with_labels(o,&
                 !   &[pnv,rpd],'ai',2,[pnv,rpd],'bj',2,'OMEGA AIBJ - Before')

                 add_contrib1 = ((PS1 .and. ((i<j.and.oidx1(1,3)==idx(1)).or.(i>j.and.oidx1(1,3)==idx(2)))  ) .or.ns==ns1 ) &
                    &.or.(.not. PS1 .and. j/=idx1(kc))
                 add_contrib2 = ((PS1 .and. ((i>j.and.oidx1(1,3)==idx(1)).or.(i<j.and.oidx1(1,3)==idx(2)))  ) .or.ns==ns1 ) &
                    &.or.(.not. PS1 .and. i/=idx1(kc))

                 if(add_contrib1)then
                    call array_reorder_2d(p10,h2,pnv,pnv,paircontrib(:,3-pair1),p10,o)
                    call array_reorder_2d(p05,h2,pnv,pnv,paircontrib(:,pair1),p10,o)
                 endif

                 !call print_tensor_unfolding_with_labels(o,&
                 !   &[pnv,rpd],'ai',2,[pnv,rpd],'bj',2,'OMEGA AIBJ')

              enddo


           enddo

        else


           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !CODE FOR RECTANGULAR RESIDUAL CONTRIBUTIONS!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


           call ass_D1to4( w1, p1, [rpd1,rpd,nv,nv] )
           if( PS1 )then
#ifdef VAR_LSDEBUG
              if(diff12 /= 1.or.nidx1/=1) then
                 call lsquit("ERROR(get_common_idx_summation_for_current_aibj): this should never occur",-1)
              endif
#endif
              do j=1,pno
                 p1(1,j,:,:) = p2(oidx1(nidx1+diff11+1,3),idx(j),:,:)
              enddo
           else
              do j=1,pno
                 do k=1,pno1
                    p1(k,j,:,:) = p2(idx1(k),idx(j),:,:)
                 enddo
              enddo
           endif
           ! transform c to \bar(c} of (ki)  and a to \bar{a} of (ij)
           call dgemm('t','t', pnv1, rpd1*rpd*nv, nv,  p10, d1, nv, w1,rpd1*rpd*nv, nul, w2, pnv1)
           call dgemm('n','n', pnv1*rpd1*rpd, pnv, nv, p10, w2, pnv1*rpd1*rpd, d, nv, nul, w1, pnv1*rpd1*rpd)
           call array_reorder_4d(p10,w1,pnv1,rpd1,rpd,pnv,[3,4,1,2],nul,w4)

           !At this point the integrals with transformed virtual indices are
           !restricted to the correct occupied space (i.e. k and j belong to the
           !corresponding amplitude spaces) and saved in the order jack :)


           !!!THE INNER CONTRACTION LOOP - BUILDING THE C INTERMEDIATE
           OneIdxSpaceLoop22: do nc2=1, p_nidx(ns)
              ! extract indices:
              ns2 = p_idx(nc2,ns)

              call check_if_contributes(ns,ns2,pno_cv,pno_S,cyc)

              if( cyc )then

                 cycle OneIdxSpaceLoop22

              endif

              call get_overlap_idx(ns,ns2,pno_cv,oidx2,nidx2,ndidx1=diff21,ndidx2=diff22)

              d2   => pno_cv(ns2)%d
              t22  => pno_t2(ns2)%elm1
              idx2 => pno_cv(ns2)%iaos
              pnv2 =  pno_cv(ns2)%ns2
              pno2 =  pno_cv(ns2)%n
              rpd2 =  pno_cv(ns2)%rpd
              PS2  =  pno_cv(ns2)%PS

              !print *,""
              !print *,"RECT CONTRIB LOOP",ns2,PS1,PS2

              !Get the integrals kdlc -> ckld and transform c and d to their
              !corresponding spaces, (iajb -> bija) 
              call ass_D1to4( govov, p2, [no, nv, no, nv]  )
              call ass_D1to4( w1,    p1, [rpd1,rpd2,nv,nv] )
              if( PS1 ) then
                 if( PS2 )then
                    p1(1,1,:,:) = p2(oidx1(nidx1+diff11+1,3),:,oidx2(nidx2+diff21+1,3),:)
                 else
                    do lc=1,pno2
                       p1(1,lc,:,:) = p2(oidx1(nidx1+diff11+1,3),:,lc,:)
                    enddo
                 endif
              else
                 if( PS2 ) then
                    do k=1,pno1
                       p1(k,1,:,:) = p2(idx1(k),:,oidx2(nidx2+diff21+1,3),:)
                    enddo
                 else
                    do lc=1,pno2
                       do k=1,pno1
                          p1(k,lc,:,:) = p2(idx1(k),:,idx2(lc),:)
                       enddo
                    enddo
                 endif
              endif

              ! transform c to \bar(c} in (ki) and d to \bar{d} in (lj)
              call dgemm('t','t',pnv1,rpd1*rpd2*nv,nv,  p10, d1, nv, w1, rpd1*rpd2*nv, nul, w3, pnv1)
              call dgemm('n','n',pnv1*rpd1*rpd2,pnv2,nv,p10, w3, pnv1*rpd1*rpd2, d2, nv, nul, w1, pnv1*rpd1*rpd2)

              !get the amplitudes in the correct order eldj -> ldje transform to a and contract to
              ! -0.5 w1(ckld) w2(ldja) += w4(ckja)
              call ass_D1to4( w3,  p3, [rpd2,pnv2,nidx2,pnv2] )
              call ass_D1to4( t22, p2, [pnv2,rpd2,pnv2,rpd2] )
              if( PS2 )then
                 if(nidx2/=1)call lsquit("ERRROROROROOROROROR",-1)
                 if( oidx2(nidx2+diff21+1,2) < oidx2(1,2) ) then
                    do b=1,pnv2
                       do a=1,pnv2
                          p3(1,b,1,a) = p2(a,1,b,1)
                       enddo
                    enddo
                 else
                    p3(1,:,1,:) = p2(:,1,:,1)
                 endif
              else
                 do b=1,pnv2
                    do j=1,nidx2
                       do a=1,pnv2
                          do i=1,pno2
                             p3(i,b,j,a) = p2(a,i,b,oidx2(j,2))
                          enddo
                       enddo
                    enddo
                 enddo
              endif

              call do_overlap_trafo(ns,ns2,2,pno_S,rpd2*pnv2*nidx2,pnv,pnv2,w3,w2,ptr=h1,ptr2=h2)
              !call print_tensor_unfolding_with_labels(w1,&
              !   &[pnv1,rpd1],'ck',2,[rpd2,pnv2],'ld',2,'inner rect mat 1')
              !call print_tensor_unfolding_with_labels(h1,&
              !   &[rpd2,pnv2],'ld',2,[nidx2,pnv],'ja',2,'inner rect mat 2')

              call dgemm('n','n', pnv1*rpd1, nidx2*pnv, rpd2*pnv2, m05, w1, pnv1*rpd1, h1, rpd2*pnv2, nul, h2, pnv1*rpd1)

              !call print_tensor_unfolding_with_labels(h2,&
              !   &[pnv1,rpd1],'ck',2,[nidx2,pnv],'ja',2,'inner rect result')

              !call print_tensor_unfolding_with_labels(w4,&
              !   &[rpd,pnv],'ja',2,[pnv1,rpd1],'ck',2,'rect mat 2 - before add')

              call ass_D1to4( h2, p3, [pnv1,rpd1,nidx2,pnv] )
              call ass_D1to4( w4, p4, [rpd,pnv,pnv1,rpd1] )
              if( PS1 )then
                 if( PS2 )then
                    do c=1,pnv1
                       do a=1,pnv
                          p4(oidx2(1,1),a,c,1) = p4(oidx2(1,1),a,c,1) + p3(c,1,1,a)
                       enddo
                    enddo
                 else
                    do c=1,pnv1
                       do a=1,pnv
                          do j=1,nidx2
                             p4(oidx2(j,1),a,c,1) = p4(oidx2(j,1),a,c,1) + p3(c,1,j,a)
                          enddo
                       enddo
                    enddo
                 endif
              else
                 if( PS2 )then
                    do k=1,pno1
                       do c=1,pnv1
                          do a=1,pnv
                             p4(oidx2(1,1),a,c,k) = p4(oidx2(1,1),a,c,k) + p3(c,k,1,a)
                          enddo
                       enddo
                    enddo
                 else
                    do k=1,pno1
                       do c=1,pnv1
                          do a=1,pnv
                             do j=1,nidx2
                                p4(oidx2(j,1),a,c,k) = p4(oidx2(j,1),a,c,k) + p3(c,k,j,a)
                             enddo
                          enddo
                       enddo
                    enddo
                 endif
              endif
              !call print_tensor_unfolding_with_labels(w4,&
              !   &[rpd,pnv],'ja',2,[pnv1,rpd1],'ck',2,'rect mat 2 - after add')

           enddo OneIdxSpaceLoop22

           !print *,""
           !print *,"EXIT RECT CONTRIB LOOP"

           !get the amplitudes, extract the necessary indices, 
           !reorder dkci -> dick :D transform to current space (bick) and do the contraction,
           !bick ja,ck^T = bija, do the permutation and addition of the contribution

           call ass_D1to4( t21, p2, [pnv1,rpd1,pnv1, rpd1] )
           call ass_D1to4( w1,  p1, [pnv1,nidx1,pnv1,rpd1] )
           if( PS1 )then

              k = oidx1(nidx1+diff11+1,3)
              i = oidx1(1,3)

              if(k<i)call array_reorder_2d(p10,p2,pnv1,pnv1,[1,2],nul,p1)
              if(k>i)call array_reorder_2d(p10,p2,pnv1,pnv1,[2,1],nul,p1)

           else

              do k=1,pno1
                 do i=1,nidx1
                    p1(:,i,:,k) = p2(:,k,:,oidx1(i,2))
                 enddo
              enddo

           endif

           call do_overlap_trafo(ns,ns1,1,pno_S, pnv,nidx1*pnv1*rpd1, pnv1,w1,w2,ptr=h1,ptr2=h2)


           !call print_tensor_unfolding_with_labels(h1,&
           !   &[pnv,nidx1],'bi',2,[pnv1,rpd1],'ck',2,'rect mat 1')
           !call print_tensor_unfolding_with_labels(w4,&
           !   &[pnv1,rpd1],'ck',2,[rpd,pnv],'ja',2,'rect mat 2')


           call dgemm('n','t', pnv*nidx1,rpd*pnv, rpd1*pnv1, m10, h1,pnv*nidx1, w4, pnv*rpd, nul, h2, pnv*nidx1)

           !call print_tensor_unfolding_with_labels(h2,&
           !   &[pnv,nidx1],'bi',2,[rpd,pnv],'ja',2,'rect result')

           !call print_tensor_unfolding_with_labels(o,&
           !   &[pnv,rpd],'ai',2,[pnv,rpd],'bj',2,'OMEGA AIBJ - Before')

           call ass_D1to4( h2, p2, [pnv,nidx1,rpd,pnv ] )
           call ass_D1to4( o,  p1, [pnv,rpd, pnv, rpd ] )
           do a=1,pnv
              do j=1,pno
                 do i=1,nidx1
                    do b=1,pnv
                       p1(a,oidx1(i,1),b,j) = p1(a,oidx1(i,1),b,j) + p2(b,i,j,a)
                       p1(a,j,b,oidx1(i,1)) = p1(a,j,b,oidx1(i,1)) + p05 * p2(b,i,j,a)
                       p1(b,j,a,oidx1(i,1)) = p1(b,j,a,oidx1(i,1)) + p2(b,i,j,a)
                       p1(b,oidx1(i,1),a,j) = p1(b,oidx1(i,1),a,j) + p05 * p2(b,i,j,a)
                    enddo
                 enddo
              enddo
           enddo
           !print *,ns,"adding",ns1," as symmetrized",pno,pno1,rpd,rpd1
           !call print_tensor_unfolding_with_labels(o,&
           !   &[pnv,rpd],'ai',2,[pnv,rpd],'bj',2,'OMEGA AIBJ')
        endif


        !if(ns==2 .and. ns1==2)then
        !   print *,"DONE"
        !   stop 0
        !endif



        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  D2 Term !!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!

        !Similar procedure as for the C2 term, just with the L integrals (which
        !could also be produced on-the-fly to reduce the memory requirements
        call ass_D1to4( w1,    p1, [nv,pno,pno1,nv] )
        call ass_D1to4( Lvoov, p2, [nv, no, no, nv] )
        do b=1,nv
           do j=1,pno1
              do i=1,pno
                 do a=1,nv
                    p1(a,i,j,b) = p2(a,idx(i), idx1(j),b)
                 enddo
              enddo
           enddo
        enddo

        ! transform c to \bar(c} of (jk)  and a to \bar{a} of (ij) and reorder
        ! to the in which it will be used later we got w4:\bar{c}k\bar{a}i
        call dgemm('t','n', pnv, pno*pno1*nv, nv,  p10, d, nv, w1, nv, nul, w2, pnv)
        call dgemm('n','n', pnv*pno*pno1, pnv1, nv, p10, w2, pnv*pno1*pno, d1, nv, nul, w1, pnv*pno*pno1)
        call array_reorder_4d( p10, w1, pnv, pno, pno1, pnv1, [4,3,1,2], nul, w4 )

        !THE INNER CONTRACTION LOOP - BUILDING THE D INTERMEDIATE
        OneIdxSpaceLoop3: do nc2=1, p_nidx(ns)
           ! extract indices:
           ns2 = p_idx(nc2,ns)

           call check_if_contributes(ns,ns2,pno_cv,pno_S,cyc)

           if( cyc )then

              cycle OneIdxSpaceLoop3

           endif

           call get_overlap_idx(ns,ns2,pno_cv,oidx2,nidx2)

           d2   => pno_cv(ns2)%d
           t22  => pno_t2(ns2)%elm1
           idx2 => pno_cv(ns2)%iaos
           pnv2 =  pno_cv(ns2)%ns2
           pno2 =  pno_cv(ns2)%n

           !Get the L integrals lfkc -> cklf and transform c and d to their
           !corresponding spaces, (iajb -> bjia) 
           call ass_D1to4( w1,    p1, [nv,pno1,pno2,nv] )
           call ass_D1to4( govov, p2, [no, nv, no,  nv] )
           do a=1,nv
              do i=1,pno2
                 do j=1,pno1
                    do b=1,nv
                       p1(b,j,i,a) = p20 * p2(idx2(i),a,idx1(j),b) - p2(idx1(j),a,idx2(i),b)
                    enddo
                 enddo
              enddo
           enddo

           ! transform c to \bar(c} in (ki) and f to \bar{d} in (lj)
           call dgemm('t','n',pnv1,pno1*pno2*nv,nv,  p10, d1, nv, w1, nv, nul, w3, pnv1)
           call dgemm('n','n',pnv1*pno1*pno2,pnv2,nv,p10, w3, pnv1*pno1*pno2, d2, nv, nul, w1, pnv1*pno1*pno2)

           !get the u amplitudes in the order eifl -> eifl  transform e to a, reorder to lfai and contract to
           ! -0.5 w1(cklf) h1(lfai) += w4(ckai)
           call ass_D1to4( w3,  p3, [pnv2,nidx2,pnv2,pno2] )
           call ass_D1to4( t22, p2, [pnv2,pno2, pnv2,pno2] )
           do j=1,pno2
              do b=1,pnv2
                 do i=1,nidx2
                    do a=1,pnv2
                       p3(a,i,b,j) = p20 * p2(a,oidx2(i,2),b,j) - p2(a,j,b,oidx2(i,2))
                    enddo
                 enddo
              enddo
           enddo

           call do_overlap_trafo(ns,ns2,1,pno_S, pnv, nidx2*pno2*pnv2, pnv2,w3,w2,ptr=h1,ptr2=h2)
           call array_reorder_4d( p10, h1, pnv, nidx2, pnv2, pno2, [4,3,1,2], nul, h2 )
           call dgemm('n','n', pnv1*pno1, pnv*nidx2, pno2*pnv2, p05, w1, pnv1*pno1, h2, pno2*pnv2, nul, h1, pnv1*pno1)


           call ass_D1to4( h1, p2, [pnv1,pno1,pnv,nidx2] )
           call ass_D1to4( w4, p4, [pnv1,pno1,pnv,pno] )
           do j=1,nidx2
              do b=1,pnv
                 do i=1,pno1
                    do a=1,pnv1
                       p4(a,i,b,oidx2(j,1)) = p4(a,i,b,oidx2(j,1)) + p2(a,i,b,j)
                    enddo
                 enddo
              enddo
           enddo

        enddo OneIdxSpaceLoop3

        !exctract amplitudes as u bjck and contract with w4 ckai
        call ass_D1to4( w1,  p1, [pnv1,nidx1,pnv1,pno1] )
        call ass_D1to4( t21, p2, [pnv1,pno1,pnv1,pno1] )
        do j=1,pno1
           do b=1,pnv1
              do i=1,nidx1
                 do a=1,pnv1
                    p1(a,i,b,j) = p20 * p2(a,oidx1(i,2),b,j) - p2(a,j,b,oidx1(i,2))
                 enddo
              enddo
           enddo
        enddo

        call do_overlap_trafo(ns,ns1,1,pno_S,pnv,nidx1*pnv1*pno1,pnv1,w1,w2,ptr=h1,ptr2=h2)
        call dgemm('n','n',pnv*nidx1,pnv*pno,pnv1*pno1,p05,h1,pnv*nidx1,w4,pnv1*pno1,nul,h2,pnv*nidx1)

        !add D2 contribution to o
        call ass_D1to4( h2, p2, [pnv,nidx1,pnv,pno] )
        call ass_D1to4( o,  p1, [pnv,pno, pnv, pno ] )
        do a=1,pnv
           do i=1,pno
              do j=1,nidx1
                 do b=1,pnv
                    p1(a,i,b,oidx1(j,1)) = p1(a,i,b,oidx1(j,1)) + p2(b,j,a,i)
                    p1(b,oidx1(j,1),a,i) = p1(b,oidx1(j,1),a,i) + p2(b,j,a,i)
                 enddo
              enddo
           enddo
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  E2 Term part 2!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !Similar procedure as for the C2 term, just nothing has to be
        !transformed for the occ-occ fock matrix --> might be constructed
        !outside the loops and saved along with the overlaps, the Foo( k, j )
        call ass_D1to2( w4,  r1, [ pno1,pno ] )
        do j=1,pno
           do i=1,pno1
              r1(i,j) = oof(idx1(i), idx(j))
           enddo
        enddo
        r1 => null()

        !THE INNER CONTRACTION LOOP - BUILDING THE E22 INTERMEDIATE
        OneIdxSpaceLoop4: do nc2=1, p_nidx(ns)
           ! extract indices:
           ns2 = p_idx(nc2,ns)

           call check_if_contributes(ns,ns2,pno_cv,pno_S,cyc)

           if( cyc )then

              cycle OneIdxSpaceLoop4

           endif

           call get_overlap_idx(ns,ns2,pno_cv,oidx2,nidx2)

           d2   => pno_cv(ns2)%d
           t22  => pno_t2(ns2)%elm1
           idx2 => pno_cv(ns2)%iaos
           pnv2 =  pno_cv(ns2)%ns2
           pno2 =  pno_cv(ns2)%n

           !Get the integrals g(kdlc) as (dklc) and transform c and d to (lj)
           !such that the order klcd is obtained
           call ass_D1to4( w1,    p1, [nv,pno1,pno2,nv] )
           call ass_D1to4( govov, p2, [no, nv, no,  nv] )
           do b=1,nv
              do j=1,pno2
                 do i=1,pno1
                    do a=1,nv
                       p1(a,i,j,b) = p2(idx1(i),a,idx2(j),b)
                    enddo
                 enddo
              enddo
           enddo
           p1 => null()
           p2 => null()
           ! transform c to \bar(c} in (lj) and d to \bar{d} in (lj)
           call dgemm('n','n',nv*pno1*pno2,pnv2,nv,  p10, w1, nv*pno1*pno2,d2,nv,nul, w3, nv*pno1*pno2)
           call dgemm('t','n',pno1*pno2*pnv2,pnv2,nv,p10, w3, nv, d2, nv, nul, w1, pno1*pno2*pnv2)

           !get the u amplitudes in the order cldj -> (lcdj) = 2 t(lcdj) - t(jcdl)
           call ass_D1to4( w3,  p3, [pno2,pnv2,pnv2,nidx2] )
           call ass_D1to4( t22, p2, [pnv2,pno2, pnv2,pno2] )
           do j=1,nidx2
              do b=1,pnv2
                 do a=1,pnv2
                    do i=1,pno2
                       p3(i,a,b,j) = p20 * p2(a,i,b,oidx2(j,2)) - p2(a,oidx2(j,2),b,i)
                    enddo
                 enddo
              enddo
           enddo

           call dgemm('n','n', pno1, nidx2, pno2*pnv2**2, p10, w1, pno1, w3, pno2*pnv2**2, nul, w2, pno1 )

           call ass_D1to2( w2, r2, [pno1,nidx2] )
           call ass_D1to2( w4, r1, [pno1,pno] )
           do j=1,nidx2
              do i=1,pno1
                 r1(i,oidx2(j,1)) = r1(i,oidx2(j,1)) + r2(i,j)
              enddo
           enddo

        enddo OneIdxSpaceLoop4

        !extract amplitudes like in C2 as aibk
        call ass_D1to4( w1,  p1, [pnv1,nidx1,pnv1,pno1] )
        call ass_D1to4( t21, p2, [pnv1,pno1,pnv1,pno1] )
        do j=1,pno1
           do b=1,pnv1
              do i=1,nidx1
                 do a=1,pnv1
                    p1(a,i,b,j) = p2(a,oidx1(i,2),b,j)
                 enddo
              enddo
           enddo
        enddo

        call do_overlap_trafo(ns,ns1,1,pno_S, pnv,nidx1*pnv1*pno1, pnv1 ,w1,w2,ptr=h1,ptr2=h2)

        call dgemm('n','n',pnv*nidx1*pnv1,pno,pno1,m10,h1,pnv*nidx1*pnv1,w4,pno1,nul,h2,pnv*nidx1*pnv1)
        call array_reorder_4d(p10,h2,pnv,nidx1,pnv1,pno,[3,4,1,2], nul, h1)

        !transform b index to the correct space
        call do_overlap_trafo(ns,ns1,1,pno_S, pnv,pno*pnv*nidx1,pnv1,h1,h2,ptr=h1)

        call ass_D1to4( h1, p2, [pnv,pno,pnv,nidx1] )
        call ass_D1to4( o,  p1, [pnv,pno, pnv, pno] )
        do i=1,nidx1
           do a=1,pnv
              do j=1,pno
                 do b=1,pnv
                    p1(a,oidx1(i,1),b,j) = p1(a,oidx1(i,1),b,j) + p2(b,j,a,i)
                    p1(b,j,a,oidx1(i,1)) = p1(b,j,a,oidx1(i,1)) + p2(b,j,a,i)
                 enddo
              enddo
           enddo
        enddo

     enddo OneIdxSpaceLoop1

     !contribution done, check the individual elements
     !call ass_D1to4(o,check_ref,[pnv,rpd,pnv,rpd])
     !if( PS )then
     !   j=2
     !   i=1
     !   do b = 1, pnv
     !      do a = 1, pnv
     !         if( abs(check_ref(a,1,b,1) - reference(a,idx(i),b,idx(j))) > 1.0E-12 )then
     !            print *,"wrong element in PNO space",ns,"el",a,b
     !            print *,"element",check_ref(a,1,b,1),reference(a,idx(i),b,idx(j)),reference(a,idx(j),b,idx(i))
     !            stop 0
     !         endif
     !      enddo
     !   enddo
     !else
     !   do j = 1, pno
     !      do b = 1, pnv
     !         do i = 1, pno
     !            do a = 1, pnv
     !               if( abs(check_ref(a,i,b,j) - reference(a,idx(i),b,idx(j))) > 1.0E-12 )then
     !                  print *,"wrong element in rectangular space",a,i,b,j,check_ref(a,i,b,j),reference(a,idx(i),b,idx(j))
     !                  stop 0
     !               endif
     !            enddo
     !         enddo
     !      enddo
     !   enddo
     !endif

  end subroutine get_common_idx_summation_for_current_aibj


  subroutine successive_4ao_mo_trafo(ao,WXYZ,WW,w,XX,x,YY,y,ZZ,z,WRKWXYZ)
     implicit none
     integer, intent(in) :: ao,w,x,y,z
     real(realk), intent(inout) :: WXYZ(ao*ao*ao*ao),WRKWXYZ(ao*ao*ao*w)
     real(realk), intent(in) :: WW(ao,w),XX(ao,x),YY(ao,y),ZZ(ao,z)
     !WXYZ(ao,ao ao ao)^T WW(ao,w)   -> WRKWXYZ (ao ao ao,w)
     call dgemm('t','n',ao*ao*ao,w,ao,1.0E0_realk,WXYZ,ao,WW,ao,0.0E0_realk,WRKWXYZ,ao*ao*ao)
     ! WRKWXYZ(ao,ao ao w)^T XX(ao,x)   -> WXYZ (ao ao w, x)
     call dgemm('t','n',ao*ao*w,x,ao,1.0E0_realk,WRKWXYZ,ao,XX,ao,0.0E0_realk,WXYZ,ao*ao*w)
     ! WXYZ(ao, ao w x)^T YY(ao,y)   -> WRKYXYX (ao w x,y)
     call dgemm('t','n',ao*w*x,y,ao,1.0E0_realk,WXYZ,ao,YY,ao,0.0E0_realk,WRKWXYZ,ao*w*x)
     ! WRKWXYZ(ao, w x y)^T ZZ(ao,z)^T   -> WXYZ (wxyz)
     call dgemm('t','n',w*x*y,z,ao,1.0E0_realk,WRKWXYZ,ao,ZZ,ao,0.0E0_realk,WXYZ,w*x*y)
  end subroutine successive_4ao_mo_trafo

  subroutine free_PNOSpaceInfo(SPINFO)
     implicit none
     type(PNOSpaceInfo), intent(inout) :: SPINFO

     if(.not.SPINFO%allocd)then
        call lsquit("ERROR(free_PNOSpaceInfo): structure not allocated, cannot free",-1)
     endif

     if(associated(SPINFO%iaos)) call mem_dealloc(SPINFO%iaos)
     if(associated(SPINFO%d   )) call mem_dealloc(SPINFO%d   )
     if(associated(SPINFO%s1  )) call mem_dealloc(SPINFO%s1  )
     if(associated(SPINFO%s2  )) call mem_dealloc(SPINFO%s2  )
     SPINFO%n    = 0
     SPINFO%ns1  = 0
     SPINFO%ns2  = 0
     SPINFO%rpd  = 0
     SPINFO%red1 = 0
     SPINFO%red2 = 0

     SPINFO%allocd = .false.

  end subroutine  free_PNOSpaceInfo


  subroutine print_tensor_unfolding_with_labels(mat,d1,l1,m1,d2,l2,m2,label)
     implicit none
     integer, intent(in) :: m1,m2
     integer, intent(in) :: d1(m1),d2(m2)
     character,intent(in) :: l1(m1),l2(m2)
     real(realk), intent(in) :: mat(*)
     character*(*), intent(in) :: label
     integer :: i,a,b,id1(m1),id2(m2),cd1,cd2

     print *,label(1:len(label))
     write (*,'("(")',advance='no') 
     cd1 = 1
     do i=1,m1
        if(i>1)write (*,'(X)',advance='no')
        write (*,'(A)',advance='no') l1(i)
        cd1 = cd1 * d1(i)
     enddo
     write (*,'(")  /  (")',advance='no') 
     cd2 = 1
     do i=1,m2
        if(i>1)write (*,'(X)',advance='no')
        write (*,'(A)',advance='no') l2(i)
        cd2 = cd2 * d2(i)
     enddo
     write (*,'(")  ")',advance='no') 

     do b = 1, cd2
        call get_midx(b,id2,d2,m2)
        write (*,'(12X,"(")',advance='no')
        do i = 1, m2
           if(i>1)write (*,'("/")',advance='no')
           write (*,'(I4)',advance='no') id2(i)
        enddo
        write (*,'(")")',advance='no')
     enddo
     write(*,*)

     do a = 1, cd1
        call get_midx(a,id1,d1,m1)
        write (*,'(" (")',advance='no') 
        do i = 1, m1
           if(i>1)write (*,'("/")',advance='no')
           write (*,'(I4)',advance='no') id1(i)
        enddo
        write (*,'(") ")',advance='no') 
        do b = 1, cd2
           write (*,'(g18.5)',advance='no') mat(a+(b-1)*cd1)
           if(b<cd2)write (*,'(",",X)',advance='no')
        enddo
        if(a<cd1)then
           write (*,'(";")')
        else
           write (*,*)
        endif
     enddo



  end subroutine print_tensor_unfolding_with_labels


end module pno_ccsd_module
