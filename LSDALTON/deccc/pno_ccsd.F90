!@file: this file contains the pno-ccsd residual routine and all the relevant
!PNO-specific routines
!author: Patrick Ettenhuber

module pno_ccsd_module
  use,intrinsic :: iso_c_binding,only:c_ptr,c_f_pointer,c_loc

  use precision
  use typedef
  use typedeftype
  use dec_typedef_module
  use screen_mod
  use tensor_interface_module
  use IntegralInterfaceDEC
  
  
  use cc_tools_module
  use dec_fragment_utils
  
  public :: get_ccsd_residual_pno_style, &
     & get_pno_trafo_matrices,      & 
     & get_pno_overlap_matrices,    &
     & successive_4ao_mo_trafo, free_PNOSpaceInfo
  
  private
  
  !SET OMP INFORMATION LOOP
  integer, parameter :: beg_array   = 1
  integer, parameter :: end_array   = 2
  integer, parameter :: nloops      = 3
  integer, parameter :: first_loop  = 1
  integer, parameter :: second_loop = 2
  integer, parameter :: third_loop  = 3
  integer, parameter :: w1_tag      = 1
  integer, parameter :: w2_tag      = 2
  integer, parameter :: w3_tag      = 3
  integer, parameter :: w4_tag      = 4
  integer, parameter :: w5_tag      = 5
  integer :: nthreads_level1_int_dir(nloops), max_nthr_int_loop
  integer,pointer :: info_omp1(:,:,:,:)


  contains


  
  !> \author Patrick Ettenhuber
  !> \date November 2013
  !> \brief this subroutine calculates the ccsd residual by transforming each
  !>        doubles amplitudes to their respective set of PNO's and then
  !>        transforming the result vector back to the reference basis. 
  subroutine get_ccsd_residual_pno_style(t1,t2,o1,o2,govov,no,nv,nb,xo,xv,yo,yv,&
        &mylsitem,fj,pno_cv,pno_S,nspaces,fock,iter,f)
     implicit none
     !ARGUMENTS
     integer, intent(in) :: no, nv, nb,iter,nspaces
     logical, intent(in) :: fj
     real(realk), intent(inout) :: t1(nv,no), t2(nv,no,nv,no)
     real(realk), intent(inout),target :: o1(nv,no)
     real(realk), intent(inout) :: o2(nv,no,nv,no), govov(no*nv*no*nv)
     real(realk), intent(in),target :: xo(nb,no), xv(nb,nv), yo(nb,no), yv(nb,nv),fock(nb,nb)
     type(lsitem), intent(inout) :: mylsitem
     type(decfrag),intent(in),optional :: f
     type(PNOSpaceInfo),intent(inout) :: pno_cv(nspaces)
     type(PNOSpaceInfo),intent(inout) :: pno_S(nspaces*(nspaces-1)/2)
     !INTERNAL VARIABLES
#ifdef MOD_UNRELEASED
     type(tensor),pointer :: pno_o2(:),pno_t2(:),pno_gvvvv(:)
     integer :: ns,c,nc,nc2
     integer(kind=8)     :: s1,   s2,   s3,    s4,   s5
     real(realk),pointer :: w1(:),w2(:),w3(:), w4(:),w5(:)
     real(realk),pointer :: p1(:,:,:,:), p2(:,:,:,:), p3(:,:,:,:), p4(:,:,:,:)
     real(realk),pointer :: h1(:), h2(:), h3(:), h4(:)
     real(realk),pointer :: r1(:,:),r2(:,:)
     real(realk),pointer :: gvvvv(:), gvovo(:), goovv(:), goovv_vvoo(:), gooov(:)
     real(realk),pointer :: Gai(:,:), Gkj(:,:)
     integer :: i, j, a, b, i_idx, la, lg, fa, fg, xa, xg
     integer(kind=8) :: o2v2
     character(tensor_MSG_LEN) :: msg
     real(realk),pointer :: d(:,:), d1(:,:), d2(:,:),t(:), t22(:), t21(:), o(:),vof(:),ovf(:)
     real(realk),pointer :: Lvoov(:)
     integer, pointer :: idx(:),idx1(:),idx2(:), p_idx(:,:), p_nidx(:), oidx1(:,:),oidx2(:,:)
     integer, pointer :: s_idx(:,:,:), s_nidx(:)
     integer :: pno,pno1,pno2,pnv,pnv1,pnv2,k,l,nidx1,nidx2,spacemax,maxvirt,maxocc,bpc,epc,ic,jc
     integer :: max_pnor,max_pnvr,beg1,beg2
     logical :: cyc,use_triangular,PS
     real(realk), pointer :: iFock(:,:), Dens(:,:)
     integer(kind=8) :: maxsize, myload
     integer :: pair,paircontribs,paircontrib(2,2),rpd
     type(tensor), pointer :: sio4(:)
     logical :: master, add_contrib
     type(pno_query_info) :: query

     !Integral stuff
     real(realk), pointer :: oof(:,:),vvf(:,:)
     integer :: alphaB,gammaB,dimAlpha,dimGamma
     integer :: dim1,dim2,dim3,MinAObatch
     integer :: iorb,nthreads
     type(int_batch) :: a_batch, g_batch
     real(realk) :: MemFree,tw,tc,tinit,tamps,tome,tbatchc,tint_dir,tfock,treord,trest,tfin
     real(realk) :: mo_time(6)
     !real(realk) :: ref(no*nv*nv*no), ref1(no*nv), u(nv,no,nv,no)
     real(realk), parameter :: p20 = 2.0E0_realk
     real(realk), parameter :: p10 = 1.0E0_realk
     real(realk), parameter :: m10 = -1.0E0_realk
     real(realk), parameter :: m05 = -0.5E0_realk
     real(realk), parameter :: p05 = 0.5E0_realk
     real(realk), parameter :: nul = 0.0E0_realk
#ifdef VAR_OMP
     logical :: nested
     logical, external :: omp_get_nested

     nested = omp_get_nested()
     call omp_set_nested(.true.)
#endif
  
     call time_start_phase(PHASE_WORK, twall = tw)
     tinit   = tw
     mo_time = 0.0E0_realk
  
     o2v2 = (i8*no**2)*nv**2
     use_triangular = .true.
     master         = .true.

     paircontribs = 2
     paircontrib(1:2,1) = [1,2]
     paircontrib(1:2,2) = [2,1]

     if(fj.and..not.present(f))call lsquit("ERROR(get_ccsd_residual_pno_style):wrong input fj without f",-1)

     ! if  we have a fragment job the basis of the atomic (pair) site has to be
     ! treated in a special way, either LO or FNO basis, this will be element 1
     ! in all the array arrays
     maxocc = 2
     if(present(f))maxocc = max(f%noccEOS,maxocc)
     maxvirt = 0
     do ns = 1, nspaces
        maxvirt = max(pno_cv(ns)%ns2,maxvirt)
     enddo

     call mem_alloc( pno_t2, nspaces  )
     call mem_alloc( pno_o2, nspaces  )
     call mem_alloc( sio4,   nspaces  )
     call mem_alloc( goovv,  o2v2     )
     call mem_alloc( Lvoov,  o2v2     )
     call mem_alloc( gooov,  no**3*nv )
     call mem_alloc( p_nidx, nspaces  )
     call mem_alloc( p_idx,  nspaces  , nspaces )
     call mem_alloc( s_nidx, no )
     call mem_alloc( s_idx,  2        , nspaces , no )
     call mem_alloc( ovf,    no*nv    )
     call mem_alloc( vof,    nv*no    )

     call mem_alloc( Gai,    nb       , no      )


     !zero the relevant quatnities
     !$OMP WORKSHARE
     Gai   = 0.0E0_realk
     goovv = 0.0E0_realk
     Lvoov = 0.0E0_realk
     gooov = 0.0E0_realk
     !$OMP END WORKSHARE

     if(DECinfo%PL>2)then
        call time_start_phase(PHASE_WORK, twall = tamps, ttot = tinit, labelttot =&
           & 'PNO: init and zeroing                :' )
     endif

     !Get all the pno amplitudes with index restrictions i<=j
     call get_pno_amplitudes(t2,pno_cv,pno_t2,nspaces,no,nv)

     if(DECinfo%PL>2)then
        call time_start_phase(PHASE_WORK, twall = tome, ttot = tamps, labelttot = &
           & 'PNO: extract and transform amplitudes:' )
     endif

     !initialize the pno_residual and the sio4 according to the allocated pno_cv
     call init_pno_residual_and_sio4(pno_cv,pno_o2,sio4,nspaces,no)

     if(DECinfo%PL>2)then
        call time_start_phase(PHASE_WORK, twall = tbatchc, ttot = tome, labelttot = &
           & 'PNO: initialize residual and sio4    :' )
     endif

     !call II_get_AbsoluteValueOcc_overlap(DECinfo%output,DECinfo%output,setting,nb,no,out)


     !==================================================
     !                  Batch construction             !
     !==================================================

     call get_batch_sizes_pno_residual(mylsitem,no,nv,nb,maxocc,maxvirt,nspaces,a_batch,&
       &g_batch,sio4,pno_cv,pno_t2,pno_o2,query)


     ! ************************************************
     ! *  Allocate matrices used in the batched loop  *
     ! ************************************************
     s1 = query%size_array(1) 
     call mem_alloc( w1, s1 )
     s2 = query%size_array(2)
     call mem_alloc( w2, s2 )
     s3 = query%size_array(3) 
     call mem_alloc( w3, s3 )
     s4 = query%size_array(4) 
     call mem_alloc( w4, s4 )
     s5 = query%size_array(5)
     call mem_alloc( w5, s5 )

     call free_query_info(query)

     if(DECinfo%PL>2)then
        call time_start_phase(PHASE_WORK, twall = tint_dir, ttot = tbatchc, labelttot = &
           & 'PNO: build batches                   :' )
     endif

     ! Do the batched interal loop
     call pno_residual_integral_direct_loop(mylsitem,w1,s1,w2,s2,w3,s3,w4,s4,w5,s5,no,nv,nb,&
        &maxocc,maxvirt,nspaces,a_batch,g_batch,sio4,pno_cv,pno_t2,pno_o2,xo,xv,yo,yv,gooov,&
        &goovv,govov,Lvoov,Gai)

     if(DECinfo%PL>2)then
        call time_start_phase(PHASE_WORK, twall = tfock, ttot = tint_dir, labelttot = &
           & 'PNO: integral direct loop            :' )
     endif

     ! Free gamma stuff
     call free_batch_info(g_batch)
     ! Free alpha stuff
     call free_batch_info(a_batch)

     call mem_dealloc( w1 )
     call mem_dealloc( w2 )
     call mem_dealloc( w3 )
     call mem_dealloc( w4 )
     call mem_dealloc( w5 )



     !SWITCH TO MO PART
     s1 = i8*nb*max(nv,no)
     call mem_alloc( w1, s1 )

     !!!!!!!!!!!!!!!!!!!
     !GET FOCK MATRICES!
     !!!!!!!!!!!!!!!!!!!
     !allocate the density matrix
     call mem_alloc( oof, no, no )
     call mem_alloc( vvf, nv, nv )

     !Transform inactive Fock matrix into the different mo subspaces
     ! -> Foo
     call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,fock,nb,0.0E0_realk,w1,no)
     call dgemm('n','n',no,no,nb,1.0E0_realk,w1,no,yo,nb,0.0E0_realk,oof,no)
     ! -> Fov
     call dgemm('n','n',no,nv,nb,1.0E0_realk,w1,no,yv,nb,0.0E0_realk,ovf,no)
     ! -> Fvo
     call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,fock,nb,0.0E0_realk,w1,nv)
     call dgemm('n','n',nv,no,nb,1.0E0_realk,w1,nv,yo,nb,0.0E0_realk,vof,nv)
     ! -> Fvv
     call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1,nv,yv,nb,0.0E0_realk,vvf,nv)

     call mem_dealloc( w1 )

     if(DECinfo%PL>2)then
        call time_start_phase(PHASE_WORK, twall = treord, ttot = tfock, labelttot = &
           & 'PNO: fock matrix construction        :' )
     endif

     Lvoov = p20 * Lvoov
     call array_reorder_4d( m10, goovv, no, no ,nv, nv, [3,2,1,4], p10, Lvoov)
     call mem_alloc( goovv_vvoo,  o2v2 )
     call array_reorder_4d( p10, goovv, no, no ,nv, nv, [3,4,1,2], nul, goovv_vvoo)
     call mem_dealloc( goovv )

     if(DECinfo%PL>2)then
        call time_start_phase(PHASE_WORK, twall = trest, ttot = treord, labelttot = &
           & 'PNO: sorting integrals for MO part   :' )
     endif

     call mem_alloc(Gkj, no, no)
     Gkj = oof

     ! D1 term
#ifdef VAR_PTR_RESHAPE
     h1(1:(i8*nv)*no) => o1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
     call c_f_pointer(c_loc(o1(1,1)),h1,[(i8*nv)*no])
#else
     call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

     h1 = vof
     h1 => null()
     ! A1 term
     call dgemm('t','n',nv,no,nb,p10,xv,nb,Gai,nb,p10,o1, nv)
     call dgemm('t','n',no,no,nb,p10,xo,nb,Gai,nb,p10,Gkj,no)
     call mem_dealloc( Gai )

     !Get pair interaction space information
     call get_pair_space_info(pno_cv,p_idx,p_nidx,s_idx,s_nidx,nspaces,no)

     spacemax = maxocc + 2 + 2

     maxsize  = max(i8*no,i8*nv)**4

     call mem_TurnONThread_Memory()
     !$OMP PARALLEL DEFAULT(NONE) PRIVATE(d,t,idx,pnv,pno,a,i,b,j,ns,pnv1,pnv2,pno1,pno2,&
     !$OMP d1,d2,t21,t22,w1,w2,w3,w4,w5,o,idx1,idx2,p1,p2,p3,p4,h1,h2,&
     !$OMP oidx1,nidx1,oidx2,nidx2,i_idx,r1,r2,cyc,& 
     !$OMP nc,nc2,rpd,PS,ic,jc,add_contrib,k,pair,l,bpc,epc) SHARED(pno_cv,pno_s,pno_t2,gvovo,goovv_vvoo,gvvvv,&
     !$OMP vvf,Lvoov,pno_o2,govov,paircontrib,paircontribs,&
     !$OMP Gkj, maxsize, nspaces, ovf,  s_idx,o1,sio4,&
     !$OMP s_nidx,gooov, no, nv, p_idx, p_nidx,spacemax) REDUCTION(+:mo_time)
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

        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  A2 Term !!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!

        !A2.1
        !Get the integral contribution, sort it first like the integrals then transform it
        !call extract_from_gvovo_transform_add(gvovo,w1,w2,d,o,idx,rpd,pno,no,pnv,nv,&
        !   &pno_cv(ns)%n,PS)

        !A2.2
        !call add_A22_contribution_simple(gvvvv,w1,w2,w3,d,t,o,pno,no,pnv,nv,pno_cv(ns)%n,PS)

        call time_start_phase(PHASE_WORK, twall = mo_time(1) )


        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  E2 Term part1, B2 !! 
        !!!!!!!!!!!!!!!!!!!!!!!!!

        call get_free_summation_for_current_aibj(no,ns,pno_cv,pno_S,pno_t2,o,&
           &w1,w2,w3,w4,w5,govov,vvf,sio4,nspaces)


        call time_start_phase(PHASE_WORK, twall = mo_time(2), ttot = mo_time(1) )
        mo_time(3) = mo_time(3) + mo_time(1)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  E2 Term part2, C2, D2 !! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Loop only over the indices which have a common index with the current
        !pair index, this could in principle also be solved with if statements
        !in the previous full loop and only doing the following work in a subset
        call get_common_idx_summation_for_current_aibj(no,ns,pno_cv,pno_S,pno_t2,&
           &o,w1,w2,w3,w4,w5,goovv_vvoo,govov,Lvoov,Gkj,p_idx,p_nidx,oidx1,oidx2,nspaces)


        call time_start_phase(PHASE_WORK, twall = mo_time(1), ttot = mo_time(2) )
        mo_time(4) = mo_time(4) + mo_time(2)

        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  B1 Term !!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !get gooov(kilc) as klic and transform c to pno basis
#ifdef VAR_PTR_RESHAPE
        p3(1:rpd,1:rpd,1:no,1:nv) => w3
        p2(1:no,1:no,1:no,1:nv) => gooov
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
        call c_f_pointer( c_loc(w3(1)),    p3, [rpd,rpd,no, nv] )
        call c_f_pointer( c_loc(gooov(1)), p2, [no, no, no, nv] )
#else
        call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

        if( PS )then


           do pair = 1, 2

              k = idx(paircontrib(1,pair))
              l = idx(paircontrib(2,pair))

              p3(1,1,:,:) = p2(k, :, l ,:)

              ! transform c such that d(c\bar{c})^T w(kli,c)^T = w1(\bar{c}kli)
              call dgemm('t','t',pnv,rpd**2*no,nv,p10,d,nv,w3,rpd**2*no,nul,w1,pnv)

              !get u from current amplitudes contract and transform back to local basis, as u(\bar{a}\bar{c}kl) from akcl
              if( k > l )then
                 call array_reorder_2d( p20, t, pnv,pnv, [2,1], nul,w2)
                 call array_reorder_2d( m10, t, pnv,pnv, [1,2], p10,w2)
              else
                 call array_reorder_2d( p20, t, pnv,pnv, [1,2], nul,w2)
                 call array_reorder_2d( m10, t, pnv,pnv, [2,1], p10,w2)
              endif

              ! carry out w2(\bar{a}\bar{c} kl) w1(\bar{c} kl i) = omega1{\bar{a}i}
              call dgemm('n','n',pnv,no,pnv*rpd**2,m10,w2,pnv,w1,pnv*rpd**2, nul, w3,pnv)
              !transform d(a\bar{a}) omega1{\bar{a} i} -> o1(a,i)
              !$OMP CRITICAL
              call dgemm('n','n',nv, no,pnv, p10,d, nv, w3, pnv, p10,o1,nv)
              !$OMP END CRITICAL
           enddo

        else

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
           call dgemm('t','t',pnv,rpd**2*no,nv,p10,d,nv,w3,rpd**2*no,nul,w1,pnv)

           !get u from current amplitudes contract and transform back to local basis, as u(\bar{a}\bar{c}kl) from akcl
           call array_reorder_4d( p20, t, pnv,rpd,pnv,rpd, [1,3,2,4], nul,w2)
           call array_reorder_4d( m10, t, pnv,rpd,pnv,rpd, [1,3,4,2], p10,w2)


           ! carry out w2(\bar{a}\bar{c} kl) w1(\bar{c} kl i) = omega1{\bar{a}i}
           call dgemm('n','n',pnv,no,pnv*rpd**2,m10,w2,pnv,w1,pnv*rpd**2, nul, w3,pnv)
           !transform d(a\bar{a}) omega1{\bar{a} i} -> o1(a,i)
           !$OMP CRITICAL
           call dgemm('n','n',nv, no,pnv, p10,d, nv, w3, pnv, p10,o1,nv)
           !$OMP END CRITICAL

        endif

        call time_start_phase(PHASE_WORK, twall = mo_time(2), ttot = mo_time(1) )
        mo_time(5) = mo_time(5) + mo_time(1)

        d   => null()
        t   => null()
        idx => null()
        o   => null()
        pnv =  0
        pno =  0 
     enddo LoopContribs
     !$OMP END DO NOWAIT

     call time_start_phase(PHASE_WORK, twall = mo_time(1) )

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
           rpd =  pno_cv(ns)%rpd
           PS  =  pno_cv(ns)%PS



            !A1 has been commented out to avoid the explicit calculation of
            !gvvov
!           !!!!!!!!!!!!!!!!!!!!!!!!!
!           !!!  A1 Term !!!!!!!!!!!!
!           !!!!!!!!!!!!!!!!!!!!!!!!!
!
!           if( PS )then
!              bpc = 3
!              epc = 2
!              do pair = 1, paircontribs
!
!                 add_contrib = .false.
!
!                 k = idx(paircontrib(1,pair))
!                 i = nc
!
!                 add_contrib = (k<i.and.i==idx(2)).or.(k>i.and.i==idx(1))  
!
!                 if(add_contrib     .and.pair==1)bpc = 1
!                 if(.not.add_contrib.and.pair==1)bpc = 2
!                 if(add_contrib     .and.pair==2)epc = 2
!                 if(.not.add_contrib.and.pair==2)epc = 1
!
!              enddo
!
!              do pair = bpc, epc
!
!                 k = idx(paircontrib(1,pair))
!                 i = nc
!
!#ifdef VAR_LSDEBUG
!                 if ( k == i )call lsquit("THIS IS WRONG AND SHOULD NEVER HAPPEN",-1)
!#endif
!
!                 !extract gvvov(adkc) as dakc and transform d and c, keep a since singles
!                 !are constructed fully -> akcd
!                 call ass_D1to4( w1,    p1, [nv, nv,rpd,nv] )
!                 call ass_D1to4( gvvov, p2, [nv, nv, no, nv] )
!                 do b=1,nv
!                    do ic=1,nv
!                       do a=1,nv
!                          p1(ic,a,1,b) = p2(a,ic, k,b)
!                       enddo
!                    enddo
!                 enddo
!
!                 call dgemm('n','n',nv*nv,pnv,nv,p10,w1,nv*nv,d,nv,nul,w2,nv*nv)
!                 call dgemm('t','n',nv*pnv,pnv,nv,p10,w2,nv,d,nv,nul,w1,nv*pnv)
!
!                 !extract amplitudes as u kcdi with i = nc2
!                 if(k>i)then
!                    call array_reorder_2d(p20,t,pnv,pnv,[2,1],nul,w3)
!                    call array_reorder_2d(m10,t,pnv,pnv,[1,2],p10,w3)
!                 else
!                    call array_reorder_2d(p20,t,pnv,pnv,[1,2],nul,w3)
!                    call array_reorder_2d(m10,t,pnv,pnv,[2,1],p10,w3)
!                 endif
!
!                 call dgemv('n',nv,pnv*pnv,p10,w1,nv,w3, 1, nul,w2,1)
!
!                 !$OMP CRITICAL
!                 o1(:,nc) = o1(:,nc) + w2(1:nv)
!                 !$OMP END CRITICAL
!              enddo
!           else
!
!              !extract gvvov(adkc) as dakc and transform d and c, keep a since singles
!              !are constructed fully -> akcd
!              call ass_D1to4( w1,    p1, [nv, nv,rpd,nv] )
!              call ass_D1to4( gvvov, p2, [nv, nv, no, nv] )
!              do b=1,nv
!                 do jc=1,pno
!                    do ic=1,nv
!                       do a=1,nv
!                          p1(ic,a,jc,b) = p2(a,ic, idx(jc),b)
!                       enddo
!                    enddo
!                 enddo
!              enddo
!              p1 => null()
!              p2 => null()
!
!              call dgemm('n','n',nv*nv*pno,pnv,nv,p10,w1,nv*nv*pno,d,nv,nul,w2,nv*nv*pno)
!              call dgemm('t','n',nv*pno*pnv,pnv,nv,p10,w2,nv,d,nv,nul,w1,nv*pno*pnv)
!
!              !extract amplitudes as u kcdi with i = nc2
!              call ass_D1to4( w3, p3, [pno,pnv,pnv,1] )
!              call ass_D1to4( t,  p2, [pnv,pno,pnv,pno] )
!              do b=1,pnv
!                 do jc=1,pno
!                    do a=1,pnv
!                       p3(jc,a,b,1) = p20 * p2(a,jc,b,i_idx) - p2(a,i_idx,b,jc)
!                    enddo
!                 enddo
!              enddo
!              p2 => null()
!              p3 => null()
!
!              call dgemv('n',nv,pno*pnv*pnv,p10,w1,nv,w3, 1, nul,w2,1)
!
!              !$OMP CRITICAL
!              o1(:,nc) = o1(:,nc) + w2(1:nv)
!              !$OMP END CRITICAL
!
!           endif

           !!!!!!!!!!!!!!!!!!!!!!!!!
           !!!  C1 Term !!!!!!!!!!!!
           !!!!!!!!!!!!!!!!!!!!!!!!!

           if( PS )then
              bpc = 3
              epc = 2
              do pair = 1, paircontribs

                 add_contrib = .false.

                 k = idx(paircontrib(2,pair))
                 i = nc

                 add_contrib = (k<i.and.i==idx(2)).or.(k>i.and.i==idx(1))  

                 if(add_contrib     .and.pair==1)bpc = 1
                 if(.not.add_contrib.and.pair==1)bpc = 2
                 if(add_contrib     .and.pair==2)epc = 2
                 if(.not.add_contrib.and.pair==2)epc = 1

              enddo

              do pair = bpc, epc
                 k = idx(paircontrib(2,pair))
                 i = nc
                 !extract amplitudes as u aick with i = nc2
                 if(i>k)then
                    call array_reorder_2d(p20,t,pnv,pnv,[2,1],nul,w3)
                    call array_reorder_2d(m10,t,pnv,pnv,[1,2],p10,w3)
                 else
                    call array_reorder_2d(p20,t,pnv,pnv,[1,2],nul,w3)
                    call array_reorder_2d(m10,t,pnv,pnv,[2,1],p10,w3)
                 endif

                 !get fock matrix in the current space and transform virtual idx
#ifdef VAR_PTR_RESHAPE
                 r1(1:rpd,1:nv) => w1
                 r2(1:no,1:nv)  => ovf
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                 call c_f_pointer( c_loc(w1(1)), r1, [rpd,nv] )
                 call c_f_pointer( c_loc(ovf(1)), r2,[ no,nv] )
#else
                 call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
                 r1(1,:) = r2(k,:)

                 call dgemm('t','t',pnv,rpd,nv,p10,d,nv,w1,rpd,nul,w2,pnv)

                 !contract amplitudes with fock matrix and transform back
                 call dgemv('n',pnv,rpd*pnv,p10,w3,pnv,w2,1,nul,w1,1)

                 !transform back
                 !$OMP CRITICAL
                 call dgemv('n',nv,pnv,p10,d,nv,w1,1,p10,o1(1,nc),1)
                 !$OMP END CRITICAL
              enddo
           else

              !extract amplitudes as u aick with i = nc2
#ifdef VAR_PTR_RESHAPE
              p3(1:pnv,1:1,1:pnv,1:rpd)    => w3
              p2(1:pnv,1:rpd,1:pnv,1:rpd)  => t
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(w3(1)), p3, [pnv,1,pnv,rpd] )
              call c_f_pointer( c_loc(t(1)),  p2, [pnv,rpd,pnv,rpd] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
              do b=1,pnv
                 do jc=1,pno
                    do a=1,pnv
                       p3(a,1,b,jc) = p20 * p2(a,i_idx,b,jc) - p2(a,jc,b,i_idx)
                    enddo
                 enddo
              enddo
              p2 => null()
              p3 => null()

              !get fock matrix in the current space and transform virtual idx
#ifdef VAR_PTR_RESHAPE
              r1(1:rpd,1:nv) => w1
              r2( 1:no,1:nv) => ovf
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(w1(1)), r1, [rpd,nv] )
              call c_f_pointer( c_loc(ovf(1)), r2,[ no,nv] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
              do b=1,nv
                 do jc=1,pno
                    r1(jc,b) = r2(idx(jc),b)
                 enddo
              enddo
              r1 => null()
              r2 => null()
              call dgemm('t','t',pnv,rpd,nv,p10,d,nv,w1,rpd,nul,w2,pnv)

              !contract amplitudes with fock matrix and transform back
              call dgemv('n',pnv,rpd*pnv,p10,w3,pnv,w2,1,nul,w1,1)

              !transform back
              !$OMP CRITICAL
              call dgemv('n',nv,pnv,p10,d,nv,w1,1,p10,o1(1,nc),1)
              !$OMP END CRITICAL
           endif


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

     call time_start_phase(PHASE_WORK, ttot = mo_time(1) )
     mo_time(6) = mo_time(6) + mo_time(1)

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

     if(DECinfo%PL>2)then
        write (*,'(" PNO: TIME MO part:")')
        write (*,'(" PNO: common contribs      :",g10.3,"s")')mo_time(3)
        write (*,'(" PNO: overlapping contribs :",g10.3,"s")')mo_time(4)
        write (*,'(" PNO: B1                   :",g10.3,"s")')mo_time(5)
        write (*,'(" PNO: C1                   :",g10.3,"s")')mo_time(6)

        call time_start_phase(PHASE_WORK, twall = tfin, ttot = trest, labelttot = &
           & 'PNO: MO part                         :' )
     endif

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
           call tensor_free( pno_t2(ns) )
           call tensor_free( pno_o2(ns) )
           call tensor_free( sio4(ns)   )
        endif

     enddo

     call mem_dealloc( pno_t2 )
     call mem_dealloc( pno_o2 )
     call mem_dealloc( sio4  )
     call mem_dealloc( goovv_vvoo )
     call mem_dealloc( Lvoov )
     call mem_dealloc( gooov )
     call mem_dealloc( p_idx )
     call mem_dealloc( p_nidx )
     call mem_dealloc( s_nidx )
     call mem_dealloc( s_idx )
     call mem_dealloc( vof )
     call mem_dealloc( ovf )
     call mem_dealloc( oof )
     call mem_dealloc( vvf )
     call mem_dealloc( Gkj )

#ifdef VAR_OMP
     call omp_set_nested(nested)
#endif
     if(DECinfo%PL>2)then
        call time_start_phase(PHASE_WORK, ttot = tfin, labelttot = &
           & 'PNO: finalization                    :' )
     endif
#endif
  end subroutine get_ccsd_residual_pno_style



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
#ifdef MOD_UNRELEASED
     integer :: i, j, c, t1,t2,dg, n1
     integer :: ns1,ns2,INFO,lwork,mindim,maxdim,red1,red2,kerdim,diag,remove
     real(realk),pointer:: s1(:,:), s2(:,:), sv(:),U(:,:), VT(:,:),work(:)
     real(realk) :: norm,thr
     real(realk),parameter :: p10 = 1.0E0_realk
     real(realk),parameter :: nul = 0.0E0_realk
     real(realk) :: time_overlap_spaces,mem_overlap_spaces,FracOfMem
     logical :: keep_pair, just_check
     integer :: allremoved, ofmindim, ofmaxdim, allocpcount

     call time_start_phase(PHASE_WORK, twall = time_overlap_spaces )

     call get_currently_available_memory(FracOfMem)
     FracOfMem = 0.1E0_realk * FracOfMem

     mem_overlap_spaces = 0.0E0_realk
     just_check         = .false.

     allocpcount= 0
     allremoved = 0
     ofmindim   = 0
     ofmaxdim   = 0
     call mem_TurnONThread_Memory()
     !$OMP PARALLEL DEFAULT(NONE) &
     !$OMP REDUCTION(+:allremoved,ofmindim,ofmaxdim,allocpcount)&
     !$OMP SHARED(pno_cv,pno_S,n,no,nv,with_svd,thr,mem_overlap_spaces,FracOfMem,just_check)&
     !$OMP PRIVATE(ns1,ns2,i,j,c,s1,s2,norm,sv,U,VT,work,remove,&
     !$OMP lwork,info,diag,kerdim,red1,red2,maxdim,mindim,dg,&
     !$OMP keep_pair)
     call init_threadmemvar()

     !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
     do i=1,n
        do j=1,n

           if(j>=i) cycle

           ! COUNT UPPER TRIANGULAR ELEMENTS WITHOUT DIAGONAL ELEMENTS
           c = (j - i + 1) + i*(i-1)/2

           !$OMP CRITICAL
           just_check = (mem_overlap_spaces/(1024.0E0_realk**3) > FracOfMem)
           !$OMP END CRITICAL

           call get_overlap_matrix_from_pno_spaces(nv,pno_cv(i),pno_cv(j),pno_S(c),with_svd,just_check)

           if(pno_S(c)%allocd)then
              if( with_svd )then
                 !$OMP CRITICAL
                 mem_overlap_spaces = mem_overlap_spaces + &
                    &(size(pno_S(c)%s1) + size(pno_S(c)%d) + size(pno_S(c)%s2) ) * 8.0E0_realk
                 !$OMP END CRITICAL
              else
                 !$OMP CRITICAL
                 mem_overlap_spaces = mem_overlap_spaces + &
                    &( size(pno_S(c)%d) ) * 8.0E0_realk
                 !$OMP END CRITICAL
              endif
           endif
        enddo
     enddo
     !$OMP END DO NOWAIT
     call collect_thread_memory()
     !$OMP END PARALLEL
     call mem_TurnOffThread_Memory()

     if(DECinfo%PL>2)then
        call time_start_phase(PHASE_WORK, ttot = time_overlap_spaces, labelttot = " PNO:&
           & getting overlap spaces:")
        write (*,'("memory requirements for pair overlap info:",g10.3," GB")')mem_overlap_spaces/(1024.0E0_realk**3)
     endif


     if(DECinfo%pno_S_on_the_fly.or.just_check)then

        DECinfo%pno_S_on_the_fly = .true.
        print *,"Switching to calculating pno overlap information on the fly"
        do i=1,n
           do j=1,n
              if(j>=i) cycle
              c = (j - i + 1) + i*(i-1)/2
              if( pno_S(c)%allocd )  call free_PNOSpaceInfo( pno_S(c) )
           enddo
        enddo

     endif

#endif
  end subroutine get_pno_overlap_matrices



  subroutine get_pno_trafo_matrices(no,nv,nb,t_mp2,cv,n,f)
     implicit none
     !ARGUMENTS
     integer, intent(in) :: no, nv, nb, n
     real(realk), intent(in) :: t_mp2(nv,no,nv,no)
     type(PNOSpaceInfo),pointer :: cv(:)
     type(decfrag),intent(in),optional :: f
     !INTERNAL
#ifdef MOD_UNRELEASED
     real(realk) :: virteival(nv),U(nv,nv),PD(nv,nv)
     integer :: i,j,oi,oj,counter, calc_parameters,det_parameters
     integer :: find_pos(no,no),maxocc,maxminocc
     logical :: doit
     logical :: fj
     real(realk), pointer :: w1(:),p1(:,:,:,:),r1(:,:)
     integer :: tid
     real(realk) :: time_pno_spaces, mem_pno_spaces
#ifdef VAR_OMP
     integer, external :: omp_get_num_threads, omp_get_thread_num
     integer :: nt_s, nt_n
#endif

     call time_start_phase(PHASE_WORK, twall = time_pno_spaces )
     mem_pno_spaces = 0.0E0_realk

     find_pos = -1
     fj = present(f)

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


     if( fj )then
        maxocc = f%noccEOS
        if(DECinfo%PNOtriangular)then
           maxminocc = 1
        else
           maxminocc = 2
        endif
     else
        if(DECinfo%PNOtriangular)then
           maxocc    = 1
           maxminocc = 1
        else
           maxocc    = 2
           maxminocc = 2
        endif
     endif



     !call mem_TurnONThread_Memory()
     !OMP PARALLEL DEFAULT(NONE) REDUCTION(+:calc_parameters,det_parameters,&
     !OMP mem_pno_spaces)&
     !OMP SHARED(no,nv,nb,n,fj,f,DECinfo,cv,find_pos,t_mp2,&
     !OMP maxocc,maxminocc)&
     !OMP PRIVATE(counter,virteival,U,PD,doit,tid,w1)
     !call init_threadmemvar()

     tid = 0
#ifdef VAR_OMP
     tid = omp_get_thread_num()
#endif

     !if( tid == 0 )then
     !   call mem_alloc(w1,2*nv*maxocc*nv*maxocc)
     !else
     !   call mem_alloc(w1,2*nv*maxminocc*nv*maxminocc)
     !endif


     if(fj)then

        if(.not.associated(f%VirtMat))then
           call lsquit("Error(get_pno_trafo_matrices)Fragment Correlation density matrix not allocated",-1)
        endif


        !$OMP MASTER
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

        mem_pno_spaces  = mem_pno_spaces  + (cv(1)%ns2**2)*8.0E0_realk

        if(.not.cv(1)%allocd)then
           call lsquit("ERROR(get_pno_trafo_matrices):EOS does not contribute&
              & according to the current threshold, skipping this fragment should be&
              & implemented",-1)
        endif
        !$OMP END MASTER

        !OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
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
                    mem_pno_spaces  = mem_pno_spaces  + (cv(counter)%ns2**2)*8.0E0_realk

                 endif
              endif
           enddo doj
        enddo doi
        !OMP END DO NOWAIT
     else
        !OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
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
                 mem_pno_spaces  = mem_pno_spaces  + (cv(counter)%ns2**2)*8.0E0_realk

              endif

           enddo dojful
        enddo doiful
        !OMP END DO NOWAIT
     endif

     !call mem_dealloc(w1)

     !call collect_thread_memory()
     !OMP END PARALLEL
     !call mem_TurnOffThread_Memory()

     print *,no**2*nv**2,"amplitudes to determine using ",calc_parameters

     call time_start_phase(PHASE_WORK, ttot = time_pno_spaces, labelttot = " PNO:&
     & getting pno spaces:")

     write (*,'("memory requirements for pno space info:",g10.3," GB")')mem_pno_spaces/(1024.0E0_realk**3)
#endif
  end subroutine get_pno_trafo_matrices



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


#ifdef MOD_UNRELEASED
  subroutine get_overlap_idx(n1,n2,cv,idx,nidx,ndidx1,ndidx2)
     implicit none
     integer,intent(in) :: n1,n2
     type(PNOSpaceInfo), intent(in) :: cv(:)
     integer, intent(out) :: idx(:,:),nidx 
     integer, intent(out), optional :: ndidx1,ndidx2
     integer :: n,nc1,nc2,pos,ncidx1,ncidx2
     logical :: yeah_this_is_one

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

           ndidx1 = cv(n1)%n - nidx
           ncidx1 = 0

           do nc1 = 1,cv(n1)%n

              yeah_this_is_one = .true.

              do n=1, nidx
                 if( nc1 == idx(n,1))then
                    yeah_this_is_one = .false.
                 endif
              enddo

              if( yeah_this_is_one )then
                 pos    = pos    + 1
                 ncidx1 = ncidx1 + 1
                 ! pos in first that equals second 
                 idx(pos,1) = nc1
                 ! pos in second that equals first - set invalid since there is none
                 idx(pos,2) = -1
                 ! the aos idx they refer to
                 idx(pos,3) = cv(n1)%iaos(nc1)
              endif
           enddo

           if( ncidx1 /= ndidx1 )then
              print *,ncidx1,ndidx1
              call lsquit("ERROR: stuff dont work",-1)
           endif

        endif

        if(present(ndidx2))then

           ndidx2 = cv(n2)%n - nidx
           ncidx2 = 0

           do nc2 = 1,cv(n2)%n

              yeah_this_is_one = .true.

              do n=1, nidx
                 if( nc2 == idx(n,2))then
                    yeah_this_is_one = .false.
                 endif
              enddo

              if( yeah_this_is_one )then
                 pos    = pos    + 1
                 ncidx2 = ncidx2 + 1
                 ! pos in first that equals second - set invalid since there is none
                 idx(pos,1) = -1
                 ! pos in second that equals first
                 idx(pos,2) = nc2
                 ! the aos idx they refer to
                 idx(pos,3) = cv(n2)%iaos(nc2)
              endif
           enddo
 
           if( ncidx2 /= ndidx2 )then
              print *,ncidx2,ndidx2
              call lsquit("ERROR: stuff dont work",-1)
           endif

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
  subroutine get_overlap_ptr(Sidx,n1,n2,pS,tr1,tr2,st,S,ldS,U,ldU,VT,ldVT,red1,red2,ns1,ns2)
     implicit none
     integer,intent(in) :: Sidx,n1,n2
     type(PNOSpaceInfo), intent(in) :: pS(:)
     character, intent(out) :: tr1,tr2
     logical, intent(out) :: st
     real(realk), pointer, intent(out) :: S(:,:)
     integer, intent(out) :: ldS
     real(realk), pointer, intent(out) :: U(:,:),VT(:,:)
     integer, intent(out) :: ldU, ldVT
     integer, intent(out),optional :: red1, red2, ns1,ns2

     !Get the overlap identificaton and transformation props

     if(n1>n2)then

        !trafo from n1 to n2 
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

        if(.not.pS(Sidx)%allocd)call lsquit("ERROR(get_overlap_ptr):overlap not allocd",-1)

     else if(n2>n1)then

        !trafo from n2 to n1
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

        if(.not.pS(Sidx)%allocd)call lsquit("ERROR(get_overlap_ptr):overlap not allocd",-1)
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
                    cntr             = cntr + 1
                    s_nidx(n1)       = cntr
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
     type(tensor), intent(in) :: pno_o2(n)
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
#ifdef VAR_PTR_RESHAPE
           w1(1:nv,1:rpd,1:nv,1:rpd) => tmp1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer(c_loc(tmp1(1)),w1,[nv,rpd,nv,rpd])
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

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
        cyc       =  .not. S(Sidx)%contributes

     else if(n2>n1)then

        !trafo from n2 to n1
        Sidx      =  (n1 - n2 + 1) + n2 * (n2 - 1 )/2
        cyc       =  .not. S(Sidx)%contributes

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
  subroutine do_overlap_trafo(nv,ns1,ns2,pos_of_overlap,CV,S,m,n,k,A,C,ptr,ptr2,pC)
     implicit none
     integer, intent(in) :: nv,pos_of_overlap,m,n,k,ns1,ns2
     type(PNOSpaceInfo),intent(inout) :: S(:),CV(:)
     real(realk),intent(in),target :: A(:)
     real(realk),intent(inout),target :: C(:)
     real(realk),pointer,optional :: ptr(:),ptr2(:)
     real(realk),optional :: pC
     real(realk),pointer :: S1(:,:)
     real(realk),pointer :: tmp1(:),tmp2(:), U(:,:), VT(:,:)
     integer :: ldS1,ldU,ldVT,Sidx,maxns,minns,cc
     logical :: skiptrafo
     character :: tr1,tr2
     real(realk), parameter :: nul = 0.0E0_realk
     real(realk), parameter :: p10 = 1.0E0_realk
     real(realk) :: preC
     integer :: tid
#ifdef VAR_OMP
     integer, external :: omp_get_thread_num
#endif


     preC = 0.0E0_realk
     if(present(pC))preC = pC

     maxns = max(ns1,ns2)
     minns = min(ns1,ns2)
     Sidx = (minns - maxns + 1) + maxns*(maxns-1)/2

     if(DECinfo%pno_S_on_the_fly.and.ns1/=ns2)then
#ifdef VAR_OMP
        tid  = omp_get_thread_num()
        Sidx = tid + 1
#endif
        call get_overlap_matrix_from_pno_spaces(nv,CV(maxns),CV(minns),S(Sidx),.true.,.false.)
     endif

     call get_overlap_ptr(Sidx,ns1,ns2,S,tr1,tr2,skiptrafo,S1,ldS1,U=U,ldU=ldU,VT=VT,ldVT=ldVT)


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

     if(DECinfo%pno_S_on_the_fly.and.ns1/=ns2)then
        if( S(Sidx)%allocd )  call free_PNOSpaceInfo( S(Sidx) )
     endif
  end subroutine do_overlap_trafo


  
  subroutine get_overlap_matrix_from_pno_spaces(nv,pno1,pno2,S12,with_svd,just_check)
     implicit none
     integer, intent(in) :: nv
     type(PNOSpaceInfo),intent(in) :: pno1,pno2
     type(PNOSpaceInfo),intent(inout) :: S12
     logical, intent(in) :: with_svd, just_check
     integer :: t1,t2,dg, n1, maxindex, minindex
     integer :: ns1,ns2,INFO,lwork,mindim,maxdim,red1,red2,kerdim,diag,remove
     real(realk),pointer:: s1(:,:), s2(:,:), sv(:),U(:,:), VT(:,:),work(:)
     real(realk) :: norm,thr
     real(realk),parameter :: p10 = 1.0E0_realk
     real(realk),parameter :: nul = 0.0E0_realk
     real(realk) :: time_overlap_spaces,mem_overlap_spaces
     logical :: keep_pair

     !FIXME: REMOVE ALLOCATION INSIDE THIS SUBROUTINE

     if( DECinfo%noPNOoverlaptrunc ) then
        thr = -1.0E0_realk * huge(thr)
     else
        thr = DECinfo%PNOoverlapthr
     endif

     ns1 = pno1%ns2
     ns2 = pno2%ns2

     S12%ns1 = ns1
     S12%ns2 = ns2

     keep_pair = ( pno1%allocd .and. pno2%allocd )

     if( keep_pair )then


        S12%s_associated = .false.
        call mem_alloc(S12%d,ns1,ns2)

        s1 => pno1%d
        s2 => pno2%d

        call dgemm('t','n',ns1,ns2,nv,p10,s1,nv,s2,nv,nul,S12%d,ns1)

        if(with_svd)then

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
           call dgesvd('A','A',ns1,ns2,S12%d,ns1,sv,U,ns1,VT,ns2,work,lwork,INFO)
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

           keep_pair = ( red1 > 0 .and. red2 > 0 )

           if(red1/=diag.or.red2/=diag)call &
              &lsquit("ERROR(get_pno_overlap_matrices)calculated wrong dimensions",-1)


           call mem_dealloc( S12%d )

           if( keep_pair )then

              S12%s_associated = .true.

              call mem_alloc( S12%s1, ns1,  red1 )
              call mem_alloc( S12%s2, red2, ns2  )
              call mem_alloc( S12%d,  red1, red2 )

              S12%red1 = red1
              S12%red2 = red2

              S12%s1 = U(:,1:red1)
              S12%s2 = VT(1:red2,:)
              S12%d  = nul

              do dg = 1, diag
                 S12%d(dg,dg) = sv(dg)
              enddo


           endif


           call mem_dealloc( U )
           call mem_dealloc( VT )
           call mem_dealloc( work )
           call mem_dealloc( sv )
        endif


        if( keep_pair ) then

           S12%n = 2
           call mem_alloc(S12%iaos,S12%n)
           maxindex = max(pno1%iaos(1),pno2%iaos(1))
           minindex = min(pno1%iaos(1),pno2%iaos(1))
           S12%iaos = [maxindex,minindex]

        endif

        S12%allocd      = keep_pair
        S12%contributes = keep_pair

        if(just_check)then

           call mem_dealloc( S12%iaos )
           call mem_dealloc( S12%s1   )
           call mem_dealloc( S12%s2   )
           call mem_dealloc( S12%d    )

           S12%allocd = .false.

        endif

     else

        S12%allocd      = .false.
        S12%contributes = .false.

     endif
  end subroutine get_overlap_matrix_from_pno_spaces

  subroutine init_pno_residual_and_sio4(cv,o2,sio4,n,no)
     implicit none
     integer, intent(in) :: n,no
     type(PNOSpaceInfo),intent(in) :: cv(n)
     type(tensor),intent(inout) :: o2(n),sio4(n)
     integer :: nn, pnv,pno,rpd

     do nn=1,n

        pnv = cv(nn)%ns2
        pno = cv(nn)%n
        rpd = cv(nn)%rpd

        if(cv(nn)%allocd)then

           call tensor_init(o2(nn), [pnv,rpd,pnv,rpd],4)
           call tensor_zero(o2(nn))


           if( cv(nn)%PS )then
              call tensor_init(sio4(nn),[no,no],2)
           else
              call tensor_init(sio4(nn),[no,no,pno,pno],4)
           endif

           call tensor_zero(sio4(nn))

        endif
     enddo
  end subroutine init_pno_residual_and_sio4

  subroutine get_pno_amplitudes(t2,cv,pno_t2,n,no,nv)
     implicit none
     integer, intent(in) :: n,no,nv
     real(realk),intent(in) :: t2(:,:,:,:)
     type(PNOSpaceInfo), intent(in) :: cv(n)
     type(tensor), intent(inout) :: pno_t2(n)
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

           call tensor_init(pno_t2(nn), [pnv,rpd,pnv,rpd],4)

#ifdef VAR_PTR_RESHAPE
           w1(1:nv,1:rpd,1:nv,1:rpd) => tmp1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer(c_loc(tmp1(1)),w1,[nv,rpd,nv,rpd])
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

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

     if(DECinfo%PL>3)then
        call print_norm(PD,i8*nv*nv,"PD matrix norm:")
     endif
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


  subroutine get_batch_sizes_pno_residual(mylsitem,no,nv,nb,maxocc,maxvirt,nspaces,a_batch,&
        &g_batch,sio4,pno_cv,pno_t2,pno_o2,query)
     implicit none
     type(lsitem), intent(inout) :: mylsitem
     integer, intent(in) :: no,nv,nb,maxocc,maxvirt,nspaces
     type(int_batch),intent(inout) :: a_batch, g_batch
     type(tensor), intent(inout) :: sio4(nspaces), pno_t2(nspaces),pno_o2(nspaces)
     type(PNOSpaceInfo),intent(inout) :: pno_cv(nspaces)
     type(pno_query_info),intent(inout) :: query
     logical :: converged
     integer :: MaxAllowedDimGamma, MaxAllowedDimAlpha,newsize,oldsize,MinAObatch
     integer :: gamma_size,alpha_size, basic,narray,ns
     real(realk) :: MemFree,used_mem, Frac_of_mem
     real(realk), pointer :: dummy1(:), dummy2(:,:)
     integer(kind=8) :: nelms
     integer :: edit,max_nthreads,max_w_per_thr(5,nloops),loop,thread
#ifdef VAR_OMP
     integer, external :: omp_get_num_threads, omp_get_max_threads
#endif

     !this should coincide with the basic number of arrays (offset_for_omp) in pno_residual_integral_direct_loop
     basic = 11
     call init_query_info(query,basic+(nspaces*5*nloops))

     call get_currently_available_memory(MemFree)

     call mem_alloc(dummy1,1)
     call mem_alloc(dummy2,1,1)

#ifdef VAR_MPI
     !call lsmpi_reduce_realk_min(MemFree,infpar%master,infpar%lg_comm)
#endif

     nthreads_level1_int_dir = 1
     max_nthr_int_loop       = 1
#ifdef VAR_OMP
     max_nthr_int_loop       = omp_get_max_threads()
#endif

     frac_of_mem = 0.8_realk*MemFree

     call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch,'R')
        
     converged  = .false.
     gamma_size = nb

     do while(.not. converged)

        call init_batch_info(mylsitem,g_batch,gamma_size,nb)

        do alpha_size = gamma_size,MinAObatch,-1

           call init_batch_info(mylsitem,a_batch,alpha_size,nb)

           query%size_array = 0

           call pno_residual_integral_direct_loop(mylsitem,dummy1,1_long,dummy1,1_long,dummy1,1_long,dummy1,1_long,&
              &dummy1,1_long,no,nv,nb,maxocc,maxvirt,nspaces,a_batch,g_batch,sio4,pno_cv,pno_t2,pno_o2,dummy2,dummy2,dummy2,dummy2,&
              & dummy1,dummy1,dummy1,dummy1,dummy2, query = query )

           !analyze query information
           nelms = 0
           do narray = 1, basic
              nelms = nelms + query%size_array(narray)
           enddo

           used_mem = (nelms*8.0E0_realk)/(1.024E3_realk**3)

           converged = ( used_mem < frac_of_mem )

           if(converged)then
              exit
           else
              call free_batch_info(a_batch)
           endif

        enddo

        if(.not.converged)then

           call free_batch_info(g_batch)

           oldsize    = gamma_size
           gamma_size = oldsize - 1

           if( gamma_size < MinAObatch )then
              call lsquit("ERROR(get_batch_sizes_pno_residual): Not enough&
                 & memory for a pno calculation",-1)
           endif

        endif

     enddo

     call mem_dealloc(dummy1)
     call mem_dealloc(dummy2)

     !do narray = 1, query%n_arrays
     !   print*,"Array",narray," of size:",query%size_array(narray)
     !enddo

     max_w_per_thr = 0
     do loop = 1, nloops
        do narray = 1,5
           do ns = 1, nspaces
              edit = basic - 1 + (ns-1) * nloops * 5 + (loop - 1)  * 5
              max_w_per_thr(narray,loop) = max(max_w_per_thr(narray,loop),query%size_array(edit+narray))
           enddo
        end do
        !print *,"query says",max(query%size_array(:)),"max is",max_w_per_thr(:,loop)
     enddo


#ifdef VAR_OMP

     nthreads_level1_int_dir(first_loop)  = omp_get_max_threads()
     nthreads_level1_int_dir(second_loop) = omp_get_max_threads()
     nthreads_level1_int_dir(third_loop)  = omp_get_max_threads()


     do narray = 1, 5

        if(max_w_per_thr(narray,first_loop)/=0.and.(narray/=1))then
           nthreads_level1_int_dir(first_loop)  =&
              & min(nthreads_level1_int_dir(first_loop), query%size_array(narray)/max_w_per_thr(narray,first_loop))
        endif
        if(max_w_per_thr(narray,second_loop)/=0.and.(narray/=3))then
           nthreads_level1_int_dir(second_loop) =&                                                      
              & min(nthreads_level1_int_dir(second_loop),query%size_array(narray)/max_w_per_thr(narray,second_loop))
        endif
        if(max_w_per_thr(narray,third_loop)/=0.and.(narray/=4.and.narray/=5))then
           nthreads_level1_int_dir(third_loop)  =&                                                      
              & min(nthreads_level1_int_dir(third_loop), query%size_array(narray)/max_w_per_thr(narray,third_loop))
        endif

     enddo
#endif

     max_nthr_int_loop = max(nthreads_level1_int_dir(first_loop),&
        & nthreads_level1_int_dir(second_loop),&
        & nthreads_level1_int_dir(third_loop))


     if(DECinfo%PL>3)then
        print *,max_nthr_int_loop,"<-max,nthreads->",nthreads_level1_int_dir
     endif

     call mem_alloc(info_omp1,2,max_nthr_int_loop,5,nloops)

     !initialize with invalid numbers
     info_omp1 = -1

     do loop = 1, nloops
        do narray = 1, 5
           do thread=1, max_nthr_int_loop
              if(thread<=nthreads_level1_int_dir(loop))then
                 info_omp1(1,thread,narray,loop) = 1 + (thread-1)*max_w_per_thr(narray,loop)
                 info_omp1(2,thread,narray,loop) = ( thread )*max_w_per_thr(narray,loop)
                 !print *,"INFO FOR thread",thread," and array ",narray," in loop ",loop
                 !print *,info_omp1(1,thread,narray,loop), info_omp1(2,thread,narray,loop)
              endif
           enddo
        enddo
     enddo


     if(DECinfo%PL>2)write(*,'("INFO: allocating ",g10.3," GB in PNO integral &
        &direct loop with #a:",I6," and #g:",I6)')used_mem,a_batch%nbatches,g_batch%nbatches

  end subroutine get_batch_sizes_pno_residual


  subroutine pno_residual_integral_direct_loop(mylsitem,w1,s1,w2,s2,w3,s3,w4,s4,w5,s5,no,nv,nb,maxocc,maxvirt,&
         & nspaces,a_batch,g_batch,sio4,pno_cv,pno_t2,pno_o2,xo,xv,yo,yv,&
         & gooov,goovv,govov,Lvoov,Gai, query )

     implicit none
     !Arguments
     type(lsitem), intent(inout) :: mylsitem
     integer(kind=8),intent(in)  :: s1,   s2,   s3,    s4,   s5
     real(realk), target, intent(inout)  :: w1(:),w2(:),w3(:), w4(:),w5(:)
     integer, intent(in) :: no, nv, nb, maxocc,maxvirt,nspaces
     type(int_batch),intent(in) :: a_batch, g_batch
     type(tensor), intent(inout) :: sio4(nspaces), pno_t2(nspaces),pno_o2(nspaces)
     type(PNOSpaceInfo),intent(inout) :: pno_cv(nspaces)
     real(realk), intent(in), target :: xo(nb,no)
     real(realk), intent(in) :: xv(:,:), yo(:,:), yv(:,:)
     real(realk), intent(inout) :: gooov(:), goovv(:),govov(:),Lvoov(:),Gai(:,:)
     type(pno_query_info),optional, intent(inout) :: query
     !internal variables
     integer(kind=8) :: myload,var_inp(4)
     integer :: ns
     integer, pointer :: idx(:)
     real(realk), pointer :: d(:,:),t(:),o(:),h1(:),h2(:),h3(:),r1(:,:)
#ifdef VAR_PTR_RESHAPE
     real(realk), contiguous, pointer :: xv_pair(:,:,:)
#else
     real(realk), pointer :: xv_pair(:,:,:)
#endif
     real(realk), pointer :: xo_pair(:,:,:),yo_pair(:,:,:),yv_pair(:,:,:)
     real(realk), pointer :: tpl(:,:), tmi(:,:)
     integer :: max_pnor, max_pnvr
     integer :: pno_comb
     integer :: pair,paircontribs,paircontrib(2,2),rpd,pno,pnv
     integer :: goffs,aoffs,tlen,tred,nor,nvr,pos1,pos2,beg1,beg2,ic,i
     integer :: alphaB,gammaB,dimAlpha,dimGamma,fa,fg,la,lg,xa,xg,a
     integer :: offset_for_omp,contract
     logical :: FoundInMem,fullRHS, doscreen, PS, EOS, add_contrib
     type(DECscreenITEM)    :: DecScreen
     character :: INTSPEC(5)
     logical :: master, I_PLUS_MINUS_DONE,SORT_INT_TO_W2_DONE, this_is_query, this_is_not_query
     real(realk) :: tw,tc,tbeg,times(9),times_in_loops(8,3)
     integer(kind=8) :: my_s1,my_s2,my_s3,my_s4,my_s5
     real(realk), pointer :: my_w1(:),my_w2(:),my_w3(:),my_w4(:),my_w5(:)
     real(realk), target  :: dummy(1)
     real(realk), parameter :: p20 = 2.0E0_realk
     real(realk), parameter :: p10 = 1.0E0_realk
     real(realk), parameter :: m10 = -1.0E0_realk
     real(realk), parameter :: m05 = -0.5E0_realk
     real(realk), parameter :: p05 = 0.5E0_realk
     real(realk), parameter :: nul = 0.0E0_realk
     integer :: tid, edit,orig_nthr
#ifdef VAR_OMP
     logical :: nested
     logical, external :: omp_get_nested
     integer, external :: omp_get_thread_num, omp_get_num_threads
     nested    = omp_get_nested()
     call omp_set_nested(.true.)
#endif

     master        = .true.
     tw            = 0.0E0_realk
     tc            = 0.0E0_realk
     times         = 0.0E0_realk
     times_in_loops= 0.0E0_realk
     this_is_query = present(query)

     this_is_not_query = .not. this_is_query

     paircontribs  = 2
     paircontrib(1:2,1) = [1,2]
     paircontrib(1:2,2) = [2,1]

     max_pnor = maxocc  * (maxocc  + 1)/2
     max_pnvr = maxvirt * (maxvirt + 1)/2

     if(this_is_query)then

        !make sure this offset is the last accessed element in this if statement
        offset_for_omp = 10 

        if(query%n_arrays /= offset_for_omp + 1 + (nspaces*5*nloops) )then
           call lsquit("ERROR(pno_residual_integral_direct_loop):query not&
              & possible, n_arrays needs to be set correctly, please check&
              & offset_for_omp",-1)
        endif

        query%size_array(6)  = i8*nb*maxocc*max_nthr_int_loop
        query%size_array(7)  = i8*nb*maxocc*max_nthr_int_loop
        query%size_array(8)  = i8*nb*maxvirt*max_nthr_int_loop
        query%size_array(9)  = i8*nb*maxvirt*max_nthr_int_loop
        query%size_array(10) = i8*max_pnvr*max_pnor*max_nthr_int_loop
        query%size_array(11) = i8*max_pnvr*max_pnor*max_nthr_int_loop

        var_inp = 0
     else
        call mem_alloc( xo_pair, nb, maxocc, max_nthr_int_loop    )
        call mem_alloc( yo_pair, nb, maxocc, max_nthr_int_loop    )
        call mem_alloc( xv_pair, nb, maxvirt, max_nthr_int_loop   )
        call mem_alloc( yv_pair, nb, maxvirt, max_nthr_int_loop   )

        call mem_alloc( tpl, max_pnor*max_pnvr, max_nthr_int_loop )
        call mem_alloc( tmi, max_pnor*max_pnvr, max_nthr_int_loop )
     endif



     ! Set integral info
     ! *****************
     INTSPEC(1)               = 'R' !R = Regular Basis set on the 1th center 
     INTSPEC(2)               = 'R' !R = Regular Basis set on the 2th center 
     INTSPEC(3)               = 'R' !R = Regular Basis set on the 3th center 
     INTSPEC(4)               = 'R' !R = Regular Basis set on the 4th center 
     INTSPEC(5)               = 'C' !C = Coulomb operator
     doscreen                 = MyLsItem%setting%scheme%cs_screen.OR.MyLsItem%setting%scheme%ps_screen

     ! ***********************************************
     ! *  precalculate the full screening matrix     *
     ! ***********************************************

     if(.not.this_is_query)then
        ! This subroutine builds the full screening matrix.
        call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting,&
           & a_batch%nbatches,g_batch%nbatches,INTSPEC,DECinfo%IntegralThreshold)
        if(mylsitem%setting%scheme%cs_screen.or.mylsitem%setting%scheme%ps_screen)then
           call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
              & nb,a_batch%nbatches,g_batch%nbatches,&
              & a_batch%batchsize,g_batch%batchsize,a_batch%batchindex,g_batch%batchindex,&
              & a_batch%batchdim,g_batch%batchdim,INTSPEC,DECinfo%output,DECinfo%output)
           call II_getBatchOrbitalScreenK(DecScreen,mylsitem%setting,&
              & nb,a_batch%nbatches,g_batch%nbatches,a_batch%batchsize,g_batch%batchsize,&
              & a_batch%batchindex,g_batch%batchindex,&
              & a_batch%batchdim,g_batch%batchdim,INTSPEC,DECinfo%output,DECinfo%output)
        endif

        call time_start_phase(PHASE_WORK, twall = tbeg )

     endif

     myload = 0

     fullRHS = g_batch%nbatches.EQ.1.AND.a_batch%nbatches.EQ.1


     !INTEGRAL DIRECT STUFF HAPPENINING HERE
     BatchGamma: do gammaB = 1,g_batch%nbatches  ! AO batches
        lg = g_batch%batchdim(gammaB)              ! Dimension of gamma batch
        fg = g_batch%batch2orb(gammaB)%orbindex(1) ! First index in gamma batch
        xg = g_batch%batch2orb(gammaB)%orbindex(lg)


        BatchAlpha: do alphaB = 1, a_batch%nbatches

           !check if the current job is to be done by current node
           !call check_job(scheme,first_round,dynamic_load,alphaB,gammaB,nbatchesAlpha,&
           !  &nbatchesGamma,tasks,tasksw,print_debug)
           !break the loop if alpha become too large, necessary to account for all
           !of the mpi and non mpi schemes, this is accounted for, because static,
           !and dynamic load balancing are enabled
           !if(alphaB>a_batch%nbatches) exit
           !print *,"JOB ",alphaB,a_batch%nbatches,gammaB,g_batch%nbatches
           if(DECinfo%PL>2.and.this_is_not_query)then
#ifdef VAR_MPI
              write (*, '("Rank ",I3," starting job (",I3,"/",I3,",",I3,"/",I3,")")') infpar%mynum,&
                 &alphaB,a_batch%nbatches,gammaB,g_batch%nbatches
#else
              write (*, '("starting job (",I3,"/",I3,",",I3,"/",I3,")")')&
                 &alphaB,a_batch%nbatches,gammaB,g_batch%nbatches
#endif
           endif

           la = a_batch%batchdim(alphaB)                              ! Dimension of alpha batch
           fa = a_batch%batch2orb(alphaB)%orbindex(1)                 ! First index in alpha batch
           xa = a_batch%batch2orb(alphaB)%orbindex(la)

 
           !short hand notation
           myload     = myload + la * lg

           I_PLUS_MINUS_DONE   = .false.
           SORT_INT_TO_W2_DONE = .false.


           if(this_is_query)then

              query%size_array(1) = max(query%size_array(1),(i8*la*nb)*lg*nb)

              var_inp = 0

           else

              call time_start_phase(PHASE_WORK, twall = times(1) )

              IF(doscreen)Mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGabKLHS(alphaB)%p
              IF(doscreen)Mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGabKRHS(gammaB)%p

              call II_GET_DECPACKED4CENTER_K_ERI(DECinfo%output,DECinfo%output, &
                 & Mylsitem%setting,w1,a_batch%batchindex(alphaB),g_batch%batchindex(gammaB),&
                 & a_batch%batchsize(alphaB),g_batch%batchsize(gammaB),la,nb,lg,nb,INTSPEC,&
                 & fullRHS,DECinfo%IntegralThreshold)

              call time_start_phase(PHASE_WORK, ttot = times(1) )
              times(3) = times(3) + times(1)

           endif


           !TODO: IMPLEMENT LOOP OVER occ occ on the outside for simple one
           !sided updates if integrals are stored in batches in occ-occ and
           !reorder such that the first trafo is to occ
           if(this_is_query)then
              query%size_array(w2_tag) = max(query%size_array(w2_tag),(i8*la*nb)*lg*nv)
              query%size_array(w3_tag) = max(query%size_array(w3_tag),(i8*lg*nv)*la*nb)
              query%size_array(w4_tag) = max(query%size_array(w4_tag),(i8*lg*nv)*la*no)
              !goovv              tag
              query%size_array(w2_tag) = max(query%size_array(w2_tag),(i8*nv*nv)*la*no)
              query%size_array(w3_tag) = max(query%size_array(w3_tag),(i8*la*nb)*lg*nv)
              !gooov and Lvoov    tag
              query%size_array(w2_tag) = max(query%size_array(w2_tag),(i8*no*no)*la*no)
              query%size_array(w3_tag) = max(query%size_array(w3_tag),(i8*nv*nv)*la*no)
           else
              call dgemm('n','n',la*nb*lg,nv,nb,p10,w1,la*nb*lg,yv,nb,nul,w2,la*nb*lg)
              call array_reorder_4d(p10,w2,la,nb,lg,nv,[3,4,1,2],nul,w3)
              call dgemm('n','n',lg*nv*la,no,nb,p10,w3,lg*nv*la,yo,nb,nul,w4,lg*nv*la)
              !goovv
              call dgemm('t','n',nv,nv*la*no,lg,p10,xv(fg,1),nb,w4,lg,nul,w2,nv)
              call array_reorder_4d(p10,w2,nv,nv,la,no,[3,4,1,2],nul,w3)
              call dgemm('t','n',no,no*nv*nv,la,p10,xo(fa,1),nb,w3,la,p10,goovv,no)
              !gooov and Lvoov
              call dgemm('t','n',no,nv*la*no,lg,p10,xo(fg,1),nb,w4,lg,nul,w2,no)
              call array_reorder_4d(p10,w2,no,nv,la,no,[3,4,1,2],nul,w3)
              call dgemm('t','n',no,no*no*nv,la,p10,xo(fa,1),nb,w3,la,p10,gooov,no)
              call dgemm('t','n',nv,no*no*nv,la,p10,xv(fa,1),nb,w3,la,p10,Lvoov,nv)
           endif

           call time_start_phase(PHASE_WORK, twall = times(1), at = times(4) )

           !$OMP PARALLEL  DEFAULT(NONE) IF(this_is_not_query)&
           !$OMP SHARED(nspaces,nthreads_level1_int_dir,&
           !$OMP pno_cv,pno_t2,pno_o2,offset_for_omp,no,nv,nb,xo,xv,yo,yv,&
           !$OMP fa,la,fg,lg,w1,w2,w3,w4,w5,query,this_is_query,this_is_not_query,&
           !$OMP sio4,info_omp1,xv_pair,xo_pair,yv_pair,yo_pair) &
           !$OMP PRIVATE(a,d,t,idx,pnv,pno,rpd,PS,o,ns,i,&
           !$OMP pno_comb,beg1,beg2,goffs,aoffs,nor,nvr,tlen,tred,my_w1,my_w2,my_w3,&
           !$OMP my_w4,my_w5,tid,EOS,edit) REDUCTION(+:times,times_in_loops)&
           !$OMP NUM_THREADS(nthreads_level1_int_dir(first_loop))

           tid = 0
#ifdef VAR_OMP
           tid = omp_get_thread_num()
#endif

           if(this_is_not_query)then
              my_w2 => w2(info_omp1(beg_array,tid+1,w2_tag,first_loop):info_omp1(end_array,tid+1,w2_tag,first_loop))
              my_w3 => w3(info_omp1(beg_array,tid+1,w3_tag,first_loop):info_omp1(end_array,tid+1,w3_tag,first_loop))
              my_w4 => w4(info_omp1(beg_array,tid+1,w4_tag,first_loop):info_omp1(end_array,tid+1,w4_tag,first_loop))

              my_w1 => null()
              my_w5 => null()
           endif

           !$OMP DO SCHEDULE(DYNAMIC)
           do ns = 1, nspaces

              if(.not.pno_cv(ns)%allocd)then

                 cycle 

              endif

              d   => pno_cv(ns)%d
              t   => pno_t2(ns)%elm1
              idx => pno_cv(ns)%iaos
              pnv =  pno_cv(ns)%ns2
              pno =  pno_cv(ns)%n
              rpd =  pno_cv(ns)%rpd
              PS  =  pno_cv(ns)%PS
              EOS =  pno_cv(ns)%is_FA_space 

              o   => pno_o2(ns)%elm1

              nor=rpd*(rpd+1)/2


              if(this_is_query)then

                 edit = offset_for_omp+(ns-1)*nloops*5+(first_loop-1)*5

                 query%size_array(edit+w2_tag) = max(query%size_array(edit+w2_tag),(i8*la*nb)*lg*rpd)
                 query%size_array(edit+w3_tag) = max(query%size_array(edit+w3_tag),(i8*la*nb)*lg*rpd)
                 query%size_array(edit+w4_tag) = max(query%size_array(edit+w4_tag),(i8*lg*rpd)*la*rpd)

                 !VOVO
                 query%size_array(edit+w2_tag) = max(query%size_array(edit+w2_tag),(i8*pnv*rpd)*la*rpd)
                 query%size_array(edit+w3_tag) = max(query%size_array(edit+w3_tag),(i8*pnv*rpd)*la*rpd)

                 !SIO4
                 query%size_array(edit+w2_tag) = max(query%size_array(edit+w2_tag),(i8*no*rpd)*la*rpd)
                 query%size_array(edit+w3_tag) = max(query%size_array(edit+w3_tag),(i8*rpd*rpd)*no*no)

              else

                 !Get the transformation matrices
                 !  xo_pair 
                 do i = 1, pno
                    do a = 1, nb
                       xo_pair(a,i,tid+1) = xo(a,idx(i))
                    enddo
                 enddo
                 !  xv_pair
                 call dgemm('n','n',nb,pnv,nv,p10,xv,nb,d,nv,nul,xv_pair(1,1,tid+1),nb)
                 !  yo_pair
                 do i = 1, pno
                    do a = 1, nb
                       yo_pair(a,i,tid+1) = yo(a,idx(i))
                    enddo
                 enddo
                 !  yv_pair
                 call dgemm('n','n',nb,pnv,nv,p10,yv,nb,d,nv,nul,yv_pair(1,1,tid+1),nb)
#ifndef VAR_OMP
                 call time_start_phase(PHASE_WORK, at = times_in_loops(1,1) )
#endif

                 !gvovo contribution
                 if( PS )then
                    beg1 = 2
                    beg2 = 1
                 else
                    beg1 = 1
                    beg2 = 1
                 endif
                 !VOVO and sio4
                 call dgemm('n','n',la*nb*lg,rpd,nb,p10,w1,la*nb*lg,yo_pair(1,beg1,tid+1),nb,nul,my_w2,la*nb*lg)
                 call array_reorder_4d(p10,my_w2,la,nb,lg,rpd,[3,4,1,2],nul,my_w3)
                 call dgemm('n','n',lg*rpd*la,rpd,nb,p10,my_w3,lg*rpd*la,yo_pair(1,beg2,tid+1),nb,nul,my_w4,lg*rpd*la)
                 !VOVO contribution
                 call dgemm('t','n',pnv,rpd*la*rpd,lg,p10,xv_pair(fg,1,tid+1),nb,my_w4,lg,nul,my_w2,pnv)
                 call array_reorder_4d(p10,my_w2,pnv,rpd,la,rpd,[3,4,1,2],nul,my_w3)
                 call dgemm('t','n',pnv,rpd*pnv*rpd,la,p10,xv_pair(fa,1,tid+1),nb,my_w3,la,p10,o,pnv)
                 !sio4
                 call dgemm('t','n',no,rpd*la*rpd,lg,p10,xo(fg,1),nb,my_w4,lg,nul,my_w2,no)
                 call array_reorder_4d(p10,my_w2,no,rpd,la,rpd,[2,4,1,3],nul,my_w3)
                 call dgemm('n','n',rpd*rpd*no,no,la,p10,my_w3,rpd*rpd*no,xo(fa,1),nb,p10,sio4(ns)%elm1,rpd*rpd*no)

#ifndef VAR_OMP
                 call time_start_phase(PHASE_WORK, at = times_in_loops(3,1) )
#endif

              endif

           enddo
           !$OMP END DO NOWAIT
           !$OMP END PARALLEL


           call time_start_phase(PHASE_WORK, ttot = times(1))
           times(5) = times(5) + times(1)

           !shift intrgrals to w2 in the order \delta \gamma \alpha \beta
           if(this_is_query)then

              query%size_array(w3_tag) = max(query%size_array(w3_tag),(i8*la*nb)*lg*nb)

           else if(.not.SORT_INT_TO_W2_DONE)then

              !\alpha \beta \gamma \delta -> \delta \gamma \alpha \beta
              call array_reorder_4d(p10,w1,la,nb,lg,nb,[4,3,1,2],nul,w3)
              SORT_INT_TO_W2_DONE = .true.

           endif

           call time_start_phase(PHASE_WORK, at = times(6), twall = times(2)  )

           !$OMP PARALLEL  DEFAULT(NONE) IF(this_is_not_query)&
           !$OMP SHARED(nspaces,paircontrib,xa,xg,Gai,nthreads_level1_int_dir,&
           !$OMP pno_cv,pno_t2,pno_o2,offset_for_omp,no,nv,nb,xo,xv,yo,yv,&
           !$OMP fa,la,fg,lg,w1,w2,w3,w4,w5,query,this_is_query,this_is_not_query,&
           !$OMP sio4,info_omp1,xv_pair,xo_pair,yv_pair,yo_pair) &
           !$OMP PRIVATE(a,d,t,idx,pnv,pno,rpd,PS,o,ns,i,h1,h2,contract,&
           !$OMP pno_comb,beg1,beg2,goffs,aoffs,nor,nvr,tlen,tred,my_w1,my_w2,my_w3,&
           !$OMP my_w4,my_w5,tid,EOS,edit,r1) REDUCTION(+:times,times_in_loops)&
           !$OMP NUM_THREADS(nthreads_level1_int_dir(second_loop))
           tid = 0
#ifdef VAR_OMP
           tid = omp_get_thread_num()
#endif
           if(this_is_not_query)then
              my_w1 => w1(info_omp1(beg_array,tid+1,w1_tag,second_loop):info_omp1(end_array,tid+1,w1_tag,second_loop))
              my_w2 => w2(info_omp1(beg_array,tid+1,w2_tag,second_loop):info_omp1(end_array,tid+1,w2_tag,second_loop))
              my_w4 => w4(info_omp1(beg_array,tid+1,w4_tag,second_loop):info_omp1(end_array,tid+1,w4_tag,second_loop))
              my_w5 => w5(info_omp1(beg_array,tid+1,w5_tag,second_loop):info_omp1(end_array,tid+1,w5_tag,second_loop))

              my_w3 => null()
           endif
           !$OMP DO SCHEDULE(DYNAMIC)
           do ns = 1, nspaces

              if(.not.pno_cv(ns)%allocd)then

                 cycle 

              endif

              d   => pno_cv(ns)%d
              t   => pno_t2(ns)%elm1
              idx => pno_cv(ns)%iaos
              pnv =  pno_cv(ns)%ns2
              pno =  pno_cv(ns)%n
              rpd =  pno_cv(ns)%rpd
              PS  =  pno_cv(ns)%PS

              o   => pno_o2(ns)%elm1

              if(this_is_query)then

                 edit = offset_for_omp+(ns-1)*nloops*5+(second_loop-1)*5

                 !CAREFUL: THIS ACCOUNTS FOR BOTH CASES PS AND NOT PS, WHEN
                 !CHANGING THE ROUTINE BELOW, THE QUERY HAS TO BE ADAPTED
                 !ACCORDINGLY
                 query%size_array(edit+w2_tag) = max(query%size_array(edit+w2_tag),(i8*rpd*lg)*la*nb)
                 query%size_array(edit+w1_tag) = max(query%size_array(edit+w1_tag),(i8*rpd*lg)*la*pnv)
                 query%size_array(edit+w2_tag) = max(query%size_array(edit+w2_tag),(i8*rpd*lg)*la*pnv)
                 query%size_array(edit+w4_tag) = max(query%size_array(edit+w4_tag),(i8*pnv*pnv)*rpd*rpd)

                 if( pnv < lg )then
                    query%size_array(edit+w1_tag) = max(query%size_array(edit+w1_tag),(i8*la*pnv)*rpd*pnv)
                 else
                    query%size_array(edit+w1_tag) = max(query%size_array(edit+w1_tag),(i8*rpd*pnv)*rpd*lg)
                 endif

                 query%size_array(edit+w5_tag) = max(query%size_array(edit+w5_tag),(i8*la*rpd))

              else
                 !Get the transformation matrices
                 !  xo_pair 
                 do i = 1, pno
                    do a = 1,nb
                       xo_pair(a,i,tid+1) = xo(a,idx(i))
                    enddo
                 enddo

                 !  yv_pair
                 call dgemm('n','n',nb,pnv,nv,p10,yv,nb,d,nv,nul,yv_pair(1,1,tid+1),nb)
#ifndef VAR_OMP
                 call time_start_phase(PHASE_WORK, at = times_in_loops(1,2)  )
#endif

                 !Get G_{\alpha i}

                 if( PS )then
#ifdef VAR_PTR_RESHAPE
                    r1(1:la,1:rpd) => my_w5
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                    call c_f_pointer(c_loc(my_w5(1)),r1,[la,rpd])
#else
                    call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
                    do pair = 1,2
                       call dgemm('t','n',rpd,lg*la*nb,nb,p10,xo_pair(1,pair,tid+1),nb,w3,nb,nul,my_w2,rpd)
                       call dgemm('n','n',rpd*lg*la,pnv,nb,p10,my_w2,rpd*lg*la,yv_pair(1,1,tid+1),nb,nul,my_w1,rpd*lg*la)
                       call array_reorder_4d(p10,my_w1,rpd,lg,la,pnv,[3,4,1,2],nul,my_w2)

                       !u{cldj}(j,dlc)
                       call array_reorder_2d(p20,t,pnv,pnv,paircontrib(:,3-pair),nul,my_w4)
                       call array_reorder_2d(m10,t,pnv,pnv,paircontrib(:,pair),  p10,my_w4)

                       !this hopefully optimizes the construction of this term
                       if( pnv < lg )then
                          call dgemm('n','n',la*pnv*rpd,pnv,lg,p10,my_w2,la*pnv*rpd,&
                             &yv_pair(fg,1,tid+1),nb,nul,my_w1,la*pnv*rpd)
                          h1 => my_w1
                          h2 => my_w4
                          contract = pnv*rpd*pnv
                       else
                          call dgemm('n','t',rpd*pnv*rpd,lg,pnv,p10,my_w4,rpd*pnv*rpd,&
                             &yv_pair(fg,1,tid+1),nb,nul,my_w1,rpd*pnv*rpd)
                          h1 => my_w2
                          h2 => my_w1
                          contract = pnv*rpd*lg
                       endif

                       call dgemm('n','t',la,rpd,contract,p10,h1,la,h2,rpd,nul,my_w5,la)


                       !$OMP CRITICAL
                       Gai(fa:xa,idx(3-pair)) = Gai(fa:xa,idx(3-pair)) + r1(:,1)
                       !$OMP END CRITICAL
                    enddo
                 else
                    call dgemm('t','n',rpd,lg*la*nb,nb,p10,xo_pair(1,1,tid+1),nb,w3,nb,nul,my_w2,rpd)
                    call dgemm('n','n',rpd*lg*la,pnv,nb,p10,my_w2,rpd*lg*la,yv_pair(1,1,tid+1),nb,nul,my_w1,rpd*lg*la)
                    call array_reorder_4d(p10,my_w1,rpd,lg,la,pnv,[3,4,1,2],nul,my_w2)

                    !u{cldj}(j,dlc)
                    call array_reorder_4d(p20,t,pnv,rpd,pnv,rpd,[4,3,2,1],nul,my_w4)
                    call array_reorder_4d(m10,t,pnv,rpd,pnv,rpd,[2,3,4,1],p10,my_w4)

                    !this hopefully optimizes the construction of this term
                    if( pnv < lg )then
                       call dgemm('n','n',la*pnv*rpd,pnv,lg,p10,my_w2,la*pnv*rpd,yv_pair(fg,1,tid+1),nb,nul,my_w1,la*pnv*rpd)
                       h1 => my_w1
                       h2 => my_w4
                       contract = pnv*rpd*pnv
                    else
                       call dgemm('n','t',rpd*pnv*rpd,lg,pnv,p10,my_w4,rpd*pnv*rpd,yv_pair(fg,1,tid+1),nb,nul,my_w1,rpd*pnv*rpd)
                       h1 => my_w2
                       h2 => my_w1
                       contract = pnv*rpd*lg
                    endif

                    call dgemm('n','t',la,rpd,contract,p10,h1,la,h2,rpd,nul,my_w5,la)

#ifdef VAR_PTR_RESHAPE
                    r1(1:la,1:rpd) => my_w5
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                    call c_f_pointer(c_loc(my_w5(1)),r1,[la,rpd])
#else
                    call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
                    !$OMP CRITICAL
                    do ic = 1, rpd
                       Gai(fa:xa,idx(ic)) = Gai(fa:xa,idx(ic)) + r1(:,ic)
                    enddo
                    !$OMP END CRITICAL

                 endif
#ifndef VAR_OMP
                 call time_start_phase(PHASE_WORK, at = times_in_loops(2,2) )
#endif
              endif

           enddo
           !$OMP END DO NOWAIT
           !$OMP END PARALLEL


           call time_start_phase(PHASE_WORK, ttot = times(2), twall = times(1) )
           times(7) = times(7) + times(2)

           !ADD THE gvvvv contribution, note that we destroy the original
           !integrals in w1 and replace them by those with restricted indices in
           !alpha and gamma -- the following code is copied from
           !get_a22_and_prepb22_terms_ex, where more info can be found
           if(fa<=fg+lg-1)then

              goffs=0
              if(fa-fg>0)goffs=fa-fg
              aoffs=0
              if(fg-fa>0)aoffs=fg-fa
              tred=0
              tlen=min(min(min(la,lg),fg+lg-fa),fa+la-fg)

              if(fa+la-1>=fg)tred= tlen*(tlen+1)/2
              if(fa>=fg.and.fg+lg-fa-tlen>0)tred=tred+(fg+lg-fa-tlen)*la
              if(fa<fg)tred=tred+lg*aoffs
              if(fa<fg.and.fa+la>fg) tred=tred+(lg-tlen)*(la-aoffs)
              if(tlen<=0)then
                 tlen=0
                 aoffs=0
                 goffs=0
                 tred=la*lg
              endif


              if(.not.I_PLUS_MINUS_DONE)then

                 if(this_is_query)then
                    query%size_array(w1_tag) = max(query%size_array(w1_tag),var_inp(2))
                    var_inp = 0

                    !note that the mem requirements for w4 and w5 are the same
                    query%size_array(w1_tag) = max(query%size_array(w1_tag),(i8*nb*nb)*la*lg)
                    query%size_array(w3_tag) = max(query%size_array(w3_tag),(i8*nb*nb)*tred)
                    query%size_array(w4_tag) = max(query%size_array(w4_tag),(i8*la*nb)*lg*nb)
                    query%size_array(w5_tag) = max(query%size_array(w5_tag),(i8*la*nb)*lg*nb)
                 else
                    call array_reorder_4d(p10,w3,nb,lg,la,nb,[3,4,2,1],nul,w1)
                    ! (w3): I[beta delta alpha gamma] <= (w1): I[alpha beta gamma delta]
                    call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w3)

                    call get_I_plusminus_le(w3,'+',fa,fg,la,lg,nb,tlen,tred,goffs,s3)
                    !(w4): I+ [delta alpha<=gamma beta] <= (w3): I+ [beta delta alpha<=gamma]
                    call array_reorder_3d(1.0E0_realk,w3,nb,nb,tred,[2,3,1],0.0E0_realk,w4)
                    ! (w2): I[beta delta alpha gamma] <= (w1): I[alpha beta gamma delta]
                    call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w3)

                    call get_I_plusminus_le(w3,'-',fa,fg,la,lg,nb,tlen,tred,goffs,s3)

                    !(w5): I- [delta alpha<=gamma beta] <= (w3): I+ [beta delta alpha<=gamma]
                    call array_reorder_3d(1.0E0_realk,w3,nb,nb,tred,[2,3,1],0.0E0_realk,w5)
                 endif

                 I_PLUS_MINUS_DONE = .true.

              endif

              call time_start_phase(PHASE_WORK, at = times(8), twall = times(1) )

              !OMP PARALLEL  DEFAULT(NONE) IF(this_is_not_query)&
              !OMP SHARED(nspaces,paircontrib,xa,xg,Gai,goffs,aoffs,tlen,tred,&
              !OMP pno_cv,pno_t2,pno_o2,offset_for_omp,no,nv,nb,xo,xv,yo,yv,tmi,tpl,&
              !OMP fa,la,fg,lg,w1,w2,w3,w4,w5,query,this_is_query,this_is_not_query,&
              !OMP sio4,info_omp1,xv_pair,xo_pair,yv_pair,yo_pair,nthreads_level1_int_dir) &
              !OMP PRIVATE(a,d,t,idx,pnv,pno,rpd,PS,o,ns,i,h1,h2,contract,&
              !OMP pno_comb,beg1,beg2,nor,nvr,my_w1,my_w2,my_w3,my_s1,my_s2,my_s3,my_s4,my_s5,&
              !OMP my_w4,my_w5,tid,EOS,edit,h3) REDUCTION(+:times,times_in_loops,tc,tw)&
              !OMP NUM_THREADS(nthreads_level1_int_dir(third_loop))
              tid = 0
#ifdef VAR_OMP
              tid = omp_get_thread_num()
#endif
              if(this_is_not_query)then

                 my_w1 => w1(info_omp1(beg_array,tid+1,w1_tag,third_loop):info_omp1(end_array,tid+1,w1_tag,third_loop))
                 my_w2 => w2(info_omp1(beg_array,tid+1,w2_tag,third_loop):info_omp1(end_array,tid+1,w2_tag,third_loop))
                 my_w3 => w3(info_omp1(beg_array,tid+1,w3_tag,third_loop):info_omp1(end_array,tid+1,w3_tag,third_loop))

                 my_s1 = info_omp1(end_array,tid+1,w1_tag,third_loop)-info_omp1(beg_array,tid+1,w1_tag,third_loop)+1
                 my_s2 = info_omp1(end_array,tid+1,w2_tag,third_loop)-info_omp1(beg_array,tid+1,w2_tag,third_loop)+1
                 my_s3 = info_omp1(end_array,tid+1,w3_tag,third_loop)-info_omp1(beg_array,tid+1,w3_tag,third_loop)+1

                 my_w4 => null()
                 my_w5 => null()

                 my_s4 = 0
                 my_s5 = 0

              endif

              !OMP DO SCHEDULE(DYNAMIC)
              do ns = 1, nspaces

                 if(.not.pno_cv(ns)%allocd)then

                    cycle 

                 endif


                 d   => pno_cv(ns)%d
                 t   => pno_t2(ns)%elm1
                 idx => pno_cv(ns)%iaos
                 pnv =  pno_cv(ns)%ns2
                 pno =  pno_cv(ns)%n
                 rpd =  pno_cv(ns)%rpd
                 PS  =  pno_cv(ns)%PS

                 o   => pno_o2(ns)%elm1

                 nor=rpd*(rpd+1)/2
                 nvr=pnv*(pnv+1)/2

                 if(this_is_query)then

                    edit = offset_for_omp+(ns-1)*nloops*5+(third_loop-1)*5

                    if(.not.PS)then
                       query%size_array(edit+w2_tag) = max(query%size_array(edit+w2_tag),(i8*pnv*rpd)*pnv*rpd)
                    endif

                    !The sizes are the same for both blocks + an -
                    query%size_array(edit+w2_tag) = max(query%size_array(edit+w2_tag),(i8*nb*tred)*pnv)
                    query%size_array(edit+w1_tag) = max(query%size_array(edit+w1_tag),(i8*tred*pnv)*pnv)
                    query%size_array(edit+w2_tag) = max(query%size_array(edit+w2_tag),(i8*tred)*nvr)
                    query%size_array(edit+w3_tag) = max(query%size_array(edit+w3_tag),(i8*nor)*tred*2)

                    my_w1 => dummy
                    my_w2 => dummy
                    my_w3 => dummy

                    my_s1 = 1
                    my_s2 = 1
                    my_s3 = 1
                    my_s4 = 0
                    my_s5 = 0
                 else

                    !Get the transformation matrices
                    !  xo_pair 
                    do i = 1, pno
                       do a = 1, nb 
                          xo_pair(a,i,tid+1) = xo(a,idx(i))
                       enddo
                    enddo
                    !  xv_pair
                    call dgemm('n','n',nb,pnv,nv,p10,xv,nb,d,nv,nul,xv_pair(1,1,tid+1),nb)
                    !  yo_pair
                    do i = 1, pno
                       do a = 1, nb
                          yo_pair(a,i,tid+1) = yo(a,idx(i))
                       enddo
                    enddo
                    !  yv_pair
                    call dgemm('n','n',nb,pnv,nv,p10,yv,nb,d,nv,nul,yv_pair(1,1,tid+1),nb)

                    call time_start_phase(PHASE_WORK, at = times_in_loops(1,3) )

                    ! get tpl and tmi
                    if( PS )then
                       h1 => t
                    else
                       call array_reorder_4d(p10,t,pnv,rpd,pnv,rpd,[1,3,2,4],nul,my_w2)
                       h1 => w2(1:pnv**2*rpd**2)
                    endif
                    call get_tpl_and_tmi_fort(h1,pnv,rpd,tpl(:,tid+1),tmi(:,tid+1))

                    call time_start_phase(PHASE_WORK, at = times_in_loops(2,3) )

                    !Plus block
                    call dgemm('n','n',nb*tred,pnv,nb,p10,w4,nb*tred,yv_pair(1,1,tid+1),nb,nul,my_w2,nb*tred)
                    !(w0):I+ [alpha<=gamma c d] = (w2):I+ [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
                    call dgemm('t','n',tred*pnv,pnv,nb,p10,my_w2,nb,yv_pair(1,1,tid+1),nb,nul,my_w1,pnv*tred)
                    !(w2):I+ [alpha<=gamma c>=d] <= (w0):I+ [alpha<=gamma c d] 
                    call get_I_cged(my_w2,my_w1,tred,pnv)
                    !(w3.1):sigma+ [alpha<=gamma i>=j] = (w2):I+ [alpha<=gamma c>=d] * (w0):t+ [c>=d i>=j]
                    call dgemm('n','n',tred,nor,nvr,p05,my_w2,tred,tpl(1,tid+1),nvr,nul,my_w3,tred)


                    !Minus block
                    call dgemm('n','n',nb*tred,pnv,nb,p10,w5,nb*tred,yv_pair(1,1,tid+1),nb,nul,my_w2,nb*tred)
                    !(w0):I- [alpha<=gamma c d] = (w2):I- [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
                    call dgemm('t','n',tred*pnv,pnv,nb,p10,my_w2,nb,yv_pair(1,1,tid+1),nb,nul,my_w1,pnv*tred)
                    !(w2):I- [alpha<=gamma c<=d] <= (w0):I- [alpha<=gamma c d] 
                    call get_I_cged(my_w2,my_w1,tred,pnv)
                    !(w3.2):sigma- [alpha<=gamma i<=j] = (w2):I- [alpha<=gamma c>=d] * (w0):t- [c>=d i>=j]
                    call dgemm('n','n',tred,nor,nvr,p05,my_w2,tred,tmi(1,tid+1),nvr,nul,my_w3(tred*nor+1),tred)

                    call time_start_phase(PHASE_WORK, twall = times(2), at = times_in_loops(3,3) )


#ifdef VAR_PTR_RESHAPE
                    h1(1:nb*no)  => xo
                    h3(1:nb*pnv) => xv_pair
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                    call c_f_pointer(c_loc(xo(1,1)),h1,[(i8*nb)*no])
                    call c_f_pointer(c_loc(xv_pair(1,1,tid+1)),h3,[(i8*nb)*pnv])
#else
                    call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

                 endif

                 if( PS )then
                    call combine_and_transform_sigma(pno_o2(ns),my_w1,my_w2,my_w3,my_s1,my_s2,my_s3,h3,h1,sio4(ns),nor,tlen,&
                       &tred,fa,fg,la,lg,rpd,pnv,nb,goffs,aoffs,4,.false.,tw,tc,rest_occ_om2=.true., act_no = no, &
                       &query = this_is_query )  
                 else
                    call combine_and_transform_sigma(pno_o2(ns),my_w1,my_w2,my_w3,my_s1,my_s2,my_s3,h3,h1,sio4(ns),nor,&
                       &tlen,tred,fa,fg,la,lg,rpd,pnv,nb,goffs,aoffs,4,.false.,tw,tc, &
                       &order=[1,3,2,4], act_no = no, sio4_ilej = .false., query = this_is_query )
                 endif

                 if( this_is_query )then

                    query%size_array(edit+w1_tag) = max(query%size_array(edit+w1_tag),my_s1)
                    query%size_array(edit+w2_tag) = max(query%size_array(edit+w2_tag),my_s2)
                    query%size_array(edit+w3_tag) = max(query%size_array(edit+w3_tag),my_s3)

                    my_w1 => null()
                    my_w2 => null()
                    my_w3 => null()

                 else

                    call time_start_phase(PHASE_WORK, ttot =  times(2))
                    times_in_loops(4,3) =  times_in_loops(4,3) + times(2)
                 endif


              enddo
              !OMP END DO NOWAIT
              !OMP END PARALLEL

              call time_start_phase( PHASE_WORK, ttot = times(1) )
              times(9) = times(9) + times(1)

           endif

        enddo BatchAlpha
     enddo BatchGamma
     
     if(.not.this_is_query)then
        ! Free integral stuff
        ! *******************
        nullify(Mylsitem%setting%LST_GAB_LHS)
        nullify(Mylsitem%setting%LST_GAB_RHS)
        call free_decscreen(DECSCREEN)

        !Free rest
        call mem_dealloc( xo_pair )
        call mem_dealloc( yo_pair )
        call mem_dealloc( xv_pair )
        call mem_dealloc( yv_pair )

        call mem_dealloc( tpl )
        call mem_dealloc( tmi )

        call mem_dealloc(info_omp1)

        call time_start_phase(PHASE_WORK, ttot = tbeg)
        if(DECinfo%PL>2)then
           write (*,'(" PNO: total time batched loop     :",1g10.3)')tbeg
           write (*,'(" PNO: times split up in            ")')
           write (*,'(" PNO: Integral calculation        :",1g10.3)')times(3)
           write (*,'(" PNO: Conventional transformation :",1g10.3)')times(4)
           write (*,'(" PNO: loop 1                      :",4g10.3)')times(5),&
              &times_in_loops(1,1),times_in_loops(2,1),times_in_loops(3,1)
           write (*,'(" PNO: sorting for Gai             :",1g10.3)')times(6)
           write (*,'(" PNO: loop 2                      :",3g10.3)')times(7),&
              &times_in_loops(1,2),times_in_loops(2,2)
           write (*,'(" PNO: sorting for B2 term         :",1g10.3)')times(8)
           write (*,'(" PNO: loop 3                      :",5g10.3)')times(9),&
              &times_in_loops(1,3),times_in_loops(2,3),times_in_loops(3,3),times_in_loops(4,3)
        endif

#ifdef VAR_OMP
        call omp_set_nested(nested)
#endif
     endif

  end subroutine pno_residual_integral_direct_loop


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
        &w1,w2,w3,w4,w5,govov,vvf,sio4,nspaces)
     implicit none
     integer, intent(in) :: no,ns,nspaces
     type(PNOSpaceInfo), intent(inout) :: pno_cv(nspaces),pno_S(nspaces*(nspaces-1)/2)
     type(tensor), intent(in) :: pno_t2(nspaces),sio4(nspaces)
     real(realk),pointer,intent(inout) :: o2_space(:)
     real(realk),pointer,intent(inout) :: w1(:),w2(:),w3(:),w4(:),w5(:)
     real(realk),intent(in),target :: govov(no**2*pno_cv(ns)%ns1**2)
     real(realk),intent(in) :: vvf(:,:)
     character :: tr11,tr12,tr21,tr22
     real(realk),pointer :: p1(:,:,:,:), p2(:,:,:,:), p3(:,:,:,:), p4(:,:,:,:),h1(:), h2(:), r1(:,:),r2(:,:),d(:,:),d1(:,:)
     real(realk),pointer :: o(:),t(:),S1(:,:), t21(:)
     real(realk),pointer :: sig(:,:,:,:),t1(:,:,:)
     logical :: skiptrafo, cyc, PS, PS1
     integer :: space, a,i,b,j, nv,pno,pnv, ns2,k,l,ilej
     integer :: pno1,pnv1,Sidx1,ldS1, rpd1,rpd,pnor, pos1,pos2,pos3
     integer :: pair,paircontrib(2,2)
     integer,pointer :: idx(:),idx1(:)
     integer(kind=8) :: o2v2
     integer, parameter :: paircontribs = 2
     real(realk) :: tw, tc
     integer :: mv(pno_cv(ns)%ns2**2/2), st
     real(realk), parameter :: nul =  0.0E0_realk
     real(realk), parameter :: p10 =  1.0E0_realk
     real(realk), parameter :: p20 =  2.0E0_realk
     real(realk), parameter :: m10 = -1.0E0_realk

     !print *," CALCULATE CONTRIBUTION", ns
     !print *,""
     !print *,""
     !print *,""
     !print *,""
     !print *,""
     !print *,""
     !print *,""



     tw = 0.0E0_realk
     tc = 0.0E0_realk

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

     !jilk = ijkl
#ifdef VAR_PTR_RESHAPE
     sig(1:rpd,1:rpd,1:no,1:no) => sio4(ns)%elm1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
     call c_f_pointer(c_loc(sio4(ns)%elm1(1)),sig,[rpd,rpd,no,no])
#else
     call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

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
#ifdef VAR_PTR_RESHAPE
        p1(1:rpd1,1:nv,1:rpd1,1:nv) => w1
        p2(1:no,  1:nv,1:no,  1:nv) => govov(1:no*nv*no*nv)
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
        call c_f_pointer( c_loc(w1(1)),    p1, [rpd1,nv,rpd1,nv] )
        call c_f_pointer( c_loc(govov(1)), p2, [no,   nv,no, nv] )
#else
        call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

        if( PS1 )then

           do pair = 1, paircontribs


              p1(1,:,1,:) = p2(idx1(paircontrib(1,pair)),:,idx1(paircontrib(2,pair)),:)


              !transform integral contribution govov(ldkc) => govov(\bar{d} k l  \bar{c}) to the space of (ij) -> w2
              call dgemm( 'n', 'n', nv, pnv, nv, p10, w1, nv, d, nv, nul, w2, nv)
              call dgemm( 't', 'n', pnv1, pnv, nv, p10, d1, nv, w2, nv, nul, w1, pnv1 )

              ! Quadratic part of the E2 term use u^{bd}_{kl} (bkdl) as b,dkl
              call array_reorder_2d( p20, t21, pnv1, pnv1, paircontrib(:,3-pair), nul, w3)
              call array_reorder_2d( m10, t21, pnv1, pnv1, paircontrib(:,pair), p10, w3)
              call do_overlap_trafo(nv,ns,ns2,1,pno_cv,pno_S,pnv,pnv1,pnv1,w3,w2,ptr=h1)

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
           call do_overlap_trafo(nv,ns,ns2,1,pno_cv,pno_S,pnv,rpd1*pnv1*rpd1,pnv1,w3,w2,ptr=h1)


           !contract amplitudes in h1 with integrals in w1 and add to w4 : -1 * h1(bdkl) w1(dlkc) += w4(bc)
           call dgemm('n','n',pnv,pnv,pnv1*rpd1**2,m10, h1,pnv,w1,pnv1*rpd1**2,p10,w5,pnv)


        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  B2 Term !!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!

        !prepare 4 occupied integral goooo for B2 term
        ! and Get the integral contribution, sort it first like the integrals then transform it, govov
#ifdef VAR_PTR_RESHAPE
        p1(1:rpd,1:rpd,1:rpd1,1:rpd1) => w3
        p3(1:rpd1,1:nv,1:rpd1,1:nv)   => w1
        p4(1:no,1:nv,1:no,1:nv)       => govov
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
        call c_f_pointer( c_loc(w3(1)),    p1, [rpd,rpd,rpd1,rpd1] )
        call c_f_pointer( c_loc(w1(1)),    p3, [rpd1,nv,rpd1,nv] )
        call c_f_pointer( c_loc(govov(1)), p4, [no,  nv,no,  nv] )
#else
        call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

        !print *," CONRTIB LOOP", ns, ns2
        !print *,""
        !print *,""

        !CODE FOR PAIR SPACES WITH RESTRICTIONS
        if( PS1 )then

           !loop over pair contributions kl and lk
           do pair = 1, paircontribs

              !!p3(1,:,1,:) = p4(idx1(paircontrib(1,pair)),:,idx1(paircontrib(2,pair)),:)

              !!!transform integral contribution
              !!!govov(kcld) => govov(\bar{c} \bar{d} k l) or lk to the space of (ij) -> w2
              !!call dgemm( 'n', 'n', nv, pnv, nv, p10, w1, nv, d, nv, nul, w3, nv )
              !!call dgemm( 't', 'n', pnv, pnv, nv, p10, d, nv, w3, nv, nul, w4, pnv )

              !!if( PS )then
              !!   p1(1,1,1,1) = p2(idx1(paircontrib(1,pair)),idx(1),idx1(paircontrib(2,pair)),idx(2))
              !!else
              !!   do j = 1,pno
              !!      do i = 1,pno
              !!         p1(i,j,1,1) = p2(idx1(paircontrib(1,pair)),idx(i),idx1(paircontrib(2,pair)),idx(j))
              !!      enddo
              !!   enddo
              !!endif

              !!!sort the amplitudes and contract cidj -> ijcd, ijcd cdkl + ijkl = ijkl
              !!call array_reorder_4d( p10, t, pnv, rpd, pnv, rpd, [2,4,1,3], nul, w1 )
              !!call dgemm( 'n', 'n', rpd**2, rpd1**2, pnv**2, p10, w1, rpd**2, w4, pnv**2, p10, w3, rpd**2 )

              !!!contract the B intermediate in w3 with the amplitudes (kl) from the
              !!!inner loop and use the overlap to transform to the omega space, ijkl klab
              !!call array_reorder_2d( p10, t21, pnv1, pnv1 , paircontrib(:,pair), nul, w1 )

              !!!print *,"pair",pair,"contirb", paircontrib(:,3-pair)
              !!call dgemm( 'n', 'n', rpd**2, pnv1**2, rpd1**2, p10, w3, rpd**2, w1, rpd1**2, dble(pair-1), w2, rpd**2 )

              do i = 1, rpd
                 do j = 1, rpd
                    w3(j+(i-1)*pno) = sig(j,i,idx1(paircontrib(1,pair)),idx1(paircontrib(2,pair)))
                 enddo
              enddo


              call array_reorder_2d( p10, t21, pnv1, pnv1, paircontrib(:,pair) , nul, w1 )

              !call print_tensor_unfolding_with_labels(w1,&
              !   &[pnv1,pnv1],'ab',2,[rpd1,rpd1],'lk',2,'AMPS P')

              !call print_tensor_unfolding_with_labels(w3,&
              !   &[rpd,rpd],'ij',2,[rpd1,rpd1],'lk',2,'INTS 1')

              call dgemm( 'n', 't', pnv1**2,rpd*rpd,1, p10, w1, pnv1**2, w3,rpd*rpd, dble(pair-1), w4, pnv1**2 )

           enddo

           !CODE FOR RECTANGULAR SPACE IN \sum_kl
        else

           !!do j=1,pno1
           !!   do i=1,pno1
           !!      p3(i,:,j,:) = p4(idx1(i),:,idx1(j),:)
           !!   enddo
           !!enddo

           !!!transform integral contribution
           !!!govov(kcld) => govov(\bar{c} \bar{d} k l) to the space of (ij) -> w2
           !!call dgemm( 'n', 'n', nv*rpd1**2, pnv, nv, p10, w1, nv*rpd1**2, d, nv, nul, w2, nv*rpd1**2 )
           !!call array_reorder_4d( p10, w2, rpd1, nv, rpd1, pnv, [2,4,1,3], nul, w1 )
           !!call dgemm( 't', 'n', pnv, rpd1**2*pnv, nv, p10, d, nv, w1, nv, nul, w4, pnv )

           !!!prepare 4 occupied integral goooo for B2 term
           !!if( PS )then
           !!   do j=1,pno1
           !!      do i=1,pno1
           !!         p1(1,1,i,j) = p2(idx1(i),idx(1),idx1(j),idx(2))
           !!      enddo
           !!   enddo
           !!else
           !!   do j=1,pno1
           !!      do i=1,pno1
           !!         do b=1,pno
           !!            do a=1,pno
           !!               p1(a,b,i,j) = p2(idx1(i),idx(a),idx1(j),idx(b))
           !!            enddo
           !!         enddo
           !!      enddo
           !!   enddo
           !!endif
           !!!sort the amplitudes and contract cidj -> ijcd, ijcd cdkl + ijkl = ijkl
           !!call array_reorder_4d( p10, t, pnv, rpd, pnv, rpd, [2,4,1,3], nul, w1 )
           !!call dgemm( 'n', 'n', rpd**2, rpd1**2, pnv**2, p10, w1, rpd**2, w4, pnv**2, p10, w3, rpd**2 )

           !!!contract the B intermediate in w3 with the amplitudes (kl) from the
           !!!inner loop and use the overlap to transform to the omega space, ijkl klab
           !!call array_reorder_4d( p10, t21, pnv1, rpd1, pnv1, rpd1, [2,4,1,3], nul, w1 )

           !!call dgemm( 'n', 'n', rpd**2, pnv1**2, rpd1**2, p10, w3, rpd**2, w1, rpd1**2, nul, w2, rpd**2 )

#ifdef VAR_PTR_RESHAPE
           p1(1:rpd,1:rpd,1:pno1,1:pno1) => w3
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer(c_loc(w3(1)),p1,[rpd,rpd,pno1,pno1])
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
           do l=1,pno1
              do k=1,pno1
                 p1(:,:,l,k) = sig(:,:,idx1(l),idx1(k))
              enddo
           enddo

           call array_reorder_4d( p10, t21, pnv1, rpd1, pnv1, rpd1, [3,1,4,2], nul, w1 )

           !call print_tensor_unfolding_with_labels(w1,&
           !   &[pnv1,pnv1],'ba',2,[rpd1,rpd1],'lk',2,'AMPS R')

           !call print_tensor_unfolding_with_labels(w3,&
           !   &[rpd,rpd],'ji',2,[rpd1,rpd1],'lk',2,'INTS R')

           call dgemm( 'n', 't', pnv1**2, rpd**2 , rpd1**2, p10, w1, pnv1**2, w3, rpd**2, nul, w4, pnv1**2 )

        endif

        !call print_tensor_unfolding_with_labels(w4,&
        !   &[pnv1,pnv1],'ba',2,[rpd,rpd],'ij',2,'res - foo')

        !Squareup and transpose subblocks
        !if( .not.PS )then
           !do j=pno,1,-1
           !   do i=j,1,-1
           !      pos1=1+((i+j*(j-1)/2)-1)*pnv1*pnv1
           !      pos2=1+(i-1)*pnv1*pnv1+(j-1)*pno*pnv1*pnv1
           !      if(j/=1) w4(pos2:pos2+pnv1*pnv1-1) = w4(pos1:pos1+pnv1*pnv1-1)
           !   enddo
           !enddo
           !do j=pno,1,-1
           !   do i=j,1,-1
           !      pos1=1+(i-1)*pnv1*pnv1+(j-1)*pno*pnv1*pnv1
           !      pos2=1+(j-1)*pnv1*pnv1+(i-1)*pno*pnv1*pnv1
           !      if(i/=j) w4(pos2:pos2+pnv1*pnv1-1) = w4(pos1:pos1+pnv1*pnv1-1)
           !   enddo
           !enddo

           !do j=pno,1,-1
           !   do i=j,1,-1
           !      pos1=1+(i-1)*pnv1*pnv1+(j-1)*pno*pnv1*pnv1
           !      call alg513(w4(pos1:pnv1*pnv1+pos1-1),pnv1,pnv1,pnv1*pnv1,mv,(pnv1*pnv1)/2,st)
           !   enddo
           !enddo
           !do j=pno,1,-1
           !   do i=j,1,-1
           !      pos1=1+((i+j*(j-1)/2)-1)*pnv1*pnv1
           !      pos2=1+(i-1)*pnv1*pnv1+(j-1)*pno*pnv1*pnv1
           !      pos3=1+(j-1)*pnv1*pnv1+(i-1)*pno*pnv1*pnv1
           !      if(i/=j)then
           !         w2(pos3:pos3+pnv1*pnv1-1) = w4(pos1:pos1+pnv1*pnv1-1)
           !         !call array_reorder_2d(p10,w4(pos1:pos1+pnv1*pnv1-1),pnv1,pnv1,[1,2],nul,w2(pos3:pos3+pnv1*pnv1-1))
           !      endif
           !      call array_reorder_2d(p10,w4(pos1:pos1+pnv1*pnv1-1),pnv1,pnv1,[2,1],nul,w2(pos2:pos2+pnv1*pnv1-1))
           !   enddo
           !enddo
           !call array_reorder_4d(p10,w2,pnv1,pnv1,rpd,rpd,[3,4,1,2],nul,w4)
           !h2 => w4
        !else
           !h2 => w4
        !endif

        ! transform back, or in the case of ns==ns2 just order correctly and add
        ! the contributions directly to the residual, and use symmetry of the
        ! term here
        !call print_tensor_unfolding_with_labels(o,&
        !   &[pnv,rpd],'ai',2,[pnv,rpd],'bj',2,'OMEGA BEFORE')

        if(ns==ns2)then
           call array_reorder_4d( p10, w4, pnv, pnv,rpd,rpd, [2,4,1,3], p10, o )
        else
           call do_overlap_trafo(nv,ns,ns2,1,pno_cv,pno_S,pnv, rpd*pnv1*rpd, pnv1,w4,w1)
           call array_reorder_4d( p10, w1, pnv, pnv1, rpd, rpd, [2,4,1,3], nul, w3 )
           call do_overlap_trafo(nv,ns,ns2,1,pno_cv,pno_S,pnv,rpd*pnv*rpd, pnv1,w3,o,pC=p10)
        endif
        !call print_tensor_unfolding_with_labels(o,&
        !   &[pnv,rpd],'ai',2,[pnv,rpd],'bj',2,'OMEGA AFTER')



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
        &w1,w2,w3,w4,w5,goovv,govov,Lvoov,oof,p_idx,p_nidx,oidx1,oidx2,nspaces)
     implicit none
     integer, intent(in) :: no,ns,nspaces
     type(PNOSpaceInfo), intent(inout) :: pno_cv(nspaces),pno_S(nspaces*(nspaces-1)/2)
     type(tensor), intent(in) :: pno_t2(nspaces)
     real(realk),pointer,intent(inout) :: o2_space(:)
     real(realk),pointer,intent(inout) :: w1(:),w2(:),w3(:),w4(:),w5(:)
     real(realk),intent(in),target :: goovv(no*no*pno_cv(ns)%ns1**2),govov(no*no*pno_cv(ns)%ns1**2)
     real(realk),intent(in),target :: Lvoov(no*no*pno_cv(ns)%ns1**2)
     real(realk),intent(in) :: oof(:,:)
     integer,intent(in)    :: p_idx(:,:),p_nidx(:)
     integer,intent(inout) :: oidx1(:,:),oidx2(:,:)
     character :: tr11,tr12,tr21,tr22,TRamp_pos1(2),TRamp_pos2(2),trh1,trh2
     real(realk),pointer :: p1(:,:,:,:), p2o(:,:,:,:), p2i(:,:,:,:), p3(:,:,:,:), p4(:,:,:,:)
     real(realk),pointer :: p5(:,:,:,:), p6(:,:,:,:)
     real(realk),pointer :: pij(:,:,:,:), pji(:,:,:,:)
     real(realk),pointer :: h1(:), h2(:), h3(:), h4(:),h1i(:),h2i(:)
     real(realk),pointer :: r1(:,:),r2(:,:),d(:,:),d1(:,:),d2(:,:)
     real(realk),pointer :: o(:),t(:),S1(:,:), t21(:), t22(:), f1(:,:,:)
     logical :: skiptrafo,cyc
     integer :: space,nv,rpd,pno,pnv,nc,nc2,nidx1,nidx2,ns1,ns2,kp(2)
     integer :: a,i,ic,b,j,jc,k,kc,l,lc,c
     integer :: rpd1,pno1,pnv1,Sidx1,ldS1
     integer :: rpd2,pno2,pnv2,Sidx2,ldS2
     integer :: pair1,pair2,paircontrib(2,2)
     integer :: ldh1,ldh2,ldh3
     integer :: i_idx,j_idx,pos_in_res,ipos_in_res
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
     !integer, intent(in) :: yep
     !real(realk), intent(in) :: reference(yep,no,yep,no)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     real(realk), pointer :: check_ref(:,:,:,:), phelp(:,:,:,:)
     integer :: diff11,diff12,diff21,diff22, id(2),kcount
     integer :: bpc,epc
     !reduce the no
     logical :: add_contrib,add_contrib1,add_contrib2,FAspace,FAspace1,PS,PS1,PS2
     real(realk) :: pref

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
        !print *,"ENTERING LOOP",ns,ns1,cyc,pno,pno_cv(ns1)%n,pno_cv(ns1)%PS

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
        !up during the loops -> please note that we use the reordered goovv, to
        !avoid cache misses here
#ifdef VAR_PTR_RESHAPE
        p2o(1:nv,1:nv,1:no,1:no) => goovv
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
        call c_f_pointer( c_loc(goovv(1)), p2o, [nv,  nv, no, no] )
#else
        call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
        if( PS )then

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !CODE FOR TRIANGULAR RESIDUAL CONTRIBUTIONS!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !associate pointers such as to write k last in this case, does not
           !change anything for the pno1==2 case but helps a lot otherwise
           !store also the part for Pij
#ifdef VAR_PTR_RESHAPE
           p1(1:nv,1:rpd,1:nv,1:rpd1) => w1
           p3(1:nv,1:rpd,1:nv,1:rpd1) => w1(nv*rpd1*rpd*nv+1:)
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer( c_loc(w1(1)),    p1, [nv,rpd,nv,rpd1] )
           call c_f_pointer( c_loc(w1(nv*rpd1*rpd*nv+1)), p3, [nv,rpd,nv,rpd1] )
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

           !FIND OUT WHICH CONTRIBUTION TO ADD AND SKIP THE OTHER, FROM THE
           !BEGINNING bpc and epc are set such that no looping occurs, both are
           !then replaced in the loops according to the memreqs, there should be
           !no case where no loop is required, so the choice of begin pair
           !contribuion (bpc) and end pair contribution (epc) should be valid in
           !the end
           bpc = 3
           epc = 2
           get_pair_C2: do pair1 = 1, paircontribs

              add_contrib = .false.

              i = idx(paircontrib(1,pair1))
              j = idx(paircontrib(2,pair1))

              if( PS1 )then


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
                 k = idx1(1)

              endif

              add_contrib = ((PS1 .and. ((i<j.and.oidx1(1,3)==idx(1)).or.(i>j.and.oidx1(1,3)==idx(2)))  ) .or.ns==ns1 ) &
                 &.or.(.not. PS1 .and. i==oidx1(1,3)) .or. add_contrib

              if(add_contrib     .and.pair1==1)bpc = 1
              if(.not.add_contrib.and.pair1==1)bpc = 2
              if(add_contrib     .and.pair1==2)epc = 2
              if(.not.add_contrib.and.pair1==2)epc = 1

           enddo get_pair_C2

           !LOOP OVER THE CONTRIBUTIONS
           pair_contrib_C2: do pair1 = bpc, epc
              

              i = idx(paircontrib(1,pair1))
              j = idx(paircontrib(2,pair1))

#ifdef VAR_PTR_RESHAPE
              if(pair1==1)then
                 p4(1:nv,1:rpd,1:nv,1:rpd1) => w1
              else if(pair1==2)then
                 p4(1:nv,1:rpd,1:nv,1:rpd1) => w1(nv*rpd1*rpd*nv+1:)
              endif
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              if(pair1==1)then
                 call c_f_pointer( c_loc(w1(1)),    p4, [nv,rpd,nv,rpd1] )
              else if(pair1==2)then
                 call c_f_pointer( c_loc(w1(nv*rpd1*rpd*nv+1)), p4, [nv,rpd,nv,rpd1] )
              endif
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

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

                 p4(:,1,:,1) = p2o(:,:,k,j)

                 ! transform c to \bar(c} of (ki)  and a to \bar{a} of (ij)
                 call dgemm('t','n', pnv, nv, nv,  p10, d, nv, p4, nv, nul, w2, pnv)
                 call dgemm('n','n', pnv, pnv1,nv, p10, w2, pnv, d1, nv, nul, w4, pnv)

              else

                 do kc=1,pno1

                    p4(:,1,:,kc) = p2o(:,:,idx1(kc),j)
                    call dgemm('t','n', pnv, nv, nv, p10, d, nv, p4(1,1,1,kc), nv, nul, w2, pnv)
                    call dgemm('n','n', pnv, pnv1,nv, p10, w2, pnv, d1, nv, nul, w4(1+(kc-1)*pnv*pnv1), pnv)

                 enddo

              endif

              !!THE INNER CONTRACTION LOOP - BUILDING THE C INTERMEDIATE
              !Here we loop over all amplitudes sharing an index with the ns
              !space, wehere l /= j if PS and all single index spaces
              OneIdxSpaceLoop21: do nc2 = 1, p_nidx(ns)

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


                 !FIXME: Cycle conditions simlify
                 if( PS2 )then

                    !ns2 refers to a pair space with restrictions

                    if(ns==ns2)then
                       l = idx2(paircontrib(1,pair1))
                    else
                       l = oidx2(nidx2+diff21+1,3)
                    endif


                    if(      ( ns /= ns1 .and. oidx1(1,3) /= i )&!.and. .not. ns == ns2) &
                       &.or. ( ns /= ns2 .and. oidx2(1,3) /= j )&!.and. .not. ns == ns1) &
                       &.and. .not. (ns == ns1 .and. ns == ns2 ) .and. PS)then
                       !print *,"CYCLE 1  ",ns,ns1,ns2,PS1,PS2,i,j,k,l,"(",idx,")(",idx1,")(",idx2,")",oidx1(1,3),oidx2(1,3)
                       cycle OneIdxSpaceLoop21
                    endif

                 else if( rpd2 == 1 )then

                    !FIXME: ns2 refers to a single space, should be the same as below
                    !and can therefore be removed 

                    l = idx2(1)

                    !cycle if l = i, which can never happen because then j = i
                    !whcih is impossible for pair space. these contributions are
                    !added through the ractangular part of the code
                    if( l == i )then
                       !print *,"CYCLE 2  ",ns,ns1,ns2,PS1,PS2,i,j,k,l,"(",idx,")(",idx1,")(",idx2,")"
                       cycle OneIdxSpaceLoop21
                    endif

                 else

                    !cycle if the overlapping index corresponds to i
                    
                    if((oidx2(1,3) == i .and. ns/=ns2) .or. (oidx1(1,3) /= i.and.ns/=ns1))then
                       cycle OneIdxSpaceLoop21
                    endif

                 endif


                 !Get the integrals kdlc -> ckld and transform c and d to their
                 !corresponding spaces, (iajb -> bija) 
#ifdef VAR_PTR_RESHAPE
                 p2i(1:no,1:nv,1:no,1:nv) => govov
                 p1(1:rpd1,1:rpd2,1:nv,1:nv) => w1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                 call c_f_pointer( c_loc(govov(1)), p2i, [no, nv, no, nv]  )
                 call c_f_pointer( c_loc(w1(1)),    p1, [rpd1,rpd2,nv,nv] )
#else
                 call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
                 if( PS1 ) then

                    if( PS2 )then

                       !find index k, k=j if ns=ns1 else the non-overlapping index
                       p1(1,1,:,:) = p2i(k,:,l,:)

                    else

                       do lc=1,pno2
                          p1(1,lc,:,:) = p2i(k,:,idx2(lc),:)
                       enddo

                    endif

                 else

                    if( PS2 ) then

                       do kc=1,pno1
                          p1(kc,1,:,:) = p2i(idx1(kc),:,l,:)
                       enddo

                    else

                       do lc=1,pno2
                          do kc=1,pno1
                             p1(kc,lc,:,:) = p2i(idx1(kc),:,idx2(lc),:)
                          enddo
                       enddo

                    endif

                 endif

                 ! transform c to \bar(c} in (ki) and d to \bar{d} in (lj)
                 call dgemm('t','t',pnv1,rpd1*rpd2*nv,nv,  p10, d1, nv, w1, rpd1*rpd2*nv, nul, w3, pnv1)
                 call dgemm('n','n',pnv1*rpd1*rpd2,pnv2,nv,p10, w3, pnv1*rpd1*rpd2, d2, nv, nul, w1, pnv1*rpd1*rpd2)

                 !get the amplitudes in the correct order eldj -> ldje transform to a and contract to
                 ! -0.5 w1(ckld) w2(ldja) += w4(ckja)
#ifdef VAR_PTR_RESHAPE
                 p2i(1:pnv2,1:rpd2,1:pnv2,1:rpd2) => t22
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                 call c_f_pointer( c_loc(t22(1)), p2i, [pnv2,rpd2,pnv2,rpd2] )
#else
                 call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
                 if( PS2 )then
                    if(l>j) call array_reorder_2d(p10,t22,pnv2,pnv2,[1,2],nul,w3)
                    if(j>l) call array_reorder_2d(p10,t22,pnv2,pnv2,[2,1],nul,w3)
                    nidx_h2 = 1
                 else
#ifdef VAR_PTR_RESHAPE
                    p3(1:rpd2,1:pnv2,1:nidx_h2,1:pnv2) => w3
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                    call c_f_pointer( c_loc(w3(1)),  p3, [rpd2,pnv2,nidx_h2,pnv2] )
#else
                    call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
                    do b=1,pnv2
                       do jc=1,nidx2
                          do a=1,pnv2
                             do ic=1,pno2
                                p3(ic,b,jc,a) = p2i(a,ic,b,oidx2(jc,2))
                             enddo
                          enddo
                       enddo
                    enddo
                    nidx_h2 = nidx2
                 endif

                 call do_overlap_trafo(nv,ns,ns2,2,pno_cv,pno_S,rpd2*pnv2*nidx_h2,pnv,pnv2,w3,w2,ptr=h1i,ptr2=h2i)

                 call dgemm('n','n', pnv1*rpd1, nidx_h2*pnv, rpd2*pnv2, m05, w1, pnv1*rpd1, h1i, rpd2*pnv2, nul, h2i, pnv1*rpd1)

#ifdef VAR_PTR_RESHAPE
                 p3(1:pnv1,1:rpd1,1:nidx_h2,1:pnv) => h2i
                 p4(1:rpd,1:pnv,1:pnv1,1:rpd1)     => w4
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                 call c_f_pointer( c_loc(h2i(1)), p3, [pnv1,rpd1,nidx_h2,pnv] )
                 call c_f_pointer( c_loc(w4(1)), p4, [rpd,pnv,pnv1,rpd1] )
#else
                 call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif


                 !FIXME: introduce the tensor_reorder_3/2d
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

                 endif

              enddo OneIdxSpaceLoop21

              !get the amplitudes, extract the necessary indices, 
              !reorder dkci -> dick :D transform to current space (bick) and do the contraction,
              !bick ja,ck^T = bija, do the permutation and addition of the contribution

              if( PS1 )then

                 if(k<i)call array_reorder_2d(p10,t21,pnv1,pnv1,[1,2],nul,w1)
                 if(k>i)call array_reorder_2d(p10,t21,pnv1,pnv1,[2,1],nul,w1)

              else

#ifdef VAR_PTR_RESHAPE
                 p5(1:pnv1,1:1   ,1:pnv1, 1:rpd1) => w1
                 p6(1:pnv1,1:rpd1,1:pnv1, 1:rpd1) => t21
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                 call c_f_pointer( c_loc(w1(1)),  p5, [pnv1,1   ,pnv1, rpd1] )
                 call c_f_pointer( c_loc(t21(1)), p6, [pnv1,rpd1,pnv1, rpd1] )
#else
                 call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

                 !call print_tensor_unfolding_with_labels(t21,&
                 !   &[pnv1,pno1],'bi',2,[pnv1,pno1],'ck',2,'t - pair')

                 do kc=1,pno1
                    p5(:,1,:,kc) = p6(:,kc,:,oidx1(1,2))
                 enddo

                 !call print_tensor_unfolding_with_labels(w1,&
                 !   &[pnv1,1],'bi',2,[pnv1,pno1],'ck',2,'t - pair')

              endif


              kcounter1: do kc = 1,rpd1
                 !call print_tensor_unfolding_with_labels(w1(1+(kc-1)*pnv1**2:),&
                 !   &[pnv,1],'bi',2,[pnv1,1],'ck',2,'t before- pair')

                 call do_overlap_trafo(nv,ns,ns1,1,pno_cv, pno_S, pnv,pnv1, pnv1,w1(1+(kc-1)*pnv1**2:),w2,ptr=h1,ptr2=h2)

                 !call print_tensor_unfolding_with_labels(h1,&
                 !   &[pnv,1],'bi',2,[pnv1,1],'ck',2,'t - pair')

                 !call print_tensor_unfolding_with_labels(w4(1+(kc-1)*pnv1*pnv:),&
                 !   &[1,pnv],'ja',2,[pnv1,1],'ck',2,'Int - pair')

                 call dgemm('n','t', pnv, pnv, pnv1, m10, h1, pnv, w4(1+(kc-1)*pnv1*pnv), pnv, nul, w3, pnv)

                 !call print_tensor_unfolding_with_labels(w3,&
                 !   &[pnv,1],'bi',2,[1,pnv],'ja',2,'res - pair')

                 call array_reorder_2d(p10,w3,pnv,pnv,paircontrib(:,3-pair1),p10,o)
                 call array_reorder_2d(p05,w3,pnv,pnv,paircontrib(:,pair1),p10,o)

                 !call print_tensor_unfolding_with_labels(o,&
                 !   &[pnv,1],'bi',2,[1,pnv],'ja',2,'OMEGA - pair')

              enddo kcounter1


           enddo pair_contrib_C2

        else


           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !CODE FOR RECTANGULAR RESIDUAL CONTRIBUTIONS!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef VAR_PTR_RESHAPE
           p1(1:rpd1,1:rpd,1:nv,1:nv) => w1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer( c_loc(w1(1)), p1, [rpd1,rpd,nv,nv] )
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
           if( PS1 )then
#ifdef VAR_LSDEBUG
              if(diff12 /= 1.or.nidx1/=1) then
                 call lsquit("ERROR(get_common_idx_summation_for_current_aibj): this should never occur",-1)
              endif
#endif
              do j=1,pno
                 p1(1,j,:,:) = p2o(:,:,oidx1(nidx1+diff11+1,3),idx(j))
              enddo
           else
              do j=1,pno
                 do k=1,pno1
                    p1(k,j,:,:) = p2o(:,:,idx1(k),idx(j))
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

              !Get the integrals kdlc -> ckld and transform c and d to their
              !corresponding spaces, (iajb -> bija) 
#ifdef VAR_PTR_RESHAPE
              p2i(1:no, 1:nv, 1:no, 1:nv) => govov
              p1(1:rpd1,1:rpd2,1:nv,1:nv) => w1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(govov(1)), p2i, [no, nv, no, nv]  )
              call c_f_pointer( c_loc(w1(1)),    p1, [rpd1,rpd2,nv,nv] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
              if( PS1 ) then
                 if( PS2 )then
                    p1(1,1,:,:) = p2i(oidx1(nidx1+diff11+1,3),:,oidx2(nidx2+diff21+1,3),:)
                 else
                    do lc=1,pno2
                       p1(1,lc,:,:) = p2i(oidx1(nidx1+diff11+1,3),:,idx2(lc),:)
                    enddo
                 endif
              else
                 if( PS2 ) then
                    do k=1,pno1
                       p1(k,1,:,:) = p2i(idx1(k),:,oidx2(nidx2+diff21+1,3),:)
                    enddo
                 else
                    do lc=1,pno2
                       do k=1,pno1
                          p1(k,lc,:,:) = p2i(idx1(k),:,idx2(lc),:)
                       enddo
                    enddo
                 endif
              endif


              ! transform c to \bar(c} in (ki) and d to \bar{d} in (lj)
              call dgemm('t','t',pnv1,rpd1*rpd2*nv,nv,  p10, d1, nv, w1, rpd1*rpd2*nv, nul, w3, pnv1)
              call dgemm('n','n',pnv1*rpd1*rpd2,pnv2,nv,p10, w3, pnv1*rpd1*rpd2, d2, nv, nul, w1, pnv1*rpd1*rpd2)

              !get the amplitudes in the correct order eldj -> ldje transform to a and contract to
              ! -0.5 w1(ckld) w2(ldja) += w4(ckja)
#ifdef VAR_PTR_RESHAPE
              p3(1:rpd2,1:pnv2,1:nidx2,1:pnv2) => w3
              p2i(1:pnv2,1:rpd2,1:pnv2,1:rpd2) => t22
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(w3(1)),  p3, [rpd2,pnv2,nidx2,pnv2] )
              call c_f_pointer( c_loc(t22(1)), p2i, [pnv2,rpd2,pnv2,rpd2] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
              if( PS2 )then
#ifdef VAR_LSDEBUG
                 if(nidx2/=1)call lsquit("ERRROROROROOROROROR",-1)
#endif
                 if( oidx2(nidx2+diff21+1,2) < oidx2(1,2) ) then
                    do b=1,pnv2
                       do a=1,pnv2
                          p3(1,b,1,a) = p2i(a,1,b,1)
                       enddo
                    enddo
                 else
                    p3(1,:,1,:) = p2i(:,1,:,1)
                 endif
              else
                 do b=1,pnv2
                    do j=1,nidx2
                       do a=1,pnv2
                          do i=1,pno2
                             p3(i,b,j,a) = p2i(a,i,b,oidx2(j,2))
                          enddo
                       enddo
                    enddo
                 enddo
              endif

              call do_overlap_trafo(nv,ns,ns2,2,pno_cv,pno_S,rpd2*pnv2*nidx2,pnv,pnv2,w3,w2,ptr=h1,ptr2=h2)

              call dgemm('n','n', pnv1*rpd1, nidx2*pnv, rpd2*pnv2, m05, w1, pnv1*rpd1, h1, rpd2*pnv2, nul, h2, pnv1*rpd1)

#ifdef VAR_PTR_RESHAPE
              p3(1:pnv1,1:rpd1,1:nidx2,1:pnv) => h2
              p4(1:rpd,1:pnv,1:pnv1,1:rpd1)   => w4
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(h2(1)), p3, [pnv1,rpd1,nidx2,pnv] )
              call c_f_pointer( c_loc(w4(1)), p4, [rpd,pnv,pnv1,rpd1] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
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

           enddo OneIdxSpaceLoop22

           !get the amplitudes, extract the necessary indices, 
           !reorder dkci -> dick :D transform to current space (bick) and do the contraction,
           !bick ja,ck^T = bija, do the permutation and addition of the contribution

#ifdef VAR_PTR_RESHAPE
           p2o(1:pnv1,1:rpd1,1:pnv1,1:rpd1) => t21
           p1(1:pnv1,1:nidx1,1:pnv1,1:rpd1) => w1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer( c_loc(t21(1)), p2o, [pnv1,rpd1,pnv1, rpd1] )
           call c_f_pointer( c_loc(w1(1)),  p1, [pnv1,nidx1,pnv1,rpd1] )
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

           if( PS1 )then

              k = oidx1(nidx1+diff11+1,3)
              i = oidx1(1,3)

              if(k>i)call array_reorder_2d(p10,t21,pnv1,pnv1,[2,1],nul,p1)
              if(k<i)call array_reorder_2d(p10,t21,pnv1,pnv1,[1,2],nul,p1)

           else

              do k=1,pno1
                 do i=1,nidx1
                    p1(:,i,:,k) = p2o(:,k,:,oidx1(i,2))
                 enddo
              enddo
              !call print_tensor_unfolding_with_labels(t21,&
              !       &[pnv1,pno1],'bi',2,[pnv1,pno1],'ck',2,'t - rect')

           endif

           !call print_tensor_unfolding_with_labels(w1,&
           !   &[pnv,nidx1],'bi',2,[pnv1,rpd1],'ck',2,'t before- rect')
           call do_overlap_trafo(nv,ns,ns1,1,pno_cv,pno_S, pnv,nidx1*pnv1*rpd1, pnv1,w1,w2,ptr=h1,ptr2=h2)

           !call print_tensor_unfolding_with_labels(h1,&
           !   &[pnv,nidx1],'bi',2,[pnv1,rpd1],'ck',2,'t - rect')

           !call print_tensor_unfolding_with_labels(w4,&
           !   &[rpd,pnv],'ja',2,[pnv1,rpd1],'ck',2,'Int - rect')

           call dgemm('n','t', pnv*nidx1,rpd*pnv, rpd1*pnv1, m10, h1,pnv*nidx1, w4, pnv*rpd, nul, h2, pnv*nidx1)

           !call print_tensor_unfolding_with_labels(h2,&
           !   &[pnv,nidx1],'bi',2,[rpd,pnv],'ja',2,'res - rect')

#ifdef VAR_PTR_RESHAPE
           p2o(1:pnv,1:nidx1,1:rpd,1:pnv) => h2
           p1(1:pnv,1:rpd, 1:pnv, 1:rpd ) => o
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer( c_loc(h2(1)), p2o, [pnv,nidx1,rpd,pnv ] )
           call c_f_pointer( c_loc(o(1)),  p1,  [pnv,rpd, pnv, rpd ] )
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

           do a=1,pnv
              do j=1,pno
                 do i=1,nidx1
                    do b=1,pnv
                       p1(a,oidx1(i,1),b,j) = p1(a,oidx1(i,1),b,j) + p2o(b,i,j,a)
                       p1(a,j,b,oidx1(i,1)) = p1(a,j,b,oidx1(i,1)) + p05 * p2o(b,i,j,a)
                       p1(b,j,a,oidx1(i,1)) = p1(b,j,a,oidx1(i,1)) + p2o(b,i,j,a)
                       p1(b,oidx1(i,1),a,j) = p1(b,oidx1(i,1),a,j) + p05 * p2o(b,i,j,a)
                    enddo
                 enddo
              enddo
           enddo

           !call print_tensor_unfolding_with_labels(o,&
           !   &[pnv,rpd],'bi',2,[rpd,pnv],'ja',2,'OMEGA - rect')

        endif




        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  D2 Term !!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!

        if( PS )then

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !CODE FOR TRIANGULAR RESIDUAL CONTRIBUTIONS!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           get_pair_D2: do pair1 = 1, paircontribs

              add_contrib = .false.

              i = idx(paircontrib(1,pair1))
              j = idx(paircontrib(2,pair1))

              if( PS1 )then


                 if(ns==ns1)then
                    k = j
                 else
                    k = oidx1(nidx1+diff11+1,3)
                 endif


              else

                 k = idx1(1)

              endif

              add_contrib = ((PS1 .and. ((i<j.and.oidx1(1,3)==idx(2)).or.(i>j.and.oidx1(1,3)==idx(1)))  ) .or.ns==ns1 ) &
                 &.or.(.not. PS1 .and. j==oidx1(1,3)) .or. add_contrib

              if(add_contrib     .and.pair1==1)bpc = 1
              if(.not.add_contrib.and.pair1==1)bpc = 2
              if(add_contrib     .and.pair1==2)epc = 2
              if(.not.add_contrib.and.pair1==2)epc = 1

           enddo get_pair_D2

           !LOOP OVER THE CONTRIBUTIONS
           pair_contrib_D2: do pair1 = bpc, epc

              i = idx(paircontrib(1,pair1))
              j = idx(paircontrib(2,pair1))

              !Similar procedure as for the C2 term, just with the L integrals (which
              !could also be produced on-the-fly to reduce the memory requirements
#ifdef VAR_PTR_RESHAPE
              p2o(1:nv,1:no,1:no,1:nv) => Lvoov
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(Lvoov(1)), p2o, [nv, no, no, nv] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
              ! extract and transform c to \bar(c} of (jk)  and a to \bar{a} of (ij) and reorder
              ! to the in which it will be used later we got w4:\bar{c}k\bar{a}i
              if( PS1 )then
#ifdef VAR_PTR_RESHAPE
                 p1(1:nv,1:rpd,1:rpd1,1:nv) => w1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                 call c_f_pointer( c_loc(w1(1)), p1, [nv,rpd,rpd1,nv] )
#else
                 call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

                 if(ns==ns1)then
                    k = idx1(paircontrib(1,pair1))
                 else
                    k = oidx1(nidx1+diff11+1,3)
                 endif
                 p1(:,1,1,:) = p2o(:,i,k,:)
                 call dgemm('t','n', pnv, nv,   nv, p10,  d,  nv, w1, nv, nul, w2, pnv)
                 call dgemm('n','n', pnv, pnv1, nv, p10, w2, pnv, d1, nv, nul, w4, pnv)

              else
#ifdef VAR_PTR_RESHAPE
                 p1(1:nv,1:rpd,1:nv,1:rpd1) => w1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                 call c_f_pointer( c_loc(w1(1)), p1, [nv,rpd,nv,rpd1] )
#else
                 call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
                 do kc=1,pno1
                    p1(:,1,:,kc) = p2o(:,i, idx1(kc),:)
                    call dgemm('t','n', pnv, nv,   nv, p10,  d,  nv, p1(1,1,1,kc), nv, nul, w2, pnv)
                    call dgemm('n','n', pnv, pnv1, nv, p10, w2, pnv, d1, nv, nul, w4(1+(kc-1)*pnv*pnv1), pnv)
                 enddo
              endif


              OneIdxSpaceLoop31: do nc2=1, p_nidx(ns)
                 ! extract indices:
                 ns2 = p_idx(nc2,ns)

                 call check_if_contributes(ns,ns2,pno_cv,pno_S,cyc)

                 if( cyc )then

                    cycle OneIdxSpaceLoop31

                 endif

                 call get_overlap_idx(ns,ns2,pno_cv,oidx2,nidx2,ndidx1=diff21,ndidx2=diff22)

                 d2   => pno_cv(ns2)%d
                 t22  => pno_t2(ns2)%elm1
                 idx2 => pno_cv(ns2)%iaos
                 pnv2 =  pno_cv(ns2)%ns2
                 pno2 =  pno_cv(ns2)%n
                 rpd2 =  pno_cv(ns2)%rpd
                 PS2  =  pno_cv(ns2)%PS

                 !FIXME: Cycle conditions, same as for C2 term
                 if( PS2 )then

                    !ns2 refers to a pair space with restrictions

                    if(ns==ns2)then
                       l = idx2(paircontrib(2,pair1))
                    else
                       l = oidx2(nidx2+diff21+1,3)
                    endif

                    if(      ( ns /= ns1 .and. oidx1(1,3) /= j )&
                       &.or. ( ns /= ns2 .and. oidx2(1,3) /= i )&
                       &.and. .not. (ns == ns1 .and. ns == ns2 ) .and. PS)then
                    cycle OneIdxSpaceLoop31
                 endif

              else if( rpd2 == 1 )then

                 !ns2 refers to a single space

                 l = idx2(1)

                 if( l == j )then
                    cycle OneIdxSpaceLoop31
                 endif

              else

                 !ns2 refers to a general rectangular space

                 if((oidx2(1,3) == j .and. ns/=ns2) .or. (oidx1(1,3) /= j.and.ns/=ns1))then
                    cycle OneIdxSpaceLoop31
                 endif

              endif


              !Get the L integrals lfkc -> cklf and transform c and d to their
              !corresponding spaces, (iajb -> bjia) 
#ifdef VAR_PTR_RESHAPE
              p1(1:nv,1:rpd1,1:rpd2,1:nv) => w1
              p2i(1:no,1:nv,1:no, 1:nv)   => govov
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(w1(1)),    p1, [nv,rpd1,rpd2,nv] )
              call c_f_pointer( c_loc(govov(1)), p2i, [no, nv, no,  nv] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
              if( PS1 )then
                 if( PS2 )then
                    do a=1,nv
                       do b=1,nv
                          p1(b,1,1,a) = p20 * p2i(l,a,k,b)&
                             & - p2i(k,a,l,b)
                       enddo
                    enddo
                 else
                    do a=1,nv
                       do lc=1,pno2
                          do b=1,nv
                             p1(b,1,lc,a) = p20 * p2i(idx2(lc),a,k,b)&
                                & - p2i(k,a,idx2(lc),b)
                          enddo
                       enddo
                    enddo
                 endif
              else
                 if( PS2 )then
                    do a=1,nv
                       do kc=1,pno1
                          do b=1,nv
                             p1(b,kc,1,a) = p20 * p2i(l,a,idx1(kc),b) &
                                &- p2i(idx1(kc),a,l,b)
                          enddo
                       enddo
                    enddo
                 else
                    do a=1,nv
                       do lc=1,pno2
                          do kc=1,pno1
                             do b=1,nv
                                p1(b,kc,lc,a) = p20 * p2i(idx2(lc),a,idx1(kc),b) - p2i(idx1(kc),a,idx2(lc),b)
                             enddo
                          enddo
                       enddo
                    enddo
                 endif
              endif


              ! transform c to \bar(c} in (ki) and f to \bar{d} in (lj)
              call dgemm('t','n',pnv1,rpd1*rpd2*nv,nv,  p10, d1, nv, w1, nv, nul, w3, pnv1)
              call dgemm('n','n',pnv1*rpd1*rpd2,pnv2,nv,p10, w3, pnv1*rpd1*rpd2, d2, nv, nul, w1, pnv1*rpd1*rpd2)

              !get the u amplitudes in the order eifl -> eifl  transform e to a, reorder to lfai and contract to
              ! -0.5 w1(cklf) h1(lfai) += w4(ckai)
              if( PS2 )then
                 if( l < i ) then
                    call array_reorder_2d(p20,t22,pnv2,pnv2,[2,1],nul,w3)
                    call array_reorder_2d(m10,t22,pnv2,pnv2,[1,2],p10,w3)
                 else
                    call array_reorder_2d(p20,t22,pnv2,pnv2,[1,2],nul,w3)
                    call array_reorder_2d(m10,t22,pnv2,pnv2,[2,1],p10,w3)
                 endif
                 nidx_h2 = 1
              else

#ifdef VAR_PTR_RESHAPE
                 p3(1:pnv2,1:nidx2,1:pnv2,1:rpd2) => w3
                 p2i(1:pnv2,1:rpd2,1:pnv2,1:rpd2) => t22
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                 call c_f_pointer( c_loc(w3(1)),  p3, [pnv2,nidx2,pnv2,rpd2] )
                 call c_f_pointer( c_loc(t22(1)), p2i, [pnv2,rpd2, pnv2,rpd2] )
#else
                 call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

                 do jc=1,pno2
                    do b=1,pnv2
                       do ic=1,nidx2
                          do a=1,pnv2
                             p3(a,ic,b,jc) = p20 * p2i(a,oidx2(ic,2),b,jc) - p2i(a,jc,b,oidx2(ic,2))
                          enddo
                       enddo
                    enddo
                 enddo
                 nidx_h2 = nidx2
              endif

              call do_overlap_trafo(nv,ns,ns2,1,pno_cv,pno_S, pnv, nidx_h2*rpd2*pnv2, pnv2,w3,w2,ptr=h1,ptr2=h2)

              call array_reorder_4d( p10, h1, pnv, nidx_h2, pnv2, rpd2, [4,3,1,2], nul, h2 )

              !call print_tensor_unfolding_with_labels(h2,&
              !   &[rpd2,pnv2],'lf',2,[pnv,nidx_h2],'ai',2,'u2 -- inner -- pair')

              !call print_tensor_unfolding_with_labels(w4,&
              !   &[pnv1,rpd1],'ck',2,[rpd2,pnv2],'lf',2,'Int -- inner -- pair')

              call dgemm('t','t', pnv*nidx_h2,pnv1*rpd1, rpd2*pnv2, p05, h2, rpd2*pnv2, w1, pnv1*rpd1, p10, w4, pnv*nidx_h2)



           enddo OneIdxSpaceLoop31

           !exctract amplitudes as u bjck and contract with w4 ckai

           if( PS1 )then

              if(j>k)then
                 call array_reorder_2d(p20,t21,pnv1,pnv1,[2,1],nul,w1)
                 call array_reorder_2d(m10,t21,pnv1,pnv1,[1,2],p10,w1)
              endif
              if(j<k)then
                 call array_reorder_2d(p20,t21,pnv1,pnv1,[1,2],nul,w1)
                 call array_reorder_2d(m10,t21,pnv1,pnv1,[2,1],p10,w1)
              endif

           else

#ifdef VAR_PTR_RESHAPE
              p5(1:pnv1,1:1,1:pnv1,1:rpd1)    => w1
              p6(1:pnv1,1:rpd1,1:pnv1,1:rpd1) => t21
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(w1(1)),  p5, [pnv1,1,pnv1,rpd1] )
              call c_f_pointer( c_loc(t21(1)), p6, [pnv1,rpd1,pnv1,rpd1] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

              do kc=1,pno1
                 p5(:,1,:,kc) = p20 * p6(:,oidx1(1,2),:,kc) - p6(:,kc,:,oidx1(1,2))
              enddo

           endif

           kcounter2: do kc = 1, rpd1

              call do_overlap_trafo(nv,ns,ns1,1,pno_cv,pno_S,pnv,nidx1*pnv1*rpd1,pnv1,w1(1+(kc-1)*pnv1**2:),w2,ptr=h1,ptr2=h2)

              !call print_tensor_unfolding_with_labels(h1,&
              !   &[pnv,1],'bj',2,[pnv1,1],'ck',2,'u2 pair')

              !call print_tensor_unfolding_with_labels(w4,&
              !   &[pnv1,1],'ai',2,[pnv1,1],'ck',2,'Int pair')

              call dgemm('n','t', pnv, pnv, pnv1,p05,w4(1+(kc-1)*pnv1*pnv),pnv,h1, pnv,nul,w3,pnv)

              !call print_tensor_unfolding_with_labels(h2,&
              !   &[pnv,1],'ai',2,[pnv,1],'ck',2,'res pair')

              !add D2 contribution to o

              !call print_tensor_unfolding_with_labels(o,&
              !   &[pnv,rpd],'ai',2,[pnv,rpd],'bj',2,'OMEGA pair before add')

              call array_reorder_2d(p10,w3,pnv,pnv,paircontrib(:,pair1),p10,o)

              !call print_tensor_unfolding_with_labels(o,&
              !   &[pnv,rpd],'ai',2,[pnv,rpd],'bj',2,'OMEGA pair after add')

           enddo kcounter2
        enddo pair_contrib_D2

     else


           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !CODE FOR RECTANGULAR RESIDUAL CONTRIBUTIONS!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !Similar procedure as for the C2 term, just with the L integrals (which
           !could also be produced on-the-fly to reduce the memory requirements
#ifdef VAR_PTR_RESHAPE
           p1(1:nv,1:rpd,1:rpd1,1:nv) => w1
           p2o(1:nv,1:no,1:no,1:nv)   => Lvoov
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer( c_loc(w1(1)),    p1,  [nv,rpd,rpd1,nv] )
           call c_f_pointer( c_loc(Lvoov(1)), p2o, [nv, no, no, nv] )
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

           if( PS1 )then
              do ic=1,pno
                 p1(:,ic,1,:) = p2o(:,idx(ic), oidx1(nidx1+diff11+1,3),:)
              enddo
           else
              do jc=1,pno1
                 do ic=1,pno
                    p1(:,ic,jc,:) = p2o(:,idx(ic), idx1(jc),:)
                 enddo
              enddo
           endif

           ! transform c to \bar(c} of (jk)  and a to \bar{a} of (ij) and reorder
           ! to the in which it will be used later we got w4:\bar{a}i\bar{c}k
           call dgemm('t','n', pnv, rpd*rpd1*nv, nv,  p10, d, nv, w1, nv, nul, w2, pnv)
           call dgemm('n','n', pnv*rpd*rpd1, pnv1, nv, p10, w2, pnv*rpd1*rpd, d1, nv, nul, w1, pnv*rpd*rpd1)

           call array_reorder_4d( p10, w1, pnv, rpd, rpd1, pnv1, [1,2,4,3], nul, w4 )

           !call print_tensor_unfolding_with_labels(w4,&
           !  &[pnv1,rpd1],'ck',2,[pnv,rpd],'ai',2,'D_ckai')

           !THE INNER CONTRACTION LOOP - BUILDING THE D INTERMEDIATE
           OneIdxSpaceLoop32: do nc2=1, p_nidx(ns)
              ! extract indices:
              ns2 = p_idx(nc2,ns)

              call check_if_contributes(ns,ns2,pno_cv,pno_S,cyc)

              if( cyc )then

                 cycle OneIdxSpaceLoop32

              endif

              call get_overlap_idx(ns,ns2,pno_cv,oidx2,nidx2,ndidx1=diff21,ndidx2=diff22)

              d2   => pno_cv(ns2)%d
              t22  => pno_t2(ns2)%elm1
              idx2 => pno_cv(ns2)%iaos
              pnv2 =  pno_cv(ns2)%ns2
              pno2 =  pno_cv(ns2)%n
              rpd2 =  pno_cv(ns2)%rpd
              PS2  =  pno_cv(ns2)%PS

              !Get the L integrals lfkc -> cklf and transform c and d to their
              !corresponding spaces, (iajb -> bjia) 
#ifdef VAR_PTR_RESHAPE
              p1(1:nv,1:rpd1,1:rpd2,1:nv) => w1
              p2i(1:no,1:nv,1:no,1:nv)    => govov
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(w1(1)),    p1, [nv,rpd1,rpd2,nv] )
              call c_f_pointer( c_loc(govov(1)), p2i, [no, nv, no,  nv] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
              if( PS1 )then
                 if( PS2 )then
                    do a=1,nv
                       do b=1,nv
                          p1(b,1,1,a) = p20 * p2i(oidx2(nidx2+diff21+1,3),a,oidx1(nidx1+diff11+1,3),b)&
                             & - p2i(oidx1(nidx1+diff11+1,3),a,oidx2(nidx2+diff21+1,3),b)
                       enddo
                    enddo
                 else
                    do a=1,nv
                       do lc=1,pno2
                          do b=1,nv
                             p1(b,1,lc,a) = p20 * p2i(idx2(lc),a,oidx1(nidx1+diff11+1,3),b)&
                                & - p2i(oidx1(nidx1+diff11+1,3),a,idx2(lc),b)
                          enddo
                       enddo
                    enddo
                 endif
              else
                 if( PS2 )then
                    do a=1,nv
                       do kc=1,pno1
                          do b=1,nv
                             p1(b,kc,1,a) = p20 * p2i(oidx2(nidx2+diff21+1,3),a,idx1(kc),b) &
                                &- p2i(idx1(kc),a,oidx2(nidx2+diff21+1,3),b)
                          enddo
                       enddo
                    enddo
                 else
                    do a=1,nv
                       do lc=1,pno2
                          do kc=1,pno1
                             do b=1,nv
                                p1(b,kc,lc,a) = p20 * p2i(idx2(lc),a,idx1(kc),b) - p2i(idx1(kc),a,idx2(lc),b)
                             enddo
                          enddo
                       enddo
                    enddo
                 endif
              endif

              ! transform c to \bar(c} in (ki) and f to \bar{d} in (lj)
              call dgemm('t','n',pnv1,rpd1*rpd2*nv,nv,  p10, d1, nv, w1, nv, nul, w3, pnv1)
              call dgemm('n','n',pnv1*rpd1*rpd2,pnv2,nv,p10, w3, pnv1*rpd1*rpd2, d2, nv, nul, w1, pnv1*rpd1*rpd2)

              !get the u amplitudes in the order eifl -> eifl  transform e to a, reorder to lfai and contract to
              ! -0.5 w1(cklf) h1(lfai) += w4(ckai)
              if( PS2 )then
#ifdef VAR_LSDEBUG
                 if(nidx2/=1)call lsquit("ERRROROROROOROROROR",-1)
#endif
                 if( oidx2(nidx2+diff21+1,2) < oidx2(1,2) ) then
                    call array_reorder_2d(p20,t22,pnv2,pnv2,[2,1],nul,w3)
                    call array_reorder_2d(m10,t22,pnv2,pnv2,[1,2],p10,w3)
                 else
                    call array_reorder_2d(p20,t22,pnv2,pnv2,[1,2],nul,w3)
                    call array_reorder_2d(m10,t22,pnv2,pnv2,[2,1],p10,w3)
                 endif
              else
#ifdef VAR_PTR_RESHAPE
                 p3(1:pnv2,1:nidx2,1:pnv2,1:rpd2) => w3
                 p2i(1:pnv2,1:rpd2,1:pnv2,1:rpd2) => t22
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
                 call c_f_pointer( c_loc(w3(1)),  p3, [pnv2,nidx2,pnv2,rpd2] )
                 call c_f_pointer( c_loc(t22(1)), p2i, [pnv2,rpd2, pnv2,rpd2] )
#else
                 call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
                 do j=1,pno2
                    do b=1,pnv2
                       do i=1,nidx2
                          do a=1,pnv2
                             p3(a,i,b,j) = p20 * p2i(a,oidx2(i,2),b,j) - p2i(a,j,b,oidx2(i,2))
                          enddo
                       enddo
                    enddo
                 enddo
              endif

              call do_overlap_trafo(nv,ns,ns2,1,pno_cv,pno_S, pnv, nidx2*rpd2*pnv2, pnv2,w3,w2,ptr=h1,ptr2=h2)

              call array_reorder_4d( p10, h1, pnv, nidx2, pnv2, rpd2, [4,3,1,2], nul, h2 )

              !call print_tensor_unfolding_with_labels(h2,&
              !   &[rpd2,pnv2],'lf',2,[pnv,nidx2],'ai',2,'u2 -- inner -- rect')

              !call print_tensor_unfolding_with_labels(w4,&
              !   &[pnv1,rpd1],'ck',2,[rpd2,pnv2],'lf',2,'Int -- inner -- rect')

              call dgemm('n','n', pnv1*rpd1, pnv*nidx2, rpd2*pnv2, p05, w1, pnv1*rpd1, h2, rpd2*pnv2, nul, h1, pnv1*rpd1)

              !call print_tensor_unfolding_with_labels(h2,&
              !   &[pnv,rpd1],'ck',2,[pnv,nidx2],'ai',2,'res -- inner -- rect')

#ifdef VAR_PTR_RESHAPE
              p2i(1:pnv1,1:rpd1,1:pnv,1:nidx2) => h1
              p4(1:pnv,1:rpd,1:pnv1,1:rpd1)    => w4
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(h1(1)), p2i, [pnv1,rpd1,pnv,nidx2] )
              call c_f_pointer( c_loc(w4(1)), p4, [pnv,rpd,pnv1,rpd1] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
              if( PS1 )then
                 if( PS2 )then
                    do b=1,pnv
                       do a=1,pnv1
                          p4(b,oidx2(1,1),a,1) = p4(b,oidx2(1,1),a,1) + p2i(a,1,b,1)
                       enddo
                    enddo
                 else
                    do j=1,nidx2
                       do b=1,pnv
                          do a=1,pnv1
                             p4(b,oidx2(j,1),a,1) = p4(b,oidx2(j,1),a,1) + p2i(a,1,b,j)
                          enddo
                       enddo
                    enddo
                 endif
              else
                 if( PS2 )then
                    do b=1,pnv
                       do i=1,pno1
                          do a=1,pnv1
                             p4(b,oidx2(1,1),a,i) = p4(b,oidx2(1,1),a,i) + p2i(a,i,b,1)
                          enddo
                       enddo
                    enddo
                 else
                    do j=1,nidx2
                       do b=1,pnv
                          do i=1,pno1
                             do a=1,pnv1
                                p4(b,oidx2(j,1),a,i) = p4(b,oidx2(j,1),a,i) + p2i(a,i,b,j)
                             enddo
                          enddo
                       enddo
                    enddo
                 endif
              endif

              !call print_tensor_unfolding_with_labels(w4,&
              !   &[pnv,rpd],'ai',2,[pnv1,rpd1],'ck',2,'aft -- inner -- rect')

           enddo OneIdxSpaceLoop32

           !exctract amplitudes as u bjck and contract with w4 ckai
#ifdef VAR_PTR_RESHAPE
           p1(1:pnv1,1:nidx1,1:pnv1,1:rpd1) => w1
           p2o(1:pnv1,1:rpd1,1:pnv1,1:rpd1) => t21
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer( c_loc(w1(1)),  p1, [pnv1,nidx1,pnv1,rpd1] )
           call c_f_pointer( c_loc(t21(1)), p2o, [pnv1,rpd1,pnv1,rpd1] )
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

           if( PS1 )then

              k = oidx1(nidx1+diff11+1,3)
              j = oidx1(1,3)

              if(k<j)then
                 call array_reorder_2d(p20,t21,pnv1,pnv1,[2,1],nul,w1)
                 call array_reorder_2d(m10,t21,pnv1,pnv1,[1,2],p10,w1)
              endif
              if(k>j)then
                 call array_reorder_2d(p20,t21,pnv1,pnv1,[1,2],nul,w1)
                 call array_reorder_2d(m10,t21,pnv1,pnv1,[2,1],p10,w1)
              endif

           else

              do kc=1,pno1
                 do ic=1,nidx1
                    p1(:,ic,:,kc) = p20 * p2o(:,oidx1(ic,2),:,kc) - p2o(:,kc,:,oidx1(ic,2))
                 enddo
              enddo

           endif

           call do_overlap_trafo(nv,ns,ns1,1,pno_cv,pno_S,pnv,nidx1*pnv1*rpd1,pnv1,w1,w2,ptr=h1,ptr2=h2)

           !call print_tensor_unfolding_with_labels(h1,&
           !   &[pnv,nidx1],'bj',2,[pnv1,rpd1],'ck',2,'u2 rect')

           !call print_tensor_unfolding_with_labels(w4,&
           !   &[pnv,rpd],'ai',2,[pnv1,rpd1],'ck',2,'Int rect')


           call dgemm('n','t',pnv*rpd,pnv*nidx1,pnv1*rpd1,p05,w4,pnv*rpd,h1,pnv*nidx1,nul,h2,pnv*rpd)

           !call print_tensor_unfolding_with_labels(h2,&
           !   &[pnv,rpd],'ai',2,[pnv,nidx1],'ck',2,'res rect')

           !add D2 contribution to o
#ifdef VAR_PTR_RESHAPE
           p5(1:pnv,1:rpd,1:pnv,1:nidx1) => h2
           p1(1:pnv,1:rpd,1:pnv,1:rpd)   => o
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer( c_loc(h2(1)), p5, [pnv,rpd,pnv,nidx1] )
           call c_f_pointer( c_loc(o(1)),  p1, [pnv,rpd, pnv, rpd ] )
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

           !call print_tensor_unfolding_with_labels(o,&
           !   &[pnv,rpd],'ai',2,[pnv,rpd],'bj',2,'OMEGA rect before add')

           do jc=1,nidx1
              do b=1,pnv
                 do ic=1,pno
                    do a=1,pnv
                       p1(a,ic,b,oidx1(jc,1)) = p1(a,ic,b,oidx1(jc,1)) + p5(a,ic,b,jc)
                       p1(b,oidx1(jc,1),a,ic) = p1(b,oidx1(jc,1),a,ic) + p5(a,ic,b,jc)
                    enddo
                 enddo
              enddo
           enddo

           !call print_tensor_unfolding_with_labels(o,&
           !  &[pnv,rpd],'ai',2,[pnv,rpd],'bj',2,'OMEGA rect after add')

        endif







        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  E2 Term part 2!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !Similar procedure as for the C2 term, just nothing has to be
        !transformed for the occ-occ fock matrix --> might be constructed
        !outside the loops and saved along with the overlaps, the Foo( k, j )

        if( PS )then


           bpc = 3
           epc = 2
           get_pair_E22: do pair1 = 1, paircontribs
              add_contrib = .false.

              i = idx(paircontrib(1,pair1))
              j = idx(paircontrib(2,pair1))

              if( PS1 )then


                 if(ns==ns1)then
                    k = j
                 else
                    k = oidx1(nidx1+diff11+1,3)
                 endif


              else

                 k = idx1(1)

              endif

              add_contrib = ((PS1 .and. ((i<j.and.oidx1(1,3)==idx(1)).or.(i>j.and.oidx1(1,3)==idx(2)))  ) .or.ns==ns1 ) &
                 &.or.(.not. PS1 .and. i==oidx1(1,3)) .or. add_contrib

              if(add_contrib     .and.pair1==1)bpc = 1
              if(.not.add_contrib.and.pair1==1)bpc = 2
              if(add_contrib     .and.pair1==2)epc = 2
              if(.not.add_contrib.and.pair1==2)epc = 1

           enddo get_pair_E22


           pair_contrib_E22: do pair1 = bpc, epc

              i = idx(paircontrib(1,pair1))
              j = idx(paircontrib(2,pair1))

#ifdef VAR_PTR_RESHAPE
              r1(1:rpd1,1:rpd) => w4
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(w4(1)),  r1, [ rpd1,rpd ] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
              if( PS1 )then
                 if(ns==ns1)then
                    k = idx1(paircontrib(2,pair1))
                 else
                    k = oidx1(nidx1+diff11+1,3)
                 endif

                 r1(1,1) = oof(k,j)
              else

                 do kc=1,pno1
                    r1(kc,1) = oof(idx1(kc), j)
                 enddo

              endif


              !THE INNER CONTRACTION LOOP - BUILDING THE E22 INTERMEDIATE
              !OneIdxSpaceLoop41: do nc2=1, p_nidx(ns)

              !   ! extract indices:
              !   ns2 = p_idx(nc2,ns)

              !   call check_if_contributes(ns,ns2,pno_cv,pno_S,cyc)

              !   if( cyc )then

              !      cycle OneIdxSpaceLoop41

              !   endif

              !   call get_overlap_idx(ns,ns2,pno_cv,oidx2,nidx2,ndidx1=diff21,ndidx2=diff22)

              !   d2   => pno_cv(ns2)%d
              !   t22  => pno_t2(ns2)%elm1
              !   idx2 => pno_cv(ns2)%iaos
              !   pnv2 =  pno_cv(ns2)%ns2
              !   pno2 =  pno_cv(ns2)%n
              !   rpd2 =  pno_cv(ns2)%rpd
              !   PS2  =  pno_cv(ns2)%PS

              !   !FIXME: Cycle conditions simlify
              !   if( PS2 )then

              !      !ns2 refers to a pair space with restrictions

              !      if(ns==ns2)then
              !         l = idx2(paircontrib(1,pair1))
              !      else
              !         l = oidx2(nidx2+1+1,3)
              !      endif


              !      if(      ( ns /= ns1 .and. oidx1(1,3) /= i )&!.and. .not. ns == ns2) &
              !         &.or. ( ns /= ns2 .and. oidx2(1,3) /= j )&!.and. .not. ns == ns1) &
              !         &.and. .not. (ns == ns1 .and. ns == ns2 ) .and. PS)then
              !         cycle OneIdxSpaceLoop41
              !      endif

              !   else if( rpd2 == 1 )then

              !      !FIXME: ns2 refers to a single space, should be the same as below
              !      !and can therefore be removed 

              !      l = idx2(1)

              !      !cycle if l = i, which can never happen because then j = i
              !      !whcih is impossible for pair space. these contributions are
              !      !added through the ractangular part of the code
              !      if( l == i )then
              !         !print *,"CYCLE 2  ",ns,ns1,ns2,PS1,PS2,i,j,k,l,"(",idx,")(",idx1,")(",idx2,")"
              !         cycle OneIdxSpaceLoop41
              !      endif

              !   else

              !      !cycle if the overlapping index corresponds to i

              !      if((oidx2(1,3) == i .and. ns/=ns2) .or. (oidx1(1,3) /= i.and.ns/=ns1))then
              !         cycle OneIdxSpaceLoop41
              !      endif

              !   endif

              !   !Get the integrals g(kdlc) as (dklc) and transform c and d to (lj)
              !   !such that the order klcd is obtained
              !   call ass_D1to4( w1,    p1, [nv,rpd1,rpd2,nv] )
              !   call ass_D1to4( govov, p2i, [no, nv, no,  nv] )
              !   if( PS1 )then
              !      if( PS2 )then
              !         p1(:,1,1,:) = p2i(k,:,l,:)
              !      else
              !         do lc=1,pno2
              !            p1(:,1,lc,:) = p2i(k,:,idx2(lc),:)
              !         enddo
              !      endif
              !   else
              !      if( PS2 )then
              !         do kc=1,pno1
              !            p1(:,kc,1,:) = p2i(idx1(kc),:,l,:)
              !         enddo
              !      else
              !         do lc=1,pno2
              !            do kc=1,pno1
              !               p1(:,kc,lc,:) = p2i(idx1(kc),:,idx2(lc),:)
              !            enddo
              !         enddo
              !      endif
              !   endif
              !   ! transform c to \bar(c} in (lj) and d to \bar{d} in (lj)
              !   call dgemm('n','n',nv*rpd1*rpd2,pnv2,nv,  p10, w1, nv*rpd1*rpd2,d2,nv,nul, w3, nv*rpd1*rpd2)
              !   call dgemm('t','n',rpd1*rpd2*pnv2,pnv2,nv,p10, w3, nv, d2, nv, nul, w1, rpd1*rpd2*pnv2)

              !   !get the u amplitudes in the order cldj -> (lcdj) = 2 t(lcdj) - t(jcdl)
              !   call ass_D1to4( w3,  p3, [rpd2,pnv2,pnv2,nidx2] )
              !   call ass_D1to4( t22, p2i, [pnv2,rpd2, pnv2,rpd2] )
              !   if( PS2 )then
              !      if( l > j ) then
              !         call array_reorder_2d(p20,t22,pnv2,pnv2,[2,1],nul,w3)
              !         call array_reorder_2d(m10,t22,pnv2,pnv2,[1,2],p10,w3)
              !      else
              !         call array_reorder_2d(p20,t22,pnv2,pnv2,[1,2],nul,w3)
              !         call array_reorder_2d(m10,t22,pnv2,pnv2,[2,1],p10,w3)
              !      endif
              !   else
              !      do jc=1,nidx2
              !         do ic=1,pno2
              !            p3(ic,:,:,jc) = p20 * p2i(:,ic,:,oidx2(jc,2)) - p2i(:,oidx2(jc,2),:,ic)
              !         enddo
              !      enddo
              !   endif

              !   call dgemm('n','n', rpd1, nidx2, rpd2*pnv2**2, p10, w1, rpd1, w3, rpd2*pnv2**2, nul, w2, rpd1 )

              !   call ass_D1to2( w2, r2, [rpd1,nidx2] )
              !   call ass_D1to2( w4, r1, [rpd1,rpd] )
              !   if( PS1 )then
              !      r1(1,1) = r1(1,1) + r2(1,1)
              !   else
              !      do kc=1,pno1
              !         r1(kc,1) = r1(kc,1) + r2(kc,1)
              !      enddo
              !   endif

              !enddo OneIdxSpaceLoop41


              !extract amplitudes like in C2 as aibk
#ifdef VAR_PTR_RESHAPE
              p1(1:pnv1,1:nidx1,1:pnv1,1:pno1) => w1
              p2o(1:pnv1,1:rpd1,1:pnv1,1:rpd1) => t21
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(w1(1)),  p1, [pnv1,nidx1,pnv1,pno1] )
              call c_f_pointer( c_loc(t21(1)), p2o, [pnv1,rpd1,pnv1,rpd1] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
              if( PS1 )then
                 if(i<k)call array_reorder_2d(p10,t21,pnv1,pnv1,[1,2],nul,p1)
                 if(i>k)call array_reorder_2d(p10,t21,pnv1,pnv1,[2,1],nul,p1)
              else
                 do kc=1,pno1
                    p1(:,1,:,kc) = p2o(:,oidx1(1,2),:,kc)
                 enddo
              endif

              call do_overlap_trafo(nv,ns,ns1,1,pno_cv,pno_S, pnv,pnv1*rpd1, pnv1 ,w1,w2,ptr=h1,ptr2=h2)

              call dgemm('n','n',pnv*pnv1,rpd,rpd1,m10,h1,pnv*pnv1,w4,rpd1,nul,h2,pnv*pnv1)
              call array_reorder_4d(p10,h2,pnv,1,pnv1,rpd,[3,4,1,2], nul, h1)

              !transform b index to the correct space
              call do_overlap_trafo(nv,ns,ns1,1,pno_cv,pno_S, pnv,rpd*pnv,pnv1,h1,h2,ptr=h1)

              call array_reorder_2d(p10,h1,pnv,pnv,paircontrib(:,3-pair1),p10,o)


           enddo pair_contrib_E22



        else

#ifdef VAR_PTR_RESHAPE
           r1(1:rpd1,1:rpd ) => w4
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer( c_loc(w4(1)),  r1, [ rpd1,rpd ] )
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

           if( PS1 )then
              do jc=1,pno
                 r1(1,jc) = oof(oidx1(nidx1+diff11+1,3), idx(jc))
              enddo
           else
              do jc=1,pno
                 do kc=1,pno1
                    r1(kc,jc) = oof(idx1(kc), idx(jc))
                 enddo
              enddo
           endif

           !call print_tensor_unfolding_with_labels(oof,&
           !   &[no],'k',1,[no],'j',1,'full oof')

           !THE INNER CONTRACTION LOOP - BUILDING THE E22 INTERMEDIATE
           !OneIdxSpaceLoop42: do nc2=1, p_nidx(ns)
           !   ! extract indices:
           !   ns2 = p_idx(nc2,ns)

           !   call check_if_contributes(ns,ns2,pno_cv,pno_S,cyc)

           !   if( cyc )then

           !      cycle OneIdxSpaceLoop42

           !   endif

           !   call get_overlap_idx(ns,ns2,pno_cv,oidx2,nidx2,ndidx1=diff21,ndidx2=diff22)

           !   d2   => pno_cv(ns2)%d
           !   t22  => pno_t2(ns2)%elm1
           !   idx2 => pno_cv(ns2)%iaos
           !   pnv2 =  pno_cv(ns2)%ns2
           !   pno2 =  pno_cv(ns2)%n
           !   rpd2 =  pno_cv(ns2)%rpd
           !   PS2  =  pno_cv(ns2)%PS

           !   !Get the integrals g(kdlc) as (dklc) and transform c and d to (lj)
           !   !such that the order klcd is obtained
           !   call ass_D1to4( w1,    p1, [nv,rpd1,rpd2,nv] )
           !   call ass_D1to4( govov, p2i, [no, nv, no,  nv] )
           !   if( PS1 )then
           !      if( PS2 )then
           !         p1(:,1,1,:) = p2i(oidx1(nidx1+diff11+1,3),:,oidx2(nidx2+diff21+1,3),:)
           !      else
           !         do jc=1,pno2
           !            p1(:,1,jc,:) = p2i(oidx1(nidx1+diff11+1,3),:,idx2(jc),:)
           !         enddo
           !      endif
           !   else
           !      if( PS2 ) then
           !         do ic=1,pno1
           !            p1(:,ic,1,:) = p2i(idx1(ic),:,oidx2(nidx2+diff21+1,3),:)
           !         enddo
           !      else
           !         do jc=1,pno2
           !            do ic=1,pno1
           !               p1(:,ic,jc,:) = p2i(idx1(ic),:,idx2(jc),:)
           !            enddo
           !         enddo
           !      endif
           !   endif
           !   ! transform c to \bar(c} in (lj) and d to \bar{d} in (lj)
           !   call dgemm('n','n',nv*rpd1*rpd2,pnv2,nv,  p10, w1, nv*rpd1*rpd2,d2,nv,nul, w3, nv*rpd1*rpd2)
           !   call dgemm('t','n',rpd1*rpd2*pnv2,pnv2,nv,p10, w3, nv, d2, nv, nul, w1, rpd1*rpd2*pnv2)

           !   !get the u amplitudes in the order cldj -> (lcdj) = 2 t(lcdj) - t(jcdl)
           !   call ass_D1to4( w3,  p3, [rpd2,pnv2,pnv2,nidx2] )
           !   call ass_D1to4( t22, p2i, [pnv2,rpd2, pnv2,rpd2] )
           !   if( PS2 )then
           !      if( oidx2(nidx2+diff21+1,2) > oidx2(1,2) ) then
           !         call array_reorder_2d(p20,t22,pnv2,pnv2,[2,1],nul,w3)
           !         call array_reorder_2d(m10,t22,pnv2,pnv2,[1,2],p10,w3)
           !      else
           !         call array_reorder_2d(p20,t22,pnv2,pnv2,[1,2],nul,w3)
           !         call array_reorder_2d(m10,t22,pnv2,pnv2,[2,1],p10,w3)
           !      endif
           !   else
           !      do jc=1,nidx2
           !         do ic=1,pno2
           !            p3(ic,:,:,jc) = p20 * p2i(:,ic,:,oidx2(jc,2)) - p2i(:,oidx2(jc,2),:,ic)
           !         enddo
           !      enddo
           !   endif

           !   call dgemm('n','n', rpd1, nidx2, rpd2*pnv2**2, p10, w1, rpd1, w3, rpd2*pnv2**2, nul, w2, rpd1 )

           !   call ass_D1to2( w2, r2, [rpd1,nidx2] )
           !   call ass_D1to2( w4, r1, [rpd1,rpd] )
           !   if( PS1 )then
           !      if( PS2 )then
           !         r1(1,oidx2(1,1)) = r1(1,oidx2(1,1)) + r2(1,1)
           !      else
           !         do j=1,nidx2
           !            r1(1,oidx2(j,1)) = r1(1,oidx2(j,1)) + r2(1,j)
           !         enddo
           !      endif
           !   else
           !      if( PS2 )then
           !         do i=1,pno1
           !            r1(i,oidx2(1,1)) = r1(i,oidx2(1,1)) + r2(i,1)
           !         enddo
           !      else
           !         do j=1,nidx2
           !            do i=1,pno1
           !               r1(i,oidx2(j,1)) = r1(i,oidx2(j,1)) + r2(i,j)
           !            enddo
           !         enddo
           !      endif
           !   endif

           !enddo OneIdxSpaceLoop42

           !extract amplitudes like in C2 as aibk
           if( PS1 )then
              k = oidx1(nidx1+diff11+1,3)
              i = oidx1(1,3)
              if(i<k)call array_reorder_2d(p10,t21,pnv1,pnv1,[1,2],nul,w1)
              if(i>k)call array_reorder_2d(p10,t21,pnv1,pnv1,[2,1],nul,w1)
           else
#ifdef VAR_PTR_RESHAPE
              p1(1:pnv1,1:nidx1,1:pnv1,1:pno1) => w1
              p2o(1:pnv1,1:pno1,1:pnv1,1:pno1) => t21
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
              call c_f_pointer( c_loc(w1(1)),  p1, [pnv1,nidx1,pnv1,pno1] )
              call c_f_pointer( c_loc(t21(1)), p2o, [pnv1,pno1,pnv1,pno1] )
#else
              call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
              do j=1,pno1
                 do i=1,nidx1
                    p1(:,i,:,j) = p2o(:,oidx1(i,2),:,j)
                 enddo
              enddo
           endif

           call do_overlap_trafo(nv,ns,ns1,1,pno_cv,pno_S, pnv,nidx1*pnv1*rpd1, pnv1 ,w1,w2,ptr=h1,ptr2=h2)

           !call print_tensor_unfolding_with_labels(h1,&
           !   &[pnv,nidx1,pnv],'aib',3,[rpd1],'k',1,'t2 rect')

           !call print_tensor_unfolding_with_labels(w4,&
           !   &[rpd1],'k',1,[rpd],'j',1,'F rect')

           call dgemm('n','n',pnv*nidx1*pnv1,rpd,rpd1,m10,h1,pnv*nidx1*pnv1,w4,rpd1,nul,h2,pnv*nidx1*pnv1)

           !call print_tensor_unfolding_with_labels(h2,&
           !   &[pnv,nidx1,pnv],'aib',3,[rpd],'j',1,'res rect')

           call array_reorder_4d(p10,h2,pnv,nidx1,pnv1,rpd,[3,4,1,2], nul, h1)

           !transform b index to the correct space
           call do_overlap_trafo(nv,ns,ns1,1,pno_cv,pno_S, pnv,rpd*pnv*nidx1,pnv1,h1,h2,ptr=h1)

#ifdef VAR_PTR_RESHAPE
           p2o(1:pnv,1:rpd,1:pnv,1:nidx1) => h1
           p1(1:pnv,1:rpd,1:pnv,1:rpd)    => o
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer( c_loc(h1(1)), p2o, [pnv,rpd,pnv,nidx1] )
           call c_f_pointer( c_loc(o(1)),  p1, [pnv,rpd, pnv, rpd] )
#else
           call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
           do i=1,nidx1
              do a=1,pnv
                 do j=1,pno
                    do b=1,pnv
                       p1(a,oidx1(i,1),b,j) = p1(a,oidx1(i,1),b,j) + p2o(b,j,a,i)
                       p1(b,j,a,oidx1(i,1)) = p1(b,j,a,oidx1(i,1)) + p2o(b,j,a,i)
                    enddo
                 enddo
              enddo
           enddo
        endif

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


  subroutine successive_4ao_mo_trafo_exch(ao,WXYZ,WW,w,XX,x,YY,y,ZZ,z,WRKWXYZ1,WRKWXYZ2,fa,la,fg,lg,res,query,is_this_query)
     implicit none
     integer, intent(in) :: ao,w,x,y,z,fa,la,fg,lg
     real(realk), intent(inout) :: WXYZ(:),WRKWXYZ1(:),WRKWXYZ2(:),res(w*x*y*z)
     real(realk), intent(in) :: WW(ao,w),XX(ao,x),YY(ao,y),ZZ(ao,z)
     integer(kind=8), optional, intent(inout) :: query(3)
     logical, optional, intent(in) :: is_this_query
     integer :: split_ao1, parts1, fe1, le1, leng
     logical :: qu

     qu = .false.
     if(present(query).and.present(is_this_query))qu = is_this_query
     !leng   = ao
     !parts1 = ao/leng

     if(qu)then
        ! Trafo 1
        query(1) = max(query(1),la*ao*lg*ao)
        query(2) = max(query(2),z*la*ao*lg)
        !Trafo 2
        query(3) = max(query(3),y*z*la*ao)
        !Trafo 3
        query(2) = max(query(2),x*y*z*la)
     else

        !do split_ao1 = 1, parts1

        !   fe1 = split_ao1*leng
        !   le1 = leng
        !   if(split_ao1 == parts1.and.mod(ao,leng)/=0) le1 = mod(ao,leng)

           !ZZ(ao,z)^T    WXYZ(la ao lg, ao)^T -> WRKWXYZ1 (z, la ao lg)
           call dgemm('t','t',z,la*ao*lg,ao,1.0E0_realk,ZZ,ao,WXYZ,la*ao*lg,0.0E0_realk,WRKWXYZ1,z)
           !YY(ao,y)^T WRKWXYZ1 (z la ao, lg)^T ->    WXYZ(y, z la ao)
           call dgemm('t','t',y,z*la*ao,lg,1.0E0_realk,YY(fg,1),ao,WRKWXYZ1,z*la*ao,0.0E0_realk,WRKWXYZ2,y)
           !XX(ao,x)^T    WXYZ(y z la, ao)^T   -> WRKWXYZ1 (x y z la)
           call dgemm('t','t',x,y*z*la,ao,1.0E0_realk,XX,ao,WRKWXYZ2,y*z*la,0.0E0_realk,WRKWXYZ1,x)
        !enddo

        !WW(ao,w)^T    WXYZ(x y z, la)^T   -+>  res(w x y z)
        call dgemm('t','t',w,x*y*z,la,1.0E0_realk,WW(fa,1),ao,WRKWXYZ1,x*y*z,1.0E0_realk,res,w)
     endif

  end subroutine successive_4ao_mo_trafo_exch

  subroutine init_query_info(query,n)
     implicit none
     type(pno_query_info) :: query
     integer, intent(in) :: n
     query%n_arrays = n
     call mem_alloc(query%size_array,n)
     query%size_array = 0
  end subroutine init_query_info

  subroutine free_query_info(query)
     implicit none
     type(pno_query_info) :: query
     query%n_arrays = 0
     call mem_dealloc(query%size_array)
  end subroutine free_query_info

#endif

 !subroutine not_implemented ccsd
     !!DEBUG: A2 term
     !!**************
     !u = p20*t2
     !call array_reorder_4d( m10, t2,   nv, no, nv, no, [1,4,3,2], p10, u  )
     !ref = gvovo

     !!A2.2 contribution
     !call array_reorder_4d( p10, gvvvv, nv, nv, nv, nv, [1,3,2,4], nul, w1  )
     !call array_reorder_4d( p10, t2,    nv, no, nv, no, [1,3,2,4], nul, w2  )
     !call dgemm( 'n', 'n', nv**2, no**2, nv**2, p10, w1, nv**2, w2, nv**2, nul, w3, nv**2 )
     !!call array_reorder_4d( p10, w3,    nv, nv, no, no, [1,3,2,4], p10, ref )

     !call print_norm(w3,o2v2,nnorm,.true.)
     !norm = 0.0E0_realk
     !!call print_norm(ref,o2v2,norm,.true.)
     !write(*,*)' DEBUG A2/TOT:',sqrt(nnorm),sqrt(norm)

     !!DEBUG: B2 term
     !!**************
     !!ref = 0.0E0_realk
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
     !!ref = 0.0E0_realk
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
     !ref = ref + w3(1:o2v2)


     !call print_norm(w3,o2v2,nnorm,.true.)
     !call print_norm(ref,o2v2,norm,.true.)
     !write (*,*)' DEBUG D2/TOT:',sqrt(nnorm),sqrt(norm)


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

     !ref = ref + w2(1:o2v2)

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

     !ref = ref + w2(1:o2v2)

     !call print_norm(w2,o2v2,nnorm,.true.)
     !call print_norm(ref,o2v2,norm,.true.)
     !write (*,*)' DEBUG E22/TOT:',sqrt(nnorm),sqrt(norm)


     !ref1 = vof

     !!DEBUG SINGLES A1
     !!****************
     !call array_reorder_4d( p10, u, nv, no, nv, no, [3,2,1,4], nul, w1) ! ckdi -> dkci
     !call dgemm( 'n','n',nv, no, nv**2*no, p10,gvvov,nv,w1,nv**2*no,nul,w2,nv)

     !ref1 = ref1 + w2(1:nv*no)

     !call print_norm(w2,i8*nv*no,nnorm,.true.)
     !call print_norm(ref1,i8*nv*no,norm,.true.)
     !write (*,*)' DEBUG A1/TOT:',sqrt(nnorm),sqrt(norm)


     !!DEBUG SINGLES B1
     !!****************
     !call array_reorder_4d( p10, gooov, no, no, no, nv, [1,4,3,2], nul, w1)
     !call dgemm('n','n',nv, no,nv*no**2,m10,u,nv,w1,nv*no**2,nul,w2,nv)

     !ref1 = ref1 + w2(1:nv*no)

     !call print_norm(w2,i8*nv*no,nnorm,.true.)
     !call print_norm(ref1,i8*nv*no,norm,.true.)
     !write (*,*)' DEBUG B1/TOT:',sqrt(nnorm),sqrt(norm)

     !!DEBUG SINGLES C1
     !!****************
     !call array_reorder_2d( p10, ovf, no, nv, [2,1], nul, w2)
     !call dgemv('n',no*nv,no*nv,p10,u,no*nv,w2,1, nul,w1,1)

     !ref1 = ref1 + w1(1:nv*no)

     !call print_norm(w1,i8*nv*no,nnorm,.true.)
     !call print_norm(ref1,i8*nv*no,norm,.true.)
     !write (*,*)' DEBUG C1/TOT:',sqrt(nnorm),sqrt(norm)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! ref is not written after this point!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !end subroutine not_implemented_ccsd
end module pno_ccsd_module
