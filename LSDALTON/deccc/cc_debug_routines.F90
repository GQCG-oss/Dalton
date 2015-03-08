!@file:
! This file should contain cc debug routines, noddy codes etc which we
! should be able to call with a keyword
!> \author Patrick Ettenhuber
!> \Date September 2013

module cc_debug_routines_module
   use precision
   use typedef
   use typedeftype
   use dec_typedef_module
   use tensor_interface_module
   use screen_mod
   use II_XC_interfaceModule
   use IntegralInterfaceMOD 

   ! DEC DEPENDENCIES (within deccc directory)   
   ! *****************************************
   use crop_tools_module
   use array2_simple_operations
   use array4_simple_operations
   use mp2_module
   use ccintegrals
   use ccsd_lhtr_module
   use ccsd_module
   use pno_ccsd_module
   use ccsdpt_module
   use orbital_operations
   use rpa_module


   contains
#ifdef MOD_UNRELEASED

   !> \author Marcin Ziolkowski (modified by Kasper Kristensen and Patrick
   !  Ettenhuber)
   subroutine ccsolver_debug(ccmodel,Co_f,Cv_f,fock_f,nbasis,nocc,nvirt, &
        & mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f,ccenergy, &
        & t1_final,t2_final,VOVO,longrange_singles,JOB,m2,m1,use_pnos,fraginfo)

     implicit none

     !> CC model
     integer,intent(in) :: ccmodel
     !> Number of occupied orbitals in full molecule/fragment AOS
     integer, intent(in) :: nocc
     !> Number of virtual orbitals in full molecule/fragment AOS
     integer, intent(in) :: nvirt
     !> Number of basis functions in full molecule/atomic extent
     integer, intent(in) :: nbasis
     !> Fock matrix in AO basis for fragment or full molecule
     real(realk), dimension(nbasis,nbasis), intent(in) :: fock_f
     !> Occupied MO coefficients for fragment/full molecule
     real(realk), dimension(nbasis,nocc), intent(in) :: Co_f
     !> Virtual MO coefficients for fragment/full molecule
     real(realk), dimension(nbasis,nvirt), intent(in) :: Cv_f
     !> Occ-occ block of Fock matrix in MO basis
     real(realk), dimension(nocc,nocc), intent(in) :: ppfock_f
     !> Virt-virt block of Fock matrix in MO basis
     real(realk), dimension(nvirt,nvirt), intent(in) :: qqfock_f
     real(realk),pointer :: dens(:,:)
     !> Is this a fragment job (true) or a full molecular calculation (false)
     logical, intent(in) :: fragment_job
     !> LS item information
     type(lsitem), intent(inout) :: mylsitem
     !> How much to print? ( ccPrintLevel>0 --> print info stuff)
     integer, intent(in) :: ccPrintLevel
     !> Coupled cluster energy for fragment/full molecule
     real(realk),intent(inout) :: ccenergy
     !> Final singles amplitudes
     type(array2),intent(inout) :: t1_final
     !> Final doubles amplitudes
     type(array4),intent(inout) :: t2_final
     !> Two electron integrals (a i | b j) stored as (a,i,b,j)
     type(array4),intent(inout) :: VOVO
     integer,intent(in) :: JOB
    
     type(array4),optional,intent(inout) :: m2
     type(array2),optional,intent(inout) :: m1
     logical,optional,intent(in)         :: use_pnos
     type(decfrag),optional,intent(inout)   :: fraginfo
     !> Include long-range singles effects using singles amplitudes
     !> from previous fragment calculations.
     !> IMPORTANT: If this it TRUE, then the singles amplitudes for the fragment
     !> (from previous calculations) must be stored in t1_final at input!
     logical,intent(in) :: longrange_singles
     real(realk),pointer :: Co2_d(:,:), Cv2_d(:,:),Co_d(:,:),Cv_d(:,:),focc(:),fvirt(:)
     real(realk),pointer :: ppfock_d(:,:),qqfock_d(:,:), Uocc(:,:), Uvirt(:,:)

     integer, dimension(2) :: occ_dims, virt_dims, ao2_dims, ampl2_dims
     integer, dimension(4) :: ampl4_dims
     type(array2) :: fock,Co,Cv,Co2,Cv2
     type(array2) :: ppfock,qqfock,pqfock,qpfock
     type(array4) :: gao,gmo,aibj,iajb
     type(array4), pointer :: t2(:),omega2(:)
     type(array2), pointer :: t1(:),omega1(:)
     type(array2) :: omega1_opt, t1_opt, omega1_prec
     type(array2) :: xocc,yocc,xvirt,yvirt,h1
     real(realk) :: two_norm_total, one_norm_total, one_norm1, one_norm2, &
          prev_norm
     real(realk), pointer :: B(:,:),c(:)
     integer :: iter,last_iter,i,j,k,l
     logical :: crop_ok,break_iterations,get_mult
     type(array4) :: omega2_opt, t2_opt, omega2_prec, u
     type(array2) :: ifock,delta_fock,fockguess
     type(ri) :: l_ao
     type(array2) :: ppfock_prec, qqfock_prec,t1tmp
     type(array4) :: Lmo
     real(realk) :: tcpu, twall, ttotend_cpu, ttotend_wall, ttotstart_cpu, ttotstart_wall
     real(realk) :: iter_cpu,iter_wall, sosex
     real(realk) :: t1fnorm2,m1fnorm2,t2fnorm2,m2fnorm2
     character(18) :: save_to,keep
     character(tensor_MSG_LEN) :: msg
     integer :: ii,aa, cc
     integer :: MaxSubSpace
     logical :: restart, u_pnos
     type(tensor) :: govov
     integer :: nspaces
     type(PNOSpaceInfo), pointer :: pno_cv(:), pno_S(:)


     call LSTIMER('START',ttotstart_cpu,ttotstart_wall,DECinfo%output)
     if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)
   
     u_pnos   = .false.
     if(present(use_pnos)) u_pnos = use_pnos
     get_mult = (JOB==SOLVE_MULTIPLIERS)

     ! Sanity check 1: Number of orbitals
     if( (nvirt < 1) .or. (nocc < 1) ) then
        write(DECinfo%output,*) 'Number of occupied orbitals = ', nocc
        write(DECinfo%output,*) 'Number of virtual  orbitals = ', nvirt
        call lsquit('ccsolver: Empty occupied or virtual space!',DECinfo%output)
     endif

     ! Sanity check 2: Singles amplitudes initiated appropriately
     if(longrange_singles) then
        if(.not. associated(t1_final%val)) then
           call lsquit('ccsolver: Long range singles corrections requested, &
                & but t1_final does not contain existing amplitudes!',DECinfo%output)
        end if
     end if

     if(get_mult)then

       ! Sanity check 3: if CCSD multipliers are requested, make sure that both are there
       if((present(m2).and..not.present(m1)).or.(present(m1).and..not.present(m2)))then
         call lsquit("ERROR(ccsolver_debug):requested unkown multipliers",-1)
       endif
     endif

     if(u_pnos)then
       ! Sanity check 4: if PNO use is requested, currently only ccsd is implemented, no
       !multipliers, though 
       if(.not.present(m2))then
         call lsquit("ERROR(ccsolver_debug):PNO ccsd requires mp2 amplitudes passed as m2",-1)
       endif

       if( .not.present(fraginfo).and. fragment_job )then
         call lsquit("ERROR(ccsolver_debug):PNO ccsd requires the fragment information if it is a fragment job",-1)
       endif

       !if( .not. associated(VOVO%val)) then
       !  call lsquit("ERROR(ccsolver_debug):PNO ccsd requires the VOVO integrals",-1)
       !endif
     endif



     MaxSubSpace = DECinfo%ccMaxDIIS

     ! title
     Call print_ccjob_header(ccmodel,ccPrintLevel,fragment_job,get_mult,nbasis,nocc,nvirt,MaxSubSpace,.false.,.false.,1)

     ! dimension vectors
     occ_dims   = [nbasis,nocc]
     virt_dims  = [nbasis,nvirt]
     ao2_dims   = [nbasis,nbasis]
     ampl4_dims = [nvirt,nocc,nvirt,nocc]
     ampl2_dims = [nvirt,nocc]

     ! go to a (pseudo) canonical basis
     call mem_alloc(focc,nocc)
     call mem_alloc(fvirt,nvirt)
     call mem_alloc(Co_d,nbasis,nocc)
     call mem_alloc(Cv_d,nbasis,nvirt)
     call mem_alloc(Co2_d,nbasis,nocc)
     call mem_alloc(Cv2_d,nbasis,nvirt)
     call mem_alloc(ppfock_d,nocc,nocc)
     call mem_alloc(qqfock_d,nvirt,nvirt)
     call mem_alloc(Uocc,nocc,nocc)
     call mem_alloc(Uvirt,nvirt,nvirt)


     !TRANSFORM TO A CANONICAL BASIS
     !******************************

     ! prevent if explicitly requested or if PNOs are requested

     if(DECinfo%CCSDpreventcanonical)then
       !nocc diagonalization
       Co_d    = Co_f
       Cv_d    = Cv_f
       ppfock_d = ppfock_f
       qqfock_d = qqfock_f
       Uocc     = 0.0E0_realk
       Uvirt    = 0.0E0_realk
       do ii=1,nocc
         Uocc(ii,ii) = 1.0E0_realk
       enddo

       do aa=1,nvirt
         Uvirt(aa,aa) = 1.0E0_realk
       enddo

     else

       call get_canonical_integral_transformation_matrices(nocc,nvirt,nbasis,ppfock_f,&
            &qqfock_f,Co_f,Cv_f,Co_d,Cv_d,Uocc,Uvirt,focc,fvirt)

       if(u_pnos)then
         !Do not destroy the locality of the occupied space if PNOs are used,
         !otherwise the adaption cannot happen to a confined space
         Co_d     = Co_f
         ppfock_d = ppfock_f
         do ii=1,nocc
           Uocc(ii,ii) = 1.0E0_realk
         enddo

       else

         ppfock_d = 0.0E0_realk
         do ii=1,nocc
           ppfock_d(ii,ii) = focc(ii)
         enddo

       endif

       qqfock_d = 0.0E0_realk
       do aa=1,nvirt
         qqfock_d(aa,aa) = fvirt(aa)
       enddo

       if(get_mult)then

         if(DECinfo%use_singles)then
           call local_can_trans(nocc,nvirt,nbasis,Uocc,Uvirt,vovo=t2_final%val,vo=t1_final%val)
         else
           call local_can_trans(nocc,nvirt,nbasis,Uocc,Uvirt,vovo=t2_final%val)
         endif

       elseif(u_pnos)then
         call local_can_trans(nocc,nvirt,nbasis,Uocc,Uvirt,vovo=m2%val)
       endif


     endif

     call mem_dealloc(focc)
     call mem_dealloc(fvirt)

     ! Copy MO coeffcients. It is very convenient to store them twice to handle transformation
     ! (including transposed MO matrices) efficiently. 
     Co2_d = Co_d
     Cv2_d = Cv_d


     ! create transformation matrices in array form
     Co    = array2_init(occ_dims,Co_d)
     Cv    = array2_init(virt_dims,Cv_d)
     Co2   = array2_init(occ_dims,Co2_d)
     Cv2   = array2_init(virt_dims,Cv2_d)
     fock  = array2_init(ao2_dims,fock_f)

     call mem_dealloc(Co_d)
     call mem_dealloc(Cv_d)
     call mem_dealloc(Co2_d)
     call mem_dealloc(Cv2_d)


     ! Get Fock matrix correction (for fragment and/or frozen core)
     ! ************************************************************
     ! Full molecule/frozen core: The correction corresponds to difference between actual Fock matrix
     !                            and Fock matrix where the density is made from only valence orbitals.
     ! Fragment: The correction correspond to the difference between actual Fock matrix
     !           and Fock matrix calculated from a "fragment density" determined from
     !           fragment's occupied molecular orbitals (which for frozen core includes only valence
     !           orbitals).

     ! Density corresponding to input MOs
     !call mem_alloc(dens,nbasis,nbasis)
     !call get_density_from_occ_orbitals(nbasis,nocc,Co%val,dens)
       
     !call mem_dealloc(dens)

     ! get two-electron integrals in ao
     if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a)') 'debug :: calculating AO integrals'
     ! Always in the DEBUG SOLVER: calculate full 4-dimensional AO integrals 
     call get_full_eri(mylsitem,nbasis,gao)
     if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,/)') 'debug :: AO integrals done'

     if(DECinfo%PL>1) call LSTIMER('CCSOL: INIT',tcpu,twall,DECinfo%output)
     if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)


     ! special MP2 things
     MP2Special : if(CCmodel == MODEL_MP2 .or. CCmodel == MODEL_RPA) then

        write(DECinfo%output,*)
        write(DECinfo%output,*) ' ********************  WARNING  **********************'
        write(DECinfo%output,*) 'CCsolver is called for MP2 model.'
        write(DECinfo%output,*) 'This will work fine but it is recommended to use the non-iterative'
        write(DECinfo%output,*) 'MP2_integrals_and_amplitudes_workhorse to use get the MP2 amplitudes'
        write(DECinfo%output,*)
        call get_VOVO_integrals(mylsitem,nbasis,nocc,nvirt,Cv%val,Co%val,gmo)

        ! Construct L: L_{bjai} = 2*g_{bjai} - g_{ajbi}
        Lmo = getL_simple_from_gmo(gmo)

        ppfock = array2_similarity_transformation(Co,fock,Co2,[nocc,nocc])
        qqfock = array2_similarity_transformation(Cv,fock,Cv2,[nvirt,nvirt])

        if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') ' debug :: gmo(vovo) norm  : ',gmo*gmo
     end if MP2Special

     ! get fock matrices for preconditioning
     Preconditioner : if(DECinfo%use_preconditioner .or. DECinfo%use_preconditioner_in_b) then
        if(DECinfo%precondition_with_full) then
           ppfock_prec = array2_init([nocc,nocc],ppfock_d)
           qqfock_prec = array2_init([nvirt,nvirt],qqfock_d)
        else
           ppfock_prec = array2_similarity_transformation(Co,fock,Co2,[nocc,nocc])
           qqfock_prec = array2_similarity_transformation(Cv,fock,Cv2,[nvirt,nvirt])
        end if
     end if Preconditioner
     call mem_dealloc(ppfock_d)
     call mem_dealloc(qqfock_d)


     ! allocate things
     if(DECinfo%use_singles) then
        call mem_alloc(t1,DECinfo%ccMaxIter)
        call mem_alloc(omega1,DECinfo%ccMaxIter)
     end if
     call mem_alloc(t2,DECinfo%ccMaxIter)
     call mem_alloc(omega2,DECinfo%ccMaxIter)

     ! initialize T1 matrices and fock transformed matrices for CC pp,pq,qp,qq
     if(decinfo%ccmodel /= MODEL_MP2 .and. decinfo%ccmodel /= MODEL_RPA ) then
        xocc = array2_init(occ_dims)
        yocc = array2_init(occ_dims)
        xvirt = array2_init(virt_dims)
        yvirt = array2_init(virt_dims)
        h1 = array2_init_plain(ao2_dims)
        CALL II_get_h1_mixed_full(DECinfo%output,DECinfo%output,MyLsitem%SETTING,&
             & h1%val,nbasis,nbasis,AORdefault,AORdefault)
     end if

     call mem_alloc(B,DECinfo%ccMaxIter,DECinfo%ccMaxIter)
     call mem_alloc(c,DECinfo%ccMaxIter)


     ! Get MO int. (govov) for RPA:
     if (ccmodel==MODEL_RPA) then
       call tensor_minit(govov, [nocc,nvirt,nocc,nvirt],4,local=.true.,atype='TDAR')
       call tensor_zero(govov)
       ! FIXME: Johannes if this code is still suppose to work, you should
       !        find govov in a different way, else please clean this file.
       ! Date:  Nov. 2014
       call lsquit("RPA in debug solver is missing govov",DECinfo%output)
     end if




     set_pno_info:if(u_pnos)then



       !GET THE PNO TRANSFORMATION MATRICES
       if(fragment_job)then

                   !COUNT PAIRS OUTSIDE EOS
         nspaces = ( nocc - fraginfo%noccEOS ) * ( nocc - fraginfo%noccEOS + 1) / 2 &
                   !COUNT PAIRS WITH 1 IDX IN EOS                   !EOS
                &+ fraginfo%noccEOS * ( nocc - fraginfo%noccEOS ) + 1

         fraginfo%nspaces = nspaces

         call mem_alloc( fraginfo%CLocPNO, nspaces )
         call get_pno_trafo_matrices(nocc,nvirt,nbasis,m2%val,&
         &fraginfo%CLocPNO,fraginfo%nspaces,f=fraginfo)
         pno_cv => fraginfo%CLocPNO

       else
                   !ALL PAIRS
         nspaces = nocc * ( nocc + 1 ) / 2
         call mem_alloc( pno_cv, nspaces )
         call get_pno_trafo_matrices(nocc,nvirt,nbasis,m2%val,&
         &pno_cv,nspaces,f=fraginfo)

       endif


       !do ii = 1,nspaces
       !   if(pno_cv(ii)%allocd)call print_norm(pno_cv(ii)%d,size(pno_cv(ii)%d,kind=8))
       !enddo

       !GET THE OVERLAP BETWEEN THE PNO SPACES
       call mem_alloc( pno_S , nspaces * (nspaces - 1)/2 )   
       !Get all the overlap matrices necessary
       call get_pno_overlap_matrices(nocc,nvirt,pno_cv,pno_S,nspaces,.true.)

       !print *,"-----------------------------------"
       !do ii = 1,nspaces * (nspaces - 1)/2 
       !   if(pno_S(ii)%allocd)call print_norm(pno_S(ii)%d,size(pno_S(ii)%d,kind=8))
       !enddo

     endif set_pno_info




     ! readme : the iteration sequence is universal and may be used for all
     !          iterative cc models (linear or non-linear) and is
     !          semi-independent on the storage of vectors (allocation and
     !          deallocation, etc)

     ! iterate
     break_iterations = .false.
     crop_ok          = .false.
     prev_norm        = 1.0E6_realk


     CCIteration : do iter=1,DECinfo%ccMaxIter

        if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)
        call LSTIMER('START',iter_cpu,iter_wall,DECinfo%output)

        ! remove old vectors
        RemoveOldVectors : if(iter > MaxSubSpace) then

           if(DECinfo%cc_driver_debug) then
              write(DECinfo%output,'(a,i4)') ' debug :: vector to delete : ',iter-MaxSubSpace
           end if

           if(DECinfo%use_singles) then
              call array2_free(t1(iter-MaxSubSpace))
              Call array2_free(omega1(iter-MaxSubSpace))
           end if
           call array4_free(t2(iter-MaxSubSpace))
           call array4_free(omega2(iter-MaxSubSpace))
        end if RemoveOldVectors


        ! If we want the multipliers, the t1 transformation has to be calculated
        ! with the t1_final and not the t1(iter) (so actually only once, but as
        ! it is a debug code this solution was simpler)
        if(iter == 1.and.DECinfo%use_singles .and. get_mult) then
           call getT1transformation(t1_final,xocc,xvirt,yocc,yvirt, &
           Co,Cv,Co2,Cv2)
        endif


        ! get new amplitude vectors
        GetGuessVectors : if(iter == 1) then
   
           if(DECinfo%use_singles)then
             call get_guess_vectors_simple(mylsitem,&
             &gao,get_mult,Co,Co2,Cv2,nocc,nvirt,nbasis,&
             &xocc,yvirt,restart,t2=t2(iter),t1=t1(iter),&
             &t1f=t1_final)
           else
             call get_guess_vectors_simple(mylsitem,&
             &gao,get_mult,Co,Co2,Cv2,nocc,nvirt,nbasis,&
             &xocc,yvirt,restart,t2=t2(iter))
           endif

        end if GetGuessVectors

        ! Initialize residual vectors
        if(DECinfo%use_singles)then
          omega1(iter) = array2_init(ampl2_dims)
        endif
        if(DECinfo%array4OnFile) then
           ! KK, initialize omega2(iter) using storing type 2
           omega2(iter) = array4_init(ampl4_dims,2,.true.)
        else
           omega2(iter) = array4_init(ampl4_dims)
        endif


        ! GET SINGLES:
        T1Related : if(DECinfo%use_singles) then

           ! get the T1 transformation matrices
           if(.not.get_mult)then
              call getT1transformation(t1(iter),xocc,xvirt,yocc,yvirt,Co,Cv,Co2,Cv2)
              ifock = getInactiveFock_simple(h1,gao,xocc,yocc,nocc,nbasis)
           else
              ifock = getInactiveFock_simple(h1,gao,Co,Co,nocc,nbasis)
           endif

           ! Note: If not fock_with_ri or ccsd_old, then the relevant
           ! ifock is calculated below in get_ccsd_residual_integral_direct.
           ! Long range fock matrix correction using old scheme
           ! (See comments above regarding Fock correction)
           if(iter == 1) then
              ! calculate fock correction in first iteration
              write(DECinfo%output,'(a)') 'long range fock correction requested'
              if(fragment_job) then
                 delta_fock = getFockCorrection(fock,ifock)
              else ! full molecule: correction is zero by definition
                 delta_fock= array2_init(ao2_dims)
              end if

           end if

           ! Add fock correction to to existing T1-transformed Fock matrix
           call array2_add_to(ifock,1.0E0_realk,delta_fock)

           ! readme : this should be done in a more clear way
           if(CCmodel == MODEL_CC2) then
              ! CC2
              ppfock = array2_similarity_transformation(xocc,fock,yocc,[nocc,nocc])
              qqfock = array2_similarity_transformation(xvirt,fock,yvirt,[nvirt,nvirt])
           else if(CCmodel >= MODEL_CCSD) then
              ! CCSD
              ppfock = array2_similarity_transformation(xocc,ifock,yocc,[nocc,nocc])
              qqfock = array2_similarity_transformation(xvirt,ifock,yvirt,[nvirt,nvirt])
           endif

           pqfock = array2_similarity_transformation(xocc,ifock,yvirt,[nocc,nvirt])
           qpfock = array2_similarity_transformation(xvirt,ifock,yocc,[nvirt,nocc])

           iajb = get_gmo_simple(gao,xocc,yvirt,xocc,yvirt)

        end if T1Related

        if(DECinfo%PL>1) call LSTIMER('CCIT: INIT',tcpu,twall,DECinfo%output)
        if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)



        ! MODIFY FOR NEW MODEL
        ! If you implement a new model, please insert call to your own residual routine here!
        SelectCoupledClusterModel : if(CCmodel==MODEL_MP2) then

          call getDoublesResidualMP2_simple(Omega2(iter),t2(iter),gmo,ppfock,qqfock, &
                & nocc,nvirt)

        elseif(CCmodel==MODEL_RIMP2) then

           call lsquit('MODEL_RIMP2 have no residual non iterative scheme',-1)

        elseif(CCmodel==MODEL_CC2) then
           u = get_u(t2(iter))
           call getSinglesResidualCCSD(omega1(iter),u,gao,pqfock,qpfock, &
                xocc,xvirt,yocc,yvirt,nocc,nvirt)
           call array4_free(u)

           gmo = get_gmo_simple(gao,xvirt,yocc,xvirt,yocc)
           call getDoublesResidualMP2_simple(omega2(iter),t2(iter),gmo,ppfock,qqfock, &
                nocc,nvirt)
           call array4_free(gmo)


        elseif(CCmodel==MODEL_CCSD .or. CCmodel==MODEL_CCSDpT) then  ! CCSD or CCSD(T)
  

           if(get_mult)then
              if(u_pnos)then
                call lsquit("ERROR(cc_driver_debug):only one of get_mult,&
                &and use_pnos should be true at the same time",-1)
              endif

              call get_ccsd_multipliers_simple(omega1(iter)%val,omega2(iter)%val,t1_final%val&
              &,t2_final%val,t1(iter)%val,t2(iter)%val,delta_fock%val,xocc%val,yocc%val,xvirt%val,yvirt%val&
              &,nocc,nvirt,nbasis,MyLsItem,gao_ex=gao)
           
           else if(u_pnos.and..not.DECinfo%hack2)then
              if(get_mult)then
                call lsquit("ERROR(cc_driver_debug):only one of get_mult,&
                &and use_pnos should be true at the same time",-1)
              endif

              !transform back to original basis   
              !if(.not.DECinfo%CCSDpreventcanonical)then
              !  if(DECinfo%use_singles)then
              !    call ccsolver_can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
              !    &vovo=t2(iter)%val,vo=t1(iter)%val,bo=xocc%val,bv=xvirt%val)
              !    call ccsolver_can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
              !    &vovo=omega2(iter)%val,vo=omega1(iter)%val,bo=yocc%val,bv=yvirt%val)
              !  else
              !    call ccsolver_can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
              !    &vovo=t2(iter)%val,bo=xocc%val,bv=xvirt%val)
              !    call ccsolver_can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
              !    &vovo=omega2(iter)%val,bo=yocc%val,bv=yvirt%val)
              !  endif
              !endif
     
              !FIXME: fock is defined in a wrong way !!!
              !if(iter == 1) t2(iter)%val = m2%val
              !if(.not.fragment_job)then
              !  call get_ccsd_residual_pno_style(t1(iter)%val,t2(iter)%val,omega1(iter)%val,&
              !  &omega2(iter)%val,iajb%val,nocc,nvirt,nbasis,xocc%val,xvirt%val,yocc%val,yvirt%val,mylsitem,&
              !  &fragment_job,pno_cv,pno_S,nspaces,fock%val,iter)
              !else
              !  call get_ccsd_residual_pno_style(t1(iter)%val,t2(iter)%val,omega1(iter)%val,&
              !  &omega2(iter)%val,iajb%val,nocc,nvirt,nbasis,xocc%val,xvirt%val,yocc%val,yvirt%val,mylsitem,&
              !  &fragment_job,pno_cv,pno_S,nspaces,fock%val,iter,f=fraginfo)
              !endif

              !stop 0
         
              !transform to pseudo diagonal basis for the solver
              !if(.not.DECinfo%CCSDpreventcanonical)then
              !  if(DECinfo%use_singles)then
              !    call ccsolver_local_can_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
              !    &vovo=t2(iter)%val,vo=t1(iter)%val,bo=xocc%val,bv=xvirt%val)
              !    call ccsolver_local_can_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
              !    &vovo=omega2(iter)%val,vo=omega1(iter)%val,bo=yocc%val,bv=yvirt%val)
              !  else
              !    call ccsolver_local_can_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
              !    &vovo=t2(iter)%val,bo=xocc%val,bv=xvirt%val)
              !    call ccsolver_local_can_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
              !    &vovo=omega2(iter)%val,bo=yocc%val,bv=yvirt%val)
              !  endif
              !endif

           else

             u = get_u(t2(iter))
             call getSinglesResidualCCSD(omega1(iter),u,gao,pqfock,qpfock,xocc,xvirt,yocc,yvirt,nocc,nvirt)
             aibj = get_gmo_simple(gao,yocc,xvirt,yocc,xvirt)
             call array4_reorder(aibj,[2,1,4,3])
             call getDoublesResidualCCSD_simple(omega2(iter),t2(iter),u,gao,aibj,iajb,nocc,nvirt, &
                  ppfock,qqfock,xocc,xvirt,yocc,yvirt)

             call array4_free(aibj)
             call array4_free(u)
           endif


        elseif(CCmodel==MODEL_RPA) then

           !if(get_mult)then
           !  !rpa_multipliers not yet implemented
           !  call RPA_multiplier(Omega2(iter),t2_final,t2(iter),gmo,ppfock,qqfock,nocc,nvirt)
           !else
!#ifdef VAR_MPI
!             call RPA_residualpar(Omega2(iter),t2(iter),govov%elm1,ppfock,qqfock,nocc,nvirt)
!#else
             !msg =' Norm of gmo'
             !call print_norm(govov,msg)
             call RPA_residualdeb(Omega2(iter),t2(iter),govov%elm1,ppfock,qqfock,nocc,nvirt)
!#endif
           !call lsquit('ccsolver_debug: Residual for model is not implemented!',-1)
           !  call RPA_residual(Omega2(iter),t2(iter),govov,ppfock,qqfock,nocc,nvirt)
           !endif


        else
           print *, 'MODEL = ', DECinfo%cc_models(DECinfo%ccmodel)
           call lsquit('ccsolver_debug: Residual for model is not implemented!',-1)

        end if SelectCoupledClusterModel

        if(DECinfo%PL>1) call LSTIMER('CCIT: RESIDUAL',tcpu,twall,DECinfo%output)


        ForDebug : if(DECinfo%cc_driver_debug) then

           if(DECinfo%use_singles) then
              write(DECinfo%output,'(a,f16.10)') ' debug :: t1 norm         ',t1(iter)*t1(iter)
              write(DECinfo%output,'(a,f16.10)') ' debug :: omega1 norm     ',omega1(iter)*omega1(iter)
           end if

           write(DECinfo%output,'(a,f16.10)') ' debug :: t2 norm         ',t2(iter)*t2(iter)
           write(DECinfo%output,'(a,f16.10)') ' debug :: omega2 norm     ',omega2(iter)*omega2(iter)
           write(DECinfo%output,'(a,f16.10)') ' debug :: ppfock norm     ',ppfock*ppfock
           write(DECinfo%output,'(a,f16.10)') ' debug :: qqfock norm     ',qqfock*qqfock

           if(DECinfo%use_singles) then
              write(DECinfo%output,'(a,f16.10)') ' debug :: pqfock norm     ',pqfock*pqfock
              write(DECinfo%output,'(a,f16.10)') ' debug :: qpfock norm     ',qpfock*qpfock
           end if

        end if ForDebug

        if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)

        ! calculate crop/diis matrix
        B=0.0E0_realk; c=0.0E0_realk
        do i=iter,max(iter-MaxSubSpace+1,1),-1
           do j=iter,i,-1
              if(DECinfo%use_singles) then
                 if(DECinfo%use_preconditioner_in_b) then
                    omega1_prec = precondition_singles(omega1(j),ppfock_prec,qqfock_prec)
                    omega2_prec = precondition_doubles(omega2(j),ppfock_prec,qqfock_prec)
                    B(i,j) = omega1(i)*omega1_prec + omega2(i)*omega2_prec
                    call array2_free(omega1_prec)
                    call array4_free(omega2_prec)
                 else
                    B(i,j) = omega1(i)*omega1(j) + omega2(i)*omega2(j)
                 end if
              else
                 ! just doubles
                 if(DECinfo%use_preconditioner_in_b) then
                    omega2_prec = precondition_doubles(omega2(j),ppfock_prec,qqfock_prec)
                    B(i,j) = omega2(i)*omega2_prec
                    call array4_free(omega2_prec)
                 else
                    B(i,j) = omega2(i)*omega2(j)
                    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,i4,a,i4,a,f16.10)') &
                         ' debug :: B(',i,',',j,')=',B(i,j)
                 end if
              end if
              B(j,i) = B(i,j)
           end do
        end do

        if(DECinfo%PL>1) call LSTIMER('CCIT: CROP MAT',tcpu,twall,DECinfo%output)

        ! solve crop/diis equation
        call CalculateDIIScoefficients(MaxSubSpace,DECinfo%ccMaxIter,iter,B,c, &
             DECinfo%PL>3)

        if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)

        ! mixing to get optimal
        if(DECinfo%use_singles) then
           t1_opt = array2_init(ampl2_dims)
           omega1_opt = array2_init(ampl2_dims)
        end if
        if(DECinfo%array4OnFile) then ! store array elements of file (storing type 2)
           omega2_opt = array4_init(ampl4_dims,2,.true.)
           t2_opt = array4_init(ampl4_dims,2,.true.)
        else
           omega2_opt = array4_init(ampl4_dims)
           t2_opt = array4_init(ampl4_dims)
        end if

        do i=iter,max(iter-MaxSubSpace+1,1),-1
           ! mix singles
           if(DECinfo%use_singles) then
              call array2_add_to(t1_opt,c(i),t1(i))
              call array2_add_to(omega1_opt,c(i),omega1(i))
           end if
           ! mix doubles
           call array4_add_to(t2_opt,c(i),t2(i))
           call array4_add_to(omega2_opt,c(i),omega2(i))
        end do

        if(DECinfo%PL>1) call LSTIMER('CCIT: MIXING',tcpu,twall,DECinfo%output)
        if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)

        ! if crop, put the optimal in place of trial (not for diis)
        if(DECinfo%use_crop) then
           if(DECinfo%use_singles) then
              call array2_free(t1(iter))
              call array2_free(omega1(iter))
              t1(iter) = array2_duplicate(t1_opt)
              omega1(iter) = array2_duplicate(omega1_opt)
           end if
           call array4_free(t2(iter))
           call array4_free(omega2(iter))
           if(DECinfo%array4OnFile) then
              t2(iter) = array4_duplicate_same_file(t2_opt)
              omega2(iter) = array4_duplicate_same_file(omega2_opt)
           else
              t2(iter) = array4_duplicate(t2_opt)
              omega2(iter) = array4_duplicate(omega2_opt)
           end if
        end if

        if(DECinfo%PL>1) call LSTIMER('CCIT: COPY OPT',tcpu,twall,DECinfo%output)
        if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)

        !HACK to get CCD
        if(DECinfo%CCDHACK)then
          omega1(iter)%val = 0.0E0_realk
        endif

        ! check for the convergence
        one_norm1 = 0.0E0_realk
        one_norm2 = 0.0E0_realk
        if(DECinfo%use_singles) one_norm1 = array2_norm(omega1(iter))
        one_norm2 = array4_norm(omega2(iter))
        one_norm_total = one_norm1 + one_norm2
        two_norm_total = sqrt(one_norm_total)

        ! simple crop diagnostics
        if(two_norm_total < prev_norm) then
           crop_ok=.true.
        else
           crop_ok=.false.
           write(DECinfo%output,'(a)') ' warning :: total norm was smaller in previous iteration !!! '
        end if
        prev_norm=two_norm_total

        if(DECinfo%PL>1) call LSTIMER('CCIT: CONV',tcpu,twall,DECinfo%output)
        if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)

        ! calculate the correlation energy and fragment energy
        ! MODIFY FOR NEW MODEL
        ! If you implement a new model, please insert call to energy routine here,
        ! or insert a call to get_cc_energy if your model uses the standard CC energy expression.
        if(.not. get_mult)then
           EnergyForCCmodel: if(CCmodel==MODEL_MP2.OR.CCmodel==MODEL_RIMP2) then  
              ! MP2
              ccenergy = get_mp2_energy(t2(iter),Lmo)
           elseif(CCmodel==MODEL_CC2 .or. CCmodel==MODEL_CCSD .or. CCmodel==MODEL_CCSDpT )then
              ! CC2, CCSD, or CCSD(T) (for (T) calculate CCSD contribution here)
              ccenergy = get_cc_energy(t1(iter),t2(iter),iajb,nocc,nvirt)
           elseif(CCmodel==MODEL_RPA) then
             !Here the energy is computed not in P U [bar P]
             !but in the AOS
              !ccenergy = RPA_energy(t2(iter),govov)
              !sosex = SOSEX_contribution(t2(iter),govov)
              !ccenergy=ccenergy+sosex
           else
              print *, 'MODEL = ', DECinfo%cc_models(DECinfo%ccmodel)
              call lsquit('ccsolver_debug: Energy for model is not implemented!',-1)
           end if EnergyForCCmodel
        endif


        if(DECinfo%PL>1) call LSTIMER('CCIT: ENERGY',tcpu,twall,DECinfo%output)

        ! check if this is the last iteration
        if(iter == DECinfo%ccMaxIter .or. two_norm_total < DECinfo%ccConvergenceThreshold) &
             break_iterations=.true.

        if(DECinfo%use_singles .and. (.not. break_iterations)) then
          call array4_free(iajb)
        end if

        if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)

        ! generate next trial vector if this is not the last iteration
        if(.not.break_iterations) then
           if(DECinfo%use_preconditioner) then
              if(DECinfo%use_singles) then
                 omega1_prec = precondition_singles(omega1_opt,ppfock_prec,qqfock_prec)
                 t1(iter+1) = t1_opt + omega1_prec
                 call array2_free(omega1_prec)
              end if
              omega2_prec = precondition_doubles(omega2_opt,ppfock_prec,qqfock_prec)
              t2(iter+1) = t2_opt + omega2_prec
              call array4_free(omega2_prec)
           else
              if(DECinfo%use_singles) t1(iter+1) = t1_opt + omega1_opt
              t2(iter+1) = t2_opt + omega2_opt
           end if
           
           if(.not.DECinfo%CCSDnosaferun)then
             if(DECinfo%use_singles)then
               call save_current_guess_simple(iter,get_mult,t2=t2(iter),t1=t1(iter))
             else
               call save_current_guess_simple(iter,get_mult,t2=t2(iter))
             endif
           endif
        end if

        if(DECinfo%PL>1) call LSTIMER('CCIT: NEXT VEC',tcpu,twall,DECinfo%output)

        ! delete optimals
        if(DECinfo%use_singles) then
           call array2_free(t1_opt)
           call array2_free(omega1_opt)
        end if
        if(DECinfo%array4OnFile) then
           ! Free optimial arrays BUT keep file, because the same file is used by t2(iter)
           call array4_free(t2_opt,keep=.true.)
           call array4_free(omega2_opt,keep=.true.)
        else
           call array4_free(t2_opt)
           call array4_free(omega2_opt)
        end if


        ! delete fock matrices
        if(DECinfo%use_singles) then
           call array2_free(ifock)
           call array2_free(ppfock)
           call array2_free(pqfock)
           call array2_free(qpfock)
           call array2_free(qqfock)
        end if


        call LSTIMER('CC ITERATION',iter_cpu,iter_wall,DECinfo%output)

#ifdef __GNUC__
        call flush(DECinfo%output)
#endif
    

        !Print Iter info
        !---------------
        call print_ccjob_iterinfo(iter,two_norm_total,ccenergy,get_mult,fragment_job)

        last_iter = iter
        if(break_iterations) exit

     end do CCIteration

     call LSTIMER('START',ttotend_cpu,ttotend_wall,DECinfo%output)



     ! Free memory and save final amplitudes
     ! *************************************

     ! free memory
     if (ccmodel==MODEL_RPA) call tensor_free(govov)

     ! remove rest of the singles amplitudes and residuals
     do i=last_iter,max(last_iter-MaxSubSpace+1,1),-1


        ! remove the lase files of t2 and omega2
        call array4_delete_file(omega2(i))
        call array4_delete_file(t2(i))

        if(DECinfo%use_singles) then

           ! Save final singles amplitudes
           if(i==last_iter) then
              if(get_mult)then
                   m1 = array2_duplicate(t1(last_iter))
              else
                if(longrange_singles) then ! just copy
                   call array2_copy(t1_final,t1(last_iter))
                else ! initialize and copy
                   t1_final = array2_duplicate(t1(last_iter))
                end if
              endif
           end if

           ! Free singles amplitudes and residuals
           call array2_free(t1(i))
           call array2_free(omega1(i))

        end if

        ! Free doubles residuals
        call array4_free(omega2(i))

        ! Save final double amplitudes
        if(i==last_iter) then
           if(get_mult)then
             m2       = array4_duplicate(t2(last_iter))
           else
             t2_final = array4_duplicate(t2(last_iter))
           endif
        end if

        ! Free doubles amplitudes
        call array4_free(t2(i))

     end do
     !if(.not.DECinfo%use_singles)then
     !  t1_final = array2_init([nvirt,nocc])
     !endif
     call print_norm(t2_final,t2fnorm2,.true.)
     if(DECinfo%use_singles)then
       call print_norm(t1_final,t1fnorm2,.true.)
     endif
     if(get_mult)then
        call print_norm(m2,m2fnorm2,.true.)
        if(DECinfo%use_singles)then
           call print_norm(m1,m1fnorm2,.true.)
        endif
     endif


     ! Write finalization message
     !---------------------------
     call print_ccjob_summary(break_iterations,get_mult,fragment_job,last_iter,DECinfo%use_singles, &
     &ccenergy,ttotend_wall,ttotstart_wall,ttotend_cpu,ttotstart_cpu,t1fnorm2,t2fnorm2, nm1 = m1fnorm2, nm2 = m2fnorm2)


     ! Save two-electron integrals in the order (virt,occ,virt,occ)
     if(CCmodel == MODEL_MP2 .or. CCmodel==MODEL_RPA) then
        call array4_free(lmo) ! also free lmo integrals
        VOVO = array4_duplicate(gmo)
        call array4_free(gmo)
     else
        VOVO = array4_duplicate(iajb)
        call array4_free(iajb)
        call array4_reorder(VOVO,[2,1,4,3])
     end if

     ! deallocate stuff
     if(DECinfo%use_singles) then
        call mem_dealloc(t1)
        call mem_dealloc(omega1)
     end if

     call mem_dealloc(t2)
     call mem_dealloc(omega2)

     call mem_dealloc(B)
     call mem_dealloc(c)


     ! remove fock correction
     call array2_free(delta_fock)
     call array4_free(gao)

     if(DECinfo%use_preconditioner .or. DECinfo%use_preconditioner_in_b) then
        call array2_free(ppfock_prec)
        call array2_free(qqfock_prec)
     end if

     if(DECinfo%use_singles) then
        call array2_free(h1)
        call array2_free(xocc)
        call array2_free(yocc)
        call array2_free(xvirt)
        call array2_free(yvirt)
        call array2_free(pqfock)
        call array2_free(qpfock)
     end if

     call array2_free(ppfock)
     call array2_free(qqfock)

     call array2_free(Co)
     call array2_free(Co2)
     call array2_free(Cv)
     call array2_free(Cv2)
     call array2_free(fock)


     !Free PNO information
     if(u_pnos)then

       if(.not.fragment_job)then
         do i = 1, nspaces

           if( pno_cv(i)%allocd )then
              call free_PNOSpaceInfo(pno_cv(i))
           endif

           do j = 1, i - 1
             cc = (j - i + 1) + i*(i-1)/2
             if( pno_S(cc)%allocd )  call free_PNOSpaceInfo( pno_S(cc) )
           enddo
         enddo

         call mem_dealloc( pno_cv )

       else
         do i = 1, nspaces

           if( fraginfo%CLocPNO(i)%allocd )then
              call free_PNOSpaceInfo( fraginfo%CLocPNO(i) )
           endif

           do j = 1, i - 1
             cc = (j - i + 1) + i*(i-1)/2
             if( pno_S(cc)%allocd )  call free_PNOSpaceInfo( pno_S(cc) )
           enddo
         enddo

         call mem_dealloc( fraginfo%CLocPNO )
         pno_cv => null()
       endif

       call mem_dealloc( pno_S )

     endif


     !transform back to original basis   
     if(DECinfo%use_singles)then
       call can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
       &vovo=t2_final%val,vo=t1_final%val)
       call can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
          &vovo=VOVO%val)
       if(get_mult)then
          call can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
           &vovo=m2%val,vo=m1%val)
       endif
     else
       call can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
       &vovo=t2_final%val)
       call can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
       &vovo=VOVO%val)
       if(get_mult)then
          call can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
           &vovo=m2%val)
       endif
     endif

     call mem_dealloc(Uocc)
     call mem_dealloc(Uvirt)
     !if ((ccmodel==MODEL_RPA)) call mem_dealloc(govov)


   end subroutine ccsolver_debug




   !> \brief Singles residual for CCSD model
   subroutine getSinglesResidualCCSD(omega1,u,gao,pqfock,qpfock, & 
               xocc,xvirt,yocc,yvirt,nocc,nvirt)
 
     implicit none
     type(array2), intent(inout) :: omega1
     type(array4), intent(inout) :: u,gao
     type(array2), intent(inout) :: pqfock,qpfock
     type(array2), intent(inout) :: xocc,xvirt,yocc,yvirt
     type(array2) :: tmp
     type(array4) :: qqpq, pppq
     integer, intent(in) :: nocc,nvirt
     integer :: a,i,k,l,c,d
     real(realk) :: starttime,endtime,aStart,aEnd,bStart,bEnd,cStart,cEnd, &
               dStart,dEnd
 
    
     aStart=0.0E0_realk; aEnd=0.0E0_realk
     bStart=0.0E0_realk; bEnd=0.0E0_realk
     cStart=0.0E0_realk; cEnd=0.0E0_realk
     dStart=0.0E0_realk; dEnd=0.0E0_realk
     starttime=0.0E0_realk; endtime=0.0E0_realk
 
     call cpu_time(starttime)
 
#ifdef SINGLES_EXTRA_SIMPLE    
 
     qqpq = get_gmo_simple(gao,xvirt,yvirt,xocc,yvirt)
 
     do i=1,nocc
       do a=1,nvirt
 
         do c=1,nvirt
         do k=1,nocc
         do d=1,nvirt
         omega1%val(a,i) = omega1%val(a,i) + & 
                 u%val(c,k,d,i)*qqpq%val(a,d,k,c)
         end do
         end do
         end do
 
       end do
     end do
     print *,'A1 ',omega1*omega1
 
     call array4_free(qqpq)
 
     pppq = get_gmo_simple(gao,xocc,yocc,xocc,yvirt)
 
     do i=1,nocc
       do a=1,nvirt
 
         do c=1,nvirt
         do k=1,nocc
         do l=1,nocc
 
         omega1%val(a,i) = omega1%val(a,i) - &
             u%val(a,k,c,l) * pppq%val(k,i,l,c)
 
         end do
         end do
         end do
 
       end do
     end do
     print *,'B1 ',omega1*omega1
 
     call array4_free(pppq)
 
     do i=1,nocc
       do a=1,nvirt
       
         do c=1,nvirt
           do k=1,nocc
             omega1%val(a,i) = omega1%val(a,i) + &
                 u%val(a,i,c,k) * pqfock%val(k,c)
           end do
         end do  
 
       end do
     end do
     print *,'C1 ',omega1*omega1
 
     omega1%val = omega1%val + qpfock%val
     print *,'D1 ',omega1*omega1
 
#else
 
     tmp = array2_init([nvirt,nocc])
 
     ! A1
     call cpu_time(aStart)
     qqpq = get_gmo_simple(gao,xocc,yvirt,yvirt,xvirt)
     call array4_reorder(qqpq,[2,1,3,4]) ! qqpq[ic,ab] -> qqpq[ci,ab]
     call array4_contract3(qqpq,u,tmp)
     call array4_free(qqpq)
     call array2_add_to(omega1,1.0E0_realk,tmp)
     call cpu_time(aEnd)
     if(DECinfo%cc_driver_debug) then
       print *,'A1 done, norm : ',omega1*omega1
     end if
 
     ! B1
     call cpu_time(bStart)
     pppq = get_gmo_simple(gao,xocc,yocc,xocc,yvirt)
     call array4_reorder(u,[4,3,2,1])  ! u[ai,bj] -> u[jb,ia]
     call array4_reorder(pppq,[3,4,1,2])
     call array2_zero(tmp)
     call array4_contract3(u,pppq,tmp)
     call array2_add_to(omega1,-1.0E0_realk,tmp)
     call array4_free(pppq)
     call cpu_time(bEnd)
     if(DECinfo%cc_driver_debug) then
       print *,'B1 done, norm: ',omega1*omega1
     end if
 
     ! C1
     call cpu_time(cStart)
     call array2_zero(tmp)
     call array4_reorder(u,[1,2,4,3]) ! u[jb,ia] -> u[jb,ai] 
     call array4_contract_array2(u,pqfock,tmp)
     call array4_reorder(u,[3,4,2,1])
     call array2_add_to(omega1,1.0E0_realk,tmp)
     call cpu_time(cEnd)
     if(DECinfo%cc_driver_debug) then
       print *,'C1 done, norm: ',omega1*omega1
     end if
 
     ! D1
     call cpu_time(dStart)
     call array2_add_to(omega1,1.0E0_realk,qpfock)
     call cpu_time(dEnd)
     if(DECinfo%cc_driver_debug) then
       print *,'D1 done, norm: ',omega1*omega1
     end if
 
     call array2_free(tmp)
     call cpu_time(endtime)
     if(DECinfo%PL>1) then
       write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD A1 : ',aEnd-aStart,' s'
       write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD B1 : ',bEnd-bStart,' s'
       write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD C1 : ',cEnd-cStart,' s'
       write(DECinfo%output,'(a,f16.3,a)') ' time :: CCSD D1 : ',dEnd-dStart,' s'
       write(DECinfo%output,'(a,f16.3,a)')  &
          ' time :: CCSD singles : ',endtime-starttime,' s'
     end if
 
#endif
 
     return
   end subroutine getSinglesResidualCCSD


  
   !> \brief Simple double residual for CCSD
   subroutine getDoublesResidualCCSD_simple(omega2,t2,u,gao,aibj,iajb,nocc,nvirt, &
        ppfock,qqfock,xocc,xvirt,yocc,yvirt)

     implicit none
     real(realk) :: aStart,aEnd,bStart,bEnd,cStart,cEnd, &
          dStart,dEnd,eStart,eEnd
     type(array4), intent(inout) :: omega2,t2
     type(array4), intent(inout) :: u,gao,aibj,iajb
     type(array2), intent(inout) :: ppfock, qqfock
     integer, intent(in) :: nocc,nvirt
     type(array2), intent(inout) :: xocc,xvirt,yocc,yvirt
     type(array4) :: abcd, tmp1, X
     type(array4) :: l1, l2, tmp
     type(array2) :: ppX,qqY,pptmp,qqtmp
     integer :: a,i,b,j,k,l,c,d
     character(tensor_MSG_LEN) :: msg

     aStart=0.0E0_realk; aEnd=0.0E0_realk
     bStart=0.0E0_realk; bEnd=0.0E0_realk
     cStart=0.0E0_realk; cEnd=0.0E0_realk
     dStart=0.0E0_realk; dEnd=0.0E0_realk
     eStart=0.0E0_realk; eEnd=0.0E0_realk


     ! -- A2
     call cpu_time(aStart)
     call array4_add_to(omega2,1.0E0_realk,aibj)

     abcd = get_gmo_simple(gao,xvirt,yvirt,xvirt,yvirt)
     call array4_reorder(t2,[1,3,2,4]) ! -> t2[ab,ij]
     call array4_reorder(abcd,[2,4,1,3])
     tmp1 = array4_init([nvirt,nvirt,nocc,nocc]) ! tmp1[ab,ij]
     call array4_contract2(abcd,t2,tmp1)
     call array4_reorder(tmp1,[1,3,2,4]) ! -> tmp1[ai,bj]
     call array4_add_to(omega2,1.0E0_realk,tmp1)
     call array4_free(tmp1)
     call array4_reorder(t2,[1,3,2,4]) ! -> t2[ai,bj]
     call array4_free(abcd)
     call cpu_time(aEnd)
     if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: A2 done, norm :',omega2*omega2


     ! -- B2
     call cpu_time(bStart)
     call array4_reorder(iajb,[2,4,1,3]) ! iajb[kc,ld] -> iajb[cd,kl] (iajb[ia,jb] -> iajb[ab,ij])
     call array4_reorder(t2,[1,3,2,4]) ! t2[ci,dj] -> t2[cd,ij] (t2[ai,bj] -> t2[ab,ij])
     tmp1 = array4_init([nocc,nocc,nocc,nocc])
     call array4_contract2(t2,iajb,tmp1) ! tmp1[ij,kl]
     X = get_gmo_simple(gao,xocc,yocc,xocc,yocc) ! X[ki,lj]
     call array4_reorder(X,[2,4,1,3]) ! X[ki,lj] -> X[ij,kl]
     call array4_add_to(X,1.0E0_realk,tmp1)
     call array4_free(tmp1)

     call array4_reorder(X,[3,4,1,2]) ! X[ij,kl] -> X[kl,ij]
     call array4_reorder(t2,[3,4,1,2]) ! t2[ab,kl]-> t2[kl,ab]
     tmp1 = array4_init([nvirt,nvirt,nocc,nocc])
     call array4_contract2(t2,X,tmp1) ! tmp1[ab,ij]
     call array4_reorder(tmp1,[1,3,2,4]) ! tmp1[ab,ij] -> tmp1[ai,bj]
     call array4_add_to(omega2,1.0E0_realk,tmp1)
     call array4_free(X)
     call array4_free(tmp1)
     call cpu_time(bEnd)
     if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: B2 done, norm :',omega2*omega2

     !! -- C2
     !call cpu_time(cStart)

     !! std                                  ij,ab  ->    bi,aj
     call array4_reorder(t2,[4,1,3,2]) ! t2[li,ad] -> t2[dl,ai]
     call array4_reorder(iajb,[1,4,3,2]) ! iajb[dc,kl] -> iajb[dl,kc]
     X = get_gmo_simple(gao,xocc,yocc,xvirt,yvirt) ! X[ki,ac]
     tmp1 = array4_init([nvirt,nocc,nocc,nvirt]) ! tmp1[ai,kc]
     call array4_contract2(t2,iajb,tmp1)
     call array4_reorder(tmp1,[3,2,1,4]) ! -> tmp1[ki,ac]
     call array4_add_to(X,-0.5E0_realk,tmp1)
     call array4_free(tmp1)

     !! a
     call array4_reorder(X,[4,1,3,2]) ! X[ki,ac] -> X[ck,ai]
     tmp1 = array4_init([nvirt,nocc,nvirt,nocc])
     call array4_contract2(X,t2,tmp1) ! tmp1[ai,bj]
     call array4_add_to(omega2,-0.5E0_realk,tmp1)

     !! b
     call array4_reorder(tmp1,[1,4,3,2]) ! tmp1[aj,bi] -> tmp[ai,bj]
     call array4_add_to(omega2,-1.0E0_realk,tmp1)

     !! c
     call array4_reorder(tmp1,[3,2,1,4]) ! tmp1[] -> tmp1[]
     call array4_add_to(omega2,-0.5E0_realk,tmp1)

     !! d
     call array4_reorder(tmp1,[1,4,3,2]) ! tmp[] -> tmp1[]
     call array4_add_to(omega2,-1.0E0_realk,tmp1)

     call array4_free(tmp1)
     call array4_free(X)
     call array4_reorder(t2,[3,2,1,4]) ! t2[dl,ai] -> t2[al,di]
     call cpu_time(cEnd)
     if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: C2 done, norm :',omega2*omega2

     !! -- D2
     call cpu_time(dStart)
     l1 = getL(gao,xvirt,yocc,xocc,yvirt, &
          xvirt,yvirt,xocc,yocc)
     call array4_reorder(iajb,[3,1,2,4]) ! iajb[dl,kc] -> iajb[kd,lc]
     l2 = getL(iajb)
     X = array4_init([nvirt,nocc,nocc,nvirt]) ! X[ai,kc]
     call array4_reorder(u,[4,3,1,2]) ! u[ai,dl] -> u[ld,ai]
     call array4_contract2(u,l2,X) ! X[ai,kc]
     call array4_scale(X,0.5E0_realk)
     call array4_add_to(X,1.0E0_realk,l1)

     if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: l1 norm : ',l1*l1
     if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: l1 norm : ',l2*l2

     call array4_free(l1)
     call array4_free(l2)

     !! a
     call array4_reorder(X,[3,4,1,2]) ! X[ai,kc] -> X[kc,ai]
     tmp = array4_init([nvirt,nocc,nvirt,nocc])
     call array4_contract2(X,u,tmp)
     call array4_add_to(omega2,0.5E0_realk,tmp)
     call array4_reorder(u,[3,4,2,1]) ! u[ld,ai] -> u[ai,dl]

     !! b
     call array4_reorder(tmp,[3,4,1,2]) ! tmp[bj,ai] -> tmp[ai,bj]
     call array4_add_to(omega2,0.5E0_realk,tmp)
     call array4_free(X)
     call array4_free(tmp)

     call cpu_time(dEnd)
     if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: D2 done, norm',omega2*omega2

     !! -- E2
     call cpu_time(eStart)
     qqtmp = array2_init([nvirt,nvirt])
     pptmp = array2_init([nocc,nocc])

     !! get qqY
     call array4_reorder(u,[4,3,2,1]) ! u[bk,dl] -> u[ldk,b]
     call array4_contract3(u,iajb,qqtmp)
     qqY = array2_add(1.0E0_realk,qqfock,-1.0E0_realk,qqtmp)
     call array2_free(qqtmp)

     !! get ppX
     call array4_reorder(u,[4,3,2,1]) ! u[jdl,c] -> u[cld,j]
     call array4_reorder(iajb,[4,3,2,1]) ! iajb[kdl,c] -> iajb[cld,k]
     call array4_contract3(iajb,u,pptmp)
     ppX = array2_add(1.0E0_realk,ppfock,1.0E0_realk,pptmp)
     call array2_free(pptmp)
     call array4_reorder(iajb,[4,3,2,1])

     !! 1
     call array2_transpose(qqY)
     call array4_reorder(t2,[3,4,1,2])
     tmp = array4_init([nvirt,nocc,nvirt,nocc])
     call array4_contract1(t2,qqY,tmp,.true.)
     call array4_reorder(tmp,[3,4,1,2])
     call array4_add_to(omega2,1.0E0_realk,tmp)
     call array4_free(tmp)
     call array4_reorder(t2,[3,4,1,2])

     !! 2
     call array4_contract1(t2,qqY,omega2,.false.)
     call array2_transpose(qqY)

     !! 3
     call array4_reorder(t2,[4,3,2,1])
     tmp = array4_init([nocc,nvirt,nocc,nvirt])
     call array4_contract1(t2,ppX,tmp,.true.)
     call array4_reorder(t2,[4,3,2,1])
     call array4_reorder(tmp,[4,3,2,1])
     call array4_add_to(omega2,-1.0E0_realk,tmp)
     call array4_free(tmp)

     !! 4
     call array4_reorder(t2,[2,1,3,4])
     tmp = array4_init([nocc,nvirt,nvirt,nocc])
     call array4_contract1(t2,ppX,tmp,.true.)
     call array4_reorder(t2,[2,1,3,4])
     call array4_reorder(tmp,[2,1,3,4])
     call array4_add_to(omega2,-1.0E0_realk,tmp)
     call array4_free(tmp)

     call array2_free(ppX)
     call array2_free(qqY)
     call cpu_time(eEnd)

     if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') 'debug :: E2 done, norm',omega2*omega2

   end subroutine getDoublesResidualCCSD_simple


     

    !> \brief should be a general subroutine to get the guess amplitudes when
    !starting up a CCSD or CC2 calculation and checks for files which contain
    !amplitudes of a previous calculation. If none are found the usual zero guess
    !is returned
    !> \author Patrick Ettenhuber
    !> \date December 2012
   subroutine get_guess_vectors_simple(mylsitem,gao,&
   &get_mult,Co,Co2,Cv2,no,nv,nb,xocc,yvirt,restart,t2,t1,t1f)
     implicit none
     type(lsitem),intent(inout) :: mylsitem
     !> contains the guess doubles amplitudes on output
     type(array4),intent(inout),optional :: t2
     !> contains the singles amplitudes on output
     type(array2),intent(inout), optional :: t1,t1f
     type(array2),intent(inout) :: Co,Co2,Cv2,xocc,yvirt
     type(array4),intent(inout) :: gao
     integer, intent(in)        :: no,nv,nb
     logical,intent(in)         :: get_mult
     integer                    :: d2(4), d1(2)
     type(array2) :: fockguess,t1tmp
     type(array4) :: iajb
     logical,intent(out) :: restart
     !> the filenames to check for valid singles amplitudes
     character(3):: safefilet11,safefilet12
     !> the filenames to check for valid doubles amplitudes
     character(3):: safefilet21,safefilet22
     integer :: fu_t11,fu_t12,fu_t21,fu_t22,fu_t1,fu_t2
     logical :: file_exists11,file_exists12,file_exists21,file_exists22
     logical :: file_status11,file_status12,file_status21,file_status22
     logical :: readfile1, readfile2
     integer :: saved_iter11,saved_iter12,saved_iter21,saved_iter22,iter_start
     integer :: saved_nel11,saved_nel12,saved_nel21,saved_nel22
     character(11) :: fullname11, fullname12, fullname21, fullname22
     character(tensor_MSG_LEN) :: msg
     if(get_mult)then
       safefilet21 = 'm21'
       safefilet22 = 'm22'
       safefilet11 = 'm11'
       safefilet12 = 'm12'
     else
       safefilet21 = 't21'
       safefilet22 = 't22'
       safefilet11 = 't11'
       safefilet12 = 't12'
     endif


     !print *,"CHECK INPUT",safefilet11,safefilet12,all_singles,DECinfo%use_singles
     fu_t11=111
     fu_t12=112
     fu_t21=121
     fu_t22=122

     d1 = [nv,no]
     d2 = [nv,no,nv,no]
     if(DECinfo%use_singles)then
       t1 = array2_init(d1)
     endif
     t2 = array4_init(d2)

     


      iter_start=0
      !check for safe files of the amplitudes in the current directory and read
      !them if they exist and ok
      readfile1=.false.
      readfile2=.false.
      
      !this can be skipped by input, but restart is default
      if(DECinfo%DECrestart)then
        if(DECinfo%use_singles)then
          fullname11=safefilet11//'.restart'
          fullname12=safefilet12//'.restart'

          file_status11=.false.
          INQUIRE(FILE=fullname11,EXIST=file_exists11)
          if(file_exists11)then
            file_status11=.true.
            OPEN(fu_t11,FILE=fullname11,STATUS='OLD',FORM='UNFORMATTED')
            READ(fu_t11)saved_iter11
            READ(fu_t11)saved_nel11
            if(saved_nel11/=no*nv)then
              call lsquit("ERROR(ccsolver_debug):wrong dimensions in amplitude &
              &file",DECinfo%output)
            endif
          endif
         
          file_status12=.false.
          INQUIRE(FILE=fullname12,EXIST=file_exists12)
          if(file_exists12)then
            file_status12=.true.
            OPEN(fu_t12,FILE=fullname12,STATUS='OLD',FORM='UNFORMATTED')
            READ(fu_t12)saved_iter12
            READ(fu_t12)saved_nel12
            if(saved_nel12/=no*nv)then
              call lsquit("ERROR(ccsolver_debug):wrong dimensions in amplitude &
              &file",DECinfo%output)
            endif
          endif
         
          !CHECK WHICH IS THE PREFERRED FILE TO READ
          if(file_status11.and.file_status12)then
            if(saved_iter11>saved_iter12)then
              fu_t1=fu_t11
              CLOSE(fu_t12)
              readfile1=.true.
            else
              fu_t1=fu_t12
              CLOSE(fu_t11)
              readfile1=.true.
            endif
          elseif(file_status11)then
            fu_t1=fu_t11
            readfile1=.true.
          elseif(file_status12)then
            fu_t1=fu_t12
            readfile1=.true.
          endif  
        endif

        fullname21=safefilet21//'.restart'
        fullname22=safefilet22//'.restart'
       
        file_status21=.false.
        INQUIRE(FILE=fullname21,EXIST=file_exists21)
        if(file_exists21)then
          file_status21=.true.
          OPEN(fu_t21,FILE=fullname21,STATUS='OLD',FORM='UNFORMATTED')
          READ(fu_t21)saved_iter21
          READ(fu_t21)saved_nel21
          if(saved_nel21/=no*no*nv*nv)then
            call lsquit("ERROR(ccsolver_debug):wrong dimensions in amplitude &
            &file",DECinfo%output)
          endif
        endif
       
        file_status22=.false.
        INQUIRE(FILE=fullname22,EXIST=file_exists22)
        if(file_exists22)then
          file_status22=.true.
          OPEN(fu_t22,FILE=fullname22,STATUS='OLD',FORM='UNFORMATTED')
          READ(fu_t22)saved_iter22
          READ(fu_t22)saved_nel22
          if(saved_nel22/=no*no*nv*nv)then
            call lsquit("ERROR(ccsolver_debug):wrong dimensions in amplitude &
            &file",DECinfo%output)
          endif
        endif
       
        !CHECK WHICH IS THE PREFERRED FILE TO READ
        if(file_status21.and.file_status22)then
          if(saved_iter21>saved_iter22)then
            iter_start=saved_iter21
            fu_t2=fu_t21
            CLOSE(fu_t22)
            readfile2=.true.
          else
            iter_start=saved_iter22
            fu_t2=fu_t22
            CLOSE(fu_t21)
            readfile2=.true.
          endif
          WRITE(DECinfo%output,'("RESTARTING CC CALCULATION WITH TRIAL-VECS FROM: ",I3)')iter_start
        elseif(file_status21)then
          iter_start=saved_iter21
          fu_t2=fu_t21
          WRITE(DECinfo%output,'("RESTARTING CC CALCULATION WITH TRIAL-VECS FROM: ",I3)')iter_start
          readfile2=.true.
        elseif(file_status22)then
          iter_start=saved_iter22
          fu_t2=fu_t22
          WRITE(DECinfo%output,'("RESTARTING CC CALCULATION WITH TRIAL-VECS FROM: ",I3)')iter_start
          readfile2=.true.
        else
          iter_start=1
        endif  
      endif


      if(readfile1)then
        READ(fu_t1)t1%val
        CLOSE(fu_t1)
        restart = .true.
      else
        restart = .false.
        if(get_mult)then
          fockguess=array2_init([nb,nb])
          call Get_AOt1Fock(mylsitem,t1f,fockguess,no,nv,nb,Co,Co2,Cv2)
          t1tmp = array2_similarity_transformation(xocc,fockguess,yvirt,[no,nv]) 
          call array2_free(fockguess)
          call mat_transpose(no,nv,-1.0E0_realk,t1tmp%val,0.0E0_realk,t1%val)
          call array2_free(t1tmp)
          call dscal(no*nv,2.0E0_realk,t1%val,1)
        endif 
      endif
      if(readfile2)then
        READ(fu_t2) t2%val
        CLOSE(fu_t2)
        restart = .true.
      else
        restart = .false.
        if(get_mult)then
          iajb = get_gmo_simple(gao,xocc,yvirt,xocc,yvirt)
          call array4_reorder(iajb,[2,1,4,3])
          call daxpy(no**2*nv**2,-4.0E0_realk,iajb%val,1,t2%val,1)
          call array4_reorder(iajb,[1,4,3,2])
          call daxpy(no**2*nv**2,2.0E0_realk,iajb%val,1,t2%val,1)
          call array4_reorder(iajb,[4,1,2,3])
          call array4_free(iajb)
        endif
      endif
   end subroutine get_guess_vectors_simple
   !> \brief Subroutine to save the current guess amplitudes for the next
   !iteration
   !> \author Patrick Ettenhuber
   !> \date Dezember 2012
   subroutine save_current_guess_simple(iter,get_mult,t2,t1)
    implicit none
    !> iteration number
    integer,intent(in) :: iter
    !> doubles guess amplitudes for the next iteration
    type(array4), intent(in),optional :: t2
    !> singles guess amplitudes for the next iteration
    type(array2), intent(in),optional :: t1
    logical, intent(in) :: get_mult
    !> alternating filenames for the doubles amplitudes
    character(3) :: safefilet21,safefilet22
    !> alternating filenames for the singles amplitudes
    character(3) :: safefilet11,safefilet12
    integer :: fu_t21,fu_t22
    integer :: fu_t11,fu_t12
    logical :: file_status11,file_status12,file_status21,file_status22
    logical :: all_singles
    character(tensor_MSG_LEN) :: msg
#ifdef SYS_AIX
    character(12) :: fullname11,  fullname12,  fullname21,  fullname22
    character(12) :: fullname11D, fullname12D, fullname21D, fullname22D
#else
    character(11) :: fullname11, fullname12, fullname21, fullname22
    character(11) :: fullname11D, fullname12D, fullname21D, fullname22D
#endif
    if(get_mult)then
      safefilet21 = 'm21'
      safefilet22 = 'm22'
      safefilet11 = 'm11'
      safefilet12 = 'm12'
    else
      safefilet21 = 't21'
      safefilet22 = 't22'
      safefilet11 = 't11'
      safefilet12 = 't12'
    endif
    fu_t11=111
    fu_t12=112
    fu_t21=121
    fu_t22=122
    if(DECinfo%use_singles)then
      !msg="singles norm save"
      !call print_norm(t1,msg)
#ifdef SYS_AIX
      fullname11=safefilet11//'.writing\0'
      fullname12=safefilet12//'.writing\0'
#else
      fullname11=safefilet11//'.writing'
      fullname12=safefilet12//'.writing'
#endif

      if(mod(iter,2)==1)then
        file_status11=.false. 
        OPEN(fu_t11,FILE=fullname11,STATUS='REPLACE',FORM='UNFORMATTED')
        WRITE(fu_t11)iter
        WRITE(fu_t11)t1%dims(1) * t1%dims(2)
        WRITE(fu_t11)t1%val
        file_status11=.true.
        WRITE(fu_t11)file_status11
        ENDFILE(fu_t11)
        CLOSE(fu_t11)
#ifdef SYS_AIX
        fullname11D=safefilet11//'.restart\0'
        fullname11=safefilet11//'\0'
#else
        fullname11D=safefilet11//'.restart'
#endif
       if(file_status11)call rename(fullname11,fullname11D)

      elseif(mod(iter,2)==0)then
        file_status12=.false. 
        OPEN(fu_t12,FILE=fullname12,STATUS='REPLACE',FORM='UNFORMATTED')
        WRITE(fu_t12)iter
        WRITE(fu_t12)t1%dims(1) * t1%dims(2)
        WRITE(fu_t12)t1%val
        file_status12=.true.
        WRITE(fu_t12)file_status12
        ENDFILE(fu_t12)
        CLOSE(fu_t12)
#ifdef SYS_AIX
        fullname12D=safefilet12//'.restart\0'
        fullname12=safefilet12//'\0'
#else
        fullname12D=safefilet12//'.restart'
#endif
        if(file_status12)call rename(fullname12,fullname12D)

      else
        call lsquit("ERROR(ccdriver_par):impossible iteration&
        &number)",DECinfo%output)
      endif
    endif

    !msg="doubles norm save"
    !call print_norm(t2,msg)
    fullname21=safefilet21//'.writing'
    fullname22=safefilet22//'.writing'
    if(mod(iter,2)==1)then
      file_status21=.false. 
      OPEN(fu_t21,FILE=fullname21,STATUS='REPLACE',FORM='UNFORMATTED')
      WRITE(fu_t21)iter
      WRITE(fu_t21)t2%dims(1) * t2%dims(2) * t2%dims(3) * t2%dims(4)
      WRITE(fu_t21)t2%val
      file_status21=.true. 
      WRITE(fu_t22)file_status21
      ENDFILE(fu_t21)
      CLOSE(fu_t21)
#ifdef SYS_AIX
      fullname21D=safefilet21//'.restart\0'
      fullname21=safefilet21//'\0'
#else
      fullname21D=safefilet21//'.restart'
#endif
      if(file_status21)call rename(fullname21,fullname21D)

    elseif(mod(iter,2)==0)then
      file_status22=.false. 
      OPEN(fu_t22,FILE=fullname22,STATUS='REPLACE',FORM='UNFORMATTED')
      WRITE(fu_t22)iter
      WRITE(fu_t22)t2%dims(1) * t2%dims(2) * t2%dims(3) * t2%dims(4)
      WRITE(fu_t22)t2%val
      file_status22=.true.
      WRITE(fu_t22)file_status22
      ENDFILE(fu_t22)
      CLOSE(fu_t22)
#ifdef SYS_AIX
      fullname22D=safefilet22//'.restart\0'
      fullname22=safefilet22//'\0'
#else
      fullname22D=safefilet22//'.restart'
#endif
      if(file_status22)call rename(fullname22,fullname22D)
    else
      call lsquit("ERROR(ccdriver_par):impossible iteration&
      &number)",DECinfo%output)
    endif
  end subroutine save_current_guess_simple
#else
  subroutine cc_debug_void()
     implicit none
  end subroutine cc_debug_void
#endif

end module cc_debug_routines_module
