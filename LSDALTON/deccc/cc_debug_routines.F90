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
   use tensor_type_def_module
! begin pablo
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif
! end pablo 
 

   ! DEC DEPENDENCIES (within deccc directory)   
   ! *****************************************
   use crop_tools_module
   use array2_simple_operations
   use array4_simple_operations
   use ri_simple_operations
   use mp2_module
   use ccintegrals
   use ccsd_module
   use ccsdpt_module
   use orbital_operations
   use rpa_module
   

   contains
   subroutine ccsolver_energy_multipliers(ccmodel,Co_f,Cv_f,fock_f,nbasis,nocc,nvirt, &
        &mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f,ccenergy)

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
     !> Occupied MO coefficients  for fragment/full molecule
     real(realk), dimension(nbasis,nocc), intent(inout) :: Co_f
     !> Virtual MO coefficients  for fragment/full molecule
     real(realk), dimension(nbasis,nvirt), intent(inout) :: Cv_f
     !> Occ-occ block of Fock matrix in MO basis
     real(realk), dimension(nocc,nocc), intent(inout) :: ppfock_f
     !> Virt-virt block of Fock matrix in MO basis
     real(realk), dimension(nvirt,nvirt), intent(inout) :: qqfock_f
     !> Is this a fragment job (true) or a full molecular calculation (false)
     logical, intent(in) :: fragment_job
     !> LS item information
     type(lsitem), intent(inout) :: mylsitem
     !> How much to print? ( ccPrintLevel>0 --> print info stuff)
     integer, intent(in) :: ccPrintLevel
     !> Coupled cluster energy for fragment/full molecule
     real(realk),intent(inout) :: ccenergy!,ccsdpt_e4,ccsdpt_e5,ccsdpt_tot
     type(array4) :: t2_final,VOVO!,ccsdpt_t2
     type(array2) :: t1_final!,ccsdpt_t1
     !> stuff needed for pair analysis
     type(array2) :: ccsd_mat_tot,ccsd_mat_tmp
     integer :: natoms,ncore,nocc_tot,p,pdx,i
     type(decorbital), pointer :: occ_orbitals(:)
     type(decorbital), pointer :: unocc_orbitals(:)
     logical, pointer :: orbitals_assigned(:)
     type(array4) :: mult2
     type(array2) :: mult1


     call ccsolver_debug(ccmodel,Co_f, Cv_f, fock_f, nbasis, nocc, nvirt, &
        & mylsitem, ccPrintLevel, fragment_job, ppfock_f, qqfock_f, ccenergy, &
        & t1_final, t2_final, VOVO, .false.)

     call array4_free(VOVO)

     call ccsolver_debug(ccmodel,Co_f, Cv_f, fock_f, nbasis, nocc, nvirt, &
        & mylsitem, ccPrintLevel, fragment_job, ppfock_f, qqfock_f, ccenergy, &
        & t1_final, t2_final, VOVO, .false., m2 = mult2, m1 = mult1)

     !call print_norm(mult2%val,int(i8*nvirt*nvirt*nocc*nocc,kind=8))
     !call print_norm(mult1%val,int(i8*nvirt*nocc,kind=8))

     ! Free arrays
     call array2_free(t1_final)
     call array2_free(mult1)
     call array4_free(t2_final)
     call array4_free(mult2)
     call array4_free(VOVO)

   end subroutine ccsolver_energy_multipliers

   !> \author Marcin Ziolkowski (modified by Kasper Kristensen and Patrick
   !  Ettenhuber)
   subroutine ccsolver_debug(ccmodel,Co_f,Cv_f,fock_f,nbasis,nocc,nvirt, &
        & mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f,ccenergy, &
        & t1_final,t2_final,VOVO,longrange_singles,m2,m1)

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
    
     type(array4),optional,intent(inout) :: m2
     type(array2),optional,intent(inout) :: m1
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
     character(18) :: save_to,keep
     character(ARR_MSG_LEN) :: msg
     integer :: ii,aa
     integer :: MaxSubSpace
     logical :: restart

     ! begin pablo 
     real(realk),pointer :: pack_gmo(:)
     type(MObatchInfo) :: MOinfo
     logical :: small_frag
     ! end pablo


     call LSTIMER('START',ttotstart_cpu,ttotstart_wall,DECinfo%output)
     if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)



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

     ! Sanity check 3: if CCSD multipliers are requested, make sure that both are
     ! there
     if((present(m2).and..not.present(m1)).or.(present(m1).and..not.present(m2)))then
       call lsquit("ERROR(ccsolver_debug):requested unkown multipliers",-1)
     endif


     get_mult    = (present(m2).and.present(m1))
     MaxSubSpace = DECinfo%ccMaxDIIS

     ! title
     Call print_ccjob_header(ccmodel,ccPrintLevel,fragment_job,get_mult,nbasis,nocc,nvirt,MaxSubSpace)

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

       ppfock_d = 0.0E0_realk
       qqfock_d = 0.0E0_realk

       do ii=1,nocc
         ppfock_d(ii,ii) = focc(ii)
       enddo
       do aa=1,nvirt
         qqfock_d(aa,aa) = fvirt(aa)
       enddo

       if(get_mult)then

         if(DECinfo%use_singles)then
           call ccsolver_local_can_trans(VOVO%val,t2_final%val,nocc,nvirt,Uocc,Uvirt,t1_final%val)
         else
           call ccsolver_local_can_trans(VOVO%val,t2_final%val,nocc,nvirt,Uocc,Uvirt)
         endif

       endif


     endif

     call mem_dealloc(focc)
     call mem_dealloc(fvirt)

     ! Copy MO coeffcients. It is very convenient to store them twice to handle transformation
     ! (including transposed MO matrices) efficiently. 
     Co2_d = Co_d
     Cv2_d = Cv_d

     ! create transformation matrices in array form
     Co   = array2_init(occ_dims,Co_d)
     Cv   = array2_init(virt_dims,Cv_d)
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
     call mem_alloc(dens,nbasis,nbasis)
     call get_density_from_occ_orbitals(nbasis,nocc,Co%val,dens)

     call mem_dealloc(dens)

     ! get two-electron integrals in ao
     if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a)') 'debug :: calculating AO integrals'
     ! Only calculate full 4-dimensional AO integrals 
     call get_full_eri(mylsitem,nbasis,gao)
     if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,/)') 'debug :: AO integrals done'

     ! Simulate two-electron integrals (debug mode)
     if(DECinfo%simulate_eri .or. DECinfo%fock_with_ri) then
        if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a)') &
             'debug :: calculate RI intermediate - temporary'
        l_ao = get_ao_ri_intermediate(mylsitem)
        if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a)') &
             'debug :: intermediates done'
     end if

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
     if(CCmodel /= MODEL_MP2) then
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

     ! readme : the iteration sequence is universal and may be used for all
     !          iterative cc models (linear or non-linear) and is
     !          semi-independent on the storage of vectors (allocation and
     !          deallocation, etc)

     ! iterate
     break_iterations = .false.
     crop_ok          = .false.
     prev_norm        = 1.0E6_realk


     ! begin pablo
     ! criterion will need to be improved/adjusted.
     !if (nbasis<=300) small_frag=.true. 
     !> get gmo and packed them
     if (small_frag) then
       call get_packed_gmo(mylsitem,Co%val,Cv2%val,pack_gmo, &
            & nbasis,nocc,nvirt,MOinfo)
     end if
     ! end pablo



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

           call get_guess_vectors_simple(mylsitem,t2(iter),t1(iter),&
           &t1_final,gao,get_mult,Co,Co2,Cv2,nocc,nvirt,nbasis,xocc,yvirt,restart)

        end if GetGuessVectors

        ! Initialize residual vectors
        if(DECinfo%use_singles) omega1(iter) = array2_init(ampl2_dims)
        if(DECinfo%array4OnFile) then
           ! KK, initialize omega2(iter) using storing type 2
           omega2(iter) = array4_init(ampl4_dims,2,.true.)
        else
           omega2(iter) = array4_init(ampl4_dims)
        endif

        ! get singles
        
        T1Related : if(DECinfo%use_singles) then

           ! get the T1 transformation matrices
           if(.not.get_mult)then
              call getT1transformation(t1(iter),xocc,xvirt,yocc,yvirt, &
                Co,Cv,Co2,Cv2)
           endif

           ! get inactive fock
           if(DECinfo%fock_with_ri) then
              ! Debug mode
              ifock = getInactiveFockFromRI(l_ao,xocc,yocc,h1)
           else
              ifock = getInactiveFock_simple(h1,gao,xocc,yocc,nocc,nbasis)
           end if
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

              call get_ccsd_multipliers_simple(omega1(iter)%val,omega2(iter)%val,t1_final%val&
              &,t2_final%val,t1(iter)%val,t2(iter)%val,gao,xocc%val,yocc%val,xvirt%val,yvirt%val&
              &,nocc,nvirt,nbasis,MyLsItem)
           
           ! begin pablo
           !> call CCSD code for small fragment:
           else if (small_frag) then

             ! reorder array like in Patrick's code
             call array4_reorder(t2(iter),[1,3,2,4]) ! -> t2[ab,ij]
             call array4_reorder(omega2(iter),[1,3,2,4]) ! -> om2[ab,ij]

             call get_ccsd_residual_small_frag(pack_gmo,t1(iter)%val, &
                  & t2(iter)%val,omega2(iter)%val,nbasis,nocc,nvirt,MOinfo)

             ! Calculate and print simple A2.2 residual for debug:
             call ultra_simple_a22(gao,xocc,xvirt,yocc,yvirt,t1(iter),t2(iter),nvirt,nocc)

             ! restor previous order:
             call array4_reorder(omega2(iter),[1,3,2,4]) ! -> om2[ai,bj]
             call array4_reorder(t2(iter),[1,3,2,4]) ! -> t2[ai,bj]

             call lsquit('CC iteration cannot continue because the CCSD &
                        & residual is not fully calculated',DECinfo%output)
           ! end pablo
           else

             u = get_u(t2(iter))
             call getSinglesResidualCCSD(omega1(iter),u,gao,pqfock,qpfock,xocc,xvirt,yocc,yvirt,nocc,nvirt)
             aibj = get_gmo_simple(gao,yocc,xvirt,yocc,xvirt)
             call array4_reorder(aibj,[2,1,4,3])
             call getDoublesResidualCCSD_simple(omega2(iter),t2(iter),u,gao,aibj,iajb,nocc,nvirt, &
                  ppfock,qqfock,xocc,xvirt,yocc,yvirt)

           endif
           call array4_free(aibj)
           call array4_free(u)


        elseif(CCmodel==MODEL_RPA) then

           call RPA_residual(Omega2(iter),t2(iter),gmo,ppfock,qqfock,nocc,nvirt)

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
             DECinfo%cc_driver_debug)

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
           EnergyForCCmodel: if(CCmodel==MODEL_MP2) then  
              ! MP2
              ccenergy = get_mp2_energy(t2(iter),Lmo)
           elseif(CCmodel==MODEL_CC2 .or. CCmodel==MODEL_CCSD .or. CCmodel==MODEL_CCSDpT )then
              ! CC2, CCSD, or CCSD(T) (for (T) calculate CCSD contribution here)
              ccenergy = get_cc_energy(t1(iter),t2(iter),iajb,nocc,nvirt)
           elseif(CCmodel==MODEL_RPA) then
              ccenergy = RPA_energy(t2(iter),gmo)
              sosex = SOSEX_contribution(t2(iter),gmo)
              ccenergy=ccenergy+sosex
           end if EnergyForCCmodel
        endif


        if(DECinfo%PL>1) call LSTIMER('CCIT: ENERGY',tcpu,twall,DECinfo%output)

        ! check if this is the last iteration
        if(iter == DECinfo%ccMaxIter .or. two_norm_total < DECinfo%ccConvergenceThreshold) &
             break_iterations=.true.

        if(DECinfo%use_singles .and. (.not. break_iterations) ) then
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
           
           if(.not.DECinfo%CCSDnosaferun)&
             &call save_current_guess_simple(iter,t2(iter),t1(iter),get_mult)
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
        call print_ccjob_iterinfo(iter,two_norm_total,ccenergy,get_mult)

        last_iter = iter
        if(break_iterations) exit

     end do CCIteration

     call LSTIMER('START',ttotend_cpu,ttotend_wall,DECinfo%output)



     ! Free memory and save final amplitudes
     ! *************************************


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

     ! Write finalization message
     !---------------------------
     call print_ccjob_summary(break_iterations,get_mult,fragment_job,last_iter,&
     &ccenergy,ttotend_wall,ttotstart_wall,ttotend_cpu,ttotstart_cpu,t1_final,t2_final)


     ! Save two-electron integrals in the order (virt,occ,virt,occ)
     if(CCmodel == MODEL_MP2) then
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

     if(DECinfo%simulate_eri .or. DECinfo%fock_with_ri) then
        call ri_free(l_ao)
        call ri_reset()
     end if


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


     !transform back to original basis   
     if(DECinfo%use_singles)then
       call ccsolver_can_local_trans(VOVO%val,t2_final%val,nocc,nvirt,Uocc,Uvirt,t1_final%val)
     else
       call ccsolver_can_local_trans(VOVO%val,t2_final%val,nocc,nvirt,Uocc,Uvirt)
     endif

     call mem_dealloc(Uocc)
     call mem_dealloc(Uvirt)


   end subroutine ccsolver_debug




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

   !> \author Patrick Ettenhuber
   !> \Date September 2013
   !> \brief This debug routine calculates the residual of the CCSD left
   !transformations rho1 and rho2, on input, the converged CCSD amplitudes t1f
   !and t2f must be given, furthermore a first guess for the multipliers m1 and
   !m2, the full ao integrals as array4 the t1 transformation matrices xo yo xv
   !and yv being the occupied and virtual particle and hole matrices, the number
   !of occupied orbitals no, the number of virtual orbitals nv and the number of
   !basis functions nb, MyLsItem is required to calculate the Fock matrix in
   !here
   subroutine get_ccsd_multipliers_simple(rho1,rho2,t1f,t2f,m1,m2,gao,xo,yo,xv,yv,no,nv,nb,MyLsItem)
     implicit none

     type(lsitem), intent(inout) :: MyLsItem
     real(realk),intent(inout) :: rho1(:,:),rho2(:,:,:,:)
     type(array4),intent(inout):: gao
     real(realk),intent(inout) :: t1f(:,:),t2f(:,:,:,:),m1(:,:),m2(:,:,:,:)
     real(realk),intent(in)    :: xo(:,:),yo(:,:),xv(:,:),yv(:,:)
     integer, intent(in)       :: no,nv,nb
     real(realk), pointer      :: w1(:), w2(:), w3(:), w4(:)
     real(realk), pointer      :: Lovov(:),Lvoov(:), Lovoo(:)
     real(realk), pointer      :: gvovv(:), gvooo(:), govvv(:), gooov(:), govov(:)
     real(realk), pointer      :: goovv(:)
     real(realk), pointer      :: oof(:),ovf(:),vof(:),vvf(:)
     character(ARR_MSG_LEN)    :: msg
     integer                   :: v4,o4,o3v,ov3,o2v2,ov,b2,v2,o2
     type(matrix)              :: iFock, Dens
     integer                   :: i,a,j,b,ctr
     real(realk)               :: norm,nrmt2,nrmm2

     b2   = nb*nb
     v2   = nv*nv
     o2   = no*no
     ov   = no*nv
     o2v2 = ov*ov
     ov3  = ov*v2
     o3v  = ov*o2
     o4   = o2*o2
     v4   = v2*v2

     if( DECinfo%PL>2 )then 

        write (msg,*)"t1 n**2 n"
        call print_norm(t1f,int(ov,kind=8),norm,.true.)
        print *,msg,norm,sqrt(norm)

        write (msg,*)"z1 n**2 n"
        call print_norm(m1,int(ov,kind=8),norm,.true.)
        print *,msg,norm,sqrt(norm)

        write (msg,*)"t2 n**2 n pn**2"
        call print_norm(t2f,int(o2v2,kind=8),norm,.true.)
        nrmt2 = 0.0E0_realk
        nrmm2 = 0.0E0_realk
        do j = 1, no
          do b = 1, nv
            do i = 1, no
              do a = 1,nv
                if(a+(i-1)*nv<=b+(j-1)*nv)then
                  nrmt2 = nrmt2 + t2f(a,i,b,j)**2
                  nrmm2 = nrmm2 + m2(a,i,b,j)**2
                endif
              enddo
            enddo
          enddo
        enddo
        print *,msg,norm,sqrt(norm),nrmt2
        
        write (msg,*)"z2 n**2 n pn**2 pn"
        call print_norm(m2,int(o2v2,kind=8),norm,.true.)
        print *,msg,norm,sqrt(norm),nrmm2,sqrt(nrmm2)
     endif

     call mem_alloc(w2,max(max(no,nv)*nb**3,max(max(max(max(o2v2,ov3),v4),o2*v2),o4)))

     call mem_alloc(oof,o2)
     call mem_alloc(ovf,ov)
     call mem_alloc(vof,ov)
     call mem_alloc(vvf,v2)

     call mem_alloc(govov,o2v2)
     call mem_alloc(goovv,o2v2)
     call mem_alloc(Lovov,o2v2)
     call mem_alloc(Lvoov,o2v2)


     call mem_alloc(gvovv,ov3)
     call mem_alloc(govvv,ov3)

     call mem_alloc(Lovoo,o3v)
     call mem_alloc(gvooo,o3v)
     call mem_alloc(gooov,o3v)


     !construct Ls from gs

     !govov
     call array4_read(gao)
     call successive_4ao_mo_trafo(nb,gao%val,xo,no,yv,nv,xo,no,yv,nv,w2)
     call dcopy(o2v2,gao%val,1,govov,1)

     call array_reorder_4d(2.0E0_realk,gao%val,no,nv,no,nv,[1,2,3,4],0.0E0_realk,Lovov)
     call array_reorder_4d(-1.0E0_realk,gao%val,no,nv,no,nv,[1,4,3,2],1.0E0_realk,Lovov)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gvoov
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yo,no,xo,no,yv,nv,w2)
     !Lvoov (a i j b) += gvoov (a i j b)
     call array_reorder_4d(2.0E0_realk,gao%val,nv,no,no,nv,[1,2,3,4],0.0E0_realk,Lvoov)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gvvoo
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yv,nv,xo,no,yo,no,w2)
     !write (msg,*)"gvvoo"
     !call print_norm(gao%val,int(no**2*nv**2,kind=8),msg)
     !Lvoov (a i j b) += gvvoo (a b i j)
     call array_reorder_4d(-1.0E0_realk,gao%val,nv,nv,no,no,[1,4,3,2],1.0E0_realk,Lvoov)

     call array_reorder_4d(1.0E0_realk,gao%val,nv,nv,no,no,[3,4,1,2],0.0E0_realk,goovv)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gvvov
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yv,nv,xo,no,yv,nv,w2)
     call array_reorder_4d(1.0E0_realk,gao%val,nv,nv,no,nv,[3,4,1,2],0.0E0_realk,govvv)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gvovv
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yo,no,xv,nv,yv,nv,w2)
     call dcopy(no*nv**3,gao%val,1,gvovv,1)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gooov
     call successive_4ao_mo_trafo(nb,gao%val,xo,no,yo,no,xo,no,yv,nv,w2)

     call dcopy(nv*no**3,gao%val,1,gooov,1)


     call array_reorder_4d(2.0E0_realk,gao%val,no,no,no,nv,[3,4,1,2],0.0E0_realk,Lovoo)
     call array_reorder_4d(-1.0E0_realk,gao%val,no,no,no,nv,[1,4,3,2],1.0E0_realk,Lovoo)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gvooo
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yo,no,xo,no,yo,no,w2)
     call dcopy(nv*no**3,gao%val,1,gvooo,1)

     call array4_dealloc(gao)


     !allocate the density matrix
     call mat_init(iFock,nb,nb)
     call mat_init(Dens,nb,nb)

     !calculate inactive fock matrix in ao basis
     call dgemm('n','t',nb,nb,no,1.0E0_realk,yo,nb,xo,nb,0.0E0_realk,Dens%elms,nb)
     call mat_zero(iFock)
     call dec_fock_transformation(iFock,Dens,MyLsItem,.false.)

     call ii_get_h1_mixed_full(DECinfo%output,DECinfo%output,MyLsItem%setting,&
          & Dens%elms,nb,nb,AORdefault,AORdefault)
     ! Add one- and two-electron contributions to Fock matrix
     call daxpy(b2,1.0E0_realk,Dens%elms,1,iFock%elms,1)
     !Free the density matrix
     call mat_free(Dens)

     call mem_alloc(w1,max(max(max(max(o2v2,ov3),v4),o2*v2),o4))

     !Transform inactive Fock matrix into the different mo subspaces
     ! -> Foo
     call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,iFock%elms,nb,0.0E0_realk,w1,no)
     call dgemm('n','n',no,no,nb,1.0E0_realk,w1,no,yo,nb,0.0E0_realk,oof,no)
     ! -> Fov
     call dgemm('n','n',no,nv,nb,1.0E0_realk,w1,no,yv,nb,0.0E0_realk,ovf,no)
     ! -> Fvo
     call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,iFock%elms,nb,0.0E0_realk,w1,nv)
     call dgemm('n','n',nv,no,nb,1.0E0_realk,w1,nv,yo,nb,0.0E0_realk,vof,nv)
     ! -> Fvv
     call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1,nv,yv,nb,0.0E0_realk,vvf,nv)

     call mat_free(iFock)

     call mem_alloc(w3,max(max(max(max(o2v2,ov3),v4),o2*v2),o4))
     call mem_alloc(w4,max(max(ov3,o3v),o2v2))


     !The notation in this routine is according to Halkier et. al. J. chem. phys.,
     !Vol 107. No.3 15 july 1997, and the order of terms is according to the
     !DALTON program  for simple comparison of the terms

     !SINGLES EXPRESSIONS
     !*******************

     !rho f
     !-----
     rho1 = 0.0E0_realk
     ! part1
     !sort amps(dkfj) -> dkjf
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w2)
     ! w3 : \sum_{fj} t^{df}_{kj} (d k f j) Lovvv (j f e a)
     call array_reorder_4d(2.0E0_realk,govvv,no,nv,nv,nv,[1,2,3,4],0.0E0_realk,w1)
     call array_reorder_4d(-1.0E0_realk,govvv,no,nv,nv,nv,[1,4,3,2],1.0E0_realk,w1)
     call dgemm('n','n',ov,v2,ov,1.0E0_realk,w2,ov,w1,ov,0.0E0_realk,w3,ov)
     !write(msg,*)"rho f 1"
     !call print_norm(w3,int(ov3,kind=8),msg)
     ! part2
     ! sort t2f (e j f k) -[1,4,2,3]> t2f (e k j f)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,4,2,3],0.0E0_realk,w2)
     ! sort govvv(j a d f) -[1,4,2,3]> govvv (j f a d) 
     call array_reorder_4d(1.0E0_realk,govvv,no,nv,nv,nv,[1,4,2,3],0.0E0_realk,w1)
     ! w1 : \sum_{jf} t^{ef}_{jk}(e k j f) govvv_{jadf}(j f a d)
     call dgemm('n','n',ov,v2,ov,1.0E0_realk,w2,ov,w1,ov,0.0E0_realk,w4,ov)
     !write(msg,*)"rho f 2"
     !call print_norm(w4,int(ov3,kind=8),msg)
     ! sort result w4(e k a d) -[4,2,1,3]+> w3 (d k e a)
     call array_reorder_4d(-1.0E0_realk,w4,nv,no,nv,nv,[4,2,1,3],1.0E0_realk,w3)
     !write(msg,*)"rho f 2 - added"
     !call print_norm(w3,int(ov3,kind=8),msg)
     ! part3
    
     ! sort t2f (d j f k) -[1,4,2,3]> t2f (d k j f)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,4,2,3],0.0E0_realk,w2)
     ! \sum_{jf} t^{df}_{jk} (d k j f)(still in w2) govvv(j f e a)
     call dgemm('n','n',ov,v2,ov,-1.0E0_realk,w2,ov,govvv,ov,1.0E0_realk,w3,ov)
     !write(msg,*)"rho f 3"
     !call print_norm(w3,int(ov3,kind=8),msg)

     !\sum_{dke} \zeta^{d e}_{k i}(d k e , i)^T w3(d k e a)
     call dgemm('t','n',no,nv,v2*no,1.0E0_realk,m2,v2*no,w3,v2*no,0.0E0_realk,w2,no)
     !write(msg,*)"rho f1-3(LT21I)"
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T
     call mat_transpose(no,nv,1.0E0_realk,w2,1.0E0_realk,rho1)

     !rho c - 1
     !-----
     ! part1
     ! sort \zeta^{df}_{kl} (d k f l) -[3 1 2 4]> w1(f d k l)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! sort t^{de}_{kl} (d k e l) -[1,2,4,3]> w2(d k l e)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w2)
     ! w3 : \sum_{dkl} w1 (f d k l) w2 (d k l e)
     call dgemm('n','n',nv,nv,no*no*nv,1.0E0_realk,w1,nv,w2,nv*no*no,0.0E0_realk,w3,nv)
     ! w1 : \sum_{fe} w3 (fe) L_{feia} (f e i a)
     call array_reorder_4d(2.0E0_realk,govvv,no,nv,nv,nv,[3,4,1,2],0.0E0_realk,w1)
     call array_reorder_4d(-1.0E0_realk,govvv,no,nv,nv,nv,[3,2,1,4],1.0E0_realk,w1)
     call dgemv('t',v2,ov,1.0E0_realk,w1,v2,w3,1,0.0E0_realk,w2,1)
     !write(msg,*)"rho c - 1(LT21A)"
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T
     call mat_transpose(no,nv,1.0E0_realk,w2,1.0E0_realk,rho1)
   
     if( DECinfo%PL > 2) then
        write(msg,*)"rho1 after (LT21A)"
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     !rho g 1-3
     !-----
     ! part1
     ! sort Lovoo (j f i l) -[2,1,3,4]> (f j i l)
     call array_reorder_4d(1.0E0_realk,Lovoo,no,nv,no,no,[2,1,3,4],0.0E0_realk,w2)
     ! -\sum_{jf} t^{df}_{kj} (d k f j) Lovoo (f j i l)
     call dgemm('n','n',ov,o2,ov,-1.0E0_realk,t2f,ov,w2,ov,0.0E0_realk,w3,ov)
     ! part2
     ! sort gooov (jkif) -[1,4,2,3]> (jfki)
     call array_reorder_4d(1.0E0_realk,gooov,no,no,no,nv,[1,4,2,3],0.0E0_realk,w2)
     ! sort t^{df}_{jl} (d j f l) -[1,4,2,3]> (dljf)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,4,2,3],0.0E0_realk,w1)
     ! \sum_{jf} t^{df}_{jl} gooov_{jkif}(jfki)
     call dgemm('n','n',ov,o2,ov,1.0E0_realk,w1,ov,w2,ov,0.0E0_realk,w4,ov)
     ! sort result w4(d l k i) -[1,3,4,2]+> w3 (d k i l)
     call array_reorder_4d(1.0E0_realk,w4,nv,no,no,no,[1,3,4,2],1.0E0_realk,w3)
     ! part3
     ! govoo from gooov through 3,4,1,2
     call array_reorder_4d(1.0E0_realk,gooov,no,no,no,nv,[3,4,1,2],0.0E0_realk,w2)
     !\sum_{jf} t^{df}_{jk}(d k j f)in w1 govoo_{jfil}(jfil)
     call dgemm('n','n',ov,o2,ov,1.0E0_realk,w1,ov,w2,ov,1.0E0_realk,w3,ov)

     ! sort w3 (d k i l) -[1,2,4,3]> w2 (d k l i)
     call array_reorder_4d(1.0E0_realk,w3,nv,no,no,no,[1,2,4,3],0.0E0_realk,w2)
     ! sort \zeta^{da}_{kl} (d k a l} -[3,1,2,4]> w1(a d k l)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! \sum_{dkl} \zeta^{da}_{kl} w2(dkli)
     call dgemm('n','n',nv,no,o2*nv,1.0E0_realk,w1,nv,w2,o2*nv,0.0E0_realk,w4,nv)
     !write(msg,*)"rho g terms 1-3 (21H)"
     !call print_norm(w4,int(ov,kind=8),msg)
     !rho1 += w2(i a)
     call daxpy(ov,1.0E0_realk,w4,1,rho1,1)


     !rho e
     !-----
     ! part1
     ! \sum_{dke} \zeta^{de}_{ki}(d k e, i)^T g_{dkea} (d k e a)
     call dgemm('t','n',no,nv,v2*no,1.0E0_realk,m2,v2*no,gvovv,v2*no,0.0E0_realk,w2,no)
     !write(msg,*)"rho e - 1 (21DC)" 
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T
     call mat_transpose(no,nv,1.0E0_realk,w2,1.0E0_realk,rho1)

     if( DECinfo%PL > 2) then
        write(msg,*)"singles after loop" 
        call print_norm(rho1,int(ov,kind=8),msg)
     endif


     ! f-term
     ! part4
     ! sort t^{de}_{jl} (d j e l) -[1,3,2,4]> (d e j l)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,3,2,4],0.0E0_realk,w1)
     ! sort m^{de}_{ki} (d k e i) -[2,4,1,3]> (k i d e)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[2,4,1,3],0.0E0_realk,w2)
     ! \sum_{de} m2(k i d e)  t^{de}_{jl} (d e j l)
     call dgemm('n','n',o2,o2,v2,1.0E0_realk,w2,o2,w1,v2,0.0E0_realk,w3,o2)
     ! w3 (k i j l) -> w2 (i j l k)
     call array_reorder_4d(1.0E0_realk,w3,no,no,no,no,[2,3,4,1],0.0E0_realk,w2)
    
     ! sort gooov(jkla) -[1,3,2,4]> (j l k a)
     call array_reorder_4d(1.0E0_realk,gooov,no,no,no,nv,[1,3,2,4],0.0E0_realk,w1)
     ! \sum_{jl} w3(i j l k) gooov_{jkla} (j l k a)
     call dgemm('n','n',no,nv,o2*no,1.0E0_realk,w2,no,w1,no*o2,0.0E0_realk,w3,no)
     !write(msg,*)"rho f 4(LT21G)"
     !call print_norm(w3,int(ov,kind=8),msg)
     call mat_transpose(no,nv,1.0E0_realk,w3,1.0E0_realk,rho1)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho 4.f(21G):" 
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     !rho g
     !-----
     ! part 4
     ! sort t^{fe}_{kl} (f k e l) -[2,4,1,3]> (k l f e)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[2,4,1,3],0.0E0_realk,w1)
     ! sort govvv_{iedf} -[3,4,1,2]> gvvov_{dfie} (d f i e) -[2,4,1,3]> (f e d i)
     ! : [4,2,3,1]
     call array_reorder_4d(1.0E0_realk,govvv,no,nv,nv,nv,[4,2,3,1],0.0E0_realk,w2)
     ! \sum_{ef} t^{fe}_{kl}(k l f e)  gvvov_{dfie} (f e d i)
     call dgemm('n','n',o2,ov,v2,1.0E0_realk,w1,o2,w2,v2,0.0E0_realk,w4,o2)
     ! sort result w4(k l d i) -[3,1,2,4]> w3 (d k l i)
     call array_reorder_4d(1.0E0_realk,w4,no,no,nv,no,[3,1,2,4],0.0E0_realk,w3)
     ! sort \zeta^{da}_{kl} (d k a l} -[3,1,2,4]> w1(a d k l)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! \sum_{dkl} \zeta^{da}_{kl} (a d k l) w3(dkli)
     call dgemm('n','n',nv,no,o2*nv,1.0E0_realk,w1,nv,w3,o2*nv,0.0E0_realk,w4,nv)

     !rho e
     !-----
     ! part2
     ! sort: \zeta^{da}_{kl} (d k a l) -[3,1,2,4]> w1 (a, d k l)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! sort: gvooo_{dkil} -[1,2,4,3]> (dkli)
     call array_reorder_4d(1.0E0_realk,gvooo,nv,no,no,no,[1,2,4,3],0.0E0_realk,w3)
     ! \sum_{dke} \zeta^{da}_{kl}(a, dkl) g_{dkil} (d k l i)
     call dgemm('n','n',nv,no,o2*nv,1.0E0_realk,w1,nv,w3,o2*nv,0.0E0_realk,w2,nv)
     !w4 += w2(i a)^T
     call daxpy(ov,1.0E0_realk,w2,1,w4,1)
     !write(msg,*)"4.g+2.e 21BF:"
     !call print_norm(w4,int(ov,kind=8),msg)
     call daxpy(ov,-1.0E0_realk,w4,1,rho1,1)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho 4.g +2.e(21BF):" 
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     !rho a
     !-----
     ! sort \hat{L}_{bkia} (b k i a) -> w1
     call dcopy(o2v2,Lvoov,1,w1,1)
     ! sort u2^{bd}_{kl} (b k d l) -[1,2,4,3]> (b k l d)
     call array_reorder_4d(2.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w2)
     call array_reorder_4d(-1.0E0_realk,t2f,nv,no,nv,no,[1,4,2,3],1.0E0_realk,w2)
     ! w1 : \sum_{dl} u^{bd}_{kl}(b k l d) \hat{L}_{ldia} (l d i a) + \hat{L}_{bkia} (b k i a)
     call dgemm('n','n',ov,ov,ov,1.0E0_realk,w2,ov,Lovov,ov,1.0E0_realk,w1,ov)
     ! w2 : \sum_{bk}  w1(b k , i a)^T \zeta_k^b (b k)
     call dgemv('t',ov,ov,1.0E0_realk,w1,ov,m1,1,0.0E0_realk,w2,1)
     !write(msg,*)"rho a (11A)"
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T
     call mat_transpose(no,nv,1.0E0_realk,w2,1.0E0_realk,rho1)
     
     if( DECinfo%PL > 2) then
        write(msg,*)"rho a (11A):" 
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     !rho b
     !-----
     !part2
     ! w2 : sort t(d j b k) -[3,1,4,2]> t (b d k j)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[3,1,4,2],0.0E0_realk,w2)
     ! w3 : sort \hat{L}_{kbid} (k b i d} -[3,2,4,1]> \hat{L}_{kbid} (i b d k)
     call array_reorder_4d(1.0E0_realk,Lovov,no,nv,no,nv,[3,2,4,1],0.0E0_realk,w3)
     ! w1 : sort  F_{ij} -> w1
     call dcopy(o2,oof,1,w1,1)
     ! w1 : \sum_{bdk} \hat{L}_{kbid}(i b d k) t^{db}_{jk} (b d k j) + \hat{F}_{ij}
     call dgemm('n','n',no,no,no*v2,1.0E0_realk,w3,no,w2,no*v2,1.0E0_realk,w1,no)
     ! w2 : \sum_{j} \zeta_{j}^{a} (a j) w1(ij)^T = w2
     call dgemm('n','t',nv,no,no,1.0E0_realk,m1,nv,w1,no,0.0E0_realk,w2,nv) 
     !write(msg,*)"rho b - 2"
     !call print_norm(w2,int(ov,kind=8),msg)
     call daxpy(ov,-1.0E0_realk,w2,1,rho1,1)
     !part1
     ! w2 : sort t(d l b k) -[3,2,1,4]> t (b l d k)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[3,2,1,4],0.0E0_realk,w2)
     ! w1 : sort  F_{ba} -> w1
     call dcopy(v2,vvf,1,w1,1)
     ! w1 : \sum_{bdk} t^{db}_{lk} (b l d k)\hat{L}_{ldka}(l d k a)  + \hat{F}_{ba}
     call dgemm('n','n',nv,nv,o2*nv,-1.0E0_realk,w2,nv,Lovov,o2*nv,1.0E0_realk,w1,nv)
     ! w2 : \sum_{j} w1(b a)^T \zeta_{j}^{b} (b j) = w2
     call dgemm('t','n',nv,no,nv,1.0E0_realk,w1,nv,m1,nv,0.0E0_realk,w2,nv) 
     !write(msg,*)"rho b - 1"
     !call print_norm(w2,int(ov,kind=8),msg)
     call daxpy(ov,1.0E0_realk,w2,1,rho1,1)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho b (11B):"
        call print_norm(rho1,int(ov,kind=8),msg)
     endif
     
     !rho c - part2
     !-------------
     ! w3 : \sum_{dke} t^{de}_{kl} (d k e ,l)^T \zeta^{de}_{kj} (d k e j)
     call dgemm('t','n',no,no,no*v2,1.0E0_realk,t2f,v2*no,m2,v2*no,0.0E0_realk,w3,no)
     ! w1 : \sum_{lj} w3 (lj) L_{l j i a} (l j i a)
     call array_reorder_4d(2.0E0_realk,gooov,no,no,no,nv,[1,2,3,4],0.0E0_realk,w1)
     call array_reorder_4d(-1.0E0_realk,gooov,no,no,no,nv,[3,2,1,4],1.0E0_realk,w1)
     call dgemv('t',o2,ov,1.0E0_realk,w1,o2,w3,1,0.0E0_realk,w2,1)
     !rho1 += w2(i a)^T
     call mat_transpose(no,nv,-1.0E0_realk,w2,1.0E0_realk,rho1)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho c - 2(LT21B)"
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     !rho d
     !-----
     ! part1
     ! sort \zeta^{da}_{kl} (d k a l) -[3 1 2 4]> w1(a d k l)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! sort t^{de}_{kl} (d k e l) -[1,2,4,3]> w2(d k l e)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w2)
     ! w3 : \sum_{dkl} w1 (a d k l) w2 (d k l e)
     call dgemm('n','n',nv,nv,o2*nv,1.0E0_realk,w1,nv,w2,nv*o2,0.0E0_realk,w3,nv)
     ! w1 : \sum_{e} w3 (ae) F_{ie} (i e)^T
     call dgemm('n','t',nv,no,nv,1.0E0_realk,w3,nv,ovf,no,0.0E0_realk,w2,nv)
     !write(msg,*)"rho d - 1"
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T
     call daxpy(ov,-1.0E0_realk,w2,1,rho1,1)
     ! part2
     ! w3 : \sum_{dke} zeta^{de}_{ki} (d k e ,i)^T t^{de}_{kl} (d k e ,l)
     call dgemm('t','n',no,no,no*v2,1.0E0_realk,m2,v2*no,t2f,v2*no,0.0E0_realk,w3,no)
     ! w1 : \sum_{e} w3 (il) F_{la} (l a)
     call dgemm('n','n',no,nv,no,1.0E0_realk,w3,no,ovf,no,0.0E0_realk,w2,no)
     !write(msg,*)"rho d - 2"
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T
     call mat_transpose(no,nv,-1.0E0_realk,w2,1.0E0_realk,rho1)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho d(LT2EFM)"
        call print_norm(rho1,int(ov,kind=8),msg)

        write(msg,*)"singles end" 
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     call mem_dealloc(w4)
     
   
     !DOUBLES EXPRESSIONS
     !*******************

     !rho B
     !-----
     rho2 = 0.0E0_realk
     !gvvvv
     call array4_read(gao)
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yv,nv,xv,nv,yv,nv,w2)
     ! sort gvvvv_{cadb} (cadb) -> (cd ab)
     call array_reorder_4d(1.0E0_realk,gao%val,nv,nv,nv,nv,[1,3,2,4],0.0E0_realk,w1)
     call array4_dealloc(gao)
     ! sort \zeta^{cd}_{ij} (c i d j) -> (i j c d)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[2,4,1,3],0.0E0_realk,w2)
     ! \sum_{cd} \zeta{cd}_{ij} (ij cd) gvvvv_{cadb} (cd ab)
     call dgemm('n','n',o2,v2,v2,1.0E0_realk,w2,o2,w1,v2,0.0E0_realk,w3,o2)
     !write(msg,*)"rho B"
     !call print_norm(w3,int(o2v2,kind=8),msg)
     ! order w3(i j a b) and add to residual rho2 (a i b j)
     call array_reorder_4d(0.5E0_realk,w3,no,no,nv,nv,[3,1,4,2],1.0E0_realk,rho2)

     !rho 2. H part 
     ! \sum_{c} \zeta_{j}^{c} (c,j)^T Lvvov_{cbia} (cbia)
     call array_reorder_4d(2.0E0_realk,govvv,no,nv,nv,nv,[3,4,1,2],0.0E0_realk,w2)
     call array_reorder_4d(-1.0E0_realk,govvv,no,nv,nv,nv,[3,2,1,4],1.0E0_realk,w2)
     call dgemm('t','n',no,v2*no,nv,1.0E0_realk,m1,nv,w2,nv,0.0E0_realk,w1,no)
     !write(msg,*)"rho H - 2"
     !call print_norm(w1,int(o2v2,kind=8),msg)
     ! order w1(j b i a) and add to residual rho2 (a i b j)
     call array_reorder_4d(1.0E0_realk,w1,no,nv,no,nv,[4,3,2,1],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"H(B-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif

     !rho E
     !----- 
     ! sort t^{cd}_{mn} (c m d n) -[1,3,2,4]> (c d m n)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,3,2,4],0.0E0_realk,w1)
     ! sort govov_{manb} (m a n b) -[1,3,2,4]> (m n a b)
     call array_reorder_4d(1.0E0_realk,govov,no,nv,no,nv,[1,3,2,4],0.0E0_realk,w2)
     ! \sum_{mn} t^{cd}_{mn} (c d m n) govov_{manb} (m n a b)
     call dgemm('n','n',v2,v2,o2,1.0E0_realk,w1,v2,w2,o2,0.0E0_realk,w3,v2)
     ! sort \zeta^{cd}_{ij} (c i d j) -[2,4,1,3]> (i j c d)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[2,4,1,3],0.0E0_realk,w2)
     ! \sum_{cd} \zeta^{cd}_{ij} (i j c d) w3 (c d a b)
     call dgemm('n','n',o2,v2,v2,0.5E0_realk,w2,o2,w3,v2,0.0E0_realk,w1,o2)
     !write(msg,*)"rho E"
     !call print_norm(w1,int(o2v2,kind=8),msg)
     ! order w1(i j  a b) and add to residual rho2 (a i b j)
     call array_reorder_4d(1.0E0_realk,w1,no,no,nv,nv,[3,1,4,2],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"E(AM-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif


     !rho F
     !-----
     !part1
     ! \sum_{efm} \zeta^{ef}_{mj} (e m f , j)^T t^{ef}_{mn} (e m f , n)
     call dgemm('t','n',no,no,v2*no,1.0E0_realk,m2,v2*no,t2f,v2*no,0.0E0_realk,w1,no)
     ! sort Lovov{ianb} (i a n b) -[2,1,4,3]> (a i b n)
     call array_reorder_4d(1.0E0_realk,Lovov,no,nv,no,nv,[2,1,4,3],0.0E0_realk,w2)
     ! \sum_{n} Lovov_{ianb} (aibn) w1 (j n)^T
     call dgemm('n','t',v2*no,no,no,1.0E0_realk,w2,v2*no,w1,no,0.0E0_realk,w3,v2*no)
     !write(msg,*)"rho F - 1"
     !call print_norm(w3,int(o2v2,kind=8),msg)
     !rho2 += w1(a i b j)
     call daxpy(o2v2,-1.0E0_realk,w3,1,rho2,1)
     !part2
     ! sort \zeta^{ea}_{mn} (e m a n) -[3,1,2,4]> (a e m n)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! sort \t^{ef}_{mn} (e m f n) -[1,2,4,3]> (e m n f)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w2)
     ! \sum_{emn} \zeta^{ea}_{mn}(a e m n) t^{ef}_{mn} (e m n f)
     call dgemm('n','n',nv,nv,o2*nv,1.0E0_realk,w1,nv,w2,o2*nv,0.0E0_realk,w3,nv)
     ! sort Lovov{ifjb} (i f j b) -> (f i b j)
     call array_reorder_4d(1.0E0_realk,Lovov,no,nv,no,nv,[2,1,4,3],0.0E0_realk,w2)
     ! \sum_{f} w3 (a f) Lovov_{ifjb} (f i b j) 
     call dgemm('n','n',nv,o2*nv,nv,1.0E0_realk,w3,nv,w2,nv,0.0E0_realk,w1,nv)
     !write(msg,*)"rho F - 2"
     !call print_norm(w1,int(o2v2,kind=8),msg)
     !rho2 += w1(a i b j)
     call daxpy(o2v2,-1.0E0_realk,w1,1,rho2,1)

     if( DECinfo%PL > 2) then
        write (msg,*)"F(EM-TRM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif


     !rho G
     !-----
     ! part1
     ! copy vv F to w3
     call dcopy(v2,vvf,1,w3,1)
     ! sort t^{fe}_{nm} (f n e m) -[3,2,1,4]> (e n f m)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[3,2,1,4],0.0E0_realk,w1)
     ! vv F - \sum_{fmn} t^{fe}_{nm}(e n f m) Lovov_{nfmb} (n f m b)
     call dgemm('n','n',nv,nv,o2*nv,-1.0E0_realk,w1,nv,Lovov,o2*nv,1.0E0_realk,w3,nv)
     ! sort \zeta^{ae}_{ij} (a i e j) -> (a i j e)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w1)
     ! \sum_{e} \zeta^{ae}_{ij} (a i j e) w3 (e a)
     call dgemm('n','n',o2*nv,nv,nv,1.0E0_realk,w1,o2*nv,w3,nv,0.0E0_realk,w2,o2*nv)
     !write(msg,*)"rho G - 1"
     !call print_norm(w2,int(o2v2,kind=8),msg)
     ! order w1(a i j b) and add to residual rho2 (a i b j)
     call array_reorder_4d(1.0E0_realk,w2,nv,no,no,nv,[1,2,4,3],1.0E0_realk,rho2)
     !part2
     ! copy oo F to w3
     call dcopy(o2,oof,1,w3,1)
     ! sort t^{fe}_{nm} (f n e m) -[4,3,1,2]> (m e f n)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[4,3,1,2],0.0E0_realk,w2)
     ! sort Lovov_{mejf} (mejf} -[3,1,2,4]> (j m e f)
     call array_reorder_4d(1.0E0_realk,Lovov,no,nv,no,nv,[3,1,2,4],0.0E0_realk,w1)
     ! oo F + \sum_{efm} Lovov_{mejf} (j mef) t^{fe}_{nm}(mef n) 
     call dgemm('n','n',no,no,v2*no,1.0E0_realk,w1,no,w2,v2*no,1.0E0_realk,w3,no)
     ! \sum_{n} \zeta^{ab}_{in} (a i b n) w3 (j n)^T
     call dgemm('n','t',v2*no,no,no,1.0E0_realk,m2,v2*no,w3,no,0.0E0_realk,w2,v2*no)
     !write(msg,*)"rho G - 2"
     !call print_norm(w2,int(o2v2,kind=8),msg)
     ! order w1(a i b j) and add to residual rho2 (a i b j)
     call array_reorder_4d(-1.0E0_realk,w2,nv,no,nv,no,[1,2,3,4],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"G(22E-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif


     !rho A
     !-----
     ! copy goooo[j n i m ]  -[2,4,3,1]> w1[n m i j]
     !goooo
     call array4_read(gao)
     call successive_4ao_mo_trafo(nb,gao%val,xo,no,yo,no,xo,no,yo,no,w2)
     call array_reorder_4d(1.0E0_realk,gao%val,no,no,no,no,[2,4,3,1],0.0E0_realk,w1)
     call array4_dealloc(gao)
     ! sort govov{jfie}(jfie) -[4,2,3,1]> (e f i j)
     call array_reorder_4d(1.0E0_realk,govov,no,nv,no,nv,[4,2,3,1],0.0E0_realk,w2)
     ! sort t^{fe}_{nm} (f n e m) -[2,4,3,1]> (n m e f)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[2,4,3,1],0.0E0_realk,w3)
     !\sum_{ef} t^{fe}_{nm} (n m e f) govov_{jfie} (ef ij) + goooo
     call dgemm('n','n',o2,o2,v2,1.0E0_realk,w3,o2,w2,v2,1.0E0_realk,w1,o2)
     ! sort \zeta^{ab}_{mn} (a m b n) -[1,3,4,2]> (a b n m)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[1,3,4,2],0.0E0_realk,w2)
     !\frac{1}{2} \sum_{mn} w2:\zeta^{ab}_{mn} (a b n m} w1(n m i j)
     call dgemm('n','n',v2,o2,o2,0.5E0_realk,w2,v2,w1,o2,0.0E0_realk,w3,v2)
     !write(msg,*)"rho A"
     !call print_norm(w3,int(o2v2,kind=8),msg)
     ! order and add (a b i j) to residual rho2 (a i b j)
     call array_reorder_4d(1.0E0_realk,w3,nv,nv,no,no,[1,3,2,4],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"A(22A-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif

     !rho D
     !-----
     ! copy Lvoov in w3 
     call dcopy(o2v2,Lvoov,1,w3,1)
     ! sort u2^{ef}_{mn} (emfn) -> em nf
     call array_reorder_4d(2.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w1)
     call array_reorder_4d(-1.0E0_realk,t2f,nv,no,nv,no,[1,4,2,3],1.0E0_realk,w1)
     ! \sum_{fn} u2^{ef}_{mn} (emnf) Lovov_{nfjb} (nfjb) + Lvoov_{emjb} (emjb)
     call dgemm('n','n',ov,ov,ov,1.0E0_realk,w1,ov,Lovov,ov,1.0E0_realk,w3,ov)
     ! \sum_{em} \zeta^{ae}_{im} (a i e m) w3 (e m j b)
     call dgemm('n','n',ov,ov,ov,1.0E0_realk,m2,ov,w3,ov,0.0E0_realk,w2,ov)
     !write(msg,*)"rho D"
     !call print_norm(w2,int(o2v2,kind=8),msg)
     ! order w1(a i j b) and add to residual rho2 (a i b j)
     call array_reorder_4d(1.0E0_realk,w2,nv,no,no,nv,[1,2,4,3],1.0E0_realk,rho2)
     ! order w1(a j i b) and add to residual rho2 (a i b j)
     call array_reorder_4d(-0.5E0_realk,w2,nv,no,no,nv,[1,3,4,2],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"D(22D-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif

     !rho C
     !-----
     ! sort govov(nbif) -[4,1,2,3]> (fnbi)
     call array_reorder_4d(1.0E0_realk,govov,no,nv,no,nv,[4,1,2,3],0.0E0_realk,w1)
     ! sort t^{ef}_{nm} (enfm) -[1,4,3,2]> (emfn)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,4,3,2],0.0E0_realk,w2)
     ! sort goovv(imeb) -[3,2,4,1]> (e m b i)
     call array_reorder_4d(1.0E0_realk,goovv,no,no,nv,nv,[3,2,4,1],0.0E0_realk,w3)
     ! \sum_{fn} t^{ef)_{nm}(em fn) govov_{nbif} (fn bi) + goovv
     call dgemm('n','n',ov,ov,ov,-1.0E0_realk,w2,ov,w1,ov,1.0E0_realk,w3,ov)
     ! sort make anti-u2 analog of \zeta^{ae}_{mj} (a m e j) as 2 (a j e m) + (a m e j)
     call array_reorder_4d(2.0E0_realk,m2,nv,no,nv,no,[1,4,3,2],0.0E0_realk,w2)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[1,2,3,4],1.0E0_realk,w2)
     ! \sum_{em} w2(a j e m) w3(em bi)
     call dgemm('n','n',ov,ov,ov,1.0E0_realk,w2,ov,w3,ov,0.0E0_realk,w1,ov)
     !write(msg,*)"rho C"
     !call print_norm(w1,int(o2v2,kind=8),msg)
     ! order w1(a j b i) and add to residual rho2 (a i b j)
     call array_reorder_4d(-0.5E0_realk,w1,nv,no,nv,no,[1,4,3,2],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"C(22C-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif

     !rho H
     !-----
     ! 2*\zeta_i^a * F_{jb} - \zeta_j^b F_{ib}
     do j = 1, no
       do b = 1, nv
         do i = 1, no
           do a = 1 ,nv
             rho2(a,i,b,j) = rho2(a,i,b,j) + 2*m1(a,i)*ovf(j+(b-1)*no) - m1(a,j)*ovf(i+(b-1)*no)
           enddo
         enddo
       enddo
     enddo

     if( DECinfo%PL > 2) then
        write (msg,*)"H(A12-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif


     !\sum_{k} Lovoo_{jbik} \zeta_{k}^{a} (a k)^T
     call dgemm('n','t',o2*nv,nv,no,1.0E0_realk,Lovoo,o2*nv,m1,nv,0.0E0_realk,w1,o2*nv)
     ! order w1(j b i a) and add to residual rho2 (a i b j)
     call array_reorder_4d(-1.0E0_realk,w1,no,nv,no,nv,[4,3,2,1],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
       write (msg,*)"H(C12-TERM)"
       print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif


     !END OF CONTRIBUTIONS
     !********************


     !PERMUTE RHO 2
     call array_reorder_4d(1.0E0_realk,rho2,nv,no,nv,no,[1,2,3,4],0.0E0_realk,w1)
     call array_reorder_4d(1.0E0_realk,w1,nv,no,nv,no,[3,4,1,2],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho end"
        call print_norm(rho2,int(o2v2,kind=8),msg)
     endif

     !ADD RIGHT HAND SIDES
     call array_reorder_4d(2.0E0_realk,Lovov,no,nv,no,nv,[2,1,4,3],1.0E0_realk,rho2)
     call mat_transpose(no,nv,2.0E0_realk,ovf,1.0E0_realk,rho1)

     call mem_dealloc(oof)
     call mem_dealloc(ovf)
     call mem_dealloc(vof)
     call mem_dealloc(vvf)

     call mem_dealloc(govov)
     call mem_dealloc(goovv)
     call mem_dealloc(Lovov)
     call mem_dealloc(Lvoov)


     call mem_dealloc(gvovv)
     call mem_dealloc(govvv)

     call mem_dealloc(Lovoo)
     call mem_dealloc(gvooo)
     call mem_dealloc(gooov)

     call mem_dealloc(w1)
     call mem_dealloc(w2)
     call mem_dealloc(w3)

    end subroutine get_ccsd_multipliers_simple

    !This function calculates the norm according to how it should be compared to
    !DALTON, because they construct the symmetrization for each term
    !individually and store the packed residual, furthermore they only print the
    !square of the 2 norm
    function symmetrized_packed_1norm(m,nv,no)result(nrm)
      implicit none
      integer,intent(in)     :: nv,no
      real(realk),intent(in) :: m(nv,no,nv,no)
      integer                :: a,i,b,j
      real(realk)            :: nrm

      nrm = 0.0E0_realk
      do b = 1, nv
        do j = 1, no
          do a = 1,nv
            do i = 1, no
              if(a+(i-1)*nv<=b+(j-1)*nv)then
                nrm = nrm + (m(a,i,b,j) + m(b,j,a,i) )**2
              endif
            enddo
          enddo
        enddo
      enddo

    end function symmetrized_packed_1norm
     
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

    !> \brief should be a general subroutine to get the guess amplitudes when
    !starting up a CCSD or CC2 calculation and checks for files which contain
    !amplitudes of a previous calculation. If none are found the usual zero guess
    !is returned
    !> \author Patrick Ettenhuber
    !> \date December 2012
   subroutine get_guess_vectors_simple(mylsitem,t2,t1,t1_final,gao,&
   &get_mult,Co,Co2,Cv2,no,nv,nb,xocc,yvirt,restart)
     implicit none
     type(lsitem),intent(inout) :: mylsitem
     !> contains the guess doubles amplitudes on output
     type(array4),intent(inout) :: t2
     !> contains the singles amplitudes on output
     type(array2),intent(inout) :: t1
     type(array2),intent(inout) :: t1_final,Co,Co2,Cv2,xocc,yvirt
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
     character(ARR_MSG_LEN) :: msg
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
          call Get_AOt1Fock(mylsitem,t1_final,fockguess,no,nv,nb,Co,Co2,Cv2)
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
   subroutine save_current_guess_simple(iter,t2,t1,get_mult)
    implicit none
    !> iteration number
    integer,intent(in) :: iter
    !> doubles guess amplitudes for the next iteration
    type(array4), intent(in) :: t2
    !> singles guess amplitudes for the next iteration
    type(array2), intent(in) :: t1
    logical, intent(in) :: get_mult
    !> alternating filenames for the doubles amplitudes
    character(3) :: safefilet21,safefilet22
    !> alternating filenames for the singles amplitudes
    character(3) :: safefilet11,safefilet12
    integer :: fu_t21,fu_t22
    integer :: fu_t11,fu_t12
    logical :: file_status11,file_status12,file_status21,file_status22
    logical :: all_singles
    character(ARR_MSG_LEN) :: msg
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

  ! begin pablo
  !> Purpose: calculate AO int. in batches and transform them to
  !           full MO basis (non T1-transformed)
  !           Based on Patrick's CCSD routine: get_ccsd_residual_integral_driven
  !
  !> Author:  Pablo Baudin
  !> Date:    October 2013
  !
  subroutine get_packed_gmo(MyLsItem,cmo_occ,cmo_vir,pack_gmo, &
             & nbas,nocc,nvir,PRbatchInfo)

    implicit none

    integer, intent(in) :: nbas, nocc, nvir
    real(realk), pointer, intent(out) :: pack_gmo(:)
    real(realk), intent(in) :: cmo_occ(nbas,nocc), cmo_vir(nbas,nvir)
    !> dim packed :
    integer :: i,j,p,q,r,s,pq,rs, job_index, full_ind, pack_ind
    type(array2) :: cmo_full
    type(array4) :: pqrs
    !
    ! Elementary types needed for the calculation
    real(realk), pointer :: gao(:), gmo(:)
    integer(kind=long) :: gaosize, gmosize, pack_gmosize
    integer :: iter, niter
    logical :: first_round, dynamic_load, print_debug, local
    !
    ! CHECKING and MEASURING variables
    integer(kind=long) :: maxsize64
    real(realk) :: MemFree, tcpu, twall, tcpu_start, twall_start, tcpu_end, twall_end
    real(realk) :: time_start, timewall_start 
    integer     :: scheme
    integer(kind=long) :: els2add
    !   
    ! variables used for MO BATCH
    integer :: nRbatch, nPbatch, PR_batch, dimR, dimP
    integer :: max_mem, min_mem, PR_sta, PR_end, P_sta, P_end, R_sta, R_end
    type(MObatchInfo), intent(out) :: PRbatchInfo
    ! 
    ! variables used for AO BATCH construction and INTEGRAL calculation
    integer :: alphaB, gammaB, dimAlpha, dimGamma
    integer :: dim1, dim2, dim3, K, MinAObatch
    integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd
    integer :: iorb, idx
    type(batchtoorb), pointer :: batch2orbAlpha(:)
    type(batchtoorb), pointer :: batch2orbGamma(:)
    Character            :: INTSPEC(5)
    logical :: FoundInMem, fullRHS, doscreen
    integer :: nbatches
    integer :: MaxAllowedDimAlpha, MaxActualDimAlpha, nbatchesAlpha
    integer :: MaxAllowedDimGamma, MaxActualDimGamma, nbatchesGamma
    integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), &
                      & batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchdimGamma(:), &
                      & batchsizeGamma(:), batchindexGamma(:)
    !
    ! MPI variables:
    logical :: master
    integer(kind=ls_mpik) :: mynum, nnod, ierr
    !
    ! Screening integrals stuff:
    type(DECscreenITEM) :: DecScreen

    !> LS item with information needed for integrals
    type(lsitem), intent(inout) :: MyLsItem


    ! Set default values for the path throug the routine
    ! **************************************************
    scheme           = 0
    dynamic_load     = DECinfo%dyn_load
    print_debug      = (DECinfo%PL>2)
    iter             = 1
    local            = .false.

    !> info about the scheme for packing MO integrals:
    !
    !    MOscheme = 0 => no packing of integrals;       mem = N**4
    !    MOscheme = 1 => packing using pq>=rs;          mem = N**2*(N**2+1)/2
    !    MOscheme = 2 => packing using p>=q and r>=s;   mem = [N*(N+1)]**2/4
    !    MOscheme = 3 => combination of scheme 1 and 2; mem = N*(N+1)*[N*(N+1) + 1]/8
    !

    ! Set integral info
    ! *****************
    INTSPEC(1)       = 'R' !R = Regular Basis set on the 1th center 
    INTSPEC(2)       = 'R' !R = Regular Basis set on the 2th center 
    INTSPEC(3)       = 'R' !R = Regular Basis set on the 3th center 
    INTSPEC(4)       = 'R' !R = Regular Basis set on the 4th center 
    INTSPEC(5)       = 'C' !C = Coulomb operator

    doscreen = MyLsItem%setting%scheme%cs_screen.OR. &
             & MyLsItem%setting%scheme%ps_screen

    ! Set MPI related info
    ! ********************
    master           = .true.
    mynum            = int(0,kind=ls_mpik)
    nnod             = 1
#ifdef VAR_MPI
    mynum            = infpar%lg_mynum
    master           = (infpar%lg_mynum == 0)
    nnod             = infpar%lg_nodtot
#endif

    ! Some timings
    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu_start,twall_start,DECinfo%output)
    call LSTIMER('START',time_start,timewall_start,DECinfo%output)

    ! Initialize stuff
    nullify(orb2batchAlpha)
    nullify(batchdimAlpha)
    nullify(batchsizeAlpha)
    nullify(batch2orbAlpha)
    nullify(batchindexAlpha)
    nullify(orb2batchGamma)
    nullify(batchdimGamma)
    nullify(batchsizeGamma)
    nullify(batch2orbGamma)
    nullify(batchindexGamma)
    nullify(gao)
    nullify(gmo)
    nullify(pack_gmo)
    nullify(PRbatchInfo%Pdims)
    nullify(PRbatchInfo%Rdims)
    nullify(PRbatchInfo%PStarts)
    nullify(PRbatchInfo%RStarts)
    nullify(PRbatchInfo%PR_index)
    nullify(PRbatchInfo%PR_packInd)
    

    cmo_full%dims = [nbas, nbas]
    call mem_alloc(cmo_full%val,nbas,nbas)

    ! get full MO coeficients:
    cmo_full%val(:,:nocc)       = cmo_occ
    cmo_full%val(:,nocc+1:nbas) = cmo_vir
    !
    !
    !==================================================
    !                  Batch construction             !
    !==================================================


    if (master) then
      ! Get free memory and determine maximum batch sizes
      ! -------------------------------------------------
      call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch,'R')
      call get_currently_available_memory(MemFree)
      call get_max_batch_sizes(scheme,nbas,nvir,nocc,MaxAllowedDimAlpha, &
           & MaxAllowedDimGamma,MinAObatch,DECinfo%manual_batchsizes,iter, &
           & MemFree,.true.,els2add,local)
    end if

    ! MPI: here you should start the slaves!!

    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************

    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchGamma,nbas)
    call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma, &
         & nbas,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma, &
         & nbatchesGamma,orb2BatchGamma,'R')

    if (master.and.DECinfo%PL>1) write(DECinfo%output,*) & 
               & 'BATCH: Number of Gamma batches   = ', nbatchesGamma, &
               & 'with maximum size', MaxActualDimGamma 


    ! Translate batchindex to orbital index
    ! -------------------------------------
    call mem_alloc(batch2orbGamma,nbatchesGamma)
    do idx=1,nbatchesGamma
       call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx))
       batch2orbGamma(idx)%orbindex = 0
       batch2orbGamma(idx)%norbindex = 0
    end do
    do iorb=1,nbas
       idx = orb2batchGamma(iorb)
       batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
       K = batch2orbGamma(idx)%norbindex
       batch2orbGamma(idx)%orbindex(K) = iorb
    end do


    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************

    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchAlpha,nbas)
    call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha, &
         & nbas,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha, &
         & nbatchesAlpha,orb2BatchAlpha,'R')

    if (master.and.DECinfo%PL>1) write(DECinfo%output,*) & 
        & 'BATCH: Number of Alpha batches   = ', nbatchesAlpha, &
        & 'with maximum size',MaxActualDimAlpha


    ! Translate batchindex to orbital index
    ! -------------------------------------
    call mem_alloc(batch2orbAlpha,nbatchesAlpha)
    do idx=1,nbatchesAlpha
       call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
       batch2orbAlpha(idx)%orbindex = 0
       batch2orbAlpha(idx)%norbindex = 0
    end do
    do iorb=1,nbas
       idx = orb2batchAlpha(iorb)
       batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
       K = batch2orbAlpha(idx)%norbindex
       batch2orbAlpha(idx)%orbindex(K) = iorb
    end do



    ! ****************************************************
    ! *  Allocate matrices and determine MO batch info.  *
    ! ****************************************************

    ! get minimum and maximum memory requirements for MO batches:
    !
    !      tmp1 [P_batch, beta, gammaB, delta]
    !      tmp2 [P_batch, beta, gammaB, s]
    !      tmp3 [R_batch, s, P_batch, beta]
    !      gmo  [R_batch, s, P_batch, r] 
    !
    min_mem = 2*nbas*nbas*(MaxActualdimGamma + 1)
    max_mem = 2*nbas*nbas*nbas*(MaxActualdimGamma + nbas)

    if (nnod>1) then  
      ! Prepare MO transfo. of integrals:
      call get_MO_batches_size(min_mem, max_mem, nbas, nPbatch, dimP, &
           & 2*nbas*nbas, 2*nbas*nbas*MaxActualdimGamma)
      gmosize = int(i8*nbas*nbas*dimR*dimP,kind=long)
    else
      ! then the whole integrals have to be stored in core:
      !
      !      gmo  [p, q, r, s] nbas**4
      !
      gmosize = int(i8*nbas*nbas*nbas*nbas,kind=long)
      min_mem = nbas*nbas*(2*MaxActualdimGamma + nbas*nbas + 1) 
      call get_MO_batches_size(min_mem, max_mem, nbas, nPbatch, dimP, &
           & 2*nbas*nbas, 2*nbas*nbas*MaxActualdimGamma)
    end if 

    dimR = dimP
    nRbatch = nPbatch
    ! DEBUG HACK !!!!
    if (DECInfo%cc_driver_debug) then
      dimP = 6
      dimR = 6
      nPbatch = (nbas-1)/dimP + 1
      nRbatch = (nbas-1)/dimR + 1
      print *, 'debug test dimP, dimR, nPbatch, nRbatch', &
               & dimP, dimR, nPbatch, nRbatch
    end if
    ! DEBUG HACK !!!!

    PRbatchInfo%nRbatch = nRbatch
    PRbatchInfo%nPbatch = nPbatch
    call get_MO_batches_info(PRbatchInfo, dimP, dimR, nbas)

    gaosize = int(i8*nbas*nbas*MaxActualDimAlpha*MaxActualDimGamma,kind=long)
    call mem_alloc(gmo,gmosize)
    call mem_alloc(gao,gaosize)
    gmo = 0.0E0_realk

    ! Sanity checks for matrix sizes which need to be filled:
    if (gmosize>MaxInt.or.gaosize>MaxInt) then
    call lsquit("ERROR(CCSD):matrix sizes too large, &
        & please recompile with 64bit integers",-1)
    endif


    ! *******************************************************
    ! *  This subroutine builds the full screening matrix.
    call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting, &
         & nbatches,nbatchesAlpha,nbatchesGamma,INTSPEC)

    if (mylsitem%setting%scheme%cs_screen .or. &
         & mylsitem%setting%scheme%ps_screen) then

       call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
            & nbas,nbatchesAlpha,nbatchesGamma,&
            & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)

       call II_getBatchOrbitalScreenK(DecScreen,mylsitem%setting,&
            & nbas,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
            & batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
    end if
    ! *******************************************************

    if (master) call lstimer('CCSD part A',time_start,timewall_start,DECinfo%output)

    fullRHS = (nbatchesGamma.eq.1).and.(nbatchesAlpha.eq.1)


    !**********************************
    ! Begin the loop over gamma batches
    !**********************************

    first_round = .false.
    if (dynamic_load) first_round = .true.

    BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches

       dimGamma   = batchdimGamma(gammaB)                         ! Dimension of gamma batch
       GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
       GammaEnd   = batch2orbGamma(gammaB)%orbindex(dimGamma)     ! Last index in gamma batch

    !**********************************
    ! Begin the loop over alpha batches
    !**********************************

    BatchAlpha: do alphaB = 1,nbatchesAlpha    ! AO batches
      
#ifdef VAR_MPI
       !! check if this is my job:
       !job_index = alphaB + (gammaB-1)*nbatchesAlpha
       !job_index = mod(job_index, nnod)
       !if (job_index .ne. mynum) cycle ! this is not my job
       !print *, "mynum",mynum
#endif

       dimAlpha   = batchdimAlpha(alphaB)                        ! Dimension of alpha batch
       AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)           ! First index in alpha batch
       AlphaEnd   = batch2orbAlpha(alphaB)%orbindex(dimAlpha)    ! Last index in alpha batch

       
       if (doscreen)Mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGabKLHS(alphaB)%p
       if (doscreen)Mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGabKRHS(gammaB)%p
       !
       ! Get two-electron int. in exchange form: (alphaB, beta, gammaB, delta):
       call II_GET_DECPACKED4CENTER_K_ERI(DECinfo%output,DECinfo%output, &
            & Mylsitem%setting,gao,batchindexAlpha(alphaB),batchindexGamma(gammaB), &
            & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),dimAlpha,nbas,dimGamma, &
            & nbas,nbatches,INTSPEC,fullRHS)
       call lsmpi_poke()

       ! Loop over MO batches:
       BatchPR: do PR_batch = 1, nRbatch*nPbatch

         PR_sta = PRbatchInfo%PR_index(PR_batch)
         PR_end = PRbatchInfo%PR_index(PR_batch+1) - 1
         P_sta  = PRbatchInfo%PStarts(PR_batch)
         dimP   = PRbatchInfo%Pdims(PR_batch)
         R_sta  = PRbatchInfo%RStarts(PR_batch)
         dimR   = PRbatchInfo%Rdims(PR_batch)

         call gao_to_gmo(gmo(PR_sta:PR_end), gao, cmo_full%val, nbas, &
              & AlphaStart, dimAlpha, GammaStart, dimGamma, P_sta, dimP, R_sta, dimR)

         ! GMO INT. MAY BE PACKED HERE BEFORE SUM OVER AO BATCHES 

       end do BatchPR


    end do BatchAlpha
    end do BatchGamma


    ! Free gamma stuff
    call mem_dealloc(orb2batchGamma)
    call mem_dealloc(batchdimGamma)
    call mem_dealloc(batchsizeGamma)
    call mem_dealloc(batchindexGamma)
    do idx=1,nbatchesGamma
       call mem_dealloc(batch2orbGamma(idx)%orbindex)
       batch2orbGamma(idx)%orbindex => null()
    end do
    call mem_dealloc(batch2orbGamma)

    ! Free alpha stuff
    call mem_dealloc(orb2batchAlpha)
    call mem_dealloc(batchdimAlpha)
    call mem_dealloc(batchsizeAlpha)
    call mem_dealloc(batchindexAlpha)
    do idx=1,nbatchesAlpha
       call mem_dealloc(batch2orbAlpha(idx)%orbindex)
       batch2orbAlpha(idx)%orbindex => null()
    end do
    call mem_dealloc(batch2orbAlpha)

    pack_gmosize = int(i8*nbas*nbas*nbas*nbas, kind=long) 
    call mem_alloc(pack_gmo,pack_gmosize)

    ! MO INT. SHOULD BE PACKED HERE (at least):
    call dcopy(nbas*nbas*nbas*nbas,gmo,1,pack_gmo,1)

    ! Free integrals
    call mem_dealloc(gao)
    call mem_dealloc(gmo)

  end subroutine get_packed_gmo


  !> Purpose: Transform AO int. into MO int. in batches
  !           
  !> Author:  Pablo Baudin
  !> Date:    October 2013
  subroutine gao_to_gmo(gmo, gao, cmo_full, nbas, &
            & AlphaStart, dimAlpha, GammaStart, dimGamma, &
            & PStart, dimP, RStart, dimR)

    implicit none

    real(realk), intent(inout) :: gmo(:)
    real(realk), intent(in) :: gao(dimAlpha*nbas*dimGamma*nbas), cmo_full(nbas,nbas)
    integer, intent(in) :: nbas, PStart, dimP, RStart, dimR
    integer, intent(in) :: AlphaStart, dimAlpha, GammaStart, dimGamma

    integer :: AlphaEnd, GammaEnd, Pend, Rend
    real(realk), pointer, dimension(:)   :: tmp1, tmp2, tmp3     => null()
    real(realk), pointer, dimension(:,:) :: cmo_alpha, cmo_beta  => null()
    real(realk), pointer, dimension(:,:) :: cmo_gamma, cmo_delta => null()
   

    AlphaEnd = AlphaStart+dimAlpha-1
    GammaEnd = GammaStart+dimAlpha-1
    Pend     = PStart+dimP-1
    Rend     = RStart+dimR-1
 
    ! allocation stuff:
    call mem_alloc(cmo_alpha,dimAlpha,dimP)
    call mem_alloc(cmo_beta,nbas,nbas)
    call mem_alloc(cmo_gamma,dimGamma,dimR)
    call mem_alloc(cmo_delta,nbas,nbas)

    call mem_alloc(tmp1, int(i8*nbas*nbas*dimP*dimGamma,kind=long))
    call mem_alloc(tmp2, int(i8*nbas*nbas*dimP*dimGamma,kind=long))
    call mem_alloc(tmp3, int(i8*nbas*nbas*dimP*dimR,kind=long))

    ! initialisation of transfo. matrices:
    cmo_alpha = cmo_full(AlphaStart:AlphaEnd,PStart:Pend)
    cmo_beta  = cmo_full
    cmo_gamma = cmo_full(GammaStart:GammaEnd,RStart:Rend)
    cmo_delta = cmo_full


    ! transfo 1st index => [P_batch, beta, gammaB, delta]
    call dgemm('t','n',dimP,nbas*dimGamma*nbas,dimAlpha,1.0E0_realk, &
         & cmo_alpha,dimAlpha,gao,dimAlpha,0.0E0_realk,tmp1,dimP)
    call lsmpi_poke() 

    ! transfo last index => [P_batch, beta, gammaB, s]
    call dgemm('n','n',dimP*nbas*dimGamma,nbas,nbas,1.0E0_realk, &
         & tmp1,dimP*nbas*dimGamma,cmo_delta,nbas,0.0E0_realk, &
         & tmp2,dimP*nbas*dimGamma)
    call lsmpi_poke() 

    ! reorder array => [gammaB, s, P_batch, beta]
    call array_reorder_4d(1.0E0_realk,tmp2,dimP,nbas,dimGamma,nbas, &
         & [3,4,1,2],0.0E0_realk,tmp1)
    call lsmpi_poke() 

    ! transfo 1st index => [R_batch, s, P_batch, beta]
    call dgemm('t','n',dimR,nbas*dimP*nbas,dimGamma,1.0E0_realk, &
         & cmo_gamma,dimGamma,tmp1,dimGamma,0.0E0_realk,tmp3,dimR) 
    call lsmpi_poke() 
     
    ! transfo last index => [R_batch, s, P_batch, q]
    call dgemm('n','n',dimR*nbas*dimP,nbas,nbas,1.0E0_realk, &                               
         & tmp3,dimR*nbas*dimP,cmo_beta,nbas,0.0E0_realk, &
         & tmp1,dimR*nbas*dimP)
    call lsmpi_poke() 

    ! reorder [R_batch, s, P_batch, q] => [P_batch, q, R_batch, s]
    ! sum over previous integrals: 
    call mat_transpose(nbas*dimR,nbas*dimP,1.0E0_realk,tmp1,1.0E0_realk,gmo)
    call lsmpi_poke() 


    ! free array
    call mem_dealloc(tmp1)
    call mem_dealloc(tmp2)
    call mem_dealloc(tmp3)

    call mem_dealloc(cmo_alpha)
    call mem_dealloc(cmo_beta)
    call mem_dealloc(cmo_gamma)
    call mem_dealloc(cmo_delta)

  end subroutine gao_to_gmo


  !> Purpose: Provide a balanced number of MO batches for two
  !           equivalent nested loops.
  !
  !  (IN):    max_mem: maximum memory requirement (no batches)
  !  (IN):    min_mem: minimum memory requirement (dimBatch = 1)
  !  (IN):    tot_size: total size of the loops (e.g. nbas, nvir or nocc)
  !
  !  (OUT):   dimBatch is calculated as solution of:
  !
  !           const_A*dimBatch**2 + const_B*dimBatch <= free_mem
  !
  !> Author:  Pablo Baudin
  !> Date:    October 2013
  subroutine get_MO_batches_size(min_mem, max_mem, tot_size, nBatch, &
             & dimBatch, const_A, const_B)

    implicit none
  
    integer, intent(in) :: min_mem, max_mem, tot_size
    integer, intent(in) :: const_A, const_B
    integer, intent(inout) :: nBatch, dimBatch
     
    integer :: free_mem
    real(realk) :: MemFree, deter


    call get_currently_available_memory(MemFree)
    free_mem = int(MemFree*1024**3/realk)

    ! Check than min_mem is available:
    if ((free_mem-min_mem) < 0) stop &
       & 'Insufficient memory in MO batches construction'

    ! Check that batching is needed:
    if ((free_mem-max_mem) >= 0) then
      nBatch = 1
      dimBatch = tot_size
 
    else
     ! solve: const_A*dimBatch**2 + const_B*dimBatch <= free_mem
     deter = const_B*const_B + 4.0E0_realk*const_A*free_mem

     ! check that constants are all positive:
     if (deter<0.0E0_realk) stop 'Constants A and B have to be &
        & positive in MO batches construction' 
 
     dimBatch = int((-const_B + sqrt(deter))/(const_A*2.0E0_realk)) 
     nBatch   = (tot_size-1)/dimBatch + 1
    end if

  end subroutine get_MO_batches_size


  !> Purpose: Get information regarding MO batches for two 
  !           equivalent nested loops:
  !
  !           PR_index: Starting storage position of PR batch
  !                     in the full MO int. array. (only relevant
  !                     if each process stors more than one batch)
  !           XStarts:  MO starting index of each X batch
  !           Xdims:    Number of MOs included in each X batch
  !
  !> Author:  Pablo Baudin
  !> Date:    October 2013
  subroutine get_MO_batches_info(PRbatchInfo, dimP, dimR, TotSize)

    implicit none

    integer, intent(in) :: TotSize
    integer, intent(inout) :: dimP, dimR

    type(MObatchInfo), intent(inout) :: PRbatchInfo

    integer :: PR_batch, Rbatch, Pbatch, RStart, PStart, nRb, nPb

    ! short cuts:
    nRb = PRbatchInfo%nRbatch
    nPb = PRbatchInfo%nPbatch

    ! Allocate arrays:
    call mem_alloc(PRbatchInfo%PR_index,   int(i8*(nRb*nPb+1),kind=long))
    call mem_alloc(PRbatchInfo%PR_packInd, int(i8*(nRb*nPb+1),kind=long))
    call mem_alloc(PRbatchInfo%PStarts,    int(i8*nRb*nPb,kind=long))
    call mem_alloc(PRbatchInfo%RStarts,    int(i8*nRb*nPb,kind=long))
    call mem_alloc(PRbatchInfo%Pdims,      int(i8*nRb*nPb,kind=long))
    call mem_alloc(PRbatchInfo%Rdims,      int(i8*nRb*nPb,kind=long))

    ! Initialization
    PRbatchInfo%PR_index(1) = 1 
    PR_batch = 1
    RStart = 1

    ! Loop over MO batches:
    BatchR: do Rbatch = 1, nRb

      ! get dimension of last R batch:
      if (Rbatch == nRb) dimR = TotSize - RStart + 1

      PStart = 1
      BatchP: do Pbatch = 1, nPb

        ! get dimension of last P batch:
        if (Pbatch == nPb) dimP = TotSize - PStart + 1
 
        ! Store info about this PR batch:
        PRbatchInfo%PStarts(PR_batch) = PStart
        PRbatchInfo%RStarts(PR_batch) = RStart
        PRbatchInfo%PDims(PR_batch)   = dimP
        PRbatchInfo%RDims(PR_batch)   = dimR

        ! index for next batch:
        PRbatchInfo%PR_index(PR_batch+1)  = PRbatchInfo%PR_index(PR_batch) &
                                          & + dimP*dimR*TotSize*TotSize

        PR_batch = PR_batch + 1
        PStart = PStart + dimP

      end do BatchP
      ! restore dimension of P batch
      dimP = dimR
      RStart = RStart + dimR
    end do BatchR

  end subroutine get_MO_batches_info


  !> Purpose: Read packed MO int. and calculate A2.2 contribution
  !           to the CCSD doubles residual. 
  !           Based on Patrick's routine: get_a22_and_prepb22_terms_ex
  !
  !> Author:  Pablo Baudin
  !> Date:    October 2013
  subroutine get_ccsd_residual_small_frag(pack_gmo,t1,t2,omega2, &
             & nbas,nocc,nvir,MOinfo)

    implicit none

    !> MO pack integrals; amplitudes and residuals:
    real(realk), intent(inout) :: omega2(nvir,nvir,nocc,nocc)
    real(realk), intent(in) :: pack_gmo(:), t1(nvir,nocc)
    real(realk), intent(in) :: t2(nvir,nvir,nocc,nocc) 
    integer, intent(in) :: nbas, nocc, nvir

    !> Batches info:
    type(MObatchInfo), intent(in) :: MOinfo
    integer :: ipack, ibatch, P_sta, P_end, R_sta, R_end
    integer :: dimMO, dimP, dimR, nPb, nRb, PR_batch

    !> Arrays stuff:
    real(realk), pointer :: x_ap(:), tmp0(:), tmp1(:), tmp2(:), tmp3(:)
    integer(kind=long) :: tmp0_size, tmp1_size, tmp2_size, tmp3_size
    integer :: n_ij, n_cd, n_pr, ncopy, rd, c

    !> Symmetric and Anti-Symmetric combinations of doubles amplitude:
    real(realk), pointer :: tpl(:), tmi(:)

    !> logical to perform contraction p<=r:
    logical :: pack

    !> debug:
    real(realk), external :: ddot


    ! Initialize stuff
    nullify(tpl)
    nullify(tmi)
    nullify(x_ap)
    nullify(tmp0)
    nullify(tmp1)
    nullify(tmp2)
    nullify(tmp3)


    ! Allocate memory to working arrays:
    dimMO = MOinfo%Pdims(1)
    n_ij = nocc*(nocc+1)/2
    n_cd = nvir*(nvir+1)/2

    call mem_alloc(tpl,  int(i8*n_ij*n_cd, kind=long))
    call mem_alloc(tmi,  int(i8*n_ij*n_cd, kind=long))
    call mem_alloc(x_ap, int(i8*nvir*nbas, kind=long))

    tmp0_size = max(nvir*nvir*dimMO*dimMO, n_ij*dimMO*dimMO, nvir*nvir*nocc*nocc)
    tmp1_size = max(nvir*nvir*dimMO*dimMO, nvir*dimMO*n_ij, dimMO*dimMO*n_ij)
    tmp2_size = max(n_cd*dimMO*dimMO, n_ij*dimMO*dimMO)

    tmp0_size = int(i8*tmp0_size, kind=long)
    tmp1_size = int(i8*tmp1_size, kind=long)
    tmp2_size = int(i8*tmp2_size, kind=long)
    tmp3_size = int(i8*2*dimMO*dimMO*n_ij, kind=long)

    call mem_alloc(tmp0, tmp0_size)
    call mem_alloc(tmp1, tmp1_size)
    call mem_alloc(tmp2, tmp2_size)
    call mem_alloc(tmp3, tmp3_size)


    ! Get symmetric and antisymmetric combinations of doubles amplitude:
    ! tpl and tmi => [c>=d; i<=j]
    call get_sym_and_anti(t2,nvir,nvir,.true.,nocc,nocc,.true.,tpl,tmi,.true.)

 
    ! Calculate transformation matrix:
    !  x_ap = delta_ap - t_ai delta_ip
    x_ap = 0.0E0_realk
    call daxpy(nvir*nocc, -1.0E0_realk, t1, 1, x_ap, 1)
    do c = 1, nvir
      x_ap(nvir*nocc + c + (c-1)*nvir) = 1.0E0_realk
    end do


    ! Begin loop over MO batches
    omega2 = 0.0E0_realk
    nPb = MOinfo%nPbatch
    nRb = MOinfo%nRbatch
    BatchPR: do PR_batch = 1, nRb*nPb

      P_sta = MOinfo%PStarts(PR_batch)
      dimP  = MOinfo%Pdims(PR_batch)
      R_sta = MOinfo%RStarts(PR_batch)
      dimR  = MOinfo%Rdims(PR_batch)

      ! INTEGRAL ARE NOT PACKED YET BUT SHOULD BE UNPACKED HERE:
      !
      ! get PcRd from PqRs:
      ipack = dimP*nbas*dimR*nocc + dimP*nocc + MOInfo%PR_index(PR_batch)
      ibatch = 1
      ncopy = dimP*nvir
      do rd = 1, dimR*nvir
        call dcopy(ncopy, pack_gmo(ipack), 1, tmp0(ibatch), 1)
        ibatch = ibatch + dimP*nvir
        ipack = ipack + dimP*nbas
      end do
      ! tmp0 = g[pc; rd]


      ! SYMMETRIC COMBINATION:
      !
      ! get g+(pc|rd) = g(pc|rd) + g(pd|rc)
      !     g-(pc|rd) = g(pc|rd) - g(pd|rc)
      !
      if (P_sta==R_sta.and.dimP==dimR) then
        n_pr = dimP*(dimR+1)/2
        pack = .true.
      else if (P_sta<R_sta) then
        n_pr = dimP*dimR
        pack = .false.
      else if (P_sta>R_sta) then
        ! p<=r empty:
        print *, 'batch empty start next one'
        cycle
      else
        print *, "batch not square", PR_batch, dimP, dimR
        stop 'problem with batching in CCSD A2.2 calculation'
      end if
      ! 
      ! get tmp1 = g[cd;pr] <-- g(pc|rd)
      call array_reorder_4d(1.0E0_realk,tmp0,dimP,nvir,dimR,nvir, &
                            & [2,4,1,3],0.0E0_realk,tmp1)
 
      call get_sym_and_anti(tmp1,nvir,nvir,.true.,dimP,dimR,pack, &
             & tmp2,tmp0)
      ! tmp2 = g+[c>=d; p<=r]
      ! tmp0 = g-[c>=d; p<=r]

      ! sigma+[p<=r; i<=j] = sum_[c>=d] g+[c>=d; p<=r]' t+[c>=d; i<=j]
      call dgemm('t','n', n_pr, n_ij, n_cd, 0.5E0_realk, tmp2, n_cd, &
                & tpl, n_cd, 0.0E0_realk, tmp3, n_pr)

      ! ANTI-SYMMETRIC COMBINATION:
      !
      ! sigma-[p<=r; i<=j] = sum_[c>=d] g-[c>=d; p<=r]' t-[c>=d; i<=j]
      call dgemm('t','n', n_pr, n_ij, n_cd, 0.5E0_realk, tmp0, n_cd, &
                & tmi, n_cd, 0.0E0_realk, tmp3(n_pr*n_ij+1), n_pr)


      ! Add contribution to doubles amplitude residual:
      call get_A22_MO_ccsd(nocc, nvir, dimP, dimR, n_ij, n_pr, P_sta, R_sta, &
                          & pack, tmp0, tmp1, tmp2, tmp3, x_ap, omega2)
 
    end do BatchPR

    ! Calculate norm of A2.2:
    print *, 'ccsd residual A2.2 norm:', ddot(nocc*nocc*nvir*nvir,omega2,1,omega2,1)


    ! Free arrays:
    call mem_dealloc(tpl)
    call mem_dealloc(tmi)
    call mem_dealloc(tmp0)
    call mem_dealloc(tmp1)
    call mem_dealloc(tmp2)
    call mem_dealloc(tmp3)

  end subroutine get_ccsd_residual_small_frag
  


  !> Purpose: Calculate the A2.2 part of the doubles amplitudes residual
  !
  !           sigma[p<=r; i<=j] = 0.5*(sigma+[p<=r; i<=j] + sigma-[p<=r; i<=j])
  !           sigma[p>=r; i<=j] = 0.5*(sigma+[p<=r; i<=j] - sigma-[p<=r; i<=j])
  !
  !           omega2 += sum_pr x_ap x_br sigma[p r; i<=j]
  !           omega2 += sum_pr x_bp x_ar sigma[p r; j<i]
  !
  !> Author:  Pablo Baudin
  !> Date:    October 2013
  subroutine get_A22_MO_ccsd(nocc, nvir, dimP, dimR, n_ij, n_pr, P_sta, R_sta, &
             & pack, tmp0, tmp1, tmp2, tmp3, x_ap, omega2)

    implicit none

    !> dimensions and indices for arrays:
    integer, intent(in) :: nocc, nvir, dimP, dimR, P_sta, R_sta
    !> contracted dimensions for arrays:
    integer, intent(in) :: n_ij, n_pr
    !> working arrays:
    real(realk), intent(inout) :: tmp0(:), tmp1(:), tmp2(:)
    !> working arrays with sigma+- on input:
    real(realk), intent(in) :: tmp3(:)
    !> doubles residual array:
    real(realk), intent(inout) :: omega2(nvir,nvir,nocc,nocc)
    !> transformation matrix:
    real(realk), intent(in) :: x_ap(:)
    !> sigma indices packed p<=r ??
    logical, intent(in) :: pack

    !> orbital indices for loops:
    integer :: ij, r, i, j, b, a
    !> variable indices for arrays:
    integer :: pos, pos1, pos2
    integer :: mv((nvir*nvir)/2), st


    ! Calculate sigma[p r; i<=j]
    !
    if (pack) then
      ! Index for sigma
      pos = 1
      do ij=1,n_ij
        do r=1,dimR
          ! Index for full sigma:
          pos1 = 1 + (r-1)*dimP + (ij-1)*dimP*dimR
          pos2 = r + (ij-1)*dimP*dimR
          !
          ! Get symmetric combination of sigmas:
          !
          ! sigma[p<r; i<=j] = 0.5*(tmp3.1): sigma+[p<=r; i<=j] 
          !                  + 0.5*(tmp3.2): sigma-[p<=r; i<=j]
          call dcopy(r-1, tmp3(pos), 1, tmp2(pos1), 1)
          call daxpy(r-1, 1.0E0_realk, tmp3(pos+n_pr*n_ij), 1, tmp2(pos1), 1)
          ! 
          ! Get anti-symmetric combination of sigmas:
          !
          ! sigma[p>=r; i<=j] = 0.5*(tmp3.1): sigma+[p<=r; i<=j] 
          !                   - 0.5*(tmp3.2): sigma-[p<=r; i<=j]
          call dcopy(r, tmp3(pos), 1, tmp2(pos2), dimP)
          call daxpy(r, -1.0E0_realk, tmp3(pos+n_pr*n_ij), 1, tmp2(pos2), dimP)
          !
          pos = pos + r
        end do
      end do
  
    ! IF NOT PACK: TREAT SYMMETRIC PART FIRST:
    else 
      ! Get symmetric combination of sigmas:
      !
      ! sigma[p r; i<=j] = 0.5*(tmp3.1): sigma+[p r; i<=j] 
      !                  + 0.5*(tmp3.2): sigma-[p r; i<=j]
      call dcopy(dimP*dimR*n_ij, tmp3, 1, tmp2, 1)
      call daxpy(dimP*dimR*n_ij, 1.0E0_realk, tmp3(1+n_pr*n_ij), 1, tmp2,1)
    end if

    ! tmp2 = sigma[p r; i<=j]
    !
    ! Transpose sigma matrix from sigma[p; r i<=j ] to sigma[r i<=j; p ]
    call mat_transpose(dimP, dimR*n_ij, 1.0E0_realk, tmp2, 0.0E0_realk, tmp0)
    call lsmpi_poke() 
    !
    ! Transform r -> b
    call dgemm('n','n',nvir,n_ij*dimP,dimR,1.0E0_realk,x_ap(1+(R_sta-1)*nvir),nvir, &
         & tmp0, dimR, 0.0E0_realk, tmp1, nvir)
    call lsmpi_poke()
    !
    ! Transform p -> a; order is now: sigma[a b i j]
    call dgemm('n','t',nvir,nvir*n_ij,dimP,1.0E0_realk,x_ap(1+(P_sta-1)*nvir),nvir, &
         & tmp1, nvir*n_ij, 0.0E0_realk, tmp0, nvir)
    call lsmpi_poke()


    ! Sum up sigma PR batches contributions to A2.2 part of CCSD residual:
    do j=nocc,1,-1
      do i=j,1,-1
        pos1=1+((i+j*(j-1)/2)-1)*nvir*nvir
        pos2=1+(i-1)*nvir*nvir+(j-1)*nocc*nvir*nvir
        if(j/=1) tmp0(pos2:pos2+nvir*nvir-1) = tmp0(pos1:pos1+nvir*nvir-1)
      enddo
    enddo
    do j=nocc,1,-1
      do i=j,1,-1
        pos1=1+(i-1)*nvir*nvir+(j-1)*nocc*nvir*nvir
        pos2=1+(j-1)*nvir*nvir+(i-1)*nocc*nvir*nvir
        if(i/=j) tmp0(pos2:pos2+nvir*nvir-1) = tmp0(pos1:pos1+nvir*nvir-1)
      enddo
    enddo
    do j=nocc,1,-1
      do i=j,1,-1
        pos1=1+(i-1)*nvir*nvir+(j-1)*nocc*nvir*nvir
        call alg513(tmp0(pos1:nvir*nvir+pos1-1),nvir,nvir,nvir*nvir,mv,(nvir*nvir)/2,st)
      enddo
    enddo
    call daxpy(int(nocc*nocc*nvir*nvir),1.0E0_realk,tmp0,1,omega2,1)
    call lsmpi_poke()



    ! IF NOT PACK: TREAT ANTISYMMETRIC PART:
    !
    if (.not.pack) then ! If NOT pack then first treat symmetric part:
      ! Get symmetric combination of sigmas:
      !
      ! sigma[r p; i<=j] = 0.5*(tmp3.1): sigma+[p r; i<=j] 
      !                  - 0.5*(tmp3.2): sigma-[p r; i<=j]
      call dcopy(dimP*dimR*n_ij, tmp3, 1, tmp2, 1)
      call daxpy(dimP*dimR*n_ij, -1.0E0_realk, tmp3(1+n_pr*n_ij), 1, tmp2,1)
      call array_reorder_4d(1.0E0_realk,tmp2,dimP,dimR,nvir,nvir, &
                            & [2,1,3,4],0.0E0_realk,tmp1)

      ! tmp1 = sigma[r p; i<=j]
      !
      ! Transpose sigma matrix from sigma[r; p i<=j ] to sigma[p i<=j; r ]
      call mat_transpose(dimR, dimP*n_ij, 1.0E0_realk, tmp1, 0.0E0_realk, tmp0)
      call lsmpi_poke() 
      !
      ! Transform p -> b
      call dgemm('n','n',nvir,n_ij*dimR,dimP,1.0E0_realk,x_ap(1+(P_sta-1)*nvir),nvir, &
           & tmp0, dimP, 0.0E0_realk, tmp1, nvir)
      call lsmpi_poke()
      !
      ! Transform r -> a; order is now: sigma[a b i j]
      call dgemm('n','t',nvir,nvir*n_ij,dimR,1.0E0_realk,x_ap(1+(R_sta-1)*nvir),nvir, &
           & tmp1, nvir*n_ij, 0.0E0_realk, tmp0, nvir)
      call lsmpi_poke()

      ! Sum up sigma PR batches contributions to A2.2 part of CCSD residual:
      do j=nocc,1,-1
        do i=j,1,-1
          pos1=1+((i+j*(j-1)/2)-1)*nvir*nvir
          pos2=1+(i-1)*nvir*nvir+(j-1)*nocc*nvir*nvir
          if(j/=1) tmp0(pos2:pos2+nvir*nvir-1) = tmp0(pos1:pos1+nvir*nvir-1)
        enddo
      enddo
      do j=nocc,1,-1
        do i=j,1,-1
          pos1=1+(i-1)*nvir*nvir+(j-1)*nocc*nvir*nvir
          pos2=1+(j-1)*nvir*nvir+(i-1)*nocc*nvir*nvir
          if(i/=j) tmp0(pos2:pos2+nvir*nvir-1) = tmp0(pos1:pos1+nvir*nvir-1)
        enddo
      enddo
      do j=nocc,1,-1
        do i=j,1,-1
          pos1=1+(i-1)*nvir*nvir+(j-1)*nocc*nvir*nvir
          call alg513(tmp0(pos1:nvir*nvir+pos1-1),nvir,nvir,nvir*nvir,mv,(nvir*nvir)/2,st)
        enddo
      enddo
      call daxpy(int(nocc*nocc*nvir*nvir),1.0E0_realk,tmp0,1,omega2,1)
      call lsmpi_poke()

    end if

  end subroutine get_A22_MO_ccsd


  subroutine ultra_simple_a22(gao,xocc,xvirt,yocc,yvirt,t1,t2,nvir,nocc)

    implicit none
 
    type(array4), intent(inout) :: t2 
    integer, intent(in) :: nocc,nvir
    type(array2), intent(inout) :: xocc,xvirt,yocc,yvirt,t1
    type(array4), intent(inout) :: gao
    type(array4) :: abcd, tmp1

    !> debug:
    real(realk), pointer :: x_ap(:), newxvirt(:), cmo(:,:)
    real(realk) :: norm
    integer :: a,b,c,d
    real(realk), external :: ddot


    abcd = get_gmo_simple(gao,xvirt,yvirt,xvirt,yvirt)
    call array4_reorder(abcd,[2,4,1,3]) ! (ab|cd) > [b d a c] 
    tmp1 = array4_init([nvir,nvir,nocc,nocc]) ! tmp1[ab,ij]
    call array4_contract2(abcd,t2,tmp1)
    call array4_reorder(tmp1,[1,3,2,4]) ! -> tmp1[ai,bj]
    call array4_free(abcd)
 
    print *, 'debug: ccsd residual A2.2 norm:', tmp1*tmp1
    call array4_free(tmp1)

  end subroutine ultra_simple_a22

  !> Purpose: Generalization of Patrick's subroutine to get symmetric
  !           and antisymmetric combination of arrays.
  !           dim1 and dim2 have to be identical. And contraction over
  !           the indices is done depending on pack12/34.
  !           ampli is an optional parameter needed to take care of the
  !           diagonal elements of the symmetric amplitudes in the CCSD code
  !
  !> Author:  Pablo Baudin
  !> Date:    October 2013
  subroutine get_sym_and_anti(array,dim1,dim2,pack12,dim3,dim4,pack34, &
             & sym,anti,ampli)

    implicit none

    !> dimension of array
    integer, intent(in) :: dim1, dim2, dim3, dim4
    !> pack dimensions 12/34:
    logical, intent(in) :: pack12, pack34
    !> full array on input 
    real(realk), intent(in) :: array(dim1,dim2,dim3,dim4)
    !> symmetric combination on output with p>=q; r<=s depending on pack12/pack34
    real(realk), pointer, intent(inout) :: sym(:)
    !> anti-symmetric combination on output with p>=q; r<=s depending on pack12/pack34
    real(realk), pointer, intent(inout) :: anti(:)
    !> if array on input is amplitude we may take care of the diagonal 
    !  part of the symmetric combination for the CCSD code.
    logical, optional, intent(in) :: ampli
    logical :: dodiag

    integer :: q, r, s, dim12, ind12, ind34, ncopy

    ! do not treat diagonal by default:
    dodiag = .false.
    if (present(ampli)) dodiag = ampli

    ! sanity check:
    if (dim1/=dim2) then
      call lsquit("ERROR(get_sym_and_anti): dimensions of array not &
                 & adapted symmetric and antisymmetric combinations",DECinfo%output)
    end if
    if ((pack12.and.dim1/=dim2).or.(pack34.and.dim3/=dim4)) then
      call lsquit("ERROR(get_sym_and_anti): dimensions of array not &
                 & adapted for contraction",DECinfo%output)
    endif

    ! get the dimension of the 2 first indices combined:
    if (pack12) then
      dim12 = dim1*(dim2+1)/2
    else 
      dim12 = dim1*dim2
    end if

    do s=1,dim4
      do r=1,dim3
        ind12 = 1
        if (pack34) then 
          ind34 = r + (s-1)*s/2
        else 
          ind34 = r + (s-1)*dim3
        end if

        do q=1,dim2
          if (pack12) then
            ncopy = dim1 - q + 1
          else
            ncopy = dim1
          end if

          !copy the elements cd into sym and anti array:
          call dcopy(ncopy, array(q,q,r,s), 1, sym(ind12+(ind34-1)*dim12), 1)
          call dcopy(ncopy, array(q,q,r,s), 1, anti(ind12+(ind34-1)*dim12), 1)

          !add and subtract the counterparts for sym and antisym array respectively:
          call daxpy(ncopy,-1.0E0_realk,array(q,q,r,s),dim1,anti(ind12+(ind34-1)*dim12),1)
          call daxpy(ncopy,1.0E0_realk,array(q,q,r,s),dim1,sym(ind12+(ind34-1)*dim12),1)

          if (dodiag) then
            !take care of symmetric combination diagonal elements to be only half
            if (pack12) then 
              sym(ind12+(ind34-1)*dim12)=0.5E0_realk*sym(ind12+(ind34-1)*dim12)
            else
              sym(q+ind12+(ind34-1)*dim12)=0.5E0_realk*sym(q+ind12+(ind34-1)*dim12)
            end if             
          end if
         
          ! count in the triangular part of the matrix
          if (pack12) then 
            ind12 = ind12 + dim1 - q + 1
          else
            ind12 = ind12 + dim1
          end if
        enddo
      enddo
    enddo
  end subroutine get_sym_and_anti
  ! end pablo


end module cc_debug_routines_module
