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
   use screen_mod
 

   ! DEC DEPENDENCIES (within deccc directory)   
   ! *****************************************
   use crop_tools_module
   use array2_simple_operations
   use array4_simple_operations
   use mp2_module
   use ccintegrals
   use ccsd_module
   use ccsdpt_module
   use orbital_operations
   use rpa_module
   type SpaceInfo
     integer              :: n,ns1,ns2,pno,red1,red2
     integer, pointer     :: iaos(:)
     real(realk), pointer :: d(:,:)
     real(realk), pointer :: s1(:,:),s2(:,:)
     logical              :: allocd
   end type SpaceInfo
   integer,parameter :: SOLVE_AMPLITUDES  = 1
   integer,parameter :: SOLVE_MULTIPLIERS = 2
   

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
        & t1_final, t2_final, VOVO, .false.,SOLVE_AMPLITUDES)

     call array4_free(VOVO)

     call ccsolver_debug(ccmodel,Co_f, Cv_f, fock_f, nbasis, nocc, nvirt, &
        & mylsitem, ccPrintLevel, fragment_job, ppfock_f, qqfock_f, ccenergy, &
        & t1_final, t2_final, VOVO, .false., SOLVE_MULTIPLIERS, m2 = mult2, m1 = mult1)

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
     type(decfrag),optional,intent(in)   :: fraginfo
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
     logical :: restart, u_pnos
     type(array) :: govov


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
     endif



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

     ! prevent if explicitly requested or if PNOs are requested

     if(u_pnos.or.DECinfo%CCSDpreventcanonical)then
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
           call ccsolver_local_can_trans(nocc,nvirt,nbasis,Uocc,Uvirt,vovo=t2_final%val,vo=t1_final%val)
         else
           call ccsolver_local_can_trans(nocc,nvirt,nbasis,Uocc,Uvirt,vovo=t2_final%val)
         endif

       elseif(u_pnos)then
         call ccsolver_local_can_trans(nocc,nvirt,nbasis,Uocc,Uvirt,vovo=m2%val)
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
     if(decinfo%ccmodel /= MODEL_MP2 .and.decinfo%ccmodel /= MODEL_RPA ) then
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
       govov = array_minit([nocc,nvirt,nocc,nvirt],4,local=.true.,atype='TDAR')
       call array_zero(govov)

       call  wrapper_to_get_t1_free_gmo(nbasis,nocc,nvirt,Co%val,Cv2%val,&
             & govov,ccmodel,mylsitem)
     end if


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
           endif

           ifock = getInactiveFock_simple(h1,gao,xocc,yocc,nocc,nbasis)
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
              if(u_pnos)then
                call lsquit("ERROR(cc_driver_debug):only one of get_mult,&
                &and use_pnos should be true at the same time",-1)
              endif

              call get_ccsd_multipliers_simple(omega1(iter)%val,omega2(iter)%val,t1_final%val&
              &,t2_final%val,t1(iter)%val,t2(iter)%val,gao,xocc%val,yocc%val,xvirt%val,yvirt%val&
              &,nocc,nvirt,nbasis,MyLsItem)
           
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

              if(.not.fragment_job)then
                call get_ccsd_residual_pno_style(t1(iter)%val,t2(iter)%val,omega1(iter)%val,&
                &omega2(iter)%val,nocc,nvirt,nbasis,xocc%val,xvirt%val,yocc%val,yvirt%val,mylsitem,&
                &fragment_job,m2%val,ppfock%val,qqfock%val,delta_fock%val,iter)
              else
                call get_ccsd_residual_pno_style(t1(iter)%val,t2(iter)%val,omega1(iter)%val,&
                &omega2(iter)%val,nocc,nvirt,nbasis,xocc%val,xvirt%val,yocc%val,yvirt%val,mylsitem,&
                &fragment_job,m2%val,ppfock%val,qqfock%val,delta_fock%val,iter,f=fraginfo)
              endif

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
#ifdef VAR_MPI
             call RPA_residualpar(Omega2(iter),t2(iter),govov%elm1,ppfock,qqfock,nocc,nvirt)
#else
             call RPA_residualdeb(Omega2(iter),t2(iter),govov%elm1,ppfock,qqfock,nocc,nvirt)
#endif
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
           EnergyForCCmodel: if(CCmodel==MODEL_MP2) then  
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
     if (ccmodel==MODEL_RPA) call array_free(govov)

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
     if(.not.DECinfo%use_singles)then
       t1_final = array2_init([nvirt,nocc])
     endif

     ! Write finalization message
     !---------------------------
     call print_ccjob_summary(break_iterations,get_mult,fragment_job,last_iter,&
     &ccenergy,ttotend_wall,ttotstart_wall,ttotend_cpu,ttotstart_cpu,t1_final,t2_final)


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


     !transform back to original basis   
     if(DECinfo%use_singles)then
       call ccsolver_can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
       &vovo=t2_final%val,vo=t1_final%val)
       call ccsolver_can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
       &vovo=VOVO%val)
     else
       call ccsolver_can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
       &vovo=t2_final%val)
       call ccsolver_can_local_trans(nocc,nvirt,nbasis,Uocc,Uvirt,&
       &vovo=VOVO%val)
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
     character(ARR_MSG_LEN) :: msg

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


  !> \author Patrick Ettenhuber
  !> \date November 2013
  !> \brief this subroutine calculates the ccsd residual by transforming each
  !>        doubles amplitudes to their respective set of PNO's and then
  !>        transforming the result vector back to the reference basis. 
  subroutine get_ccsd_residual_pno_style(t1,t2,o1,o2,no,nv,nb,xo,xv,yo,yv,&
             &mylsitem,fj,t_mp2,oof,vvf,ifo,iter,f)
    implicit none
    !ARGUMENTS
    integer, intent(in) :: no, nv, nb
    real(realk), intent(inout) :: t1(nv,no), t2(nv,no,nv,no)
    real(realk), intent(inout) :: o1(nv,no), o2(nv,no,nv,no)
    real(realk), intent(in) :: xo(nb,no), xv(nb,nv), yo(nb,no), yv(nb,nv),ifo(nb,nb)
    type(lsitem), intent(inout) :: mylsitem
    real(realk), intent(in) :: t_mp2(nv,no,nv,no)
    real(realk), intent(inout) :: oof(no,no),vvf(nv,nv)
    logical, intent(in) :: fj
    integer, intent(in) :: iter
    type(decfrag),intent(in),optional :: f
    !INTERNAL VARIABLES
    type(SpaceInfo),pointer :: pno_cv(:),pno_S(:)
    type(array),pointer :: pno_o2(:),pno_t2(:),pno_gvvvv(:),pno_govov(:),pno_gvovo(:)
    integer :: nspaces,ns,ns2,ns3,c,nc,nc2
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
    logical :: skiptrafo,skiptrafo2,save_gvvvv_is,with_screening,cyc,use_triangular
    real(realk), pointer :: iFock(:,:), Dens(:,:)
    integer(kind=8) :: maxsize, myload
    integer :: pair,paircontribs,paircontrib(2,2)
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
    !real(realk) :: ref(no*nv*nv*no), ref1(no*nv), u(nv,no,nv,no)
    real(realk), parameter :: p20 = 2.0E0_realk
    real(realk), parameter :: p10 = 1.0E0_realk
    real(realk), parameter :: m10 = -1.0E0_realk
    real(realk), parameter :: m05 = -0.5E0_realk
    real(realk), parameter :: p05 = 0.5E0_realk
    real(realk), parameter :: nul = 0.0E0_realk

    
    o2v2 = (i8*no**2)*nv**2
    with_screening = .true.
    use_triangular = .true.
    master         = .true.

    paircontribs = 2
    paircontrib(1:2,1) = [1,2]
    paircontrib(1:2,2) = [2,1]

    if(fj.and..not.present(f))call lsquit("ERROR(get_ccsd_residual_pno_style):wrong input fj without f",-1)

    !DETERMINE NUMBER OF SPACES TO BE CONSIDERED
    if(fj)then
                !COUNT PAIRS OUTSIDE EOS                        !COUNT PAIRS WITH 1 IDX IN EOS   !EOS
      nspaces = ( no - f%noccEOS ) * ( no - f%noccEOS + 1) / 2 + f%noccEOS * ( no - f%noccEOS ) + 1
    else
                !ALL PAIRS
      nspaces = no * ( no + 1 ) / 2
    endif

    allocate(pno_cv(nspaces))
    allocate(pno_S(nspaces*(nspaces-1)/2))
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
    !===============================================================
    !begin setting up all density matrices and finding the PNO basis
    !===============================================================
    !DEBUG
    !t2 = t_mp2
    !write(msg,*)'DEBUG t2:'
    !call print_norm(t2,(i8*nv**2)*no**2,msg)

    ! if  we have a fragment job the basis of the atomic (pair) site has to be
    ! treated in a special way, either LO or FNO basis, this will be element 1
    ! in all the array arrays
    if(fj)then
      call get_pno_trafo_matrices(no,nv,nb,t_mp2,pno_cv,nspaces,fj,save_gvvvv_is,f=f)
      spacemax = max(f%noccEOS,2)
    else
      call get_pno_trafo_matrices(no,nv,nb,t_mp2,pno_cv,nspaces,fj,save_gvvvv_is)
      spacemax = 2
    endif


    !Get all the overlap matrices necessary
    call get_pno_overlap_matrices(no,nv,pno_cv,pno_S,nspaces,with_screening)

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
    IF(mylsitem%setting%scheme%cs_screen .OR. &
         & mylsitem%setting%scheme%ps_screen)THEN
       call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
            & nb,nbatchesAlpha,nbatchesGamma,&
            & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
       call II_getBatchOrbitalScreenK(DecScreen,mylsitem%setting,&
            & nb,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
            & batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
    ENDIF

    myload = 0

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
    !!write(*,*)' DEBUG A2/TOT:',sqrt(nnorm),sqrt(norm)
    !
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
    !!write (*,*)' DEBUG B2/TOT:',sqrt(nnorm),sqrt(norm)
    !

    !!DEBUG: C2 term
    !!**************
    !call array_reorder_4d( p10, goovv, no, no, nv, nv, [4,1,2,3], nul, w1 ) ! kjac -> ckja
    !call array_reorder_4d( p10, t2,    nv, no, nv, no, [2,3,4,1], nul, w2 ) ! aldj -> ldja
    !call array_reorder_4d( p10, govov, no, nv, no, nv, [4,1,3,2], nul, w3 ) ! kdlc -> ckld
    !!C intermediate w1(ckja) -  0.5 w3(ckld) w2(ldja)
    !call dgemm( 'n', 'n', no*nv, no*nv, no*nv, m05, w3, no*nv, w2, no*nv, p10, w1, no*nv)
    !call array_reorder_4d( p10, t2,    nv, no, nv, no, [1,4,3,2], nul, w2 ) ! bkci -> bick
    !call dgemm( 'n', 'n', no*nv, no*nv, no*nv, m10, w2, no*nv, w1, no*nv, nul, w3, no*nv)
    !!USE THE SYMMETRIZED CONTRIBUTION, i.e. P_{ij}^{ab} (1+0.5P_{ij}) * w3
    !call array_reorder_4d( p10, w3, nv, no, no, nv, [4,2,1,3], nul, w2 ) ! bija -> aibj
    !call array_reorder_4d( p05, w3, nv, no, no, nv, [4,3,1,2], p10, w2 ) ! bjia -> aibj
    !call array_reorder_4d( p10, w3, nv, no, no, nv, [1,3,4,2], p10, w2 ) ! ajib -> aibj
    !call array_reorder_4d( p05, w3, nv, no, no, nv, [1,2,4,3], p10, w2 ) ! aijb -> aibj
    !
    !ref = ref + w2(1:o2v2)

    !call print_norm(w2,o2v2,nnorm,.true.)
    !call print_norm(ref,o2v2,norm,.true.)
    !!write (*,*)' DEBUG C2/TOT:',sqrt(nnorm),sqrt(norm)
  
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
    !!write (*,*)' DEBUG D2/TOT:',sqrt(nnorm),sqrt(norm)
    !

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
    !!write (*,*)' DEBUG E21/TOT:',sqrt(nnorm),sqrt(norm)

    !!part 2
    !call ass_D2to1(oof,h1,[no,no])
    !w1(1:no**2) = h1(1:no**2)
    !h1 => null()
    !call array_reorder_4d( p10, govov, no, nv, no, nv, [1,4,3,2], nul, w2)
    !call dgemm('n','n',no,no,nv*no*nv, p10, w2,no,u,nv*no*nv, p10, w1, no)
    !call dgemm('n','n',nv*no*nv,no,no, m10, t2, nv*no*nv,w1,no, nul, w2, nv*no*nv)
    !w3(1:o2v2) = w2(1:o2v2)
    !call array_reorder_4d( p10, w3, nv, no, nv, no, [3,4,1,2], p10, w2)
    !
    !ref = ref + w2(1:o2v2)

    !call print_norm(w2,o2v2,nnorm,.true.)
    !call print_norm(ref,o2v2,norm,.true.)
    !write (*,*)' DEBUG E22/TOT:',sqrt(nnorm),sqrt(norm)

    !
    !ref1 = vof

    !!DEBUG SINGLES A1
    !!****************
    !call array_reorder_4d( p10, u, nv, no, nv, no, [3,2,1,4], nul, w1) ! ckdi -> dkci
    !call dgemm( 'n','n',nv, no, nv**2*no, p10,gvvov,nv,w1,nv**2*no,nul,w2,nv)
    !
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
    !$OMP ns2,ns3,nc,nc2) SHARED(pno_cv,pno_s,pno_t2,gvovo,goovv,gvvvv,&
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

      o   => pno_o2(ns)%elm1


      !!!!!!!!!!!!!!!!!!!!!!!!!
      !!!  A2 Term !!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!

      !A2.1
      !Get the integral contribution, sort it first like the integrals then transform it
      call extract_from_gvovo_transform_add(gvovo,w1,w2,d,o,idx,pno,no,pnv,nv,pno_cv(ns)%n)
     
      !A2.2
      call add_A22_contribution_simple(gvvvv,w1,w2,w3,d,t,o,pno,no,pnv,nv,pno_cv(ns)%n)

      
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
           &o,w1,w2,w3,w4,w5,goovv,govov,Lvoov,oof,p_idx,p_nidx,oidx1,oidx2,nspaces)



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
    call backtransform_omegas(pno_o2,pno_cv,o2,nspaces,no,nv)

    !call print_norm(o2,o2v2)
    !call print_norm(o1,i8*no*nv)

    do ns = 1, nspaces

      do ns2 = 1, ns-1

        c = (ns2 - ns + 1) + ns*(ns-1)/2
        if( pno_S(c)%allocd )then
          call mem_dealloc( pno_S(c)%iaos )
          call mem_dealloc( pno_S(c)%d    )
          if( with_screening )then
            call mem_dealloc( pno_S(c)%s1   )
            call mem_dealloc( pno_S(c)%s2   )
          endif
        endif

      enddo

      if( pno_cv(ns)%allocd )then
        call mem_dealloc( pno_cv(ns)%iaos )
        call mem_dealloc( pno_cv(ns)%d    )
        call array_free( pno_t2(ns) )
        call array_free( pno_o2(ns) )
      endif

    enddo

    deallocate( pno_cv )
    deallocate( pno_S )
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

  subroutine get_overlap_idx(n1,n2,cv,idx,nidx)
    implicit none
    integer,intent(in) :: n1,n2
    type(SpaceInfo), intent(in) :: cv(:)
    integer, intent(out) :: idx(:,:),nidx 
    integer :: nc1,nc2

    if(n1/=n2)then
      nidx = 0
      do nc1=1,cv(n1)%n
        do nc2=1,cv(n2)%n
          if(cv(n1)%iaos(nc1) == cv(n2)%iaos(nc2))then
            nidx = nidx + 1
            ! pos in first that equals second 
            idx(nidx,1) = nc1
            ! pos in second that equals first
            idx(nidx,2) = nc2
            ! the aos idx they refer to
            idx(nidx,3) = cv(n2)%iaos(nc2)
          endif
        enddo
      enddo
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
    type(SpaceInfo), intent(in) :: pS(:)
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

    elseif(n2>n1)then

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
    type(SpaceInfo),intent(in) :: cv(ns)
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
                  !print *,"found",n1,n2,cv(n1)%iaos(k),cv(n2)%iaos(l)
                  cntr = cntr+1
                  p_nidx(n1)  = cntr
                  p_idx(cntr,n1) = n2
                  cycle SpaceLoop2
                endif
              enddo
            endif
          enddo SpaceLoop2
        enddo idxloop
      endif
    enddo SpaceLoop

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
    type(SpaceInfo),intent(in) :: pno_cv(n)
    real(realk), intent(inout) :: o2(nv,no,nv,no)
    integer :: ns,pno,pnv,i,j,a,b
    real(realk), pointer :: tmp1(:),tmp2(:),d(:,:),po2(:,:,:,:), w1(:,:,:,:)
    real(realk) :: one,nul
    integer, pointer :: idx(:)

    one = 1.0E0_realk
    nul = 0.0E0_realk

    call mem_alloc(tmp1,nv**2*no**2)
    call mem_alloc(tmp2,nv**2*no**2)

    do ns = 1, n

      if(pno_cv(ns)%allocd)then

        pno =  pno_cv(ns)%n
        pnv =  pno_cv(ns)%ns2
        d   => pno_cv(ns)%d
        po2 => pno_o2(ns)%elm4
        idx => pno_cv(ns)%iaos

        call dgemm( 'n', 'n', nv, pno**2*pnv, pnv, one, d, nv, po2, pnv, nul, tmp1, nv )
        
        call array_reorder_4d( one, tmp1, nv, pno, pnv, pno, [3,4,1,2], nul, tmp2)

        call dgemm( 'n', 'n', nv, pno**2*nv, pnv, one, d, nv, tmp2, pnv, nul, tmp1, nv )


        !sort the contribution back and add up, again, because we assume a
        !symmetrized contribution in pno_o2 we can add up without taking care
        call ass_D1to4(tmp1,w1,[nv,pno,nv,pno])

        if(pno/=2)then
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
          do j = 1, pno
            do b = 1, nv
              do i = j + 1, pno
                do a = 1, nv
                  o2(a,idx(i),b,idx(j)) = o2(a,idx(i),b,idx(j)) + w1(a,i,b,j)
                  o2(a,idx(j),b,idx(i)) = o2(a,idx(j),b,idx(i)) + w1(a,j,b,i)
                enddo
              enddo
            enddo
          enddo
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
    type(SpaceInfo),intent(in) :: cv(:),S(:)
    logical, intent(out) :: cyc
    integer :: Sidx

    cyc = .false.    

    !if the trafo matrix has been screened away, cycle
    if(n1>n2)then

      !trafo from n1 to n2 
      Sidx      =  (n2 - n1 + 1) + n1 * (n1 - 1 )/2
      cyc       =  .not. S(Sidx)%allocd

    elseif(n2>n1)then

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
  !overlap specifies wheter the overlap is the first or the second argument in a
  !gemm. right now the routine assumes that A is never transposed. optionally
  !the pointers ptr and ptr2 can be assoicated with the transformed data and the
  !matrix not containing relevant data. this simplifies the code in the loops a
  !bit
  !> \author Patrick Ettenhuber
  !> \date december 2013
  ! 
  ! TODO: make this routine independent of the internal mem_allocations
  subroutine do_overlap_trafo(ns1,ns2,pos_of_overlap,S,m,n,k,A,C,ptr,ptr2)
    implicit none
    integer, intent(in) :: pos_of_overlap,m,n,k,ns1,ns2
    type(SpaceInfo),intent(inout) :: S(:)
    real(realk),intent(in),target :: A(:)
    real(realk),intent(inout),target :: C(:)
    real(realk),pointer,optional :: ptr(:),ptr2(:)
    real(realk),pointer :: S1(:,:)
    real(realk),pointer :: tmp1(:),tmp2(:), U(:,:), VT(:,:)
    integer :: ldS1,ldU,ldVT
    logical :: skiptrafo
    character :: tr1,tr2
    
    
    call get_overlap_ptr(ns1,ns2,S,tr1,tr2,skiptrafo,S1,ldS1,U=U,ldU=ldU,VT=VT,ldVT=ldVT)


    if(.not. skiptrafo)then

      select case(pos_of_overlap)
      case(1)

        call mem_alloc(tmp1,ldS1 * n)
        call mem_alloc(tmp2,ldVT * n)
        if(tr2=='t')then
          call dgemm('t','n',ldS1, n,ldU,  1.0E0_realk, U,ldU,A,k,0.0E0_realk,tmp1,ldS1)
          call dgemm('t','n',ldVT,n,ldS1,1.0E0_realk, S1,ldS1,tmp1,ldS1,0.0E0_realk,tmp2,ldVT)
          call dgemm('t','n',m,n,ldVT,  1.0E0_realk, VT,ldVT,tmp2,ldVT,0.0E0_realk,C,m)
        elseif(tr2=='n')then
          call dgemm('n','n',ldVT, n,k,  1.0E0_realk, VT,ldVT,A,k,0.0E0_realk,tmp2,ldVT)
          call dgemm('n','n',ldS1,n,ldVT,1.0E0_realk, S1,ldS1,tmp2,ldVT,0.0E0_realk,tmp1,ldS1)
          call dgemm('n','n',m,n,ldS1,  1.0E0_realk, U,ldU,tmp1,ldS1,0.0E0_realk,C,m)
        else
          call lsquit("ERROR(do_overlap_trafo):this should never happen, check get_overlap_ptr",-1)
        endif
        call mem_dealloc(tmp1)
        call mem_dealloc(tmp2)
        !call get_overlap_ptr(ns1,ns2,S,tr1,tr2,skiptrafo,S1,ldS1)
        !call dgemm(tr2,'n',m,n,k,1.0E0_realk,S1,ldS1,A,k,0.0E0_realk,C,m)

      case(2)
        call mem_alloc(tmp1,ldS1 * m)
        call mem_alloc(tmp2,ldVT * m)
        if(tr1=='t')then
          call dgemm('n','t',m,ldVT,k, 1.0E0_realk, A,m,VT,ldVT,0.0E0_realk,tmp2,m)
          call dgemm('n','t',m,ldS1,ldVT,1.0E0_realk, tmp2,m,S1,ldS1,0.0E0_realk,tmp1,m)
          call dgemm('n','t',m, n,ldS1,  1.0E0_realk, tmp1,m,U,ldU,0.0E0_realk,C,m)
        elseif(tr1=='n')then
          call dgemm('n','n',m, ldS1,k,  1.0E0_realk, A,m,U,ldU,0.0E0_realk,tmp1,m)
          call dgemm('n','n',m,ldVT,ldS1,1.0E0_realk, tmp1,m,S1,ldS1,0.0E0_realk,tmp2,m)
          call dgemm('n','n',m,n,ldVT,  1.0E0_realk, tmp2,m,VT,ldVT,0.0E0_realk,C,m)
        else
          call lsquit("ERROR(do_overlap_trafo):this should never happen, check get_overlap_ptr",-1)
        endif
        call mem_dealloc(tmp1)
        call mem_dealloc(tmp2)
        !call get_overlap_ptr(ns1,ns2,S,tr1,tr2,skiptrafo,S1,ldS1)
        !call dgemm('n',tr1,m,n,k,1.0E0_realk,A,m,S1,ldS1,0.0E0_realk,C,m)

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
    type(SpaceInfo),intent(in) :: pno_cv(n)
    type(SpaceInfo),intent(inout) :: pno_S(n*(n-1)/2)
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
    type(SpaceInfo),intent(in) :: cv(n)
    type(array),intent(inout) :: o2(n)
    integer :: nn, pnv,pno

    do nn=1,n

      pnv = cv(nn)%ns2
      pno = cv(nn)%n

      if(cv(nn)%allocd)then

        o2(nn) = array_init([pnv,pno,pnv,pno],4)
      
      endif
    enddo
  end subroutine init_pno_residual

  subroutine get_pno_amplitudes(t2,cv,pno_t2,n,no,nv)
    implicit none
    integer, intent(in) :: n,no,nv
    real(realk),intent(in) :: t2(:,:,:,:)
    type(SpaceInfo), intent(in) :: cv(n)
    type(array), intent(inout) :: pno_t2(n)
    real(realk), pointer :: tmp1(:),tmp2(:)
    real(realk), pointer :: w1(:,:,:,:)
    integer :: nn, pnv, pno, a, b, i, j
    call mem_alloc(tmp1,no**2*nv**2)
    call mem_alloc(tmp2,no**2*nv**2)

    do nn=1,n

      pnv = cv(nn)%ns2
      pno = cv(nn)%n

      if(cv(nn)%allocd)then

        pno_t2(nn) = array_init([pnv,pno,pnv,pno],4)

        call ass_D1to4(tmp1,w1,[nv,pno,nv,pno])
        do j=1,pno
        do b=1,nv
          do i=1,pno
          do a=1,nv
            w1(a,i,b,j) = t2(a,cv(nn)%iaos(i),b,cv(nn)%iaos(j))
          enddo
          enddo
        enddo
        enddo

        call dgemm('t','n',pnv,pno**2*nv,nv,1.0E0_realk,cv(nn)%d,nv,w1,nv,0.0E0_realk,tmp2,pnv)
        call array_reorder_4d(1.0E0_realk,tmp2,pnv,pno,nv,pno,[3,4,1,2],0.0E0_realk,tmp1)

        !the amplitudes are symmetric, also after the transformation,
        !therefore it is not important in which order they are stored
        call dgemm('t','n',pnv,pno**2*pnv,nv,1.0E0_realk,cv(nn)%d,nv,tmp1,nv,0.0E0_realk,pno_t2(nn)%elm1,pnv)
        w1 => null()

        !To avoid double counting -> this can be removed and savings introduced in
        !the algorithm if restrictions are put on the pair indices, but a bit more
        !complicated than this
        if(pno==2)then
          do i = 1, pno
            pno_t2(nn)%elm4(:,i,:,i) = 0.0E0_realk
          enddo
        endif
      endif

    enddo

    call mem_dealloc(tmp1)
    call mem_dealloc(tmp2)
  end subroutine get_pno_amplitudes

  subroutine get_pno_trafo_matrices(no,nv,nb,t_mp2,cv,n,fj,sgvvvvis,f)
    implicit none
    !ARGUMENTS
    integer, intent(in) :: no, nv, nb, n
    real(realk), intent(in) :: t_mp2(nv,no,nv,no)
    type(SpaceInfo),pointer :: cv(:)
    logical,intent(in) :: fj
    logical,intent(out) :: sgvvvvis
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
      cv(1)%n    = f%noccEOS
      cv(1)%iaos = f%idxo
      counter = 1

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
          doit= ( find_pos(i,j)/= -1 )

          !calculate the pair density matrix, diagonalize it, truncate the
          !respective transformation matrix and save it in c
          if(doit)then
            counter = find_pos(i,j)
            call calculate_pair_density_matrix(PD,t_mp2(:,i,:,j),nv,(i==j))
            call solve_eigenvalue_problem_unitoverlap(nv,PD,virteival,U)
            call truncate_trafo_mat_from_EV(U,virteival,nv,cv(counter))
            if(cv(counter)%allocd)then
              if(i==j)then
                cv(counter)%n = 1
                call mem_alloc(cv(counter)%iaos,cv(counter)%n)
                cv(counter)%iaos = [i]
                det_parameters = det_parameters + cv(counter)%ns2*cv(counter)%ns2*cv(counter)%n**2
              else
                cv(counter)%n = 2
                call mem_alloc(cv(counter)%iaos,cv(counter)%n)
                cv(counter)%iaos = [i,j]
                det_parameters = det_parameters + cv(counter)%ns2*cv(counter)%ns2*2
              endif
              calc_parameters = calc_parameters + cv(counter)%ns2*cv(counter)%ns2*cv(counter)%n**2
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
          if(cv(counter)%allocd)then
            if(i==j)then
              cv(counter)%n = 1
              call mem_alloc(cv(counter)%iaos,cv(counter)%n)
              cv(counter)%iaos = [i]
              det_parameters = det_parameters + cv(counter)%ns2*cv(counter)%ns2*cv(counter)%n**2
            else
              cv(counter)%n = 2
              call mem_alloc(cv(counter)%iaos,cv(counter)%n)
              cv(counter)%iaos = [i,j]
              det_parameters = det_parameters + cv(counter)%ns2*cv(counter)%ns2*2
            endif
            calc_parameters = calc_parameters + cv(counter)%ns2*cv(counter)%ns2*cv(counter)%n**2
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
    print *,"full gvvvv",nv**4," vs ",pno_gvvvv_size

    sgvvvvis = (pno_gvvvv_size<nv**4)

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
    type(SpaceInfo),intent(inout) :: NU
    real(realk), intent(in),optional :: ext_thr
    !INTERNAL
    integer :: i,nn
    real(realk) :: thr

    if(DECinfo%noPNOtrunc)then
      thr = -1.0*huge(thr)
    else
      thr = DECinfo%simplePNOthr
    endif
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

    if(DECinfo%PL>2.and.present(ext_thr))write(DECinfo%output,'("The FO trafo  matrix has dims",2I4)')n,nn

   
    if(DECinfo%noPNOtrafo.and..not.present(ext_thr))then
      NU%ns1 = n
      NU%ns2 = n
    else
      NU%ns1 = n
      NU%ns2 = nn
    endif

    if(nn<=0)then

      NU%allocd = .false.

    else

      call mem_alloc(NU%d,n,nn)

      if(DECinfo%noPNOtrafo)then

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

  subroutine extract_from_gvovo_transform_add(gvovo,w1,w2,d,o,idx,pno,no,pnv,nv,n)
    implicit none
    integer,intent(in) :: pno,no,pnv,nv,n
    integer,intent(in) :: idx(n)
    real(realk),intent(in) ::gvovo(nv,no,nv,no), d(nv,pnv)
    real(realk),intent(inout) :: o(nv*nv*no*no)
    real(realk),intent(out) :: w1(nv,pno,nv,pno)
    real(realk),pointer,intent(inout) :: w2(:)
    integer :: i,j,a,b
    real(realk), parameter :: nul = 0.0E0_realk
    real(realk), parameter :: p10 = 1.0E0_realk

    if(n==2.and.DECinfo%PNOtriangular)then
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

    if(n==2.and.DECinfo%PNOtriangular)then
      !transform integral contribution, use symmetry  gvovo(aibj) => gvovo(\bar{a} i \bar{b} j)
      call dgemm( 't', 'n', pnv, nv,  nv, p10, d,  nv, w1, nv, nul, w2, pnv )
      call dgemm( 'n', 'n', pnv, pnv, nv, p10, w2, pnv, d, nv, nul, o, pnv )
    else
      !transform integral contribution, use symmetry  gvovo(aibj) => gvovo(\bar{b} j \bar{a} i)
      call dgemm( 't', 'n', pnv, pno**2*nv, nv, p10, d, nv, w1, nv, nul, w2, pnv )
      call array_reorder_4d( p10, w2, pnv, pno, nv, pno, [3,4,1,2], nul, w1 )
      call dgemm( 't', 'n', pnv, pno**2*pnv, nv, p10, d, nv, w1, nv, nul, o, pnv )
    endif

  end subroutine extract_from_gvovo_transform_add


  subroutine add_A22_contribution_simple(gvvvv,w1,w2,w3,d,t,o,pno,no,pnv,nv,n)
    implicit none
    integer,intent(in) :: pno,no,pnv,nv,n
    real(realk), intent(in)    :: gvvvv(nv**4), t(pnv**2*pno**2), d(nv,pnv)
    real(realk), intent(out)   :: w1(:),w2(:),w3(:)
    real(realk), intent(inout) :: o(pnv**2*pno**2)
    real(realk), parameter :: nul = 0.0E0_realk
    real(realk), parameter :: p10 = 1.0E0_realk
    !A2.2
    !transform to basis of space gvvvv(acbd) => gvvvv(\bar{a}\bar{c}\bar{b}\bar{d})
    call dgemm( 't', 'n', nv**3,     pnv, nv, p10, gvvvv, nv, d, nv, nul, w1, nv**3     )
    call dgemm( 't', 'n', nv**2*pnv, pnv, nv, p10, w1   , nv, d, nv, nul, w2, nv**2*pnv )
    call dgemm( 't', 'n', nv*pnv**2, pnv, nv, p10, w2   , nv, d, nv, nul, w1, nv*pnv**2 )
    call dgemm( 't', 'n', pnv**3,    pnv, nv, p10, w1   , nv, d, nv, nul, w2, pnv**3    )
    !end transformation, integrals are now in the order gvvvv(\bar{a}\bar{c}\bar{b}\bar{d})
    !reorder corresponding integrals (w1) to gvvvv(\bar{a}\bar{b}\bar{c}\bar{d})and 
    call array_reorder_4d( p10, w2, pnv, pnv, pnv, pnv, [1,3,2,4], nul, w1 )

    if(n==2.and.DECinfo%PNOtriangular)then
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
    type(SpaceInfo), intent(inout) :: pno_cv(nspaces),pno_S(nspaces*(nspaces-1)/2)
    type(array), intent(in) :: pno_t2(nspaces)
    real(realk),pointer,intent(inout) :: o2_space(:)
    real(realk),pointer,intent(inout) :: w1(:),w2(:),w3(:),w4(:),w5(:)
    real(realk),intent(in) :: goooo(:),govov(:),vvf(:,:)
    character :: tr11,tr12,tr21,tr22
    real(realk),pointer :: p1(:,:,:,:), p2(:,:,:,:), p3(:,:,:,:), p4(:,:,:,:),h1(:), h2(:), r1(:,:),r2(:,:),d(:,:),d1(:,:)
    real(realk),pointer :: o(:),t(:),S1(:,:), t21(:)
    logical :: skiptrafo, cyc
    integer :: space, a,i,b,j, nv,pno,pnv, ns2
    integer :: pno1,pnv1,Sidx1,ldS1
    integer :: pair,paircontribs,paircontrib(2,2)
    integer :: order1(4)
    integer :: suborder(2)
    integer,pointer :: idx(:),idx1(:)
    integer(kind=8) :: o2v2
    real(realk), parameter :: nul =  0.0E0_realk
    real(realk), parameter :: p10 =  1.0E0_realk
    real(realk), parameter :: p20 =  2.0E0_realk
    real(realk), parameter :: m10 = -1.0E0_realk

    nv  =  pno_cv(ns)%ns1
    d   => pno_cv(ns)%d
    idx => pno_cv(ns)%iaos
    pnv =  pno_cv(ns)%ns2
    pno =  pno_cv(ns)%n
    t   => pno_t2(ns)%elm1
    o   => o2_space

    o2v2 = (i8*no**2)*nv**2

    paircontribs = 2
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
    
    
    FullSpaceLoop1: do ns2 = 1, nspaces

      call check_if_contributes(ns,ns2,pno_cv,pno_S,cyc)
 
      if(cyc)then

        cycle FullSpaceLoop1

      endif

      d1   => pno_cv(ns2)%d
      t21  => pno_t2(ns2)%elm1
      idx1 => pno_cv(ns2)%iaos
      pnv1 =  pno_cv(ns2)%ns2
      pno1 =  pno_cv(ns2)%n


      !!!!!!!!!!!!!!!!!!!!!!!!!
      !!!  E2 Term part1!!!!!!! - quadratic contribution
      !!!!!!!!!!!!!!!!!!!!!!!!!

      !Get the integral contribution, sort it first like the integrals then transform it, govov
      call ass_D1to4( w1,    p1, [pno1,nv,pno1,nv] )
      call ass_D1to4( govov, p2, [no,   nv,no, nv] )
      do j=1,pno1
      do b=1,nv
        do i=1,pno1
        do a=1,nv
          p1(i,a,j,b) = p2(idx1(i),a,idx1(j),b)
        enddo
        enddo
      enddo
      enddo
      p1 => null()
      p2 => null()
      
      !transform integral contribution, use symmetry  govov(ldkc) => govov(\bar{d} k l  \bar{c}) to the space of (ij) -> w2
      call dgemm( 'n', 'n', nv*pno1**2, pnv, nv, p10, w1, nv*pno1**2, d, nv, nul, w2, nv*pno1**2 )
      call array_reorder_4d( p10, w2, pno1, nv, pno1, pnv, [2,3,1,4], nul, w3 )
      call dgemm( 't', 'n', pnv1, pno1**2*pnv, nv, p10, d1, nv, w3, nv, nul, w1, pnv1 )

      ! Quadratic part of the E2 term use u^{bd}_{kl} (bkdl) as b,dkl
      call array_reorder_4d( p20, t21, pnv1, pno1, pnv1, pno1, [1,3,2,4], nul, w3)
      call array_reorder_4d( m10, t21, pnv1, pno1, pnv1, pno1, [1,3,4,2], p10, w3)
      call do_overlap_trafo(ns,ns2,1,pno_S,pnv,pno1*pnv1*pno1,pnv1,w3,w2,ptr=h1)

      !contract amplitudes in h1 with integrals in w1 and add to w4 : -1 * h1(bdkl) w1(dlkc) += w4(bc)
      call dgemm('n','n',pnv,pnv,pnv1*pno1*pno1,m10, h1,pnv,w1,pnv1*pno1*pno1,p10,w5,pnv)


      !!!!!!!!!!!!!!!!!!!!!!!!!
      !!!  B2 Term !!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!
      
      !prepare 4 occupied integral goooo for B2 term
      call ass_D1to4( w3(1:pno**2*pno1**2), p1, [pno,pno,pno1,pno1] )
      call ass_D1to4( goooo, p2, [ no, no,  no,  no] )
      !Get the integral contribution, sort it first like the integrals then transform it, govov
      call ass_D1to4( w1,    p3, [pno1,nv,pno1,nv] )
      call ass_D1to4( govov, p4, [no,  nv,no,  nv] )


      !CODE FOR PAIR SPACES WITH RESTRICTIONS
      if(pno_cv(ns2)%n==2.and.DECinfo%PNOtriangular)then

        !loop over pair contributions kl and lk
        do pair = 1, paircontribs

          p3(1,:,1,:) = p4(idx1(paircontrib(1,pair)),:,idx1(paircontrib(2,pair)),:)
          !transform integral contribution
          !govov(kcld) => govov(\bar{c} \bar{d} k l) or lk to the space of (ij) -> w2
          call dgemm( 'n', 'n', nv, pnv, nv, p10, w1, nv, d, nv, nul, w3, nv )
          call dgemm( 't', 'n', pnv, pnv, nv, p10, d, nv, w3, nv, nul, w4, pnv )

          p1(1,1,1,1) = p2(idx1(paircontrib(1,pair)),idx(1),idx1(paircontrib(2,pair)),idx(2))

          !sort the amplitudes and contract cidj -> ijcd, ijcd cdkl + ijkl = ijkl
          call array_reorder_4d( p10, t, pnv, pno, pnv, pno, [2,4,1,3], nul, w1 )
          call dgemm( 'n', 'n', pno**2, pno1**2, pnv**2, p10, w1, pno**2, w4, pnv**2, p10, w3, pno**2 )

          !contract the B intermediate in w3 with the amplitudes (kl) from the
          !inner loop and use the overlap to transform to the omega space, ijkl klab
          call array_reorder_2d( p10, t21, pnv1, pnv1 , paircontrib(:,pair), nul, w1 )

          call dgemm( 'n', 'n', pno**2, pnv1**2, pno1**2, p10, w3, pno**2, w1, pno1**2, dble(pair-1), w2, pno**2 )
          
        enddo

      !CODE FOR RECTANGULAR SPACES
      else
        do b=1,nv
        do j=1,pno1
          do a=1,nv
          do i=1,pno1
            p3(i,a,j,b) = p4(idx1(i),a,idx1(j),b)
          enddo
          enddo
        enddo
        enddo

        !transform integral contribution
        !govov(kcld) => govov(\bar{c} \bar{d} k l) to the space of (ij) -> w2
        call dgemm( 'n', 'n', nv*pno1**2, pnv, nv, p10, w1, nv*pno1**2, d, nv, nul, w2, nv*pno1**2 )
        call array_reorder_4d( p10, w2, pno1, nv, pno1, pnv, [2,4,1,3], nul, w1 )
        call dgemm( 't', 'n', pnv, pno1**2*pnv, nv, p10, d, nv, w1, nv, nul, w4, pnv )

        !prepare 4 occupied integral goooo for B2 term
        !if(pno_cv(ns)%n/=2)then
        do j=1,pno1
        do i=1,pno1
          do b=1,pno
          do a=1,pno
            p1(a,b,i,j) = p2(idx1(i),idx(a),idx1(j),idx(b))
          enddo
          enddo
        enddo
        enddo
        !sort the amplitudes and contract cidj -> ijcd, ijcd cdkl + ijkl = ijkl
        call array_reorder_4d( p10, t, pnv, pno, pnv, pno, [2,4,1,3], nul, w1 )
        call dgemm( 'n', 'n', pno**2, pno1**2, pnv**2, p10, w1, pno**2, w4, pnv**2, p10, w3, pno**2 )

        !contract the B intermediate in w3 with the amplitudes (kl) from the
       !inner loop and use the overlap to transform to the omega space, ijkl klab
        call array_reorder_4d( p10, t21, pnv1, pno1, pnv1, pno1, [2,4,1,3], nul, w1 )

        call dgemm( 'n', 'n', pno**2, pnv1**2, pno1**2, p10, w3, pno**2, w1, pno1**2, nul, w2, pno**2 )

      endif

      ! transform back, or in the case of ns==ns2 just order correctly
      if(ns==ns2)then
        call array_reorder_4d( p10, w2, pno, pno, pnv, pnv, [3,1,4,2], nul, w1 )
      else
        call do_overlap_trafo(ns,ns2,2,pno_S,pno**2*pnv1, pnv, pnv1,w2,w1)
        call array_reorder_4d( p10, w1, pno, pno, pnv1, pnv, [3,1,4,2], nul, w3 )
        call do_overlap_trafo(ns,ns2,1,pno_S,pnv,pno**2*pnv, pnv1,w3,w1)
      endif
   
      ! add up the correcly ordered contributions
      o = o + w1(1:pno**2*pnv**2)
  
    enddo FullSpaceLoop1


    !Add the E21 contribution
    call array_reorder_4d( p10, t, pnv, pno, pnv, pno, [3,4,1,2], nul, w1)
    call dgemm('n','n',pnv,pno*pnv*pno,pnv,p10,w5,pnv,w1,pnv,nul,w2,pnv)
    o = o + w2(1:pnv*pno*pnv*pno)
    call array_reorder_4d( p10, w2, pnv, pno, pnv, pno, [3,4,1,2], p10, o )

  end subroutine get_free_summation_for_current_aibj


  subroutine get_common_idx_summation_for_current_aibj(no,ns,pno_cv,pno_S,pno_t2,o2_space,&
             &w1,w2,w3,w4,w5,goovv,govov,Lvoov,oof,p_idx,p_nidx,oidx1,oidx2,nspaces)
    implicit none
    integer, intent(in) :: no,ns,nspaces
    type(SpaceInfo), intent(inout) :: pno_cv(nspaces),pno_S(nspaces*(nspaces-1)/2)
    type(array), intent(in) :: pno_t2(nspaces)
    real(realk),pointer,intent(inout) :: o2_space(:)
    real(realk),pointer,intent(inout) :: w1(:),w2(:),w3(:),w4(:),w5(:)
    real(realk),intent(in) :: goovv(:),govov(:),Lvoov(:),oof(:,:)
    integer,intent(in)    :: p_idx(:,:),p_nidx(:)
    integer,intent(inout) :: oidx1(:,:),oidx2(:,:)
    character :: tr11,tr12,tr21,tr22,TRamp_pos1(2),TRamp_pos2(2),trh1,trh2
    real(realk),pointer :: p1(:,:,:,:), p2(:,:,:,:), p3(:,:,:,:), p4(:,:,:,:),h1(:), h2(:), h3(:)
    real(realk),pointer :: r1(:,:),r2(:,:),d(:,:),d1(:,:),d2(:,:)
    real(realk),pointer :: o(:),t(:),S1(:,:), t21(:), t22(:), f1(:,:,:)
    logical :: skiptrafo,cyc
    integer :: space, a,i,b,j, nv,pno,pnv,nc,nc2,nidx1,nidx2,ns2,ns3
    integer :: pno1,pnv1,Sidx1,ldS1
    integer :: pno2,pnv2,Sidx2,ldS2
    integer :: pair1,pair2,paircontrib(2,2)
    integer,parameter :: paircontribs = 2
    integer :: order1(4), ldh1,ldh2,ldh3
    integer :: suborder(2),i_idx,j_idx,pos_in_res,ipos_in_res
    integer,pointer :: idx(:),idx1(:),idx2(:)
    integer(kind=8) :: o2v2
    real(realk), parameter :: nul =  0.0E0_realk
    real(realk), parameter :: p20 =  2.0E0_realk
    real(realk), parameter :: p10 =  1.0E0_realk
    real(realk), parameter :: p05 =  0.5E0_realk
    real(realk), parameter :: m05 = -0.5E0_realk
    real(realk), parameter :: m10 = -1.0E0_realk

    nv  =  pno_cv(ns)%ns1
    d   => pno_cv(ns)%d
    idx => pno_cv(ns)%iaos
    pnv =  pno_cv(ns)%ns2
    pno =  pno_cv(ns)%n
    t   => pno_t2(ns)%elm1
    o   => o2_space

    OneIdxSpaceLoop1: do nc = 1, p_nidx(ns)
      ! extract indices:
      ns2 = p_idx(nc,ns)

      call check_if_contributes(ns,ns2,pno_cv,pno_S,cyc)

      if( cyc )then

        cycle OneIdxSpaceLoop1

      endif

      call get_overlap_idx(ns,ns2,pno_cv,oidx1,nidx1)

      d1   => pno_cv(ns2)%d
      t21  => pno_t2(ns2)%elm1
      idx1 => pno_cv(ns2)%iaos
      pnv1 =  pno_cv(ns2)%ns2
      pno1 =  pno_cv(ns2)%n


      !!!!!!!!!!!!!!!!!!!!!!!!!
      !!!  C2 Term !!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!

      !Transform the integral contribution to the space of the current amps
      !get the order kjac -> ckja, please note, that this loop is only
      !inside OneIdxSpaceLoop2 because I only use 3 working matrices, p10
      !might easily move the following part outside the loop and add stuff
      !up during the loops
      call ass_D1to4( w1,    p1, [nv,pno1,pno,nv] )
      call ass_D1to4( goovv, p2, [no, no,nv,  nv] )
      do a=1,nv
      do j=1,pno
        do i=1,pno1
        do b=1,nv
          p1(b,i,j,a) = p2(idx1(i),idx(j),a,b)
        enddo
        enddo
      enddo
      enddo

      ! transform c to \bar(c} of (ki)  and a to \bar{a} of (ij)
      call dgemm('t','n', pnv1, pno1*pno*nv, nv,  p10, d1, nv, w1, nv, nul, w2, pnv1)
      call dgemm('n','n', pnv1*pno1*pno, pnv, nv, p10, w2, pnv1*pno1*pno, d, nv, nul, w4, pnv1*pno1*pno)

      !THE INNER CONTRACTION LOOP - BUILDING THE C INTERMEDIATE
      OneIdxSpaceLoop2: do nc2=1, p_nidx(ns)
        ! extract indices:
        ns3 = p_idx(nc2,ns)

        call check_if_contributes(ns,ns3,pno_cv,pno_S,cyc)

        if( cyc )then

          cycle OneIdxSpaceLoop2

        endif

        call get_overlap_idx(ns,ns3,pno_cv,oidx2,nidx2)

        d2   => pno_cv(ns3)%d
        t22  => pno_t2(ns3)%elm1
        idx2 => pno_cv(ns3)%iaos
        pnv2 =  pno_cv(ns3)%ns2
        pno2 =  pno_cv(ns3)%n

        !Get the integrals kdlc -> ckld and transform c and d to their
        !corresponding spaces, (iajb -> bija) 
        call ass_D1to4( w1,    p1, [nv,pno1,pno2,nv] )
        call ass_D1to4( govov, p2, [no, nv,no,   nv] )
        do a=1,nv
        do j=1,pno2
          do i=1,pno1
          do b=1,nv
            p1(b,i,j,a) = p2(idx1(i),a,idx2(j),b)
          enddo
          enddo
        enddo
        enddo

        ! transform c to \bar(c} in (ki) and d to \bar{d} in (lj)
        call dgemm('t','n',pnv1,pno1*pno2*nv,nv,  p10, d1, nv, w1, nv, nul, w3, pnv1)
        call dgemm('n','n',pnv1*pno1*pno2,pnv2,nv,p10, w3, pnv1*pno1*pno2, d2, nv, nul, w1, pnv1*pno1*pno2)

        !get the amplitudes in the correct order eldj -> ldje transform to a and contract to
        ! -0.5 w1(ckld) w2(ldja) += w4(ckja)
        call ass_D1to4( w3,  p3, [pno2,pnv2,nidx2,pnv2] )
        call ass_D1to4( t22, p2, [pnv2,pno2,pnv2,pno2] )
        do b=1,pnv2
        do j=1,nidx2
          do a=1,pnv2
          do i=1,pno2
            p3(i,b,j,a) = p2(a,i,b,oidx2(j,2))
          enddo
          enddo
        enddo
        enddo

        call do_overlap_trafo(ns,ns3,2,pno_S,pno2*pnv2*nidx2,pnv,pnv2,w3,w2,ptr=h1,ptr2=h2)

        call dgemm('n','n', pnv1*pno1, nidx2*pnv, pno2*pnv2, m05, w1, pnv1*pno1, h1, pno2*pnv2, nul, h2, pnv1*pno1)

        call ass_D1to4( h2, p3, [pnv1,pno1,nidx2,pnv] )
        call ass_D1to4( w4, p4, [pnv1,pno1,pno,pnv] )
        do a=1,pnv
        do j=1,nidx2
          do i=1,pno1
          do b=1,pnv1
            p4(b,i,oidx2(j,1),a) = p4(b,i,oidx2(j,1),a) + p3(b,i,j,a)
          enddo
          enddo
        enddo
        enddo

      enddo OneIdxSpaceLoop2

      !get the amplitudes, extract the necessary indices, 
      !reorder dkci -> dick :D transform to current space (bick) and do the contraction,
      !bick ckja = bija, do the permutation and addition of the contribution
      call ass_D1to4( w1,  p1, [pnv1,nidx1,pnv1,pno1] )
      call ass_D1to4( t21, p2, [pnv1,pno1,pnv1,pno1] )
      do j=1,pno1
      do b=1,pnv1
        do i=1,nidx1
        do a=1,pnv1
          p1(a,i,b,j) = p2(a,j,b,oidx1(i,2))
        enddo
        enddo
      enddo
      enddo

      call do_overlap_trafo(ns,ns2,1,pno_S, pnv,nidx1*pnv1*pno1, pnv1,w1,w2,ptr=h1,ptr2=h2)

      call dgemm('n','n', pnv*nidx1,pno*pnv, pno1*pnv1, m10, h1,pnv*nidx1, w4, pnv1*pno1, nul, h2, pnv*nidx1)
      call ass_D1to4( h2, p2, [pnv,nidx1,pno,pnv] )
      call ass_D1to4( o,  p1, [pnv,pno, pnv, pno ] )
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
        ns3 = p_idx(nc2,ns)

        call check_if_contributes(ns,ns3,pno_cv,pno_S,cyc)

        if( cyc )then

          cycle OneIdxSpaceLoop3

        endif

        call get_overlap_idx(ns,ns3,pno_cv,oidx2,nidx2)

        d2   => pno_cv(ns3)%d
        t22  => pno_t2(ns3)%elm1
        idx2 => pno_cv(ns3)%iaos
        pnv2 =  pno_cv(ns3)%ns2
        pno2 =  pno_cv(ns3)%n

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

        call do_overlap_trafo(ns,ns3,1,pno_S, pnv, nidx2*pno2*pnv2, pnv2,w3,w2,ptr=h1,ptr2=h2)
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

      call do_overlap_trafo(ns,ns2,1,pno_S,pnv,nidx1*pnv1*pno1,pnv1,w1,w2,ptr=h1,ptr2=h2)
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
        ns3 = p_idx(nc2,ns)

        call check_if_contributes(ns,ns3,pno_cv,pno_S,cyc)

        if( cyc )then

          cycle OneIdxSpaceLoop4

        endif

        call get_overlap_idx(ns,ns3,pno_cv,oidx2,nidx2)

        d2   => pno_cv(ns3)%d
        t22  => pno_t2(ns3)%elm1
        idx2 => pno_cv(ns3)%iaos
        pnv2 =  pno_cv(ns3)%ns2
        pno2 =  pno_cv(ns3)%n

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

      call do_overlap_trafo(ns,ns2,1,pno_S, pnv,nidx1*pnv1*pno1, pnv1 ,w1,w2,ptr=h1,ptr2=h2)
    
      call dgemm('n','n',pnv*nidx1*pnv1,pno,pno1,m10,h1,pnv*nidx1*pnv1,w4,pno1,nul,h2,pnv*nidx1*pnv1)
      call array_reorder_4d(p10,h2,pnv,nidx1,pnv1,pno,[3,4,1,2], nul, h1)

      !transform b index to the correct space
      call do_overlap_trafo(ns,ns2,1,pno_S, pnv,pno*pnv*nidx1,pnv1,h1,h2,ptr=h1)

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
  end subroutine get_common_idx_summation_for_current_aibj
  
  !> Purpose: Wrapper for the RPA model: get MO integrals (non-T1 transformed)
  !
  !> Author:  Pablo Baudin
  !> Date:    January 2014
  subroutine wrapper_to_get_t1_free_gmo(nb,no,nv,Co,Cv,govov,ccmodel,mylsitem)

    implicit none

    integer, intent(in) :: nb, no, nv
    real(realk), pointer, intent(in) :: Co(:,:), Cv(:,:)
    type(array), intent(inout) :: govov
    integer, intent(in) :: ccmodel
    !> LS item with information needed for integrals
    type(lsitem), intent(inout) :: MyLsItem
     
    ! dummy arguments:
    type(array) :: pgmo_diag, pgmo_up
    type(MObatchInfo) :: MOinfo
    logical :: mo_ccsd  
  
    call get_t1_free_gmo(mo_ccsd,mylsitem,Co,Cv,govov,pgmo_diag,pgmo_up, &
                        & nb,no,nv,CCmodel,MOinfo)
 
  end subroutine wrapper_to_get_t1_free_gmo


end module cc_debug_routines_module
