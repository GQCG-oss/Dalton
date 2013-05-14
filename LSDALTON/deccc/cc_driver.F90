!> @file
!> Main CC driver (mp2, cc2 and ccsd are so far implemented)
!> \author Marcin Ziolkowski @ AU 2009,2010 teozio(at)gmail.com

!> General coupled-cluster solver for both full molecular and dec
!> calculations. This should depend only on the integral program and classes
!> related to storage of two- and four-dimensional arrays. All other parameters
!> should be passed as parameters.
module ccdriver

  use precision
  use lstiming!, only: lstimer
  use typedeftype!,only:lsitem
  use typedef
  use files!,only:lsopen,lsclose
  use memory_handling
  use dec_typedef_module
  use integralinterfaceMod
  use integralparameters
#ifdef VAR_LSMPI
  use infpar_module
#endif
  use tensor_interface_module


  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use dec_fragment_utils!,only: get_density_from_occ_orbitals
  use crop
  use array2_simple_operations
  use array4_simple_operations
  use ri_simple_operations!,only: get_ao_ri_intermediate, ri_reset,ri_init, ri_free
  use mp2_module!,only: get_VOVO_integrals
  use ccintegrals!,only:get_full_eri,getL_simple_from_gmo,&
!       & get_gmo_simple,get_h1
  use ccsd_module!,only: getDoublesResidualMP2_simple, &
!       & getDoublesResidualCCSD_simple,getDoublesResidualCCSD_simple2, &
!       & precondition_doubles,get_ccsd_residual_integral_driven,&
!       & get_ccsd_residual_integral_driven_oldarray_wrapper
#ifdef MOD_UNRELEASED
  use ccsdpt_module
!endif mod_unreleased
#endif
  use orbital_operations
  use rpa_module


  ! Interface for CC2 or CCSD energies
  interface get_cc_energy
     module procedure get_cc_energy_arrold
     module procedure get_cc_energy_arrnew
  end interface


  public :: ccsolver, ccsolver_par, fragment_ccsolver, ccsolver_justenergy,&
       & ccsolver_justenergy_pt, mp2_solver
  private

contains

  !> \brief Coupled-cluster solver, works for both full molecule
  !> fragment calculation (depending on the value of fragment_job).
  !> Returns CC energy, converged singles and doubles amplitudes, and
  !> two-electron integrals (a i | b j) stored as (a,i,b,j)
  !> \author Marcin Ziolkowski (modified by Kasper Kristensen and Patrick
  !  Ettenhuber)
  subroutine ccsolver(ypo_f,ypv_f,fock_f,nbasis,nocc,nvirt, &
       & mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f,ccenergy, &
       & t1_final,t2_final,VOVO,longrange_singles)

    implicit none

    !> Number of occupied orbitals in full molecule/fragment AOS
    integer, intent(in) :: nocc
    !> Number of virtual orbitals in full molecule/fragment AOS
    integer, intent(in) :: nvirt
    !> Number of basis functions in full molecule/atomic extent
    integer, intent(in) :: nbasis
    !> Fock matrix in AO basis for fragment or full molecule
    real(realk), dimension(nbasis,nbasis), intent(in) :: fock_f
    !> Occupied MO coefficients for fragment/full molecule
    real(realk), dimension(nbasis,nocc), intent(in) :: ypo_f
    !> Virtual MO coefficients for fragment/full molecule
    real(realk), dimension(nbasis,nvirt), intent(in) :: ypv_f
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
    !> Include long-range singles effects using singles amplitudes
    !> from previous fragment calculations.
    !> IMPORTANT: If this it TRUE, then the singles amplitudes for the fragment
    !> (from previous calculations) must be stored in t1_final at input!
    logical,intent(in) :: longrange_singles
    real(realk),pointer :: yho_d(:,:), yhv_d(:,:),ypo_d(:,:),ypv_d(:,:),focc(:),fvirt(:)
    real(realk),pointer :: ppfock_d(:,:),qqfock_d(:,:), Uocc(:,:), Uvirt(:,:)

    integer, dimension(2) :: occ_dims, virt_dims, ao2_dims, ampl2_dims
    integer, dimension(4) :: ampl4_dims
    type(array2) :: fock,ypo,ypv,yho,yhv
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
    logical :: crop_ok,break_iterations
    type(array4) :: omega2_opt, t2_opt, omega2_prec, u
    type(array2) :: ifock,delta_fock
    type(ri) :: l_ao
    type(array2) :: ppfock_prec, qqfock_prec
    type(array4) :: Lmo
    real(realk) :: tcpu, twall, ttotend_cpu, ttotend_wall, ttotstart_cpu, ttotstart_wall
    real(realk) :: iter_cpu,iter_wall, sosex
    character(18) :: save_to,keep
    character(ARR_MSG_LEN) :: msg
    integer :: ii,aa


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



    ! title
    Call print_ccjob_header(ccPrintLevel,fragment_job,nbasis,nocc,nvirt)

    ! dimension vectors
    occ_dims = [nbasis,nocc]
    virt_dims = [nbasis,nvirt]
    ao2_dims = [nbasis,nbasis]
    ampl4_dims = [nvirt,nocc,nvirt,nocc]
    ampl2_dims = [nvirt,nocc]

    ! go to a (pseudo) canonical basis
    call mem_alloc(focc,nocc)
    call mem_alloc(fvirt,nvirt)
    call mem_alloc(ypo_d,nbasis,nocc)
    call mem_alloc(ypv_d,nbasis,nvirt)
    call mem_alloc(yho_d,nbasis,nocc)
    call mem_alloc(yhv_d,nbasis,nvirt)
    call mem_alloc(ppfock_d,nocc,nocc)
    call mem_alloc(qqfock_d,nvirt,nvirt)
    call mem_alloc(Uocc,nocc,nocc)
    call mem_alloc(Uvirt,nvirt,nvirt)
    if(DECinfo%CCSDpreventcanonical)then
      !nocc diagonalization
      ypo_d   = ypo_f
      ypv_d   = ypv_f
      ppfock_d = ppfock_f
      qqfock_d = qqfock_f
      Uocc=0.0E0_realk
      Uvirt=0.0E0_realk
      do ii=1,nocc
        Uocc(ii,ii) = 1.0E0_realk
      enddo
      do aa=1,nvirt
        Uvirt(aa,aa) = 1.0E0_realk
      enddo
    else
      call get_canonical_integral_transformation_matrices(nocc,nvirt,nbasis,ppfock_f,&
           &qqfock_f,ypo_f,ypv_f,ypo_d,ypv_d,Uocc,Uvirt,focc,fvirt)
      ppfock_d = 0.0E0_realk
      qqfock_d = 0.0E0_realk
      do ii=1,nocc
        ppfock_d(ii,ii) = focc(ii)
      enddo
      do aa=1,nvirt
        qqfock_d(aa,aa) = fvirt(aa)
      enddo
    endif
    call mem_dealloc(focc)
    call mem_dealloc(fvirt)

    ! Copy MO coeffcients. It is very convenient to store them twice to handle transformation
    ! (including transposed MO matrices) efficiently. 
    yho_d = ypo_d
    yhv_d = ypv_d

    ! create transformation matrices in array form
    ypo = array2_init(occ_dims,ypo_d)
    ypv = array2_init(virt_dims,ypv_d)
    yho = array2_init(occ_dims,yho_d)
    yhv = array2_init(virt_dims,yhv_d)
    fock = array2_init(ao2_dims,fock_f)

    call mem_dealloc(ypo_d)
    call mem_dealloc(ypv_d)
    call mem_dealloc(yho_d)
    call mem_dealloc(yhv_d)



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
    call get_density_from_occ_orbitals(nbasis,nocc,ypo%val,dens)

    if((.not. DECinfo%ccsd_old)) then

       if(fragment_job) then ! fragment: calculate correction

          ifock=array2_init(ao2_dims)
          if(longrange_singles) then
             ! Get Fock matrix using singles amplitudes from previous
             ! fragment calculation, thereby effectively including long-range
             ! correlated polarization effects
             call Get_AOt1Fock(mylsitem,t1_final,ifock,nocc,nvirt,nbasis,ypo,yho,yhv)
          else
             ! Fock matrix for fragment from density made from input MOs
             call get_fock_matrix_for_dec(nbasis,dens,mylsitem,ifock,.true.)
          end if

          ! Long range Fock correction:
          delta_fock = getFockCorrection(fock,ifock)
          call array2_free(ifock)

       else 
          ! Full molecule: deltaF = F(Dcore) for frozen core (0 otherwise)
          if(DECinfo%frozencore) then
             ! Fock matrix from input MOs
             ifock=array2_init(ao2_dims)
             call get_fock_matrix_for_dec(nbasis,dens,mylsitem,ifock,.true.)
             ! Correction to actual Fock matrix
             delta_fock = getFockCorrection(fock,ifock)
             call array2_free(ifock)
          else
             delta_fock= array2_init(ao2_dims)
          end if
       end if

    end if
    call mem_dealloc(dens)

    ! get two-electron integrals in ao
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a)') 'debug :: calculating AO integrals'

    ! Only calculate full 4-dimensional AO integrals for old/debug mode
    if(DECinfo%ccsd_old) then
       call get_full_eri(mylsitem,nbasis,gao)
       if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,/)') 'debug :: AO integrals done'
    end if

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
    MP2Special : if(DECinfo%ccModel == 1 .or. DECinfo%ccModel == 5) then

       write(DECinfo%output,*)
       write(DECinfo%output,*) ' ********************  WARNING  **********************'
       write(DECinfo%output,*) 'CCsolver is called for MP2 model.'
       write(DECinfo%output,*) 'This will work fine but it is recommended to use the non-iterative'
       write(DECinfo%output,*) 'MP2_integrals_and_amplitudes_workhorse to use get the MP2 amplitudes'
       write(DECinfo%output,*)
       call get_VOVO_integrals(mylsitem,nbasis,nocc,nvirt,ypv,ypo,gmo)

       ! Construct L: L_{bjai} = 2*g_{bjai} - g_{ajbi}
       Lmo = getL_simple_from_gmo(gmo)

       ppfock = array2_similarity_transformation(ypo,fock,yho,[nocc,nocc])
       qqfock = array2_similarity_transformation(ypv,fock,yhv,[nvirt,nvirt])

       if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,f16.10)') ' debug :: gmo(vovo) norm  : ',gmo*gmo
    end if MP2Special


    ! get fock matrices for preconditioning
    Preconditioner : if(DECinfo%use_preconditioner .or. DECinfo%use_preconditioner_in_b) then
       if(DECinfo%precondition_with_full) then
          ppfock_prec = array2_init([nocc,nocc],ppfock_d)
          qqfock_prec = array2_init([nvirt,nvirt],qqfock_d)
       else
          ppfock_prec = array2_similarity_transformation(ypo,fock,yho,[nocc,nocc])
          qqfock_prec = array2_similarity_transformation(ypv,fock,yhv,[nvirt,nvirt])
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
    if(DECinfo%ccModel /= 1) then
       xocc = array2_init(occ_dims)
       yocc = array2_init(occ_dims)
       xvirt = array2_init(virt_dims)
       yvirt = array2_init(virt_dims)
       if(DECinfo%ccsd_old) then
          h1 = array2_init_plain(ao2_dims)
          CALL II_get_h1_mixed_full(DECinfo%output,DECinfo%output,MyLsitem%SETTING,&
               & h1%val,nbasis,nbasis,AORdefault,AORdefault)
       end if
       if(.not.DECinfo%ccsd_old) iajb =array4_init([nocc,nvirt,nocc,nvirt])
    end if


    call mem_alloc(B,DECinfo%ccMaxIter,DECinfo%ccMaxIter)
    call mem_alloc(c,DECinfo%ccMaxIter)

    ! readme : the iteration sequence is universal and may be used for all
    !          iterative cc models (linear or non-linear) and is
    !          semi-independent on the storage of vectors (allocation and
    !          deallocation, etc)

    ! iterate
    break_iterations = .false.
    crop_ok = .false.
    prev_norm = 1.0E6_realk




    print *
    print *, '### Starting CC iterations'
    print *, '### ----------------------'
    print '(1X,a)',  '###  Iteration     Residual norm          CC energy'

    write(DECinfo%output,*)
    write(DECinfo%output,*) '### Starting CC iterations'
    write(DECinfo%output,*) '### ----------------------'
    write(DECinfo%output,'(1X,a)')  '###  Iteration     Residual norm          CC energy'

    CCIteration : do iter=1,DECinfo%ccMaxIter

       if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)
       call LSTIMER('START',iter_cpu,iter_wall,DECinfo%output)

       ! remove old vectors
       RemoveOldVectors : if(iter > DECinfo%ccMaxDIIS) then

          if(DECinfo%cc_driver_debug) then
             write(DECinfo%output,'(a,i4)') ' debug :: vector to delete : ',iter-DECinfo%ccMaxDIIS
          end if

          if(DECinfo%use_singles) then
             call array2_free(t1(iter-DECinfo%ccMaxDIIS))
             Call array2_free(omega1(iter-DECinfo%ccMaxDIIS))
          end if
          call array4_free(t2(iter-DECinfo%ccMaxDIIS))
          call array4_free(omega2(iter-DECinfo%ccMaxDIIS))
       end if RemoveOldVectors

       ! get new amplitude vectors
       GetGuessVectors : if(iter == 1) then
          if(DECinfo%use_singles) t1(iter) = array2_init(ampl2_dims)
          if(DECinfo%array4OnFile) then
             ! Initialize t2(iter) using storing type 2
             t2(iter) = array4_init(ampl4_dims,2,.true.)
          else
             t2(iter) = array4_init(ampl4_dims)
          end if
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
          call getT1transformation(t1(iter),xocc,xvirt,yocc,yvirt, &
               ypo,ypv,yho,yhv)

          ! get inactive fock
          if(DECinfo%fock_with_ri) then
             ! Debug mode
             ifock = getInactiveFockFromRI(l_ao,xocc,yocc,h1)
          else
             if(DECinfo%ccsd_old) ifock = getInactiveFock_simple(h1,gao,xocc,yocc,nocc,nbasis)
          end if
          ! Note: If not fock_with_ri or ccsd_old, then the relevant
          ! ifock is calculated below in get_ccsd_residual_integral_direct.
          ! Long range fock matrix correction using old scheme
          ! (See comments above regarding Fock correction)
          if(DECinfo%ccsd_old) then

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
          end if

          ! readme : this should be done in a more clear way
          if(DECinfo%ccModel == 2 .and. DECinfo%ccsd_old) then
             ! CC2
             ppfock = array2_similarity_transformation(xocc,fock,yocc,[nocc,nocc])
             qqfock = array2_similarity_transformation(xvirt,fock,yvirt,[nvirt,nvirt])
          else if(DECinfo%ccModel >= 3 .and. DECinfo%ccsd_old) then
             ! CCSD
             ppfock = array2_similarity_transformation(xocc,ifock,yocc,[nocc,nocc])
             qqfock = array2_similarity_transformation(xvirt,ifock,yvirt,[nvirt,nvirt])
          endif

          if(DECinfo%ccsd_old) then
             pqfock = array2_similarity_transformation(xocc,ifock,yvirt,[nocc,nvirt])
             qpfock = array2_similarity_transformation(xvirt,ifock,yocc,[nvirt,nocc])
             iajb = get_gmo_simple(gao,xocc,yvirt,xocc,yvirt)
          else
             ppfock=array2_init([nocc,nocc])
             pqfock=array2_init([nocc,nvirt])
             qpfock=array2_init([nvirt,nocc])
             qqfock=array2_init([nvirt,nvirt])
          endif

       end if T1Related

       if(DECinfo%PL>1) call LSTIMER('CCIT: INIT',tcpu,twall,DECinfo%output)
       if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)


       !Write the vectors to file and keep only current in memory
!!$       if (iter>1) then
!!$         write (save_to,'("omega2_iter",i2.2,".data")') iter-1
!!$         call dump_array4_to_file(omega2(iter-1),save_to)
!!$         write (save_to,'("t2_iter",i2.2,".data")') iter-1
!!$         call dump_array4_to_file(t2(iter-1),save_to)
!!$         k = min(iter,DECinfo%ccMaxDIIS)
!!$         do i=iter-k+1,iter-2,1
!!$           call memory_deallocate_4d(omega2(i)%val)
!!$           call memory_deallocate_4d(t2(i)%val)
!!$           !call array4_free(omega2(i))
!!$           !call array4_free(t2(i))
!!$         enddo
!!$       endif


       ! MODIFY FOR NEW MODEL
       ! If you implement a new model, please insert call to your own residual routine here!
       SelectCoupledClusterModel : if(DECinfo%ccModel==1) then
          call getDoublesResidualMP2_simple(Omega2(iter),t2(iter),gmo,ppfock,qqfock, &
               & nocc,nvirt)
       elseif(DECinfo%ccModel==2) then
          if (DECinfo%ccsd_old) then
             !print *,"old cc2 scheme"
             u = get_u(t2(iter))
             call getSinglesResidualCCSD(omega1(iter),u,gao,pqfock,qpfock, &
                  xocc,xvirt,yocc,yvirt,nocc,nvirt)
             !call print_norm(omega1(iter)%val,nvirt,nocc)
             call array4_free(u)

             gmo = get_gmo_simple(gao,xvirt,yocc,xvirt,yocc)
             call getDoublesResidualMP2_simple(omega2(iter),t2(iter),gmo,ppfock,qqfock, &
                  nocc,nvirt)
             !call print_norm(1.0E0_realk,omega2(iter)%val,nvirt*nvirt,nocc*nocc)
             call array4_free(gmo)
          else
             !print *,"new cc2 scheme"
               call get_ccsd_residual_integral_driven_oldarray_wrapper(delta_fock,&
                  & omega2(iter),t2(iter),&
                  & fock,iajb%val,nocc,nvirt,ppfock,qqfock,pqfock,qpfock,xocc,xvirt,&
                  & yocc,yvirt,nbasis,MyLsItem, omega1(iter),iter)
          endif
       elseif(DECinfo%ccmodel==3 .or. DECinfo%ccmodel==4) then  ! CCSD or CCSD(T)
          if(DECinfo%ccsd_old) then

             u = get_u(t2(iter))
             call getSinglesResidualCCSD(omega1(iter),u,gao,pqfock,qpfock, &
                  xocc,xvirt,yocc,yvirt,nocc,nvirt)
             !aibj = get_gmo_simple(gao,xvirt,yocc,xvirt,yocc)
             aibj = get_gmo_simple(gao,yocc,xvirt,yocc,xvirt)
             call array4_reorder(aibj,[2,1,4,3])

             if (DECinfo%ccsd_expl) then
               !print *,"calling the old ccsd scheme"
               call getDoublesResidualCCSD_simple(omega2(iter),t2(iter),u,gao,aibj,iajb,nocc,nvirt, &
                    ppfock,qqfock,xocc,xvirt,yocc,yvirt)
               ! just to debug
               !          call getDoublesResidual_explicite(omega2(iter),t2(iter),u,gao,aibj,iajb,ppfock,qqfock, &
               !            xocc,xvirt,yocc,yvirt,nocc,nvirt,nbasis)
             else
               call getDoublesResidualCCSD_simple2(omega2(iter),t2(iter),u,gao,aibj,iajb,nocc,nvirt, &
                    ppfock,qqfock,xocc,xvirt,yocc,yvirt,nbasis)
             endif
          else
             call get_ccsd_residual_integral_driven_oldarray_wrapper(delta_fock,omega2(iter),&
                & t2(iter),&
                & fock,iajb%val,nocc,nvirt,ppfock,qqfock,pqfock,qpfock,xocc,xvirt,&
                & yocc,yvirt,nbasis,MyLsItem, omega1(iter),iter)
          end if
          if (DECinfo%ccsd_old) then
             call array4_free(aibj)
             call array4_free(u)
          endif

       elseif(DECinfo%ccmodel==5) then
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
       do i=iter,max(iter-DECinfo%ccMaxDIIS+1,1),-1
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
       call CalculateDIIScoefficients(DECinfo%ccMaxDIIS,DECinfo%ccMaxIter,iter,B,c, &
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

       do i=iter,max(iter-DECinfo%ccMaxDIIS+1,1),-1
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
       EnergyForCCmodel: if(DECinfo%ccmodel==1) then  
          ! MP2
          ccenergy = get_mp2_energy(t2(iter),Lmo)
       elseif(DECinfo%ccmodel==2 .or. DECinfo%ccmodel==3 .or. DECinfo%ccmodel==4 ) then
          ! CC2, CCSD, or CCSD(T) (for (T) calculate CCSD contribution here)
          ccenergy = get_cc_energy(t1(iter),t2(iter),iajb,nocc,nvirt)
       elseif(DECinfo%ccmodel==5) then
          ccenergy = RPA_energy(t2(iter),gmo)
          sosex = SOSEX_contribution(t2(iter),gmo)
          ccenergy=ccenergy+sosex
       end if EnergyForCCmodel


       if(DECinfo%PL>1) call LSTIMER('CCIT: ENERGY',tcpu,twall,DECinfo%output)

       ! check if this is the last iteration
       if(iter == DECinfo%ccMaxIter .or. two_norm_total < DECinfo%ccConvergenceThreshold) &
            break_iterations=.true.

       if(DECinfo%use_singles .and. (.not. break_iterations) ) then
         if(DECinfo%ccsd_old)call array4_free(iajb)
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
       end if
       !msg="t1+1:"
       !call print_norm(t1(iter+1)%val,nocc*nvirt,msg) 
       !msg="t2+1:"
       !call print_norm(t2(iter+1)%val,nocc*nocc*nvirt*nvirt,msg) 

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

      print '(1X,a,2X,i4,5X,g19.9,4X,g19.9)',  '### ',iter, two_norm_total,ccenergy
      write(DECinfo%output,'(1X,a,2X,i4,5X,g19.9,4X,g19.9)') &
            &   '### ',iter, two_norm_total,ccenergy
       last_iter = iter
       if(break_iterations) exit

    end do CCIteration

    call LSTIMER('START',ttotend_cpu,ttotend_wall,DECinfo%output)



    write(DECinfo%output,*)
    write(DECinfo%output,'(/,a)') '-------------------------------'
    write(DECinfo%output,'(a)')   '  Coupled-cluster job summary  '
    write(DECinfo%output,'(a,/)') '-------------------------------'
    if(break_iterations) then
       write(DECinfo%output,'(a)')     'Hooray! CC equation is solved!'
    else
       write(DECinfo%output,'(a,i4,a)')  'CC equation not solved in ', &
            & DECinfo%ccMaxIter, ' iterations!'
       call lsquit('CC equation not solved!',DECinfo%output)
    end if
    write(DECinfo%output,'(a,f16.3,a)') 'CCSOL: Total cpu time    = ',ttotend_cpu-ttotstart_cpu,' s'
    write(DECinfo%output,'(a,f16.3,a)') 'CCSOL: Total wall time   = ',ttotend_wall-ttotstart_wall,' s'

    if(fragment_job) then
       write(DECinfo%output,'(a,f16.10)')  'Frag. corr. energy = ',ccenergy
    else
       write(DECinfo%output,'(a,f16.10)')  'Corr. energy       = ',ccenergy
    end if
    write(DECinfo%output,'(a,i5)') 'Number of CC iterations =', last_iter


    ! Free memory and save final amplitudes
    ! *************************************


    ! remove rest of the singles amplitudes and residuals
    do i=last_iter,max(last_iter-DECinfo%ccMaxDIIS+1,1),-1


       ! remove the lase files of t2 and omega2
       call array4_delete_file(omega2(i))
       call array4_delete_file(t2(i))

       if(DECinfo%use_singles) then

          ! Save final singles amplitudes
          if(i==last_iter) then
             if(longrange_singles) then ! just copy
                call array2_copy(t1_final,t1(last_iter))
             else ! initialize and copy
                t1_final = array2_duplicate(t1(last_iter))
             end if
          end if

          ! Free singles amplitudes and residuals
          call array2_free(t1(i))
          call array2_free(omega1(i))

       end if

       ! Free doubles residuals
       call array4_free(omega2(i))

       ! Save final double amplitudes
       if(i==last_iter) then
          t2_final = array4_duplicate(t2(last_iter))
       end if

       ! Free doubles amplitudes
       call array4_free(t2(i))

    end do


    ! Save two-electron integrals in the order (virt,occ,virt,occ)
    if(DECinfo%ccModel == 1) then
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

    if(DECinfo%ccsd_expl .or. DECinfo%ccsd_old) then
       call array4_free(gao)
    endif

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

    call array2_free(ypo)
    call array2_free(yho)
    call array2_free(ypv)
    call array2_free(yhv)
    call array2_free(fock)
    !transform back to original basis   
    if(DECinfo%use_singles)then
      call ccsolver_can_local_trans(VOVO%val,t2_final%val,nocc,nvirt,Uocc,Uvirt,t1_final%val)
    else
      call ccsolver_can_local_trans(VOVO%val,t2_final%val,nocc,nvirt,Uocc,Uvirt)
    endif

    call mem_dealloc(Uocc)
    call mem_dealloc(Uvirt)


  end subroutine ccsolver

  !> \brief Get coupled-cluster energy by calling general ccsolver.
  !> \author Kasper Kristensen
  function ccsolver_justenergy(ypo_f,ypv_f,fock_f,nbasis,nocc,nvirt, &
       &mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f) result(ccenergy)

    implicit none

    !> Number of occupied orbitals in full molecule/fragment AOS
    integer, intent(in) :: nocc
    !> Number of virtual orbitals in full molecule/fragment AOS
    integer, intent(in) :: nvirt
    !> Number of basis functions in full molecule/atomic extent
    integer, intent(in) :: nbasis
    !> Fock matrix in AO basis for fragment or full molecule
    real(realk), dimension(nbasis,nbasis), intent(in) :: fock_f
    !> Occupied MO coefficients  for fragment/full molecule
    real(realk), dimension(nbasis,nocc), intent(inout) :: ypo_f
    !> Virtual MO coefficients  for fragment/full molecule
    real(realk), dimension(nbasis,nvirt), intent(inout) :: ypv_f
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
    real(realk) :: ccenergy!,ccsdpt_e4,ccsdpt_e5,ccsdpt_tot
    type(array4) :: t2_final,VOVO!,ccsdpt_t2
    type(array2) :: t1_final!,ccsdpt_t1

    if(.not. DECinfo%solver_par .or. DECinfo%ccModel==1)then
      call ccsolver(ypo_f,ypv_f,fock_f,nbasis,nocc,nvirt, &
         & mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f,ccenergy, &
         & t1_final,t2_final,VOVO,.false.)
    else
      call ccsolver_par(ypo_f,ypv_f,fock_f,nbasis,nocc,nvirt, &
         & mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f,ccenergy, &
         & t1_final,t2_final,VOVO,.false.)
    endif

    ! Free arrays
    call array2_free(t1_final)
    call array4_free(t2_final)
    call array4_free(VOVO)

  end function ccsolver_justenergy

#ifdef MOD_UNRELEASED

  !> \brief get ccsd(t) corrections for full molecule.
  !> \author Janus Juul Eriksen
  !> \date February 2013
  function ccsolver_justenergy_pt(MyMolecule,nbasis,nocc,nvirt,mylsitem,&
                      & ccPrintLevel,fragment_job,ypo_fc,ppfock_fc) result(ccenergy)

    implicit none

    !> full molecule information
    type(fullmolecule), intent(in) :: MyMolecule
    !> Number of occupied orbitals in full molecule/fragment AOS
    integer, intent(in) :: nocc
    !> Number of virtual orbitals in full molecule/fragment AOS
    integer, intent(in) :: nvirt
    !> Number of basis functions in full molecule/atomic extent
    integer, intent(in) :: nbasis
    !> Is this a fragment job (true) or a full molecular calculation (false)
    logical, intent(in) :: fragment_job
    !> LS item information
    type(lsitem), intent(inout) :: mylsitem
    !> How much to print? ( ccPrintLevel>0 --> print info stuff)
    integer, intent(in) :: ccPrintLevel
    !> Occupied MO coefficients  for fragment/full molecule (only used for Frozen core)
    real(realk), dimension(nbasis,nocc), intent(in),optional :: ypo_fc
    !> Occ-occ block of Fock matrix in MO basis (only used for frozen core)
    real(realk), dimension(nocc,nocc), intent(in),optional :: ppfock_fc
    !> Coupled cluster energy for full molecule
    real(realk) :: ccenergy,ccsdpt_e4,ccsdpt_e5,ccsdpt_tot
    type(array4) :: t2_final,ccsdpt_t2,VOVO
    type(array2) :: t1_final,ccsdpt_t1,ccsd_mat_tot,ccsd_mat_tmp,e4_mat_tot,e4_mat_tmp,e5_mat_tot
    integer :: natoms,ncore,nocc_tot,p,pdx,i
    !> stuff needed for pair analysis
    real(realk), pointer :: distance_table(:,:)
    type(ccorbital), pointer :: occ_orbitals(:)
    type(ccorbital), pointer :: unocc_orbitals(:)
    logical, pointer :: orbitals_assigned(:)

    ! is this a frozen core calculation or not?
    if (DECinfo%frozencore) then
       ncore = MyMolecule%ncore
       if(.not. present(ypo_fc)) then
          call lsquit('ccsolver_justenergy_pt: Occ MOs not present for frozencore!',-1)
       end if
       if(.not. present(ppfock_fc)) then
          call lsquit('ccsolver_justenergy_pt: Occ-occ Fock matrix not present for frozencore!',-1)
       end if

       if (.not. DECinfo%solver_par) then
          call ccsolver(ypo_fc,MyMolecule%ypv,MyMolecule%fock,nbasis,nocc,nvirt,&
               & mylsitem,ccPrintLevel,fragment_job,ppfock_fc,MyMolecule%qqfock,ccenergy,&
               & t1_final,t2_final,VOVO,.false.)
       else
          call ccsolver_par(ypo_fc,MyMolecule%ypv,MyMolecule%fock,nbasis,nocc,nvirt,&
               & mylsitem,ccPrintLevel,fragment_job,ppfock_fc,MyMolecule%qqfock,ccenergy,&
               & t1_final,t2_final,VOVO,.false.)
       endif

    else
       ncore = 0

       if (.not. DECinfo%solver_par) then
          call ccsolver(MyMolecule%ypo,MyMolecule%ypv,MyMolecule%fock,nbasis,nocc,nvirt,&
               & mylsitem,ccPrintLevel,fragment_job,MyMolecule%ppfock,MyMolecule%qqfock,ccenergy,&
               & t1_final,t2_final,VOVO,.false.)
       else
          call ccsolver_par(MyMolecule%ypo,MyMolecule%ypv,MyMolecule%fock,nbasis,nocc,nvirt,&
               & mylsitem,ccPrintLevel,fragment_job,MyMolecule%ppfock,MyMolecule%qqfock,ccenergy,&
               & t1_final,t2_final,VOVO,.false.)
       endif

    end if


    ! do ccsd

!    ! free integrals
!    call array4_free(VOVO)

    natoms = MyMolecule%natoms
    nocc_tot = MyMolecule%numocc

    ccsdpt_t1 = array2_init([nvirt,nocc])
    ccsdpt_t2 = array4_init([nvirt,nvirt,nocc,nocc])

    if(DECinfo%frozencore) then
       call ccsdpt_driver(nocc,nvirt,nbasis,ppfock_fc,MyMolecule%qqfock,ypo_fc,MyMolecule%ypv,mylsitem,t2_final,&
            & ccsdpt_t1,ccsdpt_t2)
    else
       call ccsdpt_driver(nocc,nvirt,nbasis,MyMolecule%ppfock,MyMolecule%qqfock,MyMolecule%ypo,&
            & MyMolecule%ypv,mylsitem,t2_final,ccsdpt_t1,ccsdpt_t2)
    end if


    ! as we want to  print out fragment and pair interaction fourth-order energy contributions,
    ! then for locality analysis purposes we need distance_table and occ_orbitals/
    ! unocc_orbitals (adapted from fragment_energy.f90)

    ! -- Calculate distance matrix
    call mem_alloc(distance_table,natoms,natoms)
    distance_table = 0.0E0_realk
    call GetDistances(distance_table,natoms,mylsitem,DECinfo%output) ! distances in atomic units

    ! -- Analyze basis and create orbitals
    call mem_alloc(occ_orbitals,nocc_tot)
    call mem_alloc(unocc_orbitals,nvirt)
    call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc_tot,nvirt,natoms, &
         & occ_orbitals,unocc_orbitals,distance_table)

    ! Orbital assignment
    call mem_alloc(orbitals_assigned,natoms)
    orbitals_assigned=.false.
    do p=1,nocc_tot
       pdx = occ_orbitals(p)%centralatom
       orbitals_assigned(pdx) = .true.
    end do
    do p=1,nvirt
       pdx = unocc_orbitals(p)%centralatom
       orbitals_assigned(pdx) = .true.
    end do

    ! reorder VOVO integrals from (a,i,b,j) to (a,b,i,j)
    call array4_reorder(VOVO,[1,3,2,4])

    ! print out ccsd fragment and pair interaction energies
    ccsd_mat_tot = array2_init([natoms,natoms])
    ccsd_mat_tmp = array2_init([natoms,natoms])

    call ccsd_energy_full(nocc,nvirt,natoms,ncore,t2_final,t1_final,VOVO,occ_orbitals,&
                           & ccsd_mat_tot%val,ccsd_mat_tmp%val)

    call print_ccsd_full(natoms,ccsd_mat_tot%val,orbitals_assigned,distance_table)

    ! release ccsd stuff
    call array2_free(ccsd_mat_tot)
    call array2_free(ccsd_mat_tmp)

    ! free integrals
    call array4_free(VOVO)

    ! now we calculate fourth-order (which are printed out in print_e4_full) and fifth-order energies
    e4_mat_tot = array2_init([natoms,natoms])
    e4_mat_tmp = array2_init([natoms,natoms])
    e5_mat_tot = array2_init([natoms,natoms])

    call ccsdpt_energy_e4_full(nocc,nvirt,natoms,ncore,t2_final,ccsdpt_t2,occ_orbitals,&
                             & e4_mat_tot%val,e4_mat_tmp%val,ccsdpt_e4)

    call ccsdpt_energy_e5_full(nocc,nvirt,natoms,ncore,t1_final,ccsdpt_t1,&
                             & occ_orbitals,unocc_orbitals,e5_mat_tot%val,ccsdpt_e5)

    ! print out the fourth- and fifth-order fragment and pair interactin energies
    call print_e4_full(natoms,e4_mat_tot%val,orbitals_assigned,distance_table)

    call print_e5_full(natoms,e5_mat_tot%val,orbitals_assigned,distance_table)

    ! release stuff
    call array2_free(e4_mat_tot)
    call array2_free(e4_mat_tmp)
    call array2_free(e5_mat_tot)
    call mem_dealloc(distance_table)
    do i=1,nocc_tot
       call orbital_free(occ_orbitals(i))
    end do
    call mem_dealloc(occ_orbitals)
    do i=1,nvirt
       call orbital_free(unocc_orbitals(i))
    end do
    call mem_dealloc(unocc_orbitals)
    call mem_dealloc(orbitals_assigned)

    ! sum up energies
    ccsdpt_tot = ccsdpt_e4 + ccsdpt_e5

    ! finally, print out total energies 
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '*****************************************************************************'
    write(DECinfo%output,'(1X,a)') '*                      Full CCSD(T) calculation is done !                   *'
    write(DECinfo%output,'(1X,a)') '*****************************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,g20.10)') 'Total CCSD correlation energy           =', ccenergy
    write(DECinfo%output,*)

    ! now update ccenergy with ccsd(t) correction
    ccenergy = ccenergy + ccsdpt_tot

    write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,g20.10)') 'The E4 doubles and triples contribution =', ccsdpt_e4
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,g20.10)') 'The E5 singles and triples contribution =', ccsdpt_e5
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,g20.10)') 'Total CCSD(T) energy contribution       =', ccsdpt_tot
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,g20.10)') 'Total CCSD(T) correlation energy        =', ccenergy
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    ! free amplitude arrays
    call array2_free(t1_final)
    call array4_free(t2_final)
    call array2_free(ccsdpt_t1)
    call array4_free(ccsdpt_t2)

  end function ccsolver_justenergy_pt
!endif mod_unreleased
#endif

  !> \brief For a given fragment, calculate singles and doubles amplitudes and
  !> two-electron integrals (a i | bj ) required for CC energy.
  !> Intended to be used for CC2 and CCSD (and NOT for MP2).
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine fragment_ccsolver(MyFragment,t1,t2,VOVO)

    implicit none

    !> Fragment info (only t1 information in MyFragment may be changed here)
    type(ccatom), intent(inout) :: MyFragment
    !> Singles amplitudes t1(a,i)
    type(array2),intent(inout) :: t1
    !> Doubles amplitudes t2(a,i,b,j)
    type(array4),intent(inout) :: t2
    !> Two electron integrals (a i | b j) stored as (a,i,b,j)
    type(array4),intent(inout) :: VOVO
    integer :: dims(2)
    real(realk) :: ccenergy

    ! Sanity check: This routine is not intended for MP2
    if(DECinfo%ccmodel == 1) then
       call lsquit('fragment_ccsolver cannot be used for MP2!',&
            & DECinfo%output)
    end if

    ! If MyFragment%t1_stored is TRUE, then we reuse the singles amplitudes
    ! from previous fragment calculations to describe long-range
    ! singles effects.
    ! In this case the fragment t1 amplitudes are stored in MyFragment%t1
    if(MyFragment%t1_stored) then
       dims(1) = MyFragment%nunoccAOS
       dims(2) = MyFragment%noccAOS
       t1 = array2_init(dims,MyFragment%t1)
    end if

    if(DECinfo%solver_par)then
      call ccsolver_par(myfragment%ypo,myfragment%ypv,&
         & myfragment%fock, myfragment%number_basis,myfragment%noccAOS,&
         & myfragment%nunoccAOS,myfragment%mylsitem,DECinfo%PL,&
         & .true.,myfragment%ppfock,myfragment%qqfock,ccenergy,t1,t2,VOVO,MyFragment%t1_stored)
    else
      call ccsolver(myfragment%ypo,myfragment%ypv,&
         & myfragment%fock, myfragment%number_basis,myfragment%noccAOS,&
         & myfragment%nunoccAOS,myfragment%mylsitem,DECinfo%PL,&
         & .true.,myfragment%ppfock,myfragment%qqfock,ccenergy,t1,t2,VOVO,MyFragment%t1_stored)
    endif

    ! Save singles amplitudes in fragment structure
    if(DECinfo%SinglesPolari) then
       call save_fragment_t1_AOSAOSamplitudes(MyFragment,t1)
    end if

  end subroutine fragment_ccsolver



  !> \brief For a given fragment, calculate combined doubles+singles amplitudes:
  !> u(a,i,b,j) = t2(a,i,b,j) + t1(a,i)*t1(b,j)
  !> and two-electron integrals (a i | bj ).
  !> IMPORTANT: The EOS indices are extracted for both the occupied and
  !> the virtual spaces, such that the output amplitudes/integrals
  !> can be used directly in the DEC hybrid scheme for determining the
  !> individual orbital contributions to the fragment energy.
  !> Intended to be used for CC2 and CCSD (and NOT for MP2).
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine CCfragment_get_EOS_amplitudes_and_integrals(MyFragment,&
       & uocc,uvirt,VOVOocc,VOVOvirt)

    implicit none

    !> Fragment info (only t1 information in MyFragment may be changed here)
    type(ccatom), intent(inout) :: MyFragment
    !> Combined doubles+singles amplitudes for occupied partitioning
    type(array4),intent(inout) :: uocc
    !> Combined doubles+singles amplitudes for virtual partitioning
    type(array4),intent(inout) :: uvirt
    !> Two electron integrals (a i | b j) stored as (a,i,b,j) for occ part.
    !> (a,b: AOS orbitals;   i,j: EOS orbitals)
    type(array4),intent(inout) :: VOVOocc
    !> Two electron integrals (a i | b j) stored as (a,i,b,j) for virt part.
    !> (a,b: EOS orbitals;   i,j: AOS orbitals)
    type(array4),intent(inout) :: VOVOvirt
    type(array2) :: t1
    type(array4) :: t2,u,VOVO


    ! Solve CC equation to calculate amplitudes and integrals
    ! *******************************************************
    ! Here all output indices in t1,t2, and VOVO are AOS indices.
    call fragment_ccsolver(MyFragment,t1,t2,VOVO)


    ! Extract EOS indices for integrals
    ! *********************************
    call array4_extract_eos_indices_both_schemes(VOVO, &
         & VOVOocc, VOVOvirt, MyFragment)
    call array4_free(VOVO)


    ! Calculate combined single+doubles amplitudes
    ! ********************************************
    ! u(a,i,b,j) = t2(a,i,b,j) + t1(a,i)*t1(b,j)
    call get_combined_SingleDouble_amplitudes(t1,t2,u)
    call array2_free(t1)
    call array4_free(t2)


    ! Extract EOS indices for amplitudes
    ! **********************************
    call array4_extract_eos_indices_both_schemes(u, &
         & uocc, uvirt, MyFragment)
    call array4_free(u)


  end subroutine CCfragment_get_EOS_amplitudes_and_integrals



  !> \brief Solve MP2 equation:
  !> RHS_{bjai} =  - sum_{c} t_{bjci} F_{ca}
  !>               - sum_{c} t_{aicj} F_{cb}
  !>               + sum_{k} t_{bjak} F_{ki}
  !>               + sum_{k} t_{aibk} F_{kj}
  !> It is assumed that RHS_{bjai} = R_{aibj} !
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine mp2_solver(nocc,nvirt,ppfock,qqfock,RHS,t2)

    implicit none
    !> Number of occupied orbitals (dimension of ppfock)
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals (dimension of qqfock)
    integer, intent(in) :: nvirt
    !> Occupied-occupied block of Fock matrix
    real(realk) :: ppfock(nocc,nocc)
    !> Unoccupied-unoccupied block of Fock matrix
    real(realk) :: qqfock(nvirt,nvirt)
    !> RHS array
    type(array4), intent(in) :: RHS
    !> Solution array
    type(array4), intent(inout) :: t2
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    if(DECinfo%array4OnFile) then ! RHS and t2 values are stored on file
       call mp2_solver_file(nocc,nvirt,ppfock,qqfock,RHS,t2)
    else ! RHS and t2 values are stored in memory
       call mp2_solver_mem(nocc,nvirt,ppfock,qqfock,RHS,t2)
    end if

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)


  end subroutine mp2_solver



  !> \brief Solve MP2 equation when RHS and t2 are values are stored in memory.
  !> See mp2_solver for details about the equation.
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine mp2_solver_mem(nocc,nvirt,ppfock,qqfock,RHS,t2)

    implicit none
    !> Number of occupied orbitals (dimension of ppfock)
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals (dimension of qqfock)
    integer, intent(in) :: nvirt
    !> Occupied-occupied block of Fock matrix
    real(realk) :: ppfock(nocc,nocc)
    !> Unoccupied-unoccupied block of Fock matrix
    real(realk) :: qqfock(nvirt,nvirt)
    !> RHS array
    type(array4), intent(in) :: RHS
    !> Solution array
    type(array4), intent(inout) :: t2
    type(array4) :: tmp1,tmp2
    real(realk),pointer :: Cocc_data(:,:), Cvirt_data(:,:), Socc(:,:), Svirt(:,:)
    real(realk),pointer :: EVocc(:), EVvirt(:)
    type(array2) :: Cocc, Cvirt
    integer :: I,J,A,B
    real(realk) :: tcpu, twall, deltaF
    integer :: dims(4), occdims(2), virtdims(2)
    ! real(realk) :: test


    ! Strategy for solving MP2 equation:
    ! 1. Find basis where Fock matrix is diagonal
    ! 2. Transform 2-electron integrals to diagonal basis
    ! 3. In diagonal basis the solution is trivial and the amplitudes are found.
    ! 4. Transform amplitudes in diagonal basis back to LCM basis.

    call LSTIMER('START',tcpu,twall,DECinfo%output)


    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Entering MP2 solver - store array values in memory'
    write(DECinfo%output,*)


    ! Sanity checks
    ! *************
    ! Check that nvirt /= 0.
    if(nvirt<1 .or. nocc<1) then
       write(DECinfo%output,*) 'Number of occupied orbitals = ', nocc
       write(DECinfo%output,*) 'Number of unoccupied orbitals = ', nvirt
       call lsquit('Error in mp2_solver: Number of orbitals is smaller than one!', DECinfo%output)
    endif

    ! Only implemented for MP2
    if(DECinfo%ccModel /= 1) then
       call lsquit('mp2_solver called with other model than MP2',DECinfo%output)
    end if




    ! Initialize stuff
    ! ****************
    dims = [nvirt,nocc,nvirt,nocc]
    occdims = [nocc,nocc]
    virtdims = [nvirt,nvirt]




    ! 1. Solve Fock eigenvalue problem - each block separately
    ! ********************************************************

    ! OCCUPIED-OCCUPIED BLOCK
    ! '''''''''''''''''''''''

    ! Eigenvectors
    call mem_alloc(Cocc_data,nocc,nocc)

    ! Eigenvalues
    call mem_alloc(EVocc,nocc)

    ! The overlap matrix is simply the unit matrix because
    ! the LCM/MLM basis is orthogonal.
    call mem_alloc(Socc,nocc,nocc)
    Socc=0.0E0_realk
    do i=1,nocc
       Socc(i,i) = 1E0_realk
    end do

    ! Solve eigenvalue problem
    call solve_eigenvalue_problem(nocc,ppfock,Socc,EVocc,Cocc_data)

    ! For later, it is convenient to keep eigenvectors in array2 form
    Cocc = array2_init(occdims,Cocc_data)

    ! Done with some matrices
    call mem_dealloc(Cocc_data)
    call mem_dealloc(Socc)




    ! VIRTUAL-VIRTUAL BLOCK
    ! '''''''''''''''''''''

    ! Eigenvectors
    call mem_alloc(Cvirt_data,nvirt,nvirt)

    ! Eigenvalues
    call mem_alloc(EVvirt,nvirt)

    ! Unit overlap for virtual space
    call mem_alloc(Svirt,nvirt,nvirt)
    Svirt=0.0E0_realk
    do i=1,nvirt
       Svirt(i,i) = 1E0_realk
    end do

    ! Solve eigenvalue problem
    call solve_eigenvalue_problem(nvirt,qqfock,Svirt,EVvirt,Cvirt_data)

    ! For later, it is convenient to keep eigenvectors in array2 form
    Cvirt = array2_init(virtdims,Cvirt_data)

    ! Done with some matrices
    call mem_dealloc(Cvirt_data)
    call mem_dealloc(Svirt)


    call LSTIMER('SOLVE: EIVAL',tcpu,twall,DECinfo%output)



    ! 2. Transform two-electron integrals to diagonal basis
    ! *****************************************************

    ! Using notation that (a,i,b,j) are LCM indices
    ! and (A,I,B,J) are indices in the diagonal basis,
    ! we want to carry out the transformations:
    ! RHS_{AIBJ} = sum_{aibj} C_{aA} C_{iI} C_{bB} C_{jJ} RHS_{aibj} (*)

    ! 1. Init temporary array
    tmp1= array4_init(dims)

    ! 2. Index A: RHS(a,i,b,j) --> tmp1(A,i,b,j)
    call array4_contract1(RHS,Cvirt,tmp1,.true.)

    ! 3. Index I: tmp1(A,i,b,j) --> tmp1(i,A,b,j) --> tmp2(I,A,b,j)
    call array4_reorder(tmp1,[2,1,3,4])
    tmp2= array4_init([nocc,nvirt,nvirt,nocc])

    call array4_contract1(tmp1,Cocc,tmp2,.true.)
    call array4_free(tmp1)

    ! 4. Index J: tmp2(I,A,b,j) --> tmp2(j,b,A,I) --> tmp1(J,b,A,I)
    call array4_reorder(tmp2,[4,3,2,1])
    tmp1= array4_init([nocc,nvirt,nvirt,nocc])
    call array4_contract1(tmp2,Cocc,tmp1,.true.)
    call array4_free(tmp2)

    ! 5. Index B: tmp1(J,b,A,I) --> tmp1(b,J,A,I) --> t2(B,J,A,I)
    call array4_reorder(tmp1,[2,1,3,4])
    t2 = array4_init(dims)
    call array4_contract1(tmp1,Cvirt,t2,.true.)
    call array4_free(tmp1)

    ! Now t2 contains the two-electron integrals (AI|BJ) = (BJ|AI) in the order (B,J,A,I).
    ! [Due to the symmetry it is not necessary to reorder t2 back to (A,I,B,J)].


    call LSTIMER('SOLVE: TRANS 1',tcpu,twall,DECinfo%output)



    ! 3. Solve MP2 equation in the diagonal basis
    ! *******************************************

    ! In the diagonal basis the solution to the MP2 equation is trivial!
    ! The equation:
    ! RHS_{bjai} =  - sum_{c} t_{bjci} F_{ca}
    !               - sum_{c} t_{aicj} F_{cb}
    !               + sum_{k} t_{bjak} F_{ki}
    !               + sum_{k} t_{aibk} F_{kj}
    !
    ! simply becomes (using t_{BJAI} = t_{AIBJ}):
    !
    ! RHS_{BJAI} =  - t_{BJAI} F_{AA}
    !               - t_{AIBJ} F_{BB}
    !               + t_{BJAI} F_{II}
    !               + t_{AIBJ} F_{JJ}
    !            =  t_{BJAI} [ F_{II} + F_{JJ} - F_{AA} - F_{BB} ]
    !
    ! In other words:
    !
    ! t_{BJAI} = RHS_{BJAI} / [ F_{II} + F_{JJ} - F_{AA} - F_{BB} ] (**)
    !

    ! Recalling that currently t_{BJAI} = RHS_{BJAI} and that the
    ! Fock matrix in the diagonal basis are the eigenvalues EVocc and EVvirt -
    ! we simply modify the t2 elements as follows:

    !    test=0E0_realk
    do I=1,nocc
       do A=1,nvirt
          do J=1,nocc
             do B=1,nvirt

                ! deltaF = F_{II} + F_{JJ} - F_{AA} - F_{BB}
                deltaF = EVocc(I) + EVocc(J) - EVvirt(A) - EVvirt(B)

                ! Sanity check
                if( abs(deltaF) < 1e-9 ) then
                   write(DECinfo%output,*) 'WARNING: SMALL NUMBERS OCCURING IN SOLVER!!!'
                   write(DECinfo%output,*) 'WARNING: SOLVER MAY BE UNSTABLE!!!'
                   write(DECinfo%output,*) 'Delta epsilon value = ', deltaF
                end if

                ! Canonical energy check
                ! test = test &
                ! & + (2E0_realk*t2%val(B,J,A,I) - t2%val(B,I,A,J))*t2%val(B,J,A,I)/deltaF

                ! t2 according to (**)
                t2%val(B,J,A,I) = t2%val(B,J,A,I)/deltaF

             end do
          end do
       end do
    end do

    call LSTIMER('SOLVE: CALC T2',tcpu,twall,DECinfo%output)




    ! 4. Transform amplitudes back to LCM/MLM basis
    ! *********************************************

    ! Since LCM/MLM and the diagonal basis are connected by a unitary
    ! transformation this basically corresponds to repeating step 2
    ! above with the transposed transformation matrices:
    !
    ! t_{aibj} = sum_{AIBJ} C_{Aa} C_{Ia} C_{Bb} C_{Jj} t_{AIBJ}

    ! 1. Transpose transformation matrices
    call array2_transpose(Cocc)
    call array2_transpose(Cvirt)

    ! 2. Index A: t(A,I,B,J) --> tmp1(a,I,B,J)
    tmp1= array4_init(dims)
    call array4_contract1(t2,Cvirt,tmp1,.true.)
    call array4_free(t2)

    ! 3. Index I: tmp1(a,I,B,J) --> tmp1(I,a,B,J) --> tmp2(i,a,B,J)
    call array4_reorder(tmp1,[2,1,3,4])
    tmp2= array4_init([nocc,nvirt,nvirt,nocc])
    call array4_contract1(tmp1,Cocc,tmp2,.true.)
    call array4_free(tmp1)

    ! 4. Index J: tmp2(i,a,B,J) --> tmp2(J,B,a,i) --> tmp1(j,B,a,i)
    call array4_reorder(tmp2,[4,3,2,1])
    tmp1 = array4_init([nocc,nvirt,nvirt,nocc])
    call array4_contract1(tmp2,Cocc,tmp1,.true.)
    call array4_free(tmp2)

    ! 5. Index B: tmp1(j,B,a,i) --> tmp1(B,j,a,i) --> t2(b,j,a,i)
    call array4_reorder(tmp1,[2,1,3,4])
    t2 = array4_init(dims)
    call array4_contract1(tmp1,Cvirt,t2,.true.)
    call array4_free(tmp1)

    call LSTIMER('SOLVE: TRANS 2',tcpu,twall,DECinfo%output)


    ! Clean up
    ! ********
    call mem_dealloc(EVocc)
    call mem_dealloc(EVvirt)
    call array2_free(Cocc)
    call array2_free(Cvirt)



  end subroutine mp2_solver_mem



  !> \brief Solve MP2 equation when RHS and t2 are values are stored on file.
  !> See mp2_solver for details about the equation.
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine mp2_solver_file(nocc,nvirt,ppfock,qqfock,RHS,t2)

    implicit none
    !> Number of occupied orbitals (dimension of ppfock)
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals (dimension of qqfock)
    integer, intent(in) :: nvirt
    !> Occupied-occupied block of Fock matrix
    real(realk) :: ppfock(nocc,nocc)
    !> Unoccupied-unoccupied block of Fock matrix
    real(realk) :: qqfock(nvirt,nvirt)
    !> RHS array, storing type 1 (see array4_init_file)
    type(array4), intent(in) :: RHS
    !> Solution array, storing type 1 (see array4_init_file)
    type(array4), intent(inout) :: t2
    type(array4) :: tmp1,tmp2
    type(array4) :: RHSaib,RHSbaj,RHStmp,t2tmp
    real(realk),pointer :: Cocc_data(:,:), Cvirt_data(:,:), Socc(:,:), Svirt(:,:)
    real(realk),pointer :: EVocc(:), EVvirt(:)
    type(array2) :: Cocc, Cvirt,CoccT,CvirtT
    integer :: I,J,A,B
    real(realk) :: tcpu, twall, deltaF
    integer :: occdims(2), virtdims(2)



    ! Strategy for solving MP2 equation:
    ! 1. Find basis where Fock matrix is diagonal
    ! 2. Transform 2-electron integrals to diagonal basis
    ! 3. In diagonal basis the solution is trivial and the amplitudes are found.
    ! 4. Transform amplitudes in diagonal basis back to LCM basis.

    ! Steps 2-4 necessarily overlap when we store array values on file.
    call LSTIMER('START',tcpu,twall,DECinfo%output)


    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Entering MP2 solver - store array values on file'
    write(DECinfo%output,*)


    ! Sanity checks
    ! *************
    ! Check that nvirt /= 0.
    if(nvirt<1 .or. nocc<1) then
       write(DECinfo%output,*) 'Number of occupied orbitals = ', nocc
       write(DECinfo%output,*) 'Number of unoccupied orbitals = ', nvirt
       call lsquit('Error in mp2_solver: Number of orbitals is smaller than one!', DECinfo%output)
    endif

    ! Only implemented for MP2
    if(DECinfo%ccModel /= 1) then
       call lsquit('mp2_solver called with other model than MP2',DECinfo%output)
    end if


    ! Initialize stuff
    ! ****************
    occdims = [nocc,nocc]
    virtdims = [nvirt,nvirt]



    ! 1. Solve Fock eigenvalue problem - each block separately
    ! ********************************************************

    ! OCCUPIED-OCCUPIED BLOCK
    ! '''''''''''''''''''''''

    ! Eigenvectors
    call mem_alloc(Cocc_data,nocc,nocc)
    ! Eigenvalues
    call mem_alloc(EVocc,nocc)

    ! The overlap matrix is simply the unit matrix because
    ! the LCM/MLM basis is orthogonal.
    call mem_alloc(Socc,nocc,nocc)
    Socc=0.0E0_realk
    do i=1,nocc
       Socc(i,i) = 1E0_realk
    end do

    ! Solve eigenvalue problem
    call solve_eigenvalue_problem(nocc,ppfock,Socc,EVocc,Cocc_data)

    ! For later, it is convenient to keep eigenvectors in array2 form
    ! and also to have the transposed matrices available
    Cocc = array2_init(occdims,Cocc_data)
    CoccT = array2_init(occdims)
    call array2_copy(CoccT,Cocc)
    call array2_transpose(CoccT)

    ! Done with some matrices
    call mem_dealloc(Cocc_data)
    call mem_dealloc(Socc)




    ! VIRTUAL-VIRTUAL BLOCK
    ! '''''''''''''''''''''

    ! Eigenvectors
    call mem_alloc(Cvirt_data,nvirt,nvirt)

    ! Eigenvalues
    call mem_alloc(EVvirt,nvirt)

    ! Unit overlap for virtual space
    call mem_alloc(Svirt,nvirt,nvirt)
    Svirt=0.0E0_realk
    do i=1,nvirt
       Svirt(i,i) = 1E0_realk
    end do

    ! Solve eigenvalue problem
    call solve_eigenvalue_problem(nvirt,qqfock,Svirt,EVvirt,Cvirt_data)

    ! For later, it is convenient to keep eigenvectors in array2 form
    ! and also to have the transposed matrices available
    Cvirt = array2_init(virtdims,Cvirt_data)
    CvirtT = array2_init(virtdims)
    call array2_copy(CvirtT,Cvirt)
    call array2_transpose(CvirtT)


    ! Done with some matrices
    call mem_dealloc(Cvirt_data)
    call mem_dealloc(Svirt)
    call LSTIMER('SOLVE: EIVAL',tcpu,twall,DECinfo%output)



    ! Transform three indices to diagonal basis (aib-->AIB)
    ! *****************************************************

    ! Using the notation that (a,i,b,j) are LCM indices
    ! and (A,I,B,J) are indices in the diagonal basis,
    ! we want to carry out the transformations:
    ! RHS_{AIBJ} = sum_{aibj} C_{aA} C_{iI} C_{bB} C_{jJ} RHS_{aibj}

    ! Since we cannot have full four-dimensional arrays in memory,
    ! we do this in steps. First, for a fixed "j" we transform the other indices:
    !
    ! RHS_{AIBj} = sum_{aib} C_{aA} C_{iI} C_{bB} RHS_{aibj}

    ! Temporary RHS where three indices are transformed, stored on file.
    RHStmp = array4_init_file([nvirt,nvirt,nocc,nocc],2,.false.)
    ! Temporary array for keeping RHS_{aibj} for a fixed j (last dimension is one).
    RHSaib = array4_init([nvirt,nocc,nvirt,1])
    call array4_open_file(RHS)
    call array4_open_file(RHStmp)


    do j=1,nocc

       ! 1. Read in RHS_{aibj} for fixed j
       call array4_read_file_type1(RHS,j,&
            & RHSaib%val(1:nvirt,1:nocc,1:nvirt,1),nvirt,nocc,nvirt)

       ! 2. Index A: RHS(a,i,b,j) --> tmp1(A,i,b,j)
       tmp1= array4_init([nvirt,nocc,nvirt,1])
       call array4_contract1(RHSaib,Cvirt,tmp1,.true.)

       ! 3. Index I: tmp1(A,i,b,j) --> tmp1(i,A,b,j) --> tmp2(I,A,b,j)
       call array4_reorder(tmp1,[2,1,3,4])
       tmp2= array4_init([nocc,nvirt,nvirt,1])
       call array4_contract1(tmp1,Cocc,tmp2,.true.)
       call array4_free(tmp1)

       ! 4. Index B: tmp2(I,A,b,j) --> tmp2(b,A,I,j) --> tmp1(B,A,I,j)
       call array4_reorder(tmp2,[3,2,1,4])
       tmp1= array4_init([nvirt,nvirt,nocc,1])
       call array4_contract1(tmp2,Cvirt,tmp1,.true.)
       call array4_free(tmp2)

       ! 5. Write to file referenced by temporary RHS array (storing type 2)
       do I=1,nocc
          call array4_write_file_type2(RHStmp,I,j,tmp1%val(1:nvirt,1:nvirt,I,1),nvirt,nvirt)
       end do
       call array4_free(tmp1)

    end do
    call array4_close_file(RHS,'KEEP')
    call array4_free(RHSaib)
    call LSTIMER('SOLVE: STEP 1',tcpu,twall,DECinfo%output)
    ! Now the file referenced by RHStmp contains RHS_{AIBj},
    ! stored in the order (B,A,I,j) using storing type 2.



    ! Transform j-->J; Solve equation in diag basis; Back transform ABJ-->abj
    ! ***********************************************************************

    ! Temporary array
    t2tmp = array4_init_file([nvirt,nvirt,nocc,nocc],2,.false.)
    call array4_open_file(t2tmp)


    I_loop: do I=1,nocc


       ! Read in RHS_{AIBj} for fixed I
       ! ------------------------------
       RHSbaj=array4_init([nvirt,nvirt,nocc,1])
       do j=1,nocc
          call array4_read_file_type2(RHStmp,I,j,&
               & RHSbaj%val(1:nvirt,1:nvirt,j,1), nvirt,nvirt)
       end do
       ! Now RHSbaj contains elements RHS_{AIBj} for a fixed I
       ! stored in the order (B,A,j,I)


       ! Transform final index j-->J
       ! ---------------------------
       ! Reorder: RHSbaj(B,A,j,I) --> RHSbaj(j,B,A,I)
       call array4_reorder(RHSbaj,[3,1,2,4])


       ! Transform: RHSbaj(j,B,A,I) --> tmp1(J,B,A,I)
       tmp1=array4_init([nocc,nvirt,nvirt,1])
       call array4_contract1(RHSbaj,Cocc,tmp1,.true.)
       call array4_free(RHSbaj)


       ! Now tmp1 contains the RHS_{BJAI} in the diagonal basis
       ! stored in the order (J,B,A,I) [fixed I].


       ! Divide by deltaF (solve equation in diagonal basis)
       ! ---------------------------------------------------
       ! In the diagonal basis the t2 solution vector is simply:
       ! t2_{BJAI} = RHS_{BJAI} / [ F_{II} + F_{JJ} - F_{AA} - F_{BB} ] (**)
       ! [See (**) in subroutine mp2_solver_mem]
       do A=1,nvirt
          do B=1,nvirt
             do J=1,nocc

                ! deltaF = F_{II} + F_{JJ} - F_{AA} - F_{BB}
                deltaF = EVocc(I) + EVocc(J) - EVvirt(A) - EVvirt(B)

                ! Sanity check
                if( abs(deltaF) < 1e-9 ) then
                   write(DECinfo%output,*) 'WARNING: SMALL NUMBERS OCCURING IN SOLVER!!!'
                   write(DECinfo%output,*) 'WARNING: SOLVER MAY BE UNSTABLE!!!'
                   write(DECinfo%output,*) 'Delta epsilon value = ', deltaF
                end if

                ! Amplitude element according to (**)
                tmp1%val(J,B,A,1) = tmp1%val(J,B,A,1)/deltaF

             end do
          end do
       end do

       ! Now tmp1 contains the t2_{BJAI} solution vector in the diagonal basis
       ! stored in the order (J,B,A,I) [fixed I], and we "just"
       ! need to transform back to the original basis.


       ! Transform back three indices: JBA --> jba
       ! -----------------------------------------
       ! Note: Use Transposed transformation matrices to back transform

       ! 1. Index j: tmp1(J,B,A,I) -->  tmp2(j,B,A,I)
       tmp2= array4_init([nocc,nvirt,nvirt,1])
       call array4_contract1(tmp1,CoccT,tmp2,.true.)
       call array4_free(tmp1)

       ! 2. Index b: tmp2(j,B,A,I) --> tmp2(B,A,j,I) --> tmp1(b,A,j,I)
       call array4_reorder(tmp2,[2,3,1,4])
       tmp1= array4_init([nvirt,nvirt,nocc,1])
       call array4_contract1(tmp2,CvirtT,tmp1,.true.)
       call array4_free(tmp2)

       ! 3. Index a: tmp1(b,A,j,I) --> tmp1(A,b,j,I) --> tmp2(a,b,j,I)
       call array4_reorder(tmp1,[2,1,3,4])
       tmp2= array4_init([nvirt,nvirt,nocc,1])
       call array4_contract1(tmp1,CvirtT,tmp2,.true.)
       call array4_free(tmp1)


       ! Write solution vector t2_{aIbj} to file using storing type 2, order: (a,b,j,I)
       ! ------------------------------------------------------------------------------

       ! Note: Now we only need to transform the last "I" index of t2 back.
       ! To do this we need to save on file such that we later can read in the full set of "I's"
       do j=1,nocc
          call array4_write_file_type2(t2tmp,j,I,&
               & tmp2%val(1:nvirt,1:nvirt,j,1), nvirt,nvirt )
       end do
       call array4_free(tmp2)

    end do I_loop

    call array4_close_file(RHStmp,'DELETE')
    call array4_free(RHStmp)
    call LSTIMER('SOLVE: STEP 2',tcpu,twall,DECinfo%output)
    ! Now the file assosicated with t2tmp contains the final amplitudes,
    ! except that the I index must be transformed back.



    ! Transform final t2 index (I-->i) and write solution vector to file
    ! ******************************************************************

    ! Final t2 solution vector
    t2 = array4_init_file([nvirt,nocc,nvirt,nocc],1,.false.)
    call array4_open_file(t2)

    do j=1,nocc

       tmp1 = array4_init([nvirt,nvirt,nocc,1])
       do I=1,nocc
          call array4_read_file_type2(t2tmp,j,I,tmp1%val(1:nvirt,1:nvirt,I,1),nvirt,nvirt)
       end do
       ! tmp1 now contains amplitudes t2_{aIbj} in the order (a,b,I,j) for fixed j.


       ! Transform final index I-->i
       ! ---------------------------

       ! tmp1(a,b,I,j) --> tmp1(I,a,b,j) --> tmp2(i,a,b,j)
       call array4_reorder(tmp1,[3,1,2,4])
       tmp2 = array4_init([nocc,nvirt,nvirt,1])
       call array4_contract1(tmp1,CoccT,tmp2,.true.)
       call array4_free(tmp1)
       ! Reorder to final storing order: tmp2(i,a,b,j) --> tmp2(a,i,b,j)
       call array4_reorder(tmp2,[2,1,3,4])


       ! Write t2_{aibj} to file for each j
       ! ----------------------------------
       call array4_write_file_type1(t2,j,&
            & tmp2%val(1:nvirt,1:nocc,1:nvirt,1),nvirt,nocc,nvirt)
       call array4_free(tmp2)

    end do

    call array4_close_file(t2tmp,'DELETE')
    call array4_free(t2tmp)
    call array4_close_file(t2,'KEEP')
    call LSTIMER('SOLVE: STEP 3',tcpu,twall,DECinfo%output)


    ! Clean up
    ! ********
    call mem_dealloc(EVocc)
    call mem_dealloc(EVvirt)
    call array2_free(Cocc)
    call array2_free(Cvirt)
    call array2_free(CoccT)
    call array2_free(CvirtT)



  end subroutine mp2_solver_file




  !> \brief Standard mp2 correlation energy
  !> \author Kasper Kristensen
  !> \return Full molecular MP2 energy
  !> \param t2 Double amplitudes
  !> \param Lmo Two-electron integrals L_{bjai} = 2*g_{bjai} - g_{ajbi}
  function get_mp2_energy(t2,Lmo) result(ecorr)

    implicit none
    type(array4), intent(in) :: Lmo,t2
    real(realk) :: ecorr

    ! Ecorr = sum_{aibj} t2_{bjai}*Lmo_{bjai}
    Ecorr=t2*Lmo

  end function get_mp2_energy




  !> \brief Coupled-cluster correlation energy
  !> \author Marcin Ziolkowski
  !> \return Full molecular CC correlation energy
  !> \param t2 Single amplitudes
  !> \param t2 Double amplitudes
  !> \param gmo Two-electron integrals (ia|jb)
  !> \param nocc Number of occupied orbitals
  !> \param nvirt Number of unoccupied orbitals
  function get_cc_energy_arrold(t1,t2,gmo,nocc,nvirt) result(ecorr)

    implicit none
    type(array2), intent(in) :: t1
    type(array4), intent(in) :: gmo,t2
    integer, intent(in) :: nocc,nvirt
    real(realk) :: ecorr,ecorr_s,ecorr_d
    integer :: a,i,b,j

    ecorr = 0.0E0_realk
    ecorr_s = 0.0E0_realk
    ecorr_d = 0.0E0_realk

    do j=1,nocc
       do b=1,nvirt
          do i=1,nocc
             do a=1,nvirt
                ecorr_d = ecorr_d + t2%val(a,i,b,j)* &
                     (2.0E0_realk*gmo%val(i,a,j,b)-gmo%val(i,b,j,a))
                ecorr_s = ecorr_s + ( t1%val(a,i)*t1%val(b,j) ) * &
                     (2.0E0_realk*gmo%val(i,a,j,b)-gmo%val(i,b,j,a))
             end do
          end do
       end do
    end do

    if(DECinfo%cc_driver_debug) then
       print *,' Singles energy : ',ecorr_s
       print *,' Doubles energy : ',ecorr_d
    end if

    ecorr = ecorr_s + ecorr_d

    return
  end function get_cc_energy_arrold

  function get_cc_energy_arrnew(t1,t2,gmo,nocc,nvirt) result(ecorr)

    implicit none
    type(array), intent(inout) :: t1
    type(array), intent(in) :: t2
    type(array), intent(inout) :: gmo
    integer, intent(in) :: nocc,nvirt
    real(realk) :: ecorr,ecorr_s,ecorr_d
    integer :: a,i,b,j

    ecorr = 0.0E0_realk
    ecorr_s = 0.0E0_realk
    ecorr_d = 0.0E0_realk

    if(t2%atype==DENSE.and.gmo%atype==DENSE.and.(t1%atype==DENSE.or.t1%atype==REPLICATED))then
      do j=1,nocc
         do b=1,nvirt
            do i=1,nocc
               do a=1,nvirt
                  ecorr_d = ecorr_d + t2%elm4(a,b,i,j)* &
                       (2.0E0_realk*gmo%elm4(i,a,j,b)-gmo%elm4(i,b,j,a))
                  ecorr_s = ecorr_s + ( t1%elm2(a,i)*t1%elm2(b,j) ) * &
                       (2.0E0_realk*gmo%elm4(i,a,j,b)-gmo%elm4(i,b,j,a))
               end do
            end do
         end do
      end do

      if(DECinfo%cc_driver_debug) then
         print *,' Singles energy : ',ecorr_s
         print *,' Doubles energy : ',ecorr_d
      end if

      ecorr = ecorr_s + ecorr_d
    elseif(t2%atype==TILED_DIST.and.gmo%atype==TILED_DIST)then
      t1%atype=REPLICATED
      call array_sync_replicated(t1)
      ecorr=get_cc_energy_parallel(t1,t2,gmo)
      t1%atype=DENSE
    endif


  end function get_cc_energy_arrnew


  !> \brief Get antisymmetrized double amplitudes
  !> \return Array4 structure with antisymmetrized double amplitudes
  function get_u(t2) result(u)

    implicit none
    type(array4), intent(inout) :: t2
    type(array4) :: u
    integer, dimension(4) :: dims
    integer :: a,i,b,j

    dims = t2%dims

#ifdef EXTRA_SIMPLE

    u = array4_init(dims)

    do j=1,dims(4)
       do b=1,dims(3)
          do i=1,dims(2)
             do a=1,dims(1)

                u%val(a,i,b,j) = 2.0E0_realk*t2%val(a,i,b,j) - t2%val(a,j,b,i)

             end do
          end do
       end do
    end do

#else

    u = array4_duplicate(t2)
    call array4_scale(u,2.0E0_realk)
    call array4_reorder(t2,[1,4,3,2])
    call array4_add_to(u,-1.0E0_realk,t2)
    call array4_reorder(t2,[1,4,3,2])

#endif

    return
  end function get_u



  !> \brief Print header and info about coupled-cluster job
  !> \author Marcin Ziolkowski
  !> \param ccPrintLevel Print level
  !> \param framgment_job Fragment job
  !> \param nbasis Number of basis functions
  !> \param nocc Number of occupied orbitals
  !> \param nvirt Number of unoccupied orbitals
  subroutine print_ccjob_header(ccPrintLevel,fragment_job,nbasis,nocc,nvirt)
    implicit none
    integer, intent(in) :: ccPrintLevel,nbasis,nocc,nvirt
    logical, intent(in) :: fragment_job

    if(ccPrintLevel > 0) then
       if(.not.fragment_job) then
          write(DECinfo%output,'(/,a)') '--------------------------'
          write(DECinfo%output,'(a)')   '  Coupled-cluster energy  '
          write(DECinfo%output,'(a,/)') '--------------------------'
          write(DECinfo%output,'(a,a)')      'Wave function    = ',DECinfo%cc_models(DECinfo%ccModel)
          write(DECinfo%output,'(a,i4)')     'MaxIter          = ',DECinfo%ccMaxIter
          write(DECinfo%output,'(a,i4)')     'Num. b.f.        = ',nbasis
          write(DECinfo%output,'(a,i4)')     'Num. occ. orb.   = ',nocc
          write(DECinfo%output,'(a,i4)')     'Num. unocc. orb. = ',nvirt
          write(DECinfo%output,'(a,e8.1e2)') 'Convergence      = ',DECinfo%ccConvergenceThreshold
          write(DECinfo%output,'(a,l1)')     'Debug mode       = ',DECinfo%cc_driver_debug
          write(DECinfo%output,'(a,i4)')     'Print level      = ',ccPrintLevel
          write(DECinfo%output,'(a,l1)')     'Use CROP         = ',DECinfo%use_crop
          write(DECinfo%output,'(a,i4)')     'CROP subspace    = ',DECinfo%ccMaxDIIS
          write(DECinfo%output,'(a,l1)')     'Preconditioner   = ',DECinfo%use_preconditioner
          write(DECinfo%output,'(a,l1)')     'Precond. B       = ',DECinfo%use_preconditioner_in_b
          write(DECinfo%output,'(a,l1)')     'Singles          = ',DECinfo%use_singles
       else
          write(DECinfo%output,'(/,a)') '  Coupled-cluster energy  -> Fragment job '
          write(DECinfo%output,'(a)')   '------------------------------------------'
          write(DECinfo%output,'(a,a,$)')       'Wave function    = ',DECinfo%cc_models(DECinfo%ccModel)
          write(DECinfo%output,'(4x,a,l1)')     'Debug mode       = ',DECinfo%cc_driver_debug
          write(DECinfo%output,'(a,i4,$)')      'MaxIter          = ',DECinfo%ccMaxIter
          write(DECinfo%output,'(5x,a,e8.1e2)') 'Convergence      = ',DECinfo%ccConvergenceThreshold
          write(DECinfo%output,'(a,i4,$)')      'Num. b.f.        = ',nbasis
          write(DECinfo%output,'(5x,a,i4)')     'Print level      = ',ccPrintLevel
          write(DECinfo%output,'(a,i4,$)')      'Num. occ. orb.   = ',nocc
          write(DECinfo%output,'(5x,a,i4)')     'CROP subspace    = ',DECinfo%ccMaxDIIS
          write(DECinfo%output,'(a,i4,$)')      'Num. unocc. orb. = ',nvirt
          write(DECinfo%output,'(5x,a,l1)')     'Preconditioner   = ',DECinfo%use_preconditioner
       end if

       ! cc parameters
       if(ccPrintLevel > 0) then
          if(fragment_job) then
             write(DECinfo%output,'(/,a,a)') &
                  '----  -------------   -------------   -------------   -------------   -------------     ------'
             write(DECinfo%output,'(a,a)') &
                  'Iter   1norm(S)        1norm(D)        2norm(S+D)      Targ-N(S+D)     Targ-energy       time  '
             write(DECinfo%output,'(a,a)') &
                  '----  -------------   -------------   -------------   -------------   -------------     ------'
          else
             write(DECinfo%output,'(/,a,a)') &
                  '----  -------------   -------------   -------------   -------------   -------------     ------'
             write(DECinfo%output,'(a,a)') &
                  'Iter   1norm(S)        1norm(D)        1norm(S+D)      2norm(S+D)      energy            time  '
             write(DECinfo%output,'(a,a)') &
                  '----  -------------   -------------   -------------   -------------   -------------     ------'
          end if
       end if
    end if

    return
  end subroutine print_ccjob_header


  !> \brief adaption of the ccsolver routine, rebuild for the 
  ! use of parallel distributed memory
  !> \author Patrick Ettenhuber (adapted from Marcin)
  subroutine ccsolver_par(ypo_f,ypv_f,fock_f,nb,no,nv, &
       & mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f,ccenergy, &
       & t1_final,t2_final,VOVO,longrange_singles)

    implicit none

    !> Number of occupied orbitals in full molecule/fragment AOS
    integer, intent(in) :: no
    !> Number of virtual orbitals in full molecule/fragment AOS
    integer, intent(in) :: nv
    !> Number of basis functions in full molecule/atomic extent
    integer, intent(in) :: nb
    !> Fock matrix in AO basis for fragment or full molecule
    real(realk), dimension(nb,nb), intent(in) :: fock_f
    !> Occupied MO coefficients for fragment/full molecule
    real(realk), dimension(nb,no), intent(in) :: ypo_f
    !> Virtual MO coefficients for fragment/full molecule
    real(realk), dimension(nb,nv), intent(in) :: ypv_f
    !> Occ-occ block of Fock matrix in MO basis
    real(realk), dimension(no,no), intent(in) :: ppfock_f
    !> Virt-virt block of Fock matrix in MO basis
    real(realk), dimension(nv,nv), intent(in) :: qqfock_f
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
    type(array) :: t1_final_work
    !> Final doubles amplitudes
    type(array4),intent(inout) :: t2_final
    type(array) :: t2_final_work
    !> Two electron integrals (a i | b j) stored as (a,i,b,j)
    type(array4),intent(inout) :: VOVO
    !> Include long-range singles effects using singles amplitudes
    !> from previous fragment calculations.
    !> IMPORTANT: If this it TRUE, then the singles amplitudes for the fragment
    !> (from previous calculations) must be stored in t1_final at input!
    logical,intent(in) :: longrange_singles
    !
    !work stuff
    real(realk),pointer :: ypo_d(:,:),ypv_d(:,:),yho_d(:,:), yhv_d(:,:),focc(:),fvirt(:)
    real(realk),pointer :: ppfock_d(:,:),qqfock_d(:,:),Uocc(:,:),Uvirt(:,:)
    integer, dimension(2) :: occ_dims, virt_dims, ao2_dims, ampl2_dims
    integer, dimension(4) :: ampl4_dims
    type(array) :: fock,ypo,ypv,yho,yhv
    type(array) :: ppfock,qqfock,pqfock,qpfock
    type(array) :: ifock,delta_fock
    type(array4) :: gao,gmo
    type(array) :: aibj,iajb
    type(array), pointer :: t2(:),omega2(:)
    type(array), pointer :: t1(:),omega1(:)
    type(array) :: omega1_opt, t1_opt, omega1_prec
    type(array) :: omega2_opt, t2_opt, omega2_prec, u
    type(array) :: xo,yo,xv,yv,h1
    type(array) :: Lmo
    !type(array2) :: xocc,yocc,xvirt,yvirt,h1
    real(realk) :: two_norm_total, one_norm_total, one_norm1, one_norm2, &
         prev_norm
    real(realk), pointer :: B(:,:),c(:)
    integer :: iter,last_iter,i,j,k,l
    logical :: crop_ok,break_iterations,saferun
    type(ri) :: l_ao
    type(array) :: ppfock_prec, qqfock_prec
    real(realk) :: tcpu, twall, ttotend_cpu, ttotend_wall, ttotstart_cpu, ttotstart_wall
    real(realk) :: iter_cpu,iter_wall
    integer :: nnodes
    real(realk), external :: ddot
    character(3) :: safefilet11,safefilet12,safefilet21,safefilet22
    !SOME DUMMIES FOR TESTING
    type(array) :: tmp
    character(ARR_MSG_LEN) :: msg
    integer :: ii,jj,aa,bb
    logical :: restart

    restart = .false.
    
    safefilet11='t11'
    safefilet12='t12'
    safefilet21='t21'
    safefilet22='t22'

    nnodes=1
#ifdef VAR_LSMPI
    nnodes=infpar%lg_nodtot
#endif

    call LSTIMER('START',ttotstart_cpu,ttotstart_wall,DECinfo%output)
    if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)


    ! Sanity check 1: Number of orbitals
    if( (nv < 1) .or. (no < 1) ) then
       write(DECinfo%output,*) 'Number of occupied orbitals = ', no
       write(DECinfo%output,*) 'Number of virtual  orbitals = ', nv
       call lsquit('ccsolver: Empty occupied or virtual space!',DECinfo%output)
    endif

    ! Sanity check 2: Singles amplitudes initiated appropriately
    if(longrange_singles) then
       if(.not. associated(t1_final%val)) then
          call lsquit('ccsolver: Long range singles corrections requested, &
               & but t1_final does not contain existing amplitudes!',DECinfo%output)
       end if
    end if


    ! go to a (pseudo) canonical basis
    call mem_alloc(focc,no)
    call mem_alloc(fvirt,nv)
    call mem_alloc(ypo_d,nb,no)
    call mem_alloc(ypv_d,nb,nv)
    call mem_alloc(yho_d,nb,no)
    call mem_alloc(yhv_d,nb,nv)
    call mem_alloc(ppfock_d,no,no)
    call mem_alloc(qqfock_d,nv,nv)
    call mem_alloc(Uocc,no,no)
    call mem_alloc(Uvirt,nv,nv)
    if(DECinfo%CCSDpreventcanonical)then
      !no diagonalization
      ypo_d   = ypo_f
      ypv_d   = ypv_f
      ppfock_d = ppfock_f
      qqfock_d = qqfock_f
      Uocc=0.0E0_realk
      Uvirt=0.0E0_realk
      do ii=1,no
        Uocc(ii,ii) = 1.0E0_realk
      enddo
      do aa=1,nv
        Uvirt(aa,aa) = 1.0E0_realk
      enddo
    else
      call get_canonical_integral_transformation_matrices(no,nv,nb,ppfock_f,qqfock_f,ypo_f,ypv_f,&
                                       & ypo_d,ypv_d,Uocc,Uvirt,focc,fvirt)
      ppfock_d = 0.0E0_realk
      qqfock_d = 0.0E0_realk
      do ii=1,no
        ppfock_d(ii,ii) = focc(ii)
      enddo
      do aa=1,nv
        qqfock_d(aa,aa) = fvirt(aa)
      enddo
    endif
    call mem_dealloc(focc)
    call mem_dealloc(fvirt)

    ! Copy MO coeffcients. It is very convenient to store them twice to handle transformation
    ! (including transposed MO matrices) efficiently. 
    yho_d = ypo_d
    yhv_d = ypv_d


    ! title
    Call print_ccjob_header(ccPrintLevel,fragment_job,nb,no,nv)
    ! dimension vectors
    occ_dims = [nb,no]
    virt_dims = [nb,nv]
    ao2_dims = [nb,nb]
    ampl4_dims = [nv,nv,no,no]
    ampl2_dims = [nv,no]

    ! create transformation matrices in array form
    ypo  = array_init(occ_dims,2)
    ypv  = array_init(virt_dims,2)
    yho  = array_init(occ_dims,2)
    yhv  = array_init(virt_dims,2)
    fock = array_init(ao2_dims,2)

    call array_convert(ypo_d,ypo)
    call array_convert(ypv_d,ypv)
    call array_convert(yho_d,yho)
    call array_convert(yhv_d,yhv)
    call array_convert(fock_f,fock)

    call mem_dealloc(ypo_d)
    call mem_dealloc(ypv_d)
    call mem_dealloc(yho_d)
    call mem_dealloc(yhv_d)
    ! Get Fock matrix correction (for fragment and/or frozen core)
    ! ************************************************************
    ! Full molecule/frozen core: The correction corresponds to difference between actual Fock matrix
    !                            and Fock matrix where the density is made from only valence orbitals.
    ! Fragment: The correction correspond to the difference between actual Fock matrix
    !           and Fock matrix calculated from a "fragment density" determined from
    !           fragment's occupied molecular orbitals (which for frozen core includes only valence
    !           orbitals).

    ! Density corresponding to input MOs
    call mem_alloc(dens,nb,nb)
    call get_density_from_occ_orbitals(nb,no,ypo%elm2,dens)

    if(fragment_job) then ! fragment: calculate correction

       ifock=array_init(ao2_dims,2)

       if(longrange_singles) then
          ! Get Fock matrix using singles amplitudes from previous
          ! fragment calculation, thereby effectively including long-range
          ! correlated polarization effects
          call Get_AOt1Fock(mylsitem,t1_final,ifock,no,nv,nb,ypo,yho,yhv)
       else
          ! Fock matrix for fragment from density made from input MOs
           call get_fock_matrix_for_dec(nb,dens,mylsitem,ifock,.true.)
       end if

       ! Long range Fock correction:
       !delta_fock = getFockCorrection(fock,ifock)
       delta_fock=array_init(ao2_dims,2)

       call array_cp_data(fock,delta_fock)
       call array_add(delta_fock,-1.0E0_realk,ifock)
       call array_free(ifock)

    else 
       ! Full molecule: deltaF = F(Dcore) for frozen core (0 otherwise)
       if(DECinfo%frozencore) then
          ! Fock matrix from input MOs
          ifock=array_init(ao2_dims,2)
          call get_fock_matrix_for_dec(nb,dens,mylsitem,ifock,.true.)
          ! Correction to actual Fock matrix
          delta_fock=array_init(ao2_dims,2)
          call array_cp_data(fock,delta_fock)
          call array_add(delta_fock,-1.0E0_realk,ifock)
          call array_free(ifock)
       else
          delta_fock=array_init(ao2_dims,2)
          call array_zero(delta_fock)
       end if
    end if

    call mem_dealloc(dens)

    ! get two-electron integrals in ao
    if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a)') 'debug :: calculating AO integrals'

    ! Only calculate full 4-dimensional AO integrals for old/debug mode

    ! Simulate two-electron integrals (debug mode)
!    if(DECinfo%simulate_eri .or. DECinfo%fock_with_ri) then
!       if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a)') &
!            'debug :: calculate RI intermediate - temporary'
!       l_ao = get_ao_ri_intermediate(mylsitem)
!       if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a)') &
!            'debug :: intermediates done'
!    end if

    if(DECinfo%PL>1) call LSTIMER('CCSOL: INIT',tcpu,twall,DECinfo%output)
    if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! get fock matrices for preconditioning
    Preconditioner : if(DECinfo%use_preconditioner .or. DECinfo%use_preconditioner_in_b) then
       ppfock_prec = array_minit_rpseudo_dense([no,no],2)
       qqfock_prec = array_minit_rpseudo_dense([nv,nv],2)

       call array_change_atype_to_rep(ppfock_prec)
       call array_change_atype_to_rep(qqfock_prec)
       if(DECinfo%precondition_with_full) then
          call array_convert(ppfock_d,ppfock_prec)
          call array_convert(qqfock_d,qqfock_prec)
       else
          tmp = array_init([nb,no],2)
          call array_contract_outer_indices_rl(1.0E0_realk,fock,yho,0.0E0_realk,tmp)
          call array_contract_outer_indices_ll(1.0E0_realk,ypo,tmp,0.0E0_realk,ppfock_prec)
          call array_free(tmp)

          tmp = array_init([nb,nv],2)
          call array_contract_outer_indices_rl(1.0E0_realk,fock,yhv,0.0E0_realk,tmp)
          call array_contract_outer_indices_ll(1.0E0_realk,ypv,tmp,0.0E0_realk,qqfock_prec)
          call array_free(tmp)
       end if
       call array_change_atype_to_d(ppfock_prec)
       call array_change_atype_to_d(qqfock_prec)
    end if Preconditioner

    call mem_dealloc(ppfock_d)
    call mem_dealloc(qqfock_d)

    ! allocate things
    if(DECinfo%use_singles) then
       call mem_alloc(t1,DECinfo%ccMaxIter)
       call mem_alloc(omega1,DECinfo%ccMaxIter)
       ppfock=array_init([no,no],2)
       pqfock=array_init([no,nv],2)
       qpfock=array_init([nv,no],2)
       qqfock=array_init([nv,nv],2)
    end if
    call mem_alloc(t2,DECinfo%ccMaxIter)
    call mem_alloc(omega2,DECinfo%ccMaxIter)


    ! initialize T1 matrices and fock transformed matrices for CC pp,pq,qp,qq
    if(DECinfo%ccModel /= 1) then
       xo = array_init(occ_dims,2)
       yo = array_init(occ_dims,2)
       xv = array_init(virt_dims,2)
       yv = array_init(virt_dims,2)
    end if
    !iajb=array_minit_tdpseudo_dense([no,nv,no,nv],4)
    iajb=array_minit_td([no,nv,no,nv],4)
    call array_zero(iajb)

    call mem_alloc(B,DECinfo%ccMaxIter,DECinfo%ccMaxIter)
    call mem_alloc(c,DECinfo%ccMaxIter)



    ! readme : the iteration sequence is universal and may be used for all
    !          iterative cc models (linear or non-linear) and is
    !          semi-independent on the storage of vectors (allocation and
    !          deallocation, etc)

    ! iterate
    break_iterations = .false.
    crop_ok = .false.
    prev_norm = 1.0E6_realk


    print *
    print *, '### Starting CC iterations'
    print *, '### ----------------------'
    print '(1X,a)',  '###  Iteration     Residual norm          CC energy'

    write(DECinfo%output,*)
    write(DECinfo%output,*) '### Starting CC iterations'
    write(DECinfo%output,*) '### ----------------------'
    write(DECinfo%output,'(1X,a)')  '###  Iteration     Residual norm          CC energy'


    CCIteration : do iter=1,DECinfo%ccMaxIter

       if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)
       call LSTIMER('START',iter_cpu,iter_wall,DECinfo%output)

       ! remove old vectors
       RemoveOldVectors : if(iter > DECinfo%ccMaxDIIS) then
          !print *,"remove old"
          if(DECinfo%cc_driver_debug) then
             write(DECinfo%output,'(a,i4)') ' debug :: vector to delete : ',iter-DECinfo%ccMaxDIIS
          end if

          if(DECinfo%use_singles) then
             call array_free_rpseudo_dense(t1(iter-DECinfo%ccMaxDIIS))
             Call array_free(omega1(iter-DECinfo%ccMaxDIIS))
             
          end if
          call array_free(t2(iter-DECinfo%ccMaxDIIS))
          call array_free(omega2(iter-DECinfo%ccMaxDIIS))
       end if RemoveOldVectors

       ! get guess amplitude vectors in the first iteration --> zero
       GetGuessVectors : if(iter == 1) then
          if(DECinfo%use_singles)then
            t1(iter) = array_minit_rpseudo_dense(ampl2_dims,2)
            t2(iter) = array_minit_tdpseudo_dense(ampl4_dims,4)
            call get_guess_vectors(restart,t2(iter),safefilet21,safefilet22,t1(iter),safefilet11,safefilet12)
          else
            t2(iter) = array_minit_tdpseudo_dense(ampl4_dims,4)
            call get_guess_vectors(restart,t2(iter),safefilet21,safefilet22)
         endif
       end if GetGuessVectors

       ! Initialize residual vectors
       if(DECinfo%use_singles)then
         omega1(iter) = array_init(ampl2_dims,2)
         call array_zero(omega1(iter))
       endif
       omega2(iter) = array_minit_td(ampl4_dims,4)
       call array_zero(omega2(iter))

       ! get singles
       T1Related : if(DECinfo%use_singles) then

          ! get the T1 transformation matrices
          call array_cp_data(yhv,yv)
          call array_cp_data(ypv,xv)
          call array_contract_outer_indices_rr(-1.0E0_realk,ypo,t1(iter),1.0E0_realk,xv)

          call array_cp_data(yho,yo)
          call array_cp_data(ypo,xo)
          call array_contract_outer_indices_rl(1.0E0_realk,yhv,t1(iter),1.0E0_realk,yo)

          ! get inactive fock
          !if(DECinfo%fock_with_ri) then
          !   ! Debug mode
          !   ifock = getInactiveFockFromRI(l_ao,xocc,yocc,h1)
          !end if


       end if T1Related
       if(DECinfo%PL>1) call LSTIMER('CCIT: INIT',tcpu,twall,DECinfo%output)
       if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)


       ! readme : get residuals, so far this solver supports only singles and doubles
       !          amplitudes (extension to higher iterative model is trivial); differences
       !          are mainly based on the set of residual vectors to be evaluated

       ! MODIFY FOR NEW MODEL
       ! If you implement a new model, please insert call to your own residual routine here!
       SelectCoupledClusterModel : if(DECinfo%ccModel==1) then
          !call getDoublesResidualMP2_simple(Omega2(iter),t2(iter),gmo,ppfock,qqfock, &
          !     & no,nv)
          call lsquit("ERROR(ccsolver_par):no mp2 implemented --> use&
          & ccsolver",DECinfo%output)
       elseif(DECinfo%ccModel==2) then
          call lsquit("ERROR(ccsolver_par):no cc2 implemented --> use&
          & ccsolver",DECinfo%output)
       else  ! CCSD or CCSD(T)

          call get_ccsd_residual_integral_driven(delta_fock%elm1,omega2(iter),t2(iter),&
             & fock%elm1,iajb,no,nv,ppfock%elm1,qqfock%elm1,pqfock%elm1,qpfock%elm1,xo%elm1,&
             & xv%elm1,yo%elm1,yv%elm1,nb,MyLsItem,omega1(iter)%elm1,iter,rest=restart)

       end if SelectCoupledClusterModel
       
       if(DECinfo%PL>1) call LSTIMER('CCIT: RESIDUAL',tcpu,twall,DECinfo%output)


       if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)

       ! calculate crop/diis matrix
       B=0.0E0_realk
       do i=iter,max(iter-DECinfo%ccMaxDIIS+1,1),-1
          do j=iter,i,-1
             if(DECinfo%use_singles) then
                if(DECinfo%use_preconditioner_in_b) then

                   omega1_prec = precondition_singles(omega1(j),ppfock_prec,qqfock_prec)
                   omega2_prec = precondition_doubles(omega2(j),ppfock_prec,qqfock_prec)
                   B(i,j) = array_ddot(omega1(i),omega1_prec) 
                   B(i,j) = B(i,j) + array_ddot(omega2(i),omega2_prec)
                   call array_free(omega1_prec)
                   call array_free(omega2_prec)
                else
                   !B(i,j) = d_omega1(i)*d_omega1(j) + d_omega2(i)*d_omega2(j)
                end if
            else
                ! just doubles
                if(DECinfo%use_preconditioner_in_b) then
                   omega2_prec = precondition_doubles(omega2(j),ppfock_prec,qqfock_prec)
                   B(i,j) = array_ddot(omega2(i),omega2_prec)
                   call array_free(omega2_prec)
                else
                  print *,"STOP, NOT YET IMPLEMENTED(ccsolver_par)"
                  stop 0
                  ! B(i,j) = ddot(omega2(i)%nelms,omega2(i)%elm1,1,omega2(j)%elm1,1)
                  ! if(DECinfo%cc_driver_debug) write(DECinfo%output,'(a,i4,a,i4,a,f16.10)') &
                  !      ' debug :: B(',i,',',j,')=',B(i,j)
                end if
             end if
             B(j,i) = B(i,j)
          end do
       end do
       !msg="DIIS mat, new"
       !call print_norm(B,DECinfo%ccMaxIter*DECinfo%ccMaxIter,msg)

       if(DECinfo%PL>1) call LSTIMER('CCIT: CROP MAT',tcpu,twall,DECinfo%output)
       ! solve crop/diis equation
       c=0.0E0_realk
       call CalculateDIIScoefficients(DECinfo%ccMaxDIIS,DECinfo%ccMaxIter,iter,B,c, &
            DECinfo%cc_driver_debug)


       if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)


       
       ! mixing omega to get optimal
       if(DECinfo%use_singles) then
          t1_opt     = array_init(ampl2_dims,2)
          omega1_opt = array_init(ampl2_dims,2)
       end if
       omega2_opt  = array_minit_td(ampl4_dims,4)
       call array_zero(omega2_opt)
       t2_opt = array_minit_td(ampl4_dims,4)
       call array_zero(t2_opt)
       do i=iter,max(iter-DECinfo%ccMaxDIIS+1,1),-1
          ! mix singles
          if(DECinfo%use_singles) then
            call array_add(omega1_opt,c(i),omega1(i))
            call array_add(t1_opt,c(i),t1(i))
          end if
          ! mix doubles
          call array_add(omega2_opt,c(i),omega2(i))
          call array_add(t2_opt,c(i),t2(i))
       end do


       if(DECinfo%PL>1) call LSTIMER('CCIT: MIXING',tcpu,twall,DECinfo%output)
       if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)

       ! if crop, put the optimal in place of trial (not for diis)
       if(DECinfo%use_crop) then
          if(DECinfo%use_singles) then
             call array_cp_data(omega1_opt,omega1(iter))
             call array_cp_data(t1_opt,t1(iter))
          end if
          call array_cp_data(omega2_opt,omega2(iter))
          call array_cp_data(t2_opt,t2(iter))
       end if

       if(DECinfo%PL>1) call LSTIMER('CCIT: COPY OPT',tcpu,twall,DECinfo%output)
       if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)

       ! check for the convergence
       one_norm1 = 0.0E0_realk
       one_norm2 = 0.0E0_realk
       if(DECinfo%use_singles)call print_norm(omega1(iter),one_norm1,.true.)
       call print_norm(omega2(iter),one_norm2,.true.)
       one_norm_total = one_norm1 + one_norm2
       two_norm_total = sqrt(one_norm_total)
       !if(iter==3)then
       !  print*,"SETTING TWONORM TO QUIT";two_norm_total=0.9E-5_realk
       !endif
       ! simple crop diagnostics
       if(two_norm_total < prev_norm) then
          crop_ok=.true.
       else
          crop_ok=.false.
          write(DECinfo%output,'(a)') ' warning :: total norm was smaller in previous iteration !!! '
       end if
       prev_norm=two_norm_total
       ! check if this is the last iteration
       if(iter == DECinfo%ccMaxIter .or. two_norm_total < DECinfo%ccConvergenceThreshold) &
            break_iterations=.true.


       if(DECinfo%PL>1) call LSTIMER('CCIT: CONV',tcpu,twall,DECinfo%output)
       if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)


       ! calculate the correlation energy and fragment energy
       ! MODIFY FOR NEW MODEL
       ! If you implement a new model, please insert call to energy routine here,
       ! or insert a call to get_cc_energy if your model uses the standard CC energy expression.
       ! Note: This routine uses massive parallelization, if you do not want your model to use
       ! this you do not need to make modifications here.
       EnergyForCCmodel: if(DECinfo%ccmodel==1) then  
          ! MP2
          call lsquit("ERROR(ccsolver_par):CCD/MP2 energy not yet implemented",-1)
       elseif(DECinfo%ccmodel==2 .or. DECinfo%ccmodel==3 .or. DECinfo%ccmodel==4 ) then
          ! CC2, CCSD, or CCSD(T) (for (T) calculate CCSD contribution here)
          ccenergy = get_cc_energy(t1(iter),t2(iter),iajb,no,nv)
       end if EnergyForCCmodel

       if(DECinfo%PL>1) call LSTIMER('CCIT: ENERGY',tcpu,twall,DECinfo%output)
       if(DECinfo%PL>1) call LSTIMER('START',tcpu,twall,DECinfo%output)

       ! generate next trial vector if this is not the last iteration
       if(.not.break_iterations) then
          if(DECinfo%use_preconditioner) then
             if(DECinfo%use_singles) then
                omega1_prec = precondition_singles(omega1_opt,ppfock_prec,qqfock_prec)
                t1(iter+1) = array_minit_rpseudo_dense(ampl2_dims,2)
                call array_cp_data(t1_opt,t1(iter+1))
                call array_add(t1(iter+1),1.0E0_realk,omega1_prec)
                call array_free(omega1_prec)
             end if
             omega2_prec = precondition_doubles(omega2_opt,ppfock_prec,qqfock_prec)
             t2(iter+1) = array_minit_tdpseudo_dense(ampl4_dims,4)
             call array_cp_data(t2_opt,t2(iter+1))
             call array_add(t2(iter+1),1.0E0_realk,omega2_prec)
             call array_free(omega2_prec)
          else
             if(DECinfo%use_singles)then
                t1(iter+1) = array_minit_rpseudo_dense(ampl2_dims,2)
                call array_cp_data(t1_opt,t1(iter+1))
                call array_add(t1(iter+1),1.0E0_realk,omega1_opt)
             endif
             t2(iter+1) = array_minit_tdpseudo_dense(ampl4_dims,4)
             call array_cp_data(t2_opt,t2(iter+1))
             call array_add(t2(iter+1),1.0E0_realk,omega2_opt)
          end if

          !if DECinfo%CCSDsaferun option is set, make sure data is in dense
          if(DECinfo%CCSDsaferun)then
            if(DECinfo%use_singles)then
              call save_current_guess(iter,t2(iter+1),safefilet21,safefilet22,&
              &t1(iter+1),safefilet11,safefilet12)
            else
              call save_current_guess(iter,t2(iter+1),safefilet21,safefilet22)
            endif
          endif
       end if

       if(DECinfo%PL>1) call LSTIMER('CCIT: NEXT VEC',tcpu,twall,DECinfo%output)


       ! delete optimals
       if(DECinfo%use_singles) then
          call array_free(t1_opt)
          call array_free(omega1_opt)
       end if
       call array_free(t2_opt)
       call array_free(omega2_opt)

       call LSTIMER('CC ITERATION',iter_cpu,iter_wall,DECinfo%output)

#ifdef __GNUC__
       call flush(DECinfo%output)
#endif

       print '(1X,a,2X,i4,5X,g19.9,4X,g19.9)',  '### ',iter, two_norm_total,ccenergy
       write(DECinfo%output,'(1X,a,2X,i4,5X,g19.9,4X,g19.9)') &
            &   '### ',iter, two_norm_total,ccenergy
       last_iter = iter
       if(break_iterations) exit
         
   end do CCIteration

   call LSTIMER('START',ttotend_cpu,ttotend_wall,DECinfo%output)


    write(DECinfo%output,*)
    write(DECinfo%output,'(/,a)') '-------------------------------'
    write(DECinfo%output,'(a)')   '  Coupled-cluster job summary  '
    write(DECinfo%output,'(a,/)') '-------------------------------'
    if(break_iterations) then
       write(DECinfo%output,'(a)')     'Hooray! CC equation is solved!'
    else
       write(DECinfo%output,'(a,i4,a)')  'CC equation not solved in ', &
            & DECinfo%ccMaxIter, ' iterations!'
       call lsquit('CC equation not solved!',DECinfo%output)
    end if
    write(DECinfo%output,'(a,f16.3,a)') 'CCSOL: Total cpu time    = ',ttotend_cpu-ttotstart_cpu,' s'
    write(DECinfo%output,'(a,f16.3,a)') 'CCSOL: Total wall time   = ',ttotend_wall-ttotstart_wall,' s'

    if(fragment_job) then
       write(DECinfo%output,'(a,f16.10)')  'Frag. corr. energy = ',ccenergy
    else
       write(DECinfo%output,'(a,f16.10)')  'Corr. energy       = ',ccenergy
    end if
    write(DECinfo%output,'(a,i5)') 'Number of CC iterations =', last_iter


    ! Free memory and save final amplitudes
    ! *************************************

    ! remove rest of the singles amplitudes and residuals
    do i=last_iter,max(last_iter-DECinfo%ccMaxDIIS+1,1),-1
       if(DECinfo%use_singles) then
          ! Save final singles amplitudes
          if(i==last_iter) then
             if(.not.longrange_singles) then ! just copy
                t1_final = array2_init(ampl2_dims)
             end if
             call dcopy(int(t1(i)%nelms),t1(last_iter)%elm1,1,t1_final%val,1)
          end if

          ! Free singles amplitudes and residuals
          call array_free_rpseudo_dense(t1(i))
          call array_free(omega1(i))

       end if


       ! Save final double amplitudes
       if(i==last_iter) then
          t2_final = array4_init([nv,no,nv,no])
          call array_cp_tiled2dense(t2(last_iter),.true.)
          call array_reorder_4d(1.0E0_realk,t2(last_iter)%elm1,nv,nv,no,no,[1,3,2,4],0.0E0_realk,t2_final%val)
          call array_change_atype_to_td(t2(last_iter))
       end if

       ! Free doubles residuals
       call array_free_tdpseudo_dense(omega2(i))
       ! Free doubles amplitudes
       call array_free_tdpseudo_dense(t2(i))

    end do


    ! Save two-electron integrals in the order (virt,occ,virt,occ)
    if(DECinfo%ccModel == 1) then
            print *,"not implemented"
            stop 0
       !call array4_free(lmo) ! also free lmo integrals
       !VOVO = array4_duplicate(gmo)
       !call array4_free(gmo)
    else
       VOVO = array4_init([no,nv,no,nv])
       call array_convert(iajb,VOVO%val)
       call array_free(iajb)
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
    call array_free(delta_fock)


    if(DECinfo%simulate_eri .or. DECinfo%fock_with_ri) then
       call ri_free(l_ao)
       call ri_reset()
    end if


    if(DECinfo%use_preconditioner .or. DECinfo%use_preconditioner_in_b) then
       call array_free_rpseudo_dense(ppfock_prec)
       call array_free_rpseudo_dense(qqfock_prec)
    end if

    if(DECinfo%use_singles) then
       !call array2_free(h1)
       call array_free(xo)
       call array_free(yo)
       call array_free(xv)
       call array_free(yv)
       call array_free(pqfock)
       call array_free(qpfock)
    end if

    call array_free(ppfock)
    call array_free(qqfock)

    call array_free(ypo)
    call array_free(yho)
    call array_free(ypv)
    call array_free(yhv)
    call array_free(fock)


    !transform back to original basis   
    if(DECinfo%use_singles)then
      call ccsolver_can_local_trans(VOVO%val,t2_final%val,no,nv,Uocc,Uvirt,t1_final%val)
    else
      call ccsolver_can_local_trans(VOVO%val,t2_final%val,no,nv,Uocc,Uvirt)
    endif

    call mem_dealloc(Uocc)
    call mem_dealloc(Uvirt)

#ifdef MOD_UNRELEASED
    call array4_print_statistics(DECinfo%output)
    call array_print_mem_info(DECinfo%output,.true.,.false.)
#endif

  end subroutine ccsolver_par

  !> \brief: transform ccsd_doubles, ccsdpt_singles and ccsdpt_doubles from canonical to local basis
  !> \author: Patrick Ettenhuber adapted from Janus Juul Eriksen
  !> \date: April 2013
  !> \param: t2, gvovo, t1, no and nv are nocc and nvirt, respectively, 
  !<         and U_occ and U_virt are unitary matrices from canonical --> local basis
  subroutine ccsolver_can_local_trans(gvovo,t2,no,nv,Uocc,Uvirt,t1)

    implicit none
    !> integers
    integer, intent(in) :: no, nv
    !> ccsd_doubles and ccsdpt_doubles
    real(realk), intent(inout) :: t2(nv*nv*no*no), gvovo(nv*nv*no*no)
    !> unitary transformation matrices
    real(realk), intent(inout) :: Uocc(no*no), Uvirt(nv*nv)
    !> ccsdpt_singles
    real(realk), intent(inout),optional :: t1(nv*no)
    !> temp array2 and array4 structures
    real(realk),pointer :: tmp(:)

    call mem_alloc(tmp,nv*nv*no*no)

    ! (a,i,b,j) are local basis indices and (A,I,B,J) refer to the canonical basis.
    ! on input t2 and gvovo are ordered AIBJ and t1 AI

    call successive_xyxy_trafo(nv,no,t2,Uvirt,Uocc,tmp)

    !successive transformation of gvovo:
    call successive_xyxy_trafo(nv,no,gvovo,Uvirt,Uocc,tmp)

    !if t1 trafo has to be done as well
    if(present(t1))then
      !U(a,A) t(AI)    -> t(aI)
      call dgemm('n','n',nv,no,nv,1.0E0_realk,Uvirt,nv,t1,nv,0.0E0_realk,tmp,nv)
      ! tmp(aI) U(i,I)^T   -> t(ai)
      call dgemm('n','t',nv,no,no,1.0E0_realk,tmp,nv,Uocc,no,0.0E0_realk,t1,nv)
    endif

    call mem_dealloc(tmp)
  end subroutine ccsolver_can_local_trans

  subroutine successive_xyxy_trafo(x,y,XYXY,XX,YY,WRKYXYX)
    implicit none
    integer, intent(in) :: x,y
    real(realk), intent(inout) :: XYXY(x*y*x*y),WRKYXYX(y*x*y*x)
    real(realk), intent(in) :: XX(x,x),YY(y,y)
    !XYXY(X,YXY)^T XX(x,X)^T   -> WRKYXYX (YXY,x)
    call dgemm('t','t',x*y*y,x,x,1.0E0_realk,XYXY,x,XX,x,0.0E0_realk,WRKYXYX,x*y*y)
    ! WRKYXYX(Y,XYx)^T YY(y,Y)^T   -> XYXY (XYx,y)
    call dgemm('t','t',x*x*y,y,y,1.0E0_realk,WRKYXYX,y,YY,y,0.0E0_realk,XYXY,x*x*y)
    ! XYXY(X,Yxy)^T XX(x,X)^T   -> WRKYXYX (Yxy,x)
    call dgemm('t','t',x*y*y,x,x,1.0E0_realk,XYXY,x,XX,x,0.0E0_realk,WRKYXYX,x*y*y)
    ! WRKYXYX(Y,xyx)^T YY(y,Y)^T   -> XYXY (xyxy)
    call dgemm('t','t',x*x*y,y,y,1.0E0_realk,WRKYXYX,y,YY,y,0.0E0_realk,XYXY,x*x*y)
  end subroutine successive_xyxy_trafo

  !> \brief should be a general subroutine to get the guess amplitudes when
  !starting up a CCSD or CC2 calculation and checks for files which contain
  !amplitudes of a previous calculation. If none are found the usual zero guess
  !is returned
  !> \author Patrick Ettenhuber
  !> \date December 2012
  subroutine get_guess_vectors(restart,t2,safefilet21,safefilet22,t1,safefilet11,safefilet12)
    implicit none
    logical,intent(out) :: restart
    !> contains the guess doubles amplitudes on output
    type(array),intent(inout) :: t2
    !> the filenames to check for valid doubles amplitudes
    character(3),intent(in) :: safefilet21,safefilet22
    !> contains the singles amplitudes on output
    type(array),intent(inout),optional :: t1
    !> the filenames to check for valid singles amplitudes
    character(3),intent(in), optional :: safefilet11,safefilet12
    integer :: no,nv
    integer :: fu_t11,fu_t12,fu_t21,fu_t22,fu_t1,fu_t2
    logical :: file_exists11,file_exists12,file_exists21,file_exists22
    logical :: file_status11,file_status12,file_status21,file_status22
    logical :: readfile1, readfile2
    integer :: saved_iter11,saved_iter12,saved_iter21,saved_iter22,iter_start
    integer :: saved_nel11,saved_nel12,saved_nel21,saved_nel22
    logical :: all_singles
    character(11) :: fullname11, fullname12, fullname21, fullname22
    character(ARR_MSG_LEN) :: msg
    all_singles=present(t1).and.present(safefilet11).and.present(safefilet12)
    !print *,"CHECK INPUT",safefilet11,safefilet12,all_singles,DECinfo%use_singles
    fu_t11=111
    fu_t12=112
    fu_t21=121
    fu_t22=122
    nv=t2%dims(1)
    no=t2%dims(3)


    iter_start=0
    !check for safe files of the amplitudes in the current directory and read
    !them if they exist and ok
    readfile1=.false.
    readfile2=.false.
    
    !this can be skipped by input, but restart is default
    if(.not.DECinfo%CCSDno_restart)then
      if(DECinfo%use_singles.and.all_singles)then
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
            call lsquit("ERROR(ccsolver_par):wrong dimensions in amplitude &
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
            call lsquit("ERROR(ccsolver_par):wrong dimensions in amplitude &
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
          call lsquit("ERROR(ccsolver_par):wrong dimensions in amplitude &
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
          call lsquit("ERROR(ccsolver_par):wrong dimensions in amplitude &
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
      READ(fu_t1)t1%elm1
      CLOSE(fu_t1)
      restart = .true.
    else
      call array_zero(t1)
    endif
    if(readfile2)then
      READ(fu_t2) t2%elm1
      CLOSE(fu_t2)
      restart = .true.
    else
      call array_zero(t2)
    endif
  end subroutine get_guess_vectors

  !> \brief Subroutine to save the current guess amplitudes for the next
  !iteration
  !> \author Patrick Ettenhuber
  !> \date Dezember 2012
  subroutine save_current_guess(iter,t2,safefilet21,safefilet22,&
  &t1,safefilet11,safefilet12)
    implicit none
    !> iteration number
    integer,intent(in) :: iter
    !> doubles guess amplitudes for the next iteration
    type(array), intent(in) :: t2
    !> alternating filenames for the doubles amplitudes
    character(3),intent(in) :: safefilet21,safefilet22
    !> singles guess amplitudes for the next iteration
    type(array), intent(in), optional :: t1
    !> alternating filenames for the singles amplitudes
    character(3),intent(in), optional :: safefilet11,safefilet12
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
    all_singles=present(t1).and.present(safefilet11).and.present(safefilet12)
    fu_t11=111
    fu_t12=112
    fu_t21=121
    fu_t22=122
    if(DECinfo%use_singles.and.all_singles)then
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
        WRITE(fu_t11)t1%nelms
        WRITE(fu_t11)t1%elm1
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
        WRITE(fu_t12)t1%nelms
        WRITE(fu_t12)t1%elm1
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
      WRITE(fu_t21)t2%nelms
      WRITE(fu_t21)t2%elm1
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
      WRITE(fu_t22)t2%nelms
      WRITE(fu_t22)t2%elm1
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
  end subroutine save_current_guess


end module ccdriver

! LocalWords:  iter
