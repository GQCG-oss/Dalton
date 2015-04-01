!> @file
!> DEC-CCSD(T) routines
!> \brief: ccsd(t) module
!> \author: Janus Juul Eriksen
!> \date: 2012-2014, Aarhus
module ccsdpt_module

#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif
  use precision
  use dec_typedef_module
  use memory_handling
  use lstiming!, only: lstimer
  use screen_mod!, only: DECscreenITEM
  use BUILDAOBATCH
  use typedeftype!, only: Lsitem,lssetting
  use screen_mod!,only: free_decscreen, DECSCREENITEM
  use IntegralInterfaceDEC!, only: II_precalc_DECScreenMat,&
!       & II_getBatchOrbitalScreen, II_GET_DECPACKED4CENTER_J_ERI
  use IntegralInterfaceMOD
  use IchorErimoduleHost
  use Fundamental, only: bohr_to_angstrom
  use tensor_interface_module
  use lspdm_tensor_operations_module
  use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
#ifdef VAR_OPENACC
  use openacc
#endif
#if defined(VAR_CUDA) || defined(VAR_OPENACC)
  use gpu_interfaces
#endif

  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
#ifdef VAR_MPI
  use decmpi_module
#endif
  use dec_workarounds_module
  use crop_tools_module
  use cc_tools_module
  use dec_fragment_utils
  use array2_simple_operations
  use array3_simple_operations
  use array4_simple_operations
  
#ifdef MOD_UNRELEASED
  public :: ccsdpt_driver,ccsdpt_energy_e4_frag,ccsdpt_energy_e5_frag,&
       & ccsdpt_energy_e4_pair, ccsdpt_energy_e5_pair, ccsdpt_energy_e5_ddot, &
       & ccsdpt_decnp_e4_frag, ccsdpt_decnp_e5_frag
  private
#endif

contains

#ifdef MOD_UNRELEASED

  !> \brief: driver routine for dec-ccsd(t)
  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  subroutine ccsdpt_driver(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,mylsitem,vovo,ccsd_doubles,&
                         & ccsdpt_singles,print_frags,abc,ccsdpt_doubles,e4)

    implicit none

    !> nocc, nvirt, and nbasis for fragment or full molecule
    integer, intent(in) :: nocc, nvirt, nbasis
    !> ppfock and qqfock for fragment or full molecule
    real(realk), intent(in) :: ppfock(:,:), qqfock(:,:)
    !> mo coefficents for occ and virt space for fragment or full molecule
    real(realk), intent(in) :: Co(:,:), Cv(:,:)
    !> mylsitem for fragment or full molecule
    type(lsitem), intent(inout) :: mylsitem
    !> ccsd doubles amplitudes
    type(tensor), intent(inout) :: ccsd_doubles
    !> incoming vovo integrals
    type(tensor), intent(inout) :: vovo
    !> input for the actual triples computation
    type(tensor),intent(inout) :: ccsdpt_singles
    logical :: print_frags,abc
    type(tensor),intent(inout),optional :: ccsdpt_doubles
    real(realk),optional :: e4
    !> 2-el integrals
    ! ijk scheme
    type(tensor) :: ovoo ! integrals (AI|JK) in the order (J,A,I,K)
    ! vvvo is of type DENSE, if this is a serial calculation, and TILED_DIST,
    ! if this is a parallel calculation
    type(tensor) :: vvvo ! integrals (AI|BC) in the order (C,B,A,I)
    ! abc scheme
    type(tensor) :: ooov ! integrals (AI|JK) in the order (I,J,K,A)
    ! vovv is of type DENSE, if this is a serial calculation, and TILED_DIST,
    ! if this is a parallel calculation
    type(tensor) :: vovv ! integrals (AI|BC) in the order (B,I,A,C)
    integer :: nodtotal
    integer :: ijk_nbuffs,abc_nbuffs,abc_tile_size
    !> orbital energies
    real(realk), pointer :: eivalocc(:), eivalvirt(:)
    !> MOs and unitary transformation matrices
    type(array2) :: C_can_occ, C_can_virt, Uocc, Uvirt
    !> dimensions
    integer, dimension(2) :: occdims, virtdims, virtoccdims,occAO,virtAO
    integer, dimension(3) :: dims_aaa
    integer, dimension(4) :: dims_iaai, dims_aaii
    logical :: master
    type(tensor) :: ccsdpt_doubles_2
#ifdef VAR_OPENACC
    !> device type
    integer(acc_device_kind) :: acc_device_type
#endif
    real(realk) :: tcpu,twall

    call time_start_phase(PHASE_WORK)

!#ifdef VAR_OPENACC
!
!    ! probe for device type
!    acc_device_type = acc_get_device_type()
!
!    ! initialize the device
!    call acc_init(acc_device_type)
!
!#endif

    master = .true.
    nodtotal = 1
#ifdef VAR_MPI
    master = (infpar%lg_mynum .eq. infpar%master)
    nodtotal = infpar%lg_nodtot
#endif

    if (master) then

       call LSTIMER('START',tcpu,twall,DECinfo%output)

       write(DECinfo%output,*) ''
       write(DECinfo%output,*) ''
       write(DECinfo%output,*) '=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*'
       write(DECinfo%output,*) '        Inside the CCSD(T) driver routine.        '
       write(DECinfo%output,*) '*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*='
       write(DECinfo%output,*) ''

    endif

    ! init dimensions
    occdims     = [nocc,nocc]
    virtdims    = [nvirt,nvirt]
    virtoccdims = [nvirt,nocc]
    dims_iaai   = [nocc,nvirt,nvirt,nocc]
    dims_aaii   = [nvirt,nvirt,nocc,nocc]
    dims_aaa    = [nvirt,nvirt,nvirt]
    occAO       = [nbasis,nocc]
    virtAO      = [nbasis,nvirt]

    if (print_frags) then

       !Zero to be able to sum up 
       call tensor_zero(ccsdpt_singles)
       call tensor_zero(ccsdpt_doubles)

    else

       !Zero to be able to sum up
       call tensor_zero(ccsdpt_singles)
 
       if (present(e4)) then

          e4 = 0.0E0_realk

       else

          call lsquit('print_frags == .false., but e4 is missing... aborting.',DECinfo%output) 

       endif

    endif

    call mem_alloc(eivalocc,nocc)
    call mem_alloc(eivalvirt,nvirt)
    C_can_occ  = array2_init(occAO)
    C_can_virt = array2_init(virtAO)

    if (master) then
      ! *************************************
      ! get arrays for transforming integrals
      ! *************************************
      ! C_can_occ, C_can_virt:  MO coefficients for canonical basis
      ! Uocc, Uvirt: unitary transformation matrices for canonical --> local basis (and vice versa)
      ! note: Uocc and Uvirt have indices (local,canonical)

      Uocc       = array2_init(occdims)
      Uvirt      = array2_init(virtdims)
      call get_canonical_integral_transformation_matrices(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,&
                         & C_can_occ%val,C_can_virt%val,Uocc%val,Uvirt%val,eivalocc,eivalvirt)

      ! ************************************************************
      ! transform vovo and ccsd doubles amplitudes to diagonal basis
      ! ************************************************************
      if (abc) then

         call local_can_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,oovv=vovo%elm1)
         call local_can_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,oovv=ccsd_doubles%elm1)

      else

         call local_can_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,vvoo=vovo%elm1)
         call local_can_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,vvoo=ccsd_doubles%elm1)

      endif

    end if

#ifdef VAR_MPI
    call time_start_phase(PHASE_COMM)

    ! bcast the JOB specifier and distribute data to all the slaves within local group
    waking_the_slaves: if ((nodtotal .gt. 1) .and. master) then

       ! slaves are in lsmpi_slave routine (or corresponding dec_mpi_slave) and are now awaken
       call ls_mpibcast(CCSDPTSLAVE,infpar%master,infpar%lg_comm)

       ! distribute ccsd doubles and fragment or full molecule quantities to the slaves
       call mpi_communicate_ccsdpt_calcdata(nocc,nvirt,nbasis,vovo%elm4,ccsd_doubles%elm4,&
                                          & mylsitem,print_frags,abc)

    end if waking_the_slaves
#endif

    call ccsdpt_info(nbasis,nocc,nvirt,print_frags,abc,ijk_nbuffs,abc_nbuffs,abc_tile_size,nodtotal)

#ifdef VAR_MPI
    ! Communicate important information:
    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    call ls_mpi_buffer(eivalocc,nocc,infpar%master)
    call ls_mpi_buffer(eivalvirt,nvirt,infpar%master)
    call ls_mpi_buffer(C_can_occ%val,nbasis,nocc,infpar%master)
    call ls_mpi_buffer(C_can_virt%val,nbasis,nvirt,infpar%master)
    call ls_mpi_buffer(ijk_nbuffs,infpar%master)
    call ls_mpi_buffer(abc_nbuffs,infpar%master)
    call ls_mpi_buffer(abc_tile_size,infpar%master)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

    call time_start_phase(PHASE_WORK)
#endif

    ! ********************************************
    ! get vo³ and v³o integrals in proper sequence
    ! ********************************************
    ! note: the integrals are calculated in canonical basis

    if (abc) then

       call get_CCSDpT_integrals_abc(mylsitem,nbasis,nocc,nvirt,C_can_occ%val,C_can_virt%val,ooov,vovv,abc_tile_size)

    else

       call get_CCSDpT_integrals_ijk(mylsitem,nbasis,nocc,nvirt,C_can_occ%val,C_can_virt%val,ovoo,vvvo)

    endif

    write(DECinfo%output,*) ''
    write(DECinfo%output,*) ''
    write(DECinfo%output,*) '=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*='
    write(DECinfo%output,*) '        Done with CCSD(T) integrals        '
    write(DECinfo%output,*) '*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*'
    write(DECinfo%output,*) ''

    ! release occ and virt canonical MOs
    call array2_free(C_can_occ)
    call array2_free(C_can_virt)


    ! ********************************
    ! begin actual triples calculation
    ! ********************************

    ! in all comments in the below, we employ the notation of eqs. (14.6.60) [with (i,j,k)/(a,c,d)]
    ! and (14.6.64).

    ! objective is three-fold:
    ! 1) calculate triples amplitudes, collect in array3 structures, trip_*** [canonical basis]
    ! 2) calculate ^{*}T^{a}_{i} and ^{*}T^{ab}_{ij} amplitudes in array2 and array4 structures, 
    !    ccsdpt_singles and ccsdpt_doubles [canonical basis]
    !    here: ccsdpt_doubles_2 is a temp array towards the generation of ccsdpt_doubles
    ! 3) transform ccsd_doubles, ccsdpt_singles and ccsdpt_doubles into local basis [local basis]

    ! *****************************************************
    ! ***************** trip generation *******************
    ! *****************************************************

    ! init ccsdpt_doubles_2 array structure.
    ! we merge ccsdpt_doubles and ccsdpt_doubles_2 at the end into ccsdpt_doubles. 
    ! we have dimensioned ccsdpt_doubles as dims_aaii and ccsdpt_doubles_2 as dims_iaai 
    ! in order to load in data consecutive in memory inside ccsdpt_contract_21 
    ! and ccsdpt_contract_22, respectively.
    if (print_frags) then

       if (abc) then

          call tensor_init(ccsdpt_doubles_2, [nvirt,nocc,nocc,nvirt],4)

       else

          call tensor_init(ccsdpt_doubles_2, [nocc,nvirt,nvirt,nocc],4)
 
       endif

       call tensor_zero(ccsdpt_doubles_2)

    endif

    !************************************************************!
    ! here: the main  (t) loop: this is where the magic happens! !
    !************************************************************!

#ifdef VAR_MPI

    call time_start_phase(PHASE_WORK)

    if (abc) then

       ! the parallel version of the abc-loop
       if (print_frags) then

          call abc_loop_par(nocc,nvirt,ooov%elm1,vovo%elm1,vovv,ccsd_doubles%elm1,&
                          & eivalocc,eivalvirt,nodtotal,abc_nbuffs,abc_tile_size,&
                          & ccsdpt_singles%elm1,ccsdpt_doubles%elm1,ccsdpt_doubles_2%elm1)

       else

          call abc_loop_par(nocc,nvirt,ooov%elm1,vovo%elm1,vovv,ccsd_doubles%elm1,&
                          & eivalocc,eivalvirt,nodtotal,abc_nbuffs,abc_tile_size,&
                          & ccsdpt_singles%elm1,e4=e4)

       endif

    else

       ! the parallel version of the ijk-loop
       if (print_frags) then
   
          call ijk_loop_par(nocc,nvirt,ovoo%elm1,vovo%elm1,vvvo,ccsd_doubles%elm1,&
                          & eivalocc,eivalvirt,nodtotal,ijk_nbuffs,&
                          & ccsdpt_singles%elm1,ccsdpt_doubles%elm1,ccsdpt_doubles_2%elm1)
   
       else
   
          call ijk_loop_par(nocc,nvirt,ovoo%elm1,vovo%elm1,vvvo,ccsd_doubles%elm1,&
                          & eivalocc,eivalvirt,nodtotal,ijk_nbuffs,&
                          & ccsdpt_singles%elm1,e4=e4)
   
       endif

    endif

    call time_start_phase(PHASE_WORK)

#else

    if (abc) then

       ! the serial version of the abc-loop
       if (print_frags) then

          call abc_loop_ser(nocc,nvirt,ooov%elm1,vovo%elm1,vovv,ccsd_doubles%elm1,&
                          & eivalocc,eivalvirt,ccsdpt_singles%elm1,&
                          & ccsdpt_doubles%elm1,ccsdpt_doubles_2%elm1)

       else

          call abc_loop_ser(nocc,nvirt,ooov%elm1,vovo%elm1,vovv,ccsd_doubles%elm1,&
                          & eivalocc,eivalvirt,ccsdpt_singles%elm1,e4=e4)

       endif

    else

       ! the serial version of the ijk-loop
       if (print_frags) then
   
          call ijk_loop_ser(nocc,nvirt,ovoo%elm1,vovo%elm1,vvvo%elm1,ccsd_doubles%elm1,&
                          & eivalocc,eivalvirt,ccsdpt_singles%elm1,&
                          & ccsdpt_doubles%elm1,ccsdpt_doubles_2%elm1)
   
       else
   
          call ijk_loop_ser(nocc,nvirt,ovoo%elm1,vovo%elm1,vvvo%elm1,ccsd_doubles%elm1,&
                          & eivalocc,eivalvirt,ccsdpt_singles%elm1,e4=e4)
   
       endif

    endif

#endif

    ! *******************************************
    ! *********** done w/ main loop *************
    ! *******************************************

#ifdef VAR_MPI

    ! here, synchronize all procs
    call time_start_phase(PHASE_IDLE)
    call lsmpi_barrier(infpar%lg_comm)

    ! reduce singles and doubles arrays into that residing on the master
    reducing_to_master: if (nodtotal .gt. 1) then

       call time_start_phase(PHASE_COMM)

       call lsmpi_local_reduction(ccsdpt_singles%elm1,ccsdpt_singles%nelms,infpar%master)

       if (print_frags) then

          call lsmpi_local_reduction(ccsdpt_doubles%elm1,ccsdpt_doubles%nelms,infpar%master)
          call lsmpi_local_reduction(ccsdpt_doubles_2%elm1,ccsdpt_doubles_2%nelms,infpar%master)

       else

          call lsmpi_local_reduction(e4,infpar%master)

       endif

       call time_start_phase(PHASE_WORK)

    end if reducing_to_master

    ! release stuff located on slaves
    releasing_the_slaves: if ((nodtotal .gt. 1) .and. .not. master) then

       call time_start_phase(PHASE_WORK)

       ! release stuff initialized herein
       if (print_frags) call tensor_free(ccsdpt_doubles_2) 
       call mem_dealloc(eivalocc)
       call mem_dealloc(eivalvirt)

       ! release o^3v and v^3o integrals
       if (abc) then

          call tensor_free(ooov)
          call tensor_free(vovv)

       else

          call tensor_free(ovoo)
          call tensor_free(vvvo)

       endif

       ! now, release the slaves  
       return

    end if releasing_the_slaves

    call time_start_phase(PHASE_WORK)

#endif

    ! now everything resides on the master...

    if (print_frags) then

       ! collect ccsdpt_doubles and ccsdpt_doubles_2 into ccsdpt_doubles array structure
       ! ccsdpt_doubles(a,b,i,j) = ccsdpt_doubles(a,b,i,j) + ccsdpt_doubles_2(j,a,b,i) (*)
       ! (*) here, ccsdpt_doubles_2 is simultaneously reordered as (j,a,b,i) --> (a,b,i,j)
       call array_reorder_4d(1.0E0_realk,ccsdpt_doubles_2%elm1,ccsdpt_doubles_2%dims(1),&
                                  &ccsdpt_doubles_2%dims(2),ccsdpt_doubles_2%dims(3),ccsdpt_doubles_2%dims(4),&
                                  &[2,3,4,1],1.0E0_realk,ccsdpt_doubles%elm1)
   
       ! release ccsdpt_doubles_2 array structure
       call tensor_free(ccsdpt_doubles_2)

    endif

    ! release o^3v and v^3o integrals
    if (abc) then

       call tensor_free(ooov)
       call tensor_free(vovv)

    else

       call tensor_free(ovoo)
       call tensor_free(vvvo)

    endif

    ! *************************************************
    ! ***** do canonical --> local transformation *****
    ! *************************************************

    if (print_frags) then

       if (abc) then

          call can_local_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,oovv=ccsdpt_doubles%elm1,ov=ccsdpt_singles%elm1)
          call can_local_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,oovv=ccsd_doubles%elm1)
          call can_local_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,oovv=vovo%elm1)

       else

          call can_local_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,vvoo=ccsdpt_doubles%elm1,vo=ccsdpt_singles%elm1)
          call can_local_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,vvoo=ccsd_doubles%elm1)
          call can_local_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,vvoo=vovo%elm1)

       endif

    else

       if (abc) then

          call can_local_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,ov=ccsdpt_singles%elm1)

       else

          call can_local_trans(nocc,nvirt,nbasis,Uocc%val,Uvirt%val,vo=ccsdpt_singles%elm1)

       endif

    endif

    ! now, release Uocc and Uvirt
    call array2_free(Uocc)
    call array2_free(Uvirt)

    ! clean up
    call mem_dealloc(eivalocc)
    call mem_dealloc(eivalvirt)

!#ifdef VAR_OPENACC
!
!    ! shut down the device
!    call acc_shutdown(acc_device_type)
!
!#endif

    if (master) call LSTIMER('CCSDPT_DRIVER (TOTAL)',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine ccsdpt_driver


#ifdef VAR_MPI
  !> \brief: main ijk-loop (parallel version)
  !> \author: Janus Juul Eriksen
  !> \date: january 2014
  subroutine ijk_loop_par(nocc,nvirt,ovoo,vvoo,vvvo,ccsd_doubles,&
                        & eivalocc,eivalvirt,nodtotal,nbuffs,ccsdpt_singles,&
                        & ccsdpt_doubles,ccsdpt_doubles_2,e4)

    implicit none

    !> nocc,nvirt
    integer, intent(in)      :: nocc,nvirt
    !> 2-el integrals
    real(realk), dimension(nocc,nvirt,nocc,nocc) :: ovoo ! integrals (AI|JK) in the order (J,A,I,K)
    real(realk), dimension(nvirt,nvirt,nocc,nocc) :: vvoo ! integrals (AI|BJ) in the order (A,B,I,J)
    type(tensor), intent(inout)  :: vvvo ! integrals (AI|BC) in the order (C,B,A,I)
    !> ccsd doubles amplitudes
    real(realk), dimension(nvirt,nvirt,nocc,nocc) :: ccsd_doubles
    ! o*v^2 portions of ccsd_doubles
    real(realk), pointer, dimension(:,:,:) :: ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k
    !> triples amplitudes and 3d work array
    real(realk), pointer, dimension(:,:,:) :: trip_tmp, trip_ampl
    !> ccsd(t) intermediates 
    real(realk), dimension(nvirt,nocc) :: ccsdpt_singles
    real(realk), dimension(nvirt,nvirt,nocc,nocc),optional :: ccsdpt_doubles
    real(realk), dimension(nocc,nvirt,nvirt,nocc),optional :: ccsdpt_doubles_2
    real(realk),optional :: e4
    logical :: full_no_frags
    !> orbital energies
    real(realk), intent(inout)  :: eivalocc(nocc), eivalvirt(nvirt)
    integer, intent(in) :: nbuffs
    !> job distribution
    real(realk), pointer, dimension(:) :: vvvo_pdm_i,vvvo_pdm_j,vvvo_pdm_k ! v^3 tiles from cbai
    real(realk), pointer, dimension(:,:) :: vvvo_pdm_buff      ! buffers to prefetch tiles
    integer :: b_size,njobs,nodtotal,ij,ij_count,i_old,j_old
    integer, pointer :: ij_array(:),jobs(:)
    !> loop integers
    integer :: i,j,k,idx,ij_type,tuple_type
    !> ij loop and k loop buffer handling
    integer :: ibuf, jbuf, kbuf
    integer,pointer :: tiles_in_buf(:)
    integer(kind=ls_mpik), pointer :: req(:)
    logical,pointer :: needed(:)
    integer :: num_ids,m
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind), pointer, dimension(:) :: async_id
    integer(kind=acc_device_kind) :: acc_device_type
#else
    integer, pointer, dimension(:) :: async_id
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat
    ! timings
    real(realk) :: tcpu,twall,time_pt_ijk,time_pt_ijk_min,time_pt_ijk_max
    real(realk) :: time_trip,time_efull,time_driv,time_preload
    real(realk) :: time_trip_min,time_efull_min,time_driv_min,time_preload_min
    real(realk) :: time_trip_max,time_efull_max,time_driv_max,time_preload_max
    real(realk) :: time_trip_tot,time_efull_tot,time_driv_tot,time_preload_tot
    real(realk) :: phase_cntrs(nphases)
    real(realk) :: flushing_time, flushing_time_min, flushing_time_max
    real(realk) :: unlock_time,   unlock_time_min,   unlock_time_max
    real(realk) :: waiting_time,  waiting_time_min,  waiting_time_max
    real(realk) :: time_w_min, time_w_max
    real(realk) :: time_c_min, time_c_max
    real(realk) :: time_i_min, time_i_max

    ! init timings
    unlock_time   = time_lsmpi_win_unlock
    waiting_time  = time_lsmpi_wait
    flushing_time = time_lsmpi_win_flush
    time_trip_tot = 0.0E0_realk; time_preload_tot = 0.0E0_realk; time_efull_tot = 0.0E0_realk; time_driv_tot = 0.0E0_realk
    call time_start_phase( PHASE_WORK, twall = time_pt_ijk )
    call time_phases_get_current(current_wt=phase_cntrs)
    if (infpar%lg_mynum .eq. infpar%master) call LSTIMER('START',tcpu,twall,DECinfo%output)

    full_no_frags = .false.

    if (present(e4)) full_no_frags = .true.

    call time_start_phase(PHASE_WORK)

    call mem_alloc( vvvo_pdm_buff, nvirt**3, nbuffs )
    call mem_alloc( needed,       nbuffs )
    call mem_alloc( tiles_in_buf, nbuffs )
    call mem_alloc( req,          nbuffs )
    if(alloc_in_dummy)call tensor_lock_wins(vvvo,'s',all_nodes=.true.)

    ! init ccsd_doubles_help_arrays
    call mem_alloc(ccsd_doubles_portions_i,nocc,nvirt,nvirt)
    call mem_alloc(ccsd_doubles_portions_j,nocc,nvirt,nvirt)
    call mem_alloc(ccsd_doubles_portions_k,nocc,nvirt,nvirt)

    ! init triples tuples structure
    call mem_alloc(trip_ampl,nvirt,nvirt,nvirt)
    ! init 3d wrk array
    call mem_alloc(trip_tmp,nvirt,nvirt,nvirt)

    ! create job distribution list
    ! first, determine common batch size from number of tasks and nodes
    ! in the ij matrix, njobs is the number of elements in the lower triangular matrix
    ! always an even number [ n(n+1) is always an even number ]
    njobs = int((nocc**2 + nocc)/2)
    b_size = int(njobs/nodtotal)

    ! ij_array stores all jobs for composite ij indices in descending order
    call mem_alloc(ij_array,njobs)
    ! init list (one more than b_size since mod(njobs,nodtotal) is not necessearily zero
    call mem_alloc(jobs,b_size + 1)

    ! create ij_array
    call create_ij_tensor_ccsdpt(njobs,nocc,ij_array)
    ! fill the list
    call job_distrib_ccsdpt(b_size,njobs,ij_array,jobs)

    ! release ij_array
    call mem_dealloc(ij_array)

    ! now follows the main loop

    ! a note on the mpi scheme.
    ! since we (in a dec picture) often have many nodes compared to nocc, we explicitly collapse the i- and j-loop.
    ! by doing this, we are guaranteed that all nodes participate.
    ! the composite index ij is incremented in the collapsed loop, and we may calculate i and j from ij.

    ! init ij and i_old/j_old
    ij           = 0
    i_old        = 0
    j_old        = 0
    ij_type      = 0
    needed       = .false.
    tiles_in_buf = -1

    ! set async handles. if we are not using gpus, just set them to arbitrary negative numbers
    ! handle 1: ccsd_doubles
    ! handle 2: vvvo and ovoo integrals
    ! handle 3: vvoo integrals and ccsdpt_doubles / ccsdpt_doubles_2 intermediates
    ! handle 4: triples amplitudes
    ! handle 5: energy evaluation
    num_ids = 5
    call mem_alloc(async_id,num_ids)

#ifdef VAR_OPENACC

    if (DECinfo%acc_sync) then
       async_id = acc_async_sync
    else
       do m = 1,num_ids
          async_id(m) = int(m,kind=acc_handle_kind)
       enddo
    endif

#else

    if (DECinfo%acc_sync) then
       async_id = 0
    else
       do m = 1,num_ids
          async_id(m) = -m
       enddo
    endif

#endif

#ifdef VAR_CUBLAS

    ! initialize the CUBLAS context
    stat = cublasCreate_v2(cublas_handle)

#endif

!$acc wait

!$acc enter data create(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k)&
!$acc& copyin(eivalvirt,ccsdpt_singles,e4) if(full_no_frags)
!
!$acc enter data create(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k)&
!$acc& copyin(eivalvirt,ccsdpt_singles) if(.not. full_no_frags)

!$acc wait

 ijrun_par: do ij_count = 1,b_size + 1

               ! get value of ij from job disttribution list
               ij = jobs(ij_count)

               ! no more jobs to be done? otherwise leave the loop
               if (ij .lt. 0) exit

               ! calculate i and j from composite ij value
               call calc_i_leq_j(ij,nocc,i,j)

               ! has the i and j index changed?
               if ((i .eq. i_old) .and. (j .ne. j_old)) then

                  ij_type = 1

               else if ((i .ne. i_old) .and. (j .eq. j_old)) then

                  ij_type = 2

               else ! (i .ne. i_old) .and. (j .ne. j_old))

                  ij_type = 3

               end if

               !FIND i in buffer
               call assoc_ptr_to_buf(i,vvvo,nbuffs,tiles_in_buf,needed,vvvo_pdm_i,vvvo_pdm_buff,ibuf,req)

               !FIND j in buffer
               call assoc_ptr_to_buf(j,vvvo,nbuffs,tiles_in_buf,needed,vvvo_pdm_j,vvvo_pdm_buff,jbuf,req)

               ! select ij combination
               TypeOf_ij_combi: select case(ij_type)

               case(1) ! i .gt. j

!$acc enter data copyin(ccsd_doubles(:,:,:,i),ccsd_doubles(:,:,:,j)) async(async_id(1))

#ifdef VAR_OPENACC
                  call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,j),nvirt,nvirt,&
                          & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_j,async_id(1))
#else
                  call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,j),nvirt,nvirt,&
                          & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_j)
#endif

                  ! get the j'th v^3 tile only
                  call time_start_phase(PHASE_COMM)

                  if( alloc_in_dummy )then

                     call lsmpi_wait(req(jbuf))

                  else

                     if(vvvo%lock_set(j)) call tensor_unlock_win(vvvo,j)

                  endif

                  call time_start_phase(PHASE_WORK)

!$acc enter data copyin(vvvo_pdm_i,vvvo_pdm_j,&
!$acc& ovoo(:,:,i,j),ovoo(:,:,j,i)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,i,j),vvoo(:,:,j,i)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,i,j),vvoo(:,:,j,i),&
!$acc& ccsdpt_doubles_2(:,:,:,i),ccsdpt_doubles_2(:,:,:,j),&
!$acc& ccsdpt_doubles(:,:,i,j),ccsdpt_doubles(:,:,j,i)) async(async_id(3)) if(.not. full_no_frags)

                  ! store j index
                  j_old = j

               case(2)

                  if (i .eq. j) then

! if we can use ccsd_doubles(:,:,:,j), ccsd_doubles_portions_j, vvvo_pdm_j, and ccsdpt_doubles_2(:,:,:,j),
! then there is no need for this (although we need ccsd_doubles(:,:,:,j)

!$acc enter data copyin(ccsd_doubles(:,:,:,i)) async(async_id(1))

#ifdef VAR_OPENACC
                     call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,i),nvirt,nvirt,&
                             & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i,async_id(1))
#else
                     call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,i),nvirt,nvirt,&
                             & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i)
#endif

                     ! get the i'th v^3 tile
                     call time_start_phase(PHASE_COMM)

                     if( alloc_in_dummy )then

                        call lsmpi_wait(req(ibuf))

                     else

                        if(vvvo%lock_set(i)) call tensor_unlock_win(vvvo,i)

                     endif


                     call time_start_phase(PHASE_WORK)

!$acc enter data copyin(vvvo_pdm_i,&
!$acc& ovoo(:,:,i,j)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,i,j)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,i,j),&
!$acc& ccsdpt_doubles_2(:,:,:,i),&
!$acc& ccsdpt_doubles(:,:,i,j)) async(async_id(3)) if(.not. full_no_frags)

                  else ! i .gt. j

!$acc enter data copyin(ccsd_doubles(:,:,:,i),ccsd_doubles(:,:,:,j)) async(async_id(1))

#ifdef VAR_OPENACC
                     call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,i),nvirt,nvirt,&
                             & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i,async_id(1))
#else
                     call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,i),nvirt,nvirt,&
                             & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i)
#endif

                     ! get the i'th v^3 tile only
                     call time_start_phase(PHASE_COMM)

                     if( alloc_in_dummy )then

                        call lsmpi_wait(req(ibuf))

                     else

                        if(vvvo%lock_set(i)) call tensor_unlock_win(vvvo,i)

                     endif

                     call time_start_phase(PHASE_WORK)

!$acc enter data copyin(vvvo_pdm_i,vvvo_pdm_j,&
!$acc& ovoo(:,:,i,j),ovoo(:,:,j,i)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,i,j),vvoo(:,:,j,i)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,i,j),vvoo(:,:,j,i),&
!$acc& ccsdpt_doubles_2(:,:,:,i),ccsdpt_doubles_2(:,:,:,j),&
!$acc& ccsdpt_doubles(:,:,i,j),ccsdpt_doubles(:,:,j,i)) async(async_id(3)) if(.not. full_no_frags)

                  end if

                  ! store i index
                  i_old = i

               case(3)

                  if (i .eq. j) then

! if we can use ccsd_doubles(:,:,:,j), ccsd_doubles_portions_j, vvvo_pdm_j, and ccsdpt_doubles_2(:,:,:,j),
! then these should be used here instead of the ones for 'i'

!$acc enter data copyin(ccsd_doubles(:,:,:,i)) async(async_id(1))

#ifdef VAR_OPENACC
                     call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,i),nvirt,nvirt,&
                             & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i,async_id(1))
#else
                     call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,i),nvirt,nvirt,&
                             & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i)
#endif

                     ! get the i'th v^3 tile
                     call time_start_phase(PHASE_COMM)

                     if( alloc_in_dummy )then

                        call lsmpi_wait(req(ibuf))

                     else

                        if(vvvo%lock_set(i)) call tensor_unlock_win(vvvo,i)

                     endif


                     call time_start_phase(PHASE_WORK)

!$acc enter data copyin(vvvo_pdm_i,&
!$acc& ovoo(:,:,i,j)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,i,j)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,i,j),&
!$acc& ccsdpt_doubles_2(:,:,:,i),&
!$acc& ccsdpt_doubles(:,:,i,j)) async(async_id(3)) if(.not. full_no_frags)

                  else ! i .gt. j

!$acc enter data copyin(ccsd_doubles(:,:,:,i),ccsd_doubles(:,:,:,j)) async(async_id(1))

#ifdef VAR_OPENACC
                     call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,i),nvirt,nvirt,&
                             & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i,async_id(1))
                     call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,j),nvirt,nvirt,&
                             & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_j,async_id(1))
#else
                     call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,i),nvirt,nvirt,&
                             & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i)
                     call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,j),nvirt,nvirt,&
                             & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_j)
#endif

                     ! get the i'th and j'th v^3 tiles
                     call time_start_phase(PHASE_COMM)

                     if( alloc_in_dummy )then

                        call lsmpi_wait(req(ibuf))

                     else

                        if(vvvo%lock_set(i)) call tensor_unlock_win(vvvo,i)

                     endif


                     if( alloc_in_dummy )then

                        call lsmpi_wait(req(jbuf))

                     else

                        if(vvvo%lock_set(j)) call tensor_unlock_win(vvvo,j)

                     endif

                     call time_start_phase(PHASE_WORK)

!$acc enter data copyin(vvvo_pdm_i,vvvo_pdm_j,&
!$acc& ovoo(:,:,i,j),ovoo(:,:,j,i)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,i,j),vvoo(:,:,j,i)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,i,j),vvoo(:,:,j,i),&
!$acc& ccsdpt_doubles_2(:,:,:,i),ccsdpt_doubles_2(:,:,:,j),&
!$acc& ccsdpt_doubles(:,:,i,j),ccsdpt_doubles(:,:,j,i)) async(async_id(3)) if(.not. full_no_frags)

                  end if

                  ! store i and j indices
                  i_old = i
                  j_old = j

               end select TypeOf_ij_combi

               needed(ibuf) = .true.
               needed(jbuf) = .true.


        krun_par: do k=1,j

                     ! select type of tuple
                     tuple_type = -1

                     !FIND k in buffer
                     call assoc_ptr_to_buf(k,vvvo,nbuffs,tiles_in_buf,needed,vvvo_pdm_k,vvvo_pdm_buff,kbuf,req)

                     if ((i .eq. j) .and. (j .eq. k)) then

                        ! i == j == k
                        ! this always gives zero contribution
                        cycle

                     else if ((i .eq. j) .and. (j .gt. k)) then

                        ! i == j > k
                        tuple_type = 1

!$acc enter data copyin(ccsd_doubles(:,:,:,k)) async(async_id(1))

#ifdef VAR_OPENACC
                        call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,k),nvirt,nvirt,&
                                & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_k,async_id(1))
#else
                        call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,k),nvirt,nvirt,&
                                & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_k)
#endif

                        ! get the k'th tile
                        call time_start_phase(PHASE_COMM)

                        if( alloc_in_dummy )then

                           call lsmpi_wait(req(kbuf))

                        else

                           if(vvvo%lock_set(k)) call tensor_unlock_win(vvvo,k)

                        endif

                        call time_start_phase(PHASE_WORK)

!$acc enter data copyin(vvvo_pdm_k,&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,i,k),vvoo(:,:,k,i)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,i,k),vvoo(:,:,k,i),&
!$acc& ccsdpt_doubles_2(:,:,:,k),&
!$acc& ccsdpt_doubles(:,:,i,k),ccsdpt_doubles(:,:,k,i)) async(async_id(3)) if(.not. full_no_frags)

                     else if ((i .gt. j) .and. (j .eq. k)) then

                        ! i > j == k
                        tuple_type = 2

!$acc enter data copyin(ovoo(:,:,j,k)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,j,k)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,j,k),&
!$acc& ccsdpt_doubles(:,:,j,k)) async(async_id(3)) if(.not. full_no_frags)

                     else

                        ! i > j > k
                        tuple_type = 3

!$acc enter data copyin(ccsd_doubles(:,:,:,k)) async(async_id(1))

#ifdef VAR_OPENACC
                        call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,k),nvirt,nvirt,&
                                & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_k,async_id(1))
#else
                        call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,k),nvirt,nvirt,&
                                & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_k)
#endif

                        ! get the k'th tile
                        call time_start_phase(PHASE_COMM)

                        if( alloc_in_dummy )then

                           call lsmpi_wait(req(kbuf))

                        else

                           if(vvvo%lock_set(k)) call tensor_unlock_win(vvvo,k)

                        endif

                        call time_start_phase(PHASE_WORK)

!$acc enter data copyin(vvvo_pdm_k,&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i),ovoo(:,:,j,k),ovoo(:,:,k,j)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,i,k),vvoo(:,:,k,i),vvoo(:,:,j,k),vvoo(:,:,k,j)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,i,k),vvoo(:,:,k,i),vvoo(:,:,j,k),vvoo(:,:,k,j),&
!$acc& ccsdpt_doubles_2(:,:,:,k),&
!$acc& ccsdpt_doubles(:,:,i,k),ccsdpt_doubles(:,:,k,i),&
!$acc& ccsdpt_doubles(:,:,j,k),ccsdpt_doubles(:,:,k,j)) async(async_id(3)) if(.not. full_no_frags)

                     end if

                     needed(kbuf) = .true.

                     call time_start_phase(PHASE_WORK, twall = time_preload )
                     call preload_tiles_in_bg_buf(vvvo,jobs,b_size,nocc,i,j,k,ij_count,nbuffs,needed,tiles_in_buf,vvvo_pdm_buff,req)
                     call time_start_phase(PHASE_WORK, ttot = time_preload )
                     time_preload_tot = time_preload_tot + time_preload

                     ! generate tuple(s)
                     TypeOfTuple_par_ijk: select case(tuple_type)

                     case(1)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                        call time_start_phase(PHASE_WORK, twall = time_trip )
                        call trip_generator_ijk_case1(i,k,nocc,nvirt,ccsd_doubles(:,:,:,i),ccsd_doubles(:,:,:,k),&
                                                & ccsd_doubles_portions_i,ccsd_doubles_portions_k,&
                                                & vvvo_pdm_i,vvvo_pdm_k,&
                                                & ovoo(:,:,i,i),ovoo(:,:,i,k),ovoo(:,:,k,i),&
                                                & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                        call time_start_phase(PHASE_WORK, ttot = time_trip )
                        time_trip_tot = time_trip_tot + time_trip

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,k)) async(async_id(1))

                        if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(vvvo_pdm_k,&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                           call time_start_phase(PHASE_WORK, twall = time_efull )
                           call ccsdpt_energy_full_ijk_case1(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                & vvoo(:,:,i,i),vvoo(:,:,i,k),vvoo(:,:,k,i),&
                                                & ccsdpt_singles(:,i),ccsdpt_singles(:,k),&
                                                & e4,async_id,num_ids,cublas_handle)
                           call time_start_phase(PHASE_WORK, ttot = time_efull )
                           time_efull_tot = time_efull_tot + time_efull

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,k),vvoo(:,:,k,i)) async(async_id(3))
 
                        else

                           call time_start_phase(PHASE_WORK, twall = time_driv )
#ifdef VAR_OPENACC
                           call trip_denom_ijk_acc(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else
                           call trip_denom_ijk_cpu(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif   

!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))
 
                           call ccsdpt_driver_ijk_case1(i,k,nocc,nvirt,vvoo(:,:,i,i),vvoo(:,:,i,k),vvoo(:,:,k,i),&
                                                & ovoo(:,:,i,i),ovoo(:,:,i,k),ovoo(:,:,k,i),&
                                                & vvvo_pdm_i,vvvo_pdm_k,&
                                                & ccsdpt_singles(:,i),ccsdpt_singles(:,k),&
                                                & ccsdpt_doubles(:,:,i,i),ccsdpt_doubles(:,:,i,k),&
                                                & ccsdpt_doubles(:,:,k,i),ccsdpt_doubles_2(:,:,:,i),&
                                                & ccsdpt_doubles_2(:,:,:,k),trip_tmp,trip_ampl,&
                                                & async_id,num_ids,cublas_handle)
                           call time_start_phase(PHASE_WORK, ttot = time_driv )
                           time_driv_tot = time_driv_tot + time_driv

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(vvvo_pdm_k,&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,k),vvoo(:,:,k,i))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,k),&
!$acc& ccsdpt_doubles(:,:,i,k),ccsdpt_doubles(:,:,k,i)) async(async_id(3))

                        endif

                     case(2)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                        call time_start_phase(PHASE_WORK, twall = time_trip )
                        call trip_generator_ijk_case2(i,j,nocc,nvirt,ccsd_doubles(:,:,:,i),ccsd_doubles(:,:,:,j),&
                                                & ccsd_doubles_portions_i,ccsd_doubles_portions_j,&
                                                & vvvo_pdm_i,vvvo_pdm_j,&
                                                & ovoo(:,:,i,j),ovoo(:,:,j,i),ovoo(:,:,j,j),&
                                                & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                        call time_start_phase(PHASE_WORK, ttot = time_trip )
                        time_trip_tot = time_trip_tot + time_trip

                        if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(ovoo(:,:,j,k)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                           call time_start_phase(PHASE_WORK, twall = time_efull )
                           call ccsdpt_energy_full_ijk_case2(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                & vvoo(:,:,i,j),vvoo(:,:,j,i),vvoo(:,:,j,j),&
                                                & ccsdpt_singles(:,i),ccsdpt_singles(:,j),&
                                                & e4,async_id,num_ids,cublas_handle)
                           call time_start_phase(PHASE_WORK, ttot = time_efull )
                           time_efull_tot = time_efull_tot + time_efull

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,j,k)) async(async_id(3))

                        else

                           call time_start_phase(PHASE_WORK, twall = time_driv )   
#ifdef VAR_OPENACC
                           call trip_denom_ijk_acc(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else
                           call trip_denom_ijk_cpu(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif

!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))

                           call ccsdpt_driver_ijk_case2(i,j,nocc,nvirt,vvoo(:,:,i,j),vvoo(:,:,j,i),vvoo(:,:,j,j),&
                                                & ovoo(:,:,i,j),ovoo(:,:,j,i),ovoo(:,:,j,j),&
                                                & vvvo_pdm_i,vvvo_pdm_j,&
                                                & ccsdpt_singles(:,i),ccsdpt_singles(:,j),&
                                                & ccsdpt_doubles(:,:,i,j),ccsdpt_doubles(:,:,j,i),&
                                                & ccsdpt_doubles(:,:,j,j),ccsdpt_doubles_2(:,:,:,i),&
                                                & ccsdpt_doubles_2(:,:,:,j),trip_tmp,trip_ampl,&
                                                & async_id,num_ids,cublas_handle)
                           call time_start_phase(PHASE_WORK, ttot = time_driv )
                           time_driv_tot = time_driv_tot + time_driv

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(ovoo(:,:,j,k)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,j,k))&
!$acc& copyout(ccsdpt_doubles(:,:,j,k)) async(async_id(3))

                        endif

                     case(3)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                        call time_start_phase(PHASE_WORK, twall = time_trip )
                        call trip_generator_ijk_case3(i,j,k,nocc,nvirt,ccsd_doubles(:,:,:,i),ccsd_doubles(:,:,:,j),&
                                                & ccsd_doubles(:,:,:,k),ccsd_doubles_portions_i,ccsd_doubles_portions_j,&
                                                & ccsd_doubles_portions_k,vvvo_pdm_i,&
                                                & vvvo_pdm_j,vvvo_pdm_k,&
                                                & ovoo(:,:,i,j),ovoo(:,:,i,k),ovoo(:,:,j,i),&
                                                & ovoo(:,:,j,k),ovoo(:,:,k,i),ovoo(:,:,k,j),&
                                                & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)
                        call time_start_phase(PHASE_WORK, ttot = time_trip )
                        time_trip_tot = time_trip_tot + time_trip

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,k)) async(async_id(1))

                        if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(vvvo_pdm_k,&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i),ovoo(:,:,j,k),ovoo(:,:,k,j)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                           call time_start_phase(PHASE_WORK, twall = time_efull )
                           call ccsdpt_energy_full_ijk_case3(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                & vvoo(:,:,i,j),vvoo(:,:,i,k),vvoo(:,:,j,i),&
                                                & vvoo(:,:,j,k),vvoo(:,:,k,i),vvoo(:,:,k,j),&
                                                & ccsdpt_singles(:,i),ccsdpt_singles(:,j),ccsdpt_singles(:,k),&
                                                & e4,async_id,num_ids,cublas_handle)
                           call time_start_phase(PHASE_WORK, ttot = time_efull )
                           time_efull_tot = time_efull_tot + time_efull

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,k),vvoo(:,:,k,i),vvoo(:,:,j,k),vvoo(:,:,k,j)) async(async_id(3))

                        else

                           call time_start_phase(PHASE_WORK, twall = time_driv )
#ifdef VAR_OPENACC 
                           call trip_denom_ijk_acc(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else
                           call trip_denom_ijk_cpu(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif

!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))

                           call ccsdpt_driver_ijk_case3(i,j,k,nocc,nvirt,vvoo(:,:,i,j),vvoo(:,:,i,k),&
                                                & vvoo(:,:,j,i),vvoo(:,:,j,k),vvoo(:,:,k,i),&
                                                & vvoo(:,:,k,j),ovoo(:,:,i,j),ovoo(:,:,i,k),&
                                                & ovoo(:,:,j,i),ovoo(:,:,j,k),ovoo(:,:,k,i),ovoo(:,:,k,j),&
                                                & vvvo_pdm_i,vvvo_pdm_j,vvvo_pdm_k,&
                                                & ccsdpt_singles(:,i),ccsdpt_singles(:,j),ccsdpt_singles(:,k),&
                                                & ccsdpt_doubles(:,:,i,j),ccsdpt_doubles(:,:,i,k),&
                                                & ccsdpt_doubles(:,:,j,i),ccsdpt_doubles(:,:,j,k),&
                                                & ccsdpt_doubles(:,:,k,i),ccsdpt_doubles(:,:,k,j),&
                                                & ccsdpt_doubles_2(:,:,:,i),ccsdpt_doubles_2(:,:,:,j),&
                                                & ccsdpt_doubles_2(:,:,:,k),trip_tmp,trip_ampl,&
                                                & async_id,num_ids,cublas_handle)
                           call time_start_phase(PHASE_WORK, ttot = time_driv )
                           time_driv_tot = time_driv_tot + time_driv

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(vvvo_pdm_k,&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i),ovoo(:,:,j,k),ovoo(:,:,k,j)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,k),vvoo(:,:,k,i),vvoo(:,:,j,k),vvoo(:,:,k,j))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,k),&
!$acc& ccsdpt_doubles(:,:,i,k),ccsdpt_doubles(:,:,k,i),&
!$acc& ccsdpt_doubles(:,:,j,k),ccsdpt_doubles(:,:,k,j)) async(async_id(3))

                        endif

                     end select TypeOfTuple_par_ijk

                     needed(kbuf) = .false.

                  end do krun_par

                  needed(ibuf) = .false.
                  needed(jbuf) = .false.

               if (j .eq. i) then

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,i)) async(async_id(1))

                  if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(vvvo_pdm_i,&
!$acc& ovoo(:,:,i,j)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,j)) async(async_id(3))

                  else

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(vvvo_pdm_i,&
!$acc& ovoo(:,:,i,j)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,j))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,i),&
!$acc& ccsdpt_doubles(:,:,i,j)) async(async_id(3))

                  endif

               else ! i .gt. j

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,i),ccsd_doubles(:,:,:,j)) async(async_id(1))

                  if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(vvvo_pdm_i,vvvo_pdm_j,&
!$acc& ovoo(:,:,i,j),ovoo(:,:,j,i)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,j),vvoo(:,:,j,i)) async(async_id(3))

                  else

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(vvvo_pdm_i,vvvo_pdm_j,&
!$acc& ovoo(:,:,i,j),ovoo(:,:,j,i)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,j),vvoo(:,:,j,i))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,i),ccsdpt_doubles_2(:,:,:,j),&
!$acc& ccsdpt_doubles(:,:,i,j),ccsdpt_doubles(:,:,j,i)) async(async_id(3))

                  endif

               end if

            end do ijrun_par

    call time_start_phase(PHASE_WORK)

!$acc wait

!$acc exit data delete(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k,&
!$acc& eivalvirt)&
!$acc& copyout(ccsdpt_singles,e4) if(full_no_frags)
!
!$acc exit data delete(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k,&
!$acc& eivalvirt) copyout(ccsdpt_singles) if(.not. full_no_frags)

!$acc wait

    if (alloc_in_dummy) call tensor_unlock_wins(vvvo,all_nodes=.true.)

#ifdef VAR_CUBLAS

    ! Destroy the CUBLAS context
    stat = cublasDestroy_v2(cublas_handle)

#endif

    ! release async handles array
    call mem_dealloc(async_id)

    ! release ccsd_doubles_help_arrays
    call mem_dealloc(ccsd_doubles_portions_i)
    call mem_dealloc(ccsd_doubles_portions_j)
    call mem_dealloc(ccsd_doubles_portions_k)

    ! release pdm work arrays and job list
    call mem_dealloc(vvvo_pdm_buff)
    call mem_dealloc(needed)
    call mem_dealloc(req)
    call mem_dealloc(tiles_in_buf)
    call mem_dealloc(jobs)

    ! release triples ampl structures
    call mem_dealloc(trip_ampl)
    call mem_dealloc(trip_tmp)

    call time_phases_get_diff(current_wt=phase_cntrs)
    call time_start_phase( PHASE_WORK, ttot = time_pt_ijk )

    ! timings
    if (DECinfo%PL .gt. 2) then

       ! minima
       time_pt_ijk_min    = time_pt_ijk
       time_preload_min   = time_preload_tot
       time_trip_min      = time_trip_tot
       time_efull_min     = time_efull_tot
       time_driv_min      = time_driv_tot

       call lsmpi_reduce_realk_min( time_pt_ijk_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_trip_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_efull_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_driv_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_preload_min   , infpar%master, infpar%lg_comm )

       ! maxima
       time_pt_ijk_max    = time_pt_ijk
       time_preload_max   = time_preload_tot
       time_trip_max      = time_trip_tot
       time_efull_max     = time_efull_tot
       time_driv_max      = time_driv_tot

       call lsmpi_reduce_realk_max( time_pt_ijk_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_trip_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_efull_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_driv_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_preload_max   , infpar%master, infpar%lg_comm )

       ! reductions
       call lsmpi_local_reduction( time_pt_ijk   , infpar%master )
       call lsmpi_local_reduction( time_trip_tot   , infpar%master )
       call lsmpi_local_reduction( time_efull_tot   , infpar%master )
       call lsmpi_local_reduction( time_driv_tot   , infpar%master )
       call lsmpi_local_reduction( time_preload_tot   , infpar%master )

       ! unlock, waiting, and flushing
       unlock_time   = time_lsmpi_win_unlock - unlock_time
       waiting_time  = time_lsmpi_wait       - waiting_time
       flushing_time = time_lsmpi_win_flush  - flushing_time
   
       ! minima
       unlock_time_min    = unlock_time
       waiting_time_min   = waiting_time
       flushing_time_min  = flushing_time

       call lsmpi_reduce_realk_min( unlock_time_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( waiting_time_min  , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( flushing_time_min , infpar%master, infpar%lg_comm )

       ! maxima
       unlock_time_max    = unlock_time
       waiting_time_max   = waiting_time
       flushing_time_max  = flushing_time
   
       call lsmpi_reduce_realk_max( unlock_time_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( waiting_time_max  , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( flushing_time_max , infpar%master, infpar%lg_comm )

       ! reductions 
       call lsmpi_local_reduction( unlock_time   , infpar%master )
       call lsmpi_local_reduction( waiting_time  , infpar%master )
       call lsmpi_local_reduction( flushing_time , infpar%master )

       ! work, communication, and idle times
       ! minima
       time_w_min = phase_cntrs( PHASE_WORK_IDX )
       time_c_min = phase_cntrs( PHASE_COMM_IDX )
       time_i_min = phase_cntrs( PHASE_IDLE_IDX )

       call lsmpi_reduce_realk_min( time_w_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_c_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_i_min , infpar%master, infpar%lg_comm )

       ! maxima
       time_w_max = phase_cntrs( PHASE_WORK_IDX )
       time_c_max = phase_cntrs( PHASE_COMM_IDX )
       time_i_max = phase_cntrs( PHASE_IDLE_IDX )
   
       call lsmpi_reduce_realk_max( time_w_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_c_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_i_max , infpar%master, infpar%lg_comm )

       ! reductions 
       call lsmpi_local_reduction(phase_cntrs,nphases,infpar%master)

       if (infpar%lg_mynum .eq. infpar%master) then

          write(*,'("CCSD(T)-ijk time_trip                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_trip_max  ,  time_trip_tot   /dble(nodtotal),time_trip_min  ,  time_trip_tot   / time_pt_ijk
          write(*,'("CCSD(T)-ijk time_preload             ",g10.3,g10.3,g10.3,g10.3)')&
             & time_preload_max  ,  time_preload_tot   /dble(nodtotal),time_preload_min  ,  time_preload_tot   / time_pt_ijk
          write(*,'("CCSD(T)-ijk time_efull               ",g10.3,g10.3,g10.3,g10.3)')&
             & time_efull_max  ,  time_efull_tot   /dble(nodtotal),time_efull_min  ,  time_efull_tot   / time_pt_ijk
          write(*,'("CCSD(T)-ijk time_driv                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_driv_max  ,  time_driv_tot   /dble(nodtotal),time_driv_min  ,  time_driv_tot   / time_pt_ijk
          write(*,'("CCSD(T)-ijk time in lsmpi_win_unlock ",g10.3,g10.3,g10.3,g10.3)')&
             & unlock_time_max, unlock_time/dble(nodtotal),unlock_time_min,          unlock_time   / time_pt_ijk
          write(*,'("CCSD(T)-ijk time in lsmpi_wait       ",g10.3,g10.3,g10.3,g10.3)')&
             & waiting_time_max,   waiting_time  /dble(nodtotal),waiting_time_min,   waiting_time  / time_pt_ijk
          write(*,'("CCSD(T)-ijk time in lsmpi_win_flush  ",g10.3,g10.3,g10.3,g10.3)')&
             & flushing_time_max,  flushing_time /dble(nodtotal),flushing_time_min,  flushing_time / time_pt_ijk
          write(*,'("CCSD(T)-ijk time WORK                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_w_max,phase_cntrs(PHASE_WORK_IDX)/dble(nodtotal),time_w_min,phase_cntrs(PHASE_WORK_IDX)/time_pt_ijk
          write(*,'("CCSD(T)-ijk time COMM                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_c_max,phase_cntrs(PHASE_COMM_IDX)/dble(nodtotal),time_c_min,phase_cntrs(PHASE_COMM_IDX)/time_pt_ijk
          write(*,'("CCSD(T)-ijk time IDLE                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_i_max,phase_cntrs(PHASE_IDLE_IDX)/dble(nodtotal),time_i_min,phase_cntrs(PHASE_IDLE_IDX)/time_pt_ijk

       endif

    endif

    if (infpar%lg_mynum .eq. infpar%master) call LSTIMER('IJK_LOOP_PAR',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine ijk_loop_par
#endif


  subroutine preload_tiles_in_bg_buf(vvvo,jobs,b_size,nocc,current_i,current_j,current_k,current_ij_count,nbuffs,needed,&
        &tiles_in_buf,vvvo_pdm_buff,req)
     implicit none
     type(tensor), intent(inout) :: vvvo
     integer, intent(in) :: b_size,current_i, current_j, current_k, current_ij_count, nbuffs,nocc
     integer, intent(in) :: jobs(b_size+1)
     logical, intent(inout) :: needed(nbuffs)
     integer, intent(inout) :: tiles_in_buf(nbuffs)
     real(realk), pointer, intent(inout) :: vvvo_pdm_buff(:,:)
     integer(kind=ls_mpik),intent(inout) :: req(nbuffs)

     integer :: i_test,j_test,k_test,i_search_buf
     integer :: ibuf_test, jbuf_test, kbuf_test, ij_count_test, ij_test
     logical :: keep_looping,ij_done, new_i_needed, new_j_needed, new_k_needed,found
     integer :: ts
     integer(kind=ls_mpik) :: mode
#ifdef VAR_MPI
     mode = MPI_MODE_NOCHECK

     !set testing integers
     i_test        = current_i
     j_test        = current_j
     k_test        = current_k
     ij_count_test = current_ij_count
     ij_done       = .false.
     keep_looping  = (count(needed)<nbuffs)

     !load next bunch of tiles needed
     fill_buffer: do while(keep_looping)

        !break condition
        keep_looping = (count(needed)<nbuffs.and..not.(ij_done .and. .not. k_test<=j_test))

        if(k_test<=j_test)then

           !Load the next k tile
           call check_if_new_instance_needed(k_test,tiles_in_buf,nbuffs,new_k_needed,set_needed=needed)

           !load k
           if( new_k_needed )then
              !find pos in buff
              call find_free_pos_in_buf(needed,nbuffs,kbuf_test,found)

              if(found)then

                 if( .not.alloc_in_dummy ) call tensor_lock_win(vvvo,k_test,'s',assert=mode)
                 call get_tile_dim(ts,vvvo,k_test)
                 if( alloc_in_dummy )then
                    call tensor_get_tile(vvvo,k_test,vvvo_pdm_buff(:,kbuf_test),ts,&
                       &lock_set=.true.,req=req(kbuf_test))
                 else
                    call tensor_get_tile(vvvo,k_test,vvvo_pdm_buff(:,kbuf_test),ts,&
                       &lock_set=.true.,flush_it=.true.)
                 endif
                 needed(kbuf_test)       = .true.
                 tiles_in_buf(kbuf_test) = k_test

              endif


           endif

        else

           !Load the next i and j tiles

           ij_count_test = ij_count_test + 1

           if(ij_count_test<=b_size)then

              !is incremented by one at the end of the loop
              !therefore we set it to 0 here
              k_test = 0

              ij_test = jobs(ij_count_test)

              call calc_i_leq_j(ij_test,nocc,i_test,j_test)


              call check_if_new_instance_needed(j_test,tiles_in_buf,nbuffs,new_j_needed,set_needed=needed)

              !load new j
              if( new_j_needed )then
                 !find pos in buff
                 call find_free_pos_in_buf(needed,nbuffs,jbuf_test,found)

                 if(found)then

                    if( .not. alloc_in_dummy ) call tensor_lock_win(vvvo,j_test,'s',assert=mode)
                    call get_tile_dim(ts,vvvo,j_test)
                    if(alloc_in_dummy)then
                       call tensor_get_tile(vvvo,j_test,vvvo_pdm_buff(:,jbuf_test),ts,&
                          &lock_set=.true.,req=req(jbuf_test))
                    else
                       call tensor_get_tile(vvvo,j_test,vvvo_pdm_buff(:,jbuf_test),ts,&
                          &lock_set=.true.,flush_it=.true.)
                    endif
                    needed(jbuf_test)       = .true.
                    tiles_in_buf(jbuf_test) = j_test

                 endif

              endif

              call check_if_new_instance_needed(i_test,tiles_in_buf,nbuffs,new_i_needed,set_needed=needed)

              !load new i
              if( new_i_needed )then

                 !find pos in buff
                 call find_free_pos_in_buf(needed,nbuffs,ibuf_test,found)

                 if(found)then

                    if( .not. alloc_in_dummy ) call tensor_lock_win(vvvo,i_test,'s',assert=mode)
                    call get_tile_dim(ts,vvvo,i_test)
                    if(alloc_in_dummy )then
                       call tensor_get_tile(vvvo,i_test,vvvo_pdm_buff(:,ibuf_test),ts,&
                          &lock_set=.true.,req=req(ibuf_test))
                    else
                       call tensor_get_tile(vvvo,i_test,vvvo_pdm_buff(:,ibuf_test),ts,&
                          &lock_set=.true.,flush_it=.true.)
                    endif

                    needed(ibuf_test)       = .true.
                    tiles_in_buf(ibuf_test) = i_test

                 endif

              endif

           else

              ij_done = .true.

           endif


        endif

        k_test = k_test + 1

     enddo fill_buffer
#endif
  end subroutine preload_tiles_in_bg_buf


  !> \brief: main ijk-loop (serial version)
  !> \author: Janus Juul Eriksen
  !> \date: january 2014
  subroutine ijk_loop_ser(nocc,nvirt,ovoo,vvoo,vvvo,ccsd_doubles,&
                        & eivalocc,eivalvirt,ccsdpt_singles,&
                        & ccsdpt_doubles,ccsdpt_doubles_2,e4)

    implicit none

    !> nocc,nvirt
    integer, intent(in) :: nocc,nvirt
    !> 2-el integrals
    real(realk), dimension(nocc,nvirt,nocc,nocc) :: ovoo ! integrals (AI|JK) in the order (J,A,I,K)
    real(realk), dimension(nvirt,nvirt,nocc,nocc) :: vvoo ! integrals (AI|BJ) in the order (A,B,I,J)
    real(realk), dimension(nvirt,nvirt,nvirt,nocc) :: vvvo ! integrals (AI|BC) in the order (C,B,A,I)
    !> ccsd doubles amplitudes
    real(realk), dimension(nvirt,nvirt,nocc,nocc) :: ccsd_doubles
    ! o*v^2 portions of ccsd_doubles
    real(realk), pointer, dimension(:,:,:) :: ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k
    !> triples amplitudes and 3d work array
    real(realk), pointer, dimension(:,:,:) :: trip_tmp, trip_ampl
    !> ccsd(t) intermediates
    real(realk), dimension(nvirt,nocc) :: ccsdpt_singles
    real(realk), dimension(nvirt,nvirt,nocc,nocc),optional :: ccsdpt_doubles
    real(realk), dimension(nocc,nvirt,nvirt,nocc),optional :: ccsdpt_doubles_2
    real(realk),optional :: e4
    logical :: full_no_frags
    !> orbital energiesi
    real(realk), intent(inout)  :: eivalocc(nocc), eivalvirt(nvirt)
    !> loop integers
    integer :: i,j,k,tuple_type
    !> async handles
    integer :: num_ids,m
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind), pointer, dimension(:) :: async_id
    integer(kind=acc_device_kind) :: acc_device_type
#else
    integer, pointer, dimension(:) :: async_id
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat
    real(realk) :: tcpu,twall,norm

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    full_no_frags = .false.

    if (present(e4)) full_no_frags = .true.

    ! init ccsd_doubles_help_arrays
    call mem_alloc(ccsd_doubles_portions_i,nocc,nvirt,nvirt)
    call mem_alloc(ccsd_doubles_portions_j,nocc,nvirt,nvirt)
    call mem_alloc(ccsd_doubles_portions_k,nocc,nvirt,nvirt)

    ! init triples tuples structure
    call mem_alloc(trip_ampl,nvirt,nvirt,nvirt)
    ! init 3d wrk array
    call mem_alloc(trip_tmp,nvirt,nvirt,nvirt)

    ! set async handles. if we are not using gpus, just set them to arbitrary negative numbers
    ! handle 1: ccsd_doubles
    ! handle 2: vvvo and ovoo integrals
    ! handle 3: vvoo integrals and ccsdpt_doubles / ccsdpt_doubles_2 intermediates
    ! handle 4: triples amplitudes
    ! handle 5: energy evaluation 
    num_ids = 5
    call mem_alloc(async_id,num_ids)

#ifdef VAR_OPENACC

    if (DECinfo%acc_sync) then
       async_id = acc_async_sync
    else
       do m = 1,num_ids
          async_id(m) = int(m,kind=acc_handle_kind)
       enddo
    endif

#else

    if (DECinfo%acc_sync) then
       async_id = 0
    else
       do m = 1,num_ids
          async_id(m) = -m
       enddo
    endif

#endif

#ifdef VAR_CUBLAS

    ! initialize the CUBLAS context
    stat = cublasCreate_v2(cublas_handle)

#endif

!$acc wait

!$acc enter data create(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k)&
!$acc& copyin(eivalvirt,ccsdpt_singles,e4) if(full_no_frags)
!
!$acc enter data create(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k)&
!$acc& copyin(eivalvirt,ccsdpt_singles) if(.not. full_no_frags)

!$acc wait

 irun_ser: do i=2,nocc ! i == j == k == 1 gives zero contribution

!$acc enter data copyin(ccsd_doubles(:,:,:,i)) async(async_id(1))

#ifdef VAR_OPENACC
             call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,i),nvirt,nvirt,&
                     & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i,async_id(1))
#else
             call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,i),nvirt,nvirt,&
                     & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i)
#endif

!$acc enter data copyin(vvvo(:,:,:,i)) async(async_id(2))

!$acc enter data copyin(ccsdpt_doubles_2(:,:,:,i)) async(async_id(3)) if(.not. full_no_frags)

    jrun_ser: do j=1,i

                 if (j .eq. i) then 

!$acc enter data copyin(ovoo(:,:,i,j)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,i,j)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,i,j),&
!$acc& ccsdpt_doubles(:,:,i,j)) async(async_id(3)) if(.not. full_no_frags)

                 else ! i .gt. j

!$acc enter data copyin(ccsd_doubles(:,:,:,j)) async(async_id(1))

#ifdef VAR_OPENACC
                    call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,j),nvirt,nvirt,&
                            & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_j,async_id(1))
#else
                    call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,j),nvirt,nvirt,&
                            & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_j)
#endif

!$acc enter data copyin(vvvo(:,:,:,j),&
!$acc& ovoo(:,:,i,j),ovoo(:,:,j,i)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,i,j),vvoo(:,:,j,i)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,i,j),vvoo(:,:,j,i),&
!$acc& ccsdpt_doubles_2(:,:,:,j),&
!$acc& ccsdpt_doubles(:,:,i,j),ccsdpt_doubles(:,:,j,i)) async(async_id(3)) if(.not. full_no_frags)

                 end if

       krun_ser: do k=1,j

                    ! select type of tuple
                    tuple_type = -1

                    if ((i .eq. j) .and. (j .eq. k)) then

                       ! i == j == k
                       ! this always gives zero contribution
                       cycle

                    else if ((i .eq. j) .and. (j .gt. k)) then

                       ! i == j > k
                       tuple_type = 1

!$acc enter data copyin(ccsd_doubles(:,:,:,k)) async(async_id(1))

!$acc enter data copyin(vvvo(:,:,:,k),&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,i,k),vvoo(:,:,k,i)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,i,k),vvoo(:,:,k,i),&
!$acc& ccsdpt_doubles_2(:,:,:,k),&
!$acc& ccsdpt_doubles(:,:,i,k),ccsdpt_doubles(:,:,k,i)) async(async_id(3)) if(.not. full_no_frags)

                    else if ((i .gt. j) .and. (j .eq. k)) then

                       ! i > j == k
                       tuple_type = 2

!$acc enter data copyin(ovoo(:,:,j,k)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,j,k)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,j,k),&
!$acc& ccsdpt_doubles(:,:,j,k)) async(async_id(3)) if(.not. full_no_frags)

                    else

                       ! i > j > k 
                       tuple_type = 3

!$acc enter data copyin(ccsd_doubles(:,:,:,k)) async(async_id(1))

!$acc enter data copyin(vvvo(:,:,:,k),&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i),ovoo(:,:,j,k),ovoo(:,:,k,j)) async(async_id(2))

!$acc enter data copyin(vvoo(:,:,i,k),vvoo(:,:,k,i),vvoo(:,:,j,k),vvoo(:,:,k,j)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(vvoo(:,:,i,k),vvoo(:,:,k,i),vvoo(:,:,j,k),vvoo(:,:,k,j),&
!$acc& ccsdpt_doubles_2(:,:,:,k),&
!$acc& ccsdpt_doubles(:,:,i,k),ccsdpt_doubles(:,:,k,i),&
!$acc& ccsdpt_doubles(:,:,j,k),ccsdpt_doubles(:,:,k,j)) async(async_id(3)) if(.not. full_no_frags)

                    end if

                    if ((tuple_type .eq. 1) .or. (tuple_type .eq. 3)) then

#ifdef VAR_OPENACC
                       call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,k),nvirt,nvirt,&
                               & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_k,async_id(1))
#else
                       call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,k),nvirt,nvirt,&
                               & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_k)
#endif

                    end if
    
                    ! generate tuple(s)
                    TypeOfTuple_ser_ijk: select case(tuple_type)
    
                    case(1)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                       call trip_generator_ijk_case1(i,k,nocc,nvirt,ccsd_doubles(:,:,:,i),ccsd_doubles(:,:,:,k),&
                                               & ccsd_doubles_portions_i,ccsd_doubles_portions_k,&
                                               & vvvo(:,:,:,i),vvvo(:,:,:,k),&
                                               & ovoo(:,:,i,i),ovoo(:,:,i,k),ovoo(:,:,k,i),&
                                               & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,k)) async(async_id(1))

                       if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(vvvo(:,:,:,k),&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                          call ccsdpt_energy_full_ijk_case1(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                               & vvoo(:,:,i,i),vvoo(:,:,i,k),vvoo(:,:,k,i),&
                                               & ccsdpt_singles(:,i),ccsdpt_singles(:,k),&
                                               & e4,async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,k),vvoo(:,:,k,i)) async(async_id(3))
 
                       else

#ifdef VAR_OPENACC 
                          call trip_denom_ijk_acc(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else
                          call trip_denom_ijk_cpu(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif

!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))
 
                          call ccsdpt_driver_ijk_case1(i,k,nocc,nvirt,vvoo(:,:,i,i),vvoo(:,:,i,k),vvoo(:,:,k,i),&
                                               & ovoo(:,:,i,i),ovoo(:,:,i,k),ovoo(:,:,k,i),&
                                               & vvvo(:,:,:,i),vvvo(:,:,:,k),&
                                               & ccsdpt_singles(:,i),ccsdpt_singles(:,k),&
                                               & ccsdpt_doubles(:,:,i,i),ccsdpt_doubles(:,:,i,k),&
                                               & ccsdpt_doubles(:,:,k,i),ccsdpt_doubles_2(:,:,:,i),&
                                               & ccsdpt_doubles_2(:,:,:,k),trip_tmp,trip_ampl,&
                                               & async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(vvvo(:,:,:,k),&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,k),vvoo(:,:,k,i))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,k),&
!$acc& ccsdpt_doubles(:,:,i,k),ccsdpt_doubles(:,:,k,i)) async(async_id(3))

                       endif

                    case(2)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                       call trip_generator_ijk_case2(i,j,nocc,nvirt,ccsd_doubles(:,:,:,i),ccsd_doubles(:,:,:,j),&
                                               & ccsd_doubles_portions_i,ccsd_doubles_portions_j,&
                                               & vvvo(:,:,:,i),vvvo(:,:,:,j),&
                                               & ovoo(:,:,i,j),ovoo(:,:,j,i),ovoo(:,:,j,j),&
                                               & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

                       if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(ovoo(:,:,j,k)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                          call ccsdpt_energy_full_ijk_case2(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                               & vvoo(:,:,i,j),vvoo(:,:,j,i),vvoo(:,:,j,j),&
                                               & ccsdpt_singles(:,i),ccsdpt_singles(:,j),&
                                               & e4,async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,j,k)) async(async_id(3))

                       else

#ifdef VAR_OPENACC       
                          call trip_denom_ijk_acc(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else    
                          call trip_denom_ijk_cpu(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif

!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))

                          call ccsdpt_driver_ijk_case2(i,j,nocc,nvirt,vvoo(:,:,i,j),vvoo(:,:,j,i),vvoo(:,:,j,j),&
                                               & ovoo(:,:,i,j),ovoo(:,:,j,i),ovoo(:,:,j,j),&
                                               & vvvo(:,:,:,i),vvvo(:,:,:,j),&
                                               & ccsdpt_singles(:,i),ccsdpt_singles(:,j),&
                                               & ccsdpt_doubles(:,:,i,j),ccsdpt_doubles(:,:,j,i),&
                                               & ccsdpt_doubles(:,:,j,j),ccsdpt_doubles_2(:,:,:,i),&
                                               & ccsdpt_doubles_2(:,:,:,j),trip_tmp,trip_ampl,&
                                               & async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(ovoo(:,:,j,k)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,j,k))&
!$acc& copyout(ccsdpt_doubles(:,:,j,k)) async(async_id(3))

                       endif

                    case(3)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                       call trip_generator_ijk_case3(i,j,k,nocc,nvirt,ccsd_doubles(:,:,:,i),ccsd_doubles(:,:,:,j),&
                                               & ccsd_doubles(:,:,:,k),ccsd_doubles_portions_i,ccsd_doubles_portions_j,&
                                               & ccsd_doubles_portions_k,vvvo(:,:,:,i),&
                                               & vvvo(:,:,:,j),vvvo(:,:,:,k),&
                                               & ovoo(:,:,i,j),ovoo(:,:,i,k),ovoo(:,:,j,i),&
                                               & ovoo(:,:,j,k),ovoo(:,:,k,i),ovoo(:,:,k,j),&
                                               & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle)

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,k)) async(async_id(1))

                       if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(vvvo(:,:,:,k),&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i),ovoo(:,:,j,k),ovoo(:,:,k,j)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                          call ccsdpt_energy_full_ijk_case3(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                               & vvoo(:,:,i,j),vvoo(:,:,i,k),vvoo(:,:,j,i),&
                                               & vvoo(:,:,j,k),vvoo(:,:,k,i),vvoo(:,:,k,j),&
                                               & ccsdpt_singles(:,i),ccsdpt_singles(:,j),ccsdpt_singles(:,k),&
                                               & e4,async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,k),vvoo(:,:,k,i),vvoo(:,:,j,k),vvoo(:,:,k,j)) async(async_id(3))

                       else

#ifdef VAR_OPENACC       
                          call trip_denom_ijk_acc(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else
                          call trip_denom_ijk_cpu(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif    

!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))

                          call ccsdpt_driver_ijk_case3(i,j,k,nocc,nvirt,vvoo(:,:,i,j),vvoo(:,:,i,k),vvoo(:,:,j,i),&
                                               & vvoo(:,:,j,k),vvoo(:,:,k,i),vvoo(:,:,k,j),ovoo(:,:,i,j),&
                                               & ovoo(:,:,i,k),ovoo(:,:,j,i),ovoo(:,:,j,k),ovoo(:,:,k,i),&
                                               & ovoo(:,:,k,j),vvvo(:,:,:,i),vvvo(:,:,:,j),vvvo(:,:,:,k),&
                                               & ccsdpt_singles(:,i),ccsdpt_singles(:,j),ccsdpt_singles(:,k),&
                                               & ccsdpt_doubles(:,:,i,j),ccsdpt_doubles(:,:,i,k),&
                                               & ccsdpt_doubles(:,:,j,i),ccsdpt_doubles(:,:,j,k),&
                                               & ccsdpt_doubles(:,:,k,i),ccsdpt_doubles(:,:,k,j),&
                                               & ccsdpt_doubles_2(:,:,:,i),ccsdpt_doubles_2(:,:,:,j),&
                                               & ccsdpt_doubles_2(:,:,:,k),trip_tmp,trip_ampl,&
                                               & async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(vvvo(:,:,:,k),&
!$acc& ovoo(:,:,i,k),ovoo(:,:,k,i),ovoo(:,:,j,k),ovoo(:,:,k,j)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,k),vvoo(:,:,k,i),vvoo(:,:,j,k),vvoo(:,:,k,j))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,k),&
!$acc& ccsdpt_doubles(:,:,i,k),ccsdpt_doubles(:,:,k,i),&
!$acc& ccsdpt_doubles(:,:,j,k),ccsdpt_doubles(:,:,k,j)) async(async_id(3))

                       endif

                    end select TypeOfTuple_ser_ijk

                 end do krun_ser

                 if (j .eq. i) then

                    if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(ovoo(:,:,i,j)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,j)) async(async_id(3))

                    else

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(ovoo(:,:,i,j)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,j))&
!$acc& copyout(ccsdpt_doubles(:,:,i,j)) async(async_id(3))

                    endif

                 else ! i .gt. j

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,j)) async(async_id(1))

                    if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(vvvo(:,:,:,j),&
!$acc& ovoo(:,:,i,j),ovoo(:,:,j,i)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,j),vvoo(:,:,j,i)) async(async_id(3))

                    else

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(vvvo(:,:,:,j),&
!$acc& ovoo(:,:,i,j),ovoo(:,:,j,i)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(vvoo(:,:,i,j),vvoo(:,:,j,i))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,j),&
!$acc& ccsdpt_doubles(:,:,i,j),ccsdpt_doubles(:,:,j,i)) async(async_id(3))

                    endif

                 end if
 
              end do jrun_ser

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,i)) async(async_id(1))

              if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(vvvo(:,:,:,i)) async(async_id(2))

              else

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(vvvo(:,:,:,i)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data copyout(ccsdpt_doubles_2(:,:,:,i)) async(async_id(3))

              endif

           end do irun_ser

!$acc wait

!$acc exit data delete(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k,&
!$acc& eivalvirt)&
!$acc& copyout(ccsdpt_singles,e4) if(full_no_frags)
!
!$acc exit data delete(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k,&
!$acc& eivalvirt) copyout(ccsdpt_singles) if(.not. full_no_frags)

!$acc wait

#ifdef VAR_CUBLAS

    ! Destroy the CUBLAS context
    stat = cublasDestroy_v2(cublas_handle)

#endif

    ! release async handles array
    call mem_dealloc(async_id)

    ! release ccsd_doubles_help_arrays
    call mem_dealloc(ccsd_doubles_portions_i)
    call mem_dealloc(ccsd_doubles_portions_j)
    call mem_dealloc(ccsd_doubles_portions_k)

    ! release triples ampl structures
    call mem_dealloc(trip_ampl)
    call mem_dealloc(trip_tmp)

    call LSTIMER('IJK_LOOP_SER',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine ijk_loop_ser


#ifdef VAR_MPI
  !> \brief: main abc-loop (mpi version)
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine abc_loop_par(nocc,nvirt,ooov,oovv,vovv,ccsd_doubles,&
                        & eivalocc,eivalvirt,nodtotal,nbuffs,tile_size,ccsdpt_singles,&
                        & ccsdpt_doubles,ccsdpt_doubles_2,e4)

    implicit none

    !> nocc,nvirt
    integer, intent(in) :: nocc,nvirt
    !> 2-el integrals
    real(realk), dimension(nocc,nocc,nocc,nvirt) :: ooov ! integrals (AI|JK) in the order (K,I,J,A)
    real(realk), dimension(nocc,nocc,nvirt,nvirt) :: oovv ! integrals (AI|BJ) in the order (I,J,A,B)
    type(tensor), intent(inout)  :: vovv ! integrals (AI|BC) in the order (B,I,A,C)
    real(realk), pointer, dimension(:) :: vovv_pdm_a,vovv_pdm_b,vovv_pdm_c ! ov^2*tile_size tiles from vovv
    real(realk), pointer, dimension(:,:) :: vovv_pdm_buff      ! buffers to prefetch tiles
    integer, intent(inout) :: nodtotal, tile_size
    !> ccsd doubles amplitudes
    real(realk), dimension(nocc,nocc,nvirt,nvirt) :: ccsd_doubles
    ! v*o^2 portions of ccsd_doubles
    real(realk), pointer, dimension(:,:,:) :: ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c
    !> triples amplitudes and 3d work array
    real(realk), pointer, dimension(:,:,:) :: trip_tmp, trip_ampl
    !> ccsd(t) intermediates
    real(realk), dimension(nocc,nvirt) :: ccsdpt_singles
    real(realk), dimension(nocc,nocc,nvirt,nvirt), optional :: ccsdpt_doubles
    real(realk), dimension(nvirt,nocc,nocc,nvirt), optional :: ccsdpt_doubles_2
    real(realk),optional :: e4
    logical :: full_no_frags
    !> orbital energiesi
    real(realk), intent(inout)  :: eivalocc(nocc), eivalvirt(nvirt) 
    !> loop integers
    integer :: a,b,c,tuple_type,counter
    integer :: a_tile,b_tile,c_tile,a_pos,b_pos,c_pos,a_count,b_count,c_count
    integer :: a_buf,b_buf,c_buf
    integer :: total_num_tiles,nelms,tile_size_tmp_a,tile_size_tmp_b,tile_size_tmp_c
    !> preloading
    integer, intent(in) :: nbuffs
    integer,pointer :: tiles_in_buf(:)
    integer(kind=ls_mpik), pointer :: req(:)
    logical,pointer :: needed(:)
    !> async handles
    integer :: num_ids,m
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind), pointer, dimension(:) :: async_id
    integer(kind=acc_device_kind) :: acc_device_type
#else
    integer, pointer, dimension(:) :: async_id
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat
    ! timings
    real(realk) :: tcpu,twall,time_pt_abc,time_pt_abc_min,time_pt_abc_max
    real(realk) :: time_trip,time_efull,time_driv,time_preload
    real(realk) :: time_trip_min,time_efull_min,time_driv_min,time_preload_min
    real(realk) :: time_trip_max,time_efull_max,time_driv_max,time_preload_max
    real(realk) :: time_trip_tot,time_efull_tot,time_driv_tot,time_preload_tot
    real(realk) :: phase_cntrs(nphases)
    real(realk) :: flushing_time, flushing_time_min, flushing_time_max
    real(realk) :: unlock_time,   unlock_time_min,   unlock_time_max
    real(realk) :: waiting_time,  waiting_time_min,  waiting_time_max
    real(realk) :: time_w_min, time_w_max
    real(realk) :: time_c_min, time_c_max
    real(realk) :: time_i_min, time_i_max

    ! init timings
    unlock_time   = time_lsmpi_win_unlock
    waiting_time  = time_lsmpi_wait
    flushing_time = time_lsmpi_win_flush
    time_trip_tot = 0.0E0_realk; time_preload_tot = 0.0E0_realk; time_efull_tot = 0.0E0_realk; time_driv_tot = 0.0E0_realk
    call time_start_phase( PHASE_WORK, twall = time_pt_abc )
    call time_phases_get_current(current_wt=phase_cntrs)
    if (infpar%lg_mynum .eq. infpar%master) call LSTIMER('START',tcpu,twall,DECinfo%output)

    full_no_frags = .false.

    if (present(e4)) full_no_frags = .true.

    call time_start_phase(PHASE_WORK)

    ! alloc and init stuff for preloading
    call mem_alloc( vovv_pdm_buff, nocc*nvirt**2*tile_size, nbuffs )
    call mem_alloc( needed,       nbuffs )
    call mem_alloc( tiles_in_buf, nbuffs )
    call mem_alloc( req,          nbuffs )
#ifdef VAR_MPI
    if (alloc_in_dummy) call tensor_lock_wins(vovv,'s',all_nodes=.true.)
#endif
    needed       = .false.
    tiles_in_buf = -1

    ! init ccsd_doubles_help_arrays
    call mem_alloc(ccsd_doubles_portions_a,nvirt,nocc,nocc)
    call mem_alloc(ccsd_doubles_portions_b,nvirt,nocc,nocc)
    call mem_alloc(ccsd_doubles_portions_c,nvirt,nocc,nocc)

    ! init triples tuples structure
    call mem_alloc(trip_ampl,nocc,nocc,nocc)
    ! init 3d wrk array
    call mem_alloc(trip_tmp,nocc,nocc,nocc)

    ! set async handles. if we are not using gpus, just set them to arbitrary negative numbers
    ! handle 1: ccsd_doubles
    ! handle 2: vovv and ooov integrals
    ! handle 3: oovv integrals and ccsdpt_doubles / ccsdpt_doubles_2 intermediates
    ! handle 4: triples amplitudes
    ! handle 5: energy evaluation
    num_ids = 5
    call mem_alloc(async_id,num_ids)

#ifdef VAR_OPENACC

    if (DECinfo%acc_sync) then
       async_id = acc_async_sync
    else
       do m = 1,num_ids
          async_id(m) = int(m,kind=acc_handle_kind)
       enddo
    endif

#else

    if (DECinfo%acc_sync) then
       async_id = 0
    else
       do m = 1,num_ids
          async_id(m) = -m
       enddo
    endif

#endif

#ifdef VAR_CUBLAS

    ! initialize the CUBLAS context
    stat = cublasCreate_v2(cublas_handle)

#endif

    total_num_tiles = vovv%ntiles

    a_count = 0
    b_count = 0
    c_count = 0

    tile_size_tmp_a = 0
    tile_size_tmp_b = 0
    tile_size_tmp_c = 0

    counter = total_num_tiles - infpar%lg_mynum

!$acc wait

!$acc enter data create(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c)&
!$acc& copyin(eivalocc,ccsdpt_singles,e4) if(full_no_frags)
!
!$acc enter data create(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c)&
!$acc& copyin(eivalocc,ccsdpt_singles) if(.not. full_no_frags)

!$acc wait

    do a_tile = total_num_tiles,1,-1

       a_pos = (a_tile - 1) * tile_size + 1

       if (a_tile .ne. counter) then

          cycle

       else

          counter = counter - nodtotal

       endif

       call get_tileinfo_nels_fromarr8(nelms,vovv,i8*a_tile)
       tile_size_tmp_a = int(nelms/(nocc*nvirt**2))

       !FIND a in buffer
       call assoc_ptr_to_buf(a_tile,vovv,nbuffs,tiles_in_buf,needed,vovv_pdm_a,vovv_pdm_buff,a_buf,req)

       call time_start_phase(PHASE_COMM)

       if( alloc_in_dummy )then

          call lsmpi_wait(req(a_buf))

       else

          if(vovv%lock_set(a_tile)) call tensor_unlock_win(vovv,a_tile)

       endif

       needed(a_buf) = .true.

       call time_start_phase(PHASE_WORK)

!!$acc enter data copyin(ptr_pdm_a) async(async_id(2))
!$acc enter data copyin(vovv_pdm_a) async(async_id(2))

       do b_tile = 1,a_tile

          b_pos = (b_tile - 1) * tile_size + 1

          call get_tileinfo_nels_fromarr8(nelms,vovv,i8*b_tile) 
          tile_size_tmp_b = int(nelms/(nocc*nvirt**2))

          !FIND b in buffer
          call assoc_ptr_to_buf(b_tile,vovv,nbuffs,tiles_in_buf,needed,vovv_pdm_b,vovv_pdm_buff,b_buf,req)
   
          call time_start_phase(PHASE_COMM)
   
          if( alloc_in_dummy )then
   
             call lsmpi_wait(req(b_buf))
   
          else
   
             if(vovv%lock_set(b_tile)) call tensor_unlock_win(vovv,b_tile)
   
          endif

          needed(b_buf) = .true.
   
          call time_start_phase(PHASE_WORK)

!!$acc enter data copyin(ptr_pdm_b) async(async_id(2))
!$acc enter data copyin(vovv_pdm_b) async(async_id(2)) if(b_tile .ne. a_tile) 

          do c_tile = 1,b_tile

             c_pos = (c_tile - 1) * tile_size + 1

             call get_tileinfo_nels_fromarr8(nelms,vovv,i8*c_tile)
             tile_size_tmp_c = int(nelms/(nocc*nvirt**2))

             !FIND c in buffer
             call assoc_ptr_to_buf(c_tile,vovv,nbuffs,tiles_in_buf,needed,vovv_pdm_c,vovv_pdm_buff,c_buf,req)
   
             call time_start_phase(PHASE_COMM)
   
             if( alloc_in_dummy )then
   
                call lsmpi_wait(req(c_buf))
   
             else
  
                if(vovv%lock_set(c_tile)) call tensor_unlock_win(vovv,c_tile)
   
             endif

             needed(c_buf) = .true.
   
             call time_start_phase(PHASE_WORK)

!!$acc enter data copyin(ptr_pdm_c) async(async_id(2))
!$acc enter data copyin(vovv_pdm_c) async(async_id(2)) if(c_tile .ne. b_tile)

             call time_start_phase(PHASE_WORK, twall = time_preload )
             call preload_tiles_in_bg_buf_abc(vovv,a_tile,b_tile,c_tile,nbuffs,needed,tiles_in_buf,vovv_pdm_buff,req)
             call time_start_phase(PHASE_WORK, ttot = time_preload )
             time_preload_tot = time_preload_tot + time_preload

! ##########################

             do a = a_pos,a_pos + tile_size_tmp_a - 1

                a_count = a_count + 1

!$acc enter data copyin(ccsd_doubles(:,:,:,a)) async(async_id(1))

#ifdef VAR_OPENACC
                call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,a),nocc,nocc,&
                        & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_a,async_id(1))
#else
                call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,a),nocc,nocc,&
                        & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_a)
#endif

!$acc enter data copyin(ooov(:,:,:,a)) async(async_id(2))

!$acc enter data copyin(ccsdpt_doubles_2(:,:,:,a)) async(async_id(3)) if(.not. full_no_frags)

                do b = b_pos,b_pos + tile_size_tmp_b - 1
         
                   b_count = b_count + 1

                   if (b .gt. a) then

                      b_count = 0
                      cycle

                   endif

                   if (b .eq. a) then

!$acc enter data copyin(oovv(:,:,a,b)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(oovv(:,:,a,b),&
!$acc& ccsdpt_doubles(:,:,a,b)) async(async_id(3)) if(.not. full_no_frags)

                   else ! a .gt. b

!$acc enter data copyin(ccsd_doubles(:,:,:,b)) async(async_id(1))

#ifdef VAR_OPENACC
                      call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,b),nocc,nocc,&
                              & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_b,async_id(1))
#else
                      call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,b),nocc,nocc,&
                              & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_b)
#endif

!$acc enter data copyin(ooov(:,:,:,b)) async(async_id(2))

!$acc enter data copyin(oovv(:,:,a,b),oovv(:,:,b,a)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(oovv(:,:,a,b),oovv(:,:,b,a),&
!$acc& ccsdpt_doubles_2(:,:,:,b),&
!$acc& ccsdpt_doubles(:,:,a,b),ccsdpt_doubles(:,:,b,a)) async(async_id(3)) if(.not. full_no_frags)

                   endif

                   do c = c_pos,c_pos + tile_size_tmp_c - 1

                      c_count = c_count + 1

                      if ((c .gt. b) .or. (c .gt. a)) then

                         c_count = 0
                         cycle
   
                      endif

                      ! select type of tuple
                      tuple_type = -1
         
                      if ((a .eq. b) .and. (b .eq. c)) then
         
                         ! a == b == c
                         ! this always gives zero contribution

                         if (c_count .eq. tile_size_tmp_c) c_count = 0
                         cycle
         
                      else if ((a .eq. b) .and. (b .gt. c)) then
         
                         ! a == b > c
                         tuple_type = 1

!$acc enter data copyin(ccsd_doubles(:,:,:,c)) async(async_id(1))

!$acc enter data copyin(ooov(:,:,:,c)) async(async_id(2))

!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a),&
!$acc& ccsdpt_doubles_2(:,:,:,c),&
!$acc& ccsdpt_doubles(:,:,a,c),ccsdpt_doubles(:,:,c,a)) async(async_id(3)) if(.not. full_no_frags)
         
                      else if ((a .gt. b) .and. (b .eq. c)) then
         
                         ! a > b == c
                         tuple_type = 2

!$acc enter data copyin(oovv(:,:,b,c)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(oovv(:,:,b,c),&
!$acc& ccsdpt_doubles(:,:,b,c)) async(async_id(3)) if(.not. full_no_frags)
         
                      else
         
                         ! a > b > c 
                         tuple_type = 3

!$acc enter data copyin(ccsd_doubles(:,:,:,c)) async(async_id(1))

!$acc enter data copyin(ooov(:,:,:,c)) async(async_id(2))

!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b),&
!$acc& ccsdpt_doubles_2(:,:,:,c),&
!$acc& ccsdpt_doubles(:,:,a,c),ccsdpt_doubles(:,:,c,a),&
!$acc& ccsdpt_doubles(:,:,b,c),ccsdpt_doubles(:,:,c,b)) async(async_id(3)) if(.not. full_no_frags)
         
                      end if
         
                      if ((tuple_type .eq. 1) .or. (tuple_type .eq. 3)) then

#ifdef VAR_OPENACC
                         call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,c),nocc,nocc,&
                                 & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_c,async_id(1))
#else
                         call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,c),nocc,nocc,&
                                 & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_c)
#endif
         
                      end if

                      ! generate tuple(s)
                      TypeOfTuple_par_abc: select case(tuple_type)

                      case(1)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                         call time_start_phase(PHASE_WORK, twall = time_trip )         
                         call trip_generator_abc_case1(a,c,nocc,nvirt,ccsd_doubles(:,:,:,a),ccsd_doubles(:,:,:,c),&
                                                 & ccsd_doubles_portions_a,ccsd_doubles_portions_c,&
                                                 & ooov(:,:,:,a),ooov(:,:,:,c),vovv,&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                                 & a_count,c_count,tile_size_tmp_a,tile_size_tmp_c,&
                                                 & vovv_pdm_a,vovv_pdm_c)
                         call time_start_phase(PHASE_WORK, ttot = time_trip )
                         time_trip_tot = time_trip_tot + time_trip

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,c)) async(async_id(1))

                         if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,c)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                            call time_start_phase(PHASE_WORK, twall = time_efull )
                            call ccsdpt_energy_full_abc_case1(a,a,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                  & oovv(:,:,a,a),oovv(:,:,a,c),oovv(:,:,c,a),&
                                                  & ccsdpt_singles(:,a),ccsdpt_singles(:,c),e4,async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_efull )
                            time_efull_tot = time_efull_tot + time_efull


!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a)) async(async_id(3))
         
                         else

                            call time_start_phase(PHASE_WORK, twall = time_driv )
#ifdef VAR_OPENACC
                            call trip_denom_abc_acc(a,a,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else
                            call trip_denom_abc_cpu(a,a,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif

!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))
            
                            call ccsdpt_driver_abc_case1(a,c,nocc,nvirt,oovv(:,:,a,a),oovv(:,:,a,c),oovv(:,:,c,a),vovv,&
                                                 & ooov(:,:,:,a),ooov(:,:,:,c),&
                                                 & ccsdpt_singles(:,a),ccsdpt_singles(:,c),&
                                                 & ccsdpt_doubles(:,:,a,a),ccsdpt_doubles(:,:,a,c),&
                                                 & ccsdpt_doubles(:,:,c,a),ccsdpt_doubles_2(:,:,:,a),&
                                                 & ccsdpt_doubles_2(:,:,:,c),trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                                 & a_count,c_count,tile_size_tmp_a,tile_size_tmp_c,&
                                                 & vovv_pdm_a,vovv_pdm_c)
                            call time_start_phase(PHASE_WORK, ttot = time_driv )
                            time_driv_tot = time_driv_tot + time_driv

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,c)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,c),&
!$acc& ccsdpt_doubles(:,:,a,c),ccsdpt_doubles(:,:,c,a)) async(async_id(3))

                         endif

                      case(2)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                         call time_start_phase(PHASE_WORK, twall = time_trip )         
                         call trip_generator_abc_case2(a,b,nocc,nvirt,ccsd_doubles(:,:,:,a),ccsd_doubles(:,:,:,b),&
                                                 & ccsd_doubles_portions_a,ccsd_doubles_portions_b,&
                                                 & ooov(:,:,:,a),ooov(:,:,:,b),vovv,&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                                 & a_count,b_count,tile_size_tmp_a,tile_size_tmp_b,&
                                                 & vovv_pdm_a,vovv_pdm_b)
                         call time_start_phase(PHASE_WORK, ttot = time_trip )
                         time_trip_tot = time_trip_tot + time_trip

                         if (full_no_frags) then

! this is different...
!$acc wait(async_id(4)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                            call time_start_phase(PHASE_WORK, twall = time_efull )
                            call ccsdpt_energy_full_abc_case2(a,b,b,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                  & oovv(:,:,a,b),oovv(:,:,b,a),oovv(:,:,b,b),&
                                                  & ccsdpt_singles(:,a),ccsdpt_singles(:,b),e4,async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_efull )
                            time_efull_tot = time_efull_tot + time_efull

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,b,c)) async(async_id(3))

                         else
 
                            call time_start_phase(PHASE_WORK, twall = time_driv )
#ifdef VAR_OPENACC            
                            call trip_denom_abc_acc(a,b,b,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else
                            call trip_denom_abc_cpu(a,b,b,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif
            
!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))
            
                            call ccsdpt_driver_abc_case2(a,b,nocc,nvirt,oovv(:,:,a,b),oovv(:,:,b,a),oovv(:,:,b,b),vovv,&
                                                 & ooov(:,:,:,a),ooov(:,:,:,b),&
                                                 & ccsdpt_singles(:,a),ccsdpt_singles(:,b),&
                                                 & ccsdpt_doubles(:,:,a,b),ccsdpt_doubles(:,:,b,a),&
                                                 & ccsdpt_doubles(:,:,b,b),ccsdpt_doubles_2(:,:,:,a),&
                                                 & ccsdpt_doubles_2(:,:,:,b),trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                                 & a_count,b_count,tile_size_tmp_a,tile_size_tmp_b,&
                                                 & vovv_pdm_a,vovv_pdm_b)
                            call time_start_phase(PHASE_WORK, ttot = time_driv )
                            time_driv_tot = time_driv_tot + time_driv

! this is different...
!$acc wait(async_id(5)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,b,c))&
!$acc& copyout(ccsdpt_doubles(:,:,b,c)) async(async_id(3))
         
                         endif

                      case(3)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                         call time_start_phase(PHASE_WORK, twall = time_trip )         
                         call trip_generator_abc_case3(a,b,c,nocc,nvirt,ccsd_doubles(:,:,:,a),ccsd_doubles(:,:,:,b),&
                                                 & ccsd_doubles(:,:,:,c),ccsd_doubles_portions_a,ccsd_doubles_portions_b,&
                                                 & ccsd_doubles_portions_c,ooov(:,:,:,a),&
                                                 & ooov(:,:,:,b),ooov(:,:,:,c),vovv,&
                                                 & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                                 & a_count,b_count,c_count,tile_size_tmp_a,tile_size_tmp_b,tile_size_tmp_c,&
                                                 & vovv_pdm_a,vovv_pdm_b,vovv_pdm_c)
                         call time_start_phase(PHASE_WORK, ttot = time_trip )
                         time_trip_tot = time_trip_tot + time_trip

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,c)) async(async_id(1))

                         if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,c)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                            call time_start_phase(PHASE_WORK, twall = time_efull )
                            call ccsdpt_energy_full_abc_case3(a,b,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                                  & oovv(:,:,a,b),oovv(:,:,a,c),oovv(:,:,b,a),&
                                                  & oovv(:,:,b,c),oovv(:,:,c,a),oovv(:,:,c,b),&
                                                  & ccsdpt_singles(:,a),ccsdpt_singles(:,b),ccsdpt_singles(:,c),&
                                                  & e4,async_id,num_ids,cublas_handle)
                            call time_start_phase(PHASE_WORK, ttot = time_efull )
                            time_efull_tot = time_efull_tot + time_efull

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b)) async(async_id(3))

                         else

                            call time_start_phase(PHASE_WORK, twall = time_driv ) 
#ifdef VAR_OPENACC            
                            call trip_denom_abc_acc(a,b,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else
                            call trip_denom_abc_cpu(a,b,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif

!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))            
            
                            call ccsdpt_driver_abc_case3(a,b,c,nocc,nvirt,oovv(:,:,a,b),oovv(:,:,a,c),oovv(:,:,b,a),&
                                                 & oovv(:,:,b,c),oovv(:,:,c,a),oovv(:,:,c,b),vovv,&
                                                 & ooov(:,:,:,a),ooov(:,:,:,b),ooov(:,:,:,c),&
                                                 & ccsdpt_singles(:,a),ccsdpt_singles(:,b),ccsdpt_singles(:,c),&
                                                 & ccsdpt_doubles(:,:,a,b),ccsdpt_doubles(:,:,a,c),&
                                                 & ccsdpt_doubles(:,:,b,a),ccsdpt_doubles(:,:,b,c),&
                                                 & ccsdpt_doubles(:,:,c,a),ccsdpt_doubles(:,:,c,b),&
                                                 & ccsdpt_doubles_2(:,:,:,a),ccsdpt_doubles_2(:,:,:,b),&
                                                 & ccsdpt_doubles_2(:,:,:,c),trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                                 & a_count,b_count,c_count,tile_size_tmp_a,tile_size_tmp_b,tile_size_tmp_c,&
                                                 & vovv_pdm_a,vovv_pdm_b,vovv_pdm_c)
                            call time_start_phase(PHASE_WORK, ttot = time_driv )
                            time_driv_tot = time_driv_tot + time_driv

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,c)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,c),&
!$acc& ccsdpt_doubles(:,:,a,c),ccsdpt_doubles(:,:,c,a),&
!$acc& ccsdpt_doubles(:,:,b,c),ccsdpt_doubles(:,:,c,b)) async(async_id(3))
         
                         endif

                      end select TypeOfTuple_par_abc

                      if (c_count .eq. tile_size_tmp_c) c_count = 0
 
                   end do ! end c loop 

                   if (b_count .eq. tile_size_tmp_b) b_count = 0

                   if (b .eq. a) then
         
                      if (full_no_frags) then

! this is different
!$acc wait(async_id(4)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,b)) async(async_id(3))

                      else

! this is different
!$acc wait(async_id(5)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,b))&
!$acc& copyout(ccsdpt_doubles(:,:,a,b)) async(async_id(3))

                      endif
         
                   else ! a .gt. b

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,b)) async(async_id(1))

                      if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,b)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,b),oovv(:,:,b,a)) async(async_id(3))

                      else

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,b)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,b),oovv(:,:,b,a))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,b),&
!$acc& ccsdpt_doubles(:,:,a,b),ccsdpt_doubles(:,:,b,a)) async(async_id(3))

                      endif
         
                   endif

                end do ! end b loop
          
                if (a_count .eq. tile_size_tmp_a) a_count = 0

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,a)) async(async_id(1))

                if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,a)) async(async_id(2))

                else

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,a)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data copyout(ccsdpt_doubles_2(:,:,:,a)) async(async_id(3))

                endif

             end do ! end a loop

! ##########################

!!$acc exit data delete(ptr_pdm_c) async(async_id(2))
!$acc exit data delete(vovv_pdm_c) async(async_id(2))

          needed(c_buf) = .false.

          end do ! end c_tile loop

!!$acc exit data delete(ptr_pdm_b) async(async_id(2))
!$acc exit data delete(vovv_pdm_b) async(async_id(2))

       needed(b_buf) = .false.

       end do ! end b_tile loop

!!$acc exit data delete(ptr_pdm_a) async(async_id(2))
!$acc exit data delete(vovv_pdm_a) async(async_id(2))

    needed(a_buf) = .false.

    end do ! end a_tile loop

    call time_start_phase(PHASE_WORK)

!$acc wait

!$acc exit data delete(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c,&
!$acc& eivalocc)&
!$acc& copyout(ccsdpt_singles,e4) if(full_no_frags)
!
!$acc exit data delete(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c,&
!$acc& eivalocc) copyout(ccsdpt_singles) if(.not. full_no_frags)

!$acc wait

    if (alloc_in_dummy) call tensor_unlock_wins(vovv,all_nodes=.true.)

#ifdef VAR_CUBLAS

    ! Destroy the CUBLAS context
    stat = cublasDestroy_v2 ( cublas_handle )

#endif

    ! release async handles array
    call mem_dealloc(async_id)

    ! release ccsd_doubles_help_arrays
    call mem_dealloc(ccsd_doubles_portions_a)
    call mem_dealloc(ccsd_doubles_portions_b)
    call mem_dealloc(ccsd_doubles_portions_c)

    ! release preloading stuff
    call mem_dealloc(vovv_pdm_buff)
    call mem_dealloc(needed)
    call mem_dealloc(req)
    call mem_dealloc(tiles_in_buf)

    ! release triples ampl structures
    call mem_dealloc(trip_ampl)
    call mem_dealloc(trip_tmp)

    call time_phases_get_diff(current_wt=phase_cntrs)
    call time_start_phase( PHASE_WORK, ttot = time_pt_abc )

    ! timings
    if (DECinfo%PL .gt. 2) then

       ! minima
       time_pt_abc_min    = time_pt_abc
       time_preload_min   = time_preload_tot
       time_trip_min      = time_trip_tot
       time_efull_min     = time_efull_tot
       time_driv_min      = time_driv_tot

       call lsmpi_reduce_realk_min( time_pt_abc_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_trip_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_efull_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_driv_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_preload_min   , infpar%master, infpar%lg_comm )

       ! maxima
       time_pt_abc_max    = time_pt_abc
       time_preload_max   = time_preload_tot
       time_trip_max      = time_trip_tot
       time_efull_max     = time_efull_tot
       time_driv_max      = time_driv_tot

       call lsmpi_reduce_realk_max( time_pt_abc_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_trip_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_efull_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_driv_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_preload_max   , infpar%master, infpar%lg_comm )

       ! reductions
       call lsmpi_local_reduction( time_pt_abc   , infpar%master )
       call lsmpi_local_reduction( time_trip_tot   , infpar%master )
       call lsmpi_local_reduction( time_efull_tot   , infpar%master )
       call lsmpi_local_reduction( time_driv_tot   , infpar%master )
       call lsmpi_local_reduction( time_preload_tot   , infpar%master )

       ! unlock, waiting, and flushing
       unlock_time   = time_lsmpi_win_unlock - unlock_time
       waiting_time  = time_lsmpi_wait       - waiting_time
       flushing_time = time_lsmpi_win_flush  - flushing_time
   
       ! minima
       unlock_time_min    = unlock_time
       waiting_time_min   = waiting_time
       flushing_time_min  = flushing_time

       call lsmpi_reduce_realk_min( unlock_time_min   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( waiting_time_min  , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( flushing_time_min , infpar%master, infpar%lg_comm )

       ! maxima
       unlock_time_max    = unlock_time
       waiting_time_max   = waiting_time
       flushing_time_max  = flushing_time
   
       call lsmpi_reduce_realk_max( unlock_time_max   , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( waiting_time_max  , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( flushing_time_max , infpar%master, infpar%lg_comm )

       ! reductions 
       call lsmpi_local_reduction( unlock_time   , infpar%master )
       call lsmpi_local_reduction( waiting_time  , infpar%master )
       call lsmpi_local_reduction( flushing_time , infpar%master )

       ! work, communication, and idle times
       ! minima
       time_w_min = phase_cntrs( PHASE_WORK_IDX )
       time_c_min = phase_cntrs( PHASE_COMM_IDX )
       time_i_min = phase_cntrs( PHASE_IDLE_IDX )

       call lsmpi_reduce_realk_min( time_w_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_c_min , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_min( time_i_min , infpar%master, infpar%lg_comm )

       ! maxima
       time_w_max = phase_cntrs( PHASE_WORK_IDX )
       time_c_max = phase_cntrs( PHASE_COMM_IDX )
       time_i_max = phase_cntrs( PHASE_IDLE_IDX )
   
       call lsmpi_reduce_realk_max( time_w_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_c_max , infpar%master, infpar%lg_comm )
       call lsmpi_reduce_realk_max( time_i_max , infpar%master, infpar%lg_comm )

       ! reductions 
       call lsmpi_local_reduction(phase_cntrs,nphases,infpar%master)

       if (infpar%lg_mynum .eq. infpar%master) then

          write(*,'("CCSD(T)-abc time_trip                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_trip_max  ,  time_trip_tot   /dble(nodtotal),time_trip_min  ,  time_trip_tot   / time_pt_abc
          write(*,'("CCSD(T)-abc time_preload             ",g10.3,g10.3,g10.3,g10.3)')&
             & time_preload_max  ,  time_preload_tot   /dble(nodtotal),time_preload_min  ,  time_preload_tot   / time_pt_abc
          write(*,'("CCSD(T)-abc time_efull               ",g10.3,g10.3,g10.3,g10.3)')&
             & time_efull_max  ,  time_efull_tot   /dble(nodtotal),time_efull_min  ,  time_efull_tot   / time_pt_abc
          write(*,'("CCSD(T)-abc time_driv                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_driv_max  ,  time_driv_tot   /dble(nodtotal),time_driv_min  ,  time_driv_tot   / time_pt_abc
          write(*,'("CCSD(T)-abc time in lsmpi_win_unlock ",g10.3,g10.3,g10.3,g10.3)')&
             & unlock_time_max, unlock_time/dble(nodtotal),unlock_time_min,          unlock_time   / time_pt_abc
          write(*,'("CCSD(T)-abc time in lsmpi_wait       ",g10.3,g10.3,g10.3,g10.3)')&
             & waiting_time_max,   waiting_time  /dble(nodtotal),waiting_time_min,   waiting_time  / time_pt_abc
          write(*,'("CCSD(T)-abc time in lsmpi_win_flush  ",g10.3,g10.3,g10.3,g10.3)')&
             & flushing_time_max,  flushing_time /dble(nodtotal),flushing_time_min,  flushing_time / time_pt_abc
          write(*,'("CCSD(T)-abc time WORK                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_w_max,phase_cntrs(PHASE_WORK_IDX)/dble(nodtotal),time_w_min,phase_cntrs(PHASE_WORK_IDX)/time_pt_abc
          write(*,'("CCSD(T)-abc time COMM                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_c_max,phase_cntrs(PHASE_COMM_IDX)/dble(nodtotal),time_c_min,phase_cntrs(PHASE_COMM_IDX)/time_pt_abc
          write(*,'("CCSD(T)-abc time IDLE                ",g10.3,g10.3,g10.3,g10.3)')&
             & time_i_max,phase_cntrs(PHASE_IDLE_IDX)/dble(nodtotal),time_i_min,phase_cntrs(PHASE_IDLE_IDX)/time_pt_abc

       endif

    endif

    if (infpar%lg_mynum .eq. infpar%master) call LSTIMER('ABC_LOOP_PAR',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine abc_loop_par
#endif

#ifdef VAR_MPI
  subroutine preload_tiles_in_bg_buf_abc(vovv,cur_a_tile,cur_b_tile,cur_c_tile,nbuffs,needed,tiles_in_buf,vovv_pdm_buff,req)

     implicit none
     type(tensor), intent(inout) :: vovv
     ! current a, b, and c tiles
     integer, intent(in) :: cur_a_tile, cur_b_tile, cur_c_tile
     ! number of buffers
     integer, intent(in) :: nbuffs
     logical, intent(inout) :: needed(nbuffs)
     integer, intent(inout) :: tiles_in_buf(nbuffs)
     real(realk), pointer, intent(inout) :: vovv_pdm_buff(:,:)
     integer(kind=ls_mpik),intent(inout) :: req(nbuffs)

     integer :: a_tile, b_tile, c_tile, a_buf, b_buf, c_buf
     logical :: new_a_needed, new_b_needed, new_c_needed, found
     integer :: ts
     integer(kind=ls_mpik) :: mode

     mode = MPI_MODE_NOCHECK

     a_loop_preload: do a_tile = cur_a_tile,1,-1

        do b_tile = cur_b_tile,a_tile

           do c_tile = cur_c_tile,b_tile

              if (count(needed) .lt. nbuffs) exit a_loop_preload

              ! load the next c tile
              call check_if_new_instance_needed(c_tile,tiles_in_buf,nbuffs,new_c_needed,set_needed=needed)
   
              ! load c
              if( new_c_needed )then
                 ! find pos in buff
                 call find_free_pos_in_buf(needed,nbuffs,c_buf,found)
   
                 if (found) then
   
                    if( .not. alloc_in_dummy ) call tensor_lock_win(vovv,c_tile,'s',assert=mode)
                    call get_tile_dim(ts,vovv,c_tile)
                    if( alloc_in_dummy )then
                       call tensor_get_tile(vovv,c_tile,vovv_pdm_buff(:,c_buf),ts,&
                          &lock_set=.true.,req=req(c_buf))
                    else
                       call tensor_get_tile(vovv,c_tile,vovv_pdm_buff(:,c_buf),ts,&
                          &lock_set=.true.,flush_it=.true.)
                    endif
                    needed(c_buf)       = .true.
                    tiles_in_buf(c_buf) = c_tile
   
                 endif

              endif

           enddo

           if (count(needed) .lt. nbuffs) exit a_loop_preload

           ! load the next b tile
           call check_if_new_instance_needed(b_tile,tiles_in_buf,nbuffs,new_b_needed,set_needed=needed)

           ! load b
           if( new_b_needed )then
              ! find pos in buff
              call find_free_pos_in_buf(needed,nbuffs,b_buf,found)

              if (found) then

                 if( .not. alloc_in_dummy ) call tensor_lock_win(vovv,b_tile,'s',assert=mode)
                 call get_tile_dim(ts,vovv,b_tile)
                 if( alloc_in_dummy )then
                    call tensor_get_tile(vovv,b_tile,vovv_pdm_buff(:,b_buf),ts,&
                       &lock_set=.true.,req=req(b_buf))
                 else
                    call tensor_get_tile(vovv,b_tile,vovv_pdm_buff(:,b_buf),ts,&
                       &lock_set=.true.,flush_it=.true.)
                 endif
                 needed(b_buf)       = .true.
                 tiles_in_buf(b_buf) = b_tile

              endif

           endif

        enddo

        if (count(needed) .lt. nbuffs) exit a_loop_preload

        ! load the next a tile
        call check_if_new_instance_needed(a_tile,tiles_in_buf,nbuffs,new_a_needed,set_needed=needed)

        ! load a
        if( new_a_needed )then
           ! find pos in buff
           call find_free_pos_in_buf(needed,nbuffs,a_buf,found)

           if (found) then

              if( .not. alloc_in_dummy ) call tensor_lock_win(vovv,a_tile,'s',assert=mode)
              call get_tile_dim(ts,vovv,a_tile)
              if( alloc_in_dummy )then
                 call tensor_get_tile(vovv,a_tile,vovv_pdm_buff(:,a_buf),ts,&
                    &lock_set=.true.,req=req(a_buf))
              else
                 call tensor_get_tile(vovv,a_tile,vovv_pdm_buff(:,a_buf),ts,&
                    &lock_set=.true.,flush_it=.true.)
              endif
              needed(a_buf)       = .true.
              tiles_in_buf(a_buf) = a_tile

           endif

        endif

     enddo a_loop_preload

  end subroutine preload_tiles_in_bg_buf_abc
#endif


  !> \brief: main abc-loop (serial version)
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine abc_loop_ser(nocc,nvirt,ooov,oovv,vovv,ccsd_doubles,&
                        & eivalocc,eivalvirt,ccsdpt_singles,&
                        & ccsdpt_doubles,ccsdpt_doubles_2,e4)

    implicit none

    !> nocc,nvirt
    integer, intent(in) :: nocc,nvirt
    !> 2-el integrals
    real(realk), dimension(nocc,nocc,nocc,nvirt) :: ooov ! integrals (AI|JK) in the order (K,I,J,A)
    real(realk), dimension(nocc,nocc,nvirt,nvirt) :: oovv ! integrals (AI|BJ) in the order (I,J,A,B)
    type(tensor), intent(inout) :: vovv ! integrals (AI|BC) in the order (B,I,A,C)
    !> ccsd doubles amplitudes
    real(realk), dimension(nocc,nocc,nvirt,nvirt) :: ccsd_doubles
    ! v*o^2 portions of ccsd_doubles
    real(realk), pointer, dimension(:,:,:) :: ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c
    !> triples amplitudes and 3d work array
    real(realk), pointer, dimension(:,:,:) :: trip_tmp, trip_ampl
    !> ccsd(t) intermediates
    real(realk), dimension(nocc,nvirt) :: ccsdpt_singles
    real(realk), dimension(nocc,nocc,nvirt,nvirt), optional :: ccsdpt_doubles
    real(realk), dimension(nvirt,nocc,nocc,nvirt), optional :: ccsdpt_doubles_2
    real(realk),optional :: e4
    logical :: full_no_frags
    !> orbital energiesi
    real(realk), intent(inout)  :: eivalocc(nocc), eivalvirt(nvirt)
    !> loop integers
    integer :: a,b,c,tuple_type
    !> async handles
    integer :: num_ids,m
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind), pointer, dimension(:) :: async_id
    integer(kind=acc_device_kind) :: acc_device_type
#else
    integer, pointer, dimension(:) :: async_id
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat
    real(realk) :: tcpu,twall

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    full_no_frags = .false.

    if (present(e4)) full_no_frags = .true.

    ! init ccsd_doubles_help_arrays
    call mem_alloc(ccsd_doubles_portions_a,nvirt,nocc,nocc)
    call mem_alloc(ccsd_doubles_portions_b,nvirt,nocc,nocc)
    call mem_alloc(ccsd_doubles_portions_c,nvirt,nocc,nocc)

    ! init triples tuples structure
    call mem_alloc(trip_ampl,nocc,nocc,nocc)
    ! init 3d wrk array
    call mem_alloc(trip_tmp,nocc,nocc,nocc)

    ! set async handles. if we are not using gpus, just set them to arbitrary negative numbers
    ! handle 1: ccsd_doubles
    ! handle 2: vovv and ooov integrals
    ! handle 3: oovv integrals and ccsdpt_doubles / ccsdpt_doubles_2 intermediates
    ! handle 4: triples amplitudes
    ! handle 5: energy evaluation
    num_ids = 5
    call mem_alloc(async_id,num_ids)

#ifdef VAR_OPENACC

    if (DECinfo%acc_sync) then
       async_id = acc_async_sync
    else
       do m = 1,num_ids
          async_id(m) = int(m,kind=acc_handle_kind)
       enddo
    endif

#else

    if (DECinfo%acc_sync) then
       async_id = 0
    else
       do m = 1,num_ids
          async_id(m) = -m
       enddo
    endif

#endif

#ifdef VAR_CUBLAS

    ! initialize the CUBLAS context
    stat = cublasCreate_v2(cublas_handle)

#endif

!$acc wait

!$acc enter data create(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c)&
!$acc& copyin(eivalocc,ccsdpt_singles,e4) if(full_no_frags)
!
!$acc enter data create(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c)&
!$acc& copyin(eivalocc,ccsdpt_singles) if(.not. full_no_frags)

!$acc wait

    do a = 2,nvirt ! a == b == c == 1 gives zero contribution

!$acc enter data copyin(ccsd_doubles(:,:,:,a)) async(async_id(1))

#ifdef VAR_OPENACC
             call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,a),nocc,nocc,&
                     & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_a,async_id(1))
#else
             call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,a),nocc,nocc,&
                     & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_a)
#endif

!$acc enter data copyin(ooov(:,:,:,a)) async(async_id(2))

!$acc enter data copyin(ccsdpt_doubles_2(:,:,:,a)) async(async_id(3)) if(.not. full_no_frags)

       do b=1,a

          if (b .eq. a) then

!$acc enter data copyin(vovv%elm4(:,:,a,b)) async(async_id(2))

!$acc enter data copyin(oovv(:,:,a,b)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(oovv(:,:,a,b),&
!$acc& ccsdpt_doubles(:,:,a,b)) async(async_id(3)) if(.not. full_no_frags)

          else ! a .gt. b

!$acc enter data copyin(ccsd_doubles(:,:,:,b)) async(async_id(1))

#ifdef VAR_OPENACC
             call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,b),nocc,nocc,&
                     & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_b,async_id(1))
#else
             call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,b),nocc,nocc,&
                     & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_b)
#endif

!$acc enter data copyin(ooov(:,:,:,b),&
!$acc& vovv%elm4(:,:,a,b),vovv%elm4(:,:,b,a)) async(async_id(2))

!$acc enter data copyin(oovv(:,:,a,b),oovv(:,:,b,a)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(oovv(:,:,a,b),oovv(:,:,b,a),&
!$acc& ccsdpt_doubles_2(:,:,:,b),&
!$acc& ccsdpt_doubles(:,:,a,b),ccsdpt_doubles(:,:,b,a)) async(async_id(3)) if(.not. full_no_frags)

          endif

          do c=1,b

             ! select type of tuple
             tuple_type = -1

             if ((a .eq. b) .and. (b .eq. c)) then

                ! a == b == c
                ! this always gives zero contribution
                cycle

             else if ((a .eq. b) .and. (b .gt. c)) then

                ! a == b > c
                tuple_type = 1

!$acc enter data copyin(ccsd_doubles(:,:,:,c)) async(async_id(1))

!$acc enter data copyin(ooov(:,:,:,c),&
!$acc& vovv%elm4(:,:,a,c),vovv%elm4(:,:,c,a)) async(async_id(2))

!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a),&
!$acc& ccsdpt_doubles_2(:,:,:,c),&
!$acc& ccsdpt_doubles(:,:,a,c),ccsdpt_doubles(:,:,c,a)) async(async_id(3)) if(.not. full_no_frags)

             else if ((a .gt. b) .and. (b .eq. c)) then

                ! a > b == c
                tuple_type = 2

!$acc enter data copyin(vovv%elm4(:,:,b,c)) async(async_id(2))

!$acc enter data copyin(oovv(:,:,b,c)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(oovv(:,:,b,c),&
!$acc& ccsdpt_doubles(:,:,b,c)) async(async_id(3)) if(.not. full_no_frags)

             else

                ! a > b > c 
                tuple_type = 3

!$acc enter data copyin(ccsd_doubles(:,:,:,c)) async(async_id(1))

!$acc enter data copyin(ooov(:,:,:,c),&
!$acc& vovv%elm4(:,:,a,c),vovv%elm4(:,:,c,a),vovv%elm4(:,:,b,c),vovv%elm4(:,:,c,b)) async(async_id(2))

!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b)) async(async_id(3)) if(full_no_frags)
!
!$acc enter data copyin(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b),&
!$acc& ccsdpt_doubles_2(:,:,:,c),&
!$acc& ccsdpt_doubles(:,:,a,c),ccsdpt_doubles(:,:,c,a),&
!$acc& ccsdpt_doubles(:,:,b,c),ccsdpt_doubles(:,:,c,b)) async(async_id(3)) if(.not. full_no_frags)

             end if

             if ((tuple_type .eq. 1) .or. (tuple_type .eq. 3)) then

#ifdef VAR_OPENACC
                call array_reorder_3d_acc(1.0E0_realk,ccsd_doubles(:,:,:,c),nocc,nocc,&
                        & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_c,async_id(1))
#else
                call array_reorder_3d(1.0E0_realk,ccsd_doubles(:,:,:,c),nocc,nocc,&
                        & nvirt,[3,2,1],0.0E0_realk,ccsd_doubles_portions_c)
#endif

             end if

             ! generate tuple(s)
             TypeOfTuple_ser_abc: select case(tuple_type)

             case(1)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                call trip_generator_abc_case1(a,c,nocc,nvirt,ccsd_doubles(:,:,:,a),ccsd_doubles(:,:,:,c),&
                                        & ccsd_doubles_portions_a,ccsd_doubles_portions_c,&
                                        & ooov(:,:,:,a),ooov(:,:,:,c),vovv,&
                                        & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                        & a,c,1,1)

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,c)) async(async_id(1))

                if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,c),&
!$acc& vovv%elm4(:,:,a,c),vovv%elm4(:,:,c,a)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_energy_full_abc_case1(a,a,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                         & oovv(:,:,a,a),oovv(:,:,a,c),oovv(:,:,c,a),&
                                         & ccsdpt_singles(:,a),ccsdpt_singles(:,c),e4,async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a)) async(async_id(3))

                else

#ifdef VAR_OPENACC   
                   call trip_denom_abc_acc(a,a,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else
                   call trip_denom_abc_cpu(a,a,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif

!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_driver_abc_case1(a,c,nocc,nvirt,oovv(:,:,a,a),oovv(:,:,a,c),oovv(:,:,c,a),vovv,&
                                        & ooov(:,:,:,a),ooov(:,:,:,c),&
                                        & ccsdpt_singles(:,a),ccsdpt_singles(:,c),&
                                        & ccsdpt_doubles(:,:,a,a),ccsdpt_doubles(:,:,a,c),&
                                        & ccsdpt_doubles(:,:,c,a),ccsdpt_doubles_2(:,:,:,a),&
                                        & ccsdpt_doubles_2(:,:,:,c),trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                        & a,c,1,1)

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,c),&
!$acc& vovv%elm4(:,:,a,c),vovv%elm4(:,:,c,a)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,c),&
!$acc& ccsdpt_doubles(:,:,a,c),ccsdpt_doubles(:,:,c,a)) async(async_id(3))

                endif

             case(2)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                call trip_generator_abc_case2(a,b,nocc,nvirt,ccsd_doubles(:,:,:,a),ccsd_doubles(:,:,:,b),&
                                        & ccsd_doubles_portions_a,ccsd_doubles_portions_b,&
                                        & ooov(:,:,:,a),ooov(:,:,:,b),vovv,&
                                        & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                        & a,b,1,1)

                if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(vovv%elm4(:,:,b,c)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_energy_full_abc_case2(a,b,b,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                         & oovv(:,:,a,b),oovv(:,:,b,a),oovv(:,:,b,b),&
                                         & ccsdpt_singles(:,a),ccsdpt_singles(:,b),e4,async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,b,c)) async(async_id(3))

                else

#ifdef VAR_OPENACC   
                   call trip_denom_abc_acc(a,b,b,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else   
                   call trip_denom_abc_cpu(a,b,b,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif

!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_driver_abc_case2(a,b,nocc,nvirt,oovv(:,:,a,b),oovv(:,:,b,a),oovv(:,:,b,b),vovv,&
                                        & ooov(:,:,:,a),ooov(:,:,:,b),&
                                        & ccsdpt_singles(:,a),ccsdpt_singles(:,b),&
                                        & ccsdpt_doubles(:,:,a,b),ccsdpt_doubles(:,:,b,a),&
                                        & ccsdpt_doubles(:,:,b,b),ccsdpt_doubles_2(:,:,:,a),&
                                        & ccsdpt_doubles_2(:,:,:,b),trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                        & a,b,1,1)

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(vovv%elm4(:,:,b,c)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,b,c))&
!$acc& copyout(ccsdpt_doubles(:,:,b,c)) async(async_id(3))

                endif

             case(3)

!$acc wait(async_id(1),async_id(2),async_id(5)) async(async_id(4))

                call trip_generator_abc_case3(a,b,c,nocc,nvirt,ccsd_doubles(:,:,:,a),ccsd_doubles(:,:,:,b),&
                                        & ccsd_doubles(:,:,:,c),ccsd_doubles_portions_a,ccsd_doubles_portions_b,&
                                        & ccsd_doubles_portions_c,ooov(:,:,:,a),&
                                        & ooov(:,:,:,b),ooov(:,:,:,c),vovv,&
                                        & trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                        & a,b,c,1,1,1)

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,c)) async(async_id(1))

                if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,c),&
!$acc& vovv%elm4(:,:,a,c),vovv%elm4(:,:,c,a),vovv%elm4(:,:,b,c),vovv%elm4(:,:,c,b)) async(async_id(2))

!$acc wait(async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_energy_full_abc_case3(a,b,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,trip_tmp,&
                                         & oovv(:,:,a,b),oovv(:,:,a,c),oovv(:,:,b,a),oovv(:,:,b,c),oovv(:,:,c,a),oovv(:,:,c,b),&
                                         & ccsdpt_singles(:,a),ccsdpt_singles(:,b),ccsdpt_singles(:,c),&
                                         & e4,async_id,num_ids,cublas_handle)

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b)) async(async_id(3))

                else

#ifdef VAR_OPENACC  
                   call trip_denom_abc_acc(a,b,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl,async_id(4))
#else
                   call trip_denom_abc_cpu(a,b,c,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
#endif

!$acc wait(async_id(2),async_id(3),async_id(4)) async(async_id(5))

                   call ccsdpt_driver_abc_case3(a,b,c,nocc,nvirt,oovv(:,:,a,b),oovv(:,:,a,c),oovv(:,:,b,a),&
                                        & oovv(:,:,b,c),oovv(:,:,c,a),oovv(:,:,c,b),vovv,&
                                        & ooov(:,:,:,a),ooov(:,:,:,b),ooov(:,:,:,c),&
                                        & ccsdpt_singles(:,a),ccsdpt_singles(:,b),ccsdpt_singles(:,c),&
                                        & ccsdpt_doubles(:,:,a,b),ccsdpt_doubles(:,:,a,c),&
                                        & ccsdpt_doubles(:,:,b,a),ccsdpt_doubles(:,:,b,c),&
                                        & ccsdpt_doubles(:,:,c,a),ccsdpt_doubles(:,:,c,b),&
                                        & ccsdpt_doubles_2(:,:,:,a),ccsdpt_doubles_2(:,:,:,b),&
                                        & ccsdpt_doubles_2(:,:,:,c),trip_tmp,trip_ampl,async_id,num_ids,cublas_handle,&
                                        & a,b,c,1,1,1)

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,c),&
!$acc& vovv%elm4(:,:,a,c),vovv%elm4(:,:,c,a),vovv%elm4(:,:,b,c),vovv%elm4(:,:,c,b)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,c),oovv(:,:,c,a),oovv(:,:,b,c),oovv(:,:,c,b))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,c),&
!$acc& ccsdpt_doubles(:,:,a,c),ccsdpt_doubles(:,:,c,a),&
!$acc& ccsdpt_doubles(:,:,b,c),ccsdpt_doubles(:,:,c,b)) async(async_id(3))

                endif

             end select TypeOfTuple_ser_abc

          end do

          if (b .eq. a) then

             if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(vovv%elm4(:,:,a,b)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,b)) async(async_id(3))

             else

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(vovv%elm4(:,:,a,b)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,b))&
!$acc& copyout(ccsdpt_doubles(:,:,a,b)) async(async_id(3))

             endif

          else ! a .gt. b

!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,b)) async(async_id(1))

             if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,b),&
!$acc& vovv%elm4(:,:,a,b),vovv%elm4(:,:,b,a)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,b),oovv(:,:,b,a)) async(async_id(3))

             else

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,b),&
!$acc& vovv%elm4(:,:,a,b),vovv%elm4(:,:,b,a)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data delete(oovv(:,:,a,b),oovv(:,:,b,a))&
!$acc& copyout(ccsdpt_doubles_2(:,:,:,b),&
!$acc& ccsdpt_doubles(:,:,a,b),ccsdpt_doubles(:,:,b,a)) async(async_id(3))

             endif

          endif

       end do
 
!$acc wait(async_id(4)) async(async_id(1))
!$acc exit data delete(ccsd_doubles(:,:,:,a)) async(async_id(1))

       if (full_no_frags) then

!$acc wait(async_id(4)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,a)) async(async_id(2))

       else

!$acc wait(async_id(5)) async(async_id(2))
!$acc exit data delete(ooov(:,:,:,a)) async(async_id(2))

!$acc wait(async_id(5)) async(async_id(3))
!$acc exit data copyout(ccsdpt_doubles_2(:,:,:,a)) async(async_id(3))

       endif

    end do 

!$acc wait

!$acc exit data delete(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c,&
!$acc& eivalocc)&
!$acc& copyout(ccsdpt_singles,e4) if(full_no_frags)
!
!$acc exit data delete(trip_tmp,trip_ampl,&
!$acc& ccsd_doubles_portions_a,ccsd_doubles_portions_b,ccsd_doubles_portions_c,&
!$acc& eivalocc) copyout(ccsdpt_singles) if(.not. full_no_frags)

!$acc wait

#ifdef VAR_CUBLAS

    ! Destroy the CUBLAS context
    stat = cublasDestroy_v2 ( cublas_handle )

#endif

    ! release async handles array
    call mem_dealloc(async_id)

    ! release ccsd_doubles_help_arrays
    call mem_dealloc(ccsd_doubles_portions_a)
    call mem_dealloc(ccsd_doubles_portions_b)
    call mem_dealloc(ccsd_doubles_portions_c)

    ! release triples ampl structures
    call mem_dealloc(trip_ampl)
    call mem_dealloc(trip_tmp)

    call LSTIMER('ABC_LOOP_SER',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine abc_loop_ser


  subroutine ccsdpt_energy_full_ijk_case1(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & vvoo_tile_12,vvoo_tile_13,vvoo_tile_31,&
                                         & ccsdpt_singles_1,ccsdpt_singles_3,e4,async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: o1,o2,o3,no,nv
    !> trip arrays
    real(realk), dimension(nv,nv,nv), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4
    !> ccsd(t) singles amplitudes
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_3
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_13, vvoo_tile_31
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp e4 energy
    real(realk) :: e4_tmp, e4_tmp1, e4_tmp2, e4_tmp3
    !> ddot
    real(realk), external :: ddot

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! for explanations on the calls to ccsdpt_contract_ijk_11/12,
    ! see the ccsdpt_driver_ijk_case1 routine 

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,i8*nv**3)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,2.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & 2.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case1 - 1) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = 2.0E0_realk * ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_ijk_11(o1,o1,o3,nv,no,vvoo_tile_13,vvoo_tile_31,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_ijk_12(o1,o1,o3,nv,no,vvoo_tile_12,vvoo_tile_12,&
                 & ccsdpt_singles_3,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,-1.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & -1.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case1 - 2) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,-1.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & -1.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case1 - 3) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_ijk_11(o3,o1,o1,nv,no,vvoo_tile_12,vvoo_tile_12,&
                 & ccsdpt_singles_3,trip_ampl,.true.,handle,cublas_handle)
    call ccsdpt_contract_ijk_12(o3,o1,o1,nv,no,vvoo_tile_13,vvoo_tile_31,&
                 & ccsdpt_singles_1,trip_ampl,.true.,handle,cublas_handle)

#ifndef VAR_OPENACC
    e4 = e4 + e4_tmp
#endif

  end subroutine ccsdpt_energy_full_ijk_case1


  subroutine ccsdpt_energy_full_abc_case1(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & oovv_tile_12,oovv_tile_13,oovv_tile_31,&
                                         & ccsdpt_singles_1,ccsdpt_singles_3,e4,async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: v1,v2,v3,no,nv
    !> trip arrays
    real(realk), dimension(no,no,no), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4
    !> ccsd(t) singles amplitudes
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_3
    !> tiles of vvoo integrals
    real(realk), dimension(no,no) :: oovv_tile_12, oovv_tile_13, oovv_tile_31
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp e4 energy
    real(realk) :: e4_tmp, e4_tmp1, e4_tmp2, e4_tmp3
    !> ddot
    real(realk), external :: ddot

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! for explanations on the calls to ccsdpt_contract_abc_11/12,
    ! see the ccsdpt_driver_abc_case1 routine 

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,i8*no**3)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,2.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & 2.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case1 - 1) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = 2.0E0_realk * ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_abc_11(v1,v1,v3,nv,no,oovv_tile_13,oovv_tile_31,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_abc_12(v1,v1,v3,nv,no,oovv_tile_12,oovv_tile_12,&
                 & ccsdpt_singles_3,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,-1.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & -1.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case1 - 2) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,-1.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & -1.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case1 - 3) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_abc_11(v3,v1,v1,nv,no,oovv_tile_12,oovv_tile_12,&
                 & ccsdpt_singles_3,trip_ampl,.true.,handle,cublas_handle)
    call ccsdpt_contract_abc_12(v3,v1,v1,nv,no,oovv_tile_13,oovv_tile_31,&
                 & ccsdpt_singles_1,trip_ampl,.true.,handle,cublas_handle)

#ifndef VAR_OPENACC
    e4 = e4 + e4_tmp
#endif

  end subroutine ccsdpt_energy_full_abc_case1


  subroutine ccsdpt_energy_full_ijk_case2(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & vvoo_tile_12,vvoo_tile_21,vvoo_tile_23,&
                                         & ccsdpt_singles_1,ccsdpt_singles_2,e4,async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: o1,o2,o3,no,nv
    !> trip arrays
    real(realk), dimension(nv,nv,nv), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4
    !> ccsd(t) singles amplitudes
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_2
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_21, vvoo_tile_23
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp e4 energy
    real(realk) :: e4_tmp, e4_tmp1, e4_tmp2, e4_tmp3
    !> ddot
    real(realk), external :: ddot

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! for explanations on the calls to ccsdpt_contract_ijk_11/12,
    ! see the ccsdpt_driver_ijk_case2 routine 

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,i8*nv**3)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,2.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & 2.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case2 - 1) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = 2.0E0_realk * ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_ijk_11(o1,o2,o2,nv,no,vvoo_tile_23,vvoo_tile_23,&
                 & ccsdpt_singles_1,trip_ampl,.true.,handle,cublas_handle)
    call ccsdpt_contract_ijk_12(o1,o2,o2,nv,no,vvoo_tile_21,vvoo_tile_12,&
                 & ccsdpt_singles_2,trip_ampl,.true.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,-1.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & -1.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case2 - 2) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_ijk_11(o2,o2,o1,nv,no,vvoo_tile_21,vvoo_tile_12,&
                 & ccsdpt_singles_2,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_ijk_12(o2,o2,o1,nv,no,vvoo_tile_23,vvoo_tile_23,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,-1.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & -1.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case2 - 3) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

#ifndef VAR_OPENACC
    e4 = e4 + e4_tmp
#endif

  end subroutine ccsdpt_energy_full_ijk_case2


  subroutine ccsdpt_energy_full_abc_case2(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & oovv_tile_12,oovv_tile_21,oovv_tile_23,&
                                         & ccsdpt_singles_1,ccsdpt_singles_2,e4,async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: v1,v2,v3,no,nv
    !> trip arrays
    real(realk), dimension(no,no,no), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4
    !> ccsd(t) singles amplitudes
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_2
    !> tiles of vvoo integrals
    real(realk), dimension(no,no) :: oovv_tile_12, oovv_tile_21, oovv_tile_23
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp e4 energy
    real(realk) :: e4_tmp, e4_tmp1, e4_tmp2, e4_tmp3
    !> ddot
    real(realk), external :: ddot

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! for explanations on the calls to ccsdpt_contract_abc_11/12,
    ! see the ccsdpt_driver_abc_case2 routine 

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,i8*no**3)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,2.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & 2.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case2 - 1) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = 2.0E0_realk * ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_abc_11(v1,v2,v2,nv,no,oovv_tile_23,oovv_tile_23,&
                 & ccsdpt_singles_1,trip_ampl,.true.,handle,cublas_handle)
    call ccsdpt_contract_abc_12(v1,v2,v2,nv,no,oovv_tile_21,oovv_tile_12,&
                 & ccsdpt_singles_2,trip_ampl,.true.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,-1.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & -1.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case2 - 2) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_abc_11(v2,v2,v1,nv,no,oovv_tile_21,oovv_tile_12,&
                 & ccsdpt_singles_2,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_abc_12(v2,v2,v1,nv,no,oovv_tile_23,oovv_tile_23,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,-1.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & -1.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case2 - 3) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

#ifndef VAR_OPENACC
    e4 = e4 + e4_tmp
#endif

  end subroutine ccsdpt_energy_full_abc_case2


  subroutine ccsdpt_energy_full_ijk_case3(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & vvoo_tile_12,vvoo_tile_13,vvoo_tile_21,vvoo_tile_23,vvoo_tile_31,vvoo_tile_32,&
                                         & ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3,e4,async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: o1,o2,o3,no,nv
    !> trip arrays
    real(realk), dimension(nv,nv,nv), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4
    !> ccsd(t) singles amplitudes
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_13, vvoo_tile_21
    real(realk), dimension(nv,nv) :: vvoo_tile_23, vvoo_tile_31, vvoo_tile_32
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp e4 energy
    real(realk) :: e4_tmp, e4_tmp1, e4_tmp2, e4_tmp3, e4_tmp4, e4_tmp5, e4_tmp6
    !> ddot
    real(realk), external :: ddot

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! for explanations on the calls to ccsdpt_contract_ijk_11/12,
    ! see the ccsdpt_driver_ijk_case3 routine 

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,i8*nv**3)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,8.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & 8.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case3 - 1) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = 4.0E0_realk * ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_ijk_11(o1,o2,o3,nv,no,vvoo_tile_23,vvoo_tile_32,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_ijk_12(o1,o2,o3,nv,no,vvoo_tile_21,vvoo_tile_12,&
                 & ccsdpt_singles_3,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[2,3,1],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,2.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & 2.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case3 - 2) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp + ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_ijk_11(o2,o3,o1,nv,no,vvoo_tile_31,vvoo_tile_13,&
                 & ccsdpt_singles_2,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_ijk_12(o2,o3,o1,nv,no,vvoo_tile_32,vvoo_tile_23,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[3,1,2],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,2.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & 2.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case3 - 3) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp + ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_ijk_11(o3,o1,o2,nv,no,vvoo_tile_12,vvoo_tile_21,&
                 & ccsdpt_singles_3,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_ijk_12(o3,o1,o2,nv,no,vvoo_tile_13,vvoo_tile_31,&
                 & ccsdpt_singles_2,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[3,2,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[3,2,1],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,-4.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & -4.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case3 - 4) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - 2.0E0_realk * ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_ijk_11(o3,o2,o1,nv,no,vvoo_tile_21,vvoo_tile_12,&
                 & ccsdpt_singles_3,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_ijk_12(o3,o2,o1,nv,no,vvoo_tile_23,vvoo_tile_32,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[1,3,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[1,3,2],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,-4.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & -4.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case3 - 5) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - 2.0E0_realk * ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_ijk_11(o1,o3,o2,nv,no,vvoo_tile_32,vvoo_tile_23,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_ijk_12(o1,o3,o2,nv,no,vvoo_tile_31,vvoo_tile_13,&
                 & ccsdpt_singles_2,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,[2,1,3],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,[2,1,3],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_ijk_acc(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_ijk_cpu(o1,o2,o3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,nv**3,-4.0E0_realk,trip_tmp,1,trip_ampl,nv**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(nv**3,kind=4),&
                          & -4.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(nv**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_ijk_case3 - 6) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - 2.0E0_realk * ddot(nv**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_ijk_11(o2,o1,o3,nv,no,vvoo_tile_13,vvoo_tile_31,&
                 & ccsdpt_singles_2,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_ijk_12(o2,o1,o3,nv,no,vvoo_tile_12,vvoo_tile_21,&
                 & ccsdpt_singles_3,trip_ampl,.false.,handle,cublas_handle)

#ifndef VAR_OPENACC
    e4 = e4 + 2.0E0_realk * e4_tmp
#endif

  end subroutine ccsdpt_energy_full_ijk_case3


  subroutine ccsdpt_energy_full_abc_case3(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,trip_tmp,&
                                         & oovv_tile_12,oovv_tile_13,oovv_tile_21,oovv_tile_23,oovv_tile_31,oovv_tile_32,&
                                         & ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3,e4,async_idx,num_idxs,cublas_handle)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: v1,v2,v3,no,nv
    !> trip arrays
    real(realk), dimension(no,no,no), target, intent(inout) :: trip_ampl,trip_tmp
    !> orbital energies
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    !> e4 energy
    real(realk), target, intent(inout) :: e4
    !> ccsd(t) singles amplitudes
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3
    !> tiles of vvoo integrals
    real(realk), dimension(no,no) :: oovv_tile_12, oovv_tile_13, oovv_tile_21
    real(realk), dimension(no,no) :: oovv_tile_23, oovv_tile_31, oovv_tile_32
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle
    !> temp e4 energy
    real(realk) :: e4_tmp, e4_tmp1, e4_tmp2, e4_tmp3, e4_tmp4, e4_tmp5, e4_tmp6
    !> ddot
    real(realk), external :: ddot

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! for explanations on the calls to ccsdpt_contract_abc_11/12,
    ! see the ccsdpt_driver_abc_case3 routine 

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC)
    call assign_in_subblocks(trip_tmp,'=',trip_ampl,i8*no**3)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_tmp = trip_ampl
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,8.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & 8.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case3 - 1) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = 4.0E0_realk * ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_abc_11(v1,v2,v3,nv,no,oovv_tile_23,oovv_tile_32,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_abc_12(v1,v2,v3,nv,no,oovv_tile_21,oovv_tile_12,&
                 & ccsdpt_singles_3,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[2,3,1],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,2.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & 2.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case3 - 2) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp + ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_abc_11(v2,v3,v1,nv,no,oovv_tile_31,oovv_tile_13,&
                 & ccsdpt_singles_2,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_abc_12(v2,v3,v1,nv,no,oovv_tile_32,oovv_tile_23,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[3,1,2],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,2.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & 2.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case3 - 3) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp + ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_abc_11(v3,v1,v2,nv,no,oovv_tile_12,oovv_tile_21,&
                 & ccsdpt_singles_3,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_abc_12(v3,v1,v2,nv,no,oovv_tile_13,oovv_tile_31,&
                 & ccsdpt_singles_2,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[3,2,1],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[3,2,1],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,-4.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & -4.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case3 - 4) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - 2.0E0_realk * ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_abc_11(v3,v2,v1,nv,no,oovv_tile_21,oovv_tile_12,&
                 & ccsdpt_singles_3,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_abc_12(v3,v2,v1,nv,no,oovv_tile_23,oovv_tile_32,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[1,3,2],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[1,3,2],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,-4.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & -4.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case3 - 5) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - 2.0E0_realk * ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_abc_11(v1,v3,v2,nv,no,oovv_tile_32,oovv_tile_23,&
                 & ccsdpt_singles_1,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_abc_12(v1,v3,v2,nv,no,oovv_tile_31,oovv_tile_13,&
                 & ccsdpt_singles_2,trip_ampl,.false.,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,[2,1,3],0.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,[2,1,3],0.0E0_realk,trip_ampl)
#endif

#ifdef VAR_OPENACC
    call trip_denom_abc_acc(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl,handle)
#else
    call trip_denom_abc_cpu(v1,v2,v3,no,nv,eigenocc,eigenvirt,trip_ampl)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    call dgemm_acc_openacc_async(handle,'n','n',1,1,no**3,-4.0E0_realk,trip_tmp,1,trip_ampl,no**3,1.0E0_realk,e4,1)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_tmp,trip_ampl,e4)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(1,kind=4),int(1,kind=4),int(no**3,kind=4),&
                          & -4.0E0_realk,c_loc(trip_tmp),int(1,kind=4),c_loc(trip_ampl),int(no**3,kind=4),&
                          & 1.0E0_realk,c_loc(e4),int(1,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_energy_full_abc_case3 - 6) = ",stat
!       stop
!    end if

#endif
#else
    e4_tmp = e4_tmp - 2.0E0_realk * ddot(no**3,trip_tmp,1,trip_ampl,1)
#endif

    call ccsdpt_contract_abc_11(v2,v1,v3,nv,no,oovv_tile_13,oovv_tile_31,&
                 & ccsdpt_singles_2,trip_ampl,.false.,handle,cublas_handle)
    call ccsdpt_contract_abc_12(v2,v1,v3,nv,no,oovv_tile_12,oovv_tile_21,&
                 & ccsdpt_singles_3,trip_ampl,.false.,handle,cublas_handle)

#ifndef VAR_OPENACC
    e4 = e4 + 2.0E0_realk * e4_tmp
#endif

  end subroutine ccsdpt_energy_full_abc_case3


  subroutine ccsdpt_energy_e5_ddot(no,nv,ccsdpt_singles,ccsd_singles,e5)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: no,nv
    !> ccsd(t) and ccsd singles amplitudes
    real(realk), dimension(nv,no), intent(inout) :: ccsdpt_singles,ccsd_singles
    !> e5 energy
    real(realk), intent(inout) :: e5
    !> ddot
    real(realk), external :: ddot

    e5 = 2.0E0_realk * ddot(no*nv,ccsdpt_singles,1,ccsd_singles,1)

  end subroutine ccsdpt_energy_e5_ddot


  !> \brief: create ij_array for ccsd(t)
  !> \author: Janus Juul Eriksen
  !> \date: july 2013
  subroutine create_ij_tensor_ccsdpt(njobs,no,ij_array)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: njobs,no
    !> ij_array
    integer, dimension(njobs), intent(inout) :: ij_array
    !> integers
    integer :: counter,offset,fill_1,fill_2

    ! since i .ge. j, the composite ij indices will make up a lower triangular matrix.
    ! for each ij, k (where j .ge. k) jobs have to be carried out.
    ! thus, the largest jobs for a given i-value will be those that have the largest j-value,
    ! i.e. the largest jobs will be those for which the ij index appears near the diagonal.
    ! as the value of j specifies how large a given job is, we fill up the ij_array with jobs
    ! for j-values in descending order.

    ! the below is the lower triangular part of the ij (5*5) matrix written in row-major order

    ! ||   1              ||
    ! ||   2  3           ||
    ! ||   4  5  6        ||
    ! ||   7  8  9 10     ||
    ! ||  11 12 13 14 15  ||

    ! examples of ij --> i,j conversion
    ! - ij index 15 corresponds to (i,j)=(5,5) and thus to k=1,2,3,4,5
    ! - ij index 9  corresponds to (i,j)=(4,3) and thus to k=1,2,3
    ! - ij index 11  corresponds to (i,j)=(5,1) and thus to k=1

    ! we want ij_array to look like this
    ! (15 , 14 , 10 , 13 , 9 , 6 , 12 , 8 , 5 , 3 , 11 , 7 , 4 , 2 , 1)

    ! counter specifies the index of ij_array
    counter = 1

    do fill_1 = 0,no-1

       ! zero the offset
       offset = 0

       if (fill_1 .eq. 0) then

          ! this is largest possible job, i.e., the (no,no)-th entry in the ij matrix
          ij_array(counter) = njobs
          ! increment counter
          counter = counter + 1

       else

          do fill_2 = 0,fill_1

             if (fill_2 .eq. 0) then

                ! this is the largest i-value, for which we have to do k number of jobs, 
                ! that is, we are at the no'th row essentially moving from right towards left.
                ij_array(counter) = njobs - fill_1
                ! increment counter
                counter = counter + 1

             else

                ! we loop through the i-values keeping the j-value (and k-range) fixed
                ! we thus loop from i == no up towards the diagonal of the lower triangular matrix
                offset = offset + (no - fill_2)
                ! 'njobs - fill_1' gives the current column, while 'offset' moves us up through the rows
                ! while staying below or on the diagonal.(still row-major numbering)
                ij_array(counter) = njobs - fill_1 - offset
                ! increment counter
                counter = counter + 1

             end if

          end do

       end if

    end do

  end subroutine create_ij_tensor_ccsdpt


  !> \brief: make job distribution list for ccsd(t)
  !> \author: Janus Juul Eriksen
  !> \date: july 2013
  subroutine job_distrib_ccsdpt(b_size,njobs,index_array,jobs)

    implicit none

    !> batch size (without remainder contribution) and njobs 
    integer, intent(in) :: b_size,njobs
    !> index_array
    integer, dimension(njobs), intent(inout) :: index_array
    !> jobs array
    integer, dimension(b_size+1), intent(inout) :: jobs
    !> integers
    integer :: nodtotal,fill,fill_sum

#ifdef VAR_MPI

    nodtotal = infpar%lg_nodtot

    ! fill the jobs array with composite index values stored in index_array.
    ! there are njobs jobs in total.

    ! the below algorithm distributes the jobs evenly among the nodes.

    do fill = 0,b_size

       fill_sum = infpar%lg_mynum + 1 + fill*nodtotal

       if (fill_sum .le. njobs) then

          jobs(fill + 1) = index_array(fill_sum) 

       else

          ! fill jobs array with negative number such that this number won't appear for any value of the composite index
          jobs(fill + 1) = -1

       end if

    end do

#endif

  end subroutine job_distrib_ccsdpt



  !> \brief: generator for triples amplitudes, case(1)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_ijk_case1(oindex1,oindex3,no,nv,ccsd_doubles_1,ccsd_doubles_3,&
                                & ccsd_doubles_portions_1,ccsd_doubles_portions_3,&
                                & vvvo_tile_1,vvvo_tile_3,ovoo_tile_11,&
                                & ovoo_tile_13,ovoo_tile_31,trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex3,no,nv
    !> nv**2 tiles of ccsd_doubles
    real(realk), dimension(nv,nv,no) :: ccsd_doubles_1, ccsd_doubles_3
    !> no*nv**2 tiles of ccsd_doubles
    real(realk), dimension(no,nv,nv) :: ccsd_doubles_portions_1,ccsd_doubles_portions_3
    !> tiles of ovoo 2-el integrals
    real(realk), dimension(no,nv) :: ovoo_tile_11, ovoo_tile_13, ovoo_tile_31
    !> tiles of vvvo 2-el integrals
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_1, vvvo_tile_3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs 
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! iik,iki
    call trip_amplitudes_ijk_virt(oindex1,oindex1,oindex3,no,nv,ccsd_doubles_1(:,:,oindex1),&
                            & vvvo_tile_3,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex1,oindex3,oindex1,no,nv,ccsd_doubles_portions_1,&
                            & ovoo_tile_13,trip_tmp,handle,cublas_handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC)
    call assign_in_subblocks(trip_ampl,'=',trip_tmp,i8*nv**3)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_ampl = trip_tmp
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
#endif

    ! kii,iik
    call trip_amplitudes_ijk_virt(oindex3,oindex1,oindex1,no,nv,ccsd_doubles_3(:,:,oindex1),&
                            & vvvo_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex1,oindex1,oindex3,no,nv,ccsd_doubles_portions_1,&
                            & ovoo_tile_31,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

    ! iki.kii
    call trip_amplitudes_ijk_virt(oindex1,oindex3,oindex1,no,nv,ccsd_doubles_1(:,:,oindex3),&
                            & vvvo_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex3,oindex1,oindex1,no,nv,ccsd_doubles_portions_3,&
                            & ovoo_tile_11,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_ijk_case1


  !> \brief: generator for triples amplitudes, case(1)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_abc_case1(vindex1,vindex3,no,nv,ccsd_doubles_1,ccsd_doubles_3,&
                                & ccsd_doubles_portions_1,ccsd_doubles_portions_3,&
                                & ooov_tile_1,ooov_tile_3,vovv,&
                                & trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle,&
                                & a,c,tile_size_a,tile_size_c,vovv_tile_1,vovv_tile_3)

    implicit none

    !> a, b, a, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex3,no,nv
    integer, intent(in) :: a,c,tile_size_a,tile_size_c
    !> no**2 tiles of ccsd_doubles
    real(realk), dimension(no,no,nv) :: ccsd_doubles_1, ccsd_doubles_3
    !> nv*no**2 tiles of ccsd_doubles
    real(realk), dimension(nv,no,no) :: ccsd_doubles_portions_1,ccsd_doubles_portions_3
    !> vovv integrals
    type(tensor), intent(inout)  :: vovv
    !> tiles of vovv 2-el integrals
    real(realk), dimension(nv,no,nv,tile_size_a), intent(inout), optional :: vovv_tile_1
    real(realk), dimension(nv,no,nv,tile_size_c), intent(inout), optional :: vovv_tile_3
    !> tiles of ooov 2-el integrals
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_1, ooov_tile_3
    !> triples amplitude and work array
    real(realk), dimension(no,no,no) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! aac,aca
    call trip_amplitudes_abc_occ(vindex1,vindex1,vindex3,no,nv,ccsd_doubles_1(:,:,vindex1),&
                            & ooov_tile_3,trip_tmp,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex1,vindex3,vindex1,no,nv,ccsd_doubles_portions_1,&
                            & vovv_tile_3(:,:,vindex1,c),trip_tmp,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex1,vindex3,vindex1,no,nv,ccsd_doubles_portions_1,&
                            & vovv%elm4(:,:,vindex1,vindex3),trip_tmp,handle,cublas_handle)
#endif

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC)
    call assign_in_subblocks(trip_ampl,'=',trip_tmp,i8*no**3)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_ampl = trip_tmp
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
#endif

    ! caa,aac
    call trip_amplitudes_abc_occ(vindex3,vindex1,vindex1,no,nv,ccsd_doubles_3(:,:,vindex1),&
                            & ooov_tile_1,trip_tmp,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex1,vindex1,vindex3,no,nv,ccsd_doubles_portions_1,&
                            & vovv_tile_1(:,:,vindex3,a),trip_tmp,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex1,vindex1,vindex3,no,nv,ccsd_doubles_portions_1,&
                            & vovv%elm4(:,:,vindex3,vindex1),trip_tmp,handle,cublas_handle)
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

    ! aca.caa
    call trip_amplitudes_abc_occ(vindex1,vindex3,vindex1,no,nv,ccsd_doubles_1(:,:,vindex3),&
                            & ooov_tile_1,trip_tmp,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex3,vindex1,vindex1,no,nv,ccsd_doubles_portions_3,&
                            & vovv_tile_1(:,:,vindex1,a),trip_tmp,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex3,vindex1,vindex1,no,nv,ccsd_doubles_portions_3,&
                            & vovv%elm4(:,:,vindex1,vindex1),trip_tmp,handle,cublas_handle)
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_abc_case1


  !> \brief: generator for triples amplitudes, case(2)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_ijk_case2(oindex1,oindex2,no,nv,ccsd_doubles_1,ccsd_doubles_2,&
                                & ccsd_doubles_portions_1,ccsd_doubles_portions_2,&
                                & vvvo_tile_1,vvvo_tile_2,ovoo_tile_12,&
                                & ovoo_tile_21,ovoo_tile_22,trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,no,nv
    !> nv**2 tiles of ccsd_doubles
    real(realk), dimension(nv,nv,no) :: ccsd_doubles_1, ccsd_doubles_2
    !> no*nv**2 tiles of ccsd_doubles
    real(realk), dimension(no,nv,nv) :: ccsd_doubles_portions_1,ccsd_doubles_portions_2
    !> tiles of ovoo 2-el integrals
    real(realk), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_21, ovoo_tile_22
    !> tiles of vvvo 2-el integrals
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_1, vvvo_tile_2
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! ijj.jji
    call trip_amplitudes_ijk_virt(oindex1,oindex2,oindex2,no,nv,ccsd_doubles_1(:,:,oindex2),&
                            & vvvo_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex2,oindex2,oindex1,no,nv,ccsd_doubles_portions_2,&
                            & ovoo_tile_12,trip_tmp,handle,cublas_handle)

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC)
    call assign_in_subblocks(trip_ampl,'=',trip_tmp,i8*nv**3)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_ampl = trip_tmp
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
#else 
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
#endif

    ! jij,ijj
    call trip_amplitudes_ijk_virt(oindex2,oindex1,oindex2,no,nv,ccsd_doubles_2(:,:,oindex1),&
                            & vvvo_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex1,oindex2,oindex2,no,nv,ccsd_doubles_portions_1,&
                            & ovoo_tile_22,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
#endif 

    ! jji,jij
    call trip_amplitudes_ijk_virt(oindex2,oindex2,oindex1,no,nv,ccsd_doubles_2(:,:,oindex2),&
                            & vvvo_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex2,oindex1,oindex2,no,nv,ccsd_doubles_portions_2,&
                            & ovoo_tile_21,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_ijk_case2


  !> \brief: generator for triples amplitudes, case(2)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_abc_case2(vindex1,vindex2,no,nv,ccsd_doubles_1,ccsd_doubles_2,&
                                & ccsd_doubles_portions_1,ccsd_doubles_portions_2,&
                                & ooov_tile_1,ooov_tile_2,vovv,&
                                & trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle,&
                                & a,b,tile_size_a,tile_size_b,vovv_tile_1,vovv_tile_2)

    implicit none

    !> a, b, a, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex2,no,nv
    integer, intent(in) :: a,b,tile_size_a,tile_size_b
    !> no**2 tiles of ccsd_doubles
    real(realk), dimension(no,no,nv) :: ccsd_doubles_1, ccsd_doubles_2
    !> nv*no**2 tiles of ccsd_doubles
    real(realk), dimension(nv,no,no) :: ccsd_doubles_portions_1,ccsd_doubles_portions_2
    !> vovv integrals
    type(tensor), intent(inout)  :: vovv
    !> tiles of vovv 2-el integrals
    real(realk), dimension(nv,no,nv,tile_size_a), intent(inout), optional :: vovv_tile_1
    real(realk), dimension(nv,no,nv,tile_size_b), intent(inout), optional :: vovv_tile_2
    !> tiles of ooov 2-el integrals
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_1, ooov_tile_2
    !> triples amplitude and work array
    real(realk), dimension(no,no,no) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! abb.bba
    call trip_amplitudes_abc_occ(vindex1,vindex2,vindex2,no,nv,ccsd_doubles_1(:,:,vindex2),&
                            & ooov_tile_2,trip_tmp,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex2,vindex2,vindex1,no,nv,ccsd_doubles_portions_2,&
                            & vovv_tile_2(:,:,vindex1,b),trip_tmp,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex2,vindex2,vindex1,no,nv,ccsd_doubles_portions_2,&
                            & vovv%elm4(:,:,vindex1,vindex2),trip_tmp,handle,cublas_handle)
#endif

#if defined(VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN) && !defined(VAR_OPENACC)
    call assign_in_subblocks(trip_ampl,'=',trip_tmp,i8*no**3)
#else
!$acc kernels present(trip_ampl,trip_tmp) async(handle)
    trip_ampl = trip_tmp
!$acc end kernels
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
#endif

    ! bab,abb
    call trip_amplitudes_abc_occ(vindex2,vindex1,vindex2,no,nv,ccsd_doubles_2(:,:,vindex1),&
                            & ooov_tile_2,trip_tmp,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex1,vindex2,vindex2,no,nv,ccsd_doubles_portions_1,&
                            & vovv_tile_2(:,:,vindex2,b),trip_tmp,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex1,vindex2,vindex2,no,nv,ccsd_doubles_portions_1,&
                            & vovv%elm4(:,:,vindex2,vindex2),trip_tmp,handle,cublas_handle)
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
#endif

    ! bba,bab
    call trip_amplitudes_abc_occ(vindex2,vindex2,vindex1,no,nv,ccsd_doubles_2(:,:,vindex2),&
                            & ooov_tile_1,trip_tmp,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex2,vindex1,vindex2,no,nv,ccsd_doubles_portions_2,&
                            & vovv_tile_1(:,:,vindex2,a),trip_tmp,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex2,vindex1,vindex2,no,nv,ccsd_doubles_portions_2,&
                            & vovv%elm4(:,:,vindex2,vindex1),trip_tmp,handle,cublas_handle)
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_abc_case2


  !> \brief: generator for triples amplitudes, case(3)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_ijk_case3(oindex1,oindex2,oindex3,no,nv,ccsd_doubles_1,ccsd_doubles_2,&
                                & ccsd_doubles_3,ccsd_doubles_portions_1,ccsd_doubles_portions_2,&
                                & ccsd_doubles_portions_3,vvvo_tile_1,vvvo_tile_2,vvvo_tile_3,&
                                & ovoo_tile_12,ovoo_tile_13,ovoo_tile_21,ovoo_tile_23,ovoo_tile_31,&
                                & ovoo_tile_32,trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,oindex3,no,nv
    !> nv**2 tiles of ccsd_doubles
    real(realk), dimension(nv,nv,no) :: ccsd_doubles_1, ccsd_doubles_2, ccsd_doubles_3
    !> no*nv**2 tiles of ccsd_doubles
    real(realk), dimension(no,nv,nv) :: ccsd_doubles_portions_1,ccsd_doubles_portions_2,ccsd_doubles_portions_3
    !> tiles of ovoo 2-el integrals
    real(realk), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_13, ovoo_tile_21
    real(realk), dimension(no,nv) :: ovoo_tile_23, ovoo_tile_31, ovoo_tile_32
    !> tiles of vvvo 2-el integrals 
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_1, vvvo_tile_2, vvvo_tile_3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! ijk.jki
    call trip_amplitudes_ijk_virt(oindex1,oindex2,oindex3,no,nv,ccsd_doubles_1(:,:,oindex2),&
                            & vvvo_tile_3,trip_ampl,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex2,oindex3,oindex1,no,nv,ccsd_doubles_portions_2,&
                            & ovoo_tile_13,trip_ampl,handle,cublas_handle)

    ! jik,ikj
    call trip_amplitudes_ijk_virt(oindex2,oindex1,oindex3,no,nv,ccsd_doubles_2(:,:,oindex1),&
                            & vvvo_tile_3,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex1,oindex3,oindex2,no,nv,ccsd_doubles_portions_1,&
                            & ovoo_tile_23,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
#endif

    ! kij,ijk
    call trip_amplitudes_ijk_virt(oindex3,oindex1,oindex2,no,nv,ccsd_doubles_3(:,:,oindex1),&
                            & vvvo_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex1,oindex2,oindex3,no,nv,ccsd_doubles_portions_1,&
                            & ovoo_tile_32,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
#endif

    ! jki,kij
    call trip_amplitudes_ijk_virt(oindex2,oindex3,oindex1,no,nv,ccsd_doubles_2(:,:,oindex3),&
                            & vvvo_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex3,oindex1,oindex2,no,nv,ccsd_doubles_portions_3,&
                            & ovoo_tile_21,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
#endif

    ! ikj,kji
    call trip_amplitudes_ijk_virt(oindex1,oindex3,oindex2,no,nv,ccsd_doubles_1(:,:,oindex3),&
                            & vvvo_tile_2,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex3,oindex2,oindex1,no,nv,ccsd_doubles_portions_3,&
                            & ovoo_tile_12,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
#endif 

    ! kji,jik
    call trip_amplitudes_ijk_virt(oindex3,oindex2,oindex1,no,nv,ccsd_doubles_3(:,:,oindex2),&
                            & vvvo_tile_1,trip_tmp,handle,cublas_handle)
    call trip_amplitudes_ijk_occ(oindex2,oindex1,oindex3,no,nv,ccsd_doubles_portions_2,&
                            & ovoo_tile_31,trip_tmp,handle,cublas_handle)

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_ijk_case3


  !> \brief: generator for triples amplitudes, case(3)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_abc_case3(vindex1,vindex2,vindex3,no,nv,ccsd_doubles_1,ccsd_doubles_2,&
                                & ccsd_doubles_3,ccsd_doubles_portions_1,ccsd_doubles_portions_2,&
                                & ccsd_doubles_portions_3,ooov_tile_1,ooov_tile_2,ooov_tile_3,vovv,&
                                & trip_tmp,trip_ampl,async_idx,num_idxs,cublas_handle,&
                                & a,b,c,tile_size_a,tile_size_b,tile_size_c,vovv_tile_1,vovv_tile_2,vovv_tile_3)

    implicit none

    !> a, b, a, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex2,vindex3,no,nv
    integer, intent(in) :: a,b,c,tile_size_a,tile_size_b,tile_size_c
    !> no**2 tiles of ccsd_doubles
    real(realk), dimension(no,no,nv) :: ccsd_doubles_1, ccsd_doubles_2, ccsd_doubles_3
    !> nv*no**2 tiles of ccsd_doubles
    real(realk), dimension(nv,no,no) :: ccsd_doubles_portions_1,ccsd_doubles_portions_2,ccsd_doubles_portions_3
    !> vovv integrals
    type(tensor), intent(inout)  :: vovv
    !> tiles of vovv 2-el integrals
    real(realk), dimension(nv,no,nv,tile_size_a), intent(inout), optional :: vovv_tile_1
    real(realk), dimension(nv,no,nv,tile_size_b), intent(inout), optional :: vovv_tile_2
    real(realk), dimension(nv,no,nv,tile_size_c), intent(inout), optional :: vovv_tile_3
    !> tiles of ooov 2-el integrals 
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_1, ooov_tile_2, ooov_tile_3
    !> triples amplitude and work array
    real(realk), dimension(no,no,no) :: trip_tmp, trip_ampl
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(4)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! abc.bca
    call trip_amplitudes_abc_occ(vindex1,vindex2,vindex3,no,nv,ccsd_doubles_1(:,:,vindex2),&
                            & ooov_tile_3,trip_ampl,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex2,vindex3,vindex1,no,nv,ccsd_doubles_portions_2,&
                            & vovv_tile_3(:,:,vindex1,c),trip_ampl,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex2,vindex3,vindex1,no,nv,ccsd_doubles_portions_2,&
                            & vovv%elm4(:,:,vindex1,vindex3),trip_ampl,handle,cublas_handle)
#endif

    ! cab,abc
    call trip_amplitudes_abc_occ(vindex3,vindex1,vindex2,no,nv,ccsd_doubles_3(:,:,vindex1),&
                            & ooov_tile_2,trip_tmp,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex1,vindex2,vindex3,no,nv,ccsd_doubles_portions_1,&
                            & vovv_tile_2(:,:,vindex3,b),trip_tmp,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex1,vindex2,vindex3,no,nv,ccsd_doubles_portions_1,&
                            & vovv%elm4(:,:,vindex3,vindex2),trip_tmp,handle,cublas_handle)
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
#endif

    ! bca,cab
    call trip_amplitudes_abc_occ(vindex2,vindex3,vindex1,no,nv,ccsd_doubles_2(:,:,vindex3),&
                            & ooov_tile_1,trip_tmp,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex3,vindex1,vindex2,no,nv,ccsd_doubles_portions_3,&
                            & vovv_tile_1(:,:,vindex2,a),trip_tmp,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex3,vindex1,vindex2,no,nv,ccsd_doubles_portions_3,&
                            & vovv%elm4(:,:,vindex2,vindex1),trip_tmp,handle,cublas_handle)
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
#endif

    ! acb,cba
    call trip_amplitudes_abc_occ(vindex1,vindex3,vindex2,no,nv,ccsd_doubles_1(:,:,vindex3),&
                            & ooov_tile_2,trip_tmp,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex3,vindex2,vindex1,no,nv,ccsd_doubles_portions_3,&
                            & vovv_tile_2(:,:,vindex1,b),trip_tmp,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex3,vindex2,vindex1,no,nv,ccsd_doubles_portions_3,&
                            & vovv%elm4(:,:,vindex1,vindex2),trip_tmp,handle,cublas_handle)
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
#endif

    ! bac,acb
    call trip_amplitudes_abc_occ(vindex2,vindex1,vindex3,no,nv,ccsd_doubles_2(:,:,vindex1),&
                            & ooov_tile_3,trip_tmp,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex1,vindex3,vindex2,no,nv,ccsd_doubles_portions_1,&
                            & vovv_tile_3(:,:,vindex2,c),trip_tmp,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex1,vindex3,vindex2,no,nv,ccsd_doubles_portions_1,&
                            & vovv%elm4(:,:,vindex2,vindex3),trip_tmp,handle,cublas_handle)
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
#endif

    ! cba,bac
    call trip_amplitudes_abc_occ(vindex3,vindex2,vindex1,no,nv,ccsd_doubles_3(:,:,vindex2),&
                            & ooov_tile_1,trip_tmp,handle,cublas_handle)
#ifdef VAR_MPI
    call trip_amplitudes_abc_virt(vindex2,vindex1,vindex3,no,nv,ccsd_doubles_portions_2,&
                            & vovv_tile_1(:,:,vindex3,a),trip_tmp,handle,cublas_handle)
#else
    call trip_amplitudes_abc_virt(vindex2,vindex1,vindex3,no,nv,ccsd_doubles_portions_2,&
                            & vovv%elm4(:,:,vindex3,vindex1),trip_tmp,handle,cublas_handle)
#endif

#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl,handle)
#else
    call array_reorder_3d(1.0E0_realk,trip_tmp,no,no,no,&
                        & [3,2,1],1.0E0_realk,trip_ampl)
#endif

  end subroutine trip_generator_abc_case3


  !> \brief: driver routine for contractions in case(1) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_ijk_case1(oindex1,oindex3,no,nv,vvoo_tile_12,vvoo_tile_13,vvoo_tile_31,&
                            & ovoo_tile_12,ovoo_tile_13,ovoo_tile_31,&
                            & vvvo_tile_o1,vvvo_tile_o3,&
                            & ccsdpt_singles_1,ccsdpt_singles_3,&
                            & ccsdpt_doubles_12,ccsdpt_doubles_13,ccsdpt_doubles_31,&
                            & ccsdpt_doubles_2_1,ccsdpt_doubles_2_3,wrk_3d,trip,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex3,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv) :: ccsdpt_doubles_12,ccsdpt_doubles_13,ccsdpt_doubles_31
    real(realk), dimension(no,nv,nv) :: ccsdpt_doubles_2_1,ccsdpt_doubles_2_3
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_3
    !> tiles of ovoo integrals
    real(realk), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_13, ovoo_tile_31
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_13, vvoo_tile_31
    !> tiles of vvvo 2-el integrals determined by incomming oindex1,oindex3
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o1
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : iik --132--> iki --231--> kii
    ! in 211/212 : kii ........ kii ........ iik
    ! in 221/222 : iki ........ iik ........ iki

    do idx = 1,3

       if (idx .eq. 1) then ! iik

          ! calculate contribution to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex1,oindex1,oindex3,nv,no,vvoo_tile_13,vvoo_tile_31,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex1,oindex3,nv,no,vvoo_tile_12,vvoo_tile_12,&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex3,oindex1,oindex1,nv,no,&
                           & ccsdpt_doubles_31,&
                           & wrk_3d,trip,vvvo_tile_o1,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex3,oindex1,oindex1,nv,no,&
                           & ccsdpt_doubles_12,&
                           & wrk_3d,trip,vvvo_tile_o3,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex3,oindex1,no,nv,ovoo_tile_31,ovoo_tile_13,&
                           & ccsdpt_doubles_2_1,trip,.true.,handle,cublas_handle)

       else if (idx .eq. 2) then ! kii

          ! iki: this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_iki trip amplitudes and therefore end up
          ! canceling each other when added to ccsdpt_singles

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_12,&
                           & trip,wrk_3d,vvvo_tile_o3,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_31,&
                           & trip,wrk_3d,vvvo_tile_o1,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex3,oindex1,oindex1,no,nv,ovoo_tile_12,ovoo_tile_12,&
                           & ccsdpt_doubles_2_3,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex3,oindex1,oindex1,no,nv,ovoo_tile_31,&
                           & ccsdpt_doubles_2_1,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 3) then ! iki

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,3,1],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,3,1],0.0E0_realk,trip)
#endif

          call ccsdpt_contract_ijk_11(oindex3,oindex1,oindex1,nv,no,vvoo_tile_12,vvoo_tile_12,&
                       & ccsdpt_singles_3,&
                       & trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex3,oindex1,oindex1,nv,no,vvoo_tile_13,vvoo_tile_31,&
                       & ccsdpt_singles_1,&
                       & trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles_13,&
                           & wrk_3d,trip,vvvo_tile_o1,.true.,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex1,oindex3,no,nv,ovoo_tile_13,ovoo_tile_31,&
                           & ccsdpt_doubles_2_1,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex1,oindex3,no,nv,ovoo_tile_12,&
                           & ccsdpt_doubles_2_3,trip,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_ijk_case1


  !> \brief: driver routine for contractions in case(1) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine ccsdpt_driver_abc_case1(vindex1,vindex3,no,nv,oovv_tile_12,oovv_tile_13,oovv_tile_31,vovv,&
                            & ooov_tile_v1,ooov_tile_v3,&
                            & ccsdpt_singles_1,ccsdpt_singles_3,&
                            & ccsdpt_doubles_12,ccsdpt_doubles_13,ccsdpt_doubles_31,&
                            & ccsdpt_doubles_2_1,ccsdpt_doubles_2_3,wrk_3d,trip,async_idx,num_idxs,cublas_handle,&
                            & a,c,tile_size_a,tile_size_c,vovv_tile_1,vovv_tile_3)

    implicit none

    !> a, b, c, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex3,no,nv
    integer, intent(in) :: a,c,tile_size_a,tile_size_c
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(no,no) :: ccsdpt_doubles_12,ccsdpt_doubles_13,ccsdpt_doubles_31
    real(realk), dimension(nv,no,no) :: ccsdpt_doubles_2_1,ccsdpt_doubles_2_3
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_3
    !> vovv integrals
    type(tensor), intent(inout)  :: vovv
    !> tiles of vovv integrals
    real(realk), dimension(nv,no,nv,tile_size_a), intent(inout), optional :: vovv_tile_1
    real(realk), dimension(nv,no,nv,tile_size_c), intent(inout), optional :: vovv_tile_3 
    !> tiles of oovv integrals
    real(realk), dimension(no,no) :: oovv_tile_12, oovv_tile_13, oovv_tile_31
    !> tiles of ooov 2-el integrals determined by incomming vindex1,vindex3
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v1
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v3
    !> triples amplitude and work array
    real(realk), dimension(no,no,no) :: trip, wrk_3d
    !> loop integer
    integer :: idx
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : aac --132--> aca --231--> caa
    ! in 211/212 : caa ........ caa ........ aac
    ! in 221/222 : aca ........ aac ........ aca

    do idx = 1,3

       if (idx .eq. 1) then ! aac

          ! calculate contribution to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex1,vindex1,vindex3,nv,no,oovv_tile_13,oovv_tile_31,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex1,vindex3,nv,no,oovv_tile_12,oovv_tile_12,&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex1,vindex3,vindex1,no,nv,vovv_tile_1(:,:,vindex3,a),vovv_tile_3(:,:,vindex1,c),&
                           & ccsdpt_doubles_2_1,trip,.true.,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex1,vindex3,vindex1,no,nv,vovv%elm4(:,:,vindex3,vindex1),vovv%elm4(:,:,vindex1,vindex3),&
                           & ccsdpt_doubles_2_1,trip,.true.,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex3,vindex1,vindex1,nv,no,&
                           & ccsdpt_doubles_31,&
                           & wrk_3d,trip,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex3,vindex1,vindex1,nv,no,&
                           & ccsdpt_doubles_12,&
                           & wrk_3d,trip,ooov_tile_v3,handle,cublas_handle)

       else if (idx .eq. 2) then ! caa

          ! aca: this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_aca trip amplitudes and therefore end up
          ! canceling each other when added to ccsdpt_singles

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex3,vindex1,vindex1,no,nv,vovv_tile_1(:,:,vindex1,a),vovv_tile_1(:,:,vindex1,a),&
                           & ccsdpt_doubles_2_3,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex1,vindex1,no,nv,vovv_tile_1(:,:,vindex3,a),&
                           & ccsdpt_doubles_2_1,wrk_3d,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex3,vindex1,vindex1,no,nv,vovv%elm4(:,:,vindex1,vindex1),vovv%elm4(:,:,vindex1,vindex1),&
                           & ccsdpt_doubles_2_3,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex1,vindex1,no,nv,vovv%elm4(:,:,vindex3,vindex1),&
                           & ccsdpt_doubles_2_1,wrk_3d,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_12,&
                           & trip,wrk_3d,ooov_tile_v3,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_31,&
                           & trip,wrk_3d,ooov_tile_v1,handle,cublas_handle)

       else if (idx .eq. 3) then ! aca

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,3,1],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,3,1],0.0E0_realk,trip)
#endif

          call ccsdpt_contract_abc_11(vindex3,vindex1,vindex1,nv,no,oovv_tile_12,oovv_tile_12,&
                       & ccsdpt_singles_3,&
                       & trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex3,vindex1,vindex1,nv,no,oovv_tile_13,oovv_tile_31,&
                       & ccsdpt_singles_1,&
                       & trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex1,vindex1,vindex3,no,nv,vovv_tile_3(:,:,vindex1,c),vovv_tile_1(:,:,vindex3,a),&
                           & ccsdpt_doubles_2_1,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex1,vindex3,no,nv,vovv_tile_1(:,:,vindex1,a),&
                           & ccsdpt_doubles_2_3,trip,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex1,vindex1,vindex3,no,nv,vovv%elm4(:,:,vindex1,vindex3),vovv%elm4(:,:,vindex3,vindex1),&
                           & ccsdpt_doubles_2_1,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex1,vindex3,no,nv,vovv%elm4(:,:,vindex1,vindex1),&
                           & ccsdpt_doubles_2_3,trip,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex3,vindex1,nv,no,&
                           & ccsdpt_doubles_13,&
                           & wrk_3d,trip,ooov_tile_v1,.true.,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_abc_case1


  !> \brief: driver routine for contractions in case(2) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_ijk_case2(oindex1,oindex2,no,nv,vvoo_tile_12,vvoo_tile_21,vvoo_tile_23,&
                            & ovoo_tile_12,ovoo_tile_21,ovoo_tile_23,&
                            & vvvo_tile_o1,vvvo_tile_o2,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,&
                            & ccsdpt_doubles_12,ccsdpt_doubles_21,ccsdpt_doubles_23,&
                            & ccsdpt_doubles_2_1,ccsdpt_doubles_2_2,wrk_3d,trip,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv) :: ccsdpt_doubles_12,ccsdpt_doubles_21,ccsdpt_doubles_23
    real(realk), dimension(no,nv,nv) :: ccsdpt_doubles_2_1,ccsdpt_doubles_2_2
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_2
    !> tiles of ovoo integrals
    real(realk), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_21, ovoo_tile_23
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_21, vvoo_tile_23
    !> tiles of vvvo 2-el integrals determined by incomming oindex1,oindex2
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o1
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o2
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : ijj --312--> jij --312--> jji
    ! in 211/212 : jij ........ jji ........ ijj
    ! in 221/222 : jji ........ ijj ........ jij

    do idx = 1,3

       if (idx .eq. 1) then
  
          ! calculate contributions to ccsdpt_singles:
 
          call ccsdpt_contract_ijk_11(oindex1,oindex2,oindex2,nv,no,vvoo_tile_23,vvoo_tile_23,&
                       & ccsdpt_singles_1,&
                       & trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex2,oindex2,nv,no,vvoo_tile_21,vvoo_tile_12,&
                       & ccsdpt_singles_2,&
                       & trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles_21,&
                           & wrk_3d,trip,vvvo_tile_o2,.true.,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex2,oindex1,no,nv,ovoo_tile_21,ovoo_tile_12,&
                           & ccsdpt_doubles_2_2,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex2,oindex2,oindex1,no,nv,ovoo_tile_23,&
                           & ccsdpt_doubles_2_1,trip,handle,cublas_handle)

       else if (idx .eq. 2) then
   
          ! this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_jij trip amplitudes and therefore end up
          ! canceling each other when added to T_star

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,wrk_3d)
#endif 

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_23,&
                           & trip,wrk_3d,vvvo_tile_o1,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex2,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_12,&
                           & trip,wrk_3d,vvvo_tile_o2,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex2,oindex2,no,nv,ovoo_tile_23,ovoo_tile_23,&
                           & ccsdpt_doubles_2_1,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex2,oindex2,no,nv,ovoo_tile_12,&
                           & ccsdpt_doubles_2_2,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 3) then

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:
   
          call ccsdpt_contract_ijk_11(oindex2,oindex2,oindex1,nv,no,vvoo_tile_21,vvoo_tile_12,&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex2,oindex2,oindex1,nv,no,vvoo_tile_23,vvoo_tile_23,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
   
          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex2,oindex2,nv,no,&
                           & ccsdpt_doubles_12,&
                           & wrk_3d,trip,vvvo_tile_o2,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex2,oindex2,nv,no,&
                           & ccsdpt_doubles_23,&
                           & wrk_3d,trip,vvvo_tile_o1,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex1,oindex2,no,nv,ovoo_tile_12,ovoo_tile_21,&
                           & ccsdpt_doubles_2_2,trip,.true.,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_ijk_case2


  !> \brief: driver routine for contractions in case(2) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine ccsdpt_driver_abc_case2(vindex1,vindex2,no,nv,oovv_tile_12,oovv_tile_21,oovv_tile_23,vovv,&
                            & ooov_tile_v1,ooov_tile_v2,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,&
                            & ccsdpt_doubles_12,ccsdpt_doubles_21,ccsdpt_doubles_23,&
                            & ccsdpt_doubles_2_1,ccsdpt_doubles_2_2,wrk_3d,trip,async_idx,num_idxs,cublas_handle,&
                            & a,b,tile_size_a,tile_size_b,vovv_tile_1,vovv_tile_2)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex2,no,nv
    integer, intent(in) :: a,b,tile_size_a,tile_size_b
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(no,no) :: ccsdpt_doubles_12,ccsdpt_doubles_21,ccsdpt_doubles_23
    real(realk), dimension(nv,no,no) :: ccsdpt_doubles_2_1,ccsdpt_doubles_2_2
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_2
    !> vovv integrals
    type(tensor), intent(inout)  :: vovv
    !> tiles of vovv integrals
    real(realk), dimension(nv,no,nv,tile_size_a), intent(inout), optional :: vovv_tile_1
    real(realk), dimension(nv,no,nv,tile_size_b), intent(inout), optional :: vovv_tile_2
    !> tiles of oovv integrals
    real(realk), dimension(no,no) :: oovv_tile_12, oovv_tile_21, oovv_tile_23
    !> tiles of ooov 2-el integrals determined by incomming vindex1,vindex2
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v1
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v2
    !> triples amplitude and work array
    real(realk), dimension(no,no,no) :: trip, wrk_3d
    !> loop integer
    integer :: idx
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : abb --312--> bab --312--> bba
    ! in 211/212 : bab ........ bba ........ abb
    ! in 221/222 : bba ........ abb ........ bab

    do idx = 1,3

       if (idx .eq. 1) then
  
          ! calculate contributions to ccsdpt_singles:
 
          call ccsdpt_contract_abc_11(vindex1,vindex2,vindex2,nv,no,oovv_tile_23,oovv_tile_23,&
                       & ccsdpt_singles_1,&
                       & trip,.true.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex2,vindex2,nv,no,oovv_tile_21,oovv_tile_12,&
                       & ccsdpt_singles_2,&
                       & trip,.true.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex2,vindex2,vindex1,no,nv,vovv_tile_1(:,:,vindex2,a),vovv_tile_2(:,:,vindex1,b),&
                           & ccsdpt_doubles_2_2,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex2,vindex1,no,nv,vovv_tile_2(:,:,vindex2,b),&
                           & ccsdpt_doubles_2_1,trip,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex2,vindex2,vindex1,no,nv,vovv%elm4(:,:,vindex2,vindex1),vovv%elm4(:,:,vindex1,vindex2),&
                           & ccsdpt_doubles_2_2,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex2,vindex1,no,nv,vovv%elm4(:,:,vindex2,vindex2),&
                           & ccsdpt_doubles_2_1,trip,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex1,vindex2,nv,no,&
                           & ccsdpt_doubles_21,&
                           & wrk_3d,trip,ooov_tile_v2,.true.,handle,cublas_handle)

       else if (idx .eq. 2) then
   
          ! this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_bab trip amplitudes and therefore end up
          ! canceling each other when added to T_star

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[3,1,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[3,1,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex1,vindex2,vindex2,no,nv,vovv_tile_2(:,:,vindex2,b),vovv_tile_2(:,:,vindex2,b),&
                           & ccsdpt_doubles_2_1,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex2,vindex2,no,nv,vovv_tile_2(:,:,vindex1,b),&
                           & ccsdpt_doubles_2_2,wrk_3d,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex1,vindex2,vindex2,no,nv,vovv%elm4(:,:,vindex2,vindex2),vovv%elm4(:,:,vindex2,vindex2),&
                           & ccsdpt_doubles_2_1,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex2,vindex2,no,nv,vovv%elm4(:,:,vindex1,vindex2),&
                           & ccsdpt_doubles_2_2,wrk_3d,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_23,&
                           & trip,wrk_3d,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex2,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_12,&
                           & trip,wrk_3d,ooov_tile_v2,handle,cublas_handle)

       else if (idx .eq. 3) then

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[3,1,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[3,1,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:
   
          call ccsdpt_contract_abc_11(vindex2,vindex2,vindex1,nv,no,oovv_tile_21,oovv_tile_12,&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex2,vindex2,vindex1,nv,no,oovv_tile_23,oovv_tile_23,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
   
          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex2,vindex1,vindex2,no,nv,vovv_tile_2(:,:,vindex1,b),vovv_tile_1(:,:,vindex2,a),&
                           & ccsdpt_doubles_2_2,trip,.true.,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex2,vindex1,vindex2,no,nv,vovv%elm4(:,:,vindex1,vindex2),vovv%elm4(:,:,vindex2,vindex1),&
                           & ccsdpt_doubles_2_2,trip,.true.,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex2,vindex2,nv,no,&
                           & ccsdpt_doubles_12,&
                           & wrk_3d,trip,ooov_tile_v2,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex2,vindex2,nv,no,&
                           & ccsdpt_doubles_23,&
                           & wrk_3d,trip,ooov_tile_v1,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_abc_case2


  !> \brief: driver routine for contractions in case(3) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_ijk_case3(oindex1,oindex2,oindex3,no,nv,&
                            & vvoo_tile_12, vvoo_tile_13, vvoo_tile_21,&
                            & vvoo_tile_23, vvoo_tile_31, vvoo_tile_32,&
                            & ovoo_tile_12, ovoo_tile_13, ovoo_tile_21,&
                            & ovoo_tile_23, ovoo_tile_31, ovoo_tile_32,&
                            & vvvo_tile_o1,vvvo_tile_o2,vvvo_tile_o3,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3,&
                            & ccsdpt_doubles_12,ccsdpt_doubles_13,ccsdpt_doubles_21,&
                            & ccsdpt_doubles_23,ccsdpt_doubles_31,ccsdpt_doubles_32,&
                            & ccsdpt_doubles_2_1,ccsdpt_doubles_2_2,ccsdpt_doubles_2_3,&
                            & wrk_3d,trip,async_idx,num_idxs,cublas_handle)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,oindex3,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv) :: ccsdpt_doubles_12,ccsdpt_doubles_13,ccsdpt_doubles_21
    real(realk), dimension(nv,nv) :: ccsdpt_doubles_23,ccsdpt_doubles_31,ccsdpt_doubles_32
    real(realk), dimension(no,nv,nv) :: ccsdpt_doubles_2_1,ccsdpt_doubles_2_2,ccsdpt_doubles_2_3
    real(realk), dimension(nv) :: ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3
    !> tiles of ovoo integrals
    real(realk), dimension(no,nv) :: ovoo_tile_12, ovoo_tile_13, ovoo_tile_21
    real(realk), dimension(no,nv) :: ovoo_tile_23, ovoo_tile_31, ovoo_tile_32
    !> tiles of vvoo integrals
    real(realk), dimension(nv,nv) :: vvoo_tile_12, vvoo_tile_13, vvoo_tile_21
    real(realk), dimension(nv,nv) :: vvoo_tile_23, vvoo_tile_31, vvoo_tile_32
    !> tiles of vvvo 2-el integrals determined by incomming oindex1,oindex2,oindex3
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o1
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o2
    real(realk), dimension(nv,nv,nv), intent(inout) :: vvvo_tile_o3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. 
    ! in order to minimize the number of reorderings needed to be performed, 
    ! in order to take optimal advantage of the symmetry of the amplitudes, 
    ! AND to finish all work that involves int_virt_tile_o3 first such that this may be updated (gpu),
    ! we carry out the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : ijk --213--> jik --132--> jki --321--> ikj --213--> kij --132--> kji
    ! in 211/212 : kij ........ ikj ........ ijk ........ kji ........ jki ........ jik
    ! in 221/222 : jki ........ kji ........ kij ........ jik ........ ijk ........ ikj

    do idx = 1,6

       if (idx .eq. 1) then ! ijk

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex1,oindex2,oindex3,nv,no,vvoo_tile_23,vvoo_tile_32,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex2,oindex3,nv,no,vvoo_tile_21,vvoo_tile_12,&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex3,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles_31,&
                           & wrk_3d,trip,vvvo_tile_o2,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex3,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles_21,&
                           & wrk_3d,trip,vvvo_tile_o3,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex3,oindex1,no,nv,ovoo_tile_31,ovoo_tile_13,&
                           & ccsdpt_doubles_2_2,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex2,oindex3,oindex1,no,nv,ovoo_tile_23,&
                           & ccsdpt_doubles_2_1,trip,handle,cublas_handle)

       else if (idx .eq. 2) then ! kij

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex2,oindex1,oindex3,nv,no,vvoo_tile_13,vvoo_tile_31,&
                       & ccsdpt_singles_2,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex2,oindex1,oindex3,nv,no,vvoo_tile_12,vvoo_tile_21,&
                       & ccsdpt_singles_3,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex3,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_32,&
                           & trip,wrk_3d,vvvo_tile_o1,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex3,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles_12,&
                           & trip,wrk_3d,vvvo_tile_o3,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex3,oindex2,no,nv,ovoo_tile_32,ovoo_tile_23,&
                           & ccsdpt_doubles_2_1,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex3,oindex2,no,nv,ovoo_tile_13,&
                           & ccsdpt_doubles_2_2,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 3) then ! jki

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex2,oindex3,oindex1,nv,no,vvoo_tile_31,vvoo_tile_13,&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex2,oindex3,oindex1,nv,no,vvoo_tile_32,vvoo_tile_23,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex2,oindex3,nv,no,&
                           & ccsdpt_doubles_12,&
                           & wrk_3d,trip,vvvo_tile_o3,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex2,oindex3,nv,no,&
                           & ccsdpt_doubles_32,&
                           & wrk_3d,trip,vvvo_tile_o1,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex3,oindex1,oindex2,no,nv,ovoo_tile_12,ovoo_tile_21,&
                           & ccsdpt_doubles_2_3,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex3,oindex1,oindex2,no,nv,ovoo_tile_31,&
                           & ccsdpt_doubles_2_2,trip,handle,cublas_handle)

       else if (idx .eq. 4) then ! kji

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,2,1],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,2,1],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex1,oindex3,oindex2,nv,no,vvoo_tile_32,vvoo_tile_23,&
                       & ccsdpt_singles_1,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex1,oindex3,oindex2,nv,no,vvoo_tile_31,vvoo_tile_13,&
                       & ccsdpt_singles_2,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_21,&
                           & trip,wrk_3d,vvvo_tile_o3,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex2,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles_31,&
                           & trip,wrk_3d,vvvo_tile_o2,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex3,oindex2,oindex1,no,nv,ovoo_tile_21,ovoo_tile_12,&
                           & ccsdpt_doubles_2_3,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex3,oindex2,oindex1,no,nv,ovoo_tile_32,&
                           & ccsdpt_doubles_2_1,wrk_3d,handle,cublas_handle)

       else if (idx .eq. 5) then ! ikj

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex3,oindex1,oindex2,nv,no,vvoo_tile_12,vvoo_tile_21,&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex3,oindex1,oindex2,nv,no,vvoo_tile_13,vvoo_tile_31,&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex2,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles_23,&
                           & wrk_3d,trip,vvvo_tile_o1,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex2,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles_13,&
                           & wrk_3d,trip,vvvo_tile_o2,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex1,oindex2,oindex3,no,nv,ovoo_tile_23,ovoo_tile_32,&
                           & ccsdpt_doubles_2_1,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex1,oindex2,oindex3,no,nv,ovoo_tile_12,&
                           & ccsdpt_doubles_2_3,trip,handle,cublas_handle)

       else if (idx .eq. 6) then ! jik

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_ijk_11(oindex3,oindex2,oindex1,nv,no,vvoo_tile_21,vvoo_tile_12,&
                       & ccsdpt_singles_3,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_12(oindex3,oindex2,oindex1,nv,no,vvoo_tile_23,vvoo_tile_32,&
                       & ccsdpt_singles_1,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_ijk_211(oindex1,oindex3,oindex2,nv,no,&
                           & ccsdpt_doubles_13,&
                           & trip,wrk_3d,vvvo_tile_o2,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_212(oindex1,oindex3,oindex2,nv,no,&
                           & ccsdpt_doubles_23,&
                           & trip,wrk_3d,vvvo_tile_o1,handle,cublas_handle)

          ! now do occ part:

          call ccsdpt_contract_ijk_221(oindex2,oindex1,oindex3,no,nv,ovoo_tile_13,ovoo_tile_31,&
                           & ccsdpt_doubles_2_2,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_ijk_222(oindex2,oindex1,oindex3,no,nv,ovoo_tile_21,&
                           & ccsdpt_doubles_2_3,wrk_3d,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_ijk_case3


  !> \brief: driver routine for contractions in case(3) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine ccsdpt_driver_abc_case3(vindex1,vindex2,vindex3,no,nv,&
                            & oovv_tile_12, oovv_tile_13, oovv_tile_21,&
                            & oovv_tile_23, oovv_tile_31, oovv_tile_32,vovv,&
                            & ooov_tile_v1,ooov_tile_v2,ooov_tile_v3,&
                            & ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3,&
                            & ccsdpt_doubles_12,ccsdpt_doubles_13,ccsdpt_doubles_21,&
                            & ccsdpt_doubles_23,ccsdpt_doubles_31,ccsdpt_doubles_32,&
                            & ccsdpt_doubles_2_1,ccsdpt_doubles_2_2,ccsdpt_doubles_2_3,&
                            & wrk_3d,trip,async_idx,num_idxs,cublas_handle,&
                            & a,b,c,tile_size_a,tile_size_b,tile_size_c,&
                            & vovv_tile_1, vovv_tile_2, vovv_tile_3)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: vindex1,vindex2,vindex3,no,nv
    integer, intent(in) :: a,b,c,tile_size_a,tile_size_b,tile_size_c
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(no,no) :: ccsdpt_doubles_12,ccsdpt_doubles_13,ccsdpt_doubles_21
    real(realk), dimension(no,no) :: ccsdpt_doubles_23,ccsdpt_doubles_31,ccsdpt_doubles_32
    real(realk), dimension(nv,no,no) :: ccsdpt_doubles_2_1,ccsdpt_doubles_2_2,ccsdpt_doubles_2_3
    real(realk), dimension(no) :: ccsdpt_singles_1,ccsdpt_singles_2,ccsdpt_singles_3
    !> vovv integrals
    type(tensor), intent(inout)  :: vovv
    !> tiles of vovv integrals
    real(realk), dimension(nv,no,nv,tile_size_a), intent(inout), optional :: vovv_tile_1
    real(realk), dimension(nv,no,nv,tile_size_b), intent(inout), optional :: vovv_tile_2
    real(realk), dimension(nv,no,nv,tile_size_c), intent(inout), optional :: vovv_tile_3
    !> tiles of oovv integrals
    real(realk), dimension(no,no) :: oovv_tile_12, oovv_tile_13, oovv_tile_21
    real(realk), dimension(no,no) :: oovv_tile_23, oovv_tile_31, oovv_tile_32
    !> tiles of ooov 2-el integrals determined by incomming vindex1,vindex2,vindex3
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v1
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v2
    real(realk), dimension(no,no,no), intent(inout) :: ooov_tile_v3
    !> triples amplitude and work array
    real(realk), dimension(no,no,no) :: trip, wrk_3d
    !> loop integer
    integer :: idx
    integer, intent(in) :: num_idxs
    integer*4 :: stat
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx(num_idxs), handle
#ifdef VAR_PGF90
    integer*4, external :: acc_set_cuda_stream
#endif
#else
    integer :: async_idx(num_idxs), handle
#endif
    type(c_ptr) :: cublas_handle

    handle = async_idx(5)

#ifdef VAR_CUBLAS
    stat = acc_set_cuda_stream(handle,cublas_handle)
#endif

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. 
    ! in order to minimize the number of reorderings needed to be performed, 
    ! in order to take optimal advantage of the symmetry of the amplitudes, 
    ! AND to finish all work that involves int_occ_tile_v3 first such that this may be updated (gpu),
    ! we carry out the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : abc --213--> bac --132--> bca --321--> acb --213--> cab --132--> cba
    ! in 211/212 : cab ........ acb ........ abc ........ cba ........ bca ........ bac
    ! in 221/222 : bca ........ cba ........ cab ........ bac ........ abc ........ acb

    do idx = 1,6

       if (idx .eq. 1) then ! abc

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex1,vindex2,vindex3,nv,no,oovv_tile_23,oovv_tile_32,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex2,vindex3,nv,no,oovv_tile_21,oovv_tile_12,&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex2,vindex3,vindex1,no,nv,vovv_tile_1(:,:,vindex3,a),vovv_tile_3(:,:,vindex1,c),&
                           & ccsdpt_doubles_2_2,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex3,vindex1,no,nv,vovv_tile_3(:,:,vindex2,c),&
                           & ccsdpt_doubles_2_1,trip,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex2,vindex3,vindex1,no,nv,vovv%elm4(:,:,vindex3,vindex1),vovv%elm4(:,:,vindex1,vindex3),&
                           & ccsdpt_doubles_2_2,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex3,vindex1,no,nv,vovv%elm4(:,:,vindex2,vindex3),&
                           & ccsdpt_doubles_2_1,trip,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex3,vindex1,vindex2,nv,no,&
                           & ccsdpt_doubles_31,&
                           & wrk_3d,trip,ooov_tile_v2,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex3,vindex1,vindex2,nv,no,&
                           & ccsdpt_doubles_21,&
                           & wrk_3d,trip,ooov_tile_v3,handle,cublas_handle)

       else if (idx .eq. 2) then ! cab

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[2,1,3],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[2,1,3],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex2,vindex1,vindex3,nv,no,oovv_tile_13,oovv_tile_31,&
                       & ccsdpt_singles_2,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex2,vindex1,vindex3,nv,no,oovv_tile_12,oovv_tile_21,&
                       & ccsdpt_singles_3,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex1,vindex3,vindex2,no,nv,vovv_tile_2(:,:,vindex3,b),vovv_tile_3(:,:,vindex2,c),&
                           & ccsdpt_doubles_2_1,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex3,vindex2,no,nv,vovv_tile_3(:,:,vindex1,c),&
                           & ccsdpt_doubles_2_2,wrk_3d,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex1,vindex3,vindex2,no,nv,vovv%elm4(:,:,vindex3,vindex2),vovv%elm4(:,:,vindex2,vindex3),&
                           & ccsdpt_doubles_2_1,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex3,vindex2,no,nv,vovv%elm4(:,:,vindex1,vindex3),&
                           & ccsdpt_doubles_2_2,wrk_3d,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex3,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_32,&
                           & trip,wrk_3d,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex3,vindex2,vindex1,nv,no,&
                           & ccsdpt_doubles_12,&
                           & trip,wrk_3d,ooov_tile_v3,handle,cublas_handle)

       else if (idx .eq. 3) then ! bca

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[1,3,2],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[1,3,2],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex2,vindex3,vindex1,nv,no,oovv_tile_31,oovv_tile_13,&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex2,vindex3,vindex1,nv,no,oovv_tile_32,oovv_tile_23,&
                       & ccsdpt_singles_1,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex3,vindex1,vindex2,no,nv,vovv_tile_2(:,:,vindex1,b),vovv_tile_1(:,:,vindex2,a),&
                           & ccsdpt_doubles_2_3,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex1,vindex2,no,nv,vovv_tile_1(:,:,vindex3,a),&
                           & ccsdpt_doubles_2_2,trip,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex3,vindex1,vindex2,no,nv,vovv%elm4(:,:,vindex1,vindex2),vovv%elm4(:,:,vindex2,vindex1),&
                           & ccsdpt_doubles_2_3,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex1,vindex2,no,nv,vovv%elm4(:,:,vindex3,vindex1),&
                           & ccsdpt_doubles_2_2,trip,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex2,vindex3,nv,no,&
                           & ccsdpt_doubles_12,&
                           & wrk_3d,trip,ooov_tile_v3,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex2,vindex3,nv,no,&
                           & ccsdpt_doubles_32,&
                           & wrk_3d,trip,ooov_tile_v1,handle,cublas_handle)

       else if (idx .eq. 4) then ! cba

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[3,2,1],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[3,2,1],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex1,vindex3,vindex2,nv,no,oovv_tile_32,oovv_tile_23,&
                       & ccsdpt_singles_1,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex1,vindex3,vindex2,nv,no,oovv_tile_31,oovv_tile_13,&
                       & ccsdpt_singles_2,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex3,vindex2,vindex1,no,nv,vovv_tile_1(:,:,vindex2,a),vovv_tile_2(:,:,vindex1,b),&
                           & ccsdpt_doubles_2_3,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex2,vindex1,no,nv,vovv_tile_2(:,:,vindex3,b),&
                           & ccsdpt_doubles_2_1,wrk_3d,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex3,vindex2,vindex1,no,nv,vovv%elm4(:,:,vindex2,vindex1),vovv%elm4(:,:,vindex1,vindex2),&
                           & ccsdpt_doubles_2_3,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex3,vindex2,vindex1,no,nv,vovv%elm4(:,:,vindex3,vindex2),&
                           & ccsdpt_doubles_2_1,wrk_3d,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_21,&
                           & trip,wrk_3d,ooov_tile_v3,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex2,vindex1,vindex3,nv,no,&
                           & ccsdpt_doubles_31,&
                           & trip,wrk_3d,ooov_tile_v2,handle,cublas_handle)

       else if (idx .eq. 5) then ! acb

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,1,3],0.0E0_realk,trip,handle)
#else
          call array_reorder_3d(1.0E0_realk,wrk_3d,no,no,&
                           & no,[2,1,3],0.0E0_realk,trip)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex3,vindex1,vindex2,nv,no,oovv_tile_12,oovv_tile_21,&
                       & ccsdpt_singles_3,&
                       & trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex3,vindex1,vindex2,nv,no,oovv_tile_13,oovv_tile_31,&
                       & ccsdpt_singles_2,&
                       & trip,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex1,vindex2,vindex3,no,nv,vovv_tile_3(:,:,vindex2,c),vovv_tile_2(:,:,vindex3,b),&
                           & ccsdpt_doubles_2_1,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex2,vindex3,no,nv,vovv_tile_2(:,:,vindex1,b),&
                           & ccsdpt_doubles_2_3,trip,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex1,vindex2,vindex3,no,nv,vovv%elm4(:,:,vindex2,vindex3),vovv%elm4(:,:,vindex3,vindex2),&
                           & ccsdpt_doubles_2_1,trip,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex1,vindex2,vindex3,no,nv,vovv%elm4(:,:,vindex1,vindex2),&
                           & ccsdpt_doubles_2_3,trip,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex2,vindex3,vindex1,nv,no,&
                           & ccsdpt_doubles_23,&
                           & wrk_3d,trip,ooov_tile_v1,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex2,vindex3,vindex1,nv,no,&
                           & ccsdpt_doubles_13,&
                           & wrk_3d,trip,ooov_tile_v2,handle,cublas_handle)

       else if (idx .eq. 6) then ! bac

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array

#ifdef VAR_OPENACC
          call array_reorder_3d_acc(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d,handle)
#else
          call array_reorder_3d(1.0E0_realk,trip,no,no,&
                           & no,[1,3,2],0.0E0_realk,wrk_3d)
#endif

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_abc_11(vindex3,vindex2,vindex1,nv,no,oovv_tile_21,oovv_tile_12,&
                       & ccsdpt_singles_3,&
                       & wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_12(vindex3,vindex2,vindex1,nv,no,oovv_tile_23,oovv_tile_32,&
                       & ccsdpt_singles_1,&
                       & wrk_3d,.false.,handle,cublas_handle)

          ! calculate contributions to ccsdpt_doubles (virt part):

#ifdef VAR_MPI
          call ccsdpt_contract_abc_211(vindex2,vindex1,vindex3,no,nv,vovv_tile_3(:,:,vindex1,c),vovv_tile_1(:,:,vindex3,a),&
                           & ccsdpt_doubles_2_2,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex1,vindex3,no,nv,vovv_tile_1(:,:,vindex2,a),&
                           & ccsdpt_doubles_2_3,wrk_3d,handle,cublas_handle)
#else
          call ccsdpt_contract_abc_211(vindex2,vindex1,vindex3,no,nv,vovv%elm4(:,:,vindex1,vindex3),vovv%elm4(:,:,vindex3,vindex1),&
                           & ccsdpt_doubles_2_2,wrk_3d,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_212(vindex2,vindex1,vindex3,no,nv,vovv%elm4(:,:,vindex2,vindex1),&
                           & ccsdpt_doubles_2_3,wrk_3d,handle,cublas_handle)
#endif

          ! now do occ part:

          call ccsdpt_contract_abc_221(vindex1,vindex3,vindex2,nv,no,&
                           & ccsdpt_doubles_13,&
                           & trip,wrk_3d,ooov_tile_v2,.false.,handle,cublas_handle)
          call ccsdpt_contract_abc_222(vindex1,vindex3,vindex2,nv,no,&
                           & ccsdpt_doubles_23,&
                           & trip,wrk_3d,ooov_tile_v1,handle,cublas_handle)

       end if

    end do

  end subroutine ccsdpt_driver_abc_case3


  !> \brief: transform ccsd doubles from local to canonical basis
  !> \author: Janus Juul Eriksen
  !> \date: september 2012
  !> \param: ccsd_t2, no and nv are nocc and nvirt, respectively, and U_occ and U_virt
  !          are unitary matrices from local --> canonical basis
  subroutine ccsdpt_local_can_trans(ccsd_t2_arr,no,nv,U_occ,U_virt)

    implicit none
    !> ccsd doubles
    type(tensor), intent(inout) :: ccsd_t2_arr
    !> unitary transformation matrices
    type(array2), intent(inout) :: U_occ, U_virt
    !> integers
    integer, intent(in) :: no, nv
    !> temp array4 structures
    type(array4) :: tmp1, tmp2
    type(array4) :: ccsd_t2

    !FIXME: THIS IS A CRAPPY HACK, USING A LOT OF MEMORY
    ccsd_t2       = array4_init(ccsd_t2_arr%dims)
    ccsd_t2%val   = ccsd_t2_arr%elm4
    call deassoc_ptr_arr(ccsd_t2_arr) 

    ! (a,i,b,j) are local basis indices and (A,I,B,J) refer to the canonical basis.
    ! we want to carry out the transformation:
    ! T^{AB}_{IJ} = sum_{aibj} U_{aA} U_{iI} U_{bB} U_{jJ} T^{ab}_{ij}

    ! 1. Init temporary arrays, dims = aiai
    tmp1 = array4_init_standard([nv,no,nv,no])

    ! 2. 1st index: doub_ampl(a,i,b,j) --> tmp1(A,i,b,j)
    call array4_contract1(ccsd_t2,U_virt,tmp1,.true.)
    call array4_free(ccsd_t2)

    ! 3. 2nd index: tmp1(A,i,b,j) --> tmp1(i,A,b,j) --> tmp2(I,A,b,j)
    call array4_reorder(tmp1,[2,1,3,4])
    tmp2 = array4_init_standard([no,nv,nv,no])
    call array4_contract1(tmp1,U_occ,tmp2,.true.)
    call array4_free(tmp1)

    ! 4. 3rd index: tmp2(I,A,b,j) --> tmp2(j,b,A,I) --> tmp1(J,b,A,I)
    call array4_reorder(tmp2,[4,3,2,1])
    tmp1 = array4_init_standard([no,nv,nv,no])
    call array4_contract1(tmp2,U_occ,tmp1,.true.)
    call array4_free(tmp2)

    ! 5. 4th index: tmp1(J,b,A,I) --> tmp1(b,J,A,I) --> doub_ampl(B,J,A,I) = doub_ampl(A,I,B,J)
    call array4_reorder(tmp1,[2,1,3,4])
    ccsd_t2 = array4_init_standard([nv,no,nv,no])
    call array4_contract1(tmp1,U_virt,ccsd_t2,.true.)
    call array4_free(tmp1)

    ccsd_t2_arr%dims = ccsd_t2%dims        
    call assoc_ptr_arr(ccsd_t2_arr) 
    ccsd_t2_arr%elm4 = ccsd_t2%val  
    call array4_free(ccsd_t2)

  end subroutine ccsdpt_local_can_trans


  !> \brief: transform ccsd_doubles, ccsdpt_singles and ccsdpt_doubles from canonical to local basis
  !> \author: Janus Juul Eriksen
  !> \date: september 2012
  !> \param: ccsd_t2, ccsdpt_t1, ccsdpt_t2, no and nv are nocc and nvirt, respectively, 
  !<         and U_occ and U_virt are unitary matrices from canonical --> local basis
  subroutine ccsdpt_can_local_trans(ccsd_t2_arr,ccsdpt_t1_arr,ccsdpt_t2_arr,no,nv,U_occ,U_virt)

    implicit none
    !> ccsdpt_singles
    type(tensor), intent(inout) :: ccsdpt_t1_arr
    !> ccsd_doubles and ccsdpt_doubles
    type(tensor), intent(inout) :: ccsd_t2_arr, ccsdpt_t2_arr
    !> unitary transformation matrices
    type(array2), intent(inout) :: U_occ, U_virt
    !> integers
    integer, intent(in) :: no, nv
    !> temp array2 and array4 structures
    type(array2) :: tmp0
    type(array4) :: tmp1, tmp2, tmp3, tmp4
    type(array2) :: ccsdpt_t1
    type(array4) :: ccsd_t2, ccsdpt_t2

    !FIXME: THIS IS A CRAPPY HACK, USING A LOT OF MEMORY
    ccsdpt_t1     = array2_init(ccsdpt_t1_arr%dims)
    ccsdpt_t1%val = ccsdpt_t1_arr%elm2
    ccsdpt_t2     = array4_init(ccsdpt_t2_arr%dims)
    ccsdpt_t2%val = ccsdpt_t2_arr%elm4
    ccsd_t2       = array4_init(ccsd_t2_arr%dims)
    ccsd_t2%val   = ccsd_t2_arr%elm4
    
    call deassoc_ptr_arr(ccsdpt_t1_arr) 
    call deassoc_ptr_arr(ccsdpt_t2_arr) 
    call deassoc_ptr_arr(ccsd_t2_arr) 

    ! (a,i,b,j) are local basis indices and (A,I,B,J) refer to the canonical basis.

    ! 1a. init temporary array4s, tmp1 and tmp3
    tmp1 = array4_init_standard([nv,nv,no,no])
    tmp3 = array4_init_standard([nv,nv,no,no])
    ! 1b. init temporary array2, tmp0
    tmp0 = array2_init_plain([nv,no])

    ! 2. 1st index:
    ! ccsdpt_t2(A,B,I,J) --> tmp1(a,B,I,J)
    ! ccsd_t2(B,A,J,I) --> tmp3(b,A,J,I)
    call array4_contract1(ccsdpt_t2,U_virt,tmp1,.true.)
    call array4_contract1(ccsd_t2,U_virt,tmp3,.true.)
    ! ccsdpt_t1(A,I) --> tmp0(a,I)
    call array2_matmul(U_virt,ccsdpt_t1,tmp0,'t','n',1.0E0_realk,0.0E0_realk)

    ! free ccsdpt_doubles and ccsd_doubles
    call array4_free(ccsdpt_t2)
    call array4_free(ccsd_t2)

    ! 3. 2nd index:
    ! tmp1(a,B,I,J) --> tmp1(B,a,I,J) --> tmp2(b,a,I,J)
    ! tmp3(b,A,J,I) --> tmp3(A,b,J,I) --> tmp4(a,b,J,I) 
    ! tmp0(a,I) --> ccsdpt_t1(a,i)
    call array4_reorder(tmp1,[2,1,3,4])
    call array4_reorder(tmp3,[2,1,3,4])

    ! init temporary array4s, tmp2 and tmp4
    tmp2 = array4_init_standard([nv,nv,no,no])
    tmp4 = array4_init_standard([nv,nv,no,no])

    ! transformation time - ccsdpt_doubles and ccsd_doubles case
    call array4_contract1(tmp1,U_virt,tmp2,.true.)
    call array4_contract1(tmp3,U_virt,tmp4,.true.)
    ! ccsdpt_singles case
    call array2_matmul(tmp0,U_occ,ccsdpt_t1,'n','n',1.0E0_realk,0.0E0_realk)

    ! free tmp1 and tmp3
    call array4_free(tmp1)
    call array4_free(tmp3)
    ! free tmp0
    call array2_free(tmp0)

    ! 4. 3rd index:
    ! tmp2(b,a,I,J) --> tmp2(J,I,a,b) --> tmp1(j,I,a,b)
    ! tmp4(a,b,J,I) --> tmp4(I,J,b,a) --> tmp3(i,J,b,a)
    call array4_reorder(tmp2,[4,3,2,1])
    call array4_reorder(tmp4,[4,3,2,1])

    ! init temporary array4s, tmp1 and tmp3, once again
    tmp1 = array4_init_standard([no,no,nv,nv])
    tmp3 = array4_init_standard([no,no,nv,nv])

    ! transformation time
    call array4_contract1(tmp2,U_occ,tmp1,.true.)
    call array4_contract1(tmp4,U_occ,tmp3,.true.)

    ! free tmp2 and tmp4
    call array4_free(tmp2)
    call array4_free(tmp4)

    ! 5. 4th index:
    ! tmp1(j,I,a,b) --> tmp1(I,j,a,b) --> ccsdpt_doubles(i,j,a,b)
    ! tmp3(i,J,b,a) --> tmp3(J,i,b,a) --> ccsd_doubles(j,i,b,a)
    call array4_reorder(tmp1,[2,1,3,4])
    call array4_reorder(tmp3,[2,1,3,4])

    ! init ccsdpt_t2 and ccsd_t2 array4s once again
    ccsdpt_t2 = array4_init_standard([no,no,nv,nv])
    ccsd_t2 = array4_init_standard([no,no,nv,nv])

    ! transformation time
    call array4_contract1(tmp1,U_occ,ccsdpt_t2,.true.)
    call array4_contract1(tmp3,U_occ,ccsd_t2,.true.)

    ! free tmp1 and tmp3
    call array4_free(tmp1)
    call array4_free(tmp3)

    ccsdpt_t1_arr%dims    =  ccsdpt_t1%dims      
    ccsdpt_t2_arr%dims    =  ccsdpt_t2%dims      
    ccsd_t2_arr%dims      =  ccsd_t2%dims        

    call assoc_ptr_arr(ccsdpt_t1_arr) 
    call assoc_ptr_arr(ccsdpt_t2_arr) 
    call assoc_ptr_arr(ccsd_t2_arr) 

    ccsdpt_t1_arr%elm2 = ccsdpt_t1%val
    ccsdpt_t2_arr%elm4 = ccsdpt_t2%val
    ccsd_t2_arr%elm4   = ccsd_t2%val  

    call array2_free(ccsdpt_t1)
    call array4_free(ccsdpt_t2)
    call array4_free(ccsd_t2)

  end subroutine ccsdpt_can_local_trans


  !> \brief: create VIRTUAL part of a triples amplitude ([a,b,c] tuple) for a fixed [i,j,k] tuple, that is, t^{***}_{ijk}
  !          saved as an array3 structure (amplitudes)
  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  !> \param: oindex1, oindex2, and oindex3 are the three occupied indices of the outer loop in the ccsd(t) driver
  !> \param: no and nv are nocc and nvirt, respectively
  !> \param: doub_ampl are ccsd ampltidues, t^{ab}_{ij}
  !> \param: int_virt is a v^3 part of cbai of driver routine
  !> \param: trip holds the triples tuple [a,b,c], that is, of the size (virt)³ kept in memory
  subroutine trip_amplitudes_ijk_virt(oindex1,oindex2,oindex3,no,nv,doub_ampl_v2,int_virt_tile,trip,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), dimension(nv,nv,nv), target, intent(in) :: int_virt_tile
    real(realk), dimension(nv,nv), target, intent(in) :: doub_ampl_v2
    real(realk), dimension(nv,nv,nv), target, intent(inout) :: trip
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(doub_ampl_v2,int_virt_tile,trip)
    call dgemm_acc_openacc_async(async_idx,'t','n',nv,nv**2,nv,1.0E0_realk,doub_ampl_v2,nv,int_virt_tile,nv,&
!    call dgemm_acc('t','n',nv,nv2,nv,1.0E0_realk,doub_ampl_v2,nv,int_virt_tile,nv,&
                      & 0.0E0_realk,trip,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(doub_ampl_v2,int_virt_tile,trip)
    stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(nv,kind=4),int(nv**2,kind=4),int(nv,kind=4),&
                          & 1.0E0_realk,c_loc(doub_ampl_v2),int(nv,kind=4),c_loc(int_virt_tile),int(nv,kind=4),&
                          & 0.0E0_realk,c_loc(trip),int(nv,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (trip_amplitudes_ijk_virt) = ",stat
!       stop
!    end if

#endif
#else
    call dgemm('t','n',nv,nv**2,nv,1.0E0_realk,doub_ampl_v2,nv,int_virt_tile,nv,&
                   & 0.0E0_realk,trip,nv)
#endif

  end subroutine trip_amplitudes_ijk_virt


  !> \brief: create OCCUPIED part of a triples amplitude ([i,j,k] tuple) for a fixed [a,b,c] tuple, that is, t^{***}_{abc}
  !          saved as an array3 structure (amplitudes)
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine trip_amplitudes_abc_occ(vindex1,vindex2,vindex3,no,nv,doub_ampl_o2,int_occ_tile,trip,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, no, nv
    real(realk), dimension(no,no,no), target, intent(in) :: int_occ_tile
    real(realk), dimension(no,no), target, intent(in) :: doub_ampl_o2
    real(realk), dimension(no,no,no), target, intent(inout) :: trip
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(doub_ampl_o2,int_occ_tile,trip)
    call dgemm_acc_openacc_async(async_idx,'t','n',no,no**2,no,-1.0E0_realk,doub_ampl_o2,no,int_occ_tile,no,&
                   & 0.0E0_realk,trip,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(doub_ampl_o2,int_occ_tile,trip)
    stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(no,kind=4),int(no**2,kind=4),int(no,kind=4),&
                          & -1.0E0_realk,c_loc(doub_ampl_o2),int(no,kind=4),c_loc(int_occ_tile),int(no,kind=4),&
                          & 0.0E0_realk,c_loc(trip),int(no,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (trip_amplitudes_abc_occ) = ",stat
!       stop
!    end if

#endif
#else
    call dgemm('t','n',no,no**2,no,-1.0E0_realk,doub_ampl_o2,no,int_occ_tile,no,&
                   & 0.0E0_realk,trip,no)
#endif

  end subroutine trip_amplitudes_abc_occ


  !> \brief: create OCCUPIED part of a triples amplitude ([a,b,c] tuple) for a fixed [i,j,k] tuple, that is, t^{***}_{ijk}
  !          saved as an array3 structure (amplitudes)
  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  !> \param: oindex1, oindex2, and oindex3 are the three occupied indices of the outer loop in the ccsd(t) driver
  !> \param: no and nv are nocc and nvirt, respectively
  !> \param: doub_ampl are ccsd ampltidues, t^{ab}_{ij}
  !> \param: int_occ is a ov part of jaik of driver routine
  !> \param: trip holds the triples tuple [c,a,b], that is, of the size (virt)³ kept in memory
  subroutine trip_amplitudes_ijk_occ(oindex1,oindex2,oindex3,no,nv,doub_ampl_ov2,int_occ_portion,trip,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), dimension(no,nv), target, intent(in) :: int_occ_portion
    real(realk), dimension(no,nv,nv), target, intent(in) :: doub_ampl_ov2
    real(realk), dimension(nv,nv,nv), target, intent(inout) :: trip
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(int_occ_portion,doub_ampl_ov2,trip)
    call dgemm_acc_openacc_async(async_idx,'t','n',nv,nv**2,no,-1.0E0_realk,int_occ_portion,no,doub_ampl_ov2,no,&
!    call dgemm_acc('t','n',nv,nv2,no,-1.0E0_realk,int_occ_portion,no,doub_ampl_ov2,no,&
                      & 1.0E0_realk,trip,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(int_occ_portion,doub_ampl_ov2,trip)
    stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(nv,kind=4),int(nv**2,kind=4),int(no,kind=4),&
                          & -1.0E0_realk,c_loc(int_occ_portion),int(no,kind=4),c_loc(doub_ampl_ov2),int(no,kind=4),&
                          & 1.0E0_realk,c_loc(trip),int(nv,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (trip_amplitudes_ijk_occ) = ",stat
!       stop
!    end if

#endif
#else
    call dgemm('t','n',nv,nv**2,no,-1.0E0_realk,int_occ_portion,no,doub_ampl_ov2,no,&
                   & 1.0E0_realk,trip,nv)
#endif

  end subroutine trip_amplitudes_ijk_occ


  !> \brief: create VIRTUAL part of a triples amplitude ([i,j,k] tuple) for a fixed [a,b,c] tuple, that is, t^{***}_{abc}
  !          saved as an array3 structure (amplitudes)
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  subroutine trip_amplitudes_abc_virt(vindex1,vindex2,vindex3,no,nv,doub_ampl_vo2,int_virt_portion,trip,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, no, nv
    real(realk), dimension(nv,no), target, intent(in) :: int_virt_portion
    real(realk), dimension(nv,no,no), target, intent(in) :: doub_ampl_vo2
    real(realk), dimension(no,no,no), target, intent(inout) :: trip
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(int_virt_portion,doub_ampl_vo2,trip)
    call dgemm_acc_openacc_async(async_idx,'t','n',no,no**2,nv,1.0E0_realk,int_virt_portion,nv,doub_ampl_vo2,nv,&
                   & 1.0E0_realk,trip,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(int_virt_portion,doub_ampl_vo2,trip)
    stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(no,kind=4),int(no**2,kind=4),int(nv,kind=4),&
                          & 1.0E0_realk,c_loc(int_virt_portion),int(nv,kind=4),c_loc(doub_ampl_vo2),int(nv,kind=4),&
                          & 1.0E0_realk,c_loc(trip),int(no,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (trip_amplitudes_abc_virt) = ",stat
!       stop
!    end if

#endif
#else
    call dgemm('t','n',no,no**2,nv,1.0E0_realk,int_virt_portion,nv,doub_ampl_vo2,nv,&
                   & 1.0E0_realk,trip,no)
#endif

  end subroutine trip_amplitudes_abc_virt


  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  !> \param: oindex1, oindex2, and oindex3 are the three occupied indices of the outer loop in the ccsd(t) driver
  !> \param: no and nv are nocc and nvirt, respectively
  !> \param: eigenocc and eigenvirt are vectors containing occupied and virtual orbital energies, respectively
  !> \param: amplitudes are the final triples amplitude tuple [a,b,c], that is, of the size (virt)³ kept in memory
  subroutine trip_denom_ijk_cpu(oindex1,oindex2,oindex3,no,nv,eigenocc,eigenvirt,trip)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    real(realk), dimension(nv,nv,nv), intent(inout) :: trip
    !> temporary quantities
    integer :: a, b, c
    real(realk) :: e_orb_occ

    ! at first, calculate the sum of the three participating occupied orbital energies, as this
    ! is a constant for the three incomming occupied indices

    e_orb_occ = eigenocc(oindex1) + eigenocc(oindex2) + eigenocc(oindex3)

!$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,b,c),SHARED(nv,trip,eigenvirt,e_orb_occ)
    do a=1,nv
       do b=1,nv
          do c=1,nv

                  trip(c,b,a) = trip(c,b,a) / (e_orb_occ - eigenvirt(a) - eigenvirt(b) - eigenvirt(c))

          end do
       end do
    end do
!$OMP END PARALLEL DO

  end subroutine trip_denom_ijk_cpu


#ifdef VAR_OPENACC
  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  !> \param: oindex1, oindex2, and oindex3 are the three occupied indices of the outer loop in the ccsd(t) driver
  !> \param: no and nv are nocc and nvirt, respectively
  !> \param: eigenocc and eigenvirt are vectors containing occupied and virtual orbital energies, respectively
  !> \param: amplitudes are the final triples amplitude tuple [a,b,c], that is, of the size (virt)³ kept in memory
  subroutine trip_denom_ijk_acc(oindex1,oindex2,oindex3,no,nv,eigenocc,eigenvirt,trip,async_idx)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    real(realk), dimension(nv,nv,nv), intent(inout) :: trip
    !> temporary quantities
    integer :: a, b, c
    integer(kind=acc_handle_kind) :: async_idx
    real(realk) :: e_orb_occ

    ! at first, calculate the sum of the three participating occupied orbital energies, as this
    ! is a constant for the three incomming occupied indices

    e_orb_occ = eigenocc(oindex1) + eigenocc(oindex2) + eigenocc(oindex3)

!$acc parallel present(trip,eigenvirt) firstprivate(nv,e_orb_occ) &
!$acc& private(a,b,c) async(async_idx)
!$acc loop gang
    do a=1,nv
!$acc loop worker
       do b=1,nv
!$acc loop vector
          do c=1,nv

                  trip(c,b,a) = trip(c,b,a) / (e_orb_occ - eigenvirt(a) - eigenvirt(b) - eigenvirt(c))

          end do
!$acc end loop
       end do
!$acc end loop
    end do
!$acc end loop
!$acc end parallel

  end subroutine trip_denom_ijk_acc
#endif


  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  !> \param: vindex1, vindex2, and vindex3 are the three virtual indices of the outer loop in the ccsd(t) driver
  !> \param: no and nv are nocc and nvirt, respectively
  !> \param: eigenocc and eigenvirt are vectors containing occupied and virtual orbital energies, respectively
  !> \param: amplitudes are the final triples amplitude tuple [i,j,k], that is, of the size (occ)³ kept in memory
  subroutine trip_denom_abc_cpu(vindex1,vindex2,vindex3,no,nv,eigenocc,eigenvirt,trip)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, no, nv
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    real(realk), dimension(no,no,no), intent(inout) :: trip
    !> temporary quantities
    integer :: i, j, k
    real(realk) :: e_orb_virt

    ! at first, calculate the sum of the three participating virtual orbital energies, as this
    ! is a constant for the three incomming occupied indices

    e_orb_virt = eigenvirt(vindex1) + eigenvirt(vindex2) + eigenvirt(vindex3)

!$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,j,k),SHARED(no,trip,eigenocc,e_orb_virt)
    do i=1,no
       do j=1,no
          do k=1,no

                  trip(k,j,i) = trip(k,j,i) / (eigenocc(i) + eigenocc(j) + eigenocc(k) - e_orb_virt)

          end do
       end do
    end do
!$OMP END PARALLEL DO

  end subroutine trip_denom_abc_cpu


#ifdef VAR_OPENACC
  !> \author: Janus Juul Eriksen
  !> \date: april 2014
  !> \param: vindex1, vindex2, and vindex3 are the three virtual indices of the outer loop in the ccsd(t) driver
  !> \param: no and nv are nocc and nvirt, respectively
  !> \param: eigenocc and eigenvirt are vectors containing occupied and virtual orbital energies, respectively
  !> \param: amplitudes are the final triples amplitude tuple [i,j,k], that is, of the size (occ)³ kept in memory
  subroutine trip_denom_abc_acc(vindex1,vindex2,vindex3,no,nv,eigenocc,eigenvirt,trip,async_idx)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, no, nv
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    real(realk), dimension(no,no,no), intent(inout) :: trip
    !> temporary quantities
    integer :: i, j, k
    integer(kind=acc_handle_kind) :: async_idx
    real(realk) :: e_orb_virt

    ! at first, calculate the sum of the three participating virtual orbital energies, as this
    ! is a constant for the three incomming occupied indices

    e_orb_virt = eigenvirt(vindex1) + eigenvirt(vindex2) + eigenvirt(vindex3)

!$acc parallel present(trip,eigenocc) firstprivate(no,e_orb_virt)&
!$acc& private(i,j,k) async(async_idx)
!$acc loop gang
    do i=1,no
!$acc loop worker
       do j=1,no
!$acc loop vector
          do k=1,no

                  trip(k,j,i) = trip(k,j,i) / (eigenocc(i) + eigenocc(j) + eigenocc(k) - e_orb_virt)

          end do
!$acc end loop
       end do
!$acc end loop
    end do
!$acc end loop
!$acc end parallel

  end subroutine trip_denom_abc_acc
#endif


  !> brief: do the first of the two contraction over 'cdkl' (here: 'c' and 'd', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.63) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine. int_normal is abij of driver.
  !> nv is nvirt and T_star is ccsdpt_singles of driver. trip_ampl is the triples amplitude array.
  subroutine ccsdpt_contract_ijk_11(oindex1,oindex2,oindex3,nv,no,int_normal_23,int_normal_32,T_star_o1,&
                              & trip_ampl,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv), target :: T_star_o1 ! T_star(:,oinedx1)
    real(realk), dimension(nv,nv), target :: int_normal_23, int_normal_32
    real(realk), dimension(nv,nv,nv), target :: trip_ampl
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_ijk_11: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal and we thus only need to contract with the coulumb part. 

       ! now contract coulumb term over both indices
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,int_normal_23,T_star_o1)
       call dgemm_acc_openacc_async(async_idx,'n','n',nv,1,nv**2,&
!       call dgemm_acc('n','n',nv,1,nv2,&
                & 1.0E0_realk,trip_ampl,nv,int_normal_23,nv**2,1.0E0_realk,T_star_o1,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,int_normal_23,T_star_o1)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(nv,kind=4),int(1,kind=4),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(trip_ampl),int(nv,kind=4),c_loc(int_normal_23),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o1),int(nv,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_ijk_11 - case 0) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',nv,1,nv**2,&
                & 1.0E0_realk,trip_ampl,nv,int_normal_23,nv**2,1.0E0_realk,T_star_o1,nv)
#endif

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,int_normal_23,int_normal_32,T_star_o1)
       call dgemm_acc_openacc_async(async_idx,'n','n',nv,1,nv**2,&
!       call dgemm_acc('n','n',nv,1,nv2,&
                & 2.0E0_realk,trip_ampl,nv,int_normal_23,nv**2,1.0E0_realk,T_star_o1,nv)
       call dgemm_acc_openacc_async(async_idx,'n','n',nv,1,nv**2,&
!       call dgemm_acc('n','n',nv,1,nv2,&
                & -1.0E0_realk,trip_ampl,nv,int_normal_32,nv**2,1.0E0_realk,T_star_o1,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,int_normal_23,int_normal_32,T_star_o1)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(nv,kind=4),int(1,kind=4),int(nv**2,kind=4),&
                             & 2.0E0_realk,c_loc(trip_ampl),int(nv,kind=4),c_loc(int_normal_23),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o1),int(nv,kind=4))
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(nv,kind=4),int(1,kind=4),int(nv**2,kind=4),&
                             & -1.0E0_realk,c_loc(trip_ampl),int(nv,kind=4),c_loc(int_normal_32),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o1),int(nv,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_ijk_11 - case 1) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',nv,1,nv**2,&
                & 2.0E0_realk,trip_ampl,nv,int_normal_23,nv**2,1.0E0_realk,T_star_o1,nv)
       call dgemm('n','n',nv,1,nv**2,&
                & -1.0E0_realk,trip_ampl,nv,int_normal_32,nv**2,1.0E0_realk,T_star_o1,nv)
#endif


    end select TypeofContraction_ijk_11

  end subroutine ccsdpt_contract_ijk_11


  !> brief: do the first of the two contraction over 'cdkl' (here: 'c' and 'd', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.63) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_11(vindex1,vindex2,vindex3,nv,no,int_normal_23,int_normal_32,T_star_v1,&
                              & trip_ampl,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, nv, no
    real(realk), dimension(no), target :: T_star_v1 ! T_star(:,oinedx1)
    real(realk), dimension(no,no), target :: int_normal_23, int_normal_32
    real(realk), dimension(no,no,no), target :: trip_ampl
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_abc_11: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal and we thus only need to contract with the coulumb part. 
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,int_normal_23,T_star_v1)
       call dgemm_acc_openacc_async(async_idx,'n','n',no,1,no**2,&
                & 1.0E0_realk,trip_ampl,no,int_normal_23,no**2,1.0E0_realk,T_star_v1,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,int_normal_23,T_star_v1)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(no,kind=4),int(1,kind=4),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(trip_ampl),int(no,kind=4),c_loc(int_normal_23),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v1),int(no,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_abc_11 (case 0)) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',no,1,no**2,&
                & 1.0E0_realk,trip_ampl,no,int_normal_23,no**2,1.0E0_realk,T_star_v1,no)
#endif

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,int_normal_23,int_normal_32,T_star_v1)
       call dgemm_acc_openacc_async(async_idx,'n','n',no,1,no**2,&
                & 2.0E0_realk,trip_ampl,no,int_normal_23,no**2,1.0E0_realk,T_star_v1,no)
       call dgemm_acc_openacc_async(async_idx,'n','n',no,1,no**2,&
                & -1.0E0_realk,trip_ampl,no,int_normal_32,no**2,1.0E0_realk,T_star_v1,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,int_normal_23,int_normal_32,T_star_v1)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(no,kind=4),int(1,kind=4),int(no**2,kind=4),&
                             & 2.0E0_realk,c_loc(trip_ampl),int(no,kind=4),c_loc(int_normal_23),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v1),int(no,kind=4))
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(no,kind=4),int(1,kind=4),int(no**2,kind=4),&
                             & -1.0E0_realk,c_loc(trip_ampl),int(no,kind=4),c_loc(int_normal_32),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v1),int(no,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_abc_11 (case 1)) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',no,1,no**2,&
                & 2.0E0_realk,trip_ampl,no,int_normal_23,no**2,1.0E0_realk,T_star_v1,no)
       call dgemm('n','n',no,1,no**2,&
                & -1.0E0_realk,trip_ampl,no,int_normal_32,no**2,1.0E0_realk,T_star_v1,no)
#endif

    end select TypeofContraction_abc_11

  end subroutine ccsdpt_contract_abc_11


  !> brief: do the second of the two contraction over 'cdkl' (here: 'c' and 'd', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.63) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine. int_normal is abij of driver.
  !> nv is nvirt and T_star is ccsdpt_singles of driver. trip_ampl is the triples amplitude array.
  subroutine ccsdpt_contract_ijk_12(oindex1,oindex2,oindex3,nv,no,int_normal_21,int_normal_12,T_star_o3,&
                              & trip_ampl,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv), target :: T_star_o3 ! T_star(:,oinedx3)
    real(realk), dimension(nv,nv), target :: int_normal_21, int_normal_12
    real(realk), dimension(nv,nv,nv), target :: trip_ampl
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_ijk_12: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal
       ! and we thus only need to contract with (-1)*coulumb part. 

       ! now contract coulumb term over both indices
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,int_normal_21,T_star_o3)
       call dgemm_acc_openacc_async(async_idx,'n','n',nv,1,nv**2,&
!       call dgemm_acc('n','n',nv,1,nv2,&
             & -1.0E0_realk,trip_ampl,nv,int_normal_21,nv**2,1.0E0_realk,T_star_o3,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,int_normal_21,T_star_o3)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(nv,kind=4),int(1,kind=4),int(nv**2,kind=4),&
                             & -1.0E0_realk,c_loc(trip_ampl),int(nv,kind=4),c_loc(int_normal_21),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o3),int(nv,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_ijk_12 - case 0) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',nv,1,nv**2,&
                & -1.0E0_realk,trip_ampl,nv,int_normal_21,nv**2,1.0E0_realk,T_star_o3,nv)
#endif

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,int_normal_21,int_normal_12,T_star_o3)
       call dgemm_acc_openacc_async(async_idx,'n','n',nv,1,nv**2,&
!       call dgemm_acc('n','n',nv,1,nv2,&
             & -2.0E0_realk,trip_ampl,nv,int_normal_21,nv**2,1.0E0_realk,T_star_o3,nv)
       call dgemm_acc_openacc_async(async_idx,'n','n',nv,1,nv**2,&
!       call dgemm_acc('n','n',nv,1,nv2,&
             & 1.0E0_realk,trip_ampl,nv,int_normal_12,nv**2,1.0E0_realk,T_star_o3,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,int_normal_21,int_normal_12,T_star_o3)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(nv,kind=4),int(1,kind=4),int(nv**2,kind=4),&
                             & -2.0E0_realk,c_loc(trip_ampl),int(nv,kind=4),c_loc(int_normal_21),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o3),int(nv,kind=4))
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(nv,kind=4),int(1,kind=4),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(trip_ampl),int(nv,kind=4),c_loc(int_normal_12),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o3),int(nv,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_ijk_12 - case 1) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',nv,1,nv**2,&
                & -2.0E0_realk,trip_ampl,nv,int_normal_21,nv**2,1.0E0_realk,T_star_o3,nv)
       call dgemm('n','n',nv,1,nv**2,&
                & 1.0E0_realk,trip_ampl,nv,int_normal_12,nv**2,1.0E0_realk,T_star_o3,nv)
#endif

    end select TypeofContraction_ijk_12

  end subroutine ccsdpt_contract_ijk_12


  !> brief: do the second of the two contraction over 'cdkl' (here: 'c' and 'd', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.63) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_12(vindex1,vindex2,vindex3,nv,no,int_normal_21,int_normal_12,T_star_v3,&
                              & trip_ampl,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, nv, no
    real(realk), dimension(no), target :: T_star_v3 ! T_star(:,oinedx3)
    real(realk), dimension(no,no), target :: int_normal_21, int_normal_12
    real(realk), dimension(no,no,no), target :: trip_ampl
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_abc_12: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal
       ! and we thus only need to contract with (-1)*coulumb part. 

       ! now contract coulumb term over both indices
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,int_normal_21,T_star_v3)
       call dgemm_acc_openacc_async(async_idx,'n','n',no,1,no**2,&
                & -1.0E0_realk,trip_ampl,no,int_normal_21,no**2,1.0E0_realk,T_star_v3,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,int_normal_21,T_star_v3)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(no,kind=4),int(1,kind=4),int(no**2,kind=4),&
                             & -1.0E0_realk,c_loc(trip_ampl),int(no,kind=4),c_loc(int_normal_21),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v3),int(no,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_abc_12 (case 0)) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',no,1,no**2,&
                & -1.0E0_realk,trip_ampl,no,int_normal_21,no**2,1.0E0_realk,T_star_v3,no)
#endif

    case(1)

       ! now contract coulumb term over both indices, then contract exchange term over both indices
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,int_normal_21,int_normal_12,T_star_v3)
       call dgemm_acc_openacc_async(async_idx,'n','n',no,1,no**2,&
                & -2.0E0_realk,trip_ampl,no,int_normal_21,no**2,1.0E0_realk,T_star_v3,no)
       call dgemm_acc_openacc_async(async_idx,'n','n',no,1,no**2,&
                & 1.0E0_realk,trip_ampl,no,int_normal_12,no**2,1.0E0_realk,T_star_v3,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,int_normal_21,int_normal_12,T_star_v3)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(no,kind=4),int(1,kind=4),int(no**2,kind=4),&
                             & -2.0E0_realk,c_loc(trip_ampl),int(no,kind=4),c_loc(int_normal_21),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v3),int(no,kind=4))
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(no,kind=4),int(1,kind=4),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(trip_ampl),int(no,kind=4),c_loc(int_normal_12),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v3),int(no,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_abc_12 (case 1)) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',no,1,no**2,&
                & -2.0E0_realk,trip_ampl,no,int_normal_21,no**2,1.0E0_realk,T_star_v3,no)
       call dgemm('n','n',no,1,no**2,&
                & 1.0E0_realk,trip_ampl,no,int_normal_12,no**2,1.0E0_realk,T_star_v3,no)
#endif

    end select TypeofContraction_abc_12

  end subroutine ccsdpt_contract_abc_12


  !> brief: do the first of the two contractions over 'cdk' (here: 'cd', 'k' is the summation in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_1 of driver. trip_ampl is the triples amplitude array.
  !> int_virt_tile is a v^3 tile determined by driver occ index
  !> tmp_g is a 3d work array
  subroutine ccsdpt_contract_ijk_211(oindex1,oindex2,oindex3,nv,no,&
       & T_star_o1o2,tmp_g,trip_ampl,int_virt_tile,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv,nv), target :: T_star_o1o2 ! T_star(:,:,oindex1,oindex2)
    real(realk), dimension(nv,nv,nv), target :: tmp_g,trip_ampl
    real(realk), dimension(nv,nv,nv), target, intent(in) :: int_virt_tile
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type, i,j,k
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_211: select case(contraction_type)

    case(0)

       ! note: here we collect contract over L_{dkbc} and g_{dkbc} in one go.

       ! reorder to obtain coulumb term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g)
#endif

       ! now contract coulumb term over 2 first indices
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,tmp_g,T_star_o1o2)
       call dgemm_acc_openacc_async(async_idx,'t','n',nv,nv,nv**2,&
!       call dgemm_acc('t','n',nv,nv,nv2,&
            & 1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,tmp_g,T_star_o1o2)
       stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(nv,kind=4),int(nv,kind=4),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(trip_ampl),int(nv**2,kind=4),c_loc(tmp_g),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o1o2),int(nv,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_ijk_211 - case 0-1) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('t','n',nv,nv,nv**2,&
            & 1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv)
#endif

       ! reorder to obtain exchange term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_virt_tile,nv,nv,nv,[3,1,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[3,1,2],0.0E0_realk,tmp_g)
#endif

       ! now contract exchange term over 2 first indices2
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,tmp_g,T_star_o1o2)
       call dgemm_acc_openacc_async(async_idx,'t','n',nv,nv,nv**2,&
!       call dgemm_acc('t','n',nv,nv,nv2,&
            -1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,tmp_g,T_star_o1o2)
       stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(nv,kind=4),int(nv,kind=4),int(nv**2,kind=4),&
                             & -1.0E0_realk,c_loc(trip_ampl),int(nv**2,kind=4),c_loc(tmp_g),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o1o2),int(nv,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_ijk_211 - case 0-2) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('t','n',nv,nv,nv**2,&
            -1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv)
#endif

    case(1)

       ! note: here we contract over L_{dkbc}.

       ! reorder to obtain coulumb term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g)
#endif

       ! now contract coulumb term over 2 first indices
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,tmp_g,T_star_o1o2)
       call dgemm_acc_openacc_async(async_idx,'t','n',nv,nv,nv**2,&
!       call dgemm_acc('t','n',nv,nv,nv2,&
            2.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,tmp_g,T_star_o1o2)
       stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(nv,kind=4),int(nv,kind=4),int(nv**2,kind=4),&
                             & 2.0E0_realk,c_loc(trip_ampl),int(nv**2,kind=4),c_loc(tmp_g),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o1o2),int(nv,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_ijk_211 - case 1-1) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('t','n',nv,nv,nv**2,&
            2.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv)
#endif

       ! reorder to obtain exchange term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_virt_tile,nv,nv,nv,[3,1,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[3,1,2],0.0E0_realk,tmp_g)
#endif

       ! now contract exchange term over 2 first indices
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,tmp_g,T_star_o1o2)
       call dgemm_acc_openacc_async(async_idx,'t','n',nv,nv,nv**2,&
!       call dgemm_acc('t','n',nv,nv,nv2,&
            -1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,tmp_g,T_star_o1o2)
       stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(nv,kind=4),int(nv,kind=4),int(nv**2,kind=4),&
                             & -1.0E0_realk,c_loc(trip_ampl),int(nv**2,kind=4),c_loc(tmp_g),int(nv**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o1o2),int(nv,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_ijk_211 - case 1-2) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('t','n',nv,nv,nv**2,&
            -1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o1o2,nv)
#endif

    end select TypeofContraction_211

  end subroutine ccsdpt_contract_ijk_211


  !> brief: do the first of the two contractions over 'cdk' (here: 'cd', 'k' is the summation in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_211(vindex1,vindex2,vindex3,no,nv,&
                               & int_virt_23,int_virt_32,T_star_v1,trip_ampl,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, no, nv
    real(realk), dimension(nv,no,no), target :: T_star_v1 ! T_star(:,:,:,vindex1)
    real(realk), dimension(nv,no), target :: int_virt_23, int_virt_32
    real(realk), dimension(no,no,no), target :: trip_ampl
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 221 and 222 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 221 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_abc_211: select case(contraction_type)

    case(0)

       ! now contract coulumb term over first index, then contract exchange term over first index 
       ! for this special case, we only have to subtract one coulumb term
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(int_virt_32,int_virt_23,trip_ampl,T_star_v1)
       call dgemm_acc_openacc_async(async_idx,'n','n',nv,no**2,no,1.0E0_realk,int_virt_32,nv,&
                      & trip_ampl,no,1.0E0_realk,T_star_v1,nv)
       call dgemm_acc_openacc_async(async_idx,'n','n',nv,no**2,no,-1.0E0_realk,int_virt_23,nv,&
                      & trip_ampl,no,1.0E0_realk,T_star_v1,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(int_virt_32,int_virt_23,trip_ampl,T_star_v1)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(nv,kind=4),int(no**2,kind=4),int(no,kind=4),&
                             & 1.0E0_realk,c_loc(int_virt_32),int(nv,kind=4),c_loc(trip_ampl),int(no,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v1),int(nv,kind=4))
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(nv,kind=4),int(no**2,kind=4),int(no,kind=4),&
                             & -1.0E0_realk,c_loc(int_virt_23),int(nv,kind=4),c_loc(trip_ampl),int(no,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v1),int(nv,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_abc_211 (case 0)) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',nv,no**2,no,1.0E0_realk,int_virt_32,nv,&
                      & trip_ampl,no,1.0E0_realk,T_star_v1,nv)
       call dgemm('n','n',nv,no**2,no,-1.0E0_realk,int_virt_23,nv,&
                      & trip_ampl,no,1.0E0_realk,T_star_v1,nv)
#endif

    case(1)
 
       ! now contract coulumb term over first index, next contract exchange term over first index
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(int_virt_32,int_virt_23,trip_ampl,T_star_v1)
       call dgemm_acc_openacc_async(async_idx,'n','n',nv,no**2,no,2.0E0_realk,int_virt_32,nv,&
                      & trip_ampl,no,1.0E0_realk,T_star_v1,nv)
       call dgemm_acc_openacc_async(async_idx,'n','n',nv,no**2,no,-1.0E0_realk,int_virt_23,nv,&
                      & trip_ampl,no,1.0E0_realk,T_star_v1,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(int_virt_32,int_virt_23,trip_ampl,T_star_v1)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(nv,kind=4),int(no**2,kind=4),int(no,kind=4),&
                             & 2.0E0_realk,c_loc(int_virt_32),int(nv,kind=4),c_loc(trip_ampl),int(no,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v1),int(nv,kind=4))
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(nv,kind=4),int(no**2,kind=4),int(no,kind=4),&
                             & -1.0E0_realk,c_loc(int_virt_23),int(nv,kind=4),c_loc(trip_ampl),int(no,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v1),int(nv,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_abc_211 (case 1)) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',nv,no**2,no,2.0E0_realk,int_virt_32,nv,&
                      & trip_ampl,no,1.0E0_realk,T_star_v1,nv)
       call dgemm('n','n',nv,no**2,no,-1.0E0_realk,int_virt_23,nv,&
                      & trip_ampl,no,1.0E0_realk,T_star_v1,nv)
#endif

    end select TypeofContraction_abc_211

  end subroutine ccsdpt_contract_abc_211


  !> brief: do the second of the two contractions over 'cdk' (here: 'cd', 'k' is the summation in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_1 of driver. trip_ampl is the triples amplitude array.
  !> int_virt_tile is a v^3 tile determined by driver occ index
  !> tmp_g is a 3d work array
  subroutine ccsdpt_contract_ijk_212(oindex1,oindex2,oindex3,nv,no,&
       & T_star_o3o2,tmp_g,trip_ampl,int_virt_tile,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv,nv), target :: T_star_o3o2 ! T_star(:,:,oindex3,oindex2)
    real(realk), dimension(nv,nv,nv), target, intent(in) :: int_virt_tile
    real(realk), dimension(nv,nv,nv), target :: tmp_g,trip_ampl
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! reorder to obtain coulumb term 
#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
    call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g)
#endif

    ! now contract coulumb term over 2 first indices
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,tmp_g,T_star_o3o2)
    call dgemm_acc_openacc_async(async_idx,'t','n',nv,nv,nv**2,&
!    call dgemm_acc('t','n',nv,nv,nv2,&
         & -1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o3o2,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,tmp_g,T_star_o3o2)
    stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(nv,kind=4),int(nv,kind=4),int(nv**2,kind=4),&
                          & -1.0E0_realk,c_loc(trip_ampl),int(nv**2,kind=4),c_loc(tmp_g),int(nv**2,kind=4),&
                          & 1.0E0_realk,c_loc(T_star_o3o2),int(nv,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_contract_ijk_212) = ",stat
!       stop
!    end if

#endif
#else
    call dgemm('t','n',nv,nv,nv**2,&
         & -1.0E0_realk,trip_ampl,nv**2,tmp_g,nv**2,1.0E0_realk,T_star_o3o2,nv)
#endif

  end subroutine ccsdpt_contract_ijk_212


  !> brief: do the second of the two contractions over 'cdk' (here: 'cd', 'k' is the summation in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_212(vindex1,vindex2,vindex3,no,nv,int_virt_12,T_star_v3,trip_ampl,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, no, nv
    real(realk), dimension(nv,no,no), target :: T_star_v3 ! T_star(:,:,:,vindex3)
    real(realk), dimension(nv,no), target :: int_virt_12
    real(realk), dimension(no,no,no), target :: trip_ampl
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(int_virt_12,trip_ampl,T_star_v3)
    call dgemm_acc_openacc_async(async_idx,'n','n',nv,no**2,no,-1.0E0_realk,int_virt_12,nv,&
                   & trip_ampl,no,1.0E0_realk,T_star_v3,nv)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(int_virt_12,trip_ampl,T_star_v3)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(nv,kind=4),int(no**2,kind=4),int(no,kind=4),&
                          & -1.0E0_realk,c_loc(int_virt_12),int(nv,kind=4),c_loc(trip_ampl),int(no,kind=4),&
                          & 1.0E0_realk,c_loc(T_star_v3),int(nv,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_contract_abc_212) = ",stat
!       stop
!    end if

#endif
#else
    call dgemm('n','n',nv,no**2,no,-1.0E0_realk,int_virt_12,nv,&
                   & trip_ampl,no,1.0E0_realk,T_star_v3,nv)
#endif

  end subroutine ccsdpt_contract_abc_212


  !> brief: do the first of the two contractions over 'ckl' (here: 'c', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_2 of driver. trip_ampl is the triples amplitud array.
  subroutine ccsdpt_contract_ijk_221(oindex1,oindex2,oindex3,no,nv,&
                               & int_occ_23,int_occ_32,T_star_o1,trip_ampl,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), dimension(no,nv,nv), target :: T_star_o1 ! T_star(:,:,:,oindex1)
    real(realk), dimension(no,nv), target :: int_occ_23, int_occ_32
    real(realk), dimension(nv,nv,nv), target :: trip_ampl
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 221 and 222 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 221 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_221: select case(contraction_type)

    case(0)

       ! now contract coulumb term over first index, then contract exchange term over first index 
       ! for this special case, we only have to subtract one coulumb term
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(int_occ_32,int_occ_23,trip_ampl,T_star_o1)
       call dgemm_acc_openacc_async(async_idx,'n','n',no,nv**2,nv,-1.0E0_realk,int_occ_32,no,&
!       call dgemm_acc('n','n',no,nv2,nv,-1.0E0_realk,int_occ_32,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star_o1,no)
       call dgemm_acc_openacc_async(async_idx,'n','n',no,nv**2,nv,1.0E0_realk,int_occ_23,no,&
!       call dgemm_acc('n','n',no,nv2,nv,1.0E0_realk,int_occ_23,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star_o1,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(int_occ_32,int_occ_23,trip_ampl,T_star_o1)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(no,kind=4),int(nv**2,kind=4),int(nv,kind=4),&
                             & -1.0E0_realk,c_loc(int_occ_32),int(no,kind=4),c_loc(trip_ampl),int(nv,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o1),int(no,kind=4))
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(no,kind=4),int(nv**2,kind=4),int(nv,kind=4),&
                             & 1.0E0_realk,c_loc(int_occ_23),int(no,kind=4),c_loc(trip_ampl),int(nv,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o1),int(no,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_ijk_221 - case 0) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',no,nv**2,nv,-1.0E0_realk,int_occ_32,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star_o1,no)
       call dgemm('n','n',no,nv**2,nv,1.0E0_realk,int_occ_23,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star_o1,no)
#endif

    case(1)
 
       ! now contract coulumb term over first index, next contract exchange term over first index
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(int_occ_32,int_occ_23,trip_ampl,T_star_o1)
       call dgemm_acc_openacc_async(async_idx,'n','n',no,nv**2,nv,-2.0E0_realk,int_occ_32,no,&
!       call dgemm_acc('n','n',no,nv2,nv,-2.0E0_realk,int_occ_32,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star_o1,no)
       call dgemm_acc_openacc_async(async_idx,'n','n',no,nv**2,nv,1.0E0_realk,int_occ_23,no,&
!       call dgemm_acc('n','n',no,nv2,nv,1.0E0_realk,int_occ_23,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star_o1,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(int_occ_32,int_occ_23,trip_ampl,T_star_o1)
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(no,kind=4),int(nv**2,kind=4),int(nv,kind=4),&
                             & -2.0E0_realk,c_loc(int_occ_32),int(no,kind=4),c_loc(trip_ampl),int(nv,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o1),int(no,kind=4))
       stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(no,kind=4),int(nv**2,kind=4),int(nv,kind=4),&
                             & 1.0E0_realk,c_loc(int_occ_23),int(no,kind=4),c_loc(trip_ampl),int(nv,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_o1),int(no,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_ijk_221 - case 1) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('n','n',no,nv**2,nv,-2.0E0_realk,int_occ_32,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star_o1,no)
       call dgemm('n','n',no,nv**2,nv,1.0E0_realk,int_occ_23,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star_o1,no)
#endif

    end select TypeofContraction_221

  end subroutine ccsdpt_contract_ijk_221


  !> brief: do the first of the two contractions over 'ckl' (here: 'c', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_221(vindex1,vindex2,vindex3,nv,no,&
       & T_star_v1v2,tmp_g,trip_ampl,int_occ_tile,special,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, nv, no
    real(realk), dimension(no,no), target :: T_star_v1v2 ! T_star(:,:,vindex1,vindex2)
    real(realk), dimension(no,no,no), target :: tmp_g,trip_ampl
    real(realk), dimension(no,no,no), target, intent(in) :: int_occ_tile
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    TypeofContraction_abc_221: select case(contraction_type)

    case(0)

       ! reorder to obtain coulumb term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,tmp_g,T_star_v1v2)
       call dgemm_acc_openacc_async(async_idx,'t','n',no,no,no**2,&
            & -1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,tmp_g,T_star_v1v2)
       stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(no,kind=4),int(no,kind=4),int(no**2,kind=4),&
                             & -1.0E0_realk,c_loc(trip_ampl),int(no**2,kind=4),c_loc(tmp_g),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v1v2),int(no,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_abc_221 (case 0-1)) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('t','n',no,no,no**2, &
            & -1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no)
#endif

       ! reorder to obtain exchange term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_occ_tile,no,no,no,[3,1,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_occ_tile,no,no,no,[3,1,2],0.0E0_realk,tmp_g)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,tmp_g,T_star_v1v2)
       call dgemm_acc_openacc_async(async_idx,'t','n',no,no,no**2, &
            1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,tmp_g,T_star_v1v2)
       stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(no,kind=4),int(no,kind=4),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(trip_ampl),int(no**2,kind=4),c_loc(tmp_g),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v1v2),int(no,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_abc_221 (case 0-2)) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('t','n',no,no,no**2,&
            1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no)
#endif

    case(1)

       ! reorder to obtain coulumb term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,tmp_g,T_star_v1v2)
       call dgemm_acc_openacc_async(async_idx,'t','n',no,no,no**2,&
            -2.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,tmp_g,T_star_v1v2)
       stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(no,kind=4),int(no,kind=4),int(no**2,kind=4),&
                             & -2.0E0_realk,c_loc(trip_ampl),int(no**2,kind=4),c_loc(tmp_g),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v1v2),int(no,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_abc_221 (case 1-1)) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('t','n',no,no,no**2,&
            -2.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no)
#endif

       ! reorder to obtain exchange term
#ifdef VAR_OPENACC
       call array_reorder_3d_acc(1.0E0_realk,int_occ_tile,no,no,no,[3,1,2],0.0E0_realk,tmp_g,async_idx)
#else
       call array_reorder_3d(1.0E0_realk,int_occ_tile,no,no,no,[3,1,2],0.0E0_realk,tmp_g)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,tmp_g,T_star_v1v2)
       call dgemm_acc_openacc_async(async_idx,'t','n',no,no,no**2,&
            1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,tmp_g,T_star_v1v2)
       stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(no,kind=4),int(no,kind=4),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(trip_ampl),int(no**2,kind=4),c_loc(tmp_g),int(no**2,kind=4),&
                             & 1.0E0_realk,c_loc(T_star_v1v2),int(no,kind=4))
!$acc end host_data

!       if (stat .ne. 0 ) then
!          print *, "stat (ccsdpt_contract_abc_221 (case 1-2)) = ",stat
!          stop
!       end if

#endif
#else
       call dgemm('t','n',no,no,no**2, &
            1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v1v2,no)
#endif

    end select TypeofContraction_abc_221

  end subroutine ccsdpt_contract_abc_221


  !> brief: do the second of the two contractions over 'ckl' (here: 'c', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_2 of driver. trip_ampl is the triples amplitud array.
  subroutine ccsdpt_contract_ijk_222(oindex1,oindex2,oindex3,no,nv,int_occ_12,T_star_o3,trip_ampl,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), dimension(no,nv,nv), target :: T_star_o3 ! T_star(:,:,:,oindex3)
    real(realk), dimension(no,nv), target :: int_occ_12
    real(realk), dimension(nv,nv,nv), target :: trip_ampl
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! contract coulumb term over first index
#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(int_occ_12,trip_ampl,T_star_o3)
    call dgemm_acc_openacc_async(async_idx,'n','n',no,nv**2,nv,1.0E0_realk,int_occ_12,no,&
!    call dgemm_acc('n','n',no,nv2,nv,1.0E0_realk,int_occ_12,no,&
                   & trip_ampl,nv,1.0E0_realk,T_star_o3,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(int_occ_12,trip_ampl,T_star_o3)
    stat = cublasDgemm_v2(cublas_handle,int(0,kind=4),int(0,kind=4),int(no,kind=4),int(nv**2,kind=4),int(nv,kind=4),&
                          & 1.0E0_realk,c_loc(int_occ_12),int(no,kind=4),c_loc(trip_ampl),int(nv,kind=4),&
                          & 1.0E0_realk,c_loc(T_star_o3),int(no,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_contract_ijk_222) = ",stat
!       stop
!    end if

#endif
#else
    call dgemm('n','n',no,nv**2,nv,1.0E0_realk,int_occ_12,no,&
                   & trip_ampl,nv,1.0E0_realk,T_star_o3,no)
#endif

  end subroutine ccsdpt_contract_ijk_222


  !> brief: do the second of the two contractions over 'ckl' (here: 'c', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: april 2014
  subroutine ccsdpt_contract_abc_222(vindex1,vindex2,vindex3,nv,no,&
       & T_star_v3v2,tmp_g,trip_ampl,int_occ_tile,async_idx,cublas_handle)

    implicit none
    !> input
    integer, intent(in) :: vindex1, vindex2, vindex3, nv, no
    real(realk), dimension(no,no), target :: T_star_v3v2 ! T_star(:,:,vindex3,vindex2)
    real(realk), dimension(no,no,no), target, intent(in) :: int_occ_tile
    real(realk), dimension(no,no,no), target :: tmp_g,trip_ampl
#ifdef VAR_OPENACC
    integer(kind=acc_handle_kind) :: async_idx
#else
    integer :: async_idx
#endif
    type(c_ptr) :: cublas_handle
    integer*4 :: stat

    ! reorder
#ifdef VAR_OPENACC
    call array_reorder_3d_acc(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g,async_idx)
#else
    call array_reorder_3d(1.0E0_realk,int_occ_tile,no,no,no,[1,3,2],0.0E0_realk,tmp_g)
#endif

#ifdef VAR_OPENACC
#if defined(VAR_CRAY) && !defined(VAR_CUBLAS)
!$acc host_data use_device(trip_ampl,tmp_g,T_star_v3v2)
    call dgemm_acc_openacc_async(async_idx,'t','n',no,no,no**2,&
         & 1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v3v2,no)
!$acc end host_data
#elif defined(VAR_CUBLAS)

!$acc host_data use_device(trip_ampl,tmp_g,T_star_v3v2)
    stat = cublasDgemm_v2(cublas_handle,int(1,kind=4),int(0,kind=4),int(no,kind=4),int(no,kind=4),int(no**2,kind=4),&
                          & 1.0E0_realk,c_loc(trip_ampl),int(no**2,kind=4),c_loc(tmp_g),int(no**2,kind=4),&
                          & 1.0E0_realk,c_loc(T_star_v3v2),int(no,kind=4))
!$acc end host_data

!    if (stat .ne. 0 ) then
!       print *, "stat (ccsdpt_contract_abc_222) = ",stat
!       stop
!    end if

#endif
#else
    call dgemm('t','n',no,no,no**2,&
         & 1.0E0_realk,trip_ampl,no**2,tmp_g,no**2,1.0E0_realk,T_star_v3v2,no)
#endif

  end subroutine ccsdpt_contract_abc_222


  !> \brief: calculate E[5] contribution to fragment decnp-ccsd(t) energy correction
  !> \author: Pablo Baudin (from Janus Eriksen and Kasper Kristensen)
  !> \date: Mar 2015
  subroutine ccsdpt_decnp_e5_frag(MyFragment,ccsd_singles,ccsdpt_singles)

    implicit none

    !> fragment info
    type(decfrag), intent(inout) :: MyFragment
    ! ccsd and ccsd(t) singles amplitudes
    type(tensor), intent(inout) :: ccsd_singles, ccsdpt_singles
    real(realk) :: energy_tmp, ccsdpt_e5
    logical :: SEC_occ(MyFragment%noccAOS), SEC_virt(MyFragment%nvirtAOS)
    integer :: noccAOS,nvirtAOS,i,a
    !> which partitioning schemes?
    logical :: do_occ, do_virt

    noccAOS = MyFragment%noccAOS
    nvirtAOS = MyFragment%nvirtAOS

    ! Sanity check
    if(MyFragment%nEOSatoms/=1) then
       print *, 'nEOSatoms ',MyFragment%nEOSatoms
       call lsquit('ccsdpt_decnp_e5_frag called with wrong number of EOS atoms!',-1)
    end if

    ! Determine which occ and virt orbitals to include in energy contributions
    ! based on SECONDARY assignment.
    call secondary_assigning(MyFragment,SEC_occ,SEC_virt)


    ! ***********************
    !   do E[5] energy part
    ! ***********************

    ! init energy reals to be on the safe side.
    ! note: OccEnergyPT and VirtEnergyPT have been initialized in the e4 routine.
    MyFragment%energies(FRAGMODEL_OCCpT5) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_VIRTpT5) = 0.0E0_realk

    do_occ = .false.
    do_virt = .false.
    if ((.not. DECinfo%OnlyOccPart) .and. (.not. DECinfo%OnlyVirtPart)) then
       do_occ = .true.
       do_virt = .true.
    else if (DECinfo%OnlyOccPart .and. (.not. DECinfo%OnlyVirtPart)) then
       do_occ = .true.
    else if (DECinfo%OnlyVirtPart .and. (.not. DECinfo%OnlyOccPart)) then
       do_virt = .true.
    end if

    ! *******************************
    ! do occupied partitioning scheme
    ! *******************************

    if (do_occ) then
       ! init temp energy
       ccsdpt_e5 = 0.0E0_realk

       iloop: do i=1,noccAOS
          ! Only include contribution if consistent with secondary assignment
          if(.not. SEC_occ(i)) cycle iloop
          do a=1,nvirtAOS

             energy_tmp = ccsd_singles%elm2(a,i) * ccsdpt_singles%elm2(a,i)
             ccsdpt_e5 = ccsdpt_e5 + energy_tmp

          end do
       end do iloop
       MyFragment%energies(FRAGMODEL_OCCpT5) = 2.0E0_realk * ccsdpt_e5

       ! insert into occ. part. scheme part
       MyFragment%energies(FRAGMODEL_OCCpT) = MyFragment%energies(FRAGMODEL_OCCpT) &
          & + MyFragment%energies(FRAGMODEL_OCCpT5)
    end if

    ! *********************************
    ! do virtupied partitioning scheme
    ! *********************************

    if (do_virt) then 
       ! init temp energy
       ccsdpt_e5 = 0.0E0_realk

       do i=1,noccAOS
          aloop: do a=1,nvirtAOS
             if(.not. SEC_virt(a)) cycle aloop

             energy_tmp = ccsd_singles%elm2(a,i) * ccsdpt_singles%elm2(a,i)
             ccsdpt_e5 = ccsdpt_e5 + energy_tmp

          end do aloop
       end do
       MyFragment%energies(FRAGMODEL_VIRTpT5) = 2.0E0_realk * ccsdpt_e5

       ! insert into virt. part. scheme part
       MyFragment%energies(FRAGMODEL_VIRTpT) = MyFragment%energies(FRAGMODEL_VIRTpT) &
          & + MyFragment%energies(FRAGMODEL_VIRTpT5)
    end if

  end subroutine ccsdpt_decnp_e5_frag


  !> \brief: calculate E[5] contribution to single fragment ccsd(t) energy correction
  !> \author: Janus Eriksen and Kasper Kristensen
  !> \date: september 2012
  subroutine ccsdpt_energy_e5_frag(MyFragment,ccsd_singles,ccsdpt_singles)

    implicit none

    !> fragment info
    type(decfrag), intent(inout) :: MyFragment
    ! ccsd and ccsd(t) singles amplitudes
    type(tensor), intent(inout) :: ccsd_singles, ccsdpt_singles
    real(realk) :: energy_tmp, ccsdpt_e5
    logical :: SEC_occ(MyFragment%noccAOS), SEC_virt(MyFragment%nvirtAOS)
    integer :: noccAOS,nvirtAOS,i,a

    noccAOS = MyFragment%noccAOS
    nvirtAOS = MyFragment%nvirtAOS

    ! Sanity check
    if(MyFragment%nEOSatoms/=1) then
       print *, 'nEOSatoms ',MyFragment%nEOSatoms
       call lsquit('ccsdpt_energy_e5_frag called with wrong number of EOS atoms!',-1)
    end if

    ! Determine which occ and virt orbitals to include in energy contributions
    ! based on SECONDARY assignment.
    call secondary_assigning(MyFragment,SEC_occ,SEC_virt)


    ! ***********************
    !   do E[5] energy part
    ! ***********************

    ! init energy reals to be on the safe side.
    ! note: OccEnergyPT and VirtEnergyPT have been initialized in the e4 routine.
    MyFragment%energies(FRAGMODEL_OCCpT5) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_VIRTpT5) = 0.0E0_realk

    ! init temp energy
    ccsdpt_e5 = 0.0E0_realk

    iloop: do i=1,noccAOS
       ! Only include contribution if consistent with secondary assignment
       if(.not. SEC_occ(i)) cycle iloop
       aloop: do a=1,nvirtAOS
          if(.not. SEC_virt(a)) cycle aloop

          energy_tmp = ccsd_singles%elm2(a,i) * ccsdpt_singles%elm2(a,i)
          ccsdpt_e5 = ccsdpt_e5 + energy_tmp

       end do aloop
    end do iloop
    MyFragment%energies(FRAGMODEL_OCCpT5) = 2.0E0_realk * ccsdpt_e5

    ! insert into occ. part. scheme part
    MyFragment%energies(FRAGMODEL_OCCpT) = MyFragment%energies(FRAGMODEL_OCCpT) &
         & + MyFragment%energies(FRAGMODEL_OCCpT5)



    ! *********************************
    ! do virtupied partitioning scheme
    ! *********************************

    ! singles contribution is the same as in occupied partitioning scheme
    MyFragment%energies(FRAGMODEL_VIRTpT) = MyFragment%energies(FRAGMODEL_VIRTpT) &
         & + MyFragment%energies(FRAGMODEL_OCCpT5)
    ! insert into virt_e5 part
    MyFragment%energies(FRAGMODEL_VIRTpT5) = MyFragment%energies(FRAGMODEL_VIRTpT5) &
         & + MyFragment%energies(FRAGMODEL_OCCpT5)

  end subroutine ccsdpt_energy_e5_frag


  !> \brief: calculate E[5] contribution to pair fragment ccsd(t) energy correction
  !> \author: Janus Eriksen
  !> \date: september 2012
  subroutine ccsdpt_energy_e5_pair(PairFragment,ccsd_singles,ccsdpt_singles)

    implicit none

    !> fragment info
    type(decfrag), intent(inout) :: PairFragment
    ! ccsd and ccsd(t) singles amplitudes
    type(tensor), intent(inout) :: ccsd_singles, ccsdpt_singles
    real(realk) :: energy_tmp, ccsdpt_e5
    logical :: SEC_occ(PairFragment%noccAOS), SEC_virt(PairFragment%nvirtAOS)
    integer :: noccAOS,nvirtAOS,i,a,atomi,atoma

    noccAOS = PairFragment%noccAOS
    nvirtAOS = PairFragment%nvirtAOS

    ! Sanity check
    if(PairFragment%nEOSatoms/=2) then
       print *, 'nEOSatoms ',PairFragment%nEOSatoms
       call lsquit('ccsdpt_energy_e5_pair called with wrong number of EOS atoms!',-1)
    end if

    ! Determine which occ and virt orbitals to include in energy contributions
    ! based on SECONDARY assignment.
    call secondary_assigning(PairFragment,SEC_occ,SEC_virt)


    ! ***********************
    !   do E[5] energy part
    ! ***********************

    ! init energy reals to be on the safe side.
    ! note: OccEnergyPT and VirtEnergyPT have been initialized in the e4 routine.
    PairFragment%energies(FRAGMODEL_OCCpT5) = 0.0E0_realk
    PairFragment%energies(FRAGMODEL_VIRTpT5) = 0.0E0_realk

    ! init temp energy
    ccsdpt_e5 = 0.0E0_realk

    iloop: do i=1,noccAOS
       ! Only include contribution if consistent with secondary assignment
       if(.not. SEC_occ(i)) cycle iloop
       atomi = PairFragment%occAOSorb(i)%secondaryatom
       aloop: do a=1,nvirtAOS
          if(.not. SEC_virt(a)) cycle aloop
          atoma = PairFragment%virtAOSorb(a)%secondaryatom

          ! Only include if atomi/=atoma
          ! (the atomi=atoma contributions were included for atomic fragments)
          if(atomi/=atoma) then
             energy_tmp = ccsd_singles%elm2(a,i) * ccsdpt_singles%elm2(a,i)
             ccsdpt_e5 = ccsdpt_e5 + energy_tmp
          end if

       end do aloop
    end do iloop
    PairFragment%energies(FRAGMODEL_OCCpT5) = 2.0E0_realk * ccsdpt_e5

    ! insert into occ. part. scheme part
    PairFragment%energies(FRAGMODEL_OCCpT) = PairFragment%energies(FRAGMODEL_OCCpT) &
         & + PairFragment%energies(FRAGMODEL_OCCpT5)

    ! *********************************
    ! do virtupied partitioning scheme
    ! *********************************

    ! singles contribution is the same as in occupied partitioning scheme
    PairFragment%energies(FRAGMODEL_VIRTpT) = PairFragment%energies(FRAGMODEL_VIRTpT) &
         & + PairFragment%energies(FRAGMODEL_OCCpT5)
    ! insert into virt_e5 part
    PairFragment%energies(FRAGMODEL_VIRTpT5) = PairFragment%energies(FRAGMODEL_VIRTpT5) &
         & + PairFragment%energies(FRAGMODEL_OCCpT5)


  end subroutine ccsdpt_energy_e5_pair



  !> \brief: calculate E[4] contribution to single fragment ccsd(t) energy correction
  !> \author: Janus Eriksen
  !> \date: september 2012
  subroutine ccsdpt_energy_e4_frag(MyFragment,ccsd_doubles,ccsdpt_doubles,&
                             & occ_contribs,virt_contribs,fragopt_pT)

    implicit none

    !> fragment info
    type(decfrag), intent(inout) :: MyFragment
    ! ccsd and ccsd(t) doubles amplitudes
    type(tensor), intent(inout) :: ccsd_doubles, ccsdpt_doubles
    !> is this called from inside the ccsd(t) fragment optimization routine?
    logical, optional, intent(in) :: fragopt_pT
    !> incomming orbital contribution vectors
    real(realk), intent(inout) :: occ_contribs(MyFragment%noccAOS), virt_contribs(MyFragment%nvirtAOS)
    !> integers
    integer :: nocc_eos, nocc_aos, nvirt_eos, nvirt_aos, i,j,a,b, i_eos, j_eos, a_eos, b_eos
    !> energy reals
    real(realk) :: energy_tmp, energy_res_cou, energy_res_exc
    !> which partitioning schemes?
    logical :: do_occ, do_virt

    ! init dimensions
    nocc_eos = MyFragment%noccEOS
    nvirt_eos = MyFragment%nvirtEOS
    nocc_aos = MyFragment%noccAOS
    nvirt_aos = MyFragment%nvirtAOS

    ! **************************************************************
    ! ************** do energy for single fragment *****************
    ! **************************************************************

    ! ***********************
    !   do E[4] energy part
    ! ***********************

    ! init energy reals to be on the safe side
    ! note: OccEnergyPT and VirtEnergyPT is also initialized from in here
    !       as this (e4) routine is called before the e5 routine
    MyFragment%energies(FRAGMODEL_OCCpT) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_VIRTpT) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_OCCpT4) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_VIRTpT4) = 0.0E0_realk

    do_occ = .false.
    do_virt = .false.

    if ((.not. DECinfo%OnlyOccPart) .and. (.not. DECinfo%OnlyVirtPart)) then

       do_occ = .true.
       do_virt = .true.

    else if (DECinfo%OnlyOccPart .and. (.not. DECinfo%OnlyVirtPart)) then

       do_occ = .true.

    else if (DECinfo%OnlyVirtPart .and. (.not. DECinfo%OnlyOccPart)) then

       do_virt = .true.

    end if

    ! *******************************
    ! do occupied partitioning scheme
    ! *******************************

    if (do_occ) then

       energy_res_cou = 0.0E0_realk
       energy_res_exc = 0.0E0_realk
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,i_eos,j,j_eos,a,b,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc_eos,nvirt_aos,MyFragment),&
       !$OMP REDUCTION(+:energy_res_cou),REDUCTION(+:virt_contribs)
       do j=1,nocc_eos
          j_eos = MyFragment%idxo(j)
          do i=1,nocc_eos
             i_eos = MyFragment%idxo(i)
   
             do b=1,nvirt_aos
                do a=1,nvirt_aos
   
                   energy_tmp = 4.0E0_realk * ccsd_doubles%elm4(a,b,i_eos,j_eos) &
                                  & * ccsdpt_doubles%elm4(a,b,i_eos,j_eos)
                   energy_res_cou = energy_res_cou + energy_tmp
   
                   ! update contribution from aos orbital a
                   virt_contribs(a) = virt_contribs(a) + energy_tmp
   
                   ! update contribution from aos orbital b 
                   ! (only if different from aos orbital a to avoid double counting)
                   if (a .ne. b) virt_contribs(b) = virt_contribs(b) + energy_tmp
   
                end do
             end do
   
          end do
       end do
       !$OMP END PARALLEL DO
   
       ! reorder from (a,b,i,j) to (a,b,j,i)
       call tensor_reorder(ccsd_doubles,[1,2,4,3])
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,i_eos,j,j_eos,a,b,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc_eos,nvirt_aos,MyFragment),&
       !$OMP REDUCTION(+:energy_res_exc),REDUCTION(+:virt_contribs)
       do j=1,nocc_eos
          j_eos = MyFragment%idxo(j)
          do i=1,nocc_eos
             i_eos = MyFragment%idxo(i)
   
             do b=1,nvirt_aos
                do a=1,nvirt_aos
   
                   energy_tmp = 2.0E0_realk * ccsd_doubles%elm4(a,b,i_eos,j_eos) &
                                  & * ccsdpt_doubles%elm4(a,b,i_eos,j_eos)
                   energy_res_exc = energy_res_exc - energy_tmp
   
                   ! update contribution from aos orbital a
                   virt_contribs(a) = virt_contribs(a) - energy_tmp
   
                   ! update contribution from aos orbital b 
                   ! (only if different from aos orbital a to avoid double counting)
                   if (a .ne. b) virt_contribs(b) = virt_contribs(b) - energy_tmp
   
                end do
             end do
   
          end do
       end do
       !$OMP END PARALLEL DO
   
       !get total fourth--order energy contribution
       MyFragment%energies(FRAGMODEL_OCCpT4) = energy_res_cou + energy_res_exc
   
       ! insert into occ. part. scheme part
       MyFragment%energies(FRAGMODEL_OCCpT) = MyFragment%energies(FRAGMODEL_OCCpT) &
         & + MyFragment%energies(FRAGMODEL_OCCpT4)

    end if

    ! *********************************
    ! do virtupied partitioning scheme
    ! *********************************

    if (do_virt) then

       if (do_occ .and. do_virt) then

          ! initially, reorder ccsd_doubles and ccsdpt_doubles
          ! ccsd_doubles from from (a,b,j,i) sequence to (j,i,a,b) sequence
          ! ccsdpt_doubles from from (a,b,i,j) sequence to (i,j,a,b) sequence
          call tensor_reorder(ccsd_doubles,[3,4,1,2])
          call tensor_reorder(ccsdpt_doubles,[3,4,1,2])

       else

          ! initially, reorder ccsd_doubles and ccsdpt_doubles
          ! ccsd_doubles from from (a,b,i,j) sequence to (j,i,a,b) sequence
          ! ccsdpt_doubles from from (a,b,i,j) sequence to (i,j,a,b) sequence
          call tensor_reorder(ccsd_doubles,[4,3,1,2])
          call tensor_reorder(ccsdpt_doubles,[3,4,1,2])

       end if

       energy_res_cou = 0.0E0_realk
       energy_res_exc = 0.0E0_realk
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,a_eos,b,b_eos,i,j,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nvirt_eos,nocc_aos,MyFragment),&
       !$OMP REDUCTION(+:energy_res_exc),REDUCTION(+:occ_contribs)
       do b=1,nvirt_eos
          b_eos = MyFragment%idxu(b)
          do a=1,nvirt_eos
             a_eos = MyFragment%idxu(a)
   
             do j=1,nocc_aos
                do i=1,nocc_aos
   
                   energy_tmp = 2.0E0_realk * ccsd_doubles%elm4(i,j,a_eos,b_eos) &
                                  & * ccsdpt_doubles%elm4(i,j,a_eos,b_eos)
                   energy_res_exc = energy_res_exc - energy_tmp
   
                   ! update contribution from aos orbital i
                   occ_contribs(i) = occ_contribs(i) - energy_tmp
   
                   ! update contribution from aos orbital j 
                   ! (only if different from aos orbital i to avoid double counting)
                   if (i .ne. j) occ_contribs(j) = occ_contribs(j) - energy_tmp
   
                end do
             end do
   
          end do
       end do
       !$OMP END PARALLEL DO
   
       ! reorder form (j,i,a,b) to (i,j,a,b)
       call tensor_reorder(ccsd_doubles,[2,1,3,4])
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,a_eos,b,b_eos,i,j,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nvirt_eos,nocc_aos,MyFragment),&
       !$OMP REDUCTION(+:energy_res_cou),REDUCTION(+:occ_contribs)
       do b=1,nvirt_eos
          b_eos = MyFragment%idxu(b)
          do a=1,nvirt_eos
             a_eos = MyFragment%idxu(a)
   
             do j=1,nocc_aos
                do i=1,nocc_aos
   
                   energy_tmp = 4.0E0_realk * ccsd_doubles%elm4(i,j,a_eos,b_eos) &
                                  & * ccsdpt_doubles%elm4(i,j,a_eos,b_eos)
                   energy_res_cou = energy_res_cou + energy_tmp
   
                   ! update contribution from aos orbital i
                   occ_contribs(i) = occ_contribs(i) + energy_tmp
   
                   ! update contribution from aos orbital j 
                   ! (only if different from aos orbital i to avoid double counting)
                   if (i .ne. j) occ_contribs(j) = occ_contribs(j) + energy_tmp
   
                end do
             end do
   
          end do
       end do
       !$OMP END PARALLEL DO
   
       !get total fourth--order energy contribution
       MyFragment%energies(FRAGMODEL_VIRTpT4) = energy_res_cou + energy_res_exc
   
       ! insert into virt. part. scheme part
       MyFragment%energies(FRAGMODEL_VIRTpT) = MyFragment%energies(FRAGMODEL_VIRTpT) &
          & + MyFragment%energies(FRAGMODEL_VIRTpT4)

    end if

    ! ******************************
    !   done with E[4] energy part
    ! ******************************

    ! ************************************************************************
    !   as we need to reuse the ccsd doubles in the fragment optimization,
    !   we here reorder back into (a,i,b,j) sequence IF fragopt_pT == .true. 
    ! ************************************************************************

    if (present(fragopt_pT)) then

       if (do_occ .and. (.not. do_virt)) then

          ! reorder from (a,b,j,i) to (a,i,b,j)
          if (fragopt_pT) call tensor_reorder(ccsd_doubles,[1,4,2,3])

       else if (do_virt .and. (.not. do_virt)) then

          ! reorder from (i,j,a,b) to (a,i,b,j)
          if (fragopt_pT) call tensor_reorder(ccsd_doubles,[3,1,4,2])

       else if (do_occ .and. do_virt) then

          ! reorder from (i,j,a,b) to (a,i,b,j)
          if (fragopt_pT) call tensor_reorder(ccsd_doubles,[3,1,4,2])

       end if

    end if

    ! *******************************************************************
    ! ************** done w/ energy for single fragment *****************
    ! *******************************************************************

  end subroutine ccsdpt_energy_e4_frag


  !> \brief: calculate E[4] contribution to fragment decnp-ccsd(t) energy correction
  !> \author: Pablo Baudin (based on Janus's routine)
  !> \date: Mar 2015
  subroutine ccsdpt_decnp_e4_frag(MyFragment,ccsd_doubles,ccsdpt_doubles,&
                             & occ_contribs,virt_contribs,fragopt_pT)

    implicit none

    !> fragment info
    type(decfrag), intent(inout) :: MyFragment
    ! ccsd and ccsd(t) doubles amplitudes
    type(tensor), intent(inout) :: ccsd_doubles, ccsdpt_doubles
    !> is this called from inside the ccsd(t) fragment optimization routine?
    logical, optional, intent(in) :: fragopt_pT
    !> incomming orbital contribution vectors
    real(realk), intent(inout) :: occ_contribs(MyFragment%noccAOS), virt_contribs(MyFragment%nvirtAOS)
    !> integers
    integer :: nocc_eos, nocc_aos, nvirt_eos, nvirt_aos, i,j,a,b, i_eos, j_eos, a_eos, b_eos
    !> energy reals
    real(realk) :: energy_tmp, energy_res_cou, energy_res_exc
    !> which partitioning schemes?
    logical :: do_occ, do_virt

    ! init dimensions
    nocc_eos = MyFragment%noccEOS
    nvirt_eos = MyFragment%nvirtEOS
    nocc_aos = MyFragment%noccAOS
    nvirt_aos = MyFragment%nvirtAOS

    ! **************************************************************
    ! ************** do energy for single fragment *****************
    ! **************************************************************

    ! ***********************
    !   do E[4] energy part
    ! ***********************

    ! init energy reals to be on the safe side
    ! note: OccEnergyPT and VirtEnergyPT is also initialized from in here
    !       as this (e4) routine is called before the e5 routine
    MyFragment%energies(FRAGMODEL_OCCpT) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_VIRTpT) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_OCCpT4) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_VIRTpT4) = 0.0E0_realk

    do_occ = .false.
    do_virt = .false.

    if ((.not. DECinfo%OnlyOccPart) .and. (.not. DECinfo%OnlyVirtPart)) then

       do_occ = .true.
       do_virt = .true.

    else if (DECinfo%OnlyOccPart .and. (.not. DECinfo%OnlyVirtPart)) then

       do_occ = .true.

    else if (DECinfo%OnlyVirtPart .and. (.not. DECinfo%OnlyOccPart)) then

       do_virt = .true.

    end if

    ! *******************************
    ! do occupied partitioning scheme
    ! *******************************

    if (do_occ) then

       energy_res_cou = 0.0E0_realk
       energy_res_exc = 0.0E0_realk
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,i_eos,j,a,b,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc_eos,nocc_aos,nvirt_aos,MyFragment),&
       !$OMP REDUCTION(+:energy_res_cou),REDUCTION(+:virt_contribs)
       do j=1,nocc_aos
          do i=1,nocc_eos
             i_eos = MyFragment%idxo(i)
   
             do b=1,nvirt_aos
                do a=1,nvirt_aos
   
                   energy_tmp = 4.0E0_realk * ccsd_doubles%elm4(a,b,i_eos,j) &
                                  & * ccsdpt_doubles%elm4(a,b,i_eos,j)
                   energy_res_cou = energy_res_cou + energy_tmp
   
                   ! update contribution from aos orbital a
                   virt_contribs(a) = virt_contribs(a) + energy_tmp
   
                   ! update contribution from aos orbital b 
                   ! (only if different from aos orbital a to avoid double counting)
                   if (a .ne. b) virt_contribs(b) = virt_contribs(b) + energy_tmp
   
                end do
             end do
   
          end do
       end do
       !$OMP END PARALLEL DO
   
       ! reorder from (a,b,i,j) to (a,b,j,i)
       call tensor_reorder(ccsd_doubles,[1,2,4,3])
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,i_eos,j,a,b,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc_eos,nocc_aos,nvirt_aos,MyFragment),&
       !$OMP REDUCTION(+:energy_res_exc),REDUCTION(+:virt_contribs)
       do j=1,nocc_aos
          do i=1,nocc_eos
             i_eos = MyFragment%idxo(i)
   
             do b=1,nvirt_aos
                do a=1,nvirt_aos
   
                   energy_tmp = 2.0E0_realk * ccsd_doubles%elm4(a,b,i_eos,j) &
                                  & * ccsdpt_doubles%elm4(a,b,i_eos,j)
                   energy_res_exc = energy_res_exc - energy_tmp
   
                   ! update contribution from aos orbital a
                   virt_contribs(a) = virt_contribs(a) - energy_tmp
   
                   ! update contribution from aos orbital b 
                   ! (only if different from aos orbital a to avoid double counting)
                   if (a .ne. b) virt_contribs(b) = virt_contribs(b) - energy_tmp
   
                end do
             end do
   
          end do
       end do
       !$OMP END PARALLEL DO
   
       !get total fourth--order energy contribution
       MyFragment%energies(FRAGMODEL_OCCpT4) = energy_res_cou + energy_res_exc
   
       ! insert into occ. part. scheme part
       MyFragment%energies(FRAGMODEL_OCCpT) = MyFragment%energies(FRAGMODEL_OCCpT) &
         & + MyFragment%energies(FRAGMODEL_OCCpT4)

    end if

    ! *********************************
    ! do virtupied partitioning scheme
    ! *********************************

    if (do_virt) then

       if (do_occ .and. do_virt) then

          ! initially, reorder ccsd_doubles and ccsdpt_doubles
          ! ccsd_doubles from from (a,b,j,i) sequence to (j,i,a,b) sequence
          ! ccsdpt_doubles from from (a,b,i,j) sequence to (i,j,a,b) sequence
          call tensor_reorder(ccsd_doubles,[3,4,1,2])
          call tensor_reorder(ccsdpt_doubles,[3,4,1,2])

       else

          ! initially, reorder ccsd_doubles and ccsdpt_doubles
          ! ccsd_doubles from from (a,b,i,j) sequence to (j,i,a,b) sequence
          ! ccsdpt_doubles from from (a,b,i,j) sequence to (i,j,a,b) sequence
          call tensor_reorder(ccsd_doubles,[4,3,1,2])
          call tensor_reorder(ccsdpt_doubles,[3,4,1,2])

       end if

       energy_res_cou = 0.0E0_realk
       energy_res_exc = 0.0E0_realk
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,a_eos,b,i,j,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nvirt_eos,nvirt_aos,nocc_aos,MyFragment),&
       !$OMP REDUCTION(+:energy_res_exc),REDUCTION(+:occ_contribs)
       do b=1,nvirt_aos
          do a=1,nvirt_eos
             a_eos = MyFragment%idxu(a)
   
             do j=1,nocc_aos
                do i=1,nocc_aos
   
                   energy_tmp = 2.0E0_realk * ccsd_doubles%elm4(i,j,a_eos,b) &
                                  & * ccsdpt_doubles%elm4(i,j,a_eos,b)
                   energy_res_exc = energy_res_exc - energy_tmp
   
                   ! update contribution from aos orbital i
                   occ_contribs(i) = occ_contribs(i) - energy_tmp
   
                   ! update contribution from aos orbital j 
                   ! (only if different from aos orbital i to avoid double counting)
                   if (i .ne. j) occ_contribs(j) = occ_contribs(j) - energy_tmp
   
                end do
             end do
   
          end do
       end do
       !$OMP END PARALLEL DO
   
       ! reorder form (j,i,a,b) to (i,j,a,b)
       call tensor_reorder(ccsd_doubles,[2,1,3,4])
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,a_eos,b,i,j,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nvirt_eos,nvirt_aos,nocc_aos,MyFragment),&
       !$OMP REDUCTION(+:energy_res_cou),REDUCTION(+:occ_contribs)
       do b=1,nvirt_aos
          do a=1,nvirt_eos
             a_eos = MyFragment%idxu(a)
   
             do j=1,nocc_aos
                do i=1,nocc_aos
   
                   energy_tmp = 4.0E0_realk * ccsd_doubles%elm4(i,j,a_eos,b) &
                                  & * ccsdpt_doubles%elm4(i,j,a_eos,b)
                   energy_res_cou = energy_res_cou + energy_tmp
   
                   ! update contribution from aos orbital i
                   occ_contribs(i) = occ_contribs(i) + energy_tmp
   
                   ! update contribution from aos orbital j 
                   ! (only if different from aos orbital i to avoid double counting)
                   if (i .ne. j) occ_contribs(j) = occ_contribs(j) + energy_tmp
   
                end do
             end do
   
          end do
       end do
       !$OMP END PARALLEL DO
   
       !get total fourth--order energy contribution
       MyFragment%energies(FRAGMODEL_VIRTpT4) = energy_res_cou + energy_res_exc
   
       ! insert into virt. part. scheme part
       MyFragment%energies(FRAGMODEL_VIRTpT) = MyFragment%energies(FRAGMODEL_VIRTpT) &
          & + MyFragment%energies(FRAGMODEL_VIRTpT4)

    end if

    ! ******************************
    !   done with E[4] energy part
    ! ******************************

    ! ************************************************************************
    !   as we need to reuse the ccsd doubles in the fragment optimization,
    !   we here reorder back into (a,i,b,j) sequence IF fragopt_pT == .true. 
    ! ************************************************************************

    if (present(fragopt_pT)) then

       if (do_occ .and. (.not. do_virt)) then

          ! reorder from (a,b,j,i) to (a,i,b,j)
          if (fragopt_pT) call tensor_reorder(ccsd_doubles,[1,4,2,3])

       else if (do_virt .and. (.not. do_virt)) then

          ! reorder from (i,j,a,b) to (a,i,b,j)
          if (fragopt_pT) call tensor_reorder(ccsd_doubles,[3,1,4,2])

       else if (do_occ .and. do_virt) then

          ! reorder from (i,j,a,b) to (a,i,b,j)
          if (fragopt_pT) call tensor_reorder(ccsd_doubles,[3,1,4,2])

       end if

    end if

    ! *******************************************************************
    ! ************** done w/ energy for single fragment *****************
    ! *******************************************************************

  end subroutine ccsdpt_decnp_e4_frag


  !> \brief: calculate E[4] contribution to pair fragment ccsd(t) energy correction
  !> \author: Janus Eriksen
  !> \date: september 2012
  subroutine ccsdpt_energy_e4_pair(Fragment1,Fragment2,PairFragment,ccsd_doubles,ccsdpt_doubles)

    implicit none

    !> fragment # 1 in the pair fragment
    type(decfrag),intent(in) :: Fragment1
    !> fragment # 2 in the pair fragment
    type(decfrag),intent(in) :: Fragment2
    !> pair fragment info
    type(decfrag), intent(inout) :: PairFragment
    ! ccsd and ccsd(t) doubles amplitudes
    type(tensor), intent(inout) :: ccsd_doubles, ccsdpt_doubles
    ! logical pointers for keeping hold of which pairs are to be handled
    logical, pointer :: dopair_occ(:,:), dopair_virt(:,:)
    !> integers
    integer :: nocc_eos, nocc_aos, nvirt_eos, nvirt_aos, i,j,a,b, i_eos, j_eos, a_eos, b_eos
    !> temporary energy arrays
    type(array2) :: energy_interm_cou, energy_interm_exc, energy_interm_ccsdpt
    !> energy reals
    real(realk) :: energy_tmp, energy_res_cou, energy_res_exc  
    !> which partitioning schemes?
    logical :: do_occ, do_virt

    ! init dimensions
    nocc_eos = PairFragment%noccEOS
    nvirt_eos = PairFragment%nvirtEOS
    nocc_aos = PairFragment%noccAOS
    nvirt_aos = PairFragment%nvirtAOS

    do_occ = .false.
    do_virt = .false.

    if ((.not. DECinfo%OnlyOccPart) .and. (.not. DECinfo%OnlyVirtPart)) then

       do_occ = .true.
       do_virt = .true.
       call mem_alloc(dopair_occ,nocc_eos,nocc_eos)
       call mem_alloc(dopair_virt,nvirt_eos,nvirt_eos)
       call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
       call which_pairs_virt(Fragment1,Fragment2,PairFragment,dopair_virt)

    else if (DECinfo%OnlyOccPart .and. (.not. DECinfo%OnlyVirtPart)) then

       do_occ = .true.
       call mem_alloc(dopair_occ,nocc_eos,nocc_eos)
       call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)

    else if (DECinfo%OnlyVirtPart .and. (.not. DECinfo%OnlyOccPart)) then

       do_virt = .true.
       call mem_alloc(dopair_virt,nvirt_eos,nvirt_eos)
       call which_pairs_virt(Fragment1,Fragment2,PairFragment,dopair_virt)

    end if

    ! *************************************************************
    ! ************** do energy for pair fragments *****************
    ! *************************************************************

    ! ***********************
    !   do E[4] energy part
    ! ***********************

    ! init energy reals to be on the safe side
    ! note: OccEnergyPT and VirtEnergyPT is also initialized from in here
    !       as this (e4) routine is called before the e5 routine
    PairFragment%energies(FRAGMODEL_OCCpT) = 0.0E0_realk
    PairFragment%energies(FRAGMODEL_VIRTpT) = 0.0E0_realk
    PairFragment%energies(FRAGMODEL_OCCpT4) = 0.0E0_realk
    PairFragment%energies(FRAGMODEL_VIRTpT4) = 0.0E0_realk

    ! *******************************
    ! do occupied partitioning scheme
    ! *******************************

    if (do_occ) then

       energy_res_cou = 0.0E0_realk
       energy_res_exc = 0.0E0_realk
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,i_eos,j,j_eos,a,b,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc_eos,nvirt_aos,&
       !$OMP PairFragment,dopair_occ),REDUCTION(+:energy_res_cou)
       do j=1,nocc_eos
          j_eos = PairFragment%idxo(j)
          do i=1,nocc_eos
             i_eos = PairFragment%idxo(i)
   
             if (.not. dopair_occ(i,j)) cycle 
   
             do b=1,nvirt_aos
                do a=1,nvirt_aos
   
                   energy_tmp = ccsd_doubles%elm4(a,b,i_eos,j_eos) &
                              & * ccsdpt_doubles%elm4(a,b,i_eos,j_eos)
                   energy_res_cou = energy_res_cou + energy_tmp
   
                end do
             end do
   
          end do
       end do
       !$OMP END PARALLEL DO
   
       ! reorder from (a,b,i,j) to (a,b,j,i)
       call tensor_reorder(ccsd_doubles,[1,2,4,3])
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,i_eos,j,j_eos,a,b,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc_eos,nvirt_aos,&
       !$OMP PairFragment,dopair_occ),REDUCTION(+:energy_res_exc)
       do j=1,nocc_eos
          j_eos = PairFragment%idxo(j)
          do i=1,nocc_eos
             i_eos = PairFragment%idxo(i)
   
             if (.not. dopair_occ(i,j)) cycle
   
             do b=1,nvirt_aos
                do a=1,nvirt_aos
   
                   energy_tmp = ccsd_doubles%elm4(a,b,i_eos,j_eos) &
                              & * ccsdpt_doubles%elm4(a,b,i_eos,j_eos)
                   energy_res_exc = energy_res_exc + energy_tmp
   
                end do
             end do
   
          end do
       end do
       !$OMP END PARALLEL DO
   
       ! get total fourth--order energy contribution
       PairFragment%energies(FRAGMODEL_OCCpT4) = 4.0E0_realk * energy_res_cou &
          & - 2.0E0_realk * energy_res_exc
   
       ! insert into occ. part. scheme part
       PairFragment%energies(FRAGMODEL_OCCpT) = PairFragment%energies(FRAGMODEL_OCCpT) &
          &+ PairFragment%energies(FRAGMODEL_OCCpT4)

    end if

    ! *********************************
    ! do virtupied partitioning scheme
    ! *********************************

    if (do_virt) then

       if (do_occ .and. do_virt) then

          ! initially, reorder ccsd_doubles and ccsdpt_doubles
          ! ccsd_doubles from from (a,b,j,i) sequence to (j,i,a,b) sequence
          ! ccsdpt_doubles from from (a,b,i,j) sequence to (i,j,a,b) sequence
          call tensor_reorder(ccsd_doubles,[3,4,1,2])
          call tensor_reorder(ccsdpt_doubles,[3,4,1,2])

       else

          ! initially, reorder ccsd_doubles and ccsdpt_doubles
          ! ccsd_doubles from from (a,b,i,j) sequence to (j,i,a,b) sequence
          ! ccsdpt_doubles from from (a,b,i,j) sequence to (i,j,a,b) sequence
          call tensor_reorder(ccsd_doubles,[4,3,1,2])
          call tensor_reorder(ccsdpt_doubles,[3,4,1,2])

       end if

       energy_res_cou = 0.0E0_realk
       energy_res_exc = 0.0E0_realk
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,a_eos,b,b_eos,i,j,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nvirt_eos,nocc_aos,&
       !$OMP PairFragment,dopair_virt),REDUCTION(+:energy_res_exc)
       do b=1,nvirt_eos
          b_eos = PairFragment%idxu(b)
          do a=1,nvirt_eos
             a_eos = PairFragment%idxu(a)
   
             if (.not. dopair_virt(a,b)) cycle
       
             do j=1,nocc_aos
                do i=1,nocc_aos
   
                   energy_tmp = ccsd_doubles%elm4(i,j,a_eos,b_eos) &
                              & * ccsdpt_doubles%elm4(i,j,a_eos,b_eos)
                   energy_res_exc = energy_res_exc + energy_tmp
   
                end do
             end do
       
          end do
       end do
       !$OMP END PARALLEL DO
   
       ! reorder form (j,i,a,b) to (i,j,a,b)
       call tensor_reorder(ccsd_doubles,[2,1,3,4])
   
       !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,a_eos,b,b_eos,i,j,energy_tmp),&
       !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nvirt_eos,nocc_aos,&
       !$OMP PairFragment,dopair_virt),REDUCTION(+:energy_res_cou)
       do b=1,nvirt_eos
          b_eos = PairFragment%idxu(b)
          do a=1,nvirt_eos
             a_eos = PairFragment%idxu(a)
   
             if (.not. dopair_virt(a,b)) cycle
   
             do j=1,nocc_aos
                do i=1,nocc_aos
   
                   energy_tmp = ccsd_doubles%elm4(i,j,a_eos,b_eos) &
                              & * ccsdpt_doubles%elm4(i,j,a_eos,b_eos)
                   energy_res_cou = energy_res_cou + energy_tmp
   
                end do
             end do
   
          end do
       end do
       !$OMP END PARALLEL DO
   
       ! get total fourth--order energy contribution
       PairFragment%energies(FRAGMODEL_VIRTpT4) = 4.0E0_realk * energy_res_cou &
          &- 2.0E0_realk * energy_res_exc
   
       ! insert into virt. part. scheme part
       PairFragment%energies(FRAGMODEL_VIRTpT) = PairFragment%energies(FRAGMODEL_VIRTpT) &
          & + PairFragment%energies(FRAGMODEL_VIRTpT4)

    end if

    ! ******************************
    !   done with E[4] energy part
    ! ******************************

    ! now release logical pair arrays
    if (do_occ .and. do_virt) then

      call mem_dealloc(dopair_occ)
      call mem_dealloc(dopair_virt)

    else if (do_occ .and. (.not. do_virt)) then

      call mem_dealloc(dopair_occ)

       else if (do_virt .and. (.not. do_occ)) then

      call mem_dealloc(dopair_virt)

    end if

    ! ******************************************************************
    ! ************** done w/ energy for pair fragments *****************
    ! ******************************************************************

  end subroutine ccsdpt_energy_e4_pair


  !> \brief Get MO integrals for ijk-CCSD(T) (in canonical basis), see integral storing order below.
  !> \author Janus Eriksen and Kasper Kristensen
  !> \date September-October 2012
  subroutine get_CCSDpT_integrals_ijk(MyLsitem,nbasis,nocc,nvirt,Cocc,Cvirt,ovoo,vvvo)

    implicit none

    !> Integral info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of basis functions
    integer,intent(in) :: nbasis
    !> Number of occupied orbitals
    integer,intent(in) :: nocc
    !> Number of virtual orbitals
    integer,intent(in) :: nvirt
    !> Occupied MO coefficients
    real(realk), dimension(nbasis,nocc),intent(in) :: Cocc
    !> Virtual MO coefficients
    real(realk), dimension(nbasis,nvirt),intent(in) :: Cvirt
    ! ovoo: Integrals (AI|JK) in the order (J,A,I,K)
    type(tensor), intent(inout) :: ovoo
    ! vvvo: Integrals (AI|BC) in the order (C,B,A,I)
    type(tensor), intent(inout) :: vvvo
    integer :: gammadim, alphadim,iorb
    integer :: alphaB,gammaB,dimAlpha,dimGamma,idx
    real(realk),pointer :: tmp1(:),tmp2(:),tmp3(:)
    integer(kind=long) :: size1,size2,size3
    integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd,m,k,n,i,j,dims(4),order(4)
    logical :: FullRHS,doscreen
    real(realk) :: tcpu, twall
    real(realk),pointer :: CoccT(:,:), CvirtT(:,:)
    integer :: MaxActualDimAlpha,nbatchesAlpha
    integer :: MaxActualDimGamma,nbatchesGamma
    type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
    type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
    integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint
    logical :: MoTrans, NoSymmetry,SameMol
    integer, pointer :: orb2batchAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchsizeGamma(:), batchindexGamma(:)
    TYPE(DECscreenITEM)   :: DecScreen
    type(batchtoorb), pointer :: batch2orbAlpha(:)
    type(batchtoorb), pointer :: batch2orbGamma(:)
    integer, pointer :: batchdimAlpha(:), batchdimGamma(:)
    ! distribution stuff needed for mpi parallelization
    integer, pointer :: distribution(:)
    Character            :: intSpec(5)
    integer :: myload,first_el_i_block
    logical :: master
    integer(kind=long) :: o3v,v3
    real(realk), pointer :: dummy2(:)
    integer(kind=ls_mpik) :: mode,dest,nel2t, wi_idx
    integer :: p,pos
    call time_start_phase(PHASE_WORK)

    o3v           = nocc*nocc*nocc*nvirt
    v3            = nvirt**3

#ifdef VAR_MPI

    master = (infpar%lg_mynum .eq. infpar%master)
    if (master) call LSTIMER('START',tcpu,twall,DECinfo%output)

#else

    master = .true.
    call LSTIMER('START',tcpu,twall,DECinfo%output)

#endif

    ! Set integral info
    ! *****************
    INTSPEC(1)='R' !R = Regular Basis set on the 1th center 
    INTSPEC(2)='R' !R = Regular Basis set on the 2th center 
    INTSPEC(3)='R' !R = Regular Basis set on the 3th center 
    INTSPEC(4)='R' !R = Regular Basis set on the 4th center 
    INTSPEC(5)='C' !C = Coulomb operator
    if (DECinfo%useichor) then
       iprint = 0           !print level for Ichor Integral code
       MoTrans = .FALSE.    !Do not transform to MO basis! 
       NoSymmetry = .FALSE. !Use Permutational Symmetry! 
       SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
       !Determine the full number of AO batches - not to be confused with the batches of AOs
       !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
       iAO = 1
       call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
    else
       ! Integral screening?
       doscreen = mylsitem%setting%scheme%cs_screen .or. mylsitem%setting%scheme%ps_screen
    endif
    ! allocate arrays to update during integral loop 
    ! **********************************************
    
    ! note 1: this must be done before call to get_optimal_batch_sizes_ccsdpt_integrals

    ! Integrals (AI|KJ) in the order (J,A,I,K)
    dims = [nocc,nvirt,nocc,nocc]
    call tensor_init(ovoo, dims,4)
    call tensor_zero(ovoo)

    ! Integrals (AB|IC) in the order (C,B,A,I)
    dims = [nvirt,nvirt,nvirt,nocc]
#ifdef VAR_MPI
    mode   = MPI_MODE_NOCHECK

    call tensor_ainit(vvvo,dims,4,tdims=[nvirt,nvirt,nvirt,1],atype="TDAR")
    call tensor_zero_tiled_dist(vvvo)

#else

    call tensor_init(vvvo, dims,4)
    call tensor_zero(vvvo)

#endif

    ! For efficiency when calling dgemm, save transposed matrices
    call mem_alloc(CoccT,nocc,nbasis)
    call mem_alloc(CvirtT,nvirt,nbasis)
    call mat_transpose(nbasis,nocc,1.0E0_realk,Cocc,0.0E0_realk,CoccT)
    call mat_transpose(nbasis,nvirt,1.0E0_realk,Cvirt,0.0E0_realk,CvirtT)

    ! Determine optimal batchsizes and corresponding sizes of arrays
    call get_optimal_batch_sizes_ccsdpt_integrals(mylsitem,nbasis,nocc,nvirt,alphadim,gammadim,&
         & size1,size2,size3,.true.,.false.,1)


    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************
    if (DECinfo%useichor) then
       iAO = 4 !Gamma is the 4. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the gammadim, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',gammadim,&
            & nbatchesGamma,DECinfo%output)
       call mem_alloc(AOGammabatchinfo,nbatchesGamma)
       !Construct the batches of AOS based on the gammadim, the requested
       !size of the AO batches - gammadim must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimGamma must be less og equal to gammadim
       call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',gammadim,&
            & nbatchesGamma,AOGammabatchinfo,MaxActualDimGamma,DECinfo%output)
    else
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchGamma,nbasis)
       call build_batchesofAOS(DECinfo%output,mylsitem%setting,gammadim,&
            & nbasis,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,&
            & nbatchesGamma,orb2BatchGamma,'R')
    endif
    if(master.and.DECinfo%PL>1)write(*,*) 'BATCH: Number of Gamma batches   = ', nbatchesGamma

    if (.not. DECinfo%useichor) then
       ! Translate batchindex to orbital index
       ! -------------------------------------
       call mem_alloc(batch2orbGamma,nbatchesGamma)
   
       do idx=1,nbatchesGamma
   
          call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx) )
          batch2orbGamma(idx)%orbindex = 0
          batch2orbGamma(idx)%norbindex = 0
   
       end do
   
       do iorb=1,nbasis
   
          idx = orb2batchGamma(iorb)
          batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
          K = batch2orbGamma(idx)%norbindex
          batch2orbGamma(idx)%orbindex(K) = iorb
   
       end do
    endif

    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************

    if (DECinfo%useichor) then
       iAO = 3 !Alpha is the 3. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the alphadim, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',alphadim,&
            & nbatchesAlpha,DECinfo%output)
       call mem_alloc(AOAlphabatchinfo,nbatchesAlpha)
       !Construct the batches of AOS based on the alphadim, the requested
       !size of the AO batches - alphadim must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimAlpha must be less og equal to alphadim
       call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',alphadim,&
            & nbatchesAlpha,AOAlphabatchinfo,MaxActualDimAlpha,DECinfo%output)
    else
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchAlpha,nbasis)
       call build_batchesofAOS(DECinfo%output,mylsitem%setting,alphadim,&
            & nbasis,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,&
            & nbatchesAlpha,orb2BatchAlpha,'R')
    endif

    if(master.and.DECinfo%PL>1)write(*,*) 'BATCH: Number of Alpha batches   = ', nbatchesAlpha

    if (.not. DECinfo%useichor) then
       ! Translate batchindex to orbital index
       ! -------------------------------------
       call mem_alloc(batch2orbAlpha,nbatchesAlpha)
   
       do idx=1,nbatchesAlpha
   
          call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
          batch2orbAlpha(idx)%orbindex = 0
          batch2orbAlpha(idx)%norbindex = 0
   
       end do
   
       do iorb=1,nbasis
   
          idx = orb2batchAlpha(iorb)
          batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
          K = batch2orbAlpha(idx)%norbindex
          batch2orbAlpha(idx)%orbindex(K) = iorb
   
       end do
    endif


    if (DECinfo%useichor) then
       !Calculate Screening integrals 
       SameMOL = .TRUE. !Specifies same molecule on all centers 
       call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
    else
       call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mylsitem%setting,&
               & nbatchesAlpha,nbatchesGamma,INTSPEC,DECinfo%IntegralThreshold)
       if (doscreen) then
   
          call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
               & nbasis,nbatchesAlpha,nbatchesGamma,&
               & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
               & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
   
       end if
    endif

    FullRHS = (nbatchesGamma .eq. 1) .and. (nbatchesAlpha .eq. 1)


    ! Allocate array for AO integrals
    ! *******************************
    call mem_alloc(tmp1,size1)
    call mem_alloc(tmp2,size2)
    call mem_alloc(tmp3,size3)

#ifdef VAR_MPI

    ! alloc distribution array
    nullify(distribution)
    call mem_alloc(distribution,nbatchesGamma*nbatchesAlpha)

    ! init distribution
    distribution = 0
    myload = 0
    if (DECinfo%useichor) then
       call mem_alloc(batchdimAlpha,nbatchesAlpha)
       do idx=1,nbatchesAlpha
          batchdimAlpha(idx) = AOAlphabatchinfo(idx)%dim 
       enddo
       call mem_alloc(batchdimGamma,nbatchesGamma)
       do idx=1,nbatchesGamma
          batchdimGamma(idx) = AOGammabatchinfo(idx)%dim 
       enddo
    endif
    call distribute_mpi_jobs(distribution,nbatchesAlpha,nbatchesGamma,&
    &batchdimAlpha,batchdimGamma,myload,infpar%lg_nodtot,infpar%lg_mynum)
    if (DECinfo%useichor) then
       call mem_dealloc(batchdimAlpha)
       call mem_dealloc(batchdimGamma)
    endif

#endif

    ! Start looping over gamma and alpha batches and calculate integrals
    ! ******************************************************************

    BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
       if (DECinfo%useichor) then
          dimGamma = AOGammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
          GammaStart = AOGammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
          GammaEnd = AOGammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
          AOGammaStart = AOGammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch
          AOGammaEnd = AOGammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch
       else
          dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
          GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
          GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch
       endif

       BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
          if (DECinfo%useichor) then
             dimAlpha = AOAlphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
             AlphaStart = AOAlphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
             AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
             AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
             AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
          else
             dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
             AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
             AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch
          endif

#ifdef VAR_MPI

          ! distribute tasks
          if (distribution((alphaB-1)*nbatchesGamma+gammaB) .ne. infpar%lg_mynum) then

             cycle BatchAlpha

          end if

          if(DECinfo%PL>2)write (*, '("Rank(T) ",I3," starting job (",I3,"/",I3,",",I3,"/",I3,")")')&
             &infpar%lg_mynum,alphaB,nbatchesAlpha,gammaB,nbatchesGamma

#endif

          ! Get (beta delta | alphaB gammaB) integrals using (beta,delta,alphaB,gammaB) ordering
          ! ************************************************************************************

          if (DECinfo%useichor) then
             call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,nbasis,nbasis,dimAlpha,dimGamma,&
                  & tmp1,INTSPEC,FULLRHS,1,nAObatches,1,nAObatches,AOAlphaStart,AOAlphaEnd,&
                  & AOGammaStart,AOGammaEnd,MoTrans,nbasis,nbasis,dimAlpha,dimGamma,NoSymmetry,DECinfo%IntegralThreshold)
          else
             if (doscreen) mylsitem%setting%LST_GAB_LHS => DECSCREEN%masterGabLHS
             if (doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGab(alphaB,gammaB)%p
   
   
             call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, &
                  & mylsitem%setting,tmp1,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
                  & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),nbasis,nbasis,dimAlpha,dimGamma,&
                  & FullRHS,INTSPEC,DECinfo%IntegralThreshold)
          endif

          ! tmp2(delta,alphaB,gammaB;A) = sum_{beta} [tmp1(beta;delta,alphaB,gammaB)]^T Cvirt(beta,A)
          m = nbasis*dimGamma*dimAlpha
          k = nbasis
          n = nvirt
          call dgemm('T','N',m,n,k,1.0E0_realk,tmp1,k,Cvirt,k,0.0E0_realk,tmp2,m)

          ! tmp3(B;alphaB,gammaB,A) = sum_{delta} CvirtT(B,delta) tmp2(delta;alphaB,gammaB,A)
          m = nvirt
          k = nbasis
          n = dimAlpha*dimGamma*nvirt
          call dgemm('N','N',m,n,k,1.0E0_realk,CvirtT,m,tmp2,k,0.0E0_realk,tmp3,m)

          ! tmp1(I;,alphaB,gammaB,A) = sum_{delta} CoccT(I,delta) tmp2(delta,alphaB,gammaB,A)
          m = nocc
          k = nbasis
          n = dimAlpha*dimGamma*nvirt
          call dgemm('N','N',m,n,k,1.0E0_realk,CoccT,m,tmp2,k,0.0E0_realk,tmp1,m)

          ! Reorder: tmp1(I,alphaB;gammaB,A) --> tmp2(gammaB,A;I,alphaB)
          m = nocc*dimAlpha
          n = dimGamma*nvirt
          call mat_transpose(m,n,1.0E0_realk,tmp1,0.0E0_realk,tmp2)

          ! tmp1(J;A,I,alphaB) = sum_{gamma in gammaB} CoccT(J,gamma) tmp2(gamma,A,I,alphaB)
          m = nocc
          k = dimGamma
          n = nvirt*nocc*dimAlpha
          call dgemm('N','N',m,n,k,1.0E0_realk,CoccT(1,GammaStart),nocc,tmp2,k,0.0E0_realk,tmp1,m)

          ! ovoo(J,A,I;K) += sum_{alpha in alphaB} tmp1(J,A,I,alpha) Cocc(alpha,K)
          m = nvirt*nocc**2
          k = dimAlpha
          n = nocc
          call dgemm('N','N',m,n,k,1.0E0_realk,tmp1,m,Cocc(AlphaStart,1),nbasis,1.0E0_realk,ovoo%elm1,m)

          ! Reorder: tmp3(B,alphaB;gammaB,A) --> tmp1(gammaB,A;B,alphaB)
          m = nvirt*dimAlpha
          n = dimGamma*nvirt
          call mat_transpose(m,n,1.0E0_realk,tmp3,0.0E0_realk,tmp1)

          ! tmp3(C;A,B,alphaB) = sum_{gamma in gammaB} CvirtT(C,gamma) tmp1(gamma,A,B,alphaB)
          m = nvirt
          k = dimGamma
          n = dimAlpha*nvirt**2
          call dgemm('N','N',m,n,k,1.0E0_realk,CvirtT(1,GammaStart),nvirt,tmp1,k,0.0E0_realk,tmp3,m)

          m = nvirt**3
          k = dimAlpha
          n = 1
#ifdef VAR_MPI
          ! reorder tmp1 and do vvvo(B,A,C,I) += sum_{i in IB} tmp1(B,A,C,i)
          do i=1,nocc

             ! tmp1(C,A,B,i) = sum_{alpha in alphaB} tmp3(C,A,B,alpha) Cocc(alpha,i)
             !call dgemm('N','N',m,n,k,1.0E0_realk,tmp3,m,Cocc(AlphaStart,i),nbasis,0.0E0_realk,tmp1,m)
             call dgemv('N',m,k,1.0E0_realk,tmp3,m,Cocc(AlphaStart,i),1,0.0E0_realk,tmp1,1)

             ! *** tmp1 corresponds to (AB|iC) in Mulliken notation. Noting that the v³o integrals
             ! are normally written as g_{AIBC}, we may also write this Mulliken integral (with substitution
             ! of dummy indices A=B, B=C, and C=A) as (BC|IA). In order to align with the vvvo order of
             ! ccsd(t) driver routine, we reorder as:
             ! (BC|IA) --> (CB|AI), i.e., tmp1(C,A,B,i) = ABCI(A,B,C,i) (norm. notat.) --> 
             !                                            tmp1(C,B,A,i) (norm. notat.) = tmp1(B,A,C,i) (notat. herein)
             ! 
             ! next, we accumulate
             ! vvvo(B,A,C,I) += sum_{i in IB} tmp1(B,A,C,i)

             call array_reorder_3d(1.0E0_realk,tmp1,nvirt,nvirt,nvirt,[3,2,1],0.0E0_realk,tmp2)

             call time_start_phase(PHASE_COMM)

#ifdef VAR_HAVE_MPI3
             call tensor_lock_win(vvvo,i,'s')
#endif
             !call tensor_accumulate_tile(vvvo,i,tmp2,nvirt**3,lock_set=.true.,flush_it=.true.)

             call get_residence_of_tile(dest,i,vvvo,idx_on_node = pos)

             if( alloc_in_dummy )then
                wi_idx = 1
                p      = pos - 1
             else
                wi_idx = i 
                p      = 0
             endif

             do first_el_i_block=1,v3,MAX_SIZE_ONE_SIDED
#ifndef VAR_HAVE_MPI3
                call tensor_lock_win(vvvo,i,'s',assert=mode)
#endif
                nel2t=MAX_SIZE_ONE_SIDED
                if(((v3-first_el_i_block)<MAX_SIZE_ONE_SIDED).and.&
                   &(mod(v3-first_el_i_block+1,i8*MAX_SIZE_ONE_SIDED)/=0))&
                   &nel2t=int(mod(v3,i8*MAX_SIZE_ONE_SIDED),kind=ls_mpik)


                call lsmpi_acc(tmp2(first_el_i_block:first_el_i_block+nel2t-1),nel2t,p+first_el_i_block,dest,vvvo%wi(wi_idx))

#ifdef VAR_HAVE_MPI3
                call lsmpi_win_flush(vvvo%wi(wi_idx),rank=dest,local=.true.)
#else
                call tensor_unlock_win(vvvo,i)
#endif
             enddo

#ifdef VAR_HAVE_MPI3
             call tensor_unlock_win(vvvo,i)
#endif
             call time_start_phase(PHASE_WORK)

          end do

#else

          do i=1,nocc

             ! for description, see mpi section above
             !call dgemm('N','N',m,n,k,1.0E0_realk,tmp3,m,Cocc(AlphaStart,i),nbasis,0.0E0_realk,tmp1,m)
             call dgemv('N',m,k,1.0E0_realk,tmp3,m,Cocc(AlphaStart,i),1,0.0E0_realk,tmp1,1)

             call array_reorder_3d(1.0E0_realk,tmp1,nvirt,nvirt,nvirt,[3,2,1],1.0E0_realk,vvvo%elm4(:,:,:,i))

          end do

#endif

       end do BatchAlpha
    end do BatchGamma

#ifdef VAR_MPI

    if (infpar%lg_nodtot .gt. 1) then

#ifdef VAR_PTR_RESHAPE
       dummy2(1:(i8*nocc*nvirt)*nocc*nocc) => ovoo%elm1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
       call c_f_pointer(c_loc(ovoo%elm1(1)),dummy2,[(i8*nocc*nvirt)*nocc*nocc])      
#else
       call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
       
       call time_start_phase(PHASE_IDLE)
       call lsmpi_barrier(infpar%lg_comm)

       ! now, reduce o^3v integrals onto master
       call time_start_phase(PHASE_COMM)
       call lsmpi_allreduce(dummy2,o3v,infpar%lg_comm) 
       call time_start_phase(PHASE_WORK)

    end if

    ! dealloc distribution array
    call mem_dealloc(distribution)

#endif

    ! free stuff
    ! **********
    call mem_dealloc(tmp1)
    call mem_dealloc(tmp2)
    call mem_dealloc(tmp3)
    call mem_dealloc(CoccT)
    call mem_dealloc(CvirtT)
    if (DECinfo%useichor) then
       call FREE_SCREEN_ICHORERI()
       call mem_dealloc(AOGammabatchinfo)
       call mem_dealloc(AOAlphabatchinfo)
    else
       call free_decscreen(DECSCREEN)
       call mem_dealloc(orb2batchGamma)
       call mem_dealloc(batchdimGamma)
       call mem_dealloc(batchsizeGamma)
       call mem_dealloc(batchindexGamma)
       do idx=1,nbatchesGamma
          call mem_dealloc(batch2orbGamma(idx)%orbindex)
       end do
       call mem_dealloc(batch2orbGamma)
       call mem_dealloc(orb2batchAlpha)
       call mem_dealloc(batchdimAlpha)
       call mem_dealloc(batchsizeAlpha)
       call mem_dealloc(batchindexAlpha)
       do idx=1,nbatchesAlpha
          call mem_dealloc(batch2orbAlpha(idx)%orbindex)
          batch2orbAlpha(idx)%orbindex => null()
       end do
       call mem_dealloc(batch2orbAlpha)
       nullify(mylsitem%setting%LST_GAB_LHS)
       nullify(mylsitem%setting%LST_GAB_RHS)
    endif
    if (master) call LSTIMER('CCSD(T) INT (IJK)',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine get_CCSDpT_integrals_ijk


  !> \brief Get MO integrals for abc-CCSD(T) (in canonical basis), see integral storing order below.
  !> \author Janus Eriksen and Kasper Kristensen
  !> \date September 2014
  subroutine get_CCSDpT_integrals_abc(MyLsitem,nbasis,nocc,nvirt,Cocc,Cvirt,ooov,vovv,tile_size)

    implicit none

    !> Integral info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of basis functions
    integer,intent(in) :: nbasis
    !> Number of occupied orbitals
    integer,intent(in) :: nocc
    !> Number of virtual orbitals
    integer,intent(in) :: nvirt
    !> Occupied MO coefficients
    real(realk), dimension(nbasis,nocc),intent(in) :: Cocc
    !> Virtual MO coefficients
    real(realk), dimension(nbasis,nvirt),intent(in) :: Cvirt
    ! Integrals (AI|JK) in the order (I,J,K,A)
    type(tensor), intent(inout) :: ooov
    ! Integrals (AI|BC) in the order (B,I,A,C)
    type(tensor), intent(inout) :: vovv
    integer, intent(inout) :: tile_size
    integer :: gammadim, alphadim,iorb,tile
    integer :: alphaB,gammaB,dimAlpha,dimGamma,idx
    real(realk),pointer :: tmp1(:),tmp2(:),tmp3(:)
    integer(kind=long) :: size1,size2,size3
    integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd,m,k,n,i,j,c,dims(4),order(4)
    logical :: FullRHS,doscreen
    real(realk) :: tcpu, twall, vovv_norm
    real(realk),pointer :: CoccT(:,:), CvirtT(:,:)
    integer :: MaxActualDimAlpha,nbatchesAlpha
    integer :: MaxActualDimGamma,nbatchesGamma
    type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
    type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
    integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint
    logical :: MoTrans, NoSymmetry,SameMol
    integer, pointer :: orb2batchAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchsizeGamma(:), batchindexGamma(:)
    TYPE(DECscreenITEM)   :: DecScreen
    type(batchtoorb), pointer :: batch2orbAlpha(:)
    type(batchtoorb), pointer :: batch2orbGamma(:)
    integer, pointer :: batchdimAlpha(:), batchdimGamma(:)

    ! distribution stuff needed for mpi parallelization
    integer, pointer :: distribution(:)
    Character            :: intSpec(5)
    integer :: myload,first_el_c_block,nelms,tile_size_tmp,total_num_tiles
    logical :: master
    integer(kind=long) :: o3v,v3,ov2
    real(realk), pointer :: dummy2(:)
    integer(kind=ls_mpik) :: mode,dest,nel2t, wi_idx
    integer :: p,pos
    call time_start_phase(PHASE_WORK)

    o3v           = nocc*nocc*nocc*nvirt
    v3            = nvirt**3
    ov2           = nocc*nvirt**2

#ifdef VAR_MPI

    master = (infpar%lg_mynum .eq. infpar%master)
    if (master) call LSTIMER('START',tcpu,twall,DECinfo%output)

#else

    master = .true.
    call LSTIMER('START',tcpu,twall,DECinfo%output)

#endif

    ! Set integral info
    ! *****************
    INTSPEC(1)='R' !R = Regular Basis set on the 1th center 
    INTSPEC(2)='R' !R = Regular Basis set on the 2th center 
    INTSPEC(3)='R' !R = Regular Basis set on the 3th center 
    INTSPEC(4)='R' !R = Regular Basis set on the 4th center 
    INTSPEC(5)='C' !C = Coulomb operator
    if (DECinfo%useichor) then
       iprint = 0           !print level for Ichor Integral code
       MoTrans = .FALSE.    !Do not transform to MO basis! 
       NoSymmetry = .FALSE. !Use Permutational Symmetry! 
       SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
       !Determine the full number of AO batches - not to be confused with the batches of AOs
       !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched
       iAO = 1
       call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
    else
       ! Integral screening?
       doscreen = mylsitem%setting%scheme%cs_screen .or. mylsitem%setting%scheme%ps_screen
    endif

    ! allocate arrays to update during integral loop 
    ! **********************************************
    
    ! note 1: this must be done before call to get_optimal_batch_sizes_ccsdpt_integrals

    ! ooov: Integrals (AI|KJ) in the order (I,J,K,A)
    dims = [nocc,nocc,nocc,nvirt]
    call tensor_init(ooov, dims,4)
    call tensor_zero(ooov)

    ! vovv: Integrals (AB|IC) in the order (B,I,A,C)
    dims = [nvirt,nocc,nvirt,nvirt]

#ifdef VAR_MPI

    mode   = MPI_MODE_NOCHECK

    call tensor_ainit(vovv,dims,4,tdims=[nvirt,nocc,nvirt,tile_size],atype="TDAR")
    call tensor_zero_tiled_dist(vovv)

#else

    call tensor_init(vovv, dims,4)
    call tensor_zero(vovv)

#endif

    ! For efficiency when calling dgemm, save transposed matrices
    call mem_alloc(CoccT,nocc,nbasis)
    call mem_alloc(CvirtT,nvirt,nbasis)
    call mat_transpose(nbasis,nocc,1.0E0_realk,Cocc,0.0E0_realk,CoccT)
    call mat_transpose(nbasis,nvirt,1.0E0_realk,Cvirt,0.0E0_realk,CvirtT)

    ! Determine optimal batchsizes and corresponding sizes of arrays
    call get_optimal_batch_sizes_ccsdpt_integrals(mylsitem,nbasis,nocc,nvirt,alphadim,gammadim,&
         & size1,size2,size3,.true.,.true.,tile_size)


    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************
    if (DECinfo%useichor) then
       iAO = 4 !Gamma is the 4. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the gammadim, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',gammadim,&
            & nbatchesGamma,DECinfo%output)
       call mem_alloc(AOGammabatchinfo,nbatchesGamma)
       !Construct the batches of AOS based on the gammadim, the requested
       !size of the AO batches - gammadim must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimGamma must be less og equal to gammadim
       call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',gammadim,&
            & nbatchesGamma,AOGammabatchinfo,MaxActualDimGamma,DECinfo%output)
    else
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchGamma,nbasis)
       call build_batchesofAOS(DECinfo%output,mylsitem%setting,gammadim,&
            & nbasis,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,&
            & nbatchesGamma,orb2BatchGamma,'R')
     endif

    if(master.and.DECinfo%PL>1)write(*,*) 'BATCH: Number of Gamma batches   = ', nbatchesGamma

    if (.not. DECinfo%useichor) then
       ! Translate batchindex to orbital index
       ! -------------------------------------
       call mem_alloc(batch2orbGamma,nbatchesGamma)
   
       do idx=1,nbatchesGamma
   
          call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx) )
          batch2orbGamma(idx)%orbindex = 0
          batch2orbGamma(idx)%norbindex = 0
   
       end do
   
       do iorb=1,nbasis
   
          idx = orb2batchGamma(iorb)
          batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
          K = batch2orbGamma(idx)%norbindex
          batch2orbGamma(idx)%orbindex(K) = iorb
   
       end do
    endif

    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************
    if (DECinfo%useichor) then
       iAO = 3 !Alpha is the 3. Center of the 4 center two electron coulomb integral
       !Determine how many batches of AOS based on the alphadim, the requested
       !size of the AO batches. iAO is the center that the batching should occur on. 
       !'R'  !Specifies that it is the Regular AO basis that should be batched 
       call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',alphadim,&
            & nbatchesAlpha,DECinfo%output)
       call mem_alloc(AOAlphabatchinfo,nbatchesAlpha)
       !Construct the batches of AOS based on the alphadim, the requested
       !size of the AO batches - alphadim must be unchanged since the call 
       !to determine_Ichor_nbatchesofAOS
       !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
       !So MaxActualDimAlpha must be less og equal to alphadim
       call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',alphadim,&
            & nbatchesAlpha,AOAlphabatchinfo,MaxActualDimAlpha,DECinfo%output)
    else
       ! Orbital to batch information
       ! ----------------------------
       call mem_alloc(orb2batchAlpha,nbasis)
       call build_batchesofAOS(DECinfo%output,mylsitem%setting,alphadim,&
            & nbasis,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,&
            & nbatchesAlpha,orb2BatchAlpha,'R')
     endif

    if(master.and.DECinfo%PL>1)write(*,*) 'BATCH: Number of Alpha batches   = ', nbatchesAlpha

    if (.not. DECinfo%useichor) then
       ! Translate batchindex to orbital index
       ! -------------------------------------
       call mem_alloc(batch2orbAlpha,nbatchesAlpha)
   
       do idx=1,nbatchesAlpha
   
          call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
          batch2orbAlpha(idx)%orbindex = 0
          batch2orbAlpha(idx)%norbindex = 0
   
       end do
   
       do iorb=1,nbasis
   
          idx = orb2batchAlpha(iorb)
          batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
          K = batch2orbAlpha(idx)%norbindex
          batch2orbAlpha(idx)%orbindex(K) = iorb
   
       end do
    endif

    if (DECinfo%useichor) then
       !Calculate Screening integrals 
       SameMOL = .TRUE. !Specifies same molecule on all centers 
       call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
    else
       call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mylsitem%setting,&
               & nbatchesAlpha,nbatchesGamma,INTSPEC,DECinfo%IntegralThreshold)

       if (doscreen) then
   
          call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
               & nbasis,nbatchesAlpha,nbatchesGamma,&
               & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
               & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
   
       end if
    endif

    FullRHS = (nbatchesGamma .eq. 1) .and. (nbatchesAlpha .eq. 1)


    ! Allocate array for AO integrals
    ! *******************************
    call mem_alloc(tmp1,size1)
    call mem_alloc(tmp2,size2)
    call mem_alloc(tmp3,size3)

#ifdef VAR_MPI

    ! alloc distribution array
    nullify(distribution)
    call mem_alloc(distribution,nbatchesGamma*nbatchesAlpha)

    ! init distribution
    distribution = 0
    myload = 0
    if (DECinfo%useichor) then
       call mem_alloc(batchdimAlpha,nbatchesAlpha)
       do idx=1,nbatchesAlpha
          batchdimAlpha(idx) = AOAlphabatchinfo(idx)%dim 
       enddo
       call mem_alloc(batchdimGamma,nbatchesGamma)
       do idx=1,nbatchesGamma
          batchdimGamma(idx) = AOGammabatchinfo(idx)%dim 
       enddo
    endif
    call distribute_mpi_jobs(distribution,nbatchesAlpha,nbatchesGamma,&
    &batchdimAlpha,batchdimGamma,myload,infpar%lg_nodtot,infpar%lg_mynum)
    if (DECinfo%useichor) then
       call mem_dealloc(batchdimAlpha)
       call mem_dealloc(batchdimGamma)
    endif

#endif

    ! Start looping over gamma and alpha batches and calculate integrals
    ! ******************************************************************

    BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
       if (DECinfo%useichor) then
          dimGamma = AOGammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
          GammaStart = AOGammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
          GammaEnd = AOGammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
          AOGammaStart = AOGammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch
          AOGammaEnd = AOGammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch
       else
          dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
          GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
          GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch
       endif

       BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
          if (DECinfo%useichor) then
             dimAlpha = AOAlphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
             AlphaStart = AOAlphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
             AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
             AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
             AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
          else
             dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
             AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
             AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch
          endif

#ifdef VAR_MPI

          ! distribute tasks
          if (distribution((alphaB-1)*nbatchesGamma+gammaB) .ne. infpar%lg_mynum) then

             cycle BatchAlpha

          end if

          if(DECinfo%PL>2)write (*, '("Rank(T) ",I3," starting job (",I3,"/",I3,",",I3,"/",I3,")")')&
             &infpar%lg_mynum,alphaB,nbatchesAlpha,gammaB,nbatchesGamma

#endif

          if (DECinfo%useichor) then
             call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,nbasis,nbasis,dimAlpha,dimGamma,&
                  & tmp1,INTSPEC,FULLRHS,1,nAObatches,1,nAObatches,AOAlphaStart,AOAlphaEnd,&
                  & AOGammaStart,AOGammaEnd,MoTrans,nbasis,nbasis,dimAlpha,dimGamma,NoSymmetry,&
                  & DECinfo%IntegralThreshold)
          else
             if (doscreen) mylsitem%setting%LST_GAB_LHS => DECSCREEN%masterGabLHS
             if (doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGab(alphaB,gammaB)%p
   
   
             ! Get (beta delta | alphaB gammaB) integrals using (beta,delta,alphaB,gammaB) ordering
             ! ************************************************************************************
             call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, &
                  & mylsitem%setting,tmp1,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
                  & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),nbasis,nbasis,dimAlpha,dimGamma,&
                  & FullRHS,INTSPEC,DECinfo%IntegralThreshold)
          endif
          ! tmp2(delta,alphaB,gammaB;I) = sum_{beta} [tmp1(beta;delta,alphaB,gammaB)]^T Cocc(beta,I)
          m = nbasis*dimGamma*dimAlpha
          k = nbasis
          n = nocc
          call dgemm('T','N',m,n,k,1.0E0_realk,tmp1,k,Cocc,k,0.0E0_realk,tmp2,m)

          ! tmp3(J;alphaB,gammaB,I) = sum_{delta} CoccT(J,delta) tmp2(delta;alphaB,gammaB,I)
          m = nocc
          k = nbasis
          n = dimAlpha*dimGamma*nocc
          call dgemm('N','N',m,n,k,1.0E0_realk,CoccT,m,tmp2,k,0.0E0_realk,tmp3,m)

          ! tmp1(A;alphaB,gammaB,I) = sum_{delta} CvirtT(A,delta) tmp2(delta,alphaB,gammaB,I)
          m = nvirt
          k = nbasis
          n = dimAlpha*dimGamma*nocc
          call dgemm('N','N',m,n,k,1.0E0_realk,CvirtT,m,tmp2,k,0.0E0_realk,tmp1,m)

          ! Reorder: tmp3(J,alphaB;gammaB,I) --> tmp2(gammaB,I;J,alphaB)
          m = nocc*dimAlpha
          n = dimGamma*nocc
          call mat_transpose(m,n,1.0E0_realk,tmp3,0.0E0_realk,tmp2)

          ! tmp3(K;I,J,alphaB) = sum_{gamma in gammaB} CoccT(K,gamma) tmp2(gamma,I,J,alphaB)
          m = nocc
          k = dimGamma
          n = dimAlpha*nocc**2
          call dgemm('N','N',m,n,k,1.0E0_realk,CoccT(1,GammaStart),m,tmp2,k,0.0E0_realk,tmp3,m)

          ! ooov(K,I,J;A) += sum_{alpha in alphaB} tmp3(K,I,J,alpha) Cvirt(alpha,A)
          m = nocc**3
          k = dimAlpha
          n = nvirt
          call dgemm('N','N',m,n,k,1.0E0_realk,tmp3,m,Cvirt(AlphaStart,1),nbasis,1.0E0_realk,ooov%elm1,m)

          ! Reorder: tmp1(A,alphaB;gammaB,I) --> tmp2(gammaB,I;A,alphaB)
          m = nvirt*dimAlpha
          n = dimGamma*nocc
          call mat_transpose(m,n,1.0E0_realk,tmp1,0.0E0_realk,tmp2)

          ! tmp1(B;I,A,alphaB) = sum_{gamma in gammaB} CvirtT(B,gamma) tmp2(gamma,I,A,alphaB)
          m = nvirt
          k = dimGamma
          n = nvirt*nocc*dimAlpha
          call dgemm('N','N',m,n,k,1.0E0_realk,CvirtT(1,GammaStart),m,tmp2,k,0.0E0_realk,tmp1,m)

#ifdef VAR_MPI

          ! mpi   : 1) tmp2(B,I,A,tile) = sum_{alpha in alphaB} tmp1(B,I,A,alpha) Cvirt(alpha,tile)
          !         2) vovv(B,I,A,C) += sum_{tile in CB} tmp2(B,I,A,tile)
          ! serial: vovv(B,I,A,C) += sum_{alpha in alphaB} tmp1(B,I,A,alpha) Cvirt(alpha,C)
          m = nocc*nvirt**2
          k = dimAlpha

          total_num_tiles = vovv%ntiles
          tile = 0

          do c=1,nvirt,tile_size

             tile = tile + 1

             call get_tileinfo_nels_fromarr8(nelms,vovv,i8*tile)
             tile_size_tmp = nelms/(nocc*nvirt**2)

             n = tile_size_tmp

             ! tmp2(B,I,A,tile) = sum_{alpha in alphaB} tmp1(B,I,A,alpha) Cvirt(alpha,tile)
             call dgemm('N','N',m,n,k,1.0E0_realk,tmp1,m,Cvirt(AlphaStart,c),nbasis,0.0E0_realk,tmp2,m)

             call time_start_phase(PHASE_COMM)
#ifdef VAR_HAVE_MPI3
             call tensor_lock_win(vovv,tile,'s',assert=mode)
#endif
             call get_residence_of_tile(dest,tile,vovv, idx_on_node = pos)

             if( alloc_in_dummy )then
                wi_idx = 1
                p      = pos - 1
             else
                wi_idx = tile
                p      = 0
             endif


             do first_el_c_block=1,ov2*tile_size_tmp,MAX_SIZE_ONE_SIDED
#ifndef VAR_HAVE_MPI3
                call tensor_lock_win(vovv,tile,'s',assert=mode)
#endif

                nel2t=MAX_SIZE_ONE_SIDED
                if(((ov2*tile_size_tmp-first_el_c_block)<MAX_SIZE_ONE_SIDED).and.&
                   &(mod(ov2*tile_size_tmp-first_el_c_block+1,i8*MAX_SIZE_ONE_SIDED)/=0))&
                   &nel2t=int(mod(ov2*tile_size_tmp,i8*MAX_SIZE_ONE_SIDED),kind=ls_mpik)

                call lsmpi_acc(tmp2(first_el_c_block:first_el_c_block+nel2t-1),nel2t,p+first_el_c_block,dest,vovv%wi(wi_idx))

#ifdef VAR_HAVE_MPI3
                call lsmpi_win_flush(vovv%wi(wi_idx),rank=dest,local=.true.)
#else
                call tensor_unlock_win(vovv,tile)
#endif
             enddo

#ifdef VAR_HAVE_MPI3
             call tensor_unlock_win(vovv,tile)
#endif
             call time_start_phase(PHASE_WORK)

          end do

#else

          m = nocc*nvirt**2
          k = dimAlpha
          n = nvirt
          call dgemm('N','N',m,n,k,1.0E0_realk,tmp1,m,Cvirt(AlphaStart,1),nbasis,1.0E0_realk,vovv%elm1,m)

#endif

       end do BatchAlpha
    end do BatchGamma

#ifdef VAR_MPI

    if (infpar%lg_nodtot .gt. 1) then

#ifdef VAR_PTR_RESHAPE
       dummy2(1:(i8*nocc*nocc)*nocc*nvirt) => ooov%elm1
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
       call c_f_pointer(c_loc(ooov%elm1(1)),dummy2,[(i8*nocc*nocc)*nocc*nvirt])
#else
       call lsquit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
       
       call time_start_phase(PHASE_IDLE)
       call lsmpi_barrier(infpar%lg_comm)

       ! now, reduce o^3v integrals onto master
       call time_start_phase(PHASE_COMM)
       call lsmpi_allreduce(dummy2,o3v,infpar%lg_comm) 
       call time_start_phase(PHASE_WORK)

    end if

    ! dealloc distribution array
    call mem_dealloc(distribution)

#endif

    ! free stuff
    ! **********
    call mem_dealloc(tmp1)
    call mem_dealloc(tmp2)
    call mem_dealloc(tmp3)

    if (DECinfo%useichor) then
       call FREE_SCREEN_ICHORERI()
       call mem_dealloc(AOGammabatchinfo)
       call mem_dealloc(AOAlphabatchinfo)
    else
       call free_decscreen(DECSCREEN)
       call mem_dealloc(orb2batchGamma)
       call mem_dealloc(batchdimGamma)
       call mem_dealloc(batchsizeGamma)
       call mem_dealloc(batchindexGamma)
       do idx=1,nbatchesGamma
          call mem_dealloc(batch2orbGamma(idx)%orbindex)
       end do
       call mem_dealloc(batch2orbGamma)
       call mem_dealloc(orb2batchAlpha)
       call mem_dealloc(batchdimAlpha)
       call mem_dealloc(batchsizeAlpha)
       call mem_dealloc(batchindexAlpha)
       do idx=1,nbatchesAlpha
          call mem_dealloc(batch2orbAlpha(idx)%orbindex)
          batch2orbAlpha(idx)%orbindex => null()
       end do
       call mem_dealloc(batch2orbAlpha)
       nullify(mylsitem%setting%LST_GAB_LHS)
       nullify(mylsitem%setting%LST_GAB_RHS)
    endif
    call mem_dealloc(CoccT)
    call mem_dealloc(CvirtT)

    ! finally, reorder ooov(K,I,J,A) --> ooov(I,J,K,A)
    call tensor_reorder(ooov,[2,3,1,4])

    if (master) call LSTIMER('CCSD(T) INT (ABC)',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine get_CCSDpT_integrals_abc


  !> \brief Get optimal batch sizes to be used in get_CCSDpT_integrals
  !> using the available memory.
  !> \author Kasper Kristensen & Janus Eriksen
  !> \date September 2011, rev. October 2012
  subroutine get_optimal_batch_sizes_ccsdpt_integrals(mylsitem,nbasis,nocc,nvirt,alphadim,gammadim,&
        & size1,size2,size3,adapt_to_nnodes,abc,tile_size)

     implicit none

     !> Integral info
     type(lsitem), intent(inout) :: mylsitem
     !> Number of AO basis functions
     integer,intent(in) :: nbasis
     !> Number of occupied (AOS) orbitals
     integer,intent(in) :: nocc
     !> Number of virt (AOS) orbitals
     integer,intent(in) :: nvirt
     !> Max size for AO alpha batch
     integer,intent(inout) :: alphadim
     !> Max size for AO gamma batch
     integer,intent(inout) :: gammadim
     !> Dimension of temporary array 1
     integer(kind=long),intent(inout) :: size1
     !> Dimension of temporary array 2
     integer(kind=long),intent(inout) :: size2
     !> Dimension of temporary array 3
     integer(kind=long),intent(inout) :: size3
     !> choose to split if more nodes are available than necessary
     logical,intent(in) :: adapt_to_nnodes
     !> is this for the abc partitioning?
     logical, intent(in) :: abc
     !> tile_size for abc partitioning
     integer, intent(in) :: tile_size
     !> memory reals
     real(realk) :: MemoryNeeded, MemoryAvailable
     integer :: MaxAObatch, MinAOBatch, AlphaOpt, GammaOpt,alpha,gamma,iAO
     integer(kind=ls_mpik) :: nnod,me
     logical :: master
     ! Memory currently available
     ! **************************
     call get_currently_available_memory(MemoryAvailable)
     ! Note: We multiply by 95 % to be on the safe side!
     MemoryAvailable = 0.95*MemoryAvailable

     nnod = 1
     me   = 0
#ifdef VAR_MPI
     nnod = infpar%lg_nodtot
     me   = infpar%lg_mynum
     call lsmpi_reduce_realk_min(MemoryAvailable,infpar%master,infpar%lg_comm)
#endif

     master = (me == 0)



     if ( master ) then

        ! Maximum and minimum possible batch sizes
        ! ****************************************

        ! The largest possible AO batch is the number of basis functions
        MaxAObatch = nbasis

        ! The smallest possible AO batch depends on the basis set
        ! (More precisely, if all batches are made as small as possible, then the
        !  call below determines the largest of these small batches).
        if (DECinfo%useichor) then
           !Determine the minimum allowed AObatch size MinAObatch
           !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
           !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch = 6 (the 2*(Px,Py,Pz))
           !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
           !'R'  !Specifies that it is the Regular AO basis that should be batched
           iAO = 1 !the center that the batching should occur on.  
           call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
        else
           call determine_maxBatchOrbitalsize(DECinfo%output,mylsitem%setting,MinAObatch,'R')
        endif

        ! Initialize batch sizes to be the minimum possible and then start increasing sizes below
        AlphaDim=MinAObatch
        GammaDim=MinAObatch

        GammaOpt = 0
        AlphaOpt = 0

        ! Gamma batch size
        ! =================================
        GammaLoop: do gamma = MaxAObatch,MinAOBatch,-1

           call get_max_arraysizes_for_ccsdpt_integrals(alphadim,gamma,nbasis,nocc,nvirt,&
              & size1,size2,size3,MemoryNeeded,abc,tile_size)

           if(MemoryNeeded < MemoryAvailable .or. (gamma==minAObatch) ) then
              if(adapt_to_nnodes)then
                 if( (nbasis/gamma)*(nbasis/MinAOBatch) > nnod * 3 )then

                    GammaOpt = gamma
                    exit GammaLoop

                 endif
              else

                 GammaOpt = gamma
                 exit GammaLoop

              endif
           end if

        end do GammaLoop

        if (GammaOpt .eq. 0) then

           GammaOpt = GammaDim

        endif 

        ! If gamma batch size was set manually we use that value instead
        if(DECinfo%ccsdGbatch/=0) then

           write(DECinfo%output,*) 'Gamma batch size was set manually, use that value instead!'
           GammaOpt=DECinfo%ccsdGbatch

        end if 

        ! The optimal gamma batch size is GammaOpt.
        ! We now find the maximum possible gamma batch size smaller than or equal to GammaOpt
        ! and store this number in gammadim.
        call determine_MaxOrbitals(DECinfo%output,mylsitem%setting,GammaOpt,gammadim,'R')


        ! Largest possible alpha batch size
        ! =================================
        AlphaLoop: do alpha = MaxAObatch,MinAOBatch,-1

           call get_max_arraysizes_for_ccsdpt_integrals(alpha,gammadim,nbasis,nocc,nvirt,&
              & size1,size2,size3,MemoryNeeded,abc,tile_size)

           if(MemoryNeeded < MemoryAvailable .or. (alpha==minAObatch) ) then

              if( adapt_to_nnodes  )then

                 if( (nbasis/GammaOpt)*(nbasis/alpha) > nnod * 3)then

                    AlphaOpt = alpha
                    exit AlphaLoop

                 endif
              else

                 AlphaOpt = alpha
                 exit AlphaLoop

              endif
           end if

        end do AlphaLoop

        if (AlphaOpt .eq. 0) then

           AlphaOpt = AlphaDim

        endif

        ! If alpha batch size was set manually we use that value instead
        if(DECinfo%ccsdAbatch/=0) then

           write(DECinfo%output,*) 'Alpha batch size was set manually, use that value instead!'
           AlphaOpt=DECinfo%ccsdAbatch

        end if

        ! The optimal alpha batch size is AlphaOpt.
        ! We now find the maximum possible alpha batch size smaller than or equal to AlphaOpt
        ! and store this number in alphadim.
        call determine_MaxOrbitals(DECinfo%output,mylsitem%setting,AlphaOpt,alphadim,'R')

        ! Print out and sanity check
        ! ==========================

        write(DECinfo%output,*) '======================================================================='
        write(DECinfo%output,*) '                     CCSD(T) INTEGRALS: MEMORY SUMMARY                 '
        if (abc) then
           write(DECinfo%output,*) '                             ABC partitioning                          '
        else
           write(DECinfo%output,*) '                             IJK partitioning                          '
        endif
        write(DECinfo%output,*) '======================================================================='
        write(DECinfo%output,*)
        write(DECinfo%output,*) 'To be on the safe side we use only 95% of the estimated available memory'
        write(DECinfo%output,*)
        write(DECinfo%output,'(1X,a,g10.3)') '95% of available memory (GB)            =', MemoryAvailable
        write(DECinfo%output,*)
        write(DECinfo%output,'(1X,a,i8)')    'Number of atomic basis functions        =', nbasis
        write(DECinfo%output,'(1X,a,i8)')    'Number of occupied orbitals             =', nocc
        write(DECinfo%output,'(1X,a,i8)')    'Number of virtual  orbitals             =', nvirt
        write(DECinfo%output,'(1X,a,i8)')    'Maximum alpha batch dimension           =', alphadim
        write(DECinfo%output,'(1X,a,i8)')    'Maximum gamma batch dimension           =', gammadim
        write(DECinfo%output,'(1X,a,g14.3)') 'Size of tmp array 1                     =', size1*realk*1.0E-9_realk
        write(DECinfo%output,'(1X,a,g14.3)') 'Size of tmp array 2                     =', size2*realk*1.0E-9_realk
        write(DECinfo%output,'(1X,a,g14.3)') 'Size of tmp array 3                     =', size3*realk*1.0E-9_realk
        write(DECinfo%output,*)
        write(DECinfo%output,*)

     endif

#ifdef VAR_MPI
     call ls_mpibcast(Gammadim,infpar%master,infpar%lg_comm)
     call ls_mpibcast(Alphadim,infpar%master,infpar%lg_comm)
#endif

     ! Sanity check
     call get_max_arraysizes_for_ccsdpt_integrals(alphadim,gammadim,nbasis,nocc,nvirt,&
        & size1,size2,size3,MemoryNeeded,abc,tile_size)

     if(MemoryNeeded > MemoryAvailable) then
        write(DECinfo%output,*) 'Requested/available memory: ', MemoryNeeded, MemoryAvailable
        call lsquit('CCSD(T) integrals: Insufficient memory!',-1)
     end if


  end subroutine get_optimal_batch_sizes_ccsdpt_integrals



  !> \brief Get sizes of temporary arrays used in CCSD(T) integral routine (get_CCSDpT_integrals)
  !> with the chosen AO batch sizes.
  !> NOTE: If get_CCSDpT_integrals is modified, this routine must be changed accordingly!
  !> \author Kasper Kristensen & Janus Eriksen
  !> \date September 2011, rev. October 2012
  subroutine get_max_arraysizes_for_ccsdpt_integrals(alphadim,gammadim,nbasis,nocc,nvirt,&
                     & size1,size2,size3,mem,abc,tile_size)
    implicit none
    !> Max size for AO alpha batch
    integer,intent(in) :: alphadim
    !> Max size for AO gamma batch
    integer,intent(in) :: gammadim
    !> Number of AO basis functions
    integer,intent(in) :: nbasis
    !> Number of occupied (AOS) orbitals
    integer,intent(in) :: nocc
    !> Number of virt (AOS) orbitals
    integer,intent(in) :: nvirt
    !> Dimension of temporary array 1
    integer(kind=long),intent(inout) :: size1
    !> Dimension of temporary array 2
    integer(kind=long),intent(inout) :: size2
    !> Dimension of temporary array 3
    integer(kind=long),intent(inout) :: size3
    !> Tot size of temporary arrays (in GB)
    real(realk), intent(inout) :: mem
    !> is this for the abc partitioning?
    logical, intent(in) :: abc
    !> tle_size only relevant for abc partitioning
    integer, intent(in) :: tile_size 
    real(realk) :: GB
    integer(kind=long) :: tmpI
    GB = 1073741824.0E0_realk ! 1 GB
    ! Array sizes needed in get_CCSDpT_integrals are checked and the largest one is found
 
    ! Tmp array 1
    if (abc) then

       size1 = i8*alphadim*gammadim*nbasis*nbasis
       tmpI = i8*alphadim*gammadim*nocc*nvirt
       size1 = max(size1,tmpI)
! this one is big
       tmpI = i8*alphadim*nocc*nvirt**2
       size1 = max(size1,tmpI)

    else

       size1 = i8*alphadim*gammadim*nbasis*nbasis
       tmpI = i8*nvirt**2*gammadim*alphadim
       size1 = max(size1,tmpI)
       tmpI = i8*nvirt*nocc*gammadim*alphadim
       size1 = max(size1,tmpI)
       tmpI = i8*nvirt*nocc**2*alphadim
       size1 = max(size1,tmpI)
       tmpI = i8*nvirt**3
       size1 = max(size1,tmpI)

    endif
  
    ! tmp array 2
    if (abc) then

       size2 = i8*alphadim*gammadim*nbasis*nocc
       tmpI = i8*alphadim*gammadim*nocc**2
       size2 = max(size2,tmpI)
       tmpI = i8*alphadim*gammadim*nocc*nvirt
       size2 = max(size2,tmpI)
#ifdef VAR_MPI
       tmpI = i8*nocc*nvirt**2*tile_size
       size2 = max(size2,tmpI)
#endif

    else

       size2 = i8*alphadim*gammadim*nbasis*nvirt
       tmpI = i8*alphadim*gammadim*nvirt*nocc
       size2 = max(size2,tmpI)
       tmpI = i8*nvirt**3
       size2 = max(size2,tmpI)

    endif
  
    ! Tmp array3
    if (abc) then

       size3 = i8*alphadim*gammadim*nocc**2
       tmpI = i8*alphadim*nocc**3
       size3 = max(size3,tmpI)
 
    else

       size3 = i8*alphadim*gammadim*nvirt**2
! this one is big
       tmpI = i8*alphadim*nvirt**3
       size3 = max(size3,tmpI)

    endif

    ! Size = size1+size2+size3,  convert to GB
    mem = realk*(size1+size2+size3)/GB


  end subroutine get_max_arraysizes_for_ccsdpt_integrals

  subroutine ccsdpt_info(nbasis,nocc,nvirt,print_frags,abc,ijk_nbuffs,abc_nbuffs,abc_tile_size,nodtotal)

      use iso_c_binding
      implicit none

      integer, intent(in) :: nbasis,nocc,nvirt
      logical, intent(in) :: print_frags,abc
      integer, intent(in) :: nodtotal
      logical :: ijk,manual_ijk,manual_abc_1,manual_abc_2,gpu
      integer, intent(inout) :: ijk_nbuffs,abc_nbuffs,abc_tile_size
      integer :: num_gpu,me
      logical :: master
#ifdef VAR_OPENACC
      integer(kind=acc_device_kind) :: acc_device_type
#endif
      logical :: acc_async
      real(realk) :: free_cpu ! in gb
      integer(c_size_t) :: total_gpu,free_gpu ! in bytes
      real(realk), parameter :: gb =  1073741824.0E0_realk ! 1 GB
      integer, parameter :: ijk_default = 1000000, abc_default = 1000000

      ijk_nbuffs = 0
      abc_nbuffs = 0
      abc_tile_size = 0
      manual_abc_1 = .false.
      manual_abc_2 = .false.
      manual_ijk = .false.
      gpu = .false.
      num_gpu = 0
      free_cpu = 0.0E0_realk
      total_gpu = 0
      free_gpu = 0
      acc_async = .true.

      if (DECinfo%acc_sync) acc_async = .false.

      call get_currently_available_memory(free_cpu)

      me   = 0
#ifdef VAR_MPI
      me   = infpar%lg_mynum
      call lsmpi_reduce_realk_min(free_cpu,infpar%master,infpar%lg_comm)
#endif

      master = (me .eq. 0)

      if (master) then

#ifdef VAR_OPENACC

         acc_device_type = acc_get_device_type() 
         num_gpu = acc_get_num_devices(acc_device_type)

#ifdef VAR_CUDA
         call get_dev_mem(total_gpu,free_gpu)
#endif

#endif

         if (num_gpu .gt. 0) gpu = .true.
   
         if (abc) then
            ijk = .false.
         else
            ijk = .true.
         endif
   
         if (ijk) then
   
            if (DECinfo%ijk_nbuffs .lt. ijk_default) then
    
                manual_ijk = .true.
                ijk_nbuffs = DECinfo%ijk_nbuffs
   
            endif
   
            if (manual_ijk) then
   
               if (ijk_nbuffs .lt. 3) call lsquit('manually set ijk_nbuffs (NBUFFS_IJK) .lt. 3 - aborting...',DECinfo%output)
   
            else
   
!               ! here; determine ijk_nbuffs based on available cpu/gpu memory
!               if (ijk_nbuffs .eq. ijk_default) call new_ijk_nbuffs_routine()
               ijk_nbuffs = 6
   
            endif
   
         else ! abc == .true.
   
            if (DECinfo%abc_nbuffs .lt. abc_default) then
   
                manual_abc_1 = .true.
                abc_nbuffs = DECinfo%abc_nbuffs
   
            endif
   
            if (manual_abc_1) then
   
               if (abc_nbuffs .lt. 3) call lsquit('manually set abc_nbuffs (NBUFFS_ABC) .lt. 3 - aborting...',DECinfo%output)
   
            else
   
!               ! here; determine ijk_nbuffs based on available cpu/gpu memory
!               if (abc_nbuffs .eq. abc_default) call new_abc_nbuffs_routine()
               abc_nbuffs = 6
   
            endif
   
            if (DECinfo%abc_tile_size .lt. abc_default) then
   
               manual_abc_2 = .true.
               abc_tile_size = DECinfo%abc_tile_size
   
            endif
   
            if (manual_abc_2) then
   
               if (abc_tile_size .lt. 1) call lsquit('manually set tile size (.ABC_TILE) .lt. 1 - aborting...',DECinfo%output)
   
            else
   
               call abc_tile_size_routine(nocc,nvirt,print_frags,free_cpu,nodtotal,abc_nbuffs,abc_tile_size)
   
            endif
   
            if (abc_tile_size .gt. nvirt) call lsquit('manually set tile size (.ABC_TILE) .gt. nvirt - aborting...',DECinfo%output)
   
         endif
   
         write(DECinfo%output,'(/,a)') '-----------------------------'
         write(DECinfo%output,'(a)')   '      CCSD(T) information    '
         write(DECinfo%output,'(a,/)') '-----------------------------'
#ifdef VAR_MPI
         write(DECinfo%output,'(a,i4)')     'Number of nodes in lg  = ',nodtotal
#endif
         write(DECinfo%output,'(a,l4)')     'Print frag. energies   = ',print_frags
         write(DECinfo%output,'(a,l4)')     'IJK partitioning       = ',ijk
         if (ijk) then
            if (manual_ijk) then
               write(DECinfo%output,'(a,i4)')     'Input # IJK buffers    = ',ijk_nbuffs
            else
               write(DECinfo%output,'(a,i4)')     '# IJK buffers          = ',ijk_nbuffs
            endif
         endif
         write(DECinfo%output,'(a,l4)')     'ABC partitioning       = ',abc
         if (abc) then
            if (manual_abc_1) then
               write(DECinfo%output,'(a,i4)')     'Input # ABC buffers    = ',abc_nbuffs
            else
               write(DECinfo%output,'(a,i4)')     '# ABC buffers          = ',abc_nbuffs
            endif
            if (manual_abc_2) then
               write(DECinfo%output,'(a,i4)')     'Input ABC tile size    = ',abc_tile_size
            else
               write(DECinfo%output,'(a,i4)')     'ABC tile size          = ',abc_tile_size
            endif
         endif
         write(DECinfo%output,'(a,g11.4)')     'Free CPU memory (GB)   = ',free_cpu
         write(DECinfo%output,'(a,l4)')     'Are we using GPUs?     = ',gpu
         if (gpu) then
            write(DECinfo%output,'(a,l4)')     'Asynchronous OpenACC?  = ',acc_async
            write(DECinfo%output,'(a,i4)')     'Number of GPUs         = ',num_gpu
            write(DECinfo%output,'(a,g11.4)')     'Total GPU memory (GB)  = ',total_gpu / gb
            write(DECinfo%output,'(a,g11.4)')     'Free GPU memory (GB)   = ',free_gpu / gb
            print *,'Total GPU memory (Bytes) = ',total_gpu
            print *,'Free GPU memory  (Bytes) = ',free_gpu
         endif
         write(DECinfo%output,*)
         write(DECinfo%output,*)

      endif

  end subroutine ccsdpt_info

  subroutine abc_tile_size_routine(nocc,nvirt,print_frags,free_cpu,nodtotal,abc_nbuffs,abc_tile_size)

      implicit none

      integer, intent(in) :: nocc,nvirt,nodtotal,abc_nbuffs
      real(realk), intent(in) :: free_cpu
      logical, intent(in) :: print_frags
      integer, intent(inout) :: abc_tile_size
      real(realk) :: mem_avail_start,mem_accum_tmp,mem_est_avail_tmp,mem_est_avail,mem_vovv_pdm,mem_vovv_local
      integer(kind=long) :: ccsd_doubles,ccsd_doubles_portions,ooov,oovv,vovv_total
      integer(kind=long) :: eivalocc,eivalvirt,trip_ampls,ccsdpt_singles,ccsdpt_doubles
      integer(kind=long) :: vovv_pdm,vovv_local
      integer(kind=long) :: mem_int_tmp
      integer(kind=long) :: max_abc_tile_size,ts
      integer :: remainder_1,remainder_2,num_tiles_tot,num_tiles_node
      real(realk), parameter :: GB = 1073741824.0E0_realk ! 1 GB = 1024.**3

      ! available memory - note that we multiply by 95 % to be on the safe side!
      mem_avail_start = 0.95 * free_cpu

      ! ccsd quantities
      ccsd_doubles = i8*nocc**2*nvirt**2
      ccsd_doubles_portions = i8*3*nocc*nvirt**2

      ! local integrals
      ooov = i8*nocc**3*nvirt
      oovv = i8*nocc**2*nvirt**2

      ! total distributed integrals
      vovv_total = i8*nocc*nvirt**3

      ! orbital energies
      eivalocc = i8*nocc
      eivalvirt = i8*nvirt

      ! triples amplitudes and temp array
      trip_ampls = i8*2*nocc**3

      ! ccsdpt intermediates
      ccsdpt_singles = i8*nocc*nvirt
      ccsdpt_doubles = i8*2*nocc**2*nvirt**2

      ! temp sum of integer elements
      mem_int_tmp = ccsd_doubles + &
                  & ccsd_doubles_portions + &
                  & ooov + &
                  & oovv + &
                  & eivalocc + &
                  & eivalvirt + &
                  & trip_ampls

      if (print_frags) mem_int_tmp = mem_int_tmp + ccsdpt_singles + ccsdpt_doubles 

      ! estimate the memory needed for allocations
      mem_accum_tmp = realk*mem_int_tmp / GB

      ! estimate available memory AFTER allocations, but BEFORE vovv allocation
      mem_est_avail_tmp = mem_avail_start - mem_accum_tmp
      if (mem_est_avail_tmp .lt. 0.0E0_realk) call lsquit('mem_est_avail_tmp .lt. 0 GB - aborting...',DECinfo%output) 

      ! max tile_size
      max_abc_tile_size = int(nvirt / nodtotal)

      do ts = max_abc_tile_size,1,-1

         ! how many tiles are there in total?
         num_tiles_tot = int(nvirt / ts)

         ! modulo_1
         remainder_1 = mod(nvirt,ts)

         ! update total number of tiles
         if (remainder_1 .gt. 0) num_tiles_tot = num_tiles_tot + 1

         ! how many tiles per node in PDM?
         num_tiles_node = int(num_tiles_tot / nodtotal)

         ! modulo_2
         remainder_2 = mod(num_tiles_tot,nodtotal)

         ! update number of tiles per node
         if (remainder_2 .gt. 0) num_tiles_node = num_tiles_node + 1
   
         ! calculate the PDM memory requirements for the given tile_size
         vovv_pdm = i8*num_tiles_node*(nocc*nvirt**2*ts)
         mem_vovv_pdm = realk*vovv_pdm / GB

         ! calculate the local memory requirements for the given tile_size and the given number of tiles
         vovv_local = i8*abc_nbuffs*(nocc*nvirt**2*ts)
         mem_vovv_local = realk*vovv_local / GB

         ! estimate available memory AFTER vovv allocation
         mem_est_avail = mem_est_avail_tmp - (mem_vovv_pdm + mem_vovv_local)

         if ((mem_est_avail .lt. 0.0E0_realk)) then

            if (ts .eq. 1) then

               print *,'ts                = ',ts
               print *,'num_tiles_tot     = ',num_tiles_tot
               print *,'num_tiles_node    = ',num_tiles_node
               print *,'mem_avail_start   = ',mem_avail_start
               print *,'mem_est_avail_tmp = ',mem_est_avail_tmp
               print *,'mem_vovv_pdm      = ',mem_vovv_pdm
               print *,'mem_vovv_local    = ',mem_vovv_local
               print *,'mem_est_avail     = ',mem_est_avail

               call lsquit('mem_est_avail .lt. 0 GB for smallest possible abc_tile_size - aborting...',DECinfo%output)

            endif

            cycle

         else

            abc_tile_size = ts

            return

         endif

      enddo

  end subroutine abc_tile_size_routine

!endif mod_unreleased
#endif

  subroutine dummy_ccsdpt_routine()

  end subroutine dummy_ccsdpt_routine

end module ccsdpt_module

#ifdef MOD_UNRELEASED

  !> \brief slaves enter here from lsmpi_slave (or dec_lsmpi_slave) and need to get to work 
  !> \author Janus Juul Eriksen
  !> \date x-mas 2012
#ifdef VAR_MPI

  subroutine ccsdpt_slave()

  use infpar_module
  use lsmpi_type
  use decmpi_module

  use precision
  use dec_typedef_module
  use memory_handling
  use lstiming!, only: lstimer
  use typedeftype, only: Lsitem,lssetting

  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
  use array2_simple_operations, only: array2_init_plain,array2_free 
  use array4_simple_operations, only: array4_init_standard,array4_free
  use atomic_fragment_operations
  use ccsdpt_module, only: ccsdpt_driver

    implicit none
    integer :: nocc, nvirt,nbasis
    real(realk), pointer :: ppfock(:,:), qqfock(:,:), Co(:,:), Cv(:,:)
    type(tensor) :: ccsdpt_t1
    type(tensor) :: vovo,ccsd_t2, ccsdpt_t2
    real(realk) :: ccsdpt_e4
    type(lsitem) :: mylsitem
    logical :: print_frags,abc

    abc = .false.
    print_frags = .false.

    call time_start_phase(PHASE_COMM)

    ! call ccsd(t) data routine in order to receive data from master
    call mpi_communicate_ccsdpt_calcdata(nocc,nvirt,nbasis,vovo%elm4,ccsd_t2%elm4,mylsitem,print_frags,abc)

    !FIXME: split MPI messages!!!!!!!!!!
    ! init and receive vovo and ccsd_doubles array structures
    if (abc) then

       call tensor_init(vovo,[nocc,nocc,nvirt,nvirt],4)
       call tensor_init(ccsd_t2, [nocc,nocc,nvirt,nvirt],4)

       call ls_mpibcast(vovo%elm4,nocc,nocc,nvirt,nvirt,infpar%master,infpar%lg_comm)
       call ls_mpibcast(ccsd_t2%elm4,nocc,nocc,nvirt,nvirt,infpar%master,infpar%lg_comm)

    else

       call tensor_init(vovo, [nvirt,nvirt,nocc,nocc],4)
       call tensor_init(ccsd_t2, [nvirt,nvirt,nocc,nocc],4)

       call ls_mpibcast(vovo%elm4,nvirt,nvirt,nocc,nocc,infpar%master,infpar%lg_comm)
       call ls_mpibcast(ccsd_t2%elm4,nvirt,nvirt,nocc,nocc,infpar%master,infpar%lg_comm)

    endif

    if (print_frags) then
 
       ! init ccsd(t) singles and ccsd(t) doubles
       if (abc) then

          call tensor_init(ccsdpt_t1, [nocc,nvirt],2)
          call tensor_init(ccsdpt_t2, [nocc,nocc,nvirt,nvirt],4)
          
       else

          call tensor_init(ccsdpt_t1, [nvirt,nocc],2)
          call tensor_init(ccsdpt_t2, [nvirt,nvirt,nocc,nocc],4)

       endif

    else

       ! init ccsd(t) singles
       if (abc) then

          call tensor_init(ccsdpt_t1, [nocc,nvirt],2)

       else

          call tensor_init(ccsdpt_t1, [nvirt,nocc],2)

       endif
       ccsdpt_e4 = 0.0E0_realk

    endif

    call time_start_phase(PHASE_WORK)

    ! now enter the ccsd(t) driver routine
    if (print_frags) then

       call ccsdpt_driver(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,mylsitem,vovo,ccsd_t2,&
                               & ccsdpt_t1,print_frags,abc,ccsdpt_doubles=ccsdpt_t2)

    else

       call ccsdpt_driver(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,mylsitem,vovo,ccsd_t2,&
                               & ccsdpt_t1,print_frags,abc,e4=ccsdpt_e4)

    endif

    call time_start_phase(PHASE_WORK)

    ! now, release all amplitude arrays, both ccsd and ccsd(t)
    call tensor_free(vovo)
    call tensor_free(ccsd_t2)

    if (print_frags) then

       call tensor_free(ccsdpt_t1)
       call tensor_free(ccsdpt_t2)

    else

       call tensor_free(ccsdpt_t1)

    endif

    call ls_free(mylsitem)

  end subroutine ccsdpt_slave

#endif
!endif mod_unreleased
#endif
