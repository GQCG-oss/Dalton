!> @file
!> DEC-CCSD(T) routines
!> \brief: ccsd(t) module
!> \author: Janus Juul Eriksen
!> \date: 2012-2014, Aarhus
module ccsdpt_module
  use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer, c_size_t

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
  use tensor_basic_module
  use background_buffer_module
#ifdef VAR_OPENACC
  use openacc
#endif
#if defined(VAR_CUDA) || defined(VAR_OPENACC)
  use gpu_interfaces
#endif

  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
#ifdef VAR_REAL_SP
  use sp_ccsdpt_kernels_module
#endif
  use ccsdpt_kernels_module
  use decmpi_module
  use dec_workarounds_module
  use crop_tools_module
  use cc_tools_module
  use dec_fragment_utils
  use array2_simple_operations
  use array3_simple_operations
  use array4_simple_operations
  
#ifdef MOD_UNRELEASED
  public :: ccsdpt_driver,ccsdpt_info,ccsdpt_energy_e5_frag,&
       & ccsdpt_energy_e5_pair, ccsdpt_energy_e5_ddot, &
       & ccsdpt_decnp_e4_frag, ccsdpt_decnp_e5_frag
  private
#endif

contains

#ifdef MOD_UNRELEASED

  !> \brief: driver routine for dec-ccsd(t)
  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  subroutine ccsdpt_driver(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,mylsitem,vovo_in,ccsd_doubles_in,&
                         & print_frags,abc,ccsdpt_singles,ccsdpt_doubles,e4,e5,ccsd_singles)

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
    type(tensor), intent(inout) :: ccsd_doubles_in
    type(tensor) :: ccsd_doubles
    !> incoming vovo integrals
    type(tensor), intent(inout) :: vovo_in
    type(tensor) :: vovo
    !> input for the actual triples computation
    logical :: print_frags,abc
    type(tensor),intent(inout),optional :: ccsdpt_singles
    type(tensor),intent(inout),optional :: ccsdpt_doubles
    real(realk),optional :: e4,e5
    type(tensor),intent(inout),optional :: ccsd_singles
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
    integer :: ijk_nbuffs,abc_nbuffs,ijk_tile_size,abc_tile_size
    !> orbital energies
    real(realk), pointer :: eivalocc(:), eivalvirt(:)
    !> MOs and unitary transformation matrices
    type(tensor) :: C_can_occ, C_can_virt, Uocc, Uvirt
    !> dimensions
    integer, dimension(2) :: occdims, virtdims, virtoccdims,occAO,virtAO
    integer, dimension(3) :: dims_aaa
    integer, dimension(4) :: dims_iaai, dims_aaii
    logical :: master, use_bg
#ifdef VAR_OPENACC
    !> device type
    integer(acc_device_kind) :: acc_device_type
#endif
    real(realk) :: tcpu,twall

    call time_start_phase(PHASE_WORK)
    
    use_bg = mem_is_background_buf_init()

    ! some quit statements
    if (nocc .gt. nvirt) call lsquit('CCSD(T) with nocc .gt. nvirt has not been implemented...',DECinfo%output)
    if (print_frags .and. DECinfo%pt_hack) call lsquit('print_frags .and. .PT_HACK is not allowed...',DECinfo%output) 
    if (print_frags .and. DECinfo%pt_single_prec) &
                     & call lsquit('print_frags .and. .PT_SINGLE_PREC is not allowed...',DECinfo%output)
#ifndef VAR_REAL_SP
    if (DECinfo%pt_single_prec) call lsquit('.PT_SINGLE_PREC only works for single prec. builds! re-compile...',DECinfo%output)
#endif

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

       if (present(e4) .and. present(e5)) then

          if (present(ccsd_singles)) then

             e4 = 0.0E0_realk
             e5 = 0.0E0_realk

          else

             call lsquit('print_frags == .false., but ccsd t1 is missing... aborting.',DECinfo%output)

          endif

       else

          call lsquit('print_frags == .false., but either e4 or e5 are missing... aborting.',DECinfo%output) 

       endif

    endif

    call mem_alloc(eivalocc,nocc)
    call mem_alloc(eivalvirt,nvirt)

    if (master) then
      ! *************************************
      ! get arrays for transforming integrals
      ! *************************************
      ! C_can_occ, C_can_virt:  MO coefficients for canonical basis
      ! Uocc, Uvirt: unitary transformation matrices for canonical --> local basis (and vice versa)
      ! note: Uocc and Uvirt have indices (local,canonical)

      call tensor_init(Uocc  ,occdims,  2, bg=use_bg )
      call tensor_init(Uvirt ,virtdims, 2, bg=use_bg )

      call tensor_init(C_can_occ  ,occAO , 2, bg=use_bg )
      call tensor_init(C_can_virt ,virtAO, 2, bg=use_bg )

      if (.not. DECinfo%pt_hack) then

         call get_canonical_integral_transformation_matrices(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,&
                            & C_can_occ%elm2,C_can_virt%elm2,Uocc%elm2,Uvirt%elm2,eivalocc,eivalvirt)

         if (.not. print_frags) then         

            if (abc) then
   
               call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,ov=ccsd_singles%elm1)
   
            else
   
               call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,vo=ccsd_singles%elm1)
   
            endif

         endif

      else

         call random_seed()
         call random_number(C_can_occ%elm2)
         call dscal(nbasis*nocc,1.0E-2_realk,C_can_occ%elm2,1)
         call random_seed()
         call random_number(C_can_virt%elm2)
         call dscal(nbasis*nvirt,1.0E-2_realk,C_can_virt%elm2,1)
         call random_seed()
         call random_number(eivalocc)
         call dscal(nocc,-1.0E-2_realk,eivalocc,1)
         call random_seed()
         call random_number(eivalvirt)
         call dscal(nvirt,1.0E-2_realk,eivalvirt,1)

      endif

    else !Slave only allocate

      call tensor_init(C_can_occ  ,occAO , 2, bg=use_bg )
      call tensor_init(C_can_virt ,virtAO, 2, bg=use_bg )

    endif

!#ifdef VAR_MPI
!      call tensor_minit(Uo,[nocc,nocc],2,local=.false.,atype="TDPD")
!      call tensor_minit(Uv,[nvirt,nvirt],2,local=.false.,atype="TDPD")
!      call tensor_convert(Uocc%val,Uo)
!      call tensor_convert(Uvirt%val,Uv)
!#endif

#ifdef VAR_MPI
    call time_start_phase(PHASE_COMM)

    if (nodtotal .gt. 1) then

       ! bcast the JOB specifier and distribute data to all the slaves within local group
       waking_the_slaves_info: if (master) then
   
          ! slaves are in lsmpi_slave routine (or corresponding dec_mpi_slave) and are now awaken
          call ls_mpibcast(CCSDPTSLAVE_INFO,infpar%master,infpar%lg_comm)
   
          call ccsdpt_info(nbasis,nocc,nvirt,print_frags,abc,ijk_nbuffs,abc_nbuffs,ijk_tile_size,abc_tile_size,nodtotal)
   
       end if waking_the_slaves_info

    else

       call ccsdpt_info(nbasis,nocc,nvirt,print_frags,abc,ijk_nbuffs,abc_nbuffs,ijk_tile_size,abc_tile_size,nodtotal)

    endif

    call time_start_phase(PHASE_WORK)
#else
    call ccsdpt_info(nbasis,nocc,nvirt,print_frags,abc,ijk_nbuffs,abc_nbuffs,ijk_tile_size,abc_tile_size,nodtotal)
#endif

    ! convert the ccsd doubles and the vovo from the distribution used in ccsd to the one we want to use for (t)
    if (master) call convert_ccsd_and_vovo(nocc,nvirt,nbasis,Uocc,Uvirt,ccsd_doubles_in,vovo_in,ccsd_doubles,vovo,&
                                         & nodtotal,abc,ijk_tile_size,abc_tile_size,use_bg)

#ifdef VAR_MPI
    call time_start_phase(PHASE_COMM)

    if (nodtotal .gt. 1) then

       ! bcast the JOB specifier and distribute data to all the slaves within local group
       waking_the_slaves_work: if (master) then
   
          ! slaves are in lsmpi_slave routine (or corresponding dec_mpi_slave) and are now awaken
          call ls_mpibcast(CCSDPTSLAVE_WORK,infpar%master,infpar%lg_comm)
   
          ! distribute ccsd doubles and fragment or full molecule quantities to the slaves
          if (print_frags) then

             call mpi_communicate_ccsdpt_calcdata(nocc,nvirt,nbasis,vovo,ccsd_doubles,&
                                                & mylsitem,print_frags,abc)

          else

             call mpi_communicate_ccsdpt_calcdata(nocc,nvirt,nbasis,vovo,ccsd_doubles,&
                                                & mylsitem,print_frags,abc,ccsd_singles)

          endif
   
       end if waking_the_slaves_work
   
       call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
       call ls_mpi_buffer(eivalocc,nocc,infpar%master)
       call ls_mpi_buffer(eivalvirt,nvirt,infpar%master)
       call ls_mpi_buffer(C_can_occ%elm2,nbasis,nocc,infpar%master)
       call ls_mpi_buffer(C_can_virt%elm2,nbasis,nvirt,infpar%master)
       call ls_mpi_buffer(ijk_nbuffs,infpar%master)
       call ls_mpi_buffer(abc_nbuffs,infpar%master)
       call ls_mpi_buffer(ijk_tile_size,infpar%master)
       call ls_mpi_buffer(abc_tile_size,infpar%master)
       call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

       if( .not. master )then
          ccsd_doubles = ccsd_doubles_in
          vovo         = vovo_in
       endif

    endif

    call time_start_phase(PHASE_WORK)
#endif

    ! ********************************************
    ! get vo³ and v³o integrals in proper sequence
    ! ********************************************
    ! note: the integrals are calculated in canonical basis

    if (abc) then

       call get_CCSDpT_integrals_abc(mylsitem,nbasis,nocc,nvirt,C_can_occ%elm2,C_can_virt%elm2,ooov,vovv,abc_tile_size)

    else

       call get_CCSDpT_integrals_ijk(mylsitem,nbasis,nocc,nvirt,C_can_occ%elm2,C_can_virt%elm2,ovoo,vvvo,ijk_tile_size)

    endif

    write(DECinfo%output,*) ''
    write(DECinfo%output,*) ''
    write(DECinfo%output,*) '=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*='
    write(DECinfo%output,*) '        Done with CCSD(T) integrals        '
    write(DECinfo%output,*) '*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*'
    write(DECinfo%output,*) ''

    ! release occ and virt canonical MOs
    call tensor_free(C_can_virt)
    call tensor_free(C_can_occ)

    ! ********************************
    ! begin actual triples calculation
    ! ********************************

    ! in all comments in the below, we employ the notation of eqs. (14.6.60) [with (i,j,k)/(a,c,d)]
    ! and (14.6.64).

    ! objective is three-fold:
    ! 1) calculate triples amplitudes, collect in array3 structures, trip_*** [canonical basis]
    ! 2) calculate ^{*}T^{a}_{i} and ^{*}T^{ab}_{ij} amplitudes in array2 and array4 structures, 
    !    ccsdpt_singles and ccsdpt_doubles [canonical basis]
    ! 3) transform ccsd_doubles, ccsdpt_singles and ccsdpt_doubles into local basis [local basis]

    ! *****************************************************
    ! ***************** trip generation *******************
    ! *****************************************************

    !************************************************************!
    ! here: the main (t) loop: this is where the magic happens! !
    !************************************************************!

#ifdef VAR_MPI

    call time_start_phase(PHASE_WORK)

    if (abc) then

       ! the parallel version of the abc-loop
       if (print_frags) then

          if (master) then

             if (nodtotal .gt. 1) then

                call abc_loop_par(nocc,nvirt,ooov%elm1,vovo,vovv,ccsd_doubles,&
                                & eivalocc,eivalvirt,nodtotal,abc_nbuffs,abc_tile_size,&
                                & ccsdpt_singles%elm1,ccsdpt_doubles%elm1)

             else

                call abc_loop_ser(nocc,nvirt,ooov%elm1,vovo%elm1,vovv%elm1,ccsd_doubles%elm1,&
                                & eivalocc,eivalvirt,ccsdpt_singles%elm1,ccsdpt_doubles%elm1)

             endif

          else

             call abc_loop_par(nocc,nvirt,ooov%elm1,vovo_in,vovv,ccsd_doubles_in,&
                             & eivalocc,eivalvirt,nodtotal,abc_nbuffs,abc_tile_size,&
                             & ccsdpt_singles%elm1,ccsdpt_doubles%elm1)

          endif

       else

          if (master) then

             if (nodtotal .gt. 1) then

                if (DECinfo%pt_single_prec) then

#ifdef VAR_REAL_SP
                   call sp_abc_loop_par(nocc,nvirt,ooov%elm1,vovo,vovv,ccsd_doubles,&
                                   & eivalocc,eivalvirt,nodtotal,abc_nbuffs,abc_tile_size,&
                                   & e4=e4,e5=e5,t1=ccsd_singles%elm1)
#endif

                else

                   call abc_loop_par(nocc,nvirt,ooov%elm1,vovo,vovv,ccsd_doubles,&
                                   & eivalocc,eivalvirt,nodtotal,abc_nbuffs,abc_tile_size,&
                                   & e4=e4,e5=e5,t1=ccsd_singles%elm1)

                endif

             else

                if (DECinfo%pt_single_prec) then

#ifdef VAR_REAL_SP
                   call sp_abc_loop_ser(nocc,nvirt,ooov%elm1,vovo%elm1,vovv%elm1,ccsd_doubles%elm1,&
                                   & eivalocc,eivalvirt,e4=e4,e5=e5,t1=ccsd_singles%elm1)
#endif

                else

                   call abc_loop_ser(nocc,nvirt,ooov%elm1,vovo%elm1,vovv%elm1,ccsd_doubles%elm1,&
                                   & eivalocc,eivalvirt,e4=e4,e5=e5,t1=ccsd_singles%elm1)

                endif

             endif

          else

             if (DECinfo%pt_single_prec) then

#ifdef VAR_REAL_SP
                call sp_abc_loop_par(nocc,nvirt,ooov%elm1,vovo_in,vovv,ccsd_doubles_in,&
                                & eivalocc,eivalvirt,nodtotal,abc_nbuffs,abc_tile_size,&
                                & e4=e4,e5=e5,t1=ccsd_singles%elm1)
#endif

             else

                call abc_loop_par(nocc,nvirt,ooov%elm1,vovo_in,vovv,ccsd_doubles_in,&
                                & eivalocc,eivalvirt,nodtotal,abc_nbuffs,abc_tile_size,&
                                & e4=e4,e5=e5,t1=ccsd_singles%elm1)

             endif

          endif

       endif

    else

       ! the parallel version of the ijk-loop
       if (print_frags) then
  
          if (master) then 

             if (nodtotal .gt. 1) then

                call ijk_loop_par(nocc,nvirt,ovoo%elm1,vovo,vvvo,ccsd_doubles,&
                                & eivalocc,eivalvirt,nodtotal,ijk_nbuffs,ijk_tile_size,&
                                & ccsdpt_singles%elm1,ccsdpt_doubles%elm1)

             else

                call ijk_loop_ser(nocc,nvirt,ovoo%elm1,vovo%elm1,vvvo%elm1,ccsd_doubles%elm1,&
                                & eivalocc,eivalvirt,ccsdpt_singles%elm1,ccsdpt_doubles%elm1)

             endif

          else

             call ijk_loop_par(nocc,nvirt,ovoo%elm1,vovo_in,vvvo,ccsd_doubles_in,&
                             & eivalocc,eivalvirt,nodtotal,ijk_nbuffs,ijk_tile_size,&
                             & ccsdpt_singles%elm1,ccsdpt_doubles%elm1)

          endif   

       else

          if (master) then

             if (nodtotal .gt. 1) then

                if (DECinfo%pt_single_prec) then

#ifdef VAR_REAL_SP
                   call sp_ijk_loop_par(nocc,nvirt,ovoo%elm1,vovo,vvvo,ccsd_doubles,&
                                   & eivalocc,eivalvirt,nodtotal,ijk_nbuffs,ijk_tile_size,&
                                   & e4=e4,e5=e5,t1=ccsd_singles%elm1)
#endif

                else

                   call ijk_loop_par(nocc,nvirt,ovoo%elm1,vovo,vvvo,ccsd_doubles,&
                                   & eivalocc,eivalvirt,nodtotal,ijk_nbuffs,ijk_tile_size,&
                                   & e4=e4,e5=e5,t1=ccsd_singles%elm1)

                endif

             else

                if (DECinfo%pt_single_prec) then

#ifdef VAR_REAL_SP
                   call sp_ijk_loop_ser(nocc,nvirt,ovoo%elm1,vovo%elm1,vvvo%elm1,ccsd_doubles%elm1,&
                                   & eivalocc,eivalvirt,e4=e4,e5=e5,t1=ccsd_singles%elm1)
#endif

                else

                   call ijk_loop_ser(nocc,nvirt,ovoo%elm1,vovo%elm1,vvvo%elm1,ccsd_doubles%elm1,&
                                   & eivalocc,eivalvirt,e4=e4,e5=e5,t1=ccsd_singles%elm1)

                endif

             endif

          else

             if (DECinfo%pt_single_prec) then

#ifdef VAR_REAL_SP
                call sp_ijk_loop_par(nocc,nvirt,ovoo%elm1,vovo_in,vvvo,ccsd_doubles_in,&
                                & eivalocc,eivalvirt,nodtotal,ijk_nbuffs,ijk_tile_size,&
                                & e4=e4,e5=e5,t1=ccsd_singles%elm1)
#endif

             else

                call ijk_loop_par(nocc,nvirt,ovoo%elm1,vovo_in,vvvo,ccsd_doubles_in,&
                                & eivalocc,eivalvirt,nodtotal,ijk_nbuffs,ijk_tile_size,&
                                & e4=e4,e5=e5,t1=ccsd_singles%elm1)

             endif

          endif   

       endif

    endif

    call time_start_phase(PHASE_WORK)

#else

    if (abc) then

       ! the serial version of the abc-loop
       if (print_frags) then

          call abc_loop_ser(nocc,nvirt,ooov%elm1,vovo%elm1,vovv%elm1,ccsd_doubles%elm1,&
                          & eivalocc,eivalvirt,ccsdpt_singles%elm1,ccsdpt_doubles%elm1)

       else

          if (DECinfo%pt_single_prec) then

#ifdef VAR_REAL_SP
             call sp_abc_loop_ser(nocc,nvirt,ooov%elm1,vovo%elm1,vovv%elm1,ccsd_doubles%elm1,&
                             & eivalocc,eivalvirt,e4=e4,e5=e5,t1=ccsd_singles%elm1)
#endif

          else

             call abc_loop_ser(nocc,nvirt,ooov%elm1,vovo%elm1,vovv%elm1,ccsd_doubles%elm1,&
                             & eivalocc,eivalvirt,e4=e4,e5=e5,t1=ccsd_singles%elm1)

          endif

       endif

    else

       ! the serial version of the ijk-loop
       if (print_frags) then
   
          call ijk_loop_ser(nocc,nvirt,ovoo%elm1,vovo%elm1,vvvo%elm1,ccsd_doubles%elm1,&
                          & eivalocc,eivalvirt,ccsdpt_singles%elm1,ccsdpt_doubles%elm1)
   
       else

          if (DECinfo%pt_single_prec) then

#ifdef VAR_REAL_SP
             call sp_ijk_loop_ser(nocc,nvirt,ovoo%elm1,vovo%elm1,vvvo%elm1,ccsd_doubles%elm1,&
                             & eivalocc,eivalvirt,e4=e4,e5=e5,t1=ccsd_singles%elm1)
#endif

          else

             call ijk_loop_ser(nocc,nvirt,ovoo%elm1,vovo%elm1,vvvo%elm1,ccsd_doubles%elm1,&
                             & eivalocc,eivalvirt,e4=e4,e5=e5,t1=ccsd_singles%elm1)

          endif

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
    call time_start_phase(PHASE_WORK)

#endif
    ! release o^3v and v^3o integrals
    if (abc) then

       call tensor_free(vovv)
       call tensor_free(ooov)

    else

       call tensor_free(vvvo)
       call tensor_free(ovoo)

    endif

    ! now everything resides on the master...
    call tensor_free(ccsd_doubles)
    call tensor_free(vovo)

#ifdef VAR_MPI

    ! reduce singles and doubles arrays into that residing on the master
    reducing_to_master: if (nodtotal .gt. 1) then

       call time_start_phase(PHASE_COMM)

       if (print_frags) then

          call lsmpi_local_reduction(ccsdpt_singles%elm1,ccsdpt_singles%nelms,infpar%master)
          call lsmpi_local_reduction(ccsdpt_doubles%elm1,ccsdpt_doubles%nelms,infpar%master)

       else

          call lsmpi_local_reduction(e4,infpar%master)
          call lsmpi_local_reduction(e5,infpar%master)

       endif

       call time_start_phase(PHASE_WORK)

    end if reducing_to_master


    ! release stuff located on slaves
    releasing_the_slaves: if ((nodtotal .gt. 1) .and. .not. master) then


       ! release stuff initialized herein
       call mem_dealloc(eivalocc)
       call mem_dealloc(eivalvirt)


       ! now, release the slaves  
       return

    end if releasing_the_slaves

    call time_start_phase(PHASE_WORK)

#endif


    ! *************************************************
    ! ***** do canonical --> local transformation *****
    ! *************************************************

    if (print_frags) then

       if (abc) then

          call can_local_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,oovv=ccsdpt_doubles%elm1,ov=ccsdpt_singles%elm1)

       else

          call can_local_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,vvoo=ccsdpt_doubles%elm1,vo=ccsdpt_singles%elm1)

       endif

    else

       if (.not. DECinfo%pt_hack) then

          if (abc) then
  
             call can_local_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,ov=ccsd_singles%elm1)
   
          else
   
             call can_local_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,vo=ccsd_singles%elm1)
   
          endif

       endif

    endif

    ! now, release Uocc and Uvirt
    call tensor_free(Uvirt)
    call tensor_free(Uocc)

    ! clean up
    call mem_dealloc(eivalocc)
    call mem_dealloc(eivalvirt)

    if (master) call LSTIMER('CCSDPT_DRIVER (TOTAL)',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

  end subroutine ccsdpt_driver


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
    type(tensor), intent(in) :: ccsd_singles, ccsdpt_singles
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
    type(tensor), intent(in) :: ccsd_singles, ccsdpt_singles
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

  !> \brief Get MO integrals for ijk-CCSD(T) (in canonical basis), see integral storing order below.
  !> \author Janus Eriksen and Kasper Kristensen
  !> \date September-October 2012
  subroutine get_CCSDpT_integrals_ijk(MyLsitem,nbasis,nocc,nvirt,Cocc,Cvirt,ovoo,vvvo,tile_size)

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
    integer, intent(inout) :: tile_size
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
    integer, pointer      :: distribution(:)
    type(c_ptr)           :: distributionc
    integer(kind=ls_mpik) :: distributionw
    Character            :: intSpec(5)
    integer :: myload,first_el_i_block,nelms,tile_size_tmp,total_num_tiles,tile,ats1,ats2
    logical :: master, local
    integer(kind=long) :: o3v,v3
    real(realk), pointer :: dummy2(:)
    integer(kind=ls_mpik) :: mode,dest,nel2t, wi_idx, lg_me, nodtotal
    integer :: p,pos, batch, old_gammaB, dpos
    !> use background buffering to avoid memory fragmentation problems?
    logical :: use_bg_buf, first_round, dynamic_load
    call time_start_phase(PHASE_WORK)

    o3v          = nocc*nocc*nocc*nvirt
    v3           = nvirt**3
    use_bg_buf   = mem_is_background_buf_init()
    dynamic_load = DECinfo%dyn_load

#ifdef VAR_MPI

    nodtotal = infpar%lg_nodtot
    master = (infpar%lg_mynum .eq. infpar%master)
    lg_me  = infpar%lg_mynum
    if (master) call LSTIMER('START',tcpu,twall,DECinfo%output)

#else

    master = .true.
    lg_me  = 0
    call LSTIMER('START',tcpu,twall,DECinfo%output)

#endif

    if (DECinfo%pt_hack2) then

#ifdef VAR_MPI
       mode   = MPI_MODE_NOCHECK

       if (nodtotal .gt. 1) then

          ! Integrals (AI|KJ) in the order (J,A,I,K)
          dims = [nocc,nvirt,nocc,nocc]
          call tensor_init(ovoo, dims,4,bg=use_bg_buf)

          ! Integrals (AB|IC) in the order (C,B,A,I)
          dims = [nvirt,nvirt,nvirt,nocc]   
          call tensor_ainit(vvvo,dims,4,tdims=[nvirt,nvirt,nvirt,tile_size],atype="TDAR",bg=use_bg_buf)

          call tensor_random(ovoo)
          call tensor_random(vvvo)

          call tensor_scale(ovoo,1.0E-4_realk)
          call tensor_scale(vvvo,1.0E-4_realk)

       elseif (nodtotal .eq. 1) then

          dims = [nocc,nvirt,nocc,nocc]
          call tensor_init(ovoo, dims,4,bg=use_bg_buf)

          dims = [nvirt,nvirt,nvirt,nocc]
          call tensor_init(vvvo, dims,4,bg=use_bg_buf)

          call tensor_random(ovoo)
          call tensor_random(vvvo)

          call dscal(nvirt*nocc**3,1.0E-4_realk,ovoo%elm1,1)
          call dscal(nocc*nvirt**3,1.0E-4_realk,vvvo%elm1,1)

       endif
#else

       dims = [nocc,nvirt,nocc,nocc]
       call tensor_init(ovoo, dims,4,bg=use_bg_buf)

       dims = [nvirt,nvirt,nvirt,nocc]
       call tensor_init(vvvo, dims,4,bg=use_bg_buf)

       call tensor_random(ovoo)
       call tensor_random(vvvo)

       call dscal(nvirt*nocc**3,1.0E-4_realk,ovoo%elm1,1)
       call dscal(nocc*nvirt**3,1.0E-4_realk,vvvo%elm1,1)

#endif

       if (master) call LSTIMER('CCSD(T) INT (IJK - .PT_HACK2)',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

       return

    endif

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
    call tensor_init(ovoo, dims,4,bg=use_bg_buf)
    call tensor_zero(ovoo)

    ! Integrals (AB|IC) in the order (C,B,A,I)
    dims = [nvirt,nvirt,nvirt,nocc]
    local = .true.
#ifdef VAR_MPI
    if (infpar%lg_nodtot .gt. 1) then
       mode   = MPI_MODE_NOCHECK
       local = .false.
    endif
#endif

    call tensor_ainit(vvvo,dims,4,tdims=[nvirt,nvirt,nvirt,tile_size],atype="TDAR",local=local,bg=use_bg_buf)
    call tensor_zero(vvvo)

    ! For efficiency when calling dgemm, save transposed matrices
    if(use_bg_buf)then
       call mem_pseudo_alloc(CoccT,int(i8*nocc,kind=8),int(i8*nbasis,kind=8))
       call mem_pseudo_alloc(CvirtT,int(i8*nvirt,kind=8),int(i8*nbasis,kind=8))
    else
       call mem_alloc(CoccT,nocc,nbasis)
       call mem_alloc(CvirtT,nvirt,nbasis)
    endif
    call mat_transpose(nbasis,nocc,1.0E0_realk,Cocc,0.0E0_realk,CoccT)
    call mat_transpose(nbasis,nvirt,1.0E0_realk,Cvirt,0.0E0_realk,CvirtT)

    ! Determine optimal batchsizes and corresponding sizes of arrays
    call get_optimal_batch_sizes_ccsdpt_integrals(mylsitem,nbasis,nocc,nvirt,alphadim,gammadim,&
         & size1,size2,size3,.true.,.false.,tile_size)


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
    if(use_bg_buf)then
       call mem_pseudo_alloc(tmp1,size1)
       call mem_pseudo_alloc(tmp2,size2)
       call mem_pseudo_alloc(tmp3,size3)
    else
       call mem_alloc(tmp1,size1)
       call mem_alloc(tmp2,size2)
       call mem_alloc(tmp3,size3)
    endif

#ifdef VAR_MPI
    if(dynamic_load) then
       call mem_alloc( distribution, distributionc, 1)

       distribution = 0
       if(infpar%lg_mynum == 0) distribution(1) = infpar%lg_nodtot+1

       call lsmpi_win_create(distribution,distributionw,1,infpar%lg_comm)
#ifdef VAR_HAVE_MPI3
       call lsmpi_win_lock_all(distributionw,ass=mode)
#endif

    else


       ! alloc distribution array
       nullify(distribution)
       call mem_alloc(distribution,nbatchesGamma*nbatchesAlpha)
       ! init distribution
       distribution = 0
       myload = 0

       if (infpar%lg_nodtot .gt. 1) then

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

       endif
    endif

#endif

    ! Start looping over gamma and alpha batches and calculate integrals
    ! ******************************************************************

    old_gammaB = -1
    batch      =  0

    first_round=.false.
    if(dynamic_load)then
       first_round = .true.
       batch       = lg_me + 1
    endif


    BatchLoop: do while(batch <= nbatchesGamma*nbatchesAlpha) ! AO batches


       !check if the current job is to be done by current node
       call check_job(batch,first_round,dynamic_load,alphaB,gammaB,nbatchesAlpha,&
          &nbatchesGamma,distribution, distributionw, DECinfo%PL>2)

       !exit the loop
       if(batch > nbatchesGamma*nbatchesAlpha ) exit BatchLoop


       ! If the new gamma is different form the old gamma batch
       NewGammaBatch: if( gammaB /= old_gammaB )then
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
       endif NewGammaBatch

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

       if(DECinfo%ccsolverskip)then
          if(first_round)then
          call random_number(tmp1)
          call random_number(tmp2)
          call random_number(tmp3)
          call random_number(ovoo%elm1)
          first_round=.false.
          endif
       else

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

#ifdef VAR_MPI

       if (infpar%lg_nodtot .gt. 1) then

          m = nvirt**3
          k = dimAlpha
          total_num_tiles = vvvo%ntiles
          tile = 0

          ! adapt comment to IJK scheme!!!
          ! mpi   : 1) tmp2(B,I,A,tile) = sum_{alpha in alphaB} tmp1(B,I,A,alpha) Cvirt(alpha,tile)
          !         2) vovv(B,I,A,C) += sum_{tile in CB} tmp2(B,I,A,tile)
          ! serial: vovv(B,I,A,C) += sum_{alpha in alphaB} tmp1(B,I,A,alpha) Cvirt(alpha,C)

          ! reorder tmp1 and do vvvo(B,A,C,I) += sum_{i in IB} tmp1(B,A,C,i)
          do i=1,nocc,tile_size

             tile = tile + 1

             call get_tile_dim(nelms,vvvo,i8*tile)
             tile_size_tmp = nelms/((i8*nvirt)*(i8*nvirt**2))

             n = tile_size_tmp

             ! tmp1(C,A,B,tile) = sum_{alpha in alphaB} tmp3(C,A,B,alpha) Cocc(alpha,tile)
             if( n == 1)then
                call dgemv('N',m,k,1.0E0_realk,tmp3,m,Cocc(AlphaStart,i),1,0.0E0_realk,tmp1,1)
             else
                call dgemm('N','N',m,n,k,1.0E0_realk,tmp3,m,Cocc(AlphaStart,i),nbasis,0.0E0_realk,tmp1,m)
             endif

             ! *** tmp1 corresponds to (AB|iC) in Mulliken notation. Noting that the v³o integrals
             ! are normally written as g_{AIBC}, we may also write this Mulliken integral (with substitution
             ! of dummy indices A=B, B=C, and C=A) as (BC|IA). In order to align with the vvvo order of
             ! ccsd(t) driver routine, we reorder as:
             ! (BC|IA) --> (CB|AI), i.e., tmp1(C,A,B,tile) = ABCI(A,B,C,tile) (norm. notat.) --> 
             !                                            tmp1(C,B,A,tile) (norm. notat.) = tmp1(B,A,C,tile) (notat. herein)
             ! 
             ! next, we accumulate
             ! vvvo(B,A,C,I) += sum_{tile in IB} tmp1(B,A,C,tile)

             call array_reorder_4d(1.0E0_realk,tmp1,nvirt,nvirt,nvirt,tile_size_tmp,[3,2,1,4],0.0E0_realk,tmp2)

             call time_start_phase(PHASE_COMM)
#ifdef VAR_HAVE_MPI3
             call tensor_lock_win(vvvo,tile,'s')
#endif

             call get_residence_of_tile(vvvo,tile,dest, dpos, pos, wi_idx)

             p = pos - 1

             do first_el_i_block=1,v3*tile_size_tmp,MAX_SIZE_ONE_SIDED
#ifndef VAR_HAVE_MPI3
                call tensor_lock_win(vvvo,tile,'s',assert=mode)
#endif
                nel2t=MAX_SIZE_ONE_SIDED
                if(((v3*tile_size_tmp-first_el_i_block)<MAX_SIZE_ONE_SIDED).and.&
                   &(mod(v3*tile_size_tmp-first_el_i_block+1,i8*MAX_SIZE_ONE_SIDED)/=0))&
                   &nel2t=int(mod(v3*tile_size_tmp,i8*MAX_SIZE_ONE_SIDED),kind=ls_mpik)


                call lsmpi_acc(tmp2(first_el_i_block:first_el_i_block+nel2t-1),nel2t,p+first_el_i_block,dest,vvvo%wi(wi_idx))

#ifdef VAR_HAVE_MPI3
                call lsmpi_win_flush(vvvo%wi(wi_idx),rank=dest,local=.true.)
#else
                call tensor_unlock_win(vvvo,tile)
#endif
             enddo

#ifdef VAR_HAVE_MPI3
             call tensor_unlock_win(vvvo,tile)
#endif
             call time_start_phase(PHASE_WORK)

          end do

       else

          do i=1,nocc

             m = nvirt**3
             k = dimAlpha
             ! for description, see mpi section above
             call dgemv('N',m,k,1.0E0_realk,tmp3,m,Cocc(AlphaStart,i),1,0.0E0_realk,tmp1,1)

             call array_reorder_3d(1.0E0_realk,tmp1,nvirt,nvirt,nvirt,[3,2,1],1.0E0_realk,vvvo%elm4(:,:,:,i))

          end do

       endif

#else

       do i=1,nocc

          m = nvirt**3
          k = dimAlpha
          ! for description, see mpi section above
          call dgemv('N',m,k,1.0E0_realk,tmp3,m,Cocc(AlphaStart,i),1,0.0E0_realk,tmp1,1)

          call array_reorder_3d(1.0E0_realk,tmp1,nvirt,nvirt,nvirt,[3,2,1],1.0E0_realk,vvvo%elm4(:,:,:,i))

       end do

#endif
       endif

    end do BatchLoop

#ifdef VAR_MPI

    if (infpar%lg_nodtot .gt. 1) then

#ifdef VAR_PTR_RESHAPE
       dummy2(1:(i8*nocc*nvirt)*nocc*nocc) => ovoo%elm1(:)
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
    if(dynamic_load)then
#ifdef VAR_HAVE_MPI3
       call lsmpi_win_unlock_all(distributionw)
#endif
       call lsmpi_win_free(distributionw)
       call mem_dealloc(distribution,distributionc)
    else
       call mem_dealloc(distribution)
    endif

#endif

    if(DECinfo%ccsolverskip)then
       ats1 = vvvo%access_type
       ats2 = ovoo%access_type
       vvvo%access_type = AT_ALL_ACCESS
       ovoo%access_type = AT_ALL_ACCESS
       call tensor_random(vvvo)
       call tensor_random(ovoo)
       vvvo%access_type = ats1
       ovoo%access_type = ats2
    endif

    ! free stuff
    ! **********
    if( use_bg_buf)then
       call mem_pseudo_dealloc(tmp3)
       call mem_pseudo_dealloc(tmp2)
       call mem_pseudo_dealloc(tmp1)
       call mem_pseudo_dealloc(CvirtT)
       call mem_pseudo_dealloc(CoccT)
    else
       call mem_dealloc(tmp3)
       call mem_dealloc(tmp2)
       call mem_dealloc(tmp1)
       call mem_dealloc(CvirtT)
       call mem_dealloc(CoccT)
    endif
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
    logical :: master,local
    integer(kind=long) :: o3v,v3,ov2
    real(realk), pointer :: dummy2(:)
    integer(kind=ls_mpik) :: mode,dest,nel2t, wi_idx, nodtotal
    integer :: p,pos,dpos
    !> use background buffering to avoid memory fragmentation problems?
    logical :: use_bg_buf
    call time_start_phase(PHASE_WORK)

    o3v        = nocc*nocc*nocc*nvirt
    v3         = nvirt**3
    ov2        = nocc*nvirt**2
    use_bg_buf = mem_is_background_buf_init()

#ifdef VAR_MPI

    nodtotal = infpar%lg_nodtot
    master = (infpar%lg_mynum .eq. infpar%master)
    if (master) call LSTIMER('START',tcpu,twall,DECinfo%output)

#else

    master = .true.
    call LSTIMER('START',tcpu,twall,DECinfo%output)

#endif

    if (DECinfo%pt_hack2) then

#ifdef VAR_MPI
       mode   = MPI_MODE_NOCHECK

       if (nodtotal .gt. 1) then

          ! ooov: Integrals (AI|KJ) in the order (I,J,K,A)
          dims = [nocc,nocc,nocc,nvirt]
          call tensor_init(ooov, dims,4,bg=use_bg_buf)

          ! vovv: Integrals (AB|IC) in the order (B,I,A,C)
          dims = [nvirt,nocc,nvirt,nvirt]
          call tensor_ainit(vovv,dims,4,tdims=[nvirt,nocc,nvirt,tile_size],atype="TDAR",bg=use_bg_buf)

          call tensor_random(ooov)
          call tensor_random(vovv)

          call tensor_scale(ooov,1.0E-4_realk)
          call tensor_scale(vovv,1.0E-4_realk)

       elseif (nodtotal .eq. 1) then

          ! ooov: Integrals (AI|KJ) in the order (I,J,K,A)
          dims = [nocc,nocc,nocc,nvirt]
          call tensor_init(ooov, dims,4,bg=use_bg_buf)

          ! vovv: Integrals (AB|IC) in the order (B,I,A,C)
          dims = [nvirt,nocc,nvirt,nvirt]
          call tensor_init(vovv,dims,4,bg=use_bg_buf)

          call tensor_random(ooov)
          call tensor_random(vovv)

          call dscal(nvirt*nocc**3,1.0E-4_realk,ooov%elm1,1)
          call dscal(nocc*nvirt**3,1.0E-4_realk,vovv%elm1,1)

       endif
#else

       ! ooov: Integrals (AI|KJ) in the order (I,J,K,A)
       dims = [nocc,nocc,nocc,nvirt]
       call tensor_init(ooov, dims,4,bg=use_bg_buf)

       ! vovv: Integrals (AB|IC) in the order (B,I,A,C)
       dims = [nvirt,nocc,nvirt,nvirt]
       call tensor_init(vovv,dims,4,bg=use_bg_buf)

       call tensor_random(ooov)
       call tensor_random(vovv)

       call dscal(nvirt*nocc**3,1.0E-4_realk,ooov%elm1,1)
       call dscal(nocc*nvirt**3,1.0E-4_realk,vovv%elm1,1)

#endif

       if (master) call LSTIMER('CCSD(T) INT (ABC - .PT_HACK2)',tcpu,twall,DECinfo%output,FORCEPRINT=.true.)

       return

    endif

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
    call tensor_init(ooov, dims,4,bg=use_bg_buf)
    call tensor_zero(ooov)

    ! vovv: Integrals (AB|IC) in the order (B,I,A,C)
    dims  = [nvirt,nocc,nvirt,nvirt]
    local = .true.

#ifdef VAR_MPI
    if (infpar%lg_nodtot .gt. 1) then
       mode   = MPI_MODE_NOCHECK
       local  = .false.
    endif
#endif

    call tensor_ainit(vovv,dims,4,tdims=[nvirt,nocc,nvirt,tile_size],atype="TDAR",local=local,bg=use_bg_buf)
    call tensor_zero(vovv)

    ! For efficiency when calling dgemm, save transposed matrices
    if(use_bg_buf)then
       call mem_pseudo_alloc(CoccT,int(i8*nocc,kind=8),int(i8*nbasis,kind=8))
       call mem_pseudo_alloc(CvirtT,int(i8*nvirt,kind=8),int(i8*nbasis,kind=8))
    else
       call mem_alloc(CoccT,nocc,nbasis)
       call mem_alloc(CvirtT,nvirt,nbasis)
    endif
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
    if( use_bg_buf )then
       call mem_pseudo_alloc(tmp1,size1)
       call mem_pseudo_alloc(tmp2,size2)
       call mem_pseudo_alloc(tmp3,size3)
    else
       call mem_alloc(tmp1,size1)
       call mem_alloc(tmp2,size2)
       call mem_alloc(tmp3,size3)
    endif

#ifdef VAR_MPI

    if (infpar%lg_nodtot .gt. 1) then

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

          if (infpar%lg_nodtot .gt. 1) then

             ! distribute tasks
             if (distribution((alphaB-1)*nbatchesGamma+gammaB) .ne. infpar%lg_mynum) then
   
                cycle BatchAlpha
   
             end if
   
             if(DECinfo%PL>2)write (*, '("Rank(T) ",I3," starting job (",I3,"/",I3,",",I3,"/",I3,")")')&
                &infpar%lg_mynum,alphaB,nbatchesAlpha,gammaB,nbatchesGamma

          endif

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

          if (infpar%lg_nodtot .gt. 1) then

             ! mpi   : 1) tmp2(B,I,A,tile) = sum_{alpha in alphaB} tmp1(B,I,A,alpha) Cvirt(alpha,tile)
             !         2) vovv(B,I,A,C) += sum_{tile in CB} tmp2(B,I,A,tile)
             ! serial: vovv(B,I,A,C) += sum_{alpha in alphaB} tmp1(B,I,A,alpha) Cvirt(alpha,C)
             m = nocc*nvirt**2
             k = dimAlpha
   
             total_num_tiles = vovv%ntiles
             tile = 0
   
             do c=1,nvirt,tile_size
   
                tile = tile + 1
   
                call get_tile_dim(nelms,vovv,i8*tile)
                tile_size_tmp = nelms/(nocc*(i8*nvirt**2))
   
                n = tile_size_tmp
   
                ! tmp2(B,I,A,tile) = sum_{alpha in alphaB} tmp1(B,I,A,alpha) Cvirt(alpha,tile)
                call dgemm('N','N',m,n,k,1.0E0_realk,tmp1,m,Cvirt(AlphaStart,c),nbasis,0.0E0_realk,tmp2,m)
   
                call time_start_phase(PHASE_COMM)
#ifdef VAR_HAVE_MPI3
                call tensor_lock_win(vovv,tile,'s',assert=mode)
#endif
                call get_residence_of_tile(vovv,tile,dest, dpos, pos, wi_idx)

                p      = pos - 1
   
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

          elseif (infpar%lg_nodtot .eq. 1) then 

             m = nocc*nvirt**2
             k = dimAlpha
             n = nvirt
             call dgemm('N','N',m,n,k,1.0E0_realk,tmp1,m,Cvirt(AlphaStart,1),nbasis,1.0E0_realk,vovv%elm1,m)

          endif

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
       dummy2(1:(i8*nocc*nocc)*nocc*nvirt) => ooov%elm1(:)
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

       ! dealloc distribution array
       call mem_dealloc(distribution)

    end if

#endif

    ! free stuff
    ! **********
    if( use_bg_buf )then
       call mem_pseudo_dealloc(tmp3)
       call mem_pseudo_dealloc(tmp2)
       call mem_pseudo_dealloc(tmp1)
       call mem_pseudo_dealloc(CvirtT)
       call mem_pseudo_dealloc(CoccT)
    else
       call mem_dealloc(tmp3)
       call mem_dealloc(tmp2)
       call mem_dealloc(tmp1)
       call mem_dealloc(CvirtT)
       call mem_dealloc(CoccT)
    endif

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
     !> tile_size
     integer, intent(in) :: tile_size
     !> memory reals
     real(realk) :: MemoryNeeded, MemoryAvailable
     integer :: MaxAObatch, MinAOBatch, AlphaOpt, GammaOpt,alpha,gamma,iAO
     integer(kind=ls_mpik) :: nnod,me
     logical :: master, use_bg_buf
     ! Memory currently available
     ! **************************
     use_bg_buf = mem_is_background_buf_init()
     if(use_bg_buf)then
        MemoryAvailable = (mem_get_bg_buf_free()*8.0E0_realk)/1024.0E0_realk**3
     else
        call get_currently_available_memory(MemoryAvailable)
     endif
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
    !> tle_size
    integer, intent(in) :: tile_size 
    real(realk) :: GB
    integer(kind=long) :: tmpI
    GB = 1073741824.0E0_realk ! 1 GB
    ! Array sizes needed in get_CCSDpT_integrals are checked and the largest one is found
 
    ! Tmp array 1
    if (abc) then

       size1 = (i8*alphadim*gammadim)*(i8*nbasis**2)
       tmpI = (i8*alphadim*gammadim)*(i8*nocc*nvirt)
       size1 = max(size1,tmpI)
! this one is big
       tmpI = (i8*alphadim*nocc)*(i8*nvirt**2)
       size1 = max(size1,tmpI)

    else

       size1 = (i8*alphadim*gammadim)*(i8*nbasis**2)
       tmpI = (i8*nvirt**2)*(i8*gammadim*alphadim)
       size1 = max(size1,tmpI)
       tmpI = (i8*nvirt*nocc)*(i8*gammadim*alphadim)
       size1 = max(size1,tmpI)
       tmpI = (i8*nvirt)*(i8*nocc**2*alphadim)
       size1 = max(size1,tmpI)
#ifdef VAR_MPI
       tmpI = (i8*nvirt)*(i8*nvirt**2)*tile_size
       size1 = max(size1,tmpI)
#else
       tmpI = (i8*nvirt)*(i8*nvirt**2)
       size1 = max(size1,tmpI)
#endif

    endif
  
    ! tmp array 2
    if (abc) then

       size2 = (i8*alphadim*gammadim)*(i8*nbasis*nocc)
       tmpI = (i8*alphadim*gammadim)*(i8*nocc**2)
       size2 = max(size2,tmpI)
       tmpI = (i8*alphadim*gammadim)*(i8*nocc*nvirt)
       size2 = max(size2,tmpI)
#ifdef VAR_MPI
       tmpI = (i8*nvirt**2)*(i8*nocc*tile_size)
       size2 = max(size2,tmpI)
#endif

    else

       size2 = (i8*alphadim*gammadim)*(i8*nbasis*nvirt)
       tmpI = (i8*alphadim*gammadim)*(i8*nvirt*nocc)
       size2 = max(size2,tmpI)
#ifdef VAR_MPI
       tmpI = (i8*nvirt)*(i8*nvirt**2)*tile_size
       size2 = max(size2,tmpI)
#else
       tmpI = (i8*nvirt)*(i8*nvirt**2)
       size2 = max(size2,tmpI)
#endif

    endif
  
    ! Tmp array3
    if (abc) then

       size3 = (i8*alphadim*gammadim)*nocc**2
       tmpI = alphadim*(i8*nocc**3)
       size3 = max(size3,tmpI)
 
    else

       size3 = (i8*alphadim*gammadim)*(i8*nvirt**2)
! this one is big
       tmpI = (i8*alphadim*nvirt)*(i8*nvirt**2)
       size3 = max(size3,tmpI)

    endif

    ! Size = size1+size2+size3,  convert to GB
    mem = realk*(size1+size2+size3)/GB


  end subroutine get_max_arraysizes_for_ccsdpt_integrals

  subroutine ccsdpt_info(nbasis,nocc,nvirt,print_frags,abc,ijk_nbuffs,abc_nbuffs,ijk_tile_size,abc_tile_size,nodtotal)

      implicit none

      integer, intent(in) :: nbasis,nocc,nvirt
      logical, intent(in) :: print_frags,abc
      integer, intent(in) :: nodtotal
      logical :: ijk,manual_ijk_1,manual_ijk_2,manual_abc_1,manual_abc_2,gpu
      integer, intent(inout) :: ijk_nbuffs,abc_nbuffs,ijk_tile_size,abc_tile_size
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
      logical :: use_bg_buf

      ijk_nbuffs = 0
      abc_nbuffs = 0
      ijk_tile_size = 0
      abc_tile_size = 0
      manual_abc_1 = .false.; manual_abc_2 = .false.
      manual_ijk_1 = .false.; manual_ijk_2 = .false.
      gpu = .false.
      num_gpu = 0
      free_cpu = 0.0E0_realk
      total_gpu = 0
      free_gpu = 0
      acc_async = .true.

      if (DECinfo%acc_sync) acc_async = .false.

      use_bg_buf = mem_is_background_buf_init()
      if(use_bg_buf)then
         free_cpu = (mem_get_bg_buf_free()*8.0E0_realk)/1024.0E0_realk**3
      else
         call get_currently_available_memory(free_cpu)
      endif

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
    
                manual_ijk_1 = .true.
                ijk_nbuffs = DECinfo%ijk_nbuffs
   
            endif
   
            if (manual_ijk_1) then
   
               if (ijk_nbuffs .lt. 1) call lsquit('manually set ijk_nbuffs (NBUFFS_IJK) .lt. 1 - aborting...',DECinfo%output)
   
            else
   
!               ! here; determine ijk_nbuffs based on available cpu/gpu memory
!               if (ijk_nbuffs .eq. ijk_default) call new_ijk_nbuffs_routine()
               ijk_nbuffs = 2
   
            endif
   
            if (DECinfo%ijk_tile_size .lt. ijk_default) then

               manual_ijk_2 = .true.
               ijk_tile_size = DECinfo%ijk_tile_size

            endif

            if (manual_ijk_2) then

               if (ijk_tile_size .lt. 1) call lsquit('manually set tile size (.IJK_TILE) .lt. 1 - aborting...',DECinfo%output)

            else

               call ijk_tile_size_routine(nocc,nvirt,print_frags,free_cpu,nodtotal,ijk_nbuffs,ijk_tile_size)

            endif

            if (ijk_tile_size .gt. nocc) call lsquit('manually set tile size (.IJK_TILE) .gt. nocc - aborting...',DECinfo%output)

         else ! abc == .true.
   
            if (DECinfo%abc_nbuffs .lt. abc_default) then
   
                manual_abc_1 = .true.
                abc_nbuffs = DECinfo%abc_nbuffs
   
            endif
   
            if (manual_abc_1) then
   
               if (abc_nbuffs .lt. 1) call lsquit('manually set abc_nbuffs (NBUFFS_ABC) .lt. 1 - aborting...',DECinfo%output)
   
            else
   
!               ! here; determine ijk_nbuffs based on available cpu/gpu memory
!               if (abc_nbuffs .eq. abc_default) call new_abc_nbuffs_routine()
               abc_nbuffs = 2
   
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
            if (manual_ijk_1) then
               write(DECinfo%output,'(a,i4)')     'Input # IJK buffers    = ',ijk_nbuffs
            else
               write(DECinfo%output,'(a,i4)')     'Number of IJK buffers  = ',ijk_nbuffs
            endif
            if (manual_ijk_2) then
               write(DECinfo%output,'(a,i4)')     'Input IJK tile size    = ',ijk_tile_size
            else
               write(DECinfo%output,'(a,i4)')     'IJK tile size          = ',ijk_tile_size
            endif
         endif
         write(DECinfo%output,'(a,l4)')     'ABC partitioning       = ',abc
         if (abc) then
            if (manual_abc_1) then
               write(DECinfo%output,'(a,i4)')     'Input # ABC buffers    = ',abc_nbuffs
            else
               write(DECinfo%output,'(a,i4)')     'Number of ABC buffers  = ',abc_nbuffs
            endif
            if (manual_abc_2) then
               write(DECinfo%output,'(a,i4)')     'Input ABC tile size    = ',abc_tile_size
            else
               write(DECinfo%output,'(a,i4)')     'ABC tile size          = ',abc_tile_size
            endif
         endif
         write(DECinfo%output,'(a,l4)')     'Single. prec. calc.?   = ',DECinfo%pt_single_prec
         write(DECinfo%output,'(a,g11.4)')     'Free CPU memory (GB)   = ',free_cpu
         write(DECinfo%output,'(a,l4)')     'Are we using GPUs?     = ',gpu
         if (gpu) then
            write(DECinfo%output,'(a,l4)')     'Asynchronous OpenACC?  = ',acc_async
            write(DECinfo%output,'(a,i4)')     'Number of GPUs         = ',num_gpu
            write(DECinfo%output,'(a,g11.4)')     'Total GPU memory (GB)  = ',total_gpu / gb
            write(DECinfo%output,'(a,g11.4)')     'Free GPU memory (GB)   = ',free_gpu / gb
         endif
         write(DECinfo%output,*)
         write(DECinfo%output,*)

      endif

  end subroutine ccsdpt_info


  subroutine ijk_tile_size_routine(nocc,nvirt,print_frags,free_cpu,nodtotal,ijk_nbuffs,ijk_tile_size)

      implicit none

      integer, intent(in) :: nocc,nvirt,nodtotal,ijk_nbuffs
      real(realk), intent(in) :: free_cpu
      logical, intent(in) :: print_frags
      integer, intent(inout) :: ijk_tile_size
      real(realk) :: mem_avail_start,mem_accum_tmp,mem_est_avail_tmp,mem_est_avail
      real(realk) :: mem_vvvo_pdm,mem_vvvo_local,mem_vvoo_pdm,mem_vvoo_local,mem_ccsd_pdm,mem_ccsd_local
      integer(kind=long) :: ovoo,ccsd_total,vvoo_total,vvvo_total
      integer(kind=long) :: eivalocc,eivalvirt,trip_ampls,ccsdpt_singles,ccsdpt_doubles
      integer(kind=long) :: vvvo_pdm,vvvo_local,vvoo_pdm,vvoo_local,ccsd_pdm,ccsd_local
      integer(kind=long) :: mem_int_tmp
      integer(kind=long) :: max_ijk_tile_size,ts
      integer :: remainder_1,remainder_2,num_tiles_tot,num_tiles_node
      real(realk), parameter :: GB = 1073741824.0E0_realk ! 1 GB = 1024.**3

      ! available memory - note that we multiply by 95 % to be on the safe side!
      mem_avail_start = 0.95 * free_cpu

      ! total distributed ccsd doubles ampls
      ccsd_total = (i8*nocc**2)*(i8*nvirt**2)

      ! total distributed integrals
      vvvo_total = (i8*nocc*nvirt)*(i8*nvirt**2)
      vvoo_total = (i8*nocc**2)*(i8*nvirt**2)

      ! local integrals
      ovoo = (i8*nocc**2)*(i8*nocc*nvirt)

      ! orbital energies
      eivalocc = i8*nocc
      eivalvirt = i8*nvirt

      ! triples amplitudes and temp array
      trip_ampls = 2*(i8*nvirt)*(i8*nvirt**2)

      ! ccsdpt intermediates
      ccsdpt_singles = i8*nocc*nvirt
      ccsdpt_doubles = (i8*nocc**2)*(i8*nvirt**2)

      ! temp sum of integer elements
      mem_int_tmp = ovoo + &
                  & eivalocc + &
                  & eivalvirt + &
                  & trip_ampls

      if (print_frags) mem_int_tmp = mem_int_tmp + ccsdpt_singles + ccsdpt_doubles

      ! estimate the memory needed for allocations
      mem_accum_tmp = realk*mem_int_tmp / GB

      ! estimate available memory AFTER allocations, but BEFORE ccsd t2, vvoo, and vvvo allocations
      mem_est_avail_tmp = mem_avail_start - mem_accum_tmp
      if (mem_est_avail_tmp .lt. 0.0E0_realk) call lsquit('mem_est_avail_tmp .lt. 0 GB (IJK) - aborting...',DECinfo%output)

      ! max tile_size
      max_ijk_tile_size = int(nocc / nodtotal)

      do ts = max_ijk_tile_size,1,-1

         ! how many tiles are there in total?
         num_tiles_tot = int(nocc / ts)

         ! modulo_1
         remainder_1 = mod(nocc,ts)

         ! update total number of tiles
         if (remainder_1 .gt. 0) num_tiles_tot = num_tiles_tot + 1

         ! how many tiles per node in PDM?
         num_tiles_node = int(num_tiles_tot / nodtotal)

         ! modulo_2
         remainder_2 = mod(num_tiles_tot,nodtotal)

         ! update number of tiles per node
         if (remainder_2 .gt. 0) num_tiles_node = num_tiles_node + 1

         ! calculate the PDM memory requirements for the given tile_size
         vvvo_pdm = num_tiles_node*(i8*nvirt)*(i8*nvirt**2)*ts
         mem_vvvo_pdm = realk*vvvo_pdm / GB
         vvoo_pdm = num_tiles_node*(i8*nvirt**2)*ts**2
         mem_vvoo_pdm = realk*vvoo_pdm / GB
         ccsd_pdm = num_tiles_node*nocc*(i8*nvirt**2)*ts
         mem_ccsd_pdm = realk*ccsd_pdm / GB

         ! calculate the local memory requirements for the given tile_size and the given number of tiles
         vvvo_local = 3*ijk_nbuffs*(i8*nvirt)*(i8*nvirt**2)*ts
         mem_vvvo_local = realk*vvvo_local / GB
         vvoo_local = 6*ijk_nbuffs*(i8*nvirt**2)*ts**2
         mem_vvoo_local = realk*vvoo_local / GB
         ccsd_local = 3*ijk_nbuffs*nocc*(i8*nvirt**2)*ts
         mem_ccsd_local = realk*ccsd_local / GB

         ! estimate available memory AFTER pdm allocation
         mem_est_avail = mem_est_avail_tmp - &
                       & ((mem_vvvo_pdm + mem_vvvo_local) + &
                       & (mem_vvoo_pdm + mem_vvoo_local) + &
                       & (mem_ccsd_pdm + mem_ccsd_local))

         if ((mem_est_avail .lt. 0.0E0_realk)) then

            if (ts .eq. 1) then

               print *,'ts                = ',ts
               print *,'num_tiles_tot     = ',num_tiles_tot
               print *,'num_tiles_node    = ',num_tiles_node
               print *,'mem_avail_start   = ',mem_avail_start
               print *,'mem_est_avail_tmp = ',mem_est_avail_tmp
               print *,'mem_vvvo_pdm      = ',mem_vvvo_pdm
               print *,'mem_vvvo_local    = ',mem_vvvo_local
               print *,'mem_vvoo_pdm      = ',mem_vvoo_pdm
               print *,'mem_vvoo_local    = ',mem_vvoo_local
               print *,'mem_ccsd_pdm      = ',mem_ccsd_pdm
               print *,'mem_ccsd_local    = ',mem_ccsd_local
               print *,'mem_est_avail     = ',mem_est_avail

               call lsquit('mem_est_avail .lt. 0 GB for smallest possible ijk_tile_size - aborting...',DECinfo%output)

            endif

            cycle

         else

            ijk_tile_size = ts

            return

         endif

      enddo

  end subroutine ijk_tile_size_routine


  subroutine abc_tile_size_routine(nocc,nvirt,print_frags,free_cpu,nodtotal,abc_nbuffs,abc_tile_size)

      implicit none

      integer, intent(in) :: nocc,nvirt,nodtotal,abc_nbuffs
      real(realk), intent(in) :: free_cpu
      logical, intent(in) :: print_frags
      integer, intent(inout) :: abc_tile_size
      real(realk) :: mem_avail_start,mem_accum_tmp,mem_est_avail_tmp,mem_est_avail
      real(realk) :: mem_vovv_pdm,mem_vovv_local,mem_oovv_pdm,mem_oovv_local,mem_ccsd_pdm,mem_ccsd_local
      integer(kind=long) :: ooov,ccsd_total,oovv_total,vovv_total
      integer(kind=long) :: eivalocc,eivalvirt,trip_ampls,ccsdpt_singles,ccsdpt_doubles
      integer(kind=long) :: vovv_pdm,vovv_local,oovv_pdm,oovv_local,ccsd_pdm,ccsd_local
      integer(kind=long) :: mem_int_tmp
      integer(kind=long) :: max_abc_tile_size,ts
      integer :: remainder_1,remainder_2,num_tiles_tot,num_tiles_node
      real(realk), parameter :: GB = 1073741824.0E0_realk ! 1 GB = 1024.**3

      ! available memory - note that we multiply by 95 % to be on the safe side!
      mem_avail_start = 0.95 * free_cpu

      ! total distributed ccsd doubles ampls
      ccsd_total = (i8*nocc**2)*(i8*nvirt**2)

      ! total distributed integrals
      vovv_total = (i8*nocc*nvirt)*(i8*nvirt**2)
      oovv_total = (i8*nocc**2)*(i8*nvirt**2)

      ! local integrals
      ooov = (i8*nocc**3)*nvirt

      ! orbital energies
      eivalocc = i8*nocc
      eivalvirt = i8*nvirt

      ! triples amplitudes and temp array
      trip_ampls = 2*(i8*nocc**3)

      ! ccsdpt intermediates
      ccsdpt_singles = i8*nocc*nvirt
      ccsdpt_doubles = (i8*nocc**2)*(i8*nvirt**2)

      ! temp sum of integer elements
      mem_int_tmp = ooov + &
                  & eivalocc + &
                  & eivalvirt + &
                  & trip_ampls

      if (print_frags) mem_int_tmp = mem_int_tmp + ccsdpt_singles + ccsdpt_doubles 

      ! estimate the memory needed for allocations
      mem_accum_tmp = realk*mem_int_tmp / GB

      ! estimate available memory AFTER allocations, but BEFORE ccsd t2, oovv, and vovv allocations
      mem_est_avail_tmp = mem_avail_start - mem_accum_tmp
      if (mem_est_avail_tmp .lt. 0.0E0_realk) call lsquit('mem_est_avail_tmp .lt. 0 GB (ABC) - aborting...',DECinfo%output) 

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
         vovv_pdm = num_tiles_node*nocc*(i8*nvirt**2)*ts
         mem_vovv_pdm = realk*vovv_pdm / GB
         oovv_pdm = num_tiles_node*(i8*nocc**2)*ts**2
         mem_oovv_pdm = realk*oovv_pdm / GB
         ccsd_pdm = num_tiles_node*(i8*nvirt*nocc**2)*ts
         mem_ccsd_pdm = realk*ccsd_pdm / GB

         ! calculate the local memory requirements for the given tile_size and the given number of tiles
         vovv_local = 3*abc_nbuffs*nocc*(i8*nvirt**2)*ts
         mem_vovv_local = realk*vovv_local / GB
         oovv_local = 6*abc_nbuffs*(i8*nocc**2)*ts**2
         mem_oovv_local = realk*oovv_local / GB
         ccsd_local = 3*abc_nbuffs*(i8*nvirt*nocc**2)*ts
         mem_ccsd_local = realk*ccsd_local / GB

         ! estimate available memory AFTER pdm allocation
         mem_est_avail = mem_est_avail_tmp - & 
                       & ((mem_vovv_pdm + mem_vovv_local) + &
                       & (mem_oovv_pdm + mem_oovv_local) + &
                       & (mem_ccsd_pdm + mem_ccsd_local))

         if ((mem_est_avail .lt. 0.0E0_realk)) then

            if (ts .eq. 1) then

               print *,'ts                = ',ts
               print *,'num_tiles_tot     = ',num_tiles_tot
               print *,'num_tiles_node    = ',num_tiles_node
               print *,'mem_avail_start   = ',mem_avail_start
               print *,'mem_est_avail_tmp = ',mem_est_avail_tmp
               print *,'mem_vovv_pdm      = ',mem_vovv_pdm
               print *,'mem_vovv_local    = ',mem_vovv_local
               print *,'mem_oovv_pdm      = ',mem_oovv_pdm
               print *,'mem_oovv_local    = ',mem_oovv_local
               print *,'mem_ccsd_pdm      = ',mem_ccsd_pdm
               print *,'mem_ccsd_local    = ',mem_ccsd_local
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

  subroutine convert_ccsd_and_vovo(nocc,nvirt,nbasis,Uocc,Uvirt,ccsd_doubles_in,vovo_in,ccsd_doubles,vovo,&
                                 & nodtotal,abc,ijk_tile_size,abc_tile_size,use_bg)

    implicit none

    integer, intent(in) :: nocc,nvirt,nbasis,nodtotal
    logical, intent(in) :: abc,use_bg
    integer, intent(in) :: ijk_tile_size,abc_tile_size
    type(tensor), intent(inout) :: Uocc,Uvirt
    type(tensor), intent(inout) :: ccsd_doubles_in,vovo_in
    type(tensor), intent(inout) :: ccsd_doubles,vovo
    ! tmp tensors
    type(tensor) :: tmp_tensor_1,tmp_tensor_2
    logical :: transform_bg

    if (DECinfo%pt_hack) then

#ifdef VAR_MPI
       if (nodtotal .gt. 1) then

          if (abc) then

             call tensor_minit(vovo,[nocc,nocc,nvirt,nvirt],4,&
                &tdims=[nocc,nocc,abc_tile_size,abc_tile_size],atype='TDAR',bg=use_bg)
             call tensor_minit(ccsd_doubles,[nocc,nocc,nvirt,nvirt],4,&
                &tdims=[nocc,nocc,nvirt,abc_tile_size],atype='TDAR',bg=use_bg)

          else

             call tensor_minit(vovo,[nvirt,nvirt,nocc,nocc],4,&
                &tdims=[nvirt,nvirt,ijk_tile_size,ijk_tile_size],atype='TDAR',bg=use_bg)
             call tensor_minit(ccsd_doubles,[nvirt,nvirt,nocc,nocc],4,&
                &tdims=[nvirt,nvirt,nocc,ijk_tile_size],atype='TDAR',bg=use_bg)

          endif

          call tensor_random(ccsd_doubles)
          call tensor_random(vovo)
          call tensor_scale(ccsd_doubles,1.0E-2_realk)
          call tensor_scale(vovo,1.0E-2_realk)

       elseif (nodtotal .eq. 1) then

          if (abc) then

             call tensor_init(vovo,[nocc,nocc,nvirt,nvirt],4,bg=use_bg)
             call tensor_init(ccsd_doubles,[nocc,nocc,nvirt,nvirt],4,bg=use_bg)

          else

             call tensor_init(vovo,[nvirt,nvirt,nocc,nocc],4,bg=use_bg)
             call tensor_init(ccsd_doubles,[nvirt,nvirt,nocc,nocc],4,bg=use_bg)

          endif

          call tensor_random(ccsd_doubles)
          call tensor_random(vovo)
          call dscal(nocc**2*nvirt**2,1.0E-2_realk,ccsd_doubles%elm1,1)
          call dscal(nocc**2*nvirt**2,1.0E-2_realk,vovo%elm1,1)

       endif
#else
       if (abc) then

          call tensor_init(vovo,[nocc,nocc,nvirt,nvirt],4,bg=use_bg)
          call tensor_init(ccsd_doubles,[nocc,nocc,nvirt,nvirt],4,bg=use_bg)

       else

          call tensor_init(vovo,[nvirt,nvirt,nocc,nocc],4)
          call tensor_init(ccsd_doubles,[nvirt,nvirt,nocc,nocc],4)

       endif

       call tensor_random(vovo)
       call tensor_random(ccsd_doubles)
       call dscal(nocc**2*nvirt**2,1.0E-2_realk,ccsd_doubles%elm1,1)
       call dscal(nocc**2*nvirt**2,1.0E-2_realk,vovo%elm1,1)

#endif

    else

#ifdef VAR_MPI
       if (nodtotal .gt. 1) then

          if (abc) then

             call tensor_minit(vovo,[nocc,nocc,nvirt,nvirt],4,tdims=[nocc,nocc,abc_tile_size,abc_tile_size],atype='TDAR',bg=use_bg)
             transform_bg = use_bg .and.  (mem_get_bg_buf_free()>=(i8*nocc*nocc)*nvirt*nvirt)
             call tensor_init(tmp_tensor_1,[nocc,nocc,nvirt,nvirt],4,bg=transform_bg)
             call tensor_cp_data(vovo_in,tmp_tensor_1,order=[2,4,1,3])
             call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,oovv=tmp_tensor_1%elm1)
             call tensor_cp_data(tmp_tensor_1,vovo)
             call tensor_free(tmp_tensor_1)
             call tensor_minit(ccsd_doubles,[nocc,nocc,nvirt,nvirt],4,tdims=[nocc,nocc,nvirt,abc_tile_size],atype='TDAR',bg=use_bg)
             transform_bg = use_bg .and.  (mem_get_bg_buf_free()>=(i8*nocc*nocc)*nvirt*nvirt)
             call tensor_init(tmp_tensor_2,[nocc,nocc,nvirt,nvirt],4,bg=transform_bg)
             call tensor_cp_data(ccsd_doubles_in,tmp_tensor_2,order=[2,4,3,1])
             call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,oovv=tmp_tensor_2%elm1)
             call tensor_cp_data(tmp_tensor_2,ccsd_doubles)
             call tensor_free(tmp_tensor_2)

          else

             call tensor_minit(vovo,[nvirt,nvirt,nocc,nocc],4,&
                &tdims=[nvirt,nvirt,ijk_tile_size,ijk_tile_size],atype='TDAR',bg=use_bg)
             transform_bg = use_bg .and.  (mem_get_bg_buf_free()>=(i8*nocc*nocc)*nvirt*nvirt)
             call tensor_init(tmp_tensor_1,[nvirt,nvirt,nocc,nocc],4,bg=transform_bg)
             call tensor_cp_data(vovo_in,tmp_tensor_1,order=[1,3,2,4])
             call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,vvoo=tmp_tensor_1%elm1)
             call tensor_cp_data(tmp_tensor_1,vovo)
             call tensor_free(tmp_tensor_1)
             call tensor_minit(ccsd_doubles,[nvirt,nvirt,nocc,nocc],4,&
                &tdims=[nvirt,nvirt,nocc,ijk_tile_size],atype='TDAR',bg=use_bg)
             transform_bg = use_bg .and.  (mem_get_bg_buf_free()>=(i8*nocc*nocc)*nvirt*nvirt)
             call tensor_init(tmp_tensor_2,[nvirt,nvirt,nocc,nocc],4,bg=transform_bg)
             call tensor_cp_data(ccsd_doubles_in,tmp_tensor_2,order=[1,3,4,2])
             call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,vvoo=tmp_tensor_2%elm1)
             call tensor_cp_data(tmp_tensor_2,ccsd_doubles)
             call tensor_free(tmp_tensor_2)

          endif

       elseif (nodtotal .eq. 1) then

          if (abc) then

             call tensor_init(vovo,[nocc,nocc,nvirt,nvirt],4,bg=use_bg)
             call array_reorder_4d(1.0E0_realk,vovo_in%elm1,nvirt,nocc,nvirt,nocc,[2,4,1,3],0.0E0_realk,vovo%elm1)
             call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,oovv=vovo%elm1)
             call tensor_init(ccsd_doubles,[nocc,nocc,nvirt,nvirt],4,bg=use_bg)
             call array_reorder_4d(1.0E0_realk,ccsd_doubles_in%elm1,nvirt,nocc,nvirt,nocc,[2,4,3,1],0.0E0_realk,ccsd_doubles%elm1)
             call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,oovv=ccsd_doubles%elm1)

          else

             call tensor_init(vovo,[nvirt,nvirt,nocc,nocc],4,bg=use_bg)
             call array_reorder_4d(1.0E0_realk,vovo_in%elm1,nvirt,nocc,nvirt,nocc,[1,3,2,4],0.0E0_realk,vovo%elm1)
             call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,vvoo=vovo%elm1)
             call tensor_init(ccsd_doubles,[nvirt,nvirt,nocc,nocc],4,bg=use_bg)
             call array_reorder_4d(1.0E0_realk,ccsd_doubles_in%elm1,nvirt,nocc,nvirt,nocc,[1,3,4,2],0.0E0_realk,ccsd_doubles%elm1)
             call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,vvoo=ccsd_doubles%elm1)

          endif

       endif
#else
       if (abc) then

          call tensor_init(vovo,[nocc,nocc,nvirt,nvirt],4,bg=use_bg)
          call array_reorder_4d(1.0E0_realk,vovo_in%elm1,nvirt,nocc,nvirt,nocc,[2,4,1,3],0.0E0_realk,vovo%elm1)
          call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,oovv=vovo%elm1)
          call tensor_init(ccsd_doubles,[nocc,nocc,nvirt,nvirt],4,bg=use_bg)
          call array_reorder_4d(1.0E0_realk,ccsd_doubles_in%elm1,nvirt,nocc,nvirt,nocc,[2,4,3,1],0.0E0_realk,ccsd_doubles%elm1)
          call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,oovv=ccsd_doubles%elm1)

       else

          call tensor_init(vovo,[nvirt,nvirt,nocc,nocc],4,bg=use_bg)
          call array_reorder_4d(1.0E0_realk,vovo_in%elm1,nvirt,nocc,nvirt,nocc,[1,3,2,4],0.0E0_realk,vovo%elm1)
          call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,vvoo=vovo%elm1)
          call tensor_init(ccsd_doubles,[nvirt,nvirt,nocc,nocc],4,bg=use_bg)
          call array_reorder_4d(1.0E0_realk,ccsd_doubles_in%elm1,nvirt,nocc,nvirt,nocc,[1,3,4,2],0.0E0_realk,ccsd_doubles%elm1)
          call local_can_trans(nocc,nvirt,nbasis,Uocc%elm2,Uvirt%elm2,vvoo=ccsd_doubles%elm1)

       endif
#endif

    endif

  end subroutine convert_ccsd_and_vovo

!endif mod_unreleased
#endif

  subroutine dummy_ccsdpt_routine()

  end subroutine dummy_ccsdpt_routine

end module ccsdpt_module

#ifdef MOD_UNRELEASED

#ifdef VAR_MPI

  subroutine ccsdpt_slave_info()

  use infpar_module
  use lsmpi_type
  use decmpi_module
  use ccsdpt_module, only: ccsdpt_info

    implicit none
    integer :: nocc, nvirt,nbasis
    logical :: print_frags,abc
    integer :: ijk_nbuffs,abc_nbuffs,ijk_tile_size,abc_tile_size,nodtotal

    ! none of the variable entering here are used - thus, none are initialized
    call ccsdpt_info(nbasis,nocc,nvirt,print_frags,abc,ijk_nbuffs,abc_nbuffs,ijk_tile_size,abc_tile_size,nodtotal)

  end subroutine ccsdpt_slave_info

  !> \brief slaves enter here from lsmpi_slave (or dec_lsmpi_slave) and need to get to work 
  !> \author Janus Juul Eriksen
  !> \date x-mas 2012
  subroutine ccsdpt_slave_work()

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
    type(tensor) :: vovo,ccsd_t2,ccsd_t1,ccsdpt_t1,ccsdpt_t2
    real(realk) :: ccsdpt_e4,ccsdpt_e5
    type(lsitem) :: mylsitem
    logical :: print_frags,abc
    integer :: ijk_nbuffs,abc_nbuffs,abc_tile_size,nodtotal

    abc = .false.
    print_frags = .false.

    call time_start_phase(PHASE_COMM)

    ! call ccsd(t) data routine in order to receive data from master
    if (print_frags) then

       call mpi_communicate_ccsdpt_calcdata(nocc,nvirt,nbasis,vovo,ccsd_t2,mylsitem,print_frags,abc)

    else

       call mpi_communicate_ccsdpt_calcdata(nocc,nvirt,nbasis,vovo,ccsd_t2,mylsitem,print_frags,abc,ccsd_t1)

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

       ! init ccsd singles
       if (abc) then

          call tensor_init(ccsd_t1, [nocc,nvirt],2)

       else

          call tensor_init(ccsd_t1, [nvirt,nocc],2)

       endif

       call ls_mpibcast(ccsd_t1%elm1,nvirt*nocc,infpar%master,infpar%lg_comm)

       ccsdpt_e4 = 0.0E0_realk
       ccsdpt_e5 = 0.0E0_realk

    endif


    call time_start_phase(PHASE_WORK)

    ! now enter the ccsd(t) driver routine
    if (print_frags) then

       call ccsdpt_driver(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,mylsitem,vovo,ccsd_t2,&
                               & print_frags,abc,ccsdpt_singles=ccsdpt_t1,ccsdpt_doubles=ccsdpt_t2)

    else

       call ccsdpt_driver(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,mylsitem,vovo,ccsd_t2,&
                               & print_frags,abc,e4=ccsdpt_e4,e5=ccsdpt_e5,ccsd_singles=ccsd_t1)

    endif


    call time_start_phase(PHASE_WORK)

    ! now, release all amplitude arrays, both ccsd and ccsd(t)
    if (print_frags) then

       call tensor_free(ccsdpt_t1)
       call tensor_free(ccsdpt_t2)

    else

       call tensor_free(ccsd_t1)

    endif

    call ls_free(mylsitem)

  end subroutine ccsdpt_slave_work

#endif
!endif mod_unreleased
#endif
