subroutine lsmpi_init(OnMaster)
#ifdef VAR_MPI
   use lsmpi_type
   use infpar_module
   use ls_env
   use matrix_operations_pdmm
#ifdef VAR_SCALAPACK  
   use matrix_operations_scalapack
#endif
   implicit none
   logical, intent(inout) :: OnMaster
   integer(kind=ls_mpik)  :: ierr
#ifdef VAR_CHEMSHELL
   integer(kind=ls_mpik)  :: lsdalton_chemsh_comm
   external lsdalton_chemsh_comm
#endif

   nLog      = 0
   nDP       = 0
   nInteger4 = 0
   nInteger8 = 0
   nCha      = 0

   if (call_mpi_init) then
#ifdef VAR_CHEMSHELL
     MPI_COMM_LSDALTON = lsdalton_chemsh_comm()
#else
     call MPI_INIT( ierr )
     MPI_COMM_LSDALTON        = MPI_COMM_WORLD
#endif
   endif

   !Set default mpi message sizes
   SPLIT_MPI_MSG      = 100000000
   MAX_SIZE_ONE_SIDED =  12500000

   lsmpi_enabled_comm_procs = .false.

   !TEST WHETHER THE MPI VARIABLES AND THE INTERNALLY DECLARED ls_mpik have the
   !same kind
   if(huge(ierr) /= huge(MPI_COMM_WORLD))then
      print *,"WARNING: THE MPI INTRINSIC VARIABLES AND THE ls_mpik ARE OF&
      & DIFFERENT KIND, YOUR JOB IS LIKELY TO FAIL",huge(ierr),huge(MPI_COMM_WORLD)
   endif

   !asynchronous progress is off per default, might be switched on with an
   !environment variable, or by spawning communication processes
   LSMPIASYNCP             = .false.
   call ls_getenv(varname="LSMPI_ASYNC_PROGRESS",leng=20,output_bool=LSMPIASYNCP)


   call MPI_COMM_GET_PARENT( infpar%parent_comm, ierr )
   call get_rank_for_comm( MPI_COMM_LSDALTON, infpar%mynum  )
   call get_size_for_comm( MPI_COMM_LSDALTON, infpar%nodtot )
   !default number of nodes to be used with scalapack - can be changed
   infpar%ScalapackGroupSize = infpar%nodtot
   infpar%ScalapackWorkaround = .FALSE.
   infpar%PDMMGroupSize = infpar%nodtot
#ifdef VAR_SCALAPACK  
   scalapack_mpi_set = .FALSE.
#endif   
   pdmm_mpi_set = .FALSE.
   infpar%master = int(0,kind=ls_mpik);

   !CHECK IF WE ARE ON THE MAIN MAIN MAIN MASTER PROCESS (NOT IN LOCAL GROUP,
   !NOT A SPAWNED PROCESS)
   OnMaster = (infpar%mynum==infpar%master .and. infpar%parent_comm == MPI_COMM_NULL) 

   ! Set rank, sizes and communcators for local groups
   ! to be identical to those for world group by default.
   call lsmpi_default_mpi_group

   !if i am a child process, set the intracommunicator to the parent proc
   if( infpar%parent_comm /= MPI_COMM_NULL )  call get_parent_child_relation

   ! Assume that there will be local jobs
   infpar%lg_morejobs=.true.



   if (infpar%mynum.ne.infpar%master) then

      !if normal slave, listen on MPI_COMM_LSDALTON
      call lsmpi_slave(MPI_COMM_LSDALTON)

      elseif( infpar%parent_comm /= MPI_COMM_NULL )then

      !if spawned process, listen on the parent-child-intracomm
      call lsmpi_slave(infpar%pc_comm)

   endif
#else
   logical, intent(inout) :: OnMaster
   !IF NOT COMPILED WITH MPI SET MASTER = TRUE
   OnMaster = .true.
#endif
end subroutine lsmpi_init 
#ifdef VAR_MPI

subroutine lsmpi_slave(comm)
   use lsparameters
   use lstiming
   use infpar_module
   use lsmpi_type
   use lsmpi_test
   use integralinterfaceMod
   use dec_driver_slave_module
   use tensor_interface_module
#ifdef VAR_SCALAPACK  
   use matrix_operations_scalapack
#endif
   implicit none
   !> Communicator from which task is to be received
   integer(kind=ls_mpik) :: comm 
   integer(kind=ls_mpik) :: ierr
   integer :: job
   logical :: stay_in_slaveroutine

   stay_in_slaveroutine = .true.

   do while(stay_in_slaveroutine)

      call time_start_phase(PHASE_IDLE)
      call ls_mpibcast(job,infpar%master,comm)
      call time_start_phase(PHASE_WORK)


      select case(job)
      case(MATRIXTY);
         call lsmpi_set_matrix_type_slave
      case(MATRIXTY2);
         call lsmpi_set_matrix_tmp_type_slave
      case(LSGETINT);
         call lsmpi_getIntegrals_Slave(comm)
      case(LSJENGIN);
         call lsmpi_jengine_Slave(comm)
      case(LSLINK);
         call lsmpi_linK_Slave(comm)
      case(DFTSETFU);
         call lsmpi_setSlaveFunc
      case(DFTADDFU);
         call lsmpi_addSlaveFunc
      case(LSMPI_IIDFTKSM);
         call lsmpi_II_DFT_KSM_Slave(comm)
      case(LSMPI_IIDFTABSVALOVERLAP);
         call lsmpi_II_DFT_ABSVALOVERLAP_Slave(comm)
      case(LSMPI_IIDFTKSME);
         call lsmpi_II_DFT_KSME_Slave(comm)
      case(IIDFTGEO);
         call lsmpi_geoderiv_molgrad_Slave(comm)
      case(IIDFTLIN);
         call lsmpi_linrsp_Slave(comm)
      case(IIDFTQRS);
         call lsmpi_quadrsp_Slave(comm)
      case(IIDFTMAG);
         call lsmpi_magderiv_Slave(comm)
      case(IIDFTMAL);
         call lsmpi_magderiv_linrsp_Slave(comm)
      case(IIDFTGKS);
         call lsmpi_geoderiv_F_Slave(comm)
      case(IIDFTGLR);
         call lsmpi_geoderiv_G_Slave(comm)
      case(IISCREENINIT);
         call II_screeninit(comm)
      case(IISCREEN);
         call II_bcast_screen(comm)
      case(IISCREENFREE);
         call II_screenfree(comm)
         ! DEC MP2 integrals and amplitudes
      case(MP2INAMP);
         call MP2_integrals_and_amplitudes_workhorse_slave
      case(RIMP2INAMP);
         call RIMP2_integrals_and_amplitudes_slave
      case(LSTHCRIMP2INAMP);
         call LSTHCRIMP2_integrals_and_amplitudes_slave
      case(RIMP2FULL);
         call full_canonical_rimp2_slave
      case(LSTHCRIMP2FULL);
!         call full_canonical_ls_thc_rimp2_slave
      case(CANONMP2FULL);
         call full_canonical_mp2_slave
      case(DEC_SETTING_TO_SLAVES);
         call set_dec_settings_on_slaves
      case(CCSDDATA);
         call ccsd_data_preparation
      case(MO_INTEGRAL_SIMPLE);
         call get_mo_integral_par_slave
      case(CCSDSLV4E2);
         call calculate_E2_and_permute_slave
!      case(RPAGETRESIDUAL);
!         call rpa_res_slave
      case(RPAGETFOCK);
         call rpa_fock_slave
#ifdef MOD_UNRELEASED
      case(CCGETGMO);
         call cc_gmo_data_slave
      case(MOCCSDDATA);
         call moccsd_data_slave
      case(CCSDPTSLAVE_INFO);
         call ccsdpt_slave_info
      case(CCSDPTSLAVE_WORK);
         call ccsdpt_slave_work
#endif
      case(SIMPLE_MP2_PAR);
         call get_simple_parallel_mp2_residual_slave
      case(ARRAYTEST);
         call get_slaves_to_tensor_test
      case(GROUPINIT);
         call init_mpi_groups_slave
         ! DEC driver - main loop
      case(DECDRIVER);
         call main_fragment_driver_slave
      case(DEFAULTGROUPS);
         call lsmpi_default_mpi_group
#ifdef VAR_SCALAPACK
      case(GRIDINIT);
         IF(scalapack_member)THEN
            call PDM_GRIDINIT_SLAVE
         ENDIF
      case(GRIDEXIT);
         IF(scalapack_member)THEN
            call PDM_GRIDEXIT_SLAVE
         ENDIF
      case(PDMSLAVE);
         IF(scalapack_member)THEN
            !inside the PDM_SLAVE the 
            !scalapack_comm communicator is used
            !only call PDM_SLAVE if this rank is 
            !a member of this communicator
            call PDM_SLAVE()
         ENDIF
#endif
      case(PDMMGRIDINIT);
         call PDMM_GRIDINIT_SLAVE
      case(PDMMGRIDEXIT);
         call PDMM_GRIDEXIT_SLAVE
      case(PDMA4SLV);
         call PDM_tensor_SLAVE(comm)
      case(INITSLAVETIME);
         call init_slave_timers_slave(comm)
      case(GETSLAVETIME);
         call get_slave_timers_slave(comm)
      case(GIVE_BIRTH);
         call give_birth_to_child_process

      case(LSPDM_GIVE_BIRTH);
         call lspdm_start_up_comm_procs
      case(LSMPITEST);
         call test_mpi(comm)


      case(LSMPIPRINTINFO);
         call lsmpi_print_mem_info(6,.TRUE.)

      case(SLAVES_SHUT_DOWN_CHILD)
         if(infpar%parent_comm/=MPI_COMM_NULL)then
            call lsquit("ERROR(SLAVES_SHUT_DOWN_CHILD):I am not a parent",-1)
         endif
         call shut_down_child_process

      case(LSPDM_SLAVES_SHUT_DOWN_CHILD)
         if(infpar%parent_comm/=MPI_COMM_NULL)stay_in_slaveroutine = .false.
         call lspdm_shut_down_comm_procs

      case(JOB_LSPDM_INIT_GLOBAL_BUFFER);
         call lspdm_init_global_buffer(.false.)
      case(JOB_LSPDM_FREE_GLOBAL_BUFFER);
         call lspdm_free_global_buffer(.false.)

      case(SET_GPUMAXMEM);
         call ls_mpibcast(GPUMAXMEM,infpar%master,comm)

      case(SET_SPLIT_MPI_MSG);
         call ls_mpibcast(SPLIT_MPI_MSG,infpar%master,comm)
      case(SET_MAX_SIZE_ONE_SIDED);
         call ls_mpibcast(MAX_SIZE_ONE_SIDED,infpar%master,comm)

      case(SET_TENSOR_BACKEND_TRUE);
         call tensor_set_dil_backend_true(.false.)

      case(SET_TENSOR_DEBUG_TRUE);
         call tensor_set_debug_mode_true(.false.)

      case(SET_TENSOR_ALWAYS_SYNC_TRUE);
         call tensor_set_always_sync_true(.false.)

      case(SET_TENSOR_SEG_LENGTH);
         call tensor_set_global_segment_length(0_long)

      case(INIT_BG_BUF);
         call mem_init_background_alloc_slave(comm)
      case(FREE_BG_BUF);
         call mem_free_background_alloc_slave(comm)

         !##########################################
         !########  QUIT THE SLAVEROUTINE ##########
         !##########################################

      case(CHILD_SHUT_DOWN);
         if(infpar%parent_comm==MPI_COMM_NULL)then
            call lsquit("ERROR(SHUT_DOWN_CHILD):I am not a child to be shut down",-1)
         endif
         call shut_down_child_process
         stay_in_slaveroutine = .false.

         ! Quit but there are more local jobs - used for local group handling
      case(QUITMOREJOBS); 
         infpar%lg_morejobs   = .true.
         stay_in_slaveroutine = .false.
         ! Quit there are NO more local jobs - used for local group handling
      case(QUITNOMOREJOBS); 
         infpar%lg_morejobs   = .false.
         stay_in_slaveroutine = .false.
      case(LSMPIQUIT);  
         stay_in_slaveroutine = .false.
      case default
         call free_persistent_array()
         call lsmpi_finalize(6,.FALSE.)
         stay_in_slaveroutine = .false.
      end select

   end do

end subroutine lsmpi_slave

!> Waiting position for DEC slaves to receive job to be carried out within local group,
!> i.e., the communicator is always infpar%lg_comm.
!> \author Kasper Kristensen 
!> \date April 2012
!!$    subroutine lsmpi_local_slave(job)
!!$      use infpar_module
!!$      use lsmpi_type
!!$      implicit none
!!$      !> Give last job as output, because it is convenient
!!$      !> to distinguish between two different types of quit statements:
!!$      !> job=0, quit local slave - but you may be sent back to local slave later
!!$      !> job=-1, quit local slave - there are no more jobs to be done
!!$      integer :: job
!!$print *, 'entering local slave first'
!!$100   call ls_mpibcast(job,0,infpar%lg_comm)
!!$print *, 'entering local slave - job', job
!!$      select case(job)
!!$      case(MP2INAMP);        ! MP2
!!$         call MP2_integrals_and_amplitudes_workhorse_slave
!!$      case(CCSDDATA);       ! CC2 or CCSD
!!$         call ccsd_data_preparation
!!$      case(0) ! quit local slave, but more jobs (see above)
!!$         goto 200
!!$      case(-1) ! quit local slave, no more jobs (see above)
!!$         goto 200
!!$      end select
!!$
!!$      goto 100
!!$
!!$200   return
!!$
!!$    end subroutine lsmpi_local_slave


#endif

