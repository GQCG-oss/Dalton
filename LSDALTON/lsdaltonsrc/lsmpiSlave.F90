    subroutine lsmpi_init
#ifdef VAR_MPI
    use lsmpi_type
    use infpar_module
    use ls_env
    implicit none
    integer(kind=ls_mpik) :: ierr
#ifdef VAR_CHEMSHELL
    integer(kind=ls_mpik) :: lsdalton_chemsh_comm
    external lsdalton_chemsh_comm
#endif
    
    nLog = 0
    nDP = 0
    nInteger4=0
    nInteger8=0
    nShort = 0
    nCha = 0
#ifdef VAR_CHEMSHELL
    MPI_COMM_LSDALTON = lsdalton_chemsh_comm()
#else
    call MPI_INIT( ierr )
    MPI_COMM_LSDALTON = MPI_COMM_WORLD
    !asynchronous progress is off åer default, might be switched on with an
    !environment variable
    LSMPIASYNCP = .false.
    call ls_getenv(varname="LSMPI_ASYNC_PROGRESS",leng=20,output_bool=LSMPIASYNCP)
#endif

    call MPI_COMM_RANK( MPI_COMM_LSDALTON, infpar%mynum, ierr )
    call MPI_COMM_SIZE( MPI_COMM_LSDALTON, infpar%nodtot, ierr )
    infpar%master = 0;
    
    ! Set rank, sizes and communcators for local groups
    ! to be identical to those for world group by default.
    call lsmpi_default_mpi_group
    ! Assume that there will be local jobs
    infpar%lg_morejobs=.true.

    if (infpar%mynum.ne.infpar%master) then

      call lsmpi_slave(MPI_COMM_LSDALTON)

    endif
#endif
    end subroutine lsmpi_init 
#ifdef VAR_MPI

    subroutine lsmpi_slave(comm)
      use infpar_module
      use lsmpi_type
      use lsmpi_test
      use integralinterfaceMod
      use dec_driver_slave_module
      use lspdm_tensor_operations_module,only:free_persistent_array
    implicit none
    !> Communicator from which task is to be received
    integer(kind=ls_mpik) :: comm 
    integer(kind=ls_mpik) :: ierr
    integer :: job
100   call ls_mpibcast(job,infpar%master,comm)

      select case(job)
      case(MATRIXTY);
         call lsmpi_set_matrix_type_slave
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
      case(CCSDDATA);
         call ccsd_data_preparation
      case(CCSDSLV4E2);
         call calculate_E2_and_permute_slave
#ifdef MOD_UNRELEASED
      case(CCSDPTSLAVE);
         call ccsdpt_slave
!endif mod_unreleased
#endif
      case(ARRAYTEST);
         call get_slaves_to_array_test
      case(GROUPINIT);
         call init_mpi_groups_slave
! DEC driver - main loop
      case(DECDRIVER);
         call main_fragment_driver_slave
      case(DEFAULTGROUPS);
         call lsmpi_default_mpi_group
#ifdef VAR_SCALAPACK
      case(GRIDINIT);
         call PDM_GRIDINIT_SLAVE
      case(GRIDEXIT);
         call PDM_GRIDEXIT_SLAVE
      case(PDMSLAVE);
         call PDM_SLAVE
#endif
      case(PDMA4SLV);
         call PDM_ARRAY_SLAVE
         ! Quit but there are more local jobs - used for local group handling
      case(QUITMOREJOBS); 
         infpar%lg_morejobs=.true.
         goto 200
         ! Quit there are NO more local jobs - used for local group handling
      case(QUITNOMOREJOBS); 
         infpar%lg_morejobs=.false.
         goto 200
      case(LSMPIQUIT);  
         call free_persistent_array()
         call lsmpi_finalize(6,.FALSE.)
      case(LSMPIQUITINFO);
         call free_persistent_array()
         call lsmpi_finalize(6,.TRUE.)
      case(LSMPITEST);
         call test_mpi(comm)
      case default
         call free_persistent_array()
         call lsmpi_finalize(6,.FALSE.)
      end select
     
      goto 100

      200 return

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

