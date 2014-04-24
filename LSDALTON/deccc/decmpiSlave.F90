#ifdef VAR_MPI

subroutine dec_lsmpi_slave(comm)
   use precision
   use lstiming
   use infpar_module
   use lsmpi_type
   use integralinterfaceMod!, only: II_screeninit, &
   !           & II_bcast_screen, II_screenfree
   use dec_driver_slave_module, only: main_fragment_driver_slave

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
      case(LSGETINT);
         call lsmpi_getIntegrals_Slave(comm)
      case(LSJENGIN);
         call lsmpi_jengine_Slave(comm)
      case(LSLINK);
         call lsmpi_linK_Slave(comm)
      case(DFTSETFU);
         call lsmpi_setSlaveFunc
      case(LSMPI_IIDFTKSM);
         call lsmpi_II_DFT_KSM_Slave(comm)
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
      case(ARRAYTEST);
         call get_slaves_to_array_test
         ! DEC MP2 integrals and amplitudes
      case(MP2INAMP);
         call MP2_integrals_and_amplitudes_workhorse_slave
         ! DEC MP2 RI energy
      case(MP2INAMPRI);
         call MP2_RI_EnergyContribution_slave
      case(CCSDDATA);
         call ccsd_data_preparation
      case(CCGETGMO);
         call cc_gmo_data_slave
      case(MOCCSDDATA);
         call moccsd_data_slave
      case(MO_INTEGRAL_SIMPLE);
         call get_mo_integral_par_slave
      case(CCSDSLV4E2);
         call calculate_E2_and_permute_slave
      case(RPAGETRESIDUAL);
         call rpa_res_slave
      case(RPAGETFOCK);
         call rpa_fock_slave
#ifdef MOD_UNRELEASED 
      case(CCSDPTSLAVE);
         call ccsdpt_slave
#endif
      case(SIMPLE_MP2_PAR);
         call get_simple_parallel_mp2_residual_slave
      case(GROUPINIT);
         call init_mpi_groups_slave
         ! DEC driver - main loop
      case(DECDRIVER);
         call main_fragment_driver_slave
      case(DEFAULTGROUPS);
         call lsmpi_default_mpi_group
      case(PDMA4SLV);
         call PDM_ARRAY_SLAVE(comm)
      case(INITSLAVETIME);
         call init_slave_timers_slave(comm)
      case(GETSLAVETIME);
         call get_slave_timers_slave(comm)
#ifdef VAR_SCALAPACK
      case(GRIDINIT);
         call PDM_GRIDINIT_SLAVE
      case(GRIDEXIT);
         call PDM_GRIDEXIT_SLAVE
      case(PDMSLAVE);
         call PDM_SLAVE
#endif
         ! Quit but there are more local jobs - used for local group handling
      case(QUITMOREJOBS); 
         infpar%lg_morejobs   = .true.
         stay_in_slaveroutine = .false.
         ! Quit there are NO more local jobs - used for local group handling
      case(QUITNOMOREJOBS); 
         infpar%lg_morejobs   = .false.
         stay_in_slaveroutine = .false.
      end select
   end do


end subroutine dec_lsmpi_slave

#else

!Added to avoid "has no symbols" linking warning
subroutine dec_mpi_void()
end subroutine dec_mpi_void

#endif

