#ifdef VAR_MPI

subroutine dec_lsmpi_slave(comm)
   use precision
   use lstiming
   use infpar_module
   use lsmpi_type
   use integralinterfaceMod!, only: II_screeninit, &
   !           & II_bcast_screen, II_screenfree
   use tensor_interface_module, only: lspdm_init_global_buffer,lspdm_free_global_buffer
   use dec_driver_slave_module, only: main_fragment_driver_slave
   use f12_integrals_module 

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
         call get_slaves_to_tensor_test
         ! DEC MP2 integrals and amplitudes
      case(MP2INAMP);
         call MP2_integrals_and_amplitudes_workhorse_slave
         ! DEC MP2 RI energy
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
#ifdef MOD_UNRELEASED 
      case(F12_INTEGRAL_CALCULATION);
         call get_f12_fragment_energy_slave
#endif
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
         call PDM_TENSOR_SLAVE(comm)
      case(INITSLAVETIME);
         call init_slave_timers_slave(comm)
      case(GETSLAVETIME);
         call get_slave_timers_slave(comm)

      case(JOB_LSPDM_INIT_GLOBAL_BUFFER);
         call lspdm_init_global_buffer(.false.)
      case(JOB_LSPDM_FREE_GLOBAL_BUFFER);
         call lspdm_free_global_buffer(.false.)
      case(INIT_BG_BUF);
         call mem_init_background_alloc_slave(comm)
      case(FREE_BG_BUF);
         call mem_free_background_alloc_slave(comm)

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

