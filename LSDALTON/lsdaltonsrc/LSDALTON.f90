!> @file 
!> Contains main SCF driver, some module wrappers and miscellaneous 

!> \brief Driver for stand-alone f90 linear scaling SCF.
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2008-10-26
SUBROUTINE lsdalton
  use precision
  use configurationType, only: configitem, LowAccuracyStartType, &
       & set_Low_accuracy_start_settings,revert_Low_accuracy_start_settings
  use TYPEDEFTYPE, only: lsitem
  use matrix_module
  use memory_handling, only: mem_alloc,mem_dealloc, stats_mem
  use matrix_operations, only: set_matrix_default,max_no_of_matrices, no_of_matrices, &
       & no_of_matmuls, mat_init, mat_free, mat_assign, &
       & mat_mul, mat_no_of_matmuls, mat_write_to_disk, mat_read_from_disk, mat_diag_f,&
       & mat_TrAB, mat_print, MatrixmemBuf_init, MatrixmemBuf_free, &
       & MatrixmemBuf_print
  use configuration, only: config_shutdown, config_free,scf_purify
  use lsdalton_fock_module, only: lsint_fock_data
  use init_lsdalton_mod, only: open_lsdalton_files,init_lsdalton_and_get_lsitem
  use initial_guess, only: get_initial_dens
  use scfloop_module, only: scfloop, scf_afterplay
  use lstiming, only: lstimer, init_timers, print_timers
  use ks_settings, only: ks_init_incremental_fock, ks_free_incremental_fock
  use files, only: lsopen, lsclose
  use decompMod, only: decomp_init, decomp_shutdown, decomposition, get_oao_transformed_matrices
  use matrix_util, only: save_fock_matrix_to_file, save_overlap_matrix_to_file, util_mo_to_ao_2
  use daltoninfo, only: ls_free 
  !For making orbital .plt files
  use print_moorb_grid_mod
  ! Debug and Testing
  use dal_interface, only: di_debug_general, di_debug_general2
  use extra_output, only: print_orbital_info2
  ! Profile 
  use profile_int, only: di_profile_lsint
  ! DEC 
  use DEC_typedef_module, only: DECinfo  
  ! PROPERTIES SECTION
  use lsdalton_rsp_mod, only: lsdalton_response, get_excitation_energy
  ! DYNAMICS
  use dynamics_driver, only: LS_dyn_run
  ! SOEO
  use soeo_loop, only: soeoloop, soeo_restart
  ! GEO OPTIMIZER
  use ls_optimizer_mod, only: LS_RUNOPT
  use lsmpi_type, only: lsmpi_finalize
  use lstensorMem, only: lstmem_init, lstmem_free
  use numerical_hessian, only: get_numerical_hessian
  use pbc_setup, only: set_pbc_molecules
  use molecular_hessian_mod, only: get_molecular_hessian
  use test_molecular_hessian_mod, only: test_Hessian_contributions
  use rsp_util, only: init_rsp_util
#ifdef VAR_PAPI
  use papi_module
#endif
  use integralinterfaceMod, only: II_get_overlap, II_get_h1, &
       & II_precalc_ScreenMat, II_get_GaussianGeminalFourCenter
  use dec_main_mod!, only: dec_main_prog
  implicit none
  integer             :: nbast,lupri, luerr, lucmo
  TYPE(lsitem),target :: ls
  type(configItem),target  :: config
  real(realk)         :: t1,t2,TIMSTR,TIMEND
  TYPE(Matrix)         :: F(1),D(1), CMO
  Type(Matrix), target :: H1,S
  integer             :: matmultot, lun
  REAL(REALK)         :: mx
  ! Energy
  REAL(REALK)         :: E(1),ExcitE
  logical             :: mem_monitor,do_decomp
  real(realk), allocatable :: eival(:)
  real(realk),pointer :: GGem(:,:,:,:,:)
  integer     :: lusoeo,funit
  logical     :: soeosaveexist, HFdone,OnMaster,scfpurify
  type(matrix) :: Dmo, tmp
  integer             :: nelec
  Integer             :: Natoms
  Real(realk),pointer :: geomHessian(:,:)
  type(LowAccuracyStartType)  :: LAStype
  Interface 
     subroutine optimloc(CMO,nocc,m,ls,CFG)
       use davidson_settings,only: RedSpaceItem
       use matrix_module !matrix
       use typedeftype !lsitem
       implicit none
       type(RedSpaceItem) :: CFG
       type(Matrix), target:: CMO
       TYPE(lsitem) , intent(inout) :: ls
       integer,       intent(in)    :: nocc
       integer,       intent(in)    :: m(2)
     end subroutine optimloc
  end Interface 
  OnMaster = .TRUE.
  ! Set lupri and luerr to zero here to make sure they are not unintialized for MPI slaves
  luerr=0
  lupri=0


  call lsinit_all()

#if VAR_DEBUG
  print *,        "THIS IS A DEBUG BUILD"
  write (LUPRI,*),"THIS IS A DEBUG BUILD"
  write (LUERR,*),"THIS IS A DEBUG BUILD"
#endif

  ! Open output files LSDALTON.OUT and LSDALTON.ERR
  call open_lsdalton_files(lupri,luerr)

  ! Time the whole LSdalton calculation
  call LSTIMER('START',t1,t2,LUPRI)

  ! Init LSdalton calculation and get lsitem and config structures
  call init_lsdalton_and_get_lsitem(lupri,luerr,nbast,ls,config,mem_monitor)
  ! Timing of individual steps
  CALL LSTIMER('START ',TIMSTR,TIMEND,lupri)

  IF (config%integral%debugUncontAObatch) THEN 
     call II_test_uncontAObatch(lupri,luerr,ls%setting) 
  ENDIF

  if(config%prof%doProf)then
     call di_profile_lsint(ls,config,lupri,nbast)
     call lsmpi_finalize(lupri,config%mpi_mem_monitor)
     return
  endif
  ! Vladimir Rybkin: If we do dec and optimize geometry or run dynamics 
  ! than we don't skip HF calculation
  if ( ( (config%optinfo%optimize .EQV. .TRUE.) .OR. (config%dynamics%do_dynamics .EQV. .TRUE.)) &
       & .AND. (DECinfo%doDEC)) then
     DECinfo%doHF = .TRUE.
  endif

  ! Read in already optimized HF orbitals, and localize
  OnlyLoc:  if (config%davidOrbLoc%OnlyLocalize) then
           ! read orbitals
           lun = -1
	   call mat_init(CMO,nbast,nbast)
           CALL LSOPEN(lun,'orbitals_in.u','unknown','UNFORMATTED')
	   call mat_read_from_disk(LUN,cmo,OnMaster)
           call LSclose(LUN,'KEEP')
           if (config%decomp%cfg_mlo) then
	      write(ls%lupri,*)
	      write(ls%lupri,'(a)') '*** LOCALIZING ORBITALS ***'
              call optimloc(Cmo,config%decomp%nocc,config%decomp%cfg_mlo_m,ls,config%davidOrbLoc)
           else
               call lsquit('No localization type was requested',ls%lupri) 
	   end if
           ! write localized orbitals
           lun = -1
           CALL LSOPEN(lun,'orbitals_out.u','unknown','UNFORMATTED')
           call mat_write_to_disk(lun,Cmo,OnMaster)
           call LSclose(LUN,'KEEP')
	   call mat_free(cmo)
  end if OnlyLoc

  ! Kasper K, skip Hartree-Fock related calculations for DEC calculation if requested
  ! Also skip, if we only want to localize orbitals
  SkipHF: if( (DECinfo%doDEC .and. .not. DECinfo%doHF) .or. config%davidOrbLoc%OnlyLocalize) then
     write(lupri,*)
     write(lupri,*) 'Initital Hartree-Fock calculation is skipped!'
     write(lupri,*)
     HFdone=.false.
  else
     HFdone=.true.

     call II_precalc_ScreenMat(LUPRI,LUERR,ls%SETTING)

     do_pbc: if(config%latt_config%comp_pbc) then
        CALL mat_init(S,nbast,nbast)

        CALL II_get_overlap(lupri,luerr,ls%setting,S)
        CALL mat_init(D(1),nbast,nbast)
        CALL mat_init(H1,nbast,nbast)
        CALL mat_init(F(1),nbast,nbast)
        CALL II_get_h1(lupri,luerr,ls%setting,H1)
        call get_initial_dens(H1,S,D,ls,config)

        if(.not. config%latt_config%testcase) THEN

          lsint_fock_data%ls => ls
          lsint_fock_data%H1 => H1
          lsint_fock_data%lupri = lupri
          lsint_fock_data%luerr = luerr

          call get_initial_dens(H1,S,D,ls,config)
          config%decomp%S => S 
          do_decomp = .TRUE.
!(config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. &
!               & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
!               & config%decomp%cfg_check_converged_solution .or. &
!               & config%decomp%cfg_rsp_nexcit > 0) 

          if (do_decomp) then
             call decomp_init(nbast,config%decomp)
             call decomposition(config%decomp)
          else if (config%opt%cfg_start_guess == 'TRILEVEL') then
             call mat_free(config%decomp%lcv_CMO)
             config%decomp%decompMatInit_lcv_CMO = .FALSE.
          endif

          if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then !FIXME: put this somewhere else!
             call mat_init(config%av%Fprev,nbast,nbast)
             call mat_init(config%av%Dprev,nbast,nbast)
          endif

          config%opt%cfg_max_linscf_iterations=&
                   config%latt_config%num_its_densmat
          config%opt%opt_quit=.false.
          call scfloop(H1,F,D,S,E,ls,config)

          if (do_decomp) then
             call decomp_shutdown(config%decomp)
          endif

        endif

        Call set_pbc_molecules(ls%input,ls%setting,lupri,luerr,nbast,&
        D(1),config%latt_config)!,config%lib)
        call config_shutdown(config)

     else
        !default - non PBC
        CALL mat_init(S,nbast,nbast)

        CALL II_get_overlap(lupri,luerr,ls%setting,S)
        CALL LSTIMER('*S    ',TIMSTR,TIMEND,lupri)

        IF (config%integral%debugLSlib) THEN 
           print*,'Lslib_debug'
           call Lslib_debug(lupri,luerr,ls%setting,nbast)
        ENDIF
        IF (config%integral%debugGGem) THEN 
           call mem_alloc(GGem,nbast,nbast,nbast,nbast,1) 
           call II_get_GaussianGeminalFourCenter(lupri,luerr,ls%setting,GGem,nbast,.true.) 
           call mem_dealloc(GGem) 
        ENDIF

        CALL mat_init(D(1),nbast,nbast)
        CALL mat_init(F(1),nbast,nbast)
        CALL mat_init(H1,nbast,nbast)

        ! write(lupri,*) 'QQQ New  S:',mat_trab(S,S)
        CALL II_get_h1(lupri,luerr,ls%setting,H1)
        CALL LSTIMER('*H1   ',TIMSTR,TIMEND,lupri)
        ! write(lupri,*) 'QQQ New  H1:',mat_trab(H1,H1)
        !data to pass down to fck_get_fock subroutine
        lsint_fock_data%ls => ls
        lsint_fock_data%H1 => H1
        lsint_fock_data%lupri = lupri
        lsint_fock_data%luerr = luerr

        call get_initial_dens(H1,S,D,ls,config)
        config%decomp%S => S !This is necessary because with high optimization level
        !this pointer is sometimes deassociated

        !debug integral routines
        call di_debug_general(lupri,luerr,ls,nbast,S,D(1),config%integral%debugProp)
        if (mem_monitor) then
           write(lupri,*)
           WRITE(LUPRI,'("Max no. of matrices allocated in Level 2 / get_initial_dens: ",I10)') max_no_of_matrices
           max_no_of_matrices = no_of_matrices
        endif
        !write (lupri,*) 'Start density:'
        !call MAT_PRINT(D, 1, D%nrow, 1, D%ncol, LUPRI)

        CALL LSTIMER('*START',TIMSTR,TIMEND,lupri)

        if (config%opt%cfg_incremental) then
           call ks_init_incremental_fock(nbast)
        endif
        do_decomp =.TRUE.!(config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. &
!             & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
!             & config%decomp%cfg_check_converged_solution .or. &
!             & config%decomp%cfg_rsp_nexcit > 0) 

        if (do_decomp) then
           call decomp_init(nbast,config%decomp)
           call decomposition(config%decomp)
        else if (config%opt%cfg_start_guess == 'TRILEVEL') then
           call mat_free(config%decomp%lcv_CMO)
        endif

        if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then !FIXME: put this somewhere else!
           call mat_init(config%av%Fprev,nbast,nbast)
           call mat_init(config%av%Dprev,nbast,nbast)
        endif

        inquire (file="soeosave.out", exist=soeosaveexist)
        if (config%soeoinp%cfg_restart .and. soeosaveexist) then
          call soeo_restart (config%soeoinp%cfg_unres, config%lupri, S, F(1), D(1))
        else
           if (config%soeoinp%cfg_restart) then
              write (config%lupri, *) 'WARNING: .SOEORST specified but soeosave.out does not exist'
              write (config%lupri, *) '         Making regular optimization to get matrices!'
           endif
           if (config%opt%purescf) then
              !We do one or 2 iterations with a fast and app. exchange/coulomb
              !so we do not need to be converged as hard as level 4
              config%opt%set_convergence_threshold = config%opt%cfg_convergence_threshold*300E0_realk
           endif

           IF(config%integral%LOW_ACCURACY_START)THEN
              call set_Low_accuracy_start_settings(lupri,ls,config,LAStype)
              call scfloop(H1,F,D,S,E,ls,config)
              call revert_Low_accuracy_start_settings(lupri,ls,config,LAStype)
           ENDIF

           call scfloop(H1,F,D,S,E,ls,config)

           !Level 4
           if (config%opt%purescf) then
              Write(lupri,*)'Begin Level 4'
              print*,'Begin Level 4'
              config%opt%set_convergence_threshold = config%opt%cfg_convergence_threshold
              call scf_purify(lupri,ls,config,scfpurify)
              call scfloop(H1,F,D,S,E,ls,config)
           endif
        endif

        IF(config%decomp%cfg_DumpDensRestart)THEN !default true
           ! Kasper K, save Fock and overlap matrices to file for future use
           call save_fock_matrix_to_file(F(1))
           call save_overlap_matrix_to_file(S)
        ENDIF

        if (config%opt%cfg_incremental) call ks_free_incremental_fock()

        if (config%opt%print_final_cmo) then
           call print_orbital_info2(D(1), F(1), S, 'cmo.out', config%decomp%cfg_unres, .true., config%lupri)
        endif

        !SOEO
        if (config%soeoinp%cfg_soeo) then
           write (config%lupri, *) 'Second Order Ensemble Optimization (SOEO) chosen'
           write (config%lupri, *) 'Starting SOEO section:'
           write (config%lupri, *) '================================================'
           call soeoloop (F(1), S, D(1), config%soeoinp)
           write (config%lupri, *) 'SOEO section done'
           write (config%lupri, *) '================================================'
        endif
        !END SOEO

        !debug
!        call mat_init(Cmo,nbast,nbast)
!        allocate(eival(nbast))
!        call mat_diag_f(F,S,eival,Cmo)
!        deallocate(eival)
!        call II_get_AbsoluteValue_overlap(LUPRI,LUERR,ls%SETTING,nbast,CMO,S)
!        call mat_print(S,1,S%nrow,1,S%ncol,lupri)
!        call lsquit('test done',-1)

        !lcm basis
        if (config%decomp%cfg_lcm ) then
           ! get orbitals
           call mat_init(Cmo,nbast,nbast)
           allocate(eival(nbast))
           call mat_diag_f(F(1),S,eival,Cmo)
           deallocate(eival)

           ! write CMO orbitals
           lun = -1
           CALL LSOPEN(lun,'cmo_orbitals.u','unknown','UNFORMATTED')
           call mat_write_to_disk(lun,Cmo,OnMaster)
           call LSclose(LUN,'KEEP')

           ! localize orbitals
           call leastchange_lcm(config%decomp,Cmo,config%decomp%nocc,ls)

           if (config%decomp%cfg_mlo ) then
              write(ls%lupri,'(a)')'Pred= **** LEVEL 3 ORBITAL LOCALIZATION ****'
              call optimloc(Cmo,config%decomp%nocc,config%decomp%cfg_mlo_m,ls,config%davidOrbLoc)
	      if (config%davidOrbLoc%make_orb_plot) then
                 call make_orbitalplot_file(CMO,config%davidOrbLoc,ls)
	      end if
           end if
           ! write LCM orbitals
           lun = -1
           CALL LSOPEN(lun,'lcm_orbitals.u','unknown','UNFORMATTED')
           call mat_write_to_disk(lun,Cmo,OnMaster)
           call LSclose(LUN,'KEEP')

           if (.not. config%decomp%cfg_mlo) then
              call leastchangeOrbspreadStandalone(mx,ls,Cmo,config%decomp%lupri,config%decomp%luerr)
              write(*,*) 'Orbspread standalone LCM: ', mx
           end if

           ! set basis to CMO
           call mat_assign(config%decomp%U_inv,Cmo)
           call mat_mul(config%decomp%U_inv,S,'t','n',1E0_realk,0E0_realk,config%decomp%U)
           ! Vladimir Rybkin: We free CMO unless we do a dec geometry optimization
           ! or run dynamics
           If (.not. (config%optinfo%optimize .OR. config%dynamics%do_dynamics)) then

              ! Single point DEC calculation using current HF files
              DECcalculation: IF(DECinfo%doDEC) then
                 call dec_main_prog_input(ls,F(1),D(1),S,CMO)
              endif DECcalculation
              ! free Cmo
              call mat_free(Cmo)
           Endif
           ! 
        endif
        !
        !  Debug integrals 
        !
        call di_debug_general2(lupri,luerr,ls,nbast,S,D(1))

        !
        ! Optimization
        !
        if (config%optinfo%optimize) then
           if(config%doESGopt)then
              call get_excitation_energy(ls,config,F(1),D(1),S,ExcitE,&
            & config%decomp%cfg_rsp_nexcit)
              Write(lupri,'(A,ES20.9)')'Ground state SCF Energy:',E(1)
              Write(lupri,'(A,ES20.9)')'Excitation Energy      :',ExcitE
              E(1) = E(1) + ExcitE
              Write(lupri,*)'==============================================='
              Write(lupri,'(A,ES20.9)')'Exicted state Energy   :',E(1)
              Write(lupri,*)'==============================================='
           endif
           CALL LS_runopt(E,config,H1,F,D,S,CMO,ls)
           ! Vladimir Rybkin: We free CMO if we have used them
           if (config%decomp%cfg_lcm) then
              ! free Cmo
              call mat_free(Cmo)
           Endif
        endif

        call lsdalton_response(ls,config,F(1),D(1),S)
        
        call config_shutdown(config)

        !
        ! Dynamics
        !
        if (config%dynamics%do_dynamics .EQV. .TRUE.) then
           CALL LS_dyn_run(E,config,H1,F,D,S,CMO,ls)
           ! Vladimir Rybkin: We free CMO if we have used them
           if (config%decomp%cfg_lcm) then
              ! free Cmo
              call mat_free(Cmo)
           Endif           
        endif

        !write(lupri,*) 'mem_allocated_integer, max_mem_used_integer', mem_allocated_integer, max_mem_used_integer

        ! Numerical Derivatives
        if(config%response%tasks%doNumHess .or. &
             & config%response%tasks%doNumGrad .or. &
             & config%response%tasks%doNumGradHess)then 
           nbast=D(1)%nrow
           call get_numerical_hessian(lupri,luerr,ls,nbast,config,config%response%tasks%doNumHess,&
                & config%response%tasks%doNumGrad,config%response%tasks%doNumGradHess)
        endif

        ! Analytical geometrical Hessian
        Natoms = ls%INPUT%MOLECULE%nAtoms
        IF (config%geoHessian%testContrib) THEN
           write (*,*)     'Test the Hessian contributions'
           write (lupri,*) 'Test the Hessian contributions'
           call test_Hessian_contributions(F(1),D(1),Natoms,1,ls%setting,lupri,luerr)
        ENDIF
        IF (config%geoHessian%do_geoHessian) THEN
           call mem_alloc(geomHessian,3*Natoms,3*Natoms)
           write (lupri,*) 'Calculate the Hessian'
           write (*,*)     'Calculate the Hessian'
           call get_molecular_hessian(geomHessian,Natoms,F(1),D(1),ls%setting,config%geoHessian,lupri,luerr)   
           call mem_dealloc(geomHessian)
        ENDIF

        ! PROPERTIES SECTION
        !
        if (config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. & 
             & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
             & config%decomp%cfg_check_converged_solution .or. config%decomp%cfg_rsp_nexcit > 0) then   
           call get_oao_transformed_matrices(config%decomp,F(1),D(1))
        endif
        if (do_decomp) then
           call decomp_shutdown(config%decomp)
           !call dd_shutdown(config%decomp%cfg_unres)
        endif

        CALL LSTIMER('*SCF  ',TIMSTR,TIMEND,lupri)
        WRITE(lupri,*)

        !call mat_no_of_matmuls(matmultot)
        !WRITE(lupri,'("Total no. of matmuls in SCF optimization: ",I10)') matmultot

        if ((config%opt%cfg_start_guess.eq.'TRILEVEL')&
             &.or.(config%opt%cfg_start_guess.eq.'ATOMS')&
             &.or.config%decomp%cfg_gcbasis) then
           !  CALL trilevel_shutdown
           !DEALLOCATE(ls)
        endif

     endif do_pbc
  end if SkipHF


  ! Kasper K, calculate orbital spread for projected atomic orbitals
  if(config%decomp%cfg_PAO) then
     call mat_init(Cmo,nbast,nbast)
     funit = -1
     call lsopen(funit,'cmo_orbitals.u','OLD','UNFORMATTED')
     call mat_read_from_disk(funit,cmo,OnMaster)
     call lsclose(funit,'KEEP')
     if(.not. HFDone) then ! construct overlap matrix from scratch
        call mat_init(S,nbast,nbast)
        CALL II_get_overlap(lupri,luerr,ls%setting,S)
     end if
     call get_PAOs(cmo,S,ls,lupri)
     call mat_free(cmo)
     if(.not. HFDone) call mat_free(S)
  end if

  !
  ! FINALIZE SECTION - releases memory n'stuff
  !
  if(HFdone) then  ! only if HF calculation was carried out 
     CALL mat_free(H1)
     CALL mat_free(F(1))
     CALL mat_free(D(1))
     CALL mat_free(S)
  end if

  call mat_no_of_matmuls(no_of_matmuls)
  WRITE(LUPRI,'("Total no. of matmuls used:                ",I10)') no_of_matmuls
  WRITE(LUPRI,'("Total no. of Fock/KS matrix evaluations:  ",I10)') ls%input%nfock
  if (mem_monitor) WRITE(LUPRI,'("Max no. of matrices allocated in Level 3: ",I10)') max_no_of_matrices


  ! Single point DEC calculation using HF restart files
  ! ***************************************************
  DECcalculationHFrestart: if (.not. HFdone) then
     IF(DECinfo%doDEC) then
        call dec_main_prog_file(ls)
     ENDIF
     call config_shutdown(config)
  endif DECcalculationHFrestart
  call config_free(config)

  call ls_free(ls)

  call lsfree_all()

  call stats_mem(lupri)
  !finalize MPI 
  call lsmpi_finalize(lupri,config%mpi_mem_monitor)
  call print_timers(lupri) !timings for mat operations.
  call LSTIMER('LSDALTON',t1,t2,LUPRI)
  CALL LS_TSTAMP('End simulation',LUPRI)

  CALL LSCLOSE(LUPRI,'KEEP')
  CALL LSCLOSE(LUERR,'KEEP')


END SUBROUTINE LSDALTON

SUBROUTINE lsinit_all()
  use precision
  use matrix_operations, only: MatrixmemBuf_init, set_matrix_default
  use lstensorMem, only: lstmem_init
  use rsp_util, only: init_rsp_util
  use memory_handling, only: init_globalmemvar
  use lstiming, only: init_timers
#ifdef VAR_PAPI
  use papi_module, only: mypapi_init, eventset
#endif
implicit none
  
  ! Init PAPI FLOP counting event using global parameter "eventset" stored in papi_module
#ifdef VAR_PAPI
  call mypapi_init(eventset)
#endif
  call init_globalmemvar  !initialize the global memory counters
  call set_matrix_default !initialize global matrix counters
  call init_rsp_util      !initialize response util module
  call lstmem_init
  call MatrixmemBuf_init()
  call init_timers !initialize timers
  ! MPI initialization
  call lsmpi_init

END SUBROUTINE lsinit_all

SUBROUTINE lsfree_all()
  use precision
  use matrix_operations, only: MatrixmemBuf_free
  use lstensorMem, only: lstmem_free
implicit none
  
  call lstmem_free

  call MatrixmemBuf_free()

END SUBROUTINE lsfree_all