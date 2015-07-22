!> @file 
!> Contains main SCF driver, some module wrappers and miscellaneous 


!> \brief new lsdalton routine to have a nicer stop for the slaves and not to
!run into a "non-beautiful" stop statement
!> \author Patrick Ettenhuber
!> \date 2013
SUBROUTINE LSDALTON
  use precision
  implicit none
  logical     :: OnMaster,meminfo_slaves
  integer     :: lupri, luerr
  real(realk) :: t1,t2
  
  !Set the default values
  OnMaster       = .TRUE.
  meminfo_slaves = .FALSE.
  luerr          = 0
  lupri          = 0

  ! setup the calculation 
  call lsinit_all(OnMaster,lupri,luerr,t1,t2)

  if(OnMaster)then
     ! execute the acutal calculation
     call LSDALTON_DRIVER(OnMaster,lupri,luerr,meminfo_slaves)
  else
     call LSDALTON_DRIVER_SLAVE()
  endif

  ! free everything take time and close the files
  call lsfree_all(OnMaster,lupri,luerr,t1,t2,meminfo_slaves)


END SUBROUTINE LSDALTON


!> \brief Driver for stand-alone f90 linear scaling SCF.
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2008-10-26
SUBROUTINE LSDALTON_DRIVER(OnMaster,lupri,luerr,meminfo_slaves)
  use precision
  use configurationType, only: configitem, LowAccuracyStartType, &
       & set_Low_accuracy_start_settings,revert_Low_accuracy_start_settings
  use TYPEDEFTYPE, only: lsitem
  use matrix_module
  use memory_handling, only: mem_alloc,mem_dealloc, stats_mem, print_memory_info
  use matrix_operations, only: set_matrix_default,max_no_of_matrices, no_of_matrices, &
       & no_of_matmuls, mat_init, mat_free, mat_assign,mat_scal, &
       & mat_mul, mat_no_of_matmuls, mat_write_to_disk, mat_read_from_disk, mat_diag_f,&
       & mat_TrAB, mat_print, mat_tr
  use configuration, only: config_shutdown, config_free
  use files, only: lsopen,lsclose
  use lsdalton_fock_module, only: lsint_fock_data
  use init_lsdalton_mod, only: open_lsdalton_files,init_lsdalton_and_get_lsitem,finalize_lsdalton_driver_and_free
  use initial_guess, only: get_initial_dens
  use scfloop_module, only: scfloop, scf_afterplay
  use lstiming, only: lstimer, init_timers, print_timers
  use ks_settings, only: ks_init_incremental_fock, ks_free_incremental_fock
  use decompMod, only: decomp_init, decomp_shutdown, decomposition, get_oao_transformed_matrices
  use matrix_util, only: save_fock_matrix_to_file, save_overlap_matrix_to_file, util_mo_to_ao_2,read_fock_matrix_from_file
  ! Debug and Testing
  use dal_interface, only: di_debug_general, di_debug_general2
  use extra_output, only: print_orbital_info2
  ! Profile 
#ifdef MOD_UNRELEASED
  use profile_int, only: di_profile_lsint
#endif
  ! DEC 
  use DEC_typedef_module, only: DECinfo  
  ! PROPERTIES SECTION
#ifdef VAR_RSP
  use lsdalton_rsp_mod, only: lsdalton_response, get_excitation_energy
#endif
  use response_noOpenRSP_module, only: lsdalton_response_noOpenRSP
  ! DYNAMICS
  use dynamics_driver, only: LS_dyn_run
  ! SOEO
  use soeo_loop, only: soeoloop, soeo_restart
  ! GEO OPTIMIZER
  use ls_optimizer_mod, only: LS_RUNOPT
  use InteractionEnergyMod, only: InteractionEnergy
  use ADMMbasisOptMod, only: ADMMbasisOptSub
  use lsmpi_type, only: lsmpi_finalize
  use lsmpi_op, only: TestMPIcopySetting,TestMPIcopyScreen
  use lstensorMem, only: lstmem_init, lstmem_free
#ifdef MOD_UNRELEASED
  use pbc_setup, only: set_pbc_molecules
#endif
#ifdef MOD_UNRELEASED
  use numerical_hessian, only: get_numerical_hessian
  use molecular_hessian_mod, only: get_molecular_hessian
  use test_molecular_hessian_mod, only: test_Hessian_contributions
#endif
  use rsp_util, only: init_rsp_util
  use plt_driver_module
  use HODItest_module, only: debugTestHODI
#ifdef VAR_PAPI
  use papi_module
#endif
  use integralinterfaceMod, only: II_get_overlap, II_get_h1, &
       & II_precalc_ScreenMat, II_get_GaussianGeminalFourCenter,&
       & II_get_Fock_mat
  use II_XC_interfaceModule, only: II_get_AbsoluteValue_overlap, &
       & II_get_AbsoluteValue_overlapSame
  use integralinterfaceIchorMod, only: II_Unittest_Ichor,II_Ichor_link_test
  use dec_main_mod!, only: dec_main_prog
  use optimlocMOD, only: optimloc
#ifdef HAS_PCMSOLVER
  use ls_pcm_utils, only: init_molecule
  use ls_pcm_scf, only: ls_pcm_scf_initialize, ls_pcm_scf_finalize
#endif
  implicit none
  logical, intent(in) :: OnMaster
  logical, intent(inout):: meminfo_slaves
  integer, intent(inout) :: lupri, luerr
  integer             :: nbast, lucmo
  TYPE(lsitem),target :: ls
  type(configItem),target  :: config
  real(realk)         :: t1,t2,TIMSTR,TIMEND
  TYPE(Matrix)         :: F(1),D(1), CMO
  Type(Matrix), target :: H1,S
  integer             :: matmultot, lun
  REAL(REALK)         :: mx
  ! Energy
  REAL(REALK)         :: E(1),ExcitE
  logical             :: do_decomp
  real(realk), allocatable :: eival(:)
  real(realk),pointer :: GGem(:,:,:,:,:)
  integer     :: lusoeo,funit
  logical     :: soeosaveexist, HFdone,scfpurify,skipHFpart
  type(matrix) :: Dmo, tmp
  integer             :: nelec
  Integer             :: Natoms
#ifdef MOD_UNRELEASED
  Real(realk),pointer   ::      geomHessian(:,:)
#endif
  type(matrix) :: tempm1,tempm2

  type(LowAccuracyStartType)  :: LAStype

#if VAR_LSDEBUG
  print *,       "THIS IS A DEBUG BUILD"
  write (LUPRI,*)"THIS IS A DEBUG BUILD"
  write (LUERR,*)"THIS IS A DEBUG BUILD"
#endif

  ! Init LSdalton calculation and get lsitem and config structures, and perform
  ! basic tests
  call init_lsdalton_and_get_lsitem(lupri,luerr,nbast,ls,config)

#ifdef HAS_PCMSOLVER
        !
        ! Polarizable continuum model calculation
        !
        if (config%pcm%do_pcm) then
           ! Set molecule object in ls_pcm_utils
           call init_molecule(ls%input%Molecule)
           ! Now initialize PCM
           call ls_pcm_scf_initialize(ls%setting, lupri, luerr)
           write (lupri,*) 'PCMSolver interface correctly initialized'
        end if
#endif        
  ! Timing of individual steps
  CALL LSTIMER('START ',TIMSTR,TIMEND,lupri)
  IF(config%integral%debugIchor)THEN
     call II_unittest_Ichor(LUPRI,LUERR,LS%SETTING,config%integral%debugIchorOption,config%integral%debugIchorLink)
     !the return statement leads to memory leaks but I do not care about this
     !for now atleast
     RETURN
  ENDIF
  IF(config%papitest)THEN
#ifdef VAR_PAPI
     call papi_example(LUPRI)
#endif
  ENDIF
  IF (config%integral%debugUncontAObatch) THEN 
     call II_test_uncontAObatch(lupri,luerr,ls%setting) 
     CALL LSTIMER('II_test_uncontAObatch',TIMSTR,TIMEND,lupri)
  ENDIF

#ifdef MOD_UNRELEASED
  if(config%prof%doProf)then
     call di_profile_lsint(ls,config,lupri,nbast)
     return
  endif
#endif


  ! Skip Hartree Fock part? Done when a HF calculation has already been carried out and we want to:
  ! (i)   localize orbitals
  ! (ii)  carry out DEC calculation 
  ! (iii) Construct PLT file
  if(config%davidOrbLoc%OnlyLocalize .or. (DECinfo%doDEC .and. DECinfo%HFrestart) &
       & .or. config%doplt) then
     skipHFpart=.true.
  else
     skipHFpart=.false.
  end if
  ! Special case, if we run PLT test cases then we do need to run HF calculation first.
  if(config%plt%test) then
     skipHFpart=.false.
  end if


  SkipHF: if(skipHFpart) then   ! Skip Hartree-Fock related calculations
     write(lupri,*)
     write(lupri,*) 'Hartree-Fock calculation is skipped!'
     write(lupri,*)
     HFdone=.false.

  else 
     HFdone=.true.

#ifndef VAR_MPI
     IF(config%doTestMPIcopy)THEN
        !we basicly use the MPICOPY_SETTING routine to place the setting structure
        !in the MPI buffers, deallocate ls%setting and reallocate it again using
        !the MPI buffers - we thereby test some of the functionality of the MPI
        !system. 
        !This Test would (at this moment) break PARI and other testcases because
        !ls%setting%Molecule(1) will nolonger point to ls%input%molecule
        !instead the info in ls%input%molecule will be copied to 
        !ls%setting%Molecule(1),ls%setting%Molecule(2),ls%setting%Molecule(3),...
        !so when you assume a pointer behaviour this test would make the calc crash
        !with something like  
        !Reason: Error in Molecule_free - memory previously released
        call TestMPIcopySetting(ls%SETTING)
        CALL LSTIMER('TestMPIcopySetting',TIMSTR,TIMEND,lupri)
     ENDIF
#endif     
     call II_precalc_ScreenMat(LUPRI,LUERR,ls%SETTING)
     CALL LSTIMER('II_precalc_ScreenMat',TIMSTR,TIMEND,lupri)

#ifndef VAR_MPI
     IF(config%doTestMPIcopy)THEN
        !we basicly use the MPICOPY_SCREEN routine to place the screen structure
        !in the MPI buffers, deallocate screen in screen_mod and reallocate it 
        !again using the MPI buffers - we thereby test some of the 
        !functionality of the MPI system. 
        call TestMPIcopyScreen
     ENDIF
#endif     

     CALL Print_Memory_info(lupri,'after II_precalc_ScreesMat')

#ifdef MOD_UNRELEASED
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


          call mat_scal(2._realk,D(1))

          if (do_decomp) then
             call decomp_shutdown(config%decomp)
          endif

        endif

        write(*,*) 'nearest neighbour = ',config%latt_config%nneighbour
        Call set_pbc_molecules(ls%input,ls%setting,lupri,luerr,nbast,&
        D(1),config%latt_config)!,config%lib)
        call config_shutdown(config)

     else
#endif
        !default - non PBC
        CALL mat_init(S,nbast,nbast)

        CALL II_get_overlap(lupri,luerr,ls%setting,S)
        CALL LSTIMER('*S    ',TIMSTR,TIMEND,lupri)

        CALL Print_Memory_info(lupri,'after II_get_overlap')

        IF (config%integral%debugLSlib) THEN 
           print*,'Lslib_debug'
           call Lslib_debug(lupri,luerr,ls%setting,nbast)
        ENDIF
        IF (config%integral%debugGGem) THEN 
           call mem_alloc(GGem,nbast,nbast,nbast,nbast,1) 
           call II_get_GaussianGeminalFourCenter(lupri,luerr,ls%setting,GGem,nbast,.true.) 
           call mem_dealloc(GGem) 
        ENDIF

        IF (config%doTestHodi) THEN
          call debugTestHODI(lupri,luerr,ls%setting,S,nbast,ls%INPUT%MOLECULE%nAtoms)
          CALL LSTIMER('D-HODI',TIMSTR,TIMEND,lupri)
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
        CALL Print_Memory_info(lupri,'after get_initial_dens')

        !debug integral routines
        call di_debug_general(lupri,luerr,ls,nbast,S,D(1),config%integral%debugProp)

        IF(config%integral%debugIchorLinkFull)THEN
           call II_ichor_LinK_test(lupri,luerr,ls%setting,D)
        ENDIF

        if (config%mat_mem_monitor) then
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

        CALL Print_Memory_info(lupri,'after decomposition')

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

           IF(config%integral%LOW_ACCURACY_START)THEN
              call set_Low_accuracy_start_settings(lupri,ls,config,LAStype)
              call scfloop(H1,F,D,S,E,ls,config)
              call revert_Low_accuracy_start_settings(lupri,ls,config,LAStype)
              CALL Print_Memory_info(lupri,'after Low Accuracy Start')
           ENDIF

           if(config%skipscfloop)then              
              WRITE(config%lupri,*)'The SCF Loop has been skipped!'
              WRITE(config%lupri,*)'Warning: The use of the .SKIPSCFLOOP keyword assumes that the'
              WRITE(config%lupri,*)'fock.restart exist and that it is the final Fock/Kohn-Sham matrix'
              WRITE(config%lupri,*)'of the converged density matrix in dens.restart. Use at own risk.'
              call read_fock_matrix_from_file(F(1))
           else !default
              call scfloop(H1,F,D,S,E,ls,config)
           endif
           CALL Print_Memory_info(lupri,'after scfloop')
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

        IF(config%decomp%cfg_lcm .or. config%decomp%cfg_mlo.or.DECinfo%doDEC) then
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

           IF(config%decomp%debugAbsOverlap)THEN
              call mat_init(tempm1,nbast,nbast)
              call mat_init(tempm2,nbast,nbast)
!              This only works for NOGCBASIS               
!              call mat_mul(S,CMO,'n','n',1E0_realk,0E0_realk,tempm1)
!              call mat_mul(CMO,tempm1,'t','n',1E0_realk,0E0_realk,tempm2)
!              WRITE(lupri,*)'Overlap in MO basis '
!              call mat_print(tempm2,1,nbast,1,nbast,lupri)
              call mat_free(tempm1)
              call II_get_AbsoluteValue_overlapSame(LUPRI,LUERR,ls%SETTING,nbast,nbast,CMO%elms,tempm2%elms)
!              WRITE(lupri,*)'absolute Overlap in MO basis'
!              call mat_print(tempm2,1,nbast,1,nbast,lupri)
              WRITE(lupri,'(A,F16.8)')'Trace of Absolute nummerical overlap:',mat_tr(tempm2)
              WRITE(lupri,'(A,F16.8)')'Trace of Absolute nummerical overlap with CMO:',mat_trAB(CMO,tempm2)
              WRITE(lupri,'(A,F16.8)')'Trace of Absolute nummerical overlap with S:',mat_trAB(S,tempm2)
              call II_get_AbsoluteValue_overlap(LUPRI,LUERR,ls%SETTING,nbast,nbast,nbast,CMO%elms,CMO%elms,tempm2%elms)
!              WRITE(lupri,*)'absolute Overlap in MO basis(test2)  '
!              call mat_print(tempm2,1,nbast,1,nbast,lupri)
              WRITE(lupri,'(A,F16.8)')'Trace of Absolute nummerical overlap(test2):',mat_tr(tempm2)
              WRITE(lupri,'(A,F16.8)')'Trace of Absolute nummerical overlap(test2) with CMO:',mat_trAB(CMO,tempm2)
              WRITE(lupri,'(A,F16.8)')'Trace of Absolute nummerical overlap(test2) with S:',mat_trAB(S,tempm2)
              call mat_free(tempm2)
           ENDIF
        endif

        !lcm basis: localize orbitals
        if (config%decomp%cfg_lcm)THEN
           call leastchange_lcm(config%decomp,Cmo,config%decomp%nocc,ls)
           CALL Print_Memory_info(lupri,'before LCM')
        endif

        !localize orbitals
        if (config%decomp%cfg_mlo ) then
           CALL Print_Memory_info(lupri,'before ORBITAL LOCALIZATION')
           write(ls%lupri,'(a)')'  %LOC%   '
           write(ls%lupri,'(a)')'  %LOC%   ********************************************'
           write(ls%lupri,'(a)')'  %LOC%   *       LEVEL 3 ORBITAL LOCALIZATION       *'
           write(ls%lupri,'(a)')'  %LOC%   ********************************************'
           write(ls%lupri,'(a)')'  %LOC%   '
           call optimloc(Cmo,config%decomp%nocc,config%decomp%cfg_mlo_m,ls,config%davidOrbLoc)
           if (config%davidOrbLoc%make_orb_plot) then
              call make_orbitalplot_file(CMO,config%davidOrbLoc,ls,config%plt)
           end if
           CALL Print_Memory_info(lupri,'after ORBITAL LOCALIZATION')
        endif

        if (config%decomp%cfg_lcm .or. config%decomp%cfg_mlo) then
           !write lcm to file
           lun = -1
           CALL LSOPEN(lun,'lcm_orbitals.u','unknown','UNFORMATTED')
           call mat_write_to_disk(lun,Cmo,OnMaster)
           write(ls%lupri,'(a)') '  %LOC%'
           write(ls%lupri,'(a)') '  %LOC% Localized orbitals written to lcm_orbitals.u'
           write(ls%lupri,'(a)') '  %LOC%'
           call LSclose(LUN,'KEEP')

           if (.not. config%decomp%cfg_mlo) then
              call leastchangeOrbspreadStandalone(mx,ls,Cmo,config%decomp%lupri,config%decomp%luerr)
              write(*,*) 'Orbspread standalone LCM: ', mx
           end if

           ! set basis to CMO
           call mat_assign(config%decomp%U_inv,Cmo)
           call mat_mul(config%decomp%U_inv,S,'t','n',1E0_realk,0E0_realk,config%decomp%U)
        endif

        ! Vladimir Rybkin: We free CMO unless we do a dec geometry optimization
        ! or run dynamics
        If (.not. (config%optinfo%optimize .OR. config%dynamics%do_dynamics)) then
           ! Single point DEC calculation using current HF files
           DECcalculation: IF(DECinfo%doDEC) then
              call dec_main_prog_input(ls,config,F(1),D(1),CMO,E(1))
           endif DECcalculation
           ! free Cmo
           IF(config%decomp%cfg_lcm .or. config%decomp%cfg_mlo.or.DECinfo%doDEC) then
              call mat_free(Cmo)
           ENDIF
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
#ifdef VAR_RSP
              call get_excitation_energy(ls,config,F(1),D(1),S,ExcitE,&
            & config%decomp%cfg_rsp_nexcit)
              Write(lupri,'(A,ES20.9)')'Ground state SCF Energy:',E(1)
              Write(lupri,'(A,ES20.9)')'Excitation Energy      :',ExcitE
              E(1) = E(1) + ExcitE
              Write(lupri,*)'==============================================='
              Write(lupri,'(A,ES20.9)')'Exicted state Energy   :',E(1)
              Write(lupri,*)'==============================================='
              CALL Print_Memory_info(lupri,'after get_excitation_energy')
#else
              call lsquit('Exicted state Energy requires VAR_RSP',lupri)
#endif
           endif
           CALL LS_runopt(E,config,H1,F,D,S,CMO,ls)
           CALL Print_Memory_info(lupri,'after ESG-LS_runopt')
           ! Vladimir Rybkin: We free CMO if we have used them
           if (config%decomp%cfg_lcm) then
              ! free Cmo
              call mat_free(Cmo)
           Endif
        endif

#ifndef VAR_MPI
     IF(config%doTestMPIcopy)THEN
        !we basicly use the MPICOPY_SETTING routine to place the setting structure
        !in the MPI buffers, deallocate ls%setting and reallocate it again using
        !the MPI buffers - we thereby test some of the functionality of the MPI
        !system  
        !This Test would (at this moment) break PARI and other testcases because
        !ls%setting%Molecule(1) will nolonger point to ls%input%molecule
        !instead the info in ls%input%molecule will be copied to 
        !ls%setting%Molecule(1),ls%setting%Molecule(2),ls%setting%Molecule(3),...
        !so when you assume a pointer behaviour this test would make the calc crash
        !with something like  
        !Reason: Error in Molecule_free - memory previously released
        call TestMPIcopySetting(ls%SETTING)
     ENDIF
#endif     

        if (config%InteractionEnergy) then
           CALL InteractionEnergy(E,config,H1,F,D,S,CMO,ls)           
        endif
        if(ls%input%dalton%ADMMBASISFILE)then
           call ADMMbasisOptSub(E,config,H1,F,D,S,CMO,ls)
        endif

        !PROPERTIES SECTION

        if (config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. & 
             & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
             & config%decomp%cfg_check_converged_solution .or. config%decomp%cfg_rsp_nexcit > 0) then   
           call get_oao_transformed_matrices(config%decomp,F(1),D(1))
        endif

        IF(config%response%tasks%doDipole.OR.config%response%tasks%doResponse)THEN
           CALL Print_Memory_info(lupri,'before lsdalton_response')
           IF(config%response%noOpenRSP)THEN
              !A pure LSDALTON
              call lsdalton_response_noOpenRSP(ls,config,F,D,S)
           ELSE
#ifdef VAR_RSP
              call lsdalton_response(ls,config,F(1),D(1),S)
#else
              call lsquit('lsdalton_response requires VAR_RSP defined',-1)
#endif
           ENDIF
           CALL Print_Memory_info(lupri,'after lsdalton_response')
        ENDIF

        call config_shutdown(config)

        !
        ! Dynamics
        !
        if (config%dynamics%do_dynamics .EQV. .TRUE.) then
           CALL Print_Memory_info(lupri,'before dynamics - LS_dyn_run')
           CALL LS_dyn_run(E,config,H1,F,D,S,CMO,ls)
           ! Vladimir Rybkin: We free CMO if we have used them
           if (config%decomp%cfg_lcm) then
              ! free Cmo
              call mat_free(Cmo)
           Endif           
           CALL Print_Memory_info(lupri,'after dynamics - LS_dyn_run')
        endif

        !write(lupri,*) 'mem_allocated_integer, max_mem_used_integer', mem_allocated_integer, max_mem_used_integer

#ifdef MOD_UNRELEASED
        ! Numerical Derivatives
        if(config%response%tasks%doNumHess .or. &
             & config%response%tasks%doNumGrad .or. &
             & config%response%tasks%doNumGradHess)then
           CALL Print_Memory_info(lupri,'before get_numerical_hessian')
           nbast=D(1)%nrow
           call get_numerical_hessian(lupri,luerr,ls,nbast,config,config%response%tasks%doNumHess,&
                & config%response%tasks%doNumGrad,config%response%tasks%doNumGradHess)
           CALL Print_Memory_info(lupri,'after get_numerical_hessian')
        endif

        ! Analytical geometrical Hessian
        Natoms = ls%INPUT%MOLECULE%nAtoms
        IF (config%geoHessian%testContrib) THEN
           write (*,*)     'Test the Hessian contributions'
           write (lupri,*) 'Test the Hessian contributions'
           call test_Hessian_contributions(F(1),D(1),Natoms,1,ls%setting,lupri,luerr)
        ENDIF
        IF (config%geoHessian%do_geoHessian) THEN
           CALL Print_Memory_info(lupri,'before get_molecular_hessian')
           call mem_alloc(geomHessian,3*Natoms,3*Natoms)
           write (lupri,*) 'Calculate the Hessian'
           write (*,*)     'Calculate the Hessian'
           call get_molecular_hessian(geomHessian,Natoms,F(1),D(1),ls%setting,config,lupri,luerr)   
           call mem_dealloc(geomHessian)
           CALL Print_Memory_info(lupri,'after get_molecular_hessian')
        ENDIF
#endif

        ! END OF PROPERTIES SECTION
        !
        if (do_decomp) then
           call decomp_shutdown(config%decomp)
           !call dd_shutdown(config%decomp%cfg_unres)
        endif

        CALL LSTIMER('*SCF  ',TIMSTR,TIMEND,lupri,.TRUE.)
        WRITE(lupri,*)

        !call mat_no_of_matmuls(matmultot)
        !WRITE(lupri,'("Total no. of matmuls in SCF optimization: ",I10)') matmultot

        if ((config%opt%cfg_start_guess.eq.'TRILEVEL')&
             &.or.(config%opt%cfg_start_guess.eq.'ATOMS')&
             &.or.config%decomp%cfg_gcbasis) then
           !  CALL trilevel_shutdown
           !DEALLOCATE(ls)
        endif

#ifdef MOD_UNRELEASED
     endif do_pbc
#endif
  end if SkipHF


  ! Read in already optimized HF orbitals, and localize
  OnlyLoc:  if (config%davidOrbLoc%OnlyLocalize) then
           CALL Print_Memory_info(lupri,'before LOCALIZING ORBITALS L3')
           ! read orbitals
           lun = -1
           call mat_init(CMO,nbast,nbast)
           CALL LSOPEN(lun,'orbitals_in.u','unknown','UNFORMATTED')
           call mat_read_from_disk(LUN,cmo,OnMaster)
           call LSclose(LUN,'KEEP')
           if (config%decomp%cfg_mlo) then
              write(ls%lupri,*)
              call optimloc(Cmo,config%decomp%nocc,config%decomp%cfg_mlo_m,ls,config%davidOrbLoc)
           else
               call lsquit('No localization type was requested',ls%lupri) 
           end if
           ! write localized orbitals
           lun = -1
           CALL LSOPEN(lun,'localized_orbitals.u','unknown','UNFORMATTED')
           call mat_write_to_disk(lun,Cmo,OnMaster)
           write(ls%lupri,'(a)') '  %LOC%'
           write(ls%lupri,'(a)') '  %LOC% Localized orbitals written to localized_orbitals.u'
           write(ls%lupri,'(a)') '  %LOC%'
           call LSclose(LUN,'KEEP')
           call mat_free(cmo)
           CALL Print_Memory_info(lupri,'after LOCALIZING ORBITALS L3')
  end if OnlyLoc

  ! Construct PLT file
  ConstructPLT: if(config%doplt) then
     call plt_wrapper(ls,config%plt)
  end if ConstructPLT


  ! Single point DEC calculation using HF restart files
  DECcalculationHFrestart: if ( (DECinfo%doDEC .and. DECinfo%HFrestart) ) then
     CALL Print_Memory_info(lupri,'before dec_main_prog_file')
     call dec_main_prog_file(ls,config)
     CALL Print_Memory_info(lupri,'after dec_main_prog_file')
  endif DECcalculationHFrestart


  ! Kasper K, calculate orbital spread for projected atomic orbitals
  if(config%decomp%cfg_PAO) then
     CALL Print_Memory_info(lupri,'before cfg_PAO')
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
     CALL Print_Memory_info(lupri,'after cfg_PAO')
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
  if (config%mat_mem_monitor) WRITE(LUPRI,'("Max no. of matrices allocated in Level 3: ",I10)') max_no_of_matrices

  if(.not. HFdone) then  ! ensure there's no memory leak when HF calc was skipped
     call config_shutdown(config)
  end if
  call config_free(config)

#ifdef HAS_PCMSOLVER
  if (config%pcm%do_pcm) then
     ! Now finalize PCM
     call ls_pcm_scf_finalize
     write (lupri,*) 'PCMSolver interface correctly finalized'
  end if
#endif

  ! Free all higher structures used in LSALTON_DRIVER, and which could not be
  ! initialized before the config was read
  call finalize_lsdalton_driver_and_free(lupri,luerr,ls,config,meminfo_slaves)

END SUBROUTINE LSDALTON_DRIVER

SUBROUTINE LSDALTON_DRIVER_SLAVE()
   use lsmpi_type, only: MPI_COMM_LSDALTON
   implicit none
#ifdef VAR_MPI
   call lsmpi_slave(MPI_COMM_LSDALTON)
#endif
END SUBROUTINE LSDALTON_DRIVER_SLAVE


SUBROUTINE lsinit_all(OnMaster,lupri,luerr,t1,t2)
  use precision
  use matrix_operations, only: set_matrix_default
  use init_lsdalton_mod, only: open_lsdalton_files
  use lstensorMem, only: lstmem_init
  use rsp_util, only: init_rsp_util
  use dft_memory_handling
  use memory_handling, only: init_globalmemvar
  use lstiming, only: init_timers, lstimer,  print_timers,time_start_phase,PHASE_WORK
  use tensor_interface_module,only: tensor_initialize_interface
  use GCtransMod, only: init_AO2GCAO_GCAO2AO
  use IntegralInterfaceModuleDF,only:init_IIDF_matrix
#ifdef VAR_PAPI
  use papi_module, only: mypapi_init, eventset
#endif
  use gpu_device_handling 
#ifdef VAR_ICHOR
  use IchorSaveGabMod
#endif
  use lsmpi_type,only: NullifyMPIbuffers
  implicit none
  logical, intent(inout)     :: OnMaster
  integer, intent(inout)     :: lupri, luerr
  real(realk), intent(inout) :: t1,t2
  
  !INITIALIZING TIMERS SHOULD ALWAYS BE THE FIRST CALL
  call init_timers

  ! Init PAPI FLOP counting event using global parameter "eventset" stored in papi_module
#ifdef VAR_PAPI
  call mypapi_init(eventset)
#endif

  call Init_GPU_devices   !initialize gpu(s) (acc_init) 

  call init_globalmemvar  !initialize the global memory counters
  call NullifyMPIbuffers  !initialize the MPI buffers
  call set_matrix_default !initialize global matrix counters
  call init_rsp_util      !initialize response util module
  call lstmem_init
  call setPrintDFTmem(.FALSE.)
  call init_IIDF_matrix
#ifdef VAR_ICHOR
  call InitIchorSaveGabModule()
#endif
  call init_AO2GCAO_GCAO2AO()
  ! MPI initialization
  call lsmpi_init(OnMaster)
  !tensor initialization
  call tensor_initialize_interface()

  !INIT TIMING AND FILES
  if(OnMaster)then
    call LSTIMER('START',t1,t2,LUPRI)
    call open_lsdalton_files(lupri,luerr)
  endif

  call time_start_phase(PHASE_WORK)
END SUBROUTINE lsinit_all

SUBROUTINE lsfree_all(OnMaster,lupri,luerr,t1,t2,meminfo)
  use precision
  use memory_handling, only: stats_mem
  use files, only: lsclose
  use lstiming, only: lstimer, init_timers, print_timers
  use lstensorMem, only: lstmem_free
  use tensor_interface_module ,only: tensor_finalize_interface
  use GCtransMod, only: free_AO2GCAO_GCAO2AO
  use IntegralInterfaceModuleDF,only:free_IIDF_matrix
  use dec_settings_mod, only:free_decinfo
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif
  use gpu_device_handling, only: shutdown_GPU_devices 
#ifdef VAR_ICHOR
  use IchorSaveGabMod
#endif
#ifdef VAR_SCALAPACK
  use matrix_operations_scalapack
#endif
  use matrix_operations_pdmm
  implicit none
  logical,intent(in)         :: OnMaster
  integer,intent(inout)      :: lupri,luerr
  logical,intent(inout)      :: meminfo
  real(realk), intent(inout) :: t1,t2

  call shutdown_GPU_devices   !shut down the device (acc_shutdown)
  
  IF(OnMaster)THEN
     !these routines free matrices and must be called while the 
     !slaves are still in mpislave routine
     call free_IIDF_matrix
     call free_AO2GCAO_GCAO2AO()
  ENDIF
  call lstmem_free
#ifdef VAR_ICHOR
  if(OnMaster)call FreeIchorSaveGabModule()
#endif
  

  !IF MASTER ARRIVED, CALL THE SLAVES TO QUIT AS WELL
#ifdef VAR_MPI
  if(OnMaster)call ls_mpibcast(LSMPIQUIT,infpar%master,MPI_COMM_LSDALTON)
#endif  

  call tensor_finalize_interface()
  call free_decinfo()


  if(OnMaster) call stats_mem(lupri)

#ifdef VAR_MPI
  if( infpar%parent_comm==MPI_COMM_NULL ) then

    call ls_mpibcast(meminfo,infpar%master,MPI_COMM_LSDALTON)
    if(meminfo)call lsmpi_print_mem_info(lupri,.false.)

  endif

#ifdef VAR_SCALAPACK
  IF(scalapack_mpi_set)THEN
     !free communicator 
     call LSMPI_COMM_FREE(scalapack_comm)
  ENDIF
#endif
  IF(pdmm_mpi_set)THEN
     !free communicator 
     call LSMPI_COMM_FREE(pdmm_comm)
  ENDIF

  call lsmpi_finalize(lupri,.false.)
#else
  WRITE(LUPRI,'("  Allocated MPI memory a cross all slaves:  ",i9," byte  &
       &- Should be zero - otherwise a leakage is present")') 0
  WRITE(LUPRI,'(A)')'  This is a non MPI calculation so naturally no memory is allocated on slaves!'
#endif

  if(OnMaster)then

    call print_timers(lupri) !timings for mat operations.
    call LSTIMER('LSDALTON',t1,t2,LUPRI,.TRUE.)
    CALL LS_TSTAMP('End simulation',LUPRI)

    CALL LSCLOSE(LUPRI,'KEEP')
    CALL LSCLOSE(LUERR,'KEEP')

  endif
END SUBROUTINE lsfree_all
