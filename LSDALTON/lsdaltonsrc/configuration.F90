!> @file
!> Contains configuration module. 

!> \brief Type definitions for configuration structures and routines for reading input.
!> \author S. Host and T. Kjaergaard
!> \date March 2010
module configuration
use precision
use lstiming, only: SET_LSTIME_PRINT
use configurationType, only: configitem
use profile_type, only: profileinput, prof_set_default_config
use tensor_interface_module, only: tensor_set_dil_backend_true, &
   &tensor_set_debug_mode_true, tensor_set_always_sync_true,lspdm_init_global_buffer, &
   tensor_set_global_segment_length
#ifdef MOD_UNRELEASED
use typedeftype, only: lsitem,integralconfig,geoHessianConfig
#else
use typedeftype, only: lsitem,integralconfig
#endif
use opttype, only: opt_set_default_config
use response_wrapper_type_module, only: free_mcdinputitem, &
     & alphainputitem_set_default_config, betainputitem_set_default_config, &
     & gammainputitem_set_default_config, tpainputitem_set_default_config, &
     & dtpainputitem_set_default_config, esginputitem_set_default_config, &
     & esdinputitem_set_default_config, mcdinputitem_set_default_config, &
     & rspsolveriputitem_set_default_config
use lsdalton_response_type_mod, only: rsp_tasks_set_default_config
use soeo_typedef,only: soeoinp_set_default_config
!use matrix_module!, only: matrix
use typedef, only: integral_set_default_config
use arhDensity, only: arh_set_default_config
use av_utilities, only: av_set_default_config, av_shutdown
use diagonalization, only: diag_lshift_dorth, &
     & diag_lshift_search, diag_lshift_none, diag_lshift_vanlenthe, &
     & diag_set_default_config
use decompMod, only: decomp_set_default_config
use files, only: lsopen,lsclose
use matrix_operations, only: mat_select_type, matrix_type, &
     & mtype_symm_dense, mtype_dense, &
     & mtype_unres_dense, mtype_csr, mtype_scalapack
use matrix_operations_aux, only: mat_zero_cutoff, mat_inquire_cutoff
use DEC_settings_mod, only: dec_set_default_config, config_dec_input,&
     & check_cc_input, check_dec_input
use dec_typedef_module,only: DECinfo,MODEL_MP2,MODEL_CCSDpT,MODEL_RIMP2
use optimization_input, only: optimization_set_default_config, ls_optimization_input
use ls_dynamics, only: ls_dynamics_init, ls_dynamics_input
#ifdef MOD_UNRELEASED
  use lattice_vectors, only: pbc_setup_default
#endif
use localization_input
use davidson_settings, only: davidson_default_SCF, davidson_default
use molecule_type, only: free_moleculeinfo
use readmolefile, only: read_molfile_and_build_molecule
use IntegralInterfaceMOD, only: ii_get_nucpot
use ks_settings, only: ks_free_incremental_fock
use memory_handling, only: mem_alloc,mem_dealloc
use dft_typetype
use dft_memory_handling
use plt_driver_module
#ifdef VAR_MPI
use infpar_module
use lsmpi_mod
use lsmpi_type, only: DFTSETFU,SPLIT_MPI_MSG,MAX_SIZE_ONE_SIDED
#endif
#ifdef BUILD_CGTODIFF
use cgto_diff_eri_host_interface, only: cgto_diff_eri_xfac_general
#endif
use scf_stats, only: scf_stats_arh_header
#ifdef MOD_UNRELEASED
use molecular_hessian_mod, only: geohessian_set_default_config
#endif
use xcfun_host,only: xcfun_host_init, USEXCFUN, XCFUNDFTREPORT
use ls_util,only: capitalize_string
#ifdef HAS_PCMSOLVER
use ls_pcm_config, only: pcmtype, ls_pcm_init, ls_pcm_input
#endif
use LSparameters
use iidftksmwork, only: dft_dogga_dometa

private
public :: config_set_default_config, config_read_input, config_shutdown,&
     & config_free, set_final_config_and_print
contains

!> \brief Call routines to set default values for different structures.
!> \author S. Host
!> \date March 2010
subroutine config_set_default_config(config)
implicit none
   !> Contains info, settings and data for entire calculation
   type(ConfigItem), intent(inout) :: config
  nullify(config%solver)
  allocate(config%solver)
  call arh_set_default_config(config%solver)
  nullify(config%decomp)
  allocate(config%decomp)
  call decomp_set_default_config(config%decomp)
  call integral_set_default_config(config%integral)
#ifdef VAR_XCFUN
  USEXCFUN = .TRUE.  
#else
  USEXCFUN = .FALSE.  
#endif
  call opt_set_default_config(config%opt)
  call av_set_default_config(config%av)
  !DEC calculation. The DEC info is kept as a global derived type parameter for now.
  call dec_set_default_config(config%lupri)
  call diag_set_default_config(config%diag)
  call soeoinp_set_default_config(config%soeoinp)
  ! Orbital localization settings
  call davidson_default(config%davidOrbLoc)
  ! SCF optimization w/ davidson solver settings
  call davidson_default(config%davidSCF)
  call davidson_default_SCF(config%davidSCF)
  !RESPONSE
  ! Polarizability
  call ALPHAinputitem_set_default_config(config%response%alphainput)
  ! 1st hyperpolarizability
  call BETAinputitem_set_default_config(config%response%betainput)
  ! 2nd hyperpolarizability
  call GAMMAinputitem_set_default_config(config%response%gammainput)
  ! Standard TPA
  call TPAinputitem_set_default_config(config%response%tpainput)
  ! Damped TPA
  call DTPAinputitem_set_default_config(config%response%dtpainput)
  ! ESG
  call ESGinputitem_set_default_config(config%response%esginput)
  ! ESD
  call ESDinputitem_set_default_config(config%response%esdinput)
  ! MCD
  call MCDinputitem_set_default_config(config%response%MCDinput)
  ! RSP solver
  call RSPSOLVERiputitem_set_default_config(config%response%RSPSOLVERinput)
  call rsp_tasks_set_default_config(config%response%tasks)
#ifdef VAR_RSP
  config%response%noOpenRSP = .FALSE. !Use OpenRSP module - Default
#else
  config%response%noOpenRSP = .TRUE.  !Use LSDALTON own Response module
#endif
#ifdef MOD_UNRELEASED
  ! Molecular Hessian
  call geohessian_set_default_config(config%geoHessian)
#endif
  ! geometry optimization
  call optimization_set_default_config(config%optinfo)
  ! Dynamics
  call LS_dynamics_init(config%dynamics)
#ifdef HAS_PCMSOLVER
  ! Polarizable Continuum Model
  call ls_pcm_init(config%pcm)
#endif
  !Only for testing new sparse matrix library, should be removed afterwards!
  config%sparsetest = .false.
  config%mpi_mem_monitor = .false.
  config%doDEC = .false.
  config%InteractionEnergy = .false.
  config%access_stream = .false.
  config%SameSubSystems = .false.
  config%SubSystemDensity = .false.
  config%PrintMemory = .false.
  config%doESGopt = .false.
  config%GPUMAXMEM = 2.0E0_realk
  config%noDecEnergy = .false.
  call prof_set_default_config(config%prof)
#ifdef MOD_UNRELEASED
  call pbc_setup_default(config%latt_config)
#endif
  ! PLT info
  call pltinfo_set_default_config(config%Plt)
  config%doplt         = .false.
  !F12 calc?
  config%doF12         = .false.
  config%doRIMP2       = .false.
  config%doTestMPIcopy = .false.
  config%doTestHodi    = .false.
  config%skipscfloop   = .false.
  config%papitest      = .false.
#ifdef VAR_MPI
  infpar%inputBLOCKSIZE = 0
  print*,'config_set_default_config:',infpar%inputBLOCKSIZE
#endif
end subroutine config_set_default_config

!> \brief Wrapper to routines for read input files LSDALTON.INP and MOLECULE.INP.
!> \author S. Host and T. Kjaergaard
!> \date March 2010
subroutine config_read_input(config,lupri,luerr)
implicit none
   !> Contains info, settings and data for entire calculation
   type(ConfigItem), intent(inout) :: config
   !> Logical unit number for LSDALTON.OUT
   integer, intent(in)             :: lupri
   !> Logical unit number for LSDALTON.ERR
   integer, intent(in)             :: luerr
   LOGICAL :: LSint,nonlsdaltonrun

   config%lupri          = lupri
   config%luerr          = luerr
   config%solver%lupri   = lupri
   config%decomp%lupri   = lupri
   config%decomp%luerr   = luerr
   config%av%lupri       = lupri
   config%opt%lupri      = lupri
   config%opt%luerr      = luerr
   config%diag%lupri     = lupri
   config%soeoinp%lupri  = lupri

   !read the MOLECULE.INP and set input
   call read_molfile_and_build_molecule(lupri,config%molecule,config%LIB,&
        & .FALSE.,0,config%integral%DoSpherical,config%integral%basis,&
        & config%latt_config,config%integral%atombasis)
   config%integral%nelectrons = config%molecule%nelectrons 
   config%integral%molcharge = INT(config%molecule%charge)
   !read the LSDALTON.INP and set input
   CALL read_dalton_input(LUPRI,config) 

end subroutine config_read_input

!> \brief If data has been allocated by the default setting, it should be free'd here.
!> \author S. Host
!> \date March 2010
subroutine config_shutdown(config)
implicit none
   !> Contains info, settings and data for entire calculation
   type(ConfigItem), intent(inout) :: config

   call av_shutdown(config%av)
   call free_Moleculeinfo(config%Molecule)
   call free_MCDinputitem(config%response%MCDinput)
end subroutine config_shutdown

subroutine config_free(config)
implicit none
   !> Contains info, settings and data for entire calculation
   type(ConfigItem), intent(inout) :: config
   deallocate(config%solver)
   nullify(config%solver)
   deallocate(config%decomp)
   nullify(config%decomp)

end subroutine config_free

!> \brief Read input file LSDALTON.INP and set configuration structure accordingly.
!> \author T. Kjaergaard
!> \date March 2010
SUBROUTINE read_dalton_input(LUPRI,config)
! READ THE INPUT FOR THE INTEGRAL 
use IIDFTINT, only: II_DFTsetFunc
#if defined(ENABLE_QMATRIX)
use ls_qmatrix, only: ls_qmatrix_init, &
                      ls_qmatrix_input
#endif

implicit none
!> Logical unit number for LSDALTON.OUT
INTEGER            :: LUPRI
!> Contains info, settings and data for entire calculation
type(ConfigItem), intent(inout) :: config
INTEGER            :: LUCMD !Logical unit number for the daltoninput
INTEGER            :: IDUMMY,IPOS,IPOS2,IPOS3,IPOSMU,COUNTER
character(len=80)  :: WORD,TMPWORD,camalpha,cambeta,cammu
character(len=1024):: Func
character(len=2)   :: PROMPT
LOGICAL            :: DONE,file_exists,READWORD,LSDALTON,STARTGUESS,WAVE,exchangescale
!LINSCA variables:
real(realk)        :: shift, min_density_overlap, maxratio, zero
integer            :: nvec, i
Real(realk)        :: hfweight 

WAVE = .FALSE.
exchangescale = .FALSE.

STARTGUESS = .FALSE.
Config%integral%cfg_lsdalton = .TRUE.
COUNTER = 0

INQUIRE(file='LSDALTON.INP',EXIST=file_exists) 
IF(file_exists)THEN
   LUCMD=-1
   CALL lsOPEN(LUCMD,'LSDALTON.INP','OLD','FORMATTED')
ELSE
   CALL lsQUIT('LSDALTON.INP does not exist',lupri)
ENDIF
READWORD=.TRUE.
DONE=.FALSE.
rewind(LUCMD)
DO
   IF(DONE)EXIT
   IF(READWORD) THEN
      READ (LUCMD, '(A40)') WORD
      READWORD=.TRUE.
      COUNTER = 0
   ELSE
      IF (COUNTER.GT. 1) THEN
        WRITE(LUPRI,'(1X,2A)') 'Infinite loop for input line:',WORD
        CALL lsQUIT('Infinite loop in read_dalton_input,due to wrong input',lupri)
      ENDIF
      COUNTER = COUNTER + 1
   ENDIF
!
   PROMPT = WORD(1:2)
   IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#'))CYCLE
!   IF (WORD(1:14) == '**DALTON INPUT') CYCLE
   IF (WORD(1:10) == '**GENERAL') THEN
      READWORD = .TRUE.
      CALL GENERAL_INPUT(config,readword,word,lucmd,lupri)
   ENDIF
   IF (WORD(1:10) == '**PROFILE') THEN
      READWORD = .TRUE.
      CALL PROFILE_INPUT(config%prof,readword,word,lucmd,lupri)
   ENDIF
   IF (WORD(1:10) == '**INTEGRAL') THEN
      READWORD = .TRUE.
      IF (WAVE) CALL LSQUIT('**INTEGRAL section must be defined before **WAVE FUN in the lsdalton input',lupri)
      CALL INTEGRAL_INPUT(config%integral,readword,word,lucmd,lupri)
   ENDIF
   IF (WORD(1:6) == '**WAVE') THEN
      READWORD=.TRUE.
      DO   
         IF(READWORD) THEN
            READ (LUCMD, '(A80)') WORD
            READWORD=.TRUE.
         ENDIF
         PROMPT = WORD(1:2)
         IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
         IF(PROMPT(1:1) == '*'.OR.PROMPT(1:1) == '$')THEN
            READWORD=.FALSE.
            EXIT
         ENDIF
         IF(PROMPT(1:1) .EQ. '.') THEN
            SELECT CASE(WORD(1:13)) 
            CASE ('.HF');  config%opt%calctype = config%opt%hfcalc !Hartree-Fock calc
                     config%integral%exchangeFactor = 1E0_realk
                     config%integral%dft%HFexchangeFac = 1E0_realk
            CASE ('.Hartree');  config%opt%calctype = config%opt%hfcalc !Hartree calc
                     config%integral%exchangeFactor = 0E0_realk
                     config%integral%dft%HFexchangeFac = 0E0_realk
            CASE ('.EXCHANGESCAL');
               exchangescale = .TRUE.
               READ(LUCMD,*) config%integral%exchangeFactor
               config%integral%dft%HFexchangeFac = config%integral%exchangeFactor
            CASE ('.DFT'); 
               call read_dft_input(config,lucmd,lupri)
            CASE DEFAULT
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                    & '" not recognized in **WAVE FUNCTION'
               CALL lsQUIT('Illegal keyword in **WAVE FUNCTION',lupri)
            END SELECT
         ELSE
            WRITE(LUPRI,'(2X,A)')'Requires .DFT or .HF after **WAVEFUNCTION'
            CALL lsQUIT('Requires .DFT or .HF after **WAVEFUNCTION',lupri)
         ENDIF
      ENDDO
      IF(READWORD) THEN
         READ (LUCMD, '(A80)') WORD
         READWORD=.TRUE.
      ENDIF
   ENDIF
   IF (WORD(1:10) == '*DFT INPUT') THEN
      CALL READ_WAVE_DFTINPUT(LUPRI,LUCMD,config%integral,WORD)
      READWORD=.FALSE.
!      EXIT !done with what we are interested in
   ENDIF
   IF (WORD(1:7) == '*LINSCA') THEN
      write(lupri,*) '*LINSCA under **WAVE FUNCTION has been replaced by *DENSOPT'
      call lsquit('*LINSCA under **WAVE FUNCTION has been replaced by *DENSOPT',-1)
   else if (WORD(1:8) == '*DENSOPT') THEN
      READWORD=.TRUE.
      DO   
         IF(READWORD) THEN
            READ (LUCMD, '(A40)') WORD
            READWORD=.TRUE.
         ENDIF
         PROMPT = WORD(1:2)
         IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
         IF(PROMPT(1:1) .EQ. '.') THEN
            SELECT CASE(WORD) 
            CASE('.2ND_ALL');    config%solver%cfg_2nd_order_all = .true.
                                 config%solver%cfg_do_2nd_order = .true.
                                 config%solver%set_do_2nd_order = .true.
                                 config%opt%CFG_density_method =  config%opt%CFG_F2D_DIRECT_DENS
                                 config%solver%cfg_arh_truncate = .false.
            CASE('.2ND_LOC');    config%solver%cfg_2nd_order_local = .true.
                                 config%solver%cfg_do_2nd_order = .true.
                                 config%opt%CFG_density_method =  config%opt%CFG_F2D_DIRECT_DENS
            CASE('.ARH');        config%solver%cfg_arh_truncate = .true.
                                 config%solver%cfg_arh_crop = .true.
                                 config%opt%CFG_density_method = config%opt%CFG_F2D_ARH
            CASE('.ARH FULL');   config%solver%cfg_arh_crop = .true.
                                 config%solver%cfg_arh_truncate = .false.
                                 config%opt%CFG_density_method = config%opt%CFG_F2D_ARH
            CASE('.ARH DAVID');  config%davidSCF%arh_davidson=.true.
                                 config%davidSCF%arh_lintrans = .true.
                                 config%davidSCF%precond=.true.
                                 config%opt%cfg_saveF0andD0 = .true.
                                 config%davidSCF%stepsize=0.5
                                 config%davidSCF%arh_inp_linesearch=.false.
                                 config%davidSCF%max_stepsize = config%davidSCF%stepsize
            CASE('.ARH(LS) DAVID');  config%davidSCF%arh_davidson=.true.
                                 config%davidSCF%arh_lintrans = .true.
                                 config%davidSCF%precond=.true.
                                 config%opt%cfg_saveF0andD0 = .true.
                                 config%davidSCF%stepsize=0.5
                                 config%davidSCF%arh_inp_linesearch=.true.
                                 config%davidSCF%max_stepsize = config%davidSCF%stepsize
            CASE('.ARH DEBUG');  config%davidSCF%arh_davidson_debug=.true.
            CASE('.DAVIDSON DEBUG'); config%davidSCF%debug_info =.true.
            CASE('.DAVIDSON EXTRAVEC'); config%davidSCF%arh_extravec =.true.
                                      config%davidSCF%arh_inp_extravec =.true.
            CASE('.DAVIDSON LSDEBUG'); config%davidSCF%arh_debug_linesearch =.true.
            CASE('.NOECONTINCREM');
               IF(.NOT.config%opt%cfg_saveF0andD0)THEN
                call lsquit('.NOECONTINCREM must be placed some pointer after .ARH DAVID',-1)
               ENDIF
               config%opt%cfg_saveF0andD0 = .false.
!#ifdef MOD_UNRELEASED
!            CASE('.ASYM');       config%opt%cfg_asym = .true.
!#endif
            CASE('.CHOLESKY');   config%decomp%lowdin_diagonalize = .false.; config%decomp%cholesky_decomp   = .true.
            CASE('.CONFSHIFT');  config%diag%cfg_no_conf_shift = .false.
            CASE('.CONTFAC');    READ(LUCMD,*) config%solver%cfg_arh_contract
            CASE('.CONTRAC');    READ(LUCMD,*) config%solver%cfg_arh_contract_crit
            CASE('.CONVDYN');    READ(LUCMD,*) config%opt%cfg_convdyn_type ; config%opt%cfg_convdyn = .true.
            CASE('.CONVTHR');    READ(LUCMD,*) config%opt%cfg_convergence_threshold
                                 config%opt%set_convergence_threshold = config%opt%cfg_convergence_threshold
            CASE('.DIAGHESONLY'); config%opt%cfg_diaghesonly = .true.
            CASE('.DIIS');       config%av%CFG_averaging = config%av%CFG_AVG_DIIS
            CASE('.DISK');       config%opt%cfg_queue_on_disk = .true.
            CASE('.DISKSOLVER'); config%solver%cfg_arh_disk_micro = .true.
            !CASE('.DISKQUEUE') ; config%solver%cfg_arh_disk_macro = .true. !Not active - get_from_modFIFO_disk won't work!
            CASE('.DORTH');      config%diag%CFG_lshift = diag_lshift_dorth
                                 config%av%CFG_lshift = diag_lshift_dorth
            CASE('.DUMPMAT');    config%opt%dumpmatrices = .true.
!removed keyword - no testcase - and naturally it does not seem to work. TK
!            CASE('.EDIIS');      config%av%CFG_averaging = config%av%CFG_AVG_EDIIS
            CASE('.EXPAND');     READ(LUCMD,*) config%solver%cfg_arh_expand_crit
            CASE('.EXPFAC');     READ(LUCMD,*) config%solver%cfg_arh_expand 
            CASE('.FIXSHIFT');   READ(LUCMD,*) shift 
                                 config%solver%cfg_fixed_shift_param = shift ; config%solver%cfg_fixed_shift = .true.
                                 config%diag%cfg_fixed_shift_param = shift   ; config%diag%cfg_fixed_shift = .true.
            CASE('.FLUSH');      config%av%cfg_flush_vec = .true.
            CASE('.NORESTART');  config%decomp%cfg_DumpDensRestart = .false.
            CASE('.PAO');  config%decomp%cfg_pao = .true.
            CASE('.HESONLY');    config%opt%cfg_hesonly = .true.
                                 config%decomp%cfg_check_converged_solution = .true.
                                 config%decomp%cfg_hessian_nvec = 1
            CASE('.HESVEC');     READ(LUCMD,*) config%decomp%cfg_hessian_nvec 
                                 config%decomp%cfg_check_converged_solution = .true.
            CASE('.HLMAXIT');    READ(LUCMD,*) config%decomp%cfg_homolumo_maxit
            CASE('.NO HLGAP');   config%diag%nofinalhomolumo = .true.
            CASE('.L2THR');      READ(LUCMD,*) config%opt%cfg_level2_convfactor
            CASE('.LCV');        config%decomp%cfg_lcv = .true.
            CASE('.LCM');        config%decomp%cfg_lcv = .true. ; config%decomp%cfg_lcm=.true.
            CASE('.LCVBF');      config%decomp%cfg_lcv = .true. ; config%decomp%cfg_lcvbf=.true.
            CASE('.LCMBF');      config%decomp%cfg_lcv = .true. ; config%decomp%cfg_lcvbf=.true.; config%decomp%cfg_lcm=.true.
            CASE('.LEVELSH');    ALLOCATE(config%diag%cfg_levelshifts(100)) ; config%diag%cfg_levelshifts = 0.0E0_realk
                                 READ(LUCMD,*) config%diag%cfg_nshifts,(config%diag%cfg_levelshifts(i),i=1,config%diag%cfg_nshifts)
                                 config%diag%cfg_fixed_shift = .true. ; config%diag%cfg_custom_shift = .true.
            CASE('.LINCOMB');    READ(LUCMD,*) config%opt%cfg_weight_param 
!            Both these iterative lowdin decomp are broken
!            CASE('.LWITER') ;    config%decomp%lowdin_diagonalize = .false.; config%decomp%lowdin_iterative  = .true. 
!            CASE('.LWQITER') ;   config%decomp%lowdin_diagonalize = .false.; config%decomp%lowdin_qiterative = .true.
            CASE('.MAXELM');     READ(LUCMD,*) config%solver%cfg_max_element
                                               config%solver%set_max_element = config%solver%cfg_max_element
            CASE('.MAXIT');      READ(LUCMD,*) config%opt%cfg_max_linscf_iterations
            CASE('.NOQUITMAXIT');READ(LUCMD,*) config%opt%cfg_max_linscf_iterations
                                               config%opt%opt_quit=.false.
            CASE('.MAXRATI');    READ(LUCMD,*) maxratio
                                               config%av%cfg_settings%max_dorth_ratio = maxratio
            CASE('.MAXSTEP');    READ(LUCMD,*) config%solver%cfg_max_step 
                                               config%solver%set_max_step = config%solver%cfg_max_step
            CASE('.MICTHRS');    READ(LUCMD,*) config%solver%cfg_micro_thresh
                                 config%decomp%cfg_micro_thresh = config%solver%cfg_micro_thresh
            CASE('.MICROVECS');  READ(LUCMD,*) config%solver%cfg_arh_microvecs
            CASE('.MINDAMP');    READ(LUCMD,*) config%solver%cfg_min_lshift
                                               config%diag%cfg_min_lshift = config%solver%cfg_min_lshift
            CASE('.MOCHANGE');   print*,'Keyword .MOCHANGE nolonger supported'
                                 call lsquit('Keyword .MOCHANGE nolonger supported',-1)
            CASE('.MUOPT');      config%diag%CFG_lshift = diag_lshift_search
                                 config%av%CFG_lshift = Diag_lshift_search
            CASE('.NALPHA');     read(LUCMD,*) config%decomp%nocca ; config%decomp%alpha_specified = .true.
            CASE('.NBETA');      read(LUCMD,*) config%decomp%noccb ; config%decomp%beta_specified = .true.
            CASE('.NOAV');       config%av%CFG_averaging =   config%av%CFG_AVG_none  
            CASE('.NO HLSHIFT'); config%solver%lshift_by_hlgap = .false. !Don't use the default scheme (level shift by homo lumo gap), 
                                                                         !use instead the "old" scheme developed for the Davidson algorithm
            CASE('.NEWDAMP');    config%solver%cfg_arh_newdamp = .true.
            CASE('.NVEC');       READ(LUCMD,*) NVEC; config%av%cfg_settings%max_history_size = NVEC
                                 config%av%diis_history_size = NVEC
!                                 config%av%ediis_history_size = NVEC
            CASE('.NVECDII');    READ(LUCMD,*) NVEC
                                 config%av%diis_history_size = NVEC
!                                 config%av%ediis_history_size = NVEC
            CASE('.NOPREC');     config%solver%cfg_NOPREC = .true.
                                 config%decomp%cfg_NOPREC = .true.
            CASE('.INCREM');     config%opt%cfg_incremental = .true.
            CASE('.NOSHIFT');    config%diag%cfg_lshift = diag_lshift_none
                                 config%av%cfg_lshift = diag_lshift_none
                                 config%solver%cfg_nodamp = .true.
!            CASE('.OAO');        config%decomp%cfg_do_in_oao = .true.; config%opt%cfg_start_guess = 'H1OAO'
            CASE('.OVERLAP');    READ(LUCMD,*) min_density_overlap
                                 config%av%cfg_settings%min_density_overlap &
                                     & = min_density_overlap
            !CASE('.PURIFY');     config%opt%cfg_density_method = config%opt%cfg_f2d_purification - NO LONGER SUPPORTED! /Stinne 16-08-2010
            !                     read(LUCMD,*) config%opt%cfg_purification_method
            CASE('.PRINTFINALCMO'); config%opt%print_final_cmo=.true.
!            CASE('.MATRICESINMEMORY'); config%integral%MATRICESINMEMORY=.true.
            CASE('.SKIPSCFLOOP');    
               config%diag%CFG_restart =  .TRUE.
               config%skipscfloop =  .TRUE.
            CASE('.RESTART');    config%diag%CFG_restart =  .TRUE.
            CASE('.CRASHCALC');    config%opt%crashcalc =  .TRUE.
            CASE('.TESTABSVAL');    config%decomp%debugAbsOverlap=.true.
            CASE('.PURIFYRESTARTDENSITY'); config%diag%CFG_purifyrestart =  .TRUE.
            CASE('.REDO L2');    config%diag%cfg_redo_l2 = .true.
            CASE('.TRANSFORMRESTART');    config%decomp%CFG_transformrestart =  .TRUE. 
            CASE('.RH');         config%opt%CFG_density_method =  config%opt%CFG_F2D_ROOTHAAN
            CASE('.SAFE');       config%av%CFG_safe = .true.
! obsolete keyword - noone knows what it does. Not in manual. Not in testcases
!            CASE('.SCALVIR');    config%opt%cfg_scale_virt = .true.
            !SOEO keywords
            !To-do: Collect in some read-soeo-input
            CASE('.SOEO');       config%soeoinp%cfg_soeo = .true.
            CASE('.TEST SOEO');  config%soeoinp%test = .true.
            CASE('.SOEOSAVE');   config%soeoinp%cfg_save = .true.
            CASE('.SOEORST');    config%soeoinp%cfg_restart = .true.
            CASE('.SOEOPRNT');   read(LUCMD,*) config%soeoinp%prnt
            CASE('.SOEOSPACE');  read(LUCMD,*) config%soeoinp%Nocc, config%soeoinp%Nact
                                 config%soeoinp%spaceinput = .true.
            CASE('.SOEOMATHR');  read(LUCMD,*) config%soeoinp%macrothresh
            CASE('.SOEOMITHR');  read(LUCMD,*) config%soeoinp%microthresh
            CASE('.SOEOTRUST');  read(LUCMD,*) config%soeoinp%trust
            CASE('.SOEOOCC');    config%soeoinp%occsinput = .true.
                                 read(LUCMD,*) config%soeoinp%Nfullocc
                                 read(LUCMD,*) config%soeoinp%Nfracocc
                                 if (config%soeoinp%Nfracocc <= &
                                   & size(config%soeoinp%fracoccs)) then
                                   read (LUCMD,*) config%soeoinp%fracoccs&
                                                   &(1:config%soeoinp%Nfracocc)
                                 else
                                   write (lupri,*) &
                                  &'No more than 100 fractionally occupied&
                                  & orbitals allowed'
                                   write (lupri,*) &
                                  &'Change in linears/soeo-typedef.f90 type&
                                  & soeoItem_input'
                                   call lsquit ('Nfracocc > 100',-1)
                                 endif
            CASE('.SOEOOCCB');   if (config%soeoinp%occsinput) then
                                   read(LUCMD,*) config%soeoinp%Nfulloccb
                                   read(LUCMD,*) config%soeoinp%Nfracoccb
                                   if (config%soeoinp%Nfracoccb <= &
                                     & size(config%soeoinp%fracoccsb)) then
                                     read (LUCMD,*) config%soeoinp%fracoccsb&
                                                     &(1:config%soeoinp%Nfracoccb)
                                   else
                                     write (lupri,*) &
                                    &'No more than 100 fractionally occupied&
                                    & orbitals allowed'
                                     write (lupri,*) &
                                    &'Change in linears/soeo-typedef.f90 type&
                                    & soeoItem_input'
                                     call lsquit ('Nfracoccb > 100',-1)
                                   endif
                                 endif
            CASE('.SOEOGC');     config%soeoinp%cfg_grandcan = .true.
            CASE('.SOEODIPOLE'); config%soeoinp%cfg_dipole = .true.
            !END SOEO keywords
            CASE('.SPARSETEST'); config%sparsetest = .true.
            CASE('.SPIN');       READ(LUCMD,*) config%decomp%spin
            CASE('.STABILITY');  config%decomp%cfg_check_converged_solution = .true.
                                 config%decomp%cfg_hessian_nvec = 1
            CASE('.STAB MAXIT'); READ(LUCMD,*) config%decomp%cfg_check_maxit
            CASE('.START');      READ(LUCMD,*) config%opt%cfg_start_guess 
                                 STARTGUESS = .TRUE.
            CASE('.NOATOMSTART');config%opt%add_atoms_start=.FALSE.
            CASE('.DENSELEVEL2');
               !Use Dense Matrix type in level 2 of Trilevel
               config%opt%DENSELEVEL2=.TRUE.
            CASE('.MWPURIFYATOMSTART');               
               !Perform McWeeny purification on the non idempotent Atoms Density
               config%opt%MWPURIFYATOMSTART=.TRUE.
            CASE('.UNREST');     config%decomp%cfg_unres=.true.
                                 config%integral%unres=.true.
                                 config%diag%cfg_unres=.true.
                                 config%opt%cfg_unres=.true.
                                 config%soeoinp%cfg_unres=.true.
                                 config%response%RSPsolverinput%cfg_unres = .true.
            CASE('.UNSAFE');     config%solver%cfg_arh_crop_safe = .false.
            CASE('.VanLenthe');  config%opt%CFG_density_method =  config%opt%CFG_F2D_ROOTHAAN !Diagonalization
                                 config%av%CFG_averaging = config%av%CFG_AVG_van_lenthe
                                 config%diag%CFG_lshift = diag_lshift_vanlenthe
                                 config%av%CFG_lshift = diag_lshift_vanlenthe
            CASE('.ZERO');       READ(LUCMD,*) zero
                                 call mat_zero_cutoff(zero)
            CASE('.LOW ACCURACY START'); config%integral%LOW_ACCURACY_START=.TRUE.
            CASE DEFAULT
               WRITE(config%LUPRI,*) ' Keyword ',WORD,' not recognized in read_dalton_input'
               CALL lsQUIT('Illegal keyword in *DENSOPT section',config%lupri)
            END SELECT
         ENDIF
         IF(PROMPT .EQ. '$') THEN
            WRITE(LUPRI,*) ' Keyword ',WORD,' not recognized under *DENSOPT in read_dalton_input'
            CALL lsQUIT('Illegal keyword under *DENSOPT in read_dalton_input.',config%lupri)
         ENDIF
         IF(PROMPT .EQ. '**') THEN
            READWORD=.FALSE.
            EXIT
         ENDIF
         IF(PROMPT(1:1) .EQ. '*') THEN
            READWORD=.FALSE.
            EXIT
         ENDIF
      ENDDO
   ENDIF
   IF(WORD(1:6) .EQ. '**INFO') THEN
      READWORD = .TRUE.
      call config_info_input(config,lucmd,readword,word)
   ENDIF

#ifdef MOD_UNRELEASED
   ! Geometrical Hessian input section
   IF (WORD(1:12) == '**GEOHESSIAN') THEN
      READWORD = .TRUE.
      config%geoHessian%do_geoHessian = .TRUE.
      config%geoHessian%testContrib   = .FALSE.
      config%geoHessian%DebugGen1Int  = .FALSE.
      config%geoHessian%IntPrint = 1
      call GEOHESSIAN_INPUT(config%geohessian,readword,word,lucmd,lupri)
   ENDIF
#endif


   ! KK, change from $RESPONS to **RESPONS to be consistent with other input structure.
   ResponseInput: IF (WORD(1:9) == '**RESPONS') THEN
      READWORD=.TRUE.
      call config_rsp_input(config,lucmd,readword,WORD)
   END IF ResponseInput

   ! Input for DEC calculation
   DECInput: IF (WORD(1:5) == '**DEC') THEN
      READWORD=.TRUE.
      config%doDEC = .true.
      call config_dec_input(lucmd,config%lupri,readword,word,.false.,config%doF12,config%doRIMP2)
      config%integral%PreCalcDFscreening =  config%doRIMP2
      config%integral%PreCalcF12screening =  config%doF12
   END IF DECInput

   ! Input for full molecular CC calculation
   ! (uses same setup as DEC calculation)
   CCinput: IF (WORD(1:4) == '**CC') THEN
      READWORD=.TRUE.
      config%doDEC = .true.
      call config_dec_input(lucmd,config%lupri,readword,word,.true.,config%doF12,config%doRIMP2)
      config%integral%PreCalcDFscreening =  config%doRIMP2
      config%integral%PreCalcF12screening =  config%doF12
   END IF CCinput


   ! Input for PLT plotting: Calculate density, orbital, etc. at grid points in space based on
   ! density or orbital matrix file from previous LSDALTON calculation.
   ! This will overrule all other inputs and an LSDALTON calculation as such will not be carried out!
   ! Note: The gridbox information of the plt structure is not set here but in the **GRID section.
   PLTinput: IF (WORD(1:5) == '**PLT') THEN
      READWORD=.TRUE.
      config%doPLT = .true.
      call config_plt_input(lucmd,config%lupri,readword,word,config%plt)
   END IF PLTinput

   ! Input for PLT gridbox: Set gridbox related info for PLT structure 
   ! (nX,nY,nZ,nGRIDPOINTS,deltax,deltay,deltaz,X1,Y1,Z1,gridbox_defined)
   PLTGRIDinput: IF (WORD(1:9) == '**PLTGRID') THEN
      READWORD=.TRUE.
      call config_pltgrid_input(lucmd,config%lupri,readword,word,config%plt)
   END IF PLTGRIDinput


!   
! Find optimization input section
!
   IF (WORD(1:7) .EQ. '**OPTIM') THEN  
      config%optinfo%optimize = .TRUE.
      CALL LS_optimization_input(config%optinfo,readword,word,lucmd, &
           & lupri,config%molecule%nAtoms)
   ENDIF

! Orbital localization section
!
   IF (WORD(1:7) .EQ. '**LOCAL') THEN
     READWORD=.true.
     call orbitalloc_input(lucmd,lupri,config,READWORD,word)
   ENDIF



!
! Find dynamics input section
!
   IF (WORD(1:7) .EQ. '**DYNAM') THEN
      config%dynamics%do_dynamics = .TRUE.
      CALL LS_dynamics_input(config%dynamics,readword,word,&
           & lucmd,lupri,config%molecule%NAtoms)
   ENDIF
!

!
! Find PCM input section
!
#ifdef HAS_PCMSOLVER
   IF (WORD(1:5) .EQ. '**PCM') THEN
      config%pcm%do_pcm = .true.
      call ls_pcm_input(config%pcm,readword,word,&
           & lucmd,lupri)
   ENDIF
#endif
!
   
#ifdef MOD_UNRELEASED
   IF (WORD(1:5) == '**PBC') THEN
     READWORD=.TRUE.
     !should be in MOLECULE.INP not LSDALTON.INP
     !READ(WORD(6:),*) config%latt_config%max_layer,config%latt_config%nneighbour
     config%latt_config%max_layer = 10
     config%latt_config%comp_pbc= .true.
     config%latt_config%wannier_direct= 'indirectly'
     config%latt_config%testcase= .false.
     config%latt_config%compare_elmnts= .false.
     config%latt_config%lmax=15
     config%latt_config%Tlmax=15
     config%latt_config%num_its=100
     config%latt_config%num_store=7
     config%latt_config%error=1.0d-8
     config%latt_config%num_its_densmat=3
     config%latt_config%nf=6
     config%latt_config%ndmat=6
     config%latt_config%intthr=1.0d-8
     config%latt_config%realthr = -12
     config%latt_config%read_file=.false.
     config%latt_config%store_mats=.false.
     DO
     IF(READWORD) THEN
        READ (LUCMD, '(A40)') WORD
        READWORD=.TRUE.
     ENDIF
     PROMPT = WORD(1:2)
     IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
     IF(PROMPT .EQ. '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
     IF(PROMPT(1:1) .EQ. '.') THEN
        SELECT CASE(WORD)
        CASE ('.COMP') 
          config%latt_config%compare_elmnts= .true.
        CASE ('.WANNIER'); config%latt_config%wannier_direct= 'directly'
        CASE('.STARTDENS')
          READ (LUCMD, '(I2)') config%latt_config%num_its_densmat
        CASE('.LATTICE')
           READ (LUCMD, '(I2)')config%latt_config%max_layer
        CASE('.NFIELD')
          READ (LUCMD, '(I2)') config%latt_config%nf
        CASE('.NDENSMATCUTOFF')
          READ (LUCMD, '(I2)') config%latt_config%ndmat  
        CASE('.RECLAT')
          READ (LUCMD, *) config%latt_config%nk1,config%latt_config%nk2,&
                               & config%latt_config%nk3

          if(config%latt_config%nk2 .gt. 1 ) then
            if(.not.config%latt_config%ldef%is_active(2)) then
              WRITE(*,*) 'Reciprocal vector 2 should be set to 1'
              call LSQUIT('ERROR IN RECIPROCAL VECTORS',lupri)
            endif
          endif

          if(config%latt_config%nk3 .gt. 1 ) then
            if(.not.config%latt_config%ldef%is_active(3)) then
              WRITE(*,*) 'Reciprocal vector 3 should be set to 1'
              call LSQUIT('ERROR IN RECIPROCAL VECTORS',lupri)
            endif
          endif

        CASE('.MLMAX')
          READ (LUCMD, '(I2)')config%latt_config%lmax
          config%latt_config%Tlmax = config%latt_config%lmax

        CASE('.TESTCASE')
          config%latt_config%testcase= .true.

        CASE('.DEBUGFILE')
          config%latt_config%read_file=.true.
          READ (LUCMD, '(A30)') config%latt_config%debugdensfile

        CASE ('.WRITE TO FILE') 
          config%latt_config%store_mats= .true.
          
        CASE('.CONVTHR')
          READ(LUCMD,*) config%latt_config%error

        CASE('.INTTHR')
          READ(LUCMD,*) config%latt_config%intthr

        CASE('.SCFCYCLES')
          READ(LUCMD,*) config%latt_config%num_its

        CASE('.PREVCYCLES')
          READ(LUCMD,*) config%latt_config%num_store

        CASE('.REALTHR')
          READ(LUCMD,*) config%latt_config%realthr

        CASE DEFAULT
           WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                & '" not recognized in **PBC readin.'
           print*,'Keyword ',WORD
           CALL lsQUIT('Illegal keyword in **PBC.',lupri)
        END SELECT
     ENDIF
     IF (WORD == '*END OF INPUT') THEN
        READWORD=.FALSE.
        DONE=.TRUE.
        EXIT
     ENDIF
  ENDDO

   ENDIF
#endif

#if defined(ENABLE_QMATRIX)
   ! QMatrix library
   if (WORD=='**QMATRIX') then
       config%do_qmatrix = .true.
       ! initializes the QMatrix interface
       call ls_qmatrix_init(config%ls_qmat)
       ! processes input
       READWORD = .true.
       call ls_qmatrix_input(config%ls_qmat, LUCMD, LUPRI, READWORD, WORD)
   end if
#endif

   IF (WORD == '*END OF INPUT') THEN
      DONE=.TRUE.
   ENDIF
ENDDO
!IF(LSDALTON .AND. (.NOT. STARTGUESS)) config%integral%TRILEVEL = .TRUE.
!IF(.NOT.config%integral%LINSCA)THEN
!   WRITE (LUPRI,'(/,3A,/)') ' The Keyword .LINSCA was not used in&
!   & **INTEGRAL section.'
!   CALL lsQUIT('.LINSCA WAS NOT USED.')
!ENDIF
CALL lsCLOSE(LUCMD,'KEEP')

if (config%solver%do_dft) then
   hfweight = 0.0E0_realk
   CALL II_DFTsetFunc(config%integral%dft%dftfunc,hfweight,lupri)
   !it is assumed that hfweight is set to zero and only  
   !changed if the functional require a HF weight  
   !different from zero. 
   IF (.NOT.exchangescale) THEN
     config%integral%exchangeFactor = hfweight
     config%integral%dft%HFexchangeFac = hfweight
   ENDIF
#ifdef BUILD_CGTODIFF
#ifdef VAR_MPI
   call lsquit('cgto_diff_eri not testet for MPI',-1)
#endif                     
   !call cgto_diff_eri_xfac_general(config%integral%exchangeFactor)
#endif
endif

if(config%solver%do_dft.OR.config%integral%ADMM_EXCHANGE)THEN
   call init_gridObject(config%integral%dft,config%integral%DFT%GridObject)
   call init_dftfunc(config%integral%DFT)
   IF (config%integral%CAM) THEN
      write(CAMalpha,'(G23.16)') config%integral%CAMalpha
      write(CAMbeta,'(G23.16)')  config%integral%CAMbeta
      write(CAMmu,'(G23.16)')    config%integral%CAMmu
      CAMalpha = trim(adjustl(CAMalpha))
      CAMbeta  = trim(adjustl(CAMbeta))
      CAMmu    = trim(adjustl(CAMmu))
      IF (USEXCFUN) THEN
         !IF (config%integral%DFT%DFTfuncObject().EQ.'CAMB3LYP') THEN
         write(Func,'(A8)') 'CAMCOMPX'
         !ELSE
         !  CALL LSQUIT('Todo: fix admm camcompx for camb3lyp with',-1)
         !  !Simen fix admm for non default
         !  write(Func,'(A28,A11,G22.16,A10,G22.16,A13,G22.16)') 'GGAKEY BECKEX=1 BECKECAMX=-1',&
         !               & ' CAM_ALPHA=',config%integral%CAMalpha,&
         !               & ' CAM_BETA=',config%integral%CAMbeta,&
         !               & ' RANGESEP_MU=',config%integral%CAMmu
         !ENDIF
      ELSE
         write(Func,'(A15,A18,A6,A18,A4,A18)') 'Camcompx alpha=',CAMalpha,&
           & ' beta=',CAMbeta,' mu=',CAMmu
      ENDIF
   ELSE
      Func = config%integral%admm_func
   ENDIF
   config%integral%DFT%DFTfuncObject(dftfunc_ADMML2) = Func
endif

END SUBROUTINE read_dalton_input


subroutine read_dft_input(config,lucmd,lupri)

  implicit none
  !> Contains info, settings and data for entire calculation
  type(ConfigItem), intent(inout) :: config
  character(len=80)  :: WORD,MUWORD,TMPWORD,FormatString,CAMalpha,CAMbeta,CAMmu
  character(len=1024):: XCfunString
  integer,intent(in) :: LUCMD !Logical unit number for the daltoninput
  integer,intent(in) :: LUPRI !Logical unit number for the daltonoutput file
  integer            :: ipos,ipos2,iposmu
  
  config%opt%calctype = config%opt%dftcalc !DFT calc
  config%av%CFG_SET_type = config%av%CFG_THR_dft
  config%solver%do_dft = .true.
  config%soeoinp%do_dft = .true.
  config%davidSCF%arh_dodft = .true.
  
  DO 
     READ (LUCMD, '(A80)') WORD
     IF ((WORD(1:1) .EQ. '!') .OR. (WORD(1:1) .EQ. '#')) CYCLE   
     IF (WORD(1:1) .EQ. '.' .OR. WORD(1:1) .EQ. '*') THEN
        WRITE (LUPRI,'(/A/A//A)')&
             & '--> Input error for line following .DFT',&
             & '    expected functional specification but read:',&
             & WORD
     ELSE
        call capitalize_string(WORD)
        IF ((INDEX(WORD,'LDAERF')).NE.0) THEN !LDAERF   
           config%integral%CAM=.TRUE. !Activate lr intergrals
           config%integral%CAMalpha=0E0_realk !No exact exchange
           config%integral%CAMbeta=1E0_realk !Lr HF exchange prefactor
           READ (LUCMD, '(A80)') MUWORD
           call capitalize_string(MUWORD)
           IPOS = INDEX(MUWORD,'MU')
           mu: IF (IPOS .NE. 0) THEN
              IPOS2 = INDEX(MUWORD(IPOS:),'=')
              IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 3)) THEN
                 WRITE (LUPRI,'(2X,A40)') 'Incorrect input for ERF parameters'
                 WRITE (LUPRI,'(2X,A40)') 'Format is "mu=?"'
                 CALL lsQUIT('Incorrect range-separation parameter mu',lupri)
              ELSE
                 !Set value of the range-separation parameter mu
                 READ (MUWORD((IPOS+IPOS2):80),*) config%integral%CAMmu 
                 write(LUPRI,*) ''
                 WRITE(LUPRI,*) 'This is a LDAERF functional with mu=',config%integral%CAMmu
                 !String to be passed to xcfun_host
                 !EXX=1 sets hfweight and activates exchange integral evaluation
                 XCfunString = ('LDAERF EXX=1 RANGESEP_MU='//trim(MUWORD((IPOS+IPOS2):80))) 
                 write(LUPRI,*) XCfunString
              ENDIF
           ELSE
              CALL lsQUIT('Range-separation parameter mu required when using LDAERF',lupri)
           ENDIF mu
        ELSE IF ((INDEX(WORD,'CAM')).NE.0) THEN !CAM
           config%integral%CAM=.TRUE.
           IPOS = INDEX(WORD,'ALPHA')
           IF (IPOS .NE. 0) THEN
              IPOS2 = INDEX(WORD(IPOS:),'=')
              IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
                 WRITE (LUPRI,'(2X,A40)') 'Incorrect input for CAM parameters'
                 WRITE (LUPRI,'(2X,A40)') 'Format is "alpha=?  beta=? mu=?"'
                 CALL lsQUIT('Incorrect input for alpha parameter',lupri)
              ELSE
                 READ (WORD((IPOS+IPOS2):80),*) config%integral%CAMalpha
                 IPOS = INDEX(WORD,'BETA')
                 IF (IPOS .NE. 0) THEN
                    IPOS2 = INDEX(WORD(IPOS:),'=')
                    IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 5)) THEN
                       WRITE (LUPRI,'(2X,A40)') 'Incorrect input for CAM parameters'
                       WRITE (LUPRI,'(2X,A40)') 'Format is "alpha=?  beta=? mu=?"'
                       CALL lsQUIT('Incorrect input for beta parameter',lupri)
                    ELSE
                       READ (WORD((IPOS+IPOS2):80),*) config%integral%CAMbeta
                       IPOS = INDEX(WORD,'MU')
                       IF (IPOS .NE. 0) THEN
                          IPOS2 = INDEX(WORD(IPOS:),'=')
                          IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 3)) THEN
                             WRITE (LUPRI,'(2X,A40)') 'Incorrect input for CAM parameters'
                             WRITE (LUPRI,'(2X,A40)') 'Format is "alpha=?  beta=? mu=?"'
                             CALL lsQUIT('Incorrect input for mu parameter',lupri)
                          ELSE
                             READ (WORD((IPOS+IPOS2):80),*) config%integral%CAMmu
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
              if (USEXCFUN) then
                 write(CAMalpha,'(G23.16)') config%integral%CAMalpha
                 write(CAMbeta,'(G23.16)')  config%integral%CAMbeta
                 write(CAMmu,'(G23.16)')    config%integral%CAMmu
                 CAMalpha = trim(adjustl(CAMalpha))
                 CAMbeta  = trim(adjustl(CAMbeta))
                 CAMmu    = trim(adjustl(CAMmu))
                 FormatString= "(A18,A11,A18,A10,A18,A29,A13,A18)" 
                 write(XCfunString,FormatString) 'GGAKEY BECKECAMX=1',&
                      & ' CAM_ALPHA=',CAMalpha,&
                      & ' CAM_BETA=',CAMbeta,&
                      & ' VWN5C=0.19 LYPC=0.81 EXX=1.0',&
                      & ' RANGESEP_MU=',CAMmu
                 XCfunString = trim(XCfunString)
              else
                 XCfunString = WORD
              endif
           ELSE
              config%integral%CAMalpha=0.19E0_realk
              config%integral%CAMbeta=0.46E0_realk
              config%integral%CAMmu=0.33E0_realk
              XCfunString = WORD
           ENDIF
           WRITE(LUPRI,*) 'This is a CAM functional with'
           WRITE(LUPRI,*) 'config%integral%CAMalpha',config%integral%CAMalpha
           WRITE(LUPRI,*) 'config%integral%CAMbeta',config%integral%CAMbeta
           WRITE(LUPRI,*) 'config%integral%CAMmu',config%integral%CAMmu
        ELSE
           XCfunString = WORD
        END IF
        ! Remove starting blancks
        DO IPOS=1,LEN(XCfunString)
           IF (XCfunString(IPOS:IPOS).NE." ") EXIT
        ENDDO
        ! Store functional string
        ! Waiting for user to decide whether to use xcfun or not
        config%integral%dft%dftfunc = XCfunString(IPOS:)
     ENDIF
     EXIT
  ENDDO
  
end subroutine read_dft_input

!> \brief For DEC calculations, check that DEC input is consistent with input for other parts 
!> of the code. Specifically, check input for geometry optimization orbital localization.
!> If necessary, modify config structure to comply with DEC calculation.
!> \author Kasper Kristensen
!> \date April 2013
subroutine DEC_meaningful_input(config)
  implicit none
  !> Contains info, settings and data for entire calculation
  type(ConfigItem), intent(inout) :: config
  logical :: NotAcceptedModel

  ! Only make modifications to config for DEC calculation AND if it is not
  ! a full CC calculation
  DECcalculation: if(config%doDEC) then

     if(DECinfo%ccmodel == MODEL_CCSDpT) then
        if (.not. DECinfo%full_molecular_cc) DECinfo%print_frags = .true.
     endif

     ! CCSD does not work for SCALAPACK, Hubi please fix!
     if(matrix_type==mtype_scalapack .and. (DECinfo%ccmodel/=MODEL_MP2) ) then
        call lsquit('Error in input: Coupled-cluster beyond MP2 is not implemented for .SCALAPACK!',-1)
     end if
     ! CCSD does not work for CSR, Thomas/Hubi please fix - make not matrix type
     ! subroutines
     if(config%opt%cfg_prefer_CSR .and. (DECinfo%ccmodel/=MODEL_MP2) ) then
        call lsquit('Error in input: Coupled-cluster beyond MP2 is not implemented for .CSR!',-1)
     end if
     if(DECinfo%ccmodel.EQ.MODEL_RIMP2 .AND. (.NOT. config%integral%basis(AuxBasParam)))then
        WRITE(config%LUPRI,'(/A)') &
             &     'You have specified .RIMP2 in the input but not supplied a fitting basis set'
        CALL lsquit('MP2 RI input inconsitensy: add fitting basis set',config%lupri)
     endif
     ! DEC and response do not go together right now...
     if(config%response%tasks%doResponse) then
        call lsquit('Error in input: **DEC or **CC cannot be used together with **RESPONS!',-1)
     end if

     if(config%opt%calctype == config%opt%dftcalc) then
        DECinfo%DFTreference = .TRUE.
        WRITE(config%lupri,*)' '
        WRITE(config%lupri,*)'Warning you use a Kohn-Sham Reference for DEC or CC calculation!'
        WRITE(config%lupri,*)'I hope you know what you are doing!'
        WRITE(config%lupri,*)' '
     end if

     DECinfo%GCBASIS = config%decomp%cfg_gcbasis

     ! DEC geometry optimization 
     ! *************************
     GeoOptCheck: if(config%optinfo%optimize .or. config%dynamics%do_dynamics) then
        ! Always use dynamical optimization procedure
        config%optinfo%dynopt=.true.

        ! Modify DECinfo to calculate first order properties (gradient) for MP2
        DECinfo%gradient=.true.
        DECinfo%first_order=.true.
        if (DECinfo%ccmodel /= MODEL_MP2) then
           write(DECinfo%output,*) "WARNING: DEC Geometry optimization only available for MP2"
           write(DECinfo%output,*) "WARNING: We are switching to DEC-MP2  !!!"
           DECinfo%ccmodel = MODEL_MP2
        end if 

        ! DEC restart for geometry optimizations not implemented
        if(DECinfo%HFrestart .or. DECinfo%DECrestart) then
           write(config%lupri,*) 'Warning: DEC restart not implemented for geometry optimization...'
           write(config%lupri,*) '--> I turn off DEC restart!'
           DECinfo%HFrestart=.false.
           DECinfo%DECrestart=.false.
        end if
     end if GeoOptCheck

     ! Orbital localization must be turned on for DEC calculations
     ! ***********************************************************
     ! If LCM and orbital localization has not been set in input, we enforce it here.
     ! We use the orbital variance localization function with exponent m=2 for
     ! both occupied and virtual localization.
     OrbLocCheck: if( (.not. config%decomp%cfg_mlo) .or. (.not. config%decomp%cfg_lcv) &
          &  .or. (.not. config%decomp%cfg_lcm) ) then
        !Either .LCM or .PSM/.PFM/... keyword was not used 
        IF(.NOT.DECinfo%use_canonical)THEN
           ! Turn on LCM scheme
           IF(.NOT.config%decomp%cfg_lcv)THEN
              WRITE(config%lupri,*)'Warning **DEC and **CC enforces the use of the .LCM keyword'
              config%decomp%cfg_lcv = .true.
              config%decomp%cfg_lcm=.true.
           ENDIF
           IF(.NOT.config%decomp%cfg_mlo)THEN
              ! Turn on orbital localization
              config%decomp%cfg_mlo = .true.
           
              ! use orbspread localization function and line search
              config%davidOrbLoc%orbspread=.true.
              config%davidOrbLoc%linesearch=.true.

              ! Exponent 2 for both occ and virt orbitals
              config%decomp%cfg_mlo_m(1) = 2
              config%decomp%cfg_mlo_m(2) = 2   
           ELSE
              !orbital localization already turned on by keywords.
              !Use user specified input!
           ENDIF
        else
           !No reason to Localize when Using canonical Orbitals
        endif
        ! For the release we only include DEC-MP2
#ifndef MOD_UNRELEASED
        NotAcceptedModel = DECinfo%ccmodel/=MODEL_MP2 .AND. DECinfo%ccmodel/=MODEL_RIMP2
        if(NotAcceptedModel .and. (.not. DECinfo%full_molecular_cc) ) then
           print *, 'Note that you may run a full molecular CC calculation (not linear-scaling)'
           print *, 'using the **CC section rather than the **DEC section.'
           call lsquit('DEC is currently only available for the MP2 and RI-MP2 model!',-1)
        end if
#endif
     end if OrbLocCheck

     ! Check DEC input internally
     call check_dec_input()

  end if DECcalculation


end subroutine DEC_meaningful_input

subroutine PROFILE_INPUT(profinput,readword,word,lucmd,lupri)
  implicit none
  LOGICAL,intent(inout)                :: READWORD
  TYPE(profileinput),intent(inout)   :: profinput
  character(len=80)  :: WORD
  INTEGER,intent(in) :: LUCMD !Logical unit number for the daltoninput
  INTEGER,intent(in) :: LUPRI !Logical unit number for the daltonoutput file
!
  INTEGER            :: IDUMMY,I
  character(len=2)   :: PROMPT
  profinput%doProf = .TRUE.
  DO   
     IF(READWORD) THEN
        READ (LUCMD, '(A40)') WORD
        READWORD=.TRUE.
     ENDIF
     PROMPT = WORD(1:2)
     IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
     IF(PROMPT .EQ. '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
     IF(PROMPT(1:1) .EQ. '.') THEN
        SELECT CASE(WORD) 
        CASE ('.OVERLAP'); PROFINPUT%OVERLAP = .TRUE.
        CASE ('.COULOMB'); PROFINPUT%COULOMB = .TRUE.
        CASE ('.COULOMBECONT'); PROFINPUT%COULOMBECONT = .TRUE.
        CASE ('.EXCHANGE'); PROFINPUT%EXCHANGE = .TRUE.
        CASE ('.EXCHANGEECONT'); PROFINPUT%EXCHANGEECONT = .TRUE.
        CASE ('.EXCHANGEGRAD'); PROFINPUT%EXCHANGEGRAD = .TRUE.
        CASE ('.EXCHANGE MANYD'); PROFINPUT%EXCHANGEmanyD = .TRUE.
        CASE ('.XC'); PROFINPUT%XC = .TRUE.
        CASE ('.XC LINRSP'); PROFINPUT%XCLINRSP = .TRUE.
        CASE ('.XC ENERGY'); PROFINPUT%XCENERGY = .TRUE.
        CASE ('.FOCK');  PROFINPUT%FOCK = .TRUE.
        CASE ('.NEGRAD'); PROFINPUT%NEGRAD = .TRUE.
        CASE ('.ICHORDEC'); PROFINPUT%IchorDEC = .TRUE.
        CASE ('.ICHOR'); PROFINPUT%Ichor = .TRUE.
           PROFINPUT%IchorProfDoIchor = .TRUE.
        CASE ('.ICHOR PURE THERMITE'); PROFINPUT%IchorProfDoIchor = .FALSE.;
           PROFINPUT%IchorProfDoThermite = .TRUE.;
        CASE ('.ICHOR ADD THERMITE'); PROFINPUT%IchorProfDoThermite = .TRUE.;
        CASE ('.ICHOR BASIS'); PROFINPUT%IchorProfInputBasis = .TRUE.
           DO I=1,4
              READ(LUCMD, '(A40)')WORD
              profinput%IchorProfInputBasisString(I) = WORD(1:20)
           ENDDO
        CASE DEFAULT
           WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                & '" not recognized in **PROFILE readin.'
           print*,'Keyword ',WORD
           CALL lsQUIT('Illegal keyword in **PROFILE.',lupri)
        END SELECT
     ENDIF
     IF (WORD(1:2) == '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
  ENDDO
END subroutine PROFILE_INPUT

subroutine GENERAL_INPUT(config,readword,word,lucmd,lupri)
  implicit none
  LOGICAL,intent(inout)                :: READWORD
  type(ConfigItem), intent(inout) :: config
  character(len=80)  :: WORD
  INTEGER,intent(in) :: LUCMD !Logical unit number for the daltoninput
  INTEGER,intent(in) :: LUPRI !Logical unit number for the daltonoutput file
!
  INTEGER            :: IDUMMY
  character(len=2)   :: PROMPT
  DO   
     IF(READWORD) THEN
        READ (LUCMD, '(A40)') WORD
        READWORD=.TRUE.
     ENDIF
     PROMPT = WORD(1:2)
     IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
     IF(PROMPT .EQ. '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
     IF(PROMPT(1:1) .EQ. '.') THEN
        SELECT CASE(WORD) 
        CASE ('.XCFUN')
#ifdef VAR_XCFUN
           USEXCFUN = .TRUE. 
           CONFIG%INTEGRAL%DFT%XCFUN = .TRUE.
           print*,'USEXCFUN',USEXCFUN
           write(lupri,*) 'The XCfun module is activated'
#else
           call lsquit('.XCFUN requires ENABLE_XCFUN', -1)
#endif
        CASE ('.PSFUN')
           USEXCFUN = .FALSE. 
           CONFIG%INTEGRAL%DFT%XCFUN = .FALSE.
           print*,'USEXCFUN',USEXCFUN
           write(lupri,*) 'The XCfun module is deactivated, PS functional evaluation is used instead'
        CASE('.INTERACTIONENERGY')
           !Calculated the Interaction energy 
           !using Counter Poise Correction
           config%InteractionEnergy = .true.
        CASE('.PAPITEST')
           config%papitest=.true.
        CASE('.ACCESS_STREAM')
           ! Use stream access on all files open with lsopen
           config%access_stream = .true.
           access_stream = .true.
        CASE('.SAMESUBSYSTEMS')
           config%SameSubSystems = .true.
        CASE('.SUBSYSTEMDENSITY')
           config%SubSystemDensity = .true.
        CASE('.CSR');        config%opt%cfg_prefer_CSR = .true.
        CASE('.SCALAPACK');  config%opt%cfg_prefer_SCALAPACK = .true.
        CASE('.PDMM')
#ifdef VAR_MPI
           config%opt%cfg_prefer_PDMM = .true.
           !Set tensor synchronization to always, TODO: see if this can be optimized
           call tensor_set_always_sync_true(.true.)
           !Set the background buffer on, this will use additional memory
           call lspdm_init_global_buffer(.true.)
#endif
        CASE('.PDMMBLOCKSIZE');  
#ifdef VAR_MPI
           print*,'PDMMBLOCKSIZE CHOSEN'
           READ(LUCMD,*) infpar%inputBLOCKSIZE
#endif
        CASE('.PDMMGROUPSIZE');
#ifdef VAR_MPI
           READ(LUCMD,*) infpar%PDMMGroupSize
#endif
#ifdef VAR_MPI
        CASE('.SCALAPACKGROUPSIZE');
           READ(LUCMD,*) infpar%ScalapackGroupSize
        CASE('.SCALAPACKWORKAROUND');
           infpar%ScalapackWORKAROUND=.TRUE.
        CASE('.SCALAPACKAUTOGROUPSIZE');
           infpar%ScalapackGroupSize = -1
        CASE('.SCALAPACKBLOCKSIZE');  
           print*,'SCALAPACKBLOCKSIZE CHOSEN'
           READ(LUCMD,*) infpar%inputBLOCKSIZE
#endif
        CASE('.TIME');                  call SET_LSTIME_PRINT(.TRUE.)
        CASE('.GCBASIS');               config%decomp%cfg_gcbasis    = .true. ! left for backward compatibility
        CASE('.NOGCBASIS');             config%decomp%cfg_gcbasis    = .false.
        CASE('.FORCEGCBASIS');          config%INTEGRAL%FORCEGCBASIS = .true.
        CASE('.TESTMPICOPY');           config%doTestMPIcopy         = .true.
        CASE('.TESTHODI');              config%doTestHodi            = .true.
           ! Max memory available on gpu measured in GB. By default set to 2 GB
        CASE('.GPUMAXMEM');             
           READ(LUCMD,*) config%GPUMAXMEM
           IF(config%GPUMAXMEM.LT.0.0E0_realk)THEN
              CALL LSQUIT('.GPUMAXMEM error: less than 0 GB supplied on input',-1)
           ELSEIF(config%GPUMAXMEM.GT.100.0E0_realk)THEN
              CALL LSQUIT('.GPUMAXMEM error: More than 100 GB supplied on input',-1)
           ENDIF
#ifdef VAR_MPI
           call ls_mpibcast(SET_GPUMAXMEM,infpar%master,MPI_COMM_LSDALTON)
           call ls_mpibcast(config%GPUMAXMEM,infpar%master,MPI_COMM_LSDALTON)
#endif
#ifdef VAR_MPI
        CASE('.MAX_MPI_MSG_SIZE_NEL');
           READ(LUCMD,*) SPLIT_MPI_MSG 
           call ls_mpibcast(SET_SPLIT_MPI_MSG,infpar%master,MPI_COMM_LSDALTON)
           call ls_mpibcast(SPLIT_MPI_MSG,infpar%master,MPI_COMM_LSDALTON)
        CASE('.MAX_MPI_MSG_SIZE_ONESIDED_NEL');  
           READ(LUCMD,*) MAX_SIZE_ONE_SIDED
           call ls_mpibcast(SET_MAX_SIZE_ONE_SIDED,infpar%master,MPI_COMM_LSDALTON)
           call ls_mpibcast(MAX_SIZE_ONE_SIDED,infpar%master,MPI_COMM_LSDALTON)
#endif
        CASE DEFAULT
           WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                & '" not recognized in **GENERAL readin.'
           print*,'Keyword ',WORD
           CALL lsQUIT('Illegal keyword in **GENERAL.',lupri)
        END SELECT
     ENDIF

     if (WORD(1:7) == '*TENSOR') then
        READWORD=.TRUE.
        call TENSOR_INPUT(word,LUCMD)
     endif

     IF (WORD(1:2) == '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
     IF (WORD(1:7) == '*END OF') THEN
           READWORD=.FALSE.
           EXIT
     ENDIF

  ENDDO
END subroutine GENERAL_INPUT

subroutine TENSOR_INPUT(word,lucmd)
   implicit none
   character(len=80),intent(inout)  :: word
   integer,intent(in) :: lucmd 
   integer(kind=8)    :: seg_len
   do
      read(LUCMD,'(A40)') word
      if(word(1:1) == '!' .or. word(1:1) == '#') cycle
      if(word(1:1) == '*') then ! New property or *END OF INPUT
         backspace(LUCMD)
         exit
      end if
      select case(word)
      case('.DIL_BACKEND')
         call tensor_set_dil_backend_true(.true.)
      case('.DEBUG')
         call tensor_set_debug_mode_true(.true.)
      case('.SEGMENT_LENGTH')
         read(LUCMD,*) seg_len
         call tensor_set_global_segment_length(seg_len)
      case default
         print *,"UNRECOGNIZED KEYWORD: ",word
         call lsquit("ERROR(GENERAL_INPUT): unrecognized keyword in *TENSOR section",-1)

      end select
   enddo
end subroutine TENSOR_INPUT

subroutine INTEGRAL_INPUT(integral,readword,word,lucmd,lupri)
  implicit none
  LOGICAL,intent(inout)                :: READWORD
  TYPE(integralconfig),intent(inout)   :: integral
  character(len=80)  :: WORD
  INTEGER,intent(in) :: LUCMD !Logical unit number for the daltoninput
  INTEGER,intent(in) :: LUPRI !Logical unit number for the daltonoutput file
!
  INTEGER            :: IDUMMY
  character(len=2)   :: PROMPT
  integral%LINSCA = .TRUE. !should be obsolete

  DO   
     IF(READWORD) THEN
        READ (LUCMD, '(A40)') WORD
        READWORD=.TRUE.
     ENDIF
     PROMPT = WORD(1:2)
     IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
     IF(PROMPT .EQ. '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
     IF(PROMPT(1:1) .EQ. '.') THEN
        SELECT CASE(WORD) 
        CASE ('.CONTANG'); INTEGRAL%CONTANG=.TRUE. ! Specifies that the AO-shell ordering is contracted first then 
                                                   ! angular components (for genereally contracted functions)
        CASE ('.NOGCINTEGRALTRANSFORM'); INTEGRAL%NOGCINTEGRALTRANSFORM=.TRUE.
        CASE ('.NOBQBQ'); INTEGRAL%NOBQBQ=.TRUE.
        CASE ('.FRAGMENT'); READ(LUCMD,*) INTEGRAL%numAtomsPerFragment; INTEGRAL%FRAGMENT = .TRUE.
        CASE ('.DUMP4CENTERERI');  INTEGRAL%DUMP4CENTERERI = .TRUE.
        CASE ('.2CENTERERI'); INTEGRAL%DO2CENTERERI = .TRUE.
        CASE ('.3CENTEROVL'); INTEGRAL%DO3CENTEROVL = .TRUE.
        CASE ('.4CENTERERI');  INTEGRAL%DO4CENTERERI = .TRUE.
        CASE ('.AOPRINT');  READ(LUCMD,*) INTEGRAL%AOPRINT
        CASE ('.BASPRINT');  READ(LUCMD,*) INTEGRAL%BASPRINT
        CASE ('.ForceRIMP2memReduced'); INTEGRAL%ForceRIMP2memReduced = .TRUE.
        CASE ('.DEBUGICHOR')
           INTEGRAL%DEBUGICHOR = .TRUE.
           READ(LUCMD,*) INTEGRAL%DEBUGICHORoption
        CASE ('.DEBUGICHORLINK')
           INTEGRAL%DEBUGICHORLINK = .TRUE.
        CASE ('.DEBUGICHORLINKFULL')
           INTEGRAL%DEBUGICHORLINKFULL = .TRUE.
        CASE ('.DEBUGPROP');  INTEGRAL%DEBUGPROP = .TRUE.
        CASE ('.DEBUGGEN1INT')
#ifdef BUILD_GEN1INT_LSDALTON
           INTEGRAL%DEBUGGEN1INT = .TRUE.
#else
           call lsquit('.DEBUGGEN1INT requires OpenRSP -DBUILD_GEN1INT_LSDALTON',-1)
#endif
        CASE ('.DEBUGCGTODIFF')
#ifdef BUILD_CGTODIFF
           INTEGRAL%DEBUGCGTODIFF = .TRUE.
#else
           call lsquit('.DEBUGCGTODIFF requires OpenRSP -DBUILD_CGTODIFF',-1)
#endif
        CASE ('.DEBUG PS');  INTEGRAL%PS_DEBUG = .TRUE.
        CASE ('.DEBUGEP'); INTEGRAL%DEBUGEP = .TRUE.
        CASE ('.DEBUGSCREEN'); INTEGRAL%DEBUGscreen = .TRUE.
        CASE ('.DEBUGMAGDERIV'); INTEGRAL%DEBUGMAGDERIV = .TRUE.
        CASE ('.DEBUGGEODERIVOVERLAP'); INTEGRAL%DEBUGGEODERIVOVERLAP = .TRUE.
        CASE ('.DEBUGGEODERIVKINETIC'); INTEGRAL%DEBUGGEODERIVKINETIC = .TRUE.
        CASE ('.DEBUGGEODERIVEXCHANGE'); INTEGRAL%DEBUGGEODERIVEXCHANGE = .TRUE.
        CASE ('.DEBUGGEODERIVCOULOMB'); INTEGRAL%DEBUGGEODERIVCOULOMB = .TRUE.
        CASE ('.DEBUGMAGDERIVOVERLAP'); INTEGRAL%DEBUGMAGDERIVOVERLAP = .TRUE.
        CASE ('.DEBUGNUCPOT'); INTEGRAL%DEBUGNUCPOT = .TRUE.
        CASE ('.DEBUGGGEM'); INTEGRAL%DEBUGGGEM = .TRUE.
        CASE ('.DEBUGLSLIB'); INTEGRAL%DEBUGLSlib = .TRUE.
        CASE ('.DEBUGUNCONTAOBATCH'); INTEGRAL%DEBUGUNCONTAOBATCH = .TRUE.
        CASE ('.DEBUGKINETIC'); INTEGRAL%DEBUGKINETIC = .TRUE.
        CASE ('.DEBUGOVERLAP');  INTEGRAL%DEBUGOVERLAP = .TRUE.
        CASE ('.DEBUG4CENTER');  INTEGRAL%DEBUG4CENTER = .TRUE.
        CASE ('.DEBUG4CENTERERI');  INTEGRAL%DEBUG4CENTER_ERI = .TRUE.
        CASE ('.DEBUGCCFRAGMENT');  INTEGRAL%DEBUGCCFRAGMENT = .TRUE.
        CASE ('.DEBUGDECPACKED'); INTEGRAL%DEBUGDECPACKED = .TRUE.
        CASE ('.CART-E'); INTEGRAL%HermiteEcoeff = .FALSE.
        CASE ('.INTPRINT');  READ(LUCMD,*) INTEGRAL%INTPRINT
        CASE ('.NOJENGINE'); INTEGRAL%JENGINE = .FALSE.
        CASE ('.MAXPASSES'); READ(LUCMD,*) INTEGRAL%maxpasses
        CASE ('.NOLINK'); INTEGRAL%LINK = .FALSE.

        !Line search density accelerated screening
        CASE ('.NOLSDASCREEN'); 
           !no line search density accelerated screening
           INTEGRAL%LSDASCREEN = .FALSE.           
           !default is to do line search density accelerated screening
        CASE ('.NOLSDALINK'); INTEGRAL%LSDALINK = .FALSE.
        CASE ('.NOLSDAJENGINE'); INTEGRAL%LSDAJENGINE = .FALSE.
        CASE ('.LSDACOULOMB'); INTEGRAL%LSDACOULOMB = .TRUE.
        CASE ('.LSDASCREEN_THRLOG'); 
           READ(LUCMD,*) INTEGRAL%LSDASCREEN_THRLOG

        !Density accelerated screening   
        CASE ('.DALINK'); INTEGRAL%LINK = .TRUE.
           INTEGRAL%DALINK = .TRUE.
        CASE ('.DAJENGINE'); INTEGRAL%JENGINE = .TRUE.
           INTEGRAL%DAJENGINE = .TRUE.
        CASE ('.DACOULOMB'); INTEGRAL%JENGINE = .FALSE.
           INTEGRAL%DACOULOMB = .TRUE.
        CASE ('.DASCREEN_THRLOG'); 
           READ(LUCMD,*) INTEGRAL%DASCREEN_THRLOG

        CASE ('.NSETUV'); INTEGRAL%nonSphericalETUV = .TRUE.
        CASE('.LOW RJ000 ACCURACY'); INTEGRAL%HIGH_RJ000_ACCURACY = .FALSE.
        CASE('.FTUVMAXPRIM'); READ(LUCMD,*)INTEGRAL%FTUVmaxprim
        CASE ('.LINSCAPRINT')
           READ(LUCMD,*) INTEGRAL%LINSCAPRINT
           INTEGRAL%MOLPRINT=INTEGRAL%LINSCAPRINT
           INTEGRAL%BASPRINT=INTEGRAL%LINSCAPRINT
           INTEGRAL%AOPRINT=INTEGRAL%LINSCAPRINT
!           INTEGRAL%INTPRINT=INTEGRAL%LINSCAPRINT
        CASE ('.PRINTATOMCOORD');  INTEGRAL%PRINTATOMCOORD = .TRUE.
        CASE ('.MOLPRINT');  READ(LUCMD,*) INTEGRAL%MOLPRINT
        CASE ('.NO SCREEN');  
           INTEGRAL%CS_SCREEN = .FALSE. 
           INTEGRAL%PS_SCREEN = .FALSE. 
           INTEGRAL%OD_SCREEN = .FALSE.
           INTEGRAL%MBIE_SCREEN = .FALSE.
           INTEGRAL%PARI_SCREEN = .FALSE. 
        CASE ('.THRESH')
           READ(LUCMD,*) INTEGRAL%THRESHOLD
           IF(INTEGRAL%THRESHOLD.LT.shortintCRIT)THEN
              write(lupri,'(A)')'Error: you cannot chose integral threshold less then'
              write(lupri,'(ES15.6,A)')shortintCRIT, 'due to technical reasons'
              write(lupri,'(A)')'you could use .NO SCREEN (you may have to deactivate LinK with .NOLINK)'
              write(*,'(A)')'Error: you cannot chose integral threshold less then'
              write(*,'(ES15.6,A)')shortintCRIT, 'due to technical reasons'
              write(*,'(A)')'you could use .NO SCREEN (you may have to deactivate LinK with .NOLINK)'
              call lsquit('Error in choice of integral threshold',-1)
           ENDIF
        CASE ('.DO NOT SAVE GAB');  
           INTEGRAL%saveGABtoMem = .FALSE. 
        CASE ('.NO OMP');  INTEGRAL%noOMP = .TRUE. 
        CASE ('.ICHORGPU');  INTEGRAL%IchorForceGPU = .TRUE. 
        CASE ('.ICHORCPU');  INTEGRAL%IchorForceCPU = .TRUE. 
        CASE ('.NO PASS');  INTEGRAL%DOPASS = .FALSE. 
        CASE ('.NO CS');  INTEGRAL%CS_SCREEN = .FALSE. 
        CASE ('.NO PS');  INTEGRAL%PS_SCREEN = .FALSE. 
        CASE ('.NO PARISCREEN');  INTEGRAL%PARI_SCREEN = .FALSE. 
        CASE ('.MBIE');  INTEGRAL%MBIE_SCREEN = .TRUE. 
        CASE ('.ADMM'); !Defaults to ADMM2 with B88 exchange
           IF (INTEGRAL%ADMM_EXCHANGE) THEN
             CALL LSQUIT('Illegal input under **INTEGRAL. Only one choice of ADMM basis.',lupri)
           ENDIF
           INTEGRAL%ADMM_EXCHANGE   = .TRUE.
        CASE ('.ADMM1');
           IF (INTEGRAL%ADMM_EXCHANGE) THEN
             CALL LSQUIT('Illegal input under **INTEGRAL. Only one choice of ADMM basis.',lupri)
           ENDIF
           INTEGRAL%ADMM_EXCHANGE   = .TRUE.
           INTEGRAL%ADMM1           = .TRUE.
        CASE ('.ADMM2');
           IF (INTEGRAL%ADMM_EXCHANGE) THEN
             CALL LSQUIT('Illegal input under **INTEGRAL. Only one choice of ADMM basis.',lupri)
           ENDIF
           INTEGRAL%ADMM_EXCHANGE   = .TRUE.
        CASE ('.ADMMS');
           IF (INTEGRAL%ADMM_EXCHANGE) THEN
             CALL LSQUIT('Illegal input under **INTEGRAL. Only one choice of ADMM basis.',lupri)
           ENDIF
           INTEGRAL%ADMM_EXCHANGE   = .TRUE.
           INTEGRAL%ADMMS           = .TRUE.
        CASE ('.ADMMP');
           IF (INTEGRAL%ADMM_EXCHANGE) THEN
             CALL LSQUIT('Illegal input under **INTEGRAL. Only one choice of ADMM basis.',lupri)
           ENDIF
           INTEGRAL%ADMM_EXCHANGE   = .TRUE.
           INTEGRAL%ADMMP           = .TRUE.
        CASE ('.ADMMQ');
           IF (INTEGRAL%ADMM_EXCHANGE) THEN
             CALL LSQUIT('Illegal input under **INTEGRAL. Only one choice of ADMM basis.',lupri)
           ENDIF
           INTEGRAL%ADMM_EXCHANGE   = .TRUE.
           INTEGRAL%ADMMQ           = .TRUE.
        CASE ('.ADMM-FUNC');
           READ(LUCMD,*) INTEGRAL%ADMM_FUNC
        CASE ('.ADMM-separateX'); !EXPERIMENTAL
                                  !Calculates X and XC independently, default is to calculate
                                  !X+XC with one call to dft (and x with a separate call)
           INTEGRAL%ADMM_separateX  = .TRUE.
        CASE ('.ADMM-BASIS-FILE'); ! EXPERIMENTAL
           !Write file ADMMmin.dat containing info used for constructing a basis set optimized for ADMM
           !We use PRINT_EK3 to calculate the full exchange matrix.
           !We then use the full exchange matrix to converge the standard calculation but for every 
           !iteration we print ADMM stuff:
           !Exchange energy of aux density
           !Exchange energy of full density 
           !(GGA-type) Exchange functional energy of aux density
           !(GGA-type) Exchange functional energy of full density 
           !to easily obtain this infomation we use ADMM_separateX
           IF (.NOT.INTEGRAL%ADMM_EXCHANGE) THEN
              CALL LSQUIT('Illegal input under **INTEGRAL. .ADMM-BASIS-FILE require ADMM.',lupri)
           ENDIF
           INTEGRAL%ADMMBASISFILE   = .TRUE.
           INTEGRAL%PRINT_EK3       = .TRUE. !to calculate the full exchange matrix
           INTEGRAL%ADMM_separateX  = .TRUE. !to obtain the pure Exchange functional energy of full density 
        CASE ('.ADMM-2ERI'); ! EXPERIMENTAL
           INTEGRAL%ADMM_2ERI       = .TRUE.
        CASE ('.PRINT_EK3'); ! EXPERIMENTAL
        ! calculate and print full Exchange when doing ADMM exchange approx.
        ! > Debugging purpose only
           INTEGRAL%PRINT_EK3       = .TRUE.
        CASE ('.ADMM-K-METRIC'); ! EXPERIMENTAL
        ! calculate and print the residual error in the exchange metric
        ! > Development purpose only
           INTEGRAL%ADMMexchangeMetric = .TRUE.
        CASE ('.SREXC'); 
           INTEGRAL%MBIE_SCREEN = .TRUE.
           INTEGRAL%SR_EXCHANGE = .TRUE.
           READ(LUCMD,*) INTEGRAL%CAMmu
        CASE ('.LREXC'); 
           INTEGRAL%LR_EXCHANGE = .TRUE.
        CASE ('.LREXC-DF'); 
           INTEGRAL%LR_EXCHANGE = .TRUE.
           INTEGRAL%LR_EXCHANGE_DF = .TRUE.
        CASE ('.LREXC-PARI'); 
           INTEGRAL%LR_EXCHANGE = .TRUE.
           INTEGRAL%LR_EXCHANGE_PARI = .TRUE.
        CASE ('.NOFAMILY');
           INTEGRAL%NOFAMILY = .TRUE.
        CASE ('.OVERLAP-DF-J'); INTEGRAL%OVERLAP_DF_J=.TRUE.
        CASE ('.PARI'); INTEGRAL%PARI_J=.TRUE.;INTEGRAL%PARI_K=.TRUE.
        CASE ('.PARI-J'); INTEGRAL%PARI_J=.TRUE.
        CASE ('.EASY-PARI'); INTEGRAL%SIMPLE_PARI=.TRUE.
        CASE ('.PARI-K');  INTEGRAL%PARI_K=.TRUE.
        CASE ('.MOPARI-K'); INTEGRAL%MOPARI_K=.TRUE.
        CASE ('.DF-K');    INTEGRAL%DF_K=.TRUE.
        CASE ('.NR-PARI'); INTEGRAL%NON_ROBUST_PARI=.TRUE.
        CASE ('.PARI-UNCONSTRAINED');
           INTEGRAL%PARI_CHARGE=.FALSE.
           INTEGRAL%PARI_DIPOLE=.FALSE.
        CASE ('.PARI-CHARGE');
           INTEGRAL%PARI_CHARGE=.TRUE.
           INTEGRAL%PARI_DIPOLE=.FALSE.
        CASE ('.PARI-DIPOLE');
           INTEGRAL%PARI_CHARGE=.FALSE.
           INTEGRAL%PARI_DIPOLE=.TRUE.
        CASE ('.THR_PARI');   READ(LUCMD,*) INTEGRAL%PARI_THRESHOLD
           INTEGRAL%PARI_THRESHOLD = INTEGRAL%PARI_THRESHOLD/INTEGRAL%THRESHOLD 
        CASE ('.THR_CS');   READ(LUCMD,*) INTEGRAL%CS_THRESHOLD
           !SET IN RELATION TO THE ONE THRESHOLD
           INTEGRAL%CS_THRESHOLD = INTEGRAL%CS_THRESHOLD/INTEGRAL%THRESHOLD 
        CASE ('.THR_PS');   READ(LUCMD,*) INTEGRAL%PS_THRESHOLD
           !SET IN RELATION TO THE ONE THRESHOLD
           INTEGRAL%PS_THRESHOLD = INTEGRAL%PS_THRESHOLD/INTEGRAL%THRESHOLD 
        CASE ('.TIMINGS'); INTEGRAL%TIMINGS = .TRUE.
        CASE ('.UNCONT'); INTEGRAL%UNCONT = .TRUE.
        CASE ('.NOSEGMENT'); INTEGRAL%NOSEGMENT = .TRUE.
        CASE ('.MIXEDOVERLAP'); INTEGRAL%MIXEDOVERLAP = .TRUE.
        CASE ('.NO MM_FILES'); INTEGRAL%NO_MMFILES = .TRUE.
        CASE ('.DENSFIT'); integral%DENSFIT = .TRUE.                 
        CASE ('.INTEREST'); integral%INTEREST = .TRUE.                 
        CASE ('.RUNMM');  integral%FMM = .TRUE.
        CASE ('.DECGRA'); integral%run_dec_gradient_test=.true.
        CASE ('.MEMDIST'); integral%MEMDIST = .TRUE.
        CASE DEFAULT
           WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                & '" not recognized in **INTEGRALS readin.'
           print*,'Keyword ',WORD
           CALL lsQUIT('Illegal keyword in **INTEGRAL.',lupri)
        END SELECT
     ENDIF
     IF (WORD(1:2) == '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
     IF(PROMPT(1:1) .EQ. '*')THEN
        SELECT CASE(WORD(1:7))
        CASE('*READIN'); WRITE(LUPRI,*)'*READIN is not yet implemented in new int driver'
        CASE('*ONEINT'); WRITE(LUPRI,*)'*ONEINT is not yet implemented in new int driver'
        CASE('*TWOINT'); WRITE(LUPRI,*)'*TWOINT is not yet implemented in new int driver'
        CASE('*SUPINT'); WRITE(LUPRI,*)'*SUPINT is not yet implemented in new int driver'
        CASE('*ER2INT'); WRITE(LUPRI,*)'*ER2INT is not yet implemented in new int driver'
        CASE('*SORINT'); WRITE(LUPRI,*)'*SORINT is not yet implemented in new int driver'
        CASE('*DENFIT') 
           CALL READ_INTEGRALS_DENFIT_INPUT(LUPRI,LUCMD,integral,WORD)
           READWORD=.FALSE.
        CASE('*FMM   ')
           CALL READ_INTEGRALS_FMM_INPUT(LUPRI,LUCMD,integral,WORD)
           READWORD=.FALSE.
        CASE('*FCK3  '); 
           CALL READ_INTEGRALS_FCK3_INPUT(LUPRI,LUCMD,integral,WORD)
           READWORD=.FALSE.
           WRITE(LUPRI,*)'WORD LEAVING FCK3 ',WORD
        CASE('*PERIOD'); WRITE(LUPRI,*)'*PERIOD is not yet implemented in new int driver'
        CASE('*RIFOCK'); WRITE(LUPRI,*)'*RIFOCK is not yet implemented in new int driver'
        CASE('*END OF'); 
           READWORD=.FALSE.
           RETURN
        CASE DEFAULT
           WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                & '" not recognized in **INTEGRALS'
                          CALL lsQUIT('Illegal keyword in INTEGRALS.',lupri)
        END SELECT
     ENDIF
  ENDDO
END subroutine INTEGRAL_INPUT

#ifdef MOD_UNRELEASED
!> \brief Read the **GEOHESSIAN section in the input file LSDALTON.INP
!> \author Patrick Merlot
!> \date 14/09/2012
subroutine GEOHESSIAN_INPUT(geoHessian,readword,word,lucmd,lupri)
  implicit none
  LOGICAL,intent(inout)                :: READWORD
  TYPE(geoHessianConfig),intent(inout)   :: geoHessian
  character(len=80)  :: WORD
  INTEGER,intent(in) :: LUCMD !Logical unit number for the daltoninput
  INTEGER,intent(in) :: LUPRI !Logical unit number for the daltonoutput file
  !
  INTEGER            :: IDUMMY
  character(len=2)   :: PROMPT
  !
  DO   
     IF(READWORD) THEN
        READ (LUCMD, '(A40)') WORD
        READWORD=.TRUE.
     ENDIF
     PROMPT = WORD(1:2)
     IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
     IF(PROMPT .EQ. '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
     IF(PROMPT(1:1) .EQ. '.') THEN
        SELECT CASE(WORD) 
        CASE ('.TEST'); GEOHESSIAN%testContrib=.TRUE.
#ifdef BUILD_GEN1INT_LSDALTON
           GEOHESSIAN%DebugGen1Int = .TRUE.
!#else
!           call lsquit('.DEBUGGEN1INT requires OpenRSP -DBUILD_GEN1INT_LSDALTON',-1)
#endif
        CASE ('.INTPRINT');  READ(LUCMD,*) GEOHESSIAN%IntPrint
        CASE DEFAULT
           WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                & '" not recognized in **GEOHESSIAN readin.'
           print*,'Keyword ',WORD
           CALL lsQUIT('Illegal keyword in **GEOHESSIAN.',lupri)
        END SELECT
     ENDIF
     IF (WORD(1:2) == '**') THEN
        READWORD=.FALSE.
        EXIT
     ENDIF
     IF(PROMPT(1:1) .EQ. '*')THEN
        SELECT CASE(WORD(1:7))
        CASE('*PARTIAL'); WRITE(LUPRI,*)'*PARTIAL Hessian is not yet implemented in the Hessian driver'
        CASE('*END OF'); 
           READWORD=.FALSE.
           RETURN
        CASE DEFAULT
           WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                & '" not recognized in **GEOHESSIAN'
                          CALL lsQUIT('Illegal keyword in GEOHESSIAN.',lupri)
        END SELECT
     ENDIF
  ENDDO
END subroutine GEOHESSIAN_INPUT
#endif

!> \brief Read the **INFO section in input file LSDALTON.INP and set configuration structure accordingly.
!> \author S. Host
!> \date March 2010
SUBROUTINE config_info_input(config,lucmd,readword,word)
  implicit none
  !> Contains info, settings and data for entire calculation
  type(configItem),intent(inout) :: config
  !> Logical unit number for LSDALTON.INP
  integer,intent(in) :: lucmd
  LOGICAL,intent(inout)                :: READWORD
  character(len=80),intent(inout)  :: WORD
  character(len=2)   :: PROMPT
  integer :: i
  integer(kind=8) :: PrintMemoryLowerLimit

  INFOLOOP: DO   
     IF(READWORD) THEN
        READ (LUCMD, '(A40)') WORD
        READWORD=.TRUE.
     ENDIF
     PROMPT = WORD(1:2)
     IF ((PROMPT(1:1) .EQ. '!') .OR. (PROMPT(1:1) .EQ. '#')) CYCLE
     IF(PROMPT .EQ. '**') THEN
        READWORD=.FALSE.
        EXIT INFOLOOP
     ENDIF
     IF (PROMPT(1:1) == '*') THEN
        READWORD=.FALSE.
        EXIT INFOLOOP
     ENDIF
     SELECT CASE(WORD)
     CASE('.PRINTMEMORY');  
        config%PrintMemory = .TRUE.
        call Set_MemModParamPrintMemory(config%PrintMemory,config%lupri)
     CASE('.PRINTMEMORYLOWERLIMIT'); 
        !Do not print stack everytime the Maximum Allocated memory is increased
        !only print if Maximum Allocated memory is greater then PrintMemoryLowerLimit
        !Give input in byte
        IF(config%PrintMemory)THEN
           READ(LUCMD,*) PrintMemoryLowerLimit
           call Set_PrintMemoryLowerLimit(PrintMemoryLowerLimit)
        ELSE
           call lsquit('Error .PRINTMEMORYLOWERLIMIT requires .PRINTMEMORY',-1)
        ENDIF
     CASE('.DEBUG_SCF_MEM')
        call Set_PrintSCFmemory(.TRUE.)
        call setPrintDFTmem(.TRUE.)
     CASE('.DEBUG_MPI_MEM')
        config%mpi_mem_monitor = .true.
     CASE('.DEBUG_ARH_LINTRA')
        config%solver%DEBUG_ARH_LINTRA = .true.
     CASE('.DEBUG_ARH_PRECOND')
        config%solver%DEBUG_ARH_PRECOND = .true.
     CASE('.DEBUG_CONVERT')
        config%opt%DEBUG_CONVERT = .true.
        READ(LUCMD,*) config%opt%cfg_which_conversion 
        !CASE('.DEBUG_DCHANGE')
        !     DEBUG_DCHANGE = .true.
     CASE('.DEBUG_DD')
        READ(LUCMD,*) config%solver%cfg_nits_debug
        config%decomp%cfg_check_converged_solution = .true.
        config%solver%DEBUG_DD = .true.
     CASE('.DEBUG_DD_LINTRA')
        config%solver%DEBUG_DD_LINTRA = .true.
     CASE('.DEBUG_DD_HOMOLUMO')
        config%solver%DEBUG_DD_HOMOLUMO = .true.
     CASE('.DEBUG_DIAG_REDSPACE')
        config%solver%DEBUG_DIAG_REDSPACE = .true.
     CASE('.DEBUG_DIAG_HESSIAN')
        config%opt%DEBUG_DIAG_HESSIAN = .true.
!     CASE('.DEBUG_HESSIAN')
!        config%solver%DEBUG_HESSIAN = .true.
     CASE('.DEBUG_HESSIAN_EXACT')
        config%solver%DEBUG_HESSIAN_EXACT = .true. ; config%solver%DEBUG_HESSIAN = .true.
     CASE('.DEBUG_IDEMPOTENCY')
        config%diag%DEBUG_IDEMPOTENCY = .true.
        !CASE('.DEBUG_OAO_GRADIENT')
        !     DEBUG_OAO_GRADIENT = .true.
        !CASE('.DEBUG_RH_MU_E')
        !     DEBUG_RH_MU_E = .true.
        !     READ(LUCMD,*) cfg_nits_debug,  cfg_mu_max_debug 
        !     READ(LUCMD,*) (cfg_its_debug(i),i=1,cfg_nits_debug)
        !CASE('.DEBUG_RSP')
        !     DEBUG_RSP = .true.
        !CASE('.DEBUG_RSP_LINSCA')
        !     DEBUG_RSP_LINSCA = .true.
     CASE('.INFO_CROP')
        config%solver%INFO_CROP = .true.
     CASE('.INFO_LEVELSHIFT')
        config%solver%INFO_LEVELSHIFT    = .true.
     CASE('.INFO_STABILITY')
        config%decomp%INFO_STABILITY      = .true.
     CASE('.INFO_DIIS')
        config%av%INFO_DIIS = .true.
        config%av%INFO_WEIGHT_FINAL = .true. 
     CASE('.INFO_LINEQ')
        config%solver%INFO_LINEQ = .true.
     CASE('.INFO_MATOP')
        config%opt%INFO_MATOP = .true.
        !CASE('.INFO_ORB_E')
        !     INFO_ORB_E = .true.
     CASE('.INFO_STABILITY_REDSPACE')
        config%decomp%info_stability_redspace = .true.
     CASE('.INFO_RH')
        config%diag%INFO_RH_ITERATIONS    = .true.
        config%diag%INFO_RH_MU            = .true. 
        config%diag%INFO_RH_GAP           = .true.
     CASE('.INFO_RH_DETAIL')
        !config%diag%INFO_RH_EPRED         = .true.
        config%diag%INFO_RH_GAP           = .true.
        config%diag%INFO_RH_GRADIENT      = .true.
        config%diag%INFO_RH_ITERATIONS    = .true.
        config%diag%INFO_RH_MU            = .true.
     CASE('.INFO_RSP')
        config%response%rspsolverinput%INFO_RSP = .true.
     CASE('.INFO_RSP_REDSPACE')
        config%response%rspsolverinput%INFO_RSP_REDSPACE = .true.
        !CASE('.INFO_TIME_MAT_OPERATIONS')
        !     INFO_TIME_MAT_OPERATIONS = .true.
        !     call mat_timings
     CASE('.INFO_WEIGHT')
        config%av%INFO_WEIGHT_FINAL = .true.
     CASE DEFAULT
        WRITE(config%LUPRI,*) ' Keyword ',WORD,' not recognized in config_info_input'
        CALL lsQUIT('Illegal keyword in config_info_input.',config%lupri)
     END SELECT
  ENDDO INFOLOOP
END SUBROUTINE config_info_input



!> \brief Read the **RESPONS section under *LINSCA in input file LSDALTON.INP and set configuration structure accordingly.
!> \author S. Host
!> \date March 2010
SUBROUTINE config_rsp_input(config,lucmd,readword,WORD)
  implicit none
  !> Logical for keeping track of when to read
  LOGICAL,intent(inout)                :: READWORD
  !> Contains info, settings and data for entire calculation, including response
  type(configItem),intent(inout) :: config
  !> Logical unit number for LSDALTON.INP
  integer,intent(in) :: lucmd
  character(len=80) :: word
  character(len=8) :: labels(2)
  character(len=8) :: QRlabels(3)
  integer :: i,j,k,n,nops,nlabel
  logical :: xcomp,ycomp,zcomp
!temporary variables, fix to make test cases run!
  logical :: cfg_run_pdbs, cfg_rsp_run_quadratic
  integer :: Bterm_nr=0
  integer,allocatable     :: Bterm_index(:)
  integer :: STEPS
  real(realk) :: TMP1,TMP2,STEP,TMP3

  cfg_run_pdbs = .false.
  cfg_rsp_run_quadratic = .false.
  nops = 0
  RSPLOOP: DO
     IF(READWORD) THEN
        READ (LUCMD, '(A40)') WORD
        READWORD=.TRUE.
     ENDIF
     IF ((WORD(1:1) .EQ. '!') .OR. (WORD(1:1) .EQ. '#')) CYCLE
     IF(WORD(1:2) .EQ. '**') THEN
        READWORD=.FALSE.
        EXIT RSPLOOP
     ENDIF
     IF(WORD(1:13) == '*END OF INPUT') THEN
        backspace(LUCMD)
        EXIT RSPLOOP
     END IF

     if (WORD(1:1) == '*') then
       !which type of response is wanted??
       SELECT CASE(WORD)
       CASE('*NoOpenRSP')
          config%response%noOpenRSP = .TRUE. !Use LSDALTON own Response module
       CASE('*DIPOLE')
          config%response%tasks%doDipole=.true.
       CASE('*DIPOLEMOMENTMATRIX')
          config%response%tasks%doDipoleMatrix=.true.
          config%response%tasks%doResponse=.true.
       ! Kasper K
       CASE('*ALPHA')
           config%response%tasks%doALPHA=.true.
           config%response%tasks%doResponse=.true.
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.BFREQ')
                 config%response%alphainput%real_frequencies_in_input = .true.
                 READ(LUCMD,*) config%response%alphainput%nfreq 
                 allocate(config%response%alphainput%bfreq(config%response%alphainput%nfreq))
                 read(LUCMD,*) config%response%alphainput%bfreq(1:config%response%alphainput%nfreq)
              CASE('.IMBFREQ')
                 config%response%alphainput%imag_frequencies_in_input = .true.
                 READ(LUCMD,*) config%response%alphainput%nimfreq 
                 allocate(config%response%alphainput%imbfreq(config%response%alphainput%nimfreq))
                 read(LUCMD,*) config%response%alphainput%imbfreq(1:config%response%alphainput%nimfreq)
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *ALPHA input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*BETA')
           config%response%tasks%doBETA=.true.
           config%response%tasks%doResponse=.true.
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.BFREQ')
                 config%response%betainput%real_bfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%betainput%nbfreq
                 allocate(config%response%betainput%bfreq(config%response%betainput%nbfreq))
                 read(LUCMD,*) config%response%betainput%bfreq(1:config%response%betainput%nbfreq)
              CASE('.IMBFREQ')
                 config%response%betainput%imag_bfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%betainput%nimbfreq 
                 allocate(config%response%betainput%imbfreq(config%response%betainput%nimbfreq))
                 read(LUCMD,*) config%response%betainput%imbfreq(1:config%response%betainput%nimbfreq)
              CASE('.CFREQ')
                 config%response%betainput%real_cfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%betainput%ncfreq
                 allocate(config%response%betainput%cfreq(config%response%betainput%ncfreq))
                 read(LUCMD,*) config%response%betainput%cfreq(1:config%response%betainput%ncfreq)
              CASE('.IMCFREQ')
                 config%response%betainput%imag_cfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%betainput%nimcfreq 
                 allocate(config%response%betainput%imcfreq(config%response%betainput%nimcfreq))
                 read(LUCMD,*) config%response%betainput%imcfreq(1:config%response%betainput%nimcfreq)
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *BETA input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*GAMMA')
           config%response%tasks%doGAMMA=.true.
           config%response%tasks%doResponse=.true.
           config%response%rspsolverinput%rsp_cmplxnew = .true.
           config%response%rspsolverinput%rsp_cpp = .false.
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.BFREQ')
                 config%response%gammainput%real_bfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%nbfreq
                 allocate(config%response%gammainput%bfreq(config%response%gammainput%nbfreq))
                 read(LUCMD,*) config%response%gammainput%bfreq(1:config%response%gammainput%nbfreq)
              CASE('.IMBFREQ')
                 config%response%gammainput%imag_bfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%nimbfreq 
                 allocate(config%response%gammainput%imbfreq(config%response%gammainput%nimbfreq))
                 read(LUCMD,*) config%response%gammainput%imbfreq(1:config%response%gammainput%nimbfreq)
              CASE('.CFREQ')
                 config%response%gammainput%real_cfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%ncfreq
                 allocate(config%response%gammainput%cfreq(config%response%gammainput%ncfreq))
                 read(LUCMD,*) config%response%gammainput%cfreq(1:config%response%gammainput%ncfreq)
              CASE('.IMCFREQ')
                 config%response%gammainput%imag_cfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%nimcfreq 
                 allocate(config%response%gammainput%imcfreq(config%response%gammainput%nimcfreq))
                 read(LUCMD,*) config%response%gammainput%imcfreq(1:config%response%gammainput%nimcfreq)
              CASE('.DFREQ')
                 config%response%gammainput%real_dfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%ndfreq
                 allocate(config%response%gammainput%dfreq(config%response%gammainput%ndfreq))
                 read(LUCMD,*) config%response%gammainput%dfreq(1:config%response%gammainput%ndfreq)
              CASE('.IMDFREQ')
                 config%response%gammainput%imag_dfrequencies_in_input = .true.
                 READ(LUCMD,*) config%response%gammainput%nimdfreq 
                 allocate(config%response%gammainput%imdfreq(config%response%gammainput%nimdfreq))
                 read(LUCMD,*) config%response%gammainput%imdfreq(1:config%response%gammainput%nimdfreq)
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *GAMMA input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*TPA')
           config%response%tasks%doTPA=.true.
           config%response%tasks%doResponse=.true.
           ! Sanity check 1. We only run TPA calculations if .NEXCIT have been specified!
           TPASanityCheck1: if(config%decomp%cfg_rsp_nexcit == 0) then
              write(config%lupri,*) 'Error in TPA input: To run a TPA calculation you must&
                   & define the maximum number of excited states using the .NEXCIT&
                   & keyword under **RESPONSE before specifying *TPA!'
              write(config%lupri,*) 'Example of TPA response input for the first five excited states:'
              write(config%lupri,*) '**RESPONSE'
              write(config%lupri,*) '.NEXCIT'
              write(config%lupri,*) '5'
              write(config%lupri,*) '*TPA'
                  call lsQUIT('Error in TPA input: To run a TPA calculation you must&
                 & define the maximum number of excited states using the .NEXCIT&
                 & keyword under **RESPONSE before specifying *TPA!&
                 & See end of output file for input example!', config%lupri)
           end if TPASanityCheck1
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.EXSTATES')
                 config%response%tpainput%specific_states_in_input = .true.
                 READ(LUCMD,*) config%response%tpainput%tpa_nexci
                 allocate(config%response%tpainput%ExStates(config%response%tpainput%tpa_nexci))
                 read(LUCMD,*) config%response%tpainput%ExStates(1:config%response%tpainput%tpa_nexci)
                 ! Sanity check 2 
                 do i=1,config%response%tpainput%tpa_nexci
                      TPASanityCheck2: if(config%response%tpainput%ExStates(i) > &
                                &config%decomp%cfg_rsp_nexcit) then
                   write(config%lupri,*) 'Error in .EXSTATES in *TPA input: The maximum excited state&
                          & defined by .EXSTATES must be less than or equal to the number of&
                          & excited states defined by .NEXCIT under **RESPONSE.'
                   write(config%lupri,*) 'Example of TPA response input for the four excited states 1,3, 5, and 8:'
                   write(config%lupri,*) '**RESPONSE'
                   write(config%lupri,*) '.NEXCIT'
                   write(config%lupri,*) '8'
                   write(config%lupri,*) '*TPA'
                   write(config%lupri,*) '.EXSTATES'
                   write(config%lupri,*) '4'
                   write(config%lupri,*) '1 3 5 8'
                  CALL lsQUIT('Error in .EXSTATES in *TPA input: The maximum excited state&
                     & defined by .EXSTATES must be less than or equal the number of excited states defined by&
                     & .NEXCIT under **RESPONSE!&
                     & See end of output file for input example!', config%lupri)
               end if TPASanityCheck2
                end do
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *TPA input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*DAMPED_TPA')
           config%response%tasks%doResponse=.true.
           config%response%tasks%doDTPA=.true.
           config%response%rspsolverinput%rsp_cmplxnew = .true.
           config%response%rspsolverinput%rsp_cpp = .false.
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.OPFREQ') ! One-photon frequencies
                 READ(LUCMD,*) config%response%dtpainput%nfreq
                 allocate(config%response%dtpainput%freq(config%response%dtpainput%nfreq))
                 READ(LUCMD,*) config%response%dtpainput%freq(1:config%response%dtpainput%nfreq)
              CASE('.GAMMA') ! One-photon frequencies
                 config%response%dtpainput%gamma_specified=.true.
                 READ(LUCMD,*) config%response%dtpainput%gamma
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *DAMPED_TPA input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*ESDIPOLE')
           config%response%tasks%doESD=.true.
           config%response%tasks%doResponse=.true.
           ! Sanity check 1. We only run ESD calculations if .NEXCIT have been specified!
           ESDSanityCheck1: if(config%decomp%cfg_rsp_nexcit == 0) then
              write(config%lupri,*) 'Error in excited state gradient input: &
                   & To run an ESD calculation you must&
                   & define the maximum number of excited states using the .NEXCIT&
                   & keyword under **RESPONSE before specifying *ESDIPOLE!'
              write(config%lupri,*) 'Example of ESD response input for the first five excited states:'
              write(config%lupri,*) '**RESPONSE'
              write(config%lupri,*) '.NEXCIT'
              write(config%lupri,*) '5'
              write(config%lupri,*) '*ESDIPOLE'
                  call lsQUIT('Error in excited state gradient input: To run an ESD calculation you must&
                 & define the maximum number of excited states using the .NEXCIT&
                 & keyword under **RESPONSE before specifying *ESDIPOLE!&
                 & See end of output file for input example!', config%lupri)
               end if ESDSanityCheck1
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.EXSTATES')
                 config%response%esdinput%specific_states_in_input = .true.
                 READ(LUCMD,*) config%response%esdinput%esd_nexci
                 allocate(config%response%esdinput%ExStates(config%response%esdinput%esd_nexci))
                 read(LUCMD,*) config%response%esdinput%ExStates(1:config%response%esdinput%esd_nexci)
                 ! Sanity check 2 
                 do i=1,config%response%esdinput%esd_nexci
                      ESDSanityCheck2: if(config%response%esdinput%ExStates(i) > &
                              &config%decomp%cfg_rsp_nexcit) then
                   write(config%lupri,*) 'Error in .EXSTATES in * input: The maximum excited state&
                          & defined by .EXSTATES must be less than or equal to the number of&
                          & excited states defined by .NEXCIT under **RESPONSE.'
                   write(config%lupri,*) 'Example of excited state gradient response input &
                                &for the four excited states 1,3, 5, and 8:'
                   write(config%lupri,*) '**RESPONSE'
                   write(config%lupri,*) '.NEXCIT'
                   write(config%lupri,*) '8'
                   write(config%lupri,*) '*ESDIPOLE'
                   write(config%lupri,*) '.EXSTATES'
                   write(config%lupri,*) '4'
                   write(config%lupri,*) '1 3 5 8'
                  CALL lsQUIT('Error in .EXSTATES in *ESDIPOLE input: The maximum excited state&
                     & defined by .EXSTATES must be less than or equal the number of excited states defined by&
                     & .NEXCIT under **RESPONSE!&
                     & See end of output file for input example!', config%lupri)
               end if ESDSanityCheck2
                end do
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *ESDIPOLE input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Kasper K
       CASE('*ESGRAD')
           config%response%tasks%doESG=.true.
           ! Sanity check 1. We only run ESG calculations if .NEXCIT have been specified!
           ESGSanityCheck1: if(config%decomp%cfg_rsp_nexcit == 0) then
              write(config%lupri,*) 'Error in excited state gradient input: &
                   & To run an ESG calculation you must&
                   & define the maximum number of excited states using the .NEXCIT&
                   & keyword under **RESPONSE before specifying *ESGRAD!'
              write(config%lupri,*) 'Example of ESG response input for the first five excited states:'
              write(config%lupri,*) '**RESPONSE'
              write(config%lupri,*) '.NEXCIT'
              write(config%lupri,*) '5'
              write(config%lupri,*) '*ESGRAD'
                  call lsQUIT('Error in excited state gradient input: To run an ESG calculation you must&
                 & define the maximum number of excited states using the .NEXCIT&
                 & keyword under **RESPONSE before specifying *ESGRAD!&
                 & See end of output file for input example!', config%lupri)
               end if ESGSanityCheck1
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:1) == '*') then ! New property or *END OF INPUT
                 backspace(LUCMD)
                 exit
              end if
              SELECT CASE(word)
              CASE('.EXSTATES')
                 config%response%esginput%specific_states_in_input = .true.
                 READ(LUCMD,*) config%response%esginput%esg_nexci
                 call mem_alloc(config%response%esginput%ExStates,config%response%esginput%esg_nexci)
                 read(LUCMD,*) config%response%esginput%ExStates(1:config%response%esginput%esg_nexci)
                 ! Sanity check 2 
                 do i=1,config%response%esginput%esg_nexci
                      ESGSanityCheck2: if(config%response%esginput%ExStates(i) > &
                              &config%decomp%cfg_rsp_nexcit) then
                   write(config%lupri,*) 'Error in .EXSTATES in * input: The maximum excited state&
                          & defined by .EXSTATES must be less than or equal to the number of&
                          & excited states defined by .NEXCIT under **RESPONSE.'
                   write(config%lupri,*) 'Example of excited state gradient response input &
                                &for the four excited states 1,3, 5, and 8:'
                   write(config%lupri,*) '**RESPONSE'
                   write(config%lupri,*) '.NEXCIT'
                   write(config%lupri,*) '8'
                   write(config%lupri,*) '*ESGRAD'
                   write(config%lupri,*) '.EXSTATES'
                   write(config%lupri,*) '4'
                   write(config%lupri,*) '1 3 5 8'
                  CALL lsQUIT('Error in .EXSTATES in *ESGRAD input: The maximum excited state&
                     & defined by .EXSTATES must be less than or equal the number of excited states defined by&
                     & .NEXCIT under **RESPONSE!&
                     & See end of output file for input example!', config%lupri)
               end if ESGSanityCheck2
                end do
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *ESGRAD input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       ! Thomas K and Kasper K
       CASE('*MOLGRA')
                    WRITE(config%LUPRI,*) 'Ground state molecular gradient calculations are carried out.'
                    config%response%tasks%doResponse=.true.
                    config%response%tasks%dograd = .True.
             ! Joanna K

#ifdef MOD_UNRELEASED
       CASE('*NUMHESS')
                    WRITE(config%LUPRI,*) 'Numerical Hessian calculations are carried out using the analytical gradient'
                    config%response%tasks%doNumHess = .True.
        CASE('*NUMGRAD')
                    WRITE(config%LUPRI,*) 'Numerical Gradient calculations are carried out'
                    config%response%tasks%doNumGrad = .True.
        CASE('*NUMGRADHESS')
                    WRITE(config%LUPRI,*) 'Numerical Hessian calculations are carried out using the numerical gradient'
                    config%response%tasks%doNumGradHess = .True.
#endif
        CASE('*SOLVER')
            do
               READ(LUCMD,'(A40)') word
               if(word(1:1) == '!' .or. word(1:1) == '#') cycle
               if(word(1:1) == '*') then ! New property or *END OF INPUT
                  backspace(LUCMD)
                  exit
               end if
               SELECT CASE(word)
               CASE('.COMPLEX')
                  config%response%rspsolverinput%rsp_complex = .true.
                  READ(LUCMD,*) config%response%rspsolverinput%rsp_gamma
               CASE('.CONVTHR') 
                  READ(LUCMD,*) config%response%rspsolverinput%rsp_thresh
               CASE('.SYM_SOLVER')
                  config%response%rspsolverinput%rsp_stdnew = .true.
                  config%response%rspsolverinput%rsp_cmplxnew = .true.
                  config%response%rspsolverinput%rsp_cpp = .false.
               CASE('.PAIR_SOLVER')
                  config%response%rspsolverinput%rsp_cmplxnew = .false.
                  config%response%rspsolverinput%rsp_cpp = .false.
               CASE('.MAXIT')
                  READ(LUCMD,*) config%response%rspsolverinput%rsp_maxit 
                 config%response%rspsolverinput%rsp_maxred=2*config%response%rspsolverinput%rsp_maxit 
               CASE('.MAXRED')
                  READ(LUCMD,*) config%response%rspsolverinput%rsp_maxred
               CASE('.CONVDYN')
                  READ(LUCMD,*) config%response%rspsolverinput%rsp_convdyn_type
                  config%response%rspsolverinput%rsp_convdyn =.true.
                  SELECT CASE(config%response%rspsolverinput%rsp_convdyn_type)
                  CASE('VTIGH'); config%response%rspsolverinput%rsp_conv_factor = 1.0E-4_realk
                  CASE('TIGHT'); config%response%rspsolverinput%rsp_conv_factor = 1.0E-3_realk
                  CASE('STAND'); config%response%rspsolverinput%rsp_conv_factor = 1.0E-2_realk
                  CASE('SLOPP'); config%response%rspsolverinput%rsp_conv_factor = 1.0E-1_realk
                  CASE DEFAULT
                  WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',config%response%rspsolverinput%rsp_convdyn_type,&
                                    & '" not recognized with .CONVDYN'
                  WRITE (config%LUPRI,*) 'Options are TIGHT, STANDARD, and SLOPPY.'
                  CALL lsQUIT('Illegal keyword with .CONVDYN',config%lupri)
                  END SELECT
               CASE('.QUIET')
                  config%response%rspsolverinput%rsp_quiet = .true.
               CASE('.AOPREC')
                  config%response%rspsolverinput%rsp_mo_precond = .false.
               CASE('.AOSTART')
                  config%response%rspsolverinput%rsp_mostart = .false.
               CASE('.NOPREC')
                  config%response%rspsolverinput%rsp_no_precond = .true.
#ifdef MOD_UNRELEASED
               CASE('.EXCVECLINEQ')
                  !use excitation vectors as additional initial guess for lineq
                  config%response%rspsolverinput%UseExcitationVecs = .TRUE.    
               CASE('.SVD')
                  !use SVD decomposition on the Residual in order to improve conv.
                  config%response%rspsolverinput%doSVD = .TRUE.     
                  !This is not fully tested more test on bigger systems required. 
#endif 
               CASE('.RESTEXC')
                  READ (LUCMD,*) config%response%rspsolverinput%rsp_restart_nexci
                  config%response%rspsolverinput%rsp_restart_exci = .true.
               CASE ('.NSTART');   READ(LUCMD,*) config%response%rspsolverinput%rsp_no_of_startvectors
                  config%response%rspsolverinput%rsp_startvectors = .true.  
                  config%decomp%cfg_startvectors = .TRUE.
                  config%decomp%cfg_no_of_startvectors = config%response%rspsolverinput%rsp_no_of_startvectors
               CASE('.DTHR')
                  !threshold for when excited states is considered degenerate
                  READ(LUCMD,*) config%response%rspsolverinput%degenerateTHR
               CASE DEFAULT
                  WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                       & '" not recognized in RESPONSE *SOLVER input.'
                  CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
                 
               END SELECT
            enddo


! KASPER K: IMPORTANT NOTE!!!
! As far as I know all response input keywords below 
! (except .NEXCIT and perhaps some solver thresholds) will be
! obsolote in LSDALTON and should be removed for the release.
! If anyone has objections please let me know!



       ! kasperk
       CASE('*AORESPONSE')
                    WRITE(config%LUPRI,*) 'AO-based response function &
                                   &calculations are carried out...'
                    !cfg_run_AOresponse=.true.
                    !allocate(AOR_input_temp(100))
                    !i=1
                    !do
                    !  READ(LUCMD,'(A40)') AOR_input_i
                    !    if(AOR_input_i(1:1) == '!' .or. AOR_input_i(1:1) == '#') cycle
                    !    if(AOR_input_i(1:15) == '*END AORESPONSE') exit
                    !    AOR_input_temp(i) = AOR_input_i
                    !    if(i>1) then
                    !      if(AOR_input_temp(i-1)(1:7) == '.NEXCIT') then
                    !        read(AOR_input_i, '(I4)') nexci_temp
                    !        if(config%decomp%cfg_rsp_nexcit < nexci_temp) config%decomp%cfg_rsp_nexcit = nexci_temp
                    !      endif
                    !    endif
                    !    i=i+1
                    !enddo
                    !AOR_stacksize=i-1
                    !allocate(AOR_input(AOR_stacksize))
                    !AOR_input(1:AOR_stacksize) = AOR_input_temp(1:AOR_stacksize)
                    !deallocate(AOR_input_temp)
       CASE('*QUASIMCD')
           config%response%tasks%doResponse=.true.
           config%response%tasks%doMCD=.true.
           config%response%rspsolverinput%rsp_cmplxnew = .true.
           config%response%rspsolverinput%rsp_cpp = .false.
           do
              READ(LUCMD,'(A40)') word
              if(word(1:1) == '!' .or. word(1:1) == '#') cycle
              if(word(1:15) == '*END QUASIMCD') exit
              SELECT CASE(word)
              CASE('.DEGENERATE')
                 !Degenerate states are possible! so we need to consider Aterms 
                 config%response%MCDinput%DegenerateStates = .TRUE.
                 config%response%rspsolverinput%DegenerateStates = .TRUE.
                 config%response%MCDinput%doAterms =.TRUE.    
              CASE('.DAMPEDXCOOR')
                 READ(LUCMD,*) config%response%MCDinput%nXcoor
                 call mem_alloc(config%response%MCDinput%Xcoor,config%response%MCDinput%nXcoor)
                 DO I=1,config%response%MCDinput%nXcoor
                    READ(LUCMD,*) config%response%MCDinput%Xcoor(I)
                 ENDDO
              CASE('.DAMPEDRANGE')
                 READ(LUCMD,*) TMP1
                 READ(LUCMD,*) TMP2
                 READ(LUCMD,*) STEP
                 STEPS=1
                 TMP3=TMP1
                 DO WHILE(TMP3 < TMP2)
                    STEPS=STEPS+1
                    TMP3=TMP3+STEP
                 ENDDO
                 config%response%MCDinput%nXcoor = STEPS
                 call mem_alloc(config%response%MCDinput%Xcoor,STEPS)
                 config%response%MCDinput%Xcoor(1) = TMP1
                 TMP3=TMP1
                 DO I=2,config%response%MCDinput%nXcoor
                    TMP3=TMP3+STEP
                    config%response%MCDinput%Xcoor(I) = TMP3
                 ENDDO
              CASE('.MCDEXCIT')
                 READ(LUCMD,*) config%response%MCDinput%nexci
              CASE('.MCDEXSTATES')
                 config%response%MCDinput%specific_states_in_input = .true.
                 READ(LUCMD,*) config%response%MCDinput%nMCDexci
                 allocate(config%response%MCDinput%ExStates(config%response%MCDinput%nMCDexci))
                 read(LUCMD,*) config%response%MCDinput%ExStates(1:config%response%MCDinput%nMCDexci)
                 ! Sanity check 
                 do i=1,config%response%MCDinput%nMCDexci
                  MCDSanityCheck: if(config%response%MCDinput%ExStates(i) > &
                                & config%response%MCDinput%nexci ) then
                   write(config%lupri,*) 'Error in .EXSTATES in *QUASIMCD input: The maximum excited state&
                          & defined by .EXSTATES must be less than or equal to the number of&
                          & excited states defined by .MCDEXCIT '
                   write(config%lupri,*) 'Example of MCD response input for the four excited states 1,3, 5, and 8:'
                   write(config%lupri,*) '**RESPONSE'
                   write(config%lupri,*) '*QUASIMCD'
                   write(config%lupri,*) '.MCDEXCIT'
                   write(config%lupri,*) '8'
                   write(config%lupri,*) '.MCDEXSTATES'
                   write(config%lupri,*) '4'
                   write(config%lupri,*) '1 3 5 8'
                   write(config%lupri,*) '*END QUASIMCD'
                   CALL lsQUIT('Error in .MCDEXSTATES in *QUASIMCD input: The maximum excited state&
                     & defined by .MCDEXSTATES must be less than or equal the number of excited states defined by&
                     & .MCDEXCIT!&
                     & See end of output file for input example!', config%lupri)
                  end if MCDSanityCheck
                 end do
              CASE('.NO LONDON')
                 config%response%MCDinput%london=.FALSE.
              CASE('.NO NONLONDON')
                 config%response%MCDinput%nolondon=.FALSE.                 
              CASE('.NO SIMULATE')
                 config%response%MCDinput%simulate=.FALSE.    
              CASE('.NO ATERM')
                 config%response%MCDinput%doAterms =.FALSE.    
              CASE('.NO BTERM')
                 config%response%MCDinput%doBterms =.FALSE.    
              CASE('.NO DAMPED')
                 config%response%MCDinput%dampedMCD =.FALSE.   
              CASE('.GAUSSIAN')
                 config%response%MCDinput%lorentz=.FALSE.                 
              CASE('.LINESHAPEPARAM')
                 config%response%MCDinput%useinputgamma=.TRUE.                 
                 READ(LUCMD,*) config%response%MCDinput%Gamma 
              CASE('.NVECFORPEAK')
                 READ(LUCMD,*) config%response%MCDinput%nVecForPeak
              CASE('.NSTEPS')
                 READ(LUCMD,*) config%response%MCDinput%Nsteps
              CASE DEFAULT
                 WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                      & '" not recognized in RESPONSE *QUASIMCD input.'
                 CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri)
              END SELECT
           enddo
       CASE('*LINRSP')
                    config%response%tasks%dolinrsp = .True.
                    config%response%tasks%doResponse=.true.
                    WRITE(config%LUPRI,*) 'Linear response calculations are carried out.'
                    nlabel = 2
       !SONIANEW
       CASE('*ESGRSP')
                    WRITE(config%LUPRI,'(/A)') 'Excited state gradient is available in LSDALTON&
                         & only using the *ESGRAD keyword!'
                    CALL lsQUIT('Use *ESGRAD to get the excited state gradient!',config%lupri)
                    !WRITE(LUPRI,*) 'Excited state gradient calculations are carried out'
                    !cfg_rsp_run_exgr = .true.
                    !cfg_rsp_grad_purify = .true.
                    !cfg_start_guess = 'H1DIAG'    
                    !IF(cfg_do_in_oao) cfg_start_guess = 'H1OAO'    
                    !nlabel = 1
                    !WRITE(LUPRI,'(A)')'The Calculation uses the density matrix from the previous gometry step'
                    !WRITE(LUPRI,'(A)')'as an initial guess in the new geometry step.'
                    !WRITE(LUPRI,'(A,A,A)')'If this fails the fallback is ',cfg_start_guess,' As Huckel has not been implemented'
       !SONIANEW
       CASE('*HESSMAG')
                    WRITE(config%LUPRI,*) 'Hessian/Magnetizability calcs carried out'
                    nlabel = 2
                    !cfg_rsp_run_mag = .true.
                    !!Sonia: Make separate set of options: use_eq_79 etc..
                    !!cfg_rsp_run_hes = .true.
       !THOMAS_NEW
       CASE('*SHIELD')
          config%response%tasks%doNMRshield=.true.
          config%response%tasks%doResponse=.true.
       CASE('*PDBS')
                    WRITE(config%LUPRI,*) 'Pertubation dependent basis set &
                    & calculations are carried out'
                    cfg_run_PDBS = .true.
                    if (cfg_rsp_run_quadratic) then
                       WRITE(config%LUPRI,*) 'You cannot run Quadratic&
                       & Response Calculation and a Perturbation dependent &
                       & basis set calculation at the same time' 
                       CALL lsQUIT('REMOVE *PDBS or *QUADRSP',-1)
                    endif
       CASE('*QUIET')
                    !cfg_rsp_quiet=.TRUE.
       CASE DEFAULT
          WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
             & '" not recognized in config_rsp_input.'
          !CALL lsQUIT('Illegal keyword in config_rsp_input.',config%lupri) 
       END SELECT

     else ! keywords not starting with *
       SELECT CASE(WORD)
       CASE('.NEXCIT') 
                       READ (LUCMD,*) config%decomp%cfg_rsp_nexcit
! KK, calculate one-photon absorption strengths for all
! excitation energies by default when .NEXCIT is specified.
                       config%response%tasks%doOPA=.true.
                       config%response%tasks%doResponse=.true.
       CASE('.FREQ')
                       !ALLOCATE(cfg_rsp_freq(10))
                       !READ(LUCMD,*) cfg_rsp_nfreq,(cfg_rsp_freq(i),i=1,cfg_rsp_nfreq)
                       !cfg_rsp_freqs_specified = .true.
       CASE('.IMFREQ')
                       !ALLOCATE(cfg_rsp_imfreq(10))
                       !READ(LUCMD,*) cfg_rsp_nimfreq,(cfg_rsp_imfreq(i),i=1,cfg_rsp_nimfreq)
                       !cfg_rsp_imfreqs_specified = .true.
       CASE('.DAMP_FREQ')
                       !READ(LUCMD,*) cfg_rsp_nfreq,cfg_rsp_minfreq,cfg_rsp_maxfreq
                       !ALLOCATE(cfg_rsp_freq(cfg_rsp_nfreq))
                       !if( (cfg_rsp_minfreq/=cfg_rsp_maxfreq) .and.&
                       !& (cfg_rsp_nfreq==1)) then
                       !   WRITE (LUPRI,'(/,3A,/)') ' Error. Wrong number of frequencies. '
                       !   CALL lsQUIT('Error in damp_freq')
                       !else
                       !   cfg_rsp_freq(1)=cfg_rsp_minfreq
                       !   cfg_rsp_freq(cfg_rsp_nfreq)=cfg_rsp_maxfreq
                       !   if (cfg_rsp_nfreq>2) then
                       !      k=0
                       !      do i=2,cfg_rsp_nfreq-1
                       !        k=k+1
                       !       cfg_rsp_freq(i)=cfg_rsp_minfreq+k*((cfg_rsp_maxfreq-cfg_rsp_minfreq)/(cfg_rsp_nfreq-1))
                       !        enddo
                       !        endif
                       !cfg_rsp_freqs_specified = .true.
                       !endif
       CASE('.DAMP_TPA')   !kasperk
! 1. Determine isotropically averaged damped TPA. Input:
! .DAMP_TPA
! 0
! 
! 2. Determine only specific TPA component(s), e.g. XXXX and XXYY. Input example:
! .DAMP_TPA
! 2                 ! Number of components
! XXXX              ! stored in tp_info(1)
! XXYY              ! stored in tp_info(2)
                       !cfg_rsp_run_damped_tpa = .true.
                       !READ(LUCMD,*) Ntp_comp
                       !if(Ntp_comp /= 0) then
                       !  allocate(tp_info(Ntp_comp))
                       !    do i=1,Ntp_comp
                       !      READ(LUCMD,*) tp_info(i)
                       !    enddo
                       !endif
       CASE('.OPERATOR')
                       !WRITE(LUPRI,*)
                       !WRITE(LUPRI,*) 'Response functions/transition moments will be calculated'
                       !WRITE(LUPRI,*) 'for the following (pairs of) operators:'
                       !WRITE(LUPRI,*) 
                       !do 
                       !  READ (LUCMD, '(a40)') word
                       !  if (word(1:1) == '.' .or. word(1:1) == '$') exit
                       !  if (word(1:1) == '*') exit
                       !  if (word(1:1) == '!' .or. word(1:1) == '#') cycle
                       !  nops = nops + 1
                       !  !TODO: check that it is within ..operators boundaries
                       !  !SONIANEW: initialized within given section (LINRSP or ESGRSP)
                       !  !nlabel = 2
                       !  call trim_string(word,nlabel,labels)
                       !  do i = 1,nlabel
                       !    cfg_rsp_operators(i,nops) = CFG_RSP_LABEL_TO_INT(labels(i))
                       !  enddo
                       !enddo
                       !cycle
       CASE('.DIPOLE')
                       !cfg_rsp_dipole = .true. 
       CASE('.QRSP OAO')
                       !cfg_QRSP_OAO=.true.
       CASE('.QRSP LRF')
                       !cfg_QR_do_LRF=.true.
       CASE('.ATERM')
                       !if(cfg_run_PDBS)then
                       !  cfg_run_Aterm=.true.
                       !  READ (LUCMD,*) Aterm_nr
                       !  !PDBSitems=PDBSitems+Aterm_nr
                       !  !PDBSitems set when it is examined how many degenerate
                       !  !states there is 
                       !  allocate(Aterm_index(Aterm_nr))
                       !  do i=1,Aterm_nr 
                       !     READ(LUCMD,*) Aterm_index(i)
                       !  enddo
                       !  IF(config%decomp%cfg_rsp_nexcit .EQ. 0)THEN
                       !     config%decomp%cfg_rsp_nexcit = Aterm_nr
                       !     DO i=1,Aterm_nr 
                       !        config%decomp%cfg_rsp_nexcit = MAX(config%decomp%cfg_rsp_nexcit,Aterm_index(i))
                       !     ENDDO
                       !  ENDIF
                       !endif
       CASE('.BTERM')
                       if(cfg_run_PDBS)then
                         !cfg_run_Bterm=.true.
                         READ (LUCMD,*) Bterm_nr
                         !PDBSitems=PDBSitems+Bterm_nr
                         allocate(Bterm_index(Bterm_nr))
                         do i=1,Bterm_nr 
                            READ(LUCMD,*) Bterm_index(i)
                         enddo
                         IF(config%decomp%cfg_rsp_nexcit .EQ. 0)THEN
                            config%decomp%cfg_rsp_nexcit = Bterm_nr
                            DO i=1,Bterm_nr 
                               config%decomp%cfg_rsp_nexcit = MAX(config%decomp%cfg_rsp_nexcit,Bterm_index(i))
                            ENDDO
                         ENDIF
                       endif
       CASE('.VERDET')
                       !if(cfg_run_PDBS)then
                       !  cfg_run_VERDET=.true.
                       !  READ (LUCMD,*) verdet_nr
                       !  PDBSitems=PDBSitems+verdet_nr
                       !  allocate(Verdet_freq(verdet_nr))
                       !  do i=1,verdet_nr 
                       !    READ(LUCMD,*) Verdet_freq(i)
                       !  enddo
                       !endif
       CASE('.PDBSINFO')
                       !cfg_PDBS_INFO=.true.
       CASE('.ESGRSP')
                       !cfg_run_Exgrad = .true. 
                       !READ (LUCMD,*) Exg_nr
                       !PDBSitems=PDBSitems+Exg_nr
                       !allocate(Exg_index(Exg_nr))
                       !do i=1,Exg_nr 
                       !  READ(LUCMD,*) Exg_index(i)
                       !enddo
       CASE('.RAMAN')
                       !cfg_run_RAMAN=.true.
                       !READ (LUCMD,*) raman_nr
                       !PDBSitems=PDBSitems+raman_nr
                       !allocate(Raman_freq(raman_nr))
                       !do i=1,raman_nr 
                       !  READ(LUCMD,*) Raman_freq(i)
                       !enddo
       CASE('.HBT')
                       !cfg_run_HerzbergTeller=.true.
                       !READ (LUCMD,*) HBT_nr
                       !PDBSitems=PDBSitems+HBT_nr
                       !allocate(HBT_index(HBT_nr))
                       !do i=1,HBT_nr 
                       !  READ(LUCMD,*) HBT_index(i)
                       !enddo
       CASE('.HOTFCHT')
                       !cfg_hotfcht = .true.
       !SONIANEW = select the excited state for the exc. state gradient
       CASE('.WHICHEXST') 
                       !READ (LUCMD,*) cfg_rsp_whichexst
       CASE('.RESTEXC')
                       !READ (LUCMD,*) cfg_restart_nexci; cfg_restart_exci = .true.
       !Joanna - rsp solver using conjugate gradient with optimal vectors
       CASE('.CGOP')
                      !cfg_rsp_cgop = .true.
                      !cfg_cg_truncate = .true.
                      !no_pairing=.true.
       CASE('.CGOP_NT')
                      !cfg_rsp_cgop = .true.
                      !cfg_cg_truncate = .false.
                      !no_pairing=.true.
       CASE('.CGOP_PAIR')      
                      !cfg_rsp_cgop = .true.       
       CASE('.CG_VEC')
                      !READ (LUCMD,*) cfg_red_truc
       CASE('.NO_CGOP_PREC')                 
               !cfg_cgop_prec=.false.
       !Joanna- complex response solver
       CASE('.COMPLEX')
                   READ (LUCMD,*) config%response%rspsolverinput%rsp_gamma
                         config%response%rspsolverinput%rsp_complex=.true.
                     !cfg_rsp_moprec= .true.
                     !cfg_rsp_mostart=.true.
       CASE('.COMP_OR')
                      !READ (LUCMD,*) cfg_rsp_gamma
                      !cfg_rsp_complex_or = .true.
       !FILIP:If requested, turn off the use of the density matrix from the previous geometry step
       !      as an initial guess for the new geometry step.
       !      gradients.
       CASE('.NODRESTART')
                      ! if (.not.cfg_rsp_run_grad) then
                      !    write(lupri,*)'.NODRESTART option requested, but *MOLGRA not turned on!'
                      !    call lsquit('.NODRESTART option requested, but *MOLGRA not turned on!')
                      ! endif
                      ! cfg_rsp_grad_drestart = .false.
       !FILIP:If requested, turn off McWeeny purification of the density matrix for the molecular
       !      gradients.
       CASE('.NOPURIFY')
                      ! if (.not.cfg_rsp_run_grad) then
                      !    write(lupri,*)'.NOPURIFY option requested, but *MOLGRA not turned on!'
                      !    call lsquit('.NOPURIFY option requested, but *MOLGRA not turned on!')
                      ! endif
                      ! cfg_rsp_grad_purify = .false.
       !FILIP:Change the threshold for the difference abs(2.0*trace-Nelectrons), below which 
       !      we disregard the McWeeny purified D and fall back to a standard initial guess.
       CASE('.THRNEL')
                      ! if (.not.cfg_rsp_run_grad) then
                      !    write(lupri,*)'.THRNEL option requested, but *MOLGRA not turned on!'
                      !    call lsquit('THRNEL option requested, but *MOLGRA not turned on!')
                      ! endif
                      ! READ (LUCMD,*) THRNEL
       CASE('.POLARIZ')
                      ! cfg_rsp_polariz = .true. !polarizability
       CASE('.HYPOLAR2')
                      ! cfg_rsp_hypolar2 = .true. !1st hyperpol., n+1
       CASE('.HYPOLAR')
                      ! cfg_rsp_hypolar = .true. !1st hyperpol., 2n+1
       CASE('.SECHYP3')
                      ! cfg_rsp_sechyp3 = .true. !2nd hyperpol., n+1
       CASE('.SECHYP')
                      ! cfg_rsp_sechyp = .true. !2nd hyperpol., 2n+1 (1+2+1)
       CASE('.SECHYP1')
                      ! cfg_rsp_sechyp1 = .true. !2nd hyperpol., 2n+1 (2+1+1)
       CASE('.MAGNET')
                      ! cfg_rsp_magnet = .true.   !ajt
       CASE('.EFGB')
                      ! cfg_rsp_efgb = .true.     !ajt
       CASE('.CME')
                      ! cfg_rsp_cme = .true.      !ajt
       CASE('.ROA')
                      ! cfg_rsp_roa = .true.      !ajt
       CASE('.CARS')
                      ! cfg_rsp_cars = .true.     !ajt
       CASE('.JONES')
                      ! cfg_rsp_jones = .true.    !ajt
       CASE('.VIBBETA')
                      ! cfg_rsp_vibbeta = .true.  !ajt
       CASE('.PROPTEST')
                      ! cfg_rsp_proptest = .true. !ajt
       CASE DEFAULT
          WRITE (config%LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
               & '" not recognized in config_rsp_input'
          CALL lsQUIT('Illegal keyword in **RESPONSE section',config%lupri)
       END SELECT
     endif
  ENDDO RSPLOOP

  if (config%decomp%cfg_rsp_nexcit > 0) then !Stinne
      !cfg_rsp_maxred = config%decomp%cfg_rsp_nexcit*cfg_rsp_maxred
  else
     !if (cfg_rsp_freqs_specified) then
     !   WRITE(LUPRI,*)
     !   WRITE(LUPRI, "('The response functions will be calculated for each of the ', i3, &
     !               & ' chosen frequencies.')") cfg_rsp_nfreq
     !  !kk and ajt: default allocation of zero cfg_rsp_imfreq
     !   if (.not.cfg_rsp_imfreqs_specified) then
     !      cfg_rsp_nimfreq = cfg_rsp_nfreq
     !      allocate(cfg_rsp_imfreq(cfg_rsp_nimfreq))
     !      cfg_rsp_imfreq = 0.0E0_realk
     !   endif
     !else
     !   WRITE(LUPRI,*)
     !   WRITE(LUPRI, "('No frequencies specified in input. Only static properties will be &
     !               & calculated.')") 
     !endif
  endif

  if (config%decomp%cfg_rsp_nexcit > 0) then
    !if (cfg_rsp_dipole) then
    !   do i = 1,3  !run over dip X,Y,Z
    !     nops = nops + 1
    !     cfg_rsp_operators(1,nops) = i
    !   enddo
    !endif
    !!SONIANEW: chosen excited state must be within required solutions
    !if (cfg_rsp_run_exgr) then
    !   if (cfg_rsp_whichexst > config%decomp%cfg_rsp_nexcit) then
    !      WRITE(LUPRI,*) 'You requested an excited state which has not been determined'
    !      WRITE(LUPRI,*) 'Calculation stops. Change NEXCIT or WHICHEXST and come back'
    !      STOP '*** ERROR IN ESGRSP INPUT ( WHICHEXST > NEXCIT) '
    !   elseif (cfg_rsp_whichexst == 0 ) then
    !      WRITE(LUPRI,*) 'No excited state number specified. Use default value'
    !      WRITE(LUPRI,*) 'cfg_rsp_whichexst = config%decomp%cfg_rsp_nexcit (highest state)'
    !      cfg_rsp_whichexst = config%decomp%cfg_rsp_nexcit
    !   end if
    !endif
  else
    !if (cfg_rsp_dipole) then
    !   do j = 1,3  !run over dipole_X,Y,Z - see top of file
    !     do i = j,3  !run over dipole_X,Y,Z
    !       nops = nops + 1
    !       cfg_rsp_operators(1,nops) = i
    !       cfg_rsp_operators(2,nops) = j
    !     enddo
    !   enddo
    !endif
  endif
  !number of "vectors" = number of operator pair to evaluate
  !cfg_rsp_nvec = nops  

END SUBROUTINE config_rsp_input

SUBROUTINE READ_WAVE_DFTINPUT(LUPRI,LUCMD,DALTON,WORD)
implicit none
TYPE(integralconfig)   :: DALTON
INTEGER            :: LUCMD !Logical unit number for the daltoninput
character(len=80),intent(out)  :: WORD
character(len=1)   :: PROMPT
INTEGER            :: LUPRI,n,I
character(len=80)  :: LINE
character(len=9),pointer  :: GRIDspec(:)

DO
   READ (LUCMD, '(A40)') WORD
   PROMPT = WORD(1:1)
   IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') CYCLE
   IF (PROMPT .EQ. '.') THEN
      SELECT CASE(WORD) 
      !PLEASE KEEP IN ALPHABETICAL ORDER
      CASE ('.ANGMIN'); READ(LUCMD,*) DALTON%DFT%ANGMIN
      CASE ('.ANGINT'); 
         READ(LUCMD,*) DALTON%DFT%ANGINT
         IF(DALTON%DFT%ANGINT .GT. 64) CALL lsquit('ANGINT > 64 not implemented',-1)
      CASE('.CARPAR');
         call lsquit('Cannot call DFTCARTESIANINPUT, only new integral code is compiled',-1)
      CASE ('.DFTAC');
         CALL lsQUIT(' Not implemented in the new code - if you know what this is, &
              & you could tell me and I might implement it. Thomas Kjaergaard',lupri)
         DALTON%DFT%DFTASC = .TRUE.
         DALTON%DFT%DFTPOT = .TRUE.
         READ (LUCMD,*) DALTON%DFT%DFTIPT, DALTON%DFT%DFTBR1, DALTON%DFT%DFTBR2
      CASE ('.DFTELS'); READ(LUCMD,*) DALTON%DFT%DFTELS
      CASE ('.DFTTHR'); READ(LUCMD,*) DALTON%DFT%DFTHR0,DALTON%DFT%DFTHRL, DALTON%DFT%DFTHRI, DALTON%DFT%RHOTHR
      CASE ('.MEMORY'); call setPrintDFTmem(.TRUE.)  
      CASE ('.LB94');
         DALTON%DFT%LB94=.TRUE.
      CASE ('.CS00');
         DALTON%DFT%CS00=.TRUE.
      CASE ('.CS00 ZND1');
         IF(.NOT.DALTON%DFT%CS00)&
           & CALL LSQUIT('INPUT ERROR .CS00 ZND1 must come after .CS00',-1)
         READ(LUCMD,*) DALTON%DFT%CS00ZND1
      CASE ('.CS00 ZND2');
         IF(.NOT.DALTON%DFT%CS00)&
           & CALL LSQUIT('INPUT ERROR .CS00 ZND2 must come after .CS00',-1)
         READ(LUCMD,*) DALTON%DFT%CS00ZND2
      CASE ('.CS00 SHIFT');
         DALTON%DFT%CS00=.TRUE.
         READ(LUCMD,*) DALTON%DFT%CS00SHIFT
      CASE ('.GRID TYPE'); 
         READ(LUCMD,'(A80)') LINE
         n=10 !max words should properly be 2
         allocate(GRIDspec(n))
         call TRIM_STRING(LINE,n,GRIDspec)
         DO I=1,n
            SELECT CASE(Gridspec(I)) 
               !=====================================================================
               !  These first three options are the radial grid options 
               !  (GC2,LMG,TURBO)
               !=====================================================================
            CASE('GC2')
               DALTON%DFT%RADIALGRID = 1
               DALTON%DFT%TURBO = 0
               DALTON%DFT%ZdependenMaxAng = .FALSE.
            CASE('LMG')
               DALTON%DFT%RADIALGRID = 2 !OLDDEFAULT
               DALTON%DFT%ZdependenMaxAng = .FALSE.
               DALTON%DFT%TURBO = 0
            CASE('TURBO') 
               DALTON%DFT%RADIALGRID = 3 !DEFAULT
               DALTON%DFT%ZdependenMaxAng = .TRUE.
               DALTON%DFT%TURBO = 1
               !=====================================================================
               !  These options are the partitioning options 
               !  (SSF,BECKE,BECKEORIG,BLOCK,BLOCKSSF)
               !=====================================================================
            CASE('SSF')
               DALTON%DFT%PARTITIONING = 1
               !TheGrid.partitioning = &partitioning_ssf;break;
            CASE('BECKE')
               DALTON%DFT%PARTITIONING = 2
               !TheGrid.partitioning = &partitioning_becke_corr
            CASE('BECKEORIG')
               DALTON%DFT%PARTITIONING = 3 !OLDDEFAULT
               !TheGrid.partitioning = &partitioning_becke_orig;
            CASE('BLOCK')
               DALTON%DFT%PARTITIONING = 4
               !TheGrid.partitioning = &partitioning_block;
            CASE('BLOCKSSF')
               DALTON%DFT%PARTITIONING = 5 !DEFAULT
               !TheGrid.partitioning = &partitioning_block_ssf
            CASE('CARTESIAN')
               DALTON%DFT%PARTITIONING = 6
               !gridType = GRID_TYPE_CARTESIAN; break;
            CASE DEFAULT
               WRITE (LUPRI,'(/,3A,/)') 'GRID TYPE: ',Gridspec(I),' not recognized'
               CALL lsQUIT('Illegal keyword in .GRID TYPE',lupri)
            END SELECT
         ENDDO
         deallocate(GRIDspec)
      CASE ('.NOPRUN'); DALTON%DFT%NOPRUN = .TRUE.
      CASE ('.RADINT'); READ(LUCMD,*) DALTON%DFT%RADINT
!===================================================================
      CASE ('.ULTRAC'); DALTON%DFT%RADINT = 2.15447E-7_realk; DALTON%DFT%ANGINT = 23; 
      CASE ('.COARSE'); DALTON%DFT%RADINT = 4.64159E-9_realk; DALTON%DFT%ANGINT = 29; 
      CASE ('.NORMAL'); DALTON%DFT%RADINT = 5.01187E-14_realk; DALTON%DFT%ANGINT = 35; 
      CASE ('.FINE'  ); DALTON%DFT%RADINT = 2.15443E-17_realk; DALTON%DFT%ANGINT = 47; 
      CASE ('.ULTRAF'); DALTON%DFT%RADINT = 2.15443E-17_realk; DALTON%DFT%ANGINT = 64;
!===================================================================
      CASE ('HARTRE'); 
         DALTON%DFT%DFTADD = .FALSE.
      CASE ('.GRID1' ) 
         DALTON%DFT%RADINT = 1E-5_realk; DALTON%DFT%ANGINT = 17; DALTON%DFT%ZdependenMaxAng=.TRUE.
         DALTON%DFT%TURBO = 1; DALTON%DFT%RADIALGRID = 3; DALTON%DFT%PARTITIONING = 4
      CASE ('.GRID2' ) 
         DALTON%DFT%RADINT = 2.15447E-7_realk; DALTON%DFT%ANGINT = 23; DALTON%DFT%ZdependenMaxAng=.TRUE.
         DALTON%DFT%TURBO = 1; DALTON%DFT%RADIALGRID = 3; DALTON%DFT%PARTITIONING = 4
      CASE ('.GRID3' ) 
         DALTON%DFT%RADINT = 4.64159E-9_realk; DALTON%DFT%ANGINT = 29; DALTON%DFT%ZdependenMaxAng=.TRUE.
         DALTON%DFT%TURBO = 1; DALTON%DFT%RADIALGRID = 3; DALTON%DFT%PARTITIONING = 4
      CASE ('.GRID4' ) 
         DALTON%DFT%RADINT = 5.01187E-14_realk; DALTON%DFT%ANGINT = 35; DALTON%DFT%ZdependenMaxAng=.TRUE.
         DALTON%DFT%TURBO = 1; DALTON%DFT%RADIALGRID = 3; DALTON%DFT%PARTITIONING = 4
      CASE ('.GRID5' ) 
         DALTON%DFT%RADINT = 2.15443E-17_realk; DALTON%DFT%ANGINT = 47; DALTON%DFT%ZdependenMaxAng=.TRUE.
         DALTON%DFT%TURBO = 1; DALTON%DFT%RADIALGRID = 3; DALTON%DFT%PARTITIONING = 4
      CASE ('.AOSAVE' ) 
         CALL lsQUIT('Not implemented in the new code, because I am not&
                    &convinced about the performance, convince me and I will. TK',lupri);
      CASE ('.HARDNESS' ); READ(LUCMD,*) DALTON%DFT%HRDNES
      CASE ('.DISPER' )
         DALTON%DFT%DODISP = .TRUE.
         DALTON%DFT%DODISP3 = .TRUE.
         DALTON%DFT%DO_DFTD3 = .TRUE.
         DALTON%DFT%DO_BJDAMP = .TRUE.
!AMT
      CASE ('.DFT-D2')
         DALTON%DFT%DODISP = .TRUE.
         DALTON%DFT%DODISP2 = .TRUE.
         DALTON%DFT%DO_DFTD2 = .TRUE.
      CASE ('.D2PAR')
         DALTON%DFT%DODISP = .TRUE.
         DALTON%DFT%L_INP_D2PAR = .TRUE.
         READ(LUCMD,*)DALTON%DFT%D2_s6_inp, DALTON%DFT%D2_alp_inp, DALTON%DFT%D2_rs6_inp
      CASE ('.DFT-D3')
         DALTON%DFT%DODISP = .TRUE.
         DALTON%DFT%DODISP3 = .TRUE.
         DALTON%DFT%DO_DFTD3 = .TRUE.
      CASE ('.DFT-D3BJ')
         DALTON%DFT%DODISP = .TRUE.
         DALTON%DFT%DODISP3 = .TRUE.
         DALTON%DFT%DO_DFTD3 = .TRUE.
         DALTON%DFT%DO_BJDAMP = .TRUE.
      CASE ('.3BODY')
         DALTON%DFT%DODISP = .TRUE.
         DALTON%DFT%DO_3BODY = .TRUE.
      CASE ('.D3PAR')
         DALTON%DFT%DODISP = .TRUE.
         DALTON%DFT%L_INP_D3PAR = .TRUE.
         READ(LUCMD,*) DALTON%DFT%D3_s6_inp, DALTON%DFT%D3_alp_inp, DALTON%DFT%D3_rs6_inp, & 
     &              DALTON%DFT%D3_rs18_inp, DALTON%DFT%D3_s18_inp
!AMT
      CASE DEFAULT
         WRITE (LUPRI,'(/,3A,/)') ' Keyword ',WORD,&
              & ' not recognized in *DFT INPUT'
         CALL lsQUIT('Illegal keyword in *DFT INPUT',lupri)
      END SELECT
   ENDIF
   IF (PROMPT .EQ. '*') EXIT
ENDDO

! add this including partitioning with refs and radial grid...
!WRITE (LUPRI,'(4X,A,17X,3F12.2)')&
!&        ' DFT LSint thresholds:', DALTON%DFTHR0, DALTON%DFTHRL, DALTON%DFTHRI 
!WRITE (LUPRI,'(4X,A,17X,3E12_realk.2)')&
!&        ' DFT LSint threshold for number of electrons: ', DALTON%DFTELS 
!WRITE (LUPRI,'(4X,A,F8.4,I4)')' DFT LSint radial quadrature accuracy/ang. &
!     & expansion order:',DALTON%RADINT,DALTON%ANGINT
!
END SUBROUTINE READ_WAVE_DFTINPUT

SUBROUTINE READ_INTEGRALS_DENFIT_INPUT(LUPRI,LUCMD,DALTON,word)
implicit none
TYPE(integralconfig)   :: DALTON
INTEGER            :: LUCMD !Logical unit number for the daltoninput
character(len=80),intent(out)  :: WORD
character(len=1)   :: PROMPT
INTEGER            :: LUPRI

! the old version is SUBROUTINE DFIINP

DO
  READ (LUCMD, '(A40)') WORD
  PROMPT = WORD(1:1)
  IF (PROMPT(1:1) .EQ. '!' .OR. PROMPT .EQ. '#') CYCLE
  IF (PROMPT .EQ. '.') THEN
    SELECT CASE(WORD) 
      CASE ('.DIATOM'); WRITE(LUPRI,*)'.DIATOM NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.DOBOX '); WRITE(LUPRI,*)'.DOBOX NOT IMPLEMENTED IN NEW DRIVER'  
      CASE ('.BOXORB'); WRITE(LUPRI,*)'.BOXORB NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.BUFSIZ'); WRITE(LUPRI,*)'.BUFSIZ NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.CONSTR'); WRITE(LUPRI,*)'.CONSTR NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOROBU'); WRITE(LUPRI,*)'.NOROBU NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.ATTENU'); WRITE(LUPRI,*)'.ATTENU NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.GAUDAM'); WRITE(LUPRI,*)'.GAUDAM NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NONLIN'); WRITE(LUPRI,*)'.NONLIN NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.HERMIT'); WRITE(LUPRI,*)'.HERMIT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.3CENT '); WRITE(LUPRI,*)'.3CENT  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.OAOINT'); WRITE(LUPRI,*)'.OAOINT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NO_J_M'); WRITE(LUPRI,*)'.NO_J_M NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NO_K_M'); WRITE(LUPRI,*)'.NO_K_M NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.RI-LIN'); WRITE(LUPRI,*)'.RI-LIN NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NONRMA'); WRITE(LUPRI,*)'.NONRMA NOT IMPLEMENTED IN NEW DRIVER' 
      CASE DEFAULT
                  WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                  & '" not recognized in DENFIT INTEGRAL.'
                 CALL lsQUIT('Illegal keyword in DENFIT INTEGRAL.',lupri)
    END SELECT
  ENDIF
  IF (PROMPT .EQ. '*') EXIT
ENDDO
END SUBROUTINE READ_INTEGRALS_DENFIT_INPUT

SUBROUTINE READ_INTEGRALS_FMM_INPUT(LUPRI,LUCMD,DALTON,word)
implicit none
TYPE(integralconfig)   :: DALTON
INTEGER            :: LUCMD !Logical unit number for the daltoninput
character(len=80),intent(out)  :: WORD
character(len=1)   :: PROMPT
INTEGER            :: LUPRI

! the old version is SUBROUTINE DFIINP

DO
  READ (LUCMD, '(A40)') WORD
  PROMPT = WORD(1:1)
  IF (PROMPT(1:1) .EQ. '!' .OR. PROMPT .EQ. '#') CYCLE
  IF (PROMPT .EQ. '.') THEN
    SELECT CASE(WORD(1:7)) 
      CASE ('.SKIP  '); WRITE(LUPRI,*)'.SKIP   NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.PRINT '); WRITE(LUPRI,*)'.PRINT  NOT IMPLEMENTED IN NEW DRIVER'  
      CASE ('.STOP  '); WRITE(LUPRI,*)'.STOP   NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.LMAX  '); READ(LUCMD,*) DALTON%MM_LMAX
      CASE ('.NLEVEL'); WRITE(LUPRI,*)'.NLEVEL NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.FQUAD '); WRITE(LUPRI,*)'.FQUAD  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NNQUAD'); WRITE(LUPRI,*)'.NNQUAD NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.BQUAD '); WRITE(LUPRI,*)'.BQUAD  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NLOGN '); WRITE(LUPRI,*)'.NLOGN  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.FMM   '); WRITE(LUPRI,*)'.FMM    NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.GRSRCH'); WRITE(LUPRI,*)'.GRSRCH NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.TLMAX '); READ(LUCMD,*) DALTON%MM_TLMAX
      CASE ('.UMAT  '); WRITE(LUPRI,*)'.UMAT   NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.BRFREE'); WRITE(LUPRI,*)'.BRFREE NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.GRAIN '); WRITE(LUPRI,*)'.GRAIN  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.DYNLMA'); WRITE(LUPRI,*)'.DYNLMA NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.ALLSQR'); WRITE(LUPRI,*)'.ALLSQR NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.TRSRCH'); WRITE(LUPRI,*)'.TRSRCH NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOBOXP'); WRITE(LUPRI,*)'.NOBOXP NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.PRSTAT'); WRITE(LUPRI,*)'.PRSTAT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.TCON  '); WRITE(LUPRI,*)'.TCON   NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.SCREEN'); READ(LUCMD,*) DALTON%MM_SCREEN
         DALTON%MM_SCREEN = DALTON%MM_SCREEN/DALTON%THRESHOLD
      CASE ('.SKIPNN'); WRITE(LUPRI,*)'.SKIPNN NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.MMSAVE'); WRITE(LUPRI,*)'.MMSAVE NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOONE '); DALTON%MM_NO_ONE = .TRUE.
      CASE ('.NOMMBU'); DALTON%USEBUFMM  = .FALSE.
      CASE DEFAULT
                  WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                  & '" not recognized in FMM INTEGRAL.'
                 CALL lsQUIT('Illegal keyword in FMM INTEGRAL.',lupri)
    END SELECT
  ENDIF
  IF (PROMPT .EQ. '*') EXIT
ENDDO
END SUBROUTINE READ_INTEGRALS_FMM_INPUT

SUBROUTINE READ_INTEGRALS_FCK3_INPUT(LUPRI,LUCMD,DALTON,word)
implicit none
TYPE(integralconfig)   :: DALTON
INTEGER            :: LUCMD !Logical unit number for the daltoninput
character(len=80),intent(out)  :: WORD
character(len=1)   :: PROMPT
INTEGER            :: LUPRI

! the old version is SUBROUTINE FCK3INP

DO
  READ (LUCMD, '(A40)') WORD
  PROMPT = WORD(1:1)
  IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') CYCLE
  IF (PROMPT .EQ. '.') THEN
    SELECT CASE(WORD) 
      CASE ('.FIXTHR'); WRITE(LUPRI,*)'.FIXTHR NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.THRFCK'); WRITE(LUPRI,*)'.THRFCK NOT IMPLEMENTED IN NEW DRIVER'
         READ(LUCMD, '(A40)') WORD
         WRITE(LUPRI,*)'next word',WORD
      CASE ('.PRINT'); WRITE(LUPRI,*)'.PRINT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOSPEX'); WRITE(LUPRI,*)'.NOSPEX NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOSCRP'); WRITE(LUPRI,*)'.NOSCRP NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOSKCL'); WRITE(LUPRI,*)'.NOSKCL NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOSKCS'); WRITE(LUPRI,*)'.NOSKCS NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.SCRFBT'); WRITE(LUPRI,*)'.SCRFBT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOPACK'); WRITE(LUPRI,*)'.NOPACK NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NOSKAL'); WRITE(LUPRI,*)'.NOSKAL NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.FULL_J'); WRITE(LUPRI,*)'.FULL_J  NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.UNCONT'); WRITE(LUPRI,*)'.UNCONT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.BRFREE'); WRITE(LUPRI,*)'.BRFREE NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.CHCKNA'); WRITE(LUPRI,*)'.CHCKNA NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NCSODB'); WRITE(LUPRI,*)'.NCSODB NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NCSTUV'); WRITE(LUPRI,*)'.NCSTUV NOT IMPLEMENTED IN NEW DRIVER' 

      CASE ('.TWOGRA'); WRITE(LUPRI,*)'.TWOGRA NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.HERMIT'); WRITE(LUPRI,*)'.HERMIT NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.TWOSHI'); WRITE(LUPRI,*)'.TWOSHI NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.PRPSCR'); WRITE(LUPRI,*)'.PRPSCR NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.PRPSTH'); WRITE(LUPRI,*)'.PRPSTH NOT IMPLEMENTED IN NEW DRIVER' 
      CASE ('.NUCATR'); WRITE(LUPRI,*)'.NUCATR NOT IMPLEMENTED IN NEW DRIVER'
      CASE DEFAULT
                  WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,&
                  & '" not recognized in *FCK INTEGRAL.'
                 CALL lsQUIT('Illegal keyword in FCK3 INTEGRAL.',lupri)
    END SELECT
  ENDIF
  IF (PROMPT .EQ. '*') EXIT
ENDDO

END SUBROUTINE READ_INTEGRALS_FCK3_INPUT

!> \brief Check if keywords conform and print configuration.
!> \author S. Host
!> \date March 2010
!>
!> If keywords specified in LSDALTON.INP do not conform, there are two options: \n
!> 1. Clean up, i.e. change the settings specified by the user to something
!>    meaningful. Remember to clarify this in output! E.g. \n
!>    'H1DIAG does not work well with only few saved microvectors for ARH. 
!>    Resetting max. size of subspace in ARH linear equations to', <something meaningful>. \n
!> 2. Quit, if there is no logical way to recover. \n
!> After deciding what the final configuration should be, print selected details
!> which may be useful for the user.
!>
subroutine set_final_config_and_print(lupri,config,ls,nbast)
implicit none
   !> Contains info, settings and data for entire calculation
   type(configItem), intent(inout) :: config
   !> Logical unit number for LSDALTON.OUT
   integer, intent(in)             :: lupri,nbast
   !> Object containing integral settings and molecule
   type(lsitem), intent(inout)     :: ls
!
   integer                         :: i
!   integer                         :: omp_get_num_threads
   logical                         :: file_exists,CABS_BASIS_PRESENT,Phantom
   real(realk)                     :: conv_factor, potnuc, cutoff,inverse_std_conv_factor
   CHARACTER*24, PARAMETER :: AVG_NAMES(5) = &
        &  (/ 'None                    ', &
        &     '                        ', &
        &     'DIIS                    ', &
        &     'EDIIS                   ', &
        &     'Van Lenthe modified DIIS' /)
   CHARACTER*35, PARAMETER :: F2D_NAMES(4) = (/ &
        & 'Diagonalization            ',&
        & 'Direct density optimization',&
        & 'Purification scheme        ',&
        & 'Augmented RH optimization  ' /)
   CHARACTER*35, PARAMETER :: SHIFT_NAMES(5) = (/ &
        & 'Level shifting by MO overlap       ',&
        & 'Level shifting by line search in mu',&
        & 'Level shifting by ||Dorth|| ratio  ',&
        & 'No level shifting                  ',&
        & 'Van Lenthe fixed level shifts      '/)
   integer :: nocc,nvirt,nthreads_test,nthreads
   logical :: dogga,dometa
#ifdef VAR_OMP
integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
integer, external :: OMP_GET_NESTED
!integer, external :: OMP_SET_NESTED
#endif

   write(config%lupri,*) 'Configuration:'
   write(config%lupri,*) '=============='
nthreads = 1
#ifdef VAR_OMP
!deactivates nested OpenMP behavior - this should be false per default-
!but we want to be sure.
!IF(OMP_GET_NESTED())THEN
!   CALL LSQUIT('Nested OpenMP is set to true, deactivate using OMP_NESTED=FALSE',-1)
!ENDIF
!$OMP PARALLEL SHARED(nthreads)  
!$OMP MASTER
nthreads = OMP_GET_NUM_THREADS()
IF(nthreads.GT. 1)THEN
   WRITE(lupri,'(4X,A,I3,A)')'This is an OpenMP calculation using ',nthreads,' threads.'
ELSEIF(nthreads.EQ. 1)THEN
   WRITE(lupri,'(4X,A)')'This is a Single core calculation. (no OpenMP)'
ENDIF
!$OMP END MASTER
!$OMP END PARALLEL
#else
   WRITE(lupri,'(4X,A)')'This is a Single core calculation. (no OpenMP)'
#endif

#ifdef VAR_MPI

if(mod(SPLIT_MPI_MSG,8)/=0)call lsquit("INPUT ERROR: MAX_MPI_MSG_SIZE_NEL has to be a multiple of 8",-1)


IF(nthreads.EQ.1)THEN
 IF(infpar%nodtot.GT.1)THEN
  WRITE(lupri,'(4X,A,I3,A)')'WARNING: This is a MPI calculation using ',infpar%nodtot, &
                          & ' processes, but you are only using 1 OpenMP thread'
  WRITE(lupri,'(4X,A)')     'WARNING: This is NOT recommended! LSDALTON is designed as a MPI/OpenMP hybrid code'
  WRITE(lupri,'(4X,A)')     'WARNING: It is therefore HIGHLY recommended to use the command'
  WRITE(lupri,'(4X,A)')     'WARNING: export OMP_NUM_THREADS=X'
  WRITE(lupri,'(4X,A)')     'WARNING: Where X is the number of floating point cores on your system.'
  WRITE(lupri,'(4X,A)')     'WARNING: Note due to hyper-threading, and vendors reporting number of integer cores instead of' 
  WRITE(lupri,'(4X,A)')     'WARNING: floating point units/cores, the determination of the optimal X may require some testing'
  WRITE(lupri,'(4X,A)')     'WARNING: and may easily be half what you expect'
 ENDIF
ENDIF
#ifdef VAR_INT64
WRITE(lupri,'(4X,A,I3,A)')'This is an 64 bit integer MPI calculation using ',infpar%nodtot,' processes'
#ifdef VAR_MPI_32BIT_INT
!int64,mpi & mpi32
WRITE(lupri,'(4X,A)')'linked to a 32 bit integer MPI library.'
#else
!int64,mpi nompi32
WRITE(lupri,'(4X,A)')'linked to a 64 bit integer MPI library.'
#endif
#else
!int32 mpi
WRITE(lupri,'(4X,A,I3,A)')'This is an MPI calculation using ',infpar%nodtot,' processes'
#endif
#else
!no MPI
WRITE(lupri,'(4X,A)')'This is a serial calculation (no MPI)'
#endif

config%av%diis_history_size = config%av%cfg_settings(config%av%CFG_SET_type)%max_history_size
config%av%ediis_history_size = config%av%diis_history_size
!Printing the configuration for the calculation:
!===============================================

   IF(.NOT.USEXCFUN)THEN
      if(ls%input%DO_DFT) CALL DFTREPORT(lupri) !print the functional
   ELSE
      if(ls%input%DO_DFT) CALL XCFUNDFTREPORT(lupri) !print the functional
   ENDIF
   if(ls%input%DO_DFT)THEN
      IF(config%integral%DFT%CS00)THEN
         WRITE(LUPRI,'(2X,A)')' '
         IF(ABS(config%integral%DFT%CS00SHIFT).LT.1.0E-12)THEN
            WRITE(LUPRI,'(2X,A)')'Using the Casida-Salahub correction with a Zhan-Nichols-Dixon shift'
         ELSE
            WRITE(LUPRI,'(2X,A,D14.6,A)')'Using the Casida-Salahub correction with a shift of ',&
                 & config%integral%DFT%CS00SHIFT,' a.u.'
         ENDIF
      ENDIF
   endif
   if(ls%input%DO_DFT)then
      WRITE(LUPRI,'(2X,A)')' '
      WRITE(LUPRI,'(2X,A)')'The Exchange-Correlation Grid specifications:'
      SELECT CASE(config%integral%DFT%RADIALGRID)
      case(1)
         WRITE(LUPRI,'(2X,A)')'Radial Quadrature : GC2 scheme'
         WRITE(LUPRI,'(2X,A)')'                    Gauss-Chebyshev quadrature of second kind that we' 
         WRITE(LUPRI,'(2X,A)')'                    use to generate (GC2)'
      case(2)
         WRITE(LUPRI,'(2X,A)')'Radial Quadrature : LMG scheme'
         WRITE(LUPRI,'(2X,A)')'                    Radial grid as proposed by R. Lindh, P.A. Malmqvist'
         WRITE(LUPRI,'(2X,A)')'                    and L. Gagliardi. Theor. Chem. Acc. (2001) 106, 178'
      case(3)
         WRITE(LUPRI,'(2X,A)')'Radial Quadrature : Treutler-Ahlrichs M4-T2 scheme'
         WRITE(LUPRI,'(2X,A)')'                    (J. Chem. Phys. (1995) 102, 346).'
         WRITE(LUPRI,'(2X,A)')'                    Implies also that the angular integration &
              &quality becomes Z-dependant'
      CASE DEFAULT
         WRITE (LUPRI,*)'RADIAL GRID TYPE: ',config%integral%DFT%RADIALGRID,' not recognized'
         CALL lsQUIT('Illegal value of DFT%RADIALGRID',lupri)
      END SELECT
      SELECT CASE(config%integral%DFT%PARTITIONING)
      case(1) 
         WRITE(LUPRI,'(2X,A)')'Space partitioning: Stratmann-Scuseria-Frisch partitioning scheme'
         WRITE(LUPRI,'(2X,A)')'                    Chem. Phys. Lett. (1996) vol. 213 page 257'
      case(2)
         WRITE(LUPRI,'(2X,A)')'Space partitioning: Becke partitioning scheme with atomic size correction'
         WRITE(LUPRI,'(2X,A)')'                    J. Chem. Phys. (1988) vol 88 page 2547'
      case(3)
         WRITE(LUPRI,'(2X,A)')'Space partitioning: Becke partitioning scheme without atomic size correction'
         WRITE(LUPRI,'(2X,A)')'                    J. Chem. Phys. (1988) vol 88 page 2547'
      case(4)
         WRITE(LUPRI,'(2X,A)')'Space partitioning: Becke partitioning scheme with atomic size correction'
         WRITE(LUPRI,'(2X,A)')'                    J. Chem. Phys. (1988) vol 88 page 2547'
         WRITE(LUPRI,'(2X,A)')'                    Combined with a blockwise handling of grid points'
         WRITE(LUPRI,'(2X,A)')'                    J. Chem. Phys. (2004) vol 121, page 2915.'
         WRITE(LUPRI,'(2X,A)')'                    Useful for large molecules.'
      case(5)
         WRITE(LUPRI,'(2X,A)')'Space partitioning: Stratmann-Scuseria-Frisch partitioning scheme'
         WRITE(LUPRI,'(2X,A)')'                    Chem. Phys. Lett. (1996) vol. 213 page 257,'
         WRITE(LUPRI,'(2X,A)')'                    Combined with a blockwise handling of grid points'
         WRITE(LUPRI,'(2X,A)')'                    J. Chem. Phys. (2004) vol 121, page 2915.'
         WRITE(LUPRI,'(2X,A)')'                    Useful for large molecules.'
      case(6)
         WRITE(LUPRI,'(2X,A)')'Using CARTESIAN partitioning, currently NOT working - not implemented'
         call lsquit('CARTESIAN partitioning, currently NOT working - not implemented',lupri)
      CASE DEFAULT
         WRITE (LUPRI,*)'PARTITIONING TYPE: ',config%integral%DFT%PARTITIONING,' not recognized'
         CALL lsQUIT('Illegal value of DFT%PARTITIONING',lupri)
      END SELECT
      IF(.NOT.config%integral%DFT%NOPRUN)THEN
         WRITE(config%LUPRI,*)
         WRITE (config%LUPRI,'(2X,A)')'We use grid pruning according to Mol. Phys. (1993) vol 78 page 997'
      ENDIF
      WRITE(config%LUPRI,*)
      WRITE (config%LUPRI,'(4X,A,D11.4)') ' DFT LSint Radial integration threshold:  ',config%integral%DFT%RADINT
      WRITE (config%LUPRI,'(4X,A,I3,A,I3,A)')' DFT LSint integration order range     :   [',&
           &config%integral%DFT%ANGMIN,': ',config%integral%DFT%ANGINT,']'
      WRITE (config%LUPRI,'(4X,A,I3)')' Hardness of the partioning function   : ',config%integral%DFT%HRDNES
      WRITE (config%LUPRI,'(4X,A,3D11.2)')&
           & ' DFT LSint screening thresholds        :', config%integral%DFT%DFTHR0, &
           &config%integral%DFT%DFTHRL, config%integral%DFT%DFTHRI 
      WRITE (config%LUPRI,'(4X,A,1D11.2)')&
           & ' Threshold for number of electrons     :', config%integral%DFT%DFTELS 
      WRITE (config%LUPRI,'(4X,A,D11.4)')&
           & ' The Exact Exchange Factor             :', config%integral%exchangeFactor
   endif
   if (config%opt%cfg_density_method == config%opt%cfg_f2d_arh .and. (.not. config%davidSCF%arh_davidson )) then
      WRITE(config%LUPRI,*)
      write (config%lupri,*) '  You have requested Augmented Roothaan-Hall optimization'
      write (config%lupri,*) '  => explicit averaging is turned off!'
      config%av%cfg_averaging = config%av%cfg_avg_none
      !if (config%solver%debug_arh) then
      !   write (config%lupri,*) '  Warning: Debugging Augmented Roothaan-Hall'
      !   write (config%lupri,*) '  => ARH terms turned off!'
      !   config%solver%cfg_arhterms = .false.
      !endif
      WRITE(config%LUPRI,*)
      write(config%lupri,"('Expand trust radius if ratio is larger than:          ',F6.2)") config%solver%cfg_arh_expand_crit
      write(config%lupri,"('Contract trust radius if ratio is smaller than:       ',F6.2)") config%solver%cfg_arh_contract_crit
      write(config%lupri,"('On expansion, trust radius is expanded by a factor    ',F6.2)") config%solver%cfg_arh_expand
      write(config%lupri,"('On contraction, trust radius is contracted by a factor',F6.2)") config%solver%cfg_arh_contract
      if (config%solver%cfg_arh_truncate) then
         WRITE(config%LUPRI,*)
         if (config%opt%cfg_start_guess == 'H1DIAG' .and. config%solver%cfg_arh_microvecs < 5) then
            WRITE(config%LUPRI,*) 'H1DIAG does not work well with only few saved microvectors for ARH...'
            config%solver%cfg_arh_microvecs = 5
            WRITE(config%LUPRI,*) 'Resetting max. size of subspace in ARH linear equations to:', config%solver%cfg_arh_microvecs
         else        
            WRITE(config%LUPRI,*) 'Maximum size of subspace in ARH linear equations:', config%solver%cfg_arh_microvecs
         endif
      endif
   endif

   if (config%solver%cfg_2nd_order_all) then
      WRITE(config%LUPRI,*)
      config%av%CFG_averaging = config%av%CFG_AVG_none 
      write (config%lupri,*) 'You have requested 2nd order optimization => no averaging (no DIIS)!'
      config%solver%set_do_2nd_order = .true.
   endif

   WRITE(config%LUPRI,*)
   WRITE(config%LUPRI,"('Density subspace min. method    : ',A)") AVG_NAMES(config%av%CFG_averaging)
   WRITE(config%LUPRI,"('Density optimization : ',A)") F2D_NAMES(config%opt%CFG_density_method)
   !if (config%opt%CFG_density_method == config%opt%CFG_F2D_ROOTHAAN) then
   !   cfg_rsp_mostart = .true.
   !   cfg_rsp_moprec = .true.
   !endif

!Settings concerning diagonalization and averaging:
!==================================================

   if (config%diag%cfg_lshift .ne. diag_lshift_none) then
     WRITE(config%LUPRI,"('with ',A)") SHIFT_NAMES(config%diag%CFG_lshift)
     if (config%diag%cfg_lshift ==diag_lshift_dorth) then
        WRITE(config%LUPRI,"('  where ',f6.3,' is largest accepted Dorth ratio  grep')") &
             & config%av%cfg_settings(config%av%CFG_SET_TYPE)%max_dorth_ratio
     endif
   endif

   if (config%diag%cfg_custom_shift) then
      write(config%lupri,*)
      write(config%lupri,"(i3,' levelshifts with the values')") config%diag%cfg_nshifts
      do i = 1, config%diag%cfg_nshifts
         write(config%lupri,"(f5.2)") config%diag%cfg_levelshifts(i)
      enddo
      write(config%lupri, "('have been requested.')")
   endif

   !find the maximum number of stored vectors
   config%av%cfg_settings%max_history_size = config%av%diis_history_size
   WRITE(config%LUPRI,*)
   WRITE(config%LUPRI,*) 'Maximum size of Fock/density queue in averaging:', &
      &  config%av%cfg_settings(config%av%CFG_SET_type)%max_history_size

   if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then
      ALLOCATE(config%diag%cfg_levelshifts(100))

      config%diag%cfg_fixed_shift = .true.
      config%diag%cfg_custom_shift = .true.

      if (config%opt%calctype == config%opt%dftcalc) then
         config%diag%cfg_levelshifts = 1.0E0_realk
      else
         config%diag%cfg_levelshifts = 0.3E0_realk
         do i = 1, 5
            config%diag%cfg_levelshifts(i) = 1.0E0_realk
         enddo
      endif
   endif

!Settings for unrestricted:
!==========================

   IF(MOD(config%integral%nelectrons,2) == 0)THEN
      !Even number of electrons
!      config%decomp%nocc = (config%integral%nelectrons - config%integral%molcharge)/2
      !Cecilie change 07/07 2010: Why subtract molcharge?
      ! Shouldn't it just be nelectrons / 2 ??
      ! This way it gives nocc = 3 for four electrons in H2
      config%decomp%nocc = config%integral%nelectrons/2
      config%decomp%nactive = 0

   ELSE
      !Odd number of electrons
      !Stinne change 23/4-2010: why subtract one here???
      !config%decomp%nocc = (config%integral%nelectrons - 1 - config%integral%molcharge)/2
!      config%decomp%nocc = (config%integral%nelectrons - config%integral%molcharge)/2
      !Cecilie change 07/07 2010: Same here
      config%decomp%nocc = config%integral%nelectrons/2
      config%decomp%nactive = 1
   ENDIF
   config%diag%nocc = config%decomp%nocc

   if (config%decomp%alpha_specified .or. config%decomp%beta_specified) then
      config%integral%unres =.TRUE.
      config%decomp%cfg_unres =.TRUE.
      config%diag%cfg_unres =.TRUE.
      config%opt%cfg_unres =.TRUE.
      config%soeoinp%cfg_unres = .true.
      config%diag%nocca = config%decomp%NOCCA
      config%diag%noccb = config%decomp%NOCCb
      !write(lupri,*) 'alpha_specified, beta_specified', alpha_specified, beta_specified
      if (config%decomp%alpha_specified .and. config%decomp%beta_specified) then
         if (config%decomp%nocca + config%decomp%noccb /= 2*config%decomp%nocc + config%decomp%nactive) then
            call lsquit('Nalpha + Nbeta differs from number of electrons!',config%lupri)
         endif
      else if (config%decomp%alpha_specified) then
         config%decomp%NOCCB = 2*config%decomp%NOCC + config%decomp%nactive - config%decomp%NOCCA
      else if (config%decomp%beta_specified) then
         config%decomp%NOCCA = 2*config%decomp%NOCC + config%decomp%nactive - config%decomp%NOCCB
      endif 
      config%integral%nelectrons = config%decomp%NOCCA + config%decomp%NOCCB
      WRITE(config%LUPRI,*)
      write(LUPRI,'(/,1x,a)') '--------------------------'
      write(LUPRI,'(1x,a)')   '<Unrestricted calculation>'
      write(LUPRI,'(1x,a)')   '--------------------------'
      write(LUPRI,'(1x,a,i6)')   'ALPHA spin occupancy =',config%decomp%nocca
      write(LUPRI,'(1x,a,i6,/)') 'BETA  spin occupancy =',config%decomp%noccb
      call mat_select_type(mtype_unres_dense,lupri)
   else IF(config%decomp%nactive /= 0 .or. config%decomp%cfg_unres) THEN
      !unrestricted SCF if Nelec uneven or if cfg_unres=.true.

      config%integral%unres = .true.
      config%decomp%cfg_unres = .true.
      config%diag%cfg_unres = .true.
      config%opt%cfg_unres = .true.
      config%soeoinp%cfg_unres = .true.

      config%decomp%NOCCA = config%decomp%NOCC
      config%decomp%NOCCB = config%decomp%NOCC + config%decomp%nactive

      config%diag%nocca = config%decomp%NOCCA
      config%diag%noccb = config%decomp%NOCCb
      if(config%integral%nelectrons /= 0) then
WRITE(config%LUPRI,*)
write(config%lupri,*) 'WARNING WARNING WARNING spin check commented out!!! /Stinne'

!FIXME: What should be here?? What does spin mean here??? What do we actually
!want to check???
         !if (config%decomp%spin == 0) then
         !   if (MOD(config%integral%nelectrons, 2) /= 0) then
         !      write(config%lupri,*) 'Spin:', config%decomp%spin
         !      write(config%lupri,*) 'Nelectrons:', config%integral%nelectrons
         !      write(config%lupri,*) 'MOD(Nelectrons,2):', MOD(config%integral%nelectrons, 2)
         !      call lsquit('Spin and number of electrons do not conform!')
         !   endif
         !   config%decomp%NOCCA = config%integral%nelectrons/2
         !   config%decomp%NOCCB = config%integral%nelectrons/2
         !else if (config%decomp%spin == 1) then
         !   if (MOD(config%integral%nelectrons, 2) /= 1) then
         !      write(config%lupri,*) 'Spin:', config%decomp%spin
         !      write(config%lupri,*) 'Nelectrons:', config%integral%nelectrons
         !      write(config%lupri,*) 'MOD(Nelectrons,2):', MOD(config%integral%nelectrons, 2)
         !      call lsquit('Spin and number of electrons do not conform!') 
         !   endif
         !   config%decomp%NOCCA = config%integral%nelectrons/2
         !   config%decomp%NOCCB = config%integral%nelectrons - config%decomp%NOCCA
         !else if(config%decomp%spin == 2) then
         !   if (MOD(config%integral%nelectrons,2) /= 0) then
         !      write(config%lupri,*) 'Spin:', config%decomp%spin
         !      write(config%lupri,*) 'Nelectrons:', config%integral%nelectrons
         !      write(config%lupri,*) 'MOD(Nelectrons,2):', MOD(config%integral%nelectrons, 2)
         !      call lsquit('Spin and number of electrons do not conform!') 
         !   endif
         !   config%decomp%NOCCA = config%integral%nelectrons/2 - 1
         !   config%decomp%NOCCB = config%integral%nelectrons - config%decomp%NOCCA
         !else
         !   call lsquit('Only spin = 0, 1, and 2 implemented!') 
         !endif
      else
         config%decomp%NOCCA = config%integral%nelectrons/2 
         config%decomp%NOCCB = config%integral%nelectrons/2
         if (config%decomp%spin == 0) then
            !occupations as above  
         else if (config%decomp%spin == 1) then
            if(MOD(config%integral%nelectrons, 2) /= 1) &
                 & call lsquit('Spin and Nelec do not conform!',config%lupri) 
            config%decomp%NOCCA = config%integral%nelectrons/2
            config%decomp%NOCCB = config%integral%nelectrons - config%decomp%NOCCA
         else if(config%decomp%spin == 2) then
            if (MOD(config%integral%nelectrons, 2) /= 0) &
                 & call lsquit('Spin and Nelec do not conform!',config%lupri) 
            config%decomp%NOCCA = config%integral%nelectrons/2 - 1
            config%decomp%NOCCB = config%integral%nelectrons - config%decomp%NOCCA
         else
            call lsquit('Only spin = 0, 1, and 2 implemented!',config%lupri) 
         endif
      endif
   ENDIF

   !Check XCFUN consistency for unrestricted calculations
   IF (config%decomp%cfg_unres) THEN
      !For dft calculations we require the use of XCFUN as the interface to PSFUN is flawed
      IF (config%opt%calctype.EQ.config%opt%dftcalc) THEN
       call dft_dogga_dometa(dogga,dometa)
       IF (dogga.OR.dometa) THEN
#ifndef VAR_XCFUN
        write(LUPRI,'(1x,a)')   'The current gga or meta-gga calculation is unrestricted.'
        write(LUPRI,'(1x,a)')   'Unrestricted calculations require the use of XCFUN.'
        write(LUPRI,'(1x,a)')   'Recompile the code with XCFUN turned on to run unrestriced calculations.'
        call lsquit('Unrestricted calculations require the use of XCFUN',-1)
#endif
        IF (.NOT.CONFIG%INTEGRAL%DFT%XCFUN) THEN
          write(LUPRI,'(1x,a)')   'The current gga or meta-gga calculation is unrestricted.'
          write(LUPRI,'(1x,a)')   'Unrestricted calculations require the use of XCFUN.'
          write(LUPRI,'(1x,a)')   'Remove keyword .PSFUN under **GENERAL to run unrestriced calculations.'
          CALL LSQUIT('Unrestricted calculations require the use of XCFUN',lupri)
        ENDIF
       ENDIF
      ENDIF
      WRITE(config%LUPRI,*)
      write(LUPRI,'(/,1x,a)') '--------------------------'
      write(LUPRI,'(1x,a)')   '<Unrestricted calculation>'
      write(LUPRI,'(1x,a)')   '--------------------------'
      if(config%decomp%spin == 2) &
           & write(LUPRI,'(1x,a)') 'Spin symmetry = Triplet'
      write(LUPRI,'(1x,a,i6)')   'ALPHA spin occupancy =', config%decomp%nocca
      write(LUPRI,'(1x,a,i6,/)') 'BETA  spin occupancy =', config%decomp%noccb

      !fixme: should be available for other matrix types as well
      call mat_select_type(mtype_unres_dense,lupri)

   ENDIF

!Settings concerning SCF gradient convergence threshold:
!=======================================================

   if (config%opt%cfg_convdyn) then
      !Until someone writes a better XC grid generation code, we have to use different thresholds for HF and DFT
         SELECT CASE(config%opt%cfg_convdyn_type)
         CASE('EXTRM'); conv_factor = 1.0E-8_realk
         CASE('VTIGH'); conv_factor = 1.0E-7_realk
         CASE('TIGHT'); conv_factor = 1.0E-6_realk
         CASE('STAND'); conv_factor = 1.0E-5_realk
         CASE('SLOPP'); conv_factor = 1.0E-4_realk
         CASE DEFAULT
         WRITE(config%LUPRI,*)
         WRITE (LUPRI,'(/,3A,/)') ' Keyword "',config%opt%cfg_convdyn_type,&
              & '" not recognized with .CONVDYN'
         WRITE (LUPRI,*) 'Options are TIGHT, STANDARD, and SLOPPY.'
         CALL lsQUIT('Illegal keyword with .CONVDYN',config%lupri)
         END SELECT
         inverse_std_conv_factor = 1.0E+5_realk
      IF(conv_factor*inverse_std_conv_factor.LT.0.9)THEN
         Write(config%LUPRI,'(A)')' '
         IF(.NOT.DECinfo%HFrestart)THEN
            Write(config%LUPRI,'(A)')' Due to the tightend SCF convergence threshold we also tighten the integral Threshold'
            Write(config%LUPRI,*)'with a factor:',conv_factor*inverse_std_conv_factor
            ls%input%dalton%THRESHOLD = conv_factor*inverse_std_conv_factor*ls%input%dalton%THRESHOLD
            config%integral%THRESHOLD = conv_factor*inverse_std_conv_factor*config%integral%THRESHOLD
            ls%setting%scheme%THRESHOLD = conv_factor*inverse_std_conv_factor*ls%setting%scheme%THRESHOLD
         ENDIF
      ENDIF
      if (config%decomp%cfg_unres) then
         config%solver%lshift_by_hlgap = .false. !HOMO LUMO shift not implemented for unrestricted
         config%opt%cfg_convergence_threshold = conv_factor*sqrt((config%decomp%nocca+config%decomp%noccb)*1.0E0_realk)
         config%opt%set_convergence_threshold = conv_factor*sqrt((config%decomp%nocca+config%decomp%noccb)*1.0E0_realk)
      else
         config%opt%cfg_convergence_threshold = conv_factor*sqrt(config%decomp%nocc*2.0E0_realk)
         config%opt%set_convergence_threshold = conv_factor*sqrt(config%decomp%nocc*2.0E0_realk)
      endif

      WRITE(lupri,*)
      WRITE(config%LUPRI,"('Dynamic convergence threshold for gradient: ',E10.2)") &
           & config%opt%set_convergence_threshold
   else
      WRITE(lupri,*)
      WRITE(config%LUPRI,"('Convergence threshold for gradient: ',E10.2)") &
           & config%opt%set_convergence_threshold
   endif

!Settings for HOMO-LUMO gap, Hessian eigenvalues and rsp starting guess:
!=======================================================================

   if (config%decomp%cfg_startvectors) then
      if (config%decomp%cfg_no_of_startvectors < config%decomp%cfg_rsp_nexcit) then
         WRITE(config%LUPRI,*)
         WRITE(config%LUPRI,"('WARNING: Number of start vectors smaller than number of excitation energies!')")
         WRITE(config%LUPRI,"('Resetting .NSTART to match .NEXCI')")
         config%decomp%cfg_no_of_startvectors = config%decomp%cfg_rsp_nexcit
      else if (config%decomp%cfg_no_of_startvectors < config%decomp%cfg_hessian_nvec) then
         WRITE(config%LUPRI,*)
         WRITE(config%LUPRI,"('WARNING: Number of start vectors smaller than number of requested Hessian eigenvalues!')")
         WRITE(config%LUPRI,"('Resetting .NSTART to match .HESVEC')")
         config%decomp%cfg_no_of_startvectors = config%decomp%cfg_hessian_nvec
      endif
      if (config%decomp%cfg_homolumo_maxit < 100*config%decomp%cfg_no_of_startvectors) then
         WRITE(config%LUPRI,*)
         WRITE(config%LUPRI,"('Increasing max. no. of HOMO-LUMO iterations to ', i6)") 100*config%decomp%cfg_no_of_startvectors
         config%decomp%cfg_homolumo_maxit = 100*config%decomp%cfg_no_of_startvectors
      endif
      if (config%decomp%cfg_check_maxit < 50*config%decomp%cfg_hessian_nvec) then
         WRITE(config%LUPRI,*)
         WRITE(config%LUPRI,"('Increasing max. no. of Hessian iterations to ', i6)") 50*config%decomp%cfg_hessian_nvec
         config%decomp%cfg_check_maxit = 50*config%decomp%cfg_hessian_nvec
      endif
   else if (config%decomp%cfg_rsp_nexcit > 0) then
      if (config%decomp%cfg_homolumo_maxit < 100*config%decomp%cfg_rsp_nexcit) then
         WRITE(config%LUPRI,*)
         WRITE(config%LUPRI,"('Increasing max. no. of HOMO-LUMO iterations to ', i6)") 100*config%decomp%cfg_rsp_nexcit
         config%decomp%cfg_homolumo_maxit = 100*config%decomp%cfg_rsp_nexcit
      endif
   else if (config%decomp%cfg_hessian_nvec > 1) then
      if (config%decomp%cfg_check_maxit < 50*config%decomp%cfg_hessian_nvec) then
         WRITE(config%LUPRI,*)
         WRITE(config%LUPRI,"('Increasing max. no. of Hessian iterations to ', i6)") 50*config%decomp%cfg_hessian_nvec
         config%decomp%cfg_check_maxit = 50*config%decomp%cfg_hessian_nvec
      endif
   endif

   !Check if HOMO-LUMO is really necessary for further calculations (Hessian eival or exci energies)
   ! - if not, we can print a warning instead of quitting if convergence fails
   if (config%decomp%cfg_check_converged_solution .or. config%decomp%cfg_rsp_nexcit > 0) then
      config%decomp%cfg_hlgap_needed = .true.
   endif
   !Check if (# unoccupied)*(# occupied) is less then the number of excitations asked for
   nocc = config%decomp%nocc
   nvirt = (nbast-nocc)
   if (nocc*nvirt .LT. MAX(config%decomp%cfg_rsp_nexcit,config%response%MCDinput%nexci)) then
      write(lupri,*)'The number of Occupied Orbitals  :',nocc
      write(lupri,*)'The number of Unoccupied Orbitals:',nvirt
      write(lupri,*)'The Maximum number of allowed excitation energies:',nvirt*nocc
      print*,'The Maximum number of allowed excitation energies:',nvirt*nocc,'you asked for',&
           & MAX(config%decomp%cfg_rsp_nexcit,config%response%MCDinput%nexci)
      Call lsquit("Error in Input. The maximim number of excitation energies exceed.",-1)
   endif

! Check for Cartesian basis functions 
!===============================================

   if(.NOT.ls%input%DALTON%DoSpherical .AND. .NOT.(config%opt%cfg_start_guess.EQ.'H1DIAG'))THEN
      WRITE(config%LUPRI,'(A)')'You have specified Cartesian basisfunctions in MOLECULE.INP but at this time'
      WRITE(config%LUPRI,'(A)')'it requires the H1DIAG starting guess.'
      WRITE(config%LUPRI,'(A)')' '
      WRITE(config%LUPRI,'(A)')'.START'
      WRITE(config%LUPRI,'(A)')'H1DIAG'
      WRITE(config%LUPRI,'(A)')' '
      CALL lsQUIT('Cartesian basisfunction without H1DAIG starting guess.',config%lupri)
   endif

! Check Counter Poise Input :
   IF(config%InteractionEnergy.AND.(ls%input%molecule%nSubSystems.NE.2))THEN
      Phantom = .FALSE.
      IF(ls%input%molecule%nSubSystems.EQ.3)THEN
         DO I=1,ls%input%molecule%nAtoms
            IF(ls%input%molecule%Atom(I)%SubSystemIndex.EQ.-2)Phantom = .TRUE.
         ENDDO
      ENDIF
      IF(Phantom)THEN
         WRITE(config%lupri,*)'Molecule contains Phantom atoms (basis functions without nuclei)'
         WRITE(config%lupri,*)'Since not all Phantom atoms have been assigned a SubSystemLabel'
         WRITE(config%lupri,*)'it is assumed that both Subsystems needs these Phantom atoms'          
      ELSE
         call lsquit('.INTERACTIONENERGY keyword require SubSystems in MOLECULE.INP',-1)
      ENDIF
   ENDIF

! Check integral input:
!======================

   if(config%integral%densfit .AND. (.NOT. config%integral%basis(AuxBasParam)))then
      WRITE(config%LUPRI,'(/A)') &
           &     'You have specified .DENSFIT in the dalton input but not supplied a fitting basis set'
      CALL lsQUIT('Density fitting input inconsitensy: add fitting basis set',config%lupri)
   endif
   if(config%doRIMP2 .AND. (.NOT. config%integral%basis(AuxBasParam)))then
      WRITE(config%LUPRI,'(/A)') &
           &     'Warning: You have specified .RIMP2 in the dalton input but not supplied a fitting basis set'
      CALL lsQUIT('Density fitting input inconsitensy: add fitting basis set',config%lupri)
   endif

   CABS_BASIS_PRESENT = config%integral%basis(CABBasParam)

   if(config%doF12)THEN
      if(.NOT. CABS_BASIS_PRESENT)then         
         WRITE(config%LUPRI,'(/A)') &
              &     'You have specified .F12 in the dalton input but not supplied a CABS basis set'
         CALL lsQUIT('F12 input inconsitensy: add CABS basis set',config%lupri)
      endif
      WRITE(config%LUPRI,'(A,F7.1)')'The F12 Geminal exponent for the Regular Basis = ',&
           & ls%setting%basis(1)%p%BINFO(REGBASPARAM)%GeminalScalingFactor
      IF(config%integral%ATOMBASIS)THEN
         WRITE(config%LUPRI,'(A)') 'WARNING BASIS (NOT ATOMBASIS) is required to auto detect the Geminal exponent'
      ENDIF
   endif
!ADMM basis input
   if(config%integral%ADMM_EXCHANGE) THEN
      IF (.NOT. config%integral%basis(ADMBasParam)) THEN
        WRITE(config%LUPRI,'(/A)') &
     &   'You have specified an ADMM calculation in the dalton input but not supplied a '
        WRITE(config%LUPRI,'(A)') &
     &   'ADMM fitting basis set in the molecule input as required'
        WRITE(config%LUPRI,'(/A)') &
             &     'Please read the ADMM part in the manual and supply an ADMM basis set'
        CALL lsQUIT('ADMM fitting input inconsitensy: add ADMM fitting basis set',config%lupri)
      ENDIF
   endif

   if(config%response%tasks%doResponse.AND.(config%integral%pari_J.OR.config%integral%pari_K.OR.config%integral%mopari_K))then
      WRITE(config%LUPRI,'(/A)') 'The Pari keywords do not currently work with response'
      WRITE(config%LUPRI,'(/A)') 'Please remove the Pari keywords'
      CALL lsQUIT('The Pari keywords do not currently work with response',config%lupri)
   endif

   if(config%response%tasks%doResponse.AND.config%integral%FMM)then
      WRITE(config%LUPRI,'(/A)') 'The .RUNMM keyword do not currently work with response'
      WRITE(config%LUPRI,'(/A)') 'Please remove the .RUNMM keyword'
      CALL lsQUIT('The .RUNMM keyword do not currently work with response',config%lupri)
   endif

   if(config%integral%DALINK .AND. config%opt%cfg_incremental)THEN
      WRITE(config%LUPRI,*)
      WRITE(config%lupri,*)'DalinK in combination with incremental'
      WRITE(config%lupri,*)'can be tricky so we increse threshold'
      print*,'DalinK in combination with incremental'
      print*,'can be tricky so we increse threshold'
      config%INTEGRAL%DASCREEN_THRLOG = 3 !tighten by 1000
      ls%input%dalton%DASCREEN_THRLOG = 3
      ls%setting%scheme%DASCREEN_THRLOG = 3
   endif

   if((.NOT.config%integral%LINK).AND.(.NOT.config%integral%JENGINE))THEN
      WRITE(config%LUPRI,*)
      IF(config%opt%cfg_incremental)THEN
         WRITE(config%lupri,*)'You have chosen to run without the use of Jengine or LinK'
         WRITE(config%lupri,*)'We therefore deactivate the incremental Fock matrix scheme'
         WRITE(config%lupri,*)'as it will provide no speedup - only reduced accuracy'
      ENDIF
      config%opt%cfg_incremental = .FALSE.
   endif
   
   !Note that the placement of this if statement is important, as the config%opt%cfg_incremental is subject to change
   if(config%opt%cfg_incremental)THEN
      !incremental and no change in threshold
      WRITE(config%LUPRI,*)
      WRITE(config%lupri,*)'WARNING: You have chosen to run with the incremental Fock matrix scheme'
      WRITE(config%lupri,*)'Due to the accumulated error in the incremental scheme the integral'
      WRITE(config%lupri,*)'We recommend to tighten the threshold, with' 
      WRITE(config%lupri,*)'.THRESH'
      WRITE(config%lupri,*)'1.0D-9'
      WRITE(config%lupri,*)'The scheme is first activated when the maximum element of the '
      WRITE(config%lupri,*)'differences Matrix is below 0.1'
      WRITE(config%lupri,*)'to the LSDALTON.INP file.'
   endif

   IF(config%decomp%cfg_gcbasis)THEN
      IF(config%INTEGRAL%FORCEGCBASIS)THEN
         !do nothing
         IF(ls%input%basis%BINFO(REGBASPARAM)%DunningsBasis)THEN
            WRITE(config%lupri,*)'We have detected a Dunnings Basis but the ' 
            WRITE(config%lupri,*)'FORCEGCBASIS keyword is in effect.'
         ENDIF
      ELSEIF(ls%input%basis%BINFO(REGBASPARAM)%DunningsBasis)THEN
         WRITE(config%lupri,*)'We have detected a Dunnings Basis set so we deactivate the' 
         WRITE(config%lupri,*)'use of the Grand Canonical basis, which is normally default.'
         WRITE(config%lupri,*)'The use of Grand Canonical basis can be enforced using the FORCEGCBASIS keyword' 
         config%decomp%cfg_gcbasis = .FALSE.         
      ENDIF
   ENDIF

   IF(config%integral%contang) THEN
      WRITE(config%lupri,'(A)')' '
      WRITE(config%lupri,'(A)')'Shell ordering has been set to contracted first and angular components second.'
      WRITE(config%lupri,'(A)')'To use the default angular contracted ordering remove keyword .CONTANG'
   ENDIF

   config%INTEGRAL%NOGCINTEGRALTRANSFORM = ls%input%dalton%NOGCINTEGRALTRANSFORM
   IF(config%decomp%cfg_gcbasis)THEN
      WRITE(config%lupri,'(A)')' '
      WRITE(config%lupri,'(A)')'We perform the calculation in the Grand Canonical basis'
      WRITE(config%lupri,'(A)')'(see PCCP 2009, 11, 5805-5813)'
      WRITE(config%lupri,'(A)')'To use the standard input basis use .NOGCBASIS'
      
      IF(ls%input%dalton%NOGCINTEGRALTRANSFORM)THEN
         WRITE(config%lupri,'(A)')' '
         WRITE(lupri,'(A)')'Due to the presence of the keyword (default for correlation)'
         WRITE(lupri,'(A)')'.NOGCINTEGRALTRANSFORM'
         WRITE(lupri,'(A)')'We transform the input basis to the Grand Canonical' 
         WRITE(lupri,'(A)')'basis and perform integral evaluation using this basis' 
         ls%setting%integraltransformGC = .FALSE.
      ELSE!default
         WRITE(config%lupri,'(A)')' '
         IF(.NOT.ls%input%basis%BINFO(REGBASPARAM)%Gcont)THEN
            WRITE(lupri,'(A)')'Since the input basis set is a segmented contracted basis we'
            WRITE(lupri,'(A)')'perform the integral evaluation in the more efficient'
            WRITE(lupri,'(A)')'standard input basis and then transform to the Grand '
            WRITE(lupri,'(A)')'Canonical basis, which is general contracted.'
            WRITE(lupri,'(A)')'You can force the integral evaluation in Grand '
            WRITE(lupri,'(A)')'Canonical basis by using the keyword'
            WRITE(lupri,'(A)')'.NOGCINTEGRALTRANSFORM'
            ls%setting%integraltransformGC = .TRUE.
         ELSE
            WRITE(lupri,'(A)')'Since the input basis set is a general contracted basis like'
            WRITE(lupri,'(A)')'the Grand Canonical basis, we transform the input basis to the'
            WRITE(lupri,'(A)')'Grand Canonical basis and perform integral evaluation using this basis'
            ls%setting%integraltransformGC = .FALSE.
         ENDIF
      ENDIF
   ELSE
      WRITE(config%lupri,'(A)')'We perform the calculation in the standard input basis'      
      ls%setting%integraltransformGC = .FALSE.
   ENDIF



! Write Screening Thresholds:
!======================
   WRITE(config%LUPRI,'(A)')' '
   WRITE(config%LUPRI,'(A60,ES12.4)')'The Overall Screening threshold is set to              :',config%integral%THRESHOLD
   WRITE(config%LUPRI,'(A60,ES12.4)')'The Screening threshold used for Coulomb               :',&
& config%integral%THRESHOLD*config%integral%J_THR
   WRITE(config%LUPRI,'(A60,ES12.4)')'The Screening threshold used for Exchange              :',&
&config%integral%THRESHOLD*config%integral%K_THR
   WRITE(config%LUPRI,'(A60,ES12.4)')'The Screening threshold used for One-electron operators:',&
&config%integral%THRESHOLD*config%integral%ONEEL_THR
   if(config%integral%DALINK)THEN
      WRITE(config%LUPRI,'(A)')' '
      WRITE(config%LUPRI,'(A,ES12.4)')'   DaLink have been activated, so in addition to using ',&
           & config%integral%THRESHOLD*config%integral%K_THR      
      WRITE(config%LUPRI,'(A)')'   as a screening threshold on the integrals contribution to'
      WRITE(config%LUPRI,'(A)')'   the Fock matrix, we also use a screening threshold'
      WRITE(config%LUPRI,'(A,ES12.4)')'   on the integrals contribution to the Energy:       ',&
      &config%integral%THRESHOLD*config%integral%K_THR*(1.0E+1_realk**(-config%INTEGRAL%DASCREEN_THRLOG))
   endif

!Check setting for linear equations iterative solver:
!====================================================

  config%opt%do_trustregion = &
      & (config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
      & (config%solver%cfg_2nd_order_all .or. config%solver%cfg_2nd_order_local))

   !Currently turned off, get_from_modFIFO_disk won't work!
   !if (config%solver%cfg_arh_disk_macro .and. config%opt%cfg_queue_on_disk) then
   !   write(config%lupri,*) 'It does not make sense to use .DUMP (dump previous'
   !   write(config%lupri,*) 'F, D matrices to disk while constructing new F), since'
   !   write(config%lupri,*) 'you have also requested that these should always'
   !   write(config%lupri,*) 'be stored on disk! Ignoring .DUMP !!'
   !   config%opt%cfg_queue_on_disk = .false.
   !endif

   if (config%opt%cfg_queue_on_disk) then
      WRITE(config%LUPRI,*)
      write(config%lupri,*) 'Dump previous F, D matrices to disk while constructing new F'
      write(config%lupri,*) '- memory saving, but time consuming!'
   endif

   !if (config%opt%cfg_arh_disk_macro) then !Currently turned off, get_from_modFIFO_disk won't work!
   !   write(config%lupri,*) 'Keep queue of previous F, D matrices on disk instead of in core'
   !   write(config%lupri,*) '- memory saving, but time consuming!'
   !endif

   if (config%solver%cfg_arh_disk_micro) then
      WRITE(config%LUPRI,*)
      write(config%lupri,*) 'ARH solver: Keep trial and sigma vectors on disk instead of in core'
      write(config%lupri,*) '- memory saving, but time consuming!'
   endif

   if (config%solver%cfg_arh_newdamp .and. .not. config%solver%cfg_arh_truncate) then
      WRITE(config%LUPRI,'(/A)') &
      &     'Can only use new damping scheme with .TRUNCATE'
       CALL lsQUIT('Can only use new damping scheme with .TRUNCATE',config%lupri)
   endif 

   if (config%opt%CFG_density_method == config%opt%CFG_F2D_DIRECT_DENS) then
      config%solver%cfg_arhterms = .false.
      config%solver%set_arhterms = .false.
      config%solver%cfg_arh_crop = .true.
      config%solver%cfg_arh_truncate = .false.
   endif

!Check for stuff not implemented for unrestricted:
!=================================================

   if (config%opt%cfg_start_guess=='TRILEVEL' .and. config%decomp%cfg_unres) then
      call lsquit('Sorry, trilevel starting guess not implemented for unrestricted!',config%lupri)
   endif

   if (config%solver%debug_dd) then
      if (config%decomp%cfg_unres) call lsquit('DEBUG_DD not implemented for unrestricted',config%lupri)
      call scf_stats_arh_header(config%lupri)
   endif

#ifdef HAS_PCMSOLVER
! Cross-check for polarizable continuum model
!=================================================
   if (config%pcm%do_pcm .and. config%optinfo%optimize) then
      call lsquit('Sorry, polarizable continuum model not available for geometry optimizations!',config%lupri)
   endif
   if (config%pcm%do_pcm .and. config%dynamics%do_dynamics) then
      call lsquit('Sorry, polarizable continuum model not available for dynamics!',config%lupri)
   endif
#endif   

!MKL sanity check:
!==================

   if (config%opt%cfg_prefer_CSR) then
      if (matrix_type == mtype_unres_dense) then
         call lsquit('Compressed Sparse Row (CSR) not implemented for unrestricted!',config%lupri)
      else
#ifdef VAR_MKL
         CALL mat_select_type(mtype_csr,lupri)         
         call mat_inquire_cutoff(cutoff)
         write(config%lupri,*)
         write(config%lupri, '("Using Compressed-Sparse Row matrices - zero cutoff is ", d10.2)') cutoff
#else
         call lsquit('.CSR requires MKL library and -DVAR_MKL precompiler flag',config%lupri)
#endif
      endif
      if(config%integral%densfit)then
         call lsquit('currently .CSR do not work in combination with .DENSFIT',-1)
      endif
   endif

   if (matrix_type == mtype_csr) then
      if (config%opt%CFG_density_method == config%opt%CFG_F2D_roothaan) then
         write(config%lupri,*)
         write(config%lupri,*) 'The combination of diagonalization and CSR is very inefficient!'
         write(config%lupri,*) 'Please choose another type of density optimization' 
         write(config%lupri,*) '(e.g. .ARH, .TrFD, 2ND_ALL, 2ND_LOC) and welcome back...' 
         call lsquit('Combining diagonalization and CSR is inefficient!',-1)
      endif
   endif

!OpenMP sanity check: -DVAR_OMP should be set if openMP is active
!==================
#ifndef VAR_OMP
!OpenMP should be turned off
nthreads_test = 0
!$OMP PARALLEL SHARED(nthreads_test)

!$OMP CRITICAL
nthreads_test = nthreads_test + 1
!$OMP END CRITICAL

!$OMP END PARALLEL 
IF(nthreads_test.NE.1)THEN
   print*,'OpenMP compilation inconsistency: use -DVAR_OMP when using -openmp/-fopenmp flag'
   call lsquit('OpenMP compilation inconsistency: use -DVAR_OMP',-1)
ENDIF
#endif

!PDMM sanity check:
!==================

if (config%opt%cfg_prefer_PDMM) then
   if (matrix_type == mtype_unres_dense) then
      call lsquit('PDMM not implemented for unrestricted!',config%lupri)
   else
#ifdef VAR_MPI
      WRITE(lupri,'(4X,A,I3,A)')'This is an MPI calculation using ',infpar%nodtot,' processes combined'
      WRITE(lupri,'(4X,A)')'with PDMM for memory distribution and parallelization.'
      CALL mat_select_type(mtype_pdmm,lupri,nbast)      
#else
      !no VAR_MPI
      WRITE(lupri,'(4X,A)')'This is a Standard Serial compilation.'
      WRITE(lupri,'(4X,A)')'.PDMM requires compilation using MPI.'
      print*,'This is a Standard Serial compilation.'
      print*,'.PDMM requires compilation using MPI.'
      CALL LSQUIT('PDMM requires MPI - recompile using MPI and the -DVAR_MPI flag',config%lupri)
#endif
   endif
endif

!SCALAPACK sanity check:
!==================

   if (config%opt%cfg_prefer_SCALAPACK) then
      if (matrix_type == mtype_unres_dense) then
         call lsquit('SCALAPACK not implemented for unrestricted!',config%lupri)
      else
#ifdef VAR_SCALAPACK
#ifdef VAR_MPI
         WRITE(lupri,'(4X,A,I3,A)')'This is an MPI calculation using ',infpar%nodtot,' processes combined'
         WRITE(lupri,'(4X,A)')'with SCALAPACK for memory distribution and parallelization.'
         CALL mat_select_type(mtype_scalapack,lupri,nbast)

#ifdef VAR_INT64
#ifdef VAR_MPI_32BIT_INT
         print*,'you cannot compile using a 64 bit integers, when linking to a 32 bit integer library and'
         print*,'use the 64 bit integer BLACS/SCALAPACK provided by MKL/intel'
         write(config%lupri,*)'you cannot compile using a 64 bit integers, when linking to a 32 bit integer library and'
         write(config%lupri,*)'use the 64 bit integer BLACS/SCALAPACK provided by MKL/intel'
         call lsquit('you cannot compile with VAR_INT64 and SCALAPACK and VAR_MPI32',-1)
#endif 
#endif

#else
         !VAR_SCALAPACK but no VAR_MPI
         CALL LSQUIT('SCALAPACK requires MPI - recompile using MPI and the -DVAR_MPI flag',config%lupri)
#endif
#else
         !no VAR_SCALAPACK
#ifdef VAR_MPI
         WRITE(lupri,'(4X,A,I3,A)')'This is an MPI calculation using ',infpar%nodtot,' processes.'
         call lsquit('.SCALAPACK requires -DVAR_SCALAPACK precompiler flag',config%lupri)
#else
         !no VAR_SCALAPACK and no MPI         
         WRITE(lupri,'(4X,A)')'This is a Standard Serial compilation.'
         WRITE(lupri,'(4X,A)')'.SCALAPACK requires -DVAR_SCALAPACK precompiler flag and compilation using MPI.'
         print*,'This is a Standard Serial compilation.'
         print*,'.SCALAPACK requires -DVAR_SCALAPACK precompiler flag and compilation using MPI.'
         call lsquit('.SCALAPACK requires -DVAR_SCALAPACK precompiler flag and compilation using MPI',config%lupri)
#endif
#endif
      endif
   else
      IF(config%integral%MEMDIST)call lsquit('.MEMDIST require .SCALAPACK',-1)
      IF(ls%setting%scheme%MEMDIST)call lsquit('.MEMDIST require .SCALAPACK',-1)
   endif

!Response requires the VAR_RSP
   if(config%response%tasks%doResponse) then
#ifdef VAR_RSP
      !everything is fine
#else
      call lsquit('Response Calculations require compilation with -DVAR_RSP',-1)
#endif
   end if

!Local Excited state geometry optimization check:
!================================================

   IF(config%response%tasks%doESG.AND.config%optinfo%optimize)THEN
      print*,'config%doESGopt=.TRUE.'
      WRITE(config%lupri,*)'config%doESGopt=.TRUE.'
      config%doESGopt=.TRUE.
   ENDIF
   ls%input%dalton%unres = config%decomp%cfg_unres

   CALL II_get_nucpot(lupri,lupri,ls%setting,POTNUC)
   config%opt%potnuc = POTNUC
   ls%input%potnuc = POTNUC

   IF(config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens)THEN
      config%opt%cfg_oao_gradnrm = .TRUE.
   ELSEIF(config%opt%cfg_density_method == config%opt%cfg_f2d_roothaan)THEN
      config%opt%cfg_oao_gradnrm = .FALSE.
   ELSEIF(config%opt%cfg_density_method == config%opt%cfg_f2d_arh)THEN
      config%opt%cfg_oao_gradnrm = .TRUE.
   ELSE
      config%opt%cfg_oao_gradnrm = .TRUE.
   ENDIF
   IF(config%opt%cfg_oao_gradnrm)THEN
      Write(config%lupri,*)'The SCF Convergence Criteria is applied to the gradnorm in OAO basis'
   ELSE
      Write(config%lupri,*)'The SCF Convergence Criteria is applied to the gradnorm in AO basis'
   ENDIF

   ! ORBITAL LOCALIZATION RELATED MATTERS
   if (config%decomp%cfg_mlo .and. (.not. config%decomp%cfg_lcm )) then
       write(config%lupri,*)  
       write(config%lupri,*) '    ** ORBITAL LOCALIZATION INFORMATION **'  
       write(config%lupri,*) ' Orbital localization using canonical molecular orbitals as starting guess '
       write(config%lupri,*) ' is requested. For more efficient orbital localization (particularly for large) '
       write(config%lupri,*) ' system, use options  .START/TRILEVEL and .LCM in *DENSOPT section.   '  
       write(config%lupri,*) 
   endif 
   ! Check that DEC input is consistent with geometry optimization and orbital localization.
   call DEC_meaningful_input(config)

   if(config%doDEC) then
      nocc = config%decomp%nocc
      nvirt = (nbast-nocc)
      call check_cc_input(ls,nocc,nvirt,nbast)
   endif
   write(config%lupri,*)
   write(config%lupri,*) 'End of configuration!'
   write(config%lupri,*)
   ls%setting%GPUMAXMEM = config%GPUMAXMEM
end subroutine set_final_config_and_print

SUBROUTINE TRIM_STRING(string,n,words)
  implicit none
  character(len=*), intent(in) :: string
  integer, intent(inout) :: n  !max number of wanted words on input
                               !number of found words on output
  character(len=9), intent(out) :: words(n)
  integer :: i,imax,j,k,nspaces
  logical :: inword

  do i=1,n
    words(i) = '        '
  enddo
  inword = .false.
  imax = LEN(string)
  nspaces = 0
  i = 0  !string index
  j = 0  !number of words found
  do
    i = i + 1
    if (i > imax .or. nspaces == 10) then
!      print*, 'Finished searching the string', &
!                   & 'found ',j,' words out of ',n
      n = j
      exit
    endif
    if (string(i:i) == ' ') then
      inword = .false.
      nspaces = nspaces + 1
      if (j == n) then
!        WRITE(LUPRI,*) (words(j),j=1,n)
        exit  !found all words
      endif
      cycle
    else
      nspaces = 0
      if (.not. inword) then
        !found beginning of word
        j = j + 1
        k = 1 !word index
        inword = .true.
      else
        k = k + 1
        if (k > 9) cycle  !cannot take more that 8 char pr. word
      endif
      words(j)(k:k) = string(i:i)
    endif
  enddo

END SUBROUTINE TRIM_STRING

end module configuration

#ifdef VAR_MPI
subroutine lsmpi_setmasterToSlaveFunc(WORD,hfweight)
use infpar_module
use xcfun_host,only: USEXCFUN
use lsmpi_mod
  implicit none
  character(len=1024)  :: WORD
  real(realk) :: hfweight
  call ls_mpibcast(WORD,1024,infpar%master,MPI_COMM_LSDALTON)
  call ls_mpibcast(hfweight,infpar%master,MPI_COMM_LSDALTON)
  call ls_mpibcast(USEXCFUN,infpar%master,MPI_COMM_LSDALTON)
end subroutine lsmpi_setmasterToSlaveFunc

subroutine lsmpi_setSlaveFunc()
use infpar_module
use lsmpi_mod
use xcfun_host,only: xcfun_host_init, USEXCFUN
use IIDFTINT, only: II_DFTsetFunc
use typedef
  implicit none
  character(len=1024)  :: WORD
  real(realk) :: hfweight
  integer :: ierror
  call ls_mpibcast(WORD,1024,infpar%master,MPI_COMM_LSDALTON)
  call ls_mpibcast(hfweight,infpar%master,MPI_COMM_LSDALTON)
  call ls_mpibcast(USEXCFUN,infpar%master,MPI_COMM_LSDALTON)
  call II_DFTsetFunc(WORD,hfweight,6)
end subroutine lsmpi_setSlaveFunc

subroutine lsmpi_addSlaveFunc()
use infpar_module
use lsmpi_mod
use xcfun_host,only: USEXCFUN
use IIDFTINT, only: II_DFTaddFunc
use typedef
  implicit none
  character(len=1024)  :: WORD
  real(realk) :: hfweight
  call ls_mpibcast(WORD,1024,infpar%master,MPI_COMM_LSDALTON)
  call ls_mpibcast(hfweight,infpar%master,MPI_COMM_LSDALTON)
  call ls_mpibcast(USEXCFUN,infpar%master,MPI_COMM_LSDALTON)
  call II_DFTaddFunc(WORD,hfweight,6)
end subroutine lsmpi_addSlaveFunc
#endif
