!> @file 
!> Settings info and read of input for DEC/Dalton tests (DEC=Divide-Expand-Consolidate Coupled cluster)

!> \brief Settings info for simulated DEC tests.
!> \author \latexonly Kasper Kristensen  \endlatexonly
!> \date 2010-06-16
!>
MODULE DEC_settings_mod
  use typedeftype
  use fundamental
  use precision
  use memory_handling
  use dec_typedef_module
  use dec_fragment_utils
  use ls_util
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type, only: LSMPIASYNCP
#endif

contains

  !> \brief Set default DEC settings.
  !> See explanation of parameters in type DEC_settings.
  !> \author Kasper Kristensen
  !> \date June 2010
  subroutine dec_set_default_config(output)

    implicit none
    !> Unit number for DALTON.OUT
    integer, intent(in) :: output

    DECinfo%doDEC                  = .false.

    ! Orbital-based DEC scheme 
    DECinfo%DECCO = .false.

    ! Max memory measured in GB. By default set to 2 GB
    DECinfo%memory                 = 2.0E0_realk
    DECinfo%memory_defined         = .false.
    DECinfo%use_system_memory_info = .false.

    ! -- Type of calculation
    DECinfo%full_molecular_cc = .false. ! full molecular cc
    DECinfo%print_frags = .false.
    DECinfo%simulate_full     = .false.
    DECinfo%simulate_natoms   = 1
    DECinfo%SkipReadIn        = .false.
    DECinfo%SinglesPolari     = .false.
    DECinfo%SinglesThr        = 0.2E0_realk   ! this is completely random, currently under investigation
    DECinfo%convert64to32     = .false.
    DECinfo%convert32to64     = .false.
    DECinfo%HFrestart         = .false.
    DECinfo%DECrestart        = .false.
    DECinfo%TimeBackup        = 300.0E0_realk   ! backup every 5th minute
    DECinfo%read_dec_orbitals = .false.
    DECinfo%CheckPairs        = .false.
    DECinfo%frozencore        = .false.
    DECinfo%ncalc             = 0
    call dec_set_model_names(DECinfo)


    ! -- Debug modes
    DECinfo%CRASHCALC            = .false.
    DECinfo%cc_driver_debug      = .false.
    DECinfo%CCDEBUG              = .false.
    DECinfo%manual_batchsizes    = .false.
    DECinfo%ccsdAbatch           = 0
    DECinfo%ccsdGbatch           = 0
    DECinfo%hack                 = .false.
    DECinfo%hack2                = .false.
    DECinfo%mpisplit             = 10
#ifdef VAR_MPI
    DECinfo%dyn_load             = LSMPIASYNCP
#endif
    DECinfo%force_scheme         = .false.
    DECinfo%en_mem               = 0
    DECinfo%tensor_test           = .false.
    DECinfo%reorder_test         = .false.
    DECinfo%CCSDno_restart       = .false.
    DECinfo%CCSDnosaferun        = .false.
    DECinfo%solver_par           = .false.
    DECinfo%CCSDpreventcanonical = .false.
    DECinfo%CCSD_NO_DEBUG_COMM   = .true.
    DECinfo%spawn_comm_proc      = .false.
    DECinfo%CCSDmultipliers      = .false.
    DECinfo%use_pnos             = .false.
    DECinfo%pno_S_on_the_fly     = .false.
    DECinfo%noPNOtrafo           = .false.
    DECinfo%noPNOtrunc           = .false.
    DECinfo%simplePNOthr         = 1.0E-7_realk
    DECinfo%EOSPNOthr            = 1.0E-5_realk
    DECinfo%noFAtrunc            = .false.
    DECinfo%noFAtrafo            = .false.
    DECinfo%noPNOoverlaptrunc    = .false.
    DECinfo%PNOoverlapthr        = 1.0E-5_realk
    DECinfo%PNOtriangular        = .true.
    DECinfo%CCDhack              = .false.
    DECinfo%NO_MO_CCSD           = .false.
    DECinfo%v2o2_free_solver     = .false.

    ! -- Output options 
    DECinfo%output               = output


    ! -- Orbital
    DECinfo%mulliken_threshold           = 0.01
    DECinfo%simple_mulliken_threshold    = .false.
    DECinfo%approximated_norm_threshold  = 0.1E0_realk
    DECinfo%check_lcm_orbitals           = .false.
    DECinfo%check_Occ_SubSystemLocality  = .false.
    DECinfo%force_Occ_SubSystemLocality  = .false.
    DECinfo%use_canonical                = .false.
    DECinfo%AbsorbHatoms                 = .true.  ! reassign H atoms to heavy atom neighbour
    DECinfo%mulliken                     = .false.
    DECinfo%Distance                     = .false.
    DECinfo%BoughtonPulay                = .false.
    DECinfo%FitOrbitals                  = .true.
    DECinfo%simple_orbital_threshold     = 0.05E0_realk


    ! -- Fragment
    DECinfo%MaxIter                = 20
    DECinfo%FOTlevel               = 1
    DECinfo%FOT                    = 1.0E-4_realk
    DECinfo%InclFullMolecule       = .false.
    DECinfo%PL                     = 0
    DECinfo%PurifyMOs              = .false.
    DECinfo%precondition_with_full = .false.
    DECinfo%Frag_Exp_Scheme        = 1
    DECinfo%Frag_RedOcc_Scheme        = 1
    DECinfo%Frag_RedVir_Scheme        = 1
    DECinfo%Frag_Init_Size         = 4
    DECinfo%Frag_Exp_Size          = 10
    DECinfo%frag_red_occ_thr       = 1.0  ! times FOT
    DECinfo%frag_red_virt_thr      = 0.9  ! times FOT
    DECinfo%fragadapt              = .false.
    DECinfo%only_n_frag_jobs       =  0
    DECinfo%frag_job_nr            => null()
    DECinfo%only_pair_frag_jobs    =  .false.
    ! for CC models beyond MP2 (e.g. CCSD), option to use MP2 optimized fragments
    DECinfo%fragopt_exp_model      = MODEL_MP2  ! Use MP2 fragments for expansion procedure by default
    DECinfo%fragopt_red_model      = MODEL_MP2  ! Use MP2 fragments for reduction procedure by default
    DECinfo%no_orb_based_fragopt   = .false.
    DECinfo%OnlyOccPart            = .false.
    DECinfo%OnlyVirtPart           = .false.
    ! Repeat atomic fragment calcs after fragment optimization
    DECinfo%RepeatAF               = .true.
    ! Which scheme to used for generating correlation density defining fragment-adapted orbitals
    DECinfo%CorrDensScheme         = 1

    ! -- Pair fragments
    DECinfo%pair_distance_threshold = 1000.0E0_realk/bohr_to_angstrom
    DECinfo%PairMinDist             = 3.0E0_realk/bohr_to_angstrom  ! 3 Angstrom
    DECinfo%pairFOthr               =  0.0_realk
    DECinfo%PairMP2                 = .false.
    DECinfo%PairEstimate            = .true.
    DECinfo%PairEstimateIgnore      = .false.
    DECinfo%EstimateINITradius      = 2.0E0_realk/bohr_to_angstrom
    DECinfo%EstimateInitAtom        = 1

    ! Memory use for full molecule structure
    DECinfo%fullmolecule_memory     = 0E0_realk


    ! -- CC solver options

    DECinfo%ccsd_expl               = .false.
    DECinfo%ccMaxIter               = 100
    DECinfo%ccMaxDIIS               = 3
    DECinfo%ccModel                 = MODEL_MP2 ! see parameter-list in dec_typedef.f90
    DECinfo%F12                     = .false.
    DECinfo%F12debug                = .false.
    DECinfo%SOS                     = .false.
    DECinfo%PureHydrogenDebug       = .false.
    DECinfo%StressTest              = .false.
    DECinfo%DFTreference            = .false.
    DECinfo%ccConvergenceThreshold  = 1e-5_realk
    DECinfo%CCthrSpecified          = .false.
    DECinfo%use_singles             = .false.
    DECinfo%use_preconditioner      = .true.
    DECinfo%use_preconditioner_in_b = .true.
    DECinfo%use_crop                = .true.
    DECinfo%array4OnFile            = .false.
    DECinfo%array4OnFile_specified  = .false.

    ! ccsd(t) settings
    DECinfo%abc = .false.
    DECinfo%abc_tile_size = 1

    ! First order properties
    DECinfo%first_order = .false.

    !> MP2 density matrix   
    DECinfo%density = .false.
    DECinfo%SkipFull = .false.

    !-- MP2 gradient
    DECinfo%gradient=.false.
    DECinfo%kappa_use_preconditioner=.true.
    DECinfo%kappa_use_preconditioner_in_b=.true.
    DECinfo%kappaMaxDIIS=3
    DECinfo%kappaMaxIter=100
    DECinfo%kappa_driver_debug=.false.
    DECinfo%kappaTHR=1e-4_realk
    DECinfo%EerrFactor = 0.1_realk
    DECinfo%EerrOLD = 0.0_realk

    ! -- Timings

    !> MPI (undefined by default)
    DECinfo%MPIgroupsize=0

    ! Test stuff



  end subroutine dec_set_default_config

  !> \brief Set names for models in DEC
  subroutine dec_set_model_names(DECitem)
    implicit none
    !> The DEC item
    type(decsettings),intent(inout) :: DECitem

    DECitem%cc_models(MODEL_MP2)   ='MP2     '
    DECitem%cc_models(MODEL_CC2)   ='CC2     '
    DECitem%cc_models(MODEL_CCSD)  ='CCSD    '
    DECitem%cc_models(MODEL_CCSDpT)='CCSD(T) '
    DECitem%cc_models(MODEL_RPA)   ='RPA     '
    DECitem%cc_models(MODEL_RIMP2) ='RIMP2   '

  end subroutine dec_set_model_names


  !> \brief Read the **DEC or **CC input section in LSDALTON.INP and set 
  !> configuration structure accordingly.
  !> \author Kasper Kristensen
  !> \date September 2010
  SUBROUTINE config_dec_input(input,output,readword,word,fullcalc,doF12)
    implicit none
    !> Logical for keeping track of when to read
    LOGICAL,intent(inout)                :: READWORD
    !> Logical unit number for LSDALTON.INP
    integer,intent(in) :: input
    !> Logical unit number for DALTON.OUT
    integer,intent(in) :: output
    !> Word read from input
    character(len=80),intent(inout) :: word
    character(len=80) :: myword
    !> Is this a full calculation (fullcalc=true, input **CC) 
    !> or a DEC calculation (fullcalc=false, input=**DEC)
    logical,intent(in) :: fullcalc
    !> do we do F12 calc (is a CABS basis required?)
    logical,intent(inout) :: doF12
    logical,save :: already_called = .false.
    integer :: nworkers

    ! Sanity check that this routine is only called once for either **DEC OR **CC
    if(already_called) then
       call lsquit('Error: LSDALTON.INP must contain EITHER **DEC for DEC calculation OR &
            & **CC for full molecular calculation!',-1)
    else
       ! First call to this routine
       already_called=.true.
    end if

#ifdef VAR_MPI
    ! Number of workers = Number of nodes minus master itself
    nworkers = infpar%nodtot -1
    if(nworkers<1.and..not.fullcalc) then
       call lsquit('DEC calculations using MPI require at least two MPI processes!',-1)
    end if
#endif

    ! Just to be sure, we set the default values before
    ! applying input values.
    call dec_set_default_config(output)

    ! Is this a DEC or a full calculation?
    DECinfo%full_molecular_cc = fullcalc

    ! Do DEC calculation
    DECinfo%doDEC=.true.
    DECinfo%output =output

    DO

       IF(READWORD) THEN
          READ (input, '(A80)') WORD
          call capitalize_string(word)
          READWORD=.TRUE.
       ENDIF

       IF ((WORD(1:1) .EQ. '!') .OR. (WORD(1:1) .EQ. '#')) CYCLE

       IF(WORD(1:2) .EQ. '**') THEN
          READWORD=.FALSE.
          EXIT
       ENDIF

       IF(WORD(1:13) == '*END OF INPUT') THEN
          EXIT
       END IF



       ! DEC INPUT INFO
       ! **************

       ! See explanation of DEC parameters in type DEC_settings
       DEC_INPUT_INFO: SELECT CASE(WORD)


          ! ****************************************************************************
          ! *               Keywords available for the general user                    *
          ! ****************************************************************************
          ! These keywords should be properly documented for the release.

          ! GENERAL INFO
          ! ============


          ! DEC orbital-based
       case('.DECCO')
          DECinfo%DECCO=.true.

          ! CC model
       case('.MP2') 
          call find_model_number_from_input(word, DECinfo%ccModel)
          DECinfo%use_singles = .false.  
          DECinfo%NO_MO_CCSD  = .true.
       case('.RIMP2') 
          call find_model_number_from_input(word, DECinfo%ccModel)
          DECinfo%use_singles = .false.  
          DECinfo%NO_MO_CCSD  = .true.
       case('.CC2')
          call find_model_number_from_input(word, DECinfo%ccModel)
          DECinfo%use_singles=.true. 
       case('.CCD') 
          call find_model_number_from_input(word, DECinfo%ccModel)
          DECinfo%CCDhack=.true.
          DECinfo%use_singles=.true. 
          DECinfo%solver_par=.true.
          DECinfo%NO_MO_CCSD  = .true.
       case('.CCSD')
          call find_model_number_from_input(word, DECinfo%ccModel)
          DECinfo%use_singles=.true.; DECinfo%solver_par=.true.
       case('.CCSD(T)') 
          call find_model_number_from_input(word, DECinfo%ccModel)
          DECinfo%use_singles=.true.; DECinfo%solver_par=.true.
       case('.RPA')
          call find_model_number_from_input(word, DECinfo%ccModel)
!#ifdef VAR_MPI
          DECinfo%use_singles=.false.
!#else
!          DECinfo%use_singles=.false.; DECinfo%CCDEBUG=.true.
!#endif


          ! CC SOLVER INFO
          ! ==============

          ! Save CCSD amplitudes to be able to restart full CCSD calculation
       case('.CCSDNOSAFE'); DECinfo%CCSDnosaferun=.true.

          ! Maximum number of CC iterations
       case('.CCMAXITER'); read(input,*) DECinfo%ccMaxIter 

          ! Residual norm threshold for CC amplitude equation
       case('.CCTHR') 
          read(input,*) DECinfo%ccConvergenceThreshold
          DECinfo%CCthrSpecified=.true.
          ! Number of residual vectors to save when solving CC amplitude equation
       case('.SUBSIZE'); read(input,*) DECinfo%ccMaxDIIS

          ! CCSD(T) INFO
          ! ==============
       case('.PT_ABC'); DECinfo%abc=.true.
       case('.ABC_TILE'); read(input,*) DECinfo%abc_tile_size

          ! CHOICE OF ORBITALS
          ! ==================
          ! By default orbitals from the lcm_orbitals.u file are used for DEC or full calculation.
          ! Canonical orbitals can be invoked by this keyword
       case('.CANONICAL') 
          DECinfo%use_canonical=.true.


          ! DEC CALCULATION 
          ! ===============

          ! Restart DEC calculation (only for single point calculations, not geometry optimization)
       case('.RESTART') 
          DECinfo%HFrestart=.true.           
          DECinfo%DECrestart=.true.           

          ! Use HF info generated from previous calculation but run DEC calculation from scratch
       case('.HFRESTART') 
          DECinfo%HFrestart=.true.           
          DECinfo%DECrestart=.false.           

          ! Do not absorb H atoms when assigning orbitals
       case('.NOTABSORBH'); DECinfo%AbsorbHatoms=.false.

          ! FOT
       case('.FOT')  
          read(input,*) DECinfo%FOT

          !> Correlation density for fragment-adapted orbitals, see DECsettings type definition.
       case('.CORRDENS')  
          read(input,*) DECinfo%CorrDensScheme

          ! Use frozen core approximation
       case('.FROZENCORE') 
          DECinfo%frozencore=.true.

          ! Max memory available on node measured in GB. By default set to 16 GB
       case('.MEMORY') 
          read(input,*) DECinfo%memory           
          DECinfo%memory_defined=.true.

       case('.USE_SYS_MEM_INFO') 
          DECinfo%use_system_memory_info = .true.

       case('.PAIRTHR') 
          ! Threshold in a.u.
          read(input,*) DECinfo%pair_distance_threshold
       case('.PAIRTHRANGSTROM') 
          ! Input in Angstrom
          read(input,*) DECinfo%pair_distance_threshold
          DECinfo%pair_distance_threshold=DECinfo%pair_distance_threshold/bohr_to_angstrom

          ! Calculate MP2 gradient (without necessarily doing geometry optimization)
       case('.GRADIENT') 
          DECinfo%gradient=.true.
          DECinfo%first_order=.true.

          !> Carry out MP2 density calculation (subset of gradient calculation)
       case('.DENSITY') 
          DECinfo%density=.true.
          DECinfo%first_order=.true.

          ! Threshold for residual norm of kappabar multiplier equation in first-order MP2 calculations
       case('.KAPPATHR') 
          read(input,*) DECinfo%kappaTHR



          ! ****************************************************************************
          ! *               Keywords only available for developers                     *
          ! ****************************************************************************

          ! Keywords only used for testing in release branch, not intended to be used by end-users,
          ! so on purpose there is no documentation for those in the LSDALTON manual.

       !general testing
       case('.TESTARRAY'); DECinfo%tensor_test=.true.
       case('.TESTREORDERINGS'); DECinfo%reorder_test=.true.
    
       !CCSD testing
       case('.CCSDFORCE_SCHEME'); DECinfo%force_scheme=.true.
          read(input,*) DECinfo%en_mem
       case('.CCSD_DEBUG_COMMUNICATION'); DECinfo%CCSD_NO_DEBUG_COMM   = .false.

       case('.MANUAL_BATCHSIZES') 
          DECinfo%manual_batchsizes=.true.
          read(input,*) DECinfo%ccsdAbatch, DECinfo%ccsdGbatch
       case('.MPISPLIT'); read(input,*) DECinfo%MPIsplit
       case('.INCLUDEFULLMOLECULE');DECinfo%InclFullMolecule=.true.
          ! Size of local groups in MPI scheme
       case('.MPIGROUPSIZE') 
          read(input,*) DECinfo%MPIgroupsize
#ifndef VAR_MPI
          print *, 'WARNING: You have specified MPI groupsize - but this is a serial run!'
          print *, '--> Hence, this keyword has no effect.'
          print *
#endif
       ! Those keyword should probably not be released but are necesary for some essantial test cases:
       case('.CRASHCALC') 
          DECinfo%CRASHCALC= .true.
       case('.STRESSTEST')     
          !Calculate biggest 2 atomic fragments and the biggest pair fragment
          DECinfo%StressTest  = .true.
       case('.FRAG_INIT_SIZE'); read(input,*) DECinfo%Frag_Init_Size
       case('.FRAG_EXP_SIZE'); read(input,*) DECinfo%Frag_Exp_Size
       case('.FRAG_RED_OCC_THR'); read(input,*) DECinfo%frag_red_occ_thr
       case('.FRAG_RED_VIRT_THR'); read(input,*) DECinfo%frag_red_virt_thr
       case('.PRINTFRAGS')
          ! Print fragment energies for full molecular cc calculation
          DECinfo%print_frags = .true.


#ifdef MOD_UNRELEASED

       !CCSD SPECIFIC KEYWORDS
       !**********************
       case('.CCDEBUG');                  DECinfo%CCDEBUG              = .true.
       case('.CCSOLVER_LOCAL');           DECinfo%solver_par           = .false.
       case('.CCSDDYNAMIC_LOAD');         DECinfo%dyn_load             = .true.
       case('.CCSDNODYNAMIC_LOAD');       DECinfo%dyn_load             = .false.
       case('.CCSDNO_RESTART');           DECinfo%CCSDno_restart       = .true.
       case('.SPAWN_COMM_PROC');          DECinfo%spawn_comm_proc      = .true.
       case('.CCSDMULTIPLIERS');          DECinfo%CCSDmultipliers      = .true.
       case('.PNO_S_ON_THE_FLY');         DECinfo%pno_S_on_the_fly     = .true.
       case('.USE_PNOS');                 DECinfo%use_pnos             = .true.
       case('.NOPNOTRAFO');               DECinfo%noPNOtrafo           = .true.; DECinfo%noPNOtrunc=.true.
       case('.NOPNOTRUNCATION');          DECinfo%noPNOtrunc           = .true.
       case('.NOFATRAFO');                DECinfo%noFAtrafo            = .true.; DECinfo%noFAtrunc=.true.
       case('.NOFATRUNCATION');           DECinfo%noFAtrunc            = .true.
       case('.NOPNOOVERLAPTRUNCATION');   DECinfo%noPNOoverlaptrunc    = .true.
       case('.NO_MO_CCSD');               DECinfo%NO_MO_CCSD           = .true.
       case('.PNO_DEBUG');                DECinfo%PNOtriangular        = .false.
       case('.CCSDPREVENTCANONICAL');     DECinfo%CCSDpreventcanonical = .true.
       case('.CCSDEXPL');                 DECinfo%ccsd_expl            = .true.

       case('.PNOTHR');        read(input,*) DECinfo%simplePNOthr
       case('.EOSPNOTHR');     read(input,*) DECinfo%EOSPNOthr
       case('.PNOOVERLAPTHR'); read(input,*) DECinfo%PNOoverlapthr



       !OTHER STUFF FIXME: SORT IT INTO BLOCKS
       !***********

       case('.HACK'); DECinfo%hack=.true.
       case('.HACK2'); DECinfo%hack2=.true.
       case('.V2O2_FREE_SOLVER'); DECinfo%v2o2_free_solver= .true.
       case('.READDECORBITALS'); DECinfo%read_dec_orbitals=.true.
       case('.FRAGEXPMODEL') 
          read(input,*) myword
          call find_model_number_from_input(myword,DECinfo%fragopt_exp_model)
       case('.FRAGREDMODEL') 
          read(input,*) myword
          call find_model_number_from_input(myword,DECinfo%fragopt_red_model)
       case('.TIMEBACKUP'); read(input,*) DECinfo%TimeBackup
       case('.ONLYOCCPART'); DECinfo%OnlyOccPart=.true.
       case('.ONLYVIRTPART'); DECinfo%OnlyVirtPart=.true.

       case('.F12')
          DECinfo%F12=.true.; doF12 = .TRUE.
       case('.F12DEBUG')     
          DECinfo%F12=.true.
          DECinfo%F12DEBUG=.true.
          doF12 = .TRUE.
       case('.PUREHYDROGENDEBUG')     
          DECinfo%PureHydrogenDebug       = .true.
       case('.SOSEX')
         DECinfo%SOS = .true.
       case('.NOTPREC'); DECinfo%use_preconditioner=.false.
       case('.NOTBPREC'); DECinfo%use_preconditioner_in_b=.false.
       case('.MULLIKEN'); DECinfo%mulliken=.true.
       case('.DISTANCE'); DECinfo%distance=.true.
       case('.BOUGHTONPULAY'); DECinfo%BoughtonPulay=.true.
       case('.NOTFITORBITALS'); DECinfo%FitOrbitals=.false.
       case('.SIMPLEORBITALTHRESH')
          read(input,*) DECinfo%simple_orbital_threshold
       case('.DIIS'); DECinfo%use_crop=.false.  ! use DIIS instead of CROP
       case('.MAXITER'); read(input,*) DECinfo%MaxIter
       case('.DECPRINT'); read(input,*) DECinfo%PL
       case('.MULLIKENTHR'); read(input,*) DECinfo%mulliken_threshold
       case('.CHECKPAIRS') 
          DECinfo%checkpairs=.true.
       case('.PAIRMINDIST'); read(input,*) DECinfo%PairMinDist
       case('.PAIRFOTHR'); read(input,*) DECinfo%pairFOthr
       case('.PAIRMP2'); DECinfo%PairMP2=.true.
       case('.NOTPAIRESTIMATE'); DECinfo%PairEstimate=.false.
       case('.IGNOREPAIRESTIMATE'); DECinfo%PairEstimateIgnore=.true.
       case('.ESTIMATEINITRADIUS')
          read(input,*) DECinfo%EstimateINITradius
          DECinfo%EstimateINITradius = DECinfo%EstimateINITradius/bohr_to_angstrom
       case('.ESTIMATEINITATOM')
          read(input,*) DECinfo%EstimateInitAtom
       case('.PAIRMINDISTANGSTROM')
          read(input,*) DECinfo%PairMinDist
          DECinfo%PairMinDist = DECinfo%PairMinDist/bohr_to_angstrom
       case('.PURIFICATION'); DECinfo%PurifyMOs=.true.
       case('.PRECWITHFULL'); DECinfo%precondition_with_full=.true.
       case('.SIMPLEMULLIKENTHRESH'); DECinfo%simple_mulliken_threshold=.true.
       case('.NORMTHRESH'); read(input,*) DECinfo%approximated_norm_threshold
       case('.SIMULATEFULL'); DECinfo%simulate_full=.true.
       case('.SIMULATE_NATOMS'); read(input,*) DECinfo%simulate_natoms
       case('.SKIPREADIN'); DECinfo%SkipReadIn=.true.
       case('.SINGLESPOLARI'); DECinfo%SinglesPolari=.true.
       case('.SINGLESTHR'); read(input,*) DECinfo%SinglesThr
       case('.CONVERT64TO32')
          DECinfo%convert64to32=.true.
       case('.CONVERT32TO64')
          DECinfo%convert32to64=.true.
       case('.ARRAY4ONFILE') 
          DECinfo%array4OnFile=.true.
          DECinfo%array4OnFile_specified=.true.
       case('.FRAG_EXP_SCHEME'); read(input,*) DECinfo%Frag_Exp_Scheme
       case('.FRAG_REDOCC_SCHEME'); read(input,*) DECinfo%Frag_RedOcc_Scheme
       case('.FRAG_REDVIR_SCHEME'); read(input,*) DECinfo%Frag_RedVir_Scheme
       case('.FRAGMENTADAPTED'); DECinfo%fragadapt = .true.
       case('.NO_ORB_BASED_FRAGOPT'); DECinfo%no_orb_based_fragopt = .true.
       case('.ONLY_N_JOBS')
          read(input,*)DECinfo%only_n_frag_jobs
          call mem_alloc(DECinfo%frag_job_nr,DECinfo%only_n_frag_jobs)
          read(input,*)DECinfo%frag_job_nr(1:DECinfo%only_n_frag_jobs)
       case('.ONLY_PAIR_FRAG_JOBS'); DECinfo%only_pair_frag_jobs = .true.

          ! kappabar multiplier equation
       case('.KAPPAMAXITER'); read(input,*) DECinfo%kappaMaxIter 
       case('.KAPPAMAXDIIS'); read(input,*) DECinfo%kappaMaxDIIS
       case('.KAPPA_DEBUG'); DECinfo%kappa_driver_debug=.true.
       case('.NOTKAPPAPREC'); DECinfo%kappa_use_preconditioner=.false.
       case('.NOTKAPPABPREC'); DECinfo%kappa_use_preconditioner_in_b=.false.

          ! Check that input orbitals are orthogonal (debug)
       case('.CHECKLCM'); DECinfo%check_lcm_orbitals=.true.
       case('.CHECKSUBSYSTEMLOC'); DECinfo%check_Occ_SubSystemLocality=.true.
       case('.FORCESUBSYSTEMLOC'); DECinfo%force_Occ_SubSystemLocality=.true.

          !> Collect fragment contributions to calculate full molecular MP2 density
       case('.SKIPFULL') 
          DECinfo%SkipFull=.true.
       case('.ERRORFACTOR') 
          read(input,*) DECinfo%EerrFactor
       case('.CCDRIVERDEBUG')
          DECinfo%cc_driver_debug=.true.
#endif

       CASE DEFAULT
          WRITE (output,'(/,3A,/)') ' Keyword "',WORD,&
               & '" not recognized in config_dec_input'
          CALL lsQUIT('Illegal keyword in config_dec_input',output)

       END SELECT DEC_INPUT_INFO

    ENDDO

    ! Check that input is consistent to avoid weird segmentation fault etc...
    call check_dec_input()


  END SUBROUTINE config_dec_input


  !> \brief Check that DEC input is consistent with what is currently implemented
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine check_dec_input()
    implicit none
    integer :: nodtot
    nodtot = 1
#ifdef VAR_MPI
    nodtot = infpar%nodtot
#endif


    ! DEC orbital-based - currently limited to occupied partitioning scheme
    ! and several options are not possible
    DoDECCO: if(DECinfo%DECCO) then

       ! Occupied partitioning scheme
!!$       if(.not. DECinfo%OnlyOccPart) then
!!$          print *, 'WARNING: DECCO is implemented only for occ partitioning scheme!'
!!$          print *, '--> I will only use occupied partitioning scheme.'
!!$          write(DECinfo%output,*) 'WARNING: DECCO only for occ partitioning scheme!'
!!$          write(DECinfo%output,*) '--> I will only use occupied partitioning scheme.'
!!$          DECinfo%OnlyOccPart=.true.
!!$          DECinfo%OnlyVirtPart=.false.
!!$       end if

       ! Not simulate full
       if(DECinfo%simulate_full) then
          call lsquit('DECCO not implemented for SIMULATEFULL',-1)
       end if

!!$       ! Not working for first-order properties
!!$       if(DECinfo%first_order) then
!!$          call lsquit('DECCO is not implemented for first-order properties!',DECinfo%output)
!!$       end if


       ! No stress test implemented
       if( DECinfo%StressTest ) then
          call lsquit('DECCO is not implemented for stress test!',DECinfo%output)
       end if

       if(DECinfo%SinglesPolari) then
          call lsquit('DECCO is not implemented for singles polarization effects!',DECinfo%output)
       end if

    end if DoDECCO


    ! Check that array4OnFile is only called for the cases where it is implemented
    ArraysOnFile: if(DECinfo%array4OnFile) then

       ! Only for MP2 so far
       if(DECinfo%ccModel /= MODEL_MP2) then
          call lsquit('Storing arrays on file only implemented for MP2. &
               & Suggestion: Remove .array4OnFile keyword!', DECinfo%output)
       end if

       if(DECinfo%use_singles) then
          call lsquit('Storing arrays on file not implemented for singles. &
               & Suggestion: Remove .array4OnFile keyword!', DECinfo%output)
       end if


    end if ArraysOnFile


    FirstOrderModel: if(DECinfo%ccModel /= MODEL_MP2.and.DECinfo%ccModel /= MODEL_CCSD) then


       if(DECinfo%density) then
          call lsquit('Calculation of density matrix is only implemented for MP2/CCSD!', DECinfo%output)
       end if

       if(DECinfo%gradient) then
          call lsquit('Calculation of molecular gradient is only implemented for MP2/CCSD!', DECinfo%output)
       end if

    end if FirstOrderModel



    MP2gradientCalculation: if(DECinfo%first_order) then

       if(DECinfo%full_molecular_cc) then
          call lsquit('Full calculation for MP2 gradient is implemented via the &
               & .SimulateFull keyword', DECinfo%output)
       end if

       if(DECinfo%onlyoccpart) then
          call lsquit('DEC gradient cannot be evaluated when only occupied &
               & partitioning scheme is used!',DECinfo%output)
       end if
       if(DECinfo%onlyvirtpart) then
          call lsquit('DEC gradient cannot be evaluated when only virtual &
               & partitioning scheme is used!',DECinfo%output)
       end if

       !Make sure, that if first order is specified, we calculate the CCSD
       !multipliers
       if( DECinfo%ccmodel == MODEL_CCSD)then
          DECinfo%CCSDmultipliers = .true.
       endif

    end if MP2gradientCalculation


    ! If simulate full calculation, the full molecule must be included in "fragments"
    SimulateFullCalc: if(DECinfo%simulate_full) then
       DECinfo%InclFullMolecule = .true.
    end if SimulateFullCalc

    ! Never ignores pairs when full molecule is included in fragments
    if(DECinfo%InclFullMolecule) then
       DECinfo%PairEstimateIgnore=.true.
    end if


    ! Set CC residual threshold to be 0.01*FOT
    ! - unless it was specified explicitly in the input.
    if(.not. DECinfo%CCthrSpecified) then
       DECinfo%ccConvergenceThreshold=0.01E0_realk*DECinfo%FOT
    end if

    ! Never use gradient and density at the same time (density is a subset of gradient)
    if(DECinfo%density .and. DECinfo%gradient) then
       call lsquit('Density and gradient cannot both be turned on at the same time! &
            & Note that density is a subset of a gradient calculation',DECinfo%output)
    end if

    if(DECinfo%SinglesPolari) then
       call lsquit('Full singles polarization has been temporarily disabled!',-1)
    end if

    if((.not. (DECinfo%memory_defined .or.  DECinfo%use_system_memory_info ) )&
         & .or. (DECinfo%memory_defined .and. DECinfo%use_system_memory_info ) ) then

       write(DECinfo%output,*) 'Memory not or multiply defined for **DEC or **CC calculation!'
       write(DECinfo%output,*) 'Please specify using EITHER .MEMORY keyword (in gigabytes) OR .USE_SYS_MEM_INFO'
       write(DECinfo%output,*) 'The recommended way is using .MEMORY and specifying the memory in GB'
#ifdef VAR_MPI
       write(DECinfo%output,*) 'E.g. if each MPI process has 16 GB of memory available, then use'
#else
       write(DECinfo%output,*) 'E.g. if there are 16 GB of memory available, then use'
#endif
       write(DECinfo%output,*) '.MEMORY'
       write(DECinfo%output,*) '16.0'
       write(DECinfo%output,*) ''
       call lsquit('**DEC or **CC calculation requires specification of available memory using &
            & EITHER .MEMORY OR .USE_SYS_MEM_INFO  keyword!',-1)
    end if

    ! Use purification of FOs when using fragment-adapted orbitals.
    if(DECinfo%fragadapt) then
       DECinfo%purifyMOs=.true.
    end if

    if(DECinfo%use_system_memory_info) call get_currently_available_memory(DECinfo%memory)

    ! Check in the case of a DEC calculation that the cc-restart-files are not written
    if((.not.DECinfo%full_molecular_cc).and.(.not.DECinfo%CCSDnosaferun))then
       DECinfo%CCSDnosaferun = .true.
    endif

    if( (.not.DECinfo%full_molecular_cc) .and. DECinfo%ccmodel == MODEL_MP2 .and. &
         &(    DECinfo%fragopt_exp_model == MODEL_CC2 &
         &.or. DECinfo%fragopt_red_model == MODEL_CC2 &
         &.or. DECinfo%fragopt_exp_model == MODEL_CCSD &
         &.or. DECinfo%fragopt_red_model == MODEL_CCSD &
         &.or. DECinfo%fragopt_exp_model == MODEL_CCSDpT &
         &.or. DECinfo%fragopt_red_model == MODEL_CCSDpT )                          ) then
       call lsquit('The specification of .MP2 and .FRAGEXPMODEL > .MP2 or .FRAGREDMODEL > .MP2&
            & does not make sense, please change input!',-1)

    endif

    if((.not.DECinfo%full_molecular_cc).and.DECinfo%force_scheme)then
       call lsquit("ERROR(check_dec_input):Do not use &
            &.CCSDforce_scheme in a DEC calculation",-1)
    endif

    if((DECinfo%full_molecular_cc).and.DECinfo%force_scheme.and.(DECinfo%en_mem==1.or.&
         &DECinfo%en_mem==2.or.DECinfo%en_mem==3).and.nodtot==1)then
       call lsquit("ERROR(check_dec_input):You forced a scheme in &
            &the CCSD part which is dependent on running at least 2 &
            &MPI processes with only one process",-1)
    endif


    ! Set FOTs for geometry opt.
    call set_geoopt_FOTs(DECinfo%FOT)

  end subroutine check_dec_input

  !> \brief Check that CC input is consistent with calc requirements
  !> \author Thomas Kjaergaard
  !> \date October 2014
  subroutine check_cc_input(mylsitem,nocc,nvirt,nbasis)
    implicit none
    type(lsitem),intent(inout) :: mylsitem
    integer,intent(in) :: nocc,nvirt,nbasis
    !
    real(realk) :: OO,VV,BB,AA,intMEM, solMEM,mem_required,GB
    integer     :: intstep,nthreads
    ! Number of OMP threads
#ifdef VAR_OMP
    integer, external :: OMP_GET_MAX_THREADS
    nthreads=OMP_GET_MAX_THREADS()
#else
    ! No OMP, set number of threads to one
    nthreads=1
#endif

    GB = 1.0E+9_realk !1GB
    SELECT CASE(DECinfo%ccModel)
    CASE(MODEL_MP2)
!CODE OBSOLETE DUE TO PATRICK NEW FANCY MP2 CODE
!       OO=nocc      ! Number of occupied orbitals (as real)
!       VV=nvirt     ! Number of virtual orbitals (as real)
!       ! Maximum batch dimension (as real)
!       BB=max_batch_dimension(mylsitem,nbasis)
!       AA=nbasis    ! Number of atomic orbitals (as real)       
!       call estimate_memory_for_mp2_energy(nthreads,OO,VV,AA,BB,intMEM,intStep,solMEM)
!       mem_required = max(intMEM,solMEM)
!       mem_required = mem_required + DECinfo%fullmolecule_memory
!       mem_required = nocc*nvirt*nocc*nvirt*8.0E0_realk/GB
!       IF(mem_required.GT.DECinfo%memory)THEN
!          CALL FullMemoryError(mem_required)
!          call lsquit('Memory specification too small',DECinfo%output)
!       ENDIF
!    CASE(MODEL_CC2)
!    CASE(MODEL_CCSD)
!    CASE(MODEL_CCSDpT)
    case default
    end SELECT
  end subroutine check_cc_input

  subroutine FullMemoryError(nsize)
    implicit none
    real(realk) :: nsize
    WRITE(DECinfo%output,'(A)')'Error in Memory specification. '
    WRITE(DECinfo%output,'(A)')'The memory specified using the .MEMORY keyword'
    WRITE(DECinfo%output,'(A)')'is not big enough for the calculations requirements '
    WRITE(DECinfo%output,'(A,F10.2,A)')'Requirements    :',nsize,' Gb'
    WRITE(DECinfo%output,'(A,I12,A)')  'Memory specified:',DECinfo%memory,' Gb'
  end subroutine FullMemoryError

  subroutine DEC_settings_print(DECitem,lupri)
    type(DECsettings) :: DECitem
    integer,intent(in) :: lupri

    WRITE(lupri,*) ' The DEC settings structure '
    write(lupri,*) '****************************'

    write(lupri,*) 

    write(lupri,*) 'doDEC ', DECitem%doDEC
    write(lupri,*) 'frozencore ', DECitem%frozencore
    write(lupri,*) 'full_molecular_cc ', DECitem%full_molecular_cc
    write(lupri,*) 'use_canonical ', DECitem%use_canonical
    write(lupri,*) 'simulate_full ', DECitem%simulate_full
    write(lupri,*) 'simulate_natoms ', DECitem%simulate_natoms
    write(lupri,*) 'InclFullMolecule ', DECitem%InclFullMolecule
    write(lupri,*) 'cc_models ', DECitem%cc_models
    write(lupri,*) 'ccModel ', DECitem%ccModel
    write(lupri,*) 'use_singles ', DECitem%use_singles
    write(lupri,*) 'HFrestart ', DECitem%HFrestart
    write(lupri,*) 'DECrestart ', DECitem%HFrestart
    write(lupri,*) 'TimeBackup ', DECitem%TimeBackup
    write(lupri,*) 'read_dec_orbitals ', DECitem%read_dec_orbitals
    write(lupri,*) 'memory ', DECitem%memory
    write(lupri,*) 'memory_defined ', DECitem%memory_defined
    write(lupri,*) 'fullmolecule_memory ', DECitem%fullmolecule_memory
    write(lupri,*) 'array4OnFile ', DECitem%array4OnFile
    write(lupri,*) 'array4OnFile_specified ', DECitem%array4OnFile_specified
    write(lupri,*) 'SinglesPolari ', DECitem%SinglesPolari
    write(lupri,*) 'singlesthr ', DECitem%singlesthr
    write(lupri,*) 'convert64to32 ', DECitem%convert64to32
    write(lupri,*) 'convert32to64 ', DECitem%convert32to64
    write(lupri,*) 'CCSDnosaferun ', DECitem%CCSDnosaferun
    write(lupri,*) 'solver_par ', DECitem%solver_par
    write(lupri,*) 'force_scheme ', DECitem%force_scheme
    write(lupri,*) 'dyn_load ', DECitem%dyn_load
    write(lupri,*) 'CCDEBUG ', DECitem%CCDEBUG
    write(lupri,*) 'CCSDno_restart ', DECitem%CCSDno_restart
    write(lupri,*) 'CCSDpreventcanonical ', DECitem%CCSDpreventcanonical
    write(lupri,*) 'CRASHCALC            ', DECitem%CRASHCALC
    write(lupri,*) 'cc_driver_debug ', DECitem%cc_driver_debug
    write(lupri,*) 'en_mem ', DECitem%en_mem
    write(lupri,*) 'precondition_with_full ', DECitem%precondition_with_full
    write(lupri,*) 'ccsd_expl ', DECitem%ccsd_expl
    write(lupri,*) 'ccMaxIter ', DECitem%ccMaxIter
    write(lupri,*) 'ccMaxDIIS ', DECitem%ccMaxDIIS
    write(lupri,*) 'ccConvergenceThreshold ', DECitem%ccConvergenceThreshold
    write(lupri,*) 'CCthrSpecified ', DECitem%CCthrSpecified
    write(lupri,*) 'use_preconditioner ', DECitem%use_preconditioner
    write(lupri,*) 'use_preconditioner_in_b ', DECitem%use_preconditioner_in_b
    write(lupri,*) 'use_crop ', DECitem%use_crop
#ifdef MOD_UNRELEASED    
    write(lupri,*) 'F12 ', DECitem%F12
    write(lupri,*) 'F12DEBUG ', DECitem%F12DEBUG
#endif
    write(lupri,*) 'mpisplit ', DECitem%mpisplit
    write(lupri,*) 'MPIgroupsize ', DECitem%MPIgroupsize
    write(lupri,*) 'manual_batchsizes ', DECitem%manual_batchsizes
    write(lupri,*) 'ccsdAbatch,ccsdGbatch ', DECitem%ccsdAbatch,DECitem%ccsdGbatch
    write(lupri,*) 'hack ', DECitem%hack
    write(lupri,*) 'hack2 ', DECitem%hack2
    write(lupri,*) 'SkipReadIn ', DECitem%SkipReadIn
    write(lupri,*) 'tensor_test ', DECitem%tensor_test
    write(lupri,*) 'reorder_test ', DECitem%reorder_test
    write(lupri,*) 'check_lcm_orbitals ', DECitem%check_lcm_orbitals
    write(lupri,*) 'check_Occ_SubSystemLocality ', DECitem%check_Occ_SubSystemLocality
    write(lupri,*) 'force_Occ_SubSystemLocality ', DECitem%force_Occ_SubSystemLocality
    write(lupri,*) 'PL ', DECitem%PL
    write(lupri,*) 'SkipFull ', DECitem%SkipFull
    write(lupri,*) 'output ', DECitem%output
    write(lupri,*) 'AbsorbHatoms ', DECitem%AbsorbHatoms
    write(lupri,*) 'FitOrbitals ', DECitem%FitOrbitals
    write(lupri,*) 'simple_orbital_threshold ', DECitem%simple_orbital_threshold
    write(lupri,*) 'PurifyMOs ', DECitem%PurifyMOs
    write(lupri,*) 'FragAdapt ', DECitem%FragAdapt
    write(lupri,*) 'BoughtonPulay ', DECitem%BoughtonPulay
    write(lupri,*) 'mulliken_threshold ', DECitem%mulliken_threshold
    write(lupri,*) 'simple_mulliken_threshold ', DECitem%simple_mulliken_threshold
    write(lupri,*) 'approximated_norm_threshold ', DECitem%approximated_norm_threshold
    write(lupri,*) 'mulliken ', DECitem%mulliken
    write(lupri,*) 'FOT ', DECitem%FOT
    write(lupri,*) 'MaxIter ', DECitem%MaxIter
    write(lupri,*) 'FOTlevel ', DECitem%FOTlevel
    write(lupri,*) 'Frag_Exp_Scheme ', DECitem%Frag_Exp_Scheme
    write(lupri,*) 'Frag_RedOcc_Scheme ', DECitem%Frag_RedOcc_Scheme
    write(lupri,*) 'Frag_RedVir_Scheme ', DECitem%Frag_RedVir_Scheme
    write(lupri,*) 'Frag_Init_Size ', DECitem%Frag_Init_Size
    write(lupri,*) 'Frag_Exp_Size ', DECitem%Frag_Exp_Size
    write(lupri,*) 'Frag_Red_occ_thr ', DECinfo%frag_red_occ_thr
    write(lupri,*) 'Frag_Red_virt_thr ', DECinfo%frag_red_virt_thr
    write(lupri,*) 'fragopt_exp_model ', DECitem%fragopt_exp_model
    write(lupri,*) 'fragopt_red_model ', DECitem%fragopt_red_model
    write(lupri,*) 'No_Orb_Based_FragOpt ', DECitem%no_orb_based_fragopt
    write(lupri,*) 'pair_distance_threshold ', DECitem%pair_distance_threshold
    write(lupri,*) 'PairMinDist ', DECitem%PairMinDist
    write(lupri,*) 'CheckPairs ', DECitem%CheckPairs
    write(lupri,*) 'pairFOthr ', DECitem%pairFOthr
    write(lupri,*) 'PairMP2 ', DECitem%PairMP2
    write(lupri,*) 'PairEstimate ', DECitem%PairEstimate
    write(lupri,*) 'first_order ', DECitem%first_order
    write(lupri,*) 'density ', DECitem%density
    write(lupri,*) 'gradient ', DECitem%gradient
    write(lupri,*) 'kappa_use_preconditioner ', DECitem%kappa_use_preconditioner
    write(lupri,*) 'kappa_use_preconditioner_in_b ', DECitem%kappa_use_preconditioner_in_b
    write(lupri,*) 'kappaMaxDIIS ', DECitem%kappaMaxDIIS
    write(lupri,*) 'kappaMaxIter ', DECitem%kappaMaxIter
    write(lupri,*) 'kappa_driver_debug ', DECitem%kappa_driver_debug
    write(lupri,*) 'kappaTHR ', DECitem%kappaTHR
    write(lupri,*) 'ncalc ', DECitem%ncalc
    write(lupri,*) 'EerrFactor ', DECitem%EerrFactor
    write(lupri,*) 'EerrOLD ', DECitem%EerrOLD

  end subroutine DEC_settings_print



  !> \brief Set FOTs to (possibly) be used for geometry optimization,
  !> --> sets DECinfo%GeoFOTs such that
  !> DECinfo%GeoFOTs(1) = initial FOT
  !> DECinfo%GeoFOTs(2) = (initial FOT)/10
  !> DECinfo%GeoFOTs(3) = (initial FOT)/100
  !> etc.
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine set_geoopt_FOTs(FOT)
    implicit none
    !> Initial FOT
    real(realk),intent(in) :: FOT
    integer :: i

    ! Initial FOT level is 1
    DECinfo%FOTlevel = 1

    DECinfo%GeoFOTs(1) = FOT
    do i=2,nFOTs
       DECinfo%GeoFOTs(i) = DECinfo%GeoFOTs(i-1)*0.1_realk
    end do

  end subroutine set_geoopt_FOTs

  subroutine free_decinfo()
     implicit none
     if(associated(DECinfo%frag_job_nr))call mem_dealloc(DECinfo%frag_job_nr)
  end subroutine

  !> MODIFY FOR NEW MODEL
  !> \brief For a given model input (e.g. .MP2 or .CCSD) find model number associated with input.
  !> \author Kasper Kristensen
  !> \date November 2013
  subroutine find_model_number_from_input(myword,modelnumber)
    implicit none
    !> Word read from input
    character(len=80),intent(in) :: myword
    !> Model number corresponding to input (see MODEL_* in dec_typedef.F90)
    integer,intent(inout) :: modelnumber

    SELECT CASE(MYWORD)

    case('.MP2');     modelnumber = MODEL_MP2
    case('.CC2');     modelnumber = MODEL_CC2
    case('.CCSD');    modelnumber = MODEL_CCSD
    case('.CCD');     modelnumber = MODEL_CCSD  ! effectively use CCSD where singles amplitudes are zeroed
    case('.CCSD(T)'); modelnumber = MODEL_CCSDpT
    case('.RPA');     modelnumber = MODEL_RPA
    case('.RIMP2');   modelnumber = MODEL_RIMP2
    case default
       print *, 'Model not found: ', myword
       write(DECinfo%output,*)'Model not found: ', myword
       write(DECinfo%output,*)'Models supported are:'
       write(DECinfo%output,*)'.MP2'
       write(DECinfo%output,*)'.CC2'
       write(DECinfo%output,*)'.CCSD'
       write(DECinfo%output,*)'.CCD'
       write(DECinfo%output,*)'.CCSD(T)'
       write(DECinfo%output,*)'.RPA'
       write(DECinfo%output,*)'.RIMP2'
       call lsquit('Requested model not found!',-1)
    end SELECT

  end subroutine find_model_number_from_input


end MODULE DEC_settings_mod
