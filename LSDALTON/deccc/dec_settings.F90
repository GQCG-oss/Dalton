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
  use matrix_module
  use matrix_operations
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type, only: LSMPIASYNCP
#endif

contains

  !> \brief Read the **DEC or **CC input section in LSDALTON.INP and set 
  !> configuration structure accordingly.
  !> \author Kasper Kristensen
  !> \date September 2010
  SUBROUTINE config_dec_input(input,output,readword,word,fullcalc,doF12,doRIMP2)
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
    !> do we do RIMP2 calc (is a AUX basis required?)
    logical,intent(inout) :: doRIMP2
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


       ! SNOOP
       ! =====
       case('.SNOOP') 
          ! Perform SNOOP calculation rather than DEC (will be merged at some point)
          DECinfo%SNOOP=.true.

       case('.SNOOPJUSTHF') 
          ! Just do HF calculation in SNOOP and skip correlated CC calculation?
          DECinfo%SNOOPjustHF=.true.

       case('.SNOOPTHR') 
          ! Threshold for residual norm in SNOOP HF calculations
          read(input,*) DECinfo%SNOOPthr

       case('.SNOOPMAXITER')
          ! Maximum number of iterations in SNOOP HF calculations
          read(input,*) DECinfo%SNOOPMaxIter

       case('.SNOOPMAXDIIS')
          ! Maximum number of DIIS vectors to store in SNOOP HF calculations (RH/DIIS scheme)
          read(input,*) DECinfo%SNOOPMaxDIIS
       case('.SNOOP_DEBUG')
          ! Debug prints for SNOOP
          DECinfo%SNOOPdebug=.true.

       case('.SNOOPNOTSAMESPACE')
          !> Do not use full orbital spaces for monomer calculation as defined by natural connection,
          !> rather simply do independent DEC fragment optimization for monomers.
          DECinfo%SNOOPsamespace=.false.

       case('.SNOOPLOCALIZE')
          DECinfo%SNOOPlocalize=.true.

       case('.SNOOPRESTART')
          DECinfo%SNOOPrestart=.true.
          ! Also restart HF for full molecule by default when using SNOOP restart
          DECinfo%HFrestart=.true.

       case('.SNOOPONESUB')
          read(input,*) DECinfo%SNOOPonesub
          if(DECinfo%SNOOPonesub<0) then
             call lsquit('Error in SNOOPONESUB input!',-1)
          end if

          ! CC RESPONSE
          ! ===========
       case('.EXCITATIONENERGIES')
          DECinfo%CCexci = .true.
          read(input,*) DECinfo%JacobianNumEival

       case('.JACOBIANLEFT')
          DECinfo%JacobianLHTR = .true.

       case('.JACOBIANTHR')
          read(input,*) DECinfo%JacobianThr 

       case('.JACOBIANMAXSUBSPACE')
          read(input,*) DECinfo%JacobianMaxSubspace

       case('.JACOBIANINITSUBSPACE')
          read(input,*) DECinfo%JacobianInitialSubspace

       case('.JACOBIANMAXITER')
          read(input,*) DECinfo%JacobianMaxIter

       case('.JACOBIANNOTPRECOND')
          DECinfo%JacobianPrecond = .false.

       case('.SINGLESEW1')
          DECinfo%SinglesEW1 = .true.

       case('.LW1')
          DECinfo%LW1 = .true.

       case('.P_EOM_MBPT2')
          DECinfo%P_EOM_MBPT2 = .true.


       ! GENERAL INFO
       ! ============
       case('.DECCO')
          ! DEC orbital-based
          DECinfo%DECCO=.true.

       case('.DECNP')
          ! Alternaitve DEC energy formulation with no pairs:
          DECinfo%DECNP=.true.

          !> select CC model
       case('.MP2') 
          call find_model_number_from_input(word, DECinfo%ccModel)
          DECinfo%use_singles = .false.  
          DECinfo%NO_MO_CCSD  = .true.
       case('.RIMP2') 
          call find_model_number_from_input(word, DECinfo%ccModel)
          DECinfo%use_singles = .false.  
          DECinfo%NO_MO_CCSD  = .true.
          doRIMP2 = .TRUE.
       case('.LSTHCRIMP2') 
          call find_model_number_from_input(word, DECinfo%ccModel)
          DECinfo%use_singles = .false.  
          DECinfo%NO_MO_CCSD  = .true.
          doRIMP2 = .TRUE. !we need an Aux basis 
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
          DECinfo%use_singles=.false.
          DECinfo%solver_par=.true.
       case('.SOSEX')
          call find_model_number_from_input(word, DECinfo%ccModel)
          DECinfo%use_singles=.false.
          DECinfo%solver_par=.true.
       case('.MP3') 
          call find_model_number_from_input(word, DECinfo%ccModel)
          DECinfo%use_singles = .false.  
          DECinfo%NO_MO_CCSD  = .true.

       ! CC SOLVER INFO
       ! ==============
       case('.CCSDNOSAFE')
          ! Save CCSD amplitudes to be able to restart full CCSD calculation
          DECinfo%CCSDnosaferun=.true.

       case('.CCMAXITER')
          ! Maximum number of CC iterations
          read(input,*) DECinfo%ccMaxIter 

       case('.CCTHR') 
          ! Residual norm threshold for CC amplitude equation
          read(input,*) DECinfo%ccConvergenceThreshold
          DECinfo%CCthrSpecified=.true.

       case('.SUBSIZE')
          read(input,*) DECinfo%ccMaxDIIS
          ! Number of residual vectors to save when solving CC amplitude equation


       ! CCSD(T) INFO
       ! ==============
       case('.PT_ABC'); DECinfo%abc = .true.
       case('.IJK_TILE'); read(input,*) DECinfo%ijk_tile_size
       case('.ABC_TILE'); read(input,*) DECinfo%abc_tile_size
       case('.NBUFFS_IJK'); read(input,*) DECinfo%ijk_nbuffs
       case('.NBUFFS_ABC'); read(input,*) DECinfo%abc_nbuffs
       case('.ACC_SYNC'); DECinfo%acc_sync = .true.
       case('.PT_SINGLE_PREC'); DECinfo%pt_single_prec = .true.
       case('.PT_HACK'); DECinfo%pt_hack = .true.
       case('.PT_HACK2'); DECinfo%pt_hack2 = .true.
       case('.TEST_LEN'); read(input,*) DECinfo%test_len 

       ! DEC CALCULATION 
       ! ===============
       case('.RESTART') 
          !> Restart DEC calculation (only for single point calculations, not geometry optimization)
          DECinfo%HFrestart=.true.           
          DECinfo%DECrestart=.true.           

       case('.HFRESTART') 
          !> Use HF info generated from previous calculation but run DEC calculation from scratch
          DECinfo%HFrestart=.true.           
          DECinfo%DECrestart=.false.           

       case('.ENFORCERESTART') 
          !> Enforce restart
          DECinfo%HFrestart=.true.           
          DECinfo%DECrestart=.true.           
          DECinfo%EnforceRestart=.true.           

       case('.NOTABSORBH')
          !> Do not absorb H atoms when assigning orbitals
          DECinfo%AbsorbHatoms=.false.

       case('.FOT')  
          !> Fragment Optimization Threshold 
          read(input,*) DECinfo%FOT

          ! Use frozen core approximation
       case('.FROZENCORE') 
          DECinfo%frozencore=.true.

          ! Max memory available on node measured in GB. By default set to 16 GB
       case('.MEMORY') 
          read(input,*) DECinfo%memory           
          DECinfo%memory_defined=.true.

       case('.BG_MEMORY') 
          read(input,*) DECinfo%bg_memory           

       case('.USE_SYS_MEM_INFO') 
          DECinfo%use_system_memory_info = .true.

       case('.GRADIENT') 
          ! Calculate MP2 gradient (without necessarily doing geometry optimization)
          DECinfo%gradient=.true.
          DECinfo%first_order=.true.

       case('.DENSITY') 
          !> Carry out MP2 density calculation (subset of gradient calculation)
          DECinfo%density=.true.
          DECinfo%first_order=.true.

       case('.KAPPATHR') 
          ! Threshold for residual norm of kappabar multiplier equation in first-order MP2 calculations
          read(input,*) DECinfo%kappaTHR

       case('.DECPRINT')
          ! DEC print level
          read(input,*) DECinfo%PL

       case('.MEMDEBUGPRINT')
          ! DEC print level
          DECinfo%MemDebugPrint=.true.

       case('.ONLYOCCPART')
          ! Use only occupied partitioning scheme
          DECinfo%OnlyOccPart=.true.

       case('.ONLYVIRTPART')
          ! Use only virtual partitioning scheme
          DECinfo%OnlyVirtPart=.true.



       ! CHOICE OF ORBITALS
       ! ==================
       case('.CANONICAL') 
          ! By default orbitals from the lcm_orbitals.u file are used for DEC or full calculation.
          ! Canonical orbitals can be invoked by this keyword
          DECinfo%use_canonical=.true.


       ! ****************************************************************************
       ! *               Keywords only available for developers                     *
       ! ****************************************************************************

       ! Keywords only used for testing in release branch, not intended to be used by end-users,
       ! so on purpose there is no documentation for those in the LSDALTON manual.


       !KEYWORDS FOR DEC DEBUGGING AND TESTING
       !**************************************
       case('.HACK'); DECinfo%hack=.true.
       case('.HACK2'); DECinfo%hack2=.true.
       case('.TESTARRAY'); DECinfo%tensor_test=.true.
       case('.TESTREORDERINGS'); DECinfo%reorder_test=.true.
       case('.INCLUDEFULLMOLECULE');DECinfo%InclFullMolecule=.true.
       case('.SIMULATEFULL'); DECinfo%simulate_full=.true.
       case('.SIMULATE_NATOMS'); read(input,*) DECinfo%simulate_natoms
       case('.CRASHCALC'); DECinfo%CRASHCALC=.true.
       case('.CRASHESTI'); DECinfo%CRASHESTI=.true.
       case('.PUREHYDROGENDEBUG'); DECinfo%PureHydrogenDebug=.true.
       case('.STRESSTEST')     
          !Calculate biggest 2 atomic fragments and the biggest pair fragment
          DECinfo%StressTest = .true.
       case('.PRINTFRAGS')
          ! Print fragment energies for full molecular cc calculation
          DECinfo%print_frags = .true.
       ! Check that input orbitals are orthogonal (debug)
       case('.CHECKLCM'); DECinfo%check_lcm_orbitals=.true.
       case('.CHECKSUBSYSTEMLOC'); DECinfo%check_Occ_SubSystemLocality=.true.
       case('.FORCESUBSYSTEMLOC'); DECinfo%force_Occ_SubSystemLocality=.true.


    

       !KEYWORDS FOR DEC PARALLELISM
       !****************************
       case('.TEST_FULLY_DISTRIBUTED_INTEGRALS') 
          DECinfo%test_fully_distributed_integrals=.true.
       case('.MANUAL_BATCHSIZES') 
          DECinfo%manual_batchsizes=.true.
          read(input,*) DECinfo%ccsdAbatch, DECinfo%ccsdGbatch
       case('.MANUAL_OCCBATCHSIZES') 
          DECinfo%manual_occbatchsizes=.true.
          read(input,*) DECinfo%batchOccI, DECinfo%batchOccJ
       case('.MPISPLIT'); read(input,*) DECinfo%MPIsplit
       case('.RIMPISPLIT'); read(input,*) DECinfo%RIMPIsplit
       case('.MPIGROUPSIZE') 
          read(input,*) DECinfo%MPIgroupsize
#ifndef VAR_MPI
          print *, 'WARNING: You have specified MPI groupsize - but this is a serial run!'
          print *, '--> Hence, this keyword has no effect.'
          print *
#endif
       case('.DISTRIBUTE_FULLINFO')
          DECinfo%force_distribution      = .true.
          DECinfo%distribute_fullmolecule = .true.
       case('.NOT_DISTRIBUTE_FULLINFO')
          DECinfo%force_distribution      = .true.
          DECinfo%distribute_fullmolecule = .false.


       !KEYWORDS RELATED TO FRAGMENT SPACES
       !***********************************
       case('.FRAGEXPMODEL') 
          ! CC model used in the expansion part of the fragment optimization
          read(input,*) myword
          call find_model_number_from_input(myword,DECinfo%fragopt_exp_model)

       case('.FRAGREDMODEL') 
          ! CC model used in the reduction part of the fragment optimization
          read(input,*) myword
          call find_model_number_from_input(myword,DECinfo%fragopt_red_model)

       case('.FRAG_EXP_SCHEME')
          ! Orbital list used in the expanion: (see define_frag_expansion for details)
          read(input,*) DECinfo%Frag_Exp_Scheme
          DECinfo%use_abs_overlap = (DECinfo%Frag_Exp_Scheme==3)

       case('.FRAG_REDOCC_SCHEME')
          ! Occupied orbital list used in the expanion: (see define_frag_reduction for details)
          read(input,*) DECinfo%Frag_RedOcc_Scheme

       case('.FRAG_REDVIR_SCHEME')
          ! Virtual orbital list used in the expanion: (see define_frag_reduction for details)
          read(input,*) DECinfo%Frag_RedVir_Scheme

       case('.FRAG_INIT_SIZE')
          ! Number of "average" atoms used to initialize a fragment (Excluded EOS)
          ! (one "average" atom corresponds to the average number of orbital per atom)
          read(input,*) DECinfo%Frag_Init_Size

       case('.FRAG_EXP_SIZE')
          ! Number of "average" atoms used to expand a fragment
          ! (one "average" atom corresponds to the average number of orbital per atom)
          read(input,*) DECinfo%Frag_Exp_Size

       case('.FRAG_RED_OCC')
          ! Start reducting occupied space first
          DECinfo%frag_red_occ  = .true.
       case('.FRAG_RED_VIRT')
          ! Start reducting virtual space first
          DECinfo%frag_red_virt = .true.
       case('.FRAG_RED1_THR')
          ! Threshold for convergence of first reduced space (occ or virt) 
          ! should be a number between 0 and 1 which is then multiplied to the FOT
          read(input,*) DECinfo%frag_red1_thr
       case('.FRAG_RED2_THR')
          ! Threshold for convergence of second reduced space (occ or virt)
          ! should be a number between 0 and 1 which is then multiplied to the FOT
          read(input,*) DECinfo%frag_red2_thr

       case('.FRAGMENTADAPTED')
          ! Fragment adapted orbital instead of reduction (??)
          DECinfo%fragadapt = .true.

       case('.CORRDENS')  
          !> Correlation density for fragment-adapted orbitals, see DECsettings type definition.
          read(input,*) DECinfo%CorrDensScheme

       case('.FRACOFORBSPACE_RED')
          ! set the fraction of the fully extended orbital space that is used as 
          ! tolerance in an incomplete binary search
          read(input,*) DECinfo%FracOfOrbSpace_red

       case('.FRAG_INIT_RADIUS_NO_OPT_ALL')
          ! include all orbitals for a fragment within a given radius and calculate 
          ! the fragment energies in Angstrom
          read(input,*) DECinfo%all_init_radius
          DECinfo%all_init_radius = DECinfo%all_init_radius/bohr_to_angstrom
          DECinfo%occ_init_radius = DECinfo%all_init_radius
          DECinfo%vir_init_radius = DECinfo%all_init_radius
       case('.FRAG_INIT_RADIUS_NO_OPT_OCC')
          read(input,*) DECinfo%occ_init_radius
          DECinfo%occ_init_radius = DECinfo%occ_init_radius/bohr_to_angstrom
          DECinfo%all_init_radius = 0.0E0_realk
       case('.FRAG_INIT_RADIUS_NO_OPT_VIR')
          read(input,*) DECinfo%vir_init_radius
          DECinfo%vir_init_radius = DECinfo%vir_init_radius/bohr_to_angstrom
          DECinfo%all_init_radius = 0.0E0_realk

       case('.ATOMICEXTENT')
          !Include all atomic orbitals on atoms in the fragment 
          DECinfo%AtomicExtent  = .true.

       !KEYWORDS FOR CANONICAL FULL MOLECULAR MP2
       !**************************

       case('.MPMP2') 
          ! By default a memory conserving integral direct MP2 code is used
          ! when **CC .MP2 is combined with .CANONICAL
          ! This keyword activates a MP2 version that distribute the integrals
          ! across the nodes. The code does less integral recalculation but 
          ! requires more nodes/more memory
          DECinfo%MPMP2=.true.


       !KEYWORDS FOR RIMP2 (.RIMP2) 
       !**************************
       case('.AUXATOMICEXTENT')
          !Include all atomic orbitals on all atoms in the molecule (not just fragment) 
          !maybe need to have a procedure to optimize this set of atoms
          DECinfo%AuxAtomicExtent  = .true.
       case('.NAF')
          DECinfo%NAF                      = .true.
       case('.NAFTHRESHOLD')
          read(input,*) DECinfo%NAFthreshold
       case('.RIMP2SUBGROUPSIZE')
          read(input,*) DECinfo%RIMPSubGroupSize
       case('.RIMP2PDMTENSOR')
          DECinfo%RIMP2PDMTENSOR      = .true.
       case('.RIMP2FORCEPDMCALPHA')
          DECinfo%RIMP2ForcePDMCalpha = .true.
       case('.RIMP2_TILING')
          DECinfo%RIMP2_tiling        = .true.
       case('.RIMP2_CHOL')
          DECinfo%RIMP2_lowdin        = .false.
       case('.RIMP2_LAPLACE')
          DECinfo%RIMP2_Laplace       = .true.
       case('.RIMP2_NOOMP')
          DECinfo%RIMP2_deactivateopenmp = .true.

       !KEYWORDS FOR INTEGRAL INFO
       !**************************
       case('.INTEGRALTHRESHOLD')
          read(input,*) DECinfo%IntegralThreshold
          IF(DECinfo%IntegralThreshold.LT.shortintCRIT)THEN
             write(DECinfo%output,'(A)')'Error: you cannot chose integral threshold less then'
             write(DECinfo%output,'(ES15.6,A)') shortintCRIT, 'due to technical reasons'
             write(DECinfo%output,'(A)')'you could use .NO SCREEN (you may have to deactivate LinK with .NOLINK)'
             write(*,'(A)')'Error: you cannot chose integral threshold less then'
             write(*,'(ES15.6,A)') shortintCRIT, 'due to technical reasons'
             write(*,'(A)')'you could use .NO SCREEN (you may have to deactivate LinK with .NOLINK)'
             call lsquit('Error in choice of integral threshold',-1)
          ENDIF
       !Use the Ichor Integral Code (default is Thermite Code)   
       case('.ICHOR'); DECinfo%UseIchor = .true.


       ! MEMORY HANDLING KEYWORDS
       ! **************************
       case('.BACKGROUND_BUFFER');    DECinfo%use_bg_buffer           = .true.


#ifdef MOD_UNRELEASED
       ! CCSOLVER SPECIFIC KEYWORDS
       ! **************************
       case('.CCDRIVERDEBUG');        DECinfo%cc_driver_debug         = .true.
       case('.CCSOLVER_LOCAL');       DECinfo%solver_par              = .false.
       case('.CCSDPREVENTCANONICAL'); DECinfo%CCSDpreventcanonical    = .true.
       case('.SPAWN_COMM_PROC');      DECinfo%spawn_comm_proc         = .true.
       case('.CCSDNO_RESTART');       DECinfo%CCSDno_restart          = .true.
       case('.DIIS');                 DECinfo%use_crop                = .false.
       case('.CC_TILE_SIZE_GB')
          read(input,*) DECinfo%cc_solver_tile_mem 
       case('.NOTPREC')                 
          DECinfo%use_preconditioner=.false.
          DECinfo%ccsolver_overwrite_prec = .true.
       case('.NOTBPREC')
          DECinfo%use_preconditioner_in_b=.false.
          DECinfo%ccsolver_overwrite_prec = .true.
       case('.PRECWITHFULL')
          DECinfo%precondition_with_full=.true.
          DECinfo%ccsolver_overwrite_prec = .true.
       case('.CCSOLVERSKIP')
          DECinfo%ccsolverskip = .true.
       case('.MAXITER')
          read(input,*) DECinfo%MaxIter
       case('.TENSOR_SEGMENTING_SCHEME')
          read(input,*) DECinfo%tensor_segmenting_scheme


       ! CCSD RESIDUAL SPECIFIC KEYWORDS
       ! *******************************
       case('.CCSDDYNAMIC_LOAD');         DECinfo%dyn_load             = .true.
       case('.CCSDNODYNAMIC_LOAD');       DECinfo%dyn_load             = .false.
       case('.CCSDMULTIPLIERS');          DECinfo%CCSDmultipliers      = .true.
       case('.DEBUG_MULTIPLIERS_DIRECT'); DECinfo%simple_multipler_residual = .false.
       case('.NO_MO_CCSD');               DECinfo%NO_MO_CCSD           = .true.
       case('.CCSDEXPL');                 DECinfo%ccsd_expl            = .true.
#endif
       case('.CCSDFORCE_SCHEME');         DECinfo%force_scheme         = .true.
                                          read(input,*) DECinfo%en_mem
       case('.CCSD_DEBUG_COMMUNICATION'); DECinfo%CCSD_NO_DEBUG_COMM   = .false.

          ! Stripped-down keywords
          ! **********************
       case('.NOAOFOCK'); DECinfo%noaofock   = .true.


#ifdef MOD_UNRELEASED
       ! PNO-CCSD SPECIFIC KEYWORDS
       ! **************************
       case('.USE_PNOS');                 DECinfo%use_pnos             = .true.
       case('.PNO_DEBUG');                DECinfo%PNOtriangular        = .false.
       case('.PNOTHR');        read(input,*) DECinfo%simplePNOthr
       case('.NOPNOTRAFO');               DECinfo%noPNOtrafo           = .true.; DECinfo%noPNOtrunc=.true.
       case('.NOPNOTRUNCATION');          DECinfo%noPNOtrunc           = .true.
       case('.EOSPNOTHR');     read(input,*) DECinfo%EOSPNOthr
       case('.NOFATRAFO');                DECinfo%noFAtrafo            = .true.; DECinfo%noFAtrunc=.true.
       case('.NOFATRUNCATION');           DECinfo%noFAtrunc            = .true.
       case('.PNOOVERLAPTHR'); read(input,*) DECinfo%PNOoverlapthr
       case('.NOPNOOVERLAPTRUNCATION');   DECinfo%noPNOoverlaptrunc    = .true.
       case('.PNO_S_ON_THE_FLY');         DECinfo%pno_S_on_the_fly     = .true.
#endif

       ! KEYWORDS RELATED TO F12
       ! ***********************
       case('.F12')
          DECinfo%F12=.true.; doF12 = .TRUE.
       case('.F12FRAGOPT')     
          DECinfo%F12=.true.
          DECinfo%F12fragopt=.true.
          doF12 = .TRUE.
       case('.F12DEBUG')     
          DECinfo%F12=.true.
          DECinfo%F12DEBUG=.true.
          doF12 = .TRUE.
       case('.SKIPF12SINGLES')
          DECinfo%F12SINGLES=.false.
       case('.F12CCOUPLING')     
          DECinfo%F12Ccoupling=.true.
       case('.F12SINGLESMAXITER')
          read(input,*) DECinfo%F12singlesMaxIter
       case('.F12SINGLESTHR')
          read(input,*) DECinfo%F12singlesThr
       case('.F12SINGLESMAXDIIS')
          read(input,*) DECinfo%F12singlesMaxDIIS
       case('.F12LSALL')     
          !Use Natural linear scaling algorithm to treat 
          !these terms (do not treat with DEC nor RI)
          DECinfo%NaturalLinearScalingF12Terms   = .true.
          DECinfo%NaturalLinearScalingF12TermsB1 = .true.
          DECinfo%NaturalLinearScalingF12TermsX1 = .true.
          DECinfo%NaturalLinearScalingF12TermsV1 = .true.
       case('.F12LSB1')     
          DECinfo%NaturalLinearScalingF12Terms   = .true.
          DECinfo%NaturalLinearScalingF12TermsB1 = .true.
       case('.F12LSX1')     
          DECinfo%NaturalLinearScalingF12Terms   = .true.
          DECinfo%NaturalLinearScalingF12TermsX1 = .true.
       case('.F12LSV1')     
          DECinfo%NaturalLinearScalingF12Terms   = .true.
          DECinfo%NaturalLinearScalingF12TermsV1 = .true.

       ! KEYWORDS RELATED TO PAIR FRAGMENTS AND JOB LIST
       ! ***********************************************
       case('.PAIRTHR') 
          ! Threshold in a.u.
          read(input,*) DECinfo%pair_distance_threshold
       case('.PAIRTHRANGSTROM') 
          ! Input in Angstrom
          read(input,*) DECinfo%pair_distance_threshold
          DECinfo%pair_distance_threshold=DECinfo%pair_distance_threshold/bohr_to_angstrom
       case('.CHECKPAIRS') 
          DECinfo%checkpairs=.true.
       case('.PAIRMINDIST'); read(input,*) DECinfo%PairMinDist
       case('.PAIRFOTHR'); read(input,*) DECinfo%pairFOthr
       case('.NOTPAIRESTIMATE'); DECinfo%PairEstimate=.false.
       case('.IGNOREPAIRESTIMATE'); DECinfo%PairEstimateIgnore=.true.
       case('.ESTIMATEINITRADIUS')
          read(input,*) DECinfo%EstimateINITradius
          DECinfo%EstimateINITradius = DECinfo%EstimateINITradius/bohr_to_angstrom
       case('.ESTIMATEINITATOM')
          read(input,*) DECinfo%EstimateInitAtom
       case('.PAIRESTIMATEMODEL')
          read(input,*) myword
          call find_model_number_from_input(myword, DECinfo%PairEstimateModel)
       case('.PAIRMINDISTANGSTROM')
          read(input,*) DECinfo%PairMinDist
          DECinfo%PairMinDist = DECinfo%PairMinDist/bohr_to_angstrom
       case('.NFRAGSRED')
          ! Number of reduced space to consider in multi-FOT treatment of pairs
          read(input,*) DECinfo%nFRAGSred
       case('.FOTSCALING')
          ! Factor to scale FOT by for reduced pair fragments
          read(input,*) DECinfo%FOTscaling
       case('.NO_PAIRS'); DECinfo%No_Pairs = .true.
       case('.ONLY_N_JOBS')
          read(input,*)DECinfo%only_n_frag_jobs
          call mem_alloc(DECinfo%frag_job_nr,DECinfo%only_n_frag_jobs)
          read(input,*)DECinfo%frag_job_nr(1:DECinfo%only_n_frag_jobs)
       case('.ONLY_PAIR_FRAG_JOBS'); DECinfo%only_pair_frag_jobs = .true.


       ! FIRST ORDER PROPERTIES KEYWORD
       ! ******************************
       case('.UNRELAXDENSITY') 
          ! Calculate unrelaxed density
          DECinfo%unrelaxed =.true.
          DECinfo%density =.true.
          DECinfo%first_order=.true.

       case('.SKIPFULL') 
          !> Collect fragment contributions to calculate full molecular MP2 density
          DECinfo%SkipFull=.true.

       case('.ERRORFACTOR') 
          ! scaling factor for estimated error in DEC geometry optimization
          read(input,*) DECinfo%EerrFactor

       ! kappabar multiplier equation
       case('.KAPPAMAXITER'); read(input,*) DECinfo%kappaMaxIter 
       case('.KAPPAMAXDIIS'); read(input,*) DECinfo%kappaMaxDIIS
       case('.KAPPA_DEBUG'); DECinfo%kappa_driver_debug=.true.
       case('.NOTKAPPAPREC'); DECinfo%kappa_use_preconditioner=.false.
       case('.NOTKAPPABPREC'); DECinfo%kappa_use_preconditioner_in_b=.false.


       ! DEC ORBITAL TREATMENT
       ! *********************
       case('.READDECORBITALS'); DECinfo%read_dec_orbitals=.true.
       case('.ONLY_GENERATE_DECORBS'); DECinfo%only_generate_DECorbs=.true.
       case('.MULLIKEN'); DECinfo%mulliken=.true.
       case('.DISTANCE'); DECinfo%distance=.true.
       case('.NOTFITORBITALS'); DECinfo%FitOrbitals=.false.
       case('.SIMPLEORBITALTHRESH')
          read(input,*) DECinfo%simple_orbital_threshold
       case('.PURIFICATION'); DECinfo%PurifyMOs=.true.


       ! SINGLE POLARIZATION
       ! *******************
       case('.SINGLESPOLARI'); DECinfo%SinglesPolari=.true.
       case('.SINGLESTHR'); read(input,*) DECinfo%SinglesThr


       ! KEYWORDS RELATED TO I/O
       ! ***********************
       case('.CONVERT64TO32')
          DECinfo%convert64to32=.true.
       case('.CONVERT32TO64')
          DECinfo%convert32to64=.true.
       case('.ARRAY4ONFILE') 
          DECinfo%array4OnFile=.true.
          DECinfo%array4OnFile_specified=.true.
       case('.SKIPREADIN')
          ! Skip the read-in of molecular info files dens.restart, fock.restart, lcm_orbitals.u
          DECinfo%SkipReadIn=.true.
       case('.TIMEBACKUP'); read(input,*) DECinfo%TimeBackup

       ! KEYWORDS RELATED TO TENSOR HYPER CONTRACTION (THC)
       ! ***********************
       case('.THC_GRID'); 
          read(input,*) myword
          call capitalize_string(myword)
          DECinfo%THCradint = 1.0E-6_realk
          DECinfo%THCZdependenMaxAng=.FALSE.
          IF (INDEX(myword,'VCOARSE') .NE. 0) THEN 
             DECinfo%THC_MIN_RAD_PT = 0
             DECinfo%THCangint = 2
             DECinfo%THCNOPRUN = .FALSE.
          ELSEIF (INDEX(myword,'COARSE') .NE. 0) THEN 
             DECinfo%THC_MIN_RAD_PT = 3
             DECinfo%THCangint = 2
             DECinfo%THCNOPRUN = .FALSE.
          ELSEIF (INDEX(myword,'MEDIUM') .NE. 0) THEN 
             DECinfo%THC_MIN_RAD_PT = 8
             DECinfo%THCangint = 5
          ELSEIF (INDEX(myword,'VFINE') .NE. 0) THEN 
             DECinfo%THC_MIN_RAD_PT = 16
             DECinfo%THCangint = 9             
          ELSEIF (INDEX(myword,'FINE') .NE. 0) THEN 
             DECinfo%THC_MIN_RAD_PT = 12
             DECinfo%THCangint = 5
          ENDIF
       case('.THC_PRUNE'); DECinfo%THCNOPRUN = .FALSE.
       case('.THC_DUMP'); DECinfo%THCDUMP = .TRUE.
       case('.THC_RADINT'); read(input,*) DECinfo%THCradint 
       case('.THC_MIN_RAD_PT'); read(input,*) DECinfo%THC_MIN_RAD_PT
       case('.THC_ANGINT'); read(input,*) DECinfo%THCangint
       case('.THC_HRDNES'); read(input,*) DECinfo%THCHRDNES
       case('.THC_TURBO'); read(input,*) DECinfo%THCTURBO
       case('.THC_RADIALGRID'); read(input,*) DECinfo%THCRADIALGRID
       case('.THC_NOZDEPENDENTMAXANG'); DECinfo%THCZdependenMaxAng=.FALSE.
       case('.THC_PARTITIONING'); read(input,*) DECinfo%THCPARTITIONING

       CASE DEFAULT
          WRITE (output,'(/,3A,/)') ' Keyword "',WORD,&
               & '" not recognized in config_dec_input'
          CALL lsQUIT('Illegal keyword in config_dec_input',output)

       END SELECT DEC_INPUT_INFO

    ENDDO

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

    ! Check DECNP options compatibility:
    ! ----------------------------------
    if (DECinfo%DECNP) then
       write(DECinfo%output,*)
       write(DECinfo%output,*) "WARNING: DECNP calculate no pairs explicitly!"
       write(DECinfo%output,*) "--> All pair and pair estimate related keywords &
         &will be ignored"
       DECinfo%PairEstimate=.false.
       DECinfo%PairEstimateIgnore = .true.
       DECinfo%no_pairs = .true.
       ! MODIFY FOR NEW MODEL
       ! Some models are not compatible with DECNP
       select case(DECinfo%ccmodel) 
       case (MODEL_RPA)
          call lsquit("RPA model is not compatible with DECNP yet",DECinfo%output)
       case (MODEL_SOSEX)
          call lsquit("SOSEX model is not compatible with DECNP yet",DECinfo%output)
       case (MODEL_LSTHCRIMP2)
          call lsquit("LS-THC-RI-MP2 model is not compatible with DECNP yet",DECinfo%output)
       end select

       if (DECinfo%first_order) then
          call lsquit("No first_order properties with DECNP",DECinfo%output)
       end if

       ! DECNP and SNOOP are not compatible yet
       if (DECinfo%SNOOP) then
          call lsquit("SNOOP and DECNP are not compatible yet!",DECinfo%output)
       end if

    end if


    ! Repeat atomic fragment calcs after fragment optimization if:
    ! --------------------------------------------------------
    ! - First order properties are requested
    ! - Debug calculations: Include full Molecule, simulate full molecule:
    ! - MODIFY FOR NEW CORRECTION: A corection is requested in the target CC model (F12)
    ! - The model used in fragment reduction is different from the target CC model.
    if ( DECinfo%first_order .or. DECinfo%InclFullMolecule .or. DECinfo%simulate_full .or. &
       & DECinfo%F12 .or. (DECinfo%ccmodel/=DECinfo%fragopt_red_model )) then
       DECinfo%RepeatAF=.true.
    else
       DECinfo%RepeatAF=.false.
    end if

    
    ! Local masters print less information if we use more than 100 nodes
    if( nodtot > 100 .and. DECinfo%PL<=1)then
       DECinfo%print_small_calc = .false.
    endif

    ! Reduced pairs - certain limitations
    if(DECinfo%nFRAGSred>0) then
       if(DECinfo%fragadapt) then
          call lsquit('Reduced pairs not implemented for fragment-adapted DEC!',-1)
       end if
       if(DECinfo%use_pnos) then
          call lsquit('Reduced pairs not implemented for PNOs!',-1)
       end if
    end if

    ! SNOOP - currently limited in several ways
    if(DECinfo%SNOOP) then

       ! For DEC, we currently include all pairs for SNOOP
       write(DECinfo%output,*) 'WARNING: SNOOP currently requires all pairs to be calculated!'
       write(DECinfo%output,*) '--> ignoring pair estimates and setting pair distance threshold to be huge!'
       DECinfo%pair_distance_threshold = huge(1.0_realk)
       DECinfo%PairEstimate=.false.
       DECinfo%PairEstimateIgnore = .true.
       

       if(DECinfo%SNOOPlocalize .and. DECinfo%SNOOPsamespace) then
          call lsquit('SNOOP: Monomer orbitals cannot localized when subsystems &
               & use same orbital spaces as full system!',-1)
       end if
       
       ! Only for dense matrices for now
       if(matrix_type/=mtype_dense) then
          call lsquit('SNOOP is only implemented for dense matrices!',-1)
       end if

       ! SNOOP only tested for occupied partitioning scheme
       if(.not. DECinfo%OnlyOccPart) then
          write(DECinfo%output,*) 'WARNING: SNOOP ONLY TESTED FOR OCCUPIED PART. SCHEME'
          write(DECinfo%output,*) 'WARNING: I TURN ON OCCUPIED PART. SCHEME'
          DECinfo%onlyoccpart=.true.
       end if


       ! Not hydrogen debug
       if(decinfo%PureHydrogendebug) then
          call lsquit('SNOOP not implemented for hydrogen debug',-1)
       end if
       
       ! SimulateFull will destroy subsystem assignment and thus render SNOOP meaningless
       if(DECinfo%simulate_full) then
          call lsquit('SNOOP cannot be used in connection with the SIMULATEFULL keyword!',-1)
       end if

       ! Energy contribution analysis will not work for DECCO
       if(DECinfo%DECCO) then
          call lsquit('SNOOP cannot be used in connection with the DECCO keyword!',-1)
       end if

       if(Decinfo%distribute_fullmolecule)then
          print*,"WARNING: memory distribution for the molecule type in a snoop&
             & calculation is currently not implemented -> falling back to standart"
          Decinfo%distribute_fullmolecule = .false.
       endif
       
    end if


    ! CC response - currently not implemented for DEC
    CCresponse: if(DECinfo%CCexci) then
       IF(.not. DECinfo%full_molecular_cc) then
          call lsquit('CC response is not implemented for DEC! Use **CC instead of **DEC.',-1)
       end IF
       ! For now we enforce canonical orbitals for CC response
       if(.not. DECinfo%use_canonical) then
          write(DECinfo%output,*) 'WARNING! We enforce canonical orbitals for CC response!'
          DECinfo%use_canonical = .true.
       end if
       ! P_EOM_MBPT2 only for right transformation
       if(DECinfo%P_EOM_MBPT2 .and. DECinfo%JacobianLHTR) then
          call lsquit('P_EOM_MBPT2 only for Jacobian right transformation!',-1)
       end if
    end if CCresponse


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


    IF(DECinfo%full_molecular_cc)THEN
       IF(DECinfo%print_frags.AND.DECinfo%ccModel .EQ. MODEL_RIMP2)THEN
          call lsquit('A full molecular RIMP2 calculation do not construct the amplitudes and integrals. &
               & It is therefore not possible to print the fragment energies. &
               & Suggestion: Remove .PRINTFRAGS keyword!', DECinfo%output)
       ENDIF

       if(Decinfo%distribute_fullmolecule)then
          print*,"WARNING: memory distribution for the molecule type in a full&
          & calculation is currently not implemented -> falling back to standart"
          Decinfo%distribute_fullmolecule = .false.
       endif

    ENDIF

    FirstOrderModel: if(DECinfo%ccModel /= MODEL_MP2.and.DECinfo%ccModel /= MODEL_CCSD.and.DECinfo%ccModel /= MODEL_RIMP2) then

       if(DECinfo%density) then
          call lsquit('Calculation of density matrix is only implemented for MP2/CCSD!', DECinfo%output)
       end if

       if(DECinfo%gradient) then
          call lsquit('Calculation of molecular gradient is only implemented for MP2/CCSD!', DECinfo%output)
       end if

    end if FirstOrderModel

    MP2gradientCalculation: if(DECinfo%first_order) then

       if(DECinfo%full_molecular_cc.and.DECinfo%ccmodel==MODEL_MP2) then
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
    if(.not. DECinfo%CCthrSpecified .and. (.not. DECinfo%full_molecular_cc) ) then
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

       write(DECinfo%output,*) ''
       write(DECinfo%output,*) 'Memory not defined or ambiguously defined for **DEC or **CC calculation!'
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

       if(DECinfo%use_bg_buffer)then
          call lsquit("ERROR: bg buffer not implemented for fragment adapted orbitals",-1)
       endif

    end if

    if(DECinfo%use_system_memory_info) call get_currently_available_memory(DECinfo%memory)

    if(DECinfo%use_bg_buffer.AND.(DECinfo%bg_memory<0.0E0_realk)) then
       DECinfo%bg_memory = 0.8_realk*DECinfo%memory
       write(DECinfo%output,*) ''
       write(DECinfo%output,*) 'WARNING: User did not specify the amount of memory to be used'
       write(DECinfo%output,*) '         in connection with the background buffer.'
       write(DECinfo%output,*) ''
       write(DECinfo%output,*) 'By default, 80% of the total memory will be used:'
       write(DECinfo%output,'(A,F6.3,A)') ' Total memory             = ', DECinfo%memory,   ' GB'
       write(DECinfo%output,'(A,F6.3,A)') ' Background buffer memory = ', DECinfo%bg_memory,' GB'
       write(DECinfo%output,*) ''
       write(DECinfo%output,*) 'You can specify the amount of BG buffer memory yourself:'
#ifdef VAR_MPI
       write(DECinfo%output,*) 'E.g. if each MPI process has 16 GB of memory available, '
       write(DECinfo%output,*) 'and 8 GB should be used for the background buffer, then use'
#else
       write(DECinfo%output,*) 'E.g. if there are 16 GB of memory available, '
       write(DECinfo%output,*) 'and 8 GB should be used for the background buffer, then use'
#endif
       write(DECinfo%output,*) '.BG_MEMORY'
       write(DECinfo%output,*) '8.0'
    end if

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

    ! Meaningful RI split
    if(DECinfo%RIMPIsplit<1) then
       call lsquit('RIMPISPLIT must be larger than zero!',-1)
    end if
    

    ! Set FOTs for geometry opt.
    call set_geoopt_FOTs(DECinfo%FOT)

    if(DECinfo%noaofock) then
       if(DECinfo%SinglesPolari) then
          call lsquit('Singles polarization does not work with .NOAOFOCK keyword!',-1)
       end if
       if(DECinfo%use_canonical) then
          call lsquit('NOAOFOCK keyword does not work with canonical orbitals!',-1)
       end if       
       if(DECinfo%check_lcm_orbitals) then
          call lsquit('NOAOFOCK keyword does not work with CHECKLCM keyword!',-1)
       end if
       if(DECinfo%fragadapt) then
          call lsquit('NOAOFOCK keyword does not work with fragment-adapted orbitals!',-1)
       end if
       if(DECinfo%full_molecular_cc) then
          call lsquit('NOAOFOCK keyword does not work for full molecular calculation!',-1)
       end if
       ! The fock matrix is required in the ccsolver
       select case(DECinfo%ccmodel)
       case(MODEL_CC2,MODEL_CCSD,MODEL_CCSDpT,MODEL_RPA,MODEL_SOSEX)
          call lsquit("The CC solver require the fock matrix to be stored. Remove &
             & .NOAOFOCK keyword from input.",DECinfo%output)
       end select
    end if

    if ((.not.DECinfo%full_molecular_cc) .and. DECinfo%ccmodel==MODEL_RIMP2) then
       write(DECinfo%output,*) ''
       write(DECinfo%output,*) 'WARNING: User chose RI-MP2 as the final model,'
       write(DECinfo%output,*) '         we therefore enforce RI-MP2 to be used'
       write(DECinfo%output,*) '         also in the Fragment optimization and'
       write(DECinfo%output,*) '         Pair estimates calculations.'
       write(DECinfo%output,*) ''
       DECinfo%PairEstimateModel = MODEL_RIMP2
       DECinfo%fragopt_exp_model = MODEL_RIMP2
       DECinfo%fragopt_red_model = MODEL_RIMP2
    end if

    
    ! MP3 testing
    if(DECinfo%ccmodel==MODEL_MP3) then
       if(.not.DECinfo%full_molecular_cc) then
          call lsquit('MP3 only implemented for full molecular CC!',-1)
       end if
       if(DECinfo%first_order) then
          call lsquit('No first-order properties for MP3!',-1)
       end if
       if(.not. DECinfo%use_canonical) then
          call lsquit('MP3 only implemented for canonical orbitals, insert .CANONICAL!',-1)
       end if
    end if

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

    write(lupri,*) 'SNOOP ',DECinfo%SNOOP
    write(lupri,*) 'SNOOPjustHF ', DECinfo%SNOOPjustHF
    write(lupri,*) 'SNOOPMaxDIIS ', DECinfo%SNOOPMaxDIIS
    write(lupri,*) 'SNOOPMaxIter ', DECinfo%SNOOPMaxIter
    write(lupri,*) 'SNOOPthr ', DECinfo%SNOOPthr
    write(lupri,*) 'SNOOPdebug ', DECinfo%SNOOPdebug
    write(lupri,*) 'SNOOPsamespace ', DECinfo%SNOOPsamespace
    write(lupri,*) 'SNOOPlocalize ', DECinfo%SNOOPlocalize
    write(lupri,*) 'SNOOPrestart ', DECinfo%SNOOPrestart
    write(lupri,*) 'SNOOPonesub ', DECinfo%SNOOPonesub
    write(lupri,*) 'CCexci ', DECinfo%CCexci
    write(lupri,*) 'JacobianNumEival ', DECinfo%JacobianNumEival
    write(lupri,*) 'JacobianLHTR ', DECinfo%JacobianLHTR
    write(lupri,*) 'JacobianThr ', DECinfo%JacobianThr
    write(lupri,*) 'JacobianMaxSubspace ', DECinfo%JacobianMaxSubspace
    write(lupri,*) 'JacobianInitialSubspace ', DECinfo%JacobianInitialSubspace
    write(lupri,*) 'JacobianMaxIter ', DECinfo%JacobianMaxIter
    write(lupri,*) 'JacobianPrecond ', DECinfo%JacobianPrecond
    write(lupri,*) 'SinglesEW1 ', DECinfo%SinglesEW1
    write(lupri,*) 'LW1 ', DECinfo%LW1
    write(lupri,*) 'P_EOM_MBPT2 ', DECinfo%P_EOM_MBPT2
    write(lupri,*) 'doDEC ', DECitem%doDEC
    write(lupri,*) 'DECCO ', DECitem%DECCO
    write(lupri,*) 'DECNP ', DECitem%DECNP
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
    write(lupri,*) 'only_generate_DECorbs ', DECitem%only_generate_DECorbs
    write(lupri,*) 'IntegralThreshold ', DECitem%IntegralThreshold
    write(lupri,*) 'UseIchor ', DECitem%UseIchor
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
    write(lupri,*) 'CCSDno_restart ', DECitem%CCSDno_restart
    write(lupri,*) 'CCSDpreventcanonical ', DECitem%CCSDpreventcanonical
    write(lupri,*) 'CRASHCALC            ', DECitem%CRASHCALC
    write(lupri,*) 'CRASHESTI            ', DECitem%CRASHESTI
    write(lupri,*) 'cc_driver_debug ', DECitem%cc_driver_debug
    write(lupri,*) 'use_bg_buffer ', DECitem%use_bg_buffer
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
    write(lupri,*) 'F12 ', DECitem%F12
    write(lupri,*) 'F12DEBUG ', DECitem%F12DEBUG
    write(lupri,*) 'F12singles ', DECinfo%F12singles
    write(lupri,*) 'F12fragopt ', DECitem%F12fragopt
    write(lupri,*) 'F12CCOUPLING',DECinfo%F12Ccoupling
    write(lupri,*) 'mpisplit ', DECitem%mpisplit
    write(lupri,*) 'rimpisplit ', DECitem%rimpisplit
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
    write(lupri,*) 'MemDebugPrint ', DECitem%MemDebugPrint
    write(lupri,*) 'SkipFull ', DECitem%SkipFull
    write(lupri,*) 'output ', DECitem%output
    write(lupri,*) 'AbsorbHatoms ', DECitem%AbsorbHatoms
    write(lupri,*) 'FitOrbitals ', DECitem%FitOrbitals
    write(lupri,*) 'simple_orbital_threshold ', DECitem%simple_orbital_threshold
    write(lupri,*) 'PurifyMOs ', DECitem%PurifyMOs
    write(lupri,*) 'FragAdapt ', DECitem%FragAdapt
    write(lupri,*) 'mulliken ', DECitem%mulliken
    write(lupri,*) 'FOT ', DECitem%FOT
    write(lupri,*) 'MaxIter ', DECitem%MaxIter
    write(lupri,*) 'FOTlevel ', DECitem%FOTlevel
    write(lupri,*) 'Frag_Exp_Scheme ', DECitem%Frag_Exp_Scheme
    write(lupri,*) 'Frag_RedOcc_Scheme ', DECitem%Frag_RedOcc_Scheme
    write(lupri,*) 'Frag_RedVir_Scheme ', DECitem%Frag_RedVir_Scheme
    write(lupri,*) 'Frag_Init_Size ', DECitem%Frag_Init_Size
    write(lupri,*) 'Frag_Exp_Size ', DECitem%Frag_Exp_Size
    write(lupri,*) 'Frag_Red1_thr ', DECinfo%frag_red1_thr
    write(lupri,*) 'Frag_Red2_thr ', DECinfo%frag_red2_thr
    write(lupri,*) 'Frag_Red_Occ ', DECinfo%frag_red_occ
    write(lupri,*) 'Frag_Red_Virt ', DECinfo%frag_red_virt
    write(lupri,*) 'fragopt_exp_model ', DECitem%fragopt_exp_model
    write(lupri,*) 'fragopt_red_model ', DECitem%fragopt_red_model
    write(lupri,*) 'No_Pairs ', DECitem%no_pairs
    write(lupri,*) 'pair_distance_threshold ', DECitem%pair_distance_threshold
    write(lupri,*) 'PairMinDist ', DECitem%PairMinDist
    write(lupri,*) 'CheckPairs ', DECitem%CheckPairs
    write(lupri,*) 'pairFOthr ', DECitem%pairFOthr
    write(lupri,*) 'PairEstimate ', DECitem%PairEstimate
    write(lupri,*) 'first_order ', DECitem%first_order
    write(lupri,*) 'density ', DECitem%density
    write(lupri,*) 'unrelaxed ', DECitem%unrelaxed
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
    case('.SOSEX');   modelnumber = MODEL_SOSEX
    case('.RIMP2');   modelnumber = MODEL_RIMP2
    case('.LSTHCRIMP2'); modelnumber = MODEL_LSTHCRIMP2
    case('.MP3');     modelnumber = MODEL_MP3
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
       write(DECinfo%output,*)'.SOSEX'
       write(DECinfo%output,*)'.RIMP2'
       write(DECinfo%output,*)'.LS-THC-RIMP2'       
       write(DECinfo%output,*)'.MP3'       
       call lsquit('Requested model not found!',-1)
    end SELECT

  end subroutine find_model_number_from_input


end MODULE DEC_settings_mod
