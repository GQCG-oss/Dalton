!> @file 
!> Settings info and read of input for DEC/Dalton tests (DEC=Divide-Expand-Consolidate Coupled cluster)

!> \brief Settings info for simulated DEC tests.
!> \author \latexonly Kasper Kristensen  \endlatexonly
!> \date 2010-06-16
!>
MODULE DEC_settings_mod

  use fundamental
  use precision
  use dec_typedef_module
  use ls_util

contains

  !> \brief Set default DEC settings.
  !> See explanation of parameters in type DEC_settings.
  !> \author Kasper Kristensen
  !> \date June 2010
  subroutine dec_set_default_config(output)

    implicit none
    !> Unit number for DALTON.OUT
    integer, intent(in) :: output

    DECinfo%doDEC = .false.
    ! Max memory measured in GB. By default set to 2 GB
    DECinfo%memory=2.0E0_realk
    DECinfo%memory_defined=.false.
    DECinfo%frozencore=.false.
    DECinfo%ncalc = 0

    ! -- Type of calculation
    DECinfo%full_molecular_cc=.false. ! full molecular cc
    DECinfo%mp2energydebug=.false.
    DECinfo%simulate_full=.false.
    DECinfo%simulate_natoms=1
    DECinfo%SkipReadIn=.false.
    DECinfo%SinglesPolari=.false.
    DECinfo%SinglesThr=0.2E0_realk   ! this is completely random, currently under investigation
    DECinfo%convert64to32=.false.
    DECinfo%convert32to64=.false.
    DECinfo%restart = .false.
    DECinfo%TimeBackup = 300.0E0_realk   ! backup every 5th minute
    DECinfo%read_dec_orbitals = .false.
    DECinfo%CheckPairs=.false.
    call dec_set_model_names(DECinfo)


    ! -- Debug modes
    DECinfo%cc_driver_debug=.false.
    DECinfo%ccsd_old=.false.
    DECinfo%manual_batchsizes=.false.
    DECinfo%ccsdAbatch=0
    DECinfo%ccsdGbatch=0
    DECinfo%hack=.false.
    DECinfo%hack2=.false.
    DECinfo%mpisplit=10
    DECinfo%dyn_load=.false.
    DECinfo%force_scheme=.false.
    DECinfo%en_mem=0
    DECinfo%array_test=.false.
    DECinfo%reorder_test=.false.
    DECinfo%CCSDno_restart=.false.
    DECinfo%CCSDsaferun=.false.
    DECinfo%solver_par=.false.
    DECinfo%CCSDpreventcanonical=.false.
    DECinfo%CCSD_MPICH=.false.
    DECinfo%CCDhack = .false.

    ! -- Output options 
    DECinfo%output=output


    ! -- Orbital
    DECinfo%mulliken_threshold=0.01
    DECinfo%simple_mulliken_threshold=.false.
    DECinfo%approximated_norm_threshold=0.1E0_realk
    DECinfo%check_lcm_orbitals=.false.
    DECinfo%use_canonical=.false.
    DECinfo%user_defined_orbitals=.false.
    DECinfo%AbsorbHatoms=.true.  ! reassign H atoms to heavy atom neighbour
    DECinfo%mulliken=.false.
    DECinfo%Distance=.false.
    DECinfo%BoughtonPulay=.false.
    DECinfo%FitOrbitals=.true.
    DECinfo%simple_orbital_threshold=0.05E0_realk
    DECinfo%simple_orbital_threshold_set=.false.


    ! -- Fragment
    DECinfo%MaxIter=20
    DECinfo%FOTlevel=4
    DECinfo%maxFOTlevel=8   ! if you modify this remember to modify dimension of ncalc as well!
    DECinfo%FOT=1.0E-4_realk
    DECinfo%InclFullMolecule = .false.
    DECinfo%PL=0
    DECinfo%PurifyMOs=.false.
    DECinfo%precondition_with_full=.false.
    DECinfo%HybridScheme=.false.
    DECinfo%FragmentExpansionSize = 5
    DECinfo%fragadapt=.false.
    ! for ccsd(t) calculations, option to use MP2 optimized fragments
    DECinfo%use_mp2_frag=.true.

    ! -- Pair fragments
    DECinfo%pair_distance_threshold=10.0E0_realk/bohr_to_angstrom
    DECinfo%paircut_set=.false.
    ! Pair reduction distance - default 5 Angstrom.
    DECinfo%PairReductionDistance = 5.0E0_realk/bohr_to_angstrom
    DECinfo%PairMinDist = 3.0E0_realk/bohr_to_angstrom  ! 3 Angstrom

    ! Memory use for full molecule structure
    DECinfo%fullmolecule_memory=0E0_realk


    ! -- CC solver options

    DECinfo%ccsd_expl=.false.
    DECinfo%simulate_eri= .false.
    DECinfo%fock_with_ri= .false.
    DECinfo%ccMaxIter=100
    DECinfo%ccMaxDIIS=3
    DECinfo%ccModel=1 ! 1 - MP2, 2 - CC2, 3 - CCSD, 4 - CCSD(T), 5 - RPA
    DECinfo%F12=.false.
    DECinfo%ccConvergenceThreshold=1e-5
    DECinfo%CCthrSpecified=.false.
    DECinfo%use_singles=.false.
    DECinfo%use_preconditioner=.true.
    DECinfo%use_preconditioner_in_b=.true.
    DECinfo%use_crop=.true.
    DECinfo%array4OnFile=.false.
    DECinfo%array4OnFile_specified=.false.


    ! First order properties
    DECinfo%first_order = .false.

    !> MP2 density matrix   
    DECinfo%MP2density = .false.
    DECinfo%SkipFull = .false.

    !-- MP2 gradient
    DECinfo%gradient=.false.
    DECinfo%kappa_use_preconditioner=.true.
    DECinfo%kappa_use_preconditioner_in_b=.true.
    DECinfo%kappaMaxDIIS=3
    DECinfo%kappaMaxIter=100
    DECinfo%kappa_driver_debug=.false.
    DECinfo%kappaTHR=1e-4
    DECinfo%EerrFactor = 1.0_realk
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

    DECitem%cc_models(1)='MP2     '
    DECitem%cc_models(2)='CC2     '
    DECitem%cc_models(3)='CCSD    '
    DECitem%cc_models(4)='CCSD(T) '
    DECitem%cc_models(5)='RPA     '

  end subroutine dec_set_model_names


  !> \brief Read the **DEC or **CC input section in LSDALTON.INP and set 
  !> configuration structure accordingly.
  !> \author Kasper Kristensen
  !> \date September 2010
  SUBROUTINE config_dec_input(input,output,readword,word,fullcalc)
    implicit none
    !> Logical for keeping track of when to read
    LOGICAL,intent(inout)                :: READWORD
    !> Logical unit number for LSDALTON.INP
    integer,intent(in) :: input
    !> Logical unit number for DALTON.OUT
    integer,intent(in) :: output
    !> Word read from input
    character(len=80),intent(inout) :: word
    !> Is this a full calculation (fullcalc=true, input **CC) 
    !> or a DEC calculation (fullcalc=false, input=**DEC)
    logical,intent(in) :: fullcalc
    logical,save :: already_called = .false.
    integer :: fotlevel

    ! Sanity check that this routine is only called once for either **DEC OR **CC
    if(already_called) then
       call lsquit('Error: LSDALTON.INP must contain EITHER **DEC for DEC calculation OR &
            & **CC for full molecular calculation!',-1)
    else
       ! First call to this routine
       already_called=.true.
    end if

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

          ! CC model
       case('.MP2'); DECinfo%ccModel=1; DECinfo%use_singles=.false.  ! both DEC and full calc
       case('.CC2'); DECinfo%ccModel=2; DECinfo%use_singles=.true.   ! only for full calc
       case('.CCSD'); DECinfo%ccModel=3; DECinfo%use_singles=.true.; DECinfo%solver_par=.true.  ! only for full calc


          ! CC SOLVER INFO
          ! ==============

          ! Save CCSD amplitudes to be able to restart full CCSD calculation
       case('.CCSDSAFE'); DECinfo%CCSDsaferun=.true.

          ! Maximum number of CC iterations
       case('.CCMAXITER'); read(input,*) DECinfo%ccMaxIter 

          ! Residual norm threshold for CC amplitude equation
       case('.CCTHR') 
          read(input,*) DECinfo%ccConvergenceThreshold
          DECinfo%CCthrSpecified=.true.
          ! Number of residual vectors to save when solving CC amplitude equation
       case('.SUBSIZE'); read(input,*) DECinfo%ccMaxDIIS


          ! CHOICE OF ORBITALS
          ! ==================
          ! By default canonical orbitals are used for full molecular calculation, 
          ! while local orbitals are used for DEC calculation.
          ! These default choices can be overruled by these keywords:

          ! Use canonical orbitals 
       case('.CANONICAL') 
          DECinfo%use_canonical=.true.
          DECinfo%user_defined_orbitals=.true.

          ! Do not use canonical orbitals
       case('.NOTCANONICAL') 
          DECinfo%use_canonical=.false.
          DECinfo%user_defined_orbitals=.true.



          ! DEC CALCULATION 
          ! ===============

          ! Restart DEC calculation (only for single point calculations, not geometry optimization)
       case('.RESTART'); DECinfo%restart=.true.           

          ! Do not absorb H atoms when assigning orbitals
       case('.NOTABSORBH'); DECinfo%AbsorbHatoms=.false.

          !> See description of FOT level in set_input_for_fot_level.
          !> Note that if one does not use simple orbital threshold, then
          !> one should manually choose a suitable threshold to determine
          !> the atomic extent (e.g. threshold for Boughton-Pulay procedure).
       case('.FOT')  
          read(input,*) FOTlevel
          call set_input_for_fot_level(FOTlevel)

          ! Use frozen core approximation
       case('.FROZENCORE') 
          DECinfo%frozencore=.true.

          ! Max memory available on node measured in GB. By default set to 16 GB
       case('.MEMORY') 
          read(input,*) DECinfo%memory           
          DECinfo%memory_defined=.true.

          ! Pair distance threshold
       case('.PAIRTHR') 
          ! Threshold in a.u.
          read(input,*) DECinfo%pair_distance_threshold
          DECinfo%paircut_set=.true.  ! overwrite default pair cutoff defined by .FOT
       case('.PAIRTHRANGSTROM') 
          ! Input in Angstrom
          read(input,*) DECinfo%pair_distance_threshold
          DECinfo%pair_distance_threshold=DECinfo%pair_distance_threshold/bohr_to_angstrom
          DECinfo%paircut_set=.true.  ! overwrite default pair cutoff defined by .FOT

          ! Calculate MP2 gradient (without necessarily doing geometry optimization)
       case('.GRADIENT') 
          DECinfo%gradient=.true.
          DECinfo%first_order=.true.

          !> Carry out MP2 density calculation (subset of gradient calculation)
       case('.DENSITY') 
          DECinfo%MP2density=.true.
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
       case('.TESTARRAY'); DECinfo%array_test=.true.
       case('.TESTREORDERINGS'); DECinfo%reorder_test=.true.
    
       !CCSD testing
       case('.CCSDFORCE_SCHEME'); DECinfo%force_scheme=.true.
          read(input,*) DECinfo%en_mem

       case('.MANUAL_BATCHSIZES') 
          DECinfo%manual_batchsizes=.true.
          read(input,*) DECinfo%ccsdAbatch, DECinfo%ccsdGbatch
       case('.MPISPLIT'); read(input,*) DECinfo%MPIsplit
       case('.INCLUDEFULLMOLECULE');DECinfo%InclFullMolecule=.true.
          ! Size of local groups in MPI scheme
       case('.MPIGROUPSIZE'); read(input,*) DECinfo%MPIgroupsize

#ifdef MOD_UNRELEASED
       case('.CCSD_OLD'); DECinfo%ccsd_old=.true.
       case('.CCSDSOLVER_SERIAL'); DECinfo%solver_par=.false.
       case('.CCSDDYNAMIC_LOAD'); DECinfo%dyn_load=.true.
       case('.CCSDNO_RESTART'); DECinfo%CCSDno_restart=.true.
       case('.CCSD_WITH_MPICH'); DECinfo%CCSD_MPICH=.true.
       case('.CCSDPREVENTCANONICAL'); DECinfo%CCSDpreventcanonical=.true.
       case('.CCD'); DECinfo%CCDhack=.true.;DECinfo%ccModel=3; DECinfo%use_singles=.true.; DECinfo%solver_par=.true.
       case('.HACK'); DECinfo%hack=.true.
       case('.HACK2'); DECinfo%hack2=.true.
       case('.TIMEBACKUP'); read(input,*) DECinfo%TimeBackup
       case('.READDECORBITALS'); DECinfo%read_dec_orbitals=.true.
       case('.CCSD(T)'); DECinfo%ccModel=4; DECinfo%use_singles=.true.; DECinfo%solver_par=.true.
       case('.RPA'); DECinfo%ccModel=5; DECinfo%use_singles=.false.
       case('.NOTUSEMP2FRAG') 
          DECinfo%use_mp2_frag=.false.
          !
       case('.F12'); DECinfo%F12=.true.
       case('.NOTPREC'); DECinfo%use_preconditioner=.false.
       case('.NOTBPREC'); DECinfo%use_preconditioner_in_b=.false.
       case('.MULLIKEN'); DECinfo%mulliken=.true.
       case('.DISTANCE'); DECinfo%distance=.true.
       case('.BOUGHTONPULAY'); DECinfo%BoughtonPulay=.true.
       case('.NOTFITORBITALS'); DECinfo%FitOrbitals=.false.
       case('.SIMPLEORBITALTHRESH')
          read(input,*) DECinfo%simple_orbital_threshold
          DECinfo%simple_orbital_threshold_set=.true.
       case('.DIIS'); DECinfo%use_crop=.false.  ! use DIIS instead of CROP
       case('.MAXITER'); read(input,*) DECinfo%MaxIter
       case('.DECPRINT'); read(input,*) DECinfo%PL
       case('.MULLIKENTHR'); read(input,*) DECinfo%mulliken_threshold
       case('.SKIPPAIRS') 
          DECinfo%pair_distance_threshold=0.0E0_realk
          DECinfo%paircut_set=.true.  ! overwrite default pair cutoff defined by .FOT
       case('.CHECKPAIRS') 
          DECinfo%checkpairs=.true.
       case('.PAIRREDDIST') 
          read(input,*) DECinfo%PairReductionDistance 
       case('.PAIRREDDISTANGSTROM') 
          read(input,*) DECinfo%PairReductionDistance 
          DECinfo%PairReductionDistance = DECinfo%PairReductionDistance/bohr_to_angstrom
       case('.PAIRMINDIST'); read(input,*) DECinfo%PairMinDist
       case('.PAIRMINDISTANGSTROM')
          read(input,*) DECinfo%PairMinDist
          DECinfo%PairMinDist = DECinfo%PairMinDist/bohr_to_angstrom
       case('.CCSDEXPL'); DECinfo%ccsd_expl=.true.
       case('.PURIFICATION'); DECinfo%PurifyMOs=.true.
       case('.PRECWITHFULL'); DECinfo%precondition_with_full=.true.
       case('.SIMPLEMULLIKENTHRESH'); DECinfo%simple_mulliken_threshold=.true.
       case('.NORMTHRESH'); read(input,*) DECinfo%approximated_norm_threshold
       case('.MP2DEBUG'); DECinfo%mp2energydebug=.true.
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
       case('.FRAGMENTEXPANSIONSIZE'); read(input,*) DECinfo%FragmentExpansionSize
       case('.FRAGMENTADAPTED'); DECinfo%fragadapt=.true.

          ! kappabar multiplier equation
       case('.KAPPAMAXITER'); read(input,*) DECinfo%kappaMaxIter 
       case('.KAPPAMAXDIIS'); read(input,*) DECinfo%kappaMaxDIIS
       case('.KAPPA_DEBUG'); DECinfo%kappa_driver_debug=.true.
       case('.NOTKAPPAPREC'); DECinfo%kappa_use_preconditioner=.false.
       case('.NOTKAPPABPREC'); DECinfo%kappa_use_preconditioner_in_b=.false.

          ! Check that input orbitals are orthogonal (debug)
       case('.CHECKLCM'); DECinfo%check_lcm_orbitals=.true.

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

    ! Check that array4OnFile is only called for the cases where it is implemented

    ArraysOnFile: if(DECinfo%array4OnFile) then

       ! Only for MP2 so far
       if(DECinfo%ccModel /= 1) then
          call lsquit('Storing arrays on file only implemented for MP2. &
               & Suggestion: Remove .array4OnFile keyword!', DECinfo%output)
       end if

       if(DECinfo%use_singles) then
          call lsquit('Storing arrays on file not implemented for singles. &
               & Suggestion: Remove .array4OnFile keyword!', DECinfo%output)
       end if


    end if ArraysOnFile


    BeyondMp2: if(DECinfo%ccModel /= 1) then


       if(DECinfo%MP2density) then
          call lsquit('Calculation of density matrix is only implemented for MP2!', DECinfo%output)
       end if

       if(DECinfo%gradient) then
          call lsquit('Calculation of molecular gradient is only implemented for MP2!', DECinfo%output)
       end if

       ! Turn on the occupied/virtual hybrid scheme
       DECinfo%HybridScheme=.true.

    end if BeyondMp2


    MP2gradientCalculation: if(DECinfo%first_order) then

       if(DECinfo%full_molecular_cc) then
          call lsquit('Full calculation for MP2 gradient is implemented via the &
               & .SimulateFull keyword', DECinfo%output)
       end if

    end if MP2gradientCalculation


    ! If simulate full calculation, the full molecule must be included in "fragments"
    SimulateFullCalc: if(DECinfo%simulate_full) then
       DECinfo%InclFullMolecule = .true.
    end if SimulateFullCalc


    if(DECinfo%ccmodel==4 .and. DECinfo%restart .and. (.not. DECinfo%use_mp2_frag)) then
       call lsquit('Restart option currently not implemented for CCSD(T)!',DECinfo%output)
    end if

    ! Which orbitals to use?
    ! Default full calculation: Canonical orbitals
    ! Default DEC calculation : Local orbitals
    ! --> unless user has manually specified otherwise!
    if(.not. DECinfo%user_defined_orbitals) then  
       if(DECinfo%full_molecular_cc) then
          DECinfo%use_canonical = .true.
       else
          DECinfo%use_canonical = .false.
       end if
    end if

    ! For special case of full MP2, we simply only accept canonical orbitals!
    if(DECinfo%user_defined_orbitals .and. DECinfo%ccmodel==1 .and. DECinfo%full_molecular_cc) then
       write(DECinfo%output,*) 'WARNING! You have requested a full molecular MP2 calculation using'
       write(DECinfo%output,*) 'local orbitals. This option is currently not available so I will'
       write(DECinfo%output,*) 'use canonical orbitals instead!'
       DECinfo%use_canonical = .true.
    end if


    ! Set CC residual threshold to be 0.01*FOT
    ! - unless it was specified explicitly in the input.
    if(.not. DECinfo%CCthrSpecified) then
       DECinfo%ccConvergenceThreshold=0.01E0_realk*DECinfo%FOT
    end if

    ! Only full molecular for RPA at this stage
    if(DECinfo%ccmodel==5 .and. .not. DECinfo%full_molecular_cc) then
       call lsquit('RPA only implemented for full molecule! Use **CC rather than **DEC.',-1)
    end if

    ! Never use gradient and density at the same time (density is a subset of gradient)
    if(DECinfo%MP2density .and. DECinfo%gradient) then
       call lsquit('Density and gradient cannot both be turned on at the same time! &
            & Note that density is a subset of a gradient calculation',DECinfo%output)
    end if

    if(DECinfo%SinglesPolari) then
       call lsquit('Full singles polarization has been temporarily disabled!',-1)
    end if


    ! FOs do not work with reduced pairs, set reduction distance to 1000000 to
    ! avoid it from being used in practice
    ! Also use purification of MOs.
    if(DECinfo%fragadapt) then
       DECinfo%PairReductionDistance = 1.0e6_realk
       DECinfo%purifyMOs=.true.
    end if

#ifdef RELEASE
if(.not. DECinfo%full_molecular_cc .and. DECinfo%ccmodel/=1) then
   call lsquit('Error in input: DEC scheme only implemented for MP2 model!',-1)
end if
#endif

  end subroutine check_dec_input
  
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
    write(lupri,*) 'user_defined_orbitals ', DECitem%user_defined_orbitals
    write(lupri,*) 'simulate_full ', DECitem%simulate_full
    write(lupri,*) 'simulate_natoms ', DECitem%simulate_natoms
    write(lupri,*) 'InclFullMolecule ', DECitem%InclFullMolecule
    write(lupri,*) 'cc_models ', DECitem%cc_models
    write(lupri,*) 'ccModel ', DECitem%ccModel
    write(lupri,*) 'use_singles ', DECitem%use_singles
    write(lupri,*) 'restart ', DECitem%restart
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
    write(lupri,*) 'CCSDsaferun ', DECitem%CCSDsaferun
    write(lupri,*) 'solver_par ', DECitem%solver_par
    write(lupri,*) 'force_scheme ', DECitem%force_scheme
    write(lupri,*) 'dyn_load ', DECitem%dyn_load
    write(lupri,*) 'ccsd_old ', DECitem%ccsd_old
    write(lupri,*) 'CCSDno_restart ', DECitem%CCSDno_restart
    write(lupri,*) 'CCSDpreventcanonical ', DECitem%CCSDpreventcanonical
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
    write(lupri,*) 'simulate_eri ', DECitem%simulate_eri
    write(lupri,*) 'fock_with_ri ', DECitem%fock_with_ri
    write(lupri,*) 'F12 ', DECitem%F12
    write(lupri,*) 'mpisplit ', DECitem%mpisplit
    write(lupri,*) 'MPIgroupsize ', DECitem%MPIgroupsize
    write(lupri,*) 'manual_batchsizes ', DECitem%manual_batchsizes
    write(lupri,*) 'ccsdAbatch,ccsdGbatch ', DECitem%ccsdAbatch,DECitem%ccsdGbatch
    write(lupri,*) 'hack ', DECitem%hack
    write(lupri,*) 'hack2 ', DECitem%hack2
    write(lupri,*) 'mp2energydebug ', DECitem%mp2energydebug
    write(lupri,*) 'SkipReadIn ', DECitem%SkipReadIn
    write(lupri,*) 'array_test ', DECitem%array_test
    write(lupri,*) 'reorder_test ', DECitem%reorder_test
    write(lupri,*) 'check_lcm_orbitals ', DECitem%check_lcm_orbitals
    write(lupri,*) 'PL ', DECitem%PL
    write(lupri,*) 'SkipFull ', DECitem%SkipFull
    write(lupri,*) 'output ', DECitem%output
    write(lupri,*) 'AbsorbHatoms ', DECitem%AbsorbHatoms
    write(lupri,*) 'FitOrbitals ', DECitem%FitOrbitals
    write(lupri,*) 'simple_orbital_threshold ', DECitem%simple_orbital_threshold
    write(lupri,*) 'PurifyMOs ', DECitem%PurifyMOs
    write(lupri,*) 'FragAdapt ', DECitem%FragAdapt
    write(lupri,*) 'simple_orbital_threshold_set ', DECitem%simple_orbital_threshold_set
    write(lupri,*) 'BoughtonPulay ', DECitem%BoughtonPulay
    write(lupri,*) 'mulliken_threshold ', DECitem%mulliken_threshold
    write(lupri,*) 'simple_mulliken_threshold ', DECitem%simple_mulliken_threshold
    write(lupri,*) 'approximated_norm_threshold ', DECitem%approximated_norm_threshold
    write(lupri,*) 'mulliken ', DECitem%mulliken
    write(lupri,*) 'FOT ', DECitem%FOT
    write(lupri,*) 'MaxIter ', DECitem%MaxIter
    write(lupri,*) 'FOTlevel ', DECitem%FOTlevel
    write(lupri,*) 'maxFOTlevel ', DECitem%maxFOTlevel
    write(lupri,*) 'HybridScheme ', DECitem%HybridScheme
    write(lupri,*) 'FragmentExpansionSize ', DECitem%FragmentExpansionSize
    write(lupri,*) 'use_mp2_frag ', DECitem%use_mp2_frag
    write(lupri,*) 'pair_distance_threshold ', DECitem%pair_distance_threshold
    write(lupri,*) 'paircut_set ', DECitem%paircut_set
    write(lupri,*) 'PairReductionDistance ', DECitem%PairReductionDistance
    write(lupri,*) 'PairMinDist ', DECitem%PairMinDist
    write(lupri,*) 'CheckPairs ', DECitem%CheckPairs
    write(lupri,*) 'first_order ', DECitem%first_order
    write(lupri,*) 'MP2density ', DECitem%MP2density
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



  !> \brief Set DEC parameters in DEC structure according to FOT level,
  !> e.g. FOT itself, pair cut-off, and orbital extent threshold are set here.
  !> See details inside subroutine.
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine set_input_for_fot_level(FOTlevel)
    implicit none
    !> FOT level defining precision of calculation
    integer,intent(in) :: FOTlevel

    !> FOTlevel: Defines the precision of the whole calculation:
    !>           In general FOT=10^{-FOTlevel}
    !>           So a large FOTlevel means HIGH precision.
    !>           Reasonable pair cutoff distances are associated with each FOT level.
    !>           If necessary, the paircut off is increased
    !>           in a self-adaptive black box manner during
    !>           the calculation.
    !>           It is also possible to adapt the orbital threshold to the given FOT,
    !>           however, for now we simply set the orbital threshold to 0.01 for all levels.
    !> 
    !> FOTlevel = 1: 
    !> FOT=10^{-1}, pair_distance_threshold=4Angstrom.
    !>
    !> FOTlevel = 2: 
    !> FOT=10^{-2}, pair_distance_threshold=6Angstrom.
    !> 
    !> FOTlevel = 3: 
    !> FOT=10^{-3}, pair_distance_threshold=8Angstrom.
    !> 
    !> FOTlevel = 4: 
    !> FOT=10^{-4}, pair_distance_threshold=10Angstrom.
    !>
    !> FOTlevel = 5: 
    !> FOT=10^{-5}, pair_distance_threshold=12Angstrom.
    !> 
    !> FOTlevel = 6: 
    !> FOT=10^{-6}, pair_distance_threshold=14Angstrom.
    !>
    !> FOTlevel = 7: 
    !> FOT=10^{-7}, pair_distance_threshold=16Angstrom.
    !> 
    !> FOTlevel = 8: 
    !> FOT=10^{-8}, pair_distance_threshold=18Angstrom.
    !>
    !> Default: FOTlevel=4. If FOTlevel is not 1,2,3,4,5,6,7, or 8, the program will quit.

    ! Set FOT level
    DECinfo%FOTlevel = FOTlevel

    WhatFOTlevel: SELECT CASE(DECinfo%FOTlevel)

    case(1)
       DECinfo%FOT = 1.0E-1_realk
       ! Use FOT-adapted pair cutoff and simple orbital threshold
       ! only if these were not set set manually
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=4.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.05E0_realk
       end if

    case(2)
       DECinfo%FOT = 1.0E-2_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=6.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.05E0_realk
       end if

    case(3)
       DECinfo%FOT = 1.0E-3_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=8.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.05E0_realk
       end if

    case(4)
       DECinfo%FOT = 1.0E-4_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=10.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.05E0_realk
       end if

    case(5)
       DECinfo%FOT = 1.0E-5_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=12.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.05E0_realk
       end if

    case(6)
       DECinfo%FOT = 1.0E-6_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=14.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.05E0_realk
       end if

    case(7)
       DECinfo%FOT = 1.0E-7_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=16.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.05E0_realk
       end if

    case(8)
       DECinfo%FOT = 1.0E-8_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=18.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.05E0_realk
       end if

    case default
       call lsquit('set_input_for_fot_level: FOT level must be 1,2,3,4,5,6,7, or 8!',DECinfo%output)

    end SELECT WhatFOTlevel


  end subroutine set_input_for_fot_level

end MODULE DEC_settings_mod
