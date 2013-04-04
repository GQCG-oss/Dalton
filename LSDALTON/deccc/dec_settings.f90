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

contains

  !> \brief Set default settings for DEC test cases (default: run no test cases)
  !> See explanation of parameters in type DEC_settings.
  !> \author Kasper Kristensen
  !> \date June 2010
  subroutine dec_set_default_config(output)

    implicit none
    !> Unit number for DALTON.OUT
    integer, intent(in) :: output

    DECinfo%doDEC = .false.
    DECinfo%doHF = .true.
    ! Max memory measured in GB. By default set to 2 GB
    DECinfo%memory=2.0E0_realk
    DECinfo%memory_defined=.false.
    DECinfo%frozencore=.false.
    DECinfo%ncalc = 0

    ! -- Type of calculation
    DECinfo%full_molecular_cc=.false. ! full molecular cc
    DECinfo%FullDEC=.false.
    DECinfo%simulate_full=.false.
    DECinfo%simulate_natoms=1
    DECinfo%SkipReadIn=.false.
    DECinfo%SinglesPolari=.false.
    DECinfo%SinglesThr=0.2E0_realk   ! this is completely random, currently under investigation
    DECinfo%nsubfrag=0
    DECinfo%nsubfrag2=0
    DECinfo%subfrag=.false.
    DECinfo%convert64to32=.false.
    DECinfo%convert32to64=.false.
    DECinfo%restart = .false.
    DECinfo%TimeBackup = 300.0E0_realk   ! backup every 5th minute
    DECinfo%read_dec_orbitals = .false.
    DECinfo%NoExtraPairs=.false.


    ! -- Debug modes
    DECinfo%fragmentation_debug=.false. 
    DECinfo%dec_driver_debug=.false.
    DECinfo%cc_driver_debug=.false.
    DECinfo%ccsd_old=.false.
    DECinfo%manual_batchsizes=.false.
    DECinfo%ccsdAbatch=0
    DECinfo%ccsdGbatch=0
    DECinfo%hack=.false.
    DECinfo%hack2=.false.
    DECinfo%mpidebug=.false.
    DECinfo%mpisplit=10
    DECinfo%dyn_load=.false.
    DECinfo%force_scheme=.false.
    DECinfo%en_mem=0
    DECinfo%array_test=.false.
    DECinfo%reorder_test=.false.
    DECinfo%CCSDno_restart=.false.
    DECinfo%CCSDsaferun=.false.

    ! -- Output options 
    DECinfo%output=output


    ! -- Orbital
    DECinfo%mulliken_threshold=0.01
    DECinfo%simple_mulliken_threshold=.false.
    DECinfo%approximated_norm_threshold=0.1E0_realk
    DECinfo%check_lcm_orbitals=.false.
    DECinfo%use_canonical=.false.
    DECinfo%reassignHatoms=.false.  ! reassign H atoms to heavy atom neighbour
    DECinfo%mulliken=.false.
    DECinfo%BoughtonPulay=.false.
    DECinfo%FitOrbitals=.true.
    DECinfo%simple_orbital_threshold=0.01E0_realk
    DECinfo%simple_orbital_threshold_set=.false.


    ! -- Fragment
    DECinfo%MaxIter=20
    DECinfo%FOTlevel=4
    DECinfo%maxFOTlevel=7
    DECinfo%FOT=1.0E-4_realk
    DECinfo%InclFullMolecule = .false.
    DECinfo%PL=0
    DECinfo%SkipCC=.false.
    DECinfo%NormalizeFragment=.false.
    DECinfo%precondition_with_full=.false.
    DECinfo%HybridScheme=.false.
    DECinfo%LagStepSize = 5
    DECinfo%AEstep = 3
    DECinfo%fragadapt=.false.
    ! for ccsd(t) calculations, option to use ccsd optimized fragments
    DECinfo%use_ccsd_frag=.false.

    ! -- Pair fragments
    DECinfo%pair_distance_threshold=10.0E0_realk/bohr_to_angstrom
    DECinfo%paircut_set=.false.
    ! Pair reduction distance - default 1000 Angstrom, effectively turned off.
    DECinfo%PairReductionDistance = 1000.0E0_realk/bohr_to_angstrom
    DECinfo%PairMinDist = 3.0E0_realk/bohr_to_angstrom  ! 3 Angstrom

    ! Memory use for full molecule structure
    DECinfo%fullmolecule_memory=0E0_realk


    ! -- CC solver options

    DECinfo%ccsd_expl=.false.
    DECinfo%cc_models(1)='MP2     '
    DECinfo%cc_models(2)='CC2     '
    DECinfo%cc_models(3)='CCSD    '
    DECinfo%cc_models(4)='CCSD(T) '
    DECinfo%t2_restart=.false.
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
    DECinfo%show_time=.false.
    DECinfo%timing=.false.
    DECinfo%show_memory=.false.
    DECinfo%skip_full_ao=.true.
    DECinfo%array4OnFile=.false.
    DECinfo%array4OnFile_specified=.false.
    DECinfo%AObasedCC=.false.


    ! -- Arrays
    DECinfo%zero_threshold = 1E-10_realk


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
    DECinfo%kappaTHR=1e-5
    DECinfo%EerrFactor = 1.0_realk
    DECinfo%EerrOLD = 0.0_realk

    ! -- Timings
    DECinfo%integral_time_cpu=0E0_realk
    DECinfo%integral_time_wall=0E0_realk
    DECinfo%MOintegral_time_cpu=0E0_realk
    DECinfo%MOintegral_time_wall=0E0_realk
    DECinfo%solver_time_cpu=0E0_realk
    DECinfo%solver_time_wall=0E0_realk
    DECinfo%energy_time_cpu=0E0_realk
    DECinfo%energy_time_wall=0E0_realk
    DECinfo%density_time_cpu=0E0_realk
    DECinfo%density_time_wall=0E0_realk
    DECinfo%gradient_time_cpu=0E0_realk
    DECinfo%gradient_time_wall=0E0_realk
    DECinfo%trans_time_cpu=0E0_realk
    DECinfo%trans_time_wall=0E0_realk
    DECinfo%reorder_time_cpu=0E0_realk
    DECinfo%reorder_time_wall=0E0_realk
    DECinfo%memallo_time_cpu=0E0_realk
    DECinfo%memallo_time_wall=0E0_realk

    !> Lagrangian MP2 energy
    DECinfo%lagrangian=.true.

    !> Super fragment information
    DECinfo%SF=.true.
    DECinfo%SF_maxdist = 2.5E0_realk/bohr_to_angstrom ! 2.5 angstrom
    DECinfo%SF_thr = 0.999E0_realk
    DECinfo%SimulateSF=.false.

    !> MPI (undefined by default)
    DECinfo%MPIgroupsize=0

    ! Test stuff
    DECinfo%SaveFragFile=.false.



  end subroutine dec_set_default_config


  !> \brief Read the **DEC (Divide-expand-consolidate coupled cluster) section in 
  !> input file DALTON.INP and set configuration structure accordingly.
  !> \author Kasper Kristensen
  !> \date September 2010
  SUBROUTINE config_dec_input(input,output,readword,word)
    implicit none
    !> Logical for keeping track of when to read
    LOGICAL,intent(inout)                :: READWORD
    !> Logical unit number for DALTON.INP
    integer,intent(in) :: input
    !> Logical unit number for DALTON.OUT
    integer,intent(in) :: output
    character(len=70) :: word
    integer :: fotlevel,one_force_not_enough

    ! Just to be sure, we set the default values before
    ! applying input values.
    call dec_set_default_config(output)

    ! Do DEC calculation
    DECinfo%doDEC=.true.
    DECinfo%output =output

    DO

       IF(READWORD) THEN
          READ (input, '(A40)') WORD
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

          ! NOTE: By default, we assume that the HF has already been carried out.
       case('.ccsd_old'); DECinfo%ccsd_old=.true.
       case('.CCSDsolver_par'); DECinfo%solver_par=.true.
       case('.CCSDdynamic_load'); DECinfo%dyn_load=.true.
       case('.CCSDsaferun'); DECinfo%CCSDsaferun=.true.
       case('.CCSDno_restart'); DECinfo%CCSDno_restart=.true.
       case('.manual_batchsizes'); DECinfo%manual_batchsizes=.true.;read(input,*) DECinfo%ccsdAbatch, DECinfo%ccsdGbatch
       case('.SkipHartreeFock'); DECinfo%doHF=.false.
       case('.hack'); DECinfo%hack=.true.
       case('.hack2'); DECinfo%hack2=.true.
       case('.restart'); DECinfo%restart=.true.           
       case('.TimeBackup'); read(input,*) DECinfo%TimeBackup
       case('.ReadDECorbitals'); DECinfo%read_dec_orbitals=.true.
       case('.mpidebug'); DECinfo%mpidebug=.true.
       case('.MPIsplit'); read(input,*) DECinfo%MPIsplit
       case('.ccFull'); DECinfo%full_molecular_cc=.true. 
       case('.MP2'); DECinfo%ccModel=1; DECinfo%use_singles=.false.
       case('.CC2'); DECinfo%ccModel=2; DECinfo%use_singles=.true.
       case('.CCSD'); DECinfo%ccModel=3; DECinfo%use_singles=.true.
       case('.CCSD(T)'); DECinfo%ccModel=4; DECinfo%use_singles=.true.
       case('.RPA'); DECinfo%ccModel=5; DECinfo%use_singles=.false.
       case('.useCCSDfrag'); DECinfo%use_ccsd_frag=.true.
       !
       case('.ccMaxIter'); read(input,*) DECinfo%ccMaxIter 
       case('.F12'); DECinfo%F12=.true.
       case('.ccThr') 
          read(input,*) DECinfo%ccConvergenceThreshold
          DECinfo%CCthrSpecified=.true.
       case('.Time'); DECinfo%show_time=.true. 
       case('.Timing'); DECinfo%timing=.true. 
       case('.ShowMemory'); DECinfo%show_memory=.true.
       case('.NotPrec'); DECinfo%use_preconditioner=.false.
       case('.NotBPrec'); DECinfo%use_preconditioner_in_b=.false.
       case('.Debug'); DECinfo%fragmentation_debug=.true.;
          DECinfo%dec_driver_debug=.true.;
          DECinfo%cc_driver_debug=.true.
       case('.canonical'); DECinfo%use_canonical=.true.
       case('.ReassignHatoms'); DECinfo%reassignHatoms=.true.
       case('.Mulliken'); DECinfo%mulliken=.true.
       case('.BoughtonPulay'); DECinfo%BoughtonPulay=.true.
       case('.NotFitOrbitals'); DECinfo%FitOrbitals=.false.
       case('.SimpleOrbitalThresh')
          read(input,*) DECinfo%simple_orbital_threshold
          DECinfo%simple_orbital_threshold_set=.true.
       case('.DIIS'); DECinfo%use_crop=.false.  ! use DIIS instead of CROP
       case('.SubSize'); read(input,*) DECinfo%ccMaxDIIS
       case('.restartT2'); DECinfo%t2_restart=.true.           
       case('.MaxIter'); read(input,*) DECinfo%MaxIter

          !> See description of FOT level in set_input_for_fot_level.
          !> Note that if one does not use simple orbital threshold, then
          !> one should manually choose a suitable threshold to determine
          !> the atomic extent (e.g. threshold for Boughton-Pulay procedure).
       case('.FOT')  
          read(input,*) FOTlevel
          call set_input_for_fot_level(FOTlevel)
       case('.FrozenCore') 
          DECinfo%frozencore=.true.
       case('.decPrint'); read(input,*) DECinfo%PL
       case('.mullikenThr'); read(input,*) DECinfo%mulliken_threshold
       case('.skipPairs') 
          DECinfo%pair_distance_threshold=0.0E0_realk
          DECinfo%paircut_set=.true.  ! overwrite default pair cutoff defined by .FOT
       case('.NoExtraPairs') 
          DECinfo%NoExtraPairs=.true.  
       case('.PairRedDist') 
          read(input,*) DECinfo%PairReductionDistance 
       case('.PairRedDistAngstrom') 
          read(input,*) DECinfo%PairReductionDistance 
          DECinfo%PairReductionDistance = DECinfo%PairReductionDistance/bohr_to_angstrom
       case('.pairsThr') 
          read(input,*) DECinfo%pair_distance_threshold
          DECinfo%paircut_set=.true.  ! overwrite default pair cutoff defined by .FOT
       case('.pairsThrAngstrom') ! Input in Angstrom
          read(input,*) DECinfo%pair_distance_threshold
          DECinfo%pair_distance_threshold=DECinfo%pair_distance_threshold/bohr_to_angstrom
          DECinfo%paircut_set=.true.  ! overwrite default pair cutoff defined by .FOT
       case('.PairMinDist'); read(input,*) DECinfo%PairMinDist
       case('.PairMinDistAngstrom')
          read(input,*) DECinfo%PairMinDist
          DECinfo%PairMinDist = DECinfo%PairMinDist/bohr_to_angstrom
       case('.ccsdExpl'); DECinfo%ccsd_expl=.true.
       case('.skipCC'); DECinfo%SkipCC=.true. 
       case('.NormalizeFragment'); DECinfo%NormalizeFragment=.true.
       case('.precWithFull'); DECinfo%precondition_with_full=.true.
       case('.SimpleMullikenThresh'); DECinfo%simple_mulliken_threshold=.true.
       case('.normThresh'); read(input,*) DECinfo%approximated_norm_threshold
       case('.FullDEC'); DECinfo%FullDEC=.true.
       case('.SimulateFull'); DECinfo%simulate_full=.true.
       case('.Simulate_natoms'); read(input,*) DECinfo%simulate_natoms
       case('.SkipReadIn'); DECinfo%SkipReadIn=.true.
       case('.SinglesPolari'); DECinfo%SinglesPolari=.true.
       case('.SinglesThr'); read(input,*) DECinfo%SinglesThr
       case('.SubFrag')
          read(input,*) DECinfo%nsubfrag
          DECinfo%subfrag=.true.
       case('.SubFrag2') 
          read(input,*) DECinfo%nsubfrag2
          DECinfo%subfrag=.true.
       case('.convert64to32')
          DECinfo%convert64to32=.true.
       case('.convert32to64')
          DECinfo%convert32to64=.true.
       case('.IncludeFullMolecule');DECinfo%InclFullMolecule=.true.
       case('.array4OnFile') 
          DECinfo%array4OnFile=.true.
          DECinfo%array4OnFile_specified=.true.
       case('.LagStepSize'); read(input,*) DECinfo%LagStepSize
       case('.AEstepsize'); read(input,*) DECinfo%AEstep
       case('.FragmentAdapted'); DECinfo%fragadapt=.true.
       case('.SaveFragFile'); DECinfo%SaveFragFile=.true.

          ! MP2 response
       case('.gradient') 
          DECinfo%gradient=.true.
          DECinfo%first_order=.true.
          DECinfo%SF=.true.
          DECinfo%Lagrangian=.true.
       case('.kappaMaxIter'); read(input,*) DECinfo%kappaMaxIter 
       case('.kappaMaxDIIS'); read(input,*) DECinfo%kappaMaxDIIS
       case('.kappa_debug'); DECinfo%kappa_driver_debug=.true.
       case('.NotkappaPrec'); DECinfo%kappa_use_preconditioner=.false.
       case('.NotkappaBPrec'); DECinfo%kappa_use_preconditioner_in_b=.false.

          ! integrals tests
       case('.UseFullAO'); DECinfo%skip_full_ao=.false.
       case('.CheckLCM'); DECinfo%check_lcm_orbitals=.true.

          ! Max memory measured in GB. By default set to 16 GB
       case('.Memory') 
          read(input,*) DECinfo%memory           
          DECinfo%memory_defined=.true.
          !> Super fragment calculation
       case('.NotSuperFragment'); DECinfo%SF=.false.
       case('.MaxDistSF'); read(input,*) DECinfo%SF_maxdist
       case('.MaxDistSFAngstrom') 
          read(input,*) DECinfo%SF_maxdist
          DECinfo%SF_maxdist = DECinfo%SF_maxdist/bohr_to_angstrom
       case('.OverlapThrSF'); read(input,*) DECinfo%SF_thr
       case('.SimulateSF'); DECinfo%SimulateSF=.true.
       case('.CCSDforce_scheme'); DECinfo%force_scheme=.true.
               read(input,*) DECinfo%en_mem
       case('.TESTARRAY'); DECinfo%array_test=.true.
       case('.TESTREORDERINGS'); DECinfo%reorder_test=.true.

          ! Size of local groups in MPI scheme
       case('.MPIgroupsize'); read(input,*) DECinfo%MPIgroupsize

          !> Carry out MP2 density fragment calculations
       case('.MP2density') 
          DECinfo%MP2density=.true.
          DECinfo%first_order=.true.
          DECinfo%SF=.true.
          DECinfo%Lagrangian=.true.
          !> Collect fragment contributions to calculate full molecular MP2 density
       case('.SkipFull') 
          DECinfo%SkipFull=.true.
       case('.kappaTHR') 
          read(input,*) DECinfo%kappaTHR
       case('.ErrorFactor') 
          read(input,*) DECinfo%EerrFactor
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


       if(.not. DECinfo%skip_full_ao) then
          call lsquit('Storing arrays on file not implemented for&
               & full AO arrays. Suggestion: Remove .array4OnFile keyword!', DECinfo%output)
       end if

    end if ArraysOnFile


    BeyondMp2: if(DECinfo%ccModel /= 1) then

       !if(DECinfo%skip_full_ao) then
       !   call lsquit('Coupled cluster beyond MP2 is only implemented for full AO arrays. &
       !        & Suggestion: Insert .UseFullAO keyword!', DECinfo%output)
       !end if

       if(DECinfo%MP2density) then
          call lsquit('Calculation of density matrix is only implemented for MP2!', DECinfo%output)
       end if

       if(DECinfo%gradient) then
          call lsquit('Calculation of molecular gradient is only implemented for MP2!', DECinfo%output)
       end if

       write(DECinfo%output,*) 'Coupled-cluster beyond MP2 is requested!'
       write(DECinfo%output,*) 'I turn on the occupied/virtual hybrid partitioning scheme!'
       DECinfo%Lagrangian=.true.
       DECinfo%HybridScheme=.true.

    end if BeyondMp2


    MP2gradientCalculation: if(DECinfo%gradient) then

       if(DECinfo%full_molecular_cc) then
          call lsquit('Full calculation for MP2 gradient is implemented via the &
               & .SimulateFull keyword', DECinfo%output)
       end if

    end if MP2gradientCalculation


    ! If simulate full calculation, the full molecule must be included in "fragments"
    SimulateFullCalc: if(DECinfo%simulate_full) then
       DECinfo%InclFullMolecule = .true.
    end if SimulateFullCalc

    ccsdpt_calc_tmp: if ((DECinfo%ccModel .eq. 4 .and. (.not. DECinfo%SimulateSF)) &
                         & .and. (DECinfo%ccModel .eq. 4 .and. (.not. DECinfo%full_molecular_cc))) then

       call lsquit('temporarily, dec-ccsd(t) must be done in combination &
               & with .SimulateSF or .ccFull keywords',DECinfo%output)

    end if ccsdpt_calc_tmp    

    if(DECinfo%ccmodel==4 .and. DECinfo%restart .and. (.not. DECinfo%use_ccsd_frag)) then
       call lsquit('Restart option currently not implemented for CCSD(T)!',DECinfo%output)
    end if

    ! Fragment-adapted orbitals (FAOs) only work when super fragments 
    ! are just simulated and not actually used.
    ! Also FAOs do not work with reduced pairs, set reduction distance to 1000000 to
    ! avoid it from being used in practice
    if(DECinfo%fragadapt) then
       DECinfo%simulateSF = .true.
       DECinfo%PairReductionDistance = 1.0e6_realk
    end if


    ! Set CC residual threshold to be 0.01*FOT
    ! - unless it was specified explicitly in the input.
    if(.not. DECinfo%CCthrSpecified) then
       DECinfo%ccConvergenceThreshold=0.01E0_realk*DECinfo%FOT
    end if

    ! Only full molecular for RPA at this stage
    if(DECinfo%ccmodel==5 .and. .not. DECinfo%full_molecular_cc) then
       call lsquit('RPA only implemented for full molecule! Insert .ccFull keyword.',-1)
    end if

    ! Never use gradient and density at the same time (density is a subset of gradient)
    if(DECinfo%MP2density .and. DECinfo%gradient) then
       call lsquit('Density and gradient cannot both be turned on at the same time! &
            & Note that density is a subset of a gradient calculation',DECinfo%output)
    end if

#ifndef VAR_LSMPI
if(DECinfo%restart) then
       call lsquit('DEC Restart option only possible using MPI!',DECinfo%output)
end if
#endif

if(DECinfo%frozencore .and. DECinfo%singlespolari) then
call lsquit('Frozen core not implemented with long range singles correction!',-1)
end if

#ifdef VAR_LSMPI
if(DECinfo%SinglesPolari) then
call lsquit('Full singles polarization has been temporarily disabled for MPI',-1)
end if
#endif


  end subroutine check_dec_input
  
  subroutine DEC_settings_print(DECitem,lupri)
    type(DECsettings) :: DECitem
    integer,intent(in) :: lupri
    !
    integer :: I,J,K,L,M
    WRITE(lupri,'(A)') ' The DEC settings structure'
    WRITE(lupri,'(A32,L1)')' doDEC ' ,DECitem%doDEC
    WRITE(lupri,'(A32,L1)')' doHF ' ,DECitem%doHF
    WRITE(lupri,'(A32,ES18.9)')' memory ' ,DECitem%memory
    WRITE(lupri,'(A32,L1)')' full_molecular_cc ' ,DECitem%full_molecular_cc
    WRITE(lupri,'(A32,L1)')' FullDEC ' ,DECitem%FullDEC
    WRITE(lupri,'(A32,L1)')' simulate_full ' ,DECitem%simulate_full
    WRITE(lupri,'(A32,I8)')' simulate_natoms ' ,DECitem%simulate_natoms
    WRITE(lupri,'(A32,L1)')' SkipReadIn ' ,DECitem%SkipReadIn
    WRITE(lupri,'(A32,L1)')' SinglesPolari ' ,DECitem%SinglesPolari
    WRITE(lupri,'(A32,ES18.9)')' singlesthr ' ,DECitem%singlesthr
    WRITE(lupri,'(A32,I8)')' nsubfrag ' ,DECitem%nsubfrag
    WRITE(lupri,'(A32,I8)')' nsubfrag2 ' ,DECitem%nsubfrag2
    WRITE(lupri,'(A32,L1)')' subfrag ' ,DECitem%subfrag
    WRITE(lupri,'(A32,L1)')' restart ' ,DECitem%restart
    WRITE(lupri,'(A32,ES18.9)')' TimeBackup ' ,DECitem%TimeBackup
    WRITE(lupri,'(A32,L1)')' fragmentation_debug ' ,DECitem%fragmentation_debug
    WRITE(lupri,'(A32,L1)')' dec_driver_debug ' ,DECitem%dec_driver_debug
    WRITE(lupri,'(A32,L1)')' cc_driver_debug ' ,DECitem%cc_driver_debug
    WRITE(lupri,'(A32,L1)')' ccsd_old, ' ,DECitem%ccsd_old
    WRITE(lupri,'(A32,L1)')' solver_par, ' ,DECitem%solver_par
    WRITE(lupri,'(A32,I8)')' ccsdAbatch' ,DECitem%ccsdAbatch
    WRITE(lupri,'(A32,I8)')' ccsdGbatch' ,DECitem%ccsdGbatch
    WRITE(lupri,'(A32,L1)')' hack ' ,DECitem%hack
    WRITE(lupri,'(A32,L1)')' hack2 ' ,DECitem%hack2
    WRITE(lupri,'(A32,L1)')' mpidebug ' ,DECitem%mpidebug
    WRITE(lupri,'(A32,I8)')' output ' ,DECitem%output
    WRITE(lupri,'(A32,ES18.9)')' mulliken_threshold ' ,DECitem%mulliken_threshold
    WRITE(lupri,'(A32,L1)')' simple_mulliken_threshold ' ,DECitem%simple_mulliken_threshold
    WRITE(lupri,'(A32,ES18.9)')' approximated_norm_threshold ' ,DECitem%approximated_norm_threshold
    WRITE(lupri,'(A32,L1)')' check_lcm_orbitals ' ,DECitem%check_lcm_orbitals
    WRITE(lupri,'(A32,L1)')' use_canonical ' ,DECitem%use_canonical
    WRITE(lupri,'(A32,L1)')' ReassignHatoms ' ,DECitem%ReassignHatoms
    WRITE(lupri,'(A32,L1)')' mulliken ' ,DECitem%mulliken
    WRITE(lupri,'(A32,L1)')' BoughtonPulay ' ,DECitem%BoughtonPulay
    WRITE(lupri,'(A32,L1)')' FitOrbitals ' ,DECitem%FitOrbitals
    WRITE(lupri,'(A32,ES18.9)')' simple_orbital_threshold ' ,DECitem%simple_orbital_threshold
    WRITE(lupri,'(A32,I8)')' MaxIter ' ,DECitem%MaxIter
    WRITE(lupri,'(A32,ES18.9)')' FOT ' ,DECitem%FOT
    WRITE(lupri,'(A32,I8)')' PL ' ,DECitem%PL
    WRITE(lupri,'(A32,L1)')' SkipCC ' ,DECitem%SkipCC
    WRITE(lupri,'(A32,L1)')' NormalizeFragment ' ,DECitem%NormalizeFragment
    WRITE(lupri,'(A32,L1)')' precondition_with_full ' ,DECitem%precondition_with_full
    WRITE(lupri,'(A32,L1)')' InclFullMolecule ' ,DECitem%InclFullMolecule
    WRITE(lupri,'(A32,L1)')' HybridScheme ' ,DECitem%HybridScheme
    WRITE(lupri,'(A32,I8)')' LagStepSize ' ,DECitem%LagStepSize
    WRITE(lupri,'(A32,ES18.9)')' pair_distance_threshold ' ,DECitem%pair_distance_threshold
    WRITE(lupri,'(A32,ES18.9)')' PairReductionDistance ' ,DECitem%PairReductionDistance
    WRITE(lupri,'(A32,ES18.9)')' fullmolecule_memory ' ,DECitem%fullmolecule_memory
    WRITE(lupri,'(A32,L1)')' ccsd_expl ' ,DECitem%ccsd_expl
    WRITE(lupri,'(A32,10A)')' cc_models ' ,DECitem%cc_models
    WRITE(lupri,'(A32,L1)')' t2_restart ' ,DECitem%t2_restart
    WRITE(lupri,'(A32,L1)')' simulate_eri ' ,DECitem%simulate_eri
    WRITE(lupri,'(A32,L1)')' fock_with_ri ' ,DECitem%fock_with_ri
    WRITE(lupri,'(A32,I8)')' ccMaxIter ' ,DECitem%ccMaxIter
    WRITE(lupri,'(A32,I8)')' ccMaxDIIS ' ,DECitem%ccMaxDIIS
    WRITE(lupri,'(A32,I8)')' ccModel ' ,DECitem%ccModel
    WRITE(lupri,'(A32,L1)')' F12 ' ,DECitem%F12
    WRITE(lupri,'(A32,ES18.9)')' ccConvergenceThreshold ' ,DECitem%ccConvergenceThreshold
    WRITE(lupri,'(A32,L1)')' CCthrSpecified ' ,DECitem%CCthrSpecified
    WRITE(lupri,'(A32,L1)')' use_singles ' ,DECitem%use_singles
    WRITE(lupri,'(A32,L1)')' use_preconditioner ' ,DECitem%use_preconditioner
    WRITE(lupri,'(A32,L1)')' use_preconditioner_in_b ' ,DECitem%use_preconditioner_in_b
    WRITE(lupri,'(A32,L1)')' use_crop ' ,DECitem%use_crop
    WRITE(lupri,'(A32,L1)')' show_time ' ,DECitem%show_time
    WRITE(lupri,'(A32,L1)')' timing ' ,DECitem%timing
    WRITE(lupri,'(A32,L1)')' show_memory ' ,DECitem%show_memory
    WRITE(lupri,'(A32,L1)')' skip_full_ao ' ,DECitem%skip_full_ao
    WRITE(lupri,'(A32,L1)')' array4OnFile ' ,DECitem%array4OnFile
    WRITE(lupri,'(A32,L1)')' array4OnFile_specified ' ,DECitem%array4OnFile_specified
    WRITE(lupri,'(A32,L1)')' AObasedCC ' ,DECitem%AObasedCC
    WRITE(lupri,'(A32,ES18.9)')' zero_threshold ' ,DECitem%zero_threshold
    WRITE(lupri,'(A32,ES18.9)')' integral_time_cpu ' ,DECitem%integral_time_cpu
    WRITE(lupri,'(A32,ES18.9)')' integral_time_wall ' ,DECitem%integral_time_wall
    WRITE(lupri,'(A32,ES18.9)')' MOintegral_time_cpu ' ,DECitem%MOintegral_time_cpu
    WRITE(lupri,'(A32,ES18.9)')' MOintegral_time_wall ' ,DECitem%MOintegral_time_wall
    WRITE(lupri,'(A32,ES18.9)')' solver_time_cpu ' ,DECitem%solver_time_cpu
    WRITE(lupri,'(A32,ES18.9)')' solver_time_wall ' ,DECitem%solver_time_wall
    WRITE(lupri,'(A32,ES18.9)')' energy_time_cpu ' ,DECitem%energy_time_cpu
    WRITE(lupri,'(A32,ES18.9)')' energy_time_wall ' ,DECitem%energy_time_wall
    WRITE(lupri,'(A32,ES18.9)')' density_time_cpu ' ,DECitem%density_time_cpu
    WRITE(lupri,'(A32,ES18.9)')' density_time_wall ' ,DECitem%density_time_wall
    WRITE(lupri,'(A32,ES18.9)')' gradient_time_cpu ' ,DECitem%gradient_time_cpu
    WRITE(lupri,'(A32,ES18.9)')' gradient_time_wall ' ,DECitem%gradient_time_wall
    WRITE(lupri,'(A32,ES18.9)')' trans_time_cpu ' ,DECitem%trans_time_cpu
    WRITE(lupri,'(A32,ES18.9)')' trans_time_wall ' ,DECitem%trans_time_wall
    WRITE(lupri,'(A32,ES18.9)')' reorder_time_cpu ' ,DECitem%reorder_time_cpu
    WRITE(lupri,'(A32,ES18.9)')' reorder_time_wall ' ,DECitem%reorder_time_wall
    WRITE(lupri,'(A32,ES18.9)')' memallo_time_cpu ' ,DECitem%memallo_time_cpu
    WRITE(lupri,'(A32,ES18.9)')' memallo_time_wall ' ,DECitem%memallo_time_wall
    WRITE(lupri,'(A32,L1)')' lagrangian ' ,DECitem%lagrangian
    WRITE(lupri,'(A32,L1)')' SF ' ,DECitem%SF
    WRITE(lupri,'(A32,ES18.9)')' SF_maxdist ' ,DECitem%SF_maxdist
    WRITE(lupri,'(A31,ES18.9)')' SF_THR ' ,DECitem%SF_THR
    WRITE(lupri,'(A31,L1)')' SimulateSF ' ,DECitem%SimulateSF
    WRITE(lupri,'(A31,I8)')' MPIgroupsize ' ,DECitem%MPIgroupsize
    WRITE(lupri,'(A31,L1)')' first_order ' ,DECitem%first_order
    WRITE(lupri,'(A31,L1)')' MP2density ' ,DECitem%MP2density
    WRITE(lupri,'(A31,L1)')' gradient ' ,DECitem%gradient
    WRITE(lupri,'(A31,L1)')' kappa_use_preconditioner ' ,DECitem%kappa_use_preconditioner
    WRITE(lupri,'(A31,L1)')' kappa_use_preconditioner_in_b ' ,DECitem%kappa_use_preconditioner_in_b
    WRITE(lupri,'(A31,I8)')' kappaMaxDIIS ' ,DECitem%kappaMaxDIIS
    WRITE(lupri,'(A31,I8)')' kappaMaxIter ' ,DECitem%kappaMaxIter
    WRITE(lupri,'(A31,L1)')' kappa_driver_debug ' ,DECitem%kappa_driver_debug
    WRITE(lupri,'(A31,ES18.9)')' kappaTHR ' ,DECitem%kappaTHR
    WRITE(lupri,'(A31,L1)')' SaveFragFile ' ,DECitem%SaveFragFile
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
    !>           Reasonable pair cutoff distances and 
    !>           atomic extent orbital thresholds are associated
    !>           with each FOT level.
    !>           If necessary, the paircut off is increased
    !>           in a self-adaptive black box manner during
    !>           the calculation.
    !> 
    !> FOTlevel = 1: 
    !> FOT=10^{-1}, pair_distance_threshold=4Angstrom, simple_orbital_threshold=0.1
    !>
    !> FOTlevel = 2: 
    !> FOT=10^{-2}, pair_distance_threshold=6Angstrom, simple_orbital_threshold=0.1
    !> 
    !> FOTlevel = 3: 
    !> FOT=10^{-3}, pair_distance_threshold=8Angstrom, simple_orbital_threshold=0.03
    !> 
    !> FOTlevel = 4: 
    !> FOT=10^{-4}, pair_distance_threshold=10Angstrom, simple_orbital_threshold=0.01
    !>
    !> FOTlevel = 5: 
    !> FOT=10^{-5}, pair_distance_threshold=12Angstrom, simple_orbital_threshold=0.003
    !> 
    !> FOTlevel = 6: 
    !> FOT=10^{-6}, pair_distance_threshold=14Angstrom, simple_orbital_threshold=0.001
    !>
    !> FOTlevel = 7: 
    !> FOT=10^{-7}, pair_distance_threshold=16Angstrom, simple_orbital_threshold=0.0003
    !> 
    !> Default: FOTlevel=4. If FOTlevel is not from 1 to 7, the program will quit.

    ! Set FOT level
    DECinfo%FOTlevel = FOTlevel

    WhatFOTlevel: SELECT CASE(DECinfo%FOTlevel)

    case(1)
       DECinfo%FOT = 1.0E-1_realk
       ! Use FOT-adapted pair cutoff and simple orbital threshold
       ! only if these were not set set manually
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=4.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.1E0_realk
       end if

    case(2)
       DECinfo%FOT = 1.0E-2_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=6.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.1E0_realk
       end if

    case(3)
       DECinfo%FOT = 1.0E-3_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=8.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.03E0_realk
       end if

    case(4)
       DECinfo%FOT = 1.0E-4_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=10.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.01E0_realk
       end if

    case(5)
       DECinfo%FOT = 1.0E-5_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=12.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.003E0_realk
       end if

    case(6)
       DECinfo%FOT = 1.0E-6_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=14.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.001E0_realk
       end if

    case(7)
       DECinfo%FOT = 1.0E-7_realk
       if(.not. DECinfo%paircut_set) DECinfo%pair_distance_threshold=16.0E0_realk/bohr_to_angstrom
       if(.not. DECinfo%simple_orbital_threshold_set) then
          DECinfo%simple_orbital_threshold = 0.0003E0_realk
       end if

    case default
       call lsquit('set_input_for_fot_level: FOT level must be 1,2,3,4,5,6, or 7!',DECinfo%output)

    end SELECT WhatFOTlevel


  end subroutine set_input_for_fot_level

end MODULE DEC_settings_mod
